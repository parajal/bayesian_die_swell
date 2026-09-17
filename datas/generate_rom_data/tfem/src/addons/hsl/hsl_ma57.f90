
! Copyright (C) 2006-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, routines, for solving the system of equations using the direct solver
! MA57 from the HSL library (http://www.numerical.rl.ac.uk/hsl).
! The sysmatrix must be defined as symmetric.

module hsl_ma57_m

  use kind_defs_m
  use sparse_m, only: sparsematrix_t
  use system_defs_m

  implicit none


! type definition for solver options

  type solver_options_ma57_t

    integer :: printlevel = 0  ! 0 no printing
                               ! 1 print info, rinfo
                               ! 2 print a lot of information from the solver
                               ! 3 print minimum memory requirements after
                               !   analysis and stop

    integer :: scaling = 0  ! scaling of the matrix before factorization:
                            ! ICNTL(15) of MA57. Set to 1 for scaling the
                            ! matrix.
                            ! NOTE: the default here is no scaling, although
                            !       ma57 has scaling set by default.
                            ! NOTE: this option is NOT available for the ma57
                            !       routine from the HSL2002 library as
                            !       supplied with tfem. See the HSL website
                            !       on how to obtain a recent version (free
                            !       for academics).

    integer :: niter_ref = 0  ! Maximum permitted number of steps of iterative
                              ! refinement. If niter_ref=0 no iterative
                              ! refinement is performed.

    real(dp) :: conv_rate = 0.5_dp  ! If the norm of the scaled residuals does
                                    ! not decrease by a factor of at least
                                    ! conv_rate, iterative refinement terminates

    real(dp) :: real_storage = 1._dp ! work storage allocated for reals relative
                                     ! to the minimum required according to
                                     ! the analysis step

    real(dp) :: integer_storage = 1.2_dp ! work storage allocated for integers
                                         ! relative to the minimum required
                                         ! according to the analysis step

!   Used for threshold pivoting: CNTL(1) of MA57.
    real(dp) :: pivotthreshold = 0.01_dp

  end type solver_options_ma57_t


! type definition of storage for LU factorization (MA57)

  type lu_ma57_t

    logical :: decomp = .false.  ! LU decomposition done?

!   icntl array (MA57)
    integer, allocatable, dimension(:) :: icntl

!   ifact array (MA57)
    integer, allocatable, dimension(:) :: ifact

!   fact array (MA57)
    real(dp), allocatable, dimension(:) :: fact

  end type lu_ma57_t


! interface for generic delete subroutine

  interface delete
    module procedure delete_lu_ma57
  end interface delete

! interface for generic set_solver_options subroutine

  interface set_solver_options
    module procedure set_solver_options_ma57
  end interface set_solver_options


contains


! Helper routine for setting the solver options

  subroutine set_solver_options_ma57 ( solver_options, keep, printlevel, &
    niter_ref, conv_rate, real_storage, integer_storage, pivotthreshold )

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_ma57_t), intent(inout) :: solver_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see solver_options_ma41_t
    integer, intent(in), optional :: printlevel, niter_ref
    real(dp), intent(in), optional :: conv_rate, real_storage, integer_storage,&
      pivotthreshold


    type(solver_options_ma57_t) :: solveropt


    if ( present(keep) ) then
      if ( .not. keep ) solver_options = solveropt
    else
      solver_options = solveropt
    end if

    if ( present(printlevel) ) solver_options%printlevel = printlevel
    if ( present(niter_ref) ) solver_options%niter_ref = niter_ref
    if ( present(conv_rate) ) solver_options%conv_rate = conv_rate
    if ( present(real_storage) ) solver_options%real_storage = real_storage
    if ( present(integer_storage) ) solver_options%integer_storage = &
                                                 integer_storage
    if ( present(pivotthreshold) ) solver_options%pivotthreshold = &
                                                 pivotthreshold

  end subroutine set_solver_options_ma57


! Solve the system of equations using MA57 (a direct sparse solver for
! unsymmetric matrices)

  subroutine solve_system_ma57 ( sysmatrix, rhsd, sol, lu, solver_options )

    type(sysmatrix_t), intent(in) :: sysmatrix

!   Right-hand side vector. The part corresponding to the unknowns is used as
!   the right-hand side in the system:
!
!      Suu sol_u = rhsd_u   (*)
!
!   The effect of essential boundary conditions must already have been taken
!   into account in the rhsd vector
    type(sysvector_t), intent(in) :: rhsd

!   The solution vector. The unknown part sol_u is filled with the solution of
!   the system (*), given above. The prescribed part sol_p is not modified.
    type(sysvector_t), intent(inout) :: sol

!   When present, lu is used to store the LU factorization. The next time the
!   solve routine is called with lu present only a backsubstition is done. This
!   happens even when the matrix has changed. If lu is not present or
!   has been deallocated using delete, a full factorization is performed.
    type(lu_ma57_t), intent(inout), optional :: lu

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_ma57_t), intent(in), optional :: solver_options


    integer :: numundegfd
    type(lu_ma57_t) :: lw  ! local work space when lu is not present in heading
    type(solver_options_ma57_t) :: lsolveropt ! local options


!   some testing first

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(/a/)') &
        'Error in solve_system_ma57: no system matrix structure.'
      stop
    end if

    if ( .not. sysmatrix%finalized ) then
      write(*,'(/a/)') &
        'Error in solve_system_ma57: system matrix has not been finalized.'
      stop
    end if

    if ( .not. sysmatrix%allocated_data ) then
      write(*,'(/a/)') &
        'Error in solve_system_ma57: data in system matrix not allocated.'
      stop
    end if

    if ( .not. sysmatrix%symmetric ) then
      write(*,'(/2(a/))') &
        'Error in solve_system_ma57: ', &
        '  sysmatrix must have been created with symmetric=.true.'
      stop
    end if

    if ( .not. rhsd%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_ma57: rhsd not created.'
      stop
    end if

    if ( .not. sol%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_ma57: sol not created.'
      stop
    end if

!   initialize local solver options if solver_options in heading
    if ( present(solver_options) ) lsolveropt = solver_options

!   number of unknowns

    numundegfd = sysmatrix%Suu%n

!   copy first part of rhsd to sol

    sol%u(1:numundegfd) = rhsd%u(1:numundegfd)

!   solve, solution returned in sol(1:numundegfd)

    if ( present(lu) ) then

      call solve_ma57 ( sysmatrix%Suu, sol%u(1:numundegfd), lu, lsolveropt )

    else

      call solve_ma57 ( sysmatrix%Suu, sol%u(1:numundegfd), lw, lsolveropt )

!     remove temporary storage
      call delete ( lw )

    end if

  end subroutine solve_system_ma57


! solve linear system A u = b

  subroutine solve_ma57 ( A, b, lu, so )

!   The matrix A in CSR format
    type(sparsematrix_t), intent(in) :: A

!   On entry: the right-hand side vector b.
!   On exit: the solution u = A^{-1} b.
    real(dp), intent(inout), dimension(:) :: b

!   lu is used to store the LU factorization.
!   If lu contains a LU factorization from a previous solve only a
!   backsubstition is done.  This happens even when the matrix has changed.
!   If lu has been deallocated using delete, a full factorization is performed.
    type(lu_ma57_t), intent(inout) :: lu

!   The options/parameters for the solver
    type(solver_options_ma57_t), intent(in) :: so


    logical :: allsteps
    integer :: n, ne, job, lkeep, lwork, lfact, lifact, lrhs, nrhs, row
    integer, dimension(40) :: info
    integer, dimension(:), allocatable :: keep, iwork, irn
    real(dp) :: cntl(5)
    real(dp), dimension(20) :: rinfo
    real(dp), dimension(:), allocatable :: work, x, resid


!   initialize local parameters

    allsteps = .true.
    info = 0
    rinfo = 0

!   set dimensions of problem

    n  = A%n
    ne = A%nnz

!   do all steps (analysis, factorization, solution) or only solution?

    if ( lu%decomp ) then
!     decomposition is already done
      allsteps = .false.
    end if

!   choose job

    if ( allsteps ) then

!     all steps

!     allocate memory and do the analysis

      call lu_allocate_analysis_ma57 ( lu )

      call prininfo ( 'MA57AD', 31, 15, info, rinfo, so%printlevel )

!     allocate memory for factorisation

      allocate ( lu%fact(lfact), lu%ifact(lifact) )
      allocate ( iwork(n) )

!     factorisation

      cntl(1) = so%pivotthreshold   ! threshold for pivoting

      call ma57bd ( n, ne, A%a, lu%fact, lfact, lu%ifact, lifact, lkeep, &
        keep, iwork, lu%icntl, cntl, info, rinfo )

      deallocate ( keep, iwork )

      call prininfo ( 'MA57BD', 31, 15, info, rinfo, so%printlevel )

    end if

!   solve system

    if ( so%niter_ref==0 ) then

!     standard case without iterative refinement

!     allocate local memory
      lwork = n
      allocate ( work(lwork), iwork(n) )

      job = 1
      nrhs = 1
      lrhs = n

      call ma57cd ( job, n, lu%fact, size(lu%fact), lu%ifact, size(lu%ifact), &
        nrhs, b, lrhs, work, lwork, iwork, lu%icntl, info )

      call prininfo ( 'MA57CD', 31, 15, info, rinfo, so%printlevel )

      deallocate ( work, iwork )

    else

!     use iterative refinement

!     convergence rate of iterative refinement
      cntl(3) = so%conv_rate

!     allocate local memory
      if ( so%niter_ref > 1 ) then
        job = 0
        allocate ( work(3*n) )
      else
        job = 1
        allocate ( work(n) )
      end if

      allocate ( x(n), resid(n) )
      allocate ( irn(ne) )
      allocate ( iwork(n) )

!     fill row number array to transform CSR -> coordinate format
      do row = 1, n
        irn( A%ia(row):A%ia(row+1)-1 ) = row
      end do

      call ma57dd ( job, n, ne, A%a, irn, A%ja, lu%fact, size(lu%fact), &
        lu%ifact, size(lu%ifact), b, x, resid, work, iwork, lu%icntl, cntl, &
        info, rinfo )

      b = x

      call prininfo ( 'MA57DD', 31, 15, info, rinfo, so%printlevel )

      deallocate ( work, x, resid, irn, iwork )

    end if

!   for next call of this routine: do solution only (no factorization) if
!   lu is not deleted

    lu%decomp = .true.


  contains


!   allocate memory for lu and do analysis step

    subroutine lu_allocate_analysis_ma57 ( lu )

      type (lu_ma57_t) :: lu

      integer :: row
      integer, dimension(:), allocatable :: irn

!     allocate local memory

      lkeep  = 5 * n + ne + max(n,ne) + 42 + 2 * n

      allocate ( keep(lkeep), iwork(5*n), irn(ne) )

!     allocate memory for lu

      allocate ( lu%icntl(20) )

!     default values

      call ma57id ( cntl, lu%icntl )

!     if iterative refinement is used ICNTL(9) contains the
!     maximum number of steps. Default number is 10

      lu%icntl(9) = so%niter_ref

      if ( so%printlevel == 2 ) then
!       do a lot of printing
        lu%icntl(4) = 6
        lu%icntl(5) = 3
      end if

      if ( so%scaling > 0 ) then
!       set scaling of the matrix
        lu%icntl(15) = 1
      else
!       no scaling of the matrix (default)
        lu%icntl(15) = 0
      end if

!     fill row number array to transform CSR -> coordinate format

      do row = 1, n
        irn( A%ia(row):A%ia(row+1)-1 ) = row
      end do

      call ma57ad ( n, ne, irn, A%ja, lkeep, keep, iwork, lu%icntl, info, &
        rinfo )

      deallocate ( iwork, irn )

      if ( so%printlevel == 3 ) then

!       print memory requirement and stop

        call prinmem

        stop

      end if

!     real storage space

      lfact = nint( so%real_storage * info(9) )
      lifact = nint( so%integer_storage * info(10) )

    end subroutine lu_allocate_analysis_ma57


!   print memory requirements

    subroutine prinmem

      real(dp) :: ne, info9, info10

      ne = real(A%nnz,kind=dp)
      info9 = real(info(9),kind=dp)
      info10 = real(info(10),kind=dp)

      write(*,'(/a//3(2(a,i0),a/2(a,i0,a/)/))') &
        'Minimum storage requirements after MA57 analysis:', &
        '  Matrix: ', nint(ne), ' integers ', nint(ne), ' reals', &
        '    Using 32bit integers: ', &
             nint((ne*4+ne*8)/1000000), ' MBytes ', &
        '    Using 64bit integers: ', &
             nint((ne*8+ne*8)/1000000), ' MBytes ', &
        '  Factorization: ', info(10), ' integers ', info(9), ' reals', &
        '    Using 32bit integers: ', &
             nint((info10*4+info9*8)/1000000), ' MBytes ', &
        '    Using 64bit integers: ', &
             nint((info10*8+info9*8)/1000000), ' MBytes ', &
        '  Total: ', nint(ne+info10), ' integers ', &
             nint(ne+info9), ' reals', &
        '    Using 32bit integers: ', &
          nint(((ne+info10)*4+(ne+info9)*8)/1000000), ' MBytes ', &
        '    Using 64bit integers: ', &
          nint(((ne+info10)*8+(ne+info9)*8)/1000000), ' MBytes '

      end subroutine prinmem

  end subroutine solve_ma57


! delete single lu

  subroutine delete_single_lu_ma57 ( lu )

    type (lu_ma57_t) :: lu

    if ( .not. allocated(lu%icntl) ) then
      write(*,'(/a/)') &
        'Error delete_single_lu_ma57: arrays in lu not allocated'
      stop
    end if

    deallocate ( lu%icntl, lu%ifact )
    deallocate ( lu%fact )

    lu%decomp = .false.

  end subroutine delete_single_lu_ma57


! Delete lu

  subroutine delete_lu_ma57 ( lu_ma57_1, lu_ma57_2, lu_ma57_3, &
    lu_ma57_4, lu_ma57_5 )

    type(lu_ma57_t), intent(inout) :: lu_ma57_1
    type(lu_ma57_t), intent(inout), optional :: lu_ma57_2, lu_ma57_3, &
      lu_ma57_4, lu_ma57_5

    call delete_single_lu_ma57(lu_ma57_1)
    if ( present(lu_ma57_2) ) call delete_single_lu_ma57(lu_ma57_2)
    if ( present(lu_ma57_3) ) call delete_single_lu_ma57(lu_ma57_3)
    if ( present(lu_ma57_4) ) call delete_single_lu_ma57(lu_ma57_4)
    if ( present(lu_ma57_5) ) call delete_single_lu_ma57(lu_ma57_5)

  end subroutine delete_lu_ma57


! prininfo

  subroutine prininfo ( sub, ni, nr, info, rinfo, printlevel )

    character (len=*), intent(in) :: sub
    integer, intent(in) ::  ni, nr, printlevel
    integer, intent(in), dimension(:) :: info
    real(dp), intent(in), dimension(:) :: rinfo

    if ( info(1) < 0 ) then

!     Error in sub

      write(*,'(/3a,i0/)') 'Error in ', sub, ' INFO(1) = ', info(1)

      call prinfo

      stop

    else if ( printlevel == 1 ) then

!     Print info, rinfo

      write(*,'(/3a,i0/)') 'Info after subroutine ', sub

      call prinfo

    end if

  contains

    subroutine prinfo

      integer :: i

      write(*,'(/a)') 'INFO:'

      do i = 1, ni
         write(*,'(i5,1x,i0)') i, info(i)
      end do

      write(*,'(/a)') 'RINFO:'

      do i = 1, nr
         write(*,*) i, rinfo(i)
      end do

      write(*,'()')

    end subroutine prinfo

  end subroutine prininfo

end module hsl_ma57_m
