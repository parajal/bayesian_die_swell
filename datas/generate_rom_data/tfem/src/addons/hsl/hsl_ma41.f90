
! Copyright (C) 2004-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, routines, for solving the system of equations using the direct solver
! MA41 from the HSL library (http://www.numerical.rl.ac.uk/hsl)

module hsl_ma41_m

  use kind_defs_m
  use sparse_m, only: sparsematrix_t
  use system_defs_m

  implicit none


! type definition for solver options

  type solver_options_ma41_t

    integer :: printlevel = 0  ! 0 no printing
                               ! 1 print info, rinfo
                               ! 2 print a lot of information from the solver
                               ! 3 print minimum memory requirements after
                               !   analysis and stop

    integer :: pivot_order = 0 ! Sets the pivot order = ICNTL(7) of MA41.
                               ! Has default value 0 and must be set by the
                               ! user to 1 to use the pivot order given by
                               ! the ipermu permutation array in the
                               ! sysmatrix. The ipermu array can e.g. be
                               ! set by MeTiS and can significantly speedup the
                               ! solve for large problems.

    integer :: scaling = 0     ! scaling of the matrix before factorization:
                               ! ICNTL(8) of MA41

    integer :: niter_ref = 0   ! maximum  number  of  steps  of  iterative
                               ! refinement.  If niter_ref=0, iterative
                               ! refinement is not performed

    real(dp) :: iter_ref_tol = 1.e-10_dp ! tolerance for the convergence of
                                         ! iterative refinement

    real(dp) :: real_storage = 1._dp ! work storage allocated for reals relative
                                     ! to the minimum required according to
                                     ! the analysis step

    real(dp) :: integer_storage = 1.2_dp ! work storage allocated for integers
                                         ! relative to the minimum required
                                         ! according to the analysis step

!   Used for threshold pivoting: CNTL(1). From the manual of MA41:
!   If it is nonzero, numerical pivoting will be performed. In general, a
!   larger value of CNTL(1) leads to greater fill-in but a more accurate
!   factorization. If CNTL(1) is zero, no pivoting will be performed and the
!   subroutine will fail if a zero pivot is encountered. If the matrix A is
!   diagonally dominant, then setting CNTL(1) to zero will decrease the
!   factorization time while still providing a stable decomposition. Values
!   greater than 1.0 are treated as 1.0, and less than zero as zero.
    real(dp) :: pivotthreshold = 0.01_dp

  end type solver_options_ma41_t


! type definition of storage for LU factorization (MA41)

  type lu_ma41_t

    logical :: decomp = .false.  ! LU decomposition done?

!   irn array (MA41)
    integer, allocatable, dimension(:) :: irn

!   keep array (MA41)
    integer, allocatable, dimension(:) :: keep

!   icntl array (MA41)
    integer, allocatable, dimension(:) :: icntl

!   is array (MA41)
    integer, allocatable, dimension(:) :: is

!   cntl array (MA41)
    real(dp), allocatable, dimension(:) :: cntl

!   colsca array (MA41)
    real(dp), allocatable, dimension(:) :: colsca

!   rowsca array (MA41)
    real(dp), allocatable, dimension(:) :: rowsca

!   s array (MA41)
    real(dp), allocatable, dimension(:) :: s

  end type lu_ma41_t



! Interfaces of external HSL routines

  interface
    subroutine ma41ad ( job, n, ne, irn, jcn, aspk, rhs, colsca, rowsca, keep, &
      is, maxis, s, maxs, cntl, icntl, info, rinfo )
      use kind_defs_m
      integer :: maxs, job, n,ne, maxis, is(maxis), irn(ne), jcn(ne)
      integer :: info(20), keep(50), icntl(20)
      real(dp) :: aspk(ne), rhs(n), colsca(*), rowsca(*)
      real(dp) :: rinfo(20), s(maxs), cntl(10)
    end subroutine ma41ad
    subroutine  ma41id ( cntl, icntl, keep )
      use kind_defs_m
      integer icntl(20), keep(50)
      real(dp) cntl(10)
    end subroutine ma41id
  end interface


! interface for generic delete subroutine

  interface delete
    module procedure delete_lu_ma41
  end interface delete

! interface for generic set_solver_options subroutine

  interface set_solver_options
    module procedure set_solver_options_ma41
  end interface set_solver_options


contains


! Helper routine for setting the solver options

  subroutine set_solver_options_ma41 ( solver_options, keep, printlevel, &
    scaling, pivot_order, niter_ref, iter_ref_tol, real_storage, &
    integer_storage, pivotthreshold )

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_ma41_t), intent(inout) :: solver_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see solver_options_ma41_t
    integer, intent(in), optional :: printlevel, scaling, pivot_order, niter_ref
    real(dp), intent(in), optional :: iter_ref_tol, real_storage, &
      integer_storage, pivotthreshold


    type(solver_options_ma41_t) :: solveropt


    if ( present(keep) ) then
      if ( .not. keep ) solver_options = solveropt
    else
      solver_options = solveropt
    end if

    if ( present(printlevel) ) solver_options%printlevel = printlevel
    if ( present(scaling) ) solver_options%scaling = scaling
    if ( present(pivot_order) ) solver_options%pivot_order = pivot_order
    if ( present(niter_ref) ) solver_options%niter_ref = niter_ref
    if ( present(iter_ref_tol) ) solver_options%iter_ref_tol = iter_ref_tol
    if ( present(real_storage) ) solver_options%real_storage = real_storage
    if ( present(integer_storage) ) solver_options%integer_storage = &
                                                 integer_storage
    if ( present(pivotthreshold) ) solver_options%pivotthreshold = &
                                                 pivotthreshold

  end subroutine set_solver_options_ma41


! Solve the system of equations using MA41 (a direct sparse solver for
! unsymmetric matrices)

  subroutine solve_system_ma41 ( sysmatrix, rhsd, sol, lu, solver_options )

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
    type(lu_ma41_t), intent(inout), optional :: lu

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_ma41_t), intent(in), optional :: solver_options


    integer :: numundegfd
    type(lu_ma41_t) :: lw  ! local work space when lu is not present in heading
    type(solver_options_ma41_t) :: lsolveropt ! local options


!   some testing first

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(/a/)') &
        'Error in solve_system_ma41: no system matrix structure.'
      stop
    end if

    if ( .not. sysmatrix%finalized ) then
      write(*,'(/a/)') &
        'Error in solve_system_ma41: system matrix has not been finalized.'
      stop
    end if

    if ( .not. sysmatrix%allocated_data ) then
      write(*,'(/a/)') &
        'Error in solve_system_ma41: data in system matrix not allocated.'
      stop
    end if

    if ( sysmatrix%symmetric ) then
      write(*,'(/a/)') &
        'Error in solve_system_ma41: symmetric matrix not allowed.'
      stop
    end if

    if ( .not. rhsd%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_ma41: rhsd not created.'
      stop
    end if

    if ( .not. sol%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_ma41: sol not created.'
      stop
    end if

!   initialize local solver options if solver_options in heading
    if ( present(solver_options) ) lsolveropt = solver_options

    if ( lsolveropt%pivot_order == 1 .and. .not. sysmatrix%renumber ) then
      write(*,'(/a/a/)') &
        'Error in solve_system_ma41:', &
        '  pivot_order = 1, whereas no renumbering is present in sysmatrix.'
      stop
    end if

!   number of unknowns

    numundegfd = sysmatrix%Suu%n

!   copy first part of rhsd to sol

    sol%u(1:numundegfd) = rhsd%u(1:numundegfd)

!   solve, solution returned in sol(1:numundegfd)

    if ( present(lu) ) then

      call solve_ma41 ( sysmatrix%Suu, sol%u(1:numundegfd), lu, lsolveropt, &
        sysmatrix%ipermu )

    else

      call solve_ma41 ( sysmatrix%Suu, sol%u(1:numundegfd), lw, lsolveropt, &
        sysmatrix%ipermu )

!     remove temporary storage
      call delete ( lw )

    end if

  end subroutine solve_system_ma41


! solve linear system A u = b

  subroutine solve_ma41 ( A, b, lu, so, ipermu )

!   The matrix A in CSR format
    type(sparsematrix_t), intent(in) :: A

!   On entry: the right-hand side vector b.
!   On exit: the solution u = A^{-1} b.
    real(dp), intent(inout), dimension(:) :: b

!   lu is used to store the LU factorization.
!   If lu contains a LU factorization from a previous solve only a
!   backsubstition is done.  This happens even when the matrix has changed.
!   If lu has been deallocated using delete, a full factorization is performed.
    type(lu_ma41_t), intent(inout) :: lu

!   The options/parameters for the solver
    type(solver_options_ma41_t), intent(in) :: so

!   permutation array for the pivot order
    integer, dimension(:), intent(in) :: ipermu


    logical :: allsteps
    integer :: maxis, n, ne, job, maxs
    integer, dimension(20) :: info
    real(dp) :: factor
    real(dp), dimension(20) :: rinfo


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

      if ( so%pivot_order == 1 ) then
        maxis = nint ( so%integer_storage * ( ne + 12 * n + 1 ) )
      else
        maxis = nint ( so%integer_storage * ( 2 * ne + 11 * n + 1 ) )
      end if

!     allocate memory and do the analysis

      call lu_allocate_analysis_ma41 ( lu )

      call prininfo ( 'MA41AD', 15, 9, info, rinfo, so%printlevel )

      if ( info(7) > maxis ) then
        factor = info(7) * so%integer_storage / maxis
        write(*,'(/2a/a,f10.3/)') &
              'Error: integer work storage for ', &
              'unsymmetric solver (MA41) too small.', &
              'Increase integer_storage parameter to at least ', factor
        stop
      end if

      job = 5

    else

!     only solution step

      maxis = size ( lu%is )
      maxs  = size ( lu%s )

      job = 3

    end if

!   solve system

    call ma41ad ( job, n, ne, lu%irn, A%ja, A%a, b, lu%colsca, lu%rowsca, &
                  lu%keep, lu%is, maxis, lu%s, maxs, lu%cntl, lu%icntl,   &
                  info, rinfo )

    call prininfo ( 'MA41AD', 15, 9, info, rinfo, so%printlevel )

!   for next call of this routine: do solution only (no factorization) if
!   lu is not deleted

    lu%decomp = .true.

  contains


!   allocate memory for lu and do analysis step

    subroutine lu_allocate_analysis_ma41 ( lu )

      type (lu_ma41_t) :: lu


      integer, parameter :: maxsdum = 1
      integer :: row
      real(dp) :: sdum(maxsdum)


!     allocate memory for lu

      allocate ( lu%irn(ne), lu%keep(50), lu%icntl(20), lu%is(maxis) )
      allocate ( lu%cntl(10), lu%colsca(n), lu%rowsca(n) )

!     default values

      call ma41id ( lu%cntl, lu%icntl, lu%keep )

      lu%cntl(1) = so%pivotthreshold   ! threshold for pivoting
      lu%cntl(2) = so%iter_ref_tol     ! iterative refinement tolerance

      if ( so%printlevel == 2 ) then
!       do a lot of printing
        lu%icntl(2) = 6
        lu%icntl(3) = 6
        lu%icntl(4) = 3
      end if

      lu%icntl(7) = so%pivot_order

      if ( so%pivot_order == 1 ) then
!       use user supplied pivot order in ipermu
        lu%is(1:n) = ipermu
      end if

      lu%icntl(8) = so%scaling
      lu%icntl(10) = so%niter_ref   ! number of iterative refinement steps

      job = 1

!     fill row number array to transform CSR -> coordinate format

      do row = 1, n
        lu%irn( A%ia(row):A%ia(row+1)-1 ) = row
      end do

      call ma41ad ( job, n, ne, lu%irn, A%ja, A%a, b, lu%colsca, lu%rowsca, &
                    lu%keep, lu%is, maxis, sdum, maxsdum, lu%cntl, lu%icntl,&
                    info, rinfo )

      if ( so%printlevel == 3 ) then

!       print memory requirement and stop

        call prinmem

        stop

      end if

!     real storage space

      maxs = nint( so%real_storage * info(8) )

      allocate ( lu%s(maxs) )

    end subroutine lu_allocate_analysis_ma41


!   print memory requirements

    subroutine prinmem

      real(dp) :: ne, info7, info8

      ne = real(A%nnz,kind=dp)
      info7 = real(info(7),kind=dp)
      info8 = real(info(8),kind=dp)

      write(*,'(/a//3(2(a,i0),a/2(a,i0,a/)/))') &
        'Minimum storage requirements after MA41 analysis:', &
        '  Matrix: ', nint(2*ne), ' integers ', nint(ne), ' reals', &
        '    Using 32bit integers: ', &
             nint((2*ne*4+ne*8)/1000000), ' MBytes ', &
        '    Using 64bit integers: ', &
             nint((2*ne*8+ne*8)/1000000), ' MBytes ', &
        '  Factorization: ', info(7), ' integers ', info(8), ' reals', &
        '    Using 32bit integers: ', &
             nint((info7*4+info8*8)/1000000), ' MBytes ', &
        '    Using 64bit integers: ', &
             nint((info7*8+info8*8)/1000000), ' MBytes ', &
        '  Total: ', nint(2*ne+info7), ' integers ', nint(ne+info8), ' reals', &
        '    Using 32bit integers: ', &
          nint(((2*ne+info7)*4+(ne+info8)*8)/1000000), ' MBytes ', &
        '    Using 64bit integers: ', &
          nint(((2*ne+info7)*8+(ne+info8)*8)/1000000), ' MBytes '

      end subroutine prinmem

  end subroutine solve_ma41


! delete single lu

  subroutine delete_single_lu_ma41 ( lu )

    type (lu_ma41_t) :: lu

    if ( .not. allocated(lu%irn) ) then
      write(*,'(/a/)') &
        'Error delete_single_lu_ma41: arrays in lu not allocated'
      stop
    end if

    deallocate ( lu%irn, lu%keep, lu%icntl, lu%is )
    deallocate ( lu%cntl, lu%colsca, lu%rowsca )
    deallocate ( lu%s )

    lu%decomp = .false.

  end subroutine delete_single_lu_ma41


! Delete lu

  subroutine delete_lu_ma41 ( lu_ma41_1, lu_ma41_2, lu_ma41_3, &
    lu_ma41_4, lu_ma41_5 )

    type(lu_ma41_t), intent(inout) :: lu_ma41_1
    type(lu_ma41_t), intent(inout), optional :: lu_ma41_2, lu_ma41_3, &
      lu_ma41_4, lu_ma41_5

    call delete_single_lu_ma41(lu_ma41_1)
    if ( present(lu_ma41_2) ) call delete_single_lu_ma41(lu_ma41_2)
    if ( present(lu_ma41_3) ) call delete_single_lu_ma41(lu_ma41_3)
    if ( present(lu_ma41_4) ) call delete_single_lu_ma41(lu_ma41_4)
    if ( present(lu_ma41_5) ) call delete_single_lu_ma41(lu_ma41_5)

  end subroutine delete_lu_ma41


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

end module hsl_ma41_m
