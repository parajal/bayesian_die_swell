
! Copyright (C) 2006-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, routines, for solving the system of equations using sparse direct
! methods in the PARDISO package.
!

module pardiso_m

  use glob_defs_m
  use sparse_m
  use system_defs_m

  implicit none


! type definition for solver options

  type solver_options_pardiso_t

    integer :: printlevel = 0  ! 0 no printing
                               ! 1,2 print solver info
                               ! 3 print memory requirements
                               ! 4 print number of iterative refinement steps
                               !   performed (iparm(7))

    integer :: matrixtype = 11 ! matrix type (mtype) in PARDISO:
                       !  mtype = 1 real and structurally symmetric matrix
                       !        = 2 real and symmetric positive definite matrix
                       !        = -2 real and symmetric indefinite matrix
                       !        = 11 real and unsymmetric matrix
                       !  NOTE: complex types mtype=3,4,-4,6,13 not available.


    integer :: ordering = 2 !iparm(2): Ordering scheme used (see manual PARDISO)
                            !  ordering=2 means metis ordering

    integer :: preconditionedcgs = 0 ! iparm(4); Preconditioned CGSC.
                                     ! See manual of PARDISO.

    integer :: niter_ref = 0 ! iparm(8); max number iterative refinement steps
                             ! See manual of PARDISO.

    logical :: setdefaultparms = .true. ! if .true. set default values for
                                        ! iparm(9:). If .false. iparm(9:) is
                                        ! taken from iparm in the heading of
                                        ! pardiso

  end type solver_options_pardiso_t


! type definition of storage for LU factorization (PARDISO)

  type lu_pardiso_t

    logical :: decomp = .false.  ! LU decomposition done?

    integer :: mtype = 0, n = 0

!   pt array (defined as 64bit integers, compatible with 32bit systems)
    integer(kind=di), dimension(64) :: pt = 0

!   iparm array
    integer, dimension(64) :: iparm = 0

  end type lu_pardiso_t


! interface for generic delete subroutine

  interface delete
    module procedure delete_lu_pardiso
  end interface delete


! interface for generic set_solver_options subroutine

  interface set_solver_options
    module procedure set_solver_options_pardiso
  end interface set_solver_options


contains


! Helper routine for setting the solver options

  subroutine set_solver_options_pardiso ( solver_options, keep, printlevel, &
    matrixtype, ordering, preconditionedcgs, niter_ref, setdefaultparms )

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_pardiso_t), intent(inout) :: solver_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see solver_options_pardiso_t
    integer, intent(in), optional :: printlevel, matrixtype, ordering, &
      preconditionedcgs, niter_ref
    logical, intent(in), optional :: setdefaultparms


    type(solver_options_pardiso_t) :: solveropt


    if ( present(keep) ) then
      if ( .not. keep ) solver_options = solveropt
    else
      solver_options = solveropt
    end if

    if ( present(printlevel) ) solver_options%printlevel = printlevel
    if ( present(matrixtype) ) solver_options%matrixtype = matrixtype
    if ( present(ordering) ) solver_options%ordering = ordering
    if ( present(preconditionedcgs) ) &
       solver_options%preconditionedcgs = preconditionedcgs
    if ( present(niter_ref) ) solver_options%niter_ref = niter_ref
    if ( present(setdefaultparms) ) &
       solver_options%setdefaultparms = setdefaultparms

  end subroutine set_solver_options_pardiso


! Solve the system of equations using PARDISO

  subroutine solve_system_pardiso ( sysmatrix, rhsd, sol, lu, solver_options, &
    iparm )

!   NOTE: the matrix is modified on output: columns will be sorted.
    type(sysmatrix_t), intent(inout) :: sysmatrix

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

!   When present, lu stores pointers to the internal solver data of PARDISO.
!   The next time the solve routine is called with lu present only a
!   backsubstition is done. This happens even when the matrix has changed.
!   If lu is not present or has been deallocated using delete, a full
!   factorization is performed.
    type(lu_pardiso_t), intent(inout), optional :: lu

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_pardiso_t), intent(in), optional :: solver_options

!   if solver_options%setdefaultparms is .false. then iparm(9:) in the call
!   of pardiso is taken from this array.
!   At output iparm contains the full iparm array after the PARDISO call.
    integer, dimension(:), intent(inout), optional :: iparm


    integer :: numundegfd
!   local work space when lu is not present in heading
    type(lu_pardiso_t) :: lw
    type(solver_options_pardiso_t) :: lsolveropt


!   some testing first

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(/a/)') &
        'Error in solve_system_pardiso: no system matrix structure.'
      stop
    end if

    if ( .not. sysmatrix%finalized ) then
      write(*,'(/a/)') &
        'Error in solve_system_pardiso: system matrix has not been finalized.'
      stop
    end if

    if ( .not. sysmatrix%allocated_data ) then
      write(*,'(/a/)') &
        'Error in solve_system_pardiso: data in system matrix not allocated.'
      stop
    end if

    if ( .not. rhsd%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_pardiso: rhsd not created.'
      stop
    end if

    if ( .not. sol%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_pardiso: sol not created.'
      stop
    end if

!   initialize local solver options if solver_options in heading
    if ( present(solver_options) ) lsolveropt = solver_options

    select case ( lsolveropt%matrixtype )
      case(1,11)
        if ( sysmatrix%symmetric ) then
          write(*,'(/3(a/))') &
            'Error in solve_system_pardiso: ', &
            ' for matrixtype=1 or 11 the sysmatrix must not be ', &
            ' created with symmetric=.true.'
          stop
        end if
      case(2,-2)
        if ( .not. sysmatrix%symmetric ) then
          write(*,'(/3(a/))') &
            'Error in solve_system_pardiso: ', &
            ' for matrixtype=2 or -2 the sysmatrix must be ', &
            ' created with symmetric=.true.'
          stop
        end if
      case default
        write(*,'(/a,i0/)') &
          'Error solve_system_pardiso: invalid matrixtype = ', &
          lsolveropt%matrixtype
        stop
    end select

    if ( .not. lsolveropt%setdefaultparms .and. .not. present(iparm) ) then
      write(*,'(/2(a/))') &
        'Error in solve_system_pardiso: ', &
        ' solver_options%setdefaultparms = .false. and iparm is not present'
      stop
    end if

    if ( present(iparm) ) then
      if ( size(iparm) /= 64 ) then
        write(*,'(/2(a/))') &
          'Error in solve_system_pardiso: ', &
          ' the size of the array iparm must be equal to 64'
        stop
      end if
    end if

!   number of unknowns

    numundegfd = sysmatrix%Suu%n

!   copy first part of rhsd to sol

    sol%u(1:numundegfd) = rhsd%u(1:numundegfd)

!   solve, solution returned in sol(1:numundegfd)

    if ( present(lu) ) then

      call solve_pardiso ( sysmatrix%Suu, sol%u(1:numundegfd), lu, lsolveropt, &
        iparm )

    else

      call solve_pardiso ( sysmatrix%Suu, sol%u(1:numundegfd), lw, lsolveropt, &
        iparm )

!     remove temporary storage
      call delete ( lw )

    end if

  end subroutine solve_system_pardiso


! solve linear system A u = b

  subroutine solve_pardiso ( A, b, lu, so, iparm )

!   The matrix A in CSR format
!   NOTE: the matrix is modified on output: colums will be sorted.

    type(sparsematrix_t), intent(inout) :: A

!   On entry: the right-hand side vector b.
!   On exit: the solution u = A^{-1} b.
    real(dp), intent(inout), dimension(:) :: b

!   lu is used to store the LU factorization.
!   If lu contains a LU factorization from a previous solve only a
!   backsubstition is done.  This happens even when the matrix has changed.
!   If lu has been deallocated using delete, a full factorization is performed.
    type(lu_pardiso_t), intent(inout) :: lu

!   The options/parameters for the solver
    type(solver_options_pardiso_t), intent(in) :: so

!   if so%setdefaultparms is .false. then iparm(9:) in the call
!   of pardiso is taken from this array.
!   At output iparm contains the full iparm array after the PARDISO call.
    integer, dimension(:), intent(inout), optional :: iparm


    logical :: allsteps

    real(dp) :: x(A%n)
    integer :: maxfct, mnum, phase, nrhs, error, perm, msglvl

    integer omp_get_max_threads
    external omp_get_max_threads


!   initialize local parameters

    allsteps = .true.
    maxfct = 1
    mnum = 1
    nrhs = 1
    select case ( so%printlevel )
      case(1,2)
        msglvl = 1
      case default
        msglvl = 0
    end select

!   do all steps (analysis, factorization, solution) or only solution?

    if ( lu%decomp ) then
!     decomposition is already done
      allsteps = .false.
    end if

!   choose job

    if ( allsteps ) then

!     all steps

!     sort matrix columns

      call sort_sparsematrix ( A )

!     set parameters

      lu%mtype = so%matrixtype

      lu%n  = A%n

      lu%iparm = 0

      lu%iparm(1) = 1  ! user input values for iparm
      lu%iparm(2) = so%ordering
      lu%iparm(3) = omp_get_max_threads()  ! number of procs = OMP_NUM_THREADS
      lu%iparm(4) = so%preconditionedcgs
      lu%iparm(6) = 1  ! store solution in b
      lu%iparm(8) = so%niter_ref ! iterative refinement

      if ( so%setdefaultparms ) then
!       set default values for iparm(9:)
!       pivot perturbation
        if ( so%matrixtype == 11 ) then
          lu%iparm(10) = 13
        else if ( so%matrixtype == -2 ) then
          lu%iparm(10) = 8
        end if
!       scaling vectors and improved accuracy using weighted matching
        if ( so%matrixtype == 11 ) then
          lu%iparm(11) = 1
          lu%iparm(13) = 1
        end if
        lu%iparm(18) = -1
      else
!       copy
        lu%iparm(9:) = iparm(9:)
      end if

!     allocate memory and do the analysis

      phase = 11 ! only reordering and symbolic factorization

      call pardiso ( lu%pt, maxfct, mnum, lu%mtype, phase, lu%n, A%a, A%ia, &
        A%ja, perm, nrhs, lu%iparm, msglvl, b, x, error )

      if ( error /= 0 ) then
        write ( *,'(/a/a,i0/)') &
          ' Error in solve_pardiso after analysis: ', &
          '  error = ', error
        stop
      end if

      phase = 22 ! numerical factorization

      call pardiso ( lu%pt, maxfct, mnum, lu%mtype, phase, lu%n, A%a, A%ia, &
        A%ja, perm, nrhs, lu%iparm, msglvl, b, x, error )

      if ( error /= 0 ) then
        write ( *,'(/a/a,i0/)') &
          ' Error in solve_pardiso after factorization: ', &
          '  error = ', error
        stop
      end if

      if ( so%printlevel == 3 ) then
        call prinmem
      end if

!     for next call of this routine: do solution only (no factorization) if
!     lu is not deleted

      lu%decomp = .true.

    end if

!   solve system (solution step only)

    phase = 33

    call pardiso ( lu%pt, maxfct, mnum, lu%mtype, phase, lu%n, A%a, A%ia, &
      A%ja, perm, nrhs, lu%iparm, msglvl, b, x, error )

    if ( present(iparm) ) iparm = lu%iparm

    if ( error /= 0 ) then
      write ( *,'(/a/a,i0/)') &
        ' Error in solve_pardiso after solve: ', &
        '  error = ', error
      stop
    end if

    if ( so%printlevel == 4 ) then
      write ( *,'(/a,i0/)') &
        ' Number of interative refinement steps in solve_pardiso = ', &
        lu%iparm(7)
    end if

  contains

!   print memory requirements

    subroutine prinmem

      integer :: peak
      real(dp) :: ne

      ne = real(A%nnz,kind=dp)
      peak = max(lu%iparm(15),lu%iparm(16)+lu%iparm(17))/1000


      write(*,'(/a//2(a,i0),a/2(a,i0,a/)/a,i0,a/a/a,i0,a/a,i0,a//)') &
        'Minimum storage requirements PARDISO:', &
        '  Matrix: ', nint(ne), ' integers ', nint(ne), ' reals', &
        '    Using 32bit integers: ', &
             nint((ne*4+ne*8)/1000000), ' MBytes ', &
        '    Using 64bit integers: ', &
             nint((ne*8+ne*8)/1000000), ' MBytes ', &
        '  Factorization: ', &
             peak, ' MBytes ', &
        '  Total: ', &
        '    Using 32bit integers: ', &
           nint((ne*4+ne*8)/1000000)+peak, ' MBytes ', &
        '    Using 64bit integers: ', &
           nint((ne*8+ne*8)/1000000)+peak, ' MBytes '

      end subroutine prinmem

  end subroutine solve_pardiso



! delete single lu

  subroutine delete_single_lu_pardiso ( lu )

    type (lu_pardiso_t) :: lu

    real(dp) :: ddum(1)
    integer :: maxfct, mnum, phase, nrhs, error, perm, msglvl, idum(1)

    phase = -1
    msglvl = 0

    call pardiso ( lu%pt, maxfct, mnum, lu%mtype, phase, lu%n, ddum, idum, &
      idum, perm, nrhs, lu%iparm, msglvl, ddum, ddum, error )

    lu%pt = 0
    lu%decomp = .false.
    lu%mtype = 0
    lu%iparm = 0
    lu%n = 0

  end subroutine delete_single_lu_pardiso


! Delete lu

  subroutine delete_lu_pardiso ( lu_pardiso_1, lu_pardiso_2, lu_pardiso_3, &
    lu_pardiso_4, lu_pardiso_5 )

    type(lu_pardiso_t), intent(inout) :: lu_pardiso_1
    type(lu_pardiso_t), intent(inout), optional :: lu_pardiso_2, lu_pardiso_3, &
      lu_pardiso_4, lu_pardiso_5

    call delete_single_lu_pardiso(lu_pardiso_1)
    if ( present(lu_pardiso_2) ) call delete_single_lu_pardiso(lu_pardiso_2)
    if ( present(lu_pardiso_3) ) call delete_single_lu_pardiso(lu_pardiso_3)
    if ( present(lu_pardiso_4) ) call delete_single_lu_pardiso(lu_pardiso_4)
    if ( present(lu_pardiso_5) ) call delete_single_lu_pardiso(lu_pardiso_5)

  end subroutine delete_lu_pardiso


end module pardiso_m
