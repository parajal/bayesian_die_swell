
! Copyright (C) 2022-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, routines, for solving the system of equations using the direct solver
! UMFPACK from the SuiteSparse library:
! https://people.engr.tamu.edu/davis/suitesparse.html

module umfpack_m

  use kind_defs_m
  use sparse_m, only: sparsematrix_t, sort_sparsematrix
  use system_defs_m
  use umfpack_interfaces_m

  implicit none


! size of Info and Control arrays

  integer, parameter :: UMFPACK_INFO = 90, UMFPACK_CONTROL = 20

! some index values in the Control array:

  integer, parameter :: UMFPACK_PRL = 0, UMFPACK_ORDERING = 10, &
                        UMFPACK_STRATEGY = 5, UMFPACK_PIVOT_TOLERANCE = 3, &
                        UMFPACK_SYM_PIVOT_TOLERANCE = 15, UMFPACK_SCALE = 16, &
                        UMFPACK_DROPTOL = 18, UMFPACK_IRSTEP = 7


! type definition for solver options

  type solver_options_umfpack_t

    integer :: printlevel = 1
          ! 0 no printing
          ! 1 error messages only
          ! 2 print status, Control and Info
          ! 3 print memory requirements after symbolic analysis and stop
          ! 4 print memory requirements after numerical factorization and stop

    integer :: ordering = 1  ! = Control(UMFPACK_ORDERING)
          ! 0 CHOLMOD (AMD/COLAMD then METIS)
          ! 1 AMD/COLAMD
          ! 3 METIS
          ! 4 try many orderings, pick best
          ! The following values are valid in UMFPACK, but not allowed here:
          ! 2 user-provided Qinit (not available)
          ! 5 none (not available)
          ! 6 user-provided function (not available)

    integer :: strategy = 0 ! = Control(UMFPACK_STRATEGY)
          ! 0 use sym. or unsym. strategy (default)
          ! 1 COLAMD(A), coletree postorder, not prefer diag
          ! 3 AMD(A+A'), no coletree postorder, prefer diagonal

!   pivot threshold/tolerance. See UMFPACK Userguide
    real(dp) :: &
      pivottolerance = 0.01_dp, &  ! = Control (UMFPACK_PIVOT_TOLERANCE)
      sympivottolerance = 0.001_dp ! = Control (SYM_UMFPACK_PIVOT_TOLERANCE)

    integer :: scale = 1  ! = Control(UMFPACK_SCALE)
          ! 0 no scaling
          ! 1 default: divide each row by sum (abs (row))
          ! 2 divide each row by max (abs (row))

!   drop tolerance for entries in L,U
    real(dp) :: droptol = 0.0_dp  ! = Control(UMFPACK_DROPTOL)

!   The maximum number of iterative refinement steps to attempt.
    integer :: num_iter_ref = 2  ! = Control(UMFPACK_IRSTEP)

  end type solver_options_umfpack_t


! type definition of storage for LU factorization

  type lu_umfpack_t

    logical :: decomp = .false.  ! LU decomposition done?

    real(c_double) :: Control(0:UMFPACK_CONTROL-1)

    real(c_double) :: Info(0:UMFPACK_INFO-1)

    type(c_ptr) :: Numeric

  end type lu_umfpack_t


! interface for generic delete subroutine

  interface delete
    module procedure delete_lu_umfpack
  end interface delete


! interface for generic set_solver_options subroutine

  interface set_solver_options
    module procedure set_solver_options_umfpack
  end interface set_solver_options


contains


! Helper routine for setting the solver options

  subroutine set_solver_options_umfpack ( solver_options, keep, printlevel, &
    ordering, strategy, pivottolerance, sympivottolerance, scale, droptol, &
    num_iter_ref )

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_umfpack_t), intent(inout) :: solver_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see solver_options_umfpack_t
    integer, intent(in), optional :: printlevel, ordering, strategy, &
      scale, num_iter_ref
    real(dp), intent(in), optional :: pivottolerance, sympivottolerance, &
      droptol

    type(solver_options_umfpack_t) :: solveropt


    if ( present(keep) ) then
      if ( .not. keep ) solver_options = solveropt
    else
      solver_options = solveropt
    end if

    if ( present(printlevel) ) solver_options%printlevel = printlevel
    if ( present(ordering) ) solver_options%ordering = ordering
    if ( present(strategy) ) solver_options%strategy = strategy
    if ( present(pivottolerance) ) &
                     solver_options%pivottolerance = pivottolerance
    if ( present(sympivottolerance) ) &
                     solver_options%sympivottolerance = sympivottolerance
    if ( present(scale) ) solver_options%scale = scale
    if ( present(droptol) ) solver_options%droptol = droptol
    if ( present(num_iter_ref) ) solver_options%num_iter_ref = num_iter_ref

  end subroutine set_solver_options_umfpack


! Solve the system of equations

  subroutine solve_system_umfpack ( sysmatrix, rhsd, sol, lu, solver_options )

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

!   When present, lu stores the LU decomposition.
!   The next time the solve routine is called with lu present only a
!   backsubstition is done. This happens even when the matrix has changed.
!   If lu is not present or has been deallocated using delete, a full
!   factorization is performed.
    type(lu_umfpack_t), intent(inout), optional :: lu

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_umfpack_t), intent(in), optional :: solver_options


    integer :: numundegfd
!   local work space when lu is not present in heading
    type(lu_umfpack_t) :: lw
    type(solver_options_umfpack_t) :: lsolveropt

!   some testing first

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(/a/)') &
        'Error in solve_system_umfpack: no system matrix structure.'
      stop
    end if

    if ( .not. sysmatrix%finalized ) then
      write(*,'(/a/)') &
        'Error in solve_system_umfpack: system matrix has not been finalized.'
      stop
    end if

    if ( .not. sysmatrix%allocated_data ) then
      write(*,'(/a/)') &
        'Error in solve_system_umfpack: data in system matrix not allocated.'
      stop
    end if

    if ( sysmatrix%symmetric ) then
      write(*,'(/2(a/))') &
        'Error in solve_system_umfpack: ', &
        ' the sysmatrix must not be created with symmetric=.true.'
      stop
    end if

    if ( .not. rhsd%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_umfpack: rhsd not created.'
      stop
    end if

    if ( .not. sol%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_umfpack: sol not created.'
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

      call solve_umfpack ( sysmatrix%Suu, sol%u(1:numundegfd), lu, lsolveropt )

    else

      call solve_umfpack ( sysmatrix%Suu, sol%u(1:numundegfd), lw, lsolveropt )

!     remove temporary storage
      call delete ( lw )

    end if

  end subroutine solve_system_umfpack


! solve linear system A u = b

  subroutine solve_umfpack ( A, b, lu, so )

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
    type(lu_umfpack_t), intent(inout) :: lu

!   The options/parameters for the solver
    type(solver_options_umfpack_t), intent(in) :: so


    logical :: allsteps

    integer :: status, sys

    real(c_double) :: X(A%n)
    type(c_ptr) :: Symbolic


!   initialize local parameters

    allsteps = .true.

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

!     convert matrix indices from 1-based to 0-based.

      A%ja = A%ja - 1
      A%ia = A%ia - 1

!     set default control parameters

#if I8
      call umfpack_dl_defaults ( lu%Control )
#else
      call umfpack_di_defaults ( lu%Control )
#endif

!     modify control parameters

      call modify_control

!     Optionally print Control

#if I8
      call umfpack_dl_report_control ( lu%Control )
#else
      call umfpack_di_report_control ( lu%Control )
#endif

!     symbolic analysis

#if I8
      status = umfpack_dl_symbolic ( A%n, A%n, A%ia, A%ja, A%a, Symbolic, &
        lu%Control, lu%Info )
#else
      status = umfpack_di_symbolic ( A%n, A%n, A%ia, A%ja, A%a, Symbolic, &
        lu%Control, lu%Info )
#endif

      call prinstat ( 'symbolic analysis' )

      if ( any(so%printlevel == [3,4] ) ) then
        call prinmem_symbolic
        if ( so%printlevel == 3 ) stop
      end if

!     numerical factorization

#if I8
      status = umfpack_dl_numeric ( A%ia, A%ja, A%a, Symbolic, lu%Numeric, &
        lu%Control, lu%Info )
#else
      status = umfpack_di_numeric ( A%ia, A%ja, A%a, Symbolic, lu%Numeric, &
        lu%Control, lu%Info )
#endif

      call prinstat ( 'numerical factorization' )

      if ( so%printlevel == 4 ) then
        call prinmem_numeric
        stop
      end if

!     free object Symbolic

#if I8
      call umfpack_dl_free_symbolic ( Symbolic )
#else
      call umfpack_di_free_symbolic ( Symbolic )
#endif

!     for next call of this routine: do solution only (no factorization) if
!     lu is not deleted

      lu%decomp = .true.

    else

!     convert matrix indices from 1-based to 0-based.

      A%ja = A%ja - 1
      A%ia = A%ia - 1

    end if

!   solve system (solution step only)

    sys = 2  ! solve A^T x = b to fake CSC storage

#if I8
    status = umfpack_dl_solve ( sys, A%ia, A%ja, A%a, X, b, lu%Numeric, &
      lu%Control, lu%Info )
#else
    status = umfpack_di_solve ( sys, A%ia, A%ja, A%a, X, b, lu%Numeric, &
      lu%Control, lu%Info )
#endif

    call prinstat ( 'solve' )

    b = X

!   Optionally print Info

#if I8
    call umfpack_dl_report_info ( lu%Control, lu%Info )
#else
    call umfpack_di_report_info ( lu%Control, lu%Info )
#endif

!   convert matrix indices back to 1-based

    A%ja = A%ja + 1
    A%ia = A%ia + 1

  contains

!   modify control parameters

    subroutine modify_control

!     printlevel

      if ( so%printlevel <= 2 ) then
        lu%Control(UMFPACK_PRL) = so%printlevel
      end if

!     ordering

      if ( any ( so%ordering == [0,1,3,4] ) ) then
        lu%Control(UMFPACK_ORDERING) = so%ordering
      else if ( any ( so%ordering == [2,5,6] ) ) then
        write ( *,'(/a/a,i0/)') &
          ' Error in solve_umfpack: ordering not available', &
          '  ordering = ', so%ordering
        stop
      else
        write ( *,'(/a/a,i0/)') &
          ' Error in solve_umfpack: ordering invalid', &
          '  ordering = ', so%ordering
        stop
      end if

!     strategy

      if ( any ( so%strategy == [0,1,3] ) ) then
        lu%Control(UMFPACK_STRATEGY) = so%strategy
       else
        write ( *,'(/a/a,i0/)') &
          ' Error in solve_umfpack: strategy invalid', &
          '  strategy = ', so%strategy
        stop
      end if

!     pivot threshold

      lu%Control(UMFPACK_PIVOT_TOLERANCE) = so%pivottolerance
      lu%Control(UMFPACK_SYM_PIVOT_TOLERANCE) = so%sympivottolerance

!     Scaling

      lu%Control(UMFPACK_SCALE) = so%scale

!     Drop tolerance

      lu%Control(UMFPACK_DROPTOL) = so%droptol

!     Iterative refinement

      lu%Control(UMFPACK_IRSTEP) = so%num_iter_ref

    end subroutine modify_control

!   print status and stop on error

    subroutine prinstat ( step )

      character(*), intent(in) :: step

      if ( so%printlevel == 2 .or. status /= 0 ) then
        write ( *,'(/3a,i0)') &
          'Status in solve_umfpack after ', step, ': ', status
      end if

#if I8
      call umfpack_dl_report_status ( lu%Control, status )
#else
      call umfpack_di_report_status ( lu%Control, status )
#endif

      if ( status < 0 ) then
!       Optionally print Info
#if I8
        call umfpack_dl_report_info ( lu%Control, lu%Info )
#else
        call umfpack_di_report_info ( lu%Control, lu%Info )
#endif
        stop
      end if

    end subroutine prinstat

!   print memory requirements (symbolic)

    subroutine prinmem_symbolic

      write(*,'(a/a/2(a,f0.2,a/),a,g0.2,2(/a,f0.0))') &
            'symbolic analysis:', &
            '   estimates (upper bound) for numeric LU:', &
            '   size of LU:    ', (lu%Info(20)*lu%Info(3))/2**20, ' (MB)', &
            '   memory needed: ', (lu%Info(21)*lu%Info(3))/2**20, ' (MB)', &
            '   flop count:    ', lu%Info(22), &
            '   nnz (L):       ', lu%Info(23), &
            '   nnz (U):       ', lu%Info(24)

    end subroutine prinmem_symbolic

!   print memory requirements (numeric)

    subroutine prinmem_numeric

      write(*,'(a/a/2(a,f0.2,a/),a,g0.2,2(/a,f0.0))') &
            'numerical factorization:', &
            '   actual numeric LU statistics:', &
            '   size of LU:    ', (lu%Info(40)*lu%Info(3))/2**20, ' (MB)', &
            '   memory needed: ', (lu%Info(41)*lu%Info(3))/2**20, ' (MB)', &
            '   flop count:    ', lu%Info(42), &
            '   nnz (L):       ', lu%Info(43), &
            '   nnz (U):       ', lu%Info(44)

    end subroutine prinmem_numeric

  end subroutine solve_umfpack

! delete single lu

  subroutine delete_single_lu_umfpack ( lu )

    type (lu_umfpack_t) :: lu

#if I8
    call umfpack_dl_free_numeric ( lu%Numeric )
#else
    call umfpack_di_free_numeric ( lu%Numeric )
#endif

    lu%Control = 0
    lu%Info = 0

    lu%decomp = .false.

  end subroutine delete_single_lu_umfpack


! Delete lu

  subroutine delete_lu_umfpack ( lu_umfpack_1, lu_umfpack_2, lu_umfpack_3, &
    lu_umfpack_4, lu_umfpack_5 )

    type(lu_umfpack_t), intent(inout) :: lu_umfpack_1
    type(lu_umfpack_t), intent(inout), optional :: lu_umfpack_2, lu_umfpack_3, &
      lu_umfpack_4, lu_umfpack_5

    call delete_single_lu_umfpack(lu_umfpack_1)
    if ( present(lu_umfpack_2) ) call delete_single_lu_umfpack(lu_umfpack_2)
    if ( present(lu_umfpack_3) ) call delete_single_lu_umfpack(lu_umfpack_3)
    if ( present(lu_umfpack_4) ) call delete_single_lu_umfpack(lu_umfpack_4)
    if ( present(lu_umfpack_5) ) call delete_single_lu_umfpack(lu_umfpack_5)

  end subroutine delete_lu_umfpack

end module umfpack_m
