
! Copyright (C) 2023-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, routines, for solving eigenvalue problems using the FEAST library
! package (https://www.feast-solver.org).

module feast_m

  use kind_defs_m
  use sparse_m, only: sparsematrix_t, sort_sparsematrix
  use system_defs_m

  implicit none


! type definition for solver options

  type solver_options_feast_t

    integer :: printlevel = 0  ! 0 no printing
                               ! 1 print runtime data to standard output
                               ! n<0 Write/Append comments in the file
                               !     feast<|n|>.log
                               ! Set to fpm(1)

    integer :: stop_conv = 12 ! Stopping convergence criteria in double
                              ! precision (0 to 16) epsilon = 10^{-convstop}.
                              ! Set to fpm(3)

    integer :: max_loop = 20 ! Maximum number of FEAST refinement loops allowed.
                             ! Set to fpm(4)

    integer :: conv_crit = 1 ! Convergence criteria (for solutions in the
                             ! search contour):
                             ! 0: Using relative error on the trace epsout
                             !    i.e.  epsout < epsilon
                             ! 1: Using relative residual res
                             !    i.e. maxi res(i) < epsilon
                             ! Set to fpm(6)

    integer :: only_right_eigv = -1 ! If >= 0 then:
      ! #Contours for non-Hermitian or polynomial FEAST.
      !   0: two-sided contour (compute right/left eigenvectors)
      !   1: one-sided contour (compute only right eigenvectors)
      !   2: one sided contour (left=right* eigenvectors)
      ! Set to fpm(15)

!   Options for changing the default contour ellipse
    real(dp) :: ratio = -1._dp
      ! if ratio > 0 then:
      !   Ellipse contour ratio ’vertical axis’/’horizontal axis’
      !   Set fpm(18) to nint( ratio * 100 )
    integer :: angle = 0
      ! Ellipse rotation angle in degree from vertical axis [-180:180]
      ! Set to fpm(19)

    integer :: scaling = 1 ! Matrix scaling for sparse drivers (0: No; 1: Yes).
                           ! Set to fpm(41)

    integer :: mixed_precision = 1 ! Mixed Precision for all drivers
         ! 0: use double precision linear system solvers
         ! 1: use single precision linear system solvers
         ! Set to fpm(42)

  end type solver_options_feast_t


! type definition of in/out parameters for call of FEAST routines

  type arg_feast_t

    logical :: filled = .false.  ! has data been filled by a call

!   FEAST input parameters (see Table 1 of documentation)
    integer, dimension(64) :: fpm = 0

!   Number of FEAST subspace iterations (out)
    integer :: loop = 0

!   Search subspace dimension (in/out)
!   On entry: initial guess M0 >= M
!   On exit: new suitable M0 if guess too large
    integer :: M0 = 1

!   #Eigenvalues found inside contour (out)
    integer :: M = 0

!   Error handling (out)
!   See FEAST documentation for the info codes.
    integer :: info = 0

!   Contour ellipse: eigenvalues are searched within the contour (in)
!   Emid and r are the centroid and horizontal radius of the contour ellipse
!   in the complex plane.
    complex(dp) :: Emid = (0._dp,0._dp)
    real(dp) :: r = 1._dp

!   Trace relative error |tracek − tracek−1|/ max(|Emid| + r) (out)
    real(dp) :: epsout = 0._dp

!   Eigenvalues (in/out)
!   On entry: initial guess if fpm(5)=1 (previous FEAST run)
!   On exit: Eigenvalues solutions E(1:M)
!   Remark: the E(M+1:M0) values are outside the contour
    complex(dp), dimension(:), allocatable :: E

!   Eigenvectors (N: size of the system) (in/out)
!   On entry: initial guess if fpm(5)=1 (previous FEAST run)
!   On exit: (right) Eigenvectors solutions X(1:N,1:M)
!   Remarks: -left vectors (if calculated) in X(1:N,M0+1:M0+M)
!            -if fpm(14)=1, first Q subspace on exit
    complex(dp), dimension(:,:), allocatable :: X

!   Relative residual res(1:M) (right); res(M0+1,M0+M) (left) (out)
!   See FEAST documentation for the definition.
    real(dp), dimension(:), allocatable :: res

  end type arg_feast_t


! interface for generic delete subroutine

  interface delete
    module procedure delete_arg_feast
  end interface delete

! interface for generic set_solver_options subroutine

  interface set_solver_options
    module procedure set_solver_options_feast
  end interface set_solver_options


contains


! Helper routine for setting the solver options

  subroutine set_solver_options_feast ( solver_options, keep, printlevel, &
    stop_conv, max_loop, conv_crit, only_right_eigv, ratio, angle, scaling, &
    mixed_precision )

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_feast_t), intent(inout) :: solver_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see solver_options_feast_t
    integer, intent(in), optional :: printlevel, stop_conv, max_loop, &
      conv_crit, only_right_eigv, angle, scaling, mixed_precision
    real(dp), intent(in), optional :: ratio

    type(solver_options_feast_t) :: solveropt

    if ( present(keep) ) then
      if ( .not. keep ) solver_options = solveropt
    else
      solver_options = solveropt
    end if

    if ( present(printlevel) ) solver_options%printlevel = printlevel
    if ( present(stop_conv) ) solver_options%stop_conv = stop_conv
    if ( present(max_loop) ) solver_options%max_loop = max_loop
    if ( present(conv_crit) ) solver_options%conv_crit = conv_crit
    if ( present(only_right_eigv) ) &
             solver_options%only_right_eigv = only_right_eigv
    if ( present(ratio) ) solver_options%ratio = ratio
    if ( present(angle) ) solver_options%angle = angle
    if ( present(scaling) ) solver_options%scaling = scaling
    if ( present(mixed_precision) ) &
             solver_options%mixed_precision = mixed_precision

  end subroutine set_solver_options_feast


! Solve the real general eigenvalue system (A-lambda*B)v=0 using FEAST

  subroutine solve_system_real_gen_feast ( sysmatrixA, sysmatrixB, arg, &
    solver_options )

!   system matrices A and B
!   NOTE: the matrices are modified on output: columns will be sorted.
    type(sysmatrix_t), intent(inout) :: sysmatrixA, sysmatrixB

!   Arg contains the arguments to the FEAST routine. The components of arg
!   that are "in" should have a proper value.
    type(arg_feast_t), intent(inout) :: arg

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_feast_t), intent(in), optional :: solver_options


    type(solver_options_feast_t) :: lsolveropt ! local options


!   some testing first

    if ( any ( .not. [ sysmatrixA%initialized_structure, &
                       sysmatrixB%initialized_structure ] ) ) then
      write(*,'(/a/)') &
        'Error in solve_system_feast: no system matrix structure.'
      stop
    end if

    if ( any ( .not. [ sysmatrixA%finalized, sysmatrixB%finalized ] ) ) then
      write(*,'(/a/)') &
        'Error in solve_system_feast: system matrix has not been finalized.'
      stop
    end if

    if ( any ( .not. [ sysmatrixA%allocated_data, &
                       sysmatrixB%allocated_data ] ) ) then
      write(*,'(/a/)') &
        'Error in solve_system_feast: data in system matrix not allocated.'
      stop
    end if

    if ( any ( [ sysmatrixA%symmetric, sysmatrixB%symmetric ] ) ) then
      write(*,'(/a/)') &
        'Error in solve_system_feast: symmetric matrix not allowed.'
      stop
    end if

!   initialize local solver options if solver_options in heading
    if ( present(solver_options) ) lsolveropt = solver_options

!   solve eigen system

    call solve_real_gen_feast ( sysmatrixA%Suu, sysmatrixB%Suu, arg, &
      lsolveropt )

  end subroutine solve_system_real_gen_feast


! solve real general eigen system (A -lambda B) v = 0

  subroutine solve_real_gen_feast ( A, B, arg, so )

!   The matrices A and B in CSR format
    type(sparsematrix_t), intent(inout) :: A, B

!   Arg contains the arguments to the FEAST routine. The components of arg
!   that are "in" should have a proper value.
    type(arg_feast_t), intent(inout) :: arg

!   The options/parameters for the solver
    type(solver_options_feast_t), intent(in) :: so

    integer :: n

!   initialize fpm

    call feastinit ( arg%fpm )

!   set fpm

    arg%fpm(1) = so%printlevel
    arg%fpm(3) = so%stop_conv
    arg%fpm(4) = so%max_loop
    arg%fpm(6) = so%conv_crit
    if ( so%only_right_eigv >= 0 ) arg%fpm(15) = so%only_right_eigv
    if ( so%ratio > 0 ) arg%fpm(18) = nint( so%ratio * 100 )
    arg%fpm(19) = so%angle
    arg%fpm(41) = so%scaling
    arg%fpm(42) = so%mixed_precision

!   set dimension of problem

    n  = A%n

!   allocate space

    if ( .not. arg%filled ) then

      allocate ( arg%E(arg%M0) )
      if ( arg%fpm(15) == 0 ) then
        allocate ( arg%X(n,2*arg%M0), arg%res(2*arg%M0) )
      else
        allocate ( arg%X(n,arg%M0), arg%res(arg%M0) )
      end if

    end if

!   sort matrix columns for pardiso

    call sort_sparsematrix ( A )
    call sort_sparsematrix ( B )

!   call feast

    call dfeast_gcsrgv ( n, A%a, A%ia, A%ja, B%a, B%ia, B%ja, arg%fpm, &
      arg%epsout, arg%loop, arg%Emid, arg%r, arg%M0, arg%E, arg%X, arg%M, &
      arg%res, arg%info )

!   for next call set filled to be true

    arg%filled = .true.

  end subroutine solve_real_gen_feast


! delete single lu

  subroutine delete_single_arg_feast ( arg )

    type (arg_feast_t) :: arg

    if ( .not. allocated(arg%E) ) then
      write(*,'(/a/)') &
        'Error delete_single_arg_feast: arrays in arg not allocated'
      stop
    end if

    call deall ( arg )

  contains

    subroutine deall ( arg )
      type(arg_feast_t), intent(out) :: arg
    end subroutine deall

  end subroutine delete_single_arg_feast


! Delete arg

  subroutine delete_arg_feast ( arg_1, arg_2, arg_3, arg_4, arg_5 )

    type(arg_feast_t), intent(inout) :: arg_1
    type(arg_feast_t), intent(inout), optional :: arg_2, arg_3, arg_4, arg_5

    call delete_single_arg_feast(arg_1)
    if ( present(arg_2) ) call delete_single_arg_feast(arg_2)
    if ( present(arg_3) ) call delete_single_arg_feast(arg_3)
    if ( present(arg_4) ) call delete_single_arg_feast(arg_4)
    if ( present(arg_5) ) call delete_single_arg_feast(arg_5)

  end subroutine delete_arg_feast

end module feast_m
