
! Copyright (C) 2009-2009 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, routines, for solving the system of equations using
!
!   Interative solver with diagonal (Jacobi) preconditioning
!
! Use is made of MI2X (CG, Bi-CGSTAB, GMRES) iterative solvers from the
! HSL library (http://www.cse.scitech.ac.uk/nag/hsl/)
!

module hsl_solve_diag_m

  use glob_defs_m
  use sparse_m
  use system_defs_m

  implicit none


! type definition for solver options

  type solver_options_diag_t

    integer :: printlevel = 0  !  =0 no printing
                               ! >=1 print number of iterations
                               ! >=2 print residual after convergence
                               ! >=3 print resid each cycle
                               !   4 print rsave each cycle

    integer :: type_prec = 2 ! ICNTL(3) in the call of MI2X:
                             !  CG and Bi_CGSTAB:
                             !      0 == no preconditioning
                             !    >=1 == preconditioning (diagonal scaling)
                             !  GMRES:
                             !     0 == no preconditioning
                             !     1 == left preconditioning (diagonal scaling)
                             !     2 == right preconditioning (diagonal scaling)

    integer :: maxnumiterations = -1 ! maximum number of iterations =
                                     !       ICNTL(6) in the call of MI2X:
                                     ! <=0 max number = vector length n
                                     !  >0 the maximum number of iterations

    integer :: fixednumiterations = 0 ! fixed number of iterations:
                                     ! <=0 use convergence checking with
                                     !     eps_rel and eps_abs (see below)
                                     !  >0 the fixed number of iterations

    integer :: itsolver = 5 ! The basic iterative solver:
           ! 1 CG -- Conjugate Gradient Method (MI21)
           ! 5 Bi-CGSTAB -- Bi-Conjugate Gradient Method (MI26)
           ! 8 GMRES -- Generalized Minimum RESidual method (MI24)

    integer :: mgmres = 30 !  Dimension of Krylov space for GMRES.

    real(dp) :: eps_rel = 1e-7_dp, eps_abs = 0._dp
          ! Relative tolerance and absolute tolerance, such that
          !  || residual || <= max( eps_rel * || initial residual ||, eps_abs )
          ! NOTE: eps_rel = CNTL(1) in MI2X, eps_abs = CNTL(2) in MI2X

    logical :: stoponmaxnumits = .true.
!     stop (if stoponmaxnumits=.true.) with an error
!     message or continue (if stoponmaxnumits=.false.)
!     the program if the maximum number of iterations
!     (maxnumiterations) has been reached.

  end type solver_options_diag_t


! interface for generic set_solver_options subroutine

  interface set_solver_options
    module procedure set_solver_options_diag
  end interface set_solver_options

contains


! Helper routine for setting the solver options

  subroutine set_solver_options_diag ( solver_options, keep, &
    printlevel, maxnumiterations, itsolver, mgmres, eps_rel, eps_abs, &
    stoponmaxnumits, type_prec, fixednumiterations )

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_diag_t), intent(inout) :: solver_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see solver_options_diag_t
    logical, intent(in), optional :: stoponmaxnumits
    integer, intent(in), optional :: printlevel, maxnumiterations, itsolver, &
      mgmres, type_prec, fixednumiterations
    real(dp), intent(in), optional :: eps_rel, eps_abs


    type(solver_options_diag_t) :: solveropt


    if ( present(keep) ) then
      if ( .not. keep ) solver_options = solveropt
    else
      solver_options = solveropt
    end if

    if ( present(printlevel) ) solver_options%printlevel = printlevel
    if ( present(type_prec) ) solver_options%type_prec = type_prec
    if ( present(maxnumiterations) ) solver_options%maxnumiterations = &
                                          maxnumiterations
    if ( present(fixednumiterations) ) solver_options%fixednumiterations = &
                                          fixednumiterations
    if ( present(itsolver) ) solver_options%itsolver = itsolver
    if ( present(mgmres) ) solver_options%mgmres = mgmres
    if ( present(eps_rel) ) solver_options%eps_rel = eps_rel
    if ( present(eps_abs) ) solver_options%eps_abs = eps_abs
    if ( present(stoponmaxnumits) ) solver_options%stoponmaxnumits = &
      stoponmaxnumits

  end subroutine set_solver_options_diag


! Solve the system of equations using MI2X (CG, Bi-CGSTAB, GMRES) and diagonal
! scaling

  subroutine solve_system_diag ( sysmatrix, rhsd, sol, initsol, &
    solver_options, maxnumitsreached )

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

!   If .true.: use initial estimate for solution vector given in sol
!   default = .false.
    logical, intent(in), optional :: initsol

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_diag_t), intent(in), optional :: solver_options

!   This is an OUTPUT parameter indicating whether the maximum number
!   of iterations (maxnumiterations) has been reached.
!   Note, that this parameter only makes sense in combination with
!   solver_options%stoponmaxnumits=.false.
    logical, intent(out), optional :: maxnumitsreached



    logical :: initsl, lmaxnumitsreached
    integer :: numundegfd
    type(solver_options_diag_t) :: lsolveropt


!   some testing first

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(/a/)') &
        'Error in solve_system_diag: no system matrix structure.'
      stop
    end if

    if ( .not. sysmatrix%finalized ) then
      write(*,'(/a/)') &
        'Error in solve_system_diag: system matrix has not been finalized.'
      stop
    end if

    if ( .not. sysmatrix%allocated_data ) then
      write(*,'(/a/)') &
        'Error in solve_system_diag: data in system matrix not allocated.'
      stop
    end if

    if ( sysmatrix%symmetric ) then
      write(*,'(/a/)') &
        'Error in solve_system_diag: symmetric matrix not allowed.'
      stop
    end if

    if ( .not. rhsd%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_diag: rhsd not created.'
      stop
    end if

    if ( .not. sol%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_diag: sol not created.'
      stop
    end if

!   initialize local solver options if solver_options in heading
    if ( present(solver_options) ) lsolveropt = solver_options

!   initial estimate?

    if ( present(initsol) ) then
      initsl = initsol
    else
      initsl = .false.
    end if

!   number of unknowns

    numundegfd = sysmatrix%Suu%n

    if ( .not. initsl ) then
!     clear vector
      sol%u(1:numundegfd) = 0
    end if

    call solve_diag ( sysmatrix%Suu, rhsd%u(1:numundegfd), &
      sol%u(1:numundegfd), lsolveropt, initsl, lmaxnumitsreached )

    if ( present(maxnumitsreached) ) maxnumitsreached = lmaxnumitsreached

  end subroutine solve_system_diag


! solve linear system A u = b

  subroutine solve_diag ( A, b, u, so, initsl, maxnumitsreached )

!   The matrix A in CSR format
    type(sparsematrix_t), intent(inout) :: A

!   The right-hand side vector b.
    real(dp), intent(in), dimension(:) :: b

!   The solution u = A^{-1} b.
    real(dp), intent(inout), dimension(:) :: u

!   The options/parameters for the solver
    type(solver_options_diag_t), intent(in) :: so

!   If .true.: use initial estimate for solution vector given in sol
    logical, intent(in) :: initsl

!   Indicate whether the maximum number of iterations (maxnumiterations)
!   has been reached.
    logical, intent(out) :: maxnumitsreached


    logical :: error
    real(dp) :: invdiag(size(u))

!   variables for the iterative solver

    logical :: lsave(4)
    integer :: n, nnz, info(6), icntl(8), isave(17), locy, locz, iact
    integer :: m
    real(dp) :: cntl(5), resid
    real(dp), allocatable :: w(:,:), h(:,:), rsave(:)

    external mi21id, mi21ad, mi24id, mi24ad, mi26id, mi26ad


!   initialize local parameters

    info = 0

!   set dimensions of problem

    n   = A%n
    nnz = A%nnz

!   get diagonal of matrix

    invdiag = diagonal ( A, error )

    if ( error ) then
      write(*,'(/a/)') &
        ' Error solve_diag: diagonal elements missing from matrix '
      stop
    end if

    invdiag = 1 / invdiag

!   Set control parameters

    select case ( so%itsolver )
    case(1) ! CG

      allocate ( w(n,4), rsave(6) )
      call mi21id ( icntl, cntl, isave, rsave )
      icntl(7) = 1  ! use normalized curvature for testing breakdown

    case(5) ! Bi-CGSTAB

      allocate ( w(n,8), rsave(9) )
      call mi26id ( icntl, cntl, isave, rsave )

    case(8) ! GMRES

      m = so%mgmres
      allocate ( w(n,m+7), rsave(9), h(n,m+2) )
      call mi24id ( icntl, cntl, isave, rsave, lsave )

    case default

      write(*,'(/a,i0/)') &
        ' Incorrect iterative solver: itsolver = ', so%itsolver
      stop

    end select

    icntl(3) = so%type_prec  ! preconditioning
    icntl(6) = so%maxnumiterations  ! maximum number of iterations
    if ( so%fixednumiterations > 0 ) icntl(4) = 1
    cntl(1) = so%eps_rel
    cntl(2) = so%eps_abs

    if ( initsl ) then
      icntl(5) = 1  ! initial estimate
      w(:,2) = u
    end if

!   Set right-hand side

    w(:,1) = b

!   solve system

    iact = 0

    do

      select case ( so%itsolver )
      case(1) ! CG

        call mi21ad ( iact, n, w, n, locy, locz, resid, icntl, cntl, info, &
          isave, rsave )

      case(5) ! Bi-CGSTAB

        call mi26ad ( iact, n, w, n, locy, locz, resid, icntl, cntl, info, &
          isave, rsave )

      case(8) ! GMRES

        call mi24ad ( iact, n, m, w, n, locy, locz, h, n, resid, icntl, cntl, &
          info, isave, rsave, lsave )

      case default

        call errormsg_case_default ( 'solve_diag', 'so%itsolver', &
          int_value=so%itsolver )

      end select

      select case ( iact )
        case(:-1)
!         Error occured
          if ( info(1) == -4 ) then
!           maximum of iterations reached
            if ( so%stoponmaxnumits ) then
              print *, 'Iterative solver has iterated too many times.'
              stop
            else
              maxnumitsreached = .true.
              exit
            end if
          else
            stop 'Error in iterative solver'
          end if
        case(1)
          if ( so%fixednumiterations <=0 ) then
!           convergence
            if ( so%printlevel >= 1 ) then
              write(*,'(/a,i0,a/)') &
                'Convergence after ', info(2), ' iterations.'
            end if
            if ( so%printlevel >= 2 ) then
              write(*,'(a,es16.8/)') 'L2-norm of residual = ', resid
            end if
            u = w(:,2)
            exit
          else if ( so%fixednumiterations == info(2) ) then
!           number of iterations reached
            if ( so%printlevel >= 1 ) then
              write(*,'(/a,i0,a/)') &
                'Stopped after required ', info(2), ' iterations.'
            end if
            if ( so%printlevel >= 2 ) then
              write(*,'(a,es16.8/)') 'L2-norm of residual = ', resid
            end if
            u = w(:,2)
            exit
          end if
        case(2)
!         perform matrix-vector multiplication y=Az
          w(:,locy) = smatvec ( A, w(:,locz) )
        case(3,4)
!         perform preconditioning operation y=Pz
          w(:,locy) = w(:,locz) * invdiag
        case default
          call errormsg_case_default ( 'solve_diag', 'iact', int_value=iact )
      end select

      if ( so%printlevel >= 3 ) then
        write(*,'(a,i0,a,es16.8)') 'At iteration ', info(2), ' resid = ', resid
      end if

      if ( so%printlevel >= 4 ) then
        write(*,'(a,i0,a,9es16.8)') 'At iteration ', info(2), ' rsave = ', rsave
      end if

    end do

    deallocate ( w, rsave )
    if ( so%itsolver == 8 ) deallocate ( h )

  end subroutine solve_diag

end module hsl_solve_diag_m
