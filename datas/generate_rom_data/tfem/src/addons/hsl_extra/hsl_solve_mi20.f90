
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
!   Interative solver with Algebraic Multigrid preconditioning
!
! Use is made of the MI20 preconditioner and MI2X (CG, Bi-CGSTAB, GMRES)
! iterative solvers from the HSL library (http://www.cse.scitech.ac.uk/nag/hsl/)
!
! NOTE: systems having zeros on the diagonal, like velocity-pressure systems
! and/or other constraints cannot solved using this solver interface.

module hsl_solve_mi20_m

  use glob_defs_m
  use sparse_m
  use system_defs_m
  use hsl_mi20_double
  use HSL_ZD11_double

  implicit none


! type definition for solver options

  type solver_options_mi20_t

    integer :: printlevel = 0  !  =0 no printing
                               ! >=1 print number of iterations
                               ! >=2 print residual after convergence
                               ! >=3 print resid each cycle
                               !   4 print rsave each cycle

    integer :: type_prec = 2 ! ICNTL(3) in the call of MI2X:
                             !  CG and Bi_CGSTAB:
                             !      0 == no preconditioning
                             !    >=1 == preconditioning
                             !  GMRES:
                             !      0 == no preconditioning
                             !      1 == left preconditioning
                             !      2 == right preconditioning

    integer :: maxnumiterations = -1 ! maximum number of iterations =
                                     !       ICNTL(6) in the call of MI2X:
                                     ! <=0 max number = vector length n
                                     !  >0 the maximum number of iterations

    integer :: itsolver = 5 ! The basic iterative solver:
           ! 1 CG -- Conjugate Gradient Method (MI21)
           ! 5 Bi-CGSTAB -- Bi-Conjugate Gradient Method (MI26)
           ! 8 GMRES -- Generalized Minimum RESidual method (MI24)

    integer :: mgmres = 30 !  Dimension of Krylov space for GMRES.

    real(dp) :: eps_rel = 1e-7_dp, eps_abs = 0._dp
          ! Relative tolerance and absolute tolerance, such that
          !  || residual || <= max( eps_rel * || initial residual ||, eps_abs )
          ! NOTE: eps_rel = CNTL(1) in MI2X, eps_abs = CNTL(2) in MI2X

    logical :: move_matrix = .false.  ! If move_matrix is set to .true.:
!     The HSL routine mi20 uses allocatable components to store the matrix in
!     CRS format. Tfem also uses the CRS format. Setting move_matrix=.true. will
!     move the tfem matrix (in sysmatrix) to the mi20 components to reduce
!     memory requirements.
!     The move_alloc intrinsic (fortran2003) is used for that.
!     NOTE: mi20 reorders the matrix elements.
!     NOTE: if prec is present in the solver call (see below) the user is
!       responsible for moving the matrix back to sysmatrix by supplying
!       sysmatrix as an argument to delete:
!         call delete ( prec, sysmatrix )
!       otherwise the sysmatrix becomes illegal and cannot be used or
!       even be deleted.
!     NOTE: move_matrix=.true. is incompatible with
!       use_sysmatrix_for_matvecmul = .true.

    logical :: use_sysmatrix_for_matvecmul = .false.
!     perform matrix-vector multiplication y=Az using the original sysmatrix
!     instead of the copied matrix. This might be useful if the matrix changes
!     only slightly during a time-stepping procedure and prec is present in the
!     solver call (see below) to avoid the setup of the preconditioner after
!     the second solver call (or higher). Use with care!
!     NOTE: use_sysmatrix_for_matvecmul = .true. is incompatible with
!       move_matrix=.true.

    logical :: stoponmaxnumits = .true.
!     stop (if stoponmaxnumits=.true.) with an error
!     message or continue (if stoponmaxnumits=.false.)
!     the program if the maximum number of iterations
!     (maxnumiterations) has been reached.

  end type solver_options_mi20_t


! type definition of storage for MI20

  type prec_mi20_t

    logical :: setup = .false.  ! setup done?

!   if .true. the matrix a has been created by using move_alloc to the matrix
!   in sysmatrix (which is temporarily partly unallocated).
    logical :: move_alloc_used = .false.

!   derived types for the preconditioner

    type(zd11_type) :: a
    type(mi20_data), dimension(:), allocatable :: coarse_data
    type(mi20_control) :: control
    type(mi20_info) :: info
    type(mi20_keep) :: keep
    type(ma48_control) :: ma48_cntl

  end type prec_mi20_t


! interface for generic delete subroutine

  interface delete
    module procedure delete_single_prec_mi20
  end interface delete

! interface for generic set_solver_options subroutine

  interface set_solver_options
    module procedure set_solver_options_mi20
  end interface set_solver_options

contains


! Helper routine for setting the solver options

  subroutine set_solver_options_mi20 ( solver_options, keep, move_matrix, &
    printlevel, maxnumiterations, itsolver, mgmres, eps_rel, eps_abs, &
    use_sysmatrix_for_matvecmul, stoponmaxnumits, type_prec )

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_mi20_t), intent(inout) :: solver_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see solver_options_mi20_t
    logical, intent(in), optional :: move_matrix, use_sysmatrix_for_matvecmul, &
      stoponmaxnumits
    integer, intent(in), optional :: printlevel, maxnumiterations, itsolver, &
      mgmres, type_prec
    real(dp), intent(in), optional :: eps_rel, eps_abs


    type(solver_options_mi20_t) :: solveropt


    if ( present(keep) ) then
      if ( .not. keep ) solver_options = solveropt
    else
      solver_options = solveropt
    end if

    if ( present(move_matrix) ) solver_options%move_matrix = move_matrix
    if ( present(printlevel) ) solver_options%printlevel = printlevel
    if ( present(type_prec) ) solver_options%type_prec = type_prec
    if ( present(maxnumiterations) ) solver_options%maxnumiterations = &
                                          maxnumiterations
    if ( present(itsolver) ) solver_options%itsolver = itsolver
    if ( present(mgmres) ) solver_options%mgmres = mgmres
    if ( present(eps_rel) ) solver_options%eps_rel = eps_rel
    if ( present(eps_abs) ) solver_options%eps_abs = eps_abs
    if ( present(use_sysmatrix_for_matvecmul) ) &
      solver_options%use_sysmatrix_for_matvecmul = &
                                          use_sysmatrix_for_matvecmul
    if ( present(stoponmaxnumits) ) solver_options%stoponmaxnumits = &
      stoponmaxnumits

  end subroutine set_solver_options_mi20


! Solve the system of equations using MI2X (CG, Bi-CGSTAB, GMRES) and MI20 (AMG)

  subroutine solve_system_mi20 ( sysmatrix, rhsd, sol, prec, initsol, &
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

!   If present, prec is used to store the preconditioner setup data and a copy
!   of the matrix.
!   The next time the solve routine is called with prec present no new setup
!   of the preconditioner is done and the old data is used.
!   Note, that prec contains a copy of the matrix after initial setup. This
!   copy is used in all the operations. This means that the sysmatrix in the
!   heading is basically ignored in the second (or higher) call (except if
!   solver_options%use_sysmatrix_for_matvecmul=.true.)
!   If prec is not present or has been deallocated using delete, a new setup of
!   the preconditioner is done.
!   Also if the user wants to control the preconditioner (prec%control, see
!   manual of MI20), prec must be present. Note, that in this case the user
!   must manually delete prec to destroy the setup data and the copy of the
!   matrix.
!   NOTE: if prec is present and solver_options%move_matrix is set to .true.
!     the user is responsible for moving the matrix back to the sysmatrix by
!     supplying sysmatrix as an argument to delete:
!        call delete ( prec, sysmatrix )
!     otherwise sysmatrix becomes illegal and cannot be used or even be deleted.
!   NOTE: Use with care and know what you are doing!
    type(prec_mi20_t), intent(inout), optional :: prec

!   If .true.: use initial estimate for solution vector given in sol
!   default = .false.
    logical, intent(in), optional :: initsol

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_mi20_t), intent(in), optional :: solver_options

!   This is an OUTPUT parameter indicating whether the maximum number
!   of iterations (maxnumiterations) has been reached.
!   Note, that this parameter only makes sense in combination with
!   solver_options%stoponmaxnumits=.false.
    logical, intent(out), optional :: maxnumitsreached



    logical :: initsl, setup, lmaxnumitsreached
    integer :: numundegfd
    type(prec_mi20_t) :: precw  ! local work space when prec is not present
    type(solver_options_mi20_t) :: lsolveropt


!   some testing first

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(/a/)') &
        'Error in solve_system_mi20: no system matrix structure.'
      stop
    end if

    if ( .not. sysmatrix%finalized ) then
      write(*,'(/a/)') &
        'Error in solve_system_mi20: system matrix has not been finalized.'
      stop
    end if

    if ( .not. sysmatrix%allocated_data ) then
      write(*,'(/a/)') &
        'Error in solve_system_mi20: data in system matrix not allocated.'
      stop
    end if

    if ( sysmatrix%symmetric ) then
      write(*,'(/a/)') &
        'Error in solve_system_mi20: symmetric matrix not allowed.'
      stop
    end if

    if ( .not. rhsd%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_mi20: rhsd not created.'
      stop
    end if

    if ( .not. sol%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_mi20: sol not created.'
      stop
    end if

!   initialize local solver options if solver_options in heading
    if ( present(solver_options) ) lsolveropt = solver_options

!   setup preconditioner?

    setup = .true.

    if ( present(prec) ) then
      if ( prec%setup ) then
!       setup is already done
        setup = .false.
      end if
    end if

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

!   solve, solution returned in sol(1:numundegfd)

    if ( present(prec) ) then

      call solve_mi20 ( sysmatrix%Suu, rhsd%u(1:numundegfd), &
        sol%u(1:numundegfd), prec, setup, lsolveropt, initsl, &
        lmaxnumitsreached )

!     for next call of this routine: do solution only (no setup)
      prec%setup = .true.

    else

      call solve_mi20 ( sysmatrix%Suu, rhsd%u(1:numundegfd), &
        sol%u(1:numundegfd), precw, setup, lsolveropt, initsl, &
        lmaxnumitsreached )

      call delete ( precw, sysmatrix )

    end if

    if ( present(maxnumitsreached) ) maxnumitsreached = lmaxnumitsreached

  end subroutine solve_system_mi20


! solve linear system A u = b

  subroutine solve_mi20 ( A, b, u, prec, setup, so, initsl, maxnumitsreached )

    use HSL_MC65_double

!   The matrix A in CSR format
    type(sparsematrix_t), intent(inout) :: A

!   The right-hand side vector b.
    real(dp), intent(in), dimension(:) :: b

!   The solution u = A^{-1} b.
    real(dp), intent(inout), dimension(:) :: u

!   Preconditioner data
    type(prec_mi20_t), intent(inout) :: prec

!   indicates whether new setup of preconditioner must be done
    logical, intent(in) :: setup

!   The options/parameters for the solver
    type(solver_options_mi20_t), intent(in) :: so

!   If .true.: use initial estimate for solution vector given in sol
    logical, intent(in) :: initsl

!   Indicate whether the maximum number of iterations (maxnumiterations)
!   has been reached.
    logical, intent(out) :: maxnumitsreached



    integer :: info65  ! mc65 error flag

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

    if ( setup ) then

!     setup MI20 preconditioner

!     generate local matrix

      prec%a%m = n

      if ( so%move_matrix ) then

!       move matrix

        call move_alloc ( from=A%ja, to=prec%a%col )
        call move_alloc ( from=A%a, to=prec%a%val )
        call move_alloc ( from=A%ia, to=prec%a%ptr )

        prec%move_alloc_used = .true.

      else

!       copy matrix

        allocate( prec%a%col(nnz), prec%a%val(nnz), prec%a%ptr(n+1) )

        prec%a%col = A%ja
        prec%a%val = A%a
        prec%a%ptr = A%ia

      end if

!     setup mi20

      call mi20_setup ( prec%a, prec%coarse_data, prec%keep, prec%control, &
        prec%info )

      if (prec%info%flag < 0) then
        write(*,*) "Error return from mi20_setup"
        stop
      end if

    end if

!   test matrix present

    if ( prec%move_alloc_used .and. so%use_sysmatrix_for_matvecmul ) then
      write(*,'(/a/)') &
        ' Error: matrix moved whereas use_sysmatrix_for_matvecmul=.true.'
      stop
    end if

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

        call errormsg_case_default ( 'solve_mi20', 'so%itsolver', &
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
!         convergence
          if ( so%printlevel >= 1 ) then
            write(*,'(/a,i0,a/)') &
              'Convergence after ', info(2), ' iterations.'
          end if
          if ( so%printlevel >= 2 ) then
            write(*,'(a,es16.8/)') 'L2-norm of residual = ', resid
          end if
          u = w(:,2)
          exit
        case(2)
!         perform matrix-vector multiplication y=Az
          if ( so%use_sysmatrix_for_matvecmul ) then
            w(:,locy) = smatvec ( A, w(:,locz) )
          else
            call MC65_matrix_multiply_vector ( prec%a, w(:,locz), w(:,locy), &
              info65 )
          end if
        case(3,4)
!         perform preconditioning operation y=Pz
          call mi20_precondition ( prec%a, prec%coarse_data, w(:,locz), &
            w(:,locy), prec%keep, prec%control, prec%info, prec%ma48_cntl )
           if (prec%info%flag < 0) then
             write(*,*) "Error return from mi20_precondition"
             stop
           end if
         case default
           call errormsg_case_default ( 'solve_mi20', 'iact', int_value=iact )
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

  end subroutine solve_mi20


! delete single prec

  subroutine delete_single_prec_mi20 ( prec, sysmatrix )

    type (prec_mi20_t), intent(inout) :: prec
    type (sysmatrix_t), intent(inout), optional :: sysmatrix

    if ( .not. allocated(prec%coarse_data) ) then
      write(*,'(/a/)') &
        'Error delete_single_prec_mi20: arrays in prec not allocated'
      stop
    end if

    call mi20_finalize ( prec%coarse_data, prec%keep, prec%control, prec%info )

    if ( prec%move_alloc_used ) then

!     move matrix back to sysmatrix

      if ( .not. present(sysmatrix) ) then

        write(*,'(/3(a/))') &
          'Error delete_single_prec_mi20: ', &
          '  matrix cannot be moved back to the system matrix', &
          '  sysmatrix must be added to the argument list of the routine.'
        stop

      else if ( allocated(sysmatrix%Suu%a) ) then

        write(*,'(/3(a/))') &
          'Error delete_single_prec_mi20: ', &
          '  matrix cannot be moved back to the system matrix', &
          '  sysmatrix already allocated'
        stop

      end if

      call move_alloc ( from=prec%a%col, to=sysmatrix%Suu%ja )
      call move_alloc ( from=prec%a%val, to=sysmatrix%Suu%a )
      call move_alloc ( from=prec%a%ptr, to=sysmatrix%Suu%ia )

    else

!     delete matrix

      deallocate( prec%a%col, prec%a%val, prec%a%ptr )

    end if

    prec%setup = .false.
    prec%move_alloc_used = .false.

  end subroutine delete_single_prec_mi20

end module hsl_solve_mi20_m
