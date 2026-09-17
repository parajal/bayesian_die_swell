
! Copyright (C) 2009-2009 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, routines, for solving the system of equations of the following form
!
!     [ F  B^T ] [v] = [f]
!     [ B  0   ] [p]   [g]
!
! using a preconditioned iterative solver. The preconditioner is the inverse
! of the following form:
!
!     [ F  B^T ]
!     [ 0  -S  ]
!
! where S=BF^{-1}B^T is the Schur complement. The following assumptions are
! made:
!
!   1) The action of F^{-1} can be approximated by one or more V-cycles
!      of Algebraic Multigrid preconditioning.
!   2) The action of S^{-1} can be approximated by
!
!            M_p^{-1} F_p A_p^{-1}                               (2)
!
!      where the action of M_p^{-1} can be approximated by a fixed number of
!      iterations of a diagonally scaled CG iterative method and the action
!      of A_p^{-1} can be approximated by one or more V-cycles of Algebraic
!      Multigrid preconditioning.
!
! Use is made of the MI20 preconditioner and MI2X (CG, Bi-CGSTAB, GMRES)
! iterative solvers from the HSL library (http://www.cse.scitech.ac.uk/nag/hsl/)
!
! Background:
!
!   A good approximation for the action of S^{-1} for Navier-Stokes
!   problems is [1]
!
!            M_p^{-1} F_p A_p^{-1}                               (2)
!
!   where M_p is the pressure mass matrix, A_p is the pressure diffusion matrix
!   and F_p is a pressure matrix that resembles the Navier-Stokes structure
!   for a scalar. If convection is taken implicit, F_p will be of the form [2]:
!
!       F_p = f M_p + N(p) + eta A_p
!
!   where the factor f depends on the time-integration scheme, but scales with
!   rho/deltat (rho=density), N(p) is a pressure convection term and eta is the
!   viscosity.
!
!   [1] H.C. Elman, D.J. Silvester, A.J. Wathen: "Finite Elements and
!       Fast Iterative Solvers", Oxford University Press, 2005.
!   [2] D.A. Kay, P.M. Gresho, D.F. Griffiths and D.J. Silvester:
!       "Adaptive time-stepping for incompressible flow Part II: Navier-Stokes
!       Equations", MIMS EPrint: 2008.61, The University of Manchester.
!

module hsl_solve3_mi20_m

  use glob_defs_m
  use sparse_m
  use system_defs_m
  use hsl_mi20_double
  use HSL_ZD11_double

  implicit none


! type definition for solver options

  type solver_options3_mi20_t

    integer :: printlevel = 0  !  =0 no printing
                               ! >=1 print number of iterations
                               ! >=2 print residual after convergence
                               ! >=3 print resid each cycle
                               !   4 print rsave each cycle

    integer :: type_prec = 2 ! ICNTL(3) in the call of MI2X:
                             !  Bi_CGSTAB:
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

    integer :: cg_printlevel = 0  ! Printlevel in the CG iterative solver
!     for the S^{-1} action of the preconditioner.

    integer :: cg_fixednumits = 5 ! Fixed number of iterations in the CG
!     iterative solver for the S^{-1} action of the preconditioner.

    logical :: move_matrix = .false.  ! If move_matrix is set to .true.:
!     The HSL routine mi20 uses allocatable components to store the matrix in
!     CRS format. Tfem also uses the CRS format. Setting move_matrix=.true. will
!     move the tfem sysmatrix Ap to the mi20 components to reduce
!     memory requirements.
!     The move_alloc intrinsic (fortran2003) is used for that.
!     NOTE: mi20 reorders the matrix elements.
!     NOTE: if prec is present in the solver call (see below) the user is
!       responsible for moving the matrix back to sysmatrix by supplying
!       sysmatrix as an argument to delete:
!         call delete ( prec, Ap )
!       otherwise the sysmatrix becomes illegal and cannot be used or
!       even be deleted.

  end type solver_options3_mi20_t


! type definition of storage for MI20

  type prec3_mi20_t

    logical :: setup = .false.  ! setup done?

!   if .true. the matrix Ap has been created by using move_alloc to the matrix
!   in sparsematrix Ap (which is temporarily unallocated).
    logical :: move_alloc_used = .false.

!   derived types for the preconditioner

    type(zd11_type) :: F
    type(sparsematrix_t) :: BT
    type(zd11_type) :: Ap

    type(mi20_data), dimension(:), allocatable :: coarse_data_F
    type(mi20_control) :: control_F
    type(mi20_info) :: info_F
    type(mi20_keep) :: keep_F
    type(ma48_control) :: ma48_cntl

    type(mi20_data), dimension(:), allocatable :: coarse_data_Ap
    type(mi20_control) :: control_Ap
    type(mi20_info) :: info_Ap
    type(mi20_keep) :: keep_Ap

  end type prec3_mi20_t


! interface for generic delete subroutine

  interface delete
    module procedure delete_single_prec3_mi20
  end interface delete

! interface for generic set_solver_options subroutine

  interface set_solver_options
    module procedure set_solver_options3_mi20
  end interface set_solver_options

contains


! Helper routine for setting the solver options

  subroutine set_solver_options3_mi20 ( solver_options, keep, &
    printlevel, maxnumiterations, itsolver, mgmres, eps_rel, eps_abs, &
    stoponmaxnumits, type_prec, cg_printlevel, cg_fixednumits, move_matrix )

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options3_mi20_t), intent(inout) :: solver_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see solver_options_mi20_t
    logical, intent(in), optional :: stoponmaxnumits, move_matrix
    integer, intent(in), optional :: printlevel, maxnumiterations, itsolver, &
      mgmres, type_prec, cg_printlevel, cg_fixednumits
    real(dp), intent(in), optional :: eps_rel, eps_abs


    type(solver_options3_mi20_t) :: solveropt


    if ( present(keep) ) then
      if ( .not. keep ) solver_options = solveropt
    else
      solver_options = solveropt
    end if

    if ( present(printlevel) ) solver_options%printlevel = printlevel
    if ( present(type_prec) ) solver_options%type_prec = type_prec
    if ( present(maxnumiterations) ) solver_options%maxnumiterations = &
                                          maxnumiterations
    if ( present(itsolver) ) solver_options%itsolver = itsolver
    if ( present(mgmres) ) solver_options%mgmres = mgmres
    if ( present(eps_rel) ) solver_options%eps_rel = eps_rel
    if ( present(eps_abs) ) solver_options%eps_abs = eps_abs
    if ( present(stoponmaxnumits) ) solver_options%stoponmaxnumits = &
      stoponmaxnumits
    if ( present(cg_printlevel) ) solver_options%cg_printlevel = &
      cg_printlevel
    if ( present(cg_fixednumits) ) solver_options%cg_fixednumits = &
      cg_fixednumits
    if ( present(move_matrix) ) solver_options%move_matrix = move_matrix

  end subroutine set_solver_options3_mi20


! Solve the system of equations using MI2X (Bi-CGSTAB, GMRES) and MI20 (AMG)
! and CG as preconditioners for the velocity and pressure part, respectively

  subroutine solve_system3_mi20 ( sysmatrix, rhsd, sol, ssv, ssp, &
    sysmatrix_Mp, sysmatrix_Fp, sysmatrix_Ap, prec, initsol, solver_options, &
    maxnumitsreached )

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

!   The vector subscripts for the "v" and the "p" part, respectively.
    integer, intent(in), dimension(:) :: ssv, ssp

!   The matrices Mp, Fp and Ap in sysmatrix format
    type(sysmatrix_t), intent(inout) :: sysmatrix_Mp, sysmatrix_Fp, sysmatrix_Ap

!   If present, prec is used to store:
!    1) the preconditioner setup data for the matrix F and a copy of the matrix
!    2) a copy of the matrix B^T
!    3) the preconditioner setup data for the matrix Ap and a copy of the matrix
!   The next time the solve routine is called with prec present no new setup
!   of the preconditioner is done and the old data is used.
!   Note, that prec contains a copy of the matrix F after initial setup. This
!   copy is used in all the operations. This means that the sysmatrix in the
!   heading is basically ignored in the second (or higher) call.
!   If prec is not present or has been deallocated using delete, a new setup of
!   the preconditioner is done.
!   Also if the user wants to control the preconditioner (prec%control_F and
!   prec%control_Ap, see manual of MI20), prec must be present.
!   Note, that in this case the user must manually delete prec to destroy
!   the setup data and the copy of the matrices.
!   NOTE: if prec is present and solver_options%move_matrix is set to .true.
!     the user is responsible for moving the matrix Ap back to the sysmatrix by
!     supplying sysmatrix_Ap as an argument to delete:
!        call delete ( prec, sysmatrix_Ap )
!     otherwise sysmatrix_Ap becomes illegal and cannot be used or even be
!     deleted.
!   NOTE: Use with care and know what you are doing!
    type(prec3_mi20_t), intent(inout), optional :: prec

!   If .true.: use initial estimate for solution vector given in sol
!   default = .false.
    logical, intent(in), optional :: initsol

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options3_mi20_t), intent(in), optional :: solver_options

!   This is an OUTPUT parameter indicating whether the maximum number
!   of iterations (maxnumiterations) has been reached.
!   Note, that this parameter only makes sense in combination with
!   solver_options%stoponmaxnumits=.false.
    logical, intent(out), optional :: maxnumitsreached



    logical :: initsl, setup, lmaxnumitsreached
    integer :: numundegfd
    type(prec3_mi20_t) :: precw  ! local work space when prec is not present
    type(solver_options3_mi20_t) :: lsolveropt


!   some testing first

    call check_sysm ( sysmatrix )
    call check_sysm ( sysmatrix_Mp )
    call check_sysm ( sysmatrix_Fp )
    call check_sysm ( sysmatrix_Ap )

    if ( .not. rhsd%created ) then
      write(*,'(/a/)') &
        'Error in solve_system3_mi20: rhsd not created.'
      stop
    end if

    if ( .not. sol%created ) then
      write(*,'(/a/)') &
        'Error in solve_system3_mi20: sol not created.'
      stop
    end if

    if ( size(ssv) + size(ssp) /= sysmatrix%Suu%n ) then
      write(*,'(/a/)') &
        'Error in solve_system1_mi20: ', &
        ' number of unknowns in sysmatrix /= size(ssv)+size(ssp)'
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

      call solve3_mi20 ( sysmatrix%Suu, rhsd%u(1:numundegfd), &
        sol%u(1:numundegfd), ssv, ssp, sysmatrix_Mp%Suu, sysmatrix_Fp%Suu, &
        sysmatrix_Ap%Suu, prec, setup, lsolveropt, initsl, lmaxnumitsreached )

!     for next call of this routine: do solution only (no setup)
      prec%setup = .true.

    else

      call solve3_mi20 ( sysmatrix%Suu, rhsd%u(1:numundegfd), &
        sol%u(1:numundegfd), ssv, ssp, sysmatrix_Mp%Suu, sysmatrix_Fp%Suu, &
        sysmatrix_Ap%Suu, precw, setup, lsolveropt, initsl, lmaxnumitsreached )

      call delete ( precw, sysmatrix_Ap )

    end if

    if ( present(maxnumitsreached) ) maxnumitsreached = lmaxnumitsreached

  contains

    subroutine check_sysm ( sysmatrix )

      type(sysmatrix_t), intent(inout) :: sysmatrix

      if ( .not. sysmatrix%initialized_structure ) then
        write(*,'(/a/)') &
          'Error in solve_system3_mi20: no system matrix structure.'
        stop
      end if

      if ( .not. sysmatrix%finalized ) then
        write(*,'(/a/)') &
          'Error in solve_system3_mi20: system matrix has not been finalized.'
        stop
      end if

      if ( .not. sysmatrix%allocated_data ) then
        write(*,'(/a/)') &
          'Error in solve_system3_mi20: data in system matrix not allocated.'
        stop
      end if

      if ( sysmatrix%symmetric ) then
        write(*,'(/a/)') &
          'Error in solve_system3_mi20: symmetric matrix not allowed.'
        stop
      end if

    end subroutine check_sysm

  end subroutine solve_system3_mi20


! solve linear system A u = b

  subroutine solve3_mi20 ( A, b, u, ssv, ssp, Mp, Fp, Ap, prec, setup, so, &
    initsl, maxnumitsreached )

    use hsl_solve_diag_m
    use HSL_MC65_double

!   The full system matrix A in CSR format
    type(sparsematrix_t), intent(inout) :: A

!   The right-hand side vector b.
    real(dp), intent(in), dimension(:) :: b

!   The solution u = A^{-1} b.
    real(dp), intent(inout), dimension(:) :: u

!   The vector subscripts for the "v" and the "p" part, respectively.
    integer, intent(in), dimension(:) :: ssv, ssp

!   The matrices Mp, Fp and Ap in CSR format
    type(sparsematrix_t), intent(inout) :: Mp, Fp, Ap

!   Preconditioner data
    type(prec3_mi20_t), intent(inout) :: prec

!   indicates whether new setup of preconditioner must be done
    logical, intent(in) :: setup

!   The options/parameters for the solver
    type(solver_options3_mi20_t), intent(in) :: so

!   If .true.: use initial estimate for solution vector given in sol
    logical, intent(in) :: initsl

!   Indicate whether the maximum number of iterations (maxnumiterations)
!   has been reached.
    logical, intent(out) :: maxnumitsreached


!   variables for the iterative solver

    logical :: lsave(4)
    integer :: n, nnz, info(6), icntl(8), isave(17), locy, locz, iact
    integer :: m
    real(dp) :: cntl(5), resid
    real(dp), allocatable :: w(:,:), h(:,:), rsave(:)

    logical :: ldummy
    type(sparsematrix_t) :: F
    type(solver_options_diag_t) :: so_diag
    real(dp), allocatable, dimension(:) :: wv, wv1, wp1, wp2

    external mi21id, mi21ad, mi24id, mi24ad, mi26id, mi26ad

    allocate ( wv(size(ssv)), wv1(size(ssv)), wp1(size(ssp)), wp2(size(ssp)) )

!   initialize local parameters

    info = 0

!   set dimensions of problem

    n   = A%n
    nnz = A%nnz

    if ( setup ) then

!     setup MI20 preconditioner

!     generate local matrices F and B^T

      prec%F%m = size(ssv)

      F = extract_submat ( A, ssv, ssv )

!     move matrix

      call move_alloc ( from=F%ja, to=prec%F%col )
      call move_alloc ( from=F%a, to=prec%F%val )
      call move_alloc ( from=F%ia, to=prec%F%ptr )

      prec%BT = extract_submat ( A, ssv, ssp )

!     setup mi20

      call mi20_setup ( prec%F, prec%coarse_data_F, prec%keep_F, &
        prec%control_F, prec%info_F )

      if (prec%info_F%flag < 0) then
        write(*,*) "Error return from mi20_setup for matrix F"
        stop
      end if

!     generate local matrix Ap

      prec%Ap%m = Ap%m

      if ( so%move_matrix ) then

!       move matrix

        call move_alloc ( from=Ap%ja, to=prec%Ap%col )
        call move_alloc ( from=Ap%a, to=prec%Ap%val )
        call move_alloc ( from=Ap%ia, to=prec%Ap%ptr )

        prec%move_alloc_used = .true.

      else

!       copy matrix

        allocate( prec%Ap%col(Ap%nnz), prec%Ap%val(Ap%nnz), &
                  prec%Ap%ptr(Ap%n+1) )

        prec%Ap%col = Ap%ja
        prec%Ap%val = Ap%a
        prec%Ap%ptr = Ap%ia

      end if

!     setup mi20

      call mi20_setup ( prec%Ap, prec%coarse_data_Ap, prec%keep_Ap,&
        prec%control_Ap, prec%info_Ap )

      if (prec%info_Ap%flag < 0) then
        write(*,*) "Error return from mi20_setup for matrix Ap"
        stop
      end if

    end if

!   Set control parameters

    select case ( so%itsolver )
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

!   Set CG solve for matrix S

    so_diag%printlevel = so%cg_printlevel
    so_diag%fixednumiterations = so%cg_fixednumits
    so_diag%itsolver = 1 ! CG

!   solve system

    iact = 0

    do

      select case ( so%itsolver )
      case(5) ! Bi-CGSTAB

        call mi26ad ( iact, n, w, n, locy, locz, resid, icntl, cntl, info, &
          isave, rsave )

      case(8) ! GMRES

        call mi24ad ( iact, n, m, w, n, locy, locz, h, n, resid, icntl, cntl, &
          info, isave, rsave, lsave )

      case default

        call errormsg_case_default ( 'solve3_mi20', 'so%itsolver', &
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
          w(:,locy) = smatvec ( A, w(:,locz) )
        case(3,4)

!         perform preconditioning operation y=Pz

!         action of Ap^-1 on pressure part
          call mi20_precondition ( prec%Ap, prec%coarse_data_Ap, w(ssp,locz), &
            wp2, prec%keep_Ap, prec%control_Ap, prec%info_Ap, prec%ma48_cntl )

          if (prec%info_Ap%flag < 0) then
            write(*,*) "Error return from mi20_precondition for Ap"
            stop
          end if

!         action of Fp on pressure part
          wp1 = smatvec ( Fp, wp2 )

!         action of Mp^-1 on pressure part
          call solve_diag ( Mp, wp1, wp2, so_diag, initsl=.false., &
            maxnumitsreached=ldummy )

          w(ssp,locy) = - wp2

!         add -B^Tp to velocity part
          wv1 = w(ssv,locz) - smatvec ( prec%BT, w(ssp,locy) )

!         perform action of F^-1 on velocity part
          call mi20_precondition ( prec%F, prec%coarse_data_F, wv1, &
            wv, prec%keep_F, prec%control_F, prec%info_F, prec%ma48_cntl )

          w(ssv,locy) = wv

          if (prec%info_F%flag < 0) then
            write(*,*) "Error return from mi20_precondition of F"
            stop
          end if

        case default

          call errormsg_case_default ( 'solve3_mi20', 'iact', int_value=iact )

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

    deallocate ( wv, wv1, wp1, wp2 )

  end subroutine solve3_mi20


! delete single prec

  subroutine delete_single_prec3_mi20 ( prec, sysmatrix_Ap )

    type (prec3_mi20_t), intent(inout) :: prec
    type (sysmatrix_t), intent(inout), optional :: sysmatrix_Ap

    if ( .not. allocated(prec%coarse_data_F) ) then
      write(*,'(/a/)') &
        'Error delete_single_prec_mi20: arrays in prec not allocated'
      stop
    end if

    call mi20_finalize ( prec%coarse_data_F, prec%keep_F, prec%control_F, &
      prec%info_F )

!   delete matrix F

    deallocate( prec%F%col, prec%F%val, prec%F%ptr )

!   delete matrix BT

    call delete ( prec%BT )

    call mi20_finalize ( prec%coarse_data_Ap, prec%keep_Ap, prec%control_Ap, &
      prec%info_Ap )

    if ( prec%move_alloc_used ) then

!     move matrix back to sysmatrix

      if ( .not. present(sysmatrix_Ap) ) then

        write(*,'(/3(a/))') &
          'Error delete_single_prec3_mi20: ', &
          '  matrix Ap cannot be moved back to the system matrix', &
          '  sysmatrix_Ap must be added to the argument list of the routine.'
        stop

      else if ( allocated(sysmatrix_Ap%Suu%a) ) then

        write(*,'(/3(a/))') &
          'Error delete_single_prec3_mi20: ', &
          '  matrix cannot be moved back to the system matrix', &
          '  sysmatrix_Ap already allocated'
        stop

      end if

      call move_alloc ( from=prec%Ap%col, to=sysmatrix_Ap%Suu%ja )
      call move_alloc ( from=prec%Ap%val, to=sysmatrix_Ap%Suu%a )
      call move_alloc ( from=prec%Ap%ptr, to=sysmatrix_Ap%Suu%ia )

    else

!     delete matrix

      deallocate( prec%Ap%col, prec%Ap%val, prec%Ap%ptr )

    end if

    prec%setup = .false.
    prec%move_alloc_used = .false.

  end subroutine delete_single_prec3_mi20

end module hsl_solve3_mi20_m
