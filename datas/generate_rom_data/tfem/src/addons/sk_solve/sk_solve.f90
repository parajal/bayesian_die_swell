
! Copyright (C) 2005-2016 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, routines, for solving the system of equations using iterative
! methods in Sparskit2 with various preconditioning.
!

module sk_solve_m

  use glob_defs_m
  use sparse_m
  use system_defs_m
  use timer_m

  implicit none


! detect benchmark and call toc

  logical :: benchmark_sk = .false.


! type definition for solver options

  type solver_options_sk_t

    logical :: stoponmaxmvm = .true.
                          ! stop (if stoponmaxmvm=.true.) with an error
                          ! message or continue (if stoponmaxmvm=.false.)
                          ! the program if the maximum number of
                          ! matrix-vector multiplies (maxmvm) has been reached.

    logical :: use_renumber = .false.
                          ! use the renumbered degrees of freedom in perm/iperm
                          ! of sysmatrix to reorder the rows/columns of the
                          ! matrix. The perm/ipermu arrays can e.g. be
                          ! set by MeTis.

    integer :: printlevel = 0 ! Print level:
                             !   0   Only error messages are printed
                             !   1   A little amount of information of the
                             !       iteration method is printed.
                             !       The first and last iterations are
                             !       printed
                             !   2   A maximal amount of information of the
                             !       iteration method is printed

    integer :: type_prec = 2 ! ipar(2) in the sparskit solvers:
                             !     -- status of the preconditioning:
                             !      0 == no preconditioning
                             !      1 == left preconditioning only
                             !      2 == right preconditioning only
                             !      3 == both left and right preconditioning

    integer :: maxmvm = 300  ! Parameter ipar(6) in the sparskit solvers:
                             !    Maximum number of matrix-vector multiplies.
                             !    If not a positive number the iterative
                             !    solver will run until the convergence test is
                             !    satisfied.

    integer :: itsolver = 5 ! The basic iterative solver:
                            ! (see SPARSKIT for a description)
                            !  1 CG       -- Conjugate Gradient Method
                            !  2 CGNR     -- Conjugate Gradient Method
                            !                (Normal Residual equation)
                            !  3 BCG      -- Bi-Conjugate Gradient Method
                            !  4 DBCG     -- BCG with partial pivoting
                            !  5 BCGSTAB  -- BCG stabilized
                            !  6 TFQMR    -- Transpose-Free Quasi-Minimum
                            !                Residual method
                            !  7 FOM      -- Full Orthogonalization Method
                            !  8 GMRES    -- Generalized Minimum RESidual method
                            !  9 FGMRES   -- Flexible version of Generalized
                            !                Minimum RESidual method
                            ! 10 DQGMRES  -- Direct versions of Quasi Generalize
                            !                Minimum Residual method

    integer :: mgmres = 0 !  Dimension of Krylov space for GMRES.
                          !  Used by FOM, GMRES, FGMRES, DQGMRES.

    integer :: preconditioner = 1 ! The preconditioner:
                                  ! (see SPARSKIT for a description)
                                  !  1 ILUT    : Incomplete LU factorization
                                  !              with dual truncation strategy
                                  !  2 ILUTP   : ILUT with column pivoting
                                  !  5 ILUK    : level-k ILU
                                  !  6 ILU0    : simple ILU(0) preconditioning
                                  !  7 MILU0   : MILU(0) preconditioning

    real(dp) :: prec_store = 1._dp ! Parameter for storage of the
                                 ! preconditioning matrix:
                                 ! prec_store * NNZ reals are reserved,
                                 ! where NNZ is the size of the matrix.

    real(dp) :: droptol = 1e-5_dp  ! Sets the threshold for dropping small terms
                                   ! in the preconditioning matrix.

    real(dp) :: permtol = 0.01_dp ! tolerance ratio used to determine whether
                                  ! or not to permute two columns.
                                  ! At step i columns i and j are permuted when
                                  !    abs(a(i,j))*permtol > abs(a(i,i))
                                  ! note: permtol=0 --> never permute
                                  ! good values 0.1 to 0.01.
                                  ! Only used for preconditioner=2 (ILUTP).

    real(dp) :: fillin = 1._dp ! Parameter for computing the value of the
       ! fill-in parameter lfil.
       ! The definition of lfil depend on the preconditioner:
       ! preconditoner=1 (ILUT):
       !     Each row of L and each row of U will
       !     have a maximum of lfil elements (excluding the
       !     diagonal element).
       !     The value of lfil is determined as follows:
       !     if FILLIN>0 then
       !       lfil = nint( FILLIN * NNZ / N / 2 )
       !     if FILLIN<0 then
       !       lfil = nint( - FILLIN )
       !     where NNZ is the size of the matrix and N is the
       !     number of unknowns.
       ! preconditoner=5 (ILUK):
       !     Each element whose level-of-fill exceeds lfil during the
       !     ILU process is dropped.
       !     The value of lfil is determined as follows:
       !       lfil = nint( abs(FILLIN) )
       ! for other precondioners lfil is not used.

    real(dp) :: eps_rel = 1e-7_dp, eps_abs = 0._dp
        ! Relative tolerance and absolute tolerance, such that
        !  || residual || <= eps_rel * || initial residual || + eps_abs

  end type solver_options_sk_t


! type definition of storage for ILU factorization

  type ilu_sk_t

    logical :: decomp = .false.  ! ILU decomposition done?

    integer, allocatable, dimension(:) :: iperm ! permutation array of ILUTP

    integer, allocatable, dimension(:) :: ju

    integer, allocatable, dimension(:) :: jlu

    real(dp), allocatable, dimension(:) :: alu

  end type ilu_sk_t


! interface for generic delete subroutine

  interface delete
    module procedure delete_ilu_sk
  end interface delete

! interface for generic set_solver_options subroutine

  interface set_solver_options
    module procedure set_solver_options_sk
  end interface set_solver_options

contains


! Helper routine for setting the solver options

  subroutine set_solver_options_sk ( solver_options, keep, stoponmaxmvm, &
    printlevel, type_prec, maxmvm, itsolver, mgmres, preconditioner, &
    prec_store, droptol, fillin, eps_rel, eps_abs, use_renumber, permtol )


!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_sk_t), intent(inout) :: solver_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see solver_options_sk_t
    logical, intent(in), optional :: stoponmaxmvm, use_renumber
    integer, intent(in), optional :: printlevel, type_prec, maxmvm, itsolver, &
      mgmres, preconditioner
    real(dp), intent(in), optional :: prec_store, droptol, fillin, eps_rel, &
      eps_abs, permtol


    type(solver_options_sk_t) :: solveropt


    if ( present(keep) ) then
      if ( .not. keep ) solver_options = solveropt
    else
      solver_options = solveropt
    end if

    if ( present(stoponmaxmvm) ) solver_options%stoponmaxmvm = stoponmaxmvm
    if ( present(printlevel) ) solver_options%printlevel = printlevel
    if ( present(type_prec) ) solver_options%type_prec = type_prec
    if ( present(maxmvm) ) solver_options%maxmvm = maxmvm
    if ( present(itsolver) ) solver_options%itsolver = itsolver
    if ( present(mgmres) ) solver_options%mgmres = mgmres
    if ( present(preconditioner) ) &
                              solver_options%preconditioner = preconditioner
    if ( present(prec_store) ) solver_options%prec_store = prec_store
    if ( present(droptol) ) solver_options%droptol = droptol
    if ( present(fillin) ) solver_options%fillin = fillin
    if ( present(eps_rel) ) solver_options%eps_rel = eps_rel
    if ( present(eps_abs) ) solver_options%eps_abs = eps_abs
    if ( present(use_renumber) ) solver_options%use_renumber = use_renumber
    if ( present(permtol) ) solver_options%permtol = permtol

  end subroutine set_solver_options_sk


! Solve the system of equations using sparkskit

  subroutine solve_system_sk ( sysmatrix, rhsd, sol, ilu, initsol, &
    solver_options, maxmvmreached )

    type(sysmatrix_t), intent(inout) :: sysmatrix

!   Right-hand side vector. The part corresponding to the unknowns is used as
!   the right-hand side in the system:
!
!      Suu sol_u = rhsd_u   (*)
!
!   The effect of essential boundary conditions must already have been taken
!   into account in the rhsd vector
    type(sysvector_t), intent(inout) :: rhsd

!   The solution vector. The unknown part sol_u is filled with the solution of
!   the system (*), given above. The prescribed part sol_p is not modified.
    type(sysvector_t), intent(inout) :: sol

!   If present, ilu is used to store the ILU factorization. The next time the
!   solve routine is called with ilu present no new preconditioner is build
!   and the old one is used. This happens even when the matrix has changed.
!   If ilu is not present or has been deallocated using delete, a new ILU
!   preconditioner is build.
    type(ilu_sk_t), intent(inout), optional :: ilu

!   If .true.: use initial estimate for solution vector given in sol
!   default = .false.
    logical, intent(in), optional :: initsol

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_sk_t), intent(in), optional :: solver_options

!   This is an OUTPUT parameter indicating whether the maximum number
!   of matrix-vector multiplies (maxmvm) has been reached.
!   Note, that this parameter only makes sense in combination with
!   solver_options%stoponmaxmvm=.false.
    logical, intent(out), optional :: maxmvmreached


    logical :: precond, initsl, lmaxmvmreached
    integer :: numundegfd
    type(ilu_sk_t) :: ilw  ! local work space when ilu is not present
    type(solver_options_sk_t) :: lsolveropt ! local options


!   some testing first

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(/a/)') &
        'Error in solve_system_sk: no system matrix structure.'
      stop
    end if

    if ( .not. sysmatrix%finalized ) then
      write(*,'(/a/)') &
        'Error in solve_system_sk: system matrix has not been finalized.'
      stop
    end if

    if ( .not. sysmatrix%allocated_data ) then
      write(*,'(/a/)') &
        'Error in solve_system_sk: data in system matrix not allocated.'
      stop
    end if

    if ( sysmatrix%symmetric ) then
      write(*,'(/a/)') &
        'Error in solve_system_sk: symmetric matrix not allowed.'
      stop
    end if

    if ( .not. rhsd%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_sk: rhsd not created.'
      stop
    end if

    if ( .not. sol%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_sk: sol not created.'
      stop
    end if

!   initialize local solver options if solver_options in heading
    if ( present(solver_options) ) lsolveropt = solver_options

    if ( lsolveropt%use_renumber .and. .not. sysmatrix%renumber ) then
      write(*,'(/a/a/)') &
        'Error in solve_system_sk:', &
        '  use_renumber = T, whereas no renumbering is present in sysmatrix.'
      stop
    end if


!   build preconditioner?

    precond = .true.

    if ( present(ilu) ) then
      if ( ilu%decomp ) then
!       decomposition is already done
        precond = .false.
      end if
    end if

!   inital estimate?

    if ( present(initsol) ) then
      initsl = initsol
    else
      initsl = .false.
    end if

!   number of unknowns

    numundegfd = sysmatrix%Suu%n

!   set initial estimate

    if ( .not. initsl ) then
!     clear vector
      sol%u(1:numundegfd) = 0
    end if

!   solve, solution returned in sol(1:numundegfd)

    if ( present(ilu) ) then

      call solve_sk ( sysmatrix%Suu, rhsd%u(1:numundegfd), sol%u(1:numundegfd),&
        ilu, precond, lsolveropt, lmaxmvmreached, sysmatrix%permu, &
        sysmatrix%ipermu )

!     for next call of this routine: do solution only (no incomplete LU)
      ilu%decomp = .true.

    else

      call solve_sk ( sysmatrix%Suu, rhsd%u(1:numundegfd), sol%u(1:numundegfd),&
        ilw, precond, lsolveropt, lmaxmvmreached,sysmatrix%permu, &
        sysmatrix%ipermu )

      call delete ( ilw )

    end if

    if ( present(maxmvmreached) ) maxmvmreached = lmaxmvmreached

  end subroutine solve_system_sk


! solve linear system A u = b

  subroutine solve_sk ( A, b, u, ilu, precond, so, maxmvmreached, perm, iperm )

!   The matrix A in CSR format
    type(sparsematrix_t), intent(inout) :: A

!   the right-hand side vector b.
    real(dp), intent(inout), dimension(:) :: b

!   the solution u = A^{-1} b.
    real(dp), intent(inout), dimension(:) :: u

!   ilu is used to store the ILU factorization.
    type(ilu_sk_t), intent(inout), optional :: ilu

!   indicates whether new preconditioner must be build
    logical, intent(in) :: precond

!   The options/parameters for the solver
    type(solver_options_sk_t), intent(in) :: so

!   Indicated whether the maximum number of matrix-vector multiplies (maxmvm)
!   has been reached.
    logical, intent(out) :: maxmvmreached

!   permutation arrays for the renumbering
    integer, dimension(:), intent(in) :: perm, iperm


    integer :: n, nnz, iwk, niwork, lfil, m, nwork
    integer :: ipar(16)
    integer, dimension(:), allocatable :: jw, levs
    real(dp) :: fpar(16)
    real(dp), dimension(:), allocatable :: w, wk

!   reorder?

    if ( so%use_renumber ) then

       u = u(perm)
       b = b(perm)
       call reorder_rows_sparsematrix ( iperm, A )
       call reorder_columns_sparsematrix ( iperm, A )

    end if

!   set dimensions of problem

    n  = A%n
    nnz = A%nnz

!   preconditioner build

    if ( precond ) then

!     work arrays

      select case ( so%preconditioner )
        case(1,2)
          niwork = 2 * n
          iwk = nint( so%prec_store * nnz )
        case(5)
          niwork = 3 * n
          iwk = nint( so%prec_store * nnz )
        case(6,7)
          niwork = n
          iwk = nnz + 1
        case default
          write(*,'(a,i0)') ' Error solve_sk: Invalid preconditioner: ', &
            so%preconditioner
          stop
      end select

!     allocate memory for ilu

      allocate ( ilu%ju(n), ilu%jlu(iwk), ilu%alu(iwk) )

      if ( so%preconditioner == 2 ) then
        allocate ( ilu%iperm(2*n) )
      else
        allocate ( ilu%iperm(0) )
      end if

!     allocate temporary arrays and fill-in parameter lfil

      allocate ( jw(niwork), w(n) )

      if ( so%preconditioner == 5 ) then
        allocate ( levs(iwk) )
!       fill-in parameter
        lfil = nint( abs(so%fillin) )
      else
        allocate ( levs(0) )
!       fill-in parameter
        if ( so%fillin > 0 ) then
!          relative number
           lfil = nint( so%fillin * nnz / n / 2 )
        else
           lfil = nint( abs( so%fillin ) )
        end if
      end if

!     build preconditioner

      call sk_prec ( n, A%a, A%ja, A%ia, lfil, so%droptol, ilu%alu, ilu%jlu, &
        ilu%ju, levs, iwk, w, jw, so%permtol, ilu%iperm, so%preconditioner )

      deallocate ( jw, w, levs )

    else if ( so%preconditioner == 2 ) then

!     ilutp: permuted matrix and solution

      A%ja = ilu%iperm(A%ja+n)
      u = u(ilu%iperm(1:n))

    end if

    if ( timer .and. benchmark_sk ) call toc ( 'preconditioner' )

!   set parameters

    if ( so%itsolver <= 6 ) then
       m = 0
    else if ( so%mgmres <= 1 ) then
       m = 15
    else
       m = so%mgmres
    end if

    select case ( so%itsolver )
      case(1)
        nwork = 5 * n    !CG
      case(2)
        nwork = 5 * n    !CGNR
      case(3)
        nwork = 7 * n    !BCG
      case(4)
        nwork = 11 * n   !DBCG
      case(5)
        nwork = 8 * n    !BCGSTAB
      case(6)
        nwork = 11 * n   !TFQMR
      case(7)
        nwork = (n+3)*(m+2) + (m+1)*m/2   !FOM
      case(8)
        nwork = (n+3)*(m+2) + (m+1)*m/2   !GMRES
      case(9)
        nwork = 2*n*(m+1) + (m+1)*m/2 + 3*m + 2   !FGMRES
      case(10)
        nwork = n + (m+1) * (2*n+4) !DQGMRES
      case default
        call errormsg_case_default ( 'solve_sk', 'so%itsolver', &
          int_value=so%itsolver )
    end select

!   allocate workspace

    allocate ( wk(nwork) )

!   solve system

    ipar(2) = so%type_prec
    ipar(3) = 1
    ipar(4) = nwork
    ipar(5) = m
    ipar(6) = so%maxmvm
    fpar(1) = so%eps_rel
    fpar(2) = so%eps_abs
    fpar(11) = 0  ! set initial flops to zero

    call runsol0( n, b, u, ipar, fpar, wk, A%a, A%ja, A%ia, ilu%alu, ilu%jlu, &
      ilu%ju, so%itsolver, so%printlevel, so%stoponmaxmvm, maxmvmreached )

    deallocate ( wk )

    if ( so%preconditioner == 2 ) then

!     ilutp: permuted matrix and solution

      A%ja = ilu%iperm(A%ja)
      u = u(ilu%iperm(n+1:2*n))

    end if

!   reorder?

    if ( so%use_renumber ) then

       u = u(iperm)
       b = b(iperm)
       call reorder_rows_sparsematrix ( perm, A )
       call reorder_columns_sparsematrix ( perm, A )

    end if

    if ( timer .and. benchmark_sk ) call toc ( 'matrix-vector' )

  end subroutine solve_sk


! delete single lu

  subroutine delete_single_ilu_sk ( ilu )

    type (ilu_sk_t) :: ilu

    if ( .not. allocated(ilu%ju) ) then
      write(*,'(/a/)') &
        'Error delete_single_ilu_sk: arrays in ilu not allocated'
      stop
    end if

    deallocate ( ilu%ju, ilu%jlu, ilu%alu, ilu%iperm )

    ilu%decomp = .false.

  end subroutine delete_single_ilu_sk


! Delete lu

  subroutine delete_ilu_sk ( ilu_sk_1, ilu_sk_2, ilu_sk_3, &
    ilu_sk_4, ilu_sk_5 )

    type(ilu_sk_t), intent(inout) :: ilu_sk_1
    type(ilu_sk_t), intent(inout), optional :: ilu_sk_2, ilu_sk_3, &
      ilu_sk_4, ilu_sk_5

    call delete_single_ilu_sk(ilu_sk_1)
    if ( present(ilu_sk_2) ) call delete_single_ilu_sk(ilu_sk_2)
    if ( present(ilu_sk_3) ) call delete_single_ilu_sk(ilu_sk_3)
    if ( present(ilu_sk_4) ) call delete_single_ilu_sk(ilu_sk_4)
    if ( present(ilu_sk_5) ) call delete_single_ilu_sk(ilu_sk_5)

  end subroutine delete_ilu_sk


! Compute preconditioner for sparskit

  subroutine sk_prec ( n, a, ja, ia, lfil, droptol, alu, jlu, ju, levs, &
    iwk, w, jw, permtol, iperm, preconditioner )

    integer :: n, ia(n+1), ju(n), lfil, iwk, preconditioner
    integer, dimension(:) :: ja, jlu, levs, jw, iperm
    real(dp), dimension(:) :: a, alu
    real(dp) :: w(n), droptol, permtol


    select case ( preconditioner )
      case (1,2,5)

!       ILUT, ILUTP or ILU(k)
        call prspar1 ( preconditioner, n, a, ja, ia, lfil, droptol, alu, jlu, &
          ju, levs, iwk, w, jw, permtol, iperm )

      case(6,7)

!       ILU(0) or MILU(0)
!       Note: used jw instead of iw.

        call prspar6 ( preconditioner, n, a, ja, ia, alu, jlu, ju, jw )

      case default

        write(*,'(/a,i0,a/)') 'Error in sk_prec: Preconditioner ', &
          preconditioner, ' not available.'
        stop

    end select

  end subroutine sk_prec


! Choose solver to solve problem Ax=b iteratively (SPARSKIT)

  subroutine runsol0( n, rhs, sol, ipar, fpar, wk, a, ja, ia, alu, jlu, ju, &
    itsolver, printlevel, stoponmaxmvm, maxmvmreached )

    integer :: n, ipar(16), ia(n+1), itsolver
    integer, dimension(:) :: ja, ju, jlu
    real(dp) :: fpar(16)
    real(dp), dimension(n) :: rhs, sol
    real(dp), dimension(:) :: wk, a, alu
    integer :: printlevel
    logical :: stoponmaxmvm, maxmvmreached



    external cg, cgnr, bcg, dbcg, bcgstab, tfqmr, fom, gmres, fgmres, dqgmres


    select case ( itsolver )
      case(1)
        call runsol( n, rhs, sol, ipar, fpar, wk, a, ja, ia, alu, jlu, ju, &
          cg, printlevel, stoponmaxmvm, maxmvmreached )
      case(2)
        call runsol( n, rhs, sol, ipar, fpar, wk, a, ja, ia, alu, jlu, ju, &
          cgnr, printlevel, stoponmaxmvm, maxmvmreached )
      case(3)
        call runsol( n, rhs, sol, ipar, fpar, wk, a, ja, ia, alu, jlu, ju, &
          bcg, printlevel, stoponmaxmvm, maxmvmreached )
      case(4)
        call runsol( n, rhs, sol, ipar, fpar, wk, a, ja, ia, alu, jlu, ju, &
          dbcg, printlevel, stoponmaxmvm, maxmvmreached )
      case(5)
        call runsol( n, rhs, sol, ipar, fpar, wk, a, ja, ia, alu, jlu, ju, &
          bcgstab, printlevel, stoponmaxmvm, maxmvmreached )
      case(6)
        call runsol( n, rhs, sol, ipar, fpar, wk, a, ja, ia, alu, jlu, ju, &
          tfqmr, printlevel, stoponmaxmvm, maxmvmreached )
      case(7)
        call runsol( n, rhs, sol, ipar, fpar, wk, a, ja, ia, alu, jlu, ju, &
          fom, printlevel, stoponmaxmvm, maxmvmreached )
      case(8)
        call runsol( n, rhs, sol, ipar, fpar, wk, a, ja, ia, alu, jlu, ju, &
          gmres, printlevel, stoponmaxmvm, maxmvmreached )
      case(9)
        call runsol( n, rhs, sol, ipar, fpar, wk, a, ja, ia, alu, jlu, ju, &
          fgmres, printlevel, stoponmaxmvm, maxmvmreached )
      case(10)
        call runsol( n, rhs, sol, ipar, fpar, wk, a, ja, ia, alu, jlu, ju, &
          dqgmres, printlevel, stoponmaxmvm, maxmvmreached )
      case default
        call errormsg_case_default ( 'runsol0', 'itsolver', int_value=itsolver )
    end select

  end subroutine runsol0


! Solve problem Ax=b iteratively (SPARSKIT)

  subroutine runsol ( n, rhs, sol, ipar, fpar, wk, a, ja, ia, alu, jlu, ju, &
    solver, printlevel, stoponmaxmvm, maxmvmreached )

    integer :: n, ipar(16), ia(n+1)
    integer, dimension(:) :: ja, ju, jlu
    real(dp) :: fpar(16)
    real(dp), dimension(n) :: rhs, sol
    real(dp), dimension(:) :: wk, a, alu
    integer :: printlevel
    logical :: stoponmaxmvm, maxmvmreached

    external solver

    integer :: iiter

    maxmvmreached = .false.

    iiter  = 0

    ipar(1) = 0

 10 call solver ( n, rhs, sol, ipar, fpar, wk )

    if ( ipar(7) == iiter+1 ) then
      if ( ( printlevel ==1 .and. ipar(7) == 1 ) .or.  printlevel == 2 ) then
        print *, 'Matrix-vector multiply: ', ipar(7), ' Residual: ', fpar(5)
      end if
      iiter = iiter + 1
    end if

    if ( ipar(1) == 1 ) then
       call amux ( n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
       goto 10
    else if ( ipar(1) == 2 ) then
       call atmux ( n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
       goto 10
    else if ( ipar(1) == 3 .or. ipar(1) == 5 ) then
       call lusol ( n, wk(ipar(8)), wk(ipar(9)), alu, jlu, ju)
       goto 10
    else if ( ipar(1) == 4 .or. ipar(1) == 6 ) then
       call lutsol ( n, wk(ipar(8)), wk(ipar(9)), alu, jlu, ju)
       goto 10
    else if ( ipar(1) <= 0 ) then
       if ( printlevel >= 1 .and. ipar(1) == 0 ) then
          print *, 'Iterative solver has satisfied convergence test.'
          print * ,'At matrix-vector multiply: ', ipar(7),  &
            ' Residual: ', fpar(5)
       else if ( ipar(1) == -1  ) then
         if ( stoponmaxmvm ) then
           print *, 'Iterative solver has iterated too many times.'
           stop
         else
          print * ,'At matrix-vector multiply: ', ipar(7),  &
            ' Residual: ', fpar(5)
           maxmvmreached = .true.
         end if
       else if ( ipar(1) == -2 ) then
          print *, 'Iterative solver was not given enough work space.'
          print *, 'The work space should at least have ', ipar(4), ' elements.'
          stop
       else if ( ipar(1) == -3 ) then
          print *, 'Iterative solver is facing a break-down.'
          stop
       else if ( ipar(1) /= 0 ) then
          print *, 'Iterative solver terminated. code =', ipar(1)
          stop
       end if
    end if

  end subroutine runsol


! Compute ILUT, ILUTP or ILU(k) preconditioner for sparskit

  subroutine prspar1( ipre, n, a, ja, ia, lfil, droptol, alu, jlu, ju, levs, &
    iwk, w, jw, permtol, iperm )

!   see source of routine ilut and ilutp in sparskit for a
!   description of heading parameters

    integer :: ipre, n, ja(*), ia(n+1), jlu(*), ju(n), levs(*), jw(*), lfil, &
      iwk, iperm(*)
    real(dp) :: a(*), alu(*), w(n), droptol, permtol


    integer :: ierr


    select case ( ipre )
      case(1)

!       ILUT preconditioner

        call ilut( n, a, ja, ia, lfil, droptol, alu, jlu, ju, iwk, w, jw, ierr )

      case(2)

!       ILUTP preconditioner

        call ilutp( n, a, ja, ia, lfil, droptol, permtol, n, alu, jlu, ju, &
          iwk, w, jw, iperm, ierr )

      case(5)

!       ILU(k) preconditioner

        call iluk( n, a, ja, ia, lfil, alu, jlu, ju, levs, iwk, w, jw, ierr )

      case default

        call errormsg_case_default ( 'prspar1', 'ipre', int_value=ipre )

    end select

!   deal with error code

    if ( ierr /= 0 ) then

      write(*,'(/a,i0/)') 'Preconditioning failed: error code = ', ierr

      select case ( ierr )
        case (1:)
          write(*,'(/a,i0/)') 'Zero pivot encountered at step number ', ierr
        case(-1)
          write(*,'(/3(a/))') 'Error. input matrix may be wrong.', &
                  '(The elimination process has generated a', &
                  'row in L or U whose length is >= n)'
        case(-2)
          write(*,'(/a/)') 'The matrix L overflows the array al'
        case(-3)
          write(*,'(/a/)') 'The matrix U overflows the array alu.'
        case(-4)
          write(*,'(/a/)') 'Illegal value for lfil.'
        case(-5)
          write(*,'(/a/)') 'Zero row encountered'
        case default
          call errormsg_case_default ( 'prspar1', 'ierr', int_value=ierr )
      end select

      stop

    end if

  end subroutine prspar1


! Compute ILU(0) or MILU(0) preconditioner for sparskit

  subroutine prspar6 ( ipre, n, a, ja, ia, alu, jlu, ju, iw )

    integer :: ipre, n, ja(*), ia(n+1), jlu(*), ju(n), iw(n)
    real(dp) :: a(*), alu(*)

!     see source of routine ilut in sparskit for a
!     description of heading parameters

    integer :: ierr


    select case ( ipre )
      case(6)

!       ILU(0) preconditioner

        call ilu0 ( n, a, ja, ia, alu, jlu, ju, iw, ierr )

      case(7)

!       MILU(0) preconditioner

        call milu0 ( n, a, ja, ia, alu, jlu, ju, iw, ierr )

      case default

        call errormsg_case_default ( 'prspar6', 'ipre', int_value=ipre )

      end select

!     deal with error code

      if (ierr /= 0) then

        write(*,'(/a,i0/)') 'Preconditioning failed: error code = ', ierr

        if ( ierr >= .0 ) then
           write(*,'(/a,i0/)') 'Zero pivot encountered at step number ', ierr
        end if

        stop

      end if

  end subroutine prspar6

end module sk_solve_m
