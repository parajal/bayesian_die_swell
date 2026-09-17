
! Types, routines, for solving the system of equations using interative
! methods in MKL with ILU0 or ILUT preconditioning.

module mkl_gmres_m

  use kind_defs_m
  use sparse_m
  use system_defs_m

  implicit none


! type definition for solver options

  type solver_options_mkl_t

    logical :: stoponmaxmvm = .true. ! stop (if stoponmaxmvm=.true.) with an
                                     ! error message or continue
                                     ! (if stoponmaxmvm=.false.)
                                     ! the program if the maximum number of
                                     ! iterations (maxmvm) has been reached

    integer :: printlevel = 0 ! Print level:
                              ! 0 = Only error messages are printed
                              ! 1 = A little amount of information is printed
                              ! 2 = A maximal amount of information is printed

    integer :: maxmvm = 300   ! Parameter ipar(5) in the mkl_fgmres solver:
                              ! specifies the maximum number of iterations

    integer :: mgmres = 0     ! Obsolete Parameter. Previously ipar(15) in
                              ! the mkl_fgmres solver was set to this value.
                              ! The default value for ipar(15) is used now.

    integer :: preconditioner = 1 ! The preconditioner:
                                  ! 0 = no preconditioner used
                                  ! 1 = ILU0
                                  ! 2 = ILUT

    integer :: maxfil = 0 ! Maximum fill-in, half of the bandwidth of the
                          ! ILUT preconditioner: the number of non-zero elements
                          ! in the rows cannot exceed (2*maxfil+1)

    real(dp) :: fillin = 1._dp ! Parameter for computing the used value for
                               ! maxfil if maxfil is preset to zero (maxfil=0):
                               !    maxfil = nint( FILLIN * NNZ / N / 2 )
                               ! where NNZ is the size of the matrix and N is
                               ! the number of unknowns. In words: if fillin
                               ! is near 1 the size of the ILU matrix will be
                               ! roughly equal to the size of the matrix.

    real(dp) :: droptol = 1e-3_dp ! Sets the threshold for dropping small terms
                                  ! in the ILUT preconditioning matrix

    real(dp) :: eps_rel = 1e-9_dp ! Parameter dpar(1) in the mkl_fgmres solver:
                                  ! specifies the relative tolerance

    real(dp) :: eps_abs = 1e-6_dp ! Parameter dpar(2) in the mkl_fgmres solver:
                                  ! specifies the absolute tolerance

    integer ::  zero_diag_par = 0 ! Parameter ipar(31) in the mkl_ilu0 and
                                  ! mkl_ilut:
                                  ! if a zero diagonal element occurs then:
                                  ! 0 = the calculation is stopped
                                  ! 1 = the zero value is set to
                                  !     zero_diagonal_value

    real(dp) :: zero_diag_comp = 1e-2_dp  ! Parameter dpar(31) in the mkl_ilu0
               ! and mkl_ilut_t:
               ! for ILU0: specifies the small value that is compared with
               !           the diagonal elements; if the value of the
               !           diagonal element is smaller, it is set to
               !           zero_diag_value, if zero_diag_par = 1
               ! for ILUT: specifies the value that being multiplied
               !           by the matrix row norm is assigned to those
               !           diagonal elements whose value is less than
               !           the matrix row norm multiplied by the value
               !           of the parameter droptol, if zero_diag_par = 1

    real(dp) :: zero_diag_value = 1e-10_dp ! Parameter dpar(32) in the mkl_ilu0:
                   ! specifies the value that is assigned to the diagonal
                   ! element if its value is less than zero_diag_comp

  end type solver_options_mkl_t

  integer :: ipar(128)

  real(dp) :: dpar(128)

  type ilu_mkl_t

    logical :: decomp = .false.  ! ILU decomposition done?

    integer :: preconditioner = 0

    real(dp), allocatable, dimension(:) :: alu

    integer, allocatable, dimension(:) :: ialu

    integer, allocatable, dimension(:) :: jalu

  end type ilu_mkl_t

! interface for generic delete subroutine

  interface delete
    module procedure delete_ilu_mkl
  end interface delete

! interface for generic set_solver_options subroutine

  interface set_solver_options
    module procedure set_solver_options_mkl
  end interface set_solver_options

contains


! Helper routine for setting the solver options

  subroutine set_solver_options_mkl ( solver_options, keep, stoponmaxmvm, &
    printlevel, maxmvm, mgmres, preconditioner, maxfil, droptol, eps_rel, &
    eps_abs, zero_diag_par, zero_diag_comp, zero_diag_value, fillin )

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_mkl_t), intent(inout) :: solver_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see solver_options_mkl_t
    logical, intent(in), optional :: stoponmaxmvm
    integer, intent(in), optional :: printlevel, maxmvm, mgmres, &
      preconditioner, maxfil, zero_diag_par
    real(dp), intent(in), optional :: droptol, eps_rel, eps_abs, &
      zero_diag_comp, zero_diag_value, fillin

    type(solver_options_mkl_t) :: solveropt

    if ( present(keep) ) then
      if ( .not. keep ) solver_options = solveropt
    else
      solver_options = solveropt
    end if

    if ( present(stoponmaxmvm) ) solver_options%stoponmaxmvm = stoponmaxmvm
    if ( present(printlevel) ) solver_options%printlevel = printlevel
    if ( present(maxmvm) ) solver_options%maxmvm = maxmvm
    if ( present(mgmres) ) solver_options%mgmres = mgmres
    if ( present(preconditioner) ) &
                            solver_options%preconditioner = preconditioner
    if ( present(maxfil) ) solver_options%maxfil = maxfil
    if ( present(droptol) ) solver_options%droptol = droptol
    if ( present(eps_rel) ) solver_options%eps_rel = eps_rel
    if ( present(eps_abs) ) solver_options%eps_abs = eps_abs
    if ( present(zero_diag_par) ) solver_options%zero_diag_par = zero_diag_par
    if ( present(zero_diag_comp) ) &
                              solver_options%zero_diag_comp = zero_diag_comp
    if ( present(zero_diag_value) ) &
                              solver_options%zero_diag_value = zero_diag_value
    if ( present(fillin) ) solver_options%fillin = fillin

  end subroutine set_solver_options_mkl


! Solve the system of equations using iterative solver in MKL

  subroutine solve_system_mkl_gmres ( sysmatrix, rhsd, sol, ilu, initsol, &
    solver_options, maxmvmreached )

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

!   If present, ilu is used to store the ILU factorization. The next time the
!   solve routine is called with ilu present no new preconditioner is build
!   and the old one is used. This happens even when the matrix has changed.
!   If ilu is not present or has been deallocated using delete, a new ILU
!   preconditioner is build.
    type(ilu_mkl_t), intent(inout), optional :: ilu

!   If .true.: use initial estimate for solution vector given in sol
!   default = .false.
    logical, intent(in), optional :: initsol

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_mkl_t), intent(in), optional :: solver_options

!   This is an OUTPUT parameter indicating whether maximum number
!   of matrix-vector multiplies (maxmvm) has been reached.
!   Note, that this parameter only makes sense in combination with
!   solver_options%stoponmaxmvm=.false.
    logical, intent(out), optional :: maxmvmreached


    logical :: initsl, lmaxmvmreached
    integer :: numundegfd

    type(ilu_mkl_t) :: ilw  ! local work space when ilu is not present
    type(solver_options_mkl_t) :: lsolveropt ! local options


!   some testing first

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(/a/)') &
        'Error in solve_system_mkl_gmres: no system matrix structure.'
      stop
    end if

    if ( .not. sysmatrix%finalized ) then
      write(*,'(/a/)') &
        'Error in solve_system_mkl_gmres: system matrix has not been finalized.'
      stop
    end if

    if ( .not. sysmatrix%allocated_data ) then
      write(*,'(/a/)') &
        'Error in solve_system_mkl_gmres: data in system matrix not allocated.'
      stop
    end if

    if ( sysmatrix%symmetric ) then
      write(*,'(/a/)') &
        'Error in solve_system_mkl_gmres: symmetric matrix not allowed.'
      stop
    end if

    if ( .not. rhsd%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_mkl_gmres: rhsd not created.'
      stop
    end if

    if ( .not. sol%created ) then
      write(*,'(/a/)') &
        'Error in solve_system_mkl_gmres: sol not created.'
      stop
    end if

!   initialize local solver options if solver_options in heading
    if ( present(solver_options) ) lsolveropt = solver_options

!   number of unknowns

    numundegfd = sysmatrix%Suu%n

!   initial estimate?

    if ( present(initsol) ) then
      initsl = initsol
    else
      initsl = .false.
    end if

!   set initial estimate

    if ( .not. initsl ) then
!     clear vector
      sol%u(1:numundegfd) = 0.0_dp
    end if

!   solve, solution returned in sol%u(1:numundegfd)

    if ( present(ilu) ) then

!     build preconditioner?
      if ( .not. ilu%decomp ) then

        call build_mkl_ilu ( sysmatrix%Suu, ilu, lsolveropt )

      end if

!     for next call of this routine: do solution only (no incomplete LU)
      ilu%decomp = .true.

      call solve_mkl_gmres ( sysmatrix%Suu, rhsd%u(1:numundegfd), &
        sol%u(1:numundegfd), numundegfd, lsolveropt, lmaxmvmreached, ilu )

    else if ( lsolveropt%preconditioner /= 0 ) then

!     build preconditioner
      call build_mkl_ilu ( sysmatrix%Suu, ilw, lsolveropt )

      call solve_mkl_gmres ( sysmatrix%Suu, rhsd%u(1:numundegfd), &
        sol%u(1:numundegfd), numundegfd, lsolveropt, lmaxmvmreached, ilw )

      call delete ( ilw )

    else

!     no preconditioner used

      call solve_mkl_gmres ( sysmatrix%Suu, rhsd%u(1:numundegfd), &
        sol%u(1:numundegfd), numundegfd, lsolveropt, lmaxmvmreached )

    end if

    if ( present(maxmvmreached) ) maxmvmreached = lmaxmvmreached

  end subroutine solve_system_mkl_gmres


! Build ilu preconditioner

  subroutine build_mkl_ilu ( A, ilu, so )

    type(sparsematrix_t), intent(inout) :: A
    type(ilu_mkl_t), intent(inout) :: ilu
    type(solver_options_mkl_t), intent(in) :: so

    integer :: ierr, nwork, lmaxfil

    ipar(31) = so%zero_diag_par
    dpar(31) = so%zero_diag_comp
    dpar(32) = so%zero_diag_value

!   set maxfil
    if ( so%maxfil <= 0 ) then
!     use fillin parameter
      lmaxfil = nint( so%fillin * A%nnz / A%n / 2 )
    else
      lmaxfil = so%maxfil
    end if

    call sort_sparsematrix ( A )

    if ( so%preconditioner == 1 ) then

      allocate ( ilu%alu(A%nnz) )

!     build ilu0 preconditioner
      call dcsrilu0 ( A%n, A%a, A%ia, A%ja, ilu%alu, ipar, dpar, ierr )

      ilu%preconditioner = 1

    else if ( so%preconditioner == 2 ) then

      nwork = (2*lmaxfil + 1)*A%n - lmaxfil*(lmaxfil + 1) + 1

      allocate ( ilu%alu(nwork), ilu%ialu(A%n), ilu%jalu(nwork) )

!     build ilut preconditioner
      call dcsrilut ( A%n, A%a, A%ia, A%ja, ilu%alu, ilu%ialu, ilu%jalu, &
        so%droptol, lmaxfil, ipar, dpar, ierr )

      ilu%preconditioner = 2

    end if

  end subroutine build_mkl_ilu


! Solve iteratively by using MKL_FGMRES

  subroutine solve_mkl_gmres ( A, b, u, n, so, maxmvmreached, ilu )

!   the matrix A in CSR format
    type(sparsematrix_t), intent(inout) :: A

!   the right-hand sude vector b
    real(dp), intent(in), dimension(:) :: b

!   the solution u = A^{-1} b
    real(dp), intent(out), dimension(:) :: u

!   number of unknowns
    integer :: n

!   The options/parameters for the solver
    type(solver_options_mkl_t), intent(in) :: so

!   Indicate whether the maximum number of matrix-vector multiplies (maxmvm)
!   has been reached
    logical, intent(out) :: maxmvmreached

!   ilu is used to store ILU factorization
!   if not present, no preconditioner is used
    type(ilu_mkl_t), intent(inout), optional :: ilu

    real(dp) :: trvec(n)
    real(dp), allocatable, dimension(:) :: tmp
    integer :: itercount, rci_request, ipar15

    maxmvmreached = .false.
    ipar15 = min(150,n)

    allocate (tmp(n*(2*ipar15+1)+(ipar15*(ipar15+9))/2+1))

!   Initialize the solver

    call dfgmres_init(n, u, b, rci_request, ipar, dpar, tmp)
    if ( rci_request /= 0 ) then
      print *, 'The routine dfgmres_init failed to complete the task'
      stop
    end if

!   Set solver paramaters

    ipar(5) = so%maxmvm
    ipar(9) = 1          ! residual stopping test is performed
    ipar(10) = 0         ! user-defined stopping test not performed
    ipar(12) = 1         ! zero-norm stopping test is performed
    dpar(1) = so%eps_rel
    dpar(2) = so%eps_abs

    if ( present( ilu ) ) then
      ipar(11) = 1       ! preconditioned version
    else
      ipar(11) = 0       ! non-preconditioned version
    end if

!   Check the correctness and consistency of the newly set parameters

    call dfgmres_check(n, u, b, rci_request, ipar, dpar, tmp)
    if ( rci_request /= 0 ) then
      print *, 'The routine dfgmres_check failed to complete the task'
      stop
    end if

!   Compute the solution by FGMRES solver

1   call dfgmres(n, u, b, rci_request, ipar, dpar, tmp)

    if ( rci_request == 0) then

!     If rci_request = 0, then the solution was found with the required
!     precision

      call dfgmres_get(n, u, b, rci_request, ipar, dpar, tmp, itercount)
      if ( so%printlevel >= 1 ) then
        print *, 'Iterative solver has satisfied convergence test'
        print *, 'At iteration: ', itercount, 'with residual: ', dpar(5)
      end if

    else if ( rci_request == 1 ) then

!     If rci_request = 1, then compute the vector A*TMP(IPAR(22)) and put
!     the result in vector TMP(IPAR(23))

      call mkl_dcsrgemv('N', n, A%a, A%ia, A%ja, tmp(ipar(22)), tmp(ipar(23)))

      if ( (ipar(4) > 0 .and. so%printlevel == 2) .or. &
           (ipar(4) == 1 .and. so%printlevel == 1) ) then
        print *,'At iteration: ', ipar(4)-1, 'the residual is: ', dpar(5)
      end if

      goto 1

    else if ( rci_request == 2 ) then

!     If rci_request = 2 (if ipar(10) = 1), the user should perform
!     the stopping test

!     insert user-defined stopping test here

    else if ( rci_request == 3 ) then

!     If rci_request = 3, apply the preconditioner on the vector TMP(IPAR(22))
!     and put the result in vector TMP(IPAR(23))

      if ( so%preconditioner == 1 ) then

        call mkl_dcsrtrsv('L','N','U', n, ilu%alu, A%ia, A%ja, &
          tmp(ipar(22)), trvec )
        call mkl_dcsrtrsv('U','N','N', n, ilu%alu, A%ia, A%ja, &
          trvec, tmp(ipar(23)) )

      else if ( so%preconditioner == 2 ) then

        call mkl_dcsrtrsv('L','N','U', n, ilu%alu, ilu%ialu, ilu%jalu, &
          tmp(ipar(22)), trvec )
        call mkl_dcsrtrsv('U','N','N', n, ilu%alu, ilu%ialu, ilu%jalu, &
          trvec, tmp(ipar(23)) )

      end if

      goto 1

    else if ( rci_request == 4 ) then

!     If rci_request = 4 (if ipar(12) = 1), the user should perform the
!     test for zero norm of the currently generated vector

!     insert user-defined test here

    else if ( rci_request < 0 ) then

!     If rci_request = anything else, then DFGMRES subroutine failed to
!     compute the solution vector

      if ( rci_request == -1 ) then
        if ( so%stoponmaxmvm ) then
          print *, 'Iterative solver has iterated too many times'
          stop
        else
          maxmvmreached = .true.
        end if
      else if ( rci_request == -10 ) then
        print *, 'ERROR: an attempt to divide by zero occurs'
        stop
      else if ( rci_request == -11 ) then
        print *, 'ERROR: the routine enters the infinite cycle'
        stop
      else if ( rci_request == -12 ) then
        print *, 'ERROR: errors are found in the method parameters'
        stop
      end if

    end if

  end subroutine solve_mkl_gmres


! Delete single ilu preconditioner

  subroutine delete_single_ilu_mkl ( ilu )

    type (ilu_mkl_t) :: ilu

    if ( .not. allocated(ilu%alu) ) then
      write(*,'(/a/)') &
        'Error delete_single_ilu_mkl: arrays in ilu not allocated'
      stop
    end if

    if ( ilu%preconditioner == 1 ) then

      deallocate ( ilu%alu )

    else if ( ilu%preconditioner == 2 ) then

      deallocate ( ilu%alu, ilu%ialu, ilu%jalu )

    end if

    ilu%decomp = .false.
    ilu%preconditioner = 0

  end subroutine delete_single_ilu_mkl


! Delete ilu preconditioner

  subroutine delete_ilu_mkl ( ilu_mkl_1, ilu_mkl_2, ilu_mkl_3, &
    ilu_mkl_4, ilu_mkl_5 )

    type(ilu_mkl_t), intent(inout) :: ilu_mkl_1
    type(ilu_mkl_t), intent(inout), optional :: ilu_mkl_2, ilu_mkl_3, &
      ilu_mkl_4, ilu_mkl_5

    call delete_single_ilu_mkl(ilu_mkl_1)
    if ( present(ilu_mkl_2) ) call delete_single_ilu_mkl(ilu_mkl_2)
    if ( present(ilu_mkl_3) ) call delete_single_ilu_mkl(ilu_mkl_3)
    if ( present(ilu_mkl_4) ) call delete_single_ilu_mkl(ilu_mkl_4)
    if ( present(ilu_mkl_5) ) call delete_single_ilu_mkl(ilu_mkl_5)

  end subroutine delete_ilu_mkl

end module mkl_gmres_m
