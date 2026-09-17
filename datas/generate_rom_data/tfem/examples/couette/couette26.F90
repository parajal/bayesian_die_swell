#if FEAST

! Viscoelastic problem on a unit square (b formulation)
! linear Couette flow
! DEVSS-G/SUPG
! Building generalized eigenvalue problem
!    (A-lambda B)v = 0
! for linear stability analysis.
! Matrices A and B have a different non-zero matrix structure.
! FEAST eigenvalue solver.
! NOTE: Similar to couette24, but now using the b-formulation.

program couette26

  use tfem_m
  use viscoelastic_elements_m
  use io_utils_m
  use feast_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,     &  ! Q2 velocities
    pintpl = 4,     &  ! Q1 pressures
    gintpl = 4,     &  ! Q1 gradients
    bintpl = 4,     &  ! Q1 conformation
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    physqb = 4,     & ! physical quantity nr of the conformation tensor (mode 1)
    gauss = 3,      & ! 3x3 integration of quads
    ncompb = 4,     & ! number of conformation tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 2,      & ! UCM/Oldroyd-B model
    bvariant = 1      ! b-formulation:
                      ! 1: CDT 2: symmetric
                      ! 3: Cholesky 4: Cholesky with log

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrixA, sysmatrixB
  type(sysvector_t), target :: soliter
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscriptvec_t) :: cxx, cxy, cyy
  type(subscript_t) :: velx, vely, vel, pres, solsc, gxy
  type(subscript_t), dimension(ncompb,nmodes) :: bvalcmp
  type(subscript_t), dimension(nmodes) :: bval


! variables

  integer :: &
    nx=20,             & ! number of elements in x
    ny=20,             & ! number of elements in y
    timeint1 = 8         ! first-order time integration (BDF1)

  real(dp) :: &
    eta_s = 0.1_dp,     & ! solvent viscosity
    eta_p = 1.0_dp,     & ! polymer viscosity
    lambda = 1.0_dp,    & ! relaxation time
    beta = 1.0_dp         ! upwinding parameter in the SUPG method

! FEAST interface structures

  type(solver_options_feast_t) :: so
  type(arg_feast_t) :: arg

! FEAST variables in arg

  integer :: M0 = 10       ! Search subspace dimension (initial guess)
  complex(dp):: Emid = (0._dp,0._dp) ! centroid of the contour ellipse
  real(dp) :: r = 1._dp    ! horizontal radius of the contour ellipse

! FEAST options in so

  integer :: &
    printlevel = 1, &
    stop_conv = 8, &
    max_loop = 50, &
    only_right_eigv = 1, &
    mixed_precision=0

  real(dp) :: ratio = 1.0 ! stretch circle in y-direction (ellipse)

! other local variables

  integer :: i, m
  real(dp) :: alpha, G, csteady(3), lsteady(2), eigv(2,2), scsteady(3)
  logical, dimension(3+nmodes,3+nmodes) :: physqmaskA, physqmaskB


! namelist for input of variables; read from standard input

  namelist /comppar/ nx, ny, eta_s, eta_p, lambda, beta, &
    M0, Emid, r, printlevel, stop_conv, max_loop, only_right_eigv, &
    ratio, mixed_precision

  read ( unit=*, nml=comppar )


! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus

! physqmask: set zero blocks

! Matrix A
! ========

  physqmaskA = .true.

! Gup coupling for matrix A
  physqmaskA(1:3,1:3) = reshape ( &
      [  .true., .true., .false., &
         .true., .true.,  .true., &
        .false., .true., .false. ], [3,3], order=[2,1] )

! Gup coupling with c for matrix A
  call fill_vpG_b_physqmask ( physqmaskA )

! Internal coupling of c with other modes for matrix A
  call fill_offdiag_b_physqmask ( physqmaskA )

! Matrix B
! ========

  physqmaskB = .false.

! Only diagonal blocks of c are present
  do i = 1, nmodes
    physqmaskB(physqb+i-1,physqb+i-1) = .true.
  end do


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+2*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gauss,    bintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      0,        timeint1,    ( 0, i = 23, 150 )  &
    ]

  coefficients%i(71) = bvariant  ! variant for b-formulation
  coefficients%i(83) = physqb ! physical quantity if first mode
  coefficients%i(84) = 3  ! storage of b tensor: all modes in sysvector
  coefficients%i(86) = 1  ! sysvector number for velocity in SUPG
  coefficients%i(87) = 1  ! sysvector number for iteration b tensor
  coefficients%i(89) = 1  ! exclude time derivative from b build in matrix A

  coefficients%r = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
      0._dp, 0._dp,   0._dp,   beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G,     lambda &
    ]

! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=6, nphysq=4 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [ 4,0,4,0,4,0,4,0,0,    &  ! G
                  2,2,2,2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,1,0,0,    &  ! pressure
                  ncompb,0,ncompb,0,ncompb,0,ncompb,0,0, &  ! c
                  1,1,1,1,1,1,1,1,1,    &  ! scalar, such as vorticity
                  3,3,3,3,3,3,3,3,3 ], &  ! tensor
                  [9,6] )

  input_probdef%physq = [1,2,3,(3+m,m=1,nmodes)]
  input_probdef%probnr = 1


! Dirichlet boundary conditions

! velocity on lower boundary
  call define_essential ( mesh, input_probdef, curve1=1, physq=2 )
! velocity on upper boundary
  call define_essential ( mesh, input_probdef, curve1=3, physq=2 )
! pressure in point 1
  call define_essential ( mesh, input_probdef, point=1, physq=3 )

! constraints for periodical boundary conditions

! velocities
  call define_constraint ( mesh, input_probdef, &
    physq=2, curve1=2, curve2=5, discretization='collocation', exclude=3 )

! gradients
  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, curve2=5, discretization='collocation' )

! conformation tensor (use collocation)
  do m = 1, nmodes
    call define_constraint ( mesh, input_probdef, &
      physq=physqb+m-1, curve1=2, curve2=5, discretization='collocation' )
  end do

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for solution (excluding contraint forces)

  call create_subscript ( mesh, problem, solsc, &
    physqarr=[physqgrad,physqvel,physqpress,(physqb+m-1,m=1,nmodes)] )

! create vector subscripts for the gradient xy

  call create_subscript ( mesh, problem, gxy, physqarr=[physqgrad], degfd=2 )

! create vector subscripts for the velocity

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )

! create vector subscript for the pressure

  call create_subscript ( mesh, problem, pres, physqarr=[physqpress] )

! create vector subscript for the b tensor

  do m = 1, nmodes
    call create_subscript ( mesh, problem, bval(m), physqarr=[physqb+m-1] )
    do i = 1, ncompb
      call create_subscript ( mesh, problem, bvalcmp(i,m), &
        physqarr=[physqb+m-1], degfd=i )
    end do
  end do
  call create_subscript ( mesh, problem, cxx, degfd=1, vec=6 )
  call create_subscript ( mesh, problem, cxy, degfd=2, vec=6 )
  call create_subscript ( mesh, problem, cyy, degfd=3, vec=6 )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, soliter )
  call create_sysvector ( problem, rhsd )


! fill solution vector with essential boundary conditions

  soliter%u = 0
! set gxy = 1 on domain
  soliter%u(gxy%s) = 1
! set velx = y on domain
  call fill_sysvector ( mesh, problem, soliter, &
    node1=1, node2=mesh%nnodes, physq=2, degfd=1, func=func, funcnr=1 )

! initial solution b tensor

  csteady(1) = 1 + 2 * lambda ** 2
  csteady(2) = lambda
  csteady(3) = 1

  call eig2x2 ( csteady, lsteady, eigv )
  lsteady = sqrt ( lsteady )
  call inveig2x2 ( scsteady, lsteady, eigv )

!  do m = 1, nmodes
    m = 1
    soliter%u(bvalcmp(1,1)%s) = scsteady(1) ! initial bxx
    soliter%u(bvalcmp(2,1)%s) = scsteady(2) ! initial bxy
    soliter%u(bvalcmp(3,1)%s) = scsteady(2) ! initial byx
    soliter%u(bvalcmp(4,1)%s) = scsteady(3) ! initial byy
! end do

! create system matrices

  problem%physqmask = physqmaskA
  call create_sysmatrix_structure_base ( sysmatrixA, mesh, problem, &
    usephysqmask=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrixA, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrixA )

  call create_sysmatrix_data ( sysmatrixA )

  problem%physqmask = physqmaskB
  call create_sysmatrix_structure_base ( sysmatrixB, mesh, problem, &
    usephysqmask=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrixB, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrixB )

  call create_sysmatrix_data ( sysmatrixB )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=1, nprob=1 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => soliter
  oldvectors_ve%p(1)%p => problem  ! needed for implicit_ce_supg_elem

! matrix A
! ========

  problem%physqmask = physqmaskA

! build (assemble) matrix A for gradient/velocity/pressure part

  call build_vpG ( sysmatrixA )

! build (assemble) matrix A for conformation part

  do m = 1, nmodes

    coefficients%i(85) = m   ! set mode number

    call build_b ( sysmatrixA )

  end do

! matrix B
! ========

! build (assemble) matrix B for gradient/velocity/pressure part

  call build_vpG_zero ( sysmatrixB )

! build (assemble) matrix B for conformation part

  do m = 1, nmodes

    coefficients%i(85) = m   ! set mode number

    call build_b_timederiv ( sysmatrixB )

  end do

! check matrices A and B

  call check ( sysmatrixA )
  call check ( sysmatrixB )


! change sign of A to be the Jacobian of the rhs

  sysmatrixA%Suu%a = - sysmatrixA%Suu%a

! find eigenvalues/vectors

! search circle [Emod,r] including M eigenpairs
  arg%Emid = Emid
  arg%r = r
  arg%M0 = M0 !! M0>=M

! options
  so%printlevel = printlevel
  so%stop_conv = stop_conv
  so%max_loop = max_loop
  so%only_right_eigv = only_right_eigv
  so%ratio = ratio
  so%mixed_precision = mixed_precision

  call solve_system_real_gen_feast ( sysmatrixA, sysmatrixB, arg, &
    solver_options=so )

! write data for post-processing

  print *,'FEAST output info', arg%info
  print *,'FEAST output M0', arg%M0
  print *,'Eigenvalues/Residuals (inside contour)'
  do i=1, arg%M
    print *, i, arg%E(i), arg%res(i)
  enddo

! TODO: output of eigenvectors


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( soliter, rhsd )
  call delete ( sysmatrixA, sysmatrixB )
  call delete ( oldvectors_ve )

  do m = 1, nmodes
    do i = 1, ncompb
      call delete ( bvalcmp(i,m), bval(m) )
    end do
  end do
  call delete ( solsc, gxy )
  call delete ( arg )

contains

  subroutine fill_vpG_b_physqmask ( physqmask )

    logical, dimension(:,:), intent(inout) :: physqmask

    integer :: j

    do j = 1, nmodes
      physqmask(physqgrad,physqb+j-1) = .false.
    end do
    do j = 1, nmodes
      physqmask(physqpress,physqb+j-1) = .false.
      physqmask(physqb+j-1,physqpress) = .false.
    end do

  end subroutine fill_vpG_b_physqmask

  subroutine fill_offdiag_b_physqmask ( physqmask )

    logical, dimension(:,:), intent(inout) :: physqmask

    integer :: i, j

    if ( nmodes > 1 ) then
      do i = 1, nmodes
        do j = 1, nmodes
          if ( i == j ) cycle
            physqmask(physqb+i-1,physqb+j-1) = .false.
        end do
      end do
    end if

  end subroutine fill_offdiag_b_physqmask


  subroutine build_vpG ( sysmatrix )

    type(sysmatrix_t), intent(inout) :: sysmatrix

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients, &
      physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress], &
      buildvector=.false., usephysqmask=.true. )

!   DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_elem, coefficients=coefficients, addmatvec=.true., &
      physqrow=[physqgrad,physqvel], physqcol=[physqgrad,physqvel], &
      buildvector=.false. )

!   build -div(tau) part in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_divtau, oldvectors=oldvectors_ve, &
      physqrow=[physqvel], physqcol=[(physqb+m-1,m=1,nmodes)], &
      addmatvec=.true., coefficients=coefficients, buildvector=.false. )

!   periodical condition on velocities

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients, buildvector=.false. )

!   periodical condition on gradients

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients, buildvector=.false. )

  end subroutine build_vpG


  subroutine build_b ( sysmatrix )

    type(sysmatrix_t), intent(inout) :: sysmatrix

!   build (assemble) matrix for conformation tensor mode m

!   diagonal block
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
      physqrow=[physqb+m-1], physqcol=[physqb+m-1], &
      oldvectors=oldvectors_ve, addmatvec=.true., buildvector=.false. )

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_ce_vel_supg_elem, coefficients=coefficients, &
      physqrow=[physqb+m-1], physqcol=[physqgrad,physqvel], &
      oldvectors=oldvectors_ve, addmatvec=.true., buildvector=.false. )

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problem, sysmatrix, &
      sysvector=rhsd, constraint1=2+m, elemsub=stokes_constr_node_conn, &
      addmatvec=.true., buildvector=.false. )

  end subroutine build_b

  subroutine build_vpG_zero ( sysmatrix )

    type(sysmatrix_t), intent(inout) :: sysmatrix

!   periodical condition on velocities

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=zero_constr_node_conn, &
      coefficients=coefficients, buildvector=.false. )

!   periodical condition on gradients

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=zero_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients, buildvector=.false. )

  end subroutine build_vpG_zero

  subroutine build_b_timederiv ( sysmatrix )

    type(sysmatrix_t), intent(inout) :: sysmatrix

!   build (assemble) matrix for conformation tensor mode m

!   diagonal block of c
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_ce_timederiv_supg_elem, coefficients=coefficients, &
      physqrow=[physqb+m-1], physqcol=[physqb+m-1], &
      oldvectors=oldvectors_ve, addmatvec=.true., buildvector=.false. )

!   set to zero periodical condition on conformation tensor
    call build_system_constraint ( mesh, problem, sysmatrix, &
      sysvector=rhsd, constraint1=2+m, elemsub=zero_constr_node_conn, &
      addmatvec=.true., buildvector=.false. )

  end subroutine build_b_timederiv

  function func ( nr, x )
    use kind_defs_m
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func
    func = x(2)
  end function func

end program couette26

#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on feast', &
    ' - add the feast library for linking', &
    ' - set preprocessing macro FEAST in Mdefs.mk'
end
#endif
