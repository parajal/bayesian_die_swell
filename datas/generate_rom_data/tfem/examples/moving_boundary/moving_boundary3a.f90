! Moving boundary due to viscoelastic normal stresses and gravity forces.
! 2D domain, axisymmetric with 3D velocities.
! Wall(s) and bottom rotate with an angular velocities omega.
! Perfect slip (in plane) on the walls.
! Perfect slip (in r and theta) on the bottom.
! Surface tension on fluid-air interface and fluid-wall interface.
! Viscoelastic fluid.
! DEVSS-G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! L2-projection of c=exp(s) on the finite element approximation space in the
! momentum balance.
! second-order time integration to steady state
! Similar to moving_boundary3, but with changes:
!      P1 line elements for height function.
!      Note, this is only for illustration. The standard case works best!

program moving_boundary3a

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m
  use surface_advection_elements_m
  use update_mesh_nodes1_m
  use figplot_m

  implicit none

! constants

! T: Single wall, symmetry line (planar) or center line (axisymmetric)
!    Note, that for the axisymmetric case R1=0 in the mesh generation!
! F: Two walls
  logical, parameter :: symmetric = .false.

! constants flow problem

  integer, parameter :: &
    uintpl = 8,  & ! Q2 velocities
    pintpl = 4,  & ! Q1 pressures
    gintpl = 4,  & ! Q1 gradients
    cintpl = 4,  & ! Q1 conformation
    physqg = 1,  & ! physical quantity nr of the gradients
    physqv = 2,  & ! physical quantity nr of the velocities
    physqp = 3,  & ! physical quantity nr of the pressures
    gauss = 3,   & ! 3x3 integration of quads
    gaussb = 3,  & ! 3-point integration of boundary elements
    ncompc = 6,  & ! number of conformation tensor components
    nmodes = 1,  & ! number of modes
    startm = 501   ! start of material model data

! constants surface advection

  integer, parameter :: &
    hintpl = 2,       & ! P1 height function
    vintpl = 6,       & ! P2 velocity
    ninti_sf_adv = 3, & ! number of Gauss points
    method = 1          ! discretization method 0: Galerkin, 1: SUPG


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc, input_probdefc_proj
  type(problem_t), target :: problem, problemc, problemc_proj
  type(sysmatrix_t) :: sysmatrix, sysmatrixc, sysmatrixc_proj
  type(sysvector_t), target :: sol, solm1
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, gammadot
  type(vector_t) :: tauviscous, tauviscoelastic, ctensor
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vel8, velz, velr, velt, cval
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc, solcm1, solc_proj
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc, rhsc_proj
  type(lu_ma41_t) :: luc
  type(solver_options_ma41_t) :: solver_options_u, solver_options_c
  type(lu_ma57_t) :: lu_exps_proj

! 1D height function for surface advection
  type(mesh_t) :: mesh_sf_adv
  type(input_probdef_t) :: input_probdef_sf_adv
  type(problem_t) :: problem_sf_adv
  type(sysmatrix_t) :: sysmatrix_sf_adv
  type(sysvector_t), target :: sol_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1, &
    sol_sf_adv_pred
  type(sysvector_t) :: rhsd_sf_adv
  type(oldvectors_t) :: oldvectors_sf_adv, oldvectors_sf_adv_deriv
  type(coefficients_t) :: coefficients_sf_adv
  type(vector_t), target :: velocity_sf_adv
  type(vector_t) :: height_sf_adv
  type(subscript_t) :: hgt, hgt_end

  type(vector_t), target :: meshvel
  real(dp), allocatable, dimension(:,:) :: meshcoor_initial

! ALE mesh motion problem

  type(problem_t), target :: problem_update_mesh
  real(dp), allocatable, dimension(:) :: disp

! variables

  integer :: &
    timeint1 = 1,        & ! first-order time integration, first time step
    timeint2 = 7,        & ! second-order time integration after first time step
    numtimesteps = 125,  & ! number of time steps
    logc = 1,            & ! standard scheme or log transformation
    model = 3,           & ! 2: Oldroyd-B 3: Giesekus
    coorsys = 1, & ! planar Cartesian (0) or axisymmetric (1) coor. system
    vfuncnr = 1, & ! function nr for gravity 0=none 1=rho g (-1,0)
    plot_mesh_every = 10    ! plot mesh every .. steps

  real(dp) :: &
    rho = 1.0_dp,      & ! density
    grv = 1.0_dp,      & ! gravity constant g
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 2.0_dp,   & ! relaxation time
    alphapar = 0.1_dp, & ! alpha parameter in the Giesekus model
    gamma_fa = 0.1_dp, & ! surface tension coefficient fluid-air interface
    dgamma_w = 0.01_dp,& ! surface tension difference fluid-air wrt the wall
    omega = 1.0_dp,    & ! angular rate of the inner wall
    deltat = 5.e-2_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    betah = 1.0_dp,    & ! beta parameter for SUPG of the height function
    rs_gup = 1.5_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.0_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 1.0_dp       ! integer_storage for conformation LU (HSL)

  integer :: icomp, step, i, npar
  integer :: vertices(4) = [1,3,5,7]
  integer :: nnodes, curvew
  real(dp) :: alpha, G

  real(dp), allocatable, dimension(:,:) :: meshcoor_n, meshcoor_nm1

  character(len=20) :: filename

  logical :: surface_tension = .true.  ! include surface tension


! namelist for input of variables; read from standard input

  namelist /comppar/ numtimesteps, plot_mesh_every, &
    rho, grv, eta_s, eta_p, lambda, alphapar, gamma_fa, dgamma_w, omega, &
    deltat, rs_gup, is_gup, rs_c, is_c, surface_tension

  read ( unit=*, nml=comppar )


! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus

  if ( model == 2 ) then
    npar = 2 ! Oldroyd-B
  else if ( model == 3 ) then
    npar = 3 ! Giesekus
  end if


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+npar*nmodes )

  coefficients%i = &
    [ uintpl, pintpl, 0, 0,      gintpl, &
      physqv, physqp, 0, physqg, gauss,  &
      gaussb,  cintpl, 0, 0,      0,      &
      0,      0,      model, nmodes, startm, &
      logc,   timeint1, ( 0, i = 23, 150 )  &
    ]

  coefficients%i(14) = vfuncnr
  coefficients%i(23) = coorsys
  coefficients%i(49) = 1  ! exp(s) projection =.true. for logc=1
  coefficients%i(61) = 0  ! projected G=0, direct velocity gradient=1 in CE
  coefficients%i(67) = 1  ! 3D velocity (swirl)

  coefficients%r(1:502) = &
    [ eta_s,    0._dp,   0._dp, alpha, 0._dp, &
      0._dp, 0._dp,  deltat, beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G, lambda  &
    ]

  if ( model == 3 ) then
    coefficients%r(503) = alphapar ! Giesekus
  end if

  coefficients%vfunc => vfunc

! create mesh

  call read_mesh ( mesh, filename='mesh2b.out' )

  call fill_mesh_parts ( mesh )

  if ( coorsys==1 ) then
    if ( symmetric .and. abs(mesh%coor(mesh%points(1),2)) > tiny(0._dp) ) then
      write(*,'(a/)') &
        'Error: r coordinates of the center line non-zero. Change mesh.'
      stop
    end if
    if ( mesh%coor(mesh%points(1),2) < 0._dp ) then
      write(*,'(a/)') 'Error: negative r coordinates. Change mesh.'
      stop
    end if
  end if

! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=7, nphysq=3 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 6  ! gradients
  input_probdef%vec_elementdof(1)%a(:,2) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5) = 6  ! symmetric tensor
  input_probdef%vec_elementdof(1)%a(:,6) = 9  ! unsymmetric tensor
  input_probdef%vec_elementdof(1)%a(:,7) = 2  ! vector for mesh velocity

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1


! Dirichlet boundary conditions

! bottom (perfect slip on r and theta directions)
  call define_essential ( mesh, input_probdef, &
    curve1=9, physq=physqv, degfd=[1,0,0] )
! walls (perfect slip in z-direction)
  call define_essential ( mesh, input_probdef, &
    curve1=10, physq=physqv, degfd=[0,1,1] )

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for the velocity components for monitoring

  call create_subscript ( mesh, problem, velz, physqarr=[physqv], degfd=1 )
  call create_subscript ( mesh, problem, velr, physqarr=[physqv], degfd=2 )
  call create_subscript ( mesh, problem, velt, physqarr=[physqv], degfd=3 )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=3, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a(:,1) = 0
  input_probdefc%vec_elementdof(1)%a(vertices,1) = 1  ! c
  input_probdefc%vec_elementdof(1)%a(:,2) = 1         ! scalar for plotting
  input_probdefc%vec_elementdof(1)%a(:,3) = 6         ! symmetric tensor

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

  call problem_definition ( input_probdefc, mesh, problemc )

! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemc, cval, physqarr=[1] )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, solm1, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0
! inner wall
  call fill_sysvector ( mesh, problem, sol, curve1=1, physq=physqv, degfd=3, &
    func=vel_theta, funcnr=1 )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! define some arrays and vector for the ALE mesh movement

  allocate ( meshcoor_initial(mesh%nnodes,2) )
  meshcoor_initial = mesh%coor

  allocate ( meshcoor_n(mesh%nnodes,2), meshcoor_nm1(mesh%nnodes,2) )
  meshcoor_n = mesh%coor
  meshcoor_nm1 = mesh%coor

  call create_vector ( problem, meshvel, vec=7 )
  meshvel%u = 0._dp  ! initialize


! create system vectors (solution and right-hand side) for conformation and
! initialize vectors to zero stress

  call create ( problemc, solc, solcm1, rhsc )

  if ( logc == 0 ) then ! standard
    solc(1,1)%u = 1
    solc(2,1)%u = 0
    solc(3,1)%u = 0
    solc(4,1)%u = 1
    solc(5,1)%u = 0
    solc(6,1)%u = 1
  else if ( logc == 1 ) then ! log scheme
    do i = 1, ncompc
      solc(i,1)%u = 0
    end do
  end if

! create system matrix for conformation problem

  call create_sysmatrix_structure_base ( sysmatrixc, mesh, problemc )
  call finalize_sysmatrix_structure ( sysmatrixc )

  call create_sysmatrix_data ( sysmatrixc )


! problem definition for projected "c=exp(s)" of the log conformation s

  call create_input_probdef ( mesh, input_probdefc_proj, nvec=1, nphysq=1 )

  input_probdefc_proj%vec_elementdof(1)%a(:,1) = 0
  input_probdefc_proj%vec_elementdof(1)%a(vertices,1) = 1  ! c

  input_probdefc_proj%physq = [1]
  input_probdefc_proj%probnr = 3

  call problem_definition ( input_probdefc_proj, mesh, problemc_proj )

  call create ( problemc_proj, solc_proj, rhsc_proj )

! initialize solc_proj
  solc_proj(1,1)%u = 1
  solc_proj(2,1)%u = 0
  solc_proj(3,1)%u = 0
  solc_proj(4,1)%u = 1
  solc_proj(5,1)%u = 0
  solc_proj(6,1)%u = 1


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec2=3, nprob=3, &
    nvec=1 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s(2)%p => solm1
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%s2(2)%p => solcm1
  oldvectors_ve%s2(3)%p => solc_proj
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc
  oldvectors_ve%p(3)%p => problemc_proj
  oldvectors_ve%v(1)%p => meshvel


! create subscript for velocity sampling on curve 8

  call create_subscript ( mesh, problem, vel8, physqarr=[physqv], curves=[8] )


! fill coefficients for surface advection problem

  call create_coefficients ( coefficients_sf_adv, ncoefi=100, ncoefr=50 )

  coefficients_sf_adv%i = 0
  coefficients_sf_adv%i(1:6) = [ ninti_sf_adv, 2, 0, method, 0, hintpl ]
  coefficients_sf_adv%i(10) = vintpl
  coefficients_sf_adv%i(12) = 1 ! use shapefunc for velocity for element shape

  coefficients_sf_adv%r = 0

  coefficients_sf_adv%r(4) = deltat
  coefficients_sf_adv%r(5) = betah

! create mesh for surface advection

  call curve_to_1Dmesh ( mesh, mesh%curves(8), dim1D=2, pnts=[2,5], &
    mesh1D=mesh_sf_adv )

  call fill_mesh_parts ( mesh_sf_adv )


! problem definition for surface advection

  call create_input_probdef ( mesh_sf_adv, input_probdef_sf_adv, nvec=3, &
    nphysq=1 )

  input_probdef_sf_adv%vec_elementdof(1)%a(:,1) = [1,0,1] ! height function
  input_probdef_sf_adv%vec_elementdof(1)%a(:,2) = 2 ! velocity
  input_probdef_sf_adv%vec_elementdof(1)%a(:,3) = 1 ! height function all nodes

  input_probdef_sf_adv%physq = [1]
  input_probdef_sf_adv%probnr = 4

! define problem

  call problem_definition ( input_probdef_sf_adv, mesh_sf_adv, problem_sf_adv )

  call create_sysvector ( problem_sf_adv, sol_sf_adv, rhsd_sf_adv )
  call create_sysvector ( problem_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1 )
  call create_sysvector ( problem_sf_adv, sol_sf_adv_pred )

! fill solution vector with essential boundary conditions

  sol_sf_adv%u = mesh%coor(mesh%points(2),1)
  sol_sf_adv_n%u = sol_sf_adv%u
  sol_sf_adv_pred%u = sol_sf_adv%u

! create system matrix

  call create_sysmatrix_structure ( sysmatrix_sf_adv, mesh_sf_adv, &
    problem_sf_adv )

  call create_sysmatrix_data ( sysmatrix_sf_adv )


  ! create vectors

  call create_vector ( problem_sf_adv, velocity_sf_adv, vec=2 )
  call create_vector ( problem_sf_adv, height_sf_adv, vec=3 )

  call create_oldvectors ( oldvectors_sf_adv, nsysvec=2, nvec=1 )

  oldvectors_sf_adv%s(1)%p => sol_sf_adv_n    ! corrector at n
  oldvectors_sf_adv%s(2)%p => sol_sf_adv_nm1  ! corrector at nm1

  oldvectors_sf_adv%v(1)%p => velocity_sf_adv ! advection velocity at np1

  call create_oldvectors ( oldvectors_sf_adv_deriv, nsysvec=1 )
  oldvectors_sf_adv_deriv%s(1)%p => sol_sf_adv_pred    ! prediction

! create subscript for the height values

  call create ( mesh_sf_adv, problem_sf_adv, hgt )
  call create ( mesh_sf_adv, problem_sf_adv, hgt_end, points=[1,2] )


! create and build system matrix for projection problem
! NOTE matrix remains constant and needs to be build once.

  call create_sysmatrix_structure ( sysmatrixc_proj, mesh, problemc_proj, &
    symmetric=.true. )
  call create_sysmatrix_data ( sysmatrixc_proj )



! allocate arrays for the ALE displacement problem

  allocate ( disp(mesh_sf_adv%nnodes) )


! open files

  open ( unit=11, file='height.out', status='unknown' )
  open ( unit=23, file='out', recl=300 )

! write initial mesh

  write(filename,'(a,i4.4,a)') 'mesh', 0, '.vtk'
  call write_mesh_vtk ( mesh, filename )


! time stepping

  coefficients_sf_adv%i(5) = 1 ! start with first-order scheme

  do step = 1, numtimesteps

    if ( step == 2 ) then
!     change time integration scheme at the second time step
      coefficients%i(22) = timeint2
    end if

    if ( step >= 2 ) then

      coefficients_sf_adv%i(5) = 2 ! second-order scheme

!     predict position of the surface and adapt mesh accordingly

      call update_mesh_surface_predictor

    end if


!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG


!   exps projection (for log conformation)

    if ( logc == 1 ) call solve_exps_projection


!   build implicit terms of CE with rhs in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      addmatvec=.true., coefficients=coefficients )


!   surface tension

    if ( surface_tension ) then

!     build surface integral on free surface

      coefficients%r(19) = gamma_fa

      call add_boundary_elements ( mesh, problem, rhsd, &
        elemsub=surface_tension_curve, &
        curve=8, coefficients=coefficients, physq=[physqv] )

!     build surface integral on walls

      if ( symmetric ) then
        curvew = 6  ! only upper wall
      else
        curvew = 10 ! both walls
      end if

      coefficients%r(19) = dgamma_w

      call add_boundary_elements ( mesh, problem, rhsd, &
        elemsub=surface_tension_curve, &
        curve=curvew, coefficients=coefficients, physq=[physqv] )

    end if


    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

    call copy ( sol, solm1 )


!   build (assemble) matrix and vector for conformation problem

    if ( coefficients%i(22) == timeint1 ) then

      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    else

      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem_implicit_2nd_order, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    end if


    call check ( sysmatrixc )

    call copy ( solc, solcm1 )


!   solve conformation and keep LU decomposition in the loop over components

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do icomp = 1, ncompc
      call solve_system_ma41 ( sysmatrixc, rhsc(icomp,1), solc(icomp,1), luc, &
        solver_options=solver_options_c  )
    end do

    call delete ( luc )  ! remove LU decomposition and rebuild next time step


!   solve surface advection (corrector)

    call solve_surface_height_corrector

    if ( mod(step,plot_mesh_every) == 0 ) then
      write(filename,'(a,i4.4,a)') 'mesh', step, '.vtk'
      call write_mesh_vtk ( mesh, filename )
    end if

!   write monitor data

    write(unit=23,fmt=*) &
      step, step*deltat, maxval(solc(1,1)%u(cval%s)), &
      maxval(abs(sol%u(velz%s))), maxval(abs(sol%u(velr%s))), &
      maxval(abs(sol%u(velt%s)))

  end do

  close(unit=11)
  close(unit=23)

! post-processing

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )
  oldvectors%s(1)%p => sol

  call create_vector ( problem, velocity, physq=physqv )
  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, gammadot, vec=4 )
  call create_vector ( problem, tauviscous, vec=5 )

  call extract_physvector ( mesh, problem, sol, velocity )

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=11

  call derive_vector ( mesh, problem, gammadot, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, tauviscous, elemsub=stokes_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='moving_boundary3.vtk' )

  call write_vector_vtk ( mesh, problem, filename='moving_boundary3.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true., &
    assume3D=.true. )

  call write_vector_vtk ( mesh, problem, filename='moving_boundary3.vtk', &
    dataname='velocity_vector_planar', vector=velocity, append=.true., &
    degfd=[1,2] )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='moving_boundary3.vtk', dataname='gammadot', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='moving_boundary3.vtk', &
    dataname='tauviscous', vector=tauviscous, append=.true., assume3D=.true. )

! viscoelastic

  call create_vector ( problemc, tauviscoelastic, vec=3 )

  call derive_vector ( mesh, problemc, tauviscoelastic, &
    elemsub=deriv_viscoelastic_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemc, filename='moving_boundary3.vtk', &
    dataname='tauviscoelastic', vector=tauviscoelastic, append=.true., &
    assume3D=.true. )

  call create_vector ( problemc, ctensor, vec=3 )

  call derive_vector ( mesh, problemc, ctensor, &
    elemsub=deriv_conformation_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemc, filename='moving_boundary3.vtk', &
    dataname='c', vector=ctensor, append=.true., &
    assume3D=.true. )

! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, solm1, rhsd )
  call delete ( sysmatrix )
  call delete ( vel8 )
  call delete ( velocity, pressure, gammadot, tauviscous, &
                tauviscoelastic, ctensor )
  call delete ( oldvectors, oldvectors_ve )

  deallocate ( meshcoor_initial, meshcoor_n, meshcoor_nm1 )

  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( solc, solcm1, rhsc )

  call delete ( problemc_proj )
  call delete ( input_probdefc_proj )
  call delete ( solc_proj, rhsc_proj )

  call delete ( coefficients )
  call delete ( velz, velr, velt, cval )

  call delete ( mesh_sf_adv )
  call delete ( problem_sf_adv )
  call delete ( input_probdef_sf_adv )
  call delete ( sol_sf_adv, rhsd_sf_adv )
  call delete ( sol_sf_adv_pred, sol_sf_adv_n, sol_sf_adv_nm1 )
  call delete ( sysmatrix_sf_adv )
  call delete ( oldvectors_sf_adv )
  call delete ( coefficients_sf_adv )
  call delete ( hgt, hgt_end )

! ALE mesh motion problem

  deallocate ( disp )
  if ( numtimesteps > 1 ) then
    call delete ( problem_update_mesh )
  end if

contains


! Convert curve to a 1D mesh with only a single group.

  subroutine curve_to_1Dmesh ( mesh, geometry, dim1D, pnts, mesh1D )

!   the mesh that contains the geometry
    type(mesh_t), intent(in) :: mesh

!   the geometry of the curve
    type(geometry_t), intent(in) :: geometry

!   the coordinate direction that form the 1D curve coordinates
    integer, intent(in) :: dim1D

!   create points in the new mesh from the points in mesh
    integer, dimension(:), intent(in) :: pnts

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh1D

    real(dp), dimension(size(pnts),1) :: points

!   make skeleton mesh
    call mesh_skeleton ( mesh=mesh1D, nnodes=geometry%nnodes, &
      nelem=geometry%nelem, elshape=geometry%element%elshape, &
      ndim=1, callname='curve_to_1Dmesh' )

!   fill coordinates
    mesh1D%coor(:,1) = mesh%coor(geometry%nodes,dim1D)

!   fill topology
    mesh1D%topology(1)%a = geometry%topology(:,:,1)

!   add points

    points(:,1) = mesh%coor(mesh%points(pnts),dim1D)

    call add_to_mesh ( mesh1D, points=points )

  end subroutine curve_to_1Dmesh

! velocity theta function

  function vel_theta ( nr, xin )

    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: xin
    real(dp) :: vel_theta

    select case ( nr )

    case (1)

!     rotation with omega

      vel_theta  = omega * xin(2)

    case default

      write(*,'(/a,i0/)') 'Error vel_theta: wrong function number: ', nr
      stop

    end select

  end function vel_theta

! body force function

  function vfunc ( n, nr, xin )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: xin
    real(dp), dimension(n) :: vfunc

    select case ( nr )

    case (1)

!     gravity

      vfunc(1)  = - rho * grv
      vfunc(2:) = 0

    case default

      write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
      stop

    end select

  end function vfunc


! update the (ALE) mesh based on a predictor for the surface position

  subroutine update_mesh_surface_predictor

!   predictor for surface height
    sol_sf_adv_pred%u = 2._dp*sol_sf_adv_n%u - sol_sf_adv_nm1%u

    call derive_vector ( mesh_sf_adv, problem_sf_adv, height_sf_adv, &
      elemsub=height_function_deriv, &
      coefficients=coefficients_sf_adv, oldvectors=oldvectors_sf_adv_deriv )

    meshcoor_nm1 = meshcoor_n
    meshcoor_n = mesh%coor

    disp = height_sf_adv%u - mesh%coor(mesh%curves(8)%nodes,1)

!   update nodes of the mesh (only in the first direction)

    call update_mesh_nodes ( mesh, problem_update_mesh, disp, &
      el=mesh_sf_adv%points([1,2]), eg=mesh%points([2,5]), cf=8 )

    call find_bounds_blocks ( mesh )

!   mesh velocity

    meshvel%u = reshape ( transpose ( &
        ( 1.5_dp*mesh%coor - 2*meshcoor_n + 0.5_dp*meshcoor_nm1 ) / deltat ), &
                         [2*mesh%nnodes] )

  end subroutine update_mesh_surface_predictor


! solve convection equation for the surface height (corrector)

  subroutine solve_surface_height_corrector

    use postprocessing_m

    real(dp) :: max_height, min_height, end_height(2)
    type(solver_options_ma41_t) :: solver_options_h
    real(dp), dimension(:,:), allocatable :: tmp

!   get velocities on the moving curve

    nnodes = mesh_sf_adv%nnodes
    tmp = reshape ( sol%u(vel8%s), [3,nnodes] )
    velocity_sf_adv%u = reshape ( tmp([2,1],:), [2*nnodes] )

    call build_system ( mesh_sf_adv, problem_sf_adv, sysmatrix_sf_adv, &
      rhsd_sf_adv, elemsub=surface_advection_elem, &
      oldvectors=oldvectors_sf_adv, coefficients=coefficients_sf_adv )

    call check_filled_sysmatrix ( sysmatrix_sf_adv )

    call add_effect_of_essential_to_rhs ( problem_sf_adv, sysmatrix_sf_adv, &
      sol_sf_adv, rhsd_sf_adv )


!   MA41 solver storage

    solver_options_h%integer_storage = 2.0
    solver_options_h%real_storage    = 2.0

!   solve system

    call solve_system_ma41 ( sysmatrix_sf_adv, rhsd_sf_adv, sol_sf_adv, &
                             solver_options=solver_options_h )


    call copy ( sol_sf_adv_n, sol_sf_adv_nm1 )
    call copy ( sol_sf_adv, sol_sf_adv_n )

!   compute the maximum and end height

    max_height = maxval ( sol_sf_adv%u )
    min_height = minval ( sol_sf_adv%u )
    end_height = sol_sf_adv%u(hgt_end%s(1:2))

!   write swell height

    write(11,'(i6,6es16.8)') step, step*deltat, &
                             max_height, min_height, end_height
    print '(i6,6es16.8)', step, step*deltat, &
                          max_height, min_height, end_height

  end subroutine solve_surface_height_corrector


  subroutine build_vpG

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, physqrow=[2,3], physqcol=[2,3], &
      coefficients=coefficients )

!   DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_elem, addmatvec=.true., &
      physqrow=[1,2], physqcol=[1,2], coefficients=coefficients )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[1], physqcol=[3], zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[3], physqcol=[1], zeromatvec=.true. )

  end subroutine build_vpG


! project c=exp(s) on discrete fem space

  subroutine solve_exps_projection

    type(solver_options_ma57_t) :: solver_options_ma57

    integer :: i, m

!   build vector only (matrix is constant)

    call build_system ( mesh, problemc_proj, sysmatrixc_proj, &
      m2sysvector=rhsc_proj, elemsub=exps_projection_elem, &
      oldvectors=oldvectors_ve, coefficients=coefficients )

    ! MA57 solver storage
    solver_options_ma57%integer_storage = 1.3
    solver_options_ma57%real_storage    = 1.3

!   LU decomposition is done in the first call only

    do m = 1, nmodes
      do i = 1, ncompc

        call add_effect_of_essential_to_rhs ( problemc_proj, sysmatrixc_proj, &
           solc_proj(i,m), rhsc_proj(i,m) )

        call solve_system_ma57 ( sysmatrixc_proj, rhsc_proj(i,m), &
           solc_proj(i,m), lu_exps_proj, solver_options=solver_options_ma57 )

      end do
    end do

    call delete ( lu_exps_proj )

  end subroutine solve_exps_projection

end program moving_boundary3a
