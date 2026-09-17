! Steady Stokes flow and temperature problem with viscous dissipation.
! Periodic pipe flow. Constant coefficients.
! Error compared to analytical solution for velocity and temperature.

program stress_work1

  use tfem_m
  use math_defs_m
  use stokes_elements_m
  use generalized_stokes_elements_m, only: fill_viscous_dissipation_gauss
  use convection_diffusion_elements_m
  use stress_work_functions_m
  use hsl_ma57_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,         & ! velocity interpolation
    pintpl = 2,         & ! pressure interpolation
    tintpl = uintpl,    & ! temperature interpolation
    physqv = 1,         & ! physical quantity nr of the velocity
    physqp = 2,         & ! physical quantity nr of the pressure
    gauss = 6,          & ! Gauss integration on volumes
    gaussb = 3,         & ! Gauss integration on curves
    coorsys = 1,        & ! axisymmetric coordinate system
    nz = 2,             & ! number of elements in z-direction
    nr = 10               ! number of elements in r-direction

  real(dp), parameter :: &
    eta = 10._dp,        & ! viscosity
    kappa = 5.0_dp,      & ! thermal conductivity coefficient
    Uavg = 2.0_dp,       & ! imposed averaged velocity
    T0 = 0.5_dp,         & ! imposed wall temperature
    lz = 0.6_dp,         & ! length of the domain in z-direction
    R = 3._dp,           & ! radius of the cylinder
    rs_up = 1.5_dp,      & ! real_storage gradient-velocity-pressure LU (HSL)
    is_up = 1.0_dp,      & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_t = 1.0_dp,       & ! real_storage temperature LU (HSL)
    is_t = 1.0_dp          ! integer_storage temperature LU (HSL)

! variables

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdeft
  type(problem_t) :: problem, problemt
  type(sysmatrix_t) :: sysmatrix, sysmatrixt
  type(sysvector_t), target :: sol
  type(sysvector_t) :: sol_exact, rhsd
  type(sysvector_t), target :: solt
  type(sysvector_t) :: solt_exact, rhsdt
  type(subscript_t) :: svel1, sT
  type(coefficients_t) :: coefficients, coefficientst
  type(oldvectors_t) :: oldvectors
  type(solver_options_ma57_t) :: solver_options, solver_optionst
  type(elvector_t), target :: sw_source

! Set function parameters

  Uavg_l = Uavg
  R_l = R
  eta_l = eta
  kappa_l = kappa
  T0_l = T0


! fill coefficients flow problem

  call create_coefficients ( coefficients, ncoefi=300, ncoefr=250 )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl, pintpl, 0,     0,         0,  &
      physqv, physqp, 0,     0,     gauss,  &
      gaussb ]
  coefficients%i(23)  = coorsys

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(6) = Uavg * pi * R**2


! fill coefficients temperature problem

  call create_coefficients ( coefficientst, ncoefi=100, ncoefr=50 )

  coefficientst%i = 0
  coefficientst%i(1) = tintpl
  coefficientst%i(10:11) = [ gauss, gaussb ]
  coefficientst%i(18) = 3 ! right-hand side given per element in int. points
  coefficientst%i(23)  = coorsys

  coefficientst%r(1) = kappa
  coefficientst%r(2:) = 0


! create mesh

  meshgen_options%elshape = 4 ! 6-node triangles
  meshgen_options%nx = nz
  meshgen_options%ny = nr
  meshgen_options%lx = lz
  meshgen_options%ly = R

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] ) ! curve 5

  call fill_mesh_parts ( mesh )


! problem definition for the flow problem

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,   &  ! velocity
                 1,0,1,0,1,0,   &  ! pressure
                 1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                   [6,3] )

  input_probdef%physq = [1,2]
  input_probdef%probnr = 1

  call define_essential ( mesh, input_probdef, curve1=1, physq=physqv, &
    degfd=[0,1], exclude=2 )
  call define_essential ( mesh, input_probdef, curve1=3, physq=physqv )
  call define_essential ( mesh, input_probdef, point=1, physq=physqp )

  call define_constraint ( mesh, input_probdef, physq=physqv, curve1=2, &
    nglobalc=1 )

  call define_constraint ( mesh, input_probdef, &
    physq=physqv, curve1=2, curve2=5, discretization='collocation', exclude=2 )

  call problem_definition ( input_probdef, mesh, problem )


! problem definition for the temperature problem

  call create_input_probdef ( mesh, input_probdeft, nvec=2, nphysq=1 )

  input_probdeft%vec_elementdof(1)%a =   &
      reshape ( [1,1,1,1,1,1,   &  ! temperature
                 1,1,1,1,1,1 ], &  ! scalar in in all nodes
                   [6,2] )

  input_probdeft%physq = [1]
  input_probdeft%probnr = 2

  call define_essential ( mesh, input_probdeft, curve1=3 )

  call define_constraint ( mesh, input_probdeft, &
    curve1=2, curve2=5, discretization='collocation', exclude=2 )

  call problem_definition ( input_probdeft, mesh, problemt )


! create system vectors (solution and right-hand side)

  call create ( problem, sol, rhsd )
  call create ( problemt, solt, rhsdt )

! create vector defined per element for scalar values of right-hand

  call create ( mesh, sw_source, nreal1d=1 )

! oldvectors

  call create ( oldvectors, nsysvec=1, nelvec=2 )
  oldvectors%s(1)%p => sol
  oldvectors%e(2)%p => sw_source

! set solution vectors to zero (including ess bc)

  sol%u = 0
  solt%u = 0


! create system matrix flow problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! solve velocity field

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficients  )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%real_storage=rs_up
  solver_options%integer_storage=is_up

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

! compute viscous dissipation

  call loop_over_elements ( mesh, problem, &
    elemsub=fill_viscous_dissipation_gauss, &
    coefficients=coefficients, oldvectors=oldvectors )


! create system matrix temperature problem

  call create_sysmatrix_structure_base ( sysmatrixt, mesh, problemt, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrixt, mesh, problemt )
  call finalize_sysmatrix_structure ( sysmatrixt )

  call create_sysmatrix_data ( sysmatrixt )

! solve temperature field

  call fill_sysvector ( mesh, problemt, solt, curve1=3, value=T0 )

  call build_system ( mesh, problemt, sysmatrixt, rhsdt, &
    elemsub=scalar_diffusion_elem, coefficients=coefficientst, &
    oldvectors=oldvectors )

  call build_system_constraint ( mesh, problemt, sysmatrixt, rhsdt, &
    constraint1=1, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficientst  )

  call add_effect_of_essential_to_rhs ( problemt, sysmatrixt, solt, rhsdt )

  solver_optionst%real_storage=rs_t
  solver_optionst%integer_storage=is_t

  call solve_system_ma57 ( sysmatrixt, rhsdt, solt, &
    solver_options=solver_optionst )


! postprocessing


! create exact solution for velocity

  call create_sysvector ( problem, sol_exact )

  call fill_sysvector ( mesh, problem, sol_exact, &
    node1=1, node2=mesh%nnodes, physq=physqv, degfd=1, func=func, funcnr=3 )

  call create_subscript ( mesh, problem, svel1, physqarr=[physqv], degfd=1 )

! print maximum difference of sol-solexact to standard output

  print *, 'max u_z error', maxval( abs(sol%u(svel1%s)-sol_exact%u(svel1%s)) )


! create exact solution for temperature

  call create_sysvector ( problemt, solt_exact )

  call fill_sysvector ( mesh, problemt, solt_exact, &
    node1=1, node2=mesh%nnodes, physq=1, func=func, funcnr=4 )

  call create_subscript ( mesh, problemt, st, physqarr=[1] )

! print maximum difference of solt-solt_exact to standard output

  print *, 'max T error', maxval( abs(solt%u(st%s)-solt_exact%u(st%s)) )


! write vtk

  call write_vector_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='velocity', sysvector=sol, physq=physqv )

  call write_scalar_vtk ( mesh, problemt, filename='sol.vtk', &
    dataname='temperature', sysvector=solt, append=.true. )


! delete all data including all allocated memory

  call delete ( problem, problemt )
  call delete ( input_probdef, input_probdeft )
  call delete ( mesh )
  call delete ( sol, solt, sol_exact, solt_exact )
  call delete ( rhsd, rhsdt )
  call delete ( sw_source )
  call delete ( sysmatrix, sysmatrixt )
  call delete ( coefficients )
  call delete ( oldvectors )

end program stress_work1
