! Unsteady energy balance equation (including dissipation)
!
!   gamma ( dT/dt + u . nabla T ) - alpha nabla^2 T = f
!
! gamma = rho . cp
! alpha = kappa
! f = tau : D
!
! Newtonian fluid with temperature dependent eta
! Flow computed using a flow constraint
! The initial fluid and the fluid flowing in is at temperature t0
! Zero-flux boundary conditions are assumed on the boundaries
! The fluid will heat up due to the dissipation term

program energy1

  use tfem_m
  use hsl_ma41_m
  use hsl_ma57_m
  use generalized_stokes_elements_m
  use convection_diffusion_supg_elements_m
  use io_utils_m
  use subs_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,         & ! velocity interpolation
    pintpl = 2,         & ! pressure interpolation
    tintpl = 2,         & ! temperature interpolation
    physqvel = 1,       & ! physical quantity nr of the velocity
    physqpress = 2,     & ! physical quantity nr of the pressure
    timeint1 = 1,       & ! first-order, semi-implicit Euler
    timeint2 = 2,       & ! second-order, semi-implicit Gear
    gauss = 6,          & ! 6 point integration of triangles
    gaussb = 4,         & ! 4 point integration of lines
    coorsys = 1,        & ! axisymmetric coordinate system
    nx = 20,            & ! number of elements in x-direction
    ny = 10,            & ! number of elements in y-direction
    ibeta = 2,          & ! method to compute beta in SUPG
    numtimesteps=100      ! number of time steps

  real(dp), parameter :: &
    eta0 = 10._dp,       & ! maximum viscosity
    alpha = 1._dp,       & ! thermal conductivity coefficient
    gamma = 1._dp,       & ! density times heat capacity (at constant pressure)
    flowrate = 1.6_dp,   & ! flowrate of the inflow
    t0 = 0._dp,          & ! initial and inflow temperature
    deltat = 0.02_dp,    & ! time step
    ox = 0._dp,          & ! x-coordinate of origin of the mesh
    oy = 0._dp,          & ! y-coordinate of origin of the mesh
    lx = 2._dp,          & ! length of the domain in x-direction
    ly = 1._dp,          & ! length of the domain in y-direction
    rs = 1.6_dp,         & ! real_storage velocity-pressure LU (HSL)
    is = 1.6_dp            ! integer_storage velocity-pressure LU (HSL)

! variables

  integer :: step

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdeft
  type(problem_t), target :: problem, problemt
  type(sysmatrix_t) :: sysmatrix, sysmatrixt
  type(sysvector_t), target :: sol, sol_n, sol_nm1, solhat, rhsd
  type(sysvector_t), target :: solt, solt_n, solt_nm1, rhsdt, solthat
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors
  type(solver_options_ma41_t) :: solver_optionst
  type(solver_options_ma57_t) :: solver_options
  type(vector_t) :: viscosity_post
  type(elvector_t), target :: viscosity, energy_source
  type(subscript_t) :: tval

  character(len=199) :: filename

! pass variables to module

  leta0 = eta0
  lt0   = t0

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=450, ncoefr=400 )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gaussb ]
  coefficients%i(23)  = coorsys
  coefficients%i(401) = timeint1
  coefficients%i(402) = tintpl
  coefficients%i(403) = 1 ! 0: Galerkin, 1: SUPG
  coefficients%i(404) = 1 ! include diffusion term
  coefficients%i(405) = 0 ! no slef-dependent source term
  coefficients%i(406) = 4 ! f given in integration points given by the vector
                          ! per element oldvectors%e(2)
  coefficients%i(411) = ibeta  ! method to compute beta in SUPG
  coefficients%i(251) = 3 ! eta is a function of the temperature  given by the
                          ! vector per element oldvectors%e(1)

  coefficients%r = 0
  coefficients%r(1)   = eta0
  coefficients%r(6)   = flowrate
  coefficients%r(351:353) = [ deltat, gamma, alpha ]

! create mesh

  meshgen_options%elshape = 4 ! 6-node triangles
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%ox = ox
  meshgen_options%oy = oy
  meshgen_options%lx = lx
  meshgen_options%ly = ly

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] ) ! curve 5

  call fill_mesh_parts ( mesh )

! problem definition for the flow problem

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,    &  ! pressure
                 1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [6,3] )

  input_probdef%physq = [1,2]
  input_probdef%probnr = 1

  call define_essential ( mesh, input_probdef, curve1=1, curve2=2, &
    physq=physqvel, degfd=[0,1] )
  call define_essential ( mesh, input_probdef, curve1=3, physq=physqvel, &
    degfd=[1,1] )
  call define_essential ( mesh, input_probdef, curve1=4, physq=physqvel, &
    degfd=[0,1] )

  call define_constraint ( mesh, input_probdef, physq=physqvel, curve1=5, &
    nglobalc=1 )

  call problem_definition ( input_probdef, mesh, problem )

! problem definition for the temperature problem

  call create_input_probdef ( mesh, input_probdeft, nvec=2, nphysq=1 )

  input_probdeft%vec_elementdof(1)%a =   &
      reshape ( [1,0,1,0,1,0,     &  ! temperature
                 1,1,1,1,1,1 ],  &  ! scalar in in all nodes
                 [6,2] )

  input_probdeft%physq = [1]
  input_probdeft%probnr = 2

  call define_essential ( mesh, input_probdeft, curve1=4 )

  call problem_definition ( input_probdeft, mesh, problemt )

! create system vectors (solution and right-hand side)

  call create ( problem, sol, sol_n, sol_nm1, solhat, rhsd )
  call create ( problemt, solt, solt_n, solt_nm1, solthat, rhsdt )

! create vector defined per element for eta

  call create ( mesh, viscosity, nreal1d=1 )

! create vector defined per element for scalar values of right-hand

  call create ( mesh, energy_source, nreal1d=1 )

! oldvectors

  call create ( oldvectors, nprob=2, nsysvec=5, nelvec=2, nvec=1 )
  oldvectors%p(1)%p => problem
  oldvectors%p(2)%p => problemt
  oldvectors%s(1)%p => solhat
  oldvectors%s(2)%p => solt_n
  oldvectors%s(3)%p => solt_nm1
  oldvectors%s(4)%p => solthat
  oldvectors%s(5)%p => solt
  oldvectors%e(1)%p => viscosity
  oldvectors%e(2)%p => energy_source

! initialize solution vectors

  sol%u = 0
  solt%u = 0

! fill solution vector with essential boundary conditions for velocity field

  call fill_sysvector ( mesh, problem, sol, curve1=1, curve2=2, &
    physq=physqvel, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, curve1=3, physq=physqvel, &
    value=0._dp )
  call fill_sysvector ( mesh, problem, sol, curve1=4, physq=physqvel, &
    degfd=2, value=0._dp )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! fill initial temperature field

  call fill_sysvector ( mesh, problemt, solt, node1=1, node2=mesh%nnodes, &
    value=t0 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrixt, mesh, problemt )

  call create_sysmatrix_data ( sysmatrixt )

! solve initial velocity field

  call loop_over_elements ( mesh, problem, elemsub=fill_viscosity_gauss, &
    coefficients=coefficients, oldvectors=oldvectors )

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=generalized_stokes_elem, oldvectors=oldvectors, &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%real_storage=rs
  solver_options%integer_storage=is

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

! postprocessing

  call create ( problem, viscosity_post, vec=3 )
  call create ( mesh, problemt, tval, physqarr=[1] )

  step = 0

  call postprocessing

! fill temperature boundary conditions

  call fill_sysvector ( mesh, problemt, solt, curve1=4, value=t0 )

! start time stepping

  call copy ( sol, sol_n )
  call copy ( solt, solt_n )

  do step = 1, numtimesteps

    if ( step == 2 ) coefficients%i(401) = timeint2 ! second-order

!   predict several fields that are needed in the temperature problem

    if ( step == 1 ) then ! first-order

      solhat%u  = sol_n%u
      solthat%u = solt_n%u

    else ! second-order
      solhat%u  = 2*sol_n%u - sol_nm1%u
      solthat%u = 2*solt_n%u - solt_nm1%u

    end if

!   build diffusion matrix and right-hand side

    call loop_over_elements ( mesh, problem, elemsub=fill_viscosity_gauss, &
      coefficients=coefficients, oldvectors=oldvectors )

    call loop_over_elements ( mesh, problem, &
      elemsub=fill_viscous_dissipation_gauss, &
      coefficients=coefficients, oldvectors=oldvectors )

    call build_system ( mesh, problemt, sysmatrixt, rhsdt, &
      elemsub=conv_diff_supg_elem, coefficients=coefficients, &
      oldvectors=oldvectors )

    call add_effect_of_essential_to_rhs ( problemt, sysmatrixt, solt, rhsdt )

    solver_optionst%real_storage=rs
    solver_optionst%integer_storage=is

    call solve_system_ma41 ( sysmatrixt, rhsdt, solt, &
      solver_options=solver_optionst )

!   build velocity-pressure matrix and vector

    call loop_over_elements ( mesh, problem, elemsub=fill_viscosity_gauss, &
      coefficients=coefficients, oldvectors=oldvectors )

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=generalized_stokes_elem, oldvectors=oldvectors, &
      coefficients=coefficients )

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    solver_options%real_storage=rs
    solver_options%integer_storage=is

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options )

    write(*,'(a,i0,a,es12.4,a,es16.8/)') &
      'step = ', step, ' time = ', step*deltat

!   postprocessing

    call postprocessing

!   copy old values

    call copy ( solt_n, solt_nm1 )
    call copy ( solt, solt_n )
    call copy ( sol_n, sol_nm1 )
    call copy ( sol, sol_n )

  end do


! delete all data including all allocated memory

  call delete ( problem, problemt )
  call delete ( input_probdef, input_probdeft )
  call delete ( mesh )
  call delete ( sol, sol_n, sol_nm1, solhat )
  call delete ( solt, solt_n, solt_nm1 )
  call delete ( rhsd, rhsdt )
  call delete ( sysmatrix, sysmatrixt )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( viscosity_post )
  call delete ( viscosity )

contains

! subroutine to write vtks

  subroutine postprocessing

    type(oldvectors_t) :: oldvectors_post
    type(vector_t) :: temperature

    call create ( oldvectors_post, nsysvec=1 )
    oldvectors_post%s(1)%p => solt

!   derive temperature in all nodes

    call create_vector ( problemt, temperature, vec=2 )

    call derive_vector ( mesh, problemt, temperature, &
      elemsub=derive_q, coefficients=coefficients, &
      oldvectors=oldvectors_post )

!   compute local viscosity
    viscosity_post%u = compute_viscosity ( temperature%u )

    write(filename,'(a,i4.4,a)') 'sol', step, '.vtk'
    call write_scalar_vtk ( mesh, problemt, vector=temperature, &
      filename=filename, dataname='temperature' )

    write(filename,'(a,i4.4,a)') 'velo', step, '.vtk'
    call write_vector_vtk ( mesh, problem, filename=filename, &
      dataname='velocity', sysvector=sol, physq=physqvel )

    call write_scalar_vtk ( mesh, problem, filename=filename, &
      dataname='viscosity', vector=viscosity_post, append=.true. )

  end subroutine postprocessing

end program energy1
