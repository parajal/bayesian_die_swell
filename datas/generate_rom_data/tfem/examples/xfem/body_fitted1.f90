! 2D Stokes problem on a square domain with a cylinder.
! Boundary (body) fitted mesh.

program body_fitted1

  use tfem_m
  use stokes_elements_m
  use hsl_ma57_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3             ! 3x3 Gauss integration

  real(dp), parameter :: &
    eta = 1._dp,    & ! viscosity
    flowrate = 10._dp  ! flowrate

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0
  coefficients%r(6) = flowrate

  call write_coefficients ( coefficients, filename='coefficients.out' )

! read mesh

  call read_mesh ( mesh, filename='mesh_bf.out' )

  call add_to_mesh ( mesh, curve=[-4] ) ! curve 14

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=13, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=12, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=7, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

! define constraints on curve

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=10, nglobalc=1 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=10, curve2=14, discretization='collocation', exclude=3 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

  sol%u = 0

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficients  )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

! write profile to a file for plotting with gnuplot

  call printtofile ( mesh, problem, 'vprofile_bf.out', curve=10, &
    vector=velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! derive vectors

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! write profile to a file for plotting with gnuplot

  call printtofile ( mesh, problem, 'pressure_bf.out', curve=13, &
    vector=pressure )

  call printtofile ( mesh, problem, 'vorticity_bf.out', curve=13, &
    vector=vorticity )

! write to a fig file for plotting

  !plot_options%numlevels=9
  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_bf.fig', vector=pressure )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_bf.fig', vector=vorticity )

! write binary file for reading by streamfunction

  open(unit=10,form='unformatted',file='velocity_bin_bf.out')

  write(10) velocity%u

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program body_fitted1
