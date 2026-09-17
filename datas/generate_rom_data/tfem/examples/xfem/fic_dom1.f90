! 2D Stokes problem on a square domain with a cylinder
! Fictitious domain approach with weak constraint on the cylinder.
!
program fic_dom1

  use tfem_m
  use stokes_elements_m
  use hsl_ma57_m
  use io_utils_m
  use figplot_m
  use subs1_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    nx=80,              & ! number of elements in x
    ny=80,              & ! number of elements in y
    nsubinto = 30,      & ! number of subintegration points on object
    gauss = 3             ! 3x3 Gauss integration

  real(dp), parameter :: &
    lx = 8._dp,        & ! width of domain
    ly = 8._dp,        & ! height of domain
    radius = 1._dp,    & ! radius of the cylinder
    center(2) = [ 0._dp, 0._dp ], & ! center position of the cylinder
    eta = 1._dp,       & ! viscosity
    flowrate = 10._dp  ! flowrate

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh_particle
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(sample_t) :: sample_p, sample_v, sample_vort

  integer :: i, elem, nelem

! module timer_m

  timer = .false. ! set to .true. to show cpu time output

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=300, ncoefr=250 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0
  coefficients%r(6) = flowrate

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%lx = lx
  meshgen_options%ly = ly
  meshgen_options%ox = -lx/2
  meshgen_options%oy = -ly/2

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] ) ! curve 5

! create mesh for particle boundary

  nelem = 4 * nx * radius / lx

  call mesh_skeleton ( mesh_particle, nnodes=2*nelem, nelem=nelem, elshape=2, &
    ndim=2 )

  rp = radius ! radius of the circular object
  xpc = center ! initial position of the center of the object

  call objectscoor ( 1, mesh_particle%coor )

  do elem = 1, mesh_particle%nelem
    mesh_particle%topology(1)%a(:,elem) = [ 2*elem-1, 2*elem, 2*elem + 1 ]
  end do
  mesh_particle%topology(1)%a(3,mesh_particle%nelem) = 1 ! close circle

! one object

  call add_to_mesh ( mesh, object='mesh', objectmesh=mesh_particle, &
    topology=.true., intrule=1, nsubint=nsubinto )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves_fd.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh_fd.fig' )
  plot_options%objectpointcolor=4
  plot_options%objectpointsize=0.4
  call plot_objects ( plot_options, mesh, 'mesh_fd.fig', append=.true. )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=3, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

! define constraints on the object

  call define_constraint ( mesh, input_probdef, object=1, physq=1, &
    discretization='weak', elementdof=[2,0,2] )

! define constraints on curve

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, nglobalc=1 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, curve2=5, discretization='collocation', exclude=3 )

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

  call tic

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients, buildvector=.false. )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=elementc, addmat=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=3, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficients  )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call toc ( 'build' )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call toc ( 'solve' )


! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

! write profile to a file for plotting with gnuplot

  call printtofile ( mesh, problem, 'vprofile_fd.out', curve=2, &
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

  open( unit=15, file='sample_p_fd.out' )

! sample pressure on cylinder

  call fill_sample ( mesh, problem, sample_p, ndegfd=1, object=1, &
    elemsub=stokes_sample_pressure, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_p%nnodes
    write( unit=15, fmt=* ) sample_p%coor(i,:), sample_p%u(i,:)
  end do

  close( unit=15 )

  open( unit=15, file='sample_v_fd.out', recl=300 )

! sample velocity on cylinder

  call fill_sample ( mesh, problem, sample_v, ndegfd=2, object=1, &
    elemsub=stokes_sample_velocity, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_v%nnodes
    write( unit=15, fmt=* ) sample_v%coor(i,:), sample_v%u(i,:)
  end do

  close( unit=15 )

  open( unit=15, file='sample_vort_fd.out', recl=300 )

! sample vorticity on cylinder

  coefficients%i(13)=5

  call fill_sample ( mesh, problem, sample_vort, ndegfd=1, object=1, &
    elemsub=stokes_sample_deriv, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_v%nnodes
    write( unit=15, fmt=* ) sample_vort%coor(i,:), sample_vort%u(i,:)
  end do

  close( unit=15 )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_fd.fig', vector=pressure )
  !call plot_objects ( plot_options, mesh, 'pressure_fd.fig', append=.true. )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_fd.fig', vector=vorticity )
  !call plot_objects ( plot_options, mesh, 'vorticity_fd.fig', append=.true. )

! write binary file for reading by streamfunction

  open(unit=10,form='unformatted',file='velocity_bin_fd.out')

  write(10) velocity%u

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh, mesh_particle )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( sample_p, sample_vort, sample_v )

end program fic_dom1
