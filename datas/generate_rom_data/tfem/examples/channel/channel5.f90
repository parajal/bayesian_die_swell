! 3D Stokes problem in a channel with a square cross section. Computed is the
! developed flow using only one layer of elements in the flow direction (x) and
! assuming periodical boundary conditions.
!
! Navier-slip boundary conditions.
! Dirichlet normal wall velocity.
! Collocation for periodic boundaries.

module elements_m

  use tfem_elem_m

  implicit none

contains

! Element for the constraints (connection through collocation)
! x- and y-direction only

  subroutine stokes_constr_node_conn1 ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


!   connection through collocation

    elemmat(1,:)  = [ 1._dp, 0._dp, 0._dp ]
    elemmat(2,:)  = [ 0._dp, 1._dp, 0._dp ]
    elemmat2 = - elemmat
    elemvec  = 0

  end subroutine stokes_constr_node_conn1

! Element for the constraints (connection through collocation)
! x- and z-direction only

  subroutine stokes_constr_node_conn2 ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


!   connection through collocation

    elemmat(1,:)  = [ 1._dp, 0._dp, 0._dp ]
    elemmat(2,:)  = [ 0._dp, 0._dp, 1._dp ]
    elemmat2 = - elemmat
    elemvec  = 0

  end subroutine stokes_constr_node_conn2

! Element for the constraints (connection through collocation) x-direction only

  subroutine stokes_constr_node_conn3 ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


!   connection through collocation

    elemmat(1,:)  = [ 1._dp, 0._dp, 0._dp ]
    elemmat2 = - elemmat
    elemvec  = 0

  end subroutine stokes_constr_node_conn3

end module elements_m


program channel5

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m
  use functions_m
  use elements_m

  implicit none

! constants

  integer, parameter :: &
    nx = 1,             & ! number of elements in x-direction
    ny = 10,            & ! number of elements in y-direction
    nz = 10,            & ! number of elements in z-direction
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3x3 integration of hexahedra
    vfuncnr = 1           ! function number for the wall velocity

  real(dp), parameter :: &
    lx = 0.1_dp,         & ! size in x-direction
    ly = 1._dp,          & ! size in y-direction
    lz = 1._dp,          & ! size in z-direction
    U = 1._dp,           & ! imposed average velocity
    eta = 1._dp,         & ! viscosity
    slip = 2._dp           ! slip coefficient

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure
  type(plot_options_t) :: plot_options
  type(meshgen_options_t) :: meshgen_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma57_t) :: solver_options

! variables

  real(dp) :: flowrate

  integer :: presnod(8) = [1,3,9,7,19,21,27,25]


  flowrate = U * ly * lz

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0
  coefficients%i(55) = -1  ! wall velocity = 0
  !coefficients%i(55) = 0  ! constant wall velocity
  !coefficients%i(55) = 1  ! function for wall velocity
  coefficients%i(56) = vfuncnr

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(3) = slip
  coefficients%r(6) = flowrate
  coefficients%r(20:22) = [ 0.9_dp, 0.0_dp, 0.0_dp ] ! constant wall velocity

  coefficients%vfunc1(2)%p => wall_velocity

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, &
    lz=lz, elshape=14, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! plot curves, surfaces and mesh

  plot_options%viewpoint=[1.,0.8,0.4]
  plot_options%fontsize=10
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig', surfaces=[3,4,6] )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(presnod,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, surface1=1, physq=1, &
    degfd=[0,0,1] )
  call define_essential ( mesh, input_probdef, surface1=2, physq=1, &
    degfd=[0,1,0] )
  call define_essential ( mesh, input_probdef, surface1=4, physq=1, &
    degfd=[0,1,0] )
  call define_essential ( mesh, input_probdef, surface1=6, physq=1, &
    degfd=[0,0,1] )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=1, surface1=5, nglobalc=1 )

! constraints for periodical boundary conditions

! velocities
  call define_constraint ( mesh, input_probdef, &
    physq=1, surface1=3, surface2=5, discretization='collocation', &
    excludecurves=[2,6,7,10] )

! velocities on curves
  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, curve2=4, discretization='collocation', nodedof=2, &
    exclude=3 )
  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=10, curve2=12, discretization='collocation', nodedof=2, &
    exclude=3 )
  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=6, curve2=5, discretization='collocation', nodedof=2, &
    exclude=3 )
  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=7, curve2=8, discretization='collocation', nodedof=2, &
    exclude=3 )

! velocity in points
  call define_constraint ( mesh, input_probdef, &
    physq=1, point1=2, point2=1, discretization='collocation', nodedof=1 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, point1=3, point2=4, discretization='collocation', nodedof=1 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, point1=7, point2=8, discretization='collocation', nodedof=1 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, point1=6, point2=5, discretization='collocation', nodedof=1 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

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
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=3, constraint2=4, elemsub=stokes_constr_node_conn1, &
    addmatvec=.true., coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=5, constraint2=6, elemsub=stokes_constr_node_conn2, &
    addmatvec=.true., coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=7, constraint2=10, elemsub=stokes_constr_node_conn3, &
    addmatvec=.true., coefficients=coefficients )

  call add_boundary_elements ( mesh, problem, rhsd, surface=1, &
    physq=[physqvel], sysmatrix=sysmatrix, elemsub=stokes_natboun_slip2, &
    coefficients=coefficients )

  call add_boundary_elements ( mesh, problem, rhsd, surface=2, &
    physq=[physqvel], sysmatrix=sysmatrix, elemsub=stokes_natboun_slip2, &
    coefficients=coefficients )

  call add_boundary_elements ( mesh, problem, rhsd, surface=4, &
    physq=[physqvel], sysmatrix=sysmatrix, elemsub=stokes_natboun_slip2, &
    coefficients=coefficients )

  call add_boundary_elements ( mesh, problem, rhsd, surface=6, &
    physq=[physqvel], sysmatrix=sysmatrix, elemsub=stokes_natboun_slip2, &
    coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%real_storage=1.2

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a fig file for plotting

  plot_options%printlabels=.false.

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity, surfaces=[3] )
  call plot_points_curves ( plot_options, mesh, 'velocity.fig', append=.true. )

  call plot_color_fill ( plot_options, mesh, problem, 'vx_color.fig', &
    vector=velocity, degfd=1, surfaces=[3] )
  call plot_points_curves ( plot_options, mesh, 'vx_color.fig',  &
   append=.true. )

  call plot_color_fill ( plot_options, mesh, problem, 'vy_color.fig', &
    vector=velocity, degfd=2, surfaces=[3] )
  call plot_points_curves ( plot_options, mesh, 'vy_color.fig',  &
   append=.true. )

  call plot_color_fill ( plot_options, mesh, problem, 'vz_color.fig', &
    vector=velocity, degfd=3, surfaces=[3] )
  call plot_points_curves ( plot_options, mesh, 'vz_color.fig',  &
   append=.true. )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure, surfaces=[3,4,6] )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='channel.vtk' )

  call write_vector_vtk ( mesh, problem, filename='channel.vtk', &
    dataname='velocity_vector', sysvector=sol, append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( velocity, pressure )
  call delete ( coefficients )
  call delete ( oldvectors )

end program channel5
