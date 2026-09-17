! Axisymmetrical Stokes problem: a sphere moving with velocity U on the center
! line of a cylindrical container with a bottom (at infinity). The sphere is
! fixed and the cylinder wall is moved with velocity -U (moving frame).
! At the ends fluid is moving with velocity -U as well, representing zero
! flux in the stationary frame (bottom).
! Dirichlet boundary conditions for the velocities are used.
! High-order Taylor-Hood quadrilateral elements using blended meshes.
! This problem is similar to sphere1, but now with high-order Taylor-Hood
! velocity/pressure Qp/Qp-1 elements
! Blended mesh:
!  main mesh: (p+1)^2-node high-order quadrilaterals, equidistant nodes
!  blend mesh 1: p^2-node high-order quadrilaterals, equidistant nodes

program sphere32

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 20,    & ! Qp velocities
    pintpl = 20,    & ! Qp-1 pressures
    p = 4,          & ! polynomial order
    q = p-1,        & ! polynomial order pressure
    n0 = (p+1)**2,  & ! number of nodes in main mesh element
    n1 = (q+1)**2,  & ! number of nodes in blend mesh element
    physqvel = 1,   & ! physical quantity nr of the velocities
    physqpress = 2, & ! physical quantity nr of the pressures
    gauss = p+1,    & ! number of Gauss points (1D)
    gaussb = p+1,   & ! number of Gauss points on the boundary
    inttype = 3,    & ! numerical rules for integration
    coorsys = 1       ! axisymmetric coordinate system

! definitions

  type(mesh_t) :: mesh, mesh_Qp, mesh_Qpm1
  type(mesh_t) :: mesh_plot_Qp, mesh_plot_Qpm1
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity, divergence
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

! variables

  integer :: crv

  real(dp), parameter :: &
    eta = 1._dp,  & ! viscosity
    U = 1._dp       ! velocity of the sphere

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gaussb ]
  coefficients%i(12:) = 0
  coefficients%i(23) = coorsys

  coefficients%i(40) = inttype
  coefficients%i(51) = p ! velocity polynomial order
  coefficients%i(52) = q ! pressure polynomial order

  coefficients%r(1) = eta
  coefficients%r(3:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

! read mesh

  call read_mesh_gmsh ( mesh_Qp, filename='mesh32_Qp.msh', ndim=2 )
  call read_mesh_gmsh ( mesh_Qpm1, filename='mesh32_Qpm1.msh', ndim=2 )

! curve for printing data on the centerline+sphere

  call add_to_mesh ( mesh_Qp, curve=[1,2,3,4,5,6,7,8], newnr=crv )
  call add_to_mesh ( mesh_Qpm1, curve=[1,2,3,4,5,6,7,8], newnr=crv )

  call mesh_convert ( mesh_Qp, blendmesh=mesh_Qpm1, mesh=mesh )

  call fill_mesh_parts ( mesh )

! plot points/curves for main and blend mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, filename='curves_main32.fig' )
  call plot_points_curves ( plot_options, mesh, blend=1, &
    filename='curves_blend32.fig' )

! Create plot meshes

  call mesh_convert ( mesh_Qp, mesh_plot_Qp, warn=.false., &
    elementshapes='spectraltolinear' )
  call mesh_convert ( mesh_Qpm1, mesh_plot_Qpm1, warn=.false., &
    elementshapes='spectraltolinear' )

  call fill_mesh_parts ( mesh_plot_Qp )
  call fill_mesh_parts ( mesh_plot_Qpm1 )

  call delete ( mesh_Qp, mesh_Qpm1 )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

  input_probdef%vec_elementdof(1)%a = 0           ! initialize to zero
  input_probdef%vec_elementdof(1)%a(1:n0,1) = 2   ! velocity (Qp mesh)
  input_probdef%vec_elementdof(1)%a(n0+1:n0+n1,2) = 1  ! pressure (Qpm1 mesh)
  input_probdef%vec_elementdof(1)%a(1:n0,3) = 1    ! scalar (Qp mesh)
  input_probdef%vec_elementdof(1)%a(n0+1:n0+n1,4) = 1  ! scalar (Qpm1 mesh)

  input_probdef%physq = [1,2]

! define essential boundaries

! center line
  call define_essential ( mesh, input_probdef, curve1=1, curve2=2, physq=1, &
    degfd=[0,1] )
  call define_essential ( mesh, input_probdef, curve1=7, curve2=8, physq=1, &
    degfd=[0,1] )
! sphere
  call define_essential ( mesh, input_probdef, curve1=3, curve2=6, physq=1 )
! ends and cylinder wall
  call define_essential ( mesh, input_probdef, curve1=9, curve2=14, physq=1 )
! pressure level
  call define_essential ( mesh, input_probdef, point=15, physq=2 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0

  call fill_sysvector ( mesh, problem, sol, &
    curve1=9, curve2=14, physq=1, degfd=1, value=-U )
  call fill_sysvector ( mesh, problem, sol, &
    point=15, physq=2, value=0._dp )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call delete ( sysmatrix )

! write solution to a vtk file for further post-processing

  call create_vector ( problem, velocity, physq=1 )
  call extract_physvector ( mesh, problem, sol, velocity )
  call create_vector ( problem, pressure, physq=2 )
  call extract_physvector ( mesh, problem, sol, pressure )

  call write_vector_vtk ( mesh_plot_Qp, problem, vector=velocity, &
    dataname='velocity', filename='sphere32.vtk' )

  call write_scalar_vtk ( mesh_plot_Qpm1, problem, vector=pressure, &
    dataname='pressure', filename='pressure32.vtk', &
    offset=mesh%nnodes_blend(2) )

! print solution only on center line

  call printtofile ( mesh, problem, filename='sol_on_center.out', &
    curve=crv, sysvector=sol )

! print velocity only on center line

  call printtofile ( mesh, problem, filename='velocity_on_center.out', &
    curve=crv, vector=velocity )

! print pressure only on center line

  call printtofile ( mesh, problem, filename='pressure_on_center.out', &
    curve=crv, vector=pressure )

! derivatives

  call delete ( pressure )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )
  call create_vector ( problem, divergence, vec=3 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5
  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=7
  call derive_vector ( mesh, problem, divergence, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! write derivatives to a vtk file for further post-processing

  call write_scalar_vtk ( mesh_plot_Qp, problem, vector=pressure, &
    dataname='pressure', filename='sphere32.vtk', append=.true. )

  call write_scalar_vtk ( mesh_plot_Qp, problem, vector=vorticity, &
    dataname='vorticity', filename='sphere32.vtk', append=.true. )

  call write_scalar_vtk ( mesh_plot_Qp, problem, vector=divergence, &
    dataname='divergence', filename='sphere32.vtk', append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity, divergence )
  call delete ( coefficients )
  call delete ( oldvectors )

end program sphere32
