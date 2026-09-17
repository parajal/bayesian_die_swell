! Axisymmetrical Stokes problem: a sphere moving with velocity U on the center
! line of a cylindrical container with a bottom (at infinity). The sphere is
! fixed and the cylinder wall is moved with velocity -U (moving frame).
! At the ends fluid is moving with velocity -U as well, representing zero
! flux in the stationary frame (bottom).
! Dirichlet boundary conditions for the velocities are used.
! Equal order high-order triangular elements.
! This problem is similar to sphere1, but now with high-order elements having
! equal-order continuous stabilized velocity/pressure Pp/Pp elements
!

program sphere30

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 19,    & ! Pp velocities
    pintpl = 19,    & ! Pp pressures
    pintpl1 = 21,   & ! Pp-1 pressure projection space
    p = 4,          & ! polynomial order
    q = p,          & ! polynomial order pressure
    q1 = p-1,       & ! polynomial order pressure projection space
    physqvel = 1,   & ! physical quantity nr of the velocities
    physqpress = 2, & ! physical quantity nr of the pressures
    gauss = 2*p,    & ! number of Gauss points
    gaussb = p+1,   & ! number of Gauss points on the boundary
    inttype = 3,    & ! numerical rules for integration
    coorsys = 1       ! axisymmetric coordinate system

! definitions

  type(mesh_t) :: mesh, mesh_plot
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
  coefficients%i(53) = pintpl1 ! pressure stabilization space
  coefficients%i(54) = q1 ! pressure stabilization polynomial order

  coefficients%r(1) = eta
  coefficients%r(2) = eta ! pressure stabilization parameter
  coefficients%r(3:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

! read mesh

  call read_mesh_gmsh ( mesh, filename='mesh30.msh', ndim=2 )

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

! curve for printing data on the centerline+sphere

  call add_to_mesh ( mesh, curve=[1,2,3,4,5,6,7,8], newnr=crv )

  call fill_mesh_parts ( mesh )

! Create plot mesh

  call mesh_convert ( mesh, mesh_plot, elementshapes='spectraltolinear', &
    warn=.false. )

! write mesh_plot (read by streamfunction computation)

  call write_mesh ( mesh_plot, filename='mesh_plot.out' )

  call fill_mesh_parts( mesh_plot )

! plot mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh_plot, filename='curves.fig' )
  plot_options%fontsize = 10
  call plot_mesh ( plot_options, mesh_plot, filename='mesh_plot.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 2  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar

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
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0

  call fill_sysvector ( mesh, problem, sol, &
    curve1=9, curve2=14, physq=1, degfd=1, value=-U )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem_stab, physqrow=[physqpress], physqcol=[physqpress], &
    coefficients=coefficients, buildvector=.false., addmatvec=.true. )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call delete ( sysmatrix )

  call create_vector ( problem, velocity, physq=1 )
  call extract_physvector ( mesh, problem, sol, velocity )

! post processing

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

! write to a fig file for plotting

  plot_options%xmin = -6
  plot_options%shiftbary=20

  call plot_vector ( plot_options, mesh_plot, problem, 'velocity.fig', &
    vector=velocity )

  call plot_color_contour ( plot_options, mesh_plot, problem, 'pressure.fig', &
    vector=pressure )

  plot_options%rotatecolorbar=.true.
  plot_options%shiftbary=100

  call plot_color_fill ( plot_options, mesh_plot, problem, &
    'velocity_color.fig', vector=velocity, degfd=1 )

  call plot_color_fill ( plot_options, mesh_plot, problem, &
    'vorticity_color.fig', vector=vorticity )

  call plot_color_contour ( plot_options, mesh_plot, problem, &
    'divergence.fig', vector=divergence )

  call printtofile ( mesh_plot, problem, filename='velocity_on_center.out', &
    curve=crv, vector=velocity )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh_plot, problem, vector=pressure, &
    dataname='pressure', filename='sphere30.vtk' )

  call write_vector_vtk ( mesh_plot, problem, filename='sphere30.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true. )

  call write_scalar_vtk ( mesh_plot, problem, vector=vorticity, &
    dataname='vorticity', filename='sphere30.vtk', append=.true. )

  call write_scalar_vtk ( mesh_plot, problem, vector=divergence, &
    dataname='divergence', filename='sphere30.vtk', append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh, mesh_plot )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity, divergence )
  call delete ( coefficients )
  call delete ( oldvectors )

end program sphere30
