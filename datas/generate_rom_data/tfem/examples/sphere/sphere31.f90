! Axisymmetrical Stokes problem: a sphere moving with velocity U on the center
! line of a cylindrical container with a bottom (at infinity). The sphere is
! fixed and the cylinder wall is moved with velocity -U (moving frame).
! At the ends fluid is moving with velocity -U as well, representing zero
! flux in the stationary frame (bottom).
! Dirichlet boundary conditions for the velocities are used.
! Mixed meshes with both biquadratic and bilinear quadrilateral elements.
! This problem is similar to sphere1, but now with blended meshes.
! Blended mesh:
!  main mesh: 9-node quadrilaterals
!  blend mesh 1: 4-node quadrilaterals

program sphere31

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities
    pintpl = 4,     & ! Q1 pressures
    physqvel = 1,   & ! physical quantity nr of the velocities
    physqpress = 2, & ! physical quantity nr of the pressures
    gauss = 3,      & ! 3x3 integration of quads
    gaussb = 3,     & ! number of Gauss points on the boundary
    inttype = 3,    & ! numerical rules for integration
    coorsys = 1       ! axisymmetric coordinate system

! definitions

  type(mesh_t) :: mesh, mesh_Q2, mesh_Q1
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

  coefficients%r(1) = eta
  coefficients%r(3:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

! read mesh

  call read_mesh_gmsh ( mesh_Q2, filename='mesh31_Q2.msh', ndim=2 )
  call read_mesh_gmsh ( mesh_Q1, filename='mesh31_Q1.msh', ndim=2 )

! curve for printing data on the centerline+sphere

  call add_to_mesh ( mesh_Q2, curve=[1,2,3,4,5,6,7,8], newnr=crv )
  call add_to_mesh ( mesh_Q1, curve=[1,2,3,4,5,6,7,8], newnr=crv )

  call mesh_convert ( mesh_Q2, blendmesh=mesh_Q1, mesh=mesh )

  call delete ( mesh_Q2, mesh_Q1 )

  call fill_mesh_parts ( mesh )

! plot points/curves for main and blend mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, filename='curves_main31.fig' )
  call plot_points_curves ( plot_options, mesh, blend=1, &
    filename='curves_blend31.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

  input_probdef%vec_elementdof(1)%a = 0           ! initialize to zero
  input_probdef%vec_elementdof(1)%a(1:9,1) = 2    ! velocity (Q2 mesh)
  input_probdef%vec_elementdof(1)%a(10:13,2) = 1  ! pressure (Q1 mesh)
  input_probdef%vec_elementdof(1)%a(1:9,3) = 1    ! scalar (Q2 mesh)
  input_probdef%vec_elementdof(1)%a(10:13,4) = 1  ! scalar (Q1 mesh)

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

  call write_vector_vtk ( mesh, problem, sysvector=sol, &
    dataname='velocity', filename='sphere31.vtk' )

  call write_scalar_vtk ( mesh, problem, sysvector=sol, &
    dataname='pressure', filename='pressure31.vtk', blend=1 )

! print solution only on center line

  call printtofile ( mesh, problem, filename='sol_on_center.out', &
    curve=crv, sysvector=sol )

! derivatives

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

  call write_scalar_vtk ( mesh, problem, vector=pressure, &
    dataname='pressure', filename='sphere31.vtk', append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=vorticity, &
    dataname='vorticity', filename='sphere31.vtk', append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=divergence, &
    dataname='divergence', filename='sphere31.vtk', append=.true. )

! print velocity vector only on center line

  call create_vector ( problem, velocity, physq=1 )
  call extract_physvector ( mesh, problem, sol, velocity )

  call printtofile ( mesh, problem, filename='velocity_on_center.out', &
    curve=crv, vector=velocity )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity, divergence )
  call delete ( coefficients )
  call delete ( oldvectors )

end program sphere31
