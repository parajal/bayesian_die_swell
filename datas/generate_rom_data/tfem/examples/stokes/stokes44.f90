! Axisymmetrical Stokes problem of a sphere at the centerline of a cylindrical
! container. A force is applied on the particle.
! Gmsh mesh.
! Interfacial slip on the particle boundary.

program stokes44

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m
  use subs44_m
  use limits_m

  implicit none

  integer, parameter :: &
    uintpl = 6,         & ! P2 velocities
    pintpl = 2,         & ! P1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 6,          & ! 6-point Gauss integration of triangles
    gaussb = 3,         & ! 3 point integration of boundary elements
    coorsys = 1           ! axis symmetric problem

! definitions

  type(mesh_t) :: mesh, mesh_ext, mesh_int
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol
  type(sysvector_t) :: rhsd
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients

! variables

  real(dp), parameter ::   &
    eta = 1._dp,           & ! fluid viscosity
    eta_slip = 0.1_dp        ! slip parameter

  real(dp) :: up(1)

! skip test for multiple connected elements, otherwise a mesh is plotted
! through data plots made with figplot
  TEST_FOR_MULTIPLE_SIDELEM = .false.

! fill coefficients
  call create_coefficients ( coefficients, ncoefi=600, ncoefr=501 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gaussb ]
  coefficients%i(12:) = 0
  coefficients%i(23) = coorsys

  coefficients%r(1) = eta
  coefficients%r(2:) = 0
  coefficients%r(501) = eta_slip

  force = [ 1.0_dp ]

! read mesh from gmsh output file
  call read_mesh_gmsh ( mesh_ext, filename='mesh44.msh', ndim=2, &
    physgeom=.true. )

! create particle mesh from curve=2
  call mesh_convert ( mesh_ext, mesh_int, from_curve=2 )

! merge the fluid and particle meshes
  call mesh_merge ( mesh_ext, mesh_int, mesh, nogroupmerge=.true. )

! add a new curve given by the particle mesh
  call add_to_mesh ( mesh, curvefromgroups=[2] )

  call delete ( mesh_ext, mesh_int )

  call fill_mesh_parts ( mesh )

! plot mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  plot_options%fontsize = 10
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

! group 1
  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,    &  ! pressure
                 1,1,1,1,1,1,    &  ! scalar
                 3,3,3,3,3,3 ], &  ! tensor
                 [6,4] )

! group 2 (only velocities)
  input_probdef%vec_elementdof(2)%a =   &
      reshape ( [2,2,2,    &  ! velocity
                 0,0,0,    &  ! pressure
                 0,0,0,    &  ! scalar
                 0,0,0 ], &  ! tensor
                 [3,4] )

  input_probdef%physq = [1,2]
  input_probdef%probnr = 1

! define essential boundaries

! cell boundaries
  call define_essential ( mesh, input_probdef, curve1=1, degfd=[0,1], &
    physq=physqvel )
  call define_essential ( mesh, input_probdef, curve1=3, degfd=[0,1], &
    physq=physqvel )
  call define_essential ( mesh, input_probdef, curve1=4, curve2=6, &
    physq=physqvel )

! particle ring
  call define_essential ( mesh, input_probdef, curve1=7, degfd=[0,1], &
    physq=physqvel )

! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! define constraint on the disk: rigid-body motion
  call define_constraint ( mesh, input_probdef, curve1=7, physq=physqvel, &
    discretization='collocation', nodedof=1, naddunknowns=1 )

! normal velocity coupling
  call define_constraint ( mesh, input_probdef, physq=physqvel, &
    curve1=2, curve2=7, discretization='weak', elementdof=[1,1,1] )

! define connection: slip
  call define_connection ( mesh, input_probdef, curve1=2, curve2=7, &
   discretization='weak', physq=physqvel )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)
  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill boundary conditions

! cell boundaries
  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, physq=physqvel, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=physqvel, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=4, curve2=6, physq=physqvel, value=0._dp )

! particle ring
  call fill_sysvector ( mesh, problem, sol, &
    curve1=7, physq=physqvel, degfd=2, value=0._dp )

! pressure level
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqpress, value=0._dp )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_connection ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients, elgroup1=1 )

! fill the matrix and vector for group 2 by zeros

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    coefficients=coefficients, elgroup1=2, zeromatvec=.true., addmatvec=.true. )

! rigid-body motion

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=elementc, coefficients=coefficients, &
    addmatvec=.true. )

! coupling of normal velocities

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_normal_vel, &
    coefficients=coefficients, addmatvec=.true. )

! interfacial slip

  call build_system_connection ( mesh, problem, sysmatrix, rhsd, &
    elemsub=elementcn, coefficients=coefficients, addmatvec=.true. )

  call check_filled_sysmatrix ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

! solve gradient/velocity/pressure problem

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call delete ( sysmatrix )

! get velocity
  call get_sysvector_constraint ( mesh, problem, sol, constraint=1, &
    addunknowns=.true., u=up )

  write ( *, '(1X,A,3F10.6)' ) 'velocity = ', up(1)

! plot velocity field

  call plot_color_contour ( plot_options, mesh, problem, 'u.fig', &
    physq=1, degfd=1, sysvector=sol, groups=[1] )
  call plot_color_contour ( plot_options, mesh, problem, 'v.fig', &
    physq=1, degfd=2, sysvector=sol, groups=[1] )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( coefficients )

end program stokes44
