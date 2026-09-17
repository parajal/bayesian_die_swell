! Test problem for the surface advection add-on (generic routine)
! An analytical solution on a 2D mesh (3D velocity vector) is tested.
! First-order time integration.

program surface_advection1

  use tfem_m
  use surface_advection_elements_m
  use surface_advection_functions_m
  use hsl_ma41_m
  use io_utils_m
  use figplot_m

  implicit none


! constants

  integer, parameter :: &
    elshape = 4,        & ! 6-node triangles
    hintpl = 6,         & ! interpolation for the height, P2 triangle
    gauss = 4,          & ! 4th-order integration of triangles
!    elshape = 6,        & ! 9-node quads
!    hintpl = 8,         & ! interpolation for the height, Q2 quads
!    gauss = 3,          & ! 3x3 integration of quads
    nx=20,              & ! number of elements in x
    ny=20,              & ! number of elements in y
    numtimesteps=25,    & ! number of time steps
    funcnr=1,           & ! function number for the exact height function
    vfuncnr=1             ! function number for the velocity

  real(dp), parameter :: &
    deltat = 4e-2_dp, &  ! delta t
    beta = 1.0_dp        ! beta parameter for SUPG


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: soln
  type(sysvector_t) :: sol, rhsd, solexact
  type(vector_t), target :: velocity
  type(oldvectors_t) :: oldvectors
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: so

  integer :: step
  real(dp) :: time


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = gauss
  coefficients%i(2) = 2 ! velocity given by nodal values
  coefficients%i(4) = 1 ! 0: Galerkin, 1: SUPG
  coefficients%i(5) = 1 ! time-integration
  coefficients%i(6) = hintpl ! height interpolation
  coefficients%i(9) = 3 ! numerical table for Gauss

  coefficients%r = 0
  coefficients%r(4) = deltat
  coefficients%r(5) = beta

! create mesh

  meshgen_options%elshape = elshape
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

  call plot_mesh ( plot_options, mesh, filename='mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! 3D velocity vector

  call define_essential ( mesh, input_probdef, curve1=1 )
  call define_essential ( mesh, input_probdef, curve1=4 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, soln, rhsd, solexact )

! create velocity vector

  call create_vector ( problem, velocity, vec=1 )

! oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1, nvec=1 )

  oldvectors%s(1)%p => soln
  oldvectors%v(1)%p => velocity

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )


! fill initial solution vector

  t = 0

  call fill_sysvector ( mesh, problem, soln, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=funcnr )

! time stepping

  time = 0

  do step = 1, numtimesteps

    time = step * deltat

!   fill velocity vector with vector function

    t = time

    call fill_vector ( mesh, problem, velocity, &
      node1=1, node2=mesh%nnodes, vfunc=vfunc, vfuncnr=vfuncnr )

!   fill boundary condition at new time step

    call fill_sysvector ( mesh, problem, sol, &
      curve1=1, func=func, funcnr=funcnr )
    call fill_sysvector ( mesh, problem, sol, &
      curve1=4, func=func, funcnr=funcnr )

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=surface_advection_elem, coefficients=coefficients, &
      oldvectors=oldvectors )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    so%integer_storage=1.5

    call solve_system_ma41 ( sysmatrix, rhsd, sol, solver_options=so )

    call copy ( sol, soln )

!   fill solexact

    call fill_sysvector ( mesh, problem, solexact, &
      node1=1, node2=mesh%nnodes, func=func, funcnr=funcnr )

    print *, time,  maxval( abs(sol%u-solexact%u) )

  end do


! plot solution

  call plot_color_contour ( plot_options, mesh, problem, filename='sol.fig', &
    sysvector=sol )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, soln, rhsd, solexact )
  call delete ( sysmatrix )
  call delete ( velocity )
  call delete ( coefficients )
  call delete ( oldvectors )

end program surface_advection1

