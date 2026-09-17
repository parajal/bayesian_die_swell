! Poisson problem on a unit square with Dirichlet boundary conditions
! This is the poisson1 problem with the following changes:
!  1) triangular elements
!  2) P1isoP2 shape functions
!  3) subdomain Gauss integration (using 4 subtriangles)

program poisson1b

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use problem_m
  use system_m
  use hsl_ma57_m
  use poisson_elements_m
  use poisson_functions_m
  use figplot_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 9,         & ! P1isoP2 scalar interpolation
    gauss = 3,          & ! 3 point Gauss integration on each of the four
                          ! subtriangles
    gaussb = 2,         & ! 2 point integration of boundary elements (on each
                          ! of the two subdomains)
    nx=20,              & ! number of elements in x
    ny=20,              & ! number of elements in y
    funcnr=4              ! function number for the right-hand side

  real(dp), parameter :: &
    alpha = 1._dp     ! diffusion coefficient


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd, solexact
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients
  type(solver_options_ma57_t) :: solver_options


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]

! divide the Gauss integration into four subtriangles (volume) and two line
! intervals (boundary)
  coefficients%i(32:33) = [ 2, 2 ]

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0

  coefficients%func => func

! create mesh

  meshgen_options%elshape = 4 ! use elshape=4 to create the mesh
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  mesh%element(:)%elshape = 33 ! change elshape to 33 for isoparametric mapping

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef )

!  input_probdef%elementdof(1)%a = (/1,1,1,1,1,1,1,1,1/)
  input_probdef%elementdof(1)%a = 1

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, func=func, funcnr=3 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=poisson_elem, &
    coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%integer_storage = 2

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

! create exact solution

  call create_sysvector ( problem, solexact )

  call fill_sysvector ( mesh, problem, solexact, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=3 )

! print maximum difference of sol-solexact to standard output

  print *, maxval( abs(sol%u -solexact%u) )

! plot solution

  call plot_color_contour ( plot_options, mesh, problem, filename='sol.fig', &
    sysvector=sol )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, solexact, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program poisson1b

