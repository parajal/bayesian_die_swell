#if HSL_EXTRA

! Poisson problem on a unit square with Dirichlet boundary conditions
! Iterative solver BiCGSTAB with algebraic multigrid preconditioner.

program poisson14

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use problem_m
  use system_m
  use hsl_solve_mi20_m
  use poisson_elements_m
  use poisson_functions_m
  use figplot_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3 integration of quads
    gaussb = 3,         & ! 3 point integration of boundary elements
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
  type(solver_options_mi20_t) :: solver_options


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0

  coefficients%func => func

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef )

  input_probdef%elementdof(1)%a = [1,1,1,1,1,1,1,1,1]

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, func=func, funcnr=3 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=poisson_elem, &
    coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call set_solver_options ( solver_options, printlevel=2, move_matrix=.true., &
    eps_rel=1e-7_dp, itsolver=5 )

  call solve_system_mi20 ( sysmatrix, rhsd, sol, solver_options=solver_options )

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

end program poisson14

#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on hsl_extra', &
    ' - add the libhsl3 library for linking', &
    ' - set preprocessing macro HSL_EXTRA in Mdefs.mk'
end
#endif

