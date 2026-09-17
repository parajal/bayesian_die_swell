! Poisson problem on a curved domain with Dirichlet boundary conditions and
! natural boundary conditions.
! Spectral elements.
! This problem is similar to poisson3, but now with spectral elements.

program poisson21

  use tfem_m
  use hsl_ma57_m
  use poisson_elements_m
  use poisson_functions_m
  use functions_m
  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 13,  & ! scalar interpolation
    p = 8,        & ! polynomial order
    gauss = p+1,  & ! number of Gauss points
    gaussb = p+1, & ! number of Gauss points on the boundary
    inttype = 1,  & ! Gauss-Legendre-Lobatto for integration
    nx=5,         & ! number of elements in x
    ny=5,         & ! number of elements in y
    funcnr=4,     & ! function number for the right-hand side
    vfuncnr=2       ! function number for the flux vector

  real(dp), parameter :: &
    alpha = 1._dp     ! diffusion coefficient

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh_plot
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd, solexact
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]
  coefficients%i(6) = 2  ! use a vector function for h_N
  coefficients%i(16) = vfuncnr
  coefficients%i(39) = p
  coefficients%i(40) = inttype

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0

  coefficients%func => func
  coefficients%vfunc => vfunc

! create mesh

  meshgen_options%elshape = 102
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%p = p

  meshgen_options%regionshape = 3
  meshgen_options%x2d = &
    reshape ( [ 0.0_dp, 2.0_dp, 2.0_dp, -1.0_dp,    &
                0.0_dp, 0.0_dp, 2.0_dp,  3.0_dp ], &
              [4,2] )
  meshgen_options%curved(1:4) = .true.
  meshgen_options%funcnr(1:4) = [ 1, 2, 3, 4 ]

  call quadrilateral2d ( mesh, meshgen_options, func=cfunc )

  call fill_mesh_parts ( mesh )

! Create plot mesh

  call mesh_convert ( mesh, mesh_plot, elementshapes='spectraltolinear', &
    warn=.false. )

  call fill_mesh_parts( mesh_plot )

  call plot_points_curves ( plot_options, mesh_plot, filename='curves.fig' )
  call plot_mesh ( plot_options, mesh_plot, filename='mesh_plot.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef )

  input_probdef%elementdof(1)%a = 1

  call define_essential ( mesh, input_probdef, curves=[1,3,4] )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curves=[1,3,4], func=func, funcnr=3 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=poisson_elem, &
    coefficients=coefficients )

  call add_boundary_elements ( mesh, problem, rhsd, curve=2, &
    elemsub=poisson_natboun, coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! create exact solution

  call create_sysvector ( problem, solexact )

  call fill_sysvector ( mesh, problem, solexact, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=3 )

! print maximum difference of sol-solexact to standard output

  print *, maxval( abs(sol%u -solexact%u) )

! plot solution

  plot_options%fontsize = 14

  call plot_color_contour ( plot_options, mesh_plot, problem, &
    filename='sol.fig', sysvector=sol )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh, mesh_plot )
  call delete ( sol, solexact, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program poisson21

