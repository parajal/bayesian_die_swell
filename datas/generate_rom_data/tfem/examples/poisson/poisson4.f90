! Poisson problem on a curved domain with Dirichlet boundary conditions and
! natural boundary conditions. Derived quantity.

program poisson4

  use tfem_m
  use hsl_ma57_m
  use poisson_elements_m
  use poisson_functions_m
  use functions_m
  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3 integration of quads
    gaussb = 3,         & ! 3 point integration of boundary elements
    nx=20,              & ! number of elements in x
    ny=20,              & ! number of elements in y
    funcnr=4,           & ! function number for the right-hand side
    vfuncnr=2             ! function number for the flux vector

  real(dp), parameter :: &
    alpha = 1._dp     ! diffusion coefficient

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: grad, gradexact
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(6) = 2  ! use a vector function for h_N
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]
  coefficients%i(16) = vfuncnr

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0

  coefficients%func => func
  coefficients%vfunc => vfunc

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  meshgen_options%regionshape = 3
  meshgen_options%x2d = &
    reshape ( [ 0.0_dp, 2.0_dp, 2.0_dp, -1.0_dp,    &
                0.0_dp, 0.0_dp, 2.0_dp,  3.0_dp ], &
              [4,2] )
  meshgen_options%curved(1:4) = [ .false., .true., .true., .true. ]
  meshgen_options%funcnr(1:4) = [ 1, 2, 3, 4 ]

  call quadrilateral2d ( mesh, meshgen_options, func=cfunc )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = [1,1,1,1,1,1,1,1,1]

  input_probdef%vec_elementdof(1)%a = reshape ( [2,2,2,2,2,2,2,2,2], [9,1])

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

! create grad

  call create ( problem, grad, vec=1 )

! create oldvectors

  call create ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, grad, elemsub=poisson_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! create exact solution of grad

  call create_vector ( problem, gradexact, vec=1 )

  call fill_vector ( mesh, problem, gradexact, degfd=1, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=5 )
  call fill_vector ( mesh, problem, gradexact, degfd=2, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=6 )

! print maximum difference of grad-gradexact to standard output

  print *, maxval( abs(grad%u -gradexact%u) )

! plot solution

  plot_options%fontsize = 14

  call plot_color_contour ( plot_options, mesh, problem, filename='dudx.fig', &
    degfd=1, vector=grad )
  call plot_color_contour ( plot_options, mesh, problem, filename='dudy.fig', &
    degfd=2, vector=grad )

  plot_options%scalevector = 0.05

  call plot_vector ( plot_options, mesh, problem, filename='gradu.fig', &
    vector=grad )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( grad, gradexact )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program poisson4

