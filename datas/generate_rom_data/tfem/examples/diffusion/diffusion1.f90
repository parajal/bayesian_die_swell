! Steady scalar diffusion problem on a unit square with Dirichlet boundary
! conditions and natural boundary conditions.
! Varying alpha coefficient using a function or a vector in the nodes.

program diffusion1

  use tfem_m
  use hsl_ma57_m
  use convection_diffusion_elements_m
  use convection_diffusion_functions_m
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
    funcnr=8,           & ! function number for the right-hand side
    funcnr_alpha=1,     & ! function number for the alpha coefficient
    coeff_alpha=1,      & ! how to compute alpha coefficient:
                          ! 0: constant
                          ! 1: given by a function
                          ! 2: given by nodal values in a vector
    vfuncnr=1             ! function number for the flux vector

  real(dp), parameter :: &
    alpha = 1._dp     ! diffusion coefficient


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(vector_t), target :: alpha_nodes
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(2) = coeff_alpha
  coefficients%i(3) = funcnr_alpha
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]
  coefficients%i(16) = vfuncnr

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0

  coefficients%func => func
  coefficients%vfunc => vfunc

  coefficients%func1(1)%p => alphafunc   ! scalar function for alpha


! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a(:,1) = 1  ! for alpha in the nodes

  call define_essential ( mesh, input_probdef, curve1=1 )
  call define_essential ( mesh, input_probdef, curve1=3, curve2=4 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, func=func, funcnr=7 )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, curve2=4, func=func, funcnr=7 )

  if ( coeff_alpha == 2 ) then

!   fill alpha in a vector

    call create_vector ( problem, alpha_nodes, vec=1 )

!   as an example use a function to fill alpha in the nodes
    call fill_vector ( mesh, problem, alpha_nodes, &
      node1=1, node2=mesh%nnodes, func=func, funcnr=funcnr_alpha )

    call create ( oldvectors, nvec=1 )

    oldvectors%v(1)%p => alpha_nodes

  end if

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_diffusion_elem, coefficients=coefficients, &
    oldvectors=oldvectors )

  call add_boundary_elements ( mesh, problem, rhsd, curve=2, &
    elemsub=poisson_natboun, coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! print average sol to standard output

  print *, sum( abs(sol%u) ) / sol%n

! plot solution

  call plot_color_contour ( plot_options, mesh, problem, filename='sol.fig', &
    sysvector=sol )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  if ( coeff_alpha == 2 ) call delete ( alpha_nodes )
  call delete ( sysmatrix )
  call delete ( coefficients )
  if ( coeff_alpha == 2 ) call delete ( oldvectors )

end program diffusion1

