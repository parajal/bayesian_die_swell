! Poisson problem on a unit square with Dirichlet boundary conditions
! Gmsh mesh with 6 nodes triangles.
! L2-projection of solution on the same mesh, but with a different shape
! function.

program projection5

  use tfem_m
  use hsl_ma57_m
  use convection_diffusion_elements_m
  use projection_elements_m
  use functions_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  integer, parameter ::  &
    uintpl = 6,          & ! scalar interpolation
    gauss = 6,           & ! 6 point integration of triangles
    gaussb = 3,          & ! 3 point integration of boundary elements
    uintplp = 2,         & ! scalar interpolation projection
    gaussp = 6,          & ! 6 point integration of triangles for projection
    ncompv = 1,          & ! number of components input vector for projection
    funcnr=4               ! function number for the right-hand side

  real(dp), parameter :: &
    alpha = 1._dp     ! diffusion coefficient

! definitions

  type(mesh_t), target :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd, solexact
  type(vector_t), target :: vec_p
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients

  type(input_probdef_t) :: input_probdef2
  type(problem_t) :: problem2
  type(sysmatrix_t) :: sysmatrix2
  type(sysvector_t), target :: sol2
  type(sysvector_t) :: rhsd2
  type(vector_t) :: sol2v, solexact2v
  type(coefficients_t) :: coefficients2
  type(oldvectors_t) :: oldvectors2


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0

  coefficients%func => func

! read gmsh file of a square geometry discretized with
! 6-node triangular elements

  call read_mesh_gmsh ( mesh, 'square.msh', ndim=2 )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, filename='curves.fig' )
  call plot_mesh ( plot_options, mesh, filename='mesh.fig' )

  call printinfo ( mesh, printlevel=4 )

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a = ncompv ! for transfer to projection module

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

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! create exact solution

  call create_sysvector ( problem, solexact )

  call fill_sysvector ( mesh, problem, solexact, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=3 )

! print maximum difference of sol-solexact to standard output

  print *, maxval( abs(sol%u -solexact%u) )

! plot solution

  call plot_color_contour ( plot_options, mesh, problem, filename='sol.fig', &
    sysvector=sol )

! create vector and copy data for input of projection problem

  call create ( problem, vec_p, vec=1 )

  call transfer_data ( mesh, problem, sysvector1=sol, vector2=vec_p )


! projection problem
! ------------------

! problem definition for projection

  call create_input_probdef ( mesh, input_probdef2, nvec=1 )

  input_probdef2%elementdof(1)%a = [1,0,1,0,1,0]
  input_probdef2%vec_elementdof(1)%a(:,1) = 1  ! plotting

  call problem_definition ( input_probdef2, mesh, problem2 )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem2, sol2 )
  call create_sysvector ( problem2, rhsd2 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix2, mesh, problem2, &
    symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix2 )

! fill coefficients

  call create_coefficients ( coefficients2, ncoefi=100, ncoefr=50 )

  coefficients2%i = 0
  coefficients2%i(1) = uintpl
  coefficients2%i(2) = ncompv
  coefficients2%i(3) = uintpl
  coefficients2%i(4) = uintplp
  coefficients2%i(5) = 1 ! same mesh
  coefficients2%i(10) = gaussp

  coefficients2%r = 0

! oldvectors

  call create_oldvectors ( oldvectors2, nsysvec=1, nvec=1, nprob=1 )

  oldvectors2%v(1)%p => vec_p
  oldvectors2%p(1)%p => problem

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem2, sysmatrix2, rhsd2, &
    elemsub=projection_elem, coefficients=coefficients2, &
    oldvectors=oldvectors2 )

  call add_effect_of_essential_to_rhs ( problem2, sysmatrix2, sol2, rhsd2 )

  call solve_system_ma57 ( sysmatrix2, rhsd2, sol2 )

! plotting

  call create ( problem2, sol2v, vec=1 )

  oldvectors2%s(1)%p => sol2

  call derive_vector ( mesh, problem2, sol2v, elemsub=deriv_projection, &
    coefficients=coefficients2, oldvectors=oldvectors2 )

  call create ( problem2, solexact2v, vec=1 )

  call fill_vector ( mesh, problem2, solexact2v, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=3 )

! print maximum difference of projected sol-solexact to standard output

  print *, maxval( abs(sol2v%u -solexact2v%u) )

! plot projected solution

  call plot_color_contour ( plot_options, mesh, problem2, filename='sol2.fig',&
    vector=sol2v )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, solexact, rhsd )
  call delete ( sysmatrix )
  call delete ( vec_p )
  call delete ( coefficients )

  call delete ( problem2 )
  call delete ( input_probdef2 )
  call delete ( sol2, rhsd2 )
  call delete ( sol2v, solexact2v )
  call delete ( sysmatrix2 )
  call delete ( coefficients2 )

end program projection5
