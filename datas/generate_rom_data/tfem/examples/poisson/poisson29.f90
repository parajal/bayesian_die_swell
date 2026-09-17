! Poisson problem on a unit cube with Dirichlet boundary conditions and
! natural boundary conditions
! Hybrid mesh with linear tetrahedron and pyramid elements

program poisson29

  use tfem_m
  use hsl_ma57_m
  use poisson_elements_m
  use poisson_functions_m
  use io_utils_m
  use limits_m

  implicit none


! constants

  integer, parameter :: &
    uintpl1 = 2,        & ! scalar interpolation in group 1 (tetrahedra)
    uintpl2 = 17,       & ! scalar interpolation in group 2 (pyramids)
    gauss = 4,          & ! order of integration
    gaussb = 4,         & ! order integration of boundary elements
    funcnr=17,          & ! function number for the right-hand side
    vfuncnr=5             ! function number for the flux vector

  real(dp), parameter :: &
    alpha = 1._dp     ! diffusion coefficient

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd, solexact
  type(vector_t) :: grad, gradexact
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients(2)

  WARN_ON_TEST_FOR_VARYING_SIDES = .false.

! fill coefficients

  call create_mcoefficients ( coefficients, ncoefi=100, ncoefr=50 )

! group 1
  coefficients(1)%i = 0
  coefficients(1)%i(1) = uintpl1
  coefficients(1)%i(6) = 2  ! use a vector function for h_N
  coefficients(1)%i(10:12) = [ gauss, gaussb, funcnr ]
  coefficients(1)%i(16) = vfuncnr
  coefficients(1)%i(40) = 3  ! numeric Gauss values

  coefficients(1)%r(1) = alpha
  coefficients(1)%r(2:) = 0

  coefficients(1)%func => func
  coefficients(1)%vfunc => vfunc

! group 2
  coefficients(2) = coefficients(1)
  coefficients(2)%i(1) = uintpl2


! read mesh

  call read_mesh_gmsh ( mesh, filename='cube3.msh' )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a(:,1) = 3
  input_probdef%elementdof(2)%a = 1
  input_probdef%vec_elementdof(2)%a(:,1) = 3

  call define_essential ( mesh, input_probdef, surfaces=[1,2,4,5,6] )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    surfaces=[1,2,4,5,6], func=func, funcnr=16 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=poisson_elem, &
    mcoefficients=coefficients )

  call add_boundary_elements ( mesh, problem, rhsd, surface=3, &
    elemsub=poisson_natboun, coefficients=coefficients(1) )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! create exact solution

  call create_sysvector ( problem, solexact )

  call fill_sysvector ( mesh, problem, solexact, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=16 )

! print maximum difference of sol-solexact to standard output

  print *, maxval( abs(sol%u - solexact%u) )

! create grad

  call create ( problem, grad, vec=1 )

! create oldvectors

  call create ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, grad, elemsub=poisson_deriv, &
    mcoefficients=coefficients, oldvectors=oldvectors )

! create exact solution of grad

  call create_vector ( problem, gradexact, vec=1 )

  call fill_vector ( mesh, problem, gradexact, &
    node1=1, node2=mesh%nnodes, vfunc=vfunc, vfuncnr=6 )

! print maximum difference of grad-gradexact to standard output

  print *, maxval( abs(grad%u - gradexact%u) )

! plot solution and derivatives

  call write_scalar_vtk ( mesh, problem, filename='sol29.vtk', &
    dataname='sol', sysvector=sol )

  call write_vector_vtk ( mesh, problem, filename='sol29.vtk', &
    dataname='grad', vector=grad, append=.true. )

  sol%u = sol%u - solexact%u

  call write_scalar_vtk ( mesh, problem, filename='sol29.vtk', &
    dataname='sol_error', sysvector=sol, append=.true. )

  grad%u = grad%u - gradexact%u

  call write_vector_vtk ( mesh, problem, filename='sol29.vtk', &
    dataname='grad_error', vector=grad, append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, solexact, rhsd )
  call delete ( grad, gradexact )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program poisson29
