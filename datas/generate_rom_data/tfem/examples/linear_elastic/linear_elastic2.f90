! Linear elastic problem on a rectangle. Stretching a bar.
! Mixed displacement-pressure formulation. 2D, plane strain or plane stress.

program linear_elastic2

  use tfem_m
  use hsl_ma57_m
  use linear_elastic_elements_m
  use linear_elastic_post_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 displacements
    pintpl = 4,         & ! Q1 pressures
    physqdisp = 1,      & ! physical quantity nr of the displacement
    physqpress = 2,     & ! physical quantity nr of the pressure
    gauss = 3,          & ! 3x3 integration of quads
    nx=20,              & ! number of elements in x
    ny=20                 ! number of elements in y

  real(dp), parameter :: &
    lx = 10._dp,  &  ! length of domain in x-direction
    ly = 2._dp,  &   ! length of domain in y-direction
    Emod = 2._dp, &  ! Young's modulus
    nu = 0.3_dp, &   ! Poisson's ratio
    ul = 1e-2_dp     ! displacement at x = L

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients

  logical :: uniform = .false. ! uniform deformation, no edge effects


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(2) = pintpl
  coefficients%i(3) = 1   !  0: plane strain, 1: plane stress
  coefficients%i(4) = 0   !  0: compressible, 1: incompressible
  coefficients%i(6) = physqdisp
  coefficients%i(7) = physqpress
  coefficients%i(10) = gauss
  coefficients%i(11) = gauss

  coefficients%r = 0
  coefficients%r(1) = Emod
  coefficients%r(3) = nu

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%lx = lx
  meshgen_options%ly = ly

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! displacement
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1,    &  ! scalar
                 4,4,4,4,4,4,4,4,4 ], &  ! tensor with 4 components
                 [9,4] )

  input_probdef%physq = [1,2]

  if ( uniform ) then
    call define_essential ( mesh, input_probdef, curves=[2,4], physq=1, &
      degfd=[1,0] )
    call define_essential ( mesh, input_probdef, point=1, physq=1, degfd=[0,1] )
  else
    call define_essential ( mesh, input_probdef, curves=[2,4], physq=1 )
  end if

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  if ( uniform ) then
!   uniform deformation, no edge effects
    call fill_sysvector ( mesh, problem, sol, &
      curve1=2, physq=1, degfd=1, value=ul )
    call fill_sysvector ( mesh, problem, sol, &
      curve1=4, physq=1, degfd=1, value=0._dp )
    call fill_sysvector ( mesh, problem, sol, &
      point=1, physq=1, degfd=2, value=0._dp )
  else
!   edge effects
    call fill_sysvector ( mesh, problem, sol, &
      curve1=2, physq=1, vvalue=[ul,0._dp] )
    call fill_sysvector ( mesh, problem, sol, &
      curve1=4, physq=1, value=0._dp )
  end if

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=linear_elastic_elem, coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! post processing

  call postprocessing_linear_elastic ( mesh, problem, 'le2.vtk', sol, &
    vec_scalar=3, vec_tensor=4, coefficients=coefficients, assume33=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program linear_elastic2
