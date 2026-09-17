! Linear elastic problem on a square with a hole. Stretching in x-direction
! Standard displacement formulation. 2D, plane strain or plane stress.

program linear_elastic20

  use tfem_m
  use hsl_ma57_m
  use linear_elastic_elements_m
  use linear_elastic_post_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,         & ! Q2 displacements
    physqdisp = 1,      & ! physical quantity nr of the displacements
    gauss = 7             ! 3x3 integration of quads

  real(dp), parameter :: &
    Emod = 2._dp, &  ! Young's modulus
    nu = 0.3_dp, &   ! Poisson's ratio
    ul = 1e-2_dp     ! displacement at x = L

! definitions

  type(mesh_t) :: mesh, mesh1
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients

  logical :: uniform = .true. ! uniform deformation, no edge effects


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(3) = 0   !  0: plane strain, 1: plane stress
  coefficients%i(6) = physqdisp
  coefficients%i(10) = gauss
  coefficients%i(11) = gauss

  coefficients%r = 0
  coefficients%r(1) = Emod
  coefficients%r(3) = nu

! read mesh

  call read_mesh_gmsh ( mesh1, 'hole_in_plate.msh', ndim=2 )

  call mesh_convert ( mesh1, mesh, remove_isolated_nodes=.true. )
  call delete(mesh1)

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=1 )

  input_probdef%vec_elementdof(1)%a(:,1) = 2  ! displacement
  input_probdef%vec_elementdof(1)%a(:,2) = 1  ! scalar
  input_probdef%vec_elementdof(1)%a(:,3) = 4  ! symmetric tensor

  input_probdef%physq = [1]

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

  call postprocessing_linear_elastic ( mesh, problem, 'le20.vtk', sol, &
    vec_scalar=2, vec_tensor=3, coefficients=coefficients, assume33=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program linear_elastic20
