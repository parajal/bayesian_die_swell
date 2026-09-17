! Linear elastic problem on a 3D rectangular box. Stretching a bar.
! Standard displacement formulation.

program linear_elastic5

  use tfem_m
  use hsl_ma57_m
  use linear_elastic_elements_m
  use linear_elastic_post_m
  use figplot_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    nx = 4,            & ! number of elements in x-direction
    ny = 2,            & ! number of elements in y-direction
    nz = 4,            & ! number of elements in z-direction
!    nx = 10,            & ! number of elements in x-direction
!    ny = 5,             & ! number of elements in y-direction
!    nz = 10,            & ! number of elements in z-direction
    uintpl = 8,         & ! Q2 displacements
    physqdisp = 1,      & ! physical quantity nr of the displacements
    gauss = 3             ! 3x3x3 integration of hexahedra

  real(dp), parameter :: &
    lx = 10._dp,   & ! size in x-direction
    ly = 5._dp,    & ! size in y-direction
    lz = 10._dp,   & ! size in z-direction
    Emod = 2._dp, &  ! Young's modulus
    nu = 0.3_dp, &   ! Poisson's ratio
    ul = 1e-2_dp     ! displacement at x = L

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(plot_options_t) :: plot_options
  type(meshgen_options_t) :: meshgen_options
  type(coefficients_t) :: coefficients

  logical :: uniform = .false. ! uniform deformation, no edge effects


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(6) = physqdisp
  coefficients%i(10) = gauss
  coefficients%i(11) = gauss

  coefficients%r = 0
  coefficients%r(1) = Emod
  coefficients%r(3) = nu

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, &
    lz=lz, elshape=14, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! plot curves, surfaces and mesh

  plot_options%viewpoint=[1.,0.8,0.4]
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig', surfaces=[3,4,6] )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=1 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! displacement
  input_probdef%vec_elementdof(1)%a(:,2) = 1  ! scalar
  input_probdef%vec_elementdof(1)%a(:,3) = 6  ! symmetric tensor

  input_probdef%physq = [1]

  if ( uniform ) then
    call define_essential ( mesh, input_probdef, surfaces=[3,5], physq=1, &
      degfd=[1,0,0] )
    call define_essential ( mesh, input_probdef, point=1, physq=1, &
      degfd=[0,1,1] )
  else
    call define_essential ( mesh, input_probdef, surfaces=[3,5], physq=1 )
  end if

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  if ( uniform ) then
!   uniform deformation, no edge effects
    call fill_sysvector ( mesh, problem, sol, &
      surface1=3, physq=1, degfd=1, value=ul )
    call fill_sysvector ( mesh, problem, sol, &
      surface1=5, physq=1, degfd=1, value=0._dp )
    call fill_sysvector ( mesh, problem, sol, &
      point=1, physq=1, degsfd=[2,3], value=0._dp )
  else
!   edge effects
    call fill_sysvector ( mesh, problem, sol, &
      surface1=3, physq=1, vvalue=[ul,0._dp,0._dp] )
    call fill_sysvector ( mesh, problem, sol, &
      surface1=5, physq=1, value=0._dp )
  end if

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=linear_elastic_elem, coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )


! post processing

  call postprocessing_linear_elastic ( mesh, problem, 'le5.vtk', sol, &
    vec_scalar=2, vec_tensor=3, coefficients=coefficients )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program linear_elastic5
