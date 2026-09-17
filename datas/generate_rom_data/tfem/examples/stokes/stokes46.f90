! Stokes problem on a 3D square cavity. Only half of the domain is solved
! due to symmetry conditions.
! Arbitrary rotation of mesh over the y-axis
! Boundary conditions applied using transformations.
! Optional Neumann condition in the transformed system.
! Show various possibilities for applying transformations:
! 1) A global vector/matrix valid for each node on the geometry.
! 2) A user subroutine called for each node on the geometry.
! 3) In the direction of the normal on the geometry.
! Internal subroutine as actual subroutine argument (f2008).

program stokes46

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    nx = 10,            & ! number of elements in x-direction
    ny = 5,             & ! number of elements in y-direction
    nz = 10,            & ! number of elements in z-direction
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3             ! 3x3x3 integration of hexahedra

  real(dp), parameter :: &
    lx = 1._dp,          & ! size in x-direction
    ly = 0.5_dp,         & ! size in y-direction
    lz = 1._dp,          & ! size in z-direction
    eta = 1._dp,         & ! viscosity
    angle = 30._dp         ! angle of rotation the mesh in degrees

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, dudz
  type(plot_options_t) :: plot_options
  type(meshgen_options_t) :: meshgen_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

! OPTIONS for demonstrating the transformation possibilities:

! user_subroutine: use user subroutine
  logical :: user_subroutine = .true.
! normal_vector: use direction of the normal vector?
  logical :: normal_vector = .false.
! first_direction_only: specify first direction only of the constant global
! transformation, otherwise the full constant transformation is supplied,
! both using Amatrix.
  logical :: first_direction_only = .false.
! normal_traction_free: is surface=6 traction free in normal direction?
  logical :: normal_traction_free = .true.

! variables

  integer :: presnod(8) = [1,3,9,7,19,21,27,25], i
  real(dp) :: Q(3,3), alpha


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, &
    lz=lz, elshape=14, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

! rotate mesh over angle degrees

  alpha = 4 * atan(1._dp) * angle / 180._dp

  Q(1,:) = [ cos(alpha), 0._dp, -sin(alpha) ]
  Q(2,:) = [ 0._dp     , 1._dp,  0._dp      ]
  Q(3,:) = [ sin(alpha), 0._dp,  cos(alpha) ]

  do i = 1, mesh%nnodes
    mesh%coor(i,:) = matmul ( Q, mesh%coor(i,:) )
  end do

  call fill_mesh_parts ( mesh )

! plot curves, surfaces and mesh

  plot_options%viewpoint=[1.,0.8,0.4]
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig', surfaces=[3,4,6] )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(presnod,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2]

! Define essential boundary conditions

  call define_essential ( mesh, input_probdef, surface1=1, surface2=3, physq=1 )
  call define_essential ( mesh, input_probdef, surface1=4, physq=1, &
    degfd=[0,1,0] )
  if ( normal_traction_free ) then
!   traction free in normal direction of surface 6
    call define_essential ( mesh, input_probdef, surface1=5, physq=1 )
    if ( normal_vector ) then
     call define_essential ( mesh, input_probdef, surface1=6, physq=1, &
       degfd=[0,1,1] )
    else
     call define_essential ( mesh, input_probdef, surface1=6, physq=1, &
       degfd=[1,1,0] )
    end if
  else
!   full Dirichlet on surface 6
    call define_essential ( mesh, input_probdef, surface1=5, surface2=6, &
      physq=1 )
  end if
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

! Define transformation for boundary conditions on surface 6

  if ( user_subroutine ) then
!   use user subroutine
    call define_transformation ( mesh, input_probdef, surface=6, physq=1 )
  else if ( normal_vector ) then
!   use (negative) normal vector as first direction
    call define_transformation ( mesh, input_probdef, surface=6, physq=1, &
      normalvector=-1, v2=[0._dp,1._dp,0._dp] )
  else if ( first_direction_only ) then
!   first direction = the tangential direction after transformation
    call define_transformation ( mesh, input_probdef, surface=6, physq=1, &
      Amatrix=Q(:,1:2) )
  else
!   full matrix
    call define_transformation ( mesh, input_probdef, surface=6, physq=1, &
      Amatrix=Q )
  end if

  call problem_definition ( input_probdef, mesh, problem )

! build tranformation matrix (stored in problem structure)

  if ( user_subroutine ) then
    call build_transformation_matrix ( mesh, problem, &
      coefficients=coefficients, oldvectors=oldvectors, nodesub=nodesub )
  else
    call build_transformation_matrix ( mesh, problem )
  end if

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    surface1=1, surface2=3, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=4, physq=1, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=5, surface2=6, physq=1, value=0._dp )
  if ( normal_vector ) then
    call fill_sysvector ( mesh, problem, sol, &
      surface1=6, physq=1, degfd=2, value=1._dp )
  else
    call fill_sysvector ( mesh, problem, sol, &
      surface1=6, physq=1, degfd=1, value=1._dp )
  end if
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! transform the solution vector from the local to the global system

  call transform_to_global ( problem, sol )

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, dudz, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=3
  call derive_vector ( mesh, problem, dudz, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )


! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='cavity.vtk' )

  call write_vector_vtk ( mesh, problem, filename='cavity.vtk', &
    dataname='velocity_vector', sysvector=sol, append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=dudz, filename='cavity.vtk', &
    dataname='dudz', append=.true. )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( velocity, pressure, dudz )
  call delete ( coefficients )
  call delete ( oldvectors )

contains

! user subroutine for local transformation

  subroutine nodesub ( mesh, problem, geom, node, first, last, &
    coefficients, oldvectors, nodemat )
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: geom, node
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: nodemat

    nodemat = Q

  end subroutine nodesub

end program stokes46
