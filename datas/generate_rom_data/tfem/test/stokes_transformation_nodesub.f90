module mod45_m

  use tfem_m
  implicit none

  real(dp) :: Q(2,2)

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

end module mod45_m

! Stokes problem on a unit square.
! Lid-driven cavity flow with optional Dirichlet or traction free boundary
! condition in normal direction of the lid.
! Arbitrary rotation of mesh.
! Boundary conditions applied using transformations.
! Optional Neumann condition in the transformed system.
! Show various possibilities for applying transformations:
! 1) A global vector/matrix valid for each node on the geometry.
! 2) A user subroutine called for each node on the geometry.
! 3) In the direction of the normal on the geometry.

program stokes45

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use mod45_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=20,              & ! number of elements in x
    ny=20                 ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp,  & ! viscosity
    angle = 30._dp  ! angle of rotation of the mesh in degrees

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(coefficients_t) :: coefficients
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd

! OPTIONS for demonstrating the transformation possibilities:

! user_subroutine: use user subroutine
  logical :: user_subroutine = .true.
! normal_vector: use direction of the normal vector?
  logical :: normal_vector = .false.
! first_direction_only: specify first direction only of the constant global
! transformation, otherwise the full constant transformation is supplied,
! both using Amatrix.
  logical :: first_direction_only = .false.
! normal_traction_free: is curve=3 traction free in normal direction?
  logical :: normal_traction_free = .true.

  integer :: i
  real(dp) :: alpha

  if ( user_subroutine ) then
    normal_vector = .false.
  end if

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! rotate mesh over angle degrees

  alpha = 4 * atan(1._dp) * angle / 180._dp

  Q(1,:) = [ cos(alpha), -sin(alpha) ]
  Q(2,:) = [ sin(alpha),  cos(alpha) ]

  do i = 1, mesh%nnodes
    mesh%coor(i,:) = matmul ( Q, mesh%coor(i,:) )
  end do

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,1,0,0,    &  ! pressure
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                   [9,3] )

  input_probdef%physq = [1,2]

! Define essential boundary conditions

  if ( normal_traction_free ) then
!   traction free in normal direction of curve 3
    call define_essential ( mesh, input_probdef, curve1=1, curve2=2, physq=1 )
    if ( normal_vector ) then
     call define_essential ( mesh, input_probdef, curve1=3, physq=1, &
       degfd=[0,1] )
    else
     call define_essential ( mesh, input_probdef, curve1=3, physq=1, &
       degfd=[1,0] )
    end if
    call define_essential ( mesh, input_probdef, curve1=4, physq=1 )
  else
!   full Dirichlet on curve 3
    call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )
  end if
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

! Define transformation for boundary conditions on curve 3

  if ( user_subroutine ) then
!   use user subroutine
    call define_transformation ( mesh, input_probdef, curve=3, physq=1 )
  else if ( normal_vector ) then
!   use (negative) normal vector as first direction
    call define_transformation ( mesh, input_probdef, curve=3, physq=1, &
      normalvector=-1 )
  else if ( first_direction_only ) then
!   first direction = the tangential direction after transformation (on curve=3)
    call define_transformation ( mesh, input_probdef, curve=3, physq=1, &
      Amatrix=Q(:,1:1) )
  else
!   full matrix
    call define_transformation ( mesh, input_probdef, curve=3, physq=1, &
      Amatrix=Q )
  end if

  call problem_definition ( input_probdef, mesh, problem )

! build tranformation matrix (stored in problem structure)

  if ( user_subroutine ) then
    call build_transformation_matrix ( mesh, problem, &
      coefficients=coefficients, nodesub=nodesub )
  else
    call build_transformation_matrix ( mesh, problem )
  end if

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions defined in the
! transformed system

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, physq=1, value=0._dp )
  if ( normal_vector ) then
    call fill_sysvector ( mesh, problem, sol, &
      curve1=3, physq=1, degfd=2, value=1._dp )
  else
    call fill_sysvector ( mesh, problem, sol, &
      curve1=3, physq=1, degfd=1, value=1._dp )
  end if
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! transform the solution vector from the local to the global system

  call transform_to_global ( problem, sol )

! print average sol to standard output

  print *, sum( abs(sol%u) ) / sol%n

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )

end program stokes45
