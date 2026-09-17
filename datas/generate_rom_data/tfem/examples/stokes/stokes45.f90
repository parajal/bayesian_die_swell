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
! Internal subroutine as actual subroutine argument (f2008).

program stokes45

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use io_utils_m
  use figplot_m

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
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options
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
! normal_traction_free: is curve=3 traction free in normal direction?
  logical :: normal_traction_free = .true.

  integer :: i
  real(dp) :: Q(2,2), alpha

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

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

! rotate mesh over angle degrees

  alpha = 4 * atan(1._dp) * angle / 180._dp

  Q(1,:) = [ cos(alpha), -sin(alpha) ]
  Q(2,:) = [ sin(alpha),  cos(alpha) ]

  do i = 1, mesh%nnodes
    mesh%coor(i,:) = matmul ( Q, mesh%coor(i,:) )
  end do

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

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
      coefficients=coefficients, oldvectors=oldvectors, nodesub=nodesub )
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

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a fig file for plotting

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  call plot_color_fill ( plot_options, mesh, problem, 'pressure_color.fig', &
    vector=pressure )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )

  call plot_color_fill ( plot_options, mesh, problem, 'vorticity_color.fig', &
    vector=vorticity )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

! write binary file for reading by streamfunction

  open(unit=10,form='unformatted',file='velocity_bin.out')

  write(10) velocity%u

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
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

end program stokes45
