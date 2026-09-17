! Poisson problem on a unit square with Dirichlet boundary conditions
! Line source, using a system build on an object.

module subs9_m

  use kind_defs_m
  use tfem_m

  implicit none

! global parameters

  integer :: ndim = 2, ndf = 9, nodalpb = 3

contains


! Element routine for a source on a line object:
!   f = (a + b*u)delta(n) where delta(n) is a delta function in the
! direction of the normal.

  subroutine poisson_source ( mesh, problem, ei, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(eleminfo_t), intent(in) :: ei
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, object, elemo

    real(dp) :: wg(ei%nintps), curvel(ei%nintps)
    real(dp) :: normal(ei%nintps,ndim), xig(ei%nintps)
    real(dp) :: xr(ei%nintps,ndim), phi(ei%nintps,ndf)
    real(dp) :: dxdxi(ei%nintps,ndim)
    real(dp) :: theta(ei%nintps,nodalpb), dtheta(ei%nintps,nodalpb)
    real(dp) :: x(nodalpb,ndim)

    object = ei%object
    elemo  = ei%elemo

!   reference coordinates and shape function in mesh

    xr = mesh%objects(object)%refcoor_int(ei%intps,:,elemo)

    call shape_quad_Q2 ( xr, phi )

    wg = mesh%objects(object)%wg(ei%intps)

!   shape function in object

    xig = mesh%objects(object)%xig(ei%intps,1)

    call shape_line_P2 ( xig, theta, dtheta )

    call get_coordinates_object ( mesh, elemo, x, object )

    call isoparametric_deformation_curve ( x, dtheta, dxdxi, curvel, normal )

    if ( vector ) then

      do i = 1, ndf
        elemvec(i) = coefficients%r(49) * sum ( phi(:,i) * curvel * wg )
      end do

    end if

    if ( matrix ) then

      do i = 1, ndf
        do j = 1, ndf
          elemmat(i,j) = &
            coefficients%r(50) * sum ( phi(:,i) * phi(:,j) * curvel * wg )
        end do
      end do

    end if

  end subroutine poisson_source

end module subs9_m

program poisson9

  use math_defs_m
  use tfem_m
  use hsl_ma57_m
  use poisson_elements_m
  use figplot_m
  use subs9_m

  implicit none


! constants

! rotate=.true.: rotate the object
  logical, parameter :: rotate = .true.

  integer, parameter :: &
    uintpl = 8,   & ! scalar interpolation
    gauss = 3,    & ! 3x3 integration of quads
    gaussb = 3,   & ! 3 point integration of boundary elements
    nx = 40,      & ! number of elements in x
    ny = 40,      & ! number of elements in y
    ne = 20,      & ! number of elements on source line
    nsubint = 20    ! number of subintegration intervals for rotate=.true.

  real(dp), parameter :: &
    alpha = 1._dp,   &  ! diffusion coefficient
    length = 0.5_dp, &  ! length of line source
    shiftx = 0.0_dp, &  ! shift in y of the center of the line (rotate=.true.)
    shifty = 0.0_dp, &  ! shift in x of the center of the line (rotate=.true.)
    angle = pi/8,    &  ! rotation angle (rotate=.true.)
    a = 1._dp,       &  ! coefficient in source: f = (a + b*u)*delta(n)
    b = -1._dp          ! coefficient in source: f = (a + b*u)*delta(n)


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh_source
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients

  integer :: elem, i, nn
  real(dp) :: deltax, center(2)
  real(dp), dimension(2,2) :: rotmat


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(10:11) = [ gauss, gaussb ]

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0
  coefficients%r(49:) = a
  coefficients%r(50:) = b


! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

! create mesh for line source

  nn = 2*ne+1

  call mesh_skeleton ( mesh_source, nnodes=nn, nelem=ne, elshape=2, ndim=2 )

  deltax = length/ne/2
  mesh_source%coor(:,1) = [ ( deltax*(i-1), i=1,nn ) ] + 0.5_dp - length/2
  mesh_source%coor(:,2) = 0.5_dp

  do elem = 1, mesh_source%nelem
    mesh_source%topology(1)%a(:,elem) = [ 2*elem-1, 2*elem, 2*elem + 1 ]
  end do

  if ( rotate ) then

!   rotate object

    rotmat(:,1) = [  cos(angle), sin(angle) ]
    rotmat(:,2) = [ -sin(angle), cos(angle) ]
    center = [ 0.5_dp, 0.5_dp ]

    do i = 1, mesh_source%nnodes
       mesh_source%coor(i,:) = center + [ shiftx,shifty] +  &
             matmul ( rotmat, mesh_source%coor(i,:) - center )
    end do

!   use large number of `mid-point' integration points

    call add_to_mesh ( mesh, object='mesh', objectmesh=mesh_source, &
      topology=.true., intrule=1, nsubint=nsubint )

  else

    call add_to_mesh ( mesh, object='mesh', objectmesh=mesh_source, &
      topology=.true., intrule=3, nsubint=1 )

  end if

  call fill_mesh_parts ( mesh )

! plot mesh and objects

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )
  plot_options%objectpointcolor=4
  plot_options%objectpointsize=0.4
  call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true. )

! problem definition

  call create_input_probdef ( mesh, input_probdef )

  input_probdef%elementdof(1)%a = 1

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, curve1=1, curve2=4, value=0._dp )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=poisson_elem, &
    coefficients=coefficients, buildvector=.false. )

  rhsd%u = 0

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub1=poisson_source, &
    coefficients=coefficients, object=1, addmatvec=.true. )

! this should be added only for non-zero essential bc, but it we do it anyway

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  print *, maxval( sol%u )

! plot solution

  call plot_color_contour ( plot_options, mesh, problem, filename='sol.fig', &
    sysvector=sol )
  call plot_objects ( plot_options, mesh, 'sol.fig', append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh, mesh_source )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program poisson9

