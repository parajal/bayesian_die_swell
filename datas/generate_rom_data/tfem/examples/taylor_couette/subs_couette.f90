module subs_couette_m

  use tfem_elem_m

  implicit none

  logical :: set_force = .false.
  real(dp), allocatable, dimension(:) :: xp, force
  real(dp) :: rp

contains

  subroutine elementc ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: ncurveconstr
    real(dp) :: xr(1,2)

    ncurveconstr = problem%constraints(constr)%geometry1

    xr(1,:) = mesh%coor(mesh%curves(ncurveconstr)%nodes(node),:)

!   set shape function in the point

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

!   connection through collocation

    elemmat(1,1) = 1
    elemmat(1,2) = 0

    elemmat(2,1) = 0
    elemmat(2,2) = 1

    elemmatadd(1,:) = [ -1._dp, 0._dp, -xr(1,2) + xp(2) ]
    elemmatadd(2,:) = [ 0._dp, -1._dp, xr(1,1) - xp(1) ]

  end subroutine elementc


! subroutine for freely floating particles (3D)

  subroutine elementc_3D ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: nsurfaceconstr
    real(dp) :: xr(1,3), r(3)

    nsurfaceconstr = problem%constraints(constr)%geometry1

    xr(1,:) = mesh%coor(mesh%surfaces(nsurfaceconstr)%nodes(node),:)

!   set shape function in the point

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0
      if ( set_force ) then
        elemvecadd(1:3) = force
        set_force = .false.
      end if

    end if

!   connection through collocation
!   three constraints (vectorial)
!      _   _      _     _   _
!      u - up - omega x r = 0
!
!   where r is the vector from the centerpoint of the particle to the boundary
!
!   or in components
!
!      u - up - omega_y x rz + omega_z x ry= 0
!      v - vp - omega_z x rx + omega_x x rz= 0
!      w - wp - omega_x x ry + omega_y x rx= 0
!
!   the matrix A therefore becomes
!
!     A =  [ phi    0     0]
!          [   0  phi     0]
!          [   0    0   phi]
!
!   and the matrix A_add of the additional unknowns (up,vp,wp,omega_x,
!   omega_y, omega_z)
!
!     A_add =  [ -1  0  0   0 -rz   ry ]
!              [  0 -1  0  rz   0  -rx ]
!              [  0  0 -1 -ry  rx    0 ]

    if ( matrix ) then

      elemmat(:,:) = 0._dp
      elemmat(1,1) = 1._dp
      elemmat(2,2) = 1._dp
      elemmat(3,3) = 1._dp

      r = xr(1,:) - xp(:)

      elemmatadd(1,:) = [ -1._dp,   0._dp,   0._dp,  0._dp,  -r(3),   r(2) ]
      elemmatadd(2,:) = [  0._dp,  -1._dp,   0._dp,   r(3),  0._dp,  -r(1) ]
      elemmatadd(3,:) = [  0._dp,   0._dp,  -1._dp,  -r(2),   r(1),  0._dp ]

    end if

  end subroutine elementc_3D


  subroutine constraints_object_conn_P2 ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer, parameter :: ndf = 6

    integer :: object, object2
    real(dp) :: phi(1,ndf), phi2(1,ndf), xr(1,2), xr2(1,2)

!   connection through collocation on objects

    object = problem%constraints(constr)%object
    object2 = problem%constraints(constr)%object2

    xr(1,:) = mesh%objects(object)%refcoor(node,:)
    xr2(1,:) = mesh%objects(object2)%refcoor(node,:)

    call shape_triangle_P2 ( xr, phi )
    call shape_triangle_P2 ( xr2, phi2 )

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

      elemmat(1,1:ndf)       = phi(1,:)
      elemmat(1,ndf+1:2*ndf) = 0
      elemmat(2,1:ndf)       = 0
      elemmat(2,ndf+1:2*ndf) = phi(1,:)
      elemmat2(1,1:ndf)       = -phi2(1,:)
      elemmat2(1,ndf+1:2*ndf) = 0
      elemmat2(2,1:ndf)       = 0
      elemmat2(2,ndf+1:2*ndf) = -phi2(1,:)

    end if

  end subroutine constraints_object_conn_P2


  subroutine constraints_object_conn_P1 ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer, parameter :: ndf = 3

    integer :: object, object2, i
    real(dp) :: phi(1,ndf), phi2(1,ndf), xr(1,2), xr2(1,2)

!   connection through collocation on objects

    object = problem%constraints(constr)%object
    object2 = problem%constraints(constr)%object2

    xr(1,:) = mesh%objects(object)%refcoor(node,:)
    xr2(1,:) = mesh%objects(object2)%refcoor(node,:)

    call shape_triangle_P1 ( xr, phi )
    call shape_triangle_P1 ( xr2, phi2 )

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

      elemmat = 0
      do i = 1, size(elemmat,1)
        elemmat(i,ndf*(i-1)+1:i*ndf) = phi(1,:)
      end do
      elemmat2 = 0
      do i = 1, size(elemmat2,1)
        elemmat2(i,ndf*(i-1)+1:i*ndf) = -phi2(1,:)
      end do

    end if

  end subroutine constraints_object_conn_P1


! sample value of velocity in one node of the object

  subroutine stokes_sample_velocity1d ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, uvec )

    use stokes_elements_m
    use stokes_globals_m

!   WORKAROUND: dimension of velocity is different from dimension of mesh
    integer, parameter :: ndim2=2

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: uvec


    integer :: elgrp, elem
    real(dp) :: xr(1,mesh%ndim)

    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

    allocate ( phi(1,ndf), u(ndim2*ndf), tmp(ndf,ndim2) )

!   set shape function in the point

    call set_shape_function ( shapefunc, xr, phi )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    tmp = reshape ( u, [ndf,ndim2] )

    uvec = matmul( phi(1,:), tmp )

    deallocate ( phi, u, tmp )

  end subroutine stokes_sample_velocity1d


! sample value of conformation in one node of the object (standard)

  subroutine sample_conformation_tensor_std1d ( mesh, problem, object, &
    nodeobj, coefficients, oldvectors, cval )

    use viscoelastic_elements_m
    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: cval


    integer :: elgrp, elem, j, mode
    real(dp) :: xr(1,mesh%ndim)


    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   WORKAROUND: for 1d meshes, coorsys = 2 in stokes_globals
    coorsys = 0

    call set_viscoelastic_model ( coefficients )
    call set_globals_viscoelastic ( coefficients )

!   check coefficients

    call check ( coefficients, 'sample_conformation_tensor_std', &
      indexarray=[28], minimum=[0], maximum=[nmodes] )

    allocate ( theta(1,ndfc), c(ndfc,ncompc), posc(ndfc) )

!   set shape function in the point

    call set_shape_function ( shapefuncc, xr, theta )

    mode = max ( coefficients%i(28), 1 )

    call get_sysvector ( mesh, problem, oldvectors%s2(1)%p(1,mode), &
      elgrp, elem, c(:,1), posu=posc, layer=layer )

    do j = 2, ncompc
      c(:,j) = oldvectors%s2(1)%p(j,mode)%u(posc)
    end do

    cval = matmul( theta(1,:), c )

!   delete data

    call delete ( vemodel )

    deallocate ( theta, c, posc )

  end subroutine sample_conformation_tensor_std1d


end module subs_couette_m
