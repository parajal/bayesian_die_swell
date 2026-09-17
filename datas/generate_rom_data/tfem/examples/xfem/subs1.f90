
! Some common subroutines for the circular object in fic_dom1 and xfem1.

module subs1_m

  use tfem_elem_m
  use math_defs_m

  implicit none

  integer :: ndf = 9
  integer :: ndfb = 3, nodalpb = 3
  integer :: ndflb = 2  ! number of degrees of freedom of the Lag. Mult.
  real(dp) :: rp = 1, xpc(2) = 0

contains

! objectscoor defines the coordinates of the objects

  subroutine objectscoor ( objectnr, coor )
    integer, intent(in) :: objectnr
    real(dp), dimension(:,:), intent(inout) :: coor

    integer :: np, i
    real(dp) :: p(size(coor,1))

    np = size(coor,1)
    p = [(2*pi/np*(i-1),i=1,np)]

    coor(:,1) = rp * sin(p) + xpc(1)
    coor(:,2) = rp * cos(p) + xpc(2)

  end subroutine objectscoor

! elementc is the element subroutine for the constraints on the objects

  subroutine elementc ( mesh, problem, constr, elem, node, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemmat2, elemmatadd, &
    elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: object, i, j
    real(dp) :: phi(1,ndf), xr(1,2), psi(1,ndflb), x(nodalpb,2)
    real(dp) :: theta(1,ndfb), dtheta(1,ndfb,1), dxdxi(1,2), curvel(1)


    object = problem%constraints(constr)%object

!   reference coordinates of the integration point (node) of the element (elem)

    xr(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)

!   shape function of the velocity in the integration point

    call shape_quad_Q2 ( xr, phi )

!   shape function of the Langrangian multiplier in the integration point

    call shape_line_P1 ( mesh%objects(object)%xig(node:node,1), psi )

!   shape function of the quadratic curve in the integration point

    call shape_line_P2 ( mesh%objects(object)%xig(node:node,1), theta, &
      dtheta(:,:,1) )

!   compute deformed element

    call get_coordinates_object ( mesh, elem, x, object )

    call isoparametric_deformation_curve ( x, dtheta(:,:,1), dxdxi, curvel )

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

!   two constraints (vectorial)
!
!      u  = 0
!      -    -
!
!   or in components
!
!      u = 0
!      v = 0
!
!   the matrix A therefore becomes
!
!     A =  [ phi    0 ]
!          [   0  phi ]
!
!   In the weak version the equations are multiplied by the test function psi
!   and integrated over the element.
!   NOTE: this is the contribution of a single integration point (node) in
!   a single element (elem).

    if ( matrix ) then

      elemmat = 0

      do i = 1, ndflb
        do j = 1, ndf
          elemmat(i,j) = psi(1,i) * phi(1,j) * &
                                   curvel(1) * mesh%objects(object)%wg(node)
        end do
      end do
      elemmat(ndflb+1:2*ndflb,ndf+1:2*ndf) = elemmat(1:ndflb,1:ndf)

      elemmatadd = 0

    end if

  end subroutine elementc

end module subs1_m

