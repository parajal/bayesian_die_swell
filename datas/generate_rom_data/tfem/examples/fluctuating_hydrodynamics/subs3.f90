module subs3_m

  use tfem_elem_m

  implicit none

  integer, parameter :: nparticles=1
  real(dp) :: xp(nparticles,2) = 0, rp(nparticles) = 0, dragv(2), spr_const, &
    deltat

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

      if ( first ) then

        elemvecadd(1) = -sqrt(deltat)*spr_const*(xp(1,1))+dragv(1)
        elemvecadd(2) = +dragv(2)
        elemvecadd(3) = 0.0_dp

      end if

    end if

!   connection through collocation

    elemmat(1,1) = 1
    elemmat(1,2) = 0

    elemmat(2,1) = 0
    elemmat(2,2) = 1

    elemmatadd(1,:) = [ -1._dp, 0._dp, -xr(1,2) + xp(1,2) ]
    elemmatadd(2,:) = [ 0._dp, -1._dp, xr(1,1) - xp(1,1) ]

  end subroutine elementc

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

end module subs3_m
