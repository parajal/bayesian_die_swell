module subs_m

  use kind_defs_m
  use tfem_elem_m

  implicit none

  real(dp) :: eta1l, eta2l
  real(dp), allocatable :: xp(:,:), xpn(:,:), rp(:)

  integer :: ipart

contains

! Internal element routine to fill the viscosity in the Gauss points as a
! function of the composition c

  subroutine fill_eta_gauss_c ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use generalized_stokes_elements_m
    use generalized_stokes_globals_m
    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_generalized_stokes ( mesh, coefficients, elgrp )

      allocate ( wg(ninti), xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim), etag(ninti) )
      allocate ( c_n(ndf), c_ng(ninti) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_coordinates ( x, phi, xg )

!   determine coefficient by a function of c

    call get_vector ( mesh, oldvectors%p(2)%p, oldvectors%v(2)%p, elgrp, &
      elem, c_n, layer=layer )

    c_ng = matmul ( phi, c_n )

!   function to determine the local viscosity

    etag = eta1l * ((c_ng+1.0_dp)/2.0_dp) - eta2l * ((c_ng-1.0_dp)/2.0_dp)

    call put_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, r1=etag )

    if ( last ) then

!     last element in this group

      deallocate ( wg, xig, phi, x )
      deallocate ( xg, etag )
      deallocate ( c_n, c_ng )

    end if

  end subroutine fill_eta_gauss_c

! element routine for a freely floating particle

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

    elemmatadd(1,:) = [ -1._dp, 0._dp, -xr(1,2) + xp(ipart,2) ]
    elemmatadd(2,:) = [ 0._dp, -1._dp, xr(1,1) - xp(ipart,1) ]

  end subroutine elementc

end module subs_m
