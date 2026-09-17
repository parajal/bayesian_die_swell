
! Copyright (C) 2020-2020 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! This module contains element subroutines for integration of fields of a
! linear elastic problem over the internal domain and the boundary.
! The subroutines are to be used with "integrate" and
! "integrate_boundary_elements" in the postprocessing module.
!
! The subroutines for the internal domain are:
! - linear_elastic_integrate_stress: integral of Cauchy stress over the domain
! - linear_elastic_integrate_volume: area/volume of the domain
! - linear_elastic_L2_norm: L2-norm for displacement and pressure
!
! The subroutines for boundaries (geometries) are:
!
! - linear_elastic_integrate_area_geometry: length (2D), area(3D) of geometry
! - linear_elastic_integrate_area_tensor: integral of nn, nn-I/ndim or nn-I
! - linear_elastic_integrate_stress_geometry: integral of (sigma.n)x

module linear_elastic_elements_integrate_m

  use tfem_elem_m
  use linear_elastic_set_globals_m
  use linear_elastic_elements_generic_m, &
            only: set_linear_elastic_shape_function_global, le_stress_tensor
  use linear_elastic_material_m

  implicit none

contains


! Element routine for the integral of the total stress over the domain.
! For use with integrate.

  subroutine linear_elastic_integrate_stress ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec )

    use linear_elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j
    integer, save :: ns

    if ( first ) then

!     first element on this curve

!     set globals

      call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
        maxdisp3D=1 )

      call set_linear_elastic_material ( coefficients, coorsys, lemodel )

      select case (ndim)
      case (2)
        ns = 4
      case (3)
        ns = 6
      case default
        call errormsg_case_default ( 'linear_elastic_integrate_stress', &
          'ndim', int_value=ndim )
      end select

!     allocate arrays

      allocate ( wg(ninti), detF(ninti) )
      allocate ( xg(ninti,ndim) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), psi(ninti,ndfp) )
      allocate ( x(nodalp,ndim), dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( uvector(ndf,ndim), gradu(ninti,ndim,ndim) )
      allocate ( work2(ninti,ns), u(ndim*ndf), prs(ndfp), pr(ninti) )
      allocate ( dTemp(ninti) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      if (  ndfp > 0 .and. coefficients%i(62) == 0 ) then
        call set_shape_function ( shapefuncp, xig, psi ) ! pressure
      end if

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    xg = matmul ( phi, x )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqdisp], layer=layer )

    uvector = reshape ( u, [ndf,ndim] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( disp3D == 1 ) then

      gradu(:,:,3) = 0

      if ( coorsys == 1 ) then
        where ( xg(:,2) < 1e-10_dp )
          gradu(:,2,3) = -gradu(:,3,2)       ! -du_theta / dr (r=0)
          gradu(:,3,3) = gradu(:,2,2)        ! du_r / dr (r=0)
        elsewhere
          gradu(:,2,3) = -uvector(:,3) / xg(:,2)  ! -u_theta / r
          gradu(:,3,3) = uvector(:,2) / xg(:,2)   ! u_r / r
        end where
      end if

    end if

    if ( ndfp > 0 ) then

!     mixed elements

      if ( any ( coefficients%i(62) == [1,2] ) ) then
!       global shape function for pressure
        call set_linear_elastic_shape_function_global ( shapefuncp, x, xg, psi,&
          coefficients%i(62)==1 )
      end if

      call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, prs, &
        physq=[physqpress], layer=layer )

      pr = matmul ( psi, prs )

    end if

    if ( lemodel%thermal_expansion .and. ndfp == 0 ) then

!     thermal expansion: get temperature change

      call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
        elgrp, elem, choice=coefficients%i(9), value=coefficients%r(4), &
        func=coefficients%func1(1)%p, funcnr=coefficients%i(8), x=xg, &
        indx_v=1, layer=layer, phi=phi, indx_e=1, coef=dTemp )

    end if

!   compute stress tensor

    call le_stress_tensor ( xg, work2 )

!   integrate sigma

    do i = 1, ns
      elemvec(i) = sum ( work2(:,i) * detF * wg )
    end do

    if ( last ) then

!     last element on this curve

      deallocate ( detF )
      deallocate ( phi, psi, x )
      deallocate ( uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( xg, wg, xig )
      deallocate ( work2, u, prs, pr )
      deallocate ( dTemp )

    end if

  end subroutine linear_elastic_integrate_stress


! Element routine for the area(2D) or volume (3D) of a domain.
! For use with integrate.

  subroutine linear_elastic_integrate_volume ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec )

    use linear_elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    if ( first ) then

!     first element

!     set globals

      call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
        maxdisp3D=1 )

!     allocate arrays

      allocate ( wg(ninti), detF(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), xg(ninti,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      xg = matmul ( phi, x(1:ndf,:) )
      detF = 2 * pi * xg(:,2) * detF
    end if

!   integration of over the element

    elemvec = sum ( detF * wg )

    if ( last ) then

!     last element

      deallocate ( wg, detF )
      deallocate ( xig, phi, x )
      deallocate ( dphi, F )
      deallocate ( Finv, xg )

    end if

  end subroutine linear_elastic_integrate_volume


! Element routine for the (square of the) L2 norm of the displacement
! and pressure over the domain. The exact solution for displacement and
! pressure can be given via vfunc and func in oldvectors. In that case,
! the L2-norm of the difference will be computed.
! The pressure level, to be subtracted from the pressure solution, must be
! supplied via coefficients.
! For use with integrate.

  subroutine linear_elastic_L2_norm ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec )

    use linear_elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: ip

    if ( first ) then

!     first element

!     set globals

      call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
        maxdisp3D=1 )

!     allocate arrays

      allocate ( wg(ninti), detF(ninti) )
      allocate ( work(ninti), work2(ninti,ncompu), work1(ninti), &
                 tmp(ndf,ncompu) )
      allocate ( xg(ninti,ndim) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), psi(ninti,ndfp) )
      allocate ( x(nodalp,ndim), dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim) )
      allocate ( pr(ndfp), u(ndf*ncompu) )
      allocate ( p_exact(ninti), u_exact(ninti,ncompu) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      if ( ndfp > 0 .and. coefficients%i(62) == 0 ) then
        call set_shape_function ( shapefuncp, xig, psi ) ! pressure
      end if

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    xg = matmul ( phi, x(1:ndf,:) )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( ndfp > 0 .and. any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_linear_elastic_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   get displacement

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqdisp], layer=layer )

    tmp = reshape ( u, [ndf,ncompu] )

    work2 = matmul ( phi, tmp )

    if ( coefficients%i(63) > 0 ) then
      do ip = 1, ninti
        u_exact(ip,:) = &
                   coefficients%vfunc ( ncompu, coefficients%i(63), xg(ip,:) )
      end do
    else
      u_exact = 0
    end if

!   integration of L2-norm squared over the element

    work1 = sum ( ( work2 - u_exact ) ** 2, dim=2 )

    elemvec(1) = sum ( work1 * detF * wg )

    if ( ndfp > 0 ) then

!     get pressure

      call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, pr, &
        physq=[physqpress], layer=layer )

      work = matmul ( psi, pr )

      work = work - coefficients%r(23)

      if ( coefficients%i(64) > 0 ) then
        do ip = 1, ninti
          p_exact(ip) = coefficients%func ( coefficients%i(64), xg(ip,:) )
        end do
      else
        p_exact = 0
      end if

      elemvec(2) = sum ( ( work - p_exact ) ** 2 * detF * wg )

    end if

    if ( last ) then

!     last element

      deallocate ( wg, detF )
      deallocate ( work, work2, work1, tmp )
      deallocate ( xg )
      deallocate ( xig, phi, psi )
      deallocate ( x, dphi, F )
      deallocate ( Finv )
      deallocate ( pr, u )
      deallocate ( p_exact, u_exact )

    end if

  end subroutine linear_elastic_L2_norm


! Element routine for the length(2D) or area (3D) of a geometry.
! For use with integrate_boundary_elements.

  subroutine linear_elastic_integrate_area_geometry ( mesh, problem, &
    geometry, elem, first, last, coefficients, oldvectors, elemvec )

    use linear_elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: geometry, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    if ( first ) then

!     first element on this curve

!     set globals

      call set_globals_linear_elastic_up_boun ( mesh, coefficients, &
        ndimr=mesh%ndim-1, geometry=geometry, maxdisp3D=1 )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti), xg(ninti,ndim) )
      allocate ( xig(ninti,ndim-1), x(nodalp,ndim) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, ndimr=ndim-1, &
      geometry=geometry )

    call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl )

    if ( coorsys == 1 ) then
      xg = matmul ( phi, x )
      surfl = 2 * pi * xg(:,2) * surfl
    end if

!   integrate over the element

    elemvec(1) = sum ( surfl * wg )

    if ( last ) then

!     last element on this curve

      deallocate ( wg, surfl, xg )
      deallocate ( xig, x )
      deallocate ( phi, dphi, dxdxis )

    end if

  end subroutine linear_elastic_integrate_area_geometry


! Element routine for the integral of nn (area tensor), nn-I/ndim or nn-I.
! For use with integrate_boundary_elements.

  subroutine linear_elastic_integrate_area_tensor ( mesh, problem, geometry, &
    elem, first, last, coefficients, oldvectors, elemvec )

    use linear_elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: geometry, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j
    real(dp) :: fac

    if ( first ) then

!     first element on this geometry

      call check ( coefficients, 'linear_elastic_integrate_area_tensor', &
        indexarray=[59], minimum=[0], maximum=[2] )

!     set globals

      call set_globals_linear_elastic_up_boun ( mesh, coefficients, &
        ndimr=mesh%ndim-1, geometry=geometry, maxdisp3D=1 )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti), normal(ninti,ndim) )
      allocate ( xig(ninti,ndim-1), x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( tmp(ndim,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, ndimr=ndim-1, &
      geometry=geometry )

    call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, normal )

    if ( coorsys == 1 ) then
      xg = matmul ( phi, x )
      write(*,'(/a/a/)') 'Error linear_elastic_integrate_area_tensor:', &
        ' Axisymmetric coordinates not yet implemented'
      stop
    end if

!   diagonal blocks of nn

    select case ( coefficients%i(59) )
    case(0)
      fac = 0
    case(1)
      fac = 1._dp/ndim
    case(2)
      fac = 1
    case default
      call errormsg_case_default ( 'linear_elastic_integrate_area_tensor', &
        'coefficients%i(59)', int_value=coefficients%i(59) )
    end select

    do i = 1, ndim
      tmp(i,i) = sum ( ( normal(:,i)**2 - fac ) * surfl * wg )
    end do

!   off-diagonal blocks of nn (upper triangle only)

    do i = 1, ndim-1
      do j = i+1, ndim
        tmp(i,j) = sum ( normal(:,i) * normal(:,j) * surfl * wg )
      end do
    end do

!   fill elemvec

    if ( ndim == 2 ) then

      elemvec(1) = tmp(1,1)  ! xx
      elemvec(2) = tmp(1,2)  ! xy
      elemvec(3) = tmp(2,2)  ! yy

    else if ( ndim == 3 ) then

      elemvec(1) = tmp(1,1) ! xx
      elemvec(2) = tmp(1,2) ! xy
      elemvec(3) = tmp(1,3) ! xz
      elemvec(4) = tmp(2,2) ! yy
      elemvec(5) = tmp(2,3) ! yz
      elemvec(6) = tmp(3,3) ! zz

    end if

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl, normal )
      deallocate ( xig, x, xg )
      deallocate ( phi, dphi, dxdxis )
      deallocate ( tmp )

    end if

  end subroutine linear_elastic_integrate_area_tensor

end module linear_elastic_elements_integrate_m

