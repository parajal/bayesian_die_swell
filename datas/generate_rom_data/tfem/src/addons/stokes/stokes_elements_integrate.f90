
! Copyright (C) 2014-2020 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! This module contains element subroutines for integration of fields of a Stokes
! problem over the internal domain and the boundary.
! The subroutines are to be used with "integrate" and
! "integrate_boundary_elements" in the postprocessing module.
!
! The subroutines for the internal domain are:
! - stokes_integrate_pressure: integral of the pressure over the domain
! - stokes_integrate_stress: integral of extra stress 2*eta*D over the domain
! - stokes_integrate_total_stress: integral of stress 2*eta*D-pI over the domain
! - stokes_integrate_volume: area/volume of the domain
! - stokes_L2_norm: L2-norm for velocity and pressure
!
! The subroutines for boundaries (geometries) are:
!
! - stokes_integrate_area_geometry: length (2D), area(3D) of geometry
! - stokes_integrate_area_tensor: integral of nn, nn-I/ndim or nn-I
! - stokes_integrate_pressure_geometry: integral of (pI.n)x=pnx
! - stokes_integrate_stress_geometry: integral of (2*eta*D.n)x

module stokes_elements_integrate_m

  use tfem_elem_m
  use stokes_set_globals_m
  use stokes_elements_generic_m, only: set_stokes_shape_function_global

  implicit none

  interface stokes_L2_error_norm
    module procedure stokes_L2_norm
  end interface stokes_L2_error_norm

contains


! Element routine for the integral of the pressure over the domain.
! For use with integrate.

  subroutine stokes_integrate_pressure ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec )

    use stokes_globals_m

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

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!     allocate arrays

      allocate ( wg(ninti), detF(ninti) )
      allocate ( work(ninti), xg(ninti,ndim) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), psi(ninti,ndfp) )
      allocate ( x(nodalp,ndim), dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim) )
      allocate ( u(ndfp) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      if ( coefficients%i(62) == 0 ) then
        call set_shape_function ( shapefuncp, xig, psi ) ! pressure
      end if

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    xg = matmul ( phi, x(1:ndf,:) )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   get pressure

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqpress], layer=layer )

    work = matmul ( psi, u )

!   integration of pressure over the element

    elemvec = sum ( work * detF * wg )

    if ( last ) then

!     last element

      deallocate ( wg, detF )
      deallocate ( work, xg )
      deallocate ( xig, phi, psi )
      deallocate ( x, dphi, F )
      deallocate ( Finv )
      deallocate ( u )

    end if

  end subroutine stokes_integrate_pressure


! Element routine for the integral of the stokes extra-stress over the domain.
! For use with integrate.

  subroutine stokes_integrate_stress ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j
    integer, save :: ns
    real(dp) :: eta

    if ( first ) then

!     first element on this curve

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )

      ns = ndim*(ndim+1)/2

!     allocate arrays

      allocate ( wg(ninti), detF(ninti) )
      allocate ( u(ndim*ndf) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( uvector(ndf,ndim), gradu(ninti,ndim,ndim) )
      allocate ( work2(ninti,ns) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      write(*,'(/a/a/)') 'Error stokes_integrate_stress:', &
        ' Axisymmetric coordinates not yet implemented'
      stop
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ndim] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    eta = coefficients%r(1)

    if ( ndim == 2 ) then

      work2(:,1) = 2 * eta * gradu(:,1,1)                ! D_xx
      work2(:,2) = eta * (gradu(:,1,2) + gradu(:,2,1) )  ! D_xy
      work2(:,3) = 2 * eta * gradu(:,2,2)                ! D_yy

    else

      work2(:,1) = 2 * eta * gradu(:,1,1)                ! D_xx
      work2(:,2) = eta * ( gradu(:,1,2) + gradu(:,2,1) ) ! D_xy
      work2(:,3) = eta * ( gradu(:,1,3) + gradu(:,3,1) ) ! D_xz
      work2(:,4) = 2 * eta * gradu(:,2,2)                ! D_yy
      work2(:,5) = eta * ( gradu(:,2,3) + gradu(:,3,2) ) ! D_yz
      work2(:,6) = 2 * eta * gradu(:,3,3)                ! D_zz

    end if

!   integrate: 2*eta*D

    do i = 1, ns
      elemvec(i) = sum ( work2(:,i) * detF * wg ) ! 2*eta*Dij (i,j = 1..ns)
    end do

   if ( last ) then

!     last element on this curve

      deallocate ( detF, u )
      deallocate ( phi, x )
      deallocate ( uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( work2, wg, xig )

    end if

  end subroutine stokes_integrate_stress


! Element routine for the integral of the total stress over the domain.
! For use with integrate.

  subroutine stokes_integrate_total_stress ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j
    integer, save :: ns
    real(dp) :: eta

    if ( first ) then

!     first element on this curve

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )

      ns = ndim*(ndim+1)/2

!     allocate arrays

      allocate ( wg(ninti), detF(ninti) )
      allocate ( work(ninti), xg(ninti,ndim) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), psi(ninti,ndfp) )
      allocate ( x(nodalp,ndim), dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( uvector(ndf,ndim), gradu(ninti,ndim,ndim) )
      allocate ( work2(ninti,ns), u(ndim*ndf), pr(ndfp) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      if ( coefficients%i(62) == 0 ) then
        call set_shape_function ( shapefuncp, xig, psi ) ! pressure
      end if

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      write(*,'(/a/a/)') 'Error stokes_integrate_total_stress:', &
        ' Axisymmetric coordinates not yet implemented'
      stop
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( any ( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      xg = matmul ( phi, x )
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ndim] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    eta = coefficients%r(1)

    if ( ndim == 2) then

      work2(:,1) = 2 * eta * gradu(:,1,1)                ! D_xx
      work2(:,2) = eta * (gradu(:,1,2) + gradu(:,2,1) )  ! D_xy
      work2(:,3) = 2 * eta * gradu(:,2,2)                ! D_yy

    else

      work2(:,1) = 2 * eta * gradu(:,1,1)                ! D_xx
      work2(:,2) = eta * ( gradu(:,1,2) + gradu(:,2,1) ) ! D_xy
      work2(:,3) = eta * ( gradu(:,1,3) + gradu(:,3,1) ) ! D_xz
      work2(:,4) = 2 * eta * gradu(:,2,2)                ! D_yy
      work2(:,5) = eta * ( gradu(:,2,3) + gradu(:,3,2) ) ! D_yz
      work2(:,6) = 2 * eta * gradu(:,3,3)                ! D_zz

    end if

!   get pressure

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, pr, &
      physq=[physqpress], layer=layer )

    work = matmul ( psi, pr )

    if ( ndim == 2 ) then
      work2(:,1) = work2(:,1) - work
      work2(:,3) = work2(:,3) - work
    else if ( ndim == 3 ) then
      work2(:,1) = work2(:,1) - work
      work2(:,4) = work2(:,4) - work
      work2(:,6) = work2(:,6) - work
    end if

!   integrate: -p + 2*eta*D

    do i = 1, ns
      elemvec(i) = sum ( work2(:,i) *  detF * wg )
    end do

    if ( last ) then

!     last element on this curve

      deallocate ( detF )
      deallocate ( phi, psi, x )
      deallocate ( uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( work, xg, wg, xig )
      deallocate ( work2, u, pr )

    end if

  end subroutine stokes_integrate_total_stress


! Element routine for the area(2D) or volume (3D) of a domain.
! For use with integrate.

  subroutine stokes_integrate_volume ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec )

    use stokes_globals_m

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

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

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

  end subroutine stokes_integrate_volume


! Element routine for the (square of the) L2 norm of the velocity
! and pressure over the domain.
! The exact solution for velocity and pressure can be given via vfunc and func
! in oldvectors. In that case, the L2-norm of the difference will be computed.
! The pressure level, to be subtracted from the pressure solution, must be
! supplied via coefficients.
! For use with integrate.

  subroutine stokes_L2_norm ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec )

    use stokes_globals_m

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

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

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
      if ( coefficients%i(62) == 0 ) then
        call set_shape_function ( shapefuncp, xig, psi ) ! pressure
      end if

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    xg = matmul ( phi, x(1:ndf,:) )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   get velocity

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

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

!   get pressure

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

!   integration of L2-norm squared over the element

    work1 = sum ( ( work2 - u_exact ) ** 2, dim=2 )

    elemvec(1) = sum ( work1 * detF * wg )
    elemvec(2) = sum ( ( work - p_exact ) ** 2 * detF * wg )

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

  end subroutine stokes_L2_norm


! Element routine for the length(2D) or area (3D) of a geometry.
! For use with integrate_boundary_elements.

  subroutine stokes_integrate_area_geometry ( mesh, problem, geometry, elem, &
    first, last, coefficients, oldvectors, elemvec )

    use stokes_globals_m

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

      call set_globals_stokes_vp_boun ( mesh, coefficients, &
        ndimr=mesh%ndim-1, geometry=geometry, maxvel3D=1 )

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

  end subroutine stokes_integrate_area_geometry


! Element routine for the integral of nn (area tensor), nn-I/ndim or nn-I.
! For use with integrate_boundary_elements.

  subroutine stokes_integrate_area_tensor ( mesh, problem, geometry, &
    elem, first, last, coefficients, oldvectors, elemvec )

    use stokes_globals_m

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

      call check ( coefficients, 'stokes_integrate_area_tensor', &
        indexarray=[59], minimum=[0], maximum=[2] )

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, &
        ndimr=mesh%ndim-1, geometry=geometry, maxvel3D=1 )

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
      write(*,'(/a/a/)') 'Error stokes_integrate_area_tensor:', &
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
      call errormsg_case_default ( 'stokes_integrate_area_tensor', &
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

  end subroutine stokes_integrate_area_tensor


! Element routine for the integral of (pI.n)x=pnx
! For use with integrate_boundary_elements.
! NOTE: the default direction of n is determined by the direction of
! the elements on the geometry. Use coefficients%i(36) to "flip" the normal.

  subroutine stokes_integrate_pressure_geometry ( mesh, problem, geometry, &
    elem, first, last, coefficients, oldvectors, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: geometry, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, orientation

    if ( first ) then

!     first element on this geometry

      call check ( coefficients, 'stokes_integrate_pressure_geometry', &
        indexarray=[36], minimum=[-1], maximum=[1] )

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, &
        ndimr=mesh%ndim-1, geometry=geometry )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti), normal(ninti,ndim) )
      allocate ( xig(ninti,ndim-1), x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( pr(ndf), work(ninti) )
      allocate ( tmp(ndim,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, ndimr=ndim-1, &
      geometry=geometry )

    call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, normal )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      write(*,'(/a/a/)') 'Error stokes_integrate_pressure_geometry:', &
        ' Axisymmetric coordinates not yet implemented'
      stop
    end if

!   get pressure

    call get_vector_geometry ( mesh, problem, oldvectors%v(2)%p, elem, pr, &
      ndimr=ndim-1, geometry=geometry, layer=layer )

    work = matmul ( phi, pr )

!   set normal

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    normal = real ( orientation, kind=dp ) * normal

!   evaluate in each integration point: (p*I.n)x=pnx and integrate

    do i = 1, ndim
      do j = 1, ndim
        tmp(i,j) = sum ( work * normal(:,i) * xg(:,j) * surfl * wg )
      end do
    end do

!   element vector

    elemvec = reshape ( transpose(tmp), [ndim**2] )

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl, normal )
      deallocate ( xig, x, xg )
      deallocate ( phi, dphi, dxdxis)
      deallocate ( pr, work )
      deallocate ( tmp )

    end if

  end subroutine stokes_integrate_pressure_geometry


! Element routine for the integral of (2*eta*D.n)x
! For use with integrate_boundary_elements.
! NOTE: the default direction of n is determined by the direction of
! the elements on the geometry. Use coefficients%i(36) to "flip" the normal.

  subroutine stokes_integrate_stress_geometry ( mesh, problem, geometry, &
    elem, first, last, coefficients, oldvectors, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: geometry, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, ip
    integer, save :: ns

    if ( first ) then

!     first element on this geometry

      call check ( coefficients, 'stokes_integrate_stress_geometry', &
        indexarray=[36], minimum=[-1], maximum=[1] )

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, &
        ndimr=mesh%ndim-1, geometry=geometry )

      ns = ndim*(ndim+1)/2

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti), normal(ninti,ndim) )
      allocate ( xig(ninti,ndim-1), x(nodalp,ndim), xg(ninti,ndim))
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( st(ndf*ns), work2(ndf,ns), work4(ninti,ndim), tmp(ndim,ndim) )
      allocate ( tauten(ninti,ndim,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, ndimr=ndim-1, &
      geometry=geometry )

    call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, normal )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      write(*,'(/a/a/)') 'stokes_integrate_stress_geometry:', &
        ' Axisymmetric coordinates not yet implemented'
      stop
    end if

!   stokes stress tensor

    call get_vector_geometry ( mesh, problem, oldvectors%v(1)%p, elem, st, &
     ndimr=ndim-1, geometry=geometry, layer=layer )

    work2 = reshape ( st, [ndf,ns] )

    if ( ndim == 2 ) then

      tauten(:,1,1) = matmul ( phi, work2(:,1) )
      tauten(:,1,2) = matmul ( phi, work2(:,2) )
      tauten(:,2,1) = tauten(:,1,2)
      tauten(:,2,2) = matmul ( phi, work2(:,3) )

    else if ( ndim == 3 ) then

      tauten(:,1,1) = matmul ( phi, work2(:,1) )
      tauten(:,1,2) = matmul ( phi, work2(:,2) )
      tauten(:,1,3) = matmul ( phi, work2(:,3) )
      tauten(:,2,1) = tauten(:,1,2)
      tauten(:,2,2) = matmul ( phi, work2(:,4) )
      tauten(:,2,3) = matmul ( phi, work2(:,5) )
      tauten(:,3,1) = tauten(:,1,3)
      tauten(:,3,2) = tauten(:,2,3)
      tauten(:,3,3) = matmul ( phi, work2(:,6) )

    end if

!   evaluate in each integration point: 2*eta*D.n

    do ip = 1, ninti
      work4(ip,:) = matmul( tauten(ip,:,:), normal(ip,:) )
    end do

!   evaluate in each integration point: (2*eta*D.n)x and integrate

    do i = 1, ndim
      do j = 1, ndim
        tmp(i,j) = sum ( work4(:,i) * xg(:,j) * surfl * wg )
      end do
    end do

!   element vector

    elemvec = reshape ( transpose(tmp), [ndim**2] )

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl, normal )
      deallocate ( xig, x, xg )
      deallocate ( phi, dphi, dxdxis)
      deallocate ( st, work2, work4, tmp )
      deallocate ( tauten )

    end if

  end subroutine stokes_integrate_stress_geometry

end module stokes_elements_integrate_m

