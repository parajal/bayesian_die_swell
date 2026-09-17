
! Copyright (C) 2005-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the constitutive equation using DG
!

module viscoelastic_elements_DG_2D_m

  use viscoelastic_elements_generic_m

  implicit none


contains


! Internal element routine for the constitutive equations: DG

  subroutine ce_dg_elem ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   select depending on integration scheme

    select case ( coefficients%i(25) )
      case(1) ! explicit DG, build matrix and right-hand side
        call ce_dg_elem1 ( mesh, problem, elgrp, elem, matrix, vector, &
          first, last, coefficients, oldvectors, elemmat, elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in ce_dg_elem:', &
        ' incorrect value for coefficients%i(25) = ', coefficients%i(25)
        stop
    end select

  end subroutine ce_dg_elem


! Internal element routine for the constitutive equations: DG
! Explicit time-integration

  subroutine ce_dg_elem1 ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip, m, side, elgrpnr, sidenr, elemnr
    integer :: elshape


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic ( coefficients )
      call set_globals_viscoelastic_DG ( coefficients )

!     coefficients

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

!     allocate arrays

      allocate ( wg(ninti), detF(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( dphidx(ninti,ndf,ndim) )
      allocate ( Finv(ninti,ndim,ndim) )
      allocate ( theta(ninti,ndfc), dtheta(ninti,ndfc,ndim) )
      allocate ( dthetadx(ninti,ndfc,ndim) )
      allocate ( ungradtheta(ninti,ndfc) )
      allocate ( ungradcn(ninti,ncompc,nmodes) )
      allocate ( u(ndim*ndf), tmp(ndf,ndim), gradu(ninti,ndim,ndim) )
      allocate ( xg(ninti,ndim) )

      allocate ( uvecn(ninti,ndim), cn(ndfc,ncompc) )
      allocate ( cng(ninti,ncompc,nmodes), gvecn(ninti,4+coorsys) )
      allocate ( fng(ninti,ncompc,nmodes) )

      allocate ( xigb(nintbc), wgb(nintbc) )
      allocate ( phib(nintbc,ndfb), dphib(nintbc,ndfb) )
      allocate ( xigs(nintbc,ndim,nsides), phis(nintbc,ndf,nsides) )
      allocate ( thetas(nintbc,ndfc,nsides) )
      allocate ( ugsn(nintbc,nsides), xgs(nintbc,ndim,nsides) )
      allocate ( dxdxi(nintbc,ndim), curvels(nintbc,nsides) )
      allocate ( normals(nintbc,ndim,nsides), us(nintbc,ndim) )
      allocate ( xs(nodalpb,ndim,nsides) )
      allocate ( cnjump(nintbc,ncompc,nmodes,nsides) )
      allocate ( cnside(ndfc,ncompc) )
      allocate ( work6(ndfc,ncompc,mode2-mode1+1) )

      elshape = mesh%element(elgrp)%elshape

      if ( mesh%element(elgrp)%globalshape == 'quadrilateral' .and. &
           elshape /= 6 ) then
        write(*,'(/a/a/)') 'Error in ce_dg_elem1:', &
          ' For quadrilaterals only elshape=6 is allowed '
        stop
      else if ( mesh%element(elgrp)%globalshape == 'triangle' .and. &
           elshape /= 4 .and. elshape /= 7 ) then
        write(*,'(/a/a/)') 'Error in ce_dg_elem1:', &
          ' For triangles only elshape=4 and 7 are allowed '
        stop
      end if

!     set Gauss integration for internal element and boundary integral

      call set_Gauss_integration ( globalshape, xig, wg )

      call Gauss_Legendre_line ( nintbc, xigb, wgb )

!     set shape functions

      call set_shape_function ( globalshape, xig, phi, dphi ) ! velocity
      call set_shape_function ( globalshape, xig, theta, dtheta ) ! c tensor

!     isoparametric shape of boundary
      call shape_line_P2 ( xigb, phib, dphib )

      call spread_to_sides_2D ( xigb, xigs )

      do side = 1, nsides
!       shape function on the side for velocity and c tensor
        call set_shape_function ( globalshape, xigs(:,:,side), phis(:,:,side) )
        call set_shape_function ( globalshape, xigs(:,:,side), thetas(:,:,side))
      end do

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call get_coordinates_sides ( mesh, elgrp, elem, xs )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )
    call shape_derivative ( dtheta, Finv, dthetadx )

    do side = 1, nsides
      call isoparametric_coordinates ( x, phis(:,:,side), xgs(:,:,side) )
      call isoparametric_deformation_curve ( xs(:,:,side), dphib, dxdxi, &
        curvels(:,side), normals(:,:,side) )
    end do

    if ( coorsys == 1 ) then
      curvels = 2 * pi * xgs(:,2,:) * curvels
    end if

!   get velocity vector

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel], layer=layer )

    tmp = reshape ( u, [ndf,ndim] )

    uvecn = matmul ( phi, tmp )

!   ugsn on the sides (normal velocity)

    do side = 1, nsides
      us = matmul ( phis(:,:,side), tmp )
      do ip = 1, nintbc
        ugsn(ip,side) = dot_product ( us(ip,:), normals(ip,:,side) )
      end do
    end do

!   un.grad operator

    do ip = 1, ninti
      ungradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uvecn(ip,:) )
    end do

!   compute velocity gradient vector

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), tmp )
    end do

    gvecn(:,1) = gradu(:,1,1)
    gvecn(:,2) = gradu(:,1,2)
    gvecn(:,3) = gradu(:,2,1)
    gvecn(:,4) = gradu(:,2,2)

    if ( coorsys == 1 ) then
      gvecn(:,5) = uvecn(:,2) / xg(:,2)
    end if

!   get conformation tensor at previous time step

    do m = mode1, mode2

!     single mode conformation
      call get_conformation ( mesh, problem, elgrp, elem, &
        oldvectors, cmode=cn, mode=m )

      cng(:,:,m) = matmul ( theta, cn )

!     un.grad cn term

      ungradcn(:,:,m) = matmul ( ungradtheta, cn )

!     cnjump on the sides

      do side = 1, nsides

        elgrpnr = mesh%sidelem(elgrp)%a(side,elem,1)

        if ( elgrpnr /= 0 ) then

!         side element found: get cn in side element

          elemnr = mesh%sidelem(elgrp)%a(side,elem,2)
          sidenr = mesh%sidelem(elgrp)%a(side,elem,3)

!         single mode conformation
          call get_conformation ( mesh, problem, elgrpnr, elemnr, &
            oldvectors, cmode=cnside, mode=m )

        end if

!       cnjump = cn(outside) - cn(inside) for each intgr point on this side

        do ip = 1, nintbc

          if ( ugsn(ip,side) < 0 ) then

            cnjump(ip,:,m,side) = - matmul ( thetas(ip,:,side), cn )

            if ( elgrpnr /= 0 ) then
              cnjump(ip,:,m,side) = cnjump(ip,:,m,side) &
                          + matmul ( thetas(nintbc-ip+1,:,sidenr), cnside )
            else
              cnjump(ip,:,m,side) = cnjump(ip,:,m,side) &
                + coefficients%cinflow ( ncompc, xgs(ip,:,side), ugsn(ip,side) )
            end if

          else

            cnjump(ip,:,m,side) = 0

          end if

        end do

      end do

    end do

!   matrix and vector

    if ( matrix ) then

!     mass matrix

      do i = 1, ndfc
        do j = 1, ndfc
          elemmat(i,j) = sum ( theta(:,i) * theta(:,j) * detF * wg )
        end do
      end do

    end if

    if ( vector ) then

!     viscoelastic rhs

      if ( logc ) then
        call rhs_viscoelastic_2D_log ( vemodel, gvecn, cng, fng, mode1=mode1, &
          mode2=mode2 )
      else
        call rhs_viscoelastic_2D ( vemodel, gvecn, cng, fng, mode1=mode1, &
          mode2=mode2 )
      end if

!     right-hand side of \dot c = rhs(c,t)

      do m = mode1, mode2
        do j = 1, ncompc

!         volume integral

          do i = 1, ndfc
            work6(i,j,m) = sum ( theta(:,i) * &
                      ( - ungradcn(:,j,m) + fng(:,j,m) ) * detF * wg )
          end do

!         boundary integral

          do side = 1, nsides
            do i = 1, ndfc
               work6(i,j,m) = work6(i,j,m) &
               - sum ( ugsn(:,side) * thetas(:,i,side) * cnjump(:,j,m,side) &
                           * curvels(:,side) * wgb, mask = ugsn(:,side) < 0 )
            end do
          end do

        end do
      end do

      elemvec = reshape ( work6, [ ndfc * ncompc * ( mode2 - mode1 + 1 ) ] )

    end if

    if ( last ) then

!     last element in this group

      call delete ( vemodel )

      deallocate ( wg, detF )
      deallocate ( xig, phi, x )
      deallocate ( dphi, F )
      deallocate ( dphidx )
      deallocate ( Finv )
      deallocate ( theta, dtheta )
      deallocate ( thetas )
      deallocate ( dthetadx )
      deallocate ( ungradtheta )
      deallocate ( ungradcn )
      deallocate ( u, tmp, gradu )
      deallocate ( xg )
      deallocate ( uvecn, cn )
      deallocate ( cng, gvecn )
      deallocate ( fng )
      deallocate ( xigb, wgb )
      deallocate ( phib, dphib )
      deallocate ( xigs, phis )
      deallocate ( ugsn, xgs )
      deallocate ( dxdxi, curvels )
      deallocate ( normals, us )
      deallocate ( xs )
      deallocate ( cnjump )
      deallocate ( cnside )
      deallocate ( work6 )

    end if

  end subroutine ce_dg_elem1


! set global parameters viscoelastic (internal element)

  subroutine set_globals_viscoelastic_DG ( coefficients )

    use viscoelastic_globals_m

    type(coefficients_t), intent(in) :: coefficients

    call check ( coefficients, 'set_globals_viscoelastic_DG', &
      indexarray=[24], minimum=[1] )

    nintbc = coefficients%i(24)

!   subdomain integration: not available

    call check ( coefficients, 'set_globals_viscoelastic_DG', &
      indexarray=[32,33,77], maximum=[1,1,0] )

!   interpolation: standard shape functions only using simple interface

    call check ( coefficients, 'set_globals_viscoelastic_DG', &
      indexarray=[1,2,12], maximum=[8,8,8] )

  end subroutine set_globals_viscoelastic_DG

end module viscoelastic_elements_DG_2D_m

