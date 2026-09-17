
! Copyright (C) 2006-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for stochastic models using BCF/DG
!

module bcf_elements_DG_2D_m

  use tfem_elem_m
  use stokes_set_globals_m
  use shapefunc_gauss_m
  use stochastic_models_2D_m

  implicit none


contains


! Internal element routine for the convection term in BCF/DG
! Explicit 1st order time-integration
! Note, that loop_over_elements must be used.

  subroutine bcf_conv_dg_elem ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use viscoelastic_globals_m
    use bcf_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors


    integer :: i, j, k, ip, side, elgrpnr, sidenr, elemnr
    integer :: elshape


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_stochastic_model ( coefficients )
      call set_globals_bcf ( mesh, coefficients, elgrp )
      call set_globals_bcf_DG ( mesh, coefficients, elgrp )

      nfield = stnumpar%nfield

!     allocate arrays

      allocate ( wg(ninti), detF(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( dphidx(ninti,ndf,ndim) )
      allocate ( Finv(ninti,ndim,ndim) )
      allocate ( theta(ninti,ndfq), dtheta(ninti,ndfq,ndim) )
      allocate ( dthetadx(ninti,ndfq,ndim) )
      allocate ( ungradtheta(ninti,ndfq) )
      allocate ( ungradqn(nfield,ninti,ncompq) )
      allocate ( u(ndim*ndf), tmp(ndf,ndim) )
      allocate ( xg(ninti,ndim) )

      allocate ( posq(ndfq*ncompq), posqs(ndfq*ncompq) )
      allocate ( uvecn(ninti,ndim), qn(nfield,ndfq,ncompq) )
      allocate ( qv(ncompq*ndfq) )

      allocate ( xigb(nintbc), wgb(nintbc) )
      allocate ( phib(nintbc,ndfb), dphib(nintbc,ndfb) )
      allocate ( xigs(nintbc,ndim,nsides), phis(nintbc,ndf,nsides) )
      allocate ( thetas(nintbc,ndfq,nsides) )
      allocate ( ugsn(nintbc,nsides), xgs(nintbc,ndim,nsides) )
      allocate ( dxdxi(nintbc,ndim), curvels(nintbc,nsides) )
      allocate ( normals(nintbc,ndim,nsides), us(nintbc,ndim) )
      allocate ( xs(nodalpb,ndim,nsides) )
      allocate ( qnjump(nfield,nintbc,ncompq,nsides) )
      allocate ( qnside(nfield,ndfq,ncompq) )
      allocate ( work6(nfield,ndfq,ncompq), work(ninti), work1(nintbc) )
      allocate ( work8(nintbc) )

      work1 = 0
      work8 = 0

      elshape = mesh%element(elgrp)%elshape

      if ( mesh%element(elgrp)%globalshape == 'quadrilateral' .and. &
           elshape /= 6 ) then
        write(*,'(/a/a/)') 'Error in bcf_conv_dg_elem:', &
          ' For quadrilaterals only elshape=6 is allowed '
        stop
      else if ( mesh%element(elgrp)%globalshape == 'triangle' .and. &
           elshape /= 4 .and. elshape /= 7 ) then
        write(*,'(/a/a/)') 'Error in bcf_conv_dg_elem:', &
          ' For triangles only elshape=4 and 7 are allowed '
        stop
      end if

!     set Gauss integration for internal element and boundary integral

      call set_Gauss_integration ( globalshape, xig, wg )

      call Gauss_Legendre_line ( nintbc, xigb, wgb )

!     set shape functions

      call set_shape_function ( globalshape, xig, phi, dphi ) ! velocity
      call set_shape_function_gauss ( globalshape, xig, theta, dtheta ) ! Q

!     isoparametric shape of boundary
      call shape_line_P2 ( xigb, phib, dphib )

      call spread_to_sides_2D ( xigb, xigs )

      do side = 1, nsides
!       shape function on the side for velocity and c tensor
        call set_shape_function ( globalshape, xigs(:,:,side), phis(:,:,side) )
        call set_shape_function_gauss ( globalshape, xigs(:,:,side), &
          thetas(:,:,side))
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
      elgrp, elem, u, physq=[physqvel] )

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

!   get Q-vectors at previous time step

    call get_sysvector ( mesh, problem, oldvectors%s1(1)%p(1), elgrp, elem, &
      qv, posu=posq )

    qn(1,:,:) = reshape ( qv, [ndfq,ncompq] )

    do j = 2, nfield
      qn(j,:,:) = reshape ( oldvectors%s1(1)%p(j)%u(posq), [ndfq,ncompq] )
    end do

!   un.grad qn term

    do j = 1, nfield
      ungradqn(j,:,:) = matmul ( ungradtheta, qn(j,:,:) )
    end do

!   qnjump on the sides

    do side = 1, nsides

      elgrpnr = mesh%sidelem(elgrp)%a(side,elem,1)

      if ( elgrpnr /= 0 ) then

!       side element found: get qn in side element

        elemnr = mesh%sidelem(elgrp)%a(side,elem,2)
        sidenr = mesh%sidelem(elgrp)%a(side,elem,3)

        call get_sysvector ( mesh, problem, oldvectors%s1(1)%p(1), &
          elgrpnr, elemnr, qv, posu=posqs, order='DN' )

        qnside(1,:,:) = reshape ( qv, [ndfq,ncompq] )

        do j = 2, nfield
          qnside(j,:,:) = &
                    reshape ( oldvectors%s1(1)%p(j)%u(posqs), [ndfq,ncompq] )
        end do

      end if

!     qnjump = qn(outside) - qn(inside) for each intgr point on this side

      do ip = 1, nintbc

        if ( ugsn(ip,side) < 0 ) then

          do j = 1, nfield
            qnjump(j,ip,:,side) = - matmul ( thetas(ip,:,side), qn(j,:,:) )
          end do

          if ( elgrpnr /= 0 ) then
            do j = 1, nfield
              qnjump(j,ip,:,side) = qnjump(j,ip,:,side) &
                   + matmul ( thetas(nintbc-ip+1,:,sidenr), qnside(j,:,:) )
            end do
          else
            qnjump(:,ip,:,side) = qnjump(:,ip,:,side) &
               + coefficients%qinflow ( nfield, ncompq, xgs(ip,:,side), &
                                        ugsn(ip,side) )
          end if

        else

          qnjump(:,ip,:,side) = 0

        end if

      end do

    end do


!   right-hand side of \dot Q = - u.grad Q, i.e. convection term only

!   volume integral

    do i = 1, ndfq
      work = theta(:,i) * detF * wg
      do k = 1, ncompq
        do j = 1, nfield
          work6(j,i,k) = - sum ( work * ungradqn(j,:,k) )
        end do
      end do
    end do

!   boundary integral

    do side = 1, nsides
      where ( ugsn(:,side) < 0 )
        work8 = ugsn(:,side) * curvels(:,side) * wgb
      end where
      do i = 1, ndfq
        where ( ugsn(:,side) < 0 )
          work1 = work8 * thetas(:,i,side)
        end where
        do k = 1, ncompq
          do j = 1, nfield
            work6(j,i,k) = work6(j,i,k) &
              - sum ( work1 * qnjump(j,:,k,side), mask = ugsn(:,side) < 0 )
          end do
        end do
      end do
    end do

!   put result in oldvectors%s2

    do j = 1, nfield
      oldvectors%s1(2)%p(j)%u(posq) = reshape ( work6(j,:,:), [ndfq*ncompq] )
    end do


    if ( last ) then

!     last element in this group

      call delete ( stmodel )

      deallocate ( wg, detF )
      deallocate ( xig, phi, x )
      deallocate ( dphi, F )
      deallocate ( dphidx )
      deallocate ( Finv )
      deallocate ( theta, dtheta )
      deallocate ( dthetadx )
      deallocate ( ungradtheta )
      deallocate ( ungradqn )
      deallocate ( u, tmp )
      deallocate ( xg )

      deallocate ( posq, posqs )
      deallocate ( uvecn, qn )
      deallocate ( qv )

      deallocate ( xigb, wgb )
      deallocate ( phib, dphib )
      deallocate ( xigs, phis )
      deallocate ( thetas )
      deallocate ( ugsn, xgs )
      deallocate ( dxdxi, curvels )
      deallocate ( normals, us )
      deallocate ( xs )
      deallocate ( qnjump )
      deallocate ( qnside )
      deallocate ( work6, work, work1 )
      deallocate ( work8 )

    end if

  end subroutine bcf_conv_dg_elem


! Internal element routine for the model terms in BCF/DG.
! Explicit 1st order time-integration
! Note, that loop_over_elements must be used.

  subroutine bcf_model_dg_elem ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use viscoelastic_globals_m
    use bcf_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors


    integer :: i, j
    integer :: elshape


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_stochastic_model ( coefficients )
      call set_globals_bcf ( mesh, coefficients, elgrp )
      call set_globals_bcf_DG ( mesh, coefficients, elgrp )

      nfield = stnumpar%nfield

!     allocate arrays (number of Gauss points == ndfq now!)

      ninti = ndfq

      allocate ( wg(ninti), detF(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( dphidx(ninti,ndf,ndim) )
      allocate ( Finv(ninti,ndim,ndim) )
      allocate ( u(ndim*ndf), tmp(ndf,ndim), gradu(ninti,ndim,ndim)  )
      allocate ( xg(ninti,ndim) )

      allocate ( posq(ndfq*ncompq) )
      allocate ( uvecn(ninti,ndim), qn(nfield,ncompq,ndfq) )
      allocate ( qnp1(nfield,ncompq,ndfq) )
      allocate ( qv(ncompq*ndfq), gvecn(ninti,4+coorsys) )
      allocate ( brownf(nfield,ncompq) )

      allocate ( work(ninti) )

      elshape = mesh%element(elgrp)%elshape

      if ( mesh%element(elgrp)%globalshape == 'quadrilateral' .and. &
           elshape /= 6 ) then
        write(*,'(/a/a/)') 'Error in bcf_model_dg_elem:', &
          ' For quadrilaterals only elshape=6 is allowed '
        stop
      else if ( mesh%element(elgrp)%globalshape == 'triangle' .and. &
           elshape /= 4 .and. elshape /= 7 ) then
        write(*,'(/a/a/)') 'Error in bcf_model_dg_elem:', &
          ' For triangles only elshape=4 and 7 are allowed '
        stop
      end if

!     set Gauss integration for internal element

      call set_Gauss_integration ( globalshape, xig, wg )

!     set shape functions

      call set_shape_function ( globalshape, xig, phi, dphi ) ! velocity

!     set Brownian vector

      brownf = coefficients%brownian_vector ( nfield, ncompq )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )


!   get velocity vector

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel] )

    tmp = reshape ( u, [ndf,ndim] )

    uvecn = matmul ( phi, tmp )

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

!   get Q-vectors at previous time step

    call get_sysvector ( mesh, problem, oldvectors%s1(1)%p(1), elgrp, elem, &
      qv, posu=posq )

    qn(1,:,:) = transpose ( reshape ( qv, [ndfq,ncompq] ) )

    do j = 2, nfield
      qn(j,:,:) = &
        transpose ( reshape ( oldvectors%s1(1)%p(j)%u(posq), [ndfq,ncompq] ) )
    end do

!   get convection term for Q-vectors at previous time step

    do j = 1, nfield
      qnp1(j,:,:) = &
        transpose ( reshape ( oldvectors%s1(2)%p(j)%u(posq), [ndfq,ncompq] ) )
    end do

!   divide by diagonal mass matrix

    work = detF * wg

    do i = 1, ndfq
      qnp1(:,:,i) = qnp1(:,:,i) / work(i)
    end do


!   step forward \dot Q = ....

    call ststep_2D ( stmodel, stnumpar, gvecn, qn, qnp1, brownf, addvec=.true. )


!   put result in oldvectors%s2

    do j = 1, nfield
      oldvectors%s1(2)%p(j)%u(posq) = &
                  reshape ( transpose ( qnp1(j,:,:) ), [ndfq*ncompq] )
    end do


    if ( last ) then

!     last element in this group

      call delete ( stmodel )

      deallocate ( wg, detF )
      deallocate ( xig, phi, x )
      deallocate ( dphi, F )
      deallocate ( dphidx )
      deallocate ( Finv )
      deallocate ( u, tmp, gradu )
      deallocate ( xg )

      deallocate ( posq )
      deallocate ( uvecn, qn )
      deallocate ( qnp1 )
      deallocate ( qv, gvecn )
      deallocate ( brownf )

      deallocate ( work )

    end if

  end subroutine bcf_model_dg_elem


! Internal element routine for the right-hand side (div tau) of the Stokes
! problem

  subroutine bcf_rhs_divtau ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m
    use bcf_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: j, N, ip


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_stochastic_model ( coefficients )
      call set_globals_bcf ( mesh, coefficients, elgrp )
      call set_globals_bcf_DG ( mesh, coefficients, elgrp )

      nfield = stnumpar%nfield

      allocate ( wg(ninti), detF(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( theta(ninti,ndfq) )
      allocate ( work(ninti) )
      allocate ( xg(ninti,ndim) )

      allocate ( tauvecq(ndfq,ncompt) )
      allocate ( tauvec(ninti,ncompt), tauten(ninti,ndim,ndim) )
      allocate ( posq(ndfq*ncompq) )
      allocate ( qnp1(nfield,ncompq,ndfq) )  ! transpos
      allocate ( qv(ncompq*ndfq) )


!     set Gauss integration

      call set_Gauss_integration ( globalshape, xig, wg )

!     set shape functions

      call set_shape_function ( globalshape, xig, phi, dphi ) ! velocity
      call set_shape_function_gauss ( globalshape, xig, theta ) ! Q

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )


!   get Q-vectors

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s1(2)%p(1), &
      elgrp, elem, qv, posu=posq )

    qnp1(1,:,:) = transpose( reshape ( qv, [ndfq,ncompq] ) )

    do j = 2, nfield
      qnp1(j,:,:) = &
        transpose ( reshape ( oldvectors%s1(2)%p(j)%u(posq), [ndfq,ncompq] ) )
    end do


!   stress tensor computed by projection on Q-space

    if ( logq ) then
      write(*,'(/a/a/)') 'Error in bcf_rhs_divtau: log Q not yet available'
      stop
    else
      call stress_stochastic_2D ( stmodel, qnp1, tauvecq )
    end if

    tauvec = matmul ( theta, tauvecq )


!   build equations

    if ( matrix ) then

      write(*,'(/a/a/)') 'Error in bcf_rhs_divtau: no matrix to build.', &
        'Call build_system with buildmatrix=.false.'
      stop

    end if

    if ( vector ) then

!     - (nabla v)^T:tau

      tauten(:,1,1) = tauvec(:,1)
      tauten(:,1,2) = tauvec(:,2)
      tauten(:,2,1) = tauvec(:,2)
      tauten(:,2,2) = tauvec(:,3)

      do N = 1, ndf
        do ip = 1, ninti
          work(ip) = sum ( dphidx(ip,N,:) * tauten(ip,:,1) )
        end do
        elemvec(N) = - sum ( work * detF * wg )
        do ip = 1, ninti
          work(ip) = sum ( dphidx(ip,N,:) * tauten(ip,:,2) )
        end do
        if ( coorsys == 1 ) then
          work = work + tauvec(:,4) * phi(:,N) / xg(:,2)
        end if
        elemvec(ndf+N) = - sum ( work * detF * wg )
      end do

    end if

    if ( last ) then

!     last element in this group

      call delete ( stmodel )

      deallocate ( wg, detF )
      deallocate ( xig, phi, x )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( xg )
      deallocate ( theta )
      deallocate ( work )
      deallocate ( posq )
      deallocate ( qnp1, qv )
      deallocate ( tauvecq, tauvec, tauten )

    end if

  end subroutine bcf_rhs_divtau


! compute structure tensor in all nodes

  subroutine deriv_structure_tensor ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

!   choose transformation scheme (log or standard)

    select case ( coefficients%i(454) )
      case(0) ! standard scheme
        call deriv_structure_tensor_std ( mesh, problem, elgrp, elem, &
          first, last, coefficients, oldvectors, elemvec, elemwts )
      case(1) ! log Q
!        call deriv_structure_tensor_log ( mesh, problem, elgrp, elem, &
!          first, last, coefficients, oldvectors, elemvec, elemwts )
        write(*,'(/a/a,i0/)') 'Error derive_structure_tensor:', &
        ' log Q not yet available'
        stop
      case default
        write(*,'(/a/a,i0/)') 'Error derive_structure_tensor:', &
        ' incorrect value for log parameter coefficients%i(454) = ', &
        coefficients%i(454)
        stop
    end select

  end subroutine deriv_structure_tensor


! compute structure tensor in all nodes (standard)
! (OLD VERSION: computes structure tensor in Gauss points and extrapolates)

  subroutine deriv_structure_tensor_std_old ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m
    use bcf_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: i, j


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_stochastic_model ( coefficients )
      call set_globals_bcf ( mesh, coefficients, elgrp )
      call set_globals_bcf_DG ( mesh, coefficients, elgrp )

      nfield = stnumpar%nfield

      allocate ( theta(nodalp,ndfq) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( posq(ndfq*ncompq) )
      allocate ( qg(nfield,ncompq,ndfq) )
      allocate ( qv(ncompq*ndfq) )
      allocate ( c(ndfq,ncompc) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function_gauss ( globalshape, xrnod, theta )

    end if

!   get Q-vectors

    call get_sysvector ( mesh, problem, oldvectors%s1(1)%p(1), elgrp, elem, &
      qv, posu=posq )

    qg(1,:,:) = transpose ( reshape ( qv, [ndfq,ncompq] ) )

    do j = 2, nfield
      qg(j,:,:) = &
        transpose ( reshape ( oldvectors%s1(1)%p(j)%u(posq), [ndfq,ncompq] ) )
    end do

!   compute structure tensor in the Gauss points

    do i = 1, ndfq
      call structure_tensor_2D ( stmodel, qg(:,:,i), c(i,:) )
    end do

!   interpolate to the nodal points

    elemvec = reshape ( matmul ( theta, c ), [nodalp*ncompc] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( stmodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( posq )
      deallocate ( qg, qv, c )

    end if

  end subroutine deriv_structure_tensor_std_old


! compute structure tensor in all nodes (standard)

  subroutine deriv_structure_tensor_std ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m
    use bcf_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: i, j


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_stochastic_model ( coefficients )
      call set_globals_bcf ( mesh, coefficients, elgrp )
      call set_globals_bcf_DG ( mesh, coefficients, elgrp )

      nfield = stnumpar%nfield

      allocate ( theta(nodalp,ndfq) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( posq(ndfq*ncompq) )
      allocate ( qg(nfield,ncompq,ndfq) )
      allocate ( qn(nfield,ncompq,nodalp) )
      allocate ( qv(ncompq*ndfq) )
      allocate ( c(nodalp,ncompc) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function_gauss ( globalshape, xrnod, theta )

    end if

!   get Q-vectors

    call get_sysvector ( mesh, problem, oldvectors%s1(1)%p(1), elgrp, elem, &
      qv, posu=posq )

    qg(1,:,:) = transpose ( reshape ( qv, [ndfq,ncompq] ) )

    do j = 2, nfield
      qg(j,:,:) = &
        transpose ( reshape ( oldvectors%s1(1)%p(j)%u(posq), [ndfq,ncompq] ) )
    end do

!   compute Q vectors in the nodal points by interpolation

    do j = 1, nfield
      qn(j,:,:) = matmul ( qg(j,:,:), transpose(theta) )
    end do

!   compute structure tensor in the nodal points

    do i = 1, nodalp
      call structure_tensor_2D ( stmodel, qn(:,:,i), c(i,:) )
    end do

!   element vector

    elemvec = reshape ( c, [nodalp*ncompc] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( stmodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( posq )
      deallocate ( qg, qn, qv, c )

    end if

  end subroutine deriv_structure_tensor_std


! compute determinant of structure tensor in all nodes (standard)

  subroutine deriv_determinant_std ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m
    use bcf_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: i, j


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_stochastic_model ( coefficients )
      call set_globals_bcf ( mesh, coefficients, elgrp )
      call set_globals_bcf_DG ( mesh, coefficients, elgrp )

      nfield = stnumpar%nfield

      allocate ( theta(nodalp,ndfq) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( posq(ndfq*ncompq) )
      allocate ( qg(nfield,ncompq,ndfq) )
      allocate ( qn(nfield,ncompq,nodalp) )
      allocate ( qv(ncompq*ndfq) )
      allocate ( c(nodalp,ncompc) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function_gauss ( globalshape, xrnod, theta )

    end if

!   get Q-vectors

    call get_sysvector ( mesh, problem, oldvectors%s1(1)%p(1), elgrp, elem, &
      qv, posu=posq )

    qg(1,:,:) = transpose ( reshape ( qv, [ndfq,ncompq] ) )

    do j = 2, nfield
      qg(j,:,:) = &
        transpose ( reshape ( oldvectors%s1(1)%p(j)%u(posq), [ndfq,ncompq] ) )
    end do

!   compute Q vectors in the nodal points by interpolation

    do j = 1, nfield
      qn(j,:,:) = matmul ( qg(j,:,:), transpose(theta) )
    end do

!   compute structure tensor in the nodal points

    do i = 1, nodalp
      call structure_tensor_2D ( stmodel, qn(:,:,i), c(i,:) )
    end do

!   element vector determinant

    select case (ncompc)
    case(3)
      elemvec = c(:,1) * c(:,3) - c(:,2) ** 2
    case(4)
      elemvec = ( c(:,1) * c(:,3) - c(:,2) ** 2 ) * c(:,4)
    case default
      call errormsg_case_default ( 'deriv_determinant_std', &
        'ncompc', int_value=ncompc )
    end select

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( stmodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( posq )
      deallocate ( qg, qn, qv, c )

    end if


  end subroutine deriv_determinant_std


! compute determinant of structure tensor in the Gauss points (standard)
! This routine should be called with a vector that:
!  * has a number of degrees in an element equal to the number of
!    Gauss points (=ndfq)
!  * degrees must be discontinuous across elements by defining the vector as
!    .elementwise.=true. or storing all degrees in an internal node.

  subroutine deriv_determinant_gauss_std ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m
    use bcf_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: i, j


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_stochastic_model ( coefficients )
      call set_globals_bcf ( mesh, coefficients, elgrp )
      call set_globals_bcf_DG ( mesh, coefficients, elgrp )

      nfield = stnumpar%nfield

      allocate ( theta(nodalp,ndfq) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( posq(ndfq*ncompq) )
      allocate ( qg(nfield,ncompq,ndfq) )
      allocate ( qv(ncompq*ndfq) )
      allocate ( c(ndfq,ncompc) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function_gauss ( globalshape, xrnod, theta )

    end if

!   get Q-vectors

    call get_sysvector ( mesh, problem, oldvectors%s1(1)%p(1), elgrp, elem, &
      qv, posu=posq )

    qg(1,:,:) = transpose ( reshape ( qv, [ndfq,ncompq] ) )

    do j = 2, nfield
      qg(j,:,:) = &
        transpose ( reshape ( oldvectors%s1(1)%p(j)%u(posq), [ndfq,ncompq] ) )
    end do

!   compute structure tensor in the Gauss points

    do i = 1, ndfq
      call structure_tensor_2D ( stmodel, qg(:,:,i), c(i,:) )
    end do

!   element vector determinant

    select case (ncompc)
    case(3)
      elemvec = c(:,1) * c(:,3) - c(:,2) ** 2
    case(4)
      elemvec = ( c(:,1) * c(:,3) - c(:,2) ** 2 ) * c(:,4)
    case default
      call errormsg_case_default ( 'deriv_determinant_gauss_std', &
        'ncompc', int_value=ncompc )
    end select

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( stmodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( posq )
      deallocate ( qg, qv, c )

    end if


  end subroutine deriv_determinant_gauss_std


! compute stress tensor in all nodes (standard)
! This routine computes stress in the Gauss points first and interpolates the
! stress to the nodal points of the element for printing/plotting using the
! shape function of Q

  subroutine deriv_stress_tensor_std ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m
    use bcf_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: j


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_stochastic_model ( coefficients )
      call set_globals_bcf ( mesh, coefficients, elgrp )
      call set_globals_bcf_DG ( mesh, coefficients, elgrp )

      nfield = stnumpar%nfield

      allocate ( theta(nodalp,ndfq) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( posq(ndfq*ncompq) )
      allocate ( qn(nfield,ncompq,ndfq) )
      allocate ( qv(ncompq*ndfq) )
      allocate ( tauvec(ndfq,ncompt) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function_gauss ( globalshape, xrnod, theta )

    end if

!   get Q-vectors

    call get_sysvector ( mesh, problem, oldvectors%s1(1)%p(1), elgrp, elem, &
      qv, posu=posq )

    qn(1,:,:) = transpose ( reshape ( qv, [ndfq,ncompq] ) )

    do j = 2, nfield
      qn(j,:,:) = &
        transpose ( reshape ( oldvectors%s1(1)%p(j)%u(posq), [ndfq,ncompq] ) )
    end do

!   compute stress tensor in the Gauss points

    call stress_stochastic_2D ( stmodel, qn, tauvec )

!   interpolate to the nodal points

    elemvec = reshape ( matmul ( theta, tauvec ), [nodalp*ncompt] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( stmodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( posq )
      deallocate ( qn, qv, tauvec )

    end if

  end subroutine deriv_stress_tensor_std


! compute stress tensor in all nodes (standard)
! This routine computes the Q-vectors in the Gauss points first and
! interpolates the Q-vector to the nodal points of the element.
! The nodal point values of Q are used to compute the stress in the nodal
! points for printing/plotting.

  subroutine deriv_stress_tensor_std2 ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m
    use bcf_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: j


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_stochastic_model ( coefficients )
      call set_globals_bcf ( mesh, coefficients, elgrp )
      call set_globals_bcf_DG ( mesh, coefficients, elgrp )

      nfield = stnumpar%nfield

      allocate ( theta(nodalp,ndfq) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( posq(ndfq*ncompq) )
      allocate ( qg(nfield,ncompq,ndfq) )
      allocate ( qn(nfield,ncompq,nodalp) )
      allocate ( qv(ncompq*ndfq) )
      allocate ( tauvec(nodalp,ncompt) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function_gauss ( globalshape, xrnod, theta )

    end if

!   get Q-vectors

    call get_sysvector ( mesh, problem, oldvectors%s1(1)%p(1), elgrp, elem, &
      qv, posu=posq )

    qg(1,:,:) = transpose ( reshape ( qv, [ndfq,ncompq] ) )

    do j = 2, nfield
      qg(j,:,:) = &
        transpose ( reshape ( oldvectors%s1(1)%p(j)%u(posq), [ndfq,ncompq] ) )
    end do

!   compute Q vectors in the nodal points by interpolation

    do j = 1, nfield
      qn(j,:,:) = matmul ( qg(j,:,:), transpose(theta) )
    end do

!   compute stress tensor in the nodal points

    call stress_stochastic_2D ( stmodel, qn, tauvec )

!   element vector

    elemvec = reshape ( tauvec, [nodalp*ncompt] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( stmodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( posq )
      deallocate ( qg, qn, qv, tauvec )

    end if

  end subroutine deriv_stress_tensor_std2


! set stochastic model

  subroutine set_stochastic_model ( coefficients )

    use bcf_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   set viscoelastic model (don't forget to delete stmodel again)

    integer :: p, npar_mode, npar_tot, k

    nmodes = 1

!   check size of integer coefficients

    call check ( coefficients, 'set_stochastic_model', ncoefi=453 )

!   first create model structure stmodel

    call create_stochastic_model ( model=coefficients%i(452), stmodel=stmodel, &
      coorsys=coefficients%i(23) )

!   fill material parameters from the coefficients

    npar_mode = 2 + size(stmodel%nonlin,1)
    npar_tot  = nmodes * npar_mode

!   start pointer

    call check ( coefficients, 'set_stochastic_model', indexarray=[453], &
      minimum=[451] )

    p = coefficients%i(453)

!   check size of real array

    call check ( coefficients, 'set_stochastic_model', ncoefr=p+npar_tot-1 )

    stmodel%modulus = coefficients%r(p:p+npar_tot-1:npar_mode)
    stmodel%lambda = coefficients%r(p+1:p+npar_tot-1:npar_mode)
    do k = 1, size(stmodel%nonlin,1)
      stmodel%nonlin(k,:) = &
                    coefficients%r(p+1+k:p+npar_tot-1:npar_mode)
    end do

  end subroutine set_stochastic_model


! set global parameters bcf (internal element)

  subroutine set_globals_bcf ( mesh, coefficients, elgrp )

    use bcf_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

    ncompq = stmodel%ncompq
    ncompt = stmodel%ncompt
    ncompc = ncompt ! assume, for now, equal number of components tau_p and c

    physqgrad = coefficients%i(9)

!   check size of integer array

    call check ( coefficients, 'set_globals_bcf', ncoefi=500, ncoefr=450, &
      indexarray=[454,456], minimum=[0,1], maximum=[1,3] )

    if ( coefficients%i(454) == 0 ) then
      logq = .false.
    else if ( coefficients%i(454) == 1 ) then
      logq = .true.
    end if

    stnumpar%nfield = coefficients%i(455)
    stnumpar%timeint = coefficients%i(456)
    stnumpar%tstep = coefficients%r(401)

!   interpolation: standard shape functions only, using simple interface

    call check ( coefficients, 'set_globals_bcf', &
      indexarray=[1,2,451], maximum=[8,8,8] )

!   set number of degrees of freedom for Q

    intpolq = coefficients%i(451)

    shapefuncq%globalshape = globalshape
    shapefuncq%interpolation = intpolq
    shapefuncq%numbering = 'regular'

    call set_ndf ( shapefuncq, 'set_globals_bcf', ndf=ndfq )

!   subdomain integration: not available

    call check ( coefficients, 'set_globals_bcf', &
      indexarray=[32,33,77], maximum=[1,1,0] )

  end subroutine set_globals_bcf



! set global parameters viscoelastic (internal element)

  subroutine set_globals_bcf_DG ( mesh, coefficients, elgrp )

    use bcf_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

    call check ( coefficients, 'set_globals_bcf_DG', &
      indexarray=[457], minimum=[1] )

    nintbc = coefficients%i(457)

  end subroutine set_globals_bcf_DG

end module bcf_elements_DG_2D_m

