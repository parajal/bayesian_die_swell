
! Copyright (C) 2015-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the generalized Newtonian fluid flow equation
!
!    - div( eta(gammadot_e) (nabla u+nabla u^T) ) + nabla p = f
!      div u = 0
!
! with viscosity function eta depending on the (effective) shear rate
!
!   gammadot_e = sqrt(2D:D), with D=(nabla u+nabla u^T)/2
!

module generalized_newtonian_elements_m

  use tfem_elem_m
  use stokes_elements_m
  use stokes_elements_generic_m, only: set_stokes_shape_function_global

  implicit none


contains


! Internal element routine for generalized Newtonian fluid flow

  subroutine generalized_newtonian_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use generalized_stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, M, i, j, ip, vfuncnr
    real(dp) :: S_NMii, S_NMij, L_NMi


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

    if ( first) &
          call set_globals_generalized_newtonian ( mesh, coefficients, elgrp )

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   get velocity vector

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ndim] )

    ugvector =  matmul ( phi, uvector )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

!   evaluate fields for coefficients of the generalized Newtonian model

    call evaluate_fields ( mesh, problem, coefficients, oldvectors, elgrp, &
      elem, eta_const=etaconst, eta_0=eta0, eta_inf=etainf, lambda=lambda, &
      tau_y=tauyield, x=xg, phi=phi )

!   determine coefficient eta by a generalized Newtonian model

    call compute_viscosity_gn ( coefficients, coorsys, np=ninti, xp=xg, &
      up=ugvector, gup=gradu, eta_const=etaconst, eta_0=eta0, eta_inf=etainf, &
      lambda=lambda, tau_y=tauyield, eta=etag, deta=detag )

    if ( vector ) then

!     fill vector

      vfuncnr = coefficients%i(14)

      if ( vfuncnr > 0 ) then

        do ip = 1, ninti
          fg(ip,:) = coefficients%vfunc ( ndim, vfuncnr, xg(ip,:) )
        end do

        do j = 1, ndim
          do N = 1, ndf
            elemvec( pos(N,j) ) = sum ( fg(:,j) * phi(:,N) * detF * wg )
          end do
        end do

        elemvec( posp ) = 0

      else

        elemvec = 0

      end if

    end if


    if ( matrix ) then

!     fill matrix

!     diagonal blocks

      do N = 1, ndf
        do M = N, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dphidx(ip,M,:) )
          end do
          do i = 1, ndim
            work1 = work + dphidx(:,N,i) * dphidx(:,M,i)
            if ( coorsys == 1 .and. i == 2 ) then ! axisymmetric
              work1 = work1 + 2 * phi(:,N) * phi(:,M) / xg(:,2) ** 2
            end if
            S_NMii = sum ( etag * work1 * detF * wg )
            elemmat( pos(N,i), pos(M,i) ) = S_NMii
            elemmat( pos(M,i), pos(N,i) ) = S_NMii  ! = S_MNii symmetry
          end do
        end do
      end do

!     off-diagonal blocks

      do N = 1, ndf
        do M = 1, ndf
          do i = 1, ndim
            do j = i+1, ndim
              work = dphidx(:,N,j) * dphidx(:,M,i)
              S_NMij = sum ( etag * work * detF * wg )
              elemmat( pos(N,i), pos(M,j) ) = S_NMij
              elemmat( pos(M,j), pos(N,i) ) = S_NMij  ! = S_MNji symmetry
            end do
          end do
        end do
      end do

!     velocity-pressure part

      do N = 1, ndfp
        do M = 1, ndf
          do i = 1, ndim
            if ( coorsys == 1 .and. i == 2 ) then ! axisymmetric
              work = dphidx(:,M,2) + phi(:,M) / xg(:,2)
              L_NMi = sum ( psi(:,N) * work * detF * wg )
            else
              L_NMi = sum ( psi(:,N) * dphidx(:,M,i) * detF * wg )
            end if
            elemmat( posp(N), pos(M,i) ) = - L_NMi
            elemmat( pos(M,i), posp(N) ) = - L_NMi  ! = - L_MNi symmetry
          end do
        end do
      end do

!     pressure-pressure part

      elemmat( posp, posp ) = 0

!     set continuity equation (div u=0) to zero

      if ( coefficients%i(39) == 1 ) then

        do i = 1, ndim
          elemmat( posp, pos(:,i) ) = 0
        end do

      end if

    end if


!   Newton-Raphson

    if ( coefficients%i(253) == 1 ) then

!     (nabla v)^T:(nabla u+(nabla u)^T)

      do ip = 1, ninti
        work6(ip,:,:) = matmul ( dphidx(ip,:,:), &
                                 gradu(ip,:,:) + transpose(gradu(ip,:,:)) )
      end do
      if ( coorsys == 1 ) then
        ugvector(:,2) =  matmul ( phi, uvector(:,2) )
        do N = 1, ndf
          work6(:,N,2) = work6(:,N,2) &
                               + 2 * phi(:,N) * ugvector(:,2) / xg(:,2) ** 2
        end do
      end if

      if ( vector ) then

!       nabla.u

        work = 0
        do i = 1, ndim
          work = work + gradu(:,i,i)
        end do
        if ( coorsys == 1 ) then
          work = work + ugvector(:,2) / xg(:,2)
        end if

!       q nabla.u

        elemvec(posp) = elemvec(posp) + matmul ( work * detF * wg, psi )

!       - eta*(nabla v)^T:(nabla u+(nabla u)^T)

        do j = 1, ndim
          work2(:,j) = - matmul ( etag * detF * wg, work6(:,:,j) )
        end do
        elemvec(:ndf*ndim) = elemvec(:ndf*ndim) &
                                      + reshape ( work2, [ ndf*ndim ] )

      end if

      if ( matrix ) then

!       eta'/gammadot_e *
!       (nabla v)^T:(nabla u+(nabla u)^T)(nabla u+(nabla u)^T):(nabla du)

!       diagonal

        do N = 1, ndf
          do M = N, ndf
            do i = 1, ndim
              S_NMij = sum ( detag * work6(:,N,i) * work6(:,M,i) * detF * wg )
              elemmat( pos(N,i), pos(M,i) ) = elemmat( pos(N,i), pos(M,i) ) &
                                              + S_NMij
              elemmat( pos(M,i), pos(N,i) ) = elemmat( pos(N,i), pos(M,i) )
            end do
          end do
        end do

!       off-diagonal

        do N = 1, ndf
          do M = 1, ndf
            do i = 1, ndim
              do j = i+1, ndim
                S_NMij = sum ( detag * work6(:,N,i) * work6(:,M,j) * detF * wg )
                elemmat( pos(N,i), pos(M,j) ) = elemmat( pos(N,i), pos(M,j) ) &
                                                + S_NMij
                elemmat( pos(M,j), pos(N,i) ) = elemmat( pos(N,i), pos(M,j) )
              end do
            end do
          end do
        end do

      end if

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( etaconst(ninti), tauyield(ninti) )
      allocate ( etag(ninti), eta0(ninti), etainf(ninti), lambda(ninti) )
      allocate ( work(ninti), work1(ninti) )
      allocate ( u(ndim*ndf), uvector(ndf,ndim), ugvector(ninti,ndim) )
      allocate ( gradu(ninti,ndim,ndim) )

      if ( coefficients%i(253) == 1 ) then
        allocate ( work6(ninti,ndf,ndim), detag(ninti), work2(ndf,ndim) )
      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( etaconst, tauyield )
      deallocate ( etag, eta0, etainf, lambda )
      deallocate ( work, work1 )
      deallocate ( u, uvector, gradu, ugvector )

      if ( coefficients%i(253) == 1 ) then
        deallocate ( work6, detag, work2 )
      end if

    end subroutine deallocate_arrays

  end subroutine generalized_newtonian_elem


! Internal element routine to fill the position dependent coefficient eta
! in the Gauss/nodal points according to a generalized Newtonian model.
! The result is stored in a vector defined per element (elvector).
! This element should be used together with the routine loop_over_elements.
! The resulting elvector can be used as input to generalized Stokes elements.

  subroutine fill_eta_gn ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use generalized_stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors


    integer :: j

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

    if ( first) &
          call set_globals_generalized_newtonian ( mesh, coefficients, elgrp )

!   allocate more arrays

    if ( first .and. coefficients%i(255) == 1 ) then

!     first element in this group and nodal points

      nump = nodalp

      call deallocate_arrays_nodalp

      call allocate_arrays_nodalp

      call refcoor_nodal_points ( mesh, elgrp, xig ) ! use xig for xrnod

      call set_shape_function ( shapefunc, xig, phi, dphi )

    else if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      nump = ninti

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      nump = ninti

      call allocate_arrays

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call isoparametric_coordinates ( x, phi, xg )

!   get velocity vector

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ndim] )

    ugvector =  matmul ( phi, uvector )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

!   evaluate fields for coefficients of the generalized Newtonian model

    call evaluate_fields ( mesh, problem, coefficients, oldvectors, elgrp, &
      elem, eta_const=etaconst, eta_0=eta0, eta_inf=etainf, lambda=lambda, &
      tau_y=tauyield, x=xg, phi=phi )

!   determine viscosity coefficient eta by a generalized Newtonian model

    call compute_viscosity_gn ( coefficients, coorsys, np=nump, xp=xg, &
      up=ugvector, gup=gradu, eta_const=etaconst, eta_0=eta0, eta_inf=etainf, &
      lambda=lambda, tau_y=tauyield, eta=etag )

    call put_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, r1=etag )

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(255) == 1 ) then

!     last element in this group and nodal points

      call deallocate_arrays

    else if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays_nodalp

      allocate ( wg(nodalp), xig(nodalp,ndim), phi(nodalp,ndf), detF(nodalp) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( xg(nodalp,ndim) )

      allocate ( u(ndim*ndf), uvector(ndf,ndim), gradu(nodalp,ndim,ndim) )
      allocate ( etag(nodalp), ugvector(nodalp,ndim) )
      allocate ( etaconst(nodalp), tauyield(nodalp) )
      allocate ( eta0(nodalp), etainf(nodalp), lambda(nodalp) )

    end subroutine allocate_arrays_nodalp

    subroutine deallocate_arrays_nodalp

      deallocate ( wg, xig, phi, detF )
      deallocate ( dphi, F, Finv, dphidx )
      deallocate ( xg )

    end subroutine deallocate_arrays_nodalp

    subroutine allocate_arrays

      allocate ( u(ndim*ndf), uvector(ndf,ndim), gradu(ninti,ndim,ndim) )
      allocate ( etag(ninti), ugvector(ninti,ndim) )
      allocate ( etaconst(ninti), tauyield(ninti) )
      allocate ( eta0(ninti), etainf(ninti), lambda(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( u, uvector, gradu )
      deallocate ( etaconst, tauyield )
      deallocate ( etag, ugvector, eta0, etainf, lambda )

    end subroutine deallocate_arrays

  end subroutine fill_eta_gn


! function routine for the viscosity function

  function viscosity_function ( gammadot, coefficients, eta_const, eta_0, &
    eta_inf, lambda, tau_y )

    real(dp), intent(in) :: gammadot
    type(coefficients_t), intent(in) :: coefficients
    real(dp), intent(in) :: eta_const, eta_0, eta_inf, lambda, tau_y

    real(dp) :: viscosity_function

    real(dp) :: n, a, m, r, eta_max1, eta_max2
    integer :: gnmodel, yieldmodel

!   viscosity

    gnmodel = coefficients%i(254)

    select case( gnmodel )
    case(0)
!     Constant
      viscosity_function = eta_const
    case(1)
!     Power-law
      if ( gammadot > tiny(1._dp) ) then
        m = coefficients%r(201)
        n = coefficients%r(202)
        viscosity_function = m / gammadot ** (1-n)
      else
        eta_max1 = coefficients%r(211)
        viscosity_function = eta_max1
      end if
    case(2)
!     Carreau
      n       = coefficients%r(202)
      viscosity_function = eta_inf + &
         ( eta_0 - eta_inf ) / ( 1 + (lambda*gammadot)**2 ) ** ( (1-n)/2 )
    case(3)
!     Carreau-Yasuda
      n       = coefficients%r(202)
      a       = coefficients%r(206)
      viscosity_function = eta_inf + &
         ( eta_0 - eta_inf ) / ( 1 + (lambda*gammadot)**a ) ** ( (1-n)/a )
    case default
      write(*,'(/a,i0/)') &
        'Error viscosity_function: wrong gn model number: ', gnmodel
      stop
    end select

!   yield

    yieldmodel = coefficients%i(265)

    if ( yieldmodel > 0 ) then

      select case( yieldmodel )
      case(1)
!       Bingham
        if ( gammadot > tiny(1._dp) ) then
          viscosity_function = viscosity_function + tau_y / gammadot
        else
          eta_max2 = coefficients%r(212)
          viscosity_function = viscosity_function + eta_max2
        end if
      case(2)
!       Regularized Bingham (Papanastasiou)
        r = coefficients%r(210)
        if ( gammadot > tiny(1._dp) ) then
          viscosity_function = viscosity_function + &
                        tau_y * ( 1 - exp(-r*gammadot) ) / gammadot
        else
          viscosity_function = viscosity_function + tau_y * r
        end if
      case default
        write(*,'(/a,i0/)') &
          'Error viscosity_function: wrong yield model number: ', yieldmodel
        stop
      end select

    end if

  end function viscosity_function


! function routine for the derivative of the viscosity function divided by
! gammadot: eta'/gammadot

  function viscosity_derivative_function ( gammadot, coefficients, eta_0, &
    eta_inf, lambda, tau_y )

    real(dp), intent(in) :: gammadot
    type(coefficients_t), intent(in) :: coefficients
    real(dp), intent(in) :: eta_0, eta_inf, lambda, tau_y

    real(dp) :: viscosity_derivative_function

    real(dp) :: n, a, m, r
    integer :: gnmodel, yieldmodel

    gnmodel = coefficients%i(254)

    select case( gnmodel )
    case(0)
!     Constant
      viscosity_derivative_function = 0
    case(1)
!     Power-law
      m = coefficients%r(201)
      n = coefficients%r(202)
      viscosity_derivative_function = m * (n-1) / gammadot ** (3-n)
    case(2)
!     Carreau
      n       = coefficients%r(202)
      viscosity_derivative_function = &
         ( eta_0 - eta_inf ) * (n-1) * lambda**2 &
                            / ( 1 + (lambda*gammadot)**2 ) ** ( (3-n)/2 )
    case(3)
!     Carreau-Yasuda
      n       = coefficients%r(202)
      a       = coefficients%r(206)
      viscosity_derivative_function = &
         ( eta_0 - eta_inf ) * (n-1) * lambda**a * gammadot**(a-2) &
                            / ( 1 + (lambda*gammadot)**a ) ** ( (1+a-n)/a )
    case default
      write(*,'(/a,i0/)') &
        'Error viscosity_derivative_function: wrong model number: ', gnmodel
      stop
    end select

!   yield

    yieldmodel = coefficients%i(265)

    if ( yieldmodel > 0 ) then

      select case( yieldmodel )
      case(1)
!       Bingham
        viscosity_derivative_function = &
                   viscosity_derivative_function - tau_y / gammadot ** 3
      case(2)
!       Regularized Bingham (Papanastasiou)
        r = coefficients%r(210)
        viscosity_derivative_function = viscosity_derivative_function + &
            tau_y * ( (1+r*gammadot) * exp(-r*gammadot) - 1 ) / gammadot ** 3
      case default
        write(*,'(/a,i0/)') &
          'Error viscosity_function: wrong yield model number: ', yieldmodel
        stop
      end select

    end if

  end function viscosity_derivative_function


! Routine for computing eta and deta in the required points

  subroutine compute_viscosity_gn ( coefficients, coorsys, np, xp, up, gup, &
    eta_const, eta_0, eta_inf, lambda, tau_y, eta, deta, tau )

    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: coorsys, np
    real(dp), dimension(:,:), intent(in) :: up, xp
    real(dp), dimension(:,:,:), intent(in) :: gup

    real(dp), dimension(:), intent(in) :: eta_const, eta_0, eta_inf, lambda, &
      tau_y

    real(dp), dimension(:), intent(out) :: eta

!   compute derivative if present and coefficients%i(253)==1, otherwise not.
    real(dp), dimension(:), intent(out), optional :: deta

!   compute tau = eta * gammadot if present
    real(dp), dimension(:), intent(out), optional :: tau

    integer :: ip
    real(dp) :: gammadot

    do ip = 1, np

!     second invariant

      gammadot = sum ( ( gup(ip,:,:) + transpose(gup(ip,:,:)) ) ** 2 ) / 2

      if ( coorsys == 1 ) then
!       axisymmetric
        if ( xp(ip,2) < 1e-10_dp ) then
!         use du_r / dr (r=0) for gradu_tt
          gammadot = gammadot + 2 * ( gup(ip,2,2) ) ** 2
        else
!         use u_r / r for gradu_tt
          gammadot = gammadot + 2 * ( up(ip,2) / xp(ip,2) ) ** 2
        end if
      end if

      gammadot = sqrt ( gammadot )

      eta(ip) = viscosity_function ( gammadot, coefficients, eta_const(ip), &
        eta_0(ip), eta_inf(ip), lambda(ip), tau_y(ip) )

      if ( present(deta) .and. coefficients%i(253) == 1 ) then
!        Newton-Raphson, determine factor eta'/gammadot_e
         deta(ip) = viscosity_derivative_function ( gammadot, coefficients, &
           eta_0(ip), eta_inf(ip), lambda(ip), tau_y(ip) )
      end if

      if ( present(tau) ) then
!       shear stress = von Mises stress
        tau(ip) = eta(ip) * gammadot
      end if

    end do

  end subroutine compute_viscosity_gn


! compute viscosity in all nodes

  subroutine generalized_newtonian_viscosity ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use generalized_stokes_globals_m

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
      call set_globals_generalized_newtonian ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

!   get velocity vector

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ndim] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

!   evaluate fields for coefficients in the nodes

    call evaluate_fields ( mesh, problem, coefficients, oldvectors, elgrp, &
      elem, eta_const=etaconst, eta_0=eta0, eta_inf=etainf, lambda=lambda, &
      tau_y=tauyield, x=x, phi=phi )

!   determine viscosity coefficient eta in the nodes

    call compute_viscosity_gn ( coefficients, coorsys, np=nodalp, xp=x, &
      up=uvector, gup=gradu, eta_const=etaconst, eta_0=eta0, eta_inf=etainf, &
      lambda=lambda, tau_y=tauyield, eta=etag )

    elemvec = etag

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( dphi(nodalp,ndf,ndim) )
      allocate ( etaconst(nodalp), tauyield(nodalp) )
      allocate ( etag(nodalp), eta0(nodalp), etainf(nodalp), lambda(nodalp) )
      allocate ( u(ndim*ndf), uvector(ndf,ndim) )
      allocate ( gradu(nodalp,ndim,ndim) )
      allocate ( detF(nodalp), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xrnod, phi, x )
      deallocate ( dphi )
      deallocate ( etaconst, tauyield )
      deallocate ( etag, eta0, etainf, lambda )
      deallocate ( u, uvector )
      deallocate ( gradu )
      deallocate ( detF, F )
      deallocate ( Finv, dphidx )

    end subroutine deallocate_arrays

  end subroutine generalized_newtonian_viscosity


! compute stresses in all nodes (from velocity gradients)

  subroutine generalized_newtonian_stress ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    use generalized_stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: comp, j, maxstress


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_generalized_newtonian ( mesh, coefficients, elgrp )

!     check component = coefficients%i(13)

      select case ( coorsys )
        case(0,1)
          maxstress = 5
        case(2)
          maxstress = 7
        case default
          call errormsg_case_default ( 'generalized_newtonian_stress', &
            'coorsys', int_value=coorsys )
      end select

      call check ( coefficients, 'generalized_newtonian_stress', &
        indexarray=[13], minimum=[1], maximum=[maxstress] )

!     allocate arrays

      call allocate_arrays

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ndim] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

!   evaluate fields for coefficients in the nodes

    call evaluate_fields ( mesh, problem, coefficients, oldvectors, elgrp, &
      elem, eta_const=etaconst, eta_0=eta0, eta_inf=etainf, lambda=lambda, &
      tau_y=tauyield, x=x, phi=phi )

!   determine viscosity coefficient eta in the nodes

    call compute_viscosity_gn ( coefficients, coorsys, np=nodalp, xp=x, &
      up=uvector, gup=gradu, eta_const=etaconst, eta_0=eta0, eta_inf=etainf, &
      lambda=lambda, tau_y=tauyield, eta=etag, tau=tauvm )

    comp = coefficients%i(13)

    if ( coorsys <= 1 ) then

      select case(comp)
      case(1)
        elemvec = 2 * etag * gradu(:,1,1)                ! 2 eta dudx
      case(2)
        elemvec = etag * ( gradu(:,1,2) + gradu(:,2,1) ) ! eta ( dudy + dvdx )
      case(3)
        elemvec = 2 * etag * gradu(:,2,2)                ! 2 eta dvdy
      case(4)
        where ( x(:,2) < 1e-10_dp )
          elemvec = 2 * etag * gradu(:,2,2)              ! 2 eta du_r / dr (r=0)
        elsewhere
          elemvec = 2 * etag * uvector(:,2) / x(:,2)     ! 2 eta u_r / r
        end where
      case(5)
        elemvec = tauvm                                  ! eta * gammadot
      case default
        call errormsg_case_default ( 'generalized_newtonian_stress', &
          'comp', int_value=comp )
      end select

    else if ( coorsys == 2 ) then

      select case(comp)
      case(1)
        elemvec = 2 * etag * gradu(:,1,1)                ! 2 eta dudx
      case(2)
        elemvec = etag * ( gradu(:,1,2) + gradu(:,2,1) ) ! eta ( dudy + dvdx )
      case(3)
        elemvec = etag * ( gradu(:,1,3) + gradu(:,3,1) ) ! eta ( dudz + dwdx )
      case(4)
        elemvec = 2 * etag * gradu(:,2,2)                ! 2 eta dvdy
      case(5)
        elemvec = etag * ( gradu(:,2,3) + gradu(:,3,2) ) ! eta ( dvdz + dwdy )
      case(6)
        elemvec = 2 * etag * gradu(:,3,3)                ! 2 eta dwdz
      case(7)
        elemvec = tauvm                                  ! eta * gammadot
      case default
        call errormsg_case_default ( 'generalized_newtonian_stress', &
          'comp', int_value=comp )
      end select

    end if

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ndim*ndf), uvector(ndf,ndim), gradu(nodalp,ndim,ndim) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( etaconst(nodalp), tauyield(nodalp), tauvm(nodalp) )
      allocate ( etag(nodalp), eta0(nodalp), etainf(nodalp), lambda(nodalp) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( detF )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( etaconst, tauyield, tauvm )
      deallocate ( etag, eta0, etainf, lambda )

    end subroutine deallocate_arrays

  end subroutine generalized_newtonian_stress


! generalized Newtonian stress tensor

  subroutine generalized_newtonian_stress_tensor ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use generalized_stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer, save :: maxstress
    integer :: j


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_generalized_newtonian ( mesh, coefficients, elgrp )

!     test size of elemvec

      select case ( coorsys )
        case(0,1)
          maxstress = 3 + coorsys
        case(2)
          maxstress = 6
      case default
        call errormsg_case_default ( 'generalized_newtonian_stress_tensor', &
          'coorsys', int_value=coorsys )
      end select

      if ( size(elemvec) /= nodalp*maxstress ) then
        write(*,'(/2(a/))') 'Error in generalized_newtonian_stress_tensor:', &
          ' element vector has incorrect size for a tensor'
        stop
      end if

!     allocate arrays

      call allocate_arrays

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ndim] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

!   evaluate fields for coefficients in the nodes

    call evaluate_fields ( mesh, problem, coefficients, oldvectors, elgrp, &
      elem, eta_const=etaconst, eta_0=eta0, eta_inf=etainf, lambda=lambda, &
      tau_y=tauyield, x=x, phi=phi )

!   determine viscosity coefficient eta in the nodes

    call compute_viscosity_gn ( coefficients, coorsys, np=nodalp, xp=x, &
      up=uvector, gup=gradu, eta_const=etaconst, eta_0=eta0, eta_inf=etainf, &
      lambda=lambda, tau_y=tauyield, eta=etag )

    if ( coorsys <= 1 ) then

      work2(:,1) = 2 * etag * gradu(:,1,1)                ! tau_xx
      work2(:,2) = etag * ( gradu(:,1,2) + gradu(:,2,1) ) ! tau_xy
      work2(:,3) = 2 * etag * gradu(:,2,2)                ! tau_yy

      if ( coorsys == 1 ) then
!       axisymmetric
        where ( x(:,2) < 1e-10_dp )
          work2(:,4) = 2 * etag * gradu(:,2,2)           ! 2 eta du_r / dr (r=0)
        elsewhere
          work2(:,4) = 2 * etag * uvector(:,2) / x(:,2)  ! 2 eta u_r / r
        end where
      end if

    else if ( coorsys == 2 ) then

      work2(:,1) = 2 * etag * gradu(:,1,1)                ! tau_xx
      work2(:,2) = etag * ( gradu(:,1,2) + gradu(:,2,1) ) ! tau_xy
      work2(:,3) = etag * ( gradu(:,1,3) + gradu(:,3,1) ) ! tau_xz
      work2(:,4) = 2 * etag * gradu(:,2,2)                ! tau_yy
      work2(:,5) = etag * ( gradu(:,2,3) + gradu(:,3,2) ) ! tau_yz
      work2(:,6) = 2 * etag * gradu(:,3,3)                ! tau_zz

    end if

    elemvec = reshape ( work2, [nodalp*maxstress] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ndim*ndf), uvector(ndf,ndim), gradu(nodalp,ndim,ndim) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( work2(nodalp,maxstress) )
      allocate ( etaconst(nodalp), tauyield(nodalp) )
      allocate ( etag(nodalp), eta0(nodalp), etainf(nodalp), lambda(nodalp) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( detF )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( work2 )
      deallocate ( etaconst, tauyield )
      deallocate ( etag, eta0, etainf, lambda )

    end subroutine deallocate_arrays

  end subroutine generalized_newtonian_stress_tensor


! set global parameters (internal element)

  subroutine set_globals_generalized_newtonian ( mesh, coefficients, elgrp )

    use generalized_stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

!   check size of coefficients

    call check ( coefficients, 'set_globals_generalized_newtonian', &
      ncoefi=300, ncoefr=250, indexarray=[253,254], minimum=[0,0], &
      maximum=[1,3] )

  end subroutine set_globals_generalized_newtonian


! evaluate fields for coefficients of the generalized Newtonian models

  subroutine evaluate_fields ( mesh, problem, coefficients, oldvectors, elgrp, &
    elem, eta_const, eta_0, eta_inf, lambda, tau_y, x, phi )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    integer, intent(in) :: elgrp, elem
    real(dp), dimension(:), intent(out) :: eta_0, eta_inf, lambda, &
      eta_const, tau_y
    real(dp), dimension(:,:), intent(in) :: x, phi

!   eta_0

    call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(257), value=coefficients%r(203), &
      func=coefficients%func1(2)%p, funcnr=coefficients%i(258), x=x, &
      indx_v=2, layer=coefficients%i(38), phi=phi, indx_e=3, coef=eta_0 )

!   eta_inf

    call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(259), value=coefficients%r(204), &
      func=coefficients%func1(3)%p, funcnr=coefficients%i(260), x=x, &
      indx_v=3, layer=coefficients%i(38), phi=phi, indx_e=4, coef=eta_inf )

!   lambda

    call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(261), value=coefficients%r(205), &
      func=coefficients%func1(4)%p, funcnr=coefficients%i(262), x=x, &
      indx_v=4, layer=coefficients%i(38), phi=phi, indx_e=5, coef=lambda )

!   eta_const

    call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(263), value=coefficients%r(207), &
      func=coefficients%func1(5)%p, funcnr=coefficients%i(264), x=x, &
      indx_v=5, layer=coefficients%i(38), phi=phi, indx_e=6, coef=eta_const )

!   tau_y field

    call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(266), value=coefficients%r(209), &
      func=coefficients%func1(6)%p, funcnr=coefficients%i(267), x=x, &
      indx_v=6, layer=coefficients%i(38), phi=phi, indx_e=7, coef=tau_y )

  end subroutine evaluate_fields

end module generalized_newtonian_elements_m

