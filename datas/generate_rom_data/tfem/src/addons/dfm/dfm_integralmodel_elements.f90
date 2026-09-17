
! Copyright (C) 2016-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for integral models solved with deformation fields

module dfm_integralmodel_elements_m

  use dfm_field_elements_m

  implicit none

  save

contains


! compute stress tensor for a given deformation field history (no elements).
! Note: 3D only.

  subroutine dfm_stress_tensor_point ( coefficients, F_vec, stress_ten )

    use dfm_globals_m

!   coefficients defining the model
    type(coefficients_t), intent(in) :: coefficients

!   F(tau) tensor history: F_vec(ncompf=9,ninttau,nintvaltau)
    real(dp), dimension(:,:,:), intent(in) :: F_vec

!   The stress tensor of the integral model (3D)
    real(dp), dimension(3,3), intent(out) :: stress_ten


    integer :: j, intval, field, loca

!   set globals

!   some variables in elements set by set_globals_stokes_vp
    ndim = 3  ! 3D stress tensor
    globalshape = 'quadrilateral' ! fake shape

    call set_globals_dfm ( coefficients )
    call set_globals_intmodel ( coefficients )

!   allocate arrays

    allocate ( Fvec1(ncompf,ninttau,nintvaltau) )
    allocate ( tau_gp(ninttau,nintvaltau) )
    allocate ( tau_gw(ninttau,nintvaltau) )
    allocate ( mf(ninttau,nintvaltau) )
    allocate ( ivr(ninttau,nintvaltau,2) )
    allocate ( dpf(ninttau,nintvaltau,2) )
    allocate ( phitau(ninttau,ninttau) )
    allocate ( Ften(ndim,ndim,ninttau,nintvaltau) )
    allocate ( Bten(ndim,ndim,ninttau,nintvaltau) )

!   tau Gauss points

    call dfm_tau_gauss ( coefficients, tau_gp, tau_gw, phitau )

!   memory function in the tau Gauss points

    call memory_function ( coefficients, tau_gp, mf )

    loca = coefficients%i(361)

    if ( loca == 0 ) then

!     coupled fields

!     interpolate to tau Gauss nodes

      do intval = 1, nintvaltau
        do j = 1, ncompf
          Fvec1(j,:,intval) = matmul ( F_vec(j,:,intval), transpose(phitau) )
        end do
      end do

    else if ( loca == 1 ) then

!     decoupled fields (legacy dfm)

      Fvec1(:,:,:) =  F_vec

    end if

!   compute stress tensor

!   F tensor

    Ften = reshape( Fvec1(:,:,:), [ndim,ndim,ninttau,nintvaltau] )

    do intval = 1, nintvaltau
      do field = 1, ninttau
        Ften(:,:,field,intval) = transpose(Ften(:,:,field,intval) )
      end do
    end do

!   stress tensor

    call stress_tensor_intmodel ( coefficients, stress_ten )

    deallocate ( Fvec1 )
    deallocate ( tau_gp )
    deallocate ( tau_gw )
    deallocate ( mf )
    deallocate ( ivr )
    deallocate ( dpf )
    deallocate ( phitau )
    deallocate ( Ften, Bten )

  end subroutine dfm_stress_tensor_point


! compute stress tensor in all nodes

  subroutine dfm_deriv_stress_tensor ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    use dfm_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    real(dp) :: stressten(3,3)
    integer :: j, intval, field, loca, node, i, k


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_dfm ( coefficients )
      call set_globals_intmodel ( coefficients )

      call check ( coefficients, 'dfm_deriv_stress_tensor', &
        indexarray=[361], minimum=[0], maximum=[1] )

!     allocate arrays

      allocate ( theta(nodalp,ndff) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( posf(ndff*ninttau) )
      allocate ( fn(ndff,ncompf) )
      allocate ( fu(ndff*ninttau), fn2(ndff,ncompf,ninttau) )
      allocate ( Fvec(nodalp,ncompf,ninttau,nintvaltau) )
      allocate ( tau_gp(ninttau,nintvaltau) )
      allocate ( tau_gw(ninttau,nintvaltau) )
      allocate ( mf(ninttau,nintvaltau) )
      allocate ( ivr(ninttau,nintvaltau,2) )
      allocate ( dpf(ninttau,nintvaltau,2) )
      allocate ( phitau(ninttau,ninttau) )
      allocate ( stressvec(nodalp,ncompt) )
      allocate ( Ften(ndim,ndim,ninttau,nintvaltau) )
      allocate ( Bten(ndim,ndim,ninttau,nintvaltau) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function ( shapefuncf, xrnod, theta )

!     tau Gauss points

      call dfm_tau_gauss ( coefficients, tau_gp, tau_gw, phitau )

!     memory function in the tau Gauss points

      call memory_function ( coefficients, tau_gp, mf )

    end if

    loca = coefficients%i(361)

    if ( loca == 0 ) then

!     coupled fields

      do intval = 1, nintvaltau

        call get_sysvector ( mesh, oldvectors%p(2)%p, &
          oldvectors%s2(1)%p(1,intval), elgrp, elem, fu, posu=posf, &
          layer=layer )

        fn2(:,1,:) = reshape( fu, [ndff,ninttau] )

        do j = 2, ncompf
          fn2(:,j,:) = &
                reshape( oldvectors%s2(1)%p(j,intval)%u(posf), [ndff,ninttau] )
        end do

!       interpolate to tau Gauss nodes

        do j = 1, ncompf
          fn2(:,j,:) = matmul ( fn2(:,j,:), transpose(phitau) )
        end do

        do field = 1, ninttau
          Fvec(:,:,field,intval) = matmul( theta, fn2(:,:,field) )
        end do

      end do

    else if ( loca == 1 ) then

!     decoupled fields (legacy dfm)

      do intval = 1, nintvaltau
        do field = 1, ninttau

          call get_sysvector ( mesh,  oldvectors%p(2)%p, &
            oldvectors%s3(1)%p(1,field,intval), elgrp, elem, fn(:,1), &
            posu=posf(1:ndff), layer=layer )

          do j = 2, ncompf
            fn(:,j) = oldvectors%s3(1)%p(j,field,intval)%u(posf(1:ndff))
          end do

          Fvec(:,:,field,intval) = matmul( theta, fn )

        end do
      end do
    end if

!   compute stress tensor for all nodes

    do node = 1, nodalp

!     F tensor

      Ften = reshape( Fvec(node,:,:,:), [ndim,ndim,ninttau,nintvaltau] )

      do intval = 1, nintvaltau
        do field = 1, ninttau
          Ften(:,:,field,intval) = transpose(Ften(:,:,field,intval) )
        end do
      end do

!     stress tensor

      call stress_tensor_intmodel ( coefficients, stressten )

      k = 0
      do j = 1, ndim
        do i = j, ndim ! symmetric tensor
          k = k + 1
          stressvec(node,k) = stressten(i,j)
        end do
      end do

      if ( ndim == 2 .and. ncompt == 4 ) then
        stressvec(node,4) = stressten(3,3)
      end if

    end do

    elemvec = reshape ( stressvec, [nodalp*ncompt] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( posf )
      deallocate ( fn )
      deallocate ( fu, fn2 )
      deallocate ( Fvec )
      deallocate ( tau_gp )
      deallocate ( tau_gw )
      deallocate ( mf )
      deallocate ( ivr )
      deallocate ( dpf )
      deallocate ( phitau )
      deallocate ( stressvec )
      deallocate ( Ften, Bten )

    end if

  end subroutine dfm_deriv_stress_tensor


! Stress projection on the standard discretized space for F.

  subroutine stress_projection_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use dfm_globals_m

    implicit none

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    real(dp) :: stressten(3,3)
    integer :: i, j, ip, loca, intval, field, k


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_dfm_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_globals_intmodel ( coefficients )

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

!   Some integral model data

    if ( first ) then

!     tau Gauss points

      call dfm_tau_gauss ( coefficients, tau_gp, tau_gw, phitau )

!     memory function in the tau Gauss points

      call memory_function ( coefficients, tau_gp, mf )

    end if


    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if


    if ( matrix ) then

      do i = 1, ndff
        do j = 1, ndff
          elemmat(i,j) = sum ( theta(:,i) * theta(:,j) * detF * wg )
        end do
      end do

    end if


    if ( vector ) then

!     get fields

      loca = coefficients%i(361)

      if ( loca == 0 ) then

!       coupled fields

        do intval = 1, nintvaltau

          call get_sysvector ( mesh, oldvectors%p(2)%p, &
            oldvectors%s2(2)%p(1,intval), &
            elgrp, elem, fu, posu=posf, layer=layer )

          fn2(:,1,:) = reshape( fu, [ndff,ninttau] )

          do j = 2, ncompf
            fn2(:,j,:) = &
                reshape( oldvectors%s2(2)%p(j,intval)%u(posf), [ndff,ninttau] )
          end do

!         interpolate to tau Gauss nodes

          do j = 1, ncompf
            fn2(:,j,:) = matmul ( fn2(:,j,:), transpose(phitau) )
          end do

          do field = 1, ninttau
            Fvec(:,:,field,intval) = matmul( theta, fn2(:,:,field) )
          end do

        end do

      else if ( loca == 1 ) then

!       decoupled fields (legacy dfm)

        do intval = 1, nintvaltau
          do field = 1, ninttau

            call get_sysvector ( mesh, oldvectors%p(2)%p, &
              oldvectors%s3(1)%p(1,field,intval), elgrp, elem, fn(:,1), &
              posu=posf(1:ndff), layer=layer )

            do j = 2, ncompf
              fn(:,j) = oldvectors%s3(1)%p(j,field,intval)%u(posf(1:ndff))
            end do

            Fvec(:,:,field,intval) = matmul( theta, fn )

          end do
        end do

      end if

!     compute stress tensor for all integration points

      do ip = 1, ninti

!       F tensor

        Ften = reshape( Fvec(ip,:,:,:), [ndim,ndim,ninttau,nintvaltau] )

        do intval = 1, nintvaltau
          do field = 1, ninttau
            Ften(:,:,field,intval) = transpose(Ften(:,:,field,intval) )
          end do
        end do

!       stress tensor

        call stress_tensor_intmodel ( coefficients, stressten )

        k = 0
        do j = 1, ndim
          do i = j, ndim ! symmetric tensor
            k = k + 1
            stressvec(ip,k) = stressten(i,j)
          end do
        end do

        if ( ndim == 2 .and. ncompt == 4 ) then
          stressvec(ip,4) = stressten(3,3)
        end if

      end do

      do j = 1, ncompt
        do i = 1, ndff
          work2(i,j) = sum ( theta(:,i) * stressvec(:,j) * detF * wg )
        end do
      end do

      elemvec = reshape ( work2, [ ndff * ncompt ] )

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_dfm_elem ( last, coefficients )

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

      allocate ( posf(ndff*ninttau) )
      allocate ( fn(ndff,ncompf) )
      allocate ( fu(ndff*ninttau), fn2(ndff,ncompf,ninttau) )
      allocate ( Fvec(ninti,ncompf,ninttau,nintvaltau) )
      allocate ( tau_gp(ninttau,nintvaltau) )
      allocate ( tau_gw(ninttau,nintvaltau) )
      allocate ( mf(ninttau,nintvaltau) )
      allocate ( ivr(ninttau,nintvaltau,2) )
      allocate ( dpf(ninttau,nintvaltau,2) )
      allocate ( stressvec(ninti,ncompt) )
      allocate ( work2(ndff,ncompt) )
      allocate ( Ften(ndim,ndim,ninttau,nintvaltau) )
      allocate ( Bten(ndim,ndim,ninttau,nintvaltau) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( posf )
      deallocate ( fn )
      deallocate ( fu, fn2 )
      deallocate ( Fvec )
      deallocate ( tau_gp )
      deallocate ( tau_gw )
      deallocate ( mf )
      deallocate ( ivr )
      deallocate ( dpf )
      deallocate ( stressvec )
      deallocate ( work2 )
      deallocate ( Ften )
      deallocate ( Bten )

    end subroutine deallocate_arrays

  end subroutine stress_projection_elem


! compute projected stress tensor in all nodes

  subroutine dfm_deriv_stress_tensor_proj ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use dfm_globals_m

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
      call set_globals_dfm ( coefficients )
      call set_globals_intmodel ( coefficients )

      call check ( coefficients, 'dfm_deriv_stress_tensor_proj', &
        indexarray=[361], minimum=[0], maximum=[1] )

!     allocate arrays

      allocate ( theta(nodalp,ndff) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( s(ndff,ncompt), pos_s(ndff) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function ( shapefuncf, xrnod, theta )

    end if

    call get_sysvector ( mesh, oldvectors%p(3)%p, oldvectors%s1(1)%p(1), &
      elgrp, elem, s(:,1), posu=pos_s, layer=layer )

    do j = 2, ncompt
      s(:,j) = oldvectors%s1(1)%p(j)%u(pos_s)
    end do

    elemvec = reshape ( matmul ( theta, s ), [nodalp*ncompt] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( s, pos_s )

    end if

  end subroutine dfm_deriv_stress_tensor_proj


! memory function m(tau) in tau points (age axis).
! NOTE: set_globals_intmodel needs to be called before this routine is called.
! NOTE: dfm_globals_m is used to transfer the following variables: spectrum,
!       nummodes

  subroutine memory_function ( coefficients, tau, memf )

    use dfm_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   a matrix of tau points where the memory function must be computed
    real(dp), dimension(:,:), intent(in) :: tau

!   the memory function m(tau) for a matrix of points in tau (age) space.
    real(dp), dimension(:,:), intent(out) :: memf

    integer :: m
    real(dp), dimension(nummodes) :: Gmd, lmd, nmd

    select case ( spectrum )
    case(0)  ! discrete spectrum
      Gmd = coefficients%ra(2)%a(1:nummodes)
      lmd = coefficients%ra(3)%a(1:nummodes)
      memf = 0
      do m = 1, nummodes
        memf = memf + Gmd(m)/lmd(m) * exp(-tau/lmd(m))
      end do
    case(1)  ! Mittag-Leffler
      write(*,'(/a/a/)') 'Error in memory_function:', &
        ' Mittag-Leffler memory function not yet available '
      stop
    case(2)  ! discrete modes + Carreau model modes for G(tau)
             ! G(t) = sum_k Gk*exp(-t/lk) + sum_m Gm/(1+(t/lm)^2)^(nm/2)
      Gmd = coefficients%ra(2)%a(1:nummodes)
      lmd = coefficients%ra(3)%a(1:nummodes)
      nmd = coefficients%ra(4)%a(1:nummodes)
      memf = 0
      do m = 1, numdismodes ! first numdismodes modes are discrete
        memf = memf + Gmd(m)/lmd(m) * exp(-tau/lmd(m))
      end do
      do m = numdismodes+1, nummodes
        memf = memf + Gmd(m)/lmd(m) * &
                 nmd(m) * tau/lmd(m) / ( 1 + (tau/lmd(m))**2 ) ** (nmd(m)/2+1)
      end do
    case default
      call errormsg_case_default ( 'memory_function', 'spectrum', &
        int_value=spectrum )
    end select

  end subroutine memory_function


! derivative of memory function dm/dtau in tau points (age axis).
! NOTE: set_globals_intmodel needs to be called before this routine is called.
! NOTE: dfm_globals_m is used to transfer the following variables: spectrum,
!       nummodes

  subroutine memory_function_derivative ( coefficients, tau, dmemf )

    use dfm_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   a matrix of tau points where the memory function must be computed
    real(dp), dimension(:,:), intent(in) :: tau

!   the derivative of memory function dm/dtau(tau) for a matrix of points in
!   tau (age) space.
    real(dp), dimension(:,:), intent(out) :: dmemf

    integer :: m
    real(dp), dimension(nummodes) :: Gmd, lmd, nmd

    select case ( spectrum )
    case(0)  ! discrete spectrum
      Gmd = coefficients%ra(2)%a(1:nummodes)
      lmd = coefficients%ra(3)%a(1:nummodes)
      dmemf = 0
      do m = 1, nummodes
        dmemf = dmemf - Gmd(m)/lmd(m)**2 * exp(-tau/lmd(m))
      end do
    case(1)  ! Mittag-Leffler
      write(*,'(/a/a/)') 'Error in memory_function:', &
        ' Mittag-Leffler memory function not yet available '
      stop
    case(2)  ! discrete modes + Carreau model modes for G(tau)
             ! G(t) = sum_k Gk*exp(-t/lk) + sum_m Gm/(1+(t/lm)^2)^(nm/2)
      Gmd = coefficients%ra(2)%a(1:nummodes)
      lmd = coefficients%ra(3)%a(1:nummodes)
      nmd = coefficients%ra(4)%a(1:nummodes)
      dmemf = 0
      do m = 1, numdismodes ! first numdismodes modes are discrete
        dmemf = dmemf - Gmd(m)/lmd(m)**2 * exp(-tau/lmd(m))
      end do
      do m = numdismodes+1, nummodes
        dmemf = dmemf + Gmd(m)/lmd(m)**2 * &
                 nmd(m) * ( 1 - (nmd(m)+1)*(tau/lmd(m))**2) &
                    / ( 1 + (tau/lmd(m))**2 ) ** (nmd(m)/2+2)
      end do
    case default
      call errormsg_case_default ( 'memory_function_derivative', &
        'spectrum', int_value=spectrum )
    end select

  end subroutine memory_function_derivative


! damping function in tau points (age axis).
! NOTE: set_globals_intmodel needs to be called before this routine is called.
! NOTE: dfm_globals_m is used to transfer the following variables: dampingf
! NOTE: We assume h1+h2=1 in equilibrium.

  subroutine damping_function ( coefficients, invar, dampf )

    use dfm_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   the invariants I1=invar(:,:,1), I2=invar(:,:,2) for a matrix
!   of points in tau (age) space.
    real(dp), dimension(:,:,:), intent(in) :: invar

!   the damping functions h1=dampf(:,:,1), h2=dampf(:,:,2) for a matrix
!   of points in tau (age) space.
    real(dp), dimension(:,:,:), intent(out) :: dampf

    real(dp) :: a, alpha, beta

    select case ( dampingf )
    case(1)  ! Lodge rubberlike liquid
      dampf(:,:,1) = 1
      dampf(:,:,2) = 0
    case(2)  ! used by McKinley in fractional K-BKZ
      a = coefficients%r(303)
      dampf(:,:,1) = 1 / ( 1 + a*(invar(:,:,1)-3) )
      dampf(:,:,2) = 0
    case(3)  ! PSM
      alpha = coefficients%r(304)
      beta  = coefficients%r(305)
      dampf(:,:,1) = alpha / &
                     ( alpha - 3 + beta*invar(:,:,1) + (1-beta)*invar(:,:,2) )
      dampf(:,:,2) = 0
    case default
      call errormsg_case_default ( 'damping_function', 'dampingf', &
        int_value=dampingf )
    end select

  end subroutine damping_function


! derivative of damping function in tau points (age axis).
! NOTE: set_globals_intmodel needs to be called before this routine is called.
! NOTE: dfm_globals_m is used to transfer the following variables: dampingf
! NOTE: We assume h1+h2=1 in equilibrium.

  subroutine damping_function_derivative ( coefficients, invar, ddampf )

    use dfm_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   the invariants I1=invar(:,:,1), I2=invar(:,:,2) for a matrix
!   of points in tau (age) space.
    real(dp), dimension(:,:,:), intent(in) :: invar

!   the derivatives of damping functions
!     dhi/dIj=ddampf(:,:,i,j)
!   for a matrix of points in tau (age) space.
    real(dp), dimension(:,:,:,:), intent(out) :: ddampf

    real(dp) :: a, alpha, beta

    select case ( dampingf )
    case(1)  ! Lodge rubberlike liquid
      ddampf = 0
    case(2)  ! used by McKinley in fractional K-BKZ
      a = coefficients%r(303)
      ddampf(:,:,1,1) = - a / ( 1 + a*(invar(:,:,1)-3) ) ** 2
      ddampf(:,:,1,2) = 0
      ddampf(:,:,2,:) = 0
    case(3)  ! PSM
      alpha = coefficients%r(304)
      beta  = coefficients%r(305)
      ddampf(:,:,1,1) = - alpha * beta / &
                ( alpha - 3 + beta*invar(:,:,1) + (1-beta)*invar(:,:,2) ) ** 2
      ddampf(:,:,1,2) = - alpha * (1-beta) / &
                ( alpha - 3 + beta*invar(:,:,1) + (1-beta)*invar(:,:,2) ) ** 2
      ddampf(:,:,2,:) = 0
    case default
      call errormsg_case_default ( 'damping_function_derivative', &
        'dampingf', int_value=dampingf )
    end select

  end subroutine damping_function_derivative


! stress tensor and the (optional) double stress tensor and the (optional)
! relaxation tensor for integral models.
! NOTE: set_globals_intmodel needs to be called before this routine is called.
! NOTE: dfm_globals_m is used to transfer the following variables:
!       ndim, dpf, ddpf, nintval, ninttau, Bten, Ften, ivr, spectrum, nummodes,
!       dampingf, mf, dmf, tau_gw, ncompt
! NOTE: We assume h1+h2=1 in equilibrium.

  subroutine stress_tensor_intmodel ( coefficients, stressten, relaxstressten, &
    dblstressten )

    use dfm_globals_m

    type(coefficients_t), intent(in) :: coefficients
    real(dp), dimension(3,3), intent(out) :: stressten
    real(dp), dimension(3,3), optional, intent(out) :: relaxstressten
    real(dp), dimension(ndim,ndim,ndim,ndim), optional, intent(out) :: &
      dblstressten

    integer :: intval, field, i, j, k ,m
    real(dp) :: I1, I2

!   B tensor = F.F^T and invariants

    do intval = 1, nintvaltau
      do field = 1, ninttau

!       B
        Bten(:,:,field,intval) = matmul ( Ften(:,:,field,intval), &
                      transpose(Ften(:,:,field,intval) ) ) ! F.F^T
!Btensor Use the following line instead of the two lines above if Ften already
!Btensor defines the B-tensor.
!Btensor         Bten = Ften

!       invariants

!       I1 = tr B
        I1 = 0
        do i = 1, ndim
          I1 = I1 + Bten(i,i,field,intval)
        end do
        if ( ndim == 2 ) I1 = I1 + 1
        ivr(field,intval,1) = I1
!       I2= (I1-tr B^2)/2
        I2 = ( I1**2 - sum ( Bten(:,:,field,intval) ** 2 ) ) / 2
        if ( ndim == 2 ) I2 = I2 - 0.5_dp
        ivr(field,intval,2) = I2

      end do
    end do

!   damping function

    call damping_function ( coefficients, ivr, dpf )

!   stress tensor

    stressten = 0

    do j = 1, ndim
      do i = j, ndim ! symmetric tensor: only upper triangle
        stressten(i,j) = &
                    sum ( mf * dpf(:,:,1) * Bten(i,j,:,:) * tau_gw )
      end do
    end do

    do j = 1, ndim
      do i = j+1, ndim ! symmetric tensor: transpose
        stressten(j,i) = stressten(i,j)
      end do
    end do

    if ( ndim == 2 .and. ncompt == 4 ) then
      stressten(3,3) = sum ( mf * dpf(:,:,1) * tau_gw ) ! zz
    end if

!   optional relaxation stress tensor

    if ( present(relaxstressten) ) then

!     relaxation stress tensor

      relaxstressten = 0

      do j = 1, ndim
        do i = j, ndim ! symmetric tensor: only upper triangle
          relaxstressten(i,j) = &
                      sum ( dmf * dpf(:,:,1) * Bten(i,j,:,:) * tau_gw )
        end do
      end do

      do j = 1, ndim
        do i = j+1, ndim ! symmetric tensor: transpose
          relaxstressten(j,i) = relaxstressten(i,j)
        end do
      end do

      if ( ndim == 2 .and. ncompt == 4 ) then
        relaxstressten(3,3) = sum ( dmf * dpf(:,:,1) * tau_gw ) ! zz
      end if

    end if

!   optional double stress tensor

    if ( present(dblstressten) ) then

!     damping function

      call damping_function_derivative ( coefficients, ivr, ddpf )

!     double stress tensor (fourth-order tensor)
!     TODO: use symmetry for efficiency

      dblstressten = 0

      do j = 1, ndim
        do i = 1, ndim
          do k = 1, ndim
            do m = 1, ndim
              dblstressten(i,j,k,m) = 2 * sum ( mf * &
                      ddpf(:,:,1,1) * Bten(i,j,:,:) * Bten(k,m,:,:) * tau_gw )
            end do
          end do
        end do
      end do

    end if

  end subroutine stress_tensor_intmodel


! Internal element routine for implicit terms of an integral model in the
! momentum balance (convection implicit). For first-order and second-order
! integration.

  subroutine divtau_implicit_intmodel_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use dfm_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    logical :: validm
    integer :: i5, i6, i7, ip
    integer :: N, M, i, j, id
    integer :: loca
    real(dp) :: deltat
    real(dp), dimension(1,1) :: mf0, tau0
    real(dp), dimension(1,1,2) :: invar0, dpf0
    real(dp), dimension(3,3) :: stressten, relaxstressten


!   check coefficients

    call check ( coefficients, 'divtau_implicit_intmodel_elem', &
      indexarray=[48], minimum=[0], maximum=[1] )


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_devssg_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_dfm_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_globals_intmodel ( coefficients )

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

!   Some integral model data

    if ( first ) then

!     tau Gauss points

      call dfm_tau_gauss ( coefficients, tau_gp, tau_gw, phitau )

!     memory function in the tau Gauss points

      call memory_function ( coefficients, tau_gp, mf )
      call memory_function_derivative ( coefficients, tau_gp, dmf )

!     zero values m(0), h1(3,3), h2(3,3)

      tau0 = 0
      call memory_function ( coefficients, tau0, mf0 )
      invar0 = 3
      call damping_function ( coefficients, invar0, dpf0 )

!     mf0h0 = m(0) ( h1(3,3) - h2(3,3) )
      mf0h0 = mf0(1,1) * ( dpf0(1,1,1) - dpf0(1,1,2) )

    end if

!   Some checking

    if ( first ) then

!     first element in this group

      validm = intmodel==0 .and. any( dampingf==[1,2] )

      if ( .not. validm ) then
        write(*,'(/a/a/4(a,i0)/)') 'Error in divtau_implicit_intmodel_elem:', &
          ' Model has not been implemented:', &
          ' intmodel = ', intmodel, ' spectrum = ', spectrum, &
          ' nummodes = ', nummodes, ' dampingf = ', dampingf
        stop
      end if

      if ( size(elemmat,2) /= ndim*ndf ) then
        write(*,'(/3(a/))') 'Error in divtau_implicit_intmodel_elem:', &
          ' the degrees of freedom in physqcol ', &
          ' must be (u), i.e. the velocity only'
        stop
      end if

    end if


!   start of the element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    call shape_derivative ( dtheta, Finv, dthetadx )


!   get projected tau (stress tensor) and compute nabla tau

    call get_projected_stress

!   get deformation gradient fields at previous time step

    call get_fn

!   stress (integral)

    call compute_stress_integral


!   set time step

    deltat = coefficients%r(301)


!   get mesh velocity at the current time

    if ( coefficients%i(48) == 1 ) then

!     standard ALE

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, u )

      tmp = reshape ( u, [ndf,ndim] )

      uvecmeshnp1 = matmul ( phi, tmp )

    end if


!   set tau_n in tensor format

    if ( coorsys <= 1 ) then

      taunten(:,1,1) = stressvec(:,1)
      taunten(:,1,2) = stressvec(:,2)
      taunten(:,2,1) = taunten(:,1,2)
      taunten(:,2,2) = stressvec(:,3)

    else if ( coorsys == 2 ) then

      taunten(:,1,1) = stressvec(:,1)
      taunten(:,1,2) = stressvec(:,2)
      taunten(:,1,3) = stressvec(:,3)
      taunten(:,2,1) = taunten(:,1,2)
      taunten(:,2,2) = stressvec(:,4)
      taunten(:,2,3) = stressvec(:,5)
      taunten(:,3,1) = taunten(:,1,3)
      taunten(:,3,2) = taunten(:,2,3)
      taunten(:,3,3) = stressvec(:,6)

    end if

!   set grad tau in tensor format

    if ( coorsys <= 1 ) then

      do i = 1, ndim
        gradtauten(:,1,1,i) = gradtaun(:,1,i)
        gradtauten(:,1,2,i) = gradtaun(:,2,i)
        gradtauten(:,2,1,i) = gradtauten(:,1,2,i)
        gradtauten(:,2,2,i) = gradtaun(:,3,i)
      end do

    else if ( coorsys == 2 ) then

      do i = 1, ndim
        gradtauten(:,1,1,i) = gradtaun(:,1,i)
        gradtauten(:,1,2,i) = gradtaun(:,2,i)
        gradtauten(:,1,3,i) = gradtaun(:,3,i)
        gradtauten(:,2,1,i) = gradtauten(:,1,2,i)
        gradtauten(:,2,2,i) = gradtaun(:,4,i)
        gradtauten(:,2,3,i) = gradtaun(:,5,i)
        gradtauten(:,3,1,i) = gradtauten(:,1,3,i)
        gradtauten(:,3,2,i) = gradtauten(:,2,3,i)
        gradtauten(:,3,3,i) = gradtaun(:,6,i)
      end do

    end if


    if ( matrix ) then

      do M = 1, ndf

!       dphidx_m^M tau_mi

        do ip = 1, ninti
          work2(ip,:) = matmul ( dphidx(ip,M,:), taunten(ip,:,:) )
        end do

        do N = 1, ndf

!         dphidx_k^N dphidx_m^M tau_km

          do ip = 1, ninti
            work1(ip) = dot_product ( dphidx(ip,N,:), work2(ip,:) )
          end do

!         dphidx_k^N dc_ki/dx_j

          do ip = 1, ninti
            do id = 1, ndim
              work6(ip,:,id) = matmul ( dphidx(ip,N,:), gradtauten(ip,:,:,id) )
            end do
          end do

          do i = 1, ndim
            do j = 1, ndim

!             dphidx_k^N * dphidx_m^N * c_mi * delta_jk

              work = dphidx(:,N,j) * work2(:,i)

!             dphidx_k^N * dphidx_m^N * c_km * delta_ij

              if ( i == j ) then
                work = work + work1
              end if

!             - phi_M dphidx_k^N dc_ki/dx_j

              work = work - phi(:,M) * work6(:,i,j)

              work5(N,M,i,j) = sum ( work * detF * wg )

            end do
          end do

        end do

      end do

      if ( I1dep .or. I2dep ) then

!       contribution of the double stress tensor to the matrix

        do N = 1, ndf
          do M = 1, ndf
            do i = 1, ndim
              do j = 1, ndim

!               dphidx_k^N * K_kijm * dphidx_m^M

                do ip = 1, ninti
                  work1(ip) = dot_product ( dphidx(ip,N,:), &
                       matmul( doublestressten(ip,:,i,j,:), dphidx(ip,M,:) ) )
                end do

                work5(N,M,i,j) = work5(N,M,i,j) + sum ( work1 * detF * wg )

              end do
            end do
          end do
        end do

      end if

      work5 = deltat * work5


      if ( coorsys == 1 ) then

!       axisymmetric

        work1 = 2 * stressvec(:,4) / xg(:,2) ** 2

        do N = 1, ndf
          do M = 1, ndf
!           phi^N 2tau_tt/r^2 phi^M
            work4(N,M) = sum ( phi(:,N) * work1 * phi(:,M) * detF * wg )
          end do
        end do

        work4 = deltat * work4

        do i = 1, ndim
          work2(:,i) = gradtaun(:,4,i) / xg(:,2)
        end do

        do N = 1, ndf
          do M = 1, ndf
            do j = 1, ndim
!             - phi^N phi^M dc_tt/dx_j 1/r
              work7(N,M,j) = &
                - sum ( phi(:,N) * work2(:,j) * phi(:,M) * detF * wg )
            end do
          end do
        end do

        work7 = deltat * work7

      end if

!     put into element matrix

      i5 = ndf    ! u
      i6 = 2*ndf  ! v
      i7 = 3*ndf  ! w

      if ( coorsys == 1 ) then
        work5(:,:,2,2) = work5(:,:,2,2) + work4
        work5(:,:,2,:) = work5(:,:,2,:) + work7
      end if

      if ( coorsys <= 1 ) then

        elemmat( 1:i5, 1:i5 )       = work5(:,:,1,1)
        elemmat( 1:i5, i5+1:i6 )    = work5(:,:,1,2)
        elemmat( i5+1:i6, 1:i5 )    = work5(:,:,2,1)
        elemmat( i5+1:i6, i5+1:i6 ) = work5(:,:,2,2)

      else if ( coorsys == 2 ) then

        elemmat( 1:i5, 1:i5 )       = work5(:,:,1,1)
        elemmat( 1:i5, i5+1:i6 )    = work5(:,:,1,2)
        elemmat( 1:i5, i6+1:i7 )    = work5(:,:,1,3)
        elemmat( i5+1:i6, 1:i5 )    = work5(:,:,2,1)
        elemmat( i5+1:i6, i5+1:i6 ) = work5(:,:,2,2)
        elemmat( i5+1:i6, i6+1:i7 ) = work5(:,:,2,3)
        elemmat( i6+1:i7, 1:i5 )    = work5(:,:,3,1)
        elemmat( i6+1:i7, i5+1:i6 ) = work5(:,:,3,2)
        elemmat( i6+1:i7, i6+1:i7 ) = work5(:,:,3,3)

      end if

    end if

    if ( vector ) then

!     - (nabla v)^T: ( rhs )

      rhs = stressvec + deltat * relaxstressvec

!     add mesh velocity term to the right-hand side

      if ( coefficients%i(48) == 1 ) then

        do ip = 1, ninti
          rhs(ip,:) = rhs(ip,:) &
                           + matmul ( gradtaun(ip,:,:), uvecmeshnp1(ip,:) )
        end do

      end if

      if ( coorsys <= 1 ) then

        rhsten(:,1,1) = rhs(:,1) + mf0h0 * deltat
        rhsten(:,1,2) = rhs(:,2)
        rhsten(:,2,1) = rhsten(:,1,2)
        rhsten(:,2,2) = rhs(:,3) + mf0h0 * deltat

      else if ( coorsys == 2 ) then

        rhsten(:,1,1) = rhs(:,1) + mf0h0 * deltat
        rhsten(:,1,2) = rhs(:,2)
        rhsten(:,1,3) = rhs(:,3)
        rhsten(:,2,1) = rhsten(:,1,2)
        rhsten(:,2,2) = rhs(:,4) + mf0h0 * deltat
        rhsten(:,2,3) = rhs(:,5)
        rhsten(:,3,1) = rhsten(:,1,3)
        rhsten(:,3,2) = rhsten(:,2,3)
        rhsten(:,3,3) = rhs(:,5) + mf0h0 * deltat

      end if

      if ( coorsys == 1 ) then
!       work1 = rhsten_tt
        work1 = rhs(:,4) + mf0h0 * deltat
      end if

!     - (nabla v)^T:tau

      do ip = 1, ninti
        work9(ip,:,:) = matmul ( dphidx(ip,:,:), rhsten(ip,:,:) )
      end do
      if ( coorsys == 1 ) then
        do N = 1, ndf
          work9(:,N,2) = work9(:,N,2) + work1 * phi(:,N) / xg(:,2)
        end do
      end if
      do j = 1, ndim
        work10(:,j) = - matmul ( detF * wg, work9(:,:,j) )
      end do
      elemvec = reshape ( work10, [ ndf*ndim ] )

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_devssg_elem ( last, coefficients )

    call unset_dfm_elem ( last, coefficients )


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

      allocate ( work(ninti), work2(ninti,ndim) )
      allocate ( tmp(ndf,ndim) )
      allocate ( work4(ndf,ndf), work1(ninti) )
      allocate ( work5(ndf,ndf,ndim,ndim) )
      allocate ( work6(ninti,ndim,ndim) )
      allocate ( work7(ndf,ndf,ndim) )
      allocate ( work9(ninti,ndf,ndim), work10(ndf,ndim) )
      allocate ( gradtaun(ninti,ncompt,ndim) )
      allocate ( taun(ninti,ncompt) )

      allocate ( u(ndf*ndim), uvecmeshnp1(ninti,ndim) )

      allocate ( posf(ndff*ninttau) )
      allocate ( fn(ndff,ncompf) )
      allocate ( fu(ndff*ninttau), fn2(ndff,ncompf,ninttau) )
      allocate ( Fvec(ninti,ncompf,ninttau,nintvaltau) )
      allocate ( tau_gp(ninttau,nintvaltau) )
      allocate ( tau_gw(ninttau,nintvaltau) )
      allocate ( mf(ninttau,nintvaltau) )
      allocate ( dmf(ninttau,nintvaltau) )
      allocate ( ivr(ninttau,nintvaltau,2) )
      allocate ( dpf(ninttau,nintvaltau,2) )
      allocate ( stressvec(ninti,ncompt) )
      allocate ( relaxstressvec(ninti,ncompt) )
      allocate ( Ften(ndim,ndim,ninttau,nintvaltau) )
      allocate ( Bten(ndim,ndim,ninttau,nintvaltau) )

      allocate ( s(ndff,ncompt), pos_s(ndff) )

      allocate ( taunten(ninti,ndim,ndim) )
      allocate ( gradtauten(ninti,ndim,ndim,ndim) )
      allocate ( rhsten(ninti,ndim,ndim) )
      allocate ( rhs(ninti,ncompt) )

      if ( I1dep .or. I2dep ) then
        allocate ( ddpf(ninttau,nintvaltau,2,2) )
        allocate ( doublestressten(ninti,ndim,ndim,ndim,ndim) )
      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work, work2 )
      deallocate ( tmp )
      deallocate ( work4, work1 )
      deallocate ( work5 )
      deallocate ( work6 )
      deallocate ( work7 )
      deallocate ( work9, work10 )
      deallocate ( gradtaun )
      deallocate ( taun )

      deallocate ( u, uvecmeshnp1 )

      deallocate ( posf )
      deallocate ( fn )
      deallocate ( fu, fn2 )
      deallocate ( Fvec )
      deallocate ( tau_gp )
      deallocate ( tau_gw )
      deallocate ( mf )
      deallocate ( dmf )
      deallocate ( ivr )
      deallocate ( dpf )
      deallocate ( stressvec )
      deallocate ( relaxstressvec )
      deallocate ( Ften )
      deallocate ( Bten )

      deallocate ( s, pos_s )

      deallocate ( taunten )
      deallocate ( gradtauten )
      deallocate ( rhsten )
      deallocate ( rhs )

      if ( I1dep .or. I2dep ) then
        deallocate ( ddpf )
        deallocate ( doublestressten )
      end if

    end subroutine deallocate_arrays


!   get Fn

    subroutine get_fn

      integer :: intval, j, field

      loca = coefficients%i(361)

      if ( loca == 0 ) then

!       coupled fields

        do intval = 1, nintvaltau

          call get_sysvector ( mesh, oldvectors%p(2)%p, &
            oldvectors%s2(2)%p(1,intval), &
            elgrp, elem, fu, posu=posf, layer=layer )

          fn2(:,1,:) = reshape( fu, [ndff,ninttau] )

          do j = 2, ncompf
            fn2(:,j,:) = &
                reshape( oldvectors%s2(2)%p(j,intval)%u(posf), [ndff,ninttau] )
          end do

!         interpolate to tau Gauss nodes

          do j = 1, ncompf
            fn2(:,j,:) = matmul ( fn2(:,j,:), transpose(phitau) )
          end do

          do field = 1, ninttau
            Fvec(:,:,field,intval) = matmul( theta, fn2(:,:,field) )
          end do

        end do

      else if ( loca == 1 ) then

!       decoupled fields (legacy dfm)

        do intval = 1, nintvaltau
          do field = 1, ninttau

            call get_sysvector ( mesh, oldvectors%p(2)%p, &
              oldvectors%s3(1)%p(1,field,intval), elgrp, elem, fn(:,1), &
              posu=posf(1:ndff), layer=layer )

            do j = 2, ncompf
              fn(:,j) = oldvectors%s3(1)%p(j,field,intval)%u(posf(1:ndff))
            end do

            Fvec(:,:,field,intval) = matmul( theta, fn )

          end do
        end do

      end if

    end subroutine get_fn


    subroutine get_projected_stress

      integer :: i, j

      call get_sysvector ( mesh, oldvectors%p(3)%p, oldvectors%s1(1)%p(1), &
        elgrp, elem, s(:,1), posu=pos_s, layer=layer )

      do j = 2, ncompt
        s(:,j) = oldvectors%s1(1)%p(j)%u(pos_s)
      end do

      taun = matmul ( theta, s )

!     grad taun term

      do i = 1, ndim
        gradtaun(:,:,i) = matmul ( dthetadx(:,:,i), s )
      end do

    end subroutine get_projected_stress


!   stress integral in all integration points

    subroutine compute_stress_integral

      integer :: ip, intval, field, i, j, k
!     compute stress tensor components for all integration points

      do ip = 1, ninti

!       F tensor

        Ften = reshape( Fvec(ip,:,:,:), [ndim,ndim,ninttau,nintvaltau] )

        do intval = 1, nintvaltau
          do field = 1, ninttau
            Ften(:,:,field,intval) = transpose(Ften(:,:,field,intval) )
          end do
        end do

        if ( I1dep .or. I2dep ) then

!         stress tensor and double stress tensor

          call stress_tensor_intmodel ( coefficients, stressten, &
            relaxstressten, doublestressten(ip,:,:,:,:) )

        else

!         stress tensor

          call stress_tensor_intmodel ( coefficients, stressten, &
            relaxstressten )

        end if

        k = 0
        do j = 1, ndim
          do i = j, ndim ! symmetric tensor
            k = k + 1
            stressvec(ip,k) = stressten(i,j)
            relaxstressvec(ip,k) = relaxstressten(i,j)
          end do
        end do

        if ( ndim == 2 .and. ncompt == 4 ) then
          stressvec(ip,4) = stressten(3,3)
          relaxstressvec(ip,4) = relaxstressten(3,3)
        end if

      end do

    end subroutine compute_stress_integral

  end subroutine divtau_implicit_intmodel_elem


! tau values and weights of the Gauss points (age axis).
! NOTE: set_globals_dfm needs to be called before this routine is called.

  subroutine dfm_tau_gauss ( coefficients, tau_gauss_points, &
    tau_gauss_weights, phibasistau )

    use dfm_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   array of size ninttau x nintvaltau, that give the age Gauss points and
!   weights
    real(dp), dimension(:,:), intent(out) :: tau_gauss_points, tau_gauss_weights

!   basis functions in tau direction (regular nodes)
    real(dp), dimension(:,:), optional, intent(out) :: phibasistau

    integer :: m

!   allocate arrays

    allocate ( phitaunod(ninttau,2) )
    allocate ( taug2(ninttau), wtaug2(ninttau) )

!   Gauss points

    call Gauss_Legendre_line ( ninttau, taug2, wtaug2 )

!   Linear shape functions for the age coordinate

    call shape_line_P1 ( taug2, phitaunod )

!   compute tau gauss points and weights based

    do m = 1, nintvaltau
      tau_gauss_points(:,m) = matmul(phitaunod,coefficients%ra(1)%a(m:m+1))
      tau_gauss_weights(:,m) = &
           wtaug2 * ( coefficients%ra(1)%a(m+1) - coefficients%ra(1)%a(m) ) / 2
    end do

    if ( present(phibasistau) ) then

!     basis functions in tau direction (regular nodes)

      select case ( intpoltau )
      case (1)  ! P1
        call shape_line_P1 ( taug2, phibasistau )
      case (2)  ! P2
        call shape_line_P2 ( taug2, phibasistau )
      case default
        call errormsg_case_default ( 'dfm_tau_gauss', 'intpoltau', &
          int_value=intpoltau )
      end select

    end if

!   deallocate arrays

    deallocate ( phitaunod, taug2, wtaug2 )

  end subroutine dfm_tau_gauss


! set global parameters integral models

  subroutine set_globals_intmodel ( coefficients )

    use dfm_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   integral model

    call check ( coefficients, 'set_globals_intmodel', &
      indexarray=[363,364,366], minimum=[0,0,0], maximum=[0,2,3] )
    call check ( coefficients, 'set_globals_intmodel', &
      indexarray=[365], minimum=[0] )

    intmodel = coefficients%i(363)
    spectrum = coefficients%i(364)
    nummodes = get_coefficient ( coefficients, index=365, default=1 )
    dampingf = get_coefficient ( coefficients, index=366, default=1 )

    call check ( coefficients, 'set_globals_intmodel', &
      indexarray=[367], minimum=[0], maximum=[nummodes] )

    numdismodes = coefficients%i(367)

    ncompt = ( ndim + 1 ) * ndim / 2
    if ( ndim == 2 ) ncompt = ncompt + 1  ! zz-component

    if ( intmodel == 0 ) then

!     separable Rivlin-Sawyers

      if ( spectrum == 0 ) then
!       discrete spectrum
        if ( size(coefficients%ra) < 3 ) then
          write(*,'(/a/a/)') 'Error in set_globals_intmodel:', &
            '  size of ra component in coefficients must be at least 3'
          stop
        end if
        if ( min(size(coefficients%ra(2)%a),size(coefficients%ra(3)%a) ) &
                < nummodes ) then
          write(*,'(/a/a/a,i0/)') 'Error in set_globals_intmodel:', &
            '  size of arrays ra(2)%a and ra(3)%a in coefficients ', &
            '  must be at least equal to the number of discrete modes = ', &
            nummodes
          stop
        end if
      else if ( spectrum == 2 ) then
!       Carreau spectrum
        if ( size(coefficients%ra) < 4 ) then
          write(*,'(/a/a/)') 'Error in set_globals_intmodel:', &
            '  size of ra component in coefficients must be at least 4'
          stop
        end if
        if ( minval([size(coefficients%ra(2)%a),size(coefficients%ra(3)%a),&
                     size(coefficients%ra(4)%a)] ) < nummodes ) then
          write(*,'(/a/a/a,i0/)') 'Error in set_globals_intmodel:', &
            '  size of arrays ra(2)%a, ra(3)%a and ra(4)%a in coefficients ', &
            '  must be at least equal to the number of discrete modes = ', &
            nummodes
          stop
        end if
      end if

      select case ( dampingf )
      case(1)  ! UCM
        I1dep = .false.
        I2dep = .false.
        Binv = .false.
      case(2)  ! used by McKinley in fractional K-BKZ
        I1dep = .true.
        I2dep = .false.
        Binv = .false.
      case(3)  ! PSM
        I1dep = .true.
        I2dep = .true.
        Binv = .false.
      case default
        call errormsg_case_default ( 'set_globals_intmodel', 'dampingf', &
          int_value=dampingf )
      end select

    else

!     non-separable Rivlin-Sawyers

      I1dep = .true.
      I2dep = .true.
      Binv = .true.

    end if

  end subroutine set_globals_intmodel

end module dfm_integralmodel_elements_m

