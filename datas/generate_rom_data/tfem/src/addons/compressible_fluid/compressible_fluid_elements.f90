
! Copyright (C) 2012-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the following constitutive equation for the pressure:
!    .
!    p/K + div u = 0
!
! Starting from the reference state (p=p0,J=1), the solution (in Lagrangian
! sense) is given by:
!
!    p = p0 - K ln J
!
! with:
!   K the compression modulus
!   p0 the reference pressure
!   J the Jacobian (=det F) of the deformation tensor
!

module compressible_fluid_elements_m

  use tfem_elem_m
  use stokes_set_globals_m

  implicit none

contains


! Internal element routine for the pressure constitutive equation
! semi-implicit version (with prediction of convection velocity)

  subroutine pressure_ce_elem ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use compressible_fluid_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, N, M, ip, timeint, isvn, isvnp1
    real(dp) :: L_NMi, Kbulk, deltat, fac1, fac2

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'pressure_ce_elem', ncoefi=350, ncoefr=300, &
        indexarray=[48,62,301,304], minimum=[0,0,1,0], maximum=[1,0,2,1] )

    end if

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

!   coefficients

    timeint = coefficients%i(301)
    isvn    = get_coefficient ( coefficients, index=302, default=1 )
    isvnp1  = get_coefficient ( coefficients, index=303, default=1 )
    deltat  = coefficients%r(251)
    Kbulk   = coefficients%r(252)

!   derivative of shape function for pressure

    call set_shape_function ( shapefuncp, xig, psi, dpsi )

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )
    call shape_derivative ( dpsi, Finv, dpsidx )

!   get velocity vector for convection term

    if ( coefficients%i(304) == 0 ) then

!     specified by user

      call get_sysvector ( mesh, problem, oldvectors%s(isvnp1)%p, &
        elgrp, elem, u, physq=[physqvel], layer=layer )

      tmp = reshape ( u, [ndf,ndim] )

    else if ( coefficients%i(304) == 1 ) then

!     computed on element level

      if ( timeint == 1 ) then

!       first-order, semi-implicit Euler, prediction of unp1=un

        call get_sysvector ( mesh, problem, oldvectors%s(isvn)%p, &
          elgrp, elem, u, physq=[physqvel], layer=layer )

        tmp = reshape ( u, [ndf,ndim] )

      else if ( timeint == 2 ) then

!       second-order, semi-implicit Gear (BDF2), prediction of unp1=2*un-unm1

        call get_sysvector ( mesh, problem, oldvectors%s(isvn)%p, &
          elgrp, elem, u, physq=[physqvel], layer=layer )

        tmp = 2 * reshape ( u, [ndf,ndim] )

        call get_sysvector ( mesh, problem, oldvectors%s(isvn+1)%p, &
          elgrp, elem, u, physq=[physqvel], layer=layer )

        tmp = tmp - reshape ( u, [ndf,ndim] )

      end if

    end if

    uvechat = matmul ( phi, tmp )

!   Eulerian frame or ALE

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity for convection term

      if ( coefficients%i(304) == 0 ) then

!       specified by user

        call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )

        tmp = reshape ( u, [ndf,ndim] )

      else if ( coefficients%i(304) == 1 ) then

!       computed on element level

        if ( timeint == 1 ) then

!         first-order, semi-implicit Euler, prediction of unp1=un

          call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )

          tmp = reshape ( u, [ndf,ndim] )

        else if ( timeint == 2 ) then

!         second-order, semi-implicit Gear (BDF2), prediction of unp1=2*un-unm1

          call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )

          tmp = 2 * reshape ( u, [ndf,ndim] )

          call get_vector ( mesh, problem, oldvectors%v(2)%p, elgrp, elem, u )

          tmp = tmp - reshape ( u, [ndf,ndim] )

        end if

      end if

      uvecmeshhat = matmul ( phi, tmp )

!     (un-ugrid).grad operator

      do ip = 1, ninti
        ungradpsi(ip,:) = &
                  matmul ( dpsidx(ip,:,:), uvechat(ip,:) - uvecmeshhat(ip,:) )
      end do

    else

!     Eulerian frame

      uvecmeshhat = 0

!     un.grad operator

      do ip = 1, ninti
        ungradpsi(ip,:) = matmul ( dpsidx(ip,:,:), uvechat(ip,:) )
      end do

    end if

!   get pressure at tn

    call get_sysvector ( mesh, problem, oldvectors%s(isvn)%p, &
      elgrp, elem, p, physq=[physqpress], layer=layer )

    png = matmul ( psi, p )

    if ( timeint == 2 ) then

!     get pressure at time step n-1

      call get_sysvector ( mesh, problem, oldvectors%s(isvn+1)%p, &
        elgrp, elem, p, physq=[physqpress], layer=layer )

      pnm1g = matmul ( psi, p )

    end if

    fac1 = 1._dp / Kbulk
    fac2 = 1._dp / deltat

    if ( vector ) then

!     right-hand side

      if ( timeint == 1 ) then

!       first-order, semi-implicit Euler

        do N = 1, ndfp
          elemvec(N) = fac1 * fac2 * sum ( psi(:,N) * png * detF * wg )
        end do

      else if ( timeint == 2 ) then

!       second-order, semi-implicit Gear (BDF2)

        do N = 1, ndfp
          elemvec(N) = fac1 * fac2 * sum ( psi(:,N) * &
              ( 2*png - pnm1g/2 ) * detF * wg )
        end do

      end if

    end if

    if ( matrix ) then

!     pressure-pressure part

      if ( timeint == 1 ) then

        do N = 1, ndfp
          do M = 1, ndfp
            elemmat(N,posp(M)) = fac1 * sum ( psi(:,N) * &
                ( fac2 * psi(:,M) + ungradpsi(:,M) ) * detF * wg )
          end do
        end do

      else if ( timeint == 2 ) then

        do N = 1, ndfp
          do M = 1, ndfp
            elemmat(N,posp(M)) = fac1 * sum ( psi(:,N) * &
                ( 3*fac2/2 * psi(:,M) + ungradpsi(:,M) ) * detF * wg )
          end do
        end do

      end if

!     pressure-velocity part

      do N = 1, ndfp
        do M = 1, ndf
          do i = 1, ndim
            if ( coorsys == 1 .and. i == 2 ) then ! axisymmetric
              work = dphidx(:,M,2) + phi(:,M) / xg(:,2)
              L_NMi = sum ( psi(:,N) * work * detF * wg )
            else
              L_NMi = sum ( psi(:,N) * dphidx(:,M,i) * detF * wg )
            end if
            elemmat( N, pos(M,i) ) = L_NMi
          end do
        end do
      end do

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

      allocate ( work(ninti), dpsi(ninti,ndfp,ndim), dpsidx(ninti,ndfp,ndim) )
      allocate ( ungradpsi(ninti,ndfp) )
      allocate ( tmp(ndf,ndim), u(ndim*ndf), uvechat(ninti,ndim) )
      allocate ( uvecmeshhat(ninti,ndim) )
      allocate ( p(ndfp), png(ninti), pnm1g(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work, dpsi, dpsidx )
      deallocate ( ungradpsi )
      deallocate ( tmp, u, uvechat )
      deallocate ( uvecmeshhat )
      deallocate ( p, png, pnm1g )

    end subroutine deallocate_arrays

  end subroutine pressure_ce_elem


! support old first-order interface, which was used in the examples

  subroutine pressure_ce_elem1 ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    call pressure_ce_elem ( mesh, problem, elgrp, elem, matrix, vector, &
      first, last, coefficients, oldvectors, elemmat, elemvec )

  end subroutine pressure_ce_elem1


! Internal element routine for the pressure constitutive equation
! fully implicit version (for differences delta p, delta u)

  subroutine implicit_pressure_ce_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use compressible_fluid_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, N, M, ip, timeint, isvn, isv_iter
    real(dp) :: Kbulk, deltat, fac1, fac2

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'implicit_pressure_ce_elem', ncoefi=350, &
        ncoefr=300, indexarray=[48,62,301,306], minimum=[0,0,3,0], &
        maximum=[1,0,4,1] )

    end if

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

!   coefficients

    timeint  = coefficients%i(301)
    isvn     = get_coefficient ( coefficients, index=302, default=1 )
    isv_iter = get_coefficient ( coefficients, index=305, default=1 )
    deltat   = coefficients%r(251)
    Kbulk    = coefficients%r(252)

!   derivative of shape function for pressure

    call set_shape_function ( shapefuncp, xig, psi, dpsi )

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )
    call shape_derivative ( dpsi, Finv, dpsidx )

!   get velocity vector iteration

    call get_sysvector ( mesh, problem, oldvectors%s(isv_iter)%p, &
      elgrp, elem, u, physq=[physqvel], layer=layer )

    tmp = reshape ( u, [ndf,ndim] )

    uvec_iter = matmul ( phi, tmp )

!   div u_iter

    if ( coorsys == 1 ) then ! axisymmetric
      divgiter = uvec_iter(:,2) / xg(:,2)
    else
      divgiter = 0
    end if
    do i = 1, ndim
      divgiter = divgiter + matmul ( dphidx(:,:,i), tmp(:,i) )
    end do

!   Eulerian frame or ALE

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity for convection term

      if ( coefficients%i(304) == 0 ) then

!       specified by user

        call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )

        tmp = reshape ( u, [ndf,ndim] )

      else if ( coefficients%i(304) == 1 ) then

!       computed on element level

        if ( timeint == 3 ) then

!         first-order BDF1

          call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )

          tmp = reshape ( u, [ndf,ndim] )

        else if ( timeint == 4 ) then

!         second-order BDF2

          call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )

          tmp = 2 * reshape ( u, [ndf,ndim] )

          call get_vector ( mesh, problem, oldvectors%v(2)%p, elgrp, elem, u )

          tmp = tmp - reshape ( u, [ndf,ndim] )

        end if

      end if

      uvecmeshhat = matmul ( phi, tmp )

!     (u-ugrid).grad operator

      do ip = 1, ninti
        unp1gradpsi(ip,:) = &
                  matmul ( dpsidx(ip,:,:), uvec_iter(ip,:) - uvecmeshhat(ip,:) )
      end do

    else

!     Eulerian frame

      uvecmeshhat = 0

!     u.grad operator

      do ip = 1, ninti
        unp1gradpsi(ip,:) = matmul ( dpsidx(ip,:,:), uvec_iter(ip,:) )
      end do

    end if

!   get iteration pressure

    call get_sysvector ( mesh, problem, oldvectors%s(isv_iter)%p, &
      elgrp, elem, p, physq=[physqpress], layer=layer )

    png_iter = matmul ( psi, p )

!   nabla p
    do i = 1, ndim
      gradpiter(:,i) = matmul ( dpsidx(:,:,i), p )
    end do

!   u.nabla p
    unp1gradpiter = matmul ( unp1gradpsi, p )

    if ( coefficients%i(306) == 1 ) then

!     no time-derivative included

      png = 0

    else

!     get pressure at tn

      call get_sysvector ( mesh, problem, oldvectors%s(isvn)%p, &
        elgrp, elem, p, physq=[physqpress], layer=layer )

      png = matmul ( psi, p )

    end if

    if ( timeint == 4 ) then

      if ( coefficients%i(306) == 1 ) then

!       no time-derivative included

        pnm1g = 0

      else

!       get pressure at time step n-1

        call get_sysvector ( mesh, problem, oldvectors%s(isvn+1)%p, &
          elgrp, elem, p, physq=[physqpress], layer=layer )

        pnm1g = matmul ( psi, p )

      end if

    end if

    fac1 = 1._dp / Kbulk
    if ( coefficients%i(306) == 1 ) then
      fac2 = 0  ! exclude time-derivative
    else
      fac2 = 1._dp / deltat
    end if

    if ( vector ) then

!     right-hand side

      if ( timeint == 3 ) then

!       first-order, BDF1

        do N = 1, ndfp
          elemvec(N) = - sum ( psi(:,N) * ( &
            fac1 * ( fac2 * ( png_iter - png ) + unp1gradpiter ) + divgiter &
                                          ) * detF * wg )
        end do

      else if ( timeint == 4 ) then

!       second-order, BDF2

        do N = 1, ndfp
          elemvec(N) = - sum ( psi(:,N) * ( &
            fac1 * ( fac2 * ( 3*png_iter/2 - 2*png + pnm1g/2 ) &
                          + unp1gradpiter ) + divgiter &
                                          ) * detF * wg )
        end do

      end if

    end if

    if ( matrix ) then

!     pressure-pressure part

      if ( timeint == 3 ) then

!       first-order, BDF1

        do N = 1, ndfp
          do M = 1, ndfp
            elemmat(N,posp(M)) = fac1 * sum ( psi(:,N) * &
                ( fac2 * psi(:,M) + unp1gradpsi(:,M) ) * detF * wg )
          end do
        end do

      else if ( timeint == 4 ) then

!       second-order, BDF2

        do N = 1, ndfp
          do M = 1, ndfp
            elemmat(N,posp(M)) = fac1 * sum ( psi(:,N) * &
                ( 3*fac2/2 * psi(:,M) + unp1gradpsi(:,M) ) * detF * wg )
          end do
        end do

      end if

!     pressure-velocity part (delta of u.nabla p/K + div u with respect to u)

      do N = 1, ndfp
        do M = 1, ndf
          do i = 1, ndim
            if ( coorsys == 1 .and. i == 2 ) then ! axisymmetric
              work = dphidx(:,M,2) + phi(:,M) / xg(:,2)
            else
              work = dphidx(:,M,i)
            end if
            elemmat( N, pos(M,i) ) = sum ( psi(:,N) * &
                      ( fac1 * phi(:,M) * gradpiter(:,i) + work ) * detF * wg )
          end do
        end do
      end do

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

      allocate ( work(ninti), dpsi(ninti,ndfp,ndim), dpsidx(ninti,ndfp,ndim) )
      allocate ( unp1gradpsi(ninti,ndfp) )
      allocate ( unp1gradpiter(ninti) )
      allocate ( gradpiter(ninti,ndim) )
      allocate ( tmp(ndf,ndim), u(ndim*ndf), uvec_iter(ninti,ndim) )
      allocate ( divgiter(ninti) )
      allocate ( uvecmeshhat(ninti,ndim) )
      allocate ( p(ndfp), png(ninti), pnm1g(ninti), png_iter(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work, dpsi, dpsidx )
      deallocate ( unp1gradpsi )
      deallocate ( unp1gradpiter )
      deallocate ( gradpiter )
      deallocate ( tmp, u, uvec_iter )
      deallocate ( divgiter )
      deallocate ( uvecmeshhat )
      deallocate ( p, png, pnm1g, png_iter )

    end subroutine deallocate_arrays

  end subroutine implicit_pressure_ce_elem


! Element routine for the testing the integral pressure constraint.
! For use with integrate.

  subroutine integrate_pressure_function ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec )

    use compressible_fluid_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: intnr
    real(dp) :: Kmod, p0, V0, p_level

    if ( first ) then

!     first element

!     check coefficients

      call check ( coefficients, 'integrate_pressure_function', ncoefi=350, &
        ncoefr=300, indexarray=[62,308], minimum=[0,1], maximum=[0,5] )

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!     allocate arrays

      allocate ( wg(ninti), detF(ninti) )
      allocate ( work(ninti), xg(ninti,ndim) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), psi(ninti,ndfp) )
      allocate ( x(nodalp,ndim), dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim) )
      allocate ( p(ndfp) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      call set_shape_function ( shapefuncp, xig, psi ) ! pressure

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    xg = matmul ( phi, x(1:ndf,:) )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

!   get pressure

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, p, &
      physq=[physqpress], layer=layer )

    work = matmul ( psi, p )

!   integration of pressure function over the element

    intnr = coefficients%i(308)
    Kmod = coefficients%r(252)
    p0 = coefficients%r(253)
    V0 = coefficients%r(254)

    if ( any( intnr == [1,2] ) ) then
      elemvec = sum ( exp((work-p0)/Kmod) * detF * wg )
      if ( first ) elemvec = elemvec - V0
    else if ( any( intnr == [3,4] ) ) then
      elemvec = sum ( ( exp((work-p0)/Kmod) - 1 ) * detF * wg )
    else if ( intnr == 5 ) then
      p_level = coefficients%r(255)
      elemvec = sum ( ( work - p_level ) * detF * wg )
    end if

    if ( any( intnr == [2,4] ) ) elemvec = Kmod * elemvec

    if ( last ) then

!     last element

      deallocate ( wg, detF )
      deallocate ( work, xg )
      deallocate ( xig, phi, psi )
      deallocate ( x, dphi, F )
      deallocate ( Finv )
      deallocate ( p )

    end if

  end subroutine integrate_pressure_function


! Internal element routine for the contribution of the integral pressure
! constraint to the pressure-pressure part of the Jacobian matrix.

  subroutine pp_constraint_Jacobian_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use compressible_fluid_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: N, M, isv_iter, intnr
    real(dp) :: Kmod, p0

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'pp_constraint_Jacobian_elem', &
        ncoefi=350, ncoefr=300, indexarray=[307], minimum=[1] )
      call check ( coefficients, 'pp_constraint_Jacobian_elem', &
        indexarray=[62,308], minimum=[0,1], maximum=[0,4] )

    end if

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

!   coefficients

    isv_iter = get_coefficient ( coefficients, index=305, default=1 )
    intnr = coefficients%i(308)
    Kmod = coefficients%r(252)
    p0 = coefficients%r(253)

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

!   get iteration pressure

    call get_sysvector ( mesh, problem, oldvectors%s(isv_iter)%p, &
      elgrp, elem, p, physq=[physqpress], layer=layer )

    pg_iter = matmul ( psi, p )

!   get iteration Lagrange multiplier

    call get_sysvector_constraint ( mesh, problem, oldvectors%s(isv_iter)%p, &
      constraint=coefficients%i(307), u=lagmul )

    if ( vector ) then

      write(*,'(/a/a/)') 'Error in pp_constraint_Jacobian_elem: ', &
        ' no vector to build. Call build_system with buildvector=.false.'
      stop

    end if

    if ( matrix ) then

!     pressure-pressure part

      work = lagmul(1) / Kmod ** 2 * exp( (pg_iter-p0)/Kmod )

      do N = 1, ndfp
        do M = 1, ndfp
          elemmat(N,M) = sum ( psi(:,N) * psi(:,M) * work * detF * wg )
        end do
      end do

      if ( any( intnr == [2,4] ) ) elemmat = Kmod * elemmat

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

      allocate ( work(ninti) )
      allocate ( lagmul(1) )
      allocate ( p(ndfp), pg_iter(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work )
      deallocate ( lagmul )
      deallocate ( p, pg_iter )

    end subroutine deallocate_arrays

  end subroutine pp_constraint_Jacobian_elem


! Constraint on an elementset (having current volume V) of
!
!    \int_V exp((p-p0)/K) dV - V0 = 0
!
! where V0 is the reference volume of the volume V.
!
  subroutine constraint_pressure_int_elementset ( mesh, problem, constr, &
    elem, node, matrix, vector, first, last, coefficients, oldvectors, &
    elemmat, elemmat2, elemmatadd, elemvec, elemvecadd )

    use compressible_fluid_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: elgrp, i
    integer :: isv_iter, intnr
    real(dp) :: Kmod, p0, V0, p_level

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'constraint_pressure_int_elementset', &
        ncoefi=350, ncoefr=300, indexarray=[62,308], minimum=[0,1], &
        maximum=[0,5] )

    end if

    elgrp = node ! node is used to transfer element group number

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors, maxvel3D=1 )

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

!   coefficients

    isv_iter = get_coefficient ( coefficients, index=305, default=1 )
    intnr = coefficients%i(308)
    Kmod = coefficients%r(252)
    p0 = coefficients%r(253)
    V0 = coefficients%r(254)

!   get iteration pressure

    call get_sysvector ( mesh, problem, oldvectors%s(isv_iter)%p, &
      elgrp, elem, p, physq=[physqpress], layer=layer )

    pg_iter = matmul ( psi, p )

    if ( intnr /= 5 ) work = exp( (pg_iter-p0)/Kmod )

!   vector

    if ( vector ) then

      if ( any( intnr == [1,2] ) ) then
        elemvec = - sum ( work * detF * wg )
        if ( first ) elemvec = elemvec + V0
      else if ( any( intnr == [3,4] ) ) then
        elemvec = - sum ( ( work - 1 ) * detF * wg )
      else if ( intnr == 5 ) then
        p_level = coefficients%r(255)
        elemvec = - sum ( ( pg_iter - p_level ) * detF * wg )
      end if

      if ( any( intnr == [2,4] ) ) elemvec = Kmod * elemvec

    end if

!   matrix

    if ( matrix ) then

      if ( intnr == 5 ) then
        do i = 1, ndfp
          elemmat(1,i) = sum ( psi(:,i) * detF * wg )
        end do
      else
        do i = 1, ndfp
          elemmat(1,i) = sum ( work * psi(:,i) * detF * wg ) / Kmod
        end do
      end if

      if ( any( intnr == [2,4] ) ) elemmat = Kmod * elemmat

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

      allocate ( work(ninti) )
      allocate ( p(ndfp), pg_iter(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work )
      deallocate ( p, pg_iter )

    end subroutine deallocate_arrays

  end subroutine constraint_pressure_int_elementset

end module compressible_fluid_elements_m
