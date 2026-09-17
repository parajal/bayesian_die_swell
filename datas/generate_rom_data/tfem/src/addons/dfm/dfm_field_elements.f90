
! Copyright (C) 2016-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for solving deformation fields with SUPG stabilization
!

module dfm_field_elements_m

  use tfem_elem_m
  use stokes_set_globals_m
  use devss_set_globals_m

  implicit none

  save

contains


! Internal element routine for the deformation fields using SUPG
! Decoupled fields in tau-Gauss points using Gauss basis functions.
! Legacy method from 2001 paper with DG in space replaced by SUPG
! and solved for F instead of B.

  subroutine dfm_elem1 ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use dfm_globals_m
    use supg_utils_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    logical :: checkzero
    integer :: i, j, ip, m, method, k, n
    real(dp) :: deltat, fac, gamma0, alpha0, alpha1, alpha2, esize

!   check coefficients

    call check ( coefficients, 'dfm_elem1', &
      indexarray=[48,355], minimum=[0,0], maximum=[1,1] )

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    if ( coefficients%i(357) == 0 ) then
!     Projected gradient
      call set_devssg_elem ( mesh, problem, elgrp, elem, first, &
        last, coefficients, oldvectors )
    end if

    call set_dfm_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, gauss_basis=.true. )


!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of the element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      write(*,'(/a/a/)') 'Error in dfm_elem1:', &
                         '  axisymmetric element not yet implemented'
      stop
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( coefficients%i(357) == 1 ) then
      call shape_derivative ( dphi, Finv, dphidx )
    end if

    call shape_derivative ( dtheta, Finv, dthetadx )


!   get velocity vector

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel], layer=layer )

    tmp = reshape ( u, [ndf,ndim] )

    uvec = matmul ( phi, tmp )


!   get velocity gradient vector

    if ( coefficients%i(357) == 0 ) then

!     Projected gradient

      select case ( coefficients%i(356) )
      case(0)
!       DEVSS: in solution vector
        call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
          elgrp, elem, g, physq=[physqgrad], layer=layer )
      case(1)
!       separate vector
        call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(3)%p, &
          elgrp, elem, g, layer=layer )
      case default
        call errormsg_case_default ( 'dfm_elem1', 'coefficients%i(356)', &
          int_value=coefficients%i(356) )
      end select

      work2 = reshape ( g, [ndfg,ncompg] )

      gvec(:,1:ncompg) = matmul ( zeta, work2 )

      do ip = 1, ninti
        gradu(ip,:,:) = &
                       transpose ( reshape ( gvec(ip,1:ncompg), [ndim,ndim] ) )
      end do

    else if ( coefficients%i(357) == 1 ) then

!     direct velocity gradient

!     get velocity vector

      call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
        elgrp, elem, u, physq=[physqvel], layer=layer )

      tmp = reshape ( u, [ndf,ndim] )

      do j = 1, ndim
        gradu(:,:,j) = matmul ( dphidx(:,:,j), tmp )
      end do

      do ip = 1, ninti
        gvec(ip,1:ndim**2) = &
                       reshape ( transpose ( gradu(ip,:,:) ), [ndim*ndim] )
      end do

    end if

    if ( coorsys == 1 ) then
      gvec(:,5) = uvec(:,2) / xg(:,2)
    end if


!   convection operator

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, u )

      tmp = reshape ( u, [ndf,ndim] )

      uvecmesh = matmul ( phi, tmp )

!     (u-ugrid).grad operator

      do ip = 1, ninti
        ugradtheta(ip,:) = &
                  matmul ( dthetadx(ip,:,:), uvec(ip,:) - uvecmesh(ip,:) )
      end do

    else

!     Eulerian frame

      uvecmesh = 0

!     un.grad operator

      do ip = 1, ninti
        ugradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uvec(ip,:) )
      end do

    end if


!   get deformation fields at previous time steps

    do m = 1, nintvaltau
       do n = 1, ninttau

!       deformation fields at time tn

        call get_sysvector ( mesh, problem, oldvectors%s3(1)%p(1,n,m), &
          elgrp, elem, fn(:,1), posu=posf, layer=layer )

        do j = 2, ncompf
          fn(:,j) = oldvectors%s3(1)%p(j,n,m)%u(posf)
        end do

        fng(:,:,n,m) = matmul ( theta, fn )

        if ( timeint >= 2 ) then

!         deformation fields at time tn-1

          do j = 1, ncompf
            fnm1(:,j) = oldvectors%s3(2)%p(j,n,m)%u(posf)
          end do

          fnm1g(:,:,n,m) = matmul ( theta, fnm1 )

        end if

        if ( timeint >= 3 ) then

!         deformation fields at time tn-2

          do j = 1, ncompf
            fnm2(:,j) = oldvectors%s3(3)%p(j,n,m)%u(posf)
          end do

          fnm2g(:,:,n,m) = matmul ( theta, fnm2 )

        end if

      end do
    end do

!   time integration

    select case ( timeint )
    case(1) ! Euler

      fhatg = fng
      gamma0 = 1; alpha0 = 1; alpha1 = 0; alpha2 = 0

    case(2) ! second-order with prediction

      fhatg = 2 * fng - fnm1g
      gamma0 = 1.5_dp; alpha0 = 2; alpha1 = -0.5_dp; alpha2 = 0

    case(3) ! third-order with prediction

      fhatg = 3 * fng - 3 * fnm1g + fnm2g
      gamma0 = 11._dp/6; alpha0 = 3; alpha1 = -1.5_dp; alpha2 = 1._dp/3

    case default

      call errormsg_case_default ( 'dfm_elem1', 'timeint', int_value=timeint )

    end select


!   coefficients

    method = coefficients%i(355)

    if ( method == 1 ) then

!     SUPG

!     compute h/U

      checkzero = coefficients%i(368) == 1

      if ( checkzero ) esize = sum ( detF * wg )

      call supg_hU ( ndim, globalshape, hoverU, htype=2, hlocation=2, &
        Uscaling=3, uvec=uvec(:,1:ndim), Finv=Finv, esize=esize, &
        Uglobal=coefficients%r(306), checkzero=checkzero )

      tauSUPG = coefficients%r(302) * hoverU / 2  ! upwinding parameter

    end if

    deltat = coefficients%r(301)
    fac = 1._dp / deltat


!   element matrix

    if ( matrix ) then

      do i = 1, ndff
        do j = 1, ndff
          elemmat(i,j) = sum ( ( theta(:,i) + tauSUPG * ugradtheta(:,i) ) * &
                 ( fac * gamma0 * theta(:,j) +  ugradtheta(:,j) ) * detF * wg )
        end do
      end do

    end if


!   element vector

    if ( vector ) then

!     inverse mapping of reference interval on real interval

      do j = 1, nintvaltau
        fac1(j) = 2 / ( coefficients%ra(1)%a(j+1) - coefficients%ra(1)%a(j) )
      end do

!     jump between boundary condition for tau=0 (F=I) and first interval

      if ( ndim == 2 ) then
        do ip = 1, ninti
          Fjump(ip,:,1) = [ 1, 0, 0, 1 ]
        end do
      else if ( ndim == 3 ) then
        do ip = 1, ninti
          Fjump(ip,:,1) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ]
        end do
      end if
      do ip = 1, ninti
        Fjump(ip,:,1) = matmul( fhatg(ip,:,:,1), phitaulr(1,:) ) - Fjump(ip,:,1)
      end do

!     jump between intervals at tau_j, j=1,...,N

      do j = 2, nintvaltau
        do ip = 1, ninti
          Fjump(ip,:,j) = matmul( fhatg(ip,:,:,j-1), phitaulr(2,:) )
          Fjump(ip,:,j) = matmul( fhatg(ip,:,:,j), phitaulr(1,:) ) &
                          - Fjump(ip,:,j)
        end do
      end do

!     Coupling between fields (in age tau): -(dF/dtau)^k-[F]_j*phi(tau^j)/w_k

      do k = 1, ninttau
        do j = 1, nintvaltau
          do ip = 1, ninti
            rhsg(ip,:,k,j) = - fac1(j) * &
                   ( matmul( fhatg(ip,:,:,j), dphitau(k,:) ) + &
                           Fjump(ip,:,j) * phitaulr(1,k) / wtaug(k) )
          end do
        end do
      end do

!     L.F using matrix multiplication

      do n = 1, ninttau
        do m = 1, nintvaltau
          do ip = 1, ninti
            fmat = transpose( reshape( fhatg(ip,1:ndim**2,n,m), [ndim,ndim] ) )
            rhsg(ip,1:ndim**2,n,m) = rhsg(ip,1:ndim**2,n,m) + &
                      reshape( transpose( &
                            matmul( gradu(ip,:,:), fmat ) ), [ndim*ndim] )
          end do
        end do
      end do

!     time discretization

      select case ( timeint )

      case(1) ! Euler

        rhsg = rhsg + fac * fng

      case(2) ! second-order with prediction

        rhsg = rhsg + fac * ( alpha0 * fng + alpha1 * fnm1g )

      case(3) ! third-order with prediction

        rhsg = rhsg + fac * ( alpha0 * fng + alpha1 * fnm1g + alpha2 * fnm2g )

      case default

        call errormsg_case_default ( 'dfm_elem1', 'timeint', int_value=timeint )

      end select

!     spatial discretization

      do n = 1, ninttau
        do m = 1, nintvaltau
          do k = 1, ncompf
            do i = 1, ndff
              work5(i,k,n,m) = &
                sum ( ( theta(:,i) + tauSUPG * ugradtheta(:,i) ) * &
                          rhsg(:,k,n,m) * detF * wg )
            end do
          end do
        end do
      end do

      elemvec = reshape ( work5, [ ndff * ncompf * nfields ] )

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    if ( coefficients%i(357) == 0 ) then
!     Projected gradient
      call unset_devssg_elem ( last, coefficients )
    end if

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

      allocate ( u(ndim*ndf) )
      allocate ( gradu(ninti,ndim,ndim) )

      allocate ( posf(ndff) )
      allocate ( uvec(ninti,ndim) )
      allocate ( tmp(ndf,ndim) )

      if ( coefficients%i(357) == 0 ) then
        allocate ( g(ncompg*ndfg), work2(ndfg,ncompg) )
      end if

      allocate ( uvecmesh(ninti,ndim) )

      if ( coorsys == 1 ) then
        allocate ( gvec(ninti,ndim**2+1) )
      else
        allocate ( gvec(ninti,ndim**2) )
      end if

      allocate ( fn(ndff,ncompf), fng(ninti,ncompf,ninttau,nintvaltau) )
      if ( timeint >= 2 ) then
        allocate ( fnm1(ndff,ncompf), fnm1g(ninti,ncompf,ninttau,nintvaltau) )
      end if
      if ( timeint >= 3 ) then
        allocate ( fnm2(ndff,ncompf), fnm2g(ninti,ncompf,ninttau,nintvaltau) )
      end if
      allocate ( fhatg(ninti,ncompf,ninttau,nintvaltau) )
      allocate ( fmat(ndim,ndim), rhsg(ninti,ncompf,ninttau,nintvaltau) )
      allocate ( work5(ndff,ncompf,ninttau,nintvaltau) )
      allocate ( Fjump(ninti,ncompf,nintvaltau) )
      allocate ( fac1(nintvaltau) )

      allocate ( tauSUPG(ninti), hoverU(ninti) )
      allocate ( ugradtheta(ninti,ndff) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( u )
      deallocate ( gradu )

      deallocate ( posf )
      deallocate ( uvec )
      deallocate ( tmp )

      if ( coefficients%i(357) == 0 ) then
        deallocate ( g, work2 )
      end if

      deallocate ( uvecmesh )

      deallocate ( gvec )

      deallocate ( fn, fng )
      if ( timeint >= 2 ) then
        deallocate ( fnm1, fnm1g )
      end if
      if ( timeint >= 3 ) then
        deallocate ( fnm2, fnm2g )
      end if
      deallocate ( fhatg )
      deallocate ( fmat, rhsg )
      deallocate ( work5 )
      deallocate ( Fjump )
      deallocate ( fac1 )

      deallocate ( tauSUPG, hoverU )
      deallocate ( ugradtheta )

    end subroutine deallocate_arrays

  end subroutine dfm_elem1


! Internal element routine for the deformation fields using SUPG.
! Decoupled fields in tau-Gauss points. Legacy method from 2001 paper with
! DG in space replaced by SUPG and solved for F instead of B.
! Coupling of fields on interval level for convection in tau direction.
! Fully implicit tau convection with standard basis functions.
! NOTE: The fields in a tau (age) interval become coupled.
!       The number of degrees in the nodes must be set accordingly.
! NOTE: The first oldvectors%s2(1) must contain the "updated" solution vector,
!       assuming the intervals are solved subsequently starting from interval 1.

  subroutine dfm_elem2 ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use dfm_globals_m
    use supg_utils_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    logical :: checkzero
    integer :: i, j, ip, m, method, k, n, i1, i2, k1, k2
    real(dp) :: deltat, fac, gamma0, alpha0, alpha1, alpha2, esize


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    if ( coefficients%i(357) == 0 ) then
!     Projected gradient
      call set_devssg_elem ( mesh, problem, elgrp, elem, first, &
        last, coefficients, oldvectors )
    end if

    call set_dfm_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

!   check coefficients

    call check ( coefficients, 'dfm_elem2', indexarray=[48,355,362], &
      minimum=[0,0,1], maximum=[1,1,nintvaltau] )


!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of the element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      write(*,'(/a/a/)') 'Error in dfm_elem2:', &
                         '  axisymmetric element not yet implemented'
      stop
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( coefficients%i(357) == 1 ) then
      call shape_derivative ( dphi, Finv, dphidx )
    end if

    call shape_derivative ( dtheta, Finv, dthetadx )


!   get velocity vector

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel], layer=layer )

    tmp = reshape ( u, [ndf,ndim] )

    uvec = matmul ( phi, tmp )


!   get velocity gradient vector

    if ( coefficients%i(357) == 0 ) then

!     Projected gradient

      select case ( coefficients%i(356) )
      case(0)
!       DEVSS: in solution vector
        call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
          elgrp, elem, g, physq=[physqgrad], layer=layer )
      case(1)
!       separate vector
        call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(3)%p, &
          elgrp, elem, g, layer=layer )
      case default
        call errormsg_case_default ( 'dfm_elem2', 'coefficients%i(356)', &
          int_value=coefficients%i(356) )
      end select

      work2 = reshape ( g, [ndfg,ncompg] )

      gvec(:,1:ncompg) = matmul ( zeta, work2 )

      do ip = 1, ninti
        gradu(ip,:,:) = &
                       transpose ( reshape ( gvec(ip,1:ncompg), [ndim,ndim] ) )
      end do

    else if ( coefficients%i(357) == 1 ) then

!     direct velocity gradient

!     get velocity vector

      call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
        elgrp, elem, u, physq=[physqvel], layer=layer )

      tmp = reshape ( u, [ndf,ndim] )

      do j = 1, ndim
        gradu(:,:,j) = matmul ( dphidx(:,:,j), tmp )
      end do

      do ip = 1, ninti
        gvec(ip,1:ndim**2) = &
                       reshape ( transpose ( gradu(ip,:,:) ), [ndim*ndim] )
      end do

    end if

    if ( coorsys == 1 ) then
      gvec(:,5) = uvec(:,2) / xg(:,2)
    end if


!   convection operator

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, u )

      tmp = reshape ( u, [ndf,ndim] )

      uvecmesh = matmul ( phi, tmp )

!     (u-ugrid).grad operator

      do ip = 1, ninti
        ugradtheta(ip,:) = &
                  matmul ( dthetadx(ip,:,:), uvec(ip,:) - uvecmesh(ip,:) )
      end do

    else

!     Eulerian frame

      uvecmesh = 0

!     un.grad operator

      do ip = 1, ninti
        ugradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uvec(ip,:) )
      end do

    end if


!   get deformation fields at previous time steps

    i2 = coefficients%i(362) ! the interval to do implicitly

!   deformation fields at time tn

    call get_sysvector ( mesh, problem, oldvectors%s2(2)%p(1,i2), &
      elgrp, elem, fu, posu=posf, layer=layer )

    fn2(:,1,:) = reshape( fu, [ndff,ninttau] )

    do j = 2, ncompf
      fn2(:,j,:) = reshape ( oldvectors%s2(2)%p(j,i2)%u(posf), [ndff,ninttau] )
    end do

    do n = 1, ninttau
      fn2g(:,:,n) = matmul ( theta, fn2(:,:,n) )
    end do

    if ( timeint >= 2 ) then

!     deformation fields at time tn-1

      do j = 1, ncompf
        fn2m1(:,j,:) = &
               reshape ( oldvectors%s2(3)%p(j,i2)%u(posf), [ndff,ninttau] )
      end do

      do n = 1, ninttau
        fn2m1g(:,:,n) = matmul ( theta, fn2m1(:,:,n) )
      end do

    end if

    if ( timeint >= 3 ) then

!     deformation fields at time tn-2

      do j = 1, ncompf
        fn2m2(:,j,:) = &
               reshape ( oldvectors%s2(4)%p(j,i2)%u(posf), [ndff,ninttau] )
      end do

      do n = 1, ninttau
        fn2m2g(:,:,n) = matmul ( theta, fn2m2(:,:,n) )
      end do

    end if

!   time integration

    select case ( timeint )
    case(1) ! Euler

      fhat2g = fn2g
      gamma0 = 1; alpha0 = 1; alpha1 = 0; alpha2 = 0

    case(2) ! second-order with prediction

      fhat2g = 2 * fn2g - fn2m1g
      gamma0 = 1.5_dp; alpha0 = 2; alpha1 = -0.5_dp; alpha2 = 0

    case(3) ! third-order with prediction

      fhat2g = 3 * fn2g - 3 * fn2m1g + fn2m2g
      gamma0 = 11._dp/6; alpha0 = 3; alpha1 = -1.5_dp; alpha2 = 1._dp/3

    case default

      call errormsg_case_default ( 'dfm_elem2', 'timeint', int_value=timeint )

    end select

    if ( i2 /= 1 ) then

!     get F at tn+1 already computed in interval before this one

      i1 = i2 - 1   ! the interval before i2

!     deformation fields at time tn+1

      do j = 1, ncompf
        fnp1(:,j,:) = &
               reshape ( oldvectors%s2(1)%p(j,i1)%u(posf), [ndff,ninttau] )
      end do

      do n = 1, ninttau
        fnp1g(:,:,n) = matmul ( theta, fnp1(:,:,n) )
      end do

    else

      i1 = i2

    end if

!   coefficients

    method = coefficients%i(355)

    if ( method == 1 ) then

!     SUPG

!     compute h/U

      checkzero = coefficients%i(368) == 1

      if ( checkzero ) esize = sum ( detF * wg )

      call supg_hU ( ndim, globalshape, hoverU, htype=2, hlocation=2, &
        Uscaling=3, uvec=uvec(:,1:ndim), Finv=Finv, esize=esize, &
        Uglobal=coefficients%r(306), checkzero=checkzero )

      tauSUPG = coefficients%r(302) * hoverU / 2  ! upwinding parameter

    end if

    deltat = coefficients%r(301)
    fac = 1._dp / deltat

!   mapping of reference interval on real interval

    do j = i1, i2
      fac1(j) = ( coefficients%ra(1)%a(j+1) - coefficients%ra(1)%a(j) ) / 2
    end do

!   m_ij for the tau part

    do i = 1, ninttau
      do j = 1, ninttau
        mmatrix(i,j) = fac1(i2) * sum ( phitau(:,i) * phitau(:,j) * wtaug )
      end do
    end do

!   element matrix

    if ( matrix ) then

!     n_ij for the tau part

      do i = 1, ninttau
        do j = 1, ninttau
          nmatrix(i,j) = sum ( phitau(:,i) * dphitau(:,j) * wtaug ) &
                             + phitaulr(1,i) * phitaulr(1,j)
        end do
      end do

!     work arrays for the spatial part

      do i = 1, ndff
        do j = 1, ndff
          work11(i,j) = sum ( ( theta(:,i) + tauSUPG * ugradtheta(:,i) ) * &
                 ( fac * gamma0 * theta(:,j) +  ugradtheta(:,j) ) * detF * wg )
          work10(i,j) = sum ( ( theta(:,i) + tauSUPG * ugradtheta(:,i) ) &
                                 * theta(:,j) * detF * wg )
        end do
      end do

      do i = 1, ninttau
        do j = 1, ninttau
          k1 = (i-1)*ndff
          k2 = (j-1)*ndff
          do k = 1, ndff
            do m = 1, ndff
              elemmat(k1+k,k2+m) = mmatrix(i,j) * work11(k,m)  &
                                        + nmatrix(i,j) * work10(k,m)
            end do
          end do
        end do
      end do

    end if


!   element vector

    if ( vector ) then

      if ( i2 == 1 ) then

!       jump between boundary condition for tau=0 (F=I) and first interval

        if ( ndim == 2 ) then
          do ip = 1, ninti
            Fjump2(ip,:) = [ 1, 0, 0, 1 ]
          end do
        else if ( ndim == 3 ) then
          do ip = 1, ninti
            Fjump2(ip,:) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ]
          end do
        end if
        Fjump2 = - Fjump2

      else

!       jump between intervals at tau_j, j=1,...,N

        do ip = 1, ninti
          Fjump2(ip,:) = - matmul( fnp1g(ip,:,:), phitaulr(2,:) )
        end do

      end if

!     Coupling between fields (in age tau): phi(tau^j) * F_(j-1)

      do i = 1, ninttau
        do ip = 1, ninti
          rhs2g(ip,:,i) = - phitaulr(1,i) * Fjump2(ip,:)
        end do
      end do

!     L.F using matrix multiplication

      do n = 1, ninttau
        do ip = 1, ninti
          fmat = transpose( reshape( fhat2g(ip,1:ndim**2,n), [ndim,ndim] ) )
          rhs1g(ip,1:ndim**2,n) = reshape( transpose( &
                          matmul( gradu(ip,:,:), fmat ) ), [ndim*ndim] )
!Btensor Use the following lines instead of the line above for B instead of F,
!Btensor i.e. use L.B + B.L^T for the rhs.
!Btensor               matmul( gradu(ip,:,:), fmat ) + &
!Btensor               matmul( fmat, transpose(gradu(ip,:,:)) ) ), [ndim*ndim] )
        end do
      end do

!     time discretization

      select case ( timeint )

      case(1) ! Euler

        rhs1g = rhs1g + fac * fn2g

      case(2) ! second-order with prediction

        rhs1g = rhs1g + fac * ( alpha0 * fn2g + alpha1 * fn2m1g )

      case(3) ! third-order with prediction

        rhs1g = rhs1g + &
                 fac * ( alpha0 * fn2g + alpha1 * fn2m1g + alpha2 * fn2m2g )

      case default

        call errormsg_case_default ( 'dfm_elem2', 'timeint', int_value=timeint )

      end select

      do ip = 1, ninti
        rhs1g(ip,:,:) = matmul ( rhs1g(ip,:,:), mmatrix )
      end do

!     spatial discretization

      do k = 1, ncompf
        do n = 1, ninttau
          do i = 1, ndff
            work6(i,n,k) = &
              sum ( ( theta(:,i) + tauSUPG * ugradtheta(:,i) ) * &
                    ( rhs1g(:,k,n) + rhs2g(:,k,n) ) * detF * wg )
          end do
        end do
      end do

      elemvec = reshape ( work6, [ ndff * ninttau * ncompf ] )

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    if ( coefficients%i(357) == 0 ) then
!     Projected gradient
      call unset_devssg_elem ( last, coefficients )
    end if

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

      allocate ( u(ndim*ndf) )
      allocate ( gradu(ninti,ndim,ndim) )

      allocate ( posf(ndff*ninttau) )
      allocate ( uvec(ninti,ndim) )
      allocate ( tmp(ndf,ndim) )

      if ( coefficients%i(357) == 0 ) then
        allocate ( g(ncompg*ndfg), work2(ndfg,ncompg) )
      end if

      allocate ( uvecmesh(ninti,ndim) )

      if ( coorsys == 1 ) then
        allocate ( gvec(ninti,ndim**2+1) )
      else
        allocate ( gvec(ninti,ndim**2) )
      end if

      allocate ( fu(ndff*ninttau) )
      allocate ( fn2(ndff,ncompf,ninttau), &
                 fn2g(ninti,ncompf,ninttau) )
      if ( timeint >= 2 ) then
        allocate ( fn2m1(ndff,ncompf,ninttau) )
        allocate ( fn2m1g(ninti,ncompf,ninttau) )
      end if
      if ( timeint >= 3 ) then
        allocate ( fn2m2(ndff,ncompf,ninttau) )
        allocate ( fn2m2g(ninti,ncompf,ninttau) )
      end if
      allocate ( fnp1(ndff,ncompf,ninttau) )
      allocate ( fhat2g(ninti,ncompf,ninttau) )
      allocate ( fnp1g(ninti,ncompf,ninttau) )
      allocate ( fmat(ndim,ndim), rhs1g(ninti,ncompf,ninttau) )
      allocate ( rhs2g(ninti,ncompf,ninttau) )
      allocate ( work6(ndff,ninttau,ncompf) )
      allocate ( work10(ndff,ndff), work11(ndff,ndff) )
      allocate ( Fjump2(ninti,ncompf) )
      allocate ( fac1(nintvaltau) )
      allocate ( mmatrix(ninttau,ninttau), nmatrix(ninttau,ninttau) )

      allocate ( tauSUPG(ninti), hoverU(ninti) )
      allocate ( ugradtheta(ninti,ndff) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( u )
      deallocate ( gradu )

      deallocate ( posf )
      deallocate ( uvec )
      deallocate ( tmp )

      if ( coefficients%i(357) == 0 ) then
        deallocate ( g, work2 )
      end if

      deallocate ( uvecmesh )

      deallocate ( gvec )

      deallocate ( fu )
      deallocate ( fn2, fn2g )
      if ( timeint >= 2 ) then
        deallocate ( fn2m1, fn2m1g )
      end if
      if ( timeint >= 3 ) then
        deallocate ( fn2m2, fn2m2g )
      end if
      deallocate ( fnp1 )
      deallocate ( fhat2g )
      deallocate ( fnp1g )
      deallocate ( fmat, rhs1g )
      deallocate ( rhs2g )
      deallocate ( work6 )
      deallocate ( work10, work11 )
      deallocate ( Fjump2 )
      deallocate ( fac1 )
      deallocate ( mmatrix, nmatrix )

      deallocate ( tauSUPG, hoverU )
      deallocate ( ugradtheta )

    end subroutine deallocate_arrays

  end subroutine dfm_elem2


! compute deformation field component in all nodes

  subroutine dfm_deriv_field ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    use dfm_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: comp, intval, field, loca


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_dfm ( coefficients )

!     check coefficients

      call check ( coefficients, 'dfm_deriv_field', &
        indexarray=[358,359,360,361], minimum=[1,1,1,0], &
        maximum=[ncompf,nintvaltau,ninttau,1] )

!     allocate arrays

      allocate ( theta(nodalp,ndff) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( u(ndff) )
      allocate ( fu(ndff*ninttau) )
      allocate ( work2(ndff,ninttau) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function ( shapefuncf, xrnod, theta )

    end if

    comp = coefficients%i(358)
    intval = coefficients%i(359)
    field = coefficients%i(360)
    loca = coefficients%i(361)

    if ( loca == 0 ) then

!     coupled fields

      call get_sysvector ( mesh, problem, oldvectors%s2(1)%p(comp,intval), &
        elgrp, elem, fu, layer=layer )

      work2 = reshape( fu, [ndff,ninttau] )

      elemvec = matmul ( theta, work2(:,field) )

    else if ( loca == 1 ) then

!     decoupled fields (legacy dfm)

      call get_sysvector ( mesh, problem, &
        oldvectors%s3(1)%p(comp,field,intval), elgrp, elem, u, layer=layer )

      elemvec = matmul ( theta, u )

    end if

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( u )
      deallocate ( fu )
      deallocate ( work2 )

    end if

  end subroutine dfm_deriv_field


! compute deformation field tensor in all nodes

  subroutine dfm_deriv_field_tensor ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    use dfm_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: j, intval, field, loca


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_dfm ( coefficients )

!     check coefficients

      call check ( coefficients, 'dfm_deriv_field_tensor', &
        indexarray=[359,360,361], minimum=[1,1,0], &
        maximum=[nintvaltau,ninttau,1] )

!     allocate arrays

      allocate ( theta(nodalp,ndff) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( posf(ndff*ninttau) )
      allocate ( fn(ndff,ncompf) )
      allocate ( fu(ndff*ninttau), fn2(ndff,ncompf,ninttau) )
      allocate ( work2(nodalp,ncompf) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function ( shapefuncf, xrnod, theta )

    end if

    intval = coefficients%i(359)
    field = coefficients%i(360)
    loca = coefficients%i(361)

    if ( loca == 0 ) then

!     coupled fields

      call get_sysvector ( mesh, problem, oldvectors%s2(1)%p(1,intval), &
        elgrp, elem, fu, posu=posf, layer=layer )

      fn2(:,1,:) = reshape( fu, [ndff,ninttau] )

      do j = 2, ncompf
        fn2(:,j,:) = &
              reshape( oldvectors%s2(2)%p(j,intval)%u(posf), [ndff,ninttau] )
      end do

      work2 = matmul( theta, fn2(:,:,field) )

    else if ( loca == 1 ) then

!     decoupled fields (legacy dfm)

      call get_sysvector ( mesh, problem, oldvectors%s3(1)%p(1,field,intval), &
        elgrp, elem, fn(:,1), posu=posf(1:ndff), layer=layer )

      do j = 2, ncompf
        fn(:,j) = oldvectors%s3(1)%p(j,field,intval)%u(posf(1:ndff))
      end do

      work2 = matmul( theta, fn )

    end if

    elemvec = reshape ( work2, [nodalp*ncompf] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( posf )
      deallocate ( fn )
      deallocate ( fu, fn2 )
      deallocate ( work2 )

    end if

  end subroutine dfm_deriv_field_tensor


! compute deformation field derivatives in all nodes

  subroutine dfm_deriv ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    use dfm_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: i, j, intval, field, loca


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_dfm ( coefficients )

!     check coefficients

      call check ( coefficients, 'dfm_deriv', &
        indexarray=[358,359,360,361], minimum=[1,1,1,0], &
        maximum=[2,nintvaltau,ninttau,1] )

!     allocate arrays

      allocate ( theta(nodalp,ndff) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( posf(ndff*ninttau) )
      allocate ( fn(ndff,ncompf) )
      allocate ( fu(ndff*ninttau), fn2(ndff,ncompf,ninttau) )
      allocate ( work6(nodalp,ndim,ndim) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function ( shapefuncf, xrnod, theta )

    end if

    intval = coefficients%i(359)
    field = coefficients%i(360)
    loca = coefficients%i(361)

    if ( loca == 0 ) then

!     coupled fields

      call get_sysvector ( mesh, problem, oldvectors%s2(1)%p(1,intval), &
        elgrp, elem, fu, posu=posf, layer=layer )

      fn2(:,1,:) = reshape( fu, [ndff,ninttau] )

      do j = 2, ncompf
        fn2(:,j,:) = &
              reshape( oldvectors%s2(1)%p(j,intval)%u(posf), [ndff,ninttau] )
      end do

      work6 = reshape( matmul( theta, fn2(:,1:ndim**2,field) ), &
                                                     [nodalp,ndim,ndim] )

    else if ( loca == 1 ) then

!     decoupled fields (legacy dfm)

      call get_sysvector ( mesh, problem, oldvectors%s3(1)%p(1,field,intval), &
        elgrp, elem, fn(:,1), posu=posf(1:ndff), layer=layer )

      do j = 2, ncompf
        fn(:,j) = oldvectors%s3(1)%p(j,field,intval)%u(posf(1:ndff))
      end do

      work6 = reshape( matmul( theta, fn(:,1:ndim**2) ), [nodalp,ndim,ndim] )

    end if

    do i = 1, nodalp
      work6(i,:,:) = transpose ( work6(i,:,:) )
    end do

    select case ( coefficients%i(358) )
    case (1)

!     determinant of F

      if ( ndim == 2 ) then

        elemvec =  work6(:,1,1)*work6(:,2,2) - work6(:,2,1)*work6(:,1,2)

      else if ( ndim == 3 ) then

        elemvec =  work6(:,1,1)*work6(:,2,2)*work6(:,3,3) &
                  -work6(:,1,1)*work6(:,2,3)*work6(:,3,2) &
                  -work6(:,2,1)*work6(:,1,2)*work6(:,3,3) &
                  +work6(:,2,1)*work6(:,1,3)*work6(:,3,2) &
                  +work6(:,3,1)*work6(:,1,2)*work6(:,2,3) &
                  -work6(:,3,1)*work6(:,1,3)*work6(:,2,2)

      end if

    case (2)

!     B=F.F^T

      do i = 1, nodalp
        work6(i,:,:) = matmul ( work6(i,:,:), transpose(work6(i,:,:) ) ) ! F.F^T
        work6(i,:,:) = transpose ( work6(i,:,:) ) ! correct storage in elemvec
      end do

      elemvec = reshape ( work6, [ nodalp*ndim**2 ] )

    case default

      call errormsg_case_default ( 'dfm_elem2', 'coefficients%i(358)', &
        int_value=coefficients%i(358) )

    end select

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( posf )
      deallocate ( fn )
      deallocate ( fu, fn2 )
      deallocate ( work6 )

    end if

  end subroutine dfm_deriv


! sample value of deformation fields in one node of the object

  subroutine dfm_sample ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, uvect )

    use dfm_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: uvect


    integer :: elgrp, elem, intval, field, loca, j
    real(dp), dimension(:,:), allocatable :: xr


    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

!   set globals

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )
    call set_globals_dfm ( coefficients )

!   allocate arrays

    allocate ( theta(1,ndff) )
    allocate ( xr(1,ndim) )
    allocate ( posf(ndff*ninttau) )
    allocate ( fn(ndff,ncompf) )
    allocate ( fu(ndff*ninttau), fn2(ndff,ncompf,ninttau) )
    allocate ( work6(ncompf,ninttau,nintvaltau) )

    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

!   set shape function in the point

    call set_shape_function ( shapefuncf, xr, theta )

    loca = coefficients%i(361)

    if ( loca == 0 ) then

!     coupled fields

      do intval = 1, nintvaltau

        call get_sysvector ( mesh, problem, oldvectors%s2(1)%p(1,intval), &
          elgrp, elem, fu, posu=posf, layer=layer )

        fn2(:,1,:) = reshape( fu, [ndff,ninttau] )

        do j = 2, ncompf
          fn2(:,j,:) = &
                reshape( oldvectors%s2(1)%p(j,intval)%u(posf), [ndff,ninttau] )
        end do

        do j = 1, ncompf
          work6(j,:,intval) = matmul( theta(1,:), fn2(:,j,:) )
        end do

      end do

    else if ( loca == 1 ) then

!     decoupled fields (legacy dfm)

      do intval = 1, nintvaltau
        do field = 1, ninttau

          call get_sysvector ( mesh, problem, &
            oldvectors%s3(1)%p(1,field,intval), elgrp, elem, fn(:,1), &
            posu=posf(1:ndff), layer=layer )

          do j = 2, ncompf
            fn(:,j) = oldvectors%s3(1)%p(j,field,intval)%u(posf(1:ndff))
          end do

          work6(:,field,intval) = matmul( theta(1,:), fn )

        end do
      end do

    end if

    uvect = reshape ( work6, [ ncompf*ninttau*nintvaltau ] )

    deallocate ( theta )
    deallocate ( xr )
    deallocate ( posf )
    deallocate ( fn )
    deallocate ( fu, fn2 )
    deallocate ( work6 )

  end subroutine dfm_sample


! Age (tau) nodal points array

  subroutine dfm_tau_nodal_points ( coefficients, tau_nodal_points )

    use dfm_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   array of size ninttau x nintvaltau, that give the age points of the
!   "nodal points" in tau direction
    real(dp), dimension(:,:), intent(out) :: tau_nodal_points


    integer :: m


!   deformation fields

    call check ( coefficients, 'dfm_tau_nodal_points', &
      indexarray=[352], minimum=[1], maximum=[2] )

    intpoltau = coefficients%i(352)
    nintvaltau = coefficients%i(353)
    ninttau = intpoltau + 1

!   allocate arrays

    allocate ( phitaunod(ninttau,2) )

!   Shape functions

    select case ( intpoltau )
    case (1)  ! P1
      call shape_line_P1 ( [-1._dp,1._dp], phitaunod )
    case (2)  ! P2
      call shape_line_P1 ( [-1._dp,0._dp,1._dp], phitaunod )
    case default
      call errormsg_case_default ( 'dfm_tau_nodal_points', &
        'intpoltau', int_value=intpoltau )
    end select

!   compute tau points

    do m = 1, nintvaltau
      tau_nodal_points(:,m) = matmul(phitaunod,coefficients%ra(1)%a(m:m+1))
    end do

!   deallocate arrays

    deallocate ( phitaunod )

  end subroutine dfm_tau_nodal_points


! Interpolate deformation fields from Gauss points to equidistant (nodal) points
! Only for decoupled fields in tau-Gauss points.

  subroutine dfm_interpolate_to_tau_nodes ( coefficients, fval, &
    tau_nodal_points, sol_f, sol_f1 )

    use dfm_globals_m
    use shapefunc_gauss_m

    type(coefficients_t), intent(in) :: coefficients

!   the subscript giving the degrees of freedom that need to be interpolated.
    type(subscript_t), intent(in) :: fval

!   array of size ninttau x nintvaltau, that give the age points of the
!   interpolated fields
    real(dp), dimension(:,:), intent(out) :: tau_nodal_points

!   sol_f: input. If sol_f1 is present output is written to sol_f1, otherwise
!   sol_f is overwritten with the new data.
    type(sysvector_t), dimension(:,:,:), intent(inout) :: sol_f
    type(sysvector_t), dimension(:,:,:), intent(inout), optional :: sol_f1


    integer :: j, comp, k, i, m


    call set_globals_dfm ( coefficients )

!   allocate arrays

    allocate ( phitau(ninttau,ninttau), phitaunod(ninttau,2), work(ninttau) )

!   Shape functions

    select case ( intpoltau )
    case (1)  ! P1
      call shape_line_gauss_P1 ( [-1._dp,1._dp], phitau )
      call shape_line_P1 ( [-1._dp,1._dp], phitaunod )
    case (2)  ! P2
      call shape_line_gauss_P2 ( [-1._dp,0._dp,1._dp], phitau )
      call shape_line_P1 ( [-1._dp,0._dp,1._dp], phitaunod )
    case default
      call errormsg_case_default ( 'dfm_interpolate_to_tau_nodes', &
        'intpoltau', int_value=intpoltau )
    end select

!   copy

    if ( present(sol_f1) ) call copy ( sol_f, sol_f1 )

!   interpolate and compute tau points

    do m = 1, nintvaltau
      do comp = 1, ncompf
        do k = 1, size(fval%s)
          do i = 1, ninttau
            work(i) = 0
            do j = 1, ninttau
              work(i) = work(i) + phitau(i,j)*sol_f(comp,j,m)%u(fval%s(k))
            end do
          end do
          if ( present(sol_f1) ) then
            do i = 1, ninttau
              sol_f1(comp,i,m)%u(fval%s(k)) = work(i)
            end do
          else
            do i = 1, ninttau
              sol_f(comp,i,m)%u(fval%s(k)) = work(i)
            end do
          end if
        end do
      end do
      tau_nodal_points(:,m) = matmul(phitaunod,coefficients%ra(1)%a(m:m+1))
    end do

!   deallocate arrays

    deallocate ( phitau, phitaunod )

  end subroutine dfm_interpolate_to_tau_nodes


! set global parameters deformation fields (internal element)

  subroutine set_globals_dfm ( coefficients )

    use dfm_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   deformation fields

    call check ( coefficients, 'set_globals_dfm', ncoefi=400, ncoefr=350, &
      indexarray=[354,356,357,352,368], minimum=[0,0,0,1,0], &
      maximum=[3,1,1,2,1] )

    ncompf = ndim ** 2

    intpoltau = coefficients%i(352)
    nintvaltau = coefficients%i(353)
    ninttau = intpoltau + 1
    nfields = nintvaltau * ninttau

!   set number of degrees of freedom for F

    intpolf = coefficients%i(351)

    shapefuncf%globalshape = globalshape
    shapefuncf%interpolation = intpolf
    shapefuncf%numbering = 'regular'

    call set_ndf ( shapefuncf, 'set_globals_dfm', ndf=ndff )

    timeint = get_coefficient ( coefficients, index=354, default=intpoltau+1 )

  end subroutine set_globals_dfm


! Preamble for the deformation field element

  subroutine set_dfm_elem ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, gauss_basis )

    use dfm_globals_m
    use shapefunc_gauss_m
    use set_optional_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
!   Gauss basis functions. default=.false.
    logical, intent(in), optional :: gauss_basis

    logical:: lgauss_basis


    lgauss_basis = set_optional ( variable=gauss_basis, default=.false. )


!   spatial discretization

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

!     set globals

      call set_globals_dfm ( coefficients )

!     allocate arrays

      allocate ( theta(ninti,ndff), dtheta(ninti,ndff,ndim) )
      allocate ( dthetadx(ninti,ndff,ndim) )

!     set shape function F tensor
      call set_shape_function ( shapefuncf, xig, theta, dtheta )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      if ( first ) then

!       first element in this group

!       set globals

        call set_globals_dfm ( coefficients )

      end if

!     allocate arrays

      allocate ( theta(ninti,ndff), dtheta(ninti,ndff,ndim) )
      allocate ( dthetadx(ninti,ndff,ndim) )

!     set shape function F tensor
      call set_shape_function ( shapefuncf, xig, theta, dtheta )

    end if

!   time discretization

    if ( first ) then

!     first element in this group

!     allocate arrays

      allocate ( phitau(ninttau,ninttau), dphitau(ninttau,ninttau) )
      allocate ( phitaulr(2,ninttau) )
      allocate ( taug(ninttau), wtaug(ninttau) )

!     Gauss points

      call Gauss_Legendre_line ( ninttau, taug, wtaug )

      if ( lgauss_basis ) then

!       Shape functions with nodes in the Gauss points

        select case ( intpoltau )
        case (1)  ! P1
          call shape_line_gauss_P1 ( taug, phitau, dphitau )
          call shape_line_gauss_P1 ( [-1._dp,1._dp], phitaulr )
        case (2)  ! P2
          call shape_line_gauss_P2 ( taug, phitau, dphitau )
          call shape_line_gauss_P2 ( [-1._dp,1._dp], phitaulr )
        case default
          call errormsg_case_default ( 'set_dfm_elem', 'intpoltau', &
            int_value=intpoltau )
        end select

      else

!       Shape functions with regular nodes

        select case ( intpoltau )
        case (1)  ! P1
          call shape_line_P1 ( taug, phitau, dphitau )
          call shape_line_P1 ( [-1._dp,1._dp], phitaulr )
        case (2)  ! P2
          call shape_line_P2 ( taug, phitau, dphitau )
          call shape_line_P2 ( [-1._dp,1._dp], phitaulr )
        case default
          call errormsg_case_default ( 'set_dfm_elem', 'intpoltau', &
            int_value=intpoltau )
        end select

      end if

    end if

  end subroutine set_dfm_elem


! Unset the preamble for the deformation fields element
! (deallocate arrays allocated in set_... )

  subroutine unset_dfm_elem ( last, coefficients )

    use dfm_globals_m

    logical, intent(in) :: last
    type(coefficients_t), intent(in) :: coefficients

!   deallocate memory

!   spatial discretization

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      deallocate ( theta, dtheta, dthetadx )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      if ( last ) then

!       last element in this group

      end if

      deallocate ( theta, dtheta, dthetadx )

    end if

!   time discretization

    if ( last ) then

!     last element in this group

      deallocate ( phitau, dphitau )
      deallocate ( phitaulr )
      deallocate ( taug, wtaug )

    end if

  end subroutine unset_dfm_elem

end module dfm_field_elements_m

