! Copyright (C) 2009-2009 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for three-phase diffuse-interface equation

module diffuse_interface_tp_coupled_elements_generic_m

  use tfem_elem_m
  use stokes_set_globals_m

  implicit none

contains


! Internal element routine for the diffuse-interface equation
! The solution vector consists of (c1, c2, mu1, mu2).
!
!          dc1dt + u.nabla c1 - div( M(c1,c2) nabla c1 )  = 0
!          dc2dt + u.nabla c2 - div( M(c1,c2) nabla c2 )  = 0
!    (...) * c1 + (...) * c2 + kappa nabla^2 c1 -kappa/2 nabla^2 c2 + mu1 = 0
!    (...) * c1 + (...) * c2 + kappa nabla^2 c2 -kappa/2 nabla^2 c1 + mu2 = 0
!
!   where (...) depends on the choice for the ternary free energy.
!
!   standard ternary model
!         f(c1,c2) = 1/4[alpha * c1 * c2 + beta * c1^2 * c2^2 * c3 ^2 -
!                        gamm * c1 * c2 * c3 ]
!
!         c3 = 1 - c1 - c2
!         alpha, beta and gamm are passed with coefficients
!
!  The mobility can be taken constant using coefficients%i(151) = 0
!  else it is a function of composition.

  subroutine diffuse_interface_tp_coupled_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip, i1, i2, i3, i4, non_lin_diffusion
    real(dp) :: alpha, beta, gamm, Mcoef, kappa
    real(dp) :: deltat


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


!   some tests

    if ( first ) then

      call check ( coefficients, 'diffuse_interface_tp_coupled_elem', &
        ncoefi=200, ncoefr=150 )

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!
!   oldvectors_di%s(1) = sol
!   oldvectors_di%s(2) = soldin
!   oldvectors_di%s(3) = soldi
!   oldvectors_di%p(1) = problem
!   oldvectors_di%p(2) = problemdi

!   get velocity vector

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, elgrp, &
      elem,u, physq=[physqvel], layer=layer )

    tmp = reshape ( u, [ndf,ndim] )

    u_n = matmul ( phi, tmp )

!   un.grad operator

    do ip = 1, ninti
      ungradphi(ip,:) = matmul ( dphidx(ip,:,:), u_n(ip,:) )
    end do


!   get c_1n variable (concentration) at previous time step for first phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, c_n, physq=[1] , layer=layer )

    c_ng = matmul ( phi, c_n )

!   get c_1iter variable (concentration) at previous iteration for first phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(3)%p, elgrp, &
      elem, citer, physq=[1] , layer=layer )

    citer_g = matmul ( phi, citer )

!   get c_2 variable (concentration) at previous time step for the second phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, c2_n, physq=[2] , layer=layer )

    c2_ng = matmul ( phi, c2_n )

!   get c_2 variable (concentration) at previous iteration for second phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(3)%p, elgrp, &
      elem, c2iter, physq=[2] , layer=layer )

    c2iter_g = matmul ( phi, c2iter )


!   coefficients

    deltat = coefficients%r(101)
    Mcoef  = coefficients%r(102)
    alpha  = coefficients%r(103)
    beta   = coefficients%r(104)
    kappa  = coefficients%r(105)
    gamm   = coefficients%r(107)

    non_lin_diffusion = coefficients%i(151)

!   pointers in unknown vector

    i1 = ndf
    i2 = 2*ndf
    i3 = 3*ndf
    i4 = 4*ndf

    if ( matrix ) then

!     diffusion matrix

      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,i,:) * dphidx(ip,j,:) )
          end do
          S(i,j) = sum ( work * detF * wg )
          S(j,i) = S(i,j)
        end do
      end do

      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = ( citer_g(ip) * ( 1.0_dp - citer_g(ip)) )  * &
                             sum( dphidx(ip,i,:) * dphidx(ip,j,:) )
          end do
          S1(i,j) = sum ( work * detF * wg )
          S1(j,i) = S1(i,j)
        end do
      end do


      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) =  citer_g(ip) * c2iter_g(ip) * &
                                      sum ( dphidx(ip,i,:) * dphidx(ip,j,:) )
          end do
          S2(i,j) = sum ( work * detF * wg )
          S2(j,i) = S2(i,j)
        end do
      end do

      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) =  ( c2iter_g(ip) * ( 1.0_dp- c2iter_g(ip)) ) * &
                               sum ( dphidx(ip,i,:) * dphidx(ip,j,:) )
          end do
          S3(i,j) = sum ( work * detF * wg )
          S3(j,i) = S3(i,j)
        end do
      end do

!     mass matrix

      do i = 1, ndf
        do j = i, ndf
          Mmat(i,j) = sum ( phi(:,i) * phi(:,j) * detF * wg )
          Mmat(j,i) = Mmat(i,j)
        end do
      end do

!    (\partial f  / \partial c_1) = work * c_1  in chemical potential

      work =  - 0.5_dp * alpha * c2iter_g**2 - 0.5_dp * beta &
              + 1.5_dp * citer_g * beta - citer_g**2 * beta  &
              + c2iter_g * beta -1.5_dp * citer_g * c2iter_g * beta &
              - c2iter_g**2 * beta - 0.5_dp * c2iter_g * gamm

      do i = 1, ndf
        do j = i, ndf
          LLbm(i,j) = sum ( phi(:,i) * work * phi(:,j) * detF * wg )
          LLbm(j,i) = LLbm(i,j)
        end do
      end do

!    (\partial f  / \partial c_1) = work * c_2  in chemical potential

      work =   0.5_dp * c2iter_g * beta - 0.5_dp * c2iter_g**2 * beta &
             + 0.25_dp * gamm  -0.25_dp * c2iter_g * gamm


      do i = 1, ndf
        do j = i, ndf
          LLbm1(i,j) = sum ( phi(:,i) * work * phi(:,j) * detF * wg )
          LLbm1(j,i) = LLbm1(i,j)
        end do
      end do

!    (\partial f  / \partial c_2) = work * c_1 in chemical potential

      work =  - 0.5_dp * citer_g * c2iter_g * alpha &
              + 0.5_dp * citer_g * beta - 0.5_dp * citer_g**2 * beta  &
              + c2iter_g * beta - citer_g * c2iter_g * beta  &
              - 1.5_dp * c2iter_g**2 * beta + 0.25_dp * gamm &
               - 0.25_dp * citer_g * gamm - 0.5_dp * c2iter_g * gamm

      do i = 1, ndf
        do j = i, ndf
          LLbm2(i,j) = sum ( phi(:,i) * work * phi(:,j) * detF * wg )
          LLbm2(j,i) = LLbm2(i,j)
        end do
      end do

!    (\partial f  / \partial c_2) = work * c_2 in chemical potential

      work =    -0.5_dp * beta  + 1.5_dp * c2iter_g * beta - c2iter_g**2 * beta

      do i = 1, ndf
        do j = i, ndf
          LLbm3(i,j) = sum ( phi(:,i) * work * phi(:,j) * detF * wg )
          LLbm3(j,i) = LLbm3(i,j)
        end do
      end do

!     advection term

      do i = 1, ndf
        do j = 1, ndf
          LUbm(i,j) = sum ( phi(:,i) * ungradphi(:,j) * detF * wg )
        end do
      end do


!     Fill element matrix
!
!     S : nabla^2
!     S1: nabla  . c1 * (1 - c1) nabla
!     S2: nabla  . c1 * c2 nabla
!     S3: nabla  . c2 * (1 - c2) nabla

! Mauri, Molin, Anderson:

      elemmat(    1:i1,    1:i1 ) = LUbm + Mmat / deltat
      elemmat(    1:i1, i1+1:i2 ) = 0.0_dp
      elemmat(    1:i1, i2+1:i3 ) = Mcoef * S1
      elemmat(    1:i1, i3+1:i4 ) = Mcoef * S2

      elemmat( i1+1:i2,    1:i1 ) = 0.0_dp
      elemmat( i1+1:i2, i1+1:i2 ) = LUbm + Mmat / deltat
      elemmat( i1+1:i2, i2+1:i3 ) = Mcoef * S2
      elemmat( i1+1:i2, i3+1:i4 ) = Mcoef * S3

      elemmat( i2+1:i3,    1:i1 ) = LLbm  - kappa * S
      elemmat( i2+1:i3, i1+1:i2 ) = LLbm1 - (kappa / 2.0_dp) * S
      elemmat( i2+1:i3, i2+1:i3 ) = Mmat
      elemmat( i2+1:i3, i3+1:i4 ) = 0.0_dp

      elemmat( i3+1:i4,    1:i1 ) = LLbm2 -(kappa / 2.0_dp) * S
      elemmat( i3+1:i4, i1+1:i2 ) = LLbm3 - kappa * S
      elemmat( i3+1:i4, i2+1:i3 ) = 0.0_dp
      elemmat( i3+1:i4, i3+1:i4 ) = Mmat

! Lowengrub:

      if (non_lin_diffusion == 0) then
        elemmat(    1:i1, i2+1:i3 ) = Mcoef * S
        elemmat(    1:i1, i3+1:i4 ) = 0.0_dp

        elemmat( i1+1:i2, i2+1:i3 ) = 0.0_dp
        elemmat( i1+1:i2, i3+1:i4 ) = Mcoef * S
      end if

    end if

    if ( vector ) then

!     right-hand side of the c-equation

      do i = 1, ndf
        work(i) = sum ( phi(:,i) * c_ng * detF * wg ) / deltat
      end do

      do i = 1, ndf
        work1(i) = sum ( phi(:,i) * c2_ng * detF * wg ) / deltat
      end do

!     Fill element vector

      elemvec(1:i1)    = work  ! c1
      elemvec(i1+1:i2) = work1 ! c2
      elemvec(i2+1:i3) = 0
      elemvec(i3+1:i4) = 0

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

      allocate ( ungradphi(ninti,ndf) )
      allocate ( u(ndim*ndf), c_n(ndf), c_ng(ninti) )
      allocate ( c2_n(ndf), c2_ng(ninti) )
      allocate ( citer(ndf), citer_g(ninti) )
      allocate ( c2iter(ndf), c2iter_g(ninti) )
      allocate ( work(ninti), work1(ndf) )
      allocate ( u_n(ninti,ndim), tmp(ndf,ndim) )
      allocate ( Mmat(ndf,ndf), S(ndf,ndf), LLbm(ndf,ndf), LUbm(ndf,ndf) )
      allocate ( S1(ndf,ndf), S2(ndf,ndf), S3(ndf,ndf) )
      allocate ( LLbm1(ndf,ndf), LLbm2(ndf,ndf), LLbm3(ndf,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( ungradphi )
      deallocate ( u, c_n, c_ng )
      deallocate ( c2_n, c2_ng )
      deallocate ( citer, citer_g )
      deallocate ( c2iter, c2iter_g )
      deallocate ( work, work1 )
      deallocate ( u_n, tmp )
      deallocate ( Mmat, S, LLbm, LUbm )
      deallocate ( S1, S2, S3 )
      deallocate ( LLbm1, LLbm2, LLbm3 )


    end subroutine deallocate_arrays


  end subroutine diffuse_interface_tp_coupled_elem


! Internal element routine for the diffuse-interface equation
! using the free energy function by Boyer and Lapuerta (2006)
! (https://doi.org/10.1051/m2an:2006028)
!
!          dc1dt + u.nabla c1 - div( Mcoef  nabla c1 )  = 0
!          dc2dt + u.nabla c2 - div( Mcoef2 nabla c2 )  = 0
!    (...) * c1 + (...) * c2 + kappa1 nabla^2 c1 + mu1 = 0
!    (...) * c1 + (...) * c2 + kappa2 nabla^2 c2 + mu2 = 0
!
!   where (...) depends on the choice for the ternary free energy.
!
!  Boyer and Lapuerta model
!         f(c1,c2,c3) = 1/2 [ alpha * c1^2 * (1-c1)^2 +
!                             beta  * c2^2 * (1-c2)^2 +
!                             gamm  * c3^2 * (1-c3)^2 ]
!                     + 3 * omega * c1^2 * c2^2 * c3^2
!         where c3 = 1 - c1 - c2
!
!         alpha, beta, gamm, omega, Mcoef, kappa, kappa2 and kappa3 are
!         passed with coefficients, while Mcoef2 and Mcoef3 are computed
!         within the element routine by the relationship
!         M1 * kappa = M2 * kappa2 = M3 * kappa3

  subroutine diffuse_interface_tp2_coupled_elem ( mesh, problem, elgrp, &
    elem, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemvec )

    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, ip, i1, i2, i3, i4
    real(dp) :: alpha, beta, gamm, omega
    real(dp) :: Mcoef, Mcoef2, Mcoef3, Mcoefs, kappa, kappa2, kappa3
    real(dp) :: deltat

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

!   some tests

    if ( first ) then

      call check ( coefficients, 'diffuse_interface_tp2_coupled_elem', &
        ncoefi=200, ncoefr=150 )

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   get velocity vector at previous time step

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel], layer=layer )

    tmp = reshape ( u, [ndf,ndim] )

    u_n = matmul ( phi, tmp )

!   un.grad operator

    do ip = 1, ninti
      ungradphi(ip,:) = matmul ( dphidx(ip,:,:), u_n(ip,:) )
    end do

!   get c_1n variable at previous time step for first phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, c_n, physq=[1] , layer=layer )

    c_ng = matmul ( phi, c_n )

!   get c_1iter variable at previous iteration for first phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(3)%p, elgrp, &
      elem, citer, physq=[1] , layer=layer )

    citer_g = matmul ( phi, citer )

!   get c_2n variable at previous time step for the second phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, c2_n, physq=[2] , layer=layer )

    c2_ng = matmul ( phi, c2_n )

!   get c_2iter variable at previous iteration for second phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(3)%p, elgrp, &
      elem, c2iter, physq=[2] , layer=layer )

    c2iter_g = matmul ( phi, c2iter )

!   coefficients

    deltat = coefficients%r(101)
    Mcoef  = coefficients%r(102)
    alpha  = coefficients%r(103)
    beta   = coefficients%r(104)
    kappa  = coefficients%r(105)
    gamm   = coefficients%r(107)
    kappa2 = coefficients%r(109)
    kappa3 = coefficients%r(110)
    omega  = coefficients%r(111)

!   compute Mcoef2 and Mcoef3 from the relationship
!   Mcoef * kappa = Mcoef2 * kapp2 = Mcoef3 * kappa3

    Mcoef2 = Mcoef * kappa / kappa2
    Mcoef3 = Mcoef * kappa / kappa3

    Mcoefs = Mcoef + Mcoef2 + Mcoef3

!   pointers in unknown vector

    i1 = ndf
    i2 = 2*ndf
    i3 = 3*ndf
    i4 = 4*ndf

    if ( matrix ) then

!     diffusion matrix

      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,i,:) * dphidx(ip,j,:) )
          end do
          S(i,j) = sum ( work * detF * wg )
          S(j,i) = S(i,j)
        end do
      end do

!     mass matrix

      do i = 1, ndf
        do j = i, ndf
          Mmat(i,j) = sum ( phi(:,i) * phi(:,j) * detF * wg )
          Mmat(j,i) = Mmat(i,j)
        end do
      end do

!     work * c1 part in mu1

      work = - Mcoef2 / Mcoefs * ( alpha * ( 1._dp - citer_g ) &
                                  * ( 1._dp - 2._dp * citer_g ) &
                                  + 6._dp * omega * c2iter_g  &
                                  * ( 1._dp - citer_g - c2iter_g)**2  &
                                  * ( c2iter_g - citer_g ) ) &
             - Mcoef3 / Mcoefs * ( alpha * ( 1._dp - citer_g ) &
                                  * ( 1._dp - 2._dp * citer_g ) &
                                - gamm * ( 1._dp - citer_g - c2iter_g ) &
                              * ( 2._dp * citer_g + 2._dp * c2iter_g - 1._dp ) &
                                - 6._dp * omega * c2iter_g**2 &
                                * ( 1._dp - citer_g - c2iter_g ) &
                                * ( 2._dp * citer_g + c2iter_g - 1._dp ) )

      do i = 1, ndf
        do j = i, ndf
          LLbm(i,j) = sum ( phi(:,i) * work * phi(:,j) * detF * wg )
          LLbm(j,i) = LLbm(i,j)
        end do
      end do

!     work * c2 part in mu1

      work = Mcoef2 / Mcoefs * ( beta * ( 1._dp - c2iter_g ) &
                                 * ( 1._dp - 2._dp * c2iter_g) ) &
           + Mcoef3 / Mcoefs * ( gamm * (1._dp - citer_g - c2iter_g ) &
                             * ( 2._dp * citer_g + 2._dp * c2iter_g - 1._dp ) )

      do i = 1, ndf
        do j = i, ndf
          LLbm1(i,j) = sum ( phi(:,i) * work * phi(:,j) * detF * wg )
          LLbm1(j,i) = LLbm1(i,j)
        end do
      end do

!     work * c1 part in mu2

      work = Mcoef / Mcoefs * ( alpha * ( 1._dp - citer_g ) &
                                 * ( 1._dp - 2._dp * citer_g) ) &
           + Mcoef3 / Mcoefs * ( gamm * (1._dp - citer_g - c2iter_g ) &
                             * ( 2._dp * citer_g + 2._dp * c2iter_g - 1._dp ) )

      do i = 1, ndf
        do j = i, ndf
          LLbm2(i,j) = sum ( phi(:,i) * work * phi(:,j) * detF * wg )
          LLbm2(j,i) = LLbm2(i,j)
        end do
      end do

!     work * c2 part in mu2

      work = - Mcoef / Mcoefs * ( beta * ( 1._dp - c2iter_g ) &
                                  * ( 1._dp - 2._dp * c2iter_g ) &
                                  - 6._dp * omega * citer_g  &
                                  * ( 1._dp - citer_g - c2iter_g)**2  &
                                  * ( c2iter_g - citer_g ) ) &
             - Mcoef3 / Mcoefs * ( beta * ( 1._dp - c2iter_g ) &
                                  * ( 1._dp - 2._dp * c2iter_g ) &
                                - gamm * ( 1._dp - citer_g - c2iter_g ) &
                              * ( 2._dp * citer_g + 2._dp * c2iter_g - 1._dp ) &
                                + 6._dp * omega * citer_g**2 &
                                * ( 1._dp - citer_g - c2iter_g ) &
                                * ( 1._dp - citer_g - 2._dp * c2iter_g ) )

      do i = 1, ndf
        do j = i, ndf
          LLbm3(i,j) = sum ( phi(:,i) * work * phi(:,j) * detF * wg )
          LLbm3(j,i) = LLbm3(i,j)
        end do
      end do

!     advection term

      do i = 1, ndf
        do j = 1, ndf
          LUbm(i,j) = sum ( phi(:,i) * ungradphi(:,j) * detF * wg )
        end do
      end do


!     Fill element matrix

      elemmat(    1:i1,    1:i1 ) = LUbm + Mmat / deltat
      elemmat(    1:i1, i1+1:i2 ) = 0.0_dp
      elemmat(    1:i1, i2+1:i3 ) = Mcoef * S
      elemmat(    1:i1, i3+1:i4 ) = 0.0_dp

      elemmat( i1+1:i2,    1:i1 ) = 0.0_dp
      elemmat( i1+1:i2, i1+1:i2 ) = LUbm + Mmat / deltat
      elemmat( i1+1:i2, i2+1:i3 ) = 0.0_dp
      elemmat( i1+1:i2, i3+1:i4 ) = Mcoef2 * S

      elemmat( i2+1:i3,    1:i1 ) = LLbm - kappa * S
      elemmat( i2+1:i3, i1+1:i2 ) = LLbm1
      elemmat( i2+1:i3, i2+1:i3 ) = Mmat
      elemmat( i2+1:i3, i3+1:i4 ) = 0.0_dp

      elemmat( i3+1:i4,    1:i1 ) = LLbm2
      elemmat( i3+1:i4, i1+1:i2 ) = LLbm3 - kappa2 * S
      elemmat( i3+1:i4, i2+1:i3 ) = 0.0_dp
      elemmat( i3+1:i4, i3+1:i4 ) = Mmat

    end if

    if ( vector ) then

!     right-hand side of the c-equation

      do i = 1, ndf
        work1(i) = sum ( phi(:,i) * c_ng * detF * wg ) / deltat
        work8(i) = sum ( phi(:,i) * c2_ng * detF * wg ) / deltat
      end do

!     Fill element vector

      elemvec(1:i1)    = work1  ! c1
      elemvec(i1+1:i2) = work8  ! c2
      elemvec(i2+1:i3) = 0.0_dp  ! mu1
      elemvec(i3+1:i4) = 0.0_dp  ! mu2

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

      allocate ( ungradphi(ninti,ndf) )
      allocate ( u(ndim*ndf), c_n(ndf), c_ng(ninti) )
      allocate ( c2_n(ndf), c2_ng(ninti) )
      allocate ( citer(ndf), citer_g(ninti) )
      allocate ( c2iter(ndf), c2iter_g(ninti) )
      allocate ( work(ninti), work1(ndf) )
      allocate ( work8(ndf) )
      allocate ( u_n(ninti,ndim), tmp(ndf,ndim) )
      allocate ( Mmat(ndf,ndf), S(ndf,ndf), LLbm(ndf,ndf), LUbm(ndf,ndf) )
      allocate ( LLbm1(ndf,ndf), LLbm2(ndf,ndf), LLbm3(ndf,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( ungradphi )
      deallocate ( u, c_n, c_ng )
      deallocate ( c2_n, c2_ng )
      deallocate ( citer, citer_g )
      deallocate ( c2iter, c2iter_g )
      deallocate ( work, work1 )
      deallocate ( work8 )
      deallocate ( u_n, tmp )
      deallocate ( Mmat, S, LLbm, LUbm )
      deallocate ( LLbm1, LLbm2, LLbm3 )

    end subroutine deallocate_arrays

  end subroutine diffuse_interface_tp2_coupled_elem


! Internal element routine for the rhs in the momentum balance

  subroutine rhs_tp_coupled_mugradc ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j
    real(dp) :: rho

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


!   some tests

    if ( first ) then

      call check ( coefficients, 'rhs_tp_coupled_mugradc', ncoefi=200, &
        ncoefr=150 )

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )


!   get c variable (concentration) for the first phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, c_n, physq=[1] ,layer=layer )

!   grad c term for the first phase

    do i = 1, ndim
      gradc(:,i) = matmul ( dphidx(:,:,i), c_n )
    end do

!   get c variable (concentration) for the second phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, c2_n, physq=[2] ,layer=layer )

!   grad c term for the second phase

    do i = 1, ndim
      gradc2(:,i) = matmul ( dphidx(:,:,i), c2_n )
    end do

!   get mu variable (chemical potential)) for the first phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, mu, physq=[3] ,layer=layer )

    mu_g = matmul ( phi, mu )

!   get mu variable (chemical potential)) for the second phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, mu2, physq=[4] ,layer=layer )

    mu2_g = matmul ( phi, mu2 )

!   coefficients

    rho  = coefficients%r(106)

    if ( vector ) then

!     term rho (mu1 * gradc1 + mu2 * gradc2)

      do j = 1, ndim
        do i = 1, ndf
          tmp(i,j) = rho * sum ( phi(:,i) * mu_g * gradc(:,j) * detF * wg ) + &
                     rho * sum ( phi(:,i) * mu2_g * gradc2 (:,j) *detF * wg )
        end do
      end do

      elemvec = reshape ( tmp, [ndim*ndf] )

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

      allocate ( tmp(ndf,ndim) )
      allocate ( c_n(ndf), gradc(ninti,ndim) )
      allocate ( mu(ndf), mu_g(ninti) )
      allocate ( c2_n(ndf), gradc2(ninti,ndim) )
      allocate ( mu2(ndf), mu2_g(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( tmp )
      deallocate ( c_n, gradc )
      deallocate ( mu, mu_g )
      deallocate ( c2_n, gradc2 )
      deallocate ( mu2, mu2_g )

    end subroutine deallocate_arrays

  end subroutine rhs_tp_coupled_mugradc


! Internal element routine for the rhs in the momentum balance
! In Boyer and Lapuerta model, the rhs is
!   mu1 * grad c1 + mu2 * grad c2 + mu3 * grad c3
! where mu3 can be obtained from the relationship
!   Mcoef * mu1 + Mcoef2 * mu2 + Mcoef3 * mu3 = 0
! and grad c3 = - grad c1 - grad c2

  subroutine rhs_tp2_coupled_mugradc ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j
    real(dp) :: rho, Mcoef, Mcoef2, Mcoef3
    real(dp) :: kappa, kappa2, kappa3

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

!   some tests

    if ( first ) then

      call check ( coefficients, 'rhs_tp2_coupled_mugradc', ncoefi=200, &
        ncoefr=150 )

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   get c variable (concentration) for the first phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, c_n, physq=[1] ,layer=layer )

!   grad c term for the first phase

    do i = 1, ndim
      gradc(:,i) = matmul ( dphidx(:,:,i), c_n )
    end do

!   get c variable (concentration) for the second phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, c2_n, physq=[2] ,layer=layer )

!   grad c term for the second phase

    do i = 1, ndim
      gradc2(:,i) = matmul ( dphidx(:,:,i), c2_n )
    end do

!   get mu variable (chemical potential)) for the first phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, mu, physq=[3] ,layer=layer )

    mu_g = matmul ( phi, mu )

!   get mu variable (chemical potential)) for the second phase

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, mu2, physq=[4] ,layer=layer )

    mu2_g = matmul ( phi, mu2 )

!   coefficients

    Mcoef  = coefficients%r(102)
    kappa  = coefficients%r(105)
    rho    = coefficients%r(106)
    kappa2 = coefficients%r(109)
    kappa3 = coefficients%r(110)

!   compute Mcoef2 and Mcoef3 from the relationship
!   Mcoef * kappa = Mcoef2 * kapp2 = Mcoef3 * kappa3

    Mcoef2 = Mcoef * kappa / kappa2
    Mcoef3 = Mcoef * kappa / kappa3

    if ( vector ) then

!     term rho (mu1 * grad c1 + mu2 * grad hc2 + mu3 * grad c3)

      do j = 1, ndim
        do i = 1, ndf
          tmp(i,j) = rho * sum ( phi(:,i) * ( mu_g + Mcoef / Mcoef3 * mu_g &
                        + Mcoef2 / Mcoef3 * mu2_g ) * gradc(:,j) * detF * wg ) &
                   + rho * sum ( phi(:,i) * ( mu2_g + Mcoef / Mcoef3 * mu_g &
                         + Mcoef2 / Mcoef3 * mu2_g ) * gradc2(:,j) * detF *wg )
        end do
      end do

      elemvec = reshape ( tmp, [ndim*ndf] )

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

      allocate ( tmp(ndf,ndim) )
      allocate ( c_n(ndf), gradc(ninti,ndim) )
      allocate ( mu(ndf), mu_g(ninti) )
      allocate ( c2_n(ndf), gradc2(ninti,ndim) )
      allocate ( mu2(ndf), mu2_g(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( tmp )
      deallocate ( c_n, gradc )
      deallocate ( mu, mu_g )
      deallocate ( c2_n, gradc2 )
      deallocate ( mu2, mu2_g )

    end subroutine deallocate_arrays

  end subroutine rhs_tp2_coupled_mugradc


end module diffuse_interface_tp_coupled_elements_generic_m
