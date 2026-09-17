!
! Copyright (C) 2004-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system (e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Element routines for the interfacial boundary condition:
!             __
!  t| + t|  = \/ . tau
!    d    m     s     s
!
! using different interfacial rheological models.
!
! This module contains stuff that only works for curves (in 2D space).

module interfacial_rheology_elements_2D_curve_m

  use tfem_elem_m
  use interfacial_rheology_elements_generic_m

  implicit none

contains

! Boundary element for surface tension on a curve (implicit) with viscous,
! elastic, or viscoelastic (quasi-linear Kelvin-Voigt) interface conditions.
! NOTE: only in combination with Lagrangian interface tracking

  subroutine interfacial_rheology_curve ( mesh, problem, curve, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use interfacial_rheology_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, M, i, j, ip, i1, i2, imodel, surfimpl
    real(dp) :: deltat, gammac, kappac, Kc, muc, Gc

    if ( first ) then

!     first element on this curve

!     set globals

      call set_globals_interfacial_bc ( mesh, coefficients, curve=curve )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( xig(ninti,ndim-1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim), xng(ninti,ndim) )
      allocate ( xar(nodalp*ndim), xn(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,ndim-1), dxdxis(ninti,ndim,ndim-1) )
      allocate ( gi_up(ninti,ndim,ndim-1), dxndxis(ninti,ndim,ndim-1) )
      allocate ( gij_up(ninti,ndim-1,ndim-1) )
      allocate ( work(ninti) )

      allocate ( x0(nodalp,ndim) )
      allocate ( dxdxis0(ninti,ndim,ndim-1) )
      allocate ( surfl0(ninti) )
      allocate ( xg0(ninti,ndim) )
      allocate ( gi_up0(ninti,ndim,ndim-1) )

      allocate ( Fs(ninti,ndim,ndim), Bs(ninti,ndim,ndim) )
      allocate ( Jinv(ninti), Jdet(ninti), trBs(ninti) )
      allocate ( work_upup(ninti,ndim,ndim), work_upup2(ninti,ndim,ndim) )
      allocate ( work_gijup(ninti), work_updown(ninti,ndim,ndim) )
      allocate ( work2(ndim,ndim), work3(ndim,ndim), work4(ndim,ndim) )
      allocate ( work8(ndim,ndim) )
      allocate ( work_I(ninti), work_B(ninti), work_axi(ninti) )
      allocate ( Suu(ndf,ndf), Svv(ndf,ndf), Suv(ndf,ndf) )
      allocate ( work5(ndf,ndim), work6(ndf,ndim), work7(ndf,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

!   get interface position and shape at tn

    if ( coefficients%i(12) == 1 ) then

      call get_vector ( oldvectors%m(1)%p, oldvectors%p(1)%p, &
        oldvectors%v(1)%p, elgrp=1, elem=elem, u=xar )

      xn = reshape ( xar, [nodalp,ndim] )

      call isoparametric_deformation_curved ( xn, dphi, dxndxis, surfl )

      call isoparametric_coordinates ( xn, phi, xng )

    end if

!   get interface position and shape at tn+1

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )

    call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, &
      gij_up=gij_up, gi_up=gi_up )

    call isoparametric_coordinates ( x, phi, xg )

    if ( coorsys == 1 ) then
      surfl = 2 * pi * xg(:,2) * surfl
    end if

    if ( any ( coefficients%i(11) == [(i, i=2,5)] ) ) then

!     get interface position and shape at t0

      call get_coordinates_geometry ( oldvectors%m(2)%p, elem, x0, curve=curve )

      call isoparametric_deformation_curved ( x0, dphi, dxdxis0, surfl0, &
        gi_up=gi_up0 )

      call isoparametric_coordinates ( x0, phi, xg0 )

!     deformation gradient and Finger tensor

      do i = 1, ndim
        do j = 1, ndim
          Fs(:,i,j) = dxdxis(:,i,1)*gi_up0(:,j,1)  ! Deformation gradient tensor
        end do
      end do

      do ip = 1, ninti
        Bs(ip,:,:) = matmul ( Fs(ip,:,:), transpose(Fs(ip,:,:)) )  ! Finger
      end do

      if ( coorsys == 1 ) surfl0 = 2 * pi * xg0(:,2) * surfl0

      Jinv = surfl0/surfl
      Jdet = surfl/surfl0

    end if

!   surface tension coefficient

    deltat = coefficients%r(1)    ! time step
    imodel = coefficients%i(11)   ! interface model number
    surfimpl = coefficients%i(12) ! implicit or explicit surface tension

    gammac = coefficients%r(2) ! surface tension coefficient
    kappac = coefficients%r(3) ! surface dilatational viscosity
    muc    = coefficients%r(4) ! surface shear viscosity
    Kc     = coefficients%r(5) ! surface dilatational elasticity
    Gc     = coefficients%r(6) ! surface shear elasticity

    i1 = ndf
    i2 = 2*ndf

    if ( matrix ) then

    elemmat = 0

!   implicit surface tension term
    if ( surfimpl == 1 ) then
      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphi(ip,i,:) * &
                       matmul ( gij_up(ip,:,:), dphi(ip,j,:) ) )
          end do
          elemmat(i,j) = sum ( work * surfl * wg )
          elemmat(j,i) = elemmat(i,j)
        end do
      end do

      do i = 1+ndf, 2*ndf
        do j = i, 2*ndf
          do ip = 1, ninti
            work(ip) = sum ( dphi(ip,i-ndf,:) * &
                       matmul ( gij_up(ip,:,:), dphi(ip,j-ndf,:) ) )
          end do
          if ( coorsys == 1 ) then
            work = work + phi(:,i-ndf) * phi(:,j-ndf) / xg(:,2)**2
          end if
          elemmat(i,j) = sum ( work * surfl * wg )
          elemmat(j,i) = elemmat(i,j)
        end do
      end do

        elemmat = gammac * deltat * elemmat

    end if

    select case ( imodel )

    case (1,3,5) ! Boussinesq-Scriven, Kelvin-Voigt

!     rate-of-deformation and divergence terms
      do N = 1, ndf
        do M = N, ndf
          do ip = 1, ninti
            work_gijup(ip) = sum ( dphi(ip,N,:) * &
                  matmul ( gij_up(ip,:,:), dphi(ip,M,:) ) )
          end do

          do i = 1, ndim
            do j = 1, ndim
              work_upup(:,i,j) = sum ( dphi(:,N,:) * gi_up(:,i,:), 2 ) * &
                                 sum ( dphi(:,M,:) * gi_up(:,j,:), 2 )
              work_upup2(:,i,j) = sum ( dphi(:,N,:) * gi_up(:,j,:), 2 ) * &
                                  sum ( dphi(:,M,:) * gi_up(:,i,:), 2 )
              work_updown(:,i,j) = sum ( gi_up(:,i,:) * dxdxis(:,j,:), 2 )
            end do
          end do

          if ( coorsys == 1 ) then
            work_upup(:,1,2) =   sum ( dphi(:,N,:) * gi_up(:,1,:), 2 ) * &
                               ( sum ( dphi(:,M,:) * gi_up(:,2,:), 2 ) + &
                                 phi(:,M) / xg(:,2) )
            work_upup(:,2,1) = ( sum ( dphi(:,N,:) * gi_up(:,2,:), 2 ) + &
                                 phi(:,N) / xg(:,2) ) * &
                                 sum ( dphi(:,M,:) * gi_up(:,1,:), 2 )
            work_upup(:,2,2) = ( sum ( dphi(:,N,:) * gi_up(:,2,:), 2 ) + &
                                 phi(:,N) / xg(:,2) ) * &
                               ( sum ( dphi(:,M,:) * gi_up(:,2,:), 2 ) + &
                                 phi(:,M) / xg(:,2) )
          end if

          do i = 1, ndim
            do j = 1, ndim
              work2(i,j) = sum ( work_upup(:,i,j) * surfl * wg )
              work3(i,j) = sum ( work_upup2(:,i,j) * surfl * wg )
              work4(i,j) = sum ( work_gijup * work_updown(:,i,j) * &
                                  surfl * wg )
            end do
          end do

          if ( coorsys == 1 ) then
            work_axi = phi(:,N) * phi(:,M) / xg(:,2)**2
            work3(2,2) = work3(2,2) + sum ( work_axi * surfl * wg )
            work4(2,2) = work4(2,2) + sum ( work_axi * surfl * wg )
          end if

          work8 = ( kappac - muc ) * work2 + muc * ( work3 + work4 )

          Suu(N,M) = work8(1,1)
          Suu(M,N) = Suu(N,M)
          Suv(N,M) = work8(1,2)
          Suv(M,N) = work8(2,1)
          Svv(N,M) = work8(2,2)
          Svv(M,N) = Svv(N,M)
        end do
      end do

      elemmat(    1:i1,    1:i1 ) = elemmat( 1:i1, 1:i1 ) + Suu
      elemmat(    1:i1, i1+1:i2 ) = elemmat( 1:i1, i1+1:i2 ) + Suv
      elemmat( i1+1:i2,    1:i1 ) = elemmat( i1+1:i2, 1:i1 ) + transpose(Suv)
      elemmat( i1+1:i2, i1+1:i2 ) = elemmat( i1+1:i2, i1+1:i2 ) + Svv

    case(0,2,4) ! Constant surface tension or generalized Hooke's law
    case default
      write(*,'(a,a,i0,a,a)') &
        'Error in interfacial_rheology_curve: ', &
        'cannot build system matrix for interface model = ', imodel, &
        '. Choose different interface model number or call build_system ', &
        'with buildmatrix=.false.'
      stop
    end select

    end if

    if ( vector ) then

!   surface tension

    do i = 1, ndim
      do N = 1, ndf
      if ( surfimpl == 1 ) then
        do ip = 1, ninti
          work_I(ip) = sum ( dphi(ip,N,:) * matmul ( gij_up(ip,:,:), &
            dxndxis(ip,i,:) ) )
        end do
      else
        work_I = dphi(:,N,1) * gi_up(:,i,1)
      end if
      work5(N,i) = sum ( work_I * surfl * wg )
      end do
    end do

    if ( coorsys == 1 ) then
      do N = 1, ndf
        if ( surfimpl == 1 ) then
          work_I = phi(:,N) * xng(:,2) / xg(:,2)**2
        else
          work_I = phi(:,N) / xg(:,2)
        end if
        work5(N,2) = work5(N,2) + sum ( work_I * surfl * wg )
      end do
    end if

    elemvec = reshape ( - gammac * work5, [ndf*ndim] )

    select case ( imodel )

    case (0,1) ! Surface tension only or Boussinesq-Scriven
    case (2:5) ! Generalized Hooke's law, Kelvin-Voigt

      if ( any ( imodel == [2,3] ) ) then

!       Neo-Hookean model

        do i = 1, ndim
          do N = 1, ndf
            work_I = dphi(:,N,1) * gi_up(:,i,1)
            work_B = dphi(:,N,1) * sum ( gi_up(:,:,1) * Bs(:,:,i), 2 )
            work5(N,i) = sum ( work_I * surfl * wg )
            work6(N,i) = sum ( log(Jdet) * work_I * surfl * wg )
            work7(N,i) = sum ( Jinv * work_B * surfl * wg )
          end do
        end do

        if ( coorsys == 1 ) then
          do N = 1, ndf
            work_I = phi(:,N) / xg(:,2)
            work_B = work_I * ( xg(:,2)/xg0(:,2) )**2
            work5(N,2) = work5(N,2) + sum ( work_I * surfl * wg )
            work6(N,2) = work6(N,2) + sum ( log(Jdet) * work_I * surfl * wg )
            work7(N,2) = work7(N,2) + sum ( Jinv * work_B * surfl * wg )
          end do
        end if

      else

!       Hütter-Tervoort model

        do ip = 1, ninti
          trBs(ip) = Bs(ip,1,1) + Bs(ip,2,2)
          if ( coorsys == 1 ) trBs(ip) = trBs(ip) + ( xg(ip,2)/xg0(ip,2) )**2
        end do
        trBs = trBs / 2._dp

        do i = 1, ndim
          do N = 1, ndf
            work_I = dphi(:,N,1) * gi_up(:,i,1)
            work_B = dphi(:,N,1) * sum ( gi_up(:,:,1) * Bs(:,:,i), 2 )
            work5(N,i) = sum ( Jinv**2 * trBs * work_I * surfl * wg )
            work6(N,i) = sum ( Jinv * log(Jdet) * work_I * surfl * wg )
            work7(N,i) = sum ( Jinv**2 * work_B * surfl * wg )
          end do
        end do

        if ( coorsys == 1 ) then
          do N = 1, ndf
            work_I = phi(:,N) / xg(:,2)
            work_B = work_I * ( xg(:,2)/xg0(:,2) )**2
            work5(N,2) = work5(N,2) + &
              sum ( Jinv**2 * trBs * work_I * surfl * wg )
            work6(N,2) = work6(N,2) + &
              sum ( Jinv * log(Jdet) * work_I * surfl * wg )
            work7(N,2) = work7(N,2) + sum ( Jinv**2 * work_B * surfl * wg )
          end do
        end if

      end if

      elemvec = elemvec + reshape ( ( - Kc * work6 &
                                      - Gc * work7 &
                                      + Gc * work5 ), [ndf*ndim] )

    case default
      write(*,'(a,a,i0,a,a)') &
        'Error in interfacial_rheology_curve: ', &
        'cannot build system vector for interface model = ', imodel, &
        '. Choose different interface model number or call build_system ', &
        'with buildvector=.false.'
      stop
    end select

    end if

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl )
      deallocate ( xig, phi, x )
      deallocate ( xg, xng, xar, xn )
      deallocate ( dphi, dxdxis )
      deallocate ( gi_up, dxndxis )
      deallocate ( gij_up )
      deallocate ( work )

      deallocate ( x0 )
      deallocate ( dxdxis0 )
      deallocate ( surfl0 )
      deallocate ( xg0 )
      deallocate ( gi_up0 )

      deallocate ( Fs, Bs )
      deallocate ( Jinv, Jdet, trBs )
      deallocate ( work_gijup, work_upup )
      deallocate ( work_upup2, work_updown )
      deallocate ( work2, work3, work4, work8 )
      deallocate ( work_I, work_B, work_axi )
      deallocate ( Suu, Svv, Suv )
      deallocate ( work5, work6, work7 )

    end if

  end subroutine interfacial_rheology_curve

end module interfacial_rheology_elements_2D_curve_m
