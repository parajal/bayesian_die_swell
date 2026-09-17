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
! This module contains stuff that only works for surfaces (in 3D space).

module interfacial_rheology_elements_3D_surface_m

  use tfem_elem_m
  use interfacial_rheology_elements_generic_m

  implicit none

contains

! Boundary element for surface tension on a curve (implicit) with viscous,
! elastic, or viscoelastic (quasi-linear Kelvin-Voigt) interface conditions.
! NOTE: only in combination with Lagrangian interface tracking

  subroutine interfacial_rheology_surface ( mesh, problem, surface, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use interfacial_rheology_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: surface, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, M, ip, i, j, k, i1, i2, i3, imodel, surfimpl
    real(dp) :: gammac, deltat, kappac, Kc, muc, Gc


    if ( first ) then

!     first element on this surface

!     set globals

      call set_globals_interfacial_bc ( mesh, coefficients, &
        surface=surface )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xar(nodalp*ndim), xn(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )
      allocate ( gi_up(ninti,ndim,2), dxndxis(ninti,ndim,2) )
      allocate ( gij_up(ninti,2,2) )
      allocate ( work(ninti) )

      allocate ( x0(nodalp,ndim) )
      allocate ( dxdxis0(ninti,ndim,2) )
      allocate ( surfl0(ninti) )
      allocate ( gi_up0(ninti,ndim,2) )

      allocate ( Fs(ninti,ndim,ndim), Bs(ninti,ndim,ndim) )
      allocate ( Jinv(ninti), Jdet(ninti) )
      allocate ( work_upup(ninti,ndim,ndim), trBs(ninti) )
      allocate ( work_gijup(ninti), work_updown(ninti,ndim,ndim) )
      allocate ( work2(ndim,ndim), work3(ndim,ndim), work4(ndim,ndim) )
      allocate ( work5(ndim,ndim), work6(ndf,ndim), work7(ndf,ndim) )
      allocate ( work8(ndf,ndim), work_I(ninti), work_B(ninti) )
      allocate ( Suu(ndf,ndf), Suv(ndf,ndf), Suw(ndf,ndf) )
      allocate ( Svv(ndf,ndf), Svw(ndf,ndf), Sww(ndf,ndf) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

!   get interface position and shape at tn

    if ( coefficients%i(12) == 1 ) then

      call get_vector ( oldvectors%m(1)%p, oldvectors%p(1)%p, &
        oldvectors%v(1)%p, elgrp=1, elem=elem, u=xar )

      xn = reshape ( xar, [nodalp,ndim] )

      call isoparametric_deformation_surface ( xn, dphi, dxndxis, surfl )

    end if

!   get interface position and shape at tn+1

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, &
      gij_up=gij_up, gi_up=gi_up )

!   compute deformation gradient Fs and Finger tensor Bs for elastic models

    if ( any ( coefficients%i(11) == [(i, i=2,5)] ) ) then

!     get interface position and shape at t0

      call get_coordinates_geometry ( oldvectors%m(2)%p, elem, x0, &
        surface=surface )
      call isoparametric_deformation_surface ( x0, dphi, dxdxis0, surfl0, &
        gi_up=gi_up0 )

!     Deformation gradient tensor

      do i = 1, ndim
        do j = 1, ndim
          Fs(:,i,j) = sum(dxdxis(:,i,:)*gi_up0(:,j,:),2)
        end do
      end do

!     Finger tensor

      do ip = 1, ninti
        Bs(ip,:,:) = matmul ( Fs(ip,:,:), transpose(Fs(ip,:,:)) )
      end do

!     Determinant of Fs

      Jinv = surfl0/surfl
      Jdet = surfl/surfl0

    end if

!   material and numerical parameters

    deltat = coefficients%r(1)    ! time step
    imodel = coefficients%i(11)   ! interface model number
    surfimpl = coefficients%i(12) ! implicit or explicit surface tension

    gammac = coefficients%r(2) ! surface tension
    kappac = coefficients%r(3) ! surface dilatational viscosity
    muc    = coefficients%r(4) ! surface shear viscosity
    Kc     = coefficients%r(5) ! surface dilatational elasticity
    Gc     = coefficients%r(6) ! surface shear elasticity

    i1 = ndf
    i2 = 2*ndf
    i3 = 3*ndf

    if ( matrix ) then

    elemmat = 0

!   implicit surface tension term

    if ( surfimpl == 1 ) then

      do k = 1, 3
        do i = (k-1)*ndf+1, k*ndf
          do j = i, k*ndf
            do ip = 1, ninti
              work(ip) = sum ( dphi(ip,i-(k-1)*ndf,:) * &
                         matmul ( gij_up(ip,:,:), dphi(ip,j-(k-1)*ndf,:) ) )
            end do
            elemmat(i,j) = sum ( work * surfl * wg )
            elemmat(j,i) = elemmat(i,j)
          end do
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
              work_updown(:,i,j) = sum ( gi_up(:,i,:) * dxdxis(:,j,:), 2 )
            end do
          end do

          do i = 1, ndim
            do j = 1, ndim
              work2(i,j) = sum ( work_upup(:,i,j) * surfl * wg )
              work3(i,j) = sum ( work_upup(:,j,i) * surfl * wg )
              work4(i,j) = sum ( work_gijup * work_updown(:,i,j) * &
                                  surfl * wg )
            end do
          end do

          work5 = ( kappac - muc ) * work2 + muc * ( work3 + work4 )
          Suu(N,M) = work5(1,1)
          Suu(M,N) = Suu(N,M)
          Suv(N,M) = work5(1,2)
          Suv(M,N) = work5(2,1)
          Suw(N,M) = work5(1,3)
          Suw(M,N) = work5(3,1)
          Svv(N,M) = work5(2,2)
          Svv(M,N) = Svv(N,M)
          Svw(N,M) = work5(2,3)
          Svw(M,N) = work5(3,2)
          Sww(N,M) = work5(3,3)
          Sww(M,N) = Sww(N,M)
        end do
      end do

      elemmat(    1:i1,    1:i1 ) = elemmat( 1:i1, 1:i1 ) + Suu
      elemmat(    1:i1, i1+1:i2 ) = elemmat( 1:i1, i1+1:i2 ) + Suv
      elemmat( i1+1:i2,    1:i1 ) = elemmat( i1+1:i2, 1:i1 ) + transpose(Suv)
      elemmat(    1:i1, i2+1:i3 ) = elemmat( 1:i1, i2+1:i3 ) + Suw
      elemmat( i2+1:i3,    1:i1 ) = elemmat( i2+1:i3, 1:i1 ) + transpose(Suw)
      elemmat( i1+1:i2, i1+1:i2 ) = elemmat( i1+1:i2, i1+1:i2 ) + Svv
      elemmat( i1+1:i2, i2+1:i3 ) = elemmat( i1+1:i2, i2+1:i3 ) + Svw
      elemmat( i2+1:i3, i1+1:i2 ) = elemmat( i2+1:i3, i1+1:i2 ) + transpose(Svw)
      elemmat( i2+1:i3, i2+1:i3 ) = elemmat( i2+1:i3, i2+1:i3 ) + Sww

    case(0,2,4) ! Constant surface tension or generalized Hooke's law
    case default
      write(*,'(a,a,i0,a,a)') &
        'Error in interfacial_rheology_surface: ', &
        'cannot build system matrix for interface model = ', imodel, &
        '. Choose different interface model number or call build_system ', &
        'with buildmatrix=.false.'
      stop
    end select

    end if

    if ( vector ) then

!   constant surface tension

    do i = 1, ndim
      do N = 1, ndf
      if ( surfimpl == 1 ) then
        do ip = 1, ninti
          work_I(ip) = sum ( dphi(ip,N,:) * matmul ( gij_up(ip,:,:), &
            dxndxis(ip,i,:) ) )
        end do
      else
        work_I = sum(dphi(:,N,:) * gi_up(:,i,:),2)
      end if
      work6(N,i) = sum ( work_I * surfl * wg )
      end do
    end do

    elemvec = reshape ( - gammac * work6, [ndf*ndim] )

    select case ( imodel )
    case (0,1) ! Surface tension only or Boussinesq-Scriven
    case (2:5) ! Generalized Hooke's law, Kelvin-Voigt
      do i = 1, ndim
        do N = 1, ndf
          work_I = sum(dphi(:,N,:)*gi_up(:,i,:),2)
          do ip = 1, ninti
            work_B(ip) = dot_product ( matmul ( gi_up(ip,:,:), &
              dphi(ip,N,:) ), Bs(ip,:,i) )
          end do
          if ( any ( imodel == [2,3] ) ) then
            work6(N,i) = sum ( work_I * surfl * wg )
            work7(N,i) = sum ( log(Jdet) * work_I * surfl * wg )
            work8(N,i) = sum ( Jinv * work_B * surfl * wg )
          else
            do ip = 1, ninti
              trBs(ip) = 0.5_dp * ( Bs(ip,1,1) + Bs(ip,2,2) + Bs(ip,3,3) )
            end do
            work6(N,i) = sum ( Jinv**2 * trBs * work_I * surfl * wg )
            work7(N,i) = sum ( Jinv * log(Jdet) * work_I * surfl * wg )
            work8(N,i) = sum ( Jinv**2 * work_B * surfl * wg )
          end if
        end do
      end do

      elemvec = elemvec + reshape ( ( - Kc * work7 &
                                      - Gc * work8 &
                                      + Gc * work6 ), [ndf*ndim] )
    case default
      write(*,'(a,a,i0,a,a)') &
        'Error in interfacial_rheology_surface: ', &
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
      deallocate ( xar, xn )
      deallocate ( dphi, dxdxis )
      deallocate ( gi_up, dxndxis )
      deallocate ( gij_up )
      deallocate ( work )

      deallocate ( x0 )
      deallocate ( dxdxis0 )
      deallocate ( surfl0 )
      deallocate ( gi_up0 )

      deallocate ( Fs, Bs )
      deallocate ( Jinv, Jdet )
      deallocate ( work_gijup, work_upup )
      deallocate ( work_updown, trBs )
      deallocate ( work2, work3, work4, work5 )
      deallocate ( work_I, work_B )
      deallocate ( Suu, Suv, Suw, Svv, Svw, Sww )
      deallocate ( work6, work7, work8 )

    end if

  end subroutine interfacial_rheology_surface


! Compute interfacial stress in all nodes

  subroutine interfacial_stress_ten ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    use interfacial_rheology_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    real(dp) :: gammac, kappac, Kc, muc, Gc
    integer :: i, j, k, imodel

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interfacial_stress ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

!     shape function in the nodes

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

!   interface model nr and material parameters

    imodel = coefficients%i(11)
    gammac = coefficients%r(2)
    kappac = coefficients%r(3)
    muc = coefficients%r(4)
    Kc = coefficients%r(5)
    Gc = coefficients%r(6)

!   position

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, &
     gi_up=gi_up )

    if ( any ( imodel == [(i, i=2,5)] ) ) then
      call get_coordinates ( oldvectors%m(1)%p, elgrp, elem, x0 )
      call isoparametric_deformation_surface ( x0, dphi, dxdxis0, surfl0, &
       gi_up=gi_up0 )
    end if

!   interfacial tension
    do k = 1, nodalp
      do i = 1, ndim
        do j = 1, ndim
          work9(k,i,j) = dot_product ( gi_up(k,i,:), dxdxis(k,j,:) )
        end do
      end do
    end do
    work12 = gammac*work9

!   viscous stress
    select case ( imodel )
    case ( 0,2,4 )
    case ( 1,3,5 )
!     velocities in nodal points
      call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )

      uvec = reshape ( u, [ndf,ndim] )

!     div_s u
      do i = 1, ndim-1
        work10(:,:,i) = matmul ( dphi (:,:,i), uvec )
      end do
      do i = 1, nodalp
        work(i) = sum ( [(dot_product ( gi_up(i,j,:), &
          work10(i,j,:) ), j=1,ndim )] )
      end do

!     grad_s u
      do k = 1, nodalp
        do i = 1, ndim
          do j = 1, ndim
            work11(k,i,j) = dot_product ( gi_up(k,i,:), work10(k,j,:) )
          end do
        end do
      end do
!     total viscous stress
      do i = 1, nodalp
        work12(i,:,:) = work12(i,:,:) + &
                        (kappac-muc)*work(i)*work9(i,:,:) + &
                        muc*( matmul(work11(i,:,:), work9(i,:,:)) + &
                              matmul(work9(i,:,:), transpose(work11(i,:,:))) )
      end do
    case default
      print *, 'Erorr in interfacial_stress_ten: interface model imodel = ', &
        imodel, ' not available.'
      stop
    end select

!   elastic stress
    select case ( imodel )
    case ( 0,1 )
    case ( 2:5 )
!     Determinant of F_s
      work1 = surfl / surfl0

!     total elastic stress
      do k = 1, nodalp
        do i = 1, ndim
          do j = 1, ndim
!           F_s
            work11(k,i,j) = dot_product ( dxdxis(k,i,:), gi_up0(k,j,:) )
          end do
        end do
!       B_s
        work11(k,:,:) = matmul ( work11(k,:,:), transpose(work11(k,:,:)) )
      end do

      if ( any ( imodel == [2,3] ) ) then
        do i = 1, nodalp
          work12(i,:,:) = work12(i,:,:) + &
                          Kc*log(work1(i))*work9(i,:,:) + &
                          Gc*( 1._dp/work1(i)*work11(i,:,:) - work9(i,:,:) )
        end do
      else
        do i = 1, nodalp
          work(i) = 0.5_dp * ( work11(i,1,1) + work11(i,2,2) + work11(i,3,3) )
          work12(i,:,:) = work12(i,:,:) + &
                          Kc*log(work1(i))/work1(i) * work9(i,:,:) + &
                          Gc / work1(i)**2 * ( work11(i,:,:) - &
                                               work(i) * work9(i,:,:) )
        end do
      end if
    case default
      call errormsg_case_default ( 'interfacial_stress_ten', &
        'imodel', int_value=imodel )
    end select

!   interfacial stress tensor

    do i = 1, nodalp
      elemvec(i:i+5*nodalp:nodalp) = [work12(i,1,1:3), work12(i,2,2:3), &
        work12(i,3,3)]
    end do

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xrnod(nodalp,ndim-1), x0(nodalp,ndim), u(ndf*ndim) )
      allocate ( phi(nodalp,ndf), dphi(nodalp,ndf,ndim-1) )
      allocate ( x(nodalp,ndim), uvec(ndf,ndim) )
      allocate ( dxdxis(nodalp,ndim,ndim-1), surfl(nodalp) )
      allocate ( dxdxis0(nodalp,ndim,ndim-1), surfl0(nodalp) )
      allocate ( gi_up(nodalp,ndim,ndim-1), gi_up0(nodalp,ndim,ndim-1) )
      allocate ( work(nodalp), work1(nodalp) )
      allocate ( work10(nodalp,ndim,ndim-1), work9(nodalp,ndim,ndim) )
      allocate ( work11(nodalp,ndim,ndim), work12(nodalp,ndim,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xrnod, x0, u )
      deallocate ( phi, dphi )
      deallocate ( x, uvec )
      deallocate ( dxdxis, surfl )
      deallocate ( dxdxis0, surfl0 )
      deallocate ( gi_up, gi_up0 )
      deallocate ( work, work1 )
      deallocate ( work9, work10 )
      deallocate ( work11, work12 )

    end subroutine deallocate_arrays

  end subroutine interfacial_stress_ten


! Compute interfacial stress in all nodes

  subroutine interfacial_Bs_ten ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    use interfacial_rheology_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    integer :: i, j, k

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interfacial_stress ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

!     shape function in the nodes

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

!   position

    call get_coordinates ( mesh, elgrp, elem, x )
    call get_coordinates ( oldvectors%m(1)%p, elgrp, elem, x0 )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, &
     gi_up=gi_up )
    call isoparametric_deformation_surface ( x0, dphi, dxdxis0, surfl0, &
     gi_up=gi_up0 )

!   B_s tensor
    do k = 1, nodalp
      do i = 1, ndim
        do j = 1, ndim
!         F_s
          work3(i,j) = dot_product ( dxdxis(k,i,:), gi_up0(k,j,:) )
        end do
      end do
!     B_s
      work3 = matmul ( work3, transpose(work3) )
      elemvec(k:k+5*nodalp:nodalp) = [work3(1,1:3), work3(2,2:3), &
        work3(3,3)]
    end do

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xrnod(nodalp,ndim-1), x0(nodalp,ndim) )
      allocate ( phi(nodalp,ndf), dphi(nodalp,ndf,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( dxdxis(nodalp,ndim,ndim-1), surfl(nodalp) )
      allocate ( dxdxis0(nodalp,ndim,ndim-1), surfl0(nodalp) )
      allocate ( gi_up(nodalp,ndim,ndim-1), gi_up0(nodalp,ndim,ndim-1) )
      allocate ( work3(ndim,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xrnod, x0 )
      deallocate ( phi, dphi )
      deallocate ( x )
      deallocate ( dxdxis, surfl )
      deallocate ( dxdxis0, surfl0 )
      deallocate ( gi_up, gi_up0 )
      deallocate ( work3 )

    end subroutine deallocate_arrays

  end subroutine interfacial_Bs_ten

! Compute part of interfacial stress in all nodes

  subroutine interfacial_stress_scal ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    use interfacial_rheology_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    integer :: i, j, k, pchoice

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interfacial_stress ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

!     shape function in the nodes

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

!   position

    call get_coordinates ( mesh, elgrp, elem, x )
    call get_coordinates ( oldvectors%m(1)%p, elgrp, elem, x0 )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, &
     gi_up=gi_up, normal=normal )
    call isoparametric_deformation_surface ( x0, dphi, dxdxis0, surfl0 )


!   velocities in nodal points

    call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )

    uvec = reshape ( u, [ndf,ndim] )

!   du/dxi_i
    do i = 1, ndim-1
      work9(:,:,i) = matmul ( dphi (:,:,i), uvec )
    end do

!   1/2 * L_s = 1/2 * grad_s u * I_s
    do k = 1, nodalp
      do i = 1, ndim
        do j = 1, ndim
!         I_s
          tmp(i,j) = dot_product ( gi_up(k,i,:), dxdxis(k,j,:) )
!         grad_s u
          work3(i,j) = dot_product ( gi_up(k,i,:), work9(k,j,:) )
        end do
      end do
      work10(k,:,:) = 0.5_dp * matmul(work3, tmp)
    end do

    pchoice = coefficients%i(13)

    select case ( pchoice )
    case ( 1 )
!     div_s u
      do i = 1, nodalp
        elemvec(i) = sum ( [(2._dp * work10(i,j,j), j=1,ndim )] )
      end do
    case ( 2 )
!     rate of strain
      do i = 1, nodalp
        elemvec(i) = sqrt ( 2._dp * sum ( ( work10(i,:,:) + &
          transpose ( work10(i,:,:) ) ) ** 2 ) )
      end do
    case ( 3 )
!     surface vorticity
      do i = 1, nodalp
        work3 = work10(i,:,:) - transpose ( work10(i,:,:) )
        work = [ work3(2,3), work3(3,1), work3(1,2) ]
        elemvec(i) = 2._dp * dot_product ( work, normal(i,:) )
      end do
    case ( 4 )
!     determinant of F_s
      elemvec = surfl / surfl0
    case default
      print *, 'Error in interfacial_stress_scal: set ', &
        'coefficients%i(13) to 1-4.'
      stop
    end select

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xrnod(nodalp,ndim-1), x0(nodalp,ndim), u(ndf*ndim) )
      allocate ( phi(nodalp,ndf), dphi(nodalp,ndf,ndim-1) )
      allocate ( x(nodalp,ndim), uvec(ndf,ndim) )
      allocate ( dxdxis(nodalp,ndim,ndim-1), surfl(nodalp) )
      allocate ( dxdxis0(nodalp,ndim,ndim-1), surfl0(nodalp) )
      allocate ( gi_up(nodalp,ndim,ndim-1) )
      allocate ( work(ndim), work3(ndim,ndim) )
      allocate ( tmp(ndim,ndim), work9(nodalp,ndim,ndim-1) )
      allocate ( work10(nodalp,ndim,ndim), normal(nodalp,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xrnod, x0, u )
      deallocate ( phi, dphi )
      deallocate ( x, uvec )
      deallocate ( dxdxis, surfl )
      deallocate ( dxdxis0, surfl0 )
      deallocate ( gi_up )
      deallocate ( work, work3 )
      deallocate ( tmp, work9 )
      deallocate ( work10, normal )

    end subroutine deallocate_arrays

  end subroutine interfacial_stress_scal

! Compute the eigenvalues and -vectors of a symmetric interfacial tensor
! stored in a vector. Output the eigenvalues and -vectors
! in separate vectors.

  subroutine eig ( a, lambda1, lambda2, eigv1, eigv2 )

#if NO_LIBHSL
    use eig2D3D_m, only: eig3x3
#else
    use hsl_ea23_m, only: eig3x3
#endif

    real(dp), dimension(:), intent(in) :: a
    real(dp), dimension(:), intent(out) :: lambda1, lambda2, eigv1, eigv2

    real(dp), dimension(6,size(a)/6) :: b
    real(dp) :: lambda(3), eigv(3,3)
    integer :: np, i, j, ssl(3), ssv(3), ssl1(1), ssl2(1), ssl3(1)

    np = size(a)/6
    b = reshape ( a, [6, np] )

    do i = 1, np

      call eig3x3 ( b(:,i), lambda, eigv )

      ssl = [(j, j=1,size(lambda))]
      ssv = ssl + (i - 1) * size(lambda)
      ssl1 = minloc(abs(lambda))
      ssl2 = maxloc(lambda, ssl /= ssl1(1) )
      ssl3 = pack(ssl, ssl /= ssl1(1) .and. ssl /= ssl2(1) )

      lambda1(i) = lambda(ssl2(1))
      eigv1(ssv) = eigv(:,ssl2(1))
      lambda2(i) = lambda(ssl3(1))
      eigv2(ssv) = eigv(:,ssl3(1))

    end do

  end subroutine eig

end module interfacial_rheology_elements_3D_surface_m
