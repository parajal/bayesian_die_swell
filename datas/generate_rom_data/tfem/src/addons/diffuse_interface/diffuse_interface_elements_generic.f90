
! Copyright (C) 2006-2014 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the diffuse-interface equation
!
!         dcdt + u.nabla c - div( M nabla c )  = 0
!    alpha c - beta c^3 + kappa nabla^2 c + mu = 0
!
! The solution vector consists of (c,mu).
!

module diffuse_interface_elements_generic_m

  use tfem_elem_m
  use stokes_set_globals_m

  implicit none


contains


! Internal element routine for the diffuse-interface equation
! (Classical first-order implementation).

  subroutine diffuse_interface_elem ( mesh, problem, elgrp, elem, matrix, &
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


    integer :: i, j, ip, i1, i2
    real(dp) :: alpha, beta, Mcoef, kappa
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

      call check ( coefficients, 'diffuse_interface_elem', ncoefi=200, &
        ncoefr=150, indexarray=[48,151], minimum=[0,0], maximum=[0,0]  )

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )


!   get velocity vector

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel], layer=layer )

    tmp = reshape ( u, [ndf,ndim] )

    u_n = matmul ( phi, tmp )

!   un.grad operator

    do ip = 1, ninti
      ungradphi(ip,:) = matmul ( dphidx(ip,:,:), u_n(ip,:) )
    end do

!   get c variable (concentration) at previous time step

    call get_sysvector ( mesh, problem, oldvectors%s(2)%p, elgrp, elem, c_n, &
      physq=[1], layer=layer )

    c_ng = matmul ( phi, c_n )

!   get c variable (concentration) at previous iteration

    call get_sysvector ( mesh, problem, oldvectors%s(3)%p, elgrp, elem, citer, &
      physq=[1], layer=layer )

    citer_g = matmul ( phi, citer )

!   coefficients

    deltat = coefficients%r(101)
    Mcoef  = coefficients%r(102)
    alpha  = coefficients%r(103)
    beta   = coefficients%r(104)
    kappa  = coefficients%r(105)

!   pointers in unknown vector

    i1 = ndf
    i2 = 2*ndf

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

!     left lower block matrix of GL-term in chemical potential

      work = alpha - beta * citer_g ** 2  ! Picard iteration

      do i = 1, ndf
        do j = i, ndf
          LLbm(i,j) = sum ( phi(:,i) * work * phi(:,j) * detF * wg )
          LLbm(j,i) = LLbm(i,j)
        end do
      end do

!     left lower block matrix (diffusion term)

      LLbm = LLbm - kappa * S

!     left upper block matrix (advection term)

      do i = 1, ndf
        do j = 1, ndf
          LUbm(i,j) = sum ( phi(:,i) * ungradphi(:,j) * detF * wg )
        end do
      end do

!     left upper block matrix (time derivative term)

      LUbm = LUbm + Mmat / deltat

!     Fill element matrix

      elemmat(    1:i1,    1:i1 ) = LUbm
      elemmat(    1:i1, i1+1:i2 ) = Mcoef * S
      elemmat( i1+1:i2,    1:i1 ) = LLbm
      elemmat( i1+1:i2, i1+1:i2 ) = Mmat

    end if

    if ( vector ) then

!     right-hand side of the c-equation

      do i = 1, ndf
        work1(i) = sum ( phi(:,i) * c_ng * detF * wg ) / deltat
      end do

!     Fill element vector

      elemvec(1:i1)    = work1  ! c
      elemvec(i1+1:i2) = 0      ! mu

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
      allocate ( citer(ndf), citer_g(ninti) )
      allocate ( work(ninti), work1(ndf) )
      allocate ( u_n(ninti,ndim), tmp(ndf,ndim) )
      allocate ( Mmat(ndf,ndf), S(ndf,ndf), LLbm(ndf,ndf), LUbm(ndf,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( ungradphi )
      deallocate ( u, c_n, c_ng )
      deallocate ( citer, citer_g )
      deallocate ( work, work1 )
      deallocate ( u_n, tmp )
      deallocate ( Mmat, S, LLbm, LUbm )

    end subroutine deallocate_arrays

  end subroutine diffuse_interface_elem


! Internal element routine for the diffuse-interface equation
! (first and second order, ALE and tALE).

  subroutine diffuse_interface_elem1 ( mesh, problem, elgrp, elem, matrix, &
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

    integer :: i, j, ip, i1, i2, time_int
    real(dp) :: alpha, beta, Mcoef, kappa
    real(dp) :: deltat

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

 !  allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   some tests

    if ( first ) then

      call check ( coefficients, 'diffuse_interface_elem1', ncoefi=200, &
        ncoefr=150, indexarray=[151,152], minimum=[0,0], maximum=[2,2]  )

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   coefficients time stepping

    time_int = get_coefficient ( coefficients, index=151, default=1 )
    deltat = coefficients%r(101)

!   get velocity, mu and c vectors at previous time step

    if ( coefficients%i(48) == 2 ) then

!     temporary ALE

!     compute reference and real coordinates in the temporary ALE meshes
!     NOTE: groups not yet used

      select case ( time_int )

      case (1)

!       first-order

        call set_globals_stokes_tALE ( oldvectors%m(1)%p, oldvectors%m(2)%p )

      case (2)

!       second-order

        call set_globals_stokes_tALE ( oldvectors%m(1)%p, oldvectors%m(2)%p, &
          oldvectors%m(3)%p, oldvectors%m(4)%p )

      case default

        call errormsg_case_default ( 'diffuse_interface_elem1', 'time_int', &
          int_value=time_int )

      end select

      call check ( coefficients, 'diffuse_interface_elem1', &
        indexarray=[153], minimum=[0], maximum=[1] )

!     get variables at previous time step(s) in the integration points

      do ip = 1, ninti

!       get velocity at time step n

        call get_sysvector ( oldvectors%m(2)%p, oldvectors%p(1)%p, &
           oldvectors%s(1)%p, grpelm_n(ip,1), grpelm_n(ip,2), u, &
           physq=[physqvel], layer=layer )

        tmp = reshape ( u, [ndf,ndim] )

        u_n(ip,:) = matmul ( phi_n(ip,:), tmp )

!       get c variable (concentration) at time step n

        call get_sysvector ( oldvectors%m(2)%p, oldvectors%p(2)%p, &
          oldvectors%s(2)%p, grpelm_n(ip,1), grpelm_n(ip,2), c_n, &
          physq=[1], layer=layer )

        c_ng(ip) = dot_product ( phi_n(ip,:), c_n )

        if ( time_int == 2 ) then

!         get velocity at time step n-1

          call get_sysvector ( oldvectors%m(4)%p, oldvectors%p(1)%p, &
            oldvectors%s(5)%p, grpelm_nm1(ip,1), grpelm_nm1(ip,2), u, &
            physq=[physqvel], layer=layer )

          tmp = reshape ( u, [ndf,ndim] )

          u_nm1(ip,:) = matmul ( phi_nm1(ip,:), tmp )

!         get c variable (concentration) at time step n-1

          call get_sysvector ( oldvectors%m(4)%p, oldvectors%p(2)%p, &
            oldvectors%s(4)%p, grpelm_nm1(ip,1), grpelm_nm1(ip,2), c_nm1, &
            physq=[1], layer=layer )

          c_nm1g(ip) = dot_product ( phi_nm1(ip,:), c_nm1 )

        end if

      end do

    else

!     Eulerian frame or ALE

      call check ( coefficients, 'diffuse_interface_elem1', &
        indexarray=[153], minimum=[0], maximum=[2] )

      if ( any(coefficients%i(153)==[0,1]) ) then

!       get velocity at time step n

        call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
          elgrp, elem, u, physq=[physqvel], layer=layer )

        tmp = reshape ( u, [ndf,ndim] )

        u_n = matmul ( phi, tmp )

      else if ( coefficients%i(153) == 2 ) then

!       get velocity from previous iteration step (u_i)

        call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(3)%p, &
          elgrp, elem, u, physq=[physqvel], layer=layer )

        tmp = reshape ( u, [ndf,ndim] )

        u_iter = matmul ( phi, tmp )

      end if

!     get c variable (concentration) at time step n

      call get_sysvector ( mesh, problem, oldvectors%s(2)%p, elgrp, elem, c_n, &
        physq=[1], layer=layer )

      c_ng = matmul ( phi, c_n )

      if ( time_int == 2 ) then ! higher order

        if ( any(coefficients%i(153)==[0,1]) ) then

!         get velocity at time step n-1

          call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(5)%p, &
            elgrp, elem, u, physq=[physqvel], layer=layer )

          tmp = reshape ( u, [ndf,ndim] )

          u_nm1 = matmul ( phi, tmp )

        end if

!       get c variable at time step n-1

        call get_sysvector ( mesh, problem, oldvectors%s(4)%p, elgrp, elem, &
          c_nm1, physq=[1], layer=layer )

        c_nm1g = matmul ( phi, c_nm1 )

      end if

    end if

!   get c variable (concentration) at previous iteration

    call get_sysvector ( mesh, problem, oldvectors%s(3)%p, elgrp, elem, &
      citer, physq=[1], layer=layer )

    citer_g = matmul ( phi, citer )

!   get mesh velocity for ALE or tALE

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity
      allocate( u_mesh(ninti,ndim) )
      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, u )
      tmp = reshape ( u, [ndf,ndim] )
      u_mesh = matmul ( phi, tmp )

    else if ( coefficients%i(48) == 2 ) then

!     tALE

!     get mesh velocity
      allocate( u_mesh(ninti,ndim) )
      select case ( time_int )
      case (1)
!       first-order
        u_mesh = ( xg - xg_n ) / deltat
      case (2)
!       second-order
        u_mesh = ( 1.5_dp*xg - 2.0_dp*xg_n + 0.5_dp*xg_nm1 ) / deltat
      case default
        call errormsg_case_default ( 'diffuse_interface_elem1', &
          'time_int', int_value=time_int )
      end select

    end if

!   determine ugradphi term

    if ( any(coefficients%i(48) == [1,2])  ) then ! ALE and tALE

      if ( any(coefficients%i(153)==[0,1]) ) then

!       (uhat-ugrid).grad operator

        select case ( time_int )
        case (1) ! First order Euler
          do ip = 1, ninti
            ungradphi(ip,:) = &
                      matmul ( dphidx(ip,:,:), u_n(ip,:) - u_mesh(ip,:) )
          end do
        case (2) ! Second order Gear
          do ip = 1, ninti
            ungradphi(ip,:) = &
                      matmul ( dphidx(ip,:,:), 2*u_n(ip,:) - u_nm1(ip,:) &
                               - u_mesh(ip,:) )
          end do
        case default
          call errormsg_case_default ( 'diffuse_interface_elem1', &
            'time_int', int_value=time_int )
        end select

      else if ( coefficients%i(153) == 2 ) then

!       (u_i-ugrid).grad operator

        do ip = 1, ninti
          ungradphi(ip,:) = &
                    matmul ( dphidx(ip,:,:), u_iter(ip,:) - u_mesh(ip,:) )
        end do

      end if

      deallocate( u_mesh )

    else

!     Eulerian frame

      if ( any(coefficients%i(153)==[0,1]) ) then

!       uhat.grad operator

        select case ( time_int )
        case (1) ! First order Euler
          do ip = 1, ninti
            ungradphi(ip,:) = matmul ( dphidx(ip,:,:), u_n(ip,:) )
          end do
        case (2) ! Second order Gear
          do ip = 1, ninti
            ungradphi(ip,:) = &
                      matmul ( dphidx(ip,:,:), 2*u_n(ip,:) - u_nm1(ip,:) )
          end do
        case default
          call errormsg_case_default ( 'diffuse_interface_elem1', &
            'time_int', int_value=time_int )
        end select

      else if ( coefficients%i(153) == 2 ) then

!       u_i.grad operator

        do ip = 1, ninti
          ungradphi(ip,:) = matmul ( dphidx(ip,:,:), u_iter(ip,:) )
        end do

      end if

    end if

!   material coefficients

    Mcoef  = coefficients%r(102)
    alpha  = coefficients%r(103)
    beta   = coefficients%r(104)
    kappa  = coefficients%r(105)

!   pointers in unknown vector

    i1 = ndf
    i2 = 2*ndf

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

!     left lower block matrix of GL-term in chemical potential

      if ( any ( coefficients%i(152) == [0,1] ) ) then
        work = alpha - beta * citer_g ** 2  ! Picard iteration
      else
        work = alpha - 3 * beta * citer_g ** 2  ! Newton iteration
      end if

      do i = 1, ndf
        do j = i, ndf
          LLbm(i,j) = sum ( phi(:,i) * work * phi(:,j) * detF * wg )
          LLbm(j,i) = LLbm(i,j)
        end do
      end do

!     left lower block matrix (diffusion term)

      LLbm = LLbm - kappa * S

!     left upper block matrix (advection term)

      do i = 1, ndf
        do j = 1, ndf
          LUbm(i,j) = sum ( phi(:,i) * ungradphi(:,j) * detF * wg )
        end do
      end do

!     left upper block matrix (time derivative term)

      select case ( time_int )
      case (1) ! Euler
        LUbm = LUbm + Mmat / deltat
      case (2) ! Gear
        LUbm = LUbm + (3*Mmat/2) / deltat
      case default
        call errormsg_case_default ( 'diffuse_interface_elem1', &
          'time_int', int_value=time_int )
      end select

!     Fill element matrix

      elemmat(    1:i1,    1:i1 ) = LUbm
      elemmat(    1:i1, i1+1:i2 ) = Mcoef * S
      elemmat( i1+1:i2,    1:i1 ) = LLbm
      elemmat( i1+1:i2, i1+1:i2 ) = Mmat

    end if


    if ( vector ) then

!     right-hand side of the c-equation

      select case ( time_int )
      case (1) ! Euler
        do i = 1, ndf
          work1(i) = sum ( phi(:,i) * c_ng * detF * wg ) / deltat
        end do
      case (2) ! Gear
        do i = 1, ndf
          work1(i) = sum ( phi(:,i) * ( 2*c_ng - c_nm1g/2 ) * detF * wg ) &
           / deltat
        end do
      case default
        call errormsg_case_default ( 'diffuse_interface_elem1', &
          'time_int', int_value=time_int )
      end select

!     Fill element vector

      elemvec(1:i1) = work1  ! c
      if ( any ( coefficients%i(152) == [0,1] ) ) then ! Picard iteration
        elemvec(i1+1:i2) = 0      ! mu
      else ! Newton iteration
        do i = 1, ndf
          work1(i) = - sum ( phi(:,i) * 2 * beta * (citer_g ** 3) * detF * wg )
        end do
        elemvec(i1+1:i2) = work1 ! mu
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

      allocate ( ungradphi(ninti,ndf) )
      allocate ( u(ndim*ndf), c_n(ndf), c_ng(ninti) )
      allocate ( c_nm1(ndf), c_nm1g(ninti) )
      allocate ( citer(ndf), citer_g(ninti) )
      allocate ( work(ninti), work1(ndf) )
      allocate ( tmp(ndf,ndim) )
      allocate ( Mmat(ndf,ndf), S(ndf,ndf), LLbm(ndf,ndf), LUbm(ndf,ndf) )
      if ( any(coefficients%i(153)==[0,1]) ) then
        allocate ( u_n(ninti,ndim), u_nm1(ninti,ndim) )
      else if ( coefficients%i(153) == 2 ) then
        allocate ( u_iter(ninti,ndim) )
      end if

      if ( coefficients%i(48) == 2 ) then

!       temporary ALE scheme; allocate additional arrays

        allocate ( xg_n(ninti,ndim), grpelm_n(ninti,2) )
        allocate ( xig_n(ninti,ndim), phi_n(ninti,ndf) )
        allocate ( xg_nm1(ninti,ndim), grpelm_nm1(ninti,2) )
        allocate ( xig_nm1(ninti,ndim), phi_nm1(ninti,ndf) )

        allocate ( dphi_n(ninti,ndf,ndim), F_n(ninti,ndim,ndim) )
        allocate ( Finv_n(ninti,ndim,ndim), detF_n(ninti) )

      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( ungradphi )
      deallocate ( u, c_n, c_ng )
      deallocate ( c_nm1, c_nm1g )
      deallocate ( citer, citer_g )
      deallocate ( work, work1 )
      deallocate ( tmp )
      deallocate ( Mmat, S, LLbm, LUbm )
      if ( any(coefficients%i(153)==[0,1]) ) then
        deallocate ( u_n, u_nm1 )
      else if ( coefficients%i(153) == 2 ) then
        deallocate ( u_iter )
      end if

      if ( coefficients%i(48) == 2 ) then

!       temporary ALE scheme; deallocate additional arrays

        deallocate ( xg_n, grpelm_n )
        deallocate ( xig_n, phi_n )
        deallocate ( xg_nm1, grpelm_nm1 )
        deallocate ( xig_nm1, phi_nm1 )

        deallocate ( dphi_n, F_n, Finv_n, detF_n )

      end if

    end subroutine deallocate_arrays

  end subroutine diffuse_interface_elem1


! Internal element routine for the rhs in the momentum balance

  subroutine rhs_mugradc ( mesh, problem, elgrp, elem, matrix, &
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

      call check ( coefficients, 'diffuse_mugradc', ncoefi=200, &
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


!   get c variable (concentration)

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, c_n, physq=[1], layer=layer )

!   grad c term

    do i = 1, ndim
      gradc(:,i) = matmul ( dphidx(:,:,i), c_n )
    end do

!   get mu variable (chemical potential))

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, mu, physq=[2], layer=layer )

    mu_g = matmul ( phi, mu )


!   coefficients

    rho  = coefficients%r(106)

    if ( vector ) then

!     term rho mu * gradc

      do j = 1, ndim
        do i = 1, ndf
          tmp(i,j) = rho * sum ( phi(:,i) * mu_g * gradc(:,j) * detF * wg )
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

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( tmp )
      deallocate ( c_n, gradc )
      deallocate ( mu, mu_g )

    end subroutine deallocate_arrays

  end subroutine rhs_mugradc


! Element routine for the divergence of the capillary stress in the momentum
! balance

  subroutine rhs_divstressdi ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, N, ip
    real(dp) :: kappa, rho


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )


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

      call check ( coefficients, 'rhs_divstressdi', ncoefi=200, &
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


!   coefficients

    kappa  = coefficients%r(105)
    rho    = coefficients%r(106)


!   get c variable (concentration)

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, c_n, physq=[1], layer=layer )


!   grad c term

    do i = 1, ndim
      gradc(:,i) = matmul ( dphidx(:,:,i), c_n )
    end do


!   the -gradcgradc tensor

    do i = 1, ndim
      do j = i, ndim
        gradcgradc(:,i,j) = - gradc(:,i) * gradc(:,j)
        if ( i /= j ) gradcgradc(:,j,i) = gradcgradc(:,i,j)
      end do
    end do


!   the isotropic part |gradc|^2

    gradc_sq = sum ( gradc**2._dp, 2 )


!   combine the gradcgradc and |gradc|^2 I tensors to form the capillary
!   stress tensor

    tau_c = rho * kappa * gradcgradc
    do i = 1, ndim
      tau_c(:,i,i) = tau_c(:,i,i) + rho * kappa * gradc_sq
    end do


!   build equations

    if ( matrix ) then

      elemmat = 0

    end if

    if ( vector ) then

!     - (nabla v)^T:stress_di

      do ip = 1, ninti
        work6(ip,:,:) = matmul ( dphidx(ip,:,:), tau_c(ip,:,:) )
      end do
      if ( coorsys == 1 ) then
        do N = 1, ndf
          work6(:,N,2) = work6(:,N,2) + &
                                 rho * kappa * gradc_sq * phi(:,N) / xg(:,2)
        end do
      end if
      do j = 1, ndim
        work2(:,j) = matmul ( detF * wg, work6(:,:,j) )
      end do
      elemvec = - reshape ( work2, [ ndf*ndim ] )

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

      allocate ( work6(ninti,ndf,ndim), work2(ndf,ndim) )
      allocate ( c_n(ndf), gradc(ninti,ndim) )
      allocate ( gradcgradc(ninti,ndim,ndim) )
      allocate ( tau_c(ninti,ndim,ndim) )
      allocate ( gradc_sq(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work6, work2 )
      deallocate ( c_n, gradc )
      deallocate ( gradcgradc )
      deallocate ( tau_c )
      deallocate ( gradc_sq )

    end subroutine deallocate_arrays

  end subroutine rhs_divstressdi


! element routine for the -cgradmu term on the right hand side of the momemtum
! balance

  subroutine rhs_mincgradmu ( mesh, problem, elgrp, elem, matrix, &
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

      call check ( coefficients, 'rhs_mincgradmu', ncoefi=200, &
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


!   get c variable (concentration)

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, c_n, physq=[1], layer=layer )

    c_ng = matmul ( phi, c_n )


!   get mu variable (chemical potential))

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, mu, physq=[2], layer=layer )


!   grad mu term

    do i = 1, ndim
      gradmu(:,i) = matmul ( dphidx(:,:,i), mu )
    end do


!   coefficients

    rho  = coefficients%r(106)

    if ( vector ) then

!     term rho -c * grad mu

      do j = 1, ndim
        do i = 1, ndf
          tmp(i,j) = - rho * sum ( phi(:,i) * c_ng * gradmu(:,j) * detF * wg )
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
      allocate ( c_n(ndf), gradmu(ninti,ndim) )
      allocate ( mu(ndf), c_ng(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( tmp )
      deallocate ( c_n, gradmu )
      deallocate ( mu, c_ng )

    end subroutine deallocate_arrays

  end subroutine rhs_mincgradmu


! subroutine to sample the (c,mu)

  subroutine diffuse_interface_sample ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, cmuvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: cmuvec


    integer :: elgrp, elem
    real(dp) :: xr(1,mesh%ndim)


    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

    allocate ( phi(1,ndf), u(2*ndf), tmp(ndf,2) )

!   set shape function in the point

    call set_shape_function ( shapefunc, xr, phi )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      layer=layer )

    tmp = reshape ( u, [ndf,2] )

    cmuvec = matmul( phi(1,:), tmp )

    deallocate ( phi, u, tmp )

  end subroutine diffuse_interface_sample


! subroutine to sample the concentration (c)

  subroutine diffuse_interface_sample_c ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, cval )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: cval


    integer :: elgrp, elem
    real(dp) :: xr(1,mesh%ndim)


    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

    allocate ( phi(1,ndf), u(ndf) )

!   set shape function in the point

    call set_shape_function ( shapefunc, xr, phi )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[1], layer=layer )

    cval = matmul( phi, u )

    deallocate ( phi, u )

  end subroutine diffuse_interface_sample_c


! sample the chemical potential (mu)

  subroutine diffuse_interface_sample_mu ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, muval )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: muval


    integer :: elgrp, elem
    real(dp) :: xr(1,mesh%ndim)

    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

    allocate ( phi(1,ndf), u(ndf) )

!   set shape function in the point

    call set_shape_function ( shapefunc, xr, phi )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[2], layer=layer )

    muval = matmul( phi, u )

    deallocate ( phi, u )

  end subroutine diffuse_interface_sample_mu


! Boundary element for a natural boundary on a curve for the DIM

  subroutine di_natboun_curve ( mesh, problem, curve, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j
    real(dp) :: phi_di


    if ( first ) then

!     first element on this curve

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, curve=curve )

!     allocate arrays

      allocate ( wg(ninti), fg(ninti,ndim), curvel(ninti) )
      allocate ( tmp(ndf,ndim) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,1), dxdxi(ninti,ndim) )
      allocate ( citer(ndf), citer_g(ninti) )
      allocate ( emat(ndf,ndf), evec(ndf) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )

    call isoparametric_deformation_curve ( x, dphi(:,:,1), dxdxi, curvel, &
      normal )

    call isoparametric_coordinates ( x, phi, xg )

    if ( coorsys == 1 ) then
      curvel = 2 * pi * xg(:,2) * curvel
    end if

!   get c variable (concentration) at previous iteration

    call get_sysvector_geometry ( mesh, problem, oldvectors%s(3)%p, elem, &
      citer, curve=curve, physq=[1], layer=layer )

    citer_g = matmul ( phi, citer )

!   get the wetting potential

    phi_di = coefficients%r(108)

!   fill elemvec

    elemvec= 0

    if ( any ( coefficients%i(152) == [0,1] ) ) then ! Picard iteration
      do i = 1,ndf
       evec(i) = phi_di * sum ( phi(:,i) * curvel * wg )
      end do
    else ! Newton iteration
      do i = 1,ndf
        evec(i) = ( phi_di * sum ( phi(:,i) * curvel * wg ) + &
                    phi_di * sum ( (citer_g**2) * phi(:,i) * curvel * wg ) )
      end do
    end if

    elemvec(1+ndf:2*ndf) = evec

    if ( matrix ) then

!     Robin boundary condition

    elemmat = 0

    if ( any ( coefficients%i(152) == [0,1] ) ) then ! Picard iteration
      do i = 1, ndf
        do j = i, ndf
          emat(i,j) =  phi_di * sum ( citer_g * phi(:,i) * phi(:,j) * &
            curvel * wg )
          emat(j,i) = emat(i,j)  ! symmetry
        end do
      end do
    else ! Newton iteration
      do i = 1, ndf
        do j = i, ndf
          emat(i,j) = 2 * phi_di * sum ( citer_g * phi(:,i) * phi(:,j) * &
            curvel * wg )
          emat(j,i) = emat(i,j)  ! symmetry
        end do
      end do
    end if

    elemmat ( 1+ndf:2*ndf, 1:ndf) = emat

    end if

    if ( last ) then

!     last element on this curve

      deallocate ( wg, fg, curvel )
      deallocate ( tmp )
      deallocate ( normal )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi, dxdxi )
      deallocate ( citer, citer_g )
      deallocate ( emat, evec )

    end if

  end subroutine di_natboun_curve


! Boundary element for a natural boundary on a surface for the DIM

  subroutine di_natboun_surface ( mesh, problem, surface, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: surface, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j
    real(dp) :: phi_di


    if ( first ) then

!     first element on this curve

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, surface=surface )

!     allocate arrays

      allocate ( wg(ninti), fg(ninti,ndim), surfl(ninti) )
      allocate ( tmp(ndf,ndim) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )
      allocate ( citer(ndf), citer_g(ninti) )
      allocate ( emat(ndf,ndf), evec(ndf) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, &
      normal )

    call isoparametric_coordinates ( x, phi, xg )

    if ( coorsys == 1 ) then
      curvel = 2 * pi * xg(:,2) * curvel
    end if

!   get c variable (concentration) at previous iteration

    call get_sysvector_geometry ( mesh, problem, oldvectors%s(3)%p, elem, &
      citer, surface=surface, physq=[1], layer=layer )

    citer_g = matmul ( phi, citer )

!   get the wetting potential

    phi_di = coefficients%r(108)

!   fill elemvec

    elemvec= 0

    if ( any ( coefficients%i(152) == [0,1] ) ) then ! Picard iteration
      do i = 1,ndf
       evec(i) = phi_di * sum ( phi(:,i) * surfl * wg )
      end do
    else ! Newton iteration
      do i = 1,ndf
        evec(i) = ( phi_di * sum ( phi(:,i) * surfl * wg ) + &
                    phi_di * sum ( (citer_g**2) * phi(:,i) * surfl * wg ) )
      end do
    end if

    elemvec(1+ndf:2*ndf) = evec

    if ( matrix ) then

!     Robin boundary condition

    elemmat = 0

    if ( any ( coefficients%i(152) == [0,1] ) ) then ! Picard iteration
      do i = 1, ndf
        do j = i, ndf
          emat(i,j) =  phi_di * sum ( citer_g * phi(:,i) * phi(:,j) * &
            surfl * wg )
          emat(j,i) = emat(i,j)  ! symmetry
        end do
      end do
    else ! Newton iteration
      do i = 1, ndf
        do j = i, ndf
          emat(i,j) = 2 * phi_di * sum ( citer_g * phi(:,i) * phi(:,j) * &
            surfl * wg )
          emat(j,i) = emat(i,j)  ! symmetry
        end do
      end do
    end if

    elemmat ( 1+ndf:2*ndf, 1:ndf) = emat

    end if

    if ( last ) then

!     last element on this curve

      deallocate ( wg, fg, surfl )
      deallocate ( tmp )
      deallocate ( normal )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi, dxdxis )
      deallocate ( citer, citer_g )
      deallocate ( emat, evec )

    end if

  end subroutine di_natboun_surface



end module diffuse_interface_elements_generic_m
