
! Copyright (C) 2007-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Elements for flow with inertia.
! Use in combination with the stokes elements for the viscous and pressure
! terms or viscoelastic elements for the extra stress.
! Elements can be used in 1D, 2D and 3D.
! Also 3D velocities on a 2D domain is available.

module inertia_elements_m

  use tfem_elem_m

  implicit none

contains


! Convection term: u dot grad u.
! Picard iteration: u^{j} dot grad u^{j+1}
! For the first step a starting routine is needed

  subroutine inertia_elem_conv_picard ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use inertia_globals_m
    use stokes_set_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, k, ip
    real(dp) :: rho


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


!   Some tests

    if ( first ) then

      call check ( coefficients, 'inertia_elem_conv_picard', ncoefi=250, &
        ncoefr=200, indexarray=[48], minimum=[0], maximum=[1] )

    end if


!   map reference element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

!   axisymmetric

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

!   get the density
    rho = coefficients%r(151)

!   get the velocity at iteration j
    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, &
      physq=[physqvel], u=uvec, layer=layer )

    do i = 1, ncompu
      uvecn(:,i) = matmul( phi, uvec( pos(:,i) ) )
    end do

!   Eulerian frame or ALE

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity at tn

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, u=um, layer=layer )

      tmp = reshape ( um, [ndf,ndim] )

      uvecmeshn = matmul ( phi, tmp )

!     compute (u^{j}-um^{j}) dot grad u^{j+1} term
      do ip = 1, ninti
        ungradphi(ip,:) = &
                matmul ( dphidx(ip,:,:) , uvecn(ip,1:ndim) - uvecmeshn(ip,:) )
      end do

    else

!     Eulerian frame

!     compute u^{j} dot grad u^{j+1} term
      do ip = 1, ninti
        ungradphi(ip,:) = matmul ( dphidx(ip,:,:) , uvecn(ip,1:ndim) )
      end do

    end if

    if ( vector ) then
      elemvec = 0
    end if

    if ( matrix ) then

      elemmat = 0

      do j = 1, ndf
        do k = 1, ndf
          Suu(j,k) = rho * sum( phi(:,j) * ungradphi(:,k) * wg * detF )
        end do
      end do

!     Fill diagonal blocks

      do i = 1, ncompu
!       subblock i,i
        elemmat( pos(:,i), pos(:,i) ) = Suu
      end do

      if ( coorsys == 1 .and. vel3D == 1 ) then

!       cylindrical, axisymmetric with 3D velocities

!       Fill off-diagonal blocks

        do j = 1, ndf
          do k = 1, ndf
            Suu(j,k) = rho * &
                  sum( phi(:,j) * uvecn(:,3) / xg(:,2) * phi(:,k) * wg * detF )
          end do
        end do

!       Fill r,thetha and theta,r block

!       - v_r * u_theta^j * u_theta^j+1 / r
        elemmat( pos(:,2), pos(:,3) ) = - Suu
!       v_theta * u_theta^j * u_r^j+1 / r
        elemmat( pos(:,3), pos(:,2) ) = Suu

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

      allocate ( Suu(ndf,ndf) )
      allocate ( uvec(ndf*ncompu), um(ndf*ndim) )
      allocate ( uvecn(ninti,ncompu), uvecmeshn(ninti,ndim) )
      allocate ( ungradphi(ninti,ndf), tmp(ndf,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( Suu )
      deallocate ( uvec, um )
      deallocate ( uvecn, uvecmeshn )
      deallocate ( ungradphi, tmp )

    end subroutine deallocate_arrays

  end subroutine inertia_elem_conv_picard


! Convection term: u dot grad u.
! Newton iteration:
! u^{j} dot grad u^{j+1} + u^{j+1} dot grad u^{j} - u^{j} dot grad u^{j}
! for the first step a starting routine is needed

  subroutine inertia_elem_conv_newton ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use inertia_globals_m
    use stokes_set_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, k, l, ip
    real(dp) :: rho


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


!   Some tests

    if ( first ) then

      call check ( coefficients, 'inertia_elem_conv_newton', ncoefi=250, &
        ncoefr=200, indexarray=[48], minimum=[0], maximum=[1] )

    end if

!   map reference element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

!   axisymmetric

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

!   get density
    rho = coefficients%r(151)

!   get the velocity at t_{i}

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, &
      physq=[physqvel], u=uvec, layer=layer )

    do i = 1, ncompu
      uvecn(:,i) = matmul( phi, uvec( pos(:,i) ) )
    end do

!   dudx = velocity gradient L (only non-zero d/dx terms).

    do i = 1, ncompu
      do j = 1, ndim
        dudx(:,i,j) = matmul( dphidx(:,:,j), uvec( pos(:,i) ) )
      end do
    end do

!   Eulerian frame or ALE

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity at tn

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, u=um, layer=layer )

      tmp = reshape ( um, [ndf,ndim] )

      uvecmeshn = matmul ( phi, tmp )

!     compute (u^{j}-um^{j}) dot grad u^{j+1} term
      do ip = 1, ninti
        ungradphi(ip,:) = &
                 matmul ( dphidx(ip,:,:) , uvecn(ip,1:ndim) - uvecmeshn(ip,:) )
      end do

    else

!     Eulerian frame

!     compute u^{j} dot grad u^{j+1} term
      do ip = 1, ninti
        ungradphi(ip,:) = matmul ( dphidx(ip,:,:) , uvecn(ip,1:ndim) )
      end do

    end if

!   build vector
!   -u^{j} dot grad u^{j} and put it in the rhs (cancels minus sign!)

    if ( vector ) then

      do i = 1, ncompu
        do ip = 1, ninti
          work(ip) = dot_product ( uvecn(ip,1:ndim), dudx(ip,i,:) )
        end do
        do j = 1, ndf
          elemvec( pos(j,i) ) = rho * sum ( phi(:,j) * work * wg * detF )
        end do
      end do

      if ( coorsys == 1 .and. vel3D == 1 ) then

!       cylindrical, axisymmetric with 3D velocities

!       - v_r * u_theta^j * 2 / r

        do j = 1, ndf
          work1(j) = rho * &
                  sum( phi(:,j) * uvecn(:,3) ** 2 / xg(:,2) * wg * detF )
        end do
        elemvec( pos(:,2) ) = elemvec( pos(:,2) ) - work1

!       v_theta * u_theta^j * u_r^j / r

        do j = 1, ndf
          work1(j) = rho * &
                sum( phi(:,j) * uvecn(:,3) * uvecn(:,2) / xg(:,2) * wg * detF )
        end do
        elemvec( pos(:,3) ) = elemvec( pos(:,3) ) + work1

      end if

    end if

!   build matrix

    if ( matrix ) then

      elemmat = 0

!     off-diagonal blocks

      do i = 1, ncompu
        do l = 1, ndim
          if ( i /= l ) then
            do j = 1, ndf
              do k = 1, ndf
                Suv(j,k) = &
                    rho * sum( phi(:,j) * phi(:,k) * dudx(:,i,l) * wg * detF )
              end do
            end do
!           subblock i,l
            elemmat( pos(:,i), pos(:,l) ) = Suv
          end if
        end do
      end do

!     diagonal blocks

      do i = 1, ndim
!       v_i, u_i components
        do j = 1, ndf
          do k = 1, ndf
            Suu(j,k) = rho * sum( phi(:,j) * &
              ( ungradphi(:,k) + phi(:,k) * dudx(:,i,i) ) * wg * detF )
          end do
        end do
!       subblock i,i
        elemmat( pos(:,i), pos(:,i) ) = Suu
      end do

      if ( vel3D == 1 ) then

!       diagonal block for 3D velocities on a 2D domain

        do j = 1, ndf
          do k = 1, ndf
            Suu(j,k) = rho * sum( phi(:,j) * ungradphi(:,k) * wg * detF )
          end do
        end do
!       subblock theta,theta
        elemmat( pos(:,3), pos(:,3) ) = Suu

        if ( coorsys == 1 ) then

!         cylindrical, axisymmetric with 3D velocities

          do j = 1, ndf
            do k = 1, ndf
              Suu(j,k) = rho * &
                  sum( phi(:,j) * uvecn(:,3) / xg(:,2) * phi(:,k) * wg * detF )
            end do
          end do

          do j = 1, ndf
            do k = 1, ndf
              Sww(j,k) = rho * &
                  sum( phi(:,j) * uvecn(:,2) / xg(:,2) * phi(:,k) * wg * detF )
            end do
          end do

!         Fill r,thetha block

!         - v_r * 2 * u_theta^j * u_theta^j+1 / r
          elemmat( pos(:,2), pos(:,3) ) = &
                                      elemmat( pos(:,2), pos(:,3) ) - 2 * Suu

!         Fill theta,r and theta,theta block

!         v_theta * u_theta^j * u_r^j+1 / r
          elemmat( pos(:,3), pos(:,2) ) = elemmat( pos(:,3), pos(:,2) ) + Suu
!         v_theta * u_theta^j+1 * u_r^j / r
          elemmat( pos(:,3), pos(:,3) ) = elemmat( pos(:,3), pos(:,3) ) + Sww

        end if

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

      allocate ( Suu(ndf,ndf), Suv(ndf,ndf) )
      allocate ( uvecn(ninti,ncompu), uvecmeshn(ninti,ndim), uvec(ndf*ncompu) )
      allocate ( ungradphi(ninti,ndf), work(ninti), dudx(ninti,ncompu,ndim) )
      allocate ( tmp(ndf,ndim), um(ndf*ndim) )

      if ( coorsys == 1 .and. vel3D == 1 ) then
        allocate ( Sww(ndf,ndf), work1(ndf) )
      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( Suu, Suv )
      deallocate ( uvecn, uvecmeshn, uvec )
      deallocate ( ungradphi, work, dudx )
      deallocate ( tmp, um )

      if ( coorsys == 1 .and. vel3D == 1 ) then
        deallocate ( Sww, work1 )
      end if

    end subroutine deallocate_arrays

  end subroutine inertia_elem_conv_newton


! Unsteady term (time-derivative), discretised as:
! du/dt = (u_{n+1} - u_{n})/ Delta t
! Use additional element routines for creating the total discretized problem.
! Use the proper prefactors by setting factormat and factorvec in the heading
! of build_system() to create various high-order time integration schemes.

  subroutine inertia_elem_dudt ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m
    use stokes_set_globals_m
    use inertia_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, k
    real(dp) :: rho, deltat, fac, temp


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


!   Some tests

    if ( first ) then

      call check ( coefficients, 'inertia_elem_dudt', ncoefi=250, &
        ncoefr=200 )

    end if


!   map standard element
    call get_coordinates ( mesh, elgrp, elem, x )
    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

!   axisymmetric
    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

!   density
    rho = coefficients%r(151)

!   time step
    deltat = coefficients%r(8)
    fac = 1 / deltat

!   get the velocity vector
    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, &
      physq=[physqvel], u=uvec, layer=layer )

    do i = 1, ncompu
      uvecn(:,i) = matmul( phi, uvec( pos(:,i) ) )
    end do

    if ( vector ) then

!     rhs: rho * un / delta

      do i = 1, ncompu
        do j = 1, ndf
          elemvec ( pos(j,i) ) = &
                 rho * fac * sum ( phi(:,j) * uvecn(:,i) * wg * detF )
        end do
      end do

    end if

    if ( matrix ) then

!     Compute mass matrix

!     all matrix components in the offdiagonal blocks are zero, only
!     diagonal blocks are non-zero

      elemmat = 0

!     fill diagonal blocks

      do j = 1, ndf
        do k = j, ndf

          temp = rho * fac * sum ( phi(:,j) * phi(:,k) * wg * detF )

          do i = 1, ncompu
            elemmat( pos(j,i), pos(k,i) ) = temp
            elemmat( pos(k,i), pos(j,i) ) = temp ! symmetry
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

      allocate ( uvecn(ninti,ncompu), uvec(ndf*ncompu) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( uvecn, uvec )

    end subroutine deallocate_arrays

  end subroutine inertia_elem_dudt


! The non-linear term
!    - un dot grad un
! for a given velocity vector un
! Builds the vector (right-hand side) only.
! Use the proper prefactor by setting factorvec in the heading of build_system()
! to create higher-order time integration schemes.

  subroutine inertia_elem_ungradun ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m
    use stokes_set_globals_m
    use inertia_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, ip
    real(dp) :: rho


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


!   Some tests

    if ( first ) then

      call check ( coefficients, 'inertia_elem_ungradun', ncoefi=250, &
        ncoefr=200, indexarray=[48], minimum=[0], maximum=[1] )

    end if


!   map standard element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

!   axisymmetric
    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   get density
    rho = coefficients%r(151)

!   get the velocity vector
    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, &
      physq=[physqvel], u=uvec, layer=layer )

    do i = 1, ncompu
      uvecn(:,i) = matmul( phi, uvec( pos(:,i) ) )
    end do

!   dudx = velocity gradient L (only non-zero d/dx terms).

    do i = 1, ncompu
      do j = 1, ndim
        dudx(:,i,j) = matmul( dphidx(:,:,j), uvec( pos(:,i) ) )
      end do
    end do

!   Eulerian frame or ALE

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity at tn

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, u=um, layer=layer )

      tmp = reshape ( um, [ndf,ndim] )

      uvecmeshn = matmul ( phi, tmp )

!     compute (u^{n}-um^{n}) dot grad u^{n} term

      do ip = 1, ninti
        ungradun(ip,:) = &
                    matmul ( dudx(ip,:,:), uvecn(ip,1:ndim) - uvecmeshn(ip,:) )
      end do

    else

!     Eulerian frame

!     compute u^{n} dot grad u^{n} term

      do ip = 1, ninti
        ungradun(ip,:) = matmul ( dudx(ip,:,:), uvecn(ip,1:ndim) )
      end do

    end if


    if ( vector ) then

!     - u dot grad u

      do i = 1, ncompu
        do j = 1, ndf
          elemvec ( pos(j,i) ) = &
                        - rho * sum ( phi(:,j) * ungradun(:,i) * wg * detF )
        end do
      end do

      if ( coorsys == 1 .and. vel3D == 1 ) then

!       cylindrical, axisymmetric with 3D velocities

!       + v_r * u_theta^j * 2 / r

        do j = 1, ndf
          work1(j) = rho * &
                  sum( phi(:,j) * uvecn(:,3) ** 2 / xg(:,2) * wg * detF )
        end do
        elemvec( pos(:,2) ) = elemvec( pos(:,2) ) + work1

!       - v_theta * u_theta^j * u_r^j / r

        do j = 1, ndf
          work1(j) = rho * &
                sum( phi(:,j) * uvecn(:,3) * uvecn(:,2) / xg(:,2) * wg * detF )
        end do
        elemvec( pos(:,3) ) = elemvec( pos(:,3) ) - work1

      end if

    end if

    if ( matrix ) then

      write(*,'(/3a/)') 'Error in inertia_elem_ungradun:', &
        'Only the vector is built by this element.', &
        'Please use other routines for building the matrix.'
      stop

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

      allocate ( uvecn(ninti,ncompu), uvec(ndf*ncompu) )
      allocate ( dudx(ninti,ncompu,ndim), ungradun(ninti,ncompu) )
      allocate ( uvecmeshn(ninti,ndim) )
      allocate ( tmp(ndf,ndim), um(ndf*ndim) )

      if ( coorsys == 1 .and. vel3D == 1 ) then
        allocate ( work1(ndf) )
      end if

   end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( uvecn, uvec )
      deallocate ( dudx, ungradun )
      deallocate ( uvecmeshn )
      deallocate ( tmp, um )

      if ( coorsys == 1 .and. vel3D == 1 ) then
        deallocate ( work1 )
      end if

   end subroutine deallocate_arrays

  end subroutine inertia_elem_ungradun


! Convection term: uhat dot grad u with 2nd order prediction:
!    uhat = 2 * u_n - u_{n-1}
! For the first step a starting routine is needed.

  subroutine inertia_elem_conv_pred_2nd_order ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use inertia_globals_m
    use stokes_set_globals_m
    use limits_m, only: BLOCK_IMPROPER_ALE

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, k, ip
    real(dp) :: rho


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


!   Some tests

    if ( first ) then

      if ( coefficients%i(48) == 1 .and. BLOCK_IMPROPER_ALE ) then
        write(*,'(/6(a/))') 'Error in inertia_elem_conv_pred_2nd_order:', &
          ' Combining this routine with a moving mesh (ALE) is not ', &
          ' allowed. The effect of the moving mesh on the change ', &
          ' in nabla operator for different times within a time step ', &
          ' has not been taken into account yet.', &
          ' Set BLOCK_IMPROPER_ALE=.false. in limits_m to bypass this block.'
        stop
      end if

      call check ( coefficients, 'inertia_elem_conv_pred_2nd_order', &
        ncoefi=250, ncoefr=200, indexarray=[48], minimum=[0], &
        maximum=[1] )

    end if


!   map reference element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

!   axisymmetric
    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

!   get the density
    rho = coefficients%r(151)

!   get the velocity at iteration tn
    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, &
      physq=[physqvel], u=uvec, layer=layer )

    do i = 1, ncompu
      uvecn(:,i) = matmul( phi, uvec( pos(:,i) ) )
    end do

!   get the velocity at iteration tn-1
    call get_sysvector ( mesh, problem, oldvectors%s(2)%p, elgrp, elem, &
      physq=[physqvel], u=uvec, layer=layer )

    do i = 1, ncompu
      uvecnm1(:,i) = matmul( phi, uvec( pos(:,i) ) )
    end do


!   Eulerian frame or ALE

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity at tn

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, u=um, layer=layer )

      tmp = reshape ( um, [ndf,ndim] )

      uvecmeshn = matmul ( phi, tmp )

!     get mesh velocity at tn-1

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(2)%p, elgrp, &
        elem, u=um, layer=layer )

      tmp = reshape ( um, [ndf,ndim] )

      uvecmeshnm1 = matmul ( phi, tmp )

!     compute (2*u_n-u_{n-1}) dot grad u_{n+1} term
      do ip = 1, ninti
        uhatgradphi(ip,:) = matmul ( dphidx(ip,:,:), &
                               2 * ( uvecn(ip,1:ndim) - uvecmeshn(ip,:) ) &
                                 - ( uvecnm1(ip,1:ndim) - uvecmeshnm1(ip,:) ) )
      end do

    else

!     Eulerian frame

!     compute (2*u_n-u_{n-1}) dot grad u_{n+1} term
      do ip = 1, ninti
        uhatgradphi(ip,:) = matmul ( dphidx(ip,:,:), &
                               2 * uvecn(ip,1:ndim) - uvecnm1(ip,:) )
      end do

    end if


    if ( vector ) then
      elemvec = 0
    end if

    if ( matrix ) then

      elemmat = 0

      do j = 1, ndf
        do k =1, ndf
          Suu(j,k) = rho * sum( phi(:,j) * uhatgradphi(:,k) * wg * detF )
        end do
      end do

      do i = 1, ncompu
!       subblock i,i
        elemmat( pos(:,i), pos(:,i) ) = Suu
      end do

      if ( coorsys == 1 .and. vel3D == 1 ) then

!       cylindrical, axisymmetric with 3D velocities

!       Fill off-diagonal blocks

        do j = 1, ndf
          do k = 1, ndf
            Suu(j,k) = rho * sum( phi(:,j) * &
                     ( 2 * uvecn(:,3) - uvecnm1(:,3) ) / xg(:,2) * &
                                  phi(:,k) * wg * detF )
          end do
        end do

!       Fill r,thetha and theta,r block

!       - v_r * u_theta^j * u_theta^j+1 / r
        elemmat( pos(:,2), pos(:,3) ) = - Suu
!       v_theta * u_theta^j * u_r^j+1 / r
        elemmat( pos(:,3), pos(:,2) ) = Suu

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

      allocate ( Suu(ndf,ndf) )
      allocate ( uvecn(ninti,ncompu), uvecmeshn(ninti,ndim), uvec(ndf*ncompu) )
      allocate ( uvecnm1(ninti,ncompu), uvecmeshnm1(ninti,ndim) )
      allocate ( uhatgradphi(ninti,ndf), tmp(ndf,ndim), um(ndf*ndim) )


    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( Suu )
      deallocate ( uvecn, uvecmeshn, uvec )
      deallocate ( uvecnm1, uvecmeshnm1 )
      deallocate ( uhatgradphi, tmp, um )

    end subroutine deallocate_arrays

  end subroutine inertia_elem_conv_pred_2nd_order

end module inertia_elements_m
