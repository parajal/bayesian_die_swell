
! Copyright (C) 2005-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the Stokes equation
!
!    - div( eta(nabla u+nabla u^T) ) + nabla p = f
!      div u = 0
!
! This module contains stuff that only works in 3D or only for surfaces.
!

module stokes_elements_3D_surface_m

  use tfem_elem_m
  use stokes_set_globals_m
  use stokes_elements_generic_m, only: set_stokes_shape_function_global

  implicit none

contains


! Internal element routine for the Stokes Equation (3D)

  subroutine stokes_elem_3D ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, M, ip, i1, i2, i3, i4, vfuncnr, j, i
    real(dp) :: eta


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


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    call shape_derivative ( dphi, Finv, dphidx )

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if


!   pointers in unknowns

    i1 = ndf
    i2 = 2*ndf
    i3 = 3*ndf
    i4 = 3*ndf+ndfp

    if ( vector ) then

      vfuncnr = coefficients%i(14)

      if ( vfuncnr > 0 ) then

        do ip = 1, ninti
          fg(ip,:) = coefficients%vfunc ( ndim, vfuncnr, xg(ip,:) )
        end do

        do j = 1, ndim
          do N = 1, ndf
            tmp(N,j) = sum ( fg(:,j) * phi(:,N) * detF * wg )
          end do
        end do

        elemvec(1:i3) = reshape ( tmp, [ndim*ndf] )

        elemvec(i3+1:i4) = 0

      else

        elemvec(1:i4) = 0

      end if

    end if

    if ( matrix ) then

      eta = coefficients%r(1)

      do N = 1, ndf
        do M = N, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dphidx(ip,M,:) )
          end do
          work1 = work + dphidx(:,N,1) * dphidx(:,M,1)
          Suu(N,M) = eta * sum ( work1 * detF * wg )
          Suu(M,N) = Suu(N,M)
          work1 = work + dphidx(:,N,2) * dphidx(:,M,2)
          Svv(N,M) = eta * sum ( work1 * detF * wg )
          Svv(M,N) = Svv(N,M)
          work1 = work + dphidx(:,N,3) * dphidx(:,M,3)
          Sww(N,M) = eta * sum ( work1 * detF * wg )
          Sww(M,N) = Sww(N,M)
        end do
      end do

      do N = 1, ndf
        do M = 1, ndf
          work = dphidx(:,N,2) * dphidx(:,M,1)
          Suv(N,M) = eta * sum ( work * detF * wg )
          work = dphidx(:,N,3) * dphidx(:,M,1)
          Suw(N,M) = eta * sum ( work * detF * wg )
          work = dphidx(:,N,3) * dphidx(:,M,2)
          Svw(N,M) = eta * sum ( work * detF * wg )
        end do
      end do

      do N = 1, ndfp
        do M = 1, ndf
          Lu(N,M) = sum ( psi(:,N) * dphidx(:,M,1) * detF * wg )
          Lv(N,M) = sum ( psi(:,N) * dphidx(:,M,2) * detF * wg )
          Lw(N,M) = sum ( psi(:,N) * dphidx(:,M,3) * detF * wg )
        end do
      end do

      elemmat(    1:i1,    1:i1 ) = Suu
      elemmat(    1:i1, i1+1:i2 ) = Suv
      elemmat(    1:i1, i2+1:i3 ) = Suw
      elemmat( i1+1:i2,    1:i1 ) = transpose(Suv)
      elemmat( i1+1:i2, i1+1:i2 ) = Svv
      elemmat( i1+1:i2, i2+1:i3 ) = Svw
      elemmat( i2+1:i3,    1:i1 ) = transpose(Suw)
      elemmat( i2+1:i3, i1+1:i2 ) = transpose(Svw)
      elemmat( i2+1:i3, i2+1:i3 ) = Sww
      elemmat( i3+1:i4,    1:i1 ) = -Lu
      elemmat( i3+1:i4, i1+1:i2 ) = -Lv
      elemmat( i3+1:i4, i2+1:i3 ) = -Lw
      elemmat(    1:i1, i3+1:i4 ) = -transpose(Lu)
      elemmat( i1+1:i2, i3+1:i4 ) = -transpose(Lv)
      elemmat( i2+1:i3, i3+1:i4 ) = -transpose(Lw)
      elemmat( i3+1:i4, i3+1:i4)  = 0

!     set continuity equation (div u=0) to zero

      if ( coefficients%i(39) == 1 ) then

        do i = 1, ndim
          elemmat( posp, pos(:,i) ) = 0
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

      allocate ( tmp(ndf,ndim) )
      allocate ( work(ninti), work1(ninti) )
      allocate ( Suu(ndf,ndf), Suv(ndf,ndf), Svv(ndf,ndf) )
      allocate ( Suw(ndf,ndf), Svw(ndf,ndf), Sww(ndf,ndf) )
      allocate ( Lu(ndfp,ndf), Lv(ndfp,ndf), Lw(ndfp,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( tmp )
      deallocate ( work, work1 )
      deallocate ( Suu, Suv, Svv )
      deallocate ( Suw, Svw, Sww )
      deallocate ( Lu, Lv, Lw )

    end subroutine deallocate_arrays

  end subroutine stokes_elem_3D


! Boundary element for a natural boundary on a surface for the Stokes equation

  subroutine stokes_natboun_surface ( mesh, problem, surface, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: surface, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, ip, alt(2,3), funcnr, vfuncnr, direction, j

    alt(:,1) = [2,3]
    alt(:,2) = [1,3]
    alt(:,3) = [1,2]

    if ( first ) then

!     first element on this surface

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, surface=surface )

!     allocate arrays

      allocate ( wg(ninti), fg(ninti,ndim), surfl(ninti) )
      allocate ( tmp(ndf,ndim) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, normal )

    call isoparametric_coordinates ( x, phi, xg )

    vfuncnr = coefficients%i(15)
    funcnr = coefficients%i(16)
    direction = coefficients%i(17)

    if ( funcnr > 0 ) then
      do ip = 1, ninti
        fg(ip,direction) = coefficients%func ( funcnr, xg(ip,:) )
        fg(ip,alt(:,direction)) = 0
      end do
    else if ( vfuncnr > 0 ) then
      do ip = 1, ninti
        fg(ip,:) = coefficients%vfunc ( ndim, vfuncnr, xg(ip,:) )
      end do
    else
      write(*,'(/2(a/))') 'Error in stokes_natboun_surface:', &
        ' either funcnr > 0 or vfuncnr > 0 '
      stop
    end if

    do j = 1, ndim
      do N = 1, ndf
        tmp(N,j) = sum ( fg(:,j) * phi(:,N) * surfl * wg )
      end do
    end do

    elemvec = reshape ( tmp, [ndim*ndf] )

    if ( matrix ) then
      write(*,'(/3(a/))') 'Error in stokes_natboun_surface:', &
        ' matrix == .true. Make sure add_boundary_elements is called ', &
        ' with buildmatrix=.false.'
      stop
    end if

    if ( last ) then

!     last element on this surface

      deallocate ( wg, fg, surfl, normal )
      deallocate ( xig, phi, x )
      deallocate ( tmp )
      deallocate ( xg )
      deallocate ( dphi, dxdxis )

    end if

  end subroutine stokes_natboun_surface


! Element for the flowrate constraint on a surface

  subroutine stokes_constr_flowr_surface ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


    integer :: surface, i, j
    real(dp) :: flowrate, delta_p


    surface = problem%constraints(constr)%geometry1

!   imposed flow rate

    if ( first ) then

!     first element on this surface

      call set_globals_stokes_vp_boun ( mesh, coefficients, surface=surface, &
        maxvel3D=1 )

!     allocate arrays

      allocate ( normal(ninti,ncompu) ) ! note ncompu to allow vel3D
      allocate ( wg(ninti), surfl(ninti) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )
      allocate ( tmp(ndf,ncompu) ) ! note ncompu to allow vel3D

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, normal )

    if ( vector ) then

      if ( size(elemvecadd) == 0 ) then

!       specify flowrate

        if ( first ) then
!         first element: rhs's of equations
!         specify flowrate
          flowrate = coefficients%r(6)
          elemvec(1) = flowrate
        else
          elemvec(1) = 0
        end if

      else if ( size(elemvecadd) == 1 ) then

!       specify pressure difference

        elemvec(1) = 0

        if ( first ) then
!         first element: rhs's of equations
!         specify pressure difference
          delta_p = coefficients%r(7)
          elemvecadd(1) = delta_p
        else
          elemvecadd(1) = 0
        end if

      end if

    end if

    if ( matrix ) then

      do i = 1, ndf
        do j = 1, ncompu
          tmp(i,j) = sum ( phi(:,i) * surfl * wg * normal(:,j) )
        end do
      end do

      elemmat(1,:) = reshape ( tmp, [ndf*ncompu] )

      if ( size(elemvecadd) == 1 ) then

!       specify pressure difference

        if ( first ) then
!         first element: matrix B
          elemmatadd(1,1) = -1
        else
!         not first element: matrix B = 0
          elemmatadd(1,1) = 0
        end if

      end if

    end if

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl, normal )
      deallocate ( xig, phi, x )
      deallocate ( dphi, dxdxis )
      deallocate ( tmp )

    end if

  end subroutine stokes_constr_flowr_surface


! Element for the constraints (connection through elements)

  subroutine stokes_constr_elem_conn_surface ( mesh, problem, constr, elem, &
    node, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


    integer :: surface, i, j


!   connection through elements (weak)

    surface = problem%constraints(constr)%geometry1

    if ( first ) then

!     first element on this surface

      call set_globals_stokes_vp_boun ( mesh, coefficients, surface=surface )
      call set_globals_stokes_l_boun ( coefficients )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( psi(ninti,ndfl) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )
      allocate ( tmp(ndfl,ndf) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      call set_shape_function ( shapefuncl, xig, psi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl )

    do i = 1, ndfl
      do j = 1, ndf
        tmp(i,j) = sum ( psi(:,i) * phi(:,j) * surfl * wg )
      end do
    end do

    elemmat = 0
    elemmat(1:ndfl,1:ndf) = tmp
    elemmat(ndfl+1:2*ndfl,ndf+1:2*ndf) = tmp
    elemmat(2*ndfl+1:3*ndfl,2*ndf+1:3*ndf) = tmp

    elemmat2 = - elemmat
    elemvec  = 0

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl )
      deallocate ( xig, phi, x )
      deallocate ( psi )
      deallocate ( dphi, dxdxis )
      deallocate ( tmp )

    end if

  end subroutine stokes_constr_elem_conn_surface


! constraint for \int p dGamma=0 on a surface

  subroutine stokes_zero_pressure_int_surface ( mesh, problem, constr, elem, &
    node, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: surface, i


    surface = problem%constraints(constr)%geometry1

    if ( first ) then

!     first element on this surface

      call set_globals_stokes_vp_boun ( mesh, coefficients, surface=surface, &
        maxvel3D=1 )

      if ( ndfp /= size(elemmat,2) ) then
        write(*,'(/a/2(a,i0)/a/)') &
          'Error in stokes_zero_pressure_int_surface:', &
          ' ndfp = ', ndfp, ' whereas size(elemmat,2) = ', size(elemmat,2) , &
          ' possible cause: discontinuous pressure element '
        stop
      end if

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )
      allocate ( xg(ninti,ndim) )
      allocate ( psi(ninti,ndfp) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      call set_shape_function ( shapefuncp, xig, psi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl )

!   vector

    if ( vector ) then

      elemvec = 0

    end if

!   matrix

    if ( matrix ) then

      do i = 1, ndfp
        elemmat(1,i) = sum ( psi(:,i) * surfl * wg )
      end do

    end if

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl )
      deallocate ( xig, phi, x )
      deallocate ( dphi, dxdxis )
      deallocate ( xg )
      deallocate ( psi )

    end if

  end subroutine stokes_zero_pressure_int_surface


! constraint for \int p dGamma=0 on a volume

  subroutine stokes_zero_pressure_int_volume ( mesh, problem, constr, elem, &
    node, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: volume, i


    volume = problem%constraints(constr)%geometry1

    if ( first ) then

!     first element on this volume

      call set_globals_stokes_vp_boun ( mesh, coefficients, volume=volume )

!     allocate arrays

      allocate ( wg(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( psi(ninti,ndfp), xg(ninti,ndim) )
      allocate ( detF(ninti) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( coefficients%i(62) == 0 ) then
        call set_shape_function ( shapefuncp, xig, psi )
      end if

    end if

    call get_coordinates_geometry ( mesh, elem, x, volume=volume )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x, phi, xg )

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   vector

    if ( vector ) then

      elemvec = 0

    end if

!   matrix

    if ( matrix ) then

      do i = 1, ndfp
        elemmat(1,i) = sum ( psi(:,i) * detF * wg )
      end do

    end if

    if ( last ) then

!     last element on this volume

      deallocate ( wg )
      deallocate ( xig, phi, x )
      deallocate ( psi, xg )
      deallocate ( detF )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )

    end if

  end subroutine stokes_zero_pressure_int_volume


! Element for the flow rate (for integrate_boundary_elements in postprocessing)

  subroutine stokes_flowrate_surface ( mesh, problem, surface, elem, first, &
    last, coefficients, oldvectors, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: surface, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: ip


!   compute flow rate

    if ( first ) then

!     first element on this surface

      call set_globals_stokes_vp_boun ( mesh, coefficients, surface=surface, &
        maxvel3D=1 )

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )
      allocate ( normal(ninti,ncompu) ) ! note ncompu to allow vel3D
      allocate ( u(ncompu*ndf), tmp(ndf,ncompu), uvector(ninti,ncompu) )
      allocate ( work(ninti) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, surface=surface  )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, normal )

    call get_sysvector_geometry ( mesh, problem, oldvectors%s(1)%p, elem, u, &
      surface=surface, physq=[physqvel], layer=layer )

    tmp = reshape ( u, [ndf,ncompu] )

    uvector = matmul ( phi, tmp )

!   normal velocity

    do ip = 1, ninti
      work(ip) = dot_product ( uvector(ip,:), normal(ip,:) )
    end do

    elemvec(1) = sum ( work * surfl * wg )

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl )
      deallocate ( xig, phi, x )
      deallocate ( dphi, dxdxis )
      deallocate ( normal )
      deallocate ( u, tmp, uvector )
      deallocate ( work )

    end if

  end subroutine stokes_flowrate_surface


! Element for the drag force using integrate_boundary_elements in
! the postprocessing_m module.
!
! oldvectors:
!  v(1) = "viscous stress_tensor" derived from subroutine "stokes_stress_tensor"
!  v(2) = "pressure" derived from subroutine "stokes_pressure"

  subroutine stokes_drag_surface ( mesh, problem, surface, elem, first, last, &
    coefficients, oldvectors, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: surface, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, ip


    if ( first ) then

!     first element on this surface

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, surface=surface )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )
      allocate ( st(ndf*6), pr(ndf) )
      allocate ( work2(ndf,6), work(ninti), work4(ninti,ndim) )
      allocate ( tauten(ninti,ndim,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, normal )

!   viscous stress tensor

    call get_vector_geometry ( mesh, problem, oldvectors%v(1)%p, elem, st, &
      surface=surface, layer=layer )

    work2 = reshape ( st, [ndf,6] )

    tauten(:,1,1) = matmul ( phi, work2(:,1) )
    tauten(:,1,2) = matmul ( phi, work2(:,2) )
    tauten(:,1,3) = matmul ( phi, work2(:,3) )
    tauten(:,2,1) = tauten(:,1,2)
    tauten(:,2,2) = matmul ( phi, work2(:,4) )
    tauten(:,2,3) = matmul ( phi, work2(:,5) )
    tauten(:,3,1) = tauten(:,1,3)
    tauten(:,3,2) = tauten(:,2,3)
    tauten(:,3,3) = matmul ( phi, work2(:,6) )

!   pressure

    call get_vector_geometry ( mesh, problem, oldvectors%v(2)%p, elem, pr, &
      surface=surface, layer=layer )

    work = matmul ( phi, pr )

!   Add pressure to stress tensor

    tauten(:,1,1) = tauten(:,1,1) - work
    tauten(:,2,2) = tauten(:,2,2) - work
    tauten(:,3,3) = tauten(:,3,3) - work

!   traction force in each integration point: {t}=[tau]{n}

    do ip = 1, ninti
      work4(ip,:) = matmul( tauten(ip,:,:), normal(ip,:) )
    end do

!   integrate traction to force

    do i = 1, ndim
      elemvec(i) = sum ( work4(:,i) * surfl * wg )
    end do

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl )
      deallocate ( normal )
      deallocate ( xig, phi, x )
      deallocate ( dphi, dxdxis )
      deallocate ( st, pr )
      deallocate ( work2, work, work4 )
      deallocate ( tauten )

    end if

  end subroutine stokes_drag_surface


! Boundary element for surface tension on a surface

  subroutine surface_tension_surface ( mesh, problem, surface, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: surface, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, ip, j
    real(dp) :: gammac


    if ( first ) then

!     first element on this surface

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, surface=surface )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( tmp(ndf,ndim) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )
      allocate ( gi_up(ninti,ndim,2) )
      allocate ( work(ninti) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, &
      gi_up=gi_up )

!   surface tension coefficient

    gammac = coefficients%r(19)

    if ( matrix ) then
      write(*,'(/3(a/))') 'Error in surface_tension_surface:', &
        ' matrix == .true. Make sure add_boundary_elements is called ', &
        ' with buildmatrix=.false.'
      stop
    end if

    if ( vector ) then

!     - (nabla v)^T:Gamma*(I-nn) = - v^k_j Gamma dphi_k/dxi_i g^i_j

      do j = 1, ndim
        do N = 1, ndf
          do ip = 1, ninti
            work(ip) = dot_product ( dphi(ip,N,:), gi_up(ip,j,:) )
          end do
          tmp(N,j) = - gammac * sum ( work * surfl * wg )
        end do
      end do

      elemvec = reshape ( tmp, [ndf*ndim] )

    end if

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl )
      deallocate ( xig, phi, x )
      deallocate ( tmp )
      deallocate ( dphi, dxdxis )
      deallocate ( gi_up )
      deallocate ( work )

    end if

  end subroutine surface_tension_surface


! Boundary element for surface tension on surface that depends linearly
! on surfactant concentration

  subroutine surface_tension_surfactants ( mesh, problem, surface, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: surface, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, ip, j
    real(dp) :: beta, gammac0, cs0


    if ( first ) then

!     first element on this surface

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, surface=surface )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( tmp(ndf,ndim) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )
      allocate ( gi_up(ninti,ndim,2) )
      allocate ( work(ninti) )
      allocate ( csg(ninti), gammacg(ninti) )
      allocate ( car(ndf) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, &
      gi_up=gi_up )

!   concentration c at tn

    call get_sysvector ( oldvectors%m(1)%p, oldvectors%p(1)%p, &
      oldvectors%s(1)%p, elgrp=1, elem=elem, u=car )

    csg = matmul( phi, car )

!   beta physicochemical constant

    beta = coefficients%r(25)

!   initial surface tension

    gammac0 = coefficients%r(26)

!   initial concentration

    cs0 = coefficients%r(27)

!   surface tension as a function of surfactant concentration
!       Gamma(c) = Gamma0 / ( 1 - beta ) * ( 1 - c/c0 * beta )

    gammacg = gammac0 / ( 1._dp - beta ) * ( 1._dp - csg/cs0 * beta )

    if ( matrix ) then
      write(*,'(/3(a/))') 'Error in surface_tension_surfactants:', &
        ' matrix == .true. Make sure add_boundary_elements is called ', &
        ' with buildmatrix=.false.'
      stop
    end if

    if ( vector ) then

!     - (nabla v)^T:Gamma*(I-nn) = - v^k_j Gamma dphi_k/dxi_i g^i_j

      do j = 1, ndim
        do N = 1, ndf
          do ip = 1, ninti
            work(ip) = dot_product ( dphi(ip,N,:), gi_up(ip,j,:) )
          end do
          tmp(N,j) = - sum ( gammacg * work * surfl * wg )
        end do
      end do

      elemvec = reshape ( tmp, [ndf*ndim] )

    end if

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl )
      deallocate ( xig, phi, x )
      deallocate ( tmp )
      deallocate ( dphi, dxdxis )
      deallocate ( gi_up )
      deallocate ( work )
      deallocate ( csg, gammacg )
      deallocate ( car )

    end if

  end subroutine surface_tension_surfactants

end module stokes_elements_3D_surface_m

