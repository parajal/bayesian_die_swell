
! Copyright (C) 2005-2025 Martien A. Hulsen
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
! For constant eta, these equations reduce to (Laplace form):
!
!    - div( eta nabla u ) + nabla p = f
!      div u = 0
!
! This module contains the generic stuff (works for 2D and 3D).

module stokes_elements_generic_m

  use tfem_elem_m

  implicit none

contains


! Internal element routine for the mass matrix without density.

  subroutine stokes_mass_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, M, i
    real(dp) :: fac


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors, maxvel3D=1 )

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if


    if ( vector ) then

      elemvec = 0

    end if


    if ( matrix ) then

!     fill matrix

      elemmat = 0

!     fill diagonal blocks

      do N = 1, ndf
        do M = N, ndf

          fac = sum ( phi(:,N) * phi(:,M) * wg * detF )

          do i = 1, ncompu
            elemmat( pos(N,i), pos(M,i) ) = fac
            elemmat( pos(M,i), pos(N,i) ) = fac ! symmetry
          end do

        end do
      end do

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )


  end subroutine stokes_mass_elem


! Internal element routine for the Stokes equation in Laplace form.

  subroutine stokes_Laplace_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, M, i, j, ip, vfuncnr
    real(dp) :: S_NMii, L_NMi
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


    if ( vector ) then

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

!     viscosity constant

      eta = coefficients%r(1)

!     fill matrix

      elemmat = 0

!     fill diagonal blocks

      do N = 1, ndf
        do M = N, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dphidx(ip,M,:) )
          end do
          do i = 1, ndim
            work1 = work
            if ( coorsys == 1 .and. i == 2 ) then ! axisymmetric
              work1 = work1 + 2 * phi(:,N) * phi(:,M) / xg(:,2) ** 2
            end if
            S_NMii = eta * sum ( work1 * detF * wg )
            elemmat( pos(N,i), pos(M,i) ) = S_NMii
            elemmat( pos(M,i), pos(N,i) ) = S_NMii  ! = S_MNii symmetry
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

      allocate ( work(ninti), work1(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work, work1 )

    end subroutine deallocate_arrays

  end subroutine stokes_Laplace_elem


! Internal element routine for the right-hand side of the momentum balance:
!
!   ... = div sigma
!
! with sigma = -p I + 2 eta D, with the velocity and pressure known vectors.

  subroutine stokes_rhs_divsigma ( mesh, problem, elgrp, elem, matrix, vector, &
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


    real(dp) :: eta
    integer :: i, j, N, ip


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


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if


!   get velocity and compute gradients

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( vel3D == 1 ) then

!     3D velocities on a 2D domain

      if ( coorsys == 0 ) then
!       Cartesian coordinates
        gradu(:,:,3) = 0
      else if ( coorsys == 1 ) then
        ugvector(:,2:3) =  matmul ( phi, uvector(:,2:3) )
        gradu(:,1,3) = 0
        gradu(:,2,3) = - ugvector(:,3) / xg(:,2) ! - u_theta / r
        gradu(:,3,3) = ugvector(:,2) / xg(:,2)   ! u_r / r
      end if

    end if


!   get pressure

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, pr, &
      physq=[physqpress], layer=layer )

    press = matmul ( psi, pr )


!   compute Cauchy stress tensor sigma

    eta = coefficients%r(1)

    do i = 1, ncompu
!     diagonal
      sigmaten(:,i,i) = - press + 2 * eta * gradu(:,i,i)
      do j = i+1, ncompu
!       off-diagonal
        sigmaten(:,i,j) = eta * ( gradu(:,i,j) + gradu(:,j,i) )
        sigmaten(:,j,i) = sigmaten(:,i,j) ! symmetry
      end do
    end do

    if ( coorsys == 1 .and. vel3d == 0 ) then
!     axisymmetric
      ugvector(:,2) =  matmul ( phi, uvector(:,2) )
      sigmatt = - press + 2 * eta * ugvector(:,2) / xg(:,2) ! -p + 2 eta u_r / r
    end if


!   build equations

    if ( matrix ) then

      write(*,'(/a/a/)') 'Error in rhs_divsigma: no matrix to build.', &
        'Call build_system with buildmatrix=.false.'
      stop

    end if

    if ( vector ) then

!     - (nabla v)^T:sigma

      do ip = 1, ninti
        work6(ip,:,:) = matmul ( dphidx(ip,:,:), sigmaten(ip,1:ndim,:) )
      end do
      if ( coorsys == 1 .and. vel3D == 0 ) then
        do N = 1, ndf
          work6(:,N,2) = work6(:,N,2) + sigmatt * phi(:,N) / xg(:,2)
        end do
      else if ( coorsys == 1 .and. vel3D == 1 ) then
        do N = 1, ndf
!         v_r/r * sigma_theta,theta
          work6(:,N,2) = work6(:,N,2) + sigmaten(:,3,3) * phi(:,N) / xg(:,2)
!         - v_theta/r * sigma_theta,r
          work6(:,N,3) = work6(:,N,3) - sigmaten(:,3,2) * phi(:,N) / xg(:,2)
        end do
      end if
      do j = 1, ncompu
        work2(:,j) = - matmul ( detF * wg, work6(:,:,j) )
      end do
      elemvec = reshape ( work2, [ ndf*ncompu ] )

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

      allocate ( work6(ninti,ndf,ncompu), work2(ndf,ncompu) )
      allocate ( u(ncompu*ndf), uvector(ndf,ncompu) )
      allocate ( gradu(ninti,ncompu,ncompu) )
      allocate ( pr(ndfp), press(ninti) )
      allocate ( sigmatt(ninti), sigmaten(ninti,ncompu,ncompu) )
      allocate ( ugvector(ninti,ncompu) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work6, work2 )
      deallocate ( u, uvector, gradu )
      deallocate ( pr, press )
      deallocate ( sigmatt, sigmaten )
      deallocate ( ugvector )

    end subroutine deallocate_arrays

  end subroutine stokes_rhs_divsigma


! Internal element routine for the right-hand side (q,divu)

  subroutine stokes_rhs_divu ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j


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

    call shape_derivative ( dphi, Finv, dphidx )

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   get velocity vector

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    if ( coorsys == 1 ) ugvector(:,2) =  matmul ( phi, uvector(:,2) )

    do j = 1, ndim
      gradu(:,j,j) = matmul ( dphidx(:,:,j), uvector(:,j) )
    end do

    if ( vector ) then

!     nabla.u

      work = 0
      do i = 1, ndim
        work = work + gradu(:,i,i)
      end do
      if ( coorsys == 1 ) then
        work = work + ugvector(:,2) / xg(:,2)
      end if

!     q nabla.u

      elemvec = matmul ( work * detF * wg, psi )

    end if

    if ( matrix ) then

      write(*,'(/a/a/)') 'Error in stokes_rhs_divu:', &
        ' No matrix to build. Call build_system with buildmatrix=.false.'
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

      allocate ( work(ninti) )
      allocate ( u(ndf*ncompu), uvector(ndf,ncompu) )
      allocate ( gradu(ninti,ncompu,ncompu) )
      if ( coorsys == 1 ) allocate ( ugvector(ninti,ncompu) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work )
      deallocate ( u, uvector, gradu )
      if ( coorsys == 1 ) deallocate ( ugvector )

    end subroutine deallocate_arrays

  end subroutine stokes_rhs_divu


! Element for the constraints (connection through collocation)

  subroutine stokes_constr_node_conn ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


    integer :: i

!   connection through collocation

    elemmat  = 0
    do i = 1, size(elemmat,1)
      elemmat(i,i) = 1
    end do
    elemmat2 = - elemmat
    elemvec  = 0

  end subroutine stokes_constr_node_conn


! Element for constraints with zero matrix (connection through collocation)

  subroutine zero_constr_node_conn ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

!   connection through collocation

    elemmat  = 0
    elemmat2 = 0
    elemvec  = 0

  end subroutine zero_constr_node_conn


! open boundary element using an object

  subroutine stokes_open_boundary ( mesh, problem, eleminfo, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(eleminfo_t), intent(in) :: eleminfo
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   variables

    integer :: elgrp, elem, object, elemo
    integer :: i, j, N, M, ip, orientation
    real(dp) :: eta


!   set globals fluid element

    elgrp  = eleminfo%elgrp

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals object

    object = eleminfo%object

    call set_globals_stokes_object ( mesh, coefficients, object )

    elem   = eleminfo%elem
    elemo  = eleminfo%elemo
    ninti  = eleminfo%nintps  ! reset ninti to the number in nintps

!   allocate arrays

    allocate ( pos(ndf,ndim), posp(ndfp) )
    allocate ( xig(ninti,ndim), x(nodalpo,ndim), xel(ndf,ndim) )
    allocate ( wg(ninti), curvel(ninti), normal(ninti,ndim) )
    allocate ( phi(ninti,ndf), dxdxi(ninti,ndim), psi(ninti,ndfp) )
    allocate ( detF(ninti), dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
    allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
    allocate ( theta(ninti,nodalpo), dtheta(ninti,nodalpo,ndim-1) )
    allocate ( work(ninti), work1(ninti), work2(ninti,ndim*ndf) )
    allocate ( xigo(ninti,ndim-1), dxdxis(ninti,ndim,2) )
    allocate ( xg(ninti,ndim) )

!   position arrays
    pos  = reshape ( [ ( i, i = 1, ndf*ndim ) ], [ ndf, ndim ] )
    posp = [ ( i, i = 1, ndfp ) ] + ndf*ndim

!   integration points and weighting factors

    xig  = mesh%objects(object)%refcoor_int(eleminfo%intps,:,elemo)
    xigo = mesh%objects(object)%xig(eleminfo%intps,:)
    wg   = mesh%objects(object)%wg(eleminfo%intps)

!   set shape functions

    call set_shape_function ( shapefunc, xig, phi, dphi ) ! velocity
    if ( coefficients%i(62) == 0 ) then
      call set_shape_function ( shapefuncp, xig, psi ) ! pressure
    end if
    call set_shape_function ( shapefunco, xigo, theta, dtheta ) ! shape object

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, xel )

    call isoparametric_deformation ( xel, dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

!   compute deformed object element

    call get_coordinates_object ( mesh, elemo, x, object )

    xg =  matmul ( theta, x )

    if ( ndim == 2 ) then
      call isoparametric_deformation_curve ( x, dtheta(:,:,1), dxdxi, curvel, &
        normal )
      if ( coorsys == 1 ) then
        curvel = 2 * pi * xg(:,2) * curvel
      end if
    else if ( ndim == 3 ) then
      call isoparametric_deformation_surface ( x, dtheta, dxdxis, curvel, &
        normal )
    end if

!   global shape function for pressure

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, xel, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   set normal

    call check ( coefficients, 'stokes_open_boundary', indexarray=[36], &
      minimum=[-1], maximum=[1] )

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    normal = real ( orientation, kind=dp ) * normal

!   vector

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     viscous traction term

      eta = coefficients%r(1)

!     diagonal blocks

      do N = 1, ndf
        do M = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,M,:) * normal(ip,:) )
          end do
          do i = 1, ndim
            work1 = phi(:,N) * ( work + dphidx(:,M,i) * normal(:,i) )
            elemmat( pos(N,i), pos(M,i) ) = - eta * sum ( work1 * curvel * wg )
          end do
        end do
      end do

!     off-diagonal blocks

      do N = 1, ndf
        do M = 1, ndf
          do i = 1, ndim
            do j = 1, ndim
              if ( i == j ) cycle
              work = phi(:,N) * dphidx(:,M,i) * normal(:,j)
              elemmat( pos(N,i), pos(M,j) ) = - eta * sum ( work * curvel * wg )
            end do
          end do
        end do
      end do

!     velocity-pressure part

      do N = 1, ndf
        do M = 1, ndfp
          do i = 1, ndim
            elemmat( pos(N,i), posp(M) ) = &
                        sum ( phi(:,N) * psi(:,M) * normal(:,i) * curvel * wg )
          end do
        end do
      end do

!     pressure-pressure part

      elemmat( posp, posp ) = 0

!     pressure-velocity part

      do i = 1, ndim
        elemmat( posp, pos(:,i) ) = 0
      end do

    end if

!   deallocate arrays

    deallocate ( pos, posp )
    deallocate ( xig, x, xel )
    deallocate ( wg, curvel, normal )
    deallocate ( phi, dxdxi, psi )
    deallocate ( detF, dphi, F )
    deallocate ( Finv, dphidx )
    deallocate ( theta, dtheta )
    deallocate ( work, work1, work2 )
    deallocate ( xigo, dxdxis )
    deallocate ( xg )

  end subroutine stokes_open_boundary


! constraint for \int p dGamma=0 on an object

  subroutine stokes_zero_pressure_int_object ( mesh, problem, constr, elem, &
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

    integer :: object, i, elgrp, elemnr


    object = problem%constraints(constr)%object

!   element group number and element number in the fluid

    elgrp = mesh%objects(object)%grpelm_int(node,1,elem)
    elemnr = mesh%objects(object)%grpelm_int(node,2,elem)

!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

!   set globals object

    call set_globals_stokes_object ( mesh, coefficients, object )

!   allocate arrays

    allocate ( xig(1,ndim), x(nodalpo,ndim), xel(ndf,ndim) )
    allocate ( wg(1), curvel(1) )
    allocate ( phi(1,ndf), dxdxi(1,ndim), psi(1,ndfp) )
    allocate ( detF(1), dphi(1,ndf,ndim), F(1,ndim,ndim) )
    allocate ( Finv(1,ndim,ndim), dphidx(1,ndf,ndim) )
    allocate ( theta(1,nodalpo), dtheta(1,nodalpo,ndim-1) )
    allocate ( xigo(1,ndim-1), dxdxis(1,ndim,2) )

!   reference coordinates of the integration point (node) of the element (elem)

    xig(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)

!   set Gauss integration and shape function

    call set_shape_function ( shapefunc, xig, phi, dphi ) ! velocity
    if ( coefficients%i(62) == 0 ) then
      call set_shape_function ( shapefuncp, xig, psi ) ! pressure
    else if ( any( coefficients%i(62) == [1,2] ) ) then
      write(*,'(2a)') ' Error in stokes_zero_pressure_int_object: ', &
                      'global interpolation not implemented'
    end if
    call set_shape_function ( shapefunco, xig, theta, dtheta ) ! shape object

!   integration rule

    xigo(1,:) = mesh%objects(object)%xig(node,:)
    wg = mesh%objects(object)%wg(node)

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elemnr, xel )

    call isoparametric_deformation ( xel, dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

!   compute deformed object element

    call get_coordinates_object ( mesh, elem, x, object )

    if ( shapefunco%globalshape == 'line' ) then
      call isoparametric_deformation_curve ( x, dtheta(:,:,1), dxdxi, curvel )
    else if ( shapefunco%globalshape == 'triangle' .or. &
              shapefunco%globalshape == 'quadrilateral' ) then
      call isoparametric_deformation_surface ( x, dtheta, dxdxis, curvel )
    else
      write(*,'(/a/2a/)') 'Error in stokes_zero_pressure_int_object: ', &
        ' invalid globalshape: ', shapefunco%globalshape
      stop
    end if

!   vector

    if ( vector ) then

      elemvec = 0

    end if

!   matrix

    if ( matrix ) then

      do i = 1, ndfp
        elemmat(1,i) = psi(1,i) * curvel(1) * wg(1)
      end do

    end if

!   deallocate arrays

    deallocate ( xig, x, xel )
    deallocate ( wg, curvel )
    deallocate ( phi, dxdxi, psi )
    deallocate ( detF, dphi, F )
    deallocate ( Finv, dphidx )
    deallocate ( theta, dtheta )
    deallocate ( xigo, dxdxis )

  end subroutine stokes_zero_pressure_int_object


! constraint for \int p dGamma=0 on an elementset

  subroutine stokes_zero_pressure_int_elementset ( mesh, problem, constr, &
    elem, node, matrix, vector, first, last, coefficients, oldvectors, &
    elemmat, elemmat2, elemmatadd, elemvec, elemvecadd )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: elgrp, i


    elgrp = node ! node is used to transfer element group number

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors, maxvel3D=1 )

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

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

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

  end subroutine stokes_zero_pressure_int_elementset


! compute velocity in nodal points (when some nodes have zero velocity)

  subroutine stokes_velocity ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

!     allocate arrays

      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf) )
      allocate ( u(ncompu*ndf), tmp(ndf,ncompu) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi )

    end if

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    tmp = reshape ( u, [ndf,ncompu] )

    elemvec = reshape ( matmul ( phi, tmp ), [nodalp*ncompu] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( xrnod, phi )
      deallocate ( u, tmp )

    end if

  end subroutine stokes_velocity


! compute velocity derivatives, vorticity, divergence, ...

  subroutine stokes_deriv ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: comp, j, maxderiv, ip


    if ( first ) then

!     first element in this group

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

!     check component = coefficients%i(13)

      if ( coorsys <= 1 .and. vel3D == 0 ) then
        maxderiv = 8
      else
        maxderiv = 11
      end if

      call check ( coefficients, 'stokes_deriv', indexarray=[13], &
        minimum=[1], maximum=[maxderiv] )

!     allocate arrays

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ncompu*ndf), uvector(ndf,ncompu) )
      allocate ( gradu(nodalp,ncompu,ncompu) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( vel3D == 1 ) then

      gradu(:,:,3) = 0

      if ( coorsys == 1 ) then
        where ( x(:,2) < 1e-10_dp )
          gradu(:,2,3) = -gradu(:,3,2)       ! -du_theta / dr (r=0)
          gradu(:,3,3) = gradu(:,2,2)        ! du_r / dr (r=0)
        elsewhere
          gradu(:,2,3) = -uvector(:,3) / x(:,2)   ! -u_theta / r
          gradu(:,3,3) = uvector(:,2) / x(:,2)   ! u_r / r
        end where
      end if

    end if

    comp = coefficients%i(13)

    if ( coorsys <= 1 .and. vel3D == 0 ) then

!     2D

      select case(comp)
      case(1)
        elemvec = gradu(:,1,1) ! dudx
      case(2)
        elemvec = gradu(:,1,2) ! dudy
      case(3)
        elemvec = gradu(:,2,1) ! dvdx
      case(4)
        elemvec = gradu(:,2,2) ! dvdy
      case(5)
        elemvec = gradu(:,2,1) - gradu(:,1,2) ! dvdx - dudy = vorticity
      case(6)
        if ( coorsys == 0 ) then
          write(*,'(/2(a/))') 'Error in stokes_deriv:', &
            ' gradu in theta direction not applicable for coorsys=0 '
          stop
        end if
        where ( x(:,2) < 1e-10_dp )
          elemvec = gradu(:,2,2)             ! du_r / dr (r=0)
        elsewhere
          elemvec = uvector(:,2) / x(:,2)    ! u_r / r
        end where
      case(7)
        elemvec = gradu(:,1,1) + gradu(:,2,2) ! dudx + dudy = divergence
        if ( coorsys == 1 ) then ! axisymmetric
          where ( x(:,2) < 1e-10_dp )
            elemvec = elemvec + gradu(:,2,2)             ! du_r / dr (r=0)
          elsewhere
            elemvec = elemvec + uvector(:,2) / x(:,2)    ! u_r / r
          end where
        end if
      case(8) ! effective strain rate = sqrt(2II_D)
        do ip = 1, nodalp
!         second invariant
          elemvec(ip) = &
               sum ( ( gradu(ip,:,:) + transpose(gradu(ip,:,:)) ) ** 2 ) / 2
          if ( coorsys == 1 ) then
!           axisymmetric
            if ( x(ip,2) < 1e-10_dp ) then
!             use du_r / dr (r=0) for gradu_tt
              elemvec(ip) = elemvec(ip) + 2 * ( gradu(ip,2,2) ) ** 2
            else
!             use u_r / r for gradu_tt
              elemvec(ip) = elemvec(ip) + 2 * ( uvector(ip,2) / x(ip,2) ) ** 2
            end if
          end if
        end do
        elemvec = sqrt ( elemvec )
      case default
        call errormsg_case_default ( 'stokes_deriv', &
          'comp', int_value=comp )
      end select

    else if ( coorsys == 2 .or. vel3D == 1 ) then

!     3D or vel3D = 1

      select case(comp)
      case(1)
        elemvec = gradu(:,1,1) ! dudx
      case(2)
        elemvec = gradu(:,1,2) ! dudy
      case(3)
        elemvec = gradu(:,1,3) ! dudz
      case(4)
        elemvec = gradu(:,2,1) ! dvdx
      case(5)
        elemvec = gradu(:,2,2) ! dvdy
      case(6)
        elemvec = gradu(:,2,3) ! dvdz
      case(7)
        elemvec = gradu(:,3,1) ! dwdx
      case(8)
        elemvec = gradu(:,3,2) ! dwdy
      case(9)
        elemvec = gradu(:,3,3) ! dwdz
      case(10)
        elemvec = gradu(:,1,1) + gradu(:,2,2) + gradu(:,3,3) ! divergence
      case(11) ! effective strain rate = sqrt(2II_D)
        do ip = 1, nodalp
!         2 * second invariant
          elemvec(ip) = &
               sum ( ( gradu(ip,:,:) + transpose(gradu(ip,:,:)) ) ** 2 ) / 2
        end do
        elemvec = sqrt ( elemvec )
      case default
        call errormsg_case_default ( 'stokes_deriv', &
          'comp', int_value=comp )
      end select

    end if

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( detF )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )

    end if

  end subroutine stokes_deriv


! D tensor

  subroutine stokes_D_tensor ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer, save :: maxcomp
    integer :: j


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

!     test size of elemvec

      if ( coorsys <= 1 .and. vel3D == 0 ) then
        maxcomp = 3 + coorsys
      else
        maxcomp = 6
      end if

      if ( size(elemvec) /= nodalp*maxcomp ) then
        write(*,'(/2(a/))') 'Error in stokes_D_tensor:', &
          ' element vector has incorrect size for a symmetric tensor'
        stop
      end if

!     allocate arrays

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ncompu*ndf), uvector(ndf,ncompu) )
      allocate ( gradu(nodalp,ncompu,ncompu) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( work2(nodalp,maxcomp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( vel3D == 1 ) then

      gradu(:,:,3) = 0

      if ( coorsys == 1 ) then
        where ( x(:,2) < 1e-10_dp )
          gradu(:,2,3) = -gradu(:,3,2)       ! -du_theta / dr (r=0)
          gradu(:,3,3) = gradu(:,2,2)        ! du_r / dr (r=0)
        elsewhere
          gradu(:,2,3) = -uvector(:,3) / x(:,2)   ! -u_theta / r
          gradu(:,3,3) = uvector(:,2) / x(:,2)   ! u_r / r
        end where
      end if

    end if

    if ( coorsys <= 1 .and. vel3D == 0 ) then

      work2(:,1) = gradu(:,1,1)                        ! D_xx
      work2(:,2) = ( gradu(:,1,2) + gradu(:,2,1) ) / 2 ! D_xy
      work2(:,3) = gradu(:,2,2)                        ! D_yy

      if ( coorsys == 1 ) then
!       axisymmetric
        where ( x(:,2) < 1e-10_dp )
          work2(:,4) = gradu(:,2,2)           ! du_r / dr (r=0)
        elsewhere
          work2(:,4) = uvector(:,2) / x(:,2)  ! u_r / r
        end where
      end if

    else if ( coorsys == 2 .or. vel3D == 1 ) then

!     3D or vel3D = 1

      work2(:,1) = gradu(:,1,1)                        ! D_xx
      work2(:,2) = ( gradu(:,1,2) + gradu(:,2,1) ) / 2 ! D_xy
      work2(:,3) = ( gradu(:,1,3) + gradu(:,3,1) ) / 2 ! D_xz
      work2(:,4) = gradu(:,2,2)                        ! D_yy
      work2(:,5) = ( gradu(:,2,3) + gradu(:,3,2) ) / 2 ! D_yz
      work2(:,6) = gradu(:,3,3)                        ! D_zz

    end if

    elemvec = reshape ( work2, [nodalp*maxcomp] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( detF )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( work2 )

    end if

  end subroutine stokes_D_tensor


! velocity gradient tensor

  subroutine stokes_gradu_tensor ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer, save :: maxcomp
    integer :: j


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

!     test size of elemvec

      if ( coorsys <= 1 .and. vel3D == 0 ) then
        maxcomp = 4 + coorsys
      else
        maxcomp = 9
      end if

      if ( size(elemvec) /= nodalp*maxcomp ) then
        write(*,'(/2(a/))') 'Error in stokes_gradu_tensor:', &
          ' element vector has incorrect size for an unsymmetric tensor'
        stop
      end if

!     allocate arrays

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ncompu*ndf), uvector(ndf,ncompu) )
      allocate ( gradu(nodalp,ncompu,ncompu) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( work2(nodalp,maxcomp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( vel3D == 1 ) then

      gradu(:,:,3) = 0

      if ( coorsys == 1 ) then
        where ( x(:,2) < 1e-10_dp )
          gradu(:,2,3) = -gradu(:,3,2)       ! -du_theta / dr (r=0)
          gradu(:,3,3) = gradu(:,2,2)        ! du_r / dr (r=0)
        elsewhere
          gradu(:,2,3) = -uvector(:,3) / x(:,2)   ! -u_theta / r
          gradu(:,3,3) = uvector(:,2) / x(:,2)   ! u_r / r
        end where
      end if

    end if

    if ( coorsys <= 1 .and. vel3D == 0 ) then

      work2(:,1) = gradu(:,1,1) ! dudx
      work2(:,2) = gradu(:,1,2) ! dudy
      work2(:,3) = gradu(:,2,1) ! dvdx
      work2(:,4) = gradu(:,2,2) ! dvdy

      if ( coorsys == 1 ) then
!       axisymmetric
        where ( x(:,2) < 1e-10_dp )
          work2(:,5) = gradu(:,2,2)           ! du_r / dr (r=0)
        elsewhere
          work2(:,5) = uvector(:,2) / x(:,2)  ! u_r / r
        end where
      end if

    else if ( coorsys == 2 .or. vel3D == 1 ) then

!    3D or vel3D = 1

      work2(:,1) = gradu(:,1,1) ! dudx
      work2(:,2) = gradu(:,1,2) ! dudy
      work2(:,3) = gradu(:,1,3) ! dudz
      work2(:,4) = gradu(:,2,1) ! dvdx
      work2(:,5) = gradu(:,2,2) ! dvdy
      work2(:,6) = gradu(:,2,3) ! dvdz
      work2(:,7) = gradu(:,3,1) ! dwdx
      work2(:,8) = gradu(:,3,2) ! dwdy
      work2(:,9) = gradu(:,3,3) ! dwdz

    end if

    elemvec = reshape ( work2, [nodalp*maxcomp] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( detF )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( work2 )

    end if

  end subroutine stokes_gradu_tensor


! compute persistence-of-straining tensor
!
! P_pq =
!  ( eigval(p) - eigval(q) ) * eigvec(p) dot W dot eigvec(q)
!                  .
!  + eigvec(p) dot D dot eigvec(q)
!
! where flow is (quasi-)steady thus material derivative of D
! is reduced to u dot nabla D.
!
! Note:
!   P_pq represents matrix components of P in principal coordinate
!   of D tensor.
!   P_pq matrix is symmetric with zero diagonal.
!   Get D tensor from oldvector for computation of u dot nabla D.

  subroutine stokes_P_tensor ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use stokes_globals_m

#if NO_LIBHSL
    use eig2D3D_m
#else
    use eig2D3D_m, only: eig2x2
    use hsl_ea23_m, only: eig3x3
#endif

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    integer :: j
    integer, save :: maxcomp, maxcompD, maxcompW

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

!     test size of elemvec

      if ( coorsys <= 1 .and. vel3D == 0 ) then
        maxcomp  = 1
        maxcompD = 3 + coorsys
        maxcompW = 1
      else
        maxcomp  = 3
        maxcompD = 6
        maxcompW = 3
      end if

      if ( size(elemvec) /= nodalp*maxcomp ) then
        write(*,'(/2(a/))') 'Error in stokes_P_tensor:', &
          ' element vector has incorrect size for P tensor'
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

!   get velocity vector

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

!   velocity at nodal points

    ugvector = matmul ( phi, uvector )

!   compute velocity gradient at nodal points

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( vel3D == 1 ) then

      gradu(:,:,3) = 0

      if ( coorsys == 1 ) then
        where ( x(:,2) < 1e-10_dp )
          gradu(:,2,3) = -gradu(:,3,2)       ! -du_theta / dr (r=0)
          gradu(:,3,3) = gradu(:,2,2)        ! du_r / dr (r=0)
        elsewhere
          gradu(:,2,3) = -uvector(:,3) / x(:,2)   ! -u_theta / r
          gradu(:,3,3) = uvector(:,2) / x(:,2)   ! u_r / r
        end where
      end if

    end if

!   compute rotation-rate tensor

    if ( coorsys <= 1 .and. vel3D == 0 ) then

      work2(:,1) = ( gradu(:,1,2) - gradu(:,2,1) ) / 2 ! W_xy

    else if ( coorsys == 2 .or. vel3D == 1 ) then

      work2(:,1) = ( gradu(:,1,2) - gradu(:,2,1) ) / 2 ! W_xy
      work2(:,2) = ( gradu(:,1,3) - gradu(:,3,1) ) / 2 ! W_xz
      work2(:,3) = ( gradu(:,2,3) - gradu(:,3,2) ) / 2 ! W_yz

    end if

!   get strain-rate tensor (D) components

    call get_vector ( mesh, problem, oldvectors%v(3)%p, elgrp, elem, &
      work, layer=layer )

    Dmat = reshape ( work, [ndf,maxcompD] )

!   D components at nodal points

    work10 = matmul ( phi, Dmat )

!   Eulerian frame or ALE

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity at tn

      call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, &
        elem, um )

      umesh = reshape ( um, [ndf,ndim] )

      u_mesh = matmul ( phi, umesh )

!     (u-ugrid) dot nabla operator

      do j = 1, nodalp
        tmp(j,:) = &
               matmul ( dphidx(j,:,:), ugvector(j,1:ndim) - u_mesh(j,:) )
      end do

    else

!     Eulerian frame

      u_mesh = 0

!     u dot nabla operator

      do j = 1, nodalp
        tmp(j,:) = matmul ( dphidx(j,:,:), ugvector(j,1:ndim) )
      end do

    end if

!   u dot nabla D at nodal points

    work4 = matmul ( tmp, Dmat )

!   eigenvalues and eigenvectors of D at nodal points

    do j = 1, nodalp

      if ( coorsys <= 1 .and. vel3D == 0 ) then

        call eig2x2 ( work10(j,1:3), eigval(j,:), eigvec(j,:,:) )

      else if ( coorsys == 2 .or. vel3D == 1 ) then

        call eig3x3 ( work10(j,:), eigval(j,:), eigvec(j,:,:) )

      end if

    end do

!   persistence-of-straining: work11(nodalp,maxcomp)

    if ( coorsys <= 1 .and. vel3D == 0 ) then

!     compute 12-component

      work11(:,1) = ( eigval(:,1) - eigval(:,2) ) * &
                      work2(:,1) * &
                      ( eigvec(:,1,1) * eigvec(:,2,2) - &
                        eigvec(:,2,1) * eigvec(:,1,2) ) &
                  + &
                    ( eigvec(:,1,1) * work4(:,1) + &
                      eigvec(:,2,1) * work4(:,2) ) * eigvec(:,1,2) &
                  + ( eigvec(:,1,1) * work4(:,2) + &
                      eigvec(:,2,1) * work4(:,3) ) * eigvec(:,2,2)

    else if ( coorsys == 2 .or. vel3D == 1 ) then

!     compute 12-component

      work11(:,1) =  ( eigval(:,1) - eigval(:,2) ) * &
                     ( work2(:,1) * &
                       ( eigvec(:,1,1) * eigvec(:,2,2) - &
                         eigvec(:,2,1) * eigvec(:,1,2) ) + &
                       work2(:,2) * &
                       ( eigvec(:,1,1) * eigvec(:,3,2) - &
                         eigvec(:,3,1) * eigvec(:,1,2) ) + &
                       work2(:,3) * &
                       ( eigvec(:,2,1) * eigvec(:,3,2) - &
                         eigvec(:,3,1) * eigvec(:,2,2) ) ) &
                   + ( eigvec(:,1,1) * work4(:,1) + &
                       eigvec(:,2,1) * work4(:,2) + &
                       eigvec(:,3,1) * work4(:,3) ) * eigvec(:,1,2) &
                   + ( eigvec(:,1,1) * work4(:,2) + &
                       eigvec(:,2,1) * work4(:,4) + &
                       eigvec(:,3,1) * work4(:,5) ) * eigvec(:,2,2) &
                   + ( eigvec(:,1,1) * work4(:,3) + &
                       eigvec(:,2,1) * work4(:,5) + &
                       eigvec(:,3,1) * work4(:,6) ) * eigvec(:,3,2)

!     compute 13-component

      work11(:,2) =  ( eigval(:,1) - eigval(:,3) ) * &
                     ( work2(:,1) * &
                       ( eigvec(:,1,1) * eigvec(:,2,3) - &
                         eigvec(:,2,1) * eigvec(:,1,3) ) + &
                       work2(:,2) * &
                       ( eigvec(:,1,1) * eigvec(:,3,3) - &
                         eigvec(:,3,1) * eigvec(:,1,3) ) + &
                       work2(:,3) * &
                       ( eigvec(:,2,1) * eigvec(:,3,3) - &
                         eigvec(:,3,1) * eigvec(:,2,3) ) ) &
                   + ( eigvec(:,1,1) * work4(:,1) + &
                       eigvec(:,2,1) * work4(:,2) + &
                       eigvec(:,3,1) * work4(:,3) ) * eigvec(:,1,3) &
                   + ( eigvec(:,1,1) * work4(:,2) + &
                       eigvec(:,2,1) * work4(:,4) + &
                       eigvec(:,3,1) * work4(:,5) ) * eigvec(:,2,3) &
                   + ( eigvec(:,1,1) * work4(:,3) + &
                       eigvec(:,2,1) * work4(:,5) + &
                       eigvec(:,3,1) * work4(:,6) ) * eigvec(:,3,3)

!     compute 23-component

      work11(:,3) =  ( eigval(:,2) - eigval(:,3) ) * &
                     ( work2(:,1) * &
                       ( eigvec(:,1,2) * eigvec(:,2,3) - &
                         eigvec(:,2,2) * eigvec(:,1,3) ) + &
                       work2(:,2) * &
                       ( eigvec(:,1,2) * eigvec(:,3,3) - &
                         eigvec(:,3,2) * eigvec(:,1,3) ) + &
                       work2(:,3) * &
                       ( eigvec(:,2,2) * eigvec(:,3,3) - &
                         eigvec(:,3,2) * eigvec(:,2,3) ) ) &
                   + ( eigvec(:,1,2) * work4(:,1) + &
                       eigvec(:,2,2) * work4(:,2) + &
                       eigvec(:,3,2) * work4(:,3) ) * eigvec(:,1,3) &
                   + ( eigvec(:,1,2) * work4(:,2) + &
                       eigvec(:,2,2) * work4(:,4) + &
                       eigvec(:,3,2) * work4(:,5) ) * eigvec(:,2,3) &
                   + ( eigvec(:,1,2) * work4(:,3) + &
                       eigvec(:,2,2) * work4(:,5) + &
                       eigvec(:,3,2) * work4(:,6) ) * eigvec(:,3,3)

    end if

    elemvec = reshape ( work11, [nodalp*maxcomp] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ncompu*ndf), uvector(ndf,ncompu), ugvector(nodalp,ncompu) )
      allocate ( gradu(nodalp,ncompu,ncompu), work(ndf*maxcompD) )
      allocate ( Dmat(ndf,maxcompD) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( work2(nodalp,maxcompW) )
      allocate ( tmp(nodalp,ndf), work4(nodalp,maxcompD) )
      allocate ( work10(nodalp,maxcompD))
      allocate ( eigval(nodalp,ncompu) )
      allocate ( eigvec(nodalp,ncompu,ncompu) )
      allocate ( work11(nodalp,maxcomp) )
      allocate ( um(ndf*ndim), u_mesh(nodalp,ndim), umesh(ndf,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( detF )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, ugvector, gradu, work )
      deallocate ( Dmat )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( work2 )
      deallocate ( tmp, work4, work10 )
      deallocate ( eigval, eigvec )
      deallocate ( work11 )
      deallocate ( um, u_mesh, umesh )

    end subroutine deallocate_arrays

  end subroutine stokes_P_tensor


! compute a flow classification parameter
! based on a persistence-of-straining concept
!     sqrt ( tr P^2 / 2 )
!    --------------------
!          tr D^2
! Note:
!   Get D tensor from oldvectors while P tensor is computed
!   in the routine in the same way as the routine stokes_P_tensor.

  subroutine stokes_P_scalar ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use stokes_globals_m

#if NO_LIBHSL
    use eig2D3D_m
#else
    use eig2D3D_m, only: eig2x2
    use hsl_ea23_m, only: eig3x3
#endif

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    integer :: j
    integer, save :: maxcomp, maxcompD, maxcompW
    real(dp) :: trD2

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

!     test size of elemvec

      if ( coorsys <= 1 .and. vel3D == 0 ) then
        maxcomp  = 1
        maxcompD = 3 + coorsys
        maxcompW = 1
      else
        maxcomp  = 3
        maxcompD = 6
        maxcompW = 3
      end if

      if ( size(elemvec) /= nodalp ) then
        write(*,'(/2(a/))') 'Error in stokes_P_scalar:', &
          ' element vector has incorrect size for a scalar'
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

!   get velocity vector

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

!   velocity at nodal points

    ugvector = matmul ( phi, uvector )

!   compute velocity gradient at nodal points

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( vel3D == 1 ) then

      gradu(:,:,3) = 0

      if ( coorsys == 1 ) then
        where ( x(:,2) < 1e-10_dp )
          gradu(:,2,3) = -gradu(:,3,2)       ! -du_theta / dr (r=0)
          gradu(:,3,3) = gradu(:,2,2)        ! du_r / dr (r=0)
        elsewhere
          gradu(:,2,3) = -uvector(:,3) / x(:,2)   ! -u_theta / r
          gradu(:,3,3) = uvector(:,2) / x(:,2)   ! u_r / r
        end where
      end if

    end if

!   compute rotation-rate tensor

    if ( coorsys <= 1 .and. vel3D == 0 ) then

      work2(:,1) = ( gradu(:,1,2) - gradu(:,2,1) ) / 2 ! W_xy

    else if ( coorsys == 2 .or. vel3D == 1 ) then

      work2(:,1) = ( gradu(:,1,2) - gradu(:,2,1) ) / 2 ! W_xy
      work2(:,2) = ( gradu(:,1,3) - gradu(:,3,1) ) / 2 ! W_xz
      work2(:,3) = ( gradu(:,2,3) - gradu(:,3,2) ) / 2 ! W_yz

    end if

!   get strain-rate tensor (D) components

    call get_vector ( mesh, problem, oldvectors%v(3)%p, elgrp, elem, &
      work, layer=layer )

    Dmat = reshape ( work, [ndf,maxcompD] )

!   D components at nodal points

    work10 = matmul ( phi, Dmat )

!   Eulerian frame or ALE

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity at tn

      call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, &
        elem, um )

      umesh = reshape ( um, [ndf,ndim] )

      u_mesh = matmul ( phi, umesh )

!     (u-ugrid) dot nabla operator

      do j = 1, nodalp
        tmp(j,:) = &
               matmul ( dphidx(j,:,:), ugvector(j,1:ndim) - u_mesh(j,:) )
      end do

    else

!     Eulerian frame

      u_mesh = 0

!     u dot nabla operator

      do j = 1, nodalp
        tmp(j,:) = matmul ( dphidx(j,:,:), ugvector(j,1:ndim) )
      end do

    end if

!   u dot nabla D at nodal points

    work4 = matmul ( tmp, Dmat )

!   eigenvalues and eigenvectors of D at nodal points

    do j = 1, nodalp

      if ( coorsys <= 1 .and. vel3D == 0 ) then

        call eig2x2 ( work10(j,1:3), eigval(j,:), eigvec(j,:,:) )

      else if ( coorsys == 2 .or. vel3D == 1 ) then

        call eig3x3 ( work10(j,:), eigval(j,:), eigvec(j,:,:) )

      end if

    end do

!   persistence-of-straining at nodal points

    if ( coorsys <= 1 .and. vel3D == 0 ) then

!     12-component

      work11(:,1) = ( eigval(:,1) - eigval(:,2) ) * &
                      work2(:,1) * &
                      ( eigvec(:,1,1) * eigvec(:,2,2) - &
                        eigvec(:,2,1) * eigvec(:,1,2) ) &
                  + &
                    ( eigvec(:,1,1) * work4(:,1) + &
                      eigvec(:,2,1) * work4(:,2) ) * eigvec(:,1,2) &
                  + ( eigvec(:,1,1) * work4(:,2) + &
                      eigvec(:,2,1) * work4(:,3) ) * eigvec(:,2,2)

    else if ( coorsys == 2 .or. vel3D == 1 ) then

!     compute 12-component

      work11(:,1) =  ( eigval(:,1) - eigval(:,2) ) * &
                     ( work2(:,1) * &
                       ( eigvec(:,1,1) * eigvec(:,2,2) - &
                         eigvec(:,2,1) * eigvec(:,1,2) ) + &
                       work2(:,2) * &
                       ( eigvec(:,1,1) * eigvec(:,3,2) - &
                         eigvec(:,3,1) * eigvec(:,1,2) ) + &
                       work2(:,3) * &
                       ( eigvec(:,2,1) * eigvec(:,3,2) - &
                         eigvec(:,3,1) * eigvec(:,2,2) ) ) &
                   + ( eigvec(:,1,1) * work4(:,1) + &
                       eigvec(:,2,1) * work4(:,2) + &
                       eigvec(:,3,1) * work4(:,3) ) * eigvec(:,1,2) &
                   + ( eigvec(:,1,1) * work4(:,2) + &
                       eigvec(:,2,1) * work4(:,4) + &
                       eigvec(:,3,1) * work4(:,5) ) * eigvec(:,2,2) &
                   + ( eigvec(:,1,1) * work4(:,3) + &
                       eigvec(:,2,1) * work4(:,5) + &
                       eigvec(:,3,1) * work4(:,6) ) * eigvec(:,3,2)

!     compute 13-component

      work11(:,2) =  ( eigval(:,1) - eigval(:,3) ) * &
                     ( work2(:,1) * &
                       ( eigvec(:,1,1) * eigvec(:,2,3) - &
                         eigvec(:,2,1) * eigvec(:,1,3) ) + &
                       work2(:,2) * &
                       ( eigvec(:,1,1) * eigvec(:,3,3) - &
                         eigvec(:,3,1) * eigvec(:,1,3) ) + &
                       work2(:,3) * &
                       ( eigvec(:,2,1) * eigvec(:,3,3) - &
                         eigvec(:,3,1) * eigvec(:,2,3) ) ) &
                   + ( eigvec(:,1,1) * work4(:,1) + &
                       eigvec(:,2,1) * work4(:,2) + &
                       eigvec(:,3,1) * work4(:,3) ) * eigvec(:,1,3) &
                   + ( eigvec(:,1,1) * work4(:,2) + &
                       eigvec(:,2,1) * work4(:,4) + &
                       eigvec(:,3,1) * work4(:,5) ) * eigvec(:,2,3) &
                   + ( eigvec(:,1,1) * work4(:,3) + &
                       eigvec(:,2,1) * work4(:,5) + &
                       eigvec(:,3,1) * work4(:,6) ) * eigvec(:,3,3)

!     compute 23-component

      work11(:,3) =  ( eigval(:,2) - eigval(:,3) ) * &
                     ( work2(:,1) * &
                       ( eigvec(:,1,2) * eigvec(:,2,3) - &
                         eigvec(:,2,2) * eigvec(:,1,3) ) + &
                       work2(:,2) * &
                       ( eigvec(:,1,2) * eigvec(:,3,3) - &
                         eigvec(:,3,2) * eigvec(:,1,3) ) + &
                       work2(:,3) * &
                       ( eigvec(:,2,2) * eigvec(:,3,3) - &
                         eigvec(:,3,2) * eigvec(:,2,3) ) ) &
                   + ( eigvec(:,1,2) * work4(:,1) + &
                       eigvec(:,2,2) * work4(:,2) + &
                       eigvec(:,3,2) * work4(:,3) ) * eigvec(:,1,3) &
                   + ( eigvec(:,1,2) * work4(:,2) + &
                       eigvec(:,2,2) * work4(:,4) + &
                       eigvec(:,3,2) * work4(:,5) ) * eigvec(:,2,3) &
                   + ( eigvec(:,1,2) * work4(:,3) + &
                       eigvec(:,2,2) * work4(:,5) + &
                       eigvec(:,3,2) * work4(:,6) ) * eigvec(:,3,3)

    end if

!   flow classification parameter at nodal points

    if ( coorsys <= 1 .and. vel3D == 0 ) then

      do j = 1, nodalp

        trD2 = work10(j,1)**2 + &
             2*work10(j,2)**2 + &
               work10(j,3)**2

        if ( coorsys == 1 ) then ! axisymmetric
          trD2 = trD2 + work10(j,4)**2
        end if

        if ( trD2 == 0._dp ) then
          work8(j) = 1.e+10
        else
          work8(j) = sqrt ( work11(j,1)**2 ) / trD2
        end if

      end do

    else if ( coorsys == 2 .or. vel3D == 1 ) then

      do j = 1, nodalp

        trD2 = work10(j,1)**2 + &
             2*work10(j,2)**2 + &
             2*work10(j,3)**2 + &
               work10(j,4)**2 + &
             2*work10(j,5)**2 + &
               work10(j,6)**2

        if ( trD2 == 0._dp ) then
          work8(j) = 1.e+10
        else
          work8(j) = sqrt ( work11(j,1)**2 &
                          + work11(j,2)**2 &
                          + work11(j,3)**2 ) &
                    / trD2
        end if

      end do

    end if

    elemvec = work8

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ncompu*ndf), uvector(ndf,ncompu), ugvector(nodalp,ncompu) )
      allocate ( gradu(nodalp,ncompu,ncompu), work(ndf*maxcompD) )
      allocate ( Dmat(ndf,maxcompD) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( work2(nodalp,maxcompW) )
      allocate ( tmp(nodalp,ndf), work4(nodalp,maxcompD) )
      allocate ( work10(nodalp,maxcompD))
      allocate ( eigval(nodalp,ncompu) )
      allocate ( eigvec(nodalp,ncompu,ncompu) )
      allocate ( work11(nodalp,maxcomp) )
      allocate ( work8(nodalp) )
      allocate ( um(ndf*ndim), u_mesh(nodalp,ndim), umesh(ndf,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( detF )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, ugvector, gradu, work )
      deallocate ( Dmat )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( work2 )
      deallocate ( tmp, work4, work10 )
      deallocate ( eigval, eigvec )
      deallocate ( work11 )
      deallocate ( work8 )
      deallocate ( um, u_mesh, umesh )

    end subroutine deallocate_arrays

  end subroutine stokes_P_scalar


! compute stresses in all nodes (from velocity gradients)

  subroutine stokes_stress ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: comp, j, maxstress
    real(dp) :: eta


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

!     check component = coefficients%i(13)

      if ( coorsys <= 1 .and. vel3D == 0 ) then
        maxstress = 3 + coorsys
      else
        maxstress = 6
      end if

      call check ( coefficients, 'stokes_stress', indexarray=[13], &
        minimum=[1], maximum=[maxstress] )

!     allocate arrays

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ncompu*ndf), uvector(ndf,ncompu) )
      allocate ( gradu(nodalp,ncompu,ncompu) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( vel3D == 1 ) then

      gradu(:,:,3) = 0

      if ( coorsys == 1 ) then
        where ( x(:,2) < 1e-10_dp )
          gradu(:,2,3) = -gradu(:,3,2)       ! -du_theta / dr (r=0)
          gradu(:,3,3) = gradu(:,2,2)        ! du_r / dr (r=0)
        elsewhere
          gradu(:,2,3) = -uvector(:,3) / x(:,2)  ! -u_theta / r
          gradu(:,3,3) = uvector(:,2) / x(:,2)   ! u_r / r
        end where
      end if

    end if

    eta = coefficients%r(1)

    comp = coefficients%i(13)

    if ( coorsys <= 1 .and. vel3D == 0 ) then

      select case(comp)
      case(1)
        elemvec = 2 * eta * gradu(:,1,1)                ! 2 eta dudx
      case(2)
        elemvec = eta * ( gradu(:,1,2) + gradu(:,2,1) ) ! eta ( dudy + dvdx )
      case(3)
        elemvec = 2 * eta * gradu(:,2,2)                ! 2 eta dvdy
      case(4)
        where ( x(:,2) < 1e-10_dp )
          elemvec = 2 * eta * gradu(:,2,2)              ! 2 eta du_r / dr (r=0)
        elsewhere
          elemvec = 2 * eta * uvector(:,2) / x(:,2)    ! 2 eta u_r / r
        end where
      case default
        call errormsg_case_default ( 'stokes_stress', &
          'comp', int_value=comp )
      end select

    else if ( coorsys == 2 .or. vel3D == 1 ) then

      select case(comp)
      case(1)
        elemvec = 2 * eta * gradu(:,1,1)                ! 2 eta dudx
      case(2)
        elemvec = eta * ( gradu(:,1,2) + gradu(:,2,1) ) ! eta ( dudy + dvdx )
      case(3)
        elemvec = eta * ( gradu(:,1,3) + gradu(:,3,1) ) ! eta ( dudz + dwdx )
      case(4)
        elemvec = 2 * eta * gradu(:,2,2)                ! 2 eta dvdy
      case(5)
        elemvec = eta * ( gradu(:,2,3) + gradu(:,3,2) ) ! eta ( dvdz + dwdy )
      case(6)
        elemvec = 2 * eta * gradu(:,3,3)                ! 2 eta dwdz
      case default
        call errormsg_case_default ( 'stokes_stress', &
          'comp', int_value=comp )
      end select

    end if

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( detF )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )

    end if

  end subroutine stokes_stress


! stokes stress tensor (for use with drag computations)

  subroutine stokes_stress_tensor ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer, save :: maxstress
    integer :: j
    real(dp) :: eta


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

!     test size of elemvec

      if ( coorsys <= 1 .and. vel3D == 0 ) then
        maxstress = 3 + coorsys
      else
        maxstress = 6
      end if

      if ( size(elemvec) /= nodalp*maxstress ) then
        write(*,'(/2(a/))') 'Error in stokes_stress_tensor:', &
          ' element vector has incorrect size for a symmetric tensor'
        stop
      end if

!     allocate arrays

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ncompu*ndf), uvector(ndf,ncompu) )
      allocate ( gradu(nodalp,ncompu,ncompu) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( work2(nodalp,maxstress) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( vel3D == 1 ) then

      gradu(:,:,3) = 0

      if ( coorsys == 1 ) then
        where ( x(:,2) < 1e-10_dp )
          gradu(:,2,3) = -gradu(:,3,2)       ! -du_theta / dr (r=0)
          gradu(:,3,3) = gradu(:,2,2)        ! du_r / dr (r=0)
        elsewhere
          gradu(:,2,3) = -uvector(:,3) / x(:,2)  ! -u_theta / r
          gradu(:,3,3) = uvector(:,2) / x(:,2)   ! u_r / r
        end where
      end if

    end if

    eta = coefficients%r(1)

    if ( coorsys <= 1 .and. vel3D == 0 ) then

      work2(:,1) = 2 * eta * gradu(:,1,1)                ! tau_xx
      work2(:,2) = eta * ( gradu(:,1,2) + gradu(:,2,1) ) ! tau_xy
      work2(:,3) = 2 * eta * gradu(:,2,2)                ! tau_yy

      if ( coorsys == 1 ) then
!       axisymmetric
        where ( x(:,2) < 1e-10_dp )
          work2(:,4) = 2 * eta * gradu(:,2,2)            ! 2 eta du_r / dr (r=0)
        elsewhere
          work2(:,4) = 2 * eta * uvector(:,2) / x(:,2)  ! 2 eta u_r / r
        end where
      end if

    else if ( coorsys == 2 .or. vel3D == 1 ) then

      work2(:,1) = 2 * eta * gradu(:,1,1)                ! tau_xx
      work2(:,2) = eta * ( gradu(:,1,2) + gradu(:,2,1) ) ! tau_xy
      work2(:,3) = eta * ( gradu(:,1,3) + gradu(:,3,1) ) ! tau_xz
      work2(:,4) = 2 * eta * gradu(:,2,2)                ! tau_yy
      work2(:,5) = eta * ( gradu(:,2,3) + gradu(:,3,2) ) ! tau_yz
      work2(:,6) = 2 * eta * gradu(:,3,3)                ! tau_zz

    end if

    elemvec = reshape ( work2, [nodalp*maxstress] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( detF )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( work2 )

    end if

  end subroutine stokes_stress_tensor



! compute pressures in all nodes

  subroutine stokes_pressure ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    if ( first ) then

!     first element in this group

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

      allocate ( psi(nodalp,ndfp), x(nodalp,ndim) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf) )
      allocate ( u(ndfp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      if ( coefficients%i(62) == 0 ) then
        call set_shape_function ( shapefuncp, xrnod, psi )
      end if

    end if

    if ( any ( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call get_coordinates ( mesh, elgrp, elem, x )
      call set_stokes_shape_function_global ( shapefuncp, x, x, psi, &
        coefficients%i(62)==1 )
    end if

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqpress], layer=layer )

    elemvec = matmul ( psi, u )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( psi, x )
      deallocate ( xrnod, phi )
      deallocate ( u )

    end if

  end subroutine stokes_pressure


! sample value of velocity in one node of the object

  subroutine stokes_sample_velocity ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, uvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: uvec


    integer :: elgrp, elem
    real(dp) :: xr(1,mesh%ndim)


    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

    call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

    allocate ( phi(1,ndf), u(ncompu*ndf), tmp(ndf,ncompu) )

!   set shape function in the point

    call set_shape_function ( shapefunc, xr, phi )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    tmp = reshape ( u, [ndf,ncompu] )

    uvec = matmul( phi(1,:), tmp )

    deallocate ( phi, u, tmp )

  end subroutine stokes_sample_velocity


! sample value of pressure in one node of the object

  subroutine stokes_sample_pressure ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, p )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: p


    integer :: elgrp, elem
    real(dp) :: xr(1,mesh%ndim)


    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

    call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

    allocate ( psi(1,ndfp), u(ndfp) )

!   set shape function in the point

    if ( any ( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      write(*,'(2a)') ' Error in stokes_sample_pressure:', &
                      ' global interpolation not implemented'
    end if

    call set_shape_function ( shapefuncp, xr, psi )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqpress], layer=layer )

    p = matmul( psi, u )

    deallocate ( psi, u )

  end subroutine stokes_sample_pressure


! sample value of velocity gradients in one node of the object

  subroutine stokes_sample_deriv ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, dvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: dvec


    integer :: elgrp, elem, comp, j, maxderiv, ip
    real(dp) :: xr(1,mesh%ndim)


    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

    call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

!   check component = coefficients%i(13)

    if ( coorsys <= 1 .and. vel3D == 0 ) then
      maxderiv = 8
    else
      maxderiv = 11
    end if

    call check ( coefficients, 'stokes_sample_deriv', indexarray=[13], &
      minimum=[1], maximum=[maxderiv] )

!   allocate arrays

    allocate ( detF(1) )
    allocate ( phi(1,ndf), x(nodalp,ndim) )
    allocate ( u(ncompu*ndf), uvector(ndf,ncompu), gradu(1,ncompu,ncompu) )
    allocate ( dphi(1,ndf,ndim), F(1,ndim,ndim) )
    allocate ( Finv(1,ndim,ndim), dphidx(1,ndf,ndim) )
    allocate ( xg(1,ndim), ugvector(1,ncompu) )

!   set shape function in the point

    call set_shape_function ( shapefunc, xr, phi, dphi )

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )
    ugvector = matmul ( phi, uvector )
    xg = matmul ( phi, x )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( vel3D == 1 ) then

      gradu(:,:,3) = 0

      if ( coorsys == 1 ) then
        where ( xg(:,2) < 1e-10_dp )
          gradu(:,2,3) = -gradu(:,3,2)       ! -du_theta / dr (r=0)
          gradu(:,3,3) = gradu(:,2,2)        ! du_r / dr (r=0)
        elsewhere
          gradu(:,2,3) = -ugvector(:,3) / xg(:,2)  ! -u_theta / r
          gradu(:,3,3) = ugvector(:,2) / xg(:,2)   ! u_r / r
        end where
      end if

    end if

    comp = coefficients%i(13)

    if ( coorsys <= 1 .and. vel3D == 0 ) then

!     2D

      select case(comp)
      case(1)
        dvec = gradu(:,1,1) ! dudx
      case(2)
        dvec = gradu(:,1,2) ! dudy
      case(3)
        dvec = gradu(:,2,1) ! dvdx
      case(4)
        dvec = gradu(:,2,2) ! dvdy
      case(5)
        dvec = gradu(:,2,1) - gradu(:,1,2) ! dvdx - dudy = vorticity
      case(6)
        if ( coorsys == 0 ) then
          write(*,'(/2(a/))') 'Error in stokes_deriv:', &
            ' gradu in theta direction not applicable for coorsys=0 '
          stop
        end if
        where ( xg(:,2) < 1e-10_dp )
          dvec = gradu(:,2,2)              ! du_r / dr (r=0)
        elsewhere
          dvec = ugvector(:,2) / xg(:,2)    ! u_r / r
        end where
      case(7)
        dvec = gradu(:,1,1) + gradu(:,2,2) ! dudx + dudy = divergence
        if ( coorsys == 1 ) then ! axisymmetric
          where ( xg(:,2) < 1e-10_dp )
            dvec = dvec + gradu(:,2,2)              ! du_r / dr (r=0)
          elsewhere
            dvec = dvec + ugvector(:,2) / xg(:,2)    ! u_r / r
          end where
        end if
      case(8) ! effective strain rate = sqrt(2II_D)
!       second invariant
        ip = 1
        dvec(ip) = &
             sum ( ( gradu(ip,:,:) + transpose(gradu(ip,:,:)) ) ** 2 ) / 2
        if ( coorsys == 1 ) then
!         axisymmetric
          if ( xg(ip,2) < 1e-10_dp ) then
!           use du_r / dr (r=0) for gradu_tt
            dvec(ip) = dvec(ip) + 2 * ( gradu(ip,2,2) ) ** 2
          else
!           use u_r / r for gradu_tt
            dvec(ip) = dvec(ip) + 2 * ( ugvector(ip,2) / xg(ip,2) ) ** 2
          end if
        end if
        dvec = sqrt ( dvec )
      case default
        call errormsg_case_default ( 'stokes_sample_deriv', &
          'comp', int_value=comp )
      end select

    else if ( coorsys == 2 .or. vel3D == 1 ) then

!     3D or vel3D = 1

      select case(comp)
      case(1)
        dvec = gradu(:,1,1) ! dudx
      case(2)
        dvec = gradu(:,1,2) ! dudy
      case(3)
        dvec = gradu(:,1,3) ! dudz
      case(4)
        dvec = gradu(:,2,1) ! dvdx
      case(5)
        dvec = gradu(:,2,2) ! dvdy
      case(6)
        dvec = gradu(:,2,3) ! dvdz
      case(7)
        dvec = gradu(:,3,1) ! dwdx
      case(8)
        dvec = gradu(:,3,2) ! dwdy
      case(9)
        dvec = gradu(:,3,3) ! dwdz
      case(10)
        dvec = gradu(:,1,1) + gradu(:,2,2) + gradu(:,3,3) ! divergence
      case(11) ! effective strain rate = sqrt(2II_D)
        ip = 1
!       2 * second invariant
        dvec(ip) = &
               sum ( ( gradu(ip,:,:) + transpose(gradu(ip,:,:)) ) ** 2 ) / 2
        dvec = sqrt ( dvec )
      case default
        call errormsg_case_default ( 'stokes_sample_deriv', &
          'comp', int_value=comp )
      end select

    end if

    deallocate ( detF )
    deallocate ( phi, x )
    deallocate ( u, uvector, gradu )
    deallocate ( dphi, F )
    deallocate ( Finv, dphidx )
    deallocate ( xg, ugvector )

  end subroutine stokes_sample_deriv


! surface tension on an object

  subroutine surface_tension_object ( mesh, problem, eleminfo, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(eleminfo_t), intent(in) :: eleminfo
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   variables

    integer :: elgrp, elem, object, elemo
    integer :: i, j, N, ip, nc
    real(dp) :: gammac


!   set globals fluid element

    elgrp  = eleminfo%elgrp

    call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

!   set globals object

    object = eleminfo%object

    call set_globals_stokes_object ( mesh, coefficients, object )

    elem   = eleminfo%elem
    elemo  = eleminfo%elemo
    ninti  = eleminfo%nintps  ! reset ninti to the number in nintps

!   allocate arrays

    allocate ( xig(ninti,ndim), x(nodalpo,ndim), xel(ndf,ndim) )
    allocate ( wg(ninti), curvel(ninti), normal(ninti,ndim) )
    allocate ( phi(ninti,ndf), dxdxi(ninti,ndim) )
    allocate ( detF(ninti), dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
    allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
    allocate ( theta(ninti,nodalpo), dtheta(ninti,nodalpo,ndim-1) )
    allocate ( xigo(ninti,ndim-1), dxdxis(ninti,ndim,2) )
    allocate ( xg(ninti,ndim), work2(ndf,ndim) )
    allocate ( work6(ninti,ndf,ndim), tauten(ninti,ndim,ndim) )

!   integration points and weighting factors

    xig  = mesh%objects(object)%refcoor_int(eleminfo%intps,:,elemo)
    xigo = mesh%objects(object)%xig(eleminfo%intps,:)
    wg   = mesh%objects(object)%wg(eleminfo%intps)

!   set shape functions

    call set_shape_function ( shapefunc, xig, phi, dphi ) ! velocity
    call set_shape_function ( shapefunco, xigo, theta, dtheta ) ! shape object

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, xel )

    call isoparametric_deformation ( xel, dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

!   compute deformed object element

    call get_coordinates_object ( mesh, elemo, x, object )

    xg =  matmul ( theta, x )

    if ( ndim == 2 ) then
      call isoparametric_deformation_curve ( x, dtheta(:,:,1), dxdxi, curvel, &
        normal )
      if ( coorsys == 1 ) then
        curvel = 2 * pi * xg(:,2) * curvel
      end if
    else if ( ndim == 3 ) then
      call isoparametric_deformation_surface ( x, dtheta, dxdxis, curvel, &
        normal )
    end if

!   surface tension coefficient

    gammac = coefficients%r(19)


    if ( matrix ) then
      write(*,'(/a/a/)') 'Error in surface_tension_object:', &
        ' No matrix to build. Call build_system with buildmatrix=.false.'
      stop
      elemmat = 0._dp
    end if


    if ( vector ) then

!     - (nabla v)^T:tau, tau=Gamma*(I-nn)

!     diagonal components of tau

      do i = 1, ndim
        tauten(:,i,i) = 1 - normal(:,i)**2
      end do

!     off-diagonal components of tau

      do i = 1, ndim-1
        do j = i+1, ndim
          tauten(:,i,j) = - normal(:,i)*normal(:,j)
          tauten(:,j,i) = tauten(:,i,j)
        end do
      end do

      do ip = 1, ninti
        work6(ip,:,:) = matmul ( dphidx(ip,:,:), tauten(ip,:,:) )
      end do
      if ( coorsys == 1 ) then
        do N = 1, ndf
          work6(:,N,2) = work6(:,N,2) + phi(:,N) / xg(:,2)
        end do
      end if
      do j = 1, ndim
        work2(:,j) = - gammac * matmul ( curvel * wg, work6(:,:,j) )
      end do

      nc = ndf*ndim

      if ( vel3D == 0 ) then
        elemvec = reshape ( work2, [nc] )
      else
        elemvec(1:nc) = reshape ( work2, [nc] )
        elemvec(nc+1:) = 0
      end if

    end if

!   deallocate arrays

    deallocate ( xig, x, xel )
    deallocate ( wg, curvel, normal )
    deallocate ( phi, dxdxi )
    deallocate ( detF, dphi, F )
    deallocate ( Finv, dphidx )
    deallocate ( theta, dtheta )
    deallocate ( xigo, dxdxis )
    deallocate ( xg, work2 )
    deallocate ( work6, tauten )

  end subroutine surface_tension_object


! removed routine because coefficient etaslip has moved

  subroutine stokes_natboun_slip ( mesh, problem, geometry, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: geometry, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    write(*,'(4(/a)/)') 'Error in stokes_natboun_slip:', &
      '  Element stokes_natboun_slip has been removed. Please, use the', &
      '  routine stokes_natboun_slip2 instead. NOTE: the etaslip coefficient', &
      '  has been moved from position 19 to position 3 of coefficients%r.'
    stop

  end subroutine stokes_natboun_slip


! Navier slip boundary conditions for geometries (curves,surfaces).
! For use with add_boundary_elements.

  subroutine stokes_natboun_slip2 ( mesh, problem, geometry, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: geometry, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, ip, N, M
    real(dp) :: etaslip, S_NMij, S_NMii

    if ( first ) then

!     first element on this geometry

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, &
        ndimr=mesh%ndim-1, geometry=geometry )

!     allocate arrays

      allocate ( wg(ninti), fg(ninti,ndim), curvel(ninti) )
      allocate ( tmp(ndf,ndim) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,ndim-1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim-1), dxdxis(ninti,ndim,ndim-1) )
      allocate ( work6(ninti,ndim,ndim) )
      allocate ( work2(ninti,ndim) )
      allocate ( pos(ndf,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

!     position arrays
      pos  = reshape ( [ ( i, i = 1, ndf*ndim ) ], [ ndf, ndim ] )

    end if

    call get_coordinates_geometry ( mesh, elem, x, ndimr=ndim-1, &
      geometry=geometry )

    call isoparametric_deformation_curved ( x, dphi, dxdxis, curvel, &
      normal )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      curvel = 2 * pi * xg(:,2) * curvel
    end if

    etaslip = coefficients%r(3)

!   diagonal blocks of I-nn

    do i = 1, ndim
      work6(:,i,i) = 1 - normal(:,i)**2
    end do

!   off-diagonal blocks of I-nn

    do i = 1, ndim-1
      do j = i+1, ndim
        work6(:,i,j) = - normal(:,i)*normal(:,j)
        work6(:,j,i) = work6(:,i,j)
      end do
    end do

    if ( vector ) then

      if ( coefficients%i(55) >= 0 ) then

!       determine prescribed wall velocity

        call evaluate_vector_coefficient_simple ( choice=coefficients%i(55), &
          value=coefficients%r(20:19+ndim), vfunc=coefficients%vfunc1(2)%p, &
          vfuncnr=coefficients%i(56), x=x, coef=fg )

!       (I - nn).uwall
        do ip = 1, ninti
          work2(ip,:) =  matmul ( work6(ip,:,:), fg(ip,:) )
        end do

        do j = 1, ndim
          do N = 1, ndf
            tmp(N,j) = etaslip * sum ( work2(:,j) * phi(:,N) * curvel * wg )
          end do
        end do

        elemvec = reshape ( tmp, [ndim*ndf] )

      else

        elemvec = 0

      end if

    end if

    if ( matrix ) then

!     diagonal blocks

      do N = 1, ndf
        do M = N, ndf
          do i = 1, ndim
            S_NMii = etaslip * &
                    sum ( work6(:,i,i) * phi(:,N) * phi(:,M) * curvel * wg )
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
              S_NMij = etaslip * &
                    sum ( work6(:,i,j) * phi(:,N) * phi(:,M) * curvel * wg )
              elemmat( pos(N,i), pos(M,j) ) = S_NMij
              elemmat( pos(M,j), pos(N,i) ) = S_NMij  ! = S_MNji symmetry
            end do
          end do
        end do
      end do

    end if

    if ( last ) then

!     last element on this geometry

      deallocate ( wg, fg, curvel, normal )
      deallocate ( tmp )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi, dxdxis )
      deallocate ( work6 )
      deallocate ( work2 )
      deallocate ( pos )

    end if

  end subroutine stokes_natboun_slip2


! Normal traction boundary condition for geometries (curves,surfaces).
! For use with add_boundary_elements.

  subroutine stokes_natboun_normal ( mesh, problem, geometry, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: geometry, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, N, nc

    if ( first ) then

!     first element on this geometry

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, &
        ndimr=mesh%ndim-1, geometry=geometry, maxvel3D=1 )

!     allocate arrays

      allocate ( wg(ninti), trac(ninti), curvel(ninti) )
      allocate ( tmp(ndf,ndim) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,ndim-1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim-1), dxdxis(ninti,ndim,ndim-1) )
      allocate ( pos(ndf,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

!     position arrays
      pos  = reshape ( [ ( i, i = 1, ndf*ndim ) ], [ ndf, ndim ] )

    end if

    call get_coordinates_geometry ( mesh, elem, x, ndimr=ndim-1, &
      geometry=geometry )

    call isoparametric_deformation_curved ( x, dphi, dxdxis, curvel, &
      normal )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      curvel = 2 * pi * xg(:,2) * curvel
    end if

    if ( vector ) then

      if ( coefficients%i(65) >= 0 ) then

!       determine prescribed traction

        call evaluate_scalar_coefficient_simple ( choice=coefficients%i(65), &
          value=coefficients%r(24), func=coefficients%func, &
          funcnr=coefficients%i(66), x=x, coef=trac )

        do j = 1, ndim
          do N = 1, ndf
            tmp(N,j) = sum ( trac * normal(:,j) * phi(:,N) * curvel * wg )
          end do
        end do

        elemvec = reshape ( tmp, [ndim*ndf] )

        nc = ndf*ndim

        if ( vel3D == 0 ) then
          elemvec = reshape ( tmp, [nc] )
        else
          elemvec(1:nc) = reshape ( tmp, [nc] )
          elemvec(nc+1:) = 0
        end if

      else

        elemvec = 0

      end if

    end if


    if ( matrix ) then
      write(*,'(/a/a/)') 'Error in stokes_natboun_normal:', &
        ' No matrix to build. Call build_system with buildmatrix=.false.'
      stop
      elemmat = 0._dp
    end if


    if ( last ) then

!     last element on this geometry

      deallocate ( wg, trac, curvel, normal )
      deallocate ( tmp )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi, dxdxis )
      deallocate ( pos )

    end if

  end subroutine stokes_natboun_normal


! Intersect and find reference and real coordinates on temporary ALE meshes
! of the integration points in xg. One step for use with first-order and
! (optional) two previous steps for use with second-order time integration.
! NOTE: data transfer is through global variables in stokes_globals_m that
! need to be properly allocated.

  subroutine set_globals_stokes_tALE ( mesh_ALE_np1, mesh_n, mesh_ALE_n, &
    mesh_nm1, deform_n, groups )

    use stokes_globals_m, only: nodalp, ndim, ninti, xg, xg_n, xg_nm1, &
      grpelm_n, xig_n, phi_n, grpelm_nm1, xig_nm1, phi_nm1, shapefunc, &
      dphi_n, detF_n, F_n, Finv_n

    use set_optional_m

!   meshes used:
!     mesh_ALE_np1: temporary ALE mesh for tn->tn+1
!     mesh_n:       computational mesh at tn
    type(mesh_t), intent(in) :: mesh_ALE_np1, mesh_n

!   optional meshes used (second-order time stepping):
!     mesh_ALE_n:   temporary ALE mesh for tn-1->tn
!     mesh_nm1:     computational mesh at tn-1
    type(mesh_t), intent(in), optional :: mesh_ALE_n, mesh_nm1

!   If present and .true. compute the deformation from t_n -> tn+1
!   default = .false.
    logical, intent(in), optional :: deform_n

!   If present: include element groups for intersection:
!   group(elgrp) = .true.  : include group
!   group(elgrp) = .false. : exclude group
!   Default: include all groups
    logical, dimension(:), intent(in), optional :: groups

    logical :: ldeform_n
    integer :: ip
    real(dp) :: xe(nodalp,ndim)


    ldeform_n = set_optional ( variable=deform_n, default=.false. )


!   first step

!   compute reference coordinates in mesh_ALE_np1

    call find_refcoor_points ( mesh_ALE_np1, groups, coor=xg, &
      grpelm=grpelm_n, refcoor=xig_n )

    if ( any ( grpelm_n == 0 ) ) then
      write(*,'(/a/a/)') 'Error set_globals_stokes_tALE:', &
        ' No reference coordinates found in mesh_ALE_np1'
      stop
    end if

    if ( mesh_ALE_np1%nelgrp > 1 .and. &
         any ( mesh_ALE_np1%element(grpelm_n(:,1))%elshape /= &
                                 mesh_ALE_np1%element(1)%elshape ) ) then
      write(*,'(/a/a/)') 'Error set_globals_stokes_tALE:', &
        ' Different element shapes not implemented. '
      stop
    end if


!   set shapefunction for step tn->tn+1

    if ( ldeform_n ) then
      call set_shape_function ( shapefunc, xig_n, phi_n, dphi_n )
    else
      call set_shape_function ( shapefunc, xig_n, phi_n )
    end if

!   compute coordinates
!   NOTE: element can be different for each integration point

    do ip = 1, ninti
      call get_coordinates ( mesh_n, grpelm_n(ip,1), grpelm_n(ip,2), xe )
      xg_n(ip,:) = matmul ( phi_n(ip,:), xe )
      if ( ldeform_n ) then
        call isoparametric_deformation ( xe, dphi_n(ip:ip,:,:), F_n(ip:ip,:,:),&
          Finv_n(ip:ip,:,:), detF_n(ip:ip) )
      end if
    end do


!   second step

    if ( present(mesh_ALE_n) ) then

!     compute reference coordinates in mesh_ALE_n

      call find_refcoor_points ( mesh_ALE_n, groups, coor=xg_n, &
        grpelm=grpelm_nm1, refcoor=xig_nm1 )

      if ( any ( grpelm_nm1 == 0 ) ) then
        write(*,'(/a/a/)') 'Error set_globals_stokes_tALE:', &
          ' No reference coordinates found in mesh_ALE_n'
        stop
      end if

      if ( mesh_ALE_n%nelgrp > 1 .and. &
           any ( mesh_ALE_n%element(grpelm_nm1(:,1))%elshape /= &
                                   mesh_ALE_n%element(1)%elshape ) ) then
        write(*,'(/a/a/)') 'Error set_globals_stokes_tALE:', &
          ' Different element shapes not implemented. '
        stop
      end if


!     set shapefunction for step tn-1->tn

      call set_shape_function ( shapefunc, xig_nm1, phi_nm1 )

!     compute coordinates
!     NOTE: element can be different for each integration point

      do ip = 1, ninti
        call get_coordinates ( mesh_nm1, grpelm_nm1(ip,1), grpelm_nm1(ip,2), &
          xe )
        xg_nm1(ip,:) = matmul ( phi_nm1(ip,:), xe )
      end do

    end if

  end subroutine set_globals_stokes_tALE


! set global parameters (internal element)

  subroutine set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D )

    use stokes_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp
!   the maximum value for vel3D allowed. default=0, i.e. no vel3D available.
    integer, intent(in), optional :: maxvel3D

    integer :: lmaxvel3D, ndimr

    lmaxvel3D = set_optional ( variable=maxvel3D, default=0 )

!   check size of coefficients

    call check ( coefficients, 'set_globals_stokes_vp', ncoefi=100, &
      ncoefr=50, indexarray=[23,34,37,62,67], minimum=[0,0,0,0,0], &
      maximum=[1,1,2,2,lmaxvel3D], wncoefi=150, wncoefr=100 )

    ndim = mesh%element(elgrp)%ndim
    nodalp = mesh%element(elgrp)%numnod
    physqvel = coefficients%i(6)
    physqpress = coefficients%i(7)
    layer = coefficients%i(38)
    if ( ndim == 3 ) then
      coorsys = 2
      vel3D = 0
      ncompu = ndim
    else
      coorsys = coefficients%i(23)
      vel3D = coefficients%i(67)
      ncompu = ndim + vel3D
    end if
    nsides = mesh%element(elgrp)%numsides
    nodalpb = mesh%element(elgrp)%sidnumnod
    globalshape = mesh%element(elgrp)%globalshape

!   set number of degrees of freedom velocity

    intpol = coefficients%i(1)

    ndimr = mesh%element(elgrp)%ndimr
    if ( intpol == 13 .and. any(mesh%element(elgrp)%p(:ndimr,2) /= 1 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_stokes_vp:', &
        ' For spectral elements GLL nodal distribution is required. '
      stop
    else if ( intpol == 20 .and. &
                            any(mesh%element(elgrp)%p(:ndimr,2) /= 0 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_stokes_vp:', &
        ' For high-order elements equidistant nodal distribution is required. '
      stop
    end if

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(51)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, 'set_globals_stokes_vp', ndf=ndf, ndfb=ndfb )

    if ( coefficients%i(34) == 0 .and. ndf /= nodalp ) then
      write(*,'(/5(a/))') 'Error in set_globals_stokes_vp:', &
        ' Number of degrees of freedom of the velocity shape function (ndf) ', &
        ' is different from the number of nodal points (nodalp) ', &
        ' If you really want ndf /= nodalp, set coefficients%i(34) == 1,', &
        ' but you are on your own!! '
      stop
    end if

!   set number of degrees of freedom pressure

    intpolp = coefficients%i(2)

    shapefuncp%globalshape = globalshape
    shapefuncp%interpolation = intpolp
    shapefuncp%numbering = 'regular'
    shapefuncp%p = coefficients%i(52)
    shapefuncp%spec_eval = 'gauss'
!   needed for spectral quads/hexahedra only to call Pp_at_GLL routine:
    shapefuncp%intrule = shapefunc%p + 1

    call set_ndf ( shapefuncp, 'set_globals_stokes_vp', ndf=ndfp )

!   set integration

    inttype = coefficients%i(40)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_stokes_vp:', &
        ' For spectral elements Gauss-Legendre-Lobatto integration ', &
        ' needs to be specified. '
      stop
    end if

    if ( coefficients%i(77) == 1 .or. SET_GAUSS_BY_ORDER ) then
      intrule = set_intrule ( globalshape, inttype, order=coefficients%i(10) )
      intrule2 = set_intrule2 ( globalshape, inttype, order=coefficients%i(75) )
    else
      intrule = coefficients%i(10)
      intrule2 = coefficients%i(75)
    end if
    nsubint = get_coefficient ( coefficients, index=32, default=1 )

    if ( globalshape == 'prism' .and. intrule2 == 0 ) then
      write(*,'(/a/3a/)') 'Error in set_globals_stokes_vp:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape
      stop
    end if

    if ( globalshape == 'pyramid' .and. inttype == 2 .and. intrule2 == 0 ) then
      write(*,'(/a/a/3a,i0/)') 'Error in set_globals_stokes_vp:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape, ' and inttype = ', inttype
      stop
    end if

    gauss%globalshape = globalshape
    gauss%intrule = intrule
    gauss%intrule2 = intrule2
    gauss%nsubint = nsubint
    gauss%inttype = inttype

    call set_ninti ( gauss, ninti )

    if ( intpol == 13 .and. ninti /= ndf ) then
      write(*,'(/5(a/))') 'Error in set_globals_stokes_vp:', &
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_stokes_vp


! set global parameters (boundary element)

  subroutine set_globals_stokes_vp_boun ( mesh, coefficients, curve, surface, &
    volume, ndimr, geometry, maxvel3D )

    use stokes_globals_m
    use set_optional_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in), optional :: curve, surface, volume
    integer, intent(in), optional :: ndimr, geometry
!   the maximum value for vel3D allowed. default=0, i.e. no vel3D available.
    integer, intent(in), optional :: maxvel3D

    integer :: lcurve, lsurface, lvolume
    integer :: lmaxvel3D

    lmaxvel3D = set_optional ( variable=maxvel3D, default=0 )

!   check size of coefficients

    call check ( coefficients, 'set_globals_stokes_vp_boun', ncoefi=100, &
      ncoefr=50, indexarray=[23,67], minimum=[0,0], &
      maximum=[1,lmaxvel3D], wncoefi=150, wncoefr=100 )

    physqvel = coefficients%i(6)
    physqpress = coefficients%i(7)
    layer = coefficients%i(38)

!   traditional interface (legacy)

    lcurve = set_optional ( variable=curve, default=0 )
    lsurface = set_optional ( variable=surface, default=0 )
    lvolume = set_optional ( variable=volume, default=0 )

!   interface with dimension of reference space of geometry

    if ( present(ndimr) .and. present(geometry) ) then
      select case (ndimr)
        case(1); lcurve = geometry
        case(2); lsurface = geometry
        case(3); lvolume = geometry
      case default
        call errormsg_case_default ( 'set_globals_stokes_vp_boun', &
          'ndimr', int_value=ndimr )
      end select
    end if

    if ( lcurve > 0 ) then
      ndim = mesh%curves(lcurve)%ndim
      if ( ndim == 3 ) then
        coorsys = 2
        vel3D = 0
        ncompu = ndim
      else
        coorsys = coefficients%i(23)
        vel3D = coefficients%i(67)
        ncompu = ndim + vel3D
      end if
      globalshape = mesh%curves(lcurve)%element%globalshape
      nodalp = mesh%curves(lcurve)%element%numnod
    else if ( lsurface > 0 ) then
      ndim = mesh%surfaces(lsurface)%ndim
      if ( ndim == 3 ) then
        coorsys = 2
        vel3D = 0
        ncompu = ndim
      else
        coorsys = coefficients%i(23)
        vel3D = coefficients%i(67)
        ncompu = ndim + vel3D
      end if
      globalshape = mesh%surfaces(lsurface)%element%globalshape
      nodalp = mesh%surfaces(lsurface)%element%numnod
    else if ( lvolume > 0 ) then
      ndim = mesh%volumes(lvolume)%ndim
      coorsys = 2 ! always 3D
      globalshape = mesh%volumes(lvolume)%element%globalshape
      nodalp = mesh%volumes(lvolume)%element%numnod
    end if

!   set number of degrees of freedom velocity

    intpol = coefficients%i(1)

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(51)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, 'set_globals_stokes_vp_boun', ndf=ndf )

    if ( ndf /= nodalp ) then
      write(*,'(/5(a/))') 'Error in set_globals_stokes_vp_boun:', &
        ' Number of degrees of freedom of the shape function (ndf) ', &
        ' is different from the number of nodal points in the element', &
        ' (nodalp) ', &
        ' This possibility (ndf /= nodalp) is not available.'
      stop
    end if

!   set number of degrees of freedom pressure

    intpolp = coefficients%i(2)

    shapefuncp%globalshape = globalshape
    shapefuncp%interpolation = intpolp
    shapefuncp%numbering = 'regular'
    shapefuncp%p = coefficients%i(52)
    shapefuncp%spec_eval = 'gauss'
!   needed for spectral quads only to call Pp_at_GLL routine:
    shapefuncp%intrule = shapefunc%p + 1

    call set_ndf ( shapefuncp, 'set_globals_stokes_vp_boun', ndf=ndfp )

!   set integration

    inttype = coefficients%i(40)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_stokes_vp_boun:', &
        ' For spectral elements Gauss-Legendre-Lobatto integration ', &
        ' needs to be specified. '
      stop
    end if

    if ( coefficients%i(77) == 1 .or. SET_GAUSS_BY_ORDER ) then
      intrule = set_intrule ( globalshape, inttype, order=coefficients%i(11) )
      intrule2 = set_intrule2 ( globalshape, inttype, order=coefficients%i(76) )
    else
      intrule = coefficients%i(11)
      intrule2 = coefficients%i(76)
    end if
    nsubint = get_coefficient ( coefficients, index=33, default=1 )

    if ( globalshape == 'prism' .and. intrule2 == 0 ) then
      write(*,'(/a/3a/)') 'Error in set_globals_stokes_vp_boun:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape
      stop
    end if

    if ( globalshape == 'pyramid' .and. inttype == 2 .and. intrule2 == 0 ) then
      write(*,'(/a/a/3a,i0/)') 'Error in set_globals_stokes_vp_boun:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape, ' and inttype = ', inttype
      stop
    end if

    gauss%globalshape = globalshape
    gauss%intrule = intrule
    gauss%intrule2 = intrule2
    gauss%nsubint = nsubint
    gauss%inttype = inttype

    call set_ninti ( gauss, ninti )

    if ( intpol == 13 .and. ninti /= ndf ) then
      write(*,'(/5(a/))') 'Error in set_globals_stokes_vp_boun:', &
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_stokes_vp_boun


! set global parameters (boundary Langrangian multiplier)

  subroutine set_globals_stokes_l_boun ( coefficients )

    use stokes_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   Lagrange multiplier interpolation

!   set number of degrees of freedom

    call set_ndf ( shapefunc, name_of_routine='set_globals_stokes_l_boun', &
      ndfl=ndfl, intpoll=intpoll )

    shapefuncl = shapefunc
    shapefuncl%interpolation = intpoll

  end subroutine set_globals_stokes_l_boun


! set global parameters (object)

  subroutine set_globals_stokes_object ( mesh, coefficients, object )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: object

    integer :: intpolo
    character (len=13) :: globalshapeo

!   check object

    if ( .not. mesh%objects(object)%topol ) then
      write(*,'(/a/a,i0/)') 'Error in set_globals_stokes_vp_object:', &
        ' Object does not have a topology. Object = ', object
      stop
    end if

    nodalpo = mesh%objects(object)%element%numnod

!   set shapefunction for object element (for geometrical shape).

    globalshapeo = mesh%objects(object)%element%globalshape
    intpolo = coefficients%i(35)

    shapefunco%globalshape = globalshapeo
    shapefunco%interpolation = intpolo
    shapefunco%numbering = 'standard'

  end subroutine set_globals_stokes_object


! Preamble for the stokes element

  subroutine set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, maxvel3D )

    use stokes_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    integer, intent(in), optional :: maxvel3D


    integer :: i, i1(3)


    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D )

!     allocate arrays

      allocate ( pos(ndf,ncompu), posp(ndfp) )
      allocate ( x(nodalp,ndim) )

      allocate ( xig(ninti,ndim), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim), psi(ninti,ndfp) )

      allocate ( detF(ninti), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( fg(ninti,ncompu), xg(ninti,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      if ( coefficients%i(62) == 0 ) then
        call set_shape_function ( shapefuncp, xig, psi )
      end if

!     position arrays
      pos  = reshape ( [ ( i, i = 1, ndf*ncompu ) ], [ ndf, ncompu ] )
      posp = [ ( i, i = 1, ndfp ) ] + ndf*ncompu

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      if ( first ) then

!       first element in this group

!       set globals

        call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D )

!       allocate arrays

        allocate ( pos(ndf,ncompu), posp(ndfp) )
        allocate ( x(nodalp,ndim) )

!       position arrays
        pos  = reshape ( [ ( i, i = 1, ndf*ncompu ) ], [ ndf, ncompu ] )
        posp = [ ( i, i = 1, ndfp ) ] + ndf*ncompu

      end if

!     now choose between elvector and user subroutines

      if ( coefficients%i(37) == 1 ) then

!       get element values of intrule, intrule2 and nsubint

        if ( gauss%globalshape == 'prism' .or. &
             gauss%globalshape == 'pyramid' .and. gauss%inttype == 2 ) then

!         intrule2 required
          call get_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, i1=i1 )

          if ( coefficients%i(77) == 1 .or. SET_GAUSS_BY_ORDER ) then
            intrule2 = set_intrule2 ( gauss%globalshape, gauss%inttype, &
                                      order=i1(3) )
          else
            intrule2 = i1(3)
          end if

          gauss%intrule2 = intrule2

        else

          call get_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, i1=i1(1:2) )

        end if

        if ( coefficients%i(77) == 1 .or. SET_GAUSS_BY_ORDER ) then
          intrule = set_intrule ( gauss%globalshape, gauss%inttype, &
                                  order=i1(1) )
        else
          intrule = i1(1)
        end if
        nsubint = i1(2)

        gauss%intrule = intrule
        gauss%nsubint = nsubint

        call set_ninti ( gauss, ninti )

!       allocate arrays

        allocate ( xig(ninti,ndim), wg(ninti) )

!       set Gauss integration

        call set_Gauss_integration ( gauss, xig, wg )

      else if ( coefficients%i(37) == 2 ) then

!       set ninti

        call coefficients%set_ninti_user ( mesh, problem, elgrp, elem, first, &
          last, coefficients, oldvectors )

!       allocate arrays

        allocate ( xig(ninti,ndim), wg(ninti) )

!       set Gauss integration

        call coefficients%set_Gauss_integration_user ( mesh, problem, &
          elgrp, elem, first, last, coefficients, oldvectors )

      end if

!     allocate arrays

      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim), psi(ninti,ndfp) )

      allocate ( detF(ninti), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( fg(ninti,ncompu), xg(ninti,ndim) )

!     set shape function

      call set_shape_function ( shapefunc, xig, phi, dphi )
      if ( coefficients%i(62) == 0 ) then
        call set_shape_function ( shapefuncp, xig, psi )
      end if

    end if

  end subroutine set_stokes_elem


! Unset the preamble for the stokes element
! (deallocate arrays allocated in set_... )

  subroutine unset_stokes_elem ( last, coefficients )

    use stokes_globals_m

    logical, intent(in) :: last
    type(coefficients_t), intent(in) :: coefficients

!   deallocate memory

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      deallocate ( pos, posp )
      deallocate ( x )
      deallocate ( xig, wg )
      deallocate ( phi, dphi, psi )

      deallocate ( detF, F )
      deallocate ( Finv, dphidx )
      deallocate ( fg, xg )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      if ( last ) then

!       last element in this group

        deallocate ( pos, posp )
        deallocate ( x )

      end if

      deallocate ( xig, wg )
      deallocate ( phi, dphi, psi )

      deallocate ( detF, F )
      deallocate ( Finv, dphidx )
      deallocate ( fg, xg )

    end if

  end subroutine unset_stokes_elem


! set shapefunction based on scaled global coordinates

  subroutine set_stokes_shape_function_global ( shapefunc, x, xg, phi, scaling )

    type(shapefunc_t), intent(in) :: shapefunc

!   coordinates of the nodes of the element
!   x(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: x

!   Global coordinates where phi must be computed:
!   xg(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xg

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   perform scaling of the shape function
    logical, intent(in) :: scaling

    integer :: ndim, ninti, i
    real(dp), allocatable, dimension(:) :: xmin, xmax
    real(dp), allocatable, dimension(:,:) :: xg_scale

    if ( scaling ) then

      ndim = size(x,2)
      ninti = size(xg,1)

      allocate ( xmin(ndim), xmax(ndim), xg_scale(ninti,ndim) )

      xmin = maxval ( x, dim=1 )  ! maximum in all coordinate directions
      xmax = minval ( x, dim=1 )  ! minimum in all coordinate directions

!     scale global coordinates:
!       - coordinates relative to the midpoint of the element (xmax+xmin)/2
!       - scale with half the extend xmax - xmin
      do i = 1, ndim
        xg_scale(:,i) = &
                  ( 2 * xg(:,i) - xmax(i) - xmin(i) ) / ( xmax(i) - xmin(i) )
      end do

      call set_shape_function_global ( shapefunc, xg_scale, phi )

      deallocate ( xmin, xmax, xg_scale )

    else

      call set_shape_function_global ( shapefunc, xg, phi )

    end if

  end subroutine set_stokes_shape_function_global

end module stokes_elements_generic_m

