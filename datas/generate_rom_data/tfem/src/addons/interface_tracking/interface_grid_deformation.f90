! Copyright (C) 2012-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system (e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Routines for performing interface grid deformation inspired by:

!   Decheng Wan and Stefan Turek, "Fictitious boundary and moving mesh methods
!   for the numerical simulation of rigid particulate flows", J. Comp. Phys.,
!   222 (2007) 28-56

!   M. Grajewski, M. Kester, and S. Turek., "Mathematical and numerical analysis
!   of a robust and efficient grid deformation method in the finite element
!   context", SIAM Journal on Scientific Computing, 31 (2009) 1539-1557

! Like in the original method a Poisson problem is solved for determining a
! grid velocity field by taking a gradient of the field. The main
! differences are:
!  1) The poisson problem is solved on the curvilinear system (curve,
!     surface).
!  2) The grid velocity field is determined by the surface gradient,
!     which is a tangential velocity.
!  3) The grid velocity is used to move the grid within a time step,
!     whereas the original method uses the grid velocity in a virtual
!     time stepping scheme to create the new grid from a static regular
!     grid.

module interface_grid_deformation_m

  use tfem_elem_m
  use interface_tracking_elements_m, only: set_globals_interface_tracking

  implicit none

contains


! Element routine to compute the g area function

  subroutine g_area_deriv ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use interface_tracking_globals_m

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

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( intpolv > 0 ) then
        call set_shape_function ( shapefuncv, xig, phiv, dphiv )
      else
        phiv = phi
        dphiv = dphi
      end if

    end if

!   get coordinates of nodes of interface

    call get_coordinates ( mesh, elgrp, elem, x )

    if ( isoshape == 0 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl )

    else if ( isoshape == 1 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphiv, dxdxis, surfl )

    end if

    elemvec = sum ( surfl * wg )

    elemwts = 1

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x )
      deallocate ( dxdxis )
      deallocate ( surfl )

    end subroutine deallocate_arrays

  end subroutine g_area_deriv


! Element routine to scale a function 1/f of a function f on a surface Gamma
! (given in the nodes):
!
!  c_f \int_Gamma 1/f dx = Gamma
!
! Computed are the elemvec(1) = integral 1/f, elemvec(2) = (area of) Gamma
! For use with "integrate" in postprocessing.

  subroutine integrate_inverse_elem ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec )

    use interface_tracking_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( intpolv > 0 ) then
        call set_shape_function ( shapefuncv, xig, phiv, dphiv )
      else
        phiv = phi
        dphiv = dphi
      end if

    end if

!   get coordinates of nodes of interface

    call get_coordinates ( mesh, elgrp, elem, x )

!   get f in the nodes

    call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )

    if ( isoshape == 0 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl )

!     1/f by interpolation of f from the nodes

      work = 1 / matmul ( phi, u )

    else if ( isoshape == 1 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphiv, dxdxis, surfl )

!     1/f by interpolation of f from the nodes

      work = 1 / matmul ( phiv, u )

    end if

!   integration of 1/f over the element

    elemvec(1) = sum ( work * surfl * wg )

!   area

    elemvec(2) = sum ( surfl * wg )

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )
      allocate ( work(ninti) )
      allocate ( u(max(ndf,ndfv)) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x )
      deallocate ( dxdxis )
      deallocate ( surfl )
      deallocate ( work )
      deallocate ( u )

    end subroutine deallocate_arrays

  end subroutine integrate_inverse_elem


! poisson element on a surface with additional elemvec 1/f-1/g:
!         2
!  - nabla  c = 1/f-1/g
!         s
!
  subroutine poisson_gd_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use interface_tracking_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( intpolv > 0 ) then
        call set_shape_function ( shapefuncv, xig, phiv, dphiv )
      else
        phiv = phi
        dphiv = dphi
      end if

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    if ( isoshape == 0 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, normal, &
        gij_up=gij_up )

      if ( vector ) then

!       1/f

        call get_vector ( mesh, problem, oldvectors%v(2)%p, elgrp, elem, u )
        work = 1/ matmul ( phi, u )

!       1/g

        call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )
        work = work - 1 / matmul ( phi, u )

      end if

    else if ( isoshape == 1 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphiv, dxdxis, surfl, normal, &
        gij_up=gij_up )

      if ( vector ) then

!       1/f

        call get_vector ( mesh, problem, oldvectors%v(2)%p, elgrp, elem, u )
        work = 1/ matmul ( phiv, u )

!       1/g

        call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )
        work = work - 1 / matmul ( phiv, u )

      end if

    end if

    if ( vector ) then

      elemvec = matmul ( work * surfl * wg, phi )

    end if

    if ( matrix ) then

      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphi(ip,i,:) * &
                               matmul ( gij_up(ip,:,:), dphi(ip,j,:) ) )

          end do
          elemmat(i,j) = sum ( work * surfl * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

    end if

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )
      allocate ( normal(ninti,ndim) )
      allocate ( work(ninti) )
      allocate ( gij_up(ninti,ndim-1,ndim-1) )
      allocate ( u(max(ndf,ndfv)) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x )
      deallocate ( dxdxis )
      deallocate ( surfl )
      deallocate ( normal )
      deallocate ( work )
      deallocate ( gij_up )
      deallocate ( u )

    end subroutine deallocate_arrays

  end subroutine poisson_gd_elem


! Curvature element on a surface/curve mesh.
! This is based on the surface tension element with surface tension = 1.
!
  subroutine curvature_gd_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use interface_tracking_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip, N

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( intpolv > 0 ) then
        call set_shape_function ( shapefuncv, xig, phiv, dphiv )
      else
        phiv = phi
        dphiv = dphi
      end if

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    if ( isoshape == 0 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, &
        gi_up=gi_up )

!     coordinates of integration points
      xg = matmul ( phi, x )

    else if ( isoshape == 1 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphiv, dxdxis, surfl, &
        gi_up=gi_up )

!     coordinates of integration points
      xg = matmul ( phiv, x )

    end if

    if ( coorsys == 1 ) then
      surfl = 2 * pi * xg(:,2) * surfl
    end if

    if ( vector ) then

!     - (nabla v)^T:(I-nn) = - v^k_j dphi_k/dxi_i g^i_j

      do j = 1, ndim
        do N = 1, ndf
          do ip = 1, ninti
            work(ip) = dot_product ( dphi(ip,N,:), gi_up(ip,j,:) )
          end do
          tmp(N,j) = sum ( work * surfl * wg )
        end do
      end do

      if ( coorsys == 1 ) then
        do N = 1, ndf
          work = phi(:,N) / xg(:,2)
          tmp(N,2) = tmp(N,2) + sum ( work * surfl * wg )
        end do
      end if

      elemvec = - reshape ( tmp, [ndf*ndim] )

    end if

    if ( matrix ) then

      do i = 1, ndf
        do j = i, ndf
          elemmat(i,j) = sum ( phi(:,i) * phi(:,j) * surfl * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

    end if

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )
      allocate ( work(ninti) )
      allocate ( gi_up(ninti,ndim,ndim-1) )
      allocate ( tmp(ndf,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x )
      deallocate ( xg )
      deallocate ( dxdxis )
      deallocate ( surfl )
      deallocate ( work )
      deallocate ( gi_up )
      deallocate ( tmp )

    end subroutine deallocate_arrays

  end subroutine curvature_gd_elem

end module interface_grid_deformation_m
