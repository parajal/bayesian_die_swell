
! Copyright (C) 2005-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

!
! Internal streamfunction element routine using the Poisson Equation:
!
!   - nabla^2 psi = omega
!
! where omega is derived from the velocity vector.

! Boundary element for a natural boundary on a curve for the streamfunction
! equation (Poisson equation):
!
!    dpsidn = -v * nx + u * ny   (= -tangential velocity)
!
! since dpsidx = -v, dpsidy = u.
!

! Use Dirichlet in a single point to set level of psi.


module streamfunction_elements_m

  use tfem_elem_m
  use stokes_set_globals_m

  implicit none


contains


! Internal streamfunction element routine using the Poisson Equation:
!
!   - nabla^2 psi = omega
!

  subroutine streamfunction_elem ( mesh, problem, elgrp, elem, matrix, vector, &
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


    integer :: i, j, ip


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!     allocate arrays

      allocate ( wg(ninti), detF(ninti), work(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( dudy(ninti), dvdx(ninti), u(2*nodalp), tmp(nodalp,2) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    call shape_derivative ( dphi, Finv, dphidx )

!   build equations

    if ( vector ) then

      call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u, &
        layer=layer )

      tmp = reshape ( u, [nodalp,ndim] )

      dudy = matmul ( dphidx(:,:,2), tmp(1:ndf,1) )
      dvdx = matmul ( dphidx(:,:,1), tmp(1:ndf,2) )

      if ( coorsys == 1 ) then
        work = matmul ( phi, u(1:ndf) )
        elemvec = 2 * pi * &
            matmul ( ( xg(:,2) * (dvdx - dudy) - work ) * detF * wg, phi )
      else
        elemvec = matmul ( (dvdx - dudy) * detF * wg, phi )
      end if

    end if

    if ( matrix ) then

      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,i,:) * dphidx(ip,j,:) )
          end do
          elemmat(i,j) = sum ( work * detF * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

    end if

    if ( last ) then

!     last element in this group

      deallocate ( Finv, dphidx )
      deallocate ( dphi, F )
      deallocate ( xg )
      deallocate ( xig, phi, x, work )
      deallocate ( wg, detF )
      deallocate ( dudy, dvdx, u, tmp )

    end if

  end subroutine streamfunction_elem


! Boundary element for a natural boundary on a curve for the streamfunction
! equation (Poisson equation):
!
!    dpsidn = -v * nx + u * ny
!
! since dpsidx = -v, dpsidy = u.
!

  subroutine streamfunction_natboun_curve ( mesh, problem, curve, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: ip


    if ( first ) then

!     first element on this curve

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, curve )

!     allocate arrays

      allocate ( wg(ninti), work(ninti), curvel(ninti), normal(ninti,ndim) )
      allocate ( xig(ninti,1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,1), dxdxi(ninti,ndim) )
      allocate ( u(2*ndf), gradpsi(ninti,ndim), tmp(ndf,2) )

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

    call get_vector_geometry ( mesh, problem, oldvectors%v(1)%p, elem, u, &
      curve=curve, layer=layer )

    tmp = reshape ( u, [ndf,2] )

    gradpsi(:,1) = - matmul ( phi, tmp(:,2) )
    gradpsi(:,2) =   matmul ( phi, tmp(:,1) )

    do ip = 1, ninti
      work(ip) = dot_product ( normal(ip,:), gradpsi(ip,:) )
    end do

    elemvec = matmul ( work * curvel * wg, phi )

    if ( last ) then

!     last element on this curve

      deallocate ( wg, work, curvel, normal )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi, dxdxi )
      deallocate ( u, gradpsi, tmp )

    end if

  end subroutine streamfunction_natboun_curve


! compute streamfunction in nodal points (when some nodes are missing)

  subroutine streamfunction_deriv ( mesh, problem, elgrp, elem, first, last, &
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

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )

      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf) )
      allocate ( u(ndf) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefunc, xrnod, phi )

    end if

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      layer=layer )

    elemvec = matmul ( phi, u(1:ndf) )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( xrnod, phi )
      deallocate ( u )

    end if

  end subroutine streamfunction_deriv

end module streamfunction_elements_m

