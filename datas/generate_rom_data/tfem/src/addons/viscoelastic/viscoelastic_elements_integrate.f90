
! Copyright (C) 2014-2024 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! This module contains element subroutines for integration of fields of a
! viscoelastic problem over the internal domain and the boundary.
! The subroutines are to be used with "integrate" and
! "integrate_boundary_elements" in the postprocessing module.
!
! The subroutines for the internal domain are:
!
! - viscoelastic_integrate_tau_tensor: integral of tau
!
! The subroutines for boundaries (geometries) are:
!
! - viscoelastic_integrate_tau_geometry: integral of (tau.n)x over a geometry
!

module viscoelastic_elements_integrate_m

  use stokes_elements_integrate_m
  use viscoelastic_elements_generic_m

  implicit none

contains


! Element routine for the integral of the viscoelastic extra-stress over the
! domain. For use with integrate.

  subroutine viscoelastic_integrate_tau_tensor ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: m, comp

    if ( first ) then

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic ( coefficients )

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

!     allocate arrays

      allocate ( wg(ninti) )
      allocate ( x(nodalp,ndim), xig(ninti,ndim) )
      allocate ( F(ninti,ndim,ndim), Finv(ninti,ndim,ndim), detF(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim), theta(ninti,ndfc) )
      allocate ( c(ndfc,ncomp), cg(ninti,ncomp,nmodes) )
      allocate ( tauvec(ninti,ncompt), tauten(ninti,ncompu,ncompu) )

!     set Gauss integration

      call set_Gauss_integration ( gauss, xig, wg )

!     set shape functions

      call set_shape_function ( shapefunc, xig, phi, dphi ) ! velocity
      call set_shape_function ( shapefuncc, xig, theta )    ! c tensor

      if ( vemcompressible ) then
        allocate ( psi(ninti,ndfp) )
        call set_shape_function ( shapefuncp, xig, psi ) ! pressure
        allocate ( pr(ndfp), press(ninti) )
        call create_vemopt ( vemopt, dep_J=.true., np=ninti )
      else
        call create_vemopt ( vemopt )
      end if

    end if

!   geometry dependent on space dimension

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      write(*,'(/a/a/)') 'Error viscoelastic_integrate_stress:', &
        ' Axisymmetric coordinates not yet implemented'
      stop
    end if

!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure ( mesh, oldvectors%p(1)%p, elgrp, elem, &
        coefficients, oldvectors%s(1)%p, vemopt )

    end if

!   get conformation tensor

    do m = mode1, mode2

!     single mode conformation
      call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
        cmode=c, mode=m )

      cg(:,:,m) = matmul ( theta, c )

    end do

    call stress_tensor_viscoelastic ( cg, vemopt=vemopt )

    do comp = 1, ncompt
      elemvec(comp) = sum( tauvec(:,comp) * detF * wg )
    end do

    if ( last ) then

!     last element on this curve

      call delete ( vemodel )

      deallocate ( wg )
      deallocate ( x, xig )
      deallocate ( F, Finv, detF )
      deallocate ( phi, dphi, theta )
      deallocate ( c, cg )
      deallocate ( tauvec, tauten )

      if ( vemcompressible ) then
        deallocate ( psi )
        deallocate ( pr, press )
      end if

      call delete ( vemopt )

    end if

  end subroutine viscoelastic_integrate_tau_tensor


! Element routine for the integral of the (tau.n)x over the geometry
! For use with integrate_boundary_elements.

  subroutine viscoelastic_integrate_tau_geometry ( mesh, problem, geometry, &
    elem, first, last, coefficients, oldvectors, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: geometry, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: ip, m, i, j

    if ( first ) then

!     first element on this geometry

      call check ( coefficients, 'viscoelastic_integrate_stress_geometry', &
        indexarray=[36], minimum=[-1], maximum=[1] )

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, &
        ndimr=mesh%ndim-1, geometry=geometry )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic_boun ( coefficients )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,ndim-1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim-1), dxdxis(ninti,ndim,ndim-1) )
      allocate ( work(ndf*ncomp), c(ndf,ncomp), cg(ninti,ncomp,nmodes) )
      allocate ( work4(ninti,ndim), tmp(ndim,ndim) )
      allocate ( tauten(ninti,ndim,ndim), tauvec(ninti,ncompt) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( globalshape, xig, wg )

      call set_shape_function ( globalshape, xig, phi, dphi )

      if ( vemcompressible ) then
        allocate ( psi(ninti,ndfp) )
        call set_shape_function ( shapefuncp, xig, psi )
        allocate ( pr(ndfp), press(ninti) )
        call create_vemopt ( vemopt, dep_J=.true., np=ninti )
      else
        call create_vemopt ( vemopt )
      end if

    end if

!   geometry dependent on space dimension

    call get_coordinates_geometry ( mesh, elem, x, ndimr=ndim-1, &
      geometry=geometry )

    call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, normal )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      write(*,'(/a/a/)') 'viscoelastic_integrate_stress_geometry:', &
        ' Axisymmetric coordinates not yet implemented'
      stop
    end if

!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure_geometry ( mesh, problem, elem, coefficients, &
        oldvectors%s(1)%p, vemopt, ndimr=ndim-1, geometry=geometry )

    end if

!   get conformation tensor

    do m = mode1, mode2

      call get_vector_geometry ( mesh, problem, oldvectors%v1(1)%p(m), elem, &
        work, ndimr=ndim-1, geometry=geometry, layer=layer )

      c = reshape ( work, [ndf,ncomp] )

      cg(:,:,m) = matmul ( phi, c )

    end do

!   get stress tensor
!   WARNING: since oldvectors%v1(1) contains the conformation tensor
!   derived from "deriv_conformation_tensor" subroutine, the C tensor values
!   are already in non-log formulation, independently of logc = 0 or logc = 1.
!   So, independently of the formulation used, the "stress_viscoelastic_2D/3D"
!   has to called to get the stress tensor

    if ( ndim == 2 ) then

      call stress_viscoelastic_2D ( vemodel, cg, tauvec, mode1=mode1, &
        mode2=mode2, vemopt=vemopt )

      tauten(:,1,1) = tauvec(:,1)
      tauten(:,1,2) = tauvec(:,2)
      tauten(:,2,1) = tauvec(:,2)
      tauten(:,2,2) = tauvec(:,3)

    else if ( ndim == 3 ) then

      call stress_viscoelastic_3D ( vemodel, cg, tauvec, mode1=mode1, &
        mode2=mode2, vemopt=vemopt )

      tauten(:,1,1) = tauvec(:,1)
      tauten(:,1,2) = tauvec(:,2)
      tauten(:,1,3) = tauvec(:,3)
      tauten(:,2,1) = tauvec(:,2)
      tauten(:,2,2) = tauvec(:,4)
      tauten(:,2,3) = tauvec(:,5)
      tauten(:,3,1) = tauvec(:,3)
      tauten(:,3,2) = tauvec(:,5)
      tauten(:,3,3) = tauvec(:,6)

    end if

!   traction force in each integration point: t=tau.n

    do ip = 1, ninti
      work4(ip,:) = matmul( tauten(ip,:,:), normal(ip,:) )
    end do

!   evaluate in each integration point: (tau.n)x = tx and integrate

    do i = 1, ndim
      do j = 1, ndim
        tmp(i,j) = sum ( work4(:,i) * xg(:,j) * surfl * wg )
      end do
    end do

!   element vector

    elemvec = reshape ( transpose(tmp), [ndim**2] )

    if ( last ) then

!     last element on this curve

      call delete ( vemodel )

      deallocate ( wg, surfl )
      deallocate ( normal )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi, dxdxis )
      deallocate ( work, c, cg )
      deallocate ( work4, tmp )
      deallocate ( tauten, tauvec )

      if ( vemcompressible ) then
        deallocate ( psi )
        deallocate ( pr, press )
      end if

      call delete ( vemopt )

    end if

  end subroutine viscoelastic_integrate_tau_geometry

end module viscoelastic_elements_integrate_m
