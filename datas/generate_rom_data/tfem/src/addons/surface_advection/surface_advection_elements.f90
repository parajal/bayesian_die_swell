
! Copyright (C) 2011-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the surface advection.
!
! Surfaces are described by a single height function H(t,x) in 2D flow or
! H(t,x,y) in 3D flow.
!
! The condition that the surface is material (fluid particles remain on the
! surface) leads to the kinematical condition:
!
!   dH      dH
!   -- + ux -- = uy    for 2D flow
!   dt      dx
!
! or
!
!   dH      dH      dH
!   -- + ux -- + uy -- = uz    for 3D flow
!   dt      dx      dy
!
! where u=(ux,uy) or u=(ux,uy,uz) the velocity vector at the surface.
! Note, that (x,H) or (x,y,H) form an orthogonal system.
!
! In polar coordinates a material free surface can be described by
! the kinematical condition:
!
!   dH   u_theta  dH
!   -- +   --     --    = u_r   for 2D flow
!   dt     r    dtheta
!
! Note, that the height function H is defined as the distance from the
! origin of the polar coordinate system, i.e. H=r. The velocity vector
! that is supplied to the element is (u_theta,u_r), i.e. the theta direction
! is the first component.
!
! For 3D flow the height function in polar coordinates reads:
!
! dH   u_theta dH           dH
! -- +  ----   ---   + u_z  -- = u_r  for 3D flow
! dt      r   dtheta        dz
!
! Note, that the height function H is defined as the distance from the
! origin of the polar coordinate system, i.e. H=r. The velocity vector
! that is supplied to the element is (u_z,u_theta,u_r).
!
! For material lines that swell in 2 directions the condition that the
! lines are material leads to the kinematical condition:
!
!  dH      dH
!  -- + ux -- = u2D  for 3D flow
!  dt      dx
!
! where u2D=(uy,uz) the velocity at the line and H=(hy,hz).

module surface_advection_elements_m

  use tfem_elem_m
  use supg_utils_m

  implicit none

contains


! Internal element routine (generic element for 1D and 2D).

  subroutine surface_advection_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use surface_advection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   select depending on integration scheme

    select case ( coefficients%i(5) )
      case(1) ! first-order
        call surface_advection_elem1 ( mesh, problem, elgrp, elem, &
          matrix, vector, first, last, coefficients, oldvectors, elemmat, &
          elemvec )
      case(2) ! second-order
        call surface_advection_elem2 ( mesh, problem, elgrp, elem, &
          matrix, vector, first, last, coefficients, oldvectors, elemmat, &
          elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in surface_advection_elem:', &
        ' incorrect time integration scheme  = ', coefficients%i(5)
        stop
    end select

  end subroutine surface_advection_elem


! Internal element routine (generic element for 1D and 2D).
! First-order time integration.

  subroutine surface_advection_elem1 ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use surface_advection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    logical :: checkzero
    integer :: method, i, j, ip
    real(dp) :: deltat, beta, h, esize, Uglobal
    real(dp), parameter :: eps=3*tiny(1._dp)

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_surface_advection ( mesh, coefficients, elgrp )

      if ( ndim > 2 ) then
        write(*,'(/2(a/))') 'Error in surface_advection_elem1:', &
          ' Space dimension of mesh ndim must be 1 or 2'
        stop
      end if

      ncompv = ndim + 1

!     allocate arrays

      call allocate_arrays

!     Gauss rule and shape function

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
      call isoparametric_deformation ( x, dphi, F, Finv, detF )
      xg = matmul ( phi, x )
    else if ( isoshape == 1 ) then
      call isoparametric_deformation ( x, dphiv, F, Finv, detF )
      xg = matmul ( phiv, x )
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   velocity vector at tn+1

    call evaluate_vector_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(2), value=coefficients%r(1:ncompv), &
      vfunc=coefficients%vfunc, vfuncnr=coefficients%i(3), x=xg, &
      indx_v=1, layer=0, phi=phiv, indx_e=1, coef=uvecg )

!   height function at tn

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u=Hn )

    Hng = matmul ( phi, Hn )

!   u.gradphi operator

    do ip = 1, ninti
      ugradphi(ip,:) = matmul ( dphidx(ip,:,:), uvecg(ip,1:ndim) )
    end do

!   coefficients

    method = coefficients%i(4)

    if ( method == 1 ) then

!     SUPG

      checkzero = coefficients%i(14) == 1
      Uglobal = coefficients%r(6)

      if ( globalshape == 'line' ) then

!       compute h/U for a line element

        h = abs(x(nodalp,1)-x(1,1))

        beta = coefficients%r(5)

        if ( checkzero .and. any ( abs(uvecg(:,1)) < eps ) ) then
          if ( Uglobal < eps ) then
            write(*,'(/a/a/)') 'Error in surface_advection_elem1:', &
              ' Uglobal negative or near zero '
            stop
          else
            tau = beta * h / 2 / Uglobal
          end if
        else
          tau = beta * h / 2 / abs(uvecg(:,1))
        end if

      else

!       compute h/U

        if ( checkzero ) esize = sum ( detF * wg )

        call supg_hU ( ndim, globalshape, hoverU, htype=2, hlocation=2, &
          Uscaling=3, uvec=uvecg(:,1:ndim), Finv=Finv, esize=esize, &
          Uglobal=Uglobal, checkzero=checkzero )

        beta = coefficients%r(5)

        tau = beta * hoverU / 2  ! upwinding parameter

      end if

    end if

    deltat = coefficients%r(4)


    if ( matrix ) then

      if ( method == 0 ) then

!       Galerkin
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( phi(:,i) * &
              ( phi(:,j) / deltat + ugradphi(:,j) ) * detF * wg )
          end do
        end do

      else if ( method == 1 ) then

!       SUPG
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( ( phi(:,i) + tau * ugradphi(:,i)) * &
              ( phi(:,j) / deltat + ugradphi(:,j) ) * detF * wg )
          end do
        end do

      end if

    end if

    if ( vector ) then

      if ( method == 0 ) then

!       Galerkin
        do i = 1, ndf
          elemvec(i) = &
            sum ( phi(:,i) * ( Hng / deltat + uvecg(:,ndim+1) ) * detF * wg )
        end do

      else if ( method == 1 ) then

!       SUPG
        do i = 1, ndf
          elemvec(i) = sum ( ( phi(:,i) +  tau * ugradphi(:,i) ) * &
                             ( Hng / deltat + uvecg(:,ndim+1) ) * detF * wg )
        end do

      end if

    end if

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim) )
      allocate ( x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( F(ninti,ndim,ndim), Finv(ninti,ndim,ndim) )
      allocate ( fg(ninti), detF(ninti) )
      allocate ( dphidx(ninti,ndf,ndim) )
      allocate ( Hn(ndf), Hng(ninti) )
      allocate ( uvecg(ninti,ncompv) )
      allocate ( tau(ninti), hoverU(ninti) )
      allocate ( ugradphi(ninti,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x, xg )
      deallocate ( F, Finv )
      deallocate ( fg, detF )
      deallocate ( dphidx )
      deallocate ( Hn, Hng )
      deallocate ( uvecg )
      deallocate ( tau, hoverU )
      deallocate ( ugradphi )

    end subroutine deallocate_arrays

  end subroutine surface_advection_elem1


! Internal element routine (generic element for 1D and 2D).
! Second-order time integration.

  subroutine surface_advection_elem2 ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use surface_advection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    logical :: checkzero
    integer :: method, i, j, ip
    real(dp) :: deltat, beta, h, esize, Uglobal
    real(dp), parameter :: eps=3*tiny(1._dp)


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_surface_advection ( mesh, coefficients, elgrp )

      if ( ndim > 2 ) then
        write(*,'(/2(a/))') 'Error in surface_advection_elem2:', &
          ' Space dimension of mesh ndim must be 1 or 2'
        stop
      end if

      ncompv = ndim + 1

!     allocate arrays

      call allocate_arrays

!     Gauss rule and shape function

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
      call isoparametric_deformation ( x, dphi, F, Finv, detF )
      xg = matmul ( phi, x )
    else if ( isoshape == 1 ) then
      call isoparametric_deformation ( x, dphiv, F, Finv, detF )
      xg = matmul ( phiv, x )
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   velocity vector at tn+1

    call evaluate_vector_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(2), value=coefficients%r(1:ncompv), &
      vfunc=coefficients%vfunc, vfuncnr=coefficients%i(3), x=xg, &
      indx_v=1, layer=0, phi=phiv, indx_e=1, coef=uvecg )

!   height function at tn

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u=Hn )

    Hng = matmul ( phi, Hn )

!   height function at tn-1

    call get_sysvector ( mesh, problem, oldvectors%s(2)%p, elgrp, elem, u=Hnm1 )

    Hnm1g = matmul ( phi, Hnm1 )

!   u.gradphi operator

    do ip = 1, ninti
      ugradphi(ip,:) = matmul ( dphidx(ip,:,:), uvecg(ip,1:ndim) )
    end do

!   coefficients

    method = coefficients%i(4)

    if ( method == 1 ) then

!     SUPG

      checkzero = coefficients%i(14) == 1
      Uglobal = coefficients%r(6)

      if ( globalshape == 'line' ) then

!       compute h/U for a line element

        h = abs(x(nodalp,1)-x(1,1))

        beta = coefficients%r(5)

        if ( checkzero .and. any ( abs(uvecg(:,1)) < eps ) ) then
          if ( Uglobal < eps ) then
            write(*,'(/a/a/)') 'Error in surface_advection_elem2:', &
              ' Uglobal negative or near zero '
            stop
          else
            tau = beta * h / 2 / Uglobal
          end if
        else
          tau = beta * h / 2 / abs(uvecg(:,1))
        end if

      else

!       compute h/U

        if ( checkzero ) esize = sum ( detF * wg )

        call supg_hU ( ndim, globalshape, hoverU, htype=2, hlocation=2, &
          Uscaling=3, uvec=uvecg(:,1:ndim), Finv=Finv, esize=esize, &
          Uglobal=Uglobal, checkzero=checkzero )

        beta = coefficients%r(5)

        tau = beta * hoverU / 2  ! upwinding parameter

      end if

    end if

    deltat = coefficients%r(4)

    if ( matrix ) then

      if ( method == 0 ) then

!       Galerkin
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( phi(:,i) * &
              ( 1.5_dp * phi(:,j) / deltat + ugradphi(:,j) ) * detF * wg )
          end do
        end do

      else if ( method == 1 ) then

!       SUPG
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( ( phi(:,i) + tau * ugradphi(:,i)) * &
              ( 1.5_dp * phi(:,j) / deltat + ugradphi(:,j) ) * detF * wg )
          end do
        end do

      end if

    end if

    if ( vector ) then

      if ( method == 0 ) then

!       Galerkin
        do i = 1, ndf
          elemvec(i) = sum ( phi(:,i) * &
             ( ( 2 * Hng - 0.5_dp * Hnm1g ) / deltat + uvecg(:,ndim+1) ) &
                          * detF * wg )
        end do

      else if ( method == 1 ) then

!       SUPG
        do i = 1, ndf
          elemvec(i) = sum ( ( phi(:,i) +  tau * ugradphi(:,i) ) * &
             ( ( 2 * Hng - 0.5_dp * Hnm1g ) / deltat + uvecg(:,ndim+1) ) &
                          * detF * wg )
        end do

      end if

    end if

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim) )
      allocate ( x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( F(ninti,ndim,ndim), Finv(ninti,ndim,ndim) )
      allocate ( fg(ninti), detF(ninti) )
      allocate ( dphidx(ninti,ndf,ndim) )
      allocate ( Hn(ndf), Hng(ninti) )
      allocate ( Hnm1(ndf), Hnm1g(ninti) )
      allocate ( uvecg(ninti,ncompv) )
      allocate ( tau(ninti), hoverU(ninti) )
      allocate ( ugradphi(ninti,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x, xg )
      deallocate ( F, Finv )
      deallocate ( fg, detF )
      deallocate ( dphidx )
      deallocate ( Hn, Hng )
      deallocate ( Hnm1, Hnm1g )
      deallocate ( uvecg )
      deallocate ( tau, hoverU )
      deallocate ( ugradphi )

    end subroutine deallocate_arrays

  end subroutine surface_advection_elem2


! element routine to calculate positions of material lines

  subroutine surface_material_line_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use surface_advection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   select depending on integration scheme

    select case ( coefficients%i(5) )
      case(1) ! first-order
        call surface_material_line_elem1 ( mesh, problem, elgrp, elem, &
          matrix, vector, first, last, coefficients, oldvectors, elemmat, &
          elemvec )
      case(2) ! second-order
        call surface_material_line_elem2 ( mesh, problem, elgrp, elem, &
          matrix, vector, first, last, coefficients, oldvectors, elemmat, &
          elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in surface_material_line_elem:', &
        ' incorrect time integration scheme  = ', coefficients%i(5)
        stop
    end select

  end subroutine surface_material_line_elem


! Internal element routine material lines.
! First-order time integration.

  subroutine surface_material_line_elem1 ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use surface_advection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    logical :: checkzero
    integer :: method, i, j, ip
    integer, parameter :: nrhsdline = 2
    real(dp) :: deltat, beta, h, Uglobal
    real(dp), parameter :: eps=3*tiny(1._dp)


    if ( first ) then

!     first element in this group

!     set globals

      ndim = mesh%element(elgrp)%ndim
      nodalp = mesh%element(elgrp)%numnod
      ninti = coefficients%i(1)

      call set_globals_surface_advection ( mesh, coefficients, elgrp )

      if ( ndim > 1 ) then
        write(*,'(/2(a/))') 'Error in surface_material_line_elem1:', &
          ' Space dimension of mesh ndim must be 1'
        stop
      end if

      ncompv = 3

!     allocate arrays

      call allocate_arrays

!     Gauss rule and shape function

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
      call isoparametric_deformation ( x, dphi, F, Finv, detF )
      xg = matmul ( phi, x )
    else if ( isoshape == 1 ) then
      call isoparametric_deformation ( x, dphiv, F, Finv, detF )
      xg = matmul ( phiv, x )
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   velocity vector at tn+1

    call evaluate_vector_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(2), value=coefficients%r(1:ncompv), &
      vfunc=coefficients%vfunc, vfuncnr=coefficients%i(3), x=xg, &
      indx_v=1, layer=0, phi=phiv, indx_e=1, coef=uvecg )

!   height function at tn

    do i = 1, nrhsdline
      call get_sysvector ( mesh, problem, oldvectors%s1(1)%p(i), elgrp, elem, &
                           u=Hn )
      Hngl(:,i) = matmul ( phi, Hn )
    end do

!   u.gradphi operator

    do ip = 1, ninti
      ugradphi(ip,:) = matmul ( dphidx(ip,:,:), uvecg(ip,1:ndim) )
    end do

!   coefficients

    method = coefficients%i(4)

    if ( method == 1 ) then

!     SUPG

      checkzero = coefficients%i(14) == 1
      Uglobal = coefficients%r(6)

!     compute h/U for a line element

      h = abs(x(nodalp,1)-x(1,1))

      beta = coefficients%r(5)

      if ( checkzero .and. any ( abs(uvecg(:,1)) < eps ) ) then
        if ( Uglobal < eps ) then
          write(*,'(/a/a/)') 'Error in surface_material_line_elem1:', &
            ' Uglobal negative or near zero '
          stop
        else
          tau = beta * h / 2 / Uglobal
        end if
      else
        tau = beta * h / 2 / abs(uvecg(:,1))
      end if

    end if

    deltat = coefficients%r(4)


    if ( matrix ) then

      if ( method == 0 ) then

!       Galerkin
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( phi(:,i) * &
              ( phi(:,j) / deltat + ugradphi(:,j) ) * detF * wg )
          end do
        end do

      else if ( method == 1 ) then

!       SUPG
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( ( phi(:,i) + tau * ugradphi(:,i)) * &
              ( phi(:,j) / deltat + ugradphi(:,j) ) * detF * wg )
          end do
        end do

      end if

    end if

    if ( vector ) then

      if ( method == 0 ) then

!       Galerkin
        do j = 1,nrhsdline
          do i = 1, ndf
            work(i,j) = &
              sum ( phi(:,i) * ( Hngl(:,j) / deltat + uvecg(:,j+1) ) * &
                    detF * wg )
          end do
        end do

      else if ( method == 1 ) then

!       SUPG
        do j = 1,nrhsdline
          do i = 1, ndf
            work(i,j) = sum ( ( phi(:,i) +  tau * ugradphi(:,i) ) * &
                             ( Hngl(:,j) / deltat + uvecg(:,j+1) ) * &
                               detF * wg )
          end do
        end do

      end if

      elemvec = reshape ( work, [ndf * nrhsdline] )

    end if

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim) )
      allocate ( x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( F(ninti,ndim,ndim), Finv(ninti,ndim,ndim) )
      allocate ( fg(ninti), detF(ninti) )
      allocate ( dphidx(ninti,ndf,ndim) )
      allocate ( Hn(ndf) )
      allocate ( Hngl(ninti,nrhsdline) )
      allocate ( work(ndf,nrhsdline) )
      allocate ( uvecg(ninti,ncompv) )
      allocate ( tau(ninti), hoverU(ninti) )
      allocate ( ugradphi(ninti,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x, xg )
      deallocate ( F, Finv )
      deallocate ( fg, detF )
      deallocate ( dphidx )
      deallocate ( Hn )
      deallocate ( Hngl )
      deallocate ( work )
      deallocate ( uvecg )
      deallocate ( tau, hoverU )
      deallocate ( ugradphi )

    end subroutine deallocate_arrays

  end subroutine surface_material_line_elem1


! Internal element routine material line.
! Second-order time integration.

  subroutine surface_material_line_elem2 ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use surface_advection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    logical :: checkzero
    integer :: method, i, j, ip
    integer, parameter :: nrhsdline=2
    real(dp) :: deltat, beta, h, Uglobal
    real(dp), parameter :: eps=3*tiny(1._dp)


    if ( first ) then

!     first element in this group

!     set globals

      ndim = mesh%element(elgrp)%ndim
      nodalp = mesh%element(elgrp)%numnod
      ninti = coefficients%i(1)

      call set_globals_surface_advection ( mesh, coefficients, elgrp )

      if ( ndim > 1 ) then
        write(*,'(/2(a/))') 'Error in surface_material_line_elem2:', &
          ' Space dimension of mesh ndim must be 1'
        stop
      end if

      ncompv = 3

!     allocate arrays

      call allocate_arrays

!     Gauss rule and shape function

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
      call isoparametric_deformation ( x, dphi, F, Finv, detF )
      xg = matmul ( phi, x )
    else if ( isoshape == 1 ) then
      call isoparametric_deformation ( x, dphiv, F, Finv, detF )
      xg = matmul ( phiv, x )
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   velocity vector at tn+1

    call evaluate_vector_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(2), value=coefficients%r(1:ncompv), &
      vfunc=coefficients%vfunc, vfuncnr=coefficients%i(3), x=xg, &
      indx_v=1, layer=0, phi=phiv, indx_e=1, coef=uvecg )

!   height function at tn

    do i = 1,nrhsdline

      call get_sysvector ( mesh, problem, oldvectors%s1(1)%p(i), elgrp, &
                           elem, u=Hn )

      Hngl(:,i) = matmul ( phi, Hn )

    end do

!   height function at tn-1

    do i = 1, nrhsdline

      call get_sysvector ( mesh, problem, oldvectors%s1(2)%p(i), elgrp, &
                           elem, u=Hnm1 )

      Hnm1gl(:,i) = matmul ( phi, Hnm1 )

    end do

!   u.gradphi operator

    do ip = 1, ninti
      ugradphi(ip,:) = matmul ( dphidx(ip,:,:), uvecg(ip,1:ndim) )
    end do

!   coefficients

    method = coefficients%i(4)

    if ( method == 1 ) then

!     SUPG

      checkzero = coefficients%i(14) == 1
      Uglobal = coefficients%r(6)

!     compute h/U for a line element

      h = abs(x(nodalp,1)-x(1,1))

      beta = coefficients%r(5)

      if ( checkzero .and. any ( abs(uvecg(:,1)) < eps ) ) then
        if ( Uglobal < eps ) then
          write(*,'(/a/a/)') 'Error in surface_material_line_elem2:', &
            ' Uglobal negative or near zero '
          stop
        else
          tau = beta * h / 2 / Uglobal
        end if
      else
        tau = beta * h / 2 / abs(uvecg(:,1))
      end if

    end if

    deltat = coefficients%r(4)

    if ( matrix ) then

      if ( method == 0 ) then

!       Galerkin
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( phi(:,i) * &
              ( 1.5_dp * phi(:,j) / deltat + ugradphi(:,j) ) * detF * wg )
          end do
        end do

      else if ( method == 1 ) then

!       SUPG
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( ( phi(:,i) + tau * ugradphi(:,i)) * &
              ( 1.5_dp * phi(:,j) / deltat + ugradphi(:,j) ) * detF * wg )
          end do
        end do

      end if

    end if

    if ( vector ) then

      if ( method == 0 ) then

!       Galerkin
        do j = 1,nrhsdline
          do i = 1, ndf
            work(i,j) = sum ( phi(:,i) * &
                   ( ( 2 * Hngl(:,j) - 0.5_dp * Hnm1gl(:,j) ) / deltat + &
                       uvecg(:,j+1) ) * detF * wg )
          end do
        end do

      else if ( method == 1 ) then

!       SUPG
        do j = 1,nrhsdline
          do i = 1, ndf
            work(i,j) = sum ( ( phi(:,i) +  tau * ugradphi(:,i) ) * &
             ( ( 2 * Hngl(:,j) - 0.5_dp * Hnm1gl(:,j) ) / deltat + &
                 uvecg(:,j+1) ) * detF * wg )
          end do
        end do

      end if

      elemvec = reshape ( work, [ndf * 1 * nrhsdline ] )

    end if

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim) )
      allocate ( x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( F(ninti,ndim,ndim), Finv(ninti,ndim,ndim) )
      allocate ( fg(ninti), detF(ninti) )
      allocate ( dphidx(ninti,ndf,ndim) )
      allocate ( Hn(ndf))
      allocate ( Hnm1(ndf))
      allocate ( Hngl(ninti,nrhsdline), Hnm1gl(ninti,nrhsdline) )
      allocate ( work(ndf,nrhsdline) )
      allocate ( uvecg(ninti,ncompv) )
      allocate ( tau(ninti), hoverU(ninti) )
      allocate ( ugradphi(ninti,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x, xg )
      deallocate ( F, Finv )
      deallocate ( fg, detF )
      deallocate ( dphidx )
      deallocate ( Hn )
      deallocate ( Hnm1 )
      deallocate ( Hngl, Hnm1gl )
      deallocate ( work )
      deallocate ( uvecg )
      deallocate ( tau, hoverU )
      deallocate ( ugradphi )

    end subroutine deallocate_arrays

  end subroutine surface_material_line_elem2


! Internal element routine for polar coordinates

  subroutine surface_advection_polar_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use surface_advection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   select depending on integration scheme

    select case ( coefficients%i(5) )
      case(1) ! first-order
        call surface_advection_polar_elem1 ( mesh, problem, elgrp, elem, &
          matrix, vector, first, last, coefficients, oldvectors, elemmat, &
          elemvec )
      case(2) ! second-order
        call surface_advection_polar_elem2 ( mesh, problem, elgrp, elem, &
          matrix, vector, first, last, coefficients, oldvectors, elemmat, &
          elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in surface_advection_polar_elem:', &
        ' incorrect time integration scheme  = ', coefficients%i(5)
        stop
    end select

  end subroutine surface_advection_polar_elem


! Internal element routine for polar coordinates
! First-order time integration.

  subroutine surface_advection_polar_elem1 ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use surface_advection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    logical :: checkzero
    integer :: method, i, j, ip
    real(dp) :: deltat, beta, h, esize, Uglobal
    real(dp), parameter :: eps=3*tiny(1._dp)

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_surface_advection ( mesh, coefficients, elgrp )

      if ( ndim /= 1 ) then
        write(*,'(/2(a/))') 'Error in surface_advection_polar_elem1:', &
          ' Space dimension of mesh ndim must be 1'
        stop
      end if

      ncompv = ndim + 1

!     allocate arrays

      call allocate_arrays

!     Gauss rule and shape function

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
      call isoparametric_deformation ( x, dphi, F, Finv, detF )
      xg = matmul ( phi, x )
    else if ( isoshape == 1 ) then
      call isoparametric_deformation ( x, dphiv, F, Finv, detF )
      xg = matmul ( phiv, x )
    end if


    call shape_derivative ( dphi, Finv, dphidx )

!   velocity vector at tn+1

    call evaluate_vector_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(2), value=coefficients%r(1:ncompv), &
      vfunc=coefficients%vfunc, vfuncnr=coefficients%i(3), x=xg, &
      indx_v=1, layer=0, phi=phiv, indx_e=1, coef=uvecg )

!   height function at tn

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u=Hn )

    Hng = matmul ( phi, Hn )

    do ip = 1,ninti
      uvecg(ip,1) = uvecg(ip,ndim) / Hng(ip)  ! uvecg(ip,1) = u_theta / r
    end do

!   u.gradphi operator

    do ip = 1, ninti
      ugradphi(ip,:) = matmul ( dphidx(ip,:,:), uvecg(ip,1:ndim) )
    end do

!   coefficients

    method = coefficients%i(4)

    if ( method == 1 ) then

!     SUPG

      checkzero = coefficients%i(14) == 1
      Uglobal = coefficients%r(6)

      if ( globalshape == 'line' ) then

!       compute h/U for a line element

        h = abs(x(nodalp,1)-x(1,1))

        beta = coefficients%r(5)

        if ( checkzero .and. any ( abs(uvecg(:,1)) < eps ) ) then
          if ( Uglobal < eps ) then
            write(*,'(/a/a/)') 'Error in surface_advection_elem1:', &
              ' Uglobal negative or near zero '
            stop
          else
            tau = beta * h / 2 / Uglobal
          end if
        else
          tau = beta * h / 2 / abs(uvecg(:,1))
        end if

      else

!       compute h/U

        if ( checkzero ) esize = sum ( detF * wg )

        call supg_hU ( ndim, globalshape, hoverU, htype=2, hlocation=2, &
          Uscaling=3, uvec=uvecg(:,1:ndim), Finv=Finv, esize=esize, &
          Uglobal=Uglobal, checkzero=checkzero )

        beta = coefficients%r(5)

        tau = beta * hoverU / 2  ! upwinding parameter

      end if

    end if

    deltat = coefficients%r(4)

    if ( matrix ) then

      if ( method == 0 ) then

!       Galerkin
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( phi(:,i) * &
              ( phi(:,j) / deltat + ugradphi(:,j) ) * detF * wg )
          end do
        end do

      else if ( method == 1 ) then

!       SUPG
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( ( phi(:,i) + tau * ugradphi(:,i)) * &
              ( phi(:,j) / deltat + ugradphi(:,j) ) * detF * wg )
          end do
        end do

      end if

    end if

    if ( vector ) then

      if ( method == 0 ) then

!       Galerkin
        do i = 1, ndf
          elemvec(i) = &
            sum ( phi(:,i) * ( Hng / deltat + uvecg(:,ndim+1) ) * detF * wg )
        end do

      else if ( method == 1 ) then

!       SUPG
        do i = 1, ndf
          elemvec(i) = sum ( ( phi(:,i) +  tau * ugradphi(:,i) ) * &
                             ( Hng / deltat + uvecg(:,ndim+1) ) * detF * wg )
        end do

      end if

    end if

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim) )
      allocate ( x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( F(ninti,ndim,ndim), Finv(ninti,ndim,ndim) )
      allocate ( fg(ninti), detF(ninti) )
      allocate ( dphidx(ninti,ndf,ndim) )
      allocate ( Hn(ndf), Hng(ninti) )
      allocate ( uvecg(ninti,ncompv) )
      allocate ( tau(ninti), hoverU(ninti) )
      allocate ( ugradphi(ninti,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x, xg )
      deallocate ( F, Finv )
      deallocate ( fg, detF )
      deallocate ( dphidx )
      deallocate ( Hn, Hng )
      deallocate ( uvecg )
      deallocate ( tau, hoverU )
      deallocate ( ugradphi )

    end subroutine deallocate_arrays

  end subroutine surface_advection_polar_elem1


! Internal element routine (generic element for 1D and 2D).
! Second-order time integration.

  subroutine surface_advection_polar_elem2 ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use surface_advection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    logical :: checkzero
    integer :: method, i, j, ip
    real(dp) :: deltat, beta, esize, h, Uglobal
    real(dp), parameter :: eps=3*tiny(1._dp)

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_surface_advection ( mesh, coefficients, elgrp )

      if ( ndim /= 1 ) then
        write(*,'(/2(a/))') 'Error in surface_advection_polar_elem2:', &
          ' Space dimension of mesh ndim must be 1'
        stop
      end if

      ncompv = ndim + 1

!     allocate arrays

      call allocate_arrays

!     Gauss rule and shape function

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
      call isoparametric_deformation ( x, dphi, F, Finv, detF )
      xg = matmul ( phi, x )
    else if ( isoshape == 1 ) then
      call isoparametric_deformation ( x, dphiv, F, Finv, detF )
      xg = matmul ( phiv, x )
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   velocity vector at tn+1

    call evaluate_vector_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(2), value=coefficients%r(1:ncompv), &
      vfunc=coefficients%vfunc, vfuncnr=coefficients%i(3), x=xg, &
      indx_v=1, layer=0, phi=phiv, indx_e=1, coef=uvecg )

!   height function at tn

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u=Hn )

    Hng = matmul ( phi, Hn )

!   height function at tn-1

    call get_sysvector ( mesh, problem, oldvectors%s(2)%p, elgrp, elem, u=Hnm1 )

    Hnm1g = matmul ( phi, Hnm1 )

    do ip = 1,ninti
!     uvecg(ip,1) = u_theta / r
      uvecg(ip,1) = uvecg(ip,ndim) / ( 2*Hng(ip)-Hnm1g(ip) )
    end do

!   u.gradphi operator

    do ip = 1, ninti
      ugradphi(ip,:) = matmul ( dphidx(ip,:,:), uvecg(ip,1:ndim) )
    end do

!   coefficients

    method = coefficients%i(4)

    if ( method == 1 ) then

!     SUPG

      checkzero = coefficients%i(14) == 1
      Uglobal = coefficients%r(6)

!     compute h/U for a line element

      if ( globalshape == 'line' ) then

!       compute h/U for a line element

        h = abs(x(nodalp,1)-x(1,1))

        beta = coefficients%r(5)

        if ( checkzero .and. any ( abs(uvecg(:,1)) < eps ) ) then
          if ( Uglobal < eps ) then
            write(*,'(/a/a/)') 'Error in surface_advection_elem1:', &
              ' Uglobal negative or near zero '
            stop
          else
            tau = beta * h / 2 / Uglobal
          end if
        else
          tau = beta * h / 2 / abs(uvecg(:,1))
        end if

      else

!       compute h/U

        if ( checkzero ) esize = sum ( detF * wg )

        call supg_hU ( ndim, globalshape, hoverU, htype=2, hlocation=2, &
          Uscaling=3, uvec=uvecg(:,1:ndim), Finv=Finv, esize=esize, &
          Uglobal=Uglobal, checkzero=checkzero )

        beta = coefficients%r(5)

        tau = beta * hoverU / 2  ! upwinding parameter

      end if

    end if

    deltat = coefficients%r(4)

    if ( matrix ) then

      if ( method == 0 ) then

!       Galerkin
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( phi(:,i) * &
              ( 1.5_dp * phi(:,j) / deltat + ugradphi(:,j) ) * detF * wg )
          end do
        end do

      else if ( method == 1 ) then

!       SUPG
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( ( phi(:,i) + tau * ugradphi(:,i)) * &
              ( 1.5_dp * phi(:,j) / deltat + ugradphi(:,j) ) * detF * wg )
          end do
        end do

      end if

    end if

    if ( vector ) then

      if ( method == 0 ) then

!       Galerkin
        do i = 1, ndf
          elemvec(i) = sum ( phi(:,i) * &
             ( ( 2 * Hng - 0.5_dp * Hnm1g ) / deltat + uvecg(:,ndim+1) ) &
                          * detF * wg )
        end do

      else if ( method == 1 ) then

!       SUPG
        do i = 1, ndf
          elemvec(i) = sum ( ( phi(:,i) +  tau * ugradphi(:,i) ) * &
             ( ( 2 * Hng - 0.5_dp * Hnm1g ) / deltat + uvecg(:,ndim+1) ) &
                          * detF * wg )
        end do

      end if

    end if

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim) )
      allocate ( x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( F(ninti,ndim,ndim), Finv(ninti,ndim,ndim) )
      allocate ( fg(ninti), detF(ninti) )
      allocate ( dphidx(ninti,ndf,ndim) )
      allocate ( Hn(ndf), Hng(ninti) )
      allocate ( Hnm1(ndf), Hnm1g(ninti) )
      allocate ( uvecg(ninti,ncompv) )
      allocate ( tau(ninti), hoverU(ninti) )
      allocate ( ugradphi(ninti,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x, xg )
      deallocate ( F, Finv )
      deallocate ( fg, detF )
      deallocate ( dphidx )
      deallocate ( Hn, Hng )
      deallocate ( Hnm1, Hnm1g )
      deallocate ( uvecg )
      deallocate ( tau, hoverU )
      deallocate ( ugradphi )

    end subroutine deallocate_arrays

  end subroutine surface_advection_polar_elem2


! Compute height function in all nodes

  subroutine height_function_deriv (  mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    use surface_advection_globals_m

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

      call set_globals_surface_advection ( mesh, coefficients, elgrp )

      if ( ndim > 2 ) then
        write(*,'(/2(a/))') 'Error in height_function_deriv:', &
          ' Space dimension of mesh ndim must be 1 or 2'
        stop
      end if

      ncompv = ndim + 1

!     allocate arrays

      call allocate_arrays

!     shape function in the nodes

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefunc, xrnod, phi )

    end if

!   get height function

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u=Hf )

    elemvec = matmul ( phi, Hf )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xrnod(nodalp,ndim) )
      allocate ( phi(nodalp,ndf) )
      allocate ( Hf(ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xrnod )
      deallocate ( phi )
      deallocate ( Hf )

    end subroutine deallocate_arrays

  end subroutine height_function_deriv


! set global parameters

  subroutine set_globals_surface_advection ( mesh, coefficients, elgrp )

    use surface_advection_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

    integer :: ndimr

!   check size of coefficients

    call check ( coefficients, 'set_globals_surface_advection', ncoefi=100, &
      ncoefr=50, indexarray=[12,14], minimum=[0,0], maximum=[1,1] )

    ndim = mesh%element(elgrp)%ndim
    nodalp = mesh%element(elgrp)%numnod
    globalshape = mesh%element(elgrp)%globalshape
    isoshape = coefficients%i(12)

!   set number of degrees of freedom for height function

    intpol = coefficients%i(6)

    ndimr = mesh%element(elgrp)%ndimr
    if ( intpol == 13 .and. any(mesh%element(elgrp)%p(:ndimr,2) /= 1 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_surface_advection:', &
        ' For spectral elements GLL nodal distribution is required. '
      stop
    else if ( intpol == 20 .and. &
                            any(mesh%element(elgrp)%p(:ndimr,2) /= 0 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_surface_advection:', &
        ' For high-order elements equidistant nodal distribution is required. '
      stop
    end if

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(7)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, 'set_globals_surface_advection', ndf=ndf )

!   set number of degrees of freedom for velocity

    intpolv = coefficients%i(10)

    if ( intpolv > 0 ) then

      shapefuncv%globalshape = globalshape
      shapefuncv%interpolation = intpolv
      shapefuncv%numbering = 'standard'
      shapefuncv%p = coefficients%i(11)
      shapefuncv%spec_eval = 'gauss'

      call set_ndf ( shapefuncv, 'set_globals_surface_advection', ndf=ndfv )

    else

     shapefuncv = shapefunc
     ndfv = ndf

    end if

!   set integration

    inttype = coefficients%i(9)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_surface_advection:', &
        ' For spectral elements Gauss-Legendre-Lobatto integration ', &
        ' needs to be specified. '
      stop
    end if

    if ( coefficients%i(13) == 1 .or. SET_GAUSS_BY_ORDER ) then
      intrule = set_intrule ( globalshape, inttype, order=coefficients%i(1) )
    else
      intrule = coefficients%i(1)
    end if
    nsubint = get_coefficient ( coefficients, index=8, default=1 )

    gauss%globalshape = globalshape
    gauss%intrule = intrule
    gauss%nsubint = nsubint
    gauss%inttype = inttype

    call set_ninti ( gauss, ninti )

    if ( intpol == 13 .and. ninti /= ndf ) then
      write(*,'(/5(a/))') 'Error in set_globals_surface_advection:', &
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_surface_advection

end module surface_advection_elements_m
