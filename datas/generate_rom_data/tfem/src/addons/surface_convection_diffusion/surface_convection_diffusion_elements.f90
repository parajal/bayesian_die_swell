
! Copyright (C) 2017-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Element routines to calculate concentration on a curve/surface according to
!
!    dc                __       __           __       __
!    -- +  ( u - u ) . \/ c + ( \/ . u ) c - \/ . ( D \/ c ) = 0
!    dt      -   -m      s        s  -         s        s
!
! where
!
!  d/dt the grid time derivative
!
!  c the concentration
!
!  u the velocity vector field
!  -
!  u the grid velocity field
!  -m
!
!  D the diffusion coefficient


module surface_convection_diffusion_elements_m

  use tfem_elem_m

  implicit none

contains

  subroutine surface_convection_diffusion_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use surface_convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, k, m, ip, timeint
    real(dp) :: Dcoef, deltat


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_surface_convection_diffusion ( mesh, coefficients, &
        elgrp )

!     check

      call check ( coefficients, 'surface_convection_diffusion_elem', &
        indexarray=[6], minimum=[1], maximum=[2] )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

!   get coordinates of nodes of interface

    call get_coordinates ( mesh, elgrp, elem, x )

    xg = matmul ( phi, x )

!   surface geometry

    call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, &
      gij_up=gij_up, gi_up=gi_up )

    if ( coorsys == 1 ) then
      surfl = 2 * pi * xg(:,2) * surfl
    end if

!   time integration scheme

    timeint = coefficients%i(6)

!   concentration at tn

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, &
      u=car )

    cg = matmul( phi, car )

    if ( timeint == 2 ) then

      call get_sysvector ( mesh, problem, oldvectors%s(2)%p, elgrp, elem, &
        u=car )

      cm1g = matmul( phi, car )

    end if

!   velocity at tn

    call get_vector ( mesh, oldvectors%p(2)%p, oldvectors%v(1)%p, elgrp, &
      elem, u=uar )

    un = reshape ( uar, [ndf,ndim] )
    ug = matmul ( phi, un )

!   get grid velocities u
!                       -m

    call get_vector ( mesh, oldvectors%p(2)%p, oldvectors%v(2)%p, elgrp, &
      elem, u=uar )

    uvecmg = matmul ( phi, reshape ( uar, [ndf,ndim] ) )

!   diffusion coefficient

    Dcoef = coefficients%r(4)

!   set time step

    deltat = coefficients%r(5)


    if ( matrix ) then

      elemmat = 0

      do i = 1, ndim
        do m = 1, ndf
          do ip = 1, ninti
            grad(ip,m,i) = sum ( dphi(ip,m,:) * gi_up(ip,i,:) )
          end do
        end do
      end do
      do m = 1, ndf
        do ip = 1, ninti
          ugrad(ip,m) = sum (grad(ip,m,:) * ( ug(ip,:) - uvecmg(ip,:) ) )
        end do
      end do

      do i = 1, ndim-1
        dphiu(:,:,i) = matmul ( dphi(:,:,i), un )
      end do
      do i = 1, ndim
        do ip = 1, ninti
          gdphiu(ip,i) = sum ( dphiu(ip,i,:) * gi_up(ip,i,:) )
        end do
      end do


!     dc/dt term

      if ( timeint == 1 ) then

        do m = 1, ndf
          do k = 1, ndf
            do ip = 1, ninti
              work(ip) =  phi(ip,k) * phi(ip,m)
            end do
            elemmat(k,m) = sum ( work * surfl * wg ) / deltat
          end do
        end do

      else if ( timeint == 2 ) then

        do m = 1, ndf
          do k = 1, ndf
            do ip = 1, ninti
              work(ip) =  1.5_dp * phi(ip,k) * phi(ip,m)
            end do
            elemmat(k,m) = sum ( work * surfl * wg ) / deltat
          end do
        end do

      end if

!     change in area
!       __
!     ( \/ . u ) c
!        s   -

      do m = 1, ndf
        do k = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( gdphiu(ip,:)) * phi(ip,k) * phi(ip,m)
          end do
          if ( coorsys == 1 ) then
            work = work + ug(:,2) * phi(:,k) * phi(:,m) / xg(:,2)
          end if
          elemmat(k,m) = elemmat(k,m) + sum ( work * surfl * wg )
        end do
      end do

!     convection term
!                __
!     ( u - u ). \/ c
!       -   -m     s

      do m = 1, ndf
        do k = 1, ndf
          elemmat(k,m) = elemmat(k,m) + &
                                  sum (ugrad(:,m) * phi(:,k) * surfl * wg )
        end do
      end do

!     diffusion term
!       __      __
!     - \/ . (D \/ c )
!        s       s
!     where D is the diffusion coefficient

      do m = 1, ndf
        do k = 1, ndf
          do ip = 1, ninti
            work(ip) = &
                  sum ( dphi(ip,m,:) * matmul( gij_up(ip,:,:), dphi(ip,k,:) ) )
          end do
          elemmat(k,m) = elemmat(k,m) + Dcoef * sum ( work * surfl * wg )
        end do
      end do

    end if

    if ( vector ) then

      elemvec = 0

      if ( timeint == 1 ) then

        do k = 1, ndf
          elemvec(k) = sum ( phi(:,k) * cg * surfl * wg ) / deltat
        end do

      else if ( timeint == 2 ) then

        do k = 1, ndf
          elemvec(k) = sum ( phi(:,k) * &
                         (2.0_dp * cg - 0.5_dp * cm1g)  * surfl * wg ) / deltat
        end do

      end if

    end if

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( un(ndf,ndim) )
      allocate ( car(ndf), cg(ninti), cm1g(ninti) )
      allocate ( uar(ndf*ndim) )
      allocate ( surfl(ninti), work(ninti) )
      allocate ( gi_up(ninti,ndim,ndim-1) )
      allocate ( gij_up(ninti,ndim-1,ndim-1) )
      allocate ( ug (ninti, ndim ), uvecmg(ninti, ndim) )
      allocate ( grad ( ninti, ndf, ndim ) )
      allocate ( ugrad ( ninti, ndf ) )
      allocate ( dphiu ( ninti, ndim ,ndim-1 ) )
      allocate ( gdphiu ( ninti, ndim  ) )
      allocate ( xg ( ninti, ndim ) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( wg )
      deallocate ( xig )
      deallocate ( phi, dphi )
      deallocate ( x )
      deallocate ( dxdxis )
      deallocate ( un )
      deallocate ( car, cg, cm1g )
      deallocate ( uar )
      deallocate ( surfl, work )
      deallocate ( gi_up )
      deallocate ( gij_up )
      deallocate ( ug, uvecmg )
      deallocate ( grad )
      deallocate ( ugrad )
      deallocate ( dphiu, gdphiu)
      deallocate ( xg )

    end subroutine deallocate_arrays

  end subroutine surface_convection_diffusion_elem


! integrate concentration

  subroutine integrate_surface_convection_diffusion ( mesh, problem, &
    elgrp, elem, first, last, coefficients, oldvectors, elemvec )

    use surface_convection_diffusion_globals_m

!   input/output
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

      call set_globals_surface_convection_diffusion ( mesh, coefficients, &
        elgrp )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

!   get coordinates of nodes of interface

    call get_coordinates ( mesh, elgrp, elem, x )

!   surface geometry

    call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl )

    if ( coorsys == 1 ) then
      xg = matmul ( phi, x )
      surfl = 2 * pi * xg(:,2) * surfl
    end if

!   get concentration

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, &
      u=car )

    work = matmul ( phi, car )

    elemvec = sum ( work * wg * surfl )

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( work(ninti) )
      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( xg ( ninti, ndim ) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )
      allocate ( car(ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work )
      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( x )
      deallocate ( xg )
      deallocate ( dxdxis )
      deallocate ( surfl )
      deallocate ( car )

    end subroutine deallocate_arrays

  end subroutine integrate_surface_convection_diffusion


! area of the surface

  subroutine integrate_surface_area_elem ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec )

    use surface_convection_diffusion_globals_m

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

      call set_globals_surface_convection_diffusion ( mesh, coefficients, &
        elgrp )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

!   get coordinates of nodes of interface

    call get_coordinates ( mesh, elgrp, elem, x )

!   surface geometry

    call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl )

    if ( coorsys == 1 ) then
      xg = matmul ( phi, x )
      surfl = 2 * pi * xg(:,2) * surfl
    end if

!   area

    elemvec = sum ( surfl * wg )

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( xg ( ninti, ndim ) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( x )
      deallocate ( xg )
      deallocate ( dxdxis )
      deallocate ( surfl )

    end subroutine deallocate_arrays

  end subroutine integrate_surface_area_elem


! set global parameters

  subroutine set_globals_surface_convection_diffusion ( mesh, coefficients, &
    elgrp )

    use surface_convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

    integer :: ndimr

!   check size of coefficients

    call check ( coefficients, 'set_globals_surface_convection_diffusion', &
      ncoefi=100, ncoefr=50, indexarray=[7], minimum=[0], maximum=[1]  )

    ndim = mesh%element(elgrp)%ndim
    nodalp = mesh%element(elgrp)%numnod
    globalshape = mesh%element(elgrp)%globalshape
    if ( ndim == 2 ) then
      coorsys = coefficients%i(7)
    else
      coorsys = 2
    end if

!   set number of degrees of freedom velocity

    intpol = coefficients%i(1)

    ndimr = mesh%element(elgrp)%ndimr
    if ( intpol == 13 .and. any(mesh%element(elgrp)%p(:ndimr,2) /= 1 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_surface_convection_diffusion:',&
        ' For spectral elements GLL nodal distribution is required. '
      stop
    else if ( intpol == 20 .and. &
                            any(mesh%element(elgrp)%p(:ndimr,2) /= 0 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_surface_convection_diffusion:',&
        ' For high-order elements equidistant nodal distribution is required. '
      stop
    end if

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(2)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, 'set_globals_surface_convection_diffusion', &
      ndf=ndf )

!   set integration

    inttype = coefficients%i(5)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_surface_convection_diffusion:',&
        ' For spectral elements Gauss-Legendre-Lobatto integration ', &
        ' needs to be specified. '
    end if

    intrule = coefficients%i(3)
    nsubint = get_coefficient ( coefficients, index=4, default=1 )

    gauss%globalshape = globalshape
    gauss%intrule = intrule
    gauss%nsubint = nsubint
    gauss%inttype = inttype

    call set_ninti ( gauss, ninti )

    if ( intpol == 13 .and. ninti /= ndf ) then
      write(*,'(/5(a/))') 'Error in set_globals_surface_convection_diffusion:',&
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_surface_convection_diffusion

end module surface_convection_diffusion_elements_m
