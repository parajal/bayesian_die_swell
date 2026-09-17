
! Copyright (C) 2011-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for embedded boundary or interface conditions for the
! viscoelastic flow equation using an eltree to divide the domain and/or
! introduce a boundary or interface.

module viscoelastic_elements_embedded_m

  use tfem_elem_m
  use stokes_set_globals_m
  use stokes_elements_embedded_boundary_m, only: set_stokes_eltree, &
      unset_stokes_eltree
  use viscoelastic_elements_generic_m

  implicit none

contains


! Element routine on an eltree for open boundary conditions in the implicit
! stress formulation

  subroutine implicit_stress_open_boundary_eltree ( mesh, problem, elgrp, &
    elem, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: ip, orientation
    integer :: i, j

    integer :: i1, i2, i3, N, M
    integer :: md
    real(dp) :: deltat


!   check coefficients

    call check ( coefficients, 'implicit_stress_open_boundary_eltree', &
      indexarray=[36,48], minimum=[-1,0], maximum=[1,1] )

    cstorage = get_coefficient ( coefficients, index=84, default=1 )

!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

    call set_viscoelastic_model ( coefficients )

    call set_globals_viscoelastic ( coefficients )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( theta(ninti,ndfc) )
    allocate ( dtheta(ninti,ndfc,ndim), dthetadx(ninti,ndfc,ndim) )
    allocate ( c(ndfc,ncomp), cn(ndfc,ncompc), cng(ninti,ncompc,nmodes) )
    allocate ( gradcn(ninti,ncompc,ndim,nmodes) )
    allocate ( cten(ninti,ndim,ndim), gradcten(ninti,ndim,ndim,ndim) )
    allocate ( hten(ninti,ndim,ndim), rhsd(ninti,ncompc,nmodes) )
    allocate ( Gmod(mode2-mode1+1) )
    allocate ( work2(ndf,ndim), work5(ndf,ndf,ndim,ndim) )
    allocate ( work10(ndf,ndf), work11(ninti,ndim) )
    allocate ( work(ninti), work1(ninti), work8(ninti) )
    allocate ( work4(ninti,ndim), tmp(ndf,ndim) )

    allocate ( u(ndf*ndim), uvecmeshnp1(ninti,ndim) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    call set_shape_function ( shapefuncc, xig, theta, dtheta ) ! c tensor

    call shape_derivative ( dtheta, Finv, dthetadx )

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
    end do

!   set orientation

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    dan = real ( orientation, kind=dp ) * dan


!   get conformation tensor

    do md = mode1, mode2

      if ( exps .or. cproj ) then

!       projection of c=exp(s) or c=b.b^T has been performed in
!       a separate problem

!       single mode conformation
        call get_conformation ( mesh, oldvectors%p(3)%p, elgrp, elem, &
          oldvectors, cmode=cn, mode=md, isv=3, cst=1 )

      else

!       standard case

!       single mode conformation
        call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
          oldvectors, cmode=c, mode=md )

        if ( logc ) then
!         transform sn to cn
!         NOTE: we are transforming the nodal values sn -> cn directly and
!         define cn in the element with the same interpolation as sn.
!         It is better to do a global least-square projection of cn = exp(sn)
!         to define the nodal point values of cn. See exps=.true.
          if ( coorsys <= 1 ) then
            call conformation_2D_log ( vemodel, c, cn )
          else if ( coorsys == 2 ) then
            call conformation_3D_log ( vemodel, c, cn )
          end if
        else if ( bten ) then
!         transform bn to cn
!         NOTE: we are transforming the nodal values bn -> cn directly and
!         define cn in the element with the same interpolation as bn. Maybe it
!         is better to do a global least-square projection of cn = bn.bn^T to
!         define the nodal point values of cn. See cproj=.true.
          if ( coorsys <= 1 ) then
            call conformation_2D_b ( c, cn, vemodel%bvariant )
          else if ( coorsys == 2 ) then
            call conformation_3D_b ( c, cn, vemodel%bvariant )
          end if
        else
          cn = c
        end if

      end if

!     c in Gauss points

      cng(:,:,md) = matmul ( theta, cn )

!     grad cn term

      do i = 1, ndim
        gradcn(:,:,i,md) = matmul ( dthetadx(:,:,i), cn )
      end do

    end do

!   get mesh velocity at the current time

    if ( coefficients%i(48) == 1 ) then

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, u )

      tmp = reshape ( u, [ndf,ndim] )

      uvecmeshnp1 = matmul ( phi, tmp )

    end if

!   some parameters

    deltat = coefficients%r(8)

    Gmod = vemodel%modulus(mode1:mode2)


    if ( matrix ) then

!     set c in tensor format

      if ( coorsys <= 1 ) then

        cten(:,1,1) = matmul ( cng(:,1,mode1:mode2), Gmod )
        cten(:,1,2) = matmul ( cng(:,2,mode1:mode2), Gmod )
        cten(:,2,1) = cten(:,1,2)
        cten(:,2,2) = matmul ( cng(:,3,mode1:mode2), Gmod )

      else if ( coorsys == 2 ) then

        cten(:,1,1) = matmul ( cng(:,1,mode1:mode2), Gmod )
        cten(:,1,2) = matmul ( cng(:,2,mode1:mode2), Gmod )
        cten(:,1,3) = matmul ( cng(:,3,mode1:mode2), Gmod )
        cten(:,2,1) = cten(:,1,2)
        cten(:,2,2) = matmul ( cng(:,4,mode1:mode2), Gmod )
        cten(:,2,3) = matmul ( cng(:,5,mode1:mode2), Gmod )
        cten(:,3,1) = cten(:,1,3)
        cten(:,3,2) = cten(:,2,3)
        cten(:,3,3) = matmul ( cng(:,6,mode1:mode2), Gmod )

      end if

!     set gradc in tensor format

      if ( coorsys <= 1 ) then

        do i = 1, ndim
          gradcten(:,1,1,i) = matmul ( gradcn(:,1,i,mode1:mode2), Gmod )
          gradcten(:,1,2,i) = matmul ( gradcn(:,2,i,mode1:mode2), Gmod )
          gradcten(:,2,1,i) = gradcten(:,1,2,i)
          gradcten(:,2,2,i) = matmul ( gradcn(:,3,i,mode1:mode2), Gmod )
        end do

      else if ( coorsys == 2 ) then

        do i = 1, ndim
          gradcten(:,1,1,i) = matmul ( gradcn(:,1,i,mode1:mode2), Gmod )
          gradcten(:,1,2,i) = matmul ( gradcn(:,2,i,mode1:mode2), Gmod )
          gradcten(:,1,3,i) = matmul ( gradcn(:,3,i,mode1:mode2), Gmod )
          gradcten(:,2,1,i) = gradcten(:,1,2,i)
          gradcten(:,2,2,i) = matmul ( gradcn(:,4,i,mode1:mode2), Gmod )
          gradcten(:,2,3,i) = matmul ( gradcn(:,5,i,mode1:mode2), Gmod )
          gradcten(:,3,1,i) = gradcten(:,1,3,i)
          gradcten(:,3,2,i) = gradcten(:,2,3,i)
          gradcten(:,3,3,i) = matmul ( gradcn(:,6,i,mode1:mode2), Gmod )
        end do

      end if


      do N = 1, ndf
        do M = 1, ndf

          do i = 1, ndim
            do j = 1, ndim

              do ip = 1, ninti
                ! -( v, (u dot gradc) dot n )
                work(ip) = - dot_product ( gradcten(ip,i,:,j), dan(ip,:) )

                ! ( v, (c dot gradu) dot n )
                work1(ip) = dot_product ( cten(ip,i,:), dphidx(ip,M,:) )

              end do

              work5(N,M,i,j) = sum ( phi(:,N) * phi(:,M) * work(:) ) &
                + sum ( phi(:,N) * work1(:) * dan(:,j) )

            end do
          end do


          ! ( v, (gradu^T dot c) dot n )
          do ip = 1, ninti
            work11(ip,:) = matmul ( cten(ip,:,:), dan(ip,:) )
            work8(ip) = dot_product ( work11(ip,:), dphidx(ip,M,:) )
          end do

          work10(N,M) = sum ( phi(:,N) * work8 )

        end do
      end do

      work5 = deltat * work5
      work10 = deltat * work10

      ! element matrix

      i1 = ndf    ! u
      i2 = 2*ndf  ! v
      i3 = 3*ndf  ! w

      if ( coorsys <= 1 ) then

        elemmat( 1:i1, 1:i1 )       = work5(:,:,1,1) + work10(:,:)
        elemmat( 1:i1, i1+1:i2 )    = work5(:,:,1,2)
        elemmat( i1+1:i2, 1:i1 )    = work5(:,:,2,1)
        elemmat( i1+1:i2, i1+1:i2 ) = work5(:,:,2,2) + work10(:,:)

      else if ( coorsys == 2 ) then

        elemmat( 1:i1, 1:i1 )       = work5(:,:,1,1) + work10(:,:)
        elemmat( 1:i1, i1+1:i2 )    = work5(:,:,1,2)
        elemmat( 1:i1, i2+1:i3 )    = work5(:,:,1,3)
        elemmat( i1+1:i2, 1:i1 )    = work5(:,:,2,1)
        elemmat( i1+1:i2, i1+1:i2 ) = work5(:,:,2,2) + work10(:,:)
        elemmat( i1+1:i2, i2+1:i3 ) = work5(:,:,2,3)
        elemmat( i2+1:i3, 1:i1 )    = work5(:,:,3,1)
        elemmat( i2+1:i3, i1+1:i2 ) = work5(:,:,3,2)
        elemmat( i2+1:i3, i2+1:i3 ) = work5(:,:,3,3) + work10(:,:)

      end if

      elemmat = -elemmat

    end if


    if ( vector ) then

!     - (nabla v)^T: ( G h )

!     relaxation term of CE only

      if ( coorsys <= 1 ) then
        call rhs_viscoelastic_relax_2D ( vemodel, cng, rhsd, mode1=mode1, &
          mode2=mode2)
      else if ( coorsys == 2 ) then
        call rhs_viscoelastic_relax_3D ( vemodel, cng, rhsd, mode1=mode1, &
          mode2=mode2)
      end if

!     add mesh velocity term to the right-hand side

      if ( coefficients%i(48) == 1 ) then

        do md = 1, nmodes
          do ip = 1, ninti
            rhsd(ip,:,md) = rhsd(ip,:,md) &
                          + matmul ( gradcn(ip,:,:,md), uvecmeshnp1(ip,:) )
          end do
        end do

      end if

      if ( coorsys <= 1 ) then

        hten(:,1,1) = matmul ( cng(:,1,mode1:mode2) + deltat * ( &
            rhsd(:,1,mode1:mode2) ) - 1, Gmod )
        hten(:,1,2) = matmul ( cng(:,2,mode1:mode2) + deltat * ( &
            rhsd(:,2,mode1:mode2) ), Gmod )
        hten(:,2,1) = hten(:,1,2)
        hten(:,2,2) = matmul ( cng(:,3,mode1:mode2) + deltat * ( &
            rhsd(:,3,mode1:mode2) ) - 1, Gmod )

      else if ( coorsys == 2 ) then

        hten(:,1,1) = matmul ( cng(:,1,mode1:mode2) + deltat * ( &
            rhsd(:,1,mode1:mode2) ) - 1, Gmod )
        hten(:,1,2) = matmul ( cng(:,2,mode1:mode2) + deltat * ( &
            rhsd(:,2,mode1:mode2) ), Gmod )
        hten(:,1,3) = matmul ( cng(:,3,mode1:mode2) + deltat * ( &
            rhsd(:,3,mode1:mode2) ), Gmod )
        hten(:,2,1) = hten(:,1,2)
        hten(:,2,2) = matmul ( cng(:,4,mode1:mode2) + deltat * ( &
            rhsd(:,4,mode1:mode2) ) - 1, Gmod )
        hten(:,2,3) = matmul ( cng(:,5,mode1:mode2) + deltat * ( &
            rhsd(:,5,mode1:mode2) ), Gmod )
        hten(:,3,1) = hten(:,1,3)
        hten(:,3,2) = hten(:,2,3)
        hten(:,3,3) = matmul ( cng(:,6,mode1:mode2) + deltat * ( &
            rhsd(:,6,mode1:mode2) ) - 1, Gmod )

      end if

!     -(v, hten.n)

      do ip = 1, ninti
        work4(ip,:) = matmul( hten(ip,:,:), dan(ip,:) )
      end do

      do i = 1, ndf
        do j = 1, ndim
          work2(i,j) = sum ( phi(:,i) * work4(:,j) )
        end do
      end do
      elemvec = reshape ( work2, [ ndf*ndim ] )

    end if

!   deallocate arrays

    call delete ( vemodel )

    deallocate ( theta, dtheta, dthetadx )
    deallocate ( c, cn, cng )
    deallocate ( gradcn )
    deallocate ( cten, gradcten )
    deallocate ( hten, rhsd )
    deallocate ( Gmod )
    deallocate ( work2, work5 )
    deallocate ( work10, work11 )
    deallocate ( work, work1, work8 )
    deallocate ( work4, tmp )

    deallocate ( u, uvecmeshnp1 )

    call unset_stokes_eltree

  end subroutine implicit_stress_open_boundary_eltree


! viscoelastic stress contribution to the drag

  subroutine integration_viscoelastic_drag_eltree ( mesh, problem, elgrp, &
     elem, first, last, coefficients, oldvectors, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: ip, i, m, orientation

!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   set globals viscoelastic

    call set_viscoelastic_model ( coefficients )

    call set_globals_viscoelastic ( coefficients )

    cstorage = get_coefficient ( coefficients, index=84, default=1 )

!   allocate arrays

    allocate ( theta(ninti,ndfc) )
    allocate ( c(ndfc,ncomp), cg(ninti,ncomp,nmodes) )
    allocate ( tauten(ninti,ndim,ndim), tauvec(ninti,ncompt) )
    allocate ( work4(ninti,ndim) )

!   set shape functions

    call set_shape_function ( shapefunc, xig, phi, dphi ) ! velocity
    call set_shape_function ( shapefuncc, xig, theta ) ! c tensor

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
    end do

!   set orientation

    call check ( coefficients, 'integration_viscous_drag_eltree', &
      indexarray=[36], minimum=[-1], maximum=[1] )

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    dan = real ( orientation, kind=dp ) * dan


!   get conformation tensor

    do m = mode1, mode2

!     single mode conformation
      call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
        cmode=c, mode=m )

      cg(:,:,m) = matmul ( theta, c )

    end do

    call stress_tensor_viscoelastic ( cg )

!   evaluate traction force in each integration point: {t}=[tau_p]{n}

    do ip = 1, ninti
      work4(ip,:) = matmul( tauten(ip,:,:), dan(ip,:) )
    end do

!   integrate

    do i = 1, ndim
      elemvec(i) = sum ( work4(:,i) )
    end do

!   deallocate arrays

    call delete ( vemodel )

    deallocate ( theta )
    deallocate ( c, cg )
    deallocate ( tauten, tauvec )
    deallocate ( work4 )

    call unset_stokes_eltree

  end subroutine integration_viscoelastic_drag_eltree


end module viscoelastic_elements_embedded_m
