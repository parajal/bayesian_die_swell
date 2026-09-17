
! Copyright (C) 2020-2020 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Element routines for the linear elastic equation (Navier)
!
!    - div( lambda div u I + mu (nabla u+nabla u^T) ) = f
!
! in displacement formulation or
!
!    - div( mu (nabla u+nabla u^T) ) + nabla p = f
!      p/lambda + div u = 0
!
! in mixed displacement-pressure formulation.
!
! If thermal expansion is included the equations become
!
!    - div( lambda div u I + mu (nabla u+nabla u^T) + 3*K*alpha*DT I ) = f
!
! in displacement formulation or
!
!    - div( mu (nabla u+nabla u^T) ) + nabla p = f
!      p/lambda + div u - alpha*DT*(1+nu)/nu  = 0
!
! in mixed displacement-pressure formulation.
!
! Here (lambda,mu) are the Lame parameters and K is the bulk modulus
!
!               nu * E                   E              E
!  lambda = --------------,  mu = G = -------,  K = ---------
!           (1+nu)(1-2*nu)            2(1+nu)       3(1-2*nu)
!
! where E is Young's modulus and nu Poisson's ratio. Note, G is the shear
! modulus.
!
! This module contains stuff that only works in 3D or only for surfaces.
!

module linear_elastic_elements_3D_surface_m

  use tfem_elem_m
  use linear_elastic_set_globals_m
  use linear_elastic_elements_generic_m, &
                  only: set_linear_elastic_shape_function_global
  use linear_elastic_material_m

  implicit none

contains


! Internal element routine for the linear elastic equation (3D)

  subroutine linear_elastic_elem_3D ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use linear_elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, M, ip, i1, i2, i3, i4, vfuncnr, j
    real(dp) :: mu, lambda, facKmod, alpha, fac


!   set globals, gauss, shapefunctions, ...

    call set_linear_elastic_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

    call set_linear_elastic_material ( coefficients, coorsys, lemodel )

!   mixed elements needed for incompressible materials
    if ( lemodel%incompressible .and. ndfp == 0 ) &
      call errormsg_notmixed_incompressible ( 'linear_elastic_elem_3D' )


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

    if ( ndfp > 0 .and. any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_linear_elastic_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   set material parameters

    if ( .not. lemodel%incompressible ) lambda = lemodel%lambda
    mu = lemodel%mu
    if ( lemodel%thermal_expansion) then
      alpha = lemodel%alpha
      if ( .not. lemodel%incompressible ) facKmod = 3 * lemodel%Kmod
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

      if ( lemodel%thermal_expansion ) then

!       thermal expansion: get temperature change

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(9), value=coefficients%r(4), &
          func=coefficients%func1(1)%p, funcnr=coefficients%i(8), x=xg, &
          indx_v=1, layer=layer, phi=phi, indx_e=1, coef=dTemp )

        if ( ndfp > 0 ) then

!         mixed (u,p) elements

          fac = alpha * (1+1/lemodel%nu)
          do N = 1, ndfp
            fp(N) = - fac * sum ( psi(:,N) * dTemp * detF * wg )
          end do

          elemvec(i3+1:i4) = fp

        else

!         standard displacement elements

          do j = 1, ndim
            do N = 1, ndf
              work = dphidx(:,N,j)
              tmp(N,j) = facKmod * alpha * sum ( work * dTemp * detF * wg )
            end do
          end do

          elemvec(1:i3) = elemvec(1:i3) + reshape ( tmp, [ndim*ndf] )

        end if

      end if

    end if

    if ( matrix ) then

      do N = 1, ndf
        do M = N, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dphidx(ip,M,:) )
          end do
          work1 = work + dphidx(:,N,1) * dphidx(:,M,1)
          Suu(N,M) = mu * sum ( work1 * detF * wg )
          Suu(M,N) = Suu(N,M)
          work1 = work + dphidx(:,N,2) * dphidx(:,M,2)
          Svv(N,M) = mu * sum ( work1 * detF * wg )
          Svv(M,N) = Svv(N,M)
          work1 = work + dphidx(:,N,3) * dphidx(:,M,3)
          Sww(N,M) = mu * sum ( work1 * detF * wg )
          Sww(M,N) = Sww(N,M)
        end do
      end do

      do N = 1, ndf
        do M = 1, ndf
          work = dphidx(:,N,2) * dphidx(:,M,1)
          Suv(N,M) = mu * sum ( work * detF * wg )
          work = dphidx(:,N,3) * dphidx(:,M,1)
          Suw(N,M) = mu * sum ( work * detF * wg )
          work = dphidx(:,N,3) * dphidx(:,M,2)
          Svw(N,M) = mu * sum ( work * detF * wg )
        end do
      end do

      if ( ndfp > 0 ) then

!       mixed (u,p) elements

!       pressure-pressure part

        if ( lemodel%incompressible ) then
          Cmat = 0
        else
          do N = 1, ndfp
            do M = N, ndfp
              Cmat(N,M) = sum ( psi(:,N) * psi(:,M) * detF * wg )
              Cmat(M,N) = Cmat(N,M)
            end do
          end do
          Cmat = Cmat / lambda
        end if

!       displacement-pressure part

        do N = 1, ndfp
          do M = 1, ndf
            Lu(N,M) = sum ( psi(:,N) * dphidx(:,M,1) * detF * wg )
            Lv(N,M) = sum ( psi(:,N) * dphidx(:,M,2) * detF * wg )
            Lw(N,M) = sum ( psi(:,N) * dphidx(:,M,3) * detF * wg )
          end do
        end do

      else

!       standard displacement elements

        do N = 1, ndf
          do M = N, ndf
            work = dphidx(:,N,1) * dphidx(:,M,1)
            Suu(N,M) = Suu(N,M) + lambda * sum ( work * detF * wg )
            Suu(M,N) = Suu(N,M)
            work = dphidx(:,N,2) * dphidx(:,M,2)
            Svv(N,M) = Svv(N,M) + lambda * sum ( work * detF * wg )
            Svv(M,N) = Svv(N,M)
            work = dphidx(:,N,3) * dphidx(:,M,3)
            Sww(N,M) = Sww(N,M) + lambda * sum ( work * detF * wg )
            Sww(M,N) = Sww(N,M)
          end do
        end do

        do N = 1, ndf
          do M = 1, ndf
            work = dphidx(:,N,1) * dphidx(:,M,2)
            Suv(N,M) = Suv(N,M) + lambda * sum ( work * detF * wg )
            work = dphidx(:,N,1) * dphidx(:,M,3)
            Suw(N,M) = Suw(N,M) + lambda * sum ( work * detF * wg )
            work = dphidx(:,N,2) * dphidx(:,M,3)
            Svw(N,M) = Svw(N,M) + lambda * sum ( work * detF * wg )
          end do
        end do

      end if

      elemmat(    1:i1,    1:i1 ) = Suu
      elemmat(    1:i1, i1+1:i2 ) = Suv
      elemmat(    1:i1, i2+1:i3 ) = Suw
      elemmat( i1+1:i2,    1:i1 ) = transpose(Suv)
      elemmat( i1+1:i2, i1+1:i2 ) = Svv
      elemmat( i1+1:i2, i2+1:i3 ) = Svw
      elemmat( i2+1:i3,    1:i1 ) = transpose(Suw)
      elemmat( i2+1:i3, i1+1:i2 ) = transpose(Svw)
      elemmat( i2+1:i3, i2+1:i3 ) = Sww

      if ( ndfp > 0 ) then
        elemmat( i3+1:i4,    1:i1 ) = -Lu
        elemmat( i3+1:i4, i1+1:i2 ) = -Lv
        elemmat( i3+1:i4, i2+1:i3 ) = -Lw
        elemmat(    1:i1, i3+1:i4 ) = -transpose(Lu)
        elemmat( i1+1:i2, i3+1:i4 ) = -transpose(Lv)
        elemmat( i2+1:i3, i3+1:i4 ) = -transpose(Lw)
        elemmat( i3+1:i4, i3+1:i4)  = -Cmat
      end if

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_linear_elastic_elem ( last, coefficients )

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
      allocate ( Cmat(ndfp,ndfp) )
      allocate ( dTemp(ninti), fp(ndfp) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( tmp )
      deallocate ( work, work1 )
      deallocate ( Suu, Suv, Svv )
      deallocate ( Suw, Svw, Sww )
      deallocate ( Lu, Lv, Lw )
      deallocate ( Cmat )
      deallocate ( dTemp, fp )

    end subroutine deallocate_arrays

  end subroutine linear_elastic_elem_3D


! Boundary element for a natural boundary on a surface for the linear elastic
! equation

  subroutine linear_elastic_natboun_surface ( mesh, problem, surface, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use linear_elastic_globals_m

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

      call set_globals_linear_elastic_up_boun ( mesh, coefficients, &
        surface=surface )

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
      write(*,'(/2(a/))') 'Error in linear_elastic_natboun_surface:', &
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
      write(*,'(/3(a/))') 'Error in linear_elastic_natboun_surface:', &
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

  end subroutine linear_elastic_natboun_surface

end module linear_elastic_elements_3D_surface_m

