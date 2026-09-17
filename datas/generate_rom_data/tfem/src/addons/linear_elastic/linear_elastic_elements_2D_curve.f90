
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
! modulus. For the plane-stress case (2D, Cartesian), lambda and 3*K needs to be
! replaced by lambda_bar and 2*K_bar with
!
!                nu * E                  E
!  lambda_bar = ---------,    K_bar = -------
!               (1-nu**2)             2(1-nu)
!
! This module contains stuff that only works in 2D or only for curves.

! Note, that routines work in 2D with displacements in 3D:
!   * warping, where the 2D domain is a cross-section.
!   * torsion, where the 2D domain is a (z,r) cross-section (axisymmetric).

module linear_elastic_elements_2D_curve_m

  use tfem_elem_m
  use linear_elastic_set_globals_m
  use linear_elastic_elements_generic_m, &
    only: set_linear_elastic_shape_function_global
  use linear_elastic_material_m

  implicit none

contains


! Internal element routine for the linear elastic Equation
! (2D and axisymmetrical)

  subroutine linear_elastic_elem_2D ( mesh, problem, elgrp, elem, &
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


    integer :: N, M, ip, i1, i2, i3, vfuncnr, j
    real(dp) :: mu, lambda, facKmod, alpha, fac

!   set globals, gauss, shapefunctions, ...

    call set_linear_elastic_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

    call set_linear_elastic_material ( coefficients, coorsys, lemodel )

!   mixed elements needed for incompressible materials
    if ( lemodel%incompressible .and. ndfp == 0 ) &
      call errormsg_notmixed_incompressible ( 'linear_elastic_elem_2D' )

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

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( ndfp > 0 .and. any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_linear_elastic_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   set material parameters

    if ( lemodel%plane_stress ) then
      lambda = lemodel%lambda_bar
    else if ( .not. lemodel%incompressible ) then
      lambda = lemodel%lambda
    end if
    mu = lemodel%mu
    if ( lemodel%thermal_expansion) then
      alpha = lemodel%alpha
      if ( lemodel%plane_stress ) then
        facKmod = 2 * lemodel%Kmod_bar
      else if ( .not. lemodel%incompressible ) then
        facKmod = 3 * lemodel%Kmod
      end if
    end if

!   pointers in unknowns

    i1 = ndf
    i2 = 2*ndf
    i3 = 2*ndf+ndfp

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

        elemvec(1:i2) = reshape ( tmp, [ndim*ndf] )

        elemvec(i2+1:i3) = 0

      else

        elemvec(1:i3) = 0

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

          elemvec(i2+1:i3) = fp

        else

!         standard displacement elements

          do N = 1, ndf
            work = dphidx(:,N,1)
            tmp(N,1) = facKmod * alpha * sum ( work * dTemp * detF * wg )
          end do

          if ( coorsys == 1 ) then

            do N = 1, ndf
              work = dphidx(:,N,2) + phi(:,N) / xg(:,2)
              tmp(N,2) = facKmod * alpha * sum ( work * dTemp * detF * wg )
            end do

          else

            do N = 1, ndf
              work = dphidx(:,N,2)
              tmp(N,2) = facKmod * alpha * sum ( work * dTemp * detF * wg )
            end do

          end if

          elemvec(1:i2) = elemvec(1:i2) + reshape ( tmp, [ndim*ndf] )

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
          if ( coorsys == 1 ) then
            work1 = work1 + 2 * phi(:,N) * phi(:,M) / xg(:,2) ** 2
          end if
          Svv(N,M) = mu * sum ( work1 * detF * wg )
          Svv(M,N) = Svv(N,M)
        end do
      end do

      do N = 1, ndf
        do M = 1, ndf
          work = dphidx(:,N,2) * dphidx(:,M,1)
          Suv(N,M) = mu * sum ( work * detF * wg )
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
            if ( coorsys == 1 ) then
              work = dphidx(:,M,2) + phi(:,M) / xg(:,2)
              Lv(N,M) = sum ( psi(:,N) * work * detF * wg )
            else
              Lv(N,M) = sum ( psi(:,N) * dphidx(:,M,2) * detF * wg )
            end if
          end do
        end do

      else

!       standard displacement elements

        do N = 1, ndf
          do M = N, ndf
            work = dphidx(:,N,1) * dphidx(:,M,1)
            Suu(N,M) = Suu(N,M) + lambda * sum ( work * detF * wg )
            Suu(M,N) = Suu(N,M)
          end do
        end do

        if ( coorsys == 1 ) then

          do N = 1, ndf
            do M = N, ndf
              work = ( dphidx(:,N,2) + phi(:,N) / xg(:,2) ) * &
                     ( dphidx(:,M,2) + phi(:,M) / xg(:,2) )
              Svv(N,M) = Svv(N,M) + lambda * sum ( work * detF * wg )
              Svv(M,N) = Svv(N,M)
            end do
          end do

          do N = 1, ndf
            do M = 1, ndf
              work = dphidx(:,N,1) * ( dphidx(:,M,2) + phi(:,M) / xg(:,2) )
              Suv(N,M) = Suv(N,M) + lambda * sum ( work * detF * wg )
            end do
          end do

        else

          do N = 1, ndf
            do M = N, ndf
              work = dphidx(:,N,2) * dphidx(:,M,2)
              Svv(N,M) = Svv(N,M) + lambda * sum ( work * detF * wg )
              Svv(M,N) = Svv(N,M)
            end do
          end do

          do N = 1, ndf
            do M = 1, ndf
              work = dphidx(:,N,1) * dphidx(:,M,2)
              Suv(N,M) = Suv(N,M) + lambda * sum ( work * detF * wg )
            end do
          end do

        end if

      end if

      elemmat(    1:i1,    1:i1 ) = Suu
      elemmat(    1:i1, i1+1:i2 ) = Suv
      elemmat( i1+1:i2,    1:i1 ) = transpose(Suv)
      elemmat( i1+1:i2, i1+1:i2 ) = Svv

      if ( ndfp > 0 ) then
        elemmat( i2+1:i3,    1:i1 ) = -Lu
        elemmat( i2+1:i3, i1+1:i2 ) = -Lv
        elemmat(    1:i1, i2+1:i3 ) = -transpose(Lu)
        elemmat( i1+1:i2, i2+1:i3 ) = -transpose(Lv)
        elemmat( i2+1:i3, i2+1:i3)  = -Cmat
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
      allocate ( Lu(ndfp,ndf), Lv(ndfp,ndf), Cmat(ndfp,ndfp) )
      allocate ( dTemp(ninti), fp(ndfp) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( tmp )
      deallocate ( work, work1 )
      deallocate ( Suu, Suv, Svv )
      deallocate ( Lu, Lv, Cmat )
      deallocate ( dTemp, fp )

    end subroutine deallocate_arrays

  end subroutine linear_elastic_elem_2D


! Internal element routine for the linear elastic Equation
! (2D with 3D displacements)
!
  subroutine linear_elastic_elem_2D_disp3D ( mesh, problem, elgrp, elem, &
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
      coefficients, oldvectors, maxdisp3D=1 )

    call set_linear_elastic_material ( coefficients, coorsys, lemodel )

!   disp3D not set
    if ( disp3D == 0 ) &
                call errormsg_notdisp3D ( 'linear_elastic_elem_2D_disp3D' )

!   mixed elements needed for incompressible materials
    if ( lemodel%incompressible .and. ndfp == 0 ) &
      call errormsg_notmixed_incompressible ( 'linear_elastic_elem_2D_disp3D' )


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

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

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
          fg(ip,:) = coefficients%vfunc ( ncompu, vfuncnr, xg(ip,:) )
        end do

        do j = 1, ncompu
          do N = 1, ndf
            tmp(N,j) = sum ( fg(:,j) * phi(:,N) * detF * wg )
          end do
        end do

        elemvec(1:i3) = reshape ( tmp, [ncompu*ndf] )

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

          do N = 1, ndf
            work = dphidx(:,N,1)
            fuv(N,1) = facKmod * alpha * sum ( work * dTemp * detF * wg )
          end do

          if ( coorsys == 1 ) then

            do N = 1, ndf
              work = dphidx(:,N,2) + phi(:,N) / xg(:,2)
              fuv(N,2) = facKmod * alpha * sum ( work * dTemp * detF * wg )
            end do

          else

            do N = 1, ndf
              work = dphidx(:,N,2)
              fuv(N,2) = facKmod * alpha * sum ( work * dTemp * detF * wg )
            end do

          end if

          elemvec(1:i2) = elemvec(1:i2) + reshape ( fuv, [ndim*ndf] )

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
          if ( coorsys == 1 ) then
            work1 = work1 + 2 * phi(:,N) * phi(:,M) / xg(:,2) ** 2
          end if
          Svv(N,M) = mu * sum ( work1 * detF * wg )
          Svv(M,N) = Svv(N,M)
          work1 = work
          if ( coorsys == 1 ) then
            work1 = work1 + phi(:,N) * phi(:,M) / xg(:,2) ** 2  &
                     - phi(:,N) * dphidx(:,M,2) / xg(:,2)  &
                     - phi(:,M) * dphidx(:,N,2) / xg(:,2)
          end if
          Sww(N,M) = mu * sum ( work1 * detF * wg )
          Sww(M,N) = Sww(N,M)
        end do
      end do

      do N = 1, ndf
        do M = 1, ndf
          work = dphidx(:,N,2) * dphidx(:,M,1)
          Suv(N,M) = mu * sum ( work * detF * wg )
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
            if ( coorsys == 1 ) then
              work = dphidx(:,M,2) + phi(:,M) / xg(:,2)
              Lv(N,M) = sum ( psi(:,N) * work * detF * wg )
            else
              Lv(N,M) = sum ( psi(:,N) * dphidx(:,M,2) * detF * wg )
            end if
          end do
        end do

      else

!       standard displacement elements

        do N = 1, ndf
          do M = N, ndf
            work = dphidx(:,N,1) * dphidx(:,M,1)
            Suu(N,M) = Suu(N,M) + lambda * sum ( work * detF * wg )
            Suu(M,N) = Suu(N,M)
          end do
        end do

        if ( coorsys == 1 ) then

          do N = 1, ndf
            do M = N, ndf
              work = ( dphidx(:,N,2) + phi(:,N) / xg(:,2) ) * &
                     ( dphidx(:,M,2) + phi(:,M) / xg(:,2) )
              Svv(N,M) = Svv(N,M) + lambda * sum ( work * detF * wg )
              Svv(M,N) = Svv(N,M)
            end do
          end do

          do N = 1, ndf
            do M = 1, ndf
              work = dphidx(:,N,1) * ( dphidx(:,M,2) + phi(:,M) / xg(:,2) )
              Suv(N,M) = Suv(N,M) + lambda * sum ( work * detF * wg )
            end do
          end do

        else

          do N = 1, ndf
            do M = N, ndf
              work = dphidx(:,N,2) * dphidx(:,M,2)
              Svv(N,M) = Svv(N,M) + lambda * sum ( work * detF * wg )
              Svv(M,N) = Svv(N,M)
            end do
          end do

          do N = 1, ndf
            do M = 1, ndf
              work = dphidx(:,N,1) * dphidx(:,M,2)
              Suv(N,M) = Suv(N,M) + lambda * sum ( work * detF * wg )
            end do
          end do

        end if

      end if

      elemmat(    1:i1,    1:i1 ) = Suu
      elemmat(    1:i1, i1+1:i2 ) = Suv
      elemmat(    1:i1, i2+1:i3 ) = 0
      elemmat( i1+1:i2,    1:i1 ) = transpose(Suv)
      elemmat( i1+1:i2, i1+1:i2 ) = Svv
      elemmat( i1+1:i2, i2+1:i3 ) = 0
      elemmat( i2+1:i3,    1:i1 ) = 0
      elemmat( i2+1:i3, i1+1:i2 ) = 0
      elemmat( i2+1:i3, i2+1:i3 ) = Sww
      elemmat( i3+1:i4,    1:i1 ) = -Lu
      elemmat( i3+1:i4, i1+1:i2 ) = -Lv
      elemmat( i3+1:i4, i2+1:i3 ) = 0
      elemmat(    1:i1, i3+1:i4 ) = -transpose(Lu)
      elemmat( i1+1:i2, i3+1:i4 ) = -transpose(Lv)
      elemmat( i2+1:i3, i3+1:i4 ) = 0
      elemmat( i3+1:i4, i3+1:i4)  = -Cmat

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

      allocate ( tmp(ndf,ncompu) )
      allocate ( work(ninti), work1(ninti) )
      allocate ( Suu(ndf,ndf), Suv(ndf,ndf), Svv(ndf,ndf) )
      allocate ( Sww(ndf,ndf) )
      allocate ( Lu(ndfp,ndf), Lv(ndfp,ndf), Cmat(ndfp,ndfp) )
      allocate ( dTemp(ninti), fp(ndfp), fuv(ndf,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( tmp )
      deallocate ( work, work1 )
      deallocate ( Suu, Suv, Svv )
      deallocate ( Sww )
      deallocate ( Lu, Lv, Cmat )
      deallocate ( dTemp, fp, fuv )

    end subroutine deallocate_arrays

  end subroutine linear_elastic_elem_2D_disp3D


! Boundary element for a natural boundary on a curve for the linear elastic
! equation

  subroutine linear_elastic_natboun_curve ( mesh, problem, curve, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use linear_elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, ip, funcnr, vfuncnr, direction, j


    if ( first ) then

!     first element on this curve

!     set globals

      call set_globals_linear_elastic_up_boun ( mesh, coefficients, &
        curve=curve, maxdisp3D=1 )

!     allocate arrays

      allocate ( wg(ninti), fg(ninti,ncompu), curvel(ninti) )
      allocate ( tmp(ndf,ncompu) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,1), dxdxi(ninti,ndim) )

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

    vfuncnr = coefficients%i(15)
    funcnr = coefficients%i(16)
    direction = coefficients%i(17)

    if ( funcnr > 0 ) then
      fg = 0
      do ip = 1, ninti
        fg(ip,direction) = coefficients%func ( funcnr, xg(ip,:) )
      end do
    else if ( vfuncnr > 0 ) then
      do ip = 1, ninti
        fg(ip,:) = coefficients%vfunc ( ncompu, vfuncnr, xg(ip,:) )
      end do
    else
      write(*,'(/2(a/))') 'Error in linear_elastic_natboun_curve:', &
        ' either funcnr > 0 or vfuncnr > 0 '
      stop
    end if

    do j = 1, ncompu
      do N = 1, ndf
        tmp(N,j) = sum ( fg(:,j) * phi(:,N) * curvel * wg )
      end do
    end do

    elemvec = reshape ( tmp, [ncompu*ndf] )

    if ( matrix ) then
      write(*,'(/3(a/))') 'Error in linear_elastic_natboun_curve:', &
        ' matrix == .true. Make sure add_boundary_elements is called ', &
        ' with buildmatrix=.false.'
      stop
    end if

    if ( last ) then

!     last element on this curve

      deallocate ( wg, fg, curvel, normal )
      deallocate ( tmp )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi, dxdxi )

    end if

  end subroutine linear_elastic_natboun_curve

end module linear_elastic_elements_2D_curve_m
