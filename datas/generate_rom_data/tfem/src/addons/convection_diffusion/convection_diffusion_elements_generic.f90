
! Copyright (C) 2007-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the convection-diffusion-reaction equation:
!
!          dc
!  gamma ( -- + u.nabla c ) - nabla ( alpha nabla c ) + beta c = f
!          dt
!
! where u is a velocity vector and alpha, beta and gamma are coefficients
! that may depend on position.
! The coefficient alpha can be either a scalar or a tensor.

module convection_diffusion_elements_generic_m

  use tfem_elem_m
  use poisson_elements_m

  implicit none

contains


! Internal element routine for the scalar gamma*dc/dt term, discretized using
! dc/dt = (c_{n+1} - c_{n})/ Delta t
! Use the proper prefactors by setting factormat and factorvec in the heading
! of build_system() to create various high-order time integration schemes.

  subroutine scalar_dcdt_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j
    real(dp) :: deltat, fac


!   set globals, gauss, shapefunctions, ...

    call set_convection_diffusion_elem ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )


!   allocate more arrays

    if ( first .and. coefficients%i(4) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

!   determine coefficient gamma

    call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(28), value=coefficients%r(8), &
      func=coefficients%func1(3)%p, funcnr=coefficients%i(29), x=xg, indx_v=5, &
      layer=layer, phi=phi, indx_e=5, coef=gammag )

!   time step

    deltat = coefficients%r(3)
    fac = 1 / deltat

    if ( vector ) then

!     get the scalar at the old time step

      call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, &
        u=cn, layer=layer )

      cng = matmul ( phi, cn )

!     rhs = gamma * cn / delta

      do i = 1, ndf
        elemvec(i) = fac * sum ( gammag * phi(:,i) * cng * detF * wg )
      end do

    end if

    if ( matrix ) then

!     mass matrix

      do i = 1, ndf
        do j = i, ndf
          elemmat(i,j) = fac * sum ( gammag * phi(:,i) * phi(:,j) * detF * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_convection_diffusion_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(4) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if


  contains

    subroutine allocate_arrays

      allocate ( gammag(ninti) )
      allocate ( cn(ndf), cng(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( gammag )
      deallocate ( cn, cng )

    end subroutine deallocate_arrays

  end subroutine scalar_dcdt_elem


! Internal element routine for the scalar gamma*u.nabla c term (matrix only)

  subroutine scalar_ugradc_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip


!   set globals, gauss, shapefunctions, ...

    call set_convection_diffusion_elem ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )


!   allocate more arrays

    if ( first .and. coefficients%i(4) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   determine coefficient gamma

    call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(28), value=coefficients%r(8), &
      func=coefficients%func1(3)%p, funcnr=coefficients%i(29), x=xg, indx_v=5, &
      layer=layer, phi=phi, indx_e=5, coef=gammag )

!   determine coefficient u

    call evaluate_vector_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(24), value=coefficients%r(4:3+ndim), &
      vfunc=coefficients%vfunc1(1)%p, vfuncnr=coefficients%i(25), x=xg, &
      indx_v=3, layer=layer, phi=phi, indx_e=3, coef=uvecg )

!   u.grad operator

    do ip = 1, ninti
      ugradphi(ip,:) = matmul ( dphidx(ip,:,:), uvecg(ip,:) )
    end do

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     implicit matrix

      do i = 1, ndf
        do j = 1, ndf
          elemmat(i,j) = sum ( gammag * phi(:,i) * ugradphi(:,j) * detF * wg )
        end do
      end do

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_convection_diffusion_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(4) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if


  contains

    subroutine allocate_arrays

      allocate ( ugradphi(ninti,ndf) )
      allocate ( gammag(ninti), uvecg(ninti,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( ugradphi )
      deallocate ( gammag, uvecg )

    end subroutine deallocate_arrays

  end subroutine scalar_ugradc_elem


! Internal element routine for the scalar -gamma*u.nabla cn term (vector only)

  subroutine scalar_ugradcn_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, ip


!   set globals, gauss, shapefunctions, ...

    call set_convection_diffusion_elem ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )


!   allocate more arrays

    if ( first .and. coefficients%i(4) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   determine coefficient gamma

    call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(28), value=coefficients%r(8), &
      func=coefficients%func1(3)%p, funcnr=coefficients%i(29), x=xg, indx_v=5, &
      layer=layer, phi=phi, indx_e=5, coef=gammag )

!   determine coefficient u

    call evaluate_vector_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(24), value=coefficients%r(4:3+ndim), &
      vfunc=coefficients%vfunc1(1)%p, vfuncnr=coefficients%i(25), x=xg, &
      indx_v=3, layer=layer, phi=phi, indx_e=3, coef=uvecg )

!   u.grad operator

    do ip = 1, ninti
      ugradphi(ip,:) = matmul ( dphidx(ip,:,:), uvecg(ip,:) )
    end do

!   get the scalar at the old time step

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, &
      u=cn, layer=layer )

!   u.grad cn

    ugradcn = matmul ( ugradphi, cn )

    if ( vector ) then

!     explicit vector

      do i = 1, ndf
        elemvec(i) = - sum ( gammag * phi(:,i) * ugradcn * detF * wg )
      end do

    end if

    if ( matrix ) then

      elemmat = 0

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_convection_diffusion_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(4) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if


  contains

    subroutine allocate_arrays

      allocate ( cn(ndf), ugradphi(ninti,ndf), ugradcn(ninti) )
      allocate ( gammag(ninti), uvecg(ninti,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( cn, ugradphi, ugradcn )
      deallocate ( gammag, uvecg )

    end subroutine deallocate_arrays

  end subroutine scalar_ugradcn_elem


! Internal element routine for the scalar beta*c term. (matrix only)

  subroutine scalar_betac_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use convection_diffusion_globals_m

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

    call set_convection_diffusion_elem ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )


!   allocate more arrays

    if ( first .and. coefficients%i(4) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

!   determine coefficient beta

    call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(26), value=coefficients%r(7), &
      func=coefficients%func1(2)%p, funcnr=coefficients%i(27), x=xg, indx_v=4, &
      layer=layer, phi=phi, indx_e=4, coef=betag )

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     mass matrix

      do i = 1, ndf
        do j = i, ndf
          elemmat(i,j) = sum ( betag * phi(:,i) * phi(:,j) * detF * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_convection_diffusion_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(4) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if


  contains

    subroutine allocate_arrays

      allocate ( betag(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( betag )

    end subroutine deallocate_arrays

  end subroutine scalar_betac_elem


! Internal element routine for the scalar -beta*cn term. (vector only)

  subroutine scalar_betacn_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i


!   set globals, gauss, shapefunctions, ...

    call set_convection_diffusion_elem ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )


!   allocate more arrays

    if ( first .and. coefficients%i(4) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

!   determine coefficient beta

    call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(26), value=coefficients%r(7), &
      func=coefficients%func1(2)%p, funcnr=coefficients%i(27), x=xg, indx_v=4, &
      layer=layer, phi=phi, indx_e=4, coef=betag )

    if ( vector ) then

!     get the scalar at the old time step

      call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, &
        u=cn, layer=layer )

      cng = matmul ( phi, cn )

!     rhs = - beta * cn

      do i = 1, ndf
        elemvec(i) = - sum ( betag * phi(:,i) * cng * detF * wg )
      end do

    end if

    if ( matrix ) then

      elemmat = 0

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_convection_diffusion_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(4) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if


  contains

    subroutine allocate_arrays

      allocate ( betag(ninti) )
      allocate ( cn(ndf), cng(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( betag )
      deallocate ( cn, cng )

    end subroutine deallocate_arrays

  end subroutine scalar_betacn_elem


! Internal element routine for the scalar diffusion equation with varying
! coefficient alpha and (optionally) varying Gauss rule per element;
!   - nabla ( alpha nabla c ) = f
! Coefficient alpha is a scalar.

  subroutine scalar_diffusion_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip, funcnr


!   set globals, gauss, shapefunctions, ...

    call set_convection_diffusion_elem ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )


!   allocate more arrays

    if ( first .and. coefficients%i(4) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( vector ) then

      if ( coefficients%i(12) > 0 ) then

!       function number specified (old method)

        funcnr = coefficients%i(12)

        do ip = 1, ninti
          fg(ip) = coefficients%func ( funcnr, xg(ip,:) )
        end do

        do i = 1, ndf
          elemvec(i) = sum ( fg * phi(:,i) * detF * wg )
        end do

      else if ( coefficients%i(18) > 0 ) then

!       evaluate right-hand side

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(18), value=0._dp, &
          func=coefficients%func, funcnr=coefficients%i(19), x=xg, indx_v=2, &
          layer=layer, phi=phi, indx_e=2, coef=fg )

        do i = 1, ndf
          elemvec(i) = sum ( fg * phi(:,i) * detF * wg )
        end do

      else

        elemvec = 0

      end if

    end if

    if ( matrix ) then

!     determine coefficient alpha

      call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
        elgrp, elem, choice=coefficients%i(2), value=coefficients%r(1), &
        func=coefficients%func1(1)%p, funcnr=coefficients%i(3), x=xg, indx_v=1,&
        layer=layer, phi=phi, indx_e=1, coef=alphag )

!     fill matrix

      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = alphag(ip) * sum ( dphidx(ip,i,:) * dphidx(ip,j,:) )
          end do
          elemmat(i,j) = sum ( work * detF * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_convection_diffusion_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(4) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if


  contains

    subroutine allocate_arrays

      allocate ( work(ninti) )
      allocate ( alphag(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work )
      deallocate ( alphag )

    end subroutine deallocate_arrays

  end subroutine scalar_diffusion_elem


! Internal element routine for the scalar diffusion term with varying
! coefficient alpha and (optionally) varying Gauss rule per element:
!    ... = nabla ( alpha nabla cn ) + ..
! Vector (right-hand side) only with cn substituted.
! Coefficient alpha is a scalar.

  subroutine scalar_diffusion_rhs_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, ip


!   set globals, gauss, shapefunctions, ...

    call set_convection_diffusion_elem ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )


!   allocate more arrays

    if ( first .and. coefficients%i(4) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   get the scalar at the old time step

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, &
      u=cn, layer=layer )

    do i = 1, ndim
      gradcn(:,i) = matmul ( dphidx(:,:,i), cn )
    end do

!   determine coefficient alpha

    call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(2), value=coefficients%r(1), &
      func=coefficients%func1(1)%p, funcnr=coefficients%i(3), x=xg, indx_v=1, &
      layer=layer, phi=phi, indx_e=1, coef=alphag )

    if ( vector ) then

      do i = 1, ndf
        do ip = 1, ninti
          work(ip) = alphag(ip) * sum ( dphidx(ip,i,:) * gradcn(ip,:) )
        end do
        elemvec(i) = - sum ( work * detF * wg )
      end do

    end if

    if ( matrix ) then

      elemmat = 0

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_convection_diffusion_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(4) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if


  contains

    subroutine allocate_arrays

      allocate ( work(ninti) )
      allocate ( alphag(ninti) )
      allocate ( gradcn(ninti,ndim), cn(ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work )
      deallocate ( alphag )
      deallocate ( gradcn, cn )

    end subroutine deallocate_arrays

  end subroutine scalar_diffusion_rhs_elem


! Internal element routine for the scalar diffusion equation with varying
! coefficient alpha and (optionally) varying Gauss rule per element;
!   - nabla ( alpha nabla c ) = f
! Coefficient alpha is a tensor.

  subroutine scalar_diffusion2_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip, funcnr


!   set globals, gauss, shapefunctions, ...

    call set_convection_diffusion_elem ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )


!   allocate more arrays

    if ( first .and. coefficients%i(4) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( vector ) then

      if ( coefficients%i(12) > 0 ) then

!       function number specified (old method)

        funcnr = coefficients%i(12)

        do ip = 1, ninti
          fg(ip) = coefficients%func ( funcnr, xg(ip,:) )
        end do

        do i = 1, ndf
          elemvec(i) = sum ( fg * phi(:,i) * detF * wg )
        end do

      else if ( coefficients%i(18) > 0 ) then

!       evaluate right-hand side

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(18), value=0._dp, &
          func=coefficients%func, funcnr=coefficients%i(19), x=xg, indx_v=2, &
          layer=layer, phi=phi, indx_e=2, coef=fg )

        do i = 1, ndf
          elemvec(i) = sum ( fg * phi(:,i) * detF * wg )
        end do

      else

        elemvec = 0

      end if

    end if


    if ( matrix ) then

!     determine coefficient alpha

      call evaluate_tensor_coefficient ( mesh, problem, oldvectors, &
        elgrp, elem, choice=coefficients%i(2), &
        value=reshape(coefficients%r(9:8+ndim**2), [ndim,ndim]), &
        tfunc=coefficients%tfunc1(1)%p, tfuncnr=coefficients%i(3), x=xg, &
        indx_v=1, layer=layer, phi=phi, indx_e=1, coef=alphatg )

!     fill matrix

      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,i,:) * &
                               matmul ( alphatg(ip,:,:), dphidx(ip,j,:) ) )
          end do
          elemmat(i,j) = sum ( work * detF * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_convection_diffusion_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(4) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if


  contains

    subroutine allocate_arrays

      allocate ( work(ninti) )
      allocate ( alphatg(ninti,ndim,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work )
      deallocate ( alphatg )

    end subroutine deallocate_arrays

  end subroutine scalar_diffusion2_elem


! Internal element routine for the scalar diffusion term with varying
! coefficient alpha and (optionally) varying Gauss rule per element:
!    ... = nabla ( alpha nabla cn ) + ..
! Vector (right-hand side) only with cn substituted.
! Coefficient alpha is a tensor.

  subroutine scalar_diffusion2_rhs_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, ip


!   set globals, gauss, shapefunctions, ...

    call set_convection_diffusion_elem ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )


!   allocate more arrays

    if ( first .and. coefficients%i(4) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   get the scalar at the old time step

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, &
      u=cn, layer=layer )

    do i = 1, ndim
      gradcn(:,i) = matmul ( dphidx(:,:,i), cn )
    end do

!   determine coefficient alpha

    call evaluate_tensor_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(2), &
      value=reshape(coefficients%r(9:8+ndim**2), [ndim,ndim]), &
      tfunc=coefficients%tfunc1(1)%p, tfuncnr=coefficients%i(3), x=xg, &
      indx_v=1, layer=layer, phi=phi, indx_e=1, coef=alphatg )

    if ( vector ) then

      do i = 1, ndf
        do ip = 1, ninti
          work(ip) = sum ( dphidx(ip,i,:) * &
                            matmul ( alphatg(ip,:,:), gradcn(ip,:) ) )
        end do
        elemvec(i) = - sum ( work * detF * wg )
      end do

    end if

    if ( matrix ) then

      elemmat = 0

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_convection_diffusion_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(4) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if


  contains

    subroutine allocate_arrays

      allocate ( work(ninti) )
      allocate ( alphatg(ninti,ndim,ndim) )
      allocate ( gradcn(ninti,ndim), cn(ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work )
      deallocate ( alphatg )
      deallocate ( gradcn, cn )

    end subroutine deallocate_arrays

  end subroutine scalar_diffusion2_rhs_elem


! set global parameters (internal element)

  subroutine set_globals_convection_diffusion ( mesh, coefficients, elgrp )

    use poisson_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

!   check size of coefficients

    call check ( coefficients, 'set_globals_convection_diffusion', ncoefi=100, &
      ncoefr=50, indexarray=[2,4], minimum=[0,0], maximum=[3,2] )

  end subroutine set_globals_convection_diffusion


! Preamble for the convection diffusion element

  subroutine set_convection_diffusion_elem ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors )

    use convection_diffusion_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors

    integer :: i1(3)

    if ( first .and. coefficients%i(4) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

!     set globals

      call set_globals_poisson ( mesh, coefficients, elgrp )
      call set_globals_convection_diffusion ( mesh, coefficients, elgrp )

!     allocate arrays

      allocate ( x(nodalp,ndim) )

      allocate ( xig(ninti,ndim), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim) )

      allocate ( detF(ninti), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( fg(ninti), xg(ninti,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      if ( first ) then

!       first element in this group

!       set globals

        call set_globals_poisson ( mesh, coefficients, elgrp )
        call set_globals_convection_diffusion ( mesh, coefficients, elgrp )

        allocate ( x(nodalp,ndim) )

      end if

!     now choose between elvector and user subroutines

      if ( coefficients%i(4) == 1 ) then

!       get element values of intrule, intrule2 and nsubint

        if ( gauss%globalshape == 'prism' .or. &
             gauss%globalshape == 'pyramid' .and. gauss%inttype == 2 ) then

!         intrule2 required
          call get_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, i1=i1 )

          if ( coefficients%i(50) == 1 .or. SET_GAUSS_BY_ORDER ) then
            intrule2 = set_intrule2 ( gauss%globalshape, gauss%inttype, &
                                      order=i1(3) )
          else
            intrule2 = i1(3)
          end if

          gauss%intrule2 = intrule2

        else

          call get_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, i1=i1(1:2) )

        end if

        if ( coefficients%i(50) == 1 .or. SET_GAUSS_BY_ORDER ) then
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

!       set Gauss integration and shape function

        call set_Gauss_integration ( gauss, xig, wg )

      else if ( coefficients%i(4) == 2 ) then

!       set ninti

        call coefficients%set_ninti_user ( mesh, problem, elgrp, elem, first, &
          last, coefficients, oldvectors )

!       allocate arrays

        allocate ( xig(ninti,ndim), wg(ninti) )

!       set Gauss integration

        call coefficients%set_Gauss_integration_user ( mesh, problem, elgrp, &
          elem, first, last, coefficients, oldvectors )

      end if

      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim) )

      allocate ( detF(ninti), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( fg(ninti), xg(ninti,ndim) )

!     set shape function

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

  end subroutine set_convection_diffusion_elem


! Unset the preamble for the convection-diffusion element
! (deallocate arrays allocated in set_... )

  subroutine unset_convection_diffusion_elem ( last, coefficients )

    use convection_diffusion_globals_m

    logical, intent(in) :: last
    type(coefficients_t), intent(in) :: coefficients

!   deallocate memory

    if ( last .and. coefficients%i(4) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      deallocate ( x )

      deallocate ( xig, wg )
      deallocate ( phi, dphi )

      deallocate ( detF, F )
      deallocate ( Finv, dphidx )
      deallocate ( fg, xg )

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      if ( last ) then

!       last element in this group

        deallocate ( x )

      end if

      deallocate ( xig, wg )
      deallocate ( phi, dphi )

      deallocate ( detF, F )
      deallocate ( Finv, dphidx )
      deallocate ( fg, xg )

    end if

  end subroutine unset_convection_diffusion_elem

end module convection_diffusion_elements_generic_m

