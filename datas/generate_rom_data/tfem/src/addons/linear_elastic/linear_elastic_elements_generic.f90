
! Copyright (C) 2020-2021 Martien A. Hulsen
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
! This module contains the generic stuff (works for 2D and 3D).

module linear_elastic_elements_generic_m

  use tfem_elem_m
  use linear_elastic_material_m

  implicit none

contains


! Element for the constraints (connection through collocation)

  subroutine linear_elastic_constr_node_conn ( mesh, problem, constr, elem, &
    node, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

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

  end subroutine linear_elastic_constr_node_conn


! compute displacement in nodal points (when some nodes have zero displacement)

  subroutine linear_elastic_displacement ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use linear_elastic_globals_m

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

      call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
        maxdisp3D=1 )

!     allocate arrays

      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf) )
      allocate ( u(ncompu*ndf), tmp(ndf,ncompu) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi )

    end if

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqdisp], layer=layer )

    tmp = reshape ( u, [ndf,ncompu] )

    elemvec = reshape ( matmul ( phi, tmp ), [nodalp*ncompu] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( xrnod, phi )
      deallocate ( u, tmp )

    end if

  end subroutine linear_elastic_displacement


! compute displacement derivatives

  subroutine linear_elastic_deriv ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use linear_elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: comp, j, maxderiv


    if ( first ) then

!     first element in this group

      call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
        maxdisp3D=1 )

      call set_linear_elastic_material ( coefficients, coorsys, lemodel )

!     check component = coefficients%i(13)

      maxderiv = 1

      call check ( coefficients, 'linear_elastic_deriv', indexarray=[13], &
        minimum=[1], maximum=[maxderiv] )

!     allocate arrays

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ncompu*ndf), uvector(ndf,ncompu) )
      allocate ( gradu(nodalp,ncompu,ncompu) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( dTemp(nodalp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqdisp], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( disp3D == 1 ) then

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

    comp = coefficients%i(13)

    if ( coorsys <= 1 .and. disp3D == 0 ) then

!     2D

      if ( lemodel%plane_stress .and. lemodel%thermal_expansion ) then

!       thermal expansion: get temperature change

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(9), value=coefficients%r(4), &
          func=coefficients%func1(1)%p, funcnr=coefficients%i(8), x=x, &
          indx_v=1, layer=layer, phi=phi, indx_e=1, coef=dTemp )

      end if

      select case(comp)
      case(1)

        elemvec = gradu(:,1,1) + gradu(:,2,2) ! dudx + dudy = divergence

        if ( lemodel%plane_stress ) then

!         plane stress

          elemvec = elemvec &
             - ( gradu(:,1,1) + gradu(:,2,2) ) * lemodel%nu / (1-lemodel%nu)

          if ( lemodel%thermal_expansion ) then

!           thermal expansion

            elemvec = elemvec &
                   + lemodel%alpha * dTemp * (1+lemodel%nu) / (1-lemodel%nu)

          end if

        else if ( coorsys == 1 ) then ! axisymmetric

          where ( x(:,2) < 1e-10_dp )
            elemvec = elemvec + gradu(:,2,2)             ! du_r / dr (r=0)
          elsewhere
            elemvec = elemvec + uvector(:,2) / x(:,2)    ! u_r / r
          end where

        end if

      case default

        call errormsg_case_default ( 'linear_elastic_deriv', &
          'comp', int_value=comp )

      end select

    else if ( coorsys == 2 .or. disp3D == 1 ) then

!     3D or disp3D = 1

      select case(comp)
      case(1)
        elemvec = gradu(:,1,1) + gradu(:,2,2) + gradu(:,3,3) ! divergence
      case default
        call errormsg_case_default ( 'linear_elastic_deriv', &
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
      deallocate ( dTemp )

    end if

  end subroutine linear_elastic_deriv


! epsilon strain tensor

  subroutine linear_elastic_strain_tensor ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use linear_elastic_globals_m

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

      call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
        maxdisp3D=1 )

      call set_linear_elastic_material ( coefficients, coorsys, lemodel )

!     test size of elemvec

      if ( coorsys <= 1 .and. disp3D == 0 ) then
        maxcomp = 4
      else
        maxcomp = 6
      end if

      if ( size(elemvec) /= nodalp*maxcomp ) then
        write(*,'(/2(a/))') 'Error in linear_elastic_strain_tensor:', &
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
      allocate ( dTemp(nodalp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqdisp], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( disp3D == 1 ) then

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

    if ( coorsys <= 1 .and. disp3D == 0 ) then

      work2(:,1) = gradu(:,1,1)                        ! epsilon_xx
      work2(:,2) = ( gradu(:,1,2) + gradu(:,2,1) ) / 2 ! epsilon_xy
      work2(:,3) = gradu(:,2,2)                        ! epsilon_yy

      if ( lemodel%plane_stress ) then

!       plane stress
        work2(:,4) = &
             - ( gradu(:,1,1) + gradu(:,2,2) ) * lemodel%nu / (1-lemodel%nu)

        if ( lemodel%thermal_expansion ) then

!         thermal expansion

          call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
            elgrp, elem, choice=coefficients%i(9), value=coefficients%r(4), &
            func=coefficients%func1(1)%p, funcnr=coefficients%i(8), x=x, &
            indx_v=1, layer=layer, phi=phi, indx_e=1, coef=dTemp )

            work2(:,4) = work2(:,4) &
                  + lemodel%alpha * dTemp * (1+lemodel%nu) / (1-lemodel%nu)

        end if

      else if ( coorsys == 1 ) then

!       axisymmetric

        where ( x(:,2) < 1e-10_dp )
          work2(:,4) = gradu(:,2,2)           ! du_r / dr (r=0)
        elsewhere
          work2(:,4) = uvector(:,2) / x(:,2)  ! u_r / r
        end where

      else

!       plain strain

         work2(:,4) = 0

      end if

    else if ( coorsys == 2 .or. disp3D == 1 ) then

!     3D or disp3D = 1

      work2(:,1) = gradu(:,1,1)                        ! epsilon_xx
      work2(:,2) = ( gradu(:,1,2) + gradu(:,2,1) ) / 2 ! epsilon_xy
      work2(:,3) = ( gradu(:,1,3) + gradu(:,3,1) ) / 2 ! epsilon_xz
      work2(:,4) = gradu(:,2,2)                        ! epsilon_yy
      work2(:,5) = ( gradu(:,2,3) + gradu(:,3,2) ) / 2 ! epsilon_yz
      work2(:,6) = gradu(:,3,3)                        ! epsilon_zz

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
      deallocate ( dTemp )

    end if

  end subroutine linear_elastic_strain_tensor


! displacement gradient tensor

  subroutine linear_elastic_gradu_tensor ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use linear_elastic_globals_m

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

      call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
        maxdisp3D=1 )

      call set_linear_elastic_material ( coefficients, coorsys, lemodel )

!     test size of elemvec

      if ( coorsys <= 1 .and. disp3D == 0 ) then
        maxcomp = 5
      else
        maxcomp = 9
      end if

      if ( size(elemvec) /= nodalp*maxcomp ) then
        write(*,'(/2(a/))') 'Error in linear_elastic_gradu_tensor:', &
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
      allocate ( dTemp(nodalp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqdisp], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( disp3D == 1 ) then

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

    if ( coorsys <= 1 .and. disp3D == 0 ) then

      work2(:,1) = gradu(:,1,1) ! dudx
      work2(:,2) = gradu(:,1,2) ! dudy
      work2(:,3) = gradu(:,2,1) ! dvdx
      work2(:,4) = gradu(:,2,2) ! dvdy

      if ( lemodel%plane_stress ) then

!       plane stress

         work2(:,5) = &
             - ( gradu(:,1,1) + gradu(:,2,2) ) * lemodel%nu / (1-lemodel%nu)

        if ( lemodel%thermal_expansion ) then

!         thermal expansion

          call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
            elgrp, elem, choice=coefficients%i(9), value=coefficients%r(4), &
            func=coefficients%func1(1)%p, funcnr=coefficients%i(8), x=x, &
            indx_v=1, layer=layer, phi=phi, indx_e=1, coef=dTemp )

            work2(:,5) = work2(:,5) &
                  + lemodel%alpha * dTemp * (1+lemodel%nu) / (1-lemodel%nu)

        end if

      else if ( coorsys == 1 ) then

!       axisymmetric

        where ( x(:,2) < 1e-10_dp )
          work2(:,5) = gradu(:,2,2)           ! du_r / dr (r=0)
        elsewhere
          work2(:,5) = uvector(:,2) / x(:,2)  ! u_r / r
        end where

      else

!       plane strain

        work2(:,5) = 0

      end if

    else if ( coorsys == 2 .or. disp3D == 1 ) then

!    3D or disp3D = 1

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
      deallocate ( dTemp )

    end if

  end subroutine linear_elastic_gradu_tensor


! linear_elastic stress deriv

  subroutine linear_elastic_stress_deriv ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use linear_elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer, save :: maxstress
    integer :: j, maxderiv, comp

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
        maxdisp3D=1 )

      call set_linear_elastic_material ( coefficients, coorsys, lemodel )

!     test size of elemvec

      maxderiv = 2

      call check ( coefficients, 'linear_elastic_stress_deriv', &
        indexarray=[13], minimum=[1], maximum=[maxderiv] )

!     set maxstress

      if ( coorsys <= 1 .and. disp3D == 0 ) then
        maxstress = 4
      else
        maxstress = 6
      end if

!     allocate arrays

      allocate ( detF(nodalp), work(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ncompu*ndf), uvector(ndf,ncompu) )
      allocate ( gradu(nodalp,ncompu,ncompu) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( work2(nodalp,maxstress) )
      allocate ( psi(nodalp,ndfp), prs(ndfp), pr(nodalp) )
      allocate ( dTemp(nodalp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )
      if ( ndfp > 0 .and. coefficients%i(62) == 0 ) then
        call set_shape_function ( shapefuncp, xrnod, psi )
      end if

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqdisp], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( disp3D == 1 ) then

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

    if ( ndfp > 0 ) then

!     mixed elements

      if ( any ( coefficients%i(62) == [1,2] ) ) then
!       global shape function for pressure
        call set_linear_elastic_shape_function_global ( shapefuncp, x, x, psi, &
          coefficients%i(62)==1 )
      end if

      call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, prs, &
        physq=[physqpress], layer=layer )

      pr = matmul ( psi, prs )

    end if

    if ( lemodel%thermal_expansion .and. ndfp == 0 ) then

!     thermal expansion: get temperature change

      call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
        elgrp, elem, choice=coefficients%i(9), value=coefficients%r(4), &
        func=coefficients%func1(1)%p, funcnr=coefficients%i(8), x=x, &
        indx_v=1, layer=layer, phi=phi, indx_e=1, coef=dTemp )

    end if

!   compute stress tensor

    call le_stress_tensor ( x, work2 )

!   deriv

    comp = coefficients%i(13)

    if ( coorsys <= 1 .and. disp3D == 0 ) then

!     2D

      select case(comp)
      case(1) ! hydrostatic stress tr(sigma)/3
        elemvec =  ( work2(:,1) + work2(:,3) + work2(:,4) ) / 3
      case(2) ! von Mises stress
        elemvec = ( ( work2(:,1) - work2(:,3) ) ** 2 + &
                    ( work2(:,3) - work2(:,4) ) ** 2 + &
                    ( work2(:,4) - work2(:,1) ) ** 2 + &
                    6 * work2(:,2) ** 2 )
        elemvec = sqrt( elemvec / 2 )
      case default
        call errormsg_case_default ( 'linear_elastic_stress_deriv', &
          'comp', int_value=comp )
      end select

    else if ( coorsys == 2 .or. disp3D == 1 ) then

!     3D or disp3D = 1

      select case(comp)
      case(1) ! hydrostatic stress tr(sigma)/3
        elemvec =  ( work2(:,1) + work2(:,4) + work2(:,6) ) / 3
      case(2) ! von Mises stress
        elemvec = ( ( work2(:,1) - work2(:,4) ) ** 2 + &
                    ( work2(:,4) - work2(:,6) ) ** 2 + &
                    ( work2(:,6) - work2(:,1) ) ** 2 + &
                    6 * work2(:,2) ** 2 + &
                    6 * work2(:,3) ** 2 + &
                    6 * work2(:,5) ** 2 )
        elemvec = sqrt( elemvec / 2 )
      case default
        call errormsg_case_default ( 'linear_elastic_stress_deriv', &
          'comp', int_value=comp )
      end select

    end if

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( detF, work )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( work2 )
      deallocate ( psi, prs, pr )
      deallocate ( dTemp )

    end if

  end subroutine linear_elastic_stress_deriv


! linear_elastic stress tensor

  subroutine linear_elastic_stress_tensor ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use linear_elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer, save :: maxstress
    integer :: j

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
        maxdisp3D=1 )

      call set_linear_elastic_material ( coefficients, coorsys, lemodel )

!     test size of elemvec

      if ( coorsys <= 1 .and. disp3D == 0 ) then
        maxstress = 4
      else
        maxstress = 6
      end if

      if ( size(elemvec) /= nodalp*maxstress ) then
        write(*,'(/2(a/))') 'Error in linear_elastic_stress_tensor:', &
          ' element vector has incorrect size for a symmetric tensor'
        stop
      end if

!     allocate arrays

      allocate ( detF(nodalp), work(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ncompu*ndf), uvector(ndf,ncompu) )
      allocate ( gradu(nodalp,ncompu,ncompu) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( work2(nodalp,maxstress) )
      allocate ( psi(nodalp,ndfp), prs(ndfp), pr(nodalp) )
      allocate ( dTemp(nodalp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )
      if ( ndfp > 0 .and. coefficients%i(62) == 0 ) then
        call set_shape_function ( shapefuncp, xrnod, psi )
      end if

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqdisp], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( disp3D == 1 ) then

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

    if ( ndfp > 0 ) then

!     mixed elements

      if ( any ( coefficients%i(62) == [1,2] ) ) then
!       global shape function for pressure
        call set_linear_elastic_shape_function_global ( shapefuncp, x, x, psi, &
          coefficients%i(62)==1 )
      end if

      call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, prs, &
        physq=[physqpress], layer=layer )

      pr = matmul ( psi, prs )

    end if

    if ( lemodel%thermal_expansion .and. ndfp == 0 ) then

!     thermal expansion: get temperature change

      call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
        elgrp, elem, choice=coefficients%i(9), value=coefficients%r(4), &
        func=coefficients%func1(1)%p, funcnr=coefficients%i(8), x=x, &
        indx_v=1, layer=layer, phi=phi, indx_e=1, coef=dTemp )

    end if

!   compute stress tensor

    call le_stress_tensor ( x, work2 )

    elemvec = reshape ( work2, [nodalp*maxstress] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( detF, work )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( work2 )
      deallocate ( psi, prs, pr )
      deallocate ( dTemp )

    end if

  end subroutine linear_elastic_stress_tensor


! Compute stress tensor for linear elastic material.
! NOTE, that arguments are supplied via module linear_elastic_globals_m:
!     input: lemodel, coorsys, disp3D, gradu, uvector, ndfp, work, pr, dTemp

  subroutine le_stress_tensor ( xs, sigma )

    use linear_elastic_globals_m

!   the coordinates where the stress tensor is calculated
    real(dp), dimension(:,:), intent(in) :: xs

!   output of the stress tensor
    real(dp), dimension(:,:), intent(out) :: sigma

    real(dp) :: lambda, mu, facKmod, alpha

!   set material parameters

    if ( lemodel%plane_stress ) then
      lambda = lemodel%lambda_bar
      facKmod = 2 * lemodel%Kmod_bar
    else if ( .not. lemodel%incompressible ) then
      lambda = lemodel%lambda
      facKmod = 3 * lemodel%Kmod
    end if
    mu = lemodel%mu
    if ( lemodel%thermal_expansion) alpha = lemodel%alpha

!   first the 2*mu*epsilon part of the stress tensor

    if ( coorsys <= 1 .and. disp3D == 0 ) then

      sigma(:,1) = 2 * mu * gradu(:,1,1)                ! xx
      sigma(:,2) = mu * ( gradu(:,1,2) + gradu(:,2,1) ) ! xy
      sigma(:,3) = 2 * mu * gradu(:,2,2)                ! yy

      if ( coorsys == 1 ) then
!       axisymmetric
        where ( xs(:,2) < 1e-10_dp )
          sigma(:,4) = 2 * mu * gradu(:,2,2)           ! 2 mu du_r / dr (r=0)
        elsewhere
          sigma(:,4) = 2 * mu * uvector(:,2) / xs(:,2)  ! 2 mu u_r / r
        end where
      end if

    else if ( coorsys == 2 .or. disp3D == 1 ) then

      sigma(:,1) = 2 * mu * gradu(:,1,1)                ! xx
      sigma(:,2) = mu * ( gradu(:,1,2) + gradu(:,2,1) ) ! xy
      sigma(:,3) = mu * ( gradu(:,1,3) + gradu(:,3,1) ) ! xz
      sigma(:,4) = 2 * mu * gradu(:,2,2)                ! yy
      sigma(:,5) = mu * ( gradu(:,2,3) + gradu(:,3,2) ) ! yz
      sigma(:,6) = 2 * mu * gradu(:,3,3)                ! zz

    end if

!   Now add the "pressure" part. For mixed (u,p) elements this part is in the
!   pressure and must be computed separately.

    if ( ndfp == 0 ) then

!     compute pressure part directly from the strain

      if ( coorsys <= 1 .and. disp3D == 0 ) then

        work = gradu(:,1,1) + gradu(:,2,2)

        if ( coorsys == 1 ) then
!         axisymmetric
          where ( xs(:,2) < 1e-10_dp )
            work = work + gradu(:,2,2)
          elsewhere
            work = work + uvector(:,2) / xs(:,2)
          end where
        end if

        sigma(:,1) = sigma(:,1) + lambda * work ! xx
        sigma(:,3) = sigma(:,3) + lambda * work ! yy

        if ( lemodel%plane_stress ) then
          sigma(:,4) = 0
        else if ( coorsys == 1 ) then
  !       axisymmetric
          sigma(:,4) = sigma(:,4) + lambda * work
        else
  !       plane strain
          sigma(:,4) = lambda * work
        end if

        if ( lemodel%thermal_expansion ) then

!         add thermal expansion term

          sigma(:,1) = sigma(:,1) - facKmod * alpha * dTemp ! xx
          sigma(:,3) = sigma(:,3) - facKmod * alpha * dTemp ! yy

          if ( .not. lemodel%plane_stress ) then
  !         axisymmetric or plane strain
            sigma(:,4) = sigma(:,4) - facKmod * alpha * dTemp
          end if

        end if

      else if ( coorsys == 2 .or. disp3D == 1 ) then

        work = gradu(:,1,1) + gradu(:,2,2) + gradu(:,3,3)

        sigma(:,1) = sigma(:,1) + lambda * work
        sigma(:,4) = sigma(:,4) + lambda * work
        sigma(:,6) = sigma(:,6) + lambda * work

        if ( lemodel%thermal_expansion ) then

!         add thermal expansion term

          sigma(:,1) = sigma(:,1) - facKmod * alpha * dTemp ! xx
          sigma(:,4) = sigma(:,4) - facKmod * alpha * dTemp ! yy
          sigma(:,6) = sigma(:,6) - facKmod * alpha * dTemp ! zz

        end if

      end if

    else if ( ndfp > 0 ) then

!     mixed (u,p) elements

      if ( coorsys <= 1 .and. disp3D == 0 ) then

        sigma(:,1) = sigma(:,1) - pr ! xx
        sigma(:,3) = sigma(:,3) - pr ! yy

        if ( lemodel%plane_stress ) then
          sigma(:,4) = 0
        else if ( coorsys == 1 ) then
  !       axisymmetric
          sigma(:,4) = sigma(:,4) - pr
        else
  !       plane strain
          sigma(:,4) = - pr
        end if

      else if ( coorsys == 2 .or. disp3D == 1 ) then

        sigma(:,1) = sigma(:,1) - pr
        sigma(:,4) = sigma(:,4) - pr
        sigma(:,6) = sigma(:,6) - pr

      end if

    end if

  end subroutine le_stress_tensor


! compute pressures in all nodes

  subroutine linear_elastic_pressure ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use linear_elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    if ( first ) then

!     first element in this group

      call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
        maxdisp3D=1 )

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
      call set_linear_elastic_shape_function_global ( shapefuncp, x, x, psi, &
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

  end subroutine linear_elastic_pressure


! sample value of displacement in one node of the object

  subroutine linear_elastic_sample_displacement ( mesh, problem, object, &
    nodeobj, coefficients, oldvectors, uvec )

    use linear_elastic_globals_m

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

    call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
      maxdisp3D=1 )

    allocate ( phi(1,ndf), u(ncompu*ndf), tmp(ndf,ncompu) )

!   set shape function in the point

    call set_shape_function ( shapefunc, xr, phi )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqdisp], layer=layer )

    tmp = reshape ( u, [ndf,ncompu] )

    uvec = matmul( phi(1,:), tmp )

    deallocate ( phi, u, tmp )

  end subroutine linear_elastic_sample_displacement


! sample value of pressure in one node of the object

  subroutine linear_elastic_sample_pressure ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, p )

    use linear_elastic_globals_m

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

    call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
      maxdisp3D=1 )

    allocate ( psi(1,ndfp), u(ndfp) )

!   set shape function in the point

    if ( any ( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      write(*,'(2a)') ' Error in linear_elastic_sample_pressure:', &
                      ' global interpolation not implemented'
    end if

    call set_shape_function ( shapefuncp, xr, psi )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqpress], layer=layer )

    p = matmul( psi, u )

    deallocate ( psi, u )

  end subroutine linear_elastic_sample_pressure


! sample value of displacement gradients in one node of the object

  subroutine linear_elastic_sample_deriv ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, dvec )

    use linear_elastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: dvec


    integer :: elgrp, elem, comp, j, maxderiv
    real(dp) :: xr(1,mesh%ndim)


    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

    call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
      maxdisp3D=1 )

!   check component = coefficients%i(13)

    if ( coorsys <= 1 .and. disp3D == 0 ) then
      maxderiv = 8
    else
      maxderiv = 11
    end if

    call check ( coefficients, 'linear_elastic_sample_deriv', &
      indexarray=[13], minimum=[1], maximum=[maxderiv] )

!   allocate arrays

    allocate ( detF(1) )
    allocate ( phi(1,ndf), x(nodalp,ndim) )
    allocate ( u(ncompu*ndf), uvector(ndf,ncompu), gradu(1,ncompu,ncompu) )
    allocate ( dphi(1,ndf,ndim), F(1,ndim,ndim) )
    allocate ( Finv(1,ndim,ndim), dphidx(1,ndf,ndim) )
    allocate ( xg(1,ndim), ugvector(1,ncompu) )
    allocate ( dTemp(1) )

!   set shape function in the point

    call set_shape_function ( shapefunc, xr, phi, dphi )

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqdisp], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )
    ugvector = matmul ( phi, uvector )
    xg = matmul ( phi, x )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( disp3D == 1 ) then

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

    if ( coorsys <= 1 .and. disp3D == 0 ) then

!     2D

      if ( comp == 9 .and. lemodel%plane_stress .and. &
           lemodel%thermal_expansion ) then

!       thermal expansion: get temperature change

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(9), value=coefficients%r(4), &
          func=coefficients%func1(1)%p, funcnr=coefficients%i(8), x=x, &
          indx_v=1, layer=layer, phi=phi, indx_e=1, coef=dTemp )

      end if

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
          write(*,'(/2(a/))') 'Error in linear_elastic_deriv:', &
            ' gradu in theta direction not applicable for coorsys=0 '
          stop
        end if

        where ( xg(:,2) < 1e-10_dp )
          dvec = gradu(:,2,2)               ! du_r / dr (r=0)
        elsewhere
          dvec = ugvector(:,2) / xg(:,2)    ! u_r / r
        end where

      case(7)

        dvec = gradu(:,1,1) + gradu(:,2,2) ! dudx + dudy = divergence

        if ( coorsys == 1 ) then ! axisymmetric
          where ( xg(:,2) < 1e-10_dp )
            dvec = dvec + gradu(:,2,2)               ! du_r / dr (r=0)
          elsewhere
            dvec = dvec + ugvector(:,2) / xg(:,2)    ! u_r / r
          end where
        end if

      case(8) ! effective strain

        write(*,'(/2(a/))') 'Error in linear_elastic_deriv:', &
          ' case 8: effective strain not yet available'
        stop

      case(9) ! plane stress epsilon_zz

        if ( lemodel%plane_stress ) then

          dvec = &
              - ( gradu(:,1,1) + gradu(:,2,2) ) * lemodel%nu / (1-lemodel%nu)

          if ( lemodel%thermal_expansion ) then
!           thermal expansion
            dvec = dvec &
                    + lemodel%alpha * dTemp * (1+lemodel%nu) / (1-lemodel%nu)
          end if

        else

          write(*,'(/2(a/))') 'Error in linear_elastic_deriv:', &
            ' case 9: epsilon_zz not available when plane_stress not set.'
          stop

        end if

      case default

        call errormsg_case_default ( 'linear_elastic_sample_deriv', &
          'comp', int_value=comp )

      end select

    else if ( coorsys == 2 .or. disp3D == 1 ) then

!     3D or disp3D = 1

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
      case(11) ! effective strain = sqrt(2II_D)
        write(*,'(/2(a/))') 'Error in linear_elastic_sample_deriv:', &
          ' case 11: effective strain not yet available'
        stop
      case default
        call errormsg_case_default ( 'linear_elastic_sample_deriv', &
          'comp', int_value=comp )
      end select

    end if

    deallocate ( detF )
    deallocate ( phi, x )
    deallocate ( u, uvector, gradu )
    deallocate ( dphi, F )
    deallocate ( Finv, dphidx )
    deallocate ( xg, ugvector )

  end subroutine linear_elastic_sample_deriv


! Normal traction boundary condition for geometries (curves,surfaces).
! For use with add_boundary_elements.

  subroutine linear_elastic_natboun_normal ( mesh, problem, geometry, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use linear_elastic_globals_m

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

      call set_globals_linear_elastic_up_boun ( mesh, coefficients, &
        ndimr=mesh%ndim-1, geometry=geometry, maxdisp3D=1 )

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

        if ( disp3D == 0 ) then
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
      write(*,'(/a/a/)') 'Error in linear_elastic_natboun_normal:', &
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

  end subroutine linear_elastic_natboun_normal


! set global parameters (internal element)

  subroutine set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
    maxdisp3D )

    use linear_elastic_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp
!   the maximum value for disp3D allowed. default=0, i.e. no disp3D available.
    integer, intent(in), optional :: maxdisp3D

    integer :: lmaxdisp3D

    lmaxdisp3D = set_optional ( variable=maxdisp3D, default=0 )

!   check size of coefficients

    call check ( coefficients, 'set_globals_linear_elastic_up', ncoefi=100, &
      ncoefr=50, indexarray=[23,34,37,62,67], minimum=[0,0,0,0,0], &
      maximum=[1,1,2,2,lmaxdisp3D] )

    ndim = mesh%element(elgrp)%ndim
    nodalp = mesh%element(elgrp)%numnod
    physqdisp = coefficients%i(6)
    physqpress = coefficients%i(7)
    layer = coefficients%i(38)
    if ( ndim == 3 ) then
      coorsys = 2
      disp3D = 0
      ncompu = ndim
    else
      coorsys = coefficients%i(23)
      disp3D = coefficients%i(67)
      ncompu = ndim + disp3D
    end if
    nsides = mesh%element(elgrp)%numsides
    nodalpb = mesh%element(elgrp)%sidnumnod
    globalshape = mesh%element(elgrp)%globalshape

!   set number of degrees of freedom displacement

    intpol = coefficients%i(1)

    if ( intpol == 13 .and. any(mesh%element(elgrp)%p(:ndim,2) /= 1 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_linear_elastic_up:', &
        ' For spectral elements GLL nodal distribution is required. '
      stop
    else if ( intpol == 20 .and. &
                            any(mesh%element(elgrp)%p(:ndim,2) /= 0 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_linear_elastic_up:', &
        ' For high-order elements equidistant nodal distribution is required. '
      stop
    end if

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(51)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, 'set_globals_linear_elastic_up', ndf=ndf, &
      ndfb=ndfb )

    if ( coefficients%i(34) == 0 .and. ndf /= nodalp ) then
      write(*,'(/5(a/))') 'Error in set_globals_linear_elastic_up:', &
        ' Number of degrees of freedom of the displacement shape function ', &
        ' (ndf) is different from the number of nodal points (nodalp) ', &
        ' If you really want ndf /= nodalp, set coefficients%i(34) == 1,', &
        ' but you are on your own!! '
      stop
    end if

!   set number of degrees of freedom pressure

    intpolp = coefficients%i(2)

    if ( intpolp > 0 ) then

      shapefuncp%globalshape = globalshape
      shapefuncp%interpolation = intpolp
      shapefuncp%numbering = 'regular'
      shapefuncp%p = coefficients%i(52)
      shapefuncp%spec_eval = 'gauss'
!     needed for spectral quads/hexahedra only to call Pp_at_GLL routine:
      shapefuncp%intrule = shapefunc%p + 1

      call set_ndf ( shapefuncp, 'set_globals_linear_elastic_up', ndf=ndfp )

    else

      ndfp = 0

    end if

!   set integration

    inttype = coefficients%i(40)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_linear_elastic_up:', &
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
      write(*,'(/a/3a/)') 'Error in set_globals_linear_elastic_up:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape
      stop
    end if

    if ( globalshape == 'pyramid' .and. inttype == 2 .and. intrule2 == 0 ) then
      write(*,'(/a/a/3a,i0/)') 'Error in set_globals_linear_elastic_up:', &
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
      write(*,'(/5(a/))') 'Error in set_globals_linear_elastic_up:', &
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_linear_elastic_up


! set global parameters (boundary element)

  subroutine set_globals_linear_elastic_up_boun ( mesh, coefficients, &
    curve, surface, volume, ndimr, geometry, maxdisp3D )

    use linear_elastic_globals_m
    use set_optional_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in), optional :: curve, surface, volume
    integer, intent(in), optional :: ndimr, geometry
!   the maximum value for disp3D allowed. default=0, i.e. no disp3D available.
    integer, intent(in), optional :: maxdisp3D

    integer :: lcurve, lsurface, lvolume
    integer :: lmaxdisp3D

    lmaxdisp3D = set_optional ( variable=maxdisp3D, default=0 )

!   check size of coefficients

    call check ( coefficients, 'set_globals_linear_elastic_up_boun', &
      ncoefi=100, ncoefr=50, indexarray=[23,67], minimum=[0,0], &
      maximum=[1,lmaxdisp3D] )


    physqdisp = coefficients%i(6)
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
          call errormsg_case_default ( 'set_globals_linear_elastic_up_boun', &
            'ndimr', int_value=ndimr )
      end select
    end if

    if ( lcurve > 0 ) then
      ndim = mesh%curves(lcurve)%ndim
      if ( ndim == 3 ) then
        coorsys = 2
        disp3D = 0
        ncompu = ndim
      else
        coorsys = coefficients%i(23)
        disp3D = coefficients%i(67)
        ncompu = ndim + disp3D
      end if
      globalshape = mesh%curves(lcurve)%element%globalshape
      nodalp = mesh%curves(lcurve)%element%numnod
    else if ( lsurface > 0 ) then
      ndim = mesh%surfaces(lsurface)%ndim
      if ( ndim == 3 ) then
        coorsys = 2
        disp3D = 0
        ncompu = ndim
      else
        coorsys = coefficients%i(23)
        disp3D = coefficients%i(67)
        ncompu = ndim + disp3D
      end if
      globalshape = mesh%surfaces(lsurface)%element%globalshape
      nodalp = mesh%surfaces(lsurface)%element%numnod
    else if ( lvolume > 0 ) then
      ndim = mesh%volumes(lvolume)%ndim
      coorsys = 2 ! always 3D
      globalshape = mesh%volumes(lvolume)%element%globalshape
      nodalp = mesh%volumes(lvolume)%element%numnod
    end if

!   set number of degrees of freedom displacement

    intpol = coefficients%i(1)

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(51)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, 'set_globals_linear_elastic_up_boun', ndf=ndf )

    if ( ndf /= nodalp ) then
      write(*,'(/5(a/))') 'Error in set_globals_linear_elastic_up_boun:', &
        ' Number of degrees of freedom of the shape function (ndf) ', &
        ' is different from the number of nodal points in the element', &
        ' (nodalp) ', &
        ' This possibility (ndf /= nodalp) is not available.'
      stop
    end if

!   set number of degrees of freedom pressure

    intpolp = coefficients%i(2)

    if ( intpolp > 0 ) then

      shapefuncp%globalshape = globalshape
      shapefuncp%interpolation = intpolp
      shapefuncp%numbering = 'regular'
      shapefuncp%p = coefficients%i(52)
      shapefuncp%spec_eval = 'gauss'
!     needed for spectral quads only to call Pp_at_GLL routine:
      shapefuncp%intrule = shapefunc%p + 1

      call set_ndf ( shapefuncp, 'set_globals_linear_elastic_up_boun', &
        ndf=ndfp )

    else

      ndfp = 0

    end if

!   set integration

    inttype = coefficients%i(40)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_linear_elastic_up_boun:', &
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
      write(*,'(/a/3a/)') 'Error in set_globals_linear_elastic_up_boun:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape
      stop
    end if

    if ( globalshape == 'pyramid' .and. inttype == 2 .and. intrule2 == 0 ) then
      write(*,'(/a/a/3a,i0/)') 'Error in set_globals_linear_elastic_up_boun:', &
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
      write(*,'(/5(a/))') 'Error in set_globals_linear_elastic_up_boun:', &
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_linear_elastic_up_boun


! set global parameters (boundary Langrangian multiplier)

  subroutine set_globals_linear_elastic_l_boun ( coefficients )

    use linear_elastic_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   Lagrange multiplier interpolation

!   set number of degrees of freedom

    call set_ndf ( shapefunc, &
      name_of_routine='set_globals_linear_elastic_l_boun', &
      ndfl=ndfl, intpoll=intpoll )

    shapefuncl = shapefunc
    shapefuncl%interpolation = intpoll

  end subroutine set_globals_linear_elastic_l_boun


! Preamble for the linear_elastic element

  subroutine set_linear_elastic_elem ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, maxdisp3D )

    use linear_elastic_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    integer, intent(in), optional :: maxdisp3D


    integer :: i, i1(3)


    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

!     set globals

      call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
        maxdisp3D )

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
      if ( ndfp > 0 .and. coefficients%i(62) == 0 ) then
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

        call set_globals_linear_elastic_up ( mesh, coefficients, elgrp, &
          maxdisp3D )

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
      if ( ndfp > 0 .and. coefficients%i(62) == 0 ) then
        call set_shape_function ( shapefuncp, xig, psi )
      end if

    end if

  end subroutine set_linear_elastic_elem


! Unset the preamble for the linear_elastic element
! (deallocate arrays allocated in set_... )

  subroutine unset_linear_elastic_elem ( last, coefficients )

    use linear_elastic_globals_m

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

  end subroutine unset_linear_elastic_elem


! set shapefunction based on scaled global coordinates

  subroutine set_linear_elastic_shape_function_global ( shapefunc, x, xg, &
    phi, scaling )

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

  end subroutine set_linear_elastic_shape_function_global

end module linear_elastic_elements_generic_m

