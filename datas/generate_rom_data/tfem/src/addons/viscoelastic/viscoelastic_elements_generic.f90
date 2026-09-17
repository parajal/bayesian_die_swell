
! Copyright (C) 2005-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! (Generic) Element routines for viscoelastic flow
!

module viscoelastic_elements_generic_m

  use tfem_elem_m
  use stokes_set_globals_m
  use devss_set_globals_m
  use viscoelastic_models_m

  implicit none

contains


! Internal element routine for the right-hand side (div tau) of the momentum
! balance

  subroutine rhs_divtau ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: j, N, ip, m


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, gammap=.true. )

!   coefficients

    if ( first ) then
      cstorage = get_coefficient ( coefficients, index=84, default=1 )
    end if

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of the element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )


!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure ( mesh, problem, elgrp, elem, coefficients, &
        oldvectors%s(1)%p, vemopt )

    end if

!   get conformation tensors

    do m = mode1, mode2

!     single mode conformation
      call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
        oldvectors, cmode=c, mode=m )

      cg(:,:,m) = matmul ( theta, c )

    end do

!   build equations

    if ( matrix ) then

      write(*,'(/a/a/)') 'Error in rhs_divtau: no matrix to build.', &
        'Call build_system with buildmatrix=.false.'
      stop

    end if

    if ( vector ) then

!     stress tensor: use of global variables; result in tauten

      call stress_tensor_viscoelastic ( cg, vemopt=vemopt )

!     - (nabla v)^T:tau

      do ip = 1, ninti
        work6(ip,:,:) = matmul ( dphidx(ip,:,:), tauten(ip,1:ndim,:) )
      end do
      if ( coorsys == 1 .and. vel3D == 0 ) then
        do N = 1, ndf
!         v_r/r * tau_theta,theta
          work6(:,N,2) = work6(:,N,2) + tauvec(:,4) * phi(:,N) / xg(:,2)
        end do
      else if ( coorsys == 1 .and. vel3D == 1 ) then
        do N = 1, ndf
!         v_r/r * tau_theta,theta
          work6(:,N,2) = work6(:,N,2) + tauvec(:,6) * phi(:,N) / xg(:,2)
!         - v_theta/r * tau_theta,r
          work6(:,N,3) = work6(:,N,3) - tauvec(:,5) * phi(:,N) / xg(:,2)
        end do
      end if
      do j = 1, ncompu
        work2(:,j) = - matmul ( detF * wg, work6(:,:,j) )
      end do
      elemvec = reshape ( work2, [ ndf*ncompu ] )

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )

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

      allocate ( work6(ninti,ndf,ncompu), work2(ndf,ncompu) )
      allocate ( c(ndfc,ncomp), cg(ninti,ncomp,nmodes) )
      allocate ( tauvec(ninti,ncompt), tauten(ninti,ncompu,ncompu) )

      if ( vemcompressible ) then
        allocate ( pr(ndfp), press(ninti) )
        call create_vemopt ( vemopt, dep_J=.true., np=ninti )
      else
        call create_vemopt ( vemopt )
      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work6, work2 )
      deallocate ( c, cg )
      deallocate ( tauvec, tauten )

      if ( vemcompressible ) then
        deallocate ( pr, press )
      end if

      call delete ( vemopt )

    end subroutine deallocate_arrays

  end subroutine rhs_divtau


! Boundary element routine for the right-hand side (tau.n) in the momentum
! balance (for curve). Needed for open and internal periodic boundary conditions
! NOTE: only for continuous viscoelastic stresses!

  subroutine rhs_tau_n_curve ( mesh, problem, curve, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: j, N, m, ip


    if ( first ) then

!     first element on this curve

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, curve=curve )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic_boun ( coefficients )

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

!     allocate arrays

      allocate ( wg(ninti), curvel(ninti) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,1), dxdxi(ninti,ndim) )
      allocate ( work4(ninti,ndim), work2(ndf,ndim) )
      allocate ( c(ndfc,ncomp), cg(ninti,ncomp,nmodes) )
      allocate ( tauvec(ninti,ncompt), tauten(ninti,ndim,ndim) )
      allocate ( theta(ninti,ndfc) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      call set_shape_function ( shapefuncc, xig, theta )

      if ( vemcompressible ) then
        allocate ( psi(ninti,ndfp) )
        call set_shape_function ( shapefuncp, xig, psi )
        allocate ( pr(ndfp), press(ninti) )
        call create_vemopt ( vemopt, dep_J=.true., np=ninti )
      else
        call create_vemopt ( vemopt )
      end if

    end if

!   start of the element

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )

    call isoparametric_deformation_curve ( x, dphi(:,:,1), dxdxi, curvel, &
      normal )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x, phi, xg )
      curvel = 2 * pi * xg(:,2) * curvel
    end if

!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure_geometry ( mesh, problem, elem, coefficients, &
        oldvectors%s(1)%p, vemopt, curve=curve )

    end if

!   get conformation tensors

    do m = mode1, mode2

      call get_conformation_geometry ( mesh, oldvectors%p(2)%p, elem, &
        oldvectors, curve=curve, mode=m, cmode=c )

      cg(:,:,m) = matmul ( theta, c )

    end do

!   build equations

    if ( matrix ) then

      write(*,'(/3(a/))') 'Error in rhs_tau_n_curve: no matrix to build.', &
        ' Do not call add_boundary_elements with buildmatrix=.true.', &
        ' or set buildmatrix=.false.'
      stop

    end if

    if ( vector ) then

!     stress tensor: use of global variables; result in tauten

      call stress_tensor_viscoelastic ( cg, vemopt=vemopt )

!     traction force in each integration point: {t}=[tau]{n}

      do ip = 1, ninti
        work4(ip,:) = matmul( tauten(ip,:,:), normal(ip,:) )
      end do

      do j = 1, ndim
        do N = 1, ndf
          work2(N,j) = sum ( work4(:,j) * phi(:,N) * curvel * wg )
        end do
      end do

      elemvec = reshape ( work2, [ ndf*ndim ] )

    end if

    if ( last ) then

!     last element on this curve

      call delete ( vemodel )

      deallocate ( wg, curvel )
      deallocate ( normal )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi, dxdxi )
      deallocate ( work4, work2 )
      deallocate ( c, cg )
      deallocate ( tauvec, tauten )
      deallocate ( theta )

      if ( vemcompressible ) then
        deallocate ( psi )
        deallocate ( pr, press )
      end if

      call delete ( vemopt )

    end if

  end subroutine rhs_tau_n_curve


! Boundary element routine for the right-hand side (tau.n) in the momentum
! balance (for surface). Needed for open and internal periodic boundary
! conditions
! NOTE: only for continuous viscoelastic stresses!

  subroutine rhs_tau_n_surface ( mesh, problem, surface, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: surface, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: j, N, m, ip


    if ( first ) then

!     first element on this surface

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, surface=surface )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic_boun ( coefficients )

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,ndim-1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim-1) )
      allocate ( work4(ninti,ndim), work2(ndf,ndim) )
      allocate ( c(ndfc,ncomp), cg(ninti,ncomp,nmodes) )
      allocate ( tauvec(ninti,ncompt), tauten(ninti,ndim,ndim) )
      allocate ( theta(ninti,ndfc) )
      allocate ( dxdxis(ninti,ndim,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      call set_shape_function ( shapefuncc, xig, theta )

      if ( vemcompressible ) then
        allocate ( psi(ninti,ndfp) )
        call set_shape_function ( shapefuncp, xig, psi )
        allocate ( pr(ndfp), press(ninti) )
        call create_vemopt ( vemopt, dep_J=.true., np=ninti )
      else
        call create_vemopt ( vemopt )
      end if

    end if

!   start of the element

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, normal )

!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure_geometry ( mesh, problem, elem, coefficients, &
        oldvectors%s(1)%p, vemopt, surface=surface )

    end if

!   get conformation tensors

    do m = mode1, mode2

      call get_conformation_geometry ( mesh, oldvectors%p(2)%p, elem, &
        oldvectors, surface=surface, mode=m, cmode=c )

      cg(:,:,m) = matmul ( theta, c )

    end do

!   build equations

    if ( matrix ) then

      write(*,'(/3(a/))') 'Error in rhs_tau_n_surface: no matrix to build.', &
        ' Do not call add_boundary_elements with buildmatrix=.true.', &
        ' or set buildmatrix=.false.'
      stop

    end if

    if ( vector ) then

!     stress tensor: use of global variables; result in tauten

      call stress_tensor_viscoelastic ( cg, vemopt=vemopt )

!     traction force in each integration point: {t}=[tau]{n}

      do ip = 1, ninti
        work4(ip,:) = matmul( tauten(ip,:,:), normal(ip,:) )
      end do

      do j = 1, ndim
        do N = 1, ndf
          work2(N,j) = sum ( work4(:,j) * phi(:,N) * surfl * wg )
        end do
      end do

      elemvec = reshape ( work2, [ ndf*ndim ] )

    end if

    if ( last ) then

!     last element on this surface

      call delete ( vemodel )

      deallocate ( wg, surfl )
      deallocate ( normal )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi )
      deallocate ( work4, work2 )
      deallocate ( c, cg )
      deallocate ( tauvec, tauten )
      deallocate ( theta )
      deallocate ( dxdxis )

      if ( vemcompressible ) then
        deallocate ( psi )
        deallocate ( pr, press )
      end if

      call delete ( vemopt )

    end if

  end subroutine rhs_tau_n_surface


! Internal element routine for implicit terms of the CE in the momentum
! balance (convection implicit). For first-order and second-order integration.

  subroutine divtau_implicit_ce_elem_c ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i5, i6, i7, ip, md
    integer :: N, M, i, j, id
    real(dp) :: deltat

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_devssg_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1 )

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   first element in this group

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'divtau_implicit_ce_elem_c', &
        indexarray=[48], minimum=[0], maximum=[2] )

      if ( varpar ) then
!       check coefficients for variable coefficients
        call check ( coefficients, 'divtau_implicit_ce_elem_c', &
          indexarray=[68,69,70,81], minimum=[0,0,0,0], maximum=[3,3,3,3] )
      end if

      select case ( vemodel%model )
        case(4,8:12,14,17,20:21,23:)
          write(*,'(/a/a,i0/)') 'Error in divtau_implicit_ce_elem_c:', &
            ' Model has not yet been implemented: model = ', vemodel%model
          stop
        case default
!         allowed models
      end select

      if ( size(elemmat,2) /= ncompu*ndf ) then
        write(*,'(/3(a/))') 'Error in divtau_implicit_ce_elem_c:', &
          ' the degrees of freedom in physqcol ', &
          ' must be (u), i.e. the velocity only'
        stop
      end if

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

    end if

!   start of the element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    call shape_derivative ( dtheta, Finv, dthetadx )


!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure ( mesh, problem, elgrp, elem, coefficients, &
        oldvectors%s(1)%p, vemoptn )

    end if


!   get conformation tensor at previous time step

    if ( coefficients%i(48) == 2 ) then

!     temporary ALE

      if ( logc .and. .not. exps ) then
        write(*,'(/3(a/))') 'Error get_cn_tALE:', &
          ' for tALE a log parameter coefficients%i(21) = 1 must be ', &
          ' combined with L2 projection (coefficients%i(49) = 1). '
        stop
      else if ( bten .and. .not. cproj ) then
        write(*,'(/3(a/))') 'Error get_cn_tALE:', &
          ' for tALE a btensor parameter coefficients%i(71) >= 1 must be ', &
          ' combined with L2 projection (coefficients%i(72) = 1). '
        stop
      end if

!     compute reference and real coordinates in the temporary ALE meshes
!     NOTE: groups not yet used.

      select case ( coefficients%i(22) )

      case (1,3)

!       first-order

        call set_globals_stokes_tALE ( oldvectors%m(1)%p, oldvectors%m(2)%p, &
          deform_n=.true. )

      case (2,4:7)

!       second-order
        call set_globals_stokes_tALE ( oldvectors%m(1)%p, oldvectors%m(2)%p, &
          oldvectors%m(3)%p, oldvectors%m(4)%p, deform_n=.true. )

      case default

        call errormsg_case_default ( 'divtau_implicit_ce_elem_c', &
          'coefficients%i(22)', int_value=coefficients%i(22) )

      end select

!     shape functions for the conformation tensor

      call set_shape_function ( shapefuncc, xig_n, theta_n, dtheta_n )
      call shape_derivative ( dtheta_n, Finv_n, dthetadx_n )

!     get old value of c (cn)

      call get_cn_tALE

    else

!     standard case

!     get old value of c (cn)

      if ( c_direct ) then
        call get_cn_standard2
      else
        call get_cn_standard1
      end if

    end if

!   set time step

    deltat = coefficients%r(8)

!   get mesh velocity at the current time

    if ( coefficients%i(48) == 1 ) then

!     standard ALE

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, u )

      tmp = reshape ( u, [ndf,ndim] )

      uvecmeshnp1 = matmul ( phi, tmp )

    else if ( coefficients%i(48) == 2 ) then

!     temporary ALE

      select case ( coefficients%i(22) )
      case (1,3)
!       first-order
        uvecmeshnp1 = ( xg - xg_n ) / deltat
      case (2,4:7)
!       second-order
        uvecmeshnp1 = ( 1.5_dp*xg - 2.0_dp*xg_n + 0.5_dp*xg_nm1 ) / deltat
      case default
        call errormsg_case_default ( 'divtau_implicit_ce_elem_c', &
          'coefficients%i(22)', int_value=coefficients%i(22) )
      end select

    end if

    if ( varpar ) then

!     determine the variable material parameters in mvemodel

      call evaluate_mve ( mesh, problem, elem, elgrp, coefficients, &
        oldvectors, modulus=.true., lambda=.true., nonlin=.true., alam=.true. )

!     copy arrays from mvemodel for easier handling

      do ip = 1, ninti
        Gmod_g(ip,:) = mvemodel(ip)%modulus(mode1:mode2)
        lambda_g(ip,:) = mvemodel(ip)%lambda(mode1:mode2)
        nonlin_g(ip,:,:) = mvemodel(ip)%nonlin(:,mode1:mode2)
      end do

    else

!     fixed coefficients

      do ip = 1, ninti
        Gmod_g(ip,:) = vemodel%modulus(mode1:mode2)
        lambda_g(ip,:) = vemodel%lambda(mode1:mode2)
        nonlin_g(ip,:,:) = vemodel%nonlin(:,mode1:mode2)
      end do

    end if

!   lambda modified G modulus

    if ( coverl ) then
!     c/lambda added to both sides of the constitutive equation
      if ( coefficients%i(58) > 0 ) then ! add c/non-linear parameter
        lambda_small_g = nonlin_g(:,coefficients%i(58),:)
      else ! add c/lambda
        lambda_small_g = lambda_g
      end if
      Gmodl_g = Gmod_g / ( 1 + deltat / lambda_small_g )
    else
      Gmodl_g = Gmod_g
    end if

    if ( matrix ) then

!     set c in tensor format

      if ( coorsys <= 1 .and. vel3D == 0 ) then

        do ip = 1, ninti
          cten(ip,1,1) = dot_product ( cng(ip,1,mode1:mode2), Gmodl_g(ip,:) )
          cten(ip,1,2) = dot_product ( cng(ip,2,mode1:mode2), Gmodl_g(ip,:) )
          cten(ip,2,1) = cten(ip,1,2)
          cten(ip,2,2) = dot_product ( cng(ip,3,mode1:mode2), Gmodl_g(ip,:) )
        end do

      else if ( coorsys == 2 .or. vel3D == 1 ) then

        do ip = 1, ninti
          cten(ip,1,1) = dot_product ( cng(ip,1,mode1:mode2), Gmodl_g(ip,:) )
          cten(ip,1,2) = dot_product ( cng(ip,2,mode1:mode2), Gmodl_g(ip,:) )
          cten(ip,1,3) = dot_product ( cng(ip,3,mode1:mode2), Gmodl_g(ip,:) )
          cten(ip,2,1) = cten(ip,1,2)
          cten(ip,2,2) = dot_product ( cng(ip,4,mode1:mode2), Gmodl_g(ip,:) )
          cten(ip,2,3) = dot_product ( cng(ip,5,mode1:mode2), Gmodl_g(ip,:) )
          cten(ip,3,1) = cten(ip,1,3)
          cten(ip,3,2) = cten(ip,2,3)
          cten(ip,3,3) = dot_product ( cng(ip,6,mode1:mode2), Gmodl_g(ip,:) )
        end do

      end if

!     set gradc in tensor format

      if ( coorsys <= 1 .and. vel3D == 0 ) then

        do j = 1, ndim
          do ip = 1, ninti
            gradcten(ip,1,1,j) = &
                      dot_product ( gradcn(ip,1,j,mode1:mode2), Gmodl_g(ip,:) )
            gradcten(ip,1,2,j) = &
                      dot_product ( gradcn(ip,2,j,mode1:mode2), Gmodl_g(ip,:) )
            gradcten(ip,2,1,j) = gradcten(ip,1,2,j)
            gradcten(ip,2,2,j) = &
                      dot_product ( gradcn(ip,3,j,mode1:mode2), Gmodl_g(ip,:) )
          end do
        end do

      else if ( coorsys == 2 .or. vel3D == 1 ) then

        do j = 1, ndim
          do ip = 1, ninti
            gradcten(ip,1,1,j) = &
                      dot_product ( gradcn(ip,1,j,mode1:mode2), Gmodl_g(ip,:) )
            gradcten(ip,1,2,j) = &
                      dot_product ( gradcn(ip,2,j,mode1:mode2), Gmodl_g(ip,:) )
            gradcten(ip,1,3,j) = &
                      dot_product ( gradcn(ip,3,j,mode1:mode2), Gmodl_g(ip,:) )
            gradcten(ip,2,1,j) = gradcten(ip,1,2,j)
            gradcten(ip,2,2,j) = &
                      dot_product ( gradcn(ip,4,j,mode1:mode2), Gmodl_g(ip,:) )
            gradcten(ip,2,3,j) = &
                      dot_product ( gradcn(ip,5,j,mode1:mode2), Gmodl_g(ip,:) )
            gradcten(ip,3,1,j) = gradcten(ip,1,3,j)
            gradcten(ip,3,2,j) = gradcten(ip,2,3,j)
            gradcten(ip,3,3,j) = &
                      dot_product ( gradcn(ip,6,j,mode1:mode2), Gmodl_g(ip,:) )
          end do
        end do

      end if

      do M = 1, ndf

!       dphidx_m^M c_mi

        do ip = 1, ninti
          work2(ip,:) = matmul ( dphidx(ip,M,:), cten(ip,1:ndim,:) )
        end do

        do N = 1, ndf

!         dphidx_k^N dphidx_m^M c_km

          do ip = 1, ninti
            work1(ip) = dot_product ( dphidx(ip,N,:), work2(ip,1:ndim) )
          end do

!         dphidx_k^N dc_ki/dx_j

          do ip = 1, ninti
            do id = 1, ndim
              work6(ip,:,id) = &
                          matmul ( dphidx(ip,N,:), gradcten(ip,1:ndim,:,id) )
            end do
          end do

          do i = 1, ncompu
            do j = 1, ndim

!             dphidx_k^N * dphidx_m^N * c_mi * delta_jk

              work = dphidx(:,N,j) * work2(:,i)

!             dphidx_k^N * dphidx_m^N * c_km * delta_ij

              if ( i == j ) then
                work = work + work1
              end if

!             - phi_M dphidx_k^N dc_ki/dx_j

              work = work - phi(:,M) * work6(:,i,j)

              work5(N,M,i,j) = sum ( work * detF * wg )

            end do
          end do

          if ( vel3D == 1 ) then
            work5(N,M,1:2,3) = 0
            work5(N,M,3,3) = sum ( work1 * detF * wg )
          end if

        end do

      end do

      work5 = deltat * work5


      if ( coorsys == 1 .and. vel3D == 0 ) then

!       axisymmetric

        do ip = 1, ninti
          work1(ip) = 2 * dot_product ( cng(ip,4,mode1:mode2), &
                                          Gmodl_g(ip,:) ) / xg(ip,2) ** 2
        end do

        do N = 1, ndf
          do M = 1, ndf
!           phi^N 2c_tt/r^2 phi^M
            work4(N,M) = sum ( phi(:,N) * work1 * phi(:,M) * detF * wg )
          end do
        end do

        work5(:,:,2,2) = work5(:,:,2,2) + deltat * work4

        do i = 1, ndim
          do ip = 1, ninti
            work11(ip,i) = dot_product ( gradcn(ip,4,i,mode1:mode2), &
                                                Gmodl_g(ip,:) ) / xg(ip,2)
         end do
       end do

        do N = 1, ndf
          do M = 1, ndf
            do j = 1, ndim
!             - phi^N phi^M dc_tt/dx_j 1/r
              work7(N,M,j) = &
                - sum ( phi(:,N) * work11(:,j) * phi(:,M) * detF * wg )
            end do
          end do
        end do

        work5(:,:,2,:) = work5(:,:,2,:) + deltat * work7

      else if ( coorsys == 1 .and. vel3D == 1 ) then

!       axisymmetric with 3D velocities

!       additional terms from (nabla v)^T:(-u.nabla c + L.c + c.L^T)

!       z-theta (all terms lead to zero after summation)

!       r-theta (note: some terms lead to zero after summation)

        do N = 1, ndf
          do M = 1, ndf
            work4(N,M) = &
                2 * sum ( phi(:,N) * cten(:,1,3) * dphidx(:,M,1) / xg(:,2) &
                            * detF * wg ) &
              + 2 * sum ( phi(:,N) * cten(:,2,3) * dphidx(:,M,2) / xg(:,2) &
                            * detF * wg ) &
              - 2 * sum ( phi(:,N) * cten(:,2,3) * phi(:,M) / xg(:,2) ** 2 &
                            * detF * wg )
          end do
        end do

        work5(:,:,2,3) = work5(:,:,2,3) + deltat * work4

!       theta-r

        do N = 1, ndf
          do M = 1, ndf
            work4(N,M) = &
                    sum ( dphidx(:,N,1) * cten(:,1,3) * phi(:,M) / xg(:,2) &
                            * detF * wg ) &
                  + sum ( dphidx(:,N,2) * cten(:,2,3) * phi(:,M) / xg(:,2) &
                            * detF * wg ) &
                  - sum ( phi(:,N) * cten(:,2,3) * phi(:,M) / xg(:,2) ** 2 &
                            * detF * wg ) &
                  - sum ( phi(:,N) * cten(:,1,3) * dphidx(:,M,1) / xg(:,2) &
                            * detF * wg ) &
                  - sum ( phi(:,N) * cten(:,2,3) * dphidx(:,M,2) / xg(:,2) &
                            * detF * wg )
          end do
        end do

        work5(:,:,3,2) = work5(:,:,3,2) + deltat * work4

!       theta-theta

        do N = 1, ndf
          do M = 1, ndf
            work4(N,M) = &
                  - sum ( dphidx(:,N,2) * cten(:,3,3) * phi(:,M) / xg(:,2) &
                            * detF * wg ) &
                  + sum ( phi(:,N) * cten(:,3,3) * phi(:,M) / xg(:,2) ** 2 &
                            * detF * wg ) &
                  - sum ( phi(:,N) * cten(:,1,2) * dphidx(:,M,1) / xg(:,2) &
                            * detF * wg ) &
                  - sum ( phi(:,N) * cten(:,2,2) * dphidx(:,M,2) / xg(:,2) &
                            * detF * wg ) &
                  + sum ( phi(:,N) * ( cten(:,2,2) - cten(:,3,3) ) &
                           * phi(:,M) / xg(:,2) ** 2 * detF * wg ) &
                  - sum ( dphidx(:,N,1) * cten(:,1,2) * phi(:,M) / xg(:,2) &
                            * detF * wg ) &
                  - sum ( dphidx(:,N,2) * ( cten(:,2,2) - cten(:,3,3) ) &
                           * phi(:,M) / xg(:,2) * detF * wg )
          end do
        end do

        work5(:,:,3,3) = work5(:,:,3,3) + deltat * work4

!       r-r

        do N = 1, ndf
          do M = 1, ndf
!           phi^N 2c_tt/r^2 phi^M
            work4(N,M) = &
                2 * sum ( phi(:,N) * cten(:,3,3) * phi(:,M) / xg(:,2) ** 2 &
                            * detF * wg )
          end do
        end do

        work5(:,:,2,2) = work5(:,:,2,2) + deltat * work4

!       + phi^N phi^M dc_rt/dx_j 1/r

        do N = 1, ndf
          do M = 1, ndf
            do j = 1, ndim
              work7(N,M,j) = &
                   sum ( phi(:,N) * gradcten(:,2,3,j) * phi(:,M) / xg(:,2) &
                                * detF * wg )
            end do
          end do
        end do

        work5(:,:,3,1:ndim) = work5(:,:,3,1:ndim) + deltat * work7

!       - phi^N phi^M dc_tt/dx_j 1/r

        do N = 1, ndf
          do M = 1, ndf
            do j = 1, ndim
              work7(N,M,j) = &
                 - sum ( phi(:,N) * gradcten(:,3,3,j) * phi(:,M) / xg(:,2) &
                                * detF * wg )
            end do
          end do
        end do

        work5(:,:,2,1:ndim) = work5(:,:,2,1:ndim) + deltat * work7

      end if

!     put into element matrix

      i5 = ndf    ! u
      i6 = 2*ndf  ! v
      i7 = 3*ndf  ! w

      if ( coorsys <= 1. .and. vel3D == 0 ) then

        elemmat( 1:i5, 1:i5 )       = work5(:,:,1,1)
        elemmat( 1:i5, i5+1:i6 )    = work5(:,:,1,2)
        elemmat( i5+1:i6, 1:i5 )    = work5(:,:,2,1)
        elemmat( i5+1:i6, i5+1:i6 ) = work5(:,:,2,2)

      else if ( coorsys == 2 .or. vel3D == 1 ) then

        elemmat( 1:i5, 1:i5 )       = work5(:,:,1,1)
        elemmat( 1:i5, i5+1:i6 )    = work5(:,:,1,2)
        elemmat( 1:i5, i6+1:i7 )    = work5(:,:,1,3)
        elemmat( i5+1:i6, 1:i5 )    = work5(:,:,2,1)
        elemmat( i5+1:i6, i5+1:i6 ) = work5(:,:,2,2)
        elemmat( i5+1:i6, i6+1:i7 ) = work5(:,:,2,3)
        elemmat( i6+1:i7, 1:i5 )    = work5(:,:,3,1)
        elemmat( i6+1:i7, i5+1:i6 ) = work5(:,:,3,2)
        elemmat( i6+1:i7, i6+1:i7 ) = work5(:,:,3,3)

      end if

    end if

    if ( vector ) then

!     - (nabla v)^T: ( G h )

!     relaxation term of CE only

      if ( varpar ) then

!       variable parameters

        if ( coorsys <= 1 .and. vel3D == 0 ) then
          call rhs_viscoelastic_relax_2D ( vemodel, cng, rhsd, mode1=mode1, &
            mode2=mode2, mvemodel=mvemodel, vemopt=vemoptn )
        else if ( coorsys == 2 .or. vel3D == 1 ) then
          call rhs_viscoelastic_relax_3D ( vemodel, cng, rhsd, mode1=mode1, &
            mode2=mode2, mvemodel=mvemodel, vemopt=vemoptn )
        end if

      else

!       fixed parameters

        if ( coorsys <= 1 .and. vel3D == 0 ) then
          call rhs_viscoelastic_relax_2D ( vemodel, cng, rhsd, mode1=mode1, &
            mode2=mode2, vemopt=vemoptn )
        else if ( coorsys == 2 .or. vel3D == 1 ) then
          call rhs_viscoelastic_relax_3D ( vemodel, cng, rhsd, mode1=mode1, &
            mode2=mode2, vemopt=vemoptn )
        end if

      end if

!     add mesh velocity term to the right-hand side

      if ( any( coefficients%i(48) == [1,2] ) ) then

        do md = 1, nmodes
          do ip = 1, ninti
            rhsd(ip,:,md) = rhsd(ip,:,md) &
                             + matmul ( gradcn(ip,:,:,md), uvecmeshnp1(ip,:) )
          end do
        end do

      end if


      if ( coorsys <= 1 .and. vel3D == 0 ) then

        do ip = 1, ninti
          hten(ip,1,1) = &
            dot_product ( cng(ip,1,mode1:mode2) - 1, Gmod_g(ip,:) ) + &
            deltat * dot_product ( rhsd(ip,1,mode1:mode2),Gmodl_g(ip,:) )
          hten(ip,1,2) = dot_product ( cng(ip,2,mode1:mode2), Gmod_g(ip,:) ) + &
            deltat * dot_product ( rhsd(ip,2,mode1:mode2),Gmodl_g(ip,:) )
          hten(ip,2,1) = hten(ip,1,2)
          hten(ip,2,2) = &
            dot_product ( cng(ip,3,mode1:mode2) - 1, Gmod_g(ip,:) ) + &
            deltat * dot_product ( rhsd(ip,3,mode1:mode2),Gmodl_g(ip,:) )
        end do

      else if ( coorsys == 2 .or. vel3D == 1 ) then

        do ip = 1, ninti
          hten(ip,1,1) = &
            dot_product ( cng(ip,1,mode1:mode2) - 1, Gmod_g(ip,:) ) + &
            deltat * dot_product ( rhsd(ip,1,mode1:mode2), Gmodl_g(ip,:) )
          hten(ip,1,2) = dot_product ( cng(ip,2,mode1:mode2), Gmod_g(ip,:) ) &
            + deltat * dot_product ( rhsd(ip,2,mode1:mode2), Gmodl_g(ip,:) )
          hten(ip,1,3) = dot_product ( cng(ip,3,mode1:mode2), Gmod_g(ip,:) ) &
            + deltat * dot_product ( rhsd(ip,3,mode1:mode2), Gmodl_g(ip,:) )
          hten(ip,2,1) = hten(ip,1,2)
          hten(ip,2,2) = &
            dot_product ( cng(ip,4,mode1:mode2) - 1, Gmod_g(ip,:) ) + &
            deltat * dot_product ( rhsd(ip,4,mode1:mode2), Gmodl_g(ip,:) )
          hten(ip,2,3) = dot_product ( cng(ip,5,mode1:mode2), Gmod_g(ip,:) ) &
            + deltat * dot_product ( rhsd(ip,5,mode1:mode2), Gmodl_g(ip,:) )
          hten(ip,3,1) = hten(ip,1,3)
          hten(ip,3,2) = hten(ip,2,3)
          hten(ip,3,3) = &
            dot_product ( cng(ip,6,mode1:mode2) - 1, Gmod_g(ip,:) ) + &
            deltat * dot_product ( rhsd(ip,6,mode1:mode2), Gmodl_g(ip,:) )
        end do

      end if

      if ( coorsys == 1 .and. vel3D == 0 ) then
!       work1 = hten_tt
        do ip = 1, ninti
          work1(ip) = dot_product ( cng(ip,4,mode1:mode2) - 1, Gmod_g(ip,:) ) &
                + deltat * dot_product ( rhsd(ip,4,mode1:mode2), Gmodl_g(ip,:) )
        end do
      end if

!     - (nabla v)^T:tau

      do ip = 1, ninti
        work9(ip,:,:) = matmul ( dphidx(ip,:,:), hten(ip,1:ndim,:) )
      end do
      if ( coorsys == 1 .and. vel3D == 0 ) then
        do N = 1, ndf
          work9(:,N,2) = work9(:,N,2) + work1 * phi(:,N) / xg(:,2)
        end do
      else if ( coorsys == 1 .and. vel3D == 1 ) then
        do N = 1, ndf
!         v_r/r * h_theta,theta
          work9(:,N,2) = work9(:,N,2) + hten(:,3,3) * phi(:,N) / xg(:,2)
!         - v_theta/r * h_theta,r
          work9(:,N,3) = work9(:,N,3) - hten(:,3,2) * phi(:,N) / xg(:,2)
        end do
      end if
      do j = 1, ncompu
        work10(:,j) = - matmul ( detF * wg, work9(:,:,j) )
      end do
      elemvec = reshape ( work10, [ ndf*ncompu ] )

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_devssg_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )


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

      allocate ( work(ninti), work2(ninti,ncompu), work11(ninti,ndim) )
      allocate ( tmp(ndf,ndim) )
      allocate ( work4(ndf,ndf), work1(ninti) )
      allocate ( work5(ndf,ndf,ncompu,ncompu) )
      allocate ( work6(ninti,ncompu,ndim) )
      allocate ( work7(ndf,ndf,ndim) )
      allocate ( work9(ninti,ndf,ncompu), work10(ndf,ncompu) )
      allocate ( gradcn(ninti,ncompc,ndim,nmodes) )

      allocate ( u(ndf*ndim), uvecmeshnp1(ninti,ndim) )

      allocate ( c(ndfc,ncomp), cn(ndfc,ncompc), cng(ninti,ncompc,nmodes) )
      allocate ( sg(ninti,ncompc) )
      allocate ( cten(ninti,ncompu,ncompu) )
      allocate ( gradcten(ninti,ncompu,ncompu,ndim) )
      allocate ( hten(ninti,ncompu,ncompu) )
      allocate ( rhsd(ninti,ncompc,nmodes) )
      allocate ( Gmod_g(ninti,mode2-mode1+1) )
      allocate ( Gmodl_g(ninti,mode2-mode1+1) )
      allocate ( lambda_g(ninti,mode2-mode1+1) )
      allocate ( nonlin_g(ninti,size(vemodel%nonlin,1),mode2-mode1+1) )

      if ( coefficients%i(48) == 2 ) then

!       temporary ALE scheme; allocate additional arrays

        allocate ( xg_n(ninti,ndim), grpelm_n(ninti,2) )
        allocate ( xig_n(ninti,ndim), phi_n(ninti,ndf) )
        allocate ( xg_nm1(ninti,ndim), grpelm_nm1(ninti,2) )
        allocate ( xig_nm1(ninti,ndim), phi_nm1(ninti,ndf) )

        allocate ( dphi_n(ninti,ndf,ndim), F_n(ninti,ndim,ndim) )
        allocate ( Finv_n(ninti,ndim,ndim), detF_n(ninti) )
        allocate ( theta_n(ninti,ndfc), dtheta_n(ninti,ndfc,ndim) )
        allocate ( dthetadx_n(ninti,ndfc,ndim) )

      end if

      if ( coverl ) then
        allocate ( lambda_small_g(ninti,mode2-mode1+1) )
      end if

      if ( varpar ) then
        allocate ( mvemodel(ninti) )
        mvemodel = vemodel
      end if

      if ( vemcompressible ) then
        allocate ( pr(ndfp), press(ninti) )
        call create_vemopt ( vemoptn, dep_J=.true., np=ninti )
      else
        call create_vemopt ( vemoptn )
      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      integer :: ip

      deallocate ( work, work2, work11 )
      deallocate ( tmp )
      deallocate ( work4, work1 )
      deallocate ( work5 )
      deallocate ( work6 )
      deallocate ( work7 )
      deallocate ( work9, work10 )
      deallocate ( gradcn )

      deallocate ( u, uvecmeshnp1 )

      deallocate ( c, cn, cng )
      deallocate ( sg )
      deallocate ( cten )
      deallocate ( gradcten )
      deallocate ( hten )
      deallocate ( rhsd )
      deallocate ( Gmod_g )
      deallocate ( Gmodl_g )
      deallocate ( lambda_g )
      deallocate ( nonlin_g )

      if ( coefficients%i(48) == 2 ) then

!       temporary ALE scheme; deallocate additional arrays

        deallocate ( xg_n, grpelm_n )
        deallocate ( xig_n, phi_n )
        deallocate ( xg_nm1, grpelm_nm1 )
        deallocate ( xig_nm1, phi_nm1 )

        deallocate ( dphi_n, F_n, Finv_n, detF_n )
        deallocate ( theta_n, dtheta_n )
        deallocate ( dthetadx_n )

      end if

      if ( coverl ) then
        deallocate ( lambda_small_g )
      end if

      if ( varpar ) then
        do ip = 1, ninti
          call delete ( mvemodel(ip) )
        end do
        deallocate ( mvemodel )
      end if

      if ( vemcompressible ) then
        deallocate ( pr, press )
      end if

      call delete ( vemoptn )

    end subroutine deallocate_arrays


!   Get cn for the standard case (ALE or Euler)
!   This routine computes both cn and gradcn in the integration points
!   from nodal values.

    subroutine get_cn_standard1

      integer :: md, i

      do md = mode1, mode2

        if ( exps .or. cproj ) then

!         projection of c=exp(s) or c=b.b^T has been performed in
!         a separate problem

!         single mode conformation
          call get_conformation ( mesh, oldvectors%p(3)%p, elgrp, elem, &
            oldvectors, cmode=cn, mode=md, isv=3, cst=1 )

        else

!         standard case

!         single mode conformation
          call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
            oldvectors, cmode=c, mode=md )

          if ( logc ) then
!           transform sn to cn
!           NOTE: we are transforming the nodal values sn -> cn directly and
!           define cn in the element with the same interpolation as sn. Maybe it
!           is better to do a global least-square projection of cn = exp(sn) to
!           define the nodal point values of cn. See exps=.true.
            if ( coorsys <= 1 .and. vel3D == 0 ) then
              call conformation_2D_log ( vemodel, c, cn )
            else if ( coorsys == 2 .or. vel3D == 1 ) then
              call conformation_3D_log ( vemodel, c, cn )
            end if
          else if ( bten ) then
!           transform bn to cn
!           NOTE: we are transforming the nodal values bn -> cn directly and
!           define cn in the element with the same interpolation as bn. Maybe it
!           is better to do a global least-square projection of cn = bn.bn^T to
!           define the nodal point values of cn. See cproj=.true.
            if ( coorsys <= 1 .and. vel3D == 0 ) then
              call conformation_2D_b ( c, cn, vemodel%bvariant )
            else if ( coorsys == 2 .or. vel3D == 1 ) then
              call conformation_3D_b ( c, cn, vemodel%bvariant )
            end if
          else
            cn = c
          end if

        end if

!       c in Gauss points

        cng(:,:,md) = matmul ( theta, cn )

!       grad cn term

        do i = 1, ndim
          gradcn(:,:,i,md) = matmul ( dthetadx(:,:,i), cn )
        end do

      end do

    end subroutine get_cn_standard1


!   get cn for the standard case (ALE or Euler)
!   This routine computes cn directly in the Gauss integration points
!   and gradcn from nodal values.

    subroutine get_cn_standard2

      integer :: md, i

      do md = mode1, mode2

!       get s, b or c in the nodes

!       single mode conformation
        call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
          oldvectors, cmode=c, mode=md )

!       s, b or c in Gauss points

        sg = matmul ( theta, c )

!       c in Gauss points directly from exp(s), b.b^T or copy

        if ( logc ) then
!         transform s to c
          if ( coorsys <= 1 .and. vel3D == 0 ) then
            call conformation_2D_log ( vemodel, sg, cng(:,:,md) )
          else if ( coorsys == 2 .or. vel3D == 1 ) then
            call conformation_3D_log ( vemodel, sg, cng(:,:,md) )
          end if
        else if ( bten ) then
!         transform b to c
          if ( coorsys <= 1 .and. vel3D == 0 ) then
            call conformation_2D_b ( sg, cng(:,:,md), vemodel%bvariant )
          else if ( coorsys == 2 .or. vel3D == 1 ) then
            call conformation_3D_b ( sg, cng(:,:,md), vemodel%bvariant )
          end if
        else
!         copy
          cng(:,:,md) = sg
        end if

!       gradc

!       first compute c in the nodes

        if ( exps .or. cproj ) then

!         projection of c=exp(s) or c=b.b^T has been performed in
!         a separate problem

!         single mode conformation
          call get_conformation ( mesh, oldvectors%p(3)%p, elgrp, elem, &
            oldvectors, cmode=cn, mode=md, isv=3, cst=1 )

        else

!         direct nodal method (no projection)

          if ( logc ) then
!           transform sn to cn
!           NOTE: we are transforming the nodal values sn -> cn directly and
!           define cn in the element with the same interpolation as sn. Maybe it
!           is better to do a global least-square projection of cn = exp(sn) to
!           define the nodal point values of cn. See exps=.true.
            if ( coorsys <= 1 .and. vel3D == 0 ) then
              call conformation_2D_log ( vemodel, c, cn )
            else if ( coorsys == 2 .or. vel3D == 1 ) then
              call conformation_3D_log ( vemodel, c, cn )
            end if
          else if ( bten ) then
!           transform bn to cn
!           NOTE: we are transforming the nodal values bn -> cn directly and
!           define cn in the element with the same interpolation as bn. Maybe it
!           is better to do a global least-square projection of cn = bn.bn^T to
!           define the nodal point values of cn. See cproj=.true.
            if ( coorsys <= 1 .and. vel3D == 0 ) then
              call conformation_2D_b ( c, cn, vemodel%bvariant )
            else if ( coorsys == 2 .or. vel3D == 1 ) then
              call conformation_3D_b ( c, cn, vemodel%bvariant )
            end if
          else
            cn = c
          end if

        end if

!       grad cn term

        do i = 1, ndim
          gradcn(:,:,i,md) = matmul ( dthetadx(:,:,i), cn )
        end do

      end do

    end subroutine get_cn_standard2


!   get cn for the temporary ALE case

    subroutine get_cn_tALE

      integer :: ip, md, i

      do ip = 1, ninti

        do md = mode1, mode2

          if ( exps .or. cproj ) then

!           projection of c=exp(s) or c=b.b^T has been performed in
!           a separate problem

!           single mode conformation
            call get_conformation ( oldvectors%m(2)%p, oldvectors%p(3)%p, &
              grpelm_n(ip,1), grpelm_n(ip,2), oldvectors, cmode=cn, mode=md, &
              isv=3, cst=1 )

          else

!           standard case

!           single mode conformation
            call get_conformation ( oldvectors%m(2)%p, oldvectors%p(2)%p, &
              grpelm_n(ip,1), grpelm_n(ip,2), oldvectors, cmode=cn, mode=md )

          end if

!         c in Gauss points

          cng(ip,:,md) = matmul ( theta_n(ip,:), cn )

!         grad cn term

          do i = 1, ndim
            gradcn(ip,:,i,md) = matmul ( dthetadx_n(ip,:,i), cn )
          end do

        end do

      end do

    end subroutine get_cn_tALE

  end subroutine divtau_implicit_ce_elem_c


! Element routine on an object for open boundary conditions in the implicit
! stress formulation

  subroutine implicit_stress_open_boundary ( mesh, problem, eleminfo, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(eleminfo_t), intent(in) :: eleminfo
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: elgrp, elem, object, elemo
    integer :: ip, orientation
    integer :: i, j

    integer :: i1, i2, i3, N, M
    integer :: md
    real(dp) :: deltat


!   check coefficients

    call check ( coefficients, 'implicit_stress_open_boundary', &
      indexarray=[36,48], minimum=[-1,0], maximum=[1,1] )


!   set globals fluid element

    elgrp  = eleminfo%elgrp

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

    call set_viscoelastic_model ( coefficients )

    call set_globals_viscoelastic ( coefficients )

!   coefficients

    cstorage = get_coefficient ( coefficients, index=84, default=1 )

!   set globals object

    object = eleminfo%object

    call set_globals_stokes_object ( mesh, coefficients, object )

    elem   = eleminfo%elem
    elemo  = eleminfo%elemo
    ninti  = eleminfo%nintps  ! reset ninti to the number in nintps


!   allocate arrays

    allocate ( xig(ninti,ndim), x(nodalpo,ndim), xel(ndf,ndim) )
    allocate ( wg(ninti), curvel(ninti), normal(ninti,ndim) )
    allocate ( phi(ninti,ndf), dxdxi(ninti,ndim), theta(ninti,ndfc) )
    allocate ( detF(ninti), dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
    allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
    allocate ( theta1(ninti,nodalpo), dtheta1(ninti,nodalpo,ndim-1) )
    allocate ( xigo(ninti,ndim-1), dxdxis(ninti,ndim,2) )

    allocate ( dtheta(ninti,ndfc,ndim), dthetadx(ninti,ndfc,ndim) )
    allocate ( c(ndfc,ncomp), cn(ndfc,ncompc), cng(ninti,ncompc,nmodes) )
    allocate ( sg(ninti,ncompc) )
    allocate ( gradcn(ninti,ncompc,ndim,nmodes) )
    allocate ( cten(ninti,ndim,ndim), gradcten(ninti,ndim,ndim,ndim) )
    allocate ( hten(ninti,ndim,ndim), rhsd(ninti,ncompc,nmodes) )
    allocate ( Gmod(mode2-mode1+1), Gmodl(mode2-mode1+1) )
    allocate ( work2(ndf,ndim), work5(ndf,ndf,ndim,ndim) )
    allocate ( work10(ndf,ndf), work11(ninti,ndim) )
    allocate ( work(ninti), work1(ninti), work8(ninti) )
    allocate ( work4(ninti,ndim), tmp(ndf,ndim) )
    allocate ( xg(ninti,ndim) )

    allocate ( u(ndf*ndim), uvecmeshnp1(ninti,ndim) )

    if ( coverl ) then
      allocate ( lambda_small(mode2-mode1+1) )
    end if

!   integration points and weighting factors

    xig  = mesh%objects(object)%refcoor_int(eleminfo%intps,:,elemo)
    xigo = mesh%objects(object)%xig(eleminfo%intps,:)
    wg   = mesh%objects(object)%wg(eleminfo%intps)


!   set shape functions

    call set_shape_function ( shapefunc, xig, phi, dphi ) ! velocity
    call set_shape_function ( shapefuncc, xig, theta, dtheta ) ! c tensor
    call set_shape_function ( shapefunco, xigo, theta1, dtheta1 ) ! shape object

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, xel )

    call isoparametric_deformation ( xel, dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )
    call shape_derivative ( dtheta, Finv, dthetadx )


!   compute deformed object element

    call get_coordinates_object ( mesh, elemo, x, object )

    xg =  matmul ( theta1, x )

    if ( ndim == 2 ) then
      call isoparametric_deformation_curve ( x, dtheta1(:,:,1), dxdxi, curvel, &
        normal )
      if ( coorsys == 1 ) then
        curvel = 2 * pi * xg(:,2) * curvel
      end if
    else if ( ndim == 3 ) then
      call isoparametric_deformation_surface ( x, dtheta1, dxdxis, curvel, &
        normal )
    end if

!   set normal

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    normal = real ( orientation, kind=dp ) * normal


!   get conformation tensor

    if ( c_direct ) then
      call get_cn_standard2
    else
      call get_cn_standard1
    end if


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

!   lambda modified G modulus

    if ( coverl ) then
!     c/lambda added to both sides of the constitutive equation
      if ( coefficients%i(58) > 0 ) then
        lambda_small = vemodel%nonlin(coefficients%i(58),mode1:mode2)
      else
        lambda_small = vemodel%lambda(mode1:mode2)
      end if
      Gmodl = Gmod / ( 1 + deltat / lambda_small )
    else
      Gmodl = Gmod
    end if

    if ( matrix ) then

!     set c in tensor format

      if ( coorsys <= 1 ) then

        cten(:,1,1) = matmul ( cng(:,1,mode1:mode2), Gmodl )
        cten(:,1,2) = matmul ( cng(:,2,mode1:mode2), Gmodl )
        cten(:,2,1) = cten(:,1,2)
        cten(:,2,2) = matmul ( cng(:,3,mode1:mode2), Gmodl )

      else if ( coorsys == 2 ) then

        cten(:,1,1) = matmul ( cng(:,1,mode1:mode2), Gmodl )
        cten(:,1,2) = matmul ( cng(:,2,mode1:mode2), Gmodl )
        cten(:,1,3) = matmul ( cng(:,3,mode1:mode2), Gmodl )
        cten(:,2,1) = cten(:,1,2)
        cten(:,2,2) = matmul ( cng(:,4,mode1:mode2), Gmodl )
        cten(:,2,3) = matmul ( cng(:,5,mode1:mode2), Gmodl )
        cten(:,3,1) = cten(:,1,3)
        cten(:,3,2) = cten(:,2,3)
        cten(:,3,3) = matmul ( cng(:,6,mode1:mode2), Gmodl )

      end if

!     set gradc in tensor format

      if ( coorsys <= 1 ) then

        do i = 1, ndim
          gradcten(:,1,1,i) = matmul ( gradcn(:,1,i,mode1:mode2), Gmodl )
          gradcten(:,1,2,i) = matmul ( gradcn(:,2,i,mode1:mode2), Gmodl )
          gradcten(:,2,1,i) = gradcten(:,1,2,i)
          gradcten(:,2,2,i) = matmul ( gradcn(:,3,i,mode1:mode2), Gmodl )
        end do

      else if ( coorsys == 2 ) then

        do i = 1, ndim
          gradcten(:,1,1,i) = matmul ( gradcn(:,1,i,mode1:mode2), Gmodl )
          gradcten(:,1,2,i) = matmul ( gradcn(:,2,i,mode1:mode2), Gmodl )
          gradcten(:,1,3,i) = matmul ( gradcn(:,3,i,mode1:mode2), Gmodl )
          gradcten(:,2,1,i) = gradcten(:,1,2,i)
          gradcten(:,2,2,i) = matmul ( gradcn(:,4,i,mode1:mode2), Gmodl )
          gradcten(:,2,3,i) = matmul ( gradcn(:,5,i,mode1:mode2), Gmodl )
          gradcten(:,3,1,i) = gradcten(:,1,3,i)
          gradcten(:,3,2,i) = gradcten(:,2,3,i)
          gradcten(:,3,3,i) = matmul ( gradcn(:,6,i,mode1:mode2), Gmodl )
        end do

      end if


      do N = 1, ndf
        do M = 1, ndf

          do i = 1, ndim
            do j = 1, ndim

              do ip = 1, ninti
                ! -( v, (u dot gradc) dot n )
                work(ip) = - dot_product ( gradcten(ip,i,:,j), normal(ip,:) )

                ! ( v, (c dot gradu) dot n )
                work1(ip) = dot_product ( cten(ip,i,:), dphidx(ip,M,:) )

              end do

              work5(N,M,i,j) = sum ( phi(:,N) * phi(:,M) * work(:) &
                * curvel * wg ) &
                + sum ( phi(:,N) * work1(:) * normal(:,j) * curvel * wg )

            end do
          end do


          ! ( v, (gradu^T dot c) dot n )
          do ip = 1, ninti
            work11(ip,:) = matmul ( cten(ip,:,:), normal(ip,:) )
            work8(ip) = dot_product ( work11(ip,:), dphidx(ip,M,:) )
          end do

          work10(N,M) = sum ( phi(:,N) * work8 * curvel * wg )

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
          mode2=mode2 )
      else if ( coorsys == 2 ) then
        call rhs_viscoelastic_relax_3D ( vemodel, cng, rhsd, mode1=mode1, &
          mode2=mode2 )
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

        hten(:,1,1) = matmul ( cng(:,1,mode1:mode2) - 1, Gmod ) &
                         + deltat * matmul ( rhsd(:,1,mode1:mode2), Gmodl )
        hten(:,1,2) = matmul ( cng(:,2,mode1:mode2), Gmod ) &
                         + deltat * matmul ( rhsd(:,2,mode1:mode2), Gmodl )
        hten(:,2,1) = hten(:,1,2)
        hten(:,2,2) = matmul ( cng(:,3,mode1:mode2) - 1, Gmod ) &
                         + deltat * matmul ( rhsd(:,3,mode1:mode2), Gmodl )

      else if ( coorsys == 2 ) then

        hten(:,1,1) = matmul ( cng(:,1,mode1:mode2) - 1, Gmod ) &
                         + deltat * matmul ( rhsd(:,1,mode1:mode2), Gmodl )
        hten(:,1,2) = matmul ( cng(:,2,mode1:mode2), Gmod ) &
                         + deltat * matmul ( rhsd(:,2,mode1:mode2), Gmodl )
        hten(:,1,3) = matmul ( cng(:,3,mode1:mode2), Gmod ) &
                         + deltat * matmul ( rhsd(:,3,mode1:mode2), Gmodl )
        hten(:,2,1) = hten(:,1,2)
        hten(:,2,2) = matmul ( cng(:,4,mode1:mode2) - 1, Gmod ) &
                         + deltat * matmul ( rhsd(:,4,mode1:mode2), Gmodl )
        hten(:,2,3) = matmul ( cng(:,5,mode1:mode2), Gmod ) &
                         + deltat * matmul ( rhsd(:,5,mode1:mode2), Gmodl )
        hten(:,3,1) = hten(:,1,3)
        hten(:,3,2) = hten(:,2,3)
        hten(:,3,3) = matmul ( cng(:,6,mode1:mode2) - 1, Gmod ) &
                         + deltat * matmul ( rhsd(:,6,mode1:mode2), Gmodl )

      end if

!     -(v, hten.n)

      do ip = 1, ninti
        work4(ip,:) = matmul( hten(ip,:,:), normal(ip,:) )
      end do

      do i = 1, ndf
        do j = 1, ndim
          work2(i,j) = sum ( phi(:,i) * work4(:,j) * curvel * wg )
        end do
      end do
      elemvec = reshape ( work2, [ ndf*ndim ] )

    end if


!   deallocate arrays

    call delete ( vemodel )

    deallocate ( xig, x, xel )
    deallocate ( wg, curvel, normal )
    deallocate ( phi, dxdxi, theta )
    deallocate ( detF, dphi, F )
    deallocate ( Finv, dphidx )
    deallocate ( theta1, dtheta1 )
    deallocate ( xigo, dxdxis )
    deallocate ( xg )

    deallocate ( dtheta, dthetadx )
    deallocate ( c, cn, cng )
    deallocate ( sg )
    deallocate ( gradcn )
    deallocate ( cten, gradcten )
    deallocate ( hten, rhsd )
    deallocate ( Gmod, Gmodl )
    deallocate ( work2, work5 )
    deallocate ( work10, work11 )
    deallocate ( work, work1, work8 )
    deallocate ( work4, tmp )

    deallocate ( u, uvecmeshnp1 )

    if ( coverl ) then
      deallocate ( lambda_small )
    end if

  contains

!   Get cn for the standard case (ALE or Euler)
!   This routine computes both cn and gradcn in the integration points
!   from nodal values.

    subroutine get_cn_standard1

      integer :: md, i

      do md = mode1, mode2

        if ( exps .or. cproj ) then

!         projection of c=exp(s) or c=b.b^T has been performed in
!         a separate problem

!         single mode conformation
          call get_conformation ( mesh, oldvectors%p(3)%p, elgrp, elem, &
            oldvectors, cmode=cn, mode=md, isv=3, cst=1 )

        else

!         standard case

!         single mode conformation
          call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
            oldvectors, cmode=c, mode=md )

          if ( logc ) then
!           transform sn to cn
!           NOTE: we are transforming the nodal values sn -> cn directly and
!           define cn in the element with the same interpolation as sn. Maybe it
!           is better to do a global least-square projection of cn = exp(sn) to
!           define the nodal point values of cn. See exps=.true.
            if ( coorsys <= 1 .and. vel3D == 0 ) then
              call conformation_2D_log ( vemodel, c, cn )
            else if ( coorsys == 2 .or. vel3D == 1 ) then
              call conformation_3D_log ( vemodel, c, cn )
            end if
          else if ( bten ) then
!           transform bn to cn
!           NOTE: we are transforming the nodal values bn -> cn directly and
!           define cn in the element with the same interpolation as bn. Maybe it
!           is better to do a global least-square projection of cn = bn.bn^T to
!           define the nodal point values of cn. See cproj=.true.
            if ( coorsys <= 1 .and. vel3D == 0 ) then
              call conformation_2D_b ( c, cn, vemodel%bvariant )
            else if ( coorsys == 2 .or. vel3D == 1 ) then
              call conformation_3D_b ( c, cn, vemodel%bvariant )
            end if
          else
            cn = c
          end if

        end if

!       c in Gauss points

        cng(:,:,md) = matmul ( theta, cn )

!       grad cn term

        do i = 1, ndim
          gradcn(:,:,i,md) = matmul ( dthetadx(:,:,i), cn )
        end do

      end do

    end subroutine get_cn_standard1


!   get cn for the standard case (ALE or Euler)
!   This routine computes cn directly in the Gauss integration points
!   and gradcn from nodal values.

    subroutine get_cn_standard2

      integer :: md, i

      do md = mode1, mode2

!       get s, b or c in the nodes

!       single mode conformation
        call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
          oldvectors, cmode=c, mode=md )

!       s, b or c in Gauss points

        sg = matmul ( theta, c )

!       c in Gauss points directly from exp(s), b.b^T or copy

        if ( logc ) then
!         transform s to c
          if ( coorsys <= 1 .and. vel3D == 0 ) then
            call conformation_2D_log ( vemodel, sg, cng(:,:,md) )
          else if ( coorsys == 2 .or. vel3D == 1 ) then
            call conformation_3D_log ( vemodel, sg, cng(:,:,md) )
          end if
        else if ( bten ) then
!         transform b to c
          if ( coorsys <= 1 .and. vel3D == 0 ) then
            call conformation_2D_b ( sg, cng(:,:,md), vemodel%bvariant )
          else if ( coorsys == 2 .or. vel3D == 1 ) then
            call conformation_3D_b ( sg, cng(:,:,md), vemodel%bvariant )
          end if
        else
!         copy
          cng(:,:,md) = sg
        end if

!       gradc

!       first compute c in the nodes

        if ( exps .or. cproj ) then

!         projection of c=exp(s) or c=b.b^T has been performed in
!         a separate problem

!         single mode conformation
          call get_conformation ( mesh, oldvectors%p(3)%p, elgrp, elem, &
            oldvectors, cmode=cn, mode=md, isv=3, cst=1 )

        else

!         direct nodal method (no projection)

          if ( logc ) then
!           transform sn to cn
!           NOTE: we are transforming the nodal values sn -> cn directly and
!           define cn in the element with the same interpolation as sn. Maybe it
!           is better to do a global least-square projection of cn = exp(sn) to
!           define the nodal point values of cn. See exps=.true.
            if ( coorsys <= 1 .and. vel3D == 0 ) then
              call conformation_2D_log ( vemodel, c, cn )
            else if ( coorsys == 2 .or. vel3D == 1 ) then
              call conformation_3D_log ( vemodel, c, cn )
            end if
          else if ( bten ) then
!           transform bn to cn
!           NOTE: we are transforming the nodal values bn -> cn directly and
!           define cn in the element with the same interpolation as bn. Maybe it
!           is better to do a global least-square projection of cn = bn.bn^T to
!           define the nodal point values of cn. See cproj=.true.
            if ( coorsys <= 1 .and. vel3D == 0 ) then
              call conformation_2D_b ( c, cn, vemodel%bvariant )
            else if ( coorsys == 2 .or. vel3D == 1 ) then
              call conformation_3D_b ( c, cn, vemodel%bvariant )
            end if
          else
            cn = c
          end if

        end if

!       grad cn term

        do i = 1, ndim
          gradcn(:,:,i,md) = matmul ( dthetadx(:,:,i), cn )
        end do

      end do

    end subroutine get_cn_standard2

  end subroutine implicit_stress_open_boundary


! L2 projection of exp(s) on the standard discretized space spanned by the
! shape functions

  subroutine exps_projection_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    implicit none

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, m


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1  )

!   coefficients

    if ( first ) then
      cstorage = get_coefficient ( coefficients, index=84, default=1 )
    end if

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( vector ) then

!     get conformation tensor

      do m = mode1, mode2

!       single mode conformation
        call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
          oldvectors, cmode=snod, mode=m )

        ! s in Gauss points
        s = matmul ( theta, snod )

        ! c = exp(s)
        if ( coorsys <= 1 .and. vel3D == 0 ) then
          call conformation_2D_log ( vemodel, s, c )
        else if ( coorsys == 2 .or. vel3D == 1 ) then
          call conformation_3D_log ( vemodel, s, c )
        end if

        ! c in Gauss points
        cng(:,:,m) = c

      end do

    end if


    if ( matrix ) then

      do i = 1, ndfc
        do j = 1, ndfc
          elemmat(i,j) = sum ( theta(:,i) * theta(:,j) * detF * wg )
        end do
      end do

    end if


    if ( vector ) then

      do m = mode1, mode2
        do j = 1, ncompc
          do i = 1, ndfc
            work6(i,j,m) = sum ( theta(:,i) * cng(:,j,m) * detF * wg )
          end do
        end do
      end do

      elemvec = reshape ( work6, [ ndfc * ncompc * ( mode2 - mode1 + 1 ) ] )

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )

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

      allocate ( snod(ndfc,ncomp), s(ninti,ncomp) )
      allocate ( c(ninti,ncompc), cng(ninti,ncompc,nmodes) )
      allocate ( work6(ndfc,ncompc,mode2-mode1+1) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( snod, s )
      deallocate ( c, cng )
      deallocate ( work6 )

    end subroutine deallocate_arrays

  end subroutine exps_projection_elem



! L2 projection of c=b*b^T on the standard discretized space spanned by
! the shape functions

  subroutine c_projection_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    implicit none

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, m


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1  )

!   coefficients

    if ( first ) then
      cstorage = get_coefficient ( coefficients, index=84, default=1 )
    end if

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( vector ) then

!     get conformation tensor

      do m = mode1, mode2

!       single mode conformation
        call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
          oldvectors, cmode=bnod, mode=m )

        ! b in Gauss points
        b = matmul ( theta, bnod )

        ! c = b.b^T
        if ( coorsys <= 1 .and. vel3D == 0 ) then
          call conformation_2D_b ( b, c, vemodel%bvariant )
        else if ( coorsys == 2 .or. vel3D == 1 ) then
          call conformation_3D_b ( b, c, vemodel%bvariant )
        end if

        ! c in Gauss points
        cng(:,:,m) = c

      end do

    end if


    if ( matrix ) then

      do i = 1, ndfc
        do j = 1, ndfc
          elemmat(i,j) = sum ( theta(:,i) * theta(:,j) * detF * wg )
        end do
      end do

    end if


    if ( vector ) then

      do m = mode1, mode2
        do j = 1, ncompc
          do i = 1, ndfc
            work6(i,j,m) = sum ( theta(:,i) * cng(:,j,m) * detF * wg )
          end do
        end do
      end do

      elemvec = reshape ( work6, [ ndfc * ncompc * ( mode2 - mode1 + 1 ) ] )

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )

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

      allocate ( bnod(ndfc,ncomp), b(ninti,ncomp) )
      allocate ( c(ninti,ncompc), cng(ninti,ncompc,nmodes) )
      allocate ( work6(ndfc,ncompc,mode2-mode1+1) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( bnod, b )
      deallocate ( c, cng )
      deallocate ( work6 )

    end subroutine deallocate_arrays

  end subroutine c_projection_elem


! L2 projection of b'=sqrt(b*b^T) on the standard discretized space
! spanned by the shape functions

  subroutine b_projection_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    implicit none

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, m


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1  )

!   coefficients

    if ( first ) then
      cstorage = get_coefficient ( coefficients, index=84, default=1 )
    end if

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( vector ) then

!     get contravariant deformation tensor

      do m = mode1, mode2

!       single mode conformation
        call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
          oldvectors, cmode=bnod, mode=m )

        ! b in Gauss points
        b = matmul ( theta, bnod )

        ! b' = sqrt(b*b^T)
        if ( coorsys <= 1 .and. vel3D == 0 ) then
          call sqrtc_2D_b ( b, c, bvariant=vemodel%bvariant )
        else if ( coorsys == 2 .or. vel3D == 1 ) then
          call sqrtc_3D_b ( b, c, bvariant=vemodel%bvariant )
        end if

        ! b' in Gauss points
        cng(:,:,m) = c

      end do

    end if


    if ( matrix ) then

      do i = 1, ndfc
        do j = 1, ndfc
          elemmat(i,j) = sum ( theta(:,i) * theta(:,j) * detF * wg )
        end do
      end do

    end if


    if ( vector ) then

      do m = mode1, mode2
        do j = 1, ncompc
          do i = 1, ndfc
            work6(i,j,m) = sum ( theta(:,i) * cng(:,j,m) * detF * wg )
          end do
        end do
      end do

      elemvec = reshape ( work6, [ ndfc * ncompc * ( mode2 - mode1 + 1 ) ] )

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )

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

      allocate ( bnod(ndfc,ncomp), b(ninti,ncomp) )
      allocate ( c(ninti,ncompc), cng(ninti,ncompc,nmodes) )
      allocate ( work6(ndfc,ncompc,mode2-mode1+1) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( bnod, b )
      deallocate ( c, cng )
      deallocate ( work6 )

    end subroutine deallocate_arrays

  end subroutine b_projection_elem


! compute conformation tensor component in all nodes

  subroutine deriv_conformation ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

!   choose transformation scheme

    if ( coefficients%i(21) == 0 .and. coefficients%i(71) == 0  ) then
      ! standard scheme
      call deriv_conformation_std ( mesh, problem, elgrp, elem, first, last, &
        coefficients, oldvectors, elemvec, elemwts )
    else
!     transformation (log or b-tensor)
      call deriv_conformation_trn ( mesh, problem, elgrp, elem, first, last, &
        coefficients, oldvectors, elemvec, elemwts )
    end if

  end subroutine deriv_conformation


! compute conformation tensor component in all nodes (standard)

  subroutine deriv_conformation_std ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: comp, mode


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic ( coefficients, maxvar=1 )

!     check coefficients

      call check ( coefficients, 'deriv_conformation_std', &
        indexarray=[13,28], minimum=[1,0], maximum=[ncomp,nmodes] )

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

      allocate ( theta(nodalp,ndfc) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( cnod(ndfc,ncomp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function ( shapefuncc, xrnod, theta )

    end if

    mode = max ( coefficients%i(28), 1 )

!   single mode conformation
    call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
      cmode=cnod, mode=mode )

    comp = coefficients%i(13)

    elemvec = matmul ( theta, cnod(:,comp) )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( vemodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( cnod )

    end if

  end subroutine deriv_conformation_std


! compute conformation tensor in all nodes (transform version)

  subroutine deriv_conformation_trn ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    integer :: comp, mode


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic ( coefficients, maxvar=1 )

      if ( .not. ( logc .or. bten ) ) then
        write(*,'(/a/a/a/)') 'Error in deriv_conformation_trn:', &
        ' This routine can only be used for the log or b-tensor.', &
        ' Use the routine deriv_conformation_std instead.'
        stop
      end if

!     check coefficients

      call check ( coefficients, 'deriv_conformation_trn', &
        indexarray=[13,28], minimum=[1,0], maximum=[ncompc,nmodes] )

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

      allocate ( theta(nodalp,ndfc) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( s(ndfc,ncomp), snod(nodalp,ncomp), cnod(nodalp,ncompc) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function ( shapefuncc, xrnod, theta )

    end if

!   get log conformation tensor

    mode = max ( coefficients%i(28), 1 )

!   single mode conformation
    call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
      cmode=s, mode=mode )

    snod = matmul ( theta, s )

    if ( logc ) then
      ! c = exp(s)
      if ( coorsys <= 1 .and. vel3D == 0 ) then
        call conformation_2D_log ( vemodel, snod, cnod )
      else if ( coorsys == 2 .or. vel3D == 1 ) then
        call conformation_3D_log ( vemodel, snod, cnod )
      end if
    else if ( bten ) then
      ! c = b.b^T
      if ( coorsys <= 1 .and. vel3D == 0 ) then
        call conformation_2D_b ( snod, cnod, vemodel%bvariant )
      else if ( coorsys == 2 .or. vel3D == 1 ) then
        call conformation_3D_b ( snod, cnod, vemodel%bvariant )
      end if
    end if

    comp = coefficients%i(13)

    elemvec = cnod(:,comp)

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( vemodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( s, snod, cnod )

    end if

  end subroutine deriv_conformation_trn


! compute conformation tensor in all nodes

  subroutine deriv_conformation_tensor ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

!   choose transformation scheme

    if ( coefficients%i(21) == 0 .and. coefficients%i(71) == 0  ) then
      ! standard scheme
      call deriv_conformation_tensor_std ( mesh, problem, elgrp, elem, &
        first, last, coefficients, oldvectors, elemvec, elemwts )
    else
!     transformation (log or b-tensor)
      call deriv_conformation_tensor_trn ( mesh, problem, elgrp, elem, &
        first, last, coefficients, oldvectors, elemvec, elemwts )
    end if

  end subroutine deriv_conformation_tensor


! compute conformation tensor in all nodes (standard)

  subroutine deriv_conformation_tensor_std ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: mode


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic ( coefficients, maxvar=1, gammap=.true. )

!     check coefficients

      call check ( coefficients, 'deriv_conformation_tensor_std', &
        indexarray=[28], minimum=[0], maximum=[nmodes] )

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

      allocate ( theta(nodalp,ndfc) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( c(ndfc,ncomp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function ( shapefuncc, xrnod, theta )

    end if

    mode = max ( coefficients%i(28), 1 )

!   single mode conformation
    call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
      cmode=c, mode=mode )

    elemvec = reshape ( matmul ( theta, c ), [nodalp*ncomp] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( vemodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( c )

    end if

  end subroutine deriv_conformation_tensor_std


! compute conformation tensor in all nodes (tranform version)

  subroutine deriv_conformation_tensor_trn ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    integer :: mode


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic ( coefficients, maxvar=1, gammap=.true. )

      if ( .not. ( logc .or. bten ) ) then
        write(*,'(/a/a/a/)') 'Error in deriv_conformation_tensor_trn:', &
        ' This routine can only be used for the log or b-tensor.', &
        ' Use the routine deriv_conformation_std instead.'
        stop
      end if

!     check coefficients

      call check ( coefficients, 'deriv_conformation_tensor_trn', &
        indexarray=[28], minimum=[0], maximum=[nmodes] )

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

      allocate ( theta(nodalp,ndfc) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( s(ndfc,ncomp), snod(nodalp,ncomp), cnod(nodalp,ncompc) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function ( shapefuncc, xrnod, theta )

    end if

!   get log conformation tensor

    mode = max ( coefficients%i(28), 1 )

!   single mode conformation
    call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
      cmode=s, mode=mode )

    snod = matmul ( theta, s )

    if ( logc ) then
      ! c = exp(s)
      if ( coorsys <= 1 .and. vel3D == 0 ) then
        call conformation_2D_log ( vemodel, snod, cnod )
      else if ( coorsys == 2 .or. vel3D == 1 ) then
        call conformation_3D_log ( vemodel, snod, cnod )
      end if
    else if ( bten ) then
      ! c = b.b^T
      if ( coorsys <= 1 .and. vel3D == 0 ) then
        call conformation_2D_b ( snod, cnod, vemodel%bvariant )
      else if ( coorsys == 2 .or. vel3D == 1 ) then
        call conformation_3D_b ( snod, cnod, vemodel%bvariant )
      end if
    end if

    elemvec = reshape ( cnod, [nodalp*ncompc] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( vemodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( s, snod, cnod )

    end if

  end subroutine deriv_conformation_tensor_trn


! sample value of conformation tensor in one node of the object

  subroutine sample_conformation_tensor ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, cval )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: cval

!   choose transformation scheme

    if ( coefficients%i(21) == 0 .and. coefficients%i(71) == 0  ) then
!     standard scheme
      call sample_conformation_tensor_std ( mesh, problem, object, nodeobj, &
        coefficients, oldvectors, cval )
    else
!     transformation (log or b-tensor)
      call sample_conformation_tensor_trn ( mesh, problem, object, nodeobj, &
        coefficients, oldvectors, cval )
    end if

  end subroutine sample_conformation_tensor


! sample value of conformation in one node of the object (standard)

  subroutine sample_conformation_tensor_std ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, cval )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: cval


    integer :: elgrp, elem, mode
    real(dp) :: xr(1,mesh%ndim)


    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

    call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )
    call set_viscoelastic_model ( coefficients )
    call set_globals_viscoelastic ( coefficients, maxvar=1, gammap=.true. )

!   check coefficients

    call check ( coefficients, 'sample_conformation_tensor_std', &
      indexarray=[28], minimum=[0], maximum=[nmodes] )

    cstorage = get_coefficient ( coefficients, index=84, default=1 )

    allocate ( theta(1,ndfc), c(ndfc,ncomp) )

!   set shape function in the point

    call set_shape_function ( shapefuncc, xr, theta )

    mode = max ( coefficients%i(28), 1 )

!   single mode conformation
    call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
      cmode=c, mode=mode )

    cval = matmul( theta(1,:), c )

!   delete data

    call delete ( vemodel )

    deallocate ( theta, c )

  end subroutine sample_conformation_tensor_std


! sample value of conformation in one node of the object (standard)

  subroutine sample_conformation_tensor_trn ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, cval )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: cval


    integer :: elgrp, elem, mode
    real(dp) :: xr(1,mesh%ndim)


    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

    call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )
    call set_viscoelastic_model ( coefficients )
    call set_globals_viscoelastic ( coefficients, maxvar=1, gammap=.true. )

!   check coefficients

    call check ( coefficients, 'sample_conformation_tensor_trn', &
      indexarray=[28], minimum=[0], maximum=[nmodes] )

    cstorage = get_coefficient ( coefficients, index=84, default=1 )

    allocate ( theta(1,ndfc), s(ndfc,ncomp) )
    allocate ( snod(1,ncomp), cnod(1,ncompc) )

!   set shape function in the point

    call set_shape_function ( shapefuncc, xr, theta )

    mode = max ( coefficients%i(28), 1 )

!   single mode conformation
    call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
      cmode=s, mode=mode )

    snod = matmul( theta, s )

    if ( logc ) then
      ! c = exp(s)
      if ( coorsys <= 1 .and. vel3D == 0 ) then
        call conformation_2D_log ( vemodel, snod, cnod )
      else if ( coorsys == 2 .or. vel3D == 1 ) then
        call conformation_3D_log ( vemodel, snod, cnod )
      end if
    else if ( bten ) then
      ! c = b.b^T
      if ( coorsys <= 1 .and. vel3D == 0 ) then
        call conformation_2D_b ( snod, cnod, vemodel%bvariant )
      else if ( coorsys == 2 .or. vel3D == 1 ) then
        call conformation_3D_b ( snod, cnod, vemodel%bvariant )
      end if
    end if

    cval = cnod(1,:)

!   delete data

    call delete ( vemodel )

    deallocate ( theta, s )
    deallocate ( snod, cnod )

  end subroutine sample_conformation_tensor_trn


! viscoelastic stress (scalar components) in all nodes

  subroutine deriv_viscoelastic_stress_scalar ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: m, ip, comp
    real(dp) :: chpar, tau_e, mu, tau_y


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic ( coefficients, maxvar=1, gammap=.true. )

!     check coefficients

      call check ( coefficients, 'deriv_viscoelastic_stress_scalar', &
        indexarray=[13,28,95], minimum=[1,0,0], maximum=[3,nmodes,nmodes] )

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

      allocate ( theta(nodalp,ndfc) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( c(ndfc,ncomp) )
      allocate ( cg(nodalp,ncomp,nmodes) )
      allocate ( tauvec(nodalp,ncompt), tauten(nodalp,ncompu,ncompu) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function ( shapefuncc, xrnod, theta )

      if ( varpar ) then
        allocate ( phi(nodalp,ndf) )
        call set_shape_function ( shapefunc, xrnod, phi )
        allocate ( mvemodel(nodalp) )
        mvemodel = vemodel
      end if

      if ( vemcompressible ) then
        allocate ( psi(nodalp,ndfp) )
        call set_shape_function ( shapefuncp, xrnod, psi )
        allocate ( pr(ndfp), press(nodalp) )
        call create_vemopt ( vemopt, dep_J=.true., np=nodalp )
      else
        call create_vemopt ( vemopt )
      end if

    end if

!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure ( mesh, oldvectors%p(1)%p, elgrp, elem, &
        coefficients, oldvectors%s(1)%p, vemopt )

    end if

!   get conformation tensors

    do m = mode1, mode2

!     single mode conformation
      call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
        cmode=c, mode=m )

      cg(:,:,m) = matmul ( theta, c )

    end do

    if ( varpar ) then

!     determine the variable material parameters in mvemodel

      call evaluate_mve ( mesh, problem, elem, elgrp, coefficients, &
        oldvectors, modulus=.true., nonlin=.true. )

    end if

    comp = coefficients%i(13)

!   stress tensor: use of global variables; result in tauvec

    if ( coefficients%i(28) == 0 ) then
!     total stress
      call stress_tensor_viscoelastic ( cg, vemopt=vemopt )
      if ( comp == 3 ) m = coefficients%i(95) ! mode number for SRM DP model
    else
!     stress of one mode
      call stress_tensor_viscoelastic ( cg, mode=coefficients%i(28), &
        vemopt=vemopt )
      if ( comp == 3 ) m = coefficients%i(28) ! mode number for SRM DP model
    end if

    if ( coorsys <= 1. .and. vel3D == 0 ) then

!     2D

      select case ( comp )
      case(1) ! von mises equivalent shear stress
        do ip = 1, nodalp
          elemvec(ip) = vonmises_2D ( tauvec(ip,:) )
        end do
      case(2) ! trace of stress tensor
        do ip = 1, nodalp
          elemvec(ip) = trace ( tauvec(ip,:) )
        end do
      case(3) ! regime in Saramito DP model
        do ip = 1, nodalp
          if ( varpar ) then
            tau_y = mvemodel(ip)%nonlin(1,m)
            mu = mvemodel(ip)%nonlin(2,m)
          else
            tau_y = vemodel%nonlin(1,m)
            mu = vemodel%nonlin(2,m)
          end if
          tau_e = vonmises_2D ( tauvec(ip,:) )
          chpar = mu * trace ( tauvec(ip,:) ) - 3 * tau_y
          if ( chpar <= - 3 * tau_e ) then
!           regime I: sticking
            elemvec(ip) = 1
          else if ( chpar >= 2 * mu**2 * tau_e ) then
!           regime III: loosing contact
            elemvec(ip) = 3
          else
!           regime II: sliding
            elemvec(ip) = 2
          end if
        end do
      case default
        call errormsg_case_default ( 'deriv_viscoelastic_stress_scalar', &
          'comp', int_value=comp )
      end select

    else if ( coorsys == 2 .or. vel3D == 1 ) then

!     3D

      select case ( comp )
      case(1) ! von mises equivalent shear stress
        do ip = 1, nodalp
          elemvec(ip) = vonmises_3D ( tauvec(ip,:) )
        end do
      case(2) ! trace of stress tensor
        do ip = 1, nodalp
          elemvec(ip) = trace ( tauvec(ip,:) )
        end do
      case(3) ! regime in Saramito DP model
        do ip = 1, nodalp
          if ( varpar ) then
            tau_y = mvemodel(ip)%nonlin(1,m)
            mu = mvemodel(ip)%nonlin(2,m)
          else
            tau_y = vemodel%nonlin(1,m)
            mu = vemodel%nonlin(2,m)
          end if
          tau_e = vonmises_3D ( tauvec(ip,:) )
          chpar = mu * trace ( tauvec(ip,:) ) - 3 * tau_y
          if ( chpar <= - 3 * tau_e ) then
!           regime I: sticking
            elemvec(ip) = 1
          else if ( chpar >= 2 * mu**2 * tau_e ) then
!           regime III: loosing contact
            elemvec(ip) = 3
          else
!           regime II: sliding
            elemvec(ip) = 2
          end if
        end do
      case default
        call errormsg_case_default ( 'deriv_viscoelastic_stress_scalar', &
          'comp', int_value=comp )
      end select

    end if

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( vemodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( c )
      deallocate ( cg )
      deallocate ( tauvec, tauten )

      if ( varpar ) then
        deallocate ( phi )
        do ip = 1, nodalp
          call delete ( mvemodel(ip) )
        end do
        deallocate ( mvemodel )
      end if

      if ( vemcompressible ) then
        deallocate ( psi )
        deallocate ( pr, press )
      end if

      call delete ( vemopt )

    end if

  contains

    function trace ( tau )

      real(dp), dimension(:), intent(in) :: tau
      real(dp) :: trace

      if ( ncompt == 3 ) then
        trace = tau(1) + tau(3)
      else if ( ncompt == 4 ) then
        trace = tau(1) + tau(3) + tau(4)
      else if ( ncompt == 6 ) then
        trace = tau(1) + tau(4) + tau(6)
      else
        write(*,'(/2(a/),a,i0,/)') &
          'Error in deriv_viscoelastic_stress_scalar:', &
          '  incorrect value for computing trace of stress tensor', &
          '  ncompt = ', ncompt
        stop
      end if

    end function trace

  end subroutine deriv_viscoelastic_stress_scalar


! viscoelastic stress tensor in all nodes

  subroutine deriv_viscoelastic_stress_tensor ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: m, ip


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic ( coefficients, maxvar=1, gammap=.true. )

!     check coefficients

      call check ( coefficients, 'deriv_viscoelastic_stress_tensor', &
        indexarray=[28], minimum=[0], maximum=[nmodes] )

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

!     test size of elemvec

      if ( size(elemvec) /= nodalp*ncompt ) then
        write(*,'(/2(a/))') 'Error in deriv_viscoelastic_stress_tensor:', &
          ' element vector has incorrect size for a stress tensor'
        stop
      end if

      allocate ( theta(nodalp,ndfc) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( c(ndfc,ncomp) )
      allocate ( cg(nodalp,ncomp,nmodes) )
      allocate ( tauvec(nodalp,ncompt), tauten(nodalp,ncompu,ncompu) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function ( shapefuncc, xrnod, theta )

      if ( varpar ) then
        allocate ( phi(nodalp,ndf) )
        call set_shape_function ( shapefunc, xrnod, phi )
        allocate ( mvemodel(nodalp) )
        mvemodel = vemodel
      end if

      if ( vemcompressible ) then
        allocate ( psi(nodalp,ndfp) )
        call set_shape_function ( shapefuncp, xrnod, psi )
        allocate ( pr(ndfp), press(nodalp) )
        call create_vemopt ( vemopt, dep_J=.true., np=nodalp )
      else
        call create_vemopt ( vemopt )
      end if

    end if

!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure ( mesh, oldvectors%p(1)%p, elgrp, elem, &
        coefficients, oldvectors%s(1)%p, vemopt )

    end if

!   get conformation tensors

    do m = mode1, mode2

!     single mode conformation
      call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
        cmode=c, mode=m )

      cg(:,:,m) = matmul ( theta, c )

    end do

    if ( varpar ) then

!     determine the variable material parameters in mvemodel

      call evaluate_mve ( mesh, problem, elem, elgrp, coefficients, &
        oldvectors, modulus=.true., nonlin=.true. )

    end if

!   stress tensor: use of global variables; result in tauvec

    if ( coefficients%i(28) == 0 ) then
!     total stress
      call stress_tensor_viscoelastic ( cg, vemopt=vemopt )
    else
!     stress of one mode
      call stress_tensor_viscoelastic ( cg, mode=coefficients%i(28), &
        vemopt=vemopt )
    end if

    elemvec = reshape ( tauvec, [nodalp*ncompt] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( vemodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( c )
      deallocate ( cg )
      deallocate ( tauvec, tauten )

      if ( varpar ) then
        deallocate ( phi )
        do ip = 1, nodalp
          call delete ( mvemodel(ip) )
        end do
        deallocate ( mvemodel )
      end if

      if ( vemcompressible ) then
        deallocate ( psi )
        deallocate ( pr, press )
      end if

      call delete ( vemopt )

    end if

  end subroutine deriv_viscoelastic_stress_tensor


! Internal element routine to fill the viscoelastic stress work in the Gauss
! points. The result is stored in a vector defined per element (elvector).
! This element should be used together with the routine loop_over_elements.
! The resulting elvector can be used as input to other elements based on
! the energy equation.

  subroutine fill_viscoelastic_stress_work_gauss ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors

    logical :: standard_c
    integer :: j, ip, m, pos_elvec


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1, gammap=.true. )

!   position of elvector depends on the number of non-linear+alam parameters

    pos_elvec = 3 + max(1,size(vemodel%nonlin,1)+size(vemodel%alam,1))

!   coefficients

    if ( first ) then
      cstorage = get_coefficient ( coefficients, index=84, default=1 )
    end if

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

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    call isoparametric_coordinates ( x, phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   get velocity vector

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    ugvector = matmul ( phi, uvector )

!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure ( mesh,oldvectors%p(1)%p, elgrp, elem, &
        coefficients, oldvectors%s(1)%p, vemopt )

    end if

!   gradu

    do j = 1, ndim
      gradu(:,1:ncompu,j) = matmul ( dphidx(:,:,j), uvector )
    end do

    if ( coorsys == 1 .and. vel3D == 0 ) then

!     axisymmetric

      gradu(:,:,3) = 0
      gradu(:,3,:) = 0

      where ( xg(:,2) < 1e-10_dp)
        gradu(:,3,3) = gradu(:,2,2)              ! du_r / dr (r=0)
      elsewhere
        gradu(:,3,3) = ugvector(:,2) / xg(:,2)   ! u_r / r
      end where

    else if ( vel3D == 1 ) then

!     three velocity components

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

!   rate-of-deformation tensor D

    do ip = 1, ninti
      Dten(ip,:,:) = ( gradu(ip,:,:) + transpose(gradu(ip,:,:)) ) / 2
    end do

!   get conformation tensors

    standard_c = coefficients%i(74) == 1 .and. ( exps .or. cproj )

    do m = mode1, mode2

      if ( standard_c ) then

!       projection of c=exp(s) or c=b.b^T has been performed in
!       a separate problem and used here as the conformation tensor

!       single mode conformation
        call get_conformation ( mesh, oldvectors%p(3)%p, elgrp, elem, &
          oldvectors, cmode=c, mode=m, isv=2, cst=1 )

      else

!       single mode conformation
        call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
          cmode=c, mode=m )

      end if

      cg(:,:,m) = matmul ( theta, c )

    end do

    if ( varpar ) then

!     determine the variable material parameters in mvemodel

      call evaluate_mve ( mesh, problem, elem, elgrp, coefficients, &
        oldvectors, modulus=.true., nonlin=.true. )

    end if

!   stress tensor: use of global variables; result in tauten and tauvec

    call stress_tensor_viscoelastic ( cg, standard_c, vemopt=vemopt )

    if ( coorsys == 1 .and. vel3D == 0 ) then
!     axisymmetric
      tauten(:,:,3) = 0
      tauten(:,3,:) = 0
      tauten(:,3,3) = tauvec(:,4)
    end if

!   stress work tau : D

    do ip = 1, ninti
      work(ip) = sum ( tauten(ip,:,:) * Dten(ip,:,:) )
    end do

    call put_elvector ( mesh, oldvectors%e(pos_elvec)%p, elgrp, elem, r1=work )

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )


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

      allocate ( u(ncompu*ndf), ugvector(ninti,ncompu) )
      allocate ( uvector(ndf,ncompu) )
      allocate ( work(ninti) )
      allocate ( c(ndfc,ncomp) )
      allocate ( cg(ninti,ncomp,nmodes) )
      allocate ( tauvec(ninti,ncompt) )

      if ( coorsys <= 1 .and. vel3D == 0 ) then
        allocate ( gradu(ninti,ndim+coorsys,ndim+coorsys) )
        allocate ( Dten(ninti,ndim+coorsys,ndim+coorsys) )
        allocate ( tauten(ninti,ndim+coorsys,ndim+coorsys) )
      else
        allocate ( gradu(ninti,ncompu,ncompu) )
        allocate ( Dten(ninti,ncompu,ncompu) )
        allocate ( tauten(ninti,ncompu,ncompu) )
      end if

      if ( varpar ) then
        allocate ( mvemodel(ninti) )
        mvemodel = vemodel
      end if

      if ( vemcompressible ) then
        allocate ( pr(ndfp), press(ninti) )
        call create_vemopt ( vemopt, dep_J=.true., np=ninti )
      else
        call create_vemopt ( vemopt )
      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      integer :: ip

      deallocate ( u, ugvector )
      deallocate ( uvector )
      deallocate ( gradu )
      deallocate ( Dten )
      deallocate ( tauten )
      deallocate ( work )
      deallocate ( c )
      deallocate ( cg )
      deallocate ( tauvec )

      if ( varpar ) then
        do ip = 1, ninti
          call delete ( mvemodel(ip) )
        end do
        deallocate ( mvemodel )
      end if

      if ( vemcompressible ) then
        deallocate ( pr, press )
      end if

      call delete ( vemopt )

    end subroutine deallocate_arrays

  end subroutine fill_viscoelastic_stress_work_gauss


! Element for the drag force on a CURVE due to the polymer stress using
! integrate_boundary_elements in the postprocessing_m module.
!
! oldvectors:
!  v1(1)%p(1:nmodes) = conformation tensors derived from subroutine
!                "deriv_conformation_tensor_std" for _both_ the standard
!                and the log conformation. Thus in the log formulation
!                they contain the nodal values of logc.

  subroutine viscoelastic_drag_curve ( mesh, problem, curve, elem, first, &
    last, coefficients, oldvectors, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, ip, m


    if ( first ) then

!     first element on this curve/surface

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, curve=curve, &
        maxvel3D=1 )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic_boun ( coefficients, gammap=.true. )

!     allocate arrays

      allocate ( wg(ninti), curvel(ninti) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,ndim-1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim-1), dxdxi(ninti,ndim) )
      allocate ( work(ndf*ncomp), c(ndf,ncomp), cg(ninti,ncomp,nmodes) )
      allocate ( work4(ninti,ncompu) )
      allocate ( tauten(ninti,ncompu,ncompu), tauvec(ninti,ncompt) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( vemcompressible ) then
        allocate ( psi(ninti,ndfp) )
        call set_shape_function ( shapefuncp, xig, psi )
        allocate ( pr(ndfp), press(ninti) )
        call create_vemopt ( vemopt, dep_J=.true., np=ninti )
      else
        call create_vemopt ( vemopt )
      end if

    end if

!   geometry

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )

    call isoparametric_deformation_curve ( x, dphi(:,:,1), dxdxi, curvel, &
      normal )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x, phi, xg )
      curvel = 2 * pi * xg(:,2) * curvel
    end if

!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure_geometry ( mesh, problem, elem, coefficients, &
        oldvectors%s(1)%p, vemopt, curve=curve )

    end if

!   conformation tensor

    do m = mode1, mode2

      call get_vector_geometry ( mesh, problem, oldvectors%v1(1)%p(m), elem, &
        work, curve=curve, layer=layer )

      c = reshape ( work, [ndf,ncomp] )

      cg(:,:,m) = matmul ( phi, c )

    end do

!   stress tensor: use of global variables; result in tauten

    call stress_tensor_viscoelastic ( cg, vemopt=vemopt )

!   traction force in each integration point: {t}=[tau]{n}

    do ip = 1, ninti
      work4(ip,:) = matmul( tauten(ip,:,1:ndim), normal(ip,:) )
    end do

!   integrate traction to force

    do i = 1, ncompu
      elemvec(i) = sum ( work4(:,i) * curvel * wg )
    end do

    if ( last ) then

!     last element on this curve

      call delete ( vemodel )

      deallocate ( wg, curvel )
      deallocate ( normal )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi, dxdxi )
      deallocate ( work, c, cg )
      deallocate ( work4 )
      deallocate ( tauten, tauvec )

      if ( vemcompressible ) then
        deallocate ( psi )
        deallocate ( pr, press )
      end if

      call delete ( vemopt )

    end if

  end subroutine viscoelastic_drag_curve


! Element for the drag force on a SURFACE due to the polymer stress using
! integrate_boundary_elements in the postprocessing_m module.
!
! oldvectors:
!  v1(1)%p(1:nmodes) = conformation tensors derived from subroutine
!                "deriv_conformation_tensor_std" for _both_ the standard
!                and the log conformation. Thus in the log formulation
!                they contain the nodal values of logc.

  subroutine viscoelastic_drag_surface ( mesh, problem, surface, elem, first, &
    last, coefficients, oldvectors, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: surface, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, ip, m


    if ( first ) then

!     first element on this surface

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, surface=surface )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic_boun ( coefficients, gammap=.true. )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,ndim-1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim-1) )
      allocate ( work(ndf*ncomp), c(ndf,ncomp), cg(ninti,ncomp,nmodes) )
      allocate ( work4(ninti,ndim) )
      allocate ( tauten(ninti,ndim,ndim), tauvec(ninti,ncompt) )
      allocate ( dxdxis(ninti,ndim,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( vemcompressible ) then
        allocate ( psi(ninti,ndfp) )
        call set_shape_function ( shapefuncp, xig, psi )
        allocate ( pr(ndfp), press(ninti) )
        call create_vemopt ( vemopt, dep_J=.true., np=ninti )
      else
        call create_vemopt ( vemopt )
      end if

    end if

!   geometry

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, normal )

!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure_geometry ( mesh, problem, elem, coefficients, &
        oldvectors%s(1)%p, vemopt, surface=surface )

    end if

!   conformation tensor

    do m = mode1, mode2

      call get_vector_geometry ( mesh, problem, oldvectors%v1(1)%p(m), elem, &
        work, surface=surface, layer=layer )

      c = reshape ( work, [ndf,ncomp] )

      cg(:,:,m) = matmul ( phi, c )

    end do

!   stress tensor: use of global variables; result in tauten

    call stress_tensor_viscoelastic ( cg, vemopt=vemopt )

!   traction force in each integration point: {t}=[tau]{n}

    do ip = 1, ninti
      work4(ip,:) = matmul( tauten(ip,:,:), normal(ip,:) )
    end do

!   integrate traction to force

    do i = 1, ndim
      elemvec(i) = sum ( work4(:,i) * surfl * wg )
    end do

    if ( last ) then

!     last element on this surface

      call delete ( vemodel )

      deallocate ( wg, surfl )
      deallocate ( normal )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi )
      deallocate ( work, c, cg )
      deallocate ( work4 )
      deallocate ( tauten, tauvec )
      deallocate ( dxdxis )

      if ( vemcompressible ) then
        deallocate ( psi )
        deallocate ( pr, press )
      end if

      call delete ( vemopt )

    end if

  end subroutine viscoelastic_drag_surface


! stress tensor viscoelastic
! NOTE: this is not a general usable subroutine, but more like an inlined
!       procedure. Most in/out variables are in a global module and need to
!       be defined/allocated.

  subroutine stress_tensor_viscoelastic ( cvec, standard_c, mode, vemopt )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      coorsys, vel3D, varpar, logc, bten, vemodel, mode1, mode2, mvemodel, &
!     output:
      tauvec, tauten

!   (log, b or standard) conformation tensor in vector format
    real(dp), dimension(:,:,:), intent(in) :: cvec

!   If present and .true. cvec contains the standard c tensor (non-transformed)
!   irrespective of the values of logc, bten.
!   default=.false.
    logical, optional, intent(in) :: standard_c

!   If present the stress in computed for the given mode only.
    integer, optional, intent(in) :: mode

!   if present: additional optional parameters (see type description)
    type(vemopt_t), intent(inout), optional :: vemopt


    integer :: m1, m2
    logical :: lstandard_c

    lstandard_c = set_optional ( variable=standard_c, default=.false. )

    if ( present(mode) ) then
      m1 = mode; m2 = mode
    else
      m1 = mode1; m2 = mode2
    end if

    if ( coorsys <= 1 .and. vel3D == 0 ) then

!     2D and axisymmetric

      if ( varpar ) then

!       variable coefficients

        if ( lstandard_c ) then
          call stress_viscoelastic_2D ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        else if ( logc ) then
          call stress_viscoelastic_2D_log ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        else if ( bten ) then
          call stress_viscoelastic_2D_b ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        else
          call stress_viscoelastic_2D ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        end if

      else

!       fixed coefficients

        if ( lstandard_c ) then
          call stress_viscoelastic_2D ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, vemopt=vemopt )
        else if ( logc ) then
          call stress_viscoelastic_2D_log ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, vemopt=vemopt )
        else if ( bten ) then
          call stress_viscoelastic_2D_b ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, vemopt=vemopt )
        else
          call stress_viscoelastic_2D ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, vemopt=vemopt )
        end if

      end if

      tauten(:,1,1) = tauvec(:,1)
      tauten(:,1,2) = tauvec(:,2)
      tauten(:,2,1) = tauvec(:,2)
      tauten(:,2,2) = tauvec(:,3)

    else if ( coorsys == 2 .or. vel3d == 1 ) then

!     3D

      if ( varpar ) then

!       variable coefficients

        if ( lstandard_c ) then
          call stress_viscoelastic_3D ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        else if ( logc ) then
          call stress_viscoelastic_3D_log ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        else if ( bten ) then
          call stress_viscoelastic_3D_b ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        else
          call stress_viscoelastic_3D ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        end if

      else

!       fixed coefficients

        if ( lstandard_c ) then
          call stress_viscoelastic_3D ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, vemopt=vemopt )
        else if ( logc ) then
          call stress_viscoelastic_3D_log ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, vemopt=vemopt )
        else if ( bten ) then
          call stress_viscoelastic_3D_b ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, vemopt=vemopt )
        else
          call stress_viscoelastic_3D ( vemodel, cvec, tauvec, mode1=m1, &
            mode2=m2, vemopt=vemopt )
        end if

      end if

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

  end subroutine stress_tensor_viscoelastic


! Jacobian of the stress tensor viscoelastic (single mode)
! NOTE: this is not general usable subroutine, but more like an inlined
!       procedure. Most in/out variables are in a global module and need to
!       be defined/allocated.

  subroutine dstress_tensor_viscoelastic ( cvec, standard_c, mode, vemopt )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      coorsys, vel3D, varpar, logc, bten, vemodel, mvemodel, &
!     output:
      dtauvec

!   (log, b or standard) conformation tensor in vector format
    real(dp), dimension(:,:,:), intent(in) :: cvec

!   If present and .true. cvec contains the standard c tensor (non-transformed)
!   irrespective of the values of logc, bten.
!   default=.false.
    logical, optional, intent(in) :: standard_c

!   the mode number
    integer, intent(in) :: mode

!   if present: additional optional parameters (see type description)
    type(vemopt_t), intent(inout), optional :: vemopt


    logical :: lstandard_c

    lstandard_c = set_optional ( variable=standard_c, default=.false. )

    if ( logc .and. .not. lstandard_c ) then
      write(*,'(/2(a/))') 'Error in dstress_tensor_viscoelastic:', &
        ' Jacobian for log c formulation not available.'
      stop
    end if

    if ( coorsys <= 1 .and. vel3D == 0 ) then

!     2D and axisymmetric

      if ( varpar ) then

!       variable coefficients

        if ( lstandard_c ) then
          call dstress_viscoelastic_2D ( vemodel, cvec, dtauvec, mode, &
            mvemodel=mvemodel, vemopt=vemopt )
        else if ( bten ) then
          call dstress_viscoelastic_2D_b ( vemodel, cvec, dtauvec, mode, &
            mvemodel=mvemodel, vemopt=vemopt )
        else
          call dstress_viscoelastic_2D ( vemodel, cvec, dtauvec, mode, &
            mvemodel=mvemodel, vemopt=vemopt )
        end if

      else

!       fixed coefficients

        if ( lstandard_c ) then
          call dstress_viscoelastic_2D ( vemodel, cvec, dtauvec, mode, &
            vemopt=vemopt )
        else if ( bten ) then
          call dstress_viscoelastic_2D_b ( vemodel, cvec, dtauvec, mode, &
            vemopt=vemopt )
        else
          call dstress_viscoelastic_2D ( vemodel, cvec, dtauvec, mode, &
            vemopt=vemopt )
        end if

      end if

    else if ( coorsys == 2 .or. vel3d == 1 ) then

!     3D

      if ( varpar ) then

!       variable coefficients

        if ( lstandard_c ) then
          call dstress_viscoelastic_3D ( vemodel, cvec, dtauvec, mode, &
            mvemodel=mvemodel, vemopt=vemopt )
        else if ( bten ) then
          call dstress_viscoelastic_3D_b ( vemodel, cvec, dtauvec, mode, &
            mvemodel=mvemodel, vemopt=vemopt )
        else
          call dstress_viscoelastic_3D ( vemodel, cvec, dtauvec, mode, &
            mvemodel=mvemodel, vemopt=vemopt )
        end if

      else

!       fixed coefficients

        if ( lstandard_c ) then
          call dstress_viscoelastic_3D ( vemodel, cvec, dtauvec, mode, &
            vemopt=vemopt )
        else if ( bten ) then
          call dstress_viscoelastic_3D_b ( vemodel, cvec, dtauvec, mode, &
            vemopt=vemopt )
        else
          call dstress_viscoelastic_3D ( vemodel, cvec, dtauvec, mode, &
            vemopt=vemopt )
        end if

      end if

    end if

  end subroutine dstress_tensor_viscoelastic


! derivative of the stress tensor viscoelastic with respect to J
! NOTE: this is not a general usable subroutine, but more like an inlined
!       procedure. Most in/out variables are in a global module and need to
!       be defined/allocated.

  subroutine dstressdJ_tensor_viscoelastic ( cvec, standard_c, mode, vemopt )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      coorsys, vel3D, varpar, logc, bten, vemodel, mode1, mode2, mvemodel, &
!     output:
      dtauvecdJ, dtautendJ

!   (log, b or standard) conformation tensor in vector format
    real(dp), dimension(:,:,:), intent(in) :: cvec

!   If present and .true. cvec contains the standard c tensor (non-transformed)
!   irrespective of the values of logc, bten.
!   default=.false.
    logical, optional, intent(in) :: standard_c

!   If present the stress in computed for the given mode only.
    integer, optional, intent(in) :: mode

!   If present: additional optional parameters (see type description)
    type(vemopt_t), intent(inout), optional :: vemopt

    integer :: m1, m2
    logical :: lstandard_c

    if ( logc ) then
      write(*,'(/2(a/))') 'Error in dstressdJ_tensor_viscoelastic:', &
        ' derivative with respect to J for log c formulation not available.'
      stop
    end if

    lstandard_c = set_optional ( variable=standard_c, default=.false. )

    if ( present(mode) ) then
      m1 = mode; m2 = mode
    else
      m1 = mode1; m2 = mode2
    end if

    if ( coorsys <= 1 .and. vel3D == 0 ) then

!     2D and axisymmetric

      if ( varpar ) then

!       variable coefficients

        if ( lstandard_c ) then
          call dstressdJ_viscoelastic_2D ( vemodel, cvec, dtauvecdJ, &
            mode1=m1, mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        else if ( bten ) then
          call dstressdJ_viscoelastic_2D_b ( vemodel, cvec, dtauvecdJ, &
            mode1=m1, mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        else
          call dstressdJ_viscoelastic_2D ( vemodel, cvec, dtauvecdJ, &
            mode1=m1, mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        end if

      else

!       fixed coefficients

        if ( lstandard_c ) then
          call dstressdJ_viscoelastic_2D ( vemodel, cvec, dtauvecdJ, &
            mode1=m1, mode2=m2, vemopt=vemopt )
        else if ( bten ) then
          call dstressdJ_viscoelastic_2D_b ( vemodel, cvec, dtauvecdJ, &
            mode1=m1, mode2=m2, vemopt=vemopt )
        else
          call dstressdJ_viscoelastic_2D ( vemodel, cvec, dtauvecdJ, &
            mode1=m1, mode2=m2, vemopt=vemopt )
        end if

      end if

      dtautendJ(:,1,1) = dtauvecdJ(:,1)
      dtautendJ(:,1,2) = dtauvecdJ(:,2)
      dtautendJ(:,2,1) = dtauvecdJ(:,2)
      dtautendJ(:,2,2) = dtauvecdJ(:,3)

    else if ( coorsys == 2 .or. vel3d == 1 ) then

!     3D

      if ( varpar ) then

!       variable coefficients

        if ( lstandard_c ) then
          call dstressdJ_viscoelastic_3D ( vemodel, cvec, dtauvecdJ, &
            mode1=m1, mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        else if ( bten ) then
          call dstressdJ_viscoelastic_3D_b ( vemodel, cvec, dtauvecdJ, &
            mode1=m1, mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        else
          call dstressdJ_viscoelastic_3D ( vemodel, cvec, dtauvecdJ, &
            mode1=m1, mode2=m2, mvemodel=mvemodel, vemopt=vemopt )
        end if

      else

!       fixed coefficients

        if ( lstandard_c ) then
          call dstressdJ_viscoelastic_3D ( vemodel, cvec, dtauvecdJ, &
            mode1=m1, mode2=m2, vemopt=vemopt )
        else if ( bten ) then
          call dstressdJ_viscoelastic_3D_b ( vemodel, cvec, dtauvecdJ, &
            mode1=m1, mode2=m2, vemopt=vemopt )
        else
          call dstressdJ_viscoelastic_3D ( vemodel, cvec, dtauvecdJ, &
            mode1=m1, mode2=m2, vemopt=vemopt )
        end if

      end if

      dtautendJ(:,1,1) = dtauvecdJ(:,1)
      dtautendJ(:,1,2) = dtauvecdJ(:,2)
      dtautendJ(:,1,3) = dtauvecdJ(:,3)
      dtautendJ(:,2,1) = dtauvecdJ(:,2)
      dtautendJ(:,2,2) = dtauvecdJ(:,4)
      dtautendJ(:,2,3) = dtauvecdJ(:,5)
      dtautendJ(:,3,1) = dtauvecdJ(:,3)
      dtautendJ(:,3,2) = dtauvecdJ(:,5)
      dtautendJ(:,3,3) = dtauvecdJ(:,6)

    end if

  end subroutine dstressdJ_tensor_viscoelastic


! right-hand side of viscoelastic model
! NOTE: this is not a general usable subroutine, but more like an inlined
!       procedure. Most in/out variables are in a global module and need to
!       be defined/allocated.

  subroutine rhs_viscoelastic ( gvec, cvec, rhs, vemopt )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      coorsys, vel3D, varpar, logc, bten, vemodel, mode1, mode2, mvemodel

    real(dp), dimension(:,:), intent(in) :: gvec
    real(dp), dimension(:,:,:), intent(in) :: cvec
    real(dp), dimension(:,:,:), intent(out) :: rhs
    type(vemopt_t), optional, intent(inout) :: vemopt

    real(dp), dimension(size(gvec,1),3,3) :: gradv

    if ( coorsys <= 1 .and. vel3D == 0 ) then

!     2D and axisymmetric

      if ( varpar ) then

!       variable coefficients

        if ( logc ) then
          call rhs_viscoelastic_2D_log ( vemodel, gvec, cvec, rhs, &
            mode1, mode2, mvemodel=mvemodel, vemopt=vemopt )
        else if ( bten ) then
          call rhs_viscoelastic_2D_b ( vemodel, gvec, cvec, rhs, &
            mode1, mode2, mvemodel=mvemodel, vemopt=vemopt )
        else
          call rhs_viscoelastic_2D ( vemodel, gvec, cvec, rhs, &
            mode1, mode2, mvemodel=mvemodel, vemopt=vemopt )
        end if

      else

!       fixed coefficients

        if ( logc ) then
          call rhs_viscoelastic_2D_log ( vemodel, gvec, cvec, rhs, &
            mode1, mode2, vemopt=vemopt )
        else if ( bten ) then
          call rhs_viscoelastic_2D_b ( vemodel, gvec, cvec, rhs, &
            mode1, mode2, vemopt=vemopt )
        else
          call rhs_viscoelastic_2D ( vemodel, gvec, cvec, rhs, &
            mode1, mode2, vemopt=vemopt )
        end if

      end if

    else if ( coorsys == 2 .or. vel3D == 1 ) then

!     3D

      gradv(:,1,1) = gvec(:,1)
      gradv(:,1,2) = gvec(:,2)
      gradv(:,1,3) = gvec(:,3)
      gradv(:,2,1) = gvec(:,4)
      gradv(:,2,2) = gvec(:,5)
      gradv(:,2,3) = gvec(:,6)
      gradv(:,3,1) = gvec(:,7)
      gradv(:,3,2) = gvec(:,8)
      gradv(:,3,3) = gvec(:,9)

      if ( varpar ) then

!       variable coefficients

        if ( logc ) then
          call rhs_viscoelastic_3D_log ( vemodel, gradv, cvec, rhs, &
            mode1, mode2, mvemodel=mvemodel, vemopt=vemopt )
        else if ( bten ) then
          call rhs_viscoelastic_3D_b ( vemodel, gradv, cvec, rhs, &
            mode1, mode2, mvemodel=mvemodel, vemopt=vemopt )
        else
          call rhs_viscoelastic_3D ( vemodel, gradv, cvec, rhs, &
            mode1, mode2, mvemodel=mvemodel, vemopt=vemopt )
        end if

      else

!       fixed coefficients

        if ( logc ) then
          call rhs_viscoelastic_3D_log ( vemodel, gradv, cvec, rhs, &
            mode1, mode2, vemopt=vemopt )
        else if ( bten ) then
          call rhs_viscoelastic_3D_b ( vemodel, gradv, cvec, rhs, &
            mode1, mode2, vemopt=vemopt )
        else
          call rhs_viscoelastic_3D ( vemodel, gradv, cvec, rhs, &
            mode1, mode2, vemopt=vemopt )
        end if

      end if

    end if

  end subroutine rhs_viscoelastic


! evaluate the modulus in mvemodel

  subroutine evaluate_modulus_mve ( mesh, problem, elem, elgrp, &
    coefficients, oldvectors )

    use viscoelastic_globals_m

    type(problem_t), intent(in) :: problem
    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    integer, intent(in) :: elem, elgrp

    integer :: ip, i, np
    real(dp) :: tmpar(size(mvemodel)), Gnod(size(phi,2))

    np = size(mvemodel)

    select case ( coefficients%i(68) )

      case(0) ! modulus is constant

        do ip = 1, np
          mvemodel(ip)%modulus(mode1:mode2) = vemodel%modulus(mode1:mode2)
        end do

      case(1) ! modulus is given by function

        do i = mode1, mode2
          do ip = 1, np
            mvemodel(ip)%modulus(i) = coefficients%func1(1)%p ( i, x(ip,:) )
          end do
        end do

      case(2) ! modulus is given by nodal point values

        do i = mode1, mode2
          call get_vector ( mesh, oldvectors%p(1)%p, &
            oldvectors%v1(1)%p(i), elgrp, elem, Gnod, &
            layer=coefficients%i(38) )
          do ip = 1, np
            mvemodel(ip)%modulus(i) = dot_product ( phi(ip,:), Gnod )
          end do
        end do

      case(3) ! modulus is given by values per element

        do i = mode1, mode2
          call get_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, &
            r1=tmpar, nr=i )
          do ip = 1, np
            mvemodel(ip)%modulus(i) = tmpar(ip)
          end do
        end do

      case default

        call errormsg_case_default ( 'evaluate_modulus_mve', &
          'coefficients%i(68)', int_value=coefficients%i(68) )

    end select

  end subroutine evaluate_modulus_mve


! evaluate the relaxation time in mvemodel

  subroutine evaluate_lambda_mve ( mesh, problem, elem, elgrp, &
    coefficients, oldvectors )

    use viscoelastic_globals_m

    type(problem_t), intent(in) :: problem
    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    integer, intent(in) :: elem, elgrp

    integer :: ip, i, np
    real(dp) :: tmpar(size(mvemodel)), lambdanod(size(phi,2))

    np = size(mvemodel)

    select case ( coefficients%i(69) )

      case(0) ! relaxation time is constant

        do ip = 1, np
          mvemodel(ip)%lambda(mode1:mode2) = vemodel%lambda(mode1:mode2)
        end do

      case(1) ! relaxation time is given by function

        do i = mode1, mode2
          do ip = 1, np
            mvemodel(ip)%lambda(i) = coefficients%func1(2)%p ( i, x(ip,:) )
          end do
        end do

      case(2) ! relaxation time is given by nodal point values

        do i = mode1, mode2
          call get_vector ( mesh, oldvectors%p(1)%p, &
            oldvectors%v1(2)%p(i), elgrp, elem, lambdanod, &
            layer=coefficients%i(38) )
          do ip = 1, np
            mvemodel(ip)%lambda(i) = dot_product ( phi(ip,:), lambdanod )
          end do
        end do

      case(3) ! relaxation time is given by values per element

        do i = mode1, mode2
          call get_elvector ( mesh, oldvectors%e(2)%p, elgrp, elem, &
            r1=tmpar, nr=i )
          do ip = 1, np
            mvemodel(ip)%lambda(i) = tmpar(ip)
          end do
        end do

      case default

        call errormsg_case_default ( 'evaluate_lambda_mve', &
          'coefficients%i(69)', int_value=coefficients%i(69) )

    end select

  end subroutine evaluate_lambda_mve


! evaluate the non-linear parameters in mvemodel

  subroutine evaluate_nonlin_mve ( mesh, problem, elem, elgrp, &
    coefficients, oldvectors )

    use viscoelastic_globals_m

    type(problem_t), intent(in) :: problem
    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    integer, intent(in) :: elem, elgrp

    integer :: ip, j, numnlnpar, i, np
    real(dp) :: tmpar(size(mvemodel)), nonlinnod(size(phi,2))

    np = size(mvemodel)

!   get the number of non-linear parameters

    numnlnpar = size(vemodel%nonlin,1)

    if ( numnlnpar == 0 ) return

    select case ( coefficients%i(70) )

      case(0) ! non-linear parameters are constant

        do ip = 1, np
          do j = 1, numnlnpar
            mvemodel(ip)%nonlin(j,mode1:mode2) = vemodel%nonlin(j,mode1:mode2)
          end do
        end do

      case(1) ! non-linear parameters given by function

        do i = mode1, mode2
          do j = 1, numnlnpar
            do ip = 1, np
              mvemodel(ip)%nonlin(j,i) = coefficients%func1(2+j)%p( i, x(ip,:) )
            end do
          end do
        end do

      case(2) ! non-linear parameters given by nodal point values

        do i = mode1, mode2
          do j = 1, numnlnpar
            call get_vector ( mesh, oldvectors%p(1)%p, &
              oldvectors%v1(2+j)%p(i), elgrp, elem, nonlinnod, &
              layer=coefficients%i(38) )
            do ip = 1, np
              mvemodel(ip)%nonlin(j,i) = dot_product ( phi(ip,:), nonlinnod )
            end do
          end do
        end do

      case(3) ! non-linear parameters given by values per element

        do i = mode1, mode2
          do j = 1, numnlnpar
            call get_elvector ( mesh, oldvectors%e(2+j)%p, elgrp, elem, &
              r1=tmpar, nr=i )
            do ip = 1, np
              mvemodel(ip)%nonlin(j,i) = tmpar(ip)
            end do
          end do
        end do

      case default

        call errormsg_case_default ( 'evaluate_nonlin_mve', &
          'coefficients%i(70)', int_value=coefficients%i(70) )

    end select

  end subroutine evaluate_nonlin_mve


! evaluate the adapted lambda parameters in mvemodel

  subroutine evaluate_alam_mve ( mesh, problem, elem, elgrp, &
    coefficients, oldvectors )

    use viscoelastic_globals_m

    type(problem_t), intent(in) :: problem
    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    integer, intent(in) :: elem, elgrp

    integer :: ip, j, numnlnpar, numalampar, i, np
    real(dp) :: tmpar(size(mvemodel)), alamnod(size(phi,2))

    np = size(mvemodel)
    numnlnpar = size(vemodel%nonlin,1)

!   get the number of adapted lambda parameters

    numalampar = size(vemodel%alam,1)

    if ( numalampar == 0 ) return

    select case ( coefficients%i(81) )

      case(0) ! adapted lambda parameters are constant

        do ip = 1, np
          do j = 1, numalampar
            mvemodel(ip)%alam(j,mode1:mode2) = vemodel%alam(j,mode1:mode2)
          end do
        end do

      case(1) ! adapted lambda parameters given by function

        do i = mode1, mode2
          do j = 1, numalampar
            do ip = 1, np
              mvemodel(ip)%alam(j,i) = &
                  coefficients%func1(2+numnlnpar+j)%p( i, x(ip,:) )
            end do
          end do
        end do

      case(2) ! adapted lambda parameters given by nodal point values

        do i = mode1, mode2
          do j = 1, numalampar
            call get_vector ( mesh, oldvectors%p(1)%p, &
              oldvectors%v1(2+numnlnpar+j)%p(i), elgrp, elem, alamnod, &
              layer=coefficients%i(38) )
            do ip = 1, np
              mvemodel(ip)%alam(j,i) = dot_product ( phi(ip,:), alamnod )
            end do
          end do
        end do

      case(3) ! adapted lambda parameters given by values per element

        do i = mode1, mode2
          do j = 1, numalampar
            call get_elvector ( mesh, oldvectors%e(2+numnlnpar+j)%p, &
              elgrp, elem, r1=tmpar, nr=i )
            do ip = 1, np
              mvemodel(ip)%alam(j,i) = tmpar(ip)
            end do
          end do
        end do

      case default

        call errormsg_case_default ( 'evaluate_alam_mve', &
          'coefficients%i(81)', int_value=coefficients%i(81) )

    end select

  end subroutine evaluate_alam_mve


! master routine for evaluate mvemodel

  subroutine evaluate_mve ( mesh, problem, elem, elgrp, coefficients, &
    oldvectors, modulus, lambda, nonlin, alam )

    type(problem_t), intent(in) :: problem
    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    integer, intent(in) :: elem, elgrp
!   if present indicate the filling of variable material parameters in mvemodel
!   default=.false.
    logical, optional, intent(in) :: modulus, lambda, nonlin, alam

    logical :: lmodulus, llambda, lnonlin, lalam

    lmodulus = set_optional ( variable=modulus, default=.false. )
    llambda = set_optional ( variable=lambda, default=.false. )
    lnonlin = set_optional ( variable=nonlin, default=.false. )
    lalam = set_optional ( variable=alam, default=.false. )

    if ( lmodulus ) then

!     determine the relaxation time in mvemodel

      call evaluate_modulus_mve ( mesh, problem, elem, elgrp, coefficients, &
        oldvectors )

    end if

    if ( llambda ) then

!     determine the relaxation time in mvemodel

      call evaluate_lambda_mve ( mesh, problem, elem, elgrp, coefficients, &
        oldvectors )

    end if

    if ( lnonlin ) then

!     determine the non-linear parameters in mvemodel

      call evaluate_nonlin_mve ( mesh, problem, elem, elgrp, coefficients, &
        oldvectors )

    end if

    if ( lalam ) then

!     determine the adapted lambda parameters in mvemodel

      call evaluate_alam_mve ( mesh, problem, elem, elgrp, coefficients, &
        oldvectors )

    end if

  end subroutine evaluate_mve


! get conformation tensor for a single mode

! NOTE: this is not a general usable subroutine, but more like an inlined
!       procedure. Some in/out variables are in a global module and need to
!       be defined and/or allocated.

  subroutine get_conformation ( mesh, problem, elgrp, elem, oldvectors, mode, &
   cmode, isv, cst )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      cstorage, layer, ndfc, physqc

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    type(oldvectors_t), intent(in) :: oldvectors

!   mode number
    integer, intent(in) :: mode

!   conformation for a single mode
    real(dp), dimension(:,:), intent(out) :: cmode

!   index in the oldvectors for the sysvector
    integer, optional, intent(in) :: isv

!   override the cstorage as specified in the globals module
    integer, optional, intent(in) :: cst

    integer :: isvl, j, cstl, ncmp, posc(ndfc)
    real(dp), allocatable, dimension(:) :: cm


    isvl = set_optional ( variable=isv, default=1 )
    cstl = set_optional ( variable=cst, default=cstorage )
    ncmp = size(cmode,2)
    if ( cstl > 1 ) allocate(cm(ndfc*ncmp))

    select case ( cstl )

    case(1) ! each component in separate sysvector

      call get_sysvector ( mesh, problem, oldvectors%s2(isvl)%p(1,mode), &
        elgrp, elem, cmode(:,1), posu=posc, layer=layer )

      do j = 2, ncmp
        cmode(:,j) = oldvectors%s2(isvl)%p(j,mode)%u(posc)
      end do

    case(2) ! each mode in separate sysvector

      call get_sysvector ( mesh, problem, oldvectors%s1(isvl)%p(mode), &
        elgrp, elem, cm, layer=layer )

      cmode = reshape ( cm, [ ndfc, ncmp ] )

    case(3) ! all components/all modes in sysvector, physical quantites

      call get_sysvector ( mesh, problem, oldvectors%s(isvl)%p, &
        elgrp, elem, cm, physq=[physqc+mode-1], layer=layer )

      cmode = reshape ( cm, [ ndfc, ncmp ] )

    case default

      call errormsg_case_default ( 'get_conformation', &
        'cstl', int_value=cstl )

    end select

  end subroutine get_conformation


! get conformation tensor for a single mode on a geometry

! NOTE: this is not a general usable subroutine, but more like an inlined
!       procedure. Some in/out variables are in a global module and need to
!       be defined and/or allocated.

  subroutine get_conformation_geometry ( mesh, problem, elem, oldvectors, &
    curve, surface, volume, ndimr, geometry, mode, cmode, isv, cst )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      cstorage, layer, ndfc, physqc

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elem
    type(oldvectors_t), intent(in) :: oldvectors

!   the geometry number: one of curve, surface, volume must be present (legacy)
    integer, intent(in), optional :: curve, surface, volume

!   dimension of reference space (ndimr) and geometry number
    integer, intent(in), optional :: ndimr, geometry

!   mode number
    integer, intent(in) :: mode

!   conformation for a single mode
    real(dp), dimension(:,:), intent(out) :: cmode

!   index in the oldvectors for the sysvector
    integer, optional, intent(in) :: isv

!   override the cstorage as specified in the globals module
    integer, optional, intent(in) :: cst

    integer :: isvl, j, cstl, ncmp, posc(ndfc)
    real(dp), allocatable, dimension(:) :: cm


    isvl = set_optional ( variable=isv, default=1 )
    cstl = set_optional ( variable=cst, default=cstorage )
    ncmp = size(cmode,2)
    if ( cstl > 1 ) allocate(cm(ndfc*ncmp))

    select case ( cstl )

    case(1) ! each component in separate sysvector

      call get_sysvector_geometry ( mesh, problem, &
        oldvectors%s2(isvl)%p(1,mode), elem, cmode(:,1), &
        curve, surface, volume, ndimr, geometry, &
        posu=posc, layer=layer )

      do j = 2, ncmp
        cmode(:,j) = oldvectors%s2(isvl)%p(j,mode)%u(posc)
      end do

    case(2) ! each mode in separate sysvector

      call get_sysvector_geometry ( mesh, problem, &
        oldvectors%s1(isvl)%p(mode), elem, cm, &
        curve, surface, volume, ndimr, geometry, &
        posu=posc, layer=layer )

      cmode = reshape ( cm, [ ndfc, ncmp ] )

    case(3) ! all components/all modes in sysvector, physical quantites

      call get_sysvector_geometry ( mesh, problem, &
        oldvectors%s(isvl)%p, elem, cm, &
        curve, surface, volume, ndimr, geometry, &
        physq=[physqc+mode-1], posu=posc, layer=layer )

      cmode = reshape ( cm, [ ndfc, ncmp ] )

    case default

      call errormsg_case_default ( 'get_conformation_geometry', &
        'cstl', int_value=cstl )

    end select

  end subroutine get_conformation_geometry


! get J from the pressure in the system vector

! NOTE: this is not a general usable subroutine, but more like an inlined
!       procedure. Some in/out variables are in a global module and need to
!       be defined and/or allocated.

  subroutine get_J_from_pressure ( mesh, problem, elgrp, elem, coefficients, &
    sysvector, vemopt, vemopt_gammap )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      shapefuncp, x, xg, physqpress, layer, &
!     output:
      pr, press, &
!     input/output:
      psi

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    type(coefficients_t), intent(in) :: coefficients
    type(sysvector_t), intent(in) :: sysvector
    type(vemopt_t), optional, intent(inout) :: vemopt
    type(vemopt_gammap_t), optional, intent(inout) :: vemopt_gammap

    real(dp) :: Kmod, p0

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   get pressure and compute relative change in volume J

    call get_sysvector ( mesh, problem, sysvector, elgrp, elem, pr, &
      physq=[physqpress], layer=layer )

    press = matmul ( psi, pr )

    Kmod = coefficients%r(29)
    p0 = coefficients%r(30)

    if ( present(vemopt) ) then
      vemopt%J = exp( -(press - p0) / Kmod )
    else if ( present(vemopt_gammap) ) then
      vemopt_gammap%J = exp( -(press - p0) / Kmod )
    end if

  end subroutine get_J_from_pressure


! get J from the pressure in the system vector on a geometry

! NOTE: this is not a general usable subroutine, but more like an inlined
!       procedure. Some in/out variables are in a global module and need to
!       be defined and/or allocated.

  subroutine get_J_from_pressure_geometry ( mesh, problem, elem, coefficients, &
    sysvector, vemopt, curve, surface, volume, ndimr, geometry )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      psi, physqpress, layer, &
!     output:
      pr, press

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elem
    type(coefficients_t), intent(in) :: coefficients
    type(sysvector_t), intent(in) :: sysvector
    type(vemopt_t), intent(inout) :: vemopt

!   the geometry number: one of curve, surface, volume must be present (legacy)
    integer, intent(in), optional :: curve, surface, volume

!   dimension of reference space (ndimr) and geometry number
    integer, intent(in), optional :: ndimr, geometry

    real(dp) :: Kmod, p0

!   get pressure and compute relative change in volume J

    call get_sysvector_geometry ( mesh, problem, sysvector, elem, pr, &
      curve, surface, volume, ndimr, geometry, physq=[physqpress], &
      layer=layer )

    press = matmul ( psi, pr )

    Kmod = coefficients%r(29)
    p0 = coefficients%r(30)

    vemopt%J = exp( -(press - p0) / Kmod )

  end subroutine get_J_from_pressure_geometry


! get the equivalent plastic strain

! NOTE: this is not a general usable subroutine, but more like an inlined
!       procedure. Some in/out variables are in a global module and need to
!       be defined and/or allocated.

  subroutine get_gammap ( mesh, problem, elgrp, elem, oldvectors, &
   gp, isv, gpst )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      gpstorage, layer, physqgammap

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    type(oldvectors_t), intent(in) :: oldvectors

!   equivalent plastic strain
    real(dp), dimension(:), intent(out) :: gp

!   index in the oldvectors for the sysvector
    integer, optional, intent(in) :: isv

!   override the gpstorage as specified in the globals module
    integer, optional, intent(in) :: gpst

    integer :: isvl, gpstl

    isvl = set_optional ( variable=isv, default=1 )
    gpstl = set_optional ( variable=gpst, default=gpstorage )

    select case ( gpstl )

    case(1) ! stored in a separate sysvector

      call get_sysvector ( mesh, problem, oldvectors%s(isvl)%p, &
        elgrp, elem, gp, layer=layer )

    case(2) ! multiple physical quantities in sysvector

      call get_sysvector ( mesh, problem, oldvectors%s(isvl)%p, &
        elgrp, elem, gp, physq=[physqgammap], layer=layer )

    case default

      call errormsg_case_default ( 'get_gammap', &
        'gpstl', int_value=gpstl )

    end select

  end subroutine get_gammap


! get conformation tensor for a single mode on a geometry

! NOTE: this is not a general usable subroutine, but more like an inlined
!       procedure. Some in/out variables are in a global module and need to
!       be defined and/or allocated.

  subroutine get_gammap_geometry ( mesh, problem, elem, oldvectors, &
    curve, surface, volume, ndimr, geometry, gp, isv, gpst )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      gpstorage, layer, physqgammap

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elem
    type(oldvectors_t), intent(in) :: oldvectors

!   the geometry number: one of curve, surface, volume must be present (legacy)
    integer, intent(in), optional :: curve, surface, volume

!   dimension of reference space (ndimr) and geometry number
    integer, intent(in), optional :: ndimr, geometry

!   equivalent plastic strain
    real(dp), dimension(:), intent(out) :: gp

!   index in the oldvectors for the sysvector
    integer, optional, intent(in) :: isv

!   override the gpstorage as specified in the globals module
    integer, optional, intent(in) :: gpst

    integer :: isvl, gpstl

    isvl = set_optional ( variable=isv, default=1 )
    gpstl = set_optional ( variable=gpst, default=gpstorage )

    select case ( gpstl )

    case(1) ! stored in a separate sysvector

      call get_sysvector_geometry ( mesh, problem, oldvectors%s(isvl)%p, &
        elem, gp, curve, surface, volume, ndimr, geometry, &
        layer=layer )

    case(2) ! multiple physical quantities in sysvector

      call get_sysvector_geometry ( mesh, problem, oldvectors%s(isvl)%p, &
        elem, gp, curve, surface, volume, ndimr, geometry, &
        physq=[physqgammap], layer=layer )

    case default

      call errormsg_case_default ( 'get_gammap_geometry', &
        'gpstl', int_value=gpstl )

    end select

  end subroutine get_gammap_geometry


! Right-hand side for plastic strain evolution
! NOTE: this is not a general usable subroutine, but more like an inlined
!       procedure. Most in/out variables are in a global module and need to
!       be defined/allocated.

  subroutine rhs_plastic_strain ( cvec, gammap, rhs_gammap, vemopt_gammap )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      coorsys, vel3D, varpar, logc, bten, vemodel, mode1, mode2, mvemodel

    real(dp), dimension(:,:,:), intent(in) :: cvec
    real(dp), dimension(:), intent(in) :: gammap
    real(dp), dimension(:), intent(out) :: rhs_gammap
    type(vemopt_gammap_t), optional, intent(inout) :: vemopt_gammap

    if ( coorsys <= 1 .and. vel3D == 0 ) then

!     2D and axisymmetric

      if ( varpar ) then

!       variable coefficients

        if ( logc ) then
          call rhs_plastic_strain_2D_log ( vemodel, cvec, gammap, rhs_gammap, &
            mode1, mode2, mvemodel=mvemodel, vemopt_gammap=vemopt_gammap )
        else if ( bten ) then
          call rhs_plastic_strain_2D_b ( vemodel, cvec, gammap, rhs_gammap, &
            mode1, mode2, mvemodel=mvemodel, vemopt_gammap=vemopt_gammap )
        else
          call rhs_plastic_strain_2D ( vemodel, cvec, gammap, rhs_gammap, &
            mode1, mode2, mvemodel=mvemodel, vemopt_gammap=vemopt_gammap )
        end if

      else

!       fixed coefficients

        if ( logc ) then
          call rhs_plastic_strain_2D_log ( vemodel, cvec, gammap, rhs_gammap, &
            mode1, mode2, vemopt_gammap=vemopt_gammap )
        else if ( bten ) then
          call rhs_plastic_strain_2D_b ( vemodel, cvec, gammap, rhs_gammap, &
            mode1, mode2, vemopt_gammap=vemopt_gammap )
        else
          call rhs_plastic_strain_2D ( vemodel, cvec, gammap, rhs_gammap, &
            mode1, mode2, vemopt_gammap=vemopt_gammap )
        end if

      end if

    else if ( coorsys == 2 .or. vel3D == 1 ) then

!     3D

      if ( varpar ) then

!       variable coefficients

        if ( logc ) then
          call rhs_plastic_strain_3D_log ( vemodel, cvec, gammap, rhs_gammap, &
            mode1, mode2, mvemodel=mvemodel, vemopt_gammap=vemopt_gammap )
        else if ( bten ) then
          call rhs_plastic_strain_3D_b ( vemodel, cvec, gammap, rhs_gammap, &
            mode1, mode2, mvemodel=mvemodel, vemopt_gammap=vemopt_gammap )
        else
          call rhs_plastic_strain_3D ( vemodel, cvec, gammap, rhs_gammap, &
            mode1, mode2, mvemodel=mvemodel, vemopt_gammap=vemopt_gammap )
        end if

      else

!       fixed coefficients

        if ( logc ) then
          call rhs_plastic_strain_3D_log ( vemodel, cvec, gammap, rhs_gammap, &
            mode1, mode2, vemopt_gammap=vemopt_gammap )
        else if ( bten ) then
          call rhs_plastic_strain_3D_b ( vemodel, cvec, gammap, rhs_gammap, &
            mode1, mode2, vemopt_gammap=vemopt_gammap )
        else
          call rhs_plastic_strain_3D ( vemodel, cvec, gammap, rhs_gammap, &
            mode1, mode2, vemopt_gammap=vemopt_gammap )
        end if

      end if

    end if

  end subroutine rhs_plastic_strain


! set global parameters viscoelastic (internal element)

  subroutine set_globals_viscoelastic ( coefficients, maxvar, singlemode, &
    gammap )

    use viscoelastic_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   the maximum value for coefficients%i(73) (varying coefficients) allowed.
!   default=0, i.e. constant coefficients
    integer, intent(in), optional :: maxvar

!   allow a single mode to be set by coefficients%i(85)
!   default=.false.
    logical, intent(in), optional :: singlemode

!   the use of equivalent plastic strain (gammap) has been implemented
!   default=.false.
    logical, intent(in), optional :: gammap

    integer :: lmaxvar
    logical :: lsinglemode, lgammap

    lmaxvar = set_optional ( variable=maxvar, default=0 )
    lsinglemode = set_optional ( variable=singlemode, default=.false. )
    lgammap = set_optional ( variable=gammap, default=.false. )

    nmodes = vemodel%nmodes
    ncompc = vemodel%ncompc
    ncompt = vemodel%ncompt
    vemcompressible = vemodel%compressible
    vemgammap = vemodel%gammap

    physqc = coefficients%i(83)
    physqgammap = coefficients%i(97)

    call check ( coefficients, 'set_globals_viscoelastic', &
      indexarray=[21,26,27,29,31,49,57,58,60,61,71,72,73,78,79,84,92], &
      minimum=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], &
      maximum=[1,nmodes,nmodes,3,4,1,1,size(vemodel%nonlin,1),1,1,4,1,&
               lmaxvar,1,2,3,1] )

    if ( lsinglemode ) then
      call check ( coefficients, 'set_globals_viscoelastic', &
        indexarray=[85], minimum=[0], maximum=[nmodes] )
    end if

    if ( vemgammap .and. .not. lgammap ) then
      write(*,'(/2(a/))') 'Error in set_globals_viscoelastic:', &
        ' equivalent plastic strain (gammap) has not been implemented.'
      stop
    end if

    logc = coefficients%i(21) == 1
    exps = logc .and. coefficients%i(49) == 1

    bten = coefficients%i(71) >= 1
    cproj = bten .and. coefficients%i(72) == 1

    if ( logc .and. bten ) then
      write(*,'(/2(a/))') 'Error in set_globals_viscoelastic:', &
        ' both b-tensor and log c formulation not possible.'
      stop
    end if

    if ( bten ) then
      ncomp = vemodel%ncompb   ! convected c will contain the b-tensor
    else
      ncomp = vemodel%ncompc   ! convected c will contain the c-tensor
    end if

    c_direct = get_coefficient ( coefficients, index=79, default=1 ) == 2

!   set more integer parameters

    if ( lsinglemode .and. coefficients%i(85) > 0 ) then
!     Set a single mode number
      mode1 = coefficients%i(85)
      mode2 = mode1  ! single mode
    else
      mode1 = get_coefficient ( coefficients, index=26, default=1 )
      mode2 = get_coefficient ( coefficients, index=27, default=nmodes )
    end if
    nsubm = max(mode2-mode1+1,0)

    htype     = get_coefficient ( coefficients, index=29, default=2 )
    Uscaling  = get_coefficient ( coefficients, index=31, default=3 )

    if ( Uscaling == 1 ) then
      write(*,'(/2(a/))') 'Error in set_globals_viscoelastic:', &
        ' Uscaling == 1 not available for viscoelastic elements '
      stop
    end if

!   small relaxation times

    coverl = coefficients%i(57) == 1

!   varying viscoelastic coefficients

    varpar = coefficients%i(73) == 1

!   set number of degrees of freedom for c

    intpolc = coefficients%i(12)

    shapefuncc%globalshape = globalshape
    shapefuncc%interpolation = intpolc
    shapefuncc%numbering = 'regular'

    call set_ndf ( shapefuncc, 'set_globals_viscoelastic', ndf=ndfc )

  end subroutine set_globals_viscoelastic


! set global parameters viscoelastic (boundary element)

  subroutine set_globals_viscoelastic_boun ( coefficients, maxvar, &
    singlemode, gammap )

    use viscoelastic_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   the maximum value for coefficients%i(73) (varying coefficients) al lowed.
!   default=0, i.e. constant coefficients
    integer, intent(in), optional :: maxvar

!   allow a single mode to be set by coefficients%i(85)
!   default=.false.
    logical, intent(in), optional :: singlemode

!   the use of equivalent plastic strain (gammap) has been implemented
!   default=.false.
    logical, intent(in), optional :: gammap

    integer :: lmaxvar
    logical :: lsinglemode, lgammap

    lmaxvar = set_optional ( variable=maxvar, default=0 )
    lsinglemode = set_optional ( variable=singlemode, default=.false. )
    lgammap = set_optional ( variable=gammap, default=.false. )

    nmodes = vemodel%nmodes
    ncompc = vemodel%ncompc
    ncompt = vemodel%ncompt
    vemcompressible = vemodel%compressible
    vemgammap = vemodel%gammap

    physqc = coefficients%i(83)
    physqgammap = coefficients%i(97)

    call check ( coefficients, 'set_globals_viscoelastic_boun', &
      indexarray=[21,26,27,49,57,58,60,61,71,72,73,79,84], &
      minimum=[0,0,0,0,0,0,0,0,0,0,0,0,0], &
      maximum=[1,nmodes,nmodes,1,1,size(vemodel%nonlin,1),1,1,4,1,lmaxvar,2,3] )

    if ( lsinglemode ) then
      call check ( coefficients, 'set_globals_viscoelastic_boun', &
        indexarray=[85], minimum=[0], maximum=[nmodes] )
    end if

    if ( vemgammap .and. .not. lgammap ) then
      write(*,'(/2(a/))') 'Error in set_globals_viscoelastic_boun:', &
        ' equivalent plastic strain (gammap) has not been implemented.'
      stop
    end if

    logc = coefficients%i(21) == 1
    exps = logc .and. coefficients%i(49) == 1

    bten = coefficients%i(71) >= 1
    cproj = bten .and. coefficients%i(72) == 1

    if ( logc .and. bten ) then
      write(*,'(/2(a/))') 'Error in set_globals_viscoelastic:', &
        ' both b-tensor and log c formulation not possible.'
      stop
    end if

    if ( bten ) then
      ncomp = vemodel%ncompb   ! convected c will contain the b-tensor
    else
      ncomp = vemodel%ncompc   ! convected c will contain the c-tensor
    end if

    c_direct = get_coefficient ( coefficients, index=79, default=1 ) == 2

!   set modes to be computed

    if ( lsinglemode .and. coefficients%i(85) > 0 ) then
!     Set a single mode number
      mode1 = coefficients%i(85)
      mode2 = mode1  ! single mode
    else
      mode1 = get_coefficient ( coefficients, index=26, default=1 )
      mode2 = get_coefficient ( coefficients, index=27, default=nmodes )
    end if
    nsubm = max(mode2-mode1+1,0)

!   small relaxation times

    coverl = coefficients%i(57) == 1

!   varying viscoelastic coefficients

    varpar = coefficients%i(73) == 1

!   set number of degrees of freedom for c

    intpolc = coefficients%i(12)

    shapefuncc%globalshape = globalshape
    shapefuncc%interpolation = intpolc
    shapefuncc%numbering = 'regular'

    call set_ndf ( shapefuncc, 'set_globals_viscoelastic_boun', ndf=ndfc )

  end subroutine set_globals_viscoelastic_boun


! set viscoelastic model

  subroutine set_viscoelastic_model ( coefficients )

    use viscoelastic_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   set viscoelastic model (don't forget to delete vemodel again)

    integer :: smp, emp, np_mode, nq_mode, nr_mode, npar_mode, npar_tot
    integer :: flowtype, ns_mode

    if ( coorsys <= 1 .and. vel3D == 0 ) then
      flowtype = coorsys
    else if ( coorsys == 2 .or. vel3D == 1 ) then
      flowtype = 2 ! use 3D model
    end if

!   first create model structure vemodel

    call create_viscoelastic_model ( model=coefficients%i(18), &
      vemodel=vemodel, nmodes=coefficients%i(19), flowtype=flowtype, &
      bvariant=coefficients%i(71), alam_model=coefficients%i(80), &
      vmglobal=coefficients%i(92)==1, alamJ_model=coefficients%i(96), &
      norlxmode=coefficients%i(93), alam_gammap_model=coefficients%i(98), &
      gammap_mode=coefficients%i(99) )

!   fill material parameters from the coefficients

    np_mode = size(vemodel%nonlin,1)
    nq_mode = size(vemodel%alam,1)
    nr_mode = size(vemodel%alamJ,1)
    ns_mode = size(vemodel%alam_gammap,1)
    npar_mode = 2 + np_mode + nq_mode + nr_mode + ns_mode
    npar_tot  = vemodel%nmodes * npar_mode

!   material parameters

    if ( coefficients%i(20) > 0 ) then

!     material parameters in real array r

      call parameters_in_r

    else if ( coefficients%i(20) < 0 ) then

!     material parameters in 2D real array ra2

      call parameters_in_ra2

    else

      write(*,'(/3(a/))') &
        'Error set_viscoelastic_model:', &
        ' The start index of the viscoelastic material parameters ', &
        ' coefficients%i(20) has not been set.'
      stop

    end if

  contains

!   material parameters in real array r

    subroutine parameters_in_r

      integer :: k

!     start pointer

      smp = coefficients%i(20)   ! start of material parameters

      if ( smp < 501 .or. smp > 1000 ) then
        write(*,'(/a/a,i0/,4(a/))') &
          'Error set_viscoelastic_model:', &
          ' The start index of the viscoelastic material parameters = ', smp, &
          ' This is outside the free range of 501:1000', &
          ' To avoid any clashes with data from other modules, it is', &
          ' required to move the parameters either to the free range, ', &
          ' or use a 2D array in coefficients%ra2 (use coefficients%i(20)<0).'
        stop
      end if

      emp = smp + npar_tot - 1   ! end of material parameters

!     check size of real array

      call check ( coefficients, 'set_viscoelastic_model', ncoefr=emp )

      if ( emp > 1000 ) then
        write(*,'(/a/a,i0/,4(a/))') &
          'Error set_viscoelastic_model:', &
          ' The end index of the viscoelastic material parameters = ', emp, &
          ' This is beyond the free range of 501:1000', &
          ' To avoid any clashes with data from other modules, it is', &
          ' required to the parameters either fully within the free range, ', &
          ' or use a 2D array in coefficients%ra2 (use coefficients%i(20)<0).'
        stop
      end if

      vemodel%modulus = coefficients%r(smp:emp:npar_mode)
      vemodel%lambda = coefficients%r(smp+1:emp:npar_mode)
      do k = 1, np_mode
        vemodel%nonlin(k,:) = coefficients%r(smp+1+k:emp:npar_mode)
      end do
      do k = 1, nq_mode
        vemodel%alam(k,:) = coefficients%r(smp+1+np_mode+k:emp:npar_mode)
      end do
      do k = 1, nr_mode
        vemodel%alamJ(k,:) = &
                  coefficients%r(smp+1+np_mode+nq_mode+k:emp:npar_mode)
      end do
      do k = 1, ns_mode
        vemodel%alam_gammap(k,:) = &
                  coefficients%r(smp+1+np_mode+nq_mode+nr_mode+k:emp:npar_mode)
      end do

    end subroutine parameters_in_r


!   material parameters in 2D real array ra2

    subroutine parameters_in_ra2

      integer :: m, ib, ie, sh(2)

!     start pointer

      m = - coefficients%i(20)   ! start of material parameters

      if ( m > size(coefficients%ra2) ) then
        write(*,'(/a/a,i0/a,i0/)') &
          'Error set_viscoelastic_model:', &
          ' The index for the viscoelastic material parameters = ', m, &
          ' which is larger than the size of coefficients%ra2 = ', &
          size(coefficients%ra2)
        stop
      else
        sh = shape(coefficients%ra2(m)%a)
        if ( any( sh /= [npar_mode,vemodel%nmodes] ) ) then
          write(*,'(/a/a,2(a,i0,a,i0,a/))') &
            'Error set_viscoelastic_model:', &
            ' The shape of ra2 for the viscoelastic material parameters is', &
            ' (', sh(1), ',', sh(2), '),', &
            ' whereas it must be (', npar_mode, ',', vemodel%nmodes, ').'
          stop
        end if
      end if

      vemodel%modulus = coefficients%ra2(m)%a(1,:)
      vemodel%lambda = coefficients%ra2(m)%a(2,:)
      ib = 3; ie = 2 + np_mode
      vemodel%nonlin = coefficients%ra2(m)%a(ib:ie,:)
      ib = ie+1; ie = ie + nq_mode
      vemodel%alam = coefficients%ra2(m)%a(ib:ie,:)
      ib = ie+1; ie = ie + nr_mode
      vemodel%alamJ = coefficients%ra2(m)%a(ib:ie,:)
      ib = ie+1; ie = ie + ns_mode
      vemodel%alam_gammap = coefficients%ra2(m)%a(ib:ie,:)

    end subroutine parameters_in_ra2

  end subroutine set_viscoelastic_model


! Preamble for the viscoelastic element

  subroutine set_viscoelastic_elem ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, maxvar, singlemode, gammap )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    integer, intent(in), optional :: maxvar
    logical, intent(in), optional :: singlemode, gammap


    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

!     set globals

      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic ( coefficients, maxvar, singlemode, &
        gammap )

!     allocate arrays

      allocate ( theta(ninti,ndfc), dtheta(ninti,ndfc,ndim) )
      allocate ( dthetadx(ninti,ndfc,ndim) )

!     set shape function c tensor
      call set_shape_function ( shapefuncc, xig, theta, dtheta )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      if ( first ) then

!       first element in this group

!       set globals

        call set_viscoelastic_model ( coefficients )
        call set_globals_viscoelastic ( coefficients, maxvar, singlemode, &
          gammap )

      end if

!     allocate arrays

      allocate ( theta(ninti,ndfc), dtheta(ninti,ndfc,ndim) )
      allocate ( dthetadx(ninti,ndfc,ndim) )

!     set shape function c tensor
      call set_shape_function ( shapefuncc, xig, theta, dtheta )

    end if

  end subroutine set_viscoelastic_elem


! Unset the preamble for the viscoelastic element
! (deallocate arrays allocated in set_... )

  subroutine unset_viscoelastic_elem ( last, coefficients )

    use viscoelastic_globals_m

    logical, intent(in) :: last
    type(coefficients_t), intent(in) :: coefficients

!   deallocate memory

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call delete ( vemodel )

      deallocate ( theta, dtheta, dthetadx )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      if ( last ) then

!       last element in this group

        call delete ( vemodel )

      end if

      deallocate ( theta, dtheta, dthetadx )

    end if

  end subroutine unset_viscoelastic_elem

end module viscoelastic_elements_generic_m
