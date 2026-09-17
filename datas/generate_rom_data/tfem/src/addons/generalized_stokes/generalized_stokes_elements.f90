
! Copyright (C) 2007-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the generalized Stokes flow equation
!
!    - div( eta(nabla u+nabla u^T) ) + nabla p = f
!      div u = 0
!
! with varying viscosity eta
!

module generalized_stokes_elements_m

  use tfem_elem_m
  use stokes_elements_m
  use stokes_elements_generic_m, only: set_stokes_shape_function_global
  use generalized_newtonian_elements_m
  use ziggurat_m

  implicit none


contains


! Internal element routine for the generalized Stokes equation with varying
! coefficient eta and (optionally) varying Gauss rule per element.

  subroutine generalized_stokes_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use generalized_stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, M, i, j, ip, vfuncnr
    real(dp) :: S_NMii, S_NMij, L_NMi


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

    if ( first) &
            call set_globals_generalized_stokes ( mesh, coefficients, elgrp )

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

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if


    if ( vector ) then

      vfuncnr = coefficients%i(14)

      if ( vfuncnr > 0 ) then

        do ip = 1, ninti
          fg(ip,:) = coefficients%vfunc ( ndim, vfuncnr, xg(ip,:) )
        end do

        do j = 1, ndim
          do N = 1, ndf
            elemvec( pos(N,j) ) = sum ( fg(:,j) * phi(:,N) * detF * wg )
          end do
        end do

        elemvec( posp ) = 0

      else

        elemvec = 0

      end if

    end if


    if ( matrix ) then

!     determine viscosity coefficient eta

      call evaluate_eta ( mesh, problem, coefficients, oldvectors, elgrp, &
        elem, etag, xg, phi )

!     fill matrix

!     diagonal blocks

      do N = 1, ndf
        do M = N, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dphidx(ip,M,:) )
          end do
          do i = 1, ndim
            work1 = work + dphidx(:,N,i) * dphidx(:,M,i)
            if ( coorsys == 1 .and. i == 2 ) then ! axisymmetric
              work1 = work1 + 2 * phi(:,N) * phi(:,M) / xg(:,2) ** 2
            end if
            S_NMii = sum ( etag * work1 * detF * wg )
            elemmat( pos(N,i), pos(M,i) ) = S_NMii
            elemmat( pos(M,i), pos(N,i) ) = S_NMii  ! = S_MNii symmetry
          end do
        end do
      end do

!     off-diagonal blocks

      do N = 1, ndf
        do M = 1, ndf
          do i = 1, ndim
            do j = i+1, ndim
              work = dphidx(:,N,j) * dphidx(:,M,i)
              S_NMij = sum ( etag * work * detF * wg )
              elemmat( pos(N,i), pos(M,j) ) = S_NMij
              elemmat( pos(M,j), pos(N,i) ) = S_NMij  ! = S_MNji symmetry
            end do
          end do
        end do
      end do

!     velocity-pressure part

      do N = 1, ndfp
        do M = 1, ndf
          do i = 1, ndim
            if ( coorsys == 1 .and. i == 2 ) then ! axisymmetric
              work = dphidx(:,M,2) + phi(:,M) / xg(:,2)
              L_NMi = sum ( psi(:,N) * work * detF * wg )
            else
              L_NMi = sum ( psi(:,N) * dphidx(:,M,i) * detF * wg )
            end if
            elemmat( posp(N), pos(M,i) ) = - L_NMi
            elemmat( pos(M,i), posp(N) ) = - L_NMi  ! = - L_MNi symmetry
          end do
        end do
      end do

!     pressure-pressure part

      elemmat( posp, posp ) = 0

!     set continuity equation (div u=0) to zero

      if ( coefficients%i(39) == 1 ) then

        do i = 1, ndim
          elemmat( posp, pos(:,i) ) = 0
        end do

      end if

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

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

      allocate ( etag(ninti) )
      allocate ( work(ninti), work1(ninti) )
      allocate ( u(ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( etag )
      deallocate ( work, work1 )
      deallocate ( u )

    end subroutine deallocate_arrays

  end subroutine generalized_stokes_elem


! Internal element routine to fill the viscous dissipation term in the Gauss
! points. The result is stored in a vector defined per element (elvector).
! This element should be used together with the routine loop_over_elements.
! The resulting elvector can be used as input to other elements based on
! the energy equation.

  subroutine fill_viscous_dissipation_gauss ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors )

    use generalized_stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors

    integer :: j, ip


    call set_stokes_elem ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors, maxvel3D=1 )

    if ( first) &
            call set_globals_generalized_stokes ( mesh, coefficients, elgrp )

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

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    ugvector = matmul ( phi, uvector )

!   determine viscosity coefficient eta

    call evaluate_eta ( mesh, problem, coefficients, oldvectors, elgrp, &
      elem, etag, xg, phi )

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

!   extra stress tensor tau

    do ip = 1, ninti
      tauten(ip,:,:) = 2 * etag(ip) * Dten(ip,:,:)
    end do

!   viscous dissipation ( tau : D )

    do ip = 1, ninti
      work(ip) = sum ( tauten(ip,:,:) * Dten(ip,:,:) )
    end do

    call put_elvector ( mesh, oldvectors%e(2)%p, elgrp, elem, r1=work )


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

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
      allocate ( etag(ninti) )
      if ( coorsys <= 1 .and. vel3D == 0 ) then
        allocate ( gradu(ninti,ndim+coorsys,ndim+coorsys) )
        allocate ( Dten(ninti,ndim+coorsys,ndim+coorsys) )
        allocate ( tauten(ninti,ndim+coorsys,ndim+coorsys) )
      else
        allocate ( gradu(ninti,ncompu,ncompu) )
        allocate ( Dten(ninti,ncompu,ncompu) )
        allocate ( tauten(ninti,ncompu,ncompu) )
      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( u, ugvector )
      deallocate ( uvector )
      deallocate ( gradu )
      deallocate ( Dten )
      deallocate ( tauten )
      deallocate ( work )
      deallocate ( etag )

    end subroutine deallocate_arrays

  end subroutine fill_viscous_dissipation_gauss


! viscosity coefficient eta

  subroutine evaluate_eta ( mesh, problem, coefficients, oldvectors, elgrp, &
    elem, eta, x, phi )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    integer, intent(in) :: elgrp, elem
    real(dp), dimension(:), intent(out) :: eta
    real(dp), dimension(:,:), intent(in) :: x, phi

    call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
      elgrp, elem, choice=coefficients%i(251), value=coefficients%r(1), &
      func=coefficients%func1(1)%p, funcnr=coefficients%i(252), x=x, &
      indx_v=1, layer=coefficients%i(38), phi=phi, indx_e=1, coef=eta )

  end subroutine evaluate_eta


! Internal element routine for the right-hand side of the generalized Stokes
! equation for fluctuating hydrodynamics.
! Note, that the delta correlation in time has not yet been accounted for.

  subroutine generalized_stokes_rhs_fh ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use generalized_stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: j, ip, numstress
    real(dp) :: kT


!   set gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

    if ( first) &
            call set_globals_generalized_stokes ( mesh, coefficients, elgrp )

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

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      write(*,'(/a/a/)') 'Error in generalized_stokes_rhs_fh:', &
        ' Axisymmetric elements not applicable. '
      stop
    end if

    call shape_derivative ( dphi, Finv, dphidx )


    if ( vector ) then

!     determine viscosity coefficient eta

      call evaluate_eta ( mesh, problem, coefficients, oldvectors, elgrp, &
        elem, etag, xg, phi )

!     stochastic stress tensor: stensor

      call generate_brownf ( brownf, coefficients )

      if ( ndim == 2 ) then
!       2D
        tauten(:,1,1) = sqrt(2._dp) * brownf(:,1)
        tauten(:,1,2) = brownf(:,2)
        tauten(:,2,1) = tauten(:,1,2)
        tauten(:,2,2) = sqrt(2._dp) * brownf(:,3)
      else if ( ndim == 3 ) then
!       3D
        tauten(:,1,1) = sqrt(2._dp) * brownf(:,1)
        tauten(:,1,2) = brownf(:,2)
        tauten(:,1,3) = brownf(:,3)
        tauten(:,2,1) = tauten(:,1,2)
        tauten(:,2,2) = sqrt(2._dp) * brownf(:,4)
        tauten(:,2,3) = brownf(:,5)
        tauten(:,3,1) = tauten(:,1,3)
        tauten(:,3,2) = tauten(:,2,3)
        tauten(:,3,3) = sqrt(2._dp) * brownf(:,6)
      end if

      kT = coefficients%r(208)

      work = sqrt ( 2 * kT * etag * detF * wg )

!     - (nabla v)^T:s

      do ip = 1, ninti
        work6(ip,:,:) = matmul ( dphidx(ip,:,:), tauten(ip,:,:) )
      end do
      do j = 1, ndim
        work2(:,j) = - matmul ( work, work6(:,:,j) )
      end do
      elemvec = reshape ( work2, [ ndf*ndim ] )

    end if


    if ( matrix ) then

      write(*,'(/a/a/)') 'Error in generalized_stokes_rhs_fh:', &
        'No matrix to build. Call build_system with buildmatrix=.false.'
      stop

    end if


!   unset gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

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

      select case ( coorsys )
        case(0)
          numstress = 3
        case(2)
          numstress = 6
        case default
          call errormsg_case_default ( 'generalized_stokes_rhs_fh', &
            'coorsys', int_value=coorsys )
      end select

      allocate ( etag(ninti) )
      allocate ( work(ninti), work1(ninti) )
      allocate ( work6(ninti,ndf,ndim), work2(ndf,ndim) )
      allocate ( tauten(ninti,ndim,ndim), u(ndf) )
      allocate ( brownf(ninti,numstress) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( etag )
      deallocate ( work, work1 )
      deallocate ( work6, work2 )
      deallocate ( tauten, u )
      deallocate ( brownf )

    end subroutine deallocate_arrays

  end subroutine generalized_stokes_rhs_fh


! routine for the random vector brownf

  subroutine generate_brownf ( brownf, coefficients )

    real(dp), dimension(:,:), intent(out) :: brownf
    type(coefficients_t), intent(in) :: coefficients

    integer :: i, j, ibrown

    ibrown = coefficients%i(256)

    select case ( ibrown )

    case(1)

!     create random numbers with uniform distribution, variance 1
!     use fortran intrinsic function

      call random_number ( brownf )

      brownf = ( 2 * brownf - 1 ) * sqrt(3._dp)

    case(2)

!     create random numbers with uniform distribution, variance 1
!     use ziggurat

      do i = 1, size(brownf,1)
        do j = 1, size(brownf,2)
          brownf(i,j) = uni()
        end do
      end do

      brownf = ( 2 * brownf - 1 ) * sqrt(3._dp)

    case(3)

!     create random numbers with normal distribution, variance 1
!     use ziggurat

      do i = 1, size(brownf,1)
        do j = 1, size(brownf,2)
          brownf(i,j) = rnor()
        end do
      end do

    case default

      write(*,'(/a,i0/)') &
        'Error generate_brownf: invalid random number generator: ', ibrown
      stop

    end select

  end subroutine generate_brownf


! compute viscosity in all nodes

  subroutine generalized_stokes_viscosity ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use generalized_stokes_globals_m

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

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_generalized_stokes ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

!   determine viscosity coefficient eta in the nodes

    call evaluate_eta ( mesh, problem, coefficients, oldvectors, elgrp, &
      elem, etag, x, phi )

    elemvec = etag

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ndf), etag(nodalp) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xrnod, phi, x )
      deallocate ( u, etag )

    end subroutine deallocate_arrays

  end subroutine generalized_stokes_viscosity


! compute stresses in all nodes (from velocity gradients)

  subroutine generalized_stokes_stress ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    use generalized_stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: comp, j, maxstress


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_generalized_stokes ( mesh, coefficients, elgrp )

!     check component = coefficients%i(13)

      select case ( coorsys )
        case(0,1)
          maxstress = 3 + coorsys
        case(2)
          maxstress = 6
        case default
          call errormsg_case_default ( 'generalized_stokes_stress', &
            'coorsys', int_value=coorsys )
      end select

      call check ( coefficients, 'generalized_stokes_stress', &
        indexarray=[13], minimum=[1], maximum=[maxstress] )

!     allocate arrays

      call allocate_arrays

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ndim] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

!   determine viscosity coefficient eta in the nodes

    call evaluate_eta ( mesh, problem, coefficients, oldvectors, elgrp, &
      elem, etag, x, phi )

    comp = coefficients%i(13)

    if ( coorsys <= 1 ) then

      select case(comp)
      case(1)
        elemvec = 2 * etag * gradu(:,1,1)                ! 2 eta dudx
      case(2)
        elemvec = etag * ( gradu(:,1,2) + gradu(:,2,1) ) ! eta ( dudy + dvdx )
      case(3)
        elemvec = 2 * etag * gradu(:,2,2)                ! 2 eta dvdy
      case(4)
        where ( x(:,2) < 1e-10_dp )
          elemvec = 2 * etag * gradu(:,2,2)              ! 2 eta du_r / dr (r=0)
        elsewhere
          elemvec = 2 * etag * uvector(:,2) / x(:,2)     ! 2 eta u_r / r
        end where
      case default
        call errormsg_case_default ( 'generalized_stokes_stress', &
          'comp', int_value=comp )
      end select

    else if ( coorsys == 2 ) then

      select case(comp)
      case(1)
        elemvec = 2 * etag * gradu(:,1,1)                ! 2 eta dudx
      case(2)
        elemvec = etag * ( gradu(:,1,2) + gradu(:,2,1) ) ! eta ( dudy + dvdx )
      case(3)
        elemvec = etag * ( gradu(:,1,3) + gradu(:,3,1) ) ! eta ( dudz + dwdx )
      case(4)
        elemvec = 2 * etag * gradu(:,2,2)                ! 2 eta dvdy
      case(5)
        elemvec = etag * ( gradu(:,2,3) + gradu(:,3,2) ) ! eta ( dvdz + dwdy )
      case(6)
        elemvec = 2 * etag * gradu(:,3,3)                ! 2 eta dwdz
      case default
        call errormsg_case_default ( 'generalized_stokes_stress', &
          'comp', int_value=comp )
      end select

    end if

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ndim*ndf), uvector(ndf,ndim), gradu(nodalp,ndim,ndim) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( etag(nodalp) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( detF )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( etag )

    end subroutine deallocate_arrays

  end subroutine generalized_stokes_stress


! stokes stress tensor (for use with drag computations)

  subroutine generalized_stokes_stress_tensor ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use generalized_stokes_globals_m

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

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_generalized_stokes ( mesh, coefficients, elgrp )

!     test size of elemvec

      select case ( coorsys )
        case(0,1)
          maxstress = 3 + coorsys
        case(2)
          maxstress = 6
        case default
          call errormsg_case_default ( 'generalized_stokes_stress_tensor', &
            'coorsys', int_value=coorsys )
      end select

      if ( size(elemvec) /= nodalp*maxstress ) then
        write(*,'(/2(a/))') 'Error in generalized_stokes_stress_tensor:', &
          ' element vector has incorrect size for a tensor'
        stop
      end if

!     allocate arrays

      call allocate_arrays

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ndim] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

!   determine viscosity coefficient eta in the nodes

    call evaluate_eta ( mesh, problem, coefficients, oldvectors, elgrp, &
      elem, etag, x, phi )

    if ( coorsys <= 1 ) then

      work2(:,1) = 2 * etag * gradu(:,1,1)                ! tau_xx
      work2(:,2) = etag * ( gradu(:,1,2) + gradu(:,2,1) ) ! tau_xy
      work2(:,3) = 2 * etag * gradu(:,2,2)                ! tau_yy

      if ( coorsys == 1 ) then
!       axisymmetric
        where ( x(:,2) < 1e-10_dp )
          work2(:,4) = 2 * etag * gradu(:,2,2)           ! 2 eta du_r / dr (r=0)
        elsewhere
          work2(:,4) = 2 * etag * uvector(:,2) / x(:,2)  ! 2 eta u_r / r
        end where
      end if

    else if ( coorsys == 2 ) then

      work2(:,1) = 2 * etag * gradu(:,1,1)                ! tau_xx
      work2(:,2) = etag * ( gradu(:,1,2) + gradu(:,2,1) ) ! tau_xy
      work2(:,3) = etag * ( gradu(:,1,3) + gradu(:,3,1) ) ! tau_xz
      work2(:,4) = 2 * etag * gradu(:,2,2)                ! tau_yy
      work2(:,5) = etag * ( gradu(:,2,3) + gradu(:,3,2) ) ! tau_yz
      work2(:,6) = 2 * etag * gradu(:,3,3)                ! tau_zz

    end if

    elemvec = reshape ( work2, [nodalp*maxstress] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( u(ndim*ndf), uvector(ndf,ndim), gradu(nodalp,ndim,ndim) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( work2(nodalp,maxstress) )
      allocate ( etag(nodalp) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( detF )
      deallocate ( xrnod, phi, x )
      deallocate ( u, uvector, gradu )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )
      deallocate ( work2 )
      deallocate ( etag )

    end subroutine deallocate_arrays

  end subroutine generalized_stokes_stress_tensor


! set global parameters (internal element)

  subroutine set_globals_generalized_stokes ( mesh, coefficients, elgrp )

    use generalized_stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

!   check size of coefficients

    call check ( coefficients, 'set_globals_generalized_stokes', ncoefi=300, &
      ncoefr=250, indexarray=[251,255], minimum=[0,0], &
      maximum=[3,1] )

  end subroutine set_globals_generalized_stokes

end module generalized_stokes_elements_m

