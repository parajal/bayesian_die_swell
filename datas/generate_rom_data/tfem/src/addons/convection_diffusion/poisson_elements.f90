
! Copyright (C) 2007-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the Poisson equation:
!
!    - alpha nabla^2 u = f
!

module poisson_elements_m

  use tfem_elem_m

  implicit none


contains


! Internal element routine for the Poisson Equation
! ( with constant coefficient alpha).

  subroutine poisson_elem ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use poisson_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip, funcnr
    real(dp) :: alpha


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_poisson ( mesh, coefficients, elgrp )

      allocate ( wg(ninti), fg(ninti), detF(ninti), work(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    call isoparametric_coordinates ( x, phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( vector ) then

      funcnr = coefficients%i(12)

      if ( funcnr > 0 ) then

        do ip = 1, ninti
          fg(ip) = coefficients%func ( funcnr, xg(ip,:) )
        end do

        do i = 1, ndf
          elemvec(i) = sum ( fg * phi(:,i) * detF * wg )
        end do

      else

        elemvec = 0

      end if

    end if

    if ( matrix ) then

      alpha = coefficients%r(1)

      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = alpha * sum ( dphidx(ip,i,:) * dphidx(ip,j,:) )
          end do
          elemmat(i,j) = sum ( work * detF * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

    end if

    if ( last ) then

!     last element in this group

      deallocate ( Finv, dphidx )
      deallocate ( dphi, F )
      deallocate ( xg )
      deallocate ( xig, phi, x, work )
      deallocate ( wg, fg, detF )

    end if

  end subroutine poisson_elem


! Boundary element for a natural boundary for the Poisson equation

  subroutine poisson_natboun ( mesh, problem, geom, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: geom, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    select case ( mesh%ndim )
      case(2) ! 2D and axisymmetric
        call poisson_natboun_curve ( mesh, problem, geom, elem, matrix, &
          vector, first, last, coefficients, oldvectors, elemmat, elemvec )
      case(3) ! 3D
        call poisson_natboun_surface ( mesh, problem, geom, elem, matrix, &
          vector, first, last, coefficients, oldvectors, elemmat, elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in poisson_natboun:', &
        ' poisson element not available for mesh%ndim = ', mesh%ndim
        stop
    end select

  end subroutine poisson_natboun


! Boundary element for a natural boundary on a curve for the Poisson equation

  subroutine poisson_natboun_curve ( mesh, problem, curve, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use poisson_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip, orientation, funcnr, vfuncnr


    if ( first ) then

!     first element on this curve

      call check ( coefficients, 'poisson_natboun_curve', indexarray=[6], &
        minimum=[0], maximum=[2] )

!     set globals

      call set_globals_poisson_boun ( mesh, coefficients, curve )

      allocate ( wg(ninti), fg(ninti), curvel(ninti), normal(ninti,ndim) )
      allocate ( xig(ninti,1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( sigma(ninti) )

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

    funcnr = coefficients%i(15)
    vfuncnr = coefficients%i(16)

    if ( funcnr > 0 .or. coefficients%i(6) == 1 ) then

      do ip = 1, ninti
        fg(ip) = coefficients%func ( funcnr, xg(ip,:) )
      end do

    else if ( vfuncnr > 0 .or. coefficients%i(6) == 2 ) then

      call check ( coefficients, 'poisson_natboun_curve', indexarray=[17], &
        minimum=[-1], maximum=[1] )

      orientation = get_coefficient ( coefficients, index=17, default=1 )

      normal = real ( orientation, kind=dp ) * normal

      do ip = 1, ninti
        fg(ip) = &
          dot_product ( normal(ip,:), &
                             coefficients%vfunc ( ndim, vfuncnr, xg(ip,:) ) )
      end do

    else if ( coefficients%i(6) == 0 ) then

      fg = coefficients%r(23)

    end if

    do i = 1, ndf
      elemvec(i) = - sum ( fg * phi(:,i) * curvel * wg )
    end do

    if ( matrix ) then

!     Robin boundary condition

      if ( coefficients%i(5) == 1 ) then

        do ip = 1, ninti
          sigma(ip) = coefficients%func1(7)%p ( coefficients%i(7), xg(ip,:) )
        end do

      else

        sigma = coefficients%r(2)

      end if

      do i = 1, ndf
        do j = i, ndf
          elemmat(i,j) = sum ( sigma * phi(:,i) * phi(:,j) * curvel * wg )
          elemmat(j,i) = elemmat(i,j)  ! symmetry
        end do
      end do

    end if

    if ( last ) then

!     last element on this curve

      deallocate ( wg, fg, curvel, normal )
      deallocate ( xig, phi, x )
      deallocate ( sigma )
      deallocate ( xg )
      deallocate ( dphi, dxdxi )

    end if

  end subroutine poisson_natboun_curve


! Boundary element for a natural boundary on a surface for the Poisson equation

  subroutine poisson_natboun_surface ( mesh, problem, surface, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use poisson_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: surface, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, ip, orientation, funcnr, vfuncnr


    if ( first ) then

!     first element on this surface

      call check ( coefficients, 'poisson_natboun_surface', indexarray=[6], &
        minimum=[0], maximum=[2] )

!     set globals

      call set_globals_poisson_boun ( mesh, coefficients, surface=surface )

      allocate ( wg(ninti), fg(ninti), surfl(ninti), normal(ninti,ndim) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( sigma(ninti) )

      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl, normal )

    call isoparametric_coordinates ( x, phi, xg )

    funcnr = coefficients%i(15)
    vfuncnr = coefficients%i(16)

    if ( funcnr > 0 .or. coefficients%i(6) == 1 ) then

      do ip = 1, ninti
        fg(ip) = coefficients%func ( funcnr, xg(ip,:) )
      end do

    else if ( vfuncnr > 0 .or. coefficients%i(6) == 2 ) then

      call check ( coefficients, 'poisson_natboun_surface', indexarray=[17], &
        minimum=[-1], maximum=[1] )

      orientation = get_coefficient ( coefficients, index=17, default=1 )

      normal = real ( orientation, kind=dp ) * normal

      do ip = 1, ninti
        fg(ip) = &
          dot_product ( normal(ip,:), &
                           coefficients%vfunc ( ndim, vfuncnr, xg(ip,:) ) )
      end do

    else if ( coefficients%i(6) == 0 ) then

      fg = coefficients%r(23)

    end if

    do i = 1, ndf
      elemvec(i) = - sum ( fg * phi(:,i) * surfl * wg )
    end do

    if ( matrix ) then

!     Robin boundary condition

      if ( coefficients%i(5) == 1 ) then

        do ip = 1, ninti
          sigma(ip) = coefficients%func1(7)%p ( coefficients%i(7), xg(ip,:) )
        end do

      else

        sigma = coefficients%r(2)

      end if

      do i = 1, ndf
        do j = i, ndf
          elemmat(i,j) = sum ( sigma * phi(:,i) * phi(:,j) * surfl * wg )
          elemmat(j,i) = elemmat(i,j)  ! symmetry
        end do
      end do

    end if

    if ( last ) then

!     last element on this surface

      deallocate ( wg, fg, surfl, normal )
      deallocate ( xig, phi, x )
      deallocate ( sigma )
      deallocate ( xg )
      deallocate ( dphi, dxdxis )

    end if

  end subroutine poisson_natboun_surface


! compute derivatives

  subroutine poisson_deriv ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use poisson_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: j


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_poisson ( mesh, coefficients, elgrp )

      allocate ( detF(nodalp) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf), x(nodalp,ndim) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( tmp(nodalp,ndim), u(ndf) )

!     set reference coordinates and shape function

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefunc, xrnod, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      layer=layer )

    do j = 1, ndim
      tmp(:,j) = matmul ( dphidx(:,:,j), u )
    end do

    elemvec = reshape ( tmp, [ ndim*nodalp ] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( detF )
      deallocate ( xrnod, phi, x )
      deallocate ( tmp, u )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )

    end if

  end subroutine poisson_deriv


! Element for the constraints (connection through collocation)

  subroutine poisson_node_conn ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

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

  end subroutine poisson_node_conn


! Element for the constraints (connection through elements)

  subroutine poisson_constr_elem_conn ( mesh, problem, constr, elem, &
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

    select case ( mesh%ndim )
      case(2) ! 2D and axisymmetric
        call poisson_constr_elem_conn_curve ( mesh, problem, constr, elem, &
          node, matrix, vector, first, last, coefficients, oldvectors, &
          elemmat, elemmat2, elemmatadd, elemvec, elemvecadd )
      case(3) ! 3D
        call poisson_constr_elem_conn_surface ( mesh, problem, constr, elem, &
          node, matrix, vector, first, last, coefficients, oldvectors, &
          elemmat, elemmat2, elemmatadd, elemvec, elemvecadd )
      case default
        write(*,'(/a/a,i0/)') 'Error in poisson_constr_elem_conn:', &
        ' poisson element not available for mesh%ndim = ', mesh%ndim
        stop
    end select

  end subroutine poisson_constr_elem_conn


! Element for weak connection through elements for a curve

  subroutine poisson_constr_elem_conn_curve ( mesh, problem, constr, elem, &
    node, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    use poisson_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


    integer :: curve, i, j


!   connection through elements (weak)

    curve = problem%constraints(constr)%geometry1

    if ( first ) then

!     first element on this surface

      call set_globals_poisson_boun ( mesh, coefficients, curve=curve )
      call set_globals_poisson_l_boun ( coefficients )

      allocate ( wg(ninti), curvel(ninti) )
      allocate ( xig(ninti,1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( psi(ninti,ndfl) )
      allocate ( dphi(ninti,ndf,1), dxdxi(ninti,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      call set_shape_function ( shapefuncl, xig, psi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, curve=curve  )

    call isoparametric_deformation_curve ( x, dphi(:,:,1), dxdxi, curvel )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x, phi, xg )
      curvel = 2 * pi * xg(:,2) * curvel
    end if

    do i = 1, ndfl
      do j = 1, ndf
        elemmat(i,j) = sum ( psi(:,i) * phi(:,j) * curvel * wg )
      end do
    end do

    elemmat2 = - elemmat
    elemvec  = 0

    if ( last ) then

!     last element on this surface

      deallocate ( wg, curvel )
      deallocate ( xig, phi, x )
      deallocate ( psi )
      deallocate ( dphi, dxdxi )

    end if

  end subroutine poisson_constr_elem_conn_curve


! Element for weak connection through elements for a surface

  subroutine poisson_constr_elem_conn_surface ( mesh, problem, constr, elem, &
    node, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    use poisson_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


    integer :: surface, i, j


!   connection through elements (weak)

    surface = problem%constraints(constr)%geometry1

    if ( first ) then

!     first element on this surface

      call set_globals_poisson_boun ( mesh, coefficients, surface=surface )
      call set_globals_poisson_l_boun ( coefficients )

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( psi(ninti,ndfl) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      call set_shape_function ( shapefuncl, xig, psi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, surface=surface  )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl )

    do i = 1, ndfl
      do j = 1, ndf
        elemmat(i,j) = sum ( psi(:,i) * phi(:,j) * surfl * wg )
      end do
    end do

    elemmat2 = - elemmat
    elemvec  = 0

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl )
      deallocate ( xig, phi, x )
      deallocate ( psi )
      deallocate ( dphi, dxdxis )

    end if

  end subroutine poisson_constr_elem_conn_surface


! sample value of scalar in one node of the object

  subroutine poisson_sample ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, uvec )

    use poisson_globals_m

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

    call set_globals_poisson ( mesh, coefficients, elgrp )

    allocate ( phi(1,ndf), u(ndf) )

!   set shape function in the point

    call set_shape_function ( shapefunc, xr, phi )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      layer=layer )

    uvec = matmul( phi, u )

    deallocate ( phi, u )

  end subroutine poisson_sample


! sample gradient of scalar in one node of the object

  subroutine poisson_sample_deriv ( mesh, problem, object, nodeobj, &
    coefficients, oldvectors, uvec )

    use poisson_globals_m

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

    call set_globals_poisson ( mesh, coefficients, elgrp )

    allocate ( detF(1) )
    allocate ( phi(1,ndf), x(nodalp,ndim) )
    allocate ( u(ndf) )
    allocate ( dphi(1,ndf,ndim), F(1,ndim,ndim) )
    allocate ( Finv(1,ndim,ndim), dphidx(1,ndf,ndim) )

!   set shape function in the point

    call set_shape_function ( shapefunc, xr, phi, dphi )

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      layer=layer )

    uvec = matmul ( u, dphidx(1,:,:) )

    deallocate ( detF )
    deallocate ( phi, x )
    deallocate ( u )
    deallocate ( dphi, F )
    deallocate ( Finv, dphidx )

  end subroutine poisson_sample_deriv


! set global parameters (internal element)

  subroutine set_globals_poisson ( mesh, coefficients, elgrp )

    use poisson_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

    integer :: ndimr

!   check size of coefficients

    call check ( coefficients, 'set_globals_poisson', ncoefi=100, &
      ncoefr=50, indexarray=[23], minimum=[0], maximum=[1] )

    layer = coefficients%i(38)

    ndim = mesh%element(elgrp)%ndim
    nodalp = mesh%element(elgrp)%numnod
    if ( ndim == 2 ) then
      coorsys = coefficients%i(23)
    else
      coorsys = 2
    end if
    nsides = mesh%element(elgrp)%numsides
    nodalpb = mesh%element(elgrp)%sidnumnod
    globalshape = mesh%element(elgrp)%globalshape

!   set number of degrees of freedom

    intpol = coefficients%i(1)

    ndimr = mesh%element(elgrp)%ndimr
    if ( intpol == 13 .and. any(mesh%element(elgrp)%p(:ndimr,2) /= 1 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_poisson:', &
        ' For spectral elements GLL nodal distribution is required. '
      stop
    else if ( intpol == 20 .and. &
                            any(mesh%element(elgrp)%p(:ndimr,2) /= 0 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_poisson:', &
        ' For high-order elements equidistant nodal distribution is required. '
      stop
    end if

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(39)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, name_of_routine='set_globals_poisson', ndf=ndf )

    if ( ndf /= nodalp ) then
      write(*,'(/4(a/))') 'Error in set_globals_poisson:', &
        ' Number of degrees of freedom of the shape function (ndf) ', &
        ' is different from the number of nodal points in the ', &
        ' element (nodalp). This possibility (ndf /= nodalp) is not available.'
      stop
    end if

!   set integration

    inttype = coefficients%i(40)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_poisson:', &
        ' For spectral elements Gauss-Legendre-Lobatto integration ', &
        ' needs to be specified. '
      stop
    end if

    if ( coefficients%i(50) == 1 .or. SET_GAUSS_BY_ORDER ) then
      intrule = set_intrule ( globalshape, inttype, order=coefficients%i(10) )
      intrule2 = set_intrule2 ( globalshape, inttype, order=coefficients%i(20) )
    else
      intrule = coefficients%i(10)
      intrule2 = coefficients%i(20)
    end if
    nsubint = get_coefficient ( coefficients, index=32, default=1 )

    if ( globalshape == 'prism' .and. intrule2 == 0 ) then
      write(*,'(/a/3a/)') 'Error in set_globals_poisson:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape
      stop
    end if

    if ( globalshape == 'pyramid' .and. inttype == 2 .and. intrule2 == 0 ) then
      write(*,'(/a/a/3a,i0/)') 'Error in set_globals_poisson:', &
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
      write(*,'(/5(a/))') 'Error in set_globals_poisson:', &
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_poisson


! set global parameters (boundary element)

  subroutine set_globals_poisson_boun ( mesh, coefficients, curve, surface, &
    volume, ndimr, geometry )

    use poisson_globals_m
    use set_optional_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in), optional :: curve, surface, volume
    integer, intent(in), optional :: ndimr, geometry

    integer :: lcurve, lsurface, lvolume


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
          call errormsg_case_default ( 'set_globals_poisson_boun', &
            'ndimr', int_value=ndimr )
      end select
    end if

    if ( lcurve > 0 ) then
      ndim = mesh%curves(lcurve)%ndim
      if ( ndim == 3 ) then
        coorsys = 2
      else
        coorsys = coefficients%i(23)
      end if
      globalshape = mesh%curves(lcurve)%element%globalshape
      nodalp = mesh%curves(lcurve)%element%numnod
    else if ( lsurface > 0 ) then
      ndim = mesh%surfaces(lsurface)%ndim
      if ( ndim == 3 ) then
        coorsys = 2
      else
        coorsys = coefficients%i(23)
      end if
      globalshape = mesh%surfaces(lsurface)%element%globalshape
      nodalp = mesh%surfaces(lsurface)%element%numnod
    else if ( lvolume > 0 ) then
      ndim = mesh%volumes(lvolume)%ndim
      coorsys = 2 ! always 3D
      globalshape = mesh%volumes(lvolume)%element%globalshape
      nodalp = mesh%volumes(lvolume)%element%numnod
    end if

!   set number of degrees of freedom

    intpol = coefficients%i(1)

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(39)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, name_of_routine='set_globals_poisson_boun', &
      ndf=ndf )

    if ( ndf /= nodalp ) then
      write(*,'(/5(a/))') 'Error in set_globals_poisson_boun:', &
        ' Number of degrees of freedom of the shape function (ndf) ', &
        ' is different from the number of nodal points in the element', &
        ' (nodalp) ', &
        ' This possibility (ndf /= nodalp) is not available.'
      stop
    end if

!   set integration

    inttype = coefficients%i(40)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_poisson_boun:', &
        ' For spectral elements Gauss-Legendre-Lobatto integration ', &
        ' needs to be specified. '
      stop
    end if

    if ( coefficients%i(50) == 1 .or. SET_GAUSS_BY_ORDER ) then
      intrule = set_intrule ( globalshape, inttype, order=coefficients%i(11) )
      intrule2 = set_intrule2 ( globalshape, inttype, order=coefficients%i(21) )
    else
      intrule = coefficients%i(11)
      intrule2 = coefficients%i(21)
    end if
    nsubint = get_coefficient ( coefficients, index=33, default=1 )

    if ( globalshape == 'prism' .and. intrule2 == 0 ) then
      write(*,'(/a/3a/)') 'Error in set_globals_poisson_boun:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape
      stop
    end if

    if ( globalshape == 'pyramid' .and. inttype == 2 .and. intrule2 == 0 ) then
      write(*,'(/a/a/3a,i0/)') 'Error in set_globals_poisson_boun:', &
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
      write(*,'(/5(a/))') 'Error in set_globals_poisson_boun:', &
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_poisson_boun


! set global parameters (boundary Langrangian multiplier)

  subroutine set_globals_poisson_l_boun ( coefficients )

    use poisson_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   Lagrange multiplier interpolation

!   set number of degrees of freedom

    call set_ndf ( shapefunc, name_of_routine='set_globals_poisson_l_boun', &
      ndfl=ndfl, intpoll=intpoll )

    shapefuncl = shapefunc
    shapefuncl%interpolation = intpoll

  end subroutine set_globals_poisson_l_boun

end module poisson_elements_m

