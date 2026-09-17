
! Copyright (C) 2007-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the DEVSS and DEVSS-G formulations
!
! DEVSS:
!
!   ... - 2 * alpha * div ( D(u) - E ) = ....
!
!                           D(u) - E   = 0
!
! Note: the last equation is multiplied by 2*alpha in the implementation to
! get a symmetric system matrix.
!
! DEVSS-G:
!
!   ... - alpha * div ( nabla u - G^T ) = ....
!
!                       nabla u - G^T   = 0
!
! Note: the last equation is multiplied by alpha in the implementation to
! get a symmetric system matrix.
!
! This module contains the generic stuff

module devss_elements_generic_m

  use tfem_elem_m
  use stokes_set_globals_m

  implicit none

contains


! Element routine for the right-hand side of the momentum balance:
!
!   ... = alpha * div ( nabla u - G^T )
!
! with u and G known, for example in a Newton-Raphson iteration scheme.

  subroutine devssg_rhs_div ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    real(dp) :: alpha
    integer :: j, ip


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors, maxvel3D=1 )

    call set_devssg_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

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

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )


!   get velocity and compute gradients

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do


!   get gradient

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, gtmp, &
      physq=[physqgrad], layer=layer )

    tmp = reshape ( gtmp, [ndfg,ncompg] )

    work4 = matmul ( zeta, tmp )

    do ip = 1, ninti
      ggrad(ip,:,:) = reshape ( work4(ip,:), [ncompu,ndim], order=[2,1] )
    end do


!   build equations

    if ( matrix ) then

      write(*,'(/a/a/)') 'Error in devssg_rhs_div: no matrix to build.', &
        'Call build_system with buildmatrix=.false.'
      stop

    end if

    if ( vector ) then

!     - alpha * (nabla v)^T: ( nabla u - G^T )

      alpha = coefficients%r(4)

      do ip = 1, ninti
        work6(ip,:,:) = alpha * matmul ( dphidx(ip,:,:), &
                                 transpose ( gradu(ip,:,:) - ggrad(ip,:,:) ) )
      end do

      do j = 1, ncompu
        work2(:,j) = - matmul ( detF * wg, work6(:,:,j) )
      end do

      elemvec = reshape ( work2, [ ndf*ncompu ] )

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_devssg_elem ( last, coefficients )

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
      allocate ( work4(ninti,ncompg) )
      allocate ( u(ncompu*ndf), uvector(ndf,ncompu) )
      allocate ( gradu(ninti,ncompu,ndim) )
      allocate ( gtmp(ndfg*ncompg), tmp(ndfg,ncompg) )
      allocate ( ggrad(ninti,ncompu,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work6, work2 )
      deallocate ( work4 )
      deallocate ( u, uvector, gradu )
      deallocate ( gtmp, tmp )
      deallocate ( ggrad )

    end subroutine deallocate_arrays

  end subroutine devssg_rhs_div


! Element routine for computing G from the velocity: (H, G - grad v^T) = 0

  subroutine gradient_from_velocity_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

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

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_devssg_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )


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

!   get velocity and compute gradients

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ncompu] )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

!   build equations

    if ( vector ) then

      if ( ndim == 2 .and. vel3D == 0 ) then

        do i = 1, ndfg
          work2(i,1) = sum ( zeta(:,i) * gradu(:,1,1) * detF * wg )
          work2(i,2) = sum ( zeta(:,i) * gradu(:,1,2) * detF * wg )
          work2(i,3) = sum ( zeta(:,i) * gradu(:,2,1) * detF * wg )
          work2(i,4) = sum ( zeta(:,i) * gradu(:,2,2) * detF * wg )
        end do

      else if ( ndim == 2 .and. vel3D == 1 ) then

        do i = 1, ndfg
          work2(i,1) = sum ( zeta(:,i) * gradu(:,1,1) * detF * wg )
          work2(i,2) = sum ( zeta(:,i) * gradu(:,1,2) * detF * wg )
          work2(i,3) = sum ( zeta(:,i) * gradu(:,2,1) * detF * wg )
          work2(i,4) = sum ( zeta(:,i) * gradu(:,2,2) * detF * wg )
          work2(i,5) = sum ( zeta(:,i) * gradu(:,3,1) * detF * wg )
          work2(i,6) = sum ( zeta(:,i) * gradu(:,3,2) * detF * wg )
        end do

      else

        do i = 1, ndfg
          work2(i,1) = sum ( zeta(:,i) * gradu(:,1,1) * detF * wg )
          work2(i,2) = sum ( zeta(:,i) * gradu(:,1,2) * detF * wg )
          work2(i,3) = sum ( zeta(:,i) * gradu(:,1,3) * detF * wg )
          work2(i,4) = sum ( zeta(:,i) * gradu(:,2,1) * detF * wg )
          work2(i,5) = sum ( zeta(:,i) * gradu(:,2,2) * detF * wg )
          work2(i,6) = sum ( zeta(:,i) * gradu(:,2,3) * detF * wg )
          work2(i,7) = sum ( zeta(:,i) * gradu(:,3,1) * detF * wg )
          work2(i,8) = sum ( zeta(:,i) * gradu(:,3,2) * detF * wg )
          work2(i,9) = sum ( zeta(:,i) * gradu(:,3,3) * detF * wg )
        end do

      end if

      elemvec = reshape ( work2, [ ndfg*ncompg ] )

    end if

    if ( matrix ) then

      elemmat = 0

!     mass matrix

      do i = 1, ndfg
        do j = 1, ndfg
          elemmat(i,j) = sum ( zeta(:,i) * zeta(:,j) * detF * wg )
        end do
      end do

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_devssg_elem ( last, coefficients )

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

      allocate ( u(ncompu*ndf), uvector(ndf,ncompu) )
      allocate ( gradu(ninti,ncompu,ndim), work2(ndfg, ncompg) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( u, uvector )
      deallocate ( gradu, work2 )

    end subroutine deallocate_arrays

  end subroutine gradient_from_velocity_elem


! derive element for getting the projected velocity gradient in all nodes.

  subroutine deriv_gradient_tensor ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    integer :: i

    if ( first ) then

!     first element in this group

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

      call set_globals_devssg ( mesh, coefficients, elgrp )

      physqgrad = coefficients%i(9)

      allocate ( zeta(nodalp,ndfg) )
      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf) )
      allocate ( u(ndfg*ncompg), tmp(ndfg,ncompg) )
      allocate ( gradu(nodalp,ncompu,ncompu), work2(nodalp,ncompg) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefuncg, xrnod, zeta )

    end if

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqgrad], layer=layer )

    tmp = reshape ( u, [ndfg,ncompg] )

    if ( vel3D == 1 ) then

      work2 = matmul ( zeta, tmp )

      gradu(:,1,1) = work2(:,1)
      gradu(:,1,2) = work2(:,2)
      gradu(:,2,1) = work2(:,3)
      gradu(:,2,2) = work2(:,4)
      gradu(:,3,1) = work2(:,5)
      gradu(:,3,2) = work2(:,6)
      gradu(:,:,3) = 0

      do i = 1, nodalp
         gradu(i,:,:) = transpose(gradu(i,:,:)) ! make row -> column
      end do

      elemvec = reshape ( gradu, [nodalp*ncompu**2] )

    else

      elemvec = reshape ( matmul ( zeta, tmp ), [nodalp*ncompg] )

    end if

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( zeta )
      deallocate ( xrnod, phi )
      deallocate ( u, tmp )
      deallocate ( gradu, work2 )

    end if

  end subroutine deriv_gradient_tensor


! set global parameters DEVSS (internal element)

  subroutine set_globals_devss ( mesh, coefficients, elgrp )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

    if ( mesh%ndim == 2 ) then
      ncompe = 3
    else if ( mesh%ndim == 3 ) then
      ncompe = 6
    end if

    physqgrad = coefficients%i(9)

!   set number of degrees of freedom for E

    intpole = coefficients%i(4)

    shapefunce%globalshape = globalshape
    shapefunce%interpolation = intpole
    shapefunce%numbering = 'regular'

    call set_ndf ( shapefunce, 'set_globals_devss', ndf=ndfe )

  end subroutine set_globals_devss


! set global parameters DEVSS-G (internal element)

  subroutine set_globals_devssg ( mesh, coefficients, elgrp )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

    if ( mesh%ndim == 2 .and. vel3D == 0 ) then
      ncompg = 4
    else if ( mesh%ndim == 2 .and. vel3D == 1 ) then
      ncompg = 6
    else if ( mesh%ndim == 3 ) then
      ncompg = 9
    end if

    physqgrad = coefficients%i(9)

!   set number of degrees of freedom for G

    intpolg = coefficients%i(5)

    shapefuncg%globalshape = globalshape
    shapefuncg%interpolation = intpolg
    shapefuncg%numbering = 'regular'

    call set_ndf ( shapefuncg, 'set_globals_devssg', ndf=ndfg )

  end subroutine set_globals_devssg


! Preamble for the devss element

  subroutine set_devss_elem ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors


    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

!     set globals

      call set_globals_devss ( mesh, coefficients, elgrp )

!     allocate arrays

      allocate ( zeta(ninti,ndfe) )

!     set shape function
      call set_shape_function ( shapefunce, xig, zeta )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      if ( first ) then

!       first element in this group

!       set globals

        call set_globals_devss ( mesh, coefficients, elgrp )

      end if

!     allocate arrays

      allocate ( zeta(ninti,ndfe) )

!     set shape function
      call set_shape_function ( shapefunce, xig, zeta )

    end if

  end subroutine set_devss_elem


! Unset the preamble for the devss element
! (deallocate arrays allocated in set_... )

  subroutine unset_devss_elem ( last, coefficients )

    use stokes_globals_m

    logical, intent(in) :: last
    type(coefficients_t), intent(in) :: coefficients

!   deallocate memory

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      deallocate ( zeta )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      deallocate ( zeta )

    end if

  end subroutine unset_devss_elem


! Preamble for the devssg element

  subroutine set_devssg_elem ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors


    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

!     set globals

      call set_globals_devssg ( mesh, coefficients, elgrp )

!     allocate arrays

      allocate ( zeta(ninti,ndfg) )

!     set shape function
      call set_shape_function ( shapefuncg, xig, zeta )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      if ( first ) then

!       first element in this group

!       set globals

        call set_globals_devssg ( mesh, coefficients, elgrp )

      end if

!     allocate arrays

      allocate ( zeta(ninti,ndfg) )

!     set shape function
      call set_shape_function ( shapefuncg, xig, zeta )

    end if

  end subroutine set_devssg_elem


! Unset the preamble for the devss element
! (deallocate arrays allocated in set_... )

  subroutine unset_devssg_elem ( last, coefficients )

    use stokes_globals_m

    logical, intent(in) :: last
    type(coefficients_t), intent(in) :: coefficients

!   deallocate memory

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      deallocate ( zeta )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      deallocate ( zeta )

    end if

  end subroutine unset_devssg_elem

end module devss_elements_generic_m

