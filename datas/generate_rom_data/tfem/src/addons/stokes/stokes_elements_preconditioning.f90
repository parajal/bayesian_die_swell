
! Copyright (C) 2009-2017 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the pressure preconditioning (using an approximate
! Schur complement) as decribed in the book:
! H. Elman, D. Silvester, A. Wathen: "Finite elements and fast iterative
! solvers", Oxford University Press, 2005.

module stokes_elements_preconditioning_m

  use tfem_elem_m
  use stokes_set_globals_m

  implicit none


contains


! Internal element routine for the pressure mass matrix Mp

  subroutine pressure_mass_elem ( mesh, problem, elgrp, elem, matrix, &
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

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( matrix ) then

!     mass matrix

      do i = 1, ndfp
        do j = i, ndfp
          elemmat(i,j) = sum ( psi(:,i) * psi(:,j) * detF * wg )
          elemmat(j,i) = elemmat(i,j)
        end do
      end do

    end if

    if ( vector ) then

      elemvec = 0

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )


  end subroutine pressure_mass_elem


! Internal element routine for the pressure diffusion matrix A_p

  subroutine pressure_diffusion_elem ( mesh, problem, elgrp, elem, matrix, &
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


    integer :: i, j, ip


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

!   derivative of shape function for pressure

    call set_shape_function ( shapefuncp, xig, psi, dpsi )


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dpsi, Finv, dpsidx )

    if ( matrix ) then

!     diffusion matrix

      do i = 1, ndfp
        do j = i, ndfp
          do ip = 1, ninti
            work(ip) = sum ( dpsidx(ip,i,:) * dpsidx(ip,j,:) )
          end do
          elemmat(i,j) = sum ( work * detF * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

    end if

    if ( vector ) then

      elemvec = 0

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

      allocate ( work(ninti), dpsi(ninti,ndfp,ndim), dpsidx(ninti,ndfp,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work, dpsi, dpsidx )

    end subroutine deallocate_arrays


  end subroutine pressure_diffusion_elem


! Internal element routine for the pressure convection matrix N_p

  subroutine pressure_convection_elem ( mesh, problem, elgrp, elem, matrix, &
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


    integer :: i, j, ip


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

!   derivative of shape function for pressure

    call set_shape_function ( shapefuncp, xig, psi, dpsi )


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dpsi, Finv, dpsidx )


!   get velocity vector

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, w, physq=[physqvel], layer=layer )

    tmp = reshape ( w, [ndf,ndim] )

    wvec = matmul ( phi, tmp )

!   w.grad operator

    do ip = 1, ninti
      wgradpsi(ip,:) = matmul ( dpsidx(ip,:,:), wvec(ip,:) )
    end do


    if ( matrix ) then

!     convection matrix

      do i = 1, ndfp
        do j = 1, ndfp
          elemmat(i,j) = sum ( psi(:,i) * wgradpsi(:,j) * detF * wg )
        end do
      end do

    end if

    if ( vector ) then

      elemvec = 0

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

      allocate ( work(ninti), dpsi(ninti,ndfp,ndim), dpsidx(ninti,ndfp,ndim) )
      allocate ( w(ndf*ndim), tmp(ndf,ndim), wvec(ninti,ndim) )
      allocate ( wgradpsi(ninti,ndfp) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work, dpsi, dpsidx )
      deallocate ( w, tmp, wvec )
      deallocate ( wgradpsi )

    end subroutine deallocate_arrays


  end subroutine pressure_convection_elem


end module stokes_elements_preconditioning_m
