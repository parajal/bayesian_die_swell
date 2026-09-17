
! Copyright (C) 2005-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the mixed Stokes equation:
!
!               nabla p - div( tau )  = ...   (1)
!       div u                         = 0     (2)
!    nabla u+nabla u^T  - tau / 2*eta = 0     (3)
!
! Note that here only the "div tau" term in Eq. (1) and Eq.(3) are build here.
! Others terms/equations must be build with, for example, the Stokes elements.

module mixed_stokes_elements_2D_m

  use tfem_elem_m
  use stokes_set_globals_m

  implicit none

contains


! Internal element routine for the mixed-Stokes Equation

  subroutine mixed_stokes_elem ( mesh, problem, elgrp, elem, matrix, vector, &
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


    integer :: i, j, i1, i2, i3, i4, i5
    real(dp) :: eta_mix


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_mixed_stokes ( mesh, coefficients, elgrp )

!     allocate arrays

      allocate ( wg(ninti), detF(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( theta(ninti,ndft) )
      allocate ( Bu11(ndf,ndft), Bu12(ndf,ndft), Amat(ndft,ndft) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )
      call set_shape_function ( shapefunct, xig, theta )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    call shape_derivative ( dphi, Finv, dphidx )

!   build equations

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

      do i = 1, ndf
        do j = 1, ndft
          Bu11(i,j) = sum ( dphidx(:,i,1) * theta(:,j) * detF * wg )
          Bu12(i,j) = sum ( dphidx(:,i,2) * theta(:,j) * detF * wg )
        end do
      end do

      do i = 1, ndft
        do j = 1, ndft
          Amat(i,j) = sum ( theta(:,i) * theta(:,j) * detF * wg )
        end do
      end do

      eta_mix = coefficients%r(5)

      Amat = - Amat / 2 / eta_mix

!     pointers in unknowns

      i1 = ndft     ! tau11
      i2 = 2*ndft   ! tau12
      i3 = 3*ndft   ! tau22
      i4 = i3 + ndf    ! u
      i5 = i3 + 2*ndf  ! v

!     (s,tau)

      elemmat(    1:i1,    1:i1 ) = Amat
      elemmat(    1:i1, i1+1:i3 ) = 0
      elemmat( i1+1:i2,    1:i1 ) = 0
      elemmat( i1+1:i2, i1+1:i2 ) = 2*Amat
      elemmat( i1+1:i2, i2+1:i3 ) = 0
      elemmat( i2+1:i3,    1:i2 ) = 0
      elemmat( i2+1:i3, i2+1:i3 ) = Amat

!     ( (nabla v)^T, tau )

      elemmat( i3+1:i4,    1:i1 ) = Bu11
      elemmat( i3+1:i4, i1+1:i2 ) = Bu12
      elemmat( i3+1:i4, i2+1:i3 ) = 0
      elemmat( i4+1:i5,    1:i1 ) = 0
      elemmat( i4+1:i5, i1+1:i2 ) = Bu11
      elemmat( i4+1:i5, i2+1:i3 ) = Bu12

!     transposed of ( (nabla v)^T, tau )

      elemmat(    1:i3, i3+1:i5 ) = transpose( elemmat( i3+1:i5,   1:i3 ) )

!     velocity - velocity part

      elemmat( i3+1:i5, i3+1:i5 ) = 0

    end if

    if ( last ) then

!     last element in this group

      deallocate ( Finv, dphidx )
      deallocate ( dphi, F )
      deallocate ( xg )
      deallocate ( xig, phi, x )
      deallocate ( wg, detF )
      deallocate ( theta )
      deallocate ( Bu11, Bu12, Amat )

    end if

  end subroutine mixed_stokes_elem


! compute stresses in all nodes (mixed stresses)

  subroutine mixed_stokes_stress ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    integer :: comp


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_mixed_stokes ( mesh, coefficients, elgrp )

!     check coefficients%i(13)

      call check ( coefficients, 'mixed_stokes_stress', indexarray=[13], &
        minimum=[1], maximum=[3] )

      physqstress = coefficients%i(8)

!     allocate arrays

      allocate ( theta(nodalp,ndft) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( u(3*ndft), tmp(ndft,3) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     set shape function in nodal points

      call set_shape_function ( shapefunct, xrnod, theta )

    end if

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqstress], layer=layer )

    tmp = reshape ( u, [ndft,3] )

    comp = coefficients%i(13)

    elemvec = matmul ( theta, tmp(:,comp) )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( u, tmp )

    end if

  end subroutine mixed_stokes_stress


! set global parameters mixed-stokes (internal element)

  subroutine set_globals_mixed_stokes ( mesh, coefficients, elgrp )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

!   check size of coefficients

    call check ( coefficients, 'set_globals_mixed_stokes', &
      indexarray=[37], minimum=[0], maximum=[0] )

!   set number of degrees of freedom for the stress

    intpolt = coefficients%i(3)

    shapefunct%globalshape = globalshape
    shapefunct%interpolation = intpolt
    shapefunct%numbering = 'standard'
    shapefunct%p = coefficients%i(82)
    shapefunct%spec_eval = 'gauss'

    call set_ndf ( shapefunct, 'set_globals_mixed_stokes', ndf=ndft )

  end subroutine set_globals_mixed_stokes

end module mixed_stokes_elements_2D_m

