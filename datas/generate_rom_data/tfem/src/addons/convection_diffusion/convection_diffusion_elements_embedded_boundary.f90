
! Copyright (C) 2011-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for embedded boundary conditions (Dirichlet, ... ) for the
! convection-diffusion equation using an eltree to divide the domain
! and/or introduce a boundary inside an element.

module convection_diffusion_elements_embedded_boundary_m

  use tfem_elem_m
  use convection_diffusion_set_globals_m

  implicit none


contains


! open (embedded) boundary element using an eltree for
! the scalar diffusion equation with varying coefficient alpha
!   - nabla ( alpha nabla c ) = f
! Coefficient alpha is a scalar.

  subroutine scalar_diffusion_open_boundary_eltree ( mesh, problem, &
    elgrp, elem, matrix, vector, first, last, coefficients, oldvectors, &
    elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   variables

    integer :: N, M, ip, orientation


!   set globals element

    call set_globals_poisson ( mesh, coefficients, elgrp )
    call set_globals_convection_diffusion ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_convection_diffusion_eltree ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )

!   allocate arrays

    allocate ( work(ninti) )
    allocate ( alphag(ninti) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
    end do

!   set orientation

    call check ( coefficients, 'scalar_diffusion_open_boundary_eltree', &
      indexarray=[17], minimum=[-1], maximum=[1] )

    orientation = get_coefficient ( coefficients, index=17, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   vector

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     boundary flux term

      call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
        elgrp, elem, choice=coefficients%i(2), value=coefficients%r(1), &
        func=coefficients%func1(1)%p, funcnr=coefficients%i(3), x=xg, indx_v=1,&
        layer=layer, phi=phi, indx_e=1, coef=alphag )

      do N = 1, ndf
        do M = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,M,:) * dan(ip,:) )
          end do
!         alpha_g * phi(:,N) * dphidx(:,M,i) * dan(:,i)
          elemmat(N,M) = - sum ( alphag * phi(:,N)* work )
        end do
      end do

    end if

!   deallocate arrays

    deallocate ( work )
    deallocate ( alphag )

    call unset_convection_diffusion_eltree

  end subroutine scalar_diffusion_open_boundary_eltree


! open (embedded) boundary element using an eltree for
! the scalar diffusion equation with varying coefficient alpha
!   - nabla ( alpha nabla c ) = f
! Coefficient alpha is a tensor.

  subroutine scalar_diffusion2_open_boundary_eltree ( mesh, problem, &
    elgrp, elem, matrix, vector, first, last, coefficients, oldvectors, &
    elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   variables

    integer :: N, M, ip, orientation


!   set globals element

    call set_globals_poisson ( mesh, coefficients, elgrp )
    call set_globals_convection_diffusion ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_convection_diffusion_eltree ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )

!   allocate arrays

    allocate ( work(ninti) )
    allocate ( alphatg(ninti,ndim,ndim) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
    end do

!   set orientation

    call check ( coefficients, 'scalar_diffusion2_open_boundary_eltree', &
      indexarray=[17], minimum=[-1], maximum=[1] )

    orientation = get_coefficient ( coefficients, index=17, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   vector

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     boundary flux term

      call evaluate_tensor_coefficient ( mesh, problem, oldvectors, &
        elgrp, elem, choice=coefficients%i(2), &
        value=reshape(coefficients%r(9:8+ndim**2), [ndim,ndim]), &
        tfunc=coefficients%tfunc1(1)%p, tfuncnr=coefficients%i(3), x=xg, &
        indx_v=1, layer=layer, phi=phi, indx_e=1, coef=alphatg )

      do N = 1, ndf
        do M = 1, ndf
          do ip = 1, ninti
            work(ip) = &
                  sum ( matmul( dan(ip,:), alphatg(ip,:,:) ) * dphidx(ip,M,:) )
          end do
!         phi(:,N) * dan(:,i) * alpha(:,i,j) * dphidx(:,M,j) *
          elemmat(N,M) = - sum ( phi(:,N)* work )
        end do
      end do

    end if

!   deallocate arrays

    deallocate ( work )
    deallocate ( alphatg )

    call unset_convection_diffusion_eltree

  end subroutine scalar_diffusion2_open_boundary_eltree


! (transposed or Baumann-Oden) open embedded boundary element using an eltree

  subroutine scalar_diffusion_open_boundary_transposed_eltree ( mesh, problem, &
    elgrp, elem, matrix, vector, first, last, coefficients, oldvectors, &
    elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   variables

    integer :: N, M, ip
    integer :: orientation, signfluxterm
    real(dp) :: alpha


!   set globals element

    call set_globals_poisson ( mesh, coefficients, elgrp )
    call set_globals_convection_diffusion ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_convection_diffusion_eltree ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )

!   allocate arrays

    allocate ( work(ninti) )
    allocate ( cbar(ninti) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
    end do

!   set orientation

    call check ( coefficients, &
      'scalar_diffusion_open_boundary_transposed_eltree', &
      indexarray=[2,17,44], minimum=[0,-1,-1], maximum=[0,1,1] )

    orientation = get_coefficient ( coefficients, index=17, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   set sign of flux term (+1=symmetric (Nitsche),-1=Baumann-Oden)

    signfluxterm = get_coefficient ( coefficients, index=44, default=1 )

!   the constant alpha coefficient

    alpha = coefficients%r(1)

!   vector

    if ( vector ) then

      if ( coefficients%i(42) >= 0 ) then

!       determine prescribed scalar cbar

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(42), &
          value=coefficients%r(18), func=coefficients%func1(4)%p, &
          funcnr=coefficients%i(43), x=xg, indx_v=6, &
          layer=layer, phi=phi, indx_e=6, coef=cbar )

        do N = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dan(ip,:) )
          end do
          elemvec(N) = &
               - alpha * sum ( work * cbar ) * real ( signfluxterm, kind=dp )
        end do

      else

        elemvec = 0

      end if

    end if

    if ( matrix ) then

!     flux term

      do N = 1, ndf
        do M = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dan(ip,:) )
          end do
!         phi(:,M) * dphidx(:,N,i) * dan(:,i)
          elemmat(N,M) = - alpha * sum ( work * phi(:,M) ) &
                                  * real ( signfluxterm, kind=dp )
        end do
      end do

    end if

!   deallocate arrays

    deallocate ( work )
    deallocate ( cbar )

    call unset_convection_diffusion_eltree

  end subroutine scalar_diffusion_open_boundary_transposed_eltree


! Embedded (weak) Dirichlet boundary condition for a scalar using an eltree

  subroutine scalar_diffusion_embedded_dirichlet_eltree ( mesh, problem, &
    elgrp, elem, matrix, vector, first, last, coefficients, oldvectors, &
    elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

!   variables

    integer :: i, k, N, M, ip, orientation, info
    real(dp) :: kappa


!   first build scalar mass matrix (volume integration)

!   set globals, gauss, shapefunctions, ...

    call set_convection_diffusion_elem ( mesh, problem, elgrp, elem, &
      first=.true., last=.true., coefficients=coefficients, &
      oldvectors=oldvectors )

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    xg =  matmul ( phi, x(1:ndf,:) )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

!   allocate mass matrix

    allocate ( work2(ndf,ndf) )

!   fill lower triangle of mass matrix

    do N = 1, ndf
      do M = 1, N
        work2(N,M) = sum ( phi(:,N) * phi(:,M) * wg * detF )
      end do
    end do

!   Cholesky factorization of mass matrix

    call dpotrf ( 'L', ndf, work2, ndf, info )  ! Cholesky factorization

    if ( info /= 0 ) then
      write(*,'(/a/a/a,i0/2(a,i0)/)') &
        'Error in scalar_diffusion_embedded_dirichlet_eltree: ', &
        ' Cholesky factorization failed', &
        ' info = ', info, &
        ' elgrp = ', elgrp, ' elem = ', elem
      stop
    end if

!   unset volume integration (only deallocates memory)

    call unset_convection_diffusion_elem &
                  ( last=.true., coefficients=coefficients )


!   boundary integral part

!   set globals eltree

    call set_convection_diffusion_eltree ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )

!   allocate arrays

    allocate ( work6(ndf,ndf,ndim), work7(ndf,ndf,ndim) )
    allocate ( work4(ndf,ndim), cbar(ninti) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
    end do

!   set orientation (this does not have an effect)

    call check ( coefficients, &
      'scalar_diffusion_embedded_dirichlet_eltree', &
      indexarray=[17], minimum=[-1], maximum=[1] )

    orientation = get_coefficient ( coefficients, index=17, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   kappa parameter

    kappa = coefficients%r(19)

!   (phi_N,phi_M n)

    do N = 1, ndf
      do M = N, ndf
        do i = 1, ndim
          work6(N,M,i) = sum ( phi(:,N) * phi(:,M) * dan(:,i) )
          work6(M,N,i) = work6(N,M,i) ! symmetry
        end do
      end do
    end do

!   save matrix

    work7 = work6

!   M^-1 * (phi_N,phi_M n)  or [(phi_N,phi_M n) * M^-1]^T  due to symmetry

    do i = 1, ndim

      call dpotrs ( 'L', ndf, ndf, work2, ndf, work6(:,:,i), ndf, info )

      if ( info /= 0 ) then
        write(*,'(/a/4(a,i0/))') &
          'Error in scalar_diffusion_embedded_dirichlet_eltree: ', &
          ' Cholesky solve failed. dim = ', i, &
          ' info = ', info, &
          ' elgrp = ', elgrp, &
          ' elem = ', elem
        stop
      end if

    end do

!   vector

    if ( vector ) then

      if ( coefficients%i(42) >= 0 ) then

!       determine prescribed scalar cbar

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(42), &
          value=coefficients%r(18), func=coefficients%func1(4)%p, &
          funcnr=coefficients%i(43), x=xg, indx_v=6, &
          layer=layer, phi=phi, indx_e=6, coef=cbar )

!       (phi_N n, c)

        do N = 1, ndf
          do i = 1, ndim
            work4(N,i) = sum ( phi(:,N) * dan(:,i) * cbar(:) )
          end do
        end do

!       c_N (phi_N,phi_k n) * M_km^-1 . ( phi_m n, c )

        elemvec = 0
        do k = 1, ndim
          elemvec = elemvec + matmul ( work4(:,k), work6(:,:,k) )
        end do

        elemvec = kappa * elemvec

      else

        elemvec = 0

      end if

    end if

!   matrix

    if ( matrix ) then

!     (phi_N n,phi_k) * M_km^-1 . (phi_m,phi_M n)

      elemmat = 0
      do k = 1, ndim
        elemmat = elemmat + matmul ( work7(:,:,k), work6(:,:,k) )
      end do

      elemmat = kappa * elemmat

    end if

!   deallocate arrays

    deallocate ( work2 )
    deallocate ( work6, work7 )
    deallocate ( work4, cbar )

    call unset_convection_diffusion_eltree

  end subroutine scalar_diffusion_embedded_dirichlet_eltree


! Nitsche term using an eltree

  subroutine scalar_diffusion_Nitsche_eltree ( mesh, problem, elgrp, &
    elem, matrix, vector, first, last, coefficients, oldvectors, &
    elemmat, elemvec )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   variables

    integer :: N, M, ip
    real(dp) :: factor, he


!   set globals element

    call set_globals_poisson ( mesh, coefficients, elgrp )
    call set_globals_convection_diffusion ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_convection_diffusion_eltree ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )


!   allocate arrays

    allocate ( cbar(ninti), da(ninti) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
      da(ip) = sqrt(dot_product(dan(ip,:),dan(ip,:)))
    end do

!   element length scaling

    he = 1 ! no elementscaling yet

!   Nitsche parameter

    factor = coefficients%r(20)

!   vector

    if ( vector ) then

      if ( coefficients%i(42) >= 0 ) then

!       determine prescribed scalar cbar

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(42), &
          value=coefficients%r(18), func=coefficients%func1(4)%p, &
          funcnr=coefficients%i(43), x=xg, indx_v=6, &
          layer=layer, phi=phi, indx_e=6, coef=cbar )

        do N = 1, ndf
          elemvec(N) = sum ( phi(:,N) * cbar * da )
        end do

        elemvec = factor * elemvec / he

      else

        elemvec = 0

      end if

    end if

    if ( matrix ) then

      do N = 1, ndf
        do M = N, ndf
          elemmat(N,M) = sum ( phi(:,N) * phi(:,M) * da )
          elemmat(M,N) = elemmat(N,M) ! symmetry
        end do
      end do

      elemmat = factor * elemmat / he

    end if

!   deallocate arrays

    deallocate ( cbar, da )

    call unset_convection_diffusion_eltree

  end subroutine scalar_diffusion_Nitsche_eltree


! set preamble for integration on the interface defined by an eltree

  subroutine set_convection_diffusion_eltree ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors )

    use convection_diffusion_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors

!   set alias for eltree

    eltree1 => oldvectors%ea(1)%p(elgrp,elem)%p

    if ( .not. associated(eltree1) ) then
      write(*,'(/a/a/2(a,i0)/)') 'Error in set_convection_diffusion_eltree: ', &
        ' no eltree', &
        ' elgrp = ', elgrp, ' elem = ', elem
      stop
    end if

!   Gauss for surface integration

    inttype = coefficients%i(40)

    if ( ndim == 2 ) then
      gauss_ie%globalshape = 'line'
    else if ( ndim == 3 ) then
      gauss_ie%globalshape = 'triangle'
    end if

    if ( coefficients%i(50) == 1 .or. SET_GAUSS_BY_ORDER ) then
      intrule_ie = set_intrule ( gauss_ie%globalshape, inttype, &
                                 order=coefficients%i(41) )
    else
      intrule_ie = coefficients%i(41)
    end if

    gauss_ie%intrule = intrule_ie
    gauss_ie%inttype = inttype

    call set_ninti ( gauss_ie, ninti_ie )

    allocate ( w_ie(ninti_ie), x_ie(ninti_ie,2) )

!   Gauss points and weights of a single interface element

    call set_Gauss_integration ( gauss_ie, x_ie, w_ie )

!   reset ninti

    ninti = number_of_interface_integration_points ( eltree1, ninti=ninti_ie )

!   allocate arrays

    allocate ( xig(ninti,ndim), x(nodalp,ndim) )
    allocate ( wg(ninti), wng(ninti,ndim) )
    allocate ( phi(ninti,ndf) )
    allocate ( detF(ninti), dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
    allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
    allocate ( xg(ninti,ndim) )
    allocate ( dan(ninti,ndim) )

!   integration points and weighting factors for the interface

    call interface_integration_points ( eltree1, xig, wg, wng, ninti=ninti_ie, &
      xe=x_ie, we=w_ie )

!   set shape function

    call set_shape_function ( shapefunc, xig, phi, dphi )

  end subroutine set_convection_diffusion_eltree


! unset preamble for integration on the interface defined by an eltree

  subroutine unset_convection_diffusion_eltree ( )

    use convection_diffusion_globals_m

!   deallocate arrays

    deallocate ( w_ie, x_ie )

    deallocate ( xig, x )
    deallocate ( wg, wng )
    deallocate ( phi )
    deallocate ( detF, dphi, F )
    deallocate ( Finv, dphidx )
    deallocate ( xg )
    deallocate ( dan )

  end subroutine unset_convection_diffusion_eltree

end module convection_diffusion_elements_embedded_boundary_m

