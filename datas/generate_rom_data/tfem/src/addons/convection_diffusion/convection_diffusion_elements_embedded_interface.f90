
! Copyright (C) 2011-2012 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for embedded interface conditions for the
! convection-diffusion equation using an eltree to divide the domain
! and/or introduce an interface.

module convection_diffusion_elements_embedded_interface_m

  use tfem_elem_m
  use convection_diffusion_elements_embedded_boundary_m, only: &
      set_convection_diffusion_eltree, unset_convection_diffusion_eltree
  use convection_diffusion_set_globals_m

  implicit none


contains


! embedded interface flux element using an eltree for
! the scalar diffusion equation with varying coefficient alpha
!   - nabla ( alpha nabla c ) = f
! Coefficient alpha is a scalar.
! NOTE: order='DN' in build_system, which is the default if physical quantities
! have been defined.

  subroutine scalar_diffusion_flux_interface_eltree ( mesh, problem, &
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

    integer :: N, M, ip, orientation, clayer, i


!   set globals element

    call set_globals_poisson ( mesh, coefficients, elgrp )
    call set_globals_convection_diffusion ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_convection_diffusion_eltree ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )

!   allocate arrays

    allocate ( work(ninti), alphag(ninti), jbar(ninti) )
    allocate ( da(ninti) )
    allocate ( pos1(ndf), pos2(ndf) )

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

!   set orientation (normal directs outside of the contributing side)

    call check ( coefficients, 'scalar_diffusion_flux_interface_eltree', &
      indexarray=[17,49], minimum=[-1,0], maximum=[1,2] )

    clayer = get_coefficient ( coefficients, index=49, default=1 )

    select case ( clayer )
    case(1)
      pos1 = [(i,i=1,ndf)]
      pos2 = pos1 + ndf
    case(2)
      pos2 = [(i,i=1,ndf)]
      pos1 = pos2 + ndf
    case default
      call errormsg_case_default ( 'scalar_diffusion_flux_interface_eltree', &
        'clayer', int_value=clayer )
    end select

    orientation = get_coefficient ( coefficients, index=17, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   vector

    if ( vector ) then

      if ( coefficients%i(47) >= 0 ) then

!       determine prescribed scalar jbar

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(47), &
          value=coefficients%r(22), func=coefficients%func1(6)%p, &
          funcnr=coefficients%i(48), x=xg, indx_v=8, &
          layer=layer, phi=phi, indx_e=8, coef=jbar )

        do N = 1, ndf
          elemvec(pos1(N)) = sum ( phi(:,N) * jbar * da )
        end do

        elemvec(pos2) = 0

      else

        elemvec = 0

      end if

    end if

    if ( matrix ) then

!     interface flux term

      call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
        elgrp, elem, choice=coefficients%i(2), value=coefficients%r(1), &
        func=coefficients%func1(1)%p, funcnr=coefficients%i(3), x=xg, indx_v=1,&
        layer=layer, phi=phi, indx_e=1, coef=alphag )

!     contribution to the first side

      do N = 1, ndf
        do M = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,M,:) * dan(ip,:) )
          end do
!         alpha_g * phi(:,N) * dphidx(:,M,i) * dan(:,i)
          elemmat(pos1(N),pos1(M)) = - sum ( alphag * phi(:,N)* work )
        end do
      end do

!     contribution to the second side

      elemmat(pos2,pos1) = - elemmat(pos1,pos1)

!     fill zero to the contribution from the other side (needs separate build)

      elemmat(:,pos2) = 0

    end if

!   deallocate arrays

    deallocate ( work, alphag, jbar )
    deallocate ( da )
    deallocate ( pos1, pos2 )

    call unset_convection_diffusion_eltree

  end subroutine scalar_diffusion_flux_interface_eltree


! embedded interface flux element using an eltree for
! the scalar diffusion equation with varying coefficient alpha
!   - nabla ( alpha nabla c ) = f
! Coefficient alpha is a tensor.
! NOTE: order='DN' in build_system, which is the default if physical quantities
! have been defined.

  subroutine scalar_diffusion2_flux_interface_eltree ( mesh, problem, &
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

    integer :: N, M, ip, orientation, clayer, i


!   set globals element

    call set_globals_poisson ( mesh, coefficients, elgrp )
    call set_globals_convection_diffusion ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_convection_diffusion_eltree ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )

!   allocate arrays

    allocate ( work(ninti) )
    allocate ( alphatg(ninti,ndim,ndim) )
    allocate ( da(ninti) )
    allocate ( jbar(ninti), pos1(ndf), pos2(ndf) )

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

!   set orientation (normal directs outside of the contributing side)

    call check ( coefficients, 'scalar_diffusion2_flux_interface_eltree', &
      indexarray=[17,49], minimum=[-1,0], maximum=[1,2] )

    clayer = get_coefficient ( coefficients, index=49, default=1 )

    select case ( clayer )
    case(1)
      pos1 = [(i,i=1,ndf)]
      pos2 = pos1 + ndf
    case(2)
      pos2 = [(i,i=1,ndf)]
      pos1 = pos2 + ndf
    case default
      call errormsg_case_default ( 'scalar_diffusion2_flux_interface_eltree', &
        'clayer', int_value=clayer )
    end select

    orientation = get_coefficient ( coefficients, index=17, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   vector

    if ( vector ) then

      if ( coefficients%i(47) >= 0 ) then

!       determine prescribed scalar jbar

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(47), &
          value=coefficients%r(22), func=coefficients%func1(6)%p, &
          funcnr=coefficients%i(48), x=xg, indx_v=8, &
          layer=layer, phi=phi, indx_e=8, coef=jbar )

        do N = 1, ndf
          elemvec(pos1(N)) = sum ( phi(:,N) * jbar * da )
        end do

        elemvec(pos2) = 0

      else

        elemvec = 0

      end if

    end if

    if ( matrix ) then

!     boundary flux term

      call evaluate_tensor_coefficient ( mesh, problem, oldvectors, &
        elgrp, elem, choice=coefficients%i(2), &
        value=reshape(coefficients%r(9:8+ndim**2), [ndim,ndim]), &
        tfunc=coefficients%tfunc1(1)%p, tfuncnr=coefficients%i(3), x=xg, &
        indx_v=1, layer=layer, phi=phi, indx_e=1, coef=alphatg )

!     contribution to the first side

      do N = 1, ndf
        do M = 1, ndf
          do ip = 1, ninti
            work(ip) = &
                  sum ( matmul( dan(ip,:), alphatg(ip,:,:) ) * dphidx(ip,M,:) )
          end do
!         phi(:,N) * dan(:,i) * alpha(:,i,j) * dphidx(:,M,j) *
          elemmat(pos1(N),pos1(M)) = - sum ( phi(:,N)* work )
        end do
      end do

!     contribution to the second side

      elemmat(pos2,pos1) = - elemmat(pos1,pos1)

!     fill zero to the contribution from the other side (needs separate build)

      elemmat(:,pos2) = 0

    end if

!   deallocate arrays

    deallocate ( work )
    deallocate ( alphatg )
    deallocate ( jbar, pos1, pos2 )
    deallocate ( da )

    call unset_convection_diffusion_eltree

  end subroutine scalar_diffusion2_flux_interface_eltree


! (transposed or Baumann-Oden) embedded interface element using an eltree
! NOTE: order='DN' in build_system, which is the default if physical quantities
! have been defined.

  subroutine scalar_diffusion_flux_interface_transposed_eltree ( mesh, problem,&
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

    integer :: N, M, ip, clayer, i
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
    allocate ( ibar(ninti) )
    allocate ( pos1(ndf), pos2(ndf) )

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

!   set orientation (normal directs outside of the contributing side)

    call check ( coefficients, &
      'scalar_diffusion_flux_interface_transposed_eltree', &
      indexarray=[2,17,44,49], minimum=[0,-1,-1,0], maximum=[0,1,1,2] )

    clayer = get_coefficient ( coefficients, index=49, default=1 )

    select case ( clayer )
    case(1)
      pos1 = [(i,i=1,ndf)]
      pos2 = pos1 + ndf
    case(2)
      pos2 = [(i,i=1,ndf)]
      pos1 = pos2 + ndf
    case default
      call errormsg_case_default ( &
        'scalar_diffusion_flux_interface_transposed_eltree', &
        'clayer', int_value=clayer )
    end select

    orientation = get_coefficient ( coefficients, index=17, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   set sign of flux term (+1=symmetric (Nitsche),-1=Baumann-Oden)

    signfluxterm = get_coefficient ( coefficients, index=44, default=1 )

!   the constant alpha coefficient

    alpha = coefficients%r(1)

!   vector

    if ( vector ) then

      if ( coefficients%i(45) >= 0 ) then

!       determine prescribed scalar ibar

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(45), &
          value=coefficients%r(21), func=coefficients%func1(5)%p, &
          funcnr=coefficients%i(46), x=xg, indx_v=7, &
          layer=layer, phi=phi, indx_e=7, coef=ibar )

        do N = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dan(ip,:) )
          end do
          elemvec(pos1(N)) = &
                 alpha * sum ( work * ibar ) * real ( signfluxterm, kind=dp )
        end do

        elemvec(pos2) = 0

      else

        elemvec = 0

      end if

    end if

    if ( matrix ) then

!     flux term

!     contribution of the first side

      do N = 1, ndf
        do M = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dan(ip,:) )
          end do
!         phi(:,M) * dphidx(:,N,i) * dan(:,i)
          elemmat(pos1(N),pos1(M)) = - alpha * sum ( work * phi(:,M) ) &
                                         * real ( signfluxterm, kind=dp )
        end do
      end do

!     contribution of the second side

      elemmat(pos1,pos2) = - elemmat(pos1,pos1)

!     fill zero the contribution to the other side (needs separate build)

      elemmat(pos2,:) = 0

    end if

!   deallocate arrays

    deallocate ( work )
    deallocate ( ibar )
    deallocate ( pos1, pos2 )

    call unset_convection_diffusion_eltree

  end subroutine scalar_diffusion_flux_interface_transposed_eltree


! Embedded interface condition for scalar using an eltree
! NOTE: order='DN' in build_system, which is the default if physical quantities
! have been defined.

  subroutine scalar_diffusion_embedded_interface_eltree ( mesh, problem, &
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
        'Error in scalar_diffusion_embedded_interface_eltree: ', &
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
    allocate ( work4(ndf,ndim), ibar(ninti) )

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
      'scalar_diffusion_embedded_interface_eltree', &
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
          'Error in scalar_diffusion_embedded_interface_eltree: ', &
          ' Cholesky solve failed. dim = ', i, &
          ' info = ', info, &
          ' elgrp = ', elgrp, &
          ' elem = ', elem
        stop
      end if

    end do

!   vector

    if ( vector ) then

      if ( coefficients%i(45) >= 0 ) then

!       determine prescribed scalar ibar

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(45), &
          value=coefficients%r(21), func=coefficients%func1(5)%p, &
          funcnr=coefficients%i(46), x=xg, indx_v=7, &
          layer=layer, phi=phi, indx_e=7, coef=ibar )

!       (phi_N n, ibar)

        do N = 1, ndf
          do i = 1, ndim
            work4(N,i) = sum ( phi(:,N) * dan(:,i) * ibar(:) )
          end do
        end do

!       c_N (phi_N,phi_k n) * M_km^-1 . ( phi_m n, ibar )

        elemvec(1:ndf) = 0
        do k = 1, ndim
          elemvec(1:ndf) = elemvec(1:ndf) + matmul ( work4(:,k), work6(:,:,k) )
        end do

        elemvec(1:ndf) = - kappa * elemvec(1:ndf)
        elemvec(ndf+1:2*ndf) = - elemvec(1:ndf)

      else

        elemvec = 0

      end if

    end if

!   matrix

    if ( matrix ) then

!     (phi_N n,phi_k) * M_km^-1 . (phi_m,phi_M n)

      elemmat(1:ndf,1:ndf) = 0
      do k = 1, ndim
        elemmat(1:ndf,1:ndf) = elemmat(1:ndf,1:ndf) &
                                  + matmul ( work7(:,:,k), work6(:,:,k) )
      end do

      elemmat(1:ndf,1:ndf) = kappa * elemmat(1:ndf,1:ndf)
      elemmat(ndf+1:2*ndf,1:ndf) = - elemmat(1:ndf,1:ndf)
      elemmat(:,ndf+1:2*ndf) = - elemmat(:,1:ndf)

    end if

!   deallocate arrays

    deallocate ( work2 )
    deallocate ( work6, work7 )
    deallocate ( work4, ibar )

    call unset_convection_diffusion_eltree

  end subroutine scalar_diffusion_embedded_interface_eltree


! Nitsche term using an eltree
! NOTE: order='DN' in build_system, which is the default if physical quantities
! have been defined.

  subroutine scalar_diffusion_Nitsche_interface_eltree ( mesh, problem, elgrp, &
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

    allocate ( ibar(ninti), da(ninti) )

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

      if ( coefficients%i(45) >= 0 ) then

!       determine prescribed scalar ibar

        call evaluate_scalar_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(45), &
          value=coefficients%r(21), func=coefficients%func1(5)%p, &
          funcnr=coefficients%i(46), x=xg, indx_v=7, &
          layer=layer, phi=phi, indx_e=7, coef=ibar )

        do N = 1, ndf
          elemvec(N) = sum ( phi(:,N) * ibar * da )
        end do

        elemvec(1:ndf) = - factor * elemvec(1:ndf) / he
        elemvec(ndf+1:2*ndf) = - elemvec(1:ndf)

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

      elemmat(1:ndf,1:ndf) = factor * elemmat(1:ndf,1:ndf) / he
      elemmat(ndf+1:2*ndf,1:ndf) = - elemmat(1:ndf,1:ndf)
      elemmat(:,ndf+1:2*ndf) = - elemmat(:,1:ndf)

    end if

!   deallocate arrays

    deallocate ( ibar, da )

    call unset_convection_diffusion_eltree

  end subroutine scalar_diffusion_Nitsche_interface_eltree


! Embedded interface convective flux element using an eltree
! Lesaint-Raviart type of upwinding along the interface.
! This element builds the contribution from both sides at once.
! NOTE: the coefficient gamma and the normal velocity un of the interface are
! determined from one side only and it is assumed that gamma*un is continuous
! across the interface.
! NOTE: the normal vector must be directed from the first layer to the second
! layer.
! NOTE: order='DN' in build_system, which is the default if physical quantities
! have been defined.

  subroutine scalar_ugradc_interface_eltree ( mesh, problem, &
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

    allocate ( gammag(ninti), uvecg(ninti,ndim) )
    allocate ( ung(ninti), work(ninti) )


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

!   set orientation (normal directs outside of the contributing side)

    call check ( coefficients, 'scalar_ugradc_interface_eltree', &
      indexarray=[17], minimum=[-1], maximum=[1] )

    orientation = get_coefficient ( coefficients, index=17, default=1 )

    dan = real ( orientation, kind=dp ) * dan

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

!   normal velocity component

    do ip = 1, ninti
      ung(ip) = dot_product ( uvecg(ip,:), dan(ip,:) )
    end do

!   vector

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     contribution from the first side

      do N = 1, ndf
        do M = 1, ndf
          elemmat(N,M) = &
                  - sum ( phi(:,N) * gammag * ung * phi(:,M), mask = ung < 0 )
          elemmat(N+ndf,M) = &
                  - sum ( phi(:,N) * gammag * ung * phi(:,M), mask = ung > 0 )
        end do
      end do

!     second side

      elemmat(:,ndf+1:2*ndf) = - elemmat(:,1:ndf)

    end if

!   deallocate arrays

    deallocate ( gammag, uvecg )
    deallocate ( ung, work )

    call unset_convection_diffusion_eltree

  end subroutine scalar_ugradc_interface_eltree


! Embedded interface convective flux element using an eltree
! Lesaint-Raviart type of upwinding along the interface.
! This element builds the contribution from one side only and needs to be called
! separately for the second side.
! NOTE: order='DN' in build_system, which is the default if physical quantities
! have been defined.

  subroutine scalar_ugradc_interface_eltree2 ( mesh, problem, &
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

    integer :: N, M, ip, orientation, clayer, i


!   set globals element

    call set_globals_poisson ( mesh, coefficients, elgrp )
    call set_globals_convection_diffusion ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_convection_diffusion_eltree ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors )

!   allocate arrays

    allocate ( gammag(ninti), uvecg(ninti,ndim) )
    allocate ( ung(ninti), work(ninti) )
    allocate ( pos1(ndf), pos2(ndf) )


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

!   set orientation (normal directs outside of the contributing side)

    call check ( coefficients, 'scalar_ugradc_interface_eltree2', &
      indexarray=[17,49], minimum=[-1,0], maximum=[1,2] )

    clayer = get_coefficient ( coefficients, index=49, default=1 )

    select case ( clayer )
    case(1)
      pos1 = [(i,i=1,ndf)]
      pos2 = pos1 + ndf
    case(2)
      pos2 = [(i,i=1,ndf)]
      pos1 = pos2 + ndf
    case default
      call errormsg_case_default ( 'scalar_ugradc_interface_eltree2', &
        'clayer', int_value=clayer )
    end select

    orientation = get_coefficient ( coefficients, index=17, default=1 )

    dan = real ( orientation, kind=dp ) * dan

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

!   normal velocity component

    do ip = 1, ninti
      ung(ip) = dot_product ( uvecg(ip,:), dan(ip,:) )
    end do

!   vector

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     contributing side adds to the side depending on the sign of ung

      do N = 1, ndf
        do M = 1, ndf
          elemmat(pos1(N),pos1(M)) = &
                  - sum ( phi(:,N) * gammag * ung * phi(:,M), mask = ung < 0 )
          elemmat(pos2(N),pos1(M)) = &
                  - sum ( phi(:,N) * gammag * ung * phi(:,M), mask = ung > 0 )
        end do
      end do

!     fill zero to the contribution from the other side (needs separate build)

      elemmat(:,pos2) = 0

    end if

!   deallocate arrays

    deallocate ( gammag, uvecg )
    deallocate ( ung, work )
    deallocate ( pos1, pos2 )

    call unset_convection_diffusion_eltree

  end subroutine scalar_ugradc_interface_eltree2

end module convection_diffusion_elements_embedded_interface_m

