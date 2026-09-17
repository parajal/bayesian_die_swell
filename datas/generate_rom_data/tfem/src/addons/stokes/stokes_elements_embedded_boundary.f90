
! Copyright (C) 2010-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for embedded boundary conditions (Dirichlet, surface
! tension, ... ) for the Stokes equation using an eltree to divide the
! domain and/or introduce a boundary inside an element.

module stokes_elements_embedded_boundary_m

  use tfem_elem_m
  use stokes_set_globals_m
  use stokes_elements_generic_m, only: set_stokes_shape_function_global

  implicit none

contains


! open (embedded) boundary element using an eltree

  subroutine stokes_open_boundary_eltree ( mesh, problem, elgrp, elem, matrix, &
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


!   variables

    integer :: i, j, N, M, ip, orientation
    real(dp) :: eta


!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( work(ninti), work1(ninti) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
    end do

!   set orientation

    call check ( coefficients, 'stokes_open_boundary_eltree', &
      indexarray=[36], minimum=[-1], maximum=[1] )

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   vector

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     viscous traction term

      eta = coefficients%r(1)

!     diagonal blocks

      do N = 1, ndf
        do M = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,M,:) * dan(ip,:) )
          end do
          do i = 1, ndim
            work1 = phi(:,N) * ( work + dphidx(:,M,i) * dan(:,i) )
            elemmat( pos(N,i), pos(M,i) ) = - eta * sum ( work1 )
          end do
        end do
      end do

!     off-diagonal blocks

      do N = 1, ndf
        do M = 1, ndf
          do i = 1, ndim
            do j = 1, ndim
              if ( i == j ) cycle
              work = phi(:,N) * dphidx(:,M,i) * dan(:,j)
              elemmat( pos(N,i), pos(M,j) ) = - eta * sum ( work )
            end do
          end do
        end do
      end do

!     velocity-pressure part

      do N = 1, ndf
        do M = 1, ndfp
          do i = 1, ndim
            elemmat( pos(N,i), posp(M) ) = &
                        sum ( phi(:,N) * psi(:,M) * dan(:,i) )
          end do
        end do
      end do

!     pressure-pressure part

      elemmat( posp, posp ) = 0

!     pressure-velocity part

      do i = 1, ndim
        elemmat( posp, pos(:,i) ) = 0
      end do

    end if

!   deallocate arrays

    deallocate ( work, work1 )

    call unset_stokes_eltree

  end subroutine stokes_open_boundary_eltree


! (transposed or Baumann-Oden) open embedded boundary element using an eltree

  subroutine stokes_open_boundary_transposed_eltree ( mesh, problem, elgrp, &
    elem, matrix, vector, first, last, coefficients, oldvectors, &
    elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   variables

    integer :: i, j, N, M, ip
    integer :: orientation, signviscousterm, signpressureterm
    real(dp) :: eta


!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( work(ninti), work1(ninti) )
    allocate ( ubar(ninti,ndim) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
    end do

!   set orientation

    call check ( coefficients, 'stokes_open_boundary_transposed_eltree', &
      indexarray=[36,44,45], minimum=[-1,-1,-1], maximum=[1,1,1] )

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   set sign of viscous term (+1=symmetric (Nitsche),-1=Baumann-Oden)

    signviscousterm = get_coefficient ( coefficients, index=44, default=1 )

!   set sign of pressure term (+1,-1=sign, 0=pressure term is absent)

    signpressureterm = get_coefficient ( coefficients, index=45, default=0 )

!   viscosity

    eta = coefficients%r(1)

!   vector

    if ( vector ) then

      if ( coefficients%i(42) >= 0 ) then

!       determine prescribed velocity ubar

        call evaluate_vector_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(42), &
          value=coefficients%r(12:11+ndim), vfunc=coefficients%vfunc1(1)%p, &
          vfuncnr=coefficients%i(43), x=xg, indx_v=1, &
          layer=layer, phi=phi, indx_e=1, coef=ubar )

!       velocity part

        do N = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dan(ip,:) )
          end do
          do ip = 1, ninti
            work1(ip) = sum ( dphidx(ip,N,:) * ubar(ip,:) )
          end do
          do i = 1, ndim
            elemvec( pos(N,i) ) = &
               - eta * ( sum ( dan(:,i) * work1 ) + sum ( work * ubar(:,i) ) ) &
                    * real ( signviscousterm, kind=dp )
          end do
        end do

        if ( signpressureterm /= 0 ) then

!         pressure part

          do ip = 1, ninti
            work1(ip) = sum ( dan(ip,:) * ubar(ip,:) )
          end do

          do N = 1, ndfp
            elemvec( posp(N) ) = &
              sum ( psi(:,N) * work1 ) * real ( signpressureterm, kind=dp )
          end do

        else

          elemvec( posp ) = 0

        end if

      else

        elemvec = 0

      end if

    end if

    if ( matrix ) then

!     viscous traction term

!     diagonal blocks

      do N = 1, ndf
        do M = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dan(ip,:) )
          end do
          do i = 1, ndim
            work1 = phi(:,M) * ( work + dphidx(:,N,i) * dan(:,i) )
            elemmat( pos(N,i), pos(M,i) ) = &
                   - eta * sum ( work1 ) * real ( signviscousterm, kind=dp )
          end do
        end do
      end do

!     off-diagonal blocks

      do N = 1, ndf
        do M = 1, ndf
          do i = 1, ndim
            do j = 1, ndim
              if ( i == j ) cycle
              work = phi(:,M) * dphidx(:,N,j) * dan(:,i)
              elemmat( pos(N,i), pos(M,j) ) = &
                    - eta * sum ( work ) * real ( signviscousterm, kind=dp )
            end do
          end do
        end do
      end do

!     velocity-pressure part

      do i = 1, ndim
        elemmat( pos(:,i), posp ) = 0
      end do

!     pressure-pressure part

      elemmat( posp, posp ) = 0

      if ( signpressureterm /= 0 ) then

!       pressure-velocity part

        do N = 1, ndfp
          do M = 1, ndf
            do j = 1, ndim
              elemmat( posp(N), pos(M,j) ) = &
                     sum ( psi(:,N) * phi(:,M) * dan(:,j) ) &
                        * real ( signpressureterm, kind=dp )
            end do
          end do
        end do

      else

        do j = 1, ndim
          elemmat( posp, pos(:,j) ) = 0
        end do

      end if

    end if

!   deallocate arrays

    deallocate ( work, work1 )
    deallocate ( ubar )

    call unset_stokes_eltree

  end subroutine stokes_open_boundary_transposed_eltree


! Embedded (weak) Dirichlet boundary condition for velocity using an eltree

  subroutine stokes_embedded_dirichlet_eltree ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

!   variables

    integer :: i, j, k, N, M, ip, orientation, info
    real(dp) :: kappa, tr


!   first build scalar mass matrix (volume integration)

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first=.true., &
      last=.true., coefficients=coefficients, oldvectors=oldvectors )

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
        'Error in stokes_embedded_dirichlet_eltree: ', &
        ' Cholesky factorization failed', &
        ' info = ', info, &
        ' elgrp = ', elgrp, ' elem = ', elem
      stop
    end if

!   unset volume integration (only deallocates memory)

    call unset_stokes_elem ( last=.true., coefficients=coefficients )


!   boundary integral part

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( work4(ndf,ndim), work5(ndf,ndf,ndim,ndim) )
    allocate ( work6(ndf,ndf,ndim), work7(ndf,ndf,ndim) )
    allocate ( work9(ndf,ndim,ndim), ubar(ninti,ndim) )

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

    call check ( coefficients, 'stokes_embedded_dirichlet_eltree', &
      indexarray=[36], minimum=[-1], maximum=[1] )

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   kappa parameter

    kappa = coefficients%r(11)

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
        write(*,'(/a/4(a,i0/))') 'Error in stokes_embedded_dirichlet_eltree: ',&
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

!       determine prescribed velocity ubar

        call evaluate_vector_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(42), &
          value=coefficients%r(12:11+ndim), vfunc=coefficients%vfunc1(1)%p, &
          vfuncnr=coefficients%i(43), x=xg, indx_v=1, &
          layer=layer, phi=phi, indx_e=1, coef=ubar )

!       (phi_N, nu + un)

        do N = 1, ndf
          do i = 1, ndim
            do j = 1, ndim
              work9(N,i,j) = sum ( phi(:,N) * ( dan(:,i) * ubar(:,j) + &
                                                ubar(:,i) * dan(:,j) ) )
            end do
          end do
        end do

!       v_N (phi_N n,phi_k) * M_km^-1 : ( phi_m, nu+un )

        do i = 1, ndim
          work4(:,i) = 0
          do k = 1, ndim
            work4(:,i) = work4(:,i) + matmul ( work9(:,k,i), work6(:,:,k) )
          end do
        end do

        elemvec = kappa * reshape ( work4, [ ndf*ndim ] )

      else

        elemvec = 0

      end if

    end if

!   matrix

    if ( matrix ) then

!     (phi_N n,phi_k) * M_km^-1 * (phi_m,phi_M n)

      do i = 1, ndim
        do j = 1, ndim
          work5(:,:,i,j) = matmul ( work7(:,:,i), work6(:,:,j) )
        end do
      end do

!     fill matrix

      do N = 1, ndf
        do M = 1, ndf

!         first part (full matrix part)

          do i = 1, ndim
            do j = 1, ndim
              elemmat( pos(N,i), pos(M,j) ) = work5(N,M,j,i)
            end do
          end do

!         second part (diagonal block matrix part)

          tr = 0
          do k = 1, ndim
            tr = tr + work5(N,M,k,k)
          end do

          do i = 1, ndim
            elemmat( pos(N,i), pos(M,i) ) = elemmat( pos(N,i), pos(M,i) ) + tr
          end do

        end do
      end do

      elemmat = kappa * elemmat

    end if

!   deallocate arrays

    deallocate ( work2 )
    deallocate ( work4, work5 )
    deallocate ( work6, work7 )
    deallocate ( work9, ubar )

    call unset_stokes_eltree

  end subroutine stokes_embedded_dirichlet_eltree


! Nitsche term using an eltree

  subroutine stokes_Nitsche_eltree ( mesh, problem, elgrp, &
    elem, matrix, vector, first, last, coefficients, oldvectors, &
    elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   variables

    integer :: i, N, M, ip
    real(dp) :: factor, he


!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( ubar(ninti,ndim), da(ninti) )
    allocate ( work2(ndf,ndf) )

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

    factor = coefficients%r(15)

!   vector

    if ( vector ) then

      if ( coefficients%i(42) >= 0 ) then

!       determine prescribed velocity ubar

        call evaluate_vector_coefficient ( mesh, problem, oldvectors, &
          elgrp, elem, choice=coefficients%i(42), &
          value=coefficients%r(12:11+ndim), vfunc=coefficients%vfunc1(1)%p, &
          vfuncnr=coefficients%i(43), x=xg, indx_v=1, &
          layer=layer, phi=phi, indx_e=1, coef=ubar )

        do N = 1, ndf
          do i = 1, ndim
            elemvec( pos(N,i) ) = sum ( phi(:,N) * ubar(:,i) * da )
          end do
        end do

        elemvec = factor * elemvec / he

      else

        elemvec = 0

      end if

    end if

    if ( matrix ) then

!     diagonal blocks

      do N = 1, ndf
        do M = N, ndf
          work2(N,M) = sum ( phi(:,N) * phi(:,M) * da )
          work2(M,N) = work2(N,M) ! symmetry
        end do
      end do

      work2 = factor * work2 / he

      elemmat = 0

      do i = 1, ndim
        elemmat( pos(:,i), pos(:,i) ) = work2
      end do

    end if

!   deallocate arrays

    deallocate ( ubar, da )
    deallocate ( work2 )

    call unset_stokes_eltree

  end subroutine stokes_Nitsche_eltree


! surface tension

  subroutine surface_tension_eltree ( mesh, problem, elgrp, elem, matrix, &
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


!   variables

    integer :: i, j, N, ip
    real(dp) :: gammac


!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( da(ninti) )
    allocate ( normal(ninti,ndim), work2(ndf,ndim) )
    allocate ( work6(ninti,ndf,ndim), tauten(ninti,ndim,ndim) )

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
      normal(ip,:) = dan(ip,:) / da(ip)
    end do

!   surface tension coefficient

    gammac = coefficients%r(19)

    if ( matrix ) then
      write(*,'(/a/a/)') 'Error in surface_tension_eltree:', &
        ' No matrix to build. Call build_system with buildmatrix=.false.'
      stop
      elemmat = 0._dp
    end if

!   vector

    if ( vector ) then

!     - (nabla v)^T:tau, tau=Gamma*(I-nn)

!     diagonal components of tau

      do i = 1, ndim
        tauten(:,i,i) = 1 - normal(:,i)**2
      end do

!     off-diagonal components of tau

      do i = 1, ndim-1
        do j = i+1, ndim
          tauten(:,i,j) = - normal(:,i)*normal(:,j)
          tauten(:,j,i) = tauten(:,i,j)
        end do
      end do

      do ip = 1, ninti
        work6(ip,:,:) = matmul ( dphidx(ip,:,:), tauten(ip,:,:) )
      end do
      if ( coorsys == 1 ) then
        do N = 1, ndf
          work6(:,N,2) = work6(:,N,2) + phi(:,N) / xg(:,2)
        end do
      end if
      do j = 1, ndim
        work2(:,j) = - gammac * matmul ( da, work6(:,:,j) )
      end do
      elemvec = reshape ( work2, [ ndf*ndim ] )

    end if

!   deallocate arrays

    deallocate ( da )
    deallocate ( normal, work2 )
    deallocate ( work6, tauten )

    call unset_stokes_eltree

  end subroutine surface_tension_eltree


! pressure contribution to the drag

  subroutine integration_pressure_drag_eltree ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: ip, i, orientation

!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( u(ndfp), work(ninti) )
    allocate ( work4(ninti,ndim) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
    end do

!   set orientation

    call check ( coefficients, 'integration_pressure_drag_eltree', &
      indexarray=[36], minimum=[-1], maximum=[1] )

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    dan = real ( orientation, kind=dp ) * dan


!   get pressure

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqpress], layer=layer )

    work = matmul ( psi, u )

!   evaluate in each integration point: -p*I.n = -p*n

    do ip = 1, ninti
      work4(ip,:) = - work(ip) * dan(ip,:)
    end do

!   integrate

    do i = 1, ndim
      elemvec(i) = sum ( work4(:,i) )
    end do

!   deallocate arrays

    deallocate ( u, work )
    deallocate ( work4 )

    call unset_stokes_eltree

  end subroutine integration_pressure_drag_eltree


! viscous stress contribution to the drag

  subroutine integration_viscous_drag_eltree ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: ip, i, j, orientation
    real(dp) :: eta

!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( u(ndim*ndf), uvector(ndf,ndim), gradu(ninti,ndim,ndim) )
    allocate ( tauten(ninti,ndim,ndim) )
    allocate ( work4(ninti,ndim) )

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

    call check ( coefficients, 'integration_viscous_drag_eltree', &
      indexarray=[36], minimum=[-1], maximum=[1] )

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    dan = real ( orientation, kind=dp ) * dan


!   get velocity

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ndim] )

!   compute velocity gradient

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

!   compute viscous stress tensor: tau = 2*eta*D

    eta = coefficients%r(1)

    if ( coorsys <= 1 ) then

      tauten(:,1,1) = 2 * eta * gradu(:,1,1)                ! tau_xx
      tauten(:,1,2) = eta * ( gradu(:,1,2) + gradu(:,2,1) ) ! tau_xy
      tauten(:,2,1) = tauten(:,1,2)
      tauten(:,2,2) = 2 * eta * gradu(:,2,2)                ! tau_yy

    else if ( coorsys == 2 ) then

      tauten(:,1,1) = 2 * eta * gradu(:,1,1)                ! tau_xx
      tauten(:,1,2) = eta * ( gradu(:,1,2) + gradu(:,2,1) ) ! tau_xy
      tauten(:,1,3) = eta * ( gradu(:,1,3) + gradu(:,3,1) ) ! tau_xz
      tauten(:,2,1) = tauten(:,1,2)
      tauten(:,2,2) = 2 * eta * gradu(:,2,2)                ! tau_yy
      tauten(:,2,3) = eta * ( gradu(:,2,3) + gradu(:,3,2) ) ! tau_yz
      tauten(:,3,1) = tauten(:,1,3)
      tauten(:,3,2) = tauten(:,2,3)
      tauten(:,3,3) = 2 * eta * gradu(:,3,3)                ! tau_zz

    end if

!   evaluate in each integration point: 2*eta*D.n

    do ip = 1, ninti
      work4(ip,:) = matmul( tauten(ip,:,:), dan(ip,:) )
    end do

!   integrate

    do i = 1, ndim
      elemvec(i) = sum ( work4(:,i) )
    end do

!   deallocate arrays

    deallocate ( u, uvector, gradu )
    deallocate ( tauten )
    deallocate ( work4 )

    call unset_stokes_eltree

  end subroutine integration_viscous_drag_eltree


! set preamble for integration on the interface defined by an eltree

  subroutine set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use stokes_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors

!   variables

    integer :: i

!   set alias for eltree

    eltree1 => oldvectors%ea(1)%p(elgrp,elem)%p

    if ( .not. associated(eltree1) ) then
      write(*,'(/a/a/2(a,i0)/)') 'Error in set_stokes_eltree: ', &
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

    if ( coefficients%i(77) == 1 .or. SET_GAUSS_BY_ORDER ) then
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

    allocate ( pos(ndf,ndim), posp(ndfp) )
    allocate ( xig(ninti,ndim), x(nodalp,ndim) )
    allocate ( wg(ninti), wng(ninti,ndim) )
    allocate ( phi(ninti,ndf), psi(ninti,ndfp) )
    allocate ( detF(ninti), dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
    allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
    allocate ( xg(ninti,ndim) )
    allocate ( dan(ninti,ndim) )

!   position arrays
    pos  = reshape ( [ ( i, i = 1, ndf*ndim ) ], [ ndf, ndim ] )
    posp = [ ( i, i = 1, ndfp ) ] + ndf*ndim

!   integration points and weighting factors for the interface

    call interface_integration_points ( eltree1, xig, wg, wng, ninti=ninti_ie, &
      xe=x_ie, we=w_ie )

!   set shape functions

    call set_shape_function ( shapefunc, xig, phi, dphi ) ! velocity
    if ( coefficients%i(62) == 0 ) then
      call set_shape_function ( shapefuncp, xig, psi ) ! pressure
    end if

  end subroutine set_stokes_eltree


! unset preamble for integration on the interface defined by an eltree

  subroutine unset_stokes_eltree ( )

    use stokes_globals_m

!   deallocate arrays

    deallocate ( w_ie, x_ie )

    deallocate ( pos, posp )
    deallocate ( xig, x )
    deallocate ( wg, wng )
    deallocate ( phi, psi )
    deallocate ( detF, dphi, F )
    deallocate ( Finv, dphidx )
    deallocate ( xg )
    deallocate ( dan )

  end subroutine unset_stokes_eltree

end module stokes_elements_embedded_boundary_m

