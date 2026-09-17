
! Copyright (C) 2011-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for embedded particle boundary conditions for the
! Stokes equation using an eltree to divide the domain and/or introduce
! an interface.

module stokes_elements_embedded_particle_m

  use tfem_elem_m
  use stokes_elements_embedded_boundary_m, only: set_stokes_eltree, &
      unset_stokes_eltree
  use stokes_set_globals_m

  implicit none

contains


! open (embedded) boundary element using an eltree.
! The elemsub5 interface of build_system_constraint is used.
! A global constraint on an elementset must be used with Lagrange multipliers
! acting as unknown particle velocity and particle rotation rate.
! The diagonal block S is also filled (regular open boundary term).
! No additional unknowns. No diagonal block C.
! NOTE: arguments changed: elem,node -> elem1,elgrp

  subroutine stokes_open_boundary_particle_eltree ( mesh, problem, constr, &
    elem1, elgrp, matrix, vector, first, last, coefficients, oldvectors, &
    elemmatsdiag, elemmatsdiag2, elemmat, elemmatt, elemmat2, elemmat2t, &
    elemmatdiag, elemmatadd, elemmataddt, elemmatdiagadd, elemvecf, &
    elemvecf2, elemvec, elemvecadd )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem1, elgrp
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmatsdiag, elemmatsdiag2
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmatt
    real(dp), intent(out), dimension(:,:) :: elemmat2, elemmat2t
    real(dp), intent(out), dimension(:,:) :: elemmatadd, elemmataddt
    real(dp), intent(out), dimension(:,:) :: elemmatdiag, elemmatdiagadd
    real(dp), intent(out), dimension(:) :: elemvecf, elemvecf2
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


!   variables

    integer :: elementset1, elem
    integer :: i, j, N, M, ip, orientation, signopenterm
    real(dp) :: eta, xp(3)

!   find real element number

    elementset1 = problem%constraints(constr)%elementset1

    if ( elementset1 == 0 ) then
      write(*,'(/a/a/)') &
        'Error in stokes_open_boundary_particle_eltree: ', &
        ' constraint not defined on an elementset'
      stop
    end if

    elem = mesh%elementsets(elementset1)%elements(elgrp)%a(elem1)

!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( work(ninti), work1(ninti) )
    allocate ( work6(ndim,ndim,ndim*ndf+ndfp) )

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

    call check ( coefficients, 'stokes_open_boundary_particle_eltree', &
      indexarray=[36,46], minimum=[-1,-1], maximum=[1,1] )

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   set sign of transposed open boundary term
!   (+1=symmetric (Nitsche), 0=term is absent,-1=Baumann-Oden)

    signopenterm = get_coefficient ( coefficients, index=46, default=0 )

!   particle position

    xp(:ndim) = coefficients%r(16:15+ndim)

!   vector

    if ( vector ) then

      elemvec = 0
      elemvecf = 0

    end if

    if ( matrix ) then

!     viscous traction term

      eta = coefficients%r(1)

!     diagonal blocks of S

      do N = 1, ndf
        do M = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,M,:) * dan(ip,:) )
          end do
          do i = 1, ndim
            work1 = phi(:,N) * ( work + dphidx(:,M,i) * dan(:,i) )
            elemmatsdiag( pos(N,i), pos(M,i) ) = - eta * sum ( work1 )
          end do
        end do
      end do

!     off-diagonal blocks of S

      do N = 1, ndf
        do M = 1, ndf
          do i = 1, ndim
            do j = 1, ndim
              if ( i == j ) cycle
              work = phi(:,N) * dphidx(:,M,i) * dan(:,j)
              elemmatsdiag( pos(N,i), pos(M,j) ) = - eta * sum ( work )
            end do
          end do
        end do
      end do

!     velocity-pressure part of S

      do N = 1, ndf
        do M = 1, ndfp
          do i = 1, ndim
            elemmatsdiag( pos(N,i), posp(M) ) = &
                            sum ( phi(:,N) * psi(:,M) * dan(:,i) )
          end do
        end do
      end do

!     pressure-pressure part of S

      elemmatsdiag( posp, posp ) = 0

!     pressure-velocity part of S

      do i = 1, ndim
        elemmatsdiag( posp, pos(:,i) ) = 0
      end do

    end if

!   rest of the matrix

    if ( matrix .and. coorsys == 1 ) then

!     axisymmetric: only translation along z-direction

      elemmat(1,:) = - sum( elemmatsdiag(pos(:,1),:), dim=1 )

      if ( signopenterm /= 0 ) then
        elemmatsdiag = elemmatsdiag + &
                       transpose(elemmatsdiag) * real ( signopenterm, kind=dp )
        elemmatt = transpose(elemmat) * real ( signopenterm, kind=dp )
      else
        elemmatt = 0  ! A^T=0
      end if

      elemmatdiag = 0

    else if ( matrix ) then

!     now the part of A regarding the translational velocity U (force=0):
!      -(-V,mu(nabla u+nabla U^T).n-pn)_Gamma
!     We use the partition of unity of the shape function: sum_k phi_k = 1, to
!     derive it from
!      -(v,mu(nabla u+nabla U^T).n-pn)_Gamma

      do i = 1, ndim
        elemmat(i,:) = - sum( elemmatsdiag(pos(:,i),:), dim=1 )
      end do

!     now the part of A regarding the rotational rate omega (torque=0):
!      -(-chi x (x-xp),mu(nabla u+nabla U^T).n-pn)_Gamma
!     We use the identity for isoparametric elements : sum_k x_k phi_k = x, to
!     derive it from
!      -(v,mu(nabla u+nabla U^T).n-pn)_Gamma

      do i = 1, ndim
        do j = 1, ndim
          work6(i,j,:) = matmul( x(:,i) - xp(i), elemmatsdiag(pos(:,j),:) )
        end do
      end do

!     take outer product
      if ( ndim == 2 .and. coorsys == 0 ) then
        elemmat(3,:) = - ( work6(1,2,:) - work6(2,1,:) )
      else if ( ndim == 3 ) then
        elemmat(4,:) = - ( work6(2,3,:) - work6(3,2,:) )
        elemmat(5,:) = - ( work6(3,1,:) - work6(1,3,:) )
        elemmat(6,:) = - ( work6(1,2,:) - work6(2,1,:) )
      end if

      if ( signopenterm /= 0 ) then
        elemmatsdiag = elemmatsdiag + &
                       transpose(elemmatsdiag) * real ( signopenterm, kind=dp )
        elemmatt = transpose(elemmat) * real ( signopenterm, kind=dp )
      else
        elemmatt = 0  ! A^T=0
      end if

      elemmatdiag = 0

    end if

!   deallocate arrays

    deallocate ( work, work1 )
    deallocate ( work6 )

    call unset_stokes_eltree

  end subroutine stokes_open_boundary_particle_eltree


! Embedded (weak) boundary element for rigid particle velocity using an eltree.
! The elemsub4 interface of build_system_constraint is used.
! A global constraint on an elementset must be used with Lagrange multipliers
! acting as unknown particle velocity and particle rotation rate.
! The diagonal block S is also filled (regular velocity stabilization term).
! The diagonal block C is also filled, so diagonal_block=.true. must be used.
! No additional unknowns.
! NOTE: arguments changed: elem,node -> elem1,elgrp

  subroutine stokes_embedded_particle_eltree ( mesh, problem, constr, elem1, &
    elgrp, matrix, vector, first, last, coefficients, oldvectors, &
    elemmatsdiag, elemmatsdiag2, elemmat, elemmat2, elemmatdiag, elemmatadd, &
    elemmatdiagadd, elemvec, elemvecadd )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem1, elgrp
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmatsdiag, elemmatsdiag2
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:,:) :: elemmatdiag, elemmatdiagadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


!   variables

    integer :: elementset1, elem
    integer :: i, j, k, N, M, ip, orientation, info
    real(dp) :: kappa, tr, xp(3)


!   find real element number

    elementset1 = problem%constraints(constr)%elementset1

    if ( elementset1 == 0 ) then
      write(*,'(/a/a/)') &
        'Error in stokes_embedded_particle_eltree: ', &
        ' constraint not defined on an elementset'
      stop
    end if

    elem = mesh%elementsets(elementset1)%elements(elgrp)%a(elem1)

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
        'Error in stokes_embedded_particle_eltree: ', &
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
    allocate ( work9(ndf,ndim,ndim) )

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

    call check ( coefficients, 'stokes_embedded_particle_eltree', &
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
        write(*,'(/a/4(a,i0/))') 'Error in stokes_embedded_particle_eltree: ', &
          ' Cholesky solve failed. dim = ', i, &
          ' info = ', info, &
          ' elgrp = ', elgrp, &
          ' elem = ', elem
        stop
      end if

    end do

!   vector

    if ( vector ) then

      elemvec = 0

    end if

!   matrix

    if ( matrix ) then

!     first build matrix S

!     (phi_N n,phi_k) * M_km^-1 * (phi_m,phi_M n)

      do i = 1, ndim
        do j = 1, ndim
          work5(:,:,i,j) = matmul ( work7(:,:,i), work6(:,:,j) )
        end do
      end do

      work5 = kappa * work5

!     fill matrix

      do N = 1, ndf
        do M = 1, ndf

!         first part (full matrix part)

          do i = 1, ndim
            do j = 1, ndim
              elemmatsdiag( pos(N,i), pos(M,j) ) = work5(N,M,j,i)
            end do
          end do

!         second part (diagonal block matrix part)

          tr = 0
          do k = 1, ndim
            tr = tr + work5(N,M,k,k)
          end do

          do i = 1, ndim
            elemmatsdiag( pos(N,i), pos(M,i) ) = &
                            elemmatsdiag( pos(N,i), pos(M,i) ) + tr
          end do

        end do
      end do

!     velocity-pressure part of S

      do i = 1, ndim
        elemmatsdiag( pos(:,i), posp ) = 0
      end do

!     pressure-pressure part of S

      elemmatsdiag( posp, posp ) = 0

!     pressure-velocity part of S

      do i = 1, ndim
        elemmatsdiag( posp, pos(:,i) ) = 0
      end do

    end if

!   deallocate arrays

    deallocate ( work2 )
    deallocate ( work4, work5 )
    deallocate ( work6, work7 )
    deallocate ( work9 )

!   particle position

    xp(:ndim) = coefficients%r(16:15+ndim)

    if ( matrix .and. coorsys == 1 ) then

!     axisymmetric: only translation along z-direction

      elemmat(1,:) = - sum( elemmatsdiag(pos(:,1),:), dim=1 )

!     now the columns of C regarding the translational velocity U:
!     We use the partition of unity of the shape function: sum_k phi_k = 1, to
!     derive it from the matrix A

      elemmatdiag(:,1) = - sum( elemmat(:,pos(:,1)), dim=2 )

    else if ( matrix ) then

!     now the part of A regarding the translational velocity U:
!     We use the partition of unity of the shape function: sum_k phi_k = 1, to
!     derive it from the matrix S

      do i = 1, ndim
        elemmat(i,:) = - sum( elemmatsdiag(pos(:,i),:), dim=1 )
      end do

!     now the part of A regarding the rotational rate omega
!     We use the identity for isoparametric elements : sum_k x_k phi_k = x, to
!     derive it from the matrix S

      allocate ( work6(ndim,ndim,ndim*ndf+ndfp) )

      do i = 1, ndim
        do j = 1, ndim
          work6(i,j,:) = matmul( x(:,i) - xp(i), elemmatsdiag(pos(:,j),:) )
        end do
      end do

!     take outer product
      if ( ndim == 2 .and. coorsys == 0 ) then
        elemmat(3,:) = - ( work6(1,2,:) - work6(2,1,:) )
      else if ( ndim == 3 ) then
        elemmat(4,:) = - ( work6(2,3,:) - work6(3,2,:) )
        elemmat(5,:) = - ( work6(3,1,:) - work6(1,3,:) )
        elemmat(6,:) = - ( work6(1,2,:) - work6(2,1,:) )
      end if

      deallocate ( work6 )

!     now the columns of C regarding the translational velocity U:
!     We use the partition of unity of the shape function: sum_k phi_k = 1, to
!     derive it from the matrix A

      do j = 1, ndim
        elemmatdiag(:,j) = - sum( elemmat(:,pos(:,j)), dim=2 )
      end do

      allocate ( work6(ndim,ndim,size(elemmat,1)) )

!     now the columns of C regarding the rotational rate omega
!     We use the identity for isoparametric elements : sum_k x_k phi_k = x, to
!     derive it from the matrix A

      do i = 1, ndim
        do j = 1, ndim
          work6(i,j,:) = matmul( elemmat(:,pos(:,j)), x(:,i) - xp(i) )
        end do
      end do

!     take outer product
      if ( ndim == 2 .and. coorsys == 0 ) then
        elemmatdiag(:,3) = - ( work6(1,2,:) - work6(2,1,:) )
      else if ( ndim == 3 ) then
        elemmatdiag(:,4) = - ( work6(2,3,:) - work6(3,2,:) )
        elemmatdiag(:,5) = - ( work6(3,1,:) - work6(1,3,:) )
        elemmatdiag(:,6) = - ( work6(1,2,:) - work6(2,1,:) )
      end if

      deallocate ( work6 )

    end if

    call unset_stokes_eltree

  end subroutine stokes_embedded_particle_eltree


! Nitsche boundary element for rigid particle velocity using an eltree.
! The elemsub4 interface of build_system_constraint is used.
! A global constraint on an elementset must be used with Lagrange multipliers
! acting as unknown particle velocity and particle rotation rate.
! The diagonal block S is also filled (regular velocity stabilization term).
! The diagonal block C is also filled, so diagonal_block=.true. must be used.
! No additional unknowns.
! NOTE: arguments changed: elem,node -> elem1,elgrp

  subroutine stokes_Nitsche_particle_eltree ( mesh, problem, constr, elem1, &
    elgrp, matrix, vector, first, last, coefficients, oldvectors, &
    elemmatsdiag, elemmatsdiag2, elemmat, elemmat2, elemmatdiag, elemmatadd, &
    elemmatdiagadd, elemvec, elemvecadd )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem1, elgrp
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmatsdiag, elemmatsdiag2
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:,:) :: elemmatdiag, elemmatdiagadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


!   variables

    integer :: elementset1, elem
    integer :: i, j, N, M, ip
    real(dp) :: factor, he, xp(3)


!   find real element number

    elementset1 = problem%constraints(constr)%elementset1

    if ( elementset1 == 0 ) then
      write(*,'(/a/a/)') &
        'Error in stokes_Nitsche_particle_eltree: ', &
        ' constraint not defined on an elementset'
      stop
    end if

    elem = mesh%elementsets(elementset1)%elements(elgrp)%a(elem1)

!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( da(ninti) )
    allocate ( work2(ndf,ndf) )

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
      da(ip) = sqrt(dot_product(dan(ip,:),dan(ip,:)))
    end do

!   element length scaling

    he = 1 ! no elementscaling yet

!   Nitsche parameter

    factor = coefficients%r(15)

!   vector

    if ( vector ) then

      elemvec = 0

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

      elemmatsdiag = 0

      do i = 1, ndim
        elemmatsdiag( pos(:,i), pos(:,i) ) = work2
      end do

    end if

!   deallocate arrays

    deallocate ( da )
    deallocate ( work2 )

!   particle position

    xp(:ndim) = coefficients%r(16:15+ndim)

    if ( matrix .and. coorsys == 1 ) then

!     axisymmetric: only translation along z-direction

      elemmat(1,:) = - sum( elemmatsdiag(pos(:,1),:), dim=1 )

!     now the columns of C regarding the translational velocity U:
!     We use the partition of unity of the shape function: sum_k phi_k = 1, to
!     derive it from the matrix A

      elemmatdiag(:,1) = - sum( elemmat(:,pos(:,1)), dim=2 )

    else if ( matrix ) then

!     now the part of A regarding the translational velocity U:
!     We use the partition of unity of the shape function: sum_k phi_k = 1, to
!     derive it from the matrix S

      do i = 1, ndim
        elemmat(i,:) = - sum( elemmatsdiag(pos(:,i),:), dim=1 )
      end do

!     now the part of A regarding the rotational rate omega
!     We use the identity for isoparametric elements : sum_k x_k phi_k = x, to
!     derive it from the matrix S

      allocate ( work6(ndim,ndim,ndim*ndf+ndfp) )

      do i = 1, ndim
        do j = 1, ndim
          work6(i,j,:) = matmul( x(:,i) - xp(i), elemmatsdiag(pos(:,j),:) )
        end do
      end do

!     take outer product
      if ( ndim == 2 .and. coorsys == 0 ) then
        elemmat(3,:) = - ( work6(1,2,:) - work6(2,1,:) )
      else if ( ndim == 3 ) then
        elemmat(4,:) = - ( work6(2,3,:) - work6(3,2,:) )
        elemmat(5,:) = - ( work6(3,1,:) - work6(1,3,:) )
        elemmat(6,:) = - ( work6(1,2,:) - work6(2,1,:) )
      end if

      deallocate ( work6 )

!     now the columns of C regarding the translational velocity U:
!     We use the partition of unity of the shape function: sum_k phi_k = 1, to
!     derive it from the matrix A

      do j = 1, ndim
        elemmatdiag(:,j) = - sum( elemmat(:,pos(:,j)), dim=2 )
      end do

      allocate ( work6(ndim,ndim,size(elemmat,1)) )

!     now the columns of C regarding the rotational rate omega
!     We use the identity for isoparametric elements : sum_k x_k phi_k = x, to
!     derive it from the matrix A

      do i = 1, ndim
        do j = 1, ndim
          work6(i,j,:) = matmul( elemmat(:,pos(:,j)), x(:,i) - xp(i) )
        end do
      end do

!     take outer product
      if ( ndim == 2 .and. coorsys == 0 ) then
        elemmatdiag(:,3) = - ( work6(1,2,:) - work6(2,1,:) )
      else if ( ndim == 3 ) then
        elemmatdiag(:,4) = - ( work6(2,3,:) - work6(3,2,:) )
        elemmatdiag(:,5) = - ( work6(3,1,:) - work6(1,3,:) )
        elemmatdiag(:,6) = - ( work6(1,2,:) - work6(2,1,:) )
      end if

      deallocate ( work6 )

    end if

    call unset_stokes_eltree

  end subroutine stokes_Nitsche_particle_eltree

end module stokes_elements_embedded_particle_m

