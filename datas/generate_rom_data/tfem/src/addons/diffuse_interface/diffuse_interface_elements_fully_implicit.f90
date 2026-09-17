
! Copyright (C) 2015-2015 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the diffuse-interface equation and momentum
! balance solved in one system
!
!                  dcdt + u.nabla c - div( M nabla c )  = 0
!             alpha c - beta c^3 + kappa nabla^2 c + mu = 0
!  -\nabla . (2 \eta(c) D) + \nabla g - \nabla . \tau_c = 0
!                                            \nabla . u = 0
!
! \nabla . \tau_c can be rewritten to \mu\nabla c or -c\nabla \mu
! (with terms that can be absorbed in the pressure)
!
! The solution vector consists of (c,mu,u,g).
!

module diffuse_interface_elements_fully_implicit_m

  use tfem_elem_m
  use stokes_set_globals_m

  implicit none


contains


! element routine to build the implicit convection terms in the DI-equation

  subroutine di_convection_implicit_newton_elem ( mesh, problem, elgrp, elem, &
     matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, k, ip


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )


 !  allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

!   some tests

    if ( first ) then

      call check ( coefficients, 'di_convection_implicit_newton_elem', &
        ncoefi=200, ncoefr=150, indexarray=[151,152], minimum=[0,2], &
        maximum=[2,2]  )

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   get c variable (concentration) at previous iteration step

    call get_sysvector ( mesh, problem, oldvectors%s(3)%p, elgrp, elem, &
      citer, physq=[1], layer=layer )

    citer_g = matmul ( phi, citer )

!   grad c term

    do i = 1, ndim
      gradc(:,i) = matmul ( dphidx(:,:,i), citer )
    end do

    if ( matrix ) then

!     term u_i+1 . grad c_i
!     NOTE: the mesh velocity in the u_i+1 . grad c_i term cancels with the
!     mesh velocity in the u_i . grad c_i term

      do k = 1, ndim
        do j = 1, ndf
          do i = 1, ndf
            work6(i,j,k) = sum ( phi(:,i) * gradc(:,k) * phi(:,j) * detF * wg )
          end do
        end do
      end do

      do k = 1, ndim
        elemmat(1:ndf,(k-1)*ndf+1:k*ndf) = work6(:,:,k)
      end do

    end if

    if ( vector ) then

      elemvec = 0

!     get velocity at previous iteration step

      call get_sysvector ( mesh, problem, oldvectors%s(3)%p, &
        elgrp, elem, u, physq=[physqvel], layer=layer )

      tmp = reshape ( u, [ndf,ndim] )

      u_iter = matmul ( phi, tmp )

!     term u_i . grad c_i
!     NOTE: the mesh velocity in the u_i+1 . grad c_i term cancels with the
!     mesh velocity in the u_i . grad c_i term

      do ip = 1, ninti
        ungradc(ip) = dot_product ( u_iter(ip,:), gradc(ip,:) )
      end do

      do i = 1, ndf
        elemvec(i) = sum ( phi(:,i) * ungradc * detF * wg )
      end do

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

      allocate ( ungradc(ninti), u(ndim*ndf) )
      allocate ( citer(ndf), citer_g(ninti) )
      allocate ( gradc(ninti,ndim) )
      allocate ( u_iter(ninti,ndim), tmp(ndf,ndim) )
      allocate ( work6(ndf,ndf,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( ungradc, u )
      deallocate ( citer, citer_g )
      deallocate ( gradc )
      deallocate ( u_iter, tmp )
      deallocate ( work6 )

    end subroutine deallocate_arrays

  end subroutine di_convection_implicit_newton_elem


! element routine to build the implicit terms of mugradc in the momentum
! balance

  subroutine mugradc_implicit_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, k
    real(dp) :: rho


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


!   some tests

    if ( first ) then

      call check ( coefficients, 'mugradc_implicit_elem', ncoefi=200, &
        ncoefr=150 )

      call check ( coefficients, 'mugradc_implicit_elem', ncoefi=200, &
        ncoefr=150, indexarray=[151,152], minimum=[0,0], maximum=[2,2]  )

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )


!   get c variable (concentration) at previous iteration step

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(3)%p, elgrp, &
      elem, c_i, physq=[1], layer=layer )


!   grad c term

    do i = 1, ndim
      gradc(:,i) = matmul ( dphidx(:,:,i), c_i )
    end do


!   get mu variable (chemical potential) at previous iteration step

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(3)%p, elgrp, &
      elem, mu, physq=[2], layer=layer )

    mu_g = matmul ( phi, mu )


!   coefficients

    rho  = coefficients%r(106)


    if ( matrix ) then

      elemmat = 0

!     term rho * mu_i+1 * gradc_i
      do k = 1, ndim
        do j = 1, ndf
          do i = 1, ndf
            work6(i,j,k) = &
                  -rho * sum ( phi(:,i) * gradc(:,k) * phi(:,j) * detF * wg )
          end do
        end do
      end do

      do k = 1, ndim
        elemmat((k-1)*ndf+1:k*ndf,ndf+1:2*ndf) = work6(:,:,k)
      end do

      if ( coefficients%i(152) == 2 ) then ! Newton-Raphson iteration

!       term rho * mu_i * gradc_i+1
        do k = 1, ndim
          do j = 1, ndf
            do i = 1, ndf
              work6(i,j,k) = &
                   -rho * sum ( phi(:,i) * mu_g * dphidx(:,j,k) * detF * wg )
            end do
          end do
        end do

        do k = 1, ndim
          elemmat((k-1)*ndf+1:k*ndf,1:ndf) = work6(:,:,k)
        end do

      end if

    end if


    if ( vector ) then

      elemvec = 0

      if ( coefficients%i(152) == 2 ) then ! Newton-Raphson iteration

!       term -rho mu * gradc

        do j = 1, ndim
          do i = 1, ndf
            work2(i,j) = &
                     - rho * sum ( phi(:,i) * mu_g * gradc(:,j) * detF * wg )
          end do
        end do

        elemvec = reshape ( work2, [ndim*ndf] )

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

      allocate ( work6(ndf,ndf,ndim) )
      allocate ( c_i(ndf), gradc(ninti,ndim) )
      allocate ( mu(ndf), mu_g(ninti) )
      allocate ( work2(ndf,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work6 )
      deallocate ( c_i, gradc )
      deallocate ( mu, mu_g )
      deallocate ( work2 )

    end subroutine deallocate_arrays

  end subroutine mugradc_implicit_elem


! element routine to build the implicit terms of -cgradmu in the momentum
! balance

  subroutine mincgradmu_implicit_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, k
    real(dp) :: rho


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

!   some tests

    if ( first ) then

      call check ( coefficients, 'mincgradmu_implicit_elem', ncoefi=200, &
        ncoefr=150 )

      call check ( coefficients, 'mincgradmu_implicit_elem', ncoefi=200, &
        ncoefr=150, indexarray=[151,152], minimum=[0,0], maximum=[2,2]  )

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )


!   get c variable (concentration) at previous iteration step

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(3)%p, elgrp, &
      elem, c_n, physq=[1], layer=layer )

    c_ng = matmul ( phi, c_n )


!   get mu variable (chemical potential) at previous iteration step

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(3)%p, elgrp, &
      elem, mu, physq=[2], layer=layer )


!   grad mu term

    do i = 1, ndim
      gradmu(:,i) = matmul ( dphidx(:,:,i), mu )
    end do


!   coefficients

    rho  = coefficients%r(106)

    if ( matrix ) then

      elemmat = 0

!     term rho * c_i * gradmu_i+1
      do k = 1, ndim
        do j = 1, ndf
          do i = 1, ndf
            work6(i,j,k) = &
                    rho * sum ( phi(:,i) * c_ng * dphidx(:,j,k) * detF * wg )
          end do
        end do
      end do

      do k = 1, ndim
        elemmat((k-1)*ndf+1:k*ndf,ndf+1:2*ndf) = work6(:,:,k)
      end do

      if ( coefficients%i(152) == 2 ) then ! Newton-Raphson iteration

!       term rho * c_i+1 * gradmu_i
        do k = 1, ndim
          do j = 1, ndf
            do i = 1, ndf
              work6(i,j,k) = &
                   rho * sum ( phi(:,i) * gradmu(:,k) * phi(:,j) * detF * wg )
            end do
          end do
        end do

        do k = 1, ndim
          elemmat((k-1)*ndf+1:k*ndf,1:ndf) = work6(:,:,k)
        end do

      end if

    end if


    if ( vector ) then

      elemvec = 0

      if ( coefficients%i(152) == 2 ) then ! Newton-Raphson iteration

!       term rho c * grad mu

        do j = 1, ndim
          do i = 1, ndf
            work2(i,j) = rho * sum ( phi(:,i) * c_ng * gradmu(:,j) * detF * wg )
          end do
        end do

        elemvec = reshape ( work2, [ndim*ndf] )

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

      allocate ( work6(ndf,ndf,ndim) )
      allocate ( c_n(ndf), gradmu(ninti,ndim) )
      allocate ( mu(ndf), c_ng(ninti) )
      allocate ( work2(ndf,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work6 )
      deallocate ( c_n, gradmu )
      deallocate ( mu, c_ng )
      deallocate ( work2 )

    end subroutine deallocate_arrays

  end subroutine mincgradmu_implicit_elem


! Element routine for the implicit divergence of the capillary stress in the
! momentum balance

  subroutine divstressdi_implicit_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, N, L, ip
    real(dp) :: kappa, rho
    real(dp) :: fac


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )


!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

!   some tests

    if ( first ) then

      call check ( coefficients, 'divstressdi_implicit_elem', ncoefi=200, &
        ncoefr=150 )

      call check ( coefficients, 'divstressdi_implicit_elem', ncoefi=200, &
        ncoefr=150, indexarray=[151,152], minimum=[0,0], maximum=[2,2]  )

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   coefficients

    kappa  = coefficients%r(105)
    rho    = coefficients%r(106)

!   get c variable (concentration)

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(3)%p, elgrp, &
      elem, c_n, physq=[1], layer=layer )

!   grad c term from previous iteration step

    do i = 1, ndim
      gradc(:,i) = matmul ( dphidx(:,:,i), c_n )
    end do

!   the -gradcgradc tensor

    do i = 1, ndim
      do j = i, ndim
        gradcgradc(:,i,j) = - gradc(:,i) * gradc(:,j)
        if ( i /= j ) gradcgradc(:,j,i) = gradcgradc(:,i,j)
      end do
    end do

!   the isotropic part |gradc|^2

    gradc_sq = sum ( gradc**2, 2 )

!   combine the gradcgradc and |gradc|^2 I tensors to form the capillary
!   stress tensor

    tau_c = rho * kappa * gradcgradc
    do i = 1, ndim
      tau_c(:,i,i) = tau_c(:,i,i) + rho * kappa * gradc_sq
    end do

!   build equations

    if ( matrix ) then

!     build the gradcdotgradphi and gradcgradphi structures

      do ip = 1, ninti
        do N = 1, ndf
          gradcdotgradphi(ip,N) = dot_product ( gradc(ip,:), dphidx(ip,N,:) )
        end do
      end do

      do N = 1, ndf
        do i = 1, ndim
          do j = 1, ndim
            gradcgradphi(:,N,i,j) = gradc(:,i) * dphidx(:,N,j)
          end do
        end do
      end do

!     build the isotropic part |gradc|^2 I equation

      select case ( coefficients%i(152) )
      case(0,1) ! Picard iteration
        fac = 1._dp
      case(2) ! Newton-Raphson iteration
        fac = 2._dp
      case default
        call errormsg_case_default ( 'divstressdi_implicit_elem', &
          'coefficients%i(152)', int_value=coefficients%i(152) )
      end select

      do i = 1, ndim
        do N = 1, ndf
          do L = 1, ndf
            iso_mat(ndf*(i-1)+N,L) = fac * rho * kappa * &
              sum ( dphidx(:,N,i) * gradcdotgradphi(:,L) * detF * wg )
          end do
        end do
      end do

      if ( coorsys == 1 ) then
        do N = 1,ndf
          do L = 1,ndf
!           only apply to row ndf+1:2*ndf: the r-velocities
            iso_mat(ndf+N,L) = iso_mat(ndf+N,L) + &
              fac * rho * kappa * sum ( phi(:,N) * gradcdotgradphi(:,L) * &
                                            detF * wg / xg(:,2) )
          end do
        end do
      end if

!     build the gradcgradc tensor equation

      do i = 1, ndim
        do N = 1, ndf
          do L = 1, ndf

            work1 = sum ( dphidx(:,N,:) * gradcgradphi(:,L,:,i), 2 )

            if ( coefficients%i(152) == 2 ) then ! Newton-Raphson iteration
              work1 = work1 + sum ( dphidx(:,N,:) * gradcgradphi(:,L,i,:), 2 )
            end if

            gradcgradc_mat(ndf*(i-1)+N,L) = &
                           - rho * kappa * sum ( work1 * detF * wg )

          end do
        end do
      end do

      elemmat = iso_mat + gradcgradc_mat

    end if

    if ( vector ) then

      elemvec = 0

      if ( coefficients%i(152) == 2 ) then ! Newton-Raphson iteration

!       + (nabla v)^T:stress_di

        do ip = 1, ninti
          work6(ip,:,:) = matmul ( dphidx(ip,:,:), tau_c(ip,:,:) )
        end do
        if ( coorsys == 1 ) then
          do N = 1, ndf
            work6(:,N,2) = work6(:,N,2) + &
                                rho * kappa * gradc_sq * phi(:,N) / xg(:,2)
          end do
        end if
        do j = 1, ndim
          work2(:,j) = matmul ( detF * wg, work6(:,:,j) )
        end do
        elemvec = reshape ( work2, [ ndf*ndim ] )

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

      allocate ( work6(ninti,ndf,ndim), work2(ndf,ndim) )
      allocate ( c_n(ndf), gradc(ninti,ndim) )
      allocate ( gradcdotgradphi(ninti,ndf) )
      allocate ( gradcgradphi(ninti,ndf,ndim,ndim) )
      allocate ( gradcgradc(ninti,ndim,ndim) )
      allocate ( tau_c(ninti,ndim,ndim) )
      allocate ( gradc_sq(ninti) )
      allocate ( iso_mat(ndim*ndf,ndf) )
      allocate ( gradcgradc_mat(ndim*ndf,ndf) )
      allocate ( work1(ninti) )


    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work6, work2 )
      deallocate ( c_n, gradc )
      deallocate ( gradcdotgradphi )
      deallocate ( gradcgradphi )
      deallocate ( gradcgradc )
      deallocate ( tau_c )
      deallocate ( gradc_sq )
      deallocate ( iso_mat )
      deallocate ( gradcgradc_mat )
      deallocate ( work1 )

    end subroutine deallocate_arrays

  end subroutine divstressdi_implicit_elem

end module diffuse_interface_elements_fully_implicit_m
