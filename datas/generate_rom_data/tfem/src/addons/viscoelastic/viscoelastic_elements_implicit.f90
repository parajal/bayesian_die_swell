
! Copyright (C) 2023-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for viscoelastic flow employing implicit time integration
!

module viscoelastic_elements_implicit_m

  use tfem_elem_m
  use stokes_set_globals_m
  use devss_set_globals_m
  use viscoelastic_models_m
  use viscoelastic_elements_generic_m

  implicit none

contains


! Internal element routine for the Jacobian wrt the c/b tensor and
! rhs of the - (div tau) term in the momentum balance

  subroutine implicit_divtau ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m
    use tensor_m
    use viscoelastic_elements_generic_m, only: stress_tensor_viscoelastic, &
      dstress_tensor_viscoelastic

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: j, k, l, N, m, ip, nd, ist, isv_c


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, gammap=.true. )

!   coefficients

    if ( first ) then
      cstorage = get_coefficient ( coefficients, index=84, default=3 )
    end if

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

!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure ( mesh, problem, elgrp, elem, coefficients, &
        oldvectors%s(1)%p, vemopt )

    end if


!   get conformation tensors

    isv_c = get_coefficient ( coefficients, index=87, default=1 )

    do m = mode1, mode2

!     single mode conformation
      call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
        cmode=c, mode=m, isv=isv_c )

      cg(:,:,m) = matmul ( theta, c )

    end do

!   build equations

    if ( matrix ) then

      nd = ndim*(ndim+1)/2 ! number of stress components for matrix ndim x ndim

      ist = 0 ! starting column pointer for elemmat

      do m = mode1, mode2

!       modal stress derivative use of global variables; result in dtauvec

        call dstress_tensor_viscoelastic ( cg, mode=m, vemopt=vemopt )

!       (nabla v)^T:dtau_mode/dc * delta c

        do j = 1, ncomp

          do ip = 1, ninti
            work6(ip,:,:) = &
               matmul ( dphidx(ip,:,:), &
                    vector_to_tensor2_symmetric ( ndim, dtauvec(ip,1:nd,j) ) )
          end do

          if ( coorsys == 1 .and. vel3D == 0 ) then
            do N = 1, ndf
!             v_r/r * tau_theta,theta
              work6(:,N,2) = work6(:,N,2) + dtauvec(:,4,j) * phi(:,N) / xg(:,2)
            end do
          else if ( coorsys == 1 .and. vel3D == 1 ) then
            do N = 1, ndf
!             v_r/r * tau_theta,theta
              work6(:,N,2) = work6(:,N,2) + dtauvec(:,6,j) * phi(:,N) / xg(:,2)
!             - v_theta/r * tau_theta,r
              work6(:,N,3) = work6(:,N,3) - dtauvec(:,5,j) * phi(:,N) / xg(:,2)
            end do
          end if

          do l = 1, ndfc
            work1 = detF * wg * theta(:,l)
            do k = 1, ncompu
              work7(:,k,l) = matmul ( work1 , work6(:,:,k) )
            end do
          end do

          elemmat(:,ist+1:ist+ndfc) = reshape ( work7, [ ndf*ncompu, ndfc ] )

          ist = ist + ndfc

        end do

      end do

    end if

    if ( vector ) then

!     stress tensor: use of global variables; result in tauten

      call stress_tensor_viscoelastic ( cg, vemopt=vemopt )

!     - (nabla v)^T:tau

      do ip = 1, ninti
        work6(ip,:,:) = matmul ( dphidx(ip,:,:), tauten(ip,1:ndim,:) )
      end do
      if ( coorsys == 1 .and. vel3D == 0 ) then
        do N = 1, ndf
!         v_r/r * tau_theta,theta
          work6(:,N,2) = work6(:,N,2) + tauvec(:,4) * phi(:,N) / xg(:,2)
        end do
      else if ( coorsys == 1 .and. vel3D == 1 ) then
        do N = 1, ndf
!         v_r/r * tau_theta,theta
          work6(:,N,2) = work6(:,N,2) + tauvec(:,6) * phi(:,N) / xg(:,2)
!         - v_theta/r * tau_theta,r
          work6(:,N,3) = work6(:,N,3) - tauvec(:,5) * phi(:,N) / xg(:,2)
        end do
      end if
      do j = 1, ncompu
        work2(:,j) = - matmul ( detF * wg, work6(:,:,j) )
      end do
      elemvec = reshape ( work2, [ ndf*ncompu ] )

    end if

!   deallocate more memory

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )

  contains

    subroutine allocate_arrays

      allocate ( work1(ninti) )
      allocate ( work6(ninti,ndf,ncompu), work2(ndf,ncompu) )
      allocate ( work7(ndf,ncompu,ndfc) )
      allocate ( c(ndfc,ncomp), cg(ninti,ncomp,nmodes) )
      allocate ( tauvec(ninti,ncompt), tauten(ninti,ncompu,ncompu) )
      allocate ( dtauvec(ninti,ncompt,ncomp) )

      if ( vemcompressible ) then
        allocate ( pr(ndfp), press(ninti) )
        call create_vemopt ( vemopt, dep_J=.true., np=ninti )
      else
        call create_vemopt ( vemopt )
      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work1 )
      deallocate ( work6, work2 )
      deallocate ( work7 )
      deallocate ( c, cg )
      deallocate ( tauvec, tauten )
      deallocate ( dtauvec )

      if ( vemcompressible ) then
        deallocate ( pr, press )
      end if

      call delete ( vemopt )

    end subroutine deallocate_arrays

  end subroutine implicit_divtau


! Internal element routine for the Jacobian wrt the relative change in volume J
! of the - (div tau) term in the momentum balance. No rhs.

  subroutine implicit_divtau_J ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m
    use tensor_m
    use viscoelastic_elements_generic_m, only: dstressdJ_tensor_viscoelastic

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, N, m, ip, isv_c
    real(dp) :: Kmod


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, gammap=.true. )

!   coefficients

    if ( first ) then

!     check if model is dependent on J

      if ( .not. vemcompressible ) then
        write(*,'(/6(a/))') 'Error in implicit_divtau_J:', &
          ' Element routine only applicable for models that depend on J.'
        stop
      end if

      cstorage = get_coefficient ( coefficients, index=84, default=3 )

    end if

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


!   get relative change in volume J from pressure in system vector

    Kmod = coefficients%r(29)

    call get_J_from_pressure ( mesh, problem, elgrp, elem, coefficients, &
      oldvectors%s(1)%p, vemopt )


!   get conformation tensors

    isv_c = get_coefficient ( coefficients, index=87, default=1 )

    do m = mode1, mode2

!     single mode conformation
      call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
        cmode=c, mode=m, isv=isv_c )

      cg(:,:,m) = matmul ( theta, c )

    end do

!   build equations

    if ( matrix ) then

!     derivative of stress tensor wrt to J: use of global variables;
!     result in dtautendJ

      call dstressdJ_tensor_viscoelastic ( cg, vemopt=vemopt )

!     - (nabla v)^T:dtaudJ * J/K * q

      do ip = 1, ninti
        work6(ip,:,:) = matmul ( dphidx(ip,:,:), dtautendJ(ip,1:ndim,:) )
      end do
      if ( coorsys == 1 .and. vel3D == 0 ) then
        do N = 1, ndf
!         v_r/r * dtaudJ_theta,theta
          work6(:,N,2) = work6(:,N,2) + dtauvecdJ(:,4) * phi(:,N) / xg(:,2)
        end do
      else if ( coorsys == 1 .and. vel3D == 1 ) then
        do N = 1, ndf
!         v_r/r * tau_theta,theta
          work6(:,N,2) = work6(:,N,2) + dtauvecdJ(:,6) * phi(:,N) / xg(:,2)
!         - v_theta/r * tau_theta,r
          work6(:,N,3) = work6(:,N,3) - dtauvecdJ(:,5) * phi(:,N) / xg(:,2)
        end do
      end if
      work2 = reshape ( work6, [ ninti, ndf*ncompu ] )
      do ip = 1, ninti
        work2(ip,:) = work2(ip,:) * vemopt%J(ip) / Kmod
      end do
      do j = 1, ndfp
        do i = 1, ndf*ncompu
          elemmat(i,j) = - sum ( detF * wg * work2(:,i) * psi(:,j) )
        end do
      end do

    end if

    if ( vector ) then
      elemvec = 0
    end if

!   deallocate more memory

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )

  contains

    subroutine allocate_arrays

      allocate ( c(ndfc,ncomp), cg(ninti,ncomp,nmodes) )
      allocate ( pr(ndfp), press(ninti) )
      allocate ( work6(ninti,ndf,ncompu), work2(ninti,ndf*ncompu) )
      allocate ( dtauvecdJ(ninti,ncompt), dtautendJ(ninti,ncompu,ncompu) )

      call create_vemopt ( vemopt, dep_J=.true., np=ninti )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( pr, press )
      deallocate ( c, cg )
      deallocate ( work6, work2 )
      deallocate ( dtauvecdJ, dtautendJ )

      call delete ( vemopt )

    end subroutine deallocate_arrays

  end subroutine implicit_divtau_J


! Internal element routine for the constitutive equation (multi-mode)
! Implicit time-integration with Newton-Raphson iteration.
! This is the conformation part (in Jacobian) and the rhs.

  subroutine implicit_ce_supg_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m
    use limits_m, only: BLOCK_IMPROPER_ALE, WARN_ON_INCOMPLETE_JACOBIAN
    use viscoelastic_elements_SUPG_m, only: add_ugrad_terms
    use supg_utils_m
    use viscoelastic_elements_generic_m, only: rhs_viscoelastic

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, k, l, ip, timeint, m, n, m1, n1, iv_supg
    integer, dimension(2) :: j14, j5, j3, jz, j12
    real(dp) :: beta, fac, deltat, esize, thetapar1, thetapar2
    real(dp) :: Kmod


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_devssg_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1, singlemode=.true., &
      gammap=.true. )

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'implicit_ce_supg_elem', &
        indexarray=[22,48,50,57,89,90,91], minimum=[8,0,0,0,0,0,0], &
        maximum=[11,1,1,0,1,1,1] )

      if ( coefficients%i(48) == 1 .and. any(coefficients%i(22) == [9,11]) &
           .and. BLOCK_IMPROPER_ALE ) then
        write(*,'(/6(a/))') 'Error in implicit_ce_supg_elem:', &
          ' Combining ALE with Crank-Nicolson or theta method not allowed, ', &
          ' since the change in the nabla operator due to the change ', &
          ' in mesh at different times within a single step ', &
          ' has not been taken into account yet.', &
          ' Set BLOCK_IMPROPER_ALE=.false. in limits_m to bypass this block.'
        stop
      end if

      if ( coefficients%i(89) == 1 ) then
        if ( all(coefficients%i(22)/=[8,10]) ) then
          write(*,'(/3(a/))') 'Error in implicit_ce_supg_elem:', &
            ' Excluding the time-derivate is only allowed if the ', &
            ' specified time integration scheme is BDF1 or BDF2.'
          stop
        end if
      end if

      if ( coefficients%i(90) == 1 ) then

        if ( coefficients%i(71) /= 1 ) then
          write(*,'(/3(a/))') 'Error in implicit_ce_supg_elem:', &
            ' Rotation reinitialization is only applicable for the ', &
            ' contravariant deformation tensor (CDT) formulation. '
          stop
        end if
        if ( any(coefficients%i(22) == [9,11]) ) then
          write(*,'(/a/a,i0,a/)') 'Error in implicit_ce_supg_elem:', &
            ' Rotation reinitialization for time integration method = ', &
              coefficients%i(22), ' not available.'
          stop
        end if

      end if

      if ( varpar ) then

!       check coefficients for variable coefficients
        call check ( coefficients, 'implicit_ce_supg_elem', &
          indexarray=[69,70,81], minimum=[0,0,0], maximum=[3,3,3] )

        if ( any(coefficients%i(22) == [9,11]) ) then
          write(*,'(/4(a/))') 'Error in implicit_ce_supg_elem:', &
            ' Combining variable coefficients with Crank-Nicolson ', &
            ' or theta-method not allowed, since the variable coefficients ', &
            ' are only available for a single time step.'
          stop
        end if

      end if

      if ( coefficients%i(50) == 1 .and. coefficients%i(91) == 1 ) then
        write(*,'(/3(a/))') 'Error in implicit_ce_supg_elem:', &
          ' Jacobian for SUPG test function not available if the ', &
          ' tau parameter is adjusted based on the Courant number.'
        stop
      end if

      if ( vemcompressible .and. any(coefficients%i(22) == [9,11]) &
            .and. WARN_ON_INCOMPLETE_JACOBIAN ) then

        write(*,'(/6(a/))') 'Warning in implicit_ce_supg_elem:', &
          ' The Jacobian for viscoelastic models that depend on J with', &
          ' Crank-Nicolson or theta-method is incomplete and the ', &
          ' Newton-Raphson iteration might show linear convergence.', &
          ' Set WARN_ON_INCOMPLETE_JACOBIAN=.false. in limits_m ', &
          ' to suppress this warning.'

      end if

      cstorage = get_coefficient ( coefficients, index=84, default=2 )

      if ( vemgammap ) then
        gpstorage = get_coefficient ( coefficients, index=100, default=2 )
      end if

    end if

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of element

    timeint = coefficients%i(22)

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coefficients%i(61) == 1 ) then
      call shape_derivative ( dphi, Finv, dphidx )
    end if

    call shape_derivative ( dtheta, Finv, dthetadx )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

!   get velocity vector unp1

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel], posu=posu )

    work4 = reshape ( u, [ndf,ncompu] )

    uvecnp1 = matmul ( phi, work4 )

!   get velocity vector u_supg

    iv_supg = get_coefficient ( coefficients, index=86, default=1 )

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(iv_supg)%p, &
      elgrp, elem, u, physq=[physqvel], posu=posu )

    work4 = reshape ( u, [ndf,ncompu] )

    uvec_supg = matmul ( phi, work4 )


!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      Kmod = coefficients%r(29)

      call get_J_from_pressure ( mesh, oldvectors%p(1)%p, elgrp, elem, &
        coefficients, oldvectors%s(1)%p, vemoptiter )

    end if


!   get velocity gradient vector at tn+1

    if ( coefficients%i(61) == 0 ) then

!     Projected gradient

      select case ( coefficients%i(60) )
      case(0)
!       DEVSS: in solution vector
        call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
          elgrp, elem, g, physq=[physqgrad], posu=posg )
      case(1)
!       separate vector
        call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(3)%p, &
          elgrp, elem, g, posu=posg )
      case default
        call errormsg_case_default ( 'implicit_ce_supg_elem', &
          'coefficients%i(60)', int_value=coefficients%i(60) )
      end select

      work2 = reshape ( g, [ndfg,ncompg] )

      work10 = matmul ( zeta, work2 )

      if ( vel3D == 1 ) then

        gradu(:,1,1) = work10(:,1)
        gradu(:,1,2) = work10(:,2)
        gradu(:,2,1) = work10(:,3)
        gradu(:,2,2) = work10(:,4)
        gradu(:,3,1) = work10(:,5)
        gradu(:,3,2) = work10(:,6)
        gradu(:,:,3) = 0

        do ip = 1, ninti
          gvecnp1(ip,:) = &
                       reshape ( transpose ( gradu(ip,:,:) ), [ncompu*ncompu] )
        end do

      else

        gvecnp1(:,1:ncompg) = work10

      end if

    else if ( coefficients%i(61) == 1 ) then

!     direct velocity gradient

!     get velocity vector

      call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
        elgrp, elem, u, physq=[physqvel], layer=layer )

      work4 = reshape ( u, [ndf,ncompu] )

      do j = 1, ndim
        gradu(:,:,j) = matmul ( dphidx(:,:,j), work4 )
      end do

      if ( vel3D == 1 ) then
        gradu(:,:,3) = 0
        do ip = 1, ninti
          gvecnp1(ip,:) = &
                       reshape ( transpose ( gradu(ip,:,:) ), [ncompu*ncompu] )
        end do
      else
        do ip = 1, ninti
          gvecnp1(ip,1:ndim**2) = &
                       reshape ( transpose ( gradu(ip,:,:) ), [ndim*ndim] )
        end do
      end if

    end if

    if ( coorsys == 1 .and. vel3D == 0 ) then
      gvecnp1(:,5) = uvecnp1(:,2) / xg(:,2)   ! u_r / r
    else if ( coorsys == 1 .and. vel3D == 1 ) then
      gvecnp1(:,6) = - uvecnp1(:,3) / xg(:,2) ! - u_theta / r
      gvecnp1(:,9) = uvecnp1(:,2) / xg(:,2)   ! u_r / r
    end if


!   set time step

    deltat = coefficients%r(8)


!   ALE, temporary ALE or Eulerian frame

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity at tn+1

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, um )

      tmp = reshape ( um, [ndf,ndim] )

      uvecmeshnp1 = matmul ( phi, tmp )

!     (unp1-ugridnp1).grad operator

      do ip = 1, ninti
        unp1gradtheta(ip,:) = &
            matmul ( dthetadx(ip,:,:), uvecnp1(ip,1:ndim) - uvecmeshnp1(ip,:) )
      end do

!     (u_supg-ugridnp1).grad operator

      do ip = 1, ninti
        ugradtheta_supg(ip,:) = &
         matmul ( dthetadx(ip,:,:), uvec_supg(ip,1:ndim) - uvecmeshnp1(ip,:) )
      end do

    else

!     Eulerian frame

      uvecmeshnp1 = 0

!     unp1.grad operator

      do ip = 1, ninti
        unp1gradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uvecnp1(ip,1:ndim) )
      end do

!     (u_supg-ugridnp1).grad operator

      do ip = 1, ninti
        ugradtheta_supg(ip,:) = &
                        matmul ( dthetadx(ip,:,:), uvec_supg(ip,1:ndim) )
      end do

    end if


!   get additional vectors at tn (CN and theta-method only)

    if ( any ( timeint == [9,11] ) ) then

!     get velocity vector un

      u = oldvectors%s(2)%p%u(posu)

      work4 = reshape ( u, [ndf,ncompu] )

      uvecn = matmul ( phi, work4 )

!     get relative change in volume J from pressure in system vector

      if ( vemcompressible ) then

        call get_J_from_pressure ( mesh, oldvectors%p(1)%p, elgrp, elem, &
          coefficients, oldvectors%s(2)%p, vemoptn )

      end if

!     get velocity gradient vector at tn

      if ( coefficients%i(61) == 0 ) then

!       Projected gradient

        select case ( coefficients%i(60) )
        case(0)
!         DEVSS: in solution vector
          g = oldvectors%s(2)%p%u(posg)
        case(1)
!         separate vector
          g = oldvectors%v(4)%p%u(posg)
        case default
          call errormsg_case_default ( 'implicit_ce_supg_elem', &
            'coefficients%i(60)', int_value=coefficients%i(60) )
        end select

        work2 = reshape ( g, [ndfg,ncompg] )

        work10 = matmul ( zeta, work2 )

        if ( vel3D == 1 ) then

          gradu(:,1,1) = work10(:,1)
          gradu(:,1,2) = work10(:,2)
          gradu(:,2,1) = work10(:,3)
          gradu(:,2,2) = work10(:,4)
          gradu(:,3,1) = work10(:,5)
          gradu(:,3,2) = work10(:,6)
          gradu(:,:,3) = 0

          do ip = 1, ninti
            gvecn(ip,:) = &
                      reshape ( transpose ( gradu(ip,:,:) ), [ncompu*ncompu] )
          end do

        else

          gvecn(:,1:ncompg) = work10

        end if

      else if ( coefficients%i(61) == 1 ) then

!       direct velocity gradient

!       get velocity vector

        u = oldvectors%s(2)%p%u(posu)

        work4 = reshape ( u, [ndf,ncompu] )

        do j = 1, ndim
          gradu(:,:,j) = matmul ( dphidx(:,:,j), work4 )
        end do

        if ( vel3D == 1 ) then
          gradu(:,:,3) = 0
          do ip = 1, ninti
            gvecn(ip,:) = &
                       reshape ( transpose ( gradu(ip,:,:) ), [ncompu*ncompu] )
          end do
        else
          do ip = 1, ninti
            gvecn(ip,1:ndim**2) = &
                        reshape ( transpose ( gradu(ip,:,:) ), [ndim*ndim] )
          end do
        end if

      end if

      if ( coorsys == 1 .and. vel3D == 0 ) then
        gvecn(:,5) = uvecn(:,2) / xg(:,2)   ! u_r / r
      else if ( coorsys == 1 .and. vel3D == 1 ) then
        gvecn(:,6) = - uvecn(:,3) / xg(:,2) ! - u_theta / r
        gvecn(:,9) = uvecn(:,2) / xg(:,2)   ! u_r / r
      end if

!     ALE or Eulerian frame

      if ( coefficients%i(48) == 1 ) then

!       ALE

!       get mesh velocity at tn

        call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(2)%p, elgrp, &
          elem, um )

        tmp = reshape ( um, [ndf,ndim] )

        uvecmeshn = matmul ( phi, tmp )

!       (un-ugrid).grad operator

        do ip = 1, ninti
          ungradtheta(ip,:) = &
                matmul ( dthetadx(ip,:,:), uvecn(ip,1:ndim) - uvecmeshn(ip,:) )
        end do

      else

!       Eulerian frame

!       un.grad operator

        do ip = 1, ninti
          ungradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uvecn(ip,1:ndim) )
        end do

      end if

    end if


!   get conformation tensor at previous time steps (cn, cnm1) and citer

    call get_c_standard


!   get gpiter and for CN at previous time step gpn

    if ( vemgammap ) call get_gp_standard


!   upwinding parameter and other factors

    beta = coefficients%r(9)

    if ( htype == 3 ) then
      esize = sum ( detF * wg )
    else
      esize = 1
    end if

!   thetapar

    if ( timeint == 9 ) then
      thetapar1 = 0.5_dp; thetapar2 = 0.5_dp
    else if ( timeint == 11 ) then
      thetapar1 = coefficients%r(28); thetapar2=1-thetapar1
    end if

!   compute h/U

    if ( coefficients%i(91) == 0 ) then
!     default h/U
      call supg_hU ( ndim, globalshape, hoverU, &
        htype=htype, Uscaling=Uscaling, Uglobal=coefficients%r(10), &
        uvec=uvec_supg(:,1:ndim)-uvecmeshnp1, Finv=Finv, x=x, esize=esize, &
        checkzero=coefficients%i(78)==1 )
    else if ( coefficients%i(91) == 1 ) then
!     include Jacobian of SUPG test function
      call supg_hU ( ndim, globalshape, hoverU, &
        htype=htype, Uscaling=Uscaling, Uglobal=coefficients%r(10), &
        uvec=uvec_supg(:,1:ndim)-uvecmeshnp1, Finv=Finv, x=x, esize=esize, &
        checkzero=coefficients%i(78)==1, dhoverU=dhoverU )
    end if

    if ( coefficients%i(50) == 1 ) then
!     Courant number dependent tau
      where ( deltat < hoverU )
        hoverU = deltat
      end where
    end if

    tau = beta * hoverU / 2  ! upwinding parameter

    if ( coefficients%i(89) == 1 ) then
      fac = 0  ! exclude time-derivative
    else
      fac = 1._dp / deltat
    end if

    if ( varpar ) then

!     determine the variable material parameters in mvemodel

      call evaluate_mve ( mesh, problem, elem, elgrp, coefficients, &
        oldvectors, lambda=.true., nonlin=.true., alam=.true. )

    end if

!   build equations

    if ( matrix ) then

      select case ( timeint )

      case(8)

!       implicit Euler (BDF1)

!       L^n+1.citer + citer.L^n+1T + f(citer)
        call rhs_viscoelastic ( gvecnp1, citerg, fiterg, vemoptiter )

        if ( coorsys == 1 .and. vel3D == 1 ) then
!         axisymmetric with 3D velocities: add terms due to u.nabla c
          call add_ugrad_terms ( xg, uvecnp1, citerg, fiterg, &
            dfadd_mm=vemoptiter%drhs_mm )
        end if

!       Jacobian of -rhs

        do n = mode1, mode2
          n1 = n - mode1 + 1
          do l = 1, ncomp
            do k = 1, ndfc
              do m = mode1, mode2
                m1 = m - mode1 + 1
                do j = 1, ncomp
                  do i = 1, ndfc
                    work14(i,j,m1,k,l,n1) = &
                      - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                        vemoptiter%drhs_mm(:,j,m,l,n) * theta(:,k) * detF * wg )
                  end do
                end do
              end do
            end do
          end do
        end do

!       Jacobian of material time derivative (identical for all components)

        do j = 1, ndfc
          do i = 1, ndfc
            work11(i,j) = sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
             ( fac * theta(:,j) + unp1gradtheta(:,j) ) &
                                   * detF * wg )
          end do
        end do

        if ( coefficients%i(91) == 1 ) then

!         Jacobian of SUPG test function

          do l = 1, ndim
            do k = 1, ndf
              do m = mode1, mode2
                m1 = m - mode1 + 1
                do j = 1, ncomp
                  do i = 1, ndfc
                    work3(i,j,m1,k,l) = &
                     - sum ( ( beta/2 * dhoverU(:,l) * ugradtheta_supg(:,i) &
                                 + tau * dthetadx(:,i,l) ) * phi(:,k) * &
                            ( - fac * citerg(:,j,m) + fac * cng(:,j,m) &
                                - unp1gradciter(:,j,m) + &
                                       fiterg(:,j,m) ) * detF * wg )
                  end do
                end do
              end do
            end do
          end do

        end if

        if ( vemcompressible ) then

!         include Jacobian of J

          do k = 1, ndfp
            do m = mode1, mode2
              m1 = m - mode1 + 1
              do j = 1, ncomp
                do i = 1, ndfc
                  work5(i,j,m1,k) = &
                    sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                               vemoptiter%drhsdJ(:,j,m) * &
                                  vemoptiter%J * psi(:,k) * detF * wg ) / Kmod
                end do
              end do
            end do
          end do

        end if

        if ( vemgammap ) then

!         include Jacobian of gammap

          do k = 1, ndfc
            do m = mode1, mode2
              m1 = m - mode1 + 1
              do j = 1, ncomp
                do i = 1, ndfc
                  work12(i,j,m1,k) = &
                   - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                               vemoptiter%drhsdgammap(:,j,m) * &
                                  theta(:,k) * detF * wg )
                end do
              end do
            end do
          end do

        end if

     case(9,11)

!       Crank-Nicolson (trapezoidal) or theta method

!       L^n+1.citer + citer.L^n+1T + f(citer)
        call rhs_viscoelastic ( gvecnp1, citerg, fiterg, vemoptiter )

        if ( coorsys == 1 .and. vel3D == 1 ) then
!         axisymmetric with 3D velocities: add terms due to u.nabla c
          call add_ugrad_terms ( xg, uvecnp1, citerg, fiterg, &
            dfadd_mm=vemoptiter%drhs_mm )
        end if

!       Jacobian of -rhs

        do n = mode1, mode2
          n1 = n - mode1 + 1
          do l = 1, ncomp
            do k = 1, ndfc
              do m = mode1, mode2
                m1 = m - mode1 + 1
                do j = 1, ncomp
                  do i = 1, ndfc
                    work14(i,j,m1,k,l,n1) = &
                      - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                                  thetapar1 * vemoptiter%drhs_mm(:,j,m,l,n) * &
                                    theta(:,k) * detF * wg )
                  end do
                end do
              end do
            end do
          end do
        end do

!       Jacobian of material time derivative (identical for all components)

        do j = 1, ndfc
          do i = 1, ndfc
            work11(i,j) = sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                    ( fac * theta(:,j) + thetapar1 * unp1gradtheta(:,j) )  &
                                   * detF * wg )
          end do
        end do

        if ( coefficients%i(91) == 1 ) then

!         Jacobian of SUPG test function

!         L^n.c^n + c^n.L^nT + f(c^n)
          call rhs_viscoelastic ( gvecn, cng, fng, vemopt=vemoptn )

          if ( coorsys == 1 .and. vel3D == 1 ) then
!           axisymmetric with 3D velocities: add terms due to u.nabla c
            call add_ugrad_terms ( xg, uvecn, cng, fng )
          end if

          do m = mode1, mode2
            rhsmodel(:,:,m) = thetapar1*fiterg(:,:,m) + thetapar2*fng(:,:,m)
          end do

          do l = 1, ndim
            do k = 1, ndf
              do m = mode1, mode2
                m1 = m - mode1 + 1
                do j = 1, ncomp
                  do i = 1, ndfc
                    work3(i,j,m1,k,l) = &
                     - sum ( ( beta/2 * dhoverU(:,l) * ugradtheta_supg(:,i) &
                                 + tau * dthetadx(:,i,l) ) * phi(:,k) * &
                          ( - fac * citerg(:,j,m) + fac * cng(:,j,m) &
                                 - thetapar1*unp1gradciter(:,j,m) &
                                 - thetapar2*ungradcn(:,j,m) + &
                                      rhsmodel(:,j,m) ) * detF * wg )
                  end do
                end do
              end do
            end do
          end do

        end if

        if ( vemcompressible ) work5 = 0 ! Jacobian of J not implemented

        if ( vemgammap ) then

!         include Jacobian of gammap

          do k = 1, ndfc
            do m = mode1, mode2
              m1 = m - mode1 + 1
              do j = 1, ncomp
                do i = 1, ndfc
                  work12(i,j,m1,k) = &
                   - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                              thetapar1 * vemoptiter%drhsdgammap(:,j,m) * &
                                  theta(:,k) * detF * wg )
                end do
              end do
            end do
          end do

        end if

      case(10)

!       implicit Gear

!       L^n+1.citer + citer.L^n+1T + f(citer)
        call rhs_viscoelastic ( gvecnp1, citerg, fiterg, vemoptiter )

        if ( coorsys == 1 .and. vel3D == 1 ) then
!         axisymmetric with 3D velocities: add terms due to u.nabla c
          call add_ugrad_terms ( xg, uvecnp1, citerg, fiterg, &
            dfadd_mm=vemoptiter%drhs_mm )
        end if

!       Jacobian of -rhs

        do n = mode1, mode2
          n1 = n - mode1 + 1
          do l = 1, ncomp
            do k = 1, ndfc
              do m = mode1, mode2
                m1 = m - mode1 + 1
                do j = 1, ncomp
                  do i = 1, ndfc
                    work14(i,j,m1,k,l,n1) = &
                      - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                        vemoptiter%drhs_mm(:,j,m,l,n) * theta(:,k) * detF * wg )
                  end do
                end do
              end do
            end do
          end do
        end do

!       Jacobian of material time derivative (identical for all components)

        do i = 1, ndfc
          do j = 1, ndfc
            work11(i,j) = sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
              ( 1.5_dp * fac * theta(:,j) + unp1gradtheta(:,j) ) &
                                   * detF * wg )
          end do
        end do

        if ( coefficients%i(91) == 1 ) then

!         Jacobian of SUPG test function

          do l = 1, ndim
            do k = 1, ndf
              do m = mode1, mode2
                m1 = m - mode1 + 1
                do j = 1, ncomp
                  do i = 1, ndfc
                    work3(i,j,m1,k,l) = &
                     - sum ( ( beta/2 * dhoverU(:,l) * ugradtheta_supg(:,i) &
                                 + tau * dthetadx(:,i,l) ) * phi(:,k) * &
                              ( - 1.5_dp * fac * citerg(:,j,m) + &
                                       2 * fac * cng(:,j,m) &
                                - 0.5_dp * fac * cnm1g(:,j,m) &
                                  - unp1gradciter(:,j,m) + &
                                              fiterg(:,j,m) ) * detF * wg )
                  end do
                end do
              end do
            end do
          end do

        end if

        if ( vemcompressible ) then

!         include Jacobian of J

          do k = 1, ndfp
            do m = mode1, mode2
              m1 = m - mode1 + 1
              do j = 1, ncomp
                do i = 1, ndfc
                  work5(i,j,m1,k) = &
                    sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                               vemoptiter%drhsdJ(:,j,m) * &
                                  vemoptiter%J * psi(:,k) * detF * wg ) / Kmod
                end do
              end do
            end do
          end do

        end if

        if ( vemgammap ) then

!         include Jacobian of gammap

          do k = 1, ndfc
            do m = mode1, mode2
              m1 = m - mode1 + 1
              do j = 1, ncomp
                do i = 1, ndfc
                  work12(i,j,m1,k) = &
                   - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                               vemoptiter%drhsdgammap(:,j,m) * &
                                  theta(:,k) * detF * wg )
                end do
              end do
            end do
          end do

        end if

      case default

        write(*,'(/a/a,i0/)') 'Error in implicit_ce_supg_elem:', &
        ' incorrect time integration scheme  = ', timeint
        stop

      end select

!     add diagonal blocks of material derivative

      do m1 = 1, mode2-mode1+1
        do k = 1, ncomp
          work14(:,k,m1,:,k,m1) = work14(:,k,m1,:,k,m1) + work11
        end do
      end do

!     compute column ranges for pressure(J), gammap, velocity, zero, c
!                               work5,       work12, work3,    zero, work14

      j5 = [ 1, ndfp ]
      j12 = [ 1, ndfc ]
      j3 = [ 1, ndf*ndim ]
      jz = [ ndf*ndim+1, ndf*ncompu ]
      j14 = [ 1, ndfc*ncomp*nsubm ]

      if ( coefficients%i(91) == 1 ) then
!       include Jacobian of SUPG test function: work3, zero, work14
        j14 = j14 + jz(2)
      end if

      if ( vemcompressible ) then
!       include Jacobian of J: set work5 to first block
        j12 = j12 + ndfp
        j3 = j3 + ndfp
        jz = jz + ndfp
        j14 = j14 + ndfp
      end if

      if ( vemgammap ) then
!       include Jacobian of gammap: set work12 second block
        j3 = j3 + ndfc
        jz = jz + ndfc
        j14 = j14 + ndfc
      end if

!     fill elemmat for conformation

      elemmat(:,j14(1):j14(2)) = &
                     reshape ( work14, [ ndfc*ncomp*nsubm, ndfc*ncomp*nsubm ] )

      if ( coefficients%i(91) == 1 ) then
!       include Jacobian of SUPG test function
        elemmat(:,j3(1):j3(2)) = &
                              reshape ( work3, [ ndfc*ncomp*nsubm, ndf*ndim ] )
        elemmat(:,jz(1):jz(2)) = 0
      end if

      if ( vemcompressible ) then
!       include Jacobian of J
        elemmat(:,j5(1):j5(2)) = reshape ( work5, [ ndfc*ncomp*nsubm, ndfp ] )
      end if

      if ( vemgammap ) then
!       include Jacobian of gammap
        elemmat(:,j12(1):j12(2)) = &
                                reshape ( work12, [ ndfc*ncomp*nsubm, ndfc ] )
      end if

    end if

    if ( vector ) then

      select case ( timeint )

      case(8)

!       implicit Euler

        if ( .not. matrix ) then

!         L^n+1.citer + citer.L^n+1T + f(citer)
          call rhs_viscoelastic ( gvecnp1, citerg, fiterg, vemoptiter )

          if ( coorsys == 1 .and. vel3D == 1 ) then
!           axisymmetric with 3D velocities: add terms due to u.nabla c
            call add_ugrad_terms ( xg, uvecnp1, citerg, fiterg )
          end if

        end if

        do m = mode1, mode2
          m1 = m - mode1 + 1
          do j = 1, ncomp
            do i = 1, ndfc
              work6(i,j,m1) = &
                sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                        ( - fac * citerg(:,j,m) + fac * cng(:,j,m) &
                            - unp1gradciter(:,j,m) + &
                                   fiterg(:,j,m) ) * detF * wg )
            end do
          end do
        end do

      case(9,11)

!       Crank-Nicolson (trapezoidal) or theta method

        if ( .not. matrix ) then

!         L^n+1.citer + citer.L^n+1T + f(citer)
          call rhs_viscoelastic ( gvecnp1, citerg, fiterg, vemoptiter )

          if ( coorsys == 1 .and. vel3D == 1 ) then
!           axisymmetric with 3D velocities: add terms due to u.nabla c
            call add_ugrad_terms ( xg, uvecnp1, citerg, fiterg )
          end if

        end if

        if ( .not. ( matrix .and. coefficients%i(91) == 1 ) ) then

!         L^n.c^n + c^n.L^nT + f(c^n)
          call rhs_viscoelastic ( gvecn, cng, fng, vemoptn )

          if ( coorsys == 1 .and. vel3D == 1 ) then
!           axisymmetric with 3D velocities: add terms due to u.nabla c
            call add_ugrad_terms ( xg, uvecn, cng, fng )
          end if

        end if

        do m = mode1, mode2
          rhsmodel(:,:,m) = thetapar1*fiterg(:,:,m) + thetapar2*fng(:,:,m)
          m1 = m - mode1 + 1
          do j = 1, ncomp
            do i = 1, ndfc
              work6(i,j,m1) = &
                sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                      ( - fac * citerg(:,j,m) + fac * cng(:,j,m) &
               - thetapar1*unp1gradciter(:,j,m) - thetapar2*ungradcn(:,j,m) + &
                                  rhsmodel(:,j,m) ) * detF * wg )
            end do
          end do
        end do

      case(10)

!       semi-implicit Gear with relaxation prediction

        if ( .not. matrix ) then

!         L^n+1.citer + citer.L^n+1T + f(citer)
          call rhs_viscoelastic ( gvecnp1, citerg, fiterg, vemoptiter )

          if ( coorsys == 1 .and. vel3D == 1 ) then
!           axisymmetric with 3D velocities: add terms due to u.nabla c
            call add_ugrad_terms ( xg, uvecnp1, citerg, fiterg )
          end if

        end if

        do m = mode1, mode2
          m1 = m - mode1 + 1
          do j = 1, ncomp
            do i = 1, ndfc
              work6(i,j,m1) = &
                sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                        ( - 1.5_dp * fac * citerg(:,j,m) + &
                                 2 * fac * cng(:,j,m) &
                          - 0.5_dp * fac * cnm1g(:,j,m) &
                            - unp1gradciter(:,j,m) + &
                                        fiterg(:,j,m) ) * detF * wg )
            end do
          end do
        end do

      case default

        write(*,'(/a/a,i0/)') 'Error in implicit_ce_supg_elem:', &
        ' incorrect time integration scheme  = ', timeint
        stop

      end select

      elemvec = reshape ( work6, [ ndfc * ncomp * nsubm ] )

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_devssg_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )


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

      allocate ( unp1gradtheta(ninti,ndfc), ugradtheta_supg(ninti,ndfc) )
      allocate ( ungradcn(ninti,ncomp,nmodes) )
      allocate ( u(ncompu*ndf), um(ndim*ndf), g(ncompg*ndfg) )
      allocate ( gradu(ninti,ncompu,ncompu) )
      allocate ( tmp(ndf,ndim), work2(ndfg,ncompg), work4(ndf,ncompu) )
      allocate ( work10(ninti,ncompg) )
      allocate ( posu(ncompu*ndf), posg(ncompg*ndfg) )
      allocate ( tau(ninti), hoverU(ninti) )

      allocate ( uvecnp1(ninti,ncompu), cn(ndfc,ncomp) )
      allocate ( uvec_supg(ninti,ncompu) )
      allocate ( cng(ninti,ncomp,nmodes), fng(ninti,ncomp,nmodes) )
      allocate ( uvecn(ninti,ncompu), cnm1(ndfc,ncomp) )
      allocate ( cnm1g(ninti,ncomp,nmodes), fnm1g(ninti,ncomp,nmodes) )
      allocate ( RT(ninti,ncompu,ncompu) )
      allocate ( work14(ndfc,ncomp,nsubm,ndfc,ncomp,nsubm) )
      allocate ( work6(ndfc,ncomp,nsubm) )
      allocate ( facv(ninti) )
      allocate ( work11(ndfc,ndfc) )

      if ( coorsys == 1 .and. vel3D == 0 ) then
        allocate ( gvecnp1(ninti,ndim**2+1), gvecn(ninti,ndim**2+1) )
      else if ( vel3D == 1 ) then
        allocate ( gvecnp1(ninti,ncompu**2), gvecn(ninti,ncompu**2) )
      else
        allocate ( gvecnp1(ninti,ndim**2), gvecn(ninti,ndim**2) )
      end if

      if ( coefficients%i(91) == 1 ) then
        allocate ( dhoverU(ninti,ndim), work3(ndfc,ncomp,nsubm,ndf,ndim) )
      end if

      allocate ( ungradtheta(ninti,ndfc) )
      allocate ( citer(ndfc,ncomp) )
      allocate ( citerg(ninti,ncomp,nmodes) )
      allocate ( unp1gradciter(ninti,ncomp,nmodes) )
      allocate ( fiterg(ninti,ncomp,nmodes), rhsmodel(ninti,ncomp,nmodes) )
      allocate ( uvecmeshn(ninti,ndim), uvecmeshnp1(ninti,ndim) )

      if ( varpar ) then
        allocate ( mvemodel(ninti) )
        mvemodel = vemodel
      end if

      if ( vemcompressible ) then
        allocate ( pr(ndfp), press(ninti) )
        allocate ( work5(ndfc,ncomp,nsubm,ndfp) )
      end if

      call create_vemopt ( vemoptiter, dep_J=vemcompressible, &
        dep_gammap=vemgammap, compute_drhs_mm=.true., &
        compute_drhsdJ=vemcompressible, compute_drhsdgammap=vemgammap, &
        np=ninti, ncomp=ncomp, nmodes=nmodes )

      if ( any ( coefficients%i(22) == [9,11] ) ) then
        call create_vemopt ( vemoptn, dep_J=vemcompressible, &
          dep_gammap=vemgammap, np=ninti )
      end if

      if ( vemgammap ) then
        allocate ( gpiter(ndfc), gpiterg(ninti) )
        allocate ( work12(ndfc,ncomp,nsubm,ndfc) )
        if ( any ( coefficients%i(22) == [9,11] ) ) then
          allocate ( gpn(ndfc), gpng(ninti) )
        end if
      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      integer :: ip

      deallocate ( unp1gradtheta, ugradtheta_supg )
      deallocate ( ungradcn )
      deallocate ( u, um, g )
      deallocate ( gradu )
      deallocate ( tmp, work2, work4 )
      deallocate ( work10 )
      deallocate ( posu, posg )
      deallocate ( tau, hoverU )

      deallocate ( uvecnp1, cn )
      deallocate ( uvec_supg )
      deallocate ( cng, fng )
      deallocate ( uvecn, cnm1 )
      deallocate ( cnm1g, fnm1g )
      deallocate ( RT )
      deallocate ( work14 )
      deallocate ( work6 )
      deallocate ( facv )
      deallocate ( work11 )

      deallocate ( gvecnp1, gvecn )

      if ( coefficients%i(91) == 1 ) then
        deallocate ( dhoverU, work3 )
      end if

      deallocate ( ungradtheta )
      deallocate ( citer )
      deallocate ( citerg, fiterg, rhsmodel )
      deallocate ( unp1gradciter )
      deallocate ( uvecmeshn, uvecmeshnp1 )

      if ( varpar ) then
        do ip = 1, ninti
          call delete ( mvemodel(ip) )
        end do
        deallocate ( mvemodel )
      end if

      if ( vemcompressible ) then
        deallocate ( pr, press )
        deallocate ( work5 )
      end if

      call delete ( vemoptiter )

      if ( any ( coefficients%i(22) == [9,11] ) ) call delete ( vemoptn )

      if ( vemgammap ) then
        deallocate ( gpiter, gpiterg )
        deallocate ( work12 )
        if ( any ( coefficients%i(22) == [9,11] ) ) then
          deallocate ( gpn, gpng )
        end if
      end if

    end subroutine deallocate_arrays


!   get cn and cnm1 for the standard case (ALE or Euler) and citer

    subroutine get_c_standard

      integer :: m, isv_c

      if ( coefficients%i(89) == 1 ) then

!       no time-derivative included

        cng(:,:,mode1:mode2) = 0

      else

!       get conformation tensor at time step tn

        isv_c = get_coefficient ( coefficients, index=88, default=1 )

        do m = mode1, mode2

!         single mode conformation
          call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
            cmode=cn, mode=m, isv=isv_c )

          cng(:,:,m) = matmul ( theta, cn )

          if ( coefficients%i(90) == 1 ) then

!           rotation reinitialization in integration points

            if ( timeint == 10 ) then ! BDF2
              call sqrtc_b ( cng(:,:,m), RT=RT )
            else
              call sqrtc_b ( cng(:,:,m) )
            end if

          end if

        end do

      end if

      if ( timeint == 10 ) then ! BDF2

        if ( coefficients%i(89) == 1 ) then

!         no time-derivative included

          cnm1g(:,:,mode1:mode2) = 0

        else

!         get conformation tensor at time step tn-1

          do m = mode1, mode2

!           single mode conformation
            call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
              cmode=cnm1, mode=m, isv=isv_c+1 )

            cnm1g(:,:,m) = matmul ( theta, cnm1 )

            if ( coefficients%i(90) == 1 ) then

!             rotation reinitialization in integration points

              call sqrtc_b ( cnm1g(:,:,m) )
              call rotate_b ( cnm1g(:,:,m), RT )

            end if

          end do

        end if

      end if

!     get iteration of conformation tensor

      isv_c = get_coefficient ( coefficients, index=87, default=3 )

      do m = mode1, mode2

!       single mode conformation
        call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
          cmode=citer, mode=m, isv=isv_c )

        citerg(:,:,m) = matmul ( theta, citer )

        unp1gradciter(:,:,m) = matmul ( unp1gradtheta, citer )

!       un.grad cn term (CN an theta method only)

        if ( any ( timeint == [9,11] ) ) then
          ungradcn(:,:,m) = matmul ( ungradtheta, cn )
        end if

      end do

    end subroutine get_c_standard


!   get gpiter and gpn for CN

    subroutine get_gp_standard

      integer :: isv_gp

!     get iteration of equivalent plastic strain

      isv_gp = get_coefficient ( coefficients, index=101, default=1 )

      call get_gammap ( mesh, problem, elgrp, elem, oldvectors, gp=gpiter, &
        isv=isv_gp )

      gpiterg = matmul ( theta, gpiter )

      vemoptiter%gammap = gpiterg

      if ( any ( timeint == [9,11] ) ) then

!       get equivalent plastic strain at time step tn

        isv_gp = get_coefficient ( coefficients, index=102, default=2 )

        call get_gammap ( mesh, problem, elgrp, elem, oldvectors, gp=gpn, &
          isv=isv_gp )

        gpng = matmul ( theta, gpn )

        vemoptn%gammap = gpng

      end if

    end subroutine get_gp_standard

  end subroutine implicit_ce_supg_elem


! Internal element routine for the constitutive equation (multi mode).
! Implicit time-integration with Newton-Raphson iteration.
! This is the velocity and G gradient part (in Jacobian). There is no rhs.

  subroutine implicit_ce_vel_supg_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m
    use limits_m, only: BLOCK_IMPROPER_ALE
    use supg_utils_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, k, l, n, ll, nn, ip, timeint, m, m1, iv_supg
    real(dp) :: beta, deltat, esize, thetapar1


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_devssg_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1, singlemode=.true., &
      gammap=.true. )

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'implicit_ce_vel_supg_elem', &
        indexarray=[22,48,50,57], minimum=[8,0,0,0], &
        maximum=[11,1,1,0] )

      if ( coefficients%i(48) == 1 .and. any(coefficients%i(22) == [9,11]) &
           .and. BLOCK_IMPROPER_ALE ) then
        write(*,'(/6(a/))') 'Error in implicit_ce_vel_supg_elem:', &
          ' Combining ALE with Crank-Nicolson or theta method not allowed, ', &
          ' since the change in the nabla operator due to the change ', &
          ' in mesh at different times within a single step ', &
          ' has not been taken into account yet.', &
          ' Set BLOCK_IMPROPER_ALE=.false. in limits_m to bypass this block.'
        stop
      end if

      if ( varpar ) then

!       check coefficients for variable coefficients
        call check ( coefficients, 'implicit_ce_vel_supg_elem', &
          indexarray=[69,70,81], minimum=[0,0,0], maximum=[3,3,3] )

        if ( any(coefficients%i(22) == [9,11])) then
          write(*,'(/4(a/))') 'Error in implicit_ce_vel_supg_elem:', &
            ' Combining variable coefficients with Crank-Nicolson ', &
            ' or theta-method not allowed, since the variable coefficients ', &
            ' are only available for a single time step.'
          stop
        end if

      end if

      cstorage = get_coefficient ( coefficients, index=84, default=2 )

    end if

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of element

    timeint = coefficients%i(22)

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coefficients%i(61) == 1 ) then
      call shape_derivative ( dphi, Finv, dphidx )
    end if

    call shape_derivative ( dtheta, Finv, dthetadx )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if


!   get velocity vector u_supg

    iv_supg = get_coefficient ( coefficients, index=86, default=1 )

    call get_sysvector ( mesh, problem, oldvectors%s(iv_supg)%p, &
      elgrp, elem, u, physq=[physqvel] )

    work4 = reshape ( u, [ndf,ncompu] )

    uvec_supg = matmul ( phi, work4 )


!   set time step

    deltat = coefficients%r(8)


!   ALE, temporary ALE or Eulerian frame

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity at tn+1

      call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, &
        elem, um )

      tmp = reshape ( um, [ndf,ndim] )

      uvecmeshnp1 = matmul ( phi, tmp )

!     (u_supg-ugridnp1).grad operator

      do ip = 1, ninti
        ugradtheta_supg(ip,:) = &
         matmul ( dthetadx(ip,:,:), uvec_supg(ip,1:ndim) - uvecmeshnp1(ip,:) )
      end do

    else

!     Eulerian frame

      uvecmeshnp1 = 0

!     (u_supg-ugridnp1).grad operator

      do ip = 1, ninti
        ugradtheta_supg(ip,:) = &
                        matmul ( dthetadx(ip,:,:), uvec_supg(ip,1:ndim) )
      end do

    end if


!   get conformation tensor citer and gradient of citer

    call get_c_standard


!   upwinding parameter and other factors

    beta = coefficients%r(9)

    if ( htype == 3 ) then
      esize = sum ( detF * wg )
    else
      esize = 1
    end if

!   thetapar

    if ( timeint == 9 ) then
      thetapar1 = 0.5_dp
    else if ( timeint == 11 ) then
      thetapar1 = coefficients%r(28)
    end if

!   compute h/U

    call supg_hU ( ndim, globalshape, hoverU, &
      htype=htype, Uscaling=Uscaling, Uglobal=coefficients%r(10), &
      uvec=uvec_supg(:,1:ndim)-uvecmeshnp1, Finv=Finv, x=x, esize=esize, &
      checkzero=coefficients%i(78)==1 )

    if ( coefficients%i(50) == 1 ) then
!     Courant number dependent tau
      where ( deltat < hoverU )
        hoverU = deltat
      end where
    end if

    tau = beta * hoverU / 2  ! upwinding parameter

!   build equations

    if ( matrix ) then

!     Jacobian of u.grad c

      if ( coorsys == 1 .and. vel3D == 1 ) then
!       axisymmetric with 3D velocities: add terms due to u.nabla c
        call add_swirl_terms ( xg, citerg, gradciter )
      end if

!     delta u.grad c

      do l = 1, ncompu
        do k = 1, ndf
          do m = mode1, mode2
            m1 = m - mode1 + 1
            do j = 1, ncomp
              do i = 1, ndfc
                work15(i,j,m1,k,l) = &
                sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                         phi(:,k) * gradciter(:,j,l,m) * detF * wg )
              end do
            end do
          end do
        end do
      end do

!     bilinear terms

      call drhsdL_viscoelastic ( citerg, drhsdLg )

      if ( coefficients%i(61) == 0 ) then

!       G-method

        do l = 1, ndim       ! col of L = coor dir
          do n = 1, ncompu   ! row of L = vel comp
            ll = (n-1)*ncompu + l    ! L component number (row-major)
            nn = (n-1)*ndim + l      ! G component number (row-major)
            do k = 1, ndfg
              do m = mode1, mode2
                m1 = m - mode1 + 1
                do j = 1, ncomp
                  do i = 1, ndfc
                    work3(i,j,m1,k,nn) = &
                       - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                                  zeta(:,k) * drhsdLg(:,j,ll,m) * detF * wg )
                  end do
                end do
              end do
            end do
          end do
        end do

      else if ( coefficients%i(61) == 1 ) then

!       direct velocity gradient

        do l = 1, ndim       ! col of L = coor dir
          do n = 1, ncompu   ! row of L = vel comp
            ll = (n-1)*ncompu + l  ! L component number (row-major)
            do k = 1, ndf
              do m = mode1, mode2
                m1 = m - mode1 + 1
                do j = 1, ncomp
                  do i = 1, ndfc
                    work15(i,j,m1,k,n) = work15(i,j,m1,k,n) - &
                      sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                               dphidx(:,k,l) * drhsdLg(:,j,ll,m) * detF * wg )
                  end do
                end do
              end do
            end do
          end do
        end do

      end if

!     Special cases, contribution to velocity columns only

      if ( coorsys == 1 .and. vel3D == 0 ) then

!       u_r / r

        do k = 1, ndf
          do m = mode1, mode2
            m1 = m - mode1 + 1
            do j = 1, ncomp
              do i = 1, ndfc
                work15(i,j,m1,k,2) = work15(i,j,m1,k,2) - &
                  sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                           phi(:,k) / xg(:,2) * drhsdLg(:,j,5,m) * detF * wg )
              end do
            end do
          end do
        end do

      else if ( coorsys == 1 .and. vel3D == 1 ) then

        do k = 1, ndf
          do m = mode1, mode2
            m1 = m - mode1 + 1
            do j = 1, ncomp
              do i = 1, ndfc
!              - u_theta / r
                work15(i,j,m1,k,3) = work15(i,j,m1,k,3) + &
                    sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                             phi(:,k) / xg(:,2) * drhsdLg(:,j,6,m) * detF * wg )
!               u_r / r
                work15(i,j,m1,k,2) = work15(i,j,m1,k,2) - &
                    sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                             phi(:,k) / xg(:,2) * drhsdLg(:,j,9,m) * detF * wg )
              end do
            end do
          end do
        end do

      end if

!     factor for timeint

      if ( any ( timeint == [9,11] ) ) then

!       Crank-Nicolson (trapezoidal) or theta method

        work15 = thetapar1 * work15
        if ( coefficients%i(61) == 0 ) work3 = thetapar1 * work3

      end if

      if ( coefficients%i(61) == 0 ) then

!       G-method

        elemmat(:,:ndfg*ncompg) = &
                         reshape ( work3, [ ndfc*ncomp*nsubm, ndfg*ncompg ] )
        elemmat(:,ndfg*ncompg+1:) = &
                         reshape ( work15, [ ndfc*ncomp*nsubm, ndf*ncompu ] )

      else if ( coefficients%i(61) == 1 ) then

!       direct velocity gradient

        elemmat = reshape ( work15, [ ndfc*ncomp*nsubm, ndf*ncompu ] )

      end if

    end if

    if ( vector ) then

      elemvec = 0

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_devssg_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )


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

      allocate ( ugradtheta_supg(ninti,ndfc) )
      allocate ( u(ncompu*ndf), um(ndim*ndf) )
      allocate ( tmp(ndf,ndim), work4(ndf,ncompu) )
      allocate ( tau(ninti), hoverU(ninti) )

      allocate ( uvec_supg(ninti,ncompu) )
      allocate ( facv(ninti) )

      allocate ( citer(ndfc,ncomp) )
      allocate ( citerg(ninti,ncomp,nmodes) )
      allocate ( uvecmeshn(ninti,ndim), uvecmeshnp1(ninti,ndim) )

      if ( varpar ) then
        allocate ( mvemodel(ninti) )
        mvemodel = vemodel
      end if

      allocate ( gradciter(ninti,ncomp,ncompu,nmodes) )
      allocate ( work15(ndfc,ncomp,nsubm,ndf,ncompu) )
      allocate ( work3(ndfc,ncomp,nsubm,ndfg,ncompg) )
      if ( coorsys == 1 .and. vel3D == 0 ) then
        allocate ( drhsdLg(ninti,ncomp,ndim**2+1,nmodes) )
      else if ( vel3D == 1 ) then
        allocate ( drhsdLg(ninti,ncomp,ncompu**2,nmodes) )
      else
        allocate ( drhsdLg(ninti,ncomp,ndim**2,nmodes) )
      end if


    end subroutine allocate_arrays

    subroutine deallocate_arrays

      integer :: ip

      deallocate ( ugradtheta_supg )
      deallocate ( u, um )
      deallocate ( tmp, work4 )
      deallocate ( tau, hoverU )

      deallocate ( uvec_supg )
      deallocate ( facv )

      deallocate ( citer )
      deallocate ( citerg )
      deallocate ( uvecmeshn, uvecmeshnp1 )

      if ( varpar ) then
        do ip = 1, ninti
          call delete ( mvemodel(ip) )
        end do
        deallocate ( mvemodel )
      end if

      deallocate ( gradciter )
      deallocate ( work15 )
      deallocate ( work3 )
      deallocate ( drhsdLg )

    end subroutine deallocate_arrays


!   get citer

    subroutine get_c_standard

      integer :: m, i, isv_c

!     get iteration of conformation tensor

      isv_c = get_coefficient ( coefficients, index=87, default=1 )

      do m = mode1, mode2

!       single mode conformation
        call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
          cmode=citer, mode=m, isv=isv_c )

        citerg(:,:,m) = matmul ( theta, citer )

!       grad citer term for mode m (up to ndim components)

        do i = 1, ndim
          gradciter(:,:,i,m) = matmul ( dthetadx(:,:,i), citer )
        end do
        gradciter(:,:,ndim+1:ncompu,m) = 0

      end do

    end subroutine get_c_standard

  end subroutine implicit_ce_vel_supg_elem


! cylindrical, axisymmetric with 3D velocities: add terms due to u.nabla c

  subroutine add_swirl_terms ( xg, cvec, dfadd )

    use viscoelastic_globals_m, only: mode1, mode2, bten, facv

    real(dp), dimension(:,:), intent(in) :: xg
    real(dp), dimension(:,:,:), intent(in) :: cvec
    real(dp), dimension(:,:,:,:), intent(inout) :: dfadd

    integer :: m

    facv = - 1 / xg(:,2)  ! - u_theta / r

    if ( bten ) then

!     b-tensor formulation

      do m = mode1, mode2
        dfadd(:,2,3,m) = dfadd(:,2,3,m) - facv * cvec(:,3,m)     ! - b_z,theta
        dfadd(:,3,3,m) = dfadd(:,3,3,m) + facv * cvec(:,2,m)     ! + b_z,r
        dfadd(:,4,3,m) = dfadd(:,4,3,m) - facv * cvec(:,7,m)     ! - b_theta,z
        dfadd(:,5,3,m) = dfadd(:,5,3,m) &
             - facv * ( cvec(:,6,m) + cvec(:,8,m) ) ! - b_r,theta - b_theta,r
        dfadd(:,6,3,m) = dfadd(:,6,3,m) &
             + facv * ( cvec(:,5,m) - cvec(:,9,m) ) ! + b_rr - b_theta,theta
        dfadd(:,7,3,m) = dfadd(:,7,3,m) + facv * cvec(:,4,m)     ! + b_r,z
        dfadd(:,8,3,m) = dfadd(:,8,3,m) &
             + facv * ( cvec(:,5,m) - cvec(:,9,m) ) ! + b_rr - b_theta,theta
        dfadd(:,9,3,m) = dfadd(:,9,3,m) &
             + facv * ( cvec(:,6,m) + cvec(:,8,m) ) ! + b_r,theta + b_theta,r
      end do

    else

!     c-tensor formulation, including logc

      do m = mode1, mode2
        dfadd(:,2,3,m) = dfadd(:,2,3,m) - facv * cvec(:,3,m)     ! - c_z,theta
        dfadd(:,3,3,m) = dfadd(:,3,3,m) + facv * cvec(:,2,m)     ! + c_z,r
        dfadd(:,4,3,m) = dfadd(:,4,3,m) - 2 * facv * cvec(:,5,m) ! - 2 c_r,theta
        dfadd(:,5,3,m) = dfadd(:,5,3,m) &
                  + facv * ( cvec(:,4,m) - cvec(:,6,m) ) ! + c_rr-c_theta,theta
        dfadd(:,6,3,m) = dfadd(:,6,3,m) + 2 * facv * cvec(:,5,m) ! + 2 c_r,theta
      end do

    end if

  end subroutine add_swirl_terms


! derivative of right-hand of viscoelastic model side wrt L
! NOTE: this is not a general usable subroutine, but more like an inlined
!       procedure. Most in/out variables are in a global module and need to
!       be defined/allocated.

  subroutine drhsdL_viscoelastic ( cvec, drhsdL )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      coorsys, vel3D, varpar, logc, bten, vemodel, mode1, mode2, mvemodel

    real(dp), dimension(:,:,:), intent(in) :: cvec
    real(dp), dimension(:,:,:,:), optional, intent(out) :: drhsdL

    if ( logc ) then
      write(*,'(/2(a/))') 'Error in drhsdL_viscoelastic:', &
        ' derivative drhsdL for log c formulation not available.'
      stop
    end if

    if ( coorsys <= 1 .and. vel3D == 0 ) then

!     2D and axisymmetric

      if ( varpar ) then

!       variable coefficients

        if ( bten ) then
          call drhsdL_viscoelastic_2D_b ( vemodel, cvec, drhsdL, &
            mode1=mode1, mode2=mode2, mvemodel=mvemodel )
        else
          call drhsdL_viscoelastic_2D ( vemodel, cvec, drhsdL, &
            mode1=mode1, mode2=mode2, mvemodel=mvemodel )
        end if

      else

!       fixed coefficients

        if ( bten ) then
          call drhsdL_viscoelastic_2D_b ( vemodel, cvec, drhsdL, &
            mode1=mode1, mode2=mode2 )
        else
          call drhsdL_viscoelastic_2D ( vemodel, cvec, drhsdL, &
            mode1=mode1, mode2=mode2 )
        end if

      end if

    else if ( coorsys == 2 .or. vel3D == 1 ) then

!     3D

      if ( varpar ) then

!       variable coefficients

        if ( bten ) then
          call drhsdL_viscoelastic_3D_b ( vemodel, cvec, drhsdL, &
            mode1=mode1, mode2=mode2, mvemodel=mvemodel )
        else
          call drhsdL_viscoelastic_3D ( vemodel, cvec, drhsdL, &
            mode1=mode1, mode2=mode2, mvemodel=mvemodel )
        end if

      else

!       fixed coefficients

        if ( bten ) then
          call drhsdL_viscoelastic_3D_b ( vemodel, cvec, drhsdL, &
            mode1=mode1, mode2=mode2 )
        else
          call drhsdL_viscoelastic_3D ( vemodel, cvec, drhsdL, &
            mode1=mode1, mode2=mode2 )
        end if

      end if

    end if

  end subroutine drhsdL_viscoelastic


! Internal element routine for the constitutive equation (multi-mode).
! Matrix of the time-derivative only.

  subroutine implicit_ce_timederiv_supg_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m
    use supg_utils_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, k, ip, iv_supg, m1
    real(dp) :: beta, esize


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1, singlemode=.true., &
      gammap=.true. )

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'implicit_ce_timederiv_supg_elem', &
        indexarray=[48,50,57], minimum=[0,0,0], maximum=[1,0,0] )

    end if

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

    call shape_derivative ( dtheta, Finv, dthetadx )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

!   get velocity vector u_supg

    iv_supg = get_coefficient ( coefficients, index=86, default=1 )

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(iv_supg)%p, &
      elgrp, elem, u, physq=[physqvel], posu=posu )

    work4 = reshape ( u, [ndf,ncompu] )

    uvec_supg = matmul ( phi, work4 )


!   ALE, temporary ALE or Eulerian frame

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity at tn+1

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, um )

      tmp = reshape ( um, [ndf,ndim] )

      uvecmeshnp1 = matmul ( phi, tmp )

!     (u_supg-ugridnp1).grad operator

      do ip = 1, ninti
        ugradtheta_supg(ip,:) = &
         matmul ( dthetadx(ip,:,:), uvec_supg(ip,1:ndim) - uvecmeshnp1(ip,:) )
      end do

    else

!     Eulerian frame

      uvecmeshnp1 = 0

!     (u_supg-ugridnp1).grad operator

      do ip = 1, ninti
        ugradtheta_supg(ip,:) = &
                        matmul ( dthetadx(ip,:,:), uvec_supg(ip,1:ndim) )
      end do

    end if


!   upwinding parameter and other factors

    beta = coefficients%r(9)

    if ( htype == 3 ) then
      esize = sum ( detF * wg )
    else
      esize = 1
    end if

!   compute h/U

    call supg_hU ( ndim, globalshape, hoverU, &
      htype=htype, Uscaling=Uscaling, Uglobal=coefficients%r(10), &
      uvec=uvec_supg(:,1:ndim)-uvecmeshnp1, Finv=Finv, x=x, esize=esize, &
      checkzero=coefficients%i(78)==1 )

    tau = beta * hoverU / 2  ! upwinding parameter

!   build matrix

    if ( matrix ) then

!     matrix of time derivative (identical for all components)

      do j = 1, ndfc
        do i = 1, ndfc
          work11(i,j) = sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                                theta(:,j) * detF * wg )
        end do
      end do

!     set diagonal blocks of material derivative

      work14 = 0

      do m1 = 1, mode2-mode1+1
        do k = 1, ncomp
          work14(:,k,m1,:,k,m1) = work11
        end do
      end do

      elemmat = reshape ( work14, [ ndfc*ncomp*nsubm, ndfc*ncomp*nsubm ] )

    end if

    if ( vector ) then

      write(*,'(/a/a,i0/)') 'Error in implicit_ce_timederiv_supg_elem:', &
        ' No vector available. Call build_system with buildvector=.false.'
      stop

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )


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

      allocate ( ugradtheta_supg(ninti,ndfc) )
      allocate ( u(ncompu*ndf), um(ndim*ndf) )
      allocate ( tmp(ndf,ndim), work4(ndf,ncompu) )
      allocate ( posu(ncompu*ndf) )
      allocate ( tau(ninti), hoverU(ninti) )

      allocate ( uvec_supg(ninti,ncompu) )
      allocate ( work14(ndfc,ncomp,nsubm,ndfc,ncomp,nsubm) )
      allocate ( work11(ndfc,ndfc) )

      allocate ( uvecmeshnp1(ninti,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( ugradtheta_supg )
      deallocate ( u, um )
      deallocate ( tmp, work4 )
      deallocate ( posu )
      deallocate ( tau, hoverU )

      deallocate ( uvec_supg )
      deallocate ( work14 )
      deallocate ( work11 )

      deallocate ( uvecmeshnp1 )

    end subroutine deallocate_arrays

  end subroutine implicit_ce_timederiv_supg_elem


! Internal element routine for the steady time-derivative in the
! constitutive equation (multi mode). This is for the b-formulation only.
! Newton-Raphson iteration: produces a Jacobian and a right-hand side for
! the conformation part (b-tensor only). Only BDF1 and BDF2 available.

  subroutine steady_ce_timederiv_supg_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m
    use supg_utils_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, k, l, ip, timeint, m, m1, iv_supg
    real(dp) :: beta, deltat, esize
    type(tdpar_t) :: tdpar

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1, singlemode=.true., &
      gammap=.true. )

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'steady_ce_timederiv_supg_elem', &
        indexarray=[22,48,50,57,91], minimum=[8,0,0,0,0], maximum=[11,1,1,0,1] )

      if ( coefficients%i(90) /= 1 ) then
        write(*,'(/2(a/))') 'Error in steady_ce_timederiv_supg_elem:', &
          ' Pointwise rotation reinitialization not set.'
        stop
      end if

      if ( coefficients%i(71) /= 1 ) then
        write(*,'(/2(a/))') 'Error in steady_ce_timederiv_supg_elem:', &
          ' Contravariant deformation tensor (CDT) formulation not set. '
        stop
      end if

      if ( all(coefficients%i(22) /= [8,10]) ) then
        write(*,'(/a/a,i0,a/)') 'Error in steady_ce_timederiv_supg_elem:', &
          ' Steady state time derivative for time integration method = ', &
            coefficients%i(22), ' not available.'
        stop
      end if

      if ( coefficients%i(50) == 1 .and. coefficients%i(91) == 1 ) then
        write(*,'(/3(a/))') 'Error in steady_ce_timederiv_supg_elem:', &
          ' Jacobian for SUPG test function not available if the ', &
          ' tau parameter is adjusted based on the Courant number.'
        stop
      end if

      cstorage = get_coefficient ( coefficients, index=84, default=2 )

    end if

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of element

    timeint = coefficients%i(22)

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coefficients%i(61) == 1 ) then
      call shape_derivative ( dphi, Finv, dphidx )
    end if

    call shape_derivative ( dtheta, Finv, dthetadx )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if


!   get velocity vector u_supg

    iv_supg = get_coefficient ( coefficients, index=86, default=1 )

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(iv_supg)%p, &
      elgrp, elem, u, physq=[physqvel], posu=posu )

    work4 = reshape ( u, [ndf,ncompu] )

    uvec_supg = matmul ( phi, work4 )


!   set timestep and tdpar

    deltat = coefficients%r(8)
    tdpar%timestep = deltat
    select case ( timeint )
    case(8)
!     implicit Euler (BDF1)
      tdpar%method = 1
    case(10)
!     implicit Gear (BDF2)
      tdpar%method = 2
    case default
      call errormsg_case_default ( 'steady_ce_timederiv_supg_elem', &
        'timeint', int_value=timeint )
    end select


!   ALE, temporary ALE or Eulerian frame

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity at tn+1

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, um )

      tmp = reshape ( um, [ndf,ndim] )

      uvecmeshnp1 = matmul ( phi, tmp )

!     (u_supg-ugridnp1).grad operator

      do ip = 1, ninti
        ugradtheta_supg(ip,:) = &
         matmul ( dthetadx(ip,:,:), uvec_supg(ip,1:ndim) - uvecmeshnp1(ip,:) )
      end do

    else

!     Eulerian frame

      uvecmeshnp1 = 0

!     (u_supg-ugridnp1).grad operator

      do ip = 1, ninti
        ugradtheta_supg(ip,:) = &
                        matmul ( dthetadx(ip,:,:), uvec_supg(ip,1:ndim) )
      end do

    end if


!   get conformation tensor at previous time steps (cn, cnm1) and citer

!   standard case (ALE or Euler)

    call get_c_standard

!   upwinding parameter and other factors

    beta = coefficients%r(9)

    if ( htype == 3 ) then
      esize = sum ( detF * wg )
    else
      esize = 1
    end if

!   compute h/U

    if ( coefficients%i(91) == 0 ) then
!     default h/U
      call supg_hU ( ndim, globalshape, hoverU, &
        htype=htype, Uscaling=Uscaling, Uglobal=coefficients%r(10), &
        uvec=uvec_supg(:,1:ndim)-uvecmeshnp1, Finv=Finv, x=x, esize=esize, &
        checkzero=coefficients%i(78)==1 )
    else if ( coefficients%i(91) == 1 ) then
!     include Jacobian of SUPG test function
      call supg_hU ( ndim, globalshape, hoverU, &
        htype=htype, Uscaling=Uscaling, Uglobal=coefficients%r(10), &
        uvec=uvec_supg(:,1:ndim)-uvecmeshnp1, Finv=Finv, x=x, esize=esize, &
        checkzero=coefficients%i(78)==1, dhoverU=dhoverU )
    end if

    if ( coefficients%i(50) == 1 ) then
!     Courant number dependent tau
      where ( deltat < hoverU )
        hoverU = deltat
      end where
    end if

    tau = beta * hoverU / 2  ! upwinding parameter

!   build equations

    if ( matrix ) then

!     lhs and Jacobian of steady time derivative in integration points
      call NRtd_steady_viscoelastic ( tdpar, bvec=citerg, Hb=fiterg, &
        dHb=dfiterg )

!     Jacobian of lhs

      work14 = 0

      do m = mode1, mode2
        m1 = m - mode1 + 1
        do l = 1, ncomp
          do k = 1, ndfc
            do j = 1, ncomp
              do i = 1, ndfc
                work14(i,j,m1,k,l,m1) = &
                  sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                              dfiterg(:,j,l,m) * theta(:,k) * detF * wg )
              end do
            end do
          end do
        end do
      end do

      if ( coefficients%i(91) == 1 ) then

!       Jacobian of SUPG test function

        do l = 1, ndim
          do k = 1, ndf
            do m = mode1, mode2
              m1 = m - mode1 + 1
              do j = 1, ncomp
                do i = 1, ndfc
                  work3(i,j,m1,k,l) = &
                      sum ( ( beta/2 * dhoverU(:,l) * ugradtheta_supg(:,i) &
                                  + tau * dthetadx(:,i,l) ) * phi(:,k) * &
                                               fiterg(:,j,m) * detF * wg )
                end do
              end do
            end do
          end do
       end do

     end if

     if ( coefficients%i(91) == 1 ) then
!      include Jacobian of SUPG test function
       elemmat(:,:ndf*ndim) = reshape ( work3, [ ndfc*ncomp*nsubm, ndf*ndim ] )
       elemmat(:,ndf*ndim+1:ndf*ncompu) = 0
       elemmat(:,ndf*ncompu+1:) = &
                     reshape ( work14, [ ndfc*ncomp*nsubm, ndfc*ncomp*nsubm ] )
     else
       elemmat = reshape ( work14, [ ndfc*ncomp*nsubm, ndfc*ncomp*nsubm ] )
     end if

    end if

    if ( vector ) then

      if ( .not. matrix ) then

!       lhs and Jacobian of steady time derivative in integration points
        call NRtd_steady_viscoelastic ( tdpar, bvec=citerg, Hb=fiterg, &
          dHb=dfiterg )

      endif

!     - lhs
      do m = mode1, mode2
        m1 = m - mode1 + 1
        do j = 1, ncomp
          do i = 1, ndfc
            work6(i,j,m1) = &
              - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                                 fiterg(:,j,m) * detF * wg )
          end do
        end do
      end do

      elemvec = reshape ( work6, [ ndfc*ncomp*nsubm ] )

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )


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

      allocate ( u(ncompu*ndf) )
      allocate ( ugradtheta_supg(ninti,ndfc) )
      allocate ( tau(ninti), hoverU(ninti) )

      allocate ( uvec_supg(ninti,ncompu) )
      allocate ( work4(ndf,ncompu) )
      allocate ( work14(ndfc,ncomp,nsubm,ndfc,ncomp,nsubm) )
      allocate ( work6(ndfc,ncomp,nsubm) )

      if ( coefficients%i(91) == 1 ) then
        allocate ( dhoverU(ninti,ndim), work3(ndfc,ncomp,nsubm,ndf,ndim) )
      end if

      allocate ( citer(ndfc,ncomp) )
      allocate ( citerg(ninti,ncomp,nmodes) )
      allocate ( fiterg(ninti,ncomp,nmodes) )
      allocate ( dfiterg(ninti,ncomp,ncomp,nmodes) )
      allocate ( uvecmeshnp1(ninti,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( u )
      deallocate ( ugradtheta_supg )
      deallocate ( tau, hoverU )

      deallocate ( uvec_supg )
      deallocate ( work4 )
      deallocate ( work14 )
      deallocate ( work6 )

      if ( coefficients%i(91) == 1 ) then
        deallocate ( dhoverU, work3 )
      end if

      deallocate ( citer )
      deallocate ( citerg, fiterg )
      deallocate ( dfiterg )
      deallocate ( uvecmeshnp1 )

    end subroutine deallocate_arrays


!   get cn and cnm1 for the standard case (ALE or Euler) and citer

    subroutine get_c_standard

      integer :: m, isv_c

!     get iteration of conformation tensor

      isv_c = get_coefficient ( coefficients, index=87, default=3 )

      do m = mode1, mode2

!       single mode conformation
        call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
          cmode=citer, mode=m, isv=isv_c )

        citerg(:,:,m) = matmul ( theta, citer )

      end do

    end subroutine get_c_standard

  end subroutine steady_ce_timederiv_supg_elem


! right-hand side of viscoelastic model
! NOTE: this is not a general usable subroutine, but more like an inlined
!       procedure. Most in/out variables are in a global module and need to
!       be defined/allocated.

  subroutine NRtd_steady_viscoelastic ( tdpar, bvec, Hb, dHb )

!   global variables:
    use viscoelastic_globals_m, only: &
!     input:
      coorsys, vel3D, vemodel, mode1, mode2

    real(dp), dimension(:,:,:), intent(in) :: bvec
    real(dp), dimension(:,:,:), intent(out) :: Hb
    real(dp), dimension(:,:,:,:), intent(out) :: dHb
    type(tdpar_t), intent(in) :: tdpar

    if ( coorsys <= 1 .and. vel3D == 0 ) then

!     2D and axisymmetric

      call NRtd_steady_viscoelastic_2D_b ( vemodel, tdpar, bvec, Hb, dHb, &
        mode1, mode2 )

    else if ( coorsys == 2 .or. vel3D == 1 ) then

!     3D

      call NRtd_steady_viscoelastic_3D_b ( vemodel, tdpar, bvec, Hb, dHb, &
        mode1, mode2 )

    end if

  end subroutine NRtd_steady_viscoelastic

end module viscoelastic_elements_implicit_m
