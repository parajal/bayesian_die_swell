
! Copyright (C) 2016-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Element routines for the convection-diffusion-reaction equation
!
!  gamma ( dc/dt + u . nabla c ) - alpha nabla^2 c + beta c = f
!
! NOTE: this routine is written for 2D/3D convection-diffusion problems
! with or without SUPG stabilization

module convection_diffusion_supg_elements_m

  use tfem_elem_m
  use stokes_set_globals_m

  implicit none

contains


! Internal element routine for the convection-diffusion-reaction equation

  subroutine conv_diff_supg_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use convection_diffusion_supg_globals_m
    use supg_utils_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip, timeint, nrhs
    real(dp) :: deltat, gamma, alpha, beta, fac
    real(dp) :: esize
    logical :: supg, diff, self

!   some tests

    if ( first ) then

      call check ( coefficients, 'conv_diff_supg_elem', &
        indexarray=[48,401,403,404,405,406,408,410,411,412,414], &
        minimum=[0,0,0,0,0,0,0,0,0,0,0], maximum=[1,2,1,1,1,4,3,4,5,1,4] )

      call check ( coefficients, 'conv_diff_supg_elem', indexarray=[413], &
        minimum=[0] )

    end if

    supg = coefficients%i(403) == 1 ! Galerkin or SUPG for spatial discetization
    diff = coefficients%i(404) == 1 ! diffusion term included ?
    self = coefficients%i(405) == 1 ! self-dependent source term included ?
    nrhs = coefficients%i(413)      ! number of RHS sides to assemble

!   set globals, gauss, shapefuncs, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, &
      first, last, coefficients, oldvectors, maxvel3D=1 )

    call set_conv_diff_supg_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

!   allocate more arrays

    if ( first .and. coefficients%i(4) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

!   coefficients

    deltat  = coefficients%r(351)
    gamma   = coefficients%r(352)
    if ( diff ) alpha   = coefficients%r(353)
    if ( self ) beta = coefficients%r(354)
    timeint = coefficients%i(401)

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg = matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    call shape_derivative ( dchi, Finv, dchidx )

!   get velocity vector

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel], layer=layer )

    work4 = reshape ( u, [ndf,ncompu] )

    u_n = matmul ( phi, work4 )

!   get mesh velocity for ALE

    if ( coefficients%i(48) == 1 ) then ! ALE

!     get mesh velocity

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, um )
      tmp = reshape ( um, [ndf,ndim] )
      umesh = matmul ( phi, tmp )

!     (un-umesh).grad operator

      do ip = 1, ninti
        ungradchi(ip,:) = matmul ( dchidx(ip,:,:), u_n(ip,1:ndim)-umesh(ip,:) )
      end do

    else

      umesh = 0

!     un.grad operator

      do ip = 1, ninti
        ungradchi(ip,:) = matmul ( dchidx(ip,:,:), u_n(ip,1:ndim) )
      end do

    end if


    if ( any ( nrhs == [0,1] ) ) then

!     get variable at time step n

      call get_sysvector ( mesh, problem, oldvectors%s(2)%p, elgrp, elem, q_n, &
        physq=[1], layer=layer )

      q_ng = matmul ( chi, q_n )

      if ( timeint == 2 ) then

!       get variable at time step n-1

        call get_sysvector ( mesh, problem, oldvectors%s(3)%p, elgrp, elem, &
          q_nm1, physq=[1], layer=layer )

        q_nm1g = matmul ( chi, q_nm1 )

      end if

    else

!     get variable at time step n

      do j = 1, nrhs

        call get_sysvector ( mesh, problem, oldvectors%s1(1)%p(j), elgrp, elem,&
          mq_n(:,j), physq=[1], layer=layer )

          mq_ng(:,j) = matmul ( chi, mq_n(:,j) )

      end do

      if ( timeint == 2 ) then

!       get variable at time step n-1

        do j = 1, nrhs

          call get_sysvector ( mesh, problem, oldvectors%s1(2)%p(j), elgrp, &
            elem, mq_nm1(:,j), physq=[1], layer=layer )

          mq_nm1g(:,j) = matmul ( chi, mq_nm1(:,j) )

        end do

      end if

    end if


!   upwinding parameter and other factors

    if ( supg ) then

      if ( htypeq == 3 ) then
        esize = sum ( detF * wg )
      else
        esize = 1
      end if

!     compute h/U and h*U

      call supg_hU ( ndim, globalshape, htype=htypeq, Uscaling=Uscalingq, &
        Uglobal=coefficients%r(357), uvec=u_n(:,1:ndim)-umesh, Finv=Finv, x=x, &
        esize=esize, hoverU=hoverU, htimesU=htimesU )

      if ( coefficients%i(412) == 1 ) then
!       Courant number dependent tau
        where ( deltat < hoverU )
          hoverU = deltat
        end where
      end if

      if ( diff ) Peh = htimesU * gamma / (2*alpha)

      tauq = supg_beta ( choice=coefficients%i(411), Peh=Peh, &
                         betaconst=coefficients%r(356) ) * hoverU / 2

    end if

!   get scalars for right-hand side

    if ( any ( nrhs == [0,1] ) ) then

      call evaluate_rhs_scalar ( mesh, problem, coefficients, oldvectors, &
        elgrp, elem, f_g, xg, u, phi )

    else

      do j = 1, nrhs

        call evaluate_rhs_scalar_multi ( mesh, problem, coefficients, &
          oldvectors, elgrp, elem, mf_g(:,j), xg, u, phi, nr=j )

      end do

    end if

    fac = 1/deltat

    if ( vector ) then

!     right-hand side

      if ( any ( timeint == [0,1] ) ) then

!       first-order, semi-implicit Euler

        if ( supg ) then ! SUPG

          if ( any ( nrhs == [0,1] ) ) then  ! assemble single RHS

            do i = 1, ndfq
              elemvec(i) = sum ( ( chi(:,i) + tauq * &
                ungradchi(:,i) ) * ( fac*gamma*q_ng + f_g ) * detF * wg )
            end do

          else  ! assemble multiple RHSs

            do j = 1, nrhs
              do i = 1, ndfq
                elemvec(ndfq*(j-1)+i) = sum ( ( chi(:,i) + tauq * &
                  ungradchi(:,i) ) * ( fac*gamma*mq_ng(:,j) + mf_g(:,j) ) * &
                  detF * wg )
              end do
            end do

          end if

        else ! Galerkin

          if ( any ( nrhs == [0,1] ) ) then ! assemble single RHS

            do i = 1, ndfq
              elemvec(i) = sum ( chi(:,i) * ( fac*gamma*q_ng + f_g ) * &
                detF * wg )
            end do

          else ! assemble multiple RHSs

            do j = 1, nrhs
              do i = 1, ndfq
                elemvec(ndfq*(j-1)+i) = sum ( chi(:,i) * &
                  ( fac*gamma*mq_ng(:,j) + mf_g(:,j) ) * detF * wg )
              end do
            end do

          end if

        end if

      else

!       second-order, semi-implicit Gear (BDF2)

        if ( supg ) then ! SUPG

          if ( any ( nrhs == [0,1] ) ) then  ! assemble single RHS

            do i = 1, ndfq
              elemvec(i) = sum ( ( chi(:,i) + &
                tauq(:) * ungradchi(:,i) ) * &
                  ( fac*gamma*2*q_ng - fac*gamma*q_nm1g/2 + f_g ) * detF * wg )
            end do

          else  ! assemble multiple RHSs

            do j = 1, nrhs
              do i = 1, ndfq
                elemvec(ndfq*(j-1)+i) = sum ( ( chi(:,i) + &
                  tauq(:) * ungradchi(:,i) ) * ( fac*gamma*2*mq_ng(:,j) - &
                  fac*gamma*mq_nm1g(:,j)/2 + mf_g(:,j) ) * detF * wg )
              end do
            end do

          end if

        else ! Galerkin

          if ( any ( nrhs == [0,1] ) ) then ! assemble single RHS

            do i = 1, ndfq
              elemvec(i) = sum ( chi(:,i) * &
                  ( fac*gamma*2*q_ng - fac*gamma*q_nm1g/2 + f_g ) * detF * wg )
            end do

          else ! assemble multiple RHSs

            do j = 1, nrhs
              do i = 1, ndfq
                elemvec(ndfq*(j-1)+i) = sum ( chi(:,i) * &
                  ( fac*gamma*2*mq_ng(:,j) - fac*gamma*mq_nm1g(:,j)/2 + &
                                                       mf_g(:,j) ) * detF * wg )
              end do
            end do

          end if

        end if

      end if

    end if

    if ( matrix ) then

!     diffusion matrix

      if ( diff ) then

        do i = 1, ndfq
          do j = i, ndfq
            do ip = 1, ninti
              work(ip) = sum ( dchidx(ip,i,:) * dchidx(ip,j,:) )
            end do
            Smat(i,j) = sum ( work * detF * wg )
            Smat(j,i) = Smat(i,j)
          end do
        end do

      end if

!     mass matrix

      if ( supg ) then ! SUPG

        do i = 1, ndfq
          do j = 1, ndfq
            Mmat(i,j) = sum ( ( chi(:,i) + tauq * ungradchi(:,i) ) * &
                        chi(:,j) * detF * wg )
          end do
        end do

      else ! Galerkin

        do i = 1, ndfq
          do j = 1, ndfq
            Mmat(i,j) = sum ( chi(:,i)  * chi(:,j) * detF * wg )
          end do
        end do

      end if

!     advection matrix

      if ( supg ) then ! SUPG

        do i = 1, ndfq
          do j = 1, ndfq
            convmat(i,j) = sum ( ( chi(:,i) + tauq * ungradchi(:,i) ) * &
                           ungradchi(:,j) * detF * wg )
          end do
        end do

      else ! Galerkin

        do i = 1, ndfq
          do j = 1, ndfq
            convmat(i,j) = sum ( chi(:,i) * ungradchi(:,j) * detF * wg )
          end do
        end do

      end if

!     fill matrix

      if ( any ( timeint == [0,1] ) ) then

!       first-order, semi-implicit Euler

        elemmat = gamma * ( fac*Mmat + convmat )

      else

!       second-order, semi-implicit Gear (BDF2)

        elemmat = gamma * ( fac*3*Mmat/2 + convmat )

      end if

!     add the diffusion and/or self-dependent source term

      if ( diff ) elemmat = elemmat + alpha * Smat
      if ( self ) elemmat = elemmat + beta * Mmat

    end if

!   unset globals, gauss, shapefuncs, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_conv_diff_supg_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(4) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( ungradchi(ninti,ndfq) )
      allocate ( work(ninti) )
      allocate ( u(ncompu*ndf), um(ndim*ndf) )
      if ( any ( nrhs == [0,1] ) ) then
        allocate ( q_n(ndfq), q_ng(ninti) )
        allocate ( q_nm1(ndfq), q_nm1g(ninti) )
        allocate ( f_g(ninti) )
      else
        allocate ( mq_n(ndfq,nrhs), mq_ng(ninti,nrhs) )
        allocate ( mq_nm1(ndfq,nrhs), mq_nm1g(ninti,nrhs) )
        allocate ( mf_g(ninti,nrhs) )
      end if
      allocate ( u_n(ninti,ncompu), tmp(ndf,ndim), work4(ndf,ncompu) )
      allocate ( Mmat(ndfq,ndfq), Smat(ndfq,ndfq) )
      allocate ( convmat(ndfq,ndfq) )
      allocate ( tauq(ninti), hoverU(ninti) )
      allocate ( htimesU(ninti) )
      allocate ( Peh(ninti) )
      allocate ( umesh(ninti,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( ungradchi )
      deallocate ( work )
      deallocate ( u, um )
      if ( any ( nrhs == [0,1] ) ) then
        deallocate ( q_n, q_ng )
        deallocate ( q_nm1, q_nm1g )
        deallocate ( f_g )
      else
        deallocate ( mq_n, mq_ng )
        deallocate ( mq_nm1, mq_nm1g )
        deallocate ( mf_g )
      end if
      deallocate ( u_n, tmp, work4 )
      deallocate ( Mmat, Smat )
      deallocate ( convmat )
      deallocate ( tauq, hoverU, htimesU )
      deallocate ( Peh )
      deallocate ( umesh )

    end subroutine deallocate_arrays

  end subroutine conv_diff_supg_elem


! compute variable in all nodes

  subroutine derive_q ( mesh, problemt, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use convection_diffusion_supg_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problemt
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    if ( first ) then

!     first element in this group

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )

      call set_globals_conv_diff_supg ( coefficients )

      allocate ( chi(nodalp,ndfq) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( q_n(ndfq) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefuncq, xrnod, chi )

    end if

    call get_sysvector ( mesh, problemt, oldvectors%s(1)%p, elgrp, elem, q_n, &
      physq=[1], layer=layer )

    elemvec = matmul ( chi, q_n )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( chi )
      deallocate ( xrnod )
      deallocate ( q_n )

    end if

  end subroutine derive_q


! Preamble for the conv_diff_supg element

  subroutine set_conv_diff_supg_elem ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors )

    use convection_diffusion_supg_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors


    if ( first .and. coefficients%i(4) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

!     set globals

      call set_globals_conv_diff_supg ( coefficients )

!     allocate arrays

      allocate ( chi(ninti,ndfq), dchi(ninti,ndfq,ndim) )
      allocate ( dchidx(ninti,ndfq,ndim) )

!     set shape functions
      call set_shape_function ( shapefuncq, xig, chi, dchi )

    else if ( coefficients%i(4) == 1 .or. coefficients%i(4) == 2 ) then

!     different Gauss rule for each element

      if ( first ) then

!       first element in this group

!       set globals

        call set_globals_conv_diff_supg ( coefficients )

      end if

!     allocate arrays

      allocate ( chi(ninti,ndfq), dchi(ninti,ndfq,ndim) )
      allocate ( dchidx(ninti,ndfq,ndim) )

!     set shape functions
      call set_shape_function ( shapefuncq, xig, chi, dchi )

    end if

  end subroutine set_conv_diff_supg_elem


! set global parameters (internal element)

  subroutine set_globals_conv_diff_supg ( coefficients )

    use convection_diffusion_supg_globals_m

    type(coefficients_t), intent(in) :: coefficients

!   check size of coefficients

    call check ( coefficients, 'set_globals_conv_diff_supg', ncoefi=450, &
      ncoefr=400 )

!   set more integer parameters

    htypeq     = get_coefficient ( coefficients, index=408, default=2 )
    Uscalingq  = get_coefficient ( coefficients, index=410, default=3 )

    if ( Uscalingq == 1 ) then
      write(*,'(/2(a/))') 'Error in set_globals_conv_diff_supg:', &
        ' Uscaling == 1 not available for convection-diffusion elements '
      stop
    end if

    if ( coefficients%i(409) > 0 ) then
      write(*,'(/4(a/))') 'Warning in set_globals_conv_diff_supg:', &
        ' Setting hlocation not available for convection-diffusion elements.', &
        ' Default value is used.', &
        ' Set coefficients%i(409)=0 to get rid of this warning.'
    end if

!   set number of degrees of freedom for q

    intpolt = coefficients%i(402)

    shapefuncq%globalshape = globalshape
    shapefuncq%interpolation = intpolt
    shapefuncq%numbering = 'regular'

    call set_ndf ( shapefuncq, 'set_globals_conv_diff_supg', ndf=ndfq )

  end subroutine set_globals_conv_diff_supg


! Unset the global parameters
! (deallocate arrays allocated in set_... )

  subroutine unset_conv_diff_supg_elem ( last, coefficients )

    use convection_diffusion_supg_globals_m

    logical, intent(in) :: last
    type(coefficients_t), intent(in) :: coefficients

!   deallocate memory

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      deallocate ( chi, dchi, dchidx )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      deallocate ( chi, dchi, dchidx )

    end if

  end subroutine unset_conv_diff_supg_elem


! evaluate the scalar RHS term (f) in the convection-diffusion equation

  subroutine evaluate_rhs_scalar ( mesh, problem, coefficients, &
    oldvectors, elgrp, elem, f, x, u, phi )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    integer, intent(in) :: elgrp, elem
    real(dp), dimension(:), intent(out) :: f, u
    real(dp), dimension(:,:), intent(in) :: x, phi

    integer :: ip, funcnr

    select case ( coefficients%i(406) )

      case(0) ! right-hand side is zero

        f = 0

      case(1) ! right-hand side is constant

        f = coefficients%r(355)

      case(2) ! right-hand side is given by function

        funcnr = coefficients%i(407)
        do ip = 1, size(f)
          f(ip) = coefficients%func1(2)%p ( funcnr, x(ip,:) )
        end do

      case(3) ! right-hand side is given by nodal point values

        call get_vector ( mesh, problem, oldvectors%v(2)%p, elgrp, elem, u, &
          layer=coefficients%i(38) )

        f = matmul ( phi, u )

      case(4) ! right-hand side is given by values in the Gauss points

        call get_elvector ( mesh, oldvectors%e(2)%p, elgrp, elem, r1=f )

      case default

        call errormsg_case_default ( 'evaluate_rhs_scalar', &
          'coefficients%i(406)', int_value=coefficients%i(406) )

    end select

  end subroutine evaluate_rhs_scalar


! evaluate the scalar RHS term (f) in the convection-diffusion equation for
! multiple RHSs

  subroutine evaluate_rhs_scalar_multi ( mesh, problem, coefficients, &
    oldvectors, elgrp, elem, f, x, u, phi, nr )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    integer, intent(in) :: elgrp, elem
    real(dp), dimension(:), intent(out) :: f, u
    real(dp), dimension(:,:), intent(in) :: x, phi
    integer, intent(in) ::  nr ! the current RHS number

    integer :: ip, funcnr

    select case ( coefficients%i(414) )

      case(0) ! right-hand side is zero

        f = 0

      case(1) ! right-hand side is constant

        if ( .not. allocated(coefficients%ra(2)%a) ) then
          write(*,'(/2(a/))') 'Error in evaluate_rhs_scalar_multi:', &
           ' coefficients%ra(2)%a not allocated'
          stop
        end if

        if ( nr > size(coefficients%ra(2)%a) ) then
          write(*,'(/2(a/))') 'Error in evaluate_rhs_scalar_multi:', &
           ' RHS number exceeds size of coefficients%ra(2)%a'
          stop
        end if

        f = coefficients%ra(2)%a(nr)

      case(2) ! right-hand side is given by function

        if ( .not. allocated(coefficients%ia(2)%a) ) then
          write(*,'(/2(a/))') 'Error in evaluate_rhs_scalar_multi:', &
           ' coefficients%ia(2)%a not allocated'
          stop
        end if

        if ( nr > size(coefficients%ia(2)%a) ) then
          write(*,'(/2(a/))') 'Error in evaluate_rhs_scalar_multi:', &
           ' RHS number exceeds size of coefficients%ia(2)%a'
          stop
        end if

        funcnr =  coefficients%ia(2)%a(nr)
        do ip = 1, size(f)
          f(ip) = coefficients%func1(2)%p ( funcnr, x(ip,:) )
        end do

      case(3) ! right-hand side is given by nodal point values

        call get_vector ( mesh, problem, oldvectors%v1(2)%p(nr), elgrp, elem,&
          u, layer=coefficients%i(38) )

        f = matmul ( phi, u )

      case(4) ! right-hand side is given by values in the Gauss points

        call get_elvector ( mesh, oldvectors%e(2)%p, elgrp, elem, r1=f, nr=nr )

      case default

        call errormsg_case_default ( 'evaluate_rhs_scalar_multi', &
          'coefficients%i(414)', int_value=coefficients%i(414) )

    end select

  end subroutine evaluate_rhs_scalar_multi

end module convection_diffusion_supg_elements_m

