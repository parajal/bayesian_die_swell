
! Copyright (C) 2025-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the evolution of equivalent plastic strain (gammap)
!

module viscoelastic_elements_gammap_m

  use tfem_elem_m
  use stokes_set_globals_m
  use devss_set_globals_m
  use viscoelastic_models_m
  use viscoelastic_elements_generic_m
  use supg_utils_m

  implicit none

contains


! Internal element routine for the equivalent plastic strain (gammap)
! semi-implicit time integratione

  subroutine gammap_supg_elem ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   select depending on integration scheme

    select case ( coefficients%i(22) )
      case(1,3) ! first-order
        call gammap_supg_elem1 ( mesh, problem, elgrp, elem, matrix, vector, &
          first, last, coefficients, oldvectors, elemmat, elemvec )
      case(2,4,5) ! second-order
        call gammap_supg_elem2 ( mesh, problem, elgrp, elem, matrix, vector, &
          first, last, coefficients, oldvectors, elemmat, elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in gammap_supg_elem:', &
        ' incorrect time integration scheme  = ', coefficients%i(22)
        stop
    end select

  end subroutine gammap_supg_elem


! Internal element routine for the equivalent plastic strain (gammap)
! First-order semi-implicit time-integration

  subroutine gammap_supg_elem1 ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip, timeint, m, isv_gp
    real(dp) :: beta, fac, deltat, esize


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1, gammap=.true. )


    if ( first ) then

!     check coefficients

      call check ( coefficients, 'gammap_supg_elem1', &
        indexarray=[48,50,90], minimum=[0,0,0], maximum=[1,1,1] )

!     in case the c/lambda terms are used, check if the modes are equal

      if ( varpar ) then
!       check coefficients for variable coefficients
        call check ( coefficients, 'gammap_supg_elem1', &
          indexarray=[69,70,81], minimum=[0,0,0], maximum=[3,3,3] )
      end if

      if ( coefficients%i(90) == 1 ) then
!       rotation reinitialization
        if ( coefficients%i(71) /= 1 ) then
          write(*,'(/3(a/))') 'Error in gammap_supg_elem1:', &
            ' Rotation reinitialization is only applicable for the ', &
            ' contravariant deformation tensor (CDT) formulation. '
          stop
        end if
!       exclude some time integration methods
        if ( coefficients%i(22) == 3 ) then
          write(*,'(/a/a,i0,a/)') 'Error in gammap_supg_elem1:', &
            ' Rotation reinitialization for time integration method = ', &
              coefficients%i(22), ' not available.'
          stop
        end if
      end if

!     set globals

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

      gpstorage = get_coefficient ( coefficients, index=100, default=1 )

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

    timeint = coefficients%i(22)

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dtheta, Finv, dthetadx )


!   get velocity vector

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel], layer=layer )

    work4 = reshape ( u, [ndf,ncompu] )

    uvecn = matmul ( phi, work4 )


!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure ( mesh, oldvectors%p(1)%p, elgrp, elem, &
        coefficients, oldvectors%s(1)%p, vemopt_gammap=vemoptn_gammap )

    end if


!   Convection operator: ALE or Eulerian frame

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity at tn

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, um )

      tmp = reshape ( um, [ndf,ndim] )

      uvecmeshn = matmul ( phi, tmp )

!     (un-ugrid).grad operator

      do ip = 1, ninti
        ungradtheta(ip,:) = &
               matmul ( dthetadx(ip,:,:), uvecn(ip,1:ndim) - uvecmeshn(ip,:) )
      end do

    else

!     Eulerian frame

      uvecmeshn = 0

!     un.grad operator

      do ip = 1, ninti
        ungradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uvecn(ip,1:ndim) )
      end do

    end if


!   get conformation tensor at previous time step

    do m = mode1, mode2

!     single mode conformation
      call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
        oldvectors, cmode=cn, mode=m )

      cng(:,:,m) = matmul ( theta, cn )

      if ( coefficients%i(90) == 1 ) then

!       rotation reinitialization in integration points

        call sqrtc_b ( cng(:,:,m) )

      end if

    end do


!   get gammap at previous time step

    isv_gp = get_coefficient ( coefficients, index=102, default=2 )

    call get_gammap ( mesh, problem, elgrp, elem, oldvectors, gp=gpn, &
      isv=isv_gp )

    gpng = matmul ( theta, gpn )

!   un.grad gpn term

    if ( timeint == 3 ) then

      ungradgpn = matmul ( ungradtheta, gpn )

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
      uvec=uvecn(:,1:ndim)-uvecmeshn, Finv=Finv, x=x, esize=esize, &
      checkzero=coefficients%i(78)==1 )

    deltat = coefficients%r(8)

    if ( coefficients%i(50) == 1 ) then
!     Courant number dependent tau
      where ( deltat < hoverU )
        hoverU = deltat
      end where
    end if

    tau = beta * hoverU / 2  ! upwinding parameter

    fac = 1._dp / deltat


!   modify material parameters

    if ( varpar ) then

!     determine the variable material parameters in mvemodel

      call evaluate_mve ( mesh, problem, elem, elgrp, coefficients, &
        oldvectors, lambda=.true., nonlin=.true., alam=.true. )

    end if


!   build matrix and vector

    if ( matrix ) then

      if ( timeint == 1 ) then

!       semi-implicit

        do i = 1, ndfc
          do j = 1, ndfc
            elemmat(i,j) = sum ( ( theta(:,i) + tau * ungradtheta(:,i) ) * &
                    ( fac * theta(:,j) +  ungradtheta(:,j) )  &
                                          * detF * wg )
          end do
        end do

      else if ( timeint == 3 ) then

!       explicit

        do i = 1, ndfc
          do j = 1, ndfc
            elemmat(i,j) = sum ( ( theta(:,i) + tau * ungradtheta(:,i) ) * &
                                      fac * theta(:,j) * detF * wg )
          end do
        end do

      end if

    end if

    if ( vector ) then

!     rhs f(cn,gpn)
      call rhs_plastic_strain ( cng, gpng, fgpng, vemoptn_gammap )

      if ( timeint == 1 ) then

!       semi-implicit

        do i = 1, ndfc
          elemvec(i) = &
            sum ( ( theta(:,i) + tau * ungradtheta(:,i) ) * &
                    ( fac * gpng + fgpng ) &
                           * detF * wg )
        end do

      else if ( timeint == 3 ) then

!       explicit

        do i = 1, ndfc
          elemvec(i) = &
            sum ( ( theta(:,i) + tau * ungradtheta(:,i) ) * &
                    ( fac * gpng - ungradgpn + fgpng ) &
                                   * detF * wg )
        end do

      end if

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

      allocate ( ungradtheta(ninti,ndfc) )
      allocate ( ungradgpn(ninti) )
      allocate ( u(ncompu*ndf), um(ndim*ndf) )

      allocate ( tmp(ndf,ndim), work2(ndfg,ncompg), work4(ndf,ncompu) )
      allocate ( tau(ninti), hoverU(ninti) )

      allocate ( uvecn(ninti,ncompu), cn(ndfc,ncomp) )
      allocate ( cng(ninti,ncomp,nmodes) )
      allocate ( uvecmeshn(ninti,ndim) )

      if ( varpar ) then
        allocate ( mvemodel(ninti) )
        mvemodel = vemodel
      end if

      if ( vemcompressible ) then
        allocate ( pr(ndfp), press(ninti) )
      end if

      allocate ( gpn(ndfc), gpng(ninti), fgpng(ninti) )

      call create_vemopt_gammap ( vemoptn_gammap, dep_J=vemcompressible, &
        np=ninti )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      integer :: ip

      deallocate ( ungradtheta )
      deallocate ( ungradgpn )
      deallocate ( u, um )

      deallocate ( tmp, work2, work4 )
      deallocate ( tau, hoverU )

      deallocate ( uvecn, cn )
      deallocate ( cng )
      deallocate ( uvecmeshn )

      if ( varpar ) then
        do ip = 1, ninti
          call delete ( mvemodel(ip) )
        end do
        deallocate ( mvemodel )
      end if

      if ( vemcompressible ) then
        deallocate ( pr, press )
      end if

      deallocate ( gpn, gpng, fgpng )

      call delete ( vemoptn_gammap )

   end subroutine deallocate_arrays

  end subroutine gammap_supg_elem1


! Internal element routine for the equivalent plastic strain (gammap)
! Second-order semi-implicit time-integration

  subroutine gammap_supg_elem2 ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m
    use limits_m, only: BLOCK_IMPROPER_ALE

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip, timeint, m, isv_gp
    real(dp) :: beta, fac, deltat, esize


    if ( first ) then

!     check coefficients

      if ( coefficients%i(48) == 1 .and. BLOCK_IMPROPER_ALE ) then
        write(*,'(/6(a/))') 'Error in gammap_supg_elem2:', &
          ' Combining this routine with a moving mesh (ALE) is not ', &
          ' allowed. The effect of the moving mesh on the change ', &
          ' in nabla operator for different times within a time step ', &
          ' has not been taken into account yet.', &
          ' Set BLOCK_IMPROPER_ALE=.false. in limits_m to bypass this block.'
        stop
      end if

      call check ( coefficients, 'gammap_supg_elem2', &
        indexarray=[48,50,90], minimum=[0,0,0], maximum=[1,1,1] )

      if ( coefficients%i(90) == 1 ) then
!       rotation reinitialization
        if ( coefficients%i(71) /= 1 ) then
          write(*,'(/3(a/))') 'Error in gammap_supg_elem1:', &
            ' Rotation reinitialization is only applicable for the ', &
            ' contravariant deformation tensor (CDT) formulation. '
          stop
        end if
!       exclude some time integration methods
        if ( any( coefficients%i(22) == [2,4] ) ) then
          write(*,'(/a/a,i0,a/)') 'Error in gammap_supg_elem2:', &
            ' Rotation reinitialization for time integration method = ', &
              coefficients%i(22), ' not available.'
          stop
        end if
      end if

!     set globals

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

      gpstorage = get_coefficient ( coefficients, index=100, default=1 )

    end if

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, gammap=.true. )

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

    call shape_derivative ( dtheta, Finv, dthetadx )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if


!   get velocity vector un

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel], posu=posu, layer=layer )

    work4 = reshape ( u, [ndf,ncompu] )

    uvecn = matmul ( phi, work4 )


!   get velocity vector unm1

    u = oldvectors%s(2)%p%u(posu)

    work4 = reshape ( u, [ndf,ncompu] )

    uvecnm1 = matmul ( phi, work4 )


!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure ( mesh, oldvectors%p(1)%p, elgrp, elem, &
        coefficients, oldvectors%s(1)%p, vemopt_gammap=vemoptn_gammap )
      call get_J_from_pressure ( mesh, oldvectors%p(1)%p, elgrp, elem, &
        coefficients, oldvectors%s(2)%p, vemopt_gammap=vemoptnm1_gammap )

    end if


!   Convection operator: ALE or Eulerian frame

    if ( coefficients%i(48) == 1 ) then

!     ALE

!     get mesh velocity at tn

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, elgrp, &
        elem, um )

      tmp = reshape ( um, [ndf,ndim] )

      uvecmeshn = matmul ( phi, tmp )

!     get mesh velocity at tn-1

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(2)%p, elgrp, &
        elem, um )

      tmp = reshape ( um, [ndf,ndim] )

      uvecmeshnm1 = matmul ( phi, tmp )

!     (un-ugrid).grad operator

      if ( any ( timeint == [2,4] ) ) then
        do ip = 1, ninti
          ungradtheta(ip,:) = &
                matmul ( dthetadx(ip,:,:), uvecn(ip,1:ndim) - uvecmeshn(ip,:) )
        end do
      end if

!     uhat.grad operator

      if ( any ( timeint == [2,5] ) ) then
        uhat = 2 * ( uvecn(:,1:ndim) - uvecmeshn ) &
                                     - ( uvecnm1(:,1:ndim) - uvecmeshnm1 )
        do ip = 1, ninti
          uhatgradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uhat(ip,:) )
        end do
      end if

!     (unm1-ugridnm1).grad operator

      if ( timeint == 4 ) then
        do ip = 1, ninti
          unm1gradtheta(ip,:) = &
            matmul ( dthetadx(ip,:,:), uvecnm1(ip,1:ndim) - uvecmeshnm1(ip,:) )
        end do
      end if

    else

!     Euler

      uvecmeshn = 0

!     un.grad operator

      if ( any ( timeint == [2,4] ) ) then
        do ip = 1, ninti
          ungradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uvecn(ip,1:ndim) )
        end do
      end if

!     uhat.grad operator

      if ( any ( timeint == [2,5] ) ) then
        uhat = 2 * uvecn(:,1:ndim) - uvecnm1(:,1:ndim)
        do ip = 1, ninti
          uhatgradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uhat(ip,:) )
        end do
      end if

!     unm1.grad operator

      if ( timeint == 4 ) then
        do ip = 1, ninti
          unm1gradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uvecnm1(ip,1:ndim) )
        end do
      end if

    end if


!   get conformation tensors

    do m = mode1, mode2

!     get conformation tensor at time step tn

!     single mode conformation
      call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
        oldvectors, cmode=cn, mode=m )

      cng(:,:,m) = matmul ( theta, cn )

!     get conformation tensor at time step tn-1

!     single mode conformation
      call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
        oldvectors, cmode=cnm1, mode=m, isv=2 )

      cnm1g(:,:,m) = matmul ( theta, cnm1 )

      if ( coefficients%i(90) == 1 ) then

!       rotation reinitialization in integration points

        call sqrtc_b ( cng(:,:,m), RT=RT )
        call sqrtc_b ( cnm1g(:,:,m) )
        call rotate_b ( cnm1g(:,:,m), RT )

      end if

    end do


!   get gammap at tn and tn-1

    isv_gp = get_coefficient ( coefficients, index=102, default=3 )

!   get gammap at tn
    call get_gammap ( mesh, problem, elgrp, elem, oldvectors, gp=gpn, &
      isv=isv_gp )

    gpng = matmul ( theta, gpn )

!   get gammap at tn-1

    call get_gammap ( mesh, problem, elgrp, elem, oldvectors, gp=gpnm1, &
      isv=isv_gp+1 )

    gpnm1g = matmul ( theta, gpnm1 )

!   un.grad gpn term

    if ( any ( timeint == [2,4] ) ) then
      ungradgpn = matmul ( ungradtheta, gpn )
    end if

!   unm1.grad gpnm1 term

    if ( timeint == 4 ) then
      unm1gradgpnm1 = matmul ( unm1gradtheta, gpnm1 )
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
      uvec=uvecn(:,1:ndim)-uvecmeshn, Finv=Finv, x=x, esize=esize, &
      checkzero=coefficients%i(78)==1 )

    deltat = coefficients%r(8)

    if ( coefficients%i(50) == 1 ) then
!     Courant number dependent tau
      where ( deltat < hoverU )
        hoverU = deltat
      end where
    end if

    tau = beta * hoverU / 2  ! upwinding parameter

    fac = 1._dp / deltat


!   build matrix and vector

    if ( matrix ) then

      select case ( timeint )

      case(2)

!       semi-implicit CN/AB2

        do i = 1, ndfc
          do j = 1, ndfc
            elemmat(i,j) = sum ( ( theta(:,i) + tau * uhatgradtheta(:,i) ) * &
                    ( fac * theta(:,j) + 0.5_dp * uhatgradtheta(:,j) )  &
                                   * detF * wg )
          end do
        end do

      case(4)

!       explicit AB2

        do i = 1, ndfc
          do j = 1, ndfc
            elemmat(i,j) = sum ( ( theta(:,i) + tau * ungradtheta(:,i) ) * &
                    ( fac * theta(:,j) )  &
                                   * detF * wg )
          end do
        end do

      case(5)

!       semi-implicit Gear

        do i = 1, ndfc
          do j = 1, ndfc
            elemmat(i,j) = sum ( ( theta(:,i) + tau * uhatgradtheta(:,i) ) * &
                    ( 1.5_dp * fac * theta(:,j) + uhatgradtheta(:,j) )  &
                                   * detF * wg )
          end do
        end do

      case default

        call errormsg_case_default ( 'gammap_supg_elem2', &
          'timeint', int_value=timeint )

      end select

    end if

    if ( vector ) then

!     rhs f(cn,gpn)
      call rhs_plastic_strain ( cng, gpng, fgpng, vemoptn_gammap )
!     rhs f(cnm1,gpnm1)
      call rhs_plastic_strain ( cnm1g, gpnm1g, fgpnm1g, vemoptnm1_gammap )

      select case ( timeint )

      case(2)

!       semi-implicit CN/AB2

        do i = 1, ndfc
          elemvec(i) = &
            sum ( ( theta(:,i) + tau * uhatgradtheta(:,i) ) * &
                    ( fac * gpng - 0.5_dp * ungradgpn + &
                      1.5_dp * fgpng - 0.5_dp * fgpnm1g )  &
                                   * detF * wg )
        end do

      case(4)

!       explicit AB2

        do i = 1, ndfc
          elemvec(i) = &
            sum ( ( theta(:,i) + tau * ungradtheta(:,i) ) * &
                    ( fac * fgpng + &
                      1.5_dp * ( - ungradgpn + fgpng ) &
                    - 0.5_dp * ( - unm1gradgpnm1 + fgpnm1g ) ) &
                                   * detF * wg )
        end do

      case(5)

!       semi-implicit Gear

        do i = 1, ndfc
          elemvec(i) = &
            sum ( ( theta(:,i) + tau * uhatgradtheta(:,i) ) * &
                    ( 2 * fac * gpng - 0.5_dp * fac * gpnm1g &
                  + 2 * fgpng - fgpnm1g )  &
                                   * detF * wg )
        end do

      case default

        call errormsg_case_default ( 'gammap_supg_elem2', &
          'timeint', int_value=timeint )

      end select

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

      allocate ( ungradtheta(ninti,ndfc) )
      allocate ( ungradgpn(ninti) )
      allocate ( unm1gradgpnm1(ninti) )
      allocate ( uhatgradtheta(ninti,ndfc) )
      allocate ( unm1gradtheta(ninti,ndfc) )
      allocate ( u(ncompu*ndf), um(ndim*ndf) )
      allocate ( tmp(ndf,ndim), work4(ndf,ncompu) )
      allocate ( posu(ncompu*ndf) )
      allocate ( tau(ninti), hoverU(ninti) )

      allocate ( uvecn(ninti,ncompu), cn(ndfc,ncomp) )
      allocate ( cng(ninti,ncomp,nmodes) )
      allocate ( RT(ninti,ncompu,ncompu) )
      allocate ( uvecnm1(ninti,ncompu), cnm1(ndfc,ncomp) )
      allocate ( cnm1g(ninti,ncomp,nmodes) )
      allocate ( uhat(ninti,ndim) )
      allocate ( uvecmeshn(ninti,ndim), uvecmeshnm1(ninti,ndim) )

      if ( vemcompressible ) then
        allocate ( pr(ndfp), press(ninti) )
      end if

      allocate ( gpn(ndfc), gpng(ninti), fgpng(ninti) )
      allocate ( gpnm1(ndfc), gpnm1g(ninti), fgpnm1g(ninti) )

      call create_vemopt_gammap ( vemoptn_gammap, dep_J=vemcompressible, &
        np=ninti )
      call create_vemopt_gammap ( vemoptnm1_gammap, dep_J=vemcompressible, &
        np=ninti )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( ungradtheta )
      deallocate ( ungradgpn )
      deallocate ( unm1gradgpnm1 )
      deallocate ( uhatgradtheta )
      deallocate ( unm1gradtheta )
      deallocate ( u, um )
      deallocate ( tmp, work4 )
      deallocate ( posu )
      deallocate ( tau, hoverU )

      deallocate ( uvecn, cn )
      deallocate ( cng )
      deallocate ( RT )
      deallocate ( uvecnm1, cnm1 )
      deallocate ( cnm1g )
      deallocate ( uhat )
      deallocate ( uvecmeshn, uvecmeshnm1 )

      if ( vemcompressible ) then
        deallocate ( pr, press )
      end if

      deallocate ( gpn, gpng, fgpng )
      deallocate ( gpnm1, gpnm1g, fgpnm1g )

      call delete ( vemoptn_gammap, vemoptnm1_gammap )

    end subroutine deallocate_arrays

  end subroutine gammap_supg_elem2


! Internal element routine for the equivalent plastic strain (gammap)
! Second-order semi-implicit time-integration

  subroutine gammap_supg_elem_implicit_2nd_order ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m
    use limits_m, only: BLOCK_IMPROPER_ALE

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip, timeint
    real(dp) :: beta, fac, deltat, esize


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1, gammap=.true. )

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'gammap_supg_elem_implicit_2nd_order', &
        indexarray=[48,50,90], minimum=[0,0,0], maximum=[1,1,1] )

      if ( coefficients%i(48) == 1 .and. coefficients%i(22) == 6 &
           .and. BLOCK_IMPROPER_ALE ) then
        write(*,'(/6(a/))') 'Error in gammap_supg_elem_implicit_2nd_order:', &
          ' Combining ALE with Crank-Nicolson not allowed, ', &
          ' since the change in the nabla operator due to the change ', &
          ' in mesh at different times within a single step ', &
          ' has not been taken into account yet.', &
          ' Set BLOCK_IMPROPER_ALE=.false. in limits_m to bypass this block.'
        stop
      end if

      if ( varpar ) then

!       check coefficients for variable coefficients
        call check ( coefficients, 'gammap_supg_elem_implicit_2nd_order', &
          indexarray=[69,70,81], minimum=[0,0,0], maximum=[3,3,3] )

        if ( coefficients%i(22) == 6 ) then
          write(*,'(/4(a/))') 'Error in gammap_supg_elem_implicit_2nd_order:', &
            ' Combining variable coefficients with Crank-Nicolson ', &
            ' not allowed, since the variable coefficients are only', &
            ' available for a single time step.'
          stop
        end if

      end if

      if ( coefficients%i(90) == 1 ) then
!       rotation reinitialization
        if ( coefficients%i(71) /= 1 ) then
          write(*,'(/3(a/))') 'Error in gammap_supg_elem_implicit_2nd_order:', &
            ' Rotation reinitialization is only applicable for the ', &
            ' contravariant deformation tensor (CDT) formulation. '
          stop
        end if
!       exclude some time integration methods
        if ( coefficients%i(22) == 6 ) then
          write(*,'(/a/a,i0,a/)') &
            'Error in gammap_supg_elem_implicit_2nd_order:', &
            ' Rotation reinitialization for time integration method = ', &
             coefficients%i(22), ' not available.'
          stop
        end if
      end if

!     set globals

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

      gpstorage = get_coefficient ( coefficients, index=100, default=1 )

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


!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure ( mesh, oldvectors%p(1)%p, elgrp, elem, &
        coefficients, oldvectors%s(1)%p, vemopt_gammap=vemoptnp1_gammap )

    end if


!   set time step

    deltat = coefficients%r(8)


!   Convection operator: ALE or Eulerian frame

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

    else

!     Eulerian frame

      uvecmeshnp1 = 0

!     unp1.grad operator

      do ip = 1, ninti
        unp1gradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uvecnp1(ip,1:ndim) )
      end do

    end if


!   get additional vectors at tn (CN only)

    if ( timeint == 6 ) then

!     get velocity vector un

      u = oldvectors%s(2)%p%u(posu)

      work4 = reshape ( u, [ndf,ncompu] )

      uvecn = matmul ( phi, work4 )

!     get relative change in volume J from pressure in system vector

      if ( vemcompressible ) then

        call get_J_from_pressure ( mesh, oldvectors%p(1)%p, elgrp, elem, &
          coefficients, oldvectors%s(2)%p, vemopt_gammap=vemoptn_gammap )

      end if

!     Convection operator: ALE or Eulerian frame

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


!   get conformation tensor at previous time steps (cn, cnm1)

    call get_c_standard


!   get gammap at previous time steps (gpn, gpnm1)

    call get_gp_standard


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
      uvec=uvecnp1(:,1:ndim)-uvecmeshnp1, Finv=Finv, x=x, esize=esize, &
      checkzero=coefficients%i(78)==1 )

    if ( coefficients%i(50) == 1 ) then
!     Courant number dependent tau
      where ( deltat < hoverU )
        hoverU = deltat
      end where
    end if

    tau = beta * hoverU / 2  ! upwinding parameter

    fac = 1._dp / deltat


!   modify material parameters

    if ( varpar ) then

!     determine the variable material parameters in mvemodel

      call evaluate_mve ( mesh, problem, elem, elgrp, coefficients, &
        oldvectors, lambda=.true., nonlin=.true., alam=.true. )

    end if


!   build matrix and vector

    if ( matrix ) then

      select case ( timeint )

      case(6)

!       semi-implicit Crank-Nicolson with relaxation prediction

        do i = 1, ndfc
          do j = 1, ndfc
            elemmat(i,j) = sum ( ( theta(:,i) + tau * unp1gradtheta(:,i) ) * &
                    ( fac * theta(:,j) + 0.5_dp * unp1gradtheta(:,j) )  &
                                   * detF * wg )
          end do
        end do

      case(7)

!       semi-implicit Gear with relaxation prediction

        do i = 1, ndfc
          do j = 1, ndfc
            elemmat(i,j) = sum ( ( theta(:,i) + tau * unp1gradtheta(:,i) ) * &
             ( 1.5_dp * fac * theta(:,j) + unp1gradtheta(:,j) ) &
                                   * detF * wg )
          end do
        end do

      case default

        write(*,'(/a/a,i0/)') 'Error in gammap_supg_elem_implicit_2nd_order:', &
        ' incorrect time integration scheme  = ', timeint
        stop

      end select

    end if

    if ( vector ) then

      select case ( timeint )

      case(6)

!       semi-implicit Crank-Nicolson with relaxation prediction

!       rhs f(chat,gphat)
        call rhs_plastic_strain ( chatg, gphatg, fgphatg, vemoptnp1_gammap )
!       rhs f(cn,gpn)
        call rhs_plastic_strain ( cng, gpng, fgpng, vemoptn_gammap )

        gprhsmodel = 0.5_dp*fgphatg + 0.5_dp*fgpng
        do i = 1, ndfc
          elemvec(i) = &
            sum ( ( theta(:,i) + tau * unp1gradtheta(:,i) ) * &
                  ( fac * gpng - 0.5_dp*ungradgpn + &
                    gprhsmodel ) * detF * wg )
        end do

      case(7)

!       semi-implicit Gear with relaxation prediction

!       rhs f(chat,gphat)
        call rhs_plastic_strain ( chatg, gphatg, fgphatg, vemoptnp1_gammap )

        do i = 1, ndfc
          elemvec(i) = &
            sum ( ( theta(:,i) + tau * unp1gradtheta(:,i) ) * &
                    ( 2 * fac * gpng - 0.5_dp * fac * gpnm1g &
                         + fgphatg ) * detF * wg )
        end do

      case default

        write(*,'(/a/a,i0/)') 'Error in gammap_supg_elem_implicit_2nd_order:', &
        ' incorrect time integration scheme  = ', timeint
        stop

      end select

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

      allocate ( unp1gradtheta(ninti,ndfc) )
      allocate ( ungradgpn(ninti) )
      allocate ( u(ncompu*ndf), um(ndim*ndf) )
      allocate ( tmp(ndf,ndim), work4(ndf,ncompu) )
      allocate ( posu(ncompu*ndf) )
      allocate ( tau(ninti), hoverU(ninti) )

      allocate ( uvecnp1(ninti,ncompu), cn(ndfc,ncomp) )
      allocate ( cng(ninti,ncomp,nmodes) )
      allocate ( uvecn(ninti,ncompu), cnm1(ndfc,ncomp) )
      allocate ( cnm1g(ninti,ncomp,nmodes) )
      allocate ( RT(ninti,ncompu,ncompu) )

      allocate ( ungradtheta(ninti,ndfc) )
      allocate ( chatg(ninti,ncomp,nmodes) )
      allocate ( fgphatg(ninti), gprhsmodel(ninti) )
      allocate ( uvecmeshn(ninti,ndim), uvecmeshnp1(ninti,ndim) )

      if ( varpar ) then
        allocate ( mvemodel(ninti) )
        mvemodel = vemodel
      end if

      if ( vemcompressible ) then
        allocate ( pr(ndfp), press(ninti) )
      end if

      allocate ( gpn(ndfc), gpng(ninti), fgpng(ninti) )
      allocate ( gpnm1(ndfc), gpnm1g(ninti), fgpnm1g(ninti) )

      call create_vemopt_gammap ( vemoptnp1_gammap, dep_J=vemcompressible, &
        np=ninti )
      if ( coefficients%i(22) == 6 ) then
        call create_vemopt_gammap ( vemoptn_gammap, dep_J=vemcompressible, &
          np=ninti )
      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      integer :: ip

      deallocate ( unp1gradtheta )
      deallocate ( ungradgpn )
      deallocate ( u, um )
      deallocate ( tmp, work4 )
      deallocate ( posu )
      deallocate ( tau, hoverU )

      deallocate ( uvecnp1, cn )
      deallocate ( cng )
      deallocate ( uvecn, cnm1 )
      deallocate ( cnm1g )
      deallocate ( RT )

      deallocate ( ungradtheta )
      deallocate ( chatg, fgphatg, gprhsmodel )
      deallocate ( uvecmeshn, uvecmeshnp1 )

      if ( varpar ) then
        do ip = 1, ninti
          call delete ( mvemodel(ip) )
        end do
        deallocate ( mvemodel )
      end if

      if ( vemcompressible ) deallocate ( pr, press )

      deallocate ( gpn, gpng, fgpng )
      deallocate ( gpnm1, gpnm1g, fgpnm1g )

      call delete ( vemoptnp1_gammap )
      if ( coefficients%i(22) == 6 ) call delete ( vemoptn_gammap )

    end subroutine deallocate_arrays


!   get cn and cnm1 for the standard case (ALE or Euler)

    subroutine get_c_standard

      integer :: m

      do m = mode1, mode2

!       get conformation tensor at time step tn

!       single mode conformation
        call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
          oldvectors, cmode=cn, mode=m )

        cng(:,:,m) = matmul ( theta, cn )

!       get conformation tensor at time step tn-1

!       single mode conformation
        call get_conformation ( mesh, oldvectors%p(2)%p, elgrp, elem, &
          oldvectors, cmode=cnm1, mode=m, isv=2 )

        cnm1g(:,:,m) = matmul ( theta, cnm1 )

        if ( coefficients%i(90) == 1 ) then

!         rotation reinitialization in integration points

          call sqrtc_b ( cng(:,:,m), RT=RT )
          call sqrtc_b ( cnm1g(:,:,m) )
          call rotate_b ( cnm1g(:,:,m), RT )

        end if

!       prediction of cn+1 (chat)

        chatg(:,:,m) = 2 * cng(:,:,m) - cnm1g(:,:,m)

      end do

    end subroutine get_c_standard


!   get gpn and gpnm1 for the standard case (ALE or Euler)

    subroutine get_gp_standard

      integer :: isv_gp

      isv_gp = get_coefficient ( coefficients, index=102, default=3 )

!     get gammap at tn

      call get_gammap ( mesh, problem, elgrp, elem, oldvectors, gp=gpn, &
        isv=isv_gp )

      gpng = matmul ( theta, gpn )

!     get gammap at tn-1

      call get_gammap ( mesh, problem, elgrp, elem, oldvectors, gp=gpnm1, &
        isv=isv_gp+1 )

      gpnm1g = matmul ( theta, gpnm1 )

!     prediction of gpn+1 (gphat)

      gphatg = 2 * gpng - gpnm1g

!     un.grad gpn term (CN scheme only)

      if ( timeint == 6 ) then
        ungradgpn = matmul ( ungradtheta, gpn )
      end if

    end subroutine get_gp_standard


  end subroutine gammap_supg_elem_implicit_2nd_order


! Internal element routine for the equivalent plastic strain (gammap)
! Implicit time-integration with Newton-Raphson iteration.
! This computes the Jacobian and the rhs.

  subroutine implicit_gammap_supg_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m
    use limits_m, only: BLOCK_IMPROPER_ALE
    use supg_utils_m
    use viscoelastic_elements_generic_m, only: rhs_plastic_strain

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, k, l, ip, timeint, n, n1, iv_supg
    integer, dimension(2) :: j12, j10, j6, jz, j2
    real(dp) :: beta, fac, deltat, esize, thetapar1, thetapar2
    real(dp) :: Kmod


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1, gammap=.true. )

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'implicit_gammap_supg_elem', &
        indexarray=[22,48,50,89,91,103], minimum=[8,0,0,0,0,0], &
        maximum=[11,1,1,1,1,1] )

      if ( coefficients%i(48) == 1 .and. any(coefficients%i(22) == [9,11]) &
           .and. BLOCK_IMPROPER_ALE ) then
        write(*,'(/6(a/))') 'Error in implicit_gammap_supg_elem:', &
          ' Combining ALE with Crank-Nicolson or theta method not allowed, ', &
          ' since the change in the nabla operator due to the change ', &
          ' in mesh at different times within a single step ', &
          ' has not been taken into account yet.', &
          ' Set BLOCK_IMPROPER_ALE=.false. in limits_m to bypass this block.'
        stop
      end if

      if ( coefficients%i(89) == 1 ) then
        if ( all(coefficients%i(22)/=[8,10]) ) then
          write(*,'(/3(a/))') 'Error in implicit_gammap_supg_elem:', &
            ' Excluding the time-derivate is only allowed if the ', &
            ' specified time integration scheme is BDF1 or BDF2.'
          stop
        end if
      end if

      if ( varpar ) then

!       check coefficients for variable coefficients
        call check ( coefficients, 'implicit_gammap_supg_elem', &
          indexarray=[69,70,81], minimum=[0,0,0], maximum=[3,3,3] )

        if ( any(coefficients%i(22) == [9,11]) ) then
          write(*,'(/4(a/))') 'Error in implicit_gammap_supg_elem:', &
            ' Combining variable coefficients with Crank-Nicolson ', &
            ' or theta-method not allowed, since the variable coefficients ', &
            ' are only available for a single time step.'
          stop
        end if

      end if

      if ( coefficients%i(50) == 1 .and. coefficients%i(91) == 1 ) then
        write(*,'(/3(a/))') 'Error in implicit_gammap_supg_elem:', &
          ' Jacobian for SUPG test function not available if the ', &
          ' tau parameter is adjusted based on the Courant number.'
        stop
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
        coefficients, oldvectors%s(1)%p, vemopt_gammap=vemoptiter_gammap )

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
          coefficients, oldvectors%s(2)%p, vemopt_gammap=vemoptn_gammap )

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


!   get citer and for CN at previous time step cn

    call get_c_standard


!   get gammap at previous time steps (gpn, gpnm1) and gpiter

    call get_gp_standard


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

!       rhs f(citer,gpiter)
        call rhs_plastic_strain ( citerg, gpiterg, fgpiterg, vemoptiter_gammap )

!       Jacobian of -rhs wrt to conformation

        do n = mode1, mode2
          n1 = n - mode1 + 1
          do l = 1, ncomp
            do k = 1, ndfc
              do i = 1, ndfc
                work12(i,k,l,n1) = &
                   - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                    vemoptiter_gammap%drhs_mm(:,l,n) * theta(:,k) * detF * wg )
              end do
            end do
          end do
        end do

!       Jacobian of material time derivative

        do j = 1, ndfc
          do i = 1, ndfc
            work11(i,j) = sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
             ( fac * theta(:,j) + unp1gradtheta(:,j) ) &
                                   * detF * wg )
          end do
        end do

!       include Jacobian of gammap

        do k = 1, ndfc
          do i = 1, ndfc
            work2(i,k) = &
                 - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                             vemoptiter_gammap%drhsdgammap * &
                                theta(:,k) * detF * wg )
          end do
        end do

        if ( coefficients%i(103) == 1 ) then

!         Jacobian of u.grad(gammap) wrt to u:  delta u.grad(gammap)

          do l = 1, ndim
            do k = 1, ndf
              do i = 1, ndfc
                work7(i,k,l) = &
                sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                         phi(:,k) * gradgpiter(:,l) * detF * wg )
              end do
            end do
          end do

        end if

        if ( coefficients%i(91) == 1 ) then

!         Jacobian of SUPG test function

          do l = 1, ndim
            do k = 1, ndf
              do i = 1, ndfc
                work6(i,k,l) = &
                     - sum ( ( beta/2 * dhoverU(:,l) * ugradtheta_supg(:,i) &
                                 + tau * dthetadx(:,i,l) ) * phi(:,k) * &
                            ( - fac * gpiterg + fac * gpng &
                                - unp1gradgpiter + fgpiterg ) * detF * wg )
              end do
            end do
          end do

        end if

        if ( vemcompressible ) then

!         include Jacobian of J

          do k = 1, ndfp
            do i = 1, ndfc
              work10(i,k) = &
                sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                          vemoptiter_gammap%drhsdJ * &
                          vemoptiter_gammap%J * psi(:,k) * detF * wg ) / Kmod
            end do
          end do

        end if

      case(9,11)

!       Crank-Nicolson (trapezoidal) or theta method

!       rhs f(citer,gpiter)
        call rhs_plastic_strain ( citerg, gpiterg, fgpiterg, vemoptiter_gammap )

!       Jacobian of -rhs wrt to conformation

        do n = mode1, mode2
          n1 = n - mode1 + 1
          do l = 1, ncomp
            do k = 1, ndfc
              do i = 1, ndfc
                work12(i,k,l,n1) = &
                      - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                        thetapar1 * vemoptiter_gammap%drhs_mm(:,l,n) * &
                          theta(:,k) * detF * wg )
              end do
            end do
          end do
        end do

!       Jacobian of material time derivative

        do j = 1, ndfc
          do i = 1, ndfc
            work11(i,j) = sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
             ( fac * theta(:,j) + thetapar1 * unp1gradtheta(:,j) ) &
                                   * detF * wg )
          end do
        end do

!       include Jacobian of gammap

        do k = 1, ndfc
          do i = 1, ndfc
            work2(i,k) = &
                 - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                             thetapar1 * vemoptiter_gammap%drhsdgammap * &
                                theta(:,k) * detF * wg )
          end do
        end do

        if ( coefficients%i(103) == 1 ) then

!         Jacobian of u.grad(gammap) wrt to u:  delta u.grad(gammap)

          do l = 1, ndim
            do k = 1, ndf
              do i = 1, ndfc
                work7(i,k,l) = &
                sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                         phi(:,k) * gradgpiter(:,l) * detF * wg )
              end do
            end do
          end do

          work7 = thetapar1 * work7

        end if

        if ( coefficients%i(91) == 1 ) then

!         Jacobian of SUPG test function

!         rhs f(cn,gpn)
          call rhs_plastic_strain ( cng, gpng, fgpng, vemoptn_gammap )

          gprhsmodel = thetapar1*fgpiterg + thetapar2*fgpng

          do l = 1, ndim
            do k = 1, ndf
              do i = 1, ndfc
                work6(i,k,l) = &
                 - sum ( ( beta/2 * dhoveru(:,l) * ugradtheta_supg(:,i) &
                             + tau * dthetadx(:,i,l) ) * phi(:,k) * &
                      ( - fac * gpiterg + fac * gpng &
                             - thetapar1*unp1gradgpiter &
                             - thetapar2*ungradgpn + &
                                  gprhsmodel ) * detF * wg )
              end do
            end do
          end do

        end if

        if ( vemcompressible ) work10= 0 ! Jacobian of j not implemented

      case(10)

!       implicit Gear

!       rhs f(citer,gpiter)
        call rhs_plastic_strain ( citerg, gpiterg, fgpiterg, vemoptiter_gammap )

!       Jacobian of -rhs wrt to conformation

        do n = mode1, mode2
          n1 = n - mode1 + 1
          do l = 1, ncomp
            do k = 1, ndfc
              do i = 1, ndfc
                work12(i,k,l,n1) = &
                   - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                    vemoptiter_gammap%drhs_mm(:,l,n) * theta(:,k) * detF * wg )
              end do
            end do
          end do
        end do

!       Jacobian of material time derivative

        do i = 1, ndfc
          do j = 1, ndfc
            work11(i,j) = sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
              ( 1.5_dp * fac * theta(:,j) + unp1gradtheta(:,j) ) &
                                   * detF * wg )
          end do
        end do

!       include Jacobian of gammap

        do k = 1, ndfc
          do i = 1, ndfc
            work2(i,k) = &
                 - sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                             vemoptiter_gammap%drhsdgammap * &
                                theta(:,k) * detF * wg )
          end do
        end do

        if ( coefficients%i(103) == 1 ) then

!         Jacobian of u.grad(gammap) wrt to u:  delta u.grad(gammap)

          do l = 1, ndim
            do k = 1, ndf
              do i = 1, ndfc
                work7(i,k,l) = &
                sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                         phi(:,k) * gradgpiter(:,l) * detF * wg )
              end do
            end do
          end do

        end if

        if ( coefficients%i(91) == 1 ) then

!         Jacobian of SUPG test function

          do l = 1, ndim
            do k = 1, ndf
              do i = 1, ndfc
                work6(i,k,l) = &
                 - sum ( ( beta/2 * dhoveru(:,l) * ugradtheta_supg(:,i) &
                             + tau * dthetadx(:,i,l) ) * phi(:,k) * &
                          ( - 1.5_dp * fac * gpiterg + &
                                   2 * fac * gpng &
                            - 0.5_dp * fac * gpnm1g &
                              - unp1gradgpiter + fgpiterg ) * detF * wg )
              end do
            end do
          end do

        end if

        if ( vemcompressible ) then

!         include Jacobian of J

          do k = 1, ndfp
            do i = 1, ndfc
              work10(i,k) = &
                sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                          vemoptiter_gammap%drhsdJ * &
                          vemoptiter_gammap%J * psi(:,k) * detF * wg ) / Kmod
            end do
          end do

        end if

      case default

        write(*,'(/a/a,i0/)') 'Error in implicit_gammap_supg_elem:', &
        ' incorrect time integration scheme  = ', timeint
        stop

      end select

!     add work11 to work2

      work2 = work2 + work11

      if ( coefficients%i(103) == 1 ) then
!       prepare contribution to velocity of delta u.grad gp
        if ( coefficients%i(91) == 1 ) then
          work6 = work6 + work7
        else
          work6 = work7
        end if
      end if

!     compute column ranges for pressure(J), velocity, zero, c,      gammap
!                               work10,      work6,    zero, work12, work2

      j10 = [ 1, ndfp ]
      j6 = [ 1, ndf*ndim ]
      jz = [ ndf*ndim+1, ndf*ncompu ]
      j12 = [ 1, ndfc*ncomp*nsubm ]
      j2 = j12(2) + [ 1, ndfc ]

      if ( coefficients%i(91) == 1 .or. coefficients%i(103) == 1 ) then
!       include Jacobian of supg test function: work6, zero, work12, work2
        j12 = j12 + jz(2)
        j2 = j2 + jz(2)
      end if

      if ( vemcompressible ) then
!       include Jacobian of J: set work10 to first block
        j6 = j6 + ndfp
        jz = jz + ndfp
        j12 = j12 + ndfp
        j2 = j2 + ndfp
      end if

!     fill elemmat for conformation

!     Jacobian wrt conformation
      elemmat(:,j12(1):j12(2)) = reshape ( work12, [ ndfc, ndfc*ncomp*nsubm ] )

!     Jacobian of gammap
      elemmat(:,j2(1):j2(2)) = work2

      if ( coefficients%i(91) == 1 .or. coefficients%i(103) == 1 ) then
!       include Jacobian of SUPG test function/deltau.grad(gammap)
        elemmat(:,j6(1):j6(2)) = reshape ( work6, [ ndfc, ndf*ndim ] )
        elemmat(:,jz(1):jz(2)) = 0
      end if

      if ( vemcompressible ) then
!       include Jacobian of J
        elemmat(:,j10(1):j10(2)) = work10
      end if

    end if

    if ( vector ) then

      if ( .not. matrix ) then

!       rhs f(citer,gpiter)
        call rhs_plastic_strain ( citerg, gpiterg, fgpiterg, &
          vemoptiter_gammap )

      end if

      select case ( timeint )

      case(8)

!       implicit Euler (BDF1)

        do i = 1, ndfc
          work1(i) = &
            sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                    ( - fac * gpiterg + fac * gpng &
                        - unp1gradgpiter + fgpiterg ) * detF * wg )
        end do

      case(9,11)

!       Crank-Nicolson (trapezoidal) or theta method

        if ( .not. ( matrix .and. coefficients%i(91) == 1 ) ) then

!         rhs f(cn,gpn)
          call rhs_plastic_strain ( cng, gpng, fgpng, vemoptn_gammap )

          gprhsmodel = thetapar1*fgpiterg + thetapar2*fgpng

        end if

        do i = 1, ndfc
          work1(i) = &
            sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                  ( - fac * gpiterg + fac * gpng &
           - thetapar1*unp1gradgpiter - thetapar2*ungradgpn + &
                              gprhsmodel ) * detF * wg )
        end do

      case(10)

!       semi-implicit Gear with relaxation prediction

        do i = 1, ndfc
          work1(i) = &
            sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                    ( - 1.5_dp * fac * gpiterg + &
                             2 * fac * gpng &
                      - 0.5_dp * fac * gpnm1g &
                        - unp1gradgpiter + fgpiterg ) * detF * wg )
        end do

      case default

        write(*,'(/a/a,i0/)') 'Error in implicit_gammap_supg_elem:', &
        ' incorrect time integration scheme  = ', timeint
        stop

      end select

      elemvec = work1

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )


!   deallocate more memory

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different gauss rule for each element

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( unp1gradtheta(ninti,ndfc), ugradtheta_supg(ninti,ndfc) )
      allocate ( u(ncompu*ndf), um(ndim*ndf) )
      allocate ( tmp(ndf,ndim), work4(ndf,ncompu) )
      allocate ( posu(ncompu*ndf) )
      allocate ( tau(ninti), hoveru(ninti) )

      allocate ( uvecnp1(ninti,ncompu), cn(ndfc,ncomp) )
      allocate ( uvec_supg(ninti,ncompu) )
      allocate ( cng(ninti,ncomp,nmodes) )
      allocate ( uvecn(ninti,ncompu) )
      allocate ( work12(ndfc,ndfc,ncomp,nsubm) )
      allocate ( work1(ndfc) )
      allocate ( work11(ndfc,ndfc) )
      if ( coefficients%i(91) == 1 .or. coefficients%i(103) == 1 ) then
        allocate ( work6(ndfc,ndf,ndim) )
      end if

      if ( coefficients%i(91) == 1 ) then
        allocate ( dhoveru(ninti,ndim) )
      end if

      allocate ( ungradtheta(ninti,ndfc) )
      allocate ( citer(ndfc,ncomp) )
      allocate ( citerg(ninti,ncomp,nmodes) )
      allocate ( uvecmeshn(ninti,ndim), uvecmeshnp1(ninti,ndim) )

      if ( varpar ) then
        allocate ( mvemodel(ninti) )
        mvemodel = vemodel
      end if

      if ( vemcompressible ) then
        allocate ( pr(ndfp), press(ninti) )
        allocate ( work10(ndfc,ndfp) )
      end if

      call create_vemopt_gammap ( vemoptiter_gammap, dep_J=vemcompressible, &
        compute_drhs_mm=.true., compute_drhsdJ=vemcompressible, &
        compute_drhsdgammap=.true., np=ninti, ncomp=ncomp, nmodes=nmodes )

      if ( any ( coefficients%i(22) == [9,11] ) ) then
        call create_vemopt_gammap ( vemoptn_gammap, dep_J=vemcompressible, &
          np=ninti )
      end if

      allocate ( gpiter(ndfc), gpiterg(ninti), fgpiterg(ninti) )
      allocate ( gpn(ndfc), gpng(ninti), fgpng(ninti) )
      allocate ( gpnm1(ndfc), gpnm1g(ninti) )
      allocate ( gprhsmodel(ninti) )
      allocate ( work2(ndfc,ndfc) )
      allocate ( unp1gradgpiter(ninti) )
      allocate ( ungradgpn(ninti) )
      if ( coefficients%i(103) == 1 ) then
        allocate ( gradgpiter(ninti,ndim) )
        allocate ( work7(ndfc,ndf,ndim) )
      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      integer :: ip

      deallocate ( unp1gradtheta, ugradtheta_supg )
      deallocate ( u, um )
      deallocate ( tmp, work4 )
      deallocate ( posu )
      deallocate ( tau, hoveru )

      deallocate ( uvecnp1, cn )
      deallocate ( uvec_supg )
      deallocate ( cng )
      deallocate ( uvecn )
      deallocate ( work12 )
      deallocate ( work1 )
      deallocate ( work11 )
      if ( coefficients%i(91) == 1 .or. coefficients%i(103) == 1 ) then
        deallocate ( work6 )
      end if

      if ( coefficients%i(91) == 1 ) then
        deallocate ( dhoveru )
      end if

      deallocate ( ungradtheta )
      deallocate ( citer )
      deallocate ( citerg )
      deallocate ( uvecmeshn, uvecmeshnp1 )

      if ( varpar ) then
        do ip = 1, ninti
          call delete ( mvemodel(ip) )
        end do
        deallocate ( mvemodel )
      end if

      if ( vemcompressible ) then
        deallocate ( pr, press )
        deallocate ( work10 )
      end if

      call delete ( vemoptiter_gammap )

      if ( any ( coefficients%i(22) == [9,11] ) ) call delete ( vemoptn_gammap )

      deallocate ( gpiter, gpiterg, fgpiterg )
      deallocate ( gpn, gpng, fgpng )
      deallocate ( gpnm1, gpnm1g )
      deallocate ( gprhsmodel )
      deallocate ( work2 )
      deallocate ( unp1gradgpiter )
      deallocate ( ungradgpn )
      if ( coefficients%i(103) == 1 ) then
        deallocate ( gradgpiter )
        deallocate ( work7 )
      end if

    end subroutine deallocate_arrays


!   get citer and for CN at previous time step cn

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

!     get conformation tensor at time step tn

      if ( any ( timeint == [9,11] ) ) then

        isv_c = get_coefficient ( coefficients, index=88, default=1 )

        do m = mode1, mode2

!         single mode conformation
          call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
            cmode=cn, mode=m, isv=isv_c )

          cng(:,:,m) = matmul ( theta, cn )

          if ( coefficients%i(90) == 1 ) then

!           rotation reinitialization in integration points

            call sqrtc_b ( cng(:,:,m) )

          end if

        end do

      end if

    end subroutine get_c_standard


!   get gpn and gpnm1 for the standard case (ALE or Euler) and gpiter

    subroutine get_gp_standard

      integer :: isv_gp, i

      if ( coefficients%i(89) == 1 ) then

!       no time-derivative included

        gpng = 0

      else

!       get equivalent plastic strain at time step tn

        isv_gp = get_coefficient ( coefficients, index=102, default=2 )

        call get_gammap ( mesh, problem, elgrp, elem, oldvectors, gp=gpn, &
          isv=isv_gp )

        gpng = matmul ( theta, gpn )

      end if

      if ( timeint == 10 ) then ! BDF2

        if ( coefficients%i(89) == 1 ) then

!         no time-derivative included

          gpnm1g = 0

        else

!         get equivalent plastic strain at time step tn-1

          call get_gammap ( mesh, problem, elgrp, elem, oldvectors, gp=gpnm1, &
            isv=isv_gp+1 )

          gpnm1g = matmul ( theta, gpnm1 )

        end if

      end if

!     get iteration of equivalent plastic strain

      isv_gp = get_coefficient ( coefficients, index=101, default=1 )

      call get_gammap ( mesh, problem, elgrp, elem, oldvectors, gp=gpiter, &
        isv=isv_gp )

      gpiterg = matmul ( theta, gpiter )

      unp1gradgpiter = matmul ( unp1gradtheta, gpiter )

!     un.grad gammapn term (cn an theta method only)

      if ( any ( timeint == [9,11] ) ) then
        ungradgpn = matmul ( ungradtheta, gpn )
      end if

      if ( coefficients%i(103) == 1 ) then

!       grad gpiter term

        do i = 1, ndim
          gradgpiter(:,i) = matmul ( dthetadx(:,:,i), gpiter )
        end do

      end if

    end subroutine get_gp_standard

  end subroutine implicit_gammap_supg_elem


! Internal element routine for the equivalent plastic strain (gammap)
! Matrix of the time-derivative only.

  subroutine implicit_gammap_timederiv_supg_elem ( mesh, problem, elgrp, elem, &
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


    integer :: i, j, ip, iv_supg
    real(dp) :: beta, esize


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1, singlemode=.true., &
      gammap=.true. )

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'implicit_gammap_timederiv_supg_elem', &
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
          elemmat(i,j) = sum ( ( theta(:,i) + tau * ugradtheta_supg(:,i) ) * &
                                        theta(:,j) * detF * wg )
        end do
      end do

    end if

    if ( vector ) then

      write(*,'(/a/a,i0/)') 'Error in implicit_gammap_timederiv_supg_elem:', &
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

      allocate ( uvecmeshnp1(ninti,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( ugradtheta_supg )
      deallocate ( u, um )
      deallocate ( tmp, work4 )
      deallocate ( posu )
      deallocate ( tau, hoverU )

      deallocate ( uvec_supg )

      deallocate ( uvecmeshnp1 )

    end subroutine deallocate_arrays

  end subroutine implicit_gammap_timederiv_supg_elem


! compute conformation tensor component in all nodes (standard)

  subroutine deriv_gammap ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    integer :: isv_gp


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp, maxvel3D=1 )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic ( coefficients, maxvar=1, gammap=.true. )

      gpstorage = get_coefficient ( coefficients, index=100, default=2 )

      allocate ( theta(nodalp,ndfc) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( gpnod(ndfc) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

!     shape function in nodal points

      call set_shape_function ( shapefuncc, xrnod, theta )

    end if

    isv_gp = get_coefficient ( coefficients, index=102, default=1 )

    call get_gammap ( mesh, problem, elgrp, elem, oldvectors, gp=gpnod, &
      isv=isv_gp )

    elemvec = matmul ( theta, gpnod)

    elemwts = 1

    if ( last ) then

!     last element in this group

      call delete ( vemodel )

      deallocate ( theta )
      deallocate ( xrnod )
      deallocate ( gpnod )

    end if

  end subroutine deriv_gammap

end module viscoelastic_elements_gammap_m
