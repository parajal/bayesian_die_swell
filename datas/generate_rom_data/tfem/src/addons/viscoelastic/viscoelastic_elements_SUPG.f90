! Copyright (C) 2005-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the constitutive equation based on SUPG
!

module viscoelastic_elements_SUPG_m

  use viscoelastic_elements_generic_m
  use supg_utils_m

  implicit none

  save

contains


! Internal element routine for the constitutive equations: SUPG, G-method

  subroutine ce_supg_elem ( mesh, problem, elgrp, elem, matrix, vector, &
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
        call ce_supg_elem1 ( mesh, problem, elgrp, elem, matrix, vector, &
          first, last, coefficients, oldvectors, elemmat, elemvec )
      case(2,4,5) ! second-order
        call ce_supg_elem2 ( mesh, problem, elgrp, elem, matrix, vector, &
          first, last, coefficients, oldvectors, elemmat, elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in ce_supg_elem:', &
        ' incorrect time integration scheme  = ', coefficients%i(22)
        stop
    end select

  end subroutine ce_supg_elem


! Internal element routine for the constitutive equations: SUPG, G-method
! First-order semi-implicit time-integration

  subroutine ce_supg_elem1 ( mesh, problem, elgrp, elem, matrix, vector, &
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

    call set_devssg_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1, gammap=.true. )


    if ( first ) then

!     check coefficients

      call check ( coefficients, 'ce_supg_elem1', &
        indexarray=[48,50,90], minimum=[0,0,0], maximum=[1,1,1] )

!     in case the c/lambda terms are used, check if the modes are equal

      if ( coefficients%i(57) == 1 .and. &
           coefficients%i(26) /= coefficients%i(27) ) then
        write(*,'(/4(a/))') 'Error in ce_supg_elem1:', &
          ' Additional c/lambda terms (coefficients%i(57)=1) requires ', &
          ' mode1=mode2 (coefficients%i(26)=coefficients%i(27) denoting ', &
          ' the current mode'
        stop
      end if

      if ( varpar ) then
!       check coefficients for variable coefficients
        call check ( coefficients, 'ce_supg_elem1', &
          indexarray=[69,70,81], minimum=[0,0,0], maximum=[3,3,3] )
      end if

      if ( coefficients%i(90) == 1 ) then
!       rotation reinitialization
        if ( coefficients%i(71) /= 1 ) then
          write(*,'(/3(a/))') 'Error in ce_supg_elem1:', &
            ' Rotation reinitialization is only applicable for the ', &
            ' contravariant deformation tensor (CDT) formulation. '
          stop
        end if
!       exclude some time integration methods
        if ( coefficients%i(22) == 3 ) then
          write(*,'(/a/a,i0,a/)') 'Error in ce_supg_elem1:', &
            ' Rotation reinitialization for time integration method = ', &
              coefficients%i(22), ' not available.'
          stop
        end if
      end if

!     set globals

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

      if ( vemgammap ) then
        gpstorage = get_coefficient ( coefficients, index=100, default=1 )
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


!   start of the element

    timeint = coefficients%i(22)

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x(1:ndf,:), phi, xg )
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( coefficients%i(61) == 1 ) then
      call shape_derivative ( dphi, Finv, dphidx )
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
        coefficients, oldvectors%s(1)%p, vemoptn )

    end if


!   get velocity gradient vector

    if ( coefficients%i(61) == 0 ) then

!     Projected gradient

      select case ( coefficients%i(60) )
      case(0)
!       DEVSS: in solution vector
        call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
          elgrp, elem, g, physq=[physqgrad], layer=layer )
      case(1)
!       separate vector
        call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(3)%p, &
          elgrp, elem, g, layer=layer )
      case default
        call errormsg_case_default ( 'ce_supg_elem1', &
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
          gvecn(ip,:) = reshape ( transpose ( gradu(ip,:,:) ), [ncompu*ncompu] )
        end do

      else

        gvecn(:,1:ncompg) = work10

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
          gvecn(ip,:) = reshape ( transpose ( gradu(ip,:,:) ), [ncompu*ncompu] )
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
      call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
        cmode=cn, mode=m )

      cng(:,:,m) = matmul ( theta, cn )

      if ( coefficients%i(90) == 1 ) then

!       rotation reinitialization in integration points

        call sqrtc_b ( cng(:,:,m) )

      end if

!     un.grad cn term

      if ( timeint == 3 ) then

        ungradcn(:,:,m) = matmul ( ungradtheta, cn )

      end if

    end do


!   get gammap at previous time step

    if ( vemgammap ) then

      isv_gp = get_coefficient ( coefficients, index=102, default=2 )

      call get_gammap ( mesh, oldvectors%p(3)%p, elgrp, elem, oldvectors, &
        gp=gpn, isv=isv_gp )

      gpng = matmul ( theta, gpn )

      vemoptn%gammap = gpng

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

      if ( coverl ) then
!       add c/lambda to both sides of the constitutive equation
        if ( coefficients%i(58) > 0 ) then ! add c/non-linear parameter
          lambda_g_mode1 = &
              [ ( mvemodel(ip)%nonlin(coefficients%i(58),mode1), ip=1,ninti ) ]
        else ! add c/lambda
          lambda_g_mode1 = [ ( mvemodel(ip)%lambda(mode1), ip=1,ninti ) ]
        end if
      end if

    else

      if ( coverl ) then
!       add c/lambda to both sides of the constitutive equation
        if ( coefficients%i(58) > 0 ) then ! add c/non-linear parameter
          lambda_g_mode1 = vemodel%nonlin(coefficients%i(58),mode1)
        else ! add c/lambda
          lambda_g_mode1 = vemodel%lambda(mode1)
        end if
      end if

    end if

    if ( coverl ) then
      fac2_g = 1._dp / lambda_g_mode1
    else
      fac2_g = 0
    end if


!   build matrix and vector

    if ( matrix ) then

      if ( timeint == 1 ) then

!       semi-implicit

        do i = 1, ndfc
          do j = 1, ndfc
            elemmat(i,j) = sum ( ( theta(:,i) + tau * ungradtheta(:,i) ) * &
                    ( ( fac + fac2_g ) * theta(:,j) +  ungradtheta(:,j) )  &
                                          * detF * wg )
          end do
        end do

      else if ( timeint == 3 ) then

!       explicit

        do i = 1, ndfc
          do j = 1, ndfc
            elemmat(i,j) = sum ( ( theta(:,i) + tau * ungradtheta(:,i) ) * &
                    ( ( fac + fac2_g ) * theta(:,j) )  &
                                          * detF * wg )
          end do
        end do

      end if

    end if

    if ( vector ) then

!     viscoelastic rhs

      call rhs_viscoelastic ( gvecn, cng, fng, vemopt=vemoptn )

      if ( coverl ) then
!       add c/lambda to both sides of the constitutive equation
        do ip = 1, ninti
          fng(ip,:,mode1) = fng(ip,:,mode1) + cng(ip,:,mode1) / &
                                                      lambda_g_mode1(ip)
        end do
      end if

      if ( coorsys == 1 .and. vel3D == 1 ) then

!       cylindrical, axisymmetric with 3D velocities: add terms due to u.nabla c

        call add_ugrad_terms ( xg, uvecn, cng, fng )

      end if

      if ( timeint == 1 ) then

!       semi-implicit

        do m = mode1, mode2
          do j = 1, ncomp
            do i = 1, ndfc
              work6(i,j,m-mode1+1) = &
                sum ( ( theta(:,i) + tau * ungradtheta(:,i) ) * &
                        ( fac * cng(:,j,m) + fng(:,j,m) )  &
                                       * detF * wg )
            end do
          end do
        end do

      else if ( timeint == 3 ) then

!       explicit

        do m = mode1, mode2
          do j = 1, ncomp
            do i = 1, ndfc
              work6(i,j,m-mode1+1) = &
                sum ( ( theta(:,i) + tau * ungradtheta(:,i) ) * &
                        ( fac * cng(:,j,m) - ungradcn(:,j,m) + fng(:,j,m) )  &
                                       * detF * wg )
            end do
          end do
        end do

      end if

      elemvec = reshape ( work6, [ ndfc * ncomp * ( mode2 - mode1 + 1 ) ] )

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

      allocate ( ungradtheta(ninti,ndfc) )
      allocate ( ungradcn(ninti,ncomp,nmodes) )
      allocate ( u(ncompu*ndf), um(ndim*ndf), g(ncompg*ndfg) )
      allocate ( gradu(ninti,ncompu,ncompu) )

      allocate ( tmp(ndf,ndim), work2(ndfg,ncompg), work4(ndf,ncompu) )
      allocate ( work10(ninti,ncompg) )
      allocate ( tau(ninti), hoverU(ninti) )
      allocate ( facv(ninti) )

      allocate ( uvecn(ninti,ncompu), cn(ndfc,ncomp) )
      allocate ( cng(ninti,ncomp,nmodes), fng(ninti,ncomp,nmodes) )
      allocate ( work6(ndfc,ncomp,mode2-mode1+1) )
      allocate ( uvecmeshn(ninti,ndim) )

      allocate ( lambda_g_mode1(ninti), fac2_g(ninti) )

      if ( coorsys == 1 .and. vel3D == 0 ) then
        allocate ( gvecn(ninti,ndim**2+1) )
      else if ( vel3D == 1 ) then
        allocate ( gvecn(ninti,ncompu**2) )
      else
        allocate ( gvecn(ninti,ndim**2) )
      end if

      if ( varpar ) then
        allocate ( mvemodel(ninti) )
        mvemodel = vemodel
      end if

      if ( vemcompressible ) then
        allocate ( pr(ndfp), press(ninti) )
      end if

      if ( vemgammap ) then
        allocate ( gpn(ndfc), gpng(ninti) )
      end if

      call create_vemopt ( vemoptn, dep_J=vemcompressible, &
        dep_gammap=vemgammap, np=ninti )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      integer :: ip

      deallocate ( ungradtheta )
      deallocate ( ungradcn )
      deallocate ( u, um, g )
      deallocate ( gradu )

      deallocate ( tmp, work2, work4 )
      deallocate ( work10 )
      deallocate ( tau, hoverU )
      deallocate ( facv )

      deallocate ( uvecn, cn )
      deallocate ( cng, fng )
      deallocate ( work6 )
      deallocate ( uvecmeshn )

      deallocate ( lambda_g_mode1, fac2_g )

      deallocate ( gvecn )

      if ( varpar ) then
        do ip = 1, ninti
          call delete ( mvemodel(ip) )
        end do
        deallocate ( mvemodel )
      end if

      if ( vemcompressible ) then
        deallocate ( pr, press )
      end if

      if ( vemgammap ) then
        deallocate ( gpn, gpng )
      end if

      call delete ( vemoptn )

   end subroutine deallocate_arrays

  end subroutine ce_supg_elem1


! Internal element routine for the constitutive equations: SUPG, G-method
! Second-order semi-implicit time-integration

  subroutine ce_supg_elem2 ( mesh, problem, elgrp, elem, matrix, vector, &
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
    real(dp) :: beta, fac, fac2, deltat, lambda_mode1, esize


    if ( first ) then

!     check coefficients

      if ( coefficients%i(48) == 1 .and. BLOCK_IMPROPER_ALE ) then
        write(*,'(/6(a/))') 'Error in ce_supg_elem2:', &
          ' Combining this routine with a moving mesh (ALE) is not ', &
          ' allowed. The effect of the moving mesh on the change ', &
          ' in nabla operator for different times within a time step ', &
          ' has not been taken into account yet.', &
          ' Set BLOCK_IMPROPER_ALE=.false. in limits_m to bypass this block.'
        stop
      end if

      call check ( coefficients, 'ce_supg_elem2', &
        indexarray=[48,50,90], minimum=[0,0,0], maximum=[1,1,1] )

      if ( coefficients%i(57) == 1 ) then
        if ( coefficients%i(22) /= 5 ) then
          write(*,'(/4(a/))') 'Error in ce_supg_elem2:', &
            'In this routine, adding of c/lambda terms (coefficients%i(57)=1)',&
            ' has only been implemented for second-order, semi-implicit', &
            ' Gear/Karniadakis (coefficients%i(22)=5).'
          stop
        end if
        if ( coefficients%i(26) /= coefficients%i(27) )  then
          write(*,'(/4(a/))') 'Error in ce_supg_elem2:', &
            ' Additional c/lambda terms (coefficients%i(57)=1) requires ', &
            ' mode1=mode2 (coefficients%i(26)=coefficients%i(27) denoting ', &
            ' the current mode'
          stop
        end if
      end if

      if ( coefficients%i(90) == 1 ) then
!       rotation reinitialization
        if ( coefficients%i(71) /= 1 ) then
          write(*,'(/3(a/))') 'Error in ce_supg_elem1:', &
            ' Rotation reinitialization is only applicable for the ', &
            ' contravariant deformation tensor (CDT) formulation. '
          stop
        end if
!       exclude some time integration methods
        if ( any( coefficients%i(22) == [2,4] ) ) then
          write(*,'(/a/a,i0,a/)') 'Error in ce_supg_elem2:', &
            ' Rotation reinitialization for time integration method = ', &
              coefficients%i(22), ' not available.'
          stop
        end if
      end if

!     set globals

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

      if ( vemgammap ) then
        gpstorage = get_coefficient ( coefficients, index=100, default=1 )
      end if

    end if

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_devssg_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

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

    if ( coefficients%i(61) == 1 ) then
      call shape_derivative ( dphi, Finv, dphidx )
    end if

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
        coefficients, oldvectors%s(1)%p, vemoptn )
      call get_J_from_pressure ( mesh, oldvectors%p(1)%p, elgrp, elem, &
        coefficients, oldvectors%s(2)%p, vemoptnm1 )

    end if


!   get velocity gradient vector at tn and tn-1

    if ( coefficients%i(61) == 0 ) then

!     Projected gradient at tn

      select case ( coefficients%i(60) )
      case(0)
!       DEVSS: in solution vector
        call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
          elgrp, elem, g, physq=[physqgrad], posu=posg, layer=layer )
      case(1)
!       separate vector
        call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(3)%p, &
          elgrp, elem, g, posu=posg, layer=layer )
      case default
        call errormsg_case_default ( 'ce_supg_elem2', &
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
          gvecn(ip,:) = reshape ( transpose ( gradu(ip,:,:) ), [ncompu*ncompu] )
        end do

      else

        gvecn(:,1:ncompg) = work10

      end if

!     Projected gradient at tn-1

      select case ( coefficients%i(60) )
      case(0)
!       DEVSS: in solution vector
        g = oldvectors%s(2)%p%u(posg)
      case(1)
!       separate vector
        g = oldvectors%v(4)%p%u(posg)
      case default
        call errormsg_case_default ( 'ce_supg_elem2', &
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
          gvecnm1(ip,:) = &
                      reshape ( transpose ( gradu(ip,:,:) ), [ncompu*ncompu] )
        end do

      else

        gvecnm1(:,1:ncompg) = work10

      end if

    else if ( coefficients%i(61) == 1 ) then

!     direct velocity gradient at tn

      u = oldvectors%s(1)%p%u(posu)

      work4 = reshape ( u, [ndf,ncompu] )

      do j = 1, ndim
        gradu(:,:,j) = matmul ( dphidx(:,:,j), work4 )
      end do

      if ( vel3D == 1 ) then
        gradu(:,:,3) = 0
        do ip = 1, ninti
          gvecn(ip,:) = reshape ( transpose ( gradu(ip,:,:) ), [ncompu*ncompu] )
        end do
      else
        do ip = 1, ninti
          gvecn(ip,1:ndim**2) = &
                       reshape ( transpose ( gradu(ip,:,:) ), [ndim*ndim] )
        end do
      end if

!     direct velocity gradient at tn-1

      u = oldvectors%s(2)%p%u(posu)

      work4 = reshape ( u, [ndf,ncompu] )

      do j = 1, ndim
        gradu(:,:,j) = matmul ( dphidx(:,:,j), work4 )
      end do

      if ( vel3D == 1 ) then
        gradu(:,:,3) = 0
        do ip = 1, ninti
          gvecnm1(ip,:) = &
                      reshape ( transpose ( gradu(ip,:,:) ), [ncompu*ncompu] )
        end do
      else
        do ip = 1, ninti
          gvecnm1(ip,1:ndim**2) = &
                       reshape ( transpose ( gradu(ip,:,:) ), [ndim*ndim] )
        end do
      end if

    end if

    if ( coorsys == 1 .and. vel3D == 0 ) then
      gvecn(:,5) = uvecn(:,2) / xg(:,2)   ! u_r / r
      gvecnm1(:,5) = uvecnm1(:,2) / xg(:,2)
    else if ( coorsys == 1 .and. vel3D == 1 ) then
      gvecn(:,6) = - uvecn(:,3) / xg(:,2) ! - u_theta / r
      gvecnm1(:,6) = - uvecnm1(:,3) / xg(:,2)
      gvecn(:,9) = uvecn(:,2) / xg(:,2)   ! u_r / r
      gvecnm1(:,9) = uvecnm1(:,2) / xg(:,2)
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
      call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
        cmode=cn, mode=m )

      cng(:,:,m) = matmul ( theta, cn )

!     get conformation tensor at time step tn-1

!     single mode conformation
      call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
        cmode=cnm1, mode=m, isv=2 )

      cnm1g(:,:,m) = matmul ( theta, cnm1 )

      if ( coefficients%i(90) == 1 ) then

!       rotation reinitialization in integration points

        call sqrtc_b ( cng(:,:,m), RT=RT )
        call sqrtc_b ( cnm1g(:,:,m) )
        call rotate_b ( cnm1g(:,:,m), RT )

      end if

!     un.grad cn term

      if ( any ( timeint == [2,4] ) ) then
        ungradcn(:,:,m) = matmul ( ungradtheta, cn )
      end if

!     unm1.grad cnm1 term

      if ( timeint == 4 ) then
        unm1gradcnm1(:,:,m) = matmul ( unm1gradtheta, cnm1 )
      end if

    end do


!   get gammap at tn and tn-1

    if ( vemgammap ) then

      isv_gp = get_coefficient ( coefficients, index=102, default=3 )

!     get gammap at tn
      call get_gammap ( mesh, oldvectors%p(3)%p, elgrp, elem, oldvectors, &
        gp=gpn, isv=isv_gp )

      gpng = matmul ( theta, gpn )

      vemoptn%gammap = gpng

!     get gammap at tn-1

      call get_gammap ( mesh, oldvectors%p(3)%p, elgrp, elem, oldvectors, &
        gp=gpnm1, isv=isv_gp+1 )

      gpnm1g = matmul ( theta, gpnm1 )

      vemoptnm1%gammap = gpnm1g

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

    if ( coverl ) then
!     add c/lambda to both sides of the constitutive equation
      if ( coefficients%i(58) > 0 ) then
        lambda_mode1 = vemodel%nonlin(coefficients%i(58),mode1)
      else
        lambda_mode1 = vemodel%lambda(mode1)
      end if
      fac2 = 1._dp / lambda_mode1
    else
      fac2 = 0
    end if


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
                    ( ( 1.5_dp * fac + fac2 ) * theta(:,j) &
                                              + uhatgradtheta(:,j) )  &
                                   * detF * wg )
          end do
        end do

      case default

        call errormsg_case_default ( 'ce_supg_elem2', &
          'timeint', int_value=timeint )

      end select

    end if

    if ( vector ) then

!     viscoelastic rhs

      call rhs_viscoelastic ( gvecn, cng, fng, vemopt=vemoptn )       ! tn
      call rhs_viscoelastic ( gvecnm1, cnm1g, fnm1g, vemopt=vemoptnm1 ) ! tn-1

      if ( coverl ) then
!       add c/lambda to both sides of the constitutive equation
        fng(:,:,mode1) = fng(:,:,mode1) + cng(:,:,mode1) / lambda_mode1
        fnm1g(:,:,mode1) = &
                    fnm1g(:,:,mode1) + cnm1g(:,:,mode1) / lambda_mode1
      end if

      if ( coorsys == 1 .and. vel3D == 1 ) then

!       cylindrical, axisymmetric with 3D velocities: add terms due to u.nabla c

        call add_ugrad_terms ( xg, uvecn, cng, fng )
        call add_ugrad_terms ( xg, uvecnm1, cnm1g, fnm1g )

      end if

      select case ( timeint )

      case(2)

!       semi-implicit CN/AB2

        do m = mode1, mode2
          do j = 1, ncomp
            do i = 1, ndfc
              work6(i,j,m-mode1+1) = &
                sum ( ( theta(:,i) + tau * uhatgradtheta(:,i) ) * &
                        ( fac * cng(:,j,m) - 0.5_dp * ungradcn(:,j,m) + &
                          1.5_dp * fng(:,j,m) - 0.5_dp * fnm1g(:,j,m) )  &
                                       * detF * wg )
            end do
          end do
        end do

      case(4)

!       explicit AB2

        do m = mode1, mode2
          do j = 1, ncomp
            do i = 1, ndfc
              work6(i,j,m-mode1+1) = &
                sum ( ( theta(:,i) + tau * ungradtheta(:,i) ) * &
                        ( fac * cng(:,j,m) + &
                          1.5_dp * ( - ungradcn(:,j,m) + fng(:,j,m) ) &
                        - 0.5_dp * ( - unm1gradcnm1(:,j,m) + fnm1g(:,j,m) ) ) &
                                       * detF * wg )
            end do
          end do
        end do

      case(5)

!       semi-implicit Gear

        do m = mode1, mode2
          do j = 1, ncomp
            do i = 1, ndfc
              work6(i,j,m-mode1+1) = &
                sum ( ( theta(:,i) + tau * uhatgradtheta(:,i) ) * &
                        ( 2 * fac * cng(:,j,m) - 0.5_dp * fac * cnm1g(:,j,m) &
                      + 2 * fng(:,j,m) - fnm1g(:,j,m) )  &
                                       * detF * wg )
            end do
          end do
        end do

      case default

        call errormsg_case_default ( 'ce_supg_elem2', &
          'timeint', int_value=timeint )

      end select

      elemvec = reshape ( work6, [ ndfc * ncomp * ( mode2 - mode1 + 1 ) ] )

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

      allocate ( ungradtheta(ninti,ndfc) )
      allocate ( ungradcn(ninti,ncomp,nmodes) )
      allocate ( unm1gradcnm1(ninti,ncomp,nmodes) )
      allocate ( uhatgradtheta(ninti,ndfc) )
      allocate ( unm1gradtheta(ninti,ndfc) )
      allocate ( u(ncompu*ndf), um(ndim*ndf), g(ncompg*ndfg) )
      allocate ( gradu(ninti,ncompu,ncompu) )
      allocate ( tmp(ndf,ndim), work2(ndfg,ncompg), work4(ndf,ncompu) )
      allocate ( work10(ninti,ncompg) )
      allocate ( posu(ncompu*ndf), posg(ncompg*ndfg) )
      allocate ( tau(ninti), hoverU(ninti) )
      allocate ( facv(ninti) )

      allocate ( uvecn(ninti,ncompu), cn(ndfc,ncomp) )
      allocate ( cng(ninti,ncomp,nmodes), fng(ninti,ncomp,nmodes) )
      allocate ( RT(ninti,ncompu,ncompu) )
      allocate ( uvecnm1(ninti,ncompu), cnm1(ndfc,ncomp) )
      allocate ( cnm1g(ninti,ncomp,nmodes), fnm1g(ninti,ncomp,nmodes) )
      allocate ( uhat(ninti,ndim) )
      allocate ( work6(ndfc,ncomp,mode2-mode1+1) )
      allocate ( uvecmeshn(ninti,ndim), uvecmeshnm1(ninti,ndim) )

      if ( coorsys == 1 .and. vel3D == 0 ) then
        allocate ( gvecn(ninti,ndim**2+1), gvecnm1(ninti,ndim**2+1) )
      else if ( vel3D == 1 ) then
        allocate ( gvecn(ninti,ncompu**2), gvecnm1(ninti,ncompu**2) )
      else
        allocate ( gvecn(ninti,ndim**2), gvecnm1(ninti,ndim**2) )
      end if

      if ( vemcompressible ) then
        allocate ( pr(ndfp), press(ninti) )
      end if

      if ( vemgammap ) then
        allocate ( gpn(ndfc), gpng(ninti) )
        allocate ( gpnm1(ndfc), gpnm1g(ninti) )
      end if

      call create_vemopt ( vemoptn, dep_J=vemcompressible, &
        dep_gammap=vemgammap, np=ninti )
      call create_vemopt ( vemoptnm1, dep_J=vemcompressible, &
        dep_gammap=vemgammap, np=ninti )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( ungradtheta )
      deallocate ( ungradcn )
      deallocate ( unm1gradcnm1 )
      deallocate ( uhatgradtheta )
      deallocate ( unm1gradtheta )
      deallocate ( u, um, g )
      deallocate ( gradu )
      deallocate ( tmp, work2, work4 )
      deallocate ( work10 )
      deallocate ( posu, posg )
      deallocate ( tau, hoverU )
      deallocate ( facv )

      deallocate ( uvecn, cn )
      deallocate ( cng, fng )
      deallocate ( RT )
      deallocate ( uvecnm1, cnm1 )
      deallocate ( cnm1g, fnm1g )
      deallocate ( uhat )
      deallocate ( work6 )
      deallocate ( uvecmeshn, uvecmeshnm1 )

      deallocate ( gvecn, gvecnm1 )

      if ( vemcompressible ) then
        deallocate ( pr, press )
      end if

      if ( vemgammap ) then
        deallocate ( gpn, gpng )
        deallocate ( gpnm1, gpnm1g )
      end if

      call delete ( vemoptn, vemoptnm1 )

    end subroutine deallocate_arrays

  end subroutine ce_supg_elem2


! Internal element routine for the constitutive equations: SUPG, G-method
! Second-order semi-implicit time-integration

  subroutine ce_supg_elem_implicit_2nd_order ( mesh, problem, elgrp, elem, &
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


    integer :: i, j, ip, timeint, m, isv_gp
    real(dp) :: beta, fac, deltat, esize


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvel3D=1 )

    call set_devssg_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors, maxvar=1, gammap=.true. )

    if ( first ) then

!     check coefficients

      call check ( coefficients, 'ce_supg_elem_implicit_2nd_order', &
        indexarray=[48,50,90], minimum=[0,0,0], maximum=[2,1,1] )

      if ( coefficients%i(48) == 1 .and. coefficients%i(22) == 6 &
           .and. BLOCK_IMPROPER_ALE ) then
        write(*,'(/6(a/))') 'Error in ce_supg_elem_implicit_2nd_order:', &
          ' Combining ALE with Crank-Nicolson not allowed, ', &
          ' since the change in the nabla operator due to the change ', &
          ' in mesh at different times within a single step ', &
          ' has not been taken into account yet.', &
          ' Set BLOCK_IMPROPER_ALE=.false. in limits_m to bypass this block.'
        stop
      end if

!     in case the c/lambda terms are used, check if the modes are equal

      if ( coefficients%i(57) == 1 .and. coefficients%i(22) /= 7 ) then
        write(*,'(/4(a/))') 'Error in ce_supg_elem_implicit_2nd_order:', &
          'In this routine, adding of c/lambda terms (coefficients%i(57)=1)', &
          ' has only been implemented for second-order, semi-implicit', &
          ' Gear/conformation prediction (coefficients%i(22)=7).'
        stop
      end if

      if ( varpar ) then

!       check coefficients for variable coefficients
        call check ( coefficients, 'ce_supg_elem_implicit_2nd_order', &
          indexarray=[69,70,81], minimum=[0,0,0], maximum=[3,3,3] )

        if ( coefficients%i(22) == 6 ) then
          write(*,'(/4(a/))') 'Error in ce_supg_elem_implicit_2nd_order:', &
            ' Combining variable coefficients with Crank-Nicolson ', &
            ' not allowed, since the variable coefficients are only', &
            ' available for a single time step.'
          stop
        end if

      end if

      if ( coefficients%i(90) == 1 ) then
!       rotation reinitialization
        if ( coefficients%i(71) /= 1 ) then
          write(*,'(/3(a/))') 'Error in ce_supg_elem_implicit_2nd_order:', &
            ' Rotation reinitialization is only applicable for the ', &
            ' contravariant deformation tensor (CDT) formulation. '
          stop
        end if
!       exclude some time integration methods
        if ( coefficients%i(22) == 6 ) then
          write(*,'(/a/a,i0,a/)') 'Error in ce_supg_elem_implicit_2nd_order:', &
            ' Rotation reinitialization for time integration method = ', &
             coefficients%i(22), ' not available.'
          stop
        end if
!       exclude tALE
        if ( coefficients%i(48) == 2 ) then
          write(*,'(/3(a/))') 'Error in ce_supg_elem_implicit_2nd_order:', &
            ' Rotation reinitialization for the tALE method ', &
            ' not available.'
          stop
        end if
      end if

      if ( vemgammap .and. coefficients%i(48) == 2 ) then
        write(*,'(/2(a/))') 'Error in ce_supg_elem_implicit_2nd_order:', &
          ' Gammap dependence not available for the tALE method '
        stop
      end if

!     set globals

      cstorage = get_coefficient ( coefficients, index=84, default=1 )

      if ( vemgammap ) then
        gpstorage = get_coefficient ( coefficients, index=100, default=1 )
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


!   get relative change in volume J from pressure in system vector

    if ( vemcompressible ) then

      call get_J_from_pressure ( mesh, oldvectors%p(1)%p, elgrp, elem, &
        coefficients, oldvectors%s(1)%p, vemoptnp1 )

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
        call errormsg_case_default ( 'ce_supg_elem_implicit_2nd_order', &
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


!   Convection operator: ALE, temporary ALE or Eulerian frame

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

    else if ( coefficients%i(48) == 2 ) then

!     temporary ALE

!     compute reference and real coordinates in the temporary ALE meshes
!     NOTE: groups not yet used.

      call set_globals_stokes_tALE ( oldvectors%m(1)%p, oldvectors%m(2)%p, &
        oldvectors%m(3)%p, oldvectors%m(4)%p, deform_n=.false. )

!     compute mesh velocity at tn+1

      uvecmeshnp1 = ( 1.5_dp*xg - 2.0_dp*xg_n + 0.5_dp*xg_nm1 ) / deltat

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
          call errormsg_case_default ( 'ce_supg_elem_implicit_2nd_order', &
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

!     Convection operatore: Eulerian frame or ALE

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

      else if ( coefficients%i(48) == 2 ) then

!       temporary ALE

        write(*,'(/a/a,i0/)') 'Error in ce_supg_elem_implicit_2nd_order:', &
        ' temporary ALE not implemented for time integration scheme = ', timeint
        stop

      else

!       Eulerian frame

!       un.grad operator

        do ip = 1, ninti
          ungradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uvecn(ip,1:ndim) )
        end do

      end if

    end if


!   get conformation tensor at previous time steps (cn, cnm1)

    if ( coefficients%i(48) == 2 ) then

!     temporary ALE

!     shape functions for the conformation tensor

      call set_shape_function ( shapefuncc, xig_n, theta_n )
      call set_shape_function ( shapefuncc, xig_nm1, theta_nm1 )

      call get_c_tALE

    else

!     standard case (Euler or ALE)

      call get_c_standard

    end if


!   get gammap at previous time steps (cn, cnm1)

    if ( vemgammap ) then

      call get_gp_standard

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

      if ( coverl ) then
!       add c/lambda to both sides of the constitutive equation
        if ( coefficients%i(58) > 0 ) then ! add c/non-linear parameter
          lambda_g_mode1 = &
              [ ( mvemodel(ip)%nonlin(coefficients%i(58),mode1), ip=1,ninti ) ]
        else ! add c/lambda
          lambda_g_mode1 = [ ( mvemodel(ip)%lambda(mode1), ip=1,ninti ) ]
        end if
      end if

    else

      if ( coverl ) then
!       add c/lambda to both sides of the constitutive equation
        if ( coefficients%i(58) > 0 ) then ! add c/non-linear parameter
          lambda_g_mode1 = vemodel%nonlin(coefficients%i(58),mode1)
        else ! add c/lambda
          lambda_g_mode1 = vemodel%lambda(mode1)
        end if
      end if

    end if

    if ( coverl ) then
      fac2_g = 1._dp / lambda_g_mode1
    else
      fac2_g = 0
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
             ( ( 1.5_dp * fac + fac2_g ) * theta(:,j) + unp1gradtheta(:,j) ) &
                                   * detF * wg )
          end do
        end do

      case default

        write(*,'(/a/a,i0/)') 'Error in ce_supg_elem_implicit_2nd_order:', &
        ' incorrect time integration scheme  = ', timeint
        stop

      end select

    end if

    if ( vector ) then

      select case ( timeint )

      case(6)

!       semi-implicit Crank-Nicolson with relaxation prediction

!       L^n+1.chat^n+1 + chat^n+1.L^n+1T + f(chat^n+1)
        call rhs_viscoelastic ( gvecnp1, chatg, fhatg, vemopt=vemoptnp1 )
!       L^n.c^n + c^n.L^nT + f(c^n)
        call rhs_viscoelastic ( gvecn, cng, fng, vemopt=vemoptn )

        if ( coorsys == 1 .and. vel3D == 1 ) then
!         axisymmetric with 3D velocities: add terms due to u.nabla c
          call add_ugrad_terms ( xg, uvecnp1, chatg, fhatg )
          call add_ugrad_terms ( xg, uvecn, cng, fng )
        end if

        do m = mode1, mode2
          rhsmodel(:,:,m) = 0.5_dp*fhatg(:,:,m) + 0.5_dp*fng(:,:,m)
          do j = 1, ncomp
            do i = 1, ndfc
              work6(i,j,m-mode1+1) = &
                sum ( ( theta(:,i) + tau * unp1gradtheta(:,i) ) * &
                      ( fac * cng(:,j,m) - 0.5_dp*ungradcn(:,j,m) + &
                        rhsmodel(:,j,m) ) * detF * wg )
            end do
          end do
        end do

      case(7)

!       semi-implicit Gear with relaxation prediction

!       L^n+1.chat^n+1 + chat^n+1.L^n+1T + f(chat^n+1)
        call rhs_viscoelastic ( gvecnp1, chatg, fhatg, vemopt=vemoptnp1 )

        if ( coorsys == 1 .and. vel3D == 1 ) then
!         axisymmetric with 3D velocities: add terms due to u.nabla c
          call add_ugrad_terms ( xg, uvecnp1, chatg, fhatg )
        end if

        if ( coverl ) then
!         add c/lambda to both sides of the constitutive equation
          do ip = 1, ninti
            fhatg(ip,:,mode1) = fhatg(ip,:,mode1) + &
                                 chatg(ip,:,mode1) / lambda_g_mode1(ip)
          end do
        end if

        do m = mode1, mode2
          do j = 1, ncomp
            do i = 1, ndfc
              work6(i,j,m-mode1+1) = &
                sum ( ( theta(:,i) + tau * unp1gradtheta(:,i) ) * &
                        ( 2 * fac * cng(:,j,m) - 0.5_dp * fac * &
                          cnm1g(:,j,m) + fhatg(:,j,m) ) * detF * wg )
            end do
          end do
        end do

      case default

        write(*,'(/a/a,i0/)') 'Error in ce_supg_elem_implicit_2nd_order:', &
        ' incorrect time integration scheme  = ', timeint
        stop

      end select

      elemvec = reshape ( work6, [ ndfc * ncomp * ( mode2 - mode1 + 1 ) ] )

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

      allocate ( unp1gradtheta(ninti,ndfc) )
      allocate ( ungradcn(ninti,ncomp,nmodes) )
      allocate ( u(ncompu*ndf), um(ndim*ndf), g(ncompg*ndfg) )
      allocate ( gradu(ninti,ncompu,ncompu) )
      allocate ( tmp(ndf,ndim), work2(ndfg,ncompg), work4(ndf,ncompu) )
      allocate ( work10(ninti,ncompg) )
      allocate ( posu(ncompu*ndf), posg(ncompg*ndfg) )
      allocate ( tau(ninti), hoverU(ninti) )

      allocate ( uvecnp1(ninti,ncompu), cn(ndfc,ncomp) )
      allocate ( cng(ninti,ncomp,nmodes), fng(ninti,ncomp,nmodes) )
      allocate ( uvecn(ninti,ncompu), cnm1(ndfc,ncomp) )
      allocate ( cnm1g(ninti,ncomp,nmodes), fnm1g(ninti,ncomp,nmodes) )
      allocate ( RT(ninti,ncompu,ncompu) )
      allocate ( work6(ndfc,ncomp,mode2-mode1+1) )
      allocate ( facv(ninti) )

      if ( coorsys == 1 .and. vel3D == 0 ) then
        allocate ( gvecnp1(ninti,ndim**2+1), gvecn(ninti,ndim**2+1) )
      else if ( vel3D == 1 ) then
        allocate ( gvecnp1(ninti,ncompu**2), gvecn(ninti,ncompu**2) )
      else
        allocate ( gvecnp1(ninti,ndim**2), gvecn(ninti,ndim**2) )
      end if

      allocate ( ungradtheta(ninti,ndfc) )
      allocate ( chatg(ninti,ncomp,nmodes) )
      allocate ( fhatg(ninti,ncomp,nmodes), rhsmodel(ninti,ncomp,nmodes) )
      allocate ( uvecmeshn(ninti,ndim), uvecmeshnp1(ninti,ndim) )

      if ( coefficients%i(48) == 2 ) then

!       temporary ALE scheme; allocate additional arrays

        allocate ( xg_n(ninti,ndim), grpelm_n(ninti,2) )
        allocate ( xig_n(ninti,ndim), phi_n(ninti,ndf) )
        allocate ( xg_nm1(ninti,ndim), grpelm_nm1(ninti,2) )
        allocate ( xig_nm1(ninti,ndim), phi_nm1(ninti,ndf) )

        allocate ( theta_n(ninti,ndfc), theta_nm1(ninti,ndfc) )

      end if

      allocate ( lambda_g_mode1(ninti), fac2_g(ninti) )

      if ( varpar ) then
        allocate ( mvemodel(ninti) )
        mvemodel = vemodel
      end if

      if ( vemcompressible ) then
        allocate ( pr(ndfp), press(ninti) )
      end if

      if ( vemgammap ) then
        allocate ( gpn(ndfc), gpng(ninti) )
        allocate ( gpnm1(ndfc), gpnm1g(ninti) )
      end if

      call create_vemopt ( vemoptnp1, dep_J=vemcompressible, &
        dep_gammap=vemgammap, np=ninti )
      if ( coefficients%i(22) == 6 ) then
        call create_vemopt ( vemoptn, dep_J=vemcompressible, &
          dep_gammap=vemgammap, np=ninti )
      end if

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      integer :: ip

      deallocate ( unp1gradtheta )
      deallocate ( ungradcn )
      deallocate ( u, um, g )
      deallocate ( gradu )
      deallocate ( tmp, work2, work4 )
      deallocate ( work10 )
      deallocate ( posu, posg )
      deallocate ( tau, hoverU )

      deallocate ( uvecnp1, cn )
      deallocate ( cng, fng )
      deallocate ( uvecn, cnm1 )
      deallocate ( cnm1g, fnm1g )
      deallocate ( RT )
      deallocate ( work6 )
      deallocate ( facv )

      deallocate ( gvecnp1, gvecn )

      deallocate ( ungradtheta )
      deallocate ( chatg, fhatg, rhsmodel )
      deallocate ( uvecmeshn, uvecmeshnp1 )

      if ( coefficients%i(48) == 2 ) then

!       temporary ALE scheme; deallocate additional arrays

        deallocate ( xg_n, grpelm_n )
        deallocate ( xig_n, phi_n )
        deallocate ( xg_nm1, grpelm_nm1 )
        deallocate ( xig_nm1, phi_nm1 )

        deallocate ( theta_n, theta_nm1 )

      end if

      deallocate ( lambda_g_mode1, fac2_g )

      if ( varpar ) then
        do ip = 1, ninti
          call delete ( mvemodel(ip) )
        end do
        deallocate ( mvemodel )
      end if

      if ( vemcompressible ) deallocate ( pr, press )

      if ( vemgammap ) then
        deallocate ( gpn, gpng )
        deallocate ( gpnm1, gpnm1g )
      end if

      call delete ( vemoptnp1 )
      if ( coefficients%i(22) == 6 ) call delete ( vemoptn )

    end subroutine deallocate_arrays


!   get cn and cnm1 for the standard case (Euler or ALE)

    subroutine get_c_standard

      integer :: m

      do m = mode1, mode2

!       get conformation tensor at time step tn

!       single mode conformation
        call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
          cmode=cn, mode=m )

        cng(:,:,m) = matmul ( theta, cn )

!       get conformation tensor at time step tn-1

!       single mode conformation
        call get_conformation ( mesh, problem, elgrp, elem, oldvectors, &
          cmode=cnm1, mode=m, isv=2 )

        cnm1g(:,:,m) = matmul ( theta, cnm1 )

        if ( coefficients%i(90) == 1 ) then

!         rotation reinitialization in integration points

          call sqrtc_b ( cng(:,:,m), RT=RT )
          call sqrtc_b ( cnm1g(:,:,m) )
          call rotate_b ( cnm1g(:,:,m), RT )

        end if

!       prediction of cn+1 (chat)

        chatg(:,:,m) = 2 * cng(:,:,m) - cnm1g(:,:,m)

!       un.grad cn term (CN scheme only)

        if ( timeint == 6 ) then
          ungradcn(:,:,m) = matmul ( ungradtheta, cn )
        end if

      end do

    end subroutine get_c_standard


!   get cn and cnm1 for the temporary ALE case

    subroutine get_c_tALE

      integer :: ip, m

      do ip = 1, ninti

        do m = mode1, mode2

!         get conformation tensor at time step tn

!         single mode conformation
          call get_conformation ( oldvectors%m(2)%p, problem, &
            grpelm_n(ip,1), grpelm_n(ip,2), oldvectors, cmode=cn, mode=m )

          cng(ip,:,m) = matmul ( theta_n(ip,:), cn )

!         get conformation tensor at time step tn-1

!         single mode conformation
          call get_conformation ( oldvectors%m(4)%p, problem, &
            grpelm_nm1(ip,1), grpelm_nm1(ip,2), oldvectors, cmode=cnm1, &
            mode=m, isv=2 )

          cnm1g(ip,:,m) = matmul ( theta_nm1(ip,:), cnm1 )

!         prediction of cn+1 (chat)

          chatg(ip,:,m) = 2 * cng(ip,:,m) - cnm1g(ip,:,m)

        end do

      end do

    end subroutine get_c_tALE


!   get gpn and gpnm1 for the standard case (Euler or ALE)

    subroutine get_gp_standard

      isv_gp = get_coefficient ( coefficients, index=102, default=3 )

!     get gammap at tn

      call get_gammap ( mesh, oldvectors%p(3)%p, elgrp, elem, oldvectors, &
        gp=gpn, isv=isv_gp )

      gpng = matmul ( theta, gpn )

      if ( coefficients%i(22) == 6 ) vemoptn%gammap = gpng

!     get gammap at tn-1

      call get_gammap ( mesh, oldvectors%p(3)%p, elgrp, elem, oldvectors, &
        gp=gpnm1, isv=isv_gp+1 )

      gpnm1g = matmul ( theta, gpnm1 )

!     prediction of gpn+1 (gphat)

      vemoptnp1%gammap = 2 * gpng - gpnm1g

    end subroutine get_gp_standard


  end subroutine ce_supg_elem_implicit_2nd_order


! cylindrical, axisymmetric with 3D velocities: add terms due to u.nabla c

  subroutine add_ugrad_terms ( xg, uvec, cvec, fadd, dfadd, dfadd_mm )

    use viscoelastic_globals_m, only: mode1, mode2, bten, facv

    real(dp), dimension(:,:), intent(in) :: xg, uvec
    real(dp), dimension(:,:,:), intent(in) :: cvec
    real(dp), dimension(:,:,:), intent(inout) :: fadd
    real(dp), dimension(:,:,:,:), optional, intent(inout) :: dfadd
    real(dp), dimension(:,:,:,:,:), optional, intent(inout) :: dfadd_mm

    integer :: m

    facv = - uvec(:,3) / xg(:,2)  ! - u_theta / r

    if ( bten ) then

!     b-tensor formulation

      do m = mode1, mode2
        fadd(:,2,m) = fadd(:,2,m) - facv * cvec(:,3,m)     ! - b_z,theta
        fadd(:,3,m) = fadd(:,3,m) + facv * cvec(:,2,m)     ! + b_z,r
        fadd(:,4,m) = fadd(:,4,m) - facv * cvec(:,7,m)     ! - b_theta,z
        fadd(:,5,m) = fadd(:,5,m) &
             - facv * ( cvec(:,6,m) + cvec(:,8,m) ) ! - b_r,theta - b_theta,r
        fadd(:,6,m) = fadd(:,6,m) &
             + facv * ( cvec(:,5,m) - cvec(:,9,m) ) ! + b_rr - b_theta,theta
        fadd(:,7,m) = fadd(:,7,m) + facv * cvec(:,4,m)     ! + b_r,z
        fadd(:,8,m) = fadd(:,8,m) &
             + facv * ( cvec(:,5,m) - cvec(:,9,m) ) ! + b_rr - b_theta,theta
        fadd(:,9,m) = fadd(:,9,m) &
             + facv * ( cvec(:,6,m) + cvec(:,8,m) ) ! + b_r,theta + b_theta,r
      end do

      if ( present(dfadd) ) then

        do m = mode1, mode2
          dfadd(:,2,3,m) = dfadd(:,2,3,m) - facv  ! - b_z,theta
          dfadd(:,3,2,m) = dfadd(:,3,2,m) + facv  ! + b_z,r
          dfadd(:,4,7,m) = dfadd(:,4,7,m) - facv  ! - b_theta,z
          dfadd(:,5,6,m) = dfadd(:,5,6,m) - facv  ! - b_r,theta - b_theta,r
          dfadd(:,5,8,m) = dfadd(:,5,8,m) - facv  !
          dfadd(:,6,5,m) = dfadd(:,6,5,m) + facv  ! + b_rr - b_theta,theta
          dfadd(:,6,9,m) = dfadd(:,6,9,m) - facv
          dfadd(:,7,4,m) = dfadd(:,7,4,m) + facv  ! + b_r,z
          dfadd(:,8,5,m) = dfadd(:,8,5,m) + facv  ! + b_rr - b_theta,theta
          dfadd(:,8,9,m) = dfadd(:,8,9,m) - facv  !
          dfadd(:,9,6,m) = dfadd(:,9,6,m) + facv  ! + b_r,theta + b_theta,r
          dfadd(:,9,8,m) = dfadd(:,9,8,m) + facv  ! + b_r,theta + b_theta,r
        end do

      end if

      if ( present(dfadd_mm) ) then

        do m = mode1, mode2
          dfadd_mm(:,2,m,3,m) = dfadd_mm(:,2,m,3,m) - facv  ! - b_z,theta
          dfadd_mm(:,3,m,2,m) = dfadd_mm(:,3,m,2,m) + facv  ! + b_z,r
          dfadd_mm(:,4,m,7,m) = dfadd_mm(:,4,m,7,m) - facv  ! - b_theta,z
          dfadd_mm(:,5,m,6,m) = dfadd_mm(:,5,m,6,m) &
                                         - facv  ! - b_r,theta - b_theta,r
          dfadd_mm(:,5,m,8,m) = dfadd_mm(:,5,m,8,m) - facv  !
          dfadd_mm(:,6,m,5,m) = dfadd_mm(:,6,m,5,m) &
                                         + facv  ! + b_rr - b_theta,theta
          dfadd_mm(:,6,m,9,m) = dfadd_mm(:,6,m,9,m) - facv
          dfadd_mm(:,7,m,4,m) = dfadd_mm(:,7,m,4,m) + facv  ! + b_r,z
          dfadd_mm(:,8,m,5,m) = dfadd_mm(:,8,m,5,m) &
                                         + facv  ! + b_rr - b_theta,theta
          dfadd_mm(:,8,m,9,m) = dfadd_mm(:,8,m,9,m) - facv  !
          dfadd_mm(:,9,m,6,m) = dfadd_mm(:,9,m,6,m) &
                                         + facv  ! + b_r,theta + b_theta,r
          dfadd_mm(:,9,m,8,m) = dfadd_mm(:,9,m,8,m) &
                                         + facv  ! + b_r,theta + b_theta,r
        end do

      end if

    else

!     c-tensor formulation, including logc

      do m = mode1, mode2
        fadd(:,2,m) = fadd(:,2,m) - facv * cvec(:,3,m)     ! - c_z,theta
        fadd(:,3,m) = fadd(:,3,m) + facv * cvec(:,2,m)     ! + c_z,r
        fadd(:,4,m) = fadd(:,4,m) - 2 * facv * cvec(:,5,m) ! - 2 c_r,theta
        fadd(:,5,m) = fadd(:,5,m) &
                  + facv * ( cvec(:,4,m) - cvec(:,6,m) ) ! + c_rr-c_theta,theta
        fadd(:,6,m) = fadd(:,6,m) + 2 * facv * cvec(:,5,m) ! + 2 c_r,theta
      end do

      if ( present(dfadd) ) then

        do m = mode1, mode2
          dfadd(:,2,3,m) = dfadd(:,2,3,m) - facv     ! - c_z,theta
          dfadd(:,3,2,m) = dfadd(:,3,2,m) + facv     ! + c_z,r
          dfadd(:,4,5,m) = dfadd(:,4,5,m) - 2 * facv ! - 2 c_r,theta
          dfadd(:,5,4,m) = dfadd(:,5,4,m) + facv     ! + c_rr-c_theta,theta
          dfadd(:,5,6,m) = dfadd(:,5,6,m) - facv     !
          dfadd(:,6,5,m) = dfadd(:,6,5,m) + 2 * facv ! + 2 c_r,theta
        end do

      end if

      if ( present(dfadd_mm) ) then

        do m = mode1, mode2
          dfadd_mm(:,2,m,3,m) = dfadd_mm(:,2,m,3,m) - facv     ! - c_z,theta
          dfadd_mm(:,3,m,2,m) = dfadd_mm(:,3,m,2,m) + facv     ! + c_z,r
          dfadd_mm(:,4,m,5,m) = dfadd_mm(:,4,m,5,m) - 2 * facv ! - 2 c_r,theta
          dfadd_mm(:,5,m,4,m) = dfadd_mm(:,5,m,4,m) &
                                             + facv     ! + c_rr-c_theta,theta
          dfadd_mm(:,5,m,6,m) = dfadd_mm(:,5,m,6,m) - facv     !
          dfadd_mm(:,6,m,5,m) = dfadd_mm(:,6,m,5,m) + 2 * facv ! + 2 c_r,theta
        end do

      end if

    end if

  end subroutine add_ugrad_terms

end module viscoelastic_elements_SUPG_m

