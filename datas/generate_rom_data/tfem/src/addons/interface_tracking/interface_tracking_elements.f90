! Copyright (C) 2012-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Element routines for tracking an interface according to
!
! .
! x = ( u - u ) . n n + u  + u
! -     -   -m    - -   -m   -t
!
! where
!
!   x   the position vector of the interface, which basically represents an
!   -   ALE description of the interface.
!
!   n   the normal vector on the interface
!   -
!   u   the velocity vector
!   -
!   u   an arbitrary vector field.
!   -m
!   u   a vector field tangential to the interface proportional to the
!   -t  surface gradient of scalar field c:
!
!           u  ~ nabla  c
!           -t     -  s
!
! Notes:
!    .
! 1) x.n = u, i.e. the interface is a material interface.
!    - -   -
! 2) the normal vector n is determined from the geometry of the interface.
!                      -
! 3) the tangential motion of the interface is given by
!
!          ( I - nn ) . u  + u
!            =   --     -m   -t
!    which contributes only to the redistribution of the nodes (ALE grid).
!                           .
! 4) if u  = u, u = 0, then x = u and the interface moves in a Lagrangian way.
!       -m   -  -t  -       -   -
!
! 5) The interface tracking equation can also be written as
!      .
!      x + ( u - u ) . (1 - n n) = u + u
!      -     -   -m     =   - -    -   -t
!    or
!      .                j
!      x + ( u - u ) . g g  = u + u
!      -     -   -m    - -j   -   -t
!    or
!      .                j      j
!      x + ( u - u ) . g dx/dxi  = u + u
!      -     -   -m    -  -        -   -t
!    or
!      .               __
!      x + ( u - u ) . \/ x  = u + u
!      -     -   -m     s -    -   -t
!
!    which is a (non-linear) convection equation for x on the surface.
!

module interface_tracking_elements_m

  use tfem_elem_m

  implicit none

contains


! Interface tracking for 2D and 3D flow.
! Solve for xdot for explicit time integration defined by the user.

  subroutine interface_tracking_elem ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use interface_tracking_globals_m

!   input/output
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, ip, k, m
    real(dp) :: fac


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     check size of coefficients

      call check ( coefficients, 'interface_tracking_elem', &
        indexarray=[10,12], minimum=[0,0], maximum=[1,1] )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( intpolv > 0 ) then
        call set_shape_function ( shapefuncv, xig, phiv, dphiv )
      else
        phiv = phi
        dphiv = dphi
      end if

    end if

!   get coordinates of nodes of interface

    call get_coordinates ( mesh, elgrp, elem, x )

    if ( isoshape == 0 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, normal, &
        gi_up=gi_up )

!     coordinates of integration points
      xg = matmul ( phi, x )

    else if ( isoshape == 1 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphiv, dxdxis, surfl, normal, &
        gi_up=gi_up )

!     coordinates of integration points
      xg = matmul ( phiv, x )

    end if

!   get velocities u
!                  -
    call evaluate_vector_coefficient ( mesh, problem, oldvectors, elgrp, elem, &
      choice=coefficients%i(6), value=coefficients%r(1:ndim), &
      vfunc=coefficients%vfunc1(1)%p, vfuncnr=coefficients%i(7), x=xg, &
      indx_v=1, layer=0, phi=phiv, indx_e=1, coef=uvecg )

!   get velocities u
!                  -m
    call evaluate_vector_coefficient ( mesh, problem, oldvectors, elgrp, elem, &
      choice=coefficients%i(8), value=coefficients%r(4:3+ndim), &
      vfunc=coefficients%vfunc1(2)%p, vfuncnr=coefficients%i(9), x=xg, &
      indx_v=2, layer=0, phi=phiv, indx_e=2, coef=uvecmg )


    if ( coefficients%i(10) == 1 ) then

!     get solution of Poisson problem on interface (c), from which deriving u
!                                                                           -t
      call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, elgrp, &
        elem, c )

!     compute tangential "diffusion" velocity of nodes u  = fac * nabla  c in
!                                                      -t         -----s
!     integration points

      fac = coefficients%r(7)

      do i = 1, ndim - 1
        work3(:,i) = matmul ( dphi(:,:,i), c )
      end do

      do ip = 1, ninti
        uvectg(ip,:) = fac * matmul ( gi_up(ip,:,:), work3(ip,:) )
      end do

      uvecg_rhs = uvecmg + uvectg  ! add vectors u_m + u_t

    else

      uvecg_rhs = uvecmg  ! only vector u_m

    end if

!   build matrix and vector to solve equation through FEM

    if ( matrix ) then

      do i = 1, ndf
        do j = 1, ndf
          elemmat(i,j) = sum ( phi(:,i) * phi(:,j) * wg * surfl )
        end do
      end do

      if ( coefficients%i(12) == 1 ) then

!       full matrix with all components

!       diagonal blocks

        do j = 2, ndim
          k = ndf * (j-1)
          elemmat(k+1:k+ndf,k+1:k+ndf) = elemmat(1:ndf,1:ndf)
        end do

!       off-diagonal blocks

        do i = 1, ndim-1
          k = ndf * (i-1)
          do j = i+1, ndim
            m = ndf * (j-1)
            elemmat(k+1:k+ndf,m+1:) = 0
            elemmat(m+1:,k+1:k+ndf) = 0  ! transposed
          end do
        end do

      end if

    end if

    if ( vector ) then

      do ip = 1, ninti
        un(ip) = dot_product(uvecg(ip,:)-uvecmg(ip,:),normal(ip,:))
      end do

      do j = 1, ndim
        k = (j-1) * ndf
        do i = 1, ndf
          elemvec(k+i) = sum ( phi(:,i) * wg * surfl * &
                         ( un * normal(:,j) + uvecg_rhs(:,j) ) )
        end do
      end do

    end if

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )
      allocate ( normal(ninti,ndim) )
      allocate ( un(ninti) )
      allocate ( xg(ninti,ndim) )
      allocate ( uvecmg(ninti,ndim), uvecg(ninti,ndim), uvectg(ninti,ndim) )
      allocate ( uvecg_rhs(ninti,ndim) )
      allocate ( c(ndf) )
      allocate ( work3(ninti,ndim-1) )
      allocate ( gi_up(ninti,ndim,ndim-1) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x )
      deallocate ( dxdxis )
      deallocate ( surfl )
      deallocate ( normal )
      deallocate ( un )
      deallocate ( xg )
      deallocate ( uvecg, uvecmg, uvectg, uvecg_rhs )
      deallocate ( c )
      deallocate ( work3 )
      deallocate ( gi_up )

    end subroutine deallocate_arrays

  end subroutine interface_tracking_elem


! Interface tracking for 2D and 3D flow. SUPG stabilization.
! Solve for xdot for explicit time integration defined by the user.

  subroutine interface_tracking_elem_supg ( mesh, problem, elgrp, &
    elem, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemvec )

    use interface_tracking_globals_m

!   input/output
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, ip, k, m
    real(dp) :: fac, h, beta


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     check size of coefficients

      call check ( coefficients, 'interface_tracking_elem_supg', &
        indexarray=[10,12], minimum=[0,0], maximum=[1,1] )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( intpolv > 0 ) then
        call set_shape_function ( shapefuncv, xig, phiv, dphiv )
      else
        phiv = phi
        dphiv = dphi
      end if

    end if

!   get coordinates of nodes of interface

    call get_coordinates ( mesh, elgrp, elem, x )

    if ( isoshape == 0 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, normal, &
        gij_down=gij_down, gi_up=gi_up )

!     coordinates of integration points
      xg = matmul ( phi, x )

    else if ( isoshape == 1 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphiv, dxdxis, surfl, normal, &
        gij_down=gij_down, gi_up=gi_up )

!     coordinates of integration points
      xg = matmul ( phiv, x )

    end if

!   get velocities u
!                  -
    call evaluate_vector_coefficient ( mesh, problem, oldvectors, elgrp, elem, &
      choice=coefficients%i(6), value=coefficients%r(1:ndim), &
      vfunc=coefficients%vfunc1(1)%p, vfuncnr=coefficients%i(7), x=xg, &
      indx_v=1, layer=0, phi=phiv, indx_e=1, coef=uvecg )

!   get velocities u
!                  -m
    call evaluate_vector_coefficient ( mesh, problem, oldvectors, elgrp, elem, &
      choice=coefficients%i(8), value=coefficients%r(4:3+ndim), &
      vfunc=coefficients%vfunc1(2)%p, vfuncnr=coefficients%i(9), x=xg, &
      indx_v=2, layer=0, phi=phiv, indx_e=2, coef=uvecmg )


    if ( coefficients%i(10) == 1 ) then

!     get solution of Poisson problem on interface (c), from which deriving u
!                                                                           -t
      call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, elgrp, &
        elem, c )

!     compute tangential "diffusion" velocity of nodes u  = fac * nabla  c in
!                                                      -t         -----s
!     integration points

      fac = coefficients%r(7)

      do i = 1, ndim - 1
        work3(:,i) = matmul ( dphi(:,:,i), c )
      end do

      do ip = 1, ninti
        uvectg(ip,:) = fac * matmul ( gi_up(ip,:,:), work3(ip,:) )
      end do

      uvectg = uvecmg + uvectg  ! add vectors u_m + u_t

    else

      uvectg = uvecmg  ! only vector u_m

    end if

!   characteristic element length h

    select case (ndim)
    case(2)
      h = sum ( wg * surfl )
    case(3)
      h = sqrt ( sum ( wg * surfl ) )
    case default
      call errormsg_case_default ( 'interface_tracking_elem_supg', &
        'ndim', int_value=ndim )
    end select

!   characteristic velocity = norm of tangential velocity

    do ip = 1, ninti

!     u^i = ( u - um ) . g^i
      work3(ip,:) = matmul ( uvecg(ip,:) - uvecmg(ip,:), gi_up(ip,:,:) )

!     U_c=sqrt(u^ig_i.u^jg_j)=sqrt(g_ij u^iu^j)
      Uc(ip) = dot_product ( work3(ip,:), &
                                    matmul( gij_down(ip,:,:), work3(ip,:) ) )
    end do

    Uc = sqrt(Uc)

!   factor for supg

    beta = coefficients%r(8)
    facs = beta * h / Uc / 2

!   (u-um).g^k dphi_i/xi^k

    do ip = 1, ninti
      ugradphi(ip,:) = matmul ( dphi(ip,:,:), work3(ip,:) )
    end do

!   build matrix and vector to solve equation through FEM

    if ( matrix ) then

      do i = 1, ndf
        do j = 1, ndf
          elemmat(i,j) = sum ( ( phi(:,i) + facs * ugradphi(:,i) ) * &
                                   phi(:,j) * wg * surfl )
        end do
      end do

      if ( coefficients%i(12) == 1 ) then

!       full matrix with all components

!       diagonal blocks

        do j = 2, ndim
          k = ndf * (j-1)
          elemmat(k+1:k+ndf,k+1:k+ndf) = elemmat(1:ndf,1:ndf)
        end do

!       off-diagonal blocks

        do i = 1, ndim-1
          k = ndf * (i-1)
          do j = i+1, ndim
            m = ndf * (j-1)
            elemmat(k+1:k+ndf,m+1:) = 0
            elemmat(m+1:,k+1:k+ndf) = 0  ! transposed
          end do
        end do

      end if

    end if

    if ( vector ) then

      do ip = 1, ninti
        un(ip) = dot_product(uvecg(ip,:)-uvecmg(ip,:),normal(ip,:))
      end do

      do j = 1, ndim
        k = (j-1) * ndf
        do i = 1, ndf
          elemvec(k+i) = sum ( ( phi(:,i) + facs * ugradphi(:,i) ) * &
                         ( un * normal(:,j) + uvectg(:,j) ) * wg * surfl )
        end do
      end do

    end if

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )
      allocate ( normal(ninti,ndim) )
      allocate ( un(ninti) )
      allocate ( xg(ninti,ndim) )
      allocate ( uvecmg(ninti,ndim), uvecg(ninti,ndim), uvectg(ninti,ndim) )
      allocate ( c(ndf) )
      allocate ( work3(ninti,ndim-1) )
      allocate ( gi_up(ninti,ndim,ndim-1) )
      allocate ( gij_down(ninti,ndim-1,ndim-1) )
      allocate ( Uc(ninti), facs(ninti) )
      allocate ( ugradphi(ninti,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x )
      deallocate ( dxdxis )
      deallocate ( surfl )
      deallocate ( normal )
      deallocate ( un )
      deallocate ( xg )
      deallocate ( uvecg, uvecmg, uvectg )
      deallocate ( c )
      deallocate ( work3 )
      deallocate ( gi_up )
      deallocate ( gij_down )
      deallocate ( Uc, facs )
      deallocate ( ugradphi )

    end subroutine deallocate_arrays

  end subroutine interface_tracking_elem_supg


! Interface tracking for 2D and 3D flow. SUPG stabilization.
! (Semi-)Implicit time integration based on Euler and Gear (2nd order) included
! on element level.

  subroutine interface_tracking_elem_supg_impl ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use interface_tracking_globals_m

!   input/output
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, ip, k, timeint, m
    real(dp) :: fac, h, beta, deltat


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     check size of coefficients

      call check ( coefficients, 'interface_tracking_elem_supg_impl', &
        indexarray=[10,11,12], minimum=[0,1,0], maximum=[1,2,1] )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( intpolv > 0 ) then
        call set_shape_function ( shapefuncv, xig, phiv, dphiv )
      else
        phiv = phi
        dphiv = dphi
      end if

    end if

!   start of the element

    timeint = coefficients%i(11)

!   get coordinates of nodes of interface

    call get_coordinates ( mesh, elgrp, elem, x )

    if ( isoshape == 0 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, &
        gij_down=gij_down, gi_up=gi_up )

!     coordinates of integration points
      xg = matmul ( phi, x )

!     position at tn

      call get_vector ( mesh, problem, oldvectors%v(3)%p, elgrp, elem, u=xar )

      xng = matmul ( phi, reshape ( xar, [ndf,ndim] ) )

      if ( timeint == 2 ) then

!       second-order time discretization

!       position at tn-1

        call get_vector ( mesh, problem, oldvectors%v(4)%p, elgrp, elem, u=xar )

        xnm1g = matmul ( phi, reshape ( xar, [ndf,ndim] ) )

      end if

    else if ( isoshape == 1 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphiv, dxdxis, surfl, &
        gij_down=gij_down, gi_up=gi_up )

!     coordinates of integration points
      xg = matmul ( phiv, x )

!     position at tn

      call get_vector ( mesh, problem, oldvectors%v(3)%p, elgrp, elem, u=xar )

      xng = matmul ( phiv, reshape ( xar, [ndfv,ndim] ) )

      if ( timeint == 2 ) then

!       second-order time discretization

!       position at tn-1

        call get_vector ( mesh, problem, oldvectors%v(4)%p, elgrp, elem, u=xar )

        xnm1g = matmul ( phiv, reshape ( xar, [ndfv,ndim] ) )

      end if

    end if

!   get velocities u
!                  -
    call evaluate_vector_coefficient ( mesh, problem, oldvectors, elgrp, elem, &
      choice=coefficients%i(6), value=coefficients%r(1:ndim), &
      vfunc=coefficients%vfunc1(1)%p, vfuncnr=coefficients%i(7), x=xg, &
      indx_v=1, layer=0, phi=phiv, indx_e=1, coef=uvecg )

!   get velocities u
!                  -m
    call evaluate_vector_coefficient ( mesh, problem, oldvectors, elgrp, elem, &
      choice=coefficients%i(8), value=coefficients%r(4:3+ndim), &
      vfunc=coefficients%vfunc1(2)%p, vfuncnr=coefficients%i(9), x=xg, &
      indx_v=2, layer=0, phi=phiv, indx_e=2, coef=uvecmg )

!   set right-hand side

    if ( coefficients%i(10) == 1 ) then

!     get solution of Poisson problem on interface (c), from which deriving u
!                                                                           -t
      call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, elgrp, &
        elem, c )

!     compute tangential "diffusion" velocity of nodes u  = fac * nabla  c in
!                                                      -t         -----s
!     integration points

      fac = coefficients%r(7)

      do i = 1, ndim - 1
        work3(:,i) = matmul ( dphi(:,:,i), c )
      end do

      do ip = 1, ninti
        uvectg(ip,:) = fac * matmul ( gi_up(ip,:,:), work3(ip,:) )
      end do

      uvecg_rhs = uvecg + uvectg

    else

      uvecg_rhs = uvecg

    end if

!   convection velocity

    uvecg_conv = uvecg - uvecmg

!   characteristic element length h

    select case (ndim)
    case(2)
      h = sum ( wg * surfl )
    case(3)
      h = sqrt ( sum ( wg * surfl ) )
    case default
      call errormsg_case_default ( 'interface_tracking_elem_supg_impl', &
        'ndim', int_value=ndim )
    end select

!   characteristic velocity = norm of tangential velocity

    do ip = 1, ninti

!     u^i = ( u - um ) . g^i
      work3(ip,:) = matmul ( uvecg_conv(ip,:), gi_up(ip,:,:) )

!     U_c=sqrt(u^ig_i.u^jg_j)=sqrt(g_ij u^iu^j)
      Uc(ip) = dot_product ( work3(ip,:), &
                                    matmul( gij_down(ip,:,:), work3(ip,:) ) )
    end do

    Uc = sqrt(Uc)

!   factor for supg

    beta = coefficients%r(8)
    facs = beta * h / Uc / 2

!   (u-um).g^k dphi_i/xi^k

    do ip = 1, ninti
      ugradphi(ip,:) = matmul ( dphi(ip,:,:), work3(ip,:) )
    end do

!   set time step

    deltat = coefficients%r(9)

!   build matrix and vector to solve equation through FEM

    if ( matrix ) then

      select case ( timeint)

      case(1) ! Euler implicit

        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( ( phi(:,i) + facs * ugradphi(:,i) ) * &
                  ( phi(:,j) / deltat + ugradphi(:,j) ) * wg * surfl )
          end do
        end do

      case(2) ! 2nd order Gear

        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( ( phi(:,i) + facs * ugradphi(:,i) ) * &
                  ( 1.5_dp * phi(:,j) / deltat + ugradphi(:,j) ) * wg * surfl )
          end do
        end do

      case default

        call errormsg_case_default ( 'interface_tracking_elem_supg_impl', &
          'timeint', int_value=timeint )

      end select

      if ( coefficients%i(12) == 1 ) then

!       full matrix with all components

!       diagonal blocks

        do j = 2, ndim
          k = ndf * (j-1)
          elemmat(k+1:k+ndf,k+1:k+ndf) = elemmat(1:ndf,1:ndf)
        end do

!       off-diagonal blocks

        do i = 1, ndim-1
          k = ndf * (i-1)
          do j = i+1, ndim
            m = ndf * (j-1)
            elemmat(k+1:k+ndf,m+1:) = 0
            elemmat(m+1:,k+1:k+ndf) = 0  ! transposed
          end do
        end do

      end if

    end if

    if ( vector ) then

      select case ( timeint)

      case(1) ! Euler implicit

        do j = 1, ndim
          k = (j-1) * ndf
          do i = 1, ndf
            elemvec(k+i) = sum ( ( phi(:,i) + facs * ugradphi(:,i) ) * &
                        ( xng(:,j) / deltat + uvecg_rhs(:,j) ) * wg * surfl )
          end do
        end do

      case(2) ! 2nd order Gear

        do j = 1, ndim
          k = (j-1) * ndf
          do i = 1, ndf
            elemvec(k+i) = sum ( ( phi(:,i) + facs * ugradphi(:,i) ) * &
                        ( ( 2 * xng(:,j) - 0.5_dp * xnm1g(:,j) ) / deltat + &
                                     uvecg_rhs(:,j) ) * wg * surfl )
          end do
        end do

      case default

        call errormsg_case_default ( 'interface_tracking_elem_supg_impl', &
          'timeint', int_value=timeint )

      end select

    end if

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( xar(max(ndf,ndfv)*ndim) )
      allocate ( xng(ninti,ndim), xnm1g(ninti,ndim) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )
      allocate ( xg(ninti,ndim) )
      allocate ( uvecmg(ninti,ndim), uvecg(ninti,ndim), uvectg(ninti,ndim) )
      allocate ( uvecg_conv(ninti,ndim), uvecg_rhs(ninti,ndim) )
      allocate ( c(ndf) )
      allocate ( work3(ninti,ndim-1) )
      allocate ( gi_up(ninti,ndim,ndim-1) )
      allocate ( gij_down(ninti,ndim-1,ndim-1) )
      allocate ( Uc(ninti), facs(ninti) )
      allocate ( ugradphi(ninti,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x )
      deallocate ( xar )
      deallocate ( xng, xnm1g )
      deallocate ( dxdxis )
      deallocate ( surfl )
      deallocate ( xg )
      deallocate ( uvecg, uvecmg, uvectg )
      deallocate ( uvecg_conv, uvecg_rhs )
      deallocate ( c )
      deallocate ( work3 )
      deallocate ( gi_up )
      deallocate ( gij_down )
      deallocate ( Uc, facs )
      deallocate ( ugradphi )

    end subroutine deallocate_arrays

  end subroutine interface_tracking_elem_supg_impl


! Interface tracking for 2D and 3D flow. Lagrange based: u  = u, and thus
!                                                        -m   -
!    .
!    x = u  + u
!    -   -    -t
!
! Solve for xdot for explicit time integration defined by the user.

  subroutine interface_tracking_elem_xdot_lagrange ( mesh, problem, elgrp, &
    elem, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemvec )

    use interface_tracking_globals_m

!   input/output
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, ip, k, m
    real(dp) :: fac


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     check size of coefficients

      call check ( coefficients, 'interface_tracking_elem_xdot_lagrange', &
        indexarray=[10,12], minimum=[0,0], maximum=[1,1] )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( intpolv > 0 ) then
        call set_shape_function ( shapefuncv, xig, phiv, dphiv )
      else
        phiv = phi
        dphiv = dphi
      end if

    end if

!   get coordinates of nodes of interface

    call get_coordinates ( mesh, elgrp, elem, x )

    if ( isoshape == 0 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, &
        gi_up=gi_up )

!     coordinates of integration points
      xg = matmul ( phi, x )

    else if ( isoshape == 1 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphiv, dxdxis, surfl, &
        gi_up=gi_up )

!     coordinates of integration points
      xg = matmul ( phiv, x )

    end if

!   get velocities u
!                  -
    call evaluate_vector_coefficient ( mesh, problem, oldvectors, elgrp, elem, &
      choice=coefficients%i(6), value=coefficients%r(1:ndim), &
      vfunc=coefficients%vfunc1(1)%p, vfuncnr=coefficients%i(7), x=xg, &
      indx_v=1, layer=0, phi=phiv, indx_e=1, coef=uvecg )

    if ( coefficients%i(10) == 1 ) then

!     get solution of Poisson problem on interface (c), from which deriving u
!                                                                           -t
      call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, elgrp, &
        elem, c )

!     compute tangential "diffusion" velocity of nodes u  = fac * nabla  c in
!                                                      -t         -----s
!     integration points

      fac = coefficients%r(7)

      do i = 1, ndim - 1
        work3(:,i) = matmul ( dphi(:,:,i), c )
      end do

      do ip = 1, ninti
        uvectg(ip,:) = fac * matmul ( gi_up(ip,:,:), work3(ip,:) )
      end do

      uvecg_rhs = uvecg + uvectg  ! add vectors u + u_t

    else

      uvecg_rhs = uvecg  ! only vector u

    end if

!   build matrix and vector to solve equation through FEM

    if ( matrix ) then

      do i = 1, ndf
        do j = 1, ndf
          elemmat(i,j) = sum ( phi(:,i) * phi(:,j) * wg * surfl )
        end do
      end do

      if ( coefficients%i(12) == 1 ) then

!       full matrix with all components

!       diagonal blocks

        do j = 2, ndim
          k = ndf * (j-1)
          elemmat(k+1:k+ndf,k+1:k+ndf) = elemmat(1:ndf,1:ndf)
        end do

!       off-diagonal blocks

        do i = 1, ndim-1
          k = ndf * (i-1)
          do j = i+1, ndim
            m = ndf * (j-1)
            elemmat(k+1:k+ndf,m+1:) = 0
            elemmat(m+1:,k+1:k+ndf) = 0  ! transposed
          end do
        end do

      end if

    end if

    if ( vector ) then

      do j = 1, ndim
        k = (j-1) * ndf
        do i = 1, ndf
          elemvec(k+i) = sum ( phi(:,i) * uvecg_rhs(:,j) * wg * surfl )
        end do
      end do

    end if

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )
      allocate ( xg(ninti,ndim) )
      allocate ( uvecg_rhs(ninti,ndim), uvecg(ninti,ndim), uvectg(ninti,ndim) )
      allocate ( c(ndf) )
      allocate ( work3(ninti,ndim-1) )
      allocate ( gi_up(ninti,ndim,ndim-1) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x )
      deallocate ( dxdxis )
      deallocate ( surfl )
      deallocate ( xg )
      deallocate ( uvecg, uvecg_rhs, uvectg )
      deallocate ( c )
      deallocate ( work3 )
      deallocate ( gi_up )

    end subroutine deallocate_arrays

  end subroutine interface_tracking_elem_xdot_lagrange


! Compute interface position in all nodes (component version)

  subroutine position_component_deriv (  mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec, elemwts )

    use interface_tracking_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

!     shape function in the nodes

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefunc, xrnod, phi )

    end if

!   position

    call get_sysvector ( mesh, problem, oldvectors%s(2)%p, elgrp, elem, u=xar )

    elemvec = matmul ( phi, xar )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xrnod(nodalp,ndim-1) )
      allocate ( phi(nodalp,ndf) )
      allocate ( xar(ndf), x(nodalp,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xrnod )
      deallocate ( phi )
      deallocate ( xar, x )

    end subroutine deallocate_arrays

  end subroutine position_component_deriv


! Compute interface position in all nodes (vector version)

  subroutine position_vector_deriv (  mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use interface_tracking_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

!     shape function in the nodes

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefunc, xrnod, phi )

    end if

!   position

    call get_sysvector ( mesh, problem, oldvectors%s(2)%p, elgrp, elem, u=xar )

    x = matmul ( phi, reshape ( xar, [ndf,ndim] ) )

    elemvec = reshape ( x, [ nodalp*ndim ]  )

    elemwts = 1

    if ( last ) then

!     last element in this group

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xrnod(nodalp,ndim-1) )
      allocate ( phi(nodalp,ndf) )
      allocate ( xar(ndf*ndim), x(nodalp,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xrnod )
      deallocate ( phi )
      deallocate ( xar, x )

    end subroutine deallocate_arrays

  end subroutine position_vector_deriv


! Interface tracking for 2D and 3D flow. Lagrange based: u  = u, and thus
!                                                        -m   -
!    .
!    x = u  + u
!    -   -    -t
!
! Explicit time integration based on Euler and Gear (2nd order) included
! on element level.

  subroutine interface_tracking_elem_lagrange ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use interface_tracking_globals_m

!   input/output
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, ip, k, timeint, m
    real(dp) :: deltat, fac


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     check size of coefficients

      call check ( coefficients, 'interface_tracking_elem_lagrange', &
        indexarray=[10,11,12], minimum=[0,1,0], maximum=[1,2,1] )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( intpolv > 0 ) then
        call set_shape_function ( shapefuncv, xig, phiv, dphiv )
      else
        phiv = phi
        dphiv = dphi
      end if

    end if

!   start of the element

    timeint = coefficients%i(11)

!   get coordinates of nodes of interface

    call get_coordinates ( mesh, elgrp, elem, x )

    if ( isoshape == 0 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, &
        gi_up=gi_up )

!     coordinates of integration points
      xg = matmul ( phi, x )

    else if ( isoshape == 1 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphiv, dxdxis, surfl, &
        gi_up=gi_up )

!     coordinates of integration points
      xg = matmul ( phiv, x )

    end if

!   get velocities u
!                  -
    call evaluate_vector_coefficient ( mesh, problem, oldvectors, elgrp, elem, &
      choice=coefficients%i(6), value=coefficients%r(1:ndim), &
      vfunc=coefficients%vfunc1(1)%p, vfuncnr=coefficients%i(7), x=xg, &
      indx_v=1, layer=0, phi=phiv, indx_e=1, coef=uvecg )

!   position at tn

    call get_vector ( mesh, problem, oldvectors%v(3)%p, elgrp, elem, u=xar )

    xng = matmul ( phi, reshape ( xar, [ndf,ndim] ) )

    if ( timeint == 2 ) then

!     second-order time discretization

!     position at tn-1

      call get_vector ( mesh, problem, oldvectors%v(4)%p, elgrp, elem, u=xar )

      xnm1g = matmul ( phi, reshape ( xar, [ndf,ndim] ) )

    end if

!   set right-hand side

    if ( coefficients%i(10) == 1 ) then

!     get solution of Poisson problem on interface (c), from which deriving u
!                                                                           -t
      call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, elgrp, &
        elem, c )

!     compute tangential "diffusion" velocity of nodes u  = fac * nabla  c in
!                                                      -t         -----s
!     integration points

      fac = coefficients%r(7)

      do i = 1, ndim - 1
        work3(:,i) = matmul ( dphi(:,:,i), c )
      end do

      do ip = 1, ninti
        uvectg(ip,:) = fac * matmul ( gi_up(ip,:,:), work3(ip,:) )
      end do

      uvecg_rhs = uvecg + uvectg

    else

      uvecg_rhs = uvecg

    end if

!   set time step

    deltat = coefficients%r(9)

!   build matrix and vector to solve equation through FEM

    if ( matrix ) then

      do i = 1, ndf
        do j = i, ndf
          elemmat(i,j) = sum ( phi(:,i) * phi(:,j) * wg * surfl )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

      select case ( timeint)

      case(1) ! Euler

        elemmat(:ndf,:ndf) = elemmat(:ndf,:ndf) / deltat

      case(2) ! 2nd order Gear

        elemmat(:ndf,:ndf) = 1.5_dp * elemmat(:ndf,:ndf) / deltat

      case default

        call errormsg_case_default ( 'interface_tracking_elem_lagrange', &
          'timeint', int_value=timeint )

      end select

      if ( coefficients%i(12) == 1 ) then

!       full matrix with all components

!       diagonal blocks

        do j = 2, ndim
          k = ndf * (j-1)
          elemmat(k+1:k+ndf,k+1:k+ndf) = elemmat(1:ndf,1:ndf)
        end do

!       off-diagonal blocks

        do i = 1, ndim-1
          k = ndf * (i-1)
          do j = i+1, ndim
            m = ndf * (j-1)
            elemmat(k+1:k+ndf,m+1:) = 0
            elemmat(m+1:,k+1:k+ndf) = 0  ! transposed
          end do
        end do

      end if

    end if

    if ( vector ) then

      select case ( timeint)

      case(1) ! Euler implicit

        do j = 1, ndim
          k = (j-1) * ndf
          do i = 1, ndf
            elemvec(k+i) = sum ( phi(:,i) * &
                        ( xng(:,j) / deltat + uvecg_rhs(:,j) ) * wg * surfl )
          end do
        end do

      case(2) ! 2nd order Gear

        do j = 1, ndim
          k = (j-1) * ndf
          do i = 1, ndf
            elemvec(k+i) = sum ( phi(:,i) * &
                        ( ( 2 * xng(:,j) - 0.5_dp * xnm1g(:,j) ) / deltat + &
                                     uvecg_rhs(:,j) ) * wg * surfl )
          end do
        end do

      case default

        call errormsg_case_default ( 'interface_tracking_elem_lagrange', &
          'timeint', int_value=timeint )

      end select

    end if

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( xar(ndf*ndim) )
      allocate ( xng(ninti,ndim), xnm1g(ninti,ndim) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )
      allocate ( xg(ninti,ndim) )
      allocate ( uvecg(ninti,ndim), uvectg(ninti,ndim) )
      allocate ( uvecg_rhs(ninti,ndim) )
      allocate ( c(ndf) )
      allocate ( work3(ninti,ndim-1) )
      allocate ( gi_up(ninti,ndim,ndim-1) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x )
      deallocate ( xar )
      deallocate ( xng, xnm1g )
      deallocate ( dxdxis )
      deallocate ( surfl )
      deallocate ( xg )
      deallocate ( uvecg, uvectg )
      deallocate ( uvecg_rhs )
      deallocate ( c )
      deallocate ( work3 )
      deallocate ( gi_up )

    end subroutine deallocate_arrays

  end subroutine interface_tracking_elem_lagrange


! Compute the following integrals on a closed interface (for 2D and 3D flow):
!
!           /        1  /            1  /
!   1)  V = | 1 dx = -  | div x dx = -  |  n.x ds
!           /        d  /     -      d  /  - -
!           V           V               S
!
!           /        1  /         2     1 /       2
!   2)  Q = | x dx = -  | grad |x| dx = - |  n |x| ds
!       -   / -      2  /       -       2 /  -  -
!           V           V                 S
!
! where d is the dimension of space (d=2 or 3), V is the area (d=2) or volume
! (d=3), Q is the first (linear) moment of area (d=2) or volume (d=3), and n
!        -                                                                 -
! is the outwardly directed unit normal on the interface.

! Note, that the center of area/volume can be computed by
!
!  x  = Q / V
!  -c   -
!
! after the integrations for V and Q have been performed.
!                                  -

  subroutine center_of_volume ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec )

    use interface_tracking_globals_m

!   input/output
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, ip, ics


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( intpolv > 0 ) then
        call set_shape_function ( shapefuncv, xig, phiv, dphiv )
      else
        phiv = phi
        dphiv = dphi
      end if

    end if

!   get coordinates of nodes of interface

    call get_coordinates ( mesh, elgrp, elem, x )

    if ( isoshape == 0 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, normal )

!     coordinates of integration points
      xg = matmul ( phi, x )

    else if ( isoshape == 1 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphiv, dxdxis, surfl, normal )

!     coordinates of integration points
      xg = matmul ( phiv, x )

    end if

!   compute integrands in integration points

    do ip = 1, ninti
      xn(ip) = dot_product(xg(ip,:),normal(ip,:))
      x2(ip) = dot_product(xg(ip,:),xg(ip,:))
    end do

    if ( coorsys == 1 ) then
       surfl = 2 * pi * xg(:,2) * surfl
       ics = 1
    else
       ics = 0
    end if

!   compute element contributions to the integral

!   volume V

    elemvec(1) = sum ( xn * wg * surfl ) / ( ndim + ics )

!   first moment Q

    do i = 1, ndim - ics
      elemvec(1+i) = sum ( x2 * normal(:,i) * wg * surfl ) / 2
    end do

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )
      allocate ( normal(ninti,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( xn(ninti), x2(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x )
      deallocate ( dxdxis )
      deallocate ( surfl )
      deallocate ( normal )
      deallocate ( xg )
      deallocate ( xn, x2 )

    end subroutine deallocate_arrays

  end subroutine center_of_volume


! Compute the following integrals on a closed interface (for 2D and 3D flow):
!
!           /        1  /            1  /
!   1)  V = | 1 dx = -  | div x dx = -  |  n.x ds
!           /        d  /     -      d  /  - -

!           /          1   /               1   /
!   2)  I = | xx dx =  -   | div xxx dx =  -   |  (n.x)xx ds
!       =   / --      d+2  /     ---      d+2  /   - - --
!           V              V                   S
!
! where d is the dimension of space (d=2 or 3), V is the area (d=2) or volume
! (d=3), I is the second moment of area (d=2) or volume (d=3), and n
!        =                                                         -
! is the outwardly directed unit normal on the interface.


  subroutine second_moment_of_area ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec )

    use interface_tracking_globals_m

!   input/output
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, k, ip


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interface_tracking ( mesh, coefficients, elgrp )

!     allocate arrays

      call allocate_arrays

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

      if ( intpolv > 0 ) then
        call set_shape_function ( shapefuncv, xig, phiv, dphiv )
      else
        phiv = phi
        dphiv = dphi
      end if

    end if

!   get coordinates of nodes of interface

    call get_coordinates ( mesh, elgrp, elem, x )

    if ( isoshape == 0 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, normal )

!     coordinates of integration points
      xg = matmul ( phi, x )

    else if ( isoshape == 1 ) then

!     surface geometry
      call isoparametric_deformation_curved ( x, dphiv, dxdxis, surfl, normal )

!     coordinates of integration points
      xg = matmul ( phiv, x )

    end if

!   compute integrands in integration points

    do ip = 1, ninti
      xn(ip) = dot_product(xg(ip,:),normal(ip,:))
    end do

    if ( coorsys == 1 ) then
      write(*,'(/2(a/))') 'Error second_moment_of_area:', &
        ' axisymmetric case has not been implemented. '
      stop
    end if

!   compute element contributions to the integral

!   volume V

    elemvec(1) = sum ( xn * wg * surfl ) / ndim

!   second moment I

    k = 1
    do i = 1, ndim
      do j= i, ndim
        elemvec(k+1) = &
                    sum ( xn * xg(:,i) * xg(:,j) * wg * surfl ) / ( ndim + 2 )
        k = k + 1
      end do
    end do

    if ( last ) then

!     deallocate arrays

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays

      allocate ( xig(ninti,ndim-1), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim-1) )
      allocate ( phiv(ninti,ndfv), dphiv(ninti,ndfv,ndim-1) )
      allocate ( x(nodalp,ndim) )
      allocate ( dxdxis(ninti,ndim,ndim-1) )
      allocate ( surfl(ninti) )
      allocate ( normal(ninti,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( xn(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( phiv, dphiv )
      deallocate ( x )
      deallocate ( dxdxis )
      deallocate ( surfl )
      deallocate ( normal )
      deallocate ( xg )
      deallocate ( xn )

    end subroutine deallocate_arrays

  end subroutine second_moment_of_area


! set global parameters

  subroutine set_globals_interface_tracking ( mesh, coefficients, elgrp )

    use interface_tracking_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

    integer :: ndimr

!   check size of coefficients

    call check ( coefficients, 'set_globals_interface_tracking', ncoefi=100, &
      ncoefr=50, indexarray=[13,16], minimum=[0,0], maximum=[1,1] )

    ndim = mesh%element(elgrp)%ndim
    nodalp = mesh%element(elgrp)%numnod
    globalshape = mesh%element(elgrp)%globalshape
    if ( ndim == 2 ) then
      coorsys = coefficients%i(13)
    else
      coorsys = 2
    end if
    isoshape = coefficients%i(16)

!   set number of degrees of freedom for position

    intpol = coefficients%i(1)

    ndimr = mesh%element(elgrp)%ndimr
    if ( intpol == 13 .and. any(mesh%element(elgrp)%p(:ndimr,2) /= 1 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_interface_tracking:', &
        ' For spectral elements GLL nodal distribution is required. '
      stop
    else if ( intpol == 20 .and. &
                            any(mesh%element(elgrp)%p(:ndimr,2) /= 0 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_interface_tracking:', &
        ' For high-order elements equidistant nodal distribution is required. '
      stop
    end if

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(2)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, 'set_globals_interface_tracking', ndf=ndf )

!   set number of degrees of freedom for velocity

    intpolv = coefficients%i(14)

    if ( intpolv > 0 ) then

      shapefuncv%globalshape = globalshape
      shapefuncv%interpolation = intpolv
      shapefuncv%numbering = 'standard'
      shapefuncv%p = coefficients%i(15)
      shapefuncv%spec_eval = 'gauss'

      call set_ndf ( shapefuncv, 'set_globals_interface_tracking', ndf=ndfv )

    else

     shapefuncv = shapefunc
     ndfv = ndf

    end if

!   set integration

    inttype = coefficients%i(5)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_interface_tracking:', &
        ' For spectral elements Gauss-Legendre-Lobatto integration ', &
        ' needs to be specified. '
      stop
    end if

    if ( coefficients%i(17) == 1 .or. SET_GAUSS_BY_ORDER ) then
      intrule = set_intrule ( globalshape, inttype, order=coefficients%i(3) )
    else
      intrule = coefficients%i(3)
    end if
    nsubint = get_coefficient ( coefficients, index=4, default=1 )

    gauss%globalshape = globalshape
    gauss%intrule = intrule
    gauss%nsubint = nsubint
    gauss%inttype = inttype

    call set_ninti ( gauss, ninti )

    if ( intpol == 13 .and. ninti /= ndf ) then
      write(*,'(/5(a/))') 'Error in set_globals_interface_tracking:', &
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_interface_tracking

end module interface_tracking_elements_m
