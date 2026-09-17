
! Copyright (C) 2019-2020 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for 3D beam structures

module structures_elements_beam3_m

  use tfem_elem_m
  use shapefunc_hermite_m

  implicit none


contains


! Internal element routine for a beam element with a linear elastic
! material law (Hooke's law). Geometrically linear. Two-dimensional.

  subroutine beam_elem3 ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use structures_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j
    real(dp) :: h, Emod, Gmod


    if ( first ) then

!     first element in this group

      call set_globals_beam3 ( mesh, coefficients, elgrp )

      if ( ndim /= 3 ) then
        write(*,'(/3(a/))') 'Error in beam_elem3:', &
        ' Mesh must be three-dimensional.'
        stop
      end if

!     allocate arrays

      allocate ( pos(2), pos1(ndim), pos2(ndim), por1(ndim), por2(ndim) )
      allocate ( posu(ndfu), posv(ndfv), posw(ndfv), post(ndft) )
      allocate ( e(ndim), nz(ndim), ny(ndim), Q(ndim,ndim), Qr(ndim,ndim) )
      allocate ( xig(ninti), wg(ninti) )
      allocate ( phi(ninti,ndfu), dphi_r(ninti,ndfu), dphi(ninti,ndfu) )
      allocate ( psi_r(ninti,ndfv), d2psi_r(ninti,ndfv) )
      allocate ( psi(ninti,ndfv), d2psi(ninti,ndfv) )
      allocate ( zeta(ninti,ndft), dzeta_r(ninti,ndft), dzeta(ninti,ndft) )
      allocate ( x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( qvecg(ninti,ndim), Ac(ninti) )
      allocate ( Izz(ninti), Iyz(ninti), Iyy(ninti), Jrr(ninti) )
      allocate ( fu(ndfu), fv(ndfv), fw(ndfv) )
      allocate ( qu(ninti), qv(ninti), qw(ninti) )
      allocate ( Suu(ndfu,ndfu), Stt(ndft,ndft) )
      allocate ( Svv(ndfv,ndfv), Svw(ndfv,ndfv), Sww(ndfv,ndfv) )

!     Gauss rule and shape functions

      call Gauss_Legendre_numeric_line ( ninti, xig, wg )

      select case ( intpolu )
      case(2)
        call shape_line_P1 ( xig, phi, dphi_r )
      case(6)
        call shape_line_P2 ( xig, phi, dphi_r )
      case default
        call errormsg_case_default ( 'beam_elem3', 'intpolu', &
          int_value=intpolu )
      end select

      select case ( intpolv )
      case(0)
        call shape_line_hermite_P3 ( xig, psi_r, d2phi=d2psi_r )
        pos = [ 2, 4 ]
      case(1)
        call shape_line_hermite_P4 ( xig, psi_r, d2phi=d2psi_r )
        pos = [ 2, 5 ]
      case default
        call errormsg_case_default ( 'beam_elem3', 'intpolv', &
          int_value=intpolv )
      end select

      select case ( intpolt )
      case(2)
        call shape_line_P1 ( xig, zeta, dzeta_r )
      case(6)
        call shape_line_P2 ( xig, zeta, dzeta_r )
      case default
        call errormsg_case_default ( 'beam_elem3', 'intpolt', &
          int_value=intpolt )
      end select

!     positions of degrees of freedom in element vector, matrix

      call set_pos_beam3

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the beam element

    e = x(nodalp,:) - x(1,:)
    h = sqrt(dot_product(e,e))
    e = e / h
    nz = coefficients%r(8:10)
    nz = nz / sqrt(dot_product(nz,nz))
    ny = cross_product( nz, e )
    Q(1,:) = e
    Q(2,:) = ny
    Q(3,:) = nz
    Qr(1,:) = nz
    Qr(2,:) = -ny
    Qr(3,:) = e

    do i = 1, ndim
      xg(:,i) = x(1,i) * ( 1 - xig ) / 2 + x(nodalp,i) * ( 1 + xig ) / 2
    end do

!   cross-sectional area

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(8), &
      value=coefficients%r(1), func=coefficients%func1(2)%p, &
      funcnr=coefficients%i(9), x=xg, coef=Ac )

!   second moment of area Izz

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(6), &
      value=coefficients%r(3), func=coefficients%func1(1)%p, &
      funcnr=coefficients%i(7), x=xg, coef=Izz )

!   second moment of area Iyy

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(14), &
      value=coefficients%r(11), func=coefficients%func1(3)%p, &
      funcnr=coefficients%i(15), x=xg, coef=Iyy )

!   product moment of area Iyz

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(16), &
      value=coefficients%r(12), func=coefficients%func1(4)%p, &
      funcnr=coefficients%i(17), x=xg, coef=Iyz )

!   polar second moment of area J

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(18), &
      value=coefficients%r(13), func=coefficients%func1(5)%p, &
      funcnr=coefficients%i(19), x=xg, coef=Jrr )

!   material properties

    Emod = coefficients%r(2)
    Gmod = coefficients%r(14)

!   change slope degrees from reference to real nodal values

    psi = psi_r
    psi(:,pos) = h * psi_r(:,pos) / 2
    d2psi = d2psi_r
    d2psi(:,pos) = h * d2psi(:,pos) / 2

!   change reference derivative to real derivative

    dphi = dphi_r * 2 / h
    d2psi = d2psi * 4 / h ** 2
    dzeta = dzeta_r * 2 / h

    if ( vector ) then

!     evaluate right-hand side

      if ( coefficients%i(10) >=0 ) then

        call evaluate_vector_coefficient_simple ( choice=coefficients%i(10), &
          value=coefficients%r(5:4+ndim), vfunc=coefficients%vfunc, &
          vfuncnr=coefficients%i(11), x=xg, coef=qvecg )

        qu = matmul(qvecg,e)
        qv = matmul(qvecg,ny)
        qw = matmul(qvecg,nz)

        do i = 1, ndfu
          fu(i) = sum ( qu * phi(:,i) * wg ) * h / 2
        end do
        do i = 1, ndfv
          fv(i) = sum ( qv * psi(:,i) * wg ) * h / 2
          fw(i) = sum ( qw * psi(:,i) * wg ) * h / 2
        end do

      else

        fu = 0; fv = 0; fw = 0

      end if

      elemvec(posu) = fu
      elemvec(posv) = fv
      elemvec(posw) = fw
      elemvec(post) = 0

!     transform end node displacements to global degrees

      elemvec(pos1) = matmul(transpose(Q),elemvec(pos1))
      elemvec(pos2) = matmul(transpose(Q),elemvec(pos2))

!     transform end node rotations to global degrees

      elemvec(por1) = matmul(transpose(Qr),elemvec(por1))
      elemvec(por2) = matmul(transpose(Qr),elemvec(por2))

    end if

    if ( matrix ) then

!     stiffness matrix

      do i = 1, ndfu
        do j = i, ndfu
          Suu(i,j) = sum ( Ac * dphi(:,i) * dphi(:,j) * wg )
          Suu(j,i) = Suu(i,j) ! symmetry
        end do
      end do

      Suu = Emod * Suu * h / 2

      do i = 1, ndfv
        do j = i, ndfv
          Svv(i,j) = sum ( Izz * d2psi(:,i) * d2psi(:,j) * wg )
          Svv(j,i) = Svv(i,j) ! symmetry
        end do
      end do

      Svv = Emod * Svv * h / 2

      do i = 1, ndfv
        do j = i, ndfv
          Svw(i,j) = sum ( Iyz * d2psi(:,i) * d2psi(:,j) * wg )
          Svw(j,i) = Svw(i,j) ! symmetry
        end do
      end do

      Svw = Emod * Svw * h / 2

      do i = 1, ndfv
        do j = i, ndfv
          Sww(i,j) = sum ( Iyy * d2psi(:,i) * d2psi(:,j) * wg )
          Sww(j,i) = Sww(i,j) ! symmetry
        end do
      end do

      Sww = Emod * Sww * h / 2

      do i = 1, ndft
        do j = i, ndft
          Stt(i,j) = sum ( Jrr * dzeta(:,i) * dzeta(:,j) * wg )
          Stt(j,i) = Stt(i,j) ! symmetry
        end do
      end do

      Stt = Gmod * Stt * h / 2

      elemmat = 0
      elemmat(posu,posu) = Suu
      elemmat(posv,posv) = Svv
      elemmat(posv,posw) = Svw
      elemmat(posw,posv) = Svw
      elemmat(posw,posw) = Sww
      elemmat(post,post) = Stt

!     transform end node displacements to global degrees

      elemmat(pos1,:) = matmul(transpose(Q),elemmat(pos1,:))
      elemmat(:,pos1) = matmul(elemmat(:,pos1),Q)
      elemmat(pos2,:) = matmul(transpose(Q),elemmat(pos2,:))
      elemmat(:,pos2) = matmul(elemmat(:,pos2),Q)

!     transform end node rotations to global degrees

      elemmat(por1,:) = matmul(transpose(Qr),elemmat(por1,:))
      elemmat(:,por1) = matmul(elemmat(:,por1),Qr)
      elemmat(por2,:) = matmul(transpose(Qr),elemmat(por2,:))
      elemmat(:,por2) = matmul(elemmat(:,por2),Qr)

    end if

    if ( last ) then

!     last element in this group

      deallocate ( pos, pos1, pos2, por1, por2 )
      deallocate ( posu, posv, posw, post )
      deallocate ( e, nz, ny, Q, Qr )
      deallocate ( xig, wg )
      deallocate ( phi, dphi_r, dphi )
      deallocate ( psi_r, d2psi_r )
      deallocate ( psi, d2psi )
      deallocate ( zeta, dzeta_r, dzeta )
      deallocate ( x, xg )
      deallocate ( qvecg, Ac )
      deallocate ( Izz, Iyz, Iyy, Jrr )
      deallocate ( fu, fv, fw )
      deallocate ( qu, qv, qw )
      deallocate ( Suu, Stt )
      deallocate ( Svv, Svw, Sww )

    end if

  end subroutine beam_elem3


! Compute derivatives for a beam element with a linear elastic
! material law (Hooke's law). Geometrically linear. Three-dimensional.

  subroutine beam_deriv3 ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use structures_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    integer :: comp
    real(dp) :: h, Emod, Gmod


    if ( first ) then

!     first element in this group

      call set_globals_beam3 ( mesh, coefficients, elgrp )

!     check component = coefficients%i(1)

      call check ( coefficients, 'beam_deriv3', indexarray=[1], &
        minimum=[1], maximum=[16] )

      allocate ( pos(2), pos1(ndim), pos2(ndim), por1(ndim), por2(ndim) )
      allocate ( posu(ndfu), posv(ndfv), posw(ndfv), post(ndft) )
      allocate ( e(ndim), nz(ndim), ny(ndim), Q(ndim,ndim), Qr(ndim,ndim) )
      allocate ( xrnod(nodalp,1) )
      allocate ( phi(nodalp,ndfu), dphi_r(nodalp,ndfu), dphi(nodalp,ndfu) )
      allocate ( psi_r(nodalp,ndfv), dpsi_r(nodalp,ndfv) )
      allocate ( d2psi_r(nodalp,ndfv), d3psi_r(nodalp,ndfv) )
      allocate ( psi(nodalp,ndfv), dpsi(nodalp,ndfv) )
      allocate ( d2psi(nodalp,ndfv), d3psi(nodalp,ndfv) )
      allocate ( zeta(nodalp,ndft), dzeta_r(nodalp,ndft), dzeta(nodalp,ndft) )
      allocate ( x(nodalp,ndim) )
      allocate ( Ac(nodalp), Izz(nodalp), dIzz(nodalp) )
      allocate ( Iyz(nodalp), dIyz(nodalp), Iyy(nodalp), dIyy(nodalp) )
      allocate ( Jrr(nodalp), sol(ndfu+2*ndfv+ndft) )
      allocate ( u(ndfu), v(ndfv), w(ndfv), t(ndft) )

!     set reference coordinates and shape function

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      select case ( intpolu )
      case(2)
        call shape_line_P1 ( xrnod(:,1), phi, dphi_r )
      case(6)
        call shape_line_P2 ( xrnod(:,1), phi, dphi_r )
      case default
        call errormsg_case_default ( 'beam_deriv3', 'intpolu', &
          int_value=intpolu )
      end select

      select case ( intpolv )
      case(0)
        call shape_line_hermite_P3 ( xrnod(:,1), psi_r, dpsi_r, d2psi_r, &
          d3psi_r )
        pos = [ 2, 4 ]
      case(1)
        call shape_line_hermite_P4 ( xrnod(:,1), psi_r, dpsi_r, d2psi_r, &
          d3psi_r )
        pos = [ 2, 5 ]
      case default
        call errormsg_case_default ( 'beam_deriv3', 'intpolv', &
          int_value=intpolv )
      end select

      select case ( intpolt )
      case(2)
        call shape_line_P1 ( xrnod(:,1), zeta, dzeta_r )
      case(6)
        call shape_line_P2 ( xrnod(:,1), zeta, dzeta_r )
      case default
        call errormsg_case_default ( 'beam_deriv3', 'intpolt', &
          int_value=intpolt )
      end select

!     positions of degrees of freedom in element vector, matrix

      call set_pos_beam3

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the beam element

    e = x(nodalp,:) - x(1,:)
    h = sqrt(dot_product(e,e))
    e = e / h
    nz = coefficients%r(8:10)
    nz = nz / sqrt(dot_product(nz,nz))
    ny = cross_product( nz, e )
    Q(1,:) = e
    Q(2,:) = ny
    Q(3,:) = nz
    Qr(1,:) = nz
    Qr(2,:) = -ny
    Qr(3,:) = e

!   cross-sectional area

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(8), &
      value=coefficients%r(1), func=coefficients%func1(2)%p, &
      funcnr=coefficients%i(9), x=x, coef=Ac )

!   second moment of area Izz

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(6), &
      value=coefficients%r(3), func=coefficients%func1(1)%p, &
      funcnr=coefficients%i(7), x=x, coef=Izz )

!   second moment of area Iyy

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(14), &
      value=coefficients%r(11), func=coefficients%func1(3)%p, &
      funcnr=coefficients%i(15), x=x, coef=Iyy )

!   product moment of area Iyz

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(16), &
      value=coefficients%r(12), func=coefficients%func1(4)%p, &
      funcnr=coefficients%i(17), x=x, coef=Iyz )

!   polar second moment of area J

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(18), &
      value=coefficients%r(13), func=coefficients%func1(5)%p, &
      funcnr=coefficients%i(19), x=x, coef=Jrr )

!   material properties

    Emod = coefficients%r(2)
    Gmod = coefficients%r(14)

!   change slope degrees from reference to real nodal values

    psi = psi_r
    psi(:,pos) = h * psi_r(:,pos) / 2
    dpsi = dpsi_r
    dpsi(:,pos) = h * dpsi_r(:,pos) / 2
    d2psi = d2psi_r
    d2psi(:,pos) = h * d2psi_r(:,pos) / 2
    d3psi = d3psi_r
    d3psi(:,pos) = h * d3psi_r(:,pos) / 2

!   change reference derivative to real derivative

    dphi  = dphi_r * 2 / h
    dpsi  = dpsi * 2 / h
    d2psi = d2psi * 4 / h ** 2
    d3psi = d3psi * 8 / h ** 3
    dzeta = dzeta_r * 2 / h

!   get solution

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, sol )

!   transform end node displacements/rotations to local degrees

    sol(pos1) = matmul(Q,sol(pos1))
    sol(pos2) = matmul(Q,sol(pos2))
    sol(por1) = matmul(Qr,sol(por1))
    sol(por2) = matmul(Qr,sol(por2))

!   extract axial and transverse displacements

    u = sol(posu)
    v = sol(posv)
    w = sol(posw)
    t = sol(post)

!   compute derivative

    comp = coefficients%i(1)

    select case(comp)
    case(1)
      elemvec = matmul(psi,v) ! transverse displacement v
    case(2)
      elemvec = matmul(dpsi,v) ! slope v'
    case(3)
      elemvec = matmul(d2psi,v) ! second derivative v''
    case(4)
!     bending moment Mz
      elemvec = Emod * ( Iyz * matmul(d2psi,w) + Izz * matmul(d2psi,v) )
    case(5)
!     shear force V
      dIzz = (Izz(nodalp)-Izz(1)) / h
      dIyz = (Iyz(nodalp)-Iyz(1)) / h
      elemvec = Emod * ( Izz * matmul(d3psi,v) + dIzz * matmul(d2psi,v) + &
                         Iyz * matmul(d3psi,w) + dIyz * matmul(d2psi,w) )
    case(6)
      elemvec = matmul(phi,u) ! axial displacement
    case(7)
      elemvec = matmul(dphi,u) ! axial strain
    case(8)
      elemvec = Emod * Ac * matmul(dphi,u) ! axial force
    case(9)
      elemvec = matmul(zeta,t) ! axial (torsional) rotation
    case(10)
      elemvec = matmul(dzeta,t) ! torsional strain
    case(11)
      elemvec = Gmod * Jrr * matmul(dzeta,t) ! torsional moment Mx
    case(12)
      elemvec = matmul(psi,w) ! transverse displacement w
    case(13)
      elemvec = matmul(dpsi,w) ! slope w'
    case(14)
      elemvec = matmul(d2psi,w) ! second derivative w''
    case(15)
!     bending moment My
      elemvec = - Emod * ( Iyy * matmul(d2psi,w) + Iyz * matmul(d2psi,v) )
    case(16)
!     shear force W
      dIyy = (Iyy(nodalp)-Iyy(1)) / h
      dIyz = (Iyz(nodalp)-Iyz(1)) / h
      elemvec = - Emod * ( Iyy * matmul(d3psi,w) + dIyy * matmul(d2psi,w) + &
                           Iyz * matmul(d3psi,v) + dIyz * matmul(d2psi,v) )
    case default
      call errormsg_case_default ( 'beam_deriv3', 'comp', int_value=comp )
    end select

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( pos, pos1, pos2, por1, por2 )
      deallocate ( posu, posv, posw, post )
      deallocate ( e, nz, ny, Q, Qr )
      deallocate ( xrnod )
      deallocate ( phi, dphi_r, dphi )
      deallocate ( psi_r, dpsi_r )
      deallocate ( d2psi_r, d3psi_r )
      deallocate ( psi, dpsi )
      deallocate ( d2psi, d3psi )
      deallocate ( zeta, dzeta_r, dzeta )
      deallocate ( x )
      deallocate ( Ac, Izz, dIzz )
      deallocate ( Iyz, dIyz, Iyy, dIyy )
      deallocate ( Jrr, sol )
      deallocate ( u, v, w, t )

    end if

  end subroutine beam_deriv3


! Displacement for a beam element in all nodes. Three-dimensional.
! Needed for elements with center node.

  subroutine beam_displacement3 ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use structures_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    real(dp) :: h

    if ( first ) then

!     first element in this group

      call set_globals_beam3 ( mesh, coefficients, elgrp )

      allocate ( pos(2), pos1(ndim), pos2(ndim), por1(ndim), por2(ndim) )
      allocate ( posu(ndfu), posv(ndfv), posw(ndfv), post(ndft) )
      allocate ( e(ndim), nz(ndim), ny(ndim), Q(ndim,ndim), Qr(ndim,ndim) )
      allocate ( xrnod(nodalp,1) )
      allocate ( phi(nodalp,ndfu) )
      allocate ( psi_r(nodalp,ndfv) )
      allocate ( psi(nodalp,ndfv) )
      allocate ( x(nodalp,ndim) )
      allocate ( sol(ndfu+2*ndfv+ndft) )
      allocate ( u(ndfu), v(ndfv), w(ndfv) )
      allocate ( work(ndim,nodalp) )

!     set reference coordinates and shape function

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      select case ( intpolu )
      case(2)
        call shape_line_P1 ( xrnod(:,1), phi )
      case(6)
        call shape_line_P2 ( xrnod(:,1), phi )
      case default
        call errormsg_case_default ( 'beam_displacement3', &
          'intpolu', int_value=intpolu )
      end select

      select case ( intpolv )
      case(0)
        call shape_line_hermite_P3 ( xrnod(:,1), psi_r )
        pos = [ 2, 4 ]
      case(1)
        call shape_line_hermite_P4 ( xrnod(:,1), psi_r )
        pos = [ 2, 5 ]
      case default
        call errormsg_case_default ( 'beam_displacement3', &
          'intpolv', int_value=intpolv )
      end select

!     positions of degrees of freedom in element vector, matrix

      call set_pos_beam3

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the beam element

    e = x(nodalp,:) - x(1,:)
    h = sqrt(dot_product(e,e))
    e = e / h
    nz = coefficients%r(8:10)
    nz = nz / sqrt(dot_product(nz,nz))
    ny = cross_product( nz, e )
    Q(1,:) = e
    Q(2,:) = ny
    Q(3,:) = nz
    Qr(1,:) = nz
    Qr(2,:) = -ny
    Qr(3,:) = e

!   change slope degrees from reference to real nodal values

    psi = psi_r
    psi(:,pos) = h * psi_r(:,pos) / 2

!   get solution

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, sol )

!   transform end node displacements to local degrees

    sol(pos1) = matmul(Q,sol(pos1))
    sol(pos2) = matmul(Q,sol(pos2))

!   transform end node rotations to local degrees

    sol(por1) = matmul(Qr,sol(por1))
    sol(por2) = matmul(Qr,sol(por2))

!   extract axial and transverse displacements

    u = sol(posu)
    v = sol(posv)
    w = sol(posw)

!   compute displacements in all nodes

    work(1,:) = matmul(phi,u)
    work(2,:) = matmul(psi,v)
    work(3,:) = matmul(psi,w)

    work = matmul( transpose(Q), work )

    elemvec = reshape( work, [ndim*nodalp] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( pos, pos1, pos2, por1, por2 )
      deallocate ( posu, posv, posw, post )
      deallocate ( e, nz, ny, Q, Qr )
      deallocate ( xrnod )
      deallocate ( phi )
      deallocate ( psi_r )
      deallocate ( psi )
      deallocate ( x )
      deallocate ( sol )
      deallocate ( u, v, w )
      deallocate ( work )

    end if

  end subroutine beam_displacement3


! Rotations for a beam element in all nodes. Three-dimensional.
! Needed for elements with center node.

  subroutine beam_rotation3 ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use structures_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    real(dp) :: h


    if ( first ) then

!     first element in this group

      call set_globals_beam3 ( mesh, coefficients, elgrp )

      allocate ( pos(2), pos1(ndim), pos2(ndim), por1(ndim), por2(ndim) )
      allocate ( posu(ndfu), posv(ndfv), posw(ndfv), post(ndft) )
      allocate ( e(ndim), nz(ndim), ny(ndim), Q(ndim,ndim), Qr(ndim,ndim) )
      allocate ( xrnod(nodalp,1) )
      allocate ( psi_r(nodalp,ndfv), dpsi_r(nodalp,ndfv) )
      allocate ( psi(nodalp,ndfv), dpsi(nodalp,ndfv) )
      allocate ( zeta(nodalp,ndft) )
      allocate ( x(nodalp,ndim) )
      allocate ( sol(ndfu+2*ndfv+ndft) )
      allocate ( v(ndfv), w(ndfv), t(ndft) )
      allocate ( work(ndim,nodalp) )

!     set reference coordinates and shape function

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      select case ( intpolv )
      case(0)
        call shape_line_hermite_P3 ( xrnod(:,1), psi_r, dpsi_r )
        pos = [ 2, 4 ]
      case(1)
        call shape_line_hermite_P4 ( xrnod(:,1), psi_r, dpsi_r )
        pos = [ 2, 5 ]
      case default
        call errormsg_case_default ( 'beam_rotation3', &
          'intpolv', int_value=intpolv )
      end select

      select case ( intpolt )
      case(2)
        call shape_line_P1 ( xrnod(:,1), zeta, dzeta_r )
      case(6)
        call shape_line_P2 ( xrnod(:,1), zeta, dzeta_r )
      case default
        call errormsg_case_default ( 'beam_rotation3', &
          'intpolt', int_value=intpolt )
      end select

!     positions of degrees of freedom in element vector, matrix

      call set_pos_beam3

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the beam element

    e = x(nodalp,:) - x(1,:)
    h = sqrt(dot_product(e,e))
    e = e / h
    nz = coefficients%r(8:10)
    nz = nz / sqrt(dot_product(nz,nz))
    ny = cross_product( nz, e )
    Q(1,:) = e
    Q(2,:) = ny
    Q(3,:) = nz
    Qr(1,:) = nz
    Qr(2,:) = -ny
    Qr(3,:) = e

!   change slope degrees from reference to real nodal values

    psi = psi_r
    psi(:,pos) = h * psi_r(:,pos) / 2
    dpsi = dpsi_r
    dpsi(:,pos) = h * dpsi_r(:,pos) / 2

!   change reference derivative to real derivative

    dpsi  = dpsi * 2 / h

!   get solution

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, sol )

!   transform end node displacements/rotations to local degrees

    sol(pos1) = matmul(Q,sol(pos1))
    sol(pos2) = matmul(Q,sol(pos2))
    sol(por1) = matmul(Qr,sol(por1))
    sol(por2) = matmul(Qr,sol(por2))

!   extract transverse displacements and torsion

    v = sol(posv)
    w = sol(posw)
    t = sol(post)

!   compute rotations in all nodes with respect to the global system

    work(1,:) = matmul(dpsi,v) ! slope v'
    work(2,:) = matmul(dpsi,w) ! slope w'
    work(3,:) = matmul(zeta,t) ! axial (torsional) rotation

    work = matmul( transpose(Qr), work )

    elemvec = reshape( work, [ndim*nodalp] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( pos, pos1, pos2, por1, por2 )
      deallocate ( posu, posv, posw, post )
      deallocate ( e, nz, ny, Q, Qr )
      deallocate ( xrnod )
      deallocate ( psi_r, dpsi_r )
      deallocate ( psi, dpsi )
      deallocate ( zeta )
      deallocate ( x )
      deallocate ( sol )
      deallocate ( v, w, t )
      deallocate ( work )

    end if

  end subroutine beam_rotation3


! Distributed load for a beam element in all nodes. Three-dimensional beam.
! For post processing purposes only.

  subroutine beam_q3 ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use structures_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    if ( first ) then

!     first element in this group

      call set_globals_beam3 ( mesh, coefficients, elgrp )

      allocate ( x(nodalp,ndim) )
      allocate ( qvecnod(nodalp,ndim) )

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   compute distributed loading

    if ( coefficients%i(10) >=0 ) then

      call evaluate_vector_coefficient_simple ( choice=coefficients%i(10), &
        value=coefficients%r(5:4+ndim), vfunc=coefficients%vfunc, &
        vfuncnr=coefficients%i(11), x=x, coef=qvecnod )

    else

      qvecnod = 0

    end if

!   set q vector in all nodes

    elemvec = reshape( transpose(qvecnod), [ndim*nodalp] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( x )
      deallocate ( qvecnod )

    end if

  end subroutine beam_q3


! set pos arrays for a 2D beam element

  subroutine set_pos_beam3

    use structures_globals_m

!   positions of degrees of freedom in element vector, matrix

    pos1 = [ 1, 2, 3 ]  ! displacements in first node
    por1 = [ 4, 5, 6 ]  ! rotations in first node

    select case ( intpolu )
    case(2)
      select case ( intpolt )
      case(2)
        select case ( intpolv )
        case(0)
          pos2 = [ 7, 8, 9 ]          ! displacements in last node
          por2 = [ 10, 11, 12 ]       ! rotations in last node
          posu = [ 1, 7 ]             ! axial displacement degrees of freedom
          posv = [ 2, 4, 8, 10 ]      ! transverse v displacement degs
          posw = [ 3, 5, 9, 11 ]      ! transverse w displacement degs
          post = [ 6, 12 ]            ! axial rotation degs
        case(1)
          pos2 = [ 9, 10, 11 ]        ! displacements in last node
          por2 = [ 12, 13, 14 ]       ! rotations in last node
          posu = [ 1, 9 ]             ! axial displacement degrees of freedom
          posv = [ 2, 4, 7, 10, 12 ]  ! transverse v displacement degs
          posw = [ 3, 5, 8, 11, 13 ]  ! transverse w displacement degs
          post = [ 6, 14 ]            ! axial rotation degs
        case default
          call errormsg_case_default ( 'set_pos_beam3', &
            'intpolv', int_value=intpolv )
        end select
      case(6)
        select case ( intpolv )
        case(0)
          pos2 = [ 8, 9, 10 ]         ! displacements in last node
          por2 = [ 11, 12, 13 ]       ! rotations in last node
          posu = [ 1, 8 ]             ! axial displacement degrees of freedom
          posv = [ 2, 4, 9, 11 ]      ! transverse v displacement degs
          posw = [ 3, 5, 10, 12 ]     ! transverse w displacement degs
          post = [ 6, 7, 13 ]         ! axial rotation degs
        case(1)
          pos2 = [ 10, 11, 12 ]       ! displacements in last node
          por2 = [ 13, 14, 15 ]       ! rotations in last node
          posu = [ 1, 10 ]            ! axial displacement degrees of freedom
          posv = [ 2, 4, 7, 11, 13 ]  ! transverse v displacement degs
          posw = [ 3, 5, 8, 12, 14 ]  ! transverse w displacement degs
          post = [ 6, 9, 15 ]         ! axial rotation degs
        case default
          call errormsg_case_default ( 'set_pos_beam3', &
            'intpolv', int_value=intpolv )
        end select
      case default
        call errormsg_case_default ( 'set_pos_beam3', &
          'intpolt', int_value=intpolt )
      end select
    case(6)
      select case ( intpolt )
      case(2)
        select case ( intpolv )
        case(0)
          pos2 = [ 8, 9, 10 ]        ! displacements in last node
          por2 = [ 11, 12, 13 ]      ! rotations in last node
          posu = [ 1, 7, 8 ]         ! axial displacement degrees of freedom
          posv = [ 2, 4, 9, 11 ]     ! transverse v displacement degs
          posw = [ 3, 5, 10, 12 ]    ! transverse w displacement degs
          post = [ 6, 13 ]           ! axial rotation degs
        case(1)
          pos2 = [ 10, 11, 12 ]      ! displacements in last node
          por2 = [ 13, 14, 15 ]      ! rotations in last node
          posu = [ 1, 7, 10 ]        ! axial displacement degrees of freedom
          posv = [ 2, 4, 8, 11, 13 ] ! transverse v displacement degs
          posw = [ 3, 5, 9, 12, 14 ] ! transverse w displacement degs
          post = [ 6, 15 ]           ! axial rotation degs
        case default
          call errormsg_case_default ( 'set_pos_beam3', &
            'intpolv', int_value=intpolv )
        end select
      case(6)
        select case ( intpolv )
        case(0)
          pos2 = [ 9, 10, 11 ]       ! displacements in last node
          por2 = [ 12, 13, 14 ]      ! rotations in last node
          posu = [ 1, 7, 9 ]         ! axial displacement degrees of freedom
          posv = [ 2, 4, 10, 12 ]    ! transverse v displacement degs
          posw = [ 3, 5, 11, 13 ]    ! transverse w displacement degs
          post = [ 6, 8, 14 ]        ! axial rotation degs
        case(1)
          pos2 = [ 11, 12, 13 ]      ! displacements in last node
          por2 = [ 14, 15, 16 ]      ! rotations in last node
          posu = [ 1, 7, 11 ]        ! axial displacement degrees of freedom
          posv = [ 2, 4, 8, 12, 14 ] ! transverse v displacement degs
          posw = [ 3, 5, 9, 13, 15 ] ! transverse w displacement degs
          post = [ 6, 10, 16 ]       ! axial rotation degs
        case default
          call errormsg_case_default ( 'set_pos_beam3', &
            'intpolv', int_value=intpolv )
        end select
      case default
        call errormsg_case_default ( 'set_pos_beam3', &
          'intpolt', int_value=intpolt )
      end select
    case default
      call errormsg_case_default ( 'set_pos_beam3', &
        'intpolu', int_value=intpolu )
    end select

  end subroutine set_pos_beam3


! set global parameters for a beam element (internal element)

  subroutine set_globals_beam3 ( mesh, coefficients, elgrp )

    use structures_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

!   check size of coefficients

    call check ( coefficients, 'set_globals_beam', ncoefi=100, &
      ncoefr=50, indexarray=[3,12,13], minimum=[0,0,0], maximum=[1,6,6] )

    ndim = mesh%element(elgrp)%ndim
    nodalp = mesh%element(elgrp)%numnod
    globalshape = mesh%element(elgrp)%globalshape

    if ( globalshape /= 'line' ) then
      write(*,'(/3(a/))') 'Error in set_globals_beam3:', &
      ' Only line elements are allowed'
      stop
    end if

!   shape functions for axial displacement

    intpolu = get_coefficient ( coefficients, index=12, default=2 )

    select case (intpolu)
    case(2)
!     P1
      ndfu = 2
    case(6)
!     P2
      ndfu = 3
    case default
      write(*,'(/a/a,i0/)') 'Error in set_globals_beam3:', &
        ' Invalid interpolation, intpolu = ', intpolu
      stop
    end select

!   shape functions for lateral displacement

    intpolv = coefficients%i(3)

    select case (intpolv)
    case(0)
!     P3 Hermite
      ndfv = 4
    case(1)
!     P4 Hermite
      ndfv = 5
    case default
      write(*,'(/a/a,i0/)') 'Error in set_globals_beam3:', &
        ' Invalid interpolation, intpolv = ', intpolv
      stop
    end select

!   shape functions for axial (torsion) rotation

    intpolt = get_coefficient ( coefficients, index=13, default=2 )

    select case (intpolt)
    case(2)
!     P1
      ndft = 2
    case(6)
!     P2
      ndft = 3
    case default
      write(*,'(/a/a,i0/)') 'Error in set_globals_beam3:', &
        ' Invalid interpolation, intpolt = ', intpolt
      stop
    end select

    ninti = coefficients%i(2)

!   check Z  vector to be non-zero

    if ( all( abs(coefficients%r(8:10)) < 10*tiny(1._dp) ) ) then
      write(*,'(/3(a/))') 'Error in set_globals_beam3:', &
      ' coefficients%r(8:10) too small and cannot define nz vector.'
      stop
    end if

  end subroutine set_globals_beam3

end module structures_elements_beam3_m

