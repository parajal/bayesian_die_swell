
! Copyright (C) 2019-2020 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for 2D beam structures

module structures_elements_beam2_m

  use tfem_elem_m
  use shapefunc_hermite_m

  implicit none


contains


! Internal element routine for a beam element with a linear elastic
! material law (Hooke's law). Geometrically linear. Two-dimensional.

  subroutine beam_elem2 ( mesh, problem, elgrp, elem, matrix, &
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
    real(dp) :: h, Emod


    if ( first ) then

!     first element in this group

      call set_globals_beam2 ( mesh, coefficients, elgrp )

      if ( ndim /= 2 ) then
        write(*,'(/3(a/))') 'Error in beam_elem2:', &
        ' Mesh must be two-dimensional.'
        stop
      end if

!     allocate arrays

      allocate ( pos(2), pos1(ndim), pos2(ndim), posu(ndfu), posv(ndfv) )
      allocate ( e(ndim), n(ndim), Q(ndim,ndim) )
      allocate ( xig(ninti), wg(ninti) )
      allocate ( phi(ninti,ndfu), dphi_r(ninti,ndfu), dphi(ninti,ndfu) )
      allocate ( psi_r(ninti,ndfv), d2psi_r(ninti,ndfv) )
      allocate ( psi(ninti,ndfv), d2psi(ninti,ndfv) )
      allocate ( x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( qvecg(ninti,ndim), Ac(ninti), Iz(ninti) )
      allocate ( fu(ndfu), fv(ndfv) )
      allocate ( qu(ninti), qv(ninti) )
      allocate ( Suu(ndfu,ndfu), Svv(ndfv,ndfv) )

!     Gauss rule and shape functions

      call Gauss_Legendre_numeric_line ( ninti, xig, wg )

      select case ( intpolu )
      case(2)
        call shape_line_P1 ( xig, phi, dphi_r )
      case(6)
        call shape_line_P2 ( xig, phi, dphi_r )
      case default
        call errormsg_case_default ( 'beam_elem2', &
          'intpolu', int_value=intpolu )
      end select

      select case ( intpolv )
      case(0)
        call shape_line_hermite_P3 ( xig, psi_r, d2phi=d2psi_r )
        pos = [ 2, 4 ]
      case(1)
        call shape_line_hermite_P4 ( xig, psi_r, d2phi=d2psi_r )
        pos = [ 2, 5 ]
      case default
        call errormsg_case_default ( 'beam_elem2', &
          'intpolv', int_value=intpolv )
      end select

!     positions of degrees of freedom in element vector, matrix

      call set_pos_beam2

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the beam element

    e = x(nodalp,:) - x(1,:)
    h = sqrt(dot_product(e,e))
    e = e / h
    n = [ -e(2), e(1) ]
    Q(1,:) = e
    Q(2,:) = n

    do i = 1, ndim
      xg(:,i) = x(1,i) * ( 1 - xig ) / 2 + x(nodalp,i) * ( 1 + xig ) / 2
    end do

!   cross-sectional area

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(8), &
      value=coefficients%r(1), func=coefficients%func1(2)%p, &
      funcnr=coefficients%i(9), x=xg, coef=Ac )

!   second moment of area

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(6), &
      value=coefficients%r(3), func=coefficients%func1(1)%p, &
      funcnr=coefficients%i(7), x=xg, coef=Iz )

!   material properties

    Emod = coefficients%r(2)

!   change slope degrees from reference to real nodal values

    psi = psi_r
    psi(:,pos) = h * psi_r(:,pos) / 2
    d2psi = d2psi_r
    d2psi(:,pos) = h * d2psi(:,pos) / 2

!   change reference derivative to real derivative

    dphi = dphi_r * 2 / h
    d2psi = d2psi * 4 / h ** 2

    if ( vector ) then

!     evaluate right-hand side

      if ( coefficients%i(10) >=0 ) then

        call evaluate_vector_coefficient_simple ( choice=coefficients%i(10), &
          value=coefficients%r(5:4+ndim), vfunc=coefficients%vfunc, &
          vfuncnr=coefficients%i(11), x=xg, coef=qvecg )

        qu = matmul(qvecg,e)
        qv = matmul(qvecg,n)

        do i = 1, ndfu
          fu(i) = sum ( qu * phi(:,i) * wg ) * h / 2
        end do
        do i = 1, ndfv
          fv(i) = sum ( qv * psi(:,i) * wg ) * h / 2
        end do

      else

        fu = 0; fv = 0

      end if

      elemvec(posu) = fu
      elemvec(posv) = fv

!     transform end node displacements to global degrees

      elemvec(pos1) = matmul(transpose(Q),elemvec(pos1))
      elemvec(pos2) = matmul(transpose(Q),elemvec(pos2))

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
          Svv(i,j) = sum ( Iz * d2psi(:,i) * d2psi(:,j) * wg )
          Svv(j,i) = Svv(i,j) ! symmetry
        end do
      end do

      Svv = Emod * Svv * h / 2

      elemmat(posu,posu) = Suu
      elemmat(posv,posv) = Svv
      elemmat(posu,posv) = 0
      elemmat(posv,posu) = 0

!     transform end node displacements to global degrees

      elemmat(pos1,:) = matmul(transpose(Q),elemmat(pos1,:))
      elemmat(:,pos1) = matmul(elemmat(:,pos1),Q)
      elemmat(pos2,:) = matmul(transpose(Q),elemmat(pos2,:))
      elemmat(:,pos2) = matmul(elemmat(:,pos2),Q)

    end if

    if ( last ) then

!     last element in this group

      deallocate ( pos, posu, pos1, pos2, posv )
      deallocate ( e, n, Q )
      deallocate ( xig, wg )
      deallocate ( phi, dphi_r, dphi )
      deallocate ( psi_r, d2psi_r )
      deallocate ( psi, d2psi )
      deallocate ( x, xg )
      deallocate ( qvecg, Ac, Iz )
      deallocate ( fu, fv )
      deallocate ( qu, qv )
      deallocate ( Suu, Svv )

    end if

  end subroutine beam_elem2


! Compute derivatives for a beam element with a linear elastic
! material law (Hooke's law). Geometrically linear. Two-dimensional.

  subroutine beam_deriv2 ( mesh, problem, elgrp, elem, first, last, &
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
    real(dp) :: h, Emod


    if ( first ) then

!     first element in this group

      call set_globals_beam2 ( mesh, coefficients, elgrp )

!     check component = coefficients%i(1)

      call check ( coefficients, 'beam_deriv2', indexarray=[1], &
        minimum=[1], maximum=[8] )

      allocate ( pos(2), pos1(ndim), pos2(ndim), posu(ndfu), posv(ndfv) )
      allocate ( e(ndim), n(ndim), Q(ndim,ndim) )
      allocate ( xrnod(nodalp,1) )
      allocate ( phi(nodalp,ndfu), dphi_r(nodalp,ndfu), dphi(nodalp,ndfu) )
      allocate ( psi_r(nodalp,ndfv), dpsi_r(nodalp,ndfv) )
      allocate ( d2psi_r(nodalp,ndfv), d3psi_r(nodalp,ndfv) )
      allocate ( psi(nodalp,ndfv), dpsi(nodalp,ndfv) )
      allocate ( d2psi(nodalp,ndfv), d3psi(nodalp,ndfv) )
      allocate ( x(nodalp,ndim) )
      allocate ( Ac(nodalp), Iz(nodalp), dIz(nodalp), sol(ndfu+ndfv) )
      allocate ( u(ndfu), v(ndfv) )

!     set reference coordinates and shape function

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      select case ( intpolu )
      case(2)
        call shape_line_P1 ( xrnod(:,1), phi, dphi_r )
      case(6)
        call shape_line_P2 ( xrnod(:,1), phi, dphi_r )
      case default
        call errormsg_case_default ( 'beam_deriv2', &
          'intpolu', int_value=intpolu )
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
        call errormsg_case_default ( 'beam_deriv2', &
          'intpolv', int_value=intpolv )
      end select

!     positions of degrees of freedom in element vector, matrix

      call set_pos_beam2

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the beam element

    e = x(nodalp,:) - x(1,:)
    h = sqrt(dot_product(e,e))
    e = e / h
    n = [ -e(2), e(1) ]
    Q(1,:) = e
    Q(2,:) = n

!   second moment of area

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(6), &
      value=coefficients%r(3), func=coefficients%func1(1)%p, &
      funcnr=coefficients%i(7), x=x, coef=Iz )

!   cross-sectional area

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(8), &
      value=coefficients%r(1), func=coefficients%func1(2)%p, &
      funcnr=coefficients%i(9), x=x, coef=Ac )

!   material properties

    Emod = coefficients%r(2)

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

!   get solution

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, sol )

!   transform end node displacements to local degrees

    sol(pos1) = matmul(Q,sol(pos1))
    sol(pos2) = matmul(Q,sol(pos2))

!   extract axial and transverse displacements

    u = sol(posu)
    v = sol(posv)

!   compute derivative

    comp = coefficients%i(1)

    select case(comp)
    case(1)
      elemvec = matmul(psi,v) ! transverse displacement
    case(2)
      elemvec = matmul(dpsi,v) ! slope
    case(3)
      elemvec = matmul(d2psi,v) ! second derivative
    case(4)
      elemvec = Emod * Iz * matmul(d2psi,v) ! bending moment
    case(5)
      dIz = (Iz(nodalp)-Iz(1)) / h
      elemvec = Emod * ( Iz * matmul(d3psi,v) + &
                         dIz * matmul(d2psi,v) ) ! shear force
    case(6)
      elemvec = matmul(phi,u) ! axial displacement
    case(7)
      elemvec = matmul(dphi,u) ! axial strain
    case(8)
      elemvec = Emod * Ac * matmul(dphi,u) ! axial force
    case default
      call errormsg_case_default ( 'beam_deriv2', 'comp', int_value=comp )
    end select

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( pos, pos1, pos2, posu, posv )
      deallocate ( e, n, Q )
      deallocate ( xrnod )
      deallocate ( phi, dphi_r, dphi )
      deallocate ( psi_r, dpsi_r )
      deallocate ( d2psi_r, d3psi_r )
      deallocate ( psi, dpsi )
      deallocate ( d2psi, d3psi )
      deallocate ( x )
      deallocate ( Ac, Iz, dIz, sol )
      deallocate ( u, v )

    end if

  end subroutine beam_deriv2


! Displacement for a beam element in all nodes. Two-dimensional.
! Needed for elements with center node.

  subroutine beam_displacement2 ( mesh, problem, elgrp, elem, first, last, &
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

      call set_globals_beam2 ( mesh, coefficients, elgrp )

      allocate ( pos(2), pos1(ndim), pos2(ndim), posu(ndfu), posv(ndfv) )
      allocate ( e(ndim), n(ndim), Q(ndim,ndim) )
      allocate ( xrnod(nodalp,1) )
      allocate ( phi(nodalp,ndfu) )
      allocate ( psi_r(nodalp,ndfv) )
      allocate ( psi(nodalp,ndfv) )
      allocate ( x(nodalp,ndim) )
      allocate ( sol(ndfu+ndfv) )
      allocate ( u(ndfu), v(ndfv) )
      allocate ( work(ndim,nodalp) )

!     set reference coordinates and shape function

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      select case ( intpolu )
      case(2)
        call shape_line_P1 ( xrnod(:,1), phi )
      case(6)
        call shape_line_P2 ( xrnod(:,1), phi )
      case default
        call errormsg_case_default ( 'beam_displacement2', &
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
        call errormsg_case_default ( 'beam_displacement2', &
          'intpolv', int_value=intpolv )
      end select

!     positions of degrees of freedom in element vector, matrix

      call set_pos_beam2

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the beam element

    e = x(nodalp,:) - x(1,:)
    h = sqrt(dot_product(e,e))
    e = e / h
    n = [ -e(2), e(1) ]
    Q(1,:) = e
    Q(2,:) = n

!   change slope degrees from reference to real nodal values

    psi = psi_r
    psi(:,pos) = h * psi_r(:,pos) / 2

!   get solution

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, sol )

!   transform end node displacements to local degrees

    sol(pos1) = matmul(Q,sol(pos1))
    sol(pos2) = matmul(Q,sol(pos2))

!   extract axial and transverse displacements

    u = sol(posu)
    v = sol(posv)

!   compute displacements in all nodes

    work(1,:) = matmul(phi,u)
    work(2,:) = matmul(psi,v)

    work = matmul( transpose(Q), work )

    elemvec = reshape( work, [ndim*nodalp] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( pos, pos1, pos2, posu, posv )
      deallocate ( e, n, Q )
      deallocate ( xrnod )
      deallocate ( phi )
      deallocate ( psi_r )
      deallocate ( psi )
      deallocate ( x )
      deallocate ( sol )
      deallocate ( u, v )
      deallocate ( work )

    end if

  end subroutine beam_displacement2


! Rotation for a beam element in all nodes. Two-dimensional.
! Needed for elements with center node.

  subroutine beam_rotation2 ( mesh, problem, elgrp, elem, first, last, &
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

      call set_globals_beam2 ( mesh, coefficients, elgrp )

      allocate ( pos(2), pos1(ndim), pos2(ndim), posu(ndfu), posv(ndfv) )
      allocate ( e(ndim), n(ndim), Q(ndim,ndim) )
      allocate ( xrnod(nodalp,1) )
      allocate ( psi_r(nodalp,ndfv), dpsi_r(nodalp,ndfv) )
      allocate ( psi(nodalp,ndfv), dpsi(nodalp,ndfv) )
      allocate ( x(nodalp,ndim) )
      allocate ( sol(ndfu+ndfv) )
      allocate ( v(ndfv) )

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
        call errormsg_case_default ( 'beam_rotation2', &
          'intpolv', int_value=intpolv )
      end select

!     positions of degrees of freedom in element vector, matrix

      call set_pos_beam2

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the beam element

    e = x(nodalp,:) - x(1,:)
    h = sqrt(dot_product(e,e))
    e = e / h
    n = [ -e(2), e(1) ]
    Q(1,:) = e
    Q(2,:) = n

!   change slope degrees from reference to real nodal values

    dpsi = dpsi_r
    dpsi(:,pos) = h * dpsi_r(:,pos) / 2

!   change reference derivative to real derivative

    dpsi  = dpsi * 2 / h

!   get solution

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, sol )

!   transform end node displacements to local degrees

    sol(pos1) = matmul(Q,sol(pos1))
    sol(pos2) = matmul(Q,sol(pos2))

!   extract transverse displacements

    v = sol(posv)

!   compute slope dvdx

    elemvec = matmul(dpsi,v) ! slope

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( pos, pos1, pos2, posu, posv )
      deallocate ( e, n, Q )
      deallocate ( xrnod )
      deallocate ( psi_r, dpsi_r )
      deallocate ( psi, dpsi )
      deallocate ( x )
      deallocate ( sol )
      deallocate ( v )

    end if

  end subroutine beam_rotation2


! Distributed load for a beam element in all nodes. Two-dimensional beam.
! For post processing purposes only.

  subroutine beam_q2 ( mesh, problem, elgrp, elem, first, last, &
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

      call set_globals_beam2 ( mesh, coefficients, elgrp )

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

  end subroutine beam_q2


! set pos arrays for a 2D beam element

  subroutine set_pos_beam2

    use structures_globals_m

!   positions of degrees of freedom in element vector, matrix

    pos1 = [ 1, 2 ]  ! displacements in first node

    select case ( intpolu )
    case(2)
      select case ( intpolv )
      case(0)
        pos2 = [ 4, 5 ]        ! displacements in last node
        posu = [ 1, 4 ]        ! axial displacement degrees of freedom
        posv = [ 2, 3, 5, 6 ]  ! transverse displacement degrees of freedom
      case(1)
        pos2 = [ 5, 6 ]
        posu = [ 1, 5 ]
        posv = [ 2, 3, 4, 6, 7 ]
      case default
        call errormsg_case_default ( 'set_pos_beam2', &
          'intpolv', int_value=intpolv )
      end select
    case(6)
      select case ( intpolv )
      case(0)
        pos2 = [ 5, 6 ]
        posu = [ 1, 4, 5 ]
        posv = [ 2, 3, 6, 7 ]
      case(1)
        pos2 = [ 6, 7 ]
        posu = [ 1, 4, 6 ]
        posv = [ 2, 3, 5, 7, 8 ]
      case default
        call errormsg_case_default ( 'set_pos_beam2', &
          'intpolv', int_value=intpolv )
      end select
    case default
      call errormsg_case_default ( 'set_pos_beam2', &
        'intpolu', int_value=intpolu )
    end select

  end subroutine set_pos_beam2


! set global parameters for a beam element (internal element)

  subroutine set_globals_beam2 ( mesh, coefficients, elgrp )

    use structures_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

!   check size of coefficients

    call check ( coefficients, 'set_globals_beam', ncoefi=100, &
      ncoefr=50, indexarray=[3,12], minimum=[0,0], maximum=[1,6] )

    ndim = mesh%element(elgrp)%ndim
    nodalp = mesh%element(elgrp)%numnod
    globalshape = mesh%element(elgrp)%globalshape

    if ( globalshape /= 'line' ) then
      write(*,'(/3(a/))') 'Error in set_globals_beam2:', &
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
      write(*,'(/a/a,i0/)') 'Error in set_globals_beam2:', &
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
      write(*,'(/a/a,i0/)') 'Error in set_globals_beam2:', &
        ' Invalid interpolation, intpolv = ', intpolv
      stop
    end select

    ninti = coefficients%i(2)

  end subroutine set_globals_beam2

end module structures_elements_beam2_m

