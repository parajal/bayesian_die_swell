
! Copyright (C) 2012-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for 1D beam structures

module structures_elements_beam1_m

  use tfem_elem_m
  use shapefunc_hermite_m

  implicit none

contains

! Internal element routine for a beam element with a linear elastic
! material law (Hooke's law). Geometrically linear. One-dimensional.

  subroutine beam_elem1 ( mesh, problem, elgrp, elem, matrix, &
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

      call set_globals_beam1 ( mesh, coefficients, elgrp )

      if ( ndim /= 1 ) then
        write(*,'(/3(a/))') 'Error in beam_elem1:', &
        ' Mesh must be one-dimensional.'
        stop
      end if

!     allocate arrays

      allocate ( pos(2) )
      allocate ( xig(ninti), wg(ninti) )
      allocate ( psi_r(ninti,ndfv), d2psi_r(ninti,ndfv) )
      allocate ( psi(ninti,ndfv), d2psi(ninti,ndfv) )
      allocate ( x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( qg(ninti), Iz(ninti) )

!     Gauss rule and shape function

      call Gauss_Legendre_numeric_line ( ninti, xig, wg )

      select case ( intpolv )
      case(0)
        call shape_line_hermite_P3 ( xig, psi_r, d2phi=d2psi_r )
        pos = [ 2, 4 ]
      case(1)
        call shape_line_hermite_P4 ( xig, psi_r, d2phi=d2psi_r )
        pos = [ 2, 5 ]
      case default
        call errormsg_case_default ( 'beam_elem1', 'intpolv', &
          int_value=intpolv )
      end select

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the beam element

    h = abs(x(nodalp,1) - x(1,1))
    xg(:,1) = x(1,1) * ( 1 - xig ) / 2 + x(nodalp,1) * ( 1 + xig ) / 2

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

    d2psi = d2psi * 4 / h ** 2

    if ( vector ) then

!     evaluate right-hand side

      if ( coefficients%i(4) >=0 ) then

        call evaluate_scalar_coefficient_simple ( choice=coefficients%i(4), &
          value=coefficients%r(4), func=coefficients%func, &
          funcnr=coefficients%i(5), x=xg, coef=qg )

        do i = 1, ndfv
          elemvec(i) = sum ( qg * psi(:,i) * wg ) * h / 2
        end do

      else

        elemvec = 0

      end if

    end if

    if ( matrix ) then

!     stiffness matrix

      do i = 1, ndfv
        do j = i, ndfv
          elemmat(i,j) = sum ( Iz * d2psi(:,i) * d2psi(:,j) * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

      elemmat = Emod * elemmat * h / 2

    end if

    if ( last ) then

!     last element in this group

      deallocate ( pos )
      deallocate ( xig, wg )
      deallocate ( psi_r, d2psi_r )
      deallocate ( psi, d2psi )
      deallocate ( x, xg )
      deallocate ( qg, Iz )

    end if

  end subroutine beam_elem1


! Compute derivatives for a beam element with a linear elastic
! material law (Hooke's law). Geometrically linear. One-dimensional.

  subroutine beam_deriv1 ( mesh, problem, elgrp, elem, first, last, &
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

      call set_globals_beam1 ( mesh, coefficients, elgrp )

!     check component = coefficients%i(1)

      call check ( coefficients, 'beam_deriv1', indexarray=[1], &
        minimum=[1], maximum=[5] )

      allocate ( pos(2), xrnod(nodalp,1) )
      allocate ( psi_r(nodalp,ndfv), dpsi_r(nodalp,ndfv) )
      allocate ( d2psi_r(nodalp,ndfv), d3psi_r(nodalp,ndfv) )
      allocate ( psi(nodalp,ndfv), dpsi(nodalp,ndfv) )
      allocate ( d2psi(nodalp,ndfv), d3psi(nodalp,ndfv) )
      allocate ( x(nodalp,ndim) )
      allocate ( Iz(nodalp), dIz(nodalp), u(ndfv) )

!     set reference coordinates and shape function

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

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
        call errormsg_case_default ( 'beam_deriv1', 'intpolv', &
          int_value=intpolv )
      end select

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the beam element

    h = abs(x(nodalp,1) - x(1,1))

!   second moment of area

    call evaluate_scalar_coefficient_simple ( choice=coefficients%i(6), &
      value=coefficients%r(3), func=coefficients%func1(1)%p, &
      funcnr=coefficients%i(7), x=x, coef=Iz )

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

    dpsi  = dpsi * 2 / h
    d2psi = d2psi * 4 / h ** 2
    d3psi = d3psi * 8 / h ** 3

!   get solution

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u )

!   compute derivative

    comp = coefficients%i(1)

    select case(comp)
    case(1)
      elemvec = matmul(psi,u) ! vertical displacement
    case(2)
      elemvec = matmul(dpsi,u) ! slope
    case(3)
      elemvec = matmul(d2psi,u) ! second derivative
    case(4)
      elemvec = Emod * Iz * matmul(d2psi,u) ! bending moment
    case(5)
      dIz = (Iz(nodalp)-Iz(1))/(x(nodalp,1)-x(1,1))
      elemvec = Emod * ( Iz * matmul(d3psi,u) + &
                         dIz * matmul(d2psi,u) ) ! shear force
    case default
      call errormsg_case_default ( 'beam_deriv1', 'comp', int_value=comp )
    end select

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( pos, xrnod )
      deallocate ( psi_r, dpsi_r )
      deallocate ( d2psi_r, d3psi_r )
      deallocate ( psi, dpsi )
      deallocate ( d2psi, d3psi )
      deallocate ( x )
      deallocate ( Iz, dIz, u )

    end if

  end subroutine beam_deriv1


! Displacement for a beam element in all nodes. One-dimensional beam.
! Needed for elements with center node.

  subroutine beam_displacement1 ( mesh, problem, elgrp, elem, first, last, &
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

      call set_globals_beam1 ( mesh, coefficients, elgrp )

      allocate ( pos(2), xrnod(nodalp,1) )
      allocate ( psi_r(nodalp,ndfv) )
      allocate ( psi(nodalp,ndfv) )
      allocate ( x(nodalp,ndim) )
      allocate ( sol(ndfv) )
      allocate ( work(2,nodalp) )

!     set reference coordinates and shape function

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      select case ( intpolv )
      case(0)
        call shape_line_hermite_P3 ( xrnod(:,1), psi_r )
        pos = [ 2, 4 ]
      case(1)
        call shape_line_hermite_P4 ( xrnod(:,1), psi_r )
        pos = [ 2, 5 ]
      case default
        call errormsg_case_default ( 'beam_displacement1', &
          'intpolv', int_value=intpolv )
      end select

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the beam element

    h = abs(x(nodalp,1) - x(1,1))

!   change slope degrees from reference to real nodal values

    psi = psi_r
    psi(:,pos) = h * psi_r(:,pos) / 2

!   get solution

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, sol )

!   compute displacements in all nodes

    work(1,:) = 0
    work(2,:) = matmul(psi,sol)

    elemvec = reshape( work, [2*nodalp] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( pos, xrnod )
      deallocate ( psi_r )
      deallocate ( psi )
      deallocate ( x )
      deallocate ( sol )
      deallocate ( work )

    end if

  end subroutine beam_displacement1


! Distributed load for a beam element in all nodes. One-dimensional beam.
! For post processing purposes only.

  subroutine beam_q1 ( mesh, problem, elgrp, elem, first, last, &
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

      call set_globals_beam1 ( mesh, coefficients, elgrp )

      allocate ( x(nodalp,ndim) )
      allocate ( work(2,nodalp) )
      allocate ( qnod(nodalp) )

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   compute distributed loading

    if ( coefficients%i(4) >=0 ) then

      call evaluate_scalar_coefficient_simple ( choice=coefficients%i(4), &
        value=coefficients%r(4), func=coefficients%func, &
        funcnr=coefficients%i(5), x=x, coef=qnod )

    else

      qnod = 0

    end if

!   set q vector in all nodes

    work(1,:) = 0
    work(2,:) = qnod

    elemvec = reshape( work, [2*nodalp] )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( x )
      deallocate ( work )
      deallocate ( qnod )

    end if

  end subroutine beam_q1


! set global parameters for a beam element (internal element)

  subroutine set_globals_beam1 ( mesh, coefficients, elgrp )

    use structures_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

!   check size of coefficients

    call check ( coefficients, 'set_globals_beam', ncoefi=100, &
      ncoefr=50, indexarray=[3], minimum=[0], maximum=[1] )

    ndim = mesh%element(elgrp)%ndim
    nodalp = mesh%element(elgrp)%numnod
    globalshape = mesh%element(elgrp)%globalshape

    if ( globalshape /= 'line' ) then
      write(*,'(/3(a/))') 'Error in set_globals_beam1:', &
      ' Only line elements are allowed'
      stop
    end if

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
      write(*,'(/a/a,i0/)') 'Error in set_globals_beam1:', &
        ' Invalid interpolation, intpolv = ', intpolv
      stop
    end select

    ninti = coefficients%i(2)

  end subroutine set_globals_beam1

end module structures_elements_beam1_m

