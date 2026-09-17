
! Copyright (C) 2012-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for truss structures.

module structures_elements_truss_m

  use tfem_elem_m

  implicit none

contains


! Internal element routine for a truss element with a linear elastic
! material law (Hooke's law). Geometrically linear.

  subroutine truss_elem1 ( mesh, problem, elgrp, elem, matrix, &
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


    integer :: i
    real(dp) :: l, Emod, A


    if ( first ) then

!     first element in this group

      call set_globals_truss ( mesh, coefficients, elgrp )

      allocate ( pos1(ndim), pos2(ndim) )
      allocate ( x(nodalp,ndim), e(ndim), emat(ndim,ndim) )
      allocate ( fe(ndim) )

      pos1 = [(i,i=1,ndim)]
      pos2 = pos1 + ndim

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the truss element

    e = x(nodalp,:) - x(1,:)
    l = sqrt(dot_product(e,e))
    e = e / l

    A = coefficients%r(1)

!   material properties

    Emod = coefficients%r(2)

    if ( vector ) then

!     evaluate right-hand side

      if ( coefficients%i(10) == 0 ) then

!       distributed load (gravity): split total force 50-50 to the end nodes

        fe = coefficients%r(5:4+ndim) * l / 2

        elemvec(pos1) = fe
        elemvec(pos2) = fe

      else

       elemvec = 0

      end if

    end if

    if ( matrix ) then

!     stiffness matrix

      do i = 1, ndim
        emat(i,:) = e(i) * e
      end do

      elemmat(pos1,pos1) = emat
      elemmat(pos1,pos2) = - emat
      elemmat(pos2,pos1) = - emat
      elemmat(pos2,pos2) = emat

      elemmat = Emod * A / l * elemmat

    end if

    if ( last ) then

!     last element in this group

      deallocate ( pos1, pos2 )
      deallocate ( x, e, emat )
      deallocate ( fe )

    end if

  end subroutine truss_elem1


! compute derivatives for truss

  subroutine truss_deriv1 ( mesh, problem, elgrp, elem, first, last, &
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
    real(dp) :: l, Emod, A, epsilon


    if ( first ) then

!     first element in this group

      call set_globals_truss ( mesh, coefficients, elgrp )

!     check component = coefficients%i(1)

      call check ( coefficients, 'truss_deriv1', indexarray=[1], &
        minimum=[1], maximum=[3] )

      allocate ( x(nodalp,ndim), e(ndim) )
      allocate ( uvec(ndim,ndf), u(ndim*ndf) )

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   geometry of the truss element

    e = x(nodalp,:) - x(1,:)
    l = sqrt(dot_product(e,e))
    e = e / l

    A = coefficients%r(1)

!   material properties

    Emod = coefficients%r(2)

!   get displacements

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u )

    uvec = reshape ( u, [ndim,ndf] )

!   strain

    epsilon = dot_product ( e, uvec(:,2) - uvec(:,1) ) / l

!   compute derivative

    comp = coefficients%i(1)

    select case(comp)
    case(1)
      elemvec = epsilon  ! strain
    case(2)
      elemvec = Emod * epsilon  ! stress
    case(3)
      elemvec = Emod * A * epsilon  ! normal force
    case default
      call errormsg_case_default ( 'truss_deriv1', 'comp', int_value=comp )
    end select

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( x, e )
      deallocate ( uvec, u )

    end if

  end subroutine truss_deriv1


! Distributed load for a truss element in all nodes.
! For post processing purposes only.

  subroutine truss_q1 ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use structures_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    integer :: i

    if ( first ) then

!     first element in this group

      call set_globals_truss ( mesh, coefficients, elgrp )

      allocate ( x(nodalp,ndim) )
      allocate ( qvecnod(nodalp,ndim) )

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

!   compute distributed loading

    if ( coefficients%i(10) == 0 ) then

!     distributed load

      do i = 1, nodalp
        qvecnod(i,:) = coefficients%r(5:4+ndim)
      end do

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

  end subroutine truss_q1


! set global parameters for a truss element (internal element)

  subroutine set_globals_truss ( mesh, coefficients, elgrp )

    use structures_globals_m

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

    ndim = mesh%element(elgrp)%ndim
    nodalp = mesh%element(elgrp)%numnod
    globalshape = mesh%element(elgrp)%globalshape

    if ( globalshape /= 'line' ) then
      write(*,'(/3(a/))') 'Error in set_globals_truss:', &
      ' Only line elements are allowed'
      stop
    end if

    ndf = 2

  end subroutine set_globals_truss

end module structures_elements_truss_m

