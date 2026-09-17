
! Copyright (C) 2010-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! subdivide basic element shapes (line, triangle, quadrilateral, tetrahedron and
! hexahedron) into lines, triangles or tetrahedrons such that all elements are
! on either side of the levelset=0 point, line or surface.
! This module creates the levelset=0 geometry.
! The levelset function is defined by the values in the vertex points and is
! linearly interpolated along the sides of the element.

module intmesh_m

  use glob_defs_m
  use set_optional_m

  implicit none


! simple interface mesh type

  type imesh_t

    integer :: ndim     = 0        ! space dimension
    integer :: nnodes   = 0        ! number of nodal points
    integer :: nelem    = 0        ! number of elements

!   shape of the element
!      elshape=0 : 1 node "point" element
!      elshape=1 : 2 node line element
!      elshape=3 : 3 node triangle
    integer :: elshape  = 0
    integer :: elnumnod   = 0     ! number of nodal points of a single element

!   connections of elements to nodal points
!   topology(node,elem) gives the global node number of the local node
!   in element elem
    integer, allocatable, dimension(:,:) :: topology

!   coordinates coor(nnodes,ndim)
    real(dp), allocatable, dimension(:,:) :: coor

!   jacobian(elem): jacobian of each element
!   Note: for elshape=0 the jacobian will be zero.
    real(dp), allocatable, dimension(:) :: jacobian

!   jacobiannormal(elem,ndim): jacobian * normal vector of each element.
!   The vector is pointing in the direction of increasing levelset function.
!   (from - to + ).
!   Note: for elshape=0 jacobiannormal will be zero.
    real(dp), allocatable, dimension(:,:) :: jacobiannormal

  end type imesh_t


! interface for generic create subroutine

  interface create
    module procedure create_imesh
  end interface create

! interface for generic delete subroutine

  interface delete
    module procedure delete_imesh
  end interface delete

! interface for generic volume subroutine

  interface area
    module procedure area_imesh
  end interface area

! interface for generic number_of_integration_points subroutine

  interface number_of_integration_points
    module procedure number_of_integration_points_imesh
  end interface number_of_integration_points

! interface for generic integration_points subroutine

  interface integration_points
    module procedure integration_points_imesh
  end interface integration_points


contains


! create imesh

  subroutine create_imesh ( imesh, ndim, nnodes, nelem, elshape )

    type(imesh_t), intent(inout) :: imesh
    integer, intent(in) :: ndim, nnodes, nelem, elshape

    imesh%ndim = ndim
    imesh%nnodes = nnodes
    imesh%nelem = nelem
    imesh%elshape = elshape

    select case(elshape)
    case(0)
!     point element
      imesh%elnumnod = 1
    case(1)
!     line element
      imesh%elnumnod = 2
    case(3)
!     triangle
      imesh%elnumnod = 3
    case default
      write(*,'(/a,i0/)') &
        'Error create_imesh: invalid elshape = ', elshape
      stop
    end select

    allocate(imesh%topology(imesh%elnumnod,nelem))
    allocate(imesh%coor(nnodes,ndim))

  end subroutine create_imesh


! delete imesh

  subroutine delete_imesh ( imesh )

    type(imesh_t), intent(out) :: imesh

  end subroutine delete_imesh


! subdivide a line element into line elements according to the levelset
! function s

  subroutine intmesh_line_element ( x, s, imesh )

!   coordinates of the end points
    real(dp), dimension(:,:), intent(in) :: x

!   levelset function of the end points
    real(dp), dimension(:), intent(in) :: s

!   the "interface" mesh
    type(imesh_t), intent(inout), optional :: imesh


!   counter of zero or positive level set
    integer :: countz, countp
!   the indices for the generic nodes 1, 2
    integer :: i1, i2
!   the "other" node
    integer, parameter :: io(2)=[2,1]
!   interpolation xi
    real(dp) :: xi
!   all nodes negative of positive
    logical :: alln, allp


    alln = all ( s < 0._dp )
    allp = all ( s > 0._dp )
    countz = count ( s == 0._dp )

    if ( alln .or. allp .or. countz == 2 ) then

!     one line element only

!     interface is absent for this case

      call create ( imesh, ndim=size(x,2), nnodes=0, nelem=0, elshape=0 )

    else if ( countz == 1 ) then

!     one line element only; degenerate case
!     We have the following generic case:
!
!        1-----------2
!
!     The nodes 1 and 2 correspond with nodes of the original line element.
!     The node 1 is on the zero levelset and node 2 has a positive or
!     negative levelset.
!
      countp = count ( s > 0._dp )

      if ( countp == 1 ) then

!       one positive vertex: interface is in this element

        i1 = maxval ( minloc (s) )
        i2 = io(i1)
        if ( s(i2) <= 0._dp ) stop 'internal error intmesh_line_element'

!       a single interface point (point 1)

        call create ( imesh, ndim=size(x,2), nnodes=1, nelem=1, elshape=0 )

        imesh%coor(1,:) = x(i1,:)

!       topology

        imesh%topology(1,1) = 1

      else if ( countp == 0 ) then

!       one negative vertex: interface is absent for this case

        call create ( imesh, ndim=size(x,2), nnodes=0, nelem=0, elshape=0 )

      end if

    else

!     Two line elements.
!
!     We have the following generic case:
!
!        1-------3-------2
!
!     The nodes 1 and 2 correspond with nodes of the original line element
!     The node 3 is on the zero levelset which is found by linear interpolation
!     between the levelset values in the vertex nodes.
!

!     linear interpolate to find the point on the zero levelset

      xi = s(1) / ( s(1) - s(2) )

!     a single interface point (point 3)

      call create ( imesh, ndim=size(x,2), nnodes=1, nelem=1, elshape=0 )

      imesh%coor(1,:) = (1-xi) * x(1,:) + xi * x(2,:)

!     topology

      imesh%topology(1,1) = 1

    end if

  end subroutine intmesh_line_element


! subdivide a triangle into subtriangles according to the levelset function s

  subroutine intmesh_triangle ( x, s, imesh )

!   coordinates of the vertices
    real(dp), dimension(:,:), intent(in) :: x

!   levelset function of the vertices
    real(dp), dimension(:), intent(in) :: s

!   the "interface" mesh
    type(imesh_t), intent(inout), optional :: imesh


!   the indices for the generic nodes 1, 2, 3
    integer :: i1, i2, i3
!   the "other" two nodes
    integer, parameter :: io(2,3)=reshape([2,3,3,1,1,2],[2,3])
!   number of negative values of the levelset
    integer :: countn, countz
!   interpolation xi between the sides 1-2 and 1-3
    real(dp) :: xi12, xi13, xi23
!   all vertices negative of positive
    logical :: alln, allp, cutvertex


    alln = all ( s < 0._dp )
    allp = all ( s > 0._dp )
    countz = count ( s == 0._dp )
    countn = count ( s < 0._dp )
    cutvertex = countz == 1 .and. ( countn == 0 .or. countn == 2 )

    if ( alln .or. allp .or. countz == 3 .or. cutvertex ) then

!     one triangle only

!     interface is absent for this case

      call create ( imesh, ndim=size(x,2), nnodes=0, nelem=0, elshape=1 )

    else if ( countz == 2 ) then

!     one triangle element only; degenerate case
!
!                3
!               / \             !
!              /   \            !
!             /     \           !
!            /       \          !
!           /         \         !
!          /           \        !
!         /             \       !
!        1---------------2
!
!     The nodes 1, 2 and 3 correspond with vertices of the original triangle
!     (possibly permutated).
!     Here, 1 is the node that is positive or negative and 2 and 3 are
!     on the zero levelset.
!
      if ( countn == 0 ) then

!       one positive vertex: interface is in this element

        i1 = maxval ( maxloc (s) )
        i2 = io(1,i1)
        i3 = io(2,i1)

!       a single interface line

        call create ( imesh, ndim=size(x,2), nnodes=2, nelem=1, elshape=1 )

        imesh%coor(1,:) = x(i2,:)
        imesh%coor(2,:) = x(i3,:)

!       topology

        imesh%topology(:,1) = [ 2, 1 ]

      else if ( countn == 1 ) then

!       one negative vertex: interface is absent for this case

        call create ( imesh, ndim=size(x,2), nnodes=0, nelem=0, elshape=1 )

      else
        stop 'internal error intmesh_triangle 1'
      end if

    else if ( countz == 1 .and. countn == 1 ) then

!     two triangles; degenerate case
!
!     We have the following generic case:
!
!                3
!                |\              !
!                | \             !
!                |  \            !
!                |   4           !
!                |  / \          !
!                | /   \         !
!                |/     \        !
!                1-------2
!
!     The nodes 1, 2 and 3 correspond with vertices of the original triangle
!     (possibly permutated).
!     Here, 2 is the node that is positive (negative) and 3 is
!     negative (positive). The line between the nodes 1-4 is the zero levelset.
!     The position of node 4 is is found by linear interpolation between the
!     levelset values in the vertex nodes 2 and 3.
!

      i1 = maxval ( minloc(abs(s)) )

!     find the other two nodes

      i2 = io(1,i1)
      i3 = io(2,i1)

!     linear interpolate to find the points on the zero levelset

      xi23 = s(i2) / ( s(i2) - s(i3) )

      call create ( imesh, ndim=size(x,2), nnodes=2, nelem=1, elshape=1 )

      imesh%coor(1,:) = x(i1,:)
      imesh%coor(2,:) = (1-xi23) * x(i2,:) + xi23 * x(i3,:)

!     direction according to the levelset

      if ( s(i2) > 0._dp ) then

        imesh%topology(:,1) = [ 1, 2 ]

      else if ( s(i2) < 0._dp ) then

        imesh%topology(:,1) = [ 2, 1 ]

      else
        stop 'internal error intmesh_triangle 2'
      end if

    else

!     Three triangles.
!
!     We have the following generic case:
!
!                3
!               /|\             !
!              / | \            !
!             /  |  \           !
!            5   |   \          !
!           / \  |    \         !
!          /   \ |     \        !
!         /     \|      \       !
!        1-------4-------2
!
!     The nodes 1, 2 and 3 correspond with vertices of the original triangle
!     (possibly permutated).
!     Here, 1 is the node that is positive (negative) and 2 and 3 are
!     negative (positive). The line between the nodes 4-5 is the zero levelset
!     which is found by linear interpolation between the levelset values in the
!     vertex nodes.
!

      countn = count ( s < 0._dp )

      if ( countn == 1 ) then

!       one negative vertex

        i1 = maxval ( minloc (s) )

      else if ( countn == 2 ) then

!       one positive vertex

        i1 = maxval ( maxloc (s) )

      else
        stop 'internal error intmesh_triangle 3'
      end if

!     find the other two nodes

      i2 = io(1,i1)
      i3 = io(2,i1)

!     linear interpolate to find the points on the zero levelset

      xi12 = s(i1) / ( s(i1) - s(i2) )
      xi13 = s(i1) / ( s(i1) - s(i3) )

      call create ( imesh, ndim=size(x,2), nnodes=2, nelem=1, elshape=1 )

      imesh%coor(1,:) = (1-xi12) * x(i1,:) + xi12 * x(i2,:)
      imesh%coor(2,:) = (1-xi13) * x(i1,:) + xi13 * x(i3,:)

      if ( countn == 1 ) then

!       one negative vertex

        imesh%topology(:,1) = [ 1, 2 ]

      else if ( countn == 2 ) then

!       one positive vertex

        imesh%topology(:,1) = [ 2, 1 ]

      end if

    end if

  end subroutine intmesh_triangle


! subdivide a quad into subtriangles according to the levelset function s

  subroutine intmesh_quadrilateral ( x, s, imesh )

!   coordinates of the vertices
    real(dp), dimension(:,:), intent(in) :: x

!   levelset function of the vertices
    real(dp), dimension(:), intent(in) :: s

!   the "interface" mesh
    type(imesh_t), intent(inout), optional :: imesh


    type(imesh_t) :: imesh1, imesh2

!   subdivide into two triangles
!
!       4-----------3
!       |        __/|
!       |     __/   |
!       |  __/      |
!       | /         |
!       1-----------2
!

      call intmesh_triangle ( x([1,2,3],:), s([1,2,3]), imesh1 )
      call intmesh_triangle ( x([1,3,4],:), s([1,3,4]), imesh2 )

      call merge_imesh ( imesh1, imesh2, imesh )

      call delete ( imesh1 )
      call delete ( imesh2 )

  end subroutine intmesh_quadrilateral


! subdivide a tetrahedron into subtets according to the levelset function s

  subroutine intmesh_tetrahedron ( x, s, imesh )

!   coordinates of the vertices
    real(dp), dimension(:,:), intent(in) :: x

!   levelset function of the vertices
    real(dp), dimension(:), intent(in) :: s

!   the "interface" mesh
    type(imesh_t), intent(inout) :: imesh


    integer :: num, i, indx(2), ic
!   the indices for the generic nodes 1, 2, 3, 4
    integer :: i1, i2, i3, i4
!   the "other" three nodes (generic case 1))
    integer, parameter :: io1(3,4)=reshape([2,3,4,3,1,4,1,2,4,1,3,2],[3,4])
!   the current and "other" two nodes
    integer, parameter :: ic2(2,6)=reshape([1,2,2,3,3,1,4,2,1,4,3,4],[2,6])
    integer, parameter :: io2(2,6)=reshape([3,4,1,4,2,4,1,3,2,3,1,2],[2,6])
!   the current and the other nodes (degenerate case with 2 and 3 tetrahedrons).
    integer, parameter :: ic3(2,6)=reshape([1,2,2,3,3,1,4,2,1,4,3,4],[2,6])
    integer, parameter :: io3(2,6)=reshape([3,4,1,4,2,4,1,3,2,3,1,2],[2,6])
!   the "other" three nodes (degenerate case with three zeros))
    integer, parameter :: io4(3,4)=reshape([3,2,4,1,3,4,2,1,4,1,2,3],[3,4])
!   number of negative, zero values of the levelset
    integer :: countn, countz
!   interpolation xi between the sides 1-2, 1-4, 2-4, 2-3, 1-3 and 3-4
    real(dp) :: xi12, xi14, xi24, xi23, xi13, xi34
!   all vertices negative of positive
    logical :: alln, allp
!   cut-off vertex with a zero levelset
    logical :: cutvertex
!   cut-off edge with a zero levelset
    logical :: cutedge
!   debugging
    logical, parameter :: debug = .false.


    alln = all ( s < 0._dp )
    allp = all ( s > 0._dp )
    countz = count ( s == 0._dp )
    countn = count ( s < 0._dp )
    cutvertex = countz == 1 .and. ( countn == 0 .or. countn == 3 )
    cutedge = countz == 2 .and. ( countn == 0 .or. countn == 2 )

    if ( alln .or. allp .or. cutvertex .or. cutedge .or. countz == 4  ) then

!     one tetrahedron only

!     interface is absent for this case

      call create ( imesh, ndim=size(x,2), nnodes=0, nelem=0, elshape=3 )

    else if ( countz == 1 ) then

      if (debug) print *, '3 tets'

!     three tetrahedrons; degenerate case

!     We have the generic case where one vertex is zero and not all
!     others have the same sign.
!
!     The nodes 1, 2, 3 and 4 correspond with vertices of the original
!     tetrahedron (possibly permutated).
!     Here, 1 is the node with a zero levelset. Node 2 is positive (negative)
!     and 3 and 4 are negative (positive). The triangle between the nodes
!     1-5-6 is the zero levelset. The nodes 5-6 are found by linear
!     interpolation between the levelset values in the vertex nodes 2-3 and 2-4,
!     respectively

      call create ( imesh, ndim=size(x,2), nnodes=3, nelem=1, elshape=3 )

      i1 = maxval ( minloc(abs(s)) )

      if ( countn == 1 ) then

!       one negative vertex

        i2 = maxval ( minloc (s) )

      else if ( countn == 2 ) then

!       one positive vertex

        i2 = maxval ( maxloc (s) )

      end if

      indx = [ i1, i2 ]

!     test which of the twelve possible combinations we have

      do i = 1, 6
        if ( all ( indx == ic3(:,i) ) ) then
          ic = i
          exit
        end if
        if ( all ( indx == ic3(2:1:-1,i) ) ) then
          ic = -i
          exit
        end if
      end do

!     find the other two nodes

      if ( ic > 0 ) then
        i3 = io3(1,ic)
        i4 = io3(2,ic)
      else
        i3 = io3(2,-ic)
        i4 = io3(1,-ic)
      end if

!     linear interpolate to find the points on the zero levelset

      xi23 = s(i2) / ( s(i2) - s(i3) )
      xi24 = s(i2) / ( s(i2) - s(i4) )

      imesh%coor(1,:) = x(i1,:)
      imesh%coor(2,:) = (1-xi23)*x(i2,:) + xi23*x(i3,:)
      imesh%coor(3,:) = (1-xi24)*x(i2,:) + xi24*x(i4,:)

!     topology

      if ( countn == 1 ) then

!       s(i2) < 0

        imesh%topology(:,1) = [ 1, 3, 2  ]

      else if ( countn == 2 ) then

!       s(i2) > 0

        imesh%topology(:,1) = [ 1, 2, 3  ]

      else
        stop 'internal error intmesh_tetrahedron 1'
      end if

    else if ( countz == 2 ) then

      if (debug) print *, '2 tets'

!     two tetrahedron; degenerate case

!     We have the generic case where two vertices are zero and not all
!     others have the same sign.
!
!     The nodes 1, 2, 3 and 4 correspond with vertices of the original
!     tetrahedron (possibly permutated).
!     Here, 1 and 2 are the nodes having a zero levelset. Node 3 is negative
!     and 4 is negative. The triangle between the nodes 1-2-5 is the
!     zero levelset. The node 5 is found by linear interpolation
!     between the levelset values in the vertex nodes 3-4.

      call create ( imesh, ndim=size(x,2), nnodes=3, nelem=1, elshape=3 )

      i3 = maxval ( minloc(s) )
      i4 = maxval ( maxloc(s) )

      indx = [ i3, i4 ]

!     test which of the twelve possible combinations we have

      do i = 1, 6
        if ( all ( indx == ic3(:,i) ) ) then
          ic = i
          exit
        end if
        if ( all ( indx == ic3(2:1:-1,i) ) ) then
          ic = -i
          exit
        end if
      end do

!     find the other two nodes

      if ( ic > 0 ) then
        i1 = io3(1,ic)
        i2 = io3(2,ic)
      else
        i1 = io3(2,-ic)
        i2 = io3(1,-ic)
      end if

!     linear interpolate to find the point on the zero levelset

      xi34 = s(i3) / ( s(i3) - s(i4) )

      imesh%coor(1,:) = x(i1,:)
      imesh%coor(2,:) = x(i2,:)
      imesh%coor(3,:) = (1-xi34)*x(i3,:) + xi34*x(i4,:)

!     topology

      imesh%topology(:,1) = [ 1, 2, 3  ]

    else if ( countz == 3 ) then

      if (debug) print *, '1 tets'

!     one tetrahedron; degenerate case

!     We have the generic case where three vertices are zero and the
!     fourth is either positive or negative.
!
!     The nodes 1, 2, 3 and 4 correspond with vertices of the original
!     tetrahedron (possibly permutated).
!     Here, 1, 2 and 3 are the nodes having a zero levelset.
!     Node 4 is positive (negative).
!     The triangle between the nodes 1-2-3 is the zero levelset.

      if ( countn == 0 ) then

!       one positive vertex: interface is in this element

        i4 = maxval ( maxloc (s) )
        i1 = io4(1,i4)
        i2 = io4(2,i4)
        i3 = io4(3,i4)

!       a single interface triangle

        call create ( imesh, ndim=size(x,2), nnodes=3, nelem=1, elshape=3 )

        imesh%coor(1,:) = x(i1,:)
        imesh%coor(2,:) = x(i2,:)
        imesh%coor(3,:) = x(i3,:)

!       topology

        imesh%topology(:,1) = [ 1, 2, 3 ]

      else if ( countn == 1 ) then

!       one negative vertex: interface is absent for this case

        call create ( imesh, ndim=size(x,2), nnodes=3, nelem=0, elshape=3 )

      else
        stop 'internal error intmesh_tetrahedron 2'
      end if

    else if ( countn == 1 .or. countn == 3 ) then

      if (debug) print *, '4 tets'

!     Four tetrahedrons.
!
!     We have the generic case where one vertex is cut off with tetrahedron
!
!     The nodes 1, 2, 3 and 4 correspond with vertices of the original
!     tetrahedron (possibly permutated).
!     Here, 1 is the node that is positive (negative) and 2, 3 and 4 are
!     negative (positive). The triangle between the nodes 5-6-7 is the
!     zero levelset. The nodes 5-6-7 are found by linear interpolation
!     between the levelset values in the vertex nodes 1-2, 1-3 and 1-4,
!     respectively

      call create ( imesh, ndim=size(x,2), nnodes=3, nelem=1, elshape=3 )

      if ( countn == 1 ) then

!       one negative vertex

        i1 = maxval ( minloc (s) )

      else if ( countn == 3 ) then

!       one positive vertex

        i1 = maxval ( maxloc (s) )

      end if

!     find the other three nodes

      i2 = io1(1,i1)
      i3 = io1(2,i1)
      i4 = io1(3,i1)

!     linear interpolate to find the points on the zero levelset

      xi12 = s(i1) / ( s(i1) - s(i2) )
      xi13 = s(i1) / ( s(i1) - s(i3) )
      xi14 = s(i1) / ( s(i1) - s(i4) )

      imesh%coor(1,:) = (1-xi12)*x(i1,:) + xi12*x(i2,:)
      imesh%coor(2,:) = (1-xi13)*x(i1,:) + xi13*x(i3,:)
      imesh%coor(3,:) = (1-xi14)*x(i1,:) + xi14*x(i4,:)

!     inside or outside?

      if ( countn == 1 ) then

!       one negative vertex

        imesh%topology(:,1) = [ 1, 2, 3 ]

      else if ( countn == 3 ) then

!       one positive vertex

        imesh%topology(:,1) = [ 2, 1, 3 ]

      end if

    else if ( countn == 2 ) then

      if (debug) print *, '6 tets'

!     Six tetrahedrons.
!
!     We have the generic case where two vertices are positive and two are
!     negative.
!
!     The nodes 1, 2, 3 and 4 correspond with vertices of the original
!     tetrahedron (possibly permutated).
!     Here, 1 and 2 are the negative nodes and 3 and 4 are the
!     positive ones. The two triangles between the nodes 5-6-8 and
!     and 6-7-8 are the zero levelset. The nodes 5-6-7-8 are found by
!     linear interpolation between the levelset values in the vertex nodes
!     1-4, 2-4, 2-3 and 1-3, respectively.

      call create ( imesh, ndim=size(x,2), nnodes=4, nelem=2, elshape=3 )

!     find negative nodes

      num = 0
      do i = 1, 4
        if ( s(i) < 0 ) then
          num = num + 1
          indx(num) = i
        end if
      end do

!     test which of the six possible combinations we have

      do i = 1, 6
        if ( all ( indx == ic2(:,i) ) .or. &
             all ( indx == ic2(2:1:-1,i) ) ) then
          ic = i
          exit
        end if
      end do

!     the two nodes

      i1 = ic2(1,ic)
      i2 = ic2(2,ic)

!     find the other two nodes

      i3 = io2(1,ic)
      i4 = io2(2,ic)

!     linear interpolate to find the points on the zero levelset

      xi14 = s(i1) / ( s(i1) - s(i4) )
      xi24 = s(i2) / ( s(i2) - s(i4) )
      xi23 = s(i2) / ( s(i2) - s(i3) )
      xi13 = s(i1) / ( s(i1) - s(i3) )

      imesh%coor(1,:) = (1-xi14)*x(i1,:) + xi14*x(i4,:)
      imesh%coor(2,:) = (1-xi24)*x(i2,:) + xi24*x(i4,:)
      imesh%coor(3,:) = (1-xi23)*x(i2,:) + xi23*x(i3,:)
      imesh%coor(4,:) = (1-xi13)*x(i1,:) + xi13*x(i3,:)

!     topology

      imesh%topology(:,1) = [ 1, 2, 4  ]
      imesh%topology(:,2) = [ 2, 3, 4  ]

    else
        stop 'internal error intmesh_tetrahedron 3'
    end if

  end subroutine intmesh_tetrahedron


! subdivide a hexahedron into subtets according to the levelset function s

  subroutine intmesh_hexahedron ( x, s, imesh )

    use limits_m

!   coordinates of the vertices
    real(dp), dimension(:,:), intent(in) :: x

!   levelset function of the vertices
    real(dp), dimension(:), intent(in) :: s

    type(imesh_t), intent(inout) :: imesh


    type(imesh_t) :: imesh1, imesh2, imesh3


    if ( SUBDIVIDE_HEX6 ) then

!     subdivide into six tetrahedrons

      call intmesh_tetrahedron ( x([1,2,3,5],:), s([1,2,3,5]), imesh1 )
      call intmesh_tetrahedron ( x([2,3,5,7],:), s([2,3,5,7]), imesh2 )
      call merge_imesh ( imesh1, imesh2, imesh3 )
      call delete ( imesh1 )
      call delete ( imesh2 )
      call intmesh_tetrahedron ( x([2,7,5,6],:), s([2,7,5,6]), imesh1 )
      call merge_imesh ( imesh3, imesh1, imesh2 )
      call delete ( imesh3 )
      call delete ( imesh1 )
      call intmesh_tetrahedron ( x([1,5,3,8],:), s([1,5,3,8]), imesh1 )
      call merge_imesh ( imesh2, imesh1, imesh3 )
      call delete ( imesh2 )
      call delete ( imesh1 )
      call intmesh_tetrahedron ( x([3,5,7,8],:), s([3,5,7,8]), imesh1 )
      call merge_imesh ( imesh3, imesh1, imesh2 )
      call delete ( imesh3 )
      call delete ( imesh1 )
      call intmesh_tetrahedron ( x([1,3,4,8],:), s([1,3,4,8]), imesh1 )
      call merge_imesh ( imesh2, imesh1, imesh )
      call delete ( imesh2 )
      call delete ( imesh1 )

    else

!     subdivide into five tetrahedrons

      call intmesh_tetrahedron ( x([1,2,4,5],:), s([1,2,4,5]), imesh1 )
      call intmesh_tetrahedron ( x([2,3,4,7],:), s([2,3,4,7]), imesh2 )
      call merge_imesh ( imesh1, imesh2, imesh3 )
      call delete ( imesh1 )
      call delete ( imesh2 )
      call intmesh_tetrahedron ( x([5,8,7,4],:), s([5,8,7,4]), imesh1 )
      call merge_imesh ( imesh3, imesh1, imesh2 )
      call delete ( imesh3 )
      call delete ( imesh1 )
      call intmesh_tetrahedron ( x([2,5,6,7],:), s([2,5,6,7]), imesh1 )
      call merge_imesh ( imesh2, imesh1, imesh3 )
      call delete ( imesh2 )
      call delete ( imesh1 )
      call intmesh_tetrahedron ( x([2,4,5,7],:), s([2,4,5,7]), imesh1 )
      call merge_imesh ( imesh3, imesh1, imesh )
      call delete ( imesh3 )
      call delete ( imesh1 )

    end if

  end subroutine intmesh_hexahedron


! merge imesh1 and imesh2 to imesh

  subroutine merge_imesh ( imesh1, imesh2, imesh )

    type(imesh_t), intent(in) :: imesh1, imesh2
    type(imesh_t), intent(inout) :: imesh

    integer :: ndim, nnodes, nelem, elshape

    if ( imesh1%ndim /= imesh2%ndim ) then
      write(*,'(/a/)') &
        'Error merge_imesh: ndim of imesh1 different from imesh2 '
      stop
    end if
    if ( imesh1%elshape /= imesh2%elshape ) then
      write(*,'(/a/)') &
        'Error merge_imesh: elshape of imesh1 different from imesh2 '
      stop
    end if

    ndim = imesh1%ndim
    nnodes = imesh1%nnodes + imesh2%nnodes
    nelem = imesh1%nelem + imesh2%nelem
    elshape = imesh1%elshape

    call create ( imesh, ndim, nnodes, nelem, elshape )

    imesh%topology(:,1:imesh1%nelem) = imesh1%topology
    imesh%topology(:,imesh1%nelem+1:nelem) = imesh2%topology + imesh1%nnodes
    imesh%coor(1:imesh1%nnodes,:) = imesh1%coor
    imesh%coor(imesh1%nnodes+1:nnodes,:) = imesh2%coor

  end subroutine merge_imesh


! Jacobian for imesh (relative length or area of real element versus the
! reference element)

  subroutine imesh_jacobian ( imesh )

    type(imesh_t), intent(inout) :: imesh


    integer :: elem, i1, i2, i3
    real(dp), dimension(imesh%ndim) :: d1, d2, nrml

    if ( allocated(imesh%jacobian) ) then
      return  ! jacobian already computed
    else
      allocate ( imesh%jacobian(imesh%nelem) )
    end if

    select case(imesh%elshape)

    case(0)

!     point element

      imesh%jacobian = 0

    case(1)

!     line element

      if ( imesh%ndim == 1 ) then

!       1D

        do elem = 1, imesh%nelem
          i1 = imesh%topology(1,elem)
          i2 = imesh%topology(2,elem)
          imesh%jacobian(elem) = abs( imesh%coor(i2,1) - imesh%coor(i1,1) )
        end do

      else

!       2D and 3D

        do elem = 1, imesh%nelem
          i1 = imesh%topology(1,elem)
          i2 = imesh%topology(2,elem)
          d1 = imesh%coor(i2,:) - imesh%coor(i1,:)
          imesh%jacobian(elem) = sqrt(dot_product(d1,d1))
        end do

      end if

      imesh%jacobian = imesh%jacobian / 2   ! correct for length of
                                            ! reference element = 2

    case(3)

!     triangle

      if ( imesh%ndim == 2 ) then

!       2D

        do elem = 1, imesh%nelem
          i1 = imesh%topology(1,elem)
          i2 = imesh%topology(2,elem)
          i3 = imesh%topology(3,elem)
          d1 = imesh%coor(i2,:) - imesh%coor(i1,:)
          d2 = imesh%coor(i3,:) - imesh%coor(i1,:)
          imesh%jacobian(elem) = abs( d1(1) * d2(2) - d1(2) * d2(1) )
        end do

      else

!       3D

        do elem = 1, imesh%nelem
          i1 = imesh%topology(1,elem)
          i2 = imesh%topology(2,elem)
          i3 = imesh%topology(3,elem)
          d1 = imesh%coor(i2,:) - imesh%coor(i1,:)
          d2 = imesh%coor(i3,:) - imesh%coor(i1,:)
          nrml(1) = d1(2) * d2(3) - d1(3) * d2(2)
          nrml(2) = d1(3) * d2(1) - d1(1) * d2(3)
          nrml(3) = d1(1) * d2(2) - d1(2) * d2(1)
          imesh%jacobian(elem) = sqrt(dot_product(nrml,nrml))
        end do

      end if

    case default

      call errormsg_case_default ( 'imesh_jacobian', 'imesh%elshape', &
        int_value=imesh%elshape )

    end select

  end subroutine imesh_jacobian


! Jacobian*normal vector for imesh (relative length or area of real element
! versus the reference element * the unit normal vector)

  subroutine imesh_jacobiannormal ( imesh )

    type(imesh_t), intent(inout) :: imesh


    integer :: elem, i1, i2, i3
    real(dp), dimension(imesh%ndim) :: d1, d2, nrml

    if ( allocated(imesh%jacobiannormal) ) then
      return  ! jacobiannormal already computed
    else
      allocate ( imesh%jacobiannormal(imesh%nelem,imesh%ndim) )
    end if

    select case(imesh%elshape)

    case(0)

!     point element

      imesh%jacobiannormal = 0

    case(1)

!     line element

      if ( imesh%ndim == 2 ) then

!       2D

        do elem = 1, imesh%nelem
          i1 = imesh%topology(1,elem)
          i2 = imesh%topology(2,elem)
          d1 = imesh%coor(i2,:) - imesh%coor(i1,:)
          imesh%jacobiannormal(elem,:) = [ d1(2), -d1(1) ]
        end do

      else

        write(*,'(/2(a/))') &
          'Error imesh_jacobiannormal:', &
          ' jacobiannormal only available for 2D line elements'
        stop

      end if

      imesh%jacobiannormal = imesh%jacobiannormal / 2 ! correct for length of
                                                      ! reference element = 2

    case(3)

!     triangle

      if ( imesh%ndim == 3 ) then

!       3D

        do elem = 1, imesh%nelem
          i1 = imesh%topology(1,elem)
          i2 = imesh%topology(2,elem)
          i3 = imesh%topology(3,elem)
          d1 = imesh%coor(i2,:) - imesh%coor(i1,:)
          d2 = imesh%coor(i3,:) - imesh%coor(i1,:)
          nrml(1) = d1(2) * d2(3) - d1(3) * d2(2)
          nrml(2) = d1(3) * d2(1) - d1(1) * d2(3)
          nrml(3) = d1(1) * d2(2) - d1(2) * d2(1)
          imesh%jacobiannormal(elem,:) = nrml
        end do

      else

        write(*,'(/2(a/))') &
          'Error imesh_jacobiannormal:', &
          ' jacobiannormal only available for 3D triangular elements'
        stop

      end if

    case default

      call errormsg_case_default ( 'imesh_jacobiannormal', 'imesh%ndim', &
        int_value=imesh%ndim )

    end select

  end subroutine imesh_jacobiannormal


! Map reference coordinates to real coordinates

  subroutine imesh_map ( imesh, elem, xi, x )

    type(imesh_t), intent(inout) :: imesh

!   element number
    integer, intent(in) :: elem

!   reference coordinates xi(:,ndimxi)
    real(dp), dimension(:,:), intent(in) :: xi

!   real coordinates x(:,ndim)
    real(dp), dimension(:,:), intent(out) :: x


    real(dp), dimension(size(xi,1),imesh%elnumnod) :: phi


    select case(imesh%elshape)

    case(1)

!     line element

      phi(:,1) = (1-xi(:,1))/2
      phi(:,2) = (1+xi(:,1))/2

    case(3)

!     triangle

      phi(:,1) = 1 - xi(:,1) - xi(:,2)
      phi(:,2) = xi(:,1)
      phi(:,3) = xi(:,2)

    case default

      write(*,'(/a, i0/)') &
        'Error imesh_map: invalid elshape = ', imesh%elshape
      stop

    end select

    x = matmul ( phi, imesh%coor(imesh%topology(:,elem),:) )

  end subroutine imesh_map


! compute length or area

  function area_imesh ( imesh ) result(esum)

    type(imesh_t), intent(inout) :: imesh

    real(dp) :: esum

    call imesh_jacobian ( imesh )

    esum = sum ( imesh%jacobian )

!   correct for the size of the reference element

    select case ( imesh%elshape )
      case(1)
        esum = esum * 2
      case(3)
        esum = esum / 2
      case default
        call errormsg_case_default ( 'area_imesh', &
          'imesh%elshape', int_value=imesh%elshape )
    end select

  end function area_imesh


! count number of integration points defined by imesh

  function number_of_integration_points_imesh ( imesh, ninti ) result(numint)

    type(imesh_t), intent(inout) :: imesh

!   number of integration points on a subelement
    integer, intent(in) :: ninti

    integer :: numint

    numint = ninti * imesh%nelem

  end function number_of_integration_points_imesh


! integration points for imesh

  subroutine integration_points_imesh ( imesh, x, w, wn, ninti, xe, we, numint )

    type(imesh_t), intent(inout) :: imesh

!   the reference coordinates, weights and weights*unitnormal of the
!   integration points x(:,ndim), w(:), wn(:,ndim)
    real(dp), intent(inout), dimension(:,:) :: x
    real(dp), intent(inout), dimension(:) :: w
    real(dp), intent(inout), dimension(:,:) :: wn

!   number of integration points on a subelement
    integer, intent(in) :: ninti

!   the reference coordinates, weights of the integration points
!   for a single element, xe(ninti,ndimxi), we(ninti)
    real(dp), intent(inout), dimension(:,:) :: xe
    real(dp), intent(inout), dimension(:) :: we

!   if present: the number of integration points in imesh
    integer, intent(out), optional :: numint


    integer :: elem, ipntr, i, ip


!   compute jacobiannormal

    call imesh_jacobiannormal(imesh)


!   compute weights

    ipntr = 0

    do elem = 1, imesh%nelem

      call imesh_map ( imesh, elem, xe, x(ipntr+1:ipntr+ninti,:) )

      do i = 1, imesh%ndim
        wn(ipntr+1:ipntr+ninti,i) = we * imesh%jacobiannormal(elem,i)
      end do

      ipntr = ipntr + ninti

    end do

    do ip = 1, ipntr
      w(ip) = sqrt(dot_product(wn(ip,:),wn(ip,:)))
    end do

    if ( present(numint) ) numint = ipntr

  end subroutine integration_points_imesh


! create a mesh from imesh

  subroutine imesh_to_mesh ( imesh, mesh )

    use mesh_m, only: mesh_t
    use meshgen_extra_m, only: mesh_skeleton

    type(imesh_t), intent(inout) :: imesh

!   the basic mesh created
    type(mesh_t), intent(inout) :: mesh


!   make skeleton mesh

    call mesh_skeleton ( mesh, imesh%nnodes, imesh%nelem, imesh%elshape, &
      imesh%ndim )

!   fill topology and coor

    mesh%topology(1)%a = imesh%topology
    mesh%coor = imesh%coor

  end subroutine imesh_to_mesh

end module intmesh_m
