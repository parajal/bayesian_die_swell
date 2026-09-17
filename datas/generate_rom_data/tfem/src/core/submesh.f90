
! Copyright (C) 2008-2022 Martien A. Hulsen
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
! The levelset function is defined by the values in the vertex points and is
! linearly interpolated along the sides of the element.

module submesh_m

  use glob_defs_m
  use set_optional_m

  implicit none

! simple mesh type

  type smesh_t

    integer :: ndim     = 0        ! space dimension
    integer :: nnodes   = 0        ! number of nodal points
    integer :: nelem    = 0        ! number of elements

!   shape of the element
!      elshape=1 : 2 node line element
!      elshape=3 : 3 node triangle
!      elshape=11: 4 node tetrahedron
    integer :: elshape  = 0
    integer :: elnumnod   = 0     ! number of nodal points of a single element

!   connections of elements to nodal points
!   topology(node,elem) gives the global node number of the local node
!   in element elem
    integer, allocatable, dimension(:,:) :: topology

!   coordinates coor(nnodes,ndim)
    real(dp), allocatable, dimension(:,:) :: coor

!   lsign(elem)
!   defines the sign of the "levelset" in the nodes of an element
!    +1   all larger than or equal to zero (except all zero)
!     0   all are zero
!    -1   all smaller than or equal to zero (except all zero)
!   NOTE: the definition is different from eltree!
!   NOTE: the value of 0 should normally not happen and probably means
!         there is something wrong with the definition of the levelset function.
    integer, allocatable, dimension(:) :: lsign

!   jacobian(elem): jacobian of each element
    real(dp), allocatable, dimension(:) :: jacobian

  end type smesh_t


! interface for generic create subroutine

  interface create
    module procedure create_smesh
  end interface create

! interface for generic delete subroutine

  interface delete
    module procedure delete_smesh
  end interface delete

! interface for generic volume subroutine

  interface volume
    module procedure volume_smesh
  end interface volume

! interface for generic number_of_subelements subroutine

  interface number_of_subelements
    module procedure number_of_subelements_smesh
  end interface number_of_subelements

! interface for generic number_of_subelements subroutine (vector)

  interface number_of_subelements_vector
    module procedure number_of_subelements_vector_smesh
  end interface number_of_subelements_vector

! interface for generic number_of_integration_points subroutine

  interface number_of_integration_points
    module procedure number_of_integration_points_smesh
  end interface number_of_integration_points

! interface for generic integration_points subroutine

  interface integration_points
    module procedure integration_points_smesh
  end interface integration_points


contains


! create smesh

  subroutine create_smesh ( smesh, ndim, nnodes, nelem, elshape )

    type(smesh_t), intent(inout) :: smesh
    integer, intent(in) :: ndim, nnodes, nelem, elshape

    smesh%ndim = ndim
    smesh%nnodes = nnodes
    smesh%nelem = nelem
    smesh%elshape = elshape

    select case(elshape)
    case(1)
!     line element
      smesh%elnumnod = 2
    case(3)
!     triangle
      smesh%elnumnod = 3
    case(11)
!     tetrahedron
      smesh%elnumnod = 4
    case default
      write(*,'(/a,i0/)') &
        'Error create_smesh: invalid elshape = ', elshape
      stop
    end select

    allocate(smesh%topology(smesh%elnumnod,nelem))
    allocate(smesh%coor(nnodes,ndim))
    allocate(smesh%lsign(nelem))

  end subroutine create_smesh


! delete smesh

  subroutine delete_smesh ( smesh )

    type(smesh_t), intent(out) :: smesh

  end subroutine delete_smesh


! subdivide a line element into line elements according to the levelset
! function s

  subroutine submesh_line_element ( x, s, smesh )

!   coordinates of the end points
    real(dp), dimension(:,:), intent(in) :: x

!   levelset function of the end points
    real(dp), dimension(:), intent(in) :: s

!   the "volume" mesh
    type(smesh_t), intent(inout) :: smesh


!   counter of zero or positive level set
    integer :: countz, countp
!   interpolation xi
    real(dp) :: xi
!   all nodes negative of positive
    logical :: alln, allp


    alln = all ( s < 0._dp )
    allp = all ( s > 0._dp )
    countz = count ( s == 0._dp )

    if ( alln .or. allp .or. countz == 2 ) then

!     one line element only

      call create ( smesh, ndim=size(x,2), nnodes=2, nelem=1, elshape=1 )

      smesh%coor = x
      smesh%topology(:,1) = [ 1, 2 ]
      if ( alln ) then
        smesh%lsign = -1
      else if ( allp ) then
        smesh%lsign = 1
      else if ( countz == 2 ) then
        smesh%lsign = 0
      end if

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
      call create ( smesh, ndim=size(x,2), nnodes=2, nelem=1, elshape=1 )

      smesh%coor = x
      smesh%topology(:,1) = [ 1, 2 ]

      countp = count ( s > 0._dp )

      if ( countp == 1 ) then
!       one positive vertex
        smesh%lsign = 1
      else if ( countp == 0 ) then
!       one negative vertex
        smesh%lsign = -1
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

      call create ( smesh, ndim=size(x,2), nnodes=3, nelem=2, elshape=1 )

      smesh%coor(1:2,:) = x ! vertex nodes

!     linear interpolate to find the point on the zero levelset

      xi = s(1) / ( s(1) - s(2) )

      smesh%coor(3,:) = (1-xi) * x(1,:) + xi * x(2,:)

!     topology

      smesh%topology(:,1) = [ 1, 3 ]
      smesh%topology(:,2) = [ 3, 1 ]

!     inside or outside?

      if ( s(1) < 0._dp ) then

        smesh%lsign(1) = -1
        smesh%lsign(2) = 1

      else

        smesh%lsign(1) = 1
        smesh%lsign(2) = -1

      end if

    end if

  end subroutine submesh_line_element


! subdivide a triangle into subtriangles according to the levelset function s

  subroutine submesh_triangle ( x, s, smesh )

!   coordinates of the vertices
    real(dp), dimension(:,:), intent(in) :: x

!   levelset function of the vertices
    real(dp), dimension(:), intent(in) :: s

!   the "volume" mesh
    type(smesh_t), intent(inout) :: smesh


!   the indices for the generic nodes 1, 2, 3
    integer :: i1, i2, i3
!   the "other" two nodes
    integer, parameter :: io(2,3)=reshape([2,3,3,1,1,2],[2,3])
!   number of negative or zero values of the levelset
    integer :: countn, countz
!   interpolation xi between the sides 1-2 and 1-3
    real(dp) :: xi12, xi13, xi23
!   all vertices negative of positive
    logical :: alln, allp, cutvertexn, cutvertexp, cutvertex, cutedgen, &
      cutedgep, cutedge

    alln = all ( s < 0._dp )
    allp = all ( s > 0._dp )
    countz = count ( s == 0._dp )
    countn = count ( s < 0._dp )

!   degenerate case: cut one vertex with a zero levelset
    cutvertexn = countz == 1 .and. countn == 2
    cutvertexp = countz == 1 .and. countn == 0
    cutvertex = cutvertexn .or. cutvertexp

!   degenerate case: cut one edge with a zero levelset (see intmesh)
    cutedgen = countz == 2 .and. countn == 1
    cutedgep = countz == 2 .and. countn == 0
    cutedge = cutedgen .or. cutedgep

!   choose correct case

    if ( alln .or. allp .or. countz == 3 .or. cutvertex .or. cutedge ) then

!     one triangle only

      call create ( smesh, ndim=size(x,2), nnodes=3, nelem=1, elshape=3 )

      smesh%coor = x
      smesh%topology(:,1) = [ 1, 2, 3 ]
      if ( alln .or. cutvertexn .or. cutedgen ) then
        smesh%lsign = -1
      else if ( allp .or. cutvertexp .or. cutedgep ) then
        smesh%lsign = 1
      else if ( countz == 3 ) then
        smesh%lsign = 0
      end if

    else if ( countz == 1 .and. countn == 1 ) then

!     two triangles; degenerate case
!
!     We have the following generic case:
!
!                3
!                |\          !
!                | \         !
!                |  \        !
!                |   4
!                |  / \      !
!                | /   \     !
!                |/     \    !
!                1-------2
!
!     The nodes 1, 2 and 3 correspond with vertices of the original triangle
!     (possibly permutated).
!     Here, 2 is the node that is positive (negative) and 3 is
!     negative (positive). The line between the nodes 1-4 is the zero levelset.
!     The position of node 4 is is found by linear interpolation between the
!     levelset values in the vertex nodes 2 and 3.
!

      call create ( smesh, ndim=size(x,2), nnodes=4, nelem=2, elshape=3 )

      smesh%coor(1:3,:) = x ! vertex nodes

      i1 = maxval ( minloc(abs(s)) )

!     find the other two nodes

      i2 = io(1,i1)
      i3 = io(2,i1)

!     linear interpolate to find the points on the zero levelset

      xi23 = s(i2) / ( s(i2) - s(i3) )

      smesh%coor(4,:) = (1-xi23)*smesh%coor(i2,:) + xi23*smesh%coor(i3,:)

!     topology

      smesh%topology(:,1) = [ i1, i2,  4 ]
      smesh%topology(:,2) = [ i1,  4, i3 ]

!     inside or outside?

      if ( s(i2) > 0._dp ) then

!       i2 on positive side

        smesh%lsign(1) = 1
        smesh%lsign(2) = -1

      else if ( s(i2) < 0._dp ) then

!       i2 on negative side

        smesh%lsign(1) = -1
        smesh%lsign(2) = 1

      else
        stop 'internal error submesh_triangle 1'
      end if

    else

!     Three triangles.
!
!     We have the following generic case:
!
!                3
!               /|\          !
!              / | \         !
!             /  |  \        !
!            5   |   \       !
!           / \  |    \      !
!          /   \ |     \     !
!         /     \|      \    !
!        1-------4-------2
!
!     The nodes 1, 2 and 3 correspond with vertices of the original triangle
!     (possibly permutated).
!     Here, 1 is the node that is positive (negative) and 2 and 3 are
!     negative (positive). The line between the nodes 4-5 is the zero levelset
!     which is found by linear interpolation between the levelset values in the
!     vertex nodes.
!

      call create ( smesh, ndim=size(x,2), nnodes=5, nelem=3, elshape=3 )

      smesh%coor(1:3,:) = x ! vertex nodes

      countn = count ( s < 0._dp )

      if ( countn == 1 ) then

!       one negative vertex

        i1 = maxval ( minloc (s) )

      else if ( countn == 2 ) then

!       one positive vertex

        i1 = maxval ( maxloc (s) )

      else
        stop 'internal error submesh_triangle 2'
      end if

!     find the other two nodes

      i2 = io(1,i1)
      i3 = io(2,i1)

!     linear interpolate to find the points on the zero levelset

      xi12 = s(i1) / ( s(i1) - s(i2) )
      xi13 = s(i1) / ( s(i1) - s(i3) )

      smesh%coor(4,:) = (1-xi12)*smesh%coor(i1,:) + xi12*smesh%coor(i2,:)
      smesh%coor(5,:) = (1-xi13)*smesh%coor(i1,:) + xi13*smesh%coor(i3,:)

!     topology

      smesh%topology(:,1) = [ i1,  4,  5 ]
      smesh%topology(:,2) = [  4, i2, i3 ]
      smesh%topology(:,3) = [  5,  4, i3 ]

!     inside or outside?

      if ( countn == 1 ) then

!       one negative vertex

        smesh%lsign(1) = -1
        smesh%lsign(2) = 1
        smesh%lsign(3) = 1

      else if ( countn == 2 ) then

!       one positive vertex

        smesh%lsign(1) = 1
        smesh%lsign(2) = -1
        smesh%lsign(3) = -1

      end if

    end if

  end subroutine submesh_triangle


! subdivide a quad into subtriangles according to the levelset function s

  subroutine submesh_quadrilateral ( x, s, smesh )

!   coordinates of the vertices
    real(dp), dimension(:,:), intent(in) :: x

!   levelset function of the vertices
    real(dp), dimension(:), intent(in) :: s

    type(smesh_t), intent(inout) :: smesh


    type(smesh_t) :: smesh1, smesh2

!   subdivide into two triangles
!
!       4-----------3
!       |        __/|
!       |     __/   |
!       |  __/      |
!       | /         |
!       1-----------2
!
    call submesh_triangle ( x([1,2,3],:), s([1,2,3]), smesh1 )
    call submesh_triangle ( x([1,3,4],:), s([1,3,4]), smesh2 )

    call merge_smesh ( smesh1, smesh2, smesh )

    call delete ( smesh1 )
    call delete ( smesh2 )

  end subroutine submesh_quadrilateral


! subdivide a tetrahedron into subtets according to the levelset function s

  subroutine submesh_tetrahedron ( x, s, smesh )

!   coordinates of the vertices
    real(dp), dimension(:,:), intent(in) :: x

!   levelset function of the vertices
    real(dp), dimension(:), intent(in) :: s

    type(smesh_t), intent(inout) :: smesh


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
!   number of negative, zero values of the levelset
    integer :: countn, countz
!   interpolation xi between the sides 1-2, 1-4, 2-4, 2-3, 1-3 and 3-4
    real(dp) :: xi12, xi14, xi24, xi23, xi13, xi34
!   all vertices negative of positive
    logical :: alln, allp
!   cut-off vertex with a zero levelset
    logical :: cutvertex, cutvertexn, cutvertexp
!   cut-off edge with a zero levelset
    logical :: cutedge, cutedgen, cutedgep
!   cut-off face with a zero levelset
    logical :: cutface, cutfacen, cutfacep
!   debugging
    logical, parameter :: debug = .false.


    alln = all ( s < 0._dp )
    allp = all ( s > 0._dp )
    countz = count ( s == 0._dp )
    countn = count ( s < 0._dp )
    cutvertex = countz == 1 .and. ( countn == 0 .or. countn == 3 )
    cutedge = countz == 2 .and. ( countn == 0 .or. countn == 2 )

!   degenerate case: cut one vertex with a zero levelset
    cutvertexn = countz == 1 .and. countn == 3
    cutvertexp = countz == 1 .and. countn == 0
    cutvertex = cutvertexn .or. cutvertexp

!   degenerate case: cut one edge with a zero levelset
    cutedgen = countz == 2 .and. countn == 2
    cutedgep = countz == 2 .and. countn == 0
    cutedge = cutedgen .or. cutedgep

!   degenerate case: cut one face with a zero levelset (see intmesh)
    cutfacen = countz == 3 .and. countn == 1
    cutfacep = countz == 3 .and. countn == 0
    cutface = cutfacen .or. cutfacep

    if ( alln .or. allp .or. cutvertex .or. cutedge .or. cutface .or. &
         countz == 4  ) then

!     one tetrahedron only

      call create ( smesh, ndim=size(x,2), nnodes=4, nelem=1, elshape=11 )

      smesh%coor = x
      smesh%topology(:,1) = [ 1, 2, 3, 4 ]
      if ( alln .or. cutvertexn .or. cutedgen .or. cutfacen ) then
        smesh%lsign = -1
      else if ( allp .or. cutvertexp .or. cutedgep .or. cutfacep ) then
        smesh%lsign = 1
      else if ( countz == 4 ) then
        smesh%lsign = 0
      end if

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

      call create ( smesh, ndim=size(x,2), nnodes=6, nelem=3, elshape=11 )

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

      smesh%coor(1:4,:) = x ! vertex nodes
      smesh%coor(5,:) = (1-xi23)*x(i2,:) + xi23*x(i3,:)
      smesh%coor(6,:) = (1-xi24)*x(i2,:) + xi24*x(i4,:)

!     topology

      smesh%topology(:,1) = [ i1, i2,  5,  6 ]
      smesh%topology(:,2) = [ i1,  5, i3,  6 ]
      smesh%topology(:,3) = [ i1,  6, i3, i4 ]

      if ( countn == 1 ) then

!       s(i2) < 0

        smesh%lsign(1) = -1
        smesh%lsign(2:3) = 1

      else if ( countn == 2 ) then

!       s(i2) > 0

        smesh%lsign(1) = 1
        smesh%lsign(2:3) = -1

      else
        stop 'internal error submesh_tetrahedron 1'
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

      call create ( smesh, ndim=size(x,2), nnodes=5, nelem=2, elshape=11 )

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

      smesh%coor(1:4,:) = x ! vertex nodes
      smesh%coor(5,:) = (1-xi34)*x(i3,:) + xi34*x(i4,:)

!     topology

      smesh%topology(:,1) = [ i1, i2, i3,  5 ]
      smesh%topology(:,2) = [ i1, i2,  5, i4 ]

!     inside or outside?

      if ( s(i3) > 0._dp ) then

!       s(i3) > 0

        smesh%lsign(1) = 1
        smesh%lsign(2) = -1

      else if ( s(i3) < 0._dp ) then

!       s(i3) < 0

        smesh%lsign(1) = -1
        smesh%lsign(2) = 1

      else
        stop 'internal error submesh_tetrahedron 2'
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

      call create ( smesh, ndim=size(x,2), nnodes=7, nelem=4, elshape=11 )

      smesh%coor(1:4,:) = x ! vertex nodes

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

      smesh%coor(5,:) = (1-xi12)*smesh%coor(i1,:) + xi12*smesh%coor(i2,:)
      smesh%coor(6,:) = (1-xi13)*smesh%coor(i1,:) + xi13*smesh%coor(i3,:)
      smesh%coor(7,:) = (1-xi14)*smesh%coor(i1,:) + xi14*smesh%coor(i4,:)

!     topology

      smesh%topology(:,1) = [ i1,  5,  6, 7  ]
      smesh%topology(:,2) = [  5,  6,  7, i4 ]
      smesh%topology(:,3) = [ i2, i3,  5, i4 ]
      smesh%topology(:,4) = [  5, i3,  6, i4 ]

!     inside or outside?

      if ( countn == 1 ) then

!       one negative vertex

        smesh%lsign(1) = -1
        smesh%lsign(2:4) = 1

      else if ( countn == 3 ) then

!       one positive vertex

        smesh%lsign(1) = 1
        smesh%lsign(2:4) = -1

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

      call create ( smesh, ndim=size(x,2), nnodes=8, nelem=6, elshape=11 )

      smesh%coor(1:4,:) = x ! vertex nodes

!     find negative nodes

      num = 0
      do i = 1, 4
        if ( s(i) < 0._dp ) then
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

      smesh%coor(5,:) = (1-xi14)*smesh%coor(i1,:) + xi14*smesh%coor(i4,:)
      smesh%coor(6,:) = (1-xi24)*smesh%coor(i2,:) + xi24*smesh%coor(i4,:)
      smesh%coor(7,:) = (1-xi23)*smesh%coor(i2,:) + xi23*smesh%coor(i3,:)
      smesh%coor(8,:) = (1-xi13)*smesh%coor(i1,:) + xi13*smesh%coor(i3,:)

!     topology

      smesh%topology(:,1) = [ i1, i2,  8, 6  ]
      smesh%topology(:,2) = [ i2,  7,  8, 6  ]
      smesh%topology(:,3) = [ i1,  8,  5, 6  ]
      smesh%topology(:,4) = [ i3,  8,  7, 6  ]
      smesh%topology(:,5) = [ i3,  8,  6, 5  ]
      smesh%topology(:,6) = [  5,  6, i3, i4 ]

!     inside or outside?

      smesh%lsign(1:3) = -1
      smesh%lsign(4:6) = 1

    else
        stop 'internal error submesh_tetrahedron 3'
    end if

  end subroutine submesh_tetrahedron


! subdivide a hexahedron into subtets according to the levelset function s

  subroutine submesh_hexahedron ( x, s, smesh )

    use limits_m

!   coordinates of the vertices
    real(dp), dimension(:,:), intent(in) :: x

!   levelset function of the vertices
    real(dp), dimension(:), intent(in) :: s

    type(smesh_t), intent(inout) :: smesh


    type(smesh_t) :: smesh1, smesh2, smesh3


    if ( SUBDIVIDE_HEX6 ) then

!     subdivide into six tetrahedrons

      call submesh_tetrahedron ( x([1,2,3,5],:), s([1,2,3,5]), smesh1 )
      call submesh_tetrahedron ( x([2,3,5,7],:), s([2,3,5,7]), smesh2 )
      call merge_smesh ( smesh1, smesh2, smesh3 )
      call delete ( smesh1 )
      call delete ( smesh2 )
      call submesh_tetrahedron ( x([2,7,5,6],:), s([2,7,5,6]), smesh1 )
      call merge_smesh ( smesh3, smesh1, smesh2 )
      call delete ( smesh3 )
      call delete ( smesh1 )
      call submesh_tetrahedron ( x([1,5,3,8],:), s([1,5,3,8]), smesh1 )
      call merge_smesh ( smesh2, smesh1, smesh3 )
      call delete ( smesh2 )
      call delete ( smesh1 )
      call submesh_tetrahedron ( x([3,5,7,8],:), s([3,5,7,8]), smesh1 )
      call merge_smesh ( smesh3, smesh1, smesh2 )
      call delete ( smesh3 )
      call delete ( smesh1 )
      call submesh_tetrahedron ( x([1,3,4,8],:), s([1,3,4,8]), smesh1 )
      call merge_smesh ( smesh2, smesh1, smesh )
      call delete ( smesh2 )
      call delete ( smesh1 )

    else

!     subdivide into five tetrahedrons

      call submesh_tetrahedron ( x([1,2,4,5],:), s([1,2,4,5]), smesh1 )
      call submesh_tetrahedron ( x([2,3,4,7],:), s([2,3,4,7]), smesh2 )
      call merge_smesh ( smesh1, smesh2, smesh3 )
      call delete ( smesh1 )
      call delete ( smesh2 )
      call submesh_tetrahedron ( x([5,8,7,4],:), s([5,8,7,4]), smesh1 )
      call merge_smesh ( smesh3, smesh1, smesh2 )
      call delete ( smesh3 )
      call delete ( smesh1 )
      call submesh_tetrahedron ( x([2,5,6,7],:), s([2,5,6,7]), smesh1 )
      call merge_smesh ( smesh2, smesh1, smesh3 )
      call delete ( smesh2 )
      call delete ( smesh1 )
      call submesh_tetrahedron ( x([2,4,5,7],:), s([2,4,5,7]), smesh1 )
      call merge_smesh ( smesh3, smesh1, smesh )
      call delete ( smesh3 )
      call delete ( smesh1 )

    end if

  end subroutine submesh_hexahedron


! merge smesh1 and smesh2 to smesh

  subroutine merge_smesh ( smesh1, smesh2, smesh )

    type(smesh_t), intent(in) :: smesh1, smesh2
    type(smesh_t), intent(inout) :: smesh

    integer :: ndim, nnodes, nelem, elshape

    if ( smesh1%ndim /= smesh2%ndim ) then
      write(*,'(/a/)') &
        'Error merge_smesh: ndim of smesh1 different from smesh2 '
      stop
    end if
    if ( smesh1%elshape /= smesh2%elshape ) then
      write(*,'(/a/)') &
        'Error merge_smesh: elshape of smesh1 different from smesh2 '
      stop
    end if

    ndim = smesh1%ndim
    nnodes = smesh1%nnodes + smesh2%nnodes
    nelem = smesh1%nelem + smesh2%nelem
    elshape = smesh1%elshape

    call create ( smesh, ndim, nnodes, nelem, elshape )

    smesh%topology(:,1:smesh1%nelem) = smesh1%topology
    smesh%topology(:,smesh1%nelem+1:nelem) = smesh2%topology + smesh1%nnodes
    smesh%coor(1:smesh1%nnodes,:) = smesh1%coor
    smesh%coor(smesh1%nnodes+1:nnodes,:) = smesh2%coor
    smesh%lsign(1:smesh1%nelem) = smesh1%lsign
    smesh%lsign(smesh1%nelem+1:nelem) = smesh2%lsign

  end subroutine merge_smesh


! Jacobian for smesh (relative length, area or volume of real element versus the
! reference element)

  subroutine smesh_jacobian ( smesh )

    type(smesh_t), intent(inout) :: smesh


    integer :: elem, i1, i2, i3, i4
    real(dp) :: det
    real(dp), dimension(smesh%ndim) :: d1, d2, d3, nrml

    if ( allocated(smesh%jacobian) ) then
      return  ! jacobian already computed
    else
      allocate ( smesh%jacobian(smesh%nelem) )
    end if

    select case(smesh%elshape)

    case(1)

!     line element

      if ( smesh%ndim == 1 ) then

!       1D

        do elem = 1, smesh%nelem
          i1 = smesh%topology(1,elem)
          i2 = smesh%topology(2,elem)
          smesh%jacobian(elem) = abs( smesh%coor(i2,1) - smesh%coor(i1,1) )
        end do

      else

!       2D and 3D

        do elem = 1, smesh%nelem
          i1 = smesh%topology(1,elem)
          i2 = smesh%topology(2,elem)
          d1 = smesh%coor(i2,:) - smesh%coor(i1,:)
          smesh%jacobian(elem) = sqrt(dot_product(d1,d1))
        end do

      end if

      smesh%jacobian = smesh%jacobian / 2   ! correct for length of reference
                                            ! element = 2

    case(3)

!     triangle

      if ( smesh%ndim == 2 ) then

!       2D

        do elem = 1, smesh%nelem
          i1 = smesh%topology(1,elem)
          i2 = smesh%topology(2,elem)
          i3 = smesh%topology(3,elem)
          d1 = smesh%coor(i2,:) - smesh%coor(i1,:)
          d2 = smesh%coor(i3,:) - smesh%coor(i1,:)
          smesh%jacobian(elem) = abs( d1(1) * d2(2) - d1(2) * d2(1) )
        end do

      else

!       3D

        do elem = 1, smesh%nelem
          i1 = smesh%topology(1,elem)
          i2 = smesh%topology(2,elem)
          i3 = smesh%topology(3,elem)
          d1 = smesh%coor(i2,:) - smesh%coor(i1,:)
          d2 = smesh%coor(i3,:) - smesh%coor(i1,:)
          nrml(1) = d1(2) * d2(3) - d1(3) * d2(2)
          nrml(2) = d1(3) * d2(1) - d1(1) * d2(3)
          nrml(3) = d1(1) * d2(2) - d1(2) * d2(1)
          smesh%jacobian(elem) = sqrt(dot_product(nrml,nrml))
        end do

      end if

    case(11)

!     tetrahedron

      do elem = 1, smesh%nelem
        i1 = smesh%topology(1,elem)
        i2 = smesh%topology(2,elem)
        i3 = smesh%topology(3,elem)
        i4 = smesh%topology(4,elem)
        d1 = smesh%coor(i2,:) - smesh%coor(i1,:)
        d2 = smesh%coor(i3,:) - smesh%coor(i1,:)
        d3 = smesh%coor(i4,:) - smesh%coor(i1,:)
        det = d1(1)*d2(2)*d3(3) - d1(1)*d2(3)*d3(2) - d2(1)*d1(2)*d3(3) &
            + d2(1)*d1(3)*d3(2) + d3(1)*d1(2)*d2(3) - d3(1)*d1(3)*d2(2)
        smesh%jacobian(elem) = abs(det)
      end do

    case default

      call errormsg_case_default ( 'smesh_jacobian', 'smesh%elshape', &
        int_value=smesh%elshape )

    end select

  end subroutine smesh_jacobian


! Map reference coordinates to real coordinates

  subroutine smesh_map ( smesh, elem, xi, x )

    type(smesh_t), intent(inout) :: smesh

!   element number
    integer, intent(in) :: elem

!   reference coordinates xi(:,ndimxi)
    real(dp), dimension(:,:), intent(in) :: xi

!   real coordinates x(:,ndim)
    real(dp), dimension(:,:), intent(out) :: x


    real(dp), dimension(size(xi,1),smesh%elnumnod) :: phi


    select case(smesh%elshape)

    case(1)

!     line element

      phi(:,1) = (1-xi(:,1))/2
      phi(:,2) = (1+xi(:,1))/2

    case(3)

!     triangle

      phi(:,1) = 1 - xi(:,1) - xi(:,2)
      phi(:,2) = xi(:,1)
      phi(:,3) = xi(:,2)

    case(11)

!     tetrahedron

      phi(:,1) = 1 - xi(:,1) - xi(:,2) - xi(:,3)
      phi(:,2) = xi(:,1)
      phi(:,3) = xi(:,2)
      phi(:,4) = xi(:,3)

    case default

      call errormsg_case_default ( 'smesh_map', 'smesh%elshape', &
        int_value=smesh%elshape )

    end select

    x = matmul ( phi, smesh%coor(smesh%topology(:,elem),:) )

  end subroutine smesh_map


! count number of nodes in smesh

  function number_of_nodes_smesh ( smesh, lsign ) result(nnodes)

    type(smesh_t), intent(inout) :: smesh

!   if present include elements that have an smesh%lsign value given by lsign.
!   Example lsign=(/1/) will only include elements that have an smesh%lsign
!   value of 1.
!   default = all possible values (/-1,0,1/)
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign


    integer :: nnodes

    integer :: elem, node
    integer, allocatable, dimension(:) :: work

    allocate ( work(smesh%nnodes) )

    if ( present(lsign) ) then

!     find elements

      work = 0

      do elem = 1, smesh%nelem

        if ( all ( lsign /= smesh%lsign(elem) ) ) cycle

        work ( smesh%topology(:,elem) ) = 1  ! set nodes present to 1

      end do

!     find nodes

      nnodes = 0

      do node = 1, smesh%nnodes

        if ( work(node) == 0 ) cycle ! node not present

        nnodes = nnodes + 1
        work(node) = nnodes  ! index array old node -> new node
        !work1(nnodes) = node ! index array new node -> old node

      end do

    else

!     include all elements

      nnodes = smesh%nnodes

    end if

    deallocate ( work )

  end function number_of_nodes_smesh


! count number of elements in a submesh

  function number_of_subelements_smesh ( smesh, lsign ) result(numsub)

    type(smesh_t), intent(in) :: smesh

!   if present include elements that have an smesh%lsign value given by lsign.
!   Example lsign=(/1/) will only include elements that have an smesh%lsign
!   value of 1.
!   default = all possible values (/-1,0,1/)
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign

    integer :: numsub, elem

    if ( present(lsign) ) then
      numsub = 0
      do elem = 1, smesh%nelem
        if ( all( lsign /= smesh%lsign(elem) ) ) cycle
        numsub = numsub + 1
      end do
    else
      numsub = smesh%nelem
    end if

  end function number_of_subelements_smesh


! count number of elements in a submesh (vector)

  function number_of_subelements_vector_smesh ( smesh, lsign ) result(numsub)

    type(smesh_t), intent(in) :: smesh

!   count elements that have an smesh%lsign value given by lsign (vector value)
!   Example lsign=(/-1,1/) will count elements that have an smesh%lsign
!   value of -1 or 1. The result is a vector of the same size as lsign.
!   default = all possible values (/-1,0,1/)
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:) :: lsign

    integer, dimension(size(lsign)) :: numsub

    integer :: elem

    numsub = 0
    do elem = 1, smesh%nelem
      where ( lsign == smesh%lsign(elem) )
        numsub = numsub + 1
      end where
    end do

  end function number_of_subelements_vector_smesh


! compute length, area or volume

  function volume_smesh ( smesh, lsign ) result(esum)

    type(smesh_t), intent(inout) :: smesh

!   if present include elements that have an smesh%lsign value given by lsign.
!   Example lsign=(/1/) will only include elements that have an smesh%lsign
!   value of 1.
!   default = all possible values (/-1,0,1/)
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign

    real(dp) :: esum

    integer :: elem

    call smesh_jacobian ( smesh )

    if ( present(lsign) ) then
      esum = 0
      do elem = 1, smesh%nelem
        if ( all( lsign /= smesh%lsign(elem) ) ) cycle
        esum = esum + smesh%jacobian(elem)
      end do
    else
      esum = sum ( smesh%jacobian )
    end if

!   correct for the size of the reference element

    select case ( smesh%elshape )
      case(1)
        esum = esum * 2
      case(3)
        esum = esum / 2
      case(11)
        esum = esum / 6
      case default
        call errormsg_case_default ( 'volume_smesh', 'smesh%elshape', &
          int_value=smesh%elshape )
    end select

  end function volume_smesh


! count number of integration points defined by smesh

  function number_of_integration_points_smesh ( smesh, ninti, lsign, epsjac ) &
    result(numint)

    type(smesh_t), intent(inout) :: smesh

!   number of integration points on a subelement
    integer, intent(in) :: ninti

!   if present include elements that have an smesh%lsign value given by lsign.
!   Example lsign=(/1/) will only include elements that have an smesh%lsign
!   value of 1.
!   default = all possible values (/-1,0,1/)
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign

!   if present: elements with a jacobian smaller than epsjac will be omitted
    real(dp), intent(in), optional :: epsjac


    integer :: numint, elem


!   compute jacobian

    if ( present(epsjac) ) call smesh_jacobian(smesh)


!   choose depending on the presence of lsign or epsjac

    if ( present(lsign) .and. .not. present(epsjac) ) then

      numint = 0

      do elem = 1, smesh%nelem

        if ( all( lsign /= smesh%lsign(elem) ) ) cycle

        numint = numint + ninti

      end do

    else if ( present(lsign) .and. present(epsjac) ) then

      numint = 0

      do elem = 1, smesh%nelem

        if ( all( lsign /= smesh%lsign(elem) ) ) cycle

        if ( smesh%jacobian(elem) <= epsjac ) cycle

        numint = numint + ninti

      end do

    else if ( .not. present(lsign) .and. present(epsjac) ) then

      numint = 0

      do elem = 1, smesh%nelem

        if ( smesh%jacobian(elem) <= epsjac ) cycle

        numint = numint + ninti

      end do

    else

      numint = ninti * smesh%nelem

    end if

  end function number_of_integration_points_smesh


! integration points for smesh

  subroutine integration_points_smesh ( smesh, x, w, ninti, xe, we, lsign, &
    epsjac, numint )

    type(smesh_t), intent(inout) :: smesh

!   the reference coordinates and weights of the integration points
!   x(:,ndim), w(:)
    real(dp), intent(inout), dimension(:,:) :: x
    real(dp), intent(inout), dimension(:) :: w

!   number of integration points on a subelement
    integer, intent(in) :: ninti

!   the reference coordinates and weights of the integration points
!   for a single element, xe(ninti,ndimxi), we(ninti)
    real(dp), intent(inout), dimension(:,:) :: xe
    real(dp), intent(inout), dimension(:) :: we

!   if present include elements that have an smesh%lsign value given by lsign.
!   Example lsign=(/1/) will only include elements that have an smesh%lsign
!   value of 1.
!   default = all possible values (/-1,0,1/)
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign

!   if present: elements with a jacobian smaller than epsjac will be omitted
    real(dp), intent(in), optional :: epsjac

!   if present: the number of integration points in smesh
    integer, intent(out), optional :: numint


    integer :: elem, ipntr


!   compute jacobian

    call smesh_jacobian(smesh)


!   choose depending on the presence of lsign or epsjac

    if ( present(lsign) .and. .not. present(epsjac) ) then

      ipntr = 0

      do elem = 1, smesh%nelem

        if ( all( lsign /= smesh%lsign(elem) ) ) cycle

        call fill_x_w

      end do

    else if ( present(lsign) .and. present(epsjac) ) then

      ipntr = 0

      do elem = 1, smesh%nelem

        if ( all( lsign /= smesh%lsign(elem) ) ) cycle

        if ( smesh%jacobian(elem) <= epsjac ) cycle

        call fill_x_w

      end do

    else if ( .not. present(lsign) .and. present(epsjac) ) then

      ipntr = 0

      do elem = 1, smesh%nelem

        if ( smesh%jacobian(elem) <= epsjac ) cycle

        call fill_x_w

      end do

    else

      call fill_x_w

    end if

    if ( present(numint) ) numint = ipntr

  contains

    subroutine fill_x_w

      call smesh_map ( smesh, elem, xe, x(ipntr+1:ipntr+ninti,:) )

      w(ipntr+1:ipntr+ninti) = we * smesh%jacobian(elem)

      ipntr = ipntr + ninti

    end subroutine fill_x_w

  end subroutine integration_points_smesh


! create a smesh from smesh

  subroutine smesh_to_smesh ( smesh1, smesh, lsign )

!   input smesh
    type(smesh_t), intent(in) :: smesh1

!   output smesh
    type(smesh_t), intent(out) :: smesh

!   if present include elements that have an smesh%lsign value given by lsign.
!   Example lsign=(/1/) will only include elements that have an smesh%lsign
!   value of 1.
!   default = all possible values (/-1,0,1/)
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign


    integer :: nelem, nnodes, elem, kelem, node
    integer, allocatable, dimension(:) :: work, work1


    allocate ( work(smesh1%nnodes), work1(smesh1%nnodes) )


!   make smesh

    if ( present(lsign) ) then

!     find elements

      work = 0
      nelem = 0

      do elem = 1, smesh1%nelem

        if ( all ( lsign /= smesh1%lsign(elem) ) ) cycle

        nelem = nelem + 1

        work ( smesh1%topology(:,elem) ) = 1  ! set nodes present to 1

      end do

!     find nodes

      nnodes = 0

      do node = 1, smesh1%nnodes

        if ( work(node) == 0 ) cycle ! node not present

        nnodes = nnodes + 1
        work(node) = nnodes  ! index array old node -> new node
        work1(nnodes) = node ! index array new node -> old node

      end do

    else

!     include all elements

      nelem = smesh1%nelem
      nnodes = smesh1%nnodes

    end if

    call create_smesh ( smesh, smesh1%ndim, nnodes, nelem, smesh1%elshape )

!   fill topology and coor

    if ( present(lsign) ) then

!     find elements that are present

      kelem = 0

      do elem = 1, smesh1%nelem

        if ( all ( lsign /= smesh1%lsign(elem) ) ) cycle

        kelem = kelem + 1

        smesh%topology(:,kelem) = work(smesh1%topology(:,elem))

      end do

      smesh%coor = smesh1%coor(work1(1:nnodes),:)

    else

!     all elements

      smesh%topology = smesh1%topology
      smesh%coor = smesh1%coor

    end if

    deallocate ( work, work1 )

  end subroutine smesh_to_smesh


! create a mesh from smesh

  subroutine smesh_to_mesh ( smesh, mesh, lsign )

    use mesh_m, only: mesh_t
    use meshgen_extra_m, only: mesh_skeleton

    type(smesh_t), intent(inout) :: smesh

!   the basic mesh created
    type(mesh_t), intent(inout) :: mesh

!   if present include elements that have an smesh%lsign value given by lsign.
!   Example lsign=(/1/) will only include elements that have an smesh%lsign
!   value of 1.
!   default = all possible values (/-1,0,1/)
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign


    type(smesh_t) :: smesh1


!   make skeleton mesh

    call smesh_to_smesh ( smesh, smesh1, lsign )

    call mesh_skeleton ( mesh, smesh1%nnodes, smesh1%nelem, smesh1%elshape, &
      smesh1%ndim )


!   fill topology and coor

    mesh%topology(1)%a = smesh1%topology
    mesh%coor = smesh1%coor


!   remove temporary storage

    call delete(smesh1)


  end subroutine smesh_to_mesh

end module submesh_m
