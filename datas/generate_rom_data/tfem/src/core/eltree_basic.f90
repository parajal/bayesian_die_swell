
! Copyright (C) 2008-2026 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines the eltree type and some general associated routines.
! An eltree is either a binary tree (1D), quadtree (2D) or octtree (3D).
! The type of subelement is the same in the whole tree and can be line elements
! in 1D, triangles or quadrilaterals in 2D, tetrahedrons or hexahedrons in 3D.
! The root node of the tree is the standard reference element (for example
! [-1:1,-1:1] for quads. If the region of interest is different, a mapping
! function (mapcoor) must be used to map the reference element to real
! coordinates.
! An element (node) of the tree maybe subdivided into a number of subelements
! (nodes). For example for a quadtree each quad can be further subdivided into
! four subquads etc.:
!
!                         ___o___              o: subdivided node
!                        /  / \  \             +: node at the end of a branche
!                       /  +   +  \         !
!                   ___o___     ___o___
!                  /  / \  \   /  / \  \    !
!                 +  o   +  + +  +   +  +
!                 ___|___
!                /  / \  \                  !
!               +  +   +  +
!
! The nodes at the end of a branch (+'s) make up a "mesh".

module eltree_basic_m

  use glob_defs_m
  use submesh_m
  use intmesh_m
  use set_optional_m

  implicit none


! type definition of a node of an element tree

  type eltree_t

!   defines the sign of the "levelset" in the vertices of the subelement
!    +1   all larger than zero
!     0   different signs and/or at least one zero
!    -1   all less than zero
    integer :: lsign = 0

!   coor(2,ndim)
!   coordinate array defining the size and position of the element in the
!   reference element. For line, quad and hexahedron elements coor defines
!   the coordinates of the first vertex (minimal coordinates) and the diagonal
!   positioned vertex (maximal coordinates). For example a line and a quad:
!               ----x
!      x---x    |   |
!               x----
!   For triangles and tetrahedrons the coordinates represent the quad and
!   hexahedron they are part of (see conf below).
    real(dp), dimension(:,:), allocatable :: coor

!   define the element type, 1: triangle/tetrahedron
!                            2: quadrilateral/hexahedron (default)
    integer :: eltype = 2

!   conf defines the configuration of the triangle/tetrahedron in the square or
!   cubic reference element. For triangles, a square with vertices [1,2,3,4]
!   located at [(0,0),(1,0),(1,1),(0,1)] resp. is divided into two regions:
!
!         4\                   4---------------3 [1,1]
!         |  \                  \              |
!         |    \                  \  conf=2    |
!         |      \                  \          |
!         |        \                  \        |
!         |          \                  \      |
!         | conf=1     \                  \    |
!         |              \                  \  |
!   [0,0] 1---------------2                   \2
!
!
!   | conf | vertices |
!   |-----------------|
!   |  1   | [1,2,4]  |
!   |-----------------|
!   |  2   | [2,3,4]  |
!   -------------------
!
!
!   For tetrahedrons, a cube with vertices [1,2,3,4,5,6,7,8] located at
!   [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)] resp.
!   is divided into:
!   | conf |  vertices  |
!   |-------------------|
!   |  1   | [1,2,3,7]  |
!   |-------------------|
!   |  2   | [1,4,7,8]  |
!   ---------------------
!   |  3   | [1,3,4,7]  |
!   |-------------------|
!   |  4   | [1,5,6,7]  |
!   ---------------------
!   |  5   | [1,2,6,7]  |
!   |-------------------|
!   |  6   | [1,5,7,8]  |
!   ---------------------
!   The reason for choosing these tetrahedrons is that they can be subdivided
!   into 4 subtetrahedrons of the same configuration and 4 subtetrahedrons of
!   one of the other 5 possible configurations.
!   NOTE: For triangles the root element always has to have conf = 1 since this
!         yields the same reference element for isoparametric coordinates
!         For tetrahedrons, the root element always has to have conf = 2, since
!         a mapping is performed for this reference element to the reference
!         element for isoparametric coordinates.
!   NOTE2: For quads and hexahedrons 'conf' is not used.
    integer :: conf = 1

!   If associated a further subdivision of the node/subelement is defined
!   based on a further level in the tree.
!   NOTE: the current implementation uses a pointer to an array of eltrees.
!   This means that individual branches cannot be replaced/added. This makes
!   the tree rather static, but it is suitable for our purposes. The advantage
!   is that the same type can be used for binary, quad and octtrees.
    type(eltree_t), dimension(:), pointer :: subel => null()

!   If associated a further subdivision of the node/subelement is defined
!   based on the module submesh (lines, triangles or tetrahedrons).
!   This defines the volume mesh only and can be present independently from
!   imesh.
!   NOTE: this is the (volume) submesh and does not contain the interface mesh.
!   NOTE: this can only be present if:
!     1) the node is at the end of a branch
!     2) lsign=0
    type(smesh_t), pointer :: smesh => null()

!   if associated a further subdivision of the node/subelement is defined
!   based on the module intmesh (lines, triangles or tetrahedrons).
!   This defines the interface mesh only and can be present independently from
!   smesh.
!   NOTE: this does NOT contain the (volume) submesh, only the interface mesh.
!   NOTE: this can only be present if:
!     1) the node is at the end of a branch
!     2) lsign=0
    type(imesh_t), pointer :: imesh => null()

  end type eltree_t


! pointer to an eltree encapsulated in a new type to make an array of pointers
! to an eltree, like:
!   type(eltree_p), dimension(:), allocatable :: eltree_array
! a pointer p(i) is referenced by eltree_array(i)%p

  type eltree_p
    type(eltree_t), pointer :: p => null()
  end type eltree_p


! interface for generic delete subroutine

  interface delete
    module procedure delete_eltree
  end interface delete

contains


! fill a node with the coordinates of the element and optionally
! the element type and configuration

  subroutine fill_node_eltree ( eltree, coor, conf, eltype )

    type(eltree_t), intent(inout) :: eltree
    real(dp), dimension(:,:), intent(in) :: coor

    integer, optional, intent(in) :: conf, eltype

    integer :: ndim

    ndim = size(coor,2)

    allocate(eltree%coor(2,ndim))

    eltree%coor = coor
    eltree%eltype = set_optional ( variable=eltype, default=2 )
    eltree%conf = set_optional ( variable=conf, default=1 )

  end subroutine fill_node_eltree


! subdivide a node into subelements

  subroutine subdivide_node_eltree ( eltree )

    type(eltree_t), intent(inout) :: eltree

    real(dp), dimension(2,size(eltree%coor,2)) :: coor
    real(dp) :: h(size(eltree%coor,2))
    integer :: ndim, num, i, j, k

!   h_all(i) * h yields the new first coordinate of the i-th subelement
    real(dp), dimension(:,:), allocatable :: h_all

!   conf_all(i) gives the the configuration of the i-th subelement
    integer, dimension(:), allocatable :: conf_all

!   the element type number  1:triangle/tetrahedron, 2:quadrilateral/hexahedron
    integer :: eltype

!   subdivide a node into subelements

    eltype = eltree%eltype

    ndim = size(coor,2)
    num = 2**ndim

    allocate(eltree%subel(num)) ! associate the pointer for the subelements

    h = ( eltree%coor(2,:) - eltree%coor(1,:) ) / 2  ! half-size vector

!   fill coordinates of the sub nodes for different space dimension and element
!   types

    select case ( ndim )

      case(1) ! line

        do i = 0, 1
          coor(1,1) = eltree%coor(1,1) + h(1) * i
          coor(2,:) = coor(1,:) + h
          call fill_node_eltree ( eltree%subel(i+1), coor )
        end do

      case(2)

        if ( eltype == 1 ) then ! triangle

          allocate(h_all(4,2))
          allocate(conf_all(4))

!         the subdivision for triangles depends on the configuration, note that
!         not only the new coordinates have to be given, but also the
!         configuration

          select case ( eltree%conf ) ! determine the configuration
            case(1)
              h_all = reshape( [ 0, 0, 1, 0, 0, 0, 0, 1 ] , [ 4, 2 ] )
              conf_all = [ 1, 2, 1, 1 ]
            case(2)
              h_all = reshape( [ 1, 0, 1, 1, 0, 1, 1, 1 ] , [ 4, 2 ] )
              conf_all = [ 2, 2, 1, 2 ]
            case default
              call errormsg_case_default ( 'subdivide_node_eltree', &
                'eltree%conf', int_value=eltree%conf )
          end select

          do i=1,4
            coor(1,:) = eltree%coor(1,:) + h * h_all(i,:)
            coor(2,:) = coor(1,:) + h * [1, 1]
            call fill_node_eltree ( eltree%subel(i), coor, conf_all(i), eltype )
          end do

          deallocate(h_all,conf_all)

        else if ( eltype == 2 ) then ! quad

          do j = 0, 1
            do i = 0, 1
              coor(1,:) = eltree%coor(1,:) + h * [i, j]
              coor(2,:) = coor(1,:) + h
              call fill_node_eltree ( eltree%subel(i+2*j+1), coor )
            end do
          end do

        end if

      case(3)

        if ( eltype == 1 ) then ! tetrahedron

          allocate(h_all(8,3))
          allocate(conf_all(8))

!         the subdivision for tets depends on the configuration, note that
!         not only the new coordinates have to be given, but also the
!         configuration

          select case ( eltree%conf ) ! determine the configuration
            case(1)
              h_all = reshape( [ 0, 1, 1, 1, 1, 1, 1, 1, &
                                 0, 0, 0, 0, 1, 1, 1, 1, &
                                 0, 0, 0, 0, 0, 0, 0, 1 ] , [ 8, 3 ] )
              conf_all = [ 1, 2, 3, 1, 4, 1, 5, 1 ]
            case(2)
              h_all = reshape( [ 0, 0, 0, 0, 0, 0, 0, 1, &
                                 0, 1, 1, 1, 1, 1, 1, 1, &
                                 0, 0, 0, 0, 1, 1, 1, 1 ] , [ 8, 3 ] )
              conf_all = [ 2, 2, 6, 4, 2, 3, 1, 2 ]
            case(3)
              h_all = reshape( [ 0, 0, 0, 0, 1, 1, 1, 1, &
                                 0, 1, 1, 1, 1, 1, 1, 1, &
                                 0, 0, 0, 0, 0, 0, 0, 1 ] , [ 8, 3 ] )
              conf_all = [ 3, 5, 3, 1, 6, 2, 3, 3 ]
            case(4)
              h_all = reshape( [ 0, 0, 0, 0, 1, 1, 1, 1, &
                                 0, 0, 0, 0, 0, 0, 0, 1, &
                                 0, 1, 1, 1, 1, 1, 1, 1 ] , [ 8, 3 ] )
              conf_all = [ 4, 1, 5, 4, 6, 2, 4, 4 ]
            case(5)
              h_all = reshape( [ 0, 1, 1, 1, 1, 1, 1, 1, &
                                 0, 0, 0, 0, 0, 0, 0, 1, &
                                 0, 0, 0, 0, 1, 1, 1, 1 ] , [ 8, 3 ] )
              conf_all = [ 5, 6, 4, 5, 3, 1, 5, 5 ]
            case(6)
              h_all = reshape( [ 0, 0, 0, 0, 0, 0, 0, 1, &
                                 0, 0, 0, 0, 1, 1, 1, 1, &
                                 0, 1, 1, 1, 1, 1, 1, 1 ] , [ 8, 3 ] )
              conf_all = [ 6, 6, 2, 3, 6, 4, 5, 6 ]
            case default
              call errormsg_case_default ( 'subdivide_node_eltree', &
                'eltree%conf', int_value=eltree%conf )
          end select

          do i=1,8
            coor(1,:) = eltree%coor(1,:) + h * h_all(i,:)
            coor(2,:) = coor(1,:) + h * [ 1, 1, 1 ]
            call fill_node_eltree ( eltree%subel(i), coor, conf_all(i), eltype )
          end do

          deallocate(h_all,conf_all)

        else if ( eltype == 2 ) then ! hexahedron

          do k = 0, 1
            do j = 0, 1
              do i = 0, 1
                coor(1,:) = eltree%coor(1,:) + h * [ i, j, k ]
                coor(2,:) = coor(1,:) + h
                call fill_node_eltree ( eltree%subel(i+2*j+4*k+1), coor )
              end do
            end do
          end do

        end if

      case default

        call errormsg_case_default ( 'subdivide_node_eltree', 'ndim', &
          int_value=ndim )

    end select

  end subroutine subdivide_node_eltree


! Subdivide starting from a given node (coordinates are known)

  recursive subroutine subdivide ( eltree, levelset, numsplit, split_threshold,&
    numsplitmin, mapcoor, submesh, intmesh )

    use limits_m, only: ZEROLEVELSET

    type(eltree_t), intent(inout) :: eltree

!   the levelset function, where levelset=0 defines the interface where elements
!   need to be refined.
!   NOTE: only the sign (and the zero value) is actually used. Hence the
!         function can be as simple as +1 or -1 depending whether you are
!         "in" or "out". However, if submesh=.true. the levelset function must
!         be properly defined (smooth) near the value 0: at least
!         differentiable if the interface is smooth, since linear interpolation
!         is used to find the position of the interface.
!         It does not have to be a signed distance function.
!   NOTE: if mapcoor is not present the entry coordinates x are reference
!         coordinates in the element starting from the root of the tree
!   NOTE: multiple points need to be handled at the same time, with
!         x(nnodes,ndim)
    interface
      function levelset ( x )
        use kind_defs_m
        implicit none
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1)) :: levelset
      end function levelset
    end interface

!   number of tree splits in the branches of the smallest elements near the
!   levelset=0 This means that the smallest elements are m=2^numsplit smaller
!   than the reference element. An array of m (1D), m x m (2D) and
!   m x m x m (3D) elements would completely fill the reference element.
!   Examples (2D):
!     numsplit=0: no split; the element size stays the same
!     numsplit=1: one split; the smallest element size is such that
!                 2 x 2 elements would fill the reference element
!     numsplit=2: two splits; the smallest element size is such that
!                 4 x 4 elements would fill the reference element
!     numsplit=10: ten splits; the smallest element size is such that
!                 1024 x 1024 elements would fill the reference element
    integer, intent(in) :: numsplit

!   if present, this parameter determines which subelements are considered to be
!   near the interface levelset=0 and are split according to the value of
!   numsplit.
!   Subelements where all values of the levelset function phi in the vertices
!   fullfill:
!          phi>split_threshold or phi<-split_threshold
!   are considered to be "far from the interface" and are not split according to
!   the value of numsplit. All others are split.
!   default=0
    real(dp), optional, intent(in) :: split_threshold

!   minimum number of tree splits for the complete reference element.
!   The reason for including this is the following. If the levelset function has
!   the same sign in the vertex nodes, there will be no further subdividing.
!   However it is possible that there is a structure smaller than the size of
!   the element (for example an a circle smaller than the element or a local
!   "river" going through the element). In that case you can demand that
!   subdividing is done anyway to "find" the structure.
!   default = 0 (no imposed split).
    integer, intent(in), optional :: numsplitmin

!   if present map reference coordinates to real coordinates
!   NOTE: multiple points need to be handled at the same time, with
!   x(nnodes,ndim)
    optional :: mapcoor
    interface
      function mapcoor ( x )
        use kind_defs_m
        implicit none
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1),size(x,2)) :: mapcoor
      end function mapcoor
    end interface

!   if present and .true. a subelement at the end of branch that
!   is intersected by the zero levelset interface is further subdivided
!   using the module submesh (subdivision into lines, triangles, tetrahedrons
!   which are aligned with the zero levelset.
!   This parameter defines the volume mesh.
!   NOTE: the levelset must be (smooth) near the value 0: at least
!   differentiable
!         if the interface is smooth, since linear interpolation is used to find
!         the position of the interface.
!         It does not have to be a signed distance function.
!   default = .false.
    logical, intent(in), optional :: submesh

!   if present and .true. a subelement at the end of branch that
!   is intersected by the zero levelset interface is further subdivided
!   using the module intmesh (subdivision into lines, triangles, tetrahedrons
!   which are aligned with the zero levelset.
!   This parameter defines the interface mesh.
!   NOTE: the levelset must be (smooth) near the value 0: at least
!         differentiable if the interface is smooth, since linear
!         interpolation is used to find the position of the interface.
!         It does not have to be a signed distance function.
!   default = .false.
    logical, intent(in), optional :: intmesh

!   Subdivide starting from a given node (coordinates are known)
!   The starting node is usually the root of the tree, but any node that still
!   has no subdivions associated can be used.

!   levelset value in the vertices of the (sub)element
    real(dp), dimension(:), allocatable :: lvlset

!   reference coordinates of the vertices of the (sub)element
    real(dp), dimension(:,:), allocatable :: x

    real(dp) :: lepsmin, href, leps, lthres
    real(dp), dimension(size(eltree%coor,2)) :: h
    integer :: ndim, num, inum, i, j, k, lnumsplitmin
    logical :: lsubmesh, lintmesh
!   the element type number  1:triangle/tetrahedron, 2:quadrilateral/hexahedron
    integer :: eltype


    eltype = eltree%eltype

!   size vector of the element

    h = eltree%coor(2,:) - eltree%coor(1,:)

    ndim = size(eltree%coor,2)

!   determine the size of the reference element and number of vertices

    select case ( ndim )

      case(1) ! line

        href = 2
        allocate( lvlset(2) )
        allocate( x(2,1) )

      case(2)

        if ( eltype == 1 ) then ! triangle
          href = 1
          allocate( lvlset(3) )
          allocate( x(3,2) )
        else if ( eltype == 2 ) then ! quad
          href = 2
          allocate( lvlset(4) )
          allocate( x(4,2) )
        end if

      case(3)

        if ( eltype == 1 ) then ! tetrahedron
          href = 1
          allocate( lvlset(4) )
          allocate( x(4,3) )
        else if ( eltype == 2 ) then ! hexahedron
          href = 2
          allocate( lvlset(8) )
          allocate( x(8,3) )
        end if

      case default

        call errormsg_case_default ( 'subdivide', 'ndim', int_value=ndim )

    end select

!   find reference coordinates in the vertices of the element
!   NOTE: for triangles and tets the numbering of the vertices is
!         not arbitrary. The numbering is such that the vertices [1,2,3] and
!         [1,2,3,4] for triangles and tets resp. yields shapes with the correct
!         sense for further subdivision using intmesh and submesh.

    select case ( ndim )

      case(1) ! line

        x = eltree%coor

      case(2)

        if ( eltype == 2 ) then ! quad

          do j = 0, 1
            do i = 0, 1
              x(i+2*j+1,:) = eltree%coor(1,:) + h * [ i, j ]
            end do
          end do

        else if ( eltype == 1 ) then ! triangle

          select case ( eltree%conf )
            case(1)
              x(1,:) = eltree%coor(1,:) + h * [ 0, 0 ]
              x(2,:) = eltree%coor(1,:) + h * [ 1, 0 ]
              x(3,:) = eltree%coor(1,:) + h * [ 0, 1 ]
            case(2)
              x(1,:) = eltree%coor(1,:) + h * [ 1, 0 ]
              x(2,:) = eltree%coor(1,:) + h * [ 1, 1 ]
              x(3,:) = eltree%coor(1,:) + h * [ 0, 1 ]
            case default
              call errormsg_case_default ( 'subdivide', 'ndim', int_value=ndim )
          end select

        end if

      case(3)

        if ( eltype == 2 ) then ! hexahedron

          do k = 0, 1
            do j = 0, 1
              do i = 0, 1
                x(i+2*j+4*k+1,:) = eltree%coor(1,:) + h * [ i, j, k ]
              end do
            end do
          end do

        else if ( eltype == 1 ) then ! tetrahedron

!         NOTE: because the reference shape for tets in the eltree is different
!               than the isoparametric one, an additional mapping must be done.

          select case ( eltree%conf ) ! determine the configuration
            case(1)
              x(1,:) = mapref ( eltree%coor(1,:) + h * [ 0, 0, 0] )
              x(2,:) = mapref ( eltree%coor(1,:) + h * [ 1, 0, 0] )
              x(3,:) = mapref ( eltree%coor(1,:) + h * [ 1, 1, 0] )
              x(4,:) = mapref ( eltree%coor(2,:) )
            case(2)
              x(1,:) = mapref ( eltree%coor(1,:) + h * [ 0, 0, 0] )
              x(2,:) = mapref ( eltree%coor(2,:) )
              x(3,:) = mapref ( eltree%coor(1,:) + h * [ 0, 1, 0] )
              x(4,:) = mapref ( eltree%coor(1,:) + h * [ 0, 1, 1] )
            case(3)
              x(1,:) = mapref ( eltree%coor(1,:) + h * [ 0, 0, 0] )
              x(2,:) = mapref ( eltree%coor(1,:) + h * [ 1, 1, 0] )
              x(3,:) = mapref ( eltree%coor(1,:) + h * [ 0, 1, 0] )
              x(4,:) = mapref ( eltree%coor(2,:) )
            case(4)
              x(1,:) = mapref ( eltree%coor(1,:) + h * [ 0, 0, 0] )
              x(2,:) = mapref ( eltree%coor(1,:) + h * [ 1, 0, 1] )
              x(3,:) = mapref ( eltree%coor(2,:) )
              x(4,:) = mapref ( eltree%coor(1,:) + h * [ 0, 0, 1] )
            case(5)
              x(1,:) = mapref ( eltree%coor(1,:) + h * [ 0, 0, 0] )
              x(2,:) = mapref ( eltree%coor(1,:) + h * [ 1, 0, 0] )
              x(3,:) = mapref ( eltree%coor(2,:) )
              x(4,:) = mapref ( eltree%coor(1,:) + h * [ 1, 0, 1] )
            case(6)
              x(1,:) = mapref ( eltree%coor(1,:) + h * [ 0, 0, 0] )
              x(4,:) = mapref ( eltree%coor(2,:) )
              x(2,:) = mapref ( eltree%coor(1,:) + h * [ 0, 1, 1] )
              x(3,:) = mapref ( eltree%coor(1,:) + h * [ 0, 0, 1] )
            case default
              call errormsg_case_default ( 'subdivide', 'eltree%conf', &
                int_value=eltree%conf )
          end select

        end if

      case default

        call errormsg_case_default ( 'subdivide', 'ndim', int_value=ndim )

    end select

!   find the real global coordinates (if mapcoor is present) and determine the
!   levelset, otherwise detemine the levelset for the reference coordinates
    if ( present(mapcoor) ) then
      lvlset = levelset( mapcoor(x) )
    else
      lvlset = levelset(x)
    end if


!   cut-off small value
    if ( ZEROLEVELSET > 0._dp ) then
      where ( abs(lvlset) < ZEROLEVELSET )
        lvlset = 0
      end where
    end if

    lthres = set_optional ( variable=split_threshold, default=0._dp )
    lnumsplitmin = set_optional ( variable=numsplitmin, default=0 )
    lepsmin = href / 2**lnumsplitmin + 1e-12_dp
    leps = href / 2**numsplit + 1e-12_dp

!   Determine further subdivision

    if ( all(lvlset>lthres) .and. h(1) <= lepsmin ) then

!     element full in the + region (no further subdivision)

      eltree%lsign = 1

      deallocate(lvlset, x)

    else if ( all(lvlset<-lthres) .and. h(1) <= lepsmin ) then

!     element full in the - region (no further subdivision)

      eltree%lsign = -1
      deallocate(lvlset, x)

    else if ( h(1) <= leps .and. h(1) <= lepsmin ) then

!     element is small enough (no further subdivision along the branch)

      if ( all(lvlset>0) ) then

        eltree%lsign = 1

      else if ( all(lvlset<0) ) then

        eltree%lsign = -1

      else

        eltree%lsign = 0

        lsubmesh = set_optional ( variable = submesh, default = .false. )
        lintmesh = set_optional ( variable = intmesh, default = .false. )

        if ( lsubmesh ) then

!         further subdivide using submesh

          allocate(eltree%smesh)

          if ( ndim == 1 ) then ! line
            call submesh_line_element ( x, lvlset, eltree%smesh )
          else if ( ndim == 2 .and. eltype == 1 ) then ! triangle
            call submesh_triangle ( x, lvlset, eltree%smesh )
          else if ( ndim == 2 .and. eltype == 2 ) then ! quad
            call submesh_quadrilateral ( x([1,2,4,3],:), &
                                lvlset([1,2,4,3]), eltree%smesh )
          else if ( ndim == 3 .and. eltype == 1 ) then !tet
            call submesh_tetrahedron ( x, lvlset, eltree%smesh )
          else if ( ndim == 3 .and. eltype == 2  ) then !hexa
            call submesh_hexahedron ( x([1,2,4,3,5,6,8,7],:), &
                                lvlset([1,2,4,3,5,6,8,7]), eltree%smesh )
          end if

        end if

        if ( lintmesh ) then

!         further subdivide using intmesh

          allocate(eltree%imesh)

          if ( ndim == 1 ) then ! line
            call intmesh_line_element ( x, lvlset, eltree%imesh )
          else if ( ndim == 2 .and. eltype == 1 ) then ! triangle
            call intmesh_triangle ( x, lvlset, eltree%imesh )
          else if ( ndim == 2 .and. eltype == 2 ) then ! quad
           call intmesh_quadrilateral ( x([1,2,4,3],:), &
                             lvlset([1,2,4,3]), eltree%imesh )
          else if ( ndim == 3 .and. eltype == 1 ) then !tet
            call intmesh_tetrahedron ( x, lvlset, eltree%imesh )
          else if ( ndim == 3 .and. eltype == 2  ) then !hexa
            call intmesh_hexahedron ( x([1,2,4,3,5,6,8,7],:), &
                             lvlset([1,2,4,3,5,6,8,7]), eltree%imesh )
          end if

        end if

        deallocate(lvlset, x)

      end if

    else

!     further subdivision of the element

      if ( all(lvlset>0) ) then

        eltree%lsign = 1

      else if ( all(lvlset<0) ) then

        eltree%lsign = -1

      else

        eltree%lsign = 0

      end if

      deallocate(lvlset, x)

!     fill coordinates of the sub elements

      call subdivide_node_eltree ( eltree )

      num = 2**ndim

!     subdivide further each element of the subdivision

      do inum = 1, num
        call subdivide ( eltree%subel(inum), levelset, numsplit, &
          split_threshold, numsplitmin, mapcoor, submesh, intmesh )
      end do

    end if

  end subroutine subdivide


! delete eltree from memory

  recursive subroutine delete_eltree ( eltree )

    type(eltree_t), intent(inout) :: eltree

    integer :: num, inum, ndim

    if ( associated(eltree%subel) ) then

!     continue searching

      ndim = size(eltree%coor,2)
      num = 2**ndim

      do inum = 1, num
        call delete_eltree ( eltree%subel(inum) )
      end do

      deallocate(eltree%subel)

      deallocate(eltree%coor)

    else

!     end of the branch

      deallocate(eltree%coor)

      if ( associated(eltree%smesh) ) then
        call delete ( eltree%smesh )
        deallocate(eltree%smesh)
      end if

      if ( associated(eltree%imesh) ) then
        call delete ( eltree%imesh )
        deallocate(eltree%imesh)
      end if

    end if

  end subroutine delete_eltree


! mapping for tetrahedrons from the eltree reference region to the tfem
! reference region. The mapping is performed by two deformation gradient
! matrices:
! F1 = [1  0  0]   and F2 = [ 1  0  0]
!      [0  1 -1]            [ 0  1  0]
!      [0  0  1]            [-1  0 -1]
!
! NOTE: this is the mapping for conf = 2.

  function mapref ( x )

    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x,1)) :: mapref

!   F2 . (F1 . x) is written out explicitly

    mapref = [ x(1), x(2)-x(3), x(3)-x(1) ]

  end function mapref

end module eltree_basic_m
