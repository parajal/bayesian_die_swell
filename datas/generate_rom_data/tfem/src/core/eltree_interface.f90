
! Copyright (C) 2008-2012 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines routines associated to the eltree type with respect to the "interface"

module eltree_interface_m

  use eltree_basic_m
  use meshgen_m

  implicit none


! interface for generic area subroutine

  interface area
    module procedure area_interface_eltree
  end interface area

! interface for generic number_of_interface_elements

  interface number_of_interface_elements
    module procedure number_of_interface_elements_eltree
  end interface number_of_interface_elements

! interface for generic number_of_interface_integration_points subroutine

  interface number_of_interface_integration_points
    module procedure number_of_interface_integration_points_eltree
  end interface number_of_interface_integration_points

! interface for generic interface_integration_points subroutine

  interface interface_integration_points
    module procedure interface_integration_points_eltree
  end interface interface_integration_points

contains


! count number of interface_elements in a tree

  recursive function number_of_interface_elements_eltree ( eltree ) &
    result(numsub)

    type(eltree_t), intent(in) :: eltree

    integer :: numsub

    integer :: num, inum, ndim

    numsub = 0

    if ( associated(eltree%subel) ) then

!     continue searching

      ndim = size(eltree%coor,2)
      num = 2**ndim

      do inum = 1, num
        numsub = numsub + &
                    number_of_interface_elements_eltree ( eltree%subel(inum) )
      end do

    else

!     end of the branch

      if ( associated(eltree%imesh) ) then

!       further subdivision using intmesh

        numsub = eltree%imesh%nelem

      else

!       standard case

        numsub = 0

      end if

    end if

  end function number_of_interface_elements_eltree


! count number of interface mesh nodes in a tree

  recursive function number_of_interface_mesh_nodes ( eltree ) result(numnodes)

    type(eltree_t), intent(inout) :: eltree

    integer :: numnodes

    integer :: num, inum, ndim


    numnodes = 0

    if ( associated(eltree%subel) ) then

!     continue searching

      ndim = size(eltree%coor,2)
      num = 2**ndim

      do inum = 1, num
        numnodes = numnodes + &
                      number_of_interface_mesh_nodes ( eltree%subel(inum) )
      end do

    else

!     end of the branch

      if ( associated(eltree%imesh) ) then

!       further subdivision using intmesh

        numnodes = eltree%imesh%nnodes

      else

        numnodes = 0

      end if

    end if

  end function number_of_interface_mesh_nodes


! compute length (1D), area (2D)
! NOTE: this is the reference area!

  recursive function area_interface_eltree ( eltree ) result(esum)

    type(eltree_t), intent(inout) :: eltree

    real(dp) :: esum

    integer :: num, inum, ndim


    ndim = size(eltree%coor,2)

    esum = 0

    if ( associated(eltree%subel) ) then

!     continue searching

      num = 2**ndim

      do inum = 1, num
        esum = esum + area_interface_eltree ( eltree%subel(inum) )
      end do

    else

!     end of the branch

      if ( associated(eltree%imesh) ) then

!       further subdivision using intmesh

        esum = area_imesh ( eltree%imesh )

      else

        esum = 0

      end if

    end if

  end function area_interface_eltree


! count number of interface integration points defined by a tree

  recursive function number_of_interface_integration_points_eltree ( eltree, &
    ninti ) result(numint)

    type(eltree_t), intent(inout) :: eltree

!   number of integration points on an eltree interface element (line,triangle).
    integer, intent(in) :: ninti

    integer :: numint


    integer :: num, inum, ndim


    numint = 0

    ndim = size(eltree%coor,2)

    if ( associated(eltree%subel) ) then

!     continue searching

      num = 2**ndim

      do inum = 1, num
        numint = numint + number_of_interface_integration_points_eltree ( &
                                                  eltree%subel(inum), ninti )
      end do

    else

!     end of the branch

      if ( associated(eltree%imesh) ) then

!       further subdivision using intmesh

        numint = number_of_integration_points_imesh ( eltree%imesh, ninti )

      else

!       standard case

        numint = 0

      end if

    end if

  end function number_of_interface_integration_points_eltree


! integration points for a node (recursively)

  recursive subroutine node_interface_integration_points ( eltree, ipnr, &
    x, w, wn, ninti, xe, we )

    type(eltree_t), intent(inout) :: eltree

!   start pointer for each new subelement
!   should be set to zero on first entry
    integer, intent(inout) :: ipnr

!   the reference coordinates, weights and weights*normal of the integration
!   points in the eltree interface elements (line,triangle)
!   x(:,ndim), w(:), wn(:,ndim)
    real(dp), intent(inout), dimension(:,:) :: x
    real(dp), intent(inout), dimension(:) :: w
    real(dp), intent(inout), dimension(:,:) :: wn

!   number of integration points on a single interface element
    integer, intent(in) :: ninti

!   the reference coordinates and weights of the integration points
!   for a single element in the eltree.
!   xe(ninti,ndimxi), we(ninti)
    real(dp), intent(inout), dimension(:,:) :: xe
    real(dp), intent(inout), dimension(:) :: we


    integer :: num, inum, ndim, numint


    ndim = size(eltree%coor,2)

    if ( associated(eltree%subel) ) then

!     continue searching

      num = 2**ndim

      do inum = 1, num
        call node_interface_integration_points ( eltree%subel(inum), ipnr, x, &
          w, wn, ninti, xe, we )
      end do

    else

!     end of the branch

      if ( associated(eltree%imesh) ) then

!       further subdivision using intmesh

!       integration points for imesh

        call integration_points_imesh ( eltree%imesh, x(ipnr+1:,:), &
          w(ipnr+1:), wn(ipnr+1:,:), ninti, xe, we, numint=numint )

        ipnr = ipnr + numint

      else

!       standard case

        numint = 0

      end if

    end if

  end subroutine node_interface_integration_points


! integration points for an eltree

  recursive subroutine interface_integration_points_eltree ( eltree, x, &
    w, wn, ninti, xe, we )

    type(eltree_t), intent(inout) :: eltree

!   the reference coordinates, weights and weights*normal of the integration
!   points in the eltree interface elements (line,triangle)
!   x(:,ndim), w(:), wn(:,ndim)
!   NOTE: these are with respect to the root reference element represented by
!         eltree. If mapcoor has been used the points and weights need to be
!         transformed as follows:
!
!         x' = mapcoor(x)
!
!         wn' = J F^{-T} wn
!
!         w' = |wn'|
!
!         where J=det(F) and F is the Jacobian matrix of mapcoor: dx'/dx.
    real(dp), intent(inout), dimension(:,:) :: x
    real(dp), intent(inout), dimension(:) :: w
    real(dp), intent(inout), dimension(:,:) :: wn

!   number of integration points on a single interface element
    integer, intent(in) :: ninti

!   the reference coordinates and weights of the integration points
!   for a single element in the eltree.
!   xe(ninti,ndimxi), we(ninti)
    real(dp), intent(inout), dimension(:,:) :: xe
    real(dp), intent(inout), dimension(:) :: we

    integer :: ipnr

    ipnr = 0

    call node_interface_integration_points ( eltree, ipnr, x, w, wn, ninti, &
      xe, we )

  end subroutine interface_integration_points_eltree


! create an interface mesh from an eltree

  subroutine eltree_to_interface_mesh ( eltree, mesh, mapcoor )

    use glob_defs_m

    type(eltree_t), intent(inout) :: eltree

!   the basic mesh created
!   NOTE: if mapcoor is not present the mesh is defined on the reference element
    type(mesh_t), intent(inout) :: mesh

!   if present map reference coordinates to real coordinates
!   NOTE: multiple points need to be handled at the same time, with
!   x(nnodes,ndim)
    optional :: mapcoor
    interface
      function mapcoor ( x ) result(m)
        use kind_defs_m
        implicit none
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1),size(x,2)) :: m
      end function mapcoor
    end interface


    integer :: nelem, nnodes, elshape, node, elem, ndim


    ndim = size(eltree%coor,2)

!   intmesh mesh

    nelem = number_of_interface_elements ( eltree )

    nnodes = number_of_interface_mesh_nodes ( eltree )

    select case ( ndim )
      case(2)
        elshape = 1  ! 2-node line element
      case(3)
        elshape = 3  ! 3-node triangular element
      case default
        call errormsg_case_default ( 'eltree_interface_mesh', 'ndim', &
          int_value=ndim )
    end select

    call mesh_skeleton ( mesh, nnodes, nelem, elshape, ndim )

    node = 0
    elem = 0

    call fill_coor_topology_mesh3 ( eltree, mesh, node, elem, mapcoor )

  end subroutine eltree_to_interface_mesh


! fill mesh with coordinates and topology according to intmesh

  recursive subroutine fill_coor_topology_mesh3 ( eltree, mesh, node, elem, &
    mapcoor )

    type(eltree_t), intent(inout) :: eltree
    type(mesh_t), intent(inout) :: mesh
    integer, intent(inout) :: node, elem

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


    integer :: num, inum, ndim


    ndim = size(eltree%coor,2)

    if ( associated(eltree%subel) ) then

!     continue searching

      num = 2**ndim

      do inum = 1, num
        call fill_coor_topology_mesh3 ( eltree%subel(inum), mesh, node, elem, &
          mapcoor )
      end do

    else

!     end of the branch

      if ( .not. associated(eltree%imesh) ) return ! only intmesh elements are
                                                   ! included

      if ( eltree%imesh%nelem > 0 ) then

!       add elements

        mesh%topology(1)%a(:,elem+1:elem+eltree%imesh%nelem) = &
                                            eltree%imesh%topology + node

        num = eltree%imesh%nnodes

        if ( present(mapcoor) ) then
          mesh%coor(node+1:node+num,:) = mapcoor ( eltree%imesh%coor )
        else
          mesh%coor(node+1:node+num,:) = eltree%imesh%coor
        end if

        node = node + eltree%imesh%nnodes
        elem = elem + eltree%imesh%nelem

      end if

    end if

  end subroutine fill_coor_topology_mesh3

end module eltree_interface_m
