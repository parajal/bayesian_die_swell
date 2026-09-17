
! Copyright (C) 2005-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines routine for objects and mesh

module meshgen_objects_m

  use glob_defs_m
  use mesh_m
  use shapefunc_m
  use array_defs_m
  use set_optional_m

  implicit none

contains


! Find reference coordinates for the points in objects

  subroutine find_refcoor_objects ( mesh, object1, object2, updaterefcoor2, &
    mapcoor, mapcoornr, onlyshapenod, blend )

    type(mesh_t), intent(inout) :: mesh

!   if these are present only objects object1,...,object2 are computed.
!   If only object1 is present one object is computed.
    integer, intent(in), optional :: object1, object2

!   if .false. the second connection (refcoor2) is not updated
!   default: updaterefcoor2 = .true.
    logical, intent(in), optional :: updaterefcoor2

!   if mapcoor is present the coordinates are mapped to a new position before
!   the intersection with the mesh is computed.
!   NOTE: the coordinates of the object (coor) and, if present, the coordinates
!   of the integration points (coor_int) are unaffected by the mapping. Only the
!   actual coordinates used for the intersection are mapped before the
!   intersection is performed.
    optional :: mapcoor
    interface
      function mapcoor ( nr, x ) result(m)
        use kind_defs_m
        implicit none
        integer, intent(in) :: nr
        real(dp), intent(in), dimension(:) :: x
        real(dp), dimension(size(x)) :: m
      end function mapcoor
    end interface

!   The number mapcoornr can be used to switch between different mappings in
!   the function mapcoor
    integer, intent(in), optional :: mapcoornr

!   optional parameter to indicate that only numshapenod nodes are
!   used in the geometrical shape of the element. This affects only elements
!   with additional internal nodes:
!      elshape=7 : 7 node triangle        (numshapenod=6)
!      elshape=9 : 5 node quadrilateral   (numshapenod=4)
!      elshape=10: 4 node triangle        (numshapenod=3)
!      elshape=15: 14 node tetrahedron    (numshapenod=10)
!      elshape=16: 15 node tetrahedron    (numshapenod=10)
!   default is onlyshapenod = .false. (=all nodes)
!   Reasons for not using all nodes of the element for the computation of the
!   shape include:
!    1) The internal nodes are not used for the actual unknowns that determine
!       the shape via an isoparametric element procedure and therefore must be
!       excluded in the geometrical shape mapping. For example elshape=10 is
!       used for a diffusion equation with unknowns only in the vertices (4
!       unknowns). The internal node is used for a different (discontinuous)
!       physical quantity.
!    2) Using the internal (bubble) functions in the shape can affect the
!       convergence of the Newton iteration for finding the reference
!       coordinates in a negative way. Excluding the internal nodes from the
!       geometrical shape mapping gives better convergence.
!       WARNING: depending on how the coordinates of the internal nodes are
!       generated the reference coordinates found can be different from the
!       reference coordinates for the full shape!
    logical, intent(in), optional :: onlyshapenod

!   the mesh chosen for the finding the reference coordinates
!   blend=0  get reference coordinates in the main mesh elements (mesh%element)
!   blend>0  get reference coordinates in the blend mesh elements given by
!            mesh%element_blend(:,blend),
!   default=0
    integer, intent(in), optional :: blend


!   Find reference coordinates for the points in objects


    logical :: upref2
    integer :: object, obj1, obj2
    integer :: ninti, ndim, elnumnod, nelem, elem
    real(dp), dimension(:,:), allocatable :: phi, x
    logical, dimension(:), allocatable :: onlynodes, onlynodes2


!   which objects?

    if ( present(object1) .and. present(object2) ) then
!     specified range of objects only
      obj1 = object1
      obj2 = object2
    else if ( present(object1) ) then
!     one object only
      obj1 = object1
      obj2 = object1
    else
!     all objects
      obj1 = 1
      obj2 = mesh%nobjects
    end if

    if ( obj1 < 1 .or. obj1 > mesh%nobjects .or.  &
         obj2 < 1 .or. obj2 > mesh%nobjects ) then

      write(*,'(/2a/2(a,i0),/a,i0/)') &
        'Error: object1 and/or object2 in the heading of ', &
        'find_refcoor_objects is ', &
        'out of range. object1 is ', obj1, '; object2 is ', obj2, &
        'whereas the number of objects is ', mesh%nobjects
      stop

    end if

!   update refcoor2 ?

    upref2 = set_optional ( variable=updaterefcoor2, default=.true. )

!   onlynodes

    if ( any( mesh%objects(obj1:obj2)%nodeset > 0 ) ) then
      allocate ( onlynodes(mesh%nnodes) )
    else
      allocate ( onlynodes(0) )
    end if

    if ( upref2 .and. any( mesh%objects(obj1:obj2)%nodeset2 > 0 ) ) then
      allocate ( onlynodes2(mesh%nnodes) )
    else
      allocate ( onlynodes2(0) )
    end if

!   compute reference coordinates in objects

    do object = obj1, obj2

      if ( mesh%objects(object)%nodeset > 0 ) then
        onlynodes = .false.
        onlynodes(mesh%nodesets(mesh%objects(object)%nodeset)%a) = .true.
      end if

      call find_refcoor_points ( mesh, mesh%objects(object)%groups, &
        mesh%objects(object)%nnodes, mesh%objects(object)%coor, &
        mesh%objects(object)%grpelm, mesh%objects(object)%refcoor, &
        onlynodes, mapcoor, mapcoornr, onlyshapenod, blend )

      if ( mesh%objects(object)%typeofobject == 2 .and. upref2 ) then

        if ( mesh%objects(object)%nodeset2 > 0 ) then
          onlynodes2 = .false.
          onlynodes2(mesh%nodesets(mesh%objects(object)%nodeset2)%a) = .true.
        end if

        call find_refcoor_points ( mesh, mesh%objects(object)%groups2, &
          mesh%objects(object)%nnodes, mesh%objects(object)%coor, &
          mesh%objects(object)%grpelm2, mesh%objects(object)%refcoor2, &
          onlynodes2, mapcoor, mapcoornr, onlyshapenod, blend )

      end if

      if ( mesh%objects(object)%intpoints ) then

!       also compute integration points intersection

        ninti = mesh%objects(object)%ninti
        ndim  = mesh%ndim
        elnumnod = mesh%objects(object)%elnumnod

        allocate ( phi(ninti,elnumnod), x(elnumnod,ndim) )

        call set_shape_function ( mesh%objects(object)%element%globalshape, &
                                  mesh%objects(object)%xig, phi )

        nelem = mesh%objects(object)%nelem

        do elem = 1, nelem

!         get nodal coordinates of element in object
          call get_coordinates_object ( mesh, elem, x, object )

!         compute integration point coordinates using an isoparametric element
          call isoparametric_coordinates ( x, phi, &
                                    mesh%objects(object)%coor_int(:,:,elem) )

          call find_refcoor_points ( mesh, mesh%objects(object)%groups, &
            ninti, mesh%objects(object)%coor_int(:,:,elem), &
            mesh%objects(object)%grpelm_int(:,:,elem), &
            mesh%objects(object)%refcoor_int(:,:,elem), &
            onlynodes, mapcoor, mapcoornr, onlyshapenod, blend )

          if ( mesh%objects(object)%typeofobject == 2 .and. upref2 ) then
            call find_refcoor_points ( mesh, mesh%objects(object)%groups2, &
              ninti, mesh%objects(object)%coor_int(:,:,elem), &
              mesh%objects(object)%grpelm2_int(:,:,elem), &
              mesh%objects(object)%refcoor2_int(:,:,elem), &
              onlynodes2, mapcoor, mapcoornr, onlyshapenod, blend )
          end if

        end do

        deallocate ( phi, x )

      end if

    end do

    deallocate ( onlynodes, onlynodes2 )

  end subroutine find_refcoor_objects


! Find reference coordinates for an array of points

  subroutine find_refcoor_points ( mesh, groups, nnodes, coor, grpelm, &
    refcoor, onlynodes, mapcoor, mapcoornr, onlyshapenod, blend )

    use limits_m, only: EPSELEMENT, SEARCH_BLOCKS_PASSTWO, EPSELEMENT2

    type(mesh_t), intent(in) :: mesh

!   If present: include element groups for intersection:
!   group(elgrp) = .true.  : include group
!   group(elgrp) = .false. : exclude group
!   Default: include all groups
    logical, dimension(:), intent(in), optional :: groups

!   If present: number of points/nodes that must be computed
!   Default: size(coor,1)
    integer, intent(in), optional :: nnodes

!   coordinates coor(nnodes,ndim) of the points
!   NOTE: if nnodes is present the actual first dimension of coor might be
!         larger than nnodes
    real(dp), dimension(:,:), intent(in) :: coor

!   group and element numbers of the points
!   grpelm(node,1) gives the group number
!   grpelm(node,2) gives the element number
!   if point is outside mesh: grpelm(node,:) = 0
!   NOTE: if nnodes is present the actual first dimension of grpelm might be
!         larger than nnodes
    integer, dimension(:,:), intent(out) :: grpelm

!   reference coordinates refcoor(nnodes,ndim) of the points
!   NOTE: if nnodes is present the actual first dimension of refcoor might be
!         larger than nnodes
    real(dp), dimension(:,:), intent(out) :: refcoor

!   if present and size(onlynodes>0) only elements with all nodes having
!   onlynodes=.true. are considered.
    logical, dimension(:), intent(in), optional :: onlynodes

!   if mapcoor is present the coordinates are mapped to a new position before
!   the intersection with the mesh is computed.
    optional :: mapcoor
    interface
      function mapcoor ( nr, x )
        use kind_defs_m
        implicit none
        integer, intent(in) :: nr
        real(dp), intent(in), dimension(:) :: x
        real(dp), dimension(size(x)) :: mapcoor
      end function mapcoor
    end interface

!   The number mapcoornr can be used to switch between different mappings in
!   the function mapcoor
    integer, intent(in), optional :: mapcoornr

!   optional parameter to indicate that only numshapenod nodes are
!   used in the geometrical shape of the element. This affects only elements
!   with additional internal nodes.
!   default is onlyshapenod = .false. (=all nodes)
    logical, intent(in), optional :: onlyshapenod

!   the mesh chosen for the finding the reference coordinates
!   blend=0  get reference coordinates in the main mesh elements (mesh%element)
!   blend>0  get reference coordinates in the blend mesh elements given by
!            mesh%element_blend(:,blend),
!   default=0
    integer, intent(in), optional :: blend


    integer :: elgrp, elem, node, lmapcoornr, nn, lnnodes
    integer :: ijk(mesh%ndim), p, lblend, GLL, elshape
    real(dp), dimension(size(coor,2)) :: xp, xr, x0
    real(dp) :: epsrefelem
    real(dp), dimension(:), allocatable :: xg
    real(dp), dimension(:,:), allocatable :: x
    logical :: found, lonlyshapenod, chonly, lgroups(mesh%nelgrp)

    lmapcoornr = set_optional ( variable=mapcoornr, default=1 )
    lonlyshapenod = set_optional ( variable=onlyshapenod, default=.false. )
    lnnodes = set_optional ( variable=nnodes, default=size(coor,1) )
    lblend = set_optional ( variable=blend, default=0 )

    if ( lblend > mesh%nblend .or. lblend < 0 ) then
      write(*,'(/2a,i0/)') 'Error find_refcoor_points: ', &
            ' Invalid value for blend = ', lblend
      stop
    end if

    if ( present(onlynodes) ) then
      chonly = size(onlynodes) > 0  ! check for onlynodes to be .true.
    else
      chonly = .false.
    end if

    if ( present(groups) ) then
      lgroups = groups
    else
      lgroups = .true.
    end if

!   compute reference coordinates of points

    do node = 1, lnnodes

      xp = coor(node,:)

!     map coordinates

      if ( present(mapcoor) ) xp = mapcoor ( lmapcoornr, xp )

      found = .false.

!     pass one

      epsrefelem = EPSELEMENT

      call body_search_blocks

      if ( SEARCH_BLOCKS_PASSTWO .and. .not. found ) then

!       pass two

        epsrefelem = EPSELEMENT2

        call body_search_blocks

      end if

      if ( found ) then
!       point is in mesh
        grpelm(node,1) = elgrp
        grpelm(node,2) = elem
        refcoor(node,:) = xr
      else
!       point is outside mesh
        grpelm(node,:) = 0
        refcoor(node,:) = 0 ! set to zero to avoid uninitialized values
      end if

    end do

  contains


!   search for element within blocks or sblocks

    subroutine body_search_blocks

      integer :: elm, block

!     choose blocks subdivision

      if ( mesh%sblocks%filled ) then

!       use structured blocks

        ijk = floor ( ( xp - mesh%sblocks%x0 ) / mesh%sblocks%dx )

        if ( all ( ijk >= 0 ) .and. all ( ijk <= mesh%sblocks%nbl-1 ) ) then

!         point within structured grid

          select case ( mesh%ndim )
          case(1)
            block = 1 + ijk(1)
          case(2)
            block = 1 + ijk(1) + ijk(2) * mesh%sblocks%nbl(1)
          case(3)
            block = 1 + ijk(1) + ijk(2) * mesh%sblocks%nbl(1) + &
                        ijk(3) * mesh%sblocks%nbl(1) * mesh%sblocks%nbl(2)
          case default
            call errormsg_case_default ( 'body_search_blocks', 'mesh%ndim', &
              int_value=mesh%ndim )
          end select

!         loop over all elements in sblock

          do elm = 1, mesh%sblocks%sblks(block)%nelem

            elgrp = mesh%sblocks%sblks(block)%elements(elm,1)
            elem  = mesh%sblocks%sblks(block)%elements(elm,2)

            if ( .not. lgroups(elgrp) ) cycle

            if ( any( xp < mesh%sblocks%sblks(block)%elbounds(elm,:,1) ) .or. &
                 any( xp > mesh%sblocks%sblks(block)%elbounds(elm,:,2) ) ) cycle

            if ( chonly ) then
              if ( .not. all(onlynodes(mesh%topology(elgrp)%a(:,elem))) ) cycle
            end if

            call body_elem_loop_find_refcoor

            if ( found ) exit

          end do

        end if

      else

!       use blocks

   blk: do block = 1, mesh%nblocks

          if ( any( xp < mesh%blocks(block)%bounds(:,1) ) .or. &
               any( xp > mesh%blocks(block)%bounds(:,2) ) ) cycle blk

          do elm = 1, mesh%blocks(block)%nelem

            elgrp = mesh%blocks(block)%elements(elm,1)
            elem  = mesh%blocks(block)%elements(elm,2)

            if ( .not. lgroups(elgrp) ) cycle

            if ( any( xp < mesh%blocks(block)%elbounds(elm,:,1) ) .or. &
                 any( xp > mesh%blocks(block)%elbounds(elm,:,2) ) ) cycle

            if ( chonly ) then
              if ( .not. all(onlynodes(mesh%topology(elgrp)%a(:,elem))) ) cycle
            end if

            call body_elem_loop_find_refcoor

            if ( found ) exit blk

          end do
        end do blk

      end if

    end subroutine body_search_blocks


!   body of the loop over elements in a block to find the reference coordinates

    subroutine body_elem_loop_find_refcoor

!     set some parameters based on lblend

      if ( lblend == 0 ) then

        call set_element_parameters_find_refcoor ( mesh%element(elgrp) )

      else if ( lblend > 0 ) then

        call set_element_parameters_find_refcoor ( &
                                           mesh%element_blend(elgrp,lblend) )

      end if

      allocate ( x(nn,size(xp)) )

!     get coordinates of nodal points

      call get_coordinates ( mesh, elgrp, elem, x, onlyshapenod, lblend )

!     find reference coordinates

      select case ( elshape )

        case(1) ! two-node line

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_line_P1_2 )

        case(2) ! three-node line

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_line_P2_2 )

        case(3) ! three-node triangle

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_triangle_P1 )

        case(4) ! six-node triangle

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_triangle_P2 )

        case(5) ! four-node quadrilateral

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_quad_Q1 )

        case(6) ! nine-node quadrilateral

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_quad_Q2 )

        case(7) ! seven-node triangle

          if ( lonlyshapenod ) then
            call find_refcoor_elem ( xp, xr, nn, x, x0, shape_triangle_P2, &
              onlyshapenod=.true. )
          else
            call find_refcoor_elem ( xp, xr, nn, x, x0, shape_triangle_P2plus )
          end if

        case(9) ! five-node quadrilateral

          if ( lonlyshapenod ) then
            call find_refcoor_elem ( xp, xr, nn, x, x0, shape_quad_Q1, &
              onlyshapenod=.true. )
          else
            call find_refcoor_elem ( xp, xr, nn, x, x0, shape_quad_Q1plus )
          end if

        case(10) ! four-node triangle

          if ( lonlyshapenod ) then
            call find_refcoor_elem ( xp, xr, nn, x, x0, shape_triangle_P1, &
              onlyshapenod=.true. )
          else
            call find_refcoor_elem ( xp, xr, nn, x, x0, shape_triangle_P1plus )
          end if

        case(11) ! four-node tetrahedron

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_tetra_P1 )

        case(12) ! ten-node tetrahedron

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_tetra_P2 )

        case(13) ! eight-node hexahedron

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_hexa_Q1 )

        case(14) ! 27-node hexahedron

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_hexa_Q2 )

        case(15) ! 14-node tetrahedron

          if ( lonlyshapenod ) then
            call find_refcoor_elem ( xp, xr, nn, x, x0, shape_tetra_P2, &
              onlyshapenod=.true. )
          else
            call find_refcoor_elem ( xp, xr, nn, x, x0, shape_tetra_P2plus14 )
          end if

        case(16) ! 15-node tetrahedron

          if ( lonlyshapenod ) then
            call find_refcoor_elem ( xp, xr, nn, x, x0, shape_tetra_P2, &
              onlyshapenod=.true. )
          else
            call find_refcoor_elem ( xp, xr, nn, x, x0, shape_tetra_P2plus15 )
          end if

        case(30) ! 8-node quadrilateral

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_quad_serendipity2 )

        case(31) ! 20-node hexahedron

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_hexa_serendipity2 )

        case(33) ! six-node macro triangle

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_triangle_P1isoP2 )

        case(34) ! nine-node macro quadrilateral

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_quad_Q1isoQ2 )

        case(35) ! ten-node macro tetrahedron

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_tetra_P1isoP2 )

        case(36) ! 27-node macro hexahedron

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_hexa_Q1isoQ2 )

        case(41) ! 6-node prism

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_prism_P1Q1 )

        case(43) ! 18-node prism

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_prism_P2Q2 )

        case(51) ! 5-node pyramid

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_pyramid_Q1P1r )

        case(53) ! 14-node pyramid

          call find_refcoor_elem ( xp, xr, nn, x, x0, shape_pyramid_Q2P2r )

        case(101,102,103,104,105,106,107) ! high-order

          if ( GLL == 1 ) then

!           GLL node distribution

!           compute Gauss Lobatto points
            allocate(xg(p+1))
            call GLL_points ( xg, p )

            select case ( elshape )

              case(101) ! GLL line

                call find_refcoor_elem ( xp, xr, nn, x, x0, &
                  shape_element3=shape_line_GLL_2, p=p, xg=xg )

              case(102) ! GLL quadrilateral

                call find_refcoor_elem ( xp, xr, nn, x, x0, &
                  shape_element3=shape_quad_GLL, p=p, xg=xg )

              case default

                write(*,'(/2a,i0/)') 'Error find_refcoor_points: ', &
                  ' element shape not available for GLL=1: ', elshape
                stop

            end select

            deallocate(xg)

          else

!           high-order equidistant

            select case ( elshape )

              case(101) ! high-order equidistant line

                call find_refcoor_elem ( xp, xr, nn, x, x0, &
                  shape_element2=shape_line_Pp_2, p=p )

              case(102) ! high-order equidistant quadrilateral

                call find_refcoor_elem ( xp, xr, nn, x, x0, &
                  shape_element2=shape_quad_Qp, p=p )

              case(103) ! high-order equidistant hexahedron

                call find_refcoor_elem ( xp, xr, nn, x, x0, &
                  shape_element2=shape_hexa_Qp, p=p )

              case(104) ! high-order equidistant triangle

                call find_refcoor_elem ( xp, xr, nn, x, x0, &
                  shape_element2=shape_triangle_Pp, p=p )

              case(105) ! macro high-order equidistant line

                call find_refcoor_elem ( xp, xr, nn, x, x0, &
                  shape_element4=shape_line_P1isoPp_2, p=p )

              case(106) ! macro high-order equidistant quadrilateral

                call find_refcoor_elem ( xp, xr, nn, x, x0, &
                  shape_element2=shape_quad_Q1isoQp, p=p )

              case(107) ! macro high-order equidistant hexahedron

                call find_refcoor_elem ( xp, xr, nn, x, x0, &
                  shape_element2=shape_hexa_Q1isoQp, p=p )

              case default

                write(*,'(/2a,i0/)') 'Error find_refcoor_points: ', &
                  ' element shape not available for high-order ', elshape
                stop

            end select

          end if

        case default

          write(*,'(/2a,i0/)') 'Error find_refcoor_points: ', &
            ' invalid element shape: ', elshape
          stop

      end select

      deallocate ( x )

!     test whether inside element

      select case ( elshape )

        case(3,4,7,10,11,12,33,35,104)

!         triangle or tetrahedron

          found = all ( xr > -epsrefelem ) .and. sum(xr) < 1._dp+epsrefelem

        case(1,2,5,6,9,13,14,30,31,34,36,101,102,103,105,106,107)

!         line, quadrilateral or hexahedron

          found = all ( abs (xr) < 1._dp+epsrefelem )

        case(41,42,43)

!         prism

          found = all ( xr(1:2) > -epsrefelem ) .and. &
                  sum(xr(1:2)) < 1._dp+epsrefelem .and. &
                  abs(xr(3)) < 1._dp+epsrefelem

        case(51,52,53)

!         pyramid

          found = all ( abs(xr(1:2)) < 1._dp-xr(3)+epsrefelem ) .and. &
                  xr(3) > -epsrefelem .and. xr(3) < 1._dp+epsrefelem

        case default

          call errormsg_case_default ( 'body_elem_loop_find_refcoor', &
            'elshape', int_value=elshape )

      end select

    end subroutine body_elem_loop_find_refcoor


!   set some parameters of the element for find_refcoor

    subroutine set_element_parameters_find_refcoor ( element )

      type(element_t), intent(in) :: element

!       number of nodes used in geometrical shape
        if ( lonlyshapenod ) then
          nn = element%numshapenod
        else
          nn = element%numnod
        end if

!       shape of the element
        elshape = element%elshape

!       order and position of points
        p = element%p(1,1)
        GLL = element%p(1,2)

!       start with initial value for xi to center of element
        x0 =  element%xc

    end subroutine set_element_parameters_find_refcoor

  end subroutine find_refcoor_points


! Find reference coordinates of a point given by it's coordinates in a
! specified element

  subroutine find_refcoor_elem ( xp, xr, nn, x, x0, shape_element, &
    onlyshapenod, shape_element2, shape_element3, shape_element4, p, xg )

    use limits_m, only: MAXITER, EPSITER, STOP_ON_NOCONV, WARN_ON_NOCONV, &
                        MAXVALXI

!   the coordinates of the point
!   xp(dim) are the coordinate of the point in direction dim.
    real(dp), intent(in), dimension(:) :: xp

!   the reference coordinates of the point
!   xr(dim) are the reference coordinates of the point in direction dim.
    real(dp), intent(out), dimension(:) :: xr

!   number of nodes used in the geometrical shape function
    integer, intent(in) :: nn

!   the coordinates of the nodal points of the element
!   x(nn,dim) are the coordinate of the point in direction dim.
    real(dp), intent(in), dimension(:,:) :: x

!   the initial reference coordinates of the point for Newton-Raphson
!   x0(dim) are the reference coordinates of the point in direction dim.
    real(dp), intent(in), dimension(:) :: x0

!   shape function routine of the element.
    optional :: shape_element
    interface
      subroutine shape_element ( xr, phi, dphi )
        use kind_defs_m
        implicit none
        real(dp), intent(in), dimension(:,:) :: xr
        real(dp), intent(out), dimension(:,:) :: phi
        real(dp), intent(out), dimension(:,:,:), optional :: dphi
      end subroutine shape_element
    end interface

!   optional parameter for call to get_coordinates
    logical, intent(in), optional :: onlyshapenod

!   alternative shape function routine of the element.
    optional :: shape_element2
    interface
      subroutine shape_element2 ( xr, p, phi, dphi )
        use kind_defs_m
        implicit none
        real(dp), intent(in), dimension(:,:) :: xr
        integer, intent(in) :: p
        real(dp), intent(out), dimension(:,:) :: phi
        real(dp), intent(out), dimension(:,:,:), optional :: dphi
      end subroutine shape_element2
    end interface

!   alternative shape function routine of the element.
    optional :: shape_element3
    interface
      subroutine shape_element3 ( xr, p, phi, dphi, xg )
        use kind_defs_m
        implicit none
        real(dp), intent(in), dimension(:,:) :: xr
        integer, intent(in) :: p
        real(dp), intent(out), dimension(:,:) :: phi
        real(dp), intent(out), dimension(:,:,:), optional :: dphi
        real(dp), intent(in), dimension(:), optional :: xg
      end subroutine shape_element3
    end interface

!   alternative shape function routine of the element.
    optional :: shape_element4
    interface
      subroutine shape_element4 ( xr, p, phi, dphi, iv )
        use kind_defs_m
        implicit none
        real(dp), intent(in), dimension(:,:) :: xr
        integer, intent(in) :: p
        real(dp), intent(out), dimension(:,:) :: phi
        real(dp), intent(out), dimension(:,:,:), optional :: dphi
        integer, intent(out), dimension(:), optional :: iv
      end subroutine shape_element4
    end interface

!   polynomial order
    integer, intent(in), optional :: p

!   Gauss points
    real(dp), intent(in), dimension(:), optional :: xg


    integer, parameter :: ninti = 1 ! number of `integration points' is just one
    integer :: iter
    real(dp) :: detF(ninti)
!   assume number of degrees == number of nodal points used nn
    real(dp), dimension(ninti,nn) :: phi
    real(dp), dimension(ninti,nn,size(xp)) :: dphi
    real(dp), dimension(ninti,size(xp),size(xp)) :: F, Finv
    real(dp), dimension(ninti,size(xp)) :: xi, xxi
    real(dp), dimension(size(xp)) :: dxi

!   set initial value for xi

    xi(1,:) = x0

    iter = 0

!   start do while loop (Newton-Raphson)

    do

      iter = iter + 1

      if ( present(shape_element) ) then
        call shape_element ( xi, phi, dphi )
      else if ( present(shape_element2) ) then
        call shape_element2 ( xi, p, phi, dphi )
      else if ( present(shape_element3) ) then
        call shape_element3 ( xi, p, phi, dphi, xg )
      else
        call shape_element4 ( xi, p, phi, dphi )
      end if

      call isoparametric_deformation ( x, dphi, F, Finv, detF )

      call isoparametric_coordinates ( x, phi, xxi )

      dxi = -matmul(Finv(1,:,:),xxi(1,:)-xp)

      xi(1,:) = xi(1,:) + dxi

!      print *, 'iter ', iter, 'xi', xi, 'dxi', maxval(abs(dxi))
      if ( maxval(abs(dxi)) < EPSITER ) exit

      if ( iter >= MAXITER .or. any( abs(xi(1,:)) > MAXVALXI ) ) then
        if ( STOP_ON_NOCONV ) then
          write(*,'(/2a,3f20.13/)') 'Error find_refcoor_elem: ', &
            ' no convergence for point ', xp
          stop
        else
          if ( WARN_ON_NOCONV ) then
            write(*,'(/2a,3f20.13/)') &
              'Warning find_refcoor_elem: ', &
              ' no convergence for point ', xp
            write(*,'(a/)') 'Skipping element and continue with searching'
          end if
          xi(1,:) = 2  ! make sure xi is outside the reference element
          exit
        end if
      end if

    end do

    xr = xi(1,:)

  end subroutine find_refcoor_elem

end module meshgen_objects_m
