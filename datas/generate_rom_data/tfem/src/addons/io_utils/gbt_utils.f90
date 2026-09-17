
! Copyright (C) 2007-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines routines for reading a mesh from the gambit neutral format
! Gambit is the mesh generator in the Fluent/Polyflow CFD package.
!
! Notes:
!
!   * if there are elements of different shapes, they are put into separate
!     element groups in the sequence given in the gambit file.
!     For example, if the mesh contains both 6-node triangles (elshape=4)
!     and 9-node quadrilaterals (elshape=6) and they are stored in the file
!     in that order, there will be two tfem element groups: group=1 with
!     elshape=4 and group=2 with elshape=6 (provided there is only a single
!     material).
!     NOTE: if the storage sequence is such that there are 6-nodes triangles,
!     then 9-node quadrilaterals and finally 6-nodes triangles again, there
!     will be three groups!
!
!   * The ELEMENT GROUP (material) blocks can be used to define element groups
!     having the same material properties. The element group numbering
!     within tfem will be: first all elements of the first material
!     then all elements within the second material etc. However if a
!     material/element group within gambit involves elements of different shape,
!     the material group will be split into different element groups according
!     to the shape. For example, if there are two material groups that both
!     contain 6-node triangles and 9-node quads (in that order) there will be
!     four tfem element groups:
!       material 1: group 1 (triangles, elsphape=4)
!                   group 2 (quadrilaterals, elsphape=6)
!       material 2: group 3 (triangles, elsphape=4)
!                   group 4 (quadrilaterals, elsphape=6)
!
!   * BOUNDARY CONDITIONS must be defined on either:
!      1. nodes. These are transfered to nodesets
!      2. elements (edges/faces). These are transfered to:
!         a) curves (edges of triangles or quadrilaterals)
!         b) surfaces (faces of tetrahedra or hexahedra)
!     where each BOUNDARY CONDITIONS block creates a single separate nodeset,
!     curve or surface.
!
!   * The local numbering of the nodes in a curve will be according to the
!     sequence of the element edges as given in the gambit file. Therefore the
!     local numbering of nodes will be in a natural sequence on the curve only
!     if the element edges have a natural consecutive sequence in the gambit
!     file. If not, parameters/options that depend on this natural sequence
!     cannot be used, for example:
!       step, exclude in define_essential, define_constraint,
!                        fill_sysvector, fill_vector
!
!   * The local numbering of the nodes in a surface will be according to the
!     sequence of the global nodal numbers on the surface.
!
!   * No points are generated. If points are needed, for example to set a
!     pressure level, use add_to_mesh to add points to the mesh after reading
!     the mesh.


module gbt_utils_m

  use kind_defs_m
  use mesh_m
  use limits_m, only: MAXPOINTS, MAXCURVES, MAXOBJECTS, MAXSURFACES, &
                      MAXBLOCKS, MAXGROUPS, MAXVOLUMES, MAXNODESETS, &
                      MAXELEMENTSETS
  use meshgen_parts_m, only: fill_element_core
  use misc_m
  use array_defs_m

  implicit none


  integer :: unit_gbt = 13 ! unit number for reading


! interface to support older interfaces

  interface read_gambit_mesh
    module procedure read_mesh_gambit
  end interface read_gambit_mesh

contains


! read mesh from a file written by gambit (neutral file format)

  subroutine read_mesh_gambit ( mesh, filename )

    type(mesh_t), intent(inout) :: mesh

!   the filename for reading the mesh from
    character (len=*), intent(in) :: filename


!   This routine reads basic mesh information (curves, surfaces, topology
!   and coordinates) from a file written by gambit (neutral file format).


    character (len=80) :: current, shifted_left
    character (len=32) :: nameboun

    integer :: numnp, nelem, ngrps, nbsets, ndfcd
    integer :: nr, ne, nelgrp, gelem, elemtype
    integer :: grp_elshape(MAXGROUPS), grp_numel(MAXGROUPS), &
               grp_numnod(MAXGROUPS)
    integer :: grp, bset, ibcode1, nvalues, nentry, itype, facenum, elemnum
    integer :: elshp, gtype, side, celem, selem, elem1, elem2
    integer :: elgrp, elem, curve, ios, i, j, surface, lnode, gnode, nodenr
    integer :: ncurves, nsurfaces
    integer :: ngp, nelgp
    integer :: node, nodeset, nnodesets

    integer, allocatable, dimension(:) :: ntype, ndp, elshapes, &
      gbt_index, elshape_geom, geometry_type, elnumnod, &
      tfem_elemnr, tfem_groupnr, gnodes, lnodes, ordern, gambit_elmnum, &
      typenr
    integer, allocatable, dimension(:,:) :: nodes

    type(int_array_1d_t), allocatable, dimension(:) :: elemnr, groupnr, &
      sidenr, tfem_sidenr, matgrelemnrs, bnodenr

    type(element_t), allocatable, dimension(:) :: element

    real(dp) :: dummy


    if ( mesh%meshgen ) then
      write(*,'(/a/)') 'Error in read_mesh_gambit: mesh already contains data'
      stop
    end if


!   open file

    open ( unit=unit_gbt, file=filename, form='formatted', iostat=ios, &
      status='old' )

    if ( ios /= 0 ) then
      write(*,'(/2a/)') 'Error in read_mesh_gambit: cannot open file ', filename
      stop
    end if


!   first line: CONTROL INFO

    read ( unit=unit_gbt, fmt='(a)' ) current
    shifted_left = adjustl ( current )
    if ( shifted_left(1:12) /= 'CONTROL INFO' ) then
      write(*,'(/3(a/))') 'Error in read_mesh_gambit: expecting CONTROL INFO', &
        current, &
        '     ^  '
      stop
    end if

!   skip five lines in the GAMBIT file

    do i = 1,5
      read ( unit=unit_gbt, fmt=* )
    end do


!   read mesh parameters

    read ( unit=unit_gbt, fmt=* ) numnp, nelem, ngrps, nbsets, ndfcd

    mesh%ndim = ndfcd
    mesh%nnodes = numnp
    mesh%nelem = nelem
    mesh%npoints = 0

!   skip a line

    read ( unit=unit_gbt, fmt=* )


!   COORDINATES

    read ( unit=unit_gbt, fmt='(a)' ) current
    shifted_left = adjustl ( current )
    if ( shifted_left(1:17) /= 'NODAL COORDINATES' ) then
      write(*,'(/3(a/))') &
        'Error in read_mesh_gambit: expecting NODAL COORDINATES', &
        current, &
        '     ^  '
      stop
    end if

    allocate( mesh%coor(mesh%nnodes,mesh%ndim) )

    do nodenr = 1, mesh%nnodes
      read ( unit=unit_gbt, fmt=* ) nr, mesh%coor(nodenr,:)
      if ( nr /= nodenr ) then
        write(*,'(2(/a)/2(a,i0))') 'Error in read_mesh_gambit:', &
          ' Node numbers not in sequence start:1 increment:1. Renumber first', &
          ' nodenr = ', nodenr, ' gambit node = ', nr
        stop
      end if
    end do

!   skip a line

    read ( unit=unit_gbt, fmt=* )


!   ELEMENT TOPOLOGY AND GROUPS

!   read header elements/cells

    read ( unit=unit_gbt, fmt='(a)' ) current
    shifted_left = adjustl ( current )
    if ( shifted_left(1:14) /= 'ELEMENTS/CELLS' ) then
      write(*,'(/3(a/))') &
        'Error in read_mesh_gambit: expecting ELEMENTS/CELLS', &
        current, &
        '     ^  '
      stop
    end if

!   allocate temporary storage for element connectivity

    allocate ( ntype(nelem), ndp(nelem), nodes(nelem,27) )

!   read element data
    do elem = 1, nelem
      read ( unit=unit_gbt, fmt=* ) ne, ntype(elem), ndp(elem), &
        nodes(elem,1:ndp(elem))
    end do

!   skip a line

    read ( unit=unit_gbt, fmt=* )

!   storage for material group elementnrs

    allocate ( matgrelemnrs(ngrps) )

!   read material element groups

    do grp = 1, ngrps

!     read header element groups

      read ( unit=unit_gbt, fmt='(a)' ) current
      shifted_left = adjustl ( current )
      if ( shifted_left(1:13) /= 'ELEMENT GROUP' ) then
        write(*,'(/3(a/))') &
          'Error in read_mesh_gambit: expecting ELEMENT GROUP', &
          current, &
          '     ^  '
        stop
      end if

!     read control info

      read ( unit=unit_gbt, fmt='(7x,i10,12x,i10)') ngp, nelgp

!     skip two lines

      read ( unit=unit_gbt, fmt=* )
      read ( unit=unit_gbt, fmt=* )

      allocate ( matgrelemnrs(grp)%a(nelgp) )

!     read group data
      read ( unit=unit_gbt, fmt='(10I8)' ) matgrelemnrs(grp)%a

!     skip a line

      read ( unit=unit_gbt, fmt=* )

    end do

    allocate ( elshapes(nelem), gambit_elmnum(nelem) )

!   determine elshape of each element

    call set_elshapes

!   translation of new global numbering -> gambit element numbers

    gelem = 0
    do grp = 1, ngrps
      do elem = 1, size(matgrelemnrs(grp)%a)
        gelem = gelem + 1
        gambit_elmnum(gelem) = matgrelemnrs(grp)%a(elem)
      end do
    end do

    if ( gelem /= nelem ) then
      write(*,'(/a/a/)') 'Error in read_mesh_gambit:', &
        ' Material groups not complete (elements missing). '
      stop
    end if

!   set groups

    call set_groups

    mesh%nelgrp = nelgrp

!   element and group info

    allocate( mesh%element(mesh%nelgrp) )
    allocate( mesh%element_blend(mesh%nelgrp,0) )
    mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel(mesh%nelgrp) )

    mesh%element(:)%ndim = mesh%ndim
    mesh%element(:)%elshape = grp_elshape(1:mesh%nelgrp)
    mesh%element(:)%numnod = grp_numnod(1:mesh%nelgrp)
    mesh%grpnumel = grp_numel(1:mesh%nelgrp)

!   set globalshape of elements

    do elgrp = 1, mesh%nelgrp
      select case ( mesh%element(elgrp)%elshape )
        case (1,2); mesh%element(elgrp)%globalshape = 'line'
        case (3,4,7,10); mesh%element(elgrp)%globalshape = 'triangle'
        case (5,6,9,30); mesh%element(elgrp)%globalshape = 'quadrilateral'
        case (13,14,17,31); mesh%element(elgrp)%globalshape = 'hexahedron'
        case (11,12,15,16,18); mesh%element(elgrp)%globalshape = 'tetrahedron'
        case default
          call errormsg_case_default ( 'read_mesh_gambit', &
            'mesh%element(elgrp)%elshape', &
            int_value=mesh%element(elgrp)%elshape )
      end select
    end do

!   set arrays for translation of gambit to tfem numbering

    allocate ( tfem_elemnr(nelem), tfem_groupnr(nelem) )

    gelem = 0
    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)
        gelem = gelem + 1
        tfem_elemnr(gambit_elmnum(gelem)) = elem
        tfem_groupnr(gambit_elmnum(gelem)) = elgrp
      end do
    end do

!   set topology

    allocate( mesh%topology(mesh%nelgrp) )

    gelem = 0

    do elgrp = 1, mesh%nelgrp

      allocate( mesh%topology(elgrp)%a(mesh%element(elgrp)%numnod, &
                            & mesh%grpnumel(elgrp) ) )
      allocate( gbt_index(mesh%element(elgrp)%numnod) )

!     compute permutation arrays for tfem -> gambit numbering of element nodes
      call fill_gbt_index ( mesh%element(elgrp)%elshape, gbt_index )

      do elem = 1, mesh%grpnumel(elgrp)
        gelem = gelem + 1
        mesh%topology(elgrp)%a(:,elem) = &
                             nodes( gambit_elmnum(gelem), gbt_index )
      end do

      deallocate( gbt_index )

    end do

!   test nodes

    allocate ( gnodes(mesh%nnodes) )

    gnodes = 0
    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)
        gnodes(mesh%topology(elgrp)%a(:,elem)) = 1
      end do
    end do

    if ( any( gnodes == 0 ) ) then
      write(*,'(/a/a,i0,a/)') &
        'Warning in read_mesh_gambit: ', &
        ' There are ', count(gnodes==0), &
        ' nodes not connected to the internal elements.'
    end if

    deallocate ( gnodes )


!   points

    allocate( mesh%points(MAXPOINTS) )

    mesh%npoints = 0


!   BOUNDARY SETS

    allocate ( bnodenr(nbsets), elemnr(nbsets), groupnr(nbsets) )
    allocate ( sidenr(nbsets), typenr(nbsets) )
    allocate ( tfem_sidenr(mesh%nelgrp), element(mesh%nelgrp) )

!   sidenodes and sidenode numbering for tfem and gambit
!   use a temporary storage of array of element_t types

    do elgrp = 1, mesh%nelgrp

      element(elgrp) = mesh%element(elgrp)

!     put available info on elements in tfem in element(:)
      call fill_element_core ( element(elgrp) )

!     set translation of gambit faces -> tfem sides
      allocate ( tfem_sidenr(elgrp)%a(element(elgrp)%numsides) )
      call fill_tfem_side ( element(elgrp)%elshape, tfem_sidenr(elgrp)%a )

    end do

!   reading boundary sets

    do bset = 1, nbsets

      read ( unit=unit_gbt, fmt='(a)' ) current
      shifted_left = adjustl ( current )
      if ( shifted_left(1:19) /= 'BOUNDARY CONDITIONS' ) then
        write(*,'(/3(a/))') &
          'Error in read_mesh_gambit: expecting BOUNDARY CONDITIONS', &
          current, &
          '     ^  '
        stop
      end if

!     some data

      read ( unit=unit_gbt, fmt=* ) nameboun, itype, nentry, nvalues, ibcode1

      typenr(bset) = itype

      if ( itype == 0 ) then

!       boundary set defined on nodes

        allocate ( elemnr(bset)%a(0), groupnr(bset)%a(0), sidenr(bset)%a(0) )
        allocate ( bnodenr(bset)%a(nentry) )

        do i = 1, nentry

!         actual boundary data
          read ( unit=unit_gbt, fmt=* ) node, ( dummy, j=1, nvalues )

          bnodenr(bset)%a(i) = node

        end do

      else if ( itype == 1 ) then

!       boundary set defined on element/faces

        allocate ( elemnr(bset)%a(nentry), groupnr(bset)%a(nentry) )
        allocate ( sidenr(bset)%a(nentry), bnodenr(bset)%a(0) )

        do i = 1, nentry

!         actual boundary data
          read ( unit=unit_gbt, fmt=* ) elemnum, elemtype, facenum, &
                                        ( dummy, j=1, nvalues )

!         convert to tfem numbering

          elemnr(bset)%a(i) = tfem_elemnr(elemnum)
          groupnr(bset)%a(i) = tfem_groupnr(elemnum)
          sidenr(bset)%a(i) = tfem_sidenr(groupnr(bset)%a(i))%a(facenum)

        end do

      end if

!       skip a line

      read ( unit=unit_gbt, fmt=* )

    end do

!   select curves and surfaces

    allocate ( elshape_geom(nbsets), geometry_type(nbsets), elnumnod(nbsets) )

    geometry_type = 0  ! initialize to zero

    call set_curves_surfaces

    nnodesets = count ( typenr==0 )
    ncurves   = count ( geometry_type == 1 )
    nsurfaces = count ( geometry_type == 2 )

!   nodesets

    if ( nnodesets > 0 ) then

      mesh%nnodesets = nnodesets

      allocate( mesh%nodesets(max(mesh%nnodesets,MAXNODESETS)) )

      nodeset = 0

      do bset = 1, nbsets

        if ( typenr(bset) /= 0 ) cycle

        nodeset = nodeset + 1

        nentry = size(bnodenr(bset)%a)

        allocate( mesh%nodesets(nodeset)%a(nentry) )

        mesh%nodesets(nodeset)%a = bnodenr(bset)%a

      end do

    else

      mesh%nnodesets = 0

      allocate ( mesh%nodesets(MAXNODESETS) )

    end if

!   curves

    if ( ncurves > 0 ) then

      mesh%ncurves = ncurves

      allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )

!     work space to bookkeep global and local nodes
      allocate ( gnodes(mesh%nnodes), ordern(mesh%nnodes), lnodes(mesh%nnodes) )

      mesh%curves(:mesh%ncurves)%ndim = mesh%ndim

      curve = 0

      do bset = 1, nbsets

        if ( geometry_type(bset) /= 1 ) cycle

        curve = curve + 1

        nentry = size(elemnr(bset)%a)

        mesh%curves(curve)%nelem = nentry

        gnodes = 0

!       fill topology

        allocate( mesh%curves(curve)%topology(elnumnod(bset),&
                 &mesh%curves(curve)%nelem,2) )

        do celem = 1, nentry

          elem = elemnr(bset)%a(celem)
          elgrp = groupnr(bset)%a(celem)
          side = sidenr(bset)%a(celem)

!         set topology of curve
          mesh%curves(curve)%topology(:,celem,2) = &
            mesh%topology(elgrp)%a(element(elgrp)%sidnod(:,side),elem)

!         mark gnode with a unique increasing number for sorting
          elem1 = ( celem - 1 ) * elnumnod(bset) + 1
          elem2 = elem1 + elnumnod(bset) - 1
          gnodes ( mesh%curves(curve)%topology(:,celem,2) ) = &
                              [ ( i, i=elem1,elem2 ) ]

        end do

!       element shapes on curves

        mesh%curves(curve)%element%globalshape = 'line'
        mesh%curves(curve)%element%elshape = elshape_geom(bset)
        mesh%curves(curve)%element%numnod = elnumnod(bset)
        mesh%curves(curve)%element%ndim = mesh%ndim

        mesh%curves(curve)%nblend = 0
        allocate(mesh%curves(curve)%element_blend(0))

        mesh%curves(curve)%nnodes = count ( gnodes /= 0 )

        mesh%curves(curve)%nnodes_blend = [0,mesh%curves(curve)%nnodes]

        allocate ( mesh%curves(curve)%nodes(mesh%curves(curve)%nnodes) )

!       fill local nodes

        call sort ( gnodes, ordern )

        lnode = 0

        do gnode = 1, mesh%nnodes
          if ( gnodes(gnode) == 0 ) cycle  ! global node not in curve
          lnode = lnode + 1 ! new local node
          mesh%curves(curve)%nodes(lnode) = ordern(gnode)
          lnodes(ordern(gnode)) = lnode ! store lnode
        end do

        if ( lnode /= mesh%curves(curve)%nnodes ) then
         print *, 'lnode, nnodes', lnode, mesh%curves(curve)%nnodes
         stop 'read_mesh_gambit: internal error'
        end if

!       fill local topology

        do elem = 1, mesh%curves(curve)%nelem
          mesh%curves(curve)%topology(:,elem,1) = &
             lnodes ( mesh%curves(curve)%topology(:,elem,2) )
        end do

!       set gnodes and lnodes back to zero

        gnodes ( mesh%nnodes-mesh%curves(curve)%nnodes+1: ) = 0
        lnodes ( mesh%curves(curve)%nnodes ) = 0

      end do

      deallocate ( gnodes, ordern, lnodes )

    else

      mesh%ncurves = 0

      allocate( mesh%curves(MAXCURVES) )

    end if

!   surfaces

    if ( nsurfaces > 0 ) then

      mesh%nsurfaces = nsurfaces

      allocate ( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )

      allocate ( gnodes(mesh%nnodes) ) ! work space to bookkeep global nodes

      mesh%surfaces(:mesh%nsurfaces)%ndim = mesh%ndim

      surface = 0

      do bset = 1, nbsets

        if ( geometry_type(bset) /= 2 ) cycle

        surface = surface + 1

        nentry = size(elemnr(bset)%a)

        mesh%surfaces(surface)%nelem = nentry

        gnodes = 0

!       fill topology

        allocate( &
         &mesh%surfaces(surface)%topology(elnumnod(bset),&
                                         &mesh%surfaces(surface)%nelem,2) )

        do selem = 1, nentry
          elem = elemnr(bset)%a(selem)
          elgrp = groupnr(bset)%a(selem)
          side = sidenr(bset)%a(selem)
          mesh%surfaces(surface)%topology(:,selem,2) = &
            mesh%topology(elgrp)%a(element(elgrp)%sidnod(:,side),elem)
          gnodes ( mesh%surfaces(surface)%topology(:,selem,2) ) = 1 ! mark gnode
        end do

!       element shapes on surfaces

        select case ( elshape_geom(bset) )
          case(3,4,7); mesh%surfaces(surface)%element%globalshape = 'triangle'
          case(5,6,30)
                  mesh%surfaces(surface)%element%globalshape = 'quadrilateral'
          case default
            call errormsg_case_default ( 'read_mesh_gambit', &
              'elshape_geom(bset)', int_value=elshape_geom(bset) )
        end select
        mesh%surfaces(surface)%element%elshape = elshape_geom(bset)
        mesh%surfaces(surface)%element%numnod = elnumnod(bset)
        mesh%surfaces(surface)%element%ndim = mesh%ndim

        mesh%surfaces(surface)%nblend = 0
        allocate(mesh%surfaces(surface)%element_blend(0))

        mesh%surfaces(surface)%nnodes = count ( gnodes == 1 )

        mesh%surfaces(surface)%nnodes_blend = [0,mesh%surfaces(surface)%nnodes]

        allocate ( mesh%surfaces(surface)%nodes(mesh%surfaces(surface)%nnodes) )

!       fill local nodes

        lnode = 0

        do gnode = 1, mesh%nnodes
          if ( gnodes(gnode) == 0 ) cycle  ! global node not in surface
          lnode = lnode + 1 ! new local node
          mesh%surfaces(surface)%nodes(lnode) = gnode
          gnodes(gnode) = lnode ! store lnode
        end do

        if ( lnode /= mesh%surfaces(surface)%nnodes ) then
         print *, 'lnode, nnodes', lnode, mesh%surfaces(surface)%nnodes
         stop 'read_mesh_gambit: internal error'
        end if

!       fill local topology

        do elem = 1, mesh%surfaces(surface)%nelem
          mesh%surfaces(surface)%topology(:,elem,1) = &
             gnodes ( mesh%surfaces(surface)%topology(:,elem,2) )
        end do

!       set gnodes back to zero

        gnodes ( mesh%surfaces(surface)%nodes ) = 0

      end do

      deallocate ( gnodes )

    else

      mesh%nsurfaces = 0

      allocate( mesh%surfaces(MAXSURFACES) )

    end if

!   volumes

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate( mesh%objects(MAXOBJECTS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

    mesh%meshgen = .true.

    close ( unit=unit_gbt )

!   deallocate memory

    deallocate ( ntype, ndp, nodes )
    do grp = 1, ngrps
      deallocate ( matgrelemnrs(grp)%a )
    end do
    deallocate ( matgrelemnrs )
    deallocate ( elshapes, gambit_elmnum )
    deallocate ( tfem_elemnr, tfem_groupnr )

    deallocate ( elshape_geom, geometry_type, elnumnod )

    do bset = 1, nbsets
      deallocate ( bnodenr(bset)%a, elemnr(bset)%a, groupnr(bset)%a )
      deallocate ( sidenr(bset)%a )
    end do
    deallocate ( bnodenr, elemnr, groupnr, sidenr )

    deallocate ( typenr )

    do elgrp = 1, mesh%nelgrp
      deallocate ( tfem_sidenr(elgrp)%a )
    end do
    deallocate ( tfem_sidenr )

    do elgrp = 1, mesh%nelgrp
      call deallocate_element_arrays_core ( element(elgrp) )
    end do
    deallocate ( element )

  contains


!   set shape number for each element

    subroutine set_elshapes

      integer :: elem

      do elem = 1, nelem
        select case (ntype(elem))
          case(1) ! Edge
            select case (ndp(elem))
              case(2)
                elshapes(elem) = 1
              case(3)
                elshapes(elem) = 2
              case default
                call errormsg_case_default ( 'set_elshapes', &
                  'ndp(elem)', int_value=ndp(elem) )
            end select
          case(2) ! Quadrilateral
            select case (ndp(elem))
              case(4)
                elshapes(elem) = 5
              case(8)
                elshapes(elem) = 30
              case(9)
                elshapes(elem) = 6
              case default
                call errormsg_case_default ( 'set_elshapes', &
                  'ndp(elem)', int_value=ndp(elem) )
            end select
          case(3) ! Triangle
            select case (ndp(elem))
              case(3)
                elshapes(elem) = 3
              case(6)
                elshapes(elem) = 4
              case(7)
                elshapes(elem) = 7
              case default
                call errormsg_case_default ( 'set_elshapes', &
                  'ndp(elem)', int_value=ndp(elem) )
            end select
          case(4) ! Brick
            select case (ndp(elem))
              case(8)
                elshapes(elem) = 13
              case(20)
                elshapes(elem) = 31
              case(27)
                elshapes(elem) = 14
              case default
                call errormsg_case_default ( 'set_elshapes', &
                  'ndp(elem)', int_value=ndp(elem) )
            end select
          case(5) ! Wedge (Prism)
            write(*,'(/a/)') 'Error in read_mesh_gambit: Wedge not supported'
            stop
          case(6) ! Tetrahedron
            select case (ndp(elem))
              case(4)
                elshapes(elem) = 11
              case(10)
                elshapes(elem) = 12
              case default
                call errormsg_case_default ( 'set_elshapes', &
                  'ndp(elem)', int_value=ndp(elem) )
            end select
          case(7) ! Pyramid
            write(*,'(/a/)') 'Error in read_mesh_gambit: Pyramid not supported'
            stop
          case default
            call errormsg_case_default ( 'set_elshapes', &
              'ntype(elem)', int_value=ntype(elem) )
        end select
      end do

    end subroutine set_elshapes


!   set group info

    subroutine set_groups

      integer :: grp, grelem

      nelgrp = 0

      do grp = 1, ngrps

        nelgrp = nelgrp + 1
        if ( nelgrp > MAXGROUPS ) stop 'Error gbt_utils: nelgrp > MAXGROUPS'
        grp_elshape(nelgrp) = elshapes(matgrelemnrs(grp)%a(1))
        grp_numel(nelgrp) = 1
        grp_numnod(nelgrp) = ndp(matgrelemnrs(grp)%a(1))

        do grelem = 2, size(matgrelemnrs(grp)%a)
          elem = matgrelemnrs(grp)%a(grelem)
          elem1 = matgrelemnrs(grp)%a(grelem-1)
          if ( elshapes(elem) /= elshapes(elem1) ) then
!           change in shape
            nelgrp = nelgrp + 1
            if ( nelgrp > MAXGROUPS ) stop 'Error gbt_utils: nelgrp > MAXGROUPS'
            grp_elshape(nelgrp) = elshapes(elem)
            grp_numel(nelgrp) = 1
            grp_numnod(nelgrp) = ndp(elem)
          else
            grp_numel(nelgrp) = grp_numel(nelgrp) + 1
          end if
        end do

      end do

    end subroutine set_groups


!   set index array for reordering nodes
!   gambit_node = gbt_index ( tfem_node )

    subroutine fill_gbt_index ( elshape, gbt_index )

      integer, intent(in) :: elshape
      integer, dimension(:), intent(out) :: gbt_index

      integer :: i

      select case ( elshape )
        case (12) ! 10 node tetrahedron
          gbt_index = [ 1, 2, 3, 5, 6, 4, 7, 8, 9, 10 ]
        case (13) ! 8 node hexahedron
          gbt_index = [ 1, 2, 4, 3, 5, 6, 8, 7 ]
        case default
          gbt_index = [ ( i, i=1,size(gbt_index) ) ]
      end select

    end subroutine fill_gbt_index


!   set index array for reordering faces
!   tfem_side = tfem_side ( gambit_face )

    subroutine fill_tfem_side ( elshape, tfem_side )

      integer, intent(in) :: elshape
      integer, dimension(:), intent(out) :: tfem_side

      integer :: i

      select case ( elshape )
        case (13,14,17,31) ! hexahedron
          tfem_side = [ 2, 3, 4, 5, 1, 6 ]
        case default
          tfem_side = [ ( i, i=1,size(tfem_side) ) ]
      end select

    end subroutine fill_tfem_side


  ! set curves and surfaces

    subroutine set_curves_surfaces

      integer :: bset, i

      do bset = 1, nbsets

        if ( typenr(bset) == 0 ) cycle ! nodal boundary conditions

        nentry = size(elemnr(bset)%a)

        do i = 1, nentry

!         select curve or surface elements

          elgrp = groupnr(bset)%a(i)

          select case ( mesh%element(elgrp)%elshape )
            case (3,5,9,10)
              elshp = 1  ! 2-node line
              gtype = 1
            case (4,6,7,30)
              elshp = 2  ! 3-node line
              gtype = 1
            case (13,17)
              elshp = 5  ! 4-node quad
              gtype = 2
            case (14)
              elshp = 6  ! 9-node quad
              gtype = 2
            case (31)
              elshp = 30  ! 8-node quad
              gtype = 2
            case (11,18)
              elshp = 3  ! 3-node triangle
              gtype = 2
            case (12)
              elshp = 4  ! 6-node triangle
              gtype = 2
            case (15,16)
              elshp = 7  ! 7-node triangle
              gtype = 2
            case default
              call errormsg_case_default ( 'read_mesh_gambit', &
                'mesh%element(elgrp)%elshape', &
                int_value=mesh%element(elgrp)%elshape )
          end select

          if ( i == 1 ) then
            elshape_geom(bset) = elshp
            geometry_type(bset) = gtype
            elnumnod(bset) = element(elgrp)%sidnumnod
          else
            if ( gtype /= geometry_type(bset) ) then
              write(*,'(/a/a,i0,a/)') 'Error in read_mesh_gambit: ', &
                ' Boundary conditions set ', bset, &
                ' defines both curves and surfaces. This is not supported.'
              stop
            end if
            if ( elshp /= elshape_geom(bset) ) then
              write(*,'(/a/a,i0/2a/)') 'Error in read_mesh_gambit: ', &
                ' Boundary conditions set ', bset, &
                ' defines more than one element shape for curves or surfaces.',&
                ' This is not allowed.'
              stop
            end if
          end if

        end do

      end do

    end subroutine set_curves_surfaces

  end subroutine read_mesh_gambit

end module gbt_utils_m
