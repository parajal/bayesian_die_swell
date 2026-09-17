
! Copyright (C) 2004-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Meshgeneration by construction:
!   - add points, curves, objects, blocks, sblocks, nblocks ...
!   - glue domains for sidelem array

module meshgen_construct_m

  use glob_defs_m
  use mesh_m
  use gauss_m
  use misc_m, only: sort
  use set_optional_m

  implicit none

  logical :: warn_add_to_mesh_after_meshgen_parts = .true.

contains


! Define new entities like curves, surfaces, volumes, objects, blocks,
! nodesets, elementsets

  subroutine add_to_mesh ( mesh, point, points, curve, matchingcurve, &
    curvefromgroups, surface, matchingsurface, surfacefromgroups, object, &
    volumefromgroups, blocks, sblocks, nodblocks, nodeset, nnodes, elementset, &
    group, elements, typeofobject, includeonlygroups, excludegroups, &
    includeonlygroups2, excludegroups2, onlynodeset, onlynodeset2, objectsub, &
    objectcoornr, coor, objectmesh, objectcurve, objectsurface, topology, &
    intrule, nsubint, blocksdomain, step, exclude, excludenodesets, &
    excludepoints, excludecurves, excludesurfaces, xmin, xmax, nodes, &
    replace, mergenodesets, replacegroup, mergeelementsets, elementsetnr, &
    matchingmesh, displacement, rotation_angle, x0, mindistance, &
    warn_mesh_parts, newnr )

    type(mesh_t), intent(inout) :: mesh

!   if present, a point is added to the mesh. The nodal point
!   will be the one that has the smallest distance to the coordinates
!   given by point. Example: point=(/1._dp,0._dp,2,_dp/) will find the nodal
!   point closest to the coordinates (1,0,2) and put the new "point" in that
!   node.
!   NOTE: includeonlygroups and excludegroups can be used
    real(dp), dimension(:), intent(in), optional :: point

!   add multiple points. Dimensions: points(npoints,ndim).
!   NOTE: includeonlygroups and excludegroups can be used
    real(dp), dimension(:,:), intent(in), optional :: points

!   if this array is present a new curve is created consisting of the existing
!   curves specified in the array.
!   For example curve=(/2,-6/) means: create a new curve consisting of
!   curve 2 and curve 6, where the latter is traversed in reverse direction.
!   The reason for the latter possibility is twofold:
!     - the nodal points may need to be in the same order as another curve to
!       use the curve in periodical boundary conditions or in mesh_merge.
!     - the curve may be used in boundary conditions where the normal vector
!       needs to be consistently defined in the same direction on the whole
!       curve.
!   The nodal points of the composite curve are all the unique nodal points
!   in the set of given curves. The connectivity of all the elements
!   in the given curves with the mesh nodal points is preserved. Thus if the
!   individual curves are connected, so will be the composite curve. If the
!   individual curves are disjunct, the composite curve will consist of disjunct
!   parts.
    integer, dimension(:), intent(in), optional :: curve

!   if this array is present a new curve is created consisting of the
!   the nodes of the curve in argument matchingcurve(1) but matching
!   the numbering and topology of matchingcurve(2).
!   NOTE: the coordinates of the nodal points of matchingcurve(2) must be
!   identical to the coordinates of the nodal points of matchingcurve(1)
!   (after optionally adding the virtual displacement displace and/or rotation
!   rotation_angle)
    integer, dimension(:), intent(in), optional :: matchingcurve

!   if this array is present a new surface is created consisting of the existing
!   surfaces specified in the array. For example surface=(/2,6/) means: create
!   a new surface consisting of surface 2 and surface 6.
!   The nodal points of the composite surface are all the unique nodal points
!   in the set of given surfaces and the connectivity of all the elements
!   in the given surfaces with the mesh nodal points is preserved. Thus if the
!   individual surfaces are connected, so will be the composite surface. If the
!   individual surfaces are disjunct, the composite surface will consist of
!   disjunct parts.
    integer, dimension(:), intent(in), optional :: surface

!   if this array is present a new surface is created consisting of the
!   the nodes of the surface in argument matchingsurface(1) but matching
!   the numbering and topology of matchingsurface(2).
!   NOTE: the coordinates of the nodal points of matchingsurface(2) must be
!   identical to the coordinates of the nodal points of matchingsurface(1)
!   (after optionally adding the virtual displacement displace and/or rotation
!   rotation_angle)
    integer, dimension(:), intent(in), optional :: matchingsurface

!   if this array is present a new curve is created consisting of the existing
!   elements in the specified element groups. For example
!       curvefromgroups=(/2,3/)
!   means: create a new curve consisting of all elements in group 2 and 3.
!   The nodal points of the curve are all the unique nodal points
!   in the set of given groups and the connectivity of all the elements
!   in the given groups with the mesh nodal points is preserved. Thus if the
!   individual groups are connected, so will be the curve. If the
!   individual groups are disjunct, the curve will consist of disjunct parts.
!   RESTRICTION: all specified groups must have the same element shape and have
!   a `line' globalshape.
    integer, dimension(:), intent(in), optional :: curvefromgroups

!   if this array is present a new surface is created consisting of the existing
!   elements in the specified element groups. For example
!       surfacefromgroups=(/2,3/)
!   means: create a new surface consisting of all elements in group 2 and 3.
!   The nodal points of the surface are all the unique nodal points
!   in the set of given groups and the connectivity of all the elements
!   in the given groups with the mesh nodal points is preserved. Thus if the
!   individual groups are connected, so will be the surface. If the
!   individual groups are disjunct, the surface will consist of disjunct parts.
!   RESTRICTION: all specified groups must have the same element shape and have
!   a `triangle' or `quadrilateral' globalshape.
    integer, dimension(:), intent(in), optional :: surfacefromgroups

!   if this array is present a new volume is created consisting of the existing
!   elements in the specified element groups. For example
!       volumefromgroups=(/2,3/)
!   means: create a new volume consisting of all elements in group 2 and 3.
!   The nodal points of the volume are all the unique nodal points
!   in the set of given groups and the connectivity of all the elements
!   in the given groups with the mesh nodal points is preserved. Thus if the
!   individual groups are connected, so will be the surface. If the
!   individual groups are disjunct, the surface will consist of disjunct parts.
!   RESTRICTION: all specified groups must have the same element shape and have
!   a `hexahedron', 'tetrahedron', `prism' or 'pyramid' globalshape.
    integer, dimension(:), intent(in), optional :: volumefromgroups

!   if present, an object is added to the mesh.
!   possibilities:
!    'coordinates': object consists of points that need to be defined with the
!                   supplied subroutine objectsub or with the array coor
!    'mesh': object consists of points that are extracted from the supplied
!            objectmesh.
!    'curve': object consists of points that are extracted from the curve
!             objectcurve. If objectmesh is present the curve is taken from the
!             mesh objectmesh, otherwise from mesh.
!    'surface': object consists of points that are extracted from the surface
!               objectsurface. If objectmesh is present the surface is taken
!               from the mesh objectmesh, otherwise from mesh.
!   NOTE: if replace is present the object will replace a previously defined
!         one.
    character(len=*), intent(in), optional :: object

!   if this array is present a new set of blocks is created. The size of the
!   array must be equal to the space dimension of the mesh. For example in
!   2D, specifying blocks=(/10,8/) creates 80 new blocks in a 10x8 equidistant
!   block structure.
!   The elements which have node 1 inside a particular block will be added to
!   that block. Only elements that do not yet belong to previously defined
!   blocks will be added. So an elements can be part of a single block only.
!   If the parameter blocksdomain has not been specified the blocks span the
!   whole current mesh.
    integer, dimension(:), intent(in), optional :: blocks

!   if this array is present a set of structured blocks (sblocks) is created.
!   The size of the array must be equal to the space dimension of the mesh.
!   A previous definition of sblocks will be removed. For example in
!   2D, specifying sblocks=(/10,8/) creates 80 structured blocks in a
!   10x8 equidistant block structure.
!   The elements which are (likely to be) partly of fully inside a particular
!   block will be added to that block. Elements can be part of multiple blocks.
!   If the parameter blocksdomain has not been specified the blocks span the
!   whole current mesh.
    integer, dimension(:), intent(in), optional :: sblocks

!   if this array is present a set of nodal blocks (nodblocks) is created.
!   The size of the array must be equal to the space dimension of the mesh.
!   A previous definition of nodblocks will be removed. For example in
!   2D, specifying nodblocks=(/10,8/) creates 80 nodal blocks in a
!   10x8 equidistant block structure.
!   Around each node a small box is created with the node at the center of the
!   box. The nodes where its box is partly or fully inside a particular
!   block will be added to that block. Nodes can be part of multiple blocks.
!   If the parameter blocksdomain has not been specified the blocks span the
!   whole current mesh.
    integer, dimension(:), intent(in), optional :: nodblocks

!   if present, a nodeset is added to the mesh.
!   possibilities:
!    'nodes': nodeset is specified by the nodes array
!    'merge': nodeset is given by the merging of the nodesets specified in
!             the array mergenodesets.
!   NOTE: if replace is present the nodeset will replace a previously defined
!         one.
    character(len=*), intent(in), optional :: nodeset

!   for nodeset='nodes' the nodes for specifying the nodeset. Nodes that are
!   included more than once will be detected and contribute only once.
!   NOTE: nodes that are out of range (<1 or>mesh%nnodes) are ignored.
    integer, dimension(:), intent(in), optional :: nodes

!   if present, an elementset is added to the mesh.
!   possibilities:
!    'elements': elementset is specified by the group number group
!                and the elements array
!    'merge': elementset is given by the merging of the elementsets specified in
!             the array mergeelementsets.
!    'nodes': add nodes to a previously defined elementset, given by the
!             argument elementsetnr
!   NOTE: if a mesh has multiple groups and elementset='elements' and
!         replacegroup is not present the "other" groups are set to zero
!         elements.
!   NOTE: if replace is present the elementset will replace a previously
!         defined one.
!   NOTE: if replacegroup is present the (group,elements) combo will
!         replace a (group,elements) combo in previously defined elementset.
    character(len=*), intent(in), optional :: elementset

!   if present the group number as used for setting an elementset
!   default=1
    integer, intent(in), optional :: group

!   if present: the elementset will not be added but a single group in
!   an existing elementset is changed.
!   The value of replacegroup gives the elementset number in the mesh.
!   The value of group gives the group number
    integer, intent(in), optional :: replacegroup

!   for elementset='elements' the elements for specifying the elementset.
!   Elements that are included more than once will be detected and
!   contribute only once.
!   NOTE: elements that are out of range (<1 or>mesh%grpnumel(group)) are
!         ignored.
!   NOTE: elements can be specified for a single group only (given by the
!         argument group). If more than one group needs to be specified,
!         the routine needs to called more than once. There are two
!         possibilities here:
!          1) calling add_to_mesh for multiple groups (one at the time)
!             and create a final elementset with elementset='merge'.
!          2) calling add_to_mesh for multiple groups (one at the time).
!             After the first call use replacegroup.
    integer, dimension(:), intent(in), optional :: elements

!   number of points (nodes) in the object.
!   NOTE: if coor is also present and nnodes>size(coor,1), enough memory will
!   be allocate for accommodating nnodes nodes.
    integer, intent(in), optional :: nnodes

!   type of object:
!     1: object is connected to the mesh once.
!     2: object is connected to the mesh twice.
!   default: typeofobject = 1
    integer, intent(in), optional :: typeofobject

!   include only the specified element groups for intersection of the object
!   with the mesh. For example includeonlygroups=(/1,3/), includes only the
!   element groups 1 and 3.
!   Can also be used together with point and points
!   Default is to include all element groups.
    integer, dimension(:), intent(in), optional :: includeonlygroups

!   exclude element groups for intersection of the object with the mesh.
!   For example excludegroups=(/2,4/), excludes element groups 2 and 4.
!   Can also be used together with point and points
!   Default is to include all element groups.
    integer, dimension(:), intent(in), optional :: excludegroups

!   include only the specified element groups for intersection of the object
!   with the mesh. For example includeonlygroups2=(/3,4/), includes only the
!   element groups 3 and 4.
!   Default is to include all element groups.
!   includeonlygroups2 is for the second connection and only makes sense for
!   typeofobject=2.
    integer, dimension(:),intent(in),  optional :: includeonlygroups2

!   exclude element groups for intersection of the object with the mesh.
!   For example excludegroups2=(/2,4/), excludes element groups 2 and 4.
!   Default is to include all element groups.
!   excludegroups2 is for the second connection and only makes sense for
!   typeofobject=2.
    integer, dimension(:), intent(in), optional :: excludegroups2

!   only elements that are fully "embedded" in the specified nodeset are
!   considered for intersection of the object with the mesh.
!   For example nodeset=2, considers only elements that have all nodes
!   in nodeset 2.
!   Default nodeset=0, meaning all elements are considered.
    integer, intent(in), optional :: onlynodeset

!   only elements that are fully "embedded" in the specified nodeset are
!   considered for intersection of the object with the mesh.
!   For example nodeset2=3, considers only elements that have all nodes
!   in nodeset 3.
!   Default nodeset2=0, meaning all elements are considered.
!   nodeset2 is for the second connection and only makes sense for
!   typeofobject=2.
    integer, intent(in), optional :: onlynodeset2

!   user subroutine for defining the coordinates of the points in the object
!   and the (function) number used in the subroutines. NOTE: that objectcoornr
!   is only a parameter supplied to the user subroutine objectsub and it is
!   up to the user what it means.
    optional :: objectsub, objectcoornr
    interface
      subroutine objectsub ( objectcoornr, coor )
        use kind_defs_m
        implicit none
        integer, intent(in) :: objectcoornr
        real(dp), dimension(:,:), intent(inout) :: coor
      end subroutine objectsub
    end interface
    integer, intent(in) :: objectcoornr

!   if present the object is made of the coordinates in coor
!   coor(nnodes,ndim). If the number of nodes in the array coor (=size(coor,1))
!   is smaller than nnodes given in the heading, the latter will be used for
!   memory allocation. In this way the actual number of nodes in the object
!   can be made dynamic without reallocating.
    real(dp), dimension(:,:), intent(in), optional :: coor

!   this mesh defines the object for object='mesh' or (when present) for
!   object='curve' or object='surface'.
!   It must contain a single element group only for object='mesh'
    type(mesh_t), intent(in), optional :: objectmesh

!   this curve defines the object for object='curve'.
!   It uses objectmesh when present, otherwise mesh.
!   step indicates which nodes on the curve are used:
!       step=0 all nodes (default)
!       step>0 nodes 1, 1+step, 1+2*step, ...
!       step<0 except nodes 1, 1-step, 1-2*step, ...
!   exclude indicates which nodes on the curves are excluded:
!       exclude=0 no excludes, all nodes are used (default)
!       exclude=1 exclude first point
!       exclude=2 exclude last point
!       exclude=3 exclude first and last point
!   NOTE: step, exclude only work for curves where the local node numbering is
!   in a natural sequence along the curve. This might not be the case for
!   curves constructed from several other curves or curves read from external
!   mesh generators.
    integer, intent(in), optional :: objectcurve, step, exclude

!   this surface defines the object for object='surface'.
!   It uses objectmesh when present, otherwise mesh.
    integer, intent(in), optional :: objectsurface

!   exclude nodesets from an object made from a mesh, curve or surface:
!   for example excludenodesets=(/1,3/) excludes nodesets 1 and 3.
!   Note that an object that is made by also excluding nodes cannot have a
!   topology.
    integer, intent(in), dimension(:), optional :: excludenodesets

!   exclude points from an object made from a mesh, curve or surface:
!   for example excludepoints=(/1,3/) excludes points 1 and 3.
!   Note that an object that is made by also excluding nodes cannot have a
!   topology.
    integer, intent(in), dimension(:), optional :: excludepoints

!   exclude curves from an object made from a mesh or surface:
!   for example excludecurves=(/1,3/) excludes curves 1 and 3.
!   Note that only the nodes that are actually on the surface are excluded.
!   Note that an object that is made by also excluding nodes cannot have a
!   topology.
    integer, intent(in), dimension(:), optional :: excludecurves

!   exclude surfaces from an object made from a mesh
!   for example excludesurfaces=(/1,3/) excludes surfaces 1 and 3.
!   Note that an object that is made by also excluding nodes cannot have a
!   topology.
    integer, intent(in), dimension(:), optional :: excludesurfaces

!   exclude nodes from an object made from a mesh, curve or surface if
!   any of the coordinates is less than the ones given by xmin or
!   any of the coordinates is larger than the ones given by xmax.
!   For example xmin=(/-1._dp,-2._dp/) excludes all points that have an
!   x-coordinate less than -1 and/or a y-coordinate less than -2.
!   Note that an object that is made by also excluding nodes cannot have a
!   topology.
    real(dp), intent(in), dimension(:), optional :: xmin, xmax

!   add topology for objects (for object='mesh', 'curve', 'surface).
!   default: topology=.false.
!   Note that an object that is made by also excluding nodes cannot have a
!   topology.
    logical, intent(in), optional :: topology

!   if present the blocks, sblocks or nodblocks are defined using the following
!   domain:
!     blocksdomain(:,1) = `lower left corner'
!     blocksdomain(:,2) = `upper right corner'
!   For example:
!      blocksdomains = reshape ( (/ 1, 0, 0, 2, 1, 1 /), (/3,2/) )
!   defines the blocks between the coordinates (1,0,0) and (2,1,1).
    real(dp), dimension(:,:), intent(in), optional :: blocksdomain

!   if present, integration points are added to an object
!   The value of the parameter intrule determines the integration rule.
!   See the type gauss_t for the definition this parameter.
!   NOTE: in order to add integration points objects must have a topology.
    integer, intent(in), optional :: intrule

!   The parameter determines the number of subdomains in the integration
!   defined by intrule.
!   See the type gauss_t for the definition this parameter.
!   Default value: nsubint=1
    integer, intent(in), optional :: nsubint

!   if present: the entity will not be added but an existing one is replaced.
!   This only works for: matchingcurve, matchingsurface, object, nodeset,
!                        elementset
!   The value of replace gives the entity number in the mesh.
    integer, intent(in), optional :: replace

!   if nodeset='merge', this array gives the nodesets to be merged
    integer, dimension(:), intent(in), optional :: mergenodesets

!   if elementset='merge', this array gives the elementsets to be merged
    integer, dimension(:), intent(in), optional :: mergeelementsets

!   if elementset='nodes', this value gives the elementset number
    integer, intent(in), optional :: elementsetnr

!   if present matchingcurve(2)/matchingsurface(2) is defined in
!   matchingmesh, otherwise in mesh.
    type(mesh_t), intent(in), optional :: matchingmesh

!   if present the coordinates of the nodes in
!   matchingcurve(2)/matchingsurface(2) are "displaced"
!   virtually with vector displacement for matching the coordinates of
!   matchingcurve(1)/matchingsurface(1).
!   NOTE: the coordinates of the nodes are not actually changed.
    real(dp), dimension(:), intent(in), optional :: displacement

!   if present the coordinates of the nodes in curve(2)/surface(2) are "rotated"
!   virtually through the angle specified in rotation_angle for matching the
!   coordinates of curve(1)/surface(1).
!   for 2D meshes:
!   rotation_angle(1): angle of rotation in the xy-plane
!   for 3D meshes:
!   rotation_angle(1): angle of rotation around x-axis
!   rotation_angle(2): angle of rotation around y-axis
!   rotation_angle(3): angle of rotation around z-axis
!   NOTE: the coordinates of the nodes are not actually changed.
!   NOTE2: in 3D rotations are applied in an order of increasing dimension, so
!   first around x-axis, then y-axis and finally z-axis
!   NOTE3: if displacement is also present, this is performed after the rotation
    real(dp), dimension(:), intent(in), optional :: rotation_angle

!   if present perform the rotation w.r.t to point x0 (of the same dimension as
!   the mesh), otherwise rotate around origin
    real(dp), dimension(:), intent(in), optional :: x0

!   if present the matching of coordinates in matchingcurve/matchingsurface
!   is based on the minimum distance instead of exact matching.
!   default=.false.
    logical, intent(in), optional :: mindistance

!   This has the same meaning as warn_add_to_mesh_after_meshgen_parts
!   default=.true.
    logical, intent(in), optional :: warn_mesh_parts

!   If present, the value of newnr gives the new entity number in the mesh.
!   This only works for: point, curve, surface
    integer, intent(out), optional :: newnr

!   define new entities like curves from existing entities


    integer :: nnodesl, typel, ncurves, nsurfaces, p, lstep, lexclude
    integer :: nsubintl, ninti, ndim, nelem, ndimxi, i, nn, newobj
    integer :: ne, gnr, elgrp
    integer, allocatable, dimension(:) :: combinednodes, combinedelements
    logical :: topol, lwarn
    logical, allocatable, dimension(:) :: excludenodes
    type(gauss_t) :: gauss


    lwarn = set_optional( variable=warn_mesh_parts, default=.true. )


    if ( .not. mesh%meshgen ) then
      write(*,'(/a/)') &
        'Error add_to_mesh: no mesh present in structure mesh '
      stop
    end if

!   test groups

    if ( present(includeonlygroups) ) then

      if ( any(includeonlygroups <= 0) .or. &
           any(includeonlygroups > mesh%nelgrp) ) then
        write(*,'(/a/a,i0/)') &
          'Error in add_to_mesh: includeonlygroups incorrect, ', &
          'values <= 0 or values > number of groups = ',  mesh%nelgrp
        stop
      end if

    end if

    if ( present(excludegroups) ) then

      if ( any(excludegroups <= 0) .or. &
           any(excludegroups > mesh%nelgrp) ) then
        write(*,'(/a/a,i0/)') &
          'Error in add_to_mesh: excludegroups incorrect, ', &
          'values <= 0 or values > number of groups = ',  mesh%nelgrp
        stop
      end if

    end if

!   check for blended meshes and its support

    call check_blend_add


!   choose for argument

    if ( present(point) ) then

!     add point

      call add_to_mesh_point ( mesh, point, includeonlygroups, &
        excludegroups, newnr )

    else if ( present(points) ) then

!     add points

      do p = 1, size(points,1)
        call add_to_mesh_point ( mesh, points(p,:), includeonlygroups, &
          excludegroups )
      end do

    else if ( present(curve) ) then

!     add curve

      call check_parts ( 'curve' )

      call add_to_mesh_curve ( mesh, curve, newnr )

    else if ( present(matchingcurve) ) then

!     add matching curve

      call check_parts ( 'matchingcurve' )

      call add_to_mesh_curve_match ( mesh, matchingcurve, nr=replace, &
        mesh2=matchingmesh, displace=displacement, rotate=rotation_angle, &
        x0=x0, mindistance=mindistance )

    else if ( present(curvefromgroups) ) then

!     add curve from element groups

      call check_parts ( 'curvefromgroups' )

      call add_to_mesh_curve_from_groups ( mesh, curvefromgroups )

    else if ( present(surface) ) then

!     add surface

      call check_parts ( 'surface' )

      call add_to_mesh_surface ( mesh, surface, newnr )

    else if ( present(surfacefromgroups) ) then

!     add surface from element groups

      call check_parts ( 'surfacefromgroups' )

      call add_to_mesh_surface_from_groups ( mesh, surfacefromgroups )

    else if ( present(matchingsurface) ) then

!     add matching surface

      call check_parts ( 'matchingsurface' )

      call add_to_mesh_surface_match ( mesh, matchingsurface, nr=replace, &
        mesh2=matchingmesh, displace=displacement, rotate=rotation_angle, &
        x0=x0, mindistance=mindistance )

    else if ( present(volumefromgroups) ) then

!     add volume from element groups

      call check_parts ( 'volumefromgroups' )

      call add_to_mesh_volume_from_groups ( mesh, volumefromgroups )

    else if ( present(nodeset) ) then

!     add nodeset

      if ( nodeset == 'nodes' ) then

!       add nodeset using specified nodes

        if ( .not. present(nodes) ) then
          write(*,'(2a/)') &
            'Error in add_to_mesh: nodes is missing from heading ', &
            'for nodeset=''nodes'''
          stop
        end if

        call add_to_mesh_nodeset_from_nodes ( mesh, nodes, nr=replace )

      else if ( nodeset == 'merge' ) then

!       add nodeset by merging existing nodesets

        if ( .not. present(mergenodesets) ) then
          write(*,'(2a/)') &
            'Error in add_to_mesh: mergenodesets is missing from heading ', &
            'for nodeset=''merge'''
          stop
        end if

        if ( any(mergenodesets <= 0) .or. &
             any(mergenodesets > mesh%nnodesets) ) then
          write(*,'(/a/a,i0/)') &
            'Error in add_to_mesh: mergenodesets incorrect, ', &
            'values <= 0 or values > number of nodesets = ',  mesh%nnodesets
          stop
        end if

        nn = sum( [ ( size(mesh%nodesets(mergenodesets(i))%a), &
                                           i=1,size(mergenodesets) ) ] )
        allocate ( combinednodes(nn) )

        combinednodes = [ ( mesh%nodesets(mergenodesets(i))%a, &
                                           i=1,size(mergenodesets) ) ]

        call add_to_mesh_nodeset_from_nodes ( mesh, combinednodes, nr=replace )

        deallocate ( combinednodes )

      else

        write(*,'(2a/)') &
          'Error: unknown nodeset in heading of add_to_mesh: ', nodeset
        stop

      end if

    else if ( present(elementset) ) then

!     add elementset

      if ( elementset == 'elements' ) then

!       add elementset using specified elements

        if ( .not. present(elements) ) then
          write(*,'(2a/)') &
            'Error in add_to_mesh: elements is missing from heading ', &
            'for elementset=''elements'''
          stop
        end if

        call add_to_mesh_elementset_from_elements ( mesh, elements, group, &
          nr=replace, gnr=replacegroup  )

      else if ( elementset == 'merge' ) then

!       add elementset by merging existing elementsets

        if ( .not. present(mergeelementsets) ) then
          write(*,'(2a/)') &
            'Error in add_to_mesh: mergeelementsets is missing from heading ', &
            'for elementset=''merge'''
          stop
        end if

        if ( any(mergeelementsets <= 0) .or. &
             any(mergeelementsets > mesh%nelementsets) ) then
          write(*,'(/a/a,i0/)') &
            'Error in add_to_mesh: mergeelementsets incorrect, ', &
            'values <= 0 or values > number of elementsets = ', &
            mesh%nelementsets
          stop
        end if

        do elgrp = 1, mesh%nelgrp

          ne = sum( [ ( mesh%elementsets(mergeelementsets(i))%grpnumel(elgrp),&
                                             i=1,size(mergeelementsets) ) ] )
          allocate ( combinedelements(ne) )

          combinedelements = &
              [ (  mesh%elementsets(mergeelementsets(i))%elements(elgrp)%a, &
                                             i=1,size(mergeelementsets) ) ]

          if ( elgrp == 1 ) then
!           create new elementset
            call add_to_mesh_elementset_from_elements &
                     ( mesh, combinedelements, nr=replace )
          else if ( elgrp > 1 ) then
!           modify group in existing elementset
            gnr = set_optional ( variable=replace, default=mesh%nelementsets )
            call add_to_mesh_elementset_from_elements &
                     ( mesh, combinedelements, group=elgrp, gnr=gnr )
          end if

          deallocate ( combinedelements )

        end do

      else if ( elementset == 'nodes' ) then

!       add nodes to the elementset

        if ( .not. present(elementsetnr) ) then
          write(*,'(2a/)') &
            'Error in add_to_mesh: elementsetnr is missing from heading ', &
            'for elementset=''nodes'''
          stop
        end if

        call add_to_mesh_add_nodes_to_elementset &
                     ( mesh, elementset=elementsetnr )

      else

        write(*,'(2a/)') &
          'Error: unknown elementset in heading of add_to_mesh: ', elementset
        stop

      end if

    else if ( present(object) ) then

!     add object

      newobj = set_optional ( variable=replace, default=mesh%nobjects+1 )

      if ( mesh%meshparts .and. warn_add_to_mesh_after_meshgen_parts &
           .and. lwarn ) then
        write(*,'(/a/a/a,i0,a/)') &
          'Warning add_to_mesh: adding an object ', &
          ' to a mesh which has been finalized by fill_mesh_parts.', &
          ' Call fill_mesh_parts_objects for object = ', newobj, &
          ' to finalize this object.'
      end if

      if ( present(typeofobject) ) then
        if ( typeofobject < 1 .or. typeofobject > 2 ) then
          write(*,'(/a,i0/)') &
            'Error in add_to_mesh: typeofobject has incorrect value = ', &
            typeofobject
          stop
        end if
        typel = typeofobject
      else
        typel = 1
      end if

      if ( object == 'coordinates' ) then

!       add object using points

        if ( present(objectsub) ) then

!         add object with points using user subroutine

          if ( .not. present(objectcoornr) ) then
            write(*,'(2a/)') &
              'Error in add_to_mesh: objectcoornr is missing from heading ', &
              'for object=''coordinates'' using a user subroutine'
            stop
          end if

          if ( .not. present(nnodes) ) then
            write(*,'(2a/)') &
              'Error in add_to_mesh: nnodes is missing from heading for ', &
              'objectsub subroutine call'
            stop
          end if

          call add_to_mesh_object_points_sub ( mesh, nnodes, typel, objectsub, &
            objectcoornr, nr=replace )

        else if ( present(coor) ) then

!         add object with points using coordinate array

          if ( present(nnodes) ) then
            nnodesl = nnodes ! reserve more nodes than currently in the object
          else
            nnodesl = size(coor,1)
          end if

          call add_to_mesh_object_points_coor ( mesh, nnodesl, typel, coor, &
            nr=replace )

        else

          write(*,'(2a/)') &
            'Error: for object=''coordinates'' in heading of add_to_mesh', &
            '  either objectsub or coor must be present'
          stop

        end if

      else if ( object == 'mesh' ) then

        if ( .not. present(objectmesh) ) then
          write(*,'(2a/)') &
            'Error in add_to_mesh: objectmesh is missing from heading for ', &
            'object=''mesh'''
          stop
        end if

        topol = set_optional ( variable=topology, default=.false. )

        allocate ( excludenodes(objectmesh%nnodes) )

        excludenodes = .false.
        call set_excludenodes ( objectmesh )

        call add_to_mesh_object_mesh ( mesh, typel, objectmesh, topol, &
          excludenodes, nr=replace  )

        deallocate ( excludenodes )

      else if ( object == 'curve' ) then

        if ( .not. present(objectcurve) ) then
          write(*,'(2a/)') &
            'Error in add_to_mesh: objectcurve is missing from heading for ', &
            'object=''curve'''
          stop
        end if

        if ( present(objectmesh) ) then
!         use objectmesh
          ncurves = objectmesh%ncurves
        else
!         use mesh
          ncurves = mesh%ncurves
        end if

        if ( objectcurve <= 0 .or. objectcurve > ncurves ) then
          write(*,'(/a/a,i0/)') &
            'Error in add_to_mesh: objectcurve incorrect, ', &
            'values <= 0 or values > number of curves = ',  ncurves
          stop
        end if

        topol = set_optional ( variable=topology, default=.false. )
        lstep = set_optional ( variable=step, default=0 )
        lexclude = set_optional ( variable=exclude, default=0 )

        if ( present(objectmesh) ) then

!         use objectmesh

          allocate ( excludenodes(objectmesh%nnodes) )

          call set_excludenodes_curve ( objectmesh%curves(objectcurve) )
          call set_excludenodes ( objectmesh )

          call add_to_mesh_object_geometry ( mesh, typel, &
            objectgeometry=objectmesh%curves(objectcurve), &
            coor=objectmesh%coor, topol=topol, excludenodes=excludenodes, &
            nr=replace  )

          deallocate ( excludenodes )

        else

!         use mesh

          allocate ( excludenodes(mesh%nnodes) )

          call set_excludenodes_curve ( mesh%curves(objectcurve) )
          call set_excludenodes ( mesh )

          call add_to_mesh_object_geometry ( mesh, typel, &
            objectgeometry=mesh%curves(objectcurve), &
            coor=mesh%coor, topol=topol, excludenodes=excludenodes, &
            nr=replace  )

          deallocate ( excludenodes )

        end if

      else if ( object == 'surface' ) then

        if ( .not. present(objectsurface) ) then
          write(*,'(2a/)') &
            'Error in add_to_mesh: objectsurface is missing from heading for', &
            ' object=''surface'''
          stop
        end if

        if ( present(objectmesh) ) then
!         use objectmesh
          nsurfaces = objectmesh%nsurfaces
        else
!         use mesh
          nsurfaces = mesh%nsurfaces
        end if

        if ( objectsurface <= 0 .or. objectsurface > nsurfaces ) then
          write(*,'(/a/a,i0/)') &
            'Error in add_to_mesh: objectsurface incorrect, ', &
            'values <= 0 or values > number of surfaces = ',  nsurfaces
          stop
        end if

        topol = set_optional ( variable=topology, default=.false. )

        if ( present(objectmesh) ) then

!         use objectmesh

          allocate ( excludenodes(objectmesh%nnodes) )

          excludenodes = .false.
          call set_excludenodes ( objectmesh )

          call add_to_mesh_object_geometry ( mesh, typel, &
            objectgeometry=objectmesh%surfaces(objectsurface), &
            coor=objectmesh%coor, topol=topol, excludenodes=excludenodes, &
            nr=replace  )

          deallocate ( excludenodes )

        else

!         use mesh

          allocate ( excludenodes(mesh%nnodes) )

          excludenodes = .false.
          call set_excludenodes ( mesh )

          call add_to_mesh_object_geometry ( mesh, typel, &
            objectgeometry=mesh%surfaces(objectsurface), &
            coor=mesh%coor, topol=topol, excludenodes=excludenodes, &
            nr=replace  )

          deallocate ( excludenodes )

        end if

      else

        write(*,'(2a/)') &
          'Error: unknown object in heading of add_to_mesh: ', object
        stop

      end if

!     set nodeset

      if ( present(onlynodeset) ) then
        if ( onlynodeset < 0 .or. onlynodeset > mesh%nnodesets ) then
          write(*,'(/a/a,i0/)') &
            'Error in add_to_mesh: onlynodeset has incorrect value = ', &
            'value < 0 or value > number of nodesets = ', mesh%nnodesets
          stop
        end if
        mesh%objects(newobj)%nodeset = onlynodeset
      end if

      if ( typel == 2 ) then

!       second intersection (typeofobject==2)

        if ( present(onlynodeset2) ) then
          if ( onlynodeset2 < 0 .or. onlynodeset2 > mesh%nnodesets ) then
            write(*,'(/a/a,i0/)') &
              'Error in add_to_mesh: onlynodeset2 has incorrect value = ', &
              'value < 0 or value > number of nodesets = ', mesh%nnodesets
            stop
          end if
          mesh%objects(newobj)%nodeset2 = onlynodeset2
        end if

      end if

      allocate ( mesh%objects(newobj)%groups(mesh%nelgrp) )

!     which groups?

      if ( present(includeonlygroups) ) then

        mesh%objects(newobj)%groups = .false.
        mesh%objects(newobj)%groups(includeonlygroups) = .true.

      else

!       default: all groups
        mesh%objects(newobj)%groups = .true.

      end if

      if ( present(excludegroups) ) then

        mesh%objects(newobj)%groups(excludegroups) = .false.

      end if

      if ( typel == 2 ) then

!       second intersection (typeofobject==2)

        allocate ( mesh%objects(newobj)%groups2(mesh%nelgrp) )

!       which groups?

        if ( present(includeonlygroups2) ) then

          if ( any(includeonlygroups2 <= 0) .or. &
               any(includeonlygroups2 > mesh%nelgrp) ) then
            write(*,'(/a/a,i0/)') &
              'Error in add_to_mesh: includeonlygroups2 incorrect, ', &
              'values <= 0 or values > number of groups = ',  mesh%nelgrp
            stop
          end if

          mesh%objects(newobj)%groups2 = .false.
          mesh%objects(newobj)%groups2(includeonlygroups2) = .true.

        else

!         default: all groups
          mesh%objects(newobj)%groups2 = .true.

        end if

        if ( present(excludegroups2) ) then

          if ( any(excludegroups2 <= 0) .or. &
               any(excludegroups2 > mesh%nelgrp) ) then
            write(*,'(/a/a,i0/)') &
              'Error in add_to_mesh: excludegroups2 incorrect, ', &
              'values <= 0 or values > number of groups = ', mesh%nelgrp
            stop
          end if

          mesh%objects(newobj)%groups2(excludegroups2) = .false.

        end if

      end if

!     integration points for objects

      if ( present(intrule) ) then

!       add integration points

        if ( .not. mesh%objects(newobj)%topol ) then
          write(*,'(2(a/),a,i0/)') &
            'Error in add_to_mesh: integration points for objects ', &
            '  are only possible for objects having a topology', &
            '  objectnr = ', newobj
          stop
        end if

        nsubintl = set_optional ( variable=nsubint, default=1 )

        gauss%globalshape = mesh%objects(newobj)%element%globalshape
        gauss%intrule = intrule
        gauss%nsubint = nsubintl

!       set number of integration points

        call set_ninti ( gauss, ninti )

        mesh%objects(newobj)%intpoints = .true.
        mesh%objects(newobj)%ninti = ninti

        select case ( gauss%globalshape )
        case ( 'line' )
          ndimxi = 1
        case ( 'triangle', 'quadrilateral' )
          ndimxi = 2
        case default
          ndimxi = 3
        end select
        ndim = mesh%ndim
        nelem = mesh%objects(newobj)%nelem

!       allocate data

        allocate ( mesh%objects(newobj)%xig(ninti,ndimxi), &
          mesh%objects(newobj)%wg(ninti), &
          mesh%objects(newobj)%coor_int(ninti,ndim,nelem), &
          mesh%objects(newobj)%refcoor_int(ninti,ndim,nelem), &
          mesh%objects(newobj)%grpelm_int(ninti,2,nelem) )

        if ( typel == 2 ) then

!         second intersection (typeofobject==2)

          allocate ( &
            mesh%objects(newobj)%refcoor2_int(ninti,ndim,nelem), &
            mesh%objects(newobj)%grpelm2_int(ninti,2,nelem) )

        end if

!       compute integration points and weights

        call set_Gauss_integration ( gauss, mesh%objects(newobj)%xig, &
           mesh%objects(newobj)%wg )

      end if

    else if ( present(blocks) ) then

!     add blocks

      if ( size(blocks) /= mesh%ndim ) then
        write(*,'(/2a,i0/)') &
          'Error in add_to_mesh: size of array blocks must be equal to the', &
          ' space dimension = ', mesh%ndim
        stop
      end if

      if ( any( blocks < 1 ) ) then
        write(*,'(/2a/)') &
          'Error in add_to_mesh: values in array blocks must be larger than', &
          ' zero'
        stop
      end if

      if ( present(blocksdomain) ) then

        if ( size(blocksdomain,1) /= mesh%ndim ) then
          write(*,'(/2a,i0/)') &
            'Error in add_to_mesh: size of the first dimension of', &
            ' blocksdomain must be equal to the space dimension = ', mesh%ndim
          stop
        end if
        if ( size(blocksdomain,2) /= 2 ) then
          write(*,'(/2a/)') &
            'Error in add_to_mesh: size of the second dimension of', &
            ' blocksdomain must be 2 '
          stop
        end if

        call add_to_mesh_blocks ( mesh, blocks, blocksdomain )

      else

        call add_to_mesh_blocks ( mesh, blocks )

      end if

    else if ( present(sblocks) ) then

!     add structured blocks

      if ( size(sblocks) /= mesh%ndim ) then
        write(*,'(/2a,i0/)') &
          'Error in add_to_mesh: size of array sblocks must be equal to the', &
          ' space dimension = ', mesh%ndim
        stop
      end if

      if ( any( sblocks < 1 ) ) then
        write(*,'(/2a/)') &
          'Error in add_to_mesh: values in array blocks must be larger than', &
          ' zero'
        stop
      end if

      if ( present(blocksdomain) ) then

        if ( size(blocksdomain,1) /= mesh%ndim ) then
          write(*,'(/2a,i0/)') &
            'Error in add_to_mesh: size of the first dimension of', &
            ' blocksdomain must be equal to the space dimension = ', mesh%ndim
          stop
        end if
        if ( size(blocksdomain,2) /= 2 ) then
          write(*,'(/2a/)') &
            'Error in add_to_mesh: size of the second dimension of', &
            ' blocksdomain must be 2 '
          stop
        end if

        call add_to_mesh_sblocks ( mesh, sblocks, blocksdomain )

      else

        call add_to_mesh_sblocks ( mesh, sblocks )

      end if

    else if ( present(nodblocks) ) then

!     add nodal blocks

      if ( size(nodblocks) /= mesh%ndim ) then
        write(*,'(/2a,i0/)') &
          'Error in add_to_mesh: size of array nodblocks must be equal to the',&
          ' space dimension = ', mesh%ndim
        stop
      end if

      if ( any( nodblocks < 1 ) ) then
        write(*,'(/2a/)') &
          'Error in add_to_mesh: values in array blocks must be larger than', &
          ' zero'
        stop
      end if

      if ( present(blocksdomain) ) then

        if ( size(blocksdomain,1) /= mesh%ndim ) then
          write(*,'(/2a,i0/)') &
            'Error in add_to_mesh: size of the first dimension of', &
            ' blocksdomain must be equal to the space dimension = ', mesh%ndim
          stop
        end if
        if ( size(blocksdomain,2) /= 2 ) then
          write(*,'(/2a/)') &
            'Error in add_to_mesh: size of the second dimension of', &
            ' blocksdomain must be 2 '
          stop
        end if

        call add_to_mesh_nodblocks ( mesh, nodblocks, blocksdomain )

      else

        call add_to_mesh_nodblocks ( mesh, nodblocks )

      end if

    else

      write(*,'(20(/a)/)') &
        'Error: could not determine type of addition to mesh', &
        'from the heading of add_to_mesh. Non of the keywords:',      &
        '  point ', &
        '  points ', &
        '  curve ', &
        '  matchingcurve ', &
        '  curvefromgroups ', &
        '  surface ', &
        '  matchingsurface ', &
        '  surfacefromgroups ', &
        '  volumefromgroups ', &
        '  nodeset ', &
        '  elementset ', &
        '  object ', &
        '  blocks ', &
        '  sblocks ', &
        '  nodblocks ', &
        'is present.'
      stop

    end if

  contains


!   check for support of blended meshes

    subroutine check_blend_add

      if ( present(point) ) then
        call check_blend ( mesh, 'add_to_mesh', keyword='point', &
          comment='Blended meshes commonly have double nodes. Fix needed.' )
      else if ( present(points) ) then
        call check_blend ( mesh, 'add_to_mesh', keyword='points', &
          comment='Blended meshes commonly have double nodes. Fix needed.' )
      else if ( present(curve) ) then
        call check_blend ( mesh, 'add_to_mesh', keyword='curve', &
          comment='Add curve to each mesh before creating blended mesh.' )
      else if ( present(curvefromgroups) ) then
        call check_blend ( mesh, 'add_to_mesh', keyword='curvefromgroups', &
          comment='Add curve to each mesh before creating blended mesh.' )
      else if ( present(surface) ) then
        call check_blend ( mesh, 'add_to_mesh', keyword='surface', comment= &
           'Add surface to each mesh before creating blended mesh.' )
      else if ( present(surfacefromgroups) ) then
        call check_blend ( mesh, 'add_to_mesh', keyword='surfacefromgroups', &
          comment= 'Add surface to each mesh before creating blended mesh.' )
      else if ( present(volumefromgroups) ) then
        call check_blend ( mesh, 'add_to_mesh', keyword='volumefromgroups', &
          comment= 'Add volume to each mesh before creating blended mesh.' )
      else if ( present(elementset) ) then
        if ( elementset == 'nodes' ) then
          call check_blend ( mesh, 'add_to_mesh', warning=.true., &
            keyword='elemenset=''nodes''', &
            comment='Only nodes from the main element are added. Fix needed.' )
        end if
      else if ( present(object) ) then
        if ( mesh%nblend > 0 .and. object /= 'coordinates' .and. &
             .not. present(objectmesh) ) then
          write(*,'(/a/a/a,i0,a/)') &
            'Error add_to_mesh: adding an object ', &
            ' to a blended mesh is only supported for object=''coordinates''', &
            ' or using an objectmesh not consisting of blended meshes '
          stop
          if ( present(objectmesh) ) then
            call check_blend ( objectmesh, 'add_to_mesh', keyword='object', &
            comment='Argument objectmesh should not consist of blended meshes.')
          end if
        end if
      end if

    end subroutine check_blend_add


!   check whether fill_mesh_parts has been called

    subroutine check_parts ( string )

      character(len=*), intent(in) :: string

      if ( mesh%meshparts ) then
        write(*,'(/2a/a/)') &
          'Error add_to_mesh: cannot add ', string, &
          ' to a mesh which has been finalized by fill_mesh_parts'
        stop
      end if

    end subroutine check_parts


!   set excludenodes for curve

    subroutine set_excludenodes_curve ( curve )

      type(geometry_t), intent(in) :: curve

      if ( lstep > 0 ) then
        excludenodes = .false.
        excludenodes ( curve%nodes(1::lstep) ) = .true.
      else if ( lstep < 0 ) then
        excludenodes = .true.
        excludenodes ( curve%nodes(1::-lstep) ) = .false.
      end if

      excludenodes = .false.
      if ( lexclude == 1 .or. lexclude == 3 ) then
        excludenodes(curve%nodes(1)) = .true.
      end if
      if ( lexclude == 2 .or. lexclude == 3 ) then
        excludenodes(curve%nodes(curve%nnodes)) = .true.
      end if

    end subroutine set_excludenodes_curve


!   set excludenodes from points, curves and surfaces

    subroutine set_excludenodes ( mesh )

      type(mesh_t), intent(in) :: mesh

      integer :: nst, crv, srf, node

      if ( present(excludenodesets) ) then
        if ( any(excludenodesets <=0) .or. &
             any(excludenodesets > mesh%nnodesets) ) then
          write(*,'(2(/a),i0/)') &
            'Error in add_to_mesh: ',&
            ' excludenodesets <=0 or excludenodesets > number of nodesets = ', &
             mesh%nnodesets
          stop
        end if
        do nst = 1, size(excludenodesets)
          excludenodes ( mesh%nodesets(excludenodesets(nst))%a ) = .true.
        end do
      end if

      if ( present(excludepoints) ) then
        if ( any(excludepoints <=0) .or. &
             any(excludepoints > mesh%npoints) ) then
          write(*,'(2(/a),i0/)') &
            'Error in add_to_mesh: ',&
            ' excludepoints <=0 or excludepoints > number of points = ', &
             mesh%npoints
          stop
        end if
        excludenodes ( mesh%points(excludepoints) ) = .true.
      end if

      if ( present(excludecurves) ) then
        if ( any(excludecurves <=0) .or. &
             any(excludecurves > mesh%ncurves) ) then
          write(*,'(2(/a),i0/)') &
            'Error in add_to_mesh: ',&
            ' excludecurves <=0 or excludecurves > number of curves = ', &
             mesh%ncurves
          stop
        end if
        do crv = 1, size(excludecurves)
          excludenodes ( mesh%curves(excludecurves(crv))%nodes ) = .true.
        end do
      end if

      if ( present(excludesurfaces) ) then
        if ( any(excludesurfaces <=0) .or. &
             any(excludesurfaces > mesh%nsurfaces) ) then
          write(*,'(2(/a),i0/)') &
            'Error in add_to_mesh: ',&
            ' excludesurfaces <=0 or excludesurfaces > number of surfaces = ', &
             mesh%nsurfaces
          stop
        end if
        do srf = 1, size(excludesurfaces)
          excludenodes ( mesh%surfaces(excludesurfaces(srf))%nodes ) = .true.
        end do
      end if

      if ( present(xmin) ) then
        if ( size(xmin) /= mesh%ndim ) then
          write(*,'(2(/a),i0/)') &
            'Error in add_to_mesh: ',&
            ' dimension of xmin must be equal to the space dimension = ', &
             mesh%ndim
          stop
        end if
        do node = 1, mesh%nnodes
          if ( all ( mesh%coor(node,:) > xmin ) ) cycle  ! in region
          excludenodes ( node ) = .true.
        end do
      end if

      if ( present(xmax) ) then
        if ( size(xmax) /= mesh%ndim ) then
          write(*,'(2(/a),i0/)') &
            'Error in add_to_mesh: ',&
            ' dimension of xmax must be equal to the space dimension = ', &
             mesh%ndim
          stop
        end if
        do node = 1, mesh%nnodes
          if ( all ( mesh%coor(node,:) < xmax ) ) cycle  ! in region
          excludenodes ( node ) = .true.
        end do
      end if

    end subroutine set_excludenodes


  end subroutine add_to_mesh


! Define new point from coordinates given.

  subroutine add_to_mesh_point ( mesh, point, includeonlygroups, &
    excludegroups, newnr )

    type(mesh_t), intent(inout) :: mesh

!   the coordinates (1:ndim) of the point
    real(dp), dimension(:), intent(in) :: point

!   include only the specified element groups
!   Default is to include all element groups.
    integer, dimension(:), intent(in), optional :: includeonlygroups

!   exclude element groups
!   Default is to include all element groups.
    integer, dimension(:), intent(in), optional :: excludegroups

!   If present, the value of newnr gives the new point number in the mesh.
    integer, intent(out), optional :: newnr


    logical :: groups(mesh%nelgrp)
    logical, allocatable, dimension(:) :: mask
    real(dp), parameter :: eps = 1e-10_dp
    integer :: newpoint, node, nodenr, elgrp, elem
    real(dp), allocatable, dimension(:) :: distance2

    if ( size(point) /= mesh%ndim ) then
      write(*,'(/a/a,i0/)') &
        'Error: the space dimension of point/points in add_to_mesh:', &
        ' must be the dimension of the space = ', mesh%ndim
      stop
    end if

    newpoint = mesh%npoints + 1

    if ( newpoint > size(mesh%points) ) then
      write(*,'(/a/a,i0/)') &
        'Error: maximum exceeded in add_to_mesh:', &
        ' pointnumber > maximum = ', size(mesh%points)
      stop
    end if

!   groups

    if ( present(includeonlygroups) ) then
      groups = .false.
      groups(includeonlygroups) = .true.
    else
!     default: all groups
      groups = .true.
    end if

    if ( present(excludegroups) ) then
      groups(excludegroups) = .false.
    end if

!   compute square distance

    allocate ( distance2(mesh%nnodes) )

    do node = 1, mesh%nnodes
      distance2(node) = &
         dot_product ( mesh%coor(node,:) - point, mesh%coor(node,:) - point )
    end do

!   check groups

    allocate ( mask(mesh%nnodes) )

    if ( all(groups) ) then
      nodenr = minval ( minloc ( distance2 ) )
    else
      mask = .false.
      do elgrp = 1, mesh%nelgrp
        if ( .not. groups(elgrp) ) cycle
        do elem = 1, mesh%grpnumel(elgrp)
          mask(mesh%topology(elgrp)%a(:,elem)) = .true.
        end do
      end do
      nodenr = minval ( minloc ( distance2, mask=mask ) )
    end if

    mesh%points(newpoint) = nodenr
    mesh%npoints = mesh%npoints + 1

    if ( present(newnr) ) newnr = mesh%npoints

    if ( distance2(nodenr) >= eps * maxval(distance2) ) then
      write(*,'(/2a/2(a,3/))') &
        'Warning: added point in add_to_mesh does not coincide with a nodal', &
        ' point.'
      write(*,'(a,3es26.14)') 'Requested coordinates = ', point
      write(*,'(a,3es26.14)') 'Actual coordinates    = ', mesh%coor(nodenr,:)
    end if

    deallocate ( mask, distance2 )

  end subroutine add_to_mesh_point


! Define new curve based on existing ones

  subroutine add_to_mesh_curve ( mesh, curve, newnr )

    type(mesh_t), intent(inout) :: mesh

!   a new curve is created consisting of the existing curves specified in the
!   array. For example curve=(/2,-6/) means: create a new curve consisting of
!   curve 2 and curve 6, where the latter is traversed in reverse direction.
!   The reason for the latter possibility is twofold:
!     - the nodal points may need to be in the same order as another curve to
!       use the curve in periodical boundary conditions or in mesh_merge.
!     - the curve may be used in boundary conditions where the normal vector
!       needs to be consistently defined in the same direction on the whole
!       curve.
!   The nodal points of the composite curve are all the unique nodal points
!   in the set of given curves and the connectivity of all the elements
!   in the given curves with the mesh nodal points is preserved. Thus if the
!   individual curves are connected, so will be the composite curve. If the
!   individual curves are disjunct, the composite curve will consist of disjunct
!   parts.
    integer, dimension(:), intent(in) :: curve

!   If present, the value of newnr gives the new curve number in the mesh.
    integer, intent(out), optional :: newnr


!   define new curves from existing entities

    integer :: crv, nelem, nnodes, elnumnod, newcurve, elem, cr, node
    integer :: lnodenr, newnnodes, nodenr
!   temporaries for nodes
    integer, dimension(:), allocatable :: work1, work3, work4
    integer, dimension(:,:), allocatable :: work2   ! temporary for topology

!   add curves

    if ( size( curve ) == 0 ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  curve array is empty'
      stop
    end if

    if ( any( curve == 0 ) .or. any( abs(curve) > mesh%ncurves ) ) then
      write(*,'(/a/a)',advance='no') &
        'Error in add_to_mesh: ', &
        '  curve array contains non-existing curves: '
      write(*,*) curve
      write(*,*)
      stop
    end if

    if ( any( mesh%curves(abs(curve))%element%numnod /= &
                   mesh%curves(abs(curve(1)))%element%numnod ) ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  curves in curve array have different element shapes'
      stop
    end if

    newcurve = mesh%ncurves + 1

    if ( newcurve > size(mesh%curves) ) then
      write(*,'(/a/a,i0/)') &
        'Error: maximum exceeded in add_to_mesh:', &
        ' curves > maximum = ', size(mesh%curves)
      stop
    end if

!   temporary for storing all the nodes in the curves (including double ones)
    allocate ( work1( sum(mesh%curves(abs(curve))%nnodes) ) )
!   temporary for storing the topology of the combined curves
    allocate ( work2( mesh%curves(1)%element%numnod, &
                    sum(mesh%curves(abs(curve))%nelem)) )

    nnodes = 0
    nelem = 0

    elnumnod = mesh%curves(abs(curve(1)))%element%numnod

!   loop over the curves

    do cr = 1, size(curve)

      crv = abs(curve(cr))

!     store all nodes and new topology

      if ( curve(cr) > 0 ) then

!       normal direction

        work1(nnodes+1:nnodes+mesh%curves(crv)%nnodes) = &
                                           mesh%curves(crv)%nodes
        work2(:,nelem+1:nelem+mesh%curves(crv)%nelem) = &
                              mesh%curves(crv)%topology(:,:,1) + nnodes

      else if ( curve(cr) < 0 ) then

!       reverse direction

!       reverse global nodes
        work1(nnodes+1:nnodes+mesh%curves(crv)%nnodes) = &
          mesh%curves(crv)%nodes(mesh%curves(crv)%nnodes:1:-1)
!       reverse topology
        work2(:,nelem+1:nelem+mesh%curves(crv)%nelem) = &
           mesh%curves(crv)%nnodes + 1  &
             - mesh%curves(crv)%topology(elnumnod:1:-1, &
                              mesh%curves(crv)%nelem:1:-1,1) + nnodes

      end if

      nnodes = nnodes + mesh%curves(crv)%nnodes
      nelem = nelem + mesh%curves(crv)%nelem

    end do

!   check for double nodes

!   work3: store at each (global) node in the curves the new local node number,
!          otherwise store a zero
!   work4: store at each local node the corresponding mesh node number

    allocate ( work3(mesh%nnodes), work4(nnodes) )

    work3 = 0

    lnodenr = 0
    do node = 1, nnodes
      nodenr = work1(node)
      if ( work3(nodenr) == 0 ) then
!       node is new
        lnodenr = lnodenr + 1
        work3(nodenr) = lnodenr
        work4(lnodenr) = nodenr
      end if
    end do

    newnnodes = lnodenr

!   update topology for new numbering

    do elem = 1, nelem
      work2(:,elem) = work3 ( work1 ( work2(:,elem) ) )
    end do

!   create curve

    mesh%curves(newcurve)%ndim = mesh%curves(abs(curve(1)))%ndim
    mesh%curves(newcurve)%nnodes = newnnodes
    mesh%curves(newcurve)%nelem = nelem
    mesh%curves(newcurve)%nblend = 0
    mesh%curves(newcurve)%element%elshape = &
                                 mesh%curves(abs(curve(1)))%element%elshape
    mesh%curves(newcurve)%element%p = mesh%curves(abs(curve(1)))%element%p
    mesh%curves(newcurve)%element%numnod = elnumnod
    mesh%curves(newcurve)%element%ndim = mesh%curves(abs(curve(1)))%ndim
    mesh%curves(newcurve)%element%globalshape = 'line'

    allocate( mesh%curves(newcurve)%element_blend(0) )
    mesh%curves(newcurve)%nnodes_blend = [0,mesh%curves(newcurve)%nnodes]

    mesh%curves(newcurve)%nodes = work4(:newnnodes)
    allocate( mesh%curves(newcurve)%topology(elnumnod,nelem,2) )
    mesh%curves(newcurve)%topology(:,:,1) = work2(:,:nelem)
    do elem = 1, nelem
      mesh%curves(newcurve)%topology(:,elem,2) = work4( work2(:,elem) )
    end do

    mesh%ncurves = mesh%ncurves + 1

    if ( present(newnr) ) newnr = mesh%ncurves

    deallocate ( work1, work2, work3, work4 )

  end subroutine add_to_mesh_curve


! Define a new curve having the same nodal points as an existing curve
! (argument curve(1)) but renumbered to match the numbering/topology of
! argument curve(2).

  subroutine add_to_mesh_curve_match ( mesh, curve, nr, mesh2, displace, &
    rotate, x0, mindistance )

    type(mesh_t), intent(inout) :: mesh

!   a new curve is created consisting of the nodes of the curve in argument
!   curve(1) but matching the numbering and topology of curve(2).
!   NOTE: the coordinates of the nodal points of curve(2) must be identical to
!   the coordinates of the nodal points of curve(1) (after optionally adding
!   the virtual displacement displace)
    integer, dimension(:), intent(in) :: curve

!   if present: replace a previously defined curve
    integer, intent(in), optional :: nr

!   if present curve(2) is defined in mesh2, otherwise in mesh.
    type(mesh_t), intent(in), optional :: mesh2

!   if present the coordinates of the nodes in curve(2) are "displaced"
!   virtually with vector displace for matching the coordinates of curve(1).
!   NOTE: the coordinates of the nodes are not actually changed.
    real(dp), dimension(:), intent(in), optional :: displace

!   if present the coordinates of the nodes in curve(2) are "rotated"
!   virtually through the angle specified in rotate for matching the coordinates
!   of curve(1).
!   for 2D meshes:
!   rotate(1): angle of rotation in the xy-plane
!   for 3D meshes:
!   rotate(1): angle of rotation around x-axis
!   rotate(2): angle of rotation around y-axis
!   rotate(3): angle of rotation around z-axis
!   NOTE: the coordinates of the nodes are not actually changed.
!   NOTE2: in 3D rotations are applied in an order of increasing dimension, so
!   first around x-axis, then y-axis and finally z-axis
!   NOTE3: if displace is also present, this is performed after the rotation
    real(dp), dimension(:), intent(in), optional :: rotate

!   if present perform the rotation w.r.t to point x0, otherwise rotate around
!   origin
    real(dp), dimension(:), intent(in), optional :: x0

!   if present the matching of coordinates is based on the minimum distance
!   default=.false.
    logical, intent(in), optional :: mindistance

    integer :: nnodes, newcurve, i
    integer, dimension(:), allocatable :: nodes1, nodes2, perm

!   add curve

    if ( size( curve ) /= 2 ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh/add_to_mesh_curve_match: ', &
        '  size curve array must be two'
      stop
    end if

    if ( present(mesh2) ) then
      call check_curves ( mesh2 )
    else
      call check_curves ( mesh )
    end if

    if ( present(displace) ) then
      if ( size( displace ) /= mesh%ndim ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_curve_match: ', &
          '  size of argument displace different from space dimension'
        stop
      end if
    end if

    if ( present(rotate) ) then
      if ( mesh%ndim == 1 ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_surface_match: ', &
          '  rotation not available for 1D meshes'
        stop
      end if
      if ( mesh%ndim == 2 .and. size(rotate) /= 1 ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_curve_match: ', &
          '  mesh%ndim == 2 .and. size(rotate) /= 1'
        stop
      end if
      if ( mesh%ndim == 3 .and. size(rotate) /= 3 ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_curve_match: ', &
          '  mesh%ndim == 3 .and. size(rotate) /= 3'
        stop
      end if
    end if

    if ( present(x0) ) then
      if ( size(x0) /= mesh%ndim ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_curve_match: ', &
          '  size(x0) /= mesh%ndim '
        stop
      end if
    end if

    nnodes = mesh%curves(curve(1))%nnodes
    allocate ( nodes1(nnodes), nodes2(nnodes), perm(nnodes) )

!   copy nodes from curve(1) and curve(2)

    nodes1 = mesh%curves(curve(1))%nodes

    if ( present(mesh2) ) then
      nodes2 = mesh2%curves(curve(2))%nodes
    else
      nodes2 = mesh%curves(curve(2))%nodes
    end if

    if ( present(nr) ) then

!     replace previously defined curve

      newcurve = nr

      if ( newcurve < 1 .or. newcurve > mesh%ncurves ) then
        write(*,'(/a/a,i0/)') &
          'Error: curve out of range in add_to_mesh_curve_match:', &
          ' curve < 1 or curve > ncurves = ', mesh%ncurves
        stop
      end if

      call delete_geometry ( mesh%curves(newcurve) )

    else

!     add new curve

      newcurve = mesh%ncurves + 1

      if ( newcurve > size(mesh%curves) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in add_to_mesh:', &
          ' curves > maximum = ', size(mesh%curves)
        stop
      end if

    end if

!   copy geometry to new curve

    if ( present(mesh2) ) then
      call copy ( mesh2%curves(curve(2)), mesh%curves(newcurve) )
    else
      call copy ( mesh%curves(curve(2)), mesh%curves(newcurve) )
    end if

!   renumber nodes
    if ( present(rotate) ) then
      call renumber_nodes_to_match_coor_rotate ( mesh, nodes1, nodes2, perm, &
        mesh2, displace, rotate, x0, mindistance )
    else
      call renumber_nodes_to_match_coor ( mesh, nodes1, nodes2, perm, mesh2, &
        displace, mindistance )
    end if

    if ( any( perm == 0 ) ) then
      write(*,'(/a/a,i0,1x,i0)') &
        'Error in add_to_mesh/add_to_mesh_curve_match: ', &
        ' coordinates of curves in curve array do not match: curve = ', curve
      stop
    end if

    mesh%curves(newcurve)%nodes = nodes1(perm)

    do i = 1, mesh%curves(newcurve)%nelem
      mesh%curves(newcurve)%topology(:,i,2) = mesh%curves(newcurve)%nodes( &
                                 mesh%curves(newcurve)%topology(:,i,1)   )
    end do

    if ( .not. present(nr) ) mesh%ncurves = mesh%ncurves + 1

    deallocate ( nodes1, nodes2, perm )

  contains

    subroutine check_curves ( mesha )

      type(mesh_t), intent(in) :: mesha

      if ( curve(1) <= 0 .or. curve(1) > mesh%ncurves .or. &
           curve(2) <= 0 .or. curve(2) > mesha%ncurves ) then
        write(*,'(/a/a,2(i0,1x))') &
          'Error in add_to_mesh/add_to_mesh_curve_match: ', &
          '  curve array contains non-existing curves: ', curve
        stop
      end if
      if ( mesh%curves(curve(1))%element%numnod /= &
                     mesha%curves(curve(2))%element%numnod ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_curve_match: ', &
          '  curves in curve array have different element shapes'
        stop
      end if
      if ( mesh%curves(curve(1))%nnodes /= &
                     mesha%curves(curve(2))%nnodes ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_curve_match: ', &
          '  curves in curve array have different number of nodes'
        stop
      end if
      if ( mesh%curves(curve(1))%nelem /= &
                     mesha%curves(curve(2))%nelem ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_curve_match: ', &
          '  curves in curve array have different number of nodes'
        stop
      end if
      if ( mesh%curves(curve(1))%ndim /= &
                     mesha%curves(curve(2))%ndim ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_curve_match: ', &
          '  curves in curve array have different space dimension'
        stop
      end if

    end subroutine check_curves

  end subroutine add_to_mesh_curve_match


! Define new curve based on existing groups

  subroutine add_to_mesh_curve_from_groups ( mesh, groups )

    type(mesh_t), intent(inout) :: mesh

!   a new curve is created from the existing element groups specified
!   in the array. For example groups=(/2,3/) means: create a new curve
!   consisting of elements from group 2 and 3.
!   The element types in the groups need to be consistent with curves
!   (line elements)
!   The nodal points of the composite curve are all the unique nodal points
!   in the set of given groups and the connectivity of all the elements
!   in the given groups with the mesh nodal points preserved. Thus if the
!   individual groups are connected, so will be the new curve. If the
!   individual groups are disjunct, the new curve will consist of disjunct
!   parts.
    integer, dimension(:), intent(in) :: groups

!   define new curves from existing groups

    integer :: grp, nelem, nnodes, elnumnod, newcurve, elem, gr, lnode
    integer :: gnode, gelem
!   temporary for marking global nodes
    integer, dimension(:), allocatable :: gnodes

!   add curves

    if ( size( groups ) == 0 ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  curvefromgroups array is empty'
      stop
    end if

    if ( any( groups <= 0 ) .or. any( groups > mesh%nelgrp ) ) then
      write(*,'(/a/a)',advance='no') &
        'Error in add_to_mesh: ', &
        '  curvefromgroups array contains non-existing element groups: '
      write(*,*) groups
      write(*,*)
      stop
    end if

    if ( any( mesh%element(groups)%elshape /= &
                   mesh%element(groups(1))%elshape ) ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  groups in curvefromgroups array have different element shapes'
      stop
    end if

    if ( mesh%element(groups(1))%globalshape /= 'line' ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  groups in curvefromgroups array not line elements'
      stop
    end if

    newcurve = mesh%ncurves + 1

    if ( newcurve > size(mesh%curves) ) then
      write(*,'(/a/a,i0/)') &
        'Error: maximum exceeded in add_to_mesh:', &
        ' curves > maximum = ', size(mesh%curves)
      stop
    end if

    allocate ( gnodes(mesh%nnodes) ) ! work space to bookkeep global nodes

    gnodes = 0

!   mark global nodes

    do gr = 1, size(groups)
      grp = groups(gr)
      do elem = 1, mesh%grpnumel(grp)
        gnodes ( mesh%topology(grp)%a(:,elem) ) = 1 ! mark gnode
      end do
    end do

    nnodes = count ( gnodes == 1 )

    allocate ( mesh%curves(newcurve)%nodes(nnodes) )

!   fill local nodes

    lnode = 0

    do gnode = 1, mesh%nnodes
      if ( gnodes(gnode) == 0 ) cycle  ! global node not in surface
      lnode = lnode + 1 ! new local node
      mesh%curves(newcurve)%nodes(lnode) = gnode
      gnodes(gnode) = lnode ! store lnode
    end do

    if ( lnode /= nnodes ) then
     print *, 'lnode, nnodes', lnode, nnodes
     stop 'add_to_mesh: internal error'
    end if

!   fill local and global topology

    elnumnod = mesh%element(groups(1))%numnod
    nelem = sum ( mesh%grpnumel(groups) )

    allocate( mesh%curves(newcurve)%topology(elnumnod,nelem,2) )

    gelem = 0

    do gr = 1, size(groups)
      grp = groups(gr)
      do elem = 1, mesh%grpnumel(grp)
         mesh%curves(newcurve)%topology(:,gelem+elem,1) = &
            gnodes ( mesh%topology(grp)%a(:,elem) )
         mesh%curves(newcurve)%topology(:,gelem+elem,2) = &
                     mesh%topology(grp)%a(:,elem)
      end do
      gelem = gelem + mesh%grpnumel(grp)
    end do

!   create curve

    mesh%curves(newcurve)%ndim = mesh%ndim
    mesh%curves(newcurve)%nnodes = nnodes
    mesh%curves(newcurve)%nelem = nelem
    mesh%curves(newcurve)%nblend = 0
    mesh%curves(newcurve)%element%elshape = &
                                 mesh%element(groups(1))%elshape
    mesh%curves(newcurve)%element%p = mesh%element(groups(1))%p
    mesh%curves(newcurve)%element%numnod = elnumnod
    mesh%curves(newcurve)%element%ndim = mesh%ndim
    mesh%curves(newcurve)%element%globalshape = 'line'

    allocate( mesh%curves(newcurve)%element_blend(0) )
    mesh%curves(newcurve)%nnodes_blend = [0,mesh%curves(newcurve)%nnodes]

    mesh%ncurves = mesh%ncurves + 1

    deallocate ( gnodes )

  end subroutine add_to_mesh_curve_from_groups


! Define new surface from existing surfaces

  subroutine add_to_mesh_surface ( mesh, surface, newnr )

    type(mesh_t), intent(inout) :: mesh

!   a new surface is created consisting of the existing surfaces specified in
!   the array. For example curve=(/2,6/) means: create a new surface consisting
!   of surface 2 and surface 6.
!   The nodal points of the composite surface are all the unique nodal points
!   in the set of given surfaces and the connectivity of all the elements
!   in the given surfaces with the mesh nodal points is preserved. Thus if the
!   individual surfaces are connected, so will be the composite surfaces. If the
!   individual surfaces are disjunct, the composite surface will consist of
!   disjunct parts.
    integer, dimension(:), intent(in) :: surface

!   If present, the value of newnr gives the new surface number in the mesh.
    integer, intent(out), optional :: newnr

!   define new surfaces from existing entities

    integer :: srf, nelem, nnodes, elnumnod, newsurface, elem, sr, node
    integer :: lnodenr, newnnodes, nodenr
!   temporaries for nodes
    integer, dimension(:), allocatable :: work1, work3, work4
    integer, dimension(:,:), allocatable :: work2   ! temporary for topology

!   add surfaces

    if ( size( surface ) == 0 ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  surface array is empty'
      stop
    end if

    if ( any( surface < 0 ) ) then
      write(*,'(/a/a/a/)') &
        'Error in add_to_mesh: ', &
        '  surface array contains negative surfaces; ', &
        '  It is not possible to change orientation of the surfaces '
      stop
    end if

    if ( any( surface == 0 ) .or. any( surface > mesh%nsurfaces ) ) then
      write(*,'(/a/a)',advance='no') &
        'Error in add_to_mesh: ', &
        '  surface array contains non-existing surfaces: '
      write(*,*) surface
      write(*,*)
      stop
    end if

    if ( any( mesh%surfaces(abs(surface))%element%numnod /= &
                   mesh%surfaces(abs(surface(1)))%element%numnod ) .or. &
         any( mesh%surfaces(abs(surface))%element%globalshape /= &
                   mesh%surfaces(abs(surface(1)))%element%globalshape ) ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  surfaces in surface array have different element shapes'
      stop
    end if

    newsurface = mesh%nsurfaces + 1

    if ( newsurface > size(mesh%surfaces) ) then
      write(*,'(/a/a,i0/)') &
        'Error: maximum exceeded in add_to_mesh:', &
        ' surfaces > maximum = ', size(mesh%surfaces)
      stop
    end if

!   temporary for storing all the nodes in the surfaces (including double ones)
    allocate ( work1( sum(mesh%surfaces(abs(surface))%nnodes) ) )
!   temporary for storing the topology of the combined surfaces
    allocate ( work2( mesh%surfaces(1)%element%numnod, &
                    sum(mesh%surfaces(abs(surface))%nelem)) )

    nnodes = 0
    nelem = 0

    elnumnod = mesh%surfaces(abs(surface(1)))%element%numnod

!   loop over the surfaces

    do sr = 1, size(surface)

      srf = abs(surface(sr))

!     store all nodes and new topology

      work1(nnodes+1:nnodes+mesh%surfaces(srf)%nnodes) = &
                                         mesh%surfaces(srf)%nodes
      work2(:,nelem+1:nelem+mesh%surfaces(srf)%nelem) = &
                            mesh%surfaces(srf)%topology(:,:,1) + nnodes

      nnodes = nnodes + mesh%surfaces(srf)%nnodes
      nelem = nelem + mesh%surfaces(srf)%nelem

    end do

!   check for double nodes

!   work3: store at each (global) node in the surfaces the new local node
!          number, otherwise store a zero
!   work4: store at each local node the corresponding mesh node number

    allocate ( work3(mesh%nnodes), work4(nnodes) )

    work3 = 0

    lnodenr = 0
    do node = 1, nnodes
      nodenr = work1(node)
      if ( work3(nodenr) == 0 ) then
!       node is new
        lnodenr = lnodenr + 1
        work3(nodenr) = lnodenr
        work4(lnodenr) = nodenr
      end if
    end do

    newnnodes = lnodenr

!   update topology for new numbering

    do elem = 1, nelem
      work2(:,elem) = work3 ( work1 ( work2(:,elem) ) )
    end do

!   create surface

    mesh%surfaces(newsurface)%ndim = mesh%surfaces(abs(surface(1)))%ndim
    mesh%surfaces(newsurface)%nnodes = newnnodes
    mesh%surfaces(newsurface)%nelem = nelem
    mesh%surfaces(newsurface)%nblend = 0
    mesh%surfaces(newsurface)%element%elshape = &
                               mesh%surfaces(abs(surface(1)))%element%elshape
    mesh%surfaces(newsurface)%element%p = &
                               mesh%surfaces(abs(surface(1)))%element%p
    mesh%surfaces(newsurface)%element%numnod = elnumnod
    mesh%surfaces(newsurface)%element%ndim = mesh%surfaces(abs(surface(1)))%ndim
    mesh%surfaces(newsurface)%element%globalshape = &
                           mesh%surfaces(abs(surface(1)))%element%globalshape

    allocate( mesh%surfaces(newsurface)%element_blend(0) )
    mesh%surfaces(newsurface)%nnodes_blend = &
                                      [0,mesh%surfaces(newsurface)%nnodes]

    allocate( mesh%surfaces(newsurface)%nodes(newnnodes) )
    mesh%surfaces(newsurface)%nodes = work4(:newnnodes)
    allocate( mesh%surfaces(newsurface)%topology(elnumnod,nelem,2) )
    mesh%surfaces(newsurface)%topology(:,:,1) = work2(:,:nelem)
    do elem = 1, nelem
      mesh%surfaces(newsurface)%topology(:,elem,2) = work4( work2(:,elem) )
    end do

    mesh%nsurfaces = mesh%nsurfaces + 1

    if ( present(newnr) ) newnr = mesh%nsurfaces

    deallocate ( work1, work2, work3, work4 )

  end subroutine add_to_mesh_surface


! Define a new surface having the same nodal points as an existing surface
! (argument surface(1)) but renumbered to match the numbering/topology of
! argument surface(2).

  subroutine add_to_mesh_surface_match ( mesh, surface, nr, mesh2, &
    displace, rotate, x0, mindistance )

    type(mesh_t), intent(inout) :: mesh

!   a new surface is created consisting of the nodes of the surface in argument
!   surface(1) but matching the numbering and topology of surface(2).
!   NOTE: the coordinates of the nodal points of surface(2) must be identical to
!   the coordinates of the nodal points of surface(1) (after optionally adding
!   the virtual displacement displace)
    integer, dimension(:), intent(in) :: surface

!   if present: replace a previously defined surface
    integer, intent(in), optional :: nr

!   if present surface(2) is defined in mesh2, otherwise in mesh.
    type(mesh_t), intent(in), optional :: mesh2

!   if present the coordinates of the nodes in surface(2) are "displaced"
!   virtually with vector displace for matching the coordinates of surface(1).
!   NOTE: the coordinates of the nodes are not actually changed.
    real(dp), dimension(:), intent(in), optional :: displace

!   if present the coordinates of the nodes in surface(2) are "rotated"
!   virtually through the angle specified in rotate for matching the coordinates
!   of surface(1).
!   for 2D meshes:
!   rotate(1): angle of rotation in the xy-plane
!   for 3D meshes:
!   rotate(1): angle of rotation around x-axis
!   rotate(2): angle of rotation around y-axis
!   rotate(3): angle of rotation around z-axis
!   NOTE: the coordinates of the nodes are not actually changed.
!   NOTE2: in 3D rotations are applied in an order of increasing dimension, so
!   first around x-axis, then y-axis and finally z-axis
!   NOTE3: if displace is also present, this is performed after the rotation
    real(dp), dimension(:), intent(in), optional :: rotate

!   if present perform the rotation w.r.t to point x0, otherwise rotate around
!   origin
    real(dp), dimension(:), intent(in), optional :: x0

!   if present the matching of coordinates is based on the minimum distance
!   default=.false.
    logical, intent(in), optional :: mindistance

    integer :: nnodes, newsurface, i
    integer, dimension(:), allocatable :: nodes1, nodes2, perm

!   add surface

    if ( size( surface ) /= 2 ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh/add_to_mesh_surface_match: ', &
        '  size surface array must be two'
      stop
    end if

    if ( present(mesh2) ) then
      call check_surfaces ( mesh2 )
    else
      call check_surfaces ( mesh )
    end if

    if ( present(displace) ) then
      if ( size( displace ) /= mesh%ndim ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_surface_match: ', &
          '  size of argument displace different from space dimension'
        stop
      end if
    end if

    if ( present(rotate) ) then
      if ( mesh%ndim == 1 ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_surface_match: ', &
          '  rotation not available for 1D meshes'
        stop
      end if
      if ( mesh%ndim == 2 .and. size(rotate) /= 1 ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_surface_match: ', &
          '  mesh%ndim == 2 .and. size(rotate) /= 1'
        stop
      end if
      if ( mesh%ndim == 3 .and. size(rotate) /= 3 ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_surface_match: ', &
          '  mesh%ndim == 3 .and. size(rotate) /= 3'
        stop
      end if
    end if

    if ( present(x0) ) then
      if ( size(x0) /= mesh%ndim ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_surface_match: ', &
          '  size(x0) /= mesh%ndim '
        stop
      end if
    end if

    nnodes = mesh%surfaces(surface(1))%nnodes
    allocate ( nodes1(nnodes), nodes2(nnodes), perm(nnodes) )

!   copy nodes from surface(1) and surface(2)

    nodes1 = mesh%surfaces(surface(1))%nodes

    if ( present(mesh2) ) then
      nodes2 = mesh2%surfaces(surface(2))%nodes
    else
      nodes2 = mesh%surfaces(surface(2))%nodes
    end if

    if ( present(nr) ) then

!     replace previously defined surface

      newsurface = nr

      if ( newsurface < 1 .or. newsurface > mesh%nsurfaces ) then
        write(*,'(/a/a,i0/)') &
          'Error: surface out of range in add_to_mesh_surface_match:', &
          ' surface < 1 or surface > nsurfaces = ', mesh%nsurfaces
        stop
      end if

      call delete_geometry ( mesh%surfaces(newsurface) )

    else

!     add new surface

      newsurface = mesh%nsurfaces + 1

      if ( newsurface > size(mesh%surfaces) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in add_to_mesh:', &
          ' surfaces > maximum = ', size(mesh%surfaces)
        stop
      end if

    end if

!   copy geometry to new surface

    if ( present(mesh2) ) then
      call copy ( mesh2%surfaces(surface(2)), mesh%surfaces(newsurface) )
    else
      call copy ( mesh%surfaces(surface(2)), mesh%surfaces(newsurface) )
    end if

!   renumber nodes

    if ( present(rotate) ) then
      call renumber_nodes_to_match_coor_rotate ( mesh, nodes1, nodes2, perm, &
        mesh2, displace, rotate, x0, mindistance )
    else
      call renumber_nodes_to_match_coor ( mesh, nodes1, nodes2, perm, mesh2, &
        displace, mindistance )
    end if

    if ( any( perm == 0 ) ) then
      write(*,'(/a/a,i0,1x,i0)') &
        'Error in add_to_mesh/add_to_mesh_surface_match: ', &
        ' coordinates of surfaces in surface array do not match: surface = ', &
        surface
      stop
    end if

    mesh%surfaces(newsurface)%nodes = nodes1(perm)

    do i = 1, mesh%surfaces(newsurface)%nelem
      mesh%surfaces(newsurface)%topology(:,i,2) = &
        mesh%surfaces(newsurface)%nodes( &
                                 mesh%surfaces(newsurface)%topology(:,i,1) )
    end do

    if ( .not. present(nr) ) mesh%nsurfaces = mesh%nsurfaces + 1

    deallocate ( nodes1, nodes2, perm )

  contains

    subroutine check_surfaces ( mesha )

      type(mesh_t), intent(in) :: mesha

      if ( surface(1) <= 0 .or. surface(1) > mesh%nsurfaces .or. &
           surface(2) <= 0 .or. surface(2) > mesha%nsurfaces ) then
        write(*,'(/a/a,2(i0,1x))') &
          'Error in add_to_mesh/add_to_mesh_surface_match: ', &
          '  surface array contains non-existing surfaces: ', surface
        stop
      end if
      if ( mesh%surfaces(surface(1))%element%numnod /= &
                     mesha%surfaces(surface(2))%element%numnod ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_surface_match: ', &
          '  surfaces in surface array have different element shapes'
        stop
      end if
      if ( mesh%surfaces(surface(1))%nnodes /= &
                     mesha%surfaces(surface(2))%nnodes ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_surface_match: ', &
          '  surfaces in surface array have different number of nodes'
        stop
      end if
      if ( mesh%surfaces(surface(1))%nelem /= &
                     mesha%surfaces(surface(2))%nelem ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_surface_match: ', &
          '  surfaces in surface array have different number of nodes'
        stop
      end if
      if ( mesh%surfaces(surface(1))%ndim /= &
                     mesha%surfaces(surface(2))%ndim ) then
        write(*,'(/a/a/)') &
          'Error in add_to_mesh/add_to_mesh_surface_match: ', &
          '  surfaces in surface array have different space dimension'
        stop
      end if

    end subroutine check_surfaces

  end subroutine add_to_mesh_surface_match


! Define new surface based on existing groups

  subroutine add_to_mesh_surface_from_groups ( mesh, groups )

    type(mesh_t), intent(inout) :: mesh

!   a new surface is created from the existing element groups specified
!   in the array. For example groups=(/2,3/) means: create a new surface
!   consisting of elements from group 2 and 3.
!   The element types in the groups need to be consistent with surfaces
!   (surface elements: triangle or quadrilateral).
!   The nodal points of the composite surface are all the unique nodal points
!   in the set of given groups and the connectivity of all the elements
!   in the given groups with the mesh nodal points preserved. Thus if the
!   individual groups are connected, so will be the new surface. If the
!   individual groups are disjunct, the new surface will consist of disjunct
!   parts.
    integer, dimension(:), intent(in) :: groups

!   define new surfaces from existing groups

    integer :: grp, nelem, nnodes, elnumnod, newsurface, elem, gr, lnode
    integer :: gnode, gelem
!   temporaries for nodes
    integer, dimension(:), allocatable :: gnodes

!   add surfaces

    if ( size( groups ) == 0 ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  surfacefromgroups array is empty'
      stop
    end if

    if ( any( groups <= 0 ) .or. any( groups > mesh%nelgrp ) ) then
      write(*,'(/a/a)',advance='no') &
        'Error in add_to_mesh: ', &
        '  surfacefromgroups array contains non-existing element groups: '
      write(*,*) groups
      write(*,*)
      stop
    end if

    if ( any( mesh%element(groups)%elshape /= &
                   mesh%element(groups(1))%elshape ) ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  groups in surfacefromgroups array have different element shapes'
      stop
    end if

    if ( all( mesh%element(groups(1))%globalshape /= &
                                 [ 'triangle     ',  'quadrilateral' ] ) ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  groups in surfacefromgroups array not surface elements'
      stop
    end if

    newsurface = mesh%nsurfaces + 1

    if ( newsurface > size(mesh%surfaces) ) then
      write(*,'(/a/a,i0/)') &
        'Error: maximum exceeded in add_to_mesh:', &
        ' surfaces > maximum = ', size(mesh%surfaces)
      stop
    end if

    allocate ( gnodes(mesh%nnodes) ) ! work space to bookkeep global nodes

    gnodes = 0

!   mark global nodes

    do gr = 1, size(groups)
      grp = groups(gr)
      do elem = 1, mesh%grpnumel(grp)
        gnodes ( mesh%topology(grp)%a(:,elem) ) = 1 ! mark gnode
      end do
    end do

    nnodes = count ( gnodes == 1 )

    allocate ( mesh%surfaces(newsurface)%nodes(nnodes) )

!   fill local nodes

    lnode = 0

    do gnode = 1, mesh%nnodes
      if ( gnodes(gnode) == 0 ) cycle  ! global node not in surface
      lnode = lnode + 1 ! new local node
      mesh%surfaces(newsurface)%nodes(lnode) = gnode
      gnodes(gnode) = lnode ! store lnode
    end do

    if ( lnode /= nnodes ) then
     print *, 'lnode, nnodes', lnode, nnodes
     stop 'add_to_mesh: internal error'
    end if

!   fill local and global topology

    elnumnod = mesh%element(groups(1))%numnod
    nelem = sum ( mesh%grpnumel(groups) )

    allocate( mesh%surfaces(newsurface)%topology(elnumnod,nelem,2) )

    gelem = 0

    do gr = 1, size(groups)
      grp = groups(gr)
      do elem = 1, mesh%grpnumel(grp)
         mesh%surfaces(newsurface)%topology(:,gelem+elem,1) = &
            gnodes ( mesh%topology(grp)%a(:,elem) )
         mesh%surfaces(newsurface)%topology(:,gelem+elem,2) = &
                     mesh%topology(grp)%a(:,elem)
      end do
      gelem = gelem + mesh%grpnumel(grp)
    end do

!   create surface

    mesh%surfaces(newsurface)%ndim = mesh%ndim
    mesh%surfaces(newsurface)%nnodes = nnodes
    mesh%surfaces(newsurface)%nelem = nelem
    mesh%surfaces(newsurface)%nblend = 0
    mesh%surfaces(newsurface)%element%elshape = &
                                 mesh%element(groups(1))%elshape
    mesh%surfaces(newsurface)%element%p = mesh%element(groups(1))%p
    mesh%surfaces(newsurface)%element%numnod = elnumnod
    mesh%surfaces(newsurface)%element%ndim = mesh%ndim
    mesh%surfaces(newsurface)%element%globalshape = &
                                 mesh%element(groups(1))%globalshape

    allocate( mesh%surfaces(newsurface)%element_blend(0) )
    mesh%surfaces(newsurface)%nnodes_blend = &
                                      [0,mesh%surfaces(newsurface)%nnodes]

    mesh%nsurfaces = mesh%nsurfaces + 1

    deallocate ( gnodes )

  end subroutine add_to_mesh_surface_from_groups


! Define new volume based on existing groups

  subroutine add_to_mesh_volume_from_groups ( mesh, groups )

    type(mesh_t), intent(inout) :: mesh

!   a new volume is created from the existing element groups specified
!   in the array. For example groups=(/2,3/) means: create a new volume
!   consisting of elements from group 2 and 3.
!   The element types in the groups need to be consistent with volumes
!   (volume elements: triangle or quadrilateral).
!   The nodal points of the composite volume are all the unique nodal points
!   in the set of given groups and the connectivity of all the elements
!   in the given groups with the mesh nodal points preserved. Thus if the
!   individual groups are connected, so will be the new volume. If the
!   individual groups are disjunct, the new volume will consist of disjunct
!   parts.
    integer, dimension(:), intent(in) :: groups

!   define new volumes from existing groups

    integer :: grp, nelem, nnodes, elnumnod, newvolume, elem, gr, lnode
    integer :: gnode, gelem
!   temporaries for nodes
    integer, dimension(:), allocatable :: gnodes

!   add volumes

    if ( size( groups ) == 0 ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  volumefromgroups array is empty'
      stop
    end if

    if ( any( groups <= 0 ) .or. any( groups > mesh%nelgrp ) ) then
      write(*,'(/a/a)',advance='no') &
        'Error in add_to_mesh: ', &
        '  volumefromgroups array contains non-existing element groups: '
      write(*,*) groups
      write(*,*)
      stop
    end if

    if ( any( mesh%element(groups)%elshape /= &
                   mesh%element(groups(1))%elshape ) ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  groups in volumefromgroups array have different element shapes'
      stop
    end if

    if ( all( mesh%element(groups(1))%globalshape /= &
               [ 'hexahedron   ', 'tetrahedron  ', 'prism        ', &
                 'pyramid      ' ] ) ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh: ', &
        '  groups in volumefromgroups array not volume elements'
      stop
    end if

    newvolume = mesh%nvolumes + 1

    if ( newvolume > size(mesh%volumes) ) then
      write(*,'(/a/a,i0/)') &
        'Error: maximum exceeded in add_to_mesh:', &
        ' volumes > maximum = ', size(mesh%volumes)
      stop
    end if

    allocate ( gnodes(mesh%nnodes) ) ! work space to bookkeep global nodes

    gnodes = 0

!   mark global nodes

    do gr = 1, size(groups)
      grp = groups(gr)
      do elem = 1, mesh%grpnumel(grp)
        gnodes ( mesh%topology(grp)%a(:,elem) ) = 1 ! mark gnode
      end do
    end do

    nnodes = count ( gnodes == 1 )

    allocate ( mesh%volumes(newvolume)%nodes(nnodes) )

!   fill local nodes

    lnode = 0

    do gnode = 1, mesh%nnodes
      if ( gnodes(gnode) == 0 ) cycle  ! global node not in volume
      lnode = lnode + 1 ! new local node
      mesh%volumes(newvolume)%nodes(lnode) = gnode
      gnodes(gnode) = lnode ! store lnode
    end do

    if ( lnode /= nnodes ) then
     print *, 'lnode, nnodes', lnode, nnodes
     stop 'add_to_mesh: internal error'
    end if

!   fill local and global topology

    elnumnod = mesh%element(groups(1))%numnod
    nelem = sum ( mesh%grpnumel(groups) )

    allocate( mesh%volumes(newvolume)%topology(elnumnod,nelem,2) )

    gelem = 0

    do gr = 1, size(groups)
      grp = groups(gr)
      do elem = 1, mesh%grpnumel(grp)
         mesh%volumes(newvolume)%topology(:,gelem+elem,1) = &
            gnodes ( mesh%topology(grp)%a(:,elem) )
         mesh%volumes(newvolume)%topology(:,gelem+elem,2) = &
                     mesh%topology(grp)%a(:,elem)
      end do
      gelem = gelem + mesh%grpnumel(grp)
    end do

!   create volume

    mesh%volumes(newvolume)%ndim = mesh%ndim
    mesh%volumes(newvolume)%nnodes = nnodes
    mesh%volumes(newvolume)%nelem = nelem
    mesh%volumes(newvolume)%nblend = 0
    mesh%volumes(newvolume)%element%elshape = &
                                 mesh%element(groups(1))%elshape
    mesh%volumes(newvolume)%element%p = mesh%element(groups(1))%p
    mesh%volumes(newvolume)%element%numnod = elnumnod
    mesh%volumes(newvolume)%element%ndim = mesh%ndim
    mesh%volumes(newvolume)%element%globalshape = &
                                 mesh%element(groups(1))%globalshape

    allocate( mesh%volumes(newvolume)%element_blend(0) )
    mesh%volumes(newvolume)%nnodes_blend = [0,mesh%volumes(newvolume)%nnodes]

    mesh%nvolumes = mesh%nvolumes + 1

    deallocate ( gnodes )

  end subroutine add_to_mesh_volume_from_groups


! Define new object based on points only. (user subroutine version)

  subroutine add_to_mesh_object_points_sub ( mesh, nnodes, typeofobject, &
    objectsub, objectcoornr, nr )

    type(mesh_t), intent(inout) :: mesh

    integer, intent(in) :: nnodes, typeofobject

!   user subroutine for defining the coordinates of the points in the object
!   and the number used in the subroutines
    interface
      subroutine objectsub ( objectcoornr, coor )
        use kind_defs_m
        implicit none
        integer, intent(in) :: objectcoornr
        real(dp), dimension(:,:), intent(inout) :: coor
      end subroutine objectsub
    end interface
    integer, intent(in) :: objectcoornr

!   if present: change a previously defined object
    integer, intent(in), optional :: nr


    integer :: newobject


    if ( present(nr) ) then

!     change previously defined object

      newobject = nr

      if ( newobject < 1 .or. newobject > mesh%nobjects ) then
        write(*,'(/a/a,i0/)') &
          'Error: object out of range in add_to_mesh_object_points_sub:', &
          ' object < 1 or object > nobjects = ', mesh%nobjects
        stop
      end if

      call delete_object ( mesh%objects(newobject) )

    else

!     add new object

      newobject = mesh%nobjects + 1

      if ( newobject > size(mesh%objects) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in add_to_mesh_object_points_sub:', &
          ' object > maximum = ', size(mesh%objects)
        stop
      end if

    end if

    mesh%objects(newobject)%typeofobject = typeofobject

    mesh%objects(newobject)%ndim = mesh%ndim
    mesh%objects(newobject)%nnodes = nnodes

    allocate ( &
      mesh%objects(newobject)%coor(nnodes,mesh%objects(newobject)%ndim), &
      mesh%objects(newobject)%refcoor(nnodes,mesh%objects(newobject)%ndim), &
      mesh%objects(newobject)%grpelm(nnodes,2) )

    if ( typeofobject == 2 ) then
      allocate ( &
        mesh%objects(newobject)%refcoor2(nnodes,mesh%objects(newobject)%ndim), &
        mesh%objects(newobject)%grpelm2(nnodes,2) )
    end if

!   get coor from user subroutine

    call objectsub ( objectcoornr, mesh%objects(newobject)%coor )

    if ( .not. present(nr) ) mesh%nobjects = mesh%nobjects + 1

  end subroutine add_to_mesh_object_points_sub


! Define new object based on points only. (coordinate array version)

  subroutine add_to_mesh_object_points_coor ( mesh, nnodes, typeofobject, &
    coor, nr )

    type(mesh_t), intent(inout) :: mesh

!   number of nodes in the object
!   if nnodes > size(coor,1) memory allocation is based on nnodes.
    integer, intent(in) :: nnodes

    integer, intent(in) :: typeofobject

!   coor(:,ndim). The number of nodes in coor can be smaller than nnodes.
!   the object created has size(coor,1) nodes.
    real(dp), dimension(:,:), intent(in) :: coor

!   if present: change a previously defined object
    integer, intent(in), optional :: nr


    integer :: newobject, ndim, nnodesalloc


    if ( present(nr) ) then

!     change previously defined object

      newobject = nr

      if ( newobject < 1 .or. newobject > mesh%nobjects ) then
        write(*,'(/a/a,i0/)') &
          'Error: object out of range in add_to_mesh_object_points_coor:', &
          ' object < 1 or object > nobjects = ', mesh%nobjects
        stop
      end if

      call delete_object ( mesh%objects(newobject) )

    else

!     add new object

      newobject = mesh%nobjects + 1

      if ( newobject > size(mesh%objects) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in add_to_mesh_object_points_coor:', &
          ' object > maximum = ', size(mesh%objects)
        stop
      end if

    end if

    if ( size(coor,2) /= mesh%ndim ) then
      write(*,'(/a/a/)') &
        'Error in add_to_mesh_object_points_coor: ', &
        ' dimension of coor array is not the same as mesh '
      stop
    end if

    mesh%objects(newobject)%typeofobject = typeofobject

    ndim = mesh%ndim
    mesh%objects(newobject)%ndim = ndim
!   actual number of nodes
    mesh%objects(newobject)%nnodes = min(size(coor,1),nnodes)

!   allocate memory based on the specified nnodes
    nnodesalloc = max(size(coor,1),nnodes)

    allocate ( &
      mesh%objects(newobject)%coor(nnodesalloc,ndim), &
      mesh%objects(newobject)%refcoor(nnodesalloc,ndim), &
      mesh%objects(newobject)%grpelm(nnodesalloc,2) )

    if ( typeofobject == 2 ) then
      allocate ( &
        mesh%objects(newobject)%refcoor2(nnodesalloc,ndim), &
        mesh%objects(newobject)%grpelm2(nnodesalloc,2) )
    end if

!   get coor from coor array

    mesh%objects(newobject)%coor(1:size(coor,1),:) = coor

!   initialize rest of coordinates to zero (to avoid uninitialized values)

    mesh%objects(newobject)%coor(size(coor,1)+1:nnodesalloc,:) = 0

    if ( .not. present(nr) ) mesh%nobjects = mesh%nobjects + 1

  end subroutine add_to_mesh_object_points_coor



! Define new object based on a mesh.

  subroutine add_to_mesh_object_mesh ( mesh, typeofobject, objectmesh, topol, &
    excludenodes, nr )

    type(mesh_t), intent(inout) :: mesh

    integer, intent(in) :: typeofobject

!   this mesh defines the object. It must contain a single element group only.
    type(mesh_t), intent(in) :: objectmesh

!   add topology yes/no?
    logical, intent( in) :: topol

    logical, dimension(:), intent(in) :: excludenodes

!   if present: change a previously defined object
    integer, intent(in), optional :: nr


    integer :: newobject, nnodes, node
    integer, allocatable, dimension(:) :: nodes


    if ( objectmesh%nelgrp > 1 ) then
      write(*,'(/a/a/)') &
        'Error add_to_mesh:', &
        ' objectmesh contains more than one element group.'
      stop
    end if

    if ( present(nr) ) then

!     change previously defined object

      newobject = nr

      if ( newobject < 1 .or. newobject > mesh%nobjects ) then
        write(*,'(/a/a,i0/)') &
          'Error: object out of range in add_to_mesh_object_mesh:', &
          ' object < 1 or object > nobjects = ', mesh%nobjects
        stop
      end if

      call delete_object ( mesh%objects(newobject) )

    else

!     add new object

      newobject = mesh%nobjects + 1

      if ( newobject > size(mesh%objects) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in add_to_mesh_object_mesh:', &
          ' object > maximum = ', size(mesh%objects)
        stop
      end if

    end if

    mesh%objects(newobject)%typeofobject = typeofobject

    mesh%objects(newobject)%ndim = objectmesh%ndim

!   collect nodes

    allocate ( nodes(objectmesh%nnodes) )

    nnodes = 0
    do node = 1, objectmesh%nnodes
      if ( excludenodes ( node ) ) cycle
      nnodes = nnodes + 1
      nodes(nnodes) = node
    end do
    mesh%objects(newobject)%nnodes = nnodes

    allocate ( &
      mesh%objects(newobject)%coor(nnodes,mesh%objects(newobject)%ndim), &
      mesh%objects(newobject)%refcoor(nnodes,mesh%objects(newobject)%ndim), &
      mesh%objects(newobject)%grpelm(nnodes,2) )

    if ( typeofobject == 2 ) then
      allocate ( &
        mesh%objects(newobject)%refcoor2(nnodes,mesh%objects(newobject)%ndim), &
        mesh%objects(newobject)%grpelm2(nnodes,2) )
    end if

!   get coor

    mesh%objects(newobject)%coor = objectmesh%coor(nodes(1:nnodes),:)

!   topology

    if ( topol ) then

      if ( nnodes /= objectmesh%nnodes ) then
        write(*,'(/a/a/a/)') &
          'Error add_to_mesh:', &
          ' mesh with excluded nodes used for ', &
          ' making objects cannot have a topology'
        stop
      end if

      mesh%objects(newobject)%topol = .true.
      mesh%objects(newobject)%nelem = objectmesh%nelem
      mesh%objects(newobject)%element%elshape = objectmesh%element(1)%elshape
      mesh%objects(newobject)%element%p = objectmesh%element(1)%p
      mesh%objects(newobject)%element%numnod = objectmesh%element(1)%numnod
      mesh%objects(newobject)%element%ndim = objectmesh%element(1)%ndim
      mesh%objects(newobject)%element%globalshape = &
                                       objectmesh%element(1)%globalshape
      mesh%objects(newobject)%elnumnod = objectmesh%element(1)%numnod

      allocate( mesh%objects(newobject)%topology(&
        &mesh%objects(newobject)%elnumnod,mesh%objects(newobject)%nelem) )

      mesh%objects(newobject)%topology = objectmesh%topology(1)%a

    end if

    if ( .not. present(nr) ) mesh%nobjects = mesh%nobjects + 1

    deallocate ( nodes )

  end subroutine add_to_mesh_object_mesh


! Define new object based on a geometry.

  subroutine add_to_mesh_object_geometry ( mesh, typeofobject, objectgeometry, &
    coor, topol, excludenodes, nr )

    type(mesh_t), intent(inout) :: mesh

    integer, intent(in) :: typeofobject

!   this geometry defines the object.
    type(geometry_t), intent(in) :: objectgeometry

!   the coordinates coor(1:nnodes,1:ndim) of the corresponding mesh of the
!   geometry
    real(dp), dimension(:,:), intent(in) :: coor

!   add topology yes/no?
    logical, intent(in) :: topol

    logical, dimension(:), intent(in) :: excludenodes

!   if present: change a previously defined object
    integer, intent(in), optional :: nr


    integer :: newobject, nnodes, node
    integer, allocatable, dimension(:) :: nodes


    if ( present(nr) ) then

!     change previously defined object

      newobject = nr

      if ( newobject < 1 .or. newobject > mesh%nobjects ) then
        write(*,'(/a/a,i0/)') &
          'Error: object out of range in add_to_mesh_object_geometry:', &
          ' object < 1 or object > nobjects = ', mesh%nobjects
        stop
      end if

      call delete_object ( mesh%objects(newobject) )

    else

!     add new object

      newobject = mesh%nobjects + 1

      if ( newobject > size(mesh%objects) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in add_to_mesh_object_geometry:', &
          ' object > maximum = ', size(mesh%objects)
        stop
      end if

    end if

    mesh%objects(newobject)%typeofobject = typeofobject

    mesh%objects(newobject)%ndim = objectgeometry%ndim

!   collect nodes

    allocate ( nodes(objectgeometry%nnodes) )

    nnodes = 0
    do node = 1, objectgeometry%nnodes
      if ( excludenodes ( objectgeometry%nodes(node) ) ) cycle
      nnodes = nnodes + 1
      nodes(nnodes) = objectgeometry%nodes(node)
    end do
    mesh%objects(newobject)%nnodes = nnodes

    allocate ( &
      mesh%objects(newobject)%coor(nnodes,mesh%objects(newobject)%ndim), &
      mesh%objects(newobject)%refcoor(nnodes,mesh%objects(newobject)%ndim), &
      mesh%objects(newobject)%grpelm(nnodes,2) )

    if ( typeofobject == 2 ) then
      allocate ( &
        mesh%objects(newobject)%refcoor2(nnodes,mesh%objects(newobject)%ndim), &
        mesh%objects(newobject)%grpelm2(nnodes,2) )
    end if

!   get coor

    mesh%objects(newobject)%coor = coor(nodes(1:nnodes),:)

!   topology

    if ( topol ) then

      if ( nnodes /= objectgeometry%nnodes ) then
        write(*,'(/a/a/a/)') &
          'Error add_to_mesh:', &
          ' geometries (curve,surface) with excluded nodes used for ', &
          ' making objects cannot have a topology'
        stop
      end if

      mesh%objects(newobject)%topol = .true.
      mesh%objects(newobject)%nelem = objectgeometry%nelem
      mesh%objects(newobject)%element%elshape = objectgeometry%element%elshape
      mesh%objects(newobject)%element%p = objectgeometry%element%p
      mesh%objects(newobject)%element%numnod = objectgeometry%element%numnod
      mesh%objects(newobject)%element%ndim = objectgeometry%element%ndim
      mesh%objects(newobject)%element%globalshape = &
                                       objectgeometry%element%globalshape
      mesh%objects(newobject)%elnumnod = objectgeometry%element%numnod

      allocate( mesh%objects(newobject)%topology(&
        &mesh%objects(newobject)%elnumnod,mesh%objects(newobject)%nelem) )

      mesh%objects(newobject)%topology = objectgeometry%topology(:,:,1)

    end if

    if ( .not. present(nr) ) mesh%nobjects = mesh%nobjects + 1

    deallocate ( nodes )

  end subroutine add_to_mesh_object_geometry


! Define new nodeset based on specified nodes

  subroutine add_to_mesh_nodeset_from_nodes ( mesh, nodes, nr )

    type(mesh_t), intent(inout) :: mesh

!   the nodes for creating the nodeset
!   NOTE: nodes that are out of range (<1 or>mesh%nnodes) are ignored.
    integer, dimension(:), intent(in) :: nodes

!   if present: change a previously defined nodeset
    integer, intent(in), optional :: nr


    integer :: newnodeset, nnodes, node
    integer, allocatable, dimension(:) :: sortnodes, newnodes

    allocate ( sortnodes(size(nodes)), newnodes(size(nodes)) )

    if ( present(nr) ) then

!     change previously defined nodeset

      newnodeset = nr

      if ( newnodeset < 1 .or. newnodeset > mesh%nnodesets ) then
        write(*,'(/a/a,i0/)') &
          'Error: nodeset out of range in add_to_mesh_nodeset_from_nodes:', &
          ' nodeset < 1 or nodeset > nnodesets = ', mesh%nnodesets
        stop
      end if

      deallocate ( mesh%nodesets(newnodeset)%a )

    else

!     add new nodeset

      newnodeset = mesh%nnodesets + 1

      if ( newnodeset > size(mesh%nodesets) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in add_to_mesh_nodeset_from_nodes:', &
          ' nodesets > maximum = ', size(mesh%nodesets)
        stop
      end if

    end if

    sortnodes = nodes
    call sort ( sortnodes )

!   collect nodes

    nnodes = 0

    do node = 1, size(nodes)

      if ( node < size(nodes) ) then
        if ( sortnodes(node+1) == sortnodes(node) ) then
!         node present more than once
          cycle
        end if
      end if

!     ignore nodes out of range
      if ( sortnodes(node) < 1 .or. sortnodes(node) > mesh%nnodes ) cycle

      nnodes = nnodes + 1

      newnodes(nnodes) = sortnodes(node)

    end do

    mesh%nodesets(newnodeset)%a = newnodes(1:nnodes)

    if ( .not. present(nr) ) mesh%nnodesets = mesh%nnodesets + 1

    deallocate ( sortnodes, newnodes )

  end subroutine add_to_mesh_nodeset_from_nodes


! Define new elementset based on specified (group,elements)

  subroutine add_to_mesh_elementset_from_elements &
                ( mesh, elements, group, nr, gnr )

    type(mesh_t), intent(inout) :: mesh

!   the elements for creating the elementset
!   NOTE: elements that are out of range (<1 or>mesh%grpnumel(group)) are
!         ignored.
    integer, dimension(:), intent(in) :: elements

!   if present the group number as used for setting an elementset
!   default=1
    integer, intent(in), optional :: group

!   if present: change a previously defined elementset
    integer, intent(in), optional :: nr

!   if present: change a group in a previously defined elementset
    integer, intent(in), optional :: gnr


    logical :: new
    integer :: newelementset, nelements, element, lgroup, elgrp
    integer, allocatable, dimension(:) :: sortelements, newelements

    allocate ( sortelements(size(elements)), newelements(size(elements)) )

    lgroup = set_optional ( variable=group, default=1 )

    if ( lgroup < 1 .or. lgroup > mesh%nelgrp ) then
      write(*,'(/a/a,i0/)') &
        'Error in add_to_mesh_elementset_from_elements: group incorrect, ', &
        'group < 1 or group > number of groups = ', mesh%nelgrp
      stop
    end if

    new = .true.

    if ( present(nr) ) then

!     change previously defined elementset

      newelementset = nr

      if ( newelementset < 1 .or. newelementset > mesh%nelementsets ) then
        write(*,'(/2a/a,i0/)') &
          'Error: elementset out of range in ;', &
          ' add_to_mesh_elementset_from_elements:', &
          ' elementset < 1 or elementset > nelementsets = ', mesh%nelementsets
        stop
      end if

      call delete_elementset ( mesh%elementsets(nr) )

    else if ( present(gnr) ) then

!     change group in previously defined elementset

      newelementset = gnr

      if ( newelementset < 1 .or. newelementset > mesh%nelementsets ) then
        write(*,'(/2a/a,i0/)') &
          'Error: elementset out of range in ', &
          ' add_to_mesh_elementset_from_elements:', &
          ' elementset < 1 or elementset > nelementsets = ', mesh%nelementsets
        stop
      end if

      mesh%elementsets(newelementset)%nodes_created = .false.
      mesh%elementsets(newelementset)%nnodes = 0
      deallocate ( mesh%elementsets(newelementset)%nodes )
      deallocate ( mesh%elementsets(newelementset)%elements(lgroup)%a )

      new = .false.

    else

!     add new elementset

      newelementset = mesh%nelementsets + 1

      if ( newelementset > size(mesh%elementsets) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in add_to_mesh_elementset_from_elements:', &
          ' elementset > maximum = ', size(mesh%elementsets)
        stop
      end if

    end if

    allocate ( mesh%elementsets(newelementset)%nodes(0) ) ! no nodes present

    if ( new ) then

!     new elementset; allocate and initialize data

      allocate ( mesh%elementsets(newelementset)%grpnumel(mesh%nelgrp) )
      allocate ( mesh%elementsets(newelementset)%elements(mesh%nelgrp) )
      do elgrp = 1, mesh%nelgrp
        if ( elgrp == lgroup ) cycle
        mesh%elementsets(newelementset)%grpnumel(elgrp) = 0
        allocate ( mesh%elementsets(newelementset)%elements(elgrp)%a(0) )
      end do

    end if

    sortelements = elements
    call sort ( sortelements )

!   collect elements

    nelements = 0

    do element = 1, size(elements)

      if ( element < size(elements) ) then
        if ( sortelements(element+1) == sortelements(element) ) then
!         element present more than once
          cycle
        end if
      end if

!     ignore elements out of range
      if ( sortelements(element) < 1 .or. &
           sortelements(element) > mesh%grpnumel(lgroup) ) cycle

      nelements = nelements + 1

      newelements(nelements) = sortelements(element)

    end do

    mesh%elementsets(newelementset)%grpnumel(lgroup) = nelements
    mesh%elementsets(newelementset)%elements(lgroup)%a = &
                                                  newelements(1:nelements)

    mesh%elementsets(newelementset)%nelem = &
                     sum(mesh%elementsets(newelementset)%grpnumel)

    if ( .not. ( present(nr) .or. present(gnr) ) ) &
                            mesh%nelementsets = mesh%nelementsets + 1

    deallocate ( sortelements, newelements )

  end subroutine add_to_mesh_elementset_from_elements


! Add nodes to elementset

  subroutine add_to_mesh_add_nodes_to_elementset ( mesh, elementset )

    type(mesh_t), intent(inout) :: mesh

!   add nodes to this elementset
    integer, intent(in) :: elementset


    integer :: elgrp, nnodes, n, nne, node, elem, elemnr
    integer, dimension(:), allocatable :: nodes, sortnodes


    if ( elementset < 1 .or. elementset > mesh%nelementsets ) then
      write(*,'(/2a/a,i0/)') &
        'Error: elementset out of range in ', &
        ' add_to_mesh_add_nodes_to_elementset:', &
        ' elementset < 1 or elementset > nelementsets = ', mesh%nelementsets
      stop
    end if

    deallocate ( mesh%elementsets(elementset)%nodes )

!   collect nodes

    n = sum(mesh%elementsets(elementset)%grpnumel * mesh%element(:)%numnod)
    allocate ( nodes(n), sortnodes(n) )

    nnodes = 0

    do elgrp = 1, mesh%nelgrp
      nne = mesh%element(elgrp)%numnod
      do elem = 1, mesh%elementsets(elementset)%grpnumel(elgrp)
        elemnr = mesh%elementsets(elementset)%elements(elgrp)%a(elem)
        nodes(nnodes+1:nnodes+nne)= mesh%topology(elgrp)%a(:,elemnr)
        nnodes = nnodes + nne
      end do
    end do

    if ( nnodes /= n ) stop 'error'

    sortnodes = nodes
    call sort ( sortnodes )

    nnodes = 0

    do node = 1, n

      if ( node < n ) then
        if ( sortnodes(node+1) == sortnodes(node) ) then
!         node present more than once
          cycle
        end if
      end if

      nnodes = nnodes + 1

      nodes(nnodes) = sortnodes(node)

    end do

    mesh%elementsets(elementset)%nodes_created = .true.
    mesh%elementsets(elementset)%nnodes = nnodes
    mesh%elementsets(elementset)%nodes = nodes(1:nnodes)

    deallocate ( nodes, sortnodes )

  end subroutine add_to_mesh_add_nodes_to_elementset


! Glue parts of the mesh together

  subroutine glue_mesh ( mesh, point1, point2, curve1, curve2 )

    type(mesh_t), intent(inout) :: mesh

!   points that are glued together
    integer, intent(in), optional :: point1, point2

!   curves that are glued together
    integer, intent(in), optional :: curve1, curve2


!   This routine glues part of the mesh together for use in the definition of
!   the array sidelem. Note that degrees of freedom are not afffected.


    integer :: npnt, pnt, np, nnodes, node
    integer, dimension(:), allocatable :: order
    integer, dimension(:,:), allocatable :: work


!   determine type of call

    if ( present(point1) ) then

!     glue points

      if ( point1 <=0 .or. point1 > mesh%npoints ) then
        write(*,'(/a/a,i0/)') &
          'Error: point1 in heading of glue_mesh has wrong value:', &
          'point1 <=0 or point1 > number of points = ', mesh%npoints
        stop
      end if

      if ( .not. present(point2) ) then
        write(*,'(/a/)') &
          'Error: point2 missing from heading of glue_mesh'
        stop
      else if ( point2 == point1 .or. &
                  point2 <=0 .or. point2 > mesh%npoints ) then
        write(*,'(/a/a,i0/)') &
          'Error: point2 in heading of glue_mesh has wrong value:', &
          'point2 == point1 or point2 <=0 or point2 > number of points = ', &
          mesh%npoints
        stop
      end if

      if ( allocated(mesh%gluepoints) ) then
!       glue points already exist
        npnt = mesh%ngluepoints+1
        allocate( work(npnt,2) )
        work(:mesh%ngluepoints,:) = mesh%gluepoints
        deallocate(mesh%gluepoints)
      else
!       glue points do not yet exist
        npnt = 1
        allocate( work(1,2) )
      end if

      work(npnt,1) = min(mesh%points(point1),mesh%points(point2))
      work(npnt,2) = max(mesh%points(point1),mesh%points(point2))

    else if ( present(curve1) ) then

!     glue curves

      if ( curve1 <=0 .or. curve1 > mesh%ncurves ) then
        write(*,'(/a/a,i0/)') &
          'Error: curve1 in heading of glue_mesh has wrong value:', &
          'curve1 <=0 or curve1 > number of curves = ', mesh%ncurves
        stop
      end if

      if ( .not. present(curve2) ) then
        write(*,'(/a/)') &
          'Error: curve2 missing from heading of glue_mesh'
        stop
      else if ( curve2 == curve1 .or. &
                  curve2 <=0 .or. curve2 > mesh%ncurves ) then
        write(*,'(/a/a,i0/)') &
          'Error: curve2 in heading of glue_mesh has wrong value:', &
          'curve2 == curve1 or curve2 <=0 or curve2 > number of curves = ', &
          mesh%ncurves
        stop
      else if ( mesh%curves(curve1)%nnodes /= mesh%curves(curve2)%nnodes ) then
        write(*,'(/a/a/)') &
          'Error: curve1 and curve2 in heading of glue_mesh have a different', &
          'number of nodes'
        stop
      end if

      nnodes = mesh%curves(curve1)%nnodes

      if ( allocated(mesh%gluepoints) ) then
!       glue points exist
        npnt = mesh%ngluepoints + nnodes
        allocate( work(npnt,2) )
        work(:mesh%ngluepoints,:) = mesh%gluepoints
        deallocate(mesh%gluepoints)
      else
!       glue points do not yet exist
        npnt = nnodes
        allocate( work(nnodes,2) )
      end if

      do node = 1, nnodes
        pnt = mesh%ngluepoints + node
        work(pnt,1) = min(mesh%curves(curve1)%nodes(node), &
                          mesh%curves(curve2)%nodes(node))
        work(pnt,2) = max(mesh%curves(curve1)%nodes(node), &
                          mesh%curves(curve2)%nodes(node))
      end do

    else

      write(*,'(/2a/4(/a)/)') &
        'Error: could not determine glue part ', &
        'from the heading of glue_mesh.', 'Non of the keywords:',      &
        '  point1 ',                                                   &
        '  curve1 ',                                                   &
        'is present.'
      stop

    end if

!   sort work array

    allocate(order(npnt))
    call sort( work(:,1), order )
    work(:,2) = work(order,2)
    deallocate(order)

!   delete double connections

    np = 0
    do pnt = 1, npnt
      if ( pnt < npnt ) then
        if ( work(pnt,1) == work(pnt+1,1) .and.  &
             work(pnt,2) == work(pnt+1,2) ) cycle ! double glue points
      end if
      np = np + 1
      work(np,:) = work(pnt,:) ! reuse work
    end do

!   create gluepoints

    allocate(mesh%gluepoints(np,2))
    mesh%gluepoints = work(:np,:)
    mesh%ngluepoints = np

    deallocate(work)

  end subroutine glue_mesh


! Define new blocks

  subroutine add_to_mesh_blocks ( mesh, blocks, blocksdomain )

    type(mesh_t), intent(inout) :: mesh

!   a new set of blocks is created. The size of the array must be equal to the
!   space dimension of the mesh. For example in 2D, specifying blocks=(/10,8/)
!   creates 80 new blocks in a 10x8 equidistant block structure.
!   The elements which have node 1 inside a particular block will be added to
!   that block. Only elements that do not yet belong yet to previously defined
!   blocks will be added.
!   If the parameter blocksdomain has not been specified the blocks span the
!   whole current mesh.
    integer, dimension(:), intent(in) :: blocks

!   if present the blocks are defined using the following domain:
!     blocksdomain(:,1) = `lower left corner'
!     blocksdomain(:,2) = `upper right corner'
!   For example:
!      blocksdomains = reshape ( (/ 1, 0, 0, 2, 1, 1 /), (/3,2/) )
!   defines the blocks between the coordinates (1,0,0) and (2,1,1).
    real(dp), dimension(:,:), intent(in), optional :: blocksdomain


    integer :: j, blk, firstblock, lastblock, blknr, elem, elgrp, nelem
    integer, dimension(mesh%ndim) :: ijk
    integer, allocatable, dimension(:,:) :: work
    real(dp), parameter :: eps = 1e-12_dp, epsmin = 1e-6_dp
    real(dp) :: dmax
    real(dp), dimension(mesh%ndim) :: delta, orig, xmin, xmax, x1


!   number of new blocks

    firstblock = mesh%nblocks + 1
    lastblock = mesh%nblocks + product(blocks)

    if ( lastblock > size(mesh%blocks) ) then
      write(*,'(/a/a,i0/)') &
        'Error: maximum exceeded in add_to_mesh:', &
        ' blocks > maximum = ', size(mesh%blocks)
      stop
    end if

!   define domain

    do j = 1, mesh%ndim
      xmin(j) = minval(mesh%coor(:,j))
      xmax(j) = maxval(mesh%coor(:,j))
    end do

    if ( present(blocksdomain) ) then

!     domain defined
      orig = blocksdomain(:,1)
      delta = ( blocksdomain(:,2) - blocksdomain(:,1) ) / blocks

    else

!     full mesh is covered
!     (take domain slightly larger to avoid truncation problems)
      orig = xmin - ( xmax - xmin ) * eps
      delta = ( 1 + 2 * eps ) * ( xmax - xmin ) / blocks

    end if

!   max size of the box
    dmax = maxval( abs(xmax - xmin) )

!   limit the edge size of the box to a small finite value to handle
!   zero sized edges like for a planar mesh in 3D.
    where ( abs(delta) < epsmin * dmax )
      orig = xmin - epsmin * dmax
      delta = 2 * epsmin * dmax / blocks
    end where

    allocate ( work(mesh%nelgrp,maxval(mesh%grpnumel)) )

    work = 0

!   indicate previously defined blocks

    do blk = 1, mesh%nblocks
      do elem = 1, mesh%blocks(blk)%nelem

        work(mesh%blocks(blk)%elements(elem,1), &
             mesh%blocks(blk)%elements(elem,2)) = -1

      end do
    end do

!   add new blocks

    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)

        if ( work(elgrp,elem) == -1 ) cycle ! element already in a block

        x1 = mesh%coor(mesh%topology(elgrp)%a(1,elem),:) ! first node of element
        ijk = floor ( ( x1 - orig ) / delta ) ! i, j, k of block

        if ( any ( ijk < 0 ) .or. any ( ijk > blocks-1 ) ) cycle !outside domain

!       fill work

        blknr = 1 + ijk(1)
        if ( mesh%ndim >= 2 ) blknr = blknr + ijk(2) * blocks(1)
        if ( mesh%ndim >= 3 ) blknr = blknr + ijk(3) * blocks(1) * blocks(2)
        work(elgrp,elem) = blknr + mesh%nblocks

      end do
    end do

!   fill blocks in block_t structure of mesh

    do blk = firstblock, lastblock

      nelem = count ( work == blk )

      allocate ( mesh%blocks(blk)%bounds(mesh%ndim,2) )
      allocate ( mesh%blocks(blk)%elbounds(nelem,mesh%ndim,2) )
      allocate ( mesh%blocks(blk)%elements(nelem,2) )

    end do

    mesh%blocks(firstblock:lastblock)%nelem = 0

    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)

        blknr = work(elgrp,elem)

        if ( blknr <= 0 ) cycle ! element already in a block (-1) or
                                ! element not found (0)

        mesh%blocks(blknr)%nelem = mesh%blocks(blknr)%nelem + 1
        mesh%blocks(blknr)%elements(mesh%blocks(blknr)%nelem,1) = elgrp
        mesh%blocks(blknr)%elements(mesh%blocks(blknr)%nelem,2) = elem

      end do
    end do

    deallocate ( work )

!   fill bounds

    call find_bounds_blocks ( mesh, firstblock, lastblock )

    mesh%nblocks = lastblock

  end subroutine add_to_mesh_blocks


! fill bounds of blocks

  subroutine find_bounds_blocks ( mesh, firstblock, lastblock )

    use limits_m, only: EPSBOUND

    type(mesh_t), intent(inout) :: mesh

!   if present the bounds in blocks firstblock:lastblock will be computed
!   otherwise all blocks will be used.
    integer, intent(in), optional :: firstblock, lastblock


    integer :: fblock, lblock, j, blk, elem, elgrp, elm
    real(dp), dimension(mesh%ndim) :: xmin, xmax


!   which blocks?

    if ( present(firstblock) .and. present(lastblock) ) then
!     specified range of blocks only
      fblock = firstblock
      lblock = lastblock
    else
!     all blocks
      fblock = 1
      lblock = mesh%nblocks
    end if


    do blk = fblock, lblock

!     fill elbounds

      do elm = 1, mesh%blocks(blk)%nelem

        elgrp = mesh%blocks(blk)%elements(elm,1)
        elem  = mesh%blocks(blk)%elements(elm,2)

        do j = 1, mesh%ndim
          xmin(j) = minval(mesh%coor(mesh%topology(elgrp)%a(:,elem),j))
          xmax(j) = maxval(mesh%coor(mesh%topology(elgrp)%a(:,elem),j))
        end do

        mesh%blocks(blk)%elbounds(elm,:,1) = xmin - EPSBOUND * ( xmax - xmin )
        mesh%blocks(blk)%elbounds(elm,:,2) = xmax + EPSBOUND * ( xmax - xmin )

      end do

!     fill bounds (from elbounds)

      do j = 1, mesh%ndim
        mesh%blocks(blk)%bounds(j,1) = minval(mesh%blocks(blk)%elbounds(:,j,1))
        mesh%blocks(blk)%bounds(j,2) = maxval(mesh%blocks(blk)%elbounds(:,j,2))
      end do

    end do

  end subroutine find_bounds_blocks


! Define structured blocks

  subroutine add_to_mesh_sblocks ( mesh, sblocks, blocksdomain )

    use limits_m, only: EPSBOUND

    type(mesh_t), intent(inout) :: mesh

!   a set of structured blocks is created.
!   The size of the array must be equal to the space dimension of the mesh.
!   A previous definition of structured blocks will be removed. For example in
!   2D, specifying sblocks=(/10,8/) creates 80 structured blocks in a
!   10x8 equidistant block structure.
!   The elements which are (likely to be) partly or fully inside a particular
!   block will be added to that block. Elements can be part of multiple blocks.
!   If the parameter blocksdomain has not been specified the blocks span the
!   whole current mesh.
    integer, dimension(:), intent(in) :: sblocks

!   if present the structured blocks are defined using the following domain:
!     blocksdomain(:,1) = `lower left corner'
!     blocksdomain(:,2) = `upper right corner'
!   For example:
!      blocksdomains = reshape ( (/ 1, 0, 0, 2, 1, 1 /), (/3,2/) )
!   defines the blocks between the coordinates (1,0,0) and (2,1,1).
    real(dp), dimension(:,:), intent(in), optional :: blocksdomain


    integer :: i, j, k, elem, elgrp, nb, blk
    integer, dimension(mesh%ndim) :: ijkmin, ijkmax
    real(dp), parameter :: eps = 1e-12_dp, epsmin = 1e-6_dp
    real(dp) :: dmax
    real(dp), dimension(mesh%ndim) :: delta, orig, xmin, xmax, xmine, xmaxe
    integer, allocatable, dimension(:) :: work


    if ( mesh%sblocks%filled ) then
!     delete previous defined structured blocks
      call delete_sblocks ( mesh%sblocks )
    end if

!   define domain

    do j = 1, mesh%ndim
      xmin(j) = minval(mesh%coor(:,j))
      xmax(j) = maxval(mesh%coor(:,j))
    end do

    if ( present(blocksdomain) ) then

!     domain defined
      orig = blocksdomain(:,1)
      delta = ( blocksdomain(:,2) - blocksdomain(:,1) ) / sblocks

    else

!     full mesh is covered
!     (take domain slightly larger to avoid truncation problems)
      orig = xmin - ( xmax - xmin ) * eps
      delta = ( 1 + 2 * eps ) * ( xmax - xmin ) / sblocks

    end if

!   max size of the box
    dmax = maxval( abs(xmax - xmin) )

!   limit the edge size of the box to a small finite value to handle
!   zero sized edges like for a planar mesh in 3D.
    where ( abs(delta) < epsmin * dmax )
      orig = xmin - epsmin * dmax
      delta = 2 * epsmin * dmax / sblocks
    end where

!   create structure sblocks

    allocate ( mesh%sblocks%nbl(mesh%ndim) )
    allocate ( mesh%sblocks%x0(mesh%ndim) )
    allocate ( mesh%sblocks%dx(mesh%ndim) )

    mesh%sblocks%nbl = sblocks
    mesh%sblocks%x0 = orig
    mesh%sblocks%dx = delta

    mesh%sblocks%nsblocks = product ( sblocks )

    allocate ( mesh%sblocks%sblks(mesh%sblocks%nsblocks) )
    allocate ( work(mesh%sblocks%nsblocks) )

!   first pass: count number of elements in each structured block

    mesh%sblocks%sblks(:)%nelem = 0  ! initialize nelem to zero

    do elgrp = 1, mesh%nelgrp
 e1:  do elem = 1, mesh%grpnumel(elgrp)

!       compute xmine, xmaxe, ijkmin, ijkmax

        call first_part

        do j = 1, mesh%ndim
          if ( ( ijkmin(j) < 0 .and. ijkmax(j) < 0 ) .or. &
               ( ijkmin(j) > sblocks(j)-1 .and. ijkmax(j) > sblocks(j)-1 ) ) &
          then
!           element fully outside sblocks domain
            cycle e1
          end if
        end do

!       compute nb and work subscript array

        call second_part

        mesh%sblocks%sblks(work(1:nb))%nelem = &
                             mesh%sblocks%sblks(work(1:nb))%nelem + 1

      end do e1
    end do

!   allocate sblocks

    do blk = 1, mesh%sblocks%nsblocks
      allocate ( mesh%sblocks%sblks(blk)%elbounds(&
                            &mesh%sblocks%sblks(blk)%nelem,mesh%ndim,2) )
      allocate ( mesh%sblocks%sblks(blk)%elements(&
                            &mesh%sblocks%sblks(blk)%nelem,2) )
    end do

!   second pass: store elements in each structured block

    mesh%sblocks%sblks(:)%nelem = 0  ! initialize nelem to zero again

    do elgrp = 1, mesh%nelgrp
 e2:  do elem = 1, mesh%grpnumel(elgrp)

!       compute xmine, xmaxe, ijkmin, ijkmax

        call first_part

        do j = 1, mesh%ndim
          if ( ( ijkmin(j) < 0 .and. ijkmax(j) < 0 ) .or. &
               ( ijkmin(j) > sblocks(j)-1 .and. ijkmax(j) > sblocks(j)-1 ) )&
          then
!           element fully outside sblocks domain
            cycle e2
          end if
        end do

!       compute nb and work subscript array

        call second_part

        mesh%sblocks%sblks(work(1:nb))%nelem = &
                             mesh%sblocks%sblks(work(1:nb))%nelem + 1

        do j = 1, nb
          mesh%sblocks%sblks(work(j))%elements(&
                         mesh%sblocks%sblks(work(j))%nelem,1) = elgrp
          mesh%sblocks%sblks(work(j))%elements(&
                         mesh%sblocks%sblks(work(j))%nelem,2) = elem
          mesh%sblocks%sblks(work(j))%elbounds(&
                         mesh%sblocks%sblks(work(j))%nelem,:,1) = xmin
          mesh%sblocks%sblks(work(j))%elbounds(&
                         mesh%sblocks%sblks(work(j))%nelem,:,2) = xmax
        end do

      end do e2
    end do

    mesh%sblocks%filled = .true.

    deallocate ( work )

  contains

!   first part in loop

    subroutine first_part

      integer :: j

!     element domain (outline)

      do j = 1, mesh%ndim
        xmine(j) = minval(mesh%coor(mesh%topology(elgrp)%a(:,elem),j))
        xmaxe(j) = maxval(mesh%coor(mesh%topology(elgrp)%a(:,elem),j))
      end do

!     make box around element somewhat larger
      xmine = xmine - EPSBOUND * ( xmaxe - xmine )
      xmaxe = xmaxe + EPSBOUND * ( xmaxe - xmine )

      ijkmin = floor ( ( xmine - orig ) / delta ) ! i, j, k of xmin
      ijkmax = floor ( ( xmaxe - orig ) / delta ) ! i, j, k of xmax

    end subroutine first_part

!   second part in loop

    subroutine second_part

!     only real blocks included

      ijkmin = max ( ijkmin, 0 )
      ijkmin = min ( ijkmin, sblocks - 1 )
      ijkmax = max ( ijkmax, 0 )
      ijkmax = min ( ijkmax, sblocks - 1 )

!     add element to all blocks covered between lower left (xmin)
!     and upper right (xmax)

      nb = product( ijkmax - ijkmin + 1 )

      select case ( mesh%ndim )
      case(1)
        work(1:nb) = [ (1+i,i=ijkmin(1),ijkmax(1)) ]
      case(2)
        work(1:nb) = [ ((1+i+j*sblocks(1),i=ijkmin(1),ijkmax(1)),j= &
                     ijkmin(2),ijkmax(2)) ]
      case(3)
        work(1:nb) = [ (((1+i+j*sblocks(1)+k*sblocks(1)*sblocks(2),&
                &i=ijkmin(1),ijkmax(1)),j=ijkmin(2),ijkmax(2)),&
                &k=ijkmin(3),ijkmax(3)) ]
      case default
        call errormsg_case_default ( 'second_part', 'mesh%ndim', &
          int_value=mesh%ndim )
      end select

    end subroutine second_part

  end subroutine add_to_mesh_sblocks


! Define nodal blocks

  subroutine add_to_mesh_nodblocks ( mesh, nodblocks, blocksdomain )

    use limits_m, only: EPSNODEBOX

    type(mesh_t), intent(inout) :: mesh

!   a set of nodal blocks (nodblocks) is created.
!   The size of the array must be equal to the space dimension of the mesh.
!   A previous definition of nodblocks will be removed. For example in
!   2D, specifying nodblocks=(/10,8/) creates 80 nodal blocks in a
!   10x8 equidistant block structure.
!   Around each node a small box is created with the node at the center of the
!   box. The nodes where its box is partly or fully inside a particular
!   block will be added to that block. Nodes can be part of multiple blocks.
!   If the parameter blocksdomain has not been specified the blocks span the
!   whole current mesh.
    integer, dimension(:), intent(in) :: nodblocks

!   if present the blocks, sblocks or nodblocks are defined using the following
!   domain:
!     blocksdomain(:,1) = `lower left corner'
!     blocksdomain(:,2) = `upper right corner'
!   For example:
!      blocksdomains = reshape ( (/ 1, 0, 0, 2, 1, 1 /), (/3,2/) )
!   defines the blocks between the coordinates (1,0,0) and (2,1,1).
    real(dp), dimension(:,:), intent(in), optional :: blocksdomain


    integer :: i, j, k, node, nb, blk
    integer, dimension(mesh%ndim) :: ijkmin, ijkmax
    real(dp), parameter :: eps = 1e-12_dp, epsmin = 1e-6_dp
    real(dp) :: dmax
    real(dp), dimension(mesh%ndim) :: delta, orig, xmin, xmax, xminb, xmaxb
    integer, allocatable, dimension(:) :: work


    if ( mesh%nodblocks%filled ) then
!     delete previous defined nodal blocks
      call delete_nodblocks ( mesh%nodblocks )
    end if

!   define domain

    do j = 1, mesh%ndim
      xmin(j) = minval(mesh%coor(:,j))
      xmax(j) = maxval(mesh%coor(:,j))
    end do

    if ( present(blocksdomain) ) then

!     domain defined
      orig = blocksdomain(:,1)
      delta = ( blocksdomain(:,2) - blocksdomain(:,1) ) / nodblocks

    else

!     full mesh is covered
!     (take domain slightly larger to avoid truncation problems)
      orig = xmin - ( xmax - xmin ) * eps
      delta = ( 1 + 2 * eps ) * ( xmax - xmin ) / nodblocks

    end if

!   max size of the box
    dmax = maxval( abs(xmax - xmin) )

!   limit the edge size of the box to a small finite value to handle
!   zero sized edges like for a planar mesh in 3D.
    where ( abs(delta) < epsmin * dmax )
      orig = xmin - epsmin * dmax
      delta = 2 * epsmin * dmax / nodblocks
    end where

!   create nodal nodblocks

    allocate ( mesh%nodblocks%nbl(mesh%ndim) )
    allocate ( mesh%nodblocks%x0(mesh%ndim) )
    allocate ( mesh%nodblocks%dx(mesh%ndim) )

    mesh%nodblocks%nbl = nodblocks
    mesh%nodblocks%x0 = orig
    mesh%nodblocks%dx = delta

    mesh%nodblocks%nnodblocks = product ( nodblocks )

    allocate ( mesh%nodblocks%nodblks(mesh%nodblocks%nnodblocks) )
    allocate ( work(mesh%nodblocks%nnodblocks) )

!   first pass: count number of nodes in each structured block

    mesh%nodblocks%nodblks(:)%nnodes = 0  ! initialize nnodes to zero

n1: do node = 1, mesh%nnodes

!       compute xminb, xmaxb, ijkmin, ijkmax

        call first_part

        do j = 1, mesh%ndim
          if ( ( ijkmin(j) < 0 .and. ijkmax(j) < 0 ) .or. &
           ( ijkmin(j) > nodblocks(j)-1 .and. ijkmax(j) > nodblocks(j)-1 ) )&
          then
!           node fully outside nodblocks domain
            cycle n1
          end if
        end do

!       compute nb and work subscript array

        call second_part

        mesh%nodblocks%nodblks(work(1:nb))%nnodes = &
                     mesh%nodblocks%nodblks(work(1:nb))%nnodes + 1

    end do n1

!   allocate nodblocks

    do blk = 1, mesh%nodblocks%nnodblocks
      allocate ( mesh%nodblocks%nodblks(blk)%nodes(&
                                &mesh%nodblocks%nodblks(blk)%nnodes) )
    end do

!   second pass: store nodes in each nodal block

    mesh%nodblocks%nodblks(:)%nnodes = 0  ! initialize nnodes to zero

n2: do node = 1, mesh%nnodes

!       compute xminb, xmaxb, ijkmin, ijkmax

        call first_part

        do j = 1, mesh%ndim
          if ( ( ijkmin(j) < 0 .and. ijkmax(j) < 0 ) .or. &
            ( ijkmin(j) > nodblocks(j)-1 .and. ijkmax(j) > nodblocks(j)-1 ) ) &
          then
!           element fully outside sblocks domain
            cycle n2
          end if
        end do

!       compute nb and work subscript array

        call second_part

        mesh%nodblocks%nodblks(work(1:nb))%nnodes = &
                             mesh%nodblocks%nodblks(work(1:nb))%nnodes + 1

        do j = 1, nb
          mesh%nodblocks%nodblks(work(j))%nodes(&
                         &mesh%nodblocks%nodblks(work(j))%nnodes) = node
        end do

    end do n2

    mesh%nodblocks%filled = .true.

    deallocate ( work )

  contains

!   first part in loop

    subroutine first_part

!     nodal box domain

      xminb = mesh%coor(node,:) - ( xmax - xmin ) * EPSNODEBOX
      xmaxb = mesh%coor(node,:) + ( xmax - xmin ) * EPSNODEBOX

      ijkmin = floor ( ( xminb - orig ) / delta ) ! i, j, k of xmin
      ijkmax = floor ( ( xmaxb - orig ) / delta ) ! i, j, k of xmax

    end subroutine first_part

!   second part in loop

    subroutine second_part

!     only real blocks included

      ijkmin = max ( ijkmin, 0 )
      ijkmin = min ( ijkmin, nodblocks - 1 )
      ijkmax = max ( ijkmax, 0 )
      ijkmax = min ( ijkmax, nodblocks - 1 )

!     add node to all blocks covered between lower left (xmin)
!     and upper right (xmax)

      nb = product( ijkmax - ijkmin + 1 )

      select case ( mesh%ndim )
      case(1)
        work(1:nb) = [ (1+i,i=ijkmin(1),ijkmax(1)) ]
      case(2)
        work(1:nb) = [ ((1+i+j*nodblocks(1),i=ijkmin(1),ijkmax(1)),&
                                &j=ijkmin(2),ijkmax(2)) ]
      case(3)
        work(1:nb) = [ (((1+i+j*nodblocks(1)+k*nodblocks(1)*nodblocks(2),&
                &i=ijkmin(1),ijkmax(1)),j=ijkmin(2),ijkmax(2)),&
                &k=ijkmin(3),ijkmax(3)) ]
      case default
        call errormsg_case_default ( 'second_part', 'mesh%ndim', &
          int_value=mesh%ndim )
      end select

    end subroutine second_part

  end subroutine add_to_mesh_nodblocks


! local renumbering for the nodes in nodes2 to match the coordinates of the
! nodes in nodes1

  subroutine renumber_nodes_to_match_coor ( mesh, nodes1, nodes2, perm, mesh2, &
    displace, mindistance )

    use limits_m, only: EPSNODEOVERLAP

    type(mesh_t), intent(in) :: mesh

!   the nodes in nodes2 need to be renumbered in order to match the coordinates
!   of nodes1
    integer, dimension(:), intent(in) :: nodes1, nodes2

!   the renumbering of nodes2 to match the coordinates of the nodes in nodes1
!   perm(i) is the new local node number for the node in nodes2(i)
!   NOTE: if perm(i) == 0, no matching node has been found.
    integer, dimension(:), intent(out) :: perm

!   if present the nodes in nodes2 are in mesh2 instead of mesh.
    type(mesh_t), intent(in), optional :: mesh2

!   if present the coordinates of the nodes in nodes2 are "displaced" virtually
!   with vector displace to match the coordinates of nodes1
!   NOTE: the coordinates of the nodes are not actually changed.
    real(dp), dimension(:), intent(in), optional :: displace

!   if present the matching of coordinates is based on the minimum distance
!   default=.false.
    logical, intent(in), optional :: mindistance

    integer :: nnodes, node1, node2, nodenr1, nodenr2
    real(dp), dimension(mesh%ndim) :: domainsize, disp
    real(dp), dimension(:), allocatable :: distance2
    logical :: lmindistance

    if ( present(displace) ) then
      disp = displace
    else
      disp = 0
    end if

    lmindistance = set_optional ( variable=mindistance, default=.false. )

    nnodes = size(nodes1)

    perm = 0 ! initialize

    domainsize = maxval(mesh%coor,dim=1) - minval(mesh%coor,dim=1)

!   if the domainsize is very small (~zero) in a direction, a matching node
!   will never be found since the distance should be smaller then
!   EPSNODEOVERLAP*domainsize. In that case, set the domain size to 1

    where ( domainsize < EPSNODEOVERLAP ) domainsize = 1

    if ( lmindistance ) then

!     based on the minimum distance

      allocate ( distance2(nnodes) )

      if ( present(mesh2) ) then

!       second mesh

        do node2 = 1, nnodes

          nodenr2 = nodes2(node2)

!         compute square distance

          do node1 = 1, nnodes

            nodenr1 = nodes1(node1)

            distance2(node1) = &
               dot_product ( mesh2%coor(nodenr2,:)+disp-mesh%coor(nodenr1,:), &
                             mesh2%coor(nodenr2,:)+disp-mesh%coor(nodenr1,:) )

          end do

!         matching node
          perm(node2) = minval ( minloc ( distance2 ) )

        end do

      else

!       one mesh

        do node2 = 1, nnodes

          nodenr2 = nodes2(node2)

!         compute square distance

          do node1 = 1, nnodes

            nodenr1 = nodes1(node1)

            distance2(node1) = &
               dot_product ( mesh%coor(nodenr2,:)+disp-mesh%coor(nodenr1,:), &
                             mesh%coor(nodenr2,:)+disp-mesh%coor(nodenr1,:) )

          end do

!         matching node
          perm(node2) = minval ( minloc ( distance2 ) )

        end do

      end if

      deallocate ( distance2 )

    else

!     based on zero distance to a node

      if ( present(mesh2) ) then

!       second mesh

        do node2 = 1, nnodes

          nodenr2 = nodes2(node2)

          do node1 = 1, nnodes

            nodenr1 = nodes1(node1)

            if ( all ( abs(mesh2%coor(nodenr2,:)+disp-mesh%coor(nodenr1,:) ) &
                                       < EPSNODEOVERLAP*domainsize ) ) then
!             matching node
              perm(node2) = node1

              exit ! leave loop

            end if

          end do

        end do

      else

!       one mesh

        do node2 = 1, nnodes

          nodenr2 = nodes2(node2)

          do node1 = 1, nnodes

            nodenr1 = nodes1(node1)

            if ( all ( abs(mesh%coor(nodenr2,:)+disp-mesh%coor(nodenr1,:) ) &
                                        < EPSNODEOVERLAP*domainsize ) ) then
!             matching node
              perm(node2) = node1

              exit ! leave loop

            end if

          end do

        end do

      end if

    end if

  end subroutine renumber_nodes_to_match_coor


! local renumbering for the nodes in nodes2 to match the coordinates of the
! nodes in nodes1 after a rotation and/or displacement

  subroutine renumber_nodes_to_match_coor_rotate ( mesh, nodes1, nodes2, perm, &
     mesh2, displace, rotate, x0, mindistance )

    use limits_m, only: EPSNODEOVERLAP

    type(mesh_t), intent(in) :: mesh

!   the nodes in nodes2 need to be renumbered in order to match the coordinates
!   of nodes1
    integer, dimension(:), intent(in) :: nodes1, nodes2

!   the renumbering of nodes2 to match the coordinates of the nodes in nodes1
!   perm(i) is the new local node number for the node in nodes2(i)
!   NOTE: if perm(i) == 0, no matching node has been found.
    integer, dimension(:), intent(out) :: perm

!   if present the nodes in nodes2 are in mesh2 instead of mesh.
    type(mesh_t), intent(in), optional :: mesh2

!   if present the coordinates of the nodes in nodes2 are "displaced" virtually
!   with vector displace to match the coordinates of nodes1
!   NOTE: the coordinates of the nodes are not actually changed.
    real(dp), dimension(:), intent(in), optional :: displace

!   if present the coordinates of the nodes in surface(2) are "rotated"
!   virtually through the angle specified in rotate for matching the coordinates
!   of surface(1).
!   for 2D meshes:
!   rotate(1): angle of rotation in the xy-plane
!   for 3D meshes:
!   rotate(1): angle of rotation around x-axis
!   rotate(2): angle of rotation around y-axis
!   rotate(3): angle of rotation around z-axis
!   NOTE: the coordinates of the nodes are not actually changed.
!   NOTE2: in 3D rotations are applied in an order of increasing dimension, so
!   first around x-axis, then y-axis and finally z-axis
!   NOTE3: if displace is also present, this is performed after the rotation
    real(dp), dimension(:), intent(in) :: rotate

!   if present perform the rotation w.r.t to point x0, otherwise rotate around
!   origin
    real(dp), dimension(:), intent(in), optional :: x0

!   if present the matching of coordinates is based on the minimum distance
!   default=.false.
    logical, intent(in), optional :: mindistance

    integer :: nnodes, node1, node2, nodenr1, nodenr2
    real(dp) :: theta
    real(dp), dimension(mesh%ndim) :: domainsize, disp, coor_new, lx0
    real(dp), dimension(mesh%ndim,mesh%ndim,mesh%ndim) :: rot_mat
    real(dp), dimension(:), allocatable :: distance2
    logical :: lmindistance

    if ( present(displace) ) then
      disp = displace
    else
      disp = 0
    end if

    if ( present(x0) ) then
      lx0 = x0
    else
      lx0 = 0
    end if

    if ( mesh%ndim == 2 ) then
!     rotate in the xy-place
      theta = rotate(1)
      rot_mat(1,:,:) = reshape ( [  cos(theta), sin(theta), &
                                   -sin(theta), cos(theta) ], [2,2] )
    else
!     rotation around x
      theta = rotate(1)
      rot_mat(1,:,:) = reshape ( [1._dp,       0._dp,      0._dp, &
                                  0._dp,  cos(theta), sin(theta), &
                                  0._dp, -sin(theta), cos(theta) ], [3,3] )

!     rotation around y
      theta = rotate(2)
      rot_mat(2,:,:) = reshape ( [ cos(theta), 0._dp, -sin(theta), &
                                        0._dp, 1._dp,       0._dp, &
                                   sin(theta), 0._dp,  cos(theta) ], [3,3] )

!     rotation around z
      theta = rotate(3)
      rot_mat(3,:,:) = reshape ( [ cos(theta), sin(theta), 0._dp, &
                                  -sin(theta), cos(theta), 0._dp, &
                                        0._dp,     0.0_dp, 1._dp ], [3,3] )
    end if

    lmindistance = set_optional ( variable=mindistance, default=.false. )

    nnodes = size(nodes1)

    perm = 0 ! initialize

    domainsize = maxval(mesh%coor,dim=1) - minval(mesh%coor,dim=1)

!   if the domainsize is very small (~zero) in a direction, a matching node
!   will never be found since the distance should be smaller then
!   EPSNODEOVERLAP*domainsize. In that case, set the domain size to 1

    where ( domainsize < EPSNODEOVERLAP ) domainsize = 1

    if ( lmindistance ) then

!     based on the minimum distance

      allocate ( distance2(nnodes) )

      if ( present(mesh2) ) then

!       second mesh

        do node2 = 1, nnodes

          nodenr2 = nodes2(node2)

!         compute square distance

          do node1 = 1, nnodes

            nodenr1 = nodes1(node1)

            coor_new = lx0 + matmul ( rot_mat(1,:,:), mesh2%coor(nodenr2,:)-lx0)
            if ( mesh%ndim == 3 ) then
              coor_new = lx0 + matmul ( rot_mat(2,:,:), coor_new-lx0 )
              coor_new = lx0 + matmul ( rot_mat(3,:,:), coor_new-lx0 )
            end if
            coor_new = coor_new + disp

            distance2(node1) = &
               dot_product ( coor_new-mesh%coor(nodenr1,:), &
                             coor_new-mesh%coor(nodenr1,:) )

          end do

!         matching node
          perm(node2) = minval ( minloc ( distance2 ) )

        end do

      else

!       one mesh

        do node2 = 1, nnodes

          nodenr2 = nodes2(node2)

!         compute square distance

          do node1 = 1, nnodes

            nodenr1 = nodes1(node1)

            coor_new = lx0 + matmul ( rot_mat(1,:,:), mesh%coor(nodenr2,:)-lx0)
            if ( mesh%ndim == 3 ) then
              coor_new = lx0 + matmul ( rot_mat(2,:,:), coor_new-lx0 )
              coor_new = lx0 + matmul ( rot_mat(3,:,:), coor_new-lx0 )
            end if
            coor_new = coor_new + disp

            distance2(node1) = &
               dot_product ( coor_new-mesh%coor(nodenr1,:), &
                             coor_new-mesh%coor(nodenr1,:) )

          end do

!         matching node
          perm(node2) = minval ( minloc ( distance2 ) )

        end do

      end if

      deallocate ( distance2 )

    else

!     based on zero distance to a node

      if ( present(mesh2) ) then

!       second mesh

        do node2 = 1, nnodes

          nodenr2 = nodes2(node2)

          do node1 = 1, nnodes

            nodenr1 = nodes1(node1)

            coor_new = lx0 + matmul ( rot_mat(1,:,:), mesh2%coor(nodenr2,:)-lx0)
            if ( mesh%ndim == 3 ) then
              coor_new = lx0 + matmul ( rot_mat(2,:,:), coor_new-lx0 )
              coor_new = lx0 + matmul ( rot_mat(3,:,:), coor_new-lx0 )
            end if
            coor_new = coor_new + disp

            if ( all ( abs(coor_new-mesh%coor(nodenr1,:) ) &
                                       < EPSNODEOVERLAP*domainsize ) ) then
!             matching node
              perm(node2) = node1

              exit ! leave loop

            end if

          end do

        end do

      else

!       one mesh

        do node2 = 1, nnodes

          nodenr2 = nodes2(node2)

          do node1 = 1, nnodes

            nodenr1 = nodes1(node1)

            coor_new = lx0 + matmul ( rot_mat(1,:,:), mesh%coor(nodenr2,:)-lx0)
            if ( mesh%ndim == 3 ) then
              coor_new = lx0 + matmul ( rot_mat(2,:,:), coor_new-lx0 )
              coor_new = lx0 + matmul ( rot_mat(3,:,:), coor_new-lx0 )
            end if
            coor_new = coor_new + disp

            if ( all ( abs(coor_new-mesh%coor(nodenr1,:) ) &
                                        < EPSNODEOVERLAP*domainsize ) ) then
!             matching node
              perm(node2) = node1

              exit ! leave loop

            end if

          end do

        end do

      end if

    end if

  end subroutine renumber_nodes_to_match_coor_rotate


end module meshgen_construct_m
