
! Copyright (C) 2004-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines mesh type and associated routines

module mesh_m

  use glob_defs_m
  use array_defs_m
  use set_optional_m

  implicit none

! type definition of an element

  type element_t

!   global shape of the element:
!     line
!     quadrilateral
!     triangle
!     hexahedron
!     tetrahedron
!     prism
!     pyramid
    character (len=13) :: globalshape = ''

!   shape of the element
!      elshape=1 : 2 node line element
!      elshape=2 : 3 node line element
!      elshape=3 : 3 node triangle
!      elshape=4 : 6 node triangle
!      elshape=5 : 4 node quadrilateral
!      elshape=6 : 9 node quadrilateral
!      elshape=7 : 7 node triangle
!      elshape=9 : 5 node quadrilateral
!      elshape=10: 4 node triangle
!      elshape=11: 4 node tetrahedron
!      elshape=12: 10 node tetrahedron
!      elshape=13: 8 node hexahedron
!      elshape=14: 27 node hexahedron
!      elshape=15: 14 node tetrahedron
!      elshape=16: 15 node tetrahedron
!      elshape=17: 9 node hexahedron
!      elshape=18: 5 node tetrahedron
!      elshape=30: 8 node quadrilateral
!      elshape=31: 20 node hexahedron
!      elshape=32: 3 node macro line element (2 linear subelements)
!      elshape=33: 6 node macro triangle (4 linear subelements)
!      elshape=34: 9 node macro quadrilateral (4 bilinear subelements)
!      elshape=35: 10 node macro tetrahedron (8 linear subelements)
!      elshape=36: 27 node macro hexahedron (8 trilinear subelements)
!      elshape=41: 6 node prism
!      elshape=42: 15 node prism
!      elshape=43: 18 node prism
!      elshape=51: 5 node pyramid
!      elshape=52: 13 node pyramid
!      elshape=53: 14 node pyramid
!      elshape=101: line high-order element
!      elshape=102: quadrilateral high-order element
!      elshape=103: hexahedron high-order element
!      elshape=104: triangle high-order element
!      elshape=105: line macro high-order element (linear subelements)
!      elshape=106: quadrilateral macro high-order element (bilinear subelms)
!      elshape=107: hexahedron macro high-order element (trilinear subelements)

    integer :: elshape  = 0

    integer :: ndim     = 0     ! space dimension
    integer :: ndimr    = 0     ! dimension of reference space
    integer :: numnod   = 0     ! number of nodal points
    integer :: numshapenod = 0  ! number of `shape' nodal points for mapping:
                                ! it is different from numnod for
                                ! elshape=7,9,10, having a value of numnod-1
                                ! because the center node is missing.
                                ! elshape=15,16, having a value of numnod=10
                                ! like the 10-node tetrahedron

!   sides: the nodes, edges or faces that form the usual connection to other
!          elements. Line elements connect using end nodes, quads and triangles
!          connect using edges and 3D elements (tet, hex) using faces.
    integer :: numsides = 0     ! number of sides
    integer :: sidnumvert = 0   ! number of vertices of a side
    integer :: sidnumnod = 0    ! number of nodes of a side
    integer :: numinternnod = 0 ! number of internal nodes
    integer :: sidnuminternnod = 0 ! number of internal nodes on a side

!   edges: the one-dimensional "boundaries". Curves basically consist of a
!          collection of edges.
    integer :: numedges = 0     ! number of edges
    integer :: edgenumnod = 0   ! number of nodes of an edge

!   faces: the two-dimensional "boundaries". Surfaces basically consist of a
!          collection of faces. Not available for line elements.
    integer :: numfaces = 0     ! number of faces
    integer :: facenumvert = 0  ! number of vertices of a face
    integer :: facenumnod = 0   ! number of nodes of a face

!   high-order elements:
!      p(i,1) : polynomial order in the direction i  (note: p>=1)
!               Note: variable polynomial order only possible for
!                     elshape=101,102,103 (or 105,106,107)
!      p(i,2) : layout of the nodes along the direction i:
!                 0: equidistant
!                 1: Gauss-Lobatto
!               Note: p(i,2)=1 only possible for
!                     elshape=101,102,103 (or 105,106,107)
    integer, dimension(3,2) :: p=0

    integer :: numsubdiv = 0    ! number of sub divisions for postprocessing
                                ! with figplot

    integer :: vtk_numnod  = 0  ! number of nodal points for VTK output
    integer :: vtk_numelm  = 0  ! number of subelements for VTK output
    integer :: vtk_elshape = 0  ! shape number for the VTK output

    integer :: tec_numnod  = 0  ! number of nodal points for Tecplot output
    integer :: tec_numelm  = 0  ! number of subelements for Tecplot output
    character (len=15) :: tec_elshape = ''  ! shape for Tecplot output

    integer :: gmsh_numnod  = 0  ! number of nodal points for gmsh output
    integer :: gmsh_elshape = 0  ! shape number for the gmsh output

    integer :: gmsh_numnod_parsed  = 0  ! number of nodal points for gmsh parsed
    integer :: gmsh_numelm_parsed = 0  ! number of subelements for gmsh parsed
    character (len=1) :: gmsh_elshape_parsed = ''  ! shape for gmsh parsed

!   vertices on the sides
!   sidvert(vert,side) gives the local node number of vertex vert on sidenr side
!   for 3D elements direction is according to outside normal
    integer, allocatable, dimension(:,:) :: sidvert

!   nodes on the sides
!   sidnod(node,side) gives the local node number of nodenr node on sidenr side
!   for 3D elements direction is according to outside normal
    integer, allocatable, dimension(:,:) :: sidnod

!   internal nodes
!   internnod(:) gives the local node numbers of internal nodes
    integer, allocatable, dimension(:) :: internnod

!   internal nodes on a side
!   sidinternnod(node,side) gives the local node numbers of internal nodes
!   on a side
    integer, allocatable, dimension(:,:) :: sidinternnod

!   vertices of the edges
!   edgevert(vert,edge) gives the local node number of vertex vert on the edge
    integer, allocatable, dimension(:,:) :: edgevert

!   nodes on the edges
!   edgenod(node,edge) gives the local node number of nodenr node on the edge
    integer, allocatable, dimension(:,:) :: edgenod

!   vertices on the faces
!   facevert(vert,face) gives the local node number of vertex vert on the face
!   In 3D the direction is always according to the outside normal of the element
    integer, allocatable, dimension(:,:) :: facevert

!   nodes on the faces
!   facenod(node,face) gives the local node number of nodenr node on the face
!   In 3D the direction is always according to the outside normal of the element
    integer, allocatable, dimension(:,:) :: facenod

!   subdivision for postprocessing with figplot:
!     ndim=2 subdivide into triangles
!   subtopol(node,elem) gives the local node numbers of the internal element
!   elem. Note: 1<=node<=3 (ndim=2), 1<=elem<=numsubdiv
    integer, allocatable, dimension(:,:) :: subtopol

!   vtk_index(node,elem) gives the local node numbers for the VTK format output
    integer, allocatable, dimension(:,:) :: vtk_index

!   tec_index(node,elem) gives the local node numbers for Tecplot output
    integer, allocatable, dimension(:,:) :: tec_index

!   gmsh_index(node) gives the local node numbers for Gmsh output
    integer, allocatable, dimension(:) :: gmsh_index

!   gmsh_index_parsed(node,elem) gives the local node numbers for gmsh parsed
!   output
    integer, allocatable, dimension(:,:) :: gmsh_index_parsed

!   array of size ndimr giving the reference coordinates of the center of
!   the element
    real(dp), allocatable, dimension(:) :: xc

  end type element_t


! type definition of standard nodes->elements and nodes->nodes connection
! data for curves, surfaces, objects ...

  type connectdata_t

    integer :: maxnodnumel = 0    ! maximum number of element connected
                                  ! to nodal points

    integer :: maxnodnumnod = 0    ! maximum number of nodal points connected
                                   ! to nodal points

!   array of length nnodes+1 storing the number of elements connected
!   to nodal points (accumulated).
!   nodnumel(1) = 0,
!   number of elements for nodal point n is
!   nodnumel(n+1) - nodnumel(n)
    integer, allocatable, dimension(:) :: nodnumel

!   elements connected to nodal points
!   number of elements for nodal point n is numel = nodnumel(n+1) - nodnumel(n)
!   the element numbers are stored in
!   nodelem( nodnumel(n) + i ), i = 1, numel
    integer, allocatable, dimension(:) :: nodelem

!   array of length nnodes+1 storing the number of nodal points
!   connected to nodal points (accumulated).
!   nodnumnod(1) = 0,
!   number of nodal points connected to nodal point n is
!   nodnumnod(n+1) - nodnumnod(n)
    integer, allocatable, dimension(:) :: nodnumnod

!   nodal points connected to nodal points
!   number of nodal points connected to nodal point n is
!   numnod = nodnumnod(n+1) - nodnumnod(n) and stored in
!   nodnod( nodnumnod(n) + i ), i = 1, numnod
    integer, allocatable, dimension(:) :: nodnod

  end type connectdata_t


! type definition for the geometric entities: curve, surface, volume

  type geometry_t

    logical :: meshparts = .false. ! has parts been filled in the routine
                                   ! fill_mesh_parts?

    integer :: ndim     = 0        ! space dimension
    integer :: nnodes   = 0        ! number of nodal points
    integer :: nelem    = 0        ! number of elements
    integer :: nblend   = 0        ! number of blended meshes

    integer :: n = 0, m = 0 ! if set to non-zero these numbers indicate
                            ! that the geometry has a regular structure
                            ! as generated in the basic tfem mesh generators
                            ! curves: sequence of n elements
                            ! surface: structured nxm elements
                            ! volume: not available

!   type of the elements, or main type of elements for nblend > 0
    type(element_t) :: element

!   type of elements in blended meshes (nblend > 0)
!   element_blend(m) gives the type of element for blend mesh m
    type(element_t), allocatable, dimension(:) :: element_blend

!   array of length nblend+2 storing the number of nodal points in the
!   geometry with respect to the main mesh and the blended meshes (accumulated).
!   nnodes_blend(1) = 0
!   number of nodal points in the geometry for blend mesh m is
!   nnodes_blend(m+2) - nnodes_blend(m+1)
!   where m=0 is the main mesh.
    integer, allocatable, dimension(:) :: nnodes_blend

!   number of nodal points in an element
!   elnumnod gives the number of nodes in the element
    integer :: elnumnod = 0

!   connections of elements to nodal points
!   topology(node,elem,1) gives the curve/surface/volume node number
!   in element elem
!   topology(node,elem,2) gives the global node number of the local node
!   in element elem
    integer, allocatable, dimension(:,:,:) :: topology

!   array of length nblend+2 storing the number of nodes in the element
!   in the geometry with respect to the main mesh and the blended meshes
!   (accumulated).
!   numnodtop(1) = 0
!   number of nodal points in the element for the geometry of blend mesh m is
!   numnodtop(m+2) - numnodtop(m+1)
!   where m=0 is the main mesh.
    integer, allocatable, dimension(:) :: numnodtop

!   nodal points of the curve, surface or volume
!   nodes(:) gives the nodal point numbers on the curve/surface/volume
!   for each local node number.
    integer, allocatable, dimension(:) :: nodes

!   nodal connection data nodes->elements and nodes->nodes connection
    type(connectdata_t) :: nc

!   number of mesh elements connected to the elements of the geometry
!   Note that nummeshelem(elem)=size(meshelem(elem)%a,2)
!   Not available for volumes.
    integer, allocatable, dimension(:) :: nummeshelem

!   mesh elements connected to the elements of the geometry
!     meshelem(elem)%a(1,i) = elgrpnr
!     meshelem(elem)%a(2,i) = elemnr
!     meshelem(elem)%a(3,i) = efnr
!   with
!        elem: the element number in the geometry
!     elgrpnr: the group number of the connecting mesh element
!      elemnr: the element number of the connecting mesh element
!        efnr: the edge number (for curves) or face number (for surfaces) of
!              the connecting mesh element
!           i: the i^th connected element. There is no limit to the number of
!              connected elements, but it is reasonable to expect that
!              the number of edges for curves in 2D and in the number of faces
!              for surfaces in 3D is either 1 (outer boundary) or 2 (inner
!              boundary).
!   Not available for volumes.
    type(int_array_2d_t), allocatable, dimension(:) :: meshelem

  end type geometry_t


! type definition of object
!
! An object is an entity defined by the coordinates coor in the object
! definition. So the coordinates always need to be defined. It is basically the
! object. The points are connected to the mesh by intersection.
!
! Optionally an object can have a finite element topology on its points,
! using these as the nodal points.
!
! If an object is used to define constraints, the Lagrange multipliers are
! defined in the points/nodes of the object. If there is no topology the
! Lagrange multipliers can only be used for collocation. If there is a
! topology these can be the nodal values of the finite element mesh.
!
! For fluid-solid interaction the object can be extended with multiple
! connections of the points in the object through intersection.

  type object_t

    logical :: meshparts = .false. ! has parts been filled in the routine
                                   ! fill_mesh_parts?

    logical :: topol = .false. ! if .true. points in the object have a
                               ! topology. i.e. a finite element structure

    logical :: intpoints = .false. ! the elements in the object have
                                   ! integration points defined (for use with
                                   ! weak constraints).

    integer :: typeofobject = 0 ! type of object:
                                ! 1: object is connected to the mesh once.
                                !    The connection is defined by refcoor and
                                !    corresponding group and elements in grpelm.
                                ! 2: object is connected to the mesh twice.
                                !    The first connection is defined by refcoor
                                !    and corresponding group and elements in
                                !    grpelm.
                                !    The second connection is defined by
                                !    refcoor2 and corresponding group and
                                !    elements in grpelm2.

    integer :: ndim     = 0        ! space dimension
    integer :: nnodes   = 0        ! number of nodes/points in the object

!   coordinates coor(nnodes,ndim) of the points
    real(dp), allocatable, dimension(:,:) :: coor


!   FIRST INTERSECTION ( typeofobject=1 and 2 ):

!   logical array of size mesh%nelgrp indicating whether the object must be
!   intersected with groups.
!   groups(igrp) = .true.: yes, intersect with group igrp.
!   groups(igrp) = .false.: no intersection with group igrp.
    logical, allocatable, dimension(:) :: groups

!   if >0 the intersection must be in elements that are fully "embedded" in
!   nodes of nodeset, i.e. all nodes in the elements must be in nodeset
    integer :: nodeset = 0

!   reference coordinates refcoor(nnodes,ndim) of the points
    real(dp), allocatable, dimension(:,:) :: refcoor

!   group and element numbers of the points of the object
!   grpelm(node,1) gives the group number
!   grpelm(node,2) gives the element number
!   if point is outside mesh: grpelm(node,:) = 0
    integer, allocatable, dimension(:,:) :: grpelm


!   SECOND INTERSECTION ( typeofobject=2 ):

!   logical array of size mesh%nelgrp indicating whether the object must be
!   intersected with groups.
!   groups2(igrp) = .true.: yes, intersect with group igrp.
!   groups2(igrp) = .false.: no intersection with group igrp.
    logical, allocatable, dimension(:) :: groups2

!   if >0 the intersection must be in elements that are fully "embedded" in
!   nodes of nodeset2, i.e. all nodes in the elements must be in nodeset2
    integer :: nodeset2 = 0

!   reference coordinates refcoor2(nnodes,ndim) of the points
    real(dp), allocatable, dimension(:,:) :: refcoor2

!   group and element numbers of the points of the object
!   grpelm2(node,1) gives the group number
!   grpelm2(node,2) gives the element number
!   if point is outside mesh: grpelm2(node,:) = 0
    integer, allocatable, dimension(:,:) :: grpelm2


!   the following is only for topol=.true.:

    integer :: nelem    = 0        ! number of elements

!   type of the elements
    type(element_t) :: element

!   number of nodal points in an element
!   elnumnod gives the number of nodes in the element
    integer :: elnumnod =0

!   connections of elements to nodal points
!   topology(node,elem) gives the nodal point number in element elem
    integer, allocatable, dimension(:,:) :: topology

!   nodal connection data nodes->elements and nodes->nodes connection
    type(connectdata_t) :: nc


!   the following is only for intpoints=.true.:
!   NOTE: the "real" coordinates of the integration points are defined by
!   isoparametric interpolation of the nodal coordinates.

    integer :: ninti = 0        ! number of integration points in an element

!   reference coordinates xig(ninti,ndimxi) of the integration points in a
!   single element as determined from the integration rule over a single element
!   NOTE: ndimxi=1 for lines, ndimxi=2 for triangles and quadrilaterals and
!         ndimxi=3 for all others (tetrahedrons, hexahedrons, ... ).
    real(dp), allocatable, dimension(:,:) :: xig

!   weights wg(ninti) of the integration points in a single element as
!   determined from the integration rule over a single element
    real(dp), allocatable, dimension(:) :: wg

!   coordinates coor_int(ninti,ndim,elem) of all the integration points.
    real(dp), allocatable, dimension(:,:,:) :: coor_int

!   FIRST INTERSECTION ( typeofobject=1 and 2 ):

!   reference coordinates refcoor_int(ninti,ndim,elem) of all the integration
!   points. NOTE: these are reference coordinates in elements of the mesh!
    real(dp), allocatable, dimension(:,:,:) :: refcoor_int

!   group and element numbers of the integration points of the object
!   grpelm_int(ip,1,elem) gives the group number
!   grpelm_int(ip,2,elem) gives the element number
!   if point is outside mesh: grpelm_int(ip,:,elem) = 0
!   NOTE: elem is the element number in the object and ip is the integration
!   point number within that element.
    integer, allocatable, dimension(:,:,:) :: grpelm_int

!   SECOND INTERSECTION ( typeofobject=2 ):

!   reference coordinates refcoor2_int(ninti,ndim,elem) of all the integration
!   points. NOTE: these are reference coordinates in elements of the mesh!
    real(dp), allocatable, dimension(:,:,:) :: refcoor2_int

!   group and element numbers of the integration points of the object
!   grpelm_int(ip,1,elem) gives the group number
!   grpelm_int(ip,2,elem) gives the element number
!   if point is outside mesh: grpelm_int(ip,:,elem) = 0
!   NOTE: elem is the element number in the object and ip is the integration
!   point number within that element.
    integer, allocatable, dimension(:,:,:) :: grpelm2_int

  end type object_t


! type definition of a block

  type block_t

    integer :: nelem = 0  ! number of elements in this block

!   bounds(dim,1): minimum coordinate of ALL elements in block in direction dim
!   bounds(dim,2): maximum coordinate of ALL elements in block in direction dim
!   NOTE: bounds is such that all the elements are included, but it is not
!   a minimal domain that 'just fits'.
!   NOTE: the bounds are computed using minval and maxval intrinsics. This
!   means that if nelem=0: bounds(dim,1) is the largest possible value and
!   bounds(dim,2) is the smallest possible value.
    real(dp), allocatable, dimension(:,:) :: bounds

!   elbounds(elem,dim,1), elem=1,nelem: minimum coordinate of element elem in
!                                       block in direction dim
!   elbounds(elem,dim,2), elem=1,nelem: maximum coordinate of element elem in
!                                       block in direction dim
!   NOTE: elbounds is such that an element is fully included, but it is not
!   a minimal domain that 'just fits' around an element.
    real(dp), allocatable, dimension(:,:,:) :: elbounds

!   elements(elem,1)), elem=1,nelem: group number
!   elements(elem,2)), elem=1,nelem: element number
    integer, allocatable, dimension(:,:) :: elements

  end type block_t


! type definition of a structured block

  type sblock_t

    integer :: nelem = 0 ! number of elements in this sblock

!   elbounds(elem,dim,1), elem=1,nelem: minimum coordinate of element elem in
!                                       sblock in direction dim
!   elbounds(elem,dim,2), elem=1,nelem: maximum coordinate of element elem in
!                                       sblock in direction dim
!   NOTE: elbounds is such that an element is fully included, but it is not
!   a minimal domain that 'just fits' around an element.
    real(dp), allocatable, dimension(:,:,:) :: elbounds

!   elements(elem,1)), elem=1,nelem: group number
!   elements(elem,2)), elem=1,nelem: element number
    integer, allocatable, dimension(:,:) :: elements

  end type sblock_t


! type definition of sblocks (structured blocks)

  type sblocks_t

!   logical to indicate whether the structured blocks have been filled
    logical :: filled = .false.

!   number of structured blocks. This is equal to size(sblks) if filled=.true.
    integer :: nsblocks  = 0

!   nbl(dim) gives the number of structured blocks in coordinate direction dim
    integer, allocatable, dimension(:) :: nbl

!   x0(dim): start coordinate of the structured grid of blocks in direction dim
!   dx(dim): width of a sblock in direction dim
!   Note: the grid lines separating the sblocks are given by (for 3D):
!    x = x0 + [i,j,k]*dx with i=0,..,nbl(1), j=0,..,nbl(2), k=0,..,nbl(3)
    real(dp), allocatable, dimension(:) :: x0, dx

!   sblks(blk): gives the information of sblock blk
!   Notes:
!     * number of sblocks = nbl(1) * nbl(2) * nbl(3) = size(sblks)
!     * blk = i + ( j - 1 ) * nbl(1) + ( k - 1 ) * nbl(1) * nbl(2)
!   Each sblock contains a number of elements that are "candidates" for
!   intersection for a point having its coordinates within the sblock.
!   Note, an element can be part of more than one sblock.
    type(sblock_t), allocatable, dimension(:) :: sblks

  end type sblocks_t


! type definition of a nodal block

  type nodblock_t

    integer :: nnodes  ! number of nodes in this nodblock

!   nodes(node)), node=1,nnodes
    integer, allocatable, dimension(:) :: nodes

  end type nodblock_t


! type definition of nodblocks (nodal blocks)

  type nodblocks_t

!   logical to indicate whether the nodal blocks have been filled
    logical :: filled = .false.

!   number of nodal blocks. This is equal to size(nodblks) if filled=.true.
    integer :: nnodblocks  = 0

!   nbl(dim) gives the number of nodal blocks in coordinate direction dim
    integer, allocatable, dimension(:) :: nbl

!   x0(dim): start coordinate of the nodal blocks in direction dim
!   dx(dim): width of a nodal block in direction dim
!   Note: the grid lines separating the nodblocks are given by (for 3D):
!    x = x0 + [i,j,k]*dx with i=0,..,nbl(1), j=0,..,nbl(2), k=0,..,nbl(3)
    real(dp), allocatable, dimension(:) :: x0, dx

!   nodblks(blk): gives the information of nodblock blk
!   Notes:
!     * number of nodblocks = nbl(1) * nbl(2) * nbl(3) = size(nodblks)
!     * blk = i + ( j - 1 ) * nbl(1) + ( k - 1 ) * nbl(1) * nbl(2)
!   Each nodblock contains a number of nodes that are "candidates" for
!   overlapping for a point having its coordinates within the nodblock.
!   Note, a node can be part of more than one nodblock.
    type(nodblock_t), allocatable, dimension(:) :: nodblks

  end type nodblocks_t


! type definition of an element set

  type elementset_t

    logical :: nodes_created = .false.  ! is the array nodes created?

!   number of elements in the elementset
    integer :: nelem = 0

!   number of nodal points the elements are connected to (nnodes) and the
!   nodal point numbers (nodes) the elements are connected to.
!   NOTE: these are optional and requires an additional call to add_to_mesh!
    integer :: nnodes = 0
    integer, allocatable, dimension(:) :: nodes

!   number of elements given for element group
!   grpnumel(elgrp) gives the number of elements given for group elgrp
    integer, allocatable, dimension(:) :: grpnumel

!   elements(elgrp)%a(:) the element numbers of the element set for group elgrp
    type(int_array_1d_t), allocatable, dimension(:) :: elements

  end type elementset_t


! multilevel elements

   type melements_t

    integer :: level = 1           ! The refinement level of the element

    integer, allocatable, dimension(:) :: selem    ! Sibling elements

    integer, allocatable, dimension(:,:) :: selem3D ! Sibling elements in 3D

    integer, allocatable, dimension(:) :: pelem    ! Parent elements

    integer :: reftype = 0         ! Type of refinement :
                                   ! reftype = 0 : No refinement
                                   ! reftype = 1 : Regular(Red) refinement
                                   ! reftype = 2 : Unregular(Green) refinment
                                   ! reftype = 3 : Refined to keep level
                                   !               differance of side elements
                                   !               =< 1

    logical :: leaf = .true.       ! True if element has no siblings

    real(dp) :: lside = 0          ! Longest side length

    integer :: mark = 0            ! Mark for refinement
                                   ! mark = 0 : Do not refine
                                   ! mark = 1 : Regular refine
                                   ! mark = 2 : Green refine
                                   ! mark = 3 : Unrefine

  end type melements_t


! array of multilevel elements
  type me_array_1d_t
    type(melements_t), allocatable, dimension(:) :: a
  end type me_array_1d_t


! type definition of mesh

  type mesh_t

    logical :: renumber = .false.   ! nodal point renumbering
    logical :: meshgen = .false.    ! has basic mesh (topolopy, coordinates, &
                                    ! points, curves, surfaces) been generated?
    logical :: meshparts = .false.  ! has parts been filled in the routine
                                    ! fill_mesh_parts?
    logical :: multlvlref = .false. ! has the mesh been refined?
    integer :: ndim     = 0         ! space dimension
    integer :: nnodes   = 0         ! number of nodal points
    integer :: nelem    = 0         ! number of elements
    integer :: nelgrp   = 0         ! number of element groups
    integer :: nblend   = 0         ! number of blend meshes
    integer :: npoints  = 0         ! number of points
    integer :: ncurves  = 0         ! number of curves
    integer :: nsurfaces  = 0       ! number of surfaces
    integer :: nvolumes  = 0        ! number of volumes
    integer :: maxnodnumnod = 0     ! maximum number of nodal points connected
                                    ! to nodal points
    integer :: ngluepoints = 0      ! number of glue points
    integer :: nobjects  = 0        ! number of objects
    integer :: nblocks  = 0         ! number of blocks
    integer :: nnodesets  = 0       ! number of set of nodes
    integer :: nelementsets  = 0    ! number of set of elements
    integer :: maxelnumnod = 0      ! maximum number of nodal points in an
                                    ! element

!   array of length mesh%nnodes containing renumbered nodal points
!   only present if renumber==.true.
!   new node i has old node number nodperm(i)
!   the new numbering is only used for renumbering the degrees of freedom
!   for matrix partitioning (problem%degfdperm).
    integer, allocatable, dimension(:) :: nodperm
!
!   type of elements. For nblend > 0 this is the main type of elements.
!   element(elgrp) gives the type of an element in group elgrp
    type(element_t), allocatable, dimension(:) :: element

!   type of elements in the blended meshes for nblend > 0
!   element_blend(elgrp,m) gives the type of an element in group elgrp
!   for blend mesh m
    type(element_t), allocatable, dimension(:,:) :: element_blend

!   array of length nblend+2 storing the number of nodal points in the
!   mesh with respect to the blended meshes (accumulated).
!   nnodes_blend(1) = 0
!   number of nodal points in the mesh for blend mesh m is
!   nnodes_blend(m+2) - nnodes_blend(m+1)
!   where m=0 is the main mesh.
    integer, allocatable, dimension(:) :: nnodes_blend

!   number of nodal points in an element
!   elnumnod(elgrp) gives the number of nodes in the element for group elgrp
    integer, allocatable, dimension(:) :: elnumnod

!   number of elements in an element group
!   grpnumel(elgrp) gives the number of elements in group elgrp
    integer, allocatable, dimension(:) :: grpnumel

!   connections of elements to nodal points
!   topology(elgrp)%a(node,elem) gives the global node number of the local node
!   in element elem within group elgrp.
    type(int_array_2d_t), allocatable, dimension(:) :: topology

!   array of shape (nelgrp,nblend+2) storing the number of nodes in the
!   elements with respect to the blended meshes (accumulated).
!   numnodtop(:,1) = 0
!   number of nodal points in the element for blend mesh m is
!   numnodtop(:,m+2) - numnodtop(:,m+1)
!   where m=0 is the main mesh.
    integer, allocatable, dimension(:,:) :: numnodtop

!   array of length mesh%nnodes+1 storing the number of elements connected
!   to nodal points (accumulated).
!   nodnumel(1) = 0,
!   number of elements for nodal point n is
!   nodnumel(n+1) - nodnumel(n)
    integer, allocatable, dimension(:) :: nodnumel

!   elements connected to nodal points
!   number of elements for nodal point n is numel = nodnumel(n+1) - nodnumel(n)
!   the group numbers of the elements are stored in
!   nodelem( nodnumel(n) + i, 1 ), i = 1, numel
!   the element numbers (relative to the group) are stored in
!   nodelem( nodnumel(n) + i, 2 ), i = 1, numel
    integer, allocatable, dimension(:,:) :: nodelem

!   array of length mesh%nnodes+1 storing the number of nodal points
!   connected to nodal points (accumulated).
!   nodnumnod(1) = 0,
!   number of nodal points connected to nodal point n is
!   nodnumnod(n+1) - nodnumnod(n)
    integer, allocatable, dimension(:) :: nodnumnod

!   nodal points connected to nodal points
!   number of nodal points connected to nodal point n is
!   numnod = nodnumnod(n+1) - nodnumnod(n) and stored in
!   nodnod( nodnumnod(n) + i ), i = 1, numnod
    integer, allocatable, dimension(:) :: nodnod

!   points(pnt) gives the nodal point number of point pnt
    integer, allocatable, dimension(:) :: points

!   curves(crv) gives the information of curve crv
    type(geometry_t), allocatable, dimension(:) :: curves

!   surface(srf) gives the information of surface srf
    type(geometry_t), allocatable, dimension(:) :: surfaces

!   volume(vlm) gives the information of volume vlm
    type(geometry_t), allocatable, dimension(:) :: volumes

!   gluepoints gives the nodal points that are `glued' together.
!   this is used to glue together mesh parts for the filling of sidelem.
!   gluepoint(p,1) is glued to gluepoint(p,2) and vice versa.
!   properties: 1) gluepoint(p,1) < gluepoint(p,2)
!               2) gluepoint(:,1) is sorted
    integer, allocatable, dimension(:,:) :: gluepoints

!   elements connected to the sides of elements
!   sidelem(elgrp)%a(side,elem,1) = elgrpnr
!   sidelem(elgrp)%a(side,elem,2) = elemnr
!   sidelem(elgrp)%a(side,elem,3) = sidenr
!   sidelem(elgrp)%a(side,elem,4) = orient
!   with
!       elgrp: the group number of the current element
!        side: the side number of the current element
!        elem: the element number of the current element
!     elgrpnr: the group number of the adjacent element
!      elemnr: the element number of the adjacent element
!      sidenr: the side number of the adjacent element
!      orient: the value of |orient| is:
!                 the orientation of the adjacent element.
!                 This is the number of rotations (cyclic) that has to be made
!                 over the vertices of the adjacent element + 1, such that the
!                 first vertex of the side is the same on both sides.
!                 So: 1<=|orient|<=number of vertices. For 1D it is always 1,
!                 for 2D it is always 2 (if numbering is consistent).
!              the sign of orient means:
!                 orient >0: normal (internal) mesh connection
!                 orient <0: connection generated by gluepoints
!   if there is no adjacent element: elemnr=elgrpnr=sidenr=orient=0
    type(int_array_3d_t), allocatable, dimension(:) :: sidelem

!   objects(obj) gives the information of object obj
    type(object_t), allocatable, dimension(:) :: objects

!   coordinates coor(nnodes,ndim)
    real(dp), allocatable, dimension(:,:) :: coor

!   blocks(blk): gives the information of block blk
!   This is an ordering of the elements of a mesh into blocks of elements.
    type(block_t), allocatable, dimension(:) :: blocks

!   nodesets(set)%a(:) the global node numbers of the nodeset set
    type(int_array_1d_t), allocatable, dimension(:) :: nodesets

!   elementsets(set) the elementset set
    type(elementset_t), allocatable, dimension(:) :: elementsets

!   structured blocks. Note, that structured blocks are optional.
    type(sblocks_t) :: sblocks

!   nodal blocks. Note, that nodal blocks are optional.
    type(nodblocks_t) :: nodblocks

!   multilevel elements for refinement
    type(melements_t), allocatable, dimension(:) :: melements

  end type mesh_t


! interface for generic delete subroutine

  interface delete
    module procedure delete_mesh
  end interface delete

! interface for generic check subroutine

  interface check
    module procedure check_mesh
  end interface check

! interface for generic copy subroutine

  interface copy
    module procedure copy_element_basic, copy_element_array_basic, &
                     copy_element_array2_basic, &
                     copy_geometry_basic, copy_geometry_array_basic, &
                     copy_object_basic, copy_object_array_basic, &
                     copy_block, copy_block_array, &
                     copy_sblock, copy_sblocks, &
                     copy_nodblock, copy_nodblocks, &
                     copy_mesh_basic
  end interface copy

contains


! check mesh

  subroutine check_mesh ( mesh, name_of_routine )

    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: name_of_routine

    if ( .not. mesh%meshparts ) then
      if ( .not. mesh%meshgen ) then
        write(*,'(/3a/)') &
          'Error ', name_of_routine, ': mesh has not been generated '
        stop
      else
        write(*,'(/3a/3(a/))') &
          'Error ', name_of_routine, ': mesh is incomplete ', &
          ' After meshgeneration or reading a mesh from a file the routine ', &
          ' fill_mesh_parts must be called to fill the remaining parts ', &
          ' of the mesh '
        stop
      end if
    end if

  end subroutine check_mesh


! check for consistent nblend values in mesh

  subroutine check_nblend ( mesh, name_of_routine )

    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: name_of_routine

    if ( any ( mesh%curves(:mesh%ncurves)%nblend /= mesh%nblend ) .or. &
         any ( mesh%surfaces(:mesh%nsurfaces)%nblend /= mesh%nblend ) .or. &
         any ( mesh%volumes(:mesh%nvolumes)%nblend /= mesh%nblend ) ) then
      write(*,'(/3a)') &
        'Error ', name_of_routine, ': inconsistent values for nblend in mesh'
      stop
    end if

  end subroutine check_nblend


! check for blended meshes in mesh and give an error (or optional warning)

  subroutine check_blend ( mesh, name_of_routine, warning, keyword, comment )

    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: name_of_routine
    logical, intent(in), optional :: warning
    character(len=*), intent(in), optional :: keyword, comment

!   check if blended meshes are present

    if ( mesh%nblend > 0 ) then

      if ( set_optional ( variable=warning, default=.false. ) ) then
        write(*,'(4a)') &
          'Warning ', name_of_routine, ': mesh contains blended meshes.', &
          ' Not yet supported. Use at your own risk.'
        if ( present(keyword) ) then
          write(*,'(2a)') &
            'Keyword: ', keyword
        end if
        if ( present(comment) ) then
          write(*,'(2a)') &
            'Comment: ', comment
        end if
      else
        write(*,'(4a)') &
          'Error ', name_of_routine, ': mesh contains blended meshes.', &
          ' Not yet supported.'
        if ( present(keyword) ) then
          write(*,'(2a)') &
            'Keyword: ', keyword
        end if
        if ( present(comment) ) then
          write(*,'(2a)') &
            'Comment: ', comment
        end if
        stop
      end if

    end if

  end subroutine check_blend


! check blocks

  subroutine check_blocks ( mesh )

    type(mesh_t), intent(in) :: mesh

    integer :: nelemb

    if ( mesh%nblocks > 0 ) then
      nelemb = sum(mesh%blocks(:mesh%nblocks)%nelem)
      if ( nelemb /= mesh%nelem ) then
        write(*,'(/a/2(a,i0/))') &
          'Warning check_blocks: ', &
          '  Number of elements in blocks = ', nelemb, &
          '  Number of elements in mesh = ', mesh%nelem
      end if
    end if

  end subroutine check_blocks


! Get coordinates of element given by (elgrp,elem)

  subroutine get_coordinates ( mesh, elgrp, elem, x, onlyshapenod, blend, &
    nodesx )

    type(mesh_t), intent(in) :: mesh

!   element group elgrp and element nr elem of the element
    integer, intent(in) :: elgrp, elem

!   the nodal coordinates in this element
!   x(node,dim) is the coordinate of nodal point node in dim direction
    real(dp), intent(out), dimension(:,:) :: x

!   optional parameter to indicate that only numshapenod nodes are
!   extracted. This affects only the elements with additional internal nodes.
!   default is onlyshapenod = .false. (=all nodes)
    logical, intent(in), optional :: onlyshapenod

!   the mesh chosen for the coordinates
!   blend=-1 get all coordinates of the element
!   blend=0  get coordinates of the element of the main mesh (mesh%element)
!   blend>0  get coordinates of the element of the mesh given by
!            mesh%element_blend(:,blend),
!   default=0
    integer, intent(in), optional :: blend

!   if present it contains the node numbers of the coordinates x
    integer, dimension(:), intent(out), optional :: nodesx

    logical :: onlysh
    integer :: s, lblend, nodes(size(x,1))


    onlysh = set_optional ( variable = onlyshapenod, default=.false. )
    lblend = set_optional ( variable = blend, default=0 )

    if ( lblend == 0 ) then

!     get coordinates in (main) mesh

      if ( onlysh ) then
        nodes = mesh%topology(elgrp)%a(:mesh%element(elgrp)%numshapenod,elem)
      else
        nodes = mesh%topology(elgrp)%a(:mesh%element(elgrp)%numnod,elem)
      end if

    else if ( lblend > 0 ) then

!     get coordinates of element in blend mesh lblend

      s = mesh%numnodtop(elgrp,lblend+1)

      if ( onlysh ) then
        nodes = &
           mesh%topology(elgrp)%a(s+1:s+mesh%element(elgrp)%numshapenod,elem)
      else
        nodes = &
          mesh%topology(elgrp)%a(s+1:mesh%numnodtop(elgrp,lblend+2),elem)
      end if

    else

!     get all coordinates of the element

      nodes = mesh%topology(elgrp)%a(:,elem)

    end if

    x = mesh%coor(nodes,:)

    if ( present(nodesx) ) nodesx = nodes

  end subroutine get_coordinates


! Get coordinates of nodes on the sides of an element given by (elgrp,elem)

  subroutine get_coordinates_sides ( mesh, elgrp, elem, x, blend, nodesx )

    type(mesh_t), intent(in) :: mesh

!   element group elgrp and element nr elem of the element
    integer, intent(in) :: elgrp, elem

!   the coordinates in this element of the side nodes
!   x(node,dim,side) is the coordinate of nodal point node of sidenr
!   side in direction dim
    real(dp), intent(out), dimension(:,:,:) :: x

!   the mesh chosen for the coordinates
!   blend=0  get coordinates of the element of the main mesh (mesh%element)
!   blend>0  get coordinates of the element of the mesh given by
!            mesh%element_blend(:,blend),
!   default=0
    integer, intent(in), optional :: blend

!   if present it contains the node numbers of the coordinates x
!   nodesx(node,side) is the node number of nodal point node of sidenr side
    integer, dimension(:,:), intent(out), optional :: nodesx

    integer :: side, lblend, nodes(size(x,1),size(x,3))

    lblend = set_optional ( variable = blend, default=0 )

    if ( lblend == 0 ) then

      do side = 1, mesh%element(elgrp)%numsides
        nodes(:,side) = &
              mesh%topology(elgrp)%a(mesh%element(elgrp)%sidnod(:,side),elem)
        x(:,:,side) = mesh%coor(nodes(:,side),:)
      end do

    else

      do side = 1, mesh%element_blend(elgrp,lblend)%numsides
        nodes(:,side) = &
              mesh%topology(elgrp)%a(mesh%numnodtop(elgrp,lblend+1) + &
                        mesh%element_blend(elgrp,lblend)%sidnod(:,side),elem)
        x(:,:,side) = mesh%coor(nodes(:,side),:)
      end do

    end if

    if ( present(nodesx) ) nodesx = nodes

  end subroutine get_coordinates_sides


! Get coordinates of nodes on a single side of an element given by (elgrp,elem)

  subroutine get_coordinates_side ( mesh, elgrp, elem, side, x, blend, nodesx )

    type(mesh_t), intent(in) :: mesh

!   element group elgrp, element nr elem, and side number of the element
    integer, intent(in) :: elgrp, elem, side

!   the coordinates in this element of the sidenr side
!   x(node,dim) is the coordinate of nodal point node of sidenr
!   side in direction dim
    real(dp), intent(out), dimension(:,:) :: x

!   the mesh chosen for the coordinates
!   blend=0  get coordinates of the element of the main mesh (mesh%element)
!   blend>0  get coordinates of the element of the mesh given by
!            mesh%element_blend(:,blend),
!   default=0
    integer, intent(in), optional :: blend

!   if present it contains the node numbers of the coordinates x
    integer, dimension(:), intent(out), optional :: nodesx

    integer :: lblend, nodes(size(x,1))

    lblend = set_optional ( variable = blend, default=0 )

    if ( lblend == 0 ) then

      nodes = mesh%topology(elgrp)%a(mesh%element(elgrp)%sidnod(:,side),elem)

    else

      nodes = mesh%topology(elgrp)%a(mesh%numnodtop(elgrp,lblend+1) &
              + mesh%element_blend(elgrp,lblend)%sidnod(:,side),elem)

    end if

    x = mesh%coor(nodes,:)

    if ( present(nodesx) ) nodesx = nodes

  end subroutine get_coordinates_side


! Get coordinates of element elem on a geometry given by curve, surface,
! volume or by (ndimr,geometry)

  subroutine get_coordinates_geometry ( mesh, elem, x, curve, surface, volume, &
    ndimr, geometry, blend, nodesx )

    type(mesh_t), intent(in) :: mesh

!   element nr elem on the curve, surface or volume
    integer, intent(in) :: elem

!   the coordinates in this element of the curve/surface
!   x(node,dim) is the coordinate of nodal point node of this element in
!   direction dim.
    real(dp), intent(out), dimension(:,:) :: x

!   curve/surface/volume number (legacy)
    integer, intent(in), optional :: curve, surface, volume

!   dimension of reference space (ndimr) and geometry number
    integer, intent(in), optional :: ndimr, geometry

!   the mesh chosen for the coordinates
!   blend=-1 get all coordinates of the element
!   blend=0  get coordinates of the element of the main mesh (geometry%element)
!   blend>0  get coordinates of the element of the mesh given by
!            geometry%element_blend(:,blend),
!   default=0
    integer, intent(in), optional :: blend

!   if present it contains the node numbers of the coordinates x
    integer, dimension(:), intent(out), optional :: nodesx

    integer :: lcurve, lsurface, lvolume, lblend, nodes(size(x,1))

!   traditional interface (legacy)

    lcurve = set_optional ( variable=curve, default=0 )
    lsurface = set_optional ( variable=surface, default=0 )
    lvolume = set_optional ( variable=volume, default=0 )
    lblend = set_optional ( variable = blend, default=0 )

!   interface with dimension of reference space of geometry

    if ( present(ndimr) .and. present(geometry) ) then
      select case (ndimr)
        case(1); lcurve = geometry
        case(2); lsurface = geometry
        case(3); lvolume = geometry
        case default
          call errormsg_case_default ( 'get_coordinates_geometry', 'ndimr', &
            int_value=ndimr )
      end select
    end if

    if ( lblend == 0 ) then

      if ( lcurve > 0 ) then
        nodes = mesh%curves(lcurve)%topology(&
                            :mesh%curves(lcurve)%element%numnod,elem,2)
      else if ( lsurface > 0 ) then
        nodes = mesh%surfaces(lsurface)%topology(&
                        :mesh%surfaces(lsurface)%element%numnod,elem,2)
      else if ( lvolume > 0 ) then
        nodes = mesh%volumes(lvolume)%topology(&
                          :mesh%volumes(lvolume)%element%numnod,elem,2)
      end if

    else if ( lblend > 0 ) then

      if ( lcurve > 0 ) then
        nodes = mesh%curves(lcurve)%topology( &
                   mesh%curves(lcurve)%numnodtop(lblend+1)+1: &
                          mesh%curves(lcurve)%numnodtop(lblend+2),elem,2)
      else if ( lsurface > 0 ) then
        nodes = mesh%surfaces(lsurface)%topology( &
                   mesh%surfaces(lsurface)%numnodtop(lblend+1)+1: &
                      mesh%surfaces(lsurface)%numnodtop(lblend+2),elem,2)
      else if ( lvolume > 0 ) then
        nodes = mesh%volumes(lvolume)%topology( &
                   mesh%volumes(lvolume)%numnodtop(lblend+1)+1: &
                        mesh%volumes(lvolume)%numnodtop(lblend+2),elem,2)
      end if

    else

      if ( lcurve > 0 ) then
        nodes = mesh%curves(lcurve)%topology(:,elem,2)
      else if ( lsurface > 0 ) then
        nodes = mesh%surfaces(lsurface)%topology(:,elem,2)
      else if ( lvolume > 0 ) then
        nodes = mesh%volumes(lvolume)%topology(:,elem,2)
      end if

    end if

    x = mesh%coor(nodes,:)

    if ( present(nodesx) ) nodesx = nodes

  end subroutine get_coordinates_geometry


! Get coordinates of element on an object given by (object,elem)

  subroutine get_coordinates_object ( mesh, elem, x, object )

    type(mesh_t), intent(in) :: mesh

!   element nr elem on the curve or surface
    integer, intent(in) :: elem

!   the coordinates in this element of the object
!   x(node,dim) is the coordinate of nodal point node of this element in
!   direction dim.
    real(dp), intent(out), dimension(:,:) :: x

!   object number
    integer, intent(in) :: object

    if ( .not. mesh%objects(object)%topol ) then
      write(*,'(/a/)') &
        'Error get_coordinates_object: object has no topology '
      stop
    end if

    x = mesh%objects(object)%coor(mesh%objects(object)%topology(:,elem),:)

  end subroutine get_coordinates_object


! check whether node is in elementgroups

  function node_in_groups ( mesh, groups, nodenr )

    type(mesh_t), intent(in) :: mesh
    integer, dimension(:), intent(in) :: groups
    integer, intent(in) :: nodenr
    logical :: node_in_groups

    integer :: start, numel, elm

    start = mesh%nodnumel(nodenr)
    numel = mesh%nodnumel(nodenr+1) - mesh%nodnumel(nodenr)

    node_in_groups = .false.
    do elm = 1, numel
      if ( any ( groups == mesh%nodelem(start+elm,1) ) ) then
        node_in_groups = .true.
        exit
      end if
    end do

  end function node_in_groups


! copy element1 to element: only the data needed for the basic mesh

  subroutine copy_element_basic ( element1, element )

    type(element_t), intent(in) :: element1
    type(element_t), intent(out) :: element

    element%globalshape = element1%globalshape
    element%elshape = element1%elshape
    element%p = element1%p
    element%numnod = element1%numnod
    element%ndim = element1%ndim

  end subroutine copy_element_basic


! copy element1 to element (array): only the data needed for the basic mesh

  subroutine copy_element_array_basic ( element1, element )

    type(element_t), dimension(:), intent(in) :: element1
    type(element_t), dimension(:), intent(out) :: element

    integer :: i

    do i = 1, size(element1)
      call copy_element_basic ( element1(i), element(i) )
    end do

  end subroutine copy_element_array_basic


! copy element1 to element (2D array): only the data needed for the basic mesh

  subroutine copy_element_array2_basic ( element1, element )

    type(element_t), dimension(:,:), intent(in) :: element1
    type(element_t), dimension(:,:), intent(out) :: element

    integer :: i, j

    do i = 1, size(element1,1)
      do j = 1, size(element1,2)
        call copy_element_basic ( element1(i,j), element(i,j) )
      end do
    end do

  end subroutine copy_element_array2_basic


! copy geometry1 to geometry: only the data needed for the basic mesh

  subroutine copy_geometry_basic ( geometry1, geometry, newnodenrs )

    type(geometry_t), intent(in) :: geometry1
    type(geometry_t), intent(out) :: geometry
!   if present nodes get a new number: newnr=newnodenrs(oldnr)
    integer, dimension(:), intent(in), optional :: newnodenrs

    integer :: elem

    geometry%meshparts = .false.
    geometry%ndim = geometry1%ndim
    geometry%nnodes = geometry1%nnodes
    geometry%nelem = geometry1%nelem
    geometry%nblend = geometry1%nblend
    geometry%n = geometry1%n
    geometry%m = geometry1%m

    allocate ( geometry%element_blend(geometry%nblend) )

    call copy ( geometry1%element, geometry%element )
    call copy ( geometry1%element_blend, geometry%element_blend )

    geometry%nnodes_blend = geometry1%nnodes_blend

    if ( present(newnodenrs) ) then

      geometry%nodes = newnodenrs(geometry1%nodes)
      allocate( geometry%topology(size(geometry1%topology,1),geometry%nelem,2) )
      geometry%topology(:,:,1) = geometry1%topology(:,:,1)
      do elem = 1, geometry%nelem
        geometry%topology(:,elem,2) = newnodenrs(geometry1%topology(:,elem,2))
      end do

    else

      geometry%nodes = geometry1%nodes
      allocate( geometry%topology(size(geometry1%topology,1),geometry%nelem,2) )
      geometry%topology = geometry1%topology

    end if

  end subroutine copy_geometry_basic


! copy geometry1 to geometry: only the data needed for the basic mesh

  subroutine copy_geometry_array_basic ( geometry1, geometry, newnodenrs )

    type(geometry_t), dimension(:), intent(in) :: geometry1
    type(geometry_t), dimension(:), intent(out) :: geometry
!   if present nodes get a new number: newnr=newnodenrs(oldnr)
    integer, dimension(:), intent(in), optional :: newnodenrs

    integer :: i

    do i = 1, size(geometry1)
      call copy_geometry_basic ( geometry1(i), geometry(i), newnodenrs )
    end do

  end subroutine copy_geometry_array_basic


! copy object1 to object: only the data needed for the basic mesh

  subroutine copy_object_basic ( object1, object )

    type(object_t), intent(in) :: object1
    type(object_t), intent(out) :: object

    object%meshparts = .false.
    object%topol = object1%topol
    object%intpoints = object1%intpoints
    object%typeofobject = object1%typeofobject
    object%ndim = object1%ndim
    object%nnodes = object1%nnodes

    object%coor = object1%coor

    object%groups = object1%groups
    object%refcoor = object1%refcoor
    object%grpelm = object1%grpelm
    object%nodeset = object1%nodeset

    if ( object%typeofobject == 2 ) then

      object%groups2 = object1%groups2
      object%refcoor2 = object1%refcoor2
      object%grpelm2 = object1%grpelm2
      object%nodeset2 = object1%nodeset2

    end if

    if ( object%topol ) then

      object%nelem = object1%nelem
      object%elnumnod = object1%elnumnod

      call copy_element_basic ( object1%element, object%element )

      object%topology = object1%topology

      if ( object%intpoints ) then

        object%ninti = object1%ninti

        object%xig = object1%xig
        object%wg = object1%wg
        object%coor_int = object1%coor_int

        object%refcoor_int = object1%refcoor_int
        object%grpelm_int = object1%grpelm_int

        if ( object%typeofobject == 2 ) then

          object%refcoor2_int = object1%refcoor2_int
          object%grpelm2_int = object1%grpelm2_int

        end if

      end if

    end if

  end subroutine copy_object_basic


! copy geometry1 to geometry: only the data needed for the basic mesh

  subroutine copy_object_array_basic ( object1, object )

    type(object_t), dimension(:), intent(in) :: object1
    type(object_t), dimension(:), intent(out) :: object

    integer :: i

    do i = 1, size(object1)
      call copy_object_basic ( object1(i), object(i) )
    end do

  end subroutine copy_object_array_basic


! copy block1 to block

  subroutine copy_block ( block1, block )

    type(block_t), intent(in) :: block1
    type(block_t), intent(out) :: block

    block%nelem = block1%nelem

    allocate ( block%bounds(size(block1%bounds,1),2) )
    allocate ( block%elbounds(block%nelem,size(block1%elbounds,2),2) )
    allocate ( block%elements(block%nelem,2) )
    block%bounds = block1%bounds
    block%elbounds = block1%elbounds
    block%elements = block1%elements

  end subroutine copy_block


! copy array of block1 to block

  subroutine copy_block_array ( block1, block )

    type(block_t), dimension(:), intent(in) :: block1
    type(block_t), dimension(:), intent(out) :: block

    integer :: i

    do i = 1, size(block1)
      call copy_block ( block1(i), block(i) )
    end do

  end subroutine copy_block_array


! copy sblock1 to sblock

  subroutine copy_sblock ( sblock1, sblock )

    type(sblock_t), intent(in) :: sblock1
    type(sblock_t), intent(out) :: sblock

    sblock = sblock1

  end subroutine copy_sblock


! copy sblocks1 to sblocks

  subroutine copy_sblocks ( sblocks1, sblocks )

    type(sblocks_t), intent(in) :: sblocks1
    type(sblocks_t), intent(out) :: sblocks

    if ( sblocks%filled ) then
      write(*,'(/3(a/))') &
        'Warning in copy_sblocks: ', &
        ' Structured blocks have already been defined in sblocks.', &
        ' Existing structured blocks in sblocks will be removed.'
      call delete_sblocks ( sblocks )
    end if

    if ( .not. sblocks1%filled ) return ! sblocks remain empty/undefined

    sblocks = sblocks1

  end subroutine copy_sblocks


! copy nodblock1 to nodblock

  subroutine copy_nodblock ( nodblock1, nodblock )

    type(nodblock_t), intent(in) :: nodblock1
    type(nodblock_t), intent(out) :: nodblock

    nodblock = nodblock1

  end subroutine copy_nodblock


! copy nodblocks1 to nodblocks

  subroutine copy_nodblocks ( nodblocks1, nodblocks )

    type(nodblocks_t), intent(in) :: nodblocks1
    type(nodblocks_t), intent(out) :: nodblocks

    if ( nodblocks%filled ) then
      write(*,'(/3(a/))') &
        'Warning in copy_nodblocks: ', &
        ' Structured blocks have already been defined in nodblocks.', &
        ' Existing structured blocks in nodblocks will be removed.'
      call delete_nodblocks ( nodblocks )
    end if

    if ( .not. nodblocks1%filled ) return ! nodblocks remain empty/undefined

    nodblocks = nodblocks1

  end subroutine copy_nodblocks


! copy mesh1 to mesh: only the data needed for the basic mesh

  subroutine copy_mesh_basic ( mesh1, mesh )

    use limits_m, only: MAXPOINTS, MAXCURVES, MAXOBJECTS, MAXSURFACES, &
                        MAXBLOCKS, MAXVOLUMES, MAXNODESETS, MAXELEMENTSETS

    type(mesh_t), intent(in) :: mesh1
    type(mesh_t), intent(inout) :: mesh

    integer :: nodeset, elementset

!   some testing

    if ( mesh%meshgen ) then
      write(*,'(/a/a/)') &
        'Error in copy_mesh_basic: mesh is not empty', &
        ' either use a new mesh_t variable or delete old mesh '
      stop
    end if

    if ( .not. mesh1%meshgen ) then
      write(*,'(/a/)') &
        'Error in copy_mesh_basic: mesh1 is empty'
      stop
    end if

!   copy all elements

    mesh%ndim = mesh1%ndim
    mesh%nelem = mesh1%nelem
    mesh%nnodes = mesh1%nnodes
    mesh%nelgrp = mesh1%nelgrp
    mesh%nblend = mesh1%nblend

!   internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, mesh%nblend ) )

    call copy ( mesh1%element, mesh%element )
    call copy ( mesh1%element_blend, mesh%element_blend )

    mesh%nnodes_blend = mesh1%nnodes_blend

    mesh%grpnumel = mesh1%grpnumel

    mesh%topology = mesh1%topology
    mesh%coor = mesh1%coor

!   points

    mesh%npoints = mesh1%npoints

    allocate( mesh%points(max(mesh%npoints,MAXPOINTS)) )

    mesh%points(1:mesh%npoints) = mesh1%points(1:mesh1%npoints)

!   curves

    mesh%ncurves = mesh1%ncurves

    allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )

    call copy ( mesh1%curves(:mesh1%ncurves), mesh%curves(:mesh%ncurves) )

!   surfaces

    mesh%nsurfaces = mesh1%nsurfaces

    allocate( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )

    call copy ( mesh1%surfaces(:mesh1%nsurfaces), &
                mesh%surfaces(:mesh%nsurfaces) )

!   volumes

    mesh%nvolumes = mesh1%nvolumes

    allocate( mesh%volumes(max(mesh%nvolumes,MAXVOLUMES)) )

    call copy ( mesh1%volumes(:mesh1%nvolumes), &
                mesh%volumes(:mesh%nvolumes) )

!   blocks

    mesh%nblocks = mesh1%nblocks

    allocate( mesh%blocks(max(mesh%nblocks,MAXBLOCKS)) )

    call copy ( mesh1%blocks(:mesh1%nblocks), &
                mesh%blocks(:mesh%nblocks) )

!   objects

    mesh%nobjects = mesh1%nobjects

    allocate( mesh%objects(max(mesh%nobjects,MAXOBJECTS)) )

    call copy ( mesh1%objects(:mesh1%nobjects), &
                mesh%objects(:mesh%nobjects) )

!   nodesets

    mesh%nnodesets = mesh1%nnodesets

    allocate ( mesh%nodesets(max(mesh%nnodesets,MAXNODESETS)) )
    do nodeset = 1, mesh%nnodesets
      mesh%nodesets(nodeset)%a = mesh1%nodesets(nodeset)%a
    end do

!   elementsets

    mesh%nelementsets = mesh1%nelementsets

    allocate ( mesh%elementsets(max(mesh%nelementsets,MAXELEMENTSETS)) )
    do elementset = 1, mesh%nelementsets
      mesh%elementsets(elementset)%nodes_created = &
                      mesh1%elementsets(elementset)%nodes_created
      mesh%elementsets(elementset)%nelem = &
                      mesh1%elementsets(elementset)%nelem
      mesh%elementsets(elementset)%nnodes = &
                      mesh1%elementsets(elementset)%nnodes
      mesh%elementsets(elementset)%nodes = &
                   mesh1%elementsets(elementset)%nodes
      mesh%elementsets(elementset)%grpnumel = &
                      mesh1%elementsets(elementset)%grpnumel
      mesh%elementsets(elementset)%elements = &
                     mesh1%elementsets(elementset)%elements
    end do

!   sblocks

    call copy ( mesh1%sblocks, mesh%sblocks )

!   nodblocks

    call copy ( mesh1%nodblocks, mesh%nodblocks )

!   basic mesh has been generated

    mesh%meshgen = .true.

  end subroutine copy_mesh_basic


! Delete geometry

  subroutine delete_geometry ( geometry )

    type(geometry_t), intent(out) :: geometry

  end subroutine delete_geometry


! Delete object

  subroutine delete_object ( object )

    type(object_t), intent(out) :: object

  end subroutine delete_object


! Delete block

  subroutine delete_block ( block )

    type(block_t), intent(out) :: block

  end subroutine delete_block


! Delete all blocks from mesh

  subroutine delete_blocks ( mesh )

    type(mesh_t), intent(inout) :: mesh

    integer :: block

    do block = 1, mesh%nblocks
      call delete_block ( mesh%blocks(block) )
    end do

    mesh%nblocks = 0

  end subroutine delete_blocks


! Delete sblocks from mesh

  subroutine delete_sblocks ( sblocks )

    type(sblocks_t), intent(out) :: sblocks

  end subroutine delete_sblocks


! Delete nodblocks from mesh

  subroutine delete_nodblocks ( nodblocks )

    type(nodblocks_t), intent(out) :: nodblocks

  end subroutine delete_nodblocks


! Delete elementset

  subroutine delete_elementset ( elementset )

    type(elementset_t), intent(out) :: elementset

  end subroutine delete_elementset


! Delete mesh

  subroutine delete_mesh ( mesh1, mesh2, mesh3, mesh4, mesh5, mesh6 )

    type(mesh_t), intent(inout) :: mesh1
    type(mesh_t), optional, intent(inout) :: mesh2, mesh3, mesh4, mesh5, mesh6

    call delete_single_mesh(mesh1)
    if ( present(mesh2) ) call delete_single_mesh(mesh2)
    if ( present(mesh3) ) call delete_single_mesh(mesh3)
    if ( present(mesh4) ) call delete_single_mesh(mesh4)
    if ( present(mesh5) ) call delete_single_mesh(mesh5)
    if ( present(mesh6) ) call delete_single_mesh(mesh6)

  end subroutine delete_mesh


! Delete single mesh

  subroutine delete_single_mesh ( mesh )

    type(mesh_t), intent(inout) :: mesh

    if ( .not. mesh%meshgen ) then
      write(*,'(/a/)') &
        'Error delete_mesh: no mesh present in structure mesh '
      stop
    end if

!   deallocate all allocatables and intialize to defaults
    call deall ( mesh )

  contains

    subroutine deall ( mesh )
      type(mesh_t), intent(out) :: mesh
    end subroutine deall

  end subroutine delete_single_mesh


! deallocate memory for arrays in element (core only)

  subroutine deallocate_element_arrays_core ( element )

    type(element_t), intent(inout) :: element

    deallocate( element%sidvert )
    deallocate( element%sidnod )
    deallocate( element%internnod )
    deallocate( element%sidinternnod )
    deallocate( element%edgevert )
    deallocate( element%edgenod )
    deallocate( element%facevert )
    deallocate( element%facenod )
    deallocate( element%xc )

  end subroutine deallocate_element_arrays_core


! help routine to avoid large compilation times in ifort

  subroutine copy_g ( geometry1, geometry2 )

    type(geometry_t), intent(in) :: geometry1
    type(geometry_t), intent(out) :: geometry2

    geometry2 = geometry1

  end subroutine copy_g

end module mesh_m
