
! Copyright (C) 2004-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Type definitions for the problem

module problem_defs_m

  use kind_defs_m
  use array_defs_m
  use mesh_m, only: mesh_t
  use sparse_m
  use system_defs_m
  use set_optional_m

  implicit none


! type definition for essential boundary conditions in points

  type esspoints_t

    integer :: physq = 0      ! if >0 it gives the physical quantity to be set
    integer :: layer = 0      ! if >0 it gives the layer to be set
    integer :: point = 0      ! point number

    integer, allocatable, dimension(:) :: points   ! multiple point numbers

!   degfd indicates which degrees of freedom are essential bc:
!   degfd(i) >0 => degfd i is set
!   for example degfd=(/1,0,1/) sets the first and third degree of freedom
!   to be essential.
!   Notes: - for physq >0, degfd refers to the physical quantity
!            for physq =0, degfd refers to all degrees in the nodes
!            for layer >0, degfd refers to the layer only
!            for layer =0, degfd refers to all degrees in the nodes
!          - if size(degfd) = 0, all degrees are set.
    integer, allocatable, dimension(:) :: degfd

  end type esspoints_t


! type definition for essential boundary conditions on curves

  type esscurve_t

    integer :: physq = 0      ! if >0 it gives the physical quantity to be set
    integer :: layer = 0      ! if >0 it gives the layer to be set
    integer :: first = 1      ! first curve number
    integer :: last  = 0      ! last curve number
    integer :: step  = 0      ! essential bc for nodes on the curves:
                              !   step=0 all nodes
                              !   step>0 nodes 1, 1+step, 1+2*step, ...
                              !   step<0 except nodes 1, 1-step, 1-2*step, ...
                              ! note that this is applied for all curves
                              ! _separately_ and first to last is not treated
                              ! as a single curve.
    integer :: exclude  = 0   ! exclude nodes:
                              !   exclude=0 no excludes, use all nodes
                              !   exclude=1 exclude first point of first curve
                              !   exclude=2 exclude last point of last curve
                              !   exclude=3 exclude first and last point

    integer, allocatable, dimension(:) :: curves   ! multiple curve numbers

!   exclude points:
!   for example excludepoints=(/1,3/) excludes points P1 and P3 to be essential.
!   Note that only the points that are actually on the essential boundary
!   curves are excluded.
    integer, allocatable, dimension(:) :: excludepoints

!   exclude curves:
!   for example excludecurves=(/1,3/) excludes nodes on curves C1 and C3
!   to be essential. Note that only the nodes that are actually on the
!   essential boundary curves are excluded.
    integer, allocatable, dimension(:) :: excludecurves

!   exclude surfaces:
!   for example excludesurfaces=(/1,3/) excludes nodes on surfaces S1 and S3
!   to be essential. Note that only the nodes that are actually on the
!   essential boundary curves are excluded.
    integer, allocatable, dimension(:) :: excludesurfaces

!   degfd indicates which degrees of freedom are essential bc:
!   degfd(i) >0 => degfd i is set
!   for example degfd=(/1,0,1/) sets the first and third degree of freedom
!   to be essential.
!   Notes: - for physq >0, degfd refers to the physical quantity
!            for physq =0, degfd refers to all degrees in the nodes
!            for layer >0, degfd refers to the layer only
!            for layer =0, degfd refers to all degrees in the nodes
!          - if size(degfd) = 0, all degrees are set.
    integer, allocatable, dimension(:) :: degfd

  end type esscurve_t


! type definition for essential boundary conditions on surfaces

  type esssurface_t

    integer :: physq = 0      ! if >0 it gives the physical quantity to be set
    integer :: layer = 0      ! if >0 it gives the layer to be set
    integer :: first = 1      ! first surface number
    integer :: last  = 0      ! last surface number

    integer, allocatable, dimension(:) :: surfaces   ! multiple surface numbers

!   exclude points:
!   for example excludepoints=(/1,3/) excludes points P1 and P3 to be essential.
!   Note that only the points that are actually on the essential boundary
!   surfaces are excluded.
    integer, allocatable, dimension(:) :: excludepoints

!   exclude curves:
!   for example excludecurves=(/1,3/) excludes nodes on curves C1 and C3
!   to be essential. Note that only the nodes that are actually on the surfaces
!   are excluded.
    integer, allocatable, dimension(:) :: excludecurves

!   exclude surfaces:
!   for example excludesurfaces=(/1,3/) excludes nodes on surfaces S1 and S3
!   to be essential. Note that only the nodes that are actually on the surfaces
!   having essential boundary conditions are excluded.
    integer, allocatable, dimension(:) :: excludesurfaces

!   degfd indicates which degrees of freedom are essential bc:
!   degfd(i) >0 => degfd i is set
!   for example degfd=(/1,0,1/) sets the first and third degree of freedom
!   to be essential.
!   Notes: - for physq >0, degfd refers to the physical quantity
!            for physq =0, degfd refers to all degrees in the nodes
!            for layer >0, degfd refers to the layer only
!            for layer =0, degfd refers to all degrees in the nodes
!          - if size(degfd) = 0, all degrees are set.
    integer, allocatable, dimension(:) :: degfd

  end type esssurface_t


! type definition for essential boundary conditions in elements

  type esselements_t

    integer :: physq = 0      ! if >0 it gives the physical quantity to be set
    integer :: layer = 0      ! if >0 it gives the layer to be set
    integer :: elgrp = 0      ! element group number
    integer :: elem  = 0      ! element number
    integer :: node  = 0      ! local node number in the element

!   degfd indicates which degrees of freedom are essential bc:
!   degfd(i) >0 => degfd i is set
!   for example degfd=(/1,0,1/) sets the first and third degree of freedom
!   to be essential.
!   Notes: - for physq >0, degfd refers to the physical quantity
!            for physq =0, degfd refers to all degrees in the nodes
!            for layer >0, degfd refers to the layer only
!            for layer =0, degfd refers to all degrees in the nodes
!          - if size(degfd) = 0, all degrees are set.
    integer, allocatable, dimension(:) :: degfd

  end type esselements_t


! type definition for essential boundary conditions on element groups

  type essgroups_t

    integer :: physq = 0      ! if >0 it gives the physical quantity to be set
    integer :: layer = 0      ! if >0 it gives the layer to be set
    integer :: first = 1      ! first group number
    integer :: last  = 0      ! last group number

    integer, allocatable, dimension(:) :: elgroups   ! multiple group numbers

!   exclude points:
!   for example excludepoints=(/1,3/) excludes points P1 and P3 to be essential.
!   Note that only the points that are actually in the essential groups are
!   excluded.
    integer, allocatable, dimension(:) :: excludepoints

!   exclude curves:
!   for example excludecurves=(/1,3/) excludes nodes on curves C1 and C3
!   to be essential. Note that only the nodes that are actually in the groups
!   are excluded.
    integer, allocatable, dimension(:) :: excludecurves

!   exclude surfaces:
!   for example excludesurfaces=(/1,3/) excludes nodes on surfaces S1 and S3
!   to be essential. Note that only the nodes that are actually in the groups
!   are excluded.
    integer, allocatable, dimension(:) :: excludesurfaces

!   degfd indicates which degrees of freedom are essential bc:
!   degfd(i) >0 => degfd i is set
!   for example degfd=(/1,0,1/) sets the first and third degree of freedom
!   to be essential.
!   Notes: - for physq >0, degfd refers to the physical quantity
!            for physq =0, degfd refers to all degrees in the nodes
!            for layer >0, degfd refers to the layer only
!            for layer =0, degfd refers to all degrees in the nodes
!          - if size(degfd) = 0, all degrees are set.
    integer, allocatable, dimension(:) :: degfd

  end type essgroups_t


! type definition for essential boundary conditions on nodesets

  type essnodeset_t

    integer :: physq = 0      ! if >0 it gives the physical quantity to be set
    integer :: layer = 0      ! if >0 it gives the layer to be set
    integer :: first = 1      ! first nodeset number
    integer :: last  = 0      ! last nodeset number

    integer, allocatable, dimension(:) :: nodesets   ! multiple nodesets numbers

!   degfd indicates which degrees of freedom are essential bc:
!   degfd(i) >0 => degfd i is set
!   for example degfd=(/1,0,1/) sets the first and third degree of freedom
!   to be essential.
!   Notes: - for physq >0, degfd refers to the physical quantity
!            for physq =0, degfd refers to all degrees in the nodes
!            for layer >0, degfd refers to the layer only
!            for layer =0, degfd refers to all degrees in the nodes
!          - if size(degfd) = 0, all degrees are set.
    integer, allocatable, dimension(:) :: degfd

!   exclude points:
!   for example excludepoints=(/1,3/) excludes points P1 and P3 to be essential.
!   Note that only the points that are actually in the essential nodesets are
!   excluded.
    integer, allocatable, dimension(:) :: excludepoints

!   exclude curves:
!   for example excludecurves=(/1,3/) excludes nodes on curves C1 and C3
!   to be essential. Note that only the nodes that are actually in the nodesets
!   are excluded.
    integer, allocatable, dimension(:) :: excludecurves

!   exclude surfaces:
!   for example excludesurfaces=(/1,3/) excludes nodes on surfaces S1 and S3
!   to be essential. Note that only the nodes that are actually in the nodesets
!   are excluded.
    integer, allocatable, dimension(:) :: excludesurfaces

  end type essnodeset_t


! type definition for constraints

  type constraint_t

    integer :: typeconstraint = 0    ! type of the constraint:
                              !   typeconstraint=0 constraint is not active
                              !   typeconstraint=1 distributed constraint
                              !   typeconstraint=2 global constraint

    integer :: physq1 = 0     ! if >0 it gives the physical quantity to be
                              ! put under a constraint on geometry1,
                              ! elementset1 or on the object.
                              ! If physq1=0 all unknowns are assumed to be
                              ! under constraint.

    integer :: physq2 = 0     ! if >0 it gives the physical quantity to be
                              ! put under a constraint on geometry2,
                              ! elementset2, on the "second side" of the object
                              ! or on object2.
                              ! If physq2=0 all unknowns are assumed to be
                              ! under constraint.

    integer :: layer1 = 0     ! if >0 it gives the layer to be put under a
                              ! constraint on geometry1, elementset1 or on
                              ! the object.
                              ! If layer1=0 all layers are assumed to be
                              ! under constraint.

    integer :: layer2 = 0     ! if >0 it gives the layer to be put under a
                              ! constraint on geometry2, elementset2,
                              ! on the "second side" of the object or
                              ! on object2.
                              ! If layer2=0 all layers are assumed to be
                              ! under constraint.

    integer :: errorlayer = 0 ! if layer1 > 0 or layer2 > 0 check whether
                              ! the layer is present in all nodes involved in
                              ! the constraint. Possibilities:
                              !  0: Stop the program if the check fails
                              !  1: Give a warning if the check fails
                              !  2: Do nothing (no checking, no message).

    integer :: lmnumdegfd = 1 ! The number of Lagrangian multiplier degrees of
                              ! freedom in the nodes of geometries is given by
                              ! 0: either elnumdegfd/nodenumdegfd or equal to
                              !    the number of degrees of freedom of the
                              !    constraint quantity.
                              ! 1: Apply additional rules in setting the number
                              !    of Lagrangian multiplier degrees of freedom
                              !    in the nodes of geometries in case of
                              !    layers (layer1>0 and/or layer2>0), such as
                              !    taking the number of degrees in a single
                              !    layer only and zero degrees if layer is not
                              !    present.
                              ! Note, that this only applies to geometries
                              ! and not to objects. For objects elnumdegfd or
                              ! nodenumdegfd must be used.

    integer :: typegeometry = 0  ! type of the geometry:
                              !   typegeometry=0 no geometry
                              !   typegeometry=1 point
                              !   typegeometry=2 curve
                              !   typegeometry=3 surface
                              !   typegeometry=4 volume
                              !   typegeometry=5 nodeset

    integer :: geometry1 = 0 ! geometry number to put constraint on. Lagrange
                             ! multipliers are on this geometry

    integer :: geometry2 = 0 ! geometry number on the other side of connected
                             ! geometries

    integer :: elementset1 = 0 ! elementset number to put constraint on.
                               ! Lagrange multipliers are on this elementset

    integer :: elementset2 = 0 ! elementset number on the other side of a
                               ! connected elementsets

    integer :: object = 0    ! object number to put constraint on. Lagrange
                             ! multipliers are defined in the nodes of this
                             ! object.

    integer :: object2 = 0   ! object number on the other side of connected
                             ! objects

    logical :: full2 = .false.  ! If full2=.true. then in the assembly
                                ! of the constraint equations it is
                                ! assumed that Lagrange multipliers
                                ! (constraints) are connected to _all_
                                ! degrees of freedom on geometry2.

    integer :: discretization = 0   ! type of discretization of the constraint
                                    ! and thus also the Lagrange multipliers
                                    ! for a distributed constraint:
                                    !   discretization=0 weak form
                                    !   discretization=1 collocation using
                                    !     nodes on geometry1 or object
                                    ! for a global constraint it specifies
                                    ! how the constraint is built:
                                    !   discretization=0 element by element
                                    !   discretization=1 node by node

    integer :: step = 0       ! nodes on geometry1 for typegeometry=2 (curve)
                              ! and discretization=1 collocation:
                              !   step=0 all nodes
                              !   step>0 nodes 1, 1+step, 1+2*step, ...
                              !   step<0 except nodes 1, 1-step, 1-2*step, ...

    integer :: exclude = 0    ! exclude nodes on geometry1 for
                              ! typegeometry=2 (curve) and for
                              ! discretization=1 (collocation):
                              !   exclude=0 no excludes, use all nodes
                              !   exclude=1 exclude first point of curve
                              !   exclude=2 exclude last point of curve
                              !   exclude=3 exclude first and last point

!   if .true. the diagonal block of the matrix is included in the
!   system matrix structure for the constraint part, including the additional
!   unknowns part.
!   Note, that it is not the full block, but the matrix elements that would
!   naturally be filled with zero when building the contraint point for point
!   (collocation) or element by element (weak).
    logical :: diagonal_block = .false.

!   if .true. the diagonal block of the matrix regarding the
!   the additional unknowns is included in the system matrix structure.
    logical :: diagonal_block_addunknowns = .false.

!   exclude points for a constraint on a geometry and discretization=1
!   (collocation). For example excludepoints=(/1,3/) excludes nodes in
!   points P1 and P3 from the constraint. Note that only the nodes that are
!   actually on the geometry are excluded. Note that the constraints are on
!   the first geometry (geometry1).
    integer, allocatable, dimension(:) :: excludepoints

!   exclude curves for a constraint on a geometry and discretization=1
!   (collocation). For example excludecurves=(/1,3/) excludes nodes on
!   curves C1 and C3 from the constraint. Note that only the nodes that are
!   actually on the geometry are excluded. Note that the constraints are on
!   the first geometry (geometry1).
    integer, allocatable, dimension(:) :: excludecurves

!   exclude surfaces for a constraint on a geometry and discretization=1
!   (collocation). For example excludesurfaces=(/1,3/) excludes nodes on
!   surfaces S1 and S3 from the constraint. Note that only the nodes that are
!   actually on the geometry are excluded. Note that the constraints are on
!   the first geometry (geometry1).
    integer, allocatable, dimension(:) :: excludesurfaces

!   element degrees of freedom on a geometry or object for the Lagrange
!   multiplier for discretization=0 (weak form)
!   elnumdegfd(node)
    integer, allocatable, dimension(:) :: elnumdegfd

!   number of nodal degrees of freedom on a geometry or object
!   for the Lagrange multiplier for discretization=1 (collocation)
    integer :: nodenumdegfd = 0

    integer :: nglobalc = 0   ! number of global constraints for typecontraint=2

    integer :: naddunknowns = 0   ! number of additional unknowns

!   array of length 2 storing the number of degrees of freedom for
!   the global Lagrange multiplier (accumulated).
!   globnumdegfd(1) = previous pointer in system array
!   number of degrees of freedom is
!   globnumdegfd(2) - globnumdegfd(1) == nglobalc
!   Note: this array is always filled, even if nglobalc=0
    integer, dimension(2) :: globnumdegfd = 0

!   array of length geometry%nnodes+1 storing the number of degrees of
!   freedom for the Lagrange multiplier in each nodal point (accumulated).
!   nodnumdegfd(1) = previous pointer in system array
!   number of degrees of freedom for nodal point n is
!   nodnumdegfd(n+1) - nodnumdegfd(n)
!   Note: this array is only available for a distributed constraint
!   (typeconstraint=1)
    integer, allocatable, dimension(:) :: nodnumdegfd

!   maximum number of degrees of freedom for the Lagrange multiplier
!   in the nodes
    integer :: maxnumdegfd = 0

!   array of length 2 storing the number degrees of freedom for
!   the additional constraint unknowns (accumulated).
!   addnumdegfd(1) = previous pointer in system array
!   number of degrees of freedom is
!   addnumdegfd(2) - addnumdegfd(1) == naddunknowns
!   Note: this array is always filled, even if naddunknowns=0
    integer, dimension(2) :: addnumdegfd = 0

  end type constraint_t


! type definition for connections
!
! A connection adds a local coupling between two parts of the mesh by adding
! "finite stiffness" between the parts.
! Supported are:
!  a. Points, nodesets, geometries (curves, surfaces, volumes) and elementsets.
!     In this case the elements (weak) or nodes (collocation) of the
!     are directly coupled. The assembling of the system is elementwise (weak)
!     or nodal point wise (collocation). The coupling of elements/nodes is
!     according to the sequence number defined in the geometries/elementsets.
!     The two connected points/nodesets/geometries/elementsets should
!     match regarding number of elements (weak) or number of nodes
!     (collocation). The element types of corresponding elements do not
!     have to be the same. However, points must be connected to points, curves
!     to curves, ... etc. with two exceptions:
!        1. it is possible to connect geometries to elementsets as long as
!           they match regarding number of elements (weak) or number of
!           nodes (collocation)
!        2. it is possible to connect geometries to nodesets as long as they
!           match regarding number of number of nodes (collocation).
!  b. Objects. The connection is made between the mesh elements intersected by
!     a double sided object or two separate objects. The assembling of the
!     system is always pointwise. In case of a weak connection the points are
!     the integration points in the object and for the collocated connection
!     the points are the nodes in the object.
!  c. Objects (single sided) connected to geometries, elementsets or nodesets.
!     The connection is made between the mesh elements intersected by the object
!     and the geometry, elementset or nodeset. The assembling of the
!     system is always pointwise. In case of a weak connection the points are
!     the integration points in the object and for the collocated connection
!     the points are the nodes in the object. For nodesets on the second side
!     only collocation is possible.

! Note, that a basic difference between a. (geometries/elementsets) and b.
! objects is that in the first case the connection is (usually)
! between element edges, whereas in the latter case the connection is always
! between mesh ("volume") elements. Case c. is a mix between the two.
!
! For "rigid" couplings constraints should be used.

  type connection_t

    integer :: physq1 = 0     ! if >0 it gives the physical quantity to be
                              ! connected on geometry1, elementset1 or object.
                              ! If physq1=0 all unknowns are assumed to be
                              ! connected.

    integer :: physq2 = 0     ! if >0 it gives the physical quantity to be
                              ! connected on geometry2, elementset2, the
                              ! "second side" of object or on object2.
                              ! If physq2=0 all unknowns are assumed to be
                              ! connected.

    integer :: layer1 = 0     ! if >0 it gives the layer to be connected on
                              ! geometry1, elementset1 or object.
                              ! If layer1=0 all layers are assumed to be
                              ! connected.

    integer :: layer2 = 0     ! if >0 it gives the layer to be connected
                              ! on geometry2, elementset2, the "second side" of
                              ! object or on object2.
                              ! If layer2=0 all layers are assumed to be
                              ! connected.

    integer :: typegeometry = 0  ! type of the geometry:
                              !   typegeometry=0 no geometry
                              !   typegeometry=1 point
                              !   typegeometry=2 curve
                              !   typegeometry=3 surface
                              !   typegeometry=4 volume
                              !   typegeometry=5 nodeset

    integer :: typegeometry2 = 0  ! type of the geometry on side 2:
                              !   typegeometry2=0 no geometry
                              !   typegeometry2=1 point
                              !   typegeometry2=2 curve
                              !   typegeometry2=3 surface
                              !   typegeometry2=4 volume
                              !   typegeometry2=5 nodeset

    integer :: geometry1 = 0 ! first geometry number

    integer :: geometry2 = 0 ! second geometry number

    integer :: elementset1 = 0 ! first elementset

    integer :: elementset2 = 0 ! second elementset

    integer :: object = 0    ! object number.

    integer :: object2 = 0   ! object number 2.

    integer :: discretization = 0   ! type of discretization of the connection
                                    !   discretization=0 weak form
                                    !   discretization=1 collocation using
                                    !     nodes

  end type connection_t


! type definition for local transformations
!
! Local transformations are transformations of the degrees of freedom in nodes.
! They are defined node for node, hence they are local. Transformations are
! defined as follows:
!
!    u = A u'
!
! where u are the standard degrees and u' are the transformed ones. The
! non-singular square matrix A must be defined by the user node for node if
! typetransformation=1 or a single matrix which is the same for all nodes
! if typetransformation=2.
!
! Note, that
!  a. Multiple defined transformations should not overlap. i.e. a degree
!     of freedom can be involved in a single transformation only.
!  b. Transformations do not affect the definition of element matrices.
!     These should be written using the standard degrees of freedom. The
!     transformations are applied to element matrices and vectors before
!     being assembled to the system matrix.
!  c. The definition of Dirichlet conditions is with respect to the transformed
!     degrees of freedom u'. Thus setting degfd=[1] to Dirichlet means the first
!     transformed degree of freedoms is set to Dirichlet.

  type transformation_t

    logical :: build = .false. ! has this transformation been used for the
                               ! build of the transformation matrix (only
                               ! needed for problem).

    logical :: orthogonal = .true.  ! Orthogonality of the matrix A:
                                    ! orthogonal=.true.: A is orthogonal,
                                    !    meaning A A^T = I or A^-1=A^T.
                                    ! orthogonal=.false.: A is non-orthogonal
                                    !    (general but invertable). Note, that
                                    !    the transformation from standard
                                    !    degrees to transformed degrees
                                    !    u' = A^-1 u is not yet supported
                                    !    for this case.

    integer :: typetransformation = 0    ! type of the transformation:
                   !   typetransformation=0 transformation is not active
                   !   typetransformation=1 distributed transformation, meaning
                   !        that the transformation matrix A must be supplied
                   !        node for node.
                   !   typetransformation=2 global transformation, meaning that
                   !        that the transformation matrix A is a constant for
                   !        all nodes involved in the transformation.

    integer :: normalvector = 0  ! how to build the transformation for
           ! typetransformation=1:
           !   normalvector=-1 use u1=-n on the geometry as the x'-axis
           !   normalvector=0 matrix A is fully supplied by the user.
           !   normalvector=1 use u1=n on the geometry as the x'-axis
           ! Note: normalvector/=0 can only be used for typegeometry=2 (curve)
           ! and 3 (surface).
           ! For normalvector=-1 or 1 in 2D the transformation matrix is
           ! defined as follows: the vector u1 is the vector pointing in the
           ! direction of new x'-axis. The vector u2=[-u1(2),u1(1)] is
           ! pointing in the direction of the new y'-axis. The vectors u1
           ! and u2=[-u1(2),u1(1)] are first normalized to 1: u1=u1/|u1|,
           ! u1=u2/|u2|. Then the matrix A becomes the rotation matrix
           ! A = [ u1 u2 ]
           ! This is for 2D only and only for transformation of vectors.
           ! For 3D, see the description below of the vector v2.

    real(dp), dimension(3) :: v2 = [ 0._dp, 0._dp, 1._dp ]
              ! Vector used for normalvector=-1 or 1 in 3D to define the
              ! transformation. If u1 is the unit vector in the x'-axis (normal
              ! direction), then the vector u2 pointing in the y'-direction
              ! becomes u2 = u1 x v2 and the vector u3 pointing in the
              ! z'-direction
              !    u3 = u1 x u2
              ! The vectors u1, u2 and u3 are first normalized to 1: u1=u1/|u1|,
              ! u1=u2/|u2|, u3=u3/|u3|.
              ! Then the matrix A becomes the rotation matrix
              !    A = [ u1 u2 u3 ]
              ! This is for 3D only and only for transformation of vectors.

    integer :: physq = 0      ! if >0 it gives the physical quantity to be
                              ! be transformed on the geometry.
                              ! If physq=0 all unknowns are assumed to be
                              ! transformed.

    integer :: layer = 0      ! if >0 it gives the layer to be transformed
                              ! on geometry.
                              ! If layer=0 all unknowns are assumed to be
                              ! transformed.

    integer :: errorlayer = 0 ! if layer > 0 check whether
                              ! the layer is present in all nodes involved in
                              ! the transformation. Possibilities:
                              !  0: Stop the program if the check fails
                              !  1: Give a warning if the check fails
                              !  2: Do nothing (no checking, no message).

    integer :: typegeometry = 0  ! type of the geometry:
                              !   typegeometry=0 no geometry
                              !   typegeometry=1 point
                              !   typegeometry=2 curve
                              !   typegeometry=3 surface
                              !   typegeometry=4 volume
                              !   typegeometry=5 nodeset

    integer :: geometry = 0 ! geometry number to be transformed.

    integer :: step = 0       ! nodes on geometry for typegeometry=2 (curve)
                              !   step=0 all nodes
                              !   step>0 nodes 1, 1+step, 1+2*step, ...
                              !   step<0 except nodes 1, 1-step, 1-2*step, ...

    integer :: exclude = 0    ! exclude nodes on geometry for
                              ! typegeometry=2 (curve):
                              !   exclude=0 no excludes, use all nodes
                              !   exclude=1 exclude first point of curve
                              !   exclude=2 exclude last point of curve
                              !   exclude=3 exclude first and last point

!   exclude points for the transformation.
!   For example excludepoints=(/1,3/) excludes nodes in
!   points P1 and P3 from the transformation. Note that only the nodes that are
!   actually on the geometry are excluded.
    integer, allocatable, dimension(:) :: excludepoints

!   exclude curves for the transformation.
!   For example excludecurves=(/1,3/) excludes nodes on
!   curves C1 and C3 from the transformation. Note that only the nodes that are
!   actually on the geometry are excluded.
    integer, allocatable, dimension(:) :: excludecurves

!   exclude surfaces for the transformation.
!   For example excludesurfaces=(/1,3/) excludes nodes on
!   surfaces S1 and S3 from the transformation. Note that only the nodes that
!   are actually on the geometry are excluded.
    integer, allocatable, dimension(:) :: excludesurfaces

!   The transformation matrix A for a global transformation
!   (=typetransformation=2).
    real(dp), allocatable, dimension(:,:) :: Amat_global

  end type transformation_t


! type definition for dependencies
!
! A dependency adds a coupling between two parts of the mesh by identifying
! the degrees of freedom on the first part with (a linear combination) of the
! degrees of freedom in the second part. The first part contains the dependent
! degrees of freedom, which are eliminated from the system of equations.
! In summary:
!
!    u1 = C u2 + D u_add + e
!
! where u1 are the dependent degrees of freedom, C is a matrix, u2 are the
! degrees of freedom on the second part, D is a matrix, u_add are the
! additional degrees of freedom connected to this dependency and e is a
! vector not dependent on degrees of freedom.

! If the second part (with degrees u2) is absent, additional degrees
! of freedom are the only remaining unknown degrees. If also the additional
! degrees are absent, the dependency is equivalent to imposed Dirichlet
! conditions u1=e.
!
! After assembling all dependencies, the result is a dependency between the
! dependent degrees of freedom (u_d) and the remaining degrees of freedom (u_r),
! written as follows:
!
!   u_d = E u_r + g
!
! where E is a (sparse) matrix and g is a constant vector. Due to the practical
! implementation, the dependent degrees of freedom (the first part) need to
! be defined essential. As a result the dependent degrees u_d are not a
! separate identified set of degrees in tfem. Since u_d becomes part of the
! essential degrees, the matrix E is extended to include all essential degrees
! with proper zero rows for non-dependent ones.
!
! For defining the degrees of freedom in u1, u2, the matrices C and D
! and vector e, points, nodesets, geometries (curves, surfaces, volumes),
! elementsets can be used. Additionally, for the second part objects can be
! used.
!
! Four types of dependencies are defined:
! 1) nodes
!    In this case, the nodes on both parts are coupled node for node as if the
!    two parts are two finite element meshes directly coupled. The latter case,
!    would mean C=I, but the dependencies are more general with C/=I. C can
!    even be non-square, if needed.
!    The number of nodes in both parts need to be the same. The coupling of
!    nodes is according to the sequence number defined in the geometries or
!    elementsets. The building of the matrices C, D and vector e is also nodal
!    point wise.
! 2) elements (*)
!    In this case, the nodes in part1 are coupled to all the nodes in
!    corresponding elements in part2. The building of the matrices C, D and
!    vector is also element wise.
! 3) full (*)
!    In this case the matrix C is full and all degrees in u1 are coupled to
!    all degrees in u2. The building of the matrices C, D and vector e is
!    nodal point wise.
! 4) object (*)
!    In this case part2 is given by an object (object2). The building of the
!    matrices C, D and vector e is nodal point wise (collocation).
!
! NOTES:
! - The four types are only relevant for the nodal degrees of freedom u2.
!   If only the first part has been defined, u2 is absent and
!   only additional degrees of freedom are the unknown degrees. The type is set
!   to 'nodes' and the matrix D and vector e will be build nodal point wise,
!   in that case.
!
! Two data layouts of the system matrix can be chosen:
! a) fem
!    In this case, it is assumed that the sparse matrix layout of the system
!    matrix is like the fem structure as if the two parts are coupled like
!    conforming meshes. This also means in practice that the degrees in u1
!    should not be involved in contraints and connections.
! b) general (*)
!    In this case, no special structure of the system matrix is assumed.
!    However, temporary storage of roughly the size of the system matrix is
!    needed for building the system matrix.
!
! NOTES:
! - The fem data layout can only be used in the case of "nodes" dependencies
!   without contraints or connections for u1.
! - The general data layout can be used with all types of dependencies and is
!   required for the "elements" and "full" dependencies.
!
! (*) This has not yet been implemented.

  type dependency_t

    logical :: build = .false. ! has this dependency been used for the
                               ! build of the dependency matrix (only
                               ! needed for problem)?

    integer :: typedependency = 0    ! type of the dependency:
                              !   typedependency=0 "nodes" (see intro above)
                              !   typedependency=1 "elements" (see intro above)
                              !   typedependency=2 "full" (see intro above)
                              !   typedependency=3 "object" (see intro above)

    integer :: datalayout = 0  ! Data layout:
                              !   datalayout=0 fem (see intro above)
                              !   datalayout=1 general (see intro above)

    integer :: physq1 = 0     ! if >0 it gives the physical quantity to be
                              ! involved on geometry1 or elementset1.
                              ! If physq1=0 all unknowns are assumed to be
                              ! involved.

    integer :: physq2 = 0     ! if >0 it gives the physical quantity to be
                              ! involved on geometry2 or elementset2.
                              ! If physq2=0 all unknowns are assumed to be
                              ! involved.

    integer :: layer1 = 0     ! if >0 it gives the layer to be involved on
                              ! geometry1 or elementset1.
                              ! If layer1=0 all layers are assumed to be
                              ! involved.

    integer :: layer2 = 0     ! if >0 it gives the layer to be involved
                              ! on geometry2 or elementset2.
                              ! If layer2=0 all layers are assumed to be
                              ! connected.

    integer :: typegeometry1 = 0  ! type of the geometry on part 1:
                              !   typegeometry1=0 no geometry
                              !   typegeometry1=1 point
                              !   typegeometry1=2 curve
                              !   typegeometry1=3 surface
                              !   typegeometry1=4 volume
                              !   typegeometry1=5 nodeset

    integer :: typegeometry2 = 0  ! type of the geometry on part 2:
                              !   typegeometry2=0 no geometry
                              !   typegeometry2=1 point
                              !   typegeometry2=2 curve
                              !   typegeometry2=3 surface
                              !   typegeometry2=4 volume
                              !   typegeometry2=5 nodeset

    integer :: geometry1 = 0 ! first geometry number

    integer :: geometry2 = 0 ! second geometry number

    integer :: elementset1 = 0 ! first elementset

    integer :: elementset2 = 0 ! second elementset

    integer :: object2 = 0 ! second object (there is no object1)

    integer :: step = 0       ! nodes on geometry1 for typegeometry1=2 (curve)
                              ! and typedependency=0 "nodes":
                              !   step=0 all nodes
                              !   step>0 nodes 1, 1+step, 1+2*step, ...
                              !   step<0 except nodes 1, 1-step, 1-2*step, ...

    integer :: exclude = 0    ! exclude nodes on geometry1 for
                              ! typegeometry1=2 (curve) and for
                              ! and typedependency=0 "nodes":
                              !   exclude=0 no excludes, use all nodes
                              !   exclude=1 exclude first point of curve
                              !   exclude=2 exclude last point of curve
                              !   exclude=3 exclude first and last point

!   exclude points for the dependency.
!   For example excludepoints=(/1,3/) excludes nodes in
!   points P1 and P3 from the dependency. Note that only the nodes that are
!   actually on the geometry are excluded.
    integer, allocatable, dimension(:) :: excludepoints

!   exclude curves for the dependency.
!   For example excludecurves=(/1,3/) excludes nodes on
!   curves C1 and C3 from the dependency. Note that only the nodes that are
!   actually on the geometry are excluded.
    integer, allocatable, dimension(:) :: excludecurves

!   exclude surfaces for the dependency.
!   For example excludesurfaces=(/1,3/) excludes nodes on
!   surfaces S1 and S3 from the dependency. Note that only the nodes that
!   are actually on the geometry are excluded.
    integer, allocatable, dimension(:) :: excludesurfaces

    integer :: naddunknowns = 0   ! number of additional unknowns

!   array of length 2 storing the number degrees of freedom for
!   the additional constraint unknowns (accumulated).
!   addnumdegfd(1) = previous pointer in system array
!   number of degrees of freedom is
!   addnumdegfd(2) - addnumdegfd(1) == naddunknowns
!   Note: this array is always filled, even if naddunknowns=0
    integer, dimension(2) :: addnumdegfd = 0

  end type dependency_t


! type definition for input to problem_definition

  type input_probdef_t

    logical :: created = .false.  ! has the structure been created?
    integer :: probnr = 0   ! problem number. If 0 no problem number is used.
    integer :: numesspoints = 0   ! number of points with essential bc
    integer :: numesscurves = 0   ! number of curves with essential bc
    integer :: numesssurfaces = 0 ! number of surfaces with essential bc
    integer :: numesselements = 0 ! number of elements with essential bc
    integer :: numessgroups = 0 ! number of groups with essential bc
    integer :: numessnodesets = 0 ! number of nodesets with essential bc
    integer :: numdependencies = 0 ! number of dependencies
    integer :: numconstraints = 0 ! number of constraints
    integer :: numconnections = 0 ! number of connections
    integer :: numtransformations = 0 ! number of transformations
    integer :: numinactivegroups = 0 ! number of inactive groups
    integer :: nvec = 0     ! number of vectors of special structure
    integer :: nphysq = 0   ! number physical quantities
    integer :: nphysqshifted = 0   ! number physical quantities shifted
    integer :: nelgrp = 0   ! number element groups
    integer :: numlayers = 0   ! number of layer

!   element degrees of freedom of a system vector
!   elementdof(elgrp)%a(node)
    type(int_array_1d_t), allocatable, dimension(:) :: elementdof

!   element degrees of freedom of a vector of special structure
!   vec_elementdof(elgrp)%a(node,vec)
    type(int_array_2d_t), allocatable, dimension(:) :: vec_elementdof

!   physical quantities
!   physq(i) gives the vector number of the physical quantity
    integer, allocatable, dimension(:) :: physq

!   physical quantities shifted to end of the system vector (renumber)
!   physqshifted(i) gives the physical quantity shifted
    integer, allocatable, dimension(:) :: physqshifted

!   masking physq partition. Matrix size: [nphysq,nphysq]
!   Initialized to .true. in create_input_probdef
!   The logical matrix defines the coupling between the physical quantities
!   on element level. This affects
!    1) the storage space reserved in the sparse matrix,
!    2) the assembling of the system in the build_system routine.
!   For example, assume nphysq=2 and the 2x2 logical matrix physqmask is
!   filled as follows
!          [ .true. .true. ]
!          [ .true. .false. ]
!   it means that the physical quantity 2 is not coupled with itself and the
!   corresponding block contains only zeros. By default, space is reserved
!   for this zero block in the sparse matrix. This space needs to be filled
!   with zeros, either on element level or with a separate call of
!   build_system for that partition. However, if the creation of the system
!   matrix is called with the argument usephysqmask=.true. the partitions
!   denoted with .false. are fully discarded from the matrix and the zero block
!   does not need to be filled at all. It the filling with zeros is on element
!   level, build_system needs to be called with usephysqmask=.true. as well.
!   Otherwise, an error message appears of insufficient space in the matrix.
    logical, allocatable, dimension(:,:) :: physqmask

!   layers(i) gives the nodeset defining layer i.
    integer, allocatable, dimension(:) :: layers

!   inactive groups
!   inactivegroups(i) gives the groups that are inactive.
!   An inactive group is basically the same as
!     1 setting the number of degrees of freedom in this group to zero
!       (elnumdegfd=0 and vec_elnumdegfd=0)
!     2 the routines that allow "groups" in the heading such as
!       build_system, derive_vector, ... are called with groups=(/active
!       groups/).
!   A subtle difference is though that it is now `legal' to do item 1. and
!   possible warning messages of inconsistent number of degrees between groups
!   are avoided.
!   NOTE, that all the nodal points are still there!
!   NOTE, that only element groups that are fully disconnected from other
!   groups are officially supported, since it is not allowed to have nodes that
!   have a different number of degrees of freedom in the elements connected to
!   a node. You might set allow_different_numdegfd_in_nodes = .true., but then
!   you enter messy territory and you are on your own.
    integer, allocatable, dimension(:) :: inactivegroups

!   essential boundary conditions for points
    type(esspoints_t), allocatable, dimension(:) :: esspoints

!   essential boundary conditions for curves
    type(esscurve_t), allocatable, dimension(:) :: esscurves

!   essential boundary conditions for surfaces
    type(esssurface_t), allocatable, dimension(:) :: esssurfaces

!   essential boundary conditions for elements
    type(esselements_t), allocatable, dimension(:) :: esselements

!   essential boundary conditions for groups
    type(essgroups_t), allocatable, dimension(:) :: essgroups

!   essential boundary conditions for nodesets
    type(essnodeset_t), allocatable, dimension(:) :: essnodesets

!   dependencies
    type(dependency_t), allocatable, dimension(:) :: dependencies

!   constraints
    type(constraint_t), allocatable, dimension(:) :: constraints

!   connections
    type(connection_t), allocatable, dimension(:) :: connections

!   transformations
    type(transformation_t), allocatable, dimension(:) :: transformations

  end type input_probdef_t



! type definition for problem

  type problem_t

    logical :: created = .false.  ! has the structure been created?
    logical :: orthogonal = .true. ! Orthogonality of the matrix Amat.
    logical :: buildAmat = .false. ! has Amat been build (transformations)?
    integer :: probnr = 0  ! problem number. If 0 no numbers are used.
    integer :: nelgrp = 0  ! number of element groups
    integer :: maxnoddegfd = 0 ! maximum number of degs of freedom in the nodes
    integer :: maxvecnoddegfd = 0 ! maximum number of degs of freedom in
                                  ! the nodes for all the vectors
    integer :: numessnodes = 0 ! number of nodes where essential bc are present
    integer :: numdependencies = 0 ! number of dependencies
    integer :: numconstraints = 0   ! number of constraints
    integer :: numconnections = 0   ! number of connections
    integer :: numtransformations = 0   ! number of transformations
    integer :: nvec  = 0    ! number of vectors of special structure
    integer :: nphysq = 0   ! number physical quantities
    integer :: nphysqshifted = 0   ! number physical quantities shifted

!   The number of degrees of freedom is
!       numdegfd = numnodaldegfd + numdependegfd + numconstrdegfd
!   which consists of a nodal part (from nodnumdegfd), a dependency part
!   (additional unknowns) and a constraint part (Langrangian multipliers
!   and additional unknowns).
!   After renumbering, the degrees of freedom are splitted into a unknown
!   part and a prescribed (essential) part:
!       numdegfd = numundegfd + numessdegfd
    integer :: numnodaldegfd = 0 ! number of nodal degrees of freedom
    integer :: numdependegfd = 0 ! number of dependency degrees of freedom
    integer :: numconstrdegfd = 0 ! number of constraint degrees of freedom
    integer :: numdegfd = 0  ! number of degrees of freedom
    integer :: numundegfd  = 0  ! number of unknown degrees of freedom
    integer :: numessdegfd = 0 ! number of prescribed degrees of freedom
    integer :: numinactivegroups = 0 ! number of inactive groups
    integer :: numlayers = 0 ! number of layers
    integer :: maxnodnumlayers = 0 ! maximum number of layers in the nodes
    integer :: numtransdegfd = 0 ! number of transformed degrees of freedom

!   number of degrees of freedom in each nodal point of an element
!   elnumdegfd(elgrp)%a(node) gives number of degfd in node of element
!   for element group elgrp
    type(int_array_1d_t), allocatable, dimension(:) :: elnumdegfd

!   number of degrees of freedom in each nodal point of an element for vectors
!   of special structure
!   vec_elnumdegfd(elgrp)%a(node,vec)=number of degfd in node of element
!   for element group elgrp for vector vec
    type(int_array_2d_t), allocatable, dimension(:) :: vec_elnumdegfd

!   array of length mesh%nnodes storing the layers per node:
!   nodlayers(nodenr) = bitpattern
!   bitpattern denotes the layers that are present in the node nodenr
!   for example ibset(bitpattern,2) will set layer 3.
!   Note that bitpattern=0 means _no_ layers are present.
!   NOTE: this array will be of zero size if numlayers=0
    integer, allocatable, dimension(:) :: nodlayers

!   array of length mesh%nnodes+1 storing the number of degrees of freedom
!   in each nodal point (accumulated).  nodnumdegfd(1) = 0.
!   number of degrees of freedom for nodal point n is
!   nodnumdegfd(n+1) - nodnumdegfd(n)
    integer, allocatable, dimension(:) :: nodnumdegfd

!   array of length mesh%nnodes+1 x problem%nvec storing the number of degrees
!   of freedom in each nodal point (accumulated) of vectors of special structure
!   vec_nodnumdegfd(1,vec) = 0.
!   number of degrees of freedom for nodal point n is
!   vec_nodnumdegfd(n+1,vec) - vec_nodnumdegfd(n,vec)
!   for vector of special structure vec
    integer, allocatable, dimension(:,:) :: vec_nodnumdegfd

!   number of degrees of freedom for each vector of special structure
!   vec_numdegfd(vec,1) = standard vector defined by vec_nodnumdegfd.
!   vec_numdegfd(vec,2) = vector where the components vec_elnumdegfd(:,vec)
!                      are stored per element.
    integer, allocatable, dimension(:,:) :: vec_numdegfd

!   physical quantities
!   physq(i) gives the vector number of the physical quantity
    integer, allocatable, dimension(:) :: physq

!   physical quantities shifted to end of the system vector (renumber)
!   physqshifted(i) gives the physical quantity shifted
    integer, allocatable, dimension(:) :: physqshifted

!   masking physq partition. Matrix size: [nphysq,nphysq]
!   See input_probdef_t for a description.
    logical, allocatable, dimension(:,:) :: physqmask

!   active and inactive groups
!   activegroups(i) gives the groups that are active.
!   inactivegroups(i) gives the groups that are inactive.
!   An inactive group is basically the same as
!     1 setting the number of degrees of freedom in this group to zero
!       (elnumdegfd=0 and vec_elnumdegfd=0)
!     2 the routines that allow "groups" in the heading such as
!       build_system, derive_vector, ... are called with groups=(/active
!       groups/).
!   NOTE, that all the nodal points are still there!
!   NOTE, that only element groups that are fully disconnected from other
!   groups are officially supported, since it is not allowed to have nodes that
!   have a different number of degrees of freedom in the elements connected to
!   a node. You might set allow_different_numdegfd_in_nodes = .true., but then
!   you enter messy territory and you are on your own.
    integer, allocatable, dimension(:) :: activegroups
    integer, allocatable, dimension(:) :: inactivegroups

!   nodes where essential bc are present
!   essnodes(node,1) = nodenr
!   essnodes(node,2) = bitpattern
!   bitpattern denotes the degrees of freedom that are prescribed
!   for example ibset(bitpattern,2) will set degree 3.
!   Note that bitpattern=0 means _none_ are prescribed.
    integer, allocatable, dimension(:,:) :: essnodes

!   array of dimension (problem%numdegfd,2) containing renumbered degrees
!   of freedom
!   renumbered degfd i has old        degfd number degfdperm(i,1)
!          old degfd i has renumbered degfd number degfdperm(i,2)
!   the renumbering is such that the first problem%numundegfd degrees of freedom
!   are the unknowns and the remaining numessdegfd are the prescribed degrees
!   of freedom. The renumbering is used for the storage of the system matrices
!   and the system vectors.
    integer, allocatable, dimension(:,:) :: degfdperm

!   constraints
    type(constraint_t), allocatable, dimension(:) :: constraints

!   connections
    type(connection_t), allocatable, dimension(:) :: connections

!   transformations
    type(transformation_t), allocatable, dimension(:) :: transformations

!   array of dimension problem%numtransdegfd containing the transformed
!   degrees of freedom
    integer, allocatable, dimension(:) :: degfdtrans

!   transformation matrix A.
!   Only (properly) defined for problem%numtransdegfd > 0.
    type(sparsematrix_t) :: Amat

!   dependencies
    type(dependency_t), allocatable, dimension(:) :: dependencies

  end type problem_t


! interface for generic check subroutine

  interface check
    module procedure check_problem
  end interface check

contains

! check problem

  subroutine check_problem ( problem, name_of_routine, mesh, sloppy )

    type(problem_t), intent(in) :: problem
    character(len=*), intent(in) :: name_of_routine
    type(mesh_t), intent(in), optional :: mesh
!   sloppy? Don't check nelgrp
    logical, intent(in), optional :: sloppy

    logical :: lsloppy

    lsloppy = set_optional ( variable=sloppy, default=.false. )

    if ( .not. problem%created ) then
      write(*,'(/3a/)') &
        'Error ', name_of_routine, ': problem has not been defined '
      stop
    end if

    if ( present(mesh) ) then
      if ( problem%nelgrp /= mesh%nelgrp .and. .not. lsloppy ) then
        write(*,'(/2a/a/)') &
          'Error ', name_of_routine, &
          '  Number of element groups not the same in problem and mesh '
        stop
      end if
    end if

  end subroutine check_problem


! count number of layers up until maxlayer

  elemental function count_layers ( bitpattern, maxlayer ) result(nl)

    integer, intent(in) :: bitpattern
    integer, intent(in) :: maxlayer
    integer :: nl

    integer :: i

    nl = 0

    do i = 1, maxlayer
      if ( btest(bitpattern,i-1) ) nl = nl + 1
    end do

  end function count_layers


! find layer in the given nodes

  subroutine find_layer_in_nodes ( problem, layer, nodes, numl, pnuml, lp )

    type(problem_t), intent(in) :: problem

!   the specified layer
    integer, intent(in) :: layer

!   the specified nodes
    integer, dimension(:), intent(in) :: nodes

!   if layer is present in i^th node:
!      numl(i): the number of layers in node i
!     pnuml(i): the number of layers "below" layer in node i
!   else
!     numl(i) = 0
!     pnuml(i) = 0
!   end if
    integer, dimension(:), intent(inout) :: numl, pnuml

!   if layer is present in i^th node: lp(i) = .true. else .false.
    logical, dimension(:), intent(out), optional :: lp


    logical, dimension(size(nodes)) :: wkl
    integer, dimension(size(nodes)) :: bitp

    bitp = problem%nodlayers(nodes)
    wkl = btest(problem%nodlayers(nodes),layer-1)

    if ( present(lp) ) lp(1:size(nodes)) = wkl

    where ( wkl )
!     layer present
      numl = count_layers ( bitp, problem%numlayers )
    else where
!     layer not present
      numl = 0
    end where

    where ( numl > 1 )
      pnuml = count_layers ( bitp, layer-1 )
    else where
      pnuml = 0
    end where

  end subroutine find_layer_in_nodes

end module problem_defs_m
