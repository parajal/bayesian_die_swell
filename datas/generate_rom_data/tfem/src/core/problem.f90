
! Copyright (C) 2004-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines problem type and associated routines

module problem_m

  use glob_defs_m
  use mesh_m
  use problem_defs_m
  use pos_array_m
  use set_optional_m
  use limits_m, only: MAXNUMESSPOINTS, MAXNUMESSCURVES, MAXNUMESSSURFACES, &
                      MAXNUMESSELEMENTS, MAXNUMESSGROUPS, MAXNUMESSNODESETS, &
                      MAXNUMCONSTRAINTS, MAXNUMWARNINGS, MAXNUMCONNECTIONS, &
                      MAXNUMTRANSFORMATIONS, MAXNUMDEPENDENCIES

  implicit none


! set number of bits

  integer, parameter :: NUMBEROFBITS = 32  ! DO NOT CHANGE!!


! module variables

! Allow number of degrees of freedom in nodes to be different in elements
! connected to the nodes. Setting this to .true. is not a standard
! feature and you are entering messy territory. It might work out for you
! if you know what you are doing, but you are on your own!
! YOU HAVE BEEN WARNED!
  logical :: allow_different_numdegfd_in_nodes = .false.


! interface for generic create subroutine

  interface create
    module procedure create_input_probdef
  end interface create


! interface for generic delete subroutine

  interface delete
    module procedure delete_input_probdef, delete_problem
  end interface delete


contains


! Construct input to problem definition

  subroutine create_input_probdef ( mesh, input_probdef, nvec, nphysq, &
    nphysqshifted, numinactivegroups, numlayers )

    type(mesh_t), intent(in) :: mesh
    type(input_probdef_t), intent(inout) :: input_probdef

!   number of vectors of special structure
!   if present and >0 the array input_probdef%vec_elementdof must be filled
!   by the user before entering the problem definition
    integer, intent(in), optional :: nvec

!   number of physical quantities
!   if present and >0 the array input_probdef%physq must be filled by the
!   user before entering the problem definition
    integer, intent(in), optional :: nphysq

!   number of physical quantities that are shifted to end of the system vector
!   to prevent zeros on the diagonal.
!   if present and >0 the array input_probdef%physqshifted must be filled by the
!   user before entering the problem definition
    integer, intent(in), optional :: nphysqshifted

!   number of inactive groups
!   if present and >0 the array input_probdef%inactivegroups must be filled
!   by the user before entering the problem definition
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
    integer, intent(in), optional :: numinactivegroups

!   number of layers
!   if present and >0 the array input_probdef%layers must be filled
!   by the user before entering the problem definition.
!   For example for numlayers=2:
!      input_probdef%layers = (/ 3, 4 /)
!   specifying that the two layers are defined by the nodesets 3 and 4,
!   respectively.
!   If numlayers > 0 than the degrees in the nodes are defined by the layers.
!   Each layer defines a set of degrees as given by elementdof, but only in
!   the nodes that are present in this layer. A node that belongs to two layers
!   will have double the degrees of freedom as defined by elementdof.
!   If numlayers=0, no layers are defined. This is identical to a single layer
!   defined in all nodes.
    integer, intent(in), optional :: numlayers

    integer :: elgrp, nvc


    call check ( mesh, 'create_input_probdef' )

    if ( input_probdef%created ) then
      write(*,'(/a/)') &
        'Error in create_input_probdef: input_probdef has already been created '
      stop
    end if


    input_probdef%nelgrp = mesh%nelgrp

!   allocate elementdof

    allocate( input_probdef%elementdof(mesh%nelgrp) )
    do elgrp = 1, mesh%nelgrp
      allocate( input_probdef%elementdof(elgrp)%a(mesh%elnumnod(elgrp)) )
      input_probdef%elementdof(elgrp)%a = -1  ! -1 means undefined
    end do


!   vectors of special structure

    if ( present(nvec) ) then
      input_probdef%nvec = max(nvec,0)
    else
      input_probdef%nvec = 0
    end if

    nvc = input_probdef%nvec

!   allocate vec_elementdof

    allocate( input_probdef%vec_elementdof(mesh%nelgrp) )
    do elgrp = 1, mesh%nelgrp
      allocate(&
        input_probdef%vec_elementdof(elgrp)%a(mesh%elnumnod(elgrp),nvc) )
      input_probdef%vec_elementdof(elgrp)%a = 0
    end do

!   physical quantities

    input_probdef%nphysq = 0
    input_probdef%nphysqshifted = 0

    if ( present(nphysq) ) then

      input_probdef%nphysq = max(nphysq,0)

      if ( present(nphysqshifted) ) then
        input_probdef%nphysqshifted = max(0,nphysqshifted)
      end if

    end if

    allocate(input_probdef%physq(input_probdef%nphysq))
    input_probdef%physq = 0

    allocate(input_probdef%physqshifted(input_probdef%nphysqshifted))
    input_probdef%physqshifted = 0

    allocate(&
        input_probdef%physqmask(input_probdef%nphysq,(input_probdef%nphysq)) )
    input_probdef%physqmask = .true.

!   inactive groups

    input_probdef%numinactivegroups = 0

    if ( present(numinactivegroups) ) then

      input_probdef%numinactivegroups = max(numinactivegroups,0)

    end if

    allocate(input_probdef%inactivegroups(input_probdef%numinactivegroups))
    input_probdef%inactivegroups = 0

!   layers

    input_probdef%numlayers = 0

    if ( present(numlayers) ) then

      if ( numlayers > NUMBEROFBITS ) then
        write(*,'(/a,i0/)') &
          'Error in create_input_probdef: number of layers is larger than ', &
          NUMBEROFBITS
        stop
      end if

      input_probdef%numlayers = max(numlayers,0)

    end if

    allocate(input_probdef%layers(input_probdef%numlayers))
    input_probdef%layers = 0

!   allocate esspoints

    allocate( input_probdef%esspoints(MAXNUMESSPOINTS) )

!   allocate esscurves

    allocate( input_probdef%esscurves(MAXNUMESSCURVES) )

!   allocate esssurfaces

    allocate( input_probdef%esssurfaces(MAXNUMESSSURFACES) )

!   allocate esselements

    allocate( input_probdef%esselements(MAXNUMESSELEMENTS) )

!   allocate essgroups

    allocate( input_probdef%essgroups(MAXNUMESSGROUPS) )

!   allocate essnodesets

    allocate( input_probdef%essnodesets(MAXNUMESSNODESETS) )

!   allocate constraints

    allocate( input_probdef%constraints(MAXNUMCONSTRAINTS) )

!   allocate connections

    allocate( input_probdef%connections(MAXNUMCONNECTIONS) )

!   allocate transformations

    allocate( input_probdef%transformations(MAXNUMTRANSFORMATIONS) )

!   allocate dependencies

    allocate( input_probdef%dependencies(MAXNUMDEPENDENCIES) )

    input_probdef%created = .true.

  end subroutine create_input_probdef


! Destruct (single) input to problem definition

  subroutine delete_single_input_probdef ( input_probdef, nr )

    type(input_probdef_t), intent(inout) :: input_probdef
    integer, intent(in) :: nr

    if ( .not. input_probdef%created ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error in delete_input_probdef:', &
        ' input_probdef has not been created and cannot be deleted ', &
        ' input_probdef number in heading = ', nr
      stop
    end if

!   deallocate all allocatables and intialize to defaults
   call deall ( input_probdef )

  contains

    subroutine deall ( input_probdef )
      type(input_probdef_t), intent(out) :: input_probdef
    end subroutine deall

  end subroutine delete_single_input_probdef


! Destruct input to problem definition

  subroutine delete_input_probdef ( input_probdef1, input_probdef2, &
    input_probdef3, input_probdef4, input_probdef5 )

    type(input_probdef_t), intent(inout) :: input_probdef1
    type(input_probdef_t), intent(inout), optional :: input_probdef2, &
      input_probdef3, input_probdef4, input_probdef5

    call delete_single_input_probdef(input_probdef1,1)
    if ( present(input_probdef2) ) &
      call delete_single_input_probdef(input_probdef2,2)
    if ( present(input_probdef3) ) &
      call delete_single_input_probdef(input_probdef3,3)
    if ( present(input_probdef4) ) &
      call delete_single_input_probdef(input_probdef4,4)
    if ( present(input_probdef5) ) &
      call delete_single_input_probdef(input_probdef5,5)

  end subroutine delete_input_probdef


! Definition of essential boundary conditions

  subroutine define_essential ( mesh, input_probdef, physq, layer, degfd, &
    degsfd, points, point, curves, curve1, curve2, step, exclude, group, &
    element, elnode, surfaces, surface1, surface2, elgroups, elgroup1, &
    elgroup2, nodesets, nodeset1, nodeset2, excludepoints, excludecurves, &
    excludesurfaces, num )

    type(mesh_t), intent(in) :: mesh

!   structure to store all definitions for the problem definition
    type(input_probdef_t), intent(inout) :: input_probdef

!   if physq is present physq denotes the physical quantity that must be set
!   otherwise all nodal degrees of freedom are considered.
    integer, intent(in), optional :: physq

!   if layer is present, nodal degrees of freedom to be set are restricted to
!   the given layer otherwise all nodal degrees of freedom are considered.
    integer, intent(in), optional :: layer

!   if degfd is present the array degfd indicates which degrees of freedom are
!   essential:
!       degfd(i) > 0 => degfd i is set
!   otherwise all degrees are set essential.
!   for example degfd=[1,0,1] sets the first and third degree of freedom
!   to be essential. If a degree of freedom is set that does not exist, it is
!   simply ignored. Only the first NUMBEROFBITS degrees of freedom in a node
!   can be set.
!   NOTE: if size(degfd)=0 all degrees are set, as if degfd is not present
!         if degfd=[0] no degrees are set, inactivating the essential bc.
!         if physq is present degfd refers to the physical quantity only
!         if layer is present degfd refers to the layer only
!   NOTE: degfd and degsfd are not allowed to be both present.
    integer, intent(in), dimension(:), optional :: degfd

!   if degsfd is present the array degsfd indicates which degrees of freedom are
!   essential. For example degsfd=[1,3] sets the first and third degree of
!   freedom to be essential. If a degree of freedom is set that does not exist,
!   it is simply ignored. Only the first NUMBEROFBITS degrees of freedom in a
!   node can be set.
!   NOTE: if size(degsfd)=0 all degrees are set, as if degsfd is not present
!         if degsfd=[0] no degrees are set, inactivating the essential bc.
!         if physq is present degsfd refers to the physical quantity only
!         if layer is present degsfd refers to the layer only
!   NOTE: degfd and degsfd are not allowed to be both present.
    integer, intent(in), dimension(:), optional :: degsfd

!   if point is present essential boundary conditions are set in a point
    integer, intent(in), optional :: point

!   if points is present essential boundary conditions are set in the points
!   specified. Note: points can be combines with point in one call.
    integer, intent(in), dimension(:), optional :: points

!   if curve1 is present essential boundary conditions are set along curves
!   curve1 or curve1 to curve2 if curve2 is also present.
!
!   step indicates which nodes on the curves are used:
!       step=0 all nodes
!       step>0 nodes 1, 1+step, 1+2*step, ...
!       step<0 except nodes 1, 1-step, 1-2*step, ...
!   note that this is applied for all curves _separately_ and curve1 to curve2
!   is not treated as a single curve.
!
!   exclude indicates which nodes on the curves are excluded:
!       exclude=0 no excludes, all nodes are used (default)
!       exclude=1 exclude first point of first curve
!       exclude=2 exclude last point of last curve
!       exclude=3 exclude first and last point
!
!   NOTE: step, exclude only work for curves where the local node numbering is
!   in a natural sequence along the curve. This might not be the case for
!   curves constructed from several other curves or curves read from external
!   mesh generators.
    integer, intent(in), optional :: curve1, curve2, step, exclude

!   if curves is present essential boundary conditions are set in the curves
!   specified. Note: curves can be combined with curve1 in one call.
    integer, intent(in), dimension(:), optional :: curves

!   if element is present essential boundary conditions are set in a element
    integer, intent(in), optional :: element

!   if group is present essential boundary conditions for an element in the
!   specified group is set. Default=1
    integer, intent(in), optional :: group

!   if elnode is present essential boundary conditions for an element in the
!   local node elnode is set. Default=1.
    integer, intent(in), optional :: elnode

!   if surface1 is present essential boundary conditions are set along surface
!   surface1 or surface1 to surface2 if surface2 is also present.
    integer, intent(in), optional :: surface1, surface2

!   if surfaces is present essential boundary conditions are set in the surfaces
!   specified. Note: surfaces can be combined with surface1 in one call.
    integer, intent(in), dimension(:), optional :: surfaces

!   if elgroup1 is present essential boundary conditions are set in the element
!   group elgroup1 or elgroup1 to elgroup2 if group2 is also present.
    integer, intent(in), optional :: elgroup1, elgroup2

!   if elgroups is present essential boundary conditions are set in the element
!   groups specified. Note: elgroups can be combines with elgroup1 in one call.
    integer, intent(in), dimension(:), optional :: elgroups

!   if nodeset1 is present essential boundary conditions are set along nodeset
!   nodeset1 or nodeset1 to nodeset2 if nodeset2 is also present.
    integer, intent(in), optional :: nodeset1, nodeset2

!   if nodesets is present essential boundary conditions are set in the nodesets
!   specified. Note: nodesets can be combines with nodeset1 in one call.
    integer, intent(in), dimension(:), optional :: nodesets

!   exclude points:
!   for example excludepoints=(/1,3/) excludes nodes on points P1 and P3
!   to be essential. Can be used together with the parameters curve1, curve2,
!   surface1, surface2, elgroup1, elgroup2, nodeset1, nodeset2. Note that only
!   the nodes that are actually on the curves, surfaces, groups or nodesets
!   are excluded.
    integer, intent(in), dimension(:), optional :: excludepoints

!   exclude curves:
!   for example excludecurves=(/1,3/) excludes nodes on curves C1 and C3
!   to be essential. Can be used together with the parameters curve1, curve2,
!   surface1, surface2, elgroup1, elgroup2, nodeset1, nodeset2. Note that only
!   the nodes that are actually on the curves, surfaces, groups or nodesets are
!   excluded.
    integer, intent(in), dimension(:), optional :: excludecurves

!   exclude surfaces:
!   for example excludesurfaces=(/1,3/) excludes nodes on surfaces S1 and S3
!   to be essential. Can be used together with the parameters curve1, curve2,
!   surface1, surface2, elgroup1, elgroup2, nodeset1, nodeset2. Note that only
!   the nodes that are actually on the surfaces, groups or nodesets are
!   excluded.
    integer, intent(in), dimension(:), optional :: excludesurfaces

!   if present: the sequence number of the essential condition defined,
!   i.e. the index in the arrays esspoints(:), esscurves(:), esssurfaces(:),
!   esselements(:), essgroups(:). Note, that each call of define_essential
!   assigns a sequence number to the boundary condition that is defined,
!   starting with 1 for the first call. Note also that the boundary conditions
!   involving points, curves, surfaces, elements and groups each has their own
!   sequence.
    integer, intent(out), optional :: num


!   This routine can be called multiple times:
!   MAXNUMESSPOINTS times to set essential boundary conditions in points,
!   MAXNUMESSCURVES times to set essential boundary conditions on curves,
!   MAXNUMESSSURFACES times to set essential boundary conditions on surfaces,
!   MAXNUMESSGROUPS times to set essential boundary conditions in groups
!   MAXNUMESSNODESETS times to set essential boundary conditions in nodesets,


    integer :: cg, pnt, elm, elgrp, elnod, i
    integer, allocatable, dimension(:) :: ldegfd


!   some testing

    call check ( mesh, 'define_essential' )

    if ( .not. input_probdef%created ) then
      write(*,'(/2(a/))') &
        'Error in define_essential:', &
        ' input_probdef has not been created '
      stop
    end if

    if ( present(physq) ) then
      if ( physq < 1 .or. physq > input_probdef%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq in heading of define_essential has wrong value:', &
          'physq < 1 or physq > number of physical quantities = ', &
          input_probdef%nphysq
        stop
      end if
    end if

    if ( present(layer) ) then
      if ( input_probdef%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in define_essential: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer < 1 .or. layer > input_probdef%numlayers ) then
        write(*,'(/a/a,i0/)') &
          'Error: layer in heading of define_essential has wrong value:', &
          'layer < 1 or layer > number of layers = ', &
          input_probdef%numlayers
        stop
      end if
    end if

    if ( present(excludepoints) ) then
      if ( any(excludepoints <=0) .or. any(excludepoints > mesh%npoints) ) then
        write(*,'(2(/a),i0/)') &
          'Error in define_essential: ',&
          ' excludepoints <=0 or excludepoints > number of points = ', &
           mesh%npoints
        stop
      end if
    end if

    if ( present(excludecurves) ) then
      if ( any(excludecurves <=0) .or. any(excludecurves > mesh%ncurves) ) then
        write(*,'(2(/a),i0/)') &
          'Error in define_essential: ',&
          ' excludecurves <=0 or excludecurves > number of curves = ', &
           mesh%ncurves
        stop
      end if
    end if

    if ( present(excludesurfaces) ) then
      if ( any(excludesurfaces <=0) .or. &
           any(excludesurfaces > mesh%nsurfaces) ) then
        write(*,'(2(/a),i0/)') &
          'Error in define_essential: ',&
          ' excludesurfaces <=0 or excludesurfaces > number of surfaces = ', &
           mesh%nsurfaces
        stop
      end if
    end if

    if ( present(degfd) .and. present(degsfd) ) then
      write(*,'(2(/a),i0/)') &
        'Error in define_essential: ',&
        ' degfd and degsfd are both present'
      stop
    end if

    if ( present(degsfd) ) then

!     set ldegfd

      if ( any(degsfd<=0) ) then
        write(*,'(2(/a),i0/)') &
          'Error in define_essential: degsfd <= 0 '
        stop
      end if

      if ( size(degsfd) == 0 ) then
        allocate(ldegfd(0))
      else if ( maxval(degsfd) == 0 ) then
        ldegfd = [0]
      else
        allocate(ldegfd(maxval(degsfd)))
        ldegfd = 0
        do i = 1, size(degsfd)
          if ( degsfd(i) > 0 ) ldegfd(degsfd(i)) = 1
        end do
      end if

    end if

!   determine type of call

    if ( present(point) .or. present(points) ) then

!     fill input_probdef for a point or points

      if ( present(point) ) then

        if ( point <=0 .or. point > mesh%npoints ) then
          write(*,'(/a/a,i0/)') &
            'Error: point in heading of define_essential has wrong value:', &
            'point <=0 or point > number of points = ', mesh%npoints
          stop
        end if

      end if

      if ( present(points) ) then

        if ( any( points <=0 .or. points > mesh%npoints ) ) then
          write(*,'(/a/a,i0/)') &
            'Error: points in heading of define_essential has wrong values:', &
            'points <=0 or points > number of points = ', mesh%npoints
          stop
        end if

      end if

      input_probdef%numesspoints = input_probdef%numesspoints + 1
      pnt = input_probdef%numesspoints

      if ( pnt > size(input_probdef%esspoints) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_essential:', &
          ' essential point > maximum = ', size(input_probdef%esspoints)
        stop
      end if

      if ( present(num) ) num = pnt

      if ( present(physq) ) then
        input_probdef%esspoints(pnt)%physq = physq
      else
        input_probdef%esspoints(pnt)%physq = 0
      end if

      if ( present(layer) ) then
        input_probdef%esspoints(pnt)%layer = layer
      else
        input_probdef%esspoints(pnt)%layer = 0
      end if

      if ( present(point) ) then
        input_probdef%esspoints(pnt)%point = point
      else
        input_probdef%esspoints(pnt)%point = 0
      end if

      if ( present(points) ) then
        input_probdef%esspoints(pnt)%points = points
      else
        allocate (input_probdef%esspoints(pnt)%points(0) )
      end if

      if ( present(degfd) ) then
        input_probdef%esspoints(pnt)%degfd = degfd
      else if ( present(degsfd) ) then
        input_probdef%esspoints(pnt)%degfd = ldegfd
      else
        allocate (input_probdef%esspoints(pnt)%degfd(0) )
      end if

    else if ( present(curve1) .or. present(curves) ) then

!     fill input_probdef for curves

      if ( present(curve1) ) then

        if ( curve1 <=0 .or. curve1 > mesh%ncurves ) then
          write(*,'(/a/a,i0/)') &
            'Error: curve1 in heading of define_essential has wrong value:', &
            'curve1 <=0 or curve1 > number of curves = ', mesh%ncurves
          stop
        end if

      end if

      if ( present(curves) ) then

        if ( any ( curves <=0 .or. curves > mesh%ncurves ) ) then
          write(*,'(/a/a,i0/)') &
            'Error: curves in heading of define_essential has wrong values:', &
            'curves <=0 or curves > number of curves = ', mesh%ncurves
          stop
        end if

      end if

      input_probdef%numesscurves = input_probdef%numesscurves + 1
      cg = input_probdef%numesscurves

      if ( cg > size(input_probdef%esscurves) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_essential:', &
          ' essential curve > maximum = ', size(input_probdef%esscurves)
        stop
      end if

      if ( present(num) ) num = cg

      if ( present(physq) ) then
        input_probdef%esscurves(cg)%physq = physq
      else
        input_probdef%esscurves(cg)%physq = 0
      end if

      if ( present(layer) ) then
        input_probdef%esscurves(cg)%layer = layer
      else
        input_probdef%esscurves(cg)%layer = 0
      end if

      if ( present(curves) ) then
        input_probdef%esscurves(cg)%curves = curves
      else
        allocate (input_probdef%esscurves(cg)%curves(0) )
      end if

      if ( present(degfd) ) then
        input_probdef%esscurves(cg)%degfd = degfd
      else if ( present(degsfd) ) then
        input_probdef%esscurves(cg)%degfd = ldegfd
      else
        allocate (input_probdef%esscurves(cg)%degfd(0) )
      end if

      if ( present(curve1) ) then
        input_probdef%esscurves(cg)%first = curve1
        if ( present(curve2) ) then
          if ( curve2 <=0 .or. curve2 > mesh%ncurves ) then
            write(*,'(/a/a,i0/)') &
              'Error: curve2 in heading of define_essential has wrong value:', &
              'curve2 <=0 or curve2 > number of curves = ', mesh%ncurves
            stop
          end if
          input_probdef%esscurves(cg)%last = curve2
        else
          input_probdef%esscurves(cg)%last = curve1
        end if
      end if

      if ( present(step) ) then
        input_probdef%esscurves(cg)%step = step
      end if

      if ( present(exclude) ) then
        input_probdef%esscurves(cg)%exclude = exclude
      end if

      if ( present(excludepoints) ) then
        input_probdef%esscurves(cg)%excludepoints = excludepoints
      else
        allocate (input_probdef%esscurves(cg)%excludepoints(0) )
      end if

      if ( present(excludecurves) ) then
        input_probdef%esscurves(cg)%excludecurves = excludecurves
      else
        allocate (input_probdef%esscurves(cg)%excludecurves(0) )
      end if

      if ( present(excludesurfaces) ) then
        input_probdef%esscurves(cg)%excludesurfaces = excludesurfaces
      else
        allocate (input_probdef%esscurves(cg)%excludesurfaces(0) )
      end if

    else if ( present(surface1) .or. present(surfaces) ) then

!     fill input_probdef for surfaces

      if ( present(surface1) ) then

        if ( surface1 <=0 .or. surface1 > mesh%nsurfaces ) then
          write(*,'(/a/a,i0/)') &
            'Error: surface1 in heading of define_essential has wrong value:', &
            'surface1 <=0 or surface1 > number of surfaces = ', mesh%nsurfaces
          stop
        end if

      end if

      if ( present(surfaces) ) then

        if ( any( surfaces <=0 .or. surfaces > mesh%nsurfaces ) ) then
          write(*,'(/a/a,i0/)') &
            'Error: surfaces in heading of define_essential has wrong values:',&
            'surfaces <=0 or surfaces > number of surfaces = ', mesh%nsurfaces
          stop
        end if

      end if

      input_probdef%numesssurfaces = input_probdef%numesssurfaces + 1
      cg = input_probdef%numesssurfaces

      if ( cg > size(input_probdef%esssurfaces) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_essential:', &
          ' essential surface > maximum = ', size(input_probdef%esssurfaces)
        stop
      end if

      if ( present(num) ) num = cg

      if ( present(physq) ) then
        input_probdef%esssurfaces(cg)%physq = physq
      else
        input_probdef%esssurfaces(cg)%physq = 0
      end if

      if ( present(layer) ) then
        input_probdef%esssurfaces(cg)%layer = layer
      else
        input_probdef%esssurfaces(cg)%layer = 0
      end if

      if ( present(degfd) ) then
        input_probdef%esssurfaces(cg)%degfd = degfd
      else if ( present(degsfd) ) then
        input_probdef%esssurfaces(cg)%degfd = ldegfd
      else
        allocate (input_probdef%esssurfaces(cg)%degfd(0) )
      end if

      if ( present(surfaces) ) then
        input_probdef%esssurfaces(cg)%surfaces = surfaces
      else
        allocate (input_probdef%esssurfaces(cg)%surfaces(0) )
      end if

      if ( present(surface1) ) then
        input_probdef%esssurfaces(cg)%first = surface1
        if ( present(surface2) ) then
          if ( surface2 <=0 .or. surface2 > mesh%nsurfaces ) then
            write(*,'(/a/a,i0/)') &
             'Error: surface2 in heading of define_essential has wrong value:',&
             'surface2 <=0 or surface2 > number of surfaces = ', mesh%nsurfaces
            stop
          end if
          input_probdef%esssurfaces(cg)%last = surface2
        else
          input_probdef%esssurfaces(cg)%last = surface1
        end if
      end if

      if ( present(excludepoints) ) then
        input_probdef%esssurfaces(cg)%excludepoints = excludepoints
      else
        allocate (input_probdef%esssurfaces(cg)%excludepoints(0) )
      end if

      if ( present(excludecurves) ) then
        input_probdef%esssurfaces(cg)%excludecurves = excludecurves
      else
        allocate (input_probdef%esssurfaces(cg)%excludecurves(0) )
      end if

      if ( present(excludesurfaces) ) then
        input_probdef%esssurfaces(cg)%excludesurfaces = excludesurfaces
      else
        allocate (input_probdef%esssurfaces(cg)%excludesurfaces(0) )
      end if

    else if ( present(element) ) then

!     fill input_probdef for an element

      if ( present(group) ) then
        if ( group <=0 .or. group > mesh%nelgrp ) then
          write(*,'(/a/a,i0/)') &
            'Error: group in heading of define_essential has wrong value:', &
            'group <=0 or group > number of groups = ', mesh%nelgrp
          stop
        end if
        elgrp = group
      else
        elgrp = 1
      end if

      if ( element <=0 .or. element > mesh%grpnumel(elgrp) ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error: element in heading of define_essential has wrong value:', &
          'element <=0 or element > number of elements ', &
          mesh%grpnumel(elgrp), ' in element group ', elgrp
        stop
      end if

      if ( present(elnode) ) then
        if ( elnode <=0 .or. elnode > mesh%elnumnod(elgrp) ) then
          write(*,'(/a/2(a,i0)/)') &
            'Error: elnode in heading of define_essential has wrong value:', &
            'elnode <=0 or elnode > number of nodes in element = ', &
            mesh%elnumnod(elgrp), ' in element group ', elgrp
          stop
        end if
        elnod = elnode
      else
        elnod = 1
      end if

      input_probdef%numesselements = input_probdef%numesselements + 1
      elm = input_probdef%numesselements

      if ( elm > size(input_probdef%esselements) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_essential:', &
          ' essential elements > maximum = ', size(input_probdef%esselements)
        stop
      end if

      if ( present(num) ) num = elm

      if ( present(physq) ) then
        input_probdef%esselements(elm)%physq = physq
      else
        input_probdef%esselements(elm)%physq = 0
      end if

      if ( present(layer) ) then
        input_probdef%esselements(elm)%layer = layer
      else
        input_probdef%esselements(elm)%layer = 0
      end if

      input_probdef%esselements(elm)%elgrp = elgrp
      input_probdef%esselements(elm)%elem = element
      input_probdef%esselements(elm)%node = elnod

      if ( present(degfd) ) then
        input_probdef%esselements(elm)%degfd = degfd
      else if ( present(degsfd) ) then
        input_probdef%esselements(elm)%degfd = ldegfd
      else
        allocate (input_probdef%esselements(elm)%degfd(0) )
      end if

    else if ( present(elgroup1) .or. present(elgroups) ) then

!     fill input_probdef for groups

      if ( present(elgroup1) ) then

        if ( elgroup1 <=0 .or. elgroup1 > mesh%nelgrp ) then
          write(*,'(/a/a,i0/)') &
            'Error: elgroup1 in heading of define_essential has wrong value:', &
            'elgroup1 <=0 or elgroup1 > number of groups = ', mesh%nelgrp
          stop
        end if

      end if

      if ( present(elgroups) ) then

        if ( any ( elgroups <=0 .or. elgroups > mesh%nelgrp ) ) then
          write(*,'(/a/a,i0/)') &
            'Error: elgroups in heading of define_essential has wrong values:',&
            'elgroups <=0 or elgroups > number of groups = ', mesh%nelgrp
          stop
        end if

      end if

      input_probdef%numessgroups = input_probdef%numessgroups + 1
      cg = input_probdef%numessgroups

      if ( cg > size(input_probdef%essgroups) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_essential:', &
          ' essential group > maximum = ', size(input_probdef%essgroups)
        stop
      end if

      if ( present(num) ) num = cg

      if ( present(physq) ) then
        input_probdef%essgroups(cg)%physq = physq
      else
        input_probdef%essgroups(cg)%physq = 0
      end if

      if ( present(layer) ) then
        input_probdef%essgroups(cg)%layer = layer
      else
        input_probdef%essgroups(cg)%layer = 0
      end if

      if ( present(elgroups) ) then
        input_probdef%essgroups(cg)%elgroups = elgroups
      else
        allocate (input_probdef%essgroups(cg)%elgroups(0) )
      end if

      if ( present(degfd) ) then
        input_probdef%essgroups(cg)%degfd = degfd
      else if ( present(degsfd) ) then
        input_probdef%essgroups(cg)%degfd = ldegfd
      else
        allocate (input_probdef%essgroups(cg)%degfd(0) )
      end if

      if ( present(elgroup1) ) then
        input_probdef%essgroups(cg)%first = elgroup1
        if ( present(elgroup2) ) then
          if ( elgroup2 <=0 .or. elgroup2 > mesh%nelgrp ) then
            write(*,'(/a/a,i0/)') &
             'Error: elgroup2 in heading of define_essential has wrong value:',&
             'elgroup2 <=0 or elgroup2 > number of groups = ', mesh%nelgrp
            stop
          end if
          input_probdef%essgroups(cg)%last = elgroup2
        else
          input_probdef%essgroups(cg)%last = elgroup1
        end if
      end if

      if ( present(excludepoints) ) then
        input_probdef%essgroups(cg)%excludepoints = excludepoints
      else
        allocate (input_probdef%essgroups(cg)%excludepoints(0) )
      end if

      if ( present(excludecurves) ) then
        input_probdef%essgroups(cg)%excludecurves = excludecurves
      else
        allocate (input_probdef%essgroups(cg)%excludecurves(0) )
      end if

      if ( present(excludesurfaces) ) then
        input_probdef%essgroups(cg)%excludesurfaces = excludesurfaces
      else
        allocate (input_probdef%essgroups(cg)%excludesurfaces(0) )
      end if

    else if ( present(nodeset1) .or. present(nodesets) ) then

!     fill input_probdef for nodesets

      if ( present(nodeset1) ) then

        if ( nodeset1 <=0 .or. nodeset1 > mesh%nnodesets ) then
          write(*,'(/a/a,i0/)') &
            'Error: nodeset1 in heading of define_essential has wrong value:', &
            'nodeset1 <=0 or nodeset1 > number of nodesets = ', mesh%nnodesets
          stop
        end if

      end if

      if ( present(nodesets) ) then

        if ( any ( nodesets <=0 .or. nodesets > mesh%nnodesets ) ) then
          write(*,'(/a/a,i0/)') &
            'Error: nodesets in heading of define_essential has wrong values:',&
            'nodesets <=0 or nodesets > number of nodesets = ', mesh%nnodesets
          stop
        end if

      end if

      input_probdef%numessnodesets = input_probdef%numessnodesets + 1
      cg = input_probdef%numessnodesets

      if ( cg > size(input_probdef%essnodesets) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_essential:', &
          ' essential nodeset > maximum = ', size(input_probdef%essnodesets)
        stop
      end if

      if ( present(num) ) num = cg

      if ( present(physq) ) then
        input_probdef%essnodesets(cg)%physq = physq
      else
        input_probdef%essnodesets(cg)%physq = 0
      end if

      if ( present(layer) ) then
        input_probdef%essnodesets(cg)%layer = layer
      else
        input_probdef%essnodesets(cg)%layer = 0
      end if

      if ( present(degfd) ) then
        input_probdef%essnodesets(cg)%degfd = degfd
      else if ( present(degsfd) ) then
        input_probdef%essnodesets(cg)%degfd = ldegfd
      else
        allocate (input_probdef%essnodesets(cg)%degfd(0) )
      end if

      if ( present(nodesets) ) then
        input_probdef%essnodesets(cg)%nodesets = nodesets
      else
        allocate (input_probdef%essnodesets(cg)%nodesets(0) )
      end if

      if ( present(nodeset1) ) then
        input_probdef%essnodesets(cg)%first = nodeset1
        if ( present(nodeset2) ) then
          if ( nodeset2 <=0 .or. nodeset2 > mesh%nnodesets ) then
            write(*,'(/a/a,i0/)') &
             'Error: nodeset2 in heading of define_essential has wrong value:',&
             'nodeset2 <=0 or nodeset2 > number of nodesets = ', mesh%nnodesets
            stop
          end if
          input_probdef%essnodesets(cg)%last = nodeset2
        else
          input_probdef%essnodesets(cg)%last = nodeset1
        end if
      end if

      if ( present(excludepoints) ) then
        input_probdef%essnodesets(cg)%excludepoints = excludepoints
      else
        allocate (input_probdef%essnodesets(cg)%excludepoints(0) )
      end if

      if ( present(excludecurves) ) then
        input_probdef%essnodesets(cg)%excludecurves = excludecurves
      else
        allocate (input_probdef%essnodesets(cg)%excludecurves(0) )
      end if

      if ( present(excludesurfaces) ) then
        input_probdef%essnodesets(cg)%excludesurfaces = excludesurfaces
      else
        allocate (input_probdef%essnodesets(cg)%excludesurfaces(0) )
      end if

    else

      write(*,'(9(/a)/)') &
        'Error: could not determine type of essential boundary condition', &
        'from the heading of define_essential. None of the keywords:',     &
        '  point or points ',                                              &
        '  curve1 or curves ',                                             &
        '  surface1 or surfaces ',                                         &
        '  element ',                                                      &
        '  elgroup1 or elgroups ',                                         &
        '  nodeset1 or nodesets ',                                         &
        'is present.'
      stop

    end if

  end subroutine define_essential


! Change the essential boundary conditions
! NOTE: no defaults are assumed! If a keyword is not present the value in the
! structure will not be changed with respect to value given in define_essential.

  subroutine change_essential ( mesh, input_probdef, bctype, nr, physq, layer, &
    degfd, degsfd, points, point, curves, curve1, curve2, step, exclude, &
    group, element, elnode, surfaces, surface1, surface2, elgroups, elgroup1, &
    elgroup2, nodesets, nodeset1, nodeset2, excludepoints, excludecurves, &
    excludesurfaces )

    type(mesh_t), intent(in) :: mesh

!   structure to store all definitions for the problem definition
    type(input_probdef_t), intent(inout) :: input_probdef

!   type of essential boundary condition
!     bctype='point' or 'points', 'curves', 'surfaces', 'element',
!            'groups' or 'elgroups'
!   referring to the bc's stored in the arrays esspoints(:), esscurves(:),
!   esssurfaces(:), esselements(:), essgroups(:), respectively.
    character(len=*), intent(in) :: bctype

!   the essential boundary condition sequence number, i.e. the index in the
!   arrays esspoints(:), esscurves(:), esssurfaces(:), esselements(:),
!   essgroups(:). Note, that each call of define_essential assigns a sequence
!   number to the boundary condition that is defined, starting with 1 for the
!   first call. Note also that the boundary conditions involving points, curves,
!   surfaces, elements and groups each has their own sequence. Hence, the first
!   call of define_essential for curves and the first call of define_essential
!   for points will both have a sequence number of 1.
    integer, intent(in) :: nr

!   if physq is present physq denotes the physical quantity that must be set
!   otherwise all nodal degrees of freedom are considered
    integer, intent(in), optional :: physq

!   if layer is present, nodal degrees of freedom to be set are restricted to
!   the given layer otherwise all nodal degrees of freedom are considered.
    integer, intent(in), optional :: layer

!   if degfd is present the array degfd indicates which degrees of freedom are
!   essential:
!       degfd(i) > 0 => degfd i is set
!   for example degfd=[1,0,1] sets the first and third degree of freedom
!   to be essential. If a degree of freedom is set that does not exist, it is
!   simply ignored. Only the first NUMBEROFBITS degrees of freedom in a node
!   can be set.
!   NOTE: if size(degfd)=0 all degrees are set
!         if degfd=[0] no degrees are set, inactivating the essential bc.
!         if physq is present degfd refers to the physical quantity only
!         if layer is present degfd refers to the layer only
!   NOTE: degfd and degsfd are not allowed to be both present.
    integer, intent(in), dimension(:), optional :: degfd

!   if degsfd is present the array degsfd indicates which degrees of freedom are
!   essential. For example degsfd=[1,3] sets the first and third degree of
!   freedom to be essential. If a degree of freedom is set that does not exist,
!   it is simply ignored. Only the first NUMBEROFBITS degrees of freedom in a
!   node can be set.
!   NOTE: if size(degsfd)=0 all degrees are set, as if degsfd is not present
!         if degsfd=[0] no degrees are set, inactivating the essential bc.
!         if physq is present degsfd refers to the physical quantity only
!         if layer is present degsfd refers to the layer only
!   NOTE: degfd and degsfd are not allowed to be both present.
    integer, intent(in), dimension(:), optional :: degsfd

!   if point is present essential boundary conditions are set in a point
    integer, intent(in), optional :: point

!   if points is present essential boundary conditions are set in the points
!   specified. Note: points can be combines with point in one call.
    integer, intent(in), dimension(:), optional :: points

!   if curve1 and/or curve2 is present essential boundary conditions are set
!   along curves curve1 to curve2.
!
!   step indicates which nodes on the curves are used:
!       step=0 all nodes
!       step>0 nodes 1, 1+step, 1+2*step, ...
!       step<0 except nodes 1, 1-step, 1-2*step, ...
!   note that this is applied for all curves _separately_ and curve1 to curve2
!   is not treated as a single curve.
!
!   exclude indicates which nodes on the curves are excluded:
!       exclude=0 no excludes, all nodes are used (default)
!       exclude=1 exclude first point of first curve
!       exclude=2 exclude last point of last curve
!       exclude=3 exclude first and last point
!
!   NOTE: step, exclude only work for curves where the local node numbering is
!   in a natural sequence along the curve. This might not be the case for
!   curves constructed from several other curves or curves read from external
!   mesh generators.
    integer, intent(in), optional :: curve1, curve2, step, exclude

!   if curves is present essential boundary conditions are set in the curves
!   specified. Note: curves can be combines with curve1 in one call.
    integer, intent(in), dimension(:), optional :: curves

!   if element is present essential boundary conditions are set in a element
    integer, intent(in), optional :: element

!   if group is present essential boundary conditions for an element in the
!   specified group is set.
    integer, intent(in), optional :: group

!   if elnode is present essential boundary conditions for an element in the
!   local node elnode is set. Default=1.
    integer, intent(in), optional :: elnode

!   if surface1 is present essential boundary conditions are set along surface
!   surface1 or surface1 to surface2 if surface2 is also present.
    integer, intent(in), optional :: surface1, surface2

!   if surfaces is present essential boundary conditions are set in the surfaces
!   specified. Note: surfaces can be combines with surface1 in one call.
    integer, intent(in), dimension(:), optional :: surfaces

!   if elgroup1 is present essential boundary conditions are set in the element
!   group elgroup1 or elgroup1 to elgroup2 if group2 is also present.
    integer, intent(in), optional :: elgroup1, elgroup2

!   if elgroups is present essential boundary conditions are set in the element
!   groups specified. Note: elgroups can be combines with elgroup1 in one call.
    integer, intent(in), dimension(:), optional :: elgroups

!   if nodeset1 is present essential boundary conditions are set along nodeset
!   nodeset1 or nodeset1 to nodeset2 if nodeset2 is also present.
    integer, intent(in), optional :: nodeset1, nodeset2

!   if nodesets is present essential boundary conditions are set in the nodesets
!   specified. Note: nodesets can be combines with nodeset1 in one call.
    integer, intent(in), dimension(:), optional :: nodesets

!   exclude points:
!   for example excludepoints=(/1,3/) excludes nodes on points P1 and P3
!   to be essential. Can be used together with the parameters curve1, curve2,
!   surface1, surface2, elgroup1, elgroup2, nodeset1, nodeset2. Note that only
!   the nodes that are actually on the curves, surfaces, groups or nodesets
!   are excluded.
    integer, intent(in), dimension(:), optional :: excludepoints

!   exclude curves:
!   for example excludecurves=(/1,3/) excludes nodes on curves C1 and C3
!   to be essential. Can be used together with the parameters curve1, curve2,
!   surface1, surface2, elgroup1, elgroup2, nodeset1, nodeset2. Note that only
!   the nodes that are actually on the curves, surfaces, groups or nodesets are
!   excluded.
    integer, intent(in), dimension(:), optional :: excludecurves

!   exclude surfaces:
!   for example excludesurfaces=(/1,3/) excludes nodes on surfaces S1 and S3
!   to be essential. Can be used together with the parameters curve1, curve2,
!   surface1, surface2, elgroup1, elgroup2, nodeset1, nodeset2. Note that only
!   the nodes that are actually on the surfaces, groups or nodesets are
!   excluded.
    integer, intent(in), dimension(:), optional :: excludesurfaces


    integer :: elgrp, i
    integer, allocatable, dimension(:) :: ldegfd

!   some testing

    call check ( mesh, 'change_essential' )

    if ( .not. input_probdef%created ) then
      write(*,'(/2(a/))') &
        'Error in change_essential:', &
        ' input_probdef has not been created '
      stop
    end if

    if ( present(physq) ) then
      if ( physq < 1 .or. physq > input_probdef%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq in heading of change_essential has wrong value:', &
          'physq < 1 or physq > number of physical quantities = ', &
          input_probdef%nphysq
        stop
      end if
    end if

    if ( present(layer) ) then
      if ( input_probdef%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in change_essential: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer < 1 .or. layer > input_probdef%numlayers ) then
        write(*,'(/a/a,i0/)') &
          'Error: layer in heading of change_essential has wrong value:', &
          'layer < 1 or layer > number of layers = ', &
          input_probdef%numlayers
        stop
      end if
    end if

    if ( present(excludepoints) ) then
      if ( any(excludepoints <=0) .or. any(excludepoints > mesh%npoints) ) then
        write(*,'(2(/a),i0/)') &
          'Error in change_essential: ',&
          ' excludepoints <=0 or excludepoints > number of points = ', &
           mesh%npoints
        stop
      end if
    end if

    if ( present(excludecurves) ) then
      if ( any(excludecurves <=0) .or. any(excludecurves > mesh%ncurves) ) then
        write(*,'(2(/a),i0/)') &
          'Error in change_essential: ',&
          ' excludecurves <=0 or excludecurves > number of curves = ', &
           mesh%ncurves
        stop
      end if
    end if

    if ( present(excludesurfaces) ) then
      if ( any(excludesurfaces <=0) .or. &
           any(excludesurfaces > mesh%nsurfaces) ) then
        write(*,'(2(/a),i0/)') &
          'Error in change_essential: ',&
          ' excludesurfaces <=0 or excludesurfaces > number of surfaces = ', &
           mesh%nsurfaces
        stop
      end if
    end if

    if ( present(degfd) .and. present(degsfd) ) then
      write(*,'(2(/a),i0/)') &
        'Error in change_essential: ',&
        ' degfd and degsfd are both present'
      stop
    end if

    if ( present(degsfd) ) then

!     set ldegfd

      if ( any(degsfd<=0) ) then
        write(*,'(2(/a),i0/)') &
          'Error in change_essential: degsfd <= 0 '
        stop
      end if

      if ( size(degsfd) == 0 ) then
        allocate(ldegfd(0))
      else if ( maxval(degsfd) == 0 ) then
        ldegfd = [0]
      else
        allocate(ldegfd(maxval(degsfd)))
        ldegfd = 0
        do i = 1, size(degsfd)
          if ( degsfd(i) > 0 ) ldegfd(degsfd(i)) = 1
        end do
      end if

    end if

!   determine type of call

    if ( bctype == 'point' .or. bctype == 'points'  ) then

!     change input_probdef for a point

      if ( nr < 1 .or. nr > input_probdef%numesspoints ) then
        write(*,'(/a/a,i0/)') &
          'Error: invalid nr in change_essential:', &
          ' nr  < 1 or > numesspoints = ', input_probdef%numesspoints
        stop
      end if

      if ( present(physq) ) then
        input_probdef%esspoints(nr)%physq = physq
      end if

      if ( present(layer) ) then
        input_probdef%esspoints(nr)%layer = layer
      end if

      if ( present(point) ) then
        if ( point <=0 .or. point > mesh%npoints ) then
          write(*,'(/a/a,i0/)') &
            'Error: point in heading of change_essential has wrong value:', &
            'point <=0 or point > number of points = ', mesh%npoints
          stop
        end if
        input_probdef%esspoints(nr)%point = point
      else
        input_probdef%esspoints(nr)%point = 0
      end if

      deallocate (input_probdef%esspoints(nr)%points)

      if ( present(points) ) then
        if ( any ( points <=0 .or. points > mesh%npoints ) ) then
          write(*,'(/a/a,i0/)') &
            'Error: points in heading of change_essential has wrong values:', &
            'points <=0 or points > number of points = ', mesh%npoints
          stop
        end if
        input_probdef%esspoints(nr)%points = points
      else
        allocate (input_probdef%esspoints(nr)%points(0) )
      end if

      if ( present(degfd) ) then
        input_probdef%esspoints(nr)%degfd = degfd
      else if ( present(degsfd) ) then
        input_probdef%esspoints(nr)%degfd = ldegfd
      end if

    else if ( bctype == 'curves' ) then

!     change input_probdef for curves

      if ( nr < 1 .or. nr > input_probdef%numesscurves ) then
        write(*,'(/a/a,i0/)') &
          'Error: invalid nr in change_essential:', &
          ' nr < 1 or > numesscurves = ', input_probdef%numesscurves
        stop
      end if

      if ( present(physq) ) then
        input_probdef%esscurves(nr)%physq = physq
      end if

      if ( present(layer) ) then
        input_probdef%esscurves(nr)%layer = layer
      end if

      if ( present(degfd) ) then
        input_probdef%esscurves(nr)%degfd = degfd
      else if ( present(degsfd) ) then
        input_probdef%esscurves(nr)%degfd = ldegfd
      end if

      if ( present(curves) ) then
        if ( any ( curves <=0 .or. curves > mesh%ncurves ) ) then
          write(*,'(/a/a,i0/)') &
            'Error: curves in heading of change_essential has wrong values:', &
            'curves <=0 or curves > number of curves = ', mesh%ncurves
          stop
        end if
        input_probdef%esscurves(nr)%curves = curves
      end if

      if ( present(curve1) ) then
        if ( curve1 <=0 .or. curve1 > mesh%ncurves ) then
          write(*,'(/a/a,i0/)') &
            'Error: curve1 in heading of change_essential has wrong value:', &
            'curve1 <=0 or curve1 > number of curves = ', mesh%ncurves
          stop
        end if
        input_probdef%esscurves(nr)%first = curve1
      end if

      if ( present(curve2) ) then
        if ( curve2 <=0 .or. curve2 > mesh%ncurves ) then
          write(*,'(/a/a,i0/)') &
            'Error: curve2 in heading of change_essential has wrong value:', &
            'curve2 <=0 or curve2 > number of curves = ', mesh%ncurves
          stop
        end if
        input_probdef%esscurves(nr)%last = curve2
      end if

      if ( present(step) ) then
        input_probdef%esscurves(nr)%step = step
      end if

      if ( present(exclude) ) then
        input_probdef%esscurves(nr)%exclude = exclude
      end if

      if ( present(excludepoints) ) then
        input_probdef%esscurves(nr)%excludepoints = excludepoints
      end if

      if ( present(excludecurves) ) then
        input_probdef%esscurves(nr)%excludecurves = excludecurves
      end if

      if ( present(excludesurfaces) ) then
        input_probdef%esscurves(nr)%excludesurfaces = excludesurfaces
      end if

    else if ( bctype == 'surfaces' ) then

!     change input_probdef for surfaces

      if ( nr < 1 .or. nr > input_probdef%numesssurfaces ) then
        write(*,'(/a/a,i0/)') &
          'Error: invalid nr in change_essential:', &
          ' nr < 1 or > numesssurfaces = ', input_probdef%numesssurfaces
        stop
      end if

      if ( present(physq) ) then
        input_probdef%esssurfaces(nr)%physq = physq
      end if

      if ( present(layer) ) then
        input_probdef%esssurfaces(nr)%layer = layer
      end if

      if ( present(degfd) ) then
        input_probdef%esssurfaces(nr)%degfd = degfd
      else if ( present(degsfd) ) then
        input_probdef%esssurfaces(nr)%degfd = ldegfd
      end if

      if ( present(surfaces) ) then
        if ( any ( surfaces <=0 .or. surfaces > mesh%nsurfaces ) ) then
          write(*,'(/a/a,i0/)') &
            'Error: surfaces in heading of change_essential has wrong values:',&
            'surfaces <=0 or surfaces > number of surfaces = ', mesh%nsurfaces
          stop
        end if
        input_probdef%esssurfaces(nr)%surfaces = surfaces
      end if

      if ( present(surface1) ) then
        if ( surface1 <=0 .or. surface1 > mesh%nsurfaces ) then
          write(*,'(/a/a,i0/)') &
            'Error: surface1 in heading of change_essential has wrong value:', &
            'surface1 <=0 or surface1 > number of surfaces = ', mesh%nsurfaces
          stop
        end if
        input_probdef%esssurfaces(nr)%first = surface1
      end if

      if ( present(surface2) ) then
        if ( surface2 <=0 .or. surface2 > mesh%nsurfaces ) then
          write(*,'(/a/a,i0/)') &
            'Error: surface2 in heading of change_essential has wrong value:', &
            'surface2 <=0 or surface2 > number of surfaces = ', mesh%nsurfaces
          stop
        end if
        input_probdef%esssurfaces(nr)%last = surface2
      end if

      if ( present(excludepoints) ) then
        input_probdef%esssurfaces(nr)%excludepoints = excludepoints
      end if

      if ( present(excludecurves) ) then
        input_probdef%esssurfaces(nr)%excludecurves = excludecurves
      end if

      if ( present(excludesurfaces) ) then
        input_probdef%esssurfaces(nr)%excludesurfaces = excludesurfaces
      end if

    else if ( bctype == 'element' ) then

!     change input_probdef for an element

      if ( nr < 1 .or. nr > input_probdef%numesselements ) then
        write(*,'(/a/a,i0/)') &
          'Error: invalid nr in change_essential:', &
          ' nr < 1 or > numesselements = ', input_probdef%numesselements
        stop
      end if

      if ( present(group) ) then
        if ( group <=0 .or. group > mesh%nelgrp ) then
          write(*,'(/a/a,i0/)') &
            'Error: group in heading of change_essential has wrong value:', &
            'group <=0 or group > number of groups = ', mesh%nelgrp
          stop
        end if
        input_probdef%esselements(nr)%elgrp = elgrp
        elgrp = group
      else
        elgrp = input_probdef%esselements(nr)%elgrp
      end if

      if ( present(element) ) then
        if ( element <= 0 .or. element > mesh%grpnumel(elgrp) ) then
          write(*,'(/a/2(a,i0)/)') &
            'Error: element in heading of change_essential has wrong value:', &
            'element <=0 or element > number of elements ', &
            mesh%grpnumel(elgrp), ' in element group ', elgrp
          stop
        end if
        input_probdef%esselements(nr)%elem = element
      end if

      if ( present(elnode) ) then
        if ( elnode <=0 .or. elnode > mesh%elnumnod(elgrp) ) then
          write(*,'(/a/2(a,i0)/)') &
            'Error: elnode in heading of change_essential has wrong value:', &
            'elnode <=0 or elnode > number of nodes in element = ', &
            mesh%elnumnod(elgrp), ' in element group ', elgrp
          stop
        end if
        input_probdef%esselements(nr)%node = elnode
      end if

      if ( present(physq) ) then
        input_probdef%esselements(nr)%physq = physq
      end if

      if ( present(layer) ) then
        input_probdef%esselements(nr)%layer = layer
      end if

      if ( present(degfd) ) then
        input_probdef%esselements(nr)%degfd = degfd
      else if ( present(degsfd) ) then
        input_probdef%esselements(nr)%degfd = ldegfd
      end if

    else if ( bctype == 'groups' .or. bctype == 'elgroups' ) then

!     change input_probdef for groups

      if ( nr < 1 .or. nr > input_probdef%numessgroups ) then
        write(*,'(/a/a,i0/)') &
          'Error: invalid nr in change_essential:', &
          ' nr < 1 or > numessgroups = ', input_probdef%numessgroups
        stop
      end if

      if ( present(physq) ) then
        input_probdef%essgroups(nr)%physq = physq
      end if

      if ( present(layer) ) then
        input_probdef%essgroups(nr)%layer = layer
      end if

      if ( present(degfd) ) then
        input_probdef%essgroups(nr)%degfd = degfd
      else if ( present(degsfd) ) then
        input_probdef%essgroups(nr)%degfd = ldegfd
      end if

      if ( present(elgroups) ) then
        if ( any ( elgroups <=0 .or. elgroups > mesh%nelgrp ) ) then
          write(*,'(/a/a,i0/)') &
            'Error: elgroups in heading of change_essential has wrong values:',&
            'elgroups <=0 or elgroups > number of groups = ', mesh%nelgrp
          stop
        end if
        input_probdef%essgroups(nr)%elgroups = elgroups
      end if

      if ( present(elgroup1) ) then
        if ( elgroup1 <=0 .or. elgroup1 > mesh%nelgrp ) then
          write(*,'(/a/a,i0/)') &
            'Error: elgroup1 in heading of change_essential has wrong value:', &
            'elgroup1 <=0 or elgroup1 > number of groups = ', mesh%nelgrp
          stop
        end if
        input_probdef%essgroups(nr)%first = elgroup1
      end if

      if ( present(elgroup2) ) then
        if ( elgroup2 <=0 .or. elgroup2 > mesh%nelgrp ) then
          write(*,'(/a/a,i0/)') &
            'Error: elgroup2 in heading of change_essential has wrong value:', &
            'elgroup2 <=0 or elgroup2 > number of groups = ', mesh%nelgrp
          stop
        end if
        input_probdef%essgroups(nr)%last = elgroup2
      end if

      if ( present(excludepoints) ) then
        input_probdef%essgroups(nr)%excludepoints = excludepoints
      end if

      if ( present(excludecurves) ) then
        input_probdef%essgroups(nr)%excludecurves = excludecurves
      end if

      if ( present(excludesurfaces) ) then
        input_probdef%essgroups(nr)%excludesurfaces = excludesurfaces
      end if

    else if ( bctype == 'nodesets' ) then

!     change input_probdef for nodesets

      if ( nr < 1 .or. nr > input_probdef%numessnodesets ) then
        write(*,'(/a/a,i0/)') &
          'Error: invalid nr in change_essential:', &
          ' nr < 1 or > numessnodesetss = ', input_probdef%numessnodesets
        stop
      end if

      if ( present(physq) ) then
        input_probdef%essnodesets(nr)%physq = physq
      end if

      if ( present(layer) ) then
        input_probdef%essnodesets(nr)%layer = layer
      end if

      if ( present(degfd) ) then
        input_probdef%essnodesets(nr)%degfd = degfd
      else if ( present(degsfd) ) then
        input_probdef%essnodesets(nr)%degfd = ldegfd
      end if

      if ( present(nodesets) ) then
        if ( any ( nodesets <=0 .or. nodesets > mesh%nnodesets ) ) then
          write(*,'(/a/a,i0/)') &
            'Error: nodesets in heading of change_essential has wrong values:',&
            'nodesets <=0 or nodesets > number of nodesets = ', mesh%nnodesets
          stop
        end if
        input_probdef%essnodesets(nr)%nodesets = nodesets
      end if

      if ( present(nodeset1) ) then
        if ( nodeset1 <=0 .or. nodeset1 > mesh%nnodesets ) then
          write(*,'(/a/a,i0/)') &
            'Error: nodeset1 in heading of change_essential has wrong value:', &
            'nodeset1 <=0 or nodeset1 > number of nodesets = ', mesh%nnodesets
          stop
        end if
        input_probdef%essnodesets(nr)%first = nodeset1
      end if

      if ( present(nodeset2) ) then
        if ( nodeset2 <=0 .or. nodeset2 > mesh%nnodesets ) then
          write(*,'(/a/a,i0/)') &
            'Error: nodeset2 in heading of change_essential has wrong value:', &
            'nodeset2 <=0 or nodeset2 > number of nodesets = ', mesh%nnodesets
          stop
        end if
        input_probdef%essnodesets(nr)%last = nodeset2
      end if

      if ( present(excludepoints) ) then
        input_probdef%essnodesets(nr)%excludepoints = excludepoints
      end if

      if ( present(excludecurves) ) then
        input_probdef%essnodesets(nr)%excludecurves = excludecurves
      end if

      if ( present(excludesurfaces) ) then
        input_probdef%essnodesets(nr)%excludesurfaces = excludesurfaces
      end if

    else

      write(*,'(9(/a)/)') &
        'Error: could not determine type of essential boundary condition', &
        ' from the heading of change_essential.', &
        ' Keyword bctype not equal to one of:',      &
        ' ''point or points''',                                              &
        ' ''curves''',                                                       &
        ' ''surfaces''',                                                     &
        ' ''element''',                                                      &
        ' ''groups or elgroups''',                                           &
        ' ''nodesets'''
      stop

    end if

  end subroutine change_essential


! Definition of constraints

  subroutine define_constraint ( mesh, input_probdef, physq, physq1, physq2, &
    layer, layer1, layer2, point1, point2, curve1, curve2, surface1, surface2, &
    volume, elementset1, elementset2, nodeset1, nodeset2, object, object2, &
    step, exclude, excludepoints, excludecurves, excludesurfaces, full2, &
    discretization, elementdof, nodedof, naddunknowns, nglobalc, &
    diagonal_block, diagonal_block_addunknowns, lmnumdegfd, num )

    type(mesh_t), intent(in) :: mesh

!   structure to store all definitions for the problem definition
    type(input_probdef_t), intent(inout) :: input_probdef

!   if physq is present physq denotes the physical quantity that must be
!   constrained otherwise all nodal degrees of freedom are considered
!   Note that this parameter does not specify the number of constraints, only
!   which unknowns are _involved_ in the constraints.
    integer, intent(in), optional :: physq

!   if physq1 and/or physq2 is present the physical quantity that must be
!   constrained can be different for geometry1 and geometry2, elementset1
!   and elementset2, the "two sides" of an object or object1 and object2.
!   Note that this parameter does not specify the number of constraints, only
!   which unknowns are _involved_ in the constraints.
    integer, intent(in), optional :: physq1, physq2

!   if layer is present layer denotes the layer that must be constrained
!   otherwise all nodal degrees of freedom are considered
!   Note that this parameter does not specify the number of constraints, only
!   which unknowns are _involved_ in the constraints. Only if lmnumdegfd>0 this
!   parameter also influences the number of constraints.
    integer, intent(in), optional :: layer

!   if layer1 and/or layer2 is present the layer that must be constrained
!   can be different for geometry1 and geometry2, elemenset1 and elementset2,
!   the "two sides" of an object or object1 and object2.
!   Note that this parameter does not specify the number of constraints, only
!   which unknowns are _involved_ in the constraints. Only if lmnumdegfd>0 this
!   parameter also influences the number of constraints.
    integer, intent(in), optional :: layer1, layer2

!   if point1 is present define a constraint in the point point1.
!   The Lagrange multipliers are defined in point1.
!   If point2 is also present a "link" of point point1 to point point2
!   is defined. Only discretization='collocation' is available for points.
    integer, intent(in), optional :: point1, point2

!   if curve1 is present define a constraint on the curve curve1.
!   The Lagrange multipliers are defined on curve1.
!   If curve2 is also present a "link" of curve curve1 to curve curve2
!   is defined.  Note that the direction of the curves are important: they are
!   connected starting from the beginning of the curves to the end of the
!   curves as they are defined in the mesh. Define new curves from existing
!   curves (using the routine add_to_mesh) when they don't match.
    integer, intent(in), optional :: curve1, curve2

!   if surface1 is present define a constraint on the surface surface1.
!   The Lagrange multipliers are defined on surface1.
!   if surface2 is also present a "link" of surface surface1 to surface
!   surface2 is defined. Note that the definition of the surfaces are important:
!   they are connected node for node, assuming that the topology is identical.
    integer, intent(in), optional :: surface1, surface2

!   if volume is present define a constraint on the volume volume.
    integer, intent(in), optional :: volume

!   if elementset1 is present define a constraint on the elementset elementset1.
!   The Lagrange multipliers are defined on elementset1.
!   if elementset2 is also present a "link" of elementset1 to elementset2
!   is defined. Note that the definition of the elementsets are important:
!   they are connected node for node, assuming that the topology is identical.
    integer, intent(in), optional :: elementset1, elementset2

!   if nodeset1 is present define a constraint on the nodeset nodeset1.
!   The Lagrange multipliers are defined on nodeset1.
!   if nodeset2 is also present a "link" of nodeset nodeset1 to nodeset
!   nodeset2 is defined. Note that the definition of the nodesets are important:
!   they are connected node for node. Only discretization='collocation' is
!   available for nodesets.
    integer, intent(in), optional :: nodeset1, nodeset2

!   if object is present define a constraint on the object object.
!   The Lagrange multipliers are defined on object.
!   if object2 is also present a "link" of object object to object
!   object2 is defined.
    integer, intent(in), optional :: object, object2

!   step specifies which nodes on curve1 are used for collocation of a
!   distributed constraint (discretization='collocation'):
!      step=0 all nodes
!      step>0 nodes 1, 1+step, 1+2*step, ...
!      step<0 except nodes 1, 1-step, 1-2*step, ...
!
!   exclude specifies which nodes on curve1 are excluded from collocation
!   of a distributed constraint (discretization='collocation'):
!      exclude=0 no excludes, all nodes are used (default)
!      exclude=1 exclude first point of curve
!      exclude=2 exclude last point of curve
!      exclude=3 exclude first and last point
!
!   NOTE: step, exclude only work for curves where the local node numbering is
!   in a natural sequence along the curve. This might not be the case for
!   curves constructed from several other curves or curves read from external
!   mesh generators.
    integer, intent(in), optional :: step, exclude

!   exclude points from a constraint on a geometry for
!   discretization='collocation'.
!   for example excludepoints=(/1,3/) excludes nodes in points P1 and P3
!   from the constraint. Note that only the nodes that are
!   actually on the geometry are excluded.
!   Note that the constraints are on the first geometry
!   (point1,curve1,surface1,volume,nodeset1)
!   and that it only makes sense that the points intersect with the geometry.
    integer, intent(in), dimension(:), optional :: excludepoints

!   exclude curves from a constraint on a geometry for
!   discretization='collocation'.
!   for example excludecurves=(/1,3/) excludes nodes on curves C1 and C3
!   from the constraint. Note that only the nodes that are
!   actually on the geometry are excluded.
!   Note that the constraints are on the first geometry
!   (point1,curve1,surface1,volume,nodeset1)
!   and that it only makes sense that the curves intersect with the geometry.
    integer, intent(in), dimension(:), optional :: excludecurves

!   exclude surfaces from a constraint on a geometry for
!   discretization='collocation'.
!   for example excludesurfaces=(/1,3/) excludes nodes on surfaces S1 and S3
!   from the constraint. Note that only the nodes that are
!   actually on the geometry are excluded.
!   Note that the constraints are on the first geometry
!   (point1,curve1,surface1,volume,nodeset1)
!   and that it only makes sense that the surfaces intersect with the geometry.
    integer, intent(in), dimension(:), optional :: excludesurfaces

!   full2=.true. means that the Lagrange multipliers (defined on curve1 or
!   surface1) are connected to _all_ degrees of freedom on curve2/surface2.
!   This can be used to build incompatible connections using mortars.
!
!   full2=.false. (default if not present) means that the Lagrange
!   multipliers are connected to the degrees of freedom on curve2/surface2
!   elementwise (discretization='weak') or
!   nodalwise (discretization='collocation').
    logical, intent(in), optional :: full2

!   for a distributed constraint (nglobalc is absent) discretization gives the
!   type of discretization of the Lagrange multiplier:
!       discretization='weak' the constraint is imposed in a weak form.
!       discretization='collocation' the constraint is imposed in a strong form
!                      in particular points: the collocation points.
!
!   for a global constraint (nglobalc is present) discretization gives the
!   type of building the constraint from smaller pieces:
!       discretization='weak' the constraint is build by element assembly
!       discretization='collocation' the constraint is build by nodal assembly.
!   note that step and exclude cannot be used here.
!
!   the default is discretization='weak', except for points and nodesets, where
!   discretization='collocation' is the only possibility
    character(len=*), intent(in), optional :: discretization

!   elementdof gives the degrees of freedom of the Lagrange multiplier in
!   the nodes of an element for discretization='weak'.
!   For example elementdof=(/1,0,1/) in a three node element.
!   NOTE: this is only for entities having no groups (geometries,objects).
    integer, intent(in), dimension(:), optional :: elementdof

!   if discretization='collocation' the number of Lagrange multipliers in
!   the nodes are given by nodedof. If nodedof=0 or absent from the heading
!   the number Lagrange multipliers is taken equal to the number of degrees of
!   freedom of the quantity under constraint.
!   Note that for objects this is a required parameter.
    integer, intent(in), optional :: nodedof

!   if naddunknowns is present there are additional unknowns in the constraint.
!   The number of unknowns is naddunknowns. This is typically used when the
!   imposed value still depends on a few parameters. For example when we
!   impose a rigid body rotation but the velocity and rotation are still
!   unknown.
    integer, intent(in), optional :: naddunknowns

!   if nglobalc is present, the constraint has global Lagrange multipliers
!   (not attached to elements or nodes but to curves, surfaces, objects etc).
!   The number of Lagrange multipliers (constraints) is nglobalc.
!
!   Although global constraints are not coupled to nodes or elements, the
!   constraints are build by smaller pieces (elements or nodes). For example,
!   if the constraint is an integral over a curve the integral is build
!   element by element (discretization=`weak') or node by node
!   (discretization=`collocation').
!
!   if nglobalc is absent the constraint is distributed (along a curve, object,
!   surface etc.)
    integer, intent(in), optional :: nglobalc

!   if present and .true. the diagonal block of the matrix is included in the
!   system matrix structure for the constraint part, including the additional
!   unknowns part.
!   Note, that it is not the full block, but the matrix elements that would
!   naturally be filled with zero when building the contraint point for point
!   (collocation) or element by element (weak).
!   default=.false.
    logical, intent(in), optional :: diagonal_block

!   if present and .true. the diagonal block of the matrix regarding the
!   the additional unknowns is included in the system matrix structure.
!   default=.false.
    logical, intent(in), optional :: diagonal_block_addunknowns

!   if present:
!   The number of Lagrangian multiplier degrees of freedom in the nodes of
!   geometries is given by
!   lmnumdegfd=0: either elementdof/nodedof or equal to the number of degrees
!                 of freedom of the constraint quantity.
!   lmnumdegfd=1: Apply additional rules in setting the number of Lagrangian
!                 multiplier degrees of freedom in the nodes of geometries in
!                 case of layers (layer1>0 and/or layer2>0), such as taking the
!                 number of degrees in a single layer only and zero degrees if
!                 layer is not present.
!   Note, that this only applies to geometries and not to objects. For objects
!   elementdof or nodedof must be used.
!   Default=1
    integer, intent(in), optional :: lmnumdegfd

!   if present: the constraint number defined.
    integer, intent(out), optional :: num

!
!   This routine can be called MAXNUMCONSTRAINTS times to set constraints
!

    logical, dimension(:), allocatable :: la1, la2
    integer :: cg, elnumnod


!   do some testing first

    call check ( mesh, 'define_constraint' )

    if ( .not. input_probdef%created ) then
      write(*,'(/2(a/))') &
        'Error in define_constraint:', &
        ' input_probdef has not been created '
      stop
    end if

    if ( present(nglobalc) ) then
!     global constraint
      if ( present(elementdof) .or. present(nodedof) &
           .or. present(step) .or. present(exclude)  &
           .or. present(excludepoints) &
           .or. present(excludecurves) .or. present(excludesurfaces) ) then
        write(*,'(2(/a))') &
          'Error: nglobalc in the heading of define_constraint', &
          'cannot be combined with the following parameters:'
        if ( present(elementdof) )     write(*,'(4x,a)') 'elementdof'
        if ( present(nodedof) )        write(*,'(4x,a)') 'nodedof'
        if ( present(step) )           write(*,'(4x,a)') 'step'
        if ( present(exclude) )        write(*,'(4x,a)') 'exclude'
        if ( present(excludepoints) )  write(*,'(4x,a)') 'excludepoints'
        if ( present(excludecurves) )  write(*,'(4x,a)') 'excludecurves'
        if ( present(excludesurfaces) )  write(*,'(4x,a)') 'excludesurfaces'
        write(*,*)
        stop
      end if
    end if

    if ( present(physq) .and. ( present(physq1) .or. present(physq2) ) ) then
      write(*,'(/2a/)') &
        'Error: physq and physq1/physq2 cannot be present at the same time', &
        ' in the heading of define_constraint'
      stop
    end if

    if ( present(physq) ) then
      if ( physq < 1 .or. physq > input_probdef%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq in heading of define_constraint has wrong value:', &
          'physq < 1 or physq > number of physical quantities = ', &
          input_probdef%nphysq
        stop
      end if
    end if

    if ( present(physq1) ) then
      if ( physq1 < 1 .or. physq1 > input_probdef%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq1 in heading of define_constraint has wrong value:', &
          'physq1 < 1 or physq1 > number of physical quantities = ', &
          input_probdef%nphysq
        stop
      end if
    end if

    if ( present(physq2) ) then
      if ( physq2 <=0 .or. physq2 > input_probdef%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq2 in heading of define_constraint has wrong value:', &
          'physq2 < 1 or physq2 > number of physical quantities = ', &
          input_probdef%nphysq
        stop
      end if
    end if

    if ( present(layer) .and. ( present(layer1) .or. present(layer2) ) ) then
      write(*,'(/2a/)') &
        'Error: layer and layer1/layer2 cannot be present at the same time', &
        ' in the heading of define_constraint'
      stop
    end if

    if ( present(layer) ) then
      if ( layer < 1 .or. layer > input_probdef%numlayers ) then
        write(*,'(/a/a,i0/)') &
          'Error: layer in heading of define_constraint has wrong value:', &
          'layer < 1 or layer > number of layers = ', &
          input_probdef%numlayers
        stop
      end if
    end if

    if ( present(layer1) ) then
      if ( layer1 < 1 .or. layer1 > input_probdef%numlayers ) then
        write(*,'(/a/a,i0/)') &
          'Error: layer1 in heading of define_constraint has wrong value:', &
          'layer1 < 1 or layer1 > number of layers = ', &
          input_probdef%numlayers
        stop
      end if
    end if

    if ( present(layer2) ) then
      if ( layer2 <=0 .or. layer2 > input_probdef%numlayers ) then
        write(*,'(/a/a,i0/)') &
          'Error: layer2 in heading of define_constraint has wrong value:', &
          'layer2 < 1 or layer2 > number of layers = ', &
          input_probdef%numlayers
        stop
      end if
    end if

    if ( present(excludepoints) ) then
      if ( any( excludepoints <=0 ) .or. &
           any( excludepoints > mesh%npoints ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: excludepoints in heading of define_constraint:', &
          ' value < 0 or value > number of points = ', &
          mesh%npoints
        stop
      end if
    end if

    if ( present(excludecurves) ) then
      if ( any( excludecurves <=0 ) .or. &
           any( excludecurves > mesh%ncurves ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: excludecurves in heading of define_constraint:', &
          ' value < 0 or value > number of curves = ', &
          mesh%ncurves
        stop
      end if
    end if

    if ( present(excludesurfaces) ) then
      if ( any( excludesurfaces <=0 ) .or. &
           any( excludesurfaces > mesh%nsurfaces ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: excludesurfaces in heading of define_constraint:', &
          ' value < 0 or value > number of surfaces = ', &
          mesh%nsurfaces
        stop
      end if
    end if

    if ( present(naddunknowns) ) then
      if ( naddunknowns < 0 ) then
        write(*,'(/2a/)') &
          'Error: naddunknowns in heading of define_constraint has a', &
          ' negative value'
        stop
      end if
    end if

    if ( present(elementdof) .and. present(nodedof) ) then
      write(*,'(/2a/)') &
        'Error: elementdof and nodedof cannot be both present in heading', &
        ' of define_constraint'
      stop
    end if

    if ( present(lmnumdegfd) ) then
      if ( .not. any ( lmnumdegfd == [0,1] ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: lmnumdegfd has wrong value:', &
          'lmnumdegfd must be either 0 or 1 = ', &
          lmnumdegfd
        stop
      end if
    end if

    if ( count( [ present(point1), present(curve1), present(surface1), &
                  present(volume), present(object), present(elementset1), &
                  present(nodeset1) ] ) > 1 ) then
      write(*,'(/2(a/))') &
        'Error: only one of point1, curve1, surface1, volume, elementset1, ', &
        ' nodeset1 or object can be present in the heading of define_constraint'
      stop
    end if

    la1 = [ present(point1), present(curve1), present(surface1), &
            present(elementset1), present(nodeset1) ]

    la2 = [ present(point2), present(curve2), present(surface2), &
            present(elementset2), present(nodeset2) ]

    if ( count(la2) > 1 ) then
      write(*,'(/2(a/))') &
        'Error: at most one of point2, curve2, surface2, elementset2,', &
        ' nodeset2 must be present in the heading of define_constraint'
      stop
    else if ( count(la2) == 1 .and. count ( la1 .and. la2 ) == 0 ) then
      write(*,'(/3(a/))') &
        'Error: second geometry (point2, curve2, surface2, elementset2,', &
        ' or nodeset2) is of a different kind than the first geometry ', &
        ' in the heading of define_constraint'
      stop
    end if

!   determine type of call

    if ( present(point1) .or. present(nodeset1) ) then

!     contraint in a point or on a nodeset

      input_probdef%numconstraints = input_probdef%numconstraints + 1
      cg = input_probdef%numconstraints

      if ( cg > size(input_probdef%constraints) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_constraint:', &
          ' constraint > maximum = ', size(input_probdef%constraints)
        stop
      end if

      if ( present(num) ) num = cg

      if ( present(point1) ) then
!       point
        if ( point1 <=0 .or. point1 > mesh%npoints ) then
          write(*,'(/a/a,i0/)') &
            'Error: point1 in heading of define_constraint has wrong value:', &
            'point1 <=0 or point1 > number of points = ', mesh%npoints
          stop
        end if
        input_probdef%constraints(cg)%typegeometry = 1
        input_probdef%constraints(cg)%geometry1 = point1
      else if ( present(nodeset1) ) then
!       surface
        if ( nodeset1 <=0 .or. nodeset1 > mesh%nnodesets ) then
          write(*,'(/a/a,i0/)') &
            'Error: nodeset1 in heading of define_constraint has wrong value:',&
            'nodeset1 <=0 or nodeset1 > number of nodesets = ', mesh%nnodesets
          stop
        end if
        input_probdef%constraints(cg)%typegeometry = 5
        input_probdef%constraints(cg)%geometry1 = nodeset1
      end if

      if ( present(nglobalc) ) then
!       global constraint
        input_probdef%constraints(cg)%typeconstraint = 2
        input_probdef%constraints(cg)%nglobalc = nglobalc
      else
!       distributed constraint
        input_probdef%constraints(cg)%typeconstraint = 1
        input_probdef%constraints(cg)%nglobalc = 0
      end if

      input_probdef%constraints(cg)%physq1 = 0
      input_probdef%constraints(cg)%physq2 = 0

      if ( present(physq) ) then
        input_probdef%constraints(cg)%physq1 = physq
        input_probdef%constraints(cg)%physq2 = physq
      end if

      if ( present(physq1) ) then
        input_probdef%constraints(cg)%physq1 = physq1
      end if

      if ( present(physq2) ) then
        input_probdef%constraints(cg)%physq2 = physq2
      end if

      input_probdef%constraints(cg)%layer1 = 0
      input_probdef%constraints(cg)%layer2 = 0

      if ( present(layer) ) then
        input_probdef%constraints(cg)%layer1 = layer
        input_probdef%constraints(cg)%layer2 = layer
      end if

      if ( present(layer1) ) then
        input_probdef%constraints(cg)%layer1 = layer1
      end if

      if ( present(layer2) ) then
        input_probdef%constraints(cg)%layer2 = layer2
      end if

      input_probdef%constraints(cg)%lmnumdegfd = &
                           set_optional( variable=lmnumdegfd, default=1 )

      if ( present(naddunknowns) ) then
        input_probdef%constraints(cg)%naddunknowns = naddunknowns
      else
        input_probdef%constraints(cg)%naddunknowns = 0
      end if

      if ( present(discretization) ) then
        if ( discretization == 'weak' ) then
          write(*,'(2(/a)/)') &
            'Error in the heading of define_constraint:', &
            ' discretization=''weak'' invalid for points or nodesets '
          stop
        else if ( discretization == 'collocation' ) then
          input_probdef%constraints(cg)%discretization = 1
        else
          write(*,'(2(/a)/)') &
            'Error: discretization in the heading of define_constraint', &
            'must have the value ''weak'' or ''collocation'''
          stop
        end if
      else
!       default: collocation constraint
        input_probdef%constraints(cg)%discretization = 1
      end if

      if ( present(elementdof) ) then
        write(*,'(2(/a)/)') &
          'Error: elementdof in the heading of define_constraint', &
          '  whereas discretization == ''collocation'''
        stop
      end if
      allocate (input_probdef%constraints(cg)%elnumdegfd(0) )

      if ( present(nodedof) ) then
        if ( nodedof < 0 ) then
          write(*,'(/2a/)') &
            'Error: nodedof in heading of define_constraint has a', &
            ' negative value'
          stop
        end if
        input_probdef%constraints(cg)%nodenumdegfd = nodedof
      else
        input_probdef%constraints(cg)%nodenumdegfd = 0
      end if

      if ( present(excludepoints) ) then
        input_probdef%constraints(cg)%excludepoints = excludepoints
      else
        allocate (input_probdef%constraints(cg)%excludepoints(0) )
      end if

      if ( present(excludecurves) ) then
        input_probdef%constraints(cg)%excludecurves = excludecurves
      else
        allocate (input_probdef%constraints(cg)%excludecurves(0) )
      end if

      if ( present(excludesurfaces) ) then
        input_probdef%constraints(cg)%excludesurfaces = excludesurfaces
      else
        allocate (input_probdef%constraints(cg)%excludesurfaces(0) )
      end if

      if ( present(point2) .and. present(point1) ) then

        if ( point2 == point1 .or. point2 <=0 &
               .or. point2 > mesh%npoints ) then
          write(*,'(/a/a,i0/)') &
            'Error: point2 in heading of define_constraint has wrong value:', &
            'point2 == point1 or point2 <=0 or point2 > number of points = ', &
            mesh%npoints
          stop
        end if

        input_probdef%constraints(cg)%geometry2 = point2
        if ( present(full2) ) then
          write(*,'(/2a/)') &
            'Error: in heading of define_constraint:', &
            ' full2 not allowed for points '
        end if
        input_probdef%constraints(cg)%full2 = .false.

      else if ( present(nodeset2) .and. present(nodeset1) ) then

        if ( nodeset2 == nodeset1 .or. nodeset2 <=0 &
               .or. nodeset2 > mesh%nnodesets ) then
          write(*,'(/a/2a,i0/)') &
            'Error: nodeset2 in heading of define_constraint has wrong value:',&
            'nodeset2 == nodeset1 or nodeset2 <=0 or', &
            ' nodeset2 > number of nodesets = ', mesh%nnodesets
          stop
        end if

        input_probdef%constraints(cg)%geometry2 = nodeset2
        if ( present(full2) ) then
          input_probdef%constraints(cg)%full2 = full2
        else
          input_probdef%constraints(cg)%full2 = .false.
        end if

        if ( .not. input_probdef%constraints(cg)%full2 ) then
          if ( size(mesh%nodesets(nodeset1)%a) /= &
               size(mesh%nodesets(nodeset2)%a) ) then
            write(*,'(/2a/)') &
              'Error: in heading of define_constraint:', &
              ' number of nodes in nodeset1 /= number of nodes in nodeset2'
            stop
          end if
        end if

      else
        if ( present(full2) ) then
          write(*,'(2(/a)/)') &
            'Error: full2 in the heading of define_constraint', &
            'whereas no curve2, surface2 or nodeset2 is defined'
          stop
        end if
        input_probdef%constraints(cg)%geometry2 = 0
        input_probdef%constraints(cg)%full2 = .false.
      end if

!     diagonal blocks

      input_probdef%constraints(cg)%diagonal_block = &
           set_optional( variable=diagonal_block, default=.false. )
      input_probdef%constraints(cg)%diagonal_block_addunknowns = &
           set_optional( variable=diagonal_block_addunknowns, default=.false. )

    else if ( present(curve1) .or. present(surface1) .or. present(volume) ) then

!     contraint on a curve, surface or volume

      input_probdef%numconstraints = input_probdef%numconstraints + 1
      cg = input_probdef%numconstraints

      if ( cg > size(input_probdef%constraints) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_constraint:', &
          ' constraint > maximum = ', size(input_probdef%constraints)
        stop
      end if

      if ( present(num) ) num = cg

      if ( present(curve1) ) then
!       curve
        if ( curve1 <=0 .or. curve1 > mesh%ncurves ) then
          write(*,'(/a/a,i0/)') &
            'Error: curve1 in heading of define_constraint has wrong value:', &
            'curve1 <=0 or curve1 > number of curves = ', mesh%ncurves
          stop
        end if
        input_probdef%constraints(cg)%typegeometry = 2
        input_probdef%constraints(cg)%geometry1 = curve1
        elnumnod = mesh%curves(curve1)%elnumnod
      else if ( present(surface1) ) then
!       surface
        if ( surface1 <=0 .or. surface1 > mesh%nsurfaces ) then
          write(*,'(/a/a,i0/)') &
            'Error: surface1 in heading of define_constraint has wrong value:',&
            'surface1 <=0 or surface1 > number of surfaces = ', mesh%nsurfaces
          stop
        end if
        input_probdef%constraints(cg)%typegeometry = 3
        input_probdef%constraints(cg)%geometry1 = surface1
        elnumnod = mesh%surfaces(surface1)%elnumnod
      else if ( present(volume) ) then
!       volume
        if ( volume <=0 .or. volume > mesh%nvolumes ) then
          write(*,'(/a/a,i0/)') &
            'Error: volume in heading of define_constraint has wrong value:',&
            'volume <=0 or volume > number of volumes = ', mesh%nvolumes
          stop
        end if
        input_probdef%constraints(cg)%typegeometry = 4
        input_probdef%constraints(cg)%geometry1 = volume
        elnumnod = mesh%volumes(volume)%elnumnod
      end if

      if ( present(nglobalc) ) then
!       global constraint
        input_probdef%constraints(cg)%typeconstraint = 2
        input_probdef%constraints(cg)%nglobalc = nglobalc
      else
!       distributed constraint
        input_probdef%constraints(cg)%typeconstraint = 1
        input_probdef%constraints(cg)%nglobalc = 0
      end if

      input_probdef%constraints(cg)%physq1 = 0
      input_probdef%constraints(cg)%physq2 = 0

      if ( present(physq) ) then
        input_probdef%constraints(cg)%physq1 = physq
        input_probdef%constraints(cg)%physq2 = physq
      end if

      if ( present(physq1) ) then
        input_probdef%constraints(cg)%physq1 = physq1
      end if

      if ( present(physq2) ) then
        input_probdef%constraints(cg)%physq2 = physq2
      end if

      input_probdef%constraints(cg)%layer1 = 0
      input_probdef%constraints(cg)%layer2 = 0

      if ( present(layer) ) then
        input_probdef%constraints(cg)%layer1 = layer
        input_probdef%constraints(cg)%layer2 = layer
      end if

      if ( present(layer1) ) then
        input_probdef%constraints(cg)%layer1 = layer1
      end if

      if ( present(layer2) ) then
        input_probdef%constraints(cg)%layer2 = layer2
      end if

      input_probdef%constraints(cg)%lmnumdegfd = &
                           set_optional( variable=lmnumdegfd, default=1 )

      if ( present(naddunknowns) ) then
        input_probdef%constraints(cg)%naddunknowns = naddunknowns
      else
        input_probdef%constraints(cg)%naddunknowns = 0
      end if

      if ( present(discretization) ) then
        if ( discretization == 'weak' ) then
          input_probdef%constraints(cg)%discretization = 0
        else if ( discretization == 'collocation' ) then
          input_probdef%constraints(cg)%discretization = 1
        else
          write(*,'(2(/a)/)') &
            'Error: discretization in the heading of define_constraint', &
            'must have the value ''weak'' or ''collocation'''
          stop
        end if
      else
!       default: weak (integral) constraint
        input_probdef%constraints(cg)%discretization = 0
      end if

      if ( present(elementdof) ) then
        if ( input_probdef%constraints(cg)%discretization == 1 ) then
          write(*,'(2(/a)/)') &
            'Error: elementdof in the heading of define_constraint', &
            'whereas discretization == ''collocation'''
          stop
        end if
        if ( size(elementdof) /= elnumnod ) then
          write(*,'(/2a,i0/a,i0/)') &
            'Error: elementdof in the heading of define_constraint has', &
            ' wrong size: ', size(elementdof), &
            ' whereas it should be ', elnumnod
          stop
        end if
        if ( any(elementdof < 0) ) then
          write(*,'(/2(a/))') &
            'Error: elementdof in heading of define_constraint has ', &
            ' one or more negative values'
          stop
        end if
        input_probdef%constraints(cg)%elnumdegfd = elementdof
      else
        if ( input_probdef%constraints(cg)%discretization == 0 .and. &
             input_probdef%constraints(cg)%typeconstraint == 1 ) then
          write(*,'(2(/a)/)') &
            'Error: elementdof not in the heading of define_constraint', &
            'whereas discretization == ''weak'''
          stop
        end if
        allocate (input_probdef%constraints(cg)%elnumdegfd(0) )
      end if

      if ( present(nodedof) ) then
        if ( input_probdef%constraints(cg)%discretization == 0 ) then
          write(*,'(2(/a)/)') &
            'Error: nodedof in the heading of define_constraint', &
            'whereas discretization == ''weak'''
          stop
        end if
        if ( nodedof < 0 ) then
          write(*,'(/2a/)') &
            'Error: nodedof in heading of define_constraint has a', &
            ' negative value'
          stop
        end if
        input_probdef%constraints(cg)%nodenumdegfd = nodedof
      else
        input_probdef%constraints(cg)%nodenumdegfd = 0
      end if

      if ( present(step) .and. present(curve1) ) then
        if ( input_probdef%constraints(cg)%discretization == 0 ) then
          write(*,'(2(/a)/)') &
            'Error: step in the heading of define_constraint', &
            'whereas discretization == ''weak'''
          stop
        end if
        input_probdef%constraints(cg)%step = step
      end if

      if ( present(exclude) .and. present(curve1) ) then
        if ( input_probdef%constraints(cg)%discretization == 0 ) then
          write(*,'(2(/a)/)') &
            'Error: exclude in the heading of define_constraint', &
            'whereas discretization == ''weak'''
          stop
        end if
        input_probdef%constraints(cg)%exclude = exclude
      end if

      if ( present(excludepoints) ) then
        if ( input_probdef%constraints(cg)%discretization == 0 ) then
          write(*,'(2(/a)/)') &
            'Error: excludepoints in the heading of define_constraint', &
            'whereas discretization == ''weak'''
          stop
        end if
        input_probdef%constraints(cg)%excludepoints = excludepoints
      else
        allocate (input_probdef%constraints(cg)%excludepoints(0) )
      end if

      if ( present(excludecurves) ) then
        if ( input_probdef%constraints(cg)%discretization == 0 ) then
          write(*,'(2(/a)/)') &
            'Error: excludecurves in the heading of define_constraint', &
            'whereas discretization == ''weak'''
          stop
        end if
        input_probdef%constraints(cg)%excludecurves = excludecurves
      else
        allocate (input_probdef%constraints(cg)%excludecurves(0) )
      end if

      if ( present(excludesurfaces) ) then
        if ( input_probdef%constraints(cg)%discretization == 0 ) then
          write(*,'(2(/a)/)') &
            'Error: excludesurfaces in the heading of define_constraint', &
            'whereas discretization == ''weak'''
          stop
        end if
        input_probdef%constraints(cg)%excludesurfaces = excludesurfaces
      else
        allocate (input_probdef%constraints(cg)%excludesurfaces(0) )
      end if

      if ( present(curve2) .and. present(curve1) ) then

        if ( curve2 == curve1 .or. curve2 <=0 &
               .or. curve2 > mesh%ncurves ) then
          write(*,'(/a/a,i0/)') &
            'Error: curve2 in heading of define_constraint has wrong value:', &
            'curve2 == curve1 or curve2 <=0 or curve2 > number of curves = ', &
            mesh%ncurves
          stop
        end if

        input_probdef%constraints(cg)%geometry2 = curve2
        if ( present(full2) ) then
          input_probdef%constraints(cg)%full2 = full2
        else
          input_probdef%constraints(cg)%full2 = .false.
        end if

        if ( .not. input_probdef%constraints(cg)%full2 ) then
          if ( mesh%curves(curve1)%nnodes /= mesh%curves(curve2)%nnodes ) then
            write(*,'(/2a/)') &
              'Error: in heading of define_constraint:', &
              ' number of nodes on curve1 /= number of nodes on curve2'
            stop
          end if
        end if

      else if ( present(surface2) .and. present(surface1) ) then

        if ( surface2 == surface1 .or. surface2 <=0 &
               .or. surface2 > mesh%nsurfaces ) then
          write(*,'(/a/2a,i0/)') &
            'Error: surface2 in heading of define_constraint has wrong value:',&
            'surface2 == surface1 or surface2 <=0 or', &
            ' surface2 > number of surfaces = ', mesh%nsurfaces
          stop
        end if

        input_probdef%constraints(cg)%geometry2 = surface2
        if ( present(full2) ) then
          input_probdef%constraints(cg)%full2 = full2
        else
          input_probdef%constraints(cg)%full2 = .false.
        end if

        if ( .not. input_probdef%constraints(cg)%full2 ) then
          if ( mesh%surfaces(surface1)%nnodes /= &
               mesh%surfaces(surface2)%nnodes ) then
            write(*,'(/2a/)') &
              'Error: in heading of define_constraint:', &
              ' number of nodes on surface1 /= number of nodes on surface2'
            stop
          end if
        end if

      else
        if ( present(full2) ) then
          write(*,'(2(/a)/)') &
            'Error: full2 in the heading of define_constraint', &
            'whereas no curve2 or surface2 is defined'
          stop
        end if
        input_probdef%constraints(cg)%geometry2 = 0
        input_probdef%constraints(cg)%full2 = .false.
      end if

!     diagonal blocks

      input_probdef%constraints(cg)%diagonal_block = &
           set_optional( variable=diagonal_block, default=.false. )
      input_probdef%constraints(cg)%diagonal_block_addunknowns = &
           set_optional( variable=diagonal_block_addunknowns, default=.false. )

    else if ( present(elementset1) ) then

!     constraint on an elementset

      input_probdef%numconstraints = input_probdef%numconstraints + 1
      cg = input_probdef%numconstraints

      if ( cg > size(input_probdef%constraints) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_constraint:', &
          ' constraint > maximum = ', size(input_probdef%constraints)
        stop
      end if

      if ( present(num) ) num = cg

      if ( elementset1 <=0 .or. elementset1 > mesh%nelementsets ) then
        write(*,'(/2a/a,i0/)') &
          'Error: elementset1 in heading of define_constraint has ', &
          'wrong value:', &
          'elementset1 <=0 or elementset1 > number of elementsets = ', &
          mesh%nelementsets
        stop
      end if

      input_probdef%constraints(cg)%elementset1 = elementset1

      if ( present(nglobalc) ) then
!       global constraint
        input_probdef%constraints(cg)%typeconstraint = 2
        input_probdef%constraints(cg)%nglobalc = nglobalc
      else
!       distributed constraint
        input_probdef%constraints(cg)%typeconstraint = 1
        input_probdef%constraints(cg)%nglobalc = 0
      end if

      input_probdef%constraints(cg)%physq1 = 0
      input_probdef%constraints(cg)%physq2 = 0

      if ( present(physq) ) then
        input_probdef%constraints(cg)%physq1 = physq
        input_probdef%constraints(cg)%physq2 = physq
      end if

      if ( present(physq1) ) then
        input_probdef%constraints(cg)%physq1 = physq1
      end if

      if ( present(physq2) ) then
        input_probdef%constraints(cg)%physq2 = physq2
      end if

      input_probdef%constraints(cg)%layer1 = 0
      input_probdef%constraints(cg)%layer2 = 0

      if ( present(layer) ) then
        input_probdef%constraints(cg)%layer1 = layer
        input_probdef%constraints(cg)%layer2 = layer
      end if

      if ( present(layer1) ) then
        input_probdef%constraints(cg)%layer1 = layer1
      end if

      if ( present(layer2) ) then
        input_probdef%constraints(cg)%layer2 = layer2
      end if

      input_probdef%constraints(cg)%lmnumdegfd = &
                           set_optional( variable=lmnumdegfd, default=1 )

      if ( present(naddunknowns) ) then
        input_probdef%constraints(cg)%naddunknowns = naddunknowns
      else
        input_probdef%constraints(cg)%naddunknowns = 0
      end if

      if ( present(discretization) ) then
        if ( discretization == 'weak' ) then
          input_probdef%constraints(cg)%discretization = 0
        else if ( discretization == 'collocation' ) then
          input_probdef%constraints(cg)%discretization = 1
        else
          write(*,'(2(/a)/)') &
            'Error: discretization in the heading of define_constraint', &
            'must have the value ''weak'' or ''collocation'''
          stop
        end if
      else
!       default: weak (integral) constraint
        input_probdef%constraints(cg)%discretization = 0
      end if

      if ( input_probdef%constraints(cg)%typeconstraint == 1 ) then
        write(*,'(5(/a)/)') &
          'Error: in define_constraint', &
          'Distributed constraints on elementsets have not yet been implemented'
        stop
      end if

!     elnumdegfd not set
      allocate (input_probdef%constraints(cg)%elnumdegfd(0) )
!     nodedof not set and of no use yet
      input_probdef%constraints(cg)%nodenumdegfd = 0

      if ( present(excludepoints) .or. present(excludecurves) &
           .or. present(excludesurfaces) ) then
        write(*,'(/2a/a/)') &
          'Warning: excludepoints/excludecurves/excludesurfaces ', &
          'in the heading of define_constraint.', &
          'Ignored for constraint on an elementset.'
      end if

      allocate (input_probdef%constraints(cg)%excludepoints(0) )
      allocate (input_probdef%constraints(cg)%excludecurves(0) )
      allocate (input_probdef%constraints(cg)%excludesurfaces(0) )

      if ( present(elementset2) .and. present(elementset1) ) then

        if ( elementset2 == elementset1 .or. elementset2 <=0 &
               .or. elementset2 > mesh%nelementsets ) then
          write(*,'(/2a/2a,i0/)') &
            'Error: elementset2 in heading of define_constraint ', &
            'has wrong value:', &
            'elementset2 == elementset1 or elementset2 <=0 or ', &
            ' elementset2 > number of elementsets = ', &
            mesh%nelementsets
          stop
        end if

        input_probdef%constraints(cg)%elementset2 = elementset2

        if ( mesh%elementsets(elementset1)%nnodes /= &
                    mesh%elementsets(elementset2)%nnodes ) then
          write(*,'(/2a/)') &
            'Error: in heading of define_constraint:', &
            ' number of nodes on elementset1 /= number of nodes on elementset2'
          stop
        end if

        if ( mesh%elementsets(elementset1)%nelem /= &
                    mesh%elementsets(elementset2)%nelem ) then
          write(*,'(/a/2a/)') &
            'Error: in heading of define_constraint:', &
            ' number of elements in elementset1 /= number of ', &
            'elements in elementset2'
          stop
        end if

      else
        input_probdef%constraints(cg)%elementset2 = 0
      end if

      input_probdef%constraints(cg)%full2 = .false.

!     diagonal blocks

      input_probdef%constraints(cg)%diagonal_block = &
           set_optional( variable=diagonal_block, default=.false. )
      input_probdef%constraints(cg)%diagonal_block_addunknowns = &
           set_optional( variable=diagonal_block_addunknowns, default=.false. )

    else if ( present(object) ) then

!     contraint on an object

      if ( object <=0 .or. object > mesh%nobjects ) then
        write(*,'(/a/a,i0/)') &
          'Error: object in heading of define_constraint has wrong value:', &
          'object <=0 or object > number of objects = ', mesh%nobjects
        stop
      end if

      input_probdef%numconstraints = input_probdef%numconstraints + 1
      cg = input_probdef%numconstraints

      if ( cg > size(input_probdef%constraints) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_constraint:', &
          ' constraint > maximum = ', size(input_probdef%constraints)
        stop
      end if

      if ( present(num) ) num = cg

      input_probdef%constraints(cg)%object = object

      if ( present(nglobalc) ) then
!       global constraint
        input_probdef%constraints(cg)%typeconstraint = 2
        input_probdef%constraints(cg)%nglobalc = nglobalc
      else
!       distributed constraint
        input_probdef%constraints(cg)%typeconstraint = 1
        input_probdef%constraints(cg)%nglobalc = 0
      end if

      input_probdef%constraints(cg)%physq1 = 0
      input_probdef%constraints(cg)%physq2 = 0

      if ( present(physq) ) then
        input_probdef%constraints(cg)%physq1 = physq
        input_probdef%constraints(cg)%physq2 = physq
      end if

      if ( present(physq1) ) then
        input_probdef%constraints(cg)%physq1 = physq1
      end if

      if ( present(physq2) ) then
        input_probdef%constraints(cg)%physq2 = physq2
      end if

      input_probdef%constraints(cg)%layer1 = 0
      input_probdef%constraints(cg)%layer2 = 0

      if ( present(layer) ) then
        input_probdef%constraints(cg)%layer1 = layer
        input_probdef%constraints(cg)%layer2 = layer
      end if

      if ( present(layer1) ) then
        input_probdef%constraints(cg)%layer1 = layer1
      end if

      if ( present(layer2) ) then
        input_probdef%constraints(cg)%layer2 = layer2
      end if

      input_probdef%constraints(cg)%lmnumdegfd = &
                           set_optional( variable=lmnumdegfd, default=1 )

      if ( present(naddunknowns) ) then
        input_probdef%constraints(cg)%naddunknowns = naddunknowns
      else
        input_probdef%constraints(cg)%naddunknowns = 0
      end if

      if ( present(discretization) ) then
        if ( discretization == 'weak' ) then
          input_probdef%constraints(cg)%discretization = 0
        else if ( discretization == 'collocation' ) then
          input_probdef%constraints(cg)%discretization = 1
        else
          write(*,'(2(/a)/)') &
            'Error: discretization in the heading of define_constraint ', &
            'must have the value ''weak'' or ''collocation'''
          stop
        end if
      else
!       default: weak (integral) constraint
        input_probdef%constraints(cg)%discretization = 0
      end if

      if ( input_probdef%constraints(cg)%discretization == 0 .and. &
           .not. mesh%objects(object)%intpoints ) then
        write(*,'(/a/a,i0/a/)') &
          'Error define_constraint: ', &
          ' a weak constraint has been defined on object ', object, &
          ' however it does not contain integration points.'
        stop
      end if

      if ( present(elementdof) ) then
        if ( input_probdef%constraints(cg)%discretization == 1 ) then
          write(*,'(2(/a)/)') &
            'Error: elementdof in the heading of define_constraint', &
            'whereas discretization == ''collocation'''
          stop
        end if
        if ( size(elementdof) /= mesh%objects(object)%elnumnod ) then
          write(*,'(/2a,i0/a,i0/)') &
            'Error: elementdof in the heading of define_constraint has', &
            ' wrong size: ', size(elementdof), &
            ' whereas it should be ', mesh%objects(object)%elnumnod
          stop
        end if
        if ( any(elementdof < 0) ) then
          write(*,'(/2(a/))') &
            'Error: elementdof in heading of define_constraint has ', &
            ' one or more negative values'
          stop
        end if
        input_probdef%constraints(cg)%elnumdegfd = elementdof
      else
        if ( input_probdef%constraints(cg)%discretization == 0 .and. &
             input_probdef%constraints(cg)%typeconstraint == 1 ) then
          write(*,'(2(/a)/)') &
            'Error: elementdof not in the heading of define_constraint', &
            'whereas discretization == ''weak'''
          stop
        end if
        allocate (input_probdef%constraints(cg)%elnumdegfd(0) )
      end if

      if ( present(nodedof) ) then
        if ( input_probdef%constraints(cg)%discretization == 0 ) then
          write(*,'(2(/a)/)') &
            'Error: nodedof in the heading of define_constraint', &
            'whereas discretization == ''weak'''
          stop
        end if
        if ( nodedof < 0 ) then
          write(*,'(/2a/)') &
            'Error: nodedof in heading of define_constraint has a', &
            ' negative value'
          stop
        end if
        input_probdef%constraints(cg)%nodenumdegfd = nodedof
      else
        if ( input_probdef%constraints(cg)%discretization == 1 .and. &
             input_probdef%constraints(cg)%typeconstraint == 1 ) then
          write(*,'(2(/a)/)') &
            'Error: nodedof not in the heading of define_constraint', &
            'for objects, whereas discretization == ''collocation'''
          stop
        end if
        input_probdef%constraints(cg)%nodenumdegfd = 0
      end if

      allocate (input_probdef%constraints(cg)%excludepoints(0) )
      allocate (input_probdef%constraints(cg)%excludecurves(0) )
      allocate (input_probdef%constraints(cg)%excludesurfaces(0) )

      input_probdef%constraints(cg)%full2 = .false.

      if ( present(object2) ) then

!       connection object

        if ( object2 <=0 .or. object2 > mesh%nobjects ) then
          write(*,'(/a/a,i0/)') &
            'Error: object2 in heading of define_constraint has wrong value:', &
            'object2 <=0 or object2 > number of objects = ', mesh%nobjects
          stop
        end if

        if ( mesh%objects(object)%typeofobject == 2 ) then
          write(*,'(/a/a/)') &
            'Error: object in heading of define_constraint is two-sided,', &
            'whereas also object2 is present: this is incompatible.'
          stop
        end if

        if ( input_probdef%constraints(cg)%discretization == 0 ) then
!         some testing for weak constraints
          if ( mesh%objects(object)%nelem /= mesh%objects(object2)%nelem .or. &
               mesh%objects(object)%elnumnod /= mesh%objects(object2)%elnumnod &
               .or. mesh%objects(object)%element%globalshape /= &
                    mesh%objects(object2)%element%globalshape ) then
            write(*,'(/a/a,i0/a/a,i0,a/)') &
              'Error define_constraint: ', &
              ' a weak constraint has been defined on object ', object, &
              ' however the number of elements or element type in the', &
              ' connecting object ', object2, ' is not the same.'
            stop
          end if
          if ( .not. mesh%objects(object2)%intpoints ) then
            write(*,'(/a/a,i0/a,i0,a/)') &
              'Error define_constraint: ', &
              ' a weak constraint has been defined on object ', object, &
              ' however the connecting object ', object2, &
              ' does not contain integration points.'
            stop
          end if
          if ( mesh%objects(object)%ninti /= mesh%objects(object2)%ninti  ) then
            write(*,'(/a/a,i0/a/a,i0,a/)') &
              'Error define_constraint: ', &
              ' a weak constraint has been defined on object ', object, &
              ' however the number of integration points in the', &
              ' connecting object ', object2, ' is not the same.'
            stop
          end if
        end if

        input_probdef%constraints(cg)%object2 = object2

        if ( mesh%objects(object)%nnodes /= &
             mesh%objects(object2)%nnodes ) then
          write(*,'(/2a/)') &
            'Error: in heading of define_constraint:', &
            ' number of nodes on object /= number of nodes on object2'
          stop
        end if

      else

        input_probdef%constraints(cg)%object2 = 0

      end if

!     diagonal blocks

      input_probdef%constraints(cg)%diagonal_block = &
           set_optional( variable=diagonal_block, default=.false. )
      input_probdef%constraints(cg)%diagonal_block_addunknowns = &
           set_optional( variable=diagonal_block_addunknowns, default=.false. )

!     some testing of keywords
      if ( present(step) .or. present(exclude) .or. present(excludecurves) &
         .or. present(curve2) .or. present(surface2) .or. present(full2) &
         .or. present(excludesurfaces) .or. present(elementset2) &
         .or. present(nodeset2) .or. present(point2) &
         .or.  present(excludepoints) ) then
        write(*,'(2(/a))') &
          'Error: object in the heading of define_constraint', &
          'cannot be combined with the following parameters:'
        if ( present(step) )           write(*,'(4x,a)') 'step'
        if ( present(exclude) )        write(*,'(4x,a)') 'exclude'
        if ( present(excludepoints) )  write(*,'(4x,a)') 'excludepoints'
        if ( present(excludecurves) )  write(*,'(4x,a)') 'excludecurves'
        if ( present(excludesurfaces) )  write(*,'(4x,a)') 'excludesurfaces'
        if ( present(point2) )         write(*,'(4x,a)') 'point2'
        if ( present(curve2) )         write(*,'(4x,a)') 'curve2'
        if ( present(surface2) )       write(*,'(4x,a)') 'surface2'
        if ( present(elementset2) )    write(*,'(4x,a)') 'elementset2'
        if ( present(nodeset2) )       write(*,'(4x,a)') 'nodeset2'
        if ( present(full2) )          write(*,'(4x,a)') 'full2'
        write(*,*)
        stop
      end if

    else

      write(*,'(7(/a)/)') &
        'Error: could not determine type of constraint from the heading', &
        'of define_constraint. None of the keywords:',      &
        '  point1 ',                                                       &
        '  curve1 ',                                                       &
        '  surface1 ',                                                     &
        '  volume ',                                                       &
        '  elementset1 ',                                                  &
        '  nodeset1 ',                                                     &
        '  object ',                                                       &
        'is present.'
      stop

    end if

  end subroutine define_constraint


! Definition of connections

  subroutine define_connection ( mesh, input_probdef, point1, point2, curve1, &
    curve2, surface1, surface2, volume1, volume2, elementset1, elementset2, &
    nodeset1, nodeset2, object, object2, physq, physq1, physq2, layer, &
    layer1, layer2, discretization, num )

    type(mesh_t), intent(in) :: mesh

!   structure to store all definitions for the problem definition
    type(input_probdef_t), intent(inout) :: input_probdef

!   if both are present a connection between the points is defined.
!   Only discretization='collocation' is available for points.
    integer, intent(in), optional :: point1, point2

!   If both curve1, curve2 are present a connection between the curves is
!   defined. If only curve1 is present, elementset2 or nodeset2 must be present.
!   If only curve2 is present, object must be present.
!   Note that the definition of the curves is important: they are
!   connected according to either the node sequence (collocation) or
!   element sequence (weak).
    integer, intent(in), optional :: curve1, curve2

!   If both surface1, surface2 are present a connection between the surfaces is
!   defined. If only surface1 is present, elementset2 or nodeset2 must be
!   present. If only surface2 is present, object must be present.
!   Note that the definition of the surfaces is important: they are
!   connected according to either the node sequence (collocation) or
!   element sequence (weak).
    integer, intent(in), optional :: surface1, surface2

!   If both volume1, volume2 are present a connection between the volumes is
!   defined. If only volume1 is present, elementset2 or nodeset2 must be
!   present. If only volume2 is present, object must be present.
!   Note that the definition of the volumes is important: they are
!   connected according to either the node sequence (collocation) or
!   element sequence (weak).
    integer, intent(in), optional :: volume1, volume2

!   If both elementset1, elementset2 are present a connection between the
!   elementsets is defined. If only elementset2 is present, curve1, surface1,
!   volume1 or object must be present. Only elementset1 present is not allowed.
!   Note that the definition of the elementsets is important: they are
!   connected according to either the node sequence (collocation) or
!   element sequence (weak).
    integer, intent(in), optional :: elementset1, elementset2

!   If both nodeset1, nodeset2 are present a connection between the
!   nodesets is defined. If only nodeset2 is present, curve1, surface1,
!   volume1 or object must be present.
!   Note that the definition of the nodesets is important: they are connected
!   node for node.
!   Only discretization='collocation' is available for nodesets.
    integer, intent(in), optional :: nodeset1, nodeset2

!   object is the object on which the connection is defined for a
!   two-sided object or the first object for a connection.
    integer, intent(in), optional :: object

!   if object2 is also present a connection of object to object2 is defined
    integer, intent(in), optional :: object2

!   if physq is present physq denotes the physical quantity that must be
!   connected otherwise all nodal degrees of freedom are considered.
    integer, intent(in), optional :: physq

!   if physq1 and/or physq2 is present the physical quantity that must be
!   connected can be different for the "two sides".
    integer, intent(in), optional :: physq1, physq2

!   if layer is present layer denotes the layer that must be
!   connected otherwise all nodal degrees of freedom are considered.
!   NOTE: only a connection to elements that have all nodes within a layer is
!   supported.
    integer, intent(in), optional :: layer

!   if layer1 and/or layer2 is present the layer that must be
!   connected can be different for the "two sides".
!   NOTE: only a connection to elements that have all nodes within a layer is
!   supported.
    integer, intent(in), optional :: layer1, layer2

!   discretization gives the type of discretization of connection
!       discretization='weak' the connection is imposed in a weak form. The
!                      assembling is done element by element.
!       discretization='collocation' the connection is imposed in a strong form
!                      in particular points: the collocation points.
!   the default is discretization='weak'
    character(len=*), intent(in), optional :: discretization

!   if present: the connection number defined.
    integer, intent(out), optional :: num

!   WARNING:
!   It is assumed that the connections are made between degrees of freedom that
!   are outside the scope of the standard FEM structure.  For example, it is not
!   allowed to make connections between elements that have common nodes,
!   because then there is already a connection because of the standard
!   FEM structure.
!
!   This routine can be called MAXNUMCONNECTIONS times to set connections
!

    integer :: cg


!   do some testing first

    call check ( mesh, 'define_connection' )

    if ( .not. input_probdef%created ) then
      write(*,'(/2(a/))') &
        'Error in define_connection:', &
        ' input_probdef has not been created '
      stop
    end if

    if ( present(physq) .and. ( present(physq1) .or. present(physq2) ) ) then
      write(*,'(/2a/)') &
        'Error: physq and physq1/physq2 cannot be present at the same time', &
        ' in the heading of define_connection'
      stop
    end if

    if ( present(physq) ) then
      if ( physq < 1 .or. physq > input_probdef%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq in heading of define_connection has wrong value:', &
          'physq < 1 or physq > number of physical quantities = ', &
          input_probdef%nphysq
        stop
      end if
    end if

    if ( present(physq1) ) then
      if ( physq1 < 1 .or. physq1 > input_probdef%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq1 in heading of define_connection has wrong value:', &
          'physq1 < 1 or physq1 > number of physical quantities = ', &
          input_probdef%nphysq
        stop
      end if
    end if

    if ( present(physq2) ) then
      if ( physq2 <=0 .or. physq2 > input_probdef%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq2 in heading of define_connection has wrong value:', &
          'physq2 < 1 or physq2 > number of physical quantities = ', &
          input_probdef%nphysq
        stop
      end if
    end if

    if ( present(layer) .and. ( present(layer1) .or. present(layer2) ) ) then
      write(*,'(/2a/)') &
        'Error: layer and layer1/layer2 cannot be present at the same time', &
        ' in the heading of define_connection'
      stop
    end if

    if ( present(layer) ) then
      if ( layer < 1 .or. layer > input_probdef%numlayers ) then
        write(*,'(/a/a,i0/)') &
          'Error: layer in heading of define_connection has wrong value:', &
          'layer < 1 or layer > number of layers = ', &
          input_probdef%numlayers
        stop
      end if
    end if

    if ( present(layer1) ) then
      if ( layer1 < 1 .or. layer1 > input_probdef%numlayers ) then
        write(*,'(/a/a,i0/)') &
          'Error: layer1 in heading of define_connection has wrong value:', &
          'layer1 < 1 or layer1 > number of layers = ', &
          input_probdef%numlayers
        stop
      end if
    end if

    if ( present(layer2) ) then
      if ( layer2 <=0 .or. layer2 > input_probdef%numlayers ) then
        write(*,'(/a/a,i0/)') &
          'Error: layer2 in heading of define_connection has wrong value:', &
          'layer2 < 1 or layer2 > number of layers = ', &
          input_probdef%numlayers
        stop
      end if
    end if

    if ( count( [ present(point1), present(curve1), present(surface1), &
                  present(volume1), present(elementset1), &
                  present(nodeset1), present(object) ] ) /= 1 ) then
      write(*,'(/2(a/))') &
        'Error: one of point1, curve1, surface1, elementset1, nodeset1', &
        ' or object must be present in the heading of define_connection'
      stop
    end if

    if ( present(object) ) then
      if ( present(point2) ) then
        write(*,'(/3(a/))') &
          'Error: point2 ', &
          ' must not be present in the heading of define_connection', &
          ' together with object in the heading '
        stop
      end if
    else
      if ( count( [ present(point2), present(curve2), present(surface2), &
                    present(volume2), present(elementset2), &
                    present(nodeset2) ] ) /= 1 ) then
        write(*,'(/2(a/))') &
          'Error: one of point2, curve2, surface2, elementset2, nodeset2', &
          ' must be present in the heading of define_connection'
        stop
      end if
    end if

!   determine type of call

    if ( present(point1) .or. present(nodeset1) ) then

!     connection in a point or on a nodeset

      input_probdef%numconnections = input_probdef%numconnections + 1
      cg = input_probdef%numconnections

      if ( cg > size(input_probdef%connections) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_connection:', &
          ' connection > maximum = ', size(input_probdef%connections)
        stop
      end if

      if ( present(num) ) num = cg

      if ( present(point1) ) then

!       point

        if ( .not. present(point2) ) then
          write(*,'(/a/a/)') &
            'Error: point1 in heading of define_connection requires point2', &
            '  to be present as well'
          stop
        end if

        if ( point1 <=0 .or. point1 > mesh%npoints ) then
          write(*,'(/a/a,i0/)') &
            'Error: point1 in heading of define_connection has wrong value:', &
            'point1 <=0 or point1 > number of points = ', mesh%npoints
          stop
        end if

        if ( point2 == point1 .or. point2 <=0 &
               .or. point2 > mesh%npoints ) then
          write(*,'(/a/a,i0/)') &
            'Error: point2 in heading of define_connection has wrong value:', &
            'point2 == point1 or point2 <=0 or point2 > number of points = ', &
            mesh%npoints
          stop
        end if

        input_probdef%connections(cg)%typegeometry = 1
        input_probdef%connections(cg)%geometry1 = point1
        input_probdef%connections(cg)%typegeometry2 = 1
        input_probdef%connections(cg)%geometry2 = point2

      else if ( present(nodeset1) ) then

!       nodeset

        if ( .not. present(nodeset2) ) then
          write(*,'(/a/a/)') &
            'Error: nodeset1 in heading of define_connection requires', &
            '  nodeset2 to be present as well'
          stop
        end if

        if ( nodeset1 <=0 .or. nodeset1 > mesh%nnodesets ) then
          write(*,'(/a/a,i0/)') &
            'Error: nodeset1 in heading of define_connection has wrong value:',&
            'nodeset1 <=0 or nodeset1 > number of nodesets = ', mesh%nnodesets
          stop
        end if

        if ( nodeset2 == nodeset1 .or. nodeset2 <=0 &
               .or. nodeset2 > mesh%nnodesets ) then
          write(*,'(/a/2a,i0/)') &
            'Error: nodeset2 in heading of define_connection has wrong value:',&
            'nodeset2 == nodeset1 or nodeset2 <=0 or', &
            ' nodeset2 > number nodesets = ', mesh%nnodesets
          stop
        end if

        input_probdef%connections(cg)%typegeometry = 5
        input_probdef%connections(cg)%geometry1 = nodeset1
        input_probdef%connections(cg)%typegeometry2 = 5
        input_probdef%connections(cg)%geometry2 = nodeset2

      end if

      if ( present(discretization) ) then
        if ( discretization == 'weak' ) then
          write(*,'(2(/a)/)') &
            'Error in the heading of define_connection:', &
            ' discretization=''weak'' invalid for points or nodesets '
          stop
        else if ( discretization == 'collocation' ) then
          input_probdef%connections(cg)%discretization = 1
        else
          write(*,'(2(/a)/)') &
            'Error: discretization in the heading of define_connection', &
            'must have the value ''weak'' or ''collocation'''
          stop
        end if
      else
!       default: collocation connection
        input_probdef%connections(cg)%discretization = 1
      end if

    else if ( present(curve1) .or. present(surface1) .or. &
              present(volume1) .or. present(elementset1) ) then

!     connection on a curve, surface, volume or elementset

      input_probdef%numconnections = input_probdef%numconnections + 1
      cg = input_probdef%numconnections

      if ( cg > size(input_probdef%connections) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_connection:', &
          ' connection > maximum = ', size(input_probdef%connections)
        stop
      end if

      if ( present(num) ) num = cg

      if ( present(curve1) ) then

!       curve

        if ( count( [ present(curve2), present(elementset2), &
                      present(nodeset2) ] ) == 0 ) then
          write(*,'(/a/a/)') &
            'Error: curve1 in heading of define_connection requires', &
            '  either curve2, elementset2 or nodeset2 to be present as well'
          stop
        end if

        if ( curve1 <=0 .or. curve1 > mesh%ncurves ) then
          write(*,'(/a/a,i0/)') &
            'Error: curve1 in heading of define_connection has wrong value:', &
            'curve1 <=0 or curve1 > number curves = ', mesh%ncurves
          stop
        end if

        input_probdef%connections(cg)%typegeometry = 2
        input_probdef%connections(cg)%geometry1 = curve1

        if ( present(curve2) ) then

          if ( curve2 == curve1 .or. curve2 <=0 &
                 .or. curve2 > mesh%ncurves ) then
            write(*,'(/a/2a,i0/)') &
              'Error: curve2 in heading of define_connection has wrong value:',&
              'curve2 == curve1 or curve2 <=0 or', &
              ' curve2 > number curves = ', mesh%ncurves
            stop
          end if

          input_probdef%connections(cg)%typegeometry2 = 2
          input_probdef%connections(cg)%geometry2 = curve2

        else if ( present(elementset2) ) then

          if ( elementset2 <=0 .or. elementset2 > mesh%nelementsets ) then
            write(*,'(/2a/a,i0/)') &
              'Error: elementset2 in heading of define_connection has ', &
              ' wrong value: ', &
              ' elementset1 <=0 or elementset1 > number elementsets = ', &
              mesh%nelementsets
            stop
          end if

          input_probdef%connections(cg)%geometry2 = 0
          input_probdef%connections(cg)%typegeometry2 = 0
          input_probdef%connections(cg)%elementset2 = elementset2

        else if ( present(nodeset2) ) then

          if ( nodeset2 <=0 .or. nodeset2 > mesh%nnodesets ) then
            write(*,'(/2a/a,i0/)') &
              'Error: nodeset2 in heading of define_connection ', &
              ' has wrong value:',&
              ' nodeset2 <=0 or nodeset2 > number nodesets = ', &
              mesh%nnodesets
            stop
          end if

          input_probdef%connections(cg)%typegeometry2 = 5
          input_probdef%connections(cg)%geometry2 = nodeset2

        end if

      else if ( present(surface1) ) then

!       surface

        if ( count( [ present(surface2), present(elementset2), &
                      present(nodeset2) ] ) == 0 ) then
          write(*,'(/a/a/)') &
            'Error: surface1 in heading of define_connection requires', &
            '  either surface2, elementset2 or nodeset2 to be present as well'
          stop
        end if

        if ( surface1 <=0 .or. surface1 > mesh%nsurfaces ) then
          write(*,'(/a/a,i0/)') &
            'Error: surface1 in heading of define_connection has wrong value:',&
            'surface1 <=0 or surface1 > number surfaces = ', mesh%nsurfaces
          stop
        end if

        input_probdef%connections(cg)%typegeometry = 3
        input_probdef%connections(cg)%geometry1 = surface1

        if ( present(surface2) ) then

          if ( surface2 == surface1 .or. surface2 <=0 &
                 .or. surface2 > mesh%nsurfaces ) then
            write(*,'(/2a/2a,i0/)') &
              'Error: surface2 in heading of define_connection has ', &
              'wrong value:', 'surface2 == surface1 or surface2 <=0 or', &
              ' surface2 > number surfaces = ', mesh%nsurfaces
            stop
          end if

          input_probdef%connections(cg)%typegeometry2 = 3
          input_probdef%connections(cg)%geometry2 = surface2

        else if ( present(elementset2) ) then

          if ( elementset2 <=0 .or. elementset2 > mesh%nelementsets ) then
            write(*,'(/2a/a,i0/)') &
              'Error: elementset2 in heading of define_connection has ', &
              ' wrong value: ', &
              ' elementset1 <=0 or elementset1 > number elementsets = ', &
              mesh%nelementsets
            stop
          end if

          input_probdef%connections(cg)%typegeometry2 = 0
          input_probdef%connections(cg)%geometry2 = 0
          input_probdef%connections(cg)%elementset2 = elementset2

        else if ( present(nodeset2) ) then

          if ( nodeset2 <=0 .or. nodeset2 > mesh%nnodesets ) then
            write(*,'(/2a/a,i0/)') &
              'Error: nodeset2 in heading of define_connection ', &
              ' has wrong value:',&
              ' nodeset2 <=0 or nodeset2 > number nodesets = ', &
              mesh%nnodesets
            stop
          end if

          input_probdef%connections(cg)%typegeometry2 = 5
          input_probdef%connections(cg)%geometry2 = nodeset2

        end if

      else if ( present(volume1) ) then

!       volume

        if ( count( [ present(volume2), present(elementset2), &
                      present(nodeset2) ] ) == 0 ) then
          write(*,'(/a/a/)') &
            'Error: volume1 in heading of define_connection requires', &
            '  either volume2, elementset2 or nodeset2 to be present as well'
          stop
        end if

        if ( volume1 <=0 .or. volume1 > mesh%nvolumes ) then
          write(*,'(/a/a,i0/)') &
            'Error: volume1 in heading of define_connection has wrong value:',&
            'volume1 <=0 or volume1 > number volumes = ', mesh%nvolumes
          stop
        end if

        input_probdef%connections(cg)%typegeometry = 4
        input_probdef%connections(cg)%geometry1 = volume1

        if ( present(volume2) ) then

          if ( volume2 == volume1 .or. volume2 <=0 &
                 .or. volume2 > mesh%nvolumes ) then
            write(*,'(/2a/2a,i0/)') &
              'Error: volume2 in heading of define_connection has ', &
              'wrong value:', 'volume2 == volume1 or volume2 <=0 or', &
              ' volume2 > number volumes = ', mesh%nvolumes
            stop
          end if

          input_probdef%connections(cg)%typegeometry2 = 4
          input_probdef%connections(cg)%geometry2 = volume2

        else if ( present(elementset2) ) then

          if ( elementset2 <=0 .or. elementset2 > mesh%nelementsets ) then
            write(*,'(/2a/a,i0/)') &
              'Error: elementset2 in heading of define_connection has ', &
              ' wrong value: ', &
              ' elementset2 <=0 or elementset2 > number elementsets = ', &
              mesh%nelementsets
            stop
          end if

          input_probdef%connections(cg)%typegeometry2 = 0
          input_probdef%connections(cg)%geometry2 = 0
          input_probdef%connections(cg)%elementset2 = elementset2

        else if ( present(nodeset2) ) then

          if ( nodeset2 <=0 .or. nodeset2 > mesh%nnodesets ) then
            write(*,'(/2a/a,i0/)') &
              'Error: nodeset2 in heading of define_connection ', &
              ' has wrong value:',&
              ' nodeset2 <=0 or nodeset2 > number nodesets = ', &
              mesh%nnodesets
            stop
          end if

          input_probdef%connections(cg)%typegeometry2 = 5
          input_probdef%connections(cg)%geometry2 = nodeset2

        end if

      else if ( present(elementset1) ) then

!       elementset

        if ( .not. present(elementset2) ) then
          write(*,'(/a/a/)') &
            'Error: elementset1 in heading of define_connection requires', &
            '  elementset2 must be present as well'
          stop
        end if

        if ( elementset1 <=0 .or. elementset1 > mesh%nelementsets ) then
          write(*,'(/2a/a,i0/)') &
            'Error: elementset1 in heading of define_connection', &
            ' has wrong value:', &
            'elementset1 <=0 or elementset1 > number elementsets = ', &
            mesh%nelementsets
          stop
        end if

        input_probdef%connections(cg)%typegeometry = 0
        input_probdef%connections(cg)%elementset1 = elementset1

        if ( elementset2 == elementset1 .or. elementset2 <=0 &
               .or. elementset2 > mesh%nelementsets ) then
          write(*,'(/3a,i0/)') &
            'Error: elementset2 in heading of define_connection has ', &
            'wrong value:', 'elementset2 == elementset1 or ', &
            'elementset2 <=0 or elementset2 > number elementsets = ', &
            mesh%nelementsets
          stop
        end if

        input_probdef%connections(cg)%typegeometry2 = 0
        input_probdef%connections(cg)%elementset2 = elementset2

      end if

!     discretization for curve1, surface1, volume1 or elementset1

      if ( present(nodeset2) ) then

        if ( present(discretization) ) then
          if ( discretization == 'weak' ) then
            write(*,'(2(/a)/)') &
              'Error in the heading of define_connection:', &
              ' discretization=''weak'' invalid for nodesets '
            stop
          else if ( discretization == 'collocation' ) then
            input_probdef%connections(cg)%discretization = 1
          else
            write(*,'(2(/a)/)') &
              'Error: discretization in the heading of define_connection', &
              'must have the value ''collocation'' for nodesets'
            stop
          end if
        else
!         default: collocated connection
          input_probdef%connections(cg)%discretization = 1
        end if

      else

        if ( present(discretization) ) then
          if ( discretization == 'weak' ) then
            input_probdef%connections(cg)%discretization = 0
          else if ( discretization == 'collocation' ) then
            input_probdef%connections(cg)%discretization = 1
          else
            write(*,'(2(/a)/)') &
              'Error: discretization in the heading of define_connection', &
              'must have the value ''weak'' or ''collocation'''
            stop
          end if
        else
!         default: weak (integral) connection
          input_probdef%connections(cg)%discretization = 0
        end if

      end if

    else if ( present(object) ) then

!     connection on an object

      if ( object <=0 .or. object > mesh%nobjects ) then
        write(*,'(/a/a,i0/)') &
          'Error: object in heading of define_connection has wrong value:', &
        '  object <=0 or object > number of objects = ', mesh%nobjects
        stop
      end if

      input_probdef%numconnections = input_probdef%numconnections + 1
      cg = input_probdef%numconnections

      if ( cg > size(input_probdef%connections) ) then
        write(*,'(/a/a,i0/)') &
          'Error: maximum exceeded in define_connection:', &
          ' connection > maximum = ', size(input_probdef%connections)
        stop
      end if

      if ( present(num) ) num = cg

      input_probdef%connections(cg)%object = object

!     discretization for object

      if ( present(nodeset2) ) then

        if ( present(discretization) ) then
          if ( discretization == 'weak' ) then
            write(*,'(2(/a)/)') &
              'Error in the heading of define_connection:', &
              ' discretization=''weak'' invalid for nodesets '
            stop
          else if ( discretization == 'collocation' ) then
            input_probdef%connections(cg)%discretization = 1
          else
            write(*,'(2(/a)/)') &
              'Error: discretization in the heading of define_connection', &
              'must have the value ''collocation'' for nodesets'
            stop
          end if
        else
!         default: collocated connection
          input_probdef%connections(cg)%discretization = 1
        end if

      else

        if ( present(discretization) ) then
          if ( discretization == 'weak' ) then
            input_probdef%connections(cg)%discretization = 0
          else if ( discretization == 'collocation' ) then
            input_probdef%connections(cg)%discretization = 1
          else
            write(*,'(2(/a)/)') &
              'Error: discretization in the heading of define_connection ', &
              'must have the value ''weak'' or ''collocation'''
            stop
          end if
        else
!         default: weak (integral) connection
          input_probdef%connections(cg)%discretization = 0
        end if

      end if

      if ( input_probdef%connections(cg)%discretization == 0 .and. &
           .not. mesh%objects(object)%intpoints ) then
        write(*,'(/a/a,i0/a/)') &
          'Error define_connection: ', &
          ' a weak connection has been defined on object ', object, &
          ' however it does not contain integration points.'
        stop
      end if

      if ( present(object2) ) then

!       connection object

        if ( mesh%objects(object)%typeofobject == 2 ) then
          write(*,'(/a/a/)') &
            'Error: object in heading of define_connection is two-sided,', &
            'whereas also object2 is present: this is incompatible.'
          stop
        end if

        if ( object2 < 1 .or. object2 > mesh%nobjects ) then
          write(*,'(/a/a,i0/)') &
            'Error: object2 in heading of define_connection has wrong value:', &
            'object2 < 1 or object2 > number of objects = ', mesh%nobjects
          stop
        end if

        input_probdef%connections(cg)%object2 = object2
        input_probdef%connections(cg)%typegeometry2 = 0
        input_probdef%connections(cg)%geometry2 = 0
        input_probdef%connections(cg)%elementset2 = 0

      else if ( present(curve2) ) then

!       connecting curve

        if ( mesh%objects(object)%typeofobject == 2 ) then
          write(*,'(/a/a/)') &
            'Error: object in heading of define_connection is two-sided,', &
            'whereas also curve2 is present: this is incompatible.'
          stop
        end if

        if ( curve2 <=0 .or. curve2 > mesh%ncurves ) then
          write(*,'(/a/a,i0/)') &
            'Error: curve2 in heading of define_connection has wrong value:',&
            ' curve2 <=0 or curve2 > number curves = ', mesh%ncurves
          stop
        end if

        input_probdef%connections(cg)%object2 = 0
        input_probdef%connections(cg)%typegeometry2 = 2
        input_probdef%connections(cg)%geometry2 = curve2
        input_probdef%connections(cg)%elementset2 = 0

      else if ( present(surface2) ) then

!       connecting surface

        if ( mesh%objects(object)%typeofobject == 2 ) then
          write(*,'(/a/a/)') &
            'Error: object in heading of define_connection is two-sided,', &
            'whereas also surface2 is present: this is incompatible.'
          stop
        end if

        if ( surface2 <=0 .or. surface2 > mesh%nsurfaces ) then
          write(*,'(/a/a,i0/)') &
            'Error: surface2 in heading of define_connection has wrong value:',&
            ' surface2 <=0 or surface2 > number surfaces = ', mesh%nsurfaces
          stop
        end if

        input_probdef%connections(cg)%object2 = 0
        input_probdef%connections(cg)%typegeometry2 = 3
        input_probdef%connections(cg)%geometry2 = surface2
        input_probdef%connections(cg)%elementset2 = 0

      else if ( present(volume2) ) then

!       connecting volume

        if ( mesh%objects(object)%typeofobject == 2 ) then
          write(*,'(/a/a/)') &
            'Error: object in heading of define_connection is two-sided,', &
            'whereas also volume2 is present: this is incompatible.'
          stop
        end if

        if ( volume2 <=0 .or. volume2 > mesh%nvolumes ) then
          write(*,'(/a/a,i0/)') &
            'Error: volume2 in heading of define_connection has wrong value:',&
            ' volume2 <=0 or volume2 > number volumes = ', mesh%nvolumes
          stop
        end if

        input_probdef%connections(cg)%object2 = 0
        input_probdef%connections(cg)%typegeometry2 = 4
        input_probdef%connections(cg)%geometry2 = volume2
        input_probdef%connections(cg)%elementset2 = 0

      else if ( present(elementset2) ) then

        if ( mesh%objects(object)%typeofobject == 2 ) then
          write(*,'(/a/a/)') &
            'Error: object in heading of define_connection is two-sided,', &
            'whereas also elementset2 is present: this is incompatible.'
          stop
        end if

        if ( elementset2 <=0 .or. elementset2 > mesh%nelementsets ) then
          write(*,'(/2a/a,i0/)') &
            'Error: elementset2 in heading of define_connection has ', &
            ' wrong value: ', &
            ' elementset2 <=0 or elementset2 > number elementsets = ', &
            mesh%nelementsets
          stop
        end if

        input_probdef%connections(cg)%object2 = 0
        input_probdef%connections(cg)%typegeometry2 = 0
        input_probdef%connections(cg)%geometry2 = 0
        input_probdef%connections(cg)%elementset2 = elementset2

      else if ( present(nodeset2) ) then

        if ( mesh%objects(object)%typeofobject == 2 ) then
          write(*,'(/a/a/)') &
            'Error: object in heading of define_connection is two-sided,', &
            'whereas also nodeset2 is present: this is incompatible.'
          stop
        end if

        if ( nodeset2 <=0 .or. nodeset2 > mesh%nnodesets ) then
          write(*,'(/2a/a,i0/)') &
            'Error: nodeset2 in heading of define_connection has ', &
            ' wrong value: ', &
            ' nodeset2 <=0 or nodeset2 > number nodesets = ', &
            mesh%nnodesets
          stop
        end if

        input_probdef%connections(cg)%object2 = 0
        input_probdef%connections(cg)%typegeometry2 = 5
        input_probdef%connections(cg)%geometry2 = nodeset2
        input_probdef%connections(cg)%elementset2 = 0

      else

        input_probdef%connections(cg)%object2 = 0
        input_probdef%connections(cg)%typegeometry2 = 0
        input_probdef%connections(cg)%geometry2 = 0
        input_probdef%connections(cg)%elementset2 = 0

      end if

    end if

    input_probdef%connections(cg)%physq1 = 0
    input_probdef%connections(cg)%physq2 = 0

    if ( present(physq) ) then
      input_probdef%connections(cg)%physq1 = physq
      input_probdef%connections(cg)%physq2 = physq
    end if

    if ( present(physq1) ) then
      input_probdef%connections(cg)%physq1 = physq1
    end if

    if ( present(physq2) ) then
      input_probdef%connections(cg)%physq2 = physq2
    end if

    input_probdef%connections(cg)%layer1 = 0
    input_probdef%connections(cg)%layer2 = 0

    if ( present(layer) ) then
      input_probdef%connections(cg)%layer1 = layer
      input_probdef%connections(cg)%layer2 = layer
    end if

    if ( present(layer1) ) then
      input_probdef%connections(cg)%layer1 = layer1
    end if

    if ( present(layer2) ) then
      input_probdef%connections(cg)%layer2 = layer2
    end if

  end subroutine define_connection


! Definition of transformations

  subroutine define_transformation ( mesh, input_probdef, physq, layer, &
    point, curve, surface, volume, nodeset, step, exclude, excludepoints, &
    excludecurves, excludesurfaces, normalvector, v2, Amatrix, orthogonal, num )

    use math_defs_m

    type(mesh_t), intent(in) :: mesh

!   structure to store all definitions for the problem definition
    type(input_probdef_t), intent(inout) :: input_probdef

!   if physq is present physq denotes the physical quantity that must be
!   transformed otherwise all nodal degrees of freedom are assumed.
    integer, intent(in), optional :: physq

!   if layer is present layer denotes the layer that must be transformed
!   otherwise all nodal degrees of freedom are assumed
    integer, intent(in), optional :: layer

!   if point is present define a transformation in the point point.
!   Only normalvector=0 is available for a point.
    integer, intent(in), optional :: point

!   if curve is present define a transformation on the curve curve.
    integer, intent(in), optional :: curve

!   if surface is present define a transformation on the surface surface.
    integer, intent(in), optional :: surface

!   if volume is present define a transformation on the volume volume.
!   Only normalvector=0 is available for volumes.
    integer, intent(in), optional :: volume

!   if nodeset is present define a transformation on the nodeset nodeset.
!   Only normalvector=0 is available for nodesets.
    integer, intent(in), optional :: nodeset

!   step specifies which nodes on curve are used for the transformation
!      step=0 all nodes
!      step>0 nodes 1, 1+step, 1+2*step, ...
!      step<0 all except nodes 1, 1-step, 1-2*step, ...
!
!   exclude specifies which nodes on curve are excluded from the transformation
!      exclude=0 no excludes, all nodes are used (default)
!      exclude=1 exclude first point of curve
!      exclude=2 exclude last point of curve
!      exclude=3 exclude first and last point
!
!   NOTE: step, exclude only work for curves where the local node numbering is
!   in a natural sequence along the curve. This might not be the case for
!   curves constructed from several other curves or curves read from external
!   mesh generators.
    integer, intent(in), optional :: step, exclude

!   exclude points from a transformation on a geometry.
!   For example excludepoints=(/1,3/) excludes nodes in points P1 and P3
!   from the transformation. Note that only the nodes that are
!   actually on the geometry are excluded.
    integer, intent(in), dimension(:), optional :: excludepoints

!   exclude curves from a transformation on a geometry.
!   For example excludecurves=(/1,3/) excludes nodes on curves C1 and C3
!   from the transformation. Note that only the nodes that are
!   actually on the geometry are excluded.
    integer, intent(in), dimension(:), optional :: excludecurves

!   exclude surfaces from a transformation on a geometry.
!   For example excludesurfaces=(/1,3/) excludes nodes on surfaces S1 and S3
!   from the transformation. Note that only the nodes that are
!   actually on the geometry are excluded.
    integer, intent(in), dimension(:), optional :: excludesurfaces

!   for a distributed transformation (Amatrix is absent) normalvector determines
!   how the transformation should be build:
!     normalvector=-1 use -n on the geometry as the x'-axis
!     normalvector=0 matrix A is fully supplied by the user.
!     normalvector=1 use n on the geometry as the x'-axis

!   for a global transformation (Amatrix is present) normalvector is not used
!   and should not be present.
!
!   the default is normalvector=0
    integer, intent(in), optional :: normalvector

!   Vector used for normalvector=-1 or 1 in 3D to define the transformation.
!   If u1 is the unit vector in the x'-axis (normal direction), then the vector
!   u2 pointing in the y'-direction becomes u2 = u1 x v2 and the vector u3
!   pointing in the z'-direction u3 = u1 x u2. The vectors u1, u2 and u3 are
!   first normalized to 1: u1=u1/|u1|, u1=u2/|u2|, u3=u3/|u3|. Then the matrix
!   A becomes the rotation matrix A = [ u1 u2 u3 ]
    real(dp), intent(in), dimension(3), optional :: v2

!   if Amatrix is present, the transformation is global and the transformation
!   matrix A is the same for all nodes. How A is set, depends on the shape of
!   Amatrix:
!      (2,1) : The vector u1=Amatrix(:,1) is the vector pointing in the
!              direction of new x'-axis. The vector u2=[-u1(2),u1(1)] is
!              pointing in the direction of the new y'-axis. The vectors u1
!              and u2=[-u1(2),u1(1)] are first normalized to 1: u1=u1/|u1|,
!              u1=u2/|u2|. Then the matrix A becomes the rotation matrix
!                 A = [ u1 u2 ]
!              This is for 2D only. The number of degrees of freedom involved
!              in the transformation per node should be 2.
!      (3,2) : The vector u1=Amatrix(:,1) is the vector pointing in the
!              direction of new x'-axis. The vector u2=u1 x Amatrix(:,2) is
!              pointing in the direction of the new y'-axis and the vector
!              u3=u1 x u2 is pointing in the direction of the new z'-axis.
!              The vectors u1, u2 and u3 are first normalized to 1: u1=u1/|u1|,
!              u1=u2/|u2|, u3=u3/|u3|.
!              Then the matrix A becomes the rotation matrix
!                 A = [ u1 u2 u3 ]
!              This is for 3D only. The number of degrees of freedom involved
!              per node should be 3.
!      (n,n) : Square matrix: A=Amatrix. The value of n should match the number
!              of degrees of freedom involved in the transformation per node.
!
!   if Amatrix is absent the transformation is distributed
    real(dp), intent(in), dimension(:,:), optional :: Amatrix

!   Orthogonality of the matrix Amatrix:
!     orthogonal=.true.: A is orthogonal, meaning A A^T = I or A^-1=A^T.
!     orthogonal=.false: A is non-orthogonal (general but invertable).
!                        Note, that the transformation from standard degrees
!                        to transformed degrees
!                           u' = A^-1 u
!                        is not yet supported for orthogonal=.false
!   The default is orthogonal=.true.
    logical, intent(in), optional :: orthogonal

!   if present: the transformation number defined.
    integer, intent(out), optional :: num

!
!   This routine can be called MAXNUMTRANSFORMATIONS times to set
!   transformations
!

    logical :: lorthogonal
    integer :: cg, n, m
    real(dp) :: lv
    real(dp), allocatable, dimension(:) :: u1, u2, u3


!   do some testing first

    call check ( mesh, 'define_transformation' )

    if ( .not. input_probdef%created ) then
      write(*,'(/2(a/))') &
        'Error in define_transformation:', &
        ' input_probdef has not been created '
      stop
    end if

    if ( present(Amatrix) .and. present(normalvector) ) then
      write(*,'(/2(a/))') &
        'Error: Amatrix and normalvector cannot be both present ', &
        ' in the heading of define_transformation'
        stop
    end if

    if ( present(physq) ) then
      if ( physq < 1 .or. physq > input_probdef%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq in heading of define_transformation has wrong value:', &
          ' physq < 1 or physq > number of physical quantities = ', &
          input_probdef%nphysq
        stop
      end if
    end if

    if ( present(layer) ) then
      if ( layer < 1 .or. layer > input_probdef%numlayers ) then
        write(*,'(/a/a,i0/)') &
          'Error: layer in heading of define_transformation has wrong value:', &
          ' layer < 1 or layer > number of layers = ', &
          input_probdef%numlayers
        stop
      end if
    end if

    if ( present(excludepoints) ) then
      if ( any( excludepoints <=0 ) .or. &
           any( excludepoints > mesh%npoints ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: excludepoints in heading of define_transformation:', &
          ' value < 0 or value > number of points = ', &
          mesh%npoints
        stop
      end if
    end if

    if ( present(excludecurves) ) then
      if ( any( excludecurves <=0 ) .or. &
           any( excludecurves > mesh%ncurves ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: excludecurves in heading of define_transformation:', &
          ' value < 0 or value > number of curves = ', &
          mesh%ncurves
        stop
      end if
    end if

    if ( present(excludesurfaces) ) then
      if ( any( excludesurfaces <=0 ) .or. &
           any( excludesurfaces > mesh%nsurfaces ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: excludesurfaces in heading of define_transformation:', &
          ' value < 0 or value > number of surfaces = ', &
          mesh%nsurfaces
        stop
      end if
    end if

    if ( count( [ present(point), present(curve), present(surface), &
                  present(volume), present(nodeset) ] ) > 1 ) then
      write(*,'(/2(a/))') &
        'Error: only one of point, curve, surface, volume, nodeset', &
        ' can be present in the heading of define_transformation'
      stop
    end if

    lorthogonal = set_optional ( variable=orthogonal, default=.true. )

!   set transformation

    input_probdef%numtransformations = input_probdef%numtransformations + 1
    cg = input_probdef%numtransformations

    if ( cg > size(input_probdef%transformations) ) then
      write(*,'(/a/a,i0/)') &
        'Error: maximum exceeded in define_transformations:', &
        ' transformations > maximum = ', size(input_probdef%transformations)
      stop
    end if

    if ( present(num) ) num = cg

    if ( present(point) ) then

!     point

      if ( point <=0 .or. point > mesh%npoints ) then
        write(*,'(/a/a,i0/)') &
          'Error: point in heading of define_transformation has wrong value:', &
          'point <=0 or point > number of points = ', mesh%npoints
        stop
      end if
      input_probdef%transformations(cg)%typegeometry = 1
      input_probdef%transformations(cg)%geometry = point

    else if ( present(nodeset) ) then

!     nodeset

      if ( nodeset <=0 .or. nodeset > mesh%nnodesets ) then
        write(*,'(/2a/a,i0/)') &
          'Error: nodeset in heading of define_transformation', &
          ' has wrong value:',&
          'nodeset <=0 or nodeset > number of nodesets = ', mesh%nnodesets
        stop
      end if
      input_probdef%transformations(cg)%typegeometry = 5
      input_probdef%transformations(cg)%geometry = nodeset

    else if ( present(curve) ) then

!     curve

      if ( curve <=0 .or. curve > mesh%ncurves ) then
        write(*,'(/a/a,i0/)') &
          'Error: curve in heading of define_transformations wrong value:', &
          'curve <=0 or curve1 > number of curves = ', mesh%ncurves
        stop
      end if
      input_probdef%transformations(cg)%typegeometry = 2
      input_probdef%transformations(cg)%geometry = curve

    else if ( present(surface) ) then

!     surface

      if ( surface <=0 .or. surface > mesh%nsurfaces ) then
        write(*,'(/a/a,i0/)') &
          'Error: surface in heading of define_transformations wrong value:',&
          'surface <=0 or surface > number of surfaces = ', mesh%nsurfaces
        stop
      end if
      input_probdef%transformations(cg)%typegeometry = 3
      input_probdef%transformations(cg)%geometry = surface

    else if ( present(volume) ) then

!     volume

      if ( volume <=0 .or. volume > mesh%nvolumes ) then
        write(*,'(/a/a,i0/)') &
          'Error: volume in heading of define_transformations wrong value:',&
          'volume <=0 or volume > number of volumes = ', mesh%nvolumes
        stop
      end if
      input_probdef%transformations(cg)%typegeometry = 4
      input_probdef%transformations(cg)%geometry = volume

    else

      write(*,'(7(/a)/)') &
        'Error: could not determine geometry from the heading', &
        'of define_transformation. None of the keywords:',      &
        '  point ',                                                       &
        '  curve ',                                                       &
        '  surface ',                                                     &
        '  volume ',                                                      &
        '  nodeset ',                                                     &
        'is present.'
      stop

    end if

!   Amatrix

    if ( present(Amatrix) ) then

!     global transformation

      n = size(Amatrix,1)
      m = size(Amatrix,2)

      if ( all(shape(Amatrix)==[2,1]) .or. all(shape(Amatrix)==[3,2]) ) then
        if ( .not. lorthogonal ) then
          write(*,'(/2a/)') &
            'Error: shape Amatrix = [2,1] or [3,2] and orthogonal', &
            ' specified .false. '
          stop
        end if
      else if ( n /= m ) then
        write(*,'(/a/)') &
          'Error: shape of Amatrix must be [2,1], [3,2] or square'
        stop
      end if
      input_probdef%transformations(cg)%typetransformation = 2
      allocate(input_probdef%transformations(cg)%Amat_global(n,n))

      if ( all(shape(Amatrix)==[2,1]) ) then

!       2D

        allocate ( u1(2), u2(2) )

        lv = sqrt ( dot_product ( Amatrix(:,1), Amatrix(:,1) ) )
        u1 = Amatrix(:,1) / lv
        u2 = [ -u1(2), u1(1) ]

        input_probdef%transformations(cg)%Amat_global(:,1) = u1
        input_probdef%transformations(cg)%Amat_global(:,2) = u2

        deallocate ( u1, u2 )

      else if ( all(shape(Amatrix)==[3,2]) ) then

!       3D

        allocate ( u1(3), u2(3), u3(3) )

        lv = sqrt ( dot_product ( Amatrix(:,1), Amatrix(:,1) ) )
        u1 = Amatrix(:,1) / lv
        u2 = cross_product ( u1, Amatrix(:,2) )
        lv = sqrt ( dot_product ( u2, u2 ) )
        u2 = u2 / lv
        u3 = cross_product ( u1, u2 )

        input_probdef%transformations(cg)%Amat_global(:,1) = u1
        input_probdef%transformations(cg)%Amat_global(:,2) = u2
        input_probdef%transformations(cg)%Amat_global(:,3) = u3

        deallocate ( u1, u2, u3 )

      else

!       square matrix given

        input_probdef%transformations(cg)%Amat_global = Amatrix

      end if

    else

!     distributed transformation

      input_probdef%transformations(cg)%typetransformation = 1
      allocate(input_probdef%transformations(cg)%Amat_global(0,0))

    end if

!   orthogonality of the transformation matrix

    input_probdef%transformations(cg)%orthogonal = lorthogonal

!   normalvector

    if ( present(normalvector) ) then
      if ( normalvector == 0 ) then
        input_probdef%transformations(cg)%normalvector = 0
      else if ( any ( normalvector == [-1,1] ) ) then
        if ( count( [ present(curve), present(surface) ] ) /= 1 ) then
          write(*,'(/2(a/))') &
            'Error: for normalvector=-1 or 1 one of curve, surface', &
            ' must be present in the heading of define_transformation'
          stop
        end if
        input_probdef%transformations(cg)%normalvector = normalvector
      else
        write(*,'(/a/a/)') &
          'Error: normalvector in heading of define_transformation:', &
          ' must be either -1, 0 or 1.'
        stop
      end if
    else
      input_probdef%transformations(cg)%normalvector = 0
    end if

!   v2

    if ( present(v2) ) then
      if ( all ( normalvector /= [-1,1] ) ) then
        write(*,'(/2(a/))') &
          'Error: v2 present and normalvector not -1 or 1 ', &
          ' in the heading of define_transformation'
        stop
      else if ( mesh%ndim /= 3 ) then
        write(*,'(/a/a,i0/)') &
          'Error: v2 present and dimension of the mesh is not equal to 3', &
          mesh%ndim
        stop
      end if
      input_probdef%transformations(cg)%v2 = v2
    else
      input_probdef%transformations(cg)%v2 = [ 0._dp, 0._dp, 1._dp ]
    end if

!   physq

    input_probdef%transformations(cg)%physq = 0
    if ( present(physq) ) then
      input_probdef%transformations(cg)%physq = physq
    end if

!   layer

    input_probdef%transformations(cg)%layer = 0
    if ( present(layer) ) then
      input_probdef%transformations(cg)%layer = layer
    end if

!   exclude

    if ( present(step) .and. present(curve) ) then
      input_probdef%transformations(cg)%step = step
    else
      input_probdef%transformations(cg)%step = 0
    end if

    if ( present(exclude) .and. present(curve) ) then
      input_probdef%transformations(cg)%exclude = exclude
    else
      input_probdef%transformations(cg)%exclude = 0
    end if

    if ( present(excludepoints) ) then
      input_probdef%transformations(cg)%excludepoints = excludepoints
    else
      allocate (input_probdef%transformations(cg)%excludepoints(0) )
    end if

    if ( present(excludecurves) ) then
      input_probdef%transformations(cg)%excludecurves = excludecurves
    else
      allocate (input_probdef%transformations(cg)%excludecurves(0) )
    end if

    if ( present(excludesurfaces) ) then
      input_probdef%transformations(cg)%excludesurfaces = excludesurfaces
    else
      allocate (input_probdef%transformations(cg)%excludesurfaces(0) )
    end if

  end subroutine define_transformation


! Definition of dependencies

  subroutine define_dependency ( mesh, input_probdef, point1, point2, curve1, &
    curve2, surface1, surface2, volume1, volume2, elementset1, elementset2, &
    nodeset1, nodeset2, object2, physq, physq1, physq2, layer, layer1, layer2, &
    step, exclude, excludepoints, excludecurves, excludesurfaces, &
    naddunknowns, typedependency, datalayout, num )

    type(mesh_t), intent(in) :: mesh

!   structure to store all dependency definitions for the problem definition
    type(input_probdef_t), intent(inout) :: input_probdef

!   if both are present a dependency between the points is defined.
!   If only point1 is present, naddunknowns > 0 is required.
!   Only typedependency='nodes' or 'object' is available for points.
    integer, intent(in), optional :: point1, point2

!   If both curve1, curve2 are present a dependency between the curves is
!   defined. If only curve1 is present, naddunknowns > 0 is required.
!   Note, that the degrees of freedom on curve1 are the dependent degrees of
!   freedom and should be defined essential.
!   Note, that the definition of the curves is important: they are
!   connected according to either the node sequence or element sequence
!   (typedependency=elements).
    integer, intent(in), optional :: curve1, curve2

!   If both surface1, surface2 are present a dependency between the surfaces is
!   defined. If only surface1 is present, naddunknowns > 0 is required.
!   Note, that the degrees of freedom on surface1 are the dependent degrees of
!   freedom and should be defined essential.
!   Note, that the definition of the surfaces is important: they are
!   connected according to either the node sequence or element sequence
!   (typedependency=elements).
    integer, intent(in), optional :: surface1, surface2

!   If both volume1, volume2 are present a dependency between the volumes is
!   defined. If only volume1 is present, naddunknowns > 0 is required.
!   Note, that the degrees of freedom on volume1 are the dependent degrees of
!   freedom and should be defined essential.
!   Note, that the definition of the volumes is important: they are
!   connected according to either the node sequence or element sequence
!   (typedependency=elements).
    integer, intent(in), optional :: volume1, volume2

!   If both elementset1, elementset2 are present a dependency between the
!   elementsets is defined.
!   If only elementset1 is present, naddunknowns > 0 is required.
!   Note, that the degrees of freedom on elementset1 are the dependent degrees
!   of freedom and should be defined essential.
!   Note, that the definition of the elementsets is important: they are
!   connected according to either the node sequence or element sequence
!   (typedependency=elements).
    integer, intent(in), optional :: elementset1, elementset2

!   If both nodeset1, nodeset2 are present a dependency between the nodesets is
!   defined. If only nodeset1 is present naddunknowns > 0 is required.
!   Note that the definition of the nodesets is important: they are connected
!   node for node.
!   Only typedependency='nodes' or 'object' is available for nodesets.
    integer, intent(in), optional :: nodeset1, nodeset2

!   If object2 is present a dependency between the geometry1 and object2 is
!   defined. Note that the definition of the object is important: the geometry
!   and the object are connected node for node (collocation).
!   Only typedependency='object' is available for object2.
    integer, intent(in), optional :: object2

!   if physq is present physq denotes the physical quantity that is involved
!   in the dependency. Otherwise all nodal degrees of freedom are considered.
    integer, intent(in), optional :: physq

!   if physq1 and/or physq2 is present the physical quantity that is involved
!   in the dependency can be different for the two parts.
    integer, intent(in), optional :: physq1, physq2

!   if layer is present layer denotes the layer that is involved
!   in the dependency. Otherwise all nodal degrees of freedom are considered.
    integer, intent(in), optional :: layer

!   if layer1 and/or layer2 is present the layer that is involved
!   in the dependency can be different for the two parts.
    integer, intent(in), optional :: layer1, layer2

!   step specifies which nodes on curve are used for the dependency
!      step=0 all nodes
!      step>0 nodes 1, 1+step, 1+2*step, ...
!      step<0 except nodes 1, 1-step, 1-2*step, ...
!
!   exclude specifies which nodes on curve are excluded from the dependency
!      exclude=0 no excludes, all nodes are used (default)
!      exclude=1 exclude first point of curve
!      exclude=2 exclude last point of curve
!      exclude=3 exclude first and last point
!
!   NOTE: step, exclude only work for curves where the local node numbering is
!   in a natural sequence along the curve. This might not be the case for
!   curves constructed from several other curves or curves read from external
!   mesh generators.
!   NOTE: the nodes excluded are with respect to part 1.
    integer, intent(in), optional :: step, exclude

!   exclude points from a dependency on a geometry.
!   For example excludepoints=(/1,3/) excludes nodes in points P1 and P3
!   from the dependency. Note that only the nodes that are
!   actually on the geometry are excluded.
!   NOTE: the nodes excluded are with respect to part 1.
    integer, intent(in), dimension(:), optional :: excludepoints

!   exclude curves from a dependency on a geometry.
!   For example excludecurves=(/1,3/) excludes nodes on curves C1 and C3
!   from the dependency. Note that only the nodes that are
!   actually on the geometry are excluded.
!   NOTE: the nodes excluded are with respect to part 1.
    integer, intent(in), dimension(:), optional :: excludecurves

!   exclude surfaces from a dependency on a geometry.
!   For example excludesurfaces=(/1,3/) excludes nodes on surfaces S1 and S3
!   from the dependency. Note that only the nodes that are
!   actually on the geometry are excluded.
!   NOTE: the nodes excluded are with respect to part 1.
    integer, intent(in), dimension(:), optional :: excludesurfaces

!   if naddunknowns is present there are additional unknowns in the dependency.
!   The number of additional unknowns is naddunknowns. This is typically used
!   if the dependent unknowns still depend on a few parameters. For example
!   when we impose a rigid body rotation but the velocity and rotation are still
!   unknown.
    integer, intent(in), optional :: naddunknowns

!   typedependency gives the type of type of the dependency:
!       typedependency='nodes' the dependency is node for node.
!       typedependency='elements' the dependency is elementwise.
!       typedependency='full' the dependency is between all degrees.
!   the default is typedependency='nodes' except if object2 is present.
!   typedependency is ignored.
    character(len=*), intent(in), optional :: typedependency

!   datalayout gives the type of the data layout of the system matrix:
!       datalayout='fem' the data layout is according to the fem structure
!       datalayout='general' the data layout is general
!   the default is datalayout='fem'. This is only relevant for
!   typedependency='nodes', since it is reset to datalayout='general' for the
!   other typedependency settings.
    character(len=*), intent(in), optional :: datalayout

!   if present: the dependency number defined.
    integer, intent(out), optional :: num

!   WARNING:
!   For datalayout='fem' it is assumed that AFTER the "elimination" of the
!   degrees of freedom on part1 a standard FEM structure exist between the
!   remaining degrees of freedom. Memory is reserved for that in the system
!   matrix. It is however NOT allowed that the degrees of freedom in part1 and
!   part2 are already connected by the FEM structure. For example, it is not
!   allowed to make dependencies between degrees in elements that have
!   common nodes, because then there is already a connection because of
!   the standard FEM structure.
!
!   This routine can be called MAXNUMDEPENDENCIES times to set dependencies,
!   however be careful with overlapping dependencies. For example, assume you
!   have three degrees of freedom u1, u2 and u3. Imposing
!        u1 = u2
!        u1 = u3
!   in two separate dependencies will lead to the following dependency:
!        u1 = u2 + u3
!   because the dependency matrix A is built using standard fem assembling.
!   However, the following dependencies
!        u2 = u1
!        u3 = u1
!   or
!        u2 = u1
!        u3 = u2
!   will do what you probably want (remove u2 and u3 and retain u1).

    logical :: yesfem
    logical, dimension(:), allocatable :: la1, la2
    integer :: cg


!   do some testing first

    call check ( mesh, 'define_dependency' )

    if ( .not. input_probdef%created ) then
      write(*,'(/2(a/))') &
        'Error in define_dependency:', &
        ' input_probdef has not been created '
      stop
    end if

    if ( present(physq) .and. ( present(physq1) .or. present(physq2) ) ) then
      write(*,'(/2a/)') &
        'Error: physq and physq1/physq2 cannot be present at the same time', &
        ' in the heading of define_dependency'
      stop
    end if

    if ( present(physq) ) then
      if ( physq < 1 .or. physq > input_probdef%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq in heading of define_dependency has wrong value:', &
          'physq < 1 or physq > number of physical quantities = ', &
          input_probdef%nphysq
        stop
      end if
    end if

    if ( present(physq1) ) then
      if ( physq1 < 1 .or. physq1 > input_probdef%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq1 in heading of define_dependency has wrong value:', &
          'physq1 < 1 or physq1 > number of physical quantities = ', &
          input_probdef%nphysq
        stop
      end if
    end if

    if ( present(physq2) ) then
      if ( physq2 <=0 .or. physq2 > input_probdef%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq2 in heading of define_dependency has wrong value:', &
          'physq2 < 1 or physq2 > number of physical quantities = ', &
          input_probdef%nphysq
        stop
      end if
    end if

    if ( present(layer) .and. ( present(layer1) .or. present(layer2) ) ) then
      write(*,'(/2a/)') &
        'Error: layer and layer1/layer2 cannot be present at the same time', &
        ' in the heading of define_dependency'
      stop
    end if

    if ( present(layer) ) then
      if ( layer < 1 .or. layer > input_probdef%numlayers ) then
        write(*,'(/a/a,i0/)') &
          'Error: layer in heading of define_dependency has wrong value:', &
          'layer < 1 or layer > number of layers = ', &
          input_probdef%numlayers
        stop
      end if
    end if

    if ( present(layer1) ) then
      if ( layer1 < 1 .or. layer1 > input_probdef%numlayers ) then
        write(*,'(/a/a,i0/)') &
          'Error: layer1 in heading of define_dependency has wrong value:', &
          'layer1 < 1 or layer1 > number of layers = ', &
          input_probdef%numlayers
        stop
      end if
    end if

    if ( present(layer2) ) then
      if ( layer2 <=0 .or. layer2 > input_probdef%numlayers ) then
        write(*,'(/a/a,i0/)') &
          'Error: layer2 in heading of define_dependency has wrong value:', &
          'layer2 < 1 or layer2 > number of layers = ', &
          input_probdef%numlayers
        stop
      end if
    end if

    if ( present(excludepoints) ) then
      if ( any( excludepoints <=0 ) .or. &
           any( excludepoints > mesh%npoints ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: excludepoints in heading of define_dependency:', &
          ' value < 0 or value > number of points = ', &
          mesh%npoints
        stop
      end if
    end if

    if ( present(excludecurves) ) then
      if ( any( excludecurves <=0 ) .or. &
           any( excludecurves > mesh%ncurves ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: excludecurves in heading of define_dependency:', &
          ' value < 0 or value > number of curves = ', &
          mesh%ncurves
        stop
      end if
    end if

    if ( present(excludesurfaces) ) then
      if ( any( excludesurfaces <=0 ) .or. &
           any( excludesurfaces > mesh%nsurfaces ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: excludesurfaces in heading of define_dependency:', &
          ' value < 0 or value > number of surfaces = ', &
          mesh%nsurfaces
        stop
      end if
    end if

    if ( present(naddunknowns) ) then
      if ( naddunknowns < 0 ) then
        write(*,'(/2a/)') &
          'Error: naddunknowns in heading of define_dependency has a', &
          ' negative value'
        stop
      end if
    end if

    if ( present(typedependency) ) then
      if ( .not. any ( typedependency == &
                    [ 'nodes   ', 'elements', 'full    ' ] ) ) then
        write(*,'(/a/a/)') &
          'Error: typedependency in heading of define_dependency must be ', &
          ' one of ''nodes'', ''elements'' or ''full''.'
        stop
      end if
    end if

    if ( present(datalayout) ) then
      if ( .not. any ( datalayout == [ 'fem    ', 'general' ] ) ) &
      then
        write(*,'(/a/a/)') &
          'Error: datalayout in heading of define_dependency must be ', &
          ' either ''fem'' or ''general''.'
        stop
      end if
    end if

    la1 = [ present(point1), present(curve1), present(surface1), &
            present(volume1), present(elementset1), present(nodeset1) ]

    la2 =  [ present(point2), present(curve2), present(surface2), &
             present(volume2), present(elementset2), present(nodeset2) ]

    if ( count(la1) /= 1 ) then
      write(*,'(/2(a/))') &
        'Error: one of point1, curve1, surface1, volume1, elementset1 or ', &
        ' nodeset1 must be present in the heading of define_dependency'
      stop
    end if

    if ( count(la2) > 1 ) then
      write(*,'(/2(a/))') &
        'Error: at most one of point2, curve2, surface2, volume2, nodeset2, ', &
        ' or elementset2 must be present in the heading of define_dependency'
      stop
    else if ( count(la2) == 1 .and. count ( la1 .and. la2 ) == 0 ) then
      write(*,'(/3(a/))') &
        'Error: second geometry (point2, curve2, surface2, volume2, nodeset2,',&
        ' or elementset2) is of a different kind than the first geometry ', &
        ' in the heading of define_dependency'
      stop
    else if ( present(object2) ) then
      if ( count(la2) /= 0 ) then
        write(*,'(/2(a/))') &
          'Error define_dependency: point2, curve2, surface2, volume2,', &
          '  nodeset2 or elementset2, cannot be combined with object2'
        stop
      else if ( present(elementset1) ) then
        write(*,'(/2(a/))') &
          'Error define_dependency:', &
          '  elementset1 cannot be combined with object2'
        stop
      end if
      if ( present(typedependency) ) then
        if ( typedependency /= 'object' ) then
          write(*,'(2(/a)/)') &
            'Error in the heading of define_dependency:', &
            ' typedependency /= ''object'' and object2 present'
          stop
        end if
      end if
    else if ( count(la2) == 0 .and. present(typedependency) ) then
      if ( typedependency /= 'nodes' ) then
        write(*,'(2(/a)/)') &
          'Error in the heading of define_dependency:', &
          ' typedependency /= ''nodes'' and only the first part is defined'
        stop
      end if
    end if

!   new dependency

    input_probdef%numdependencies = input_probdef%numdependencies + 1
    cg = input_probdef%numdependencies

    if ( cg > size(input_probdef%dependencies) ) then
      write(*,'(/a/a,i0/)') &
        'Error: maximum exceeded in define_dependency:', &
        ' dependency > maximum = ', size(input_probdef%dependencies)
      stop
    end if

    if ( present(num) ) num = cg

!   set data layout

    yesfem = .false.

    if ( present(datalayout) ) then
      select case (datalayout)
      case('fem')
        yesfem = .true.
        input_probdef%dependencies(cg)%datalayout = 0
      case('general')
        input_probdef%dependencies(cg)%datalayout = 1
      case default
        call errormsg_case_default ( 'define_dependency', 'datalayout', &
          char_value=datalayout )
      end select
    else
!     default setting, overruled by typedependency /= nodes
      input_probdef%dependencies(cg)%datalayout = 0
    end if

!   set type of dependency

    if ( present(typedependency) ) then
      select case (typedependency)
      case('nodes')
        input_probdef%dependencies(cg)%typedependency = 0
      case('elements')
        input_probdef%dependencies(cg)%typedependency = 1
      case('full')
        input_probdef%dependencies(cg)%typedependency = 2
      case default
        call errormsg_case_default ( 'define_dependency', 'typedependency', &
          char_value=typedependency )
      end select
    else
!     default setting, overruled by object2
      input_probdef%dependencies(cg)%typedependency = 0
    end if

!   determine type of call

    if ( present(point1) .or. present(nodeset1) ) then

!     dependency in a point or on a nodeset

      if ( present(typedependency) ) then
        if ( typedependency == 'elements' ) then
          write(*,'(2(/a)/)') &
          'Error in the heading of define_dependency:', &
          ' typedependency = ''elements'' not allowed for points or nodesets'
          stop
        end if
      end if

      if ( present(point1) ) then

!       point

        if ( point1 <=0 .or. point1 > mesh%npoints ) then
          write(*,'(/a/a,i0/)') &
            'Error: point1 in heading of define_dependency has wrong value:', &
            'point1 <=0 or point1 > number of points = ', mesh%npoints
          stop
        end if

        input_probdef%dependencies(cg)%typegeometry1 = 1
        input_probdef%dependencies(cg)%geometry1 = point1

        if ( present(point2) ) then

          if ( point2 == point1 .or. point2 <=0 &
                 .or. point2 > mesh%npoints ) then
            write(*,'(/a/a,i0/)') &
            'Error: point2 in heading of define_dependency has wrong value:', &
            'point2 == point1 or point2 <=0 or point2 > number of points = ', &
              mesh%npoints
            stop
          end if

          input_probdef%dependencies(cg)%typegeometry2 = 1
          input_probdef%dependencies(cg)%geometry2 = point2

        else if ( present(object2) ) then

          if ( object2 <=0 .or. object2 > mesh%nobjects ) then
            write(*,'(/a/a,i0/)') &
              'Error define_dependency: object2 has wrong value:',&
              'object2 <=0 or object2 > number of objects = ', mesh%nobjects
            stop
          end if

          input_probdef%dependencies(cg)%typegeometry2 = 0
          input_probdef%dependencies(cg)%geometry2 = 0
          input_probdef%dependencies(cg)%object2 = object2
          input_probdef%dependencies(cg)%typedependency = 3

        else

          if ( .not. present(naddunknowns) ) then

            write(*,'(/a/a/)') &
              'Error: since point2 is absent in heading of define_dependency:',&
              '  naddunknowns must be present'
            stop

          else if ( naddunknowns == 0 ) then

            write(*,'(/2a/)') &
              'Error: naddunknowns in heading of define_dependency', &
              ' must have a value larger than zero if point2 is absent.'
            stop

          end if

          input_probdef%dependencies(cg)%typegeometry2 = 0
          input_probdef%dependencies(cg)%geometry2 = 0
          input_probdef%dependencies(cg)%typedependency = 0

        end if

      else if ( present(nodeset1) ) then

!       nodeset

        if ( nodeset1 <=0 .or. nodeset1 > mesh%nnodesets ) then
          write(*,'(/a/a,i0/)') &
            'Error: nodeset1 in heading of define_dependency has wrong value:',&
            'nodeset1 <=0 or nodeset1 > number of nodesets = ', mesh%nnodesets
          stop
        end if

        input_probdef%dependencies(cg)%typegeometry1 = 5
        input_probdef%dependencies(cg)%geometry1 = nodeset1

        if ( present(nodeset2) ) then

          if ( nodeset2 == nodeset1 .or. nodeset2 <=0 &
                 .or. nodeset2 > mesh%nnodesets ) then
            write(*,'(/2a/2a,i0/)') &
              'Error: nodeset2 in heading of define_dependency has ', &
              'wrong value:',&
              'nodeset2 == nodeset1 or nodeset2 <=0 or', &
              ' nodeset2 > number nodesets = ', mesh%nnodesets
            stop
          end if

          input_probdef%dependencies(cg)%typegeometry2 = 5
          input_probdef%dependencies(cg)%geometry2 = nodeset2

        else if ( present(object2) ) then

          if ( object2 <=0 .or. object2 > mesh%nobjects ) then
            write(*,'(/a/a,i0/)') &
              'Error define_dependency: object2 has wrong value:',&
              'object2 <=0 or object2 > number of objects = ', mesh%nobjects
            stop
          end if

          input_probdef%dependencies(cg)%typegeometry2 = 0
          input_probdef%dependencies(cg)%geometry2 = 0
          input_probdef%dependencies(cg)%object2 = object2
          input_probdef%dependencies(cg)%typedependency = 3

        else

          if ( .not. present(naddunknowns) ) then

            write(*,'(/a/a/)') &
              'Error: since nodeset2 is absent in heading of', &
              ' define_dependency: naddunknowns must be present'
            stop

          else if ( naddunknowns == 0 ) then

            write(*,'(/2a/)') &
              'Error: naddunknowns in heading of define_dependency ', &
              ' must have a value larger than zero if nodeset2 is absent.'
            stop

          end if

          input_probdef%dependencies(cg)%typegeometry2 = 0
          input_probdef%dependencies(cg)%geometry2 = 0
          input_probdef%dependencies(cg)%typedependency = 0

        end if

      end if

    else if ( present(curve1) .or. present(surface1) .or. &
              present(volume1) .or. present(elementset1) ) then

!     dependency on a curve, surface, volume or elementset

      if ( present(curve1) ) then

!       curve

        if ( curve1 <=0 .or. curve1 > mesh%ncurves ) then
          write(*,'(/a/a,i0/)') &
            'Error: curve1 in heading of define_dependency has wrong value:', &
            'curve1 <=0 or curve1 > number of curves = ', mesh%ncurves
          stop
        end if

        input_probdef%dependencies(cg)%typegeometry1 = 2
        input_probdef%dependencies(cg)%geometry1 = curve1

        if ( present(curve2) ) then

          if ( curve2 == curve1 .or. curve2 <=0 &
                 .or. curve2 > mesh%ncurves ) then
            write(*,'(/a/a,i0/)') &
            'Error: curve2 in heading of define_dependency has wrong value:', &
            'curve2 == curve1 or curve2 <=0 or curve2 > number of curves = ', &
              mesh%ncurves
            stop
          end if

          input_probdef%dependencies(cg)%typegeometry2 = 2
          input_probdef%dependencies(cg)%geometry2 = curve2

        else if ( present(object2) ) then

          if ( object2 <=0 .or. object2 > mesh%nobjects ) then
            write(*,'(/a/a,i0/)') &
              'Error define_dependency: object2 has wrong value:',&
              'object2 <=0 or object2 > number of objects = ', mesh%nobjects
            stop
          end if

          input_probdef%dependencies(cg)%typegeometry2 = 0
          input_probdef%dependencies(cg)%geometry2 = 0
          input_probdef%dependencies(cg)%object2 = object2
          input_probdef%dependencies(cg)%typedependency = 3

        else

          if ( .not. present(naddunknowns) ) then

            write(*,'(/a/a/)') &
              'Error: since curve2 is absent in heading of define_dependency:',&
              '  naddunknowns must be present'
            stop

          else if ( naddunknowns == 0 ) then

            write(*,'(/2a/)') &
              'Error: naddunknowns in heading of define_dependency', &
              ' must have a value larger than zero if curve2 is absent.'
            stop

          end if

          input_probdef%dependencies(cg)%typegeometry2 = 0
          input_probdef%dependencies(cg)%geometry2 = 0
          input_probdef%dependencies(cg)%typedependency = 0

        end if

      else if ( present(surface1) ) then

!       surface

        if ( surface1 <=0 .or. surface1 > mesh%nsurfaces ) then
          write(*,'(/a/a,i0/)') &
            'Error: surface1 in heading of define_dependency has wrong value:',&
            'surface1 <=0 or surface1 > number of surfaces = ', mesh%nsurfaces
          stop
        end if

        input_probdef%dependencies(cg)%typegeometry1 = 3
        input_probdef%dependencies(cg)%geometry1 = surface1

        if ( present(surface2) ) then

          if ( surface2 == surface1 .or. surface2 <=0 &
                 .or. surface2 > mesh%nsurfaces ) then
            write(*,'(/a/2a,i0/)') &
            'Error: surface2 in heading of define_dependency has wrong value:',&
            'surface2 == surface1 or surface2 <=0 or surface2 > number of ', &
            'surfaces = ', &
              mesh%nsurfaces
            stop
          end if

          input_probdef%dependencies(cg)%typegeometry2 = 3
          input_probdef%dependencies(cg)%geometry2 = surface2

        else if ( present(object2) ) then

          if ( object2 <=0 .or. object2 > mesh%nobjects ) then
            write(*,'(/a/a,i0/)') &
              'Error define_dependency: object2 has wrong value:',&
              'object2 <=0 or object2 > number of objects = ', mesh%nobjects
            stop
          end if

          input_probdef%dependencies(cg)%typegeometry2 = 0
          input_probdef%dependencies(cg)%geometry2 = 0
          input_probdef%dependencies(cg)%object2 = object2
          input_probdef%dependencies(cg)%typedependency = 3

        else

          if ( .not. present(naddunknowns) ) then

            write(*,'(/a/a/)') &
              'Error: since surface2 is absent in heading of', &
              'define_dependency: naddunknowns must be present'
            stop

          else if ( naddunknowns == 0 ) then

            write(*,'(/a/a/)') &
              'Error: naddunknowns in heading of define_dependency', &
              ' must have a value larger than zero if surface2 is absent.'
            stop

          end if

          input_probdef%dependencies(cg)%typegeometry2 = 0
          input_probdef%dependencies(cg)%geometry2 = 0
          input_probdef%dependencies(cg)%typedependency = 0

        end if

      else if ( present(volume1) ) then

!       volume

        if ( volume1 <=0 .or. volume1 > mesh%nvolumes ) then
          write(*,'(/a/a,i0/)') &
            'Error: volume1 in heading of define_dependency has wrong value:',&
            'volume1 <=0 or volume1 > number of volumes = ', mesh%nvolumes
          stop
        end if

        input_probdef%dependencies(cg)%typegeometry1 = 4
        input_probdef%dependencies(cg)%geometry1 = volume1

        if ( present(volume2) ) then

          if ( volume2 == volume1 .or. volume2 <=0 &
                 .or. volume2 > mesh%nvolumes ) then
            write(*,'(/a/2a,i0/)') &
            'Error: volume2 in heading of define_dependency has wrong value:',&
            'volume2 == volume1 or volume2 <=0 or volume2 > number of ', &
            'volumes = ', &
              mesh%nvolumes
            stop
          end if

          input_probdef%dependencies(cg)%typegeometry2 = 4
          input_probdef%dependencies(cg)%geometry2 = volume2

        else if ( present(object2) ) then

          if ( object2 <=0 .or. object2 > mesh%nobjects ) then
            write(*,'(/a/a,i0/)') &
              'Error define_dependency: object2 has wrong value:',&
              'object2 <=0 or object2 > number of objects = ', mesh%nobjects
            stop
          end if

          input_probdef%dependencies(cg)%typegeometry2 = 0
          input_probdef%dependencies(cg)%geometry2 = 0
          input_probdef%dependencies(cg)%object2 = object2
          input_probdef%dependencies(cg)%typedependency = 3

        else

          if ( .not. present(naddunknowns) ) then

            write(*,'(/a/a/)') &
              'Error: since volume2 is absent in heading of', &
              'define_dependency: naddunknowns must be present'
            stop

          else if ( naddunknowns == 0 ) then

            write(*,'(/a/a/)') &
              'Error: naddunknowns in heading of define_dependency ', &
              ' must have a value larger than zero if volume2 is absent.'
            stop

          end if

          input_probdef%dependencies(cg)%typegeometry2 = 0
          input_probdef%dependencies(cg)%geometry2 = 0
          input_probdef%dependencies(cg)%typedependency = 0

        end if

      else if ( present(elementset1) ) then

!       elementset

        if ( elementset1 <=0 .or. elementset1 > mesh%nelementsets ) then
          write(*,'(/2a/a,i0/)') &
            'Error: elementset1 in heading of define_dependency has wrong', &
            ' value:',&
            'elementset1 <=0 or elementset1 > number of elementsets = ', &
            mesh%nelementsets
          stop
        end if

        input_probdef%dependencies(cg)%typegeometry1 = 0
        input_probdef%dependencies(cg)%geometry1 = 0
        input_probdef%dependencies(cg)%elementset1 = elementset1

        if ( present(elementset2) ) then

          if ( elementset2 == elementset1 .or. elementset2 <=0 &
                 .or. elementset2 > mesh%nelementsets ) then
            write(*,'(/2a/2a,i0/)') &
            'Error: elementset2 in heading of define_dependency has wrong', &
            ' value:',&
            'elementset2 == elementset1 or elementset2 <=0 or ', &
            'elementset2 > number of elementsets = ', &
              mesh%nelementsets
            stop
          end if

          input_probdef%dependencies(cg)%typegeometry2 = 0
          input_probdef%dependencies(cg)%geometry2 = 0
          input_probdef%dependencies(cg)%elementset2 = elementset2

        else

          if ( .not. present(naddunknowns) ) then

            write(*,'(/a/a/)') &
              'Error: since elementset2 is absent in heading of', &
              'define_dependency: naddunknowns must be present'
            stop

          else if ( naddunknowns == 0 ) then

            write(*,'(/a/a/)') &
              'Error: naddunknowns in heading of define_dependency', &
              ' must have a value larger than zero if elementset2 is absent.'
            stop

          end if

          input_probdef%dependencies(cg)%typegeometry2 = 0
          input_probdef%dependencies(cg)%geometry2 = 0
          input_probdef%dependencies(cg)%elementset2 = 0
          input_probdef%dependencies(cg)%typedependency = 0

        end if

      end if

    end if

!   set data layout
    if ( any ( input_probdef%dependencies(cg)%typedependency == [1,2,3] ) ) then
      if ( yesfem ) write(*,'(/a/)') &
        'Warning define_dependency: datalayout set to ''general'''
      input_probdef%dependencies(cg)%datalayout = 1 ! general datalayout
    end if

    input_probdef%dependencies(cg)%physq1 = 0
    input_probdef%dependencies(cg)%physq2 = 0

    if ( present(physq) ) then
      input_probdef%dependencies(cg)%physq1 = physq
      input_probdef%dependencies(cg)%physq2 = physq
    end if

    if ( present(physq1) ) then
      input_probdef%dependencies(cg)%physq1 = physq1
    end if

    if ( present(physq2) ) then
      input_probdef%dependencies(cg)%physq2 = physq2
    end if

    input_probdef%dependencies(cg)%layer1 = 0
    input_probdef%dependencies(cg)%layer2 = 0

    if ( present(layer) ) then
      input_probdef%dependencies(cg)%layer1 = layer
      input_probdef%dependencies(cg)%layer2 = layer
    end if

    if ( present(layer1) ) then
      input_probdef%dependencies(cg)%layer1 = layer1
    end if

    if ( present(layer2) ) then
      input_probdef%dependencies(cg)%layer2 = layer2
    end if

!   exclude

    if ( present(step) .and. present(curve1) ) then
      input_probdef%dependencies(cg)%step = step
    else
      input_probdef%dependencies(cg)%step = 0
    end if

    if ( present(exclude) .and. present(curve1) ) then
      input_probdef%dependencies(cg)%exclude = exclude
    else
      input_probdef%dependencies(cg)%exclude = 0
    end if

    if ( present(excludepoints) ) then
      input_probdef%dependencies(cg)%excludepoints = excludepoints
    else
      allocate (input_probdef%dependencies(cg)%excludepoints(0) )
    end if

    if ( present(excludecurves) ) then
      input_probdef%dependencies(cg)%excludecurves = excludecurves
    else
      allocate (input_probdef%dependencies(cg)%excludecurves(0) )
    end if

    if ( present(excludesurfaces) ) then
      input_probdef%dependencies(cg)%excludesurfaces = excludesurfaces
    else
      allocate (input_probdef%dependencies(cg)%excludesurfaces(0) )
    end if

    if ( present(naddunknowns) ) then
      input_probdef%dependencies(cg)%naddunknowns = naddunknowns
    else
      input_probdef%dependencies(cg)%naddunknowns = 0
    end if

  end subroutine define_dependency


! Problem definition

  subroutine problem_definition ( input_probdef, mesh, problem )

    type(input_probdef_t), intent(in) :: input_probdef
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(inout) :: problem

    call check ( mesh, 'problem_definition' )

    if ( .not. input_probdef%created ) then
      write(*,'(/2(a/))') &
        'Error in problem_definition:', &
        ' input_probdef has not been created '
      stop
    end if

    if ( problem%created ) then
      write(*,'(/2(a/))') &
        'Error in problem_definition:', &
        ' problem has already been defined '
      stop
    end if

    problem%probnr = input_probdef%probnr
    problem%nelgrp = mesh%nelgrp

!   Part of the problem definition concerning the number of degrees of freedom
!   in the nodal points of the mesh

    call problem_definition_numdegfd ( input_probdef, mesh, problem )

!   Part of the problem definition concerning the number of degrees of freedom
!   of the vectors of special structure

    call problem_definition_vec_numdegfd ( input_probdef, mesh, problem )

!   Part of the problem definition concerning the essential boundary conditions

    call problem_definition_essnodes ( input_probdef, mesh, problem )

!   Part of the problem definition concerning the dependencies

    call problem_definition_dependencies ( input_probdef, mesh, problem )

!   Part of the problem definition concerning the connections

    call problem_definition_constraints ( input_probdef, mesh, problem )

!   Part of the problem definition concerning the connections

    call problem_definition_connections ( input_probdef, mesh, problem )

!   Set number of degrees of freedom

    problem%numdegfd = problem%numnodaldegfd + problem%numdependegfd + &
                       problem%numconstrdegfd

!   Part of the problem definition concerning the renumbering

    call problem_definition_renumber ( mesh, problem )

!   Part of the problem definition concerning the transformations

    call problem_definition_transformations ( input_probdef, mesh, problem )

    problem%created = .true.

  end subroutine problem_definition


! Part of the problem definition concerning the number of degrees of freedom

  subroutine problem_definition_numdegfd ( input_probdef, mesh, problem )

    type(input_probdef_t), intent(in) :: input_probdef
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(inout) :: problem


    integer :: elgrp, node, ndegfd, ndegfd_prev, nodenr, elem, numwar
    integer :: ig, ag, vec, lay, numl
    logical :: agroups(mesh%nelgrp)


!   do some testing first

!   test elementdof

    do elgrp = 1, mesh%nelgrp

      if ( any( input_probdef%elementdof(elgrp)%a < -1 ) ) then
        write (*,'(/a/a,i0,a)') &
          'Error in problem_definition_numdegfd:', &
          ' input_probdef%elementdof for element group ', elgrp, &
          ' contains negative numbers:'
        write (*,'(20(i0,1x))') input_probdef%elementdof(elgrp)%a
        stop
      end if

    end do

!   test vec_elementdof

    do elgrp = 1, mesh%nelgrp
      do vec = 1, input_probdef%nvec

        if ( any( input_probdef%vec_elementdof(elgrp)%a(:,vec) < 0 ) ) then
          write (*,'(/a/a,i0/a,i0,a)') &
            'Error in problem_definition_numdegfd:', &
            ' input_probdef%vec_elementdof for element group ', elgrp, &
            ' contains negative numbers for vector nr ', vec, ' :'
          write (*,'(20(i0,1x))') input_probdef%vec_elementdof(elgrp)%a(:,vec)
          stop
        end if

      end do
    end do

!   test inactive groups

    if ( input_probdef%numinactivegroups > 0 ) then

!     groups>=1  <=nelgrp

      if ( any( input_probdef%inactivegroups < 1 ) .or. &
           any( input_probdef%inactivegroups > mesh%nelgrp ) ) then
          write(*,'(2(/a)/)') &
            'Error in problem_definition_numdegfd: ',&
            ' groups in array input_probdef%inactivegroups(:) are out of range.'
          stop
      end if

    end if

    if ( input_probdef%nphysq > 0 ) then

!     physq>= ?1

      if ( minval(input_probdef%physq) < 1 ) then
          write(*,'(3(/a)/)') &
            'Error in problem_definition_numdegfd: ',&
            ' a vector number in the physical quantities array', &
            ' input_probdef%physq(:) must be at least have a value of 1.'
          stop
      end if

!     test number of vectors available for physical quantities

      if ( input_probdef%nvec < maxval(input_probdef%physq) ) then
        write(*,'(/a/a/2a/)') &
          'Error in problem_definition_numdegfd: the number of vectors of ',&
          ' special structure (nvec) must at least be equal to the maximum',&
          ' vector number in the physical quantities array', &
          ' input_probdef%physq(:)'
        stop
      end if

!     test whether physical quantities fit in elementvector

      do elgrp = 1, mesh%nelgrp

!       skip inactive groups
        if ( any( elgrp == input_probdef%inactivegroups ) ) cycle

        if ( all( input_probdef%elementdof(elgrp)%a /= -1 ) ) then
!         user has specified element degrees of freedom
          do node = 1, mesh%elnumnod(elgrp)
            if ( input_probdef%elementdof(elgrp)%a(node) <            &
                 sum ( input_probdef%vec_elementdof(elgrp)%a(node,    &
                        &input_probdef%physq ) ) &
               ) then
              write (*,'(/4(a/),2(a,i0)/)') &
                'Error in problem_definition_numdegfd: the sum of the number', &
                ' of degrees of freedom in the physical quantities is larger', &
                ' than the number of degrees of freedom in the nodes of the',  &
                ' element.', &
                ' element group ', elgrp, ' element node ', node
              stop
            end if
          end do
        end if

      end do

    else

!     test elementdof

      do elgrp = 1, mesh%nelgrp

        if ( any( input_probdef%elementdof(elgrp)%a == -1 ) ) then
          write (*,'(/a/a,i0,a/a)') &
            'Error in problem_definition_numdegfd:', &
            ' input_probdef%elementdof for element group ', elgrp, &
            ' not filled', &
            ' one or more values are still -1 (initialized value):'
          write (*,'(20(i0,1x))') input_probdef%elementdof(elgrp)%a
          stop
        end if

      end do

    end if

!   test shifted physq

    if ( input_probdef%nphysqshifted > 0 ) then

!     physq>=1  <=nphys

      if ( minval(input_probdef%physqshifted) < 1 .or. &
           maxval(input_probdef%physqshifted) > input_probdef%nphysq ) then
          write(*,'(3(/a)/)') &
            'Error in problem_definition_numdegfd: ',&
            ' physical quantity numbers in array', &
            ' input_probdef%physqshifted(:) are out of range.'
          stop
      end if

    end if

!   test layers

    if ( input_probdef%numlayers > 0 ) then

!     layer>=1  <=nnodesets

      if ( minval(input_probdef%layers) < 1 .or. &
           maxval(input_probdef%layers) > mesh%nnodesets ) then
          write(*,'(3(/a)/)') &
            'Error in problem_definition_numdegfd: ',&
            ' nodesets in array', &
            ' input_probdef%layers(:) are out of range.'
          stop
      end if

    end if

!   physq

    problem%nphysq = input_probdef%nphysq
    problem%physq = input_probdef%physq

!   physq shifted

    problem%nphysqshifted = input_probdef%nphysqshifted
    problem%physqshifted = input_probdef%physqshifted

!   physq mask

    problem%physqmask = input_probdef%physqmask

!   layers

    problem%numlayers = input_probdef%numlayers

    if ( input_probdef%numlayers > 0 ) then
      allocate(problem%nodlayers(mesh%nnodes))
      problem%nodlayers = 0
      do lay = 1, input_probdef%numlayers
        problem%nodlayers(mesh%nodesets(input_probdef%layers(lay))%a) = &
         ibset(problem%nodlayers(mesh%nodesets(input_probdef%layers(lay))%a),&
                lay-1)
      end do
      problem%maxnodnumlayers = &
         maxval ( count_layers ( problem%nodlayers, problem%numlayers ) )
    else
      allocate(problem%nodlayers(0)) ! zero size array
      problem%maxnodnumlayers = 0
    end if

!   active and inactive groups

    agroups = .true.
    agroups(input_probdef%inactivegroups) = .false.
    problem%numinactivegroups = mesh%nelgrp - count ( agroups )
    allocate(problem%inactivegroups(problem%numinactivegroups))
    allocate(problem%activegroups(mesh%nelgrp-problem%numinactivegroups))
    ag =  0
    ig =  0
    do elgrp = 1, mesh%nelgrp
      if ( agroups(elgrp) ) then
!       active group
        ag = ag + 1
        problem%activegroups(ag) = elgrp
      else
!       inactive group
        ig = ig + 1
        problem%inactivegroups(ig) = elgrp
      end if
    end do

!   fill elnumdegfd

    allocate( problem%elnumdegfd(problem%nelgrp) )
    do elgrp = 1, problem%nelgrp
      allocate( problem%elnumdegfd(elgrp)%a(mesh%elnumnod(elgrp)) )
      if ( all( input_probdef%elementdof(elgrp)%a /= -1 ) ) then
!       user has specified element degrees of freedom
        problem%elnumdegfd(elgrp)%a = input_probdef%elementdof(elgrp)%a
      else
!       user has not specified element degrees of freedom: element degrees of
!       freedom are composed of physical quantities only
        problem%elnumdegfd(elgrp)%a = sum(&
          &input_probdef%vec_elementdof(elgrp)%a(:,problem%physq),dim=2)
      end if
      if ( .not. agroups(elgrp) ) then
!       group not active: zero degrees?
        if ( any( problem%elnumdegfd(elgrp)%a /= 0 ) ) then
          write (*,'(/2(a/),a,i0/)') &
            'Error in problem_definition_numdegfd: inactive groups must ', &
            ' have zero degrees of freedom', &
            ' element group ', elgrp
          print *, problem%elnumdegfd(elgrp)%a
          stop
        end if
      end if
    end do

!   fill nodnumdegfd

    allocate ( problem%nodnumdegfd(mesh%nnodes+1) )

    numwar = 0

    problem%nodnumdegfd(1)  = 0
    problem%nodnumdegfd(2:) = -1  ! set to -1 to detect first node

    do elgrp = 1, problem%nelgrp

      do elem = 1, mesh%grpnumel(elgrp)
        do node = 1, mesh%elnumnod(elgrp)

          nodenr = mesh%topology(elgrp)%a(node,elem)
          ndegfd_prev = problem%nodnumdegfd ( nodenr + 1 )

          if ( problem%numlayers > 0 ) then
!           layers have been defined
            numl = count_layers ( problem%nodlayers(nodenr), problem%numlayers )
            ndegfd = problem%elnumdegfd(elgrp)%a(node) * numl
          else
!           no layers
            ndegfd = problem%elnumdegfd(elgrp)%a(node)
          end if

          if ( ndegfd_prev == -1 ) then

!           first time this node
            problem%nodnumdegfd ( nodenr + 1 ) = ndegfd

          else if ( ndegfd /= ndegfd_prev ) then

!           number of degrees of freedom different

            if ( allow_different_numdegfd_in_nodes ) then

!             it is allowed by the user

              numwar = numwar + 1

              if ( numwar <= MAXNUMWARNINGS ) then
                write(*,'(/a,i0/a/)') &
                  'Warning: different number of degrees of freedom at &
                  &nodal point ', nodenr, &
                  'Using the maximum value.'
              end if

              if ( numwar == MAXNUMWARNINGS ) then
                write(*,'(/a/)') 'Further warnings will be suppressed'
              end if

              if ( ndegfd > ndegfd_prev ) then
                problem%nodnumdegfd ( nodenr + 1 ) = ndegfd
              end if

            else

!             not allowed (default)

              write(*,'(/a,i0/)') &
                'Error: different number of degrees of freedom at &
                &nodal point ', nodenr
              stop

            end if

          end if

        end do
      end do

    end do

!   set isolated or inactive nodes to zero degrees of freedom
    where ( problem%nodnumdegfd == -1 )
      problem%nodnumdegfd = 0
    end where

    problem%maxnoddegfd = maxval(problem%nodnumdegfd)

!   accumulate nodnumdegfd

    do nodenr = 1, mesh%nnodes

      problem%nodnumdegfd ( nodenr + 1 ) = &
           problem%nodnumdegfd ( nodenr + 1 ) + problem%nodnumdegfd ( nodenr )

    end do

    problem%numnodaldegfd = problem%nodnumdegfd(mesh%nnodes+1)

  end subroutine problem_definition_numdegfd


! Part of the problem definition concerning the number of degrees of freedom
! of the vectors of special structure

  subroutine problem_definition_vec_numdegfd ( input_probdef, mesh, problem )

    type(input_probdef_t), intent(in) :: input_probdef
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(inout) :: problem

    integer :: elgrp, node, ndegfd, ndegfd_prev, nodenr, elem, vec, numwar, grp
    integer :: w(mesh%nelgrp), numl


!   test inactive groups

    if ( input_probdef%numinactivegroups > 0 ) then

!     zero degrees?

      do grp = 1, input_probdef%numinactivegroups
        elgrp = input_probdef%inactivegroups(grp)
        if ( any( input_probdef%vec_elementdof(elgrp)%a /= 0 ) ) then
          write (*,'(/2(a/),a,i0/)') &
            'Error in problem_definition_vec_numdegfd: inactive groups must ', &
            ' have zero degrees of freedom in the vectors', &
            ' element group ', elgrp
          stop
        end if
      end do

    end if

    problem%nvec = input_probdef%nvec

!   fill vec_elnumdegfd

    problem%vec_elnumdegfd = input_probdef%vec_elementdof

!   fill vec_nodnumdegfd

    allocate ( problem%vec_numdegfd(problem%nvec,2) )
    allocate ( problem%vec_nodnumdegfd(mesh%nnodes+1,problem%nvec) )

    problem%maxvecnoddegfd = 0
    problem%vec_nodnumdegfd(1,:)  = 0
    problem%vec_nodnumdegfd(2:,:) = -1  ! set to -1 to detect first node

    do vec = 1, problem%nvec

      numwar = 0

      do elgrp = 1, problem%nelgrp

!       skip inactive groups
        if ( any( elgrp == problem%inactivegroups ) ) cycle

        do elem = 1, mesh%grpnumel(elgrp)
          do node = 1, mesh%elnumnod(elgrp)

            nodenr = mesh%topology(elgrp)%a(node,elem)
            ndegfd_prev = problem%vec_nodnumdegfd ( nodenr + 1, vec )

            if ( problem%numlayers > 0 ) then
!             layers have been defined
              numl = &
                  count_layers ( problem%nodlayers(nodenr), problem%numlayers )
              ndegfd = problem%vec_elnumdegfd(elgrp)%a(node,vec) * numl
            else
!             no layers
              ndegfd = problem%vec_elnumdegfd(elgrp)%a(node,vec)
            end if

            if ( ndegfd_prev == -1 ) then

!             first time this node
              problem%vec_nodnumdegfd ( nodenr + 1, vec ) = ndegfd

            else if ( ndegfd /= ndegfd_prev ) then

!             number of degrees of freedom different

              if ( allow_different_numdegfd_in_nodes ) then

!               it is allowed by the user

                numwar = numwar + 1

                if ( numwar <= MAXNUMWARNINGS ) then
                  write(*,'(/a,i0,a,i0/a/)') &
                    'Warning: different number of degrees of freedom at &
                    &nodal point ', nodenr, ' vector of special structure ', &
                    vec, 'Using the maximum value.'
                end if

                if ( numwar == MAXNUMWARNINGS ) then
                  write(*,'(/a/)') 'Further warnings will be suppressed'
                end if

                if ( ndegfd > ndegfd_prev ) then
                  problem%vec_nodnumdegfd ( nodenr + 1, vec ) = ndegfd
                end if

              else

!               not allowed (default)

                write(*,'(/a,i0,a,i0/)') &
                  'Error: different number of degrees of freedom at &
                  &nodal point ', nodenr, ' vector of special structure ', &
                  vec
                stop

              end if

            end if

          end do
        end do

      end do

      where ( problem%vec_nodnumdegfd(:,vec) == -1 )
!       set isolated nodes to zero degrees of freedom
        problem%vec_nodnumdegfd(:,vec) = 0
      end where

      problem%maxvecnoddegfd = max( problem%maxvecnoddegfd,&
                                    maxval(problem%vec_nodnumdegfd(:,vec)) )

!     accumulate vec_nodnumdegfd

      do nodenr = 1, mesh%nnodes

        problem%vec_nodnumdegfd ( nodenr + 1, vec ) = &
                  problem%vec_nodnumdegfd ( nodenr + 1, vec )  &
                     + problem%vec_nodnumdegfd ( nodenr, vec )

      end do

      problem%vec_numdegfd(vec,1) = problem%vec_nodnumdegfd(mesh%nnodes+1,vec)

      do elgrp = 1, problem%nelgrp
        w(elgrp) = mesh%grpnumel(elgrp) * &
                     sum( problem%vec_elnumdegfd(elgrp)%a(:,vec) )
      end do
      problem%vec_numdegfd(vec,2) = sum(w)

    end do

  end subroutine problem_definition_vec_numdegfd


! Part of the problem definition concerning the essential boundary conditions

  subroutine problem_definition_essnodes ( input_probdef, mesh, problem )

    use misc_m, only: sort

    type(input_probdef_t), intent(in) :: input_probdef
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(inout) :: problem

    integer :: worklength, curve, node, elemnr
    integer :: esspnt, esscrv, esssrf, essgrp, essnst
    integer :: fcrv, lcrv, fsrf, lsrf, fgrp, lgrp, fnst, lnst
    integer :: nn, nn0, sz, step, exclude, numnod, nod, surface
    integer :: numec, numes, crvnum, srfnum, elgrp
    integer :: nodeset, i
    integer :: pos(problem%maxnoddegfd)
    integer, allocatable, dimension(:) :: nodes, bitpn, order, work

!   determine worklength for storing nodes for essential bc

    worklength = 0

    if ( input_probdef%numesspoints > 0 ) then

      do esspnt = 1, input_probdef%numesspoints
        if ( input_probdef%esspoints(esspnt)%point > 0 ) then
          worklength = worklength + 1
        end if
        worklength = worklength + size(input_probdef%esspoints(esspnt)%points)
      end do

    end if

    if ( input_probdef%numesscurves > 0 ) then

      do esscrv = 1, input_probdef%numesscurves

         fcrv = input_probdef%esscurves(esscrv)%first
         lcrv = input_probdef%esscurves(esscrv)%last

         do curve = fcrv, lcrv
           worklength = worklength + size(mesh%curves(curve)%nodes)
         end do

         do i = 1, size(input_probdef%esscurves(esscrv)%curves)
           curve = input_probdef%esscurves(esscrv)%curves(i)
           worklength = worklength + size(mesh%curves(curve)%nodes)
         end do

      end do

    end if

    if ( input_probdef%numesssurfaces > 0 ) then

      do esssrf = 1, input_probdef%numesssurfaces

         fsrf = input_probdef%esssurfaces(esssrf)%first
         lsrf = input_probdef%esssurfaces(esssrf)%last

         do surface = fsrf, lsrf
           worklength = worklength + size(mesh%surfaces(surface)%nodes)
         end do

         do i = 1, size(input_probdef%esssurfaces(esssrf)%surfaces)
           surface = input_probdef%esssurfaces(esssrf)%surfaces(i)
           worklength = worklength + size(mesh%surfaces(surface)%nodes)
         end do

      end do

    end if

    if ( input_probdef%numesselements > 0 ) then

      worklength = worklength + input_probdef%numesselements

    end if

    if ( input_probdef%numessgroups > 0 ) then

      do essgrp = 1, input_probdef%numessgroups

         fgrp = input_probdef%essgroups(essgrp)%first
         lgrp = input_probdef%essgroups(essgrp)%last

         do elgrp = fgrp, lgrp
           worklength = worklength + mesh%grpnumel(elgrp) * mesh%elnumnod(elgrp)
         end do

         do i = 1, size(input_probdef%essgroups(essgrp)%elgroups)
           elgrp = input_probdef%essgroups(essgrp)%elgroups(i)
           worklength = worklength + mesh%grpnumel(elgrp) * mesh%elnumnod(elgrp)
         end do

      end do

    end if

    if ( input_probdef%numessnodesets > 0 ) then

      do essnst = 1, input_probdef%numessnodesets

         fnst = input_probdef%essnodesets(essnst)%first
         lnst = input_probdef%essnodesets(essnst)%last

         do nodeset = fnst, lnst
           worklength = worklength + size(mesh%nodesets(nodeset)%a)
         end do

         do i = 1, size(input_probdef%essnodesets(essnst)%nodesets)
           nodeset = input_probdef%essnodesets(essnst)%nodesets(i)
           worklength = worklength + size(mesh%nodesets(nodeset)%a)
         end do

      end do

    end if

!   fill essnodes

    allocate ( nodes(worklength), bitpn(worklength), order(worklength) )

    nn = 0

!   points

    if ( input_probdef%numesspoints > 0 ) then

      do esspnt = 1, input_probdef%numesspoints

        if ( input_probdef%esspoints(esspnt)%point > 0 ) then

          nodes(nn+1) = mesh%points(input_probdef%esspoints(esspnt)%point)
          bitpn(nn+1) = bitpattern ( nodes(nn+1),       &
                input_probdef%esspoints(esspnt)%physq,  &
                input_probdef%esspoints(esspnt)%layer,  &
                input_probdef%esspoints(esspnt)%degfd )

          nn = nn + 1

        end if

        do i = 1, size(input_probdef%esspoints(esspnt)%points)

          nodes(nn+1) = mesh%points(input_probdef%esspoints(esspnt)%points(i))
          bitpn(nn+1) = bitpattern ( nodes(nn+1),       &
              input_probdef%esspoints(esspnt)%physq,  &
              input_probdef%esspoints(esspnt)%layer,  &
              input_probdef%esspoints(esspnt)%degfd )

          nn = nn + 1

        end do

      end do

    end if

!   curves

    if ( input_probdef%numesscurves > 0 ) then

!     create workspace for nodes on curves

      allocate ( work(mesh%nnodes) )
      work = 0

      do esscrv = 1, input_probdef%numesscurves

!       set work for excludes

        call set_work ( mesh, &
          input_probdef%esscurves(esscrv)%excludepoints, &
          input_probdef%esscurves(esscrv)%excludecurves, &
          input_probdef%esscurves(esscrv)%excludesurfaces )

        fcrv = input_probdef%esscurves(esscrv)%first
        lcrv = input_probdef%esscurves(esscrv)%last

        nn0 = nn

        do curve = fcrv, lcrv

          sz = size(mesh%curves(curve)%nodes)
          nodes(nn+1:nn+sz) = mesh%curves(curve)%nodes

          step = input_probdef%esscurves(esscrv)%step

          if ( step == 0 ) then

!           all nodes
            do node = nn+1, nn+sz
              if ( work(nodes(node)) == 0 ) then
!               include node
                bitpn(node) = bitpattern ( nodes(node),       &
                      input_probdef%esscurves(esscrv)%physq,  &
                      input_probdef%esscurves(esscrv)%layer,  &
                      input_probdef%esscurves(esscrv)%degfd )
              else
!               node on an excluded point or curve
                bitpn(node) = 0
              end if
            end do

          else if ( step > 0 ) then

!           set bitpattern for nodes given by step
            bitpn(nn+1:nn+sz) = 0
            do node = nn+1, nn+sz, step
              if ( work(nodes(node)) == 0 ) then
!               include node
                bitpn(node) = bitpattern ( nodes(node),       &
                      input_probdef%esscurves(esscrv)%physq,  &
                      input_probdef%esscurves(esscrv)%layer,  &
                      input_probdef%esscurves(esscrv)%degfd )
              else
!               node on an excluded point or curve
                bitpn(node) = 0
              end if
            end do

          else if ( step < 0 ) then

!           exclude points: set bitpattern to zero given by step
            do node = nn+1, nn+sz
              if ( work(nodes(node)) == 0 ) then
!               include node
                bitpn(node) = bitpattern ( nodes(node),       &
                      input_probdef%esscurves(esscrv)%physq,  &
                      input_probdef%esscurves(esscrv)%layer,  &
                      input_probdef%esscurves(esscrv)%degfd )
              else
!               node on an excluded point or curve
                bitpn(node) = 0
              end if
            end do
            bitpn(nn+1:nn+sz:-step) = 0

          end if

          nn = nn + sz

        end do

        do i = 1, size(input_probdef%esscurves(esscrv)%curves)

          curve = input_probdef%esscurves(esscrv)%curves(i)

          sz = size(mesh%curves(curve)%nodes)
          nodes(nn+1:nn+sz) = mesh%curves(curve)%nodes

!         all nodes
          do node = nn+1, nn+sz
            if ( work(nodes(node)) == 0 ) then
!             include node
              bitpn(node) = bitpattern ( nodes(node),       &
                    input_probdef%esscurves(esscrv)%physq,  &
                    input_probdef%esscurves(esscrv)%layer,  &
                    input_probdef%esscurves(esscrv)%degfd )
            else
!             node on an excluded point or curve
              bitpn(node) = 0
            end if
          end do

          nn = nn + sz

        end do

        if ( nn == nn0 ) cycle  ! no nodes for this esscrv

        exclude = input_probdef%esscurves(esscrv)%exclude

        if ( exclude == 1 .or. exclude == 3 ) then

!         exclude first point
          bitpn(nn0+1) = 0

        end if

        if ( exclude == 2 .or. exclude == 3 ) then

!         exclude last point
          bitpn(nn) = 0

        end if

!       reset work to zero

        call unset_work ( mesh, &
          input_probdef%esscurves(esscrv)%excludepoints, &
          input_probdef%esscurves(esscrv)%excludecurves, &
          input_probdef%esscurves(esscrv)%excludesurfaces )

      end do

      deallocate ( work )

    end if

!   surfaces

    if ( input_probdef%numesssurfaces > 0 ) then

!     create workspace for nodes on surfaces

      allocate ( work(mesh%nnodes) )
      work = 0

      do esssrf = 1, input_probdef%numesssurfaces

!       set work for excludes

        call set_work ( mesh, &
          input_probdef%esssurfaces(esssrf)%excludepoints, &
          input_probdef%esssurfaces(esssrf)%excludecurves, &
          input_probdef%esssurfaces(esssrf)%excludesurfaces )

        fsrf = input_probdef%esssurfaces(esssrf)%first
        lsrf = input_probdef%esssurfaces(esssrf)%last

        nn0 = nn

        do surface = fsrf, lsrf

          sz = size(mesh%surfaces(surface)%nodes)
          nodes(nn+1:nn+sz) = mesh%surfaces(surface)%nodes

!         set nodes for this surface
          do node = nn+1, nn+sz
            if ( work(nodes(node)) == 0 ) then
!             include node
              bitpn(node) = bitpattern ( nodes(node),       &
                    input_probdef%esssurfaces(esssrf)%physq,  &
                    input_probdef%esssurfaces(esssrf)%layer,  &
                    input_probdef%esssurfaces(esssrf)%degfd )
            else
!             node on an excluded point, curve, surface
              bitpn(node) = 0
            end if
          end do

          nn = nn + sz

        end do

        do i = 1, size(input_probdef%esssurfaces(esssrf)%surfaces)

          surface = input_probdef%esssurfaces(esssrf)%surfaces(i)

          sz = size(mesh%surfaces(surface)%nodes)
          nodes(nn+1:nn+sz) = mesh%surfaces(surface)%nodes

!         set nodes for this surface
          do node = nn+1, nn+sz
            if ( work(nodes(node)) == 0 ) then
!             include node
              bitpn(node) = bitpattern ( nodes(node),       &
                    input_probdef%esssurfaces(esssrf)%physq,  &
                    input_probdef%esssurfaces(esssrf)%layer,  &
                    input_probdef%esssurfaces(esssrf)%degfd )
            else
!             node on an excluded point, curve, surface
              bitpn(node) = 0
            end if
          end do

          nn = nn + sz

        end do

!       reset work to zero

        call unset_work ( mesh, &
          input_probdef%esssurfaces(esssrf)%excludepoints, &
          input_probdef%esssurfaces(esssrf)%excludecurves, &
          input_probdef%esssurfaces(esssrf)%excludesurfaces )

      end do

      deallocate ( work )

    end if


!   elements

    if ( input_probdef%numesselements > 0 ) then

      do elemnr = 1, input_probdef%numesselements

        nodes(nn+1) = mesh%topology(input_probdef%esselements(elemnr)%elgrp)%&
                      &a(input_probdef%esselements(elemnr)%node,&
                      &  input_probdef%esselements(elemnr)%elem)
        bitpn(nn+1) = bitpattern ( nodes(nn+1),      &
              input_probdef%esselements(elemnr)%physq,  &
              input_probdef%esselements(elemnr)%layer,  &
              input_probdef%esselements(elemnr)%degfd )

        nn = nn + 1

      end do

    end if


!   groups

    if ( input_probdef%numessgroups > 0 ) then

!     create workspace for nodes on groups

      allocate ( work(mesh%nnodes) )
      work = 0

      do essgrp = 1, input_probdef%numessgroups

!       set work for excludes

        call set_work ( mesh, &
          input_probdef%essgroups(essgrp)%excludepoints, &
          input_probdef%essgroups(essgrp)%excludecurves, &
          input_probdef%essgroups(essgrp)%excludesurfaces )

        fgrp = input_probdef%essgroups(essgrp)%first
        lgrp = input_probdef%essgroups(essgrp)%last

        nn0 = nn

        do elgrp = fgrp, lgrp

          sz = mesh%grpnumel(elgrp)*mesh%elnumnod(elgrp)
          nodes(nn+1:nn+sz) = reshape ( mesh%topology(elgrp)%a, [sz] )

!         set nodes for this group
          do node = nn+1, nn+sz
            if ( work(nodes(node)) == 0 ) then
!             include node
              bitpn(node) = bitpattern ( nodes(node),       &
                    input_probdef%essgroups(essgrp)%physq,  &
                    input_probdef%essgroups(essgrp)%layer,  &
                    input_probdef%essgroups(essgrp)%degfd )
            else
!             node on an excluded curve
              bitpn(node) = 0
            end if
          end do

          nn = nn + sz

        end do

        do i = 1, size(input_probdef%essgroups(essgrp)%elgroups)

          elgrp = input_probdef%essgroups(essgrp)%elgroups(i)

          sz = mesh%grpnumel(elgrp)*mesh%elnumnod(elgrp)
          nodes(nn+1:nn+sz) = reshape ( mesh%topology(elgrp)%a, [sz] )

!         set nodes for this group
          do node = nn+1, nn+sz
            if ( work(nodes(node)) == 0 ) then
!             include node
              bitpn(node) = bitpattern ( nodes(node),       &
                    input_probdef%essgroups(essgrp)%physq,  &
                    input_probdef%essgroups(essgrp)%layer,  &
                    input_probdef%essgroups(essgrp)%degfd )
            else
!             node on an excluded curve
              bitpn(node) = 0
            end if
          end do

          nn = nn + sz

        end do

!       reset work to zero

        call unset_work ( mesh, &
          input_probdef%essgroups(essgrp)%excludepoints, &
          input_probdef%essgroups(essgrp)%excludecurves, &
          input_probdef%essgroups(essgrp)%excludesurfaces )

      end do

      deallocate ( work )

    end if


!   nodesets

    if ( input_probdef%numessnodesets > 0 ) then

!     create workspace for nodes on nodesets

      allocate ( work(mesh%nnodes) )
      work = 0

      do essnst = 1, input_probdef%numessnodesets

!       set work for excludes

        call set_work ( mesh, &
          input_probdef%essnodesets(essnst)%excludepoints, &
          input_probdef%essnodesets(essnst)%excludecurves, &
          input_probdef%essnodesets(essnst)%excludesurfaces )

        fnst = input_probdef%essnodesets(essnst)%first
        lnst = input_probdef%essnodesets(essnst)%last

        nn0 = nn

        do nodeset = fnst, lnst

          sz = size(mesh%nodesets(nodeset)%a)
          nodes(nn+1:nn+sz) = mesh%nodesets(nodeset)%a

!         set nodes for this nodeset
          do node = nn+1, nn+sz
            if ( work(nodes(node)) == 0 ) then
!             include node
              bitpn(node) = bitpattern ( nodes(node),       &
                    input_probdef%essnodesets(essnst)%physq,  &
                    input_probdef%essnodesets(essnst)%layer,  &
                    input_probdef%essnodesets(essnst)%degfd )
            else
!             node on an excluded curve/surface
              bitpn(node) = 0
            end if
          end do

          nn = nn + sz

        end do

        do i = 1, size(input_probdef%essnodesets(essnst)%nodesets)

          nodeset = input_probdef%essnodesets(essnst)%nodesets(i)

          sz = size(mesh%nodesets(nodeset)%a)
          nodes(nn+1:nn+sz) = mesh%nodesets(nodeset)%a

!         set nodes for this nodeset
          do node = nn+1, nn+sz
            if ( work(nodes(node)) == 0 ) then
!             include node
              bitpn(node) = bitpattern ( nodes(node),       &
                    input_probdef%essnodesets(essnst)%physq,  &
                    input_probdef%essnodesets(essnst)%layer,  &
                    input_probdef%essnodesets(essnst)%degfd )
            else
!             node on an excluded curve/surface
              bitpn(node) = 0
            end if
          end do

          nn = nn + sz

        end do

!       reset work to zero

        call unset_work ( mesh, &
          input_probdef%essnodesets(essnst)%excludepoints, &
          input_probdef%essnodesets(essnst)%excludecurves, &
          input_probdef%essnodesets(essnst)%excludesurfaces )

      end do

      deallocate ( work )

    end if


!   sort nodes and permute corresponding bitpatterns

    call sort ( nodes, order )
    bitpn = bitpn(order)

!   extract nodes

    numnod = 0

    do node = 1, size(nodes)

      if ( node < size(nodes) ) then
        if ( nodes(node+1) == nodes(node) ) then
!         node present more than once
          bitpn(node+1) = ior(bitpn(node),bitpn(node+1))
          cycle
        end if
      end if

      if ( bitpn(node) /= 0 ) then
!       only count nodes that really have degrees set
        numnod = numnod + 1
      end if

    end do

    problem%numessnodes = numnod

    allocate ( problem%essnodes(numnod,2) )

    nod = 0

    do node = 1, size(nodes)

      if ( node < size(nodes) ) then
        if ( nodes(node+1) == nodes(node) ) cycle ! node present more than once
      end if

      if ( bitpn(node) /= 0 ) then
!       node that sets essential bc has been found
        nod = nod + 1
        problem%essnodes(nod,1) = nodes(node)
        problem%essnodes(nod,2) = bitpn(node)
      end if

    end do

    deallocate ( nodes, bitpn, order )

  contains


!   set work for excludes

    subroutine set_work ( mesh, excludepoints, excludecurves, excludesurfaces )

      type(mesh_t), intent(in) :: mesh
      integer, dimension(:), intent(in) :: excludepoints, excludecurves, &
        excludesurfaces

      integer :: curve, surface

!     set work for points
      work ( mesh%points(excludepoints) ) = 1

!     set work for curves
      numec = size ( excludecurves )
      do curve = 1, numec
        crvnum = excludecurves(curve)
        work ( mesh%curves(crvnum)%nodes ) = 1
      end do

!     set work for surfaces
      numes = size ( excludesurfaces )
      do surface = 1, numes
        srfnum = excludesurfaces(surface)
        work ( mesh%surfaces(srfnum)%nodes ) = 1
      end do

    end subroutine set_work


!   unset work for excludes

    subroutine unset_work ( mesh, excludepoints, excludecurves, &
      excludesurfaces )

      type(mesh_t), intent(in) :: mesh
      integer, dimension(:), intent(in) :: excludepoints, excludecurves, &
        excludesurfaces

      integer :: curve, surface

!     reset work to zero
      work ( mesh%points(excludepoints) ) = 0
      do curve = 1, numec
        crvnum = excludecurves(curve)
        work ( mesh%curves(crvnum)%nodes ) = 0
      end do
      do surface = 1, numes
        srfnum = excludesurfaces(surface)
        work ( mesh%surfaces(srfnum)%nodes ) = 0
      end do

    end subroutine unset_work


!   set bitpattern for a single node

    function bitpattern ( nodenr, physq, layer, degfd )

      integer, intent(in) :: nodenr, physq, layer
      integer, intent(in), dimension(:) :: degfd
      integer :: bitpattern

      integer :: deg, ndegfd

      bitpattern = 0

!     get positions of degrees

      if ( physq > 0 ) then
        call pos_array_local_node ( problem, nodenr, ndegfd, pos, [physq], &
          layer )
      else
        call pos_array_local_node ( problem, nodenr, ndegfd, pos, layer=layer )
      end if

!     limit the position to the word size

      where ( pos(:ndegfd) > bit_size(bitpattern) )
        pos(:ndegfd) = 0
      end where

!     now set the degrees

      if ( size(degfd) > 0 ) then

!       set specified degrees

        do deg = 1, min(size(degfd),ndegfd)
          if ( degfd(deg) /= 0 .and. pos(deg) > 0 ) then
            bitpattern = ibset(bitpattern,pos(deg)-1)
          end if
        end do

      else

!       set all degrees

        do deg = 1, ndegfd
          if ( pos(deg) > 0 ) then
            bitpattern = ibset(bitpattern,pos(deg)-1)
          end if
        end do

      end if

    end function bitpattern

  end subroutine problem_definition_essnodes


! Part of the problem definition concerning the dependencies

  subroutine problem_definition_dependencies ( input_probdef, mesh, problem )

    type(input_probdef_t), intent(in) :: input_probdef
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(inout) :: problem

    integer :: dep, ip, elementset1, elementset2

!   fill dependencies

    problem%numdependencies = input_probdef%numdependencies

    allocate( problem%dependencies(input_probdef%numdependencies) )

    do dep = 1, problem%numdependencies

!     copy complete structure
      problem%dependencies(dep) = input_probdef%dependencies(dep)

      if ( problem%dependencies(dep)%elementset1 > 0 ) then
        elementset1 = problem%dependencies(dep)%elementset1
        if ( .not. mesh%elementsets(elementset1)%nodes_created ) then
          write(*,'(/a/a,i0,a/a/)') &
            'Error problem_definition_dependencies: ', &
            ' elementset ', elementset1, ' does not contain nodes.', &
            ' Add nodes to the elementset with add_to_mesh '
          stop
        end if
        if ( mesh%elementsets(elementset1)%nnodes == 0 ) then
          write(*,'(/a/a,i0,a/a/)') &
            'Error problem_definition_dependencies: ', &
            ' elementset ', elementset1, ' has zero number of nodes.', &
            ' Empty elementset? '
          stop
        end if
      end if

      if ( problem%dependencies(dep)%elementset2 > 0 ) then
        elementset2 = problem%dependencies(dep)%elementset2
        if ( .not. mesh%elementsets(elementset2)%nodes_created ) then
          write(*,'(/a/a,i0,a/a/)') &
            'Error problem_definition_dependencies: ', &
            ' elementset ', elementset2, ' does not contain nodes.', &
            ' Add nodes to the elementset with add_to_mesh '
          stop
        end if
        if ( mesh%elementsets(elementset2)%nnodes == 0 ) then
          write(*,'(/a/a,i0,a/a/)') &
            'Error problem_definition_dependencies: ', &
            ' elementset ', elementset2, ' has zero number of nodes.', &
            ' Empty elementset? '
          stop
        end if
     end if

    end do

!   create new arrays

    ip = problem%numnodaldegfd

    do dep = 1, problem%numdependencies

!     fill addnumdegfd

      problem%dependencies(dep)%addnumdegfd(1) = ip
      ip = ip + problem%dependencies(dep)%naddunknowns
      problem%dependencies(dep)%addnumdegfd(2) = ip

    end do

    problem%numdependegfd = ip - problem%numnodaldegfd

  end subroutine problem_definition_dependencies


! Part of the problem definition concerning the constraints

  subroutine problem_definition_constraints ( input_probdef, mesh, problem )

    type(input_probdef_t), intent(in) :: input_probdef
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(inout) :: problem


    integer :: constr, ip, elem, node, ndegfd, ndegfd_prev, nodenr, nnodes
    integer :: nelem, elnumnod, step, exclude, vec, object, elgrp
    integer :: curve, numec, crvnum, surface, numes, srfnum
    integer :: elementset1, elementset2, nodeset, point, llayer, chlayer
    integer :: nodenr1, lmnumdegfd
    type(geometry_t) :: geometry
    integer, allocatable, dimension(:) :: work

!   fill constraints

    problem%numconstraints = input_probdef%numconstraints

    allocate( problem%constraints(input_probdef%numconstraints) )

    do constr = 1, problem%numconstraints

!     copy complete structure
      problem%constraints(constr) = input_probdef%constraints(constr)

      if ( problem%constraints(constr)%elementset1 > 0 ) then
        elementset1 = problem%constraints(constr)%elementset1
        if ( .not. mesh%elementsets(elementset1)%nodes_created ) then
          write(*,'(/a/a,i0,a/a/)') &
            'Error problem_definition_constraints: ', &
            ' elementset ', elementset1, ' does not contain nodes.', &
            ' Add nodes to the elementset with add_to_mesh '
          stop
        end if
        if ( mesh%elementsets(elementset1)%nnodes == 0 ) then
          write(*,'(/a/a,i0,a/a/)') &
            'Error problem_definition_constraints: ', &
            ' elementset ', elementset1, ' has zero number of nodes.', &
            ' Empty elementset? '
          stop
        end if
      end if

      if ( problem%constraints(constr)%elementset2 > 0 ) then
        elementset2 = problem%constraints(constr)%elementset2
        if ( .not. mesh%elementsets(elementset2)%nodes_created ) then
          write(*,'(/a/a,i0,a/a/)') &
            'Error problem_definition_constraints: ', &
            ' elementset ', elementset2, ' does not contain nodes.', &
            ' Add nodes to the elementset with add_to_mesh '
          stop
        end if
        if ( mesh%elementsets(elementset2)%nnodes == 0 ) then
          write(*,'(/a/a,i0,a/a/)') &
            'Error problem_definition_constraints: ', &
            ' elementset ', elementset2, ' has zero number of nodes.', &
            ' Empty elementset? '
          stop
        end if
     end if

    end do

!   create new arrays

    ip = problem%numnodaldegfd + problem%numdependegfd

    do constr = 1, problem%numconstraints

      lmnumdegfd = input_probdef%constraints(constr)%lmnumdegfd

!     fill globnumdegfd

      problem%constraints(constr)%globnumdegfd(1) = ip
      ip = ip + problem%constraints(constr)%nglobalc
      problem%constraints(constr)%globnumdegfd(2) = ip

!     fill nodnumdegfd (distributed Lagr, mult only)

      if ( problem%constraints(constr)%typeconstraint == 1 ) then

!       fill nodnumdegfd

        if ( problem%constraints(constr)%geometry1 > 0 ) then
!         constraint on geometry
          if ( problem%constraints(constr)%typegeometry == 1 ) then
!           constraint in point; fill fake geometry
            point = problem%constraints(constr)%geometry1
            geometry%nnodes = 1
            geometry%nodes = [mesh%points(point)]
          else if ( problem%constraints(constr)%typegeometry == 2 ) then
!           constraint on curve
            call copy_g ( mesh%curves(problem%constraints(constr)%geometry1), &
              geometry )
          else if ( problem%constraints(constr)%typegeometry == 3 ) then
!           constraint on surface
            call copy_g ( mesh%surfaces(problem%constraints(constr)%geometry1),&
              geometry )
          else if ( problem%constraints(constr)%typegeometry == 4 ) then
!           constraint on volume
            call copy_g ( mesh%volumes(problem%constraints(constr)%geometry1),&
              geometry )
          else if ( problem%constraints(constr)%typegeometry == 5 ) then
!           constraint on nodeset; fill fake geometry
            nodeset = problem%constraints(constr)%geometry1
            geometry%nnodes = size(mesh%nodesets(nodeset)%a)
            geometry%nodes = mesh%nodesets(nodeset)%a
          end if
          nnodes = geometry%nnodes
          nelem  = geometry%nelem
          elnumnod = geometry%elnumnod
        else if ( problem%constraints(constr)%object > 0 ) then
!         constraint on object
          object = problem%constraints(constr)%object
          nnodes = mesh%objects(object)%nnodes
          nelem  = mesh%objects(object)%nelem
          elnumnod = mesh%objects(object)%elnumnod
        end if

        allocate ( problem%constraints(constr)%nodnumdegfd(nnodes+1) )

        problem%constraints(constr)%nodnumdegfd = 0

        if ( problem%constraints(constr)%discretization == 0 ) then

!         weak form of constraint

          do elem = 1, nelem
            do node = 1, elnumnod

              ndegfd = problem%constraints(constr)%elnumdegfd(node)

              if ( problem%constraints(constr)%geometry1 > 0 ) then
!               constraint on geometry (curve or surface)
                nodenr = geometry%topology(node,elem,1)
                if ( lmnumdegfd == 1 ) then
!                 check layers (only layer1 and not layer2!!)
                  llayer = input_probdef%constraints(constr)%layer1
                  if ( llayer > 0 ) then
                    nodenr1 = geometry%topology(node,elem,2)
                    call adjust_for_layer_A ! set ndegfd=0 if no layer present
                  end if
                end if
              else if ( problem%constraints(constr)%object > 0 ) then
!               constraint on object
                nodenr = mesh%objects(object)%topology(node,elem)
              end if

              ndegfd_prev = problem%constraints(constr)%nodnumdegfd(nodenr+1)

              if ( ndegfd_prev == 0 ) then

!               first time this node
                problem%constraints(constr)%nodnumdegfd(nodenr+1) = ndegfd

              else if ( ndegfd /= ndegfd_prev ) then

!               number of degrees of freedom different
                write(*,'(/a/2(a,i0)/)') &
                  'Warning: different number of degrees of freedom ', &
                  ' (constraints) at nodal point ', nodenr, &
                  ' of constraint nr', constr

                if ( ndegfd > ndegfd_prev ) then
                  problem%constraints(constr)%nodnumdegfd(nodenr+1) = ndegfd
                end if

              end if

            end do
          end do

        else if ( problem%constraints(constr)%discretization == 1 ) then

!         collocation of constraint

          if ( problem%constraints(constr)%step > 0 ) then
!           only nodes with increment step
            step = problem%constraints(constr)%step
          else
!           all nodes
            step = 1
          end if

          do node = 1, nnodes, step

            if ( problem%constraints(constr)%geometry1 > 0 ) then

!             constraint on curve, surface or volume

              nodenr = geometry%nodes(node)

              if ( problem%constraints(constr)%nodenumdegfd > 0 ) then
!               user specified number of degrees of freedom
                ndegfd = problem%constraints(constr)%nodenumdegfd
                chlayer = 1
              else if ( problem%constraints(constr)%physq1 > 0 ) then
!               number of degrees of freedom = constraint physq1 quantity
                vec = problem%physq(problem%constraints(constr)%physq1)
                ndegfd = problem%vec_nodnumdegfd( nodenr+1, vec ) &
                             - problem%vec_nodnumdegfd ( nodenr, vec )
                chlayer = 2
              else if ( problem%constraints(constr)%physq2 > 0 ) then
!               number of degrees of freedom = constraint physq2 quantity
                vec = problem%physq(problem%constraints(constr)%physq2)
                ndegfd = problem%vec_nodnumdegfd( nodenr+1, vec ) &
                             - problem%vec_nodnumdegfd ( nodenr, vec )
                chlayer = 3
              else
!               number of degrees of freedom = constraint quantity
                ndegfd = problem%nodnumdegfd( nodenr+1 ) &
                             - problem%nodnumdegfd ( nodenr )
                chlayer = 2
              end if

              if ( lmnumdegfd == 1 ) then

!               check layers and adjust ndegfd heuristically

                select case (chlayer)
                  case(1) ! no degrees if layer not present
                    llayer = input_probdef%constraints(constr)%layer1
                    nodenr1 = nodenr
                    if ( llayer > 0 ) call adjust_for_layer_A
                  case(2) ! single layer degrees; no degrees if layer
                          ! not present; check only layer1
                    llayer = input_probdef%constraints(constr)%layer1
                    if ( llayer > 0 ) call adjust_for_layer_B
                  case(3) ! single layer degrees; no degrees if layer
                          ! not present; check only layer2
                    llayer = input_probdef%constraints(constr)%layer2
                    if ( llayer > 0 ) call adjust_for_layer_B
                  case default
                    call errormsg_case_default ( &
                      'problem_definition_constraints', 'chlayer', &
                      int_value=chlayer )
                end select

              end if

            else if ( problem%constraints(constr)%object > 0 ) then

!             constraint on object
              elgrp = mesh%objects(object)%grpelm(node,1)

              if ( elgrp == 0 ) cycle   ! no degrees connected

!             user specified number of degrees of freedom
              ndegfd = problem%constraints(constr)%nodenumdegfd

            end if

            problem%constraints(constr)%nodnumdegfd(node+1) = ndegfd

          end do

          if ( problem%constraints(constr)%step < 0 ) then
!           exclude nodes with increment -step
            step = - problem%constraints(constr)%step
            do node = 1, nnodes, step
              problem%constraints(constr)%nodnumdegfd(node+1) = 0
            end do
          end if

          exclude = problem%constraints(constr)%exclude
          if ( exclude == 1 .or. exclude == 3 ) then
!           exclude first node
            problem%constraints(constr)%nodnumdegfd(2) = 0
          end if
          if ( exclude == 2 .or. exclude == 3 ) then
!           exclude last node
            problem%constraints(constr)%nodnumdegfd(nnodes+1) = 0
          end if

          if ( size(problem%constraints(constr)%excludepoints) > 0 .or. &
               size(problem%constraints(constr)%excludecurves) > 0 .or. &
               size(problem%constraints(constr)%excludesurfaces) > 0 ) then
!           exclude nodes in points, on curves and/or surfaces
            allocate ( work(mesh%nnodes) )
            work = 0
!           set work
            work ( &
                mesh%points ( problem%constraints(constr)%excludepoints ) ) = 1
            numec = size ( problem%constraints(constr)%excludecurves )
!           set work
            do curve = 1, numec
              crvnum = problem%constraints(constr)%excludecurves(curve)
              work ( mesh%curves(crvnum)%nodes ) = 1
            end do
            numes = size ( problem%constraints(constr)%excludesurfaces )
!           set work
            do surface = 1, numes
              srfnum = problem%constraints(constr)%excludesurfaces(surface)
              work ( mesh%surfaces(srfnum)%nodes ) = 1
            end do
!           delete nodes
            do node = 1, nnodes
              if ( work( geometry%nodes(node) ) == 1 ) then
!               exclude node
                problem%constraints(constr)%nodnumdegfd(node+1) = 0
              end if
            end do
            deallocate ( work )
          end if

        else

          write(*,'(/a/)') &
            'Internal error in problem_definition_constraint'
          stop

        end if

        problem%constraints(constr)%maxnumdegfd = &
          maxval ( problem%constraints(constr)%nodnumdegfd )

!       accumulate nodnumdegfd

        problem%constraints(constr)%nodnumdegfd(1) = ip

        do nodenr = 1, nnodes

          problem%constraints(constr)%nodnumdegfd ( nodenr + 1 ) = &
               problem%constraints(constr)%nodnumdegfd ( nodenr + 1 ) + &
                    problem%constraints(constr)%nodnumdegfd ( nodenr )

        end do

        ip = problem%constraints(constr)%nodnumdegfd ( nnodes + 1 )

      end if

!     fill addnumdegfd

      problem%constraints(constr)%addnumdegfd(1) = ip
      ip = ip + problem%constraints(constr)%naddunknowns
      problem%constraints(constr)%addnumdegfd(2) = ip

!     remove memory from geometry (only necessary for typegeometry= 1 or 5)

      if ( allocated(geometry%nodes) ) deallocate(geometry%nodes)

    end do

    problem%numconstrdegfd = ip - problem%numnodaldegfd - problem%numdependegfd

  contains


!   modify number of degrees for layer (zero degrees if layer not present)
!   NOTE: nodenr1 is the node number

    subroutine adjust_for_layer_A

      integer, dimension(1) :: numl, pnuml

      call find_layer_in_nodes ( problem, llayer, nodes=[nodenr1], &
          numl=numl, pnuml=pnuml )

      if ( numl(1) == 0 ) ndegfd = 0

    end subroutine adjust_for_layer_A


!   modify number of degrees for layer
!   (degrees for one layer only and zero degrees if layer not present)
!   NOTE: nodenr is the node number

    subroutine adjust_for_layer_B

      integer, dimension(1) :: numl, pnuml

      call find_layer_in_nodes ( problem, llayer, nodes=[nodenr], &
          numl=numl, pnuml=pnuml )

      if ( numl(1) > 1 ) then
        ndegfd = ndegfd / numl(1)
      else if ( numl(1) == 0 ) then
        ndegfd = 0
      end if

    end subroutine adjust_for_layer_B

  end subroutine problem_definition_constraints


! Part of the problem definition concerning the connections

  subroutine problem_definition_connections ( input_probdef, mesh, problem )

    type(input_probdef_t), intent(in) :: input_probdef
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(inout) :: problem

    logical :: equalnnodes, equalnelements
    integer :: conn, geometry1, geometry2, elementset1, elementset2, &
               typegeometry, typegeometry2, object, object2


!   fill connections

    problem%numconnections = input_probdef%numconnections

    allocate( problem%connections(input_probdef%numconnections) )

    do conn = 1, problem%numconnections

!     copy structure
      problem%connections(conn) = input_probdef%connections(conn)

!     some basic testing

      geometry1 = problem%connections(conn)%geometry1
      geometry2 = problem%connections(conn)%geometry2
      typegeometry = problem%connections(conn)%typegeometry
      typegeometry2 = problem%connections(conn)%typegeometry2
      elementset1 = problem%connections(conn)%elementset1
      elementset2 = problem%connections(conn)%elementset2
      object = problem%connections(conn)%object
      object2 = problem%connections(conn)%object2

      if ( input_probdef%connections(conn)%discretization == 0 ) then

!       some testing for weak connections

        if ( elementset1 > 0 ) then

!         weak elementset

          if ( mesh%elementsets(elementset1)%nelem /= &
               mesh%elementsets(elementset2)%nelem ) then

            write(*,'(/a/2(a,i0)/2a/)') &
              'Error problem_definition_connections: ', &
              ' elementset ', elementset1, ' and elementset2 ', elementset2, &
              ' have different number of elements. ', &
              ' Not allowed for a weak connection'
            stop

          end if

        else if ( geometry1 > 0 ) then

!         weak geometry

          if ( geometry2 > 0 ) then

!           connected geometry

            equalnelements = .true.

            select case ( typegeometry )
            case(2); equalnelements = mesh%curves(geometry1)%nelem == &
                                      mesh%curves(geometry2)%nelem
            case(3); equalnelements = mesh%surfaces(geometry1)%nelem == &
                                      mesh%surfaces(geometry2)%nelem
            case(4); equalnelements = mesh%volumes(geometry1)%nelem == &
                                      mesh%volumes(geometry2)%nelem
            case default
              call errormsg_case_default ( 'problem_definition_connections', &
                'typegeometry', int_value=typegeometry )
            end select

            if ( .not. equalnelements ) then
                write(*,'(/a/2(a,i0)/2a/)') &
                  'Error problem_definition_connections: ', &
                  ' geometry ', geometry1, ' and geometry ', geometry2, &
                  ' have different number of elements.', &
                  ' Not allowed for a weak connection'
                stop
            end if

          else if ( elementset2 > 0 ) then

!           connected elementset

            equalnelements = .true.

            select case ( typegeometry )
            case(2); equalnelements = mesh%curves(geometry1)%nelem == &
                                      mesh%elementsets(elementset2)%nelem
            case(3); equalnelements = mesh%surfaces(geometry1)%nelem == &
                                      mesh%elementsets(elementset2)%nelem
            case(4); equalnelements = mesh%volumes(geometry1)%nelem == &
                                      mesh%elementsets(elementset2)%nelem
            case default
              call errormsg_case_default ( 'problem_definition_connections', &
                'typegeometry', int_value=typegeometry )
            end select

            if ( .not. equalnelements ) then
                write(*,'(/a/2(a,i0)/2a/)') &
                  'Error problem_definition_connections: ', &
                  ' geometry ', geometry1, ' and elementset ', elementset2, &
                  ' have different number of elements.', &
                  ' Not allowed for a weak connection'
              stop
            end if

          end if

        else if ( object > 0 .and. object2 > 0 ) then

!         weak object (two of them)

          if ( mesh%objects(object)%nelem /= mesh%objects(object2)%nelem .or. &
               mesh%objects(object)%elnumnod /= mesh%objects(object2)%elnumnod &
               .or. mesh%objects(object)%element%globalshape /= &
                    mesh%objects(object2)%element%globalshape ) then
            write(*,'(/a/a,i0/a/a,i0,a/)') &
              'Error problem_definition_connections: ', &
              ' a weak connection has been defined on object ', object, &
              ' however the number of elements or element type in the', &
              ' connecting object ', object2, ' is not the same.'
            stop
          end if
          if ( .not. mesh%objects(object2)%intpoints ) then
            write(*,'(/a/a,i0/a,i0,a/)') &
              'Error problem_definition_connections: ', &
              ' a weak connection has been defined on object ', object, &
              ' however the connecting object ', object2, &
              ' does not contain integration points.'
            stop
          end if
          if ( mesh%objects(object)%ninti /= mesh%objects(object2)%ninti  ) then
            write(*,'(/a/a,i0/a/a,i0,a/)') &
              'Error problem_definition_connections: ', &
              ' a weak connection has been defined on object ', object, &
              ' however the number of integration points in the', &
              ' connecting object ', object2, ' is not the same.'
            stop
          end if

        else if ( object > 0 .and. geometry2 > 0 ) then

!         connected geometry

          equalnelements = .true.

          select case ( typegeometry2 )
          case(2); equalnelements = mesh%objects(object)%nelem == &
                                    mesh%curves(geometry2)%nelem
          case(3); equalnelements = mesh%objects(object)%nelem == &
                                    mesh%surfaces(geometry2)%nelem
          case(4); equalnelements = mesh%objects(object)%nelem == &
                                    mesh%volumes(geometry2)%nelem
          case default
            call errormsg_case_default ( 'problem_definition_connections', &
              'typegeometry', int_value=typegeometry )
          end select

          if ( .not. equalnelements ) then
              write(*,'(/a/2(a,i0)/2a/)') &
                'Error problem_definition_connections: ', &
                ' object ', object, ' and geometry ', geometry2, &
                ' have different number of elements.', &
                ' Not allowed for a weak connection'
            stop
          end if

        else if ( object > 0 .and. elementset2 > 0 ) then

!         connected elementset

          equalnelements = mesh%objects(object)%nelem == &
                                    mesh%elementsets(elementset2)%nelem

          if ( .not. equalnelements ) then
              write(*,'(/a/2(a,i0)/2a/)') &
                'Error problem_definition_connections: ', &
                ' object ', object, ' and elementset ', elementset2, &
                ' have different number of elements.', &
                ' Not allowed for a weak connection'
            stop
          end if

        end if

      else if ( problem%connections(conn)%discretization == 1 ) then

!       collocation: test on equal number of nodes

        if ( elementset1 > 0 ) then

!         collocated elementset

          if ( .not. mesh%elementsets(elementset1)%nodes_created .or. &
               .not. mesh%elementsets(elementset2)%nodes_created ) then
            write(*,'(/a/2(a,i0),a/a/)') &
              'Error problem_definition_connections: ', &
              ' elementsets ', elementset1, 'and/or ', elementset2, &
              ' do not contain nodes.', &
              ' Add nodes to the elementset with add_to_mesh '
            stop
          end if

          if ( mesh%elementsets(elementset1)%nnodes == 0 .or. &
               mesh%elementsets(elementset2)%nnodes == 0 ) then
            write(*,'(/a/2(a,i0),a/a/)') &
              'Error problem_definition_connections: ', &
              ' elementsets ', elementset1, ' and/or', elementset2, &
              ' has zero number of nodes.', &
              ' Empty elementset? '
            stop
          end if

          if ( mesh%elementsets(elementset1)%nnodes /= &
               mesh%elementsets(elementset2)%nnodes ) then

            write(*,'(/a/2(a,i0)/a/)') &
              'Error problem_definition_connections: ', &
              ' elementset ', elementset1, ' and elementset2 ', elementset2, &
              ' have different number of nodes. Not allowed for collocation'
            stop

          end if

        else if ( geometry1 > 0 ) then

!         collocated geometry

          if ( geometry2 > 0 ) then

!           geometry

            equalnnodes = .true.

            select case ( typegeometry )
            case(2); equalnnodes = mesh%curves(geometry1)%nnodes == &
                                   mesh%curves(geometry2)%nnodes
            case(3); equalnnodes = mesh%surfaces(geometry1)%nnodes == &
                                   mesh%surfaces(geometry2)%nnodes
            case(4); equalnnodes = mesh%volumes(geometry1)%nnodes == &
                                   mesh%volumes(geometry2)%nnodes
            case(5); equalnnodes = size(mesh%nodesets(geometry1)%a) == &
                                   size(mesh%nodesets(geometry2)%a)
            case default
              call errormsg_case_default ( 'problem_definition_connections', &
                'typegeometry', int_value=typegeometry )
            end select

            if ( .not. equalnnodes ) then
                write(*,'(/a/2(a,i0)/a/)') &
                  'Error problem_definition_connections: ', &
                  ' geometry1 ', geometry1, ' and geometry2 ', geometry2, &
                  ' have different number of nodes. Not allowed for collocation'
                stop
            end if

          else if ( elementset2 > 0 ) then

!           connected elementset

            if ( .not. mesh%elementsets(elementset2)%nodes_created ) then
              write(*,'(/a/a,i0,a/a/)') &
                'Error problem_definition_connections: ', &
                ' elementset ', elementset2, ' does not contain nodes.', &
                ' Add nodes to the elementset with add_to_mesh '
              stop
            end if
            if ( mesh%elementsets(elementset2)%nnodes == 0 ) then
              write(*,'(/a/a,i0,a/a/)') &
                'Error problem_definition_connections: ', &
                ' elementset ', elementset2, ' has zero number of nodes.', &
                ' Empty elementset? '
              stop
            end if

            equalnnodes = .true.

            select case ( typegeometry )
            case(2); equalnnodes = mesh%curves(geometry1)%nnodes == &
                                   mesh%elementsets(elementset2)%nnodes
            case(3); equalnnodes = mesh%surfaces(geometry1)%nnodes == &
                                   mesh%elementsets(elementset2)%nnodes
            case(4); equalnnodes = mesh%volumes(geometry1)%nnodes == &
                                   mesh%elementsets(elementset2)%nnodes
            case default
              call errormsg_case_default ( 'problem_definition_connections', &
                'typegeometry', int_value=typegeometry )
            end select

            if ( .not. equalnnodes ) then
                write(*,'(/a/2(a,i0)/2a/)') &
                  'Error problem_definition_connections: ', &
                  ' geometry ', geometry1, ' and elementset ', elementset2, &
                  ' have different number of nodes.', &
                  ' Not allowed for a collocated connection'
              stop
            end if

          end if

        else if ( object > 0 .and. object2 > 0 ) then

!         collocated object (two of them)

          if ( mesh%objects(object)%nnodes /= &
               mesh%objects(object2)%nnodes ) then
            write(*,'(/a/2(a,i0)/)') &
              'Error: in heading of problem_definition_connections:', &
              ' number of nodes on object ', object, &
              ' /= number of nodes on object2', object2
            stop
          end if

        else if ( object > 0 .and. geometry2 > 0 ) then

!         collocated object with geometry

          equalnnodes = .true.

          select case ( typegeometry2 )
          case(2); equalnnodes = mesh%objects(object)%nnodes == &
                                 mesh%curves(geometry2)%nnodes
          case(3); equalnnodes = mesh%objects(object)%nnodes == &
                                 mesh%surfaces(geometry2)%nnodes
          case(4); equalnnodes = mesh%objects(object)%nnodes == &
                                 mesh%volumes(geometry2)%nnodes
          case(5); equalnnodes = mesh%objects(object)%nnodes == &
                                 size(mesh%nodesets(geometry2)%a)
          case default
            call errormsg_case_default ( 'problem_definition_connections', &
              'typegeometry2', int_value=typegeometry2 )
          end select

          if ( .not. equalnnodes ) then
              write(*,'(/a/2(a,i0)/a/)') &
                'Error problem_definition_connections: ', &
                ' object ', object, ' and geometry2 ', geometry2, &
                ' have different number of nodes. Not allowed for collocation'
            stop
          end if

        else if ( object > 0 .and. elementset2 > 0 ) then

!         collocated object with elementset

          if ( .not. mesh%elementsets(elementset2)%nodes_created ) then
            write(*,'(/a/a,i0,a/a/)') &
              'Error problem_definition_connections: ', &
              ' elementset ', elementset2, ' does not contain nodes.', &
              ' Add nodes to the elementset with add_to_mesh '
            stop
          end if
          if ( mesh%elementsets(elementset2)%nnodes == 0 ) then
            write(*,'(/a/a,i0,a/a/)') &
              'Error problem_definition_connections: ', &
              ' elementset ', elementset2, ' has zero number of nodes.', &
              ' Empty elementset? '
            stop
          end if

          equalnnodes = mesh%objects(object)%nnodes == &
                                 mesh%elementsets(elementset2)%nnodes

          if ( .not. equalnnodes ) then
            write(*,'(/a/2(a,i0)/2a/)') &
                'Error problem_definition_connections: ', &
                ' object ', object, ' and elementset ', elementset2, &
                ' have different number of nodes.', &
                ' Not allowed for a collocated connection'
            stop
          end if

        end if

      end if

    end do

  end subroutine problem_definition_connections


! Part of the problem definition concerning the renumbering

  subroutine problem_definition_renumber( mesh, problem )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(inout) :: problem

    integer :: node, ndegfd, nodenr, bp, nndof, sp1, sp2, deg
    integer :: nndund, numess, degfd, newnodenr, i, nndund2, physq, physqnr
    integer, allocatable, dimension(:) :: work1, work2, work3
    integer, allocatable, dimension(:,:) :: work


!   renumbering for the matrix partitioning

    allocate ( problem%degfdperm(problem%numdegfd,2) )

!   make a work array to indicate which nodes have ess bc (for speed);

    allocate ( work(mesh%nnodes,2) )

!   a work space work1 to store the degrees of freedom for ess bc

    allocate ( work1( problem%numessnodes * problem%maxnoddegfd ) )

!   a work space work2 to store the unknown degrees of freedom with a tag for
!   shifting (if sign=-1).

    allocate ( work2( problem%numdegfd ) )

!   a work space work3 to store the degrees of freedom in a node with a tag for
!   shifting (if value is 1).

    allocate ( work3( problem%maxnoddegfd ) )

    work = 0

    do node = 1, problem%numessnodes
      work(problem%essnodes(node,1),1) = 1   ! nodes
      work(problem%essnodes(node,1),2) = problem%essnodes(node,2)  ! bit pattern
    end do

!   do the filling

    nndund  = 0
    nndund2 = 0
    numess  = 0

    do newnodenr = 1, mesh%nnodes

      if ( mesh%renumber ) then
        nodenr = mesh%nodperm(newnodenr)
      else
        nodenr = newnodenr
      end if

!     number of degrees of freedom in the node

      ndegfd = problem%nodnumdegfd ( nodenr + 1 )  &
                  - problem%nodnumdegfd ( nodenr )

      work3(1:ndegfd) = 0

      if ( problem%nphysqshifted > 0 ) then

!       degrees renumbering by shifting physq to the end of the vector
!       fill work3 array: 0=degree not shifted, 1=degree shifted

        do physq = 1, problem%nphysqshifted

          physqnr = problem%physqshifted(physq)

          bp = sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physqnr-1)) &
                         - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physqnr-1)) )
          nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physqnr)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physqnr))

          work3(bp+1:bp+nndof) = 1 ! degrees must be shifted

        end do

      end if

      if ( work(nodenr,1) == 1 ) then

!       node with ess bc

        do degfd = 1, min(ndegfd, NUMBEROFBITS)

          if ( btest(work(nodenr,2),degfd-1) ) then

!           prescribed degree of freedom
            numess = numess + 1
            work1(numess) = problem%nodnumdegfd ( nodenr ) + degfd

          else

!           unknown degree of freedom
            nndund = nndund + 1
            work2(nndund) = problem%nodnumdegfd ( nodenr ) + degfd
            if ( work3(degfd) == 1 ) then
!             tag degree for shifting
              work2(nndund) = - work2(nndund)
              nndund2 = nndund2 + 1
            end if

          end if

        end do

!       remaining degrees

        do degfd = NUMBEROFBITS + 1, ndegfd

!         unknown degree of freedom
          nndund = nndund + 1
          work2(nndund) = problem%nodnumdegfd ( nodenr ) + degfd
          if ( work3(degfd) == 1 ) then
!           tag degree for shifting
            work2(nndund) = - work2(nndund)
            nndund2 = nndund2 + 1
          end if

        end do

      else

!       all degrees of freedom are unknowns for this node

        do degfd = 1, ndegfd

          nndund = nndund + 1
          work2(nndund) = problem%nodnumdegfd ( nodenr ) + degfd
          if ( work3(degfd) == 1 ) then
!           tag degree for shifting
            work2(nndund) = - work2(nndund)
            nndund2 = nndund2 + 1
          end if

        end do

      end if

    end do

!   collect unknown degrees of freedom

    if ( problem%nphysqshifted > 0 ) then

!     shift physical degrees of freedom

      sp1 = 0  ! pointer in first part (degfds that stay)
      sp2 = nndund - nndund2  ! pointer in second part (degfds that are shifted)

      do deg = 1, nndund
        if ( work2(deg) > 0 ) then
!         degree of freedom stays here
          sp1 = sp1 + 1
          problem%degfdperm(sp1,1) = work2(deg)
        else
!         shift degree of freedom to the end
          sp2 = sp2 + 1
          problem%degfdperm(sp2,1) = - work2(deg)
        end if
      end do

      if ( sp1 /= nndund - nndund2 ) stop 'error sp1'
      if ( sp2 /= nndund ) stop 'error sp2'

    else

      problem%degfdperm(1:nndund,1) = work2(1:nndund)

    end if

!   the unknown part: unknowns in nodes + dependency + constraint unknowns

    problem%numundegfd = nndund + problem%numdependegfd + &
                                  problem%numconstrdegfd

!   dependency and constraint degrees live at the end of the system vector

    problem%degfdperm(nndund+1:problem%numundegfd,1) = &
                [ (problem%numnodaldegfd+i, i=1, problem%numdependegfd + &
                                                 problem%numconstrdegfd) ]

!   the essential part shifted to the end in the renumbered vector

    problem%numessdegfd = numess

    problem%degfdperm(problem%numdegfd-numess+1:,1) = work1(1:numess)

    deallocate ( work, work1, work2, work3 )

!   the inverse permutation array

    problem%degfdperm(problem%degfdperm(:,1),2) = &
                             [ (i, i=1, problem%numdegfd) ]

  end subroutine problem_definition_renumber


! Part of the problem definition concerning the transformations

  subroutine problem_definition_transformations ( input_probdef, mesh, problem )

    type(input_probdef_t), intent(in) :: input_probdef
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(inout) :: problem


    integer :: geometry, trans, layer, nnodes, node, nodenr, dof, step, &
      exclude, numec, crvnum, numes, srfnum, ndof, physq, n
!   workdeg: for each degree of freedom the number of degrees coupled
!            for transformation
    integer, allocatable, dimension(:) :: workdeg
!   work: indicate nodes that need to be excluded from the transformation
    logical, allocatable, dimension(:) :: work
!   nodes: array with nodes to be considered for transformation
    integer, allocatable, dimension(:) :: nodes, pos

    allocate ( workdeg(problem%numdegfd), work(mesh%nnodes) )

    workdeg = 0
    work = .false.

!   fill transformations

    problem%numtransformations = input_probdef%numtransformations

    allocate( problem%transformations(input_probdef%numtransformations) )

    if ( problem%numtransformations == 0 ) then

!     no transformations

      problem%numtransdegfd = 0
      allocate ( problem%degfdtrans(0) )

!     create dummy Amat to save memory
      call create_dummy_Amat

      return

    end if

    do trans = 1, problem%numtransformations

!     copy complete structure
      problem%transformations(trans) = input_probdef%transformations(trans)

!     inactive transformation?
      if ( problem%transformations(trans)%typetransformation == 0 ) cycle

      physq = problem%transformations(trans)%physq
      layer = problem%transformations(trans)%layer
      geometry = problem%transformations(trans)%geometry

!     exclude nodes

      call exclude_nodes ( .true. )

!     get nodes

      call get_nodes

!     max number of degrees in the nodes

      call ndof_nodes

      allocate ( pos(ndof) )

!     loop over nodes

      do node = 1, nnodes

        nodenr = nodes(node)

        if ( work(nodenr) ) cycle  ! node excluded

!       check layer

        call check_layer

!       get degrees positions

        if ( physq > 0 ) then
!         physical quantity specified
          call pos_array_node ( problem, nodenr, dof, pos, &
            physqarr=[physq], layer=layer )
        else
          call pos_array_node ( problem, nodenr, dof, pos, layer=layer )
        end if

!       check overlap

        if ( any ( workdeg(pos(1:dof)) > 0 ) ) then
          write(*,'(3(/a)/2(a,i0))') &
            'Error in problem_definition_transformations: ', &
            ' Overlap in transformation: degrees in node', &
            ' already being transformed.', &
            ' transformation = ', trans, &
            ' nodenr = ', nodenr
          if ( physq > 0 ) write(*,'(a,i0)') &
            ' physq = ', physq
          if ( layer > 0 ) write(*,'(a,i0/)') &
            ' layer = ', layer
          write(*,*)
          stop
        end if

!       check size of Amat_global

        if ( problem%transformations(trans)%typetransformation == 2 .and. &
             size( problem%transformations(trans)%Amat_global,1) /= dof ) then
          write(*,'(3(/a)/2(a,i0))') &
            'Error in problem_definition_transformations: ', &
            ' Size of global transformation matrix does not match', &
            ' the number of the degrees being transformed in the node.', &
            ' transformation = ', trans, &
            ' nodenr = ', nodenr
          if ( physq > 0 ) write(*,'(a,i0)') &
            ' physq = ', physq
          if ( layer > 0 ) write(*,'(a,i0)') &
            ' layer = ', layer
          write(*,*)
          stop
        end if

!       fill degrees

        workdeg(pos(1:dof)) = dof

      end do

!     set back exclude nodes

      call exclude_nodes ( .false. )

      deallocate ( nodes, pos )

    end do

!   create array degfdtrans and sparse matrix Amat

    problem%numtransdegfd = count ( workdeg /= 0 )
    allocate ( problem%degfdtrans(problem%numtransdegfd) )

    if ( problem%numtransdegfd > 0 ) then
!     create Amat
      call create_Amat
    else
!     create dummy Amat to save memory
      call create_dummy_Amat
    end if

    deallocate ( workdeg, work )

  contains


!  set or unset nodes that need to be excluded from transformations

    subroutine exclude_nodes ( val )

      logical, intent(in) :: val

      integer :: curve, surface, k

!     set work

      if ( problem%transformations(trans)%typegeometry == 2 ) then
!       curve
        nnodes = mesh%curves(geometry)%nnodes
        step = problem%transformations(trans)%step
        if ( step > 1 ) then
!         exclude nodes in between the ones with increment step
          do k = 2, step
            work(mesh%curves(geometry)%nodes(k:nnodes:step)) = val
          end do
        else if ( step < 0 ) then
!         exclude nodes with increment -step
          work(mesh%curves(geometry)%nodes(1:nnodes:-step)) = val
        end if
        exclude = problem%transformations(trans)%exclude
        if ( exclude == 1 .or. exclude == 3 ) then
!         exclude first node
          work(mesh%curves(geometry)%nodes(1)) = val
        end if
        if ( exclude == 2 .or. exclude == 3 ) then
!         exclude last node
          work(mesh%curves(geometry)%nodes(nnodes)) = val
        end if
      end if
      work ( &
        mesh%points ( problem%transformations(trans)%excludepoints ) ) = val
      numec = size ( problem%transformations(trans)%excludecurves )
      do curve = 1, numec
        crvnum = problem%transformations(trans)%excludecurves(curve)
        work ( mesh%curves(crvnum)%nodes ) = val
      end do
      numes = size ( problem%transformations(trans)%excludesurfaces )
      do surface = 1, numes
        srfnum = problem%transformations(trans)%excludesurfaces(surface)
        work ( mesh%surfaces(srfnum)%nodes ) = val
      end do

    end subroutine exclude_nodes


!   determine number of degrees in the nodes (max)

    subroutine ndof_nodes

      integer :: node, nodenr

      ndof = 0

      if ( physq > 0 ) then
        do node = 1, nnodes
          nodenr = nodes(node)
          ndof = max( ndof, problem%vec_nodnumdegfd(nodenr+1,&
                                    &problem%physq(physq)) &
                      - problem%vec_nodnumdegfd(nodenr,&
                                    &problem%physq(physq)) )
        end do
      else
        do node = 1, nnodes
          nodenr = nodes(node)
          ndof = max( ndof, problem%nodnumdegfd(nodenr+1) &
                      - problem%nodnumdegfd(nodenr) )
        end do
      end if

    end subroutine ndof_nodes

    subroutine check_layer

      if ( layer > 0 .and. problem%transformations(trans)%errorlayer < 2 ) then
        if ( .not. btest(problem%nodlayers(nodenr),layer-1) ) then
          select case ( problem%transformations(trans)%errorlayer )
          case(0)
            write(*,'(3(/a)/3(a,i0)/)') &
              'Error in problem_definition_transformations: ', &
              ' Node in transformation not fully within layer.', &
              ' Stop because errorlayer = 0 ', &
              ' transformation = ', trans, ' layer = ', layer, &
              ' nodenr = ', nodenr
            stop
          case(1)
            write(*,'(3(/a)/3(a,i0)/)') &
              'Warning in problem_definition_transformations: ', &
              ' Node in transformation not fully within layer.', &
              ' Warning because errorlayer = 1 ', &
              ' transformation = ', trans, ' layer = ', layer, &
              ' nodenr = ', nodenr
          case default
            call errormsg_case_default ( 'check_layer', &
              'problem%transformations(trans)%errorlayer ', &
              int_value=problem%transformations(trans)%errorlayer  )
          end select
        end if
      end if

    end subroutine check_layer

!   get the nodes for the transformation

    subroutine get_nodes

      select case ( problem%transformations(trans)%typegeometry )
      case(1) ! transformation in point
        nnodes = 1
        nodes = [ mesh%points(geometry) ]
      case(2)
!       transformation on curve
        nnodes = mesh%curves(geometry)%nnodes
        nodes = mesh%curves(geometry)%nodes
      case(3)
!       transformation on surface
        nnodes = mesh%surfaces(geometry)%nnodes
        nodes = mesh%surfaces(geometry)%nodes
      case(4)
!       transformation on volume
        nnodes = mesh%volumes(geometry)%nnodes
        nodes = mesh%volumes(geometry)%nodes
      case(5)
!       transformation on nodeset
        nnodes = size(mesh%nodesets(geometry)%a)
        nodes = mesh%nodesets(geometry)%a
        case default
          call errormsg_case_default ( 'get_nodes', &
            'problem%transformations(trans)%typegeometry', &
            int_value=problem%transformations(trans)%typegeometry )
      end select

    end subroutine get_nodes

    subroutine create_Amat

      integer :: i, j

      n = problem%numdegfd
      problem%Amat%n = n
      problem%Amat%m = n
      allocate ( problem%Amat%ia(n+1) )

!     scan workdeg

      j = 0
      do i = 1, problem%numdegfd
        if ( workdeg(i) == 0 ) then
!         degree not transformed (value of 1 on the diagonal)
          problem%Amat%ia(i+1) = 1
        else
!         degree transformed
          j = j + 1  ! degree in transformation
          problem%degfdtrans(j) = i
          problem%Amat%ia(i+1) = workdeg(i)
        end if
      end do

!     accumulate Amat%ia and compute nnz

      problem%Amat%ia(1) = 1
      do i = 1, n
        problem%Amat%ia(i+1) = problem%Amat%ia(i+1) + problem%Amat%ia(i)
      end do
      problem%Amat%nnz = problem%Amat%ia(problem%Amat%n+1) - 1

!     allocate data of Amat

      allocate( problem%Amat%a(problem%Amat%nnz) )
      allocate( problem%Amat%ja(problem%Amat%nnz) )

!     set orthogonality of Amat

      if ( all( problem%transformations&
                           &(1:problem%numtransformations)%orthogonal) ) then
        problem%orthogonal = .true.
      else
        problem%orthogonal = .false.
      end if

    end subroutine create_Amat


!   create dummy Amat (can only be deleted). Should have been the identity
!   matrix, but this is not stored to save memory for problems without
!   transformations.

    subroutine create_dummy_Amat

      problem%Amat%n = 0
      problem%Amat%m = 0
      problem%Amat%nnz = 0
      allocate( problem%Amat%ia(0) )
      allocate( problem%Amat%a(0) )
      allocate( problem%Amat%ja(0) )
      problem%orthogonal = .true.

    end subroutine create_dummy_Amat

  end subroutine problem_definition_transformations


! Delete single problem

  subroutine delete_single_problem ( problem, nr )

    type(problem_t), intent(inout) :: problem
    integer, intent(in) :: nr

    if ( .not. problem%created ) then
      write(*,'(/a/a,i0/)') &
        'Error delete_problem: problem has not been defined ', &
        ' problem number in heading = ', nr
      stop
    end if

!   deallocate all allocatables and intialize to defaults
    call deall ( problem )

  contains

    subroutine deall ( problem )
      type(problem_t), intent(out) :: problem
    end subroutine deall

  end subroutine delete_single_problem


! Delete problem

  subroutine delete_problem ( problem1, problem2, problem3, problem4, problem5 )

    type(problem_t), intent(inout) :: problem1
    type(problem_t), intent(inout), optional :: problem2, problem3, problem4, &
      problem5

    call delete_single_problem(problem1,1)
    if ( present(problem2) ) call delete_single_problem(problem2,2)
    if ( present(problem3) ) call delete_single_problem(problem3,3)
    if ( present(problem4) ) call delete_single_problem(problem4,4)
    if ( present(problem5) ) call delete_single_problem(problem5,5)

  end subroutine delete_problem

end module problem_m
