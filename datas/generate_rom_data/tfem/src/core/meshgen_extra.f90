
! Copyright (C) 2004-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines extra routines for mesh generation:
!   - mesh_merge: create a new mesh from two other meshes.
!   - mesh_convert: convert a mesh to a new mesh with different basic properties
!   - mesh_change: change a basic mesh (in place).
!   - mesh_skeleton: create a skeleton mesh, i.e. a mesh without coordinates
!                    and topology filled in (mesh is incomplete).

module meshgen_extra_m

  use glob_defs_m
  use mesh_m
  use array_defs_m
  use limits_m, only: MAXPOINTS, MAXCURVES, MAXOBJECTS, MAXSURFACES, &
                      MAXBLOCKS, MAXVOLUMES, MAXNODESETS, MAXELEMENTSETS
  use set_optional_m

  implicit none

contains


! Create a new mesh from two other meshes

  subroutine mesh_merge ( mesh1, mesh2, mesh, point1, point2, points1, &
    points2, curve1, curve2, curves1, curves2, surface1, surface2, surfaces1, &
    surfaces2, mergegroups, mergegroup1, mergegroup2, deletepoints1, &
    deletepoints2, deletecurves1, deletecurves2, deletesurfaces1, &
    deletesurfaces2, nogroupmerge, warn )

!   the two meshes that need to be merged
    type(mesh_t), intent(in) :: mesh1, mesh2

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh

!   points in the corresponding meshes that are the same
    integer, intent(in), optional :: point1, point2

!   multiple points in the corresponding meshes that are the same
!   NOTE: points2(i) corresponds to points1(i)
    integer, intent(in), dimension(:), optional :: points1, points2

!   curves in the corresponding meshes that are the same
!   if curve2 < 0, the curve number is |curve2| and the curve is traversed in
!   the opposite way
    integer, intent(in), optional :: curve1, curve2

!   multiple curves in the corresponding meshes that are the same
!   if curves2(i) < 0, the curve number is |curves2(i)| and the curve
!   is traversed in the opposite way
!   NOTE: curves2(i) corresponds to curves1(i)
    integer, intent(in), dimension(:), optional :: curves1, curves2

!   surfaces in the corresponding meshes that are the same
    integer, intent(in), optional :: surface1, surface2

!   multiple surfaces in the corresponding meshes that are the same
!   NOTE: surfaces2(i) corresponds to surfaces1(i)
    integer, intent(in), dimension(:), optional :: surfaces1, surfaces2

!   groups that need to be merged into single group
!   if not present, groups are not merged and new groups are created
!   mergegroups(:,1) = group in mesh1 becomes a single group in mesh together
!   with mergegroups(:,2) = group in mesh2.
    integer, dimension(:,:), intent(in), optional :: mergegroups

!   groups that need to be merged into single group
!   simplified interface for a single group couple.
!   if not present, groups are not merged and new groups are created
!   mergegroup1 = group in mesh1 becomes a single group in mesh together
!   with mergegroup2 = group in mesh2. If more than one couple of groups needs
!   to be merged, use mergegroups.
!   NOTE: if only mergegroup1 is present, mergegroup2=1 is assumed.
    integer, intent(in), optional :: mergegroup1, mergegroup2

!   arrays to indicate which points of mesh1 and mesh2 must not be included
!   in the new merged mesh
!   NOTE: double points are removed automatically already
    integer, dimension(:), intent(in), optional :: deletepoints1, deletepoints2

!   arrays to indicate which curves of mesh1 and mesh2 must not be included
!   in the new merged mesh
!   NOTE: double curves are removed automatically already
    integer, dimension(:), intent(in), optional :: deletecurves1, deletecurves2

!   arrays to indicate which surfaces of mesh1 and mesh2 must not be included
!   in the new merged mesh
!   NOTE: double surfaces are removed automatically already
    integer, dimension(:), intent(in), optional :: deletesurfaces1, &
      deletesurfaces2

!   do not merge groups if there is only a single group in both mesh1 and
!   mesh2. The default behavior is nogroupmerge=.false., this means that if
!   both mesh1 and mesh2 have one group of elements that are compatible,
!   one single group of elements in the new mesh is created. This is to avoid
!   that in the most common situation extra parameters in the heading are
!   needed to create a single group.
    logical, intent(in), optional :: nogroupmerge

!   Give warnings. By setting to .false. the warnings are suppressed.
!   default=.true.
    logical, intent(in), optional :: warn

!   This routine merges two meshes. Optionally points, curves and/or surfaces
!   can be identified to be the same. Only the basic structure as generated
!   in the basic meshgeneration routines such as quadrilateral2d is copied.
!   The arrays/structures generated in fill_mesh_parts still needs to done.


    logical :: mergegrps, singlegrpm, lwarn
    integer :: nn, nodenr, elgrp, sp, point, curve, crv, elem, grp, i
    integer :: numdelp1, numdelp2, numdelc1, numdelc2, numdels1, numdels2
    integer :: reduction, mergegrp2
    integer :: surface, srf, pnt, pnt2, crv1, crv2, srf1, srf2
    integer, dimension(:), allocatable :: work, gr1, gr2
    integer, dimension(:), allocatable :: delpoints1, delpoints2
    integer, dimension(:), allocatable :: delcurves1, delcurves2
    integer, dimension(:), allocatable :: delsurfaces1, delsurfaces2


!   initialize

    mergegrps = .false.
    singlegrpm = .true.

!   set parameters

    if ( present(nogroupmerge) ) then
      singlegrpm = .not. nogroupmerge
    end if

!   some testing

    if ( mesh%meshgen ) then
      write(*,'(/a/a/)') &
        'Error in mesh_merge: mesh is not empty', &
        ' either use a new mesh_t variable or delete old mesh '
      stop
    end if

    lwarn = set_optional ( variable=warn, default=.true. )

    if ( lwarn ) then

      if ( mesh1%nvolumes > 0 ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_merge: mesh1 contains volumes. ', &
          'Merging of volumes has not been implemented.'
      end if
      if ( mesh2%nvolumes > 0 ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_merge: mesh2 contains volumes. ', &
          'Merging of volumes has not been implemented.'
      end if

      if ( mesh1%nobjects > 0 ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_merge: mesh1 contains objects. ', &
          'Merging of objects has not been implemented.'
      end if
      if ( mesh2%nobjects > 0 ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_merge: mesh2 contains objects. ', &
          'Merging of objects has not been implemented.'
      end if

      if ( mesh1%nblocks > 0 ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_merge: mesh1 contains blocks. ', &
          'Merging of blocks has not been implemented.'
      end if
      if ( mesh2%nblocks > 0 ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_merge: mesh2 contains blocks. ', &
          'Merging of blocks has not been implemented.'
      end if

      if ( mesh1%sblocks%filled ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_merge: mesh1 contains structured blocks. ', &
          'Merging of structured blocks has not been implemented.'
      end if
      if ( mesh2%sblocks%filled ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_merge: mesh2 contains structured blocks. ', &
          'Merging of structured blocks has not been implemented.'
      end if

    end if

    if ( .not. mesh1%meshgen ) then
      write(*,'(/a/)') &
        'Error in mesh_merge: mesh1 is empty'
      stop
    else if ( .not. mesh2%meshgen ) then
      write(*,'(/a/)') &
        'Error in mesh_merge: mesh2 is empty'
      stop
    else if ( mesh1%ndim /= mesh2%ndim ) then
      write(*,'(/a/)') &
        'Error in mesh_merge: ndim different in mesh1 and mesh2'
      stop
    end if

    call check_blend ( mesh1, 'mesh_merge', keyword='mesh1', &
      comment='Merge the meshes before making blended meshes.' )
    call check_blend ( mesh2, 'mesh_merge', keyword='mesh2', &
      comment='Merge the meshes before making blended meshes.' )

    if ( present(mergegroups) ) then

      if ( present(mergegroup1) .or. present(mergegroup2) ) then
        write(*,'(/a/a/)') &
          'Error in mesh_merge: mergegroup1, mergegroup2 can not be present', &
          ' if mergegroups is present'
        stop
      end if

      if ( size(mergegroups,dim=1) > min(mesh1%nelgrp,mesh2%nelgrp) ) then
        write(*,'(/a/a,i0,2a,i0/)') &
          'Error in mesh_merge: size mergegroups incorrect: ', &
          ' first dimension is ', size(mergegroups,dim=1), &
          ' whereas it should be smaller than ', &
          'min(mesh1%nelgrp,mesh2%nelgrp) = ', min(mesh1%nelgrp,mesh2%nelgrp)
        stop
      end if

      if ( size(mergegroups,dim=2) /= 2 ) then
        write(*,'(/a/a,i0,a/)') &
          'Error in mesh_merge: size mergegroups incorrect: ', &
          ' second dimension is ', size(mergegroups,dim=2), &
          ' whereas it should be 2.'
        stop
      end if

      if ( any(mergegroups <= 0) .or. any(mergegroups(:,1) > mesh1%nelgrp) &
            .or. any(mergegroups(:,2) > mesh2%nelgrp) ) then
        write(*,'(/a,2(/a,i0)/)') &
          'Error in mesh_merge: mergegroups has incorrect values ', &
          'values <= 0 or values > number of groups of mesh1 = ', &
           mesh1%nelgrp, &
          'values <= 0 or values > number of groups of mesh2 = ', mesh2%nelgrp
        stop
      end if

      if ( any(mesh1%element(mergegroups(:,1))%elshape /= &
               mesh2%element(mergegroups(:,2))%elshape ) .or. &
           any(mesh1%element(mergegroups(:,1))%numnod /= &
               mesh2%element(mergegroups(:,2))%numnod ) ) then
        write(*,'(/a/a/)') &
          'Error in mesh_merge: mergegroups has incompatible groups', &
          ' element shapes and/or number of nodes are different'
        stop
      end if

      allocate( gr1(mesh1%nelgrp), gr2(mesh2%nelgrp) )
      gr1 = 0
      gr1(mergegroups(:,1)) = mergegroups(:,2)
      gr2 = 0
      gr2(mergegroups(:,2)) = mergegroups(:,1)

      mergegrps = .true.
      reduction = size(mergegroups,dim=1)

    else if ( present(mergegroup1) ) then

      if ( mergegroup1 <= 0 .or. mergegroup1 > mesh1%nelgrp ) then
        write(*,'(/a/a,i0/)') &
          'Error in mesh_merge: mergegroup1 has an incorrect value ', &
          'value <= 0 or value > number of groups of mesh1 = ', mesh1%nelgrp
        stop
      end if

      if ( present(mergegroup2) ) then
        if ( mergegroup2 <= 0 .or. mergegroup2 > mesh2%nelgrp ) then
          write(*,'(/a/a,i0/)') &
            'Error in mesh_merge: mergegroup2 has an incorrect value ', &
            'value <= 0 or value > number of groups of mesh2 = ', mesh2%nelgrp
          stop
        end if
        mergegrp2 = mergegroup2
      else
        mergegrp2 = 1   ! default
      end if

      if ( mesh1%element(mergegroup1)%elshape /= &
           mesh2%element(mergegrp2)%elshape .or. &
           mesh1%element(mergegroup1)%numnod /= &
           mesh2%element(mergegrp2)%numnod  ) then
        write(*,'(/3(a/))') &
          'Error in mesh_merge: mergegroup1 and mergegroup2 ', &
          'are incompatible:', &
          ' element shapes and/or number of nodes are different'
        stop
      end if

      allocate( gr1(mesh1%nelgrp), gr2(mesh2%nelgrp) )
      gr1 = 0
      gr1(mergegroup1) = mergegrp2
      gr2 = 0
      gr2(mergegrp2) = mergegroup1

      mergegrps = .true.
      reduction = 1

    else if ( mesh1%nelgrp == 1 .and. mesh2%nelgrp == 1 .and. singlegrpm ) then

!     single group: merge if allowed

      if ( mesh1%element(1)%elshape /= &
           mesh2%element(1)%elshape .or. &
           mesh1%element(1)%numnod /= &
           mesh2%element(1)%numnod  ) then
        if ( lwarn ) write(*,'(/3(a/))') &
          'Warning in mesh_merge: elements in mesh1 and mesh2 ', &
          'are incompatible:', &
          ' two groups are created in the new mesh '
      else

        allocate( gr1(1), gr2(1) )
        gr1 = 1
        gr2 = 1

        mergegrps = .true.
        reduction = 1

      end if

    end if

!   work storage for new numbering of new nodes for mesh2

    allocate(work(mesh2%nnodes))

    work = 0

!   arrays for deleting points when merging

    allocate(delpoints1(mesh1%npoints))
    delpoints1 = 0
    allocate(delpoints2(mesh2%npoints))
    delpoints2 = 0

!   arrays for deleting curves when merging

    allocate(delcurves1(mesh1%ncurves))
    delcurves1 = 0
    allocate(delcurves2(mesh2%ncurves))
    delcurves2 = 0

!   arrays for deleting surfaces when merging

    allocate(delsurfaces1(mesh1%nsurfaces))
    delsurfaces1 = 0
    allocate(delsurfaces2(mesh2%nsurfaces))
    delsurfaces2 = 0


!   determine merged geometrical objects

!   POINTS

    if ( present(point1) ) then

!     merge points

      if ( point1 <=0 .or. point1 > mesh1%npoints ) then
        write(*,'(/a/a,i0/)') &
          'Error: point1 in heading of mesh_merge has wrong value:', &
          'point1 <=0 or point1 > number of points = ', mesh1%npoints
        stop
      end if

      if ( .not. present(point2) ) then
        write(*,'(/a/)') &
          'Error: point2 missing from heading of mesh_merge'
        stop
      else if ( point2 <=0 .or. point2 > mesh2%npoints ) then
        write(*,'(/a/a,i0/)') &
          'Error: point2 in heading of mesh_merge has wrong value:', &
          'point2 <=0 or point2 > number of points = ', &
          mesh2%npoints
        stop
      end if

      delpoints2(point2) = 1 ! always delete point2

!     node in point2 can be removed == node in point1

      work(mesh2%points(point2)) = mesh1%points(point1)

    end if

!   multiple points

    if ( present(points1) ) then

!     merge points

      if ( any(points1 <=0 ) .or. any(points1 > mesh1%npoints) ) then
        write(*,'(/a/a,i0/)') &
          'Error: points1 in heading of mesh_merge has wrong value:', &
          'points1 <=0 or points1 > number of points = ', mesh1%npoints
        stop
      end if

      if ( .not. present(points2) ) then
        write(*,'(/a/)') &
          'Error: points2 missing from heading of mesh_merge'
        stop
      else if ( size(points2) /= size(points1) ) then
        write(*,'(/a/)') &
          'Error: points2 and points1 in the heading of mesh_merge must', &
          ' have the same dimension.'
        stop
      else if ( any(points2 == 0) .or. any(abs(points2) > mesh2%npoints) ) then
        write(*,'(/a/a,i0/)') &
          'Error: points2 in heading of mesh_merge has wrong value:', &
          ' points2 ==0 or points2 > number of points = ', &
          mesh2%npoints
        stop
      end if

      delpoints2(points2) = 1 ! always delete points in points2

!     node in points2 can be removed == node in points1

      work(mesh2%points(points2)) = mesh1%points(points1)

    end if

!   delete points?

    if ( present(deletepoints1) ) then
      if ( any(deletepoints1 <= 0) .or. any(deletepoints1 > mesh1%npoints) &
        ) then
        write(*,'(/a/a,i0/)') &
          'Error: deletepoints1 in heading of mesh_merge has wrong values:', &
          'point <=0 or point > number of points = ', mesh1%npoints
        stop
      end if
      delpoints1(deletepoints1) = 1
    end if

    if ( present(deletepoints2) ) then
      if ( any(deletepoints2 <= 0) .or. any(deletepoints2 > mesh2%npoints) &
        ) then
        write(*,'(/a/a,i0/)') &
          'Error: deletepoints2 in heading of mesh_merge has wrong values:', &
          'point <=0 or point > number of points = ', mesh2%npoints
        stop
      end if
      delpoints2(deletepoints2) = 1
    end if

    numdelp1 = count( delpoints1 == 1 )
    numdelp2 = count( delpoints2 == 1 )

!   CURVES

    if ( present(curve1) ) then

!     merge curves

      if ( curve1 <=0 .or. curve1 > mesh1%ncurves ) then
        write(*,'(/a/a,i0/)') &
          'Error: curve1 in heading of mesh_merge has wrong value:', &
          'curve1 <=0 or curve1 > number of curves = ', mesh1%ncurves
        stop
      end if

      if ( .not. present(curve2) ) then
        write(*,'(/a/)') &
          'Error: curve2 missing from heading of mesh_merge'
        stop
      else if ( curve2 == 0 .or. abs(curve2) > mesh2%ncurves ) then
        write(*,'(/a/a,i0/)') &
          'Error: curve2 in heading of mesh_merge has wrong value:', &
          ' curve2 ==0 or curve2 > number of curves = ', &
          mesh2%ncurves
        stop
      else if ( mesh1%curves(curve1)%nnodes /= &
                     mesh2%curves(abs(curve2))%nnodes ) then
        write(*,'(/a/a/)') &
          'Error: curve1 and curve2 in heading of mesh_merge have a ', &
          'different number of nodes'
        stop
      end if

      delcurves2(abs(curve2)) = 1 ! always delete curve2

!     nodes on curve2 can be removed == nodes on curve1

      if ( curve2 > 0 ) then
        work(mesh2%curves(curve2)%nodes) = mesh1%curves(curve1)%nodes
      else
        work(mesh2%curves(-curve2)%nodes(mesh2%curves(-curve2)%nnodes:1:-1)) &
                      = mesh1%curves(curve1)%nodes
      end if

    end if

!   multiple curves

    if ( present(curves1) ) then

!     merge curves

      if ( any(curves1 <=0 ) .or. any(curves1 > mesh1%ncurves) ) then
        write(*,'(/a/a,i0/)') &
          'Error: curves1 in heading of mesh_merge has wrong value:', &
          'curves1 <=0 or curves1 > number of curves = ', mesh1%ncurves
        stop
      end if

      if ( .not. present(curves2) ) then
        write(*,'(/a/)') &
          'Error: curves2 missing from heading of mesh_merge'
        stop
      else if ( size(curves2) /= size(curves1) ) then
        write(*,'(/a/)') &
          'Error: curves2 and curves1 in the heading of mesh_merge must', &
          ' have the same dimension.'
        stop
      else if ( any(curves2 == 0) .or. any(abs(curves2) > mesh2%ncurves) ) then
        write(*,'(/a/a,i0/)') &
          'Error: curves2 in heading of mesh_merge has wrong value:', &
          ' curves2 ==0 or curves2 > number of curves = ', &
          mesh2%ncurves
        stop
      else if ( any(mesh1%curves(curves1)%nnodes /= &
                     mesh2%curves(abs(curves2))%nnodes) ) then
        write(*,'(/a/a/)') &
          'Error: curves1 and curves2 in heading of mesh_merge have a ', &
          'different number of nodes'
        stop
      end if

      delcurves2(abs(curves2)) = 1 ! always delete curves in curves2

!     nodes on curve2 can be removed == nodes on curve1

      do crv = 1, size(curves2)
        crv1 = curves1(crv)
        crv2 = curves2(crv)
        if ( crv2 > 0 ) then
          work(mesh2%curves(crv2)%nodes) = mesh1%curves(crv1)%nodes
        else
          work(mesh2%curves(-crv2)%nodes(mesh2%curves(-crv2)%nnodes:1:-1)) &
                        = mesh1%curves(crv1)%nodes
        end if
      end do

    end if

!   delete curves?

    if ( present(deletecurves1) ) then
      if ( any(deletecurves1 <= 0) .or. any(deletecurves1 > mesh1%ncurves) &
        ) then
        write(*,'(/a/a,i0/)') &
          'Error: deletecurves1 in heading of mesh_merge has wrong values:', &
          'curve <=0 or curve > number of curves = ', mesh1%ncurves
        stop
      end if
      delcurves1(deletecurves1) = 1
    end if

    if ( present(deletecurves2) ) then
      if ( any(deletecurves2 <= 0) .or. any(deletecurves2 > mesh2%ncurves) &
        ) then
        write(*,'(/a/a,i0/)') &
          'Error: deletecurves2 in heading of mesh_merge has wrong values:', &
          'curve <=0 or curve > number of curves = ', mesh2%ncurves
        stop
      end if
      delcurves2(deletecurves2) = 1
    end if

    numdelc1 = count( delcurves1 == 1 )
    numdelc2 = count( delcurves2 == 1 )


!   SURFACES

    if ( present(surface1) ) then

!     merge surfaces

      if ( surface1 <=0 .or. surface1 > mesh1%nsurfaces ) then
        write(*,'(/a/a,i0/)') &
          'Error: surface1 in heading of mesh_merge has wrong value:', &
          'surface1 <=0 or surface1 > number of surfaces = ', mesh1%nsurfaces
        stop
      end if

      if ( .not. present(surface2) ) then
        write(*,'(/a/)') &
          'Error: surface2 missing from heading of mesh_merge'
        stop
      else if ( surface2 <= 0 .or. surface2 > mesh2%nsurfaces ) then
        write(*,'(/a/a,i0/)') &
          'Error: surface2 in heading of mesh_merge has wrong value:', &
          ' surface2 <=0 or surface2 > number of surfaces = ', &
          mesh2%nsurfaces
        stop
      else if ( mesh1%surfaces(surface1)%nnodes /= &
                     mesh2%surfaces(abs(surface2))%nnodes ) then
        write(*,'(/a/a/)') &
          'Error: surface1 and surface2 in heading of mesh_merge have a ', &
          'different number of nodes'
        stop
      end if

      delsurfaces2(surface2) = 1 ! always delete surface2

!     nodes on surface2 can be removed == nodes on surface1

      work(mesh2%surfaces(surface2)%nodes) = mesh1%surfaces(surface1)%nodes

    end if

!   multiple surfaces

    if ( present(surfaces1) ) then

!     merge surfaces

      if ( any(surfaces1 <=0 ) .or. any(surfaces1 > mesh1%nsurfaces) ) then
        write(*,'(/a/a,i0/)') &
          'Error: surfaces1 in heading of mesh_merge has wrong value:', &
          'surfaces1 <=0 or surfaces1 > number of surfaces = ', mesh1%nsurfaces
        stop
      end if

      if ( .not. present(surfaces2) ) then
        write(*,'(/a/)') &
          'Error: surfaces2 missing from heading of mesh_merge'
        stop
      else if ( size(surfaces2) /= size(surfaces1) ) then
        write(*,'(/a/)') &
          'Error: surfaces2 and surfaces1 in the heading of mesh_merge must', &
          ' have the same dimension.'
        stop
      else if ( any(surfaces2 == 0) .or. &
                any(abs(surfaces2) > mesh2%nsurfaces) ) then
        write(*,'(/a/a,i0/)') &
          'Error: surfaces2 in heading of mesh_merge has wrong value:', &
          ' surfaces2 ==0 or surfaces2 > number of surfaces = ', &
          mesh2%nsurfaces
        stop
      else if ( any(mesh1%surfaces(surfaces1)%nnodes /= &
                     mesh2%surfaces(abs(surfaces2))%nnodes) ) then
        write(*,'(/a/a/)') &
          'Error: surfaces1 and surfaces2 in heading of mesh_merge have a ', &
          'different number of nodes'
        stop
      end if

      delsurfaces2(abs(surfaces2)) = 1 ! always delete surfaces in surfaces2

!     nodes on surface2 can be removed == nodes on surface1

      do srf = 1, size(surfaces2)
        srf1 = surfaces1(srf)
        srf2 = surfaces2(srf)
        work(mesh2%surfaces(srf2)%nodes) = mesh1%surfaces(srf1)%nodes
      end do

    end if

!   delete surfaces?

    if ( present(deletesurfaces1) ) then
      if ( any(deletesurfaces1 <= 0) .or. &
                    any(deletesurfaces1 > mesh1%nsurfaces) ) then
        write(*,'(/2a/a,i0/)') &
          'Error: deletesurfaces1 in heading of mesh_merge has', &
          ' wrong values:', &
          'surface <=0 or surface > number of surfaces = ', mesh1%nsurfaces
        stop
      end if
      delsurfaces1(deletesurfaces1) = 1
    end if

    if ( present(deletesurfaces2) ) then
      if ( any(deletesurfaces2 <= 0) .or. &
                  any(deletesurfaces2 > mesh2%nsurfaces) ) then
        write(*,'(/2a/a,i0/)') &
          'Error: deletesurfaces2 in heading of mesh_merge has', &
          ' wrong values:', &
          'surface <=0 or surface > number of surfaces = ', mesh2%nsurfaces
        stop
      end if
      delsurfaces2(deletesurfaces2) = 1
    end if

    numdels1 = count( delsurfaces1 == 1 )
    numdels2 = count( delsurfaces2 == 1 )


    nn = mesh1%nnodes

!   new node numbers of mesh2
    do nodenr = 1, mesh2%nnodes
      if ( work(nodenr) > 0 ) cycle   ! has new node number from mesh1
      nn = nn + 1
      work(nodenr) = nn
    end do

!   fill new mesh

    mesh%ndim = mesh1%ndim
    mesh%nelem = mesh1%nelem + mesh2%nelem
    mesh%nnodes = nn

    if ( mergegrps ) then

!     merge some groups of mesh1 and mesh2

      mesh%nelgrp = mesh1%nelgrp + mesh2%nelgrp - reduction

      allocate( mesh%element ( mesh%nelgrp ) )
      allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
      mesh%nnodes_blend = [0,mesh%nnodes]
      allocate( mesh%grpnumel ( mesh%nelgrp ) )
      allocate( mesh%topology ( mesh%nelgrp ) )

!     Fill element and topology

!     copy elements mesh1 and merged groups form mesh2
      do elgrp = 1, mesh1%nelgrp
        call copy ( mesh1%element(elgrp), mesh%element(elgrp) )
        if ( gr1(elgrp) > 0 ) then
!         merge group form mesh2
          mesh%grpnumel(elgrp) = mesh1%grpnumel(elgrp) + &
                                 mesh2%grpnumel(gr1(elgrp))
          allocate( mesh%topology(elgrp)%a(mesh%element(elgrp)%numnod,&
                                           mesh%grpnumel(elgrp)) )
          mesh%topology(elgrp)%a(:,1:mesh1%grpnumel(elgrp)) = &
                                            mesh1%topology(elgrp)%a
          do elem = 1, mesh2%grpnumel(gr1(elgrp))
            mesh%topology(elgrp)%a(:,mesh1%grpnumel(elgrp)+elem) = &
                                  work(mesh2%topology(gr1(elgrp))%a(:,elem))
          end do
        else
!         copy group form mesh1
          mesh%grpnumel(elgrp) = mesh1%grpnumel(elgrp)
          mesh%topology(elgrp)%a = mesh1%topology(elgrp)%a
        end if
      end do

!     copy remaining elements mesh2
      sp = mesh1%nelgrp
      grp = 0
      do elgrp = 1, mesh2%nelgrp
        if ( gr2(elgrp) > 0 ) cycle  ! group is already merged
        grp = grp + 1
        call copy ( mesh2%element(elgrp), mesh%element(sp+grp) )
        mesh%grpnumel(sp+grp) = mesh2%grpnumel(elgrp)
        allocate( mesh%topology(sp+grp)%a(mesh2%element(elgrp)%numnod,&
                                         mesh2%grpnumel(elgrp)) )
        do elem = 1, mesh2%grpnumel(elgrp)
          mesh%topology(sp+grp)%a(:,elem) = &
                work(mesh2%topology(grp)%a(:,elem))
        end do
      end do

      deallocate ( gr1, gr2 )

    else

!     create separate groups for mesh1 and mesh2

      mesh%nelgrp = mesh1%nelgrp + mesh2%nelgrp

      allocate( mesh%element ( mesh%nelgrp ) )
      allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
      mesh%nnodes_blend = [0,mesh%nnodes]
      allocate( mesh%grpnumel ( mesh%nelgrp ) )
      allocate( mesh%topology ( mesh%nelgrp ) )

!     Fill element and topology

!     copy elements mesh1
      do elgrp = 1, mesh1%nelgrp
        call copy ( mesh1%element(elgrp), mesh%element(elgrp) )
        mesh%grpnumel(elgrp) = mesh1%grpnumel(elgrp)
        mesh%topology(elgrp)%a = mesh1%topology(elgrp)%a
      end do

!     copy elements mesh2
      sp = mesh1%nelgrp
      do elgrp = 1, mesh2%nelgrp
        call copy ( mesh2%element(elgrp), mesh%element(sp+elgrp) )
        mesh%grpnumel(sp+elgrp) = mesh2%grpnumel(elgrp)
        allocate( mesh%topology(sp+elgrp)%a(mesh2%element(elgrp)%numnod,&
                                         mesh2%grpnumel(elgrp)) )
        do elem = 1, mesh2%grpnumel(elgrp)
          mesh%topology(sp+elgrp)%a(:,elem) = &
                work(mesh2%topology(elgrp)%a(:,elem))
        end do
      end do

    end if

!   fill points

!   remove all double points

    do pnt2 = 1, mesh2%npoints
      if ( delpoints2(pnt2) == 1 ) cycle ! already deleted
      if ( mesh2%points(pnt2) == 0 ) cycle ! point not connected
!     remove double points
      if ( any( mesh1%points(:mesh1%npoints) == &
                                      work(mesh2%points(pnt2)) ) ) then
        delpoints2(pnt2) = 1
        numdelp2 = numdelp2 + 1
      end if
    end do

    mesh%npoints = mesh1%npoints + mesh2%npoints - numdelp1 - numdelp2
    sp = mesh1%npoints - numdelp1

    allocate(mesh%points(max(size(mesh1%points),mesh%npoints)))

!   copy points mesh1

    pnt = 0

    do point = 1, mesh1%npoints

      if ( delpoints1(point) == 1 ) cycle ! no copy of point

      pnt = pnt + 1

      mesh%points(pnt) = mesh1%points(point)

    end do

!   copy points mesh2

    pnt = 0

    do point = 1, mesh2%npoints

      if ( delpoints2(point) == 1 ) cycle ! no copy of point

      pnt = pnt + 1

      mesh%points(sp+pnt) = work(mesh2%points(point)) ! use new node numbers

    end do

!   fill curves

!   remove all double curves

    do crv2 = 1, mesh2%ncurves
      if ( delcurves2(crv2) == 1 ) cycle ! already deleted
!     remove double curves
      do crv1 = 1, mesh1%ncurves
        if ( mesh1%curves(crv1)%nnodes /= mesh2%curves(crv2)%nnodes ) cycle
        if ( all( mesh1%curves(crv1)%nodes == &
                         work(mesh2%curves(crv2)%nodes) ) .or. &
             all( mesh1%curves(crv1)%nodes(mesh1%curves(crv1)%nnodes:1:-1) &
                               == work(mesh2%curves(crv2)%nodes) ) ) then
          delcurves2(crv2) = 1
          numdelc2 = numdelc2 + 1
        end if
      end do
    end do

    mesh%ncurves = mesh1%ncurves + mesh2%ncurves - numdelc1 - numdelc2
    sp = mesh1%ncurves - numdelc1

    allocate(mesh%curves(max(size(mesh1%curves),mesh%ncurves)))

!   copy curves mesh1
    crv = 0

    do curve = 1, mesh1%ncurves

      if ( delcurves1(curve) == 1 ) cycle ! no copy of curve

      crv = crv + 1

      call copy ( mesh1%curves(curve), mesh%curves(crv) )

    end do

!   copy curves mesh2
    crv = 0

    do curve = 1, mesh2%ncurves

      if ( delcurves2(curve) == 1 ) cycle ! no copy of curve

      crv = crv + 1

!     copy nodes, topology with new node numbers

      call copy ( mesh2%curves(curve), mesh%curves(sp+crv), work )

    end do


!   fill surfaces

!   remove all double surfaces

    do srf2 = 1, mesh2%nsurfaces
      if ( delsurfaces2(srf2) == 1 ) cycle ! already deleted
!     remove double surfaces
      do srf1 = 1, mesh1%nsurfaces
        if ( mesh1%surfaces(srf1)%nnodes /= mesh2%surfaces(srf2)%nnodes ) cycle
        if ( all( mesh1%surfaces(srf1)%nodes == &
                         work(mesh2%surfaces(srf2)%nodes) ) ) then
          delsurfaces2(srf2) = 1
          numdels2 = numdels2 + 1
        end if
      end do
    end do

!   remove surfaces from mesh1 and mesh2
    mesh%nsurfaces = mesh1%nsurfaces + mesh2%nsurfaces - numdels1 - numdels2
    sp = mesh1%nsurfaces - numdels1

    allocate(mesh%surfaces(max(size(mesh1%surfaces),mesh%nsurfaces)))

!   copy surfaces mesh1
    srf = 0

    do surface = 1, mesh1%nsurfaces

      if ( delsurfaces1(surface) == 1 ) cycle ! no copy of surface

      srf = srf + 1

      call copy ( mesh1%surfaces(surface), mesh%surfaces(srf) )

    end do

!   copy surfaces mesh2
    srf = 0

    do surface = 1, mesh2%nsurfaces

      if ( delsurfaces2(surface) == 1 ) cycle ! no copy of surface

      srf = srf + 1

!     copy nodes, topology with new node numbers

      call copy ( mesh2%surfaces(surface), mesh%surfaces(sp+srf), work )

    end do

!   coordinates

    allocate(mesh%coor(nn,mesh%ndim))

!   coordinates mesh1
    mesh%coor(:mesh1%nnodes,:) = mesh1%coor

!   coordinates mesh2
    do i = 1, mesh%ndim
      where ( work > mesh1%nnodes )
        mesh%coor(work,i) = mesh2%coor(:,i)
      end where
    end do

!   volumes

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    mesh%meshgen = .true.

    deallocate(work)
    deallocate(delpoints1,delpoints2)
    deallocate(delcurves1,delcurves2,delsurfaces1,delsurfaces2)

  end subroutine mesh_merge


! Convert a mesh to a new basic mesh with different properties
! NOTE: a call with no options just makes a copy of the basic mesh.

  subroutine mesh_convert ( mesh1, mesh, coordinates, mirror, orientation, &
    surfaces, elementshapes, userelmesh, remove_double_nodes, excludepoints, &
    excludecurves, excludesurfaces, xmin, xmax, remove_isolated_nodes, &
    blendmesh, blend, only_groups, from_curve, from_surface, from_volume, warn )

!   the mesh that needs to be converted
    type(mesh_t), intent(in) :: mesh1

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh

!   convert coordinates of the mesh:
!   size(coordinates) gives the space dimension of the new mesh.
!   Note, that a change in dimension is only possible for:
!     line elements: space dimension can be 1, 2 or 3
!     triangles and quadrilateral, space dimension can be 2 or 3
!   |coordinates(i)| gives the old coordinate direction of the new ith
!     direction. A value of 0 means setting all coordinate values to zero
!     in the ith space direction, which is only allowed for lines, triangles
!     and quads.
!   sign(coordinates(i)) identifies whether the new ith direction is in
!     the positive of negative direction of the old coordinates
!   Examples:
!      coordinates=(/2,3,1/): cyclic rotation of the axis
!      coordinates=(/1,2,-3/): mirroring with respect to the (x,y) plane
!      coordinates=(/1,2,0/), 3D mesh: projection on the (x,y) plane
!      coordinates=(/1,2,0/), 2D mesh: convert 2D mesh to 3D with z=0.
!      coordinates=(/1,2/), 3D mesh: projection on the (x,y) plane and
!                                    convert to a 2D mesh
    integer, dimension(:), intent(in), optional :: coordinates

!   For argument 'coordinates' convert the topology of the elements in addition
!   to converting the coordinates of the nodes, assuming that a mirror
!   operation has been performed on the coordinates. Using mirror=.true. the
!   topology (sequence of the nodes) becomes standard again and the sign of
!   detF is positive.
    logical, intent(in), optional :: mirror

!   Change the orientation of the specified surfaces.
!   This option is meant for surfaces made by the routine hexahedron
!   and does not work for meshes from other sources.
!   The sign of orientation denotes the "view direction", where positive
!   means from the usual positive direction (as defined in the userguide
!   for all surfaces S1,...,S6) and negative means looking from the "flip"
!   side, where the first curve stays at the "bottom". Within that view,
!   usual or flipped, |orientation| means the following orientation:
!
!         ___________             ___________
!        |           |           |           |
!        |           |           |           |
!      ^ |           |           |           | ^
!    y | |           |           |           | | x
!      | |___________|           |___________| |
!         --> x                        y <---
!
!     |orientation| = 1         |orientation| = 2
!
!
!              <--- x             --> y
!         ___________             ___________
!        |           | |       | |           |
!        |           | | y   x | |           |
!        |           | v       v |           |
!        |           |           |           |
!        |___________|           |___________|
!
!     |orientation| = 3         |orientation| = 4
!
!   Hence orientation=1 means: no change.
!   Example: orientation=(/2,-1,-4/), surfaces=(/3,4,6/) will change the
!   orientation (nodal numbering and topology) of the surfaces 3, 4 and 6,
!   using orientation=2, -1 and -3, respectively.
!
    integer, dimension(:), intent(in), optional :: orientation, surfaces

!   Convert (some) element shapes to other element shapes:
!   Possibilities:
!
!     'lagrange' : convert serendipity elements (elshape=30,31) to
!                  Lagrangian interpolated elements (elshape=6,14).
!                  All other elements are just copied.
!
!     'serendipity' : convert Lagrange elements (elshape=6,14) to
!                     serendipity elements (elshape=30,31).
!                     All other elements are just copied.
!
!     'extend' : convert "standard" elements to "extended" elements. The
!                following elshapes are converted:
!                    before   after
!                      3       10
!                      4        7
!                      5        9
!                     11       18
!                     12       16
!                     13       17
!                 All other elements are just copied.
!
!     'extend14' : convert "standard" element with elshape=12 (10-node
!                  tetrahedron) to "extended" element with elshape=15 (14-node
!                  tetrahedron).
!                  All other elements are just copied.
!
!     'reduce' : convert "extended" elements to "standard" elements. The
!                following elshapes are converted:
!                    before   after
!                     10        3
!                      7        4
!                      9        5
!                     18       11
!                     15       12
!                     16       12
!                     17       13
!                 All other elements are just copied.
!
!     'spectraltolinear' : convert high-order (elshape=101,102,103,104)
!                and macro high-order elements (elshape=105,106,107)
!                to linear elements (elshape=1,5,13,3).
!                All other elements are just copied.
!
!     'quadratictolinear' : convert "quadratic" elements to "linear" elements
!                by removing all nodes, except the vertices.
!                The following elshapes are converted:
!                    before   after
!                      2        1
!                      4        3
!                      6        5
!                      7        3
!                     12       11
!                     14       13
!                     15       11
!                     30        5
!                     31       13
!                     42       41
!                     43       41
!                     52       51
!                     53       51
!                 All elements in the groups must be of "quadratic" shape.
!
    character(len=*), intent(in), optional :: elementshapes

!   Convert mesh1 to mesh with some elements removed and/or replaced with
!   elmesh using the user routine userelmesh.
!   NOTE: points/curves/surfaces are removed.
!   NOTE: the purpose of this option is basically for plotting xfem data
!   Heading parameters routine for element level:
!   Output for element (elgrp, elem):
!     indicator:
!       indicator = 0 : copy element
!       indicator = 1 : delete element
!       indicator = 2 : replace mesh with elmesh
!     elmesh: mesh for indicator=2
    optional :: userelmesh
    interface
      subroutine userelmesh ( mesh, elgrp, elem, indicator, elmesh )
        use kind_defs_m
        use mesh_m, only: mesh_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: elgrp, elem
        integer, intent(out) :: indicator
        type(mesh_t), intent(inout) :: elmesh
      end subroutine userelmesh
    end interface

!   By setting to .true. the nodes having the same coordinates are
!   replaced by a single node.
!   NOTE: this option requires nodblocks to be set in mesh1 using the
!   routine add_to_mesh.
!   default=.false.
    logical, intent(in), optional :: remove_double_nodes

!   exclude nodes in points for removal with remove_double_nodes=.true.
!   For example excludepoints=(/1,3/) excludes the nodes given by
!   points 1 and 3.
    integer, intent(in), dimension(:), optional :: excludepoints

!   exclude nodes on curves for removal with remove_double_nodes=.true.
!   For example excludecurves=(/1,3/) excludes the nodes given by
!   curves 1 and 3.
    integer, intent(in), dimension(:), optional :: excludecurves

!   exclude nodes on surfaces for removal with remove_double_nodes=.true.
!   For example excludesurfaces=(/1,3/) excludes the nodes given by
!   surfaces 1 and 3.
    integer, intent(in), dimension(:), optional :: excludesurfaces

!   exclude nodes if any of the coordinates is less than the ones given
!   by xmin or any of the coordinates is larger than the ones given by
!   xmax for removal with remove_double_nodes=.true.
!   For example xmin=(/-1._dp,-2._dp/) excludes all nodes that have an
!   x-coordinate less than -1 and/or a y-coordinate less than -2.
    real(dp), intent(in), dimension(:), optional :: xmin, xmax

!   If set to .true. a mesh is converted to a new basic mesh with
!   isolated nodes removed. Isolated nodes are nodes that are not connected
!   to the internal mesh.
!   Note, that geometries (points, curves, surfaces) that are connected to
!   isolated nodes are removed as well.
!   Note, that also points connected to nodal point zero (sepran) are removed.
!   default=.false.
    logical, intent(in), optional :: remove_isolated_nodes

!   Convert mesh1 to a mesh with blendmesh added as blended mesh.
!   The mesh blendmesh should be basically the same as mesh1, but with
!   different shapes for the elements.
    type(mesh_t), intent(in), optional :: blendmesh

!   Convert mesh1 to a mesh with a single shape of elements for each group.
!   blend=0 extracts the main mesh given by the elements in mesh1%element,
!   blend>0 extracts the blend mesh given by the elements in
!           mesh1%element_blend(:,blend),
    integer, intent(in), optional :: blend

!   Convert mesh1 to a mesh with only the specified groups present.
!   For example: only_groups=[1,3].
!   Note, that geometries (points, curves, surfaces) which have all nodes within
!   the specified groups are preserved. Others are discarded. The numbering of
!   the geometries might be different.
    integer, intent(in), dimension(:), optional :: only_groups

!   Convert mesh1 to a mesh created from a single curve in mesh1.
!   For example: from_curve = 3
!   Note that there will only be a mesh consisting of a single group
!   of line elements and no points, curves, surfaces and/or volumes.
    integer, intent(in), optional :: from_curve

!   Convert mesh1 to a mesh made from a single surface in mesh1.
!   For example: from_surface = 3
!   Note that there will only be a mesh consisting of a single group
!   of surface elements and no points, curves, surfaces and/or volumes.
    integer, intent(in), optional :: from_surface

!   Convert mesh1 to a mesh made from a single volume in mesh1.
!   For example: from_volume = 3
!   Note that there will only be a mesh consisting of a single group
!   of volume elements and no points, curves, surfaces and/or volumes.
    integer, intent(in), optional :: from_volume

!   Give warnings. By setting to .false. the warnings are suppressed.
!   default=.true.
    logical, intent(in), optional :: warn


    logical :: extend, lwarn, lremove_double_nodes, lremove_isolated_nodes
    integer :: numkey, elgrp, grp


!   some testing

    if ( mesh%meshgen ) then
      write(*,'(/a/a/)') &
        'Error in mesh_convert: mesh is not empty', &
        ' either use a new mesh_t variable or delete old mesh '
      stop
    end if

    lwarn = set_optional ( variable=warn, default=.true. )

    if ( lwarn ) then

      if ( mesh1%nvolumes > 0 ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_convert: mesh1 contains volumes. ', &
          'Converting of volumes has not been implemented.'
      end if

      if ( mesh1%nobjects > 0 ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_convert: mesh1 contains objects. ', &
          'Converting of objects has not been implemented.'
      end if

      if ( mesh1%nblocks > 0 ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_convert: mesh1 contains blocks. ', &
          'Converting of blocks has not been implemented.'
      end if

      if ( mesh1%nnodesets > 0 ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_convert: mesh1 contains nodesets. ', &
          'Converting of nodesets has not been implemented.'
      end if

      if ( mesh1%nelementsets > 0 ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_convert: mesh1 contains elementsets. ', &
          'Converting of elementsets has not been implemented.'
      end if

      if ( mesh1%sblocks%filled ) then
        write(*,'(/2(a/))') &
          'Warning in mesh_convert: mesh1 contains structured blocks. ', &
          'Converting of structured blocks has not been implemented.'
      end if

    end if

    if ( .not. mesh1%meshgen ) then
      write(*,'(/a/)') &
        'Error in mesh_convert: mesh1 is empty'
      stop
    end if

    numkey = 0

!   convert shape of elements

    if ( present(elementshapes) ) then

      call check_blend ( mesh1, 'mesh_convert', keyword='elementshapes' )

      if ( all ( elementshapes /= &
        [ 'lagrange         ', 'serendipity      ', &
          'extend           ', 'extend14         ', 'reduce           ', &
          'spectraltolinear ', 'quadratictolinear' ] ) ) then
        write(*,'(20(/a)/)') &
          'Error in mesh_convert: ', &
          ' heading parameter elementshape must be one of ', &
          '   lagrange ', &
          '   serendipity ', &
          '   extend ', &
          '   extend14 ', &
          '   reduce ', &
          '   spectraltolinear', &
          '   quadratictolinear'
        stop
      else if ( any ( elementshapes == [ 'lagrange   ', 'extend     ', &
        'extend14   ' ] ) ) then
        extend = .true.
      else if ( any ( elementshapes == [ 'serendipity      ', &
                    'reduce           ', 'quadratictolinear' ] ) ) then
        extend = .false.
      end if

      numkey = numkey + 1

    end if

!   convert coordinates

    if ( present(coordinates) ) then

!     test

      if ( any ( abs(coordinates) > mesh1%ndim ) ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' heading parameter coordinates is out of range'
        stop
      end if

      if ( size(coordinates) /= mesh1%ndim .or. any ( coordinates == 0 ) ) then

!       change in dimension or zero setting of coordinates

        do elgrp = 1, mesh1%nelgrp

          select case ( mesh1%element(elgrp)%elshape )
            case (1,2,32,101,105) ! line
              if ( size(coordinates) == 0 .or. size(coordinates) > 3 ) then
                write(*,'(2(/a)/)') &
                  'Error in mesh_convert: ', &
                  ' the space dimension of line elements should be 1, 2 or 3'
                stop
              end if
            case (3:10,30,33,34,102,104,106) ! triangles or quads
              if ( size(coordinates) <= 1 .or. size(coordinates) > 3 ) then
                write(*,'(3(/a)/)') &
                  'Error in mesh_convert: ', &
                  ' the space dimension of triangular or quadrilateral', &
                  ' elements should be 2 or 3'
                stop
              end if
            case default
              write(*,'(3(/a)/)') &
                'Error in mesh_convert: ', &
                ' the space dimension of 3D elements cannot be changed', &
                ' nor can the coordinates be set to zero'
              stop
          end select

        end do

      end if

      numkey = numkey + 1

    else if ( present(mirror) ) then

      write(*,'(2(/a)/)') &
        'Error in mesh_convert: ', &
        ' argument mirror can only used if argument coordinates is present.'
        stop

    end if

!   convert orientation of surfaces

    if ( present(orientation) ) then

!     test

      call check_blend ( mesh1, 'mesh_convert', keyword='orientation' )

      if ( .not. present(surfaces) ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' additional heading parameter surfaces is missing for orientation'
        stop
      end if

      if ( size(orientation) /= size(surfaces) ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' the size of orientation is different from the size of surfaces'
        stop
      end if

      if ( any ( abs(orientation) > 4 ) .or. any ( orientation == 0 ) ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' heading parameter orientation is out of range (lt -4, 0 or gt 4)'
        stop
      end if

      if ( any ( surfaces <= 0 ) .or. any( surfaces > mesh1%nsurfaces ) ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' heading parameter surfaces is out of range'
        stop
      end if

      numkey = numkey + 1

    end if

!   convert to user mesh

    if ( present(userelmesh) ) then
      call check_blend ( mesh1, 'mesh_convert', keyword='userelmesh' )
      numkey = numkey + 1
    end if

    if ( present(remove_double_nodes) ) then
      numkey = numkey + 1
    end if

    if ( present(remove_isolated_nodes) ) then
      numkey = numkey + 1
    end if

    if ( present(only_groups) ) then

      call check_blend ( mesh1, 'mesh_convert', keyword='only_groups' )

      if ( any ( only_groups <= 0 ) .or. &
           any ( only_groups > mesh1%nelgrp ) ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' heading argument only_groups is out of range'
        stop
      end if

      do grp = 1, size(only_groups)
        if ( count ( only_groups == only_groups(grp) ) > 1 ) then
          write(*,'(2(/a),i0/)') &
            'Error in mesh_convert: ', &
            ' only_groups contains identical group numbers, elgroup = ', &
             only_groups(grp)
          stop
        end if
      end do

      numkey = numkey + 1

    end if

    if ( present(from_curve) ) then

      call check_blend ( mesh1, 'mesh_convert', keyword='from_curve' )

      if ( from_curve <= 0 .or. &
           from_curve > mesh1%ncurves ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' heading argument from_curves is out of range'
        stop
      end if

      numkey = numkey + 1

    end if

    if ( present(from_surface) ) then

      call check_blend ( mesh1, 'mesh_convert', keyword='from_surface' )

      if ( from_surface <= 0 .or. &
           from_surface > mesh1%nsurfaces ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' heading argument from_surfaces is out of range'
        stop
      end if

      numkey = numkey + 1

    end if

    if ( present(from_volume) ) then

      call check_blend ( mesh1, 'mesh_convert', keyword='from_volume' )

      if ( from_volume <= 0 .or. &
           from_volume > mesh1%nvolumes ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' heading argument from_volumes is out of range'
        stop
      end if

      numkey = numkey + 1

    end if

    if ( present(blendmesh) ) then

      if ( any( [ mesh1%ndim, mesh1%nelgrp, mesh1%ncurves, &
                  mesh1%nsurfaces, mesh1%nvolumes ] /= &
                [ blendmesh%ndim, blendmesh%nelgrp, blendmesh%ncurves, &
                  blendmesh%nsurfaces, blendmesh%nvolumes ] ) ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' mesh1 and blendmesh are incompatible '
        stop
      end if

      if ( any( mesh1%grpnumel /= blendmesh%grpnumel ) ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' mesh1 and blendmesh have different number of elements '
        stop
      end if

      if ( any( mesh1%element(:)%globalshape /= &
                     blendmesh%element(:)%globalshape ) ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' mesh1 and blendmesh have different globalshapes '
        stop
      end if

      numkey = numkey + 1

    end if

    if ( present(blend) ) then

      if ( blend < 0 .or. &
           blend > mesh1%nblend ) then
        write(*,'(2(/a)/)') &
          'Error in mesh_convert: ', &
          ' heading argument blend is out of range'
        stop
      end if

      numkey = numkey + 1

    end if

!   Too many keywords?

    if ( numkey > 1 ) then
      write(*,'(/4(a/))') &
        'Error in mesh_convert: ', &
        ' heading parameters from the following list cannot be combined: ', &
        '   elementshapes, coordinates, orientation, userelmesh, ', &
        '   remove_double_nodes, remove_isolated_nodes, only_groups ', &
        '   from_curve, from_surface, from_volume, blendmesh, blend '
      stop
    end if

    lremove_double_nodes = &
             set_optional ( variable=remove_double_nodes, default=.false. )

    lremove_isolated_nodes = &
             set_optional ( variable=remove_isolated_nodes, default=.false. )

    if ( present(elementshapes) ) then

!     convert element shapes

      if ( elementshapes == 'spectraltolinear' ) then

!       convert high-order elements to linear elements

        call mesh_convert_spectraltolinear ( mesh1, mesh )

      else if ( elementshapes == 'quadratictolinear' ) then

!       convert quadratic elements to linear elements

        call mesh_convert_quadratictolinear ( mesh1, mesh )

      else if ( extend ) then

!       add nodes to element shapes

        call mesh_convert_add_elnodes ( mesh1, mesh, elementshapes )

      else

!       remove nodes from element shapes

        call mesh_convert_remove_elnodes ( mesh1, mesh, elementshapes )

      end if

    else if ( present(coordinates) ) then

!     convert coordinates

      call mesh_convert_coordinates ( mesh1, mesh, coordinates, mirror )

    else if ( present(orientation) ) then

!     convert orientation of surfaces

      call mesh_convert_orientation ( mesh1, mesh, orientation, surfaces )

    else if ( present(userelmesh) ) then

!     convert orientation of surfaces

      call mesh_convert_userelmesh ( mesh1, mesh, userelmesh )

    else if ( lremove_double_nodes ) then

!     remove multiple nodes

      call check_blend ( mesh1, 'mesh_convert', &
        keyword='remove_double_nodes=.true.', comment=&
        'Remove double nodes from meshes before creating the blended mesh.' )

      if ( .not. mesh1%nodblocks%filled ) then
        write(*,'(/3(a/))') &
          'Error in mesh_convert: ', &
          '  for remove_double_nodes = .true. the input mesh must, ', &
          '  contain nodblocks.'
        stop
      end if

      call mesh_convert_remove_double_nodes ( mesh1, mesh, excludepoints, &
        excludecurves, excludesurfaces, xmin, xmax )

    else if ( lremove_isolated_nodes ) then

!     remove isolated nodes

      call check_blend ( mesh1, 'mesh_convert', &
        keyword='remove_isolated_nodes=.true.', comment=&
        'Remove isolated nodes from meshes before creating the blended mesh.' )

      call mesh_convert_remove_isolated_nodes ( mesh1, mesh )

    else if ( present(only_groups) ) then

!     convert to a mesh with only the specified groups present

      call mesh_convert_only_groups ( mesh1, mesh, only_groups )

    else if ( present(from_curve) ) then

!     convert curve to a mesh with only a single group.

      call mesh_convert_from_geometry ( mesh1, mesh1%curves(from_curve), mesh )

    else if ( present(from_surface) ) then

!     convert surface to a mesh with only a single group.

      call mesh_convert_from_geometry ( mesh1, &
                                        mesh1%surfaces(from_surface), mesh )

    else if ( present(from_volume) ) then

!     convert volume to a mesh with only a single group.

      call mesh_convert_from_geometry ( mesh1, &
                                        mesh1%volumes(from_volume), mesh )

    else if ( present(blendmesh) ) then

!     add mesh from blendmesh as blend mesh

      call mesh_convert_add_blend ( mesh1, blendmesh, mesh )

    else if ( present(blend) ) then

!     reduce mesh to single shape mesh

      call mesh_convert_single_shape ( mesh1, blend, mesh )

    else

!     copy all elements

      call mesh_convert_copy ( mesh1, mesh )

    end if

  end subroutine mesh_convert


! Copy a mesh to a new basic mesh with identical properties

  subroutine mesh_convert_copy ( mesh1, mesh )

!   the mesh that needs to be copied
    type(mesh_t), intent(in) :: mesh1

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh


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

    allocate( mesh%points(MAXPOINTS) )

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

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    mesh%meshgen = .true.

  end subroutine mesh_convert_copy


! Copy a mesh to a new mesh where all serendipity elements have been
! converted to Langrangian elements or "standard" elements to "extended"
! elements by adding nodes to the element shapes

  subroutine mesh_convert_add_elnodes ( mesh1, mesh, elementshapes )

    use shapefunc_standard_m

!   the mesh that needs to be converted
    type(mesh_t), intent(in) :: mesh1

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh

!   Convert (some) element shapes to other element shapes:
!   Possibilities:
!     'lagrange' : convert serendipity elements (elshape=30,31) to
!                  Lagrangian interpolated elements (elshape=6,14).
!                  All other elements are just copied.
!     'extend' : convert "standard" elements to "extended" elements. The
!                following elshapes are converted:
!                    before   after
!                      3       10
!                      4        7
!                      5        9
!                     11       18
!                     12       16
!                     13       17
!                 All other elements are just copied.
!     'extend14' : convert "standard" element with elshape=12 (10-node
!                  tetrahedron) to "extended" element with elshape=15 (14-node
!                  tetrahedron).
!                  All other elements are just copied.
    character(len=*), intent(in) :: elementshapes


    logical :: serendip, extend, extend14
    integer :: elgrp, surface, numnod, elem, i, elshape
    integer :: nnr(mesh1%nnodes)! new node numbering for each old node in mesh1
    integer :: nnnodes ! number of nodes in new mesh
    integer, dimension(maxval(mesh1%element(:)%numnod)) :: nnod3


!   copy some data

    mesh%ndim = mesh1%ndim
    mesh%nelem = mesh1%nelem
    mesh%nelgrp = mesh1%nelgrp
    mesh%nblend = 0

!   internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology ( mesh%nelgrp ) )

    mesh%element(:)%ndim = mesh1%element(:)%ndim
    mesh%element(:)%globalshape = mesh1%element(:)%globalshape
    mesh%grpnumel(:) = mesh1%grpnumel(:)

!   loop over groups

    nnr = 0
    nnnodes = 0
    serendip = .false.
    extend = .false.
    extend14 = .false.

    do elgrp = 1, mesh%nelgrp

      elshape = mesh1%element(elgrp)%elshape

      if ( elshape == 30 .and. elementshapes == 'lagrange' ) then

!       8-node quad -> 9-node quad

        serendip = .true.

        call add_nodes_group ( newelshape=6, numnod=9, &
          po=[1,2,3,4,5,6,7,8], pc=9 )

      else if ( elshape == 31 .and. elementshapes == 'lagrange' ) then

!       20-node hex -> 27-node hex

        if ( .not. mesh1%meshparts ) then
          write(*,'(/3(a/))') &
            'Error in mesh_convert_add_elnodes: ', &
            ' For conversion of 20-node to 27-node hexahedron ', &
            ' meshparts of mesh1 needs to be filled.'
          stop
        end if

        serendip = .true.

        call add_nodes_group ( newelshape=14, numnod=27, &
          po=[1,2,3,4,6,7,8,9,10,12,16,18,19,20,21,22,24,25,26,27], &
          pf=[5,11,15,17,13,23], pc=14 )

      else if ( elshape == 3 .and. elementshapes == 'extend' ) then

!       3-node triangle -> 4-node triangle

        extend = .true.

        call add_nodes_group ( newelshape=10, numnod=4, po=[1,2,3], pc=4 )

      else if ( elshape == 4 .and. elementshapes == 'extend' ) then

!       6-node triangle -> 7-node triangle

        extend = .true.

        call add_nodes_group ( newelshape=7, numnod=7, po=[1,2,3,4,5,6], pc=7 )

      else if ( elshape == 5 .and. elementshapes == 'extend' ) then

!       4-node quad -> 5-node quad

        extend = .true.

        call add_nodes_group ( newelshape=9, numnod=5, po=[1,2,3,4], pc=5 )

      else if ( elshape == 11 .and. elementshapes == 'extend' ) then

!       4-node tetrahedron -> 5-node tetrahedron

        extend = .true.

        call add_nodes_group ( newelshape=18, numnod=5, po=[1,2,3,4], pc=5 )

      else if ( elshape == 12 .and. elementshapes == 'extend' ) then

!       10-node tetrahedron -> 15-node tetrahedron

        if ( .not. mesh1%meshparts ) then
          write(*,'(/3(a/))') &
            'Error in mesh_convert_add_elnodes: ', &
            ' For conversion of 10-node to 15-node tetrahedron ', &
            ' meshparts of mesh1 needs to be filled.'
          stop
        end if

        extend = .true.

        call add_nodes_group ( newelshape=16, numnod=15, &
          po=[1,2,3,4,5,6,7,8,9,10], pf=[11,12,13,14], pc=15 )

      else if ( elshape == 12 .and. elementshapes == 'extend14' ) then

!       10-node tetrahedron -> 14-node tetrahedron

        if ( .not. mesh1%meshparts ) then
          write(*,'(/3(a/))') &
            'Error in mesh_convert_add_elnodes: ', &
            ' For conversion of 10-node to 14-node tetrahedron ', &
            ' meshparts of mesh1 needs to be filled.'
          stop
        end if

        extend14 = .true.

        call add_nodes_group ( newelshape=15, numnod=14, &
          po=[1,2,3,4,5,6,7,8,9,10], pf=[11,12,13,14] )

      else if ( elshape == 13 .and. elementshapes == 'extend' ) then

!       8-node hexahedron -> 9-node hexahedron

        extend = .true.

        call add_nodes_group ( newelshape=17, numnod=9, &
                               po=[1,2,3,4,5,6,7,8], pc=9 )

      else

!       copy elements with new numbering

        numnod = mesh1%element(elgrp)%numnod

        mesh%element(elgrp)%elshape = mesh1%element(elgrp)%elshape
        mesh%element(elgrp)%p = mesh1%element(elgrp)%p
        mesh%element(elgrp)%numnod = numnod

        allocate( mesh%topology(elgrp)%a(numnod,mesh%grpnumel(elgrp)) )

!       loop over elements in this group

        do elem = 1, mesh%grpnumel(elgrp)

          nnod3(1:numnod) = nnr(mesh1%topology(elgrp)%a(:,elem))

          do i = 1, numnod
            if ( nnod3(i) == 0 ) then
  !            new node does not exist
               nnnodes = nnnodes + 1
               nnod3(i) = nnnodes
               nnr(mesh1%topology(elgrp)%a(i,elem)) = nnnodes
            end if
          end do

          mesh%topology(elgrp)%a(:,elem) = nnod3(1:numnod)

        end do

      end if

    end do

    if ( any( nnr == 0 ) ) then

      write(*,'(/a/a,i0,a/a/)') &
        'Warning in mesh_convert_add_elnodes: ', &
        ' There are ', count(nnr==0), &
        ' nodes not connected to the internal elements.' , &
        ' Continuing anyway.'

      do i = 1, mesh1%nnodes
        if ( nnr(i) == 0 ) then
          nnnodes = nnnodes + 1
          nnr(i) = nnnodes
        end if
      end do

    end if

    if ( .not. serendip .and. elementshapes == 'lagrange' ) then
      write(*,'(2(/a)/)') &
        'Warning in mesh_convert_add_elnodes: ', &
        ' No serendipity elements (8-node quad or 20-node hexahedron) found.'
    else if ( .not. extend .and. elementshapes == 'extend' ) then
      write(*,'(2(/a)/)') &
        'Warning in mesh_convert_add_elnodes: ', &
        ' No elements found that can be extended'
    else if ( .not. extend14 .and. elementshapes == 'extend14' ) then
      write(*,'(2(/a)/)') &
        'Warning in mesh_convert_add_elnodes: ', &
        ' No elements found that can be extended to 14 node tetrahedrons'
    end if

    mesh%nnodes = nnnodes

    mesh%nnodes_blend = [0,mesh%nnodes]

    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   copy old coordinates

    mesh%coor(nnr,:) = mesh1%coor

!   compute new coordinates

    do elgrp = 1, mesh%nelgrp

      elshape = mesh1%element(elgrp)%elshape

      if ( elshape == 30 .and. elementshapes == 'lagrange' ) then

!       8-node quad -> 9-node quad

        call add_coordinates ( oldnumnod=8,&
          xrnodc = reshape ( [ 0._dp, 0._dp ], [1,2] ), &
          pc=9, shape_element=shape_quad_serendipity2 )

      else if ( elshape == 31 .and. elementshapes == 'lagrange' ) then

!       20-node hex -> 27-node hex

        call add_coordinates ( oldnumnod=20,&
          xrnodf = reshape ( [  0._dp,  0._dp,  1._dp,  0._dp, -1._dp,  0._dp, &
                                0._dp, -1._dp,  0._dp,  1._dp,  0._dp,  0._dp, &
                               -1._dp,  0._dp,  0._dp,  0._dp,  0._dp,  1._dp  &
                             ], [6,3] ), &
          xrnodc = reshape ( [ 0._dp, 0._dp, 0._dp ], [1,3] ), &
          pf=[5,11,15,17,13,23], pc=14, &
          shape_element=shape_hexa_serendipity2 )

      else if ( elshape == 3 .and. elementshapes == 'extend' ) then

!       3-node triangle -> 4-node triangle

        call add_coordinates ( oldnumnod=3,&
          xrnodc = reshape ( [ 1._dp/3, 1._dp/3 ], [1,2] ), &
          pc=4, shape_element=shape_triangle_P1 )

      else if ( elshape == 4 .and. elementshapes == 'extend' ) then

!       6-node triangle -> 7-node triangle

        call add_coordinates ( oldnumnod=6,&
          xrnodc = reshape ( [ 1._dp/3, 1._dp/3 ], [1,2] ), &
          pc=7, shape_element=shape_triangle_P2 )

      else if ( elshape == 5 .and. elementshapes == 'extend' ) then

!       4-node quad -> 5-node quad

        call add_coordinates ( oldnumnod=4,&
          xrnodc = reshape ( [ 0._dp, 0._dp ], [1,2] ), &
          pc=5, shape_element=shape_quad_Q1 )

      else if ( elshape == 11 .and. elementshapes == 'extend' ) then

!       4-node tetrahedron -> 5-node tetrahedron

        call add_coordinates ( oldnumnod=4,&
          xrnodc = reshape ( [ 1._dp/4, 1._dp/4, 1._dp/4 ], [1,3] ), &
          pc=5, shape_element=shape_tetra_P1 )

     else if ( elshape == 12 .and. elementshapes == 'extend' ) then

!       10-node tetrahedron -> 15-node tetrahedron

        call add_coordinates ( oldnumnod=10,&
          xrnodf = reshape ( [ 1._dp/3, 1._dp/3, 1._dp/3, 0._dp, &
                               1._dp/3, 0._dp,   1._dp/3, 1._dp/3, &
                               0._dp,   1._dp/3, 1._dp/3, 1._dp/3 ], &
                               [4,3] ), &
          xrnodc = reshape ( [ 1._dp/4, 1._dp/4, 1._dp/4 ], [1,3] ), &
          pf=[11,12,13,14], pc=15, &
          shape_element=shape_tetra_P2 )

      else if ( elshape == 12 .and. elementshapes == 'extend14' ) then

!       10-node tetrahedron -> 14-node tetrahedron

        call add_coordinates ( oldnumnod=10,&
          xrnodf = reshape ( [ 1._dp/3, 1._dp/3, 1._dp/3, 0._dp, &
                               1._dp/3, 0._dp,   1._dp/3, 1._dp/3, &
                               0._dp,   1._dp/3, 1._dp/3, 1._dp/3 ],  &
                               [4,3] ), &
          pf=[11,12,13,14], &
          shape_element=shape_tetra_P2 )

      else if ( elshape == 13 .and. elementshapes == 'extend' ) then

!       8-node hexahedron -> 9-node hexahedron

        call add_coordinates ( oldnumnod=8,&
          xrnodc = reshape ( [ 0._dp, 0._dp, 0._dp ], [1,3] ), &
          pc=9, shape_element=shape_hexa_Q1 )

      end if

    end do

!   points (transfer to new numbering with nnr).

    mesh%npoints = mesh1%npoints

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1:mesh%npoints) = nnr(mesh1%points(1:mesh1%npoints))

!   curves (transfer to new numbering with nnr).

    mesh%ncurves = mesh1%ncurves

    allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )

    call copy ( mesh1%curves(:mesh1%ncurves), mesh%curves(:mesh%ncurves), nnr )

    mesh%curves(:mesh%ncurves)%n = 0  ! cannot be used for changing orientation

!   surfaces

    mesh%nsurfaces = mesh1%nsurfaces

    allocate( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )

!   this is only for the case that a surface is converted, otherwise there is
!   copy of this info twice
    mesh%surfaces(:mesh%nsurfaces)%ndim = mesh1%surfaces(:mesh1%nsurfaces)%ndim
    do surface = 1, mesh%nsurfaces
      call copy ( mesh1%surfaces(surface)%element, &
                  mesh%surfaces(surface)%element )
    end do
    mesh%surfaces(:mesh%nsurfaces)%nelem = &
                                    mesh1%surfaces(:mesh1%nsurfaces)%nelem

    do surface = 1, mesh%nsurfaces

      elshape = mesh1%surfaces(surface)%element%elshape

      if ( elshape == 30 .and. elementshapes == 'lagrange' ) then

!       surface element is 8-node quadrilateral: convert to 9-node

        call add_nodes_surface ( oldnumnod=8, newelshape=6, numnod=9, &
          xrnodc = reshape ( [ 0._dp, 0._dp ], [1,2] ), &
          po=[1,2,3,4,5,6,7,8], pc=9, shape_element=shape_quad_serendipity2 )

      else if ( elshape == 4 .and. &
        ( elementshapes == 'extend' .or.  elementshapes == 'extend14' ) ) then

!       surface element is 6-node triangle: convert to 7-node

        call add_nodes_surface ( oldnumnod=6, newelshape=7, numnod=7, &
          xrnodc = reshape ( [ 1._dp/3, 1._dp/3 ], [1,2] ), &
          po=[1,2,3,4,5,6], pc=7, shape_element=shape_triangle_P2 )

      else

!       copy existing surface (transfer to new numbering with nnr).

        call copy ( mesh1%surfaces(surface), mesh%surfaces(surface), nnr )

      end if

    end do

!   blend

    mesh%surfaces(:mesh%nsurfaces)%nblend = 0
    do surface = 1, mesh%nsurfaces
      allocate ( mesh%surfaces(surface)%element_blend(0) )
      mesh%surfaces(surface)%nnodes_blend = [0,mesh%surfaces(surface)%nnodes]
    end do

!   cannot be used for changing orientation
    mesh%surfaces(:mesh%nsurfaces)%n = 0
    mesh%surfaces(:mesh%nsurfaces)%m = 0

!   volumes

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    mesh%meshgen = .true.

  contains


!   inline routine for adding nodes to elements in a single group

    subroutine add_nodes_group ( newelshape, numnod, po, pf, pc )

      integer, intent(in) :: newelshape ! elshape of the new element
      integer, intent(in) :: numnod ! number of nodes in the new element

!     position in the new element of the old nodes (vertices, midside edges)
      integer, intent(in), dimension(:) :: po

!     position in the new element of the mid face nodes
      integer, intent(in), dimension(:), optional :: pf
!     position in the new element of the center node
      integer, intent(in), optional :: pc

      integer :: nnodnew(numnod) ! new node numbering in element
      integer :: face, elem, elgrpnr, elemnr, sidenr, i

      mesh%element(elgrp)%elshape = newelshape
      mesh%element(elgrp)%numnod = numnod

      allocate( mesh%topology(elgrp)%a(numnod,mesh%grpnumel(elgrp)) )

!     loop over elements in this group

      do elem = 1, mesh%grpnumel(elgrp)

!       vertices and midside edges
        nnodnew(po) = nnr( mesh1%topology(elgrp)%a(:,elem) )

        if ( present(pf) ) then

!         midside faces

          nnodnew(pf) = 0

          do face = 1, size(pf)

            elgrpnr = mesh1%sidelem(elgrp)%a(face,elem,1)
            elemnr = mesh1%sidelem(elgrp)%a(face,elem,2)
            sidenr = mesh1%sidelem(elgrp)%a(face,elem,3)

            if ( elgrpnr /= 0 ) then

!             connecting element

              if ( mesh%element(elgrpnr)%elshape /= newelshape ) then
                write(*,'(/2(a/))') &
                  'Error in mesh_convert_add_elnodes: ', &
                  ' Connecting element is not conforming.'
                stop
              end if

              if ( elgrpnr < elgrp .or. &
                   ( elgrpnr == elgrp .and. elemnr < elem ) ) then
!               element already with numbering
                nnodnew(pf(face)) = mesh%topology(elgrpnr)%a(pf(sidenr),elemnr)
                if ( nnodnew(pf(face)) == 0 ) &
                             stop 'internal error nnodnew(pf(face)) == 0 '
              end if

            end if

          end do

        end if

!       vertices and midside edges

        do i = 1, size(po)
          if ( nnodnew(po(i)) == 0 ) then
!           new node number does not exist
            nnnodes = nnnodes + 1
            nnodnew(po(i)) = nnnodes
            nnr(mesh1%topology(elgrp)%a(i,elem)) = nnnodes
          end if
        end do

        if ( present(pf) ) then

!         midside faces

          do face = 1, size(pf)
            if ( nnodnew(pf(face)) == 0 ) then
!             new node number does not exist
              nnnodes = nnnodes + 1
              nnodnew(pf(face)) = nnnodes
            end if
          end do

        end if

        if ( present(pc) ) then

!         center node

          nnnodes = nnnodes + 1
          nnodnew(pc) = nnnodes

          mesh%topology(elgrp)%a(:,elem) = nnodnew

        end if

      end do

    end subroutine add_nodes_group


!   inline routine for coordinates of the new nodes

    subroutine add_coordinates ( oldnumnod, xrnodf, xrnodc, pf, pc, &
      shape_element )

      integer, intent(in) :: oldnumnod ! number of nodes in the old element

!     reference coordinates of face nodes
      real(dp), intent(in), dimension(:,:), optional :: xrnodf

!     reference coordinates of center node
      real(dp), intent(in), dimension(:,:), optional :: xrnodc

!     position in the new element of the mid face nodes
      integer, intent(in), dimension(:), optional :: pf

!     position in the new element of the center node
      integer, intent(in), optional :: pc

!     shape function routine of the old element
      interface
        subroutine shape_element ( xr, phi, dphi )
          use kind_defs_m
          implicit none
          real(dp), intent(in), dimension(:,:) :: xr
          real(dp), intent(out), dimension(:,:) :: phi
          real(dp), intent(out), dimension(:,:,:), optional :: dphi
        end subroutine shape_element
      end interface


      integer :: elem, nodc
      integer, dimension(:), allocatable :: nodf
      real(dp) :: phic(1,oldnumnod), x(oldnumnod,mesh%ndim)
      real(dp), dimension(:,:), allocatable :: phif


!     mid face nodes

      if ( present(pf) ) then

        allocate ( nodf(size(pf)), phif(size(pf),oldnumnod) )

        call shape_element ( xrnodf, phif )

        do elem = 1, mesh%grpnumel(elgrp)

          x = mesh1%coor(mesh1%topology(elgrp)%a(:,elem),:) ! old coordinates

          nodf = mesh%topology(elgrp)%a(pf,elem) ! mid face nodes

          mesh%coor(nodf,:) = matmul ( phif, x )

        end do

        deallocate ( nodf, phif )

      end if

!     center node

      if ( present(pc) ) then

        call shape_element ( xrnodc, phic )

        do elem = 1, mesh%grpnumel(elgrp)

          x = mesh1%coor(mesh1%topology(elgrp)%a(:,elem),:) ! old coordinates

          nodc = mesh%topology(elgrp)%a(pc,elem) ! center node

          mesh%coor(nodc,:) = matmul ( phic(1,:), x )

        end do

      end if

    end subroutine add_coordinates


!   inline routine for adding a center node to the elements in a single surface

    subroutine add_nodes_surface ( oldnumnod, newelshape, numnod, xrnodc, &
      po, pc, shape_element )

      integer, intent(in) :: oldnumnod ! number of nodes in the old element
      integer, intent(in) :: newelshape ! elshape of the new element
      integer, intent(in) :: numnod ! number of nodes in the new element

!     reference coordinates of center node
      real(dp), intent(in), dimension(:,:), optional :: xrnodc

!     position in the new element of the old nodes (vertices, midside edges)
      integer, intent(in), dimension(:) :: po

!     position in the new element of the center node
      integer, intent(in), optional :: pc

!     shape function routine of the old element
      interface
        subroutine shape_element ( xr, phi, dphi )
          use kind_defs_m
          implicit none
          real(dp), intent(in), dimension(:,:) :: xr
          real(dp), intent(out), dimension(:,:) :: phi
          real(dp), intent(out), dimension(:,:,:), optional :: dphi
        end subroutine shape_element
      end interface


      integer :: elem, nelem, nnodes1, node, nodenr
      real(dp), dimension(:), allocatable :: distance2
      real(dp) :: phic(1,oldnumnod)
      real(dp) :: x(oldnumnod,3), xc(3)


!     old number of nodes on the surface
      nnodes1 = mesh1%surfaces(surface)%nnodes

      nelem = mesh1%surfaces(surface)%nelem

!     new number of nodes on the surface
      mesh%surfaces(surface)%nnodes = nnodes1 + nelem

      mesh%surfaces(surface)%element%elshape = newelshape
      mesh%surfaces(surface)%element%numnod = numnod

      allocate( mesh%surfaces(surface)%nodes(mesh%surfaces(surface)%nnodes) )

!     new node numbers for existing surface nodes
      mesh%surfaces(surface)%nodes(1:nnodes1) = &
                                  nnr(mesh1%surfaces(surface)%nodes)

!     center node must be added

      call shape_element ( xrnodc, phic )

      allocate( mesh%surfaces(surface)%topology(numnod,nelem,2) )

!     copy topology of existing surface nodes
      mesh%surfaces(surface)%topology(po,:,1) = &
                                   mesh1%surfaces(surface)%topology(:,:,1)
      do elem = 1, mesh%surfaces(surface)%nelem
        mesh%surfaces(surface)%topology(po,elem,2) = &
                          nnr(mesh1%surfaces(surface)%topology(po,elem,2))
      end do

      allocate ( distance2(mesh%nnodes) )

      do elem = 1, nelem

!        existing coordinates coordinates of the center node
         x = mesh1%coor(mesh1%surfaces(surface)%topology(:,elem,2),:)

         xc = matmul ( phic(1,:), x )

!        compute square distance
!        (here we could use sidenr for surfaces to avoid this!!)

         do node = 1, mesh%nnodes
           distance2(node) = &
               dot_product ( mesh%coor(node,:) - xc, mesh%coor(node,:) - xc )
         end do

         nodenr = minval ( minloc ( distance2 ) ) ! minimum distance

         mesh%surfaces(surface)%nodes ( nnodes1 + elem ) = nodenr

         mesh%surfaces(surface)%topology(pc,elem,1) = nnodes1 + elem
         mesh%surfaces(surface)%topology(pc,elem,2) = nodenr

      end do

      deallocate ( distance2 )

    end subroutine add_nodes_surface

  end subroutine mesh_convert_add_elnodes


! Copy a mesh to a new mesh where
!  1) all Lagrangian quadratic elements have been converted to
!     Serendipity elements or
!  2) "extended" elements have been reduced to "standard" elements or

  subroutine mesh_convert_remove_elnodes ( mesh1, mesh, elementshapes )

    use shapefunc_standard_m

!   the mesh that needs to be converted
    type(mesh_t), intent(in) :: mesh1

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh

!   Convert (some) element shapes to other element shapes:
!   Possibilities:
!     'serendipity' : convert Lagrange elements (elshape=6,14) to
!                     serendipity elements (elshape=30,31).
!                     All other elements are just copied.
!     'reduce' : convert "extended" elements to "standard" elements. The
!                following elshapes are converted:
!                    before   after
!                     10        3
!                      7        4
!                      9        5
!                     18       11
!                     15       12
!                     16       12
!                     17       13
!                 All other elements are just copied.
    character(len=*), intent(in) :: elementshapes

    logical :: lagrange, reduce

    integer :: elgrp, curve, surface, numnod, elem, i, j, elshape

    integer :: nnnodes ! number of nodes in new mesh
!   nnr: new node numbering for each old node in mesh1
    integer, allocatable, dimension(:) :: nnr
    integer, allocatable, dimension(:) :: zeropoints
    integer, dimension(maxval(mesh1%element(:)%numnod)) :: nnod3

!   copy some data

    mesh%ndim = mesh1%ndim
    mesh%nelem = mesh1%nelem
    mesh%nelgrp = mesh1%nelgrp
    mesh%nblend = 0

!   internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology ( mesh%nelgrp ) )

    mesh%element(:)%ndim = mesh1%element(:)%ndim
    mesh%element(:)%globalshape = mesh1%element(:)%globalshape
    mesh%grpnumel(:) = mesh1%grpnumel(:)

!   loop over groups

    allocate ( nnr(mesh1%nnodes) )

    nnr = 0
    nnnodes = 0
    lagrange = .false.
    reduce = .false.

    do elgrp = 1, mesh%nelgrp

      elshape = mesh1%element(elgrp)%elshape

      if ( elshape == 6 .and. elementshapes == 'serendipity' ) then

!       9-node quad -> 8-node quad

        lagrange = .true.

        call remove_nodes_group ( newelshape=30, numnod=8, &
          pn=[1,2,3,4,5,6,7,8] )

      else if ( elshape == 14 .and. elementshapes == 'serendipity' ) then

!       27-node hex -> 20-node hex

        lagrange = .true.

        call remove_nodes_group ( newelshape=31, numnod=20, &
          pn=[1,2,3,4,6,7,8,9,10,12,16,18,19,20,21,22,24,25,26,27] )

      else if ( elshape == 10 .and. elementshapes == 'reduce' ) then

!       4-node triangle -> 3-node triangle

        reduce = .true.

        call remove_nodes_group ( newelshape=3, numnod=3, pn=[1,2,3] )

      else if ( elshape == 7 .and. elementshapes == 'reduce' ) then

!       7-node triangle -> 6-node triangle

        reduce = .true.

        call remove_nodes_group ( newelshape=4, numnod=6, pn=[1,2,3,4,5,6] )

      else if ( elshape == 9 .and. elementshapes == 'reduce' ) then

!       5-node quad -> 4-node quad

        reduce = .true.

        call remove_nodes_group ( newelshape=5, numnod=4, pn=[1,2,3,4] )

      else if ( elshape == 18 .and. elementshapes == 'reduce' ) then

!       5-node tetrahedron -> 4-node tetrahedron

        reduce = .true.

        call remove_nodes_group ( newelshape=11, numnod=4, pn=[1,2,3,4] )

      else if ( elshape == 15 .and. elementshapes == 'reduce' ) then

!       14-node tetrahedron -> 10-node tetrahedron

        reduce = .true.

        call remove_nodes_group ( newelshape=12, numnod=10, &
          pn=[1,2,3,4,5,6,7,8,9,10] )

      else if ( elshape == 16 .and. elementshapes == 'reduce' ) then

!       15-node tetrahedron -> 10-node tetrahedron

        reduce = .true.

        call remove_nodes_group ( newelshape=12, numnod=10, &
          pn=[1,2,3,4,5,6,7,8,9,10] )

      else if ( elshape == 17 .and. elementshapes == 'reduce' ) then

!       9-node hexahedron -> 8-node hexahedron

        reduce = .true.

        call remove_nodes_group ( newelshape=13, numnod=8, &
          pn=[1,2,3,4,5,6,7,8] )

      else

!       copy elements with new numbering

        numnod = mesh1%element(elgrp)%numnod

        mesh%element(elgrp)%elshape = mesh1%element(elgrp)%elshape
        mesh%element(elgrp)%p = mesh1%element(elgrp)%p
        mesh%element(elgrp)%numnod = numnod

        allocate( mesh%topology(elgrp)%a(numnod,mesh%grpnumel(elgrp)) )

!       loop over elements in this group

        do elem = 1, mesh%grpnumel(elgrp)

          nnod3(1:numnod) = nnr(mesh1%topology(elgrp)%a(:,elem))

          do i = 1, numnod
            if ( nnod3(i) == 0 ) then
  !            new node does not exist
               nnnodes = nnnodes + 1
               nnod3(i) = nnnodes
               nnr(mesh1%topology(elgrp)%a(i,elem)) = nnnodes
            end if
          end do

          mesh%topology(elgrp)%a(:,elem) = nnod3(1:numnod)

        end do

      end if

    end do

    if ( .not. lagrange .and. elementshapes == 'serendipity' ) then
      write(*,'(3(/a)/)') &
        'Warning in mesh_convert_remove_elnodes: ', &
        ' No quadratic Lagrange elements found ', &
        '(9-node quad or 27-node hexahedron).'
    else if ( .not. reduce .and. elementshapes == 'reduce' ) then
      write(*,'(2(/a)/)') &
        'Warning in mesh_convert_remove_elnodes: ', &
        ' No elements found that can be reduced'
    end if

    mesh%nnodes = nnnodes

    mesh%nnodes_blend = [0,mesh%nnodes]

    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   copy old coordinates

    do i = 1, mesh1%nnodes
      if ( nnr(i) > 0 ) then
        mesh%coor(nnr(i),:) = mesh1%coor(i,:)
      end if
    end do

!   points (transfer to new numbering with nnr).

    mesh%npoints = mesh1%npoints

    allocate( mesh%points(max(mesh%npoints,MAXPOINTS)) )

    mesh%points(1:mesh%npoints) = nnr(mesh1%points(1:mesh1%npoints))

    if ( any ( mesh%points(1:mesh%npoints) == 0 ) ) then
      write(*,'(/3(a/),a)', advance='no') &
        'Warning in mesh_convert_remove_elnodes: ', &
        ' The converted mesh contains user points that are not connected ', &
        ' to nodal points. Do not try to refer to these points in any way.', &
        ' The point numbers are: '
      allocate ( zeropoints(count(mesh%points(1:mesh%npoints) == 0)) )
      j = 0
      do i = 1, mesh%npoints
        if ( mesh%points(i) /= 0 ) cycle
        j = j + 1
        zeropoints(j) = i
      end do
      write (*,*) zeropoints
      write (*,*)
      deallocate ( zeropoints )
    end if

!   curves (transfer to new numbering with nnr).

    mesh%ncurves = mesh1%ncurves

    allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )

    call copy ( mesh1%curves(:mesh1%ncurves), mesh%curves(:mesh%ncurves), nnr )

    mesh%curves(:mesh%ncurves)%n = 0  ! cannot be used for changing orientation

    do curve = 1, mesh%ncurves

      if ( any ( mesh%curves(curve)%nodes == 0 ) ) then
        write(*,'(/3(a/),a,i0)') 'Warning in mesh_convert_remove_elnodes: ', &
          ' A converted curve contains points that are not connected ', &
          ' to nodal points. Do not use this curve in any way.', &
          ' The curve number is: ', curve

      end if

    end do

!   surfaces

    mesh%nsurfaces = mesh1%nsurfaces

    allocate( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )

!   this is only for the case that a surface is converted, otherwise there is
!   copy of this info twice
    do surface = 1, mesh%nsurfaces
      call copy ( mesh1%surfaces(surface)%element, &
                  mesh%surfaces(surface)%element )
    end do
    mesh%surfaces(:mesh%nsurfaces)%ndim = mesh1%surfaces(:mesh1%nsurfaces)%ndim
    mesh%surfaces(:mesh%nsurfaces)%nelem = &
                               mesh1%surfaces(:mesh1%nsurfaces)%nelem

    do surface = 1, mesh%nsurfaces

      elshape = mesh1%surfaces(surface)%element%elshape

      if ( elshape == 6 .and. elementshapes == 'serendipity' ) then

!       surface element is 9-node quadrilateral: convert to 8-node

        call remove_nodes_geometry ( mesh1%surfaces(surface), &
          mesh%surfaces(surface), newelshape=30, numnod=8, &
          pn=[1,2,3,4,5,6,7,8], nnr=nnr )

      else if ( elshape == 7 .and. elementshapes == 'reduce' ) then

!       surface element is 7-node triangle: convert to 6-node

        call remove_nodes_geometry ( mesh1%surfaces(surface), &
          mesh%surfaces(surface), newelshape=4, numnod=6, &
          pn=[1,2,3,4,5,6], nnr=nnr )

      else

!       copy existing surface (transfer to new numbering with nnr).

        call copy ( mesh1%surfaces(:mesh1%nsurfaces), &
                    mesh%surfaces(:mesh%nsurfaces), nnr )

      end if

    end do

!   blend

    mesh%surfaces(:mesh%nsurfaces)%nblend = 0
    do surface = 1, mesh%nsurfaces
      allocate ( mesh%surfaces(surface)%element_blend(0) )
      mesh%surfaces(surface)%nnodes_blend = [0,mesh%surfaces(surface)%nnodes]
    end do

!   cannot be used for changing orientation
    mesh%surfaces(:mesh%nsurfaces)%n = 0
    mesh%surfaces(:mesh%nsurfaces)%m = 0

!   volumes

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    mesh%meshgen = .true.

    deallocate ( nnr )

  contains

!   inline routine for removing nodes from elements in a single group

    subroutine remove_nodes_group ( newelshape, numnod, pn )

      integer, intent(in) :: newelshape ! elshape of the new element
      integer, intent(in) :: numnod ! number of nodes in the new element

!     position in the old element of the new nodes (vertices, midside edges)
      integer, intent(in), dimension(:) :: pn

      integer :: nnodnew(numnod) ! new node numbering in element
      integer :: elem, i


      mesh%element(elgrp)%elshape = newelshape
      mesh%element(elgrp)%numnod = numnod

      allocate( mesh%topology(elgrp)%a(numnod,mesh%grpnumel(elgrp)) )

!     loop over elements in this group

      do elem = 1, mesh%grpnumel(elgrp)

!       vertices and midside edges
        nnodnew = nnr( mesh1%topology(elgrp)%a(pn,elem) )

        do i = 1, numnod
          if ( nnodnew(i) == 0 ) then
!           new node number does not exist
            nnnodes = nnnodes + 1
            nnodnew(i) = nnnodes
            nnr(mesh1%topology(elgrp)%a(pn(i),elem)) = nnnodes
          end if
        end do

        mesh%topology(elgrp)%a(:,elem) = nnodnew

      end do

    end subroutine remove_nodes_group

  end subroutine mesh_convert_remove_elnodes


! Copy a mesh to a new mesh where "quadratic" elements have been lowered
! to linear elements by removing all nodes, except the vertices.
! The following elshapes are converted:
!          before   after
!            2        1
!            4        3
!            6        5
!            7        3
!           12       11
!           14       13
!           15       11
!           30        5
!           31       13
!           42       41
!           43       41
!           52       51
!           53       51
! All elements in the groups must be of "quadratic" shape.

  subroutine mesh_convert_quadratictolinear ( mesh1, mesh )

    use shapefunc_standard_m

!   the mesh that needs to be converted
    type(mesh_t), intent(in) :: mesh1

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh

    integer :: elgrp, curve, surface, i, j, elshape, elem

    integer :: nnnodes ! number of nodes in new mesh
!   nnr: new node numbering for each old node in mesh1
    integer, allocatable, dimension(:) :: nnr
    integer, allocatable, dimension(:) :: zeropoints

!   copy some data

    mesh%ndim = mesh1%ndim
    mesh%nelem = mesh1%nelem
    mesh%nelgrp = mesh1%nelgrp
    mesh%nblend = 0

!   internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology ( mesh%nelgrp ) )

    mesh%element(:)%ndim = mesh1%element(:)%ndim
    mesh%element(:)%globalshape = mesh1%element(:)%globalshape
    mesh%grpnumel(:) = mesh1%grpnumel(:)

!   loop over groups

    allocate ( nnr(mesh1%nnodes) )

    nnr = 0
    nnnodes = 0

    do elgrp = 1, mesh%nelgrp

      select case ( mesh1%element(elgrp)%elshape )

      case(2) ! 3-node line -> 2 node line

        call remove_nodes_group ( newelshape=1, numnod=2, pn=[1,3] )

      case(4,7) ! 6 or 7-node triangle -> 3-node triangle

        call remove_nodes_group ( newelshape=3, numnod=3, pn=[1,3,5] )

      case(6,30) ! 8 or 9-node quad -> 4-node quad

        call remove_nodes_group ( newelshape=5, numnod=4, pn=[1,3,5,7] )

      case(12,15) ! 10 or 14-node tetrahedron -> 4-node tetrahedron

        call remove_nodes_group ( newelshape=11, numnod=4, pn=[1,3,5,10] )

      case(14) ! 27-node hexahedron -> 8-node hexahedron

        call remove_nodes_group ( newelshape=13, numnod=8, &
          pn=[1,3,9,7,19,21,27,25] )

      case(31) ! 20-node hexahedron -> 8-node hexahedron

        call remove_nodes_group ( newelshape=13, numnod=8, &
          pn=[1,3,8,6,13,15,20,18] )

      case(42) ! 15-node prism -> 6-node prism

        call remove_nodes_group ( newelshape=41, numnod=6, pn=[1,3,5,10,12,14] )

      case(43) ! 18-node prism -> 6-node prism

        call remove_nodes_group ( newelshape=41, numnod=6, pn=[1,3,5,13,15,17] )

      case(52) ! 13-node pyramid -> 5-node pyramid

        call remove_nodes_group ( newelshape=51, numnod=5, pn=[1,3,5,7,13] )

      case(53) ! 14-node pyramid -> 5-node pyramid

        call remove_nodes_group ( newelshape=51, numnod=5, pn=[1,3,5,7,14] )

      case default

        write(*,'(2(a/a,i0/))') &
          'Error in mesh_convert_quadratictolinear: ', &
          ' Element group ', elgrp, &
          ' does not have a quadratic element shape.', &
          ' elshape = ', mesh1%element(elgrp)%elshape
        stop

      end select

    end do

!   new numbering (keep original sequence)

    nnnodes = 0

    do i = 1, mesh1%nnodes
      if ( nnr(i) == 0 ) cycle
      nnnodes = nnnodes + 1
      nnr(i) = nnnodes
    end do

!   topology with new numbering

    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)
        mesh%topology(elgrp)%a(:,elem) = nnr(mesh%topology(elgrp)%a(:,elem))
      end do
    end do

    mesh%nnodes = nnnodes

    mesh%nnodes_blend = [0,mesh%nnodes]

    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   copy old coordinates

    do i = 1, mesh1%nnodes
      if ( nnr(i) > 0 ) then
        mesh%coor(nnr(i),:) = mesh1%coor(i,:)
      end if
    end do

!   points (transfer to new numbering with nnr).

    mesh%npoints = mesh1%npoints

    allocate( mesh%points(max(mesh%npoints,MAXPOINTS)) )

    mesh%points(1:mesh%npoints) = nnr(mesh1%points(1:mesh1%npoints))

    if ( any ( mesh%points(1:mesh%npoints) == 0 ) ) then
      write(*,'(/3(a/),a)', advance='no') &
        'Warning in mesh_convert_quadratictolinear: ', &
        ' The converted mesh contains user points that are not connected ', &
        ' to nodal points. Do not try to refer to these points in any way.', &
        ' The point numbers are: '
      allocate ( zeropoints(count(mesh%points(1:mesh%npoints) == 0)) )
      j = 0
      do i = 1, mesh%npoints
        if ( mesh%points(i) /= 0 ) cycle
        j = j + 1
        zeropoints(j) = i
      end do
      write (*,*) zeropoints
      write (*,*)
      deallocate ( zeropoints )
    end if

!   curves

    mesh%ncurves = mesh1%ncurves

    allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )

    do curve = 1, mesh%ncurves
      call copy ( mesh1%curves(curve)%element, mesh%curves(curve)%element )
    end do
    mesh%curves(:mesh%ncurves)%ndim = mesh1%curves(:mesh1%ncurves)%ndim
    mesh%curves(:mesh%ncurves)%nelem = &
                               mesh1%curves(:mesh1%ncurves)%nelem

    do curve = 1, mesh%ncurves

      elshape = mesh1%curves(curve)%element%elshape

      select case ( elshape )

      case(2)

!       curve element is 3-node line: convert to 2-node

        call remove_nodes_geometry ( mesh1%curves(curve), &
          mesh%curves(curve), newelshape=1, numnod=2, pn=[1,3], nnr=nnr )

      case default

        write(*,'(2(a/a,i0/))') &
          'Error in mesh_convert_quadratictolinear: ', &
          ' Curve ', curve, &
          ' does not have a quadratic element shape.', &
          ' elshape = ', elshape
        stop

      end select

    end do

!   blend

    mesh%curves(:mesh%ncurves)%nblend = 0
    do curve = 1, mesh%ncurves
      allocate ( mesh%curves(curve)%element_blend(0) )
      mesh%curves(curve)%nnodes_blend = [0,mesh%curves(curve)%nnodes]
    end do

!   cannot be used for changing orientation
    mesh%curves(:mesh%ncurves)%n = 0

!   surfaces

    mesh%nsurfaces = mesh1%nsurfaces

    allocate( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )

    do surface = 1, mesh%nsurfaces
      call copy ( mesh1%surfaces(surface)%element, &
                  mesh%surfaces(surface)%element )
    end do
    mesh%surfaces(:mesh%nsurfaces)%ndim = mesh1%surfaces(:mesh1%nsurfaces)%ndim
    mesh%surfaces(:mesh%nsurfaces)%nelem = &
                               mesh1%surfaces(:mesh1%nsurfaces)%nelem

    do surface = 1, mesh%nsurfaces

      elshape = mesh1%surfaces(surface)%element%elshape

      select case ( elshape )

      case(6,30)

!       surface element is 8 or 9-node quadrilateral: convert to 4-node

        call remove_nodes_geometry ( mesh1%surfaces(surface), &
          mesh%surfaces(surface), newelshape=5, numnod=4, pn=[1,3,5,7], &
          nnr=nnr )

      case(4,7)

!       surface element is 6 or 7-node triangle: convert to 3-node

        call remove_nodes_geometry ( mesh1%surfaces(surface), &
          mesh%surfaces(surface), newelshape=3, numnod=3, pn=[1,3,5], &
          nnr=nnr )

      case default

        write(*,'(2(a/a,i0/))') &
          'Error in mesh_convert_quadratictolinear: ', &
          ' Surface ', surface, &
          ' does not have a quadratic element shape.', &
          ' elshape = ', elshape
        stop

      end select

    end do

!   blend

    mesh%surfaces(:mesh%nsurfaces)%nblend = 0
    do surface = 1, mesh%nsurfaces
      allocate ( mesh%surfaces(surface)%element_blend(0) )
      mesh%surfaces(surface)%nnodes_blend = [0,mesh%surfaces(surface)%nnodes]
    end do

!   cannot be used for changing orientation
    mesh%surfaces(:mesh%nsurfaces)%n = 0
    mesh%surfaces(:mesh%nsurfaces)%m = 0

!   volumes

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    mesh%meshgen = .true.

    deallocate ( nnr )

  contains

!   inline routine for removing nodes from elements in a single group

    subroutine remove_nodes_group ( newelshape, numnod, pn )

      integer, intent(in) :: newelshape ! elshape of the new element
      integer, intent(in) :: numnod ! number of nodes in the new element

!     position in the old element of the new nodes (vertices, midside edges)
      integer, intent(in), dimension(:) :: pn

      integer :: elem

      mesh%element(elgrp)%elshape = newelshape
      mesh%element(elgrp)%numnod = numnod

      allocate( mesh%topology(elgrp)%a(numnod,mesh%grpnumel(elgrp)) )

!     loop over elements in this group

      do elem = 1, mesh%grpnumel(elgrp)

!       nodes to keep
        mesh%topology(elgrp)%a(:,elem) = mesh1%topology(elgrp)%a(pn,elem)

!       tag nodes
        nnr( mesh1%topology(elgrp)%a(pn,elem) ) = 1

      end do

    end subroutine remove_nodes_group

  end subroutine mesh_convert_quadratictolinear


! routine for removing a nodes from the elements in a single geometry

  subroutine remove_nodes_geometry ( geometry1, geometry, newelshape, &
    numnod, pn, nnr )

    type(geometry_t), intent(in) :: geometry1 ! geometry input
    type(geometry_t), intent(inout) :: geometry ! geometry with nodes removed

    integer, intent(in) :: newelshape ! elshape of the new element
    integer, intent(in) :: numnod ! number of nodes in the new element

!   position in the new element of the old nodes (vertices, midside edges)
    integer, intent(in), dimension(:) :: pn

!   new node numbering, with
!      nnr(i) = 0, node i has been removed
!      nnr(i) /= 0, node i has new number nnr(i)
    integer, intent(in), dimension(:) :: nnr

    integer :: elem, i
    integer :: nelem, nnodes1, node, nnodes
!   work array for storing local new numbering
    integer, allocatable, dimension(:) :: work

!   old number of nodes on the geometry
    nnodes1 = geometry1%nnodes
    nelem = geometry1%nelem

    nnodes = count ( nnr(geometry1%nodes) /= 0 ) ! number of nodes on geometry
    geometry%nnodes = nnodes
    geometry%element%elshape = newelshape
    geometry%element%numnod = numnod

    allocate ( work(nnodes1) )
    allocate ( geometry%nodes(nnodes) )

!   new node numbers for existing geometry nodes
    work = 0
    node = 0
    do i = 1, nnodes1
      if ( nnr(geometry1%nodes(i)) /= 0 ) then
         node = node + 1
         geometry%nodes(node) = nnr(geometry1%nodes(i))
         work(i) = node
      end if
    end do

    allocate( geometry%topology(numnod,nelem,2) )

!   copy topology of existing geometry nodes
    do elem = 1, geometry%nelem
      geometry%topology(:,elem,1) = work(geometry1%topology(pn,elem,1))
      geometry%topology(:,elem,2) = nnr(geometry1%topology(pn,elem,2))
    end do

  end subroutine remove_nodes_geometry


! Convert mesh1 to mesh having different space dimension and/or coordinates

  subroutine mesh_convert_coordinates ( mesh1, mesh, coordinates, mirror )

!   the mesh that needs to be converted
    type(mesh_t), intent(in) :: mesh1

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh

!   convert coordinates of the mesh
    integer, dimension(:), intent(in) :: coordinates

!   if present and .true., the tolopogy is mirrored
!   default=.false.
    logical, intent(in), optional :: mirror

    logical :: lmirror
    integer :: ic, elgrp, elem, j, m, elnumnod, s
    integer, dimension(:), allocatable :: mirror_index


    lmirror = set_optional ( variable=mirror, default=.false. )

!   copy mesh

    call mesh_convert_copy ( mesh1, mesh )

!   remove coordinates

    deallocate ( mesh%coor )

!   new coordinates

    mesh%ndim = size(coordinates)
    mesh%element(1:mesh%nelgrp)%ndim = mesh%ndim
    do j = 1, mesh%nblend
      mesh%element_blend(1:mesh%nelgrp,j)%ndim = mesh%ndim
    end do
    allocate ( mesh%coor(mesh%nnodes,mesh%ndim) )
    mesh%curves(1:mesh%ncurves)%ndim = mesh%ndim
    mesh%curves(1:mesh%ncurves)%element%ndim = mesh%ndim
    mesh%surfaces(1:mesh%nsurfaces)%ndim = mesh%ndim
    mesh%surfaces(1:mesh%nsurfaces)%element%ndim = mesh%ndim

    do ic = 1, mesh%ndim

      if ( coordinates(ic) /= 0 ) then
!       copy coordinate in direction ic with a possible change in sign
        mesh%coor(:,ic) = &
                mesh1%coor(:,abs(coordinates(ic))) * sign ( 1, coordinates(ic) )
      else
!       set coordinate in direction ic to zero
        mesh%coor(:,ic) = 0
      end if

    end do

    if ( lmirror ) then

!     mirror element topology

      do elgrp = 1, mesh%nelgrp

        elnumnod = mesh%element(elgrp)%numnod

        allocate ( mirror_index(elnumnod) )

        call fill_mirror_index ( mesh%element(elgrp)%elshape, mirror_index )

        do elem = 1, mesh%grpnumel(elgrp)
          mesh%topology(elgrp)%a(1:elnumnod,elem) = &
                             mesh%topology(elgrp)%a(mirror_index,elem)
        end do

        deallocate ( mirror_index )

!       blended meshes

        s = elnumnod

        do m = 1, mesh%nblend

          elnumnod = mesh%element_blend(elgrp,m)%numnod

          allocate ( mirror_index(elnumnod) )

          call fill_mirror_index ( mesh%element_blend(elgrp,m)%elshape, &
            mirror_index )

          do elem = 1, mesh%grpnumel(elgrp)
            mesh%topology(elgrp)%a(s+1:s+elnumnod,elem) = &
                               mesh%topology(elgrp)%a(mirror_index,elem)
          end do

          deallocate ( mirror_index )

          s = s + elnumnod

        end do

      end do

    end if

  contains

!   set index array for reordering topology for mirroring
!   old_node = mirror_index ( tfem_node )

    subroutine fill_mirror_index ( elshape, mirror_index )

      integer, intent(in) :: elshape
      integer, dimension(:), intent(out) :: mirror_index

      select case ( elshape )
        case (1) ! 2 node line
          mirror_index = [ 2, 1 ]
        case (2,32) ! 3 node line
          mirror_index = [ 3, 2, 1 ]
        case (3) ! 3 node triangle
          mirror_index = [ 2, 1, 3 ]
        case (4,33) ! 6 node triangle
          mirror_index = [ 3, 2, 1, 6, 5, 4 ]
        case (5) ! 4 node quadrilateral
          mirror_index = [ 2, 1, 4, 3 ]
        case (6,34) ! 9 node quadrilateral
          mirror_index = [ 3, 2, 1, 8, 7, 6, 5, 4, 9 ]
        case (7) ! 7 node triangle
          mirror_index = [ 3, 2, 1, 6, 5, 4, 7 ]
        case (9) ! 5 node quadrilateral
          mirror_index = [ 2, 1, 4, 3, 5 ]
        case (10) ! 4 node triangle
          mirror_index = [ 2, 1, 3, 4 ]
        case (11) ! 4 node tetrahedron
          mirror_index = [ 2, 1, 3, 4 ]
        case (12,35) ! 10 node tetrahedron
          mirror_index = [ 3, 2, 1, 6, 5, 4, 8, 7, 9, 10 ]
        case (13) ! 8 node hexahedron
          mirror_index = [ 2, 1, 4, 3, 6, 5, 8, 7 ]
        case (14,36) ! 27 node hexahedron
          mirror_index = [ 3, 2, 1, 6, 5, 4, 9, 8, 7, &
                           12, 11, 10, 15, 14, 13, 18, 17, 16, &
                           21, 20, 19, 24, 23, 22, 27, 26, 25 ]
        case (15) ! 14 node tetrahedron
          mirror_index = [ 3, 2, 1, 6, 5, 4, 8, 7, 9, 10, 11, 12, 14, 13 ]
        case (16) ! 15 node tetrahedron
          mirror_index = [ 3, 2, 1, 6, 5, 4, 8, 7, 9, 10, 11, 12, 14, 13, 15 ]
        case (17) ! 9 node hexahedron
          mirror_index = [ 2, 1, 4, 3, 6, 5, 8, 7, 9 ]
        case (18) ! 5 node tetrahedron
          mirror_index = [ 2, 1, 3, 4, 5 ]
        case (30) ! 8 node quadrilateral
          mirror_index = [ 3, 2, 1, 8, 7, 6, 5, 4 ]
        case (31) ! 20 node hexahedron
          mirror_index = [ 3, 2, 1, 5, 4, 8, 7, 6, &
                           10, 9, 12, 11, &
                           15, 14, 13, 17, 16, 20, 19, 18 ]
        case (41) ! 6 node prism
          mirror_index = [ 2, 1, 3, 5, 4, 6 ]
        case (42) ! 15 node prism
          mirror_index = [ 3, 2, 1, 6, 5, 4, 8, 7, 9, 12, 11, 10, 15, 14, 13 ]
        case (43) ! 18 node prism
          mirror_index = [ 3, 2, 1, 6, 5, 4, 9, 8, 7, 12, &
                           11, 10, 15, 14, 13, 18, 17, 16 ]
        case (51) ! 5 node pyramid
          mirror_index = [ 2, 1, 4, 3, 4, 5 ]
        case (52) ! 13 node pyramid
          mirror_index = [ 3, 2, 1, 8, 7, 6, 5, 4, 10, 9, 12, 11, 13 ]
        case (53) ! 14 node pyramid
          mirror_index = [ 3, 2, 1, 8, 7, 6, 5, 4, 9, 11, 10, 13, 12, 14 ]
        case default
          write(*,'(/a/a,i0/)') &
            'Error in mesh_convert_coordinates: ', &
            ' mirror not implemented for elshape = ', elshape
          stop
      end select

    end subroutine fill_mirror_index

  end subroutine mesh_convert_coordinates


! Convert mesh1 to mesh having different space orientation for the surfaces

  subroutine mesh_convert_orientation ( mesh1, mesh, orientation, surfaces )

!   the mesh that needs to be converted
    type(mesh_t), intent(in) :: mesh1

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh

!   Change the orientation of the specified surfaces.
    integer, dimension(:), intent(in) :: orientation, surfaces


    integer :: i, ori, flp, srf


!   copy mesh

    call mesh_convert_copy ( mesh1, mesh )

    do i = 1, size(surfaces)

      ori = abs(orientation(i))
      flp = sign ( 1, orientation(i) )

      srf = surfaces(i)

      if ( mesh%surfaces(srf)%n == 0 .or. mesh%surfaces(srf)%m == 0 ) then
        write(*,'(3(/a)/a,i0/)') &
          'Error in mesh_convert_orientation: ', &
          ' No structured mesh info available (m, n). ', &
          ' Mesh not made by hexahedron or info lost in conversion.', &
          ' surface = ', srf
        stop
      end if

      select case ( mesh%surfaces(srf)%element%elshape )

        case (5) ! 4 node quad

          call mesh_convert_orientation_quad4 ( mesh%surfaces(srf), ori, flp )

        case (6) ! 9 node quad

          call mesh_convert_orientation_quad9 ( mesh%surfaces(srf), ori, flp )

        case default

          write(*,'(2(/a)/a,i0/)') &
            'Error in mesh_convert_orientation: ', &
            ' only 4-node and 9-node quads have been implemented ', &
            ' surface = ', srf
          stop

      end select

    end do

  end subroutine mesh_convert_orientation


! convert orientation of surface with quad4 elements

  subroutine mesh_convert_orientation_quad4 ( surface, orientation, flip )

    type(geometry_t), intent(inout) :: surface

    integer, intent(in) :: orientation, flip


    integer, parameter :: inpelms = 4 ! number of nodal points per element on
                                      ! surfaces

    integer :: elsnodes(inpelms), i, j, e, nn1, nn2, m, n, elem
    integer, allocatable, dimension(:) :: newnodes

    allocate ( newnodes(surface%nnodes) )

    n = surface%n
    m = surface%m

    nn1 = n + 1
    nn2 = m + 1

    if ( flip == 1 ) then

!     usual direction

      select case ( orientation )

      case(1)

        return  ! nothing to do

      case(2)

        newnodes = [ ((i+(j-1)*nn1, j=1,nn2), i=nn1,1,-1) ]

        call surtop ( n1=m, n2=n, nn1row=m+1 )

        surface%n = m
        surface%m = n

      case(3)

        newnodes = [ ((i+(j-1)*nn1, i=nn1,1,-1), j=nn2,1,-1) ]

      case(4)

        newnodes = [ ((i+(j-1)*nn1, j=nn2,1,-1), i=1,nn1) ]

        call surtop ( n1=m, n2=n, nn1row=m+1 )

        surface%n = m
        surface%m = n

      case default

        call errormsg_case_default ( 'mesh_convert_orientation_quad4', &
          'orientation', int_value=orientation )

      end select

    else if ( flip == -1 ) then

!     flip direction

      select case ( orientation )

      case(1)

        newnodes = [ ((i+(j-1)*nn1, i=nn1,1,-1), j=1,nn2) ]

      case(2)

        newnodes = [ ((i+(j-1)*nn1, j=1,nn2), i=1,nn1) ]

        call surtop ( n1=m, n2=n, nn1row=m+1 )

        surface%n = m
        surface%m = n

      case(3)

        newnodes = [ ((i+(j-1)*nn1, i=1,nn1), j=nn2,1,-1) ]

      case(4)

        newnodes = [ ((i+(j-1)*nn1, j=nn2,1,-1), i=nn1,1,-1) ]

        call surtop ( n1=m, n2=n, nn1row=m+1 )

        surface%n = m
        surface%m = n

      case default

        call errormsg_case_default ( 'mesh_convert_orientation_quad4', &
          'orientation', int_value=orientation )

      end select

    end if

!   new nodes

    surface%nodes = surface%nodes(newnodes)

!   topology of elements on surface in global nodes

    do elem = 1, surface%nelem
      surface%topology(:,elem,2) = surface%nodes(surface%topology(:,elem,1))
    end do

    deallocate ( newnodes )

  contains

!   surface topology

    subroutine surtop ( n1, n2, nn1row )

      integer, intent(in) :: n1, n2, nn1row

      integer :: nn1, nn2, i, j

      do i = 1, n1
        do j = 1, n2

          nn1 = (j-1)*nn1row + i-1 ! number of nodes before element (i,j)
          nn2 = nn1 + nn1row   ! number of nodes before element (i,j) + 1 row

          e = i + (j-1)*n1   ! element number

          elsnodes = [ (nn1+i,i=1,2), (nn2+i,i=2,1,-1) ]

          surface%topology(:,e,1) = elsnodes(1:inpelms)

        end do
      end do

    end subroutine surtop

  end subroutine mesh_convert_orientation_quad4


! convert orientation of surface with quad9 elements

  subroutine mesh_convert_orientation_quad9 ( surface, orientation, flip )

    type(geometry_t), intent(inout) :: surface

    integer, intent(in) :: orientation, flip


    integer, parameter :: inpelms = 9 ! number of nodal points per element on
                                      ! surfaces

    integer :: elsnodes(inpelms), i, j, e, nn1, nn2, m, n, elem
    integer, allocatable, dimension(:) :: newnodes

    allocate ( newnodes(surface%nnodes) )

    n = surface%n
    m = surface%m

    nn1 = 2*n+1
    nn2 = 2*m+1

    if ( flip == 1 ) then

!     usual direction

      select case ( orientation )

      case(1)

        return  ! nothing to do

      case(2)

        newnodes = [ ((i+(j-1)*nn1, j=1,nn2), i=nn1,1,-1) ]

        call surtop ( n1=m, n2=n, nn1row=m+1 )

        surface%n = m
        surface%m = n

      case(3)

        newnodes = [ ((i+(j-1)*nn1, i=nn1,1,-1), j=nn2,1,-1) ]

      case(4)

        newnodes = [ ((i+(j-1)*nn1, j=nn2,1,-1), i=1,nn1) ]

        call surtop ( n1=m, n2=n, nn1row=m+1 )

        surface%n = m
        surface%m = n

      case default

        call errormsg_case_default ( 'mesh_convert_orientation_quad9', &
          'orientation', int_value=orientation )

      end select

    else if ( flip == -1 ) then

!     flip direction

      select case ( orientation )

      case(1)

        newnodes = [ ((i+(j-1)*nn1, i=nn1,1,-1), j=1,nn2) ]

      case(2)

        newnodes = [ ((i+(j-1)*nn1, j=1,nn2), i=1,nn1) ]

        call surtop ( n1=m, n2=n, nn1row=m+1 )

        surface%n = m
        surface%m = n

      case(3)

        newnodes = [ ((i+(j-1)*nn1, i=1,nn1), j=nn2,1,-1) ]

      case(4)

        newnodes = [ ((i+(j-1)*nn1, j=nn2,1,-1), i=nn1,1,-1) ]

        call surtop ( n1=m, n2=n, nn1row=m+1 )

        surface%n = m
        surface%m = n

      case default

        call errormsg_case_default ( 'mesh_convert_orientation_quad9', &
          'orientation', int_value=orientation )

      end select

    end if

!   new nodes

    surface%nodes = surface%nodes(newnodes)

!   topology of elements on surface in global nodes

    do elem = 1, surface%nelem
      surface%topology(:,elem,2) = surface%nodes(surface%topology(:,elem,1))
    end do

    deallocate ( newnodes )

  contains

!   surface topology

    subroutine surtop ( n1, n2, nn1row )

      integer, intent(in) :: n1, n2, nn1row

      integer :: nn1, nn2, nn3, i, j

      do i = 1, n1
        do j = 1, n2

          nn1 = 2*(j-1)*nn1row + 2*(i-1) ! number of nodes before element (i,j)
          nn2 = nn1 + nn1row   ! number of nodes before element (i,j) + 1 row
          nn3 = nn2 + nn1row   ! number of nodes before element (i,j) + 2 rows

          e = i + (j-1)*n1   ! element number

          elsnodes = [ (nn1+i,i=1,3), nn2+3, (nn3+i,i=3,1,-1), (nn2+i,i=1,2) ]

          surface%topology(:,e,1) = elsnodes(1:inpelms)

        end do
      end do

    end subroutine surtop

  end subroutine mesh_convert_orientation_quad9


! Convert mesh1 to mesh with some elements removed and/or replaced with elmesh
! NOTE: points/curves/surfaces are removed.
! NOTE: the purpose of this routine is basically for plotting xfem data

  subroutine mesh_convert_userelmesh ( mesh1, mesh, userelmesh )

!   the mesh that needs to be converted
    type(mesh_t), intent(in) :: mesh1

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh

!   routine for element level:
!     indicator = 0 : copy element
!     indicator = 1 : delete element
!     indicator = 2 : replace mesh with elmesh
    interface
      subroutine userelmesh ( mesh, elgrp, elem, indicator, elmesh )
        use kind_defs_m
        use mesh_m, only: mesh_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        integer, intent(in) :: elgrp, elem
        integer, intent(out) :: indicator
        type(mesh_t), intent(inout) :: elmesh
      end subroutine userelmesh
    end interface


    integer :: i, elgrp, elem, indicator, nnodes, j
    integer, allocatable, dimension(:) :: inodes, inodes2
    type(int_array_1d_t), dimension(mesh1%nelgrp) :: ielements
    type(mesh_t) :: elmesh, meshd
    type(mesh_t), pointer :: mesha, meshb, meshc

!   first some allocation

    allocate ( inodes(mesh1%nnodes), inodes2(mesh1%nnodes) )

    do i = 1, mesh1%nelgrp
      allocate(ielements(i)%a(mesh1%grpnumel(i)))
    end do

    allocate(mesha,meshb)

    inodes = 0

!   loop over all elements

    do elgrp = 1, mesh1%nelgrp
      do elem = 1, mesh1%grpnumel(elgrp)

!       ask user what to do for this element

        call userelmesh ( mesh1, elgrp, elem, indicator, elmesh )

        ielements(elgrp)%a(elem) = indicator

        if ( indicator == 0 ) then

!         nodes of this element need to be in the mesh

          inodes(mesh1%topology(elgrp)%a(:,elem)) = 1

        else if ( indicator == 2 ) then

!         replace element with submesh

          if ( mesha%meshgen ) then

!           merge elmesh with what we already have

            call mesh_merge ( mesha, elmesh, meshb )

            call delete ( mesha )

!           move meshb to mesha

            meshc => mesha
            mesha => meshb
            meshb => meshc

          else

!           first time: copy mesh

            call mesh_convert_copy ( elmesh, mesha )

          end if

          call delete ( elmesh )

        end if

      end do
    end do

!   which nodes need to be retained?

    nnodes = 0
    do i = 1, mesh1%nnodes
      if ( inodes(i) == 1 ) then
        nnodes = nnodes + 1
        inodes(i) = nnodes
        inodes2(nnodes) = i  ! inverse
      end if
    end do

!   now copy mesh1 to meshd (only elements that need to be retained)

    meshd%ndim = mesh1%ndim
    meshd%nnodes = nnodes
    meshd%nelgrp = mesh1%nelgrp
    meshd%nblend = 0

    allocate( meshd%element ( meshd%nelgrp ) )
    allocate( meshd%element_blend ( meshd%nelgrp, 0 ) )
    meshd%nnodes_blend = [0,meshd%nnodes]

    call copy ( mesh1%element, meshd%element )

    allocate( meshd%grpnumel ( meshd%nelgrp ) )

    do elgrp = 1, meshd%nelgrp
!     only elements that need to be retained
      meshd%grpnumel(elgrp) = count ( ielements(elgrp)%a == 0 )
    end do

    meshd%nelem = sum ( meshd%grpnumel )

    allocate( meshd%topology(meshd%nelgrp) )
    do elgrp = 1, meshd%nelgrp
      allocate( meshd%topology(elgrp)%a(meshd%element(elgrp)%numnod,&
                                        meshd%grpnumel(elgrp)) )
      j = 0
      do i = 1, mesh1%grpnumel(elgrp)
        if ( ielements(elgrp)%a(i) == 0 ) then
          j = j + 1
          meshd%topology(elgrp)%a(:,j) = inodes(mesh1%topology(elgrp)%a(:,i))
        end if
      end do
    end do
    meshd%coor = mesh1%coor(inodes2(1:nnodes),:)

!   allocate memory for points, curves, surfaces and volumes

    allocate( meshd%points(MAXPOINTS), meshd%curves(MAXCURVES), &
      meshd%surfaces(MAXSURFACES), meshd%volumes(MAXVOLUMES) )

!   blocks

    allocate ( meshd%blocks(MAXBLOCKS) )

!   objects

    allocate ( meshd%objects(MAXOBJECTS) )

!   nodesets

    allocate ( meshd%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( meshd%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    meshd%meshgen = .true.

    if ( mesha%meshgen ) then

!     merge mesh with subdivided mesh

      call mesh_merge ( meshd, mesha, mesh, nogroupmerge=.true. )

    else

      call mesh_convert ( meshd, mesh )

    end if

!   clear memory

    do i = 1, mesh1%nelgrp
      deallocate(ielements(i)%a)
    end do

    call delete (meshd)
    if ( mesha%meshgen ) call delete(mesha)

    deallocate(mesha,meshb)
    deallocate ( inodes, inodes2 )

  end subroutine mesh_convert_userelmesh


! Convert a high-order mesh to a new basic mesh with linear elements

  subroutine mesh_convert_spectraltolinear ( mesh1, mesh )

!   the mesh that needs to be copied
    type(mesh_t), intent(in) :: mesh1

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh


    integer :: elgrp, curve, surface, p, q, r, i, j, k, m, n, kelem, elem, &
      mj, elnumnod
    integer, dimension(:,:), allocatable :: mat


!   internal mesh

    mesh%ndim = mesh1%ndim
    mesh%nnodes = mesh1%nnodes
    mesh%nelgrp = mesh1%nelgrp
    mesh%nblend = 0

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    mesh%nnodes_blend = [0,mesh%nnodes]

    mesh%element(:)%ndim = mesh1%element(:)%ndim
    mesh%element(:)%globalshape = mesh1%element(:)%globalshape

    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology(mesh%nelgrp) )

    do elgrp = 1, mesh%nelgrp

      select case ( mesh1%element(elgrp)%elshape )

      case(101,105) ! high-order line element

        mesh%element(elgrp)%elshape = 1
        mesh%element(elgrp)%numnod = 2
        p = mesh1%element(elgrp)%p(1,1)
        mesh%grpnumel(elgrp) = mesh1%grpnumel(elgrp) * p

        allocate( mesh%topology(elgrp)%a(mesh%element(elgrp)%numnod, &
                                         mesh%grpnumel(elgrp)) )

        kelem = 0
        do elem = 1, mesh1%grpnumel(elgrp)
           do i = 1, p
             mesh%topology(elgrp)%a(:,kelem+i) = &
                               mesh1%topology(elgrp)%a([i,i+1],elem)
           end do
           kelem = kelem + p
        end do

      case(102,106) ! high-order quadrilateral element

        mesh%element(elgrp)%elshape = 5
        mesh%element(elgrp)%numnod = 4
        p = mesh1%element(elgrp)%p(1,1)
        q = mesh1%element(elgrp)%p(2,1)
        mesh%grpnumel(elgrp) = mesh1%grpnumel(elgrp) * p * q

        allocate( mesh%topology(elgrp)%a(mesh%element(elgrp)%numnod, &
                                         mesh%grpnumel(elgrp)) )

        kelem = 0
        do elem = 1, mesh1%grpnumel(elgrp)
           do i = 1, p
             do j = 1, q
               m = i+(j-1)*(p+1)
               mesh%topology(elgrp)%a(:,kelem+i+(j-1)*p) = &
                    mesh1%topology(elgrp)%a([m,m+1,m+p+2,m+p+1],elem)
             end do
           end do
           kelem = kelem + p * q
        end do

      case(103,107) ! high-order hexahedron element

        mesh%element(elgrp)%elshape = 13
        mesh%element(elgrp)%numnod = 8
        p = mesh1%element(elgrp)%p(1,1)
        q = mesh1%element(elgrp)%p(2,1)
        r = mesh1%element(elgrp)%p(3,1)
        mesh%grpnumel(elgrp) = mesh1%grpnumel(elgrp) * p * q * r

        allocate( mesh%topology(elgrp)%a(mesh%element(elgrp)%numnod, &
                                         mesh%grpnumel(elgrp)) )

        kelem = 0
        do elem = 1, mesh1%grpnumel(elgrp)
           do i = 1, p
             do j = 1, q
               do k = 1, r
                 m = i+(j-1)*(p+1)+(k-1)*(p+1)*(q+1)
                 n = (p+1)*(q+1)
                 mesh%topology(elgrp)%a(:,kelem+i+(j-1)*p+(k-1)*p*q) = &
                    mesh1%topology(elgrp)%a(&
                         [m,m+1,m+p+2,m+p+1,m+n,m+1+n,m+p+2+n,m+p+1+n],elem)
               end do
             end do
           end do
           kelem = kelem + p * q * r
        end do

      case(104) ! high-order triangle element

        mesh%element(elgrp)%elshape = 3
        mesh%element(elgrp)%numnod = 3
        p = mesh1%element(elgrp)%p(1,1)
        mesh%grpnumel(elgrp) = mesh1%grpnumel(elgrp) * p**2

        allocate( mesh%topology(elgrp)%a(mesh%element(elgrp)%numnod, &
                                         mesh%grpnumel(elgrp)) )

        allocate ( mat(p+1,p+1) )

        kelem = 0

        do elem = 1, mesh1%grpnumel(elgrp)

!         fill mat with nodal numbers

          do j = 1, p+1
            mj= (j-1)*(2*p+4-j)/2
            do i = 1, p+2-j
              mat(i,j) = mesh1%topology(elgrp)%a(i+mj,elem)
            end do
          end do

!         first set of elements oriented with point down left

          do j = 1, p
            mj = (j-1)*(2*p+2-j)/2
            do i = 1, p+1-j
              mesh%topology(elgrp)%a(:,kelem+i+mj) = &
                                   [ mat(i,j), mat(i+1,j), mat(i,j+1) ]
            end do
          end do

!         first set of elements oriented with point up right

          do j = 1, p-1
            mj = p*(p+1)/2 + (j-1)*(2*p-j)/2
            do i = 1, p-j
              mesh%topology(elgrp)%a(:,kelem+i+mj) = &
                                   [ mat(i+1,j), mat(i+1,j+1), mat(i,j+1) ]
            end do
          end do

          kelem = kelem + p**2

        end do

      case default ! other elements

!       just copy

        mesh%element(elgrp)%elshape = mesh1%element(elgrp)%elshape
        mesh%element(elgrp)%numnod = mesh1%element(elgrp)%numnod
        mesh%element(elgrp)%p = mesh1%element(elgrp)%p
        mesh%grpnumel(elgrp) = mesh1%grpnumel(elgrp)

        mesh%topology(elgrp)%a = mesh1%topology(elgrp)%a

      end select

    end do

    mesh%nelem = sum(mesh%grpnumel)

    mesh%coor = mesh1%coor

!   points

    mesh%npoints = mesh1%npoints

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1:mesh%npoints) = mesh1%points(1:mesh1%npoints)

!   curves

    mesh%ncurves = mesh1%ncurves

    allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )

    mesh%curves(:mesh%ncurves)%ndim = mesh1%curves(:mesh1%ncurves)%ndim
    mesh%curves(:mesh%ncurves)%element%globalshape = &
                           mesh1%curves(:mesh1%ncurves)%element%globalshape
    mesh%curves(:mesh%ncurves)%element%ndim = &
                           mesh1%curves(:mesh1%ncurves)%element%ndim
    mesh%curves(:mesh%ncurves)%nnodes = mesh1%curves(:mesh1%ncurves)%nnodes

    mesh%curves(:mesh%ncurves)%nblend = 0
    do curve = 1, mesh%ncurves
      allocate ( mesh%curves(curve)%element_blend(0) )
      mesh%curves(curve)%nnodes_blend = [0,mesh%curves(curve)%nnodes]
    end do

    do curve = 1, mesh%ncurves
      mesh%curves(curve)%nodes = mesh1%curves(curve)%nodes
    end do

    do curve = 1, mesh%ncurves

      select case ( mesh1%curves(curve)%element%elshape )

      case (101,105) ! high-order line element

        elnumnod = 2
        mesh%curves(curve)%element%elshape = 1
        mesh%curves(curve)%element%numnod = 2
        p = mesh1%curves(curve)%element%p(1,1)
        mesh%curves(curve)%nelem = mesh1%curves(curve)%nelem * p
        mesh%curves(curve)%n = 0  ! regular structure is gone

        allocate( mesh%curves(curve)%topology(elnumnod,&
                                              mesh%curves(curve)%nelem,2) )

        kelem = 0
        do elem = 1, mesh1%curves(curve)%nelem
           do i = 1, p
             mesh%curves(curve)%topology(:,kelem+i,:) = &
                            mesh1%curves(curve)%topology([i,i+1],elem,:)
           end do
           kelem = kelem + p
        end do

      case default ! other standard line elements

!       just copy

        mesh%curves(curve)%element%elshape = mesh1%curves(curve)%element%elshape
        mesh%curves(curve)%element%p = mesh1%curves(curve)%element%p
        mesh%curves(curve)%element%numnod = mesh1%curves(curve)%element%numnod
        mesh%curves(curve)%nelem = mesh1%curves(curve)%nelem
        mesh%curves(curve)%n = mesh1%curves(curve)%n

        mesh%curves(curve)%topology = mesh1%curves(curve)%topology

      end select

    end do

!   surfaces

    mesh%nsurfaces = mesh1%nsurfaces

    allocate( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )

    mesh%surfaces(:mesh%nsurfaces)%ndim = mesh1%surfaces(:mesh1%nsurfaces)%ndim
    mesh%surfaces(:mesh%nsurfaces)%element%globalshape = &
                       mesh1%surfaces(:mesh1%nsurfaces)%element%globalshape
    mesh%surfaces(:mesh%nsurfaces)%element%ndim = &
                       mesh1%surfaces(:mesh1%nsurfaces)%element%ndim
    mesh%surfaces(:mesh%nsurfaces)%nnodes = &
                       mesh1%surfaces(:mesh1%nsurfaces)%nnodes

    mesh%surfaces(:mesh%nsurfaces)%nblend = 0
    do surface = 1, mesh%nsurfaces
      allocate ( mesh%surfaces(surface)%element_blend(0) )
      mesh%surfaces(surface)%nnodes_blend = [0,mesh%surfaces(surface)%nnodes]
    end do

    do surface = 1, mesh%nsurfaces
      mesh%surfaces(surface)%nodes = mesh1%surfaces(surface)%nodes
    end do


    do surface = 1, mesh%nsurfaces

      select case ( mesh1%surfaces(surface)%element%elshape )

      case(102,106) ! high-order quadrilateral element

        elnumnod = 4
        mesh%surfaces(surface)%element%elshape = 5
        mesh%surfaces(surface)%element%numnod = 4
        p = mesh1%surfaces(surface)%element%p(1,1)
        q = mesh1%surfaces(surface)%element%p(2,1)
        mesh%surfaces(surface)%nelem = mesh1%surfaces(surface)%nelem * p * q
        mesh%surfaces(surface)%n = 0  ! structure is gone
        mesh%surfaces(surface)%m = 0  ! structure is gone

        allocate( mesh%surfaces(surface)%topology(elnumnod,&
                                           mesh%surfaces(surface)%nelem,2) )

        kelem = 0
        do elem = 1, mesh1%surfaces(surface)%nelem
           do i = 1, p
             do j = 1, q
               m = i+(j-1)*(p+1)
               mesh%surfaces(surface)%topology(:,kelem+i+(j-1)*p,:) = &
                  mesh1%surfaces(surface)%topology([m,m+1,m+p+2,m+p+1],elem,:)
             end do
           end do
           kelem = kelem + p * q
        end do

      case default ! other standard surface elements

!       just copy

        mesh%surfaces(surface)%element%elshape = &
                                   mesh1%surfaces(surface)%element%elshape
        mesh%surfaces(surface)%element%p = mesh1%surfaces(surface)%element%p
        mesh%surfaces(surface)%element%numnod = &
                                   mesh1%surfaces(surface)%element%numnod
        mesh%surfaces(surface)%nelem = mesh1%surfaces(surface)%nelem
        mesh%surfaces(surface)%n = mesh1%surfaces(surface)%n
        mesh%surfaces(surface)%m = mesh1%surfaces(surface)%m

        mesh%surfaces(surface)%topology = mesh1%surfaces(surface)%topology

      end select

    end do

!   volumes

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    mesh%meshgen = .true.

  end subroutine mesh_convert_spectraltolinear


! Copy a mesh to a new mesh where all nodes that have the same coordinates
! are replaced by a single node.

  subroutine mesh_convert_remove_double_nodes ( mesh1, mesh, excludepoints, &
    excludecurves, excludesurfaces, xmin, xmax )

    use limits_m, only: EPSNODEOVERLAP

!   the mesh that needs to be converted
    type(mesh_t), intent(in) :: mesh1

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh

!   exclude nodes in points.
!   For example excludepoints=(/1,3/) excludes the nodes given by
!   points 1 and 3.
    integer, intent(in), dimension(:), optional :: excludepoints

!   exclude nodes on curves.
!   For example excludecurves=(/1,3/) excludes the nodes given by
!   curves 1 and 3.
    integer, intent(in), dimension(:), optional :: excludecurves

!   exclude nodes on surfaces.
!   For example excludesurfaces=(/1,3/) excludes the nodes given by
!   surfaces 1 and 3.
    integer, intent(in), dimension(:), optional :: excludesurfaces

!   exclude nodes if any of the coordinates is less than the ones given by
!   xmin or any of the coordinates is larger than the ones given by xmax.
!   For example xmin=(/-1._dp,-2._dp/) excludes all nodes that have an
!   x-coordinate less than -1 and/or a y-coordinate less than -2.
    real(dp), intent(in), dimension(:), optional :: xmin, xmax


    integer :: elgrp, numnod, elem
    integer :: block, node, node2, nod

    integer, allocatable, dimension(:) :: &
      work, &! work array for nodes
      nnr    ! new node numbering for each old node in mesh1
    integer :: nnnodes ! number of nodes in new mesh
    integer :: ijk(mesh1%ndim)
    real(dp) :: domainsize(mesh1%ndim)


!   allocate work arrays

    allocate ( work(mesh1%nnodes), nnr(mesh1%nnodes) )

!   copy some data

    mesh%ndim = mesh1%ndim
    mesh%nelem = mesh1%nelem
    mesh%nelgrp = mesh1%nelgrp
    mesh%nblend = 0

!   internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )

    call copy ( mesh1%element, mesh%element )

    mesh%grpnumel = mesh1%grpnumel

    allocate( mesh%topology ( mesh%nelgrp ) )

    domainsize = mesh1%nodblocks%nbl * mesh1%nodblocks%dx


!   loop over nodes

    work = 0

    do node = 1, mesh1%nnodes

      if ( work(node) /= 0 ) cycle  ! node already removed

!     find nodblock of node

      ijk = floor ( &
         ( mesh1%coor(node,:) - mesh1%nodblocks%x0 ) / mesh1%nodblocks%dx )

      if ( all ( ijk >= 0 ) .and. all ( ijk <= mesh1%nodblocks%nbl-1 ) ) then

!       point within structured grid of nodblocks

        select case ( mesh1%ndim )
        case(1)
          block = 1 + ijk(1)
        case(2)
          block = 1 + ijk(1) + ijk(2) * mesh1%nodblocks%nbl(1)
        case(3)
          block = 1 + ijk(1) + ijk(2) * mesh1%nodblocks%nbl(1) + &
                      ijk(3) * mesh1%nodblocks%nbl(1) * mesh1%nodblocks%nbl(2)
        case default
          call errormsg_case_default ( 'mesh_convert_remove_double_nodes', &
            'mesh1%ndim', int_value=mesh1%ndim )
        end select

!       loop over all nodes in nodblock

        do nod = 1, mesh1%nodblocks%nodblks(block)%nnodes

          node2 = mesh1%nodblocks%nodblks(block)%nodes(nod)

!         check only upper triangle of connection matrix
          if ( node2 <= node ) cycle

          if ( work(node2) /= 0 ) cycle  ! node already removed

          if ( all ( abs(mesh1%coor(node2,:)-mesh1%coor(node,:) ) &
                              < EPSNODEOVERLAP*domainsize ) ) then
!           remove node
            work(node2) = node  !NOTE: node < node2
          end if

        end do

      end if

    end do

!   exclude some nodes

    call excludenodes

!   generate new numbering

    nnnodes = 0

    do node = 1, mesh1%nnodes

      if ( work(node) == 0 ) then
  !     new node
        nnnodes = nnnodes + 1
        nnr(node) = nnnodes
      else
!       removed node
        nnr(node) = nnr(work(node))
      end if

    end do

!   new topology

    do elgrp = 1, mesh%nelgrp

!     copy topology with new numbering

      numnod = mesh1%element(elgrp)%numnod

      allocate( mesh%topology(elgrp)%a(numnod,mesh%grpnumel(elgrp)) )

!     loop over elements in this group

      do elem = 1, mesh%grpnumel(elgrp)

        mesh%topology(elgrp)%a(:,elem) = nnr(mesh1%topology(elgrp)%a(:,elem))

      end do

    end do

!   new coordinates

    mesh%nnodes = nnnodes

    mesh%nnodes_blend = [0,mesh%nnodes]

    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   copy old coordinates

    do node = 1, mesh1%nnodes
      if ( work(node) == 0 ) then
        mesh%coor(nnr(node),:) = mesh1%coor(node,:)
      end if
    end do

!   points (transfer to new numbering with nnr).

    mesh%npoints = mesh1%npoints

    allocate( mesh%points(max(mesh%npoints,MAXPOINTS)) )

    mesh%points(1:mesh%npoints) = nnr(mesh1%points(1:mesh1%npoints))

!   curves (transfer to new numbering with nnr).

    mesh%ncurves = mesh1%ncurves

    allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )

    call copy ( mesh1%curves(:mesh1%ncurves), mesh%curves(:mesh%ncurves), nnr )

!   surfaces (transfer to new numbering with nnr).

    mesh%nsurfaces = mesh1%nsurfaces

    allocate( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )

    call copy ( mesh1%surfaces(:mesh1%nsurfaces), &
                mesh%surfaces(:mesh%nsurfaces), nnr )

!   remove double points in points, curves and surfaces

    call mesh_change_remove_double_nodes_geometries ( mesh )

!   volumes

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    mesh%meshgen = .true.

    deallocate ( work, nnr )

  contains

!   excludenodes from points, curves and surfaces

    subroutine excludenodes

      integer :: crv, srf, node

      if ( present(excludepoints) ) then
        if ( any(excludepoints <=0) .or. &
             any(excludepoints > mesh1%npoints) ) then
          write(*,'(2(/a),i0/)') &
            'Error in mesh_convert_remove_double_nodes: ',&
            ' excludepoints <=0 or excludepoints > number of points = ', &
             mesh1%npoints
          stop
        end if
        work ( mesh1%points(excludepoints) ) = 0
      end if

      if ( present(excludecurves) ) then
        if ( any(excludecurves <=0) .or. &
             any(excludecurves > mesh1%ncurves) ) then
          write(*,'(2(/a),i0/)') &
            'Error in mesh_convert_remove_double_nodes: ',&
            ' excludecurves <=0 or excludecurves > number of curves = ', &
             mesh1%ncurves
          stop
        end if
        do crv = 1, size(excludecurves)
          work ( mesh1%curves(excludecurves(crv))%nodes ) = 0
        end do
      end if

      if ( present(excludesurfaces) ) then
        if ( any(excludesurfaces <=0) .or. &
             any(excludesurfaces > mesh1%nsurfaces) ) then
          write(*,'(2(/a),i0/)') &
            'Error in mesh_convert_remove_double_nodes: ',&
            ' excludesurfaces <=0 or excludesurfaces > number of surfaces = ', &
             mesh1%nsurfaces
          stop
        end if
        do srf = 1, size(excludesurfaces)
          work ( mesh1%surfaces(excludesurfaces(srf))%nodes ) = 0
        end do
      end if

      if ( present(xmin) ) then
        if ( size(xmin) /= mesh1%ndim ) then
          write(*,'(2(/a),i0/)') &
            'Error in mesh_convert_remove_double_nodes: ',&
            ' dimension of xmin must be equal to the space dimension = ', &
             mesh1%ndim
          stop
        end if
        do node = 1, mesh1%nnodes
          if ( all ( mesh1%coor(node,:) > xmin ) ) cycle  ! in region
          work ( node ) = 0
        end do
      end if

      if ( present(xmax) ) then
        if ( size(xmax) /= mesh1%ndim ) then
          write(*,'(2(/a),i0/)') &
            'Error in mesh_convert_remove_double_nodes: ',&
            ' dimension of xmax must be equal to the space dimension = ', &
             mesh1%ndim
          stop
        end if
        do node = 1, mesh1%nnodes
          if ( all ( mesh1%coor(node,:) < xmax ) ) cycle  ! in region
          work ( node ) = 0
        end do
      end if

    end subroutine excludenodes

  end subroutine mesh_convert_remove_double_nodes


! Convert a mesh to a new basic mesh with isolated nodes removed.
! Note, that points connected to nodal point zero (sepran) are removed as well.
! Note, that geometries (points, curves, surfaces) that are connected to
! isolated nodes are removed as well.

  subroutine mesh_convert_remove_isolated_nodes ( mesh1, mesh )

!   the mesh that needs to be converted
    type(mesh_t), intent(in) :: mesh1

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh


    integer :: elgrp, elem, nnodes, i, nremoved, np, nc, ns
    integer :: point, curve, surface
    integer, allocatable, dimension(:) :: newnum
    integer :: removed(maxval([mesh1%npoints,mesh1%ncurves,mesh1%nsurfaces]))


!   copy some data

    mesh%ndim = mesh1%ndim
    mesh%nelem = mesh1%nelem
    mesh%nelgrp = mesh1%nelgrp
    mesh%nblend = 0

!   find isolated nodes

    allocate ( newnum(mesh1%nnodes) )

    newnum = 0

    do elgrp = 1, mesh1%nelgrp
      do elem = 1, mesh1%grpnumel(elgrp)
        newnum(mesh1%topology(elgrp)%a(:,elem)) = 1
      end do
    end do

!   set new numbering

    nnodes = 0
    do i = 1, mesh1%nnodes
      if ( newnum(i) == 0 ) cycle
      nnodes = nnodes + 1  ! node found
      newnum(i) = nnodes
    end do

    mesh%nnodes = nnodes

!   internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    mesh%nnodes_blend = [0,mesh%nnodes]

    call copy ( mesh1%element, mesh%element )

    mesh%grpnumel = mesh1%grpnumel

    allocate( mesh%topology(mesh%nelgrp) )
    do elgrp = 1, mesh%nelgrp
      allocate(&
        mesh%topology(elgrp)%a(mesh%element(elgrp)%numnod,mesh%grpnumel(elgrp)))
      do elem = 1, mesh%grpnumel(elgrp)
        mesh%topology(elgrp)%a(:,elem) = newnum(mesh1%topology(elgrp)%a(:,elem))
      end do
    end do

    allocate( mesh%coor ( nnodes, mesh%ndim ) )

    do i = 1, mesh1%nnodes
      if ( newnum(i) == 0 ) cycle
      mesh%coor(newnum(i),:) = mesh1%coor(i,:)
    end do

    removed = 0
    nremoved = 0

!   points

    allocate( mesh%points(max(mesh1%npoints,MAXPOINTS)))

    np = 0

    do point = 1, mesh1%npoints
      if ( mesh1%points(point) == 0 ) then
!       point removed since zero nodal point (probably sepran)
        nremoved = nremoved + 1
        removed(nremoved) = point
      else if ( newnum(mesh1%points(point)) == 0 ) then
!       point removed since corresponding node remove
        nremoved = nremoved + 1
        removed(nremoved) = point
      else
!       point copied to new mesh
        np = np + 1
        mesh%points(np) = newnum(mesh1%points(point))
      end if
    end do

    mesh%npoints = np

    if ( nremoved /= 0 ) then
      write(*,'(/a/a)') &
        'Warning in mesh_convert_remove_isolated_nodes: ', &
        ' The following point numbers have been removed: '
      write(*,'(25(i0,1x))') removed(1:nremoved)
      removed(1:nremoved) = 0
      nremoved = 0
    end if

!   curves

    allocate( mesh%curves(max(mesh1%ncurves,MAXCURVES)) )

    nc = 0

    do curve = 1, mesh1%ncurves
      if ( any ( newnum(mesh1%curves(curve)%nodes) == 0 ) ) then
!       curve removed
        nremoved = nremoved + 1
        removed(nremoved) = curve
      else
!       curve copied to new mesh
        nc = nc + 1
        call copy ( mesh1%curves(curve), mesh%curves(nc), newnum )
      end if
    end do

    mesh%ncurves = nc

    if ( nremoved /= 0 ) then
      write(*,'(/a/a)') &
        'Warning in mesh_convert_remove_isolated_nodes: ', &
        ' The following curve numbers have been removed: '
      write(*,'(25(i0,1x))') removed(1:nremoved)
      removed(1:nremoved) = 0
      nremoved = 0
    end if

!   surfaces

    allocate( mesh%surfaces(max(mesh1%nsurfaces,MAXSURFACES)) )

    ns = 0

    do surface = 1, mesh1%nsurfaces
      if ( any ( newnum(mesh1%surfaces(surface)%nodes) == 0 ) ) then
!       surface removed
        nremoved = nremoved + 1
        removed(nremoved) = surface
      else
!       surface copied to new mesh
        ns = ns + 1
        call copy ( mesh1%surfaces(surface), mesh%surfaces(ns), newnum )
      end if
    end do

    mesh%nsurfaces = ns

    if ( nremoved /= 0 ) then
      write(*,'(/a/a)') &
        'Warning in mesh_convert_remove_isolated_nodes: ', &
        ' The following surface numbers have been removed: '
      write(*,'(25(i0,1x))') removed(1:nremoved)
      removed(1:nremoved) = 0
      nremoved = 0
    end if

!   volumes

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    mesh%meshgen = .true.

    deallocate ( newnum )

  end subroutine mesh_convert_remove_isolated_nodes


! Convert mesh1 to a mesh with only the specified groups present.
! For example: only_groups=(/1,3/).
! Note, that geometries (points, curves, surfaces) which have all nodes within
! the specified groups are preserved. Others are discarded. The numbering of
! the geometries might be different.

  subroutine mesh_convert_only_groups ( mesh1, mesh, only_groups )

!   the mesh that needs to be copied
    type(mesh_t), intent(in) :: mesh1

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh

!   the groups specified
    integer, intent(in), dimension(:), optional :: only_groups

    integer :: elgrp, elem, node, nr, i, grp


    integer, allocatable, dimension(:) :: work
    integer, dimension(mesh1%ncurves) :: curves
    integer, dimension(mesh1%nsurfaces) :: surfaces

    allocate ( work(mesh1%nnodes) )

!   nodes

    work = 0

    do grp = 1, size(only_groups)
      elgrp = only_groups(grp)
      do elem = 1, mesh1%grpnumel(elgrp)
        work(mesh1%topology(elgrp)%a(:,elem)) = 1
      end do
    end do

!   new node numbering

    nr = 0
    do node = 1, mesh1%nnodes
      if ( work(node) == 0 ) cycle
      nr = nr + 1
      work(node) = nr
    end do

!   copy elements of specified groups only

    mesh%ndim = mesh1%ndim
    mesh%nelem = sum(mesh1%grpnumel(only_groups))
    mesh%nnodes = nr
    mesh%nelgrp = size(only_groups)
    mesh%nblend = 0

!   internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    mesh%nnodes_blend = [0,mesh%nnodes]

    call copy ( mesh1%element(only_groups), mesh%element )

    mesh%grpnumel = mesh1%grpnumel(only_groups)

    allocate( mesh%topology(mesh%nelgrp) )
    do elgrp = 1, mesh%nelgrp
      allocate(&
        mesh%topology(elgrp)%a(mesh%element(elgrp)%numnod,mesh%grpnumel(elgrp)))
      do elem = 1, mesh%grpnumel(elgrp)
        mesh%topology(elgrp)%a(:,elem) = &
                           work(mesh1%topology(only_groups(elgrp))%a(:,elem))
      end do
    end do
    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )
    do node = 1, mesh1%nnodes
      if ( work(node) == 0 ) cycle
      mesh%coor(work(node),:) = mesh1%coor(node,:)
    end do

!   points

    mesh%npoints = count( work(mesh1%points(:mesh1%npoints)) /= 0 )

    allocate( mesh%points(max(mesh%npoints,MAXPOINTS)) )

    nr = 0
    do i = 1, mesh1%npoints
      if ( work(mesh1%points(i)) == 0 ) cycle
      nr = nr + 1
      mesh%points(nr) = work(mesh1%points(i))
    end do

!   curves

    mesh%ncurves = 0

    do i = 1, mesh1%ncurves
      if ( any ( work(mesh1%curves(i)%nodes) == 0 ) ) cycle
      mesh%ncurves = mesh%ncurves + 1
      curves(mesh%ncurves) = i
    end do

    allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )

    call copy ( mesh1%curves(curves(:mesh%ncurves)), &
                mesh%curves(:mesh%ncurves), work )

!   surfaces

    mesh%nsurfaces = 0

    do i = 1, mesh1%nsurfaces
      if ( any ( work(mesh1%surfaces(i)%nodes) == 0 ) ) cycle
      mesh%nsurfaces = mesh%nsurfaces + 1
      surfaces(mesh%nsurfaces) = i
    end do

    allocate( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )

    call copy ( mesh1%surfaces(surfaces(:mesh%nsurfaces)), &
                mesh%surfaces(:mesh%nsurfaces), work )

!   volumes

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    mesh%meshgen = .true.

    deallocate ( work )

  end subroutine mesh_convert_only_groups


! Convert geometry to a mesh with only a single group.

  subroutine mesh_convert_from_geometry ( mesh1, geometry, mesh )

!   the mesh that contains the geometry
    type(mesh_t), intent(in) :: mesh1

!   the geometry
    type(geometry_t), intent(in) :: geometry

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh

!   make skeleton mesh
    call mesh_skeleton ( mesh=mesh, nnodes=geometry%nnodes, &
      nelem=geometry%nelem, elshape=geometry%element%elshape, &
      ndim=geometry%ndim, callname='mesh_convert_from_geometry' )

!   fill coordinates
    mesh%coor = mesh1%coor(geometry%nodes,:)

!   fill topology
    mesh%topology(1)%a = geometry%topology(:,:,1)

  end subroutine mesh_convert_from_geometry


! Convert a mesh to a new basic mesh with a blended mesh added

  subroutine mesh_convert_add_blend ( mesh1, blendmesh, mesh )

!   the mesh to which the blend mesh needs to be added
    type(mesh_t), intent(in) :: mesh1

!   the blend mesh
    type(mesh_t), intent(in) :: blendmesh

!   the new mesh created
    type(mesh_t), intent(out) :: mesh


    integer :: elgrp, numnod, numnod1, curve, surface


!   copy all elements from mesh1 and blendmesh

    mesh%ndim = mesh1%ndim
    mesh%nelem = mesh1%nelem
    mesh%nnodes = mesh1%nnodes + blendmesh%nnodes
    mesh%nelgrp = mesh1%nelgrp
    mesh%nblend = mesh1%nblend + 1 + blendmesh%nblend

!   internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, mesh%nblend ) )
    allocate( mesh%nnodes_blend(mesh%nblend+2) )

    call copy ( mesh1%element, mesh%element )
    if ( mesh1%nblend > 0 ) then
      call copy ( mesh1%element_blend, mesh%element_blend(:,1:mesh1%nblend) )
    end if
    call copy ( blendmesh%element, mesh%element_blend(:,mesh1%nblend+1) )
    if ( blendmesh%nblend > 0 ) then
      call copy ( blendmesh%element_blend, &
                                       mesh%element_blend(:,mesh1%nblend+2:) )
    end if

    mesh%nnodes_blend(1:mesh1%nblend+2) = mesh1%nnodes_blend
    mesh%nnodes_blend(mesh1%nblend+2:) = blendmesh%nnodes_blend + mesh1%nnodes

    mesh%grpnumel = mesh1%grpnumel

    allocate( mesh%topology(mesh%nelgrp) )
    do elgrp = 1, mesh%nelgrp
      numnod = mesh%element(elgrp)%numnod + &
                        sum(mesh%element_blend(elgrp,:)%numnod)
      allocate(mesh%topology(elgrp)%a(numnod,mesh%grpnumel(elgrp)))
      numnod1 = mesh1%element(elgrp)%numnod + &
                        sum(mesh1%element_blend(elgrp,:)%numnod)
      mesh%topology(elgrp)%a(1:numnod1,:) = mesh1%topology(elgrp)%a
      mesh%topology(elgrp)%a(numnod1+1:,:) = blendmesh%topology(elgrp)%a + &
                                                                 mesh1%nnodes
    end do
    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )
    mesh%coor(1:mesh1%nnodes,:) = mesh1%coor
    mesh%coor(mesh1%nnodes+1:,:) = blendmesh%coor

!   points

    mesh%npoints = mesh1%npoints + blendmesh%npoints

    allocate( mesh%points(max(mesh%npoints,MAXPOINTS)) )

    mesh%points(1:mesh1%npoints) = mesh1%points(1:mesh1%npoints)
    mesh%points(mesh1%npoints+1:mesh%npoints) = &
                   blendmesh%points(1:blendmesh%npoints) + mesh1%nnodes

!   curves

    mesh%ncurves = mesh1%ncurves

    allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )

    do curve = 1, mesh%ncurves
      call geometry_convert_add_blend ( mesh1%curves(curve), mesh1%nnodes, &
        blendmesh%curves(curve), mesh%curves(curve) )
    end do

!   surfaces

    mesh%nsurfaces = mesh1%nsurfaces

    allocate( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )

    do surface = 1, mesh%nsurfaces
      call geometry_convert_add_blend ( mesh1%surfaces(surface), &
        mesh1%nnodes, blendmesh%surfaces(surface), mesh%surfaces(surface) )
    end do

!   volumes

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    mesh%meshgen = .true.

  end subroutine mesh_convert_add_blend


! add geometry from blend mesh to geometry

  subroutine geometry_convert_add_blend ( geometry1, mesh1nnodes, &
    blendgeometry, geometry )

!   the geometry where geometry from blend mesh need to be added
    type(geometry_t), intent(in) :: geometry1

!   the value of mesh1%nnodes of which geometry1 is a component
    integer, intent(in) :: mesh1nnodes

!   the geometry that needs to be added
    type(geometry_t), intent(in) :: blendgeometry

!   the new geometry created
    type(geometry_t), intent(out) :: geometry

    integer :: elnumnod1

    geometry%meshparts = .false.
    geometry%ndim = geometry1%ndim
    geometry%nnodes = geometry1%nnodes + blendgeometry%nnodes
    geometry%nelem = geometry1%nelem
    geometry%nblend = geometry1%nblend + 1 + blendgeometry%nblend
    geometry%n = geometry1%n
    geometry%m = geometry1%m

    allocate ( geometry%element_blend(geometry%nblend) )
    allocate ( geometry%nnodes_blend ( geometry%nblend + 2 ) )

    call copy ( geometry1%element, geometry%element )
    if ( geometry1%nblend > 0 ) then
      call copy ( geometry1%element_blend, &
                               geometry%element_blend(1:geometry1%nblend) )
    end if
    call copy ( blendgeometry%element, &
                             geometry%element_blend(geometry1%nblend+1) )
    if ( blendgeometry%nblend > 0 ) then
      call copy ( blendgeometry%element_blend, &
                               geometry%element_blend(geometry1%nblend+2:) )
    end if

    geometry%nnodes_blend(1:geometry1%nblend+2) = geometry1%nnodes_blend
    geometry%nnodes_blend(geometry1%nblend+2:) = &
                                blendgeometry%nnodes_blend + geometry1%nnodes

    allocate(geometry%nodes(geometry%nnodes))

    geometry%nodes(1:geometry1%nnodes) = geometry1%nodes
    geometry%nodes(geometry1%nnodes+1:) = blendgeometry%nodes + mesh1nnodes

    elnumnod1 = size(geometry1%topology,1)
    allocate( geometry%topology(elnumnod1+size(blendgeometry%topology,1),&
              geometry%nelem,2) )

    geometry%topology(1:elnumnod1,:,:) = geometry1%topology
    geometry%topology(elnumnod1+1:,:,1) = &
                        blendgeometry%topology(:,:,1) + geometry1%nnodes
    geometry%topology(elnumnod1+1:,:,2) = &
                        blendgeometry%topology(:,:,2) + mesh1nnodes

  end subroutine geometry_convert_add_blend


! Convert mesh1 to a mesh with a single shape of elements for each group.

  subroutine mesh_convert_single_shape ( mesh1, blend, mesh )

!   the mesh from where single shape mesh must be extracted
    type(mesh_t), intent(in) :: mesh1

!   blend=0 extracts the main mesh
!   blend>0 extracts the blend mesh with shapes in element_blend.
    integer, intent(in) :: blend

!   the new mesh created
    type(mesh_t), intent(out) :: mesh


    integer :: elgrp, curve, surface, pnt, point
    type(element_t), dimension(mesh1%nelgrp) :: element

!   choose element

    if ( blend == 0 ) then
      element = mesh1%element
    else
      element = mesh1%element_blend(:,blend)
    end if

!   basic parameters of the mesh

    mesh%ndim = mesh1%ndim
    mesh%nelem = mesh1%nelem
    mesh%nnodes = mesh1%nnodes_blend(blend+2)-mesh1%nnodes_blend(blend+1)
    mesh%nelgrp = mesh1%nelgrp
    mesh%nblend = 0

!   internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )

    mesh%nnodes_blend = [0,mesh%nnodes]

    if ( blend > 0 ) then
      call copy ( mesh1%element_blend(:,blend), mesh%element )
    else
      call copy ( mesh1%element, mesh%element )
    end if

    mesh%grpnumel = mesh1%grpnumel

    allocate( mesh%topology(mesh%nelgrp) )

    do elgrp = 1, mesh%nelgrp
      allocate(mesh%topology(elgrp)%a(element(elgrp)%numnod, &
                                      mesh%grpnumel(elgrp)))
      mesh%topology(elgrp)%a = mesh1%topology(elgrp)%a(&
        mesh1%numnodtop(elgrp,blend+1)+1:mesh1%numnodtop(elgrp,blend+2),:) &
          -mesh1%nnodes_blend(blend+1)
    end do

    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )
    mesh%coor = &
       mesh1%coor(mesh1%nnodes_blend(blend+1)+1:mesh1%nnodes_blend(blend+2),:)

!   points

    mesh%npoints = count ( mesh1%points(1:mesh1%npoints) > &
                           mesh1%nnodes_blend(blend+1) .and. &
                      mesh1%points(1:mesh1%npoints) <= &
                           mesh1%nnodes_blend(blend+2) )

    allocate( mesh%points(max(mesh%npoints,MAXPOINTS)) )

    pnt = 0
    do point = 1, mesh1%npoints
      if ( mesh1%points(point) <= mesh1%nnodes_blend(blend+1) .or. &
           mesh1%points(point)  > mesh1%nnodes_blend(blend+2) ) cycle
      pnt = pnt + 1
      mesh%points(pnt) = mesh1%points(point) - mesh1%nnodes_blend(blend+1)
    end do

!   curves

    mesh%ncurves = mesh1%ncurves

    allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )

    do curve = 1, mesh%ncurves
      call geometry_convert_single_shape ( mesh1%curves(curve), &
        mesh%curves(curve) )
    end do

!   surfaces

    mesh%nsurfaces = mesh1%nsurfaces

    allocate( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )

    do surface = 1, mesh%nsurfaces
      call geometry_convert_single_shape ( mesh1%surfaces(surface), &
        mesh%surfaces(surface) )
    end do

!   volumes

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    mesh%meshgen = .true.

  contains

!   extract geometry

    subroutine geometry_convert_single_shape ( geometry1, geometry )

!     the geometry from where single shape must be extracted
      type(geometry_t), intent(in) :: geometry1

!     the new geometry created
      type(geometry_t), intent(out) :: geometry

      type(element_t) :: element

      integer :: i

!     choose element

      if ( blend == 0 ) then
        element = geometry1%element
      else
        element = geometry1%element_blend(blend)
      end if

      geometry%meshparts = .false.
      geometry%ndim = geometry1%ndim
      geometry%nnodes = geometry1%nnodes_blend(blend+2) &
                            - geometry1%nnodes_blend(blend+1)
      geometry%nelem = geometry1%nelem
      geometry%nblend = 0
      geometry%n = geometry1%n
      geometry%m = geometry1%m

      allocate ( geometry%element_blend(0) )
      geometry%nnodes_blend = [0,geometry%nnodes]

      call copy ( element, geometry%element )

      allocate(geometry%nodes(geometry%nnodes))

      geometry%nodes = geometry1%nodes(geometry1%nnodes_blend(blend+1)+1:&
                                       geometry1%nnodes_blend(blend+2)) &
                        - mesh1%nnodes_blend(blend+1)

      allocate( geometry%topology(geometry%element%numnod,geometry%nelem,2) )

      geometry%topology = geometry1%topology(&
            [(geometry1%numnodtop(blend+1)+i,i=1,geometry%element%numnod)],:,:)

      geometry%topology(:,:,1) = geometry%topology(:,:,1) &
                                         - geometry1%nnodes_blend(blend+1)
      geometry%topology(:,:,2) = geometry%topology(:,:,2) &
                                         - mesh1%nnodes_blend(blend+1)

    end subroutine geometry_convert_single_shape

  end subroutine mesh_convert_single_shape


! Change a basic mesh (in place).
! NOTE: a call with no options leaves the mesh unchanged.

  subroutine mesh_change ( mesh, regular_to_macro, macro_to_regular, &
    remove_double_nodes_geometries, groups, blends, curves, surfaces, &
    volumes, warn )

!   the mesh that needs to be changed
    type(mesh_t), intent(inout) :: mesh

!   By setting to .true., the elements are changed to "macro" element shapes,
!   where possible, i.e. elements that consist of linear, bilinear or trilinear
!   subelements.
!   Note, a basic mesh is not required, although elements in non_basic
!   components like objects are ignored.
!   default=.false.
    logical, intent(in), optional :: regular_to_macro

!   By setting to .true., the "macro" elements are changed to "regular" element
!   shapes, where possible, i.e. elements that consist of linear, bilinear
!   or trilinear subelements becomes regular (high-order) elements.
!   Note, a basic mesh is not required, although elements in non_basic
!   components like objects are ignored.
!   default=.false.
    logical, intent(in), optional :: macro_to_regular

!   By setting to .true., nodes in points and geometries (curves, surfaces and
!   volumes) pointing to the same mesh node are replaced by a single node.
!   As a result nodes in points and geometries are unique.
!   NOTE: points are considered as a single "geometry" for this option.
!   This means that each point will contain a unique mesh node.
!   default=.false.
    logical, intent(in), optional :: remove_double_nodes_geometries

!   If present: the element groups to include.
!   Applies to: regular_to_macro and macro_to_regular.
!   For example groups[2,4] will include only element groups 2 and 4.
!   Default: all groups
    integer, dimension(:), intent(in), optional :: groups

!   If present: the blend meshes to include.
!   Applies to: regular_to_macro and macro_to_regular.
!   For example blends[0,2] will include only blend mesh 0 and 2.
!   Note: 0 means the main mesh.
!   Default: all blend meshes including the main mesh.
    integer, dimension(:), intent(in), optional :: blends

!   If present: the curves to include.
!   Applies to: regular_to_macro and macro_to_regular.
!   For example curves[2,3] will include only curves 2 and 3.
!   Default: all curves
    integer, dimension(:), intent(in), optional :: curves

!   If present: the surfaces to include.
!   Applies to: regular_to_macro and macro_to_regular.
!   For example surfaces[2,3] will include only surfaces 2 and 3.
!   Default: all surfaces
    integer, dimension(:), intent(in), optional :: surfaces

!   If present: the volumes to include.
!   Applies to: regular_to_macro and macro_to_regular.
!   For example volumes[2,3] will include only volumes 2 and 3.
!   Default: all volumes
    integer, dimension(:), intent(in), optional :: volumes

!   Give warnings. By setting to .false. the warnings are suppressed.
!   default=.true.
    logical, intent(in), optional :: warn


    logical :: lwarn, lremove_double_nodes_geometries, lregular_to_macro, &
      lmacro_to_regular
    logical :: lgroups(mesh%nelgrp), lblends(0:mesh%nblend), &
               lcurves(mesh%ncurves), lsurfaces(mesh%nsurfaces), &
               lvolumes(mesh%nvolumes)
    integer :: numkey


!   some testing

    if ( .not. mesh%meshgen ) then
      write(*,'(/a/)') &
        'Error in mesh_change: mesh is empty'
      stop
    end if

    if ( .not. ( present(macro_to_regular) .or. &
                 present(regular_to_macro) ) ) then
      if ( mesh%meshparts ) then
        write(*,'(/3(a/))') &
          'Error in mesh_change: ', &
          ' Mesh parts already filled for this mesh.', &
          ' Use mesh_convert to convert to a basic mesh first.'
        stop
      end if
    end if

    if ( present(groups) ) then
!     specified groups only
      if ( any ( groups < 1 ) .or. any ( groups > mesh%nelgrp ) ) then
        write(*,'(/2(a/),a,i0/)') &
        'Error mesh_change: argument groups out of range:', &
        ' some groups are < 1 or larger than the number of groups ', mesh%nelgrp
        stop
      end if
      lgroups = .false.
      lgroups(groups) = .true.
    else
      lgroups = .true.
    end if

    if ( present(blends) ) then
!     specified blend meshes only
      if ( any ( blends < 0 ) .or. any ( blends > mesh%nblend ) ) then
        write(*,'(/2(a/),a,i0/)') &
        'Error mesh_change: argument blends out of range:', &
        ' some are < 0 or larger than the number of blend meshes ', mesh%nblend
        stop
      end if
      lblends = .false.
      lblends(blends) = .true.
    else
      lblends = .true.
    end if

    if ( present(curves) ) then
!     specified curves only
      if ( any ( curves < 1 ) .or. any ( curves > mesh%ncurves ) ) then
        write(*,'(/2(a/),a,i0/)') &
        'Error mesh_change: argument curves out of range:', &
        ' some curves are < 1 or larger than the number of curves ', &
                                                                 mesh%ncurves
        stop
      end if
      lcurves = .false.
      lcurves(curves) = .true.
    else
      lcurves = .true.
    end if

    if ( present(surfaces) ) then
!     specified surfaces only
      if ( any ( surfaces < 1 ) .or. any ( surfaces > mesh%nsurfaces ) ) then
        write(*,'(/2(a/),a,i0/)') &
        'Error mesh_change: argument surfaces out of range:', &
        ' some surfaces are < 1 or larger than the number of surfaces ', &
                                                                 mesh%nsurfaces
        stop
      end if
      lsurfaces = .false.
      lsurfaces(surfaces) = .true.
    else
      lsurfaces = .true.
    end if

    if ( present(volumes) ) then
!     specified volumes only
      if ( any ( volumes < 1 ) .or. any ( volumes > mesh%nvolumes ) ) then
        write(*,'(/2(a/),a,i0/)') &
        'Error mesh_change: argument volumes out of range:', &
        ' some volumes are < 1 or larger than the number of volumes ', &
                                                                 mesh%nvolumes
        stop
      end if
      lvolumes = .false.
      lvolumes(volumes) = .true.
    else
      lvolumes = .true.
    end if

    lwarn = set_optional ( variable=warn, default=.true. )

    numkey = 0

!   remove double nodes from geometries

    if ( present(remove_double_nodes_geometries) ) then

      call check_blend ( mesh, 'mesh_change', &
                         keyword='remove_double_nodes_geometries' )

      numkey = numkey + 1

    end if

!   regular to macro

    if ( present(regular_to_macro) ) then

      numkey = numkey + 1

    end if

!   regular to macro

    if ( present(macro_to_regular) ) then

      numkey = numkey + 1

    end if

!   Too many keywords?

    if ( numkey > 1 ) then
      write(*,'(/5(a/))') &
        'Error in mesh_change: ', &
        ' heading parameters from the following list cannot be combined: ', &
        '   remove_double_nodes_geometries ', &
        '   regular_to_macro ', &
        '   macro_to_regular '
      stop
    end if


    lremove_double_nodes_geometries = &
      set_optional ( variable=remove_double_nodes_geometries, default=.false. )

    if ( lremove_double_nodes_geometries ) then

!     remove multiple nodes in geometries

      if ( lwarn ) then

        if ( all( [ mesh%npoints, mesh%ncurves, mesh%nsurfaces, &
                    mesh%nvolumes ] == 0 ) ) then
          write(*,'(/3(a/))') &
            'Warning in mesh_change: ', &
            '  mesh does not have any points or geometries, whereas ', &
            '  remove_double_nodes_geometries = .true. '
        end if

      end if

      call mesh_change_remove_double_nodes_geometries ( mesh )

    end if


    lregular_to_macro = &
               set_optional ( variable=regular_to_macro, default=.false. )

    if ( lregular_to_macro ) then

!     regular mesh to macro mesh

      call mesh_change_regular_to_macro ( mesh, lgroups, lblends, lcurves, &
        lsurfaces, lvolumes )

    end if


    lmacro_to_regular = &
              set_optional ( variable=macro_to_regular, default=.false. )

    if ( lmacro_to_regular ) then

!     macro mesh to regular mesh

      call mesh_change_macro_to_regular ( mesh, lgroups, lblends, lcurves, &
        lsurfaces, lvolumes )

    end if

  end subroutine mesh_change


! Change a basic mesh such that points and geometries have unique node numbers

  subroutine mesh_change_remove_double_nodes_geometries ( mesh )

    use limits_m, only: REMOVE_DOUBLE_POINTS

!   the mesh that needs to be changed
    type(mesh_t), intent(inout) :: mesh

    integer :: curve, surface, volume
    integer, allocatable, dimension(:) :: &
      work, &! work array for nodes
      tmp    ! temporary allocatable


!   allocate work array

    allocate ( work(mesh%nnodes) )

!   points

    if ( REMOVE_DOUBLE_POINTS ) then

      call remove_dbl ( mesh%points(1:mesh%npoints), tmp )
      deallocate( mesh%points )
      mesh%npoints = size(tmp)
      allocate( mesh%points(max(mesh%npoints,MAXPOINTS)) )
      mesh%points(1:mesh%npoints) = tmp
      deallocate ( tmp )

    end if

!   curves

    do curve = 1, mesh%ncurves
      call remove_nodes_geom ( mesh%curves(curve) )
    end do

!   surfaces

    do surface = 1, mesh%nsurfaces
      call remove_nodes_geom ( mesh%surfaces(surface) )
    end do

!   volumes

    do volume = 1, mesh%nvolumes
      call remove_nodes_geom ( mesh%volumes(volume) )
    end do

  contains

!   routine needed for removing double nodes from points and geometries
!   NOTE: work array needs to be allocated (size=number of nodes in the mesh).

    subroutine remove_dbl ( nodin, nodout, nnrl )

!     nodal numbers, from which doubles need to be removed
      integer, dimension(:), intent(in) :: nodin

!     nodal numbers with doubles removed
      integer, dimension(:), allocatable, intent(out) :: nodout

!     local nodal numbers after doubles removed (size==size(nodin))
      integer, dimension(:), allocatable, intent(out), optional :: nnrl

      integer :: nn, wk(size(nodin)), i

      work = 0

!     find new number of nodes and fill work array's

      nn = 0

      do i = 1, size(nodin)
        if ( work(nodin(i)) /= 0 ) cycle ! node already there
        nn = nn + 1
        work(nodin(i)) = nn
        wk(nn) = nodin(i)
      end do

      if ( present(nnrl) ) then
         allocate ( nnrl(size(nodin) ) )
         nnrl = work(nodin)
      end if

      allocate ( nodout(nn) )

      nodout = wk(1:nn)

    end subroutine remove_dbl

!   remove nodes from geometry

    subroutine remove_nodes_geom ( geometry )

      type(geometry_t), intent(inout) :: geometry

      integer :: elem
      integer, allocatable, dimension(:) :: &
        tmp,  &! temporary allocatable
        nnrl   ! new local node numbering for curves and surfaces

      call move_alloc ( geometry%nodes, tmp )
      call remove_dbl ( tmp, geometry%nodes, nnrl )
      geometry%nnodes = size(geometry%nodes)
      do elem = 1, geometry%nelem
        geometry%topology(:,elem,1) = nnrl(geometry%topology(:,elem,1))
      end do

      deallocate ( tmp, nnrl )

    end subroutine remove_nodes_geom

  end subroutine mesh_change_remove_double_nodes_geometries


! Change the elshape of a basic mesh to a macro variant with linear subelements

  subroutine mesh_change_regular_to_macro ( mesh, groups, blends, curves, &
    surfaces, volumes )

    type(mesh_t), intent(inout) :: mesh
    logical, dimension(:), intent(in) :: groups, curves, surfaces, volumes
    logical, dimension(0:), intent(in) :: blends

    integer :: elshape, elgrp, newelshape, curve, surface, volume, m

!   bulk elements

    do elgrp = 1, mesh%nelgrp

      if ( .not. groups(elgrp) ) cycle

      elshape = mesh%element(elgrp)%elshape
      newelshape = 0

      call select_elshape_grp

      if ( newelshape > 0 .and. blends(0) ) &
                                      mesh%element(elgrp)%elshape = newelshape

      do m = 1, mesh%nblend

        if ( .not. blends(m) ) cycle

        elshape = mesh%element_blend(elgrp,m)%elshape
        newelshape = 0

        call select_elshape_grp

        if ( newelshape > 0 ) mesh%element_blend(elgrp,m)%elshape = newelshape

      end do

    end do

!   curves

    do curve = 1, mesh%ncurves

      if ( .not. curves(curve) ) cycle

      elshape = mesh%curves(curve)%element%elshape
      newelshape = 0

      call select_elshape_curve

      if ( newelshape > 0 .and. blends(0) ) &
                               mesh%curves(curve)%element%elshape = newelshape

      do m = 1, mesh%nblend

        if ( .not. blends(m) ) cycle

        elshape = mesh%curves(curve)%element_blend(m)%elshape
        newelshape = 0

        call select_elshape_curve

        if ( newelshape > 0 ) &
                     mesh%curves(curve)%element_blend(m)%elshape = newelshape

      end do

    end do

!   surfaces

    do surface = 1, mesh%nsurfaces

      if ( .not. surfaces(surface) ) cycle

      elshape = mesh%surfaces(surface)%element%elshape
      newelshape = 0

      call select_elshape_surface

      if ( newelshape > 0 .and. blends(0) ) &
                            mesh%surfaces(surface)%element%elshape = newelshape

      do m = 1, mesh%nblend

        if ( .not. blends(m) ) cycle

        elshape = mesh%surfaces(surface)%element_blend(m)%elshape
        newelshape = 0

        call select_elshape_surface

        if ( newelshape > 0 ) &
                 mesh%surfaces(surface)%element_blend(m)%elshape = newelshape

      end do

    end do

!   volumes

    do volume = 1, mesh%nvolumes

      if ( .not. volumes(volume) ) cycle

      elshape = mesh%volumes(volume)%element%elshape
      newelshape = 0

      call select_elshape_volume

      if ( newelshape > 0 .and. blends(0) ) &
                            mesh%volumes(volume)%element%elshape = newelshape

      do m = 1, mesh%nblend

        if ( .not. blends(m) ) cycle

        elshape = mesh%volumes(volume)%element_blend(m)%elshape
        newelshape = 0

        call select_elshape_volume

        if ( newelshape > 0 ) &
                     mesh%volumes(volume)%element_blend(m)%elshape = newelshape

      end do

    end do

  contains

    subroutine select_elshape_grp

      select case ( elshape )
      case(2) ! 3 node line element
        newelshape = 32
      case(4) ! 6 node triangle
        newelshape = 33
      case(6) ! 9 node quadrilateral
        newelshape = 34
      case(12) ! 10 node tetrahedron
        newelshape = 35
      case(14) ! 27 node hexahedron
        newelshape = 36
      case(101) ! line high-order element
        newelshape = 105
      case(102) ! quadrilateral high-order element
        newelshape = 106
      case(103) ! hexahedron high-order element
        newelshape = 107
      case default
        call errormsg_case_default ( 'select_elshape_grp', 'elshape', &
          int_value=elshape )
      end select

    end subroutine select_elshape_grp

    subroutine select_elshape_curve

      select case ( elshape )
      case(2) ! 3 node line element
        newelshape = 32
      case(101) ! line high-order element
        newelshape = 105
      case default
        call errormsg_case_default ( 'select_elshape_grp', 'elshape', &
          int_value=elshape )
      end select

    end subroutine select_elshape_curve

    subroutine select_elshape_surface

      select case ( elshape )
      case(4) ! 6 node triangle
        newelshape = 33
      case(6) ! 9 node quadrilateral
        newelshape = 34
      case(102) ! quadrilateral high-order element
        newelshape = 106
      case default
        call errormsg_case_default ( 'select_elshape_surface', 'elshape', &
          int_value=elshape )
      end select

    end subroutine select_elshape_surface

    subroutine select_elshape_volume

      select case ( elshape )
      case(12) ! 10 node tetrahedron
        newelshape = 35
      case(14) ! 27 node hexahedron
        newelshape = 36
      case(103) ! hexahedron high-order element
        newelshape = 107
      case default
        call errormsg_case_default ( 'select_elshape_volume', 'elshape', &
          int_value=elshape )
      end select

    end subroutine select_elshape_volume

  end subroutine mesh_change_regular_to_macro


! Change the elshape of a basic mesh to a macro variant with linear subelements

  subroutine mesh_change_macro_to_regular ( mesh, groups, blends, curves, &
    surfaces, volumes )

    type(mesh_t), intent(inout) :: mesh
    logical, dimension(:), intent(in) :: groups, curves, surfaces, volumes
    logical, dimension(0:), intent(in) :: blends

    integer :: elshape, elgrp, newelshape, curve, surface, volume, m

!   bulk elements

    do elgrp = 1, mesh%nelgrp

      if ( .not. groups(elgrp) ) cycle

      elshape = mesh%element(elgrp)%elshape
      newelshape = 0

      call select_elshape_grp

      if ( newelshape > 0 .and. blends(0) ) &
                                      mesh%element(elgrp)%elshape = newelshape

      do m = 1, mesh%nblend

        if ( .not. blends(m) ) cycle

        elshape = mesh%element_blend(elgrp,m)%elshape
        newelshape = 0

        call select_elshape_grp

        if ( newelshape > 0 ) mesh%element_blend(elgrp,m)%elshape = newelshape

      end do

    end do

!   curves

    do curve = 1, mesh%ncurves

      if ( .not. curves(curve) ) cycle

      elshape = mesh%curves(curve)%element%elshape
      newelshape = 0

      call select_elshape_curve

      if ( newelshape > 0 .and. blends(0) ) &
                               mesh%curves(curve)%element%elshape = newelshape

      do m = 1, mesh%nblend

        if ( .not. blends(m) ) cycle

        elshape = mesh%curves(curve)%element_blend(m)%elshape
        newelshape = 0

        call select_elshape_curve

        if ( newelshape > 0 ) &
                     mesh%curves(curve)%element_blend(m)%elshape = newelshape

      end do

    end do

!   surfaces

    do surface = 1, mesh%nsurfaces

      if ( .not. surfaces(surface) ) cycle

      elshape = mesh%surfaces(surface)%element%elshape
      newelshape = 0

      call select_elshape_surface

      if ( newelshape > 0 .and. blends(0) ) &
                            mesh%surfaces(surface)%element%elshape = newelshape

      do m = 1, mesh%nblend

        if ( .not. blends(m) ) cycle

        elshape = mesh%surfaces(surface)%element_blend(m)%elshape
        newelshape = 0

        call select_elshape_surface

        if ( newelshape > 0 ) &
                 mesh%surfaces(surface)%element_blend(m)%elshape = newelshape

      end do

    end do

!   volumes

    do volume = 1, mesh%nvolumes

      if ( .not. volumes(volume) ) cycle

      elshape = mesh%volumes(volume)%element%elshape
      newelshape = 0

      call select_elshape_volume

      if ( newelshape > 0 .and. blends(0) ) &
                            mesh%volumes(volume)%element%elshape = newelshape

      do m = 1, mesh%nblend

        if ( .not. blends(m) ) cycle

        elshape = mesh%volumes(volume)%element_blend(m)%elshape
        newelshape = 0

        call select_elshape_volume

        if ( newelshape > 0 ) &
                     mesh%volumes(volume)%element_blend(m)%elshape = newelshape

      end do

    end do

  contains

    subroutine select_elshape_grp

      select case ( elshape )
      case(32) ! 3 node line element
        newelshape = 2
      case(33) ! 6 node triangle
        newelshape = 4
      case(34) ! 9 node quadrilateral
        newelshape = 6
      case(35) ! 10 node tetrahedron
        newelshape = 12
      case(36) ! 27 node hexahedron
        newelshape = 14
      case(105) ! line high-order element
        newelshape = 101
      case(106) ! quadrilateral high-order element
        newelshape = 102
      case(107) ! hexahedron high-order element
        newelshape = 103
      case default
        call errormsg_case_default ( 'select_elshape_grp', 'elshape', &
          int_value=elshape )
      end select

    end subroutine select_elshape_grp

    subroutine select_elshape_curve

      select case ( elshape )
      case(32) ! 3 node line element
        newelshape = 2
      case(105) ! line high-order element
        newelshape = 101
      case default
        call errormsg_case_default ( 'select_elshape_curve', 'elshape', &
          int_value=elshape )
      end select

    end subroutine select_elshape_curve

    subroutine select_elshape_surface

      select case ( elshape )
      case(33) ! 6 node triangle
        newelshape = 4
      case(34) ! 9 node quadrilateral
        newelshape = 6
      case(106) ! quadrilateral high-order element
        newelshape = 102
      case default
        call errormsg_case_default ( 'select_elshape_surface', 'elshape', &
          int_value=elshape )
      end select

    end subroutine select_elshape_surface

    subroutine select_elshape_volume

      select case ( elshape )
      case(35) ! 10 node tetrahedron
        newelshape = 12
      case(36) ! 27 node hexahedron
        newelshape = 14
      case(107) ! hexahedron high-order element
        newelshape = 103
      case default
        call errormsg_case_default ( 'select_elshape_volume', 'elshape', &
          int_value=elshape )
      end select

    end subroutine select_elshape_volume

  end subroutine mesh_change_macro_to_regular


! Generate an incomplete geometry, having no nodes and topology.
! Fill components nodes and topology to complete the basic geometry.
! Note: fill nodes and topology(:,:,1) first and then topology(:,:,2) using
!   do elem = 1, nelem
!     topology(:,elem,2) = nodes(topology(:,elem,1))
!   end do
! Don't use the geometry until these components have been filled!

  subroutine geometry_skeleton ( geometry, nnodes, nelem, elshape, ndim, &
    callname )

    type(geometry_t), intent(inout) :: geometry

!   number of nodes in the geometry
    integer, intent(in) :: nnodes

!   number of elements in the geometry
    integer, intent(in) :: nelem

!   element types in the geometry. See element_t for the possibilities
    integer, intent(in) :: elshape

!   space dimension for line, quadrilateral and triangle elements
!   default=1 for line elements
!   default=2 for quadrilateral and triangle elements
    integer, intent(in), optional :: ndim

!   name of calling routine
    character (len=*), intent(in), optional :: callname


    character (len=13) :: globalshape
    character (len=:), allocatable :: lcallname
    integer :: inpelm ! number of nodal points per element
    integer :: ndiml
    integer, dimension (53) :: inpelm_ar

    lcallname = set_optional ( variable=callname, default='geometry_skeleton' )

!   test geometry

    if ( allocated(geometry%topology) ) then
      write(*,'(/3a/2(a/))') 'Error ', lcallname, ':', &
        ' geometry is not empty', &
        ' either use a new geometry_t variable or delete old geometry '
      stop
    end if

    if ( nnodes < 0 ) then
      write(*,'(/3a/a/)') 'Error ', lcallname, ':', &
        ' nnodes is a negative number'
      stop
    end if
    if ( nelem < 0 ) then
      write(*,'(/3a/a/)') 'Error ', lcallname, ':', &
        ' nelem is a negative number'
      stop
    end if

!   set globalshape of elements

    select case ( elshape )
      case (1,2,32); globalshape = 'line'
      case (3,4,7,10,33); globalshape = 'triangle'
      case (5,6,9,30,34); globalshape = 'quadrilateral'
      case (11,12,15,16,18,35); globalshape = 'tetrahedron'
      case (13,14,17,31,36); globalshape = 'hexahedron'
      case (41,42,43); globalshape = 'prism'
      case (51,52,53); globalshape = 'pyramid'
      case default
        write(*,'(/a/a,i0/)') 'Error geometry_skeleton:', &
          ' invalid elshape = ', elshape
        stop
    end select

!   set inpelm of elements

    inpelm_ar = 0
    inpelm_ar(1:18) = &
      [ 2, 3, 3, 6, 4, 9, 7, 0, 5, 4, 4, 10, 8, 27, 14, 15, 9, 5 ]
    inpelm_ar(30:36) = [ 8, 20, 3, 6, 9, 10, 27 ]
    inpelm_ar(41:43) = [ 6, 15, 18 ]
    inpelm_ar(51:53) = [ 5, 13, 14 ]

    inpelm = inpelm_ar(elshape)

!   set ndim

    if ( globalshape == 'line' ) then
      ndiml = set_optional ( variable=ndim, default=1 )
      if ( ndiml < 1 .or. ndiml > 3 ) then
        write(*,'(/3a/a/)') 'Error ', lcallname, ':', &
          ' ndim < 1 or ndim > 3 '
        stop
      end if
    else if ( globalshape == 'triangle' .or. &
              globalshape == 'quadrilateral' )  then
      ndiml = set_optional ( variable=ndim, default=2 )
      if ( ndiml < 2 .or. ndiml > 3 ) then
        write(*,'(/3a/a/)') 'Error ', lcallname, ':', &
          ' ndim < 2 or ndim > 3 for triangles or quadrilaterals'
        stop
      end if
    else
      ndiml = 3
    end if

!   create geometry

    geometry%ndim   = ndiml
    geometry%nelem  = nelem
    geometry%nnodes = nnodes
    geometry%nblend = 0

!   allocate memory for geometry

    allocate( geometry%element_blend(0) )
    geometry%nnodes_blend = [0,nnodes]
    allocate( geometry%topology(inpelm,nelem,2) )
    allocate( geometry%nodes(nnodes) )

!   element type

    geometry%element%elshape = elshape
    geometry%element%numnod = inpelm
    geometry%element%ndim = ndiml
    geometry%element%globalshape = globalshape

!   topology (set to zero)

    geometry%topology = 0

!   nodes (set to zero)

    geometry%nodes = 0

  end subroutine geometry_skeleton


! Generate an incomplete mesh, having no coordinates and topology
! Don't use this mesh until the coordinates and topology have been filled!

  subroutine mesh_skeleton ( mesh, nnodes, nelem, elshape, ndim, callname )

    type(mesh_t), intent(inout) :: mesh

!   number of nodes in the mesh
    integer, intent(in) :: nnodes

!   number of elements in the mesh
    integer, intent(in) :: nelem

!   element types in the mesh. See element_t for the possibilities
    integer, intent(in) :: elshape

!   space dimension for line, quadrilateral and triangle elements
!   default=1 for line elements
!   default=2 for quadrilateral and triangle elements
    integer, intent(in), optional :: ndim

!   name of calling routine
    character (len=*), intent(in), optional :: callname


    character (len=13) :: globalshape
    character (len=:), allocatable :: lcallname
    integer :: inpelm ! number of nodal points per element
    integer :: ndiml
    integer, dimension (53) :: inpelm_ar

    lcallname = set_optional ( variable=callname, default='mesh_skeleton' )

!   test mesh

    if ( mesh%meshgen ) then
      write(*,'(/3a/2(a/))') 'Error ', lcallname, ':', &
        ' mesh is not empty', &
        ' either use a new mesh_t variable or delete old mesh '
      stop
    end if

    if ( nnodes < 0 ) then
      write(*,'(/3a/a/)') 'Error ', lcallname, ':', &
        ' nnodes is a negative number'
      stop
    end if
    if ( nelem < 0 ) then
      write(*,'(/3a/a/)') 'Error ', lcallname, ':', &
        ' nelem is a negative number'
      stop
    end if

!   set globalshape of elements

    select case ( elshape )
      case (1,2,32); globalshape = 'line'
      case (3,4,7,10,33); globalshape = 'triangle'
      case (5,6,9,30,34); globalshape = 'quadrilateral'
      case (11,12,15,16,18,35); globalshape = 'tetrahedron'
      case (13,14,17,31,36); globalshape = 'hexahedron'
      case (41,42,43); globalshape = 'prism'
      case (51,52,53); globalshape = 'pyramid'
      case default
        write(*,'(/a/a,i0/)') 'Error mesh_skeleton:', &
          ' invalid elshape = ', elshape
        stop
    end select

!   set inpelm of elements

    inpelm_ar = 0
    inpelm_ar(1:18) = &
      [ 2, 3, 3, 6, 4, 9, 7, 0, 5, 4, 4, 10, 8, 27, 14, 15, 9, 5 ]
    inpelm_ar(30:36) = [ 8, 20, 3, 6, 9, 10, 27 ]
    inpelm_ar(41:43) = [ 6, 15, 18 ]
    inpelm_ar(51:53) = [ 5, 13, 14 ]

    inpelm = inpelm_ar(elshape)

!   set ndim

    if ( globalshape == 'line' ) then
      ndiml = set_optional ( variable=ndim, default=1 )
      if ( ndiml < 1 .or. ndiml > 3 ) then
        write(*,'(/3a/a/)') 'Error ', lcallname, ':', &
          ' ndim < 1 or ndim > 3 '
        stop
      end if
    else if ( globalshape == 'triangle' .or. &
              globalshape == 'quadrilateral' )  then
      ndiml = set_optional ( variable=ndim, default=2 )
      if ( ndiml < 2 .or. ndiml > 3 ) then
        write(*,'(/3a/a/)') 'Error ', lcallname, ':', &
          ' ndim < 2 or ndim > 3 for triangles or quadrilaterals'
        stop
      end if
    else
      ndiml = 3
    end if

!   create mesh

    mesh%ndim   = ndiml
    mesh%nelem  = nelem
    mesh%nnodes = nnodes
    mesh%nelgrp = 1
    mesh%nblend = 0

!   allocate memory for mesh

    allocate( mesh%element(1) )
    allocate( mesh%element_blend(1,0) )
    mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel(1) )
    allocate( mesh%topology(1) )
    allocate( mesh%topology(1)%a(inpelm,nelem) )
    allocate( mesh%coor(nnodes,ndiml) )

!   element type

    mesh%element(1)%elshape = elshape
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = ndiml
    mesh%element(1)%globalshape = globalshape

!   topology (set to zero)

    mesh%grpnumel(1) = mesh%nelem

    mesh%topology(1)%a = 0

!   coordinates (set to zero)

    mesh%coor = 0

!   allocate memory for points, curves, surfaces and volumes

    allocate( mesh%points(MAXPOINTS), mesh%curves(MAXCURVES), &
      mesh%surfaces(MAXSURFACES), mesh%volumes(MAXVOLUMES) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   basic mesh has been generated

    mesh%meshgen = .true.

  end subroutine mesh_skeleton

end module meshgen_extra_m
