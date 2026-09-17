
! Copyright (C) 2004-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Meshgeneration for remaining "parts" of the mesh

module meshgen_parts_m

  use glob_defs_m
  use mesh_m, only: mesh_t, connectdata_t, element_t
  use meshgen_objects_m
  use meshgen_construct_m, only: add_to_mesh
  use misc_m, only: sort

  implicit none

contains


! Compute other parts of the mesh from the topology

  subroutine fill_mesh_parts ( mesh, onlyshapenod, blend )

    type(mesh_t), intent(inout) :: mesh

!   optional argument to indicate that only numshapenod nodes are
!   used in the geometrical shape of the element for finding the reference
!   coordinates of objects. See heading of routine find_refcoor_objects for more
!   information.
    logical, intent(in), optional :: onlyshapenod

!   optional argument to indicate which blend mesh is being used for:
!    - filling sidelem, gluepoints, edges and faces for geometries.
!    - finding reference coordinates for intersecting objects
!   default = 0 (= main mesh)
    integer, intent(in), optional :: blend


    if ( .not. mesh%meshgen ) then
      write(*,'(/a/)') &
        'Error fill_mesh_parts: no mesh present in structure mesh '
      stop
    end if

    if ( mesh%meshparts ) then
      write(*,'(/2(a/))') &
        'Error: fill_mesh_parts can only be called once for a mesh.', &
        ' Use mesh_convert to convert to a basic mesh first.'
      stop
    end if

    call check_nblend ( mesh, 'fill_mesh_parts' )

    if ( present(blend) ) then
      if ( blend < 0 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error fill_mesh_parts: invalid argument ', &
          '  blend < 0 or larger than ', mesh%nblend
        stop
      end if
    end if

!   Fill part of the mesh related to element types

    call fill_mesh_parts_elements ( mesh )

!   Fill part of the mesh related to elements connected to nodes

    call fill_mesh_parts_nodelem ( mesh )

!   Fill part of the mesh related to nodes connected to nodes

    call fill_mesh_parts_nodnod ( mesh )

!   Fill part of the mesh related to elements connected to sides

    call fill_mesh_parts_sidelem ( mesh, blend )

!   Fill part of the mesh related to elements connected to sides (gluepoints)

    call fill_mesh_parts_sidelem_glue ( mesh, blend )

!   Fill part of the mesh related to curves, surfaces and volumes

    call fill_mesh_parts_geometries ( mesh%curves(1:mesh%ncurves) )
    call fill_mesh_parts_geometries ( mesh%surfaces(1:mesh%nsurfaces) )
    call fill_mesh_parts_geometries ( mesh%volumes(1:mesh%nvolumes) )
    call fill_mesh_parts_geometries_edges ( mesh, &
                                      mesh%curves(1:mesh%ncurves), blend )
    call fill_mesh_parts_geometries_faces ( mesh, &
                                      mesh%surfaces(1:mesh%nsurfaces), blend )
    call fill_mesh_parts_geometries_elements ( mesh, &
                                      mesh%volumes(1:mesh%nvolumes) )

!   Fill part of the mesh related to blocks

    call fill_mesh_parts_blocks ( mesh )

!   Fill part of the mesh related to objects

    call fill_mesh_parts_objects ( mesh, onlyshapenod=onlyshapenod, &
                                                                 blend=blend )

    mesh%meshparts = .true.
    mesh%curves(:mesh%ncurves)%meshparts = .true.
    mesh%surfaces(:mesh%nsurfaces)%meshparts = .true.
    mesh%objects(:mesh%nobjects)%meshparts = .true.

  end subroutine fill_mesh_parts


! Fill part of the mesh related to elements

  subroutine fill_mesh_parts_elements ( mesh )

    type(mesh_t), intent(inout) :: mesh

!   Fill element type info
!   NOTE: elshape, elshape_blend, globalshape, numnod and ndim must have been
!   filled already.
!   For high-order elements also p must have been filled already.

    integer :: elgrp, m

!   fill elnumnod

    allocate ( mesh%elnumnod ( mesh%nelgrp ) )

    mesh%elnumnod = mesh%element(:)%numnod

!   fill rest of element info

    do elgrp = 1, mesh%nelgrp
      call fill_element ( mesh%element(elgrp) )
    end do

!   Blended meshes

    do m = 1, mesh%nblend

      mesh%elnumnod = mesh%elnumnod + mesh%element_blend(:,m)%numnod

!     fill rest of element info

      do elgrp = 1, mesh%nelgrp
        call fill_element ( mesh%element_blend(elgrp,m) )
      end do

    end do

    mesh%maxelnumnod = maxval ( mesh%elnumnod )

    allocate ( mesh%numnodtop(mesh%nelgrp,mesh%nblend+2) )

    mesh%numnodtop(:,1) = 0
    mesh%numnodtop(:,2) = mesh%element(:)%numnod
    do m = 1, mesh%nblend
      mesh%numnodtop(:,m+2) = &
               mesh%numnodtop(:,m+1) + mesh%element_blend(:,m)%numnod
    end do

  end subroutine fill_mesh_parts_elements


! Fill part of the mesh related to elements connected to nodes

  subroutine fill_mesh_parts_nodelem ( mesh )

    type(mesh_t), intent(inout) :: mesh

    integer :: elgrp, elem, node, nodenr
    integer, dimension(:), allocatable :: work

!   mesh

!   fill nodnumel

    allocate ( mesh%nodnumel ( mesh%nnodes + 1 ) )

    mesh%nodnumel = 0

    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)
        do node = 1, mesh%elnumnod(elgrp)
          nodenr = mesh%topology(elgrp)%a(node,elem)
          mesh%nodnumel(nodenr+1) = mesh%nodnumel(nodenr+1) + 1
        end do
      end do
    end do

!   accumulate nodnumel

    do node = 2, mesh%nnodes
      mesh%nodnumel(node+1) = mesh%nodnumel(node) + mesh%nodnumel(node+1)
    end do

!   fill nodelem

    allocate ( mesh%nodelem( mesh%nodnumel ( mesh%nnodes + 1 ), 2 ) )
    allocate ( work( mesh%nnodes ) )

    work  = 0

    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)
        do node = 1, mesh%elnumnod(elgrp)
          nodenr = mesh%topology(elgrp)%a(node,elem)
          work(nodenr) = work(nodenr) + 1
          mesh%nodelem( mesh%nodnumel(nodenr) + work(nodenr), 1 ) = elgrp
          mesh%nodelem( mesh%nodnumel(nodenr) + work(nodenr), 2 ) = elem
        end do
      end do
    end do

    deallocate ( work )

  end subroutine fill_mesh_parts_nodelem


! Fill part of the mesh related to nodes connected to nodes

  subroutine fill_mesh_parts_nodnod ( mesh )

    type(mesh_t), intent(inout) :: mesh

    integer :: elgrp, elem, node, nodenr, elm
    integer :: nn, ipwork, nod
    integer, dimension(:), allocatable :: work, work1


!   mesh

!   fill nodnumnod

    allocate ( mesh%nodnumnod ( mesh%nnodes + 1 ) )

    mesh%nodnumnod(1) = 0

!   first create enough workspace

    allocate ( work( maxval(mesh%elnumnod) * size(mesh%nodelem,1) ) )
    allocate ( work1( mesh%nnodes ) )

    ipwork = 0   ! pointer in array work
    work1  = 0

    do nodenr = 1, mesh%nnodes

      work1(nodenr) = 1  ! Exclude node itself
      nn = 0

!     fill work with nodal points from elements connected

      do elm = 1, mesh%nodnumel(nodenr+1)-mesh%nodnumel(nodenr)

        elgrp = mesh%nodelem( mesh%nodnumel(nodenr) + elm, 1 )
        elem  = mesh%nodelem( mesh%nodnumel(nodenr) + elm, 2 )

        do node = 1, mesh%elnumnod(elgrp)

          nod = mesh%topology(elgrp)%a(node,elem)

          if ( work1(nod) /= 0 ) cycle  ! node already used

          work1(nod) = 1 ! new node
          nn = nn + 1
          work(ipwork+nn) = nod

        end do

      end do

!     put work1 back to zero

      work1(nodenr) = 0
      work1(work(ipwork+1:ipwork+nn)) = 0

!     fill nodnumnod

      mesh%nodnumnod(nodenr+1) = nn
      ipwork = ipwork + nn

    end do

    deallocate ( work1 )

    mesh%maxnodnumnod = maxval(mesh%nodnumnod)

!   accumulate nodnumnod

    do node = 2, mesh%nnodes
      mesh%nodnumnod(node+1) = mesh%nodnumnod(node) + mesh%nodnumnod(node+1)
    end do

!   fill nonnod

    allocate ( mesh%nodnod( mesh%nodnumnod(mesh%nnodes+1) ) )

    mesh%nodnod = work(1:mesh%nodnumnod(mesh%nnodes+1))

    deallocate ( work )

  end subroutine fill_mesh_parts_nodnod


! Fill part of the mesh related to elements connected to sides

  subroutine fill_mesh_parts_sidelem ( mesh, blend )

    use limits_m, &
         only: TEST_FOR_MULTIPLE_SIDELEM, WARN_ON_TEST_FOR_MULTIPLE_SIDELEM, &
               WARN_ON_TEST_FOR_VARYING_SIDES

    type(mesh_t), intent(inout) :: mesh
    integer, intent(in), optional :: blend

    integer :: elgrp, elem, side, el, elgrpnr, elemnr
    integer :: node1, nodenr1
    integer :: node, nodenr, numel, vert, sidenr
    integer, dimension(:,:), allocatable :: work
    logical :: found
    type(element_t), dimension(mesh%nelgrp) :: element
    integer, dimension(mesh%nelgrp) :: topoffset

!   select blend mesh

    topoffset = 0

    if ( present(blend) ) then
      if ( blend == 0 ) then
        element = mesh%element
      else
        element = mesh%element_blend(:,blend)
        topoffset = mesh%numnodtop(:,blend+1)
      end if
    else
      element = mesh%element
    end if

!   fill sidelem

    allocate ( mesh%sidelem ( mesh%nelgrp ) )
    do elgrp = 1, mesh%nelgrp
      allocate ( mesh%sidelem(elgrp)%a(element(elgrp)%numsides, &
                                       mesh%grpnumel(elgrp),4) )
    end do

!   test special cases

    if ( any ( element(:)%globalshape == 'prism' ) .or. &
         any ( element(:)%globalshape == 'pyramid' ) ) then

!     these meshes can have varying side shapes and should be excluded for now

      if ( WARN_ON_TEST_FOR_VARYING_SIDES ) then
        write(*,'(/12(a/))') &
          'Warning fill_mesh_parts_sidelem: mesh contains ', &
          ' prism and/or pyramid shaped element.', &
          ' Computing mesh%sidelem is not implemented ', &
          ' for elements having varying side shapes.', &
          ' Set (use limits_m module) ', &
          '   WARN_ON_TEST_FOR_VARYING_SIDES = .false. ', &
          ' to stop this warning from appearing.'
      end if

      do elgrp = 1, mesh%nelgrp
        mesh%sidelem(elgrp)%a = 0
      end do

      return

    end if

    if ( TEST_FOR_MULTIPLE_SIDELEM .and. ( &
         ( any ( element(:)%globalshape == 'line' ) &
                                            .and. mesh%ndim > 1 ) .or. &
         ( any ( element(:)%globalshape == 'triangle' ) &
                                           .and. mesh%ndim == 3 ) .or. &
         ( any ( element(:)%globalshape == 'quadrilateral' ) &
                                           .and. mesh%ndim == 3 ) ) ) then

!     these meshes can have multiple connected sides and should be excluded

      if ( WARN_ON_TEST_FOR_MULTIPLE_SIDELEM ) then
        write(*,'(/12(a/))') &
          'Warning fill_mesh_parts_sidelem: mesh contains ', &
          '  1) line elements for ndim=2,3 or ', &
          '  2) triangular or quadrilateral elements for ndim=3 ', &
          ' and therefore mesh%sidelem is not computed because ', &
          ' multiple connected elements are possible. This can have ', &
          ' some unexpected side effects, such as mesh plots in ', &
          ' data plots with figplot. Set (use limits_m module) ', &
          '   WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = .false. ', &
          ' to stop this warning from appearing or ', &
          '   TEST_FOR_MULTIPLE_SIDELEM = .false. ', &
          ' to skip the test (and compute mesh%sidelem anyway), ', &
          ' if your mesh does not have multiple connected elements.'
      end if

      do elgrp = 1, mesh%nelgrp
        mesh%sidelem(elgrp)%a = 0
      end do

      return

    end if

    allocate ( work( mesh%nelgrp, maxval(mesh%grpnumel) ) )

    work = 0

!   start loop

    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)
        do side = 1, element(elgrp)%numsides

!         fill work for elements connected to the nodes on the side

          do vert = 1, element(elgrp)%sidnumvert
            node = element(elgrp)%sidvert(vert,side) + topoffset(elgrp)
            nodenr = mesh%topology(elgrp)%a(node,elem)
            numel = mesh%nodnumel(nodenr+1)-mesh%nodnumel(nodenr)
            do el = 1, numel
              elgrpnr = mesh%nodelem( mesh%nodnumel(nodenr) + el, 1 )
              elemnr  = mesh%nodelem( mesh%nodnumel(nodenr) + el, 2 )
              work(elgrpnr,elemnr) = work(elgrpnr,elemnr) + 1
            end do
          end do

!         check whether an element shows up sidnumvert times (node1 only)

          node = element(elgrp)%sidvert(1,side) + topoffset(elgrp)
          nodenr = mesh%topology(elgrp)%a(node,elem)
          numel = mesh%nodnumel(nodenr+1)-mesh%nodnumel(nodenr)
          found = .false.
          do el = 1, numel
            elgrpnr = mesh%nodelem( mesh%nodnumel(nodenr) + el, 1 )
            elemnr  = mesh%nodelem( mesh%nodnumel(nodenr) + el, 2 )
            if ( work( elgrpnr, elemnr ) == &
                       element(elgrp)%sidnumvert .and.  &
                 ( elemnr /= elem .or. elgrpnr /= elgrp ) ) then
!             element found (excluding current element)
              found = .true.
              exit
            end if
          end do

!         fill sidelem for this element

          if ( found ) then
            mesh%sidelem(elgrp)%a(side,elem,:) = [ elgrpnr, elemnr, 0, 0 ]
          else
            mesh%sidelem(elgrp)%a(side,elem,:) = 0
          end if

!         set work to zero again

          do vert = 1, element(elgrp)%sidnumvert
            node = element(elgrp)%sidvert(vert,side) + topoffset(elgrp)
            nodenr = mesh%topology(elgrp)%a(node,elem)
            numel = mesh%nodnumel(nodenr+1)-mesh%nodnumel(nodenr)
            do el = 1, numel
              elgrpnr = mesh%nodelem( mesh%nodnumel(nodenr) + el, 1 )
              elemnr  = mesh%nodelem( mesh%nodnumel(nodenr) + el, 2 )
              work(elgrpnr,elemnr) = 0
            end do
          end do

        end do
      end do
    end do

!   Fill the side number and orientation of the adjacent elements

    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)
        do side = 1, element(elgrp)%numsides

          elgrpnr = mesh%sidelem(elgrp)%a(side,elem,1)
          elemnr  = mesh%sidelem(elgrp)%a(side,elem,2)

          if ( elgrpnr /= 0 ) then

!           check ajacent element for sidenr
            do sidenr = 1, element(elgrpnr)%numsides
              if ( elgrp == mesh%sidelem(elgrpnr)%a(sidenr,elemnr,1) .and. &
                   elem == mesh%sidelem(elgrpnr)%a(sidenr,elemnr,2) ) then
                 mesh%sidelem(elgrp)%a(side,elem,3) = sidenr
                 exit
              end if
            end do

            if ( mesh%sidelem(elgrp)%a(side,elem,3) == 0 ) then
              write(*,'(/a/)') &
                'Error fill_mesh_parts_sidelem: inconsistent mesh topology'
                stop
            end if

!           determine orientation of adjacent element
            node1 = element(elgrp)%sidvert(1,side) + topoffset(elgrp)
            nodenr1 = mesh%topology(elgrp)%a(node1,elem)
            do vert = 1, element(elgrpnr)%sidnumvert
              node = element(elgrpnr)%sidvert(vert,sidenr) + topoffset(elgrpnr)
              nodenr = mesh%topology(elgrpnr)%a(node,elemnr)
              if ( nodenr == nodenr1 ) then
                mesh%sidelem(elgrp)%a(side,elem,4) = vert
                exit
              end if
            end do

          end if

        end do
      end do
    end do

    deallocate ( work )

  end subroutine fill_mesh_parts_sidelem


! Fill part of the mesh related to elements connected to sides (add gluepoints)

  subroutine fill_mesh_parts_sidelem_glue ( mesh, blend )

    type(mesh_t), intent(inout) :: mesh
    integer, intent(in), optional :: blend

    integer :: elgrp, elem, side, el, elgrpnr, elemnr
    integer :: node1, nodenr1
    integer :: node, nodenr, numel, vert, sidenr
    integer :: pnt, pnt1, pnt2, gpnt, gnode, gnodenr, gnodenr1

!   this routine adds the effect of gluepoints to sidelem

!   work1: number of nodes connected via gluepoints (accumulated)
!     length = mesh%nnodes+1
!     work1(1) = 0,
!     number of nodes for nodal point n is
!     work1(n+1) - work1(n)
!   work2: number of nodes connected via gluepoints
!     length = mesh%nnodes
!     number of nodes for nodal point n is work2(n)
!   work3: nodes connected to nodal points via gluepoints
!     number of nodes for nodal point n is numnod = work1(n+1) - work1(n)
!     or work2(n)
!     the nodal numbers numbers are stored in
!     work3( work1(n) + i ), i = 1, numnod
    integer, dimension(:), allocatable :: work1, work2, work3
    integer, dimension(:,:), allocatable :: work
    logical :: found
    type(element_t), dimension(mesh%nelgrp) :: element
    integer, dimension(mesh%nelgrp) :: topoffset

    if ( mesh%ngluepoints == 0 ) return

!   sidelem exists?

    if ( .not. allocated(mesh%sidelem) ) then
      write(*,'(/a/)') &
        'Error fill_mesh_parts_sidelem_glue: mesh%sidelem not allocated.'
      stop
    end if

!   select blend mesh

    topoffset = 0

    if ( present(blend) ) then
      if ( blend == 0 ) then
        element = mesh%element
      else
        element = mesh%element_blend(:,blend)
        topoffset = mesh%numnodtop(:,blend+1)
      end if
    else
      element = mesh%element
    end if

!   take glue points into account: create work arrays

    allocate ( work( mesh%nelgrp, maxval(mesh%grpnumel) ) )

    work = 0

    allocate (work1(mesh%nnodes+1))
    work1 = 0

!   count connections for each node

    do gpnt = 1, mesh%ngluepoints
      work1(mesh%gluepoints(gpnt,1)+1) = work1(mesh%gluepoints(gpnt,1)+1) + 1
      work1(mesh%gluepoints(gpnt,2)+1) = work1(mesh%gluepoints(gpnt,2)+1) + 1
    end do

!   accumulate and fill nodes

    do pnt = 1, mesh%nnodes
      work1(pnt+1) = work1(pnt) + work1(pnt+1)
    end do

    allocate (work2(mesh%nnodes))
    work2 = 0

    allocate ( work3(work1(mesh%nnodes+1)) )

!   generate work2 (number of nodes) and work3 (the node numbers connected)

    do gpnt = 1, mesh%ngluepoints
!     each gluepoint generates two connections: pnt1 -> pnt2, pnt2 -> pnt1
      pnt1 = mesh%gluepoints(gpnt,1)
      work2(pnt1) = work2(pnt1) + 1
      pnt2 = mesh%gluepoints(gpnt,2)
      work2(pnt2) = work2(pnt2) + 1
      work3(work1(pnt1)+work2(pnt1)) = pnt2
      work3(work1(pnt2)+work2(pnt2)) = pnt1
    end do


!   start loop

    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)
        do side = 1, element(elgrp)%numsides

!         fill work for elements connected to the glue nodes on the side

          elgrpnr = mesh%sidelem(elgrp)%a(side,elem,1)

          if ( elgrpnr == 0 ) then

!           only add to side that doesn't have a connection yet
!           check possible gluepoints to add

            do vert = 1, element(elgrp)%sidnumvert
              node = element(elgrp)%sidvert(vert,side) + topoffset(elgrp)
              nodenr = mesh%topology(elgrp)%a(node,elem)
              do gnode = 1, work2(nodenr) ! check all possible connections
                gnodenr = work3(work1(nodenr)+gnode) ! glue node found
                numel = mesh%nodnumel(gnodenr+1)-mesh%nodnumel(gnodenr)
                do el = 1, numel
                  elgrpnr = mesh%nodelem( mesh%nodnumel(gnodenr) + el, 1 )
                  elemnr  = mesh%nodelem( mesh%nodnumel(gnodenr) + el, 2 )
                  work(elgrpnr,elemnr) = work(elgrpnr,elemnr) + 1
                end do
              end do
            end do

!           check whether an element shows up sidnumvert times (node1 only)

            node = element(elgrp)%sidvert(1,side) + topoffset(elgrp)
            nodenr = mesh%topology(elgrp)%a(node,elem)
            found = .false.
  gluenode: do gnode = 1, work2(nodenr)
              gnodenr = work3(work1(nodenr)+gnode)
              numel = mesh%nodnumel(gnodenr+1)-mesh%nodnumel(gnodenr)
              do el = 1, numel
                elgrpnr = mesh%nodelem( mesh%nodnumel(gnodenr) + el, 1 )
                elemnr  = mesh%nodelem( mesh%nodnumel(gnodenr) + el, 2 )
                if ( elemnr == elem .and. elgrpnr == elgrp ) then
                  write(*,'(/2a/2(a,i0)/)') &
                    'Error fill_mesh_parts_sidelem_glue:', &
                    ' element connected to itself:', &
                    ' elementnr = ', elem, ' groupnr = ', elgrp
                  stop
                end if
                if ( work( elgrpnr, elemnr ) == &
                           element(elgrp)%sidnumvert ) then
!                 element found
                  found = .true.
                  exit gluenode
                end if
              end do
            end do gluenode

!           fill sidelem for this element

            if ( found ) then
!             indicate newly found with ivert=-1
              mesh%sidelem(elgrp)%a(side,elem,:) = [ elgrpnr, elemnr, 0, -1 ]
            else
              mesh%sidelem(elgrp)%a(side,elem,:) = 0
            end if

!           set work to zero again

            do vert = 1, element(elgrp)%sidnumvert
              node = element(elgrp)%sidvert(vert,side) + topoffset(elgrp)
              nodenr = mesh%topology(elgrp)%a(node,elem)
              do gnode = 1, work2(nodenr)
                gnodenr = work3(work1(nodenr)+gnode)
                numel = mesh%nodnumel(gnodenr+1)-mesh%nodnumel(gnodenr)
                do el = 1, numel
                  elgrpnr = mesh%nodelem( mesh%nodnumel(gnodenr) + el, 1 )
                  elemnr  = mesh%nodelem( mesh%nodnumel(gnodenr) + el, 2 )
                  work(elgrpnr,elemnr) = 0
                end do
              end do
            end do

          end if

        end do
      end do
    end do

    deallocate ( work )

!   Fill the side number and orientation of the adjacent elements

    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)
        do side = 1, element(elgrp)%numsides

          elgrpnr = mesh%sidelem(elgrp)%a(side,elem,1)
          elemnr  = mesh%sidelem(elgrp)%a(side,elem,2)

          if ( mesh%sidelem(elgrp)%a(side,elem,4) == -1 ) then

!           check ajacent element for sidenr
!           only check elements that have vert < 0
            do sidenr = 1, element(elgrpnr)%numsides
              if ( elgrp == mesh%sidelem(elgrpnr)%a(sidenr,elemnr,1) .and. &
                   elem == mesh%sidelem(elgrpnr)%a(sidenr,elemnr,2)  .and. &
                   mesh%sidelem(elgrpnr)%a(sidenr,elemnr,4) < 0 ) then
                 mesh%sidelem(elgrp)%a(side,elem,3) = sidenr
                 exit
              end if
            end do

!           determine orientation of adjacent element
            node1 = element(elgrp)%sidvert(1,side) + topoffset(elgrp)
            nodenr1 = mesh%topology(elgrp)%a(node1,elem)
 gluenode1: do gnode = 1, work2(nodenr1) ! check all connecting nodes
              gnodenr1 = work3(work1(nodenr1)+gnode)
              do vert = 1, element(elgrpnr)%sidnumvert
                node = element(elgrpnr)%sidvert(vert,sidenr) + &
                                                            topoffset(elgrpnr)
                nodenr = mesh%topology(elgrpnr)%a(node,elemnr)
                if ( nodenr == gnodenr1 ) then
                  mesh%sidelem(elgrp)%a(side,elem,4) = - vert
                  exit gluenode1
                end if
              end do
            end do gluenode1

          end if

        end do
      end do
    end do

   deallocate( work1, work2, work3 )

  end subroutine fill_mesh_parts_sidelem_glue


! Fill part of the mesh related to curves, surfaces and volumes (geometries)

  subroutine fill_mesh_parts_geometries ( geometries )

    type(geometry_t), dimension(:), intent(inout) :: geometries

    integer :: geom, maxv, ngeom, m
    integer, dimension(:), allocatable :: work, work1

    ngeom = size(geometries)

!   fill element info

    do geom = 1, ngeom

      geometries(geom)%elnumnod = geometries(geom)%element%numnod

!     Main mesh

      call fill_element ( geometries(geom)%element )

!     Blended meshes

      do m = 1, geometries(geom)%nblend

        geometries(geom)%elnumnod = geometries(geom)%elnumnod &
                        + geometries(geom)%element_blend(m)%numnod

!       fill rest of element info

        call fill_element ( geometries(geom)%element_blend(m) )

      end do

!     fill numnodtop

      allocate ( geometries(geom)%numnodtop(geometries(geom)%nblend+2) )

      geometries(geom)%numnodtop(1) = 0
      geometries(geom)%numnodtop(2) = geometries(geom)%element%numnod
      do m = 1, geometries(geom)%nblend
        geometries(geom)%numnodtop(m+2) = &
              geometries(geom)%numnodtop(m+1) &
                    + geometries(geom)%element_blend(m)%numnod
      end do

    end do

!   fill nodnumel and nodelem (of nc) for geometries

!   first create enough workspace

    allocate ( work( maxval(geometries(:)%nnodes) ) )

    do geom = 1, ngeom

      call fill_nodelem_nc ( geometries(geom)%nnodes, &
        geometries(geom)%nelem, geometries(geom)%elnumnod, &
        geometries(geom)%topology(:,:,1), geometries(geom)%nc, work )

    end do

    deallocate ( work )

!   fill nodnumnod and nodnod (of nc) for geom

!   first create enough workspace

    maxv = 0
    do geom = 1, ngeom
       maxv = max(size(geometries(geom)%nc%nodelem),maxv)
    end do

    allocate ( work( maxval(geometries(:)%elnumnod) * maxv ) )
    allocate ( work1( maxval(geometries(:)%nnodes) ) )

    do geom = 1, ngeom

      call fill_nodnod_nc ( geometries(geom)%nnodes, &
        geometries(geom)%elnumnod, geometries(geom)%topology(:,:,1), &
        geometries(geom)%nc, work, work1 )

    end do

    deallocate ( work1, work )

  end subroutine fill_mesh_parts_geometries


! Fill part of the mesh related to connected element edges for curves

  subroutine fill_mesh_parts_geometries_edges ( mesh, curves, blend )

    use limits_m, &
      only: PRINT_CURVE_ELEMENTS_NOT_ON_EDGE, WARN_ON_EDGES_AND_FACES

    type(mesh_t), intent(in) :: mesh
    type(geometry_t), dimension(:), intent(inout) :: curves
    integer, intent(in), optional :: blend

    integer :: curve, nelem, elem, vert, node, nodenr, el, numel
    integer :: elgrpnr, elemnr, numelem, nel, edgenr, edge
    integer :: vnodes(2), edgevnodes(2), ledgevnodes(2), i, j
    integer, allocatable, dimension(:) :: zerocurves
    integer, dimension(:,:), allocatable :: work
    type(element_t), dimension(mesh%nelgrp) :: element
    type(element_t), dimension(size(curves)) :: elementg
    integer, dimension(mesh%nelgrp) :: topoffset
    integer, dimension(size(curves)) :: topoffsetg

!   select blend mesh

    topoffset = 0
    topoffsetg = 0

    if ( present(blend) ) then
      if ( blend == 0 ) then
        element = mesh%element
        elementg = [ ( curves(i)%element, i = 1, size(curves) ) ]
      else
        element = mesh%element_blend(:,blend)
        topoffset = mesh%numnodtop(:,blend+1)
        elementg = [ ( curves(i)%element_blend(blend), i = 1, size(curves) ) ]
        topoffsetg = [ ( curves(i)%numnodtop(blend+1), i = 1, size(curves) ) ]
      end if
    else
      element = mesh%element
      elementg = [ ( curves(i)%element, i = 1, size(curves) ) ]
    end if

!   start find process

    allocate ( work( mesh%nelgrp, maxval(mesh%grpnumel) ) )

    work = 0

    do curve = 1, size(curves)

!     fill nummeshelem and meshelem with edges

      nelem = curves(curve)%nelem

      allocate ( curves(curve)%nummeshelem(nelem) )
      allocate ( curves(curve)%meshelem(nelem) )

!     start loop over elements

      do elem = 1, nelem

!       fill work for elements connected to the nodes on vertices of an edge

        do vert = 1, 2
          node = elementg(curve)%edgevert(vert,1) + topoffsetg(curve)
          nodenr = curves(curve)%topology(node,elem,2)
          vnodes(vert) = nodenr
          numel = mesh%nodnumel(nodenr+1)-mesh%nodnumel(nodenr)
          do el = 1, numel
            elgrpnr = mesh%nodelem( mesh%nodnumel(nodenr) + el, 1 )
            elemnr  = mesh%nodelem( mesh%nodnumel(nodenr) + el, 2 )
            work(elgrpnr,elemnr) = work(elgrpnr,elemnr) + 1
          end do
        end do

!       check whether an element shows up two times (node1 only)

        node = elementg(curve)%edgevert(1,1) + topoffsetg(curve)
        nodenr = curves(curve)%topology(node,elem,2)
        numel = mesh%nodnumel(nodenr+1)-mesh%nodnumel(nodenr)
        numelem = 0
        do el = 1, numel
          elgrpnr = mesh%nodelem( mesh%nodnumel(nodenr) + el, 1 )
          elemnr  = mesh%nodelem( mesh%nodnumel(nodenr) + el, 2 )
          if ( work( elgrpnr, elemnr ) == 2 ) then
!           element found
            numelem = numelem + 1
          end if
        end do

        curves(curve)%nummeshelem(elem) = numelem

!       fill meshelem
        allocate ( curves(curve)%meshelem(elem)%a(3,numelem) )
        nel = 0
        do el = 1, numel
          elgrpnr = mesh%nodelem( mesh%nodnumel(nodenr) + el, 1 )
          elemnr  = mesh%nodelem( mesh%nodnumel(nodenr) + el, 2 )
          if ( work( elgrpnr, elemnr ) == 2 ) then
!           element found
            nel = nel + 1
!           find edgenr
            edgenr = 0
            do edge = 1, element(elgrpnr)%numedges
              ledgevnodes = element(elgrpnr)%edgevert(:,edge) + &
                                                            topoffset(elgrpnr)
              edgevnodes = mesh%topology(elgrpnr)%a(ledgevnodes,elemnr)
              if ( all( vnodes == edgevnodes ) .or. &
                   all( vnodes == edgevnodes(2:1:-1) ) ) then
                edgenr = edge
                exit
              end if
            end do
            curves(curve)%meshelem(elem)%a(:,nel) = &
               [ elgrpnr, elemnr, edgenr ]
          end if
        end do

        if ( any ( curves(curve)%meshelem(elem)%a(3,:) == 0 ) ) then
          write(*,'(/a/a,i0/a/a,i0/)') &
            'Error fill_mesh_parts_geometries_edges:', &
            '  Inconsistent topology of curve ', curve, &
            '  No edge found.', &
            '  elem = ', elem
            stop
        end if

!       set work to zero again

        do vert = 1, 2
          node = elementg(curve)%edgevert(vert,1) + topoffsetg(curve)
          nodenr = curves(curve)%topology(node,elem,2)
          numel = mesh%nodnumel(nodenr+1)-mesh%nodnumel(nodenr)
          do el = 1, numel
            elgrpnr = mesh%nodelem( mesh%nodnumel(nodenr) + el, 1 )
            elemnr  = mesh%nodelem( mesh%nodnumel(nodenr) + el, 2 )
            work(elgrpnr,elemnr) = 0
          end do
        end do

      end do

      if ( any ( curves(curve)%nummeshelem == 0 ) .and. &
                   WARN_ON_EDGES_AND_FACES ) then
        write(*,'(/a/a,i0/a)') &
          'Warning in fill_mesh_parts_geometries_edges: ', &
          '  There are elements on curve ', curve, &
          '  that do not overlap with edges of mesh elements.'
        if ( PRINT_CURVE_ELEMENTS_NOT_ON_EDGE ) then
          write(*,'(/a,i0/a)') &
            '  The following elements on curve ', curve, &
            '  do not overlap with edges of mesh elements:'
          allocate ( zerocurves(count(curves(curve)%nummeshelem == 0)) )
          j = 0
          do i = 1, nelem
            if ( curves(curve)%nummeshelem(i) /= 0 ) cycle
            j = j + 1
            zerocurves(j) = i
          end do
          write (*,*) zerocurves
          write (*,*)
          deallocate ( zerocurves )
        end if
      end if

    end do

    deallocate ( work )

  end subroutine fill_mesh_parts_geometries_edges


! Fill part of the mesh related to connected element faces for surfaces

  subroutine fill_mesh_parts_geometries_faces ( mesh, surfaces, blend )

    use limits_m, &
      only: PRINT_SURFACE_ELEMENTS_NOT_ON_FACE, WARN_ON_EDGES_AND_FACES, &
            WARN_ON_TEST_FOR_VARYING_SIDES
    type(mesh_t), intent(in) :: mesh
    type(geometry_t), dimension(:), intent(inout) :: surfaces
    integer, intent(in), optional :: blend

    integer :: surface, nelem, elem, vert, node, nodenr, el, numel
    integer :: elgrpnr, elemnr, numelem, nel, facenr, face, nv
    integer :: vnodes(4), facevnodes(4), lfacevnodes(4), i, j
    integer, allocatable, dimension(:) :: zerosurfaces
    integer, dimension(:,:), allocatable :: work
    type(element_t), dimension(mesh%nelgrp) :: element
    type(element_t), dimension(size(surfaces)) :: elementg
    integer, dimension(mesh%nelgrp) :: topoffset
    integer, dimension(size(surfaces)) :: topoffsetg


!   test special cases

    if ( any ( mesh%element(:)%globalshape == 'prism' ) .or. &
         any ( mesh%element(:)%globalshape == 'pyramid' ) ) then

!     these meshes can have varying faces and should be excluded for now

      if ( WARN_ON_TEST_FOR_VARYING_SIDES ) then
        write(*,'(/12(a/))') &
          'Warning fill_mesh_parts_geometries_faces: mesh contains ', &
          ' prism and/or pyramid shaped element.', &
          ' Computing meshelem for surfaces is not implemented ', &
          ' for elements having varying side shapes.', &
          ' Set (use limits_m module) ', &
          '   WARN_ON_TEST_FOR_VARYING_SIDES = .false. ', &
          ' to stop this warning from appearing.'
      end if

      do surface = 1, size(surfaces)
        nelem = surfaces(surface)%nelem
        allocate( surfaces(surface)%nummeshelem(nelem) )
        surfaces(surface)%nummeshelem = 0  ! not available
        allocate ( surfaces(surface)%meshelem(nelem) )
        do elem = 1, nelem
          allocate ( surfaces(surface)%meshelem(elem)%a(3,0) )
        end do
      end do

      return

    end if

!   select blend mesh

    topoffset = 0
    topoffsetg = 0

    if ( present(blend) ) then
      if ( blend == 0 ) then
        element = mesh%element
        elementg = [ ( surfaces(i)%element, i = 1, size(surfaces) ) ]
      else
        element = mesh%element_blend(:,blend)
        topoffset = mesh%numnodtop(:,blend+1)
        elementg = &
               [ ( surfaces(i)%element_blend(blend), i = 1, size(surfaces) ) ]
        topoffsetg =  &
                 [ ( surfaces(i)%numnodtop(blend+1), i = 1, size(surfaces) ) ]
      end if
    else
      element = mesh%element
      elementg = [ ( surfaces(i)%element, i = 1, size(surfaces) ) ]
    end if

!   start find process

    allocate ( work( mesh%nelgrp, maxval(mesh%grpnumel) ) )

    work = 0

    do surface = 1, size(surfaces)

!     fill nummeshelem and meshelem with edges

      nelem = surfaces(surface)%nelem

      allocate ( surfaces(surface)%nummeshelem(nelem) )
      allocate ( surfaces(surface)%meshelem(nelem) )

!     start loop over elements

      do elem = 1, nelem

!       fill work for elements connected to the nodes on vertices of a face

        nv = elementg(surface)%facenumvert

        do vert = 1, nv
          node = elementg(surface)%facevert(vert,1) + topoffsetg(surface)
          nodenr = surfaces(surface)%topology(node,elem,2)
          vnodes(vert) = nodenr
          numel = mesh%nodnumel(nodenr+1)-mesh%nodnumel(nodenr)
          do el = 1, numel
            elgrpnr = mesh%nodelem( mesh%nodnumel(nodenr) + el, 1 )
            elemnr  = mesh%nodelem( mesh%nodnumel(nodenr) + el, 2 )
            work(elgrpnr,elemnr) = work(elgrpnr,elemnr) + 1
          end do
        end do

        call sort ( vnodes(:nv) ) ! sort for easy compare

!       check whether an element shows up nv times (node1 only)

        node = elementg(surface)%facevert(1,1) + topoffsetg(surface)
        nodenr = surfaces(surface)%topology(node,elem,2)
        numel = mesh%nodnumel(nodenr+1)-mesh%nodnumel(nodenr)
        numelem = 0
        do el = 1, numel
          elgrpnr = mesh%nodelem( mesh%nodnumel(nodenr) + el, 1 )
          elemnr  = mesh%nodelem( mesh%nodnumel(nodenr) + el, 2 )
          if ( work( elgrpnr, elemnr ) == nv ) then
!           element found
            numelem = numelem + 1
          end if
        end do

        surfaces(surface)%nummeshelem(elem) = numelem

!       fill meshelem
        allocate ( surfaces(surface)%meshelem(elem)%a(3,numelem) )
        nel = 0
        do el = 1, numel
          elgrpnr = mesh%nodelem( mesh%nodnumel(nodenr) + el, 1 )
          elemnr  = mesh%nodelem( mesh%nodnumel(nodenr) + el, 2 )
          if ( work( elgrpnr, elemnr ) == nv ) then
!           element found
            nel = nel + 1
!           find facenr
            facenr = 0
            do face = 1, element(elgrpnr)%numfaces
              lfacevnodes(:nv) = element(elgrpnr)%facevert(:,face) + &
                                                             topoffset(elgrpnr)
              facevnodes(:nv) = &
                             mesh%topology(elgrpnr)%a(lfacevnodes(:nv),elemnr)
              call sort ( facevnodes(:nv) )
              if ( all( vnodes(:nv) == facevnodes(:nv) ) ) then
                facenr = face
                exit
              end if
            end do
            surfaces(surface)%meshelem(elem)%a(:,nel) = &
               [ elgrpnr, elemnr, facenr ]
          end if
        end do

        if ( any ( surfaces(surface)%meshelem(elem)%a(3,:) == 0 ) ) then
          write(*,'(/a/a,i0/a/a,i0/)') &
            'Error fill_mesh_parts_geometries_faces:', &
            '  Inconsistent topology of surface ', surface, &
            '  No face found.', &
            '  elem = ', elem
            stop
        end if

!       set work to zero again

        do vert = 1, nv
          node = elementg(surface)%facevert(vert,1) + topoffsetg(surface)
          nodenr = surfaces(surface)%topology(node,elem,2)
          numel = mesh%nodnumel(nodenr+1)-mesh%nodnumel(nodenr)
          do el = 1, numel
            elgrpnr = mesh%nodelem( mesh%nodnumel(nodenr) + el, 1 )
            elemnr  = mesh%nodelem( mesh%nodnumel(nodenr) + el, 2 )
            work(elgrpnr,elemnr) = 0
          end do
        end do

      end do

      if ( any ( surfaces(surface)%nummeshelem == 0 ) .and. &
                 WARN_ON_EDGES_AND_FACES ) then
        write(*,'(/a/a,i0/a)') &
          'Warning in fill_mesh_parts_geometries_faces: ', &
          '  There are elements on surface ', surface, &
          '  that do not overlap with faces of mesh elements.'
        if ( PRINT_SURFACE_ELEMENTS_NOT_ON_FACE ) then
          write(*,'(/a,i0/a)') &
            '  The following elements on surface ', surface, &
            '  do not overlap with faces of mesh elements:'
          allocate ( zerosurfaces(count(surfaces(surface)%nummeshelem == 0)) )
          j = 0
          do i = 1, nelem
            if ( surfaces(surface)%nummeshelem(i) /= 0 ) cycle
            j = j + 1
            zerosurfaces(j) = i
          end do
          write (*,*) zerosurfaces
          write (*,*)
          deallocate ( zerosurfaces )
        end if
      end if

    end do

    deallocate ( work )

  end subroutine fill_mesh_parts_geometries_faces


! Fill part of the mesh related to connected elements for volumes

  subroutine fill_mesh_parts_geometries_elements ( mesh, volumes )

    type(mesh_t), intent(in) :: mesh
    type(geometry_t), dimension(:), intent(inout) :: volumes

    integer :: volume, nelem, elem

!   volumes: info on connecting elements not available yet.

    do volume = 1, size(volumes)
      nelem =volumes(volume)%nelem
      allocate ( volumes(volume)%nummeshelem(nelem) )
      volumes(volume)%nummeshelem = 0  ! not available yet for volumes
      allocate ( volumes(volume)%meshelem(nelem) )
      do elem = 1, nelem
        allocate ( volumes(volume)%meshelem(elem)%a(3,0) )  ! not available yet
                                                            ! for volumes
      end do
    end do

  end subroutine fill_mesh_parts_geometries_elements


! Fill part of the mesh related to blocks

  subroutine fill_mesh_parts_blocks ( mesh )

    use limits_m, only: NUMBLOCKSX

    type(mesh_t), intent(inout) :: mesh


    integer :: blocks(mesh%ndim)

    if ( mesh%sblocks%filled ) return  ! blocks need not be defined

    if ( mesh%nblocks == 0 ) then

!     set default number of blocks

      blocks = NUMBLOCKSX

      call add_to_mesh ( mesh, blocks=blocks )

      call check_blocks ( mesh )

    else

!     check blocks

      call check_blocks ( mesh )

    end if

  end subroutine fill_mesh_parts_blocks


! Fill part of the mesh related to objects
! Note: filled based on intersection with the main mesh for nblend>0.

  subroutine fill_mesh_parts_objects ( mesh, object1, object2, onlyshapenod, &
    blend )

    type(mesh_t), intent(inout) :: mesh

!   if these are present only objects object1,...,object2 are computed.
!   If only object1 is present one object is computed.
    integer, intent(in), optional :: object1, object2

!   optional parameter to indicate that only numshapenod nodes are
!   used in the geometrical shape of the element for finding the reference
!   coordinates of objects. See heading of routine find_refcoor_objects for more
!   information.
    logical, intent(in), optional :: onlyshapenod

!   optional blend argument
    integer, intent(in), optional :: blend

    integer :: maxv, object, obj1, obj2
    integer, dimension(:), allocatable :: work, work1


    if ( mesh%nobjects == 0 ) return

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

!   compute reference coordinates of point/nodes in the object

    call find_refcoor_objects ( mesh, object1, object2, blend=blend, &
      onlyshapenod=onlyshapenod )

!   compute topology related stuff

    if ( any( mesh%objects(obj1:obj2)%topol ) ) then

!     test on previous call

      if ( any ( mesh%objects(obj1:obj2)%element%numsides /= 0 ) ) then

        write(*,'(/a/a,i0,a,i0,a/)') &
          'Error fill_mesh_parts_objects:', &
          '  Of the objects ', obj1, ' to ', obj2, &
          ' component element is already filled'
          stop

      end if

!     fill element info

      do object = obj1, obj2

        if ( .not. mesh%objects(object)%topol ) cycle

        call fill_element ( mesh%objects(object)%element )

      end do

!     fill element info and nodnumel and nodelem (of nc) for objects

!     first create enough workspace

      allocate ( work( maxval(mesh%objects(obj1:obj2)%nnodes) ) )

      do object = obj1, obj2

        if ( .not. mesh%objects(object)%topol ) cycle

        call fill_nodelem_nc ( mesh%objects(object)%nnodes, &
          mesh%objects(object)%nelem, mesh%objects(object)%elnumnod, &
          mesh%objects(object)%topology, mesh%objects(object)%nc, work )

      end do

      deallocate ( work )

!     fill nodnumnod and nodnod (of nc) for objects

!     first create enough workspace

      maxv = 0
      do object = obj1, obj2
        if ( .not. mesh%objects(object)%topol ) cycle
        maxv = max(size(mesh%objects(object)%nc%nodelem),maxv)
      end do

      allocate ( work( maxval(mesh%objects(obj1:obj2)%elnumnod) * maxv ) )
      allocate ( work1( maxval(mesh%objects(obj1:obj2)%nnodes) ) )

      do object = obj1, obj2

        if ( .not. mesh%objects(object)%topol ) cycle

        call fill_nodnod_nc ( mesh%objects(object)%nnodes, &
          mesh%objects(object)%elnumnod, mesh%objects(object)%topology, &
          mesh%objects(object)%nc, work, work1 )

      end do

      deallocate ( work1, work )

    end if

    mesh%objects(obj1:obj2)%meshparts = .true.

  end subroutine fill_mesh_parts_objects


! subroutine for filling nodnumel and nodelem of nc

  subroutine fill_nodelem_nc ( nnodes, nelem, elnumnod, topology, nc, work )

    integer, intent(in) :: nnodes, nelem, elnumnod
    integer, dimension(:,:), intent(in) :: topology
    type(connectdata_t), intent(inout) :: nc
    integer, dimension(:), intent(inout) :: work


    integer :: elem, node, nodenr


!   fill nodnumel

    allocate ( nc%nodnumel(nnodes+1) )

    nc%nodnumel = 0

    do elem = 1, nelem
      do node = 1, elnumnod
        nodenr = topology(node,elem)
        nc%nodnumel(nodenr+1) = nc%nodnumel(nodenr+1) + 1
      end do
    end do

    nc%maxnodnumel = maxval(nc%nodnumel)

!   accumulate nodnumel

    do node = 2, nnodes
      nc%nodnumel(node+1) = nc%nodnumel(node) + nc%nodnumel(node+1)
    end do

!   fill nodelem

    allocate ( nc%nodelem( nc%nodnumel ( nnodes + 1 ) ) )

    work = 0

    do elem = 1, nelem
      do node = 1, elnumnod
        nodenr = topology(node,elem)
        work(nodenr) = work(nodenr) + 1
        nc%nodelem( nc%nodnumel(nodenr) + work(nodenr) ) = elem
      end do
    end do

  end subroutine fill_nodelem_nc


! subroutine for filling nodnumnod and nodnod of nc

  subroutine fill_nodnod_nc ( nnodes, elnumnod, topology, nc, work, work1 )

    integer, intent(in) :: nnodes, elnumnod
    integer, dimension(:,:), intent(in) :: topology
    type(connectdata_t), intent(inout) :: nc
    integer, dimension(:), intent(inout) :: work, work1


    integer :: elem, node, nodenr, elm
    integer :: nn, ipwork, nod

!   fill nodnumnod

    allocate ( nc%nodnumnod(nnodes+1) )

    nc%nodnumnod(1) = 0

    ipwork = 0   ! pointer in array work
    work1  = 0

    do nodenr = 1, nnodes

      work1(nodenr) = 1  ! Exclude node itself
      nn = 0

!     fill work with nodal points from elements connected

      do elm = 1, nc%nodnumel(nodenr+1) - nc%nodnumel(nodenr)

        elem  = nc%nodelem( nc%nodnumel(nodenr) + elm )

        do node = 1, elnumnod

          nod = topology(node,elem)

          if ( work1(nod) /= 0 ) cycle  ! node already used

          work1(nod) = 1 ! new node
          nn = nn + 1
          work(ipwork+nn) = nod

        end do

      end do

!     put work1 back to zero

      work1(nodenr) = 0
      work1(work(ipwork+1:ipwork+nn)) = 0

!     fill nodnumnod

      nc%nodnumnod(nodenr+1) = nn
      ipwork = ipwork + nn

    end do

    nc%maxnodnumnod = maxval(nc%nodnumnod)

!   accumulate nodnumnod

    do node = 2, nnodes
      nc%nodnumnod(node+1) = nc%nodnumnod(node) + nc%nodnumnod(node+1)
    end do

!   fill nonnod

    allocate ( nc%nodnod( nc%nodnumnod(nnodes+1) ) )

    nc%nodnod = work(1:nc%nodnumnod(nnodes+1))

  end subroutine fill_nodnod_nc


! Fill elements

  subroutine fill_element ( element )

    type(element_t), intent(inout) :: element

!   Fill element type info for single element group.
!   NOTE: elshape, globalshape, numnod and ndim must have been filled already.
!         For high-order elements also p must have been filled already.

!   element info for the core of TFEM

    call fill_element_core ( element )

!   element info for the figplot addon

    call fill_element_figplot ( element )

!   element info for VTK output

    call fill_element_vtk ( element )

!   element info for Tecplot output

    call fill_element_tecplot ( element )

!   element info for Gmsh output

    call fill_element_gmsh ( element )

!   element info for Gmsh parsed output

    call fill_element_gmsh_parsed ( element )

  end subroutine fill_element


! Fill elements for the core of TFEM

  subroutine fill_element_core ( element )

    type(element_t), intent(inout) :: element

!   Fill element type info for single element group.
!   NOTE: elshape, globalshape, numnod and ndim must have been filled already.
!         For high-order elements also p must have been filled already.

    integer :: i, j, k, p, q, r

    select case ( element%elshape )

      case(1) ! two-node line element

        element%ndimr = 1           ! dimension of reference space
        element%numshapenod = 2     ! number of `shape' nodal points for mapping
        element%numsides = 2        ! number of sides
        element%sidnumvert = 1      ! number of vertices of a side
        element%sidnumnod = 1       ! number of nodes of a side
        element%numinternnod = 0    ! number internal nodes
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%numedges = 1        ! number of edges
        element%edgenumnod = 2      ! number of nodes of an edge
        element%numfaces = 0        ! number of faces
        element%facenumvert = 0     ! number of vertices of a face
        element%facenumnod = 0      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = reshape ( [ 1, 2 ], [1,2] )

        element%sidnod  = reshape ( [ 1, 2 ], [1,2] )

        !element%internnod =  []

        !element%sidinternnod =  []

        element%edgevert = reshape ( [ 1, 2 ], [2,1] )

        element%edgenod  = reshape ( [ 1, 2 ], [2,1] )

        !element%facevert = []

        !element%facenod = []

        element%xc = 0

      case(2,32) ! three-node line element

        element%ndimr = 1           ! dimension of reference space
        element%numshapenod = 3     ! number of `shape' nodal points for mapping
        element%numsides = 2        ! number of sides
        element%sidnumvert = 1      ! number of vertices of a side
        element%sidnumnod = 1       ! number of nodes of a side
        element%numinternnod = 1    ! number internal nodes
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%numedges = 1        ! number of edges
        element%edgenumnod = 3      ! number of nodes of an edge
        element%numfaces = 0        ! number of faces
        element%facenumvert = 0     ! number of vertices of a face
        element%facenumnod = 0      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = reshape ( [ 1, 3 ], [1,2] )

        element%sidnod  = reshape ( [ 1, 3 ], [1,2] )

        element%internnod = [ 2 ]

        !element%sidinternnod =  []

        element%edgevert = reshape ( [ 1, 3 ], [2,1] )

        element%edgenod  = reshape ( [ 1, 2, 3 ], [3,1] )

        !element%facevert = []

        !element%facenod = []

        element%xc = 0

      case(3,10) ! three-node and four-node triangle

        element%ndimr = 2           ! dimension of reference space
        element%numshapenod = 3     ! number of `shape' nodal points for mapping
        element%numsides = 3        ! number of sides
        element%sidnumvert = 2      ! number of vertices of a side
        element%sidnumnod = 2       ! number of nodes of a side
        select case ( element%elshape )
        case(3)
          element%numinternnod = 0  ! number internal nodes
        case(10)
          element%numinternnod = 1  ! number internal nodes
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%numedges = 3        ! number of edges
        element%edgenumnod = 2      ! number of nodes of an edge
        element%numfaces = 1        ! number of faces
        element%facenumvert = 3     ! number of vertices of a face
        element%facenumnod = 3      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = reshape ( [ 1,2, 2,3, 3,1 ], [2,3] )

        element%sidnod  = reshape ( [ 1,2, 2,3, 3,1 ], [2,3] )

        select case ( element%elshape )
        case(3)
          !element%internnod = []
        case(10)
          element%internnod = [4]
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select

        !element%sidinternnod =  []

        element%edgevert = element%sidvert

        element%edgenod  = element%sidnod

        element%facevert = reshape ( [ 1, 2, 3 ], [3,1] )

        element%facenod  = reshape ( [ 1, 2, 3 ], [3,1] )

        element%xc = 1._dp/3

      case(4,7,33) ! six-node and seven-node triangle

        element%ndimr = 2           ! dimension of reference space
        element%numshapenod = 6     ! number of `shape' nodal points for mapping
        element%numsides = 3        ! number of sides
        element%sidnumvert = 2      ! number of vertices of a side
        element%sidnumnod = 3       ! number of nodes of a side
        select case ( element%elshape )
        case(4,33)
          element%numinternnod = 0  ! number internal nodes
        case(7)
          element%numinternnod = 1  ! number internal nodes
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select
        element%sidnuminternnod = 1 ! number internal nodes on a side
        element%numedges = 3        ! number of edges
        element%edgenumnod = 3      ! number of nodes of an edge
        element%numfaces = 1        ! number of faces
        element%facenumvert = 3     ! number of vertices of a face
        select case ( element%elshape )
        case(4,33)
          element%facenumnod = 6    ! number of nodes of a face
        case(7)
          element%facenumnod = 7    ! number of nodes of a face
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select

        call allocate_element_arrays

        element%sidvert = reshape ( [ 1,3, 3,5, 5,1 ], [2,3] )

        element%sidnod = reshape ( [ 1,2,3, 3,4,5, 5,6,1 ], [3,3] )

        select case ( element%elshape )
        case(4,33)
          !element%internnod =  []
        case(7)
          element%internnod =  [7]
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select

        element%sidinternnod = reshape ( [ 2, 4, 6 ], [1,3] )

        element%edgevert = element%sidvert

        element%edgenod = element%sidnod

        element%facevert = reshape ( [ 1, 3, 5 ], [3,1] )

        element%facenod = reshape ( [ (i,i=1,element%facenumnod) ], &
                                        [element%facenumnod,1] )

        element%xc = 1._dp/3

      case(5,9) ! four-node and five-node quadrilateral

        element%ndimr = 2           ! dimension of reference space
        element%numshapenod = 4     ! number of `shape' nodal points for mapping
        element%numsides = 4        ! number of sides
        element%sidnumvert = 2      ! number of vertices of a side
        element%sidnumnod = 2       ! number of nodes of a side
        select case ( element%elshape )
        case(5)
          element%numinternnod = 0  ! number internal nodes
        case(9)
          element%numinternnod = 1  ! number internal nodes
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%numedges = 4        ! number of edges
        element%edgenumnod = 2      ! number of nodes of an edge
        element%numfaces = 1        ! number of faces
        element%facenumvert = 4     ! number of vertices of a face
        select case ( element%elshape )
        case(5)
          element%facenumnod = 4      ! number of nodes of a face
        case(9)
          element%facenumnod = 5      ! number of nodes of a face
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select

        call allocate_element_arrays

        element%sidvert = reshape ( [ 1,2, 2,3, 3,4, 4,1 ], [2,4] )

        element%sidnod  = reshape ( [ 1,2, 2,3, 3,4, 4,1 ], [2,4] )

        select case ( element%elshape )
        case(5)
          !element%internnod = []
        case(9)
          element%internnod = [5]
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select

        !element%sidinternnod =  []

        element%edgevert = element%sidvert

        element%edgenod = element%sidnod

        element%facevert = reshape ( [ 1, 2, 3, 4 ], [4,1] )

        element%facenod = reshape ( [ (i,i=1,element%facenumnod) ], &
                                        [element%facenumnod,1] )

        element%xc = 0

      case(6,34) ! nine-node quadrilateral

        element%ndimr = 2           ! dimension of reference space
        element%numshapenod = 9     ! number of `shape' nodal points for mapping
        element%numsides = 4        ! number of sides
        element%sidnumvert = 2      ! number of vertices of a side
        element%sidnumnod = 3       ! number of nodes of a side
        element%numinternnod = 1    ! number internal nodes
        element%sidnuminternnod = 1 ! number internal nodes on a side
        element%numedges = 4        ! number of edges
        element%edgenumnod = 3      ! number of nodes of an edge
        element%numfaces = 1        ! number of faces
        element%facenumvert = 4     ! number of vertices of a face
        element%facenumnod = 9      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = reshape ( [ 1,3, 3,5, 5,7, 7,1 ], [2,4] )

        element%sidnod = reshape ( [ 1,2,3, 3,4,5, 5,6,7, 7,8,1 ], [3,4] )

        element%internnod =  [ 9 ]

        element%sidinternnod = reshape ( [ 2, 4, 6, 8 ], [1,4] )

        element%edgevert = element%sidvert

        element%edgenod = element%sidnod

        element%facevert = reshape ( [ 1, 3, 5, 7 ], [4,1] )

        element%facenod = reshape ( [ (i,i=1,element%facenumnod) ], &
                                        [element%facenumnod,1] )

        element%xc = 0

      case(11,18) ! four-node and five-node tetrahedron

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = 4     ! number of `shape' nodal points for mapping
        element%numsides = 4        ! number of sides
        element%sidnumvert = 3      ! number of vertices of a side
        element%sidnumnod = 3       ! number of nodes of a side
        select case ( element%elshape )
        case(11)
          element%numinternnod = 0   ! number internal nodes
        case(18)
          element%numinternnod = 1   ! number internal nodes
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%numedges = 6        ! number of edges
        element%edgenumnod = 2      ! number of nodes of an edge
        element%numfaces = 4        ! number of faces
        element%facenumvert = 3     ! number of vertices of a face
        element%facenumnod = 3      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = reshape ( [ 1,3,2, 1,2,4, 2,3,4, 3,1,4 ], [3,4] )

        element%sidnod = reshape ( [ 1,3,2, 1,2,4, 2,3,4, 3,1,4 ], [3,4] )

        select case ( element%elshape )
        case(11)
          !element%internnod = []
        case(18)
          element%internnod = [5]
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select

        ! element%sidinternnod =  []

        element%edgevert = &
                reshape ( [ 1,2, 2,3, 3,1, 1,4, 2,4, 3,4 ], [2,6] )

        element%edgenod = element%edgevert

        element%facevert = element%sidvert

        element%facenod = element%sidnod

        element%xc = 1._dp/4

      case(12,35) ! ten-node tetrahedron

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = 10    ! number of `shape' nodal points for mapping
        element%numsides = 4        ! number of sides
        element%sidnumvert = 3      ! number of vertices of a side
        element%sidnumnod = 6       ! number of nodes of a side
        element%numinternnod = 0    ! number internal nodes
        element%numedges = 6        ! number of edges
        element%edgenumnod = 3      ! number of nodes of an edge
        element%numfaces = 4        ! number of faces
        element%facenumvert = 3     ! number of vertices of a face
        element%facenumnod = 6      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = &
          reshape ( [ 1,5,3, 1,3,10, 3,5,10, 5,1,10 ], [3,4] )

        element%sidnod = reshape ( [ 1,6,5,4,3,2, 1,2,3,8,10,7, &
          3,4,5,9,10,8, 5,6,1,7,10,9 ], [6,4] )

        !element%internnod =  []

        element%sidnuminternnod = 0   ! number internal nodes on a side

        !element%sidinternnod =  []

        element%edgevert = &
                reshape ( [ 1,3, 3,5, 5,1, 1,10, 3,10, 5,10 ], [2,6] )

        element%edgenod = &
          reshape ( [ 1,2,3, 3,4,5, 5,6,1, 1,7,10, 3,8,10, 5,9,10 ], [3,6] )

        element%facevert = element%sidvert

        element%facenod = element%sidnod

        element%xc = 1._dp/4

      case(13,17) ! eight-node and nine-node hexahedron

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = 8     ! number of `shape' nodal points for mapping
        element%numsides = 6        ! number of sides
        element%sidnumvert = 4      ! number of vertices of a side
        element%sidnumnod = 4       ! number of nodes of a side
        select case ( element%elshape )
        case(13)
          element%numinternnod = 0  ! number internal nodes
        case(17)
          element%numinternnod = 1  ! number internal nodes
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%numedges = 12       ! number of edges
        element%edgenumnod = 2      ! number of nodes of an edge
        element%numfaces = 6        ! number of faces
        element%facenumvert = 4     ! number of vertices of a face
        element%facenumnod = 4      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = &
          reshape ( [ 1,4,3,2, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8 &
                     ], [4,6] )

        element%sidnod = &
          reshape ( [ 1,4,3,2, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8 &
                     ], [4,6] )

        select case ( element%elshape )
        case(13)
          !element%internnod =  []
        case(17)
          element%internnod =  [9]
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select

        !element%sidinternnod =  []

        element%edgevert = &
                reshape ( [ 1,2, 2,3, 3,4, 4,1, &
                            1,5, 2,6, 3,7, 4,8, &
                            5,6, 6,7, 7,8, 8,5 ], [2,12] )

        element%edgenod = element%edgevert

        element%facevert = element%sidvert

        element%facenod = element%sidnod

        element%xc = 0

      case(14,36) ! 27-node hexahedron

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = 27    ! number of `shape' nodal points for mapping
        element%numsides = 6        ! number of sides
        element%sidnumvert = 4      ! number of vertices of a side
        element%sidnumnod = 9       ! number of nodes of a side
        element%numinternnod = 1    ! number internal nodes
        element%sidnuminternnod = 1 ! number internal nodes on a side
        element%numedges = 12       ! number of edges
        element%edgenumnod = 3      ! number of nodes of an edge
        element%numfaces = 6        ! number of faces
        element%facenumvert = 4     ! number of vertices of a face
        element%facenumnod = 9      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = &
          reshape ( [ 1,7,9,3, 1,3,21,19, 3,9,27,21, 9,7,25,27, 7,1,19,25, &
                      19,21,27,25 ], [4,6] )

        element%sidnod = &
          reshape ( [ 1,4,7,8,9,6,3,2,5, 1,2,3,12,21,20,19,10,11, &
                      3,6,9,18,27,24,21,12,15, 9,8,7,16,25,26,27,18,17, &
                      7,4,1,10,19,22,25,16,13, 19,20,21,24,27,26,25,22,23 &
                     ], [9,6] )

        element%internnod =  [14]

        element%sidinternnod = reshape ( [ 5, 11, 15, 17, 13, 23 ], [1,6] )

        element%edgevert = &
                reshape ( [ 1,3,   3,9,   9,7,   7,1,  &
                            1,19,  3,21,  9,27,  7,25, &
                            19,21, 21,27, 27,25, 25,19 ], [2,12] )

        element%edgenod = &
                reshape ( [ 1,2,3,    3,6,9,    9,8,7,    7,4,1,  &
                            1,10,19,  3,12,21,  9,18,27,  7,16,25, &
                            19,20,21, 21,24,27, 27,26,25, 25,22,19 ], [3,12] )

        element%facevert = element%sidvert

        element%facenod = element%sidnod

        element%xc = 0

      case(15,16) ! 14-node and 15-node tetrahedron

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = 10    ! number of `shape' nodal points for mapping
        element%numsides = 4        ! number of sides
        element%sidnumvert = 3      ! number of vertices of a side
        element%sidnumnod = 7       ! number of nodes of a side
        select case ( element%elshape )
        case(15)
          element%numinternnod = 0      ! number internal nodes
        case(16)
          element%numinternnod = 1      ! number internal nodes
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select
        element%sidnuminternnod = 1 ! number internal nodes on a side
        element%numedges = 6        ! number of edges
        element%edgenumnod = 3      ! number of nodes of an edge
        element%numfaces = 4        ! number of faces
        element%facenumvert = 3     ! number of vertices of a face
        element%facenumnod = 7      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = &
          reshape ( [ 1,5,3, 1,3,10, 3,5,10, 5,1,10 ], [3,4] )

        element%sidnod = reshape ( [ 1,6,5,4,3,2,11, 1,2,3,8,10,7,12, &
          3,4,5,9,10,8,13, 5,6,1,7,10,9,14 ], [7,4] )

        select case ( element%elshape )
        case(15)
          !element%internnod =  []
        case(16)
          element%internnod =  [15]
        case default
          call errormsg_case_default ( 'fill_element_core', &
            'element%elshape', int_value=element%elshape )
        end select

        element%sidinternnod = reshape ( [ 11, 12, 13, 14 ], [1,4] )

        element%edgevert = &
                reshape ( [ 1,3, 3,5, 5,1, 1,10, 3,10, 5,10 ], [2,6] )

        element%edgenod = &
          reshape ( [ 1,2,3, 3,4,5, 5,6,1, 1,7,10, 3,8,10, 5,9,10 ], [3,6] )

        element%facevert = element%sidvert

        element%facenod = element%sidnod

        element%xc = 1._dp/4

      case(30) ! eight-node quadrilateral

        element%ndimr = 2           ! dimension of reference space
        element%numshapenod = 8     ! number of `shape' nodal points for mapping
        element%numsides = 4        ! number of sides
        element%sidnumvert = 2      ! number of vertices of a side
        element%sidnumnod = 3       ! number of nodes of a side
        element%numinternnod = 0    ! number internal nodes
        element%sidnuminternnod = 1 ! number internal nodes on a side
        element%numedges = 4        ! number of edges
        element%edgenumnod = 3      ! number of nodes of an edge
        element%numfaces = 1        ! number of faces
        element%facenumvert = 4     ! number of vertices of a face
        element%facenumnod = 8      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = reshape ( [ 1,3, 3,5, 5,7, 7,1 ], [2,4] )

        element%sidnod = reshape ( [ 1,2,3, 3,4,5, 5,6,7, 7,8,1 ], [3,4] )

        element%sidinternnod = reshape ( [ 2, 4, 6, 8 ], [1,4] )

        element%edgevert = element%sidvert

        element%edgenod = element%sidnod

        element%facevert = reshape ( [ 1, 3, 5, 7 ], [4,1] )

        element%facenod = reshape ( [ (i,i=1,element%facenumnod) ], &
                                        [element%facenumnod,1] )

        element%xc = 0

      case(31) ! 20-node hexahedron

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = 20    ! number of `shape' nodal points for mapping
        element%numsides = 6        ! number of sides
        element%sidnumvert = 4      ! number of vertices of a side
        element%sidnumnod = 8       ! number of nodes of a side
        element%numinternnod = 0    ! number internal nodes
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%numedges = 12       ! number of edges
        element%edgenumnod = 3      ! number of nodes of an edge
        element%numfaces = 6        ! number of faces
        element%facenumvert = 4     ! number of vertices of a face
        element%facenumnod = 8      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = &
          reshape ( [ 1,6,8,3, 1,3,15,13, 3,8,20,15, 8,6,18,20, 6,1,13,18, &
                       13,15,20,18 ], [4,6] )

        element%sidnod = &
          reshape ( [ 1,4,6,7,8,5,3,2, 1,2,3,10,15,14,13,9, &
                      3,5,8,12,20,17,15,10, 8,7,6,11,18,19,20,12, &
                      6,4,1,9,13,16,18,11, 13,14,15,17,20,19,18,16 ], [8,6] )

        !element%sidinternnod =  []

        element%edgevert = &
                reshape ( [ 1,3,   3,8,   8,6,   6,1,  &
                            1,13,  3,15,  8,20,  6,18, &
                            13,15, 15,20, 20,18, 18,13 ], [2,12] )

        element%edgenod = &
                reshape ( [ 1,2,3,    3,5,8,    8,7,6,    6,4,1,   &
                            1,9,13,   3,10,15,  8,12,20,  6,11,18, &
                            13,14,15, 15,17,20, 20,19,18, 18,16,13 ], [3,12] )

        element%facevert = element%sidvert

        element%facenod = element%sidnod

        element%xc = 0

      case(41) ! six-node prism

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = 6     ! number of `shape' nodal points for mapping
        element%numsides = 5        ! number of sides
        element%numinternnod = 0    ! number internal nodes
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%numedges = 9        ! number of edges
        element%edgenumnod = 2      ! number of nodes of an edge
        element%numfaces = 5        ! number of faces

!       The following not set. DON'T USE! Extension of element_t required.
        element%sidnumvert = 0      ! number of vertices of a side
        element%sidnumnod = 0       ! number of nodes of a side
        element%facenumvert = 0     ! number of vertices of a face
        element%facenumnod = 0      ! number of nodes of a face

        call allocate_element_arrays

        ! element%internnod = []

        ! element%sidinternnod =  []

        element%edgevert = &
             reshape ( [ 1,2, 2,3, 3,1, 1,4, 2,5, 3,6, 4,5, 5,6, 6,4 ], [2,9] )

        element%edgenod = element%edgevert

        element%facevert = element%sidvert

        element%facenod = element%sidnod

        element%xc = [ 1._dp/3, 1._dp/3, 0._dp ]

!       The following not set. DON'T USE! Extension of element_t required.

        ! element%sidvert
        ! element%sidnod
        ! element%facevert
        ! element%facenod

      case(42) ! fifteen-node prism

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = 15    ! number of `shape' nodal points for mapping
        element%numsides = 5        ! number of sides
        element%numinternnod = 0    ! number internal nodes
        element%numedges = 9        ! number of edges
        element%edgenumnod = 3      ! number of nodes of an edge
        element%numfaces = 5        ! number of faces

!       The following not set. DON'T USE! Extension of element_t required.
        element%sidnumvert = 0      ! number of vertices of a side
        element%sidnumnod = 0       ! number of nodes of a side
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%facenumvert = 0     ! number of vertices of a face
        element%facenumnod = 0      ! number of nodes of a face

        call allocate_element_arrays

        ! element%internnod = []

        element%edgevert = reshape ( &
             [ 1,3, 3,5, 5,1, 1,10, 3,12, 5,14, 10,12, 12,14, 14,10 ], [2,9] )

        element%edgenod = &
                reshape ( [ 1,2,3,   3,4,5,    5,6,1,    1,7,10,  &
                            3,8,12,  5,9,14,  10,11,12, 12,13,14, &
                            14,15,10 ], [3,9] )

        element%facevert = element%sidvert

        element%facenod = element%sidnod

        element%xc = [ 1._dp/3, 1._dp/3, 0._dp ]

!       The following not set. DON'T USE! Extension of element_t required.

        ! element%sidvert
        ! element%sidnod
        ! element%sidinternnod
        ! element%facevert
        ! element%facenod

      case(43) ! eighteen-node prism

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = 18    ! number of `shape' nodal points for mapping
        element%numsides = 5        ! number of sides
        element%numinternnod = 0    ! number internal nodes
        element%numedges = 9        ! number of edges
        element%edgenumnod = 3      ! number of nodes of an edge
        element%numfaces = 5        ! number of faces

!       The following not set. DON'T USE! Extension of element_t required.
        element%sidnumvert = 0      ! number of vertices of a side
        element%sidnumnod = 0       ! number of nodes of a side
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%facenumvert = 0     ! number of vertices of a face
        element%facenumnod = 0      ! number of nodes of a face

        call allocate_element_arrays

        ! element%internnod = []

        element%edgevert = reshape ( &
             [ 1,3, 3,5, 5,1, 1,13, 3,15, 5,17, 13,15, 15,17, 17,13 ], [2,9] )

        element%edgenod = &
                reshape ( [ 1,2,3,   3,4,5,    5,6,1,    1,7,13,   &
                            3,9,15,  5,11,17,  13,14,15, 15,16,17, &
                            17,18,13 ], [3,9] )

        element%facevert = element%sidvert

        element%facenod = element%sidnod

        element%xc = [ 1._dp/3, 1._dp/3, 0._dp ]

!       The following not set. DON'T USE! Extension of element_t required.

        ! element%sidvert
        ! element%sidnod
        ! element%sidinternnod
        ! element%facevert
        ! element%facenod

      case(51) ! five-node pyramid

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = 5     ! number of `shape' nodal points for mapping
        element%numsides = 5        ! number of sides
        element%numinternnod = 0    ! number internal nodes
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%numedges = 8        ! number of edges
        element%edgenumnod = 2      ! number of nodes of an edge
        element%numfaces = 5        ! number of faces

!       The following not set. DON'T USE! Extension of element_t required.
        element%sidnumvert = 0      ! number of vertices of a side
        element%sidnumnod = 0       ! number of nodes of a side
        element%facenumvert = 0     ! number of vertices of a face
        element%facenumnod = 0      ! number of nodes of a face

        call allocate_element_arrays

        ! element%internnod = []

        ! element%sidinternnod =  []

        element%edgevert = &
             reshape ( [ 1,2, 2,3, 3,4, 4,1, 1,5, 2,5, 3,5, 4,5 ], [2,8] )

        element%edgenod = element%edgevert

        element%facevert = element%sidvert

        element%facenod = element%sidnod

        element%xc = [ 0._dp, 0._dp, 1._dp/4 ]

!       The following not set. DON'T USE! Extension of element_t required.

        ! element%sidvert
        ! element%sidnod
        ! element%facevert
        ! element%facenod

      case(52) ! thirteen-node pyramid

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = 13    ! number of `shape' nodal points for mapping
        element%numsides = 5        ! number of sides
        element%numinternnod = 0    ! number internal nodes
        element%numedges = 8        ! number of edges
        element%edgenumnod = 3      ! number of nodes of an edge
        element%numfaces = 5        ! number of faces

!       The following not set. DON'T USE! Extension of element_t required.
        element%sidnumvert = 0      ! number of vertices of a side
        element%sidnumnod = 0       ! number of nodes of a side
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%facenumvert = 0     ! number of vertices of a face
        element%facenumnod = 0      ! number of nodes of a face

        call allocate_element_arrays

        ! element%internnod = []

        element%edgevert = reshape ( &
             [ 1,3, 3,5, 5,7, 7,1, 1,13, 3,13, 5,13, 7,13 ], [2,8] )

        element%edgenod = &
                reshape ( [ 1,2,3,   3,4,5,    5,6,7,    7,8,1,   &
                            1,9,13,  3,10,13,  5,11,13,  7,12,13 ], &
                            [3,8] )

        element%facevert = element%sidvert

        element%facenod = element%sidnod

        element%xc = [ 0._dp, 0._dp, 1._dp/4 ]

!       The following not set. DON'T USE! Extension of element_t required.

        ! element%sidvert
        ! element%sidnod
        ! element%sidinternnod
        ! element%facevert
        ! element%facenod

      case(53) ! fourteen-node pyramid

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = 14    ! number of `shape' nodal points for mapping
        element%numsides = 5        ! number of sides
        element%numinternnod = 0    ! number internal nodes
        element%numedges = 8        ! number of edges
        element%edgenumnod = 3      ! number of nodes of an edge
        element%numfaces = 5        ! number of faces

!       The following not set. DON'T USE! Extension of element_t required.
        element%sidnumvert = 0      ! number of vertices of a side
        element%sidnumnod = 0       ! number of nodes of a side
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%facenumvert = 0     ! number of vertices of a face
        element%facenumnod = 0      ! number of nodes of a face

        call allocate_element_arrays

        ! element%internnod = []

        element%edgevert = reshape ( &
             [ 1,3, 3,5, 5,7, 7,1, 1,14, 3,14, 5,14, 7,14 ], [2,8] )

        element%edgenod = &
                reshape ( [ 1,2,3,    3,4,5,    5,6,7,    7,8,1,   &
                            1,10,14,  3,11,14,  5,12,14,  7,13,14 ], &
                            [3,8] )

        element%facevert = element%sidvert

        element%facenod = element%sidnod

        element%xc = [ 0._dp, 0._dp, 1._dp/4 ]

!       The following not set. DON'T USE! Extension of element_t required.

        ! element%sidvert
        ! element%sidnod
        ! element%sidinternnod
        ! element%facevert
        ! element%facenod

      case(101,105) ! line high-order element (including macro variant)

        p = element%p(1,1)

        element%ndimr = 1           ! dimension of reference space
        element%numshapenod = p+1   ! number of `shape' nodal points for mapping
        element%numsides = 2        ! number of sides
        element%sidnumvert = 1      ! number of vertices of a side
        element%sidnumnod = 1       ! number of nodes of a side
        element%numinternnod = p-1  ! number internal nodes
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%numedges = 1        ! number of edges
        element%edgenumnod = p+1    ! number of nodes of an edge
        element%numfaces = 0        ! number of faces
        element%facenumvert = 0     ! number of vertices of a face
        element%facenumnod = 0      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = reshape ( [ 1, p+1 ], [1,2] )

        element%sidnod  = reshape ( [ 1, p+1 ], [1,2] )

        element%internnod = [ (i,i=2,p) ]

        !element%sidinternnod =  []

        element%edgevert = reshape ( [ 1, p+1 ], [2,1] )

        element%edgenod  = reshape ( [ (i,i=1,p+1) ], [p+1,1] )

        !element%facevert = []

        !element%facenod = []

        element%xc = 0

      case(102,106) ! quadrilateral high-order element (including macro variant)

        p = element%p(1,1)
        q = element%p(2,1)

        element%ndimr = 2           ! dimension of reference space
        element%numshapenod = (p+1)*(q+1)   ! number of `shape' nodal points
                                            ! for mapping
        element%numsides = 4                ! number of sides
        element%sidnumvert = 2              ! number of vertices of a side
        element%numinternnod = (p-1)*(q-1)  ! number internal nodes
        element%numedges = 4                ! number of edges
        element%numfaces = 1                ! number of faces
        element%facenumvert = 4             ! number of vertices of a face
        element%facenumnod = (p+1)*(q+1)    ! number of nodes of a face

!       the following are NOT available: DON'T USE THEM
        element%sidnumnod = 0       ! number of nodes of a side
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%edgenumnod = 0      ! number of nodes of an edge

        call allocate_element_arrays

        element%sidvert = reshape ( [ 1, p+1, p+1, (p+1)*(q+1), (p+1)*(q+1), &
                                       (p+1)*q+1, (p+1)*q+1, 1 ], [2,4] )

        element%internnod =  [ ((i+(p+1)*(j-1),i=2,p),j=2,q) ]

        element%edgevert = element%sidvert

        element%facevert = reshape ( [ 1, p+1, (p+1)*(q+1), (p+1)*q+1 ], &
                                       [4,1] )

        element%facenod = reshape ( [ (i,i=1,element%facenumnod) ], &
                                        [element%facenumnod,1] )

!       the following are NOT available: DON'T USE THEM
!          element%sidnod = []
!          element%sidinternnod = []
!          element%edgenod = []

        element%xc = 0

      case(103,107) ! hexahedron high-order element (including macro variant)

        p = element%p(1,1)
        q = element%p(2,1)
        r = element%p(3,1)

        element%ndimr = 3           ! dimension of reference space
        element%numshapenod = (p+1)*(q+1)*(r+1)   ! number of `shape' nodal
                                                  ! points for mapping
        element%numsides = 6                      ! number of sides
        element%sidnumvert = 4                    ! number of vertices of a side
        element%numinternnod = (p-1)*(q-1)*(r-1)  ! number internal nodes
        element%numedges = 12                     ! number of edges
        element%numfaces = 6                      ! number of faces
        element%facenumvert = 4                   ! number of vertices of a face

!       the following are NOT available: DON'T USE THEM
        element%sidnumnod = 0       ! number of nodes of a side
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%edgenumnod = 0      ! number of nodes of an edge
        element%facenumnod = 0      ! number of nodes of a face

        call allocate_element_arrays

        element%sidvert = &
          reshape ( [ 1, (p+1)*q+1, (p+1)*(q+1), p+1, &
                      1, p+1, (p+1)*((q+1)*r+1), (p+1)*(q+1)*r+1, &
                      p+1, (p+1)*(q+1), (p+1)*(q+1)*(r+1), (p+1)*((q+1)*r+1), &
                      (p+1)*(q+1), (p+1)*q+1, (p+1)*(q+1)*(r+1)-p, &
                      (p+1)*(q+1)*(r+1), &
                      (p+1)*q+1, 1, (p+1)*(q+1)*r+1, (p+1)*(q+1)*(r+1)-p, &
                      (p+1)*(q+1)*r+1, (p+1)*((q+1)*r+1), (p+1)*(q+1)*(r+1), &
                      (p+1)*(q+1)*(r+1)-p ], [4,6] )

        element%internnod = &
                    [ (((i+(p+1)*(j-1)+(p+1)*(q+1)*(k-1),i=2,p),j=2,q),k=2,r) ]

        element%edgevert = &
                reshape ( [ 1, p+1, p+1, (p+1)*(q+1), (p+1)*(q+1), (p+1)*q+1, &
                            (p+1)*q+1, 1, 1, (p+1)*(q+1)*r+1, p+1, &
                            (p+1)*((q+1)*r+1), (p+1)*(q+1), &
                            (p+1)*(q+1)*(r+1), (p+1)*q+1,(p+1)*(q+1)*(r+1)-p,&
                            (p+1)*(q+1)*r+1, (p+1)*((q+1)*r+1), &
                            (p+1)*((q+1)*r+1), (p+1)*(q+1)*(r+1), &
                            (p+1)*(q+1)*(r+1), (p+1)*(q+1)*(r+1)-p, &
                            (p+1)*(q+1)*(r+1)-p, (p+1)*(q+1)*r+1 ], [2,12] )

        element%facevert = element%sidvert

!       the following are NOT available: DON'T USE THEM
!          element%sidnod = []
!          element%sidinternnod = []
!          element%edgenod = []
!          element%facenod = []

        element%xc = 0

      case(104) ! triangle high-order element

        p = element%p(1,1)

        element%ndimr = 2           ! dimension of reference space
        element%numshapenod = (p+1)*(p+2)/2 ! number of `shape' nodal points
                                            ! for mapping
        element%numsides = 3                ! number of sides
        element%sidnumvert = 2              ! number of vertices of a side
        element%numinternnod = (p-2)*(p-1)/2  ! number internal nodes
        element%numedges = 3                ! number of edges
        element%numfaces = 1                ! number of faces
        element%facenumvert = 3             ! number of vertices of a face
        element%facenumnod = (p+1)*(p+2)/2  ! number of nodes of a face

!       the following are NOT available: DON'T USE THEM
        element%sidnumnod = 0       ! number of nodes of a side
        element%sidnuminternnod = 0 ! number internal nodes on a side
        element%edgenumnod = 0      ! number of nodes of an edge

        call allocate_element_arrays

        element%sidvert = reshape ( [ 1, p+1, p+1, (p+1)*(p+2)/2, &
                                       (p+1)*(p+2)/2, 1 ], [2,3] )

        element%internnod =  [ ((i+(2*p+4-j)*(j-1)/2,i=2,p+1-j),j=2,p-1) ]

        element%edgevert = element%sidvert

        element%facevert = reshape ( [ 1, p+1, (p+1)*(p+2)/2 ], [3,1] )

        element%facenod = reshape ( [ (i,i=1,element%facenumnod) ], &
                                        [element%facenumnod,1] )

!       the following are NOT available: DON'T USE THEM
!          element%sidnod = []
!          element%sidinternnod = []
!          element%edgenod = []

       element%xc = 1._dp/3

      case default

        write(*,'(/a,i0/)') &
          'Error fill_element: Wrong or not implemented element shape: ', &
          element%elshape
        stop

    end select

  contains

    subroutine allocate_element_arrays

      allocate( element%sidvert(element%sidnumvert,element%numsides) )
      allocate( element%sidnod(element%sidnumnod,element%numsides) )
      allocate( element%internnod(element%numinternnod) )
      allocate( element%sidinternnod(element%sidnuminternnod,element%numsides) )
      allocate( element%edgevert(2,element%numedges) )
      allocate( element%edgenod(element%edgenumnod,element%numedges) )
      allocate( element%facevert(element%facenumvert,element%numfaces) )
      allocate( element%facenod(element%facenumnod,element%numfaces) )
      allocate( element%xc(element%ndimr) )

    end subroutine allocate_element_arrays

  end subroutine fill_element_core


! Fill elements for figplot

  subroutine fill_element_figplot ( element )

    type(element_t), intent(inout) :: element

!   Fill element type info for single element group for figplot

    select case ( element%elshape )

      case(3,10) ! three-node and four-node triangle

        element%numsubdiv = 1       ! subdivision: number of triangles

        element%subtopol = reshape ( [ 1,2,3 ], [3,1] )

      case(4,7,33) ! six-node and seven-node triangle

!       figplot

        element%numsubdiv = 4       ! subdivision: number of triangles

        element%subtopol = reshape ( [ 1,2,6, 2,3,4, 2,4,6, 6,4,5 ], [3,4] )

      case(5,9) ! four-node and five-node quadrilateral

!       figplot

        element%numsubdiv = 2       ! subdivision: number of triangles

        element%subtopol = reshape ( [ 1,2,3, 1,3,4 ], [3,2] )

      case(6,34) ! nine-node quadrilateral

!       figplot

        element%numsubdiv = 8       ! subdivision: number of triangles

        element%subtopol = reshape ( [ 1,2,9, 1,9,8, 2,3,4, 2,4,9, &
                                       9,4,5, 9,5,6, 8,9,6, 8,6,7 ], [3,8] )

      case(30) ! eight-node quadrilateral

!       figplot

        element%numsubdiv = 6       ! subdivision: number of triangles

        element%subtopol = reshape ( [ 1,2,8, 2,6,8, 8,6,7, &
                                       2,3,4, 2,4,6, 4,5,6 ], [3,6] )

      case default

        element%numsubdiv = 0       ! subdivision: number of triangles

        allocate( element%subtopol(0,0) )

    end select

  end subroutine fill_element_figplot


! Fill elements VTK

  subroutine fill_element_vtk ( element )

    use limits_m, only: USE_VTK_ELEMENTS

    type(element_t), intent(inout) :: element

!   Fill element type info for single element group for VTK

    select case ( element%elshape )

      case(1) ! two-node line element

        element%vtk_numnod = 2      ! number of vtk nodal points
        element%vtk_numelm = 1      ! number of vtk subelements
        element%vtk_elshape = 3     ! shape number of vtk format output

        element%vtk_index = reshape ( [ 1, 2 ], [2,1] )

      case(2,32) ! three-node line element

        element%vtk_numnod = 3      ! number of vtk nodal points
        element%vtk_numelm = 1      ! number of vtk subelements
        element%vtk_elshape = 21    ! shape number of vtk format output

        element%vtk_index = reshape ( [ 1, 3, 2 ], [3,1] )

      case(3,10) ! three-node and four-node triangle

        element%vtk_numnod = 3      ! number of vtk nodal points
        element%vtk_numelm = 1      ! number of vtk subelements
        element%vtk_elshape = 5     ! shape number of vtk format output

        element%vtk_index = reshape ( [ 1, 2, 3 ], [3,1] )

      case(4,33) ! six-node node triangle

        element%vtk_numnod = 6      ! number of vtk nodal points
        element%vtk_numelm = 1      ! number of vtk subelements
        element%vtk_elshape = 22    ! shape number of vtk format output

        element%vtk_index = reshape ( [ 1, 3, 5, 2, 4, 6 ], [6,1] )

      case(7) ! seven-node triangle

        if ( USE_VTK_ELEMENTS ) then

          element%vtk_numnod = 7      ! number of vtk nodal points
          element%vtk_numelm = 1      ! number of vtk subelements
          element%vtk_elshape = 34    ! shape number of vtk format output

          element%vtk_index = reshape ( [ 1, 3, 5, 2, 4, 6, 7 ], [7,1] )

        else

          element%vtk_numnod = 6      ! number of vtk nodal points
          element%vtk_numelm = 1      ! number of vtk subelements
          element%vtk_elshape = 22    ! shape number of vtk format output

          element%vtk_index = reshape ( [ 1, 3, 5, 2, 4, 6 ], [6,1] )

        end if

      case(5,9) ! four-node and five-node quadrilateral

        element%vtk_numnod = 4     ! number of vtk nodal points
        element%vtk_numelm = 1     ! number of vtk subelements
        element%vtk_elshape = 9    ! shape number of vtk format output

        element%vtk_index = reshape ( [ 1, 2, 3, 4 ], [4,1] )

      case(6,34) ! nine-node quadrilateral

        if ( USE_VTK_ELEMENTS ) then

          element%vtk_numnod = 9      ! number of vtk nodal points
          element%vtk_numelm = 1      ! number of vtk subelements
          element%vtk_elshape = 28    ! shape number of vtk format output

          element%vtk_index = reshape ( [ 1,3,5,7,2,4,6,8,9 ], [9,1] )

        else

          element%vtk_numnod = 4      ! number of vtk nodal points
          element%vtk_numelm = 4      ! number of vtk subelements
          element%vtk_elshape = 9     ! shape number of vtk format output

          element%vtk_index = reshape ( [ 1,2,9,8, 2,3,4,9, &
                                          8,9,6,7, 9,4,5,6 ], [4,4] )

        end if

      case(11,18) ! four-node and five-node tetrahedron

        element%vtk_numnod = 4      ! number of vtk nodal points
        element%vtk_numelm = 1      ! number of vtk subelements
        element%vtk_elshape = 10    ! shape number of vtk format output

        element%vtk_index = reshape ( [ 1, 2, 3, 4 ], [4,1] )

      case(12,15,16,35) ! 10-node, 14-node and 15-node tetrahedron

        element%vtk_numnod = 10     ! number of vtk nodal points
        element%vtk_numelm = 1      ! number of vtk subelements
        element%vtk_elshape = 24    ! shape number of vtk format output

        element%vtk_index = reshape ( [ 1, 3, 5, 10, 2, 4, 6, 7, 8, 9 ], &
                                        [10,1] )

      case(13,17) ! eight-node and nine-node hexahedron

        element%vtk_numnod = 8      ! number of vtk nodal points
        element%vtk_numelm = 1      ! number of vtk subelements
        element%vtk_elshape = 12    ! shape number of vtk format output

        element%vtk_index = reshape ( [ 1, 2, 3, 4, 5, 6, 7, 8 ], [8,1] )

      case(14,36) ! 27-node hexahedron

        if ( USE_VTK_ELEMENTS ) then

          element%vtk_numnod = 27    ! number of vtk nodal points
          element%vtk_numelm = 1     ! number of vtk subelements
          element%vtk_elshape = 29   ! shape number of vtk format output

          element%vtk_index = reshape ( [ 1,3,9,7, 19,21,27,25, &
                                          2,6,8,4, 20,24,26,22, &
                                          10,12,18,16, 13,15,11,17,5,23, 14 ], &
                                          [27,1] )

        else

          element%vtk_numnod = 8     ! number of vtk nodal points
          element%vtk_numelm = 8     ! number of vtk subelements
          element%vtk_elshape = 12   ! shape number of vtk format output

          element%vtk_index = reshape ( [ 1,2,5,4,     10,11,14,13, &
                                          2,3,6,5,     11,12,15,14, &
                                          4,5,8,7,     13,14,17,16, &
                                          5,6,9,8,     14,15,18,17, &
                                          10,11,14,13, 19,20,23,22, &
                                          11,12,15,14, 20,21,24,23, &
                                          13,14,17,16, 22,23,26,25, &
                                          14,15,18,17, 23,24,27,26 ], [8,8] )
        end if

      case(30) ! eight-node quadrilateral

        element%vtk_numnod = 8      ! number of vtk nodal points
        element%vtk_numelm = 1      ! number of vtk subelements
        element%vtk_elshape = 23    ! shape number of vtk format output

        element%vtk_index = reshape ( [ 1,3,5,7,2,4,6,8 ], [8,1] )

      case(31) ! twenty-node hexahedron

        element%vtk_numnod = 20      ! number of vtk nodal points
        element%vtk_numelm = 1       ! number of vtk subelements
        element%vtk_elshape = 25     ! shape number of vtk format output

        element%vtk_index = reshape ( [ 1,3,8,6, 13,15,20,18, &
                                        2,5,7,4, 14,17,19,16, &
                                        9,10,12,11 ], [20,1] )

      case(41) ! six-node prism (wedge)

        element%vtk_numnod = 6       ! number of vtk nodal points
        element%vtk_numelm = 1       ! number of vtk subelements
        element%vtk_elshape = 13     ! shape number of vtk format output

        element%vtk_index = reshape ( [ 4,5,6, 1,2,3 ], [6,1] )

      case(42) ! fifteen-node prism (wedge)

        element%vtk_numnod = 15      ! number of vtk nodal points
        element%vtk_numelm = 1       ! number of vtk subelements
        element%vtk_elshape = 26     ! shape number of vtk format output

        element%vtk_index = reshape ( [ 10,12,14, 1,3,5, 11,13,15, &
                                        2,4,6, 7,8,9 ], [15,1] )

      case(43) ! eighteen-node prism (wedge)

        element%vtk_numnod = 18      ! number of vtk nodal points
        element%vtk_numelm = 1       ! number of vtk subelements
        element%vtk_elshape = 32     ! shape number of vtk format output

        element%vtk_index = reshape ( [ 13,15,17, 1,3,5, 14,16,18, &
                                        2,4,6, 7,9,11, 8,10,12 ], [18,1] )

      case(51) ! five-node pyramid

        element%vtk_numnod = 5       ! number of vtk nodal points
        element%vtk_numelm = 1       ! number of vtk subelements
        element%vtk_elshape = 14     ! shape number of vtk format output

        element%vtk_index = reshape ( [ 1,2,3,4,5 ], [5,1] )

      case(52) ! thirteen-node pyramid

        element%vtk_numnod = 13      ! number of vtk nodal points
        element%vtk_numelm = 1       ! number of vtk subelements
        element%vtk_elshape = 27     ! shape number of vtk format output

        element%vtk_index = reshape ( [ 1,3,5,7, 13, 2,4,6,8, 9,10,11,12 ], &
                                      [13,1] )

      case(53) ! fourteen-node pyramid

        element%vtk_numnod = 10      ! number of vtk nodal points
        element%vtk_numelm = 2       ! number of vtk subelements
        element%vtk_elshape = 24     ! shape number of vtk format output

        element%vtk_index = reshape ( [ 3,5,7, 14, 4,6,9, 11,12,13, &
                                        1,3,7, 14, 2,9,8, 10,11,13 ], [10,2] )

      case default ! elements not supported

        element%vtk_numnod = 0
        element%vtk_numelm = 0
        element%vtk_elshape = 0

        allocate( element%vtk_index(0,0) )

    end select

  end subroutine fill_element_vtk


! Fill elements for tecplot

  subroutine fill_element_tecplot ( element )

    type(element_t), intent(inout) :: element

!   Fill element type info for single element group for tecplot

    select case ( element%elshape )

      case(1) ! two-node line element

        element%tec_numnod = 2      ! number of tecplot nodal points
        element%tec_numelm = 1      ! number of tecplot subelements
        element%tec_elshape = 'FELINESEG' ! element shape of tecplot output

        element%tec_index = reshape ( [ 1, 2 ], [2,1] )

      case(2,32) ! three-node line element

        element%tec_numnod = 2      ! number of tecplot nodal points
        element%tec_numelm = 2      ! number of tecplot subelements
        element%tec_elshape = 'FELINESEG' ! element shape of tecplot output

        element%tec_index =  reshape ( [ 1,2, 2,3 ], [2,2] )

      case(3,10) ! three-node and four-node triangle

        element%tec_numnod = 3      ! number of tecplot nodal points
        element%tec_numelm = 1      ! number of tecplot subelements
        element%tec_elshape = 'FETRIANGLE' ! element shape of tecplot output

        element%tec_index = reshape ( [ 1,2,3 ], [3,1] )

      case(4,7,33) ! six-node and seven-node triangle

        element%tec_numnod = 3      ! number of tecplot nodal points
        element%tec_numelm = 4      ! number of tecplot subelements
        element%tec_elshape = 'FETRIANGLE' ! element shape of tecplot output

        element%tec_index = reshape ( [ 1,2,6, 2,3,4, 2,4,6, 6,4,5 ], &
                                      [3,4])

      case(5,9) ! four-node and five-node quadrilateral

        element%tec_numnod = 4      ! number of tecplot nodal points
        element%tec_numelm = 1      ! number of tecplot subelements
        element%tec_elshape = 'FEQUADRILATERAL' ! element shape of tecplot

        element%tec_index = reshape ( [ 1,2,3,4 ], [4,1] )

      case(6,34) ! nine-node quadrilateral

        element%tec_numnod = 4      ! number of tecplot nodal points
        element%tec_numelm = 4      ! number of tecplot subelements
        element%tec_elshape = 'FEQUADRILATERAL' ! element shape of tecplot

        element%tec_index = reshape ( [ 1,2,9,8, 2,3,4,9, &
                                        8,9,6,7, 9,4,5,6 ], [4,4] )

      case(11,18) ! four-node and five-node tetrahedron

        element%tec_numnod = 4      ! number of tecplot nodal points
        element%tec_numelm = 1      ! number of tecplot subelements
        element%tec_elshape = 'FETETRAHEDRON' ! element shape of tecplot

        element%tec_index = reshape ( [ 4,1,2,3 ], [4,1] )

      case(12,15,16,35) ! 10-node, 14-node and 15-node tetrahedron

        element%tec_numnod = 4      ! number of tecplot nodal points
        element%tec_numelm = 8      ! number of tecplot subelements
        element%tec_elshape = 'FETETRAHEDRON' ! element shape of tecplot

        element%tec_index =  reshape ( [ 7,1,2,6, 8,2,3,4,  &
                                         9,4,5,6, 10,7,8,9, &
                                         7,2,8,9, 7,2,9,6,  &
                                         9,2,8,4, 9,2,4,6 ], [4,8] )

      case(13,17) ! eight-node and nine-node hexahedron

        element%tec_numnod = 8      ! number of tecplot nodal points
        element%tec_numelm = 1      ! number of tecplot subelements
        element%tec_elshape = 'FEBRICK' ! element shape of tecplot

        element%tec_index =  reshape ( [ 1, 2, 3, 4, 5, 6, 7, 8 ], [8,1] )

      case(14,36) ! 27-node hexahedron

        element%tec_numnod = 8      ! number of tecplot nodal points
        element%tec_numelm = 8      ! number of tecplot subelements
        element%tec_elshape = 'FEBRICK' ! element shape of tecplot

        element%tec_index = reshape ( [ 1,2,5,4,     10,11,14,13, &
                                        2,3,6,5,     11,12,15,14, &
                                        4,5,8,7,     13,14,17,16, &
                                        5,6,9,8,     14,15,18,17, &
                                        10,11,14,13, 19,20,23,22, &
                                        11,12,15,14, 20,21,24,23, &
                                        13,14,17,16, 22,23,26,25, &
                                        14,15,18,17, 23,24,27,26 ], [8,8] )

      case(30) ! eight-node quadrilateral

        element%tec_numnod = 3      ! number of tecplot nodal points
        element%tec_numelm = 6      ! number of tecplot subelements
        element%tec_elshape = 'FETRIANGLE' ! element shape of tecplot output

        element%tec_index = reshape ( [ 1,2,8, 2,6,8, 8,6,7, &
                                        2,3,4, 2,4,6, 4,5,6 ], [3,6] )

      case default  ! element not supported

        element%tec_numnod = 0
        element%tec_numelm = 0
        element%tec_elshape = ''

        allocate( element%tec_index(0,0) )

    end select

  end subroutine fill_element_tecplot


! Fill elements Gmsh

  subroutine fill_element_gmsh ( element )

    use meshgen_gmsh_index_m

    type(element_t), intent(inout) :: element

    integer :: p
    integer, dimension(:), allocatable :: eltypes

!   Fill element type info for single element group for Gmsh

    select case ( element%elshape )

      case(1) ! two-node line element

        element%gmsh_numnod = 2      ! number of gmsh nodal points
        element%gmsh_elshape = 1     ! shape number of gmsh format output

        element%gmsh_index = [ 1, 2 ]

      case(2,32) ! three-node line element

        element%gmsh_numnod = 3      ! number of gmsh nodal points
        element%gmsh_elshape = 8     ! shape number of gmsh format output

        element%gmsh_index = [ 1, 3, 2 ]

      case(3,10) ! three-node and four-node triangle

        element%gmsh_numnod = 3      ! number of gmsh nodal points
        element%gmsh_elshape = 2     ! shape number of gmsh format output

        element%gmsh_index = [ 1, 2, 3 ]

      case(4,7,33) ! six-node and seven-node triangle

        element%gmsh_numnod = 6      ! number of gmsh nodal points
        element%gmsh_elshape = 9     ! shape number of gmsh format output

        element%gmsh_index = [ 1, 3, 5, 2, 4, 6 ]

      case(5,9) ! four-node and five-node quadrilateral

        element%gmsh_numnod = 4     ! number of gmsh nodal points
        element%gmsh_elshape = 3    ! shape number of gmsh format output

        element%gmsh_index = [ 1, 2, 3, 4 ]

      case(6,34) ! nine-node quadrilateral

        element%gmsh_numnod = 9      ! number of gmsh nodal points
        element%gmsh_elshape = 10    ! shape number of gmsh format output

        element%gmsh_index = [ 1, 3, 5, 7, 2, 4, 6, 8, 9 ]

      case(11,18) ! four-node and five-node tetrahedron

        element%gmsh_numnod = 4      ! number of gmsh nodal points
        element%gmsh_elshape = 4    ! shape number of gmsh format output

        element%gmsh_index = [ 1, 2, 3, 4 ]

      case(12,15,16,35) ! 10-node, 14-node and 15-node tetrahedron

        element%gmsh_numnod = 10     ! number of gmsh nodal points
        element%gmsh_elshape = 11    ! shape number of gmsh format output

        element%gmsh_index = [ 1, 3, 5, 10, 2, 4, 6, 7, 9, 8 ]

      case(13,17) ! eight-node and nine-node hexahedron

        element%gmsh_numnod = 8     ! number of gmsh nodal points
        element%gmsh_elshape = 5     ! shape number of gmsh format output

        element%gmsh_index = [ 1, 2, 3, 4, 5, 6, 7, 8 ]

      case(14,36) ! 27-node hexahedron

        element%gmsh_numnod = 27    ! number of gmsh nodal points
        element%gmsh_elshape = 12   ! shape number of gmsh format output

        element%gmsh_index = [ 1, 3, 9, 7, 19, 21, 27, 25, 2, &
                               4, 10, 6, 12, 8, 18, 16, 20, 22, &
                               24, 26, 5, 11, 13, 15, 17, 23, 14 ]

      case(30) ! eight-node quadrilateral

        element%gmsh_numnod = 8      ! number of gmsh nodal points
        element%gmsh_elshape = 16    ! shape number of gmsh format output

        element%gmsh_index = [ 1, 3, 5, 7, 2, 4, 6, 8 ]

      case(31) ! twenty-node hexahedron

        element%gmsh_numnod = 20      ! number of gmsh nodal points
        element%gmsh_elshape = 17     ! shape number of gmsh format output

        element%gmsh_index = [ 1, 3, 8, 6, 13, 15, 20, 18, 2, 4, &
                               9, 5, 10, 7, 12, 11, 14, 16, 17, 19 ]

      case(41) ! six-node prism

        element%gmsh_numnod = 6      ! number of gmsh nodal points
        element%gmsh_elshape = 6     ! shape number of gmsh format output

        element%gmsh_index = [ 1, 2, 3, 4, 5, 6 ]

      case(42) ! fifteen-node prism

        element%gmsh_numnod = 15      ! number of gmsh nodal points
        element%gmsh_elshape = 18     ! shape number of gmsh format output

        element%gmsh_index = [ 1,3,5, 10,12,14, 2,6,7, 4,8,9, 11,15,13 ]

      case(43) ! eighteen-node prism

        element%gmsh_numnod = 18      ! number of gmsh nodal points
        element%gmsh_elshape = 13     ! shape number of gmsh format output

        element%gmsh_index = [ 1,3,5, 13,15,17, 2,6,7, &
                               4,9,11, 14,18,16, 8,12,10 ]

      case(51) ! five-node pyramid

        element%gmsh_numnod = 5      ! number of gmsh nodal points
        element%gmsh_elshape = 7     ! shape number of gmsh format output

        element%gmsh_index = [ 1, 2, 3, 4, 5 ]

      case(52) ! thirteen-node pyramid

        element%gmsh_numnod = 13      ! number of gmsh nodal points
        element%gmsh_elshape = 19     ! shape number of gmsh format output

        element%gmsh_index = [ 1,3,5,7, 13, 2,8,9, 4,10,6, 11,12 ]

      case(53) ! fourteen-node pyramid

        element%gmsh_numnod = 14      ! number of gmsh nodal points
        element%gmsh_elshape = 14     ! shape number of gmsh format output

        element%gmsh_index = [ 1,3,5,7, 14, 2,8,10, 4,11,6, 12,13, 9 ]

      case(101,105) ! high order line element (including macro variant)

        p = element%p(1,1)

        if ( p <= 10 ) then

          element%gmsh_numnod = p+1     ! number of gmsh nodal points

          eltypes = [ 1, 8, 26, 27, 28, 62, 63, 64, 65, 66 ]

          element%gmsh_elshape = eltypes(p) ! shape number of gmsh format output

          allocate( element%gmsh_index(element%gmsh_numnod) )

          call gmsh_index_line ( p, element%gmsh_index )

        else

!         element not supported

          element%gmsh_numnod = 0
          element%gmsh_elshape = 0

          allocate( element%gmsh_index(0) )

        end if

      case(102,106) ! high order quadrilateral element (including macro variant)

        p = element%p(1,1)

        if ( p <= 10 ) then

          element%gmsh_numnod = (p+1)**2     ! number of gmsh nodal points

          eltypes = [ 3, 10, 36, 37, 38, 47, 48, 49, 50, 51 ]

          element%gmsh_elshape = eltypes(p) ! shape number of gmsh format output

          allocate( element%gmsh_index(element%gmsh_numnod) )

          call gmsh_index_quadrilateral ( p, element%gmsh_index )

        else

!         element not supported

          element%gmsh_numnod = 0
          element%gmsh_elshape = 0

          allocate( element%gmsh_index(0) )

        end if

      case(103,107) ! high order hexahedron element (including macro variant)

        p = element%p(1,1)

        if ( p <= 9 ) then

          element%gmsh_numnod = (p+1)**3     ! number of gmsh nodal points

          eltypes = [ 5, 12, 92, 93, 94, 95, 96, 97, 98 ]

          element%gmsh_elshape = eltypes(p) ! shape number of gmsh format output

          allocate( element%gmsh_index(element%gmsh_numnod) )

          call gmsh_index_hexahedron ( p, element%gmsh_index )

        else

!         element not supported

          element%gmsh_numnod = 0
          element%gmsh_elshape = 0

          allocate( element%gmsh_index(0) )

        end if

      case(104) ! high order triangle element

        p = element%p(1,1)

        if ( p <= 10 ) then

          element%gmsh_numnod = (p+1)*(p+2)/2     ! number of gmsh nodal points

          eltypes = [ 2, 9, 21, 23, 25, 42, 43, 44, 45, 46 ]

          element%gmsh_elshape = eltypes(p) ! shape number of gmsh format output

          allocate( element%gmsh_index(element%gmsh_numnod) )

          call gmsh_index_triangle ( p, element%gmsh_index )

        else

!         element not supported

          element%gmsh_numnod = 0
          element%gmsh_elshape = 0

          allocate( element%gmsh_index(0) )

        end if

      case default ! elements not supported

        element%gmsh_numnod = 0
        element%gmsh_elshape = 0

        allocate( element%gmsh_index(0) )

    end select

  end subroutine fill_element_gmsh


! Fill elements for gmsh parsed

  subroutine fill_element_gmsh_parsed ( element )

    type(element_t), intent(inout) :: element

!   Fill element type info for single element group for gmsh parsed

    select case ( element%elshape )

      case(1) ! two-node line element

        element%gmsh_numnod_parsed = 2    ! number of gmsh parsed nodal points
        element%gmsh_numelm_parsed = 1    ! number of gmsh parsed subelements
        element%gmsh_elshape_parsed = 'L' ! element shape of gmsh parsed output

        element%gmsh_index_parsed = reshape ( [ 1, 2 ], [2,1] )

      case(2,32) ! three-node line element

        element%gmsh_numnod_parsed = 2    ! number of gmsh parsed nodal points
        element%gmsh_numelm_parsed = 2    ! number of gmsh parsed subelements
        element%gmsh_elshape_parsed = 'L' ! element shape of gmsh parsed output

        element%gmsh_index_parsed =  reshape ( [ 1,2, 2,3 ], [2,2] )

      case(3,10) ! three-node and four-node triangle

        element%gmsh_numnod_parsed = 3    ! number of gmsh parsed nodal points
        element%gmsh_numelm_parsed = 1    ! number of gmsh parsed subelements
        element%gmsh_elshape_parsed = 'T' ! element shape of gmsh parsed output

        element%gmsh_index_parsed = reshape ( [ 1,2,3 ], [3,1] )

      case(4,7,33) ! six-node and seven-node triangle

        element%gmsh_numnod_parsed = 3    ! number of gmsh parsed nodal points
        element%gmsh_numelm_parsed = 4    ! number of gmsh parsed subelements
        element%gmsh_elshape_parsed = 'T' ! element shape of gmsh parsed output

        element%gmsh_index_parsed = reshape ( [ 1,2,6, 2,3,4, 2,4,6, 6,4,5 ], &
                                      [3,4])

      case(5,9) ! four-node and five-node quadrilateral

        element%gmsh_numnod_parsed = 4    ! number of gmsh parsed nodal points
        element%gmsh_numelm_parsed = 1    ! number of gmsh parsed subelements
        element%gmsh_elshape_parsed = 'Q' ! element shape of gmsh parsed output

        element%gmsh_index_parsed = reshape ( [ 1,2,3,4 ], [4,1] )

      case(6,34) ! nine-node quadrilateral

        element%gmsh_numnod_parsed = 4    ! number of gmsh parsed nodal points
        element%gmsh_numelm_parsed = 4    ! number of gmsh parsed subelements
        element%gmsh_elshape_parsed = 'Q' ! element shape of gmsh parsed output

        element%gmsh_index_parsed = reshape ( [ 1,2,9,8, 2,3,4,9, &
                                         8,9,6,7, 9,4,5,6 ], [4,4] )

      case(11,18) ! four-node and five-node tetrahedron

        element%gmsh_numnod_parsed = 4    ! number of gmsh parsed nodal points
        element%gmsh_numelm_parsed = 1    ! number of gmsh parsed subelements
        element%gmsh_elshape_parsed = 'S' ! element shape of gmsh parsed output

        element%gmsh_index_parsed = reshape ( [ 1,2,3,4 ], [4,1] )

      case(12,15,16,35) ! 10-node, 14-node and 15-node tetrahedron

        element%gmsh_numnod_parsed = 4    ! number of gmsh parsed nodal points
        element%gmsh_numelm_parsed = 8    ! number of gmsh parsed subelements
        element%gmsh_elshape_parsed = 'S' ! element shape of gmsh parsed

        element%gmsh_index_parsed = reshape ( [ 1,2,6,7, 2,3,4,8, &
                                                4,5,6,9, 7,8,9,10, &
                                                2,8,9,7, 2,9,6,7,  &
                                                2,8,4,9, 2,4,6,9 ], [4,8] )

      case(13,17) ! eight-node and nine-node hexahedron

        element%gmsh_numnod_parsed = 8    ! number of gmsh parsed nodal points
        element%gmsh_numelm_parsed = 1    ! number of gmsh parsed subelements
        element%gmsh_elshape_parsed = 'H' ! element shape of gmsh parsed

        element%gmsh_index_parsed = reshape ( [ 1,2,3,4,5,6,7,8 ], [8,1] )

      case(14,36) ! 27-node hexahedron

        element%gmsh_numnod_parsed = 8    ! number of gmsh parsed nodal points
        element%gmsh_numelm_parsed = 8    ! number of gmsh parsed subelements
        element%gmsh_elshape_parsed = 'H' ! element shape of gmsh parsed output

        element%gmsh_index_parsed = reshape ( [ 1,2,5,4,     10,11,14,13, &
                                                2,3,6,5,     11,12,15,14, &
                                                4,5,8,7,     13,14,17,16, &
                                                5,6,9,8,     14,15,18,17, &
                                                10,11,14,13, 19,20,23,22, &
                                                11,12,15,14, 20,21,24,23, &
                                                13,14,17,16, 22,23,26,25, &
                                                14,15,18,17, 23,24,27,26 ], &
                                              [8,8] )

      case(30) ! eight-node quadrilateral

        element%gmsh_numnod_parsed = 3    ! number of gmsh parsed nodal points
        element%gmsh_numelm_parsed = 6    ! number of gmsh parsed subelements
        element%gmsh_elshape_parsed = 'T' ! element shape of gmsh parsed output

        element%gmsh_index_parsed = reshape ( [ 1,2,8, 2,6,8, 8,6,7, &
                                                2,3,4, 2,4,6, 4,5,6 ], [3,6] )

      case default  ! element not supported

        element%gmsh_numnod_parsed = 0
        element%gmsh_numelm_parsed = 0
        element%gmsh_elshape_parsed = ''

        allocate( element%gmsh_index_parsed(0,0) )

    end select

  end subroutine fill_element_gmsh_parsed

end module meshgen_parts_m
