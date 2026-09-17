
! Copyright (C) 2007-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Routines, for the system of equations involving connections

module system_connection_m

  use kind_defs_m
  use mesh_m
  use problem_defs_m
  use sparse_m
  use system_defs_m
  use system_matrix_m
  use element_defs_m
  use set_optional_m
  use misc_m, only: sort

  implicit none


contains


! Create system matrix structure of connections (ia not yet accumulated)

  subroutine create_sysmatrix_structure_connection ( sysmatrix, mesh, problem, &
    symmetric )

    type(sysmatrix_t), intent(inout) :: sysmatrix
    type(mesh_t),      intent(in)  :: mesh
    type(problem_t),   intent(in)  :: problem

!   if symmetric is present and has the values .true. the matrix is assumed to
!   be symmetric and only the non-zero elements of the upper-triangle
!   (including the diagonal) of the matrix of unknowns (Suu) are stored.
!   default: symmetric=.false.
!   NOTE: this has no effect if sysmatrix has been initialized already.
    logical, intent(in), optional :: symmetric

!   This routine fills the arrays ia with the number of non-zeros in each row:
!
!       ia(i+1) = number of non-zeros in row i
!
!   ia is not yet accumulated and may still be adapted (add more non-zeros)
!   Here the structure of the connections is created.
!
!   Connections are basically
!    a. connections between elements defined by
!       1. connecting geometries or elementsets.
!       2. intersections by points (nodes or integration points) in objects.
!    b. connections between nodes in points, nodesets, geometries or
!       elementsets.
!    c. connections between elements defined by intersections by nodes in
!       objects and nodes in geometries, nodesets or elementsets.
!   Points must be connected to points, curves to curves etc., with the
!   exceptions that geometries (curves, surfaces, volumes) can be connected to
!   elementsets and nodesets, and objects can be connected to geometries
!   (curves, surfaces, volumes), nodesets and elementsets.
!
!   If there is more than one connection, you need to be careful:
!     a. If connections have identical values for BOTH physq1 and physq2,
!        (i.e. conn1%physq1=conn2%physq1 AND conn1%physq2=conn2%physq2)
!        or different values for BOTH physq1 and physq2,
!        (i.e. conn1%physq1/=conn2%physq1 AND conn1%physq2/=conn2%physq2)
!        there is no problem. In the first case the connections are simply
!        merged into basically one single connection. In the second case
!        the degrees of freedom are separated by default.
!        If layers are defined, the above also applies to layers.
!     b. If NOT a. then be sure that the connections are FULLY SEPARATED, i.e.
!        the connection must not connect to degrees of freedom within the
!        standard FEM influence range.

!   It is assumed that the elements/nodes in geometries or elementsets or
!   intersected by the objects all already contribute to the system matrix.
!   This means, that the connections internal to an element are already
!   accounted for. In this routine only the "off-diagonal" blocks are
!   accounted for. The diagonal blocks fit in the standard FEM system.

    logical :: done(problem%numconnections)
    integer :: numess, numund, conn
    integer :: nnodes, nnodes1, nnodes2, nodenr
    integer :: object, physq, elem2
    integer :: node, ndnr, i, ndnr0
    integer :: numnod, deg_r, row, deg_c1, deg_c2, deg_c, rowp
    integer :: colp, col, maxnumnodes, n1, n2, numnod1, numnod0, conn1, nn
    integer :: dof, layer, dofr

    integer, dimension(mesh%maxnodnumnod+1) :: nodes1, worksc1, worksc2
    integer, dimension(problem%maxnoddegfd) :: pos, posr
    type(int_array_1d_t), dimension(problem%numconnections) :: worknod
    integer, dimension(:), allocatable :: utype, work, nodes0, nodes, dofc
    integer, dimension(:,:), allocatable ::worknumnod, work1, posc

    allocate ( utype(problem%numdegfd), work(problem%numdegfd) )
    allocate ( worknumnod(mesh%nnodes+1,problem%numconnections) )
    allocate ( work1(mesh%nnodes,problem%numconnections) )

!   NOTE: only one-way connections are counted, i.e. connections of degrees
!   of freedom on the first side to degrees of freedom on the second side. The
!   reverse connection is taken into account by adding the "transpose" element
!   in the matrix also.

    call check ( mesh, 'create_sysmatrix_structure_connection' )
    call check ( problem, 'create_sysmatrix_structure_connection', mesh )

    if ( sysmatrix%finalized ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysmatrix_structure_connection: ', &
        ' sysmatrix has already been finalized '
      stop
    end if

!   fill temporary utype array that gives the type of the degree of freedom
!      1: unknown  3: prescribed

    numess = problem%numessdegfd
    numund = problem%numundegfd

    where ( problem%degfdperm(:,2) > numund )
      utype = 3
    else where
      utype = 1
    end where

!   initialize?

    if ( .not. sysmatrix%initialized_structure ) then

!     initialize system_matrix

      call initialize_sysmatrix_structure ( sysmatrix, numund, numess )

!     set symmetry of the matrix

      sysmatrix%symmetric = set_optional ( variable=symmetric, default=.false. )

    end if

    worknumnod = 0

!   count number of connecting nodes

    call nodes_connected ( pass=1 )

    maxnumnodes = maxval ( sum ( worknumnod, dim=2 ) )

!   accumulate worknumnod

    do node = 2, mesh%nnodes
      worknumnod(node+1,:) = worknumnod(node,:) + worknumnod(node+1,:)
    end do

!   fill worknod

    do conn = 1, problem%numconnections
      allocate ( worknod(conn)%a( worknumnod( mesh%nnodes + 1, conn ) ) )
    end do

!   fill connecting nodes

    work1 = 0

    call nodes_connected ( pass=2 )

!   allocate arrays

    allocate ( nodes0(maxnumnodes), nodes(maxnumnodes) )
    allocate ( dofc(maxnumnodes), posc(maxnumnodes,problem%maxnoddegfd) )

!   scan all nodes for degrees of freedom and look for connections

    work = 0

    do nodenr = 1, mesh%nnodes

!     cycle node if no degrees
      if ( all ( worknumnod(nodenr+1,:) - worknumnod(nodenr,:) == 0 ) ) cycle

      done = .false.

!     fill degrees already counted for in standard FEM for this node

!     nodes are current node + all connections

      nodes1(1) = nodenr
      n1 = mesh%nodnumnod(nodenr)
      n2 = mesh%nodnumnod(nodenr+1)
      numnod1 = n2 - n1
      nodes1(2:numnod1+1) = mesh%nodnod(n1+1:n2)

      do i = 1, numnod1 + 1

        node = nodes1(i)

        deg_c1 = problem%nodnumdegfd ( node ) + 1
        deg_c2 = problem%nodnumdegfd ( node + 1 )

!       store deg_c1 and deg_c2 for later use

        worksc1(i) = deg_c1
        worksc2(i) = deg_c2

      end do

!     loop over all connections with possible overlap

      do conn = 1, problem%numconnections

        if ( done(conn) ) cycle

!       row degrees

        physq = problem%connections(conn)%physq1
        layer = problem%connections(conn)%layer1

        if ( physq > 0 ) then
          call pos_array_local_node ( problem, nodenr, dof, pos, [physq], &
            layer )
        else
          call pos_array_local_node ( problem, nodenr, dof, pos, layer=layer )
        end if

!       store and change to absolute positions

        dofr = dof
        posr(1:dof) = problem%nodnumdegfd(nodenr) + pos(1:dof)

!       collect nodes

        numnod0 = worknumnod(nodenr+1,conn) - worknumnod(nodenr,conn)
        nodes0(1:numnod0) = worknod(conn)%a( &
             [ ( worknumnod(nodenr,conn) + i, i=1,numnod0 ) ] )

!       test for other connections to have the same degrees of freedom

        do conn1 = conn + 1, problem%numconnections

          if ( problem%connections(conn1)%physq1 == physq .and. &
               problem%connections(conn1)%physq2 == &
                         problem%connections(conn)%physq2 .and. &
               problem%connections(conn1)%layer1 == layer .and. &
               problem%connections(conn1)%layer2 == &
                         problem%connections(conn)%layer2  ) then

!           add nodes of connection

            nn = worknumnod(nodenr+1,conn1) - worknumnod(nodenr,conn1)
            nodes0(numnod0+1:numnod0+nn) = worknod(conn1)%a( &
                 [ ( worknumnod(nodenr,conn1) + i, i=1,nn ) ] )
            numnod0 = numnod0 + nn

            done(conn1) = .true.

          end if

        end do

        if ( dofr == 0 ) cycle ! no need for further counting

!       sort the nodes

        call sort ( nodes0(1:numnod0) )

!       find unique node numbers

        numnod = 0

        do node = 1, numnod0

          ndnr0 = nodes0 ( node )

          if ( node < numnod0 ) then
!           node present more than once?
            if ( nodes0(node+1) == ndnr0 ) cycle
          end if

          numnod = numnod + 1
          nodes(numnod) = ndnr0

        end do

        call count_matrix_entries

      end do

    end do

!   deallocate worknod

    do conn = 1, problem%numconnections
      deallocate ( worknod(conn)%a )
    end do

!   deallocate arrays

    deallocate ( nodes0, nodes )
    deallocate ( dofc, posc )

    deallocate ( utype, work )
    deallocate ( worknumnod )
    deallocate ( work1 )

  contains


!   connecting nodes of the "second side" for each node on the "first side"

    subroutine nodes_connected ( pass )

!     pass=1: only counting
!     pass=2: storing nodes
      integer, intent(in) :: pass

      integer, dimension(:,:), allocatable :: elgrpelem2
      integer :: point1, geometry1, nodeset1, elementset1, &
                 point2, geometry2, nodeset2, elementset2, &
                 typegeometry, discretization, eleme, nelem, &
                 elgrp1, elem1, typeofobject, object2, node1, &
                 typegeometry2, conn, elem, elgrp2, node, ip

      type(geometry_t) :: geometry, geometryc


!     determine nodes involved and upper bound of connecting nodes for each node

      do conn = 1, problem%numconnections

!       what type of connection?

        typegeometry = problem%connections(conn)%typegeometry
        discretization = problem%connections(conn)%discretization
        elementset1 = problem%connections(conn)%elementset1
        object = problem%connections(conn)%object

        if ( typegeometry == 1 ) then

!         connection using points

          point1 = problem%connections(conn)%geometry1
          point2 = problem%connections(conn)%geometry2

          if ( pass == 1 ) then

!           count (upper bound) of connecting nodes

            worknumnod(mesh%points(point1)+1,conn) = &
                           worknumnod(mesh%points(point1)+1,conn) + 1

          else if ( pass == 2 ) then

!           store connecting nodes

            ndnr = mesh%points(point1)
            worknod(conn)%a( worknumnod(ndnr,conn) + work1(ndnr,conn) + 1 ) = &
                                       mesh%points(point2)
            work1(ndnr,conn) = work1(ndnr,conn) + 1

          end if

        else if ( typegeometry == 5 ) then

!         connection using nodesets

          nodeset1 = problem%connections(conn)%geometry1
          nodeset2 = problem%connections(conn)%geometry2

          if ( pass == 1 ) then

!           count (upper bound) of connecting nodes

            worknumnod(mesh%nodesets(nodeset1)%a+1,conn) = &
                           worknumnod(mesh%nodesets(nodeset1)%a+1,conn) + 1

          else if ( pass == 2 ) then

!           store connecting nodes

            worknod(conn)%a( worknumnod(mesh%nodesets(nodeset1)%a,conn) &
                             + work1(mesh%nodesets(nodeset1)%a,conn) + 1 ) = &
                                       mesh%nodesets(nodeset2)%a
            work1(mesh%nodesets(nodeset1)%a,conn) = &
                               work1(mesh%nodesets(nodeset1)%a,conn) + 1

          end if

        else if ( any ( typegeometry == [ 2, 3, 4 ] ) ) then

!         connection using geometry

          geometry1 = problem%connections(conn)%geometry1

          if ( typegeometry == 2 ) then
!           connection using curves
            call copy_g ( mesh%curves(geometry1), geometry )
          else if ( typegeometry == 3 ) then
!           connection using surfaces
            call copy_g ( mesh%surfaces(geometry1), geometry )
          else if ( typegeometry == 4 ) then
!           connection using volumes
            call copy_g ( mesh%volumes(geometry1), geometry )
          end if

          geometry2 = problem%connections(conn)%geometry2
          elementset2 = problem%connections(conn)%elementset2

          if ( geometry2 > 0 ) then

!           connection using geometry on the second side

            typegeometry2 = problem%connections(conn)%typegeometry2

            if ( typegeometry2 == 2 ) then
!             connection using curves
              call copy_g ( mesh%curves(geometry2), geometryc )
            else if ( typegeometry2 == 3 ) then
!             connection using surfaces
              call copy_g ( mesh%surfaces(geometry2), geometryc )
            else if ( typegeometry2 == 4 ) then
!             connection using volumes
              call copy_g ( mesh%volumes(geometry2), geometryc )
            else if ( typegeometry2 == 5 ) then
!             connection using nodeset; fill fake geometryc
              geometryc%nnodes = size(mesh%nodesets(geometry2)%a)
              allocate(geometryc%nodes(geometryc%nnodes))
              geometryc%nodes = mesh%nodesets(geometry2)%a
            end if

            if ( discretization == 0 ) then

!             weak connection using elements in geometry1

              nnodes2 = geometryc%elnumnod

              do elem = 1, geometry%nelem

                if ( pass == 1 ) then

!                 count (upper bound) of connecting nodes

                  worknumnod(geometry%topology(:,elem,2)+1,conn) = &
                     worknumnod(geometry%topology(:,elem,2)+1,conn) + nnodes2

                else if ( pass == 2 ) then

!                 store connecting nodes

                  nnodes1 = geometry%elnumnod
                  do node = 1, nnodes1
                    ndnr = geometry%topology(node,elem,2)
                    worknod(conn)%a( [ ( worknumnod(ndnr,conn) + &
                          work1(ndnr,conn) + i, i=1,nnodes2 ) ] ) = &
                                       geometryc%topology(:,elem,2)
                  end do
                  work1(geometry%topology(:,elem,2),conn) = &
                         work1(geometry%topology(:,elem,2),conn) + nnodes2

                end if

              end do

            else if ( discretization == 1 ) then

!             collocation using nodes in geometry1

              if ( pass == 1 ) then

!               count (upper bound) of connecting nodes

                worknumnod(geometry%nodes+1,conn) = &
                               worknumnod(geometry%nodes+1,conn) + 1

              else if ( pass == 2 ) then

!               store connecting nodes

                worknod(conn)%a( worknumnod(geometry%nodes,conn) &
                                 + work1(geometry%nodes,conn) + 1 ) = &
                                           geometryc%nodes
                work1(geometry%nodes,conn) = &
                                   work1(geometry%nodes,conn) + 1

              end if

            end if

          else if ( elementset2 > 0 ) then

!           connection using elementset on the second side

            if ( discretization == 0 ) then

!             weak connection using elements in geometry1

              elem = 0

              do elgrp2 = 1, mesh%nelgrp
                do eleme = 1, mesh%elementsets(elementset2)%grpnumel(elgrp2)

                  elem2 = &
                        mesh%elementsets(elementset2)%elements(elgrp2)%a(eleme)

                  elem = elem + 1

                  nnodes2 = mesh%elnumnod(elgrp2)

                  if ( pass == 1 ) then

!                   count (upper bound) of connecting nodes

                    worknumnod(geometry%topology(:,elem,2)+1,conn) = &
                       worknumnod(geometry%topology(:,elem,2)+1,conn) + nnodes2

                  else if ( pass == 2 ) then

!                   store connecting nodes

                    nnodes1 = geometry%elnumnod
                    do node = 1, nnodes1
                      ndnr = geometry%topology(node,elem,2)
                      worknod(conn)%a( [ ( worknumnod(ndnr,conn) + &
                            work1(ndnr,conn) + i, i=1,nnodes2 ) ] ) = &
                                         mesh%topology(elgrp2)%a(:,elem2)
                    end do
                    work1(geometry%topology(:,elem,2),conn) = &
                           work1(geometry%topology(:,elem,2),conn) + nnodes2

                  end if

                end do

              end do

            else if ( discretization == 1 ) then

!             collocation using nodes in geometry1

              if ( pass == 1 ) then

!               count (upper bound) of connecting nodes

                worknumnod(geometry%nodes+1,conn) = &
                               worknumnod(geometry%nodes+1,conn) + 1

              else if ( pass == 2 ) then

!               store connecting nodes

                worknod(conn)%a( worknumnod(geometry%nodes,conn) &
                                 + work1(geometry%nodes,conn) + 1 ) = &
                                       mesh%elementsets(elementset2)%nodes
                work1(geometry%nodes,conn) = &
                                   work1(geometry%nodes,conn) + 1

              end if

            end if

          end if

        else if ( elementset1 > 0 ) then

!         connection using elementsets

          elementset2 = problem%connections(conn)%elementset2

          if ( discretization == 0 ) then

!           weak connection using elements in elementset1

!           fill help array elgrpelem2 for groups and elements in elementset2

            nelem = mesh%elementsets(elementset1)%nelem
            allocate ( elgrpelem2(nelem,2) )

            elem = 0
            do elgrp2 = 1, mesh%nelgrp
              do eleme = 1, mesh%elementsets(elementset2)%grpnumel(elgrp2)
                elem = elem + 1
                elem2 = mesh%elementsets(elementset2)%elements(elgrp2)%a(eleme)
                elgrpelem2(elem,:) = [ elgrp2, elem2 ]
              end do
            end do

!           determine connections

            elem = 0

            do elgrp1 = 1, mesh%nelgrp
              do eleme = 1, mesh%elementsets(elementset1)%grpnumel(elgrp1)

                elem1 = mesh%elementsets(elementset1)%elements(elgrp1)%a(eleme)

                elem = elem + 1

                elgrp2 = elgrpelem2(elem,1)
                elem2 = elgrpelem2(elem,2)

                nnodes2 = mesh%elnumnod(elgrp2)

                if ( pass == 1 ) then

!                 count (upper bound) of connecting nodes

                  worknumnod(mesh%topology(elgrp1)%a(:,elem1)+1,conn) = &
                     worknumnod(mesh%topology(elgrp1)%a(:,elem1)+1,conn) &
                      + nnodes2

                else if ( pass == 2 ) then

!                 store connecting nodes

                  nnodes1 = mesh%elnumnod(elgrp1)
                  do node = 1, nnodes1
                    ndnr = mesh%topology(elgrp1)%a(node,elem1)
                    worknod(conn)%a( [ ( worknumnod(ndnr,conn) + &
                          work1(ndnr,conn) + i, i=1,nnodes2 ) ] ) = &
                                       mesh%topology(elgrp2)%a(:,elem2)
                  end do
                  work1(mesh%topology(elgrp1)%a(:,elem1),conn) = &
                         work1(mesh%topology(elgrp1)%a(:,elem1),conn) + nnodes2

                end if

              end do

            end do

            deallocate ( elgrpelem2 )

          else if ( discretization == 1 ) then

!           collocation using nodes in elementset1

            if ( pass == 1 ) then

!             count (upper bound) of connecting nodes

              worknumnod(mesh%elementsets(elementset1)%nodes+1,conn) = &
                  worknumnod(mesh%elementsets(elementset1)%nodes+1,conn) + 1

            else if ( pass == 2 ) then

!             store connecting nodes

              worknod(conn)%a( &
                     worknumnod(mesh%elementsets(elementset1)%nodes,conn)&
                    + work1(mesh%elementsets(elementset1)%nodes,conn) + 1 ) = &
                                      mesh%elementsets(elementset2)%nodes
              work1(mesh%elementsets(elementset1)%nodes,conn) = &
                            work1(mesh%elementsets(elementset1)%nodes,conn) + 1

            end if

          end if

        else if ( object > 0 ) then

!         connection using object

          typeofobject = mesh%objects(object)%typeofobject
          object2 = problem%connections(conn)%object2
          geometry2 = problem%connections(conn)%geometry2
          elementset2 = problem%connections(conn)%elementset2

!         some preliminary work

          if ( geometry2 > 0 ) then

!           connect to a geometry, set geometryc

            typegeometry2 = problem%connections(conn)%typegeometry2

            if ( typegeometry2 == 2 ) then
!             connection using curves
              call copy_g ( mesh%curves(geometry2), geometryc )
            else if ( typegeometry2 == 3 ) then
!             connection using surfaces
              call copy_g ( mesh%surfaces(geometry2), geometryc )
            else if ( typegeometry2 == 4 ) then
!             connection using volumes
              call copy_g ( mesh%volumes(geometry2), geometryc )
            else if ( typegeometry2 == 5 ) then
!             connection using nodeset; fill fake geometryc
              geometryc%nnodes = size(mesh%nodesets(geometry2)%a)
              allocate(geometryc%nodes(geometryc%nnodes))
              geometryc%nodes = mesh%nodesets(geometry2)%a
            end if

          end if

          if ( discretization == 0 ) then

!           weak connection using elements in object

!           preliminary work

            if ( elementset2 > 0 ) then

!             fill help array elgrpelem2 for groups and elements in elementset2

              nelem = mesh%elementsets(elementset2)%nelem
              allocate ( elgrpelem2(nelem,2) )

              elem = 0
              do elgrp2 = 1, mesh%nelgrp
                do eleme = 1, mesh%elementsets(elementset2)%grpnumel(elgrp2)
                  elem = elem + 1
                  elem2 = &
                        mesh%elementsets(elementset2)%elements(elgrp2)%a(eleme)
                  elgrpelem2(elem,:) = [ elgrp2, elem2 ]
                end do
              end do

            end if

!           start loop

            do elem = 1, mesh%objects(object)%nelem

              do ip = 1, mesh%objects(object)%ninti

                elgrp1 = mesh%objects(object)%grpelm_int(ip,1,elem)
                elem1  = mesh%objects(object)%grpelm_int(ip,2,elem)

                if ( elgrp1 == 0 .and. pass == 1 ) then  ! no nodes connected
                  write(*,'(2(/a)/3(a,i0)/)') &
                    'Error in create_sysmatrix_structure_connection: ',&
                    ' reference coordinates in object are missing,', &
                    ' object = ', object, ' element = ', elem, &
                    ' integration point = ', ip
                  stop
                end if

                if ( typeofobject == 2 ) then

!                 two-sided connections

                  elgrp2 = mesh%objects(object)%grpelm2_int(ip,1,elem)
                  elem2  = mesh%objects(object)%grpelm2_int(ip,2,elem)

                  if ( elgrp2 == 0 .and. pass == 1 ) then  ! no nodes connected
                    write(*,'(2(/a)/3(a,i0)/)') &
                      'Error in create_sysmatrix_structure_connection: ', &
                      ' connecting reference coordinates on second side ', &
                      ' of object are missing. object = ', object, &
                      ' element = ', elem, ' integration point = ', ip
                    stop
                  end if

                  nnodes2 = mesh%elnumnod(elgrp2)

                else if ( object2 > 0 ) then

!                 connect to other object

                  elgrp2 = mesh%objects(object2)%grpelm_int(ip,1,elem)
                  elem2  = mesh%objects(object2)%grpelm_int(ip,2,elem)

                  if ( elgrp2 == 0 .and. pass == 1 ) then  ! no nodes connected
                    write(*,'(2(/a)/3(a,i0)/)') &
                      'Error in create_sysmatrix_structure_connection: ', &
                      ' connecting reference coordinates in second object',&
                      ' are missing. object = ', object2, &
                      ' element = ', elem, ' integration point = ', ip
                    stop
                  end if

                  nnodes2 = mesh%elnumnod(elgrp2)

                else if ( geometry2 > 0 ) then

!                 connect to a geometry

                  nnodes2 = geometryc%elnumnod

                else if ( elementset2 > 0 ) then

!                 connect to an elementset

                  elgrp2 = elgrpelem2(elem,1)
                  elem2 = elgrpelem2(elem,2)

                  nnodes2 = mesh%elnumnod(elgrp2)

                end if

                if ( pass == 1 ) then

!                 count (upper bound) of connecting nodes

                  worknumnod(mesh%topology(elgrp1)%a(:,elem1)+1,conn) = &
                       worknumnod(mesh%topology(elgrp1)%a(:,elem1)+1,conn) &
                             + nnodes2

                else if ( pass == 2 ) then

!                 store connecting nodes

                  nnodes1 = mesh%elnumnod(elgrp1)
                  if ( geometry2 > 0 ) then
                    do node1 = 1, nnodes1
                      ndnr = mesh%topology(elgrp1)%a(node1,elem1)
                      worknod(conn)%a( [ ( worknumnod(ndnr,conn) + &
                            work1(ndnr,conn) + i, i=1,nnodes2 ) ] ) = &
                                         geometryc%topology(:,elem,2)
                    end do
                  else ! object, object2 or elementset2
                    do node1 = 1, nnodes1
                      ndnr = mesh%topology(elgrp1)%a(node1,elem1)
                      worknod(conn)%a( [ ( worknumnod(ndnr,conn) + &
                            work1(ndnr,conn) + i, i=1,nnodes2 ) ] ) = &
                                         mesh%topology(elgrp2)%a(:,elem2)
                    end do
                  end if
                  work1(mesh%topology(elgrp1)%a(:,elem1),conn) = &
                         work1(mesh%topology(elgrp1)%a(:,elem1),conn) + nnodes2


                end if

              end do

            end do

            if ( elementset2 > 0 ) deallocate ( elgrpelem2 )

          else if ( discretization == 1 ) then

!           collocation using nodes in object

            nnodes = mesh%objects(object)%nnodes

            do node = 1, nnodes

              elgrp1 = mesh%objects(object)%grpelm(node,1)
              elem1  = mesh%objects(object)%grpelm(node,2)

              if ( elgrp1 == 0 .and. pass == 1 ) then  ! no nodes connected
                write(*,'(2(/a)/2(a,i0)/)') &
                  'Error in create_sysmatrix_structure_connection: ', &
                  ' reference coordinates in object are missing,', &
                  ' object = ', object, ' object node = ', node
                stop
              end if

              if ( mesh%objects(object)%typeofobject == 2 ) then

!               two-sided connections

                elgrp2 = mesh%objects(object)%grpelm2(node,1)
                elem2  = mesh%objects(object)%grpelm2(node,2)

                if ( elgrp2 == 0 .and. pass == 1 ) then   ! no nodes connected
                  write(*,'(2(/a)/2(a,i0)/)') &
                    'Error in create_sysmatrix_structure_connection: ', &
                    ' connecting reference coordinates on second side of',&
                    ' object are missing. object = ', object, &
                    ' object node = ', node
                  stop
                end if

                nnodes2 = mesh%elnumnod(elgrp2)

              else if ( problem%connections(conn)%object2 > 0 ) then

!               connect to other object

                object2 = problem%connections(conn)%object2
                elgrp2 = mesh%objects(object2)%grpelm(node,1)
                elem2  = mesh%objects(object2)%grpelm(node,2)

                if ( elgrp2 == 0 .and. pass == 1 ) then   ! no nodes connected
                  write(*,'(2(/a)/2(a,i0)/)') &
                    'Error in create_sysmatrix_structure_connection: ', &
                    ' connecting reference coordinates in second object', &
                    ' are missing. object = ', object2, &
                    ' object node = ', node
                  stop
                end if

                nnodes2 = mesh%elnumnod(elgrp2)

              else if ( geometry2 > 0 ) then

!               connect to a geometry

                nnodes2 = 1

              else if ( elementset2 > 0 ) then

!               connect to an elementset

                nnodes2 = 1

              end if

              if ( pass == 1 ) then

!               count (upper bound) of connecting nodes

                worknumnod(mesh%topology(elgrp1)%a(:,elem1)+1,conn) = &
                   worknumnod(mesh%topology(elgrp1)%a(:,elem1)+1,conn) + nnodes2

              else if ( pass == 2 ) then

!               store connecting nodes

                nnodes1 = mesh%elnumnod(elgrp1)
                if ( geometry2 > 0 ) then
                  do node1 = 1, nnodes1
                    ndnr = mesh%topology(elgrp1)%a(node1,elem1)
                    worknod(conn)%a( worknumnod(ndnr,conn) + &
                          work1(ndnr,conn) + 1 ) = geometryc%nodes(node)
                  end do
                else if ( elementset2 > 0 ) then
                  do node1 = 1, nnodes1
                    ndnr = mesh%topology(elgrp1)%a(node1,elem1)
                    worknod(conn)%a( worknumnod(ndnr,conn) + &
                          work1(ndnr,conn) + 1 ) = &
                                  mesh%elementsets(elementset2)%nodes(node)
                  end do
                else ! object or object2
                  do node1 = 1, nnodes1
                    ndnr = mesh%topology(elgrp1)%a(node1,elem1)
                    worknod(conn)%a( [ ( worknumnod(ndnr,conn) + &
                          work1(ndnr,conn) + i, i=1,nnodes2 ) ] ) = &
                                       mesh%topology(elgrp2)%a(:,elem2)
                  end do
                end if
                work1(mesh%topology(elgrp1)%a(:,elem1),conn) = &
                       work1(mesh%topology(elgrp1)%a(:,elem1),conn) + nnodes2

              end if

            end do

          end if

        end if

      end do

    end subroutine nodes_connected


!   count matrix entries of node connections

    subroutine count_matrix_entries

      integer :: dr, dc, i

      do dr = 1, dofr

        deg_r = posr(dr)

!       fill work with degrees already accounted for in standard FEM

        do i = 1, numnod1 + 1
          work( worksc1(i):worksc2(i) ) = 1
        end do

        row = problem%degfdperm(deg_r,2)

!       count matrix entries of node connections

        do i = 1, numnod

          node = nodes(i)

          physq = problem%connections(conn)%physq2
          layer = problem%connections(conn)%layer2

          if ( physq > 0 ) then
            call pos_array_local_node ( problem, node, dof, pos, [physq], &
              layer )
          else
            call pos_array_local_node ( problem, node, dof, pos, layer=layer )
          end if

!         column degrees: store for later use

          dofc(i) = dof
          posc(i,1:dof) = problem%nodnumdegfd(node) + pos(1:dof)

!         start loop over degrees

          do dc = 1, dofc(i)

            deg_c = posc(i,dc)

            if ( work(deg_c) /= 0 ) cycle ! degree already counted

            col = problem%degfdperm(deg_c,2)

            if ( utype(deg_r) == 1 ) then

!             unknown degree of freedom (row)

              if ( utype(deg_c) == 1 ) then

!               unknown degree of freedom (column)

                if ( sysmatrix%symmetric ) then
!                 symmetric matrix
                  if ( col > row ) then
!                   element in upper triangle
                    sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1
                  else
!                   element in lower triangle, add transposed element
                    sysmatrix%Suu%ia(col+1) = sysmatrix%Suu%ia(col+1) + 1
                  end if
                else
!                 unsymmetric matrix
                  sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1
                  sysmatrix%Suu%ia(col+1) = sysmatrix%Suu%ia(col+1) + 1
                end if

              else if ( utype(deg_c) == 3 ) then

!               prescribed degree of freedom (column)
                sysmatrix%Sup%ia(row+1) = sysmatrix%Sup%ia(row+1) + 1
                colp = col - numund
                sysmatrix%Spu%ia(colp+1) = sysmatrix%Spu%ia(colp+1) + 1

              end if

            else if ( utype(deg_r) == 3 ) then

!             prescribed degree of freedom (row)

              rowp = row - numund

              if ( utype(deg_c) == 1 ) then

!               unknown degree of freedom (column)
                sysmatrix%Sup%ia(col+1) = sysmatrix%Sup%ia(col+1) + 1
                sysmatrix%Spu%ia(rowp+1) = sysmatrix%Spu%ia(rowp+1) + 1

              else if ( utype(deg_c) == 3 ) then

!               prescribed degree of freedom (column)
                sysmatrix%Spp%ia(rowp+1) = sysmatrix%Spp%ia(rowp+1) + 1
                colp = col - numund
                sysmatrix%Spp%ia(colp+1) = sysmatrix%Spp%ia(colp+1) + 1

              end if

            end if

            work(deg_c) = 1 ! degree counted

          end do

        end do

!       set work (degrees found) back to zero

        do i = 1, numnod
          work( posc(i,1:dofc(i)) ) = 0
        end do

        do i = 1, numnod1 + 1
          work( worksc1(i):worksc2(i) ) = 0
        end do

      end do

    end subroutine count_matrix_entries

  end subroutine create_sysmatrix_structure_connection


! Assemble system matrix and system vector for connections

  subroutine build_system_connection ( mesh, problem, sysmatrix, sysvector, &
    msysvector, m2sysvector, m3sysvector, elemsub, elemsub1, elemsub2, &
    coefficients, oldvectors, order, connection1, connection2, connections, &
    buildmatrix, buildvector, addmatvec, addmat, addvec, factormat, &
    factorvec, transform )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the system to be assembled
    type(sysmatrix_t), intent(inout) :: sysmatrix

!   the right-hand side vector to be assembled
    type(sysvector_t), intent(inout), optional :: sysvector

!   multiple right-hand side vectors to be assembled
    type(sysvector_t), dimension(:), intent(inout), optional :: msysvector
    type(sysvector_t), dimension(:,:), intent(inout), optional :: m2sysvector
    type(sysvector_t), dimension(:,:,:), intent(inout), optional :: m3sysvector

!   this is the element subroutine that must be supplied by the calling routine
!   Note that the element matrices and vectors are adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
    optional elemsub
    interface
      subroutine elemsub ( mesh, problem, conn, elem, node, matrix, vector, &
        first, last, coefficients, oldvectors, elemmat11, elemmat12, &
        elemmat21, elemmat22, elemvec1, elemvec2 )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: conn, elem, node
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat11, elemmat12, &
          elemmat21, elemmat22
        real(dp), intent(out), dimension(:) :: elemvec1, elemvec2
      end subroutine elemsub
    end interface

!   this is the element subroutine that must be supplied by the calling routine
!   Note that the element matrices and vectors are adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
!   Compared to elemsub the interface is extended with elgrp1, elem1, elgrp2,
!   and elem2, which is useful for elementsets.
    optional elemsub1
    interface
      subroutine elemsub1 ( mesh, problem, conn, elem, node, elgrp1, elem1, &
        elgrp2, elem2, matrix, vector, first, last, coefficients, oldvectors, &
        elemmat11, elemmat12, elemmat21, elemmat22, elemvec1, elemvec2 )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: conn, elem, node, elgrp1, elem1, elgrp2, elem2
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat11, elemmat12, &
          elemmat21, elemmat22
        real(dp), intent(out), dimension(:) :: elemvec1, elemvec2
      end subroutine elemsub1
    end interface

!   this is the element subroutine that must be supplied by the calling routine
!   Note that the element matrices and vectors are adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
!   Compared to elemsub the interface is extended with nodenr1, nodenr2,
!   which can be useful for collocated connections.
    optional elemsub2
    interface
      subroutine elemsub2 ( mesh, problem, conn, elem, node, nodenr1, nodenr2, &
        matrix, vector, first, last, coefficients, oldvectors, elemmat11, &
        elemmat12, elemmat21, elemmat22, elemvec1, elemvec2 )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: conn, elem, node, nodenr1, nodenr2
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat11, elemmat12, &
          elemmat21, elemmat22
        real(dp), intent(out), dimension(:) :: elemvec1, elemvec2
      end subroutine elemsub2
    end interface

!   if coefficients is present, it is passed through to the element subroutine
    type(coefficients_t), intent(in), optional :: coefficients

!   if oldvectors is present, it is passed through to the element subroutine.
!   This structure can be used to supply old sysvectors/vectors from possibly
!   other problems to the element subroutines. Note that the problem structures
!   of the other problems must be included in oldvectors.
    type(oldvectors_t), intent(in), optional :: oldvectors

!   The parameter order determines the sequence of the degrees of freedom on
!   elementlevel. There are only two possibilities and affect only the two loops
!   over nodal points and degrees of freedom (within a physical quantity if
!   problem%nphysq>0):
!     order = 'ND' : the most inner loop is over the degrees of freedom
!     order = 'DN' : the most inner loop is over the nodal points
!   Specifically we have:   nphysq = 0   nphysq > 0
!         order = 'ND' :       ND           PND
!         order = 'DN' :       DN           PDN
!   NOTE: if nphys >0 the main ordering (outside loop) of the element
!   degrees of freedom are the physical quantities.
!   If problem%numlayers > 0 the degrees are stacked into layers and we have
!                           nphysq = 0   nphysq > 0
!         order = 'ND' :      NLD           PNLD
!         order = 'DN' :      LDN           PLDN
!   See also the userguide for further explanation.
!
!   For example consider an element with two nodes and two physical quanties:
!   one with two degrees (u,v) and one with a single degree p. Then we have:
!     order = 'ND' : u1, v1, u2, v2, p1, p2
!     order = 'DN' : u1, u2, v1, v2, p1, p2
!   The parameter order generates a permutation of the degrees of freedom on
!   element level, making the life of an element programmer easier,
!   but _does not affect_ the global numbering of the unknowns.
!   The (column) layout of the assembled sparse matrix can be different due to
!   the different sequence of unknowns in the assembling process.
!
!   The default of order is:
!     order = 'ND' : if no physical quantities have been defined (nphysq=0)
!     order = 'DN' : if physical quantities have been defined (nphysq>0)
    character(len=*), intent(in), optional :: order

!   if these are present the assembling takes place for connections
!   connection1,...,connection2 only. If only connection1 is present one
!   connections is assembled only.
    integer, intent(in), optional :: connection1, connection2

!   if present: the connections to be assembled
!   For example connections=(/2,4/) will assemble for connections 2 and 4.
!   Default: all connections.
    integer, dimension(:), intent(in), optional :: connections

!   if one of these is assigned the value .false. the assembling of
!   the matrix or vector is not done.
    logical, intent(in), optional :: buildmatrix, buildvector

!   if addmatvec is set to .true. the matrix and vector are not cleared before
!   the assembling and thus the element matrices and vectors are added to an
!   existing system. The default is .false. (clearing)
    logical, intent(in), optional :: addmatvec

!   if addmat is set to .true. the matrix is not cleared before
!   the assembling and thus the element matrices and vectors are added to an
!   existing system. The default is .false. (clearing)
!   NOTE: if addmatvec is present, the value of addmat will be ignored.
    logical, intent(in), optional :: addmat

!   if addvec is set to .true. the vector is not cleared before
!   the assembling and thus the element matrices and vectors are added to an
!   existing system. The default is .false. (clearing)
!   NOTE: if addmatvec is present, the value of addvec will be ignored.
    logical, intent(in), optional :: addvec

!   if present the element matrix and/or vector will be multiplied by the
!   specified factor before assembly into the system matrix/vector.
    real(dp), intent(in), optional :: factormat, factorvec

!   By setting transform=.false. the transformation matrix is not applied
!   to the element matrix and vector. This means, the element degrees of
!   freedom are defined in the local (transformed) system.
!   default = .true.
    logical, intent(in), optional :: transform


    logical :: matrix, vector, laddmat, laddvec, first, last
    logical, dimension(mesh%maxelnumnod) :: lp
    logical :: posgrpelem1, posgrpelem2, posgeomelem1, posgeomelem2, &
               posnode1, posnode2, ltransform
!   work1: array for denoting transformed degrees of freedom
    logical, allocatable, dimension(:) :: work1
    integer :: elem, nelem, ninti, ip
    integer :: numund, numess, rowg
    integer :: dof1, dof2, nrhsd, rhsd
    integer :: cnn, object, elem1, elgrp1, elem2, elgrp2, conn, object2, eleme
    integer :: physq1, physq2, ndof1, ndof2, nnodes, nodenr1, nodenr2
    integer :: lconns(problem%numconnections), lnconns, i
    integer, allocatable, dimension(:) :: pos1, pos2, w1, w2
    integer :: node, j, si, sj, layer1, layer2, k
    integer :: point1, geometry1, nodeset1, elementset1, &
               point2, geometry2, nodeset2, elementset2, &
               typegeometry, discretization, typeofobject, typegeometry2
    integer, dimension(:,:), allocatable :: elgrpelem2
    real(dp), allocatable, dimension(:,:) :: elemmat11, elemmat12, elemmat21, &
      elemmat22
    real(dp), allocatable, dimension(:) :: elemvec1, elemvec2
    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl
    type(geometry_t) :: geometry, geometryc


    call check ( mesh, 'build_system_connection' )
    call check ( problem, 'build_system_connection', mesh )


!   initialize local parameters

    nrhsd = 1

    if ( present(oldvectors) ) oldvl = oldvectors
    if ( present(coefficients) ) coeffl = coefficients

!   ordering; set defaults

    if ( present(order) ) then
      if ( all ( order /= [ 'ND', 'DN' ] ) ) then
        write(*,'(2(/a)/)') &
          'Error in build_system_connection: ', &
          ' heading parameter order must be either ''ND'', ''DN''.'
        stop
      end if
    end if

!   test elemsub

    if ( count( [ present(elemsub), present(elemsub1), &
                  present(elemsub2) ] ) /= 1 ) then
      write(*,'(3(/a)/)') &
        'Error in build_system_connection: ', &
        ' one of elemsub, elemsub1, elemsub2 ', &
        ' must be present in the heading '
      stop
    end if

!   which connections?

    if ( present(connection1) .and. present(connection2) ) then
!     specified range of connections only
      lnconns = connection2 - connection1 + 1
      lconns(1:lnconns)= [ (i,i=connection1,connection2) ]
    else if ( present(connection1) ) then
!     one connection only
      lnconns = 1
      lconns(1) = connection1
    else if ( present(connections) ) then
      lnconns = size(connections)
      lconns(1:lnconns) = connections
    else
!     all connections
      lnconns = problem%numconnections
      lconns = [ (i,i=1,lnconns) ]
    end if

    if ( any ( lconns(1:lnconns) < 1 ) .or. &
         any ( lconns(1:lnconns) > problem%numconnections ) ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error: connection1, connection2, or connections in the heading of ', &
        ' build_system_connection is out of range: ', &
        ' some are < 1 or larger than the number of connections ', &
        problem%numconnections
      stop
    end if

!   test limitation in implementation of layers

    if ( problem%numlayers > 0 ) then
      if ( .not. all( problem%connections(lconns(1:lnconns))%layer1 > 0 ) .or. &
           .not. all( problem%connections(lconns(1:lnconns))%layer2 > 0 ) ) then
        write(*,'(4(/a)/)') &
          'Error in build_system_connection: ', &
          ' Layers have been defined and layer1=0 or layer2=0.', &
          ' Layers not implemented for this routine without specifying ', &
          ' a layer on both sides of the connection.'
        stop
      end if
    end if

!   compute matrix and vector?

    if ( present(buildmatrix) ) then
      matrix = buildmatrix
    else
      matrix = .true.
    end if

    if ( present(buildvector) ) then
      vector = buildvector
    else
      vector = .true.
    end if

    if ( .not. ( matrix .or. vector ) ) return

!   clear matrix and vector?

    if ( matrix ) then

      if ( .not. sysmatrix%initialized_structure ) then
        write(*,'(/a/)') &
          'Error in build_system_connection: no system matrix structure.'
        stop
      end if

      if ( .not. sysmatrix%finalized ) then
        write(*,'(/2a/)') &
          'Error in build_system_connection: ', &
          'system matrix has not been finalized.'
        stop
      end if

      if ( .not. sysmatrix%allocated_data ) then
        write(*,'(/2a/)') &
          'Error in build_system_connection: ', &
          'data in system matrix not allocated.'
        stop
      end if

      if ( sysmatrix%Suu%m + sysmatrix%Sup%m /= problem%numdegfd ) then
        write(*,'(/a/a/a/)') &
          'Error in build_system_connection: ', &
          'number of degrees of freedom of the ', &
          'system matrix is different from the number in problem.'
        stop
      end if

      if ( present(addmatvec) ) then
        laddmat = addmatvec
      else if ( present(addmat) ) then
        laddmat = addmat
      else
        laddmat = .false.
      end if

      if ( .not. laddmat ) then

!       clear matrix
        call clear_sysmatrix ( sysmatrix )

      end if

    end if

    if ( vector ) then

      if ( present(addmatvec) ) then
        laddvec = addmatvec
      else if ( present(addvec) ) then
        laddvec = addvec
      else
        laddvec = .false.
      end if

      if ( present(sysvector) ) then

!       single right-hand side

        if ( .not. sysvector%created ) then
          write(*,'(/a/)') &
            'Error in build_system_connection: sysvector not created.'
          stop
        end if

        if ( sysvector%n /= problem%numdegfd ) then
          write(*,'(/a/a/a/)') &
            'Error in build_system_connection: ', &
            'number of degrees of freedom of the ', &
            'sysvector vector different from the number in problem.'
          stop
        end if

        if ( .not. laddvec ) then

!         clear vector
          sysvector%u = 0

        end if

      else if ( present(msysvector) ) then

!       multiple right-hand sides

        nrhsd = size(msysvector)

        if ( .not. all(msysvector%created) ) then
          write(*,'(/a/)') &
            'Error in build_system_connection: msysvector not created.'
          stop
        end if

        if ( any( msysvector%n /= problem%numdegfd ) ) then
          write(*,'(/a/a/a/)') &
            'Error in build_system_connection: ', &
            'number of degrees of freedom of the ', &
            'msysvector vector different from the number in problem.'
          stop
        end if

        if ( .not. laddvec ) then

!         clear vector
          do rhsd = 1, nrhsd
            msysvector(rhsd)%u = 0
          end do

        end if

      else if ( present(m2sysvector) ) then

!       multiple right-hand sides (2d array)

        nrhsd = size(m2sysvector)

        if ( .not. all(m2sysvector%created) ) then
          write(*,'(/a/)') &
            'Error in build_system_connection: m2sysvector not created.'
          stop
        end if

        if ( any( m2sysvector%n /= problem%numdegfd ) ) then
          write(*,'(/a/a/a/)') &
            'Error in build_system_connection: ', &
            'number of degrees of freedom of the ', &
            'm2sysvector vector different from the number in problem.'
          stop
        end if

        if ( .not. laddvec ) then

!         clear vectors
          do i = 1, size(m2sysvector,1)
            do j = 1, size(m2sysvector,2)
              m2sysvector(i,j)%u = 0
            end do
          end do

        end if

      else if ( present(m3sysvector) ) then

!       multiple right-hand sides (2d array)

        nrhsd = size(m3sysvector)

        if ( .not. all(m3sysvector%created) ) then
          write(*,'(/a/)') &
            'Error in build_system_connection: m3sysvector not created.'
          stop
        end if

        if ( any( m3sysvector%n /= problem%numdegfd ) ) then
          write(*,'(/a/a/a/)') &
            'Error in build_system_connection: ', &
            'number of degrees of freedom of the ', &
            'm3sysvector vector different from the number in problem.'
          stop
        end if

        if ( .not. laddvec ) then

!         clear vectors
          do i = 1, size(m3sysvector,1)
            do j = 1, size(m3sysvector,2)
              do k = 1, size(m3sysvector,3)
                m3sysvector(i,j,k)%u = 0
              end do
            end do
          end do

        end if

      else

        write(*,'(/a/a/)') &
          'Error in build_system_connection:', &
          ' no sysvector, msysvector, m2sysvector or m3sysvector present'
        stop

      end if

    end if

!   transformations

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in build_system_connection: ', &
        ' transformation matrix has not been build'
      stop
    end if

    ltransform = set_optional ( variable=transform, default=.true. )

    if ( problem%numtransdegfd > 0 .and. ltransform ) then

!     work array to indicate transformed degrees of freedom
      allocate ( work1(problem%numdegfd) )
      work1 = .false.
      work1(problem%degfdtrans) = .true.

    end if

!   temporary work arrays to store information on current large matrix row

    numess = problem%numessdegfd
    numund = problem%numundegfd

    allocate( w1(numund), w2(numess) )

    w1 = 0
    w2 = 0

!   start loop over all connections

    do cnn = 1, lnconns

      conn = lconns(cnn)

!     set some connection properties

      typegeometry = problem%connections(conn)%typegeometry
      discretization = problem%connections(conn)%discretization
      elementset1 = problem%connections(conn)%elementset1
      object  = problem%connections(conn)%object

      physq1 = problem%connections(conn)%physq1
      physq2 = problem%connections(conn)%physq2
      layer1 = problem%connections(conn)%layer1
      layer2 = problem%connections(conn)%layer2

!     logicals for which posarrays must be computed before assembling

      posgrpelem1 = .false.
      posgrpelem2 = .false.
      posgeomelem1 = .false.
      posgeomelem2 = .false.
      posnode1 = .false.
      posnode2 = .false.

!     what type of connection?

      if ( typegeometry == 1 ) then

!       connection using points

        posnode1 = .true.
        posnode2 = .true.

        point1 = problem%connections(conn)%geometry1
        point2 = problem%connections(conn)%geometry2

!       find number of degrees of freedom in this point

        nodenr1 = mesh%points(point1)
        nodenr2 = mesh%points(point2)

        ndof1 = ndof_nodenr ( nodenr1, physq1 )
        ndof2 = ndof_nodenr ( nodenr2, physq2 )

        if ( ndof1 == 0 .or. ndof2 == 0 ) cycle ! nothing to do here

        call allocate_elemmatvec

        first = .true.; last = .true.

        call compute_assemble_elemmatvec

        call deallocate_elemmatvec

      else if ( typegeometry == 5 ) then

!       connection using nodesets

        posnode1 = .true.
        posnode2 = .true.

        nodeset1 = problem%connections(conn)%geometry1
        nodeset2 = problem%connections(conn)%geometry2

        do node = 1, size(mesh%nodesets(nodeset1)%a)

!         find number of degrees of freedom in this node

          nodenr1 = mesh%nodesets(nodeset1)%a(node)
          nodenr2 = mesh%nodesets(nodeset2)%a(node)

          ndof1 = ndof_nodenr ( nodenr1, physq1 )
          ndof2 = ndof_nodenr ( nodenr2, physq2 )

          if ( ndof1 == 0 .or. ndof2 == 0 ) cycle ! nothing to do here

          call allocate_elemmatvec

          first = node == 1
          last  = node == size(mesh%nodesets(nodeset1)%a)

          call compute_assemble_elemmatvec

          call deallocate_elemmatvec

        end do

      else if ( any ( typegeometry == [ 2, 3, 4 ] ) ) then

!       connection using geometry

        geometry1 = problem%connections(conn)%geometry1
        geometry2 = problem%connections(conn)%geometry2
        elementset2 = problem%connections(conn)%elementset2

        if ( typegeometry == 2 ) then
!         connection using curves
          call copy_g ( mesh%curves(geometry1), geometry )
        else if ( typegeometry == 3 ) then
!         connection using surfaces
          call copy_g ( mesh%surfaces(geometry1), geometry )
        else if ( typegeometry == 4 ) then
!         connection using volumes
          call copy_g ( mesh%volumes(geometry1), geometry )
        end if

        if ( geometry2 > 0 ) then

!         connection using geometry on the second side

          typegeometry2 = problem%connections(conn)%typegeometry2

          if ( typegeometry2 == 2 ) then
!           connection using curves
            call copy_g ( mesh%curves(geometry2), geometryc )
          else if ( typegeometry2 == 3 ) then
!             connection using surfaces
            call copy_g ( mesh%surfaces(geometry2), geometryc )
          else if ( typegeometry2 == 4 ) then
!           connection using volumes
            call copy_g ( mesh%volumes(geometry2), geometryc )
          else if ( typegeometry2 == 5 ) then
!           connection using nodeset; fill fake geometryc
            geometryc%nnodes = size(mesh%nodesets(geometry2)%a)
            allocate(geometryc%nodes(geometryc%nnodes))
            geometryc%nodes = mesh%nodesets(geometry2)%a
          end if

          if ( discretization == 0 ) then

!           weak

            posgeomelem1 = .true.
            posgeomelem2 = .true.

            do elem = 1, geometry%nelem

              ndof1 = ndof_geomelem ( geometry, elem, physq1 )
              ndof2 = ndof_geomelem ( geometryc, elem, physq2 )

              if ( ndof1 == 0 .or. ndof2 == 0 ) cycle ! nothing to do here

              call allocate_elemmatvec

              first = elem == 1
              last  = elem == geometry%nelem

              call compute_assemble_elemmatvec

              call deallocate_elemmatvec

            end do

          else if ( discretization == 1 ) then

!           collocation

            posnode1 = .true.
            posnode2 = .true.

            do node = 1, geometry%nnodes

              nodenr1 = geometry%nodes(node)
              nodenr2 = geometryc%nodes(node)

              ndof1 = ndof_nodenr ( nodenr1, physq1 )
              ndof2 = ndof_nodenr ( nodenr2, physq2 )

              if ( ndof1 == 0 .or. ndof2 == 0 ) cycle ! nothing to do here

              call allocate_elemmatvec

              first = node == 1
              last  = node == geometry%nnodes

              call compute_assemble_elemmatvec

              call deallocate_elemmatvec

            end do

          end if

        else if ( elementset2 > 0 ) then

!         connection using elementset on the second side

          if ( discretization == 0 ) then

            posgeomelem1 = .true.
            posgrpelem2 = .true.

!           weak connection using elements in geometry1

            elem = 0

            do elgrp2 = 1, mesh%nelgrp
              do eleme = 1, mesh%elementsets(elementset2)%grpnumel(elgrp2)

                elem2 = mesh%elementsets(elementset2)%elements(elgrp2)%a(eleme)

                elem = elem + 1

                ndof1 = ndof_geomelem ( geometry, elem, physq1 )
                ndof2 = ndof_elgrpelem ( elgrp2, elem2, physq2 )

                if ( ndof1 == 0 .or. ndof2 == 0 ) cycle ! nothing to do here

                call allocate_elemmatvec

                first = elem == 1
                last  = elem == geometry%nelem

                call compute_assemble_elemmatvec

                call deallocate_elemmatvec

              end do
            end do

          else if ( discretization == 1 ) then

!           collocation using nodes in geometry1

            posnode1 = .true.
            posnode2 = .true.

            do node = 1, geometry%nnodes

              nodenr1 = geometry%nodes(node)
              nodenr2 = mesh%elementsets(elementset2)%nodes(node)

              ndof1 = ndof_nodenr ( nodenr1, physq1 )
              ndof2 = ndof_nodenr ( nodenr2, physq2 )

              if ( ndof1 == 0 .or. ndof2 == 0 ) cycle ! nothing to do here

              call allocate_elemmatvec

              first = node == 1
              last  = node == geometry%nnodes

              call compute_assemble_elemmatvec

              call deallocate_elemmatvec

            end do

          end if

        end if

      else if ( elementset1 > 0 ) then

!       connection using elementsets

        elementset2 = problem%connections(conn)%elementset2

        if ( discretization == 0 ) then

!         weak connection using elements in elementset1

          posgrpelem1 = .true.
          posgrpelem2 = .true.

!         fill help array elgrpelem2 for groups and elements in elementset2

          nelem = mesh%elementsets(elementset1)%nelem
          allocate ( elgrpelem2(nelem,2) )

          elem = 0
          do elgrp2 = 1, mesh%nelgrp
            do eleme = 1, mesh%elementsets(elementset2)%grpnumel(elgrp2)
              elem = elem + 1
              elem2 = mesh%elementsets(elementset2)%elements(elgrp2)%a(eleme)
              elgrpelem2(elem,:) = [ elgrp2, elem2 ]
            end do
          end do

          elem = 0

          do elgrp1 = 1, mesh%nelgrp
            do eleme = 1, mesh%elementsets(elementset1)%grpnumel(elgrp1)

              elem1 = mesh%elementsets(elementset1)%elements(elgrp1)%a(eleme)

              elem = elem + 1

              elgrp2 = elgrpelem2(elem,1)
              elem2 = elgrpelem2(elem,2)

              ndof1 = ndof_elgrpelem ( elgrp1, elem1, physq1 )
              ndof2 = ndof_elgrpelem ( elgrp2, elem2, physq2 )

              if ( ndof1 == 0 .or. ndof2 == 0 ) cycle ! nothing to do here

              call allocate_elemmatvec

              first = elem == 1
              last  = elem == nelem

              call compute_assemble_elemmatvec

              call deallocate_elemmatvec

            end do
          end do

          deallocate ( elgrpelem2 )

        else if ( discretization == 1 ) then

!         collocation using nodes in elementset1

          posnode1 = .true.
          posnode2 = .true.

          nnodes = mesh%elementsets(elementset1)%nnodes

          do node = 1, nnodes

            nodenr1 = mesh%elementsets(elementset1)%nodes(node)
            nodenr2 = mesh%elementsets(elementset2)%nodes(node)

            ndof1 = ndof_nodenr ( nodenr1, physq1 )
            ndof2 = ndof_nodenr ( nodenr2, physq2 )

            if ( ndof1 == 0 .or. ndof2 == 0 ) cycle ! nothing to do here

            call allocate_elemmatvec

            first = node == 1
            last  = node == nnodes

            call compute_assemble_elemmatvec

            call deallocate_elemmatvec

          end do

        end if

      else if ( object > 0 ) then

!       connection using object

        typeofobject = mesh%objects(object)%typeofobject
        object2 = problem%connections(conn)%object2
        geometry2 = problem%connections(conn)%geometry2
        elementset2 = problem%connections(conn)%elementset2

        posgrpelem1 = .true.

!       some preliminary work

        if ( geometry2 > 0 ) then

!         connect to a geometry, set geometryc

          typegeometry2 = problem%connections(conn)%typegeometry2

          if ( typegeometry2 == 2 ) then
!           connection using curves
            call copy_g ( mesh%curves(geometry2), geometryc )
          else if ( typegeometry2 == 3 ) then
!           connection using surfaces
            call copy_g ( mesh%surfaces(geometry2), geometryc )
          else if ( typegeometry2 == 4 ) then
!           connection using volumes
            call copy_g ( mesh%volumes(geometry2), geometryc )
          else if ( typegeometry2 == 5 ) then
!           connection using nodeset; fill fake geometryc
            geometryc%nnodes = size(mesh%nodesets(geometry2)%a)
            allocate(geometryc%nodes(geometryc%nnodes))
            geometryc%nodes = mesh%nodesets(geometry2)%a
          end if

        end if

        if ( discretization == 0 ) then

!         weak connection

!         preliminary work

          if ( elementset2 > 0 ) then

!           fill help array elgrpelem2 for groups and elements in elementset2

            nelem = mesh%elementsets(elementset2)%nelem
            allocate ( elgrpelem2(nelem,2) )

            elem = 0
            do elgrp2 = 1, mesh%nelgrp
              do eleme = 1, mesh%elementsets(elementset2)%grpnumel(elgrp2)
                elem = elem + 1
                elem2 = &
                      mesh%elementsets(elementset2)%elements(elgrp2)%a(eleme)
                elgrpelem2(elem,:) = [ elgrp2, elem2 ]
              end do
            end do

          end if

          if ( typeofobject == 2 .or. object2 > 0 .or. elementset2 > 0 ) then
            posgrpelem2 = .true.
          else if ( geometry2 > 0 ) then
            posgeomelem2 = .true.
          end if

          nelem = mesh%objects(object)%nelem
          ninti = mesh%objects(object)%ninti

!         loop over all integration points on object

          do elem = 1, nelem
            do ip = 1, ninti

              node = ip

              elgrp1 = mesh%objects(object)%grpelm_int(ip,1,elem)
              elem1  = mesh%objects(object)%grpelm_int(ip,2,elem)

              if ( typeofobject == 2 ) then

!               two-sided connections

                elgrp2 = mesh%objects(object)%grpelm2_int(ip,1,elem)
                elem2  = mesh%objects(object)%grpelm2_int(ip,2,elem)

              else if ( object2 > 0 ) then

!               connect to other object

                elgrp2 = mesh%objects(object2)%grpelm_int(ip,1,elem)
                elem2  = mesh%objects(object2)%grpelm_int(ip,2,elem)

              else if ( elementset2 > 0 ) then

!               connect to an elementset

                elgrp2 = elgrpelem2(elem,1)
                elem2 = elgrpelem2(elem,2)

              end if

              ndof1 = ndof_elgrpelem ( elgrp1, elem1, physq1 )
              if ( posgrpelem2 ) then
                ndof2 = ndof_elgrpelem ( elgrp2, elem2, physq2 )
              else if ( posgeomelem2 ) then
                ndof2 = ndof_geomelem ( geometryc, elem, physq2 )
              end if

              if ( ndof1 == 0 .or. ndof2 == 0 ) cycle ! nothing to do here

              call allocate_elemmatvec

              first = ip == 1 .and. elem == 1
              last  = ip == ninti .and. elem == nelem

              call compute_assemble_elemmatvec

              call deallocate_elemmatvec

            end do
          end do

          if ( elementset2 > 0 ) deallocate ( elgrpelem2 )

        else if ( discretization == 1 ) then

!         collocation

          if ( typeofobject == 2 .or. object2 > 0 ) then
            posgrpelem2 = .true.
          else if ( geometry2 > 0 .or. elementset2 > 0 ) then
            posnode2 = .true.
          end if

          do node = 1, mesh%objects(object)%nnodes

            elgrp1 = mesh%objects(object)%grpelm(node,1)
            elem1  = mesh%objects(object)%grpelm(node,2)

            if ( typeofobject == 2 ) then

!             two-sided connections

              elgrp2 = mesh%objects(object)%grpelm2(node,1)
              elem2  = mesh%objects(object)%grpelm2(node,2)

            else if ( object2 > 0 ) then

!             connect to other object

              elgrp2 = mesh%objects(object2)%grpelm(node,1)
              elem2  = mesh%objects(object2)%grpelm(node,2)

            else if ( geometry2 > 0 ) then

!             connect to a geometry

              nodenr2 = geometryc%nodes(node)

            else if ( elementset2 > 0 ) then

!             connect to an elementset

              nodenr2 = mesh%elementsets(elementset2)%nodes(node)

            end if

            ndof1 = ndof_elgrpelem ( elgrp1, elem1, physq1 )
            if ( posgrpelem2 ) then
              ndof2 = ndof_elgrpelem ( elgrp2, elem2, physq2 )
            else if ( posnode2 ) then
              ndof2 = ndof_nodenr ( nodenr2, physq2 )
            end if

            if ( ndof1 == 0 .or. ndof2 == 0 ) cycle ! nothing to do here

            call allocate_elemmatvec

            first = node == 1
            last  = node == mesh%objects(object)%nnodes

            call compute_assemble_elemmatvec

            call deallocate_elemmatvec

          end do

        end if

      end if

    end do

    deallocate( w1, w2 )

    if ( problem%numtransdegfd > 0. .and. ltransform ) then
      deallocate ( work1 )
    end if

  contains


!   number of degrees in a node

    function ndof_nodenr ( nodenr, physq ) result(ndof)

      integer, intent(in) :: nodenr, physq
      integer :: ndof

      if ( physq > 0 ) then
        ndof =  problem%vec_nodnumdegfd(nodenr+1,problem%physq(physq)) &
                   - problem%vec_nodnumdegfd(nodenr,problem%physq(physq))
      else
        ndof = problem%nodnumdegfd(nodenr+1) - problem%nodnumdegfd(nodenr)
      end if

    end function ndof_nodenr


!   number of degrees in an element of a geometry

    function ndof_geomelem ( geometry, elem, physq ) result(ndof)

      type(geometry_t), intent(in) :: geometry
      integer, intent(in) :: elem, physq
      integer :: ndof

      if ( physq > 0 ) then
        ndof = sum ( problem%vec_nodnumdegfd(geometry%topology(:,elem,2)+1,&
                                        &problem%physq(physq)) &
                         - problem%vec_nodnumdegfd(geometry%topology(:,elem,2),&
                                        &problem%physq(physq)) )
      else
        ndof = sum ( problem%nodnumdegfd(geometry%topology(:,elem,2)+1) &
                   - problem%nodnumdegfd(geometry%topology(:,elem,2)) )
      end if

    end function ndof_geomelem


!   number of degrees in an element of a group

    function ndof_elgrpelem ( elgrp, elem, physq ) result(ndof)

      integer, intent(in) :: elgrp, elem, physq
      integer :: ndof

      if ( physq > 0 ) then
        ndof = sum ( problem%vec_elnumdegfd(elgrp)%a(:,physq) )
      else
        ndof = sum ( problem%elnumdegfd(elgrp)%a )
      end if

    end function ndof_elgrpelem


!   allocate the element matrix and vector

    subroutine allocate_elemmatvec

!     reserve memory for element matrix and vector and positions

      allocate ( elemmat11(ndof1,ndof1), elemmat12(ndof1,ndof2) )
      allocate ( elemmat21(ndof2,ndof1), elemmat22(ndof2,ndof2) )
      allocate ( elemvec1(ndof1*nrhsd), elemvec2(ndof2*nrhsd) )
      allocate ( pos1(ndof1), pos2(ndof2) )

    end subroutine allocate_elemmatvec


!   deallocate the element matrix and vector

    subroutine deallocate_elemmatvec

!     remove memory for element matrix and vector and positions

      deallocate ( elemmat11, elemmat12, elemmat21, elemmat22 )
      deallocate ( elemvec1, elemvec2, pos1, pos2 )

    end subroutine deallocate_elemmatvec


!   compute and assemble the element matrix and vector


    subroutine compute_assemble_elemmatvec

      integer :: k, rhsd, row, i, j
      type(sparsematrix_t) :: A

!     compute element matrix and vector

      if ( present(elemsub) ) then

!       simple interface only

        call elemsub ( mesh, problem, conn, elem, node, matrix, vector, &
          first, last, coeffl, oldvl, elemmat11, elemmat12, &
          elemmat21, elemmat22, elemvec1, elemvec2 )

      else if ( present(elemsub1) ) then

!       interface includes elgrp1/elgrp2, elem1/elem2
!       useful for weak connection of elementsets.

        call elemsub1 ( mesh, problem, conn, elem, node, elgrp1, elem1, elgrp2,&
          elem2, matrix, vector, first, last, coeffl, oldvl, elemmat11, &
          elemmat12, elemmat21, elemmat22, elemvec1, elemvec2 )

      else if ( present(elemsub2) ) then

!       interface includes nodenr1, nodenr2
!       useful for collocated connections.

        call elemsub2 ( mesh, problem, conn, elem, node, nodenr1, nodenr2, &
          matrix, vector, first, last, coeffl, oldvl, elemmat11, elemmat12, &
          elemmat21, elemmat22, elemvec1, elemvec2 )

      end if

!     compute positions in large matrix/vector of element degrees of freedom

      if ( posgrpelem1 ) then

!       normal element in a group

        if ( physq1 > 0 ) then
!         physical quantity specified
          call pos_array ( mesh, problem, elgrp1, elem1, dof1, pos1, &
            physqarr=[physq1], order=order, layer=layer1, lp=lp )
        else
          call pos_array ( mesh, problem, elgrp1, elem1, dof1, pos1, &
            order=order, layer=layer1, lp=lp )
        end if

        if ( layer1 > 0 ) then
          if ( .not. all( lp(1:mesh%elnumnod(elgrp1)) ) ) then
            write(*,'(3(/a)/4(a,i0)/)') &
                'Error in build_system_connection: ', &
              ' First connecting element not fully within layer.', &
              ' Not implemented.', &
              ' connection = ', conn, ' elgrp1 = ', elgrp1, &
              ' elem1 = ', elem1, ' layer1 = ', layer1
            stop
          end if
        end if

      else if ( posgeomelem1 ) then

!       element in a geometry

        if ( physq1 > 0 ) then
!         physical quantity specified
          call pos_array_geometry ( problem, geometry, elem, dof1, pos1, &
            physqarr=[physq1], order=order, layer=layer1, lp=lp )
        else
          call pos_array_geometry ( problem, geometry, elem, dof1, pos1, &
            order=order, layer=layer1, lp=lp )
        end if

        if ( layer1 > 0 ) then
          if ( .not. all( lp(1:geometry%elnumnod) ) ) then
            write(*,'(3(/a)/4(a,i0)/)') &
                'Error in build_system_connection: ', &
              ' First connecting element not fully within layer.', &
              ' Not implemented.', &
              ' connection = ', conn, ' geometry = ', geometry1, &
              ' elem1 = ', elem, ' layer1 = ', layer1
            stop
          end if
        end if

      else if ( posnode1 ) then

!       node connect

        if ( physq1 > 0 ) then
!         physical quantity specified
          call pos_array_node ( problem, nodenr1, dof1, pos1, &
            physqarr=[physq1], layer=layer1 )
        else
          call pos_array_node ( problem, nodenr1, dof1, pos1, layer=layer1 )
        end if

      end if

      if ( posgrpelem2 ) then

!       normal element in a group

        if ( physq2 > 0 ) then
!         physical quantity specified
          call pos_array ( mesh, problem, elgrp2, elem2, dof2, pos2, &
            physqarr=[physq2], order=order, layer=layer2, lp=lp )
        else
          call pos_array ( mesh, problem, elgrp2, elem2, dof2, pos2, &
            order=order, layer=layer2, lp=lp )
        end if

        if ( layer2 > 0 ) then
          if ( .not. all( lp(1:mesh%elnumnod(elgrp2)) ) ) then
            write(*,'(3(/a)/3(a,i0)/)') &
              'Error in build_system_connection: ', &
              ' Second connecting element not fully within layer.', &
              ' Not implemented.', &
              ' elgrp2 = ', elgrp2, ' elem2 = ', elem2, ' layer2 = ', layer2
            stop
          end if
        end if

      else if ( posgeomelem2 ) then

!       element in a geometry

        if ( physq2 > 0 ) then
!         physical quantity specified
          call pos_array_geometry ( problem, geometryc, elem, dof2, pos2, &
            physqarr=[physq2], order=order, layer=layer2, lp=lp )
        else
          call pos_array_geometry ( problem, geometryc, elem, dof2, pos2, &
            order=order, layer=layer2, lp=lp )
        end if

        if ( layer2 > 0 ) then
          if ( .not. all( lp(1:geometry%elnumnod) ) ) then
            write(*,'(3(/a)/4(a,i0)/)') &
                'Error in build_system_connection: ', &
              ' First connecting element not fully within layer.', &
              ' Not implemented.', &
              ' connection = ', conn, ' geometry = ', geometry2, &
              ' elem2 = ', elem, ' layer2 = ', layer2
            stop
          end if
        end if

      else if ( posnode2 ) then

!       node connect

        if ( physq2 > 0 ) then
!         physical quantity specified
          call pos_array_node ( problem, nodenr2, dof2, pos2, &
            physqarr=[physq2], layer=layer2 )
        else
          call pos_array_node ( problem, nodenr2, dof2, pos2, layer=layer2 )
        end if

      end if

      if ( dof1 /= ndof1 .or. dof2 /= ndof2 ) then
        write(*,'(/a/a/)') &
          'Internal error in build_system_connection: ', &
          ' dof1 /= ndof1 .or. dof2 /= ndof2'
        stop
      end if

!     transformations

      if ( problem%numtransdegfd > 0 .and. ltransform ) then

         if ( any ( work1(pos1) ) ) then

!          extract transformation matrix

           A = extract_submat ( problem%Amat, pos1, pos1 )

!          transform matrices

           if ( matrix ) then
             elemmat11 = fsmatmul ( elemmat11, A )
             elemmat11 = stfmatmul ( A, elemmat11 )
             elemmat12 = stfmatmul ( A, elemmat12 )
             elemmat21 = fsmatmul ( elemmat21, A )
           end if

!          transform vector

           if ( vector ) then

!            possibly multiple right-hand side

             do rhsd = 1, nrhsd
               k = dof1*(rhsd-1)
               elemvec1(k+1:k+dof1) = stmatvec ( A, elemvec1(k+1:k+dof1) )
             end do

           end if

           call delete(A)

         end if

         if ( any ( work1(pos2) ) ) then

!          extract transformation matrix

           A = extract_submat ( problem%Amat, pos2, pos2 )

!          transform matrices

           if ( matrix ) then
             elemmat12 = fsmatmul ( elemmat12, A )
             elemmat21 = stfmatmul ( A, elemmat21 )
             elemmat22 = fsmatmul ( elemmat22, A )
             elemmat22 = stfmatmul ( A, elemmat22 )
           end if

!          transform vector

           if ( vector ) then

!            possibly multiple right-hand side

             do rhsd = 1, nrhsd
               k = dof2*(rhsd-1)
               elemvec2(k+1:k+dof2) = stmatvec ( A, elemvec2(k+1:k+dof2) )
             end do

           end if

           call delete(A)

         end if

      end if

!     add matrix

      if ( matrix ) then

!       add element matrix row by row

        if ( present(factormat) ) then
          elemmat11 = factormat * elemmat11
          elemmat12 = factormat * elemmat12
          elemmat21 = factormat * elemmat21
          elemmat22 = factormat * elemmat22
        end if

        call add_elemmat_to_sysmatrix ( sysmatrix, elemmat11, pos1, pos1, &
          w1, w2 )
        call add_elemmat_to_sysmatrix ( sysmatrix, elemmat12, pos1, pos2, &
          w1, w2 )
        call add_elemmat_to_sysmatrix ( sysmatrix, elemmat21, pos2, pos1, &
          w1, w2 )
        call add_elemmat_to_sysmatrix ( sysmatrix, elemmat22, pos2, pos2, &
          w1, w2 )

      end if

!     add vector

      if ( vector ) then

!       add element vector to large vector

        if ( present(factorvec) ) then
          elemvec1 = factorvec * elemvec1
          elemvec2 = factorvec * elemvec2
        end if

        if ( present(msysvector) ) then

!         multiple right-hand side

          do row = 1, dof1

            rowg = pos1(row) ! global row number

            do rhsd = 1, nrhsd
              msysvector(rhsd)%u(rowg) = msysvector(rhsd)%u(rowg) &
                                          + elemvec1( row + dof1*(rhsd-1) )
            end do

          end do

          do row = 1, dof2

            rowg = pos2(row) ! global row number

            do rhsd = 1, nrhsd
              msysvector(rhsd)%u(rowg) = msysvector(rhsd)%u(rowg) &
                                          + elemvec2( row + dof2*(rhsd-1) )
            end do

          end do

        else if ( present(m2sysvector) ) then

!         multiple right-hand side (2d array)

          si = size(m2sysvector,1)

          do row = 1, dof1

            rowg = pos1(row) ! global row number

            do j = 1, size(m2sysvector,2)
              do i = 1, si
                rhsd = i + si*(j-1)
                m2sysvector(i,j)%u(rowg) = m2sysvector(i,j)%u(rowg) &
                                          + elemvec1( row + dof1*(rhsd-1) )
              end do
            end do

          end do

          do row = 1, dof2

            rowg = pos2(row) ! global row number

            do j = 1, size(m2sysvector,2)
              do i = 1, si
                rhsd = i + si*(j-1)
                m2sysvector(i,j)%u(rowg) = m2sysvector(i,j)%u(rowg) &
                                          + elemvec2( row + dof2*(rhsd-1) )
              end do
            end do

          end do

        else if ( present(m3sysvector) ) then

!         multiple right-hand side (3d array)

          si = size(m3sysvector,1)
          sj = size(m3sysvector,2)

          do row = 1, dof1

            rowg = pos1(row) ! global row number

            do k = 1, size(m3sysvector,3)
              do j = 1, sj
                do i = 1, si
                  rhsd = i + si*(j-1) + si*sj*(k-1)
                  m3sysvector(i,j,k)%u(rowg) = m3sysvector(i,j,k)%u(rowg) &
                                               + elemvec1( row + dof1*(rhsd-1) )
                end do
              end do
            end do

          end do

          do row = 1, dof2

            rowg = pos2(row) ! global row number

            do k = 1, size(m3sysvector,3)
              do j = 1, sj
                do i = 1, si
                  rhsd = i + si*(j-1) + si*sj*(k-1)
                  m3sysvector(i,j,k)%u(rowg) = m3sysvector(i,j,k)%u(rowg) &
                                               + elemvec2( row + dof2*(rhsd-1) )
                end do
              end do
            end do

          end do

        else

!         single right-hand side

          do row = 1, dof1

            rowg = pos1(row) ! global row number

            sysvector%u(rowg) = sysvector%u(rowg) + elemvec1(row)

          end do

          do row = 1, dof2

            rowg = pos2(row) ! global row number

            sysvector%u(rowg) = sysvector%u(rowg) + elemvec2(row)

          end do

        end if

      end if

    end subroutine compute_assemble_elemmatvec

  end subroutine build_system_connection

end module system_connection_m
