
! Copyright (C) 2004-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines routines for input/output of mesh, problem, ...

module tfem_utils_m

  use kind_defs_m
  use array_defs_m
  use mesh_m
  use limits_m, only: MAXPOINTS, MAXCURVES, MAXOBJECTS, MAXSURFACES, &
                      MAXBLOCKS, MAXVOLUMES, MAXNODESETS, MAXELEMENTSETS
  use problem_defs_m, only: input_probdef_t
  use problem_m, only: create_input_probdef
  use element_defs_m

  implicit none

  integer :: unit_tfemu = 13 ! unit number for reading/writing

contains


! write mesh to a file

  subroutine write_mesh ( mesh, filename )

    type(mesh_t), intent(in) :: mesh

!   the filename for writing the mesh to
    character (len=*), intent(in) :: filename

!   This routine writes basic mesh information (points, curves, topology
!   and coordinates) to an intermediate file for later reading by read_mesh
!   from a different main program. The file format is therefore binary and
!   not intended for a permanent mesh storage file since it might change in
!   the future.


    integer :: elgrp, curve, object, surface, block, volume, nodeset, sblock, &
      nodblock, elementset, ios, m, nblend


    if ( .not. mesh%meshgen ) then
      write(*,'(/a/)') 'Error in write_mesh: mesh contains no data'
      stop
    end if

!   open file

    open ( unit=unit_tfemu, file=filename, form='unformatted', iostat=ios, &
      status='replace' )

    if ( ios /= 0 ) then
      write(*,'(/a,i0/)') 'Error in write_mesh: ios = ', ios
      stop
    end if

!   some scalar data

    write ( unit=unit_tfemu ) mesh%ndim, mesh%nelem, mesh%nnodes, mesh%nelgrp, &
      mesh%npoints, mesh%ncurves, mesh%nsurfaces, mesh%nvolumes, &
      mesh%nobjects, mesh%nblocks, mesh%nnodesets, mesh%renumber, &
      mesh%sblocks%filled, mesh%nodblocks%filled, mesh%nelementsets, &
      mesh%nblend

    if ( mesh%multlvlref ) then
      write(*,'(/a/a/)') 'Warning in write_mesh:', &
        ' writing of multi-level elements not yet been implemented.'
    end if

!   element

    write ( unit=unit_tfemu ) mesh%element(:)%elshape
    write ( unit=unit_tfemu ) mesh%element(:)%globalshape
    write ( unit=unit_tfemu ) mesh%element(:)%numnod
    write ( unit=unit_tfemu ) mesh%element(:)%ndim
    do elgrp = 1, mesh%nelgrp
      write ( unit=unit_tfemu ) mesh%element(elgrp)%p
    end do

    do m = 1, mesh%nblend
      write ( unit=unit_tfemu ) mesh%element_blend(:,m)%elshape
      write ( unit=unit_tfemu ) mesh%element_blend(:,m)%globalshape
      write ( unit=unit_tfemu ) mesh%element_blend(:,m)%numnod
      write ( unit=unit_tfemu ) mesh%element_blend(:,m)%ndim
      do elgrp = 1, mesh%nelgrp
        write ( unit=unit_tfemu ) mesh%element_blend(elgrp,m)%p
      end do
    end do

    write ( unit=unit_tfemu ) mesh%nnodes_blend

!   topology

    write ( unit=unit_tfemu ) mesh%grpnumel

    do elgrp = 1, mesh%nelgrp
      write ( unit=unit_tfemu ) size(mesh%topology(elgrp)%a,1)
      write ( unit=unit_tfemu ) mesh%topology(elgrp)%a
    end do

!   points

    write ( unit=unit_tfemu ) mesh%points(:mesh%npoints)

!   curves

    write ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%ndim
    write ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%nnodes
    write ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%nelem
    write ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%nblend
    write ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%n
    write ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%m

    write ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%element%globalshape
    write ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%element%numnod
    write ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%element%elshape
    write ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%element%ndim
    do curve = 1, mesh%ncurves
      write ( unit=unit_tfemu ) mesh%curves(curve)%element%p
      write ( unit=unit_tfemu ) mesh%curves(curve)%nnodes_blend
      write ( unit=unit_tfemu ) mesh%curves(curve)%nodes
      write ( unit=unit_tfemu ) size(mesh%curves(curve)%topology,1)
      write ( unit=unit_tfemu ) mesh%curves(curve)%topology
    end do

    do curve = 1, mesh%ncurves
      nblend = mesh%curves(curve)%nblend
      write ( unit=unit_tfemu ) &
                     mesh%curves(curve)%element_blend(:nblend)%globalshape
      write ( unit=unit_tfemu ) mesh%curves(curve)%element_blend(:nblend)%numnod
      write ( unit=unit_tfemu ) &
                     mesh%curves(curve)%element_blend(:nblend)%elshape
      write ( unit=unit_tfemu ) mesh%curves(curve)%element_blend(:nblend)%ndim
      do m = 1, nblend
        write ( unit=unit_tfemu ) mesh%curves(curve)%element_blend(m)%p
      end do
    end do

!   surfaces

    write ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%ndim
    write ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%nnodes
    write ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%nelem
    write ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%nblend
    write ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%n
    write ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%m
    write ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%element%globalshape
    write ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%element%numnod
    write ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%element%elshape
    write ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%element%ndim
    do surface = 1, mesh%nsurfaces
      write ( unit=unit_tfemu ) mesh%surfaces(surface)%element%p
      write ( unit=unit_tfemu ) mesh%surfaces(surface)%nnodes_blend
      write ( unit=unit_tfemu ) mesh%surfaces(surface)%nodes
      write ( unit=unit_tfemu ) size(mesh%surfaces(surface)%topology,1)
      write ( unit=unit_tfemu ) mesh%surfaces(surface)%topology
    end do

    do surface = 1, mesh%nsurfaces
      nblend = mesh%surfaces(surface)%nblend
      write ( unit=unit_tfemu ) &
                     mesh%surfaces(surface)%element_blend(:nblend)%globalshape
      write ( unit=unit_tfemu ) &
                     mesh%surfaces(surface)%element_blend(:nblend)%numnod
      write ( unit=unit_tfemu ) &
                     mesh%surfaces(surface)%element_blend(:nblend)%elshape
      write ( unit=unit_tfemu ) &
                     mesh%surfaces(surface)%element_blend(:nblend)%ndim
      do m = 1, nblend
        write ( unit=unit_tfemu ) mesh%surfaces(surface)%element_blend(m)%p
      end do
    end do

!   volumes

    write ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%ndim
    write ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%nnodes
    write ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%nelem
    write ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%nblend
    write ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%n
    write ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%m
    write ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%element%globalshape
    write ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%element%numnod
    write ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%element%elshape
    write ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%element%ndim
    do volume = 1, mesh%nvolumes
      write ( unit=unit_tfemu ) mesh%volumes(volume)%element%p
      write ( unit=unit_tfemu ) mesh%volumes(volume)%nnodes_blend
      write ( unit=unit_tfemu ) mesh%volumes(volume)%nodes
      write ( unit=unit_tfemu ) size(mesh%volumes(volume)%topology,1)
      write ( unit=unit_tfemu ) mesh%volumes(volume)%topology
    end do

    do volume = 1, mesh%nvolumes
      nblend = mesh%volumes(volume)%nblend
      write ( unit=unit_tfemu ) &
                     mesh%volumes(volume)%element_blend(:nblend)%globalshape
      write ( unit=unit_tfemu ) &
                     mesh%volumes(volume)%element_blend(:nblend)%numnod
      write ( unit=unit_tfemu ) &
                     mesh%volumes(volume)%element_blend(:nblend)%elshape
      write ( unit=unit_tfemu ) &
                     mesh%volumes(volume)%element_blend(:nblend)%ndim
      do m = 1, nblend
        write ( unit=unit_tfemu ) mesh%volumes(volume)%element_blend(m)%p
      end do
    end do

!   blocks

    write ( unit=unit_tfemu ) mesh%blocks(:mesh%nblocks)%nelem
    do block = 1, mesh%nblocks
      write ( unit=unit_tfemu ) mesh%blocks(block)%bounds
      write ( unit=unit_tfemu ) mesh%blocks(block)%elbounds
      write ( unit=unit_tfemu ) mesh%blocks(block)%elements
    end do

!   objects

    write ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%topol
    write ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%intpoints
    write ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%typeofobject
    write ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%ndim
    write ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%nnodes
    write ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%nodeset
    write ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%nodeset2
    do object = 1, mesh%nobjects
      write ( unit=unit_tfemu ) mesh%objects(object)%groups
      write ( unit=unit_tfemu ) mesh%objects(object)%coor
      if ( mesh%objects(object)%typeofobject == 2 ) then
        write ( unit=unit_tfemu ) mesh%objects(object)%groups2
      end if
      if ( mesh%objects(object)%topol ) then
        write ( unit=unit_tfemu ) mesh%objects(object)%nelem
        write ( unit=unit_tfemu ) mesh%objects(object)%elnumnod
        write ( unit=unit_tfemu ) mesh%objects(object)%element%globalshape
        write ( unit=unit_tfemu ) mesh%objects(object)%element%numnod
        write ( unit=unit_tfemu ) mesh%objects(object)%element%elshape
        write ( unit=unit_tfemu ) mesh%objects(object)%element%p
        write ( unit=unit_tfemu ) mesh%objects(object)%element%ndim
        write ( unit=unit_tfemu ) mesh%objects(object)%topology
      end if
      if ( mesh%objects(object)%intpoints ) then
        write ( unit=unit_tfemu ) mesh%objects(object)%ninti
        write ( unit=unit_tfemu ) size(mesh%objects(object)%xig,2)
        write ( unit=unit_tfemu ) mesh%objects(object)%xig
        write ( unit=unit_tfemu ) mesh%objects(object)%wg
      end if
    end do

!   nodesets

    do nodeset = 1, mesh%nnodesets
      write ( unit=unit_tfemu ) size(mesh%nodesets(nodeset)%a)
      write ( unit=unit_tfemu ) mesh%nodesets(nodeset)%a
    end do

!   elementsets

    do elementset = 1, mesh%nelementsets
      write ( unit=unit_tfemu ) mesh%elementsets(elementset)%nodes_created
      write ( unit=unit_tfemu ) mesh%elementsets(elementset)%nelem
      write ( unit=unit_tfemu ) mesh%elementsets(elementset)%nnodes
      write ( unit=unit_tfemu ) mesh%elementsets(elementset)%nodes
      write ( unit=unit_tfemu ) mesh%elementsets(elementset)%grpnumel
      do elgrp = 1, mesh%nelgrp
        write ( unit=unit_tfemu ) &
                size(mesh%elementsets(elementset)%elements(elgrp)%a)
        write ( unit=unit_tfemu ) mesh%elementsets(elementset)%elements(elgrp)%a
      end do
    end do

!   sblocks

    if ( mesh%sblocks%filled ) then

      write ( unit=unit_tfemu ) mesh%sblocks%nsblocks

      write ( unit=unit_tfemu ) mesh%sblocks%nbl
      write ( unit=unit_tfemu ) mesh%sblocks%x0
      write ( unit=unit_tfemu ) mesh%sblocks%dx

      write ( unit=unit_tfemu ) mesh%sblocks%sblks(:)%nelem

      do sblock = 1, mesh%sblocks%nsblocks
        write ( unit=unit_tfemu ) mesh%sblocks%sblks(sblock)%elbounds
        write ( unit=unit_tfemu ) mesh%sblocks%sblks(sblock)%elements
      end do

    end if

!   nodblocks

    if ( mesh%nodblocks%filled ) then

      write ( unit=unit_tfemu ) mesh%nodblocks%nnodblocks

      write ( unit=unit_tfemu ) mesh%nodblocks%nbl
      write ( unit=unit_tfemu ) mesh%nodblocks%x0
      write ( unit=unit_tfemu ) mesh%nodblocks%dx

      write ( unit=unit_tfemu ) mesh%nodblocks%nodblks(:)%nnodes

      do nodblock = 1, mesh%nodblocks%nnodblocks
        write ( unit=unit_tfemu ) mesh%nodblocks%nodblks(nodblock)%nodes
      end do

    end if

!   coordinates

    write ( unit=unit_tfemu ) mesh%coor

!   mesh renumbering

    if ( mesh%renumber ) then
      write ( unit=unit_tfemu ) mesh%nodperm
    end if

    close ( unit=unit_tfemu )

  end subroutine write_mesh


! read mesh from a file written by write_mesh

  subroutine read_mesh ( mesh, filename )

    type(mesh_t), intent(inout) :: mesh

!   the filename for reading the mesh from
    character (len=*), intent(in) :: filename

!   This routine reads basic mesh information (points, curves, surfaces,
!   topology and coordinates) from an intermediate file written by write_mesh.
!   The file format is binary and not intended for a permanent mesh storage
!   file since it might change in the future.


    integer :: elgrp, curve, object, ios, surface, block, nelem, ndimxi, &
      volume, nodeset, nnodes, sblock, nodblock, elementset, m, nblend, &
      elnumnod


    if ( mesh%meshgen ) then
      write(*,'(/a/)') 'Error in read_mesh: mesh already contains data'
      stop
    end if

!   open file

    open ( unit=unit_tfemu, file=filename, form='unformatted', iostat=ios, &
      status='old' )

    if ( ios /= 0 ) then
      write(*,'(/2a/)') 'Error in read_mesh: cannot open file ', filename
      stop
    end if

!   some scalar data

    read ( unit=unit_tfemu ) mesh%ndim, mesh%nelem, mesh%nnodes, mesh%nelgrp, &
      mesh%npoints, mesh%ncurves, mesh%nsurfaces, mesh%nvolumes, &
      mesh%nobjects, mesh%nblocks, mesh%nnodesets, mesh%renumber, &
      mesh%sblocks%filled, mesh%nodblocks%filled, mesh%nelementsets, &
      mesh%nblend

    mesh%multlvlref = .false. ! reading of multi-level el. not yet implemented

!   element

    allocate( mesh%element(mesh%nelgrp) )

    read ( unit=unit_tfemu ) mesh%element(:)%elshape
    read ( unit=unit_tfemu ) mesh%element(:)%globalshape
    read ( unit=unit_tfemu ) mesh%element(:)%numnod
    read ( unit=unit_tfemu ) mesh%element(:)%ndim
    do elgrp = 1, mesh%nelgrp
      read ( unit=unit_tfemu ) mesh%element(elgrp)%p
    end do

    allocate( mesh%element_blend(mesh%nelgrp,mesh%nblend) )

    do m = 1, mesh%nblend
      read ( unit=unit_tfemu ) mesh%element_blend(:,m)%elshape
      read ( unit=unit_tfemu ) mesh%element_blend(:,m)%globalshape
      read ( unit=unit_tfemu ) mesh%element_blend(:,m)%numnod
      read ( unit=unit_tfemu ) mesh%element_blend(:,m)%ndim
      do elgrp = 1, mesh%nelgrp
        read ( unit=unit_tfemu ) mesh%element_blend(elgrp,m)%p
      end do
    end do

    allocate ( mesh%nnodes_blend(mesh%nblend+2) )

    read ( unit=unit_tfemu ) mesh%nnodes_blend

!   topology

    allocate( mesh%grpnumel(mesh%nelgrp) )

    read ( unit=unit_tfemu ) mesh%grpnumel

    allocate( mesh%topology(mesh%nelgrp) )

    do elgrp = 1, mesh%nelgrp
      read ( unit=unit_tfemu ) elnumnod
      allocate( mesh%topology(elgrp)%a(elnumnod,mesh%grpnumel(elgrp) ) )
      read ( unit=unit_tfemu ) mesh%topology(elgrp)%a
    end do

!   points

    allocate( mesh%points(max(mesh%npoints,MAXPOINTS)) )

    read ( unit=unit_tfemu ) mesh%points(:mesh%npoints)

!   curves

    allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )

    read ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%ndim
    read ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%nnodes
    read ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%nelem
    read ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%nblend
    read ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%n
    read ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%m
    read ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%element%globalshape
    read ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%element%numnod
    read ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%element%elshape
    read ( unit=unit_tfemu ) mesh%curves(:mesh%ncurves)%element%ndim
    do curve = 1, mesh%ncurves
      read ( unit=unit_tfemu ) mesh%curves(curve)%element%p
      allocate( mesh%curves(curve)%nnodes_blend(mesh%curves(curve)%nblend+2) )
      read ( unit=unit_tfemu ) mesh%curves(curve)%nnodes_blend
      allocate( mesh%curves(curve)%nodes(mesh%curves(curve)%nnodes) )
      read ( unit=unit_tfemu ) mesh%curves(curve)%nodes
      read ( unit=unit_tfemu ) elnumnod
      allocate( mesh%curves(curve)%topology(elnumnod,&
                                                  mesh%curves(curve)%nelem,2))
      read ( unit=unit_tfemu ) mesh%curves(curve)%topology
    end do

    do curve = 1, mesh%ncurves
      nblend = mesh%curves(curve)%nblend
      allocate ( mesh%curves(curve)%element_blend(nblend) )
      read ( unit=unit_tfemu ) &
                     mesh%curves(curve)%element_blend(:nblend)%globalshape
      read ( unit=unit_tfemu ) mesh%curves(curve)%element_blend(:nblend)%numnod
      read ( unit=unit_tfemu ) &
                     mesh%curves(curve)%element_blend(:nblend)%elshape
      read ( unit=unit_tfemu ) mesh%curves(curve)%element_blend(:nblend)%ndim
      do m = 1, nblend
        read ( unit=unit_tfemu ) mesh%curves(curve)%element_blend(m)%p
      end do
    end do

!   surfaces

    allocate( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )

    read ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%ndim
    read ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%nnodes
    read ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%nelem
    read ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%nblend
    read ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%n
    read ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%m
    read ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%element%globalshape
    read ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%element%numnod
    read ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%element%elshape
    read ( unit=unit_tfemu ) mesh%surfaces(:mesh%nsurfaces)%element%ndim
    do surface = 1, mesh%nsurfaces
      read ( unit=unit_tfemu ) mesh%surfaces(surface)%element%p
      allocate( mesh%surfaces(surface)%nnodes_blend(&
                                  mesh%surfaces(surface)%nblend+2) )
      read ( unit=unit_tfemu ) mesh%surfaces(surface)%nnodes_blend
      allocate( mesh%surfaces(surface)%nodes(mesh%surfaces(surface)%nnodes) )
      read ( unit=unit_tfemu ) mesh%surfaces(surface)%nodes
      read ( unit=unit_tfemu ) elnumnod
      allocate( mesh%surfaces(surface)%topology(elnumnod,&
                                            mesh%surfaces(surface)%nelem,2) )
      read ( unit=unit_tfemu ) mesh%surfaces(surface)%topology
    end do

    do surface = 1, mesh%nsurfaces
      nblend = mesh%surfaces(surface)%nblend
      allocate ( mesh%surfaces(surface)%element_blend(nblend) )
      read ( unit=unit_tfemu ) &
                     mesh%surfaces(surface)%element_blend(:nblend)%globalshape
      read ( unit=unit_tfemu ) &
                     mesh%surfaces(surface)%element_blend(:nblend)%numnod
      read ( unit=unit_tfemu ) &
                     mesh%surfaces(surface)%element_blend(:nblend)%elshape
      read ( unit=unit_tfemu ) &
                     mesh%surfaces(surface)%element_blend(:nblend)%ndim
      do m = 1, nblend
        read ( unit=unit_tfemu ) mesh%surfaces(surface)%element_blend(m)%p
      end do
    end do

!   volumes

    allocate( mesh%volumes(max(mesh%nvolumes,MAXVOLUMES)) )

    read ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%ndim
    read ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%nnodes
    read ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%nelem
    read ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%nblend
    read ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%n
    read ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%m
    read ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%element%globalshape
    read ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%element%numnod
    read ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%element%elshape
    read ( unit=unit_tfemu ) mesh%volumes(:mesh%nvolumes)%element%ndim
    do volume = 1, mesh%nvolumes
      read ( unit=unit_tfemu ) mesh%volumes(volume)%element%p
      allocate( mesh%volumes(volume)%nnodes_blend(&
                                  mesh%volumes(volume)%nblend+2) )
      read ( unit=unit_tfemu ) mesh%volumes(volume)%nnodes_blend
      allocate( mesh%volumes(volume)%nodes(mesh%volumes(volume)%nnodes) )
      read ( unit=unit_tfemu ) mesh%volumes(volume)%nodes
      read ( unit=unit_tfemu ) elnumnod
      allocate( mesh%volumes(volume)%topology(elnumnod,&
                                              mesh%volumes(volume)%nelem,2) )
      read ( unit=unit_tfemu ) mesh%volumes(volume)%topology
    end do

    do volume = 1, mesh%nvolumes
      nblend = mesh%volumes(volume)%nblend
      allocate ( mesh%volumes(volume)%element_blend(nblend) )
      read ( unit=unit_tfemu ) &
                     mesh%volumes(volume)%element_blend(:nblend)%globalshape
      read ( unit=unit_tfemu ) &
                     mesh%volumes(volume)%element_blend(:nblend)%numnod
      read ( unit=unit_tfemu ) &
                     mesh%volumes(volume)%element_blend(:nblend)%elshape
      read ( unit=unit_tfemu ) &
                     mesh%volumes(volume)%element_blend(:nblend)%ndim
      do m = 1, nblend
        read ( unit=unit_tfemu ) mesh%volumes(volume)%element_blend(m)%p
      end do
    end do

!   blocks

    allocate( mesh%blocks(max(mesh%nblocks,MAXBLOCKS)) )

    read ( unit=unit_tfemu ) mesh%blocks(:mesh%nblocks)%nelem

    do block = 1, mesh%nblocks
      allocate ( mesh%blocks(block)%bounds(mesh%ndim,2) )
      read ( unit=unit_tfemu ) mesh%blocks(block)%bounds
      nelem = mesh%blocks(block)%nelem
      allocate ( mesh%blocks(block)%elbounds(nelem,mesh%ndim,2) )
      read ( unit=unit_tfemu ) mesh%blocks(block)%elbounds
      allocate ( mesh%blocks(block)%elements(nelem,2) )
      read ( unit=unit_tfemu ) mesh%blocks(block)%elements
    end do

!   objects

    allocate( mesh%objects(max(mesh%nobjects,MAXOBJECTS)) )

    read ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%topol
    read ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%intpoints
    read ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%typeofobject
    read ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%ndim
    read ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%nnodes
    read ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%nodeset
    read ( unit=unit_tfemu ) mesh%objects(:mesh%nobjects)%nodeset2

    do object = 1, mesh%nobjects
      allocate( mesh%objects(object)%groups(mesh%nelgrp) )
      read ( unit=unit_tfemu ) mesh%objects(object)%groups
      allocate( mesh%objects(object)%coor(mesh%objects(object)%nnodes, &
                                          mesh%objects(object)%ndim) )
      read ( unit=unit_tfemu ) mesh%objects(object)%coor
      allocate( mesh%objects(object)%refcoor(mesh%objects(object)%nnodes, &
                                              mesh%objects(object)%ndim) )
      allocate( mesh%objects(object)%grpelm(mesh%objects(object)%nnodes,2) )
      if ( mesh%objects(object)%typeofobject == 2 ) then
        allocate( mesh%objects(object)%groups2(mesh%nelgrp) )
        read ( unit=unit_tfemu ) mesh%objects(object)%groups2
        allocate( mesh%objects(object)%refcoor2(mesh%objects(object)%nnodes, &
                                                mesh%objects(object)%ndim) )
        allocate( mesh%objects(object)%grpelm2(mesh%objects(object)%nnodes,2) )
      end if
      if ( mesh%objects(object)%topol ) then
        read ( unit=unit_tfemu ) mesh%objects(object)%nelem
        read ( unit=unit_tfemu ) mesh%objects(object)%elnumnod
        read ( unit=unit_tfemu ) mesh%objects(object)%element%globalshape
        read ( unit=unit_tfemu ) mesh%objects(object)%element%numnod
        read ( unit=unit_tfemu ) mesh%objects(object)%element%elshape
        read ( unit=unit_tfemu ) mesh%objects(object)%element%p
        read ( unit=unit_tfemu ) mesh%objects(object)%element%ndim
        allocate( mesh%objects(object)%topology(mesh%objects(object)%elnumnod,&
        &mesh%objects(object)%nelem) )
        read ( unit=unit_tfemu ) mesh%objects(object)%topology
      end if
      if ( mesh%objects(object)%intpoints ) then
        read ( unit=unit_tfemu ) mesh%objects(object)%ninti
        read ( unit=unit_tfemu ) ndimxi
        allocate( mesh%objects(object)%xig(mesh%objects(object)%ninti,ndimxi) )
        allocate( mesh%objects(object)%wg(mesh%objects(object)%ninti) )
        read ( unit=unit_tfemu ) mesh%objects(object)%xig
        read ( unit=unit_tfemu ) mesh%objects(object)%wg
!       allocate more data
        allocate ( &
          mesh%objects(object)%coor_int(mesh%objects(object)%ninti,&
            mesh%objects(object)%ndim,mesh%objects(object)%nelem), &
          mesh%objects(object)%refcoor_int(mesh%objects(object)%ninti,&
            mesh%objects(object)%ndim,mesh%objects(object)%nelem), &
          mesh%objects(object)%grpelm_int(mesh%objects(object)%ninti,2,&
            mesh%objects(object)%nelem) )
        if ( mesh%objects(object)%typeofobject == 2 ) then
!         second intersection (typeofobject==2)
          allocate ( &
            mesh%objects(object)%refcoor2_int(mesh%objects(object)%ninti,&
              mesh%objects(object)%ndim,mesh%objects(object)%nelem), &
            mesh%objects(object)%grpelm2_int(mesh%objects(object)%ninti,2,&
              mesh%objects(object)%nelem) )
        end if
      end if
    end do

!   nodesets

    allocate( mesh%nodesets(max(mesh%nnodesets,MAXNODESETS)) )

    do nodeset = 1, mesh%nnodesets
      read ( unit=unit_tfemu ) nnodes
      allocate(mesh%nodesets(nodeset)%a(nnodes))
      read ( unit=unit_tfemu ) mesh%nodesets(nodeset)%a
    end do

!   elementsets

    allocate ( mesh%elementsets(max(mesh%nelementsets,MAXELEMENTSETS)) )
    do elementset = 1, mesh%nelementsets
      read ( unit=unit_tfemu ) mesh%elementsets(elementset)%nodes_created
      read ( unit=unit_tfemu ) mesh%elementsets(elementset)%nelem
      read ( unit=unit_tfemu ) nnodes
      mesh%elementsets(elementset)%nnodes = nnodes
      allocate ( mesh%elementsets(elementset)%nodes(nnodes) )
      read ( unit=unit_tfemu ) mesh%elementsets(elementset)%nodes
      allocate ( mesh%elementsets(elementset)%grpnumel(mesh%nelgrp) )
      read ( unit=unit_tfemu ) mesh%elementsets(elementset)%grpnumel
      allocate ( mesh%elementsets(elementset)%elements(mesh%nelgrp) )
      do elgrp = 1, mesh%nelgrp
        read ( unit=unit_tfemu ) nelem
        allocate ( mesh%elementsets(elementset)%elements(elgrp)%a(nelem) )
        read ( unit=unit_tfemu ) mesh%elementsets(elementset)%elements(elgrp)%a
      end do
    end do

!   sblocks

    if ( mesh%sblocks%filled ) then

      read ( unit=unit_tfemu ) mesh%sblocks%nsblocks

      allocate ( mesh%sblocks%nbl(mesh%ndim), mesh%sblocks%x0(mesh%ndim), &
                 mesh%sblocks%dx(mesh%ndim) )

      read ( unit=unit_tfemu ) mesh%sblocks%nbl
      read ( unit=unit_tfemu ) mesh%sblocks%x0
      read ( unit=unit_tfemu ) mesh%sblocks%dx

      allocate ( mesh%sblocks%sblks(mesh%sblocks%nsblocks) )

      read ( unit=unit_tfemu ) mesh%sblocks%sblks(:)%nelem
      do sblock = 1, mesh%sblocks%nsblocks
        nelem = mesh%sblocks%sblks(sblock)%nelem
        allocate ( mesh%sblocks%sblks(sblock)%elbounds(nelem,mesh%ndim,2) )
        read ( unit=unit_tfemu ) mesh%sblocks%sblks(sblock)%elbounds
        allocate ( mesh%sblocks%sblks(sblock)%elements(nelem,2) )
        read ( unit=unit_tfemu ) mesh%sblocks%sblks(sblock)%elements
      end do

    end if

!   nodblocks

    if ( mesh%nodblocks%filled ) then

      read ( unit=unit_tfemu ) mesh%nodblocks%nnodblocks

      allocate ( mesh%nodblocks%nbl(mesh%ndim), mesh%nodblocks%x0(mesh%ndim), &
                 mesh%nodblocks%dx(mesh%ndim) )

      read ( unit=unit_tfemu ) mesh%nodblocks%nbl
      read ( unit=unit_tfemu ) mesh%nodblocks%x0
      read ( unit=unit_tfemu ) mesh%nodblocks%dx

      allocate ( mesh%nodblocks%nodblks(mesh%nodblocks%nnodblocks) )

      read ( unit=unit_tfemu ) mesh%nodblocks%nodblks(:)%nnodes

      do nodblock = 1, mesh%nodblocks%nnodblocks
        nnodes = mesh%nodblocks%nodblks(nodblock)%nnodes
        allocate ( mesh%nodblocks%nodblks(nodblock)%nodes(nnodes) )
        read ( unit=unit_tfemu ) mesh%nodblocks%nodblks(nodblock)%nodes
      end do

    end if

!   coordinates

    allocate( mesh%coor(mesh%nnodes,mesh%ndim) )

    read ( unit=unit_tfemu ) mesh%coor

!   mesh renumbering

    if ( mesh%renumber ) then
      allocate ( mesh%nodperm(mesh%nnodes) )
      read ( unit=unit_tfemu ) mesh%nodperm
    end if

    mesh%meshgen = .true.

    close ( unit=unit_tfemu )

  end subroutine read_mesh


! write input_probdef to a file

  subroutine write_input_probdef ( mesh, input_probdef, filename )

    type(mesh_t), intent(in) :: mesh
    type(input_probdef_t), intent(in) :: input_probdef

!   the filename for writing input_probdef to
    character (len=*), intent(in) :: filename

!   This routine writes the input_probdef structure to an intermediate file
!   for later reading by read_input_probdef from a different main program.
!   The file format is therefore binary and not intended for a permanent
!   storage file since it might change in the future.

    integer :: elgrp, esspnt, esscrv, esssrf, essgrp, constr, esselm, conn
    integer :: essnst, ios, trans, dep


    if ( .not. input_probdef%created ) then
      write(*,'(/a/)') 'Error in write_input_probdef', &
        ' input_probdef has not been created '
      stop
    end if

!   open file

    open ( unit=unit_tfemu, file=filename, form='unformatted', iostat=ios, &
      status='replace' )

    if ( ios /= 0 ) then
      write(*,'(/a,i0/)') 'Error in write_input_probdef: ios = ', ios
      stop
    end if

!   some scalar data

    write ( unit=unit_tfemu ) input_probdef%nvec, input_probdef%nphysq, &
      input_probdef%nphysqshifted, input_probdef%numinactivegroups, &
      input_probdef%numlayers
    write ( unit=unit_tfemu ) input_probdef%probnr, input_probdef%numesspoints,&
      input_probdef%numesscurves, input_probdef%numesssurfaces, &
      input_probdef%numesselements, input_probdef%numconstraints, &
      input_probdef%numessgroups, input_probdef%numessnodesets, &
      input_probdef%numconnections, input_probdef%numtransformations, &
      input_probdef%numdependencies

!   elementdata

    do elgrp = 1, mesh%nelgrp
      write ( unit=unit_tfemu ) input_probdef%elementdof(elgrp)%a
      write ( unit=unit_tfemu ) input_probdef%vec_elementdof(elgrp)%a
    end do

!   physical quantities

    write ( unit=unit_tfemu ) input_probdef%physq, input_probdef%physqshifted
    write ( unit=unit_tfemu ) input_probdef%physqmask

!   layers

    write ( unit=unit_tfemu ) input_probdef%layers

!   inactive groups

    write ( unit=unit_tfemu ) input_probdef%inactivegroups

!   esspoints

    do esspnt = 1, input_probdef%numesspoints
      write ( unit=unit_tfemu ) input_probdef%esspoints(esspnt)%physq,   &
                          input_probdef%esspoints(esspnt)%layer,   &
                          input_probdef%esspoints(esspnt)%point,   &
                          size(input_probdef%esspoints(esspnt)%points), &
                          size(input_probdef%esspoints(esspnt)%degfd)
      write ( unit=unit_tfemu ) input_probdef%esspoints(esspnt)%points, &
                                input_probdef%esspoints(esspnt)%degfd
    end do

!   esscurves

    do esscrv = 1, input_probdef%numesscurves
      write ( unit=unit_tfemu ) input_probdef%esscurves(esscrv)%physq,   &
                          input_probdef%esscurves(esscrv)%layer,   &
                          input_probdef%esscurves(esscrv)%first,   &
                          input_probdef%esscurves(esscrv)%last,   &
                          input_probdef%esscurves(esscrv)%step,   &
                          input_probdef%esscurves(esscrv)%exclude,   &
                          size(input_probdef%esscurves(esscrv)%curves), &
                      size(input_probdef%esscurves(esscrv)%excludepoints), &
                      size(input_probdef%esscurves(esscrv)%excludecurves), &
                      size(input_probdef%esscurves(esscrv)%excludesurfaces), &
                          size(input_probdef%esscurves(esscrv)%degfd)
      write ( unit=unit_tfemu ) &
                          input_probdef%esscurves(esscrv)%curves, &
                          input_probdef%esscurves(esscrv)%excludepoints, &
                          input_probdef%esscurves(esscrv)%excludecurves, &
                          input_probdef%esscurves(esscrv)%excludesurfaces, &
                          input_probdef%esscurves(esscrv)%degfd
    end do

!   esssurfaces

    do esssrf = 1, input_probdef%numesssurfaces
      write ( unit=unit_tfemu ) input_probdef%esssurfaces(esssrf)%physq,   &
                          input_probdef%esssurfaces(esssrf)%layer,   &
                          input_probdef%esssurfaces(esssrf)%first,   &
                          input_probdef%esssurfaces(esssrf)%last,    &
                          size(input_probdef%esssurfaces(esssrf)%surfaces), &
                      size(input_probdef%esssurfaces(esssrf)%excludepoints), &
                      size(input_probdef%esssurfaces(esssrf)%excludecurves), &
                      size(input_probdef%esssurfaces(esssrf)%excludesurfaces), &
                          size(input_probdef%esssurfaces(esssrf)%degfd)
      write ( unit=unit_tfemu ) &
                          input_probdef%esssurfaces(esssrf)%surfaces, &
                          input_probdef%esssurfaces(esssrf)%excludepoints, &
                          input_probdef%esssurfaces(esssrf)%excludecurves, &
                          input_probdef%esssurfaces(esssrf)%excludesurfaces, &
                          input_probdef%esssurfaces(esssrf)%degfd
    end do

!   esselements

    do esselm = 1, input_probdef%numesselements
      write ( unit=unit_tfemu ) input_probdef%esselements(esselm)%physq,   &
                          input_probdef%esselements(esselm)%layer,   &
                          input_probdef%esselements(esselm)%elgrp,   &
                          input_probdef%esselements(esselm)%elem,   &
                          input_probdef%esselements(esselm)%node,   &
                          size(input_probdef%esselements(esselm)%degfd)
      write ( unit=unit_tfemu ) input_probdef%esselements(esselm)%degfd
    end do

!   essgroups

    do essgrp = 1, input_probdef%numessgroups
      write ( unit=unit_tfemu ) input_probdef%essgroups(essgrp)%physq,   &
                          input_probdef%essgroups(essgrp)%layer,   &
                          input_probdef%essgroups(essgrp)%first,   &
                          input_probdef%essgroups(essgrp)%last,    &
                          size(input_probdef%essgroups(essgrp)%elgroups), &
                      size(input_probdef%essgroups(essgrp)%excludepoints), &
                      size(input_probdef%essgroups(essgrp)%excludecurves), &
                      size(input_probdef%essgroups(essgrp)%excludesurfaces), &
                          size(input_probdef%essgroups(essgrp)%degfd)
      write ( unit=unit_tfemu ) &
                          input_probdef%essgroups(essgrp)%elgroups, &
                          input_probdef%essgroups(essgrp)%excludepoints, &
                          input_probdef%essgroups(essgrp)%excludecurves, &
                          input_probdef%essgroups(essgrp)%excludesurfaces, &
                          input_probdef%essgroups(essgrp)%degfd
    end do

!   essnodesets

    do essnst = 1, input_probdef%numessnodesets
      write ( unit=unit_tfemu ) input_probdef%essnodesets(essnst)%physq,   &
                          input_probdef%essnodesets(essnst)%layer,   &
                          input_probdef%essnodesets(essnst)%first,   &
                          input_probdef%essnodesets(essnst)%last,    &
                          size(input_probdef%essnodesets(essnst)%nodesets), &
                      size(input_probdef%essnodesets(essnst)%excludepoints), &
                      size(input_probdef%essnodesets(essnst)%excludecurves), &
                      size(input_probdef%essnodesets(essnst)%excludesurfaces), &
                          size(input_probdef%essnodesets(essnst)%degfd)
      write ( unit=unit_tfemu ) &
                          input_probdef%essnodesets(essnst)%nodesets, &
                          input_probdef%essnodesets(essnst)%excludepoints, &
                          input_probdef%essnodesets(essnst)%excludecurves, &
                          input_probdef%essnodesets(essnst)%excludesurfaces, &
                          input_probdef%essnodesets(essnst)%degfd
    end do

!   constraints

    do constr = 1, input_probdef%numconstraints
      write ( unit=unit_tfemu ) &
                          input_probdef%constraints(constr)%typeconstraint, &
                          input_probdef%constraints(constr)%physq1, &
                          input_probdef%constraints(constr)%physq2, &
                          input_probdef%constraints(constr)%layer1, &
                          input_probdef%constraints(constr)%layer2, &
                          input_probdef%constraints(constr)%errorlayer, &
                          input_probdef%constraints(constr)%typegeometry, &
                          input_probdef%constraints(constr)%geometry1, &
                          input_probdef%constraints(constr)%geometry2, &
                          input_probdef%constraints(constr)%elementset1, &
                          input_probdef%constraints(constr)%elementset2, &
                          input_probdef%constraints(constr)%object, &
                          input_probdef%constraints(constr)%object2, &
                          input_probdef%constraints(constr)%full2, &
                          input_probdef%constraints(constr)%discretization, &
                          input_probdef%constraints(constr)%step, &
                          input_probdef%constraints(constr)%exclude, &
                          input_probdef%constraints(constr)%nodenumdegfd, &
                          input_probdef%constraints(constr)%nglobalc, &
                          input_probdef%constraints(constr)%naddunknowns, &
                      size(input_probdef%constraints(constr)%excludepoints), &
                      size(input_probdef%constraints(constr)%excludecurves), &
                      size(input_probdef%constraints(constr)%excludesurfaces), &
                          size(input_probdef%constraints(constr)%elnumdegfd)
      write ( unit=unit_tfemu ) &
                          input_probdef%constraints(constr)%excludepoints, &
                          input_probdef%constraints(constr)%excludecurves, &
                          input_probdef%constraints(constr)%excludesurfaces, &
                          input_probdef%constraints(constr)%elnumdegfd
    end do

!   connections

    do conn = 1, input_probdef%numconnections
      write ( unit=unit_tfemu ) &
                          input_probdef%connections(conn)%physq1, &
                          input_probdef%connections(conn)%physq2, &
                          input_probdef%connections(conn)%layer1, &
                          input_probdef%connections(conn)%layer2, &
                          input_probdef%connections(conn)%typegeometry, &
                          input_probdef%connections(conn)%typegeometry2, &
                          input_probdef%connections(conn)%geometry1, &
                          input_probdef%connections(conn)%geometry2, &
                          input_probdef%connections(conn)%elementset1, &
                          input_probdef%connections(conn)%elementset2, &
                          input_probdef%connections(conn)%object, &
                          input_probdef%connections(conn)%object2, &
                          input_probdef%connections(conn)%discretization
    end do

!   transformations

    do trans = 1, input_probdef%numtransformations
      write ( unit=unit_tfemu ) &
                      input_probdef%transformations(trans)%typetransformation, &
                      input_probdef%transformations(trans)%orthogonal, &
                      input_probdef%transformations(trans)%normalvector, &
                      input_probdef%transformations(trans)%v2, &
                      input_probdef%transformations(trans)%physq, &
                      input_probdef%transformations(trans)%layer, &
                      input_probdef%transformations(trans)%errorlayer, &
                      input_probdef%transformations(trans)%typegeometry, &
                      input_probdef%transformations(trans)%geometry, &
                      input_probdef%transformations(trans)%step, &
                      input_probdef%transformations(trans)%exclude, &
                  size(input_probdef%transformations(trans)%excludepoints), &
                  size(input_probdef%transformations(trans)%excludecurves), &
                  size(input_probdef%transformations(trans)%excludesurfaces), &
                  size(input_probdef%transformations(trans)%Amat_global,1)
      write ( unit=unit_tfemu ) &
                      input_probdef%transformations(trans)%excludepoints, &
                      input_probdef%transformations(trans)%excludecurves, &
                      input_probdef%transformations(trans)%excludesurfaces, &
                      input_probdef%transformations(trans)%Amat_global
    end do

!   dependencies

     do dep = 1, input_probdef%numdependencies
      write ( unit=unit_tfemu ) &
                      input_probdef%dependencies(dep)%typedependency, &
                      input_probdef%dependencies(dep)%datalayout, &
                      input_probdef%dependencies(dep)%physq1, &
                      input_probdef%dependencies(dep)%physq2, &
                      input_probdef%dependencies(dep)%layer1, &
                      input_probdef%dependencies(dep)%layer2, &
                      input_probdef%dependencies(dep)%typegeometry1, &
                      input_probdef%dependencies(dep)%typegeometry2, &
                      input_probdef%dependencies(dep)%geometry1, &
                      input_probdef%dependencies(dep)%geometry2, &
                      input_probdef%dependencies(dep)%elementset1, &
                      input_probdef%dependencies(dep)%elementset2, &
                      input_probdef%dependencies(dep)%object2, &
                      input_probdef%dependencies(dep)%step, &
                      input_probdef%dependencies(dep)%exclude, &
                      input_probdef%dependencies(dep)%naddunknowns, &
                  size(input_probdef%dependencies(dep)%excludepoints), &
                  size(input_probdef%dependencies(dep)%excludecurves), &
                  size(input_probdef%dependencies(dep)%excludesurfaces)
      write ( unit=unit_tfemu ) &
                      input_probdef%dependencies(dep)%excludepoints, &
                      input_probdef%dependencies(dep)%excludecurves, &
                      input_probdef%dependencies(dep)%excludesurfaces
    end do

    close ( unit=unit_tfemu )

  end subroutine write_input_probdef


! read input_probdef from a file written by write_input_probdef

  subroutine read_input_probdef ( mesh, input_probdef, filename )

    type(mesh_t), intent(in) :: mesh
    type(input_probdef_t), intent(out) :: input_probdef

!   the filename for writing input_probdef to
    character (len=*), intent(in) :: filename

!   This routine reads the input_probdef structure from an intermediate file
!   that is written by write_input_probdef, possibly in a different program.

    integer :: elgrp, esspnt, esscrv, esssrf, essgrp, constr, ios, nvec, nphysq
    integer :: nphysqshifted, numinactivegroups, conn, trans, dep
    integer :: sizeexpoints, sizeexcurves, sizeexsurfaces, sizedegfd, sizeAmat
    integer :: sizepoints, sizecurves, sizesurfaces, sizeelgroups, sizenodesets
    integer :: numlayers, esselm, essnst


!   open file

    open ( unit=unit_tfemu, file=filename, form='unformatted', iostat=ios, &
      status='old' )

    if ( ios /= 0 ) then
      write(*,'(/2a/)') 'Error in read_input_probdef: cannot open file ', &
        filename
      stop
    end if

!   some scalar data

    read ( unit=unit_tfemu ) nvec, nphysq, nphysqshifted, numinactivegroups, &
      numlayers

!   create structure

    call create_input_probdef ( mesh, input_probdef, nvec, nphysq, &
      nphysqshifted, numinactivegroups, numlayers )

!   some more scalar data

    read ( unit=unit_tfemu ) input_probdef%probnr, input_probdef%numesspoints, &
      input_probdef%numesscurves, input_probdef%numesssurfaces, &
      input_probdef%numesselements, input_probdef%numconstraints, &
      input_probdef%numessgroups, input_probdef%numessnodesets, &
      input_probdef%numconnections, input_probdef%numtransformations, &
      input_probdef%numdependencies

!   elementdata

    do elgrp = 1, mesh%nelgrp
      read ( unit=unit_tfemu ) input_probdef%elementdof(elgrp)%a
      read ( unit=unit_tfemu ) input_probdef%vec_elementdof(elgrp)%a
    end do

!   physical quantities

    read ( unit=unit_tfemu ) input_probdef%physq, input_probdef%physqshifted
    read ( unit=unit_tfemu ) input_probdef%physqmask

!   layers

    read ( unit=unit_tfemu ) input_probdef%layers

!   inactive groups

    read ( unit=unit_tfemu ) input_probdef%inactivegroups

!   esspoints

    do esspnt = 1, input_probdef%numesspoints
      read ( unit=unit_tfemu ) input_probdef%esspoints(esspnt)%physq,   &
                         input_probdef%esspoints(esspnt)%layer,   &
                         input_probdef%esspoints(esspnt)%point,   &
                         sizepoints, sizedegfd
      allocate(input_probdef%esspoints(esspnt)%points(sizepoints))
      allocate(input_probdef%esspoints(esspnt)%degfd(sizedegfd))
      read ( unit=unit_tfemu ) input_probdef%esspoints(esspnt)%points, &
                               input_probdef%esspoints(esspnt)%degfd
    end do

!   esscurves

    do esscrv = 1, input_probdef%numesscurves
      read ( unit=unit_tfemu ) input_probdef%esscurves(esscrv)%physq,   &
                         input_probdef%esscurves(esscrv)%layer,   &
                         input_probdef%esscurves(esscrv)%first,   &
                         input_probdef%esscurves(esscrv)%last,    &
                         input_probdef%esscurves(esscrv)%step,    &
                         input_probdef%esscurves(esscrv)%exclude, &
                         sizecurves, &
                         sizeexpoints, sizeexcurves, sizeexsurfaces, &
                         sizedegfd
      allocate(input_probdef%esscurves(esscrv)%curves(sizecurves))
      allocate(input_probdef%esscurves(esscrv)%excludepoints(sizeexpoints))
      allocate(input_probdef%esscurves(esscrv)%excludecurves(sizeexcurves))
      allocate(&
             input_probdef%esscurves(esscrv)%excludesurfaces(sizeexsurfaces))
      allocate(input_probdef%esscurves(esscrv)%degfd(sizedegfd))
      read ( unit=unit_tfemu ) &
                         input_probdef%esscurves(esscrv)%curves, &
                         input_probdef%esscurves(esscrv)%excludepoints, &
                         input_probdef%esscurves(esscrv)%excludecurves, &
                         input_probdef%esscurves(esscrv)%excludesurfaces, &
                         input_probdef%esscurves(esscrv)%degfd
    end do

!   esssurfaces

    do esssrf = 1, input_probdef%numesssurfaces
      read ( unit=unit_tfemu ) input_probdef%esssurfaces(esssrf)%physq,   &
                         input_probdef%esssurfaces(esssrf)%layer,   &
                         input_probdef%esssurfaces(esssrf)%first,   &
                         input_probdef%esssurfaces(esssrf)%last,   &
                         sizesurfaces, &
                         sizeexpoints, sizeexcurves, sizeexsurfaces, &
                         sizedegfd
      allocate(input_probdef%esssurfaces(esssrf)%surfaces(sizesurfaces))
      allocate(input_probdef%esssurfaces(esssrf)%excludepoints(sizeexpoints))
      allocate(input_probdef%esssurfaces(esssrf)%excludecurves(sizeexcurves))
      allocate(&
             input_probdef%esssurfaces(esssrf)%excludesurfaces(sizeexsurfaces))
      allocate(input_probdef%esssurfaces(esssrf)%degfd(sizedegfd))
      read ( unit=unit_tfemu ) &
                         input_probdef%esssurfaces(esssrf)%surfaces, &
                         input_probdef%esssurfaces(esssrf)%excludepoints, &
                         input_probdef%esssurfaces(esssrf)%excludecurves, &
                         input_probdef%esssurfaces(esssrf)%excludesurfaces, &
                         input_probdef%esssurfaces(esssrf)%degfd
    end do

!   esselements

    do esselm = 1, input_probdef%numesselements
      read ( unit=unit_tfemu ) input_probdef%esselements(esselm)%physq,   &
                         input_probdef%esselements(esselm)%layer,   &
                         input_probdef%esselements(esselm)%elgrp,   &
                         input_probdef%esselements(esselm)%elem,   &
                         input_probdef%esselements(esselm)%node,   &
                         sizedegfd
      allocate(input_probdef%esselements(esselm)%degfd(sizedegfd))
      read ( unit=unit_tfemu ) input_probdef%esselements(esselm)%degfd
    end do

!   essgroups

    do essgrp = 1, input_probdef%numessgroups
      read ( unit=unit_tfemu ) input_probdef%essgroups(essgrp)%physq,   &
                         input_probdef%essgroups(essgrp)%layer,   &
                         input_probdef%essgroups(essgrp)%first,   &
                         input_probdef%essgroups(essgrp)%last,   &
                         sizeelgroups, &
                         sizeexpoints, sizeexcurves, sizeexsurfaces, &
                         sizedegfd
      allocate(input_probdef%essgroups(essgrp)%degfd(sizeelgroups))
      allocate(input_probdef%essgroups(essgrp)%excludepoints(sizeexpoints))
      allocate(input_probdef%essgroups(essgrp)%excludecurves(sizeexcurves))
      allocate(input_probdef%essgroups(essgrp)%excludesurfaces(sizeexsurfaces))
      allocate(input_probdef%essgroups(essgrp)%degfd(sizedegfd))
      read ( unit=unit_tfemu ) &
                         input_probdef%essgroups(essgrp)%elgroups, &
                         input_probdef%essgroups(essgrp)%excludepoints, &
                         input_probdef%essgroups(essgrp)%excludecurves, &
                         input_probdef%essgroups(essgrp)%excludesurfaces, &
                         input_probdef%essgroups(essgrp)%degfd
    end do

!   essnodesets

    do essnst = 1, input_probdef%numessnodesets
      read ( unit=unit_tfemu ) input_probdef%essnodesets(essnst)%physq,   &
                         input_probdef%essnodesets(essnst)%layer,   &
                         input_probdef%essnodesets(essnst)%first,   &
                         input_probdef%essnodesets(essnst)%last,   &
                         sizenodesets, &
                         sizeexpoints, sizeexcurves, sizeexsurfaces, &
                         sizedegfd
      allocate(input_probdef%essnodesets(essnst)%nodesets(sizenodesets))
      allocate(input_probdef%essnodesets(essnst)%excludepoints(sizeexpoints))
      allocate(input_probdef%essnodesets(essnst)%excludecurves(sizeexcurves))
      allocate(&
             input_probdef%essnodesets(essnst)%excludesurfaces(sizeexsurfaces))
      allocate(input_probdef%essnodesets(essnst)%degfd(sizedegfd))
      read ( unit=unit_tfemu ) &
                         input_probdef%essnodesets(essnst)%nodesets, &
                         input_probdef%essnodesets(essnst)%excludepoints, &
                         input_probdef%essnodesets(essnst)%excludecurves, &
                         input_probdef%essnodesets(essnst)%excludesurfaces, &
                         input_probdef%essnodesets(essnst)%degfd
    end do

!   constraints

    do constr = 1, input_probdef%numconstraints
      read ( unit=unit_tfemu ) &
                         input_probdef%constraints(constr)%typeconstraint, &
                         input_probdef%constraints(constr)%physq1, &
                         input_probdef%constraints(constr)%physq2, &
                         input_probdef%constraints(constr)%layer1, &
                         input_probdef%constraints(constr)%layer2, &
                         input_probdef%constraints(constr)%errorlayer, &
                         input_probdef%constraints(constr)%typegeometry, &
                         input_probdef%constraints(constr)%geometry1, &
                         input_probdef%constraints(constr)%geometry2, &
                         input_probdef%constraints(constr)%elementset1, &
                         input_probdef%constraints(constr)%elementset2, &
                         input_probdef%constraints(constr)%object, &
                         input_probdef%constraints(constr)%object2, &
                         input_probdef%constraints(constr)%full2, &
                         input_probdef%constraints(constr)%discretization, &
                         input_probdef%constraints(constr)%step, &
                         input_probdef%constraints(constr)%exclude, &
                         input_probdef%constraints(constr)%nodenumdegfd, &
                         input_probdef%constraints(constr)%nglobalc, &
                         input_probdef%constraints(constr)%naddunknowns, &
                         sizeexpoints, sizeexcurves, sizeexsurfaces, sizedegfd
      allocate(input_probdef%constraints(constr)%excludepoints(sizeexpoints))
      allocate(input_probdef%constraints(constr)%excludecurves(sizeexcurves))
      allocate(&
             input_probdef%constraints(constr)%excludesurfaces(sizeexsurfaces))
      allocate(input_probdef%constraints(constr)%elnumdegfd(sizedegfd))
      read ( unit=unit_tfemu ) &
                         input_probdef%constraints(constr)%excludepoints, &
                         input_probdef%constraints(constr)%excludecurves, &
                         input_probdef%constraints(constr)%excludesurfaces, &
                         input_probdef%constraints(constr)%elnumdegfd
    end do

!   connections

    do conn = 1, input_probdef%numconnections
      read ( unit=unit_tfemu ) &
                          input_probdef%connections(conn)%physq1, &
                          input_probdef%connections(conn)%physq2, &
                          input_probdef%connections(conn)%layer1, &
                          input_probdef%connections(conn)%layer2, &
                          input_probdef%connections(conn)%typegeometry, &
                          input_probdef%connections(conn)%typegeometry2, &
                          input_probdef%connections(conn)%geometry1, &
                          input_probdef%connections(conn)%geometry2, &
                          input_probdef%connections(conn)%elementset1, &
                          input_probdef%connections(conn)%elementset2, &
                          input_probdef%connections(conn)%object, &
                          input_probdef%connections(conn)%object2, &
                          input_probdef%connections(conn)%discretization
    end do

!   transformations

    do trans = 1, input_probdef%numtransformations
      input_probdef%transformations(trans)%build = .false.
      read ( unit=unit_tfemu ) &
                      input_probdef%transformations(trans)%typetransformation, &
                      input_probdef%transformations(trans)%orthogonal, &
                      input_probdef%transformations(trans)%normalvector, &
                      input_probdef%transformations(trans)%v2, &
                      input_probdef%transformations(trans)%physq, &
                      input_probdef%transformations(trans)%layer, &
                      input_probdef%transformations(trans)%errorlayer, &
                      input_probdef%transformations(trans)%typegeometry, &
                      input_probdef%transformations(trans)%geometry, &
                      input_probdef%transformations(trans)%step, &
                      input_probdef%transformations(trans)%exclude, &
                      sizeexpoints, sizeexcurves, sizeexsurfaces, sizeAmat
      allocate(input_probdef%transformations(trans)%excludepoints(sizeexpoints))
      allocate(input_probdef%transformations(trans)%excludecurves(sizeexcurves))
      allocate(&
           input_probdef%transformations(trans)%excludesurfaces(sizeexsurfaces))
      allocate(&
           input_probdef%transformations(trans)%Amat_global(sizeAmat,sizeAmat))
      read ( unit=unit_tfemu ) &
                      input_probdef%transformations(trans)%excludepoints, &
                      input_probdef%transformations(trans)%excludecurves, &
                      input_probdef%transformations(trans)%excludesurfaces, &
                      input_probdef%transformations(trans)%Amat_global
    end do

!   dependencies

    do dep = 1, input_probdef%numdependencies
      input_probdef%dependencies(dep)%build = .false.
      read ( unit=unit_tfemu ) &
                      input_probdef%dependencies(dep)%typedependency, &
                      input_probdef%dependencies(dep)%datalayout, &
                      input_probdef%dependencies(dep)%physq1, &
                      input_probdef%dependencies(dep)%physq2, &
                      input_probdef%dependencies(dep)%layer1, &
                      input_probdef%dependencies(dep)%layer2, &
                      input_probdef%dependencies(dep)%typegeometry1, &
                      input_probdef%dependencies(dep)%typegeometry2, &
                      input_probdef%dependencies(dep)%geometry1, &
                      input_probdef%dependencies(dep)%geometry2, &
                      input_probdef%dependencies(dep)%object2, &
                      input_probdef%dependencies(dep)%elementset1, &
                      input_probdef%dependencies(dep)%elementset2, &
                      input_probdef%dependencies(dep)%step, &
                      input_probdef%dependencies(dep)%exclude, &
                      input_probdef%dependencies(dep)%naddunknowns, &
                      sizeexpoints, sizeexcurves, sizeexsurfaces
      allocate(input_probdef%dependencies(dep)%excludepoints(sizeexpoints))
      allocate(input_probdef%dependencies(dep)%excludecurves(sizeexcurves))
      allocate(input_probdef%dependencies(dep)%excludesurfaces(sizeexsurfaces))
      read ( unit=unit_tfemu ) &
                      input_probdef%dependencies(dep)%excludepoints, &
                      input_probdef%dependencies(dep)%excludecurves, &
                      input_probdef%dependencies(dep)%excludesurfaces
    end do

    close ( unit=unit_tfemu )

  end subroutine read_input_probdef


! write coefficients to a file

  subroutine write_coefficients ( coefficients, filename )

    type(coefficients_t), intent(in) :: coefficients

!   the filename for writing the cofficients to
    character (len=*), intent(in) :: filename

!   This routine writes to an intermediate file for later reading by
!   read_coefficients from a different main program. The file format is
!   therefore binary and not intended for a permanent storage file since
!   it might change in the future.

    integer :: ios, i, ni, nr


    if ( .not. coefficients%created ) then
      write(*,'(/a/)') 'Error in write_coefficients: ', &
        ' coefficients not created '
      stop
    end if

!   open file

    open ( unit=unit_tfemu, file=filename, form='unformatted', iostat=ios, &
      status='replace' )

    if ( ios /= 0 ) then
      write(*,'(/a,i0/)') 'Error in write_coefficients: ios = ', ios
      stop
    end if

!   sizes

    write ( unit=unit_tfemu ) size(coefficients%i), size(coefficients%r)
    ni = size(coefficients%ia)
    nr = size(coefficients%ra)
    write ( unit=unit_tfemu ) ni, nr
    write ( unit=unit_tfemu ) &
      ( size(coefficients%ia(i)%a), i=1,ni ), &
      ( size(coefficients%ra(i)%a), i=1,nr )

!   data

    write ( unit=unit_tfemu ) coefficients%i
    write ( unit=unit_tfemu ) coefficients%r
    write ( unit=unit_tfemu ) &
      ( coefficients%ia(i)%a, i=1,ni ), &
      ( coefficients%ra(i)%a, i=1,nr )

    close ( unit=unit_tfemu )

  end subroutine write_coefficients


! read coefficients from a file

  subroutine read_coefficients ( coefficients, filename )

    type(coefficients_t), intent(inout) :: coefficients

!   the filename for writing the cofficients to
    character (len=*), intent(in) :: filename

    integer :: ios, ncoefi, ncoefr, ni, nr, i
    integer, dimension(:), allocatable :: ncoefia, ncoefra


    if ( coefficients%created ) then
      write(*,'(/a/)') 'Error in read_coefficients: ', &
        ' coefficients already created '
      stop
    end if

!   open file

    open ( unit=unit_tfemu, file=filename, form='unformatted', iostat=ios, &
      status='old' )

    if ( ios /= 0 ) then
      write(*,'(/2a/)') 'Error in read_coefficients: cannot open file ', &
        filename
      stop
    end if

!   read sizes

    read ( unit=unit_tfemu ) ncoefi, ncoefr
    read ( unit=unit_tfemu ) ni, nr

    allocate ( ncoefia(ni), ncoefra(nr) )

    read ( unit=unit_tfemu ) ncoefia, ncoefra

!   create structure

    call create_coefficients ( coefficients, ncoefi=ncoefi, ncoefr=ncoefr, &
      ncoefia=ncoefia, ncoefra=ncoefra )

!   data

    read ( unit=unit_tfemu ) coefficients%i
    read ( unit=unit_tfemu ) coefficients%r
    read ( unit=unit_tfemu ) &
      ( coefficients%ia(i)%a, i=1,ni ), &
      ( coefficients%ra(i)%a, i=1,nr )

    coefficients%created = .true.

    close ( unit=unit_tfemu )

  end subroutine read_coefficients


end module tfem_utils_m
