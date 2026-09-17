
! Copyright (C) 2004-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines routines for interfacing to sepran

module sep_utils_m

  use kind_defs_m
  use mesh_m
  use set_optional_m
  use limits_m, only: MAXPOINTS, MAXCURVES, MAXOBJECTS, MAXSURFACES, &
                      MAXVOLUMES, MAXBLOCKS, MAXNODESETS, MAXELEMENTSETS

  implicit none

  integer :: unit_sep = 13 ! unit number for reading

! interface to support older interfaces

  interface read_sepran_mesh
    module procedure read_mesh_sepran
  end interface read_sepran_mesh

contains


! read mesh from a file written by sepran (sepmesh)

  subroutine read_mesh_sepran ( mesh, filename, renumber, formatted )

    type(mesh_t), intent(inout) :: mesh

!   the filename for reading the mesh from
    character (len=*), intent(in) :: filename

!   if present and .true. and the sepran mesh contains renumbering, the
!   nodes are renumbered. Default = .false.
    logical, intent(in), optional :: renumber

!   if present and .true. the sepran mesh file format is formatted.
!   Default = .true.
    logical, intent(in), optional :: formatted

!   This routine reads basic mesh information (points, curves, topology
!   and coordinates) from an file written by sepran (sepmesh).


    character (len=11) :: lform
    character (len=80) :: current
    logical :: renum, lformatted
    integer :: elgrp, elem, curve, ios, i, j, surface, lnode, gnode
    integer :: work(10), npelm, meshtp, numnod, nelem
    integer :: ncurves, nsurfaces, nvolms, nnodes, inode, numdouble, node
    integer :: kelmt, kelmu, kelmj, kelmm, kelmn, ivers, ispec(5)

    integer, allocatable, dimension(:) :: zeropoints, gnodes, workc


    if ( mesh%meshgen ) then
      write(*,'(/a/)') 'Error in read_mesh_sepran: mesh already contains data'
      stop
    end if

    renum = set_optional ( variable=renumber, default=.false. )

    lformatted = set_optional ( variable=formatted, default=.true. )

    if ( lformatted ) then
      lform = 'formatted'
    else
      lform = 'unformatted'
    end if

!   open file

    open ( unit=unit_sep, file=filename, form=lform, iostat=ios, &
      status='old' )

    if ( ios /= 0 ) then
      write(*,'(/2a/)') 'Error in read_mesh_sepran: cannot open file ', filename
      stop
    end if

!   some scalar data

    if ( lformatted ) then
      read ( unit=unit_sep, fmt=* ) work
    else
      read ( unit=unit_sep ) work
    end if

    mesh%ndim = -work(1)
    npelm = work(2)
    meshtp = work(4)
    mesh%nelgrp = work(5)
    mesh%nnodes = work(6)
    mesh%nelem = work(7)
    mesh%npoints = work(8)
    ncurves = work(9)
    nsurfaces = work(10)

    if ( lformatted ) then
      read ( unit=unit_sep, fmt=* ) work
    else
      read ( unit=unit_sep ) work
    end if

    nvolms = work(1)
    kelmt = work(2)
    !nextpnt = work(3)
    !kelmk = work(4)
    kelmu = work(5)
    kelmj = work(6)
    kelmm = work(7)
    kelmn = work(8)
    !imap = work(9)
    ivers = work(10)

    if ( ivers < 7 ) then
      write(*,'(/a/)') 'Error in read_mesh_sepran: mesh file version too old'
      stop
    end if

!   skip three lines

    if ( lformatted ) then
      read ( unit=unit_sep, fmt=* )
      read ( unit=unit_sep, fmt=* )
      read ( unit=unit_sep, fmt=* )
    else
      read ( unit=unit_sep )
      read ( unit=unit_sep )
      read ( unit=unit_sep )
    end if

!   element and group info

    allocate( mesh%element(mesh%nelgrp) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel(mesh%nelgrp) )

    mesh%element(:)%ndim = mesh%ndim

    if ( mesh%nelgrp == 1 ) then

      mesh%element(1)%elshape = meshtp
      mesh%element(1)%numnod = npelm
      mesh%grpnumel(1) = mesh%nelem

    else

!     read info of elements from file

      if ( lformatted ) then
        read ( unit=unit_sep, fmt=* ) ! skip comment
      else
        read ( unit=unit_sep ) ! skip comment
      end if

      do elgrp = 1, mesh%nelgrp
        if ( lformatted ) then
          read ( unit=unit_sep, fmt=* ) work(1:4)
        else
          read ( unit=unit_sep ) work(1:4)
        end if
        mesh%element(elgrp)%elshape = work(2)
        mesh%element(elgrp)%numnod = work(1)
        mesh%grpnumel(elgrp) = work(3)
      end do

    end if

    if ( any ( mesh%element(:)%elshape == -10001 ) ) then
      write(*,'(/3(a/))') 'Error in read_mesh_sepran: ', &
        ' The sepran mesh contains connection elements.', &
        ' Remove meshconnect from the mesh input,'
      stop
    end if

!   set globalshape of elements

    do elgrp = 1, mesh%nelgrp
      select case ( mesh%element(elgrp)%elshape )
        case (1,2); mesh%element(elgrp)%globalshape = 'line'
        case (3,4,7,10); mesh%element(elgrp)%globalshape = 'triangle'
        case (5,6,9); mesh%element(elgrp)%globalshape = 'quadrilateral'
        case (13,14,17); mesh%element(elgrp)%globalshape = 'hexahedron'
        case (11,12,15,16,18); mesh%element(elgrp)%globalshape = 'tetrahedron'
        case default
          call errormsg_case_default ( 'read_mesh_sepran', &
            'mesh%element(elgrp)%elshape', &
            int_value=mesh%element(elgrp)%elshape )
      end select
    end do

!   coordinates

    allocate( mesh%coor(mesh%nnodes,mesh%ndim) )

    if ( lformatted ) then
      read ( unit=unit_sep, fmt=* ) ! skip comment
      read ( unit=unit_sep, fmt=* ) ( mesh%coor(node,:), node = 1, mesh%nnodes )
    else
      read ( unit=unit_sep ) ! skip comment
      read ( unit=unit_sep ) ( mesh%coor(node,:), node = 1, mesh%nnodes )
    end if

!   topology

    allocate( mesh%topology(mesh%nelgrp) )

    if ( lformatted ) then
      read ( unit=unit_sep, fmt=* ) ! skip comment
    else
      read ( unit=unit_sep ) ! skip comment
    end if

    do elgrp = 1, mesh%nelgrp

      numnod = mesh%element(elgrp)%numnod
      nelem = mesh%grpnumel(elgrp)

      allocate( mesh%topology(elgrp)%a(numnod,nelem) )

      if ( lformatted ) then
        read ( unit=unit_sep, fmt=* ) &
                       ( mesh%topology(elgrp)%a(:,elem), elem=1,nelem )
      else
        read ( unit=unit_sep ) &
                       ( mesh%topology(elgrp)%a(:,elem), elem=1,nelem )
      end if

    end do

!   points

    allocate( mesh%points(max(mesh%npoints,MAXPOINTS)) )

    if ( mesh%npoints > 0 ) then

      if ( lformatted ) then
        read ( unit=unit_sep, fmt=* ) ! skip comment
        read ( unit=unit_sep, fmt=* ) mesh%points(1:mesh%npoints)
      else
        read ( unit=unit_sep ) ! skip comment
        read ( unit=unit_sep ) mesh%points(1:mesh%npoints)
      end if

      if ( any ( mesh%points(1:mesh%npoints) == 0 ) ) then
        write(*,'(/3(a/),a)', advance='no') 'Warning in read_mesh_sepran: ', &
          ' The sepran mesh contains user points that are not connected ', &
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

    end if

!   curves

    if ( ncurves > 0 .and. &
         ( kelmt == 0 .or. nsurfaces > 0 .or. kelmm > 0 ) ) then

      mesh%ncurves = ncurves

      allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )

      if ( lformatted ) then
        read ( unit=unit_sep, fmt=* ) ! skip comment
      else
        read ( unit=unit_sep ) ! skip comment
      end if

      do curve = 1, mesh%ncurves

        if ( lformatted ) then
          read ( unit=unit_sep, fmt=* ) nnodes!, iinner
        else
          read ( unit=unit_sep ) nnodes!, iinner
        end if

        allocate( workc(nnodes) ) ! work space to store sepran nodes

        if ( lformatted ) then
          read ( unit=unit_sep, fmt='(10i8)' ) workc
        else
          read ( unit=unit_sep ) workc
        end if

!       count double nodes generated by composite curves (curves(c1,..,cn))
        numdouble = 0
        do node = 1, nnodes - 1
          if ( workc(node) == workc(node+1) ) numdouble = numdouble + 1
        end do

        mesh%curves(curve)%nnodes = nnodes - numdouble
        allocate( mesh%curves(curve)%nodes(mesh%curves(curve)%nnodes) )

!       fill nodes and remove double nodes along the way
        inode = 0
        do node = 1, nnodes
          if ( node < nnodes ) then
            if ( workc(node) == workc(node+1) ) cycle
          end if
          inode = inode + 1
          mesh%curves(curve)%nodes(inode) = workc(node)
        end do

        if ( inode /= mesh%curves(curve)%nnodes ) then
         print *, 'inode, nnodes', inode, mesh%curves(curve)%nnodes
         stop 'read_mesh_sepran: internal error'
        end if

        deallocate(workc)

      end do

!     Fill other information on curves in a heuristic way
!     If elnumnod/elshape varies for different curves it should be
!     filled manually

      mesh%curves(:mesh%ncurves)%ndim = mesh%ndim
      mesh%curves(:mesh%ncurves)%element%globalshape = 'line'
      mesh%curves(:mesh%ncurves)%element%ndim = mesh%ndim

      mesh%curves(:mesh%ncurves)%nblend = 0
      do curve = 1, mesh%ncurves
        allocate(mesh%curves(curve)%element_blend(0))
        mesh%curves(curve)%nnodes_blend = [0,mesh%curves(curve)%nnodes]
      end do

      if ( any ( mesh%element(:)%elshape == 1 ) .or. &
           any ( mesh%element(:)%elshape == 3 ) .or. &
           any ( mesh%element(:)%elshape == 5 ) .or. &
           any ( mesh%element(:)%elshape == 9 ) .or. &
           any ( mesh%element(:)%elshape == 10 ) .or. &
           any ( mesh%element(:)%elshape == 11 ) .or. &
           any ( mesh%element(:)%elshape == 13 ) .or. &
           any ( mesh%element(:)%elshape == 17 ) .or. &
           any ( mesh%element(:)%elshape == 18 ) ) then

        mesh%curves(:mesh%ncurves)%nelem = mesh%curves(:mesh%ncurves)%nnodes - 1
        mesh%curves(:mesh%ncurves)%element%elshape = 1
        mesh%curves(:mesh%ncurves)%element%numnod = 2

!       topology of elements on curves

        do curve = 1, mesh%ncurves
          allocate( mesh%curves(curve)%topology(&
                                  mesh%curves(curve)%element%numnod,&
                                  mesh%curves(curve)%nelem,2) )
          do elem = 1, mesh%curves(curve)%nelem
            mesh%curves(curve)%topology(:,elem,1) = [ (elem-1+i,i=1,2) ]
            mesh%curves(curve)%topology(:,elem,2) = &
               mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
          end do
        end do

      else if ( any ( mesh%element(:)%elshape == 2 ) .or. &
                any ( mesh%element(:)%elshape == 4 ) .or. &
                any ( mesh%element(:)%elshape == 6 ) .or. &
                any ( mesh%element(:)%elshape == 7 ) .or. &
                any ( mesh%element(:)%elshape == 12 ) .or. &
                any ( mesh%element(:)%elshape == 14 ) .or. &
                any ( mesh%element(:)%elshape == 15 ) .or. &
                any ( mesh%element(:)%elshape == 16 ) ) then

        mesh%curves(:mesh%ncurves)%nelem = &
                           (mesh%curves(:mesh%ncurves)%nnodes-1)/2
        mesh%curves(:mesh%ncurves)%element%elshape = 2
        mesh%curves(:mesh%ncurves)%element%numnod = 3

!       topology of elements on curves

        do curve = 1, mesh%ncurves
          allocate( mesh%curves(curve)%topology(&
                                           mesh%curves(curve)%element%numnod,&
                                           mesh%curves(curve)%nelem,2) )
          do elem = 1, mesh%curves(curve)%nelem
            mesh%curves(curve)%topology(:,elem,1) = [ (2*(elem-1)+i,i=1,3) ]
            mesh%curves(curve)%topology(:,elem,2) = &
               mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
          end do
        end do

      else

        write(*,'(/a/)') 'Error in read_mesh_sepran: incorrect elshape'
        stop

      end if

      if ( any ( mesh%curves(1:mesh%ncurves)%nnodes == 0 ) ) then
        write(*,'(/3(a/),a)', advance='no') 'Warning in read_mesh_sepran: ', &
          ' The sepran mesh contains curves that are not connected ', &
          ' to nodal points. Do not try to refer to these curves in any way.', &
          ' The curves numbers are: '
        allocate ( zeropoints(count(mesh%curves(1:mesh%ncurves)%nnodes == 0)) )
        j = 0
        do i = 1, mesh%ncurves
          if ( mesh%curves(i)%nnodes /= 0 ) cycle
          j = j + 1
          zeropoints(j) = i
        end do
        write (*,*) zeropoints
        write (*,*)
        deallocate ( zeropoints )
      end if

    else

      mesh%ncurves = 0

      allocate( mesh%curves(MAXCURVES) )

    end if

!   surfaces

    if ( nsurfaces > 0 .and. &
         ( kelmt == 0 .or. nvolms > 0 .or. kelmn > 0 ) ) then

      mesh%nsurfaces = nsurfaces

      allocate ( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )

      allocate ( gnodes(mesh%nnodes) ) ! work space to bookkeep global nodes
      gnodes = 0

      mesh%surfaces(:mesh%nsurfaces)%ndim = mesh%ndim
      mesh%surfaces(:mesh%nsurfaces)%element%ndim = mesh%ndim

      mesh%surfaces(:mesh%nsurfaces)%nblend = 0
      do surface = 1, mesh%nsurfaces
        allocate(mesh%surfaces(surface)%element_blend(0))
      end do

      if ( lformatted ) then
        read ( unit=unit_sep, fmt=* ) ! skip comment
      else
        read ( unit=unit_sep ) ! skip comment
      end if

      do surface = 1, mesh%nsurfaces

        if ( lformatted ) then
          read ( unit=unit_sep, fmt=* ) &
            mesh%surfaces(surface)%nelem, mesh%surfaces(surface)%element%numnod
        else
          read ( unit=unit_sep ) &
            mesh%surfaces(surface)%nelem, mesh%surfaces(surface)%element%numnod
        end if

!       element shapes on surfaces

        select case ( mesh%surfaces(surface)%element%numnod )
          case(3,6,7)
             mesh%surfaces(surface)%element%globalshape = 'triangle'
          case(4,9)
             mesh%surfaces(surface)%element%globalshape = 'quadrilateral'
          case default
            call errormsg_case_default ( 'read_mesh_sepran', &
              'mesh%surfaces(surface)%element%numnod', &
              int_value=mesh%surfaces(surface)%element%numnod )
        end select

        select case ( mesh%surfaces(surface)%element%numnod )
          case(3); mesh%surfaces(surface)%element%elshape = 3
          case(4); mesh%surfaces(surface)%element%elshape = 5
          case(6); mesh%surfaces(surface)%element%elshape = 4
          case(7); mesh%surfaces(surface)%element%elshape = 7
          case(9); mesh%surfaces(surface)%element%elshape = 6
          case default
            call errormsg_case_default ( 'read_mesh_sepran', &
              'mesh%surfaces(surface)%element%numnod', &
              int_value=mesh%surfaces(surface)%element%numnod )
        end select

        nelem = mesh%surfaces(surface)%nelem
        allocate( &
          mesh%surfaces(surface)%topology(&
                             mesh%surfaces(surface)%element%numnod,nelem,2) )

        if ( lformatted ) then
          read ( unit=unit_sep, fmt=* ) &
            ( mesh%surfaces(surface)%topology(:,elem,2), elem=1,nelem )
        else
          read ( unit=unit_sep ) &
            ( mesh%surfaces(surface)%topology(:,elem,2), elem=1,nelem )
        end if

        do elem = 1, mesh%surfaces(surface)%nelem
          gnodes ( mesh%surfaces(surface)%topology(:,elem,2) ) = 1 ! mark gnode
        end do

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
         stop 'read_mesh_sepran: internal error'
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

!   renumber

    if ( kelmj > 0 ) then

!     read renumbered nodes

      allocate ( mesh%nodperm(mesh%nnodes) )

      if ( lformatted ) then
        read ( unit=unit_sep, fmt=* ) ! skip comment
        read ( unit=unit_sep, fmt=* ) mesh%nodperm
      else
        read ( unit=unit_sep ) ! skip comment
        read ( unit=unit_sep ) mesh%nodperm
      end if

      if ( renum ) then
!       set renumber on
        mesh%renumber = .true.
      else
!       leave renumber off
        mesh%renumber = .false.
        deallocate ( mesh%nodperm )
      end if

    end if

!   spectral elements

    if ( kelmu > 0 ) then

!     skip lines until spec info

      do
        if ( lformatted ) then
          read ( unit=unit_sep, fmt='(a)', iostat=ios ) current
        else
          read ( unit=unit_sep, iostat=ios ) current
        end if
        if ( ios /= 0 ) then
          write(*,'(/a/a/)') 'Error in read_mesh_sepran: ', &
                             ' Cannot find //Spectral_Elements'
          stop
        end if
        if ( current(1:19) == '//Spectral_Elements' ) exit
      end do

!     read spectral info

      if ( lformatted ) then
        read ( unit=unit_sep, fmt=* ) ispec
      else
        read ( unit=unit_sep ) ispec
      end if

!     correct data read

!     internal mesh

      do elgrp = 1, mesh%nelgrp
        call correct_element ( mesh%element(elgrp) )
      end do

!     curves

      mesh%curves(:mesh%ncurves)%nelem = &
                  (mesh%curves(:mesh%ncurves)%nnodes - 1)/(ispec(2)+1)
      mesh%curves(:mesh%ncurves)%element%numnod = ispec(2) + 2
      do curve = 1, mesh%ncurves
        call correct_element ( mesh%curves(curve)%element )
      end do

!     topology of elements on curves

      do curve = 1, mesh%ncurves
        deallocate( mesh%curves(curve)%topology )
        allocate( mesh%curves(curve)%topology(&
               mesh%curves(curve)%element%numnod,mesh%curves(curve)%nelem,2) )
        do elem = 1, mesh%curves(curve)%nelem
          mesh%curves(curve)%topology(:,elem,1) = &
                        [ ((elem-1)*(ispec(2)+1)+i,i=1,ispec(2)+2) ]
          mesh%curves(curve)%topology(:,elem,2) = &
             mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
        end do
      end do

!     surfaces

      mesh%surfaces(:mesh%nsurfaces)%element%globalshape = 'quadrilateral'
      mesh%surfaces(:mesh%nsurfaces)%element%elshape = 102
      do surface = 1, mesh%nsurfaces
        mesh%surfaces(surface)%element%p(1:2,1) = ispec(2) + 1
        mesh%surfaces(surface)%element%p(1:2,2) = ispec(5)
      end do

    end if

!   volumes

    allocate ( mesh%volumes(MAXVOLUMES) )

!   blocks

    allocate( mesh%blocks(MAXBLOCKS) )

!   objects

    allocate( mesh%objects(MAXOBJECTS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

    mesh%meshgen = .true.

    close ( unit=unit_sep )

  contains

    subroutine correct_element ( element )

      type(element_t), intent(inout) :: element

      integer :: elshape

      elshape = element%elshape

      select case ( elshape )
        case (1)
          element%elshape = 101
          element%p(1,1) = ispec(2) + 1
          element%p(1,2) = ispec(5)
        case (5)
          element%elshape = 102
          element%p(1:2,1) = ispec(2) + 1
          element%p(1:2,2) = ispec(5)
        case (13)
          element%elshape = 103
          element%p(1:3,1) = ispec(2) + 1
          element%p(1:3,2) = ispec(5)
        case default
          stop 'internal error read_mesh_sepran; spectral elements'
      end select

    end subroutine correct_element

  end subroutine read_mesh_sepran

end module sep_utils_m
