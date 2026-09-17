
! Copyright (C) 2012-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines routines for:
! 1) reading a mesh from a file in the gmsh format (MSH ASCII or binary format).
! 2) writing a mesh (or geometries only) to a file in the gmsh format
!    (MSH ASCII or binary format).
! 3) writing mesh+data to a file or multiple files in the gmsh format
!    (MSH ASCII or binary format) for postprocessing using gmsh.
! 4) writing a scalar field to a gmsh parsed (.pos) format file for adaptive
!    meshing purposes.

! Gmsh is a 3D finite element grid generator with a build-in CAD engine
! and post-processor (http://www.geuz.org/gmsh). Gmsh is distributed under the
! terms of the GNU General Public License (GPL).

! Notes regarding reading a mesh:
!
!   * Basically, all the meshes generated using gmsh are three dimensional,
!     i.e. mesh%ndim=3. So, if one wants to create a 1D or 2D mesh the routine
!     read_mesh_gmsh needs to be called with the optional argument ndim=1 or
!     ndim=2, respectively. Then only the first or the first and second
!     coordinate of each node is used, respectively. Alternatively, it is
!     possible to convert the mesh after reading to a lower dimension using
!     mesh_convert with the argument coordinates.
!
!   * The msh file format basically consists of a set of nodal coordinates (3D)
!     and a list of elements connected to the nodes. Each element has an
!     element type and can be further characterized by tags.
!     At least two tags need to be present. The first tag is the physical entity
!     to which the element belongs, but can be zero (no physical entity
!     defined). The second tag is the elementary geometrical entity (in gmsh)
!     to which the element belongs. If there are more tags, they are simply
!     ignored. The element types can include one-node points.
!
!   * Both the ASCII and binary format can be read and is automatically
!     detected by the reader.
!
!   * The elements are converted into tfem points, curves, surfaces and element
!     groups depending on the value of domaintype (an argument of the
!     read_mesh_gmsh routine).
!     For domaintype=0, the default value, the maximum dimensionality of
!     the elements determines which elements are converted into element groups:
!      1) one-node elements are converted into tfem points.
!      2) line elements are converted into
!         a. element groups if there are no surface or volume elements
!         b. curves otherwise
!      3) triangular and quadrilateral elements are converted into
!         a. element groups if there are no volume elements
!         b. surfaces otherwise
!      4) tets, hexas, prisms and pyramids are converted into element groups
!     Using a value of domaintype > 0 will set the dimensionality (1, 2 or 3)
!     of the elements in the element groups manually. Elements having a
!     dimensionality larger than domaintype will be ignored. In summary:
!     domaintype=1:
!      1) one-node elements are converted into tfem points.
!      2) line elements are converted into element groups.
!     domaintype=2:
!      1) one-node elements are converted into tfem points.
!      2) line elements are converted into curves.
!      3) surface elements (triangles and quads) are converted into el. groups.
!     domaintype=3:
!      1) one-node elements are converted into tfem points.
!      2) line elements are converted into curves.
!      3) surface elements (triangles and quads) are converted surfaces.
!      4) volume elements (tets, hexas, prisms and pyramids) are converted
!         into el. groups.
!     NOTES:
!      - elements of dimensionality domaintype do not have to be present. For
!        example a mesh generated with gmsh -2 and read with domaintype=3
!        will have only points, curves and surfaces and no element groups.
!      - reading higher-order (3-5) elements has not been implemented.
!
!   * Points, curves, surfaces and element groups are treated as separate sets.
!     Within a set the sequence of the elements is important for the
!     number and composition of curves/surfaces/element groups generated.
!     The following rules apply for creating new curves/surfaces/element
!     groups in the sequence of elements listed within a set:
!      1) curves, surfaces and groups consist of a single element type.
!      2) curves and surfaces consist of elements belonging to the same
!           a. physical entity (in gmsh) if defined and if physgeom=.true,
!           b. geometrical entity (in gmsh), otherwise.
!      3) element groups consist of elements belonging to the same
!         physical entity, if defined.
!     Notes:
!     - Different element groups consisting of the same element type can be
!       created by introducing physical entities in gmsh.
!     - The actual values used for geometrical entities and/or physical entities
!       are ignored and only the changes are monitored to create new
!       curves/surfaces/groups in the sequence of elements within a set.
!     - Elements belonging to the same physical entity are merged into a single
!       curve, surface or group, whenever possible within the constraint given
!       by 1).
!     - If physical entities are defined in gmsh, only the physical entities
!       end up in the written msh file. So points/curves and surfaces need to
!       be put into physical entities as well in order to end up in the msh
!       file.
!
!   * The local numbering of the nodes in a curve or surface will be according
!     to the sequence of the elements as given in the gmsh file, except if
!     sortphys=.true. Therefore, the local numbering of nodes on a curve will
!     be in a natural sequence only if the elements have a natural consecutive
!     sequence in the gmsh file and sortphys=.false.
!     If not, parameters/options that depend on this natural sequence cannot be
!     used, for example:
!       step, exclude in define_essential, define_constraint,
!                        fill_sysvector, fill_vector
!
!
! Notes regarding writing a mesh + data:
!
!   * Both the ASCII and binary format is supported.
!
!   * The mesh coordinates and the data that is written to the ASCII file
!     format should not be used for computations, since these are written with
!     a single precision format to save space and meant for visualization only.
!     Full (double precision) accuracy can be obtained by using the binary
!     format.
!
!   * It is possible to write a mesh that consists of the geometries (points,
!     curves, surfaces) only. This mesh can be used to create a new internal
!     mesh using gmsh. This can, for example, be useful for remeshing in an ALE
!     computation.
!
!   * It is possible to write the mesh to a separate file first and write the
!     data to (a) separate file(s).
!
!   * You can identify different time steps using the optional arguments
!     timestepnr and/or time.
!
!
! Notes regarding writing a scalar field to a gmsh parsed (.pos) format file:
!
!   * Only scalar fields are supported.
!
!   * By default only "internal" element groups are written. It is possible to
!     write points, curves and surfaces, but these are off by default.
!
!   * Although parsed format files might be used for post-processing purposes,
!     the main application is for using the data as a "background mesh" for
!     setting target element sizes in adaptive mesh procedures.
!
!   * Prisms and pyramids are not supported.
!

module gmsh_utils_m

  use glob_defs_m
  use mesh_m
  use limits_m, only: MAXPOINTS, MAXCURVES, MAXOBJECTS, MAXSURFACES, &
                      MAXBLOCKS, MAXGROUPS, MAXVOLUMES, MAXNODESETS, &
                      MAXELEMENTSETS, WARN_ON_NO_DATA_IN_NODES, &
                      PRINT_NODES_NOT_CONNECTED, WARN_ON_NODES_NOT_CONNECTED
  use misc_m
  use array_defs_m
  use set_optional_m
  use problem_defs_m
  use system_vector_m
  use vector_m
  use streamio_m
  use, intrinsic :: iso_fortran_env

  implicit none


  integer :: unit_gmsh = 13 ! unit number for reading/writing

! type definition for a single refinement field

  type refinement_field_t

!   the coordinates of the refinement points
    real(dp), allocatable, dimension(:,:) :: coor

!   the number of refinement points
    integer :: npoints

!   the elements size is set to dx_fine within a distance distmin from the
!   refinement point, and to dx_coarse outside distance distmax. In between
!   distmin and distmax, the elementssize is interpolated
    real(dp) :: distmin, distmax, dx_fine, dx_coarse

  end type refinement_field_t

! type definition for multiple refinement fields

  type refinement_fields_t

!   number of refinement fields
    integer :: n = 0

!   spatial dimension of the refinement fields
    integer :: ndim = 0

!   the refinement fields (the i-th refinement field can be accessed by
!   refinement_fields%rf(i))
    type(refinement_field_t), allocatable, dimension(:) :: rf

  end type refinement_fields_t

contains


! read mesh from a file written by gmsh (msh file format)

  subroutine read_mesh_gmsh ( mesh, filename, ndim, physgeom, sortphys, &
    domaintype, highorder )

    type(mesh_t), intent(inout) :: mesh

!   the filename for reading the mesh from
    character (len=*), intent(in) :: filename

!   Optional dimension of space. The space dimension of gmsh is always ndim=3.
!   For ndim=1 only the first coordinate is read from the file.
!   For ndim=2, only the first and second coordinate is read from the file.
!   default is ndim=3
    integer, intent(in), optional :: ndim

!   If present and .true., physgeom indicates that the geometries (curves and
!   surfaces) are defined according to the physical entities (if these are
!   defined). Otherwise, the geometries are defined according to the geometrical
!   entities. Note, that this does not affect points.
!   default=.false.
    logical, intent(in), optional :: physgeom

!   If present and .true., sortphys sorts the gmsh elements according to the
!   physical entities (if these are defined). This ensures that the numbering of
!   physical entities in gmsh is preserved in tfem, provided that these are
!   numbered without "holes". This means, for example, that if there are
!   N physical entities for curves, all physical entities from 1 until N must be
!   defined. This must be true for all points, curves, surfaces, groups
!   separately.
!   Notes:
!     - This is only applied if also physgeom=.true.
!     - This also affects points. Give each point a different physical entity,
!       if the numbering of points has to be preserved.
!     - Elements are reshuffled due to the sorting algorithm, which affects the
!       local numbering of the nodes of geometries.
    logical, intent(in), optional :: sortphys

!   If present, domaintype determines the dimensionality of the mesh domain
!   (the element groups):
!     domaintype=0 the highest dimensional elements present in the mesh file
!                  are converted to element groups.
!     domaintype=1 the line elements are converted to element groups.
!     domaintype=2 the surface elements are converted to element groups.
!     domaintype=3 the volume elements are converted to element groups.
!   For more info see description at the start of the module.
!   Default = 0.
    integer, intent(in), optional :: domaintype

!   If present and .true., highorder indicates that all low order elements
!   (linear and quadratic) need to be imported as high order tfem elements, if
!   available.
!   Note, that this affects the nodal numbering of elements.
!   default=.false.
    logical, intent(in), optional :: highorder


!   This routine reads basic mesh information (points, curves, surfaces,
!   topology and coordinates) from a file written by gmsh (msh file format).



    character (len=80) :: current, shifted_left

!   limits build in

    integer, parameter :: MAXGMSHELEMTYPE = 124
    integer, parameter :: MAXELNODES = 1000
    integer, parameter :: MAXTAGS = 6

    logical :: binary, phys, lphysgeom, phys_enti, new, lsortphys, lhighorder
    integer(int32) :: nr, elm_type, nume, ntags, ne
    integer :: nelem, ios, nodenr, elem, ldomaintype, i, esl, tg2l, ess, &
      tg2s, esv, filetype, e, l, s, v, grp, crv, j, srf, maxdtype, dtype
    integer :: point, curve, surface, elgrp, lnode, mt, loc
    integer :: ndp_gmsh(MAXGMSHELEMTYPE), order_gmsh(MAXGMSHELEMTYPE), &
               tfem_elshape(MAXGMSHELEMTYPE), &
               typeofdomain(MAXGMSHELEMTYPE)

!    integer(int32), allocatable, dimension(:) :: elmtype
    integer, allocatable, dimension(:) :: elmtype, elshapes, gnodes, &
      ordern, gmsh_inv_index, tg1l, esl1, numel, ordere, work, tg1s, ess1, &
      tg1v, esv1, numes, numev, noelements
    integer(int32), allocatable, dimension(:,:) :: tags, nodes
    real(dp) :: coor(3)

!   number of nodal points for each element type in gmsh
!   (0 means does not exist in tfem)
    ndp_gmsh = 0
    ndp_gmsh(1:66) = [  2,  3,  4,  4,  8, 6,  5,   3,  6,   9, &
                       10, 27, 18, 14,  1, 8, 20,  15, 13,   0, &
                       10, 0,  15,  0, 21, 4,  5,   6,  0,   0, &
                        0, 0,   0,  0,  0, 16, 25, 36,  0,   0, &
                        0, 28, 36, 45, 55, 66, 49, 64, 81, 100, &
                      121, 0,   0,  0,  0,  0,  0,  0,  0,   0, &
                        0, 7, 8, 9, 10, 11 ]

!   order of each order element type in gmsh
!   (0 means does not exist in tfem)
    order_gmsh = 0
    order_gmsh(1:66) = [  1, 1, 1, 1, 1,  1, 1, 2,  2,  2, &
                          2, 2, 2, 2, 0,  2, 2, 2,  2,  0, &
                          3, 0, 4, 0, 5,  3, 4, 5,  0,  0, &
                          0, 0, 0, 0, 0,  3, 4, 5,  0,  0, &
                          0, 6, 7, 8, 9, 10, 6, 7,  8,  9, &
                         10, 0, 0, 0, 0,  0, 0, 0,  0,  0, &
                          0, 6, 7, 8, 9, 10 ]

!   tfem elshape for each element type in gmsh
!   (0: does not exist in tfem, -1: a point)
    tfem_elshape = 0
    tfem_elshape(1:66) = [  1,   3,   5,  11,  13,  41,  51,  2,    4,   6, &
                           12,  14,  43,  53,  -1,  30,  31, 42,   52,   0, &
                          104,   0, 104,   0, 104, 101, 101, 101,   0,   0, &
                            0,   0,   0,   0,   0, 102, 102, 102,   0,   0, &
                            0, 104, 104, 104, 104, 104, 102, 102, 102, 102, &
                          102,   0,   0,   0,   0,   0,   0,   0,   0,   0, &
                            0, 101, 101, 101, 101, 101 ]

!   domain type for each element type in gmsh
!   (0: does not exist in tfem, -1: a point,
!    1: line, 2: surface, 3: volume)
    typeofdomain = 0
    typeofdomain(1:66) = [  1, 2, 2, 3,  3, 3, 3, 1, 2, 2, &
                            3, 3, 3, 3, -1, 2, 3, 3, 3, 0, &
                            2, 0, 2, 0,  2, 1, 1, 1, 0, 0, &
                            0, 0, 0, 0,  0, 2, 2, 2, 0, 0, &
                            0, 2, 2, 2,  2, 2, 2, 2, 2, 2, &
                            2, 0, 0, 0,  0, 0, 0, 0, 0, 0, &
                            0, 1, 1, 1, 1, 1 ]

    if ( mesh%meshgen ) then
      write(*,'(/a/)') 'Error in read_mesh_gmsh: mesh already contains data'
      stop
    end if

!   set space dimension

    mesh%ndim = set_optional ( variable=ndim, default=3 )
    if ( .not. any( mesh%ndim == [ 1, 2, 3 ] ) ) then
      write(*,'(/a,i0/)') 'Incorrect space dimension = ', mesh%ndim
    end if

!   set physgeom

    lphysgeom = set_optional ( variable=physgeom, default=.false. )

!   set sortphys

    lsortphys = set_optional ( variable=sortphys, default=.false. )

!   set domaintype

    ldomaintype = set_optional ( variable=domaintype, default=0 )

    if ( .not. any( ldomaintype == [ 0, 1, 2, 3 ] ) ) then
      write(*,'(/a,i0/)') 'Incorrect domaintype = ', ldomaintype
    end if

!   set highorder

    lhighorder = set_optional ( variable=highorder, default=.false. )

    if ( lhighorder ) then

!     import low order elements as high order tfem elements

      tfem_elshape([1,8])  = 101  ! lines
      tfem_elshape([2,9])  = 104  ! triangles
      tfem_elshape([3,10]) = 102  ! quadrilaterals

    end if

!   open file

    open ( unit=unit_gmsh, file=filename, form='formatted', iostat=ios, &
      status='old' )

    if ( ios /= 0 ) then
      write(*,'(/2a/)') 'Error in read_mesh_gmsh: cannot open file ', filename
      stop
    end if


    binary = .false.

!   skip lines until $MeshFormat is found

    call skip_lines ( '$MeshFormat' )

    read ( unit=unit_gmsh, fmt='(a)' ) current

!   read version
    shifted_left = adjustl ( current )
    if ( shifted_left(1:3) /= '2.2' ) then
      if ( shifted_left(1:1) == '4' ) then
        write(*,'(/2(a/),a)') &
         'You are using gmsh version 4.0.0 or higher, which introduces', &
         'a new msh format version. You can continue using gmsh v4 after ', &
         'applying one of the fixes in the following error message.'
      end if
      write(*,'(/6(a/),2a/)') 'Error in read_mesh_gmsh:', &
        'File version not supported.', &
        'Only file version 2.2 is currently supported.', &
        'To fix this, do one of the following:', &
        ' - run with ''gmsh -format msh2'', ', &
        ' - put a line ''Mesh.MshFileVersion = 2.2;'' in the .geo file', &
        ' - put a line ''Mesh.MshFileVersion = 2.2;'' in the ', &
        '$HOME/.gmsh-options file.'
      stop
    end if

!   ASCII or binary?
    read ( unit=shifted_left(4:), fmt=* ) filetype
    binary = filetype == 1

    if ( binary ) then
!     reopen as unformatted stream
      close ( unit=unit_gmsh )
      open ( unit=unit_gmsh, file=filename, access='stream', iostat=ios, &
        form='unformatted', status='old' )
    end if


!   skip lines until $Nodes is found

    call skip_lines ( '$Nodes' )

    if ( binary ) then
      call read_line ( unit=unit_gmsh, line=current )
      read ( unit=current, fmt=* ) mesh%nnodes
    else
      read ( unit=unit_gmsh, fmt=* ) mesh%nnodes
    end if

    allocate( mesh%coor(mesh%nnodes,mesh%ndim) )

!   read nodes

    do nodenr = 1, mesh%nnodes

      if ( binary ) then
        read ( unit=unit_gmsh ) nr, mesh%coor(nodenr,:), coor(mesh%ndim+1:3)
      else
        read ( unit=unit_gmsh, fmt=* ) nr, mesh%coor(nodenr,:)
      end if

      if ( nr /= nodenr ) then
        write(*,'(2(/a)/2(a,i0))') 'Error in read_mesh_gmsh:', &
          ' Node numbers not in sequence start:1 increment:1. Renumber first', &
          ' nodenr = ', nodenr, ' gmsh node = ', nr
        stop
      end if

    end do


!   skip lines until $Elements is found

    call skip_lines ( '$Elements' )

    if ( binary ) then
      call read_line ( unit=unit_gmsh, line=current )
      read ( unit=current, fmt=* ) nelem
    else
      read ( unit=unit_gmsh, fmt=* ) nelem
    end if

!   allocate temporary storage for element connectivity

    allocate ( elmtype(nelem), tags(nelem,MAXTAGS), nodes(nelem,MAXELNODES) )

!   read element data

    if ( binary ) then

!     binary mesh file

      elem = 0

      do

        read ( unit=unit_gmsh ) elm_type, nume, ntags

        if ( ntags < 2 ) then
          write(*,'(2(/a)/)') 'Error in read_mesh_gmsh:', &
            ' Number of element tags must be at least two.'
          stop
        end if

        do e = 1, nume

          elem = elem + 1

          elmtype(elem) = elm_type

          read ( unit=unit_gmsh ) &
            ne, tags(elem,1:ntags), nodes(elem,1:ndp_gmsh(elmtype(elem)))

          if ( ndp_gmsh(elmtype(elem)) == 0 ) then
            write(*,'(/a/a,i0,a)') 'Error in read_mesh_gmsh:', &
              ' Gmsh element type ', elmtype(elem), ' not available in tfem.'
            stop

          end if

        end do

        if ( elem == nelem ) exit

      end do

    else

!     ASCII mesh file

      do elem = 1, nelem

        read ( unit=unit_gmsh, fmt=* ) ne, elmtype(elem), ntags, &
          tags(elem,1:ntags), nodes(elem,1:ndp_gmsh(elmtype(elem)))

        if ( ndp_gmsh(elmtype(elem)) == 0 ) then
          write(*,'(/a/a,i0,a)') 'Error in read_mesh_gmsh:', &
            ' Gmsh element type ', elmtype(elem), ' not available in tfem.'
          stop
        end if
        if ( ntags < 2 ) then
          write(*,'(2(/a)/)') 'Error in read_mesh_gmsh:', &
            ' Number of element tags must be at least two.'
          stop
        end if

      end do

    end if

!   check tags for physical entities

    if ( any( tags(:,1) > 0 ) ) then
!     physical entities
      if ( any( tags(:,1) == 0 ) ) then
        write(*,'(2(/a)/)') 'Error in read_mesh_gmsh:', &
          ' Physical entities defined but some have a zero number.'
        stop
      end if
      phys_enti = .true.
      mt = maxval(tags(:,1))
    else
!     no physical entities
      phys_enti = .false.
      mt = 0
      if ( lphysgeom ) then
        write(*,'(2(/a)/)') 'Warning in read_mesh_gmsh:', &
          ' Heading argument physgeom=.true and no physical entities defined.'
      end if
    end if

    phys = lphysgeom .and. phys_enti  ! physical entity for geometry

    allocate ( elshapes(nelem) )

!   determine elshape of each element

    elshapes = tfem_elshape ( elmtype )

    if ( any ( elshapes == 0 ) ) stop 'internal error 1'

!   determine the elements of the "highest" geometry
!     maxdtype=-1 point
!     maxdtype=1 line
!     maxdtype=2 surface
!     maxdtype=3 volume

    maxdtype = maxval ( typeofdomain ( elmtype ) )

    if ( maxdtype == -1 ) then
      write(*,'(2(/a)/)') 'Warning in read_mesh_gmsh:', &
        ' Mesh consists of points only.'
    end if

!   select domaintype

    if ( ldomaintype == 0 ) then

      dtype = maxdtype

    else

      dtype = ldomaintype

    end if

!   order of elements

    allocate ( ordere(nelem) )

    if ( lsortphys .and. phys ) then
!     sort elements according to physical entities
      allocate ( work(nelem) )
      work = tags(:,1)
      call sort ( work, ordere )
      deallocate ( work )
    else
!     standard order
      do i = 1, nelem
        ordere(i) = i
      end do
      if ( lsortphys .and. .not. phys_enti ) then
        write(*,'(3(/a)/)') 'Warning in read_mesh_gmsh:', &
          ' Heading argument sortphys=.true and no physical entities defined.',&
          ' sortphys ignored.'
      end if
      if ( lsortphys .and. .not. lphysgeom ) then
        write(*,'(3(/a)/)') 'Warning in read_mesh_gmsh:', &
          ' Heading argument sortphys=.true and physgeom=.false.', &
          ' sortphys ignored.'
      end if
    end if

    allocate ( tg1l(mt), esl1(mt), numel(mt) )
    allocate ( tg1s(mt), ess1(mt), numes(mt) )
    allocate ( tg1v(mt), esv1(mt), numev(mt) )

!   count number of points, curves, surfaces, groups

    mesh%npoints = 0
    mesh%ncurves = 0
    mesh%nsurfaces = 0
    mesh%nelgrp = 0

    esl = 0; l = 0; tg2l = 0
    ess = 0; s = 0; tg2s = 0
    esv = 0; v = 0

    select case ( dtype )
      case(-1) ! only points
        mesh%npoints = nelem
      case(1) ! points and line mesh
        do j = 1, nelem ! loop over all entities
          i = ordere(j)
          if ( elshapes(i) == -1 ) then
            mesh%npoints = mesh%npoints + 1
          else if ( any ( elshapes(i) == [ 1, 2, 101 ] ) ) then
            if ( phys_enti ) then ! physical entities for groups
              if ( any ( tg1l(1:l) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1l(1:l)-tags(i,1))))
                if ( elshapes(i) == esl1(loc) ) cycle  ! same group
              end if
              l=l+1; tg1l(l) = tags(i,1); esl1(l) = elshapes(i)
            else
              if ( elshapes(i) == esl ) cycle ! same group
              esl = elshapes(i)
            end if
            mesh%nelgrp = mesh%nelgrp + 1
          end if
        end do
      case(2) ! points, curves and surface mesh
        do j = 1, nelem ! loop over all entities
          i = ordere(j)
          if ( elshapes(i) == -1 ) then
            mesh%npoints = mesh%npoints + 1
          else if ( any ( elshapes(i) == [ 1, 2, 101 ] ) ) then
            if ( phys ) then ! physical entities for geometries
              if ( any ( tg1l(1:l) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1l(1:l)-tags(i,1))))
                if ( elshapes(i) == esl1(loc) ) cycle  ! same geometry
              end if
              l=l+1; tg1l(l) = tags(i,1); esl1(l) = elshapes(i)
            else
              if ( elshapes(i) == esl .and. &
                   tags(i,2) == tg2l ) cycle ! same geometry
              esl = elshapes(i); tg2l = tags(i,2)
            end if
            mesh%ncurves = mesh%ncurves + 1
          else if ( any ( elshapes(i) == [ 3, 4, 5, 6, 30, 102, 104 ] ) ) then
            if ( phys_enti ) then ! physical entities for groups
              if ( any ( tg1s(1:s) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1s(1:s)-tags(i,1))))
                if ( elshapes(i) == ess1(loc) ) cycle  ! same group
              end if
              s=s+1; tg1s(s) = tags(i,1); ess1(s) = elshapes(i)
            else
              if ( elshapes(i) == ess ) cycle ! same group
              ess = elshapes(i)
            end if
            mesh%nelgrp = mesh%nelgrp + 1
          end if
        end do
      case(3) ! points, curves, surfaces and volume mesh
        do j = 1, nelem ! loop over all entities
          i = ordere(j)
          phys = lphysgeom .and. tags(i,1) > 0  ! physical entity for geometry
          if ( elshapes(i) == -1 ) then
            mesh%npoints = mesh%npoints + 1
          else if ( any ( elshapes(i) == [ 1, 2, 101 ] ) ) then
            if ( phys ) then ! physical entities for geometries
              if ( any ( tg1l(1:l) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1l(1:l)-tags(i,1))))
                if ( elshapes(i) == esl1(loc) ) cycle  ! same geometry
              end if
              l=l+1; tg1l(l) = tags(i,1); esl1(l) = elshapes(i)
            else
              if ( elshapes(i) == esl .and. &
                   tags(i,2) == tg2l ) cycle ! same geometry
              esl = elshapes(i); tg2l = tags(i,2)
            end if
            mesh%ncurves = mesh%ncurves + 1
          else if ( any ( elshapes(i) == [ 3, 4, 5, 6, 30, 102, 104 ] ) ) then
            if ( phys ) then ! physical entities for geometries
              if ( any ( tg1s(1:s) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1s(1:s)-tags(i,1))))
                if ( elshapes(i) == ess1(loc) ) cycle  ! same geometry
              end if
              s=s+1; tg1s(s) = tags(i,1); ess1(s) = elshapes(i)
            else
              if ( elshapes(i) == ess .and. &
                   tags(i,2) == tg2s ) cycle ! same geometry
              ess = elshapes(i); tg2s = tags(i,2)
            end if
            mesh%nsurfaces = mesh%nsurfaces + 1
          else if ( any ( elshapes(i) == &
                      [ 11, 12, 13, 14, 31, 41, 42, 43, 51, 52, 53 ] ) ) then
            if ( phys_enti ) then ! physical entities for groups
              if ( any ( tg1v(1:v) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1v(1:v)-tags(i,1))))
                if ( elshapes(i) == esv1(loc) ) cycle  ! same group
              end if
              v=v+1; tg1v(v) = tags(i,1); esv1(v) = elshapes(i)
            else
              if ( elshapes(i) == esv ) cycle ! same group
              esv = elshapes(i)
            end if
            mesh%nelgrp = mesh%nelgrp + 1
          end if
        end do
      case default
        call errormsg_case_default ( 'read_mesh_gmsh', &
          'dtype', int_value=dtype )
    end select

!   allocate data for points, curves, surfaces and element groups

    allocate( mesh%points(max(mesh%npoints,MAXPOINTS)) )
    allocate( mesh%curves(max(mesh%ncurves,MAXCURVES)) )
    allocate( mesh%surfaces(max(mesh%nsurfaces,MAXSURFACES)) )
    allocate( mesh%element(mesh%nelgrp) )
    allocate( mesh%element_blend(mesh%nelgrp,0) )
    mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel(mesh%nelgrp) )

!   determine element types, number of elements and optionally polynomial order

    curve = 0
    surface = 0
    elgrp = 0
    mesh%grpnumel = 0

    esl = 0; l = 0; tg2l = 0
    ess = 0; s = 0; tg2s = 0
    esv = 0; v = 0

    select case ( dtype )
      case(1) ! points and line mesh
        do j = 1, nelem ! loop over all entities
          i = ordere(j)
          if ( any ( elshapes(i) == [ 1, 2, 101 ] ) ) then
            if ( phys_enti ) then ! physical entities for groups
              if ( any ( tg1l(1:l) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1l(1:l)-tags(i,1))))
                if ( elshapes(i) == esl1(loc) ) then
!                 same group
                  mesh%grpnumel(loc) = mesh%grpnumel(loc) + 1
                  cycle
                end if
              end if
              l = l + 1
            else
              if ( elshapes(i) == esl ) then
!               same group
                mesh%grpnumel(elgrp) = mesh%grpnumel(elgrp) + 1
                cycle
              end if
              esl = elshapes(i)
            end if
            elgrp = elgrp + 1
            mesh%element(elgrp)%elshape = elshapes(i)
            mesh%element(elgrp)%numnod = ndp_gmsh(elmtype(i))
            if ( elshapes(i) == 101 ) then
              mesh%element(elgrp)%p(1,1) = order_gmsh(elmtype(i))
            end if
            mesh%grpnumel(elgrp) = 1
          end if
        end do
      case(2) ! points, curves and surface mesh
        do j = 1, nelem ! loop over all entities
          i = ordere(j)
          if ( any ( elshapes(i) == [ 1, 2, 101 ] ) ) then
            if ( phys ) then ! physical entities for geometries
              if ( any ( tg1l(1:l) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1l(1:l)-tags(i,1))))
                if ( elshapes(i) == esl1(loc) ) then
!                 same geometry
                  mesh%curves(loc)%nelem = mesh%curves(loc)%nelem + 1
                  cycle
                end if
              end if
              l = l + 1
            else
              if ( elshapes(i) == esl .and. &
                   tags(i,2) == tg2l ) then
!               same geometry
                mesh%curves(curve)%nelem = mesh%curves(curve)%nelem + 1
                cycle
              end if
              esl = elshapes(i); tg2l = tags(i,2)
            end if
            curve = curve + 1
            mesh%curves(curve)%element%elshape = elshapes(i)
            mesh%curves(curve)%element%numnod = ndp_gmsh(elmtype(i))
            if ( elshapes(i) == 101 ) then
              mesh%curves(curve)%element%p(1,1) = order_gmsh(elmtype(i))
            end if
            mesh%curves(curve)%nelem = 1
          else if ( any ( elshapes(i) == [ 3, 4, 5, 6, 30, 102, 104 ] ) ) then
            if ( phys_enti ) then ! physical entities for groups
              if ( any ( tg1s(1:s) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1s(1:s)-tags(i,1))))
                if ( elshapes(i) == ess1(loc) ) then
!                 same group
                  mesh%grpnumel(loc) = mesh%grpnumel(loc) + 1
                  cycle
                end if
              end if
              s = s + 1
            else
              if ( elshapes(i) == ess ) then
!               same group
                mesh%grpnumel(elgrp) = mesh%grpnumel(elgrp) + 1
                cycle
              end if
              ess = elshapes(i)
            end if
            elgrp = elgrp + 1
            mesh%element(elgrp)%elshape = elshapes(i)
            mesh%element(elgrp)%numnod = ndp_gmsh(elmtype(i))
            if ( any( elshapes(i) == [ 102, 104 ] ) ) then
              mesh%element(elgrp)%p(1:2,1) = order_gmsh(elmtype(i))
            end if
            mesh%grpnumel(elgrp) = 1
          end if
        end do
      case(3) ! points, curves, surfaces and volume mesh
        do j = 1, nelem ! loop over all entities
          i = ordere(j)
          phys = lphysgeom .and. tags(i,1) > 0  ! physical entity for geometry
          if ( any ( elshapes(i) == [ 1, 2, 101 ] ) ) then
            if ( phys ) then ! physical entities for geometries
              if ( any ( tg1l(1:l) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1l(1:l)-tags(i,1))))
                if ( elshapes(i) == esl1(loc) ) then
!                 same geometry
                  mesh%curves(loc)%nelem = mesh%curves(loc)%nelem + 1
                  cycle
                end if
              end if
              l = l + 1
            else
              if ( elshapes(i) == esl .and. &
                   tags(i,2) == tg2l ) then
!               same geometry
                mesh%curves(curve)%nelem = mesh%curves(curve)%nelem + 1
                cycle
              end if
              esl = elshapes(i); tg2l = tags(i,2)
            end if
            curve = curve + 1
            mesh%curves(curve)%element%elshape = elshapes(i)
            mesh%curves(curve)%element%numnod = ndp_gmsh(elmtype(i))
            mesh%curves(curve)%nelem = 1
            if ( elshapes(i) == 101 ) then
              mesh%curves(curve)%element%p(1,1) = order_gmsh(elmtype(i))
            end if
          else if ( any ( elshapes(i) == [ 3, 4, 5, 6, 30, 102, 104 ] ) ) then
            if ( phys ) then ! physical entities for geometries
              if ( any ( tg1s(1:s) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1s(1:s)-tags(i,1))))
                if ( elshapes(i) == ess1(loc) ) then
!                 same geometry
                  mesh%surfaces(loc)%nelem = mesh%surfaces(loc)%nelem + 1
                  cycle
                end if
              end if
              s = s + 1
            else
              if ( elshapes(i) == ess .and. &
                   tags(i,2) == tg2s ) then
!               same geometry
                mesh%surfaces(surface)%nelem = mesh%surfaces(surface)%nelem + 1
                cycle
              end if
              ess = elshapes(i); tg2s = tags(i,2)
            end if
            surface = surface + 1
            mesh%surfaces(surface)%element%elshape = elshapes(i)
            mesh%surfaces(surface)%element%numnod = ndp_gmsh(elmtype(i))
            if ( any ( elshapes(i) == [ 102, 104 ] ) ) then
              mesh%surfaces(surface)%element%p(1:2,1) = order_gmsh(elmtype(i))
            end if
            mesh%surfaces(surface)%nelem = 1
          else if ( any ( elshapes(i) == &
                      [ 11, 12, 13, 14, 31, 41, 42, 43, 51, 52, 53 ] ) ) then
            if ( phys_enti ) then ! physical entities for groups
              if ( any ( tg1v(1:v) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1v(1:v)-tags(i,1))))
                if ( elshapes(i) == esv1(loc) ) then
!                 same group
                  mesh%grpnumel(loc) = mesh%grpnumel(loc) + 1
                  cycle
                end if
              end if
              v = v + 1
            else
              if ( elshapes(i) == esv ) then
!               same group
                mesh%grpnumel(elgrp) = mesh%grpnumel(elgrp) + 1
                cycle
              end if
              esv = elshapes(i)
            end if
            elgrp = elgrp + 1
            mesh%element(elgrp)%elshape = elshapes(i)
            mesh%element(elgrp)%numnod = ndp_gmsh(elmtype(i))
            mesh%grpnumel(elgrp) = 1
          end if
        end do
      case default
        call errormsg_case_default ( 'read_mesh_gmsh', 'dtype', &
          int_value=dtype )
    end select

!   set some more data

    mesh%curves(:mesh%ncurves)%ndim  = mesh%ndim
    mesh%curves(:mesh%ncurves)%element%ndim  = mesh%ndim
    mesh%surfaces(:mesh%nsurfaces)%ndim  = mesh%ndim
    mesh%surfaces(:mesh%nsurfaces)%element%ndim  = mesh%ndim
    mesh%nelem = sum ( mesh%grpnumel )
    mesh%element(:)%ndim  = mesh%ndim

    mesh%curves(:mesh%ncurves)%nblend = 0
    do curve = 1, mesh%ncurves
      allocate(mesh%curves(curve)%element_blend(0))
    end do
    mesh%surfaces(:mesh%nsurfaces)%nblend = 0
    do surface = 1, mesh%nsurfaces
      allocate(mesh%surfaces(surface)%element_blend(0))
    end do

!   set globalshape of elements

    mesh%curves(:mesh%ncurves)%element%globalshape = 'line'

    do surface = 1, mesh%nsurfaces
      select case ( mesh%surfaces(surface)%element%elshape )
        case (3,4,104); mesh%surfaces(surface)%element%globalshape = 'triangle'
        case (5,6,30,102)
                    mesh%surfaces(surface)%element%globalshape = 'quadrilateral'
        case default
          call errormsg_case_default ( 'read_mesh_gmsh', &
            'mesh%surfaces(surface)%element%elshape', &
            int_value=mesh%surfaces(surface)%element%elshape )
      end select
    end do

    do elgrp = 1, mesh%nelgrp
      select case ( mesh%element(elgrp)%elshape )
        case (1,2,101); mesh%element(elgrp)%globalshape = 'line'
        case (3,4,104); mesh%element(elgrp)%globalshape = 'triangle'
        case (5,6,30,102); mesh%element(elgrp)%globalshape = 'quadrilateral'
        case (13,14,31); mesh%element(elgrp)%globalshape = 'hexahedron'
        case (11,12); mesh%element(elgrp)%globalshape = 'tetrahedron'
        case (41,42,43); mesh%element(elgrp)%globalshape = 'prism'
        case (51,52,53); mesh%element(elgrp)%globalshape = 'pyramid'
        case default
          call errormsg_case_default ( 'read_mesh_gmsh', &
            'mesh%element(elgrp)%elshape', &
            int_value=mesh%element(elgrp)%elshape )
      end select
    end do


!   allocate topology arrays

    do curve = 1, mesh%ncurves
      allocate( mesh%curves(curve)%topology(&
                    mesh%curves(curve)%element%numnod,&
                    mesh%curves(curve)%nelem,2) )
    end do

    do surface = 1, mesh%nsurfaces
      allocate( mesh%surfaces(surface)%topology(&
                    mesh%surfaces(surface)%element%numnod,&
                    mesh%surfaces(surface)%nelem,2) )
    end do

    allocate( mesh%topology(mesh%nelgrp) )

    do elgrp = 1, mesh%nelgrp
      allocate( mesh%topology(elgrp)%a(mesh%element(elgrp)%numnod,&
                                       mesh%grpnumel(elgrp)) )
    end do


!   fill points and topology of elements

    point = 0
    curve = 0
    surface = 0
    elgrp = 0

    esl = 0; l = 0; tg2l = 0
    ess = 0; s = 0; tg2s = 0
    esv = 0; v = 0

    numel = 0; numes = 0; numev = 0

    select case ( dtype )
      case(-1) ! only points
        do i = 1, nelem ! loop over all entities
          point = point + 1
          mesh%points(point) = nodes(i,1)
        end do
      case(1) ! points and line mesh
        do j = 1, nelem ! loop over all entities
          i = ordere(j)
          if ( elshapes(i) == -1 ) then
            point = point + 1
            mesh%points(point) = nodes(i,1)
          else if ( any ( elshapes(i) == [ 1, 2, 101 ] ) ) then
            if ( phys_enti ) then ! physical entities for groups
              new = .true.
              if ( any ( tg1l(1:l) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1l(1:l)-tags(i,1))))
                if ( elshapes(i) == esl1(loc) ) then
!                 same group
                  numel(loc) = numel(loc) + 1
                  elem = numel(loc); grp = loc; new = .false.
                end if
              end if
              if ( new ) then ! new group: increment counters, reset elem to 1
                l = l + 1; elgrp = elgrp + 1; numel(elgrp) = 1
                elem = 1; grp = elgrp
              end if
            else
              if ( elshapes(i) == esl ) then
!               same group
                elem = elem + 1
              else
                elgrp = elgrp + 1
                elem = 1; grp = elgrp
                esl = elshapes(i)
              end if
            end if
            mesh%topology(grp)%a(:,elem) = &
                                    nodes(i,1:mesh%element(grp)%numnod)
          end if
        end do
      case(2) ! points, curves and surface mesh
        do j = 1, nelem ! loop over all entities
          i = ordere(j)
          if ( elshapes(i) == -1 ) then
            point = point + 1
            mesh%points(point) = nodes(i,1)
          else if ( any ( elshapes(i) == [ 1, 2, 101 ] ) ) then
            if ( phys ) then ! physical entities for geometries
              new = .true.
              if ( any ( tg1l(1:l) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1l(1:l)-tags(i,1))))
                if ( elshapes(i) == esl1(loc) ) then
!                 same geometry
                  numel(loc) = numel(loc) + 1
                  elem = numel(loc); crv = loc; new = .false.
                end if
              end if
              if ( new ) then ! new geom: increment counters, reset elem to 1
                l = l + 1; curve = curve + 1; numel(curve) = 1
                elem = 1; crv = curve
              end if
            else
              if ( elshapes(i) == esl .and. &
                   tags(i,2) == tg2l ) then
!               same geometry
                elem = elem + 1
              else
                curve = curve + 1
                elem = 1; crv = curve
                esl = elshapes(i); tg2l = tags(i,2)
              end if
            end if
            mesh%curves(crv)%topology(:,elem,2) = &
                                    nodes(i,1:mesh%curves(crv)%element%numnod)
          else if ( any ( elshapes(i) == [ 3, 4, 5, 6, 30, 102, 104 ] ) ) then
            if ( phys_enti ) then ! physical entities for groups
              new = .true.
              if ( any ( tg1s(1:s) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1s(1:s)-tags(i,1))))
                if ( elshapes(i) == ess1(loc) ) then
!                 same group
                  numes(loc) = numes(loc) + 1
                  elem = numes(loc); grp = loc; new = .false.
                end if
              end if
              if ( new ) then ! new group: increment counters, reset elem to 1
                s = s + 1; elgrp = elgrp + 1; numes(elgrp) = 1
                elem = 1; grp = elgrp
              end if
            else
              if ( elshapes(i) == ess ) then
!               same group
                elem = elem + 1
              else
                elgrp = elgrp + 1
                elem = 1; grp = elgrp
                ess = elshapes(i)
              end if
            end if
            mesh%topology(grp)%a(:,elem) = &
                                    nodes(i,1:mesh%element(grp)%numnod)
          end if
        end do
      case(3) ! points, curves, surfaces and volume mesh
        do j = 1, nelem ! loop over all entities
          i = ordere(j)
          phys = lphysgeom .and. tags(i,1) > 0  ! physical entity for geometry
          if ( elshapes(i) == -1 ) then
            point = point + 1
            mesh%points(point) = nodes(i,1)
          else if ( any ( elshapes(i) == [ 1, 2, 101 ] ) ) then
            if ( phys ) then ! physical entities for geometries
              new = .true.
              if ( any ( tg1l(1:l) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1l(1:l)-tags(i,1))))
                if ( elshapes(i) == esl1(loc) ) then
!                 same geometry
                  numel(loc) = numel(loc) + 1
                  elem = numel(loc); crv = loc; new = .false.
                end if
              end if
              if ( new ) then ! new geom: increment counters, reset elem to 1
                l = l + 1; curve = curve + 1; numel(curve) = 1
                elem = 1; crv = curve
              end if
            else
              if ( elshapes(i) == esl .and. &
                   tags(i,2) == tg2l ) then
!               same geometry
                elem = elem + 1
              else
                curve = curve + 1
                elem = 1; crv = curve
                esl = elshapes(i); tg2l = tags(i,2)
              end if
            end if
            mesh%curves(crv)%topology(:,elem,2) = &
                                    nodes(i,1:mesh%curves(crv)%element%numnod)
          else if ( any ( elshapes(i) == [ 3, 4, 5, 6, 30, 102, 104 ] ) ) then
            if ( phys ) then ! physical entities for geometries
              new = .true.
              if ( any ( tg1s(1:s) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1s(1:s)-tags(i,1))))
                if ( elshapes(i) == ess1(loc) ) then
!                 same geometry
                  numes(loc) = numes(loc) + 1
                  elem = numes(loc); srf = loc; new = .false.
                end if
              end if
              if ( new ) then ! new geom: increment counters, reset elem to 1
                s = s + 1; surface = surface + 1; numes(surface) = 1
                elem = 1; srf = surface
              end if
            else
              if ( elshapes(i) == ess .and. &
                   tags(i,2) == tg2s ) then
!               same geometry
                elem = elem + 1
              else
                surface = surface + 1
                elem = 1; srf = surface
                ess = elshapes(i); tg2s = tags(i,2)
              end if
            end if
            mesh%surfaces(srf)%topology(:,elem,2) = &
                                  nodes(i,1:mesh%surfaces(srf)%element%numnod)
          else if ( any ( elshapes(i) == &
                      [ 11, 12, 13, 14, 31, 41, 42, 43, 51, 52, 53 ] ) ) then
            if ( phys_enti ) then ! physical entities for groups
              new = .true.
              if ( any ( tg1v(1:v) == tags(i,1) ) ) then
                loc = maxval(minloc(abs(tg1v(1:v)-tags(i,1))))
                if ( elshapes(i) == esv1(loc) ) then
!                 same group
                  numev(loc) = numev(loc) + 1
                  elem = numev(loc); grp = loc; new = .false.
                end if
              end if
              if ( new ) then ! new group: increment counters, reset elem to 1
                v = v + 1; elgrp = elgrp + 1; numev(elgrp) = 1
                elem = 1; grp = elgrp
              end if
            else
              if ( elshapes(i) == esv ) then
!               same group
                elem = elem + 1
              else
                elgrp = elgrp + 1
                elem = 1; grp = elgrp
                esv = elshapes(i)
              end if
            end if
            mesh%topology(elgrp)%a(:,elem) = &
                                    nodes(i,1:mesh%element(elgrp)%numnod)
          end if
        end do
      case default
        call errormsg_case_default ( 'read_mesh_gmsh', &
          'dtype', int_value=dtype )
    end select


!   correct topology arrays for sequence of the nodes

    do curve = 1, mesh%ncurves

      allocate ( gmsh_inv_index(mesh%curves(curve)%element%numnod) )

      call fill_gmsh_inv_index ( mesh%curves(curve)%element%elshape, &
                                 mesh%curves(curve)%element%p(1,1), &
                                 gmsh_inv_index )

      do elem = 1, mesh%curves(curve)%nelem

        mesh%curves(curve)%topology(:,elem,2) = &
                mesh%curves(curve)%topology(gmsh_inv_index,elem,2)
      end do

      deallocate ( gmsh_inv_index )

    end do

    do surface = 1, mesh%nsurfaces

      allocate ( gmsh_inv_index(mesh%surfaces(surface)%element%numnod) )

      call fill_gmsh_inv_index ( mesh%surfaces(surface)%element%elshape, &
                                 mesh%surfaces(surface)%element%p(1,1), &
                                 gmsh_inv_index )

      do elem = 1, mesh%surfaces(surface)%nelem
        mesh%surfaces(surface)%topology(:,elem,2) = &
                mesh%surfaces(surface)%topology(gmsh_inv_index,elem,2)
      end do

      deallocate ( gmsh_inv_index )

    end do

    do elgrp = 1, mesh%nelgrp

      allocate ( gmsh_inv_index(mesh%element(elgrp)%numnod) )

      call fill_gmsh_inv_index ( mesh%element(elgrp)%elshape, &
                                 mesh%element(elgrp)%p(1,1), &
                                 gmsh_inv_index )

      do elem = 1, mesh%grpnumel(elgrp)
        mesh%topology(elgrp)%a(:,elem) = &
                          mesh%topology(elgrp)%a(gmsh_inv_index,elem)
      end do

      deallocate ( gmsh_inv_index )

    end do


!   work space to bookkeep global nodes

    allocate ( gnodes(mesh%nnodes), ordern(mesh%nnodes) )

    gnodes = 0

!   fill local nodes and topology in curves and surfaces

    do curve = 1, mesh%ncurves
      call fill_local_geometry ( mesh%curves(curve) )
    end do

    do surface = 1, mesh%nsurfaces
      call fill_local_geometry ( mesh%surfaces(surface) )
    end do

!   test nodes

    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)
        gnodes(mesh%topology(elgrp)%a(:,elem)) = 1
      end do
    end do

    if ( WARN_ON_NODES_NOT_CONNECTED .and. any( gnodes == 0 ) ) then

      if ( PRINT_NODES_NOT_CONNECTED ) then

        write(*,'(/a/a)') &
          'Warning in read_mesh_gmsh: ', &
          '  Nodes not connected to the internal elements:'

        allocate ( noelements(count(gnodes==0)) )

        j = 0
        do i = 1, mesh%nnodes
          if ( gnodes(i) /= 0 ) cycle
          j = j + 1
          noelements(j) = i
        end do

        write (*,*) noelements
        write (*,*)

        deallocate ( noelements )

      else

        write(*,'(/a/a,i0,a/)') &
          'Warning in read_mesh_gmsh: ', &
          '  There are ', count(gnodes==0), &
          ' nodes not connected to internal elements.'

      end if

    end if

    deallocate ( gnodes, ordern )


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

    close ( unit=unit_gmsh )

!   remove memory

    deallocate ( elshapes )
    deallocate ( elmtype, tags, nodes )
    deallocate ( ordere )

    deallocate ( tg1l, numel )
    deallocate ( tg1s, numes )
    deallocate ( tg1v, numev )


  contains


!   skip lines until key is found

    subroutine skip_lines ( key )

      character(len=*), intent(in) :: key

      do

        if ( binary ) then

          call read_line ( unit=unit_gmsh, line=current )

        else

          read ( unit=unit_gmsh, fmt='(a)', iostat=ios ) current
          if ( ios /= 0 ) then
            write(*,'(/2a/)') 'Error in read_mesh_gmsh: cannot find ', key
            stop
          end if

        end if

        shifted_left = adjustl ( current )
        if ( shifted_left(1:len(key)) == key ) exit

      end do

    end subroutine skip_lines


!   set index array for reordering nodes
!   gmsh_node = gmsh_inv_index ( tfem_node )

    subroutine fill_gmsh_inv_index ( elshape, p, gmsh_inv_index )

      use meshgen_gmsh_index_m

      integer, intent(in) :: elshape, p
      integer, dimension(:), intent(out) :: gmsh_inv_index

      integer :: i

      select case ( elshape )
        case (2) ! 3 node line element
          gmsh_inv_index = [ 1, 3, 2 ]
        case (4) ! 6 node triangle
          gmsh_inv_index = [ 1, 4, 2, 5, 3, 6 ]
        case (6) ! 9 node quadrilateral
          gmsh_inv_index = [ 1, 5, 2, 6, 3, 7, 4, 8, 9 ]
        case (12) ! 10 node tetrahedron
          gmsh_inv_index = [ 1, 5, 2, 6, 3, 7, 8, 10, 9, 4 ]
        case (14) ! 27 node hexahedron
          gmsh_inv_index = [ 1, 9, 2, 10, 21, 12, 4, 14, 3, &
                             11, 22, 13, 23, 27, 24, 16, 25, 15, &
                             5, 17, 6, 18, 26, 19, 8, 20, 7 ]
        case (30) ! 8 node quadrilateral
          gmsh_inv_index = [ 1, 5, 2, 6, 3, 7, 4, 8 ]
        case (31) ! 20 node hexahedron
          gmsh_inv_index = [ 1, 9, 2, 10, 12, 4, 14, 3, &
                             11, 13, 16, 15, &
                             5, 17, 6, 18, 19, 8, 20, 7 ]
        case (42) ! 15 node prism
          gmsh_inv_index = [ 1, 7, 2, 10, 3, 8, 9, 11, 12, 4, &
                             13, 5, 15, 6, 14 ]
        case (43) ! 18 node prism
          gmsh_inv_index = [ 1, 7, 2, 10, 3, 8, 9, 16, 11, 18, &
                             12, 17, 4, 13, 5, 15, 6, 14 ]
        case (52) ! 13 node pyramid
          gmsh_inv_index = [ 1, 6, 2, 9, 3, 11, 4, 7, 8, 10, &
                             12, 13, 5 ]
        case (53) ! 14 node pyramid
          gmsh_inv_index = [ 1, 6, 2, 9, 3, 11, 4, 7, 14, 8, &
                             10, 12, 13, 5 ]
        case (101) ! high order line element
          call gmsh_index_line ( p, gmsh_inv_index=gmsh_inv_index )
        case (102) ! high order quadrilateral element
          call gmsh_index_quadrilateral ( p, gmsh_inv_index=gmsh_inv_index )
        case (104) ! high order triangle element
          call gmsh_index_triangle ( p, gmsh_inv_index=gmsh_inv_index )
        case default
          gmsh_inv_index = [ ( i, i=1,size(gmsh_inv_index) ) ]
      end select

    end subroutine fill_gmsh_inv_index


!   fill a single geometry (curve or surface) with local nodes

    subroutine fill_local_geometry ( geometry )

      type(geometry_t), intent(inout) :: geometry

       integer :: elem1, elem2, elem, gnode

      do elem = 1, geometry%nelem
!       mark gnode with a unique increasing number for sorting
        elem1 = ( elem - 1 ) * geometry%element%numnod + 1
        elem2 = elem1 + geometry%element%numnod - 1
        gnodes ( geometry%topology(:,elem,2) ) = [ ( i, i=elem1,elem2 ) ]
      end do

      geometry%nnodes = count ( gnodes /= 0 )

      geometry%nnodes_blend = [0,geometry%nnodes]

      allocate ( geometry%nodes(geometry%nnodes) )

!     fill local nodes

      call sort ( gnodes, ordern )

      lnode = 0

      do gnode = 1, mesh%nnodes
        if ( gnodes(gnode) == 0 ) cycle  ! global node not in geometry
          lnode = lnode + 1 ! new local node
          geometry%nodes(lnode) = ordern(gnode)
          gnodes(gnode) = lnode ! store lnode
      end do

      if ( lnode /= geometry%nnodes ) then
        print *, 'lnode, nnodes', lnode, geometry%nnodes
        stop 'read_mesh_gmsh: internal error'
      end if

!     reorder gnodes

      gnodes(ordern) = gnodes

!     fill local topology

      do elem = 1, geometry%nelem
        geometry%topology(:,elem,1) = &
        gnodes ( geometry%topology(:,elem,2) )
      end do

!     set gnodes back to zero

      gnodes ( geometry%nodes ) = 0

    end subroutine fill_local_geometry

  end subroutine read_mesh_gmsh


! write a mesh to a gmsh file

  subroutine write_mesh_gmsh ( mesh, filename, groups, binary, blend )

    type(mesh_t), intent(in) :: mesh

!   the filename for writing the mesh to (it must have the .msh extension).
    character (len=*), intent(in) :: filename

!   if present: the element groups to be written
!   For example groups=(/2,4/) will write mesh for element groups 2 and 4.
!   Default: all groups are written
!   NOTE: all the nodal points of the mesh are written and only the elements
!   are affected by choosing a subset (of all groups) to be written.
    integer, dimension(:), intent(in), optional :: groups

!   if .true. write the mesh as a binary file.
!   default=.false.
!   The advantages of a binary file are:
!   a. the data is identical to the internal representation (no errors due
!      to truncation).
!   b. the occupied space is usually smaller.
!   c. reading the file can be somewhat faster.
!   A disadvantage is that the file cannot be interpreted by humans.
    logical, intent(in), optional :: binary

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
!   default=0
    integer, intent(in), optional :: blend

    integer :: i, ios, lblend
    integer, allocatable, dimension(:) :: lgroups
    logical :: lbinary

    call check ( mesh, 'write_mesh_gmsh' )

    if ( present(blend) ) then
      if ( blend < 0 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error write_mesh_gmsh: invalid argument ', &
          '  blend < 0 or larger than ', mesh%nblend
        stop
      end if
      lblend = blend
    else
      lblend = 0
    end if

!   check groups
    if ( present(groups) ) then
      if ( any ( groups < 1 ) .or. any ( groups > mesh%nelgrp ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: parameter groups in the heading of write_mesh_gmsh is ', &
          'out of range: some groups are < 1 or larger than ', mesh%nelgrp
        stop
      end if
      allocate(lgroups(size(groups)))
      lgroups = groups
    else
      allocate(lgroups(mesh%nelgrp))
      lgroups = [ (i,i=1,mesh%nelgrp) ]
    end if

    lbinary = set_optional ( variable=binary, default=.false. )

    if ( lbinary ) then

!     write binary file

      open ( unit=unit_gmsh, file=filename, access='stream', iostat=ios, &
             form='unformatted', status='replace' )

      if ( ios /= 0 ) then
        write(*,'(/a,i0/)') 'Error in write_mesh_gmsh: ios = ', ios
        stop
      end if

      call gmsh_header_binary
      call gmsh_nodes_elements_binary ( mesh, lgroups, lblend )

    else

!     write ASCII file

      open ( unit=unit_gmsh, file=filename, iostat=ios, status='replace' )

      if ( ios /= 0 ) then
        write(*,'(/a,i0/)') 'Error in write_mesh_gmsh: ios = ', ios
        stop
      end if

      call gmsh_header
      call gmsh_nodes_elements ( mesh, lgroups, lblend )

    end if

    deallocate(lgroups)

    close ( unit=unit_gmsh )

  end subroutine write_mesh_gmsh


! write geometries (points, curves, surfaces) to a gmsh file

  subroutine write_geometries_gmsh ( mesh, filename, binary, blend )

    type(mesh_t), intent(in) :: mesh

!   the filename for writing the geometries to (it must have a .msh extension).
    character (len=*), intent(in) :: filename

!   if .true. write the mesh as a binary file.
!   default=.false.
!   The advantages of a binary file are:
!   a. the data is identical to the internal representation (no errors due
!      to truncation).
!   b. the occupied space is usually smaller.
!   c. reading the file can be somewhat faster.
!   A disadvantage is that the file cannot be interpreted by humans.
    logical, intent(in), optional :: binary

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
!   default=0
    integer, intent(in), optional :: blend

    integer :: ios, lblend
    logical :: lbinary

    call check ( mesh, 'write_geometries_gmsh' )

    if ( present(blend) ) then
      if ( blend < 0 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error write_geometries_gmsh: invalid argument ', &
          '  blend < 0 or larger than ', mesh%nblend
        stop
      end if
      lblend = blend
    else
      lblend = 0
    end if

    lbinary = set_optional ( variable=binary, default=.false. )

    if ( lbinary ) then

!     write binary file

      open ( unit=unit_gmsh, file=filename, access='stream', iostat=ios, &
             form='unformatted', status='replace' )

      if ( ios /= 0 ) then
        write(*,'(/a,i0/)') 'Error in write_geometries_gmsh: ios = ', ios
        stop
      end if

      call gmsh_header_binary
      call gmsh_geometries_binary ( mesh, lblend )

    else

!     write ASCII file

      open ( unit=unit_gmsh, file=filename, iostat=ios, status='replace' )

      if ( ios /= 0 ) then
        write(*,'(/a,i0/)') 'Error in write_geometries_gmsh: ios = ', ios
        stop
      end if

      call gmsh_header
      call gmsh_geometries ( mesh, lblend )

    end if

    close ( unit=unit_gmsh )

  end subroutine write_geometries_gmsh


! gmsh header (ASCII format)

  subroutine gmsh_header

    character (len=21) :: line
    integer :: values(8)

!   generate date and time

    call date_and_time ( values=values )

    write ( line, '(2(i2.2,a),i4,a,3(i2.2,a))' ) &
      values(3), '-', values(2), '-', values(1), ', ', &
      values(5), ':', values(6), ':', values(7)

!   general header

    write ( unit=unit_gmsh, fmt='(a)' ) '$Comments'
    write ( unit=unit_gmsh, fmt='(2a)' ) 'TFEM output generated at ', line
    write ( unit=unit_gmsh, fmt='(a)' ) '$EndComments'

    write ( unit=unit_gmsh, fmt='(a)' ) '$MeshFormat'
    write ( unit=unit_gmsh, fmt='(2a)' ) '2.2 0 8'
    write ( unit=unit_gmsh, fmt='(a)' ) '$EndMeshFormat'

  end subroutine gmsh_header


! gmsh header (binary format)

  subroutine gmsh_header_binary

    character (len=21) :: line
    integer :: values(8)

!   generate date and time

    call date_and_time ( values=values )

    write ( line, '(2(i2.2,a),i4,a,3(i2.2,a))' ) &
      values(3), '-', values(2), '-', values(1), ', ', &
      values(5), ':', values(6), ':', values(7)

!   general header

    call write_line ( unit=unit_gmsh, line='$Comments' )
    call write_line ( unit=unit_gmsh, line='TFEM output generated at '//line )
    call write_line ( unit=unit_gmsh, line='$EndComments' )

    call write_line ( unit=unit_gmsh, line='$MeshFormat' )
    call write_line ( unit=unit_gmsh, line='2.2 1 8' )
    write ( unit=unit_gmsh ) 1; call write_line ( unit=unit_gmsh, line='' )
    call write_line ( unit=unit_gmsh, line='$EndMeshFormat' )

  end subroutine gmsh_header_binary


! gmsh nodes and elements (ASCII format mesh)

  subroutine gmsh_nodes_elements ( mesh, groups, blend )

    type(mesh_t), intent(in) :: mesh

!   element groups to be to be included
!   All nodes, points, curves and surfaces are written.
    integer, dimension(:), intent(in) :: groups

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
    integer, intent(in) :: blend

    integer :: elem, elgrp, i, j, numelem, grp, elemnr, point, curve, surface
    type(element_t), dimension(mesh%nelgrp) :: element
    type(element_t), dimension(mesh%ncurves) :: elementc
    type(element_t), dimension(mesh%nsurfaces) :: elements

!   choose element

    if ( blend == 0 ) then

      element = mesh%element
      do curve = 1, mesh%ncurves
        elementc(curve) = mesh%curves(curve)%element
      end do
      do surface = 1, mesh%nsurfaces
        elements(surface) = mesh%surfaces(surface)%element
      end do

    else

      element = mesh%element_blend(:,blend)
      do curve = 1, mesh%ncurves
        elementc(curve) = mesh%curves(curve)%element_blend(blend)
      end do
      do surface = 1, mesh%nsurfaces
        elements(surface) = mesh%surfaces(surface)%element_blend(blend)
      end do

    end if

!   test

    if ( any(element(groups)%gmsh_elshape == 0) ) then
      write (*,'(3(a/),2a)') 'Error gmsh_nodes_elements: ', &
        ' there are element groups that have element shapes that are', &
        ' not (yet) implemented for gmsh output.', &
        ' gmsh_elshape(:) ='
      write (*,*) element(groups)%gmsh_elshape
      stop
    end if

!   nodes

    write ( unit=unit_gmsh, fmt='(a)' ) '$Nodes'
    write ( unit=unit_gmsh, fmt='(i0)' ) &
                     mesh%nnodes_blend(blend+2)-mesh%nnodes_blend(blend+1)
    write ( unit=unit_gmsh, fmt='(i0,3es16.8)' ) &
      ( i-mesh%nnodes_blend(blend+1), mesh%coor(i,:), &
          ( 0._dp, j=mesh%ndim+1,3 ), &
          i=mesh%nnodes_blend(blend+1)+1,mesh%nnodes_blend(blend+2) )

    write ( unit=unit_gmsh, fmt='(a)' ) '$EndNodes'

    numelem = count ( mesh%points(1:mesh%npoints) > &
                           mesh%nnodes_blend(blend+1) .and. &
                      mesh%points(1:mesh%npoints) <= &
                           mesh%nnodes_blend(blend+2) ) + &
              sum ( mesh%curves(1:mesh%ncurves)%nelem ) + &
              sum ( mesh%surfaces(1:mesh%nsurfaces)%nelem ) + &
              sum ( mesh%grpnumel(groups) )

    write ( unit=unit_gmsh, fmt='(a)' ) '$Elements'
    write ( unit=unit_gmsh, fmt='(i0)' ) numelem

    elemnr = 0

!   points

    do point = 1, mesh%npoints
      if ( mesh%points(point) <= mesh%nnodes_blend(blend+1) .or. &
           mesh%points(point)  > mesh%nnodes_blend(blend+2) ) cycle
      elemnr = elemnr + 1
      write ( unit=unit_gmsh, fmt='(6(1x,i0))' ) &
         elemnr, 15, 2, 0, elemnr, mesh%points(point) &
                          -mesh%nnodes_blend(blend+1)
    end do

!   curves

    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        elemnr = elemnr + 1
        write ( unit=unit_gmsh, fmt='(11(1x,i0))' ) &
          elemnr, elementc(curve)%gmsh_elshape, &
            2, 0, curve, mesh%curves(curve)%topology(&
                    elementc(curve)%gmsh_index &
                    +mesh%curves(curve)%numnodtop(blend+1),elem,2) &
                          -mesh%nnodes_blend(blend+1)
      end do
    end do

!   surfaces

    do surface = 1, mesh%nsurfaces
      do elem = 1, mesh%surfaces(surface)%nelem
        elemnr = elemnr + 1
        write ( unit=unit_gmsh, fmt='(11(1x,i0))' ) &
          elemnr, elements(surface)%gmsh_elshape, &
            2, 0, surface, mesh%surfaces(surface)%topology(&
                         elements(surface)%gmsh_index &
                    +mesh%surfaces(surface)%numnodtop(blend+1),elem,2) &
                          -mesh%nnodes_blend(blend+1)
      end do
    end do

!   element groups

    do grp = 1, size(groups)
      elgrp = groups(grp)
      do elem = 1, mesh%grpnumel(elgrp)
        elemnr = elemnr + 1
        write ( unit=unit_gmsh, fmt='(11(1x,i0))' ) &
          elemnr, element(elgrp)%gmsh_elshape, &
            2, 0, elgrp, mesh%topology(elgrp)%a(&
              &element(elgrp)%gmsh_index+mesh%numnodtop(elgrp,blend+1),elem) &
               -mesh%nnodes_blend(blend+1)
      end do
    end do

    write ( unit=unit_gmsh, fmt='(a)' ) '$EndElements'

  end subroutine gmsh_nodes_elements


! gmsh nodes and elements (binary format mesh)

  subroutine gmsh_nodes_elements_binary ( mesh, groups, blend )

    type(mesh_t), intent(in) :: mesh

!   element groups to be to be included
!   All nodes, points, curves and surfaces are written.
    integer, dimension(:), intent(in) :: groups

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
    integer, intent(in) :: blend

    character (len=80) :: line
    integer(int32) :: i, elemnr, point, curve, surface, nume, elm_type, ntags, &
      elgrp
    integer :: elem, j, numelem, grp
    type(element_t), dimension(mesh%nelgrp) :: element
    type(element_t), dimension(mesh%ncurves) :: elementc
    type(element_t), dimension(mesh%nsurfaces) :: elements

!   choose element

    if ( blend == 0 ) then

      element = mesh%element
      do curve = 1, mesh%ncurves
        elementc(curve) = mesh%curves(curve)%element
      end do
      do surface = 1, mesh%nsurfaces
        elements(surface) = mesh%surfaces(surface)%element
      end do

    else

      element = mesh%element_blend(:,blend)
      do curve = 1, mesh%ncurves
        elementc(curve) = mesh%curves(curve)%element_blend(blend)
      end do
      do surface = 1, mesh%nsurfaces
        elements(surface) = mesh%surfaces(surface)%element_blend(blend)
      end do

    end if

!   test

    if ( any(element(groups)%gmsh_elshape == 0) ) then
      write (*,'(3(a/),2a)') 'Error gmsh_nodes_elements: ', &
        ' there are element groups that have element shapes that are', &
        ' not (yet) implemented for gmsh output.', &
        ' gmsh_elshape(:) ='
      write (*,*) element(groups)%gmsh_elshape
      stop
    end if

!   nodes

    call write_line ( unit=unit_gmsh, line='$Nodes' )
    write ( unit=line, fmt='(i0)' ) &
                     mesh%nnodes_blend(blend+2)-mesh%nnodes_blend(blend+1)
    call write_line ( unit=unit_gmsh, line=trim(line) )
    write ( unit=unit_gmsh ) &
      ( i-mesh%nnodes_blend(blend+1), mesh%coor(i,:), &
      ( 0._dp, j=mesh%ndim+1,3 ), &
          i=mesh%nnodes_blend(blend+1)+1,mesh%nnodes_blend(blend+2) )
    call write_line ( unit=unit_gmsh, line='' )
    call write_line ( unit=unit_gmsh, line='$EndNodes' )

    numelem = count ( mesh%points(1:mesh%npoints) > &
                           mesh%nnodes_blend(blend+1) .and. &
                      mesh%points(1:mesh%npoints) <= &
                           mesh%nnodes_blend(blend+2) ) + &
              sum ( mesh%curves(1:mesh%ncurves)%nelem ) + &
              sum ( mesh%surfaces(1:mesh%nsurfaces)%nelem ) + &
              sum ( mesh%grpnumel(groups) )

    call write_line ( unit=unit_gmsh, line='$Elements' )
    write ( unit=line, fmt='(i0)' ) numelem
    call write_line ( unit=unit_gmsh, line=trim(line) )

    ntags = 2

    elemnr = 0

!   points

    elm_type = 15
    nume = count ( mesh%points(1:mesh%npoints) > &
                           mesh%nnodes_blend(blend+1) .and. &
                      mesh%points(1:mesh%npoints) <= &
                           mesh%nnodes_blend(blend+2) )

    write ( unit=unit_gmsh ) elm_type, nume, ntags

    do point = 1, mesh%npoints
      if ( mesh%points(point) <= mesh%nnodes_blend(blend+1) .or. &
           mesh%points(point)  > mesh%nnodes_blend(blend+2) ) cycle
      elemnr = elemnr + 1
      write ( unit=unit_gmsh ) elemnr, 0_int32, elemnr, &
                      int(mesh%points(point)-mesh%nnodes_blend(blend+1),int32)
    end do

!   curves

    do curve = 1, mesh%ncurves

      elm_type = elementc(curve)%gmsh_elshape
      nume = mesh%curves(curve)%nelem

      write ( unit=unit_gmsh ) elm_type, nume, ntags

      do elem = 1, mesh%curves(curve)%nelem
        elemnr = elemnr + 1
        write ( unit=unit_gmsh ) elemnr, 0_int32, curve, &
            int(mesh%curves(curve)%topology(&
                elementc(curve)%gmsh_index&
                  +mesh%curves(curve)%numnodtop(blend+1),elem,2) &
                  -mesh%nnodes_blend(blend+1),int32)
      end do

    end do

!   surfaces

    do surface = 1, mesh%nsurfaces

      elm_type = elements(surface)%gmsh_elshape
      nume = mesh%surfaces(surface)%nelem

      write ( unit=unit_gmsh ) elm_type, nume, ntags

      do elem = 1, mesh%surfaces(surface)%nelem
        elemnr = elemnr + 1
        write ( unit=unit_gmsh ) elemnr, 0_int32, surface, &
          int(mesh%surfaces(surface)%topology(&
                      &elements(surface)%gmsh_index&
                        +mesh%surfaces(surface)%numnodtop(blend+1),elem,2) &
                  -mesh%nnodes_blend(blend+1),int32)
      end do

    end do

!   element groups

    do grp = 1, size(groups)

      elgrp = groups(grp)

      elm_type = element(elgrp)%gmsh_elshape
      nume = mesh%grpnumel(elgrp)

      write ( unit=unit_gmsh ) elm_type, nume, ntags

      do elem = 1, mesh%grpnumel(elgrp)
        elemnr = elemnr + 1
        write ( unit=unit_gmsh ) elemnr, 0_int32, elgrp, &
          int(mesh%topology(elgrp)%a(element(elgrp)%gmsh_index+&
                             mesh%numnodtop(elgrp,blend+1),elem)&
                                   -mesh%nnodes_blend(blend+1),int32)
      end do

    end do

    call write_line ( unit=unit_gmsh, line='' )
    call write_line ( unit=unit_gmsh, line='$EndElements' )

  end subroutine gmsh_nodes_elements_binary


! Write geometries only (ASCII format mesh)
! Only nodes on the geometries are exported to the file.

  subroutine gmsh_geometries ( mesh, blend )

    type(mesh_t), intent(in) :: mesh

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
    integer, intent(in) :: blend

    integer :: elem, i, j, numelem, elemnr, point, curve, surface
    integer :: nnodes, node, nump
    integer, allocatable, dimension(:) :: newnum
    type(element_t), dimension(mesh%ncurves) :: elementc
    type(element_t), dimension(mesh%nsurfaces) :: elements

!   choose element

    if ( blend == 0 ) then

      do curve = 1, mesh%ncurves
        elementc(curve) = mesh%curves(curve)%element
      end do
      do surface = 1, mesh%nsurfaces
        elements(surface) = mesh%surfaces(surface)%element
      end do

    else

      do curve = 1, mesh%ncurves
        elementc(curve) = mesh%curves(curve)%element_blend(blend)
      end do
      do surface = 1, mesh%nsurfaces
        elements(surface) = mesh%surfaces(surface)%element_blend(blend)
      end do

    end if

    allocate ( newnum(mesh%nnodes) )

!   set newnum (new numbering of the nodes)

    newnum = 0

!   add nodes in points

    nump = 0
    do point = 1, mesh%npoints
      if ( mesh%points(point) <= mesh%nnodes_blend(blend+1) .or. &
           mesh%points(point)  > mesh%nnodes_blend(blend+2) ) cycle
      newnum(mesh%points(point)) = 1
      nump = nump + 1
    end do

!   add nodes on curves

    do curve = 1, mesh%ncurves
      newnum(mesh%curves(curve)%nodes) = 1
    end do

!   add nodes on surfaces

    do surface = 1, mesh%nsurfaces
      newnum(mesh%surfaces(surface)%nodes) = 1
    end do

!   only include the blend mesh required

    newnum(1:mesh%nnodes_blend(blend+1)) = 0
    newnum(mesh%nnodes_blend(blend+2)+1:) = 0

!   set new numbering

    nnodes = 0
    do i = 1, mesh%nnodes
      if ( newnum(i) == 0 ) cycle
      nnodes = nnodes + 1  ! node found
      newnum(i) = nnodes
    end do


!   write nodes

    write ( unit=unit_gmsh, fmt='(a)' ) '$Nodes'
    write ( unit=unit_gmsh, fmt='(i0)' ) nnodes
    node = 0
    do i = 1, mesh%nnodes
      if ( newnum(i) == 0 ) cycle
      node = node + 1
      write ( unit=unit_gmsh, fmt='(i0,3es16.8)' ) &
                       node, mesh%coor(i,:), ( 0._dp, j=mesh%ndim+1,3 )
    end do
    write ( unit=unit_gmsh, fmt='(a)' ) '$EndNodes'

    numelem = nump + sum ( mesh%curves(1:mesh%ncurves)%nelem ) + &
                     sum ( mesh%surfaces(1:mesh%nsurfaces)%nelem )

    write ( unit=unit_gmsh, fmt='(a)' ) '$Elements'
    write ( unit=unit_gmsh, fmt='(i0)' ) numelem

    elemnr = 0

!   points

    do point = 1, mesh%npoints
      if ( mesh%points(point) <= mesh%nnodes_blend(blend+1) .or. &
           mesh%points(point)  > mesh%nnodes_blend(blend+2) ) cycle
      elemnr = elemnr + 1
      write ( unit=unit_gmsh, fmt='(6(1x,i0))' ) &
         elemnr, 15, 2, 0, elemnr, newnum(mesh%points(point))
    end do

!   curves

    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        elemnr = elemnr + 1
        write ( unit=unit_gmsh, fmt='(11(1x,i0))' ) &
          elemnr, elementc(curve)%gmsh_elshape, &
            2, 0, curve, newnum(mesh%curves(curve)%topology(&
                              elementc(curve)%gmsh_index&
                        +mesh%curves(curve)%numnodtop(blend+1),elem,2))
      end do
    end do

!   surfaces

    do surface = 1, mesh%nsurfaces
      do elem = 1, mesh%surfaces(surface)%nelem
        elemnr = elemnr + 1
        write ( unit=unit_gmsh, fmt='(11(1x,i0))' ) &
          elemnr, elements(surface)%gmsh_elshape, &
            2, 0, surface, newnum(mesh%surfaces(surface)%topology(&
                              elements(surface)%gmsh_index&
                 +mesh%surfaces(surface)%numnodtop(blend+1),elem,2))
      end do
    end do

    write ( unit=unit_gmsh, fmt='(a)' ) '$EndElements'

    deallocate ( newnum )

  end subroutine gmsh_geometries


! Write geometries only (binary format mesh)
! Only nodes on the geometries are exported to the file.

  subroutine gmsh_geometries_binary ( mesh, blend )

    type(mesh_t), intent(in) :: mesh

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
    integer, intent(in) :: blend

    character (len=80) :: line
    integer(int32) :: node, nume, elm_type, ntags, point, curve, surface, elemnr
    integer :: elem, i, j, numelem
    integer :: nnodes, nump
    integer(int32), allocatable, dimension(:) :: newnum

    type(element_t), dimension(mesh%ncurves) :: elementc
    type(element_t), dimension(mesh%nsurfaces) :: elements

!   choose element

    if ( blend == 0 ) then

      do curve = 1, mesh%ncurves
        elementc(curve) = mesh%curves(curve)%element
      end do
      do surface = 1, mesh%nsurfaces
        elements(surface) = mesh%surfaces(surface)%element
      end do

    else

      do curve = 1, mesh%ncurves
        elementc(curve) = mesh%curves(curve)%element_blend(blend)
      end do
      do surface = 1, mesh%nsurfaces
        elements(surface) = mesh%surfaces(surface)%element_blend(blend)
      end do

    end if

    allocate ( newnum(mesh%nnodes) )

!   set newnum (new numbering of the nodes)

    newnum = 0

!   add nodes in points

    nump = 0
    do point = 1, mesh%npoints
      if ( mesh%points(point) <= mesh%nnodes_blend(blend+1) .or. &
           mesh%points(point)  > mesh%nnodes_blend(blend+2) ) cycle
      newnum(mesh%points(point)) = 1
      nump = nump + 1
    end do

!   add nodes on curves

    do curve = 1, mesh%ncurves
      newnum(mesh%curves(curve)%nodes) = 1
    end do

!   add nodes on surfaces

    do surface = 1, mesh%nsurfaces
      newnum(mesh%surfaces(surface)%nodes) = 1
    end do

!   only include the blend mesh required

    newnum(1:mesh%nnodes_blend(blend+1)) = 0
    newnum(mesh%nnodes_blend(blend+2)+1:) = 0

!   set new numbering

    nnodes = 0
    do i = 1, mesh%nnodes
      if ( newnum(i) == 0 ) cycle
      nnodes = nnodes + 1  ! node found
      newnum(i) = nnodes
    end do


!   write nodes

    call write_line ( unit=unit_gmsh, line='$Nodes' )
    write ( unit=line, fmt='(i0)' ) nnodes
    call write_line ( unit=unit_gmsh, line=trim(line) )
    node = 0
    do i = 1, mesh%nnodes
      if ( newnum(i) == 0 ) cycle
      node = node + 1
      write ( unit=unit_gmsh ) node, mesh%coor(i,:), ( 0._dp, j=mesh%ndim+1,3 )
    end do
    call write_line ( unit=unit_gmsh, line='' )
    call write_line ( unit=unit_gmsh, line='$EndNodes' )

    numelem = nump + sum ( mesh%curves(1:mesh%ncurves)%nelem ) + &
                     sum ( mesh%surfaces(1:mesh%nsurfaces)%nelem )

    call write_line ( unit=unit_gmsh, line='$Elements' )
    write ( unit=line, fmt='(i0)' ) numelem
    call write_line ( unit=unit_gmsh, line=trim(line) )

    ntags = 2

    elemnr = 0

!   points

    elm_type = 15
    nume = nump

    write ( unit=unit_gmsh ) elm_type, nume, ntags

    do point = 1, mesh%npoints
      if ( mesh%points(point) <= mesh%nnodes_blend(blend+1) .or. &
           mesh%points(point)  > mesh%nnodes_blend(blend+2) ) cycle
      elemnr = elemnr + 1
      write ( unit=unit_gmsh ) elemnr, 0_int32, elemnr, &
                               newnum(mesh%points(point))
    end do

!   curves

    do curve = 1, mesh%ncurves

      elm_type = elementc(curve)%gmsh_elshape
      nume = mesh%curves(curve)%nelem

      write ( unit=unit_gmsh ) elm_type, nume, ntags

      do elem = 1, mesh%curves(curve)%nelem
        elemnr = elemnr + 1
        write ( unit=unit_gmsh ) elemnr, 0_int32, curve, &
            newnum(mesh%curves(curve)%topology(elementc(curve)%gmsh_index &
                     +mesh%curves(curve)%numnodtop(blend+1),elem,2))
      end do

    end do

!   surfaces

    do surface = 1, mesh%nsurfaces

      elm_type = elements(surface)%gmsh_elshape
      nume = mesh%surfaces(surface)%nelem

      write ( unit=unit_gmsh ) elm_type, nume, ntags

      do elem = 1, mesh%surfaces(surface)%nelem
        elemnr = elemnr + 1
        write ( unit=unit_gmsh ) elemnr, 0_int32, surface, &
          newnum(mesh%surfaces(surface)%topology(&
                            &elements(surface)%gmsh_index &
                  +mesh%surfaces(surface)%numnodtop(blend+1),elem,2))
      end do

    end do

    call write_line ( unit=unit_gmsh, line='' )
    call write_line ( unit=unit_gmsh, line='$EndElements' )

    deallocate ( newnum )

  end subroutine gmsh_geometries_binary


! write a scalar field to a gmsh file

  subroutine write_scalar_gmsh ( mesh, problem, filename, dataname, physq, &
    degfd, sysvector, vector, append, dataonly, groups, timestepnr, time, &
    binary, blend, offset )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the data to (it must have the .msh extension).
    character (len=*), intent(in) :: filename

!   the dataname for the scalar, for example 'pressure'.
    character (len=*), intent(in) :: dataname

!   the physical quantity to be written
    integer, intent(in), optional :: physq

!   the degree of freedom to be written (within the physical quantity if
!   physq is also present)
    integer, intent(in), optional :: degfd

!   sysvector to be written
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be written
    type(vector_t), intent(in), optional :: vector

!   if .true.: append data to existing file in order to avoid writing
!   coordinate and topology data more than ones.
!   default=.false.
    logical, intent(in), optional :: append

!   if .true.: write only data to existing file (no mesh). In this way the mesh
!   and the data can be written to separate files.
!   default=.false.
    logical, intent(in), optional :: dataonly

!   if present: the element groups to be written
!   For example groups=(/2,4/) will write data for element groups 2 and 4.
!   Default: all groups are written
!   NOTE: the data in all the nodal points of the mesh is written and only
!   the elements are affected by choosing a subset (of all groups) to be
!   written.
    integer, dimension(:), intent(in), optional :: groups

!   time step number/index (starting with 0)
!   default=0
    integer, intent(in), optional :: timestepnr

!   the (physical) time value associated with the data
!   default=0.0
    real(dp), intent(in), optional :: time

!   if .true. write the data as a binary file.
!   default=.false.
!   The advantages of a binary file are:
!   a. the data is identical to the internal representation (no errors due
!      to truncation).
!   b. the occupied space is usually smaller.
!   c. reading the file can be somewhat faster.
!   A disadvantage is that the file cannot be interpreted by humans.
    logical, intent(in), optional :: binary

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
!   default=0
!   Notes:
!     - blend has no effect if vector%elementwise=.true.
!     - blend should not change its value for "append" calls.
    integer, intent(in), optional :: blend

!   In case
!    1. the structure problem is created using a mesh consisting of
!       blended meshes,
!    2. the data corresponds to a blended mesh (so not to the main mesh),
!    3. a "fake" mesh is used that only contains the nodes corresponding
!       to the data
!   the value of offset should be equal to the number nodes missing from
!   the main mesh and the blended meshes "below" the current mesh. The value
!   of offest should be 0 if the "fake" mesh is based on the main mesh.
!   default = 0
!   NOTE: the value for offset is given by mesh%nnodes_blend(m+2), where
!   mesh is the blended mesh and m is the blend mesh just "below" the
!   current blend mesh (note that m=0 represents the main mesh).
!   Example: for plotting higher-order pressures in a velocity/pressure
!   Taylor-Hood formulation on its corresponding grid, offset should be
!   given the value of the number of nodes of the velocity grid.
    integer, intent(in), optional :: offset


    logical :: fileexists, fileappend, firstskip, filedataonly
    logical :: agroups(mesh%nelgrp), elementwise
    logical :: lbinary
    integer :: i, ag, elgrp, ndofu, elemnr, elem, s, numnod
    integer :: nodenr, pos(problem%maxnoddegfd), dof, deg
    integer :: ltimestepnr, lblend, loffset
    integer, allocatable, dimension(:) :: lgroups

    real(dp) :: svalue, u(problem%maxvecnoddegfd*mesh%maxelnumnod)
    real(dp) :: ltime


    call check ( mesh, 'write_scalar_gmsh' )
    call check ( problem, 'write_scalar_gmsh', mesh )

    if ( present(blend) ) then
      if ( blend < 0 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error write_scalar_gmsh: invalid argument ', &
          '  blend < 0 or larger than ', mesh%nblend
        stop
      end if
      lblend = blend
    else
      lblend = 0
    end if

    loffset = set_optional ( variable=offset, default=0 )

    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in write_scalar_gmsh: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(sysvector) ) then
      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in write_scalar_gmsh: sysvector has not been created '
        stop
      end if
      if ( problem%probnr /= sysvector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_scalar_gmsh: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in sysvector = ',  sysvector%probnr
        stop
      end if
      elementwise = .false.
    else if ( present(vector) ) then
      if ( .not. vector%created ) then
        write(*,'(/a/)') &
          'Error in write_scalar_gmsh: vector has not been created '
        stop
      end if
      if ( problem%probnr /= vector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_scalar_gmsh: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in vector = ', vector%probnr
        stop
      end if
      elementwise = vector%elementwise
    end if

!   check groups
    if ( present(groups) ) then
      if ( any ( groups < 1 ) .or. any ( groups > mesh%nelgrp ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: parameter groups in the heading of write_scalar_gmsh is ', &
          'out of range: some groups are < 1 or larger than ', mesh%nelgrp
        stop
      end if
      if ( problem%numinactivegroups == 0 ) then
        allocate(lgroups(size(groups)))
        lgroups = groups
      else
        agroups = .false.
        agroups(groups) = .true.
        agroups(problem%activegroups) = .true.
        allocate( lgroups(count(agroups)) )
        ag = 0
        do elgrp = 1, mesh%nelgrp
          if ( agroups(elgrp) ) then
            ag = ag + 1
            lgroups(ag) = elgrp
          end if
        end do
      end if
    else if ( problem%numinactivegroups > 0 ) then
      allocate( lgroups(size(problem%activegroups)) )
      lgroups = problem%activegroups
    else
      allocate(lgroups(mesh%nelgrp))
      lgroups = [ (i,i=1,mesh%nelgrp) ]
    end if

!   set deg

    if ( present(degfd) ) then
      if ( degfd < 1 ) then
        write(*,'(/a/)') &
          'Error in write_scalar_gmsh: degfd must be larger than zero'
        stop
      end if
      deg = degfd
    else
      deg = 1 ! take first degree
    end if

!   set fileappend, filedataonly

    fileappend = set_optional ( variable=append, default=.false. )
    filedataonly = set_optional ( variable=dataonly, default=.false. )

    lbinary = set_optional ( variable=binary, default=.false. )

!   open file

    inquire ( file=filename, exist=fileexists )

    if ( fileexists .and. fileappend ) then

!     add scalar point data to existing file

      if ( lbinary ) then

        open ( unit=unit_gmsh, file=filename, status='old', position='append', &
          access='stream', form='unformatted' )

      else

        open ( unit=unit_gmsh, file=filename, status='old', position='append' )

      end if

    else

!     add scalar point data to new file

      if ( lbinary ) then

        open ( unit=unit_gmsh, file=filename, status='replace', &
          access='stream', form='unformatted' )

        call gmsh_header_binary
        if ( .not. filedataonly ) &
                     call gmsh_nodes_elements_binary ( mesh, lgroups, lblend )

      else

        open ( unit=unit_gmsh, file=filename, status='replace' )

        call gmsh_header
        if ( .not. filedataonly ) &
                           call gmsh_nodes_elements ( mesh, lgroups, lblend )

      end if


    end if


!   write first part

    ltime = set_optional ( variable=time, default=0.0_dp )
    ltimestepnr = set_optional ( variable=timestepnr, default=0 )

    if ( lbinary ) then
      call write_dataheader_gmsh_binary ( elementwise, trim(dataname), &
        ltimestepnr, ltime, nfc=1, nelem=sum(mesh%grpnumel(lgroups)), &
        nnodes=mesh%nnodes_blend(lblend+2)-mesh%nnodes_blend(lblend+1) )
    else
      call write_dataheader_gmsh ( elementwise, trim(dataname), ltimestepnr, &
        ltime, nfc=1, nelem=sum(mesh%grpnumel(lgroups)), &
        nnodes=mesh%nnodes_blend(lblend+2)-mesh%nnodes_blend(lblend+1) )
    end if


!   write scalar

    if ( present(sysvector) ) then

!     sysvector in nodes

      firstskip = .true.

      do nodenr = mesh%nnodes_blend(lblend+1)+1, mesh%nnodes_blend(lblend+2)

        if ( present(physq) ) then
          call pos_array_node ( problem, loffset+nodenr, dof, pos, &
            physqarr=[physq] )
        else
          call pos_array_node ( problem, loffset+nodenr, dof, pos )
        end if

        if ( deg > dof ) then

          if ( firstskip .and. WARN_ON_NO_DATA_IN_NODES ) then
            write(*,'(/3(a/),2a/)') &
              'Warning in write_scalar_gmsh:', &
              ' there are nodes where the required unknown in the ', &
              ' sysvector is undefined. The scalar value is set to zero. ', &
              ' dataname = ', trim(dataname)
          end if
          firstskip = .false.
          svalue = 0

        else

          svalue = sysvector%u( pos(deg) )

        end if

        if ( lbinary ) then
          write ( unit=unit_gmsh ) int(nodenr,int32), svalue
        else
          write ( unit=unit_gmsh, fmt='(i0,es16.8)' ) &
                                           nodenr, real ( svalue, kind=sp )
        end if

      end do

    else if ( present(vector) .and. elementwise ) then

!     vector in nodes given per element

!     start at internal elements
      elemnr = count ( mesh%points(1:mesh%npoints) > 0 ) + &
               sum ( mesh%curves(1:mesh%ncurves)%nelem ) + &
               sum ( mesh%surfaces(1:mesh%nsurfaces)%nelem )

      do elgrp = 1, mesh%nelgrp

        if ( .not. any ( lgroups == elgrp ) ) cycle

        numnod = mesh%element(elgrp)%numnod

        s = ( deg - 1 ) * numnod ! start pointer

        do elem = 1, mesh%grpnumel(elgrp)

!         get element vector
          call get_vector ( mesh, problem, vector, elgrp, elem, u, &
            ndofu=ndofu, sloppy=.true. )

          if ( ndofu < deg * numnod ) then
            write(*,'(/2(a/),3a,i0/)') &
              'Error in write_scalar_gmsh:', &
              ' Not enough data in nodal points of element.', &
              ' Element group ', elgrp, ' dataname = ', trim(dataname)
            stop
          else if ( elem == 1 .and. mod(ndofu,numnod) /= 0 ) then
            write(*,'(/6(a/),a,i0/)') &
              'Error in write_scalar_gmsh:', &
              ' number of unknowns in the vector ', &
              ' is not a multiple of the number of nodes.', &
              ' write_scalar_gmsh assumes that the vector is defined ', &
              ' in all nodes having the same number of degrees of freedom.', &
              ' Element group ', elgrp
            stop
          end if

          elemnr = elemnr + 1

          if ( lbinary ) then
            write ( unit=unit_gmsh ) int([elemnr,numnod],int32)
            write ( unit=unit_gmsh ) u(s+mesh%element(elgrp)%gmsh_index)
          else
            write ( unit=unit_gmsh, fmt='(i0,1x,i0)' ) elemnr, numnod
            write ( unit=unit_gmsh, fmt='(10es16.8)' ) &
                real ( u(s+mesh%element(elgrp)%gmsh_index), kind=sp )
          end if

        end do

      end do

    else

!     vector in nodes

      firstskip = .true.

      do nodenr = mesh%nnodes_blend(lblend+1)+1, mesh%nnodes_blend(lblend+2)

        dof = problem%vec_nodnumdegfd(loffset+nodenr+1,vector%vec) &
                   - problem%vec_nodnumdegfd(loffset+nodenr,vector%vec)

        if ( deg > dof ) then

          if ( firstskip .and. WARN_ON_NO_DATA_IN_NODES ) then
            write(*,'(/3(a/),2a/)') &
              'Warning in write_scalar_gmsh:', &
              ' there are nodes where the required unknown in the ', &
              ' vector is undefined. The scalar value is set to zero. ', &
              ' dataname = ', trim(dataname)
          end if
          firstskip = .false.
          svalue = 0

        else

          svalue = vector%u( &
                    problem%vec_nodnumdegfd(loffset+nodenr,vector%vec) + deg )

        end if

        if ( lbinary ) then
          write ( unit=unit_gmsh ) int(nodenr,int32), svalue
        else
          write ( unit=unit_gmsh, fmt='(i0,es16.8)' ) &
                                       nodenr, real ( svalue, kind=sp )
        end if

      end do

    end if

    if ( elementwise ) then
      if ( lbinary ) then
        call write_line ( unit=unit_gmsh, line='' )
        call write_line ( unit=unit_gmsh, line='$EndElementNodeData' )
      else
        write ( unit=unit_gmsh, fmt='(a)' ) '$EndElementNodeData'
      end if
    else
      if ( lbinary ) then
        call write_line ( unit=unit_gmsh, line='' )
        call write_line ( unit=unit_gmsh, line='$EndNodeData' )
      else
        write ( unit=unit_gmsh, fmt='(a)' ) '$EndNodeData'
      end if
    end if

    deallocate ( lgroups )

    close ( unit=unit_gmsh )

  end subroutine write_scalar_gmsh


! write data header for ASCII file

  subroutine write_dataheader_gmsh ( elementwise, dataname, timestepnr, time, &
    nfc, nelem, nnodes )

    logical, intent(in) :: elementwise
    character (len=*), intent(in) :: dataname
    integer, intent(in) :: timestepnr, nfc
    real(dp), intent(in) :: time
    integer, intent(in) :: nelem, nnodes

    if ( elementwise ) then

      write ( unit=unit_gmsh, fmt='(a)' ) '$ElementNodeData'
      write ( unit=unit_gmsh, fmt='(i0)' ) 1   ! number of string tags
      write ( unit=unit_gmsh, fmt='(3a)' ) '"', dataname, '"'
      write ( unit=unit_gmsh, fmt='(i0)' ) 1   ! number of real tags
      write ( unit=unit_gmsh, fmt='(es12.4)' ) time
      write ( unit=unit_gmsh, fmt='(i0)' ) 3   ! number of integer tags
      write ( unit=unit_gmsh, fmt='(i0)' ) timestepnr   ! time step number
      write ( unit=unit_gmsh, fmt='(i0)' ) nfc   ! number of field components
      write ( unit=unit_gmsh, fmt='(i0)' ) nelem  ! number of elements

    else

      write ( unit=unit_gmsh, fmt='(a)' ) '$NodeData'
      write ( unit=unit_gmsh, fmt='(i0)' ) 1   ! number of string tags
      write ( unit=unit_gmsh, fmt='(3a)' ) '"', dataname, '"'
      write ( unit=unit_gmsh, fmt='(i0)' ) 1   ! number of real tags
      write ( unit=unit_gmsh, fmt='(es12.4)' ) time
      write ( unit=unit_gmsh, fmt='(i0)' ) 3   ! number of integer tags
      write ( unit=unit_gmsh, fmt='(i0)' ) timestepnr   ! time step number
      write ( unit=unit_gmsh, fmt='(i0)' ) nfc   ! number of field components
      write ( unit=unit_gmsh, fmt='(i0)' ) nnodes   ! number of nodes

    end if

  end subroutine write_dataheader_gmsh


! write data header for binary file

  subroutine write_dataheader_gmsh_binary ( elementwise, dataname, timestepnr,&
    time, nfc, nelem, nnodes )

    logical, intent(in) :: elementwise
    character (len=*), intent(in) :: dataname
    integer, intent(in) :: timestepnr, nfc
    real(dp), intent(in) :: time
    integer, intent(in) :: nelem, nnodes

    character (len=80) :: line


    if ( elementwise ) then

      call write_line ( unit=unit_gmsh, line='$ElementNodeData' )
      call write_line ( unit=unit_gmsh, line='1' ) ! number of string tags
      call write_line ( unit=unit_gmsh, line='"'//dataname//'"' )
      call write_line ( unit=unit_gmsh, line='1' ) ! number of real tags
      write ( unit=line, fmt='(es12.4)' ) time
      call write_line ( unit=unit_gmsh, line=trim(line) )
      call write_line ( unit=unit_gmsh, line='3' ) ! number of integer tags
      write ( unit=line, fmt='(i0)' ) timestepnr   ! time step number
      call write_line ( unit=unit_gmsh, line=trim(line) )
      write ( unit=line, fmt='(i0)' ) nfc   ! number of field components
      call write_line ( unit=unit_gmsh, line=trim(line) )
      write ( unit=line, fmt='(i0)' ) nelem  ! number of elements
      call write_line ( unit=unit_gmsh, line=trim(line) )

    else

      call write_line ( unit=unit_gmsh, line='$NodeData' )
      call write_line ( unit=unit_gmsh, line='1' ) ! number of string tags
      call write_line ( unit=unit_gmsh, line='"'//dataname//'"' )
      call write_line ( unit=unit_gmsh, line='1' ) ! number of real tags
      write ( unit=line, fmt='(es12.4)' ) time
      call write_line ( unit=unit_gmsh, line=trim(line) )
      call write_line ( unit=unit_gmsh, line='3' ) ! number of integer tags
      write ( unit=line, fmt='(i0)' ) timestepnr   ! time step number
      call write_line ( unit=unit_gmsh, line=trim(line) )
      write ( unit=line, fmt='(i0)' ) nfc   ! number of field components
      call write_line ( unit=unit_gmsh, line=trim(line) )
      write ( unit=line, fmt='(i0)' ) nnodes   ! number of nodes
      call write_line ( unit=unit_gmsh, line=trim(line) )

    end if

  end subroutine write_dataheader_gmsh_binary


! write a vector field to a gmsh file

  subroutine write_vector_gmsh ( mesh, problem, filename, dataname, physq, &
    degfd, sysvector, vector, append, dataonly, groups, timestepnr, time, &
    binary, assume2D, assume3D, blend, offset )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the data to (it must have the .msh extension).
    character (len=*), intent(in) :: filename

!   the dataname for the vector, for example 'velocity'.
    character (len=*), intent(in) :: dataname

!   the physical quantity to be written
    integer, intent(in), optional :: physq

!   the degrees of freedom to be written as a vector (within the physical
!   quantity if physq is also present).
!   Setting one of the components to zero means this component is
!   set to zero.
    integer, intent(in), dimension(:), optional :: degfd

!   sysvector to be written
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be written
    type(vector_t), intent(in), optional :: vector

!   if .true.: append data to existing file in order to avoid writing
!   coordinate and topology data more than ones.
!   default=.false.
    logical, intent(in), optional :: append

!   if .true.: write only data to existing file (no mesh). In this way the mesh
!   and the data can be written to separate files.
!   default=.false.
    logical, intent(in), optional :: dataonly

!   if present: the element groups to be written
!   For example groups=(/2,4/) will write data for element groups 2 and 4.
!   Default: all groups are written
!   NOTE: the data in all the nodal points of the mesh is written and only
!   the elements are affected by choosing a subset (of all groups) to be
!   written.
    integer, dimension(:), intent(in), optional :: groups

!   time step number/index (starting with 0)
!   default=0
    integer, intent(in), optional :: timestepnr

!   the (physical) time value associated with the data
!   default=0.0
    real(dp), intent(in), optional :: time

!   if .true. write the data as a binary file.
!   default=.false.
!   The advantages of a binary file are:
!   a. the data is identical to the internal representation (no errors due
!      to truncation).
!   b. the occupied space is usually smaller.
!   c. reading the file can be somewhat faster.
!   A disadvantage is that the file cannot be interpreted by humans.
    logical, intent(in), optional :: binary

!   Assume a 2D vector (even if mesh%ndim /= 2).
!   default=.false.
    logical, intent(in), optional :: assume2D

!   Assume a 3D vector (even if mesh%ndim /= 3).
!   default=.false.
    logical, intent(in), optional :: assume3D

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
!   default=0
!   Notes:
!     - blend has no effect if vector%elementwise=.true.
!     - blend should not change its value for "append" calls.
    integer, intent(in), optional :: blend

!   In case
!    1. the structure problem is created using a mesh consisting of
!       blended meshes,
!    2. the data corresponds to a blended mesh (so not to the main mesh),
!    3. a "fake" mesh is used that only contains the nodes corresponding
!       to the data
!   the value of offset should be equal to the number nodes missing from
!   the main mesh and the blended meshes "below" the current mesh. The value
!   of offest should be 0 if the "fake" mesh is based on the main mesh.
!   default = 0
!   NOTE: the value for offset is given by mesh%nnodes_blend(m+2), where
!   mesh is the blended mesh and m is the blend mesh just "below" the
!   current blend mesh (note that m=0 represents the main mesh).
!   Example: for plotting higher-order pressures in a velocity/pressure
!   Taylor-Hood formulation on its corresponding grid, offset should be
!   given the value of the number of nodes of the velocity grid.
    integer, intent(in), optional :: offset


    logical :: fileexists, fileappend, firstskip, filedataonly
    logical :: agroups(mesh%nelgrp), elementwise
    logical :: lbinary, lassume2D, lassume3D
    integer :: i, ag, elgrp, ndofu, elemnr, elem, s(3), numnod, j, lblend
    integer :: nodenr, pos(problem%maxnoddegfd), dof, deg(3), deg1(3)
    integer, allocatable, dimension(:) :: lgroups
    integer :: ltimestepnr, loffset

    real(dp) :: vvalue(3), u(problem%maxvecnoddegfd*mesh%maxelnumnod)
    real(dp), dimension(:,:), allocatable :: v
    real(dp) :: ltime

!   some checking

    call check ( mesh, 'write_vector_gmsh' )
    call check ( problem, 'write_vector_gmsh', mesh )

    lassume2D = set_optional ( variable=assume2D, default=.false. )
    lassume3D = set_optional ( variable=assume3D, default=.false. )

    if ( lassume2D .and. lassume3D ) then
      write(*,'(/a/a/)') &
        'Error in write_vector_gmsh: only one of assume2D and assume3D ', &
        ' can be set .true.'
      stop
    end if

    if ( present(blend) ) then
      if ( blend < 0 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error write_vector_gmsh: invalid argument ', &
          '  blend < 0 or larger than ', mesh%nblend
        stop
      end if
      lblend = blend
    else
      lblend = 0
    end if

    loffset = set_optional ( variable=offset, default=0 )

    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in write_vector_gmsh: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(sysvector) ) then
      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in write_vector_gmsh: sysvector has not been created '
        stop
      end if
      if ( problem%probnr /= sysvector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_vector_gmsh: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in sysvector = ',  sysvector%probnr
        stop
      end if
      elementwise = .false.
    else if ( present(vector) ) then
      if ( .not. vector%created ) then
        write(*,'(/a/)') &
          'Error in write_vector_gmsh: vector has not been created '
        stop
      end if
      if ( problem%probnr /= vector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_vector_gmsh:: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in vector = ', vector%probnr
        stop
      end if
      elementwise = vector%elementwise
    end if

!   check groups
    if ( present(groups) ) then
      if ( any ( groups < 1 ) .or. any ( groups > mesh%nelgrp ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: parameter groups in the heading of write_vector_gmsh is ', &
          'out of range: some groups are < 1 or larger than ', mesh%nelgrp
        stop
      end if
      if ( problem%numinactivegroups == 0 ) then
        allocate(lgroups(size(groups)))
        lgroups = groups
      else
        agroups = .false.
        agroups(groups) = .true.
        agroups(problem%activegroups) = .true.
        allocate( lgroups(count(agroups)) )
        ag = 0
        do elgrp = 1, mesh%nelgrp
          if ( agroups(elgrp) ) then
            ag = ag + 1
            lgroups(ag) = elgrp
          end if
        end do
      end if
    else if ( problem%numinactivegroups > 0 ) then
      allocate( lgroups(size(problem%activegroups)) )
      lgroups = problem%activegroups
    else
      allocate(lgroups(mesh%nelgrp))
      lgroups = [ (i,i=1,mesh%nelgrp) ]
    end if

!   Set deg

    if ( present(degfd) ) then

      if ( size(degfd) > 3 ) then
        write(*,'(/a/)') &
          'Error in write_vector_gmsh: dimension degfd larger than 3'
        stop
      end if

      if ( any(degfd < 0 ) ) then
        write(*,'(/a/)') &
          'Error in write_vector_gmsh: values in degfd must be not be negative'
        stop
      end if

      deg(1:size(degfd)) = degfd
      deg(size(degfd)+1:3) = 0

    else if ( mesh%ndim == 3 .or. lassume3D ) then

      deg = [1,2,3]

    else if ( mesh%ndim == 2 .or. lassume2D ) then

      deg = [1,2,0]

    else if ( mesh%ndim == 1 ) then

      deg = [1,0,0]

    end if

!   Set deg1 (to avoid zero index)

    where ( deg == 0 )
      deg1 = 1
    else where
      deg1 = deg
    end where

!   set fileappend, filedataonly

    fileappend = set_optional ( variable=append, default=.false. )
    filedataonly = set_optional ( variable=dataonly, default=.false. )

    lbinary = set_optional ( variable=binary, default=.false. )

!   open file

    inquire ( file=filename, exist=fileexists )

    if ( fileexists .and. fileappend ) then

!     add vector point data to existing file

      if ( lbinary ) then

        open ( unit=unit_gmsh, file=filename, status='old', position='append', &
          access='stream', form='unformatted' )

      else

        open ( unit=unit_gmsh, file=filename, status='old', position='append' )

      end if

    else

!     add vector point data to new file

      if ( lbinary ) then

        open ( unit=unit_gmsh, file=filename, status='replace', &
          access='stream', form='unformatted' )

        call gmsh_header_binary
        if ( .not. filedataonly ) &
                     call gmsh_nodes_elements_binary ( mesh, lgroups, lblend )

      else

        open ( unit=unit_gmsh, file=filename, status='replace' )

        call gmsh_header
        if ( .not. filedataonly ) &
                            call gmsh_nodes_elements ( mesh, lgroups, lblend )

      end if

    end if


!   write first part

    ltime = set_optional ( variable=time, default=0.0_dp )
    ltimestepnr = set_optional ( variable=timestepnr, default=0 )

    if ( lbinary ) then
      call write_dataheader_gmsh_binary ( elementwise, trim(dataname), &
        ltimestepnr, ltime, nfc=3, nelem=sum(mesh%grpnumel(lgroups)), &
        nnodes=mesh%nnodes_blend(lblend+2)-mesh%nnodes_blend(lblend+1) )
    else
      call write_dataheader_gmsh ( elementwise, trim(dataname), ltimestepnr, &
      ltime, nfc=3, nelem=sum(mesh%grpnumel(lgroups)), &
      nnodes=mesh%nnodes_blend(lblend+2)-mesh%nnodes_blend(lblend+1) )
    end if


!   write vector

    if ( present(sysvector) ) then

!     sysvector in nodes

      firstskip = .true.

      do nodenr = mesh%nnodes_blend(lblend+1)+1, mesh%nnodes_blend(lblend+2)

        if ( present(physq) ) then
          call pos_array_node ( problem, loffset+nodenr, dof, pos, &
            physqarr=[physq] )
        else
          call pos_array_node ( problem, loffset+nodenr, dof, pos )
        end if

        if ( any(deg > dof) ) then

          if ( firstskip .and. WARN_ON_NO_DATA_IN_NODES ) then
            write(*,'(/3(a/),2a/)') &
              'Warning in write_vector_gmsh:', &
              ' there are nodes where the required unknown in the ', &
              ' sysvector is undefined. The vector value is set to zero. ', &
              ' dataname = ', trim(dataname)
          end if

          firstskip = .false.

          where ( deg > dof .or. deg == 0 )
            vvalue = 0
          else where
            vvalue = sysvector%u( pos(deg1) )
          end where

        else

          where ( deg == 0 )
            vvalue = 0
          else where
            vvalue = sysvector%u( pos(deg1) )
          end where

        end if

        if ( lbinary ) then
          write ( unit=unit_gmsh ) int(nodenr,int32), vvalue
        else
          write ( unit=unit_gmsh, fmt='(i0,3es16.8)' ) &
                                           nodenr, real ( vvalue, kind=sp )
        end if

      end do

    else if ( present(vector) .and. elementwise ) then

!     vector in nodes given per element

!     start at internal element number
      elemnr = count ( mesh%points(1:mesh%npoints) > 0 ) + &
               sum ( mesh%curves(1:mesh%ncurves)%nelem ) + &
               sum ( mesh%surfaces(1:mesh%nsurfaces)%nelem )

      do elgrp = 1, mesh%nelgrp

        if ( .not. any ( lgroups == elgrp ) ) cycle

        numnod = mesh%element(elgrp)%numnod

        allocate ( v(3,numnod) )

        s = ( deg1 - 1 ) * numnod   ! start pointers

        do elem = 1, mesh%grpnumel(elgrp)

!         get element vector
          call get_vector ( mesh, problem, vector, elgrp, elem, u, &
            ndofu=ndofu, sloppy=.true. )

          if ( elem == 1 .and. mod(ndofu,numnod) /= 0 ) then
            write(*,'(/6(a/),a,i0/)') &
              'Error in write_vector_gmsh:', &
              ' number of unknowns in the vector ', &
              ' is not a multiple of the number of nodes.', &
              ' write_scalar_gmsh assumes that the vector is defined ', &
              ' in all nodes having the same number of degrees of freedom.', &
              ' Element group ', elgrp
            stop
          end if

!         fill output array with nodal values of vector

          do j = 1, numnod
            where ( deg == 0 )
              v(:,j) = 0
            else where
              v(:,j) = u( s + mesh%element(elgrp)%gmsh_index(j) )
            end where
          end do

          elemnr = elemnr + 1

          if ( lbinary ) then
            write ( unit=unit_gmsh ) int([elemnr,numnod],int32)
            write ( unit=unit_gmsh ) v
          else
            write ( unit=unit_gmsh, fmt='(i0,1x,i0)' ) elemnr, numnod
            write ( unit=unit_gmsh, fmt='(10es16.8)' ) v
          end if

        end do

        deallocate ( v )

      end do

    else

!     vector in nodes

      firstskip = .true.

      do nodenr = mesh%nnodes_blend(lblend+1)+1, mesh%nnodes_blend(lblend+2)

        dof = problem%vec_nodnumdegfd(loffset+nodenr+1,vector%vec) &
                 - problem%vec_nodnumdegfd(loffset+nodenr,vector%vec)

        if ( any(deg > dof) ) then

          if ( firstskip .and. WARN_ON_NO_DATA_IN_NODES ) then
            write(*,'(/3(a/),2a/)') &
              'Warning in write_vector_gmsh:', &
              ' there are nodes where the required unknown in the ', &
              ' vector is undefined. The vector value is set to zero. ', &
              ' dataname = ', trim(dataname)
          end if

          firstskip = .false.

          where ( deg > dof .or. deg == 0 )
            vvalue = 0
          else where
            vvalue = vector%u( &
                      problem%vec_nodnumdegfd(loffset+nodenr,vector%vec)+deg1 )
          end where

        else

          where ( deg == 0 )
            vvalue = 0
          else where
            vvalue = vector%u( &
                      problem%vec_nodnumdegfd(loffset+nodenr,vector%vec)+deg1 )
          end where

        end if

        if ( lbinary ) then
          write ( unit=unit_gmsh ) int(nodenr,int32), vvalue
        else
          write ( unit=unit_gmsh, fmt='(i0,3es16.8)' )  &
                                       nodenr, real ( vvalue, kind=sp )
        end if


      end do

    end if

    if ( elementwise ) then
      if ( lbinary ) then
        call write_line ( unit=unit_gmsh, line='' )
        call write_line ( unit=unit_gmsh, line='$EndElementNodeData' )
      else
        write ( unit=unit_gmsh, fmt='(a)' ) '$EndElementNodeData'
      end if
    else
      if ( lbinary ) then
        call write_line ( unit=unit_gmsh, line='' )
        call write_line ( unit=unit_gmsh, line='$EndNodeData' )
      else
        write ( unit=unit_gmsh, fmt='(a)' ) '$EndNodeData'
      end if
    end if

    deallocate ( lgroups )

    close ( unit=unit_gmsh )

  end subroutine write_vector_gmsh


! write a tensor field to a gmsh file

  subroutine write_tensor_gmsh ( mesh, problem, filename, dataname, physq, &
    degfd, sysvector, vector, append, dataonly, groups, symmetric, timestepnr, &
    time, binary, assume2D, assume3D, assume33, blend, offset )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the data to (it must have the .gmsh extension).
    character (len=*), intent(in) :: filename

!   the dataname for the vector, for example 'stress'.
    character (len=*), intent(in) :: dataname

!   the physical quantity to be written
    integer, intent(in), optional :: physq

!   the degrees of freedom to be written as a tensor (within the physical
!   quantity if physq is also present).
!   Setting one of the components to zero means this component is
!   set to zero.
    integer, intent(in), dimension(:,:), optional :: degfd

!   sysvector to be written
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be written
    type(vector_t), intent(in), optional :: vector

!   if .true.: append data to existing file in order to avoid writing
!   coordinate and topology data more than ones.
!   default=.false.
    logical, intent(in), optional :: append

!   if .true.: write only data to existing file (no mesh). In this way the mesh
!   and the data can be written to separate files.
!   default=.false.
    logical, intent(in), optional :: dataonly

!   if present: the element groups to be written
!   For example groups=(/2,4/) will write mesh for element groups 2 and 4.
!   Default: all groups are written
!   NOTE: the data in all the nodal points of the mesh is written and only
!   the elements are affected by choosing a subset (of all groups) to be
!   written.
    integer, dimension(:), intent(in), optional :: groups

!   if .true.: the tensor is symmetric
!   default=.true.
!   NOTE: a symmetric tensor is stored in tfem as upper-diagonal, row-wise.
!         a non-symmetric tensor is stored row-wise.
!         So:  2D symmetric: xx, xy, yy
!              2D unsymmetric: xx, xy, yx, yy
!              3D symmetric: xx, xy, xz, yy, yz, zz
!              3D unsymmetric: xx, xy, xz, yx, yy, yz, zx, zy, zz
!   In the gmsh output file a complete 3x3 tensor is written, whether it is 2D,
!   3D, symmetric or unsymmetric. The file ends up being bigger by the writing
!   of zeros for the components that are not actually there.
!   Note, that in gmsh the components are numbered 0,1,...,8 representing
!   the tensor components row wise.
    logical, intent(in), optional :: symmetric

!   time step number/index (starting with 0)
!   default=0
    integer, intent(in), optional :: timestepnr

!   the (physical) time value associated with the data
!   default=0.0
    real(dp), intent(in), optional :: time

!   if .true. write the data as a binary file.
!   default=.false.
!   The advantages of a binary file are:
!   a. the data is identical to the internal representation (no errors due
!      to truncation).
!   b. the occupied space is usually smaller.
!   c. reading the file can be somewhat faster.
!   A disadvantage is that the file cannot be interpreted by humans.
    logical, intent(in), optional :: binary

!   Assume a 2D tensor (even if mesh%ndim /= 2).
!   default=.false.
    logical, intent(in), optional :: assume2D

!   Assume a 3D tensor (even if mesh%ndim /= 3).
!   default=.false.
    logical, intent(in), optional :: assume3D

!   Assume a 3D tensor with the 33 diagonal component present
!   (only applicable for mesh%ndim = 2).
!   default=.false.
    logical, intent(in), optional :: assume33

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
!   default=0
!   Notes:
!     - blend has no effect if vector%elementwise=.true.
!     - blend should not change its value for "append" calls.
    integer, intent(in), optional :: blend

!   In case
!    1. the structure problem is created using a mesh consisting of
!       blended meshes,
!    2. the data corresponds to a blended mesh (so not to the main mesh),
!    3. a "fake" mesh is used that only contains the nodes corresponding
!       to the data
!   the value of offset should be equal to the number nodes missing from
!   the main mesh and the blended meshes "below" the current mesh. The value
!   of offest should be 0 if the "fake" mesh is based on the main mesh.
!   default = 0
!   NOTE: the value for offset is given by mesh%nnodes_blend(m+2), where
!   mesh is the blended mesh and m is the blend mesh just "below" the
!   current blend mesh (note that m=0 represents the main mesh).
!   Example: for plotting higher-order pressures in a velocity/pressure
!   Taylor-Hood formulation on its corresponding grid, offset should be
!   given the value of the number of nodes of the velocity grid.
    integer, intent(in), optional :: offset


    logical :: fileexists, fileappend, firstskip, lsymmetric, filedataonly
    logical :: agroups(mesh%nelgrp), elementwise
    logical :: lbinary, lassume2D, lassume3D, lassume33
    integer :: i, ag, elgrp, ndofu, elemnr, elem, s(3,3), numnod, j, lblend
    integer :: nodenr, pos(problem%maxnoddegfd), dof, deg(3,3), deg1(3,3)
    integer, allocatable, dimension(:) :: lgroups
    integer :: ltimestepnr, loffset

    real(dp) :: tvalue(3,3), u(problem%maxvecnoddegfd*mesh%maxelnumnod)
    real(dp), dimension(:,:,:), allocatable :: t
    real(dp) :: ltime

!   some checking

    call check ( mesh, 'write_tensor_gmsh' )
    call check ( problem, 'write_tensor_gmsh', mesh )

    lassume2D = set_optional ( variable=assume2D, default=.false. )
    lassume3D = set_optional ( variable=assume3D, default=.false. )
    lassume33 = set_optional ( variable=assume33, default=.false. )

    if ( count( [ lassume2D, lassume3D, lassume33 ] ) > 1 ) then
      write(*,'(/a/a/)') &
        'Error in write_tensor_gmsh: only one of ', &
        ' assume2D, assume3D, assume33 can be set .true.'
      stop
    end if

    if ( present(blend) ) then
      if ( blend < 0 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error write_tensor_gmsh: invalid argument ', &
          '  blend < 0 or larger than ', mesh%nblend
        stop
      end if
      lblend = blend
    else
      lblend = 0
    end if

    loffset = set_optional ( variable=offset, default=0 )

    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in write_tensor_gmsh: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(sysvector) ) then
      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in write_tensor_gmsh: sysvector has not been created '
        stop
      end if
      if ( problem%probnr /= sysvector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_tensor_gmsh: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in sysvector = ',  sysvector%probnr
        stop
      end if
      elementwise = .false.
    else if ( present(vector) ) then
      if ( .not. vector%created ) then
        write(*,'(/a/)') &
          'Error in write_tensor_gmsh: vector has not been created '
        stop
      end if
      if ( problem%probnr /= vector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_tensor_gmsh: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in vector = ', vector%probnr
        stop
      end if
      elementwise = vector%elementwise
    end if

!   check groups
    if ( present(groups) ) then
      if ( any ( groups < 1 ) .or. any ( groups > mesh%nelgrp ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: parameter groups in the heading of write_tensor_gmsh is ', &
          'out of range: some groups are < 1 or larger than ', mesh%nelgrp
        stop
      end if
      if ( problem%numinactivegroups == 0 ) then
        allocate(lgroups(size(groups)))
        lgroups = groups
      else
        agroups = .false.
        agroups(groups) = .true.
        agroups(problem%activegroups) = .true.
        allocate( lgroups(count(agroups)) )
        ag = 0
        do elgrp = 1, mesh%nelgrp
          if ( agroups(elgrp) ) then
            ag = ag + 1
            lgroups(ag) = elgrp
          end if
        end do
      end if
    else if ( problem%numinactivegroups > 0 ) then
      allocate( lgroups(size(problem%activegroups)) )
      lgroups = problem%activegroups
    else
      allocate(lgroups(mesh%nelgrp))
      lgroups = [ (i,i=1,mesh%nelgrp) ]
    end if

!   symmetry

    lsymmetric = set_optional ( variable=symmetric, default=.true. )

!   Set deg

    if ( present(degfd) ) then

      if ( any ( shape(degfd) > 3 ) ) then
        write(*,'(/a/)') &
          'Error in write_tensor_gmsh: dimension degfd larger than 3 '
        stop
      end if

      if ( any(degfd < 0 ) ) then
        write(*,'(/a/)') &
          'Error in write_tensor_gmsh: values in degfd must be not be negative'
        stop
      end if

      deg = 0
      deg(1:size(degfd,1),1:size(degfd,2)) = degfd

    else if ( mesh%ndim == 3 .or. lassume3D ) then

      if ( lsymmetric ) then
        deg(1,:) = [1,2,3]
        deg(2,:) = [2,4,5]
        deg(3,:) = [3,5,6]
      else
        deg(1,:) = [1,2,3]
        deg(2,:) = [4,5,6]
        deg(3,:) = [7,8,9]
      end if

    else if ( mesh%ndim == 2 .and. lassume33 ) then

      if ( lsymmetric ) then
        deg(1,:) = [1,2,0]
        deg(2,:) = [2,3,0]
        deg(3,:) = [0,0,4]
      else
        deg(1,:) = [1,2,0]
        deg(2,:) = [3,4,0]
        deg(3,:) = [0,0,5]
      end if

    else if ( mesh%ndim == 2 .or. lassume2D ) then

      if ( lsymmetric ) then
        deg(1,:) = [1,2,0]
        deg(2,:) = [2,3,0]
        deg(3,:) = [0,0,0]
      else
        deg(1,:) = [1,2,0]
        deg(2,:) = [3,4,0]
        deg(3,:) = [0,0,0]
      end if

    else if ( mesh%ndim == 1 ) then

      deg(1,:) = [1,0,0]
      deg(2,:) = [0,0,0]
      deg(3,:) = [0,0,0]

    end if

!   Set deg1 (to avoid zero index)

    where ( deg == 0 )
      deg1 = 1
    else where
      deg1 = deg
    end where

!   set fileappend, filedataonly

    fileappend = set_optional ( variable=append, default=.false. )
    filedataonly = set_optional ( variable=dataonly, default=.false. )

    lbinary = set_optional ( variable=binary, default=.false. )

!   open file

    inquire ( file=filename, exist=fileexists )

    if ( fileexists .and. fileappend ) then

!     add tensor point data to existing file

      if ( lbinary ) then

        open ( unit=unit_gmsh, file=filename, status='old', position='append', &
          access='stream', form='unformatted' )

      else

        open ( unit=unit_gmsh, file=filename, status='old', position='append' )

      end if

    else

!     add tensor point data to new file

      if ( lbinary ) then

        open ( unit=unit_gmsh, file=filename, status='replace', &
          access='stream', form='unformatted' )

        call gmsh_header_binary
        if ( .not. filedataonly ) &
                     call gmsh_nodes_elements_binary ( mesh, lgroups, lblend )

      else

        open ( unit=unit_gmsh, file=filename, status='replace' )

        call gmsh_header
        if ( .not. filedataonly ) &
                           call gmsh_nodes_elements ( mesh, lgroups, lblend )

      end if

    end if


!   write first part

    ltime = set_optional ( variable=time, default=0.0_dp )
    ltimestepnr = set_optional ( variable=timestepnr, default=0 )

    if ( lbinary ) then
      call write_dataheader_gmsh_binary ( elementwise, trim(dataname), &
        ltimestepnr, ltime, nfc=9, nelem=sum(mesh%grpnumel(lgroups)), &
        nnodes=mesh%nnodes_blend(lblend+2)-mesh%nnodes_blend(lblend+1) )
    else
      call write_dataheader_gmsh ( elementwise, trim(dataname), ltimestepnr, &
        ltime, nfc=9, nelem=sum(mesh%grpnumel(lgroups)), &
        nnodes=mesh%nnodes_blend(lblend+2)-mesh%nnodes_blend(lblend+1) )
    end if


!   write tensor

    if ( present(sysvector) ) then

!     sysvector in nodes

      firstskip = .true.

      do nodenr = mesh%nnodes_blend(lblend+1)+1, mesh%nnodes_blend(lblend+2)

        if ( present(physq) ) then
          call pos_array_node ( problem, loffset+nodenr, dof, pos, &
            physqarr=[physq] )
        else
          call pos_array_node ( problem, loffset+nodenr, dof, pos )
        end if

        if ( any(deg > dof) ) then

          if ( firstskip .and. WARN_ON_NO_DATA_IN_NODES ) then
            write(*,'(/3(a/),2a/)') &
              'Warning in write_tensor_gmsh:', &
              ' there are nodes where the required unknown in the ', &
              ' sysvector is undefined. The tensor value is set to zero. ', &
              ' dataname = ', trim(dataname)
          end if

          firstskip = .false.

          do i = 1, 3
            where ( deg(i,:) > dof .or. deg(i,:) == 0 )
              tvalue(i,:) = 0
            else where
              tvalue(i,:) = sysvector%u( pos(deg1(i,:)) )
            end where
          end do

        else

          do i = 1, 3
            where ( deg(i,:) == 0 )
              tvalue(i,:) = 0
            else where
              tvalue(i,:) = sysvector%u( pos(deg1(i,:)) )
            end where
          end do

        end if

        if ( lbinary ) then
          write ( unit=unit_gmsh ) int(nodenr,int32), transpose(tvalue)
        else
          write ( unit=unit_gmsh, fmt='(i0,9es16.8)' ) &
                             nodenr, real ( transpose(tvalue), kind=sp )
        end if

      end do

    else if ( present(vector) .and. elementwise ) then

!     vector in nodes given per element

!     start at internal element number
      elemnr = count ( mesh%points(1:mesh%npoints) > 0 ) + &
               sum ( mesh%curves(1:mesh%ncurves)%nelem ) + &
               sum ( mesh%surfaces(1:mesh%nsurfaces)%nelem )

      do elgrp = 1, mesh%nelgrp

        if ( .not. any ( lgroups == elgrp ) ) cycle

        numnod = mesh%element(elgrp)%numnod

        allocate ( t(3,3,numnod) )

        s = ( deg1 - 1 ) * numnod   ! start pointers

        do elem = 1, mesh%grpnumel(elgrp)

!         get element vector
          call get_vector ( mesh, problem, vector, elgrp, elem, u, &
            ndofu=ndofu, sloppy=.true. )

          if ( elem == 1 .and. mod(ndofu,numnod) /= 0 ) then
            write(*,'(/6(a/),a,i0/)') &
              'Error in write_tensor_gmsh:', &
              ' number of unknowns in the vector ', &
              ' is not a multiple of the number of nodes.', &
              ' write_scalar_gmsh assumes that the vector is defined ', &
              ' in all nodes having the same number of degrees of freedom.', &
              ' Element group ', elgrp
            stop
          end if

!         fill output array with nodal values of tensor

          do i = 1, 3
            do j = 1, 3
              if ( deg(i,j) == 0 ) then
                t(i,j,:) = 0
              else
                t(i,j,:) = u( s(i,j) + mesh%element(elgrp)%gmsh_index )
              end if
            end do
          end do

          elemnr = elemnr + 1

          if ( lbinary ) then
            write ( unit=unit_gmsh ) int([elemnr, numnod],int32)
            write ( unit=unit_gmsh ) ( transpose(t(:,:,j)), j=1,numnod )
          else
            write ( unit=unit_gmsh, fmt='(i0,1x,i0)' ) elemnr, numnod
            write ( unit=unit_gmsh, fmt='(10es16.8)' ) &
                               ( transpose(t(:,:,j)), j=1,numnod )
          end if

        end do

        deallocate ( t )

      end do

    else

!     vector in nodes

      firstskip = .true.

      do nodenr = mesh%nnodes_blend(lblend+1)+1, mesh%nnodes_blend(lblend+2)

        dof = problem%vec_nodnumdegfd(loffset+nodenr+1,vector%vec) &
                 - problem%vec_nodnumdegfd(loffset+nodenr,vector%vec)

        if ( any(deg > dof) ) then

          if ( firstskip .and. WARN_ON_NO_DATA_IN_NODES ) then
            write(*,'(/3(a/),2a/)') &
              'Warning in write_tensor_gmsh:', &
              ' there are nodes where the required unknown in the ', &
              ' vector is undefined. The tensor value is set to zero. ', &
              ' dataname = ', trim(dataname)
          end if

          firstskip = .false.

          do i = 1, 3
            where ( deg(i,:) > dof .or. deg(i,:) == 0 )
              tvalue(i,:) = 0
            else where
              tvalue(i,:) = vector%u( &
                problem%vec_nodnumdegfd(loffset+nodenr,vector%vec)+deg1(i,:) )
            end where
          end do

        else

          do i = 1, 3
            where ( deg(i,:) == 0 )
              tvalue(i,:) = 0
            else where
              tvalue(i,:) = vector%u( &
                problem%vec_nodnumdegfd(loffset+nodenr,vector%vec)+deg1(i,:) )
            end where
          end do

        end if

        if ( lbinary ) then
          write ( unit=unit_gmsh ) int(nodenr,int32), transpose(tvalue)
        else
          write ( unit=unit_gmsh, fmt='(i0,9es16.8)' ) &
                              nodenr, real ( transpose(tvalue), kind=sp )
        end if

      end do

    end if

    if ( elementwise ) then
      if ( lbinary ) then
        call write_line ( unit=unit_gmsh, line='' )
        call write_line ( unit=unit_gmsh, line='$EndElementNodeData' )
      else
        write ( unit=unit_gmsh, fmt='(a)' ) '$EndElementNodeData'
      end if
    else
      if ( lbinary ) then
        call write_line ( unit=unit_gmsh, line='' )
        call write_line ( unit=unit_gmsh, line='$EndNodeData' )
      else
        write ( unit=unit_gmsh, fmt='(a)' ) '$EndNodeData'
      end if
    end if

    deallocate ( lgroups )

    close ( unit=unit_gmsh )

  end subroutine write_tensor_gmsh


! write a scalar field to a parsed postprocessing (.pos) file

  subroutine write_scalar_gmsh_parsed ( mesh, problem, filename, dataname, &
    physq, degfd, sysvector, vector, points, curves, surfaces, groups )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the data to (it must have the .pos extension).
    character (len=*), intent(in) :: filename

!   the dataname for the scalar, for example 'pressure'.
    character (len=*), intent(in) :: dataname

!   the physical quantity to be written
    integer, intent(in), optional :: physq

!   the degree of freedom to be written (within the physical quantity if
!   physq is also present)
    integer, intent(in), optional :: degfd

!   sysvector to be written
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be written
    type(vector_t), intent(in), optional :: vector

!   if present: the points to be written
!   For example points=(/1,4/) will write data for points 1 and 4.
!   Default: no points are written
    integer, dimension(:), intent(in), optional :: points

!   if present: the curves to be written
!   For example curves=(/1,4/) will write data for curves 1 and 4.
!   Default: no curves are written
    integer, dimension(:), intent(in), optional :: curves

!   if present: the surfaces to be written
!   For example surfaces=(/1,4/) will write data for curves 1 and 4.
!   Default: no surfaces are written
    integer, dimension(:), intent(in), optional :: surfaces

!   if present: the element groups to be written
!   For example groups=(/2,4/) will write data for element groups 2 and 4.
!   Default: all groups are written
    integer, dimension(:), intent(in), optional :: groups


    logical :: firstskip
    logical :: lpoints(mesh%npoints), lcurves(mesh%ncurves)
    logical :: lsurfaces(mesh%nsurfaces), lgroups(mesh%nelgrp)
    integer :: i, elgrp, ndofu, elem, s, numnod, subelem, deg, j
    integer :: dof, point, nodenr, curve, surface
    integer :: pos(problem%maxnoddegfd)
    real(dp) :: svalue
    real(dp), allocatable, dimension(:) :: u

    call check ( mesh, 'write_scalar_gmsh_parsed' )
    call check_blend ( mesh, 'write_scalar_gmsh_parsed' )
    call check ( problem, 'write_scalar_gmsh_parsed', mesh )

    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in write_scalar_gmsh_parsed only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(sysvector) ) then
      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in write_scalar_gmsh_parsed: sysvector has not been created '
        stop
      end if
      if ( problem%probnr /= sysvector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_scalar_gmsh_parsed: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in sysvector = ',  sysvector%probnr
        stop
      end if
      allocate(u(problem%maxnoddegfd*mesh%maxelnumnod))
    else if ( present(vector) ) then
      if ( .not. vector%created ) then
        write(*,'(/a/)') &
          'Error in write_scalar_gmsh_parsed: vector has not been created '
        stop
      end if
      if ( problem%probnr /= vector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_scalar_gmsh_parsed: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in vector = ', vector%probnr
        stop
      end if
      allocate(u(problem%maxvecnoddegfd*mesh%maxelnumnod))
    end if

!   check points
    if ( present(points) ) then
      if ( any ( points < 1 ) .or. any ( points > mesh%npoints ) ) then
        write(*,'(/2a/a,i0/)') &
          'Error: parameter points in the heading of ', &
          'write_scalar_gmsh_parsed ', &
          'out of range: some points are < 1 or larger than ', mesh%npoints
        stop
      end if
      lpoints = .false.
      lpoints(points) = .true.
    else
      lpoints = .false.
    end if

!   check curves
    if ( present(curves) ) then
      if ( any ( curves < 1 ) .or. any ( curves > mesh%ncurves ) ) then
        write(*,'(/2a/a,i0/)') &
          'Error: parameter curves in the heading of ', &
          'write_scalar_gmsh_parsed ', &
          'out of range: some curves are < 1 or larger than ', mesh%ncurves
        stop
      end if
      lcurves = .false.
      lcurves(curves) = .true.
    else
      lcurves = .false.
    end if

!   check surfaces
    if ( present(surfaces) ) then
      if ( any ( surfaces < 1 ) .or. any ( surfaces > mesh%nsurfaces ) ) then
        write(*,'(/2a/a,i0/)') &
          'Error: parameter surfaces in the heading of ', &
          'write_scalar_gmsh_parsed ', &
          'out of range: some surfaces are < 1 or larger than ', mesh%nsurfaces
        stop
      end if
      lsurfaces = .false.
      lsurfaces(surfaces) = .true.
    else
      lsurfaces = .false.
    end if

!   check groups
    if ( present(groups) ) then
      if ( any ( groups < 1 ) .or. any ( groups > mesh%nelgrp ) ) then
        write(*,'(/2a/a,i0/)') &
          'Error: parameter groups in the heading of ', &
          'write_scalar_gmsh_parsed ', &
          'out of range: some groups are < 1 or larger than ', mesh%nelgrp
        stop
      end if
      if ( problem%numinactivegroups == 0 ) then
        lgroups = .false.
        lgroups(groups) = .true.
      else
        lgroups = .false.
        lgroups(groups) = .true.
        lgroups(problem%activegroups) = .true.
      end if
    else if ( problem%numinactivegroups > 0 ) then
      lgroups = .false.
      lgroups(problem%activegroups) = .true.
    else
      lgroups = .true.
    end if

!   set deg

    if ( present(degfd) ) then
      if ( degfd < 1 ) then
        write(*,'(/a/)') &
          'Error in write_scalar_gmsh_parsed degfd must be larger than zero'
        stop
      end if
      deg = degfd
    else
      deg = 1 ! take first degree
    end if

!   open file

    open ( unit=unit_gmsh, file=filename, status='replace' )

    write ( unit=unit_gmsh, fmt='(3a)' ) 'View "', trim(dataname), '" {'


!   write scalar point by point

    if ( present(sysvector) ) then

!     sysvector in nodes

      firstskip = .true.

      do point = 1, mesh%npoints

        if ( .not. lpoints(point) ) cycle

        nodenr = mesh%points(point)

        if ( present(physq) ) then
          call pos_array_node ( problem, nodenr, dof, pos, physqarr=[physq] )
        else
          call pos_array_node ( problem, nodenr, dof, pos )
        end if

        if ( deg > dof ) then

          if ( firstskip .and. WARN_ON_NO_DATA_IN_NODES ) then
            write(*,'(/3(a/),2a/)') &
              'Warning in write_scalar_gmsh_parsed:', &
              ' there are nodes where the required unknown in the ', &
              ' sysvector is undefined. The scalar value is set to zero. ', &
              ' dataname = ', trim(dataname)
          end if
          firstskip = .false.
          svalue = 0

        else

          svalue = sysvector%u( pos(deg) )

        end if

        write ( unit=unit_gmsh, fmt='(a)' ) 'SP('
        call write_list ( [ mesh%coor(nodenr,:), ( 0._dp, j=mesh%ndim+1,3 ) ] )
        write ( unit=unit_gmsh, fmt='(a)' ) '){'
        write ( unit=unit_gmsh, fmt='(es15.8)' ) svalue
        write ( unit=unit_gmsh, fmt='(a)' ) '};'

      end do

    end if


!   write scalar element by element for curves

    do curve = 1, mesh%ncurves

      if ( .not. lcurves(curve) ) cycle

      numnod = mesh%curves(curve)%element%numnod

      s = ( deg - 1 ) * numnod ! start pointer

      do elem = 1, mesh%curves(curve)%nelem

!       get element vector

        if ( present(sysvector) ) then
          if ( present(physq) ) then
            call get_sysvector_geometry ( mesh, problem, sysvector, elem, u, &
              curve, physq=[physq], order='DN' )
          else
            call get_sysvector_geometry ( mesh, problem, sysvector, elem, u, &
              curve, order='DN' )
          end if
        else if ( present(vector) ) then
          call get_vector_geometry ( mesh, problem, vector, elem, u, curve )
        end if

        do subelem = 1, mesh%curves(curve)%element%gmsh_numelm_parsed
          write ( unit=unit_gmsh, fmt='(a)' ) &
            'S'//mesh%curves(curve)%element%gmsh_elshape_parsed//'('
          call write_list ( &
            [ ( mesh%coor(mesh%curves(curve)%topology( &
           mesh%curves(curve)%element%gmsh_index_parsed(i,subelem),elem,2),:), &
                 ( 0._dp, j=mesh%ndim+1,3 ), &
                 i=1,mesh%curves(curve)%element%gmsh_numnod_parsed ) ] )
          write ( unit=unit_gmsh, fmt='(a)' ) '){'
          call write_list ( &
                u(s+mesh%curves(curve)%element%gmsh_index_parsed(:,subelem)) )
          write ( unit=unit_gmsh, fmt='(a)' ) '};'
        end do

      end do

    end do


!   write scalar element by element for surfaces

    do surface = 1, mesh%nsurfaces

      if ( .not. lsurfaces(surface) ) cycle

      numnod = mesh%surfaces(surface)%element%numnod

      s = ( deg - 1 ) * numnod ! start pointer

      do elem = 1, mesh%surfaces(surface)%nelem

!       get element vector

        if ( present(sysvector) ) then
          if ( present(physq) ) then
            call get_sysvector_geometry ( mesh, problem, sysvector, elem, u, &
              surface=surface, physq=[physq], order='DN' )
          else
            call get_sysvector_geometry ( mesh, problem, sysvector, elem, u, &
              surface=surface, order='DN' )
          end if
        else if ( present(vector) ) then
          call get_vector_geometry ( mesh, problem, vector, elem, u, &
            surface=surface )
        end if

        do subelem = 1, mesh%surfaces(surface)%element%gmsh_numelm_parsed
          write ( unit=unit_gmsh, fmt='(a)' ) &
            'S'//mesh%surfaces(surface)%element%gmsh_elshape_parsed//'('
          call write_list ( &
            [ ( mesh%coor(mesh%surfaces(surface)%topology( &
           mesh%surfaces(surface)%element%gmsh_index_parsed(i,subelem),&
                &elem,2),:), ( 0._dp, j=mesh%ndim+1,3 ), &
                 i=1,mesh%surfaces(surface)%element%gmsh_numnod_parsed ) ] )
          write ( unit=unit_gmsh, fmt='(a)' ) '){'
          call write_list ( &
             u(s+mesh%surfaces(surface)%element%gmsh_index_parsed(:,subelem)) )
          write ( unit=unit_gmsh, fmt='(a)' ) '};'
        end do

      end do

    end do


!   write scalar element by element for element groups

    do elgrp = 1, mesh%nelgrp

      if ( .not. lgroups(elgrp) ) cycle

      numnod = mesh%element(elgrp)%numnod

      s = ( deg - 1 ) * numnod ! start pointer

      do elem = 1, mesh%grpnumel(elgrp)

!       get element vector

        if ( present(sysvector) ) then
          if ( present(physq) ) then
            call get_sysvector ( mesh, problem, sysvector, elgrp, elem, u, &
              physq=[physq], order='DN', ndofu=ndofu, sloppy=.true. )
          else
            call get_sysvector ( mesh, problem, sysvector, elgrp, elem, u, &
              order='DN', ndofu=ndofu, sloppy=.true. )
          end if
        else if ( present(vector) ) then
          call get_vector ( mesh, problem, vector, elgrp, elem, u, &
            ndofu=ndofu, sloppy=.true. )
        end if

        if ( ndofu < deg * numnod ) then
          write(*,'(/2(a/),3a,i0/)') &
            'Error in write_scalar_gmsh_parsed:', &
            ' Not enough data in nodal points of element.', &
            ' Element group ', elgrp, ' dataname = ', trim(dataname)
          stop
        else if ( elem == 1 .and. mod(ndofu,numnod) /= 0 ) then
          write(*,'(/7(a/),a,i0/)') &
            'Error in write_scalar_gmsh_parsed:', &
            ' number of unknowns in the vector ', &
            ' is not a multiple of the number of nodes.', &
            ' write_scalar_gmsh_parsed assumes that ', &
            ' the (sys)vector is defined ', &
            ' in all nodes having the same number of degrees of freedom.', &
            ' Element group ', elgrp
          stop
        end if

        do subelem = 1, mesh%element(elgrp)%gmsh_numelm_parsed
          write ( unit=unit_gmsh, fmt='(a)' ) &
            'S'//mesh%element(elgrp)%gmsh_elshape_parsed//'('
          call write_list ( &
            [ ( mesh%coor(mesh%topology(elgrp)%a( &
                 mesh%element(elgrp)%gmsh_index_parsed(i,subelem),elem),:), &
                 ( 0._dp, j=mesh%ndim+1,3 ), &
                 i=1,mesh%element(elgrp)%gmsh_numnod_parsed ) ] )
          write ( unit=unit_gmsh, fmt='(a)' ) '){'
          call write_list ( &
                u(s+mesh%element(elgrp)%gmsh_index_parsed(:,subelem)) )
          write ( unit=unit_gmsh, fmt='(a)' ) '};'
        end do

      end do

    end do

    write ( unit=unit_gmsh, fmt='(a)' )  '};'

    deallocate ( u )

    close ( unit=unit_gmsh )

  contains

!   write a comma-separated list of a real array

    subroutine write_list ( a )
      real(dp), dimension(:), intent(in) :: a
      integer :: i
      write(unit=unit_gmsh,fmt='(es15.8,24(a,es15.8))') &
        a(1), ( ',', a(i), i=2,size(a) )
    end subroutine write_list

  end subroutine write_scalar_gmsh_parsed


! add a single refinement field to refinement_fields.
! Note: for compatibility with Gmsh, the refinement fields should be given in
! order of increasing npoints. This routine ensures the refinement field is
! such that the order is always in increasing npoints

  subroutine add_refinement_field ( refinement_fields, distmin, distmax, &
    dx_fine, dx_coarse, coor )

!   a new refinement field will be added to refinement_fields
    type(refinement_fields_t), intent(inout) :: refinement_fields

!   input parameters for the refinement field that will be added:
!   the elements size is set to dx_fine within a distance distmin from the
!   refinement points, and to dx_coarse outside distance distmax. In between
!   distmin and distmax, the element size is interpolated
    real(dp), intent(in) :: distmin, distmax, dx_fine, dx_coarse

!   the coordinates of the refinement points of the refinement field that will
!   be added (NOTE: number of points = size(coor,1), dimension = size(coor,2) )
!   the dimension should be either 2 or 3 for 2D or 3D meshes, respectively
    real(dp), intent(in) :: coor(:,:)

    type(refinement_fields_t) :: refinement_fields_temp
    integer :: i, j, n, ndim, npoints, size1, size2

    npoints = size(coor,1)
    ndim = size(coor,2)
    n = refinement_fields%n

!   first perform some checks

    if ( ndim /= 2 .and. ndim /=3 ) then
      write(*,'(/a/)') &
        'Error add_refinement_fields: ndim /= 2 or ndim /= 3'
      stop
    end if

    if ( n /= 0 .and. ndim /= refinement_fields%ndim ) then
      write(*,'(/a/a/)') &
        'Error: add_refinement_fields: new refinement field has different', &
        'ndim than existing refinement fields'
      stop
    end if

    if ( .not. allocated(refinement_fields%rf) .and. n /= 0 ) then
      write(*,'(/a/)') &
        'Error: add_refinement_fields: refinement_field not allocated'
      stop
    end if

    if ( allocated(refinement_fields%rf) ) then
      if ( n /= size(refinement_fields%rf) ) then
        write(*,'(/a/)') &
          'Error: add_refinement_fields: n /= size(refinement_fields%rf)'
        stop
      end if
    end if

    if ( n == 0 ) then

      allocate(refinement_fields%rf(1))

      refinement_fields%rf(1)%distmin = distmin
      refinement_fields%rf(1)%distmax = distmax
      refinement_fields%rf(1)%dx_fine = dx_fine
      refinement_fields%rf(1)%dx_coarse = dx_coarse

      allocate(refinement_fields%rf(1)%coor(npoints,ndim))

      refinement_fields%rf(1)%coor = coor
      refinement_fields%rf(1)%npoints = npoints

    else

!     check if refinement fields are in order of increasing size
      do i = 1, n - 1
        if ( refinement_fields%rf(i+1)%npoints < &
                                      refinement_fields%rf(i)%npoints ) then
          write(*,'(2(/a)/)') &
            'Error add_refinement_field: refinement fields not ordered in', &
            'increasing npoints'
          stop
        end if
      end do

!     move refinement_fields to a temporary variable
      call move_alloc ( refinement_fields%rf, refinement_fields_temp%rf )
      refinement_fields_temp%n = refinement_fields%n

      allocate(refinement_fields%rf(n+1))

!     find the index of where to insert the refinement field (to keep the
!     increasing order of npoints
      i = 1
      do j = 1, n
        if ( npoints <= refinement_fields_temp%rf(i)%npoints ) exit
        i = i + 1
      end do

!     insert the refinement field
      refinement_fields%rf(i)%distmin = distmin
      refinement_fields%rf(i)%distmax = distmax
      refinement_fields%rf(i)%dx_fine = dx_fine
      refinement_fields%rf(i)%dx_coarse = dx_coarse

      allocate(refinement_fields%rf(i)%coor(npoints,ndim))

      refinement_fields%rf(i)%coor = coor
      refinement_fields%rf(i)%npoints = npoints

!     copy the old refinement fields
      i = 1
      do j = 1, n + 1

        if ( .not. allocated(refinement_fields%rf(j)%coor) ) then

          size1=size(refinement_fields_temp%rf(i)%coor,1)
          size2=size(refinement_fields_temp%rf(i)%coor,2)
          allocate(refinement_fields%rf(j)%coor(size1,size2))
          refinement_fields%rf(j)%coor = &
            refinement_fields_temp%rf(i)%coor
          refinement_fields%rf(j)%distmin = &
            refinement_fields_temp%rf(i)%distmin
          refinement_fields%rf(j)%distmax = &
            refinement_fields_temp%rf(i)%distmax
          refinement_fields%rf(j)%dx_fine = &
            refinement_fields_temp%rf(i)%dx_fine
          refinement_fields%rf(j)%dx_coarse = &
            refinement_fields_temp%rf(i)%dx_coarse
          refinement_fields%rf(j)%npoints = &
            refinement_fields_temp%rf(i)%npoints

          i = i + 1

        end if

      end do

      call delete_refinement_fields ( refinement_fields_temp )

    end if

    refinement_fields%n = n + 1
    refinement_fields%ndim = ndim

  end subroutine add_refinement_field


! deallocate a refinement_fields_t

  subroutine delete_refinement_fields ( refinement_fields )

!   the refinement fields to be deleted
    type(refinement_fields_t), intent(inout) :: refinement_fields

    integer :: i, n

    n = refinement_fields%n

!   first perform some checks

    if ( .not. allocated(refinement_fields%rf) ) then
      write(*,'(/a/)') &
        'Error: delete_refinement_fields: refinement_field not allocated'
      stop
    end if

    if ( n /= size(refinement_fields%rf) ) then
      write(*,'(/a/)') &
        'Error: delete_refinement_fields: n /= size(refinement_fields%rf)'
      stop
    end if

    do i = 1, n
      deallocate ( refinement_fields%rf(i)%coor )
    end do
    deallocate ( refinement_fields%rf )

    refinement_fields%n = 0
    refinement_fields%ndim = 0

  end subroutine delete_refinement_fields


! write refinement field info to a file

  subroutine write_refinement_fields ( refinement_fields, filename )

!   the refinement fields
    type(refinement_fields_t), intent(in) :: refinement_fields

!   the filename for writing the refinement points (the file needs to be open)
    character (len=*), intent(in) :: filename

    integer :: i, j, k, n, ndim, iunit

    ndim = refinement_fields%ndim
    n = refinement_fields%n

!   first perform some checks

    if ( ndim /= 2 .and. ndim /= 3 ) then
      write(*,'(/a/)') &
        'Error write_refinement_fields: ndim /= 2 or ndim /= 3'
      stop
    end if

    if ( .not. allocated(refinement_fields%rf) ) then
      write(*,'(/a/)') &
        'Error: write_refinement_fields: refinement_field not allocated'
      stop
    end if

    if ( n /= size(refinement_fields%rf) ) then
      write(*,'(/a/)') &
        'Error: write_refinement_fields: n /= size(refinement_fields%rf)'
      stop
    end if

    inquire(file=filename, number=iunit)

    if ( iunit == -1 ) then
      write(*,'(/a/)') &
        'Error write_refinement_fields: file not open for writing'
      stop
    end if

!   write the refinement fields

    write ( iunit, '(1X,A,I0,A)' ) 'nrefinement_field = ', n, ';'

    k = 1
    do j = 1, n

      write ( iunit, '(1X,A,I0,A,I0,A)' ) 'nrefine[',j,'] = ', &
        refinement_fields%rf(j)%npoints, ';'

      if ( ndim == 2 ) then

        do i = 1,refinement_fields%rf(j)%npoints
          write ( iunit, '(1X,A,I0,A,e19.12,A)' ) 'xp_refine[',k,'] = ', &
            refinement_fields%rf(j)%coor(i,1), ';'
          write ( iunit, '(1X,A,I0,A,e19.12,A)' ) 'yp_refine[',k,'] = ', &
            refinement_fields%rf(j)%coor(i,2), ';'
          write ( iunit, '(1X,A,I0,A,e19.12,A)' ) 'zp_refine[',k,'] = ', &
            0._dp, ';'
          k=k+1
        end do

      else

        do i = 1,refinement_fields%rf(j)%npoints
          write ( iunit, '(1X,A,I0,A,e19.12,A)' ) 'xp_refine[',k,'] = ', &
            refinement_fields%rf(j)%coor(i,1), ';'
          write ( iunit, '(1X,A,I0,A,e19.12,A)' ) 'yp_refine[',k,'] = ', &
            refinement_fields%rf(j)%coor(i,2), ';'
          write ( iunit, '(1X,A,I0,A,e19.12,A)' ) 'zp_refine[',k,'] = ', &
            refinement_fields%rf(j)%coor(i,3), ';'
          k=k+1
        end do

      end if

      write ( iunit, '(1X,A,I0,A,e19.12,A)' ) 'distmin[',j,'] = ', &
        refinement_fields%rf(j)%distmin, ';'
      write ( iunit, '(1X,A,I0,A,e19.12,A)' ) 'distmax[',j,'] = ', &
        refinement_fields%rf(j)%distmax, ';'

      write ( iunit, '(1X,A,I0,A,e19.12,A)' ) 'dx_fine[',j,'] = ', &
        refinement_fields%rf(j)%dx_fine, ';'
      write ( iunit, '(1X,A,I0,A,e19.12,A)' ) 'dx_coarse[',j,'] = ', &
        refinement_fields%rf(j)% dx_coarse, ';'

    end do

   end subroutine write_refinement_fields

end module gmsh_utils_m
