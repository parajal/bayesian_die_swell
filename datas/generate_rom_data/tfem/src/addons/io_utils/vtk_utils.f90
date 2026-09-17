
! Copyright (C) 2005-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines routines for output to the Visualization Toolkit (VTK) in the legacy
! format (extension .vtk).

module vtk_utils_m

  use kind_defs_m
  use mesh_m
  use problem_defs_m
  use system_vector_m
  use vector_m
  use streamio_m
  use, intrinsic :: iso_fortran_env
  use limits_m, only: WARN_ON_NO_DATA_IN_NODES, WRITE_BINARY_VTK, &
                      WRITE_DOUBLE_VTK

  implicit none

  integer :: unit_vtk = 13 ! unit number for writing

contains


! write a mesh only to a legacy VTK file

  subroutine write_mesh_vtk ( mesh, filename, groups, double, binary, blend )

    type(mesh_t), intent(in) :: mesh

!   the filename for writing the mesh to (it must have the .vtk extension).
    character (len=*), intent(in) :: filename

!   if present: the element groups to be written
!   For example groups=(/2,4/) will write mesh for element groups 2 and 4.
!   Default: all groups are written
!   NOTE: all the nodal points of the mesh are written and only the cells
!   are affected by choosing a subset (of all groups) to be written.
    integer, dimension(:), intent(in), optional :: groups

!   write data in the "double" format (.true.) or the "float" format (.false.).
!   NOTE: files are 40-50 % larger for the double format.
!   default=WRITE_DOUBLE_VTK (set in limits_m).
    logical, intent(in), optional :: double

!   write data in the binary format (.true.) or ASCII format (.false.)
!   default=WRITE_BINARY_VTK (set in limits_m).
!   The advantages of a binary file are:
!   a. the data is identical to the internal representation.
!   b. the files are significantly smaller than ASCII vtk-files.
!   c. reading the file might be somewhat faster.
!   A disadvantage is that the file cannot be interpreted by humans.
    logical, intent(in), optional :: binary

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
!   default=0
    integer, intent(in), optional :: blend

    logical :: ldouble, lbinary
    integer :: i, ios, lblend
    integer, allocatable, dimension(:) :: lgroups

    call check ( mesh, 'write_mesh_vtk' )

    ldouble = set_optional ( variable=double, default=WRITE_DOUBLE_VTK )
    lbinary = set_optional ( variable=binary, default=WRITE_BINARY_VTK )

    if ( present(blend) ) then
      if ( blend < 0 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error write_mesh_vtk: invalid argument ', &
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
          'Error: parameter groups in the heading of write_mesh_vtk is ', &
          'out of range: some groups are < 1 or larger than ', mesh%nelgrp
        stop
      end if
      allocate(lgroups(size(groups)))
      lgroups = groups
    else
      allocate(lgroups(mesh%nelgrp))
      lgroups = [ (i,i=1,mesh%nelgrp) ]
    end if

    if ( lbinary ) then

!     write binary file

      open ( unit=unit_vtk, file=filename, access='stream', iostat=ios, &
             form='unformatted', status='replace', convert='big_endian' )

      if ( ios /= 0 ) then
        write(*,'(/a,i0/)') 'Error in write_mesh_vtk: ios = ', ios
        stop
      end if

      call vtk_header_binary
      call vtk_geometry_ug_binary ( mesh, lgroups, ldouble, lblend )

    else

!     write ASCII file

      open ( unit=unit_vtk, file=filename, status='replace', iostat=ios )

      if ( ios /= 0 ) then
        write(*,'(/a,i0/)') 'Error in write_mesh_vtk: ios = ', ios
        stop
      end if

      call vtk_header
      call vtk_geometry_ug ( mesh, lgroups, ldouble, lblend )

    end if

    deallocate(lgroups)

    close ( unit=unit_vtk )

  end subroutine write_mesh_vtk


! write a scalar field to a legacy VTK file

  subroutine write_scalar_vtk ( mesh, problem, filename, dataname, physq, &
    degfd, sysvector, vector, append, switch, groups, double, binary, &
    blend, offset )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the data to (it must have the .vtk extension).
    character (len=*), intent(in) :: filename

!   the dataname for the scalar, for example 'pressure'.
    character (len=*), intent(in) :: dataname

!   the physical quantity to be plotted
    integer, intent(in), optional :: physq

!   the degree of freedom to be plotted (within the physical quantity if
!   physq is also present)
    integer, intent(in), optional :: degfd

!   sysvector to be plotted
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be plotted
    type(vector_t), intent(in), optional :: vector

!   if .true.: append data to existing file in order to avoid writing
!   coordinate and topology data more than ones.
!   default=.false.
    logical, intent(in), optional :: append

!   The switch=.true. parameter is only needed if both point and cell data
!   formats need to written to the same file. It needs to be set when
!   "switching" formats while appending new data (append=.true.).
!   NOTE: switch can only be used once for a single data file.
!   default=.false.
    logical, intent(in), optional :: switch

!   if present: the element groups to be written
!   For example groups=(/2,4/) will write data for element groups 2 and 4.
!   Default: all groups are written
!   NOTE: the data in all the nodal points of the mesh is written and only
!   the cells are affected by choosing a subset (of all groups) to be written.
!   NOTE: groups needs to be same between different "append" calls.
    integer, dimension(:), intent(in), optional :: groups

!   write data in the "double" format (.true.) or the "float" format (.false.).
!   NOTE: files are 40-50 % larger for the double format.
!   default=WRITE_DOUBLE_VTK (set in limits_m).
    logical, intent(in), optional :: double

!   write data in the binary format (.true.) or ASCII format (.false.)
!   default=WRITE_BINARY_VTK (set in limits_m).
!   The advantages of a binary file are:
!   a. the data is identical to the internal representation.
!   b. the files are significantly smaller than ASCII vtk-files.
!   c. reading the file might be somewhat faster.
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


    character (len=80) :: line
    logical :: fileexists, fileappend, firstskip, elementwise, ldouble, lbinary
    logical :: agroups(mesh%nelgrp), switchformat
    integer :: i, ag, elgrp, elem, ndofu
    integer :: nodenr, pos(problem%maxnoddegfd), dof, deg, lblend, loffset
    integer, allocatable, dimension(:) :: lgroups

    real(dp) :: svalue, u(problem%maxvecnoddegfd*mesh%maxelnumnod)


    call check ( mesh, 'write_scalar_vtk' )
    call check ( problem, 'write_scalar_vtk', mesh )

    ldouble = set_optional ( variable=double, default=WRITE_DOUBLE_VTK )
    lbinary = set_optional ( variable=binary, default=WRITE_BINARY_VTK )
    loffset = set_optional ( variable=offset, default=0 )

    if ( present(blend) ) then
      if ( blend < 0 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error write_scalar_vtk: invalid argument ', &
          '  blend < 0 or larger than ', mesh%nblend
        stop
      end if
      lblend = blend
    else
      lblend = 0
    end if

    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in write_scalar_vtk: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(sysvector) ) then
      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in write_scalar_vtk: sysvector has not been created '
        stop
      end if
      if ( problem%probnr /= sysvector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_scalar_vtk: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in sysvector = ',  sysvector%probnr
        stop
      end if
      elementwise = .false.
    else if ( present(vector) ) then
      if ( .not. vector%created ) then
        write(*,'(/a/)') &
          'Error in write_scalar_vtk: vector has not been created '
        stop
      end if
      if ( problem%probnr /= vector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_scalar_vtk: ',&
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
          'Error: parameter groups in the heading of write_scalar_vtk is ', &
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
          'Error in write_scalar_vtk: degfd must be larger than zero'
        stop
      end if
      deg = degfd
    else
      deg = 1 ! take first degree
    end if

!   set fileappend and swithformat

    fileappend = set_optional ( variable=append, default=.false. )
    switchformat = set_optional ( variable=switch, default=.false. )

!   open file

    inquire ( file=filename, exist=fileexists )

    if ( lbinary ) then

!     binary format

      if ( fileexists .and. fileappend ) then

!       add scalar point data to existing file

        open ( unit=unit_vtk, file=filename, access='stream', status='old', &
               position='append', form='unformatted', convert='big_endian' )

!       first line of the point or cell data

        if ( switchformat ) &
            call vtk_data_format_binary ( mesh, elementwise, lgroups, lblend )

      else

        open ( unit=unit_vtk, file=filename, access='stream', &
               form='unformatted', status='replace', convert='big_endian' )

        call vtk_header_binary
        call vtk_geometry_ug_binary ( mesh, lgroups, ldouble, lblend )

!       first line of the point or cell data

        call vtk_data_format_binary ( mesh, elementwise, lgroups, lblend )

      end if

!     write first part

      if ( ldouble ) then
        write ( unit=line, fmt='(3a)' ) 'SCALARS ', trim(dataname), ' double 1'
        call write_line ( unit=unit_vtk, line=trim(line) )
      else
        write ( unit=line, fmt='(3a)' ) 'SCALARS ', trim(dataname), ' float 1'
        call write_line ( unit=unit_vtk, line=trim(line) )
      end if
      write ( unit=line, fmt='(a)' ) 'LOOKUP_TABLE default'
      call write_line ( unit=unit_vtk, line=trim(line) )

    else

!     ASCII format

      if ( fileexists .and. fileappend ) then

!       add scalar point data to existing file

        open ( unit=unit_vtk, file=filename, status='old', position='append' )

!       first line of the point or cell data

        if ( switchformat ) &
            call vtk_data_format ( mesh, elementwise, lgroups, lblend )

      else

        open ( unit=unit_vtk, file=filename, status='replace' )

        call vtk_header
        call vtk_geometry_ug ( mesh, lgroups, ldouble, lblend )

!       first line of the point or cell data

        call vtk_data_format ( mesh, elementwise, lgroups, lblend )

      end if

!     write first part

      if ( ldouble ) then
        write ( unit=unit_vtk, fmt='(3a)' ) &
                            'SCALARS ', trim(dataname), ' double 1'
      else
        write ( unit=unit_vtk, fmt='(3a)' ) &
                            'SCALARS ', trim(dataname), ' float 1'
      end if
      write ( unit=unit_vtk, fmt='(a)' ) 'LOOKUP_TABLE default'

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
              'Warning in write_scalar_vtk:', &
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
          if ( ldouble ) then
            write ( unit=unit_vtk ) svalue
          else
            write ( unit=unit_vtk ) real ( svalue, kind=sp )
          end if
        else
          if ( ldouble ) then
            write ( unit=unit_vtk, fmt='(es25.16)' ) svalue
          else
            write ( unit=unit_vtk, fmt='(es16.8)' ) real ( svalue, kind=sp )
          end if
        end if

      end do

    else if ( present(vector) .and. elementwise ) then

!     vector in nodes given per element

      do elgrp = 1, mesh%nelgrp

        if ( .not. any ( lgroups == elgrp ) ) cycle

        do elem = 1, mesh%grpnumel(elgrp)

!         get element vector
          call get_vector ( mesh, problem, vector, elgrp, elem, u, &
            order='ND', ndofu=ndofu, sloppy=.true. )

          if ( ndofu < deg ) then
            write(*,'(/2(a/),3a,i0/)') &
              'Error in write_scalar_vtk:', &
              ' Not enough data in element.', &
              ' Element group ', elgrp, ' dataname = ', trim(dataname)
            stop
          end if

          if ( lbinary ) then
            if ( ldouble ) then
              write ( unit=unit_vtk ) u(deg)
            else
              write ( unit=unit_vtk ) real ( u(deg), kind=sp )
            end if
          else
            if ( ldouble ) then
              write ( unit=unit_vtk, fmt='(10es25.16)' ) u(deg)
            else
              write ( unit=unit_vtk, fmt='(10es16.8)' ) real ( u(deg), kind=sp )
            end if
          end if

        end do

      end do

    else if ( present(vector) ) then

!     vector in nodes

      firstskip = .true.

      do nodenr = mesh%nnodes_blend(lblend+1)+1, mesh%nnodes_blend(lblend+2)

        dof = problem%vec_nodnumdegfd(loffset+nodenr+1,vector%vec) &
                 - problem%vec_nodnumdegfd(loffset+nodenr,vector%vec)

        if ( deg > dof ) then

          if ( firstskip .and. WARN_ON_NO_DATA_IN_NODES ) then
            write(*,'(/3(a/),2a/)') &
              'Warning in write_scalar_vtk:', &
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
          if ( ldouble ) then
            write ( unit=unit_vtk ) svalue
          else
            write ( unit=unit_vtk ) real ( svalue, kind=sp )
          end if
        else
          if ( ldouble ) then
            write ( unit=unit_vtk, fmt='(es25.16)' ) svalue
          else
            write ( unit=unit_vtk, fmt='(es16.8)' ) real ( svalue, kind=sp )
          end if
        end if

      end do

    end if

    deallocate ( lgroups )

    close ( unit=unit_vtk )

  end subroutine write_scalar_vtk


! write a vector field to a legacy VTK file

  subroutine write_vector_vtk ( mesh, problem, filename, dataname, physq, &
    degfd, sysvector, vector, append, switch, groups, double, binary, &
    assume2D, assume3D, blend, offset )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the data to (it must have the .vtk extension).
    character (len=*), intent(in) :: filename

!   the dataname for the vector, for example 'velocity'.
    character (len=*), intent(in) :: dataname

!   the physical quantity to be plotted
    integer, intent(in), optional :: physq

!   the degrees of freedom to be plotted as a vector (within the physical
!   quantity if physq is also present).
!   Setting one of the components to zero means this component is
!   set to zero.
    integer, intent(in), dimension(:), optional :: degfd

!   sysvector to be plotted
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be plotted
    type(vector_t), intent(in), optional :: vector

!   if .true.: append data to existing file in order to avoid writing
!   coordinate and topology data more than ones.
!   default=.false.
    logical, intent(in), optional :: append

!   The switch=.true. parameter is only needed if both point and cell data
!   formats need to written to the same file. It needs to be set when
!   "switching" formats while appending new data (append=.true.).
!   NOTE: switch can only be used once for a single data file.
!   default=.false.
    logical, intent(in), optional :: switch

!   if present: the element groups to be written
!   For example groups=(/2,4/) will write data for element groups 2 and 4.
!   Default: all groups are written
!   NOTE: the data in all the nodal points of the mesh is written and only
!   the cells are affected by choosing a subset (of all groups) to be written.
!   NOTE: groups needs to be same between different "append" calls.
    integer, dimension(:), intent(in), optional :: groups

!   write data in the "double" format (.true.) or the "float" format (.false.).
!   NOTE: files are 40-50 % larger for the double format.
!   default=WRITE_DOUBLE_VTK (set in limits_m).
    logical, intent(in), optional :: double

!   write data in the binary format (.true.) or ASCII format (.false.)
!   default=WRITE_BINARY_VTK (set in limits_m).
!   The advantages of a binary file are:
!   a. the data is identical to the internal representation.
!   b. the files are significantly smaller than ASCII vtk-files.
!   c. reading the file might be somewhat faster.
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


    character (len=80) :: line
    logical :: fileexists, fileappend, firstskip, elementwise, ldouble, lbinary
    logical :: agroups(mesh%nelgrp), switchformat, lassume2D, lassume3D
    integer :: i, ag, elgrp, elem, ndofu, lblend, loffset
    integer :: nodenr, pos(problem%maxnoddegfd), dof, deg(3), deg1(3)
    integer, allocatable, dimension(:) :: lgroups

    real(dp) :: vvalue(3), u(problem%maxvecnoddegfd*mesh%maxelnumnod)

!   some checking

    call check ( mesh, 'write_vector_vtk' )
    call check ( problem, 'write_vector_vtk', mesh )

    ldouble = set_optional ( variable=double, default=WRITE_DOUBLE_VTK )
    lbinary = set_optional ( variable=binary, default=WRITE_BINARY_VTK )
    lassume2D = set_optional ( variable=assume2D, default=.false. )
    lassume3D = set_optional ( variable=assume3D, default=.false. )
    loffset = set_optional ( variable=offset, default=0 )

    if ( lassume2D .and. lassume3D ) then
      write(*,'(/a/a/)') &
        'Error in write_vector_vtk: only one of assume2D and assume3D ', &
        ' can be set .true.'
      stop
    end if

    if ( present(blend) ) then
      if ( blend < 0 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error write_vector_vtk: invalid argument ', &
          '  blend < 0 or larger than ', mesh%nblend
        stop
      end if
      lblend = blend
    else
      lblend = 0
    end if

    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in write_vector_vtk: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(sysvector) ) then
      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in write_vector_vtk: sysvector has not been created '
        stop
      end if
      if ( problem%probnr /= sysvector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_vector_vtk: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in sysvector = ',  sysvector%probnr
        stop
      end if
      elementwise = .false.
    else if ( present(vector) ) then
      if ( .not. vector%created ) then
        write(*,'(/a/)') &
          'Error in write_vector_vtk: vector has not been created '
        stop
      end if
      if ( problem%probnr /= vector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_vector_vtk: ',&
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
          'Error: parameter groups in the heading of write_vector_vtk is ', &
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
          'Error in write_vector_vtk: dimension degfd larger than 3.'
        stop
      end if

      if ( any(degfd < 0 ) ) then
        write(*,'(/a/)') &
          'Error in write_vector_vtk: values in degfd must be not be negative'
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

!   set fileappend and switchformat

    fileappend = set_optional ( variable=append, default=.false. )
    switchformat = set_optional ( variable=switch, default=.false. )

!   open file

    inquire ( file=filename, exist=fileexists )

    if ( lbinary ) then

!     binary format

      if ( fileexists .and. fileappend ) then

!       add scalar point data to existing file

        open ( unit=unit_vtk, file=filename, access='stream', status='old', &
               position='append', form='unformatted', convert='big_endian' )

!       first line of the point or cell data

        if ( switchformat ) &
             call vtk_data_format_binary ( mesh, elementwise, lgroups, lblend )

      else

        open ( unit=unit_vtk, file=filename, access='stream', &
               form='unformatted', status='replace', convert='big_endian' )

        call vtk_header_binary
        call vtk_geometry_ug_binary ( mesh, lgroups, ldouble, lblend )

!       first line of the point or cell data

        call vtk_data_format_binary ( mesh, elementwise, lgroups, lblend )

      end if

!     write first part

      if ( ldouble ) then
        write ( unit=line, fmt='(3a)' ) 'VECTORS ', trim(dataname), ' double'
        call write_line ( unit=unit_vtk, line=trim(line) )
      else
        write ( unit=line, fmt='(3a)' ) 'VECTORS ', trim(dataname), ' float'
        call write_line ( unit=unit_vtk, line=trim(line) )
      end if

    else

!     ASCII format

      if ( fileexists .and. fileappend ) then

!       add scalar point data to existing file

        open ( unit=unit_vtk, file=filename, status='old', position='append' )

!       first line of the point or cell data

        if ( switchformat ) &
                   call vtk_data_format ( mesh, elementwise, lgroups, lblend )

      else

        open ( unit=unit_vtk, file=filename, status='replace' )

        call vtk_header
        call vtk_geometry_ug ( mesh, lgroups, ldouble, lblend )

!       first line of the point or cell data

        call vtk_data_format ( mesh, elementwise, lgroups, lblend )

      end if

!     write first part

      if ( ldouble ) then
        write ( unit=unit_vtk, fmt='(3a)' ) &
                                 'VECTORS ', trim(dataname), ' double'
      else
        write ( unit=unit_vtk, fmt='(3a)' ) &
                                 'VECTORS ', trim(dataname), ' float'
      end if

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
              'Warning in write_vector_vtk:', &
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
          if ( ldouble ) then
            write ( unit=unit_vtk ) vvalue
          else
            write ( unit=unit_vtk ) real ( vvalue, kind=sp )
          end if
        else
          if ( ldouble ) then
            write ( unit=unit_vtk, fmt='(3es25.16)' ) vvalue
          else
            write ( unit=unit_vtk, fmt='(3es16.8)' ) real ( vvalue, kind=sp )
          end if
        end if

      end do

    else if ( present(vector) .and. elementwise ) then

!     vector in nodes given per element

      do elgrp = 1, mesh%nelgrp

        if ( .not. any ( lgroups == elgrp ) ) cycle

        do elem = 1, mesh%grpnumel(elgrp)

!         get element vector
          call get_vector ( mesh, problem, vector, elgrp, elem, u, &
            order='ND', ndofu=ndofu, sloppy=.true. )

          if ( any( deg > ndofu ) ) then
            write(*,'(/2(a/),3a,i0/)') &
              'Error in write_vector_vtk:', &
              ' Not enough data in element.', &
              ' Element group ', elgrp, ' dataname = ', trim(dataname)
            stop
          end if

          where ( deg == 0 )
            vvalue = 0
          else where
            vvalue = u(deg1)
          end where

          if ( lbinary ) then
            if ( ldouble ) then
              write ( unit=unit_vtk ) vvalue
            else
              write ( unit=unit_vtk ) real ( vvalue, kind=sp )
            end if
          else
            if ( ldouble ) then
              write ( unit=unit_vtk, fmt='(3es25.16)' ) vvalue
            else
              write ( unit=unit_vtk, fmt='(3es16.8)' ) real ( vvalue, kind=sp )
            end if
          end if

        end do

      end do

    else if ( present(vector) ) then

!     vector in nodes

      firstskip = .true.

      do nodenr = mesh%nnodes_blend(lblend+1)+1, mesh%nnodes_blend(lblend+2)

        dof = problem%vec_nodnumdegfd(loffset+nodenr+1,vector%vec) &
                 - problem%vec_nodnumdegfd(loffset+nodenr,vector%vec)

        if ( any(deg > dof) ) then

          if ( firstskip .and. WARN_ON_NO_DATA_IN_NODES ) then
            write(*,'(/3(a/),2a/)') &
              'Warning in write_vector_vtk:', &
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
          if ( ldouble ) then
            write ( unit=unit_vtk ) vvalue
          else
            write ( unit=unit_vtk ) real ( vvalue, kind=sp )
          end if
        else
          if ( ldouble ) then
            write ( unit=unit_vtk, fmt='(3es25.16)' ) vvalue
          else
            write ( unit=unit_vtk, fmt='(3es16.8)' ) real ( vvalue, kind=sp )
          end if
        end if

      end do

    end if

    deallocate ( lgroups )

    close ( unit=unit_vtk )

  end subroutine write_vector_vtk


! write a tensor field to a legacy VTK file

  subroutine write_tensor_vtk ( mesh, problem, filename, dataname, physq, &
    degfd, sysvector, vector, append, switch, groups, symmetric, double, &
    binary, assume2D, assume3D, assume33, blend, offset )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the data to (it must have the .vtk extension).
    character (len=*), intent(in) :: filename

!   the dataname for the vector, for example 'stress'.
    character (len=*), intent(in) :: dataname

!   the physical quantity to be plotted
    integer, intent(in), optional :: physq

!   the degrees of freedom to be plotted as a tensor (within the physical
!   quantity if physq is also present).
!   Setting one of the components to zero means this component is
!   set to zero.
    integer, intent(in), dimension(:,:), optional :: degfd

!   sysvector to be plotted
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be plotted
    type(vector_t), intent(in), optional :: vector

!   if .true.: append data to existing file in order to avoid writing
!   coordinate and topology data more than ones.
!   default=.false.
    logical, intent(in), optional :: append

!   The switch=.true. parameter is only needed if both point and cell data
!   formats need to written to the same file. It needs to be set when
!   "switching" formats while appending new data (append=.true.).
!   NOTE: switch can only be used once for a single data file.
!   default=.false.
    logical, intent(in), optional :: switch

!   if present: the element groups to be written
!   For example groups=(/2,4/) will write mesh for element groups 2 and 4.
!   Default: all groups are written
!   NOTE: the data in all the nodal points of the mesh is written and only
!   the cells are affected by choosing a subset of the groups to be written.
!   NOTE: groups needs to be same between different "append" calls.
    integer, dimension(:), intent(in), optional :: groups

!   if .true.: the tensor is symmetric
!   default=.true.
!   NOTE: a symmetric tensor is stored in tfem as upper-diagonal, row-wise.
!         a non-symmetric tensor is stored row-wise.
!         So:  2D symmetric: xx, xy, yy
!              2D unsymmetric: xx, xy, yx, yy
!              3D symmetric: xx, xy, xz, yy, yz, zz
!              3D unsymmetric: xx, xy, xz, yx, yy, yz, zx, zy, zz
!   In the vtk output file a complete 3x3 tensor is written, whether it is 2D,
!   3D, symmetric or unsymmetric. The file ends up being bigger by the writing
!   of zeros for the components that are not actually there.
!   Note, that in paraview the components are numbered 0,1,...,8 representing
!   the tensor components row wise and can be accessed using the calculator.
    logical, intent(in), optional :: symmetric

!   write data in the "double" format (.true.) or the "float" format (.false.).
!   NOTE: files are 40-50 % larger for the double format.
!   default=WRITE_DOUBLE_VTK (set in limits_m).
    logical, intent(in), optional :: double

!   write data in the binary format (.true.) or ASCII format (.false.)
!   default=WRITE_BINARY_VTK (set in limits_m).
!   The advantages of a binary file are:
!   a. the data is identical to the internal representation.
!   b. the files are significantly smaller than ASCII vtk-files.
!   c. reading the file might be somewhat faster.
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


    character (len=80) :: line
    logical :: ldouble, lbinary
    logical :: fileexists, fileappend, firstskip, lsymmetric, elementwise
    logical :: agroups(mesh%nelgrp), switchformat
    logical :: lassume2D, lassume3D, lassume33
    integer :: i, ag, elgrp, elem, ndofu, lblend, loffset
    integer :: nodenr, pos(problem%maxnoddegfd), dof, deg(3,3), deg1(3,3)
    integer, allocatable, dimension(:) :: lgroups

    real(dp) :: tvalue(3,3), u(problem%maxvecnoddegfd*mesh%maxelnumnod)

!   some checking

    call check ( mesh, 'write_tensor_vtk' )
    call check ( problem, 'write_tensor_vtk', mesh )

    ldouble = set_optional ( variable=double, default=WRITE_DOUBLE_VTK )
    lbinary = set_optional ( variable=binary, default=WRITE_BINARY_VTK )
    lassume2D = set_optional ( variable=assume2D, default=.false. )
    lassume3D = set_optional ( variable=assume3D, default=.false. )
    lassume33 = set_optional ( variable=assume33, default=.false. )
    loffset = set_optional ( variable=offset, default=0 )

    if ( count( [ lassume2D, lassume3D, lassume33 ] ) > 1 ) then
      write(*,'(/a/a/)') &
        'Error in write_tensor_vtk: only one of ', &
        ' assume2D, assume3D, assume33 can be set .true.'
      stop
    end if

    if ( present(blend) ) then
      if ( blend < 0 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error write_tensor_vtk: invalid argument ', &
          '  blend < 0 or larger than ', mesh%nblend
        stop
      end if
      lblend = blend
    else
      lblend = 0
    end if

    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in write_tensor_vtk: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(sysvector) ) then
      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in write_tensor_vtk: sysvector has not been created '
        stop
      end if
      if ( problem%probnr /= sysvector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_tensor_vtk: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in sysvector = ',  sysvector%probnr
        stop
      end if
      elementwise = .false.
    else if ( present(vector) ) then
      if ( .not. vector%created ) then
        write(*,'(/a/)') &
          'Error in write_tensor_vtk: vector has not been created '
        stop
      end if
      if ( problem%probnr /= vector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_tensor_vtk: ',&
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
          'Error: parameter groups in the heading of write_tensor_vtk is ', &
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

    if ( present(symmetric) ) then
      lsymmetric = symmetric
    else
      lsymmetric = .true.
    end if

!   Set deg

    if ( present(degfd) ) then

      if ( any ( shape(degfd) > 3 ) ) then
        write(*,'(/a/)') &
          'Error in write_tensor_vtk: dimension degfd larger than 3.'
        stop
      end if

      if ( any(degfd < 0 ) ) then
        write(*,'(/a/)') &
          'Error in write_tensor_vtk: values in degfd must be not be negative'
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

!   set fileappend and switchformat

    fileappend = set_optional ( variable=append, default=.false. )
    switchformat = set_optional ( variable=switch, default=.false. )

!   open file

    inquire ( file=filename, exist=fileexists )

    if ( lbinary ) then

!     binary format

      if ( fileexists .and. fileappend ) then

!       add scalar point data to existing file

        open ( unit=unit_vtk, file=filename, access='stream', status='old', &
               position='append', form='unformatted', convert='big_endian' )

!       first line of the point or cell data

        if ( switchformat ) &
             call vtk_data_format_binary ( mesh, elementwise, lgroups, lblend )

      else

        open ( unit=unit_vtk, file=filename, access='stream', &
               form='unformatted', status='replace', convert='big_endian' )

        call vtk_header_binary
        call vtk_geometry_ug_binary ( mesh, lgroups, ldouble, lblend )

!       first line of the point or cell data

        call vtk_data_format_binary ( mesh, elementwise, lgroups, lblend )

      end if

!     write first part

      if ( ldouble ) then
        write ( unit=line, fmt='(3a)' ) 'TENSORS ', trim(dataname), ' double'
        call write_line ( unit=unit_vtk, line=trim(line) )
      else
        write ( unit=line, fmt='(3a)' ) 'TENSORS ', trim(dataname), ' float'
        call write_line ( unit=unit_vtk, line=trim(line) )
      end if

    else

!     ASCII format

      if ( fileexists .and. fileappend ) then

!       add scalar point data to existing file

        open ( unit=unit_vtk, file=filename, status='old', position='append' )

!       first line of the point or cell data

        if ( switchformat ) &
                 call vtk_data_format ( mesh, elementwise, lgroups, lblend )

      else

        open ( unit=unit_vtk, file=filename, status='replace' )

        call vtk_header
        call vtk_geometry_ug ( mesh, lgroups, ldouble, lblend )

!       first line of the point or cell data

        call vtk_data_format ( mesh, elementwise, lgroups, lblend )

      end if

!     write first part

      if ( ldouble ) then
        write ( unit=unit_vtk, fmt='(3a)' ) &
                                 'TENSORS ', trim(dataname), ' double'
      else
        write ( unit=unit_vtk, fmt='(3a)' ) &
                                 'TENSORS ', trim(dataname), ' float'
      end if

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
              'Warning in write_tensor_vtk:', &
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
          if ( ldouble ) then
            write ( unit=unit_vtk ) transpose(tvalue)
          else
            write ( unit=unit_vtk ) real ( transpose(tvalue), kind=sp )
          end if
        else
          if ( ldouble ) then
            write ( unit=unit_vtk, fmt='(9es25.16)' ) transpose(tvalue)
          else
            write ( unit=unit_vtk, fmt='(9es16.8)' ) &
                                       real ( transpose(tvalue), kind=sp )
          end if
        end if

      end do

    else if ( present(vector) .and. elementwise ) then

!     vector in nodes given per element

      do elgrp = 1, mesh%nelgrp

        if ( .not. any ( lgroups == elgrp ) ) cycle

        do elem = 1, mesh%grpnumel(elgrp)

!         get element vector
          call get_vector ( mesh, problem, vector, elgrp, elem, u, &
            order='ND', ndofu=ndofu, sloppy=.true. )

          if ( any( deg > ndofu ) ) then
            write(*,'(/2(a/),3a,i0/)') &
              'Error in write_tensor_vtk:', &
              ' Not enough data in element.', &
              ' Element group ', elgrp, ' dataname = ', trim(dataname)
            stop
          end if

          do i = 1, 3
            where ( deg(i,:) == 0 )
              tvalue(i,:) = 0
            else where
              tvalue(i,:) = u(deg1(i,:))
            end where
          end do

          if ( lbinary ) then
            if ( ldouble ) then
              write ( unit=unit_vtk ) transpose(tvalue)
            else
              write ( unit=unit_vtk ) real ( transpose(tvalue), kind=sp )
            end if
          else
            if ( ldouble ) then
              write ( unit=unit_vtk, fmt='(9es25.16)' ) transpose(tvalue)
            else
              write ( unit=unit_vtk, fmt='(9es16.8)' ) &
                                       real ( transpose(tvalue), kind=sp )
            end if
          end if

        end do

      end do

    else if ( present(vector) ) then

!     vector in nodes

      firstskip = .true.

      do nodenr = mesh%nnodes_blend(lblend+1)+1, mesh%nnodes_blend(lblend+2)

        dof = problem%vec_nodnumdegfd(loffset+nodenr+1,vector%vec) &
                 - problem%vec_nodnumdegfd(loffset+nodenr,vector%vec)

        if ( any(deg > dof) ) then

          if ( firstskip .and. WARN_ON_NO_DATA_IN_NODES ) then
            write(*,'(/3(a/),2a/)') &
              'Warning in write_tensor_vtk:', &
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
          if ( ldouble ) then
            write ( unit=unit_vtk ) transpose(tvalue)
          else
            write ( unit=unit_vtk ) real ( transpose(tvalue), kind=sp )
          end if
        else
          if ( ldouble ) then
            write ( unit=unit_vtk, fmt='(9es25.16)' ) transpose(tvalue)
          else
            write ( unit=unit_vtk, fmt='(9es16.8)' ) &
                                     real ( transpose(tvalue), kind=sp )
          end if
        end if

      end do

    end if

    deallocate ( lgroups )

    close ( unit=unit_vtk )

  end subroutine write_tensor_vtk


! VTK header

  subroutine vtk_header

    character (len=21) :: line
    integer :: values(8)

!   generate date and time

    call date_and_time ( values=values )

    write ( line, '(2(i2.2,a),i4,a,3(i2.2,a))' ) &
      values(3), '-', values(2), '-', values(1), ', ', &
      values(5), ':', values(6), ':', values(7)

!   general header

    write ( unit=unit_vtk, fmt='(a)' ) '# vtk DataFile Version 3.0'
    write ( unit=unit_vtk, fmt='(2a)' ) 'TFEM output generated at ', line
    write ( unit=unit_vtk, fmt='(a)' ) 'ASCII'

  end subroutine vtk_header


! VTK header (binary format)

  subroutine vtk_header_binary

    character (len=21) :: line
    integer :: values(8)

!   generate date and time

    call date_and_time ( values=values )

    write ( line, '(2(i2.2,a),i4,a,3(i2.2,a))' ) &
      values(3), '-', values(2), '-', values(1), ', ', &
      values(5), ':', values(6), ':', values(7)

!   general header

    call write_line ( unit=unit_vtk, line='# vtk DataFile Version 3.0' )
    call write_line ( unit=unit_vtk, line='TFEM output generated at '//line )
    call write_line ( unit=unit_vtk, line='BINARY' )

  end subroutine vtk_header_binary


! VTK geometry, including points and topology for unstructured grid

  subroutine vtk_geometry_ug ( mesh, groups, double, blend )

    type(mesh_t), intent(in) :: mesh

!   element groups to be to be included
!   NOTE: all the nodal points of the mesh are written and only the cells
!   are affected by choosing a subset (of all groups) to be written.
    integer, dimension(:), intent(in) :: groups

!   write data in the "double" (.true.) or "float" (.false.) format.
    logical, intent(in) :: double

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
    integer, intent(in) :: blend

    integer :: elem, elgrp, numdata, i, j, subelem, numelem, grp
    type(element_t), dimension(mesh%nelgrp) :: element

!   choose element

    if ( blend == 0 ) then
      element = mesh%element
    else
      element = mesh%element_blend(:,blend)
    end if

!   test

    if ( any(element(groups)%vtk_elshape == 0) ) then
      write (*,'(3(a/),2a)') 'Error vtk_geometry_ug: ', &
        ' there are element groups that have element shapes that are', &
        ' not (yet) implemented for VTK output.', &
        ' vtk_elshape(:) ='
      write (*,*) element(groups)%vtk_elshape
      stop
    end if

    write ( unit=unit_vtk, fmt='(a)' ) 'DATASET UNSTRUCTURED_GRID'

!   points

    if ( double ) then
      write ( unit=unit_vtk, fmt='(a,i0,a)' ) 'POINTS ', &
        mesh%nnodes_blend(blend+2)-mesh%nnodes_blend(blend+1), ' double'
      write ( unit=unit_vtk, fmt='(3es25.16)' ) &
        ( mesh%coor(i,:), ( 0._dp, j=mesh%ndim+1,3 ), &
          i=mesh%nnodes_blend(blend+1)+1,mesh%nnodes_blend(blend+2) )
    else
      write ( unit=unit_vtk, fmt='(a,i0,a)' ) 'POINTS ', &
        mesh%nnodes_blend(blend+2)-mesh%nnodes_blend(blend+1), ' float'
      write ( unit=unit_vtk, fmt='(3es16.8)' ) &
        ( mesh%coor(i,:), ( 0._dp, j=mesh%ndim+1,3 ), &
          i=mesh%nnodes_blend(blend+1)+1,mesh%nnodes_blend(blend+2) )
    end if

!   topology of the cells

    write ( unit=unit_vtk, fmt=* )
    numdata = sum ( mesh%grpnumel(groups) * element(groups)%vtk_numelm * &
                    (1+element(groups)%vtk_numnod) )
    numelem = sum ( mesh%grpnumel(groups) * element(groups)%vtk_numelm )
    write ( unit=unit_vtk, fmt='(a,i0,1x,i0)' ) 'CELLS ', numelem, numdata
    do grp = 1, size(groups)
      elgrp = groups(grp)
      do elem = 1, mesh%grpnumel(elgrp)
        do subelem = 1, element(elgrp)%vtk_numelm
          write ( unit=unit_vtk, fmt='(11(1x,i0)/3x,10(1x,i0))' ) &
            element(elgrp)%vtk_numnod, &
            mesh%topology(elgrp)%a(element(elgrp)%vtk_index(:,subelem) &
                                       +mesh%numnodtop(elgrp,blend+1),elem) &
                                       -1-mesh%nnodes_blend(blend+1)
        end do
      end do
    end do

!   type of the cells

    write ( unit=unit_vtk, fmt=* )
    write ( unit=unit_vtk, fmt='(a,i0)' ) 'CELL_TYPES ', numelem
    do grp = 1, size(groups)
      elgrp = groups(grp)
      write ( unit=unit_vtk, fmt='(20(1x,i0))' ) &
        ( element(elgrp)%vtk_elshape, &
             elem = 1, mesh%grpnumel(elgrp) * element(elgrp)%vtk_numelm )
    end do

  end subroutine vtk_geometry_ug


! VTK geometry, including points and topology for unstructured grid (binary)

  subroutine vtk_geometry_ug_binary ( mesh, groups, double, blend )

    type(mesh_t), intent(in) :: mesh

!   element groups to be to be included
!   NOTE: all the nodal points of the mesh are written and only the cells
!   are affected by choosing a subset (of all groups) to be written.
    integer, dimension(:), intent(in) :: groups

!   write data in the "double" (.true.) or "float" (.false.) format.
    logical, intent(in) :: double

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
    integer, intent(in) :: blend

    character (len=80) :: line
    integer :: elem, elgrp, numdata, i, j, subelem, numelem, grp
    type(element_t), dimension(mesh%nelgrp) :: element

!   choose element

    if ( blend == 0 ) then
      element = mesh%element
    else
      element = mesh%element_blend(:,blend)
    end if

!   test

    if ( any(element(groups)%vtk_elshape == 0) ) then
      write (*,'(3(a/),2a)') 'Error vtk_geometry_ug_binary: ', &
        ' there are element groups that have element shapes that are', &
        ' not (yet) implemented for VTK output.', &
        ' vtk_elshape(:) ='
      write (*,*) element(groups)%vtk_elshape
      stop
    end if

    call write_line ( unit=unit_vtk, line='DATASET UNSTRUCTURED_GRID' )

!   points

    if ( double ) then
      write ( unit=line, fmt='(a,i0,a)' ) 'POINTS ', &
        mesh%nnodes_blend(blend+2)-mesh%nnodes_blend(blend+1), ' double'
      call write_line ( unit=unit_vtk, line=trim(line) )
      write ( unit=unit_vtk ) &
        ( mesh%coor(i,:), ( 0._dp, j=mesh%ndim+1,3 ), &
          i=mesh%nnodes_blend(blend+1)+1,mesh%nnodes_blend(blend+2) )
    else
      write ( unit=line, fmt='(a,i0,a)' ) 'POINTS ', &
        mesh%nnodes_blend(blend+2)-mesh%nnodes_blend(blend+1), ' float'
      call write_line ( unit=unit_vtk, line=trim(line) )
      write ( unit=unit_vtk ) ( real(mesh%coor(i,:),real32), &
                         ( 0._real32, j=mesh%ndim+1,3 ), &
                   i=mesh%nnodes_blend(blend+1)+1,mesh%nnodes_blend(blend+2) )

    end if

!   topology of the cells

    call write_line ( unit=unit_vtk, line='' )
    numdata = sum ( mesh%grpnumel(groups) * element(groups)%vtk_numelm * &
                    (1+element(groups)%vtk_numnod) )
    numelem = sum ( mesh%grpnumel(groups) * element(groups)%vtk_numelm )
    write ( unit=line, fmt='(a,i0,1x,i0)' ) 'CELLS ', numelem, numdata
    call write_line ( unit=unit_vtk, line=trim(line) )
    do grp = 1, size(groups)
      elgrp = groups(grp)
      do elem = 1, mesh%grpnumel(elgrp)
        do subelem = 1, element(elgrp)%vtk_numelm
          write ( unit=unit_vtk ) &
            int(element(elgrp)%vtk_numnod,int32), &
            int(mesh%topology(elgrp)%a(element(elgrp)%vtk_index(:,subelem) &
                                        +mesh%numnodtop(elgrp,blend+1),elem) &
                                 -1-mesh%nnodes_blend(blend+1),int32)
        end do
      end do
    end do

!   type of the cells

    call write_line ( unit=unit_vtk, line='' )
    write ( unit=line, fmt='(a,i0)' ) 'CELL_TYPES ', numelem
    call write_line ( unit=unit_vtk, line=trim(line) )
    do grp = 1, size(groups)
      elgrp = groups(grp)
      write ( unit=unit_vtk ) &
        ( int(element(elgrp)%vtk_elshape,int32), &
             elem = 1, mesh%grpnumel(elgrp) * element(elgrp)%vtk_numelm )
    end do

  end subroutine vtk_geometry_ug_binary


! first line of the point of cell data

  subroutine vtk_data_format ( mesh, elementwise, groups, blend )

    type(mesh_t), intent(in) :: mesh
    logical, intent(in) :: elementwise
    integer, dimension(:), intent(in) :: groups
    integer, intent(in) :: blend

    integer :: nelem, elgrp

    if ( elementwise ) then

!     set nelem first

      nelem = 0
      do elgrp = 1, mesh%nelgrp
        if ( .not. any ( groups == elgrp ) ) cycle
        nelem = nelem + mesh%grpnumel(elgrp)
      end do

      write ( unit=unit_vtk, fmt=* )
      write ( unit=unit_vtk, fmt='(a,i0)' ) 'CELL_DATA ', nelem

    else

      write ( unit=unit_vtk, fmt=* )
      write ( unit=unit_vtk, fmt='(a,i0)' ) 'POINT_DATA ', &
                       mesh%nnodes_blend(blend+2)-mesh%nnodes_blend(blend+1)
    end if

  end subroutine vtk_data_format


! first line of the point of cell data

  subroutine vtk_data_format_binary ( mesh, elementwise, groups, blend )

    type(mesh_t), intent(in) :: mesh
    logical, intent(in) :: elementwise
    integer, dimension(:), intent(in) :: groups
    integer, intent(in) :: blend

    character (len=21) :: line
    integer :: nelem, elgrp

    if ( elementwise ) then

!     set nelem first

      nelem = 0
      do elgrp = 1, mesh%nelgrp
        if ( .not. any ( groups == elgrp ) ) cycle
        nelem = nelem + mesh%grpnumel(elgrp)
      end do

      call write_line ( unit=unit_vtk, line='' )
      write ( unit=line, fmt='(a,i0)' ) 'CELL_DATA ', nelem
      call write_line ( unit=unit_vtk, line=trim(line) )

    else

      call write_line ( unit=unit_vtk, line='' )
      write ( unit=line, fmt='(a,i0)' ) 'POINT_DATA ', &
                         mesh%nnodes_blend(blend+2)-mesh%nnodes_blend(blend+1)
      call write_line ( unit=unit_vtk, line=trim(line) )

    end if

  end subroutine vtk_data_format_binary


! write geometries of a mesh to a legacy VTK file. Optionally, the normals to
! a line/surface can be written as well
! NOTE: the normals are approximated in a crude way, these routines should
! therefore only be used for postprocessing / visualization

  subroutine write_geometry_vtk ( mesh, filename, point, curve, surface, &
    write_normals, double, binary )

    use meshgen_m
    use problem_m
    use limits_m

!   the mesh of which to write the geometry (as indicated by point, curve or
!   surface)
    type(mesh_t), intent(in) :: mesh

!   the filename for writing the mesh to (it must have the .vtk extension).
    character (len=*), intent(in) :: filename

!   the point, curve or surface number to write to a vtk-file
!   NOTE: at least one (and not more) of these arguments should be supplied
    integer, intent(in), optional :: point, curve, surface

!   optional: plot the normals to a curve (2D) or surface (3D)
    logical, intent(in), optional :: write_normals

!   write data in the "double" format (.true.) or the "float" format (.false.).
!   NOTE: files are 40-50 % larger for the double format.
!   default=WRITE_DOUBLE_VTK (set in limits_m).
    logical, intent(in), optional :: double

!   write data in the binary format (.true.) or ASCII format (.false.)
!   default=WRITE_BINARY_VTK (set in limits_m).
!   The advantages of a binary file are:
!   a. the data is identical to the internal representation.
!   b. the files are significantly smaller than ASCII vtk-files.
!   c. reading the file might be somewhat faster.
!   A disadvantage is that the file cannot be interpreted by humans.
    logical, intent(in), optional :: binary


!   temporary mesh for plotting the geometry
    type(mesh_t) :: mesh_plot

!   problem for normals
    type(input_probdef_t) :: input_probdef
    type(problem_t), target :: problem
    type(vector_t) :: normals
    type(subscriptvec_t) :: nx, ny, nz

    character (len=21) :: line
    integer :: i, values(8)
    real(dp), allocatable, dimension(:) :: normal
    real(dp), allocatable, dimension(:,:) :: bd, x
    logical :: lwrite_normals, warn_old, ldouble, lbinary

    call check_blend ( mesh, 'write_geometry_vtk' )

    lwrite_normals = set_optional ( variable=write_normals, default=.false. )
    ldouble = set_optional ( variable=double, default=WRITE_DOUBLE_VTK )
    lbinary = set_optional ( variable=binary, default=WRITE_BINARY_VTK )

!   first perform some checks

    if ( count( [present(point), present(curve), present(surface)] ) /= 1 ) then
      write(*,'(/2(a/))') &
        'Error write_geometry_vtk: supply either a point OR a curve OR', &
        'a surface'
      stop
    end if

    if ( present(point) ) then

!     generate date and time

      call date_and_time ( values=values )

      write ( line, '(2(i2.2,a),i4,a,3(i2.2,a))' ) &
        values(3), '-', values(2), '-', values(1), ', ', &
        values(5), ':', values(6), ':', values(7)

      if ( lbinary ) then

!       general header
        open ( unit=unit_vtk, file=filename, access='stream', &
               form='unformatted', status='replace', convert='big_endian' )

        call write_line ( unit=unit_vtk, line='# vtk DataFile Version 3.0' )
        call write_line ( unit=unit_vtk, &
                                       line='TFEM output generated at '//line )
        call write_line ( unit=unit_vtk, line='BINARY' )
        call write_line ( unit=unit_vtk, line='DATASET UNSTRUCTURED_GRID' )
        if ( ldouble ) then
          write ( unit=line, fmt='(a)' ) 'POINTS 1 double'
          call write_line ( unit=unit_vtk, line=trim(line) )
          write ( unit=unit_vtk ) &
            [ mesh%coor(mesh%points(point),:), ( 0._dp, i=mesh%ndim+1,3 ) ]
        else
          write ( unit=line, fmt='(a)' ) 'POINTS 1 float'
          call write_line ( unit=unit_vtk, line=trim(line) )
          write ( unit=unit_vtk ) &
            [ real(mesh%coor(mesh%points(point),:),real32), &
            ( 0._real32, i=mesh%ndim+1,3 ) ]
        end if

      else

!       general header
        open ( unit=unit_vtk, file=filename, status='replace' )

        write ( unit=unit_vtk, fmt='(a)' ) '# vtk DataFile Version 3.0'
        write ( unit=unit_vtk, fmt='(2a)' ) 'TFEM output generated at ', line
        write ( unit=unit_vtk, fmt='(a)' ) 'ASCII'
        write ( unit=unit_vtk, fmt='(a)' ) 'DATASET UNSTRUCTURED_GRID'
        if ( ldouble ) then
          write ( unit=unit_vtk, fmt='(a)' ) 'POINTS 1 double'
          write ( unit=unit_vtk, fmt='(3es25.16)' ) &
            [ mesh%coor(mesh%points(point),:), ( 0._dp, i=mesh%ndim+1,3 ) ]
        else
          write ( unit=unit_vtk, fmt='(a)' ) 'POINTS 1 float'
          write ( unit=unit_vtk, fmt='(3es16.8)' ) &
            [ real(mesh%coor(mesh%points(point),:),real32), &
            ( 0._real32, i=mesh%ndim+1,3 ) ]
        end if

      end if

      close ( unit=unit_vtk )

      return

    end if

    if ( present(curve) ) then

      call mesh_skeleton ( mesh_plot, mesh%curves(curve)%nnodes, &
        mesh%curves(curve)%nelem, mesh%curves(curve)%element%elshape, &
        mesh%curves(curve)%ndim )

      mesh_plot%coor = mesh%coor(mesh%curves(curve)%nodes,:)
      mesh_plot%topology(1)%a = mesh%curves(curve)%topology(:,:,1)

    end if

    if ( present(surface) ) then

      call mesh_skeleton ( mesh_plot, mesh%surfaces(surface)%nnodes, &
       mesh%surfaces(surface)%nelem, mesh%surfaces(surface)%element%elshape, &
       mesh%surfaces(surface)%ndim )

      mesh_plot%coor = mesh%coor(mesh%surfaces(surface)%nodes,:)
      mesh_plot%topology(1)%a = mesh%surfaces(surface)%topology(:,:,1)

    end if

    if ( lwrite_normals .and. mesh_plot%ndim == 3 .and. &
           mesh_plot%element(1)%globalshape == 'line' ) then
      write(*,'(/a/)') &
        'Warning write_geometry_vtk: skipping normal on line in 3D'
      lwrite_normals = .false.
    end if

!   add blocks manually to avoid 'domain of blocks near zero' in fill_mesh_parts

    allocate(bd(mesh%ndim,2))

    do i=1,mesh%ndim
      bd(i,1)=minval(mesh%coor(:,i))-1.e-10_dp
      bd(i,2)=maxval(mesh%coor(:,i))+1.e-10_dp
    end do

    call add_to_mesh ( mesh_plot, blocks=[(2,i=1,mesh%ndim)], &
      blocksdomain=bd )

!   temporarily disable warning for fill_mesh_parts_sidelem (in limits_m)

    warn_old = WARN_ON_TEST_FOR_MULTIPLE_SIDELEM
    WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = .false.

    call fill_mesh_parts ( mesh_plot )

    deallocate(bd)

    WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = warn_old

    if ( lwrite_normals ) then

!     create a problem to write the normals

      call create_input_probdef ( mesh_plot, input_probdef, nvec=1 )

      input_probdef%elementdof(1)%a = 1
      input_probdef%vec_elementdof(1)%a(:,1) = mesh%ndim

      call problem_definition ( input_probdef, mesh_plot, problem )

      call create ( problem, normals, vec=1 )

      call create_subscript_vector ( mesh_plot, problem, nx, degfd=1, vec=1 )
      call create_subscript_vector ( mesh_plot, problem, ny, degfd=2, vec=1 )

!     approximate the normal in each nodal point
!     NOTE: strictly speaking, the normal at a nodal point is ill-defined. Here
!     we approximate the normal for each element and add them to the nodal
!     nodal points connected to that element

      call derive_vector ( mesh_plot, problem, vector=normals, &
        elemsub=deriv_normals )

      allocate ( normal(mesh_plot%ndim) )

!     normalize the normal vectors to length 1

      if ( mesh_plot%ndim == 2 ) then

        do i = 1, size(nx%s)
          normal = [ normals%u(nx%s(i)), normals%u(ny%s(i)) ]
          normal = normal / sqrt(dot_product(normal,normal))
          normals%u(nx%s(i)) = normal(1)
          normals%u(ny%s(i)) = normal(2)
        end do

      else

        call create_subscript_vector ( mesh_plot, problem, nz, degfd=3, vec=1 )

        do i = 1, size(nx%s)
          normal = [ normals%u(nx%s(i)), normals%u(ny%s(i)), &
                                                            normals%u(nz%s(i)) ]
          normal = normal / sqrt(dot_product(normal,normal))
          normals%u(nx%s(i)) = normal(1)
          normals%u(ny%s(i)) = normal(2)
          normals%u(nz%s(i)) = normal(3)
        end do

      end if

      deallocate ( normal )
      call delete ( input_probdef )

    end if

!   write the vtks

    if ( .not. lwrite_normals ) then
      call write_mesh_vtk ( mesh_plot, filename, double=double, binary=binary )
    else
      call write_vector_vtk ( mesh_plot, problem, filename, 'normals', &
        vector=normals, double=double, binary=binary )
      call delete ( problem )
      call delete ( normals )
    end if

    call delete ( mesh_plot )

  contains

!   approximate the normal to an interface mesh

    subroutine deriv_normals ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors, elemvec, elemwts )

      use math_defs_m

      type(mesh_t), intent(in) :: mesh
      type(problem_t), intent(in) :: problem
      integer, intent(in) :: elgrp, elem
      logical, intent(in) :: first, last
      type(coefficients_t), intent(in) :: coefficients
      type(oldvectors_t), intent(in) :: oldvectors
      real(dp), intent(out), dimension(:) :: elemvec, elemwts

      integer :: numnod, pt(3), ndim

      numnod = mesh%element(1)%numnod
      ndim = mesh%ndim

      if ( first ) allocate ( x(numnod,ndim), normal(ndim) )

!     get coordinates of nodes of interface

      call get_coordinates ( mesh, elgrp, elem, x )

      if ( ndim == 2 ) then

!       approximate the normal on the line element

        select case ( mesh%element(1)%elshape )
          case (1)
            normal = [(x(2,2)-x(1,2)), -(x(2,1)-x(1,1))]
          case (2)
            normal = [(x(3,2)-x(1,2)), -(x(3,1)-x(1,1))]
          case default
            write(*,'(/a/a,i0/)')'deriv_normals: ', &
                ' this elementshape is not available: ', mesh%element(1)%elshape
            stop
        end select

      else

!       get the points to determine the normal on the surface element

        select case ( mesh%element(1)%elshape )
          case (3,10)
            pt = [1,2,3]
          case (4,7)
            pt = [1,3,5]
          case (5,9)
            pt = [1,2,4]
          case (6)
            pt = [1,3,7]
          case default
            write(*,'(/a/a,i0/)')'deriv_normals: ', &
                ' this elementshape is not available: ', mesh%element(1)%elshape
            stop
        end select

!       approximate the normal on the surface element using a cross product

        normal = cross_product ( x(pt(2),:)-x(pt(1),:), x(pt(3),:)-x(pt(1),:) )

      end if

!     add the normal to the element vector

      if ( ndim == 2 ) then
        elemvec(:size(elemvec)/2) = normal(1)
        elemvec(size(elemvec)/2+1:) = normal(2)
      else
        elemvec(:size(elemvec)/3) = normal(1)
        elemvec(size(elemvec)/3+1:2*size(elemvec)/3) = normal(2)
        elemvec(2*size(elemvec)/3+1:) = normal(3)
      end if

      elemwts = 1

      if ( last ) deallocate ( x, normal )

    end subroutine deriv_normals

  end subroutine write_geometry_vtk


! Write ALL geometries of a mesh, each to a separate legacy VTK file.

  subroutine write_geometries_vtk ( mesh, write_normals, points, curves, &
    surfaces, double, binary )

!   the mesh of which to write the geometries
    type(mesh_t), intent(in) :: mesh

!   optional: plot the normals to a curve (2D) or surface (3D)
    logical, intent(in), optional :: write_normals

!   logical whether to write points, curves or surfaces to a vtk-file
    logical, intent(in), optional :: points, curves, surfaces

!   write data in the "double" format (.true.) or the "float" format (.false.).
!   NOTE: files are 40-50 % larger for the double format.
!   default=WRITE_DOUBLE_VTK (set in limits_m).
    logical, intent(in), optional :: double

!   write data in the binary format (.true.) or ASCII format (.false.)
!   default=WRITE_BINARY_VTK (set in limits_m).
!   The advantages of a binary file are:
!   a. the data is identical to the internal representation.
!   b. the files are significantly smaller than ASCII vtk-files.
!   c. reading the file might be somewhat faster.
!   A disadvantage is that the file cannot be interpreted by humans.
    logical, intent(in), optional :: binary

    logical :: lpoints, lcurves, lsurfaces
    integer :: i
    character(len=30) :: filename


    call check_blend ( mesh, 'write_geometries_vtk' )

    lpoints = set_optional ( variable=points, default=.true. )
    lcurves = set_optional ( variable=curves, default=.true. )
    lsurfaces = set_optional ( variable=surfaces, default=.true. )


    if ( lpoints ) then
      do i = 1, mesh%npoints
        write(filename,'(a,i4.4,a)') 'point_',i,'.vtk'
        call write_geometry_vtk ( mesh, filename, point=i, double=double, &
          binary=binary )
      end do
    end if

    if ( lcurves ) then
      do i = 1, mesh%ncurves
        write(filename,'(a,i4.4,a)') 'curve_',i,'.vtk'
        if ( mesh%ndim == 2 ) then
          call write_geometry_vtk ( mesh, filename, curve=i, &
            write_normals=write_normals, double=double, binary=binary )
        else
          call write_geometry_vtk ( mesh, filename, curve=i, double=double, &
            binary=binary )
        end if
      end do
    end if

    if ( lsurfaces ) then
      do i = 1, mesh%nsurfaces
        write(filename,'(a,i4.4,a)') 'surface_',i,'.vtk'
        call write_geometry_vtk ( mesh, filename, surface=i, &
          write_normals=write_normals, double=double, binary=binary )
      end do
    end if

  end subroutine write_geometries_vtk

end module vtk_utils_m
