
! Copyright (C) 2006-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines routines for output to the Tecplot format (extension .plt).

module tec_utils_m

  use kind_defs_m
  use mesh_m
  use problem_defs_m
  use system_vector_m
  use vector_m
  use limits_m, only: WARN_ON_NO_DATA_IN_NODES

  implicit none

  integer :: unit_tec = 13 ! unit number for writing

  ! interface to support older interfaces

  interface write_objects_tec
    module procedure write_objects_tecplot
  end interface write_objects_tec

  interface write_scalar_tec
    module procedure write_scalar_tecplot
  end interface write_scalar_tec


contains


! write objects to a Tecplot file

  subroutine write_objects_tecplot ( mesh, dataname, filename, object1, &
    object2, append, topology, time )

    type(mesh_t), intent(in) :: mesh

!   the filename for writing the data to (it must have the .plt extension).
    character (len=*), intent(in) :: filename

!   the dataname for the object
    character (len=*), intent(in) :: dataname

!   if these are present only objects object1,...,object2 are plotted.
!   If only object1 is present one object is plotted.
    integer, intent(in), optional :: object1, object2

    logical, intent(in), optional :: append

    logical, intent(in), optional :: topology

    real(dp), intent(in), optional :: time

    real(dp), dimension(mesh%ndim) :: x

    integer :: object, node, obj1, obj2

    integer :: num_object, numelem

    logical :: fileexists, fileappend, firstopen, obj_topology


!   set fileappend

    if ( present(append) ) then
      fileappend = append
    else
      fileappend = .false.
    end if

!   object topology

    if ( present(topology) ) then
      obj_topology = topology          ! write object topology as well
    else
      obj_topology = .false.           ! write only object coordinates
    end if

!   open file

    inquire ( file=filename, exist=fileexists )

    if ( fileexists .and. fileappend ) then

!     add scalar point data to existing file

      firstopen = .false.

      open ( unit=unit_tec, file=filename, status='old', position='append' )

    else

      firstopen = .true.

      open ( unit=unit_tec, file=filename, status='replace' )

    end if

!   check mesh

    call check_mesh ( mesh, 'write_objects_tecplot' )
    call check_blend ( mesh, 'write_objects_tecplot' )

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

    if ( obj1 < 1 .or. obj1 > mesh%nobjects .or.  &
         obj2 < 1 .or. obj2 > mesh%nobjects ) then

      write(*,'(/2a/2(a,i0),/a,i0/)') &
        'Error: object1 and/or object2 in the heading of ', &
        'plot_objects is ', &
        'out of range. object1 is ', obj1, '; object2 is ', obj2, &
        'whereas the number of objects is ', mesh%nobjects
      stop

    end if

!   plot coordinates of object

    num_object = sum (mesh%objects(obj1:obj2)%nnodes)

    write( unit=unit_tec, fmt='(2a)' ) '# ', dataname

    if (mesh%ndim == 2) then
      write( unit=unit_tec, fmt='(a)' ) 'variables = "x", "y"'
    else if (mesh%ndim == 3) then
      write( unit=unit_tec, fmt='(a)' ) 'variables = "x", "y", "z"'
    end if

    if ( obj_topology ) then

      if ( any ( .not. mesh%objects(obj1:obj2)%topol ) ) then
        write(*,'(/a/)') &
          'Error in write_objects_tecplot: object(s) have no topology '
        stop
      end if

      if ( any ( mesh%objects(obj1:obj2)%element%tec_elshape /= &
                         mesh%objects(obj1)%element%tec_elshape ) ) then
        write(*,'(/2a/)') &
          'Error in write_objects_tecplot:', &
                       ' all objects must have the same element shape '
        stop
      end if

      numelem = sum (mesh%objects(obj1:obj2)%nelem * &
                                mesh%objects(obj1:obj2)%element%tec_numelm )
      write ( unit=unit_tec, fmt='(2(a,i0),2a)' )  &
            'zone n=', num_object, ',e=', numelem, &
            ',datapacking=point, zonetype=', &
                                      mesh%objects(obj1)%element%tec_elshape

      if (.not. firstopen) then
        write ( unit=unit_tec, fmt='(a)' )  ',connectivitysharezone=1'
      end if

    else

      write( unit=unit_tec, fmt='(a,i0,a)' ) &
                        'zone i=', num_object, ',datapacking=point'
    end if

    if ( present(time) ) then
      write ( unit=unit_tec, fmt='(a,e14.4)' ) &
                             ',strandid=15, solutiontime=', time
    end if

    do object = obj1, obj2
      do node = 1, mesh%objects(object)%nnodes
        x = mesh%objects(object)%coor(node,1:mesh%ndim)
        write ( unit=unit_tec, fmt=* ) x !, object
      end do
    end do

!   plot mesh topology for objects

    if ( obj_topology .and. firstopen ) then

      call tec_geometry_object ( mesh, obj1, obj2 )

    end if

    close (unit_tec)

  end subroutine write_objects_tecplot


! write multiple data

  subroutine write_mdata_tecplot ( mesh, problem, filename, dataname, &
    oldvectors, append, remesh, time )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the data to (it must have the .plt extension).
    character (len=*), intent(in) :: filename

!   the dataname for the scalar, for example 'pressure'.
!   NOTE: the size of the array determines the number of scalars written
    character (len=*), dimension(:), intent(in) :: dataname

!   contains the multiple data vectors
!   NOTE: the size of the array dataname determines the number of scalars
!   written. Each vector must be a scalar (single degree in each node).
    type(oldvectors_t), intent(in) :: oldvectors

!   if .true.: append data to existing file in order to avoid writing
!   coordinate and topology data more than ones.
!   default=.false.
    logical, intent(in), optional :: append

!   if .true.: write mesh coord and topology every time for each zone
!   default=.false.
    logical, intent(in), optional :: remesh

!   time for transient data

    real(dp), intent(in), optional :: time

    integer :: i, n

    logical :: fileexists, fileappend, firstopen, fileremesh

    call check_blend ( mesh, 'write_mdata_tecplot' )

!   set fileappend

    if ( present(append) ) then
      fileappend = append
    else
      fileappend = .false.
    end if

!   set fileremesh

    if ( present(remesh) ) then
      fileremesh = remesh
    else
      fileremesh = .false.
    end if

!   open file

    inquire ( file=filename, exist=fileexists )

    if ( fileexists .and. fileappend ) then

!     add scalar point data to existing file

      firstopen = .false.
      if ( fileremesh )  firstopen = .true.

      open ( unit=unit_tec, file=filename, status='old', position='append' )

      call tec_geometry_point ( mesh, firstopen, dataname=dataname, time=time )

    else

      firstopen = .true.

      open ( unit=unit_tec, file=filename, status='replace' )

      call tec_header
      call tec_geometry_point ( mesh, firstopen, dataname=dataname, time=time )

    end if

!   write first part

    do i = 1, size(dataname)

      write( unit=unit_tec, fmt=* )
      write( unit=unit_tec, fmt='(2a)' ) '# ', dataname(i)

      do n = 1, mesh%nnodes
        write(unit=unit_tec, fmt='(es16.8)') oldvectors%v(i)%p%u(n)
      end do

    end do

!   element connectivity for tecplot

    if ( firstopen ) call tec_geometry_ug ( mesh )

    if ( fileremesh ) write( unit=unit_tec, fmt=* )

    close (unit_tec)

  end subroutine write_mdata_tecplot


! write a scalar field to a Tecplot file

  subroutine write_scalar_tecplot ( mesh, problem, filename, dataname, physq, &
    degfd, sysvector, vector, append, time )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the data to (it must have the .plt extension).
    character (len=*), intent(in) :: filename

!   the dataname for the scalar, for example 'pressure'.
    character (len=*), intent(in) :: dataname

!   the physical quantity to be plotted
    integer, intent(in), optional :: physq

!   the degree of freedom to be plotted (within the physical quantity when
!   physq is also present)
    integer, intent(in), optional :: degfd

!   sysvector to be plotted
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be plotted
    type(vector_t), intent(in), optional :: vector

!   when .true.: append data to existing file in order to avoid writing
!   coordinate and topology data more than ones.
!   default=.false.
    logical, intent(in), optional :: append

!   time for transient data

    real(dp), intent(in), optional :: time

    logical :: fileexists, fileappend, firstskip, firstopen
    integer :: nodenr, pos(problem%maxnoddegfd), dof, deg

    real(dp) :: svalue


    call check_mesh ( mesh, 'write_scalar_tecplot' )
    call check_blend ( mesh, 'write_scalar_tecplot' )
    call check_problem ( problem, 'write_scalar_tecplot', mesh )

    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in write_scalar_tecplot: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(sysvector) ) then
      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in write_scalar_tecplot: sysvector has not been created '
        stop
      end if
      if ( problem%probnr /= sysvector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_scalar_tecplot: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in sysvector = ',  sysvector%probnr
        stop
      end if
    else if ( present(vector) ) then
      if ( .not. vector%created ) then
        write(*,'(/a/)') &
          'Error in write_scalar_tecplot: vector has not been created '
        stop
      end if
      if ( vector%elementwise ) then
        write(*,'(/a/a/)') &
          'Error in write_scalar_tecplot:', &
          ' vector stored elementwise not supported'
        stop
      end if
      if ( problem%probnr /= vector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in write_scalar_tecplot: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in vector = ', vector%probnr
        stop
      end if
    end if

!   set deg

    if ( present(degfd) ) then
      if ( degfd < 1 ) then
        write(*,'(/a/)') &
          'Error in write_scalar_tecplot: degfd must be larger than zero'
        stop
      end if
      deg = degfd
    else
      deg = 1 ! take first degree
    end if

!   set fileappend

    if ( present(append) ) then
      fileappend = append
    else
      fileappend = .false.
    end if

!   open file

    inquire ( file=filename, exist=fileexists )

    if ( fileexists .and. fileappend ) then

!     add scalar point data to existing file

      firstopen = .false.

      open ( unit=unit_tec, file=filename, status='old', position='append' )

      call tec_geometry_point ( mesh, firstopen, time=time )

    else

      firstopen = .true.

      open ( unit=unit_tec, file=filename, status='replace' )

      call tec_header
      call tec_geometry_point ( mesh, firstopen, time=time )

    end if

!   write first part

    write( unit=unit_tec, fmt=* )
    write( unit=unit_tec, fmt='(2a)' ) '# ', dataname

!   write scalar

    if ( present(sysvector) ) then

!     sysvector in nodes

      firstskip = .true.

      do nodenr = 1, mesh%nnodes

        if ( present(physq) ) then
          call pos_array_node ( problem, nodenr, dof, pos, physqarr=[physq] )
        else
          call pos_array_node ( problem, nodenr, dof, pos )
        end if

        if ( deg > dof ) then

          if ( firstskip .and. WARN_ON_NO_DATA_IN_NODES ) then
            write(*,'(/3(a/),2a/)') &
              'Warning in write_scalar_tecplot:', &
              ' there are nodes where the required unknown in the ', &
              ' sysvector is undefined. The scalar value is set to zero. ', &
              ' dataname = ', dataname
          end if
          firstskip = .false.
          svalue = 0

        else

          svalue = sysvector%u( pos(deg) )

        end if

        write ( unit=unit_tec, fmt='(es16.8)' ) svalue

      end do

    else if ( present(vector) ) then

!     vector in nodes

      firstskip = .true.

      do nodenr = 1, mesh%nnodes

        dof = problem%vec_nodnumdegfd(nodenr+1,vector%vec) &
                 - problem%vec_nodnumdegfd(nodenr,vector%vec)

        if ( deg > dof ) then

          if ( firstskip .and. WARN_ON_NO_DATA_IN_NODES ) then
            write(*,'(/3(a/),2a/)') &
              'Warning in write_scalar_tecplot:', &
              ' there are nodes where the required unknown in the ', &
              ' vector is undefined. The scalar value is set to zero. ', &
              ' dataname = ', dataname
          end if
          firstskip = .false.
          svalue = 0

        else

          svalue = vector%u( problem%vec_nodnumdegfd(nodenr,vector%vec) + deg )

        end if

        write ( unit=unit_tec, fmt='(es16.8)' ) svalue

      end do

    end if

!   element connectivity for tecplot

    if ( firstopen ) call tec_geometry_ug ( mesh )

  end subroutine write_scalar_tecplot


! Tecplot header

  subroutine tec_header

    character (len=21) :: line
    integer :: values(8)

!   generate date and time

    call date_and_time ( values=values )

    write ( line, '(2(i2.2,a),i4,a,3(i2.2,a))' ) &
      values(3), '-', values(2), '-', values(1), ', ', &
      values(5), ':', values(6), ':', values(7)

!   general header

    write ( unit=unit_tec, fmt='(a)' ) '# Tecplot DataFile Version 1.0'
    write ( unit=unit_tec, fmt='(2a)' ) '# TFEM output generated at ', line

  end subroutine tec_header


! Tecplot geometry, including points and topology for unstructured grid

  subroutine tec_geometry_point ( mesh, firstopen, dataname, time )

    type(mesh_t), intent(in) :: mesh

    logical, intent(in) :: firstopen

    character(len=*), dimension(:), intent(in), optional :: dataname

    real(dp), intent(in), optional :: time

    integer :: i, j, numelem

!   test

    if ( any(mesh%element(:)%tec_numnod == 0) ) then
      write (*,'(3(a/),2a)') 'Error tec_geometry_ug: ', &
        ' there are element groups that have element shapes that are', &
        ' not (yet) implemented for Tecplot output.'
      stop
    end if

    numelem = sum ( mesh%grpnumel(:) * mesh%element(:)%tec_numelm )
!
    if ( firstopen ) then

      if ( present(dataname) ) then

        if ( mesh%ndim == 2 ) then
          write ( unit=unit_tec, fmt='(a)' ) 'variables = "x", "y",'
        else if ( mesh%ndim == 3 ) then
          write ( unit=unit_tec, fmt='(a)' ) 'variables = "x", "y", "z",'
        end if

        do i = 1, size(dataname)

          if (i /= size(dataname)) then
            write( unit = unit_tec, fmt='(3a)' ) '"', trim(dataname(i)), '",'
          else
            write( unit = unit_tec, fmt='(3a)' ) '"', trim(dataname(i)), '"'
          end if

        end do

      else

        if ( mesh%ndim == 2 ) then
          write ( unit=unit_tec, fmt='(a)' ) 'variables = "x", "y", "f"'
        else if ( mesh%ndim == 3 ) then
          write ( unit=unit_tec, fmt='(a)' ) 'variables = "x", "y", "z", "f"'
        end if

      end if

      write ( unit=unit_tec, fmt='(2(a,i0),2a)' ) &
            'zone n=', mesh%nnodes, ',e=', numelem, &
            ',datapacking=block, zonetype=', mesh%element(1)%tec_elshape

      if ( present(time) ) then
        write ( unit=unit_tec, fmt='(a,e14.4)' ) &
                               ',strandid=5, solutiontime=', time
      end if

!     points

      write ( unit=unit_tec, fmt=* )
      do j = 1, mesh%ndim
        write ( unit=unit_tec, fmt='(3es16.8)' ) &
                              ( mesh%coor(i,j), i=1,mesh%nnodes )
      end do

    else

      write( unit=unit_tec, fmt=* )

      if (mesh%ndim == 2) then
        write( unit=unit_tec, fmt='(2(a,i0),4a)' ) &
          'zone n=', mesh%nnodes, ',e=', numelem, &
          ',datapacking=block, zonetype=', mesh%element(1)%tec_elshape,  &
          ',varsharelist=([1,2]=1)', &
          ',connectivitysharezone=1'
      else if (mesh%ndim == 3) then
        write( unit=unit_tec, fmt='(2(a,i0),4a)' ) &
          'zone n=', mesh%nnodes, ',e=', numelem, &
          ',datapacking=block, zonetype=', mesh%element(1)%tec_elshape,  &
          ',varsharelist=([1,2,3]=1)', &
          ',connectivitysharezone=1'
      end if

      if ( present(time) ) then
        write ( unit=unit_tec, fmt='(a,e14.4)' ) &
                             ',strandid=5, solutiontime=', time
      end if

    end if

  end subroutine tec_geometry_point


! Tecplot geometry: topology for unstructured grid

  subroutine tec_geometry_ug ( mesh )

    type(mesh_t), intent(in) :: mesh

    integer :: elem, elgrp, subelem

!   topology of the cells

    write ( unit=unit_tec, fmt=* )

    do elgrp = 1, mesh%nelgrp
      write( unit=unit_tec, fmt='(2a)' ) '# element connectivity for ', &
        mesh%element(elgrp)%tec_elshape
      do elem = 1, mesh%grpnumel(elgrp)
        do subelem = 1, mesh%element(elgrp)%tec_numelm
          write ( unit=unit_tec, fmt='(11(1x,i0)/3x,10(1x,i0))' ) &
            mesh%topology(elgrp)%a(&
                &mesh%element(elgrp)%tec_index(:,subelem),elem)
        end do
      end do
    end do

  end subroutine tec_geometry_ug


! Tecplot geometry: topology for objects

  subroutine tec_geometry_object ( mesh, obj1, obj2 )

    type(mesh_t), intent(in) :: mesh

    integer, intent(in) :: obj1, obj2

    integer :: obj, elem, subelem, node_inc

!   topology of the object mesh

    write ( unit=unit_tec, fmt=* )

    node_inc = 0

    do obj = obj1, obj2
      write( unit=unit_tec, fmt='(3a,i6)' ) '# element connectivity for ', &
        mesh%objects(obj)%element%tec_elshape, 'nelem=', mesh%objects(obj)%nelem
      do elem = 1, mesh%objects(obj)%nelem
        do subelem = 1, mesh%objects(obj)%element%tec_numelm
          write ( unit=unit_tec, fmt='(11(1x,i0)/3x,10(1x,i0))' ) &
            mesh%objects(obj)%topology(&
                &mesh%objects(obj)%element%tec_index(:,subelem),elem) + node_inc
        end do
      end do

      node_inc = node_inc + mesh%objects(obj)%nnodes

    end do

    write ( unit=unit_tec,fmt=*) ''

  end subroutine tec_geometry_object

end module tec_utils_m

