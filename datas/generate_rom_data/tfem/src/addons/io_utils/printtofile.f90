
! Copyright (C) 2006-2013 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Print data to a file

module printtofile_m

  use kind_defs_m
  use system_defs_m
  use vector_defs_m
  use mesh_m
  use problem_defs_m

  implicit none


! print blank line after one line of data in grid printing.
! Needed for gnuplot. Set to .false. for matlab and reshape data after reading.
  logical :: blank_line = .true.

! include node number for data line
  logical :: include_node_number = .false.

  integer :: unit_printtofile = 13 ! unit number for writing

  real(dp), private :: boxmax(3), boxmin(3)


contains

! print a sysvector or vector

  subroutine printtofile ( mesh, problem, filename, curve, surface, nnx, nny, &
    nnz, xmax, xmin, ymax, ymin, zmax, zmin, vector, sysvector )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the data to
    character (len=*), intent(in) :: filename

!   curve number for plotting
    integer, intent(in), optional :: curve

!   surface number for plotting
    integer, intent(in), optional :: surface

!   nnx and nny are the number of nodes on the regular grid in x and y
!   direction respectively (if nnz is not present)
!   nnx, nny and nnz are the number of nodes on the regular grid in x, y and z
!   direction respectively
    integer, intent(in), optional :: nnx, nny, nnz

!   maximum and minimum coordinates used for printing
    real(dp), intent(in), optional :: xmin, xmax, ymin, ymax, zmin, zmax

!   vector or sysvector to be printed
!   NOTE: if no vector or sysvector is present, only a nodal point number with
!   coordinates is printed (except for grid printing with nnx,nny,nnz).
    type(vector_t), intent(in), optional :: vector
    type(sysvector_t), intent(in), optional :: sysvector




!   testing

    call check ( mesh, 'printtofile' )
    call check ( problem, 'printtofile', mesh )

    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in printtofile: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(sysvector) ) then
      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in printtofile: sysvector has not been created '
        stop
      end if
      if ( problem%probnr /= sysvector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in printtofile: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in sysvector = ',  sysvector%probnr
        stop
      end if
    else if ( present(vector) ) then
      if ( .not. vector%created ) then
        write(*,'(/a/)') &
          'Error in printtofile: vector has not been created '
        stop
      end if
      if ( problem%probnr /= vector%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in printtofile: ',&
          ' problem number in problem = ', problem%probnr, &
          ' whereas problem number in vector = ', vector%probnr
        stop
      end if
    end if

!   set printing box

    boxmax = 1.0e+100_dp
    boxmin = -1.0e+100_dp

    if ( present(xmax) ) boxmax(1) = xmax
    if ( present(ymax) ) boxmax(2) = ymax
    if ( present(zmax) ) boxmax(3) = zmax
    if ( present(xmin) ) boxmin(1) = xmin
    if ( present(ymin) ) boxmin(2) = ymin
    if ( present(zmin) ) boxmin(3) = zmin


!   how to print?

    if ( present(curve) ) then

!     print on curve

      if ( present(vector) ) then
        call printtofile_geometry ( mesh, problem, mesh%curves(curve), &
          filename, vector )
      else if ( present(sysvector) ) then
        call printtofile_geometry ( mesh, problem, mesh%curves(curve), &
          filename, sysvector=sysvector )
      else
!       coordinates
        call printtofile_geometry ( mesh, problem, mesh%curves(curve), &
          filename )
      end if

    else if ( present(surface) .and. present(nnx) .and. present(nny) ) then

!     print on surface with regular grid

      if ( present(vector) ) then
        call printtofile_2dgrid ( mesh, problem, nnx, nny, filename, vector, &
          surface=surface )
      else if ( present(sysvector) ) then
        call printtofile_2dgrid ( mesh, problem, nnx, nny, filename, &
          sysvector=sysvector, surface=surface )
      else
        write(*,'(/a/a/)') &
          'Error in printtofile: one of vector or sysvector ', &
          ' must be present in the heading for regular grid on surface.'
        stop
      end if

    else if ( present(surface) ) then

!     print on surface

      if ( present(vector) ) then
        call printtofile_geometry ( mesh, problem, mesh%surfaces(surface), &
          filename, vector )
      else if ( present(sysvector) ) then
        call printtofile_geometry ( mesh, problem, mesh%surfaces(surface), &
          filename, sysvector=sysvector )
      else
!       coordinates
        call printtofile_geometry ( mesh, problem, mesh%surfaces(surface), &
          filename )
      end if

    else if ( present(nnx) .and. present(nny) .and. present(nnz) ) then

!     3D grid

      if ( present(vector) ) then
        call printtofile_3dgrid ( mesh, problem, nnx, nny, nnz, filename, &
          vector )
      else if ( present(sysvector) ) then
        call printtofile_3dgrid ( mesh, problem, nnx, nny, nnz, filename, &
          sysvector=sysvector )
      else
        write(*,'(/a/a/)') &
          'Error in printtofile: one of vector or sysvector ', &
          ' must be present in the heading for 3D regular grid.'
        stop
      end if

    else if ( present(nnx) .and. present(nny) ) then

!     2D grid

      if ( present(vector) ) then
        call printtofile_2dgrid ( mesh, problem, nnx, nny, filename, vector )
      else if ( present(sysvector) ) then
        call printtofile_2dgrid ( mesh, problem, nnx, nny, filename, &
          sysvector=sysvector )
      else
        write(*,'(/a/a/)') &
          'Error in printtofile: one of vector or sysvector ', &
          ' must be present in the heading for 2D regular grid.'
        stop
      end if

    else

!     full printing

      if ( present(vector) ) then
        call printtofile_all ( mesh, problem, filename, vector )
      else if ( present(sysvector) ) then
        call printtofile_all ( mesh, problem, filename, sysvector=sysvector )
      else
!       coordinates
        call printtofile_all ( mesh, problem, filename )
      end if

    end if

  end subroutine printtofile


! print a sysvector or vector

  subroutine printtofile_all ( mesh, problem, filename, vector, sysvector )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the data to
    character (len=*), intent(in) :: filename

    type(vector_t), intent(in), optional :: vector
    type(sysvector_t), intent(in), optional :: sysvector

    logical, parameter :: includenode = .true.
    integer :: n1, n2, nodenr, node, elem, elgrp, n
    real(dp) :: x(mesh%ndim)

    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in printtofile_all: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    n = mesh%ndim

    if ( present(vector) ) then

!     vector

      open ( unit=unit_printtofile, file=filename, status='replace', recl=300 )

      if ( vector%elementwise ) then
        n2 = 0
        do elgrp = 1, mesh%nelgrp
          do elem = 1, mesh%grpnumel(elgrp)
            do node = 1, mesh%elnumnod(elgrp)
              n1 = n2
              n2 = n1 + problem%vec_elnumdegfd(elgrp)%a(node,vector%vec)
              nodenr = mesh%topology(elgrp)%a(node,elem)
              x = mesh%coor(nodenr,:)
              if ( n2 > n1 .and. &
                       all(x>=boxmin(:n)) .and. all(x<=boxmax(:n)) ) then
                call write_line_of_data ( includenode, nodenr, x, &
                  vector%u(n1+1:n2) )
              end if
            end do
            write(unit=unit_printtofile,fmt=*)
          end do
        end do
      else
        do nodenr = 1, mesh%nnodes
          n1 = problem%vec_nodnumdegfd ( nodenr, vector%vec )
          n2 = problem%vec_nodnumdegfd ( nodenr + 1, vector%vec )
          x = mesh%coor(nodenr,:)
          if ( n2 > n1 .and. all(x>=boxmin(:n)) .and. all(x<=boxmax(:n)) ) then
            call write_line_of_data ( includenode, nodenr, x, &
              vector%u(n1+1:n2) )
          end if
        end do
      end if

      close(unit=unit_printtofile)

    else if ( present(sysvector) ) then

!     sysvector

      open ( unit=unit_printtofile, file=filename, status='replace', recl=300 )

      do nodenr = 1, mesh%nnodes
        n1 = problem%nodnumdegfd ( nodenr )
        n2 = problem%nodnumdegfd ( nodenr + 1 )
        x = mesh%coor(nodenr,:)
        if ( n2 > n1 .and. all(x>=boxmin(:n)) .and. all(x<=boxmax(:n)) ) then
         call write_line_of_data ( includenode, nodenr, x, &
           sysvector%u(problem%degfdperm(n1+1:n2,2)) )
        end if
      end do

      close(unit=unit_printtofile)

    else

!     Coordinates

      open ( unit=unit_printtofile, file=filename, status='replace', recl=300 )

      do nodenr = 1, mesh%nnodes
        x = mesh%coor(nodenr,:)
        if ( all(x>=boxmin(:n)) .and. all(x<=boxmax(:n)) ) then
         call write_line_of_data ( includenode, nodenr, x )
        end if
      end do

      close(unit=unit_printtofile)

    end if

  end subroutine printtofile_all


! print a sysvector or vector for plotting along a curve with gnuplot

  subroutine printtofile_geometry ( mesh, problem, geometry, filename, vector, &
    sysvector )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   curve number for plotting
    type(geometry_t), intent(in) :: geometry

!   the filename for writing the data to
    character (len=*), intent(in) :: filename

    type(vector_t), intent(in), optional :: vector
    type(sysvector_t), intent(in), optional :: sysvector


    logical, parameter :: includenode = .true.
    integer :: n1, n2, nodenr, node, n
    real(dp) :: x(mesh%ndim)

    n = mesh%ndim

    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in printtofile_geometry: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(vector) ) then

!     vector

      if ( vector%elementwise ) then
        write(*,'(/a/a/)') &
          'Error in printtofile_geometry: printing of data in vectors stored', &
          ' elementwise not yet implemented.'
        stop
      end if

      open ( unit=unit_printtofile, file=filename, status='replace', recl=300 )

      do node = 1, geometry%nnodes
        nodenr = geometry%nodes(node)
        x = mesh%coor(nodenr,:)
        n1 = problem%vec_nodnumdegfd ( nodenr, vector%vec )
        n2 = problem%vec_nodnumdegfd ( nodenr + 1, vector%vec )
        if ( n2 > n1 .and. all(x>=boxmin(:n)) .and. all(x<=boxmax(:n)) ) then
          call write_line_of_data ( includenode, nodenr, x, &
            vector%u(n1+1:n2) )
        end if
      end do

      close(unit=unit_printtofile)

    else if ( present(sysvector) ) then

!     sysvector

      open ( unit=unit_printtofile, file=filename, status='replace', recl=300 )

      do node = 1, geometry%nnodes
        nodenr = geometry%nodes(node)
        n1 = problem%nodnumdegfd ( nodenr )
        n2 = problem%nodnumdegfd ( nodenr + 1 )
        x = mesh%coor(nodenr,:)
        if ( n2 > n1 .and. all(x>=boxmin(:n)) .and. all(x<=boxmax(:n)) ) then
          call write_line_of_data ( includenode, nodenr, x, &
            sysvector%u(problem%degfdperm(n1+1:n2,2)) )
        end if
      end do

      close(unit=unit_printtofile)

    else if ( present(sysvector) ) then

!     coordinates only

      open ( unit=unit_printtofile, file=filename, status='replace', recl=300 )

      do node = 1, geometry%nnodes
        nodenr = geometry%nodes(node)
        x = mesh%coor(nodenr,:)
        if ( all(x>=boxmin(:n)) .and. all(x<=boxmax(:n)) ) then
          call write_line_of_data ( includenode, nodenr, x )
        end if
      end do

      close(unit=unit_printtofile)

    end if

  end subroutine printtofile_geometry


! print a sysvector or vector for plotting with printtofile (2D domain)

  subroutine printtofile_2dgrid ( mesh, problem, nnx, nny, filename, vector, &
    sysvector, surface )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   nnx and nny are the number of nodes on the regular grid in x and y
!   direction respectively
    integer, intent(in) :: nnx, nny

!   the filename for writing the data to
    character (len=*), intent(in) :: filename

    type(vector_t), intent(in), optional :: vector
    type(sysvector_t), intent(in), optional :: sysvector

!   print on this surface
    integer, intent(in), optional :: surface


!   The purpose of this routine is to write the data to a file in regular
!   square pattern for further processing with matlab, gnuplot etc.
!   This only works if the mesh is made with the internal meshgenerator for
!   quadrilateral region of TFEM.


    logical, parameter :: includenode = .false.
    integer :: i, j, n1, n2, nodenr, nn2, elnumnod, nnodes
    logical :: wdata
    real(dp) :: x(mesh%ndim)


    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in printtofile_2dgrid: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(surface) ) then
      elnumnod = mesh%surfaces(surface)%elnumnod
      nnodes = mesh%surfaces(surface)%nnodes
    else
      elnumnod = mesh%elnumnod(1) ! expecting only one group
      nnodes = mesh%nnodes
    end if

    if ( elnumnod == 5 ) then
      nn2 = nnx * nny + (nnx-1)*(nny-1)
    else
      nn2 = nnx * nny
    end if

    if ( nnodes /= nn2 ) then
      write(*,'(/a/a/)') &
        'Error in printtofile_2dgrid: number of nodes in the mesh is ', &
        ' not compatible with a grid of nnx * nny.'
      stop
    end if

    if ( present(vector) ) then

!     vector

      if ( vector%elementwise ) then
        write(*,'(/a/a/)') &
          'Error in printtofile_2dgrid: printing of data in vectors stored', &
          ' elementwise not yet implemented.'
        stop
      end if

      open ( unit=unit_printtofile, file=filename, status='replace', recl=300 )

      do j = 1, nny
        wdata = .false.
        do i = 1, nnx
          if ( elnumnod == 5 ) then
            nodenr = i + ( 2 * nnx - 1 ) * ( j - 1 )
          else
            nodenr = i + nnx * ( j - 1 )
          end if
          if ( present(surface) ) nodenr = mesh%surfaces(surface)%nodes(nodenr)
          x = mesh%coor(nodenr,:)
          n1 = problem%vec_nodnumdegfd ( nodenr, vector%vec )
          n2 = problem%vec_nodnumdegfd ( nodenr + 1, vector%vec )
          if ( n2 > n1 ) then
            call write_line_of_data ( includenode, nodenr, x, &
              vector%u(n1+1:n2) )
            wdata = .true.
          end if
        end do
        if ( wdata .and. blank_line ) write(unit=unit_printtofile,fmt=*)
      end do

      close(unit=unit_printtofile)

    else if ( present(sysvector) ) then

!     sysvector

      open ( unit=unit_printtofile, file=filename, status='replace', recl=300 )

      do j = 1, nny
        wdata = .false.
        do i = 1, nnx
          if ( elnumnod == 5 ) then
            nodenr = i + ( 2 * nnx - 1 ) * ( j - 1 )
          else
            nodenr = i + nnx * ( j - 1 )
          end if
          if ( present(surface) ) nodenr = mesh%surfaces(surface)%nodes(nodenr)
          x = mesh%coor(nodenr,:)
          n1 = problem%nodnumdegfd ( nodenr )
          n2 = problem%nodnumdegfd ( nodenr + 1 )
          if ( n2 > n1 ) then
            call write_line_of_data ( includenode, nodenr, x, &
              sysvector%u(problem%degfdperm(n1+1:n2,2)) )
            wdata = .true.
          end if
        end do
        if ( wdata .and. blank_line ) write(unit=unit_printtofile,fmt=*)
      end do

      close(unit=unit_printtofile)

    end if

  end subroutine printtofile_2dgrid


! print a sysvector or vector for plotting with printtofile (2D domain)

  subroutine printtofile_3dgrid ( mesh, problem, nnx, nny, nnz, filename, &
    vector, sysvector )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   nnx, nny and nnz are the number of nodes on the regular grid in x, y and z
!   direction respectively
    integer, intent(in) :: nnx, nny, nnz

!   the filename for writing the data to
    character (len=*), intent(in) :: filename

    type(vector_t), intent(in), optional :: vector
    type(sysvector_t), intent(in), optional :: sysvector


!   The purpose of this routine is to write the data to a file in regular
!   square pattern for further processing with matlab, gnuplot etc.
!   This only works if the mesh is made with the internal meshgenerator for
!   quadrilateral region of TFEM.


    logical, parameter :: includenode = .false.
    integer :: i, j, n1, n2, nodenr, k
    logical :: wdata
    real(dp) :: x(mesh%ndim)


    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in printtofile_3dgrid: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( mesh%nnodes /= nnx * nny * nnz ) then
      write(*,'(/a/a/)') &
        'Error in printtofile_3dgrid: number of nodes in the mesh is ', &
        ' not compatible with a grid of nnx * nny *nnz.'
      stop
    end if

    if ( present(vector) ) then

!     vector

      if ( vector%elementwise ) then
        write(*,'(/a/a/)') &
          'Error in printtofile_3dgrid: printing of data in vectors stored', &
          ' elementwise not yet implemented.'
        stop
      end if

      open ( unit=unit_printtofile, file=filename, status='replace', recl=300 )

      do k = 1, nnz
        do j = 1, nny
          wdata = .false.
          do i = 1, nnx
            nodenr = i + nnx * ( j - 1 ) + nnx * nny * ( k - 1 )
            x = mesh%coor(nodenr,:)
            n1 = problem%vec_nodnumdegfd ( nodenr, vector%vec )
            n2 = problem%vec_nodnumdegfd ( nodenr + 1, vector%vec )
            if ( n2 > n1 ) then
              call write_line_of_data ( includenode, nodenr, x, &
                vector%u(n1+1:n2) )
              wdata = .true.
            end if
          end do
          if ( wdata .and. blank_line ) write(unit=unit_printtofile,fmt=*)
        end do
      end do

      close(unit=unit_printtofile)

    else if ( present(sysvector) ) then

!     sysvector

      open ( unit=unit_printtofile, file=filename, status='replace', recl=300 )

      do k = 1, nnz
        do j = 1, nny
          wdata = .false.
          do i = 1, nnx
            nodenr = i + nnx * ( j - 1 ) + nnx * nny * ( k - 1 )
            x = mesh%coor(nodenr,:)
            n1 = problem%nodnumdegfd ( nodenr )
            n2 = problem%nodnumdegfd ( nodenr + 1 )
            if ( n2 > n1 ) then
              call write_line_of_data ( includenode, nodenr, x, &
                sysvector%u(problem%degfdperm(n1+1:n2,2)) )
              wdata = .true.
            end if
          end do
          if ( wdata .and. blank_line ) write(unit=unit_printtofile,fmt=*)
        end do
      end do

      close(unit=unit_printtofile)

    end if

  end subroutine printtofile_3dgrid


! write line of data

  subroutine write_line_of_data ( includenode, nodenr, x, u )

    logical, intent(in) :: includenode
    integer, intent(in) :: nodenr
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(in), optional :: u

    if ( present(u) ) then
!     print vector data
      if ( include_node_number .and. includenode ) then
        write ( unit=unit_printtofile, fmt=* ) nodenr, x, u
      else
        write ( unit=unit_printtofile, fmt=* ) x, u
      end if
    else
!     print coordinates
      write ( unit=unit_printtofile, fmt=* ) nodenr, x
    end if

  end subroutine write_line_of_data

end module printtofile_m
