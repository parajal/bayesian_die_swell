
! Copyright (C) 2004-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Basic meshgeneration:
!   - mesh generation of regular meshes on simple domains
!
! mesh on a quadrilateral


module meshgen_basic_quadrilateral_m

  use glob_defs_m
  use meshgen_basic_common_m

  implicit none

contains


! Simple meshgenerator for 2D quadrilateral regions.

  subroutine quadrilateral2d ( mesh, mesh_options, func )

    type(mesh_t), intent(inout) :: mesh

!   Options for generating the mesh.
!   See type definition for possibilities and defaults.
    type(meshgen_options_t), intent(in) :: mesh_options

!   This function must give the coordinates along a curved boundary as a
!   function of the local reference coordinate xr \in [0,1]. Note, that the
!   reference coordinate xr is in the positive direction of the curves as
!   defined in the meshgenerator, i.e.:
!        curve 1: xr=xi
!        curve 2: xr=eta
!        curve 3: xr=1-xi
!        curve 4: xr=1-eta
!   where (xi,eta) are local reference coordinates in the unit square as
!   evaluated on the boundary.
!   The region [0,1]x[0,1] will be mapped onto a curved region. All nodes will
!   will be mapped and thus affected, not just the boundary nodes.
!   The mapping method used is the Blending function method of Gordon and Hall
!   for quadrilaterals as described in:
!      B. Szabo and I. Babuska: "Finite Element Analysis", Wiley 1991.

    optional :: func
    interface
      function func ( nr, xr )
        use kind_defs_m
        implicit none
        integer, intent(in) :: nr
        real(dp), intent(in) :: xr
        real(dp), dimension(2) :: func
      end function func
    end interface



!   A basic mesh structure is generated (topology, coordinates etc). Other
!   parts of the mesh structure still need to be filled.


    integer :: node, curve
    real(dp) :: xi, eta
    real(dp), allocatable, dimension(:,:) :: tmp

!   test mesh

    if ( mesh%meshgen ) then
      write(*,'(/3(a/))') 'Error quadrilateral2d:', &
        ' mesh is not empty', &
        ' either use a new mesh_t variable or delete old mesh '
      stop
    end if

!   test mesh_options

    if ( mesh_options%rx < 1 .or. mesh_options%ry < 1 ) then
      write(*,'(/2(a/))') 'Error quadrilateral2d:', &
        ' rx or ry in mesh_options must not be smaller than 1'
      stop
    end if

!   fill mesh

    select case ( mesh_options%elshape )

      case(-1) ! element shape probably not set

        write(*,'(/a/)') 'Error in quadrilateral2d: element shape not set'
        stop

      case(3) ! three-node triangle

        call rectangle2d_tria3 ( mesh, mesh_options%nx, mesh_options%ny, &
          mesh_options%rx, mesh_options%ry, mesh_options%ratio(1:4), &
          mesh_options%factor(1:4) )

      case(4) ! six-node triangle

        call rectangle2d_tria6 ( mesh, mesh_options%nx, mesh_options%ny, &
          mesh_options%rx, mesh_options%ry, mesh_options%ratio(1:4), &
          mesh_options%factor(1:4) )

      case(5) ! four-node quadrilateral

        call rectangle2d_quad4 ( mesh, mesh_options%nx, mesh_options%ny, &
          mesh_options%rx, mesh_options%ry, mesh_options%ratio(1:4), &
          mesh_options%factor(1:4) )

      case(6) ! nine-node quadrilateral

        call rectangle2d_quad9 ( mesh, mesh_options%nx, mesh_options%ny, &
          mesh_options%rx, mesh_options%ry, mesh_options%ratio(1:4), &
          mesh_options%factor(1:4) )

      case(7) ! seven-node triangle

        call rectangle2d_tria7 ( mesh, mesh_options%nx, mesh_options%ny, &
          mesh_options%rx, mesh_options%ry, mesh_options%ratio(1:4), &
          mesh_options%factor(1:4) )

      case(9) ! five-node quadrilateral

        call rectangle2d_quad5 ( mesh, mesh_options%nx, mesh_options%ny, &
          mesh_options%rx, mesh_options%ry, mesh_options%ratio(1:4), &
          mesh_options%factor(1:4) )

      case(10) ! four-node triangle

        call rectangle2d_tria4 ( mesh, mesh_options%nx, mesh_options%ny, &
          mesh_options%rx, mesh_options%ry, mesh_options%ratio(1:4), &
          mesh_options%factor(1:4) )

      case(102) ! 2D high-order quad with equal order interpolation

        if ( mesh_options%p < 1 ) then
          write(*,'(/2(a/))') 'Error quadrilateral2d:', &
            ' for high-order elements p in mesh_options must be larger than 0'
          stop
        end if

        if ( .not. any( mesh_options%l == [ 0, 1 ] ) ) then
          write(*,'(/2(a/))') 'Error quadrilateral2d:', &
            ' for high-order elements l in mesh_options must be either 0 or 1.'
          stop
        end if

        call rectangle2d_quad_ho ( mesh, mesh_options%nx, mesh_options%ny, &
          mesh_options%rx, mesh_options%ry, mesh_options%ratio(1:4), &
          mesh_options%factor(1:4), mesh_options%p, mesh_options%l )

      case(104) ! 2D high-order triangle

        if ( mesh_options%p < 1 ) then
          write(*,'(/2(a/))') 'Error quadrilateral2d:', &
            ' for high-order elements p in mesh_options must be larger than 0'
          stop
        end if

        call rectangle2d_triangle_ho ( mesh, mesh_options%nx, mesh_options%ny, &
          mesh_options%rx, mesh_options%ry, mesh_options%ratio(1:4), &
          mesh_options%factor(1:4), mesh_options%p )

      case default

        write(*,'(/a,i0/)') 'Error: Wrong element shape: ', mesh_options%elshape
        stop

    end select

!   global element shape

    select case ( mesh_options%elshape )

      case(3,4,7,10,104)

        mesh%element(:)%globalshape = 'triangle'

      case(5,6,9,30,102)

        mesh%element(:)%globalshape = 'quadrilateral'

      case default

        call errormsg_case_default ( 'quadrilateral2d', &
          'mesh_options%elshape', int_value=mesh_options%elshape )

    end select


!   choose region

    select case ( mesh_options%regionshape )

      case(1) ! rectangle

        if ( any(mesh_options%curved) ) then
          write(*,'(/a/)') &
            'Warning in quadrilateral2d: curved boundaries need regionshape=3'
        end if

!       scale coordinates

        mesh%coor(:,1) = mesh%coor(:,1) * mesh_options%lx + mesh_options%ox
        mesh%coor(:,2) = mesh%coor(:,2) * mesh_options%ly + mesh_options%oy

      case(2) ! quadrilateral

        if ( any(mesh_options%curved) ) then
          write(*,'(/a/)') &
            'Warning in quadrilateral2d: curved boundaries need regionshape=3'
        end if

        allocate ( tmp(mesh%nnodes,2) )

        call straight_quadrilateral ( tmp )

        mesh%coor = tmp

        deallocate ( tmp )

      case(3) ! curved quadrilateral

        if ( .not. present(func) ) then
          write(*,'(/a/)') &
            'Error in quadrilateral2d: func not present for curved boundaries.'
          stop
        end if

        do curve = 1, 4
          if ( mesh_options%curved(curve) .and. &
               mesh_options%funcnr(curve) == -1 ) then
            write(*,'(/a,i0/)') &
              'Error in quadrilateral2d: function number not set for curve ', &
              curve
            stop
          end if
        end do

        allocate ( tmp(mesh%nnodes,2) )

!       start with straight region

        call straight_quadrilateral ( tmp )

!       curve 1

        if ( mesh_options%curved(1) ) then

!         curve 1 is curved

          do node = 1, mesh%nnodes
             xi  = mesh%coor(node,1)
             eta = mesh%coor(node,2)
             tmp(node,:) = tmp(node,:) + ( func(mesh_options%funcnr(1),xi) &
               - (1-xi) * mesh_options%x2d(1,:) &
               - xi * mesh_options%x2d(2,:) ) * (1-eta)
          end do

        end if

        if ( mesh_options%curved(2) ) then

!         curve 2 is curved

          do node = 1, mesh%nnodes
             xi  = mesh%coor(node,1)
             eta = mesh%coor(node,2)
             tmp(node,:) = tmp(node,:) + ( func(mesh_options%funcnr(2),eta) &
               - (1-eta) * mesh_options%x2d(2,:) &
               - eta * mesh_options%x2d(3,:) ) * xi
          end do

        end if

        if ( mesh_options%curved(3) ) then

!         curve 3 is curved

          do node = 1, mesh%nnodes
             xi  = mesh%coor(node,1)
             eta = mesh%coor(node,2)
             tmp(node,:) = tmp(node,:) + ( func(mesh_options%funcnr(3),1-xi) &
               - (1-xi) * mesh_options%x2d(4,:) &
               - xi * mesh_options%x2d(3,:) ) * eta
          end do

        end if

        if ( mesh_options%curved(4) ) then

!         curve 4 is curved

          do node = 1, mesh%nnodes
             xi  = mesh%coor(node,1)
             eta = mesh%coor(node,2)
             tmp(node,:) = tmp(node,:) + ( func(mesh_options%funcnr(4),1-eta) &
               - (1-eta) * mesh_options%x2d(1,:) &
               - eta * mesh_options%x2d(4,:) ) * (1-xi)
          end do

        end if

        mesh%coor = tmp

        deallocate ( tmp )

      case default

        write(*,'(/a,i0/)') &
         'Error: Wrong region shape: ', mesh_options%regionshape
        stop

    end select

!   type of elements on curves

    mesh%curves(1:mesh%ncurves)%element%globalshape = 'line'

!   allow surfaces

    allocate( mesh%surfaces(MAXSURFACES) )

!   no volumes

    allocate( mesh%volumes(0) )

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

  contains

    subroutine straight_quadrilateral ( coor )

      real(dp), dimension(:,:), intent(out) :: coor

!     x = [ x1 * (1-xi) + x2 * xi ] (1 - eta) +
!         [ x4 * (1-xi) + x3 * xi ] eta
!     with xi = coor(1), eta = coor(2) \in [0,1] from rectangle2d

      integer :: i

      do i = 1, 2
        coor(:,i) = &
         ( mesh_options%x2d(1,i) * (1 - mesh%coor(:,1)) + &
           mesh_options%x2d(2,i) * mesh%coor(:,1) ) * (1 - mesh%coor(:,2)) + &
         ( mesh_options%x2d(4,i) * ( 1 - mesh%coor(:,1) ) + &
           mesh_options%x2d(3,i) * mesh%coor(:,1) ) * mesh%coor(:,2)
      end do

    end subroutine straight_quadrilateral

  end subroutine quadrilateral2d


! Generate mesh on region [0,1]x[0,1] using triangular elements with 3 nodes.

  subroutine rectangle2d_tria3 ( mesh, n, m, r, s, ratio, factor )

    type(mesh_t), intent(out) :: mesh

!   number of elements along the sides and refinement factors
!   The ACTUAL number of elements is n*r x m*s.
    integer, intent(in) :: n, m, r, s

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, dimension(:), intent(in) :: ratio
    real(dp), dimension(:), intent(in) :: factor

!   Generate a mesh on region  [0,1]x[0,1] using triangles with 3 nodes.
!   The numbering for both nodal points and elements is row by row.
!   For example a 2x2 mesh is numbered as follows
!
!   Nodes:
!
!       7 ------- 8 ------- 9
!       |       / |       / |
!       |    /    |    /    |
!       | /       | /       |
!       4 ------- 5 ------- 6
!       |       / |       / |
!       |    /    |    /    |
!       | /       | /       |
!       1 ------- 2 ------- 3
!
!   Elements:
!
!       x ------- x ------- x
!       |  5    / |  7    / |
!       |    /    |    /    |
!       | /     6 | /    8  |
!       x ------- x ------- x
!       |  1    / |  3    / |
!       |    /    |    /    |
!       | /    2  | /    4  |
!       x ------- x ------- x
!
!  Also generated are four points and four curves:
!
!               <-----
!                  C3
!      P4 --------------------- P3
!       |                       |           Note:   C1=P1-P2
!       |                       |                   C2=P2-P3
!  |    |                       |    ^              C3=P3-P4
!  | C4 |                       | C2 |              C4=P4-P1
! \|/   |                       |    |
!       |                       |
!       |                       |
!      P1 --------------------- P2
!                  C1
!                ----->
!

    integer, parameter :: inpelm = 3 ! number of nodal points per element
    integer, parameter :: inpelmc = 2 ! number of nodal points per element on
                                      ! curves

    integer :: elnodes(inpelm), i, j, e, k
    integer :: nn1, nn2, nn1row, nn1col, node, elem, curve, nr, mr
    real(dp) :: deltax, deltay, D
    real(dp) :: deltax1, deltax2, deltax3, deltax4
    real(dp), allocatable, dimension(:) :: x1, x2, x3, x4

    allocate ( x1(1:n*r+1), x2(1:m*s+1), x3(1:n*r+1), x4(1:m*s+1) )


    nr = n*r
    mr = m*s

    nn1row = nr+1
    nn1col = mr+1

!   internal mesh

    mesh%ndim   = 2
    mesh%nelem  = 2*nr*mr
    mesh%nnodes = nn1row*nn1col
    mesh%nelgrp = 1

!   allocate memory for internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%nnodes_blend(2) ); mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology(1) )
    allocate( mesh%topology(1)%a(inpelm,mesh%nelem) )
    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   element type

    mesh%element(1)%elshape = 3
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 2

!   topology

    mesh%grpnumel(1) = mesh%nelem

    do i = 1, nr
      do j = 1, mr

        nn1 = (j-1)*nn1row + i-1 ! number of nodes before element (i,j)
        nn2 = nn1 + nn1row   ! number of nodes before element (i,j) + 1 row

        e = 2*(i-1) + 2*(j-1)*nr + 1   ! element number

        elnodes = [ nn1+1, (nn2+i,i=2,1,-1) ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

        e = e+1 ! element number

        elnodes = [ (nn1+i, i=1,2), nn2+2 ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

      end do
    end do

!   coordinates


    if ( all(ratio == 0) ) then
!     equidistant
      deltax = 1.0_dp / nr
      deltay = 1.0_dp / mr
      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*nn1row
          mesh%coor(node,1) = (i-1) * deltax
          mesh%coor(node,2) = (j-1) * deltay
        end do
      end do
    else
!     non-equidistant
      call distribute_elements ( n, x1(::r), ratio(1), factor(1) )
      call distribute_elements ( m, x2(::s), ratio(2), factor(2) )
      call distribute_elements ( n, x3(::r), ratio(3), factor(3) )
      call distribute_elements ( m, x4(::s), ratio(4), factor(4) )
!     refine
      do elem = 1, n
        k = ( elem - 1 ) * r  ! nodes before element elem
        deltax1 = (x1(k+r+1)-x1(k+1))/r
        deltax3 = (x3(k+r+1)-x3(k+1))/r
        do node = k+2, k+r
          x1(node) = x1(node-1) + deltax1
          x3(node) = x3(node-1) + deltax3
        end do
      end do
      do elem = 1, m
        k = ( elem - 1 ) * s  ! nodes before element elem
        deltax2 = (x2(k+s+1)-x2(k+1))/s
        deltax4 = (x4(k+s+1)-x4(k+1))/s
        do node = k+2, k+s
          x2(node) = x2(node-1) + deltax2
          x4(node) = x4(node-1) + deltax4
        end do
      end do
      x3 = 1._dp - x3(nn1row:1:-1)
      x4 = 1._dp - x4(nn1col:1:-1)
!     create straight lines in reference square [0,1]x[0,1]
      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*nn1row
          D = ( x1(i) - x3(i) ) * ( x4(j) - x2(j) ) - 1
          mesh%coor(node,1) = ( x4(j) * (x1(i) - x3(i)) - x1(i) ) / D
          mesh%coor(node,2) = ( x1(i) * (x4(j) - x2(j)) - x4(j) ) / D
        end do
      end do
    end if

    deallocate ( x1, x2, x3, x4 )

!   points

    mesh%npoints = 4

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn1row
    mesh%points(3) = nn1row * nn1col
    mesh%points(4) = nn1row * ( nn1col - 1 ) + 1

!   curves

    mesh%ncurves = 4

!   allocate memory for curves

    allocate( mesh%curves(MAXCURVES) )
    allocate( mesh%curves(1)%nodes(nn1row) )
    allocate( mesh%curves(2)%nodes(nn1col) )
    allocate( mesh%curves(3)%nodes(nn1row) )
    allocate( mesh%curves(4)%nodes(nn1col) )
    allocate( mesh%curves(1)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(2)%topology(inpelmc,mr,2) )
    allocate( mesh%curves(3)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(4)%topology(inpelmc,mr,2) )

!   element type

    mesh%curves(:4)%element%elshape = 1
    mesh%curves(:4)%element%numnod = inpelmc
    mesh%curves(:4)%element%ndim = 2


!   other info of curves

    mesh%curves(:4)%ndim = 2
    mesh%curves(:4)%nblend = 0
    mesh%curves(1:4)%nnodes = [ nn1row, nn1col, nn1row, nn1col ]
    mesh%curves(1:4)%nelem = [ nr, mr, nr, mr ]
    mesh%curves(1:4)%n = mesh%curves(1:4)%nelem
    do i = 1, 4
      allocate(mesh%curves(i)%element_blend(0))
      allocate(mesh%curves(i)%nnodes_blend(2) )
      mesh%curves(i)%nnodes_blend = [0,mesh%curves(i)%nnodes]
    end do

!   nodes of curves

    mesh%curves(1)%nodes = [ (i, i=1,nn1row) ]                      ! curve 1
    mesh%curves(2)%nodes = [ (j*nn1row, j=1,nn1col) ]               ! curve 2
    mesh%curves(3)%nodes = [ ((nn1col-1)*nn1row+i, i=nn1row,1,-1) ] ! curve 3
    mesh%curves(4)%nodes = [ ((j-1)*nn1row+1, j=nn1col,1,-1) ]      ! curve 4

!   topology of elements on curves

    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        mesh%curves(curve)%topology(:,elem,1) = [ (elem-1+i,i=1,2) ]
        mesh%curves(curve)%topology(:,elem,2) = &
           mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
      end do
    end do

  end subroutine rectangle2d_tria3


! Generate mesh on region [0,1]x[0,1] using triangular elements with 4 nodes.

  subroutine rectangle2d_tria4 ( mesh, n, m, r, s, ratio, factor )

    type(mesh_t), intent(out) :: mesh

!   number of elements along the sides and refinement factors
!   The ACTUAL number of elements is n*r x m*s.
    integer, intent(in) :: n, m, r, s

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, dimension(:), intent(in) :: ratio
    real(dp), dimension(:), intent(in) :: factor

!   Generate a mesh on region  [0,1]x[0,1] using triangles with 4 nodes.
!   The numbering for both nodal points and elements is row by row.
!   For example a 2x2 mesh is numbered as follows
!
!   Nodes:
!
!       15------- 16------- 17
!       | 11    / | 13    / |
!       |    /    |    /    |
!       | /    12 | /    14 |
!       8 ------- 9 ------- 10
!       |  4    / | 6     / |
!       |    /    |    /    |
!       | /    5  | /     7 |
!       1 ------- 2 ------- 3
!
!   Elements:
!
!       x ------- x ------- x
!       |  5    / |  7    / |
!       |    /    |    /    |
!       | /     6 | /    8  |
!       x ------- x ------- x
!       |  1    / |  3    / |
!       |    /    |    /    |
!       | /    2  | /    4  |
!       x ------- x ------- x
!
!  Also generated are four points and four curves:
!
!               <-----
!                  C3
!      P4 --------------------- P3
!       |                       |           Note:   C1=P1-P2
!       |                       |                   C2=P2-P3
!  |    |                       |    ^              C3=P3-P4
!  | C4 |                       | C2 |              C4=P4-P1
! \|/   |                       |    |
!       |                       |
!       |                       |
!      P1 --------------------- P2
!                  C1
!                ----->
!

    integer, parameter :: inpelm = 4 ! number of nodal points per element
    integer, parameter :: inpelmc = 2 ! number of nodal points per element on
                                      ! curves

    integer :: elnodes(inpelm), i, j, e, k
    integer :: nn1, nn2, nn3, nn1row, nn1col, node, elem, curve, nr, mr
    integer :: node1, node3
    real(dp) :: deltax, deltay, D
    real(dp) :: deltax1, deltax2, deltax3, deltax4
    real(dp), allocatable, dimension(:) :: x1, x2, x3, x4

    allocate ( x1(1:n*r+1), x2(1:m*s+1), x3(1:n*r+1), x4(1:m*s+1) )


    nr = n*r
    mr = m*s

    nn1row = nr+1
    nn1col = mr+1

!   internal mesh

    mesh%ndim   = 2
    mesh%nelem  = 2*nr*mr
    mesh%nnodes = nn1row*nn1col + 2*nr*mr
    mesh%nelgrp = 1

!   allocate memory for internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%nnodes_blend(2) ); mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology(1) )
    allocate( mesh%topology(1)%a(inpelm,mesh%nelem) )
    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   element type

    mesh%element(1)%elshape = 10
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 2

!   topology

    mesh%grpnumel(1) = mesh%nelem

    do i = 1, nr
      do j = 1, mr

        nn1 = (j-1)*(nn1row+2*nr) + i-1 ! number of nodes before element (i,j)
        nn2 = nn1 + nn1row + i-1 ! number of nodes before element (i,j) + 1 row
        nn3 = nn1 + nn1row + 2*nr !number of nodes before element (i,j) + 2 rows

        e = 2*(i-1) + 2*(j-1)*nr + 1   ! element number

        elnodes = [ nn1+1, (nn3+i,i=2,1,-1), nn2+1 ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

        e = e+1 ! element number

        elnodes = [ (nn1+i, i=1,2), nn3+2, nn2+2 ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

      end do
    end do

!   coordinates

    if ( all(ratio == 0) ) then
!     equidistant
      deltax = 1.0_dp / nr
      deltay = 1.0_dp / mr
      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*(nn1row+2*nr)
          mesh%coor(node,1) = (i-1) * deltax
          mesh%coor(node,2) = (j-1) * deltay
        end do
      end do
      do i = 1, nr
        do j = 1, mr
          node = 2*i-1 + (j-1)*(nn1row+2*nr) + nn1row
          mesh%coor(node,1) = (i-1) * deltax + deltax/3
          mesh%coor(node,2) = (j-1) * deltay + 2*deltay/3
          node = node + 1
          mesh%coor(node,1) = (i-1) * deltax + 2*deltax/3
          mesh%coor(node,2) = (j-1) * deltay + deltay/3
        end do
      end do
    else
!     non-equidistant
      call distribute_elements ( n, x1(::r), ratio(1), factor(1) )
      call distribute_elements ( m, x2(::s), ratio(2), factor(2) )
      call distribute_elements ( n, x3(::r), ratio(3), factor(3) )
      call distribute_elements ( m, x4(::s), ratio(4), factor(4) )
!     refine
      do elem = 1, n
        k = ( elem - 1 ) * r  ! nodes before element elem
        deltax1 = (x1(k+r+1)-x1(k+1))/r
        deltax3 = (x3(k+r+1)-x3(k+1))/r
        do node = k+2, k+r
          x1(node) = x1(node-1) + deltax1
          x3(node) = x3(node-1) + deltax3
        end do
      end do
      do elem = 1, m
        k = ( elem - 1 ) * s  ! nodes before element elem
        deltax2 = (x2(k+s+1)-x2(k+1))/s
        deltax4 = (x4(k+s+1)-x4(k+1))/s
        do node = k+2, k+s
          x2(node) = x2(node-1) + deltax2
          x4(node) = x4(node-1) + deltax4
        end do
      end do
      x3 = 1._dp - x3(nn1row:1:-1)
      x4 = 1._dp - x4(nn1col:1:-1)
!     create straight lines in reference square [0,1]x[0,1]
      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*(nn1row+2*nr)
          D = ( x1(i) - x3(i) ) * ( x4(j) - x2(j) ) - 1
          mesh%coor(node,1) = ( x4(j) * (x1(i) - x3(i)) - x1(i) ) / D
          mesh%coor(node,2) = ( x1(i) * (x4(j) - x2(j)) - x4(j) ) / D
        end do
      end do
      do i = 1, nr
        do j = 1, mr
          node = 2*i-1 + (j-1)*(nn1row+2*nr) + nn1row
          node1 = i + (j-1)*(nn1row+2*nr)
          node3 = i + j*(nn1row+2*nr)
          mesh%coor(node,:) = (mesh%coor(node1,:)+ &
                               mesh%coor(node3,:)+mesh%coor(node3+1,:))/3
          node = node + 1
          node3 = node3 + 1
          mesh%coor(node,:) = (mesh%coor(node1,:)+mesh%coor(node1+1,:)+&
                               mesh%coor(node3,:))/3
        end do
      end do
    end if

    deallocate ( x1, x2, x3, x4 )

!   points

    mesh%npoints = 4

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn1row
    mesh%points(3) = nn1row * nn1col + 2*nr*mr
    mesh%points(4) = nn1row * ( nn1col - 1 ) + 2*nr*mr + 1

!   curves

    mesh%ncurves = 4

!   allocate memory for curves

    allocate( mesh%curves(MAXCURVES) )
    allocate( mesh%curves(1)%nodes(nn1row) )
    allocate( mesh%curves(2)%nodes(nn1col) )
    allocate( mesh%curves(3)%nodes(nn1row) )
    allocate( mesh%curves(4)%nodes(nn1col) )
    allocate( mesh%curves(1)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(2)%topology(inpelmc,mr,2) )
    allocate( mesh%curves(3)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(4)%topology(inpelmc,mr,2) )

!   element type

    mesh%curves(:4)%element%elshape = 1
    mesh%curves(:4)%element%numnod = inpelmc
    mesh%curves(:4)%element%ndim = 2

!   other info of curves

    mesh%curves(:4)%ndim = 2
    mesh%curves(:4)%nblend = 0
    mesh%curves(1:4)%nnodes = [ nn1row, nn1col, nn1row, nn1col ]
    mesh%curves(1:4)%nelem = [ nr, mr, nr, mr ]
    mesh%curves(1:4)%n = mesh%curves(1:4)%nelem
    do i = 1, 4
      allocate(mesh%curves(i)%element_blend(0))
      allocate(mesh%curves(i)%nnodes_blend(2) )
      mesh%curves(i)%nnodes_blend = [0,mesh%curves(i)%nnodes]
    end do

!   nodes of curves

    mesh%curves(1)%nodes = [ (i, i=1,nn1row) ]                      ! curve 1
    mesh%curves(2)%nodes = [ (j*nn1row+(j-1)*2*nr, j=1,nn1col) ]    ! curve 2
!                                                                       curve 3
    mesh%curves(3)%nodes = [ ((nn1col-1)*nn1row+2*nr*mr+i, i=nn1row,1,-1) ]
    mesh%curves(4)%nodes = [ ((j-1)*(nn1row+2*nr)+1, j=nn1col,1,-1) ]! curve 4

!   topology of elements on curves

    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        mesh%curves(curve)%topology(:,elem,1) = [ (elem-1+i,i=1,2) ]
        mesh%curves(curve)%topology(:,elem,2) = &
           mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
      end do
    end do

  end subroutine rectangle2d_tria4


! Generate mesh on region [0,1]x[0,1] using triangular elements with 6 nodes.

  subroutine rectangle2d_tria6 ( mesh, n, m, r, s, ratio, factor )

    type(mesh_t), intent(out) :: mesh

!   number of elements along the sides and refinement factors
!   The ACTUAL number of elements is n*r x m*s.
    integer, intent(in) :: n, m, r, s

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, dimension(:), intent(in) :: ratio
    real(dp), dimension(:), intent(in) :: factor

!   Generate a mesh on region  [0,1]x[0,1] using triangles with 6 nodes.
!   The numbering for both nodal points and elements is row by row.
!   For example a 2x2 mesh is numbered as follows
!
!   Nodes:
!
!      21 -- 22 -- 23 -- 24 -- 25
!       |        /  |        /  |
!      16    17    18    19    20
!       | /         | /         |
!      11 -- 12 -- 13 -- 14 -- 15
!       |        /  |        /  |
!       6     7     8     9    10
!       | /         | /         |
!       1 --  2 --  3 --  4 --  5
!
!   Elements:
!
!       x --  x --  x --  x --  x
!       | 5      /  | 7      /  |
!       x     x     x     x     x
!       | /      6  | /      8  |
!       x --  x --  x --  x --  x
!       |  1     /  |  3     /  |
!       x     x     x     x     x
!       | /      2  | /      4  |
!       x --  x --  x --  x --  x
!
!  Also generated are four points and four curves:
!
!               <-----
!                  C3
!      P4 --------------------- P3
!       |                       |           Note:   C1=P1-P2
!       |                       |                   C2=P2-P3
!  |    |                       |    ^              C3=P3-P4
!  | C4 |                       | C2 |              C4=P4-P1
! \|/   |                       |    |
!       |                       |
!       |                       |
!      P1 --------------------- P2
!                  C1
!                ----->
!

    integer, parameter :: inpelm = 6 ! number of nodal points per element
    integer, parameter :: inpelmc = 3 ! number of nodal points per element on
                                      ! curves

    integer :: elnodes(inpelm), i, j, e, k
    integer :: nn1, nn2, nn3, nn1row, nn1col, node, elem, curve, nr, mr
    real(dp) :: deltax, deltay, D
    real(dp) :: deltax1, deltax2, deltax3, deltax4
    real(dp), allocatable, dimension(:) :: x1, x2, x3, x4

    allocate ( x1(1:2*n*r+1), x2(1:2*m*s+1), x3(1:2*n*r+1), x4(1:2*m*s+1) )


    nr = n*r
    mr = m*s

    nn1row = 2*nr+1
    nn1col = 2*mr+1

!   internal mesh

    mesh%ndim   = 2
    mesh%nelem  = 2*nr*mr
    mesh%nnodes = nn1row*nn1col
    mesh%nelgrp = 1

!   allocate memory for internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%nnodes_blend(2) ); mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology(1) )
    allocate( mesh%topology(1)%a(inpelm,mesh%nelem) )
    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   element type

    mesh%element(1)%elshape = 4
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 2

!   topology

    mesh%grpnumel(1) = mesh%nelem

    do i = 1, nr
      do j = 1, mr

        nn1 = 2*(j-1)*nn1row + 2*(i-1) ! number of nodes before element (i,j)
        nn2 = nn1 + nn1row   ! number of nodes before element (i,j) + 1 row
        nn3 = nn2 + nn1row   ! number of nodes before element (i,j) + 2 rows

        e = 2*i-1 + (j-1)*2*nr   ! element number

        elnodes = [ nn1+1, nn2+2, (nn3+i,i=3,1,-1), nn2+1 ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

        e = e + 1

        elnodes = [ (nn1+i,i=1,3), nn2+3, nn3+3, nn2+2 ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

      end do
    end do

!   coordinates

    if ( all(ratio == 0) ) then
!     equidistant
      deltax = 0.5_dp / nr
      deltay = 0.5_dp / mr
      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*nn1row
          mesh%coor(node,1) = (i-1) * deltax
          mesh%coor(node,2) = (j-1) * deltay
        end do
      end do
    else
!     non-equidistant
      call distribute_elements ( n, x1(1:nn1row:2*r), ratio(1), factor(1) )
      call distribute_elements ( m, x2(1:nn1col:2*s), ratio(2), factor(2) )
      call distribute_elements ( n, x3(1:nn1row:2*r), ratio(3), factor(3) )
      call distribute_elements ( m, x4(1:nn1col:2*s), ratio(4), factor(4) )
!     mid-side nodes and refine
      do elem = 1, n
        k = ( elem - 1 ) * 2*r  ! nodes before element elem
        deltax1 = (x1(k+2*r+1)-x1(k+1))/2/r
        deltax3 = (x3(k+2*r+1)-x3(k+1))/2/r
        do node = k+2, k+2*r
          x1(node) = x1(node-1) + deltax1
          x3(node) = x3(node-1) + deltax3
        end do
      end do
      do elem = 1, m
        k = ( elem - 1 ) * 2*s  ! nodes before element elem
        deltax2 = (x2(k+2*s+1)-x2(k+1))/2/s
        deltax4 = (x4(k+2*s+1)-x4(k+1))/2/s
        do node = k+2, k+2*s
          x2(node) = x2(node-1) + deltax2
          x4(node) = x4(node-1) + deltax4
        end do
      end do
      x3 = 1._dp - x3(nn1row:1:-1)
      x4 = 1._dp - x4(nn1col:1:-1)
!     create straight lines in reference square [0,1]x[0,1]
      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*nn1row
          D = ( x1(i) - x3(i) ) * ( x4(j) - x2(j) ) - 1
          mesh%coor(node,1) = ( x4(j) * (x1(i) - x3(i)) - x1(i) ) / D
          mesh%coor(node,2) = ( x1(i) * (x4(j) - x2(j)) - x4(j) ) / D
        end do
      end do
    end if

    deallocate ( x1, x2, x3, x4 )

!   points

    mesh%npoints = 4

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn1row
    mesh%points(3) = nn1row * nn1col
    mesh%points(4) = nn1row * ( nn1col - 1 ) + 1

!   curves

    mesh%ncurves = 4

!   allocate memory for curves

    allocate( mesh%curves(MAXCURVES) )
    allocate( mesh%curves(1)%nodes(nn1row) )
    allocate( mesh%curves(2)%nodes(nn1col) )
    allocate( mesh%curves(3)%nodes(nn1row) )
    allocate( mesh%curves(4)%nodes(nn1col) )
    allocate( mesh%curves(1)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(2)%topology(inpelmc,mr,2) )
    allocate( mesh%curves(3)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(4)%topology(inpelmc,mr,2) )

!   element type

    mesh%curves(:4)%element%elshape = 2
    mesh%curves(:4)%element%numnod = inpelmc
    mesh%curves(:4)%element%ndim = 2

!   other info of curves

    mesh%curves(:4)%ndim = 2
    mesh%curves(:4)%nblend = 0
    mesh%curves(1:4)%nnodes = [ nn1row, nn1col, nn1row, nn1col ]
    mesh%curves(1:4)%nelem = [ nr, mr, nr, mr ]
    mesh%curves(1:4)%n = mesh%curves(1:4)%nelem
    do i = 1, 4
      allocate(mesh%curves(i)%element_blend(0))
      allocate(mesh%curves(i)%nnodes_blend(2) )
      mesh%curves(i)%nnodes_blend = [0,mesh%curves(i)%nnodes]
    end do

!   nodes of curves

    mesh%curves(1)%nodes = [ (i, i=1,nn1row) ]                      ! curve 1
    mesh%curves(2)%nodes = [ (j*nn1row, j=1,nn1col) ]               ! curve 2
    mesh%curves(3)%nodes = [ ((nn1col-1)*nn1row+i, i=nn1row,1,-1) ] ! curve 3
    mesh%curves(4)%nodes = [ ((j-1)*nn1row+1, j=nn1col,1,-1) ]      ! curve 4

!   topology of elements on curves

    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        mesh%curves(curve)%topology(:,elem,1) = [ (2*(elem-1)+i,i=1,3) ]
        mesh%curves(curve)%topology(:,elem,2) = &
           mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
      end do
    end do

  end subroutine rectangle2d_tria6


! Generate mesh on region [0,1]x[0,1] using triangular elements with 7 nodes.

  subroutine rectangle2d_tria7 ( mesh, n, m, r, s, ratio, factor )

    type(mesh_t), intent(out) :: mesh

!   number of elements along the sides and refinement factors
!   The ACTUAL number of elements is n*r x m*s.
    integer, intent(in) :: n, m, r, s

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, dimension(:), intent(in) :: ratio
    real(dp), dimension(:), intent(in) :: factor

!   Generate a mesh on region  [0,1]x[0,1] using triangles with 7 nodes.
!   The numbering for both nodal points and elements is row by row.
!   For example a 2x2 mesh is numbered as follows
!
!   Nodes:
!
!      29 -- 30 -- 31 -- 32 -- 25
!       | 20     /  | 22     /  |
!      24    25    26    27    28
!       | /     21  | /     23  |
!      15 -- 16 -- 17 -- 18 -- 19
!       | 6      /  | 8      /  |
!       10   11    12    13    14
!       | /      7  | /      9  |
!       1 --  2 --  3 --  4 --  5
!
!   Elements:
!
!       x --  x --  x --  x --  x
!       | 5      /  | 7      /  |
!       x     x     x     x     x
!       | /      6  | /      8  |
!       x --  x --  x --  x --  x
!       |  1     /  |  3     /  |
!       x     x     x     x     x
!       | /      2  | /      4  |
!       x --  x --  x --  x --  x
!
!  Also generated are four points and four curves:
!
!               <-----
!                  C3
!      P4 --------------------- P3
!       |                       |           Note:   C1=P1-P2
!       |                       |                   C2=P2-P3
!  |    |                       |    ^              C3=P3-P4
!  | C4 |                       | C2 |              C4=P4-P1
! \|/   |                       |    |
!       |                       |
!       |                       |
!      P1 --------------------- P2
!                  C1
!                ----->
!

    integer, parameter :: inpelm = 7 ! number of nodal points per element
    integer, parameter :: inpelmc = 3 ! number of nodal points per element on
                                      ! curves

    integer :: elnodes(inpelm), i, j, e, node1, node3, k
    integer :: nn1, nn2, nn3, nn4, nn1row, nn1col, node, elem, curve, nr, mr
    real(dp) :: deltax, deltay, D
    real(dp) :: deltax1, deltax2, deltax3, deltax4
    real(dp), allocatable, dimension(:) :: x1, x2, x3, x4

    allocate ( x1(1:2*n*r+1), x2(1:2*m*s+1), x3(1:2*n*r+1), x4(1:2*m*s+1) )


    nr = n*r
    mr = m*s

    nn1row = 2*nr+1
    nn1col = 2*mr+1

!   internal mesh

    mesh%ndim   = 2
    mesh%nelem  = 2*nr*mr
    mesh%nnodes = nn1row*nn1col + 2*nr*mr
    mesh%nelgrp = 1

!   allocate memory for internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%nnodes_blend(2) ); mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology(1) )
    allocate( mesh%topology(1)%a(inpelm,mesh%nelem) )
    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   element type

    mesh%element(1)%elshape = 7
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 2

!   topology

    mesh%grpnumel(1) = mesh%nelem

    do i = 1, nr
      do j = 1, mr

!       number of nodes before element (i,j)
        nn1 = (j-1)*( 2*nn1row + 2*nr ) + 2*(i-1)
!       number of nodes before element (i,j) + 1 row
        nn2 = nn1 + nn1row
!       number of nodes before element (i,j) + 2 rows
        nn3 = nn2 + 2*nr
!       number of nodes before element (i,j) + 3 rows
        nn4 = nn3 + nn1row

        e = 2*i-1 + (j-1)*2*nr   ! element number

        elnodes = [ nn1+1, nn3+2, (nn4+i,i=3,1,-1), nn3+1, nn2+1 ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

        e = e + 1

        elnodes = [ (nn1+i,i=1,3), nn3+3, nn4+3, nn3+2, nn2+2 ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

      end do
    end do

!   coordinates

    if ( all(ratio == 0) ) then
!     equidistant
      deltax = 1.0_dp / nr
      deltay = 1.0_dp / mr
      do i = 1, nn1row
        do j = 1, mr+1
          node = i + (j-1)*(2*nn1row+2*nr)
          mesh%coor(node,1) = (i-1) * deltax/2
          mesh%coor(node,2) = (j-1) * deltay
        end do
      end do
      do i = 1, nr
        do j = 1, mr
          node = 2*i-1 + nn1row + (j-1)*(2*nn1row+2*nr)
          mesh%coor(node,1) = (i-1) * deltax + deltax/3
          mesh%coor(node,2) = (j-1) * deltay + 2*deltay/3
          node = node + 1
          mesh%coor(node,1) = (i-1) * deltax + 2*deltax/3
          mesh%coor(node,2) = (j-1) * deltay + deltay/3
        end do
      end do
      do i = 1, nn1row
        do j = 1, mr
          node = i + nn1row+2*nr + (j-1)*(2*nn1row+2*nr)
          mesh%coor(node,1) = (i-1) * deltax/2
          mesh%coor(node,2) = (j-1) * deltay + deltay/2
        end do
      end do
    else
!     non-equidistant
      call distribute_elements ( n, x1(1:nn1row:2*r), ratio(1), factor(1) )
      call distribute_elements ( m, x2(1:nn1col:2*s), ratio(2), factor(2) )
      call distribute_elements ( n, x3(1:nn1row:2*r), ratio(3), factor(3) )
      call distribute_elements ( m, x4(1:nn1col:2*s), ratio(4), factor(4) )
!     mid-side nodes and refine
      do elem = 1, n
        k = ( elem - 1 ) * 2*r  ! nodes before element elem
        deltax1 = (x1(k+2*r+1)-x1(k+1))/2/r
        deltax3 = (x3(k+2*r+1)-x3(k+1))/2/r
        do node = k+2, k+2*r
          x1(node) = x1(node-1) + deltax1
          x3(node) = x3(node-1) + deltax3
        end do
      end do
      do elem = 1, m
        k = ( elem - 1 ) * 2*s  ! nodes before element elem
        deltax2 = (x2(k+2*s+1)-x2(k+1))/2/s
        deltax4 = (x4(k+2*s+1)-x4(k+1))/2/s
        do node = k+2, k+2*s
          x2(node) = x2(node-1) + deltax2
          x4(node) = x4(node-1) + deltax4
        end do
      end do
      x3 = 1._dp - x3(nn1row:1:-1)
      x4 = 1._dp - x4(nn1col:1:-1)
!     create straight lines in reference square [0,1]x[0,1]
      do i = 1, nn1row
        do j = 1, nn1col
          if ( mod(j,2) == 0 ) then
!           even row
            node = i + (j-1)*nn1row + j*nr
          else
!           odd row
            node = i + (j-1)*(nn1row+nr)
          end if
          D = ( x1(i) - x3(i) ) * ( x4(j) - x2(j) ) - 1
          mesh%coor(node,1) = ( x4(j) * (x1(i) - x3(i)) - x1(i) ) / D
          mesh%coor(node,2) = ( x1(i) * (x4(j) - x2(j)) - x4(j) ) / D
        end do
      end do
      do i = 1, nr
        do j = 1, mr
          node = 2*i-1 + nn1row + (j-1)*(2*nn1row+2*nr)
          node1 = 2*i-1 + (j-1)*(2*nn1row+2*nr)
          node3 = 2*i-1 + j*(2*nn1row+2*nr)
          mesh%coor(node,:) = (mesh%coor(node1,:)+ &
                               mesh%coor(node3,:)+mesh%coor(node3+2,:))/3
          node = node + 1
          node3 = node3 + 2
          mesh%coor(node,:) = (mesh%coor(node1,:)+mesh%coor(node1+2,:)+&
                               mesh%coor(node3,:))/3
        end do
      end do
    end if

    deallocate ( x1, x2, x3, x4 )

!   points

    mesh%npoints = 4

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn1row
    mesh%points(3) = nn1row * nn1col + 2*nr*mr
    mesh%points(4) = nn1row * ( nn1col - 1 ) + 2*nr*mr + 1

!   curves

    mesh%ncurves = 4

!   allocate memory for curves

    allocate( mesh%curves(MAXCURVES) )
    allocate( mesh%curves(1)%nodes(nn1row) )
    allocate( mesh%curves(2)%nodes(nn1col) )
    allocate( mesh%curves(3)%nodes(nn1row) )
    allocate( mesh%curves(4)%nodes(nn1col) )
    allocate( mesh%curves(1)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(2)%topology(inpelmc,mr,2) )
    allocate( mesh%curves(3)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(4)%topology(inpelmc,mr,2) )

!   element type

    mesh%curves(:4)%element%elshape = 2
    mesh%curves(:4)%element%numnod = inpelmc
    mesh%curves(:4)%element%ndim = 2

!   other info of curves

    mesh%curves(:4)%ndim = 2
    mesh%curves(:4)%nblend = 0
    mesh%curves(1:4)%nnodes = [ nn1row, nn1col, nn1row, nn1col ]
    mesh%curves(1:4)%nelem = [ nr, mr, nr, mr ]
    mesh%curves(1:4)%n = mesh%curves(1:4)%nelem
    do i = 1, 4
      allocate(mesh%curves(i)%element_blend(0))
      allocate(mesh%curves(i)%nnodes_blend(2) )
      mesh%curves(i)%nnodes_blend = [0,mesh%curves(i)%nnodes]
    end do

!   nodes of curves

!   curve 1
    mesh%curves(1)%nodes = [ (i, i=1,nn1row) ]
!   curve 2
    mesh%curves(2)%nodes(1:nn1col:2) = [ (j*nn1row+(j-1)*nr, j=1,nn1col,2) ]
    mesh%curves(2)%nodes(2:nn1col:2) = [ (j*nn1row+j*nr, j=2,nn1col,2) ]
!   curve 3
    mesh%curves(3)%nodes = [ ((nn1col-1)*nn1row+2*nr*mr+i, i=nn1row,1,-1) ]
!   curve 4
    mesh%curves(4)%nodes(1:nn1col:2) = &
                            [ ((j-1)*(nn1row+nr)+1, j=nn1col,1,-2) ]
    mesh%curves(4)%nodes(2:nn1col:2) = &
                            [ ((j-1)*nn1row+j*nr+1, j=nn1col-1,2,-2) ]

!   topology of elements on curves

    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        mesh%curves(curve)%topology(:,elem,1) = [ (2*(elem-1)+i,i=1,3) ]
        mesh%curves(curve)%topology(:,elem,2) = &
           mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
      end do
    end do

  end subroutine rectangle2d_tria7


! Generate mesh on region [0,1]x[0,1] using quadrilateral elements with 4 nodes.

  subroutine rectangle2d_quad4 ( mesh, n, m, r, s, ratio, factor )

    type(mesh_t), intent(out) :: mesh

!   number of elements along the sides and refinement factors
!   The ACTUAL number of elements is n*r x m*s.
    integer, intent(in) :: n, m, r, s

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, dimension(:), intent(in) :: ratio
    real(dp), dimension(:), intent(in) :: factor

!   Generate a mesh on region  [0,1]x[0,1] using quads with 4 nodes.
!   The numbering for both nodal points and elements is row by row.
!   For example a 2x2 mesh is numbered as follows
!
!   Nodes:
!
!       7 ------- 8 ------- 9
!       |         |         |
!       |         |         |
!       |         |         |
!       4 ------- 5 ------- 6
!       |         |         |
!       |         |         |
!       |         |         |
!       1 ------- 2 ------- 3
!
!   Elements:
!
!       x ------- x ------- x
!       |         |         |
!       |    3    |    4    |
!       |         |         |
!       x ------- x ------- x
!       |         |         |
!       |    1    |    2    |
!       |         |         |
!       x ------- x ------- x
!
!  Also generated are four points and four curves:
!
!               <-----
!                  C3
!      P4 --------------------- P3
!       |                       |           Note:   C1=P1-P2
!       |                       |                   C2=P2-P3
!  |    |                       |    ^              C3=P3-P4
!  | C4 |                       | C2 |              C4=P4-P1
! \|/   |                       |    |
!       |                       |
!       |                       |
!      P1 --------------------- P2
!                  C1
!                ----->
!

    integer, parameter :: inpelm = 4 ! number of nodal points per element
    integer, parameter :: inpelmc = 2 ! number of nodal points per element on
                                      ! curves

    integer :: elnodes(inpelm), i, j, e, k
    integer :: nn1, nn2, nn1row, nn1col, node, elem, curve, nr, mr
    real(dp) :: deltax, deltay, D
    real(dp) :: deltax1, deltax2, deltax3, deltax4
    real(dp), allocatable, dimension(:) :: x1, x2, x3, x4

    allocate ( x1(1:n*r+1), x2(1:m*s+1), x3(1:n*r+1), x4(1:m*s+1) )


    nr = n*r
    mr = m*s

    nn1row = nr+1
    nn1col = mr+1

!   internal mesh

    mesh%ndim   = 2
    mesh%nelem  = nr*mr
    mesh%nnodes = nn1row*nn1col
    mesh%nelgrp = 1

!   allocate memory for internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%nnodes_blend(2) ); mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology(1) )
    allocate( mesh%topology(1)%a(inpelm,mesh%nelem) )
    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   element type

    mesh%element(1)%elshape = 5
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 2

!   topology

    mesh%grpnumel(1) = mesh%nelem

    do i = 1, nr
      do j = 1, mr

        nn1 = (j-1)*nn1row + i-1 ! number of nodes before element (i,j)
        nn2 = nn1 + nn1row   ! number of nodes before element (i,j) + 1 row

        e = i + (j-1)*nr   ! element number

        elnodes = [ (nn1+i,i=1,2), (nn2+i,i=2,1,-1) ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

      end do
    end do

!   coordinates

    if ( all(ratio == 0) ) then
!     equidistant
      deltax = 1.0_dp / nr
      deltay = 1.0_dp / mr
      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*nn1row
          mesh%coor(node,1) = (i-1) * deltax
          mesh%coor(node,2) = (j-1) * deltay
        end do
      end do
    else
!     non-equidistant
      call distribute_elements ( n, x1(::r), ratio(1), factor(1) )
      call distribute_elements ( m, x2(::s), ratio(2), factor(2) )
      call distribute_elements ( n, x3(::r), ratio(3), factor(3) )
      call distribute_elements ( m, x4(::s), ratio(4), factor(4) )
!     refine
      do elem = 1, n
        k = ( elem - 1 ) * r  ! nodes before element elem
        deltax1 = (x1(k+r+1)-x1(k+1))/r
        deltax3 = (x3(k+r+1)-x3(k+1))/r
        do node = k+2, k+r
          x1(node) = x1(node-1) + deltax1
          x3(node) = x3(node-1) + deltax3
        end do
      end do
      do elem = 1, m
        k = ( elem - 1 ) * s  ! nodes before element elem
        deltax2 = (x2(k+s+1)-x2(k+1))/s
        deltax4 = (x4(k+s+1)-x4(k+1))/s
        do node = k+2, k+s
          x2(node) = x2(node-1) + deltax2
          x4(node) = x4(node-1) + deltax4
        end do
      end do
      x3 = 1._dp - x3(nn1row:1:-1)
      x4 = 1._dp - x4(nn1col:1:-1)
!     create straight lines in reference square [0,1]x[0,1]
      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*nn1row
          D = ( x1(i) - x3(i) ) * ( x4(j) - x2(j) ) - 1
          mesh%coor(node,1) = ( x4(j) * (x1(i) - x3(i)) - x1(i) ) / D
          mesh%coor(node,2) = ( x1(i) * (x4(j) - x2(j)) - x4(j) ) / D
        end do
      end do
    end if

    deallocate ( x1, x2, x3, x4 )

!   points

    mesh%npoints = 4

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn1row
    mesh%points(3) = nn1row * nn1col
    mesh%points(4) = nn1row * ( nn1col - 1 ) + 1

!   curves

    mesh%ncurves = 4

!   allocate memory for curves

    allocate( mesh%curves(MAXCURVES) )
    allocate( mesh%curves(1)%nodes(nn1row) )
    allocate( mesh%curves(2)%nodes(nn1col) )
    allocate( mesh%curves(3)%nodes(nn1row) )
    allocate( mesh%curves(4)%nodes(nn1col) )
    allocate( mesh%curves(1)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(2)%topology(inpelmc,mr,2) )
    allocate( mesh%curves(3)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(4)%topology(inpelmc,mr,2) )

!   element type

    mesh%curves(:4)%element%elshape = 1
    mesh%curves(:4)%element%numnod = inpelmc
    mesh%curves(:4)%element%ndim = 2

!   other info of curves

    mesh%curves(:4)%ndim = 2
    mesh%curves(:4)%nblend = 0
    mesh%curves(1:4)%nnodes = [ nn1row, nn1col, nn1row, nn1col ]
    mesh%curves(1:4)%nelem = [ nr, mr, nr, mr ]
    mesh%curves(1:4)%n = mesh%curves(1:4)%nelem
    do i = 1, 4
      allocate(mesh%curves(i)%element_blend(0))
      allocate(mesh%curves(i)%nnodes_blend(2) )
      mesh%curves(i)%nnodes_blend = [0,mesh%curves(i)%nnodes]
    end do

!   nodes of curves

    mesh%curves(1)%nodes = [ (i, i=1,nn1row) ]                      ! curve 1
    mesh%curves(2)%nodes = [ (j*nn1row, j=1,nn1col) ]               ! curve 2
    mesh%curves(3)%nodes = [ ((nn1col-1)*nn1row+i, i=nn1row,1,-1) ] ! curve 3
    mesh%curves(4)%nodes = [ ((j-1)*nn1row+1, j=nn1col,1,-1) ]      ! curve 4

!   topology of elements on curves

    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        mesh%curves(curve)%topology(:,elem,1) = [ (elem-1+i,i=1,2) ]
        mesh%curves(curve)%topology(:,elem,2) = &
           mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
      end do
    end do

  end subroutine rectangle2d_quad4


! Generate mesh on region [0,1]x[0,1] using quadrilateral elements with 5 nodes.

  subroutine rectangle2d_quad5 ( mesh, n, m, r, s, ratio, factor )

    type(mesh_t), intent(out) :: mesh

!   number of elements along the sides and refinement factors
!   The ACTUAL number of elements is n*r x m*s.
    integer, intent(in) :: n, m, r, s

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, dimension(:), intent(in) :: ratio
    real(dp), dimension(:), intent(in) :: factor

!   Generate a mesh on region  [0,1]x[0,1] using quads with 5 nodes.
!   The numbering for both nodal points and elements is row by row.
!   For example a 2x2 mesh is numbered as follows
!
!   Nodes:
!
!       11------- 12------- 13
!       |         |         |
!       |    9    |   10    |
!       |         |         |
!       6 ------- 7 ------- 8
!       |         |         |
!       |    4    |    5    |
!       |         |         |
!       1 ------- 2 ------- 3
!
!   Elements:
!
!       x ------- x ------- x
!       |         |         |
!       |    3    |    4    |
!       |         |         |
!       x ------- x ------- x
!       |         |         |
!       |    1    |    2    |
!       |         |         |
!       x ------- x ------- x
!
!  Also generated are four points and four curves:
!
!               <-----
!                  C3
!      P4 --------------------- P3
!       |                       |           Note:   C1=P1-P2
!       |                       |                   C2=P2-P3
!  |    |                       |    ^              C3=P3-P4
!  | C4 |                       | C2 |              C4=P4-P1
! \|/   |                       |    |
!       |                       |
!       |                       |
!      P1 --------------------- P2
!                  C1
!                ----->
!

    integer, parameter :: inpelm = 5 ! number of nodal points per element
    integer, parameter :: inpelmc = 2 ! number of nodal points per element on
                                      ! curves

    integer :: elnodes(inpelm), i, j, e, k, nr, mr
    integer :: nn1, nn2, nn3, nn1row, nn1col, node, elem, curve, node1, node4
    real(dp) :: deltax, deltay, D
    real(dp) :: deltax1, deltax2, deltax3, deltax4
    real(dp), allocatable, dimension(:) :: x1, x2, x3, x4

    allocate ( x1(1:n*r+1), x2(1:m*s+1), x3(1:n*r+1), x4(1:m*s+1) )


    nr = n*r
    mr = m*s

    nn1row = n*r+1
    nn1col = m*s+1

!   internal mesh

    mesh%ndim   = 2
    mesh%nelem  = nr*mr
    mesh%nnodes = nn1row*nn1col + nr*mr
    mesh%nelgrp = 1

!   allocate memory for internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%nnodes_blend(2) ); mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology(1) )
    allocate( mesh%topology(1)%a(inpelm,mesh%nelem) )
    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   element type

    mesh%element(1)%elshape = 9
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 2

!   topology

    mesh%grpnumel(1) = mesh%nelem

    do i = 1, nr
      do j = 1, mr

        nn1 = (j-1)*(nn1row+nr) + i-1 ! number of nodes before element (i,j)
        nn2 = nn1 + nn1row   ! number of nodes before element (i,j) + 1 row
        nn3 = nn2 + nr        ! number of nodes before element (i,j) + 2 rows

        e = i + (j-1)*nr   ! element number

        elnodes = [ (nn1+i,i=1,2), (nn3+i,i=2,1,-1), nn2+1 ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

      end do
    end do

!   coordinates


    if ( all(ratio == 0) ) then
!     equidistant
      deltax = 1.0_dp / nr
      deltay = 1.0_dp / mr
      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*(nn1row+nr)
          mesh%coor(node,1) = (i-1) * deltax
          mesh%coor(node,2) = (j-1) * deltay
        end do
      end do
      do i = 1, nr
        do j = 1, mr
          node = i + (j-1)*(nn1row+nr) + nn1row
          mesh%coor(node,1) = (i-1) * deltax + deltax/2
          mesh%coor(node,2) = (j-1) * deltay + deltay/2
        end do
      end do
    else
!     non-equidistant
      call distribute_elements ( n, x1(::r), ratio(1), factor(1) )
      call distribute_elements ( m, x2(::s), ratio(2), factor(2) )
      call distribute_elements ( n, x3(::r), ratio(3), factor(3) )
      call distribute_elements ( m, x4(::s), ratio(4), factor(4) )
!     refine
      do elem = 1, n
        k = ( elem - 1 ) * r  ! nodes before element elem
        deltax1 = (x1(k+r+1)-x1(k+1))/r
        deltax3 = (x3(k+r+1)-x3(k+1))/r
        do node = k+2, k+r
          x1(node) = x1(node-1) + deltax1
          x3(node) = x3(node-1) + deltax3
        end do
      end do
      do elem = 1, m
        k = ( elem - 1 ) * s  ! nodes before element elem
        deltax2 = (x2(k+s+1)-x2(k+1))/s
        deltax4 = (x4(k+s+1)-x4(k+1))/s
        do node = k+2, k+s
          x2(node) = x2(node-1) + deltax2
          x4(node) = x4(node-1) + deltax4
        end do
      end do
      x3 = 1._dp - x3(nn1row:1:-1)
      x4 = 1._dp - x4(nn1col:1:-1)
!     create straight lines in reference square [0,1]x[0,1]
      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*(nn1row+nr)
          D = ( x1(i) - x3(i) ) * ( x4(j) - x2(j) ) - 1
          mesh%coor(node,1) = ( x4(j) * (x1(i) - x3(i)) - x1(i) ) / D
          mesh%coor(node,2) = ( x1(i) * (x4(j) - x2(j)) - x4(j) ) / D
        end do
      end do
      do i = 1, nr
        do j = 1, mr
          node = i + (j-1)*(nn1row+nr) + nn1row
          node1 = i + (j-1)*(nn1row+nr)
          node4 = i + j*(nn1row+nr)
          mesh%coor(node,:) = (mesh%coor(node1,:)+mesh%coor(node1+1,:)+ &
                               mesh%coor(node4,:)+mesh%coor(node4+1,:))/4
        end do
      end do
    end if

    deallocate ( x1, x2, x3, x4 )

!   points

    mesh%npoints = 4

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn1row
    mesh%points(3) = nn1row * nn1col + nr*mr
    mesh%points(4) = nn1row * ( nn1col - 1 ) + nr*mr + 1

!   curves

    mesh%ncurves = 4

!   allocate memory for curves

    allocate( mesh%curves(MAXCURVES) )
    allocate( mesh%curves(1)%nodes(nn1row) )
    allocate( mesh%curves(2)%nodes(nn1col) )
    allocate( mesh%curves(3)%nodes(nn1row) )
    allocate( mesh%curves(4)%nodes(nn1col) )
    allocate( mesh%curves(1)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(2)%topology(inpelmc,mr,2) )
    allocate( mesh%curves(3)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(4)%topology(inpelmc,mr,2) )

!   element type

    mesh%curves(:4)%element%elshape = 1
    mesh%curves(:4)%element%numnod = inpelmc
    mesh%curves(:4)%element%ndim = 2

!   other info of curves

    mesh%curves(:4)%ndim = 2
    mesh%curves(:4)%nblend = 0
    mesh%curves(1:4)%nnodes = [ nn1row, nn1col, nn1row, nn1col ]
    mesh%curves(1:4)%nelem = [ nr, mr, nr, mr ]
    mesh%curves(1:4)%n = mesh%curves(1:4)%nelem
    do i = 1, 4
      allocate(mesh%curves(i)%element_blend(0))
      allocate(mesh%curves(i)%nnodes_blend(2) )
      mesh%curves(i)%nnodes_blend = [0,mesh%curves(i)%nnodes]
    end do

!   nodes of curves

    mesh%curves(1)%nodes = [ (i, i=1,nn1row) ]                      ! curve 1
    mesh%curves(2)%nodes = [ (j*nn1row+(j-1)*nr, j=1,nn1col) ]       ! curve 2
!                                                                       curve 3
    mesh%curves(3)%nodes = [ ((nn1col-1)*nn1row+nr*mr+i, i=nn1row,1,-1) ]
    mesh%curves(4)%nodes = [ ((j-1)*(nn1row+nr)+1, j=nn1col,1,-1) ]  ! curve 4

!   topology of elements on curves

    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        mesh%curves(curve)%topology(:,elem,1) = [ (elem-1+i,i=1,2) ]
        mesh%curves(curve)%topology(:,elem,2) = &
           mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
      end do
    end do

  end subroutine rectangle2d_quad5


! Generate mesh on region [0,1]x[0,1] using quadrilateral elements with 9 nodes.

  subroutine rectangle2d_quad9 ( mesh, n, m, r, s, ratio, factor)

    type(mesh_t), intent(out) :: mesh

!   number of elements along the sides and refinement factors
!   The ACTUAL number of elements is n*r x m*s.
    integer, intent(in) :: n, m, r, s

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, dimension(:), intent(in) :: ratio
    real(dp), dimension(:), intent(in) :: factor

!   Generate a mesh on region  [0,1]x[0,1] using quads with 9 nodes.
!   The numbering for both nodal points and elements is row by row.
!   For example a 2x2 mesh is numbered as follows
!
!   Nodes:
!
!      21 -- 22 -- 23 -- 24 -- 25
!       |           |           |
!      16    17    18    19    20
!       |           |           |
!      11 -- 12 -- 13 -- 14 -- 15
!       |           |           |
!       6     7     8     9    10
!       |           |           |
!       1 --  2 --  3 --  4 --  5
!
!   Elements:
!
!       x --  x --  x --  x --  x
!       |           |           |
!       x     3     x     4     x
!       |           |           |
!       x --  x --  x --  x --  x
!       |           |           |
!       x     1     x     2     x
!       |           |           |
!       x --  x --  x --  x --  x
!
!  Also generated are four points and four curves:
!
!               <-----
!                  C3
!      P4 --------------------- P3
!       |                       |           Note:   C1=P1-P2
!       |                       |                   C2=P2-P3
!  |    |                       |    ^              C3=P3-P4
!  | C4 |                       | C2 |              C4=P4-P1
! \|/   |                       |    |
!       |                       |
!       |                       |
!      P1 --------------------- P2
!                  C1
!                ----->
!

    integer, parameter :: inpelm = 9 ! number of nodal points per element
    integer, parameter :: inpelmc = 3 ! number of nodal points per element on
                                      ! curves

    integer :: elnodes(inpelm), i, j, e, k
    integer :: nn1, nn2, nn3, nn1row, nn1col, node, elem, curve, nr, mr
    real(dp) :: deltax, deltay, D
    real(dp) :: deltax1, deltax2, deltax3, deltax4
    real(dp), allocatable, dimension(:) :: x1, x2, x3, x4

    allocate ( x1(1:2*n*r+1), x2(1:2*m*s+1), x3(1:2*n*r+1), x4(1:2*m*s+1) )


    nr = n*r
    mr = m*s

    nn1row = 2*nr+1
    nn1col = 2*mr+1

!   internal mesh

    mesh%ndim   = 2
    mesh%nelem  = nr*mr
    mesh%nnodes = nn1row*nn1col
    mesh%nelgrp = 1

!   allocate memory for internal mesh

    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%nnodes_blend(2) ); mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology(1) )
    allocate( mesh%topology(1)%a(inpelm,mesh%nelem) )
    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   element type

    mesh%element(1)%elshape = 6
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 2

!   topology

    mesh%grpnumel(1) = mesh%nelem

    do i = 1, nr
      do j = 1, mr

        nn1 = 2*(j-1)*nn1row + 2*(i-1) ! number of nodes before element (i,j)
        nn2 = nn1 + nn1row   ! number of nodes before element (i,j) + 1 row
        nn3 = nn2 + nn1row   ! number of nodes before element (i,j) + 2 rows

        e = i + (j-1)*nr   ! element number

        elnodes = [ (nn1+i,i=1,3), nn2+3, (nn3+i,i=3,1,-1), (nn2+i,i=1,2) ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

      end do
    end do

!   coordinates

    if ( all(ratio == 0) ) then
!     equidistant
      deltax = 0.5_dp / nr
      deltay = 0.5_dp / mr
      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*nn1row
          mesh%coor(node,1) = (i-1) * deltax
          mesh%coor(node,2) = (j-1) * deltay
        end do
      end do
    else
!     non-equidistant
      call distribute_elements ( n, x1(1:nn1row:2*r), ratio(1), factor(1) )
      call distribute_elements ( m, x2(1:nn1col:2*s), ratio(2), factor(2) )
      call distribute_elements ( n, x3(1:nn1row:2*r), ratio(3), factor(3) )
      call distribute_elements ( m, x4(1:nn1col:2*s), ratio(4), factor(4) )
!     mid-side nodes and refine
      do elem = 1, n
        k = ( elem - 1 ) * 2*r  ! nodes before element elem
        deltax1 = (x1(k+2*r+1)-x1(k+1))/2/r
        deltax3 = (x3(k+2*r+1)-x3(k+1))/2/r
        do node = k+2, k+2*r
          x1(node) = x1(node-1) + deltax1
          x3(node) = x3(node-1) + deltax3
        end do
      end do
      do elem = 1, m
        k = ( elem - 1 ) * 2*s  ! nodes before element elem
        deltax2 = (x2(k+2*s+1)-x2(k+1))/2/s
        deltax4 = (x4(k+2*s+1)-x4(k+1))/2/s
        do node = k+2, k+2*s
          x2(node) = x2(node-1) + deltax2
          x4(node) = x4(node-1) + deltax4
        end do
      end do
      x3 = 1._dp - x3(nn1row:1:-1)
      x4 = 1._dp - x4(nn1col:1:-1)
!     create straight lines in reference square [0,1]x[0,1]
      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*nn1row
          D = ( x1(i) - x3(i) ) * ( x4(j) - x2(j) ) - 1
          mesh%coor(node,1) = ( x4(j) * (x1(i) - x3(i)) - x1(i) ) / D
          mesh%coor(node,2) = ( x1(i) * (x4(j) - x2(j)) - x4(j) ) / D
        end do
      end do
    end if

    deallocate ( x1, x2, x3, x4 )

!   points

    mesh%npoints = 4

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn1row
    mesh%points(3) = nn1row * nn1col
    mesh%points(4) = nn1row * ( nn1col - 1 ) + 1

!   curves

    mesh%ncurves = 4

!   allocate memory for curves

    allocate( mesh%curves(MAXCURVES) )
    allocate( mesh%curves(1)%nodes(nn1row) )
    allocate( mesh%curves(2)%nodes(nn1col) )
    allocate( mesh%curves(3)%nodes(nn1row) )
    allocate( mesh%curves(4)%nodes(nn1col) )
    allocate( mesh%curves(1)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(2)%topology(inpelmc,mr,2) )
    allocate( mesh%curves(3)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(4)%topology(inpelmc,mr,2) )

!   element type

    mesh%curves(:4)%element%elshape = 2
    mesh%curves(:4)%element%numnod = inpelmc
    mesh%curves(:4)%element%ndim = 2

!   other info of curves

    mesh%curves(:4)%ndim = 2
    mesh%curves(:4)%nblend = 0
    mesh%curves(1:4)%nnodes = [ nn1row, nn1col, nn1row, nn1col ]
    mesh%curves(1:4)%nelem = [ nr, mr, nr, mr ]
    mesh%curves(1:4)%n = mesh%curves(1:4)%nelem
    do i = 1, 4
      allocate(mesh%curves(i)%element_blend(0))
      allocate(mesh%curves(i)%nnodes_blend(2) )
      mesh%curves(i)%nnodes_blend = [0,mesh%curves(i)%nnodes]
    end do

!   nodes of curves

    mesh%curves(1)%nodes = [ (i, i=1,nn1row) ]                      ! curve 1
    mesh%curves(2)%nodes = [ (j*nn1row, j=1,nn1col) ]               ! curve 2
    mesh%curves(3)%nodes = [ ((nn1col-1)*nn1row+i, i=nn1row,1,-1) ] ! curve 3
    mesh%curves(4)%nodes = [ ((j-1)*nn1row+1, j=nn1col,1,-1) ]      ! curve 4

!   topology of elements on curves

    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        mesh%curves(curve)%topology(:,elem,1) = [ (2*(elem-1)+i,i=1,3) ]
        mesh%curves(curve)%topology(:,elem,2) = &
           mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
      end do
    end do

  end subroutine rectangle2d_quad9


! Generate mesh on region [0,1]x[0,1] using quadrilateral high-order elements
! with (p+1)^2 nodes, where p is the polynomial order which will be equal
! in both directions

  subroutine rectangle2d_quad_ho ( mesh, n, m, r, s, ratio, factor, p, l )

    use gauss_spectral_m

    type(mesh_t), intent(out) :: mesh

!   number of elements along the sides and refinement factors
!   The ACTUAL number of elements is n*r x m*s.
    integer, intent(in) :: n, m, r, s

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, dimension(:), intent(in) :: ratio
    real(dp), dimension(:), intent(in) :: factor

!   the polynomial order
    integer, intent(in) :: p

!   layout: l=0: equal distribution l=1: Gauss-Lobatto
    integer, intent(in) :: l

!   Generate a mesh on region  [0,1]x[0,1] using quads with (p+1)^2 nodes.
!   The numbering for both nodal points and elements is row by row.
!   For example a 2x2 mesh with p=3 is numbered as follows
!
!   Nodes:
!      43 -- 44 -- 45 -- 46 -- 47 --48 --49
!       |                 |               |
!      36    37    38    39    40   41   42
!       |                 |               |
!      29    30    31    32    33   34   35
!       |                 |               |
!      22 -- 23 -- 24 -- 25 -- 26 --27 --28
!       |                 |               |
!      15    16    17    18    19   20   21
!       |                 |               |
!       8     9    10    11    12   13   14
!       |                 |               |
!       1 --  2 --  3 --  4 --  5 -- 6 -- 7
!
!   Elements:
!
!       x --  x --  x --  x --  x -- x -- x
!       |                 |               |
!       x                 x               x
!       |        3        |       4       |
!       x                 x               x
!       |                 |               |
!       x --  x --  x --  x --  x -- x -- x
!       |                 |               |
!       x                 x               x
!       |        1        |       2       |
!       x                 x               x
!       |                 |               |
!       x --  x --  x --  x --  x -- x -- x
!
!  Also generated are four points and four curves:
!
!               <-----
!                  C3
!      P4 --------------------- P3
!       |                       |           Note:   C1=P1-P2
!       |                       |                   C2=P2-P3
!  |    |                       |    ^              C3=P3-P4
!  | C4 |                       | C2 |              C4=P4-P1
! \|/   |                       |    |
!       |                       |
!       |                       |
!      P1 --------------------- P2
!                  C1
!                ----->
!

!   to get the distribution on a line with polynomial order p
    real(dp), dimension(:), allocatable  :: x

    integer :: i, j, e, k, kl, inpelm, inpelmc, k2
    integer :: nn1, nn1row, nn1col, node, elem, curve, nr, mr
    integer :: elemx, elemy
    integer, dimension(:), allocatable :: elnodes
    real(dp), dimension(:), allocatable :: xdistribution, ydistribution
    real(dp) :: deltax, deltay, D
    real(dp) :: deltax1, deltax2, deltax3, deltax4
    real(dp), allocatable, dimension(:) :: x1, x2, x3, x4

    allocate ( x1(1:p*n*r+1), x2(1:p*m*s+1), x3(1:p*n*r+1), x4(1:p*m*s+1) )


    nr = n*r
    mr = m*s
    inpelm = (p+1)**2 ! number of nodal points per element
    inpelmc = p+1     ! number of nodal points per element on curves
    allocate(x(inpelmc))
    allocate(elnodes(inpelm))

    nn1row = p*nr+1
    nn1col = p*mr+1

!   internal mesh

    mesh%ndim   = 2
    mesh%nelem  = nr*mr
    mesh%nnodes = nn1row*nn1col
    mesh%nelgrp = 1

!   allocate memory for internal mesh
    allocate( xdistribution (nn1row) )
    allocate( ydistribution (nn1col) )
    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%nnodes_blend(2) ); mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology(1) )
    allocate( mesh%topology(1)%a(inpelm,mesh%nelem) )
    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   element type

    mesh%element(1)%elshape = 102
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 2
    mesh%element(1)%p(1,1) = p
    mesh%element(1)%p(2,1) = p
    mesh%element(1)%p(1,2) = l
    mesh%element(1)%p(2,2) = l

!   topology

    mesh%grpnumel(1) = mesh%nelem

!   node numbers per element. First x-coordinate and then the y-coordinate
!   equal polynomial order in both directions.

    do i = 1, nr
      do j = 1, mr
        do kl = 1, p+1
!         number of nodes before element (i,j)
          nn1 = p*(i-1)+p*(nr*p+1)*(j-1)+(nr*p+1)*(kl-1)
          elnodes(1+(p+1)*(kl-1):(p+1)+(p+1)*(kl-1)) = [ (nn1+i,i=1,p+1) ]
        end do
        e = i + (j-1)*nr   ! element number
        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)
      end do
    end do

!   coordinates

    if ( l == 0 ) then

!     equidistant distribution on  [-1,1]

      x = -1 + 2 * [ ( i - 1, i=1,p+1 ) ] / real(p,kind=dp)

    else if ( l == 1 ) then

!     get the GLL distribution
!     with x the GLL-points (output) and p the polynomial order (input)

      call GLL_points ( x, p )

    end if

    if ( all(ratio == 0) ) then

!     equidistant mesh

!     the x values are between -1 and 1

      deltax = 1._dp / nr

      do elemx = 1, nr
        do i = 1, p
          node = (elemx-1)*p + i
          xdistribution(node) = ( (x(i)+1._dp)/2 + elemx-1 ) * deltax
         end do
      end do
      xdistribution(nn1row) = 1

      deltay = 1._dp / mr

      do elemy = 1, mr
         do i = 1 ,p
            node = (elemy-1)*p + i
            ydistribution(node)= ( (x(i)+1._dp)/2 + elemy-1 ) * deltay
         end do
      end do
      ydistribution(nn1col) = 1

      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*nn1row
          mesh%coor(node,1) = xdistribution(i)
          mesh%coor(node,2) = ydistribution(j)
        end do
      end do

    else

!     non-equidistant

      call distribute_elements ( n, x1(1:nn1row:p*r), ratio(1), factor(1) )
      call distribute_elements ( m, x2(1:nn1col:p*s), ratio(2), factor(2) )
      call distribute_elements ( n, x3(1:nn1row:p*r), ratio(3), factor(3) )
      call distribute_elements ( m, x4(1:nn1col:p*s), ratio(4), factor(4) )

!     mid-side nodes and refine

      do elem = 1, n ! x-distribution

        k = ( elem - 1 ) * p*r  ! nodes before element elem

        deltax1 = (x1(k+p*r+1)-x1(k+1))/r
        deltax3 = (x3(k+p*r+1)-x3(k+1))/r

        do i = 1, r
          k2 = k + (i - 1._dp) * p
          do node = k2+2, k2+p+1
            j=node-k2 ! which point number in reference element
            x1(node) = x1(k2+1) + deltax1 * ((x(j)+1._dp)/2._dp)
!           x(j) asks for the distribution of the j'th point,
!           with the coordinates adapted to coordinates between 0 and 1.
            x3(node) = x3(k2+1) + deltax3 * ((x(j)+1._dp)/2._dp)
          end do
        end do

      end do

      do elem = 1, m ! y-distribution

        k = ( elem - 1 ) * p*s  ! nodes before element elem

        deltax2 = (x2(k+p*s+1)-x2(k+1))/s
        deltax4 = (x4(k+p*s+1)-x4(k+1))/s

        do i = 1, s
          k2 = k + (i - 1._dp) * p
          do node = k2+2, k2+p+1
            j=node-k2
            x2(node) = x2(k2+1) + deltax2 * ((x(j)+1._dp)/2._dp)
            x4(node) = x4(k2+1) + deltax4 * ((x(j)+1._dp)/2._dp)
          end do
        end do

      end do

      x3 = 1._dp - x3(nn1row:1:-1)
      x4 = 1._dp - x4(nn1col:1:-1)

!     create straight lines in reference square [0,1]x[0,1]

      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*nn1row
          D = ( x1(i) - x3(i) ) * ( x4(j) - x2(j) ) - 1
          mesh%coor(node,1) = ( x4(j) * (x1(i) - x3(i)) - x1(i) ) / D
          mesh%coor(node,2) = ( x1(i) * (x4(j) - x2(j)) - x4(j) ) / D
        end do
      end do

    end if

    deallocate ( x1, x2, x3, x4 )

!   points

    mesh%npoints = 4

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn1row
    mesh%points(3) = nn1row * nn1col
    mesh%points(4) = nn1row * ( nn1col - 1 ) + 1

!   curves

    mesh%ncurves = 4

!   allocate memory for curves

    allocate( mesh%curves(MAXCURVES) )
    allocate( mesh%curves(1)%nodes(nn1row) )
    allocate( mesh%curves(2)%nodes(nn1col) )
    allocate( mesh%curves(3)%nodes(nn1row) )
    allocate( mesh%curves(4)%nodes(nn1col) )
    allocate( mesh%curves(1)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(2)%topology(inpelmc,mr,2) )
    allocate( mesh%curves(3)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(4)%topology(inpelmc,mr,2) )

!   element type

    mesh%curves(:4)%element%elshape = 101
    mesh%curves(:4)%element%numnod = inpelmc
    mesh%curves(:4)%element%p(1,1) = p
    mesh%curves(:4)%element%p(2,1) = l
    mesh%curves(:4)%element%ndim = 2

!   other info of curves

    mesh%curves(:4)%ndim = 2
    mesh%curves(:4)%nblend = 0
    mesh%curves(1:4)%nnodes = [ nn1row, nn1col, nn1row, nn1col ]
    mesh%curves(1:4)%nelem = [ nr, mr, nr, mr ]
    mesh%curves(1:4)%n = mesh%curves(1:4)%nelem
    do i = 1, 4
      allocate(mesh%curves(i)%element_blend(0))
      allocate(mesh%curves(i)%nnodes_blend(2) )
      mesh%curves(i)%nnodes_blend = [0,mesh%curves(i)%nnodes]
    end do

!   nodes of curves

    mesh%curves(1)%nodes = [ (i, i=1,nn1row) ]                      ! curve 1
    mesh%curves(2)%nodes = [ (j*nn1row, j=1,nn1col) ]               ! curve 2
    mesh%curves(3)%nodes = [ ((nn1col-1)*nn1row+i, i=nn1row,1,-1) ] ! curve 3
    mesh%curves(4)%nodes = [ ((j-1)*nn1row+1, j=nn1col,1,-1) ]      ! curve 4

!   topology of elements on curves
    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        mesh%curves(curve)%topology(:,elem,1) = [ (p*(elem-1)+i,i=1,p+1) ]
        mesh%curves(curve)%topology(:,elem,2) = &
           mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
      end do
    end do

  end subroutine rectangle2d_quad_ho


! Generate mesh on region [0,1]x[0,1] using triangular high-order elements
! with (p+1)*(p+2)/2 nodes, where p is the polynomial order which will be equal
! in both directions

  subroutine rectangle2d_triangle_ho ( mesh, n, m, r, s, ratio, factor, p )

    type(mesh_t), intent(out) :: mesh

!   number of elements along the sides and refinement factors
!   The ACTUAL number of elements is n*r x m*s.
    integer, intent(in) :: n, m, r, s

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, dimension(:), intent(in) :: ratio
    real(dp), dimension(:), intent(in) :: factor

!   the polynomial order
    integer, intent(in) :: p

!   Generate a mesh on region  [0,1]x[0,1] using triangles with (p+1)*(p+2)/2
!   nodes. The numbering for both nodal points and elements is row by row.
!   For example a 2x2 mesh with p=3 is numbered as follows
!
!   Nodes:
!      43 -- 44 -- 45 -- 46 -- 47 --48 --49
!       |              /  |            /  |
!      36    37    38    39    40   41   42
!       |        /        |       /       |
!      29    30    31    32    33   34   35
!       |  /              |  /            |
!      22 -- 23 -- 24 -- 25 -- 26 --27 --28
!       |             /   |            /  |
!      15    16    17    18    19   20   21
!       |        /        |       /       |
!       8     9    10    11    12   13   14
!       |  /              |  /            |
!       1 --  2 --  3 --  4 --  5 -- 6 -- 7
!
!   Elements:
!
!       x --  x --  x --  x --  x -- x -- x
!       |              /  |             / |
!       x    5      /     x    7     /    x
!       |        /        |       /       |
!       x    /      6     x    /     8    x
!       | /               | /             |
!       x --  x --  x --  x --  x -- x -- x
!       |              /  |             / |
!       x    1      /     x    3     /    x
!       |        /        |       /       |
!       x    /       2    x    /     4    x
!       | /               | /             |
!       x --  x --  x --  x --  x -- x -- x
!
!  Also generated are four points and four curves:
!
!               <-----
!                  C3
!      P4 --------------------- P3
!       |                       |           Note:   C1=P1-P2
!       |                       |                   C2=P2-P3
!  |    |                       |    ^              C3=P3-P4
!  | C4 |                       | C2 |              C4=P4-P1
! \|/   |                       |    |
!       |                       |
!       |                       |
!      P1 --------------------- P2
!                  C1
!                ----->
!

!   to get the distribution on a line with polynomial order p
    real(dp), dimension(:), allocatable  :: x

    integer :: i, j, e, k, kl, inpelm, inpelmc, k2, ii, jj, mj
    integer :: nn1, nn1row, nn1col, node, elem, curve, nr, mr
    integer :: elemx, elemy
    integer, dimension(:), allocatable :: elnodes
    integer, dimension(:,:), allocatable :: mat
    real(dp), dimension(:), allocatable :: xdistribution, ydistribution
    real(dp) :: deltax, deltay, D
    real(dp) :: deltax1, deltax2, deltax3, deltax4
    real(dp), allocatable, dimension(:) :: x1, x2, x3, x4

    allocate ( x1(1:p*n*r+1), x2(1:p*m*s+1), x3(1:p*n*r+1), x4(1:p*m*s+1) )
    allocate ( mat(p+1,p+1) )


    nr = n*r
    mr = m*s
    inpelm = (p+1)*(p+2)/2 ! number of nodal points per element
    inpelmc = p+1          ! number of nodal points per element on curves
    allocate(x(inpelmc))
    allocate(elnodes(inpelm))

    nn1row = p*nr+1
    nn1col = p*mr+1

!   internal mesh

    mesh%ndim   = 2
    mesh%nelem  = nr*mr*2
    mesh%nnodes = nn1row*nn1col
    mesh%nelgrp = 1

!   allocate memory for internal mesh
    allocate( xdistribution (nn1row) )
    allocate( ydistribution (nn1col) )
    allocate( mesh%element ( mesh%nelgrp ) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%nnodes_blend(2) ); mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel ( mesh%nelgrp ) )
    allocate( mesh%topology(1) )
    allocate( mesh%topology(1)%a(inpelm,mesh%nelem) )
    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   element type

    mesh%element(1)%elshape = 104
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 2
    mesh%element(1)%p(1,1) = p
    mesh%element(1)%p(2,1) = p
    mesh%element(1)%p(1,2) = 0
    mesh%element(1)%p(2,2) = 0

!   topology

    mesh%grpnumel(1) = mesh%nelem

!   node numbers per element.

    do i = 1, nr
      do j = 1, mr

!       fill mat with nodal numbers of the two elements together

        do kl = 1, p+1
!         number of nodes before element (i,j) at local row kl
          nn1 = p*(i-1)+nn1row*(p*(j-1)+kl-1)
          mat(kl,:) = [ (nn1+i,i=1,p+1) ]
        end do

        e = 2*i-1 + (j-1)*2*nr   ! element number

        do jj = 1, p+1
          do ii = 1, p+2-jj
            mj= (jj-1)*(2*p+4-jj)/2
            elnodes(ii+mj) = mat(p+2-ii,jj)
          end do
        end do

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

        e = e + 1

        do jj = 1, p+1
          do ii = 1, p+2-jj
            mj = (jj-1)*(2*p+4-jj)/2
            elnodes(ii+mj) = mat(ii,p+2-jj)
          end do
        end do

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

      end do
    end do

!   coordinates

!   equidistant distribution on [0,1]

    x = [ ( i - 1, i=1,p+1 ) ] / real(p,kind=dp)

    if ( all(ratio == 0) ) then

!     equidistant mesh

      deltax = 1._dp / nr

      do elemx = 1, nr
        do i = 1, p
          node = (elemx-1)*p + i
          xdistribution(node) = ( x(i) + elemx-1 ) * deltax
         end do
      end do
      xdistribution(nn1row) = 1

      deltay = 1._dp / mr

      do elemy = 1, mr
         do i = 1 ,p
            node = (elemy-1)*p + i
            ydistribution(node)= ( x(i) + elemy-1 ) * deltay
         end do
      end do
      ydistribution(nn1col) = 1

      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*nn1row
          mesh%coor(node,1) = xdistribution(i)
          mesh%coor(node,2) = ydistribution(j)
        end do
      end do

    else

!     non-equidistant

      call distribute_elements ( n, x1(1:nn1row:p*r), ratio(1), factor(1) )
      call distribute_elements ( m, x2(1:nn1col:p*s), ratio(2), factor(2) )
      call distribute_elements ( n, x3(1:nn1row:p*r), ratio(3), factor(3) )
      call distribute_elements ( m, x4(1:nn1col:p*s), ratio(4), factor(4) )

!     mid-side nodes and refine

      do elem = 1, n ! x-distribution

        k = ( elem - 1 ) * p*r  ! nodes before element elem

        deltax1 = (x1(k+p*r+1)-x1(k+1))/r
        deltax3 = (x3(k+p*r+1)-x3(k+1))/r

        do i = 1, r
          k2 = k + (i - 1._dp) * p
          do node = k2+2, k2+p+1
            j=node-k2 ! which point number in reference element
            x1(node) = x1(k2+1) + deltax1 * x(j)
            x3(node) = x3(k2+1) + deltax3 * x(j)
          end do
        end do

      end do

      do elem = 1, m ! y-distribution

        k = ( elem - 1 ) * p*s  ! nodes before element elem

        deltax2 = (x2(k+p*s+1)-x2(k+1))/s
        deltax4 = (x4(k+p*s+1)-x4(k+1))/s

        do i = 1, s
          k2 = k + (i - 1._dp) * p
          do node = k2+2, k2+p+1
            j=node-k2
            x2(node) = x2(k2+1) + deltax2 * x(j)
            x4(node) = x4(k2+1) + deltax4 * x(j)
          end do
        end do

      end do

      x3 = 1._dp - x3(nn1row:1:-1)
      x4 = 1._dp - x4(nn1col:1:-1)

!     create straight lines in reference square [0,1]x[0,1]

      do i = 1, nn1row
        do j = 1, nn1col
          node = i + (j-1)*nn1row
          D = ( x1(i) - x3(i) ) * ( x4(j) - x2(j) ) - 1
          mesh%coor(node,1) = ( x4(j) * (x1(i) - x3(i)) - x1(i) ) / D
          mesh%coor(node,2) = ( x1(i) * (x4(j) - x2(j)) - x4(j) ) / D
        end do
      end do

    end if

    deallocate ( x1, x2, x3, x4 )

!   points

    mesh%npoints = 4

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn1row
    mesh%points(3) = nn1row * nn1col
    mesh%points(4) = nn1row * ( nn1col - 1 ) + 1

!   curves

    mesh%ncurves = 4

!   allocate memory for curves

    allocate( mesh%curves(MAXCURVES) )
    allocate( mesh%curves(1)%nodes(nn1row) )
    allocate( mesh%curves(2)%nodes(nn1col) )
    allocate( mesh%curves(3)%nodes(nn1row) )
    allocate( mesh%curves(4)%nodes(nn1col) )
    allocate( mesh%curves(1)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(2)%topology(inpelmc,mr,2) )
    allocate( mesh%curves(3)%topology(inpelmc,nr,2) )
    allocate( mesh%curves(4)%topology(inpelmc,mr,2) )

!   element type

    mesh%curves(:4)%element%elshape = 101
    mesh%curves(:4)%element%numnod = inpelmc
    mesh%curves(:4)%element%p(1,1) = p
    mesh%curves(:4)%element%p(2,1) = 0 ! equidistant
    mesh%curves(:4)%element%ndim = 2

!   other info of curves

    mesh%curves(:4)%ndim = 2
    mesh%curves(:4)%nblend = 0
    mesh%curves(1:4)%nnodes = [ nn1row, nn1col, nn1row, nn1col ]
    mesh%curves(1:4)%nelem = [ nr, mr, nr, mr ]
    mesh%curves(1:4)%n = mesh%curves(1:4)%nelem
    do i = 1, 4
      allocate(mesh%curves(i)%element_blend(0))
      allocate(mesh%curves(i)%nnodes_blend(2) )
      mesh%curves(i)%nnodes_blend = [0,mesh%curves(i)%nnodes]
    end do

!   nodes of curves

    mesh%curves(1)%nodes = [ (i, i=1,nn1row) ]                      ! curve 1
    mesh%curves(2)%nodes = [ (j*nn1row, j=1,nn1col) ]               ! curve 2
    mesh%curves(3)%nodes = [ ((nn1col-1)*nn1row+i, i=nn1row,1,-1) ] ! curve 3
    mesh%curves(4)%nodes = [ ((j-1)*nn1row+1, j=nn1col,1,-1) ]      ! curve 4

!   topology of elements on curves
    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        mesh%curves(curve)%topology(:,elem,1) = [ (p*(elem-1)+i,i=1,p+1) ]
        mesh%curves(curve)%topology(:,elem,2) = &
           mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
      end do
    end do

  end subroutine rectangle2d_triangle_ho

end module meshgen_basic_quadrilateral_m
