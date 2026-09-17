
! Copyright (C) 2006-2021 Martien A. Hulsen
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
! mesh on a hexahedron
!
! LIMITATIONS:
!   - curved boundaries have not yet been implemented.


module meshgen_basic_hexahedron_m

  use glob_defs_m
  use meshgen_basic_common_m

  implicit none

contains


! Simple meshgenerator for hexahedral regions.

  subroutine hexahedron ( mesh, mesh_options )

    type(mesh_t), intent(inout) :: mesh

!   Options for generating the mesh.
!   See type definition for possibilities and defaults.
    type(meshgen_options_t), intent(in) :: mesh_options


!   A basic mesh structure is generated (topology, coordinates etc). Other
!   parts of the mesh structure still need to be filled.


    real(dp), allocatable, dimension(:,:) :: tmp

!   test mesh

    if ( mesh%meshgen ) then
      write(*,'(/3(a/))') 'Error hexahedron:', &
        ' mesh is not empty', &
        ' either use a new mesh_t variable or delete old mesh '
      stop
    end if

!   test mesh_options

    if ( mesh_options%rx < 1 .or. mesh_options%ry < 1 .or. &
         mesh_options%rz < 1 ) then
      write(*,'(/2(a/))') 'Error hexahedron:', &
        ' rx, ry or rz in mesh_options must not be smaller than 1'
      stop
    end if

!   fill mesh

    select case ( mesh_options%elshape )

      case(-1) ! element shape probably not set

        write(*,'(/a/)') 'Error in hexahedron: element shape not set'
        stop

      case(13) ! eight-node hexahedron

        call cuboid_hexa8 ( mesh, mesh_options%nx, mesh_options%ny, &
          mesh_options%nz, mesh_options%rx, mesh_options%ry, mesh_options%rz, &
          mesh_options%ratio(1:12), mesh_options%factor(1:12) )

      case(14) ! 27-node hexahedron

        call cuboid_hexa27 ( mesh, mesh_options%nx, mesh_options%ny, &
          mesh_options%nz, mesh_options%rx, mesh_options%ry, mesh_options%rz, &
          mesh_options%ratio(1:12), mesh_options%factor(1:12) )

      case default

        write(*,'(/a,i0/)') 'Error: Wrong element shape: ', mesh_options%elshape
        stop

    end select


!   global element shape

    select case ( mesh_options%elshape )

      case(13,14,17,31)

        mesh%element(:)%globalshape = 'hexahedron'

      case(15,16,18)

        mesh%element(:)%globalshape = 'tetrahedron'

      case default

        call errormsg_case_default ( 'hexahedron', 'mesh_options%elshape', &
          int_value=mesh_options%elshape )

    end select


!   choose region

    select case ( mesh_options%regionshape )

      case(4) ! cuboid

!       scale coordinates

        mesh%coor(:,1) = mesh%coor(:,1) * mesh_options%lx + mesh_options%ox
        mesh%coor(:,2) = mesh%coor(:,2) * mesh_options%ly + mesh_options%oy
        mesh%coor(:,3) = mesh%coor(:,3) * mesh_options%lz + mesh_options%oz

      case(5) ! hexahedron

        allocate ( tmp(mesh%nnodes,3) )

        call straight_hexahedron ( tmp )

        mesh%coor = tmp

        deallocate ( tmp )

      case default

        write(*,'(/a,i0/)') &
         'Error: Wrong region shape: ', mesh_options%regionshape
        stop

    end select

!   type of elements on curves

    mesh%curves(1:mesh%ncurves)%element%globalshape = 'line'

!   type of elements on surfaces

    select case ( mesh_options%elshape )

      case(13,14,17,31)

        mesh%surfaces(1:mesh%nsurfaces)%element%globalshape = 'quadrilateral'

      case(15,16,18)

        mesh%surfaces(1:mesh%nsurfaces)%element%globalshape = 'triangle'

      case default

        call errormsg_case_default ( 'hexahedron', 'mesh_options%elshape', &
          int_value=mesh_options%elshape )

    end select

!   allow volumes in 3D

    allocate( mesh%volumes(MAXVOLUMES) )

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

    subroutine straight_hexahedron ( coor )

      real(dp), dimension(:,:), intent(out) :: coor

      integer :: i

!     x = sum_k x_k phi_k(xi,eta,zeta)
!     with xi = coor(1), eta = coor(2), zeta = coor(3) from cuboid

      do i = 1, 3
        coor(:,i) = &
         mesh_options%x3d(1,i) * &
          (1 - mesh%coor(:,1)) * (1 - mesh%coor(:,2)) * (1 - mesh%coor(:,3)) + &
         mesh_options%x3d(2,i) * &
               mesh%coor(:,1)  * (1 - mesh%coor(:,2)) * (1 - mesh%coor(:,3)) + &
         mesh_options%x3d(3,i) * &
               mesh%coor(:,1)  *      mesh%coor(:,2)  * (1 - mesh%coor(:,3)) + &
         mesh_options%x3d(4,i) * &
          (1 - mesh%coor(:,1)) *      mesh%coor(:,2)  * (1 - mesh%coor(:,3)) + &
         mesh_options%x3d(5,i) * &
          (1 - mesh%coor(:,1)) * (1 - mesh%coor(:,2)) *      mesh%coor(:,3)  + &
         mesh_options%x3d(6,i) * &
               mesh%coor(:,1)  * (1 - mesh%coor(:,2)) *      mesh%coor(:,3)  + &
         mesh_options%x3d(7,i) * &
               mesh%coor(:,1)  *      mesh%coor(:,2)  *      mesh%coor(:,3)  + &
         mesh_options%x3d(8,i) * &
          (1 - mesh%coor(:,1)) *      mesh%coor(:,2)  *      mesh%coor(:,3)
      end do

    end subroutine straight_hexahedron

  end subroutine hexahedron


! Generate a mesh on region [0,1]x[0,1]x[0,1] using hexahedrons with 8 nodes.

  subroutine cuboid_hexa8 ( mesh, n1, n2, n3, r1, r2, r3, ratio, factor )

    type(mesh_t), intent(out) :: mesh

!   number of elements along the sides and refinement factors
!   The ACTUAL number of elements is n1*r1 x n2*r2 x n3*r3.
    integer, intent(in) :: n1, n2, n3
    integer, intent(in) :: r1, r2, r3

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, dimension(:), intent(in) :: ratio
    real(dp), dimension(:), intent(in) :: factor

!   Generate a mesh on region  [0,1]x[0,1]x[0,1] using hexahedrons with 8 nodes.
!
!   Also generated are eight points, twelve curves and six surfaces.
!
!   Points:
!
!        8-----------7
!       /|          /|
!      / |         / |
!     5-----------6  |
!     |  |        |  |
!     |  4--------|--3
!     | /         | /
!     |/          |/
!     1-----------2
!
!  Curves:                    Surfaces:
!    C1=P1-P2                   S1=P1-P2-P3-P4
!    C2=P2-P3                   S2=P1-P2-P6-P5
!    C3=P4-P3                   S3=P2-P3-P7-P6
!    C4=P1-P4                   S4=P4-P3-P7-P8
!    C5=P1-P5                   S5=P1-P4-P8-P5
!    C6=P2-P6                   S6=P5-P6-P7-P8
!    C7=P3-P7
!    C8=P4-P8
!    C9=P6-P6
!    C10=P6-P7
!    C11=P8-P7
!    C12=P5-P8
!
!  NOTE: all curves are directed in positive coordinate direction
!

    integer, parameter :: inpelm = 8 ! number of nodal points per element
    integer, parameter :: inpelmc = 2 ! number of nodal points per element on
                                      ! curves
    integer, parameter :: inpelms = 4 ! number of nodal points per element on
                                      ! surfaces

    integer :: elnodes(inpelm), elsnodes(inpelms), i, j, k, e
    integer :: nn1, nn2, nn3, nn4, nn1row, nn1col, nn1hgt, nn1pln, nnodes
    integer :: node, elem, curve, surface, nr1, nr2, nr3
    real(dp) :: deltax, deltay, deltaz, D, xi, eta, zeta
    real(dp) :: deltax1, deltax2, deltax3, deltax4
    real(dp), allocatable, dimension(:) :: c1, c2, c3, c4, c5, c6, c7, c8, &
      c9, c10, c11, c12
    real(dp), allocatable, dimension(:,:,:) :: s1, s2, s3, s4, s5, s6

    allocate ( c1(n1*r1+1), c2(n2*r2+1), c3(n1*r1+1), c4(n2*r2+1) )
    allocate ( c5(n3*r3+1), c6(n3*r3+1), c7(n3*r3+1), c8(n3*r3+1) )
    allocate ( c9(n1*r1+1), c10(n2*r2+1), c11(n1*r1+1), c12(n2*r2+1) )
    allocate ( s1(n1*r1+1,n2*r2+1,2), s2(n1*r1+1,n3*r3+1,2), &
               s3(n2*r2+1,n3*r3+1,2) )
    allocate ( s4(n1*r1+1,n3*r3+1,2), s5(n2*r2+1,n3*r3+1,2),  &
               s6(n1*r1+1,n2*r2+1,2) )


    nr1 = n1 * r1
    nr2 = n2 * r2
    nr3 = n3 * r3

    nn1row = nr1+1
    nn1col = nr2+1
    nn1hgt = nr3+1
    nn1pln = nn1row * nn1col
    nnodes = nn1pln * nn1hgt

!   internal mesh

    mesh%ndim   = 3
    mesh%nelem  = nr1*nr2*nr3
    mesh%nnodes = nnodes
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

    mesh%element(1)%elshape = 13
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 3

!   topology

    mesh%grpnumel(1) = mesh%nelem

    do i = 1, nr1
      do j = 1, nr2
        do k = 1, nr3

!       number of nodes before element (i,j,k)
        nn1 = (k-1)*nn1pln + (j-1)*nn1row + i-1
!       number of nodes before element (i,j,k) + 1 row
        nn2 = nn1 + nn1row
!       number of nodes before element (i,j,k) + 1 plane
        nn3 = nn1 + nn1pln
!       number of nodes before element (i,j,k) + 1 plane + 1 row
        nn4 = nn3 + nn1row

        e = i + (j-1)*nr1 + (k-1)*nr1*nr2   ! element number

        elnodes = [ (nn1+i,i=1,2), (nn2+i,i=2,1,-1), &
                    (nn3+i,i=1,2), (nn4+i,i=2,1,-1) ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

        end do
      end do
    end do

!   coordinates

    deltax = 1.0_dp / nr1
    deltay = 1.0_dp / nr2
    deltaz = 1.0_dp / nr3

    if ( all(ratio == 0) ) then
!     equidistant
      do i = 1, nn1row
        do j = 1, nn1col
          do k = 1, nn1hgt
            node = i + (j-1)*nn1row + (k-1)*nn1pln
            mesh%coor(node,1) = (i-1) * deltax
            mesh%coor(node,2) = (j-1) * deltay
            mesh%coor(node,3) = (k-1) * deltaz
          end do
        end do
      end do
    else
!     non-equidistant
      call distribute_elements ( n1, c1(::r1), ratio(1), factor(1) )
      call distribute_elements ( n2, c2(::r2), ratio(2), factor(2) )
      call distribute_elements ( n1, c3(::r1), ratio(3), factor(3) )
      call distribute_elements ( n2, c4(::r2), ratio(4), factor(4) )
      call distribute_elements ( n3, c5(::r3), ratio(5), factor(5) )
      call distribute_elements ( n3, c6(::r3), ratio(6), factor(6) )
      call distribute_elements ( n3, c7(::r3), ratio(7), factor(7) )
      call distribute_elements ( n3, c8(::r3), ratio(8), factor(8) )
      call distribute_elements ( n1, c9(::r1), ratio(9), factor(9) )
      call distribute_elements ( n2, c10(::r2), ratio(10), factor(10) )
      call distribute_elements ( n1, c11(::r1), ratio(11), factor(11) )
      call distribute_elements ( n2, c12(::r2), ratio(12), factor(12) )
!     refine
      do elem = 1, n1
        k = ( elem - 1 ) * r1  ! nodes before element elem
        deltax1 = (c1(k+r1+1)-c1(k+1))/r1
        deltax2 = (c3(k+r1+1)-c3(k+1))/r1
        deltax3 = (c9(k+r1+1)-c9(k+1))/r1
        deltax4 = (c11(k+r1+1)-c11(k+1))/r1
        do node = k+2, k+r1
          c1(node) = c1(node-1) + deltax1
          c3(node) = c3(node-1) + deltax2
          c9(node) = c9(node-1) + deltax3
          c11(node) = c11(node-1) + deltax4
        end do
      end do
      do elem = 1, n2
        k = ( elem - 1 ) * r2  ! nodes before element elem
        deltax1 = (c2(k+r2+1)-c2(k+1))/r2
        deltax2 = (c4(k+r2+1)-c4(k+1))/r2
        deltax3 = (c10(k+r2+1)-c10(k+1))/r2
        deltax4 = (c12(k+r2+1)-c12(k+1))/r2
        do node = k+2, k+r2
          c2(node) = c2(node-1) + deltax1
          c4(node) = c4(node-1) + deltax2
          c10(node) = c10(node-1) + deltax3
          c12(node) = c12(node-1) + deltax4
        end do
      end do
      do elem = 1, n3
        k = ( elem - 1 ) * r3  ! nodes before element elem
        deltax1 = (c5(k+r3+1)-c5(k+1))/r3
        deltax2 = (c6(k+r3+1)-c6(k+1))/r3
        deltax3 = (c7(k+r3+1)-c7(k+1))/r3
        deltax4 = (c8(k+r3+1)-c8(k+1))/r3
        do node = k+2, k+r3
          c5(node) = c5(node-1) + deltax1
          c6(node) = c6(node-1) + deltax2
          c7(node) = c7(node-1) + deltax3
          c8(node) = c8(node-1) + deltax4
        end do
      end do
!     create straight lines in reference square [0,1]x[0,1] of surfaces
      call quadmesh ( nn1row, nn1col, c1, c2, c3, c4, s1 )
      call quadmesh ( nn1row, nn1hgt, c1, c6, c9, c5, s2 )
      call quadmesh ( nn1col, nn1hgt, c2, c7, c10, c6, s3 )
      call quadmesh ( nn1row, nn1hgt, c3, c7, c11, c8, s4 )
      call quadmesh ( nn1col, nn1hgt, c4, c8, c12, c5, s5 )
      call quadmesh ( nn1row, nn1col, c9, c10, c11, c12, s6 )
!     generate coordinates
      do i = 1, nn1row
        do j = 1, nn1col
          do k = 1, nn1hgt
            xi   = (i-1) * deltax
            eta  = (j-1) * deltay
            zeta = (k-1) * deltaz
            node = i + (j-1)*nn1row + (k-1)*nn1pln
            mesh%coor(node,1) = s1(i,j,1) * (1-zeta) + s6(i,j,1) * zeta &
                                  + s2(i,k,1) * (1-eta) + s4(i,k,1) * eta  &
                                  - c1(i) * (1-eta)*(1-zeta) &
                                  - c3(i) * eta*(1-zeta) &
                                  - c11(i) * eta*zeta &
                                  - c9(i) * (1-eta)*zeta
            mesh%coor(node,2) = s1(i,j,2) * (1-zeta) + s6(i,j,2) * zeta &
                                  + s5(j,k,1) * (1-xi) + s3(j,k,1) * xi  &
                                  - c4(j) * (1-xi)*(1-zeta) &
                                  - c2(j) * xi*(1-zeta) &
                                  - c10(j) * xi*zeta &
                                  - c12(j) * (1-xi)*zeta
            mesh%coor(node,3) = s2(i,k,2) * (1-eta) + s4(i,k,2) * eta &
                                  + s5(j,k,2) * (1-xi) + s3(j,k,2) * xi  &
                                  - c5(k) * (1-xi)*(1-eta) &
                                  - c6(k) * xi*(1-eta) &
                                  - c7(k) * xi*eta &
                                  - c8(k) * (1-xi)*eta
          end do
        end do
      end do
    end if

    deallocate ( c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12 )
    deallocate ( s1, s2, s3, s4, s5, s6 )


!   NOTE: the rest of the routine is identical to cubic_hexa27: should be merged

!   points

    mesh%npoints = 8

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn1row
    mesh%points(3) = nn1pln
    mesh%points(4) = nn1pln - nn1row + 1
    mesh%points(5) = nnodes - nn1pln + 1
    mesh%points(6) = nnodes - nn1pln + nn1row
    mesh%points(7) = nnodes
    mesh%points(8) = nnodes - nn1row + 1

!   curves

    mesh%ncurves = 12

!   allocate memory for curves

    allocate( mesh%curves(MAXCURVES) )
    allocate( mesh%curves(1)%nodes(nn1row) )
    allocate( mesh%curves(2)%nodes(nn1col) )
    allocate( mesh%curves(3)%nodes(nn1row) )
    allocate( mesh%curves(4)%nodes(nn1col) )
    allocate( mesh%curves(5)%nodes(nn1hgt) )
    allocate( mesh%curves(6)%nodes(nn1hgt) )
    allocate( mesh%curves(7)%nodes(nn1hgt) )
    allocate( mesh%curves(8)%nodes(nn1hgt) )
    allocate( mesh%curves(9)%nodes(nn1row) )
    allocate( mesh%curves(10)%nodes(nn1col) )
    allocate( mesh%curves(11)%nodes(nn1row) )
    allocate( mesh%curves(12)%nodes(nn1col) )
    allocate( mesh%curves(1)%topology(inpelmc,nr1,2) )
    allocate( mesh%curves(2)%topology(inpelmc,nr2,2) )
    allocate( mesh%curves(3)%topology(inpelmc,nr1,2) )
    allocate( mesh%curves(4)%topology(inpelmc,nr2,2) )
    allocate( mesh%curves(5)%topology(inpelmc,nr3,2) )
    allocate( mesh%curves(6)%topology(inpelmc,nr3,2) )
    allocate( mesh%curves(7)%topology(inpelmc,nr3,2) )
    allocate( mesh%curves(8)%topology(inpelmc,nr3,2) )
    allocate( mesh%curves(9)%topology(inpelmc,nr1,2) )
    allocate( mesh%curves(10)%topology(inpelmc,nr2,2) )
    allocate( mesh%curves(11)%topology(inpelmc,nr1,2) )
    allocate( mesh%curves(12)%topology(inpelmc,nr2,2) )

!   element type

    mesh%curves(:12)%element%elshape = 1
    mesh%curves(:12)%element%numnod = inpelmc
    mesh%curves(:12)%element%ndim = 3

!   other info of curves

    mesh%curves(:12)%ndim = 3
    mesh%curves(:12)%nblend = 0
    do i = 1, 12
      mesh%curves(i)%nnodes = size(mesh%curves(i)%nodes)
      mesh%curves(i)%nelem = size(mesh%curves(i)%topology,2)
      mesh%curves(i)%n = mesh%curves(i)%nelem
      allocate(mesh%curves(i)%element_blend(0))
      allocate(mesh%curves(i)%nnodes_blend(2) )
      mesh%curves(i)%nnodes_blend = [0,mesh%curves(i)%nnodes]
    end do

!   nodes of curves

    mesh%curves(1)%nodes = [ (i, i=1,nn1row) ]
    mesh%curves(2)%nodes = [ (j*nn1row, j=1,nn1col) ]
    mesh%curves(3)%nodes = [ (nn1pln-nn1row+i, i=1,nn1row) ]
    mesh%curves(4)%nodes = [ ((j-1)*nn1row+1, j=1,nn1col) ]
    mesh%curves(5)%nodes = [ ((k-1)*nn1pln+1, k=1,nn1hgt) ]
    mesh%curves(6)%nodes = mesh%curves(5)%nodes + nn1row - 1
    mesh%curves(7)%nodes = mesh%curves(5)%nodes + nn1pln - 1
    mesh%curves(8)%nodes = mesh%curves(5)%nodes + nn1pln - nn1row
    mesh%curves(9)%nodes = mesh%curves(1)%nodes + nnodes - nn1pln
    mesh%curves(10)%nodes = mesh%curves(2)%nodes + nnodes - nn1pln
    mesh%curves(11)%nodes = mesh%curves(3)%nodes + nnodes - nn1pln
    mesh%curves(12)%nodes = mesh%curves(4)%nodes + nnodes - nn1pln

!   topology of elements on curves

    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        mesh%curves(curve)%topology(:,elem,1) = &
                                   [ ((inpelmc-1)*(elem-1)+i,i=1,inpelmc) ]
        mesh%curves(curve)%topology(:,elem,2) = &
           mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
      end do
    end do

!   surfaces

    mesh%nsurfaces = 6

!   allocate memory for surfaces

    allocate( mesh%surfaces(MAXSURFACES) )
    allocate( mesh%surfaces(1)%nodes(nn1pln) )
    allocate( mesh%surfaces(2)%nodes(nn1row*nn1hgt) )
    allocate( mesh%surfaces(3)%nodes(nn1col*nn1hgt) )
    allocate( mesh%surfaces(4)%nodes(nn1row*nn1hgt) )
    allocate( mesh%surfaces(5)%nodes(nn1col*nn1hgt) )
    allocate( mesh%surfaces(6)%nodes(nn1pln) )
    allocate( mesh%surfaces(1)%topology(inpelms,nr1*nr2,2) )
    allocate( mesh%surfaces(2)%topology(inpelms,nr1*nr3,2) )
    allocate( mesh%surfaces(3)%topology(inpelms,nr2*nr3,2) )
    allocate( mesh%surfaces(4)%topology(inpelms,nr1*nr3,2) )
    allocate( mesh%surfaces(5)%topology(inpelms,nr2*nr3,2) )
    allocate( mesh%surfaces(6)%topology(inpelms,nr1*nr2,2) )

!   element type

    mesh%surfaces(:6)%element%elshape = 5
    mesh%surfaces(:6)%element%numnod = inpelms
    mesh%surfaces(:6)%element%ndim = 3
    mesh%surfaces(:6)%nblend = 0

!   other info of surfaces

    mesh%surfaces(:6)%ndim = 3
    do i = 1, 6
      mesh%surfaces(i)%nnodes = size(mesh%surfaces(i)%nodes)
      mesh%surfaces(i)%nelem = size(mesh%surfaces(i)%topology,2)
      allocate(mesh%surfaces(i)%element_blend(0))
      allocate(mesh%surfaces(i)%nnodes_blend(2) )
      mesh%surfaces(i)%nnodes_blend = [0,mesh%surfaces(i)%nnodes]
    end do

!   surface 1

    mesh%surfaces(1)%nodes = [ (i, i=1,nn1pln) ]

    call surtop ( srf=1, n1=nr1, n2=nr2, nn1row=nn1row )

    mesh%surfaces(1)%n = nr1
    mesh%surfaces(1)%m = nr2

!   surface 2

    mesh%surfaces(2)%nodes = &
       [ ((i+(j-1)*nn1pln, i=1,nn1row), j=1,nn1hgt) ]

    call surtop ( srf=2, n1=nr1, n2=nr3, nn1row=nn1row )

    mesh%surfaces(2)%n = nr1
    mesh%surfaces(2)%m = nr3

!   surface 3

    mesh%surfaces(3)%nodes = &
       [ ((nn1row*i+(j-1)*nn1pln, i=1,nn1col), j=1,nn1hgt) ]

    call surtop ( srf=3, n1=nr2, n2=nr3, nn1row=nn1col )

    mesh%surfaces(3)%n = nr2
    mesh%surfaces(3)%m = nr3

!   surface 4

    mesh%surfaces(4)%nodes = mesh%surfaces(2)%nodes + nn1pln - nn1row
    mesh%surfaces(4)%topology(:,:,1) = mesh%surfaces(2)%topology(:,:,1)

    mesh%surfaces(4)%n = mesh%surfaces(2)%n
    mesh%surfaces(4)%m = mesh%surfaces(2)%m

!   surface 5

    mesh%surfaces(5)%nodes = mesh%surfaces(3)%nodes - nn1row + 1
    mesh%surfaces(5)%topology(:,:,1) = mesh%surfaces(3)%topology(:,:,1)

    mesh%surfaces(5)%n = mesh%surfaces(3)%n
    mesh%surfaces(5)%m = mesh%surfaces(3)%m

!   surface 6

    mesh%surfaces(6)%nodes = mesh%surfaces(1)%nodes + nnodes - nn1pln
    mesh%surfaces(6)%topology(:,:,1) = mesh%surfaces(1)%topology(:,:,1)

    mesh%surfaces(6)%n = mesh%surfaces(1)%n
    mesh%surfaces(6)%m = mesh%surfaces(1)%m

!   topology of elements on surfaces in global nodes

    do surface = 1, mesh%nsurfaces
      do elem = 1, mesh%surfaces(surface)%nelem
        mesh%surfaces(surface)%topology(:,elem,2) = &
        mesh%surfaces(surface)%nodes(mesh%surfaces(surface)%topology(:,elem,1))
      end do
    end do

  contains


!   create mesh of straight lines in reference square [0,1]x[0,1]

    subroutine quadmesh ( nn1row, nn1col, c1, c2, c3, c4, s )

      integer, intent(in) :: nn1row, nn1col
      real(dp), dimension(:), intent(in) :: c1, c2, c3, c4
      real(dp), dimension(:,:,:), intent(out) :: s

      integer :: i, j

      do i = 1, nn1row
        do j = 1, nn1col
          D = ( c1(i) - c3(i) ) * ( c4(j) - c2(j) ) - 1
          s(i,j,1) = ( c4(j) * (c1(i) - c3(i)) - c1(i) ) / D
          s(i,j,2) = ( c1(i) * (c4(j) - c2(j)) - c4(j) ) / D
        end do
      end do

    end subroutine quadmesh


!   surface topology

    subroutine surtop ( srf, n1, n2, nn1row )

      integer, intent(in) :: srf, n1, n2, nn1row

      integer :: i, j

      do i = 1, n1
        do j = 1, n2

          nn1 = (j-1)*nn1row + i-1 ! number of nodes before element (i,j)
          nn2 = nn1 + nn1row   ! number of nodes before element (i,j) + 1 row

          e = i + (j-1)*n1   ! element number

          elsnodes = [ (nn1+i,i=1,2), (nn2+i,i=2,1,-1) ]

          mesh%surfaces(srf)%topology(:,e,1) = elsnodes(1:inpelms)

        end do
      end do

    end subroutine surtop

  end subroutine cuboid_hexa8


! Generate a mesh on region [0,1]x[0,1]x[0,1] using hexahedrons with 27 nodes.

  subroutine cuboid_hexa27 ( mesh, n1, n2, n3, r1, r2, r3, ratio, factor )

    type(mesh_t), intent(out) :: mesh

!   number of elements along the sides and refinement factors
!   The ACTUAL number of elements is n1*r1 x n2*r2 x n3*r3.
    integer, intent(in) :: n1, n2, n3
    integer, intent(in) :: r1, r2, r3

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, dimension(:), intent(in) :: ratio
    real(dp), dimension(:), intent(in) :: factor

!   Generate a mesh on region [0,1]x[0,1]x[0,1] using hexahedrons with 8 nodes.
!
!   Also generated are eight points, twelve curves and six surfaces.
!
!   Points:
!
!        8-----------7
!       /|          /|
!      / |         / |
!     5-----------6  |
!     |  |        |  |
!     |  4--------|--3
!     | /         | /
!     |/          |/
!     1-----------2
!
!  Curves:                    Surfaces:
!    C1=P1-P2                   S1=P1-P2-P3-P4
!    C2=P2-P3                   S2=P1-P2-P6-P5
!    C3=P4-P3                   S3=P2-P3-P7-P6
!    C4=P1-P4                   S4=P4-P3-P7-P8
!    C5=P1-P5                   S5=P1-P4-P8-P5
!    C6=P2-P6                   S6=P5-P6-P7-P8
!    C7=P3-P7
!    C8=P4-P8
!    C9=P6-P6
!    C10=P6-P7
!    C11=P8-P7
!    C12=P5-P8
!
!  NOTE: all curves are directed in positive coordinate direction
!

    integer, parameter :: inpelm = 27 ! number of nodal points per element
    integer, parameter :: inpelmc = 3 ! number of nodal points per element on
                                      ! curves
    integer, parameter :: inpelms = 9 ! number of nodal points per element on
                                      ! surfaces

    integer :: elnodes(inpelm), elsnodes(inpelms), i, j, k, e
    integer :: nn1, nn2, nn3, nn4, nn1row, nn1col, nn1hgt, nn1pln, nnodes
    integer :: nn5, nn6, nn7, nn8, nn9
    integer :: node, elem, curve, surface, nr1, nr2, nr3
    real(dp) :: deltax, deltay, deltaz, D, xi, eta, zeta
    real(dp) :: deltax1, deltax2, deltax3, deltax4
    real(dp), allocatable, dimension(:) :: c1, c2, c3, c4, c5, c6, c7, c8, &
      c9, c10, c11, c12
    real(dp), allocatable, dimension(:,:,:) :: s1, s2, s3, s4, s5, s6

    allocate ( c1(2*n1*r1+1), c2(2*n2*r2+1), c3(2*n1*r1+1), c4(2*n2*r2+1) )
    allocate ( c5(2*n3*r3+1), c6(2*n3*r3+1), c7(2*n3*r3+1), c8(2*n3*r3+1) )
    allocate ( c9(2*n1*r1+1), c10(2*n2*r2+1), c11(2*n1*r1+1), c12(2*n2*r2+1) )
    allocate ( s1(2*n1*r1+1,2*n2*r2+1,2), s2(2*n1*r1+1,2*n3*r3+1,2), &
                s3(2*n2*r2+1,2*n3*r3+1,2) )
    allocate ( s4(2*n1*r1+1,2*n3*r3+1,2), s5(2*n2*r2+1,2*n3*r3+1,2), &
                s6(2*n1*r1+1,2*n2*r2+1,2) )


    nr1 = n1 * r1
    nr2 = n2 * r2
    nr3 = n3 * r3

    nn1row = 2*nr1+1
    nn1col = 2*nr2+1
    nn1hgt = 2*nr3+1
    nn1pln = nn1row * nn1col
    nnodes = nn1pln * nn1hgt

!   internal mesh

    mesh%ndim   = 3
    mesh%nelem  = nr1*nr2*nr3
    mesh%nnodes = nnodes
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

    mesh%element(1)%elshape = 14
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 3

!   topology

    mesh%grpnumel(1) = mesh%nelem

    do i = 1, nr1
      do j = 1, nr2
        do k = 1, nr3

!       number of nodes before element (i,j,k)
        nn1 = 2*(k-1)*nn1pln + 2*(j-1)*nn1row + 2*(i-1)
!       number of nodes before element (i,j,k) + 1 row
        nn2 = nn1 + nn1row
!       number of nodes before element (i,j,k) + 2 rows
        nn3 = nn2 + nn1row
!       number of nodes before element (i,j,k) + 1 plane
        nn4 = nn1 + nn1pln
!       number of nodes before element (i,j,k) + 1 plane + 1 row
        nn5 = nn4 + nn1row
!       number of nodes before element (i,j,k) + 1 plane + 1 row
        nn6 = nn5 + nn1row
!       number of nodes before element (i,j,k) + 2 planes
        nn7 = nn4 + nn1pln
!       number of nodes before element (i,j,k) + 2 planes + 1 row
        nn8 = nn7 + nn1row
!       number of nodes before element (i,j,k) + 2 planes + 2 rows
        nn9 = nn8 + nn1row

        e = i + (j-1)*nr1 + (k-1)*nr1*nr2   ! element number

        elnodes = [ (nn1+i,i=1,3), (nn2+i,i=1,3), (nn3+i,i=1,3), &
                    (nn4+i,i=1,3), (nn5+i,i=1,3), (nn6+i,i=1,3), &
                    (nn7+i,i=1,3), (nn8+i,i=1,3), (nn9+i,i=1,3) ]

        mesh%topology(1)%a(:,e) = elnodes(1:inpelm)

        end do
      end do
    end do

!   coordinates

    deltax = 0.5_dp / nr1
    deltay = 0.5_dp / nr2
    deltaz = 0.5_dp / nr3

    if ( all(ratio == 0) )then
!     equidistant
      do i = 1, nn1row
        do j = 1, nn1col
          do k = 1, nn1hgt
            node = i + (j-1)*nn1row + (k-1)*nn1pln
            mesh%coor(node,1) = (i-1) * deltax
            mesh%coor(node,2) = (j-1) * deltay
            mesh%coor(node,3) = (k-1) * deltaz
          end do
        end do
      end do
    else
!     non-equidistant
      call distribute_elements ( n1, c1(::2*r1), ratio(1), factor(1) )
      call distribute_elements ( n2, c2(::2*r2), ratio(2), factor(2) )
      call distribute_elements ( n1, c3(::2*r1), ratio(3), factor(3) )
      call distribute_elements ( n2, c4(::2*r2), ratio(4), factor(4) )
      call distribute_elements ( n3, c5(::2*r3), ratio(5), factor(5) )
      call distribute_elements ( n3, c6(::2*r3), ratio(6), factor(6) )
      call distribute_elements ( n3, c7(::2*r3), ratio(7), factor(7) )
      call distribute_elements ( n3, c8(::2*r3), ratio(8), factor(8) )
      call distribute_elements ( n1, c9(::2*r1), ratio(9), factor(9) )
      call distribute_elements ( n2, c10(::2*r2), ratio(10), factor(10) )
      call distribute_elements ( n1, c11(::2*r1), ratio(11), factor(11) )
      call distribute_elements ( n2, c12(::2*r2), ratio(12), factor(12) )
!     mid-side nodes and refine
      do elem = 1, n1
        k = ( elem - 1 ) * 2*r1  ! nodes before element elem
        deltax1 = (c1(k+2*r1+1)-c1(k+1))/2/r1
        deltax2 = (c3(k+2*r1+1)-c3(k+1))/2/r1
        deltax3 = (c9(k+2*r1+1)-c9(k+1))/2/r1
        deltax4 = (c11(k+2*r1+1)-c11(k+1))/2/r1
        do node = k+2, k+2*r1
          c1(node) = c1(node-1) + deltax1
          c3(node) = c3(node-1) + deltax2
          c9(node) = c9(node-1) + deltax3
          c11(node) = c11(node-1) + deltax4
        end do
      end do
      do elem = 1, n2
        k = ( elem - 1 ) * 2*r2  ! nodes before element elem
        deltax1 = (c2(k+2*r2+1)-c2(k+1))/2/r2
        deltax2 = (c4(k+2*r2+1)-c4(k+1))/2/r2
        deltax3 = (c10(k+2*r2+1)-c10(k+1))/2/r2
        deltax4 = (c12(k+2*r2+1)-c12(k+1))/2/r2
        do node = k+2, k+2*r2
          c2(node) = c2(node-1) + deltax1
          c4(node) = c4(node-1) + deltax2
          c10(node) = c10(node-1) + deltax3
          c12(node) = c12(node-1) + deltax4
        end do
      end do
      do elem = 1, n3
        k = ( elem - 1 ) * 2*r3  ! nodes before element elem
        deltax1 = (c5(k+2*r3+1)-c5(k+1))/2/r3
        deltax2 = (c6(k+2*r3+1)-c6(k+1))/2/r3
        deltax3 = (c7(k+2*r3+1)-c7(k+1))/2/r3
        deltax4 = (c8(k+2*r3+1)-c8(k+1))/2/r3
        do node = k+2, k+2*r3
          c5(node) = c5(node-1) + deltax1
          c6(node) = c6(node-1) + deltax2
          c7(node) = c7(node-1) + deltax3
          c8(node) = c8(node-1) + deltax4
        end do
      end do
!     create straight lines in reference square [0,1]x[0,1] of surfaces
      call quadmesh ( nn1row, nn1col, c1, c2, c3, c4, s1 )
      call quadmesh ( nn1row, nn1hgt, c1, c6, c9, c5, s2 )
      call quadmesh ( nn1col, nn1hgt, c2, c7, c10, c6, s3 )
      call quadmesh ( nn1row, nn1hgt, c3, c7, c11, c8, s4 )
      call quadmesh ( nn1col, nn1hgt, c4, c8, c12, c5, s5 )
      call quadmesh ( nn1row, nn1col, c9, c10, c11, c12, s6 )
!     generate coordinates
      do i = 1, nn1row
        do j = 1, nn1col
          do k = 1, nn1hgt
            xi   = (i-1) * deltax
            eta  = (j-1) * deltay
            zeta = (k-1) * deltaz
            node = i + (j-1)*nn1row + (k-1)*nn1pln
            mesh%coor(node,1) = s1(i,j,1) * (1-zeta) + s6(i,j,1) * zeta &
                                  + s2(i,k,1) * (1-eta) + s4(i,k,1) * eta  &
                                  - c1(i) * (1-eta)*(1-zeta) &
                                  - c3(i) * eta*(1-zeta) &
                                  - c11(i) * eta*zeta &
                                  - c9(i) * (1-eta)*zeta
            mesh%coor(node,2) = s1(i,j,2) * (1-zeta) + s6(i,j,2) * zeta &
                                  + s5(j,k,1) * (1-xi) + s3(j,k,1) * xi  &
                                  - c4(j) * (1-xi)*(1-zeta) &
                                  - c2(j) * xi*(1-zeta) &
                                  - c10(j) * xi*zeta &
                                  - c12(j) * (1-xi)*zeta
            mesh%coor(node,3) = s2(i,k,2) * (1-eta) + s4(i,k,2) * eta &
                                  + s5(j,k,2) * (1-xi) + s3(j,k,2) * xi  &
                                  - c5(k) * (1-xi)*(1-eta) &
                                  - c6(k) * xi*(1-eta) &
                                  - c7(k) * xi*eta &
                                  - c8(k) * (1-xi)*eta
          end do
        end do
      end do
    end if

    deallocate ( c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12 )
    deallocate ( s1, s2, s3, s4, s5, s6 )


!   NOTE: the rest of the routine is identical to cubic_hexa8: should be merged

!   points

    mesh%npoints = 8

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn1row
    mesh%points(3) = nn1pln
    mesh%points(4) = nn1pln - nn1row + 1
    mesh%points(5) = nnodes - nn1pln + 1
    mesh%points(6) = nnodes - nn1pln + nn1row
    mesh%points(7) = nnodes
    mesh%points(8) = nnodes - nn1row + 1

!   curves

    mesh%ncurves = 12

!   allocate memory for curves

    allocate( mesh%curves(MAXCURVES) )
    allocate( mesh%curves(1)%nodes(nn1row) )
    allocate( mesh%curves(2)%nodes(nn1col) )
    allocate( mesh%curves(3)%nodes(nn1row) )
    allocate( mesh%curves(4)%nodes(nn1col) )
    allocate( mesh%curves(5)%nodes(nn1hgt) )
    allocate( mesh%curves(6)%nodes(nn1hgt) )
    allocate( mesh%curves(7)%nodes(nn1hgt) )
    allocate( mesh%curves(8)%nodes(nn1hgt) )
    allocate( mesh%curves(9)%nodes(nn1row) )
    allocate( mesh%curves(10)%nodes(nn1col) )
    allocate( mesh%curves(11)%nodes(nn1row) )
    allocate( mesh%curves(12)%nodes(nn1col) )
    allocate( mesh%curves(1)%topology(inpelmc,nr1,2) )
    allocate( mesh%curves(2)%topology(inpelmc,nr2,2) )
    allocate( mesh%curves(3)%topology(inpelmc,nr1,2) )
    allocate( mesh%curves(4)%topology(inpelmc,nr2,2) )
    allocate( mesh%curves(5)%topology(inpelmc,nr3,2) )
    allocate( mesh%curves(6)%topology(inpelmc,nr3,2) )
    allocate( mesh%curves(7)%topology(inpelmc,nr3,2) )
    allocate( mesh%curves(8)%topology(inpelmc,nr3,2) )
    allocate( mesh%curves(9)%topology(inpelmc,nr1,2) )
    allocate( mesh%curves(10)%topology(inpelmc,nr2,2) )
    allocate( mesh%curves(11)%topology(inpelmc,nr1,2) )
    allocate( mesh%curves(12)%topology(inpelmc,nr2,2) )

!   element type

    mesh%curves(:12)%element%elshape = 2
    mesh%curves(:12)%element%numnod = inpelmc
    mesh%curves(:12)%element%ndim = 3

!   other info of curves

    mesh%curves(:12)%ndim = 3
    mesh%curves(:12)%nblend = 0
    do i = 1, 12
      mesh%curves(i)%nnodes = size(mesh%curves(i)%nodes)
      mesh%curves(i)%nelem = size(mesh%curves(i)%topology,2)
      mesh%curves(i)%n = mesh%curves(i)%nelem
      allocate(mesh%curves(i)%element_blend(0))
      allocate(mesh%curves(i)%nnodes_blend(2) )
      mesh%curves(i)%nnodes_blend = [0,mesh%curves(i)%nnodes]
    end do

!   nodes of curves

    mesh%curves(1)%nodes = [ (i, i=1,nn1row) ]
    mesh%curves(2)%nodes = [ (j*nn1row, j=1,nn1col) ]
    mesh%curves(3)%nodes = [ (nn1pln-nn1row+i, i=1,nn1row) ]
    mesh%curves(4)%nodes = [ ((j-1)*nn1row+1, j=1,nn1col) ]
    mesh%curves(5)%nodes = [ ((k-1)*nn1pln+1, k=1,nn1hgt) ]
    mesh%curves(6)%nodes = mesh%curves(5)%nodes + nn1row - 1
    mesh%curves(7)%nodes = mesh%curves(5)%nodes + nn1pln - 1
    mesh%curves(8)%nodes = mesh%curves(5)%nodes + nn1pln - nn1row
    mesh%curves(9)%nodes = mesh%curves(1)%nodes + nnodes - nn1pln
    mesh%curves(10)%nodes = mesh%curves(2)%nodes + nnodes - nn1pln
    mesh%curves(11)%nodes = mesh%curves(3)%nodes + nnodes - nn1pln
    mesh%curves(12)%nodes = mesh%curves(4)%nodes + nnodes - nn1pln

!   topology of elements on curves

    do curve = 1, mesh%ncurves
      do elem = 1, mesh%curves(curve)%nelem
        mesh%curves(curve)%topology(:,elem,1) = &
                                   [ ((inpelmc-1)*(elem-1)+i,i=1,inpelmc) ]
        mesh%curves(curve)%topology(:,elem,2) = &
           mesh%curves(curve)%nodes( mesh%curves(curve)%topology(:,elem,1) )
      end do
    end do

!   surfaces

    mesh%nsurfaces = 6

!   allocate memory for surfaces

    allocate( mesh%surfaces(MAXSURFACES) )
    allocate( mesh%surfaces(1)%nodes(nn1pln) )
    allocate( mesh%surfaces(2)%nodes(nn1row*nn1hgt) )
    allocate( mesh%surfaces(3)%nodes(nn1col*nn1hgt) )
    allocate( mesh%surfaces(4)%nodes(nn1row*nn1hgt) )
    allocate( mesh%surfaces(5)%nodes(nn1col*nn1hgt) )
    allocate( mesh%surfaces(6)%nodes(nn1pln) )
    allocate( mesh%surfaces(1)%topology(inpelms,nr1*nr2,2) )
    allocate( mesh%surfaces(2)%topology(inpelms,nr1*nr3,2) )
    allocate( mesh%surfaces(3)%topology(inpelms,nr2*nr3,2) )
    allocate( mesh%surfaces(4)%topology(inpelms,nr1*nr3,2) )
    allocate( mesh%surfaces(5)%topology(inpelms,nr2*nr3,2) )
    allocate( mesh%surfaces(6)%topology(inpelms,nr1*nr2,2) )

!   element type

    mesh%surfaces(:6)%element%elshape = 6
    mesh%surfaces(:6)%element%numnod = inpelms
    mesh%surfaces(:6)%element%ndim = 3
    mesh%surfaces(:6)%nblend = 0

!   other info of surfaces

    mesh%surfaces(:6)%ndim = 3
    do i = 1, 6
      mesh%surfaces(i)%nnodes = size(mesh%surfaces(i)%nodes)
      mesh%surfaces(i)%nelem = size(mesh%surfaces(i)%topology,2)
      allocate(mesh%surfaces(i)%element_blend(0))
      allocate(mesh%surfaces(i)%nnodes_blend(2) )
      mesh%surfaces(i)%nnodes_blend = [0,mesh%surfaces(i)%nnodes]
    end do

!   surface 1

    mesh%surfaces(1)%nodes = [ (i, i=1,nn1pln) ]

    call surtop ( srf=1, n1=nr1, n2=nr2, nn1row=nn1row )

    mesh%surfaces(1)%n = nr1
    mesh%surfaces(1)%m = nr2

!   surface 2

    mesh%surfaces(2)%nodes = &
       [ ((i+(j-1)*nn1pln, i=1,nn1row), j=1,nn1hgt) ]

    call surtop ( srf=2, n1=nr1, n2=nr3, nn1row=nn1row )

    mesh%surfaces(2)%n = nr1
    mesh%surfaces(2)%m = nr3

!   surface 3

    mesh%surfaces(3)%nodes = &
       [ ((nn1row*i+(j-1)*nn1pln, i=1,nn1col), j=1,nn1hgt) ]

    call surtop ( srf=3, n1=nr2, n2=nr3, nn1row=nn1col )

    mesh%surfaces(3)%n = nr2
    mesh%surfaces(3)%m = nr3

!   surface 4

    mesh%surfaces(4)%nodes = mesh%surfaces(2)%nodes + nn1pln - nn1row
    mesh%surfaces(4)%topology(:,:,1) = mesh%surfaces(2)%topology(:,:,1)

    mesh%surfaces(4)%n = mesh%surfaces(2)%n
    mesh%surfaces(4)%m = mesh%surfaces(2)%m

!   surface 5

    mesh%surfaces(5)%nodes = mesh%surfaces(3)%nodes - nn1row + 1
    mesh%surfaces(5)%topology(:,:,1) = mesh%surfaces(3)%topology(:,:,1)

    mesh%surfaces(5)%n = mesh%surfaces(3)%n
    mesh%surfaces(5)%m = mesh%surfaces(3)%m

!   surface 6

    mesh%surfaces(6)%nodes = mesh%surfaces(1)%nodes + nnodes - nn1pln
    mesh%surfaces(6)%topology(:,:,1) = mesh%surfaces(1)%topology(:,:,1)

    mesh%surfaces(6)%n = mesh%surfaces(1)%n
    mesh%surfaces(6)%m = mesh%surfaces(1)%m

!   topology of elements on surfaces in global nodes

    do surface = 1, mesh%nsurfaces
      do elem = 1, mesh%surfaces(surface)%nelem
        mesh%surfaces(surface)%topology(:,elem,2) = &
        mesh%surfaces(surface)%nodes(mesh%surfaces(surface)%topology(:,elem,1))
      end do
    end do

  contains

    subroutine surtop ( srf, n1, n2, nn1row )

      integer, intent(in) :: srf, n1, n2, nn1row

      integer :: i, j

      do i = 1, n1
        do j = 1, n2

          nn1 = 2*(j-1)*nn1row + 2*(i-1) ! number of nodes before element (i,j)
          nn2 = nn1 + nn1row   ! number of nodes before element (i,j) + 1 row
          nn3 = nn2 + nn1row   ! number of nodes before element (i,j) + 2 rows

          e = i + (j-1)*n1   ! element number

          elsnodes = [ (nn1+i,i=1,3), nn2+3, (nn3+i,i=3,1,-1), (nn2+i,i=1,2) ]

          mesh%surfaces(srf)%topology(:,e,1) = elsnodes(1:inpelms)

        end do
      end do

    end subroutine surtop


!   create mesh of straight lines in reference square [0,1]x[0,1]

    subroutine quadmesh ( nn1row, nn1col, c1, c2, c3, c4, s )

      integer, intent(in) :: nn1row, nn1col
      real(dp), dimension(:), intent(in) :: c1, c2, c3, c4
      real(dp), dimension(:,:,:), intent(out) :: s

      integer :: i, j

      do i = 1, nn1row
        do j = 1, nn1col
          D = ( c1(i) - c3(i) ) * ( c4(j) - c2(j) ) - 1
          s(i,j,1) = ( c4(j) * (c1(i) - c3(i)) - c1(i) ) / D
          s(i,j,2) = ( c1(i) * (c4(j) - c2(j)) - c4(j) ) / D
        end do
      end do

    end subroutine quadmesh


  end subroutine cuboid_hexa27


end module meshgen_basic_hexahedron_m
