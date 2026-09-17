
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
! mesh generation on a line


module meshgen_basic_line_m

  use kind_defs_m
  use meshgen_basic_common_m

  implicit none


contains


! Simple meshgenerator for 1D lines

  subroutine line1d ( mesh, mesh_options )

    type(mesh_t), intent(inout) :: mesh

!   Options for generating the mesh.
!   See type definition for possibilities and defaults.
    type(meshgen_options_t), intent(in) :: mesh_options

!   A basic mesh structure is generated (topology, coordinates etc). Other
!   parts of the mesh structure still need to be filled.


!   test mesh

    if ( mesh%meshgen ) then
      write(*,'(/3(a/))') 'Error line1d:', &
        ' mesh is not empty', &
        ' either use a new mesh_t variable or delete old mesh '
      stop
    end if

!   test mesh_options

    if ( mesh_options%rx < 1 ) then
      write(*,'(/2(a/))') 'Error line1d:', &
        ' rx in mesh_options must not be smaller than 1'
      stop
    end if

!   fill mesh

    select case ( mesh_options%elshape )

      case(-1) ! element shape probably not set

        write(*,'(/a/)') 'Error in line1d: element shape not set'
        stop

      case(1) ! 2 node line element

        call line1d_2node ( mesh, mesh_options%nx, mesh_options%rx, &
          mesh_options%ratio(1), mesh_options%factor(1) )

      case(2) ! 3 node line element

        call line1d_3node ( mesh, mesh_options%nx, mesh_options%rx, &
          mesh_options%ratio(1), mesh_options%factor(1) )

      case(101) ! spectral line element

        if ( mesh_options%p < 1 ) then
          write(*,'(/2(a/))') 'Error line1d:', &
            ' for spectral elements p in mesh_options must be larger than 0'
          stop
        end if

        if ( .not. any( mesh_options%l == [ 0, 1 ] ) ) then
          write(*,'(/2(a/))') 'Error line1d:', &
            ' for high-order elements l in mesh_options must be either 0 or 1.'
          stop
        end if

        call line1d_ho ( mesh, mesh_options%nx, mesh_options%rx, &
          mesh_options%ratio(1), mesh_options%factor(1), mesh_options%p, &
          mesh_options%l )

      case default

        write(*,'(/a,i0/)') 'Error: Wrong element shape: ', mesh_options%elshape
        stop

    end select

!   scale coordinates

    mesh%coor(:,1) = mesh%coor(:,1) * mesh_options%lx + mesh_options%ox

!   allow curves

    allocate ( mesh%curves(MAXCURVES) )

!   no surfaces and volumes

    allocate ( mesh%surfaces(0), mesh%volumes(0) )

!   objects

    allocate ( mesh%objects(MAXOBJECTS) )

!   blocks

    allocate ( mesh%blocks(MAXBLOCKS) )

!   nodesets

    allocate ( mesh%nodesets(MAXNODESETS) )

!   elementsets

    allocate ( mesh%elementsets(MAXELEMENTSETS) )

!   global element shape

    mesh%element(:)%globalshape = 'line'

!   basic mesh has been generated

    mesh%meshgen = .true.

  end subroutine line1d


! Generate mesh on region [0,1] using line elements with 2 nodes.

  subroutine line1d_2node ( mesh, n, r, ratio, factor )

    type(mesh_t), intent(out) :: mesh

!   number of elements on the interval [0,1] and refinement factor
!   The ACTUAL number of elements is n*r
    integer, intent(in) :: n, r

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, intent(in) :: ratio
    real(dp), intent(in) :: factor

!   Generate a mesh on region [0,1] using line elements with 2 nodes.
!   The numbering is straightforward:
!   For example a 4 element mesh is numbered as follows
!
!   Nodes:
!
!       1 --  2 --  3 --  4 --  5
!
!   Elements:
!
!       x --  x --  x --  x --  x
!          1     2     3     4
!
!   Also generated are two points and the end of the interval
!
!      P1 --------------------- P2
!

    integer, parameter :: inpelm = 2 ! number of nodal points per element

    integer :: elem, node, nn, k, nr
    real(dp) :: deltax
    real(dp), allocatable, dimension(:) :: x

    allocate ( x(1:n*r+1) )


    nr = n*r

    nn = nr+1

!   internal mesh

    mesh%ndim   = 1
    mesh%nelem  = nr
    mesh%nnodes = nn
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

    mesh%element(1)%elshape = 1
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 1

!   topology

    mesh%grpnumel(1) = mesh%nelem

    do elem = 1, nr
      mesh%topology(1)%a(:,elem) = [ elem, elem + 1 ]
    end do

!   coordinates

    if ( ratio == 0 ) then
!     equidistant
      deltax = 1.0_dp / nr
      do node = 1, nn
        mesh%coor(node,1) = (node-1) * deltax
      end do
    else
!     non-equidistant
      call distribute_elements ( n, x(::r), ratio, factor )
!     refine
      do elem = 1, n
        k = ( elem - 1 ) * r  ! nodes before element elem
        deltax = (x(k+r+1)-x(k+1))/r
        do node = k+2, k+r
          x(node) = x(node-1) + deltax
        end do
      end do
      mesh%coor(1:nn,1) = x
    end if

    deallocate ( x )

!   points

    mesh%npoints = 2

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn

  end subroutine line1d_2node


! Generate mesh on region [0,1] using line elements with 3 nodes.

  subroutine line1d_3node ( mesh, n, r, ratio, factor )

    type(mesh_t), intent(out) :: mesh

!   number of elements on the interval [0,1] and refinement factor
!   The ACTUAL number of elements is n*r
    integer, intent(in) :: n, r

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, intent(in) :: ratio
    real(dp), intent(in) :: factor

!   Generate a mesh on region [0,1] using line elements with 3 nodes.
!   The numbering is straightforward:
!   For example a 2 element mesh is numbered as follows
!
!   Nodes:
!
!       1 --  2 --  3 --  4 --  5
!
!   Elements:
!
!       x --  x --  x --  x --  x
!             1           2
!
!   Also generated are two points and the end of the interval
!
!      P1 --------------------- P2
!

    integer, parameter :: inpelm = 3 ! number of nodal points per element

    integer :: elem, node, nn, nr, k
    real(dp) :: deltax
    real(dp), allocatable, dimension(:) :: x

    allocate ( x(1:2*n*r+1) )


    nr = n*r

    nn = 2*nr+1

!   internal mesh

    mesh%ndim   = 1
    mesh%nelem  = nr
    mesh%nnodes = nn
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

    mesh%element(1)%elshape = 2
    mesh%element(1)%numnod = inpelm
    mesh%element(1)%ndim = 1

!   topology

    mesh%grpnumel(1) = mesh%nelem

    do elem = 1, nr
      mesh%topology(1)%a(:,elem) = [ 2*elem-1, 2*elem, 2*elem + 1 ]
    end do

!   coordinates

    if ( ratio == 0 ) then
!     equidistant
      deltax = 0.5_dp / nr
      do node = 1, nn
        mesh%coor(node,1) = (node-1) * deltax
      end do
    else
!     non-equidistant
      call distribute_elements ( n, x(1:nn:2*r), ratio, factor )
!     mid-side nodes and refine
      do elem = 1, n
        k = ( elem - 1 ) * 2*r  ! nodes before element elem
        deltax = (x(k+2*r+1)-x(k+1))/2/r
        do node = k+2, k+2*r
          x(node) = x(node-1) + deltax
        end do
      end do
      mesh%coor(1:nn,1) = x
    end if

    deallocate ( x )

!   points

    mesh%npoints = 2

!   allocate memory for points

    allocate( mesh%points(MAXPOINTS) )

    mesh%points(1) = 1
    mesh%points(2) = nn

  end subroutine line1d_3node


! Generate mesh on region [0,1] using high-order line elements

  subroutine line1d_ho ( mesh, n, r, ratio, factor, p, l )

    use gauss_spectral_m

    type(mesh_t), intent(out) :: mesh

!   number of elements on the interval [0,1] and refinement factor
!   The ACTUAL number of elements is n*r
    integer, intent(in) :: n, r

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, intent(in) :: ratio
    real(dp), intent(in) :: factor

!   polynomial order
    integer, intent(in) :: p

!   layout: l=0: equal distribution l=1: Gauss-Lobatto
    integer, intent(in) :: l

!   simple mesh generator for a high-order line element with order p
!   on the standard mesh region [0,1]
!   the following numbering is used, with for example p=3 and number of
!   elements n=2
!
!   Nodes:
!
!      1 -- 2 -- 3 -- 4 -- 5 -- 6 -- 7
!
!   Elements:
!
!      x -- x -- x -- x -- x -- x -- x
!      |______1______||______2_______|
!
!   Also the two end points are generated
!
!      P1---------------------------P2
!

    integer:: inpelm ! number of nodal points of an element
    integer:: elem, node, nn, i, nr, j
    real(dp):: deltax, ltmp
    real(dp), allocatable, dimension(:) :: x, xr

    allocate ( x(p+1), xr(0:n) )


    inpelm = p + 1  !number of integration points of a standard element

    nr = n*r    ! the actual number of elements including refinement factor r
    nn = nr*p+1 ! total number of nodes

!   internal mesh
    mesh%ndim = 1       ! dimension is 1, because it is a line
    mesh%nelem = nr     ! number of elements
    mesh%nnodes = nn    ! number of nodes
    mesh%nelgrp = 1     ! number of elementgroups (standard nelgrp=1)

!   allocate memory for internal mesh
    allocate( mesh%element (mesh%nelgrp) )
    allocate( mesh%element_blend ( mesh%nelgrp, 0 ) )
    allocate( mesh%nnodes_blend(2) ); mesh%nnodes_blend = [0,mesh%nnodes]
    allocate( mesh%grpnumel (mesh%nelgrp) )
    allocate( mesh%topology(1) )
    allocate( mesh%topology(1)%a(inpelm,mesh%nelem) )
    allocate( mesh%coor ( mesh%nnodes, mesh%ndim ) )

!   element type
!   element shape number for spectral line elements
    mesh%element(1)%elshape = 101
    mesh%element(1)%numnod = inpelm ! gives the nodal points in one element
    mesh%element(1)%ndim = 1    ! dimension is 1
    mesh%element(1)%p(1,1) = p  ! order of polynomial
    mesh%element(1)%p(1,2) = l  ! layout

!   topology
    mesh%grpnumel(1) = mesh%nelem

!   fill the mesh topology
    do elem = 1, nr
      do i = 1, p+1
        mesh%topology(1)%a(i,elem) = (elem-1)*p+i
      end do
    end do

    if ( l == 0 ) then

!     equidistant distribution on [-1,1]

      x = -1 + 2 * [ ( i - 1, i=1,p+1 ) ] / real(p,kind=dp)

    else if ( l == 1 ) then

!     get the GLL distribution
!     with x the GLL-points (output) and p the polynomial order (input)

      call GLL_points ( x, p )

    end if

    if ( ratio == 0 ) then

!     equidistante mesh

      deltax = 1._dp / nr

      do elem = 1, nr
        do i = 1, p
          node = (elem-1)*p + i
          mesh%coor(node,1)= ( (x(i)+1._dp)/2 + elem-1 ) * deltax
        end do
      end do

      mesh%coor(nn,1) = 1

    else

!     non-equidistante mesh
      call distribute_elements( n, xr(0:n), ratio, factor)

!     all nodes, including refinement
      node=0
      do elem = 1, n
        deltax=(xr(elem)-xr(elem-1))/r
        do i = 1, r ! refinement of the mesh
          ltmp = deltax*(i-1) ! length of previous elements
          do j = 1, p ! nodes in one element
            node=node+1
            mesh%coor(node,1) = xr(elem-1) + ltmp + deltax*((x(j)+1._dp)/2._dp)
          end do
        end do
      end do
      mesh%coor(nn,1) = xr(n)

    end if

    deallocate ( x, xr )

!   points

    mesh%npoints = 2

!   allocate memory for points
    allocate( mesh%points(MAXPOINTS))

    mesh%points(1) = 1
    mesh%points(2) = nn

  end subroutine line1d_ho

end module meshgen_basic_line_m
