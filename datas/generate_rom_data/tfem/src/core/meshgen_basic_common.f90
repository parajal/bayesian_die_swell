
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
! common defs and routines

module meshgen_basic_common_m

  use kind_defs_m
  use limits_m, only: MAXPOINTS, MAXCURVES, MAXOBJECTS, MAXNUMITER, &
                      MAXSURFACES, MAXVOLUMES, MAXBLOCKS, MAXNODESETS, &
                      MAXELEMENTSETS
  use mesh_m, only: mesh_t

  implicit none


! type definition of the options for meshgeneration

  type meshgen_options_t

    integer :: regionshape = 1 ! Shape of the region:
                               !  1: rectangle
                               !  2: quadrilateral
                               !  3: quadrilateral with curved edges
                               !  4: cuboid (=rectangular box)
                               !  5: hexahedron
                               !  6: hexahedron with curved faces (*)
                               !
                               !  (*) not yet available

    integer :: nx = 1      ! number of elements in x-direction
    integer :: ny = 1      ! number of elements in y-direction
    integer :: nz = 1      ! number of elements in z-direction

    integer :: rx = 1       ! refinement factors. Note that the actual number
    integer :: ry = 1       ! of elements is nx*rx, ny*ry, nz*rz
    integer :: rz = 1       ! The difference between
                            !   a) nx=n*r, rx=1
                            !   b) nx=n, rx=r
                            ! is the distribution of the elements along the
                            ! side. The value of nx is used for determining the
                            ! ratio of the element sizes. The value of rx is
                            ! is used to further refine the distribution in an
                            ! equidistant way. This difference is important if
                            ! two domains need to be connected using mortars.

    integer :: elshape = -1 ! shape number of the elements:
                            !   elshape=1 : 2 node line element
                            !   elshape=2 : 3 node line element
                            !   elshape=3 : 3 node triangle
                            !   elshape=4 : 6 node triangle
                            !   elshape=5 : 4 node quadrilateral
                            !   elshape=6 : 9 node quadrilateral
                            !   elshape=7 : 7 node triangle
                            !   elshape=9 : 5 node quadrilateral
                            !   elshape=10: 4 node triangle
                            !   elshape=13: 8 node hexahedron
                            !   elshape=14: 27 node hexahedron
                            !   elshape=101: line high-order element
                            !   elshape=102: quadrilateral high-order element
                            !   elshape=104: triangle high-order element

    integer :: p = 0        ! Polynomial order.
                            ! This parameter can be used in combination with
                            ! elshape=101, 102 and 104.
                            ! The same order p is used in the different
                            ! directions for multidimensional elements

    integer :: l = 1        ! Layout of the nodes:
                            !  0: equidistant
                            !  1: Gauss-Lobatto
                            ! This parameter can be used in combination with
                            ! elshape=101 and 102.
                            ! The same layout l is used in the different
                            ! directions for multidimensional elements
                            ! Note: for elshape=104 l is ignored, since only
                            ! equidistant is available for triangles.

!   Input parameters for a rectangular shape of the region (regionshape=1 or 4):
    real(dp) :: lx = 1.0_dp ! size of region in x-direction
    real(dp) :: ly = 1.0_dp ! size of region in y-direction
    real(dp) :: ox = 0.0_dp ! x-coordinate of lower left corner
    real(dp) :: oy = 0.0_dp ! y-coordinate of lower left corner

!   Input parameters for a rectangular shape of the region (regionshape=4):
    real(dp) :: lz = 1.0_dp ! size of region in z-direction
    real(dp) :: oz = 0.0_dp ! z-coordinate of lower left corner

!   Input parameters for a quadrilateral shape of the region
!   (regionshape=2 or 3).
!   x2d(i,j) are the coordinates of the corners (corner i, direction j).
    real(dp) :: x2d(4,2) = &
      reshape ( [ 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp,    &
                  0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp ], [4,2] )

!   Input parameters for a hexahedral shape of the region
!   (regionshape=5 or 6).
!   x3d(i,j) are the coordinates of the corners (corner i, direction j).
!   NOTE: a cuboid is mapped on to a hexahedron using the corner coordinates.
!   If the four corner nodes of a quadrilateral edge are such that they do not
!   belong to a single plane the edge will not be flat but curved!
    real(dp) :: x3d(8,3) = &
      reshape ( [ 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp,    &
                  0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp,    &
                  0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp,    &
                  0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp,    &
                  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,    &
                  1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp ], [8,3] )

!   Input for regionshape=3: curved boundaries.
!   curved(i)=.true. means curved along curve i
    logical, dimension(4) :: curved = .false.

!   Function number for the decription of the curved boundaries
!   funcnr(i) is number for curve i
    integer, dimension(4) :: funcnr = -1

!   ratio: type of mesh refinement along the curves (2D and 3D) or line (1D).
!          NOTE: that the direction of the curves is as given in the
!          meshgenerator
!     ratio = 0: equidistant mesh
!     ratio = 1: the size of the last element is factor times the first
!     ratio = 2: the size of an element is factor times the previous one
!     ratio = 3: the size of the last element is 1/factor times the first
!     ratio = 4: the size of an element is 1/factor times the previous one
!     ratio = 5: the size of the first element is factor times the size
!                of an element when elements would be equidistant distributed
!     ratio = 6: the size of the last element is factor times the size
!                of an element when elements would be equidistant distributed
!     ratio = 7: the size of the first element is factor times the length
!                of the side
!     ratio = 8: the size of the last element is factor times the length
!                of the side
!     NOTES: - ratio=0--4 are compatible with the meshgenerator of Sepran
!            - ratio=5--8 of Sepran are not available. Split your mesh in two
!              parts if you want that functionality.
!            - ratio=7--8 do not scale with an increase of the number of
!              elements like ratio=1--6. Both factor and the number of elements
!              (nx and/or ny) need to be changed to have a uniform refinement.
    integer, dimension(12)  :: ratio = 0
    real(dp), dimension(12) :: factor = 1._dp

  end type meshgen_options_t


contains


! Helper routine for setting the mesh generator options

  subroutine set_mesh_options ( mesh_options, keep, regionshape, nx, ny, nz, &
    elshape, p, l, lx, ly, lz, ox, oy, oz, x2d, x3d, curved, funcnr, ratio, &
    factor, rx, ry, rz )

!   Options for generating the mesh.
!   See type definition for possibilities and defaults.
    type(meshgen_options_t), intent(inout) :: mesh_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see mesh_options_t
    integer, intent(in), optional :: regionshape, nx, ny, nz, elshape, p, l
    integer, intent(in), optional :: rx, ry, rz, funcnr(4), ratio(:)
    real(dp), intent(in), optional :: lx, ly, lz, ox, oy, oz
    real(dp), intent(in), optional :: x2d(4,2), x3d(8,3), factor(:)
    logical, intent(in), optional :: curved(4)


    type(meshgen_options_t) :: meshopt


    if ( present(keep) ) then
      if ( .not. keep ) mesh_options = meshopt
    else
      mesh_options = meshopt
    end if

    if ( present(regionshape) ) mesh_options%regionshape = regionshape
    if ( present(nx) ) mesh_options%nx = nx
    if ( present(ny) ) mesh_options%ny = ny
    if ( present(nz) ) mesh_options%nz = nz
    if ( present(rx) ) mesh_options%rx = rx
    if ( present(ry) ) mesh_options%ry = ry
    if ( present(rz) ) mesh_options%rz = rz
    if ( present(elshape) ) mesh_options%elshape = elshape
    if ( present(p) ) mesh_options%p = p
    if ( present(l) ) mesh_options%l = l
    if ( present(funcnr) ) mesh_options%funcnr = funcnr
    if ( present(ratio) ) mesh_options%ratio(:size(ratio)) = ratio
    if ( present(factor) ) mesh_options%factor(:size(factor)) = factor
    if ( present(curved) ) mesh_options%curved = curved
    if ( present(lx) ) mesh_options%lx = lx
    if ( present(ly) ) mesh_options%ly = ly
    if ( present(lz) ) mesh_options%lz = lz
    if ( present(ox) ) mesh_options%ox = ox
    if ( present(oy) ) mesh_options%oy = oy
    if ( present(oz) ) mesh_options%oz = oz
    if ( present(x2d) ) mesh_options%x2d = x2d
    if ( present(x3d) ) mesh_options%x3d = x3d

  end subroutine set_mesh_options


! Generate n elements on the interval [0,1] using ratio and factor

  subroutine distribute_elements ( n, x, ratio, factor )

!   number of elements
    integer, intent(in) :: n

!   coordinates of the nodes
    real(dp), intent(out), dimension(0:n) :: x

!   ratio and factor determine the distribution (see meshgen_options_t)
    integer, intent(in) :: ratio
    real(dp), intent(in) :: factor

!   The interval [0:1] is divided into n elements:
!
!     dx_{i+1} = g dx_{i-1}
!                                               1 - g^n
!     1 = (1+g+g^2+g^3+....+g^{n-1}) dx_1   (= -------- dx_1)
!                                               1 - g
!   with fac = 1+g+g^2+g^3+....+g^{n-1} we have
!
!     dx_1 = 1 / fac        (=(1-g)/(1-g^n))
!     dx_n = g^{n-1} dx_1
!

    integer :: i
    real(dp) :: g ! the factor  dx_{i+1}=g dx_i
    real(dp) :: fac, dx_1, dx
    logical :: conv


    conv = .true.

    if ( factor < 0 ) then
      write(*,'(/a,i0/)') 'Error: negative factor: ', factor
      stop
    end if

    if ( n <= 1 ) then
      write(*,'(/a,i0/)') 'Error: number of elements must be at least 2'
      stop
    end if

    select case (ratio)
      case(0)
        g = 1._dp
      case(1)
        g = exp(log(factor)/(n-1))
      case(2)
        g = factor
      case(3)
        g = exp(-log(factor)/(n-1))
      case(4)
        g = 1._dp/factor
      case(5)
        dx_1 = factor / n
        call compute_factor ( n, dx_1, g, conv )
      case(6)
        dx_1 = factor / n
        call compute_factor ( n, dx_1, g, conv )
        g = 1._dp / g
      case(7)
        dx_1 = factor
        call compute_factor ( n, dx_1, g, conv )
      case(8)
        dx_1 = factor
        call compute_factor ( n, dx_1, g, conv )
        g = 1._dp / g
      case default
        write(*,'(/a,i0/)') 'Error: invalid value for ratio: ', ratio
        stop
    end select

    if ( .not. conv ) then
      write(*,'(/a/)') 'Error: no convergence in element distribution'
      stop
    end if

!   generate mesh

    fac = 1

    do i = 1, n - 1
      fac = fac + g ** i
    end do

!   size of first element

    dx_1 = 1._dp / fac

!   generate all elements

    x(0) = 0
    dx = dx_1
    do i = 1, n
      x(i) = x(i-1) + dx
      dx = g * dx
    end do

!   test whether x(n) = 1

    if ( abs(x(n) - 1._dp) > 1.e-10_dp ) then
      write(*,'(/a,i0/)') 'Error: end value x(n) /= 1: ', x(n)
      stop
    else
      x(n) = 1._dp
    end if

  end subroutine distribute_elements


! Compute factor for distribution of elements

  subroutine compute_factor ( n, dx_1, g, conv )

!   number of elements
    integer, intent(in) :: n

!   length of the first interval
    real(dp), intent(in) :: dx_1

!   computed factor between adjacent intervals: dx_{i+1} = g dx_i
    real(dp), intent(out) :: g

!   logical for indicating convergence
    logical, intent(out) :: conv


    integer :: iter, j
    real(dp) :: fac1, fac2, y, dy, dg

    conv = .true.

!   initialize factor

    if ( dx_1 * n < 1._dp ) then
      g = exp(log(1/dx_1)/(n-1))
    else
      g = exp(log(1/dx_1)/(n-1))
      g = 1/g
    end if

!   start loop

    iter = 0

    do

      fac1 = 0
      fac2 = 0
      do j = 0, n - 1
         fac1 = fac1 + g ** j
         fac2 = fac2 + j * g ** ( j - 1 )
      end do

      y  = fac1 * dx_1 - 1._dp   !  find zero of this function
      dy = fac2 * dx_1

      dg = - y / dy    ! Newton-Raphson
      g  = g + dg

      iter = iter + 1

!      print *, iter, y

      if ( iter >= MAXNUMITER ) then
!       no convergence
        conv = .false.
        exit
      else if ( abs(y) < 1e-12_dp ) then
!       convergence
        exit
      end if

    end do

  end subroutine compute_factor

end module meshgen_basic_common_m
