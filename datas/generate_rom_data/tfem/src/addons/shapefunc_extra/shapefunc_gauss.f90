
! Copyright (C) 2006-2006 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines shape functions and associated routines with Legrendre Gauss points
! as nodal points

module shapefunc_gauss_m

  use glob_defs_m
  use mesh_m, only: mesh_t

  implicit none


contains


! Line P1

  subroutine shape_line_gauss_P1 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i) with i the point in space
    real(dp), intent(in), dimension(:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j), with i the point in space
!   j the unknown
    real(dp), intent(out), dimension(:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points,
!   which are the Legrendre Gauss points:
!
!      -- 1 ----- 2 --
!
!   The reference coordinates xr are in the region [-1,1]

    real(dp) :: sqrt3

    sqrt3 = sqrt(3._dp)

    phi(:,1) = (1-sqrt3*xr)/2
    phi(:,2) = (1+sqrt3*xr)/2

    if ( present(dphi) ) then
      dphi(:,1) = -sqrt3/2
      dphi(:,2) =  sqrt3/2
    end if

  end subroutine shape_line_gauss_P1


! Line P2

  subroutine shape_line_gauss_P2 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i) with i the point in space
    real(dp), intent(in), dimension(:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j), with i the point in space
!   j the unknown
    real(dp), intent(out), dimension(:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points,
!   which are the Legrendre Gauss points:
!
!    -- 1 --- 2 --- 3 --
!
!   The reference coordinates xr are in the region [-1,1]

    real(dp) :: sqrt15

    sqrt15 = sqrt(15._dp)


    phi(:,1) = (5*xr-sqrt15)*xr/6
    phi(:,2) = 1-5*xr**2/3
    phi(:,3) = (5*xr+sqrt15)*xr/6

    if ( present(dphi) ) then
      dphi(:,1) = (10*xr-sqrt15)/6
      dphi(:,2) = -10*xr/3
      dphi(:,3) = (10*xr+sqrt15)/6
    end if

  end subroutine shape_line_gauss_P2


! Quadrilateral Q1

  subroutine shape_quad_gauss_Q1 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points,
!   which are the Legrendre Gauss points:
!
!       ---------
!     |  3     4  |
!     |           |
!     |  1     2  |
!       ---------
!
!   The reference coordinates xr=(xi,eta) are in the region [-1,1]x[-1,1].


    integer :: i, j, m
    real(dp), dimension(size(xr,1),2) :: phi1, phi2, dphi1, dphi2

    if ( present(dphi) ) then

      call shape_line_gauss_P1 ( xr(:,1), phi1, dphi1 )
      call shape_line_gauss_P1 ( xr(:,2), phi2, dphi2 )

      do i = 1, 2
        do j = 1, 2
          m = i + 2*(j-1)
          phi(:,m)    =  phi1(:,i) *  phi2(:,j)
          dphi(:,m,1) = dphi1(:,i) *  phi2(:,j)
          dphi(:,m,2) =  phi1(:,i) * dphi2(:,j)
        end do
      end do

    else

      call shape_line_gauss_P1 ( xr(:,1), phi1 )
      call shape_line_gauss_P1 ( xr(:,2), phi2 )

      do i = 1, 2
        do j = 1, 2
          m = i + 2*(j-1)
          phi(:,m) = phi1(:,i) * phi2(:,j)
        end do
      end do

    end if

  end subroutine shape_quad_gauss_Q1


! Quadrilateral Q2

  subroutine shape_quad_gauss_Q2 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points,
!   which are the Legrendre Gauss points:
!
!      ---------
!     | 7  8  9 |
!     | 4  5  6 |
!     | 1  2  3 |
!      ---------
!
!   The reference coordinates xr=(xi,eta) are in the region [-1,1]x[-1,1].


    integer :: i, j, m
    real(dp), dimension(size(xr,1),3) :: phi1, phi2, dphi1, dphi2

    if ( present(dphi) ) then

      call shape_line_gauss_P2 ( xr(:,1), phi1, dphi1 )
      call shape_line_gauss_P2 ( xr(:,2), phi2, dphi2 )

      do i = 1, 3
        do j = 1, 3
          m = i + 3 * ( j - 1 )
          phi(:,m)    =  phi1(:,i) *  phi2(:,j)
          dphi(:,m,1) = dphi1(:,i) *  phi2(:,j)
          dphi(:,m,2) =  phi1(:,i) * dphi2(:,j)
        end do
      end do

    else

      call shape_line_gauss_P2 ( xr(:,1), phi1 )
      call shape_line_gauss_P2 ( xr(:,2), phi2 )

      do i = 1, 3
        do j = 1, 3
          m = i + 3 * ( j - 1 )
          phi(:,m) = phi1(:,i) * phi2(:,j)
        end do
      end do

    end if

  end subroutine shape_quad_gauss_Q2


! shape function chooser based on shape of element and number of degrees

  subroutine set_shape_function_gauss ( globalshape, xr, phi, dphi )

!   Global shape of the element:
!     line
!     quadrilateral
!     triangle
!     hexahedron
    character (len=*), intent(in) :: globalshape

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
!   Also ndf=size(phi,2), the number of degrees of freedom, is used to
!   determine the type of shape function.
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   set shape function

    logical :: found
    integer :: ndf


    found = .true.
    ndf = size(phi,2)

!   Select on global shape of the element

    select case ( globalshape )

    case ( 'line' )

      if ( present(dphi) ) then

        select case (ndf)
          case(1)
            phi  = 1
            dphi = 0
          case(2)
            call shape_line_gauss_P1 ( xr(:,1), phi, dphi(:,:,1) )
          case(3)
            call shape_line_gauss_P2 ( xr(:,1), phi, dphi(:,:,1) )
          case default
            found = .false.
        end select

      else

        select case (ndf)
          case(1)
            phi  = 1
          case(2)
            call shape_line_gauss_P1 ( xr(:,1), phi )
          case(3)
            call shape_line_gauss_P2 ( xr(:,1), phi )
          case default
            found = .false.
        end select

      end if

    case ( 'triangle' )

      if ( present(dphi) ) then

        select case (ndf)
          case(1)
            phi  = 1
            dphi = 0
          case default
            found = .false.
        end select

      else

        select case (ndf)
          case(1)
            phi  = 1
          case default
            found = .false.
        end select

      end if

    case ( 'quadrilateral' )

      if ( present(dphi) ) then

        select case (ndf)
          case(1)
            phi  = 1
            dphi = 0
          case(4)
            call shape_quad_gauss_Q1 ( xr, phi, dphi )
          case(9)
            call shape_quad_gauss_Q2 ( xr, phi, dphi )
          case default
            found = .false.
        end select

      else

        select case (ndf)
          case(1)
            phi  = 1
          case(4)
            call shape_quad_gauss_Q1 ( xr, phi )
          case(9)
            call shape_quad_gauss_Q2 ( xr, phi )
          case default
            found = .false.
        end select

      end if

    case ( 'hexahedron' )

      if ( present(dphi) ) then

        select case (ndf)
          case(1)
            phi  = 1
            dphi = 0
          case default
            found = .false.
        end select

      else

        select case (ndf)
          case(1)
            phi  = 1
          case default
            found = .false.
        end select

      end if

    case default

      call errormsg_case_default ( 'set_shape_function_gauss', &
        'globalshape', char_value=globalshape )

    end select

    if ( .not. found ) then
      write(*,'(/a/a,i0,2a/)') 'Error in set_shape_function_gauss:', &
        ' incorrect interpolation, ndf = size(phi,2) = ', ndf, &
        ', globalshape = ', globalshape
      stop
    end if

  end subroutine set_shape_function_gauss

end module shapefunc_gauss_m
