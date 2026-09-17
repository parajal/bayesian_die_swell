
! Copyright (C) 2007-2007 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines shape functions and associated routines for hierarchic shape functions

module shapefunc_hierarchic_m

  use glob_defs_m
  use mesh_m, only: mesh_t
  use shapefunc_m

  implicit none


contains


! Line P2

  subroutine shape_line_hierarchic_P2 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i) with i the point in space
    real(dp), intent(in), dimension(:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j), with i the point in space
!   j the unknown
    real(dp), intent(out), dimension(:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points:
!
!     1 --- 2 --- 3
!
!   The reference coordinates xr are in the region [-1,1]

    phi(:,1) = (1-xr)/2
    phi(:,2) = 1-xr**2
    phi(:,3) = (1+xr)/2

    if ( present(dphi) ) then
      dphi(:,1) = -0.5_dp
      dphi(:,2) = -2*xr
      dphi(:,3) =  0.5_dp
    end if

  end subroutine shape_line_hierarchic_P2


! Quadrilateral Q2

  subroutine shape_quad_hierarchic_Q2 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points:
!
!     7 --- 6 --- 5
!     |           |
!     8     9     4
!     |           |
!     1 --- 2 --- 3
!
!   The reference coordinates xr=(xi,eta) are in the region [-1,1]x[-1,1].


    integer, parameter, dimension(3,3) :: &
      p = reshape ( [1,2,3,8,9,4,7,6,5], [3,3] )
    integer :: i, j
    real(dp), dimension(size(xr,1),3) :: phi1, phi2, dphi1, dphi2

    if ( present(dphi) ) then

      call shape_line_hierarchic_P2 ( xr(:,1), phi1, dphi1 )
      call shape_line_hierarchic_P2 ( xr(:,2), phi2, dphi2 )

      do i = 1, 3
        do j = 1, 3
          phi(:,p(i,j))    =  phi1(:,i) *  phi2(:,j)
          dphi(:,p(i,j),1) = dphi1(:,i) *  phi2(:,j)
          dphi(:,p(i,j),2) =  phi1(:,i) * dphi2(:,j)
        end do
      end do

    else

      call shape_line_hierarchic_P2 ( xr(:,1), phi1 )
      call shape_line_hierarchic_P2 ( xr(:,2), phi2 )

      do i = 1, 3
        do j = 1, 3
          phi(:,p(i,j)) =  phi1(:,i) *  phi2(:,j)
        end do
      end do

    end if

  end subroutine shape_quad_hierarchic_Q2


! Triangle P2

  subroutine shape_triangle_hierarchic_P2 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points:
!
!     5
!     | \          !
!     6   4
!     |     \      !
!     1 --2-- 3
!
!   The reference coordinates xr=(xi,eta) are in lower triangle of
!   the region [0,1]x[0,1]  (xi+eta<=1)
!

    integer :: dm

    real(dp) :: lam(size(xr,1),3), dlam(size(xr,1),3,2)


!   barycentric coordinates lambda1, lambda2, lambda3

    if ( present(dphi) ) then
      call barycentric ( xr, lam, dlam )
    else
      call barycentric ( xr, lam )
    end if

    phi(:,1) = lam(:,1)
    phi(:,2) = 4 * lam(:,1) * lam(:,2)
    phi(:,3) = lam(:,2)
    phi(:,4) = 4 * lam(:,2) * lam(:,3)
    phi(:,5) = lam(:,3)
    phi(:,6) = 4 * lam(:,3) * lam(:,1)

    if ( present(dphi) ) then

      do dm = 1, 2
        dphi(:,1,dm) = dlam(:,1,dm)
        dphi(:,2,dm) = 4 * ( dlam(:,1,dm) * lam(:,2) + &
                             lam(:,1) * dlam(:,2,dm) )
        dphi(:,3,dm) = dlam(:,2,dm)
        dphi(:,4,dm) = 4 * ( dlam(:,2,dm) * lam(:,3) + &
                             lam(:,2) * dlam(:,3,dm) )
        dphi(:,5,dm) = dlam(:,3,dm)
        dphi(:,6,dm) = 4 * ( dlam(:,3,dm) * lam(:,1) + &
                             lam(:,3) * dlam(:,1,dm) )
      end do

    end if

  end subroutine shape_triangle_hierarchic_P2


! Triangle P2+

  subroutine shape_triangle_hierarchic_P2plus ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points:
!
!     5
!     | \           !
!     6   4
!     | 7   \       !
!     1 --2-- 3
!
!   The reference coordinates xr=(xi,eta) are in lower triangle of
!   the region [0,1]x[0,1]  (xi+eta<=1)
!

    integer :: dm

    real(dp) :: lam(size(xr,1),3), dlam(size(xr,1),3,2)
    real(dp) :: bubble(size(xr,1)), dbubble(size(xr,1),2)


!   barycentric coordinates lambda1, lambda2, lambda3

    if ( present(dphi) ) then
      call barycentric ( xr, lam, dlam )
    else
      call barycentric ( xr, lam )
    end if

    bubble = lam(:,1)*lam(:,2)*lam(:,3)

    phi(:,1) = lam(:,1)
    phi(:,2) = 4 * lam(:,1) * lam(:,2)
    phi(:,3) = lam(:,2)
    phi(:,4) = 4 * lam(:,2) * lam(:,3)
    phi(:,5) = lam(:,3)
    phi(:,6) = 4 * lam(:,3) * lam(:,1)
    phi(:,7) = 27*bubble

    if ( present(dphi) ) then

      do dm = 1, 2
        dbubble(:,dm) = dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                        lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                        lam(:,1)*lam(:,2)*dlam(:,3,dm)
      end do

      do dm = 1, 2
        dphi(:,1,dm) = dlam(:,1,dm)
        dphi(:,2,dm) = 4 * ( dlam(:,1,dm) * lam(:,2) + &
                             lam(:,1) * dlam(:,2,dm) )
        dphi(:,3,dm) = dlam(:,2,dm)
        dphi(:,4,dm) = 4 * ( dlam(:,2,dm) * lam(:,3) + &
                             lam(:,2) * dlam(:,3,dm) )
        dphi(:,5,dm) = dlam(:,3,dm)
        dphi(:,6,dm) = 4 * ( dlam(:,3,dm) * lam(:,1) + &
                             lam(:,3) * dlam(:,1,dm) )
      end do

      dphi(:,7,:) = 27*dbubble

    end if

  end subroutine shape_triangle_hierarchic_P2plus


! Hexahedron Q2

  subroutine shape_hexa_hierarchic_Q2 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points:
!
!        25-----26------27
!        /|            /|
!      22 |   23      24|
!      /  16    17   /  18
!    19-----20------21  |
!     | 13|   14    | 15|
!     |   7------8--|---9
!    10  /  11     12  /
!     | 4      5    | 6
!     |/            |/
!     1------2------3
!
!   The reference coordinates xr=(xi,eta,zeta) are
!   in the region [-1,1]x[-1,1]x[-1,1].

    integer :: i, j, k, m
    real(dp), dimension(size(xr,1),3) :: phi1, phi2, phi3, dphi1, dphi2, dphi3


    if ( present(dphi) ) then

      call shape_line_hierarchic_P2 ( xr(:,1), phi1, dphi1 )
      call shape_line_hierarchic_P2 ( xr(:,2), phi2, dphi2 )
      call shape_line_hierarchic_P2 ( xr(:,3), phi3, dphi3 )

      do k = 1, 3
        do j = 1, 3
          do i = 1, 3
            m = i + 3*(j-1) + 9*(k-1)
            phi(:,m)    =  phi1(:,i) *  phi2(:,j) *  phi3(:,k)
            dphi(:,m,1) = dphi1(:,i) *  phi2(:,j) *  phi3(:,k)
            dphi(:,m,2) =  phi1(:,i) * dphi2(:,j) *  phi3(:,k)
            dphi(:,m,3) =  phi1(:,i) *  phi2(:,j) * dphi3(:,k)
          end do
        end do
      end do

    else

      call shape_line_hierarchic_P2 ( xr(:,1), phi1 )
      call shape_line_hierarchic_P2 ( xr(:,2), phi2 )
      call shape_line_hierarchic_P2 ( xr(:,3), phi3 )

      do k = 1, 3
        do j = 1, 3
          do i = 1, 3
            m = i + 3*(j-1) + 9*(k-1)
            phi(:,m) = phi1(:,i) * phi2(:,j) * phi3(:,k)
          end do
        end do
      end do

    end if

  end subroutine shape_hexa_hierarchic_Q2


! Tetrahedron P2

  subroutine shape_tetra_hierarchic_P2 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points:
!
!       10
!      /|\           !
!     / 7  9
!     / |    \       !
!     8 1 -6-- 5
!    / /     /
!    / 2  4
!    // /
!    3
!
!   The reference coordinates xr=(xi,eta,zeta) are in the lower tetrahedron of
!   the region [0,1]x[0,1]x[0,1]  (xi+eta+zeta<=1)
!

    integer :: dm

    real(dp) :: lam(size(xr,1),4), dlam(size(xr,1),4,3)


!   barycentric coordinates lambda1, lambda2, lambda3, lambda4

    if ( present(dphi) ) then
      call barycentric ( xr, lam, dlam )
    else
      call barycentric ( xr, lam )
    end if

    phi(:,1)  = lam(:,1)
    phi(:,2)  = 4 * lam(:,1) * lam(:,2)
    phi(:,3)  = lam(:,2)
    phi(:,4)  = 4 * lam(:,2) * lam(:,3)
    phi(:,5)  = lam(:,3)
    phi(:,6)  = 4 * lam(:,3) * lam(:,1)
    phi(:,7)  = 4 * lam(:,1) * lam(:,4)
    phi(:,8)  = 4 * lam(:,2) * lam(:,4)
    phi(:,9)  = 4 * lam(:,3) * lam(:,4)
    phi(:,10) = lam(:,4)

    if ( present(dphi) ) then

      do dm = 1, 3
        dphi(:,1,dm) = dlam(:,1,dm)
        dphi(:,2,dm) = 4 * ( dlam(:,1,dm) * lam(:,2) + &
                             lam(:,1) * dlam(:,2,dm) )
        dphi(:,3,dm) = dlam(:,2,dm)
        dphi(:,4,dm) = 4 * ( dlam(:,2,dm) * lam(:,3) + &
                             lam(:,2) * dlam(:,3,dm) )
        dphi(:,5,dm) = dlam(:,3,dm)
        dphi(:,6,dm) = 4 * ( dlam(:,3,dm) * lam(:,1) + &
                             lam(:,3) * dlam(:,1,dm) )
        dphi(:,7,dm) = 4 * ( dlam(:,1,dm) * lam(:,4) + &
                             lam(:,1) * dlam(:,4,dm) )
        dphi(:,8,dm) = 4 * ( dlam(:,2,dm) * lam(:,4) + &
                             lam(:,2) * dlam(:,4,dm) )
        dphi(:,9,dm) = 4 * ( dlam(:,3,dm) * lam(:,4) + &
                             lam(:,3) * dlam(:,4,dm) )
        dphi(:,10,dm) = dlam(:,4,dm)
      end do

    end if

  end subroutine shape_tetra_hierarchic_P2


! shape function chooser based on shape of element and number of degrees

  subroutine set_shape_function_hierarchic ( globalshape, xr, phi, dphi )

!   Global shape of the element:
!     line
!     quadrilateral
!     triangle
!     hexahedron
!     tetrahedron
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
          case(3)
            call shape_line_hierarchic_P2 ( xr(:,1), phi, dphi(:,:,1) )
          case default
            found = .false.
        end select

      else

        select case (ndf)
          case(3)
            call shape_line_hierarchic_P2 ( xr(:,1), phi )
          case default
            found = .false.
        end select

      end if

    case ( 'triangle' )

      if ( present(dphi) ) then

        select case (ndf)
          case(6)
            call shape_triangle_hierarchic_P2 ( xr, phi, dphi )
          case(7)
            call shape_triangle_hierarchic_P2plus ( xr, phi, dphi )
          case default
            found = .false.
        end select

      else

        select case (ndf)
          case(6)
            call shape_triangle_hierarchic_P2 ( xr, phi )
          case(7)
            call shape_triangle_hierarchic_P2plus ( xr, phi )
          case default
            found = .false.
        end select

      end if

    case ( 'quadrilateral' )

      if ( present(dphi) ) then

        select case (ndf)
          case(9)
            call shape_quad_hierarchic_Q2 ( xr, phi, dphi )
          case default
            found = .false.
        end select

      else

        select case (ndf)
          case(9)
            call shape_quad_hierarchic_Q2 ( xr, phi )
          case default
            found = .false.
        end select

      end if

    case ( 'hexahedron' )

      if ( present(dphi) ) then

        select case (ndf)
          case(27)
            call shape_hexa_hierarchic_Q2 ( xr, phi, dphi )
          case default
            found = .false.
        end select

      else

        select case (ndf)
          case(27)
            call shape_hexa_hierarchic_Q2 ( xr, phi )
          case default
            found = .false.
        end select

      end if

    case ( 'tetrahedron' )

      if ( present(dphi) ) then

        select case (ndf)
          case(10)
            call shape_tetra_hierarchic_P2 ( xr, phi, dphi )
          case default
            found = .false.
        end select

      else

        select case (ndf)
          case(10)
            call shape_tetra_hierarchic_P2 ( xr, phi )
          case default
            found = .false.
        end select

      end if

    case default

      call errormsg_case_default ( 'set_shape_function_hierarchic', &
        'globalshape', char_value=globalshape )

    end select

    if ( .not. found ) then
      write(*,'(/a/a,i0,2a/)') 'Error in set_shape_function_hierarchic:', &
        ' incorrect interpolation, ndf = size(phi,2) = ', ndf, &
        ', globalshape = ', globalshape
      stop
    end if

  end subroutine set_shape_function_hierarchic

end module shapefunc_hierarchic_m
