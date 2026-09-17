
! Copyright (C) 2007-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines shape functions and associated routines for macro elements

module shapefunc_macro_m

  use kind_defs_m
  use mesh_m, only: mesh_t

  implicit none


! interface to support two interfaces for shape_line_P1isoP2

  interface shape_line_P1isoP2
    module procedure shape_line_P1isoP2_1, shape_line_P1isoP2_2
  end interface shape_line_P1isoP2

! interface to support two interfaces for shape_line_P1isoPp

  interface shape_line_P1isoPp
    module procedure shape_line_P1isoPp_1, shape_line_P1isoPp_2
  end interface shape_line_P1isoPp

contains


! Line P1isoP2 (one dimension in xr)

  subroutine shape_line_P1isoP2_1 ( xr, phi, dphi )

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
!     1 -- 2 -- 3
!
!   The reference coordinates xr are in the region [-1,1]
!
!   The element consists of 2 linear subelements
!
!     x-1-x-2-x
!

    where ( xr < 0.0_dp )

!     subelement 1

      phi(:,1) = - xr
      phi(:,2) = 1 + xr
      phi(:,3) = 0

    else where

!     subelement 2

      phi(:,1) = 0
      phi(:,2) = 1 - xr
      phi(:,3) = xr

    end where

    if ( present(dphi) ) then

      where ( xr < 0.0_dp )

!       subelement 1

        dphi(:,1) = - 1
        dphi(:,2) = 1
        dphi(:,3) = 0

      else where

!       subelement 2

        dphi(:,1) = 0
        dphi(:,2) = - 1
        dphi(:,3) = 1

      end where

    end if

  end subroutine shape_line_P1isoP2_1


! Line P1isoP2 (two dimensions in xr)

  subroutine shape_line_P1isoP2_2 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,1) with i the point in space
    real(dp), dimension(:,:), intent(in) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j), with i the point in space
!   j the unknown
    real(dp), intent(out), dimension(:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points:
!
!     1 -- 2 -- 3
!
!   The reference coordinates xr are in the region [-1,1]
!
!   The element consists of 2 linear subelements
!
!     x-1-x-2-x
!

    where ( xr(:,1) < 0.0_dp )

!     subelement 1

      phi(:,1) = - xr(:,1)
      phi(:,2) = 1 + xr(:,1)
      phi(:,3) = 0

    else where

!     subelement 2

      phi(:,1) = 0
      phi(:,2) = 1 - xr(:,1)
      phi(:,3) = xr(:,1)

    end where

    if ( present(dphi) ) then

      where ( xr(:,1) < 0.0_dp )

!       subelement 1

        dphi(:,1) = - 1
        dphi(:,2) = 1
        dphi(:,3) = 0

      else where

!       subelement 2

        dphi(:,1) = 0
        dphi(:,2) = - 1
        dphi(:,3) = 1

      end where

    end if

  end subroutine shape_line_P1isoP2_2


! Triangle P1isoP2

  subroutine shape_triangle_P1isoP2 ( xr, phi, dphi )

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
!     |     \       !
!     1 --2-- 3
!
!   The reference coordinates xr=(xi,eta) are in lower triangle of
!   the region [0,1]x[0,1]  (xi+eta<=1)
!
!   The element consists of 4 linear subtriangles
!
!     x
!     |3\           !
!     x-- x
!     |1\4|2\       !
!     x --x-- x
!

    phi = 0

    where ( xr(:,1) + xr(:,2) < 0.5_dp )

!     subelement 1

      phi(:,1) = 1 - 2 * xr(:,1) - 2 * xr(:,2)
      phi(:,2) = 2 * xr(:,1)
      phi(:,6) = 2 * xr(:,2)

    else where ( xr(:,1) > 0.5_dp )

!     subelement 2

      phi(:,2) = 2 - 2 * xr(:,1) - 2 * xr(:,2)
      phi(:,3) = 2 * xr(:,1) - 1
      phi(:,4) = 2 * xr(:,2)

    else where ( xr(:,2) > 0.5_dp )

!     subelement 3

      phi(:,6) = 2 - 2 * xr(:,1) - 2 * xr(:,2)
      phi(:,4) = 2 * xr(:,1)
      phi(:,5) = 2 * xr(:,2) - 1

    else where

!     subelement 4

      phi(:,2) = 1 - 2 * xr(:,2)
      phi(:,4) = 2 * xr(:,1) + 2 * xr(:,2) - 1
      phi(:,6) = 1 - 2 * xr(:,1)

    end where

    if ( present(dphi) ) then

      dphi = 0

      where ( xr(:,1) + xr(:,2) < 0.5_dp )

!       subelement 1

        dphi(:,1,1) = - 2
        dphi(:,1,2) = - 2
        dphi(:,2,1) = 2
        dphi(:,6,2) = 2

      else where ( xr(:,1) > 0.5_dp )

!       subelement 2

        dphi(:,2,1) = - 2
        dphi(:,2,2) = - 2
        dphi(:,3,1) = 2
        dphi(:,4,2) = 2

      else where ( xr(:,2) > 0.5_dp )

!       subelement 3

        dphi(:,6,1) = - 2
        dphi(:,6,2) = - 2
        dphi(:,4,1) = 2
        dphi(:,5,2) = 2

      else where

!       subelement 4

        dphi(:,2,2) = - 2
        dphi(:,4,1) = 2
        dphi(:,4,2) = 2
        dphi(:,6,1) = - 2

      end where

    end if

  end subroutine shape_triangle_P1isoP2


! Quadrilateral Q1-iso-Q2

  subroutine shape_quad_Q1isoQ2 ( xr, phi, dphi )

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
!
!   The element consists of 4 bilinear subquads
!
!     x --- x --- x
!     |  3  |  4  |
!     x --- x --- x
!     |  1  |  2  |
!     x --- x --- x
!


    integer, parameter, dimension(3,3) :: &
      p = reshape ( [1,2,3,8,9,4,7,6,5], [3,3] )
    integer :: i, j
    real(dp), dimension(size(xr,1),3) :: phi1, phi2, dphi1, dphi2

    if ( present(dphi) ) then

      call shape_line_P1isoP2 ( xr(:,1), phi1, dphi1 )
      call shape_line_P1isoP2 ( xr(:,2), phi2, dphi2 )

      do i = 1, 3
        do j = 1, 3
          phi(:,p(i,j))    =  phi1(:,i) *  phi2(:,j)
          dphi(:,p(i,j),1) = dphi1(:,i) *  phi2(:,j)
          dphi(:,p(i,j),2) =  phi1(:,i) * dphi2(:,j)
        end do
      end do

    else

      call shape_line_P1isoP2 ( xr(:,1), phi1 )
      call shape_line_P1isoP2 ( xr(:,2), phi2 )

      do i = 1, 3
        do j = 1, 3
          phi(:,p(i,j)) =  phi1(:,i) *  phi2(:,j)
        end do
      end do

    end if

  end subroutine shape_quad_Q1isoQ2


! Tetrahedron P1isoP2

  subroutine shape_tetra_P1isoP2 ( xr, phi, dphi )

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
!   The element consists of 8 linear subtetrahedrons according to the topology:
!
!   T1=[ 1 2 6 7 ]
!   T2=[ 2 3 4 8 ]
!   T3=[ 4 5 6 9 ]
!   T4=[ 7 8 9 10 ]
!   T5=[ 2 8 9 7 ]
!   T6=[ 2 9 6 7 ]
!   T7=[ 2 9 8 4 ]
!   T8=[ 2 6 9 4 ]
!

    phi = 0

    where ( xr(:,1) + xr(:,2) + xr(:,3) <= 0.5_dp )

!     subelement 1

      phi(:,1) = 1 - 2 * xr(:,1) - 2 * xr(:,2) - 2 * xr(:,3)
      phi(:,2) = 2 * xr(:,1)
      phi(:,6) = 2 * xr(:,2)
      phi(:,7) = 2 * xr(:,3)

    else where ( xr(:,1) >= 0.5_dp )

!     subelement 2

      phi(:,2) = 2 - 2 * xr(:,1) - 2 * xr(:,2) - 2 * xr(:,3)
      phi(:,3) = 2 * xr(:,1) - 1
      phi(:,4) = 2 * xr(:,2)
      phi(:,8) = 2 * xr(:,3)

    else where ( xr(:,2) >= 0.5_dp )

!     subelement 3

      phi(:,4) = 2 * xr(:,1)
      phi(:,5) = 2 * xr(:,2) - 1
      phi(:,6) = 2 - 2 * xr(:,1) - 2 * xr(:,2) - 2 * xr(:,3)
      phi(:,9) = 2 * xr(:,3)

    else where ( xr(:,3) >= 0.5_dp )

!     subelement 4

      phi(:,7) = 2 - 2 * xr(:,1) - 2 * xr(:,2) - 2 * xr(:,3)
      phi(:,8) = 2 * xr(:,1)
      phi(:,9) = 2 * xr(:,2)
      phi(:,10) = 2 * xr(:,3) - 1

    else where ( xr(:,1) + xr(:,2) <= 0.5_dp .and. &
                 xr(:,1) + xr(:,3) >= 0.5_dp )

!     subelement 5

      phi(:,2) = 1 - 2 * xr(:,3)
      phi(:,8) = 2 * xr(:,1) + 2 * xr(:,3) - 1
      phi(:,9) = 2 * xr(:,2)
      phi(:,7) = 1 - 2 * xr(:,1) - 2 * xr(:,2)

    else where ( xr(:,1) + xr(:,2) <= 0.5_dp .and. &
                 xr(:,1) + xr(:,3) <= 0.5_dp )

!     subelement 6

      phi(:,2) = 2 * xr(:,1)
      phi(:,9) = 2 * xr(:,1) + 2 * xr(:,2) + 2 * xr(:,3) - 1
      phi(:,6) = 1 - 2 * xr(:,1) - 2 * xr(:,3)
      phi(:,7) = 1 - 2 * xr(:,1) - 2 * xr(:,2)

    else where ( xr(:,1) + xr(:,2) >= 0.5_dp .and. &
                 xr(:,1) + xr(:,3) >= 0.5_dp )

!     subelement 7

      phi(:,2) = 2 - 2 * xr(:,1) - 2 * xr(:,2) - 2 * xr(:,3)
      phi(:,9) = 1 - 2 * xr(:,1)
      phi(:,8) = 2 * xr(:,1) + 2 * xr(:,3) - 1
      phi(:,4) = 2 * xr(:,1) + 2 * xr(:,2) - 1

    else where

!     subelement 8

      phi(:,2) = 1 - 2 * xr(:,2)
      phi(:,6) = 1 - 2 * xr(:,1) - 2 * xr(:,3)
      phi(:,9) = 2 * xr(:,3)
      phi(:,4) = 2 * xr(:,1) + 2 * xr(:,2) - 1

    end where

    if ( present(dphi) ) then

      dphi = 0

      where ( xr(:,1) + xr(:,2) + xr(:,3) <= 0.5_dp )

!       subelement 1

        dphi(:,1,1) = - 2
        dphi(:,1,2) = - 2
        dphi(:,1,3) = - 2
        dphi(:,2,1) = 2
        dphi(:,6,2) = 2
        dphi(:,7,3) = 2

      else where ( xr(:,1) >= 0.5_dp )

!       subelement 2

        dphi(:,2,1) = - 2
        dphi(:,2,2) = - 2
        dphi(:,2,3) = - 2
        dphi(:,3,1) = 2
        dphi(:,4,2) = 2
        dphi(:,8,3) = 2

      else where ( xr(:,2) >= 0.5_dp )

!       subelement 3

        dphi(:,4,1) = 2
        dphi(:,5,2) = 2
        dphi(:,6,1) = - 2
        dphi(:,6,2) = - 2
        dphi(:,6,3) = - 2
        dphi(:,9,3) = 2

      else where ( xr(:,3) >= 0.5_dp )

!       subelement 4

        dphi(:,7,1) = - 2
        dphi(:,7,2) = - 2
        dphi(:,7,3) = - 2
        dphi(:,8,1) = 2
        dphi(:,9,2) = 2
        dphi(:,10,3) = 2

      else where ( xr(:,1) + xr(:,2) <= 0.5_dp .and. &
                   xr(:,1) + xr(:,3) >= 0.5_dp )

!       subelement 5

        dphi(:,2,3) = - 2
        dphi(:,8,1) = 2
        dphi(:,8,3) = 2
        dphi(:,9,2) = 2
        dphi(:,7,1) = - 2
        dphi(:,7,2) = - 2

      else where ( xr(:,1) + xr(:,2) <= 0.5_dp .and. &
                   xr(:,1) + xr(:,3) <= 0.5_dp )

!       subelement 6

        dphi(:,2,1) = 2
        dphi(:,9,1) = 2
        dphi(:,9,2) = 2
        dphi(:,9,3) = 2
        dphi(:,6,1) = - 2
        dphi(:,6,3) = - 2
        dphi(:,7,1) = - 2
        dphi(:,7,2) = - 2

      else where ( xr(:,1) + xr(:,2) >= 0.5_dp .and. &
                   xr(:,1) + xr(:,3) >= 0.5_dp )

!       subelement 7

        dphi(:,2,1) = - 2
        dphi(:,2,2) = - 2
        dphi(:,2,3) = - 2
        dphi(:,9,1) = - 2
        dphi(:,8,1) = 2
        dphi(:,8,3) = 2
        dphi(:,4,1) = 2
        dphi(:,4,2) = 2

      else where

!       subelement 8

        dphi(:,2,2) = - 2
        dphi(:,6,1) = - 2
        dphi(:,6,3) = - 2
        dphi(:,9,3) = 2
        dphi(:,4,1) = 2
        dphi(:,4,2) = 2

      end where

    end if

  end subroutine shape_tetra_P1isoP2


! Hexahedron Q1-iso-Q2

  subroutine shape_hexa_Q1isoQ2 ( xr, phi, dphi )

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
!          zeta
!            ^ _ eta
!            | /|
!            |/
!             ---> xi
!
!   The reference coordinates xr=(xi,eta,zeta) are
!   in the region [-1,1]x[-1,1]x[-1,1].
!
!   The element consists of 8 trilinear subhexes
!

    integer :: i, j, k
    integer, parameter, dimension(3,3,3) :: &
      p = reshape ( [(i,i=1,27)], [3,3,3] )
    real(dp), dimension(size(xr,1),3) :: phi1, phi2, phi3, dphi1, dphi2, dphi3

    if ( present(dphi) ) then

      call shape_line_P1isoP2 ( xr(:,1), phi1, dphi1 )
      call shape_line_P1isoP2 ( xr(:,2), phi2, dphi2 )
      call shape_line_P1isoP2 ( xr(:,3), phi3, dphi3 )

      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            phi(:,p(i,j,k))    =  phi1(:,i) *  phi2(:,j) *  phi3(:,k)
            dphi(:,p(i,j,k),1) = dphi1(:,i) *  phi2(:,j) *  phi3(:,k)
            dphi(:,p(i,j,k),2) =  phi1(:,i) * dphi2(:,j) *  phi3(:,k)
            dphi(:,p(i,j,k),3) =  phi1(:,i) *  phi2(:,j) * dphi3(:,k)
          end do
        end do
      end do

    else

      call shape_line_P1isoP2 ( xr(:,1), phi1 )
      call shape_line_P1isoP2 ( xr(:,2), phi2 )
      call shape_line_P1isoP2 ( xr(:,3), phi3 )

      do i = 1, 3
        do j = 1, 3
          do k = 1, 3
            phi(:,p(i,j,k)) = phi1(:,i) * phi2(:,j) * phi3(:,k)
          end do
        end do
      end do

    end if

  end subroutine shape_hexa_Q1isoQ2


! Line P1isoPp (one dimension in xr)

  subroutine shape_line_P1isoPp_1 ( xr, p, phi, dphi, iv )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i) with i the point in space
    real(dp), intent(in), dimension(:) :: xr

!   polynomial order
    integer, intent(in) :: p

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j), with i the point in space
!   j the unknown
    real(dp), intent(out), dimension(:,:), optional :: dphi

!   Interval number. Size must be equal to size of xr
    integer, intent(out), dimension(:), optional :: iv

!   The unknowns of the element are positioned at the following nodal points:
!
!     1 -- 2 -- 3 .... p -- p+1
!
!   The reference coordinates xr are in the region [-1,1]
!
!   The element consists of p linear subelements

    integer :: m
    integer, dimension(size(xr)) :: i


!   interval number

    i = int( p*(xr+1)/2 ) + 1

!   avoid boundary nodes to fall outside range

    where ( i > p )
      i = p
    end where
    where ( i < 1 )
      i = 1
    end where

    phi = 0

    do m = 1, size(xr)
      phi(m,i(m)) = i(m) - p*(1+xr(m))/2
      phi(m,i(m)+1) = 1 - i(m) + p*(1+xr(m))/2
    end do

    if ( present(dphi) ) then
      dphi = 0
      do m = 1, size(xr)
        dphi(m,i(m)) = - real(p,dp)/2
        dphi(m,i(m)+1) = real(p,dp)/2
      end do
    end if

    if ( present(iv) ) iv = i

  end subroutine shape_line_P1isoPp_1


! Line P1isoPp (two dimensions in xr)

  subroutine shape_line_P1isoPp_2 ( xr, p, phi, dphi, iv )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,1) with i the point in space
    real(dp), dimension(:,:), intent(in) :: xr

!   polynomial order
    integer, intent(in) :: p

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,1), with i the point in space
!   j the unknown
    real(dp), intent(out), dimension(:,:,:), optional :: dphi

!   Interval number. Size must be equal to size of xr
    integer, intent(out), dimension(:), optional :: iv

!   The unknowns of the element are positioned at the following nodal points:
!
!     1 -- 2 -- 3 .... p -- p+1
!
!   The reference coordinates xr are in the region [-1,1]
!
!   The element consists of p linear subelements

    integer :: m
    integer, dimension(size(xr,1)) :: i


!   interval number

    i = int( p*(xr(:,1)+1)/2 ) + 1

!   avoid boundary nodes to fall outside range

    where ( i > p )
      i = p
    end where
    where ( i < 1 )
      i = 1
    end where

    phi = 0

    do m = 1, size(xr,1)
      phi(m,i(m)) = i(m) - p*(1+xr(m,1))/2
      phi(m,i(m)+1) = 1 - i(m) + p*(1+xr(m,1))/2
    end do

    if ( present(dphi) ) then
      dphi = 0
      do m = 1, size(xr,1)
        dphi(m,i(m),1) = - real(p,dp)/2
        dphi(m,i(m)+1,1) = real(p,dp)/2
      end do
    end if

    if ( present(iv) ) iv = i

  end subroutine shape_line_P1isoPp_2


! Quadrilateral Q1-iso-Qp

  subroutine shape_quad_Q1isoQp ( xr, p, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   polynomial order
    integer, intent(in) :: p

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points:
!
! (p+1)p+1 -- ...... --- (p+1)**2
!     |                  |
!     .                  .
!     .                  .
!   p+2 --- p+3 .... --- 2p+2
!     |                  |
!     1 --- 2 .... p --- p+1
!
!  eta
!   ^
!   |
!    --> xi
!
!   The reference coordinates xr=(xi,eta) are in the region [-1,1]x[-1,1].

!   The element consists of p x p bilinear subquads
!
!     x --- x --- x  ...  x --- x --- x
!     |     |     |  ...  |p^2-1| p^2 |
!     x --- x --- x  ...  x --- x --- x
!     .     .     .       .     .     .
!     .     .     .       .     .     .
!     x --- x --- x  ...  x --- x --- x
!     | p+1 | p+2 |  ...  |2p+1 | 2p  |
!     x --- x --- x  ...  x --- x --- x
!     |  1  |  2  |  ...  | p-1 |  p  |
!     x --- x --- x  ...  x --- x --- x


    integer :: i, j, m, ii, jj
    integer, dimension(size(xr,1)) :: iv1, iv2
    integer, dimension(p+1,p+1) :: pos
    real(dp), dimension(size(xr,1),p+1) :: phi1, phi2, dphi1, dphi2

    pos = reshape ( [(i,i=1,(p+1)**2)], [p+1,p+1] )

    if ( present(dphi) ) then

      call shape_line_P1isoPp ( xr(:,1), p, phi1, dphi1, iv1 )
      call shape_line_P1isoPp ( xr(:,2), p, phi2, dphi2, iv2 )

      phi = 0
      dphi = 0

      do m = 1, size(xr,1)
        ii = iv1(m); jj = iv2(m)
        do i = ii, ii+1
          do j = jj, jj+1
            phi(m,pos(i,j))    =  phi1(m,i) *  phi2(m,j)
            dphi(m,pos(i,j),1) = dphi1(m,i) *  phi2(m,j)
            dphi(m,pos(i,j),2) =  phi1(m,i) * dphi2(m,j)
          end do
        end do
      end do

    else

      call shape_line_P1isoPp ( xr(:,1), p, phi1, iv=iv1 )
      call shape_line_P1isoPp ( xr(:,2), p, phi2, iv=iv2 )

      phi = 0

      do m = 1, size(xr,1)
        ii = iv1(m); jj = iv2(m)
        do i = ii, ii+1
          do j = jj, jj+1
            phi(m,pos(i,j)) = phi1(m,i) * phi2(m,j)
          end do
        end do
      end do

    end if

  end subroutine shape_quad_Q1isoQp


! Hexahedron Q1-iso-Qp

  subroutine shape_hexa_Q1isoQp ( xr, p, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   polynomial order
    integer, intent(in) :: p

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi

!   The unknowns of the element are positioned at the following nodal points:
!
!    (p+1)^3-p----------(p+1)^3
!        /|            /|
!       / |           / |
!      /  |          /  |
! p(p+1)^2+1--------(p+1)((p+1)p+1)
!     |   |         |   |
!     | (p+1)p+1----|--(p+1)^2
!     |  /          |  /
!     | /           | /
!     |/            |/
!     1--- ....  ---p+1
!
!          zeta
!            ^ _ eta
!            | /|
!            |/
!             ---> xi
!
!   The reference coordinates xr=(xi,eta,zeta) are
!   in the region [-1,1]x[-1,1]x[-1,1].
!
!   The element consists of p^3 trilinear subhexes

    integer :: i, j, k, m, ii, jj, kk
    integer, dimension(size(xr,1)) :: iv1, iv2, iv3
    integer, dimension(p+1,p+1,p+1) :: pos
    real(dp), dimension(size(xr,1),p+1) :: phi1, phi2, phi3, dphi1, dphi2, dphi3

    pos = reshape ( [(i,i=1,(p+1)**3)], [p+1,p+1,p+1] )

    if ( present(dphi) ) then

      call shape_line_P1isoPp ( xr(:,1), p, phi1, dphi1, iv1 )
      call shape_line_P1isoPp ( xr(:,2), p, phi2, dphi2, iv2 )
      call shape_line_P1isoPp ( xr(:,3), p, phi3, dphi3, iv3 )

      phi = 0
      dphi = 0

      do m = 1, size(xr,1)
        ii = iv1(m); jj = iv2(m); kk = iv3(m)
        do i = ii, ii+1
          do j = jj, jj+1
            do k = kk, kk+1
              phi(m,pos(i,j,k))    =  phi1(m,i) *  phi2(m,j) *  phi3(m,k)
              dphi(m,pos(i,j,k),1) = dphi1(m,i) *  phi2(m,j) *  phi3(m,k)
              dphi(m,pos(i,j,k),2) =  phi1(m,i) * dphi2(m,j) *  phi3(m,k)
              dphi(m,pos(i,j,k),3) =  phi1(m,i) *  phi2(m,j) * dphi3(m,k)
            end do
          end do
        end do
      end do

    else

      call shape_line_P1isoPp ( xr(:,1), p, phi1, iv=iv1 )
      call shape_line_P1isoPp ( xr(:,2), p, phi2, iv=iv2 )
      call shape_line_P1isoPp ( xr(:,3), p, phi3, iv=iv3 )

      phi = 0

      do m = 1, size(xr,1)
        ii = iv1(m); jj = iv2(m)
        do i = ii, ii+1
          do j = jj, jj+1
            do k = kk, kk+1
              phi(m,pos(i,j,k)) = phi1(m,i) * phi2(m,j) * phi3(m,k)
            end do
          end do
        end do
      end do

    end if

  end subroutine shape_hexa_Q1isoQp


end module shapefunc_macro_m
