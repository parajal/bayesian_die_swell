
! Copyright (C) 2004-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines standard shape functions

module shapefunc_standard_m

  use kind_defs_m

  implicit none


! interface to support two interfaces for shape_line_P1

  interface shape_line_P1
    module procedure shape_line_P1_1, shape_line_P1_2
  end interface shape_line_P1

! interface to support two interfaces for shape_line_P2

  interface shape_line_P2
    module procedure shape_line_P2_1, shape_line_P2_2
  end interface shape_line_P2


contains


! Line P1 (one dimension in xr))

  subroutine shape_line_P1_1 ( xr, phi, dphi )

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
!     1 --------- 2
!
!   The reference coordinates xr are in the region [-1,1]


    phi(:,1) = (1-xr)/2
    phi(:,2) = (1+xr)/2

    if ( present(dphi) ) then
      dphi(:,1) = -0.5_dp
      dphi(:,2) =  0.5_dp
    end if

  end subroutine shape_line_P1_1


! Line P1 (two dimensions in xr)

  subroutine shape_line_P1_2 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,1) with i the point in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,1), with i the point in space
!   j the unknown
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points:
!
!     1 --------- 2
!
!   The reference coordinates xr are in the region [-1,1]


    phi(:,1) = (1-xr(:,1))/2
    phi(:,2) = (1+xr(:,1))/2

    if ( present(dphi) ) then
      dphi(:,1,1) = -0.5_dp
      dphi(:,2,1) =  0.5_dp
    end if

  end subroutine shape_line_P1_2


! Line P2 (one dimension in xr)

  subroutine shape_line_P2_1 ( xr, phi, dphi )

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


    phi(:,1) = -(1-xr)*xr/2
    phi(:,2) = 1-xr**2
    phi(:,3) = (1+xr)*xr/2

    if ( present(dphi) ) then
      dphi(:,1) = -(1-2*xr)/2
      dphi(:,2) = -2*xr
      dphi(:,3) = (1+2*xr)/2
    end if

  end subroutine shape_line_P2_1



! Line P2 (two dimensions in xr)

  subroutine shape_line_P2_2 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,1) with i the point in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,1), with i the point in space
!   j the unknown
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points:
!
!     1 --- 2 --- 3
!
!   The reference coordinates xr are in the region [-1,1]


    phi(:,1) = -(1-xr(:,1))*xr(:,1)/2
    phi(:,2) = 1-xr(:,1)**2
    phi(:,3) = (1+xr(:,1))*xr(:,1)/2

    if ( present(dphi) ) then
      dphi(:,1,1) = -(1-2*xr(:,1))/2
      dphi(:,2,1) = -2*xr(:,1)
      dphi(:,3,1) = (1+2*xr(:,1))/2
    end if

  end subroutine shape_line_P2_2


! Quadrilateral P1

  subroutine shape_quad_P1 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following points:
!
!     -------------
!     |           |     x center point value
!     |     |     |     -- coefficient of xi  (slope)
!     |     x--   |     |  coefficient of eta  (slope)
!     |           |
!     -------------
!
!  eta
!   ^
!   |
!    --> xi
!
!   Shape functions: 1, xi, eta
!
!   The reference coordinates xr=(xi,eta) are in the region [-1,1]x[-1,1].


    phi(:,1) = 1
    phi(:,2) = xr(:,1)
    phi(:,3) = xr(:,2)

    if ( present(dphi) ) then

      dphi(:,1,:) = 0
      dphi(:,2,1) = 1
      dphi(:,2,2) = 0
      dphi(:,3,1) = 0
      dphi(:,3,2) = 1

    end if

  end subroutine shape_quad_P1


! Quadrilateral Q1

  subroutine shape_quad_Q1 ( xr, phi, dphi )

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
!     4 --------- 3
!     |           |
!     |           |
!     |           |
!     1 --------- 2
!
!  eta
!   ^
!   |
!    --> xi
!
!   The reference coordinates xr=(xi,eta) are in the region [-1,1]x[-1,1].


    integer, parameter, dimension(2,2) :: &
      p = reshape ( [1,2,4,3], [2,2] )
    integer :: i, j
    real(dp), dimension(size(xr,1),2) :: phi1, phi2, dphi1, dphi2

    if ( present(dphi) ) then

      call shape_line_P1 ( xr(:,1), phi1, dphi1 )
      call shape_line_P1 ( xr(:,2), phi2, dphi2 )

      do i = 1, 2
        do j = 1, 2
          phi(:,p(i,j))    =  phi1(:,i) *  phi2(:,j)
          dphi(:,p(i,j),1) = dphi1(:,i) *  phi2(:,j)
          dphi(:,p(i,j),2) =  phi1(:,i) * dphi2(:,j)
        end do
      end do

    else

      call shape_line_P1 ( xr(:,1), phi1 )
      call shape_line_P1 ( xr(:,2), phi2 )

      do i = 1, 2
        do j = 1, 2
          phi(:,p(i,j)) =  phi1(:,i) *  phi2(:,j)
        end do
      end do

    end if

  end subroutine shape_quad_Q1


! Quadrilateral Q1plus (bubble)

  subroutine shape_quad_Q1plus ( xr, phi, dphi )

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
!     4 --------- 3
!     |           |
!     |     5     |
!     |           |
!     1 --------- 2
!
!  eta
!   ^
!   |
!    --> xi
!
!   The reference coordinates xr=(xi,eta) are in the region [-1,1]x[-1,1].


    integer :: dm
    real(dp) :: phiQ1(size(xr,1),4)
    real(dp) :: dphiQ1(size(xr,1),4,2)
    real(dp) :: bubble(size(xr,1)), dbubble(size(xr,1),2)


    if ( present(dphi) ) then

      call shape_quad_Q1 ( xr, phiQ1, dphiQ1 )

      bubble = phiQ1(:,1)*phiQ1(:,3)

      phi(:,1) = phiQ1(:,1) - 4*bubble
      phi(:,2) = phiQ1(:,2) - 4*bubble
      phi(:,3) = phiQ1(:,3) - 4*bubble
      phi(:,4) = phiQ1(:,4) - 4*bubble
      phi(:,5) = 16*bubble

      do dm = 1, 2
        dbubble(:,dm) = dphiQ1(:,1,dm)*phiQ1(:,3) + &
                        phiQ1(:,1)*dphiQ1(:,3,dm)
      end do

      dphi(:,1,:) = dphiQ1(:,1,:) - 4*dbubble
      dphi(:,2,:) = dphiQ1(:,2,:) - 4*dbubble
      dphi(:,3,:) = dphiQ1(:,3,:) - 4*dbubble
      dphi(:,4,:) = dphiQ1(:,4,:) - 4*dbubble
      dphi(:,5,:) = 16*dbubble

    else

      call shape_quad_Q1 ( xr, phiQ1 )

      bubble = phiQ1(:,1)*phiQ1(:,3)

      phi(:,1) = phiQ1(:,1) - 4*bubble
      phi(:,2) = phiQ1(:,2) - 4*bubble
      phi(:,3) = phiQ1(:,3) - 4*bubble
      phi(:,4) = phiQ1(:,4) - 4*bubble
      phi(:,5) = 16*bubble

    end if

  end subroutine shape_quad_Q1plus


! Quadrilateral Q2

  subroutine shape_quad_Q2 ( xr, phi, dphi )

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
!  eta
!   ^
!   |
!    --> xi
!
!   The reference coordinates xr=(xi,eta) are in the region [-1,1]x[-1,1].


    integer, parameter, dimension(3,3) :: &
      p = reshape ( [1,2,3,8,9,4,7,6,5], [3,3] )
    integer :: i, j
    real(dp), dimension(size(xr,1),3) :: phi1, phi2, dphi1, dphi2

    if ( present(dphi) ) then

      call shape_line_P2 ( xr(:,1), phi1, dphi1 )
      call shape_line_P2 ( xr(:,2), phi2, dphi2 )

      do i = 1, 3
        do j = 1, 3
          phi(:,p(i,j))    =  phi1(:,i) *  phi2(:,j)
          dphi(:,p(i,j),1) = dphi1(:,i) *  phi2(:,j)
          dphi(:,p(i,j),2) =  phi1(:,i) * dphi2(:,j)
        end do
      end do

    else

      call shape_line_P2 ( xr(:,1), phi1 )
      call shape_line_P2 ( xr(:,2), phi2 )

      do i = 1, 3
        do j = 1, 3
          phi(:,p(i,j)) =  phi1(:,i) *  phi2(:,j)
        end do
      end do

    end if

  end subroutine shape_quad_Q2


! Quadrilateral serendipity (quadratic)

  subroutine shape_quad_serendipity2 ( xr, phi, dphi )

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
!     8           4
!     |           |
!     1 --- 2 --- 3
!
!  eta
!   ^
!   |
!    --> xi
!
!   The reference coordinates xr=(xi,eta) are in the region [-1,1]x[-1,1].
!
!   Based on the formulas:  (with xi0=xi_i xi, eta0=eta_i eta)
!
!      corner nodes:
!               phi = (1+xi0)(1+eta0)(xi0+eta0-1)/4
!      midside nodes xi_i=0
!               phi = (1+eta0)(1-xi^2)/2
!      midside nodes eta_i=0
!               phi = (1+xi0)(1-eta^2)/2
!

    integer :: i
    integer, parameter :: cor(4) = [ 1, 3, 5, 7 ] ! corner nodes
    integer, parameter :: mi1(2) = [ 2, 6 ] ! midside nodes xi=0
    integer, parameter :: mi2(2) = [ 8, 4 ] ! midside nodes eta=0
    integer, dimension(8) :: xii, etai ! xi, eta in the nodes
    real(dp), dimension(size(xr,1)) :: xi, eta
    real(dp), dimension(size(xr,1),8) :: xi0, eta0

    xii  = [ -1,  0,  1, 1, 1, 0, -1, -1 ]
    etai = [ -1, -1, -1, 0, 1, 1,  1,  0 ]

    xi  = xr(:,1)
    eta = xr(:,2)

    do i = 1, 8
      xi0(:,i)  = xii(i) * xi
      eta0(:,i) = etai(i) * eta
    end do

    phi(:,cor) = ( 1 + xi0(:,cor) ) * ( 1 + eta0(:,cor) ) &
                     * ( xi0(:,cor) + eta0(:,cor) - 1 ) / 4
    phi(:,mi1) = 1 + eta0(:,mi1)
    phi(:,mi2) = 1 + xi0(:,mi2)
    do i = 1, 2
      phi(:,mi1(i)) = phi(:,mi1(i)) * ( 1 - xi**2 ) / 2
      phi(:,mi2(i)) = phi(:,mi2(i)) * ( 1 - eta**2 ) / 2
    end do

    if ( present(dphi) ) then

      dphi(:,cor,1) = ( 1 + eta0(:,cor) ) * ( 2 * xi0(:,cor) + eta0(:,cor) ) / 4
      dphi(:,cor,2) = ( 1 + xi0(:,cor)) * ( xi0(:,cor) + 2 * eta0(:,cor) ) / 4
      do i = 1, 4
        dphi(:,cor(i),1) = dphi(:,cor(i),1) * xii(cor(i))
        dphi(:,cor(i),2) = dphi(:,cor(i),2) * etai(cor(i))
      end do
      dphi(:,mi1,1) = 1 + eta0(:,mi1)
      dphi(:,mi2,2) = 1 + xi0(:,mi2)
      do i = 1, 2
        dphi(:,mi1(i),1) = - dphi(:,mi1(i),1) * xi
        dphi(:,mi1(i),2) = etai(mi1(i)) * ( 1 - xi**2 ) / 2
        dphi(:,mi2(i),1) = xii(mi2(i)) * ( 1 - eta**2 ) / 2
        dphi(:,mi2(i),2) = - dphi(:,mi2(i),2) * eta
      end do

    end if

  end subroutine shape_quad_serendipity2


! barycentric coordinates

  subroutine barycentric ( xr, lambda, dlambda )

!   Reference coordinates where lambda and dlambda must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function lambda(i,j), with i the point in space and j the number
    real(dp), intent(out), dimension(:,:) :: lambda

!   derivative of the shape function: dlambda(i,j,k), with i the point in space
!   j the number and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dlambda

!   Triangles:
!   ---------
!
!   the reference region of a triangle is
!
!      1 |\               eta \in [0,1]
!      ^ |  \             xi  \in [0,1]
!  eta | |    \           xi + eta <= 1
!        |      \      !
!      0 ---------
!        0  xi -> 1
!
!   The reference coordinates xr=(xi,eta) are in lower triangle of
!   the region [0,1]x[0,1]  (xi+eta<=1)
!
!   The three barycentric coordinates are
!
!   lambda1=1-xi-eta
!   lambda2=xi
!   lambda3=eta
!
!
!   Tetrahedrons:
!   ---------
!
!   the reference region of a tetrahedron is
!
!             1 |\               eta  \in [0,1]
!             ^ |  \             xi   \in [0,1]
!        zeta | |    \           zeta \in [0,1]
!               |      \         xi + eta + zeta <= 1
!             0 ---------
!         xi  / 0 eta -> 1
!            /
!         1 /
!
!
!   The reference coordinates xr=(xi,eta,zeta) are in lower tetrahedron of
!   the region [0,1]x[0,1]x[0,1]  (xi+eta+zeta<=1)
!
!   The four barycentric coordinates are
!
!   lambda1=1-xi-eta-zeta
!   lambda2=xi
!   lambda3=eta
!   lambda4=zeta
!

    integer :: refdim

    refdim = size(xr,2)

    if ( refdim == 2 ) then

!     triangle

      lambda(:,1) = 1 - xr(:,1) - xr(:,2)
      lambda(:,2) = xr(:,1)
      lambda(:,3) = xr(:,2)

      if ( present(dlambda) ) then

        dlambda(:,1,1) = -1
        dlambda(:,1,2) = -1
        dlambda(:,2,1) = 1
        dlambda(:,2,2) = 0
        dlambda(:,3,1) = 0
        dlambda(:,3,2) = 1

      end if

    else if ( refdim == 3 ) then

!     tetrahedron

      lambda(:,1) = 1 - xr(:,1) - xr(:,2) - xr(:,3)
      lambda(:,2) = xr(:,1)
      lambda(:,3) = xr(:,2)
      lambda(:,4) = xr(:,3)

      if ( present(dlambda) ) then

        dlambda(:,1,1) = -1
        dlambda(:,1,2) = -1
        dlambda(:,1,3) = -1
        dlambda(:,2,1) = 1
        dlambda(:,2,2) = 0
        dlambda(:,2,3) = 0
        dlambda(:,3,1) = 0
        dlambda(:,3,2) = 1
        dlambda(:,3,3) = 0
        dlambda(:,4,1) = 0
        dlambda(:,4,2) = 0
        dlambda(:,4,3) = 1

      end if

    end if

  end subroutine barycentric


! Triangle P1

  subroutine shape_triangle_P1 ( xr, phi, dphi )

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
!     3
!     | \          !
!     |   \        !
!     |     \      !
!     1 ----- 2
!
!   The reference coordinates xr=(xi,eta) are in lower triangle of
!   the region [0,1]x[0,1]  (xi+eta<=1)
!


    if ( present(dphi) ) then
      call barycentric ( xr, phi, dphi )
    else
      call barycentric ( xr, phi )
    end if

  end subroutine shape_triangle_P1


! Triangle P1+

  subroutine shape_triangle_P1plus ( xr, phi, dphi )

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
!     3
!     | \           !
!     |   \         !
!     | 4   \       !
!     1 ----- 2
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

    phi(:,1) = lam(:,1) - 9*bubble
    phi(:,2) = lam(:,2) - 9*bubble
    phi(:,3) = lam(:,3) - 9*bubble
    phi(:,4) = 27*bubble

    if ( present(dphi) ) then

      do dm = 1, 2
        dbubble(:,dm) = dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                        lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                        lam(:,1)*lam(:,2)*dlam(:,3,dm)
      end do

      dphi(:,1,:) = dlam(:,1,:) - 9*dbubble
      dphi(:,2,:) = dlam(:,2,:) - 9*dbubble
      dphi(:,3,:) = dlam(:,3,:) - 9*dbubble
      dphi(:,4,:) = 27*dbubble

    end if

  end subroutine shape_triangle_P1plus


! Triangle P2

  subroutine shape_triangle_P2 ( xr, phi, dphi )

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

    phi(:,1) = lam(:,1) * ( 2*lam(:,1) - 1 )
    phi(:,2) = 4 * lam(:,1) * lam(:,2)
    phi(:,3) = lam(:,2) * ( 2*lam(:,2) - 1 )
    phi(:,4) = 4 * lam(:,2) * lam(:,3)
    phi(:,5) = lam(:,3) * ( 2*lam(:,3) - 1 )
    phi(:,6) = 4 * lam(:,3) * lam(:,1)

    if ( present(dphi) ) then

      do dm = 1, 2
        dphi(:,1,dm) = ( 4 * lam(:,1) - 1 ) * dlam(:,1,dm)
        dphi(:,2,dm) = 4 * ( dlam(:,1,dm) * lam(:,2) + &
                             lam(:,1) * dlam(:,2,dm) )
        dphi(:,3,dm) = ( 4 * lam(:,2) - 1 ) * dlam(:,2,dm)
        dphi(:,4,dm) = 4 * ( dlam(:,2,dm) * lam(:,3) + &
                             lam(:,2) * dlam(:,3,dm) )
        dphi(:,5,dm) = ( 4 * lam(:,3) - 1 ) * dlam(:,3,dm)
        dphi(:,6,dm) = 4 * ( dlam(:,3,dm) * lam(:,1) + &
                             lam(:,3) * dlam(:,1,dm) )
      end do

    end if

  end subroutine shape_triangle_P2


! Triangle P2+

  subroutine shape_triangle_P2plus ( xr, phi, dphi )

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
!     | 7   \      !
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

    phi(:,1) = lam(:,1) * ( 2*lam(:,1) - 1 ) + 3*bubble
    phi(:,2) = 4 * lam(:,1) * lam(:,2) - 12*bubble
    phi(:,3) = lam(:,2) * ( 2*lam(:,2) - 1 ) + 3*bubble
    phi(:,4) = 4 * lam(:,2) * lam(:,3) - 12*bubble
    phi(:,5) = lam(:,3) * ( 2*lam(:,3) - 1 ) + 3*bubble
    phi(:,6) = 4 * lam(:,3) * lam(:,1) - 12*bubble
    phi(:,7) = 27*bubble

    if ( present(dphi) ) then

      do dm = 1, 2
        dbubble(:,dm) = dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                        lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                        lam(:,1)*lam(:,2)*dlam(:,3,dm)
      end do

      do dm = 1, 2
        dphi(:,1,dm) = ( 4 * lam(:,1) - 1 ) * dlam(:,1,dm) + 3*dbubble(:,dm)
        dphi(:,2,dm) = 4 * ( dlam(:,1,dm) * lam(:,2) + &
                             lam(:,1) * dlam(:,2,dm) ) - 12*dbubble(:,dm)
        dphi(:,3,dm) = ( 4 * lam(:,2) - 1 ) * dlam(:,2,dm) + 3*dbubble(:,dm)
        dphi(:,4,dm) = 4 * ( dlam(:,2,dm) * lam(:,3) + &
                             lam(:,2) * dlam(:,3,dm) ) - 12*dbubble(:,dm)
        dphi(:,5,dm) = ( 4 * lam(:,3) - 1 ) * dlam(:,3,dm) + 3*dbubble(:,dm)
        dphi(:,6,dm) = 4 * ( dlam(:,3,dm) * lam(:,1) + &
                             lam(:,3) * dlam(:,1,dm) ) - 12*dbubble(:,dm)
      end do

      dphi(:,7,:) = 27*dbubble

    end if

  end subroutine shape_triangle_P2plus


! Hexahedron P1

  subroutine shape_hexa_P1 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following points:
!
!        -------------
!       /|          /|
!      / |         / |
!     -------|/----  |    x center point value
!     |  |   x--  |  |   -- coefficient of xi   (slope)
!     |  ---------|---    / coefficient of eta  (slope)
!     | /         | /     | coefficient of zeta (slope)
!     |/          |/
!     -------------
!
!          zeta
!            ^ _ eta
!            | /|
!            |/
!             ---> xi
!
!   Shape functions: 1, xi, eta, zeta
!
!   The reference coordinates xr=(xi,eta,zeta) are
!   in the region [-1,1]x[-1,1]x[-1,1].


    phi(:,1) = 1
    phi(:,2) = xr(:,1)
    phi(:,3) = xr(:,2)
    phi(:,4) = xr(:,3)

    if ( present(dphi) ) then

      dphi(:,1,:) = 0
      dphi(:,2,1) = 1
      dphi(:,2,2) = 0
      dphi(:,2,3) = 0
      dphi(:,3,1) = 0
      dphi(:,3,2) = 1
      dphi(:,3,3) = 0
      dphi(:,4,1) = 0
      dphi(:,4,2) = 0
      dphi(:,4,3) = 1

    end if

  end subroutine shape_hexa_P1


! Hexahedron Q1

  subroutine shape_hexa_Q1 ( xr, phi, dphi )

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
!          zeta
!            ^ _ eta
!            | /|
!            |/
!             ---> xi
!
!   The reference coordinates xr=(xi,eta,zeta) are
!   in the region [-1,1]x[-1,1]x[-1,1].



    integer, parameter, dimension(2,2,2) :: &
      p = reshape ( [1,2,4,3,5,6,8,7], [2,2,2] )
    integer :: i, j, k
    real(dp), dimension(size(xr,1),2) :: phi1, phi2, phi3, dphi1, dphi2, dphi3

    if ( present(dphi) ) then

      call shape_line_P1 ( xr(:,1), phi1, dphi1 )
      call shape_line_P1 ( xr(:,2), phi2, dphi2 )
      call shape_line_P1 ( xr(:,3), phi3, dphi3 )

      do i = 1, 2
        do j = 1, 2
          do k = 1, 2
            phi(:,p(i,j,k))    =  phi1(:,i) *  phi2(:,j) *  phi3(:,k)
            dphi(:,p(i,j,k),1) = dphi1(:,i) *  phi2(:,j) *  phi3(:,k)
            dphi(:,p(i,j,k),2) =  phi1(:,i) * dphi2(:,j) *  phi3(:,k)
            dphi(:,p(i,j,k),3) =  phi1(:,i) *  phi2(:,j) * dphi3(:,k)
          end do
        end do
      end do

    else

      call shape_line_P1 ( xr(:,1), phi1 )
      call shape_line_P1 ( xr(:,2), phi2 )
      call shape_line_P1 ( xr(:,3), phi3 )

      do i = 1, 2
        do j = 1, 2
          do k = 1, 2
            phi(:,p(i,j,k)) = phi1(:,i) * phi2(:,j) * phi3(:,k)
          end do
        end do
      end do

    end if

  end subroutine shape_hexa_Q1


! Hexahedron Q1 (regular numbering)

  subroutine shape_hexa_Q1_reg ( xr, phi, dphi )

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
!        7-----------8
!       /|          /|
!      / |         / |
!     5-----------6  |
!     |  |        |  |
!     |  3--------|--4
!     | /         | /
!     |/          |/
!     1-----------2
!
!          zeta
!            ^ _ eta
!            | /|
!            |/
!             ---> xi
!
!   The reference coordinates xr=(xi,eta,zeta) are
!   in the region [-1,1]x[-1,1]x[-1,1].


    integer :: i, j, k, m
    real(dp), dimension(size(xr,1),2) :: phi1, phi2, phi3, dphi1, dphi2, dphi3

    if ( present(dphi) ) then

      call shape_line_P1 ( xr(:,1), phi1, dphi1 )
      call shape_line_P1 ( xr(:,2), phi2, dphi2 )
      call shape_line_P1 ( xr(:,3), phi3, dphi3 )

      do k = 1, 2
        do j = 1, 2
          do i = 1, 2
            m = i + 2*(j-1) + 4*(k-1)
            phi(:,m)    =  phi1(:,i) *  phi2(:,j) *  phi3(:,k)
            dphi(:,m,1) = dphi1(:,i) *  phi2(:,j) *  phi3(:,k)
            dphi(:,m,2) =  phi1(:,i) * dphi2(:,j) *  phi3(:,k)
            dphi(:,m,3) =  phi1(:,i) *  phi2(:,j) * dphi3(:,k)
          end do
        end do
      end do

    else

      call shape_line_P1 ( xr(:,1), phi1 )
      call shape_line_P1 ( xr(:,2), phi2 )
      call shape_line_P1 ( xr(:,3), phi3 )

      do k = 1, 2
        do j = 1, 2
          do i = 1, 2
            m = i + 2*(j-1) + 4*(k-1)
            phi(:,m) = phi1(:,i) * phi2(:,j) * phi3(:,k)
          end do
        end do
      end do

    end if

  end subroutine shape_hexa_Q1_reg


! Hexahedron Q2

  subroutine shape_hexa_Q2 ( xr, phi, dphi )

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

    integer :: i, j, k, m
    real(dp), dimension(size(xr,1),3) :: phi1, phi2, phi3, dphi1, dphi2, dphi3


    if ( present(dphi) ) then

      call shape_line_P2 ( xr(:,1), phi1, dphi1 )
      call shape_line_P2 ( xr(:,2), phi2, dphi2 )
      call shape_line_P2 ( xr(:,3), phi3, dphi3 )

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

      call shape_line_P2 ( xr(:,1), phi1 )
      call shape_line_P2 ( xr(:,2), phi2 )
      call shape_line_P2 ( xr(:,3), phi3 )

      do k = 1, 3
        do j = 1, 3
          do i = 1, 3
            m = i + 3*(j-1) + 9*(k-1)
            phi(:,m) = phi1(:,i) * phi2(:,j) * phi3(:,k)
          end do
        end do
      end do

    end if

  end subroutine shape_hexa_Q2


! Hexahedron serendipity (quadratic)

  subroutine shape_hexa_serendipity2 ( xr, phi, dphi )

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
!        18-----19------20
!        /|            /|
!      16 |           17|
!      /  11         /  12
!    13-----14------15  |
!     |   |         |   |
!     |   6------7--|---8
!     9  /         10  /
!     | 4           | 5
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
!   Based on the formulas:  (with xi0=xi_i xi, eta0=eta_i eta)
!
!      corner nodes:
!               phi = (1+xi0)(1+eta0)(1+zeta0)(xi0+eta0+zeta0-2)/8
!      midside nodes xi_i=0
!               phi = (1+eta0)(1+zeta0)(1-xi^2)/4
!      midside nodes eta_i=0
!               phi = (1+xi0)(1+zeta0)(1-eta^2)/4
!      midside nodes zeta_i=0
!               phi = (1+xi0)(1+eta0)(1-zeta^2)/4
!

    integer :: i
    integer, parameter :: cor(8) = [ 1, 3, 6, 8, 13, 15, 18, 20 ] ! corner nodes
    integer, parameter :: mi1(4) = [ 2, 7, 14, 19 ] ! midside nodes xi=0
    integer, parameter :: mi2(4) = [ 4, 5, 16, 17 ] ! midside nodes eta=0
    integer, parameter :: mi3(4) = [ 9, 10, 11, 12 ] ! midside nodes zeta=0
    integer, dimension(20) :: xii, etai, zetai ! xi, eta in the nodes
    real(dp), dimension(size(xr,1)) :: xi, eta, zeta
    real(dp), dimension(size(xr,1),20) :: xi0, eta0, zeta0

    xii   = [ -1,  0,  1, -1,  1, -1,  0,  1, -1,  1, &
              -1,  1, -1,  0,  1, -1,  1, -1,  0,  1 ]
    etai  = [ -1, -1, -1,  0,  0, 1,  1,  1, -1, -1, &
               1,  1, -1, -1, -1, 0,  0,  1,  1,  1  ]
    zetai = [ -1, -1, -1, -1, -1, -1, -1, -1, 0,  0, &
               0,  0,  1,  1,  1,  1,  1,  1, 1,  1 ]

    xi   = xr(:,1)
    eta  = xr(:,2)
    zeta = xr(:,3)

    do i = 1, 20
      xi0(:,i)   = xii(i)*xi
      eta0(:,i)  = etai(i)*eta
      zeta0(:,i) = zetai(i)*zeta
    end do

    phi(:,cor) = ( 1 + xi0(:,cor) ) * ( 1+eta0(:,cor) ) * ( 1 + zeta0(:,cor) ) &
                    * ( xi0(:,cor) + eta0(:,cor) + zeta0(:,cor) - 2 ) / 8
    phi(:,mi1) = ( 1 + eta0(:,mi1) ) * ( 1 + zeta0(:,mi1) )
    phi(:,mi2) = ( 1 + xi0(:,mi2) ) * ( 1 + zeta0(:,mi2) )
    phi(:,mi3) = ( 1 + xi0(:,mi3)) * ( 1 + eta0(:,mi3) )
    do i = 1, 4
      phi(:,mi1(i)) = phi(:,mi1(i)) * ( 1 - xi**2 ) / 4
      phi(:,mi2(i)) = phi(:,mi2(i)) * ( 1 - eta**2 ) / 4
      phi(:,mi3(i)) = phi(:,mi3(i)) * ( 1 - zeta**2 ) / 4
    end do

    if ( present(dphi) ) then

      dphi(:,cor,1) = ( 1 + eta0(:,cor) ) * ( 1 + zeta0(:,cor) ) &
                       * ( 2 * xi0(:,cor) + eta0(:,cor) + zeta0(:,cor) - 1 ) / 8
      dphi(:,cor,2) = ( 1 + xi0(:,cor) ) * ( 1 + zeta0(:,cor) ) &
                       * ( xi0(:,cor) + 2 * eta0(:,cor) + zeta0(:,cor) - 1 ) / 8
      dphi(:,cor,3) = ( 1 + xi0(:,cor) ) * ( 1 + eta0(:,cor) ) &
                       * ( xi0(:,cor) + eta0(:,cor) + 2 * zeta0(:,cor) - 1 ) / 8
      do i = 1, 8
        dphi(:,cor(i),1) = dphi(:,cor(i),1) * xii(cor(i))
        dphi(:,cor(i),2) = dphi(:,cor(i),2) * etai(cor(i))
        dphi(:,cor(i),3) = dphi(:,cor(i),3) * zetai(cor(i))
      end do
      dphi(:,mi1,1) = ( 1 + eta0(:,mi1) ) * ( 1 + zeta0(:,mi1) )
      dphi(:,mi1,2) = 1 + zeta0(:,mi1)
      dphi(:,mi1,3) = 1 + eta0(:,mi1)
      dphi(:,mi2,1) = 1 + zeta0(:,mi2)
      dphi(:,mi2,2) = ( 1 + xi0(:,mi2) ) * ( 1 + zeta0(:,mi2) )
      dphi(:,mi2,3) = 1 + xi0(:,mi2)
      dphi(:,mi3,1) = 1 + eta0(:,mi3)
      dphi(:,mi3,2) = 1 + xi0(:,mi3)
      dphi(:,mi3,3) = ( 1 + xi0(:,mi3) ) * ( 1 + eta0(:,mi3) )
      do i = 1, 4
        dphi(:,mi1(i),1) = - dphi(:,mi1(i),1) * xi / 2
        dphi(:,mi1(i),2) = dphi(:,mi1(i),2) * etai(mi1(i)) * ( 1 - xi**2 ) / 4
        dphi(:,mi1(i),3) = dphi(:,mi1(i),3) * zetai(mi1(i)) * ( 1 - xi**2 ) / 4
        dphi(:,mi2(i),1) = dphi(:,mi2(i),1) * xii(mi2(i)) * ( 1 - eta**2 ) / 4
        dphi(:,mi2(i),2) = - dphi(:,mi2(i),2) * eta / 2
        dphi(:,mi2(i),3) = dphi(:,mi2(i),3) * zetai(mi2(i)) * ( 1 - eta**2 ) / 4
        dphi(:,mi3(i),1) = dphi(:,mi3(i),1) * xii(mi3(i)) * ( 1 - zeta**2 ) / 4
        dphi(:,mi3(i),2) = dphi(:,mi3(i),2) * etai(mi3(i)) * ( 1 - zeta**2 ) / 4
        dphi(:,mi3(i),3) = - dphi(:,mi3(i),3) * zeta / 2
      end do

    end if

  end subroutine shape_hexa_serendipity2


! Tetrahedron P1

  subroutine shape_tetra_P1 ( xr, phi, dphi )

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
!       4
!      /|\           !
!     / |  \         !
!     / |    \       !
!     / 1 ---- 3
!    / /     /
!    / /  /
!    // /
!    2
!
!   The reference coordinates xr=(xi,eta,zeta) are in the lower tetrahedron of
!   the region [0,1]x[0,1]x[0,1]  (xi+eta+zeta<=1)
!


    if ( present(dphi) ) then
      call barycentric ( xr, phi, dphi )
    else
      call barycentric ( xr, phi )
    end if

  end subroutine shape_tetra_P1


! Tetrahedron P2

  subroutine shape_tetra_P2 ( xr, phi, dphi )

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
!      /|\         !
!     / 7  9
!     / |    \     !
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

    phi(:,1)  = lam(:,1) * ( 2*lam(:,1) - 1 )
    phi(:,2)  = 4 * lam(:,1) * lam(:,2)
    phi(:,3)  = lam(:,2) * ( 2*lam(:,2) - 1 )
    phi(:,4)  = 4 * lam(:,2) * lam(:,3)
    phi(:,5)  = lam(:,3) * ( 2*lam(:,3) - 1 )
    phi(:,6)  = 4 * lam(:,3) * lam(:,1)
    phi(:,7)  = 4 * lam(:,1) * lam(:,4)
    phi(:,8)  = 4 * lam(:,2) * lam(:,4)
    phi(:,9)  = 4 * lam(:,3) * lam(:,4)
    phi(:,10) = lam(:,4) * ( 2*lam(:,4) - 1 )

    if ( present(dphi) ) then

      do dm = 1, 3
        dphi(:,1,dm) = ( 4 * lam(:,1) - 1 ) * dlam(:,1,dm)
        dphi(:,2,dm) = 4 * ( dlam(:,1,dm) * lam(:,2) + &
                             lam(:,1) * dlam(:,2,dm) )
        dphi(:,3,dm) = ( 4 * lam(:,2) - 1 ) * dlam(:,2,dm)
        dphi(:,4,dm) = 4 * ( dlam(:,2,dm) * lam(:,3) + &
                             lam(:,2) * dlam(:,3,dm) )
        dphi(:,5,dm) = ( 4 * lam(:,3) - 1 ) * dlam(:,3,dm)
        dphi(:,6,dm) = 4 * ( dlam(:,3,dm) * lam(:,1) + &
                             lam(:,3) * dlam(:,1,dm) )
        dphi(:,7,dm) = 4 * ( dlam(:,1,dm) * lam(:,4) + &
                             lam(:,1) * dlam(:,4,dm) )
        dphi(:,8,dm) = 4 * ( dlam(:,2,dm) * lam(:,4) + &
                             lam(:,2) * dlam(:,4,dm) )
        dphi(:,9,dm) = 4 * ( dlam(:,3,dm) * lam(:,4) + &
                             lam(:,3) * dlam(:,4,dm) )
        dphi(:,10,dm) = ( 4 * lam(:,4) - 1 ) * dlam(:,4,dm)
      end do

    end if

  end subroutine shape_tetra_P2


! Tetrahedron P2+(14) (14 nodes having four face bubbles)

  subroutine shape_tetra_P2plus14 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points
!   in the vertices and mid-side edges (like a 10-node tet):
!
!       10
!      /|\         !
!     / 7  9
!     / |    \     !
!     8 1 -6-- 5
!    / /     /
!    / 2  4
!    // /
!    3
!
!   and the nodal points at the center of the four faces:
!
!       |
!      /| \        !
!     / |14\       !
!     / |    \     !
!    /12 13----    NOTE: node 13 is on the face xi+eta+zeta=1
!    /  /    /
!    / / 11/
!    //  /
!     /
!
!   The reference coordinates xr=(xi,eta,zeta) are in the lower tetrahedron of
!   the region [0,1]x[0,1]x[0,1]  (xi+eta+zeta<=1)
!

    integer :: dm

    real(dp) :: lam(size(xr,1),4), dlam(size(xr,1),4,3)
    real(dp) :: b11(size(xr,1)), db11(size(xr,1),3)
    real(dp) :: b12(size(xr,1)), db12(size(xr,1),3)
    real(dp) :: b13(size(xr,1)), db13(size(xr,1),3)
    real(dp) :: b14(size(xr,1)), db14(size(xr,1),3)


!   barycentric coordinates lambda1, lambda2, lambda3, lambda4

    if ( present(dphi) ) then
      call barycentric ( xr, lam, dlam )
    else
      call barycentric ( xr, lam )
    end if

!   face bubbles

    b11 = lam(:,1)*lam(:,2)*lam(:,3)
    b12 = lam(:,1)*lam(:,2)*lam(:,4)
    b13 = lam(:,2)*lam(:,3)*lam(:,4)
    b14 = lam(:,3)*lam(:,1)*lam(:,4)

!   shape functions

    phi(:,1)  = lam(:,1) * ( 2*lam(:,1) - 1 ) + 3*b11 + 3*b12 + 3*b14
    phi(:,2)  = 4 * lam(:,1) * lam(:,2) - 12*b11 - 12*b12
    phi(:,3)  = lam(:,2) * ( 2*lam(:,2) - 1 ) + 3*b11 + 3*b12 + 3*b13
    phi(:,4)  = 4 * lam(:,2) * lam(:,3) - 12*b11 - 12*b13
    phi(:,5)  = lam(:,3) * ( 2*lam(:,3) - 1 ) + 3*b11 + 3*b13 + 3*b14
    phi(:,6)  = 4 * lam(:,3) * lam(:,1) - 12*b11 - 12*b14
    phi(:,7)  = 4 * lam(:,1) * lam(:,4) - 12*b12 - 12*b14
    phi(:,8)  = 4 * lam(:,2) * lam(:,4) - 12*b12 - 12*b13
    phi(:,9)  = 4 * lam(:,3) * lam(:,4) - 12*b13 - 12*b14
    phi(:,10) = lam(:,4) * ( 2*lam(:,4) - 1 ) + 3*b12 + 3*b13 + 3*b14
    phi(:,11) = 27*b11
    phi(:,12) = 27*b12
    phi(:,13) = 27*b13
    phi(:,14) = 27*b14

    if ( present(dphi) ) then

      do dm = 1, 3
        db11(:,dm) = dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                     lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                     lam(:,1)*lam(:,2)*dlam(:,3,dm)
        db12(:,dm) = dlam(:,1,dm)*lam(:,2)*lam(:,4) + &
                     lam(:,1)*dlam(:,2,dm)*lam(:,4) + &
                     lam(:,1)*lam(:,2)*dlam(:,4,dm)
        db13(:,dm) = dlam(:,2,dm)*lam(:,3)*lam(:,4) + &
                     lam(:,2)*dlam(:,3,dm)*lam(:,4) + &
                     lam(:,2)*lam(:,3)*dlam(:,4,dm)
        db14(:,dm) = dlam(:,3,dm)*lam(:,1)*lam(:,4) + &
                     lam(:,3)*dlam(:,1,dm)*lam(:,4) + &
                     lam(:,3)*lam(:,1)*dlam(:,4,dm)
      end do

      do dm = 1, 3
        dphi(:,1,dm) = ( 4 * lam(:,1) - 1 ) * dlam(:,1,dm) + &
                                 3*db11(:,dm) + 3*db12(:,dm) + 3*db14(:,dm)
        dphi(:,2,dm) = 4 * ( dlam(:,1,dm) * lam(:,2) + &
                             lam(:,1) * dlam(:,2,dm) ) &
                                     - 12*db11(:,dm) - 12*db12(:,dm)
        dphi(:,3,dm) = ( 4 * lam(:,2) - 1 ) * dlam(:,2,dm) + &
                                 3*db11(:,dm) + 3*db12(:,dm) + 3*db13(:,dm)
        dphi(:,4,dm) = 4 * ( dlam(:,2,dm) * lam(:,3) + &
                             lam(:,2) * dlam(:,3,dm) ) &
                                     - 12*db11(:,dm) - 12*db13(:,dm)
        dphi(:,5,dm) = ( 4 * lam(:,3) - 1 ) * dlam(:,3,dm) + &
                                 3*db11(:,dm) + 3*db13(:,dm) + 3*db14(:,dm)
        dphi(:,6,dm) = 4 * ( dlam(:,3,dm) * lam(:,1) + &
                             lam(:,3) * dlam(:,1,dm) ) &
                                     - 12*db11(:,dm) - 12*db14(:,dm)
        dphi(:,7,dm) = 4 * ( dlam(:,1,dm) * lam(:,4) + &
                             lam(:,1) * dlam(:,4,dm) ) &
                                     - 12*db12(:,dm) - 12*db14(:,dm)
        dphi(:,8,dm) = 4 * ( dlam(:,2,dm) * lam(:,4) + &
                             lam(:,2) * dlam(:,4,dm) ) &
                                     - 12*db12(:,dm) - 12*db13(:,dm)
        dphi(:,9,dm) = 4 * ( dlam(:,3,dm) * lam(:,4) + &
                             lam(:,3) * dlam(:,4,dm) ) &
                                     - 12*db13(:,dm) - 12*db14(:,dm)
        dphi(:,10,dm) = ( 4 * lam(:,4) - 1 ) * dlam(:,4,dm) + &
                                  3*db12(:,dm) + 3*db13(:,dm) + 3*db14(:,dm)
      end do

      dphi(:,11,:) = 27*db11
      dphi(:,12,:) = 27*db12
      dphi(:,13,:) = 27*db13
      dphi(:,14,:) = 27*db14

    end if

  end subroutine shape_tetra_P2plus14


! Tetrahedron P2+(15) (15 nodes having four face bubbles and one volume bubble)

  subroutine shape_tetra_P2plus15 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following nodal points
!   in the vertices and mid-side edges (like a 10-node tet):
!
!       10
!      /|\          !
!     / 7  9
!     / |    \      !
!     8 1 -6-- 5
!    / /     /
!    / 2  4
!    // /
!    3
!
!   and the nodal points at the center of the four faces:
!
!       |
!      /| \         !
!     / |14\        !
!     / |    \      !
!    /12 13----    NOTE: node 13 is on the face xi+eta+zeta=1
!    /  /    /
!    / / 11/
!    //  /
!     /
!
!   and a nodal point at the center of the volume
!
!       |
!      /| \        !
!     / |  \       !
!     / |    \     !
!    /  15 ----
!    /  /    /
!    / /   /
!    //  /
!     /
!
!   The reference coordinates xr=(xi,eta,zeta) are in the lower tetrahedron of
!   the region [0,1]x[0,1]x[0,1]  (xi+eta+zeta<=1)
!

    integer :: dm

    real(dp) :: lam(size(xr,1),4), dlam(size(xr,1),4,3)
    real(dp) :: b11(size(xr,1)), db11(size(xr,1),3)
    real(dp) :: b12(size(xr,1)), db12(size(xr,1),3)
    real(dp) :: b13(size(xr,1)), db13(size(xr,1),3)
    real(dp) :: b14(size(xr,1)), db14(size(xr,1),3)
    real(dp) :: b15(size(xr,1)), db15(size(xr,1),3)


!   barycentric coordinates lambda1, lambda2, lambda3, lambda4

    if ( present(dphi) ) then
      call barycentric ( xr, lam, dlam )
    else
      call barycentric ( xr, lam )
    end if

!   face bubbles and volume bubble

    b11 = lam(:,1)*lam(:,2)*lam(:,3)
    b12 = lam(:,1)*lam(:,2)*lam(:,4)
    b13 = lam(:,2)*lam(:,3)*lam(:,4)
    b14 = lam(:,3)*lam(:,1)*lam(:,4)
    b15 = lam(:,1)*lam(:,2)*lam(:,3)*lam(:,4)

!   shape functions

    phi(:,1)  = lam(:,1) * ( 2*lam(:,1) - 1 ) + 3*b11 + 3*b12 + 3*b14 - 4*b15
    phi(:,2)  = 4 * lam(:,1) * lam(:,2) - 12*b11 - 12*b12 + 32*b15
    phi(:,3)  = lam(:,2) * ( 2*lam(:,2) - 1 ) + 3*b11 + 3*b12 + 3*b13 - 4*b15
    phi(:,4)  = 4 * lam(:,2) * lam(:,3) - 12*b11 - 12*b13 + 32*b15
    phi(:,5)  = lam(:,3) * ( 2*lam(:,3) - 1 ) + 3*b11 + 3*b13 + 3*b14 - 4*b15
    phi(:,6)  = 4 * lam(:,3) * lam(:,1) - 12*b11 - 12*b14 + 32*b15
    phi(:,7)  = 4 * lam(:,1) * lam(:,4) - 12*b12 - 12*b14 + 32*b15
    phi(:,8)  = 4 * lam(:,2) * lam(:,4) - 12*b12 - 12*b13 + 32*b15
    phi(:,9)  = 4 * lam(:,3) * lam(:,4) - 12*b13 - 12*b14 + 32*b15
    phi(:,10) = lam(:,4) * ( 2*lam(:,4) - 1 ) + 3*b12 + 3*b13 + 3*b14 - 4*b15
    phi(:,11) = 27*b11 - 108*b15
    phi(:,12) = 27*b12 - 108*b15
    phi(:,13) = 27*b13 - 108*b15
    phi(:,14) = 27*b14 - 108*b15
    phi(:,15) = 256*b15

    if ( present(dphi) ) then

      do dm = 1, 3
        db11(:,dm) = dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                     lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                     lam(:,1)*lam(:,2)*dlam(:,3,dm)
        db12(:,dm) = dlam(:,1,dm)*lam(:,2)*lam(:,4) + &
                     lam(:,1)*dlam(:,2,dm)*lam(:,4) + &
                     lam(:,1)*lam(:,2)*dlam(:,4,dm)
        db13(:,dm) = dlam(:,2,dm)*lam(:,3)*lam(:,4) + &
                     lam(:,2)*dlam(:,3,dm)*lam(:,4) + &
                     lam(:,2)*lam(:,3)*dlam(:,4,dm)
        db14(:,dm) = dlam(:,3,dm)*lam(:,1)*lam(:,4) + &
                     lam(:,3)*dlam(:,1,dm)*lam(:,4) + &
                     lam(:,3)*lam(:,1)*dlam(:,4,dm)
        db15(:,dm) = dlam(:,1,dm)*lam(:,2)*lam(:,3)*lam(:,4) + &
                     lam(:,1)*dlam(:,2,dm)*lam(:,3)*lam(:,4) + &
                     lam(:,1)*lam(:,2)*dlam(:,3,dm)*lam(:,4) + &
                     lam(:,1)*lam(:,2)*lam(:,3)*dlam(:,4,dm)
      end do

      do dm = 1, 3
        dphi(:,1,dm) = ( 4 * lam(:,1) - 1 ) * dlam(:,1,dm) + &
               3*db11(:,dm) + 3*db12(:,dm) + 3*db14(:,dm) - 4*db15(:,dm)
        dphi(:,2,dm) = 4 * ( dlam(:,1,dm) * lam(:,2) + &
                             lam(:,1) * dlam(:,2,dm) ) &
               - 12*db11(:,dm) - 12*db12(:,dm) + 32*db15(:,dm)
        dphi(:,3,dm) = ( 4 * lam(:,2) - 1 ) * dlam(:,2,dm) + &
               3*db11(:,dm) + 3*db12(:,dm) + 3*db13(:,dm) - 4*db15(:,dm)
        dphi(:,4,dm) = 4 * ( dlam(:,2,dm) * lam(:,3) + &
                             lam(:,2) * dlam(:,3,dm) ) &
               - 12*db11(:,dm) - 12*db13(:,dm) + 32*db15(:,dm)
        dphi(:,5,dm) = ( 4 * lam(:,3) - 1 ) * dlam(:,3,dm) + &
               3*db11(:,dm) + 3*db13(:,dm) + 3*db14(:,dm) - 4*db15(:,dm)
        dphi(:,6,dm) = 4 * ( dlam(:,3,dm) * lam(:,1) + &
                             lam(:,3) * dlam(:,1,dm) ) &
               - 12*db11(:,dm) - 12*db14(:,dm) + 32*db15(:,dm)
        dphi(:,7,dm) = 4 * ( dlam(:,1,dm) * lam(:,4) + &
                             lam(:,1) * dlam(:,4,dm) ) &
               - 12*db12(:,dm) - 12*db14(:,dm) + 32*db15(:,dm)
        dphi(:,8,dm) = 4 * ( dlam(:,2,dm) * lam(:,4) + &
                             lam(:,2) * dlam(:,4,dm) ) &
               - 12*db12(:,dm) - 12*db13(:,dm) + 32*db15(:,dm)
        dphi(:,9,dm) = 4 * ( dlam(:,3,dm) * lam(:,4) + &
                             lam(:,3) * dlam(:,4,dm) ) &
               - 12*db13(:,dm) - 12*db14(:,dm) + 32*db15(:,dm)
        dphi(:,10,dm) = ( 4 * lam(:,4) - 1 ) * dlam(:,4,dm) + &
               3*db12(:,dm) + 3*db13(:,dm) + 3*db14(:,dm) - 4*db15(:,dm)
      end do

      dphi(:,11,:) = 27*db11 - 108*db15
      dphi(:,12,:) = 27*db12 - 108*db15
      dphi(:,13,:) = 27*db13 - 108*db15
      dphi(:,14,:) = 27*db14 - 108*db15
      dphi(:,15,:) = 256*db15

    end if

  end subroutine shape_tetra_P2plus15


! Prism P1

  subroutine shape_prism_P1 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following points:
!
!        |----------/|
!       /|   _____/  |
!      /____/        |
!     // |  |        |     x center point value
!     |  |  x--      |    / coefficient of xi   (slope)
!     |  | /         |   -- coefficient of eta  (slope)
!     |  ----------_/|    | coefficient of zeta (slope)
!     | /    _____/
!     |/ ___/
!     |/
!          zeta
!            ^ _ eta
!            | /|
!            |/
!             ---> xi
!
!   Shape functions: 1, xi-1/3, eta-1/3, zeta
!
!   in the region [0,1]x[0,1]x[-1,1]  (xi+eta<=1)
!   The origin is half way between the edge separating the two orthogonal
!   square faces.
!

    phi(:,1) = 1
    phi(:,2) = xr(:,1) - 1/3._dp
    phi(:,3) = xr(:,2) - 1/3._dp
    phi(:,4) = xr(:,3)

    if ( present(dphi) ) then

      dphi(:,1,:) = 0
      dphi(:,2,1) = 1
      dphi(:,2,2) = 0
      dphi(:,2,3) = 0
      dphi(:,3,1) = 0
      dphi(:,3,2) = 1
      dphi(:,3,3) = 0
      dphi(:,4,1) = 0
      dphi(:,4,2) = 0
      dphi(:,4,3) = 1

    end if

  end subroutine shape_prism_P1


! Prism P1-Q1

  subroutine shape_prism_P1Q1 ( xr, phi, dphi )

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
!        4---------_/6
!       /|   _____/  |
!      /____/        |
!     5/ |           |
!     |  |           |
!     |  |           |
!     |  1---------_/3
!     | /    _____/
!     |/ ___/
!     2/
!
!          zeta
!            ^
!            |
!            |
!            /---> eta
!           /
!         |/_
!        xi
!
!   The reference coordinates xr=(xi,eta,zeta) are
!   in the region [0,1]x[0,1]x[-1,1]  (xi+eta<=1)
!   The origin is half way between the edge 1-4.

    integer, parameter, dimension(3,2) :: &
      p = reshape ( [1,2,3,4,5,6], [3,2] )
    integer :: i, j
    real(dp), dimension(size(xr,1),3) :: phi1
    real(dp), dimension(size(xr,1),3,2) :: dphi1
    real(dp), dimension(size(xr,1),2) :: phi2, dphi2

    if ( present(dphi) ) then

      call shape_triangle_P1 ( xr(:,1:2), phi1, dphi1 )
      call shape_line_P1 ( xr(:,3), phi2, dphi2 )

      do i = 1, 3
        do j = 1, 2
          phi(:,p(i,j))    =  phi1(:,i)   *  phi2(:,j)
          dphi(:,p(i,j),1) = dphi1(:,i,1) *  phi2(:,j)
          dphi(:,p(i,j),2) = dphi1(:,i,2) *  phi2(:,j)
          dphi(:,p(i,j),3) =  phi1(:,i)   * dphi2(:,j)
        end do
      end do

    else

      call shape_triangle_P1 ( xr(:,1:2), phi1 )
      call shape_line_P1 ( xr(:,3), phi2 )

      do i = 1, 3
        do j = 1, 2
          phi(:,p(i,j)) = phi1(:,i) * phi2(:,j)
        end do
      end do

    end if

  end subroutine shape_prism_P1Q1


! Prism P2-Q2

  subroutine shape_prism_P2Q2 ( xr, phi, dphi )

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
!       13---18----_/17
!    14 /|   _____/  |
!      /____/16      |
!    15/ 7    12    11
!     |  |   10      |
!     | 8|           |
!     9  1----6----_/5
!     |2/    _____/
!     |/ ___/4
!     3/
!
!          zeta
!            ^
!            |
!            |
!            /---> eta
!           /
!         |/_
!        xi
!
!   The reference coordinates xr=(xi,eta,zeta) are
!   in the region [0,1]x[0,1]x[-1,1]  (xi+eta<=1)
!   The origin is half way between the edge 1-13 in node 7.

    integer, parameter, dimension(6,3) :: &
      p = reshape ( [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18], [6,3] )
    integer :: i, j
    real(dp), dimension(size(xr,1),6) :: phi1
    real(dp), dimension(size(xr,1),6,2) :: dphi1
    real(dp), dimension(size(xr,1),3) :: phi2, dphi2

    if ( present(dphi) ) then

      call shape_triangle_P2 ( xr(:,1:2), phi1, dphi1 )
      call shape_line_P2 ( xr(:,3), phi2, dphi2 )

      do i = 1, 6
        do j = 1, 3
          phi(:,p(i,j))    =  phi1(:,i)   *  phi2(:,j)
          dphi(:,p(i,j),1) = dphi1(:,i,1) *  phi2(:,j)
          dphi(:,p(i,j),2) = dphi1(:,i,2) *  phi2(:,j)
          dphi(:,p(i,j),3) =  phi1(:,i)   * dphi2(:,j)
        end do
      end do

    else

      call shape_triangle_P2 ( xr(:,1:2), phi1 )
      call shape_line_P2 ( xr(:,3), phi2 )

      do i = 1, 6
        do j = 1, 3
          phi(:,p(i,j)) = phi1(:,i) * phi2(:,j)
        end do
      end do

    end if

  end subroutine shape_prism_P2Q2


! Pyramid P1

  subroutine shape_pyramid_P1 ( xr, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following points:
!
!
!           /\\               !
!          /| \ \             !
!         / /  \  \           !
!        / |  | \   \         !
!       |  |  x--|    \       !
!       |  |-/---\-----|--     x center point value (at [0,0,1/4]
!      / /        |    /       / coefficient of xi   (slope)
!     | /         \   /        -- coefficient of eta  (slope)
!     |/           | /         | coefficient of zeta (slope)
!     |------------|/
!
!          zeta
!            ^ _ eta
!            | /|
!            |/
!            ----> xi
!
!   The reference coordinates xr=(xi,eta,zeta) are
!   in the region [-1,1]x[-1,1]x[0,1].
!
!   Shape functions: 1, xi, eta, zeta-1/4
!

    phi(:,1) = 1
    phi(:,2) = xr(:,1)
    phi(:,3) = xr(:,2)
    phi(:,4) = xr(:,3) - 1/4._dp

    if ( present(dphi) ) then

      dphi(:,1,:) = 0
      dphi(:,2,1) = 1
      dphi(:,2,2) = 0
      dphi(:,2,3) = 0
      dphi(:,3,1) = 0
      dphi(:,3,2) = 1
      dphi(:,3,3) = 0
      dphi(:,4,1) = 0
      dphi(:,4,2) = 0
      dphi(:,4,3) = 1

    end if

  end subroutine shape_pyramid_P1


! Pyramid Q1-P1 rational

  subroutine shape_pyramid_Q1P1r ( xr, phi, dphi )

    use limits_m, only: EPS_PYRAMID

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following points:
!
!            5
!           /\\               !
!          /| \ \             !
!         / /  \  \           !
!        / |    \   \         !
!       |  |     |    \       !
!       |  4-----\-----3
!      / /        |    /
!     | /         \   /
!     |/           | /
!     1------------2/
!
!          zeta
!            ^ _ eta
!            | /|
!            |/
!            ----> xi
!
!   The reference coordinates xr=(xi,eta,zeta) are
!   in the region [-1,1]x[-1,1]x[0,1].
!
! Bedrosian, G. Shape functions and integration formulas for three-dimensional
! finite element analysis. Int. J. Numer. Meth. Engng. 35, 95-108 (1992).
!
! See also:
! Bergot, M., Cohen, G. & Durufle, M. Higher-order Finite Elements for Hybrid
! Meshes Using New Nodal Pyramidal Elements. J Sci Comput 42, 345-381 (2010).
!

    real(dp), dimension(size(xr,1)) :: fac

    where ( abs(1-xr(:,3)) <= EPS_PYRAMID )
      fac = 1 / EPS_PYRAMID
    else where
      fac = 1 / ( 1 - xr(:,3) )  ! 1/(1-zeta)
    end where

    phi(:,1) = ( 1 - xr(:,1) - xr(:,2) - xr(:,3) + xr(:,1) * xr(:,2) * fac ) / 4
    phi(:,2) = ( 1 + xr(:,1) - xr(:,2) - xr(:,3) - xr(:,1) * xr(:,2) * fac ) / 4
    phi(:,3) = ( 1 + xr(:,1) + xr(:,2) - xr(:,3) + xr(:,1) * xr(:,2) * fac ) / 4
    phi(:,4) = ( 1 - xr(:,1) + xr(:,2) - xr(:,3) - xr(:,1) * xr(:,2) * fac ) / 4
    phi(:,5) = xr(:,3)

    if ( present(dphi) ) then

      dphi(:,1,1) = ( -1 + xr(:,2) * fac ) / 4
      dphi(:,1,2) = ( -1 + xr(:,1) * fac ) / 4
      dphi(:,1,3) = ( -1 + xr(:,1) * xr(:,2) * fac**2 ) / 4
      dphi(:,2,1) = (  1 - xr(:,2) * fac ) / 4
      dphi(:,2,2) = ( -1 - xr(:,1) * fac ) / 4
      dphi(:,2,3) = ( -1 - xr(:,1) * xr(:,2) * fac**2 ) / 4
      dphi(:,3,1) = (  1 + xr(:,2) * fac ) / 4
      dphi(:,3,2) = (  1 + xr(:,1) * fac ) / 4
      dphi(:,3,3) = ( -1 + xr(:,1) * xr(:,2) * fac**2 ) / 4
      dphi(:,4,1) = ( -1 - xr(:,2) * fac ) / 4
      dphi(:,4,2) = (  1 - xr(:,1) * fac ) / 4
      dphi(:,4,3) = ( -1 - xr(:,1) * xr(:,2) * fac**2 ) / 4
      dphi(:,5,1) = 0
      dphi(:,5,2) = 0
      dphi(:,5,3) = 1

    end if

  end subroutine shape_pyramid_Q1P1r


! Pyramid Q2-P2 rational

  subroutine shape_pyramid_Q2P2r ( xr, phi, dphi )

    use limits_m, only: EPS_PYRAMID

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   The unknowns of the element are positioned at the following points:
!
!            14
!           /\\                !
!          /| \ \              !
!        13 /  \ 12
!        / |    \   \          !
!      10  |    11    \        !
!       |  7----6\-----5
!      / /        |    /
!     | 8     9   \   4
!     |/           | /
!     1-----2------3/
!
!          zeta
!            ^ _ eta
!            | /|
!            |/
!            ----> xi
!
!   The reference coordinates xr=(xi,eta,zeta) are
!   in the region [-1,1]x[-1,1]x[0,1].
!
! Graglia, R. D. & Gheorma, I.-L. Higher order interpolatory vector bases on
! pyramidal elements. IEEE Trans. Antennas Propagat. 47, 775-782 (1999).
!
! See also:
! Bergot, M., Cohen, G. & Durufle, M. Higher-order Finite Elements for Hybrid
! Meshes Using New Nodal Pyramidal Elements. J Sci Comput 42, 345-381 (2010).
!

    real(dp), dimension(size(xr,1)) :: fac
    real(dp), dimension(size(xr,1),5:14) :: t
    real(dp), dimension(size(xr,1),5:14,3) :: dt

    where ( abs(1-xr(:,3)) <= EPS_PYRAMID )
      fac = 1 / EPS_PYRAMID
    else where
      fac = 1 / ( 1 - xr(:,3) )  ! 1/(1-zeta)
    end where

    t(:,5) = xr(:,1) ** 2
    t(:,6) = xr(:,1) * xr(:,2)
    t(:,7) = xr(:,2) ** 2
    t(:,8) = xr(:,2) * xr(:,3)
    t(:,9) = xr(:,3) ** 2
    t(:,10) = xr(:,3) * xr(:,1)
    t(:,11) = xr(:,1) * xr(:,2) * fac
    t(:,12) = xr(:,1) ** 2 * xr(:,2) * fac
    t(:,13) = xr(:,1) * xr(:,2) ** 2 * fac
    t(:,14) = ( xr(:,1) * xr(:,2) * fac ) ** 2

    phi(:,1) = - xr(:,3) / 4 + t(:,6) / 2 + ( t(:,8) + t(:,9) + t(:,10)  &
               - t(:,11) - t(:,12) - t(:,13) + t(:,14) ) / 4
    phi(:,2) = - xr(:,2) / 2 + ( t(:,7) + t(:,8) ) / 2 &
               + ( t(:,12) - t(:,14) ) / 2
    phi(:,3) = - xr(:,3) / 4 - t(:,6) / 2 + ( t(:,8) + t(:,9) - t(:,10)  &
               + t(:,11) - t(:,12) + t(:,13) + t(:,14) ) / 4
    phi(:,4) =   xr(:,1) / 2 + t(:,5) / 2 - t(:,10) / 2  &
               - ( t(:,13) + t(:,14) ) / 2
    phi(:,5) = - xr(:,3) / 4 + t(:,6) / 2 + ( - t(:,8) + t(:,9) - t(:,10)  &
               - t(:,11) + t(:,12) + t(:,13) + t(:,14) ) / 4
    phi(:,6) =   xr(:,2) / 2 + ( t(:,7) - t(:,8) ) / 2  &
               - ( t(:,12) + t(:,14) ) / 2
    phi(:,7) = - xr(:,3) / 4 - t(:,6) / 2 + ( - t(:,8) + t(:,9) + t(:,10)  &
               + t(:,11) + t(:,12) - t(:,13) + t(:,14) ) / 4
    phi(:,8) = - xr(:,1) / 2 + t(:,5) / 2 + t(:,10) / 2  &
               + ( t(:,13) - t(:,14) ) / 2
    phi(:,9) = 1 - 2 * xr(:,3) - t(:,5) - t(:,7) + t(:,9) + t(:,14)
    phi(:,10) = xr(:,3) - t(:,6) - t(:,8) - t(:,9) - t(:,10) + t(:,11)
    phi(:,11) = xr(:,3) + t(:,6) - t(:,8) - t(:,9) + t(:,10) - t(:,11)
    phi(:,12) = xr(:,3) - t(:,6) + t(:,8) - t(:,9) + t(:,10) + t(:,11)
    phi(:,13) = xr(:,3) + t(:,6) + t(:,8) - t(:,9) - t(:,10) - t(:,11)
    phi(:,14) = - xr(:,3) + 2 * t(:,9)

    if ( present(dphi) ) then

      dt(:,5,1) = 2 * xr(:,1)
      dt(:,6,1) = xr(:,2)
      dt(:,7,1) = 0
      dt(:,8,1) = 0
      dt(:,9,1) = 0
      dt(:,10,1) = xr(:,3)
      dt(:,11,1) = xr(:,2) * fac
      dt(:,12,1) = 2 * xr(:,1) * xr(:,2) * fac
      dt(:,13,1) = xr(:,2) ** 2 * fac
      dt(:,14,1) = 2 * xr(:,1) * ( xr(:,2) * fac ) ** 2

      dt(:,5,2) = 0
      dt(:,6,2) = xr(:,1)
      dt(:,7,2) = 2* xr(:,2)
      dt(:,8,2) = xr(:,3)
      dt(:,9,2) = 0
      dt(:,10,2) = 0
      dt(:,11,2) = xr(:,1) * fac
      dt(:,12,2) = xr(:,1) ** 2 * fac
      dt(:,13,2) = 2 * xr(:,1) * xr(:,2) * fac
      dt(:,14,2) = 2 * xr(:,2) * ( xr(:,1) * fac ) ** 2

      dt(:,5,3) = 0
      dt(:,6,3) = 0
      dt(:,7,3) = 0
      dt(:,8,3) = xr(:,2)
      dt(:,9,3) = 2 * xr(:,3)
      dt(:,10,3) = xr(:,1)
      dt(:,11,3) = xr(:,1) * xr(:,2) * fac ** 2
      dt(:,12,3) = xr(:,1) ** 2 * xr(:,2) * fac ** 2
      dt(:,13,3) = xr(:,1) * xr(:,2) ** 2 * fac ** 2
      dt(:,14,3) = 2 * ( xr(:,1) * xr(:,2) ) ** 2 * fac ** 3

      dphi(:,1,1) = dt(:,6,1) / 2 + ( dt(:,8,1) + dt(:,9,1) + dt(:,10,1)  &
               - dt(:,11,1) - dt(:,12,1) - dt(:,13,1) + dt(:,14,1) ) / 4
      dphi(:,2,1) = ( dt(:,7,1) + dt(:,8,1) ) / 2 &
                              + ( dt(:,12,1) - dt(:,14,1) ) / 2
      dphi(:,3,1) = - dt(:,6,1) / 2 + ( dt(:,8,1) + dt(:,9,1) - dt(:,10,1)  &
               + dt(:,11,1) - dt(:,12,1) + dt(:,13,1) + dt(:,14,1) ) / 4
      dphi(:,4,1) = 1._dp / 2 + dt(:,5,1) / 2 - dt(:,10,1) / 2  &
                              - ( dt(:,13,1) + dt(:,14,1) ) / 2
      dphi(:,5,1) = dt(:,6,1) / 2 + ( - dt(:,8,1) + dt(:,9,1) - dt(:,10,1)  &
               - dt(:,11,1) + dt(:,12,1) + dt(:,13,1) + dt(:,14,1) ) / 4
      dphi(:,6,1) = ( dt(:,7,1) - dt(:,8,1) ) / 2  &
                              - ( dt(:,12,1) + dt(:,14,1) ) / 2
      dphi(:,7,1) = - dt(:,6,1) / 2 + ( - dt(:,8,1) + dt(:,9,1) + dt(:,10,1)  &
               + dt(:,11,1) + dt(:,12,1) - dt(:,13,1) + dt(:,14,1) ) / 4
      dphi(:,8,1) = - 1._dp / 2 + dt(:,5,1) / 2 + dt(:,10,1) / 2  &
                              + ( dt(:,13,1) - dt(:,14,1) ) / 2
      dphi(:,9,1) = - dt(:,5,1) - dt(:,7,1) + dt(:,9,1) + dt(:,14,1)
      dphi(:,10,1) = - dt(:,6,1) - dt(:,8,1) - dt(:,9,1) - dt(:,10,1) &
                     + dt(:,11,1)
      dphi(:,11,1) = dt(:,6,1) - dt(:,8,1) - dt(:,9,1) + dt(:,10,1) - dt(:,11,1)
      dphi(:,12,1) = - dt(:,6,1) + dt(:,8,1) - dt(:,9,1) + dt(:,10,1) &
                     + dt(:,11,1)
      dphi(:,13,1) = dt(:,6,1) + dt(:,8,1) - dt(:,9,1) - dt(:,10,1) - dt(:,11,1)
      dphi(:,14,1) = 2 * dt(:,9,1)

      dphi(:,1,2) = dt(:,6,2) / 2 + ( dt(:,8,2) + dt(:,9,2) + dt(:,10,2)  &
               - dt(:,11,2) - dt(:,12,2) - dt(:,13,2) + dt(:,14,2) ) / 4
      dphi(:,2,2) = - 1._dp / 2 + ( dt(:,7,2) + dt(:,8,2) ) / 2 &
                                + ( dt(:,12,2) - dt(:,14,2) ) / 2
      dphi(:,3,2) = - dt(:,6,2) / 2 + ( dt(:,8,2) + dt(:,9,2) - dt(:,10,2)  &
               + dt(:,11,2) - dt(:,12,2) + dt(:,13,2) + dt(:,14,2) ) / 4
      dphi(:,4,2) = dt(:,5,2) / 2 - dt(:,10,2) / 2  &
                                - ( dt(:,13,2) + dt(:,14,2) ) / 2
      dphi(:,5,2) = dt(:,6,2) / 2 + ( - dt(:,8,2) + dt(:,9,2) - dt(:,10,2)  &
               - dt(:,11,2) + dt(:,12,2) + dt(:,13,2) + dt(:,14,2) ) / 4
      dphi(:,6,2) = 1._dp / 2 + ( dt(:,7,2) - dt(:,8,2) ) / 2  &
                                - ( dt(:,12,2) + dt(:,14,2) ) / 2
      dphi(:,7,2) = - dt(:,6,2) / 2 + ( - dt(:,8,2) + dt(:,9,2) + dt(:,10,2)  &
               + dt(:,11,2) + dt(:,12,2) - dt(:,13,2) + dt(:,14,2) ) / 4
      dphi(:,8,2) = dt(:,5,2) / 2 + dt(:,10,2) / 2  &
                                + ( dt(:,13,2) - dt(:,14,2) ) / 2
      dphi(:,9,2) = - dt(:,5,2) - dt(:,7,2) + dt(:,9,2) + dt(:,14,2)
      dphi(:,10,2) = - dt(:,6,2) - dt(:,8,2) - dt(:,9,2) - dt(:,10,2) &
                     + dt(:,11,2)
      dphi(:,11,2) = dt(:,6,2) - dt(:,8,2) - dt(:,9,2) + dt(:,10,2) - dt(:,11,2)
      dphi(:,12,2) = - dt(:,6,2) + dt(:,8,2) - dt(:,9,2) + dt(:,10,2) &
                     + dt(:,11,2)
      dphi(:,13,2) = dt(:,6,2) + dt(:,8,2) - dt(:,9,2) - dt(:,10,2) - dt(:,11,2)
      dphi(:,14,2) = 2 * dt(:,9,2)

      dphi(:,1,3) = - 1._dp / 4 + dt(:,6,3) / 2 + ( dt(:,8,3) + dt(:,9,3) &
        + dt(:,10,3) - dt(:,11,3) - dt(:,12,3) - dt(:,13,3) + dt(:,14,3) ) / 4
      dphi(:,2,3) = ( dt(:,7,3) + dt(:,8,3) ) / 2 &
                                + ( dt(:,12,3) - dt(:,14,3) ) / 2
      dphi(:,3,3) = - 1._dp / 4 - dt(:,6,3) / 2 + ( dt(:,8,3) + dt(:,9,3) &
        - dt(:,10,3) + dt(:,11,3) - dt(:,12,3) + dt(:,13,3) + dt(:,14,3) ) / 4
      dphi(:,4,3) = dt(:,5,3) / 2 - dt(:,10,3) / 2  &
                                - ( dt(:,13,3) + dt(:,14,3) ) / 2
      dphi(:,5,3) = - 1._dp / 4 + dt(:,6,3) / 2 + ( - dt(:,8,3) + dt(:,9,3) &
        - dt(:,10,3) - dt(:,11,3) + dt(:,12,3) + dt(:,13,3) + dt(:,14,3) ) / 4
      dphi(:,6,3) = ( dt(:,7,3) - dt(:,8,3) ) / 2  &
                                - ( dt(:,12,3) + dt(:,14,3) ) / 2
      dphi(:,7,3) = - 1._dp / 4 - dt(:,6,3) / 2 + ( - dt(:,8,3) + dt(:,9,3) &
        + dt(:,10,3) + dt(:,11,3) + dt(:,12,3) - dt(:,13,3) + dt(:,14,3) ) / 4
      dphi(:,8,3) = dt(:,5,3) / 2 + dt(:,10,3) / 2  &
                                + ( dt(:,13,3) - dt(:,14,3) ) / 2
      dphi(:,9,3) = - 2 - dt(:,5,3) - dt(:,7,3) + dt(:,9,3) + dt(:,14,3)
      dphi(:,10,3) = 1 - dt(:,6,3) - dt(:,8,3) - dt(:,9,3) - dt(:,10,3) &
                     + dt(:,11,3)
      dphi(:,11,3) = 1 + dt(:,6,3) - dt(:,8,3) - dt(:,9,3) + dt(:,10,3) &
                     - dt(:,11,3)
      dphi(:,12,3) = 1 - dt(:,6,3) + dt(:,8,3) - dt(:,9,3) + dt(:,10,3) &
                     + dt(:,11,3)
      dphi(:,13,3) = 1 + dt(:,6,3) + dt(:,8,3) - dt(:,9,3) - dt(:,10,3) &
                     - dt(:,11,3)
      dphi(:,14,3) = - 1 + 2 * dt(:,9,3)

    end if

  end subroutine shape_pyramid_Q2P2r

end module shapefunc_standard_m
