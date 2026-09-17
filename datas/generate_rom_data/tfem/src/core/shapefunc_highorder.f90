
! Copyright (C) 2021-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! This module defines high-order shape functions

module shapefunc_highorder_m

  use kind_defs_m
  use shapefunc_standard_m, only: barycentric

  implicit none


! interface to support two interfaces for shape_line_highorder

  interface shape_line_Pp
    module procedure shape_line_Pp_1, shape_line_Pp_2
  end interface shape_line_Pp

contains


! Routine for computing shapefunctions and derivatives for Lagrangian
! interpolation on a line

  subroutine shape_line_lagrange ( x, xp, phi, dphi )

!   coordinates for evaluation of the phi and dphi
    real(dp), dimension(:), intent(in) :: x

!   positions of the nodes for interpolation (where the unknowns are given)
    real(dp), dimension(:), intent(in) :: xp

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), dimension(:,:), intent(out) :: phi

!   derivative of the shape function: dphi(i,j), with i the point in space
!   j the unknown
    real(dp), dimension(:,:), intent(out), optional :: dphi

    integer :: n, m, i, j, k
    real(dp), allocatable, dimension(:) :: vec
    real(dp), allocatable, dimension(:,:) :: mat
    real(dp), allocatable, dimension(:,:,:) :: mat1, mat2

    m = size(x)  ! number of evaluation points
    n = size(xp) ! number of nodes

    allocate( vec(n), mat(n,n), mat1(m,n,n), mat2(m,n,n) )

!   fill mat with xp_i-xp_j and 1 on the diagonal

    do i = 1, n
      do j = i+1, n
        mat(i,j) = xp(j) - xp(i)
        mat(j,i) = -mat(i,j)
      end do
      mat(i,i) = 1
    end do

!   the scaling factors to require the Kronecker delta property

    vec = 1 / product ( mat, dim=1 )

!   fill columns of mat1 with multiple copies of x-xp_i and 1 on the diagonal

    do i = 1, n
      mat1(:,i,1) = x - xp(i)
    end do
    do j = 2, n
      mat1(:,:,j) = mat1(:,:,1)
    end do
    do i = 1, n
      mat1(:,i,i) = 1
    end do

!   shape functions

    phi = product ( mat1, dim=2 )

!   scale to one

    do k = 1, m
      phi(k,:) = phi(k,:) * vec
    end do

    if ( present(dphi) ) then

!     derivative of the shape functions

      dphi = 0

      do i = 1, n

        mat2 = mat1

        mat2(:,i,:) = 1 ! unit factor in row i
        mat2(:,:,i) = 0 ! zero factors in column i

        dphi = dphi + product ( mat2, dim=2 )

      end do

      do k = 1, m
        dphi(k,:) = dphi(k,:) * vec
      end do

    end if

  end subroutine shape_line_lagrange


! line Pp (one dimension in xr)

  subroutine shape_line_Pp_1 ( xr, p, phi, dphi )

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

!   The unknowns of the element are positioned at the following nodal points:
!
!     1 --- 2 .... p --- p+1
!
!   The reference coordinates xr are in the region [-1,1]

    integer :: i
    real(dp), allocatable, dimension(:) :: xp


    select case (p)

    case(1) ! P1

      phi(:,1) = (1-xr)/2
      phi(:,2) = (1+xr)/2

      if ( present(dphi) ) then
        dphi(:,1) = -0.5_dp
        dphi(:,2) =  0.5_dp
      end if

    case(2) ! P2

      phi(:,1) = -(1-xr)*xr/2
      phi(:,2) = 1-xr**2
      phi(:,3) = (1+xr)*xr/2

      if ( present(dphi) ) then
        dphi(:,1) = -(1-2*xr)/2
        dphi(:,2) = -2*xr
        dphi(:,3) = (1+2*xr)/2
      end if

    case(3) ! P3

      phi(:,1) = -(9*xr**2-1)*(xr-1)/16
      phi(:,2) = 9*(xr**2-1)*(3*xr-1)/16
      phi(:,3) = -9*(xr**2-1)*(3*xr+1)/16
      phi(:,4) = (xr+1)*(9*xr**2-1)/16

      if ( present(dphi) ) then
        dphi(:,1) = -(27*xr**2-18*xr-1)/16
        dphi(:,2) = 9*(9*xr**2-2*xr-3)/16
        dphi(:,3) = -9*(9*xr**2+2*xr-3)/16
        dphi(:,4) = (27*xr**2+18*xr-1)/16
      end if

    case(4) ! P4

      phi(:,1) = (4*xr**2-1)*(xr**2-xr)/6
      phi(:,2) = -4*(xr**2-1)*(2*xr**2-xr)/3
      phi(:,3) = (xr**2-1)*(4*xr**2-1)
      phi(:,4) = -4*(xr**2-1)*(2*xr**2+xr)/3
      phi(:,5) = (4*xr**2-1)*(xr**2+xr)/6

      if ( present(dphi) ) then
        dphi(:,1) = (16*xr**3-12*xr**2-2*xr+1)/6
        dphi(:,2) = -4*(8*xr**3-3*xr**2-4*xr+1)/3
        dphi(:,3) = 16*xr**3-10*xr
        dphi(:,4) = -4*(8*xr**3+3*xr**2-4*xr-1)/3
        dphi(:,5) = (16*xr**3+12*xr**2-2*xr-1)/6
      end if

    case(5) ! P5

      phi(:,1) = -(25*xr**2-9)*(25*xr**2-1)*(xr-1)/768
      phi(:,2) = 25*(xr**2-1)*(25*xr**2-1)*(5*xr-3)/768
      phi(:,3) = -25*(xr**2-1)*(25*xr**2-9)*(5*xr-1)/384
      phi(:,4) = 25*(xr**2-1)*(25*xr**2-9)*(5*xr+1)/384
      phi(:,5) = -25*(xr**2-1)*(25*xr**2-1)*(5*xr+3)/768
      phi(:,6) = (25*xr**2-9)*(25*xr**2-1)*(xr+1)/768

      if ( present(dphi) ) then
        dphi(:,1) = -(3125*xr**4-2500*xr**3-750*xr**2+500*xr+9)/768
        dphi(:,2) = 25*(625*xr**4-300*xr**3-390*xr**2+156*xr+5)/768
        dphi(:,3) = -25*(625*xr**4-100*xr**3-510*xr**2+68*xr+45)/384
        dphi(:,4) = 25*(625*xr**4+100*xr**3-510*xr**2-68*xr+45)/384
        dphi(:,5) = -25*(625*xr**4+300*xr**3-390*xr**2-156*xr+5)/768
        dphi(:,6) = (3125*xr**4+2500*xr**3-750*xr**2-500*xr+9)/768
      end if

    case(6:) ! P6 and larger

      xp = 2*[ (i,i=0,p) ]/real(p,dp) - 1

      call shape_line_lagrange ( xr, xp, phi, dphi )

    case default

      write(*,'(/a/a,i0/)') 'Error shape_line_Pp_1:', &
       ' order not available = ', p
      stop

    end select

  end subroutine shape_line_Pp_1


! Line Pp (two dimensions in xr)

  subroutine shape_line_Pp_2 ( xr, p, phi, dphi )

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


!   The unknowns of the element are positioned at the following nodal points:
!
!     1 --- 2 .... p --- p+1
!
!   The reference coordinates xr are in the region [-1,1]

    integer :: i
    real(dp), allocatable, dimension(:) :: xp


    select case (p)

    case(1) ! P1

      phi(:,1) = (1-xr(:,1))/2
      phi(:,2) = (1+xr(:,1))/2

      if ( present(dphi) ) then
        dphi(:,1,1) = -0.5_dp
        dphi(:,2,1) =  0.5_dp
      end if

    case(2) ! P2

      phi(:,1) = -(1-xr(:,1))*xr(:,1)/2
      phi(:,2) = 1-xr(:,1)**2
      phi(:,3) = (1+xr(:,1))*xr(:,1)/2

      if ( present(dphi) ) then
        dphi(:,1,1) = -(1-2*xr(:,1))/2
        dphi(:,2,1) = -2*xr(:,1)
        dphi(:,3,1) = (1+2*xr(:,1))/2
      end if

    case(3) ! P3

      phi(:,1) = -(9*xr(:,1)**2-1)*(xr(:,1)-1)/16
      phi(:,2) = 9*(xr(:,1)**2-1)*(3*xr(:,1)-1)/16
      phi(:,3) = -9*(xr(:,1)**2-1)*(3*xr(:,1)+1)/16
      phi(:,4) = (xr(:,1)+1)*(9*xr(:,1)**2-1)/16

      if ( present(dphi) ) then
        dphi(:,1,1) = -(27*xr(:,1)**2-18*xr(:,1)-1)/16
        dphi(:,2,1) = 9*(9*xr(:,1)**2-2*xr(:,1)-3)/16
        dphi(:,3,1) = -9*(9*xr(:,1)**2+2*xr(:,1)-3)/16
        dphi(:,4,1) = (27*xr(:,1)**2+18*xr(:,1)-1)/16
      end if

    case(4) ! P4

      phi(:,1) = (4*xr(:,1)**2-1)*(xr(:,1)**2-xr(:,1))/6
      phi(:,2) = -4*(xr(:,1)**2-1)*(2*xr(:,1)**2-xr(:,1))/3
      phi(:,3) = (xr(:,1)**2-1)*(4*xr(:,1)**2-1)
      phi(:,4) = -4*(xr(:,1)**2-1)*(2*xr(:,1)**2+xr(:,1))/3
      phi(:,5) = (4*xr(:,1)**2-1)*(xr(:,1)**2+xr(:,1))/6

      if ( present(dphi) ) then
        dphi(:,1,1) = (16*xr(:,1)**3-12*xr(:,1)**2-2*xr(:,1)+1)/6
        dphi(:,2,1) = -4*(8*xr(:,1)**3-3*xr(:,1)**2-4*xr(:,1)+1)/3
        dphi(:,3,1) = 16*xr(:,1)**3-10*xr(:,1)
        dphi(:,4,1) = -4*(8*xr(:,1)**3+3*xr(:,1)**2-4*xr(:,1)-1)/3
        dphi(:,5,1) = (16*xr(:,1)**3+12*xr(:,1)**2-2*xr(:,1)-1)/6
      end if

    case(5) ! P5

      phi(:,1) = -(25*xr(:,1)**2-9)*(25*xr(:,1)**2-1)*(xr(:,1)-1)/768
      phi(:,2) = 25*(xr(:,1)**2-1)*(25*xr(:,1)**2-1)*(5*xr(:,1)-3)/768
      phi(:,3) = -25*(xr(:,1)**2-1)*(25*xr(:,1)**2-9)*(5*xr(:,1)-1)/384
      phi(:,4) = 25*(xr(:,1)**2-1)*(25*xr(:,1)**2-9)*(5*xr(:,1)+1)/384
      phi(:,5) = -25*(xr(:,1)**2-1)*(25*xr(:,1)**2-1)*(5*xr(:,1)+3)/768
      phi(:,6) = (25*xr(:,1)**2-9)*(25*xr(:,1)**2-1)*(xr(:,1)+1)/768

      if ( present(dphi) ) then
        dphi(:,1,1) = -(3125*xr(:,1)**4-2500*xr(:,1)**3-750*xr(:,1)**2 &
                      +500*xr(:,1)+9)/768
        dphi(:,2,1) = 25*(625*xr(:,1)**4-300*xr(:,1)**3-390*xr(:,1)**2 &
                        +156*xr(:,1)+5)/768
        dphi(:,3,1) = -25*(625*xr(:,1)**4-100*xr(:,1)**3-510*xr(:,1)**2 &
                         +68*xr(:,1)+45)/384
        dphi(:,4,1) = 25*(625*xr(:,1)**4+100*xr(:,1)**3-510*xr(:,1)**2 &
                        -68*xr(:,1)+45)/384
        dphi(:,5,1) = -25*(625*xr(:,1)**4+300*xr(:,1)**3-390*xr(:,1)**2 &
                         -156*xr(:,1)+5)/768
        dphi(:,6,1) = (3125*xr(:,1)**4+2500*xr(:,1)**3-750*xr(:,1)**2 &
                     -500*xr(:,1)+9)/768
      end if

    case(6:) ! P6 and larger

      xp = 2*[ (i,i=0,p) ]/real(p,dp) - 1

      call shape_line_lagrange ( xr(:,1), xp, phi, dphi(:,:,1) )

    case default

      write(*,'(/a/a,i0/)') 'Error shape_line_Pp_2:', &
       ' order not available = ', p
      stop

    end select

  end subroutine shape_line_Pp_2


! Quadrilateral Qp

  subroutine shape_quad_Qp ( xr, p, phi, dphi )

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


    integer :: i, j
    integer, allocatable, dimension(:,:) :: s
    real(dp), allocatable, dimension(:,:) :: phi1, phi2, dphi1, dphi2

    allocate(phi1(size(xr,1),p+1))
    allocate(phi2,dphi1,dphi2,mold=phi1)
    s = reshape ( [(i,i=1,(p+1)**2)], [p+1,p+1] )


    if ( present(dphi) ) then

      call shape_line_Pp ( xr(:,1), p, phi1, dphi1 )
      call shape_line_Pp ( xr(:,2), p, phi2, dphi2 )

      do i = 1, p+1
        do j = 1, p+1
          phi(:,s(i,j))    =  phi1(:,i) *  phi2(:,j)
          dphi(:,s(i,j),1) = dphi1(:,i) *  phi2(:,j)
          dphi(:,s(i,j),2) =  phi1(:,i) * dphi2(:,j)
        end do
      end do

    else

      call shape_line_Pp ( xr(:,1), p, phi1 )
      call shape_line_Pp ( xr(:,2), p, phi2 )

      do i = 1, p+1
        do j = 1, p+1
          phi(:,s(i,j)) =  phi1(:,i) *  phi2(:,j)
        end do
      end do

    end if

  end subroutine shape_quad_Qp


! Hexahedron Qp

  subroutine shape_hexa_Qp ( xr, p, phi, dphi )

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


    integer :: i, j, k
    integer, allocatable, dimension(:,:,:) :: s
    real(dp), allocatable, dimension(:,:) :: phi1, phi2, phi3, dphi1, dphi2, &
      dphi3

    allocate(phi1(size(xr,1),p+1))
    allocate(phi2,phi3,dphi1,dphi2,dphi3,mold=phi1)
    s = reshape ( [(i,i=1,(p+1)**3)], [p+1,p+1,p+1] )


    if ( present(dphi) ) then

      call shape_line_Pp ( xr(:,1), p, phi1, dphi1 )
      call shape_line_Pp ( xr(:,2), p, phi2, dphi2 )
      call shape_line_Pp ( xr(:,3), p, phi3, dphi3 )

      do i = 1, p+1
        do j = 1, p+1
          do k = 1, p+1
            phi(:,s(i,j,k))    =  phi1(:,i) *  phi2(:,j) *  phi3(:,k)
            dphi(:,s(i,j,k),1) = dphi1(:,i) *  phi2(:,j) *  phi3(:,k)
            dphi(:,s(i,j,k),2) =  phi1(:,i) * dphi2(:,j) *  phi3(:,k)
            dphi(:,s(i,j,k),3) =  phi1(:,i) *  phi2(:,j) * dphi3(:,k)
          end do
        end do
      end do

    else

      call shape_line_Pp ( xr(:,1), p, phi1 )
      call shape_line_Pp ( xr(:,2), p, phi2 )
      call shape_line_Pp ( xr(:,3), p, phi3 )

      do i = 1, p+1
        do j = 1, p+1
          do k = 1, p+1
            phi(:,s(i,j,k)) = phi1(:,i) * phi2(:,j) * phi3(:,k)
          end do
        end do
      end do

    end if

  end subroutine shape_hexa_Qp


! Routine for computing shapefunctions and derivatives for Lagrangian
! interpolation on a line using barycentric coordinates.
! This is not meant to be used as a general purpose routine and is
! only to help generate shape functions on a triangle.

  subroutine shape_line_lagrange_bar ( lam1, lam2, xp1, xp2, phi, dphi1, dphi2 )

!   barycentric coordinates for evaluation of phi, dphi1, dphi2
!   Note: lam1+lam2 doesn't have to be equal to 1 in order to be able to use
!         this routine for a triangle.
    real(dp), dimension(:), intent(in) :: lam1, lam2

!   positions of the nodes along the line for interpolation (where the
!   unknowns are given) in terms of the values of lam1, lam2 respectively.
    real(dp), dimension(:), intent(in) :: xp1, xp2

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), dimension(:,:), intent(out) :: phi

!   derivative of the shape function wrt lam1 and lam2: dphi1(i,j), dphi2(i,j)
!   with i the point in space j the unknown
    real(dp), intent(out), dimension(:,:), optional :: dphi1, dphi2

    integer :: n, m, i, j, k
    real(dp), allocatable, dimension(:) :: vec
    real(dp), allocatable, dimension(:,:) :: mat
    real(dp), allocatable, dimension(:,:,:) :: mat1, mat2

    m = size(lam1)  ! number of evaluation points
    n = size(xp1)    ! number of nodes

    allocate( vec(n), mat(n,n), mat1(m,n,n), mat2(m,n,n) )

!   fill mat with xp_j-xp_i and 1 on the diagonal

    do i = 1, n
      do j = i+1, n
        mat(i,j) = xp2(j) - xp2(i)
        mat(j,i) = xp1(i) - xp1(j)
      end do
      mat(i,i) = 1
    end do

!   the scaling factors to require the Kronecker delta property

    vec = 1 / product ( mat, dim=1 )

!   fill columns of mat1 with lam1-xp1_i (lower diagonal), lam2-xp2_i
!   (upper diagonal) and 1 on the diagonal

    do i = 1, n
      mat1(:,i,1) = lam1 - xp1(i)
      mat1(:,i,n) = lam2 - xp2(i)
    end do
    do j = 2, n-1
      mat1(:,1:j-1,j) = mat1(:,1:j-1,n)
      mat1(:,j+1:n,j) = mat1(:,j+1:n,1)
    end do
    do i = 1, n
      mat1(:,i,i) = 1
    end do

!   shape functions

    phi = product ( mat1, dim=2 )

!   scale to one

    do k = 1, m
      phi(k,:) = phi(k,:) * vec
    end do

    if ( present(dphi1) .and. present(dphi2) ) then

!     derivative of the shape functions

!     derivative wrt lambda1

      dphi1 = 0

      do i = 1, n

        mat2 = mat1

        mat2(:,i,1:i-1) = 1 ! unit factors in row i wrt lambda1
        mat2(:,i,i:n) = 0 ! zero factors in row i to exclude columns of lambda2

        dphi1 = dphi1 + product ( mat2, dim=2 )

      end do

!     derivative wrt lambda2

      dphi2 = 0

      do i = 1, n

        mat2 = mat1

        mat2(:,i,1:i) = 0 ! zero factors in row i to exclude columns of lambda1
        mat2(:,i,i+1:n) = 1 ! unit factors in row i wrt lambda2

        dphi2 = dphi2 + product ( mat2, dim=2 )

      end do

      do k = 1, m
        dphi1(k,:) = dphi1(k,:) * vec
        dphi2(k,:) = dphi2(k,:) * vec
      end do

    end if

  end subroutine shape_line_lagrange_bar


! Routine for computing shapefunctions and derivatives for bubble functions
! on a triangle using barycentric coordinates.
! The purpose of this routine is to help generate shape functions on a triangle.

  subroutine shape_bubble_triangle ( lam, zl, phi, dphidlam )

!   barycentric coordinates for evaluation of phi, dphidlam
!   lam(i,j), with i the evaluation point and j=1,2,3
    real(dp), dimension(:,:), intent(in) :: lam

!   the zero level, where the bubble function = 0
    real(dp), intent(in) :: zl

!   shape function phi(i), with i the evaluation point
    real(dp), dimension(:), intent(out) :: phi

!   derivative of the shape function wrt lam: dphidlam(i,j)
!   with i the point in space j=1,2,3
    real(dp), intent(out), dimension(:,:), optional :: dphidlam

    phi = (lam(:,1)-zl)*(lam(:,2)-zl)*(lam(:,3)-zl)

    if ( present(dphidlam) ) then

      dphidlam(:,1) = (lam(:,2)-zl)*(lam(:,3)-zl)
      dphidlam(:,2) = (lam(:,3)-zl)*(lam(:,1)-zl)
      dphidlam(:,3) = (lam(:,1)-zl)*(lam(:,2)-zl)

    end if

  end subroutine shape_bubble_triangle


! Lagrange interpolation on a triangle (regular node distribution)
! Generic routine for any value of p.

  subroutine shape_triangle_lagrange_regular ( xr, p, phi, dphi )

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
! (p+1)(p+2)/2
!     |   \                     !
!     .      .
!     .         .
!     .           .
!   p+2 --- ... --- 2p+1
!     |                \        !
!     1 --- 2 .... p --- p+1
!
!  eta
!   ^
!   |
!    --> xi
!
!   The reference coordinates xr=(xi,eta) are in the region [-1,1]x[-1,1].

    integer :: i, mat(p+1,p+1), nr, j, mj, r, nne, s, l
    integer, dimension(:), allocatable :: ind, ind1, ind2, ind3
    real(dp) :: delta
    real(dp) :: lam(size(xr,1),3), dlam(size(xr,1),3,2)
    real(dp) :: philam(size(xr,1)), dphidlam(size(xr,1),3)
    real(dp), dimension(p+1) :: xpl ! local array for position along an edge
    real(dp), dimension(size(xr,1),size(phi,2)) :: phitmp, dphi1, dphi2
    real(dp), dimension(:,:), allocatable :: lamring


!   regular distribution of reference coordinates of nodal points

    delta = 1 / real(p,dp)  ! grid increment
    xpl = [ (i, i=0,p) ] * delta

!   barycentric coordinates of the points where phi, dphi must be computed

    if ( present(dphi) ) then
      call barycentric ( xr, lam, dlam )
    else
      call barycentric ( xr, lam )
    end if

!   fill mat with local nodal numbers

    do j = 1, p+1
       mj= (j-1)*(2*p+4-j)/2
       do i = 1, p+2-j
          mat(i,j) = i+mj
       end do
    end do

!   number of edge rings
    nr = (p+2)/3

!   loop over all edge rings

    do r = 1, nr

      nne = p+4-3*r  ! number of edge nodes
      s = p+3-2*r    ! end vertex index in matrix

!     edge shape factors

!     edge 1
      ind1 = mat(r:s,r)
      call edge ( 1, 2, ind1 )

!     edge 2
      allocate(ind2(nne))
      do i = 1, nne
        ind2(i) = mat(s-i+1,r+i-1)
      end do
      call edge ( 2, 3, ind2 )

!     edge 3
      ind3 = mat(r,s:r:-1)
      call edge ( 3, 1, ind3 )

      ind = [ ind1(1:nne-1), ind2(1:nne-1), ind3(1:nne-1) ]

      deallocate ( ind2 )

      if ( r == 1 ) cycle

!     contribution of factors from external rings

      allocate( lamring(3*(nne-1),3) )

      call fill_lamring

      call shape_factor_external_rings ( ner=r-1 )

      deallocate ( ind, lamring )

    end do

!   center node

    if ( mod(p,3) == 0 ) then

      ind = [ mat(nr+1,nr+1) ]
      lamring = reshape ( [ 1, 1, 1 ] / 3._dp, [1,3] )
      phi(:,ind(1)) = 1
      if ( present(dphi) ) then
        do l = 1, 2
          dphi(:,ind(1),l) = 0
        end do
      end if

      call shape_factor_external_rings ( ner=nr )

    end if

  contains

    subroutine edge ( i, j, ind )

      integer, intent(in) :: i, j
      integer, dimension(:), intent(in) :: ind

      integer :: k, l

      call shape_line_lagrange_bar ( lam(:,i), lam(:,j), xpl(s:r:-1), xpl(r:s),&
                phitmp(:,1:nne), dphi1(:,1:nne), dphi2(:,1:nne) )
      phi(:,ind) = phitmp(:,1:nne)
      if ( present(dphi) ) then
        do k = 1, nne
          do l = 1, 2
            dphi(:,ind(k),l) = dphi1(:,k) * dlam(:,i,l) + &
                               dphi2(:,k) * dlam(:,j,l)
          end do
        end do
      end if

    end subroutine edge

    subroutine fill_lamring

      real(dp), dimension(:,:), allocatable :: xr1, xr2, xr3, lam1, lam2, lam3

      allocate( xr1(nne-1,2), xr2(nne-1,2), xr3(nne-1,2) )
      xr1(:,1) = [ (r-2+i, i=1,nne-1) ] * delta
      xr1(:,2) = (r-1) * delta
      xr2(:,1) = [ (s-i, i=1,nne-1) ] * delta
      xr2(:,2) = [ (r-2+i, i=1,nne-1) ] * delta
      xr3(:,1) = (r-1) * delta
      xr3(:,2) = [ (s-i, i=1,nne-1) ] * delta
      allocate( lam1(nne-1,3), lam2(nne-1,3), lam3(nne-1,3) )
      call barycentric ( xr1, lam1 )
      call barycentric ( xr2, lam2 )
      call barycentric ( xr3, lam3 )
      lamring(1:nne-1,:) = lam1
      lamring(nne:2*nne-2,:) = lam2
      lamring(2*nne-1:3*nne-3,:) = lam3

    end subroutine fill_lamring

    subroutine shape_factor_external_rings ( ner )

      integer, intent(in) :: ner

      integer :: er, j, l
      real(dp), dimension(size(lamring,1)) :: philamring

      do er = 1, ner

        call shape_bubble_triangle ( lam, zl=(er-1)*delta, phi=philam, &
          dphidlam=dphidlam )

        call shape_bubble_triangle ( lamring, zl=(er-1)*delta, phi=philamring )

        if ( present(dphi) ) then

          do j = 1, size(lamring,1)
            do l = 1, 2
              dphi(:,ind(j),l) = ( dphi(:,ind(j),l) * philam + &
                                     phi(:,ind(j)) * ( &
                                        dphidlam(:,1) * dlam(:,1,l) + &
                                        dphidlam(:,2) * dlam(:,2,l) + &
                                        dphidlam(:,3) * dlam(:,3,l) ) ) &
                                             /  philamring(j)
            end do
          end do

        end if

        do j = 1, size(lamring,1)
          phi(:,ind(j)) = phi(:,ind(j)) * philam / philamring(j)
        end do

      end do

    end subroutine shape_factor_external_rings

  end subroutine shape_triangle_lagrange_regular


! Triangle Pp (equidistant node distribution)
! Analytic functions for p<=5, generic algorithm for p>5.

  subroutine shape_triangle_Pp ( xr, p, phi, dphi )

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
! (p+1)(p+2)/2
!     |   \                      !
!     .      .
!     .         .
!     .           .
!   p+2 --- ... --- 2p+1
!     |                \         !
!     1 --- 2 .... p --- p+1
!
!  eta
!   ^
!   |
!    --> xi
!
!   The reference coordinates xr=(xi,eta) are in the region [-1,1]x[-1,1].

    integer :: dm

    real(dp) :: lam(size(xr,1),3), dlam(size(xr,1),3,2)


!   barycentric coordinates lambda1, lambda2, lambda3

    if ( present(dphi) ) then
      call barycentric ( xr, lam, dlam )
    else
      call barycentric ( xr, lam )
    end if

    select case (p)

    case(1) ! P1

      phi = lam
      if ( present(dphi) ) then
        dphi = dlam
      end if

    case(2) ! P2

      phi(:,1) = lam(:,1) * ( 2*lam(:,1) - 1 )
      phi(:,3) = lam(:,2) * ( 2*lam(:,2) - 1 )
      phi(:,6) = lam(:,3) * ( 2*lam(:,3) - 1 )
      phi(:,2) = 4 * lam(:,1) * lam(:,2)
      phi(:,4) = 4 * lam(:,3) * lam(:,1)
      phi(:,5) = 4 * lam(:,2) * lam(:,3)

      if ( present(dphi) ) then

        do dm = 1, 2
          dphi(:,1,dm) = ( 4 * lam(:,1) - 1 ) * dlam(:,1,dm)
          dphi(:,3,dm) = ( 4 * lam(:,2) - 1 ) * dlam(:,2,dm)
          dphi(:,6,dm) = ( 4 * lam(:,3) - 1 ) * dlam(:,3,dm)
          dphi(:,2,dm) = 4 * ( dlam(:,1,dm) * lam(:,2) + &
                               lam(:,1) * dlam(:,2,dm) )
          dphi(:,5,dm) = 4 * ( dlam(:,2,dm) * lam(:,3) + &
                               lam(:,2) * dlam(:,3,dm) )
          dphi(:,4,dm) = 4 * ( dlam(:,3,dm) * lam(:,1) + &
                               lam(:,3) * dlam(:,1,dm) )
        end do

      end if

    case(3) ! P3

      phi(:,1)  = lam(:,1) * ( 3*lam(:,1) - 1 ) * ( 3*lam(:,1) - 2 ) / 2
      phi(:,4)  = lam(:,2) * ( 3*lam(:,2) - 1 ) * ( 3*lam(:,2) - 2 ) / 2
      phi(:,10) = lam(:,3) * ( 3*lam(:,3) - 1 ) * ( 3*lam(:,3) - 2 ) / 2
      phi(:,2)  = 9 * lam(:,1) * ( 3*lam(:,1) - 1 ) * lam(:,2) / 2
      phi(:,3)  = 9 * lam(:,2) * ( 3*lam(:,2) - 1 ) * lam(:,1) / 2
      phi(:,7)  = 9 * lam(:,2) * ( 3*lam(:,2) - 1 ) * lam(:,3) / 2
      phi(:,9)  = 9 * lam(:,3) * ( 3*lam(:,3) - 1 ) * lam(:,2) / 2
      phi(:,8)  = 9 * lam(:,3) * ( 3*lam(:,3) - 1 ) * lam(:,1) / 2
      phi(:,5)  = 9 * lam(:,1) * ( 3*lam(:,1) - 1 ) * lam(:,3) / 2
      phi(:,6)  = 27 * lam(:,1) * lam(:,2) * lam(:,3)

      if ( present(dphi) ) then

        do dm = 1, 2
          dphi(:,1,dm)  = ( 27*lam(:,1)**2 - 18*lam(:,1) + 2 ) / 2 * &
                            dlam(:,1,dm)
          dphi(:,4,dm)  = ( 27*lam(:,2)**2 - 18*lam(:,2) + 2 ) / 2 * &
                            dlam(:,2,dm)
          dphi(:,10,dm) = ( 27*lam(:,3)**2 - 18*lam(:,3) + 2 ) / 2 * &
                            dlam(:,3,dm)
          dphi(:,2,dm)  = 9 * ( ( 6*lam(:,1) - 1 ) * dlam(:,1,dm) * lam(:,2) + &
                             lam(:,1) * ( 3*lam(:,1) - 1 ) * dlam(:,2,dm) ) / 2
          dphi(:,3,dm)  = 9 * ( ( 6*lam(:,2) - 1 ) * dlam(:,2,dm) * lam(:,1) + &
                             lam(:,2) * ( 3*lam(:,2) - 1 ) * dlam(:,1,dm) ) / 2
          dphi(:,7,dm)  = 9 * ( ( 6*lam(:,2) - 1 ) * dlam(:,2,dm) * lam(:,3) + &
                             lam(:,2) * ( 3*lam(:,2) - 1 ) * dlam(:,3,dm) ) / 2
          dphi(:,9,dm)  = 9 * ( ( 6*lam(:,3) - 1 ) * dlam(:,3,dm) * lam(:,2) + &
                             lam(:,3) * ( 3*lam(:,3) - 1 ) * dlam(:,2,dm) ) / 2
          dphi(:,8,dm)  = 9 * ( ( 6*lam(:,3) - 1 ) * dlam(:,3,dm) * lam(:,1) + &
                             lam(:,3) * ( 3*lam(:,3) - 1 ) * dlam(:,1,dm) ) / 2
          dphi(:,5,dm)  = 9 * ( ( 6*lam(:,1) - 1 ) * dlam(:,1,dm) * lam(:,3) + &
                             lam(:,1) * ( 3*lam(:,1) - 1 ) * dlam(:,3,dm) ) / 2
          dphi(:,6,dm)  = 27 * ( dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                                 lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                                 lam(:,1)*lam(:,2)*dlam(:,3,dm) )

        end do

      end if

    case(4) ! P4

      phi(:,1)  = lam(:,1) * ( 4*lam(:,1) - 1 ) * ( 2*lam(:,1) - 1 ) * &
                   ( 4*lam(:,1) - 3 ) / 3
      phi(:,5)  = lam(:,2) * ( 4*lam(:,2) - 1 ) * ( 2*lam(:,2) - 1 ) * &
                   ( 4*lam(:,2) - 3 ) / 3
      phi(:,15) = lam(:,3) * ( 4*lam(:,3) - 1 ) * ( 2*lam(:,3) - 1 ) * &
                   ( 4*lam(:,3) - 3 ) / 3
      phi(:,2)  = 16 * lam(:,1) * ( 4*lam(:,1) - 1 ) * ( 2*lam(:,1) - 1 ) * &
                   lam(:,2) / 3
      phi(:,4)  = 16 * lam(:,2) * ( 4*lam(:,2) - 1 ) * ( 2*lam(:,2) - 1 ) * &
                   lam(:,1) / 3
      phi(:,9)  = 16 * lam(:,2) * ( 4*lam(:,2) - 1 ) * ( 2*lam(:,2) - 1 ) * &
                   lam(:,3) / 3
      phi(:,14) = 16 * lam(:,3) * ( 4*lam(:,3) - 1 ) * ( 2*lam(:,3) - 1 ) * &
                   lam(:,2) / 3
      phi(:,13) = 16 * lam(:,3) * ( 4*lam(:,3) - 1 ) * ( 2*lam(:,3) - 1 ) * &
                   lam(:,1) / 3
      phi(:,6)  = 16 * lam(:,1) * ( 4*lam(:,1) - 1 ) * ( 2*lam(:,1) - 1 ) * &
                   lam(:,3) / 3
      phi(:,3)  = 4 * lam(:,1) * ( 4*lam(:,1) - 1 ) * ( 4*lam(:,2) - 1 ) * &
                   lam(:,2)
      phi(:,12) = 4 * lam(:,2) * ( 4*lam(:,2) - 1 ) * ( 4*lam(:,3) - 1 ) * &
                   lam(:,3)
      phi(:,10) = 4 * lam(:,3) * ( 4*lam(:,3) - 1 ) * ( 4*lam(:,1) - 1 ) * &
                   lam(:,1)
      phi(:,7)  = 32 * ( 4*lam(:,1) - 1 ) * lam(:,1) * lam(:,2) * lam(:,3)
      phi(:,8)  = 32 * ( 4*lam(:,2) - 1 ) * lam(:,1) * lam(:,2) * lam(:,3)
      phi(:,11) = 32 * ( 4*lam(:,3) - 1 ) * lam(:,1) * lam(:,2) * lam(:,3)

      if ( present(dphi) ) then

        do dm = 1, 2
          dphi(:,1,dm)  = ( ( 8 * lam(:,1) - 1 ) * ( 2*lam(:,1) - 1 ) * &
                   ( 4*lam(:,1) - 3 ) +  lam(:,1) * ( 4*lam(:,1) - 1 ) * &
                   ( 16 * lam(:,1) - 10 ) ) / 3 * dlam(:,1,dm)
          dphi(:,5,dm)  = ( ( 8 * lam(:,2) - 1 ) * ( 2*lam(:,2) - 1 ) * &
                   ( 4*lam(:,2) - 3 ) +  lam(:,2) * ( 4*lam(:,2) - 1 ) * &
                   ( 16 * lam(:,2) - 10 ) ) / 3 * dlam(:,2,dm)
          dphi(:,15,dm) = ( ( 8 * lam(:,3) - 1 ) * ( 2*lam(:,3) - 1 ) * &
                   ( 4*lam(:,3) - 3 ) +  lam(:,3) * ( 4*lam(:,3) - 1 ) * &
                   ( 16 * lam(:,3) - 10 ) ) / 3 * dlam(:,3,dm)
          dphi(:,2,dm) = 16 * ( ( ( 8*lam(:,1) - 1 ) * ( 2*lam(:,1) - 1 ) + &
             lam(:,1) * ( 4*lam(:,1) - 1 ) * 2 ) * dlam(:,1,dm) * lam(:,2) + &
             lam(:,1) * ( 4*lam(:,1) - 1 ) * ( 2*lam(:,1) - 1 ) * &
                              dlam(:,2,dm) ) / 3
          dphi(:,4,dm) = 16 * ( ( ( 8*lam(:,2) - 1 ) * ( 2*lam(:,2) - 1 ) + &
             lam(:,2) * ( 4*lam(:,2) - 1 ) * 2 ) * dlam(:,2,dm) * lam(:,1) + &
             lam(:,2) * ( 4*lam(:,2) - 1 ) * ( 2*lam(:,2) - 1 ) * &
                              dlam(:,1,dm) ) / 3
          dphi(:,9,dm) = 16 * ( ( ( 8*lam(:,2) - 1 ) * ( 2*lam(:,2) - 1 ) + &
             lam(:,2) * ( 4*lam(:,2) - 1 ) * 2 ) * dlam(:,2,dm) * lam(:,3) + &
             lam(:,2) * ( 4*lam(:,2) - 1 ) * ( 2*lam(:,2) - 1 ) * &
                              dlam(:,3,dm) ) / 3
          dphi(:,14,dm) = 16 * ( ( ( 8*lam(:,3) - 1 ) * ( 2*lam(:,3) - 1 ) + &
             lam(:,3) * ( 4*lam(:,3) - 1 ) * 2 ) * dlam(:,3,dm) * lam(:,2) + &
             lam(:,3) * ( 4*lam(:,3) - 1 ) * ( 2*lam(:,3) - 1 ) * &
                              dlam(:,2,dm) ) / 3
          dphi(:,13,dm) = 16 * ( ( ( 8*lam(:,3) - 1 ) * ( 2*lam(:,3) - 1 ) + &
             lam(:,3) * ( 4*lam(:,3) - 1 ) * 2 ) * dlam(:,3,dm) * lam(:,1) + &
             lam(:,3) * ( 4*lam(:,3) - 1 ) * ( 2*lam(:,3) - 1 ) * &
                              dlam(:,1,dm) ) / 3
          dphi(:,6,dm) = 16 * ( ( ( 8*lam(:,1) - 1 ) * ( 2*lam(:,1) - 1 ) + &
             lam(:,1) * ( 4*lam(:,1) - 1 ) * 2 ) * dlam(:,1,dm) * lam(:,3) + &
             lam(:,1) * ( 4*lam(:,1) - 1 ) * ( 2*lam(:,1) - 1 ) * &
                              dlam(:,3,dm) ) / 3
          dphi(:,3,dm)  = 4 * ( ( 8*lam(:,1) - 1 ) * dlam(:,1,dm) &
                          * ( 4*lam(:,2) - 1 ) * lam(:,2) &
                          + lam(:,1) * ( 4*lam(:,1) - 1 ) * ( 8*lam(:,2) - 1 ) &
                          * dlam(:,2,dm) )
          dphi(:,12,dm) = 4 * ( ( 8*lam(:,2) - 1 ) * dlam(:,2,dm) &
                          * ( 4*lam(:,3) - 1 ) * lam(:,3) &
                          + lam(:,2) * ( 4*lam(:,2) - 1 ) * ( 8*lam(:,3) - 1 ) &
                          * dlam(:,3,dm) )
          dphi(:,10,dm) = 4 * ( ( 8*lam(:,3) - 1 ) * dlam(:,3,dm) &
                          * ( 4*lam(:,1) - 1 ) * lam(:,1) &
                          + lam(:,3) * ( 4*lam(:,3) - 1 ) * ( 8*lam(:,1) - 1 ) &
                          * dlam(:,1,dm) )
          dphi(:,7,dm)  = 32 * ( 4*dlam(:,1,dm) &
                          * lam(:,1) * lam(:,2) * lam(:,3) &
                  + ( 4*lam(:,1) - 1 ) * ( dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                                           lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                                           lam(:,1)*lam(:,2)*dlam(:,3,dm) ) )
          dphi(:,8,dm)  = 32 * ( 4*dlam(:,2,dm) &
                          * lam(:,1) * lam(:,2) * lam(:,3) &
                  + ( 4*lam(:,2) - 1 ) * ( dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                                           lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                                           lam(:,1)*lam(:,2)*dlam(:,3,dm) ) )
          dphi(:,11,dm) = 32 * ( 4*dlam(:,3,dm) &
                          * lam(:,1) * lam(:,2) * lam(:,3) &
                  + ( 4*lam(:,3) - 1 ) * ( dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                                           lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                                           lam(:,1)*lam(:,2)*dlam(:,3,dm) ) )
        end do

      end if

    case(5) ! P5

      phi(:,1)  = lam(:,1) * ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                   ( 5*lam(:,1) - 3 ) * ( 5*lam(:,1) - 4 ) / 24
      phi(:,6)  = lam(:,2) * ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                   ( 5*lam(:,2) - 3 ) * ( 5*lam(:,2) - 4 ) / 24
      phi(:,21) = lam(:,3) * ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                   ( 5*lam(:,3) - 3 ) * ( 5*lam(:,3) - 4 ) / 24
      phi(:,2)  = 25 * lam(:,1) * ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                   ( 5*lam(:,1) - 3 ) * lam(:,2) / 24
      phi(:,5)  = 25 * lam(:,2) * ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                   ( 5*lam(:,2) - 3 ) * lam(:,1) / 24
      phi(:,11) = 25 * lam(:,2) * ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                   ( 5*lam(:,2) - 3 ) * lam(:,3) / 24
      phi(:,20) = 25 * lam(:,3) * ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                   ( 5*lam(:,3) - 3 ) * lam(:,2) / 24
      phi(:,19) = 25 * lam(:,3) * ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                   ( 5*lam(:,3) - 3 ) * lam(:,1) / 24
      phi(:,7)  = 25 * lam(:,1) * ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                   ( 5*lam(:,1) - 3 ) * lam(:,3) / 24
      phi(:,3)  = 25 * lam(:,1) * ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                   ( 5*lam(:,2) - 1 ) * lam(:,2) / 12
      phi(:,4)  = 25 * lam(:,2) * ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                   ( 5*lam(:,1) - 1 ) * lam(:,1) / 12
      phi(:,15) = 25 * lam(:,2) * ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                   ( 5*lam(:,3) - 1 ) * lam(:,3) / 12
      phi(:,18) = 25 * lam(:,3) * ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                   ( 5*lam(:,2) - 1 ) * lam(:,2) / 12
      phi(:,16) = 25 * lam(:,3) * ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                   ( 5*lam(:,1) - 1 ) * lam(:,1) / 12
      phi(:,12) = 25 * lam(:,1) * ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                   ( 5*lam(:,3) - 1 ) * lam(:,3) / 12
      phi(:,8)  = 125 * ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                  lam(:,1) * lam(:,2) * lam(:,3) / 6
      phi(:,10) = 125 * ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                  lam(:,1) * lam(:,2) * lam(:,3) / 6
      phi(:,17) = 125 * ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                  lam(:,1) * lam(:,2) * lam(:,3) / 6
      phi(:,9)  = 125 * ( 5*lam(:,1) - 1 ) * ( 5*lam(:,2) - 1 ) * &
                  lam(:,1) * lam(:,2) * lam(:,3) / 4
      phi(:,14) = 125 * ( 5*lam(:,2) - 1 ) * ( 5*lam(:,3) - 1 ) * &
                  lam(:,1) * lam(:,2) * lam(:,3) / 4
      phi(:,13) = 125 * ( 5*lam(:,3) - 1 ) * ( 5*lam(:,1) - 1 ) * &
                  lam(:,1) * lam(:,2) * lam(:,3) / 4

      if ( present(dphi) ) then

        do dm = 1, 2
          dphi(:,1,dm)  = ( ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                            ( 5*lam(:,1) - 3 ) * ( 5*lam(:,1) - 4 ) + &
                              lam(:,1) * ( 50*lam(:,1) - 15 ) * &
                              ( 5*lam(:,1) - 3 ) * ( 5*lam(:,1) - 4 ) + &
                  lam(:,1) * ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                   ( 50*lam(:,1) - 35 ) ) * dlam(:,1,dm) / 24
          dphi(:,6,dm)  = ( ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                            ( 5*lam(:,2) - 3 ) * ( 5*lam(:,2) - 4 ) + &
                              lam(:,2) * ( 50*lam(:,2) - 15 ) * &
                              ( 5*lam(:,2) - 3 ) * ( 5*lam(:,2) - 4 ) + &
                  lam(:,2) * ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                   ( 50*lam(:,2) - 35 ) ) * dlam(:,2,dm) / 24
          dphi(:,21,dm) = ( ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                            ( 5*lam(:,3) - 3 ) * ( 5*lam(:,3) - 4 ) + &
                              lam(:,3) * ( 50*lam(:,3) - 15 ) * &
                              ( 5*lam(:,3) - 3 ) * ( 5*lam(:,3) - 4 ) + &
                  lam(:,3) * ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                   ( 50*lam(:,3) - 35 ) ) * dlam(:,3,dm) / 24
          dphi(:,2,dm)  = 25 * ( ( ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                   ( 5*lam(:,1) - 3 ) + lam(:,1) * ( 50*lam(:,1) - 15 ) * &
                   ( 5*lam(:,1) - 3 ) + lam(:,1) * ( 5*lam(:,1) - 1 ) * &
                   ( 5*lam(:,1) - 2 ) * 5 ) * dlam(:,1,dm) * lam(:,2) + &
                   lam(:,1) * ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                   ( 5*lam(:,1) - 3 ) * dlam(:,2,dm) ) / 24
          dphi(:,5,dm)  = 25 * ( ( ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                   ( 5*lam(:,2) - 3 ) + lam(:,2) * ( 50*lam(:,2) - 15 ) * &
                   ( 5*lam(:,2) - 3 ) + lam(:,2) * ( 5*lam(:,2) - 1 ) * &
                   ( 5*lam(:,2) - 2 ) * 5 ) * dlam(:,2,dm) * lam(:,1) + &
                   lam(:,2) * ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                   ( 5*lam(:,2) - 3 ) * dlam(:,1,dm) ) / 24
          dphi(:,11,dm)  = 25 * ( ( ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                   ( 5*lam(:,2) - 3 ) + lam(:,2) * ( 50*lam(:,2) - 15 ) * &
                   ( 5*lam(:,2) - 3 ) + lam(:,2) * ( 5*lam(:,2) - 1 ) * &
                   ( 5*lam(:,2) - 2 ) * 5 ) * dlam(:,2,dm) * lam(:,3) + &
                   lam(:,2) * ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                   ( 5*lam(:,2) - 3 ) * dlam(:,3,dm) ) / 24
          dphi(:,20,dm)  = 25 * ( ( ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                   ( 5*lam(:,3) - 3 ) + lam(:,3) * ( 50*lam(:,3) - 15 ) * &
                   ( 5*lam(:,3) - 3 ) + lam(:,3) * ( 5*lam(:,3) - 1 ) * &
                   ( 5*lam(:,3) - 2 ) * 5 ) * dlam(:,3,dm) * lam(:,2) + &
                   lam(:,3) * ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                   ( 5*lam(:,3) - 3 ) * dlam(:,2,dm) ) / 24
          dphi(:,19,dm)  = 25 * ( ( ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                   ( 5*lam(:,3) - 3 ) + lam(:,3) * ( 50*lam(:,3) - 15 ) * &
                   ( 5*lam(:,3) - 3 ) + lam(:,3) * ( 5*lam(:,3) - 1 ) * &
                   ( 5*lam(:,3) - 2 ) * 5 ) * dlam(:,3,dm) * lam(:,1) + &
                   lam(:,3) * ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                   ( 5*lam(:,3) - 3 ) * dlam(:,1,dm) ) / 24
          dphi(:,7,dm)  = 25 * ( ( ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                   ( 5*lam(:,1) - 3 ) + lam(:,1) * ( 50*lam(:,1) - 15 ) * &
                   ( 5*lam(:,1) - 3 ) + lam(:,1) * ( 5*lam(:,1) - 1 ) * &
                   ( 5*lam(:,1) - 2 ) * 5 ) * dlam(:,1,dm) * lam(:,3) + &
                   lam(:,1) * ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                   ( 5*lam(:,1) - 3 ) * dlam(:,3,dm) ) / 24
          dphi(:,3,dm)  = 25 * ( ( ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) + &
                   lam(:,1) * ( 50*lam(:,1) - 15 ) ) * dlam(:,1,dm) * &
                   ( 5*lam(:,2) - 1 ) * lam(:,2) + &
                   lam(:,1) * ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                   ( 10*lam(:,2) - 1 ) * dlam(:,2,dm) ) / 12
          dphi(:,4,dm)  = 25 * ( ( ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) + &
                   lam(:,2) * ( 50*lam(:,2) - 15 ) ) * dlam(:,2,dm) * &
                   ( 5*lam(:,1) - 1 ) * lam(:,1) + &
                   lam(:,2) * ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                   ( 10*lam(:,1) - 1 ) * dlam(:,1,dm) ) / 12
          dphi(:,15,dm) = 25 * ( ( ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) + &
                   lam(:,2) * ( 50*lam(:,2) - 15 ) ) * dlam(:,2,dm) * &
                   ( 5*lam(:,3) - 1 ) * lam(:,3) + &
                   lam(:,2) * ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                   ( 10*lam(:,3) - 1 ) * dlam(:,3,dm) ) / 12
          dphi(:,18,dm) = 25 * ( ( ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) + &
                   lam(:,3) * ( 50*lam(:,3) - 15 ) ) * dlam(:,3,dm) * &
                   ( 5*lam(:,2) - 1 ) * lam(:,2) + &
                   lam(:,3) * ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                   ( 10*lam(:,2) - 1 ) * dlam(:,2,dm) ) / 12
          dphi(:,16,dm) = 25 * ( ( ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) + &
                   lam(:,3) * ( 50*lam(:,3) - 15 ) ) * dlam(:,3,dm) * &
                   ( 5*lam(:,1) - 1 ) * lam(:,1) + &
                   lam(:,3) * ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                   ( 10*lam(:,1) - 1 ) * dlam(:,1,dm) ) / 12
          dphi(:,12,dm) = 25 * ( ( ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) + &
                   lam(:,1) * ( 50*lam(:,1) - 15 ) ) * dlam(:,1,dm) * &
                   ( 5*lam(:,3) - 1 ) * lam(:,3) + &
                   lam(:,1) * ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                   ( 10*lam(:,3) - 1 ) * dlam(:,3,dm) ) / 12
          dphi(:,8,dm)  = 125 * ( ( 50*lam(:,1) - 15 ) * dlam(:,1,dm) * &
                          lam(:,1) * lam(:,2) * lam(:,3) + &
                          ( 5*lam(:,1) - 1 ) * ( 5*lam(:,1) - 2 ) * &
                            ( dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                              lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                              lam(:,1)*lam(:,2)*dlam(:,3,dm) ) ) / 6
          dphi(:,10,dm) = 125 * ( ( 50*lam(:,2) - 15 ) * dlam(:,2,dm) * &
                          lam(:,1) * lam(:,2) * lam(:,3) + &
                          ( 5*lam(:,2) - 1 ) * ( 5*lam(:,2) - 2 ) * &
                            ( dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                              lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                              lam(:,1)*lam(:,2)*dlam(:,3,dm) ) ) / 6
          dphi(:,17,dm) = 125 * ( ( 50*lam(:,3) - 15 ) * dlam(:,3,dm) * &
                          lam(:,1) * lam(:,2) * lam(:,3) + &
                          ( 5*lam(:,3) - 1 ) * ( 5*lam(:,3) - 2 ) * &
                            ( dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                              lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                              lam(:,1)*lam(:,2)*dlam(:,3,dm) ) ) / 6
          dphi(:,9,dm)  = 125 * ( ( 5*dlam(:,1,dm) * ( 5*lam(:,2) - 1 ) + &
                          ( 5*lam(:,1) - 1 ) * 5*dlam(:,2,dm) ) * &
                          lam(:,1) * lam(:,2) * lam(:,3) + &
                          ( 5*lam(:,1) - 1 ) * ( 5*lam(:,2) - 1 ) * &
                          ( dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                            lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                            lam(:,1)*lam(:,2)*dlam(:,3,dm) ) ) / 4
          dphi(:,14,dm) = 125 * ( ( 5*dlam(:,2,dm) * ( 5*lam(:,3) - 1 ) + &
                          ( 5*lam(:,2) - 1 ) * 5*dlam(:,3,dm) ) * &
                          lam(:,1) * lam(:,2) * lam(:,3) + &
                          ( 5*lam(:,2) - 1 ) * ( 5*lam(:,3) - 1 ) * &
                          ( dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                            lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                            lam(:,1)*lam(:,2)*dlam(:,3,dm) ) ) / 4
          dphi(:,13,dm) = 125 * ( ( 5*dlam(:,3,dm) * ( 5*lam(:,1) - 1 ) + &
                          ( 5*lam(:,3) - 1 ) * 5*dlam(:,1,dm) ) * &
                          lam(:,1) * lam(:,2) * lam(:,3) + &
                          ( 5*lam(:,3) - 1 ) * ( 5*lam(:,1) - 1 ) * &
                          ( dlam(:,1,dm)*lam(:,2)*lam(:,3) + &
                            lam(:,1)*dlam(:,2,dm)*lam(:,3) + &
                            lam(:,1)*lam(:,2)*dlam(:,3,dm) ) ) / 4
        end do

      end if

    case(6:) ! P6 and higher

!     use generic routine

      call shape_triangle_lagrange_regular ( xr, p, phi, dphi )

    case default

      write(*,'(/a/a,i0/)') 'Error shape_triangle_Pp:', &
       ' order not available = ', p
      stop

    end select

  end subroutine shape_triangle_Pp


! Compute the value of the shapefunction and derivative at any reference point
! for a polynomial Pp. The shapefunction is based on Legendre polynomials and
! is discontinuous across element boundaries.
! Although the name says differently, this routine can be used on both
! triangles and quadrilaterals. On triangles the reference coordinates
! (xi,eta) \in lower triangle of [0:1,0:1] on one quarter of the functions
! on [-1:1,-1:1] is covered.
! A constant is subtracted from the higher-order even Legendre polynomials
! to achieve that in (xi,eta) = (0,0) only phi(:,1) (=1) is non-zero.

  subroutine shape_quad_Pp ( xr, p, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   polynomial order of polynomial
    integer, intent(in) :: p

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


    integer     :: i, j, k, m, np
    real(dp)    :: x, y
    real(dp), dimension(0:p) :: lx, ly, L0

    if ( p == 0 ) then

!     P0

      np = size(xr(:,1))
      phi(1:np,1) = 1._dp
      if ( present(dphi) ) dphi(1:np,1,:) = 0
      return

    else if ( present(dphi) ) then

      write(*,'(/a/)') &
        'Error in shape_quad_Pp: dphi not yet implemented for p > 1'
      stop

    end if

!   fill L0 (=Ln(0)) for even polynomial degrees

    L0(0) = 1
    do j = 2, p, 2
      L0(j) = -(j-1)*L0(j-2)/j
    end do

    np = size(xr(:,1))

    do i = 1, np

      x = xr(i,1)
      y = xr(i,2)

!     compute Legrende polynomials Ln (n=0,...,p) up to order p

      lx(0) = 1._dp
      ly(0) = 1._dp

      lx(1) = x
      ly(1) = y

      do j = 1, p-1
        lx(j+1) = ((2*j+1)*x*lx(j) - j*lx(j-1))/(j+1)
        ly(j+1) = ((2*j+1)*y*ly(j) - j*ly(j-1))/(j+1)
      end do

!     subtract Ln(0) for even order to get lx=ly=0 at the origin
!     of an element so that the first degree of freedom is the value on the
!     center of the element. Note, that orthogonality is lost.

      do j = 2, p, 2
        lx(j) = lx(j) - L0(j)
        ly(j) = ly(j) - L0(j)
      end do

!     P0

      phi(i,1) = 1._dp

!     P1,...,Pp

      m = 1
      do j = 1, p
        do k = 0, j
          phi(i,m+k+1) = lx(k) * ly(j-k)
        end do
        m = m + j + 1
      end do

    end do

  end subroutine shape_quad_Pp

end module shapefunc_highorder_m
