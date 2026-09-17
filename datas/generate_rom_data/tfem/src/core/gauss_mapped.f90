
! Copyright (C) 2011-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Mapped Gauss integration abscissas and weights (triangles and pyramids).

module gauss_mapped_m

  use kind_defs_m
  use gauss_defs_m
  use gauss_standard_numeric_m

  implicit none


contains

! Gauss Legendre in 2D defined on a triangle that is mapped from
! the region [-1,1]x[-1,1] (quad) by putting node (1,1) on top of
! (1,-1).

  subroutine Gauss_Legendre_mapped_triangle ( n, x, w )

!   the number of integration points in one dimension
    integer, intent(in) :: n

!   the reference coordinates and weights of the integration points
!   x(ninti,ndim), w(ninti)
!   where ninti = n**2
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w

    integer :: m
    real(dp) :: tmp(size(x,1))

!   the reference region of a triangle is
!
!      1 |\               eta \in [0,1]
!      ^ |  \             xi  \in [0,1]
!  eta | |    \           xi + eta <= 1
!        |      \     !
!      0 ---------
!        0  xi -> 1
!

    m = n**2 ! number of points

!   test

    if ( size(x,1) < m .or. size(w) < m ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_mapped_triangle: size of x or w too small'
      stop
    end if

!   integration on reference quad [-1,1]x[-1,1]

    call Gauss_Legendre_numeric_quad ( n, x, w )

!   modify weight by multiplying with Jacobian

    w(:m) = w(:m) * ( 1 - x(:m,1) ) / 8

!   map coordinates from quad to triangle

    tmp(:m) = ( 1 + x(:m,1) ) / 2
    x(:m,2) = ( 1 - x(:m,1) ) * ( 1 + x(:m,2) ) / 4
    x(:m,1) = tmp(:m)

  end subroutine Gauss_Legendre_mapped_triangle


! Gauss Legendre in 3D defined on a pyramid that is mapped from
! the region [-1,1]x[-1,1]x[-1,1] (unit cube) to the pyramid with
! base [-1,1]x[-1,1] by merging the upper four corners into one apex at (0,0,1).

  subroutine Gauss_Legendre_mapped_pyramid ( n, n2, x, w )

!   the number of integration points in one dimension on bottom quad
    integer, intent(in) :: n

!   the number of integration points in one dimension in the height direction
    integer, intent(in) :: n2

!   the reference coordinates and weights of the integration points
!   x(ninti,ndim), w(ninti)
!   where ninti = n**2*n2
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w

    integer :: ipk, ip, i, j, k
    real(dp), dimension(n) :: x1, w1
    real(dp), dimension(n2) :: x2, w2

!   The reference region of a pyramid is
!
!
!           /\\             !
!          /| \ \           !
!         / /  \  \         !
!        / |    \   \       !
!       |  |     |    \     !
!       |  ------\------
!      / /        |    /
!     | /         \   /
!     |/           | /
!     --------------/
!
!          zeta
!            ^ _ eta
!            | /|
!            |/
!            ----> xi
!
!   The reference coordinates xr=(xi,eta,zeta) are
!   in the region [-1,1]x[-1,1]x[0,1].
!   Notes: zeta=0 is the bottom square, zeta=1 is the apex,
!          |xi| <= 1-zeta, |eta| <= 1-zeta

!   test

    if ( size(x,1) < n**2*n2 .or. size(w) < n**2*n2 ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_mapped_pyramid: size of x or w too small'
      stop
    end if

!   integration on reference line [-1,1]

    call Gauss_Legendre_numeric_line ( n, x1, w1 )
    call Gauss_Legendre_numeric_line ( n2, x2, w2 )

!   map coordinates and weights

!   map height from [-1,1] to [0,1]
    x2 = ( 1 + x2 ) / 2  ! zeta
    w2 = w2 / 2

!   Jacobian of the mapping
    w2 = w2 * ( 1 - x2 ) ** 2

    do k = 1, n2
      ipk = n**2 * ( k - 1)
      do i = 1, n
        do j = 1, n
          ip = i + n * ( j - 1 ) + ipk
          x(ip,1) = x1(i) * ( 1 - x2(k) ) ! mapping of xi
          x(ip,2) = x1(j) * ( 1 - x2(k) ) ! mapping of eta
          x(ip,3) = x2(k)
          w(ip)   = w1(i) * w1(j) * w2(k)
        end do
      end do
    end do

  end subroutine Gauss_Legendre_mapped_pyramid


! Set number of Gauss integration points for an element

  subroutine set_ninti_mapped ( gauss, ninti )

    type(gauss_t), intent(in) :: gauss

!   number of integration points
    integer, intent(out) :: ninti


!   set integration

    if ( gauss%globalshape == 'triangle' ) then

      ninti = gauss%intrule ** 2

    else if ( gauss%globalshape == 'pyramid' ) then

      ninti = gauss%intrule ** 2 * gauss%intrule2

    else

      write(*,'(/2a/)') &
        'Error in set_ninti_mapped: not available for ', gauss%globalshape
      stop

    end if

  end subroutine set_ninti_mapped

end module gauss_mapped_m
