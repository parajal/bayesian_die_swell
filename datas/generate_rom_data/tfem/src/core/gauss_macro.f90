
! Copyright (C) 2007-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Gauss integration abscissas and weights for elements divided into smaller
! subdomains (macroelements)

module gauss_macro_m

  use kind_defs_m
  use gauss_defs_m
  use gauss_standard_m

  implicit none

contains

! Generalized Gauss Legendre in 1D defined on the interval [-1,1] using
! m equally sized subdomains

  subroutine Gauss_Legendre_line_macro ( m, xs, ws, x, w )

!   number of subdomains on [-1,1]
    integer, intent(in) :: m

!   the reference coordinates and weights of the integration points
!   on a single subdomain [-1,1], xs(i), ws(i), i=1,n
!   Note: n = size(ws)
    real(dp), intent(in), dimension(:) :: xs, ws

!   the reference coordinates and weights of the integration points
!   x(i), w(i), i=1,ni
!   where ni=n*m
    real(dp), intent(out), dimension(:) :: x, w

    integer :: n, ni, i

    n = size(ws)
    ni = n*m

!   test

    if ( size(x) < ni .or. size(w) < ni ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_line_macro: size of x or w too small'
      stop
    end if

!   Fill x and w

    x(1:ni) = [ ( ( xs + 2 * i - m - 1 ) / m, i = 1, m ) ]

    w(1:ni) = [ ( ws / m, i = 1, m ) ]

  end subroutine Gauss_Legendre_line_macro


! Generalized Gauss Legendre in 2D defined on the region [-1,1]x[-1,1] (quad)
! using m x m equally sized subdomains

  subroutine Gauss_Legendre_quad_macro ( m, xs, ws, x, w )

!   number of subdomains on [-1,1], m x m on quad [-1,1]x[-1,1]
    integer, intent(in) :: m

!   the reference coordinates and weights of the integration points
!   on a single subdomain (xs(i,j),j=1,2), ws(i), i=1,n
!   Note: n = size(ws)
    real(dp), intent(in), dimension(:,:) :: xs
    real(dp), intent(in), dimension(:) :: ws

!   the reference coordinates and weights of the integration points
!   x(i,1:ndim), w(i), i = 1, ni
!   where ni = n*m**2
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w

    real(dp) :: d
    integer :: n, ni, i, j, mt, l
    real(dp), dimension(:,:), allocatable :: xsl
    real(dp), dimension(:), allocatable :: wsl

    n = size(ws)
    mt = m**2
    ni = n*mt

!   test

    if ( size(x) < ni .or. size(w) < ni ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_quad_macro: size of x or w too small'
      stop
    end if

    xsl = ( xs + 1 ) / m - 1
    wsl = ws / mt
    d = 2._dp / m

!   Fill x and w

    l = 0
    do j = 0, m-1
      do i = 0, m-1
        x(l+1:l+n,1) = xsl(:,1) + i * d
        x(l+1:l+n,2) = xsl(:,2) + j * d
        l = l + n
      end do
    end do

    w(1:ni) = [ ( wsl, j = 1, mt ) ]

  end subroutine Gauss_Legendre_quad_macro


! Generalized Gauss Legendre in 2D defined on the triangular reference region
! Division into m^2 equally sized subtriangles

  subroutine Gauss_Legendre_triangle_macro ( m, xs, ws, x, w )

!   number of subtriangles in one dimension
    integer, intent(in) :: m

!   the reference coordinates and weights of the integration points
!   on a single subdomain (xs(i,j),j=1,2), ws(i), i=1,n
!   Note: n = size(ws)
    real(dp), intent(in), dimension(:,:) :: xs
    real(dp), intent(in), dimension(:) :: ws

!   the reference coordinates and weights of the integration points
!   x(ninti,ndim), w(ninti)
!   where ninti = n*m**2
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w


!   the reference region of a triangle is
!
!      1 |\               eta \in [0,1]
!      ^ |  \             xi  \in [0,1]
!  eta | |    \           xi + eta <= 1
!        |      \      !
!      0 ---------
!        0  xi -> 1
!
!   the reference triangle has been divided into m**2 subtriangle. On each
!   triangle a n Gauss rule is applied. For example with n=1, m=2 we have
!   4 subtriangles where in each subtriangle there is one integration point in
!   the center of gravity of the subtriangle.
!
!     | \             !
!     | x \           !
!     |----\          !
!     | \ x| \        !
!     |x \ |x \       !
!     |--------
!

    integer :: n, ni, i, j, nlow, nup, mt
    real(dp), dimension(:), allocatable :: wsl

    n = size(ws)
    mt = m ** 2
    ni = n * mt

!   test

    if ( size(x,1) < ni .or. size(w) < ni ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_triangle_macro: size of x or w too small'
      stop
    end if

    nlow = n * m * ( m + 1 ) / 2
    nup  = n * m * ( m - 1 ) / 2

!   Fill x and w

!   "lower triangles"

    x(1:nlow,1) = [ ( ( ( xs(:,1) + i - 1 ) / m, i = 1, m-j+1 ), j = 1, m ) ]
    x(1:nlow,2) = [ ( ( ( xs(:,2) + j - 1 ) / m, i = 1, m-j+1 ), j = 1, m ) ]

!   "upper triangles"

    x(nlow+1:nlow+nup,1) = &
                    [ ( ( ( - xs(:,1) + i ) / m, i = 1, m-j ), j = 1, m-1 ) ]
    x(nlow+1:nlow+nup,2) = &
                    [ ( ( ( - xs(:,2) + j ) / m, i = 1, m-j ), j = 1, m-1 ) ]

    wsl = ws / mt

    w(1:ni) = [ ( wsl, j = 1, mt ) ]

  end subroutine Gauss_Legendre_triangle_macro


! Generalized Gauss Legendre in 3D defined on the region [-1,1]x[-1,1]x[-1,1]
! using m x m x m equally sized subdomains

  subroutine Gauss_Legendre_hexahedron_macro ( m, xs, ws, x, w )

!   number of subdomains on [-1,1], m x m x m on hexahedron [-1,1]x[-1,1]x[-1,1]
    integer, intent(in) :: m

!   the reference coordinates and weights of the integration points
!   on a single subdomain (xs(i,j),j=1,3), ws(i), i=1,n
!   Note: n = size(ws)
    real(dp), intent(in), dimension(:,:) :: xs
    real(dp), intent(in), dimension(:) :: ws

!   the reference coordinates and weights of the integration points
!   x(i,1:ndim), w(i), i = 1, ni
!   where ni = n*m**3
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w

    real(dp) :: d
    integer :: n, ni, i, j, k, mt, l
    real(dp), dimension(:,:), allocatable :: xsl
    real(dp), dimension(:), allocatable :: wsl

    n = size(ws)
    mt = m**3
    ni = n*mt

!   test

    if ( size(x) < ni .or. size(w) < ni ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_hexahedron_macro: size of x or w too small'
      stop
    end if

    xsl = ( xs + 1 ) / m - 1
    wsl = ws / mt
    d = 2._dp / m

!   Fill x and w

    l = 0
    do k = 0, m-1
      do j = 0, m-1
        do i = 0, m-1
          x(l+1:l+n,1) = xsl(:,1) + i * d
          x(l+1:l+n,2) = xsl(:,2) + j * d
          x(l+1:l+n,3) = xsl(:,3) + k * d
          l = l + n
        end do
      end do
    end do

    w(1:ni) = [ ( wsl, j = 1, mt ) ]

  end subroutine Gauss_Legendre_hexahedron_macro


! Generalized Gauss Legendre in 3D defined on the tetrahedral reference region
! Division into m^2 equally sized subtetrahedrons.

  subroutine Gauss_Legendre_tetrahedron_macro ( m, xs, ws, x, w )

!   number of subtriangles in one dimension
    integer, intent(in) :: m

!   the reference coordinates and weights of the integration points
!   on a single subdomain (xs(i,j),j=1,3), ws(i), i=1,n
!   Note: n = size(ws)
    real(dp), intent(in), dimension(:,:) :: xs
    real(dp), intent(in), dimension(:) :: ws

!   the reference coordinates and weights of the integration points
!   x(i,1:ndim), w(i), i = 1, ni
!   where ni = n*m**3
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w

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
!   the reference tetrahedron has been divided into m**3 subtetrahedrons.
!   On each tetrahedron a n Gauss rule is applied. For example with n=1, m=2
!   we have 8 subtetrahedrons where in each subtetrahedron there is one
!   integration point in the center of gravity of the subtetrahedron.
!
!   NOTES:
!   - Currently only m=2 is available. For m>2 a recursive routine needs to
!     be implemented.
!   - The subdivision is not unique and the subdidision used for P1isoP2 shape
!     function has been chosen to make the integration consistent with the
!     piecewise linear shape function.

    integer :: ni, i, j, k, is(10,4), n, mt
    real(dp) :: p(10,3), phi(size(ws),4)
    real(dp), dimension(:), allocatable :: wsl


    if ( m /= 2 ) then
      write(*,'(/a,i0/a/)') &
        'Error in Gauss_Legendre_tetrahedron_macro: subdivision m = ', m, &
        ' not implemented. Currently only m=2 is available.'
      stop
    end if

    n = size(ws)
    mt = m**3
    ni = n*mt

!   test

    if ( size(x,1) < ni .or. size(w) < ni ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_tetrahedron_macro: size of x or w too small'
      stop
    end if

!   Shape function on a reference tetrahedron

    phi(:,1) = 1 - xs(:,1) - xs(:,2) - xs(:,3)
    phi(:,2) = xs(:,1)
    phi(:,3) = xs(:,2)
    phi(:,4) = xs(:,3)

!   Nodes

    p(1,:) = [ 0, 0, 0 ]
    p(2,:) = [ 0.5_dp, 0._dp, 0._dp ]
    p(3,:) = [ 1, 0, 0 ]
    p(4,:) = [ 0.5_dp, 0.5_dp, 0._dp ]
    p(5,:) = [ 0, 1, 0 ]
    p(6,:) = [ 0._dp, 0.5_dp, 0._dp ]
    p(7,:) = [ 0._dp, 0._dp, 0.5_dp ]
    p(8,:) = [ 0.5_dp, 0._dp, 0.5_dp ]
    p(9,:) = [ 0._dp, 0.5_dp, 0.5_dp ]
    p(10,:) = [ 0, 0, 1 ]

!   Topology of subelements

    is(1,:) = [ 1, 2, 6, 7 ]
    is(2,:) = [ 2, 3, 4, 8 ]
    is(3,:) = [ 4, 5, 6, 9 ]
    is(4,:) = [ 7, 8, 9, 10 ]
    is(5,:) = [ 2, 8, 9, 7 ]
    is(6,:) = [ 2, 9, 6, 7 ]
    is(7,:) = [ 2, 9, 8, 4 ]
    is(8,:) = [ 2, 6, 9, 4 ]

!   Fill x and w

    k = 0
    do i = 1, 8
      x(k+1:k+n,:) = matmul ( phi, p(is(i,:),:) )
      k = k + n
    end do

    wsl = ws / mt

    w(1:ni) = [ ( wsl, j = 1, mt ) ]

  end subroutine Gauss_Legendre_tetrahedron_macro


! Set number of Gauss integration points for an element

  subroutine set_ninti_macro ( gauss, n, ninti )

    type(gauss_t), intent(in) :: gauss

!   number of integration points of single subdomain
    integer, intent(in) :: n

!   number of integration points
    integer, intent(out) :: ninti

    integer :: m

!   local parameters

    m = gauss%nsubint

!   set integration

    if ( gauss%globalshape == 'line' ) then

      ninti = n * m

    else if ( gauss%globalshape == 'triangle' .or. &
              gauss%globalshape == 'quadrilateral' ) then

      ninti = n * m ** 2

    else if ( gauss%globalshape == 'tetrahedron' .or. &
              gauss%globalshape == 'hexahedron' ) then

      ninti = n * m ** 3

    else

      write(*,'(/2a/)') &
        'Error in set_ninti_macro: not available for ', &
        gauss%globalshape
      stop

    end if

  end subroutine set_ninti_macro

end module gauss_macro_m
