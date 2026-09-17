
! Copyright (C) 2004-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Standard Gauss integration abscissas and weights using analytical expressions

module gauss_standard_m

  use kind_defs_m
  use gauss_defs_m

  implicit none


contains


! Gauss Legendre in 1D defined on the interval [-1,1]

  subroutine Gauss_Legendre_line ( n, x, w )

!   number of integration points
    integer, intent(in) :: n

!   the reference coordinates and weights of the integration points
!   x(i), w(i), i=1,n
    real(dp), intent(out), dimension(:) :: x, w

    if ( size(x) < n .or. size(w) < n ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_line: size of x or w too small'
      stop
    end if

    select case(n)
      case(1)
        x(1) = 0
        w(1) = 2
      case(2)
        x(2) = sqrt(3._dp)/3
        x(1) = -x(2)
        w(1) = 1
        w(2) = 1
      case(3)
        x(3) = sqrt(15._dp)/5
        x(2) = 0
        x(1) = -x(3)
        w(1) = 5._dp/9
        w(2) = 8._dp/9
        w(3) = w(1)
      case(4)
        x(4) = sqrt(525+70*sqrt(30._dp))/35
        x(3) = sqrt(525-70*sqrt(30._dp))/35
        x(2) = -x(3)
        x(1) = -x(4)
        w(1) = (18-sqrt(30._dp))/36
        w(2) = (18+sqrt(30._dp))/36
        w(3) = w(2)
        w(4) = w(1)
      case(5)
        x(5) = sqrt(245+14*sqrt(70._dp))/21
        x(4) = sqrt(245-14*sqrt(70._dp))/21
        x(3) = 0
        x(2) = -x(4)
        x(1) = -x(5)
        w(1) = (322-13*sqrt(70._dp))/900
        w(2) = (322+13*sqrt(70._dp))/900
        w(3) = 128._dp/225
        w(4) = w(2)
        w(5) = w(1)
      case default
        write(*,'(/a,i0/)') &
          'Error in Gauss_Legendre_line: invalid number of points, n = ', n
        stop
    end select

  end subroutine Gauss_Legendre_line


! Gauss Legendre in 2D defined on the region [-1,1]x[-1,1] (quad)

  subroutine Gauss_Legendre_quad ( n, x, w )

!   the number of integration points in one dimension
    integer, intent(in) :: n

!   the reference coordinates and weights of the integration points
!   x(ninti,ndim), w(ninti)
!   where ninti = n**2
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w


    integer :: i, j, ip
    real(dp), dimension(n) :: x1, w1


!   test

    if ( size(x,1) < n**2 .or. size(w) < n**2 ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_quad: size of x or w too small'
      stop
    end if


    call Gauss_Legendre_line ( n, x1, w1 )

    do i = 1, n
      do j = 1, n
        ip = i + n * ( j - 1 )
        x(ip,1) = x1(i)
        x(ip,2) = x1(j)
        w(ip)   = w1(i) * w1(j)
      end do
    end do

  end subroutine Gauss_Legendre_quad


! Gauss Legendre in 2D defined on triangular reference region

  subroutine Gauss_Legendre_triangle ( n, x, w )

!   the number of integration points
    integer, intent(in) :: n

!   the reference coordinates and weights of the integration points
!   x(ninti,ndim), w(ninti)
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w


!   the reference region of a triangle is
!
!      1 |\               eta \in [0,1]
!      ^ |  \             xi  \in [0,1]
!  eta | |    \           xi + eta <= 1
!        |      \     !
!      0 ---------
!        0  xi -> 1
!

    real(dp) :: sqrt15, t, A, r, u, s, v, B, C, sqrt10, sqrtsqrt
    real(dp) :: g1, g2, g3, g4, g5, h1, h2, w1, w2, w3

!   test

    if ( size(x,1) < n .or. size(w) < n ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_triangle: size of x or w too small'
      stop
    end if

    select case (n)
      case(1)

!       first order precision, one-point Gauss (midpoint)

        t  = 1._dp / 3
        x(1,1:2) = [ t, t ]

        w(1) = .5_dp

      case(3)

!       second order precision, three-point Gauss

        t  = 1._dp / 6
        r  = 2._dp / 3
        x(1,1:2) = [ t, t ]
        x(2,1:2) = [ r, t ]
        x(3,1:2) = [ t, r ]

        w(1:3) = t

      case(6)

!       fourth order precision, six-point Gauss (Cowper, 1973)
!       (used the data from the online lecture notes of Carlos Felippa:
!        http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/.
!        see also: Felippa, C.A. (2004). A compendium of FEM integration
!        formulas for symbolic work. Engineering Computations, 21(8), 867–890.

        sqrt10 = sqrt(10._dp)
        sqrtsqrt = sqrt(38-44*sqrt(0.4_dp))

        g1=(8-sqrt10+sqrtsqrt)/18
        g2=(8-sqrt10-sqrtsqrt)/18
        h1 = 1-2*g1
        h2 = 1-2*g2
        x(1,1:2) = [ g1, g1 ]
        x(2,1:2) = [ h1, g1 ]
        x(3,1:2) = [ g1, h1 ]
        x(4,1:2) = [ g2, g2 ]
        x(5,1:2) = [ h2, g2 ]
        x(6,1:2) = [ g2, h2 ]

        sqrtsqrt = sqrt(213125-53320*sqrt10)

        w(1:3) = (620+sqrtsqrt)/7440
        w(4:6) = (620-sqrtsqrt)/7440

      case(7)

!       fifth order precision, seven-point Gauss (Hammer, Marlowe and Stroud)

        sqrt15 = sqrt(15._dp)

        t  = 1._dp / 3
        A  = 9._dp / 80
        r = ( 6 - sqrt15 ) / 21
        u = ( 6 + sqrt15 ) / 21
        s = ( 9 + 2*sqrt15 ) / 21
        v = ( 9 - 2*sqrt15 ) / 21
        B = ( 155 - sqrt15 ) / 2400
        C = ( 155 + sqrt15 ) / 2400

        x(1,1:2) = [ t, t ]
        x(2,1:2) = [ r, r ]
        x(3,1:2) = [ r, s ]
        x(4,1:2) = [ s, r ]
        x(5,1:2) = [ u, u ]
        x(6,1:2) = [ u, v ]
        x(7,1:2) = [ v, u ]

        w(1) = A
        w(2:4) = B
        w(5:7) = C

      case(12)

!       sixth order precision, 12-point Gauss
!       (used the data from the online lecture notes of Carlos Felippa:
!        http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/
!        see also: Felippa, C.A. (2004). A compendium of FEM integration
!        formulas for symbolic work. Engineering Computations, 21(8), 867–890.)

        g1=0.063089014491502228340331602870819157_dp
        g2=0.249286745170910421291638553107019076_dp
        g3=0.053145049844816947353249671631398147_dp
        g4=0.310352451033784405416607733956552153_dp
        w1=(30*g2**3*(4*g3**2+(1-2*g4)**2+4*g3*(-1+g4))+ &
            g3**2*(1-15*g4)+(-1+g4)*g4-g3*(-1+g4)*(-1+15*g4)+ &
            2*g2*(1+60*g3*g4*(-1+g3+g4))-6*g2**2*(3+10*(-1+g4)*g4+ &
            10*g3**2*(1+3*g4)+10*g3*(-1+g4)*(1+3*g4)))/ &
           (180*(g1-g2)*(-(g2*(-1+2*g2)*(-1+g3)*g3)+(-1+g3)*(g2-2*g2**2- &
            2*g3+3*g2*g3)*g4-(g2*(-1+2*g2-3*g3)+2*g3)*g4**2+ &
            2*g1**2*(g2*(-2+3*g2)+g3-g3**2+g4-g3*g4-g4**2)+ &
            g1*(-4*g2**2+(-1+g3)*g3+(-1+g3)*(1+3*g3)*g4+(1+3*g3)* &
            g4**2-2*g2*(-1+g3**2+g3*(-1+g4)+(-1+g4)*g4))))
        w2=(-1+12*(2-3*g1)*g1*w1+4*g3**2*(-1+3*w1)+4*g3*(-1+g4)*(-1+3*w1)+ &
            4*(-1+g4)*g4*(-1+3*w1))/(12*(g2*(-2+3*g2)+g3-g3**2+g4-g3*g4-g4**2))
        w3=(1-3*w1-3*w2)/6

        h1 = 1-2*g1
        h2 = 1-2*g2
        g5 = 1-g3-g4

        x(1,1:2) = [ g1, g1 ]
        x(2,1:2) = [ h1, g1 ]
        x(3,1:2) = [ g1, h1 ]
        x(4,1:2) = [ g2, g2 ]
        x(5,1:2) = [ h2, g2 ]
        x(6,1:2) = [ g2, h2 ]
        x(7,1:2) = [ g3, g4 ]
        x(8,1:2) = [ g4, g3 ]
        x(9,1:2) = [ g5, g3 ]
        x(10,1:2) = [ g5, g4 ]
        x(11,1:2) = [ g3, g5 ]
        x(12,1:2) = [ g4, g5 ]

        w(1:3) = w1 / 2
        w(4:6) = w2 / 2
        w(7:12) = w3 / 2

      case default
        write(*,'(/a,i0/)') &
          'Error in Gauss_Legendre_triangle: invalid number of points, n = ', n
        stop
    end select

  end subroutine Gauss_Legendre_triangle


! Gauss Legendre in 3D defined on the region [-1,1]x[-1,1]x[-1,1] (hexahedron)

  subroutine Gauss_Legendre_hexahedron ( n, x, w )

!   the number of integration points in one dimension
    integer, intent(in) :: n

!   the reference coordinates and weights of the integration points
!   x(ninti,ndim), w(ninti)
!   where ninti = n**3
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w


    integer :: i, j, k, ip
    real(dp), dimension(n) :: x1, w1


!   test

    if ( size(x,1) < n**3 .or. size(w) < n**3 ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_hexahedron: size of x or w too small'
      stop
    end if


    call Gauss_Legendre_line ( n, x1, w1 )

    do i = 1, n
      do j = 1, n
        do k = 1, n
          ip = i + n * ( j - 1 ) + n**2 * ( k - 1 )
          x(ip,1) = x1(i)
          x(ip,2) = x1(j)
          x(ip,3) = x1(k)
          w(ip)   = w1(i) * w1(j) * w1(k)
        end do
      end do
    end do

  end subroutine Gauss_Legendre_hexahedron


! Gauss Legendre in 3D defined on a tetrahedral reference region

  subroutine Gauss_Legendre_tetrahedron ( n, x, w )

!   the number of integration points
    integer, intent(in) :: n

!   the reference coordinates and weights of the integration points
!   x(ninti,ndim), w(ninti)
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

    real(dp) :: sqrt15, r, s, A, s1, s2, t1, t2, B1, B2, u, v, C
    real(dp) :: sqrt5, g1, g2, g3, g4, h1, h2, h3, h4, p4, w1, w2, w3
    real(dp) :: sqrt17

!   test

    if ( size(x,1) < n .or. size(w) < n ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_tetrahedron: size of x or w too small'
      stop
    end if

    select case (n)
      case(1)

!       first order precision, one point (midpoint)

        x(1,1:3) = [ 0.25_dp, 0.25_dp, 0.25_dp ]

        w(1) = 1._dp / 6

      case(4)

!       second order precision, four-point Gauss (Hammer & Stroud)

        r = ( 5 - sqrt(5._dp) ) / 20
        s = ( 5 + 3 * sqrt(5._dp) ) / 20

        x(1,1:3) = [ r, r, r ]
        x(2,1:3) = [ s, r, r ]
        x(3,1:3) = [ r, s, r ]
        x(4,1:3) = [ r, r, s ]

        w(1:4) = 1._dp / 24

      case(8)

!       third order precision, 8-point Gauss
!       (used the data from the online lecture notes of Carlos Felippa:
!        http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/
!        see also: Felippa, C.A. (2004). A compendium of FEM integration
!        formulas for symbolic work. Engineering Computations, 21(8), 867–890.)

        sqrt17 = sqrt(17._dp)

        g1=(55-3*sqrt17+sqrt(1022-134*sqrt17))/196
        g2=(55-3*sqrt17-sqrt(1022-134*sqrt17))/196
        w1=1._dp/8+sqrt((1715161837-406006699*sqrt17)/23101)/3120
        w2=1._dp/8-sqrt((1715161837-406006699*sqrt17)/23101)/3120

        h1=1-3*g1
        h2=1-3*g2

        x(1,1:3) = [ g1, g1, g1 ]
        x(2,1:3) = [ h1, g1, g1 ]
        x(3,1:3) = [ g1, h1, g1 ]
        x(4,1:3) = [ g1, g1, h1 ]
        x(5,1:3) = [ g2, g2, g2 ]
        x(6,1:3) = [ h2, g2, g2 ]
        x(7,1:3) = [ g2, h2, g2 ]
        x(8,1:3) = [ g2, g2, h2 ]

        w(1:4) = w1 / 6
        w(5:8) = w2 / 6

      case(15)

!       fifth order precision, fifteen-point Gauss (Stroud)

        sqrt15 = sqrt(15._dp)

        r  = 0.25_dp
        A  = 16._dp / 810
        s1 = ( 7 - sqrt15 ) / 34
        s2 = ( 7 + sqrt15 ) / 34
        t1 = ( 13 + 3*sqrt15 ) / 34
        t2 = ( 13 - 3*sqrt15 ) / 34
        B1 = ( 2665 + 14*sqrt15 ) / 226800
        B2 = ( 2665 - 14*sqrt15 ) / 226800
        u  = ( 10 - 2*sqrt15 ) / 40
        v  = ( 10 + 2*sqrt15 ) / 40
        C  = 20._dp / 2268

        x(1,1:3) = [ r, r, r ]
        x(2,1:3) = [ s1, s1, s1 ]
        x(3,1:3) = [ t1, s1, s1 ]
        x(4,1:3) = [ s1, t1, s1 ]
        x(5,1:3) = [ s1, s1, t1 ]
        x(6,1:3) = [ s2, s2, s2 ]
        x(7,1:3) = [ t2, s2, s2 ]
        x(8,1:3) = [ s2, t2, s2 ]
        x(9,1:3) = [ s2, s2, t2 ]
        x(10,1:3) = [ u, u, v ]
        x(11,1:3) = [ u, v, u ]
        x(12,1:3) = [ v, u, u ]
        x(13,1:3) = [ u, v, v ]
        x(14,1:3) = [ v, u, v ]
        x(15,1:3) = [ v, v, u ]

        w(1) = A
        w(2:5) = B1
        w(6:9) = B2
        w(10:15) = C

      case(24)

!       sixth order precision, 24-point Gauss
!       (used the data from the online lecture notes of Carlos Felippa:
!        http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/
!        see also: Felippa, C.A. (2004). A compendium of FEM integration
!        formulas for symbolic work. Engineering Computations, 21(8), 867–890.)

        sqrt5 = sqrt(5._dp)

        g1 = 0.214602871259152029288839219386284991_dp
        g2 = 0.040673958534611353115579448956410059_dp
        g3 = 0.322337890142275510343994470762492125_dp
        g4 = (3-sqrt5)/12
        h4 = (5+sqrt5)/12
        p4 = (1+sqrt5)/12
        w1 = (85+2*g2*(-319+9*sqrt5+624*g2)-638*g3- &
              24*g2*(-229+472*g2)*g3+96*(13+118*g2*(-1+2*g2))*g3**2+ &
              9*sqrt5*(-1+2*g3))/(13440*(g1-g2)*(g1-g3)*(3-8*g2+ &
              8*g1*(-1+2*g2)-8*g3+16*(g1+g2)*g3))
        w2 = -(85+2*g1*(-319+9*sqrt5+624*g1)-638*g3- &
              24*g1*(-229+472*g1)*g3+96*(13+118*g1*(-1+2*g1))*g3**2+ &
              9*sqrt5*(-1+2*g3))/(13440*(g1-g2)*(g2-g3)*(3-8*g2+ &
              8*g1*(-1+2*g2)-8*g3+16*(g1+g2)*g3))
        w3 = (85+2*g1*(-319+9*sqrt5+624*g1)-638*g2- &
              24*g1*(-229+472*g1)*g2+96*(13+118*g1*(-1+2*g1))*g2**2+ &
              9*sqrt5*(-1+2*g2))/(13440*(g1-g3)*(g2-g3)*(3-8*g2+ &
              8*g1*(-1+2*g2)-8*g3+16*(g1+g2)*g3))

        h1 = 1-3*g1
        h2 = 1-3*g2
        h3 = 1-3*g3

        x(1,1:3) = [ g1, g1, g1 ]
        x(2,1:3) = [ h1, g1, g1 ]
        x(3,1:3) = [ g1, h1, g1 ]
        x(4,1:3) = [ g1, g1, h1 ]
        x(5,1:3) = [ g2, g2, g2 ]
        x(6,1:3) = [ h2, g2, g2 ]
        x(7,1:3) = [ g2, h2, g2 ]
        x(8,1:3) = [ g2, g2, h2 ]
        x(9,1:3)  = [ g3, g3, g3 ]
        x(10,1:3) = [ h3, g3, g3 ]
        x(11,1:3) = [ g3, h3, g3 ]
        x(12,1:3) = [ g3, g3, h3 ]
        x(13,1:3)  = [ h4, p4, g4 ]
        x(14,1:3)  = [ p4, h4, g4 ]
        x(15,1:3)  = [ g4, h4, p4 ]
        x(16,1:3)  = [ g4, p4, h4 ]
        x(17,1:3)  = [ g4, g4, h4 ]
        x(18,1:3)  = [ g4, g4, p4 ]
        x(19,1:3)  = [ h4, g4, p4 ]
        x(20,1:3)  = [ p4, g4, h4 ]
        x(21,1:3)  = [ g4, h4, g4 ]
        x(22,1:3)  = [ g4, p4, g4 ]
        x(23,1:3)  = [ h4, g4, g4 ]
        x(24,1:3)  = [ p4, g4, g4 ]

        w(1:4) = w1 / 6
        w(5:8) = w2 / 6
        w(9:12) = w3 / 6
        w(13:24) = 27._dp/3360

      case default
        write(*,'(/2a,i0/)') &
          'Error in Gauss_Legendre_tetrahedron: ', &
          'invalid number of points, n = ', n
        stop

    end select

  end subroutine Gauss_Legendre_tetrahedron


! Gauss Legendre in 3D defined on the prism reference region

  subroutine Gauss_Legendre_prism ( n, n2, x, w )

!   the number of integration points on the triangle base
    integer, intent(in) :: n

!   the number of integration points on the height
    integer, intent(in) :: n2

!   the reference coordinates and weights of the integration points
!   x(ninti,ndim), w(ninti)
!   where ninti = n*n2
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w


    integer :: i, j, ip
    real(dp) :: x1(n,2), w1(n)
    real(dp), dimension(n2) :: x2, w2


!   test

    if ( size(x,1) < n*n2 .or. size(w) < n*n2 ) then
      write(*,'(/a/)') &
        'Error in Gauss_Legendre_prism: size of x or w too small'
      stop
    end if

    call Gauss_Legendre_triangle ( n, x1, w1 )
    call Gauss_Legendre_line ( n2, x2, w2 )

    do i = 1, n
      do j = 1, n2
        ip = i + n * ( j - 1 )
        x(ip,1) = x1(i,1)
        x(ip,2) = x1(i,2)
        x(ip,3) = x2(j)
        w(ip)   = w1(i) * w2(j)
      end do
    end do

  end subroutine Gauss_Legendre_prism


! Set number of Gauss integration points for an element

  subroutine set_ninti_standard ( gauss, ninti )

    type(gauss_t), intent(in) :: gauss

!   number of integration points
    integer, intent(out) :: ninti


!   set integration

    if ( gauss%globalshape == 'line' .or. &
         gauss%globalshape == 'triangle' .or. &
         gauss%globalshape == 'tetrahedron' ) then

      ninti = gauss%intrule

    else if ( gauss%globalshape == 'quadrilateral' ) then

      ninti = gauss%intrule ** 2

    else if ( gauss%globalshape == 'hexahedron' ) then

      ninti = gauss%intrule ** 3

    else if ( gauss%globalshape == 'prism' ) then

      ninti = gauss%intrule * gauss%intrule2

    else

      write(*,'(/a/2a/)') &
        'Error in set_ninti_standard: ', &
        '  standard Gauss integration not available for: ', gauss%globalshape
      stop

    end if

  end subroutine set_ninti_standard

end module gauss_standard_m
