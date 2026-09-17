
! This module contains subroutines for:
!  1) the Gauss-Legendre-Lobatto (GLL) points and the according weight function
!     in order to compute spectral (line) elements
!
! Based on code by A.C. Verkaik (26-11-2009)
! Based on/copied from the book of Canuto et al. 1988, Appendix C (with alpha=0
! and beta=0, because Legendre method is used) and for some theoretical
! explanation see also Appendix B2 of the book written by Karniadakis & Sherwin
! 2005 (2nd edition)

module gauss_spectral_m

  use kind_defs_m
  use gauss_defs_m

  implicit none

contains


! compute the weights for the Gauss-Legendre-Lobatto integration

  subroutine GLL_weights ( x, p, w )

!   array filled with GLL points/integration points (ranging from -1 to 1)
    real(dp), intent(in), dimension(:) :: x

!   array filled with weights of the corresponding GLLpoints/integration points
!   (added all together should give the value 2)
    real(dp), intent(out), dimension(:) :: w

!   order of the GLL polynomial and p+1 is the number of GLLpoints
    integer, intent(in) :: p

!   local variables
    integer :: j
    real(dp) :: xj, poly
    real(dp) :: dum1, dum2, dum3, dum4, dum5 ! dummy variables

    do j = 1, p+1
      xj = x(j)
      call GLL_polynomials(p, poly, dum1, dum2, dum3, dum4, dum5, xj)
!     calculation of the weights according to formula (see Canuto for example)
      w(j)= (2.0_dp/(p*(p+1))) * (1.0_dp/(poly**2._dp))
    end do

  end subroutine GLL_weights


! compute the Gauss-Legendre-Lobatto integration points

  subroutine GLL_points ( x, p )

!   x will be filled with the GLL points/coordinates (ranging from -1 to 1)
    real(dp), intent(out), dimension(:) ::  x

!   the polynomial order
    integer, intent(in) ::  p

!   declare other variables
    integer :: m, nh
    real(dp) :: det, rp, rn, a, b, rm
    real(dp) :: pnp1p, pdnp1p, pnp, pdnp, pnm1p, pdnm1
    real(dp) :: pnp1m, pdnp1m, pnm, pdnm, pnm1m
    real(dp) :: dth, cd, sd, cs, ss, y
    real(dp) :: pnp1, pdnp1, pn, pdn, pnm1
    real(dp) :: poly, pder, recsum, eps, dely, cssave

    integer :: j, k, i

!   the computation part

    m=p+1 ! number of integration points in one element
    eps = 1.e-14 ! the accuracy of the estimation of the GLL points

    call GLL_polynomials( m, pnp1p, pdnp1p, pnp, pdnp, pnm1p, pdnm1, 1.0_dp )
    call GLL_polynomials( m, pnp1m, pdnp1m, pnm, pdnm, pnm1m, pdnm1, -1.0_dp )

    det = pnp*pnm1m - pnm*pnm1p
    rp = -pnp1p
    rn = -pnp1m
    a = (rp*pnm1m - rn*pnm1p)/det
    b = (rn*pnp -rp*pnm)/det

    rm = m
    x(m) = 1.0_dp
    nh = m/2

!   now some smart first estimation method for a 'zero' point is used is used
!   'a recursion relation is set-up to obtain the initial guess for the roots
!   of the Gauss-Legendre-Lobatto polynomial'

    dth = 4.0_dp*atan(1.0_dp)/(2.0_dp*rm - 1.0_dp)
    cd = cos(2.0_dp*dth)
    sd = sin(2.0_dp*dth)
    cs = cos(dth) ! is the Chebyschev point, which is used as a first estimater
    ss = sin(dth)

!   now the roots of the polynomial are computed by polynomial deflation
!   first the first half (and symmetry is used to compute the second half)

    do j = 2, nh

      y = cs

      do k =1, 13  ! maximum number of iterations is 13, but if a certain
                   ! accuracy is reached the calculation will stop

        call GLL_polynomials ( m, pnp1, pdnp1, pn, pdn, pnm1, pdnm1, y)

        poly = pnp1 + a*pn + b*pnm1 ! this has to be equal to L'_{k+1}
        pder = pdnp1 + a*pdn + b*pdnm1 ! this has to be equal to L''_{k+1}

        recsum = 0.0_dp
        do i = 1, j-1
          recsum = recsum + 1.0_dp/(y - x(m+1-i))
        end do

!       Newton-Raphson method is used to get the estimation of one of
!       the roots points
        dely = -poly/(pder - recsum*poly)
        y = y + dely

!       check whether the desired accuracy is achieved
        if ( abs(dely) < eps ) then
!          write(*,*) 'k:', k
!          write(*,*) 'dely:', dely
          exit
        end if

      end do

!      write(*,*) 'accuracy of GLLpoint estimation:', dely

      x(m+1-j) = y   ! save the found root point
!     some trick to get the new estimation for the next root point
      cssave = cs*cd - ss*sd
      ss = cs*sd + ss*cd
      cs = cssave

    end do

    x(1)=-1.0_dp

!   symmetry is used to fill the other of Gauss-Legendre-Lobatto polynomial
!   roots (the negative values)
    do i = 2, nh
      x(i) = -x(m+1-i)
    end do

!   check whether there is a middle point and fill it with zero if
!   the number of points is odd
    if ( mod(m,2)==1) then
      x(nh+1)=0.0_dp
    end if

  end subroutine GLL_points


!  compute the jacobi polynomial and derivatives, with alpha=0 and beta=0
!  (the Legendre method) of order p, at point xi

  subroutine GLL_polynomials ( p, poly, pder, polymin1, pdermin1, polymin2, &
    pdermin2, xi)

   integer, intent(in) :: p
   real(dp), intent(in) :: xi
   real(dp), intent(out) :: poly, pder, polymin1, pdermin1, polymin2, pdermin2

!  other variables
   real(dp) :: x, polyls, pderls, polyk, pderk, psave, pdsave
   integer :: k

   x=xi
   poly = 1.0_dp
   pder = 0.0_dp

   if (p==0) return

   polyls = poly
   pderls = pder
   poly = x
   pder = 1.0_dp

   if (p==1) return

   do k=2,p
!    recursion relation for the legendre polynomial of order k
     polyk=((2*k-1)*x*poly-(k-1)*polyls)/k
!    recursion relation for the derivative of the legendre polynomial of order k
     pderk=((2*k-1)*x*pder-(k-1)*pderls+(2*k-1)*poly)/k
     psave=polyls ! dummy to save the former solution
     pdsave=pderls
     polyls=poly  ! to save the k-1 solution
     poly=polyk
     pderls=pder  ! to save the derivative of the k-1 solution
     pder=pderk
   end do

   polymin1=polyls
   pdermin1=pderls
   polymin2=psave
   pdermin2=pdsave

  end subroutine GLL_polynomials


! Gauss Legendre Lobatto in 1D defined on the interval [-1,1]

  subroutine Gauss_Legendre_Lobatto_line ( n, x, w )

!   number of integration points
    integer, intent(in) :: n

!   the reference coordinates and weights of the integration points
!   x(i), w(i), i=1,n
    real(dp), intent(out), dimension(:) :: x
    real(dp), intent(out), dimension(:), optional :: w

    integer :: p

    p = n - 1 ! polynomial order

!   calculation of the GLL points
    call GLL_points ( x, p )

    if ( present(w) ) then

!     calculation of the weights of the integration points
      call GLL_weights ( x, p, w )

    end if

  end subroutine Gauss_Legendre_Lobatto_line


! Gauss Legendre Lobatto in 2D defined on the region [-1,1]x[-1,1] (quad)

  subroutine Gauss_Legendre_Lobatto_quad ( n, x, w )

!   number of integration points = order of GLLpolynomial + 1
    integer, intent(in) :: n

!   the reference coordinates and weights of the integration points
!   x(i,ndim), w(i), i=1,ninti
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:), optional :: w

    integer :: p, i, j, ip
    real(dp), dimension(n) :: x1, w1

    p = n - 1 ! polynomial order

!   calculation of the GLLpoints
    call GLL_points ( x1, p )

    if ( present(w) ) then

!     Now calculation of the 2D points and weights
!     until now only equal order in x- and y-direction
!     Could be improved by taking symmetries into account to save CPU time.

!     calculation of the weights of the integration points (GLLpoints)
      call GLL_weights ( x1, p, w1 )

      do i = 1, n
        do j = 1, n
          ip = i + n * (j-1)
          x(ip,1) = x1(i)
          x(ip,2) = x1(j)
          w(ip)   = w1(i)*w1(j)
        end do
      end do

    else

!     Now calculation of only the 2D points
!     until now only equal order in x- and y-direction

      do i = 1, n
        do j = 1, n
          ip = i + n * (j-1)
          x(ip,1) = x1(i)
          x(ip,2) = x1(j)
        end do
      end do

    end if

  end subroutine Gauss_Legendre_Lobatto_quad


! Gauss Legendre Lobatto in 3D defined on the region [-1,1]x[-1,1]x[-1,1] (hex)

  subroutine Gauss_Legendre_Lobatto_hexahedron ( n, x, w )

!   number of integration points = order of GLLpolynomial + 1
    integer, intent(in) :: n

!   the reference coordinates and weights of the integration points
!   x(i,ndim), w(i), i=1,ninti
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w

    integer :: p, i, j, k, ip
    real(dp), dimension(n) :: x1, w1

    p = n - 1 ! polynomial order

!   calculation of the GLLpoints
    call GLL_points ( x1, p )

!   calculation of the weights of the integration points (GLLpoints)
    call GLL_weights ( x1, p, w1 )

!   Now calculation of the 3D points and weights
!   until now only equal order in x-, y- and z-direction
!   Could be improved by taking symmetries into account to save CPU time.

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

  end subroutine Gauss_Legendre_Lobatto_hexahedron


! Set number of Gauss integration points for an element

  subroutine set_ninti_spectral ( gauss, ninti )

    type(gauss_t), intent(in) :: gauss

!   number of integration points
    integer, intent(out) :: ninti


!   set integration

    if ( gauss%globalshape == 'line' ) then

      ninti = gauss%intrule

    else if ( gauss%globalshape == 'quadrilateral' ) then

      ninti = gauss%intrule ** 2

    else if ( gauss%globalshape == 'hexahedron' ) then

      ninti = gauss%intrule ** 3

    else

      write (*, '(/a/3a/)' ) &
        ' Error set_ninti_spectral: ', &
        ' gauss%globalshape == ', gauss%globalshape, ' not available '
      stop

    end if

  end subroutine set_ninti_spectral

end module gauss_spectral_m
