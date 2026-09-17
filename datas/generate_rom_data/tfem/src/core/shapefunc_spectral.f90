
! This module contains subroutines for computing
!  1) the shapefunction (is always 1 at the GLL points (xi)) and
!     the derivative of the shapefunctions in the GLL points (xi)
!  2) the shapefunction and the derivative of the shapefunction in
!     any reference point (not complete yet)
!
! Based on routines by A.C. Verkaik.
! Some is based on/copied from the book of Canuto et al. 1988, Appendix C
! (with alpha=0 and beta=0, because Legendre method is used) and for some
! theoretical explanation see also Appendix B2 of the book written by
! Karniadakis & Sherwin 2005 (2nd edition)

module shapefunc_spectral_m

  use kind_defs_m
  use gauss_spectral_m
  use set_optional_m
  use shapefunc_highorder_m, only: shape_line_lagrange

  implicit none

! interface to support two interfaces for shape_line_GLL

  interface shape_line_GLL
    module procedure shape_line_GLL_1, shape_line_GLL_2
  end interface shape_line_GLL

contains


! compute the value of the shapefunction at the GLL points (dof's are at the
! GLL points and if necessary also the derivatives at xi, which is a
! Gauss-Legendre-Lobatto (integration) point

  subroutine shape_line_at_GLL ( xr, phi, dphi, p )

!   the integration points (Gauss-Legendre-Lobatto)
    real(dp), intent(in), dimension(:) :: xr

!   value of the different shapefunctions at point xi
    real(dp), intent(out), dimension(:,:) :: phi

!   first derivative of the shapefunctions in point xi, which has to be one
!   of the GLL points
    real(dp), intent(out), dimension(:,:), optional :: dphi

!   polynomial order of Gauss-Legendre-Lobatto polynomial (optional)
!   default=size(x)-1
    integer, intent(in), optional :: p

!   N is number of GLL points, xi is the GLL point for which the
!   shapefunction and derivatives will be calculated
    integer :: N, i, j, pl
    real(dp) :: dum1, dum2, dum3, dum4 ! dummy variables
    real(dp) :: poly, pder, poly2, pder2
    real(dp) :: xi, xj

    pl = set_optional ( variable=p, default=size(xr)-1 )

    N=pl+1   ! number of GLL/integration points

!   the shapefunction phi(j) value is equal to 1 if xi=xr(j) and zero if xi
!   is equal to one of the other GLLpoints.

    do j = 1, N
      do i = 1, N
        if (i==j) then
          phi(j,i)=1.0_dp
        else
          phi(j,i)=0.0_dp
        end if
      end do
    end do

    if ( present(dphi) ) then

!     calculation of the first derivative of the shapefunction, but this only
!     holds if xi is equal to one of the GLLpoints

      do i = 1, N

        xi = xr(i)
        call GLL_polynomials ( pl, poly, pder, dum1, dum2, dum3, dum4, xi )

        do j = 1, N
          if (j==i) then
            dphi(i,j) = 0._dp
          else
            xj = xr(j)
            call GLL_polynomials ( pl, poly2, pder2, dum1, dum2, dum3, dum4, &
              xj )
            dphi(i,j) = poly/((xi-xj) * poly2)
          end if
        end do

      end do

!     adjustment of derivative values at the end points
      dphi(1,1) = -pl * N / 4.0_dp
      dphi(N,N) = pl * N / 4.0_dp

    end if

  end subroutine shape_line_at_GLL


! compute the value of the shapefunction and derivative at any reference point
! (dof's are at the GLL points)
! (one dimension in xr)

  subroutine shape_line_GLL_1 ( xr, p, phi, dphi, xg )

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

!   If present: the reference coordinates of the Gauss Legendre Lobatto points
!   xg(i) with i the point in space
!   If not present the GLL points are computed, which is more expensive.
    real(dp), intent(in), dimension(:), optional :: xg

    real(dp), dimension(p+1) :: xgl ! local array for Gauss points


!   GLL points

    if ( present(xg) ) then
      xgl = xg
    else
      call GLL_points ( xgl, p )
    end if

!   compute phi and optionally dphi by Lagrange interpolation

    call shape_line_lagrange ( xr, xgl, phi, dphi )

  end subroutine shape_line_GLL_1


! compute the value of the shapefunction and derivative at any reference point
! (dof's are at the GLL points)
! (two dimensions in xr)

  subroutine shape_line_GLL_2 ( xr, p, phi, dphi, xg )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,1) with i the point in space
    real(dp), intent(in), dimension(:,:) :: xr

!   polynomial order
    integer, intent(in) :: p

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,1), with i the point in space
!   j the unknown
    real(dp), intent(out), dimension(:,:,:), optional :: dphi

!   If present: the reference coordinates of the Gauss Legendre Lobatto points
!   xg(i) with i the point in space
!   If not present the GLL points are computed, which is more expensive.
    real(dp), intent(in), dimension(:), optional :: xg

    real(dp), dimension(p+1) :: xgl ! local array for Gauss points


!   GLL points

    if ( present(xg) ) then
      xgl = xg
    else
      call GLL_points ( xgl, p )
    end if

!   compute phi and optionally dphi by Lagrange interpolation

    call shape_line_lagrange ( xr(:,1), xgl, phi, dphi(:,:,1) )

  end subroutine shape_line_GLL_2


! compute the value of the shapefunction at the GLL points (dof's are at the
! GLL points and if necessary also the derivatives at xi, which is a
! Gauss-Legendre-Lobatto (integration) point

  subroutine shape_quad_at_GLL ( xr, p, phi, dphi )

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   polynomial order of Gauss-Legendre-Lobatto polynomial
    integer, intent(in) :: p

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


    real(dp), dimension(:,:), allocatable :: phi1, dphi1

!   n is number of GLL points in one direction
    integer :: i, j, k, l, r, s, n, np


    n=p+1   ! number of GLL/integration points in one direction

    allocate ( phi1(n,n), dphi1(n,n) )

!   calculate phi and dphi for 1D
!   phi1=shapefunction for 1 D element with order p,

    call shape_line_at_GLL ( xr(1:n,1), phi1, dphi1, p )

!   it is assumed that the polynomial order of the GLL-polynomial is equal
!   in x- and y-direction

    np = size(xr(:,1))

    do i = 1, np
      do j = 1, np
        if (i==j) then
           phi(i,j) = 1._dp
        else
           phi(i,j) = 0._dp
        end if
      end do
    end do

    if ( present(dphi) ) then

      dphi = 0

!     Only fill the non-zeros

      do i = 1, n
        do j = 1, n
          r = (j-1)*n + i
          l = j
          s = (l-1)*n
          do k = 1, n
            s = s + 1
            dphi(s,r,1) = dphi1(k,i)
          end do
          k = i
          do l = 1, n
            s = (l-1)*n + k
            dphi(s,r,2) = dphi1(l,j)
          end do
        end do
      end do

    end if

    deallocate ( phi1, dphi1 )

  end subroutine shape_quad_at_GLL


! compute the value of the shapefunction and derivative at any reference point
! (dof's are at the GLL points)

  subroutine shape_quad_GLL ( xr, p, phi, dphi, xg )

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

!   If present: the reference coordinates of the Gauss Legendre Lobatto points
!   in 1D, xg(i) with i the point in space
!   If not present the GLL points are computed, which is more expensive.
    real(dp), intent(in), dimension(:), optional :: xg

    integer :: i, j
    integer, allocatable, dimension(:,:) :: s
    real(dp), allocatable, dimension(:,:) :: phi1, phi2, dphi1, dphi2


    allocate(phi1(size(xr,1),p+1))
    allocate(phi2,dphi1,dphi2,mold=phi1)
    s = reshape ( [(i,i=1,(p+1)**2)], [p+1,p+1] )


    if ( present(dphi) ) then

      call shape_line_GLL ( xr(:,1), p, phi1, dphi1, xg )
      call shape_line_GLL ( xr(:,2), p, phi2, dphi2, xg )

      do i = 1, p+1
        do j = 1, p+1
          phi(:,s(i,j))    =  phi1(:,i) *  phi2(:,j)
          dphi(:,s(i,j),1) = dphi1(:,i) *  phi2(:,j)
          dphi(:,s(i,j),2) =  phi1(:,i) * dphi2(:,j)
        end do
      end do

    else

      call shape_line_GLL ( xr(:,1), p, phi1, xg=xg )
      call shape_line_GLL ( xr(:,2), p, phi2, xg=xg )

      do i = 1, p+1
        do j = 1, p+1
          phi(:,s(i,j)) =  phi1(:,i) *  phi2(:,j)
        end do
      end do

    end if

  end subroutine shape_quad_GLL


! Compute the value of the shapefunction and derivative at GLL reference points
! for a polynomial Pp. The shapefunction is based on Legendre polynomials.
! A constant is subtracted from the higher-order even Legendre polynomials
! to achieve that in the center only phi(:,1) (=1) is non-zero.

  subroutine shape_quad_Pp_at_GLL ( xr, p, phi, dphi, n )

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

!   number of GLL points in one dimension (assuming equal order)
    integer, intent(in), optional :: n


    integer     :: i, j, k, l, m, s
    real(dp)    :: x
    real(dp), allocatable, dimension(:) :: L0
    real(dp), allocatable, dimension(:,:) :: lx


    if ( present(dphi) ) then
      write(*,'(/a/)') &
        'Error in shape_quad_Pp: dphi not yet implemented'
      stop
    end if


    allocate ( lx(0:p,n), L0(0:p) )

!   fill L0 (=Ln(0)) for even polynomial degrees

    L0(0) = 1
    do j = 2, p, 2
      L0(j) = -(j-1)*L0(j-2)/j
    end do


!   fill lx

    do i = 1, n

      x = xr(i,1)

!     compute Legrende polynomials Ln (n=0,...,p) up to order p

      lx(0,i) = 1._dp

      lx(1,i) = x

      do j = 1, p-1
        lx(j+1,i) = ((2*j+1)*x*lx(j,i) - j*lx(j-1,i))/(j+1)
      end do

!     subtract Ln(0) for even order to get lx=0 at the middle point
!     of an element so that the first degree of freedom is the value on the
!     center of the element. Note, that orthogonality is lost.

      do j = 2, p, 2
        lx(j,i) = lx(j,i) - L0(j)
      end do

    end do


!   fill shape functions

    do i = 1, n
      do j = 1, n

        s = i + (j-1)*n

!       P0

        phi(s,1) = 1._dp

!       P1,...,Pp

        m = 1
        do l = 1, p
          do k = 0, l
            phi(s,m+k+1) = lx(k,i) * lx(l-k,j)
!            phi(s,m+k+1) = lx(l-k,i) * lx(k,j)
          end do
          m = m + l + 1
        end do

      end do
    end do

    deallocate ( lx, L0 )

  end subroutine shape_quad_Pp_at_GLL

end module shapefunc_spectral_m
