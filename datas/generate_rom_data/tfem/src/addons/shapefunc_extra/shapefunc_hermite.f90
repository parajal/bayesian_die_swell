
! Copyright (C) 2012-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines shape functions and associated routines for Hermite shapefunctions,
! for creating C^1 continuous trial functions.

module shapefunc_hermite_m

  use kind_defs_m
  use mesh_m, only: mesh_t

  implicit none


contains


! Line P3

  subroutine shape_line_hermite_P3 ( xr, phi, dphi, d2phi, d3phi )

!   Reference coordinates where phi, dphi and d2phi must be computed:
!   xr(i) with i the point in space
    real(dp), intent(in), dimension(:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j), with i the point in space
!   j the unknown
    real(dp), intent(out), dimension(:,:), optional :: dphi

!   second-order derivative of the shape function: d2phi(i,j), with i the
!   point in space j the unknown
    real(dp), intent(out), dimension(:,:), optional :: d2phi

!   third-order derivative of the shape function: d3phi(i,j), with i the
!   point in space j the unknown
    real(dp), intent(out), dimension(:,:), optional :: d3phi


!   The unknowns of the element are positioned at the following nodal points,
!
!     1 --------- 2
!
!   In both nodes there are two degrees of freedom: the field value u and the
!   derivative of the field with respect to the reference coordinate u'. The
!   sequence of the degrees of freedom on element level is:
!    [ u_1, u'_1, u_2, u'_2 ]
!
!   The reference coordinates xr are in the region [-1,1]

    phi(:,1) =  (1-xr)**2*(2+xr)/4
    phi(:,2) =  (1-xr)**2*(1+xr)/4
    phi(:,3) =  (1+xr)**2*(2-xr)/4
    phi(:,4) = -(1+xr)**2*(1-xr)/4

    if ( present(dphi) ) then
      dphi(:,1) = -3*(1-xr**2)/4
      dphi(:,2) = -(1-xr)*(1+3*xr)/4
      dphi(:,3) =  3*(1-xr**2)/4
      dphi(:,4) =  (1+xr)*(-1+3*xr)/4
    end if

    if ( present(d2phi) ) then
      d2phi(:,1) =  3*xr/2
      d2phi(:,2) = (-1+3*xr)/2
      d2phi(:,3) = -3*xr/2
      d2phi(:,4) =  (1+3*xr)/2
    end if

    if ( present(d3phi) ) then
      d3phi(:,1) =  1.5_dp
      d3phi(:,2) =  1.5_dp
      d3phi(:,3) = -1.5_dp
      d3phi(:,4) =  1.5_dp
    end if

  end subroutine shape_line_hermite_P3


! Line P4

  subroutine shape_line_hermite_P4 ( xr, phi, dphi, d2phi, d3phi )

!   Reference coordinates where phi, dphi and d2phi must be computed:
!   xr(i) with i the point in space
    real(dp), intent(in), dimension(:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j), with i the point in space
!   j the unknown
    real(dp), intent(out), dimension(:,:), optional :: dphi

!   second-order derivative of the shape function: d2phi(i,j), with i the
!   point in space j the unknown
    real(dp), intent(out), dimension(:,:), optional :: d2phi

!   third-order derivative of the shape function: d3phi(i,j), with i the
!   point in space j the unknown
    real(dp), intent(out), dimension(:,:), optional :: d3phi


!   The unknowns of the element are positioned at the following nodal points,
!
!     1 --------- 2 ---------- 3
!
!   In both end nodes there are two degrees of freedom: the field value u and
!   the derivative of the field with respect to the reference coordinate u'.
!   In the middle node there is an additional "bubble" degree of freedom.
!   The sequence of the degrees of freedom on element level is:
!    [ u_1, u'_1, b2, u_3, u'_3 ]
!
!   The reference coordinates xr are in the region [-1,1]

    real(dp) :: fac

    fac = 3._dp/16

    phi(:,1) =  (1-xr)**2*(2+xr)/4
    phi(:,2) =  (1-xr)**2*(1+xr)/4
    phi(:,3) =  fac*(1-xr**2)**2
    phi(:,4) =  (1+xr)**2*(2-xr)/4
    phi(:,5) = -(1+xr)**2*(1-xr)/4

    if ( present(dphi) ) then
      dphi(:,1) = -3*(1-xr**2)/4
      dphi(:,2) = -(1-xr)*(1+3*xr)/4
      dphi(:,3) = -4*fac*(1-xr**2)*xr
      dphi(:,4) =  3*(1-xr**2)/4
      dphi(:,5) =  (1+xr)*(-1+3*xr)/4
    end if

    if ( present(d2phi) ) then
      d2phi(:,1) =  3*xr/2
      d2phi(:,2) = (-1+3*xr)/2
      d2phi(:,3) = 4*fac*(3*xr**2-1)
      d2phi(:,4) = -3*xr/2
      d2phi(:,5) =  (1+3*xr)/2
    end if

    if ( present(d3phi) ) then
      d3phi(:,1) =  1.5_dp
      d3phi(:,2) =  1.5_dp
      d3phi(:,3) = 24*fac*xr
      d3phi(:,4) = -1.5_dp
      d3phi(:,5) =  1.5_dp
    end if

  end subroutine shape_line_hermite_P4

end module shapefunc_hermite_m
