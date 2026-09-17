
! Copyright (C) 2015-2015 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines standard shape functions based on the global coordinate system

module shapefunc_standard_global_m

  use kind_defs_m

  implicit none


contains


! P1

  subroutine shape_global_P1 ( x, phi, dphidx )

!   Global coordinates where phi and dphidx must be computed:
!   x(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: x

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function with respect to the global coordinates:
!   dphidx(i,j,k), with i the point in space j the unknown and k the
!   direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphidx


!   The unknowns of the element are positioned at x=0 and consist of the
!   the value in x=0 and the derivatives with respect to the global coordinate
!   system.

    integer :: j

    phi(:,1) = 1
    phi(:,2:) = x

    if ( present(dphidx) ) then

      dphidx = 0

      do j = 1, size(x,2)
        dphidx(:,j+1,j) = 1
      end do

    end if

  end subroutine shape_global_P1

end module shapefunc_standard_global_m
