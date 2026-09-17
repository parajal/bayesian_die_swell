
! Copyright (C) 2005-2006 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


module bcf_functions_m

  use kind_defs_m

  implicit none

  real(dp), allocatable, dimension(:,:) :: brown_bcf
  save

contains


! inflow boundary conditions for the Q-vector in DG

  function qinflow ( nfield, ncompq, x, un )

!   number of components in the function
    integer, intent(in) :: nfield, ncompq

!   coordinates of the point
    real(dp), intent(in), dimension(:) :: x

!   normal velocity (should be negative)

    real(dp), intent(in) :: un

!   the inflow value of the q tensor at the point x
    real(dp), dimension(nfield,ncompq) :: qinflow

    real(dp), parameter :: eps = 1e-14_dp

    if ( un <= -eps ) then
      write(*,*) 'Function qinflow has not been defined'
      write(*,*) 'Coordinates:', x, 'Normal velocity', un
      stop
    end if

    qinflow = 0

  end function qinflow


! Brownian vector

  function brownian_vector ( nfield, ncompq )

!   number of components in the function
    integer, intent(in) :: nfield, ncompq

!   the Brownian vector
    real(dp), dimension(nfield,ncompq) :: brownian_vector

    brownian_vector = brown_bcf

  end function brownian_vector

end module bcf_functions_m
