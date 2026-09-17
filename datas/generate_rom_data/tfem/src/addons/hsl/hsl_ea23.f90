
! Copyright (C) 2022-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

!
! eigenvalues and vectors for 3D tensors using ea23 from the HSL library
!

module hsl_ea23_m

  use kind_defs_m

  implicit none

  interface eig3x3
    module procedure eig3x3_hsl
  end interface eig3x3

contains

! Compute eigenvalues and eigenvectors for a 3x3 symmetric matrix (HSL)

  subroutine eig3x3_hsl ( c, lambda, eigv )

    real(dp), intent(in) :: c(6)
    real(dp), intent(out) :: lambda(3), eigv(3,3)

!   Compute eigenvalues and eigenvectors for a 3x3 symmetric matrix

    integer :: iw(3), irot
    real(dp) :: a(3,3)


    a(1,:) = c(1:3)
    a(2,2:3) = c(4:5)
    a(3,3) = c(6)

!   --- HSL routine (jacobi)

    irot = 20

    call ea23ad ( 3, a, 3, lambda, eigv, .true., irot, iw, 6 )

  end subroutine eig3x3_hsl

end module hsl_ea23_m

