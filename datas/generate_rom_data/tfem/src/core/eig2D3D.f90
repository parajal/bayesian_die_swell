
! Copyright (C) 2018-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

!
! eigenvalues and vectors for 2D and 3D tensors
!

module eig2D3D_m

  use kind_defs_m

  implicit none

  interface eig3x3
    module procedure eig3x3_lapack
  end interface eig3x3

contains


! Compute eigenvalues and eigenvectors for a 2x2 symmetric matrix analytically

  subroutine eig2x2 ( c, lambda, eigv )

    real(dp), intent(in) :: c(3)
    real(dp), intent(out) :: lambda(2), eigv(2,2)

!   Compute eigenvalues and eigenvectors for a 2x2 symmetric matrix analytically
!   See note 287. Derived from elx270 of dynaflow.

    real(dp) :: c11, c12, c22, D, c11mc22, sqrtD, length, trc, sqrt2


    c11 = c(1)
    c12 = c(2)
    c22 = c(3)

    c11mc22 = c11 - c22
    D = c11mc22 ** 2 + 4 * c12 ** 2
    sqrtD = sqrt(D)
    trc = c11+c22
    lambda(1) = (trc+sqrtD)/2
    lambda(2) = lambda(1)-sqrtD
    sqrt2 = sqrt(2._dp)/2

    if ( c11 > c22 ) then
      eigv(1,1) = c11mc22 + sqrtD
      eigv(2,1) = 2*c12
      length = sqrt(2*sqrtD*(sqrtD+c11mc22))
      eigv(1,1) = eigv(1,1) / length
      eigv(2,1) = eigv(2,1) / length
      eigv(1,2) = -eigv(2,1)
      eigv(2,2) = eigv(1,1)
    else if ( c11 < c22 ) then
      eigv(1,1) = 2*c12
      eigv(2,1) = -c11mc22 + sqrtD
      length = sqrt(2*sqrtD*(sqrtD-c11mc22))
      eigv(1,1) = eigv(1,1) / length
      eigv(2,1) = eigv(2,1) / length
      eigv(1,2) = -eigv(2,1)
      eigv(2,2) = eigv(1,1)
    else if ( c12 /= 0 ) then
      eigv(1,1) = sign(sqrt2,c12)
      eigv(2,1) =  sqrt2
      eigv(1,2) = -sqrt2
      eigv(2,2) = eigv(1,1)
    else
      lambda(1) = trc / 2
      lambda(2) = lambda(1)
      eigv(1,1) = 1
      eigv(2,1) = 0
      eigv(1,2) = 0
      eigv(2,2) = 1
    end if

  end subroutine eig2x2


! Compute a 2x2 symmetric matrix from eigenvalues and eigenvectors

  subroutine inveig2x2 ( c, lambda, eigv )

    real(dp), intent(out) :: c(3)
    real(dp), intent(in) :: lambda(2), eigv(2,2)

!
!   Compute symmetric tensor from eigenvalues and eigenvectors
!
!     c = V * lambda * V^T
!

    c(1) = lambda(1) * eigv(1,1)**2 + lambda(2) * eigv(1,2)**2
    c(2) = lambda(1) * eigv(1,1)*eigv(2,1) + lambda(2) * eigv(1,2)*eigv(2,2)
    c(3) = lambda(1) * eigv(2,1)**2 + lambda(2) * eigv(2,2)**2

  end subroutine inveig2x2


! Compute eigenvalues and eigenvectors for a 3x3 symmetric matrix (Lapack)

  subroutine eig3x3_lapack ( c, lambda, eigv )

    real(dp), intent(in) :: c(6)
    real(dp), intent(out) :: lambda(3), eigv(3,3)

!   Compute eigenvalues and eigenvectors for a 3x3 symmetric matrix

    integer, parameter :: lwork = 102
    integer :: info
    real(dp) :: work(lwork)

    eigv(1,:) = c(1:3)
    eigv(2,2:3) = c(4:5)
    eigv(3,3) = c(6)

!   Lapack routine

    call dsyev ( 'V', 'U', 3, eigv, 3, lambda, work, lwork, info )

    if ( info /= 0 ) stop ' Error in eig3x3_lapack:dsyev '

  end subroutine eig3x3_lapack


! Compute a 3x3 symmetric matrix from eigenvalues and eigenvectors

  subroutine inveig3x3 ( c, lambda, eigv )

    real(dp), intent(out) :: c(6)
    real(dp), intent(in) :: lambda(3), eigv(3,3)

!
!   Compute symmetric tensor from eigenvalues and eigenvectors
!
!     c = V * lambda * V^T
!

    c(1) = sum ( lambda * eigv(1,:) * eigv(1,:) )
    c(2) = sum ( lambda * eigv(1,:) * eigv(2,:) )
    c(3) = sum ( lambda * eigv(1,:) * eigv(3,:) )
    c(4) = sum ( lambda * eigv(2,:) * eigv(2,:) )
    c(5) = sum ( lambda * eigv(2,:) * eigv(3,:) )
    c(6) = sum ( lambda * eigv(3,:) * eigv(3,:) )

  end subroutine inveig3x3


end module eig2D3D_m

