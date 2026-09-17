
! Copyright (C) 2023-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


!
! All viscoelastic models in one go.
!

module viscoelastic_models_m

  use viscoelastic_models_2D_m
  use viscoelastic_models_2D_log_m
  use viscoelastic_models_2D_b_m
  use viscoelastic_models_3D_m
  use viscoelastic_models_3D_log_m
  use viscoelastic_models_3D_b_m

  implicit none

contains


! square root of conformation from the b tensor (2D, axisymmetric or 3D
! viscoelastic models). Optionally compute R from the polar decomposition
! b = V.R, with V=sqrt(b.b^T)=sqrt(c).

  subroutine sqrtc_b ( b, sqrtc, RT, bvariant )

!   INPUT
!     contravariant deformation tensor b in c = b*b^T
!     b(ip,comp) is component comp in point ip
!       ncompb=4, 2D: bxx, bxy, byx, byy
!       ncompb=5, 2D or axisymmetric:
!                bxx, bxy, byx, byy, bzz  (2D)
!                bzz, bzr, brz, brr, btt  (axisymmetric)
!       ncompb=9, 3D: bxx, bxy, bxz, byx, byy, byz, bzx, bzy, bzz
!     where ncompb=size(b,2)
!   OUTPUT
!     if present(sqrtc)
!       b is not modified
!     otherwise it contains the V-tensor in unsymmetric storage:
!       V(ip,comp) is component comp in point ip
!       ncompb=4, 2D: Vxx, Vxy, Vyx=Vxy, Vyy
!       ncompb=5, 2D or axisymmetric:
!                Vxx, Vxy, Vyx=Vxy, Vyy, Vzz  (2D)
!                Vzz, Vzr, Vrz=Vzr, Vrr, Vtt  (axisymmetric)
!       ncompb=9, 3D: Vxx, Vxy, Vxz, Vyx, Vyy, Vyz, Vzx, Vzy, Vzz
!     where ncompb=size(b,2)
    real(dp), dimension(:,:), intent(inout) :: b

!   If present
!     the symmetric V tensor = sqrt(b.b^T) = sqrt(c) in symmetric storage.
!     sqrtc(ip,comp) is component comp in point ip
!       ncompb=4, 2D: Vxx, Vxy, Vyy
!       ncompb=5, 2D or axisymmetric:
!                 Vxx, Vxy, Vyy, Vzz  (2D)
!                 Vzz, Vzr, Vrr, Vtt  (axisymmetric)
!       ncompb=9, 3D: Vxx, Vxy, Vxz, Vyy, Vyz, Vzz
!     where ncompb=size(b,2),
!   otherwise the V tensor is output via the argument b.
    real(dp), dimension(:,:), intent(out), optional :: sqrtc

!   optional (transpose of) the rotation tensor from polar decomposition
!   RT(ip,i,j), where i,j=1,nd with nd=2 for 2D, axisymmetric and nd=3 for 3D.
    real(dp), dimension(:,:,:), intent(out), optional :: RT

!   Optionally the b variant can be given.
!   Only needed for bvariant=4 (Cholesky-log)
!   default=1
    integer, intent(in), optional :: bvariant

    integer :: ncompb

    ncompb = size(b,2)

    select case ( ncompb )
    case(4,5) ! 2D or axisymmetric
      call sqrtc_2D_b ( b, sqrtc, RT, bvariant )
    case(9) ! 3D
      call sqrtc_3D_b ( b, sqrtc, RT, bvariant )
    case default
      write(*,'(a,i0)') ' Error sqrtc_b: Invalid ncompb = ', ncompb
      stop
    end select

  end subroutine sqrtc_b


! rotate contravariant deformation for 2D, axisymmetric and 3D
! viscoelastic models.

  subroutine rotate_b ( b, RT )

!   INPUT
!     contravariant deformation tensor b in c = b*b^T
!     b(ip,comp) is component comp in point ip
!       ncompb=4, 2D: bxx, bxy, byx, byy
!       ncompb=5, 2D or axisymmetric:
!                bxx, bxy, byx, byy, bzz  (2D)
!                bzz, bzr, brz, brr, btt  (axisymmetric)
!       ncompb=9, 3D: bxx, bxy, bxz, byx, byy, byz, bzx, bzy, bzz
!     where ncompb=size(b,2)
!   OUTPUT
!     components of b.R^T
    real(dp), dimension(:,:), intent(inout) :: b

!   (transpose of) the rotation tensor R from polar decomposition
!   RT(ip,i,j), where i,j=1,nd with nd=2 for 2D, axisymmetric and nd=3 for 3D.
    real(dp), dimension(:,:,:), intent(in) :: RT

    integer :: ncompb

    ncompb = size(b,2)

    select case ( ncompb )
    case(4,5) ! 2D or axisymmetric
      call rotate_2D_b ( b, RT )
    case(9) ! 3D
      call rotate_3D_b ( b, RT )
    case default
      write(*,'(a,i0)') ' Error rotate_b: Invalid ncompb = ', ncompb
      stop
    end select

  end subroutine rotate_b


! skew norm contravariant deformation tensor for 2D, axisymmetric and 3D
! viscoelastic models. Not applicable to Cholesly-log.

  function skew_norm_b ( b )

!   INPUT
!     contravariant deformation tensor b in c = b*b^T
!     b(ip,comp) is component comp in point ip
!       ncompb=4, 2D: bxx, bxy, byx, byy
!       ncompb=5, 2D or axisymmetric:
!                bxx, bxy, byx, byy, bzz  (2D)
!                bzz, bzr, brz, brr, btt  (axisymmetric)
!       ncompb=9, 3D: bxx, bxy, bxz, byx, byy, byz, bzx, bzy, bzz
!     where ncompb=size(b,2)
    real(dp), dimension(:,:), intent(in) :: b

    real(dp), dimension(size(b,1)) :: skew_norm_b

    integer :: ncompb

    ncompb = size(b,2)

    select case ( ncompb )
    case(4,5) ! 2D or axisymmetric
      skew_norm_b = skew_norm_2D_b ( b )
    case(9) ! 3D
      skew_norm_b = skew_norm_3D_b ( b )
    case default
      write(*,'(a,i0)') ' Error rotate_b: Invalid ncompb = ', ncompb
      stop
    end select

  end function skew_norm_b

end module viscoelastic_models_m
