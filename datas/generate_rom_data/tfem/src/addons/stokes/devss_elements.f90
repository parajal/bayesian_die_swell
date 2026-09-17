
! Copyright (C) 2005-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the DEVSS and DEVSS-G formulations
!
! DEVSS:
!
!   ... - 2 * alpha * div ( D(u) - E ) = ....
!
!                           D(u) - E   = 0
!
! Note: the last equation is multiplied by 2*alpha in the implementation to
! get a symmetric system matrix.
!
! DEVSS-G:
!
!   ... - alpha * div ( nabla u - G^T ) = ....
!
!                       nabla u - G^T   = 0
!
! Note: the last equation is multiplied by alpha in the implementation to
! get a symmetric system matrix.
!

module devss_elements_m

  use tfem_elem_m
  use devss_elements_generic_m
  use devss_elements_2D_m
  use devss_elements_3D_m

  implicit none

contains


! Internal element routine for the DEVSS Equation

  subroutine devss_elem ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    select case ( mesh%ndim )
      case(2) ! 2D and axisymmetric
        call devss_elem_2D ( mesh, problem, elgrp, elem, matrix, vector, &
          first, last, coefficients, oldvectors, elemmat, elemvec )
      case(3) ! 3D
        call devss_elem_3D ( mesh, problem, elgrp, elem, matrix, vector, &
          first, last, coefficients, oldvectors, elemmat, elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in devss_elem:', &
        ' devss element not available for mesh%ndim = ', mesh%ndim
        stop
    end select

  end subroutine devss_elem


! Internal element routine for the DEVSS-G Equation

  subroutine devssg_elem ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    select case ( mesh%ndim )
      case(2) ! 2D and axisymmetric
        if ( coefficients%i(67) == 1 ) then
          call devssg_elem_2D_vel3D ( mesh, problem, elgrp, elem, matrix, &
            vector, first, last, coefficients, oldvectors, elemmat, elemvec )
        else
          call devssg_elem_2D ( mesh, problem, elgrp, elem, matrix, vector, &
            first, last, coefficients, oldvectors, elemmat, elemvec )
        end if
      case(3) ! 3D
        call devssg_elem_3D ( mesh, problem, elgrp, elem, matrix, vector, &
          first, last, coefficients, oldvectors, elemmat, elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in devssg_elem:', &
        ' devss element not available for mesh%ndim = ', mesh%ndim
        stop
    end select

  end subroutine devssg_elem

end module devss_elements_m

