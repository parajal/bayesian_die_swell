
! Copyright (C) 2005-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the constitutive equation
!

module viscoelastic_elements_m

  use stokes_elements_m
  use devss_elements_m
  use viscoelastic_elements_SUPG_m
  use viscoelastic_elements_implicit_m
  use viscoelastic_elements_gammap_m
  use viscoelastic_elements_DG_2D_m
  use viscoelastic_elements_embedded_m
  use viscoelastic_elements_integrate_m

  implicit none

  save

contains

! Element for the drag force due to the polymer stress using
! integrate_boundary_elements in the postprocessing_m module.
!
! oldvectors:
!  v1(1)%p(1:nmodes) = conformation tensors derived from subroutine
!                "derive_conformation_tensor_std" for _both_ the standard
!                and the log conformation. Thus in the log formulation
!                they contain the nodal values of logc.

  subroutine viscoelastic_drag ( mesh, problem, curve, elem, first, last, &
    coefficients, oldvectors, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    select case ( mesh%ndim )
      case(2) ! 2D and axisymmetric
        call viscoelastic_drag_curve ( mesh, problem, curve, elem, first, &
          last, coefficients, oldvectors, elemvec )
      case(3) ! 3D
        call viscoelastic_drag_surface ( mesh, problem, curve, elem, first, &
          last, coefficients, oldvectors, elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in viscoelastic_drag:', &
        ' element not available for mesh%ndim = ', mesh%ndim
        stop
    end select

  end subroutine viscoelastic_drag

end module viscoelastic_elements_m

