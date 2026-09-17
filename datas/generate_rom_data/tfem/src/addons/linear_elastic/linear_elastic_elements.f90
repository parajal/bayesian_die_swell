
! Copyright (C) 2020-2020 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the linear elastic equation (Navier)
!
!    - div( lambda div u I + mu (nabla u+nabla u^T) ) = f
!
! in displacement formulation or
!
!    - div( mu (nabla u+nabla u^T) ) + nabla p = f
!      p/lambda + div u = 0
!
! in mixed displacement-pressure formulation.
! For plane stress (2D, Cartesian) lambda needs to be replaced by lambda_bar.

module linear_elastic_elements_m

  use tfem_elem_m
  use linear_elastic_elements_generic_m
  use linear_elastic_elements_2D_curve_m
  use linear_elastic_elements_3D_surface_m
  use linear_elastic_elements_stabilized_m
  use linear_elastic_elements_integrate_m

  implicit none

contains


! Internal element routine for the Stokes Equation

  subroutine linear_elastic_elem ( mesh, problem, elgrp, elem, matrix, vector, &
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
          call linear_elastic_elem_2D_disp3D ( mesh, problem, elgrp, elem, &
            matrix, vector, first, last, coefficients, oldvectors, elemmat, &
            elemvec )
        else
          call linear_elastic_elem_2D ( mesh, problem, elgrp, elem, matrix, &
            vector, first, last, coefficients, oldvectors, elemmat, elemvec )
        end if
      case(3) ! 3D
        call linear_elastic_elem_3D ( mesh, problem, elgrp, elem, &
          matrix, vector, first, last, coefficients, oldvectors, elemmat, &
          elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in linear_elastic_elem:', &
        ' linear_elastic element not available for mesh%ndim = ', mesh%ndim
        stop
    end select

  end subroutine linear_elastic_elem


! Boundary element for a natural boundary for the Stokes equation

  subroutine linear_elastic_natboun ( mesh, problem, geom, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: geom, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    select case ( mesh%ndim )
      case(2) ! 2D and axisymmetric
        call linear_elastic_natboun_curve ( mesh, problem, geom, elem, matrix, &
          vector, first, last, coefficients, oldvectors, elemmat, elemvec )
      case(3) ! 3D
        call linear_elastic_natboun_surface ( mesh, problem, geom, elem, &
          matrix, vector, first, last, coefficients, oldvectors, elemmat, &
          elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in linear_elastic_natboun:', &
        ' linear_elastic element not available for mesh%ndim = ', mesh%ndim
        stop
    end select

  end subroutine linear_elastic_natboun

end module linear_elastic_elements_m

