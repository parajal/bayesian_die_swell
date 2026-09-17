
! Copyright (C) 2005-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the Stokes equation
!
!    - div( eta(nabla u+nabla u^T) ) + nabla p = f
!      div u = 0
!

module stokes_elements_m

  use tfem_elem_m
  use stokes_elements_generic_m
  use stokes_elements_embedded_boundary_m
  use stokes_elements_embedded_interface_m
  use stokes_elements_embedded_particle_m
  use stokes_elements_2D_curve_m
  use stokes_elements_3D_surface_m
  use stokes_elements_stabilized_m
  use stokes_elements_preconditioning_m
  use stokes_elements_integrate_m

  implicit none

contains


! Internal element routine for the Stokes Equation

  subroutine stokes_elem ( mesh, problem, elgrp, elem, matrix, vector, &
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
          call stokes_elem_2D_vel3D ( mesh, problem, elgrp, elem, matrix, &
            vector, first, last, coefficients, oldvectors, elemmat, elemvec )
        else
          call stokes_elem_2D ( mesh, problem, elgrp, elem, matrix, vector, &
            first, last, coefficients, oldvectors, elemmat, elemvec )
        end if
      case(3) ! 3D
        call stokes_elem_3D ( mesh, problem, elgrp, elem, matrix, vector, &
          first, last, coefficients, oldvectors, elemmat, elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in stokes_elem:', &
        ' stokes element not available for mesh%ndim = ', mesh%ndim
        stop
    end select

  end subroutine stokes_elem


! Boundary element for a natural boundary for the Stokes equation

  subroutine stokes_natboun ( mesh, problem, geom, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

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
        call stokes_natboun_curve ( mesh, problem, geom, elem, matrix, &
          vector, first, last, coefficients, oldvectors, elemmat, elemvec )
      case(3) ! 3D
        call stokes_natboun_surface ( mesh, problem, geom, elem, matrix, &
          vector, first, last, coefficients, oldvectors, elemmat, elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in stokes_natboun:', &
        ' stokes element not available for mesh%ndim = ', mesh%ndim
        stop
    end select

  end subroutine stokes_natboun


! Element for the flowrate constraint

  subroutine stokes_constr_flowr ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    select case ( mesh%ndim )
      case(2) ! 2D and axisymmetric
        call stokes_constr_flowr_curve ( mesh, problem, constr, elem, node, &
          matrix, vector, first, last, coefficients, oldvectors, elemmat, &
          elemmat2, elemmatadd, elemvec, elemvecadd )
      case(3) ! 3D
        call stokes_constr_flowr_surface ( mesh, problem, constr, elem, node, &
          matrix, vector, first, last, coefficients, oldvectors, elemmat, &
          elemmat2, elemmatadd, elemvec, elemvecadd )
      case default
        write(*,'(/a/a,i0/)') 'Error in stokes_constr_flowr:', &
        ' stokes element not available for mesh%ndim = ', mesh%ndim
        stop
    end select


  end subroutine stokes_constr_flowr


! Element for the constraints (connection through elements)

  subroutine stokes_constr_elem_conn ( mesh, problem, constr, elem, &
    node, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    select case ( mesh%ndim )
      case(2) ! 2D and axisymmetric
        call stokes_constr_elem_conn_curve ( mesh, problem, constr, elem, &
          node, matrix, vector, first, last, coefficients, oldvectors, &
          elemmat, elemmat2, elemmatadd, elemvec, elemvecadd )
      case(3) ! 3D
        call stokes_constr_elem_conn_surface ( mesh, problem, constr, elem, &
          node, matrix, vector, first, last, coefficients, oldvectors, &
          elemmat, elemmat2, elemmatadd, elemvec, elemvecadd )
      case default
        write(*,'(/a/a,i0/)') 'Error in stokes_constr_elem_conn:', &
        ' stokes element not available for mesh%ndim = ', mesh%ndim
        stop
    end select

  end subroutine stokes_constr_elem_conn


! Element for the flow rate (for integrate_boundary_elements in postprocessing)

  subroutine stokes_flowrate ( mesh, problem, curve, elem, first, &
    last, coefficients, oldvectors, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    select case ( mesh%ndim )
      case(2) ! 2D and axisymmetric
        call stokes_flowrate_curve ( mesh, problem, curve, elem, first, &
          last, coefficients, oldvectors, elemvec )
      case(3) ! 3D
        call stokes_flowrate_surface ( mesh, problem, curve, elem, first, &
          last, coefficients, oldvectors, elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in stokes_flowrate:', &
        ' stokes element not available for mesh%ndim = ', mesh%ndim
        stop
    end select

  end subroutine stokes_flowrate


! Element for the drag force using integrate_boundary_elements in
! the postprocessing_m module.
!
! oldvectors:
!  v(1) = "viscous stress_tensor" derived from subroutine "stokes_stress_tensor"
!  v(2) = "pressure" derived from subroutine "stokes_pressure"

  subroutine stokes_drag ( mesh, problem, curve, elem, first, last, &
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
        call stokes_drag_curve ( mesh, problem, curve, elem, first, last, &
          coefficients, oldvectors, elemvec )
      case(3) ! 3D
        call stokes_drag_surface ( mesh, problem, curve, elem, first, last, &
          coefficients, oldvectors, elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in stokes_drag:', &
        ' stokes element not available for mesh%ndim = ', mesh%ndim
        stop
    end select

  end subroutine stokes_drag

end module stokes_elements_m

