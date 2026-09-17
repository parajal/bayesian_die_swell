! Copyright (C) 2012-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Solve a Poisson problem on a curve or surface:
!
!  - nabla^2_s c = 1/f - 1/g
!
! where f and g are functions given by values in the nodes:
!   f = monitor function representing the requested area of the elements
!   g = element area function
!
! The routines in this modules assume:
!    1) MA57 solver.

! NOTES:
! - Make a special version by modifying the source below to change 1).


module poisson_gd_ma57_m

  use tfem_m
  use hsl_ma57_m
  use interface_grid_deformation_m

  implicit none

contains


! routine for the Poisson problem on a closed curve or surface.

  subroutine poisson_gd_ma57 ( mesh, problem, coefficients, sysmatrix, &
    sol, f_monitor, solver_options, elementdof )

!   Mesh of the closed curve or surface
    type(mesh_t), intent(inout) :: mesh

!   The Poisson problem defined on the domain.
!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problem

!   The coefficients
    type(coefficients_t), intent(in) :: coefficients

!   The system matrix of the Poisson problem defined on the domain.
!   Is created if empty on call and kept at output
    type(sysmatrix_t), intent(inout) :: sysmatrix

!   Solution vector
!   NOTE: it goes only out.
    type(sysvector_t), intent(out) :: sol

!   The monitor function in the nodes of the mesh. The value of the monitor
!   function gives the relative size (length/area) of elements with respect to
!   each other. The absolute value is not important.
!   If it not present, it is set to 1 in all nodes.
    real(dp), dimension(:), intent(in), optional :: f_monitor

!   solver options for the MA57 solver
    type(solver_options_ma57_t), intent(in), optional :: solver_options

!   element degrees of freedom of the system vector for the poisson problem
!   elementdof(elgrp)%a(node)
    type(int_array_1d_t), dimension(:), intent(in), optional :: elementdof

!   Local definitions

    type(input_probdef_t) :: input_probdef
    type(oldvectors_t) :: oldvectors
    type(sysvector_t) :: rhsd

!   The area and monitor function (both unscaled and scaled)
    type(vector_t), target :: garea, fmon

    integer :: elgrp

    real(dp) :: resultsum(2)


!   problem definition

    if ( .not. problem%created ) then

      call create_input_probdef ( mesh, input_probdef, nvec=1 )

      if ( present(elementdof) ) then
        input_probdef%elementdof = elementdof
        do elgrp = 1, mesh%nelgrp
          input_probdef%vec_elementdof(elgrp)%a(:,1) = 1
        end do
      else
        do elgrp = 1, mesh%nelgrp
          input_probdef%elementdof(elgrp)%a = 1
          input_probdef%vec_elementdof(elgrp)%a(:,1) = 1
        end do
      end if

      call define_essential ( mesh, input_probdef, element=1, elnode=1 )

      call problem_definition ( input_probdef, mesh, problem )

      call delete ( input_probdef )

    end if

!   create system vectors (solution and right-hand side)

    call create ( problem, sol, rhsd )

!   fill solution vector with zero essential boundary conditions

    call fill_sysvector ( mesh, problem, sol, element=1, elnode=1, value=0._dp )

!   create oldvectors

    call create_oldvectors ( oldvectors, nvec=2 )


!   area function

    call create_vector ( problem, garea, vec=1 )

!   derive area function

    call derive_vector ( mesh, problem, garea, elemsub=g_area_deriv, &
      coefficients=coefficients, oldvectors=oldvectors )

    oldvectors%v(1)%p => garea

    call integrate ( mesh, problem, resultsum, &
      elemsub=integrate_inverse_elem, coefficients=coefficients, &
      oldvectors=oldvectors )

!   scaled g such that
!         /
!         | 1/g dx = area of Gamma
!         /Gamma
!   or average(1/g)=1.

    garea%u = resultsum(1)/resultsum(2)*garea%u

!   monitor function

    call create_vector ( problem, fmon, vec=1 )

    if ( present(f_monitor) ) then

      fmon%u = f_monitor

      oldvectors%v(1)%p => fmon

      call integrate ( mesh, problem, resultsum, &
        elemsub=integrate_inverse_elem, &
        coefficients=coefficients, oldvectors=oldvectors )

!     scaled f such that
!           /
!           | 1/f dx = area of Gamma
!           /Gamma
!     or average(1/f)=1.

      fmon%u = resultsum(1)/resultsum(2)*fmon%u

    else

      fmon%u = 1._dp

    end if

!   Note that
!         /
!         | ( 1/f - 1/g ) dx = 0
!         /Omega
!   for consistency with the zero flux on the boundary.


    if ( .not. sysmatrix%allocated_data ) then

!     create system matrix

      call create_sysmatrix_structure ( sysmatrix, mesh, problem, &
        symmetric=.true. )

      call create_sysmatrix_data ( sysmatrix )

    end if

    oldvectors%v(1)%p => garea
    oldvectors%v(2)%p => fmon

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=poisson_gd_elem, coefficients=coefficients, &
      oldvectors=oldvectors )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   MA57 solver

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options )

!   delete some remaining structures

    call delete ( garea, fmon )
    call delete ( rhsd )
    call delete ( oldvectors )

  end subroutine poisson_gd_ma57

end module poisson_gd_ma57_m
