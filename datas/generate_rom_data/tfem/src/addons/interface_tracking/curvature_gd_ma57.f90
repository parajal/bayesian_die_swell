! Copyright (C) 2014-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Solve a projection problem on a curve or surface given by
!
!    c = n nabla_s.n
!    -   -         -
! where n is the normal vector on the surface or curve. The rhs is partially
! integrated like the surface tension force in the momentum balance.
!
! The routines in this modules assume:
!    1) Closed surface/curve (2D/3D) or an axisymmetric problem.
!    2) MA57 solver.

! NOTES:
! - Make a special version by modifying the source below to change 1) and 2).


module curvature_gd_ma57_m

  use tfem_m
  use hsl_ma57_m
  use interface_grid_deformation_m
  use set_optional_m

  implicit none

contains


! routine for the curvature problem on a closed curve or surface.
! Note: curvature = |nabla_s.n|, which is |1/R| for a curve in 2D and
! |1/R1+1/R2| for a surface or a curve in an axisymmetric problem.
! Note, that the sign of |nabla_s.n| is not recovered.

  subroutine curvature_gd_ma57 ( mesh, problem, coefficients, sysmatrix, &
    curvature, solver_options, square, elementdof )

!   Mesh of the closed curve or surface
    type(mesh_t), intent(inout) :: mesh

!   The curvature problem defined on the domain.
!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problem

!   The coefficients
    type(coefficients_t), intent(in) :: coefficients

!   The system matrix of the curvature problem defined on the domain.
!   Is created if empty on call and kept at output
    type(sysmatrix_t), intent(inout) :: sysmatrix

!   Curvature in all nodes of the mesh.
!   NOTE: it goes only out.
    real(dp), dimension(:), intent(out) :: curvature

!   solver options for the MA57 solver
    type(solver_options_ma57_t), intent(in), optional :: solver_options

!   if square=.true. the square of the curvature is output
!   default = .false.
    logical, intent(in), optional :: square

!   element degrees of freedom of the system vector for the curvature problem
!   elementdof(elgrp)%a(node)
    type(int_array_1d_t), dimension(:), intent(in), optional :: elementdof


!   Local definitions

    type(input_probdef_t) :: input_probdef
    type(sysvector_t) :: rhsd(mesh%ndim), sol(mesh%ndim)
    type(lu_ma57_t) :: lu

    logical :: lsquare
    integer :: elgrp, i


    lsquare = set_optional ( variable=square, default=.false. )

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

      call problem_definition ( input_probdef, mesh, problem )

      call delete ( input_probdef )

    end if

!   create system vectors (solution and right-hand side)

    call create ( problem, sol, rhsd )

!   system matrix

    if ( .not. sysmatrix%allocated_data ) then

!     create system matrix

      call create_sysmatrix_structure ( sysmatrix, mesh, problem, &
        symmetric=.true. )

      call create_sysmatrix_data ( sysmatrix )

    end if

    call build_system ( mesh, problem, sysmatrix, msysvector=rhsd, &
      elemsub=curvature_gd_elem, coefficients=coefficients )

!   MA57 solver

    do i = 1, size(sol)
      call solve_system_ma57 ( sysmatrix, rhsd(i), sol(i), &
        solver_options=solver_options, lu=lu )
    end do

!   curvature

    curvature = sol(1)%u ** 2
    do i = 2, mesh%ndim
      curvature = curvature + sol(i)%u ** 2
    end do
    if ( .not. lsquare ) curvature = sqrt(curvature)

!   delete some remaining structures

    call delete ( lu )
    call delete ( rhsd, sol )

  end subroutine curvature_gd_ma57

end module curvature_gd_ma57_m
