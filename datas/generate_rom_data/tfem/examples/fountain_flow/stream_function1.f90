! Streamfunction module for fountain_flow1

! Internal streamfunction element routine using the Poisson Equation:
!
!   - nabla^2 psi = omega
!

! Boundary element for a natural boundary on a curve for the streamfunction
! equation (Poisson equation):
!
!    dpsidn = -v * nx + u * ny   (= -tangential velocity)
!
! since dpsidx = -v, dpsidy = u.
!

! Dirichlet in a single point (P1) psi=0

module stream_function1_m

  use tfem_m
  use streamfunction_elements_m
  use hsl_ma57_m

  implicit none

contains


! routine for the streamfunction

  subroutine stream_function ( mesh, problem, coefficients, sol_vel, vel, sol )

!   input/output

    type(mesh_t), intent(inout) :: mesh

!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problem

    type(coefficients_t), intent(in) :: coefficients

!   solution containing the velocity
    type(sysvector_t), intent(in) :: sol_vel

!   subscript for the velocity vector
    type(subscript_t), intent(in) :: vel

!   solution containing the streamfunction
    type(sysvector_t), intent(inout) :: sol


!   local definitions

    type(input_probdef_t) :: input_probdef
    type(sysmatrix_t) :: sysmatrix
    type(sysvector_t) :: rhsd
    type(vector_t), target :: velocity
    type(oldvectors_t) :: oldvectors
    type(solver_options_ma57_t) :: solver_options

    integer, parameter :: crvboun = 17
    integer :: i


    if ( .not. problem%created ) then

!     problem definition

      call create_input_probdef ( mesh, input_probdef, nvec=1 )

      do i = 1, mesh%nelgrp
        input_probdef%elementdof(i)%a = 1
        input_probdef%vec_elementdof(i)%a = 2
      end do

      call define_essential ( mesh, input_probdef, point=1 ) ! Dirichlet in P1

      call problem_definition ( input_probdef, mesh, problem )

      call delete ( input_probdef )

!     create solution
      call create ( problem, sol )

    end if

!   create right-hand side

    call create ( problem, rhsd )

!   fill solution vector with essential boundary conditions

    call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )

!   create system matrix

    call create_sysmatrix_structure ( sysmatrix, mesh, problem, &
      symmetric=.true. )
    call create_sysmatrix_data ( sysmatrix )

!   create velocity vector

    call create_vector ( problem, velocity, vec=1 )

    velocity%u = sol_vel%u(vel%s)

!   build (assemble) matrix and vector from elements

    call create_oldvectors ( oldvectors, nvec=1 )
    oldvectors%v(1)%p => velocity

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=streamfunction_elem, coefficients=coefficients, &
      oldvectors=oldvectors )

    call add_boundary_elements ( mesh, problem, rhsd, curve=crvboun, &
      elemsub=streamfunction_natboun_curve, coefficients=coefficients, &
      oldvectors=oldvectors )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    solver_options%integer_storage=1.4

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options )

!   delete all data including all allocated memory

    call delete ( rhsd )
    call delete ( sysmatrix )
    call delete ( velocity )
    call delete ( oldvectors )

  end subroutine stream_function

end module stream_function1_m
