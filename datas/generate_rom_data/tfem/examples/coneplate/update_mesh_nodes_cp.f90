! update mesh routine for cone plate

module update_mesh_nodes_cp_m

  use tfem_m
  use hsl_ma57_m
  use hsl_ma41_m
  use laplace_elements_ale_m

  implicit none

contains

! routine for updating the mesh nodes in the ALE scheme (displacement based)
! including free surfaces

  subroutine update_mesh_nodes ( mesh, problemx, problemy, disp1, disp2, &
    disp3, disp4, disp5 )
!   input/output

!   coordinates of the mesh are updated after calling this routine
    type(mesh_t), intent(inout) :: mesh

!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problemx, problemy

!   displacement of the interface
    real(dp), dimension(:), intent(in) :: disp1, disp2
!   displacement cone
    real(dp), dimension(:), intent(in) :: disp3, disp4, disp5

!   local definitions

    type(input_probdef_t) :: input_probdef
    type(sysmatrix_t) :: sysmatrixx, sysmatrixy
    type(sysvector_t) :: rhsdx, solx
    type(sysvector_t) :: rhsdy, soly
    type(coefficients_t) :: coefficients
    type(lu_ma41_t) :: lu_ma41
    type(solver_options_ma41_t) :: solver_options_ma41

!   constants

    integer, parameter :: &
      uintpl = 6,         & ! scalar interpolation
      gauss = 5,          & ! 8-point Gauss integration of tets
      gaussb = 5            ! 3 point integration of boundary elements

    real(dp), parameter :: &
      alpha = 1._dp     ! diffusion coefficient. Not used in laplace_elem_ale
                        ! Only for compatibility with poisson_elem.

    type(subscript_t) :: freedispx, freedispy, subsx, subsy, conedispx, &
                         conedispy, platedispy


    call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

    coefficients%i = 0
    coefficients%i(1) = uintpl
    coefficients%i(10:11) = [ gauss, gaussb ]
    coefficients%i(40) = 3

    coefficients%r(1) = alpha
    coefficients%r(2:) = 0

    if ( .not. problemx%created ) then

!     problem definition

      call create_input_probdef ( mesh, input_probdef )

      input_probdef%elementdof(1)%a = 1

      call define_essential ( mesh, input_probdef, curve1=1 )
      call define_essential ( mesh, input_probdef, curve1=2 )
      call define_essential ( mesh, input_probdef, curve1=3 )

      call problem_definition ( input_probdef, mesh, problemx )

      call delete ( input_probdef )

    end if

    call create_subscript ( mesh, problemx, subsx )
    call create_subscript ( mesh, problemx, freedispx, curves=[3] )
    call create_subscript ( mesh, problemx, conedispx, curves=[2] )

!   create solution and right-hand side

    call create ( problemx, solx, rhsdx )

!   create system matrix

    call create_sysmatrix_structure ( sysmatrixx, mesh, problemx )
    call create_sysmatrix_data ( sysmatrixx )

    if ( .not. problemy%created ) then

!     problem definition

      call create_input_probdef ( mesh, input_probdef )

      input_probdef%elementdof(1)%a = 1

      call define_essential ( mesh, input_probdef, curve1=1 )
      call define_essential ( mesh, input_probdef, curve1=2 )
      call define_essential ( mesh, input_probdef, curve1=3 )

      call problem_definition ( input_probdef, mesh, problemy )

      call delete ( input_probdef )

    end if

    call create_subscript ( mesh, problemy, subsy )
    call create_subscript ( mesh, problemy, freedispy, curves=[3] )
    call create_subscript ( mesh, problemy, platedispy, curves=[1] )
    call create_subscript ( mesh, problemy, conedispy, curves=[2] )

!   create solution and right-hand side

    call create ( problemy, soly, rhsdy )

!   create system matrix

    call create_sysmatrix_structure ( sysmatrixy, mesh, problemy )
    call create_sysmatrix_data ( sysmatrixy )


!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problemx, sysmatrixx, rhsdx, &
      elemsub=laplace_elem_ale, coefficients=coefficients )

    call build_system ( mesh, problemy, sysmatrixy, rhsdy, &
      elemsub=laplace_elem_ale, coefficients=coefficients )

!   solve for x- and y-displacement

!   fill solution vector with zero essential boundary conditions

    solx%u = 0

    solx%u(freedispx%s) = disp1
    solx%u(conedispx%s) = disp3

!   fill solution vector with essential boundary conditions

    call add_effect_of_essential_to_rhs ( problemx, sysmatrixx, solx,&
      rhsdx )

    solver_options_ma41%integer_storage = 4.3
    solver_options_ma41%real_storage    = 4.3

!   solve the system
    call solve_system_ma41 ( sysmatrixx, rhsdx, solx, lu=lu_ma41, &
      solver_options=solver_options_ma41 )

    call delete ( lu_ma41 )

!   solve for y-displacement

!   fill solution vector with zero essential boundary conditions

    soly%u = 0

!   fill solution vector with essential boundary conditions

    soly%u(freedispy%s) = disp2
    soly%u(conedispy%s) = disp4
    soly%u(platedispy%s) = disp5

    call add_effect_of_essential_to_rhs ( problemy, sysmatrixy, soly, &
       rhsdy )

!   solve the system
    call solve_system_ma41 ( sysmatrixy, rhsdy, soly, lu=lu_ma41, &
                 solver_options=solver_options_ma41 )

    call delete ( lu_ma41 )

!   update mesh
    mesh%coor(:,1) = mesh%coor(:,1) + solx%u(subsx%s)
    mesh%coor(:,2) = mesh%coor(:,2) + soly%u(subsy%s)

!   delete all data including all allocated memory

    call delete ( coefficients )
    call delete ( rhsdx, solx )
    call delete ( rhsdy, soly )
    call delete ( sysmatrixx, sysmatrixy )

  end subroutine update_mesh_nodes

end module update_mesh_nodes_cp_m

