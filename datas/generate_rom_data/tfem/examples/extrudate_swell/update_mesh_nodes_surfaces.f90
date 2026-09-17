! update mesh routine for free surfaces and material lines

module update_mesh_nodes_surfaces_m

  use tfem_m
  use hsl_ma41_m
  use hsl_ma57_m
  use laplace_elements_ale_m

  implicit none

contains

  subroutine update_mesh_nodes_expand ( mesh, problem, problemy, &
           disp1, disp2 )
!   input/output

!   coordinates of the mesh are updated after calling this routine
    type(mesh_t), intent(inout) :: mesh

!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problem, problemy

!   displacement of the particle
!    real(dp), dimension(:), intent(in) :: disp

!   displacement of the interface
    real(dp), dimension(:), intent(in) :: disp1, disp2

!   local definitions

    type(input_probdef_t) :: input_probdef
    type(sysmatrix_t) :: sysmatrix, sysmatrixy
    type(sysvector_t) :: rhsd, sol
    type(sysvector_t) :: rhsdy, soly
    type(coefficients_t) :: coefficients
    type(lu_ma41_t) :: lu_ma41
    type(solver_options_ma41_t) :: solver_options_ma41

!   constants

    integer, parameter :: &
      uintpl = 6,         & ! scalar interpolation
      gauss = 6,          & ! 8-point Gauss integration of tets
      gaussb = 3            ! 3 point integration of boundary elements

    real(dp), parameter :: &
      alpha = 1._dp     ! diffusion coefficient. Not used in laplace_elem_ale
                        ! Only for compatibility with poisson_elem.

    type(subscript_t) :: displ_1, displ_2, subsy

    call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

    coefficients%i = 0
    coefficients%i(1) = uintpl
    coefficients%i(10:11) = [ gauss, gaussb ]

    coefficients%r(1) = alpha
    coefficients%r(2:) = 0

    if ( .not. problem%created ) then

!     problem definition

      call create_input_probdef ( mesh, input_probdef )

      input_probdef%elementdof(1)%a = 1

      call define_essential ( mesh, input_probdef, curve1=1, curve2=4)

      call problem_definition ( input_probdef, mesh, problem )

      call delete ( input_probdef )

    end if

!   create solution and right-hand side

    call create ( problem, sol, rhsd )

!   create system matrix

    call create_sysmatrix_structure ( sysmatrix, mesh, problem )
    call create_sysmatrix_data ( sysmatrix )

    if ( .not. problemy%created ) then

!     problem definition

      call create_input_probdef ( mesh, input_probdef )

      input_probdef%elementdof(1)%a = 1

      call define_essential ( mesh, input_probdef, curve1=1)
      call define_essential ( mesh, input_probdef, curve1=3, curve2=4, &
        exclude=1)

      call problem_definition ( input_probdef, mesh, problemy )

      call delete ( input_probdef )

    end if

    call create_subscript ( mesh, problemy, subsy )
    call create_subscript ( mesh, problemy, displ_1, curves=[3] )
    call create_subscript ( mesh, problemy, displ_2, curves=[4] )

!   create solution and right-hand side

    call create ( problemy, soly, rhsdy )

!   create system matrix

    call create_sysmatrix_structure ( sysmatrixy, mesh, problemy )
    call create_sysmatrix_data ( sysmatrixy )


!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=laplace_elem_ale, coefficients=coefficients )

    call build_system ( mesh, problemy, sysmatrixy, rhsdy, &
      elemsub=laplace_elem_ale, coefficients=coefficients )

!   solve for x- and y-displacement

!     fill solution vector with zero essential boundary conditions

      sol%u = 0

!     fill solution vector with essential boundary conditions

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol,&
         rhsd )

    solver_options_ma41%integer_storage = 4.3
    solver_options_ma41%real_storage    = 4.3

!     solve the system
      call solve_system_ma41 ( sysmatrix, rhsd, sol, lu=lu_ma41, &
                  solver_options=solver_options_ma41 )

!    end do

    call delete ( lu_ma41 )

!   solve for y-displacement

!   fill solution vector with zero essential boundary conditions

    soly%u = 0

!   fill solution vector with essential boundary conditions

    soly%u(displ_1%s) = disp1
    soly%u(displ_2%s) = disp2

    call add_effect_of_essential_to_rhs ( problemy, sysmatrixy, soly, &
       rhsdy )

!   solve the system
    call solve_system_ma41 ( sysmatrixy, rhsdy, soly, lu=lu_ma41, &
                 solver_options=solver_options_ma41 )

    call delete ( lu_ma41 )

!   update mesh
    mesh%coor(:,1) = mesh%coor(:,1) + sol%u(problem%degfdperm(:,2))
    mesh%coor(:,2) = mesh%coor(:,2) + soly%u(subsy%s)

!   delete all data including all allocated memory

    call delete ( coefficients )
    call delete ( rhsd, sol )
    call delete ( sysmatrix )
    call delete (problem, problemy)

  end subroutine update_mesh_nodes_expand

! routine for updating the mesh nodes in the ALE scheme (displacement based)
! including free surfaces

  subroutine update_mesh_nodes_surfaces ( mesh, problem, problemy, problemz, &
           disp1, disp2, disp3, disp4 )
!   input/output

!   coordinates of the mesh are updated after calling this routine
    type(mesh_t), intent(inout) :: mesh

!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problem, problemy, problemz

!   displacement of the particle
!    real(dp), dimension(:), intent(in) :: disp

!   displacement of the interface
    real(dp), dimension(:,:), intent(in) :: disp1, disp2, disp3, disp4

!   local definitions

    type(input_probdef_t) :: input_probdef
    type(sysmatrix_t) :: sysmatrix, sysmatrixy, sysmatrixz
    type(sysvector_t) :: rhsd, sol
    type(sysvector_t) :: rhsdz, solz
    type(sysvector_t) :: rhsdy, soly
    type(coefficients_t) :: coefficients
    type(lu_ma41_t) :: lu_ma41
    type(solver_options_ma41_t) :: solver_options_ma41

!   constants

    integer, parameter :: &
      uintpl = 6,         & ! scalar interpolation
      gauss = 8,          & ! 8-point Gauss integration of tets
      gaussb = 3            ! 3 point integration of boundary elements

    real(dp), parameter :: &
      alpha = 1._dp     ! diffusion coefficient. Not used in laplace_elem_ale
                        ! Only for compatibility with poisson_elem.

    type(subscript_t) :: displ_surf1, displ_surf2, displ_surf3, displ_surf4, &
                         subsz, subsy, displ_surf1line, displ_surf2line, &
                         displ_surf3line, displ_surf4line

    call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

    coefficients%i = 0
    coefficients%i(1) = uintpl
    coefficients%i(10:11) = [ gauss, gaussb ]

    coefficients%r(1) = alpha
    coefficients%r(2:) = 0

    if ( .not. problem%created ) then

!     problem definition

      call create_input_probdef ( mesh, input_probdef )

      input_probdef%elementdof(1)%a = 1

      call define_essential ( mesh, input_probdef, surface1=1, surface2=6)

      call problem_definition ( input_probdef, mesh, problem )

      call delete ( input_probdef )

    end if

!   create solution and right-hand side

    call create ( problem, sol, rhsd )

!   create system matrix

    call create_sysmatrix_structure ( sysmatrix, mesh, problem )
    call create_sysmatrix_data ( sysmatrix )

    if ( .not. problemy%created ) then

!     problem definition

      call create_input_probdef ( mesh, input_probdef )

      input_probdef%elementdof(1)%a = 1

      call define_essential ( mesh, input_probdef, surface1=1, surface2=6)
      call define_essential ( mesh, input_probdef, surface1=8, surface2=11)

      call problem_definition ( input_probdef, mesh, problemy )

      call delete ( input_probdef )

    end if

    call create_subscript ( mesh, problemy, subsy )
    call create_subscript ( mesh, problemy, displ_surf2, surfaces=[9] )
    call create_subscript ( mesh, problemy, displ_surf1line, surfaces=[8] )
    call create_subscript ( mesh, problemy, displ_surf4, surfaces=[11] )
    call create_subscript ( mesh, problemy, displ_surf3line, surfaces=[10] )

!   create solution and right-hand side

    call create ( problemy, soly, rhsdy )

!   create system matrix

    call create_sysmatrix_structure ( sysmatrixy, mesh, problemy )
    call create_sysmatrix_data ( sysmatrixy )

    if ( .not. problemz%created ) then

!     problem definition

      call create_input_probdef ( mesh, input_probdef )

      input_probdef%elementdof(1)%a = 1

      call define_essential ( mesh, input_probdef, surface1=1, surface2=6)

      call define_essential ( mesh, input_probdef, surface1=8, surface2=11)

      call problem_definition ( input_probdef, mesh, problemz )

      call delete ( input_probdef )

    end if

    call create_subscript ( mesh, problemz, subsz )
    call create_subscript ( mesh, problemz, displ_surf1, surfaces=[8] )
    call create_subscript ( mesh, problemz, displ_surf3, surfaces=[10] )
    call create_subscript ( mesh, problemz, displ_surf2line, surfaces=[9] )
    call create_subscript ( mesh, problemz, displ_surf4line, surfaces=[11] )

!   create solution and right-hand side

    call create ( problemz, solz, rhsdz )

!   create system matrix

    call create_sysmatrix_structure_base ( sysmatrixz, mesh, problemz )
    call finalize_sysmatrix_structure ( sysmatrixz )

    call create_sysmatrix_data ( sysmatrixz )

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=laplace_elem_ale, coefficients=coefficients )

    call build_system ( mesh, problemz, sysmatrixz, rhsdz, &
      elemsub=laplace_elem_ale, coefficients=coefficients )

    call build_system ( mesh, problemy, sysmatrixy, rhsdy, &
      elemsub=laplace_elem_ale, coefficients=coefficients )

!   solve for x- and y-displacement

!     fill solution vector with zero essential boundary conditions

      sol%u = 0

!     fill solution vector with essential boundary conditions

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol,&
         rhsd )

    solver_options_ma41%integer_storage = 4.3
    solver_options_ma41%real_storage    = 4.3

!     solve the system
      call solve_system_ma41 ( sysmatrix, rhsd, sol, lu=lu_ma41, &
                  solver_options=solver_options_ma41 )

!    end do

    call delete ( lu_ma41 )

!   solve for y-displacement

!   fill solution vector with zero essential boundary conditions

    soly%u = 0

!   fill solution vector with essential boundary conditions

    soly%u(displ_surf2%s) = disp2(:,2)
    soly%u(displ_surf1line%s) = disp1(:,1)
    soly%u(displ_surf4%s) = disp4(:,2)
    soly%u(displ_surf3line%s) = disp3(:,1)

    call add_effect_of_essential_to_rhs ( problemy, sysmatrixy, soly, &
       rhsdy )

!   solve the system
    call solve_system_ma41 ( sysmatrixy, rhsdy, soly, lu=lu_ma41, &
                 solver_options=solver_options_ma41 )

    call delete ( lu_ma41 )

!   solve for z-displacement

!   fill solution vector with zero essential boundary conditions

    solz%u = 0

!   fill solution vector with essential boundary conditions

    solz%u(displ_surf1%s) = disp1(:,2)
    solz%u(displ_surf2line%s) = disp2(:,1)
    solz%u(displ_surf3%s) = disp3(:,2)
    solz%u(displ_surf4line%s) = disp4(:,1)

    call add_effect_of_essential_to_rhs ( problemz, sysmatrixz, solz, &
       rhsdz )

!   solve the system
    call solve_system_ma41 ( sysmatrixz, rhsdz, solz, lu=lu_ma41, &
                 solver_options=solver_options_ma41 )

!   update mesh
    mesh%coor(:,1) = mesh%coor(:,1) + sol%u(problem%degfdperm(:,2))
    mesh%coor(:,2) = mesh%coor(:,2) + soly%u(subsy%s)
    mesh%coor(:,3) = mesh%coor(:,3) + solz%u(subsz%s)

!   delete all data including all allocated memory

    call delete ( coefficients )
    call delete ( rhsd, sol )
    call delete ( rhsdz, solz )
    call delete ( sysmatrix )
    call delete ( sysmatrixz )
    call delete ( lu_ma41 )

  end subroutine update_mesh_nodes_surfaces

end module update_mesh_nodes_surfaces_m
