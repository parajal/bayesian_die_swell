
module update_mesh_nodes13_m

  use tfem_m
  use hsl_ma57_m
  use laplace_elements_ale_m

  implicit none

contains


! routine for updating the mesh nodes in the ALE scheme (displacement based)

  subroutine update_mesh_nodes ( mesh, problem, disp )

!   input/output

!   coordinates of the mesh are updated after calling this routine
    type(mesh_t), intent(inout) :: mesh

!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problem

!   displacement of the particle
    real(dp), dimension(:), intent(in) :: disp


!   local definitions

    type(input_probdef_t) :: input_probdef
    type(sysmatrix_t) :: sysmatrix
    type(sysvector_t), dimension(3) :: rhsd, sol
    type(coefficients_t) :: coefficients
    type(lu_ma57_t) :: lu_ma57

!   constants

    integer, parameter :: &
      uintpl = 6,         & ! scalar interpolation
      gauss = 8,          & ! 8-point Gauss integration of tets
      gaussb = 3            ! 3 point integration of boundary elements

    real(dp), parameter :: &
      alpha = 1._dp     ! diffusion coefficient. Not used in laplace_elem_ale
                        ! Only for compatibility with poisson_elem.

    integer :: i

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

      call define_essential ( mesh, input_probdef, surface1=1, surface2=6 )
      call define_essential ( mesh, input_probdef, surface1=7 )

      call problem_definition ( input_probdef, mesh, problem )

      call delete ( input_probdef )

    end if

!   create solution and right-hand side

    call create ( problem, sol, rhsd )

!   create system matrix

    call create_sysmatrix_structure ( sysmatrix, mesh, problem, &
      symmetric=.true. )
    call create_sysmatrix_data ( sysmatrix )

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, msysvector=rhsd, &
      elemsub=laplace_elem_ale, coefficients=coefficients )

!   solve

    do i = 1,3

!     fill solution vector with zero essential boundary conditions

      call fill_sysvector ( mesh, problem, sol(i), surface1=1, surface2=6, &
        value=0._dp )

!     fill solution vector with essential boundary conditions

      call fill_sysvector ( mesh, problem, sol(i), surface1=7, value=disp(i) )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol(i),&
         rhsd(i) )

!     solve the system
      call solve_system_ma57 ( sysmatrix, rhsd(i), sol(i), lu=lu_ma57 )

    end do

!   update mesh
    mesh%coor(:,1) = mesh%coor(:,1) + sol(1)%u(problem%degfdperm(:,2))
    mesh%coor(:,2) = mesh%coor(:,2) + sol(2)%u(problem%degfdperm(:,2))
    mesh%coor(:,3) = mesh%coor(:,3) + sol(3)%u(problem%degfdperm(:,2))

!   delete all data including all allocated memory

    call delete ( coefficients )
    call delete ( rhsd )
    call delete ( sysmatrix )
    call delete ( lu_ma57 )

  end subroutine update_mesh_nodes

end module update_mesh_nodes13_m
