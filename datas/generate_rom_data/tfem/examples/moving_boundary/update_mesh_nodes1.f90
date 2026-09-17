! update mesh routine for moving_boundary1

module update_mesh_nodes1_m

  use tfem_m
  use hsl_ma57_m
!  use laplace_elements_ale_m
  use poisson_elements_m
  use set_optional_m

  implicit none

contains


! routine for updating the mesh nodes in the ALE scheme (displacement based)

  subroutine update_mesh_nodes ( mesh, problem, disp, el, eg, cf, df )

!   input/output

!   coordinates of the mesh are updated after calling this routine
    type(mesh_t), intent(inout) :: mesh

!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problem

!   displacement of the moving boundary
!   disp(node)
    real(dp), dimension(:), intent(in) :: disp

!   nodal points of the ends of moving boundary:
!     el(:) local nodes
!     eg(:) global nodes
    integer, dimension(2), intent(in) :: el, eg

!   curve number of the free surface
    integer, intent(in) :: cf

!   if present, the coordinate direction of the motion of the free boundary
!   default = 1
    integer, optional, intent(in) :: df


!   local definitions

    type(input_probdef_t) :: input_probdef
    type(sysmatrix_t) :: sysmatrix
    type(sysvector_t) :: rhsd, sol
    type(coefficients_t) :: coefficients

    type(subscript_t) :: disp1, disp6, disp_free


!   constants

    integer, parameter :: &
      uintpl = 8,         & ! scalar interpolation (Q2)
      gauss = 3,          & ! 3x3-point Gauss integration of quads
      gaussb = 3            ! 3 point integration of boundary elements

    real(dp), parameter :: &
      alpha = 1._dp     ! diffusion coefficient. Not used in laplace_elem_ale
                        ! Only for compatibility with poisson_elem.

    integer :: i, dir_free


    dir_free = set_optional ( variable=df, default=1 )


    call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

    coefficients%i = 0
    coefficients%i(1) = uintpl
    coefficients%i(10:11) = [ gauss, gaussb ]

    coefficients%r(1) = alpha
    coefficients%r(2:) = 0


    if ( .not. problem%created ) then

!     problem definition

      call create_input_probdef ( mesh, input_probdef )

      do i = 1, mesh%nelgrp
        input_probdef%elementdof(i)%a = 1
      end do

      call define_essential ( mesh, input_probdef, curve1=11 )

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

    call build_system ( mesh, problem, sysmatrix, sysvector=rhsd, &
!      elemsub=laplace_elem_ale, coefficients=coefficients )
      elemsub=poisson_elem, coefficients=coefficients )

!   creation of subscripts for interface displacement

    call create_subscript ( mesh, problem, disp1, curves=[1] )
    call create_subscript ( mesh, problem, disp6, curves=[6] )
    call create_subscript ( mesh, problem, disp_free, curves=[cf] )

!   fill solution vector with essential boundary conditions

    call fill_sysvector ( mesh, problem, sol, curve1=11, value=0._dp )

!   displace the nodes on the walls proportionally
    sol%u(disp1%s) = disp(el(1)) * &
                        mesh%coor(mesh%curves(1)%nodes,dir_free) &
                        / mesh%coor(eg(1),dir_free)
    sol%u(disp6%s) = disp(el(2)) * &
                        mesh%coor(mesh%curves(6)%nodes,dir_free) &
                        / mesh%coor(eg(2),dir_free)

    sol%u(disp_free%s) = disp ! moving boundary

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve the system
    call solve_system_ma57 ( sysmatrix, rhsd, sol )

!   update mesh
    mesh%coor(:,dir_free) = mesh%coor(:,dir_free) + &
                                             sol%u(problem%degfdperm(:,2))

!   delete all data including all allocated memory

    call delete ( coefficients )
    call delete ( rhsd )
    call delete ( sysmatrix )
    call delete ( disp1, disp6, disp_free )

  end subroutine update_mesh_nodes

end module update_mesh_nodes1_m
