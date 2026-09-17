! update mesh routine for fountain_flow1

module update_mesh_nodes1_m

  use tfem_m
  use hsl_ma57_m
!  use laplace_elements_ale_m
  use poisson_elements_m

  implicit none

contains


! routine for updating the mesh nodes in the ALE scheme (displacement based)

  subroutine update_mesh_nodes ( mesh, problem, disp, el, eg )

!   input/output

!   coordinates of the mesh are updated after calling this routine
    type(mesh_t), intent(inout) :: mesh

!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problem

!   displacement of the moving boundary
!   disp(node,j), j=1,ndim
    real(dp), dimension(:,:), intent(in) :: disp

!   nodal points of the ends of moving boundary:
!     el(:) local nodes
!     eg(:) global nodes
    integer, dimension(2), intent(in) :: el, eg


!   local definitions

    type(input_probdef_t) :: input_probdef
    type(sysmatrix_t) :: sysmatrix
    type(sysvector_t), dimension(size(disp,2)) :: rhsd, sol
    type(coefficients_t) :: coefficients
    type(lu_ma57_t) :: lu

    type(subscript_t) :: disp1, disp12, disp14


!   constants

    integer, parameter :: &
      uintpl = 8,         & ! scalar interpolation (Q2)
      gauss = 3,          & ! 3x3-point Gauss integration of quads
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

      do i = 1, mesh%nelgrp
        input_probdef%elementdof(i)%a = 1
      end do

      call define_essential ( mesh, input_probdef, curve1=17 )

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
!      elemsub=laplace_elem_ale, coefficients=coefficients )
      elemsub=poisson_elem, coefficients=coefficients )

!   creation of subscripts for interface displacement

    call create_subscript ( mesh, problem, disp1, curves=[1] )
    call create_subscript ( mesh, problem, disp12, curves=[12] )
    call create_subscript ( mesh, problem, disp14, curves=[14] )

!   solve

    do i = 1, size(sol)

!     fill solution vector with essential boundary conditions

      call fill_sysvector ( mesh, problem, sol(i), curve1=17, value=0._dp )

      if ( i == 1 ) then
!       displace the wall nodes proportionally
        sol(i)%u(disp1%s) = disp(el(1),1) * mesh%coor(mesh%curves(1)%nodes,1) &
                            / mesh%coor(eg(1),1)
        sol(i)%u(disp12%s) = disp(el(2),1) * mesh%coor(mesh%curves(12)%nodes,1)&
                            / mesh%coor(eg(2),1)
      end if

      sol(i)%u(disp14%s) = disp(:,i) ! moving boundary

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol(i),&
         rhsd(i) )

!     solve the system
      call solve_system_ma57 ( sysmatrix, rhsd(i), sol(i), lu=lu )

    end do

!   update mesh
    do i = 1, size(sol)
      mesh%coor(:,i) = mesh%coor(:,i) + sol(i)%u(problem%degfdperm(:,2))
    end do

!   delete all data including all allocated memory

    call delete ( coefficients )
    call delete ( rhsd )
    call delete ( sysmatrix )
    call delete ( lu )
    call delete ( disp1, disp12, disp14 )

  end subroutine update_mesh_nodes

end module update_mesh_nodes1_m
