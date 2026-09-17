! update mesh routine for drop8

module update_mesh_nodes8_m

  use tfem_m
  use hsl_ma57_m
  use laplace_elements_ale_m
!  use poisson_elements_m

  implicit none

contains


! routine for updating the mesh nodes in the ALE scheme (displacement based)

  subroutine update_mesh_nodes ( mesh, problem, disp )

!   input/output

!   coordinates of the mesh are updated after calling this routine
    type(mesh_t), intent(inout) :: mesh

!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problem

!   displacement of the droplet interface k
!   disp(k)%a(node,j), j=1,ndim
    type(real_array_2d_t), dimension(:), intent(in) :: disp

!   local definitions

    type(input_probdef_t) :: input_probdef
    type(sysmatrix_t) :: sysmatrix
    type(sysvector_t), dimension(size(disp(1)%a,2)) :: rhsd, sol
    type(coefficients_t) :: coefficients
    type(lu_ma57_t) :: lu

    type(subscript_t), dimension(size(disp)) :: disp7, disp8


!   constants

    integer, parameter :: &
      uintpl = 6,         & ! scalar interpolation
      gauss = 4,          & ! order of Gauss integration on tetrahedrons
      gaussb = 4            ! order of Gauss integration on triangles
                            !(boundaries)

    real(dp), parameter :: &
      alpha = 1._dp     ! diffusion coefficient. Not used in laplace_elem_ale
                        ! Only for compatibility with poisson_elem.

    integer :: i, j, nobj

    nobj=size(disp)

    call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

    coefficients%i = 0
    coefficients%i(1) = uintpl
    coefficients%i(10:11) = [ gauss, gaussb ]
    coefficients%i(40) = 3  ! "numeric" integration rules

    coefficients%r(1) = alpha
    coefficients%r(2:) = 0


    if ( .not. problem%created ) then

!     problem definition

      call create_input_probdef ( mesh, input_probdef )

      do i = 1, mesh%nelgrp
        input_probdef%elementdof(i)%a = 1
      end do

      call define_essential ( mesh, input_probdef, surface1=1, &
        surface2=6+2*nobj )

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
!      elemsub=poisson_elem, coefficients=coefficients )

!   creation of subscripts for interface displacement

    do j = 1, nobj
      call create_subscript ( mesh, problem, disp7(j), surfaces=[6+j] )
      call create_subscript ( mesh, problem, disp8(j), surfaces=[6+nobj+j] )
    end do

!   solve

    do i = 1, size(sol)

!     fill solution vector with zero essential boundary conditions

      call fill_sysvector ( mesh, problem, sol(i), surface1=1, surface2=6, &
        value=0._dp )

!     fill solution vector with essential boundary conditions

      do j = 1, nobj
        sol(i)%u(disp7(j)%s) = disp(j)%a(:,i)
        sol(i)%u(disp8(j)%s) = disp(j)%a(:,i)
      end do

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

  end subroutine update_mesh_nodes

end module update_mesh_nodes8_m
