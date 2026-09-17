! update mesh routine for drop13

module update_mesh_nodes13_m

  use tfem_m
  use hsl_ma57_m
  use laplace_elements_ale_m
! use poisson_elements_m

  implicit none

contains


! routine for updating the mesh nodes in the ALE scheme (displacement based)

  subroutine update_mesh_nodes ( mesh, problem, disp )

!   input/output

!   coordinates of the mesh are updated after calling this routine
    type(mesh_t), intent(inout) :: mesh

!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problem

!   displacement of the droplet interface k:
    type(real_array_2d_t), intent(inout) :: disp
    type(real_array_2d_t), dimension(3) :: disp_symaxis

!   local definitions

    type(input_probdef_t) :: input_probdef
    type(sysmatrix_t) :: sysmatrix
    type(sysvector_t), dimension(size(disp%a,2))  :: rhsd, sol
    type(coefficients_t) :: coefficients
    type(lu_ma57_t) :: lu

    type(subscript_t):: disp5, disp6
    type(subscript_t), dimension(3) :: disp1


!   constants

    integer, parameter :: &
      uintpl = 6,         & ! scalar interpolation
      gauss = 6,          & ! 6-point Gauss integration of triangles
      gaussb = 3            ! 3 point integration of boundary elements

    real(dp), parameter :: &
      alpha = 1._dp     ! diffusion coefficient. Not used in laplace_elem_ale
                        ! Only for compatibility with poisson_elem.

    integer :: i, lastnode
    real(dp) :: length, x1, x2, x



    call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

    coefficients%i = 0
    coefficients%i(1) = uintpl
    coefficients%i(10:11) = [ gauss, gaussb ]
    coefficients%i(23) = 1 !axisymmetric

    coefficients%r(1) = alpha
    coefficients%r(2:) = 1

    lastnode = mesh%curves(7)%nnodes

    if ( .not. problem%created ) then

!     problem definition

      call create_input_probdef ( mesh, input_probdef )

      do i = 1, mesh%nelgrp
        input_probdef%elementdof(i)%a = 1
      end do

      call define_essential ( mesh, input_probdef, curve1=1, curve2=8 )

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

!   displacement of the curves on symmetry axis
    allocate(disp_symaxis(1)%a(mesh%curves(1)%nnodes,2))
    disp_symaxis(1)%a = 0._dp
    x1 = mesh%coor(mesh%curves(7)%nodes(lastnode),1)
    x2 = mesh%coor(mesh%curves(1)%nodes(1),1)

    length = abs(x2 - x1)
    do i = 2, size(disp_symaxis(1)%a,1)-1
      disp_symaxis(1)%a(i,1) = abs(mesh%coor(mesh%curves(1)%nodes(i),1) - x2) &
        / length * disp%a(lastnode,1)
    end do
    disp_symaxis(1)%a(  size(disp_symaxis(1)%a,1),1) = disp%a(lastnode,1)


    allocate(disp_symaxis(2)%a(mesh%curves(8)%nnodes,2))
    disp_symaxis(2)%a = 0._dp
    x1= mesh%coor(mesh%curves(7)%nodes(lastnode),1)
    x2= mesh%coor(mesh%curves(7)%nodes(1),1)
    length = abs(x1-x2)
    do i = 1, size(disp_symaxis(2)%a,1)
      x = mesh%coor(mesh%curves(8)%nodes(i),1)
      disp_symaxis(2)%a(i,1) = (abs( x - x2 ) * disp%a(lastnode,1) + &
        abs( x - x1 ) * disp%a(1,1)) / length
    end do

    disp_symaxis(2)%a(1,1) = disp%a(lastnode,1)
    disp_symaxis(2)%a(mesh%curves(8)%nnodes,1) = &
       disp%a(1,1)

    allocate(disp_symaxis(3)%a(mesh%curves(3)%nnodes,2))
    disp_symaxis(3)%a = 0._dp
    x1 = mesh%coor(mesh%curves(7)%nodes(1),1)
    x2 = mesh%coor(mesh%curves(3)%nodes(mesh%curves(3)%nnodes),1)
    length = abs( x2 - x1)
    do i = 2, size(disp_symaxis(3)%a,1)-1
      disp_symaxis(3)%a(i,1) = abs(mesh%coor(mesh%curves(3)%nodes(i),1)-x2) &
        / length * disp%a(1,1)
    end do
    disp_symaxis(3)%a(1,1) = disp%a(1,1)

!   creation of subscripts for interface displacement

    call create_subscript ( mesh, problem, disp6, curves=[2] )
    call create_subscript ( mesh, problem, disp5, curves=[7] )

    call create_subscript ( mesh, problem, disp1(1), curves=[1] )
    call create_subscript ( mesh, problem, disp1(2), curves=[8] )
    call create_subscript ( mesh, problem, disp1(3), curves=[3] )

!   solve
    do i = 1, size(sol)

!     fill solution vector with zero essential boundary conditions
      if ( i == 1) then
        call fill_sysvector ( mesh, problem, sol(i), curve1=4, curve2=6, &
          value=0._dp )
      else
        call fill_sysvector ( mesh, problem, sol(i), curve1=1, curve2=8, &
          excludecurves=[7,2], value=0._dp )
      end if

!     fill solution vector with essential boundary conditions

      sol(i)%u(disp1(1)%s) = disp_symaxis(1)%a(:,i)
      sol(i)%u(disp1(2)%s) = disp_symaxis(2)%a(:,i)
      sol(i)%u(disp1(3)%s) = disp_symaxis(3)%a(:,i)

      sol(i)%u(disp5%s) = disp%a(:,i)
      sol(i)%u(disp6%s) = disp%a(:,i)

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

end module update_mesh_nodes13_m
