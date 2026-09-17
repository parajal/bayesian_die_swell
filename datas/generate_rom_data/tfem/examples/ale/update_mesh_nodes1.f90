! update mesh routine for ale1

module update_mesh_nodes1_m

  use tfem_m
  use hsl_ma57_m
  use laplace_elements_ale_m

  implicit none

contains


! routine for updating the mesh nodes in the ALE scheme (velocity based)

  subroutine update_mesh_nodes ( mesh, problem, velp, dt, sol, method )

!   input/output

!   coordinates of the mesh are updated after calling this routine
    type(mesh_t), intent(inout) :: mesh

!   Is created if empty on call and kept at output
    type(problem_t), intent(inout) :: problem

!   velocity of the particle and the time step
    real(dp), intent(in) :: velp(2), dt

!   mesh velocity
    type(sysvector_t), dimension(:), intent(inout) :: sol

!   subscript

!   Time integration:
!     method=1  Euler explicit
!     method=2  Second-order Adams-Bashforth
    integer, intent(in) :: method


!   local definitions

    type(input_probdef_t) :: input_probdef
    type(sysmatrix_t) :: sysmatrix
    type(sysvector_t), dimension(size(sol)) :: rhsd, sol_old
    type(coefficients_t) :: coefficients
    type(lu_ma57_t) :: lu

!   constants

    integer, parameter :: &
      uintpl = 6,         & ! scalar interpolation
      gauss = 6,          & ! 6-point Gauss integration of triangles
      gaussb = 3            ! 3 point integration of boundary elements

    integer :: i

    real(dp), parameter :: &
      alpha = 1._dp     ! diffusion coefficient. Not used in laplace_elem_ale
                        ! Only for compatibility with poisson_elem.


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

      call define_essential ( mesh, input_probdef, curve1=1, curve2=5 )

      call problem_definition ( input_probdef, mesh, problem )

    end if

!   create old solution vectors and right-hand side

    call create ( problem, sol_old, rhsd )

    if ( .not. sol(1)%created ) then
      call create ( problem, sol )
    else
!     copy mesh velocity to sol_old
      call copy ( sol, sol_old )
    end if


!   create system matrix

    call create_sysmatrix_structure ( sysmatrix, mesh, problem, &
      symmetric=.true. )
    call create_sysmatrix_data ( sysmatrix )

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, msysvector=rhsd, &
      elemsub=laplace_elem_ale, coefficients=coefficients )

!   solve

    do i = 1, size(sol)

!     fill solution vector with zero essential boundary conditions

      call fill_sysvector ( mesh, problem, sol(i), curve1=1, curve2=4, &
        value=0._dp )

!     fill solution vector with essential boundary conditions

      call fill_sysvector ( mesh, problem, sol(i), curve1=5, value=velp(i) )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol(i),&
         rhsd(i) )

!     solve the system
      call solve_system_ma57 ( sysmatrix, rhsd(i), sol(i), lu=lu )

    end do

!   update mesh

    do i = 1, size(sol)

      if ( method==1 ) then
!       update mesh nodes with forward Euler
        mesh%coor(:,i) = mesh%coor(:,i) + sol(i)%u(problem%degfdperm(:,2))*dt
      else if ( method==2 ) then
!       update mesh nodes with 2nd order Adams-Bashforth
        mesh%coor(:,i) = mesh%coor(:,i) + &
              dt*(3.0_dp*sol(i)%u(problem%degfdperm(:,2))/2.0_dp - &
              sol_old(i)%u(problem%degfdperm(:,2))/2.0_dp)
      end if

    end do

!   delete all data including all allocated memory

    call delete ( coefficients )
    if (input_probdef%created) call delete ( input_probdef )
    call delete ( rhsd )
    call delete ( sysmatrix )
    call delete ( lu )

  end subroutine update_mesh_nodes

end module update_mesh_nodes1_m
