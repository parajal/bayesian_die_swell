module subs_m

  use kind_defs_m
  use tfem_elem_m

  implicit none

  real(dp) :: lly, lshear_rate

contains

! Element for the constraints (connection through collocation)

  subroutine stokes_constr_node_conn_shear ( mesh, problem, constr, elem, &
    node, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


    integer :: i

!   connection through collocation

    elemmat  = 0
    do i = 1, size(elemmat,1)
      elemmat(i,i) = 1
    end do
    elemmat2 = - elemmat
    elemvec  = 0
    elemvec(1) = -lshear_rate*lly

  end subroutine stokes_constr_node_conn_shear

end module subs_m
