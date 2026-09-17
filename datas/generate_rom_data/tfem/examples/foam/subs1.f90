module subs1_m

  use tfem_elem_m

  implicit none

contains


! element for force free motion in x-direction

  subroutine elementc1 ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

!   connection through collocation

    elemmat(1,:) = 0
    elemmat(1,1) = 1

    elemmatadd(1,1) = -1._dp

  end subroutine elementc1


! element for force free motion in y-direction

  subroutine elementc2 ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

!   connection through collocation

    elemmat(1,:) = 0
    elemmat(1,2) = 1

    elemmatadd(1,1) = -1._dp

  end subroutine elementc2


! element for force free motion in z-direction

  subroutine elementc3 ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

!   connection through collocation

    elemmat(1,:) = 0
    elemmat(1,3) = 1

    elemmatadd(1,1) = -1._dp

  end subroutine elementc3


! routine for output of the field data

  subroutine output_fields ( mesh, problem, coefficients, sol, post )

    use stokes_elements_m
    use io_utils_m

!   input/output

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the coefficients used for computing the derivative quantities
    type(coefficients_t), intent(in) :: coefficients

!   solution vector containing velocity and pressure
    type(sysvector_t), target, intent(in) :: sol

!   the post processing number for the output
    integer, intent(in) :: post


!   local defs

    type(vector_t) :: pressure, viscous_stress
    type(oldvectors_t) :: oldvectors

    character(len=30) :: filename


!   Creation of oldvectors structure for velocity/pressure problem

    call create_oldvectors ( oldvectors, nsysvec=1 )

!   Storage of solution vector

    oldvectors%s(1)%p => sol

!   compute pressure field

    call create_vector ( problem, pressure, vec=3 )

    call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
      coefficients=coefficients, oldvectors=oldvectors )

!   compute viscous stress field

    call create_vector ( problem, viscous_stress, vec=4 )

    call derive_vector ( mesh, problem, viscous_stress, &
      elemsub=stokes_stress_tensor, coefficients=coefficients, &
      oldvectors=oldvectors )

!   output fields

    write(filename,'(a,i4.4,a)') 'fields', post, '.vtk'

    call write_vector_vtk ( mesh, problem, filename=filename, &
      dataname='velocity', physq=1, sysvector=sol )
    call write_scalar_vtk ( mesh, problem, filename=filename, &
      dataname='pressure', vector=pressure, append=.true. )
    call write_tensor_vtk ( mesh, problem, filename=filename, &
      dataname='viscous_stress', vector=viscous_stress, append=.true. )

!   remove data

    call delete ( oldvectors )
    call delete ( pressure, viscous_stress )

  end subroutine output_fields

end module subs1_m
