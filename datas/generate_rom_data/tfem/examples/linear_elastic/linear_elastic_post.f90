! Linear elastic post processing

module linear_elastic_post_m

  use tfem_m
  use linear_elastic_elements_m
  use io_utils_m

  implicit none

contains


! Routine for post processing linear elastic problems

  subroutine postprocessing_linear_elastic ( mesh, problem, vtkfilename, sol, &
    vec_scalar, vec_tensor, coefficients, mcoefficients, assume3D, assume33 )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the vtk postprocessing output filename for writing the data to
    character (len=*), intent(in) :: vtkfilename

!   sysvector containing the solution
    type(sysvector_t), target, intent(in) :: sol

!   vector number for the scalar and tensor output
    integer, intent(in) :: vec_scalar, vec_tensor

!   if coefficients is present, it is passed through to the element subroutine
    type(coefficients_t), intent(inout), optional :: coefficients

!   array of coefficients (length=number of element groups)
!   If present only one of them is passed through to the element subroutine
!   based on the element group number
    type(coefficients_t), dimension(:), intent(inout), optional :: mcoefficients

!   Assume a 3D tensor with the 33 diagonal component present
!   (only applicable for mesh%ndim = 2).
!   default=.false.
    logical, intent(in), optional :: assume33

!   Assume a 3D tensor (even if mesh%ndim /= 3).
!   default=.false.
    logical, intent(in), optional :: assume3D


    type(vector_t) :: scalar, tensor
    type(oldvectors_t) :: oldvectors


!   post-processing

    call create_vector ( problem, scalar, vec=vec_scalar )
    call create_vector ( problem, tensor, vec=vec_tensor )

  ! create the structure oldvectors

    call create_oldvectors ( oldvectors, nsysvec=1 )

    oldvectors%s(1)%p => sol

  ! write to a vtk file for plotting

    call write_vector_vtk ( mesh, problem, vtkfilename, 'displacement', &
      physq=1, sysvector=sol, assume3D=assume3D )

    if ( coefficients%i(67) == 1 ) then

!     3D vector on a 2D mesh

      call write_vector_vtk ( mesh, problem, vtkfilename, &
        'displacement_planar', degfd=[1,2], &
        physq=1, sysvector=sol, assume3D=assume3D, append=.true. )

      call write_vector_vtk ( mesh, problem, vtkfilename, &
        'displacement_out_of_plane', degfd=[0,0,3], &
        physq=1, sysvector=sol, assume3D=assume3D, append=.true. )

    end if

    coefficients%i(13)=1

    call derive_vector ( mesh, problem, scalar, &
      elemsub=linear_elastic_deriv, &
      coefficients=coefficients, oldvectors=oldvectors )

    call write_scalar_vtk ( mesh, problem, vtkfilename, 'divergence', &
      vector=scalar, append=.true. )

    coefficients%i(13)=1

    call derive_vector ( mesh, problem, scalar, &
      elemsub=linear_elastic_stress_deriv, &
      coefficients=coefficients, oldvectors=oldvectors )

    call write_scalar_vtk ( mesh, problem, vtkfilename, 'hydrostatic_stress', &
      vector=scalar, append=.true. )

    coefficients%i(13)=2

    call derive_vector ( mesh, problem, scalar, &
      elemsub=linear_elastic_stress_deriv, &
      coefficients=coefficients, oldvectors=oldvectors )

    call write_scalar_vtk ( mesh, problem, vtkfilename, 'von_Mises_stress', &
      vector=scalar, append=.true. )

    call derive_vector ( mesh, problem, tensor, &
      elemsub=linear_elastic_strain_tensor, &
      coefficients=coefficients, oldvectors=oldvectors )

    call write_tensor_vtk ( mesh, problem, vtkfilename, 'strain_tensor', &
      vector=tensor, assume3D=assume3D, assume33=assume33, append=.true. )

    call derive_vector ( mesh, problem, tensor, &
      elemsub=linear_elastic_stress_tensor, &
      coefficients=coefficients, oldvectors=oldvectors )

    call write_tensor_vtk ( mesh, problem, vtkfilename, 'stress_tensor', &
      vector=tensor, assume3D=assume3D, assume33=assume33, append=.true. )


  ! delete all data including all allocated memory

    call delete ( scalar, tensor )
    call delete ( oldvectors )

  end subroutine postprocessing_linear_elastic

end module linear_elastic_post_m
