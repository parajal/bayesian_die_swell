module write_vtks_m

  use tfem_m
  use viscoelastic_elements_m
  use io_utils_m

  implicit none

contains

! write the data of an implicit viscoelastic simulation to .vtk files

  subroutine write_vtks ( mesh, problem, sol, coefficients_in, filename )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(sysvector_t), intent(inout), target :: sol
    type(coefficients_t), intent(in) :: coefficients_in
    character(len=*), intent(in) :: filename

    type(oldvectors_t) :: oldvectors
    type(vector_t) :: velocity, pressure, vorticity, vonmises
    type(coefficients_t) :: coefficients
    type(vector_t) :: ctensor, btensor

    logical :: axisym
    integer :: nmodes

    coefficients = coefficients_in

    nmodes = coefficients%i(19)
    axisym = coefficients%i(23) == 1

!   post-processing

    call create_vector ( problem, velocity, physq=2 )
    call create_vector ( problem, pressure, vec=4+nmodes )
    call create_vector ( problem, vonmises, vec=4+nmodes )
    call create_vector ( problem, vorticity, vec=4+nmodes )
    call create_vector ( problem, ctensor, vec=5+nmodes )
    call create_vector ( problem, btensor, vec=6+nmodes )

    call extract_physvector ( mesh, problem, sol, velocity )

!   create the structure oldvectors

    call create_oldvectors ( oldvectors, nsysvec=1 )

!   derive vectors

    oldvectors%s(1)%p => sol

    call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
      coefficients=coefficients, oldvectors=oldvectors )

    coefficients%i(13)=5

    call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
      coefficients=coefficients, oldvectors=oldvectors )

    call write_vector_vtk ( mesh, problem, filename, &
      dataname='velocity', vector=velocity )
    call write_scalar_vtk ( mesh, problem, filename, &
      dataname='pressure', vector=pressure, append=.true. )
    call write_scalar_vtk ( mesh, problem, filename, &
      dataname='vorticity', vector=vorticity, append=.true. )

!   conformation tensor

    call derive_vector ( mesh, problem, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors )

    call write_tensor_vtk ( mesh, problem, filename=filename, &
      dataname='conformation_tensor', vector=ctensor, assume33=axisym, &
      append=.true. )

!   b tensor

    call derive_vector ( mesh, problem, btensor, &
      elemsub=deriv_conformation_tensor_std, &
      coefficients=coefficients, oldvectors=oldvectors )

    call write_tensor_vtk ( mesh, problem, filename, &
      dataname='b-tensor', vector=btensor, symmetric=.false., assume33=axisym, &
      append=.true. )

    coefficients%i(13) = 1 ! component
    coefficients%i(28) = 1 ! mode

    call derive_vector ( mesh, problem, vonmises, &
      elemsub=deriv_viscoelastic_stress_scalar, &
      coefficients=coefficients, oldvectors=oldvectors )

    call write_scalar_vtk ( mesh, problem, filename, &
      dataname='vonMises', vector=vonmises, append=.true. )

  end subroutine write_vtks

end module write_vtks_m
