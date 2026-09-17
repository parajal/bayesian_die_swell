! Output fields for drop3
! The external and internal part are written to separate files to avoid
! smoothing the jumps at the interface by paraview.
! Note, that separate meshes and problems are created for the external
! and internal part and only the data corresponding to the internal and external
! part is written, respectively. This is different from using groups on the
! original problem, in which the full mesh and all the nodal data are written
! (twice: both for the internal and external datafile). In that case, the
! separation internal/external is only made by outputing only the elements
! (connectivity) from the specified groups.

module output_fields3_m

  use tfem_m
  use viscoelastic_elements_m
  use io_utils_m

  implicit none

contains


! routine for output of the field data

  subroutine output_fields ( mesh, problem, problem_c, coefficients, &
    sol, sol_c, post )

!   input/output

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem, problem_c

!   the coefficients used for computing the derivative quantities
    type(coefficients_t), intent(in), dimension(2) :: coefficients

!   solution vector containing gradient, velocity and pressure
    type(sysvector_t), target, intent(in) :: sol

!   solution vector containing conformation
    type(sysvector_t), dimension(:,:), target, intent(in) :: sol_c

!   the post processing number for the output
    integer, intent(in) :: post


!   local defs

    type(mesh_t) :: mesh_ext, mesh_int
    type(input_probdef_t) :: input_probdef_ext, input_probdef_int
    type(problem_t) :: problem_ext, problem_int
    type(vector_t) :: pressure, viscous_stress
    type(vector_t) :: conformation, viscoelastic_stress
    type(vector_t) :: velocity_ext, pressure_ext, velocity_int, pressure_int
    type(vector_t) :: viscous_stress_ext, viscous_stress_int
    type(vector_t) :: viscoelastic_stress_int, conformation_int
    type(subscript_t) :: ssvel_ext, ssvel_int
    type(subscriptvec_t) :: ssprs_ext, ssprs_int, ssten_ext, ssten_int, &
                            sstenc_int
    type(oldvectors_t) :: oldvectors, oldvectors_c

    integer :: physqvel, ndim
    character(len=30) :: filename


!   Creation of oldvectors structure for velocity/pressure problem

    call create_oldvectors ( oldvectors, nsysvec=1 )

!   Storage of solution vector

    oldvectors%s(1)%p => sol

!   compute pressure field

    call create_vector ( problem, pressure, vec=4 )

    call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
      mcoefficients=coefficients, oldvectors=oldvectors )

!   compute viscous stress field

    call create_vector ( problem, viscous_stress, vec=5 )

    call derive_vector ( mesh, problem, viscous_stress, &
      elemsub=stokes_stress_tensor, mcoefficients=coefficients, &
      oldvectors=oldvectors )

!   Creation of oldvectors structure for conformation problem

    call create_oldvectors ( oldvectors_c, nsysvec2=1 )

!   Storage of solution vector

    oldvectors_c%s2(1)%p => sol_c

!   compute viscoelastic stress field inside the droplet (element group group 2)

    call create_vector ( problem_c, viscoelastic_stress, vec=3 )

    call derive_vector ( mesh, problem_c, viscoelastic_stress, &
      elemsub=deriv_viscoelastic_stress_tensor, &
      coefficients=coefficients(2), oldvectors=oldvectors_c, groups=[2] )

!   compute conformation field inside the droplet (element group group 2)

    call create_vector ( problem_c, conformation, vec=3 )

    call derive_vector ( mesh, problem_c, conformation, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients(2), oldvectors=oldvectors_c, groups=[2] )


!   subscripts for extracting velocities, pressures and tensors

    physqvel = coefficients(1)%i(6)

    call create_subscript ( mesh, problem, ssvel_ext, groups=[1], &
      physqarr=[physqvel] )
    call create_subscript ( mesh, problem, ssprs_ext, groups=[1], vec=4 )
    call create_subscript ( mesh, problem, ssten_ext, groups=[1], vec=5 )

    call create_subscript ( mesh, problem, ssvel_int, groups=[2], &
      physqarr=[physqvel] )
    call create_subscript ( mesh, problem, ssprs_int, groups=[2], vec=4 )
    call create_subscript ( mesh, problem, ssten_int, groups=[2], vec=5 )
    call create_subscript ( mesh, problem_c, sstenc_int, groups=[2], vec=3 )

!   split mesh in external and internal mesh

    call mesh_convert ( mesh, mesh_ext, only_groups=[1], warn=.false. )
    call mesh_convert ( mesh, mesh_int, only_groups=[2], warn=.false. )

    call fill_mesh_parts ( mesh_ext )
    call fill_mesh_parts ( mesh_int )

    ndim = mesh%ndim

!   problem definition in external fluid for plotting

    call create_input_probdef ( mesh_ext, input_probdef_ext, nvec=3 )

    input_probdef_ext%elementdof(1)%a = 0  ! no need for sol
    input_probdef_ext%vec_elementdof(1)%a(:,1) = ndim  ! vector quantity
    input_probdef_ext%vec_elementdof(1)%a(:,2) = 1     ! scalar quantity
    input_probdef_ext%vec_elementdof(1)%a(:,3) = ndim*(ndim+1)/2  ! sym. tensor

    call problem_definition ( input_probdef_ext, mesh_ext, problem_ext )

    call create_vector ( problem_ext, velocity_ext, vec=1 )
    call create_vector ( problem_ext, pressure_ext, vec=2 )
    call create_vector ( problem_ext, viscous_stress_ext, vec=3 )

!   problem definition in internal fluid for plotting

    call create_input_probdef ( mesh_int, input_probdef_int, nvec=3 )

    input_probdef_int%elementdof(1)%a = 0  ! no need for sol
    input_probdef_int%vec_elementdof(1)%a(:,1) = ndim  ! vector quantity
    input_probdef_int%vec_elementdof(1)%a(:,2) = 1     ! scalar quantity
    input_probdef_int%vec_elementdof(1)%a(:,3) = ndim*(ndim+1)/2  ! sym. tensor

    call problem_definition ( input_probdef_int, mesh_int, problem_int )

    call create_vector ( problem_int, velocity_int, vec=1 )
    call create_vector ( problem_int, pressure_int, vec=2 )
    call create_vector ( problem_int, viscous_stress_int, vec=3 )
    call create_vector ( problem_int, viscoelastic_stress_int, vec=3 )
    call create_vector ( problem_int, conformation_int, vec=3 )

!   extract data

    velocity_ext%u = sol%u(ssvel_ext%s)
    pressure_ext%u = pressure%u(ssprs_ext%s)
    viscous_stress_ext%u = viscous_stress%u(ssten_ext%s)
    velocity_int%u = sol%u(ssvel_int%s)
    pressure_int%u = pressure%u(ssprs_int%s)
    viscous_stress_int%u = viscous_stress%u(ssten_int%s)
    viscoelastic_stress_int%u = viscoelastic_stress%u(sstenc_int%s)
    conformation_int%u = conformation%u(sstenc_int%s)

!   external fields

    write(filename,'(a,i4.4,a)') 'fields_ext', post, '.vtk'

    call write_vector_vtk ( mesh_ext, problem_ext, filename=filename, &
      dataname='velocity', vector=velocity_ext )
    call write_scalar_vtk ( mesh_ext, problem_ext, filename=filename, &
      dataname='pressure', vector=pressure_ext, append=.true. )
    call write_tensor_vtk ( mesh_ext, problem_ext, filename=filename, &
      dataname='viscous_stress', vector=viscous_stress_ext, append=.true. )

!   internal fields

    write(filename,'(a,i4.4,a)') 'fields_int', post, '.vtk'

    call write_vector_vtk ( mesh_int, problem_int, filename=filename, &
      dataname='velocity', vector=velocity_int )
    call write_scalar_vtk ( mesh_int, problem_int, filename=filename, &
      dataname='pressure', vector=pressure_int, append=.true. )
    call write_tensor_vtk ( mesh_int, problem_int, filename=filename, &
      dataname='viscous_stress', vector=viscous_stress_int, append=.true. )
    call write_tensor_vtk ( mesh_int, problem_int, filename=filename, &
      dataname='viscoelastic_stress', vector=viscoelastic_stress_int, &
      append=.true. )
    call write_tensor_vtk ( mesh_int, problem_int, filename=filename, &
      dataname='conformation', vector=conformation_int, append=.true. )

!   remove data

    call delete ( oldvectors, oldvectors_c )
    call delete ( pressure, viscous_stress )
    call delete ( conformation, viscoelastic_stress )
    call delete ( mesh_ext, mesh_int )
    call delete ( problem_ext, problem_int )
    call delete ( velocity_ext, pressure_ext, velocity_int, pressure_int )
    call delete ( viscous_stress_ext, viscous_stress_int )
    call delete ( viscoelastic_stress_int, conformation_int )
    call delete ( ssvel_ext, ssvel_int )
    call delete ( ssprs_ext, ssprs_int )
    call delete ( ssten_ext, ssten_int )
    call delete ( sstenc_int )

  end subroutine output_fields

end module output_fields3_m
