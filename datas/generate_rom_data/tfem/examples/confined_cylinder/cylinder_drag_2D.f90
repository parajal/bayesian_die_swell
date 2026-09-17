! 2D viscoelastic problem: flow around a cylinder confined between two walls.

! Post-processing of the drag force

program cylinder_drag

  use tfem_m
  use viscoelastic_elements_m
  use io_utils_m

  implicit none


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc
  type(problem_t) :: problem, problemc
  type(sysvector_t), target :: sol
  type(vector_t), target :: pressure, stress_tensor, conformation_tensor(1)
  type(oldvectors_t) :: oldvectors_d, oldvectors_dve
  type(coefficients_t) :: coefficients

  type(sysvector_t), dimension(3,1), target :: solc

  integer :: i, ios

  real(dp) :: dragv(2), dragve(2)


! read coefficients

  call read_coefficients ( coefficients, filename="coefficients.out" )


! read mesh

  call read_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call read_input_probdef ( mesh, input_probdef, filename='probdef.out' )

  call problem_definition ( input_probdef, mesh, problem )


! problem definition conformation tensor

  call read_input_probdef ( mesh, input_probdefc, filename='probdefc.out' )

  call problem_definition ( input_probdefc, mesh, problemc )


! create system vectors for gradient/velocity/pressure (solution)

  call create_sysvector ( problem, sol )


! create system vectors (solution ) for conformation and

  call create_m2sysvector ( problemc, solc )


! read data for post-processing

  open ( unit = 10, file='data.out', form='unformatted', iostat=ios, &
         status='old' )

  if ( ios /= 0 ) then
    write(*,'(/2a/)') 'Error: cannot open file data.out '
    stop
  end if

  read(10) sol%u
  read(10) (solc(i,1)%u, i=1,3)

  close ( unit=10 )


! fill oldvectors for stokes problem

  call create_oldvectors ( oldvectors_d, nsysvec=1, nvec=2 )

  oldvectors_d%s(1)%p => sol

! pressure

  call create_vector ( problem, pressure, vec=4 )

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors_d )

! viscous stress

  call create_vector ( problem, stress_tensor, vec=5 )

  call derive_vector ( mesh, problem, stress_tensor, &
    elemsub=stokes_stress_tensor, coefficients=coefficients, &
    oldvectors=oldvectors_d, elcurves=[20] ) ! only include surface 20

  oldvectors_d%v(1)%p => stress_tensor
  oldvectors_d%v(2)%p => pressure

! integrate drag force

  call integrate_boundary_elements ( mesh, problem, dragv, &
    elemsub=stokes_drag, curve=20, coefficients=coefficients, &
    oldvectors=oldvectors_d )

  print *, 'viscous drag = ', -dragv(1)

  call delete ( pressure, stress_tensor )


! create oldvectors for viscoelastic problem

  call create_oldvectors ( oldvectors_dve, nsysvec2=1, nvec1=1 )

  oldvectors_dve%s2(1)%p => solc

! conformation tensor in the nodes

  call create ( problemc, conformation_tensor, vec=3 )

  call derive_vector ( mesh, problemc, conformation_tensor(1), &
    elemsub=deriv_conformation_tensor_std, coefficients=coefficients, &
    oldvectors=oldvectors_dve, elcurves=[20] ) ! only include surface 20

  oldvectors_dve%v1(1)%p => conformation_tensor

! integrate drag force

  call integrate_boundary_elements ( mesh, problemc, dragve, &
    elemsub=viscoelastic_drag, curve=20, coefficients=coefficients, &
    oldvectors=oldvectors_dve )

  print *, 'viscoelastic drag = ', -dragve(1)
  print *, 'total drag = ', -dragv(1)-dragve(1)

  call delete ( conformation_tensor )

! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol )
  call delete ( oldvectors_d, oldvectors_dve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( solc )

  call delete ( coefficients )

end program cylinder_drag
