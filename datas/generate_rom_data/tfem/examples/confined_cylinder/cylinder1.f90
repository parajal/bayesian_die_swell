! 2D Stokes problem: flow around a cylinder confined between two walls.
! The entry and exit sides are connected by periodical boundary conditions
! and a constant flow rate is imposed.
! Only the upper-half of the domain is modelled and symmetry conditions are
! assumed (zero tractions).

! The average velocity is U.
! Drag computations.
!
! TFEM generated mesh.

program cylinder1

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3             ! 3x3x3 integration of bricks

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity
  type(vector_t), target :: pressure, stress_tensor
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma57_t) :: solver_options

! variables

  logical :: tensoroutput = .true.

  real(dp) :: &
    ly = 2._dp,   & ! size in y-direction
    eta = 1._dp,  & ! viscosity
    U = 1._dp       ! average velocity in the channel

  real(dp) :: drag(2), flowrate

  integer :: presnod(4) = [1,3,5,7]


  flowrate = U * ly

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(6) = flowrate

  call write_coefficients ( coefficients, filename='coefficients.out' )

! read mesh

  call read_mesh ( mesh, filename='mesh.out' )

! connect curves for periodical bc in DG

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 25

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 2  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(presnod,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar
  input_probdef%vec_elementdof(1)%a(:,4) = 3  ! symmetric tensor

  input_probdef%physq = [1,2]

! define essential boundaries

! cylinder
  call define_essential ( mesh, input_probdef, curve1=20, physq=1 )
! top (wall)
  call define_essential ( mesh, input_probdef, curve1=21, physq=1 )
! bottom (center line)
  call define_essential ( mesh, input_probdef, curve1=23, curve2=24, physq=1, &
    degfd=[0,1] )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=25, nglobalc=1 )

! constraints for periodical boundary conditions

! velocities (use weak connection)
  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=18, curve2=25, discretization='weak', &
    elementdof=[2,0,2] )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1, nvec=2 )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_elem_conn, addmatvec=.true., &
    coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%real_storage=1.1

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

  call delete ( sysmatrix )

  call create_vector ( problem, velocity, physq=1 )
  call extract_physvector ( mesh, problem, sol, velocity )

! post processing

  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, stress_tensor, vec=4 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, stress_tensor, &
    elemsub=stokes_stress_tensor, coefficients=coefficients, &
    oldvectors=oldvectors )
!    oldvectors=oldvectors, elcurves=(/20/) ) ! only include curve 20

  oldvectors%v(1)%p => stress_tensor
  oldvectors%v(2)%p => pressure

! integrate drag force

  call integrate_boundary_elements ( mesh, problem, drag, &
    elemsub=stokes_drag, curve=20, coefficients=coefficients, &
    oldvectors=oldvectors )

  print *, 'drag force = ', -drag(1)

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='cylinder.vtk' )

  call write_vector_vtk ( mesh, problem, filename='cylinder.vtk', &
    dataname='velocity_vector', sysvector=sol, append=.true. )

  if ( tensoroutput ) then

    call write_tensor_vtk ( mesh, problem, vector=stress_tensor, &
      dataname='stress_tensor', filename='cylinder.vtk', append=.true. )

  else

    call write_scalar_vtk ( mesh, problem, vector=stress_tensor, degfd=1, &
      dataname='tau_xx', filename='cylinder.vtk', append=.true. )

    call write_scalar_vtk ( mesh, problem, vector=stress_tensor, degfd=2, &
      dataname='tau_xy', filename='cylinder.vtk', append=.true. )

    call write_scalar_vtk ( mesh, problem, vector=stress_tensor, degfd=3, &
      dataname='tau_xz', filename='cylinder.vtk', append=.true. )

    call write_scalar_vtk ( mesh, problem, vector=stress_tensor, degfd=4, &
      dataname='tau_yy', filename='cylinder.vtk', append=.true. )

    call write_scalar_vtk ( mesh, problem, vector=stress_tensor, degfd=5, &
      dataname='tau_yz', filename='cylinder.vtk', append=.true. )

    call write_scalar_vtk ( mesh, problem, vector=stress_tensor, degfd=6, &
      dataname='tau_zz', filename='cylinder.vtk', append=.true. )

  end if

! write binary file for reading by streamfunction

  open(unit=10,form='unformatted',file='velocity_bin.out')

  write(10) velocity%u

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, stress_tensor )
  call delete ( coefficients )
  call delete ( oldvectors )

end program cylinder1
