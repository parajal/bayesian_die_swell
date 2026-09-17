! Stokes problem on a rectangular (2D) domain with 3D velocities (axisymmetric).
! Cylindrical coordinates with circumferential velocity included (swirl).
! Periodical conditions and imposed flow rate in z-direction.

program concentric_cylinders1

  use tfem_m
  use math_defs_m
  use hsl_ma57_m
  use stokes_elements_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  logical, parameter :: concentric = .true.  ! F: Single pipe
                                             ! T: Two concentric walls

  integer, parameter :: &
    uintpl = 8,   & ! Q2 velocities
    pintpl = 4,   & ! Q1 pressures
    physqv = 1,   & ! physical quantity nr of the velocities
    physqp = 2,   & ! physical quantity nr of the pressures
    gauss = 3,    & ! 3x3 integration of quads
    nz=2,         & ! number of elements in z
    nr=20           ! number of elements in r

  real(dp), parameter :: &
    R1 = 0.5_dp,  & ! inner radius (for concentric=.true.)
    R2 = 1._dp,   & ! outer radius
    er = 10._dp,  & ! ratio of element length z/r.
    eta = 1._dp,  & ! viscosity
    U = 1._dp,    & ! imposed average velocity in z direction
    W = 1._dp       ! Velocity in theta direction of the (outer) wall.

  integer :: vertices(4) = [1,3,5,7]

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, duthetadr, gammadot, Dtensor, Ltensor
  type(vector_t) :: tau_rtheta, tautensor
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

  real(dp) :: flowrate

  if ( concentric ) then
    flowrate = U * pi * ( R2 ** 2 - R1 ** 2 )
  else
    flowrate = U * pi * R2 ** 2
  end if

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl, pintpl, 0, 0, 0,      &
      physqv, physqp, 0, 0, gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%i(23) = 1  ! cylindrical coordinate system (axisymmetrical)
  coefficients%i(67) = 1  ! 3D velocity (swirl)

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(6) = flowrate

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nz
  meshgen_options%ny = nr
  if ( concentric ) then
    meshgen_options%oy = R1
    meshgen_options%ly = R2 - R1
  else
    meshgen_options%ly = R2
  end if
  meshgen_options%lx = er * nz * meshgen_options%ly / nr

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(vertices,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,4) = 6  ! symmetric tensor
  input_probdef%vec_elementdof(1)%a(:,5) = 9  ! unsymmetric tensor

  input_probdef%physq = [1,2]

  if ( concentric ) then
!   inner cilinder wall
    call define_essential ( mesh, input_probdef, curve1=1, physq=physqv )
  else
!   center line
    call define_essential ( mesh, input_probdef, curve1=1, physq=physqv, &
      degfd=[0,1,1], exclude=1 )
  end if
! outer wall
  call define_essential ( mesh, input_probdef, curve1=3, physq=physqv )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=physqp )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqv, curve1=2, nglobalc=1 )

! constraint for periodical boundary conditions

  if ( concentric ) then
    call define_constraint ( mesh, input_probdef, physq=physqv, &
      curve1=2, curve2=5, discretization='collocation', exclude=3 )
  else
    call define_constraint ( mesh, input_probdef, physq=physqv, &
      curve1=2, curve2=5, discretization='collocation', exclude=2 )
  end if

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=physqv, degfd=3, value=W )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqp, value=0._dp )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr_curve, &
    addmatvec=.true., coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, duthetadr, vec=3 )
  call create_vector ( problem, gammadot, vec=3 )
  call create_vector ( problem, Dtensor, vec=4 )
  call create_vector ( problem, Ltensor, vec=5 )
  call create_vector ( problem, tau_rtheta, vec=3 )
  call create_vector ( problem, tautensor, vec=4 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=8

  call derive_vector ( mesh, problem, duthetadr, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=11

  call derive_vector ( mesh, problem, gammadot, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Dtensor, elemsub=stokes_D_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Ltensor, elemsub=stokes_gradu_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, tau_rtheta, elemsub=stokes_stress, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, tautensor, elemsub=stokes_stress_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='concyl1.vtk' )

  call write_vector_vtk ( mesh, problem, filename='concyl1.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true., &
    assume3D=.true. )

  call write_scalar_vtk ( mesh, problem, vector=duthetadr, &
    filename='concyl1.vtk', dataname='duthetadr', append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='concyl1.vtk', dataname='gammadot', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='concyl1.vtk', &
    dataname='D', vector=Dtensor, append=.true., assume3D=.true. )

  call write_tensor_vtk ( mesh, problem, filename='concyl1.vtk', &
    dataname='L', vector=Ltensor, append=.true., assume3D=.true., &
    symmetric=.false. )

  call write_scalar_vtk ( mesh, problem, vector=tau_rtheta, &
    filename='concyl1.vtk', dataname='tau_rtheta', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='concyl1.vtk', &
    dataname='tau', vector=tautensor, append=.true., assume3D=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, duthetadr, gammadot, Dtensor, Ltensor )
  call delete ( tau_rtheta, tautensor )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program concentric_cylinders1
