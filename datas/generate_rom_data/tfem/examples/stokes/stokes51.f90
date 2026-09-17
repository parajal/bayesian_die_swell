! Stokes problem on a 3D square cavity. Only half of the domain is solved
! due to symmetry conditions.
! Prisms.

program stokes51

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use io_utils_m
  use limits_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 16,        & ! P2Q2 velocities
    pintpl = 15,        & ! P1Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 6,          & ! 6 point integration of triangular base
    gauss2 = 3            ! 3 point integration in the height

  real(dp), parameter :: &
    eta = 1._dp     ! viscosity

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: pressure, vorticity, dudy
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

  integer :: vertices(6) = [1,3,5,13,15,17]

  WARN_ON_TEST_FOR_VARYING_SIDES = .false.

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       0 ]
  coefficients%i(12:) = 0
  coefficients%i(75) = gauss2

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

! read mesh

  call read_mesh_gmsh ( mesh, filename='cube51.msh' )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(vertices,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, surface1=1, physq=1 )
  call define_essential ( mesh, input_probdef, surface1=2, physq=1, &
    degfd=[0,0,1] )
  call define_essential ( mesh, input_probdef, surface1=3, surface2=6, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    surface1=1, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=2, physq=1, degfd=3, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=3, surface2=6, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=6, physq=1, degfd=1, value=-1._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! post-processing

  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )
  call create_vector ( problem, dudy, vec=3 )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=2
  call derive_vector ( mesh, problem, dudy, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='stokes51.vtk' )

  call write_vector_vtk ( mesh, problem, filename='stokes51.vtk', &
    dataname='velocity_vector', sysvector=sol, append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=dudy, filename='stokes51.vtk', &
    dataname='dudy', append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=vorticity, &
    filename='stokes51.vtk', dataname='vorticity', append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( pressure, vorticity, dudy )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program stokes51
