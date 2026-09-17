! Stokes problem on a 3D square cavity. Only half of the domain is solved
! due to symmetry conditions.

program stokes11_pardiso_sym

  use tfem_m
  use pardiso_m
  use stokes_elements_m
  use io_utils_m
  use timer_m

  implicit none

! constants

  logical, parameter :: post = .false.

  integer, parameter :: &
! 2 GBytes (MA57)
!    nx = 20,            & ! number of elements in x-direction
!    ny = 10,            & ! number of elements in y-direction
!    nz = 20,            & ! number of elements in z-direction
! 9 GBytes (MA57)
!    nx = 30,            & ! number of elements in x-direction
!    ny = 15,            & ! number of elements in y-direction
!    nz = 30,            & ! number of elements in z-direction
! 2.5 GBytes (MA41)
!    nx = 18,            & ! number of elements in x-direction
!    ny = 9,             & ! number of elements in y-direction
!    nz = 18,            & ! number of elements in z-direction
! 8.5 GBytes (MA41)
    nx = 26,            & ! number of elements in x-direction
    ny = 13,            & ! number of elements in y-direction
    nz = 26,            & ! number of elements in z-direction
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3             ! 3x3x3 integration of hexahedra

  real(dp), parameter :: &
    lx = 1._dp,          & ! size in x-direction
    ly = 0.5_dp,         & ! size in y-direction
    lz = 1._dp,          & ! size in z-direction
    eta = 1._dp            ! viscosity

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, dudz
  type(meshgen_options_t) :: meshgen_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_pardiso_t) :: solver_options_pardiso

! variables

  logical :: impose_iterative_refinement = .false.
  integer :: presnod(8) = [1,3,9,7,19,21,27,25]


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

  print '(3(a,i0)/)', 'Mesh parameters: nx=', nx, ' ny=', ny, ' nz=', nz

  call tic

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, &
    lz=lz, elshape=14, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

  call toc ( 'mesh' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(presnod,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, surface1=1, surface2=3, physq=1 )
  call define_essential ( mesh, input_probdef, surface1=4, physq=1, &
    degfd=[0,1,0] )
  call define_essential ( mesh, input_probdef, surface1=5, surface2=6, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

  call problem_definition ( input_probdef, mesh, problem )

  call toc ( 'problem definition' )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    surface1=1, surface2=3, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=4, physq=1, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=5, surface2=6, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=6, physq=1, degfd=1, value=1._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  call toc ( 'vectors' )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

  call toc ( 'create_sysmatrix' )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call check ( sysmatrix )

  call toc ( 'build_system' )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options_pardiso%matrixtype=-2
  if ( impose_iterative_refinement ) then
    solver_options_pardiso%niter_ref=10
    solver_options_pardiso%printlevel=4
  end if
  call solve_system_pardiso ( sysmatrix, rhsd, sol, &
    solver_options=solver_options_pardiso  )

  call toc ( 'solve' )

  if ( post ) then

!   post-processing

    call create_vector ( problem, velocity, physq=1 )
    call create_vector ( problem, pressure, vec=3 )
    call create_vector ( problem, dudz, vec=3 )

    call extract_physvector ( mesh, problem, sol, velocity )

    oldvectors%s(1)%p => sol

    call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
      coefficients=coefficients, oldvectors=oldvectors )

    coefficients%i(13)=3
    call derive_vector ( mesh, problem, dudz, elemsub=stokes_deriv, &
      coefficients=coefficients, oldvectors=oldvectors )

!   write to a vtk file for further post-processing

    call write_scalar_vtk ( mesh, problem, vector=pressure, &
      dataname='pressure', filename='cavity11_pardiso_sym.vtk' )

    call write_vector_vtk ( mesh, problem, filename='cavity11_pardiso_sym.vtk',&
      dataname='velocity_vector', sysvector=sol, append=.true. )

    call write_scalar_vtk ( mesh, problem, vector=dudz, &
      filename='cavity11_pardiso_sym.vtk', dataname='dudz', append=.true. )

    call toc ( 'post-processing' )

    call delete ( velocity, pressure, dudz )

  end if

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program stokes11_pardiso_sym
