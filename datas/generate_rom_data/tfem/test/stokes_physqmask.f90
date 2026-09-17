! Stokes problem on a 3D square cavity. Only half of the domain is solved
! due to symmetry conditions.

program stokes11a

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    nx = 4,            & ! number of elements in x-direction
    ny = 2,             & ! number of elements in y-direction
    nz = 4,            & ! number of elements in z-direction
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

  logical, parameter :: &
    physqmask = .true. ! no p-p zero block in the system matrix

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(meshgen_options_t) :: meshgen_options
  type(coefficients_t) :: coefficients

! variables

  integer :: presnod(8) = [1,3,9,7,19,21,27,25]


  !timer = .true.

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, &
    lz=lz, elshape=14, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(presnod,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2]
  input_probdef%physqmask(2,2) = .false.
!  input_probdef%physqmask = reshape ( &
!       [ .true., .true., &
!         .true., .false. ], [2,2], order=[2,1] )

  call define_essential ( mesh, input_probdef, surface1=1, surface2=3, physq=1 )
  call define_essential ( mesh, input_probdef, surface1=4, physq=1, &
    degfd=[0,1,0] )
  call define_essential ( mesh, input_probdef, surface1=5, surface2=6, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    surfaces=[1,2,3,5,6], physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=4, physq=1, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=6, physq=1, degfd=1, value=1._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )


! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, &
    symmetric=.true., usephysqmask=physqmask )

  call create_sysmatrix_data ( sysmatrix )

  print *, 'nnz = ', sysmatrix%Suu%nnz

  !call tic

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients, usephysqmask=physqmask )

  !call toc ( 'build' )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  print *, 'sum(abs(sol%u))/sol%n'
  print *, sum(abs(sol%u))/sol%n

  !call toc ( 'solve' )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program stokes11a
