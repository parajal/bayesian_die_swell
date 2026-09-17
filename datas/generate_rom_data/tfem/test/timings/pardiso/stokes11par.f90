! Stokes problem on a 3D square cavity. Only half of the domain is solved
! due to symmetry conditions.

program stokes11

  use tfem_m
  use pardiso_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
!    nx = 50,            & ! number of elements in x-direction
!    ny = 25,            & ! number of elements in y-direction
!    nz = 50,            & ! number of elements in z-direction
    nx = 5,            & ! number of elements in x-direction
    ny = 5,            & ! number of elements in y-direction
    nz = 5,            & ! number of elements in z-direction
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
  type(plot_options_t) :: plot_options
  type(meshgen_options_t) :: meshgen_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_pardiso_t) :: so

! variables

  integer :: presnod(8) = [1,3,9,7,19,21,27,25]


  call tic

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, &
    lz=lz, elshape=14, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! plot curves, surfaces and mesh

  plot_options%viewpoint=[1.,0.8,0.4]
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig', surfaces=[3,4,6] )

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

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call toc ( 'before_solve' )

  so%printlevel=3
  call solve_system_pardiso ( sysmatrix, rhsd, sol, solver_options=so )

  call toc ( 'after_solve' )

! post-processing

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

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='cavity.vtk' )

  call write_vector_vtk ( mesh, problem, filename='cavity.vtk', &
    dataname='velocity_vector', sysvector=sol, append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=dudz, filename='cavity.vtk', &
    dataname='dudz', append=.true. )


! write to a fig file for plotting

!  plot_options%printlabels=.false.

!  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
!    vector=velocity, surfaces=(/4,6/) )
!  call plot_points_curves ( plot_options, mesh, 'velocity.fig', append=.true. )
!
!  call plot_color_fill ( plot_options, mesh, problem, 'velocity_color.fig', &
!    vector=velocity, degfd=1, surfaces=(/4/) )
!  call plot_points_curves ( plot_options, mesh, 'velocity_color.fig',  &
!   append=.true. )
!
!  call plot_color_contour ( plot_options, mesh, problem, 'dudz_contour.fig', &
!    vector=dudz, surfaces=(/3,4,6/) )
!
!  call plot_color_contour ( plot_options, mesh, problem, &
!    'pressure_contour.fig', vector=pressure, surfaces=(/3,4,6/) )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( velocity, pressure, dudz )
  call delete ( coefficients )
  call delete ( oldvectors )

  call toc ( 'end' )

end program stokes11
