! Grid deformation for where the elements are "attracted" towards a sphere.
! Periodical boundary conditions in x-direction

program grid_deformation4

  use kind_defs_m
  use tfem_m
  use meshgen_grid_deformation_m
  use deform_mesh_ma57_m
  use subs_grid_deformation4_m
  use timer_m
  use io_utils_m


  implicit none


  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3 integration of quads
    nx=20,              & ! number of elements in x
    ny=20,              & ! number of elements in y
    nz=20,              & ! number of elements in z
    max_gridsteps = 20    ! number of pseudo time steps in the grid deformation

  real(dp), parameter :: &
!   The (relative) size of the elements near the zero levelset (min_f) and far
!   away from the the zero levelset (max_f).
    min_f = 0.02,  &
    max_f = 0.2,  &
    ox = -0.5_dp, &  ! mesh domain parameters
    oy = -0.5_dp, &
    oz = -0.5_dp, &
    lx = 1._dp,   &
    ly = 1._dp,   &
    lz = 1._dp,   &
    radius=0.2,   &  ! radius of the sphere
    center(3)=[0.15_dp,0._dp,0._dp]  ! center of the sphere


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(problem_t) :: problem_gd
  type(sysmatrix_t) :: sysmatrix_gd
  type(coefficients_t) :: coefficients_gd

  integer :: obj_gd
  real(dp), allocatable, dimension(:) :: f_monitor


  timer = .false.

! set values in subs module

  lox = ox; loy = oy; loz = oz; llx = lx; lly = ly; llz = lz

  lmin_f = min_f
  lmax_f = max_f

! create mesh

  call set_mesh_options ( meshgen_options, regionshape=4, nx=nx, ny=ny, nz=nz, &
    elshape=14, lx=lx, ly=ly, lz=lz, ox=ox, oy=oy, oz=oz )

  call hexahedron ( mesh, meshgen_options )

  call tic

! use sblocks to speed up the search for interpolation points

  call add_to_mesh ( mesh, sblocks=[4*nx,4*ny,4*nz] )
  call toc('sblocks')

! object for computing the deformed mesh

  call add_to_mesh ( mesh, object='coordinates', coor=mesh%coor )
  obj_gd = mesh%nobjects

! radius and center position of an object

  rpl = radius
  cpl = center

  call fill_mesh_parts ( mesh )

  call write_mesh_vtk ( mesh, filename='mesh_initial.vtk' )


! fill coefficients

  call create_coefficients ( coefficients_gd, ncoefi=100, ncoefr=50 )

  coefficients_gd%i = 0
  coefficients_gd%i(1) = uintpl
  coefficients_gd%i(10) = gauss

  coefficients_gd%r = 0


! set monitor function f

  allocate ( f_monitor(mesh%nnodes) )

  call set_monitor_function ( mesh, f_monitor )


! perform the grid deformation by pseudo time stepping

  obj_gd = 0  ! first call: set to non-existing object

  call deform_mesh_ma57 ( mesh, obj_gd, problem_gd, coefficients_gd, &
    sysmatrix_gd, f_monitor, numgridsteps=max_gridsteps, &
    constraint_definition=constraint_definition, &
    move_object_points=move_object_points, project_boundary=project_boundary )


! update coordinates of the mesh

  mesh%coor = mesh%objects(obj_gd)%coor


! plot deformed mesh

  call write_mesh_vtk ( mesh, filename='mesh_deformed.vtk' )


! delete all data including all allocated memory

  call delete ( problem_gd )
  call delete ( mesh )
  call delete ( sysmatrix_gd )
  call delete ( coefficients_gd )

end program grid_deformation4

