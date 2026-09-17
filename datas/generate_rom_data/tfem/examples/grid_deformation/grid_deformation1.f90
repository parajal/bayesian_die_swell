! Grid deformation: elements are "attracted" towards a circle.

program grid_deformation1

  use kind_defs_m
  use tfem_m
  use meshgen_grid_deformation_m
  use deform_mesh_ma57_m
  use subs_grid_deformation1_m
  use figplot_m
  use timer_m

  implicit none

  integer, parameter :: &
    nx=25,              & ! number of elements in x
    ny=25,              & ! number of elements in y
    uintpl = 8,         & ! scalar interpolation for Poisson problem
    gauss = 3,          & ! 3x3 integration of quads
    max_gridsteps = 20    ! number of pseudo time steps in the grid deformation

  real(dp), parameter :: &
!   The (relative) size of the elements near the zero levelset (min_f) and far
!   away from the the zero levelset (max_f).
    min_f = 0.02,  &
    max_f = 0.2,  &
    ox = -0.5_dp, &  ! mesh domain parameters
    oy = -0.5_dp, &
    lx = 1._dp,   &
    ly = 1._dp,   &
    radius=0.2,   &  ! radius of the circle
    center(2)=0._dp  ! center of the circle

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(problem_t) :: problem_gd
  type(coefficients_t) :: coefficients_gd
  type(sysmatrix_t) :: sysmatrix_gd
  type(plot_options_t) :: plot_options

  integer :: nnodes, obj_ob, obj_gd
  real(dp), allocatable, dimension(:) :: f_monitor
  real(dp), allocatable, dimension(:,:) :: coor


  timer = .false.

! set values in subs module

  lox = ox; loy = oy; llx = lx; lly = ly
  lmin_f = min_f; lmax_f = max_f

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, elshape=6, &
    lx=lx, ly=ly, ox=ox, oy=oy )

  call quadrilateral2d ( mesh, meshgen_options )

  call tic

! use sblocks to speed up the search for interpolation points

  call add_to_mesh ( mesh, sblocks=[4*nx,4*ny] )
  call toc('sblocks')

! radius and center position of an object

  rpl = radius
  cpl = center

! create object for the circle boundary (only for plotting)

  nnodes = 50

  allocate ( coor(nnodes,2) )

  call objectscoor ( 1, coor )

  call add_to_mesh ( mesh, object='coordinates', coor=coor )
  obj_ob = mesh%nobjects

  deallocate ( coor )

  call fill_mesh_parts ( mesh )

  call plot_mesh ( plot_options, mesh, 'mesh_initial.fig' )
  plot_options%objectpointcolor=4
  plot_options%objectpointsize=0.4
  call plot_objects ( plot_options, mesh, 'mesh_initial.fig', append=.true., &
    object1=obj_ob )


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
    project_boundary=project_boundary, plotfigs=.true. )


! update coordinates of the mesh

  mesh%coor = mesh%objects(obj_gd)%coor


! plot deformed mesh

  call plot_mesh ( plot_options, mesh, 'mesh_deformed.fig' )
  plot_options%objectpointcolor=4
  plot_options%objectpointsize=0.4
  call plot_objects ( plot_options, mesh, 'mesh_deformed.fig', append=.true., &
    object1=obj_ob )


! delete all data including all allocated memory

  call delete ( problem_gd )
  call delete ( mesh )
  call delete ( sysmatrix_gd )
  call delete ( coefficients_gd )

end program grid_deformation1

