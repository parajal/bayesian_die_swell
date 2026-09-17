!
!  mesh build from 4 submeshes:
!
!
!     R5  -----------------------
!         |                     |
!         |                     |
!     n5  |         4           |
!         |                     |
!         |                     |
!     R4  -----------------------
!         |                     |
!         |                     |
!     n4  |         3           |
!         |                     |
!         |                     |
!     R3  -----------------------
!         |                     |
!         |                     |
!     n3  |         2           |
!         |                     |
!         |                     |
!     R2  -----------------------
!         |                     |
!         |                     |
!     n2  |         1           |
!         |                     |
!         |                     |
!     R1  -----------------------
!         0                     H
!                  n1
!

! NOTE: one single group is created. Use nogroupmerge=.true. in heading of
! mesh_merge to obtain 4 separate element groups.

program mesh_four_blocks

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use figplot_m
  use io_utils_m

  implicit none

  type(meshgen_options_t) :: mesh_options
  type(plot_options_t) :: plot_options
  type(mesh_t) :: mesh, mesh1, mesh2, mesh3

  integer :: elshape=6
  real(dp) :: R1=1, R2=3, R3=5, R4=7, R5=9
  real(dp) :: H=5
  integer :: n1=12, n2=6, n3=6, n4=6, n5=6

  real(dp) :: f1=3, f2=3, f3=3, f4=3, f5=3


! namelist for input of variables; read from standard input

  namelist /meshpar/ elshape, R1, R2, R3, R4, R5, H, n1, n2, n3, n4, n5, &
    f1, f2, f3, f4, f5

  read ( unit=*, nml=meshpar )


  plot_options%fontsize=8

! 1

 call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=n2, oy=R1, &
    lx=H, ly=R2-R1, ratio=[3,1,1,3], factor=[f1,f2,f1,f2] )

  call quadrilateral2d ( mesh1, mesh_options )
  !call plot_points_curves ( plot_options, mesh1, 'curves1.fig' )

! 2

  call set_mesh_options ( mesh_options, elshape=elshape, oy=R2, nx=n1, ny=n3, &
    lx=H, ly=R3-R2, ratio=[3,3,1,1], factor=[f1,f3,f1,f3] )

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, curve1=3, curve2=-1 )

  !call plot_points_curves ( plot_options, mesh, 'curves2.fig' )
  call delete ( mesh1, mesh2 )

! 3

  call set_mesh_options ( mesh_options, elshape=elshape, oy=R3, nx=n1, ny=n4, &
    lx=H, ly=R4-R3, ratio=[3,1,1,3], factor=[f1,f4,f1,f4] )

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh3, mesh1, mesh2, curve1=6, curve2=-1 )

  !call plot_points_curves ( plot_options, mesh, 'curves3.fig' )
  call delete ( mesh1, mesh3 )

! 4

  call set_mesh_options ( mesh_options, elshape=elshape, oy=R4, nx=n1, ny=n5, &
    lx=H, ly=R5-R4, ratio=[3,3,1,1], factor=[f1,f5,f1,f5] )

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh2, mesh1, mesh, curve1=9, curve2=-1 )

  !call plot_points_curves ( plot_options, mesh, 'curves4.fig' )
  call delete ( mesh1, mesh2 )

! add curves

! free surface
  call add_to_mesh ( mesh, curve=[2,5,8,11] ) ! curve 14 (8)
! bottom
  call add_to_mesh ( mesh, curve=[13,10,7,4] ) ! curve 15 (9)
! walls
  call add_to_mesh ( mesh, curve=[1,12] ) ! curve 16 (10)
! full boundary
  call add_to_mesh ( mesh, curve=[1,2,5,8,11,12,13,10,7,4] ) ! curve 17 (11)
! entry boundary
  call add_to_mesh ( mesh, curve=[10,7] ) ! curve 18


  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves4b.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh4b.fig' )

!  plot_options%xmin = -1
!  plot_options%xmax = 1
!  call plot_mesh ( plot_options, mesh, 'mesh2b_zoom.fig' )

  call write_mesh ( mesh, 'mesh4b.bin' )

  call printinfo ( mesh, printlevel=2 )

  call delete ( mesh )

end program mesh_four_blocks

