!
!  mesh build from 3 submeshes:
!
!
!         --------------------------------
!         |        |            |        |
!         |        |            |        |
!         |        |            |        |
!     n4  |        |            |        |
!         |   1    |     2      |   3    |
!       R |        |            |        |
!         |        |            |        |
!         |        |            |        |
!         |        |            |        |
!         |        |            |        |
!         --------------------------------
!        -L1       -L2          0        L3
!             n1         n2         n3
!
! NOTE: one groups is created. Use nogroupmerge=.true. in heading of
! mesh_merge to obtain three separate groups.

program mesh_three_blocks

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
  real(dp) :: L1=7, L2=5, L3=5
  real(dp) :: R=1
  integer :: n1=8, n2=8, n3=12, n4=12

  real(dp) :: f1=1, f2=6, f3=8, f4=8


! namelist for input of variables; read from standard input

  namelist /meshpar/ elshape, L1, L2, L3, R, n1, n2, n3, n4, f1, f2, f3, f4

  read ( unit=*, nml=meshpar )


! 1

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=n4, ox=-L1, &
    lx=L1-L2, ly=R, ratio=[1,3,3,1], factor=[f1,f4,f1,f4] )

  call quadrilateral2d ( mesh1, mesh_options )
  !call plot_points_curves ( plot_options, mesh1, 'curves1.fig' )

! 2

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=n4, &
    ox=-L2, lx=L2, ly=R, ratio=[3,3,1,1], factor=[f2,f4,f2,f4] )

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, curve1=2, curve2=-4 )

  plot_options%fontsize=8
  !call plot_points_curves ( plot_options, mesh3, 'curves2.fig' )
  call delete ( mesh1, mesh2 )

! 3

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n3, ny=n4, &
    lx=L3, ly=R, ratio=[1,3,3,1], factor=[f3,f4,f3,f4] )

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh3, mesh1, mesh, curve1=6, curve2=-4 )
  call delete ( mesh3, mesh1 )

  plot_options%fontsize=8
  !call plot_points_curves ( plot_options, mesh, 'curves3.fig' )

! add curves

! wall
  call add_to_mesh ( mesh, curve=[-3,-7] ) ! curve 11
! wall+slip surface
  call add_to_mesh ( mesh, curve=[-3,-7,-10] ) ! curve 12
! centerline
  call add_to_mesh ( mesh, curve=[1,5,8] ) ! curve 13
! invert curve 2
  call add_to_mesh ( mesh, curve=[-2] ) ! curve 14
! full boundary
  call add_to_mesh ( mesh, curve=[1,5,8,9,10,7,3,4] ) ! curve 15


  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves3b.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh3b.fig' )

  plot_options%xmin = -1
  plot_options%xmax = 1
  call plot_mesh ( plot_options, mesh, 'mesh3b_zoom.fig' )

  plot_options%xmin = -L1
  plot_options%xmax = -L1+R
  call plot_mesh ( plot_options, mesh, 'mesh3b_zoom1.fig' )

  call write_mesh ( mesh, 'mesh3b.out' )

  call printinfo ( mesh, printlevel=2 )

  call delete ( mesh )

end program mesh_three_blocks

