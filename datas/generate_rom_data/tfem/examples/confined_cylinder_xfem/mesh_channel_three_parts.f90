!
!  mesh build from 3 submeshes:
!
!         ---------------------------------
!         |        |            |         |
!         |        |            |         |
!         |        |            |         |
!    n3   |        |            |         |
!       0 |   1    |     2      |     3   |
!         |        |            |         |
!    2*H  |        |            |         |
!         |        |            |         |
!         |        |            |         |
!         |        |            |         |
!         ---------------------------------
!                        0
!             L2        2*L1         L2
!             n2         n1          n2
!

! NOTE: one single group is created. Use nogroupmerge=.true. in heading of
! mesh_merge to obtain 2 separate element groups.

program mesh_channel_three_parts

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use meshgen_extra_m
  use figplot_m
  use io_utils_m

  implicit none

  type(meshgen_options_t) :: mesh_options
  type(plot_options_t) :: plot_options
  type(mesh_t) :: mesh, mesh1, mesh2, mesh3

  integer :: elshape=6
  real(dp) :: L1=5, L2=5
  real(dp) :: H=1
  integer :: n1=8, n2=8, n3=12
  real(dp) :: f1=7

! namelist for input of variables; read from standard input

  namelist /meshpar/ L1, L2, H, n1, n2, n3, f1

  open (unit=10, status='old', file='meshinput.txt')
  read (unit=10, nml=meshpar )
  close(unit=10)

  plot_options%fontsize=8


! 1 : center refined

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=n3, &
    ox=-L1, oy=-H, lx=2*L1, ly=2*H )

  call quadrilateral2d ( mesh1, mesh_options )

  !call plot_points_curves ( plot_options, mesh1, 'curves1.fig' )

! 2 : left part

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=n3, &
    ox=-L1-L2, oy=-H, lx=L2, ly=2*H, ratio=[3,0,1,0], &
    factor=[f1,1._dp,f1,1._dp] )

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, curve1=4, curve2=-2 )

  call delete ( mesh1, mesh2 )

  !call plot_points_curves ( plot_options, mesh3, 'curves2.fig' )

! 3 : right part

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=n3, &
    ox=L1, oy=-H, lx=L2, ly=2*H, ratio=[1,0,3,0], &
    factor=[f1,1._dp,f1,1._dp] )

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh3, mesh1, mesh, curve1=2, curve2=-4 )

  call delete ( mesh1, mesh3 )


! add curves

! inflow for periodic b.c.
  call add_to_mesh ( mesh, curve=[-7] ) ! curve 11
! walls
  call add_to_mesh ( mesh, curve=[5,1,8,10,3,6] ) ! curve 12
! full boundary
  call add_to_mesh ( mesh, curve=[5,1,8,9,10,3,6,7] ) ! curve 13


  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

!  plot_options%xmin = -1
!  plot_options%xmax = 1
!  call plot_mesh ( plot_options, mesh, 'mesh_zoom.fig' )

  call write_mesh ( mesh, 'mesh_channel_3p.out' )

  call printinfo ( mesh, printlevel=2 )

  call delete ( mesh )

end program mesh_channel_three_parts

