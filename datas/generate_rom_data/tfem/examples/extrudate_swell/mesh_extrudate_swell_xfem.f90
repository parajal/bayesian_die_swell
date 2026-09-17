
!  mesh build from 3 submeshes:
!
!                  --------------
!                  |            |
!              n4  |            |
!                  |     3      |
!                S |            |
!                  |            |
!         -----------------------
!         |        |            |
!         |        |            |
!         |        |            |
!     n3  |        |            |
!         |   1    |     2      |
!       R |        |            |
!         |        |            |
!         |        |            |
!         |        |            |
!         |        |            |
!         -----------------------
!        -L1       0            L2
!             n1         n2
!

! NOTE: one single group is created. Use nogroupmerge=.true. in heading of
! mesh_merge to obtain 2 separate element groups.

program mesh_extrudate_swell_xfem

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
  real(dp) :: R=1, S=0.5
  integer :: n1=8, n2=8, n3=12, n4=5

  real(dp) :: f1=7, f2=6, f3=8, f4=7


! namelist for input of variables; read from standard input

  namelist /meshpar/ elshape, L1, L2, R, S, n1, n2, n3, n4, f1, f2, f3, f4

  read (unit=*,nml=meshpar )


! 1

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=n3, ox=-L1, &
    lx=L1, ly=R, ratio=[3,3,1,1], factor=[f1,f3,f1,f3] )

  call quadrilateral2d ( mesh1, mesh_options )
  !call plot_points_curves ( plot_options, mesh1, 'curves1.fig' )

! 2

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=n3, &
    lx=L2, ly=R, ratio=[1,3,3,1], factor=[f2,f3,f2,f3] )

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, curve1=2, curve2=-4 )

  plot_options%fontsize=8
  !call plot_points_curves ( plot_options, mesh3, 'curves2.fig' )
  call delete ( mesh1, mesh2 )

! 3

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=n4, &
    lx=L2, ly=S, oy=R, ratio=[1,1,3,3], factor=[f2,f4,f2,f4] )

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh3, mesh1, mesh, curve1=7, curve2=-1 )

  call delete ( mesh1, mesh3 )


! add curves

! outflow
  call add_to_mesh ( mesh, curve=[6,8] ) ! curve 11
! centerline
  call add_to_mesh ( mesh, curve=[1,5] ) ! curve 12
! full boundary
  call add_to_mesh ( mesh, curve=[1,5,6,8,9,10,3,4] ) ! curve 13


  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves_esx.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh_esx.fig' )

!  plot_options%xmin = -1
!  plot_options%xmax = 1
!  call plot_mesh ( plot_options, mesh, 'mesh_esx_zoom.fig' )

  call write_mesh ( mesh, 'mesh_esx.out' )

  call printinfo ( mesh, printlevel=2 )

  call delete ( mesh )

end program mesh_extrudate_swell_xfem

