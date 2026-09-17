!
!  mesh build from 2 submeshes:
!
!
!     y3  -----------------------
!         |                     |
!         |                     |
!     n3  |         2           |
!         |                     |
!         |                     |
!     y2  -----------------------
!         |                     |
!         |                     |
!     n2  |         1           |
!         |                     |
!         |                     |
!     y1  -----------------------
!         x1                    x2
!                  n1
!

! NOTE: one single group is created. Use nogroupmerge=.true. in heading of
! mesh_merge to obtain 2 separate element groups.

program mesh_two_blocks

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use figplot_m
  use io_utils_m

  implicit none

  type(meshgen_options_t) :: mesh_options
  type(plot_options_t) :: plot_options
  type(mesh_t) :: mesh, mesh1, mesh2

  logical :: rotate=.false.
  integer :: elshape=6
  real(dp) :: x1=0, x2=5
  real(dp) :: y1=1, y2=3, y3=5
  integer :: n1=12, n2=6, n3=6

  real(dp) :: f1=3, f2=3, f3=3


! namelist for input of variables; read from standard input

  namelist /meshpar/ elshape, x1, x2, y1, y2, y3, n1, n2, n3, f1, f2, f3, &
    rotate

  read ( unit=*, nml=meshpar )


! 1

 call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=n2, ox=x1, &
    oy=y1, lx=x2-x1, ly=y2-y1, ratio=[3,1,1,3], factor=[f1,f2,f1,f2] )

  call quadrilateral2d ( mesh1, mesh_options )
  !call plot_points_curves ( plot_options, mesh1, 'curves1.fig' )

! 2

  call set_mesh_options ( mesh_options, elshape=elshape, ox=x1, oy=y2, &
    nx=n1, ny=n3, lx=x2-x1, ly=y3-y2, ratio=[3,3,1,1], &
    factor=[f1,f3,f1,f3] )

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh, curve1=3, curve2=-1 )

  plot_options%fontsize=8
  !call plot_points_curves ( plot_options, mesh, 'curves2.fig' )
  call delete ( mesh1, mesh2 )

! add new curves

! left boundary
  call add_to_mesh ( mesh, curve=[2,5] ) ! curve 8
! right boundary
  call add_to_mesh ( mesh, curve=[7,4] ) ! curve 9
! bottom and top boundary
  call add_to_mesh ( mesh, curve=[1,6] ) ! curve 10
! full boundary
  call add_to_mesh ( mesh, curve=[1,2,5,6,7,4] ) ! curve 11

  if ( rotate ) then

!   rotate mesh 90 degrees clockwise

    call mesh_convert ( mesh, mesh1, coordinates=[2,-1] )
    call delete ( mesh )
    call copy ( mesh1, mesh )
    call delete ( mesh1 )

  end if

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves2b.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh2b.fig' )

!  plot_options%xmin = -1
!  plot_options%xmax = 1
!  call plot_mesh ( plot_options, mesh, 'mesh2b_zoom.fig' )

  call write_mesh ( mesh, 'mesh2b.out' )

  call printinfo ( mesh, printlevel=2 )

  call delete ( mesh )

end program mesh_two_blocks

