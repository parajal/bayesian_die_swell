module cylinder_functions_m

  use math_defs_m

  implicit none

  real(dp), parameter :: a = 1

contains

  function cfunc ( nr, xr )
    integer, intent(in) :: nr
    real(dp), intent(in) :: xr
    real(dp), dimension(2) :: cfunc

    select case(nr)
      case(1)
        cfunc = a * [ -cos(pi*xr/4), sin(pi*xr/4) ]
      case(2)
        cfunc = a * [ -cos(pi*(xr+1)/4), sin(pi*(xr+1)/4) ]
      case(3)
        cfunc = a * [ sin(pi*xr/4), cos(pi*xr/4) ]
      case(4)
        cfunc = a * [ sin(pi*(xr+1)/4), cos(pi*(xr+1)/4) ]
      case default
        write(*,'(/a,i0/)') 'Error cfunc: wrong function number: ', nr
        stop
    end select

  end function cfunc

end module cylinder_functions_m

program mesh_confined_cylinder

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use figplot_m
  use io_utils_m
  use cylinder_functions_m

  implicit none

  type(meshgen_options_t) :: mesh_options
  type(plot_options_t) :: plot_options
  type(mesh_t) :: mesh, mesh1, mesh2, mesh3
  real(dp) :: x2d(4,2), srt2

  integer :: elshape=6
  real(dp) :: L1=15, L2=2
  real(dp) :: H=2
  integer :: n1=8, n2=8, n3=12, n4=6

  real(dp) :: f1=7, f2=6, f3=8
  real(dp) :: fH=12, fdiag=15

!
!  cylinder confined between two walls build from 6 submeshes:
!
!
!         ----------------------------------------------
!         |        |  \         |        /    |        |
!         |        |   \    3   |   4   /     |        |
!         |        |     \      |      /      |        |
!         |        |      \     |     /       |        |
!         |   1    |  2    \ n4 |    /    5   |    6   |
!     H   |        |        --------          |        |
!         |        |      /           \       |        |
!         |        |     /             \      |        |
!         |        |    / n3            \     |        |
!         |        |    |               |     |        |
!         --------------                ----------------
!        -L1      -L2  -a       0       a     L2       L1
!             n1     n2
!

! NOTE: one single group is created. Use nogroupmerge=.false. in heading of
! mesh_merge to obtain 6 separate element groups.


! namelist for input of variables; read from standard input

  namelist /meshpar/ elshape, L1, L2, H, n1, n2, n3, n4, f1, f2, f3, fH, fdiag

  read ( unit=*, nml=meshpar )


  srt2 = sqrt(2._dp)/2


! 1

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=n3, ox=-L1, &
    lx=L1-L2, ly=H, ratio=[3,1,1,3], factor=[f1,fH,f1,fH] )

  call quadrilateral2d ( mesh1, mesh_options )
  call plot_points_curves ( plot_options, mesh1, 'curves1.fig' )

! 2

  x2d(1,:) = [-L2,0._dp]
  x2d(2,:) = [-a,0._dp]
  x2d(3,:) = [-a*srt2,a*srt2]
  x2d(4,:) = [-L2,H]

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=n3, &
    x2d=x2d, regionshape=3, ratio=[3,1,1,3], factor=[f2,f3,fdiag,fH] )

  mesh_options%curved(2) = .true.
  mesh_options%funcnr(2) = 1

  call quadrilateral2d ( mesh2, mesh_options, func=cfunc )

  call mesh_merge ( mesh1, mesh2, mesh3, curve1=2, curve2=-4 )

  plot_options%fontsize=8
  call plot_points_curves ( plot_options, mesh3, 'curves2.fig' )
  call delete ( mesh1, mesh2 )

! 3

  x2d(1,:) = [-a*srt2,a*srt2]
  x2d(2,:) = [0._dp,a]
  x2d(3,:) = [0._dp,H]
  x2d(4,:) = [-L2,H]

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n4, ny=n2, &
    x2d=x2d, regionshape=3, ratio=[0,1,0,3], factor=[1._dp,f2,1._dp,fdiag] )

  mesh_options%curved(1) = .true.
  mesh_options%funcnr(1) = 2

  call quadrilateral2d ( mesh2, mesh_options, func=cfunc )

  call mesh_merge ( mesh3, mesh2, mesh1, curve1=7, curve2=-4 )

  call plot_points_curves ( plot_options, mesh1, 'curves3.fig' )
  call delete ( mesh3, mesh2 )

! 4

  x2d(1,:) = [0._dp,a]
  x2d(2,:) = [a*srt2,a*srt2]
  x2d(3,:) = [L2,H]
  x2d(4,:) = [0._dp,H]

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n4, ny=n2, &
    x2d=x2d, regionshape=3, ratio=[0,1,0,3], factor=[1._dp,fdiag,1._dp,f2] )

  mesh_options%curved(1) = .true.
  mesh_options%funcnr(1) = 3

  call quadrilateral2d ( mesh2, mesh_options, func=cfunc )

  call mesh_merge ( mesh1, mesh2, mesh3, curve1=9, curve2=-4 )
  call plot_points_curves ( plot_options, mesh3, 'curves4.fig' )
  call delete ( mesh1, mesh2 )

! 5

  x2d(1,:) = [a,0._dp]
  x2d(2,:) = [L2,0._dp]
  x2d(3,:) = [L2,H]
  x2d(4,:) = [a*srt2,a*srt2]

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=n3, &
    x2d=x2d, regionshape=3, ratio=[1,1,3,3], factor=[f2,fH,fdiag,f3] )

  mesh_options%curved(4) = .true.
  mesh_options%funcnr(4) = 4

  call quadrilateral2d ( mesh1, mesh_options, func=cfunc )

  call mesh_merge ( mesh3, mesh1, mesh2, curve1=12, curve2=-3 )
  plot_options%fontsize=5
  call plot_points_curves ( plot_options, mesh2, 'curves5.fig' )
  call delete ( mesh3, mesh1 )

! 6

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=n3, ox=L2, &
    lx=L1-L2, ly=H, ratio=[1,1,3,3], factor=[f1,fH,f1,fH] )

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh2, mesh1, mesh, curve1=15, curve2=-4 )
  call plot_points_curves ( plot_options, mesh, 'curves6.fig' )
  call delete ( mesh2, mesh1 )

! add curves

! cylinder
  call add_to_mesh ( mesh, curve=[6,8,11,16] ) ! curve 20
! top
  call add_to_mesh ( mesh, curve=[19,13,10,3] ) ! curve 21
! centerline+cylinder
  call add_to_mesh ( mesh, curve=[1,5,20,14,17] ) ! curve 22
! centerline left
  call add_to_mesh ( mesh, curve=[1,5] ) ! curve 23
! centerline right
  call add_to_mesh ( mesh, curve=[14,17] ) ! curve 24

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

  call write_mesh ( mesh, 'mesh.out' )

  print *, 'number of elements: ', mesh%nelem
  print *, 'number of nodes:    ', mesh%nnodes

  call delete ( mesh )

  call set_mesh_options ( mesh_options )

end program mesh_confined_cylinder

