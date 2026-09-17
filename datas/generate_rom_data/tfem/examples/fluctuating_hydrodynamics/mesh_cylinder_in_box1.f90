!
!  cylinder confined between four walls build from 1 quarter submesh and
!  rotating once over 90 degrees and merging. Then rotating over 180 degrees
!  and merging.
!
!              |\                         !
!              |  \                       !
!              |    \                     !
!              |      \                   !
!              |        \                 !
!              |          \               !
!              |          /
!              |         /
!              |        /
!              |        |
!           -H |        |a       o
!              |        | n2
!              |        |
!              |        \                 !
!              |         \                !
!              |          \               !
!              |          /
!              |    n1  /
!              |      /
!              |   /
!              |/

! NOTE: one single group is created. Use nogroupmerge=.false. in heading of
! mesh_merge to obtain 6 separate element groups.

module cylinder_functions_m

  use math_defs_m

  implicit none

  real(dp) :: a = 1

contains

  function cfunc ( nr, xr )
    integer, intent(in) :: nr
    real(dp), intent(in) :: xr
    real(dp), dimension(2) :: cfunc

    select case(nr)
      case(1)
        cfunc = a * [ -cos(pi*xr/2-pi/4), sin(pi*xr/2-pi/4) ]
      case default
        write(*,'(/a,i0/)') 'Error cfunc: wrong function number: ', nr
        stop
    end select

  end function cfunc

end module cylinder_functions_m

program mesh_cylinder_in_box1

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
  real(dp) :: H=4
  integer :: n1=8, n2=8

  real(dp) :: f1=15


! namelist for input of variables; read from standard input

  namelist /meshpar/ elshape, a, H, n1, n2, f1

  read ( unit=*, nml=meshpar )


  srt2 = sqrt(2._dp)/2


! 1

  x2d(1,:) = [-H,-H]
  x2d(2,:) = [-a*srt2,-a*srt2]
  x2d(3,:) = [-a*srt2,a*srt2]
  x2d(4,:) = [-H,H]

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=n2, &
    x2d=x2d, regionshape=3,ratio=[3,0,1,0], factor=[f1,1._dp,f1,1._dp] )

  mesh_options%curved(2) = .true.
  mesh_options%funcnr(2) = 1

  call quadrilateral2d ( mesh1, mesh_options, func=cfunc )

  plot_options%fontsize=8
  call plot_points_curves ( plot_options, mesh1, 'curves1.fig' )

! rotate 90 deg

  call mesh_convert ( mesh1, mesh2, coordinates=[2,-1] )

  call mesh_merge ( mesh1, mesh2, mesh3, curve1=3, curve2=-1 )

  call plot_points_curves ( plot_options, mesh3, 'curves2.fig' )

  call delete ( mesh1, mesh2 )

! rotate 180 deg

  call mesh_convert ( mesh3, mesh1, coordinates=[-1,-2] )

  call mesh_merge ( mesh3, mesh1, mesh, curves1=[1,6], curves2=[-6,-1] )

  call plot_points_curves ( plot_options, mesh, 'curves3.fig' )

  call delete ( mesh3, mesh1 )

! add curves

! cylinder

  call add_to_mesh ( mesh, curve=[2,5,8,11] ) ! curve 13

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh1.fig' )

  call write_mesh ( mesh, 'mesh1.out' )

  call printinfo ( mesh )

  call delete ( mesh )

end program mesh_cylinder_in_box1

