program contraction3D

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use io_utils_m

  implicit none

  type(meshgen_options_t) :: mesh_options
!  type(plot_options_t) :: plot_options
  type(mesh_t) :: mesh, mesh1, mesh2, mesh3
  real(dp) :: x3d(8,3)

  integer, parameter :: elshape=14
  integer, parameter :: n1=2, n2=2, n3=2, n4=3, n5=3, m1=2, m2=3, nw=2

  real(dp), parameter :: factor1=0.5, factor2=0.2
  real(dp), parameter :: L1=10, L2=1, L3=1, L4=2, L5=10, W=5
  real(dp), parameter :: H1=1.5, H2=2.5_dp, H3=1.0_dp

  integer :: i, j, k

!
!  constraction build from 7 submeshes:
!
!
!         ------------------
!         |        |       |
!  m2  H2 |   2    |   4   |
!         |        |       |
!         |--------|\____  |
!         |        |     \_|---------------------------
!  m1  H1 |   1    |   3   |    |    |                |  H3
!         |        |       | 5  | 6  |      7         |
!       0 ---------------------------------------------
!        -L1      -L2      0    L3   L4               L5
!             n1      n2     n3   n4        n5
!
! NOTE: one single group is created. Use nogroupmerge=.false. in heading of
! mesh_merge to obtain 7 separate element groups.


! 1

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=nw, nz=m1, &
    ox=-L1, lx=L1-L2, ly=W, lz=H1, regionshape=4 )

  call hexahedron ( mesh1, mesh_options )

!  call plot_points_curves ( plot_options, mesh1, 'curves1.fig' )

! 2

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=nw, nz=m2, &
    ox=-L1, oz=H1, lx=L1-L2, ly=W, lz=H2, regionshape=4 )

  call hexahedron ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, surface1=6, surface2=1 )

!  call plot_points_curves ( plot_options, mesh3, 'curves2.fig' )
  call delete ( mesh1, mesh2 )

! 3

  x3d(1,:) = [-L2,0._dp,0._dp]
  x3d(2,:) = [0,0,0]
  x3d(6,:) = [0._dp,0._dp,H3]
  x3d(5,:) = [-L2,0._dp,H1]
  x3d(4,:) = [-L2,W,0._dp]
  x3d(3,:) = [0._dp,W,0._dp]
  x3d(7,:) = [0._dp,W,H3]
  x3d(8,:) = [-L2,W,H1]

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=nw, nz=m1, &
    x3d=x3d, regionshape=5, ratio=[6,0,6,0,0,6,6,0,6,0,6,0] )

  mesh_options%factor = factor1

  call hexahedron ( mesh2, mesh_options )

  call mesh_merge ( mesh3, mesh2, mesh1, surface1=3, surface2=5 )
!  plot_options%fontsize=8
!  call plot_points_curves ( plot_options, mesh1, 'curves3.fig' )
  call delete ( mesh3, mesh2 )

! 4

  x3d(1,:) = [-L2,0._dp,H1]
  x3d(2,:) = [0._dp,0._dp,H3]
  x3d(6,:) = [0._dp,0._dp,H1+H2]
  x3d(5,:) = [-L2,0._dp,H1+H2]
  x3d(4,:) = [-L2,W,H1]
  x3d(3,:) = [0._dp,W,H3]
  x3d(7,:) = [0._dp,W,H1+H2]
  x3d(8,:) = [-L2,W,H1+H2]

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=nw, nz=m2, &
    x3d=x3d, regionshape=5, ratio=[6,0,6,0,0,5,5,0,6,0,6,0] )

  mesh_options%factor = factor1
  mesh_options%factor(2) = factor2
  mesh_options%factor(10) = factor2

  call hexahedron ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, surfaces1=[16,8], &
     surfaces2=[1,5] )
!  call plot_points_curves ( plot_options, mesh3, 'curves4.fig' )
  call delete ( mesh1, mesh2 )

! 5

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n3, ny=nw, nz=m1, &
    lx=L3, ly=W, lz=H3, regionshape=4, ratio=[5,0,5,0,6,6,6,6,5,0,5,0] )

  mesh_options%factor = factor1

  call hexahedron ( mesh1, mesh_options )

  call mesh_merge ( mesh3, mesh1, mesh2, surface1=14, surface2=5 )
!  plot_options%fontsize=5
!  call plot_points_curves ( plot_options, mesh2, 'curves5.fig' )
  call delete ( mesh3, mesh1 )

! 6

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n4, ny=nw, nz=m1, &
    ox=L3, lx=L4-L3, ly=W, lz=H3, regionshape=4, &
    ratio=[0,0,0,0,6,0,0,6,0,0,0,0] )

  mesh_options%factor = factor1

  call hexahedron ( mesh1, mesh_options )

  call mesh_merge ( mesh2, mesh1, mesh3, surface1=23, surface2=5 )
!  call plot_points_curves ( plot_options, mesh3, 'curves6.fig' )
  call delete ( mesh2, mesh1 )

! 7

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n5, ny=nw, nz=m1, &
    ox=L4, lx=L5-L4, ly=W, lz=H3, regionshape=4 )

  call hexahedron ( mesh1, mesh_options )

  call mesh_merge ( mesh3, mesh1, mesh, surface1=28, surface2=5 )
!  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call delete ( mesh3, mesh1 )

  call fill_mesh_parts ( mesh )

!  call plot_mesh ( plot_options, mesh, 'mesh.fig', &
!    surfaces=(/11,20,18,25,30,35,33,34,29,24,15,19,4,9/) )

!  call write_mesh_vtk ( mesh, filename='mesh.vtk' )

  call printinfo ( mesh )

  do i = 1, mesh%ncurves
    print *, 'curve=', i
    print *, 'mesh%curves(i)%nummeshelem'
    print *, mesh%curves(i)%nummeshelem
    do j = 1, mesh%curves(i)%nelem
      print *, (mesh%curves(i)%meshelem(j)%a(1,k),k=1,mesh%curves(i)%nummeshelem(j))
      print *, (mesh%curves(i)%meshelem(j)%a(2,k),k=1,mesh%curves(i)%nummeshelem(j))
      print *, (mesh%curves(i)%meshelem(j)%a(3,k),k=1,mesh%curves(i)%nummeshelem(j))
    end do
  end do
  do i = 1, mesh%nsurfaces
    print *, 'surfaces=', i
    print *, 'mesh%surfaces(i)%nummeshelem'
    print *, mesh%surfaces(i)%nummeshelem
    do j = 1, mesh%surfaces(i)%nelem
      print *, (mesh%surfaces(i)%meshelem(j)%a(1,k),k=1,mesh%surfaces(i)%nummeshelem(j))
      print *, (mesh%surfaces(i)%meshelem(j)%a(2,k),k=1,mesh%surfaces(i)%nummeshelem(j))
      print *, (mesh%surfaces(i)%meshelem(j)%a(3,k),k=1,mesh%surfaces(i)%nummeshelem(j))
    end do
  end do

  call delete ( mesh )

end program contraction3D

