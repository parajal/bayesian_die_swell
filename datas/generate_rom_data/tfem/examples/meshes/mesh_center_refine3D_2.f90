! 3D mesh, with refinement near the middle
! similar to mesh_center_refine3D.f90, but now using add_to_mesh with
! matchingsurface to change orientation of the surfaces.

program mesh_center_refine3D_2

  use tfem_m
  use figplot_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    nx=5,               & ! number of elements in x
    ny=5,               & ! number of elements in y
    nz=5                  ! number of elements in z

! definitions

  type(mesh_t) :: mesh, mesh1, mesh2, mesh3
  type(plot_options_t) :: plot_options
  type(meshgen_options_t) :: meshgen_options

  real(dp) :: x3d(8,3)


  plot_options%viewpoint=[1.,0.8,0.4]

! create mesh

  x3d = reshape ( [ 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp,    &
                    0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp,    &
                    0.25_dp, 0.25_dp, 0.75_dp, 0.75_dp,    &
                    -0.75_dp, -0.75_dp, 1.75_dp, 1.75_dp,    &
                    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,    &
                    1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp ], [8,3] )

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, regionshape=5, &
    elshape=14, x3d=x3d )

! basic part
  call hexahedron ( mesh1, meshgen_options )

  call plot_points_curves ( plot_options, mesh1, 'curves1.fig' )

! rotate 90 degrees
  call mesh_convert ( mesh1, mesh2, coordinates=[1,3,-2] )

  mesh2%coor(:,2) = mesh2%coor(:,2) + 0.75_dp
  mesh2%coor(:,3) = mesh2%coor(:,3) + 0.25_dp

  call mesh_merge ( mesh1, mesh2, mesh3, surface1=4, surface2=2 )

  call delete ( mesh1, mesh2 )

! shift center to (0,0)
  mesh3%coor(:,2) = mesh3%coor(:,2) - 0.5_dp
  mesh3%coor(:,3) = mesh3%coor(:,3) + 0.25_dp

  call plot_points_curves ( plot_options, mesh3, 'curves2.fig' )

! rotate 180 degrees
  call mesh_convert ( mesh3, mesh1, coordinates=[1,-2,-3] )

  call mesh_merge ( mesh3, mesh1, mesh2, surfaces1=[9,2], surfaces2=[2,9] )

  call plot_points_curves ( plot_options, mesh2, 'curves3.fig' )

  call delete ( mesh3, mesh1 )

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, regionshape=4, &
    elshape=14, ly=0.5_dp, lz=0.5_dp )

! center part
  call hexahedron ( mesh3, meshgen_options )

! shift center to (0,0)
  mesh3%coor(:,2) = mesh3%coor(:,2) - 0.25_dp
  mesh3%coor(:,3) = mesh3%coor(:,3) - 0.25_dp

! change orientation of surfaces of center to match with outer part
  call add_to_mesh ( mesh3, matchingsurface=[1,12], matchingmesh=mesh2, &
                     replace=1 )
  call add_to_mesh ( mesh3, matchingsurface=[4,7], matchingmesh=mesh2, &
                     replace=4 )

  call mesh_merge ( mesh2, mesh3, mesh, surfaces1=[12,1,7,17], &
    surfaces2=[1,6,4,2] )

  call delete ( mesh3, mesh2 )

  call fill_mesh_parts ( mesh )

  call write_mesh_vtk ( mesh, filename='mesh.vtk' )

! plot curves, surfaces and mesh

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig', &
    surfaces=[6,11,3,18,21,13,8] )

  call printinfo ( mesh )

  call delete ( mesh )

end program mesh_center_refine3D_2
