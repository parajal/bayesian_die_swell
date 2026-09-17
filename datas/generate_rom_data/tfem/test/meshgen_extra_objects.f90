program meshgen_extra_objects

  use kind_defs_m
  use mesh_m
  use meshgen_m
  !use figplot_m

  implicit none

  type(meshgen_options_t) :: mesh_options

  !type(plot_options_t) :: plot_options
  type(mesh_t) :: mesh, mesh1, mesh2, mesh3

  integer :: i

  mesh_options%elshape = 6
  mesh_options%nx = 4
  mesh_options%ny = 4

  call quadrilateral2d ( mesh1, mesh_options )

  mesh_options%nx = 4
  mesh_options%ny = 4

  mesh_options%ox = 0.4_dp

  mesh_options%lx = 0.2_dp
  mesh_options%ly = 0.8_dp

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh, nogroupmerge=.true. )

  mesh_options%elshape = 5
  mesh_options%ny = 17


  call quadrilateral2d ( mesh3, mesh_options )

  call add_to_mesh ( mesh, object='mesh', objectmesh=mesh2, typeofobject=2, &
    excludegroups=[2], excludegroups2=[1] )

  call add_to_mesh ( mesh, object='curve', objectcurve=8, typeofobject=2, &
    excludegroups=[2], excludegroups2=[1] )

  call add_to_mesh ( mesh, object='curve', objectcurve=2, objectmesh=mesh3, &
    typeofobject=2, excludegroups=[2], excludegroups2=[1] )

  call delete ( mesh1, mesh2, mesh3 )

  !call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  call fill_mesh_parts ( mesh )

  !call plot_mesh ( plot_options, mesh, 'mesh.fig' )

  !plot_options%objectpointsize=0.4
  !plot_options%objectpointcolor=4
  !call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true., object1=1 )
  !plot_options%objectpointcolor=0
  !call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true., object1=2 )
  !plot_options%objectpointcolor=2
  !call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true., object1=3 )

  do i = 1, mesh%nobjects
    print *, 'object = ', i
    print *, 'mesh%objects(:)%refcoor'
    print *, mesh%objects(i)%refcoor
    print *, 'mesh%objects(:)%grpelm'
    print *, mesh%objects(i)%grpelm(:,1)
    print *, mesh%objects(i)%grpelm(:,2)
    print *, 'mesh%objects(:)%refcoor2'
    print *, mesh%objects(i)%refcoor2
    print *, 'mesh%objects(:)%grpelm2'
    print *, mesh%objects(i)%grpelm2(:,1)
    print *, mesh%objects(i)%grpelm2(:,2)
  end do

  call delete ( mesh )

end program meshgen_extra_objects
