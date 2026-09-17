program meshgen_extra_objects3

  use kind_defs_m
  use mesh_m
  use meshgen_m
!  use figplot_m
  use io_utils_m

  implicit none

  type(meshgen_options_t) :: mesh_options

!  type(plot_options_t) :: plot_options
  type(mesh_t) :: mesh, mesh1, mesh2

  integer :: i

  mesh_options%ox = 0.2_dp

  mesh_options%elshape = 6
  mesh_options%nx = 3
  mesh_options%ny = 3

  call quadrilateral2d ( mesh1, mesh_options )

  mesh_options%nx = 6
  mesh_options%ny = 6

  mesh_options%ox = 1.2_dp

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh, nogroupmerge=.true. )

!  call add_to_mesh ( mesh, object='curve', objectcurve=8, typeofobject=2, &
!    excludegroups=(/2/), excludegroups2=(/1/), topology=.true., &
!    intrule=2, nsubint=1 )
  call add_to_mesh ( mesh, object='curve', objectcurve=2, typeofobject=2, &
    excludegroups=[2], excludegroups2=[1], topology=.true., &
    intrule=2, nsubint=2 )

  call delete ( mesh1, mesh2 )

!  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  call fill_mesh_parts ( mesh )

!  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

!  plot_options%objectpointsize=0.4
!  plot_options%objectpointcolor=4
!  call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true. )

  call write_mesh ( mesh, filename='mesh2.out' )

  do i = 1, mesh%nobjects
    print *, 'object = ', i
    print *, 'mesh%objects(:)%refcoor_int'
    print *, mesh%objects(i)%refcoor_int
    print *, 'mesh%objects(:)%grpelm_int'
    print *, mesh%objects(i)%grpelm_int(:,1,:)
    print *, mesh%objects(i)%grpelm_int(:,2,:)
    print *, 'mesh%objects(:)%refcoor2'
    print *, mesh%objects(i)%refcoor2_int
    print *, 'mesh%objects(:)%grpelm2_int'
    print *, mesh%objects(i)%grpelm2_int(:,1,:)
    print *, mesh%objects(i)%grpelm2_int(:,2,:)
  end do

  call delete ( mesh )

end program meshgen_extra_objects3
