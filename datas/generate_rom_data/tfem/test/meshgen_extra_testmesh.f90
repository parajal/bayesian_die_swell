program meshgen_extra_testmesh

  use kind_defs_m
  use mesh_m
  use meshgen_m
  !use figplot_m

  implicit none

  type(meshgen_options_t) :: mesh_options
  !type(plot_options_t) :: plot_options
  type(mesh_t) :: mesh, mesh1, mesh2, mesh3
  integer :: mergegroups(1,2) = 1

  mesh_options%elshape = 6
  mesh_options%nx = 2
  mesh_options%ny = 2

  call quadrilateral2d ( mesh1, mesh_options )

  mesh_options%ox = 1._dp

  mesh_options%nx = 4

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, curve1=2, curve2=-4, &
    nogroupmerge=.true. )
  call delete ( mesh1, mesh2 )

  mesh_options%ox = 0._dp
  mesh_options%oy = 1._dp

  mesh_options%nx = 2
  mesh_options%ny = 4

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh3, mesh1, mesh, curve1=3, curve2=-1, &
    mergegroups=mergegroups )
  call delete ( mesh1, mesh3 )

  !call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  call fill_mesh_parts ( mesh )

  !call plot_mesh ( plot_options, mesh, 'mesh.fig' )

  call printa

  call delete ( mesh )

contains

  subroutine printa

  integer :: i

  print *, 'mesh%ndim'
  print *, mesh%ndim
  print *, 'mesh%nnodes'
  print *, mesh%nnodes
  print *, 'mesh%nelem'
  print *, mesh%nelem
  print *, 'mesh%nelgrp'
  print *, mesh%nelgrp
  print *, 'mesh%npoints'
  print *, mesh%npoints
  print *, 'mesh%ncurves'
  print *, mesh%ncurves
  print *, 'mesh%element(1)%elshape'
  print *, mesh%element(1)%elshape
  print *, 'mesh%element(1)%globalshape'
  print *, mesh%element(1)%globalshape
  print *, 'mesh%element(1)%ndim'
  print *, mesh%element(1)%ndim
  print *, 'mesh%element(1)%numnod'
  print *, mesh%element(1)%numnod
  print *, 'mesh%element(1)%sidnumvert'
  print *, mesh%element(1)%sidnumvert
  print *, 'mesh%element(1)%sidvert'
  print *, mesh%element(1)%sidvert
  print *, 'mesh%elnumnod'
  print *, mesh%elnumnod
  print *, 'mesh%grpnumel'
  print *, mesh%grpnumel
  print *, 'mesh%topology(1)%a'
  print *, mesh%topology(1)%a
  print *, 'mesh%nodnumel'
  print *, mesh%nodnumel
  print *, 'mesh%nodelem(:,1)'
  print *, mesh%nodelem(:,1)
  print *, 'mesh%nodelem(:,2)'
  print *, mesh%nodelem(:,2)
  print *, 'mesh%nodnumnod'
  print *, mesh%nodnumnod
  print *, 'mesh%nodnod'
  print *, mesh%nodnod
  print *, 'mesh%points(:mesh%npoints)'
  print *, mesh%points(:mesh%npoints)
  print *, 'mesh%curves(:mesh%ncurves)%ndim'
  print *, mesh%curves(:mesh%ncurves)%ndim
  print *, 'mesh%curves(:mesh%ncurves)%nnodes'
  print *, mesh%curves(:mesh%ncurves)%nnodes
  print *, 'mesh%curves(:mesh%ncurves)%nelem'
  print *, mesh%curves(:mesh%ncurves)%nelem
  print *, 'mesh%curves(:mesh%ncurves)%elnumnod'
  print *, mesh%curves(:mesh%ncurves)%elnumnod
  print *, 'mesh%curves(:mesh%ncurves)%element%numnod'
  print *, mesh%curves(:mesh%ncurves)%element%numnod
  print *, 'mesh%curves(:mesh%ncurves)%element%globalshape'
  print *, mesh%curves(:mesh%ncurves)%element%globalshape
  print *, 'mesh%curves(:mesh%ncurves)%element%elshape'
  print *, mesh%curves(:mesh%ncurves)%element%elshape
  print *, 'mesh%curves%nodes'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%nodes
  end do
  print *, 'mesh%curves%topology(:,:,1)'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%topology(:,:,1)
  end do
  print *, 'mesh%curves%topology(:,:,2)'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%topology(:,:,2)
  end do
  print *, 'mesh%curves%nc%nodnumel'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%nc%nodnumel
  end do
  print *, 'mesh%curves%nc%nodelem'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%nc%nodelem
  end do
  print *, 'mesh%curves%nc%nodnumnod'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%nc%nodnumnod
  end do
  print *, 'mesh%curves%nc%nodnod'
  do i = 1, mesh%ncurves
    print *, mesh%curves(i)%nc%nodnod
  end do
  print *, 'mesh%sidelem(1)%a(:,:,1)'
  print *, mesh%sidelem(1)%a(:,:,1)
  print *, 'mesh%sidelem(1)%a(:,:,2)'
  print *, mesh%sidelem(1)%a(:,:,2)
  print *, 'mesh%sidelem(1)%a(:,:,3)'
  print *, mesh%sidelem(1)%a(:,:,3)
  print *, 'mesh%sidelem(1)%a(:,:,4)'
  print *, mesh%sidelem(1)%a(:,:,4)
  print *, 'mesh%ngluepoints'
  print *, mesh%ngluepoints
  if ( mesh%ngluepoints > 0 ) then
  print *, 'mesh%gluepoints(:,1)'
  print *, mesh%gluepoints(:,1)
  print *, 'mesh%gluepoints(:,2)'
  print *, mesh%gluepoints(:,2)
  end if
  print *, 'mesh%coor'
  print *, mesh%coor

  end subroutine printa

end program meshgen_extra_testmesh

