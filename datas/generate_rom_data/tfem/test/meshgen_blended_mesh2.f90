program meshgen_blended_mesh2

  use tfem_m
  use io_utils_m

  implicit none

  type(mesh_t) :: mesh_P1, mesh_P2, mesh_P3, mesh1, mesh2, mesh3

!  call read_mesh_gmsh ( mesh_P1, filename='mesh_P1.msh', ndim=2 )
  call read_mesh_gmsh ( mesh_P2, filename='mesh_P2.msh', ndim=2 )
  call read_mesh_gmsh ( mesh_P3, filename='mesh_P3.msh', ndim=2 )
  call mesh_convert ( mesh_P2, mesh_P1, elementshapes='quadratictolinear' )

  print *, 'mesh_P1'
  print *
  call printinfo ( mesh_P1, printlevel=6 )
  print *, 'mesh_P2'
  print *
  call printinfo ( mesh_P2, printlevel=6 )
  print *, 'mesh_P3'
  print *
  call printinfo ( mesh_P3, printlevel=6 )

  call mesh_convert ( mesh_P3, mesh1, blendmesh=mesh_P2 )

  print *, 'mesh1'
  print *
  call printinfo ( mesh1, printlevel=6 )

  call mesh_convert ( mesh1, mesh2, blendmesh=mesh_P1 )

  print *, 'mesh2'
  print *
  call printinfo ( mesh2, printlevel=6 )

  call mesh_convert ( mesh2, mesh3, blendmesh=mesh_P3 )

  print *, 'mesh3'
  print *
  call printinfo ( mesh3, printlevel=6 )

  call fill_mesh_parts ( mesh3 )

  print *, 'mesh3%elnumnod'
  print *, mesh3%elnumnod

  call delete ( mesh_P1, mesh_P2, mesh_P3, mesh1, mesh2, mesh3 )

end program meshgen_blended_mesh2
