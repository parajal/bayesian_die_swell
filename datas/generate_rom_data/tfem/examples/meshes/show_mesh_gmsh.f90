program show_mesh_gmsh

  use meshgen_m
  use io_utils_m

  implicit none

  type(mesh_t) :: mesh

  call read_mesh_gmsh ( mesh, filename='ellipsoid_in_a_box.msh' )

  call fill_mesh_parts ( mesh )

  call write_mesh_vtk ( mesh, 'mesh.vtk' )

  call delete ( mesh )

end program show_mesh_gmsh

