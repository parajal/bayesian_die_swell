! generate a 2D mesh of a particle in a sphere

program generate_meshp

  use tfem_m
  use io_utils_m
  use math_defs_m

  implicit none

  real(dp)  :: &
    L1 = 24._dp, &
    L2 = 6._dp,  &
    H = 2._dp,   &
    a = 1._dp

  integer :: &
    nrefine = 0

  type(mesh_t) :: mesh
  real(dp) :: hcyl, hwall, hLwall, hLcen

  hcyl   = 0.08_dp * 0.5_dp**nrefine
  hwall  = 0.4_dp  * 0.5_dp**nrefine
  hLwall = 1.2_dp  * 0.5_dp**nrefine
  hLcen  = 0.08_dp * 0.5_dp**nrefine

  open ( unit=25, file='meshp.geo' )

  write ( 25, '(1X,A,F25.14,A)' ) 'hcyl = ', hcyl, ';'
  write ( 25, '(1X,A,F25.14,A)' ) 'hwall = ', hwall , ';'
  write ( 25, '(1X,A,F25.14,A)' ) 'hLwall = ', hLwall, ';'
  write ( 25, '(1X,A,F25.14,A)' ) 'hLcen = ', hLcen, ';'

  write ( 25, '(1X,A,F25.14,A)' ) 'L1 = ', L1, ';'
  write ( 25, '(1X,A,F25.14,A)' ) 'L2 = ', L2, ';'
  write ( 25, '(1X,A,F25.14,A)' ) 'H = ', H, ';'
  write ( 25, '(1X,A,F25.14,A)' ) 'a = ', a, ';'

  write ( 25, '(/1x,a)' ) 'Mesh.Algorithm = 6;'
  write ( 25, '(/1x,a)' ) 'Include "settling_sphere.igo";'

  close ( 25 )

  call execute_command_line ('gmsh -2 -order 2 -o &
    &meshp.msh meshp.geo > outputmesh.out' )

  call read_mesh_gmsh ( mesh, ndim=2, filename='meshp.msh' )

  call fill_mesh_parts ( mesh )

  call printinfo ( mesh, printlevel=2 )

  ! call write_mesh_vtk ( mesh, filename='mesh.vtk' )
  ! call write_geometries_vtk ( mesh, write_normals=.true. )

end program generate_meshp
