! generate a 3D mesh of a particle in a sphere
! NOTE: an angular cutout of the cylinder is meshed, with an opening angle
! given by the theta variable

program generate_meshp_3d

  use tfem_m
  use io_utils_m
  use math_defs_m

  implicit none

  real(dp)  :: &
    L1 = 24._dp, &
    L2 = 6._dp,  &
    H = 2._dp,   &
    a = 1._dp,   &
    theta = pi/4._dp

  integer :: &
    nrefine = 0

  type(mesh_t) :: mesh
  type(refinement_fields_t) :: refinement_fields
  real(dp), dimension(:,:), allocatable :: refinement_coor
  real(dp) :: hcyl, hwall, hLwall, hLcen
  integer :: i

  hcyl   = 0.08_dp * 0.5_dp**nrefine
  hwall  = 0.4_dp  * 0.5_dp**nrefine
  hLwall = 1.2_dp  * 0.5_dp**nrefine
  hLcen  = 0.16_dp * 0.5_dp**nrefine

  call add_refinement_field ( refinement_fields, &
    coor=reshape([0._dp,0._dp,0._dp],[1,3]), distmin=1.1*a, &
    distmax=3*a, dx_fine=hcyl, dx_coarse=hLwall )

  call add_refinement_field ( refinement_fields, &
    coor=reshape([cos(theta/2)*H,sin(theta/2)*H,0._dp],[1,3]), &
    distmin=2*a, distmax=4*a, dx_fine=hwall, dx_coarse=hLwall )

  allocate(refinement_coor(100,3))

  refinement_coor = 0
  refinement_coor(:,3) = [(-L2 + (i-1)*(L2+L2)/(100-1),i=1,100)]

  call add_refinement_field ( refinement_fields, &
    coor=refinement_coor, distmin=0.2*a, distmax=3*a, &
    dx_fine=hLcen, dx_coarse=hLwall )

  open ( unit=25, file='meshp_3d.geo' )

  call write_refinement_fields ( refinement_fields, 'meshp_3d.geo' )

  write ( 25, '(1X,A,F25.14,A)' ) 'L1 = ', L1, ';'
  write ( 25, '(1X,A,F25.14,A)' ) 'L2 = ', L2, ';'
  write ( 25, '(1X,A,F25.14,A)' ) 'H = ', H, ';'
  write ( 25, '(1X,A,F25.14,A)' ) 'a = ', a, ';'
  write ( 25, '(1X,A,F25.14,A)' ) 'angular_opening = ', theta, ';'

  write ( 25, '(/1x,a)' ) 'Include "settling_sphere_3d_occ.igo";'
  write ( 25, '(/1x,a)' ) 'Include "refinement.igo";'

  close ( 25 )

  call delete_refinement_fields ( refinement_fields )

  call execute_command_line ('gmsh -3 -order 2 -o &
    &meshp_3d.msh meshp_3d.geo > outputmesh.out' )

  call read_mesh_gmsh ( mesh, ndim=3, physgeom=.true., &
    filename='meshp_3d.msh', sortphys=.true. )

  call fill_mesh_parts ( mesh )

  call printinfo ( mesh, printlevel=2 )

  ! call write_mesh_vtk ( mesh, filename='mesh.vtk' )
  ! call write_geometries_vtk ( mesh, write_normals=.true. )

  call write_mesh_gmsh ( mesh, filename='meshp_3d.msh' )

end program generate_meshp_3d
