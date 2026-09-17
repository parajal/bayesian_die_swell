! simple illustration for the use of the module submesh (part of eltree)

program eltree2

  use tfem_m
  use eltree_m
  use figplot_m
  use gauss_m
  use io_utils_m
  implicit none

  type(smesh_t) :: smesh
  type(mesh_t) :: mesh
  type(plot_options_t) :: plot_options
  type(gauss_t) :: gauss

  integer, parameter :: intrule = 24  ! integration rule for the tet
  integer :: ninti
  real(dp) :: x(8,3), s(8)
  real(dp) :: xe(intrule,3), we(intrule)
  real(dp), allocatable :: xg(:,:), w(:)

! define the vertices of unit hexahedron

  x(1,:) = [0._dp,0._dp,0._dp]
  x(2,:) = [1._dp,0._dp,0._dp]
  x(3,:) = [1._dp,1._dp,0._dp]
  x(4,:) = [0._dp,1._dp,0._dp]
  x(5,:) = [0._dp,0._dp,1._dp]
  x(6,:) = [1._dp,0._dp,1._dp]
  x(7,:) = [1._dp,1._dp,1._dp]
  x(8,:) = [0._dp,1._dp,1._dp]

! define the levelset value in the vertices

  s = [ -1._dp, 0.1_dp, 2._dp, 0.1_dp, 0.1_dp, 2._dp, 3._dp, 2._dp ]

  call submesh_hexahedron ( x, s, smesh )

  print *, smesh%topology
  print *, transpose(smesh%coor)
  print *, smesh%lsign

! composite integration rule on the unit hexahedron (positive levelset only)

  ninti = number_of_integration_points ( smesh, ninti=size(we), lsign=[1], &
                                         epsjac=0._dp )

  allocate ( w(ninti), xg(ninti,3) )

  gauss%globalshape = 'tetrahedron'
  gauss%intrule = size(we)

  call set_Gauss_integration ( gauss, xe, we )

  call integration_points ( smesh, xg, w, ninti=size(we), xe=xe, we=we, &
                            lsign=[1], epsjac=0._dp )

! output

  print *, 'number of integration points = ', ninti
  print *, 'area = ', sum(w), &
    'integral x^2y^2z^2 dV on the positive levelset = ', &
    sum(xg(:,1)**2*xg(:,2)**2*xg(:,3)**2*w)

  call smesh_to_mesh ( smesh, mesh )
  call fill_mesh_parts ( mesh )
  call plot_mesh ( plot_options, mesh, filename="mesh.fig" )

  call write_mesh_vtk ( mesh, filename="mesh.vtk" )

  call delete(mesh)
  call smesh_to_mesh ( smesh, mesh, lsign=[-1] )
  call fill_mesh_parts ( mesh )
  call plot_mesh ( plot_options, mesh, filename="meshm1.fig" )

  call write_mesh_vtk ( mesh, filename="meshm1.vtk" )

  call delete(mesh)
  call smesh_to_mesh ( smesh, mesh, lsign=[1] )
  call fill_mesh_parts ( mesh )
  call plot_mesh ( plot_options, mesh, filename="meshp1.fig" )

  call write_mesh_vtk ( mesh, filename="meshp1.vtk" )

  call delete(smesh)
  call delete(mesh)

end program eltree2
