module functions_m

  use math_defs_m

  implicit none

  real(dp) :: rp, xp(2)

contains

  function cfunc ( nr, xr )
    integer, intent(in) :: nr
    real(dp), intent(in) :: xr
    real(dp), dimension(2) :: cfunc

    select case(nr)
      case(1)
        cfunc = [ 2*xr, 5*xr*(1-xr) ]
      case(2)
        cfunc = [ 2+xr*(1-xr), 2*xr ]
      case(3)
        cfunc = [ -1+3*(1-xr), 2+xr*(1-xr) ]
      case(4)
        cfunc = [ xr-1+xr*(1-xr), 3*(1-xr) ]
      case default
        write(*,'(/a,i0/)') 'Error cfunc: wrong function number: ', nr
        stop
    end select

  end function cfunc

! objectscoor defines the coordinates of the objects

  subroutine objectscoor ( objectnr, coor )
    integer, intent(in) :: objectnr
    real(dp), dimension(:,:), intent(inout) :: coor

    integer :: np, i
    real(dp) :: p(size(coor,1))

    np = size(coor,1)
    p = [(2*pi/np*(i-1),i=1,np)]

    coor(:,1) = rp * sin(p) + xp(1)
    coor(:,2) = rp * cos(p) + xp(2)

  end subroutine objectscoor

end module functions_m

program objects_in_mesh

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use figplot_m
  use functions_m

  implicit none

  type(meshgen_options_t) :: mesh_options
  type(plot_options_t) :: plot_options
  type(mesh_t) :: mesh, mesh1

  integer :: nx, ny

  mesh_options%elshape = 6
  mesh_options%regionshape = 3
  mesh_options%x2d = &
    reshape ( [ 0.0_dp, 2.0_dp, 2.0_dp, -1.0_dp,    &
                0.0_dp, 0.0_dp, 2.0_dp,  3.0_dp ], &
              [4,2] )
  nx = 20
  ny = 20
  mesh_options%nx = nx
  mesh_options%ny = ny
  mesh_options%curved(1) = .true.
  mesh_options%funcnr(1) = 1
  mesh_options%curved(2) = .true.
  mesh_options%funcnr(2) = 2
  mesh_options%curved(3) = .true.
  mesh_options%funcnr(3) = 3
  mesh_options%curved(4) = .true.
  mesh_options%funcnr(4) = 4

  call quadrilateral2d ( mesh, mesh_options, func=cfunc )

! add object defined by square mesh

  call set_mesh_options ( mesh_options, elshape=5, nx=3, ny=6, &
    ox=1.5_dp, oy=1.5_dp, lx=0.3_dp, ly=0.3_dp )
  call quadrilateral2d ( mesh1, mesh_options )
  call add_to_mesh ( mesh, object='mesh', objectmesh=mesh1 )
  call delete ( mesh1 )

! add object defined by points on a circle

  rp = 0.3_dp ! radius of the circular object
  xp = [0.2_dp,1.5_dp] ! center of the object

  call add_to_mesh ( mesh, object='coordinates', nnodes=20, &
    objectsub=objectscoor, objectcoornr=1 )

! other parts of the mesh

  call fill_mesh_parts ( mesh )

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )
  plot_options%objectpointcolor=4
  plot_options%objectpointsize=0.4
  call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true. )

  call delete ( mesh )

end program objects_in_mesh
