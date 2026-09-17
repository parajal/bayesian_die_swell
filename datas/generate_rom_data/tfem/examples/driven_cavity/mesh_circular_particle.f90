! module for creating a mesh for a circular particle.

module mesh_circular_particle_m

  use math_defs_m

  implicit none

  real(dp), private :: rpc = 1._dp, xpc(2) = 0

contains

! create mesh for particle domain

  subroutine create_mesh_circular_particle ( mesh_particle, elshape, nx, ny, rp, xp )

    use meshgen_m

    implicit none

    type(mesh_t), intent(inout) :: mesh_particle
    integer, intent(in) :: elshape, nx, ny
    real(dp), intent(in) :: rp, xp(2)


    type(mesh_t) :: mesh1, mesh2, mesh3
    type(meshgen_options_t) :: meshgen_options

!    integer :: elshape = 6
!    integer :: nx, ny
    real(dp) :: srt2
    real(dp) :: lx, ly, x2d(4,2)

    srt2 = sqrt(2._dp)/2._dp

!   set function parameters
    rpc = rp
    xpc = xp

!    nx = 3
!    ny = 5
    lx = rp * 0.45
    ly = rp * 0.45

    ! 1
    x2d(1,:) = [-lx,-ly] + [xp(1),xp(2)]
    x2d(2,:) = [-rp*srt2,-rp*srt2] + [xp(1),xp(2)]
    x2d(3,:) = [-rp*srt2,rp*srt2] + [xp(1),xp(2)]
    x2d(4,:) = [-lx,ly] + [xp(1),xp(2)]

    call set_mesh_options ( meshgen_options, elshape=elshape, nx=nx, ny=ny, x2d=x2d, regionshape=3 )

    meshgen_options%curved(2) = .true.
    meshgen_options%funcnr(2) = 1

    call quadrilateral2d ( mesh1, meshgen_options, func=cfunc )

    ! 2
    x2d(1,:) = [-lx,ly] + [xp(1),xp(2)]
    x2d(2,:) = [-rp*srt2,rp*srt2] + [xp(1),xp(2)]
    x2d(3,:) = [rp*srt2,rp*srt2] + [xp(1),xp(2)]
    x2d(4,:) = [lx,ly] + [xp(1),xp(2)]

    call set_mesh_options ( meshgen_options, elshape=elshape, nx=nx, ny=ny, x2d=x2d, regionshape=3 )

    meshgen_options%curved(2) = .true.
    meshgen_options%funcnr(2) = 2

    call quadrilateral2d ( mesh2, meshgen_options, func=cfunc )

    call mesh_merge ( mesh1, mesh2, mesh3, curve1=3, curve2=-1 )

    call delete ( mesh1, mesh2 )

    ! 3
    x2d(1,:) = [lx,ly] + [xp(1),xp(2)]
    x2d(2,:) = [rp*srt2,rp*srt2] + [xp(1),xp(2)]
    x2d(3,:) = [rp*srt2,-rp*srt2] + [xp(1),xp(2)]
    x2d(4,:) = [lx,-ly] + [xp(1),xp(2)]

    call set_mesh_options ( meshgen_options, elshape=elshape, nx=nx, ny=ny, x2d=x2d, regionshape=3 )

    meshgen_options%curved(2) = .true.
    meshgen_options%funcnr(2) = 3

    call quadrilateral2d ( mesh1, meshgen_options, func=cfunc )

    call mesh_merge ( mesh3, mesh1, mesh2, curve1=6, curve2=-1 )

    call delete ( mesh1, mesh3 )

    ! 4
    x2d(1,:) = [lx,-ly] + [xp(1),xp(2)]
    x2d(2,:) = [rp*srt2,-rp*srt2] + [xp(1),xp(2)]
    x2d(3,:) = [-rp*srt2,-rp*srt2] + [xp(1),xp(2)]
    x2d(4,:) = [-lx,-ly] + [xp(1),xp(2)]

    call set_mesh_options ( meshgen_options, elshape=elshape, nx=nx, ny=ny, x2d=x2d, regionshape=3)

    meshgen_options%curved(2) = .true.
    meshgen_options%funcnr(2) = 4

    call quadrilateral2d ( mesh1, meshgen_options, func=cfunc )

    call mesh_merge ( mesh2, mesh1, mesh3, curves1=[9,1], curves2=[-1,-3] )

    call delete ( mesh1, mesh2 )

    ! 5
    x2d(1,:) = [-lx,-ly] + [xp(1),xp(2)]
    x2d(2,:) = [-lx,ly] + [xp(1),xp(2)]
    x2d(3,:) = [lx,ly] + [xp(1),xp(2)]
    x2d(4,:) = [lx,-ly] + [xp(1),xp(2)]

    call set_mesh_options ( meshgen_options, elshape=elshape, nx=ny, ny=ny, x2d=x2d, regionshape=2 )

    call quadrilateral2d ( mesh1, meshgen_options )

    call mesh_merge ( mesh3, mesh1, mesh_particle, curves1=[12,10,7,4], curves2=[-4,-3,-2,-1] )

    call delete ( mesh3, mesh1 )

  end subroutine create_mesh_circular_particle


! Function for the circular boundary of the particle (used in create_particle_mesh).

  function cfunc ( nr, xr )
    integer, intent(in) :: nr
    real(dp), intent(in) :: xr
    real(dp), dimension(2) :: cfunc

    select case(nr)
      case(1)
        cfunc = rpc * [ -cos(pi*xr/2-pi/4), sin(pi*xr/2-pi/4) ] + xpc
      case(2)
        cfunc = rpc * [ sin(pi*xr/2-pi/4), cos(pi*xr/2-pi/4) ] + xpc
      case(3)
        cfunc = rpc * [ cos(pi*xr/2-pi/4), -sin(pi*xr/2-pi/4) ] + xpc
      case(4)
        cfunc = rpc * [ -sin(pi*xr/2-pi/4), -cos(pi*xr/2-pi/4) ] + xpc
      case default
        write(*,'(/a,i0/)') 'Error cfunc: wrong function number: ', nr
        stop
    end select

  end function cfunc

end module mesh_circular_particle_m
