! A curved quadrilateral region meshed using spectral elements.
! Ratio and factor is used to refine the elements towards a boundary.
! mesh_convert is used to make a mesh that can be plotted.

module functions_m

  use kind_defs_m

  implicit none

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

end module functions_m

program mesh_spectral_quad

  use tfem_m
  use figplot_m
  use functions_m

  implicit none

! constants

  integer, parameter :: &
    p  = 7,  &   ! polynomial order of spectral elements
    nx = 5,  &   ! number of elements in x-direction
    ny = 5       ! number of elements in y-direction

! definitions

  type(mesh_t) :: mesh, mesh_plot
  type(meshgen_options_t) :: meshgen_options
  type(plot_options_t) :: plot_options

! Create spectral mesh

  call set_mesh_options ( meshgen_options, elshape=102, nx=nx, ny=ny, p=p, &
    regionshape=3, ratio=[0,5,0,6], factor=[1._dp,0.1_dp,1._dp,0.1_dp], &
    funcnr=[1,2,3,4] )
  meshgen_options%curved = .true.

  call quadrilateral2d ( mesh, meshgen_options, func=cfunc )

  call fill_mesh_parts( mesh )

! Create plot mesh

  call mesh_convert ( mesh, mesh_plot, elementshapes='spectraltolinear', &
    warn=.false. )

  call fill_mesh_parts( mesh_plot )

  call plot_points_curves ( plot_options, mesh_plot, filename='curves.fig' )
  call plot_mesh ( plot_options, mesh_plot, filename='mesh_plot.fig' )

end program mesh_spectral_quad
