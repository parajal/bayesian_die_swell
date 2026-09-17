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

program poisson22

  use tfem_m
  use functions_m

  implicit none

! constants

  integer, parameter :: &
    p = 11,       & ! polynomial order
    nx=4,         & ! number of elements in x
    ny=4            ! number of elements in y

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh

  integer, parameter :: n = 10
  integer :: grpelm(n,2), nnodes, i
  real(dp) :: coor(n,2), refcoor(n,2)

! create mesh

  meshgen_options%elshape = 102
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%p = p

  meshgen_options%regionshape = 3
  meshgen_options%x2d = &
    reshape ( [ 0.0_dp, 2.0_dp, 2.0_dp, -1.0_dp,    &
                 0.0_dp, 0.0_dp, 2.0_dp,  3.0_dp ], &
              [4,2] )
  meshgen_options%curved(1:4) = [ .false., .true., .true., .true. ]
  meshgen_options%funcnr(1:4) = [ 1, 2, 3, 4 ]

  call quadrilateral2d ( mesh, meshgen_options, func=cfunc )

  call fill_mesh_parts ( mesh )

  nnodes=2
  coor(1,:) = [ 0.5_dp, 0.6_dp]
  coor(2,:) = [ 0.545_dp, 1.62_dp]
  call find_refcoor_points ( mesh, nnodes=nnodes, &
    coor=coor, grpelm=grpelm, refcoor=refcoor )

  do i = 1, nnodes
    print *, grpelm(i,:), refcoor(i,:)
  end do

! delete all data including all allocated memory

  call delete ( mesh )

end program poisson22

