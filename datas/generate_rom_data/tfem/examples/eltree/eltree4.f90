! Example program for using the eltree module
! Use eltree (quadtree/octree) together with aligned triangles/tets near the
! interface to subdivide the region.
! Compute an integral on a (the boundary of a) circle or sphere for which
! the exact solution is known

module functions_m

  use math_defs_m

  implicit none

contains

  function levelset ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1)) :: levelset

!   circle with center at (0,0) and radius 1
    levelset = sqrt(sum(x**2,dim=2))-1.0_dp

  end function levelset

  function mapcoor ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1),size(x,2)) :: mapcoor

!   map reference [-1:1]x[-1,1](x[-1,1]) on [-2:2]x[-2,2](x[-2,2])
    mapcoor = 2 * x

  end function mapcoor

end module functions_m

program eltree4

  use eltree_m
  use functions_m
  use figplot_m
  use io_utils_m
  use gauss_m
  use limits_m

  implicit none

! constants

  integer, parameter :: &
    numsplit = 4,  &   ! relative size of smallest elements near
                       ! interface 1/2^numsplit
    numsplitmin = 1, & ! relative size of all subelements will be at
                       ! least as small as 1/2^numsplitmin
    !intrule = 2, &! integration rule on one interface element of eltree
    !ndim = 2      ! dimension of space
    intrule = 3, &! integration rule on one interface element of eltree
    ndim = 3      ! dimension of space

  real(dp), parameter :: &
    split_threshold = 0._dp    ! levelset value for being "close" enough
                               ! to the interface for tree splitting

  logical, parameter :: &
    plot = .true.   ! make plots

! definitions

  type(eltree_t) :: eltree
  type(mesh_t) :: mesh, mesh_new
  type(plot_options_t) :: plot_options
  type(gauss_t) :: gauss

  integer :: ninti, nintp
  real(dp) :: coor(2,ndim), areae, intge, areac, intg, scaling, intgn(ndim)
  real(dp), allocatable :: w(:), x(:,:), wn(:,:), we(:), xe(:,:)

  logical, parameter :: &
  remove_nodes = .false. ! remove double nodes from the surface mesh
                         ! (make it conforming).

! generate a continuous levelset=0 surface in 3D (see module limits_m)

  if ( ndim == 3 ) SUBDIVIDE_HEX6 = .true.

! fill root node of the eltree

  coor(1,:) = -1._dp  ! lower corner of the reference region (-1,-1,-1)
  coor(2,:) =  1._dp  ! upper corner of the reference region (1,1,1)

  call fill_node_eltree ( eltree, coor )

! divide root reference domain into subdomains

  call subdivide ( eltree, levelset=levelset, numsplit=numsplit, &
    split_threshold=split_threshold, numsplitmin=numsplitmin, mapcoor=mapcoor, &
    intmesh=.true., submesh=.true. )

  print *, 'number of interface_elements = ', &
            number_of_interface_elements(eltree)
  print *, 'reference area of interface elements = ', &
            area(eltree)

! integration scheme on interface elements in eltree

  if ( ndim == 2 ) then
    gauss%globalshape='line'
  else if ( ndim == 3 ) then
    gauss%globalshape='triangle'
  end if
  gauss%intrule=intrule

  call set_ninti ( gauss, ninti )

  allocate ( we(ninti), xe(ninti,ndim-1) )

  call set_Gauss_integration ( gauss, xe, we )

! composite integration scheme for the interface (circle/sphere)

  nintp = number_of_interface_integration_points( eltree, ninti=ninti )

  allocate ( w(nintp), x(nintp,ndim), wn(nintp,ndim) )

  call interface_integration_points ( eltree, x, w, wn, ninti=ninti, &
                                      xe=xe, we=we )

! compute length/area, the integral of x^2 and the integral of n
! (outside normal).

  if ( ndim == 2 ) then
    scaling = 2._dp ! scaling from the mapping in mapcoor
                    ! (relative length change)
    areae = 2*pi
    intge = pi
  else if ( ndim == 3 ) then
    scaling = 4._dp ! scaling from the mapping in mapcoor (relative area change)
    areae = 4*pi
    intge = 4*pi/3
  end if

  x = mapcoor(x)

  areac = scaling * sum(w)
  intg = scaling * sum(w*x(:,1)**2)
  intgn = scaling * sum(wn)

  print *, 'computed area of circle/sphere = ', areac, &
           'relative error = ', (areac-areae)/areae
  print *, 'computed integral of x^2 on circle/sphere = ', intg, &
           'relative error = ', (intg-intge)/intge
  print *, 'computed normal circle/sphere = ', intgn

! plotting

  if ( ndim == 2 .and. plot ) then

    plot_options%boundarycolor=0

    call eltree_to_interface_mesh ( eltree, mesh, mapcoor=mapcoor )
    call fill_mesh_parts ( mesh )
    call plot_mesh ( plot_options, mesh, filename='mesh.fig' )

    call write_mesh_vtk ( mesh, filename='mesh.vtk' )

  else if ( plot ) then

    call eltree_to_interface_mesh ( eltree, mesh, mapcoor=mapcoor )
    call fill_mesh_parts ( mesh )
    call write_mesh_vtk ( mesh, filename='mesh.vtk' )

    if ( remove_nodes ) then

      print *, 'number of nodes unconforming surface mesh =', mesh%nnodes
      call add_to_mesh ( mesh, nodblocks=[10,10,10] )
      !EPSNODEBOX=1e-3_dp
      !EPSNODEOVERLAP=1e-3_dp
      call mesh_convert ( mesh, mesh_new, remove_double_nodes=.true., &
        warn=.false. )
      print *, 'number of nodes conforming surface mesh =', mesh_new%nnodes
      call fill_mesh_parts ( mesh_new )
      call write_mesh_vtk ( mesh_new, filename='mesh_new.vtk' )
      call delete ( mesh_new )

    end if

  end if

! release memory

  deallocate ( we, xe, w, x, wn )

  call delete(eltree)
  if ( plot ) call delete(mesh)

end program eltree4

