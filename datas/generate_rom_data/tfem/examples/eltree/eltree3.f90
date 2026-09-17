! Example program for using the eltree module
! Use eltree (quadtree) together with aligned triangles/tets near the
! interface to subdivide the region.
! Compute an integral on a circle or sphere for which the exact solution is
! known

module functions_m

  use math_defs_m

  implicit none

contains

  function levelset ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1)) :: levelset

!   circle with center at (0,0) and radius 1
    levelset = sqrt(sum(x**2,dim=2))-1._dp

  end function levelset

  function mapcoor ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1),size(x,2)) :: mapcoor

!   map reference [-1:1]x[-1,1](x[-1,1]) on [-2:2]x[-2,2](x[-2,2])
    mapcoor = 2 * x

  end function mapcoor

end module functions_m

program eltree3

  use eltree_m
  use functions_m
  use figplot_m
  use io_utils_m
  use gauss_m
  use limits_m

  implicit none

! constants

  integer, parameter :: &
    numsplit = 5,  &   ! relative size of smallest elements near
                       ! interface 1/2^numsplit
    numsplitmin = 1, & ! relative size of all subelements will be at least
                       ! as small as 1/2^numsplitmin
    intrule = 2, &! integration rule on one subelement of eltree
    intrule_s = 3, &! integration rule on one subelement of submesh
    ndim = 2      ! dimension of space
!    intrule_s = 4, &! integration rule on one subelement of submesh
!    ndim = 3      ! dimension of space

  real(dp), parameter :: &
    epsjac = 0._dp, & ! size of jacobian for removal of elements
    split_threshold = 1.e-2_dp    ! levelset value for being "close"
                                  ! enough to the interface for tree splitting

  logical, parameter :: &
    plot = .true.   ! make plots

! definitions

  type(eltree_t) :: eltree
  type(mesh_t) :: mesh
  type(plot_options_t) :: plot_options
  type(gauss_t) :: gauss
  type(integration_options_t) :: intopt

  integer :: ninti, nintp, nintis
  real(dp) :: coor(2,ndim), vole, intge, vol, intg, detJ
  real(dp), allocatable :: w(:), x(:,:), we(:), xe(:,:), ws(:), xs(:,:)


  if ( ndim == 3 ) SUBDIVIDE_HEX6 = .true.

! fill root node of the eltree

  coor(1,:) = -1._dp  ! lower corner of the reference region (-1,-1,-1)
  coor(2,:) =  1._dp  ! upper corner of the reference region (1,1,1)

  call fill_node_eltree ( eltree, coor )

! divide root reference domain into subdomains

  call subdivide ( eltree, levelset=levelset, numsplit=numsplit, &
    split_threshold=split_threshold, numsplitmin=numsplitmin, &
    mapcoor=mapcoor, submesh=.true. )

  print *, 'number of subelements = ', number_of_subelements(eltree)
  print *, 'number of subelements outside = ', &
            number_of_subelements(eltree,lsign=[1])
  print *, 'number of subelements crossing interface = ', &
            number_of_subelements(eltree,lsign=[0])
  print *, 'number of subelements inside = ', &
            number_of_subelements(eltree,lsign=[-1])
  print *, 'reference volume of elements outside = ', volume(eltree,lsign=[1])
  print *, 'reference volume of elements crossing interface = ', &
            volume(eltree,lsign=[0])
  print *, 'reference volume of elements inside = ', volume(eltree,lsign=[-1])

! integration scheme on subelements in eltree

  if ( ndim == 2 ) then
    gauss%globalshape='quadrilateral'
  else if ( ndim == 3 ) then
    gauss%globalshape='hexahedron'
  end if
  gauss%intrule=intrule

  call set_ninti ( gauss, ninti )

  allocate ( we(ninti), xe(ninti,ndim) )

  call set_Gauss_integration ( gauss, xe, we )

! integration scheme on subelements in submesh

  if ( ndim == 2 ) then
    gauss%globalshape='triangle'
  else if ( ndim == 3 ) then
    gauss%globalshape='tetrahedron'
  end if
  gauss%intrule=intrule_s

  call set_ninti ( gauss, nintis )

  allocate ( ws(nintis), xs(nintis,ndim) )

  call set_Gauss_integration ( gauss, xs, ws )

! composite integration scheme for the inside of circle/sphere

  intopt%epsjac = epsjac

  nintp = number_of_integration_points( eltree, lsign=[-1], &
             ninti=ninti, nintis=nintis, integration_options=intopt )

  allocate ( w(nintp), x(nintp,ndim) )

  call integration_points ( eltree, x, w, lsign=[-1], ninti=ninti, &
    xe=xe, we=we, nintis=nintis, xs=xs, ws=ws, integration_options=intopt )

! compute volume/area and the integral of r^2.

  x = mapcoor(x)

  vol = sum(w)
  intg = sum(w*(sum(x**2,dim=2)))

  if ( ndim == 2 ) then
    detJ = 4._dp ! determinant of the Jacobian of the mapping in mapcoor
                 ! (relative volume change)
    vole = pi
    intge = pi/2
  else if ( ndim == 3 ) then
    detJ = 8._dp ! determinant of the Jacobian of the mapping in mapcoor
                 ! (relative volume change)
    vole = 4*pi/3
    intge = 4*pi/5
  end if

  vol = vol * detJ
  intg = intg * detJ

  print *, 'computed volume of circle/sphere = ', vol, &
           'relative error = ', (vol-vole)/vole
  print *, 'computed integral of r^2 on circle/sphere = ', intg, &
           'relative error = ', (intg-intge)/intge

! plotting

  if ( ndim == 2 .and. plot ) then

    plot_options%boundarycolor=0

    call eltree_to_mesh ( eltree, mesh, lsign=[1], mapcoor=mapcoor )
    call fill_mesh_parts ( mesh )
    call plot_mesh ( plot_options, mesh, filename='mesh1.fig' )
    call delete(mesh)
    call eltree_to_mesh ( eltree, mesh, lsign=[-1], mapcoor=mapcoor )
    call fill_mesh_parts ( mesh )
    call plot_mesh ( plot_options, mesh, filename='mesh3.fig' )
    call delete(mesh)
    call eltree_to_mesh ( eltree, mesh, mapcoor=mapcoor )
    call fill_mesh_parts ( mesh )
    call plot_mesh ( plot_options, mesh, filename='mesh.fig' )

    call write_mesh_vtk ( mesh, filename='mesh.vtk' )

  else if ( plot ) then

    call eltree_to_mesh ( eltree, mesh, lsign=[1], mapcoor=mapcoor )
    call fill_mesh_parts ( mesh )
    call write_mesh_vtk ( mesh, filename='mesh1.vtk' )
    call delete(mesh)
    call eltree_to_mesh ( eltree, mesh, lsign=[-1], mapcoor=mapcoor )
    call fill_mesh_parts ( mesh )
    call write_mesh_vtk ( mesh, filename='mesh3.vtk' )
    call delete(mesh)
    call eltree_to_mesh ( eltree, mesh, mapcoor=mapcoor )
    call fill_mesh_parts ( mesh )
    call write_mesh_vtk ( mesh, filename='mesh.vtk' )

  end if

! release memory

  deallocate ( we, xe, w, x, ws, xs )

  call delete(eltree)
  if ( plot ) call delete(mesh)

end program eltree3
