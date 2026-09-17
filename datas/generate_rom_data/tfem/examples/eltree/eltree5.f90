! Example program for using the eltree module
! A 2D region is divided into finite elements.
! Use eltree (quadtree) together with aligned triangles near the
! interface to subdivide the elements.
! Compute an integral on a (the boundary of a) circle for which the exact
! solution is known

module functions_sd5_m

  use tfem_elem_m
  use math_defs_m

  implicit none

  save

! mesh used for mapping of reference element
  type(mesh_t), pointer :: lmesh => null()

  integer :: lelem, lelgrp
! radius and center of the cylinder
  real(dp) :: rpl = 1._dp, cpl(2) = [0._dp,0._dp]

contains


! levelset for the circle

  function levelset ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1)) :: levelset

    real(dp) :: c(size(x,1),size(x,2))

    c(:,1) = cpl(1)
    c(:,2) = cpl(2)

!   circle with center at cpl and radius rpl
    levelset = sqrt(sum((x-c)**2,dim=2)) - rpl

  end function levelset


! function for mapping reference coordinates to the real coordinates
! in building the eltree within an element.

  function mapcoor ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1),size(x,2)) :: mapcoor

    real(dp) :: phi(size(x,1),9), xnod(9,2)

    call shape_quad_Q2 ( x, phi )

    call get_coordinates ( lmesh, lelgrp, lelem, xnod )

    mapcoor = matmul ( phi, xnod ) ! isoparametric

  end function mapcoor

  function cfunc ( nr, xr )
    integer, intent(in) :: nr
    real(dp), intent(in) :: xr
    real(dp), dimension(2) :: cfunc

    select case(nr)
      case(1)
        cfunc = [ 2*(2*xr-1), -2+2.0_dp*xr*(1-xr) ]
      case(2)
        cfunc = [ 2+2.0_dp*xr*(1-xr), 2*(2*xr-1) ]
      case(3)
        cfunc = [ -2*(2*xr-1), 2-2.0_dp*xr*(1-xr) ]
      case(4)
        cfunc = [ -2+2.0_dp*xr*(1-xr), -2*(2*xr-1) ]
      case default
        write(*,'(/a,i0/)') 'Error cfunc: wrong function number: ', nr
        stop
    end select

  end function cfunc

end module functions_sd5_m


module subs_sd5_m

  use tfem_elem_m
  use eltree_m
  use convection_diffusion_set_globals_m
  use functions_sd5_m

  implicit none

  save

! user gauss is based on a quadtree using the module eltree
  type(eltree_t) :: eltree

! the levelset function in all the nodes (it is a signed distance function
! in this case of a simple cylinder).
  real(dp), dimension(:), allocatable :: d

! the integration scheme on interface elements in the eltree:
!   intrule_ie: the integration rule
!   ninti_ie: number of integration points
!   x_ie, w_ie: points and weights
  integer :: intrule_ie, ninti_ie
  real(dp), allocatable :: w_ie(:), x_ie(:,:)

! wng: weight * reference normal in the integration points on the interface
! dan: weight * normal in the integration points on the interface
!     (after mapping )
! da: weight in the integration points on the interface (after mapping)
  real(dp), allocatable :: wng(:,:), dan(:,:), da(:)

! split parameters
!  ns : numsplit parameter in subdivide for the integration scheme
!  nsmin : numsplitmin parameter in subdivide
!  ns_plot : numsplit parameter in subdivide for the plot mesh
!  spl_th: split_threshold parameter in subdivide
  integer :: ns = 2, nsmin = 1, ns_plot = 2
  real(dp) :: spl_th = 0._dp

! some logicals
  logical :: printnumel = .false., printint = .false.

contains

! Element routine for use with "integrate" in postprocessing

  subroutine integrate_interface ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec )

    use poisson_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    logical :: interface_in_element
    integer :: ip
    integer :: nod(mesh%element(elgrp)%numnod)
    real(dp) :: coor(2,2)
    type(gauss_t) :: gauss_ie

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_poisson ( mesh, coefficients, elgrp )

      allocate ( x(nodalp,ndim) )

      gauss_ie%globalshape = 'line'
      gauss_ie%intrule = intrule_ie

      call set_ninti ( gauss_ie, ninti_ie )

      allocate ( w_ie(ninti_ie), x_ie(ninti_ie,1) )

      call set_Gauss_integration ( gauss_ie, x_ie, w_ie )

    end if

!   nodal points

    nod = mesh%topology(elgrp)%a(:,elem)

!   interface in element?

    if ( all ( d(nod) > spl_th ) .or. all ( d(nod) < -spl_th ) ) then

!     all nodes far from the interface

      interface_in_element = .false.

    else

!     element possibly contains an interface, start subdivide

!     fill root node of the eltree

      coor(1,:) = -1._dp  ! lower corner of the reference region (-1,-1,-1)
      coor(2,:) =  1._dp  ! upper corner of the reference region (1,1,1)

      call fill_node_eltree ( eltree, coor )

!     divide root reference domain into subdomains

      lelgrp = elgrp
      lelem = elem

      call subdivide ( eltree, levelset=levelset, numsplit=ns, &
        split_threshold=spl_th, numsplitmin=nsmin, mapcoor=mapcoor, intmesh=.true. )

      if ( printnumel ) then
        print *, 'number of interface_elements = ', &
          number_of_interface_elements(eltree)
      end if

!     determine ninti

      ninti = number_of_interface_integration_points ( eltree, ninti=ninti_ie )

      if ( printint ) then
        print *, 'number of integration points = ', ninti
      end if

      if ( ninti > 0 ) then

        interface_in_element = .true.

      else

        interface_in_element = .false.

        call delete(eltree)

      end if

    end if

!   integration

    if ( interface_in_element ) then

!     element crosses the interface

!     allocate arrays

      allocate ( xig(ninti,ndim), wg(ninti), wng(ninti,ndim) )

!     composite integration scheme on the interface

      call interface_integration_points ( eltree, xig, wg, wng, ninti=ninti_ie,&
        xe=x_ie, we=w_ie )

      call delete(eltree)

!     allocate arrays

      allocate ( dan(ninti,ndim), da(ninti) )

      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim) )

      allocate ( detF(ninti), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( xg(ninti,ndim) )

!     set shape function

      call set_shape_function ( shapefunc, xig, phi, dphi )

!     element geometric info

      call get_coordinates ( mesh, elgrp, elem, x )

      call isoparametric_deformation ( x, dphi, F, Finv, detF )

      xg = matmul ( phi, x )

      if ( coorsys == 1 ) then
        detF = 2 * pi * xg(:,2) * detF
      end if

      call shape_derivative ( dphi, Finv, dphidx )

!     transform reference area*normal to actual area*normal: da = J F^{-T} dA

      do ip = 1, ninti
        dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
        da(ip) = sqrt(dot_product(dan(ip,:),dan(ip,:)))
      end do

!     integral over the surface

      elemvec(1) = sum(da)
      elemvec(2) = sum(xg(:,1)**2*da)
      elemvec(3:ndim+2) = sum(dan,dim=1)

!     deallocate arrays

      deallocate ( xig, wg, wng )

      deallocate ( dan, da )

      deallocate ( phi, dphi )

      deallocate ( detF, F )
      deallocate ( Finv, dphidx )
      deallocate ( xg )

    else

!     normal element inside or outside

      elemvec = 0

    end if

    if ( last ) then

!     last element in this group

      deallocate ( x )

      deallocate ( w_ie, x_ie )

    end if

  end subroutine integrate_interface


! define the mesh for each element (interface)

  subroutine userelmesh1 ( mesh, elgrp, elem, indicator, elmesh )

    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: elgrp, elem
    integer, intent(out) :: indicator
    type(mesh_t), intent(inout) :: elmesh

    logical :: interface_in_element
    integer :: nod(mesh%element(elgrp)%numnod)
    real(dp) :: coor(2,2)


    nod = mesh%topology(elgrp)%a(:,elem)

!   interface in element?

    if ( all ( d(nod) > spl_th ) .or. all ( d(nod) < -spl_th ) ) then

!     all nodes far from the interface

      interface_in_element = .false.

    else

!     element possibly contains an interface, start subdivide

!     fill root node of the eltree

      coor(1,:) = -1._dp  ! lower corner of the reference region (-1,-1,-1)
      coor(2,:) =  1._dp  ! upper corner of the reference region (1,1,1)

      call fill_node_eltree ( eltree, coor )

!     divide root reference domain into subdomains

      lelgrp = elgrp
      lelem = elem

      call subdivide ( eltree, levelset=levelset, numsplit=ns_plot, &
        split_threshold=spl_th, numsplitmin=nsmin, mapcoor=mapcoor, intmesh=.true. )

      if ( number_of_interface_elements(eltree) > 0 ) then

        interface_in_element = .true.

      else

        interface_in_element = .false.

        call delete(eltree)

      end if

    end if

!   add elements

    if ( interface_in_element ) then

!     add only interface elements

      call eltree_to_interface_mesh ( eltree, elmesh, mapcoor=mapcoor )

      indicator = 2 ! replace element with submesh

      call delete(eltree)

    else

!     fully in or out: remove element

      indicator = 1

    end if

  end subroutine userelmesh1


! define the mesh for each element (volume)

  subroutine userelmesh2 ( mesh, elgrp, elem, indicator, elmesh )

    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: elgrp, elem
    integer, intent(out) :: indicator
    type(mesh_t), intent(inout) :: elmesh

    logical :: interface_in_element
    integer :: nod(mesh%element(elgrp)%numnod)
    real(dp) :: coor(2,2)


    nod = mesh%topology(elgrp)%a(:,elem)

!   interface in element?

    if ( all ( d(nod) > spl_th ) .or. all ( d(nod) < -spl_th ) ) then

!     all nodes far from the interface

      interface_in_element = .false.

    else

!     element possibly contains an interface, start subdivide

!     fill root node of the eltree

      coor(1,:) = -1._dp  ! lower corner of the reference region (-1,-1,-1)
      coor(2,:) =  1._dp  ! upper corner of the reference region (1,1,1)

      call fill_node_eltree ( eltree, coor )

!     divide root reference domain into subdomains

      lelgrp = elgrp
      lelem = elem

      call subdivide ( eltree, levelset=levelset, numsplit=ns_plot, &
        split_threshold=spl_th, numsplitmin=nsmin, mapcoor=mapcoor, submesh=.true., &
        intmesh=.true. )

      if ( number_of_interface_elements(eltree) > 0 ) then

        interface_in_element = .true.

      else

        interface_in_element = .false.

        call delete(eltree)

      end if

    end if

!   add elements

    if ( interface_in_element ) then

!     add sub elements

      call eltree_to_mesh ( eltree, elmesh, mapcoor=mapcoor )

      indicator = 2 ! replace element with submesh

      call delete(eltree)

    else

!     fully in or out: leave element

      indicator = 0

    end if

  end subroutine userelmesh2

end module subs_sd5_m

program eltree5

  use tfem_m
  use eltree_m
  use functions_sd5_m
  use subs_sd5_m
  use figplot_m
  use io_utils_m
  use limits_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    nx=5,               & ! number of elements in x
    ny=5,               & ! number of elements in y
    numsplit = 4,       & ! relative size of smallest elements near
                          ! interface 1/2^numsplit
    numsplitmin = 1,    & ! relative size of all subelements will be
                          ! at least as small as 1/2^numsplitmin
    numsplit_plot = 4,  & ! parameter in subdivide for the plot mesh
    gausse = 2            ! integration rule on one interface element of eltree

  real(dp), parameter :: &
    radius = 1._dp,    & ! radius of the cylinder
    center(2) = [ 0.0_dp, 0.0_dp ], & ! center position of the cylinder
    split_threshold = 1.e-2_dp    ! levelset value for being "close" enough
                                  ! to the interface for tree splitting

  logical, parameter :: &
    curved= .true.   ! curved mesh domain


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t), target :: mesh
  type(mesh_t) :: mesh_plot, mesh1, mesh2
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients

  integer :: i
  real(dp) :: resultsum(4), areae, intge

! module functions_sd5_m

  rpl = radius ! radius of the circular object
  cpl = center ! initial position of the center of the object

! module stokes_usergauss_m

  intrule_ie = gausse  ! Gauss integration for the eltree subelements
  ns = numsplit ! parameter in subdivide for the integration scheme
  nsmin = numsplitmin ! parameter in subdivide
  ns_plot = numsplit_plot ! parameter in subdivide for the plot mesh
  spl_th = split_threshold

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl

  coefficients%r(1) = 0


! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  if ( curved ) then

    meshgen_options%regionshape = 3
    meshgen_options%x2d = &
      reshape ( [ -2.0_dp,  2.0_dp, 2.0_dp, -2.0_dp,    &
                   -2.0_dp, -2.0_dp, 2.0_dp,  2.0_dp ], &
                [4,2] )
    meshgen_options%curved(1:4) = .true.
    meshgen_options%funcnr(1:4) = [ 1, 2, 3, 4 ]

    call quadrilateral2d ( mesh, meshgen_options, func=cfunc )

  else

    meshgen_options%ox = -2
    meshgen_options%oy = -2
    meshgen_options%lx = 4
    meshgen_options%ly = 4

    call quadrilateral2d ( mesh, meshgen_options )

  end if

  call fill_mesh_parts ( mesh )

  lmesh => mesh  ! supply mesh to mapcoor for mapping coordinates


! problem definition

  call create_input_probdef ( mesh, input_probdef )

  input_probdef%elementdof(1)%a = 1

  call problem_definition ( input_probdef, mesh, problem )

! levelset in all nodes

  allocate ( d(mesh%nnodes) )

  d = levelset ( mesh%coor ) ! set levelset for all nodes

! integrate on boundary of circle

  call integrate ( mesh, problem, resultsum, elemsub=integrate_interface, &
    coefficients=coefficients )

  areae = 2*pi*radius
  intge = pi*radius**3

  print *, 'computed length of circle boundary = ', resultsum(1), &
    'relative error = ', (resultsum(1)-areae)/areae
  print *, 'computed integral of x^2 on circle boundary = ', resultsum(2), &
    'relative error = ', (resultsum(2)-intge)/intge
  print *, 'computed normal integral circle = ', resultsum(3:4)


! post-processing: create mesh_plot for plotting

  call mesh_convert ( mesh, mesh1, userelmesh=userelmesh1, warn=.false. )
  !call fill_mesh_parts ( mesh1 )
  !call plot_mesh ( plot_options, mesh1, filename='mesh1.fig' )
  call mesh_convert ( mesh, mesh2, userelmesh=userelmesh2, warn=.false. )
  !call fill_mesh_parts ( mesh2 )
  !call plot_mesh ( plot_options, mesh2, filename='mesh2.fig' )
  call mesh_merge ( mesh1, mesh2, mesh_plot, warn=.false. )
  !call printinfo ( mesh1 )
  !call printinfo ( mesh2 )
  call delete ( mesh1, mesh2 )

  call fill_mesh_parts ( mesh_plot )
  !call plot_mesh ( plot_options, mesh_plot, filename='mesh_plot.fig' )

  !call printinfo ( mesh_plot)

  plot_options%boundarycolor=0
  plot_options%plotboundary=.false.

  call plot_mesh ( plot_options, mesh_plot, 'mesh_sd5_plot.fig', groups=[(i,i=3,mesh_plot%nelgrp)] )
  plot_options%meshcolor=4
  call plot_mesh ( plot_options, mesh_plot, 'mesh_sd5_plot.fig', groups=[2], append=.true. )

  call write_mesh_vtk ( mesh_plot, groups=[2], filename='mesh.vtk' )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh, mesh_plot )
  call delete ( coefficients )

  deallocate ( d )

end program eltree5
