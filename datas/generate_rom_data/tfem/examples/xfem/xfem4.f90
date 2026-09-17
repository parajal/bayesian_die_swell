! 2D Stokes problem on a square domain with a cylindrical shaped drop
! eXtended FEM approach with weak constraint on the cylinder for
! a coupling of velocities on the interface between the two fluids.
! Quadtree+triangle subdivision for integration.
! Fluids described by layers.
! Plot a layer.

! functions module for xfem4

module functions_xf4_m

  use tfem_elem_m

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

end module functions_xf4_m


! module for the user gauss routines

module stokes_usergauss_m

  use tfem_elem_m
  use functions_xf4_m
  use eltree_m

  implicit none

  save

! user gauss is based on a quadtree using the add-on eltree
  type(eltree_t) :: eltree

! the levelset function in all the nodes (it is a signed distance function
! in this case of a simple cylinder).
  real(dp), dimension(:), allocatable :: d

! the standard integration scheme:
!   ninti_s: number of integration points
!   xs, ws: points and weights
! the integration scheme on subelements in the eltree:
!   intrule_e: the integration rule
!   ninti_e: number of integration points
!   xe, we: points and weights
! the integration scheme on subelements in the submesh (triangles):
!   intrule_sub: the integration rule
!   ninti_sub: number of integration points
!   xsub, wsub: points and weights
  integer :: ninti_s, intrule_e, ninti_e, intrule_sub, ninti_sub
  real(dp), allocatable :: ws(:), xs(:,:), we(:), xe(:,:), wsub(:), xsub(:,:)

! integration options for setting the composite Gauss integration on the eltree
  type(integration_options_t) :: intopt

! numsplit parameters
!  ns : numsplit parameter in subdivide for the integration scheme
!  nsmin : numsplitmin parameter in subdivide
!  ns_plot : numsplit parameter in subdivide for the plot mesh
  integer :: ns = 2, nsmin = 1, ns_plot = 2

!  epsjac: subelements in smesh (submesh) having a jacobian smaller than
!          epsjac will be omitted.
  real(dp) :: epsjac = 0._dp

! some logicals
  logical :: warn = .false., printint = .false.

! how to integrate interface elements:
!  intregion = 1 : outside region
!  intregion = 2 : inside region

  integer :: intregion = 1

! how to generate mesh:
!  meshregion = 1 : outside region
!  meshregion = 2 : inside region

  integer :: meshregion = 1

contains


! user routine for setting the number of integration points

  subroutine set_ninti_user ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors

    integer :: nod(nodalp)
    real(dp) :: coor(2,2)
    type(gauss_t) :: gausse, gausssub

    if ( first ) then

!     first element in this group

!     Gauss rule for the standard case

      ninti_s = ninti

      allocate ( ws(ninti_s), xs(ninti_s,2) )

      call set_Gauss_integration ( gauss, xs, ws )

!     Gauss rule for the subdivided elements in the eltree

      gausse = gauss
      gausse%intrule = intrule_e

      call set_ninti ( gausse, ninti_e )

      allocate ( we(ninti_e), xe(ninti_e,2) )

      call set_Gauss_integration ( gausse, xe, we )

!     Gauss rule for the subdivided elements in the submesh

      gausssub%globalshape = 'triangle'
      gausssub%intrule = intrule_sub

      call set_ninti ( gausssub, ninti_sub )

      allocate ( wsub(ninti_sub), xsub(ninti_sub,2) )

      call set_Gauss_integration ( gausssub, xsub, wsub )

    end if

!   Set ninti

    nod = mesh%topology(elgrp)%a(:,elem)

    if ( .not. ( all ( d(nod) <= 0._dp ) .or. all ( d(nod) >= 0._dp ) ) ) then

!     element crosses the interface

!     fill root node of the eltree

      coor(1,:) = -1._dp  ! lower corner of the reference region (-1,-1,-1)
      coor(2,:) =  1._dp  ! upper corner of the reference region (1,1,1)

      call fill_node_eltree ( eltree, coor )

!     divide root reference domain into subdomains

      lelgrp = elgrp
      lelem = elem

      call subdivide ( eltree, levelset=levelset, numsplit=ns, &
        numsplitmin=nsmin, mapcoor=mapcoor, submesh=.true. )

!     composite integration scheme for levelset >= 0 or levelset <=0

      intopt%epsjac = epsjac

      if ( intregion == 1 ) then
!       outside region
        ninti = number_of_integration_points ( eltree, lsign=[0,1],&
          ninti=ninti_e, nintis=ninti_sub, integration_options=intopt )
      else if ( intregion == 2 ) then
!       inside region
        ninti = number_of_integration_points ( eltree, lsign=[-1,0],&
          ninti=ninti_e, nintis=ninti_sub, integration_options=intopt )
      end if

      if ( printint ) then
        print *, 'number of integration points = ', ninti
      end if

    else

!     fully in or out

      ninti = ninti_s

    end if

  end subroutine set_ninti_user


! set the composite gauss integration rule with a user subroutine

  subroutine set_Gauss_integration_user ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors

    integer :: nod(nodalp)

    nod = mesh%topology(elgrp)%a(:,elem)

    if ( .not. ( all ( d(nod) <= 0._dp ) .or. all ( d(nod) >= 0._dp ) ) ) then

!     element crosses the interface

!     composite integration scheme for levelset >= 0 or levelset <= 0

      if ( intregion == 1 ) then
!       outside region
        call integration_points ( eltree, xig, wg, lsign=[1], ninti=ninti_e, &
          xe=xe, we=we, nintis=ninti_sub, xs=xsub, ws=wsub, &
          integration_options=intopt )
      else if ( intregion == 2 ) then
!       inside region
        call integration_points ( eltree, xig, wg, lsign=[-1], ninti=ninti_e,&
          xe=xe, we=we, nintis=ninti_sub, xs=xsub, ws=wsub, &
          integration_options=intopt )
      end if

      call delete(eltree)

    else

!     fully out: standard integration

      call set_Gauss_integration ( gauss, xig, wg )

    end if

!   delete some allocated memory

    if ( last ) then

      deallocate ( we, xe, ws, xs, wsub, xsub )

    end if

  end subroutine set_Gauss_integration_user


! define the mesh for each element

  subroutine userelmesh ( mesh, elgrp, elem, indicator, elmesh )

    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: elgrp, elem
    integer, intent(out) :: indicator
    type(mesh_t), intent(inout) :: elmesh

    integer :: nod(mesh%element(elgrp)%numnod)
    real(dp) :: coor(2,2)


    nod = mesh%topology(elgrp)%a(:,elem)

    if ( .not. ( all ( d(nod) <= 0._dp ) .or. all ( d(nod) >= 0._dp ) ) ) then

!     element crosses the interface

!     fill root node of the eltree

      coor(1,:) = -1._dp  ! lower corner of the reference region (-1,-1,-1)
      coor(2,:) =  1._dp  ! upper corner of the reference region (1,1,1)

      call fill_node_eltree ( eltree, coor )

!     divide root reference domain into subdomains

      lelgrp = elgrp
      lelem = elem

      call subdivide ( eltree, levelset=levelset, numsplit=ns_plot, &
        numsplitmin=nsmin, mapcoor=mapcoor, submesh=.true. )

!     add only elements in the outside or inside region

      if ( meshregion == 1 ) then
!       outside region
        call eltree_to_mesh ( eltree, elmesh, lsign=[1], mapcoor=mapcoor )
      else if ( meshregion == 2 ) then
!       inside region
        call eltree_to_mesh ( eltree, elmesh, lsign=[-1], mapcoor=mapcoor )
      end if

      indicator = 2 ! replace element with submesh

      call delete(eltree)

    else if ( all ( d(nod) <= 0._dp ) ) then

!     fully in

      if ( meshregion == 1 ) then
!       outside region
        indicator = 1  ! remove element
      else if ( meshregion == 2 ) then
!       inside region
        indicator = 0  ! leave element as it is
      end if

    else

!     fully out

      if ( meshregion == 1 ) then
!       outside region
        indicator = 0  ! leave element as it is
      else if ( meshregion == 2 ) then
!       inside region
        indicator = 1  ! remove element
      end if

    end if

  end subroutine userelmesh

end module stokes_usergauss_m


! the actual program

program xfem4

  use tfem_m
  use stokes_elements_m
  use hsl_ma57_m
  use io_utils_m
  use figplot_m
  use subs4_m
  use functions_xf4_m
  use stokes_usergauss_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    nx=80,              & ! number of elements in x
    ny=80,              & ! number of elements in y
    numsplit = 3,       & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit
    numsplitmin = 1,    & ! relative size of all subelements will be at least
                          ! as small as 1/2^numsplitmin
    numsplit_plot = 3,  & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit for plotting
    nsubinto = 30,      & ! number of subintegration points on object
    gausse = 3,         & ! Gauss integration for the eltree subelements
    gausssub = 6,       & ! Gauss integration for the submesh subelements
    gauss = 3             ! 3x3 Gauss integration

  real(dp), parameter :: &
    lx = 8._dp,        & ! width of domain
    ly = 8._dp,        & ! height of domain
    radius = 1._dp,    & ! radius of the cylinder
    center(2) = [ 0.045_dp, 1.010_dp ], & ! center position of the cylinder
    eta1 = 1._dp,      & ! viscosity outside fluid
    eta2 = 10._dp,     & ! viscosity inside fluid
    flowrate = 10._dp, & ! flowrate
    epsdef = 0.01_dp    ! deformation of the mesh

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t), target :: mesh
  type(mesh_t) :: mesh_particle, mesh_plot1, mesh_plot2, mesh_plot
  type(input_probdef_t) :: input_probdef, input_probdef_plot
  type(problem_t) :: problem, problem_plot
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(sample_t) :: sample_p, sample_v, sample_vort
  type(solver_options_ma57_t) :: solver_options

  integer :: i, elem, nelem, nnodes, plotobj1, plotobj2, nbx, nby

  integer, dimension(:), allocatable :: nodes

! plot nodesets
  logical, parameter :: plot_nodesets = .true.

! plot data in layers
  logical, parameter :: plot_data_in_layers = .false.


! set some parameters in modules

! module timer_m

  timer = .false. ! set to .true. to show cpu time output

! module functions1_m

  rpl = radius ! radius of the circular object
  cpl = center ! initial position of the center of the object

! module stokes_usergauss_m

  intrule_e = gausse  ! Gauss integration for the eltree subelements
  intrule_sub = gausssub  ! Gauss integration for the eltree subelements
  ns = numsplit ! parameter in subdivide for the integration scheme
  nsmin = numsplitmin ! parameter in subdivide
  ns_plot = numsplit_plot ! parameter in subdivide for the plot mesh
  epsjac = 0._dp ! subelements in smesh (submesh) having a jacobian smaller
                 ! than epsjac will be omitted.

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0
  coefficients%i(37) = 2 ! Gauss points defined by user subroutines

  coefficients%r(1) = eta1
  coefficients%r(2:) = 0
  coefficients%r(6) = flowrate

  coefficients%set_ninti_user => set_ninti_user
  coefficients%set_Gauss_integration_user => set_Gauss_integration_user

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%lx = lx
  meshgen_options%ly = ly
  meshgen_options%ox = -lx/2
  meshgen_options%oy = -ly/2

  call quadrilateral2d ( mesh, meshgen_options )

! deform the mesh a little bit to avoid zero or very small integration areas

  mesh%coor(:,1) = mesh%coor(:,1) * ( 1 - epsdef*(1-2*abs(mesh%coor(:,1))/lx) )
  mesh%coor(:,2) = mesh%coor(:,2) * ( 1 - epsdef*(1-2*abs(mesh%coor(:,2))/ly) )

  call add_to_mesh ( mesh, curve=[-4] ) ! curve 5

! create mesh for particle boundary

  nelem = 4 * nx * radius / lx

  call mesh_skeleton ( mesh_particle, nnodes=2*nelem, nelem=nelem, elshape=2, &
    ndim=2 )

  rp = radius ! radius of the circular object
  xpc = center ! initial position of the center of the object

  call objectscoor ( 1, mesh_particle%coor )

  do elem = 1, mesh_particle%nelem
    mesh_particle%topology(1)%a(:,elem) = [ 2*elem-1, 2*elem, 2*elem + 1 ]
  end do
  mesh_particle%topology(1)%a(3,mesh_particle%nelem) = 1 ! close circle

! one object

  call add_to_mesh ( mesh, object='mesh', objectmesh=mesh_particle, &
    topology=.true., intrule=2, nsubint=nsubinto )

! set levelset in the nodes

  allocate ( nodes(mesh%nnodes), d(mesh%nnodes) )

  d = levelset ( mesh%coor ) ! set levelset for all nodes

! outside nodeset for layer 1

  call find_layer_nodes ( 'out' )

  call add_to_mesh ( mesh, nodeset='nodes', nodes=nodes(1:nnodes) )

! object for plotting the nodeset

  call add_to_mesh ( mesh, object='coordinates', &
    coor=mesh%coor(nodes(1:nnodes),:) )

! inside nodeset for layer 2

  call find_layer_nodes ( 'in' )

  call add_to_mesh ( mesh, nodeset='nodes', nodes=nodes(1:nnodes) )

! object for plotting the nodeset

  call add_to_mesh ( mesh, object='coordinates', &
    coor=mesh%coor(nodes(1:nnodes),:) )

! blocks for fast searching of mesh intersections

  nbx = nint(sqrt(real(nx)))
  nby = nint(sqrt(real(ny)))

  call add_to_mesh ( mesh, blocks=[nbx,nby] )

  call fill_mesh_parts ( mesh )

  lmesh => mesh  ! supply mesh to mapcoor for mapping coordinates

  call plot_points_curves ( plot_options, mesh, 'curves_xf4.fig' )

  if ( plot_nodesets ) then

    call plot_mesh_nodeset ( nodeset=1, filename='mesh_outside_xf4.fig' )
    call plot_mesh_nodeset ( nodeset=2, filename='mesh_inside_xf4.fig' )

  else

    call plot_mesh ( plot_options, mesh, 'mesh_xf4.fig' )

  end if

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2, &
    numlayers=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]
  input_probdef%layers = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, physq=1, layer=1 )
  call define_essential ( mesh, input_probdef, curve1=3, physq=1, layer=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2, layer=1 )

! define constraint on the object for coupling velocities between the layers

  call define_constraint ( mesh, input_probdef, object=1, layer1=1, object2=1,&
    layer2=2, physq=1, discretization='weak', elementdof=[2,0,2] )

! define constraint on curve for flowrate

  call define_constraint ( mesh, input_probdef, &
    physq=1, layer=1, curve1=2, nglobalc=1 )

! define constraint on curve for periodical boundary conditions

  call define_constraint ( mesh, input_probdef, &
    physq=1, layer=1, curve1=2, curve2=5, discretization='collocation', &
    exclude=3 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

  sol%u = 0

! create the structure oldvectors

  call create ( oldvectors, nsysvec=1 )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

  call tic

! build (assemble) matrix and vector from elements

! outside fluid

  intregion = 1

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, oldvectors=oldvectors, &
    coefficients=coefficients, buildvector=.false., layer=1 )

! inside fluid

  intregion = 2
  coefficients%r(1) = eta2

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, oldvectors=oldvectors, addmat=.true., &
    coefficients=coefficients, buildvector=.false., layer=2 )

  call toc ( 'sys' )

! build missing matrix elements with zero

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    coefficients=coefficients, zeromatvec=.true., &
    addmat=.true., exclude_single_layer =.true. )

  call toc ( 'syszero' )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=elementc, addmatvec=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=3, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficients  )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call check ( sysmatrix )

  call toc ( 'build' )

  solver_options%real_storage=4._dp

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

  call toc ( 'solve' )

! post-processing

  oldvectors%s(1)%p => sol

! sampling for gnuplot

! write profile to a file for plotting with gnuplot

  call printtofile ( mesh, problem, 'vprofile_xf4.out', curve=2, sysvector=sol )

! sample pressure on cylinder

  call sample_pressure_on_cylinder ( layer=1, filename='sample_p_outside_xf4.out' )
  call sample_pressure_on_cylinder ( layer=2, filename='sample_p_inside_xf4.out' )

! sample velocity on cylinder

  call sample_velocity_on_cylinder ( layer=1, filename='sample_v_outside_xf4.out' )
  call sample_velocity_on_cylinder ( layer=2, filename='sample_v_inside_xf4.out' )

! sample vorticity on cylinder

  call sample_vorticity_on_cylinder ( layer=1, filename='sample_vort_outside_xf4.out' )
  call sample_vorticity_on_cylinder ( layer=2, filename='sample_vort_inside_xf4.out' )


! post-processing: create mesh_plot for plotting

  meshregion = 1
  call mesh_convert ( mesh, mesh_plot1, userelmesh=userelmesh, warn=.false. )
  meshregion = 2
  call mesh_convert ( mesh, mesh_plot2, userelmesh=userelmesh, warn=.false. )

  call mesh_merge ( mesh_plot1, mesh_plot2, mesh_plot )

  call fill_mesh_parts ( mesh_plot )

  call plot_mesh ( plot_options, mesh_plot, 'mesh_plot_xf4.fig' )

  call toc ( 'plot_mesh' )


! create objects for interpolation on the original mesh

  warn_add_to_mesh_after_meshgen_parts = .false. ! .false. to suppress warning
  call add_to_mesh ( mesh, object='coordinates', coor=mesh_plot1%coor )
  plotobj1 = mesh%nobjects
  call add_to_mesh ( mesh, object='coordinates', coor=mesh_plot2%coor )
  plotobj2 = mesh%nobjects
  call fill_mesh_parts_objects ( mesh, object1=plotobj1, object2=plotobj2 )

  call delete ( mesh_plot1, mesh_plot2 )

  call toc ( 'fill_mesh_parts_objects' )


! problem definition for post_processing

  call create_input_probdef ( mesh_plot, input_probdef_plot, nvec=2 )

  do i = 1, mesh_plot%nelgrp
    input_probdef_plot%elementdof(i)%a(:) = 0

    input_probdef_plot%vec_elementdof(i)%a(:,1) = 2  ! velocity
    input_probdef_plot%vec_elementdof(i)%a(:,2) = 1  ! scalar
  end do

  call problem_definition ( input_probdef_plot, mesh_plot, problem_plot )

  call create_vector ( problem_plot, velocity, vec=1 )
  call create_vector ( problem_plot, pressure, vec=2 )
  call create_vector ( problem_plot, vorticity, vec=2 )

! sample and plot pressure on mesh_plot

  call sample_pressure
  call plot_contour_on_mesh_plot ( filename='pressure_xf4.fig', vector=pressure )

! sample and plot vorticity on mesh_plot

  call sample_vorticity
  call plot_contour_on_mesh_plot ( filename='vorticity_xf4.fig', vector=vorticity )

! sample and plot velocity on mesh_plot

  call sample_velocity
  call plot_contour_on_mesh_plot ( filename='velocity_x_xf4.fig', &
    vector=velocity, degfd=1 )
  call plot_contour_on_mesh_plot ( filename='velocity_y_xf4.fig', &
    vector=velocity, degfd=2 )

  call toc ( 'fill_sample' )


! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=pressure, &
    dataname='pressure', filename='drop_in_flow_xf4.vtk' )

  call write_vector_vtk ( mesh_plot, problem_plot, &
    filename='drop_in_flow_xf4.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true. )

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=vorticity, &
    filename='drop_in_flow_xf4.vtk', &
    dataname='vorticity', append=.true. )

  call delete ( pressure, vorticity )


  if ( plot_data_in_layers ) then

!   post-processing: plot in a layer

    call create_vector ( problem, pressure, vec=3 )
    call create_vector ( problem, vorticity, vec=3 )

    call plot_in_layer ( 'velocity_xf4_layer1.fig', &
      'pressure_contour_xf4_layer1.fig', &
      'vorticity_contour_xf4_layer1.fig', layer=1 )
    call plot_in_layer ( 'velocity_xf4_layer2.fig', &
      'pressure_contour_xf4_layer2.fig', &
      'vorticity_contour_xf4_layer2.fig', layer=2 )

    call delete ( pressure, vorticity )

  end if


! delete all data including all allocated memory

  call delete ( problem, problem_plot )
  call delete ( input_probdef, input_probdef_plot )
  call delete ( mesh, mesh_plot, mesh_particle )
  call delete ( sol, rhsd )
  call delete ( velocity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( sample_p, sample_vort, sample_v )

  deallocate ( nodes, d )

contains


! find inside and outside layer nodes

  subroutine find_layer_nodes ( region )

    character(len=*), intent(in) :: region

    integer, dimension(mesh%nnodes) :: nodes1
    integer :: elem, i

    if ( region == 'in' ) then

!     set internal nodes

      where ( d <= 0._dp )
        nodes = 1
      else where
        nodes = 0
      end where

    else if ( region == 'out' ) then

!     set external nodes

      where ( d > 0._dp )
        nodes = 1
      else where
        nodes = 0
      end where

    end if

    nodes1 = nodes

!   add nodes connected to elements at the interface

    do elem = 1, mesh%grpnumel(1)
       if ( .not. all ( nodes(mesh%topology(1)%a(:,elem)) == &
                        nodes(mesh%topology(1)%a(1,elem)) ) ) then
         nodes1(mesh%topology(1)%a(:,elem)) = 1
       end if
    end do

!   count and collect the nodes

    nnodes = 0
    do i = 1, mesh%nnodes
      if ( nodes1(i) == 0 ) cycle
      nnodes = nnodes + 1
      nodes(nnodes) = i
    end do

  end subroutine find_layer_nodes


! plot mesh including the specified nodeset

  subroutine plot_mesh_nodeset ( nodeset, filename )

    integer, intent(in) :: nodeset
    character(len=*), intent(in) :: filename

    call plot_mesh ( plot_options, mesh, filename )

    plot_options%objectpointcolor=4
    plot_options%objectpointsize=0.4

    call plot_objects ( plot_options, mesh, filename, &
      append=.true., object1=1 )

    plot_options%objectpointcolor=1
    plot_options%objectpointsize=0.1

    call plot_objects ( plot_options, mesh, filename, &
      append=.true., object1=nodeset+1 )

  end subroutine plot_mesh_nodeset


! sample pressure on cylinder

  subroutine sample_pressure_on_cylinder ( layer, filename )

    integer, intent(in) :: layer
    character(len=*), intent(in) :: filename

    integer :: i

    open( unit=15, file=filename, recl=300 )

    coefficients%i(38) = layer

!   sample on cylinder (errorobject=1 to warn for points not in layer)

    call fill_sample ( mesh, problem, sample_p, ndegfd=1, object=1, &
      layer=layer, elemsub=stokes_sample_pressure, coefficients=coefficients, &
      oldvectors=oldvectors, errorobject=1 )

    do i = 1, sample_p%nnodes
      write( unit=15, fmt=* ) sample_p%coor(i,:), sample_p%u(i,:)
    end do

    close( unit=15 )

  end subroutine sample_pressure_on_cylinder


! sample velocity on cylinder

  subroutine sample_velocity_on_cylinder ( layer, filename )

    integer, intent(in) :: layer
    character(len=*), intent(in) :: filename

    integer :: i

    open( unit=15, file=filename, recl=300 )

    coefficients%i(38) = layer

!   sample on cylinder (errorobject=1 to warn for points not in layer)

    call fill_sample ( mesh, problem, sample_v, ndegfd=2, object=1, &
      layer=layer, elemsub=stokes_sample_velocity, coefficients=coefficients, &
      oldvectors=oldvectors, errorobject=1 )

    do i = 1, sample_v%nnodes
      write( unit=15, fmt=* ) sample_v%coor(i,:), sample_v%u(i,:)
    end do

    close( unit=15 )

  end subroutine sample_velocity_on_cylinder


! sample vorticity on cylinder

  subroutine sample_vorticity_on_cylinder ( layer, filename )

    integer, intent(in) :: layer
    character(len=*), intent(in) :: filename

    integer :: i

    open( unit=15, file=filename, recl=300 )

    coefficients%i(38) = layer

    coefficients%i(13)=5

!   sample on cylinder (errorobject=1 to warn for points not in layer)

    call fill_sample ( mesh, problem, sample_vort, ndegfd=1, object=1, &
      layer=layer, elemsub=stokes_sample_deriv, coefficients=coefficients, &
      oldvectors=oldvectors, errorobject=1 )

    do i = 1, sample_vort%nnodes
      write( unit=15, fmt=* ) sample_vort%coor(i,:), sample_vort%u(i,:)
    end do

    close( unit=15 )

  end subroutine sample_vorticity_on_cylinder


! sample pressure on mesh_plot

  subroutine sample_pressure

    integer :: n

!   sample outside (use layer as argument to detect points not in layer)

    coefficients%i(38) = 1

    call fill_sample ( mesh, problem, sample_p, ndegfd=1, object=plotobj1, &
      elemsub=stokes_sample_pressure, coefficients=coefficients, &
      oldvectors=oldvectors, layer=1 )

    n = sample_p%nnodes

    pressure%u(1:n) = sample_p%u(1:n,1)

!   sample inside (use layer as argument to detect points not in layer)

    coefficients%i(38) = 2

    call fill_sample ( mesh, problem, sample_p, ndegfd=1, object=plotobj2, &
      elemsub=stokes_sample_pressure, coefficients=coefficients, &
      oldvectors=oldvectors, layer=2 )

    pressure%u(n+1:) = sample_p%u(1:sample_p%nnodes,1)

  end subroutine sample_pressure


! sample vorticity on mesh_plot

  subroutine sample_vorticity

    integer :: n

!   sample outside (use layer as argument to detect points not in layer)

    coefficients%i(38) = 1

    coefficients%i(13)=5

    call fill_sample ( mesh, problem, sample_vort, ndegfd=1, object=plotobj1, &
      elemsub=stokes_sample_deriv, coefficients=coefficients, &
      oldvectors=oldvectors, layer=1 )

    n = sample_vort%nnodes

!   sample inside (use layer as argument to detect points not in layer)

    vorticity%u(1:n) = sample_vort%u(1:n,1)

    coefficients%i(38) = 2

    call fill_sample ( mesh, problem, sample_vort, ndegfd=1, object=plotobj2, &
      elemsub=stokes_sample_deriv, coefficients=coefficients, &
      oldvectors=oldvectors, layer=2 )

    vorticity%u(n+1:) = sample_vort%u(1:sample_vort%nnodes,1)

  end subroutine sample_vorticity


! sample velocity on mesh_plot

  subroutine sample_velocity

    integer :: n

!   sample outside (use layer as argument to detect points not in layer)

    coefficients%i(38) = 1

    call fill_sample ( mesh, problem, sample_v, ndegfd=2, object=plotobj1, &
      elemsub=stokes_sample_velocity, coefficients=coefficients, &
      oldvectors=oldvectors, layer=1 )

    n = 2*sample_v%nnodes

    velocity%u(1:n) = reshape ( transpose(sample_v%u), [n] )

!   sample inside (use layer as argument to detect points not in layer)

    coefficients%i(38) = 2

    call fill_sample ( mesh, problem, sample_v, ndegfd=2, object=plotobj2, &
      elemsub=stokes_sample_velocity, coefficients=coefficients, &
      oldvectors=oldvectors, layer=2 )

    velocity%u(n+1:) = reshape ( transpose(sample_v%u), [2*sample_v%nnodes] )

  end subroutine sample_velocity


! plot contour on mesh_plot

  subroutine plot_contour_on_mesh_plot ( filename, vector, degfd )

    character(len=*), intent(in) :: filename
    type(vector_t), intent(in) :: vector
    integer, intent(in), optional :: degfd

    plot_options%plotboundary2=.false.

    call plot_color_contour ( plot_options, mesh_plot, problem_plot, &
      filename=filename, vector=vector, degfd=degfd )

    call plot_boundary ( plot_options, mesh, filename=filename, append=.true. )

    plot_options%objectpointcolor=4
    plot_options%objectpointsize=0.1

    call plot_objects ( plot_options, mesh, filename=filename, &
      append=.true., object1=1 )

  end subroutine plot_contour_on_mesh_plot


! plot data in a layer

  subroutine plot_in_layer ( filename1, filename2, filename3, layer )

    character(len=*), intent(in) :: filename1, filename2, filename3
!   layer for postprocessing: 1 is outside, 2 is inside
    integer, intent(in) :: layer

    coefficients%i(38) = layer

    call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
      coefficients=coefficients, oldvectors=oldvectors, layer=layer )

    coefficients%i(13)=5

    call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
      coefficients=coefficients, oldvectors=oldvectors, layer=layer )

!   write to a fig file for plotting

    call plot_vector ( plot_options, mesh, problem, filename1, &
      sysvector=sol, physq=1, layer=layer )
    call plot_mesh ( plot_options, mesh, filename1, append=.true. )
    plot_options%objectpointcolor=4
    plot_options%objectpointsize=0.1
    call plot_objects ( plot_options, mesh, filename1, append=.true., object1=1 )
    plot_options%objectpointcolor=1
    plot_options%objectpointsize=0.1
    call plot_objects ( plot_options, mesh, filename1, append=.true., object1=layer+1 )

    call plot_color_contour ( plot_options, mesh, problem, &
      filename2, vector=pressure, layer=layer )
    call plot_mesh ( plot_options, mesh, filename2, append=.true. )
    plot_options%objectpointcolor=4
    plot_options%objectpointsize=0.1
    call plot_objects ( plot_options, mesh, filename2, append=.true., object1=1 )
    plot_options%objectpointcolor=1
    plot_options%objectpointsize=0.1
    call plot_objects ( plot_options, mesh, filename2, append=.true., object1=layer+1 )

    call plot_color_contour ( plot_options, mesh, problem, &
      filename3, vector=vorticity, layer=layer )
    call plot_mesh ( plot_options, mesh, filename3, append=.true. )
    plot_options%objectpointcolor=4
    plot_options%objectpointsize=0.1
    call plot_objects ( plot_options, mesh, filename3, append=.true., object1=1 )
    plot_options%objectpointcolor=1
    plot_options%objectpointsize=0.1
    call plot_objects ( plot_options, mesh, filename3, append=.true., object1=layer+1 )

  end subroutine plot_in_layer

end program xfem4
