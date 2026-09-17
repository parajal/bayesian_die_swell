! 2D Stokes problem on a square domain with a cylinder
! eXtended FEM approach with weak constraint on the cylinder.
! Quadtree+triangle subdivision for integration.

! functions module for xfem2

module functions_xf2_m

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

end module functions_xf2_m


! module for the user gauss routines

module stokes_usergauss_m

  use tfem_elem_m
  use functions_xf2_m
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
  logical :: warn = .false., printnumel = .false., printint = .false.

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

      if ( warn .and. number_of_subelements(eltree,lsign=[1]) <= 1 ) then
        write(*,'(a)') 'Minimal number subelements in outer region'
      end if

      if ( printnumel ) then
        print *, 'number of subelements outside = ', &
          number_of_subelements(eltree,lsign=[1])
      end if

!     composite integration scheme for levelset >= 0

      intopt%epsjac = epsjac

      ninti = number_of_integration_points ( eltree, lsign=[0,1],&
        ninti=ninti_e, nintis=ninti_sub, integration_options=intopt )

      if ( printint ) then
        print *, 'number of integration points = ', ninti
      end if

    else if ( all ( d(nod) <= 0._dp ) ) then

!     fully in

      ninti = 0

    else

!     fully out

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

!     composite integration scheme for levelset >= 0

      call integration_points ( eltree, xig, wg, lsign=[1], ninti=ninti_e, &
        xe=xe, we=we, nintis=ninti_sub, xs=xsub, ws=wsub, &
        integration_options=intopt )

      call delete(eltree)

    else if ( all ( d(nod) <= 0._dp ) ) then

!     fully in: generate fully zero elemmat and elemvec

      xig = 0; wg = 0

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

!     add only elements in the outside region

      call eltree_to_mesh ( eltree, elmesh, lsign=[1], mapcoor=mapcoor )

      indicator = 2 ! replace element with submesh

      call delete(eltree)

    else if ( all ( d(nod) <= 0._dp ) ) then

!     fully in: remove element

      indicator = 1

    else

!     fully out: leave element as it is.

      indicator = 0

    end if

  end subroutine userelmesh

end module stokes_usergauss_m


! the actual program

program xfem2

  use tfem_m
  use stokes_elements_m
  use hsl_ma57_m
  use io_utils_m
  use figplot_m
  use subs1_m
  use functions_xf2_m
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
    center(2) = [ 0._dp, 0._dp ], & ! center position of the cylinder
    eta = 1._dp,       & ! viscosity
    flowrate = 10._dp, & ! flowrate
    epsdef = 0.01_dp    ! deformation of the mesh

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t), target :: mesh
  type(mesh_t) :: mesh_particle, mesh_plot
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

  integer :: i, elem, nelem, nnodes, plotobj, nbx, nby

  integer, dimension(:), allocatable :: nodes, nodes1

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

  coefficients%r(1) = eta
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

! nodeset for Dirichlet on internal nodes

  allocate ( nodes(mesh%nnodes), d(mesh%nnodes) )

  d = levelset ( mesh%coor ) ! set levelset for all nodes

  call find_internal_zero_nodes

  call add_to_mesh ( mesh, nodeset='nodes', nodes=nodes(1:nnodes) )

! object for plotting the nodeset

  call add_to_mesh ( mesh, object='coordinates', &
    coor=mesh%coor(nodes(1:nnodes),:) )

  nbx = nint(sqrt(real(nx)))
  nby = nint(sqrt(real(ny)))

  call add_to_mesh ( mesh, blocks=[nbx,nby] )

  call fill_mesh_parts ( mesh )

  lmesh => mesh  ! supply mesh to mapcoor for mapping coordinates

  call plot_points_curves ( plot_options, mesh, 'curves_xf.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh_xf.fig' )
  plot_options%objectpointcolor=4
  plot_options%objectpointsize=0.4
  call plot_objects ( plot_options, mesh, 'mesh_xf.fig', append=.true., &
    object1=1 )
  plot_options%objectpointcolor=5
  call plot_objects ( plot_options, mesh, 'mesh_xf.fig', append=.true.,  &
    object1=2 )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=3, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )
  call define_essential ( mesh, input_probdef, nodeset1=1 )

! define constraints on the object

  call define_constraint ( mesh, input_probdef, object=1, physq=1, &
    discretization='weak', elementdof=[2,0,2] )

! define constraints on curve

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, nglobalc=1 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, curve2=5, discretization='collocation', exclude=3 )

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

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, oldvectors=oldvectors, &
    coefficients=coefficients, buildvector=.false. )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=elementc, addmat=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=3, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficients  )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call toc ( 'build' )

  solver_options%real_storage=4._dp

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

  call toc ( 'solve' )

! post-processing

  oldvectors%s(1)%p => sol

! sampling for gnuplot

! write profile to a file for plotting with gnuplot

  call printtofile ( mesh, problem, 'vprofile_xf.out', curve=2, sysvector=sol )

  open( unit=15, file='sample_p_xf.out' )

! sample pressure on cylinder

  call fill_sample ( mesh, problem, sample_p, ndegfd=1, object=1, &
    elemsub=stokes_sample_pressure, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_p%nnodes
    write( unit=15, fmt=* ) sample_p%coor(i,:), sample_p%u(i,:)
  end do

  close( unit=15 )

  open( unit=15, file='sample_v_xf.out', recl=300 )

! sample velocity on cylinder

  call fill_sample ( mesh, problem, sample_v, ndegfd=2, object=1, &
    elemsub=stokes_sample_velocity, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_v%nnodes
    write( unit=15, fmt=* ) sample_v%coor(i,:), sample_v%u(i,:)
  end do

  close( unit=15 )

  open( unit=15, file='sample_vort_xf.out', recl=300 )

! sample vorticity on cylinder

  coefficients%i(13)=5

  call fill_sample ( mesh, problem, sample_vort, ndegfd=1, object=1, &
    elemsub=stokes_sample_deriv, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_vort%nnodes
    write( unit=15, fmt=* ) sample_vort%coor(i,:), sample_vort%u(i,:)
  end do

  close( unit=15 )

! post-processing: create mesh_plot for plotting

  call mesh_convert ( mesh, mesh_plot, userelmesh=userelmesh, warn=.false. )

  call fill_mesh_parts ( mesh_plot )

  call plot_mesh ( plot_options, mesh_plot, 'mesh_xf_plot.fig' )

  call toc ( 'plot_mesh' )

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

! create an object for interpolation on the original mesh

  warn_add_to_mesh_after_meshgen_parts = .false. ! .false. to suppress warning
  call add_to_mesh ( mesh, object='coordinates', coor=mesh_plot%coor )
  plotobj = mesh%nobjects
  call fill_mesh_parts_objects ( mesh, object1=plotobj )

  call toc ( 'fill_mesh_parts_objects' )

! sample pressure

  call fill_sample ( mesh, problem, sample_p, ndegfd=1, object=plotobj, &
    elemsub=stokes_sample_pressure, coefficients=coefficients, &
    oldvectors=oldvectors )

  pressure%u = sample_p%u(:,1)

! sample velocity

  call fill_sample ( mesh, problem, sample_v, ndegfd=2, object=plotobj, &
    elemsub=stokes_sample_velocity, coefficients=coefficients, &
    oldvectors=oldvectors )

  velocity%u = reshape ( transpose(sample_v%u), [2*mesh_plot%nnodes] )

! sample vorticity

  coefficients%i(13)=5

  call fill_sample ( mesh, problem, sample_vort, ndegfd=1, object=plotobj, &
    elemsub=stokes_sample_deriv, coefficients=coefficients, &
    oldvectors=oldvectors )

  vorticity%u = sample_vort%u(:,1)

  call toc ( 'fill_sample' )

  plot_options%plotboundary2=.false.

  call plot_color_contour ( plot_options, mesh_plot, problem_plot, &
    'pressure_xf.fig', vector=pressure )
  call plot_boundary ( plot_options, mesh, 'pressure_xf.fig', append=.true. )

  call plot_color_contour ( plot_options, mesh_plot, problem_plot, &
    'vorticity_xf.fig', vector=vorticity )
  call plot_boundary ( plot_options, mesh, 'vorticity_xf.fig', append=.true. )

  call plot_color_fill ( plot_options, mesh_plot, problem_plot, &
    'vorticity_fill_xf.fig', vector=vorticity )
  call plot_boundary ( plot_options, mesh, 'vorticity_fill_xf.fig', &
    append=.true. )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=pressure, &
    dataname='pressure', filename='cylinder_xf.vtk' )

  call write_vector_vtk ( mesh_plot, problem_plot, &
    filename='cylinder_xf.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true. )

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=vorticity, &
    filename='cylinder_xf.vtk', &
    dataname='vorticity', append=.true. )

  !print *, maxval(vorticity%u), minval(vorticity%u)
  !print *, minval(sqrt(sum(sample_vort%coor**2,dim=2)))

! delete all data including all allocated memory

  call delete ( problem, problem_plot )
  call delete ( input_probdef, input_probdef_plot )
  call delete ( mesh, mesh_plot, mesh_particle )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( sample_p, sample_vort, sample_v )

  deallocate ( nodes, d )

contains


! find internal nodes that are not connected to elements at the interface

  subroutine find_internal_zero_nodes

    integer :: i, elem

!   set internal nodes

    where ( d <= 0._dp )
      nodes = 1
    else where
      nodes = 0
    end where

    allocate ( nodes1(mesh%nnodes) )

    nodes1 = nodes

!   remove nodes connected to elements at the interface

    do elem = 1, mesh%grpnumel(1)
       if ( .not. all ( nodes(mesh%topology(1)%a(:,elem)) == &
                        nodes(mesh%topology(1)%a(1,elem)) ) ) then
         nodes1(mesh%topology(1)%a(:,elem)) = 0
       end if
    end do

!   count and collect the nodes

    nnodes = 0
    do i = 1, mesh%nnodes
      if ( nodes1(i) == 0 ) cycle
      nnodes = nnodes + 1
      nodes(nnodes) = i
    end do

    deallocate ( nodes1 )

  end subroutine find_internal_zero_nodes

end program xfem2
