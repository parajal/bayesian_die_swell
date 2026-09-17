! 2D Stokes problem on a square domain with a cylinder
! eXtended FEM approach with weak constraint on the cylinder.
! Quadtree subdivision for integration.

! functions module for xfem1


! Some common subroutines for the circular object in fic_dom1 and xfem1.

module subs1_m

  use tfem_elem_m
  use math_defs_m

  implicit none

  integer :: ndf = 9
  integer :: ndfb = 3, nodalpb = 3
  integer :: ndflb = 2  ! number of degrees of freedom of the Lag. Mult.
  real(dp) :: rp = 1, xpc(2) = 0

contains

! objectscoor defines the coordinates of the objects

  subroutine objectscoor ( objectnr, coor )
    integer, intent(in) :: objectnr
    real(dp), dimension(:,:), intent(inout) :: coor

    integer :: np, i
    real(dp) :: p(size(coor,1))

    np = size(coor,1)
    p = [(2*pi/np*(i-1),i=1,np)]

    coor(:,1) = rp * sin(p) + xpc(1)
    coor(:,2) = rp * cos(p) + xpc(2)

  end subroutine objectscoor

! elementc is the element subroutine for the constraints on the objects

  subroutine elementc ( mesh, problem, constr, elem, node, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemmat2, elemmatadd, &
    elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: object, i, j
    !real(dp) :: phi(1,ndf), xr(1,2), r(2), psi(1,ndflb), x(nodalpb,2)
    real(dp) :: phi(1,ndf), xr(1,2), psi(1,ndflb), x(nodalpb,2)
    real(dp) :: theta(1,ndfb), dtheta(1,ndfb,1), dxdxi(1,2), curvel(1)


    object = problem%constraints(constr)%object

!   reference coordinates of the integration point (node) of the element (elem)

    xr(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)
    !r = mesh%objects(object)%coor_int(node,:,elem) - xpc

!   shape function of the velocity in the integration point

    call shape_quad_Q2 ( xr, phi )

!   shape function of the Langrangian multiplier in the integration point

    call shape_line_P1 ( mesh%objects(object)%xig(node:node,1), psi )

!   shape function of the quadratic curve in the integration point

    call shape_line_P2 ( mesh%objects(object)%xig(node:node,1), theta, &
      dtheta(:,:,1) )

!   compute deformed element

    call get_coordinates_object ( mesh, elem, x, object )

    call isoparametric_deformation_curve ( x, dtheta(:,:,1), dxdxi, curvel )

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

!   two constraints (vectorial)
!
!      u  = 0
!      -    -
!
!   or in components
!
!      u = 0
!      v = 0
!
!   the matrix A therefore becomes
!
!     A =  [ phi    0 ]
!          [   0  phi ]
!
!   In the weak version the equations are multiplied by the test function psi
!   and integrated over the element.
!   NOTE: this is the contribution of a single integration point (node) in
!   a single element (elem).

    if ( matrix ) then

      elemmat = 0

      do i = 1, ndflb
        do j = 1, ndf
          elemmat(i,j) = psi(1,i) * phi(1,j) * &
                                   curvel(1) * mesh%objects(object)%wg(node)
        end do
      end do
      elemmat(ndflb+1:2*ndflb,ndf+1:2*ndf) = elemmat(1:ndflb,1:ndf)

      elemmatadd = 0

    end if

  end subroutine elementc

end module subs1_m

module functions_xf1_m

  use tfem_elem_m

  implicit none

  save

  type(mesh_t), pointer :: lmesh ! mesh used for mapping of reference element
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

end module functions_xf1_m


! replacement module for the dummy one in the standard add-on

module stokes_usergauss_m

  use tfem_elem_m
  use functions_xf1_m
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
  integer :: ninti_s, intrule_e, ninti_e
  real(dp), allocatable :: ws(:), xs(:,:), we(:), xe(:,:)

! integration options for setting the composite Gauss integration on the eltree
  type(integration_options_t) :: intopt

! numsplit parameters
!  ns : numsplit parameter in subdivide for the integration scheme
!  nsmin : numsplitmin parameter in subdivide
!  ns_plot : numsplit parameter in subdivide for the plot mesh
  integer :: ns = 2, nsmin = 1, ns_plot = 2

! some logicals
  logical :: warn = .false., refine = .false., printnumel = .false., &
    printint = .false.

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
    type(gauss_t) :: gausse

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
        numsplitmin=nsmin, mapcoor=mapcoor )

      if ( refine .and. number_of_subelements(eltree,lsign=[1]) == 0 ) then

!       further refine
        call delete (eltree)
        call fill_node_eltree ( eltree, coor )
        call subdivide ( eltree, levelset=levelset, numsplit=ns+3, &
          numsplitmin=nsmin, mapcoor=mapcoor )

      end if

      if ( warn .and. number_of_subelements(eltree,lsign=[1]) == 0 ) then
        write(*,'(a)') 'Warning no integration points in outer region'
      end if

      if ( printnumel ) then
        print *, 'number of subelements outside = ', &
          number_of_subelements(eltree,lsign=[1])
        print *, 'number of subelements crossing interface = ', &
          number_of_subelements(eltree,lsign=[0])
      end if

!     composite integration scheme for levelset >= 0

      intopt%midp = .true.
      intopt%half = .true.

      ninti = number_of_integration_points ( eltree, lsign=[0,1],&
        ninti=ninti_e, integration_options=intopt )

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

      call integration_points ( eltree, xig, wg, lsign=[0,1], ninti=ninti_e, &
        xe=xe, we=we, integration_options=intopt )

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

      deallocate ( we, xe, ws, xs )

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
        numsplitmin=nsmin, mapcoor=mapcoor )

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

program xfem1

  use tfem_m
  use stokes_elements_m
  use hsl_ma57_m
  use io_utils_m
  use subs1_m
  use functions_xf1_m
  use stokes_usergauss_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    nx=10,              & ! number of elements in x
    ny=10,              & ! number of elements in y
    numsplit = 7,       & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit
    numsplitmin = 1,    & ! relative size of all subelements will be at least
                          ! as small as 1/2^numsplitmin
    numsplit_plot = 2,  & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit for plotting
    nsubinto = 10,      & ! number of subintegration points on object
    !nsubint = 1,        & ! number of subintegration points
    gausse = 2,         & ! Gauss integration for the eltree subelements
    gauss = 3             ! 3x3 Gauss integration

  real(dp), parameter :: &
    lx = 8._dp,        & ! width of domain
    ly = 8._dp,        & ! height of domain
    radius = 1._dp,    & ! radius of the cylinder
    center(2) = [ 0._dp, 0._dp ], & ! center position of the cylinder
    eta = 1._dp,       & ! viscosity
    flowrate = 10._dp, & ! flowrate
    epsdef = 0.020_dp    ! deformation of the mesh

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
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(sample_t) :: sample_p, sample_v, sample_vort
  type(solver_options_ma57_t) :: solver_options

  integer :: i, elem, nelem, nnodes, plotobj, nbx, nby

  integer, dimension(:), allocatable :: nodes, nodes1

! set some parameters in modules

! module functions1_m

  rpl = radius ! radius of the circular object
  cpl = center ! initial position of the center of the object

! module stokes_usergauss_m

  intrule_e = gausse  ! Gauss integration for the eltree subelements
  ns = numsplit ! parameter in subdivide for the integration scheme
  nsmin = numsplitmin ! parameter in subdivide
  ns_plot = numsplit_plot ! parameter in subdivide for the plot mesh

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
  mesh%coor(:,2) = mesh%coor(:,2) * ( 1 - epsdef*(1-2*abs(mesh%coor(:,2))/lx) )

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

  solver_options%real_storage=4._dp

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

! post-processing

  oldvectors%s(1)%p => sol

! sampling for gnuplot

! post-processing: create mesh_plot for plotting

  call mesh_convert ( mesh, mesh_plot, userelmesh=userelmesh, warn=.false. )

  call fill_mesh_parts ( mesh_plot )

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

  print *, maxval(vorticity%u), minval(vorticity%u)
  print *, minval(sqrt(sum(sample_vort%coor**2,dim=2)))

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

    integer :: elem, i

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

end program xfem1
