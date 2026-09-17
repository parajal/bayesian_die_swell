! 2D Stokes problem on a square domain with a cylinder
! eXtended FEM approach with embedded Dirichlet condition on the cylinder.
! Quadtree+triangle subdivision for integration.
! Prebuild eltree for all elements crossing the interface.
! Optional non-zero imposed velocity of the cylinder.
! Similar to xfem5 but now with embedded Dirichlet condition

! functions module for xfem8


! Some common subroutines for the circular object in fic_dom1 and xfem1.

module subs8_m

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

end module subs8_m

module functions_xf8_m

  use tfem_elem_m
  use eltree_m

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

end module functions_xf8_m


! replacement module for the dummy one in the standard add-on

module stokes_usergauss_m

  use tfem_elem_m
  use functions_xf8_m
  use eltree_m

  implicit none

  save

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

! epsjac: subelements in smesh (submesh) having a jacobian smaller than
!         epsjac will be omitted.
  real(dp) :: epsjac = 0._dp

! some logicals
  logical :: printnumel = .false., printint = .false.

! array of pointers to an eltree (for userelmesh only)

  type(eltree_p), dimension(:,:), pointer :: ea

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

!   set eltree pointer (alias)

    eltree1 => oldvectors%ea(1)%p(elgrp,elem)%p

!   Set ninti

    nod = mesh%topology(elgrp)%a(:,elem)

    if ( associated(eltree1) ) then

!     element crosses the interface

      if ( printnumel ) then
        print *, 'number of subelements outside = ', &
          number_of_subelements ( eltree1, lsign=[1] )
      end if

!     composite integration scheme for levelset >= 0

      intopt%epsjac = epsjac

      ninti = number_of_integration_points ( eltree1, lsign=[1],&
        ninti=ninti_e, nintis=ninti_sub, integration_options=intopt )

      if ( printint ) then
        print *, 'number of integration points = ', ninti
      end if

    else if ( all ( d(nod) <= 0._dp ) ) then

!     fully in

      ninti = 0

    else if ( all ( d(nod) >= 0._dp ) ) then

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

!   set eltree pointer

    eltree1 => oldvectors%ea(1)%p(elgrp,elem)%p

    nod = mesh%topology(elgrp)%a(:,elem)

    if ( associated(eltree1) ) then

!     element crosses the interface

!     composite integration scheme for levelset >= 0

      call integration_points ( eltree1, xig, wg, lsign=[1], ninti=ninti_e, &
        xe=xe, we=we, nintis=ninti_sub, xs=xsub, ws=wsub, &
        integration_options=intopt )

    else if ( all ( d(nod) <= 0._dp ) ) then

!     fully in: generate fully zero elemmat and elemvec

      xig = 0; wg = 0

    else if ( all ( d(nod) >= 0._dp ) ) then

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
    type(eltree_t), pointer :: eltree


!   set eltree pointer for convenience

    eltree => ea(elgrp,elem)%p

    nod = mesh%topology(elgrp)%a(:,elem)

    if ( associated(eltree) ) then

!     element crosses the interface

      lelgrp = elgrp
      lelem = elem

!     add only elements in the outside region

      call eltree_to_mesh ( eltree, elmesh, lsign=[1], mapcoor=mapcoor )

      indicator = 2 ! replace element with submesh

    else if ( all ( d(nod) <= 0._dp ) ) then

!     fully in: remove element

      indicator = 1

    else if ( all ( d(nod) >= 0._dp ) ) then

!     fully out: leave element as it is.

      indicator = 0

    end if

  end subroutine userelmesh

end module stokes_usergauss_m


! the actual program

program xfem8

  use tfem_m
  use stokes_elements_m
  use hsl_ma41_m
  use subs8_m
  use functions_xf8_m
  use stokes_usergauss_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    nx=20,              & ! number of elements in x
    ny=20,              & ! number of elements in y
    numsplit = 3,       & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit
    numsplitmin = 1,    & ! relative size of all subelements will be at least
                          ! as small as 1/2^numsplitmin
    gausse = 3,         & ! Gauss integration for the eltree subelements
    gausssub = 6,       & ! Gauss integration for the submesh subelements
    gauss_ie = 3,       & ! Gauss integration for the interface elements
    gauss = 3             ! 3x3 Gauss integration

  real(dp), parameter :: &
    lx = 8._dp,        & ! width of domain
    ly = 8._dp,        & ! height of domain
    radius = 1._dp,    & ! radius of the cylinder
    center(2) = [ 0._dp, 0._dp ], & ! center position of the cylinder
    eta = 1._dp,       & ! viscosity
    kappa = 1.0_dp,     & ! embedded Dirichlet viscosity parameter
    flowrate = 10._dp, & ! flowrate
    split_threshold = 1.e-2_dp, & ! levelset value for being "close" enough to the interface for tree splitting
    U(2) = [ 1._dp, 2._dp ], & ! imposed velocity of the cylinder if coefficients%i(42) = 0
    epsdef = 0.02_dp    ! deformation of the mesh

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t), target :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options

  integer :: nnodes, nbx, nby

  integer, dimension(:), allocatable :: nodes
  real(dp), dimension(:,:), allocatable :: coor

! array of pointers to an eltree

  type(eltree_p), target, dimension(:,:), allocatable :: eltree_array

! set some parameters in modules

! module functions1_m

  rpl = radius ! radius of the circular object
  cpl = center ! initial position of the center of the object

! module stokes_usergauss_m

  intrule_e = gausse  ! Gauss integration for the eltree subelements
  intrule_sub = gausssub  ! Gauss integration for the eltree subelements
  epsjac = 0._dp ! subelements in smesh (submesh) having a jacobian smaller
                 ! than epsjac will be omitted.

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gauss ]
  coefficients%i(12:) = 0
  coefficients%i(36) = -1 ! normal to fluid points inwards to the object
  coefficients%i(37) = 2 ! Gauss points defined by user subroutines
  coefficients%i(41) = gauss_ie ! Gauss points for interface elements
  coefficients%i(42) = -1  ! -1: impose zero velocity, 0: impose non-zero velocity U

  coefficients%r(1) = eta
  coefficients%r(2:) = 0
  coefficients%r(6) = flowrate
  coefficients%r(11) = kappa
  coefficients%r(12:13) = U

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

! create object for particle boundary (only for plotting)

  nnodes = 8 * nx * radius / lx

  allocate ( coor(nnodes,2) )

  rp = radius ! radius of the circular object
  xpc = center ! initial position of the center of the object

  call objectscoor ( 1, coor )

  call add_to_mesh ( mesh, object='coordinates', coor=coor )

  lmesh => mesh  ! supply mesh to mapcoor for mapping coordinates

! allocate array of pointers to an eltree (eltree for all elements, one group)

  allocate ( eltree_array(1,mesh%grpnumel(1)) )

! nodeset for Dirichlet on internal nodes

  allocate ( nodes(mesh%nnodes), d(mesh%nnodes) )

  d = levelset ( mesh%coor ) ! set levelset for all nodes

  call create_eltree_in_elements ( numsplit )  ! eltree for elements crossing the interface

  call create_elementset_from_eltree ! elementset for elements crossing the interface

  call find_internal_zero_nodes

  call add_to_mesh ( mesh, nodeset='nodes', nodes=nodes(1:nnodes) )

! object for plotting the nodeset

  call add_to_mesh ( mesh, object='coordinates', &
    coor=mesh%coor(nodes(1:nnodes),:) )

  nbx = nint(sqrt(real(nx)))
  nby = nint(sqrt(real(ny)))

  call add_to_mesh ( mesh, blocks=[nbx,nby] )

  call fill_mesh_parts ( mesh )

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

  call create ( oldvectors, nsysvec=1, nelta=1 )

  oldvectors%ea(1)%p => eltree_array
  ea => eltree_array ! userelmesh does not have oldvectors as an argument

! create system matrix (non-symmetric due to embedded Dirichlet)

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, oldvectors=oldvectors, &
    coefficients=coefficients, buildvector=.false. )

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_open_boundary_eltree, oldvectors=oldvectors, &
    coefficients=coefficients, addmat=.true., elementset=1 )

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    physqrow=[1], physqcol=[1], &
    elemsub=stokes_embedded_dirichlet_eltree, oldvectors=oldvectors, &
    coefficients=coefficients, addmatvec=.true., elementset=1 )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficients  )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


  call solve_system_ma41 ( sysmatrix, rhsd, sol, solver_options=solver_options )

  print *, sum(abs(sol%u))/sol%n

  call delete_eltree_in_elements ! remove eltree for all elements

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

  deallocate ( coor, nodes, d, eltree_array )

contains


! create eltree in the elements which are crossed by the interface

  subroutine create_eltree_in_elements ( numsplit )

    integer, intent(in) :: numsplit

    integer :: nod(mesh%element(1)%numnod), elem
    real(dp) :: coor(2,2)


!   loop all elements

    do elem = 1, mesh%grpnumel(1)

!     nodal points

      nod = mesh%topology(1)%a(:,elem)

      if ( all ( d(nod) > split_threshold ) .or. all ( d(nod) < -split_threshold ) ) then

!       all nodes far from the interface

        cycle

      else

!       element possibly contains an interface, start subdivide

        allocate ( eltree_array(1,elem)%p )

!       fill root node of the eltree

        coor(1,:) = -1._dp  ! lower corner of the reference region (-1,-1,-1)
        coor(2,:) =  1._dp  ! upper corner of the reference region (1,1,1)

        call fill_node_eltree ( eltree_array(1,elem)%p, coor )

!       divide root reference domain into subdomains

        lelgrp = 1
        lelem = elem

        call subdivide ( eltree_array(1,elem)%p, levelset=levelset, numsplit=numsplit, &
          split_threshold=split_threshold, numsplitmin=numsplitmin, &
          mapcoor=mapcoor, submesh=.true., intmesh=.true. )

!       check whether interface is in element

        if ( any( number_of_subelements_vector(eltree_array(1,elem)%p,lsign=[-1,1]) == 0 ) ) then

!         no interface, remove eltree

          call delete(eltree_array(1,elem)%p)  ! delete actual eltree
          deallocate( eltree_array(1,elem)%p ) ! delete pointer target (pointer becomes disassociated)

        end if

      end if

    end do

  end subroutine create_eltree_in_elements


! create elementset of elements which are crossed by the interface

  subroutine create_elementset_from_eltree ( replace )

    integer, intent(in), optional :: replace

    integer :: elements(mesh%grpnumel(1)), numelem, elem


!   loop all elements

    numelem = 0

    do elem = 1, mesh%grpnumel(1)

!     set eltree pointer for convenience

      if ( associated(eltree_array(1,elem)%p) ) then

!       element crosses the interface

        numelem = numelem + 1

        elements(numelem) = elem

      end if

    end do

!   add elementset

    call add_to_mesh ( mesh, elementset='elements', elements=elements(1:numelem), replace=replace )

  end subroutine create_elementset_from_eltree


! delete eltree in the elements

  subroutine delete_eltree_in_elements

    integer :: elem

!   loop all elements

    do elem = 1, size(eltree_array,2)

      if ( associated(eltree_array(1,elem)%p) ) then

        call delete(eltree_array(1,elem)%p) ! delete actual eltree
        deallocate(eltree_array(1,elem)%p) ! delete pointer target (pointer becomes disassociated)

      end if

    end do

  end subroutine delete_eltree_in_elements


! find internal nodes that are not connected to elements at the interface

  subroutine find_internal_zero_nodes

    integer :: elem, i

!   set internal nodes

    where ( d <= 0._dp )
      nodes = 1
    else where
      nodes = 0
    end where

!   remove nodes connected to elements at the interface

    do elem = 1, mesh%grpnumel(1)
      if ( associated(eltree_array(1,elem)%p) ) then
         nodes(mesh%topology(1)%a(:,elem)) = 0
       end if
    end do

!   count and collect the nodes (overwrite array nodes along the way)

    nnodes = 0
    do i = 1, mesh%nnodes
      if ( nodes(i) == 0 ) cycle
      nnodes = nnodes + 1
      nodes(nnodes) = i
    end do

  end subroutine find_internal_zero_nodes

end program xfem8
