! 3D Stokes problem on a cubic domain with a freely floating sphere (F=0,M=0).
! eXtended FEM approach with embedded particle condition on the sphere.
! Quadtree+triangle subdivision for integration.
! Prebuild eltree for all elements crossing the interface.
! Removal of elements from eltree if integration volume is too small.
! High-order exact integration for small integration volumes.
! Optional symmetric matrix.
! Optional Nitsche's method.
! Similar to xfem13 but now 3D.

module subs14_m

  use math_defs_m

  implicit none

  real(dp) :: rp = 1, xpc(3) = 0

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
    coor(:,3) = xpc(3)

  end subroutine objectscoor

end module subs14_m


! functions module for xfem14

module functions_xf14_m

  use tfem_elem_m
  use eltree_m

  implicit none

  save

! mesh used for mapping of reference element
  type(mesh_t), pointer :: lmesh => null()

  integer :: lelem, lelgrp
! radius and center of the sphere
  real(dp) :: rpl = 1._dp, cpl(3) = [ 0._dp, 0._dp, 0._dp ]

contains


! levelset for the sphere

  function levelset ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1)) :: levelset

    real(dp) :: c(size(x,1),size(x,2))

    c(:,1) = cpl(1)
    c(:,2) = cpl(2)
    c(:,3) = cpl(3)

!   sphere with center at cpl and radius rpl
    levelset = sqrt(sum((x-c)**2,dim=2)) - rpl

  end function levelset


! function for mapping reference coordinates to the real coordinates
! in building the eltree within an element.

  function mapcoor ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1),size(x,2)) :: mapcoor

    real(dp) :: phi(size(x,1),27), xnod(27,3)

    call shape_hexa_Q2 ( x, phi )

    call get_coordinates ( lmesh, lelgrp, lelem, xnod )

    mapcoor = matmul ( phi, xnod ) ! isoparametric

  end function mapcoor

end module functions_xf14_m


! module for the user gauss routines

module stokes_usergauss_m

  use tfem_elem_m
  use functions_xf14_m
  use eltree_m

  implicit none

  save

! the levelset function in all the nodes (it is a signed distance function
! in this case of a simple sphere).
  real(dp), dimension(:), allocatable :: d

! integration volume near interface elements
  real(dp), dimension(:), allocatable :: vol

! the standard integration scheme:
!   ninti_s: number of integration points
!   xs, ws: points and weights
! the integration scheme on subelements in the eltree:
!   intrule_e: the integration rule
!   ninti_e: number of integration points
!   xe, we: points and weights
! the integration scheme on subelements in the submesh (tetrahedrons):
!   intrule_sub: the integration rule
!   intrule_sub_small: the integration rule for small integration volumes
!   ninti_sub: number of integration points
!   ninti_sub_small: number of integration points for small integration volumes
!   xsub, wsub: points and weights
!   xsub_small, wsub_small: points and weights for small integration volumes
  integer :: ninti_s, intrule_e, ninti_e, intrule_sub, intrule_sub_small, &
             ninti_sub, ninti_sub_small
  real(dp), allocatable :: ws(:), xs(:,:), we(:), xe(:,:), wsub(:), xsub(:,:), &
                           wsub_small(:), xsub_small(:,:)

! integration options for setting the composite Gauss integration on the eltree
  type(integration_options_t) :: intopt

! some logicals
  logical :: printnumel = .false., printint = .false.

! array of pointers to an eltree (for userelmesh only)

  type(eltree_p), dimension(:,:), pointer :: ea

! treshold for small volumes

  real(dp) :: vol_small = 1e-2_dp

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
    type(gauss_t) :: gausse, gausssub, gausssub_small

    if ( first ) then

!     first element in this group

!     Gauss rule for the standard case

      ninti_s = ninti

      allocate ( ws(ninti_s), xs(ninti_s,3) )

      call set_Gauss_integration ( gauss, xs, ws )

!     Gauss rule for the subdivided elements in the eltree

      gausse = gauss
      gausse%intrule = intrule_e

      call set_ninti ( gausse, ninti_e )

      allocate ( we(ninti_e), xe(ninti_e,3) )

      call set_Gauss_integration ( gausse, xe, we )

!     Gauss rule for the subdivided elements in the submesh

      gausssub%globalshape = 'tetrahedron'
      gausssub%intrule = intrule_sub
      gausssub%inttype = 3  ! use numerical Gauss values

      call set_ninti ( gausssub, ninti_sub )

      allocate ( wsub(ninti_sub), xsub(ninti_sub,3) )

      call set_Gauss_integration ( gausssub, xsub, wsub )

!     Gauss rule for the subdivided elements in the submesh (small volumes)

      gausssub_small%globalshape = 'tetrahedron'
      gausssub_small%intrule = intrule_sub_small
      gausssub_small%inttype = 3  ! use numerical Gauss values

      call set_ninti ( gausssub_small, ninti_sub_small )

      allocate ( wsub_small(ninti_sub_small), xsub_small(ninti_sub_small,3) )

      call set_Gauss_integration ( gausssub_small, xsub_small, wsub_small )

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

      if ( vol(elem) > vol_small ) then
        ninti = number_of_integration_points ( eltree1, lsign=[1],&
          ninti=ninti_e, nintis=ninti_sub, integration_options=intopt )
      else
        ninti = number_of_integration_points ( eltree1, lsign=[1],&
          ninti=ninti_e, nintis=ninti_sub_small, integration_options=intopt )
      end if

      if ( printint ) then
        print *, 'number of integration points = ', ninti
      end if

    else if ( all ( d(nod) >= 0._dp ) ) then

!     fully out

      ninti = ninti_s

    else

!     fully in or removed element due to small integration volume

      ninti = 0

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

      if ( vol(elem) > vol_small ) then
        call integration_points ( eltree1, xig, wg, lsign=[1], &
          ninti=ninti_e, xe=xe, we=we, nintis=ninti_sub, xs=xsub, ws=wsub, &
          integration_options=intopt )
      else
        call integration_points ( eltree1, xig, wg, lsign=[1], &
          ninti=ninti_e, xe=xe, we=we, nintis=ninti_sub_small, &
          xs=xsub_small, ws=wsub_small, integration_options=intopt )
      end if

    else if ( all ( d(nod) >= 0._dp ) ) then

!     fully out: standard integration

      call set_Gauss_integration ( gauss, xig, wg )

    else

!     fully in or removed element due to small integration volume
!     generate fully zero elemmat and elemvec

      xig = 0; wg = 0

    end if

!   delete some allocated memory

    if ( last ) then

      deallocate ( we, xe, ws, xs, wsub, xsub, wsub_small, xsub_small )

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

    else if ( all ( d(nod) >= 0._dp ) ) then

!     fully out: leave element as it is.

      indicator = 0

    else

!     fully in or removed element due to small integration volume

      indicator = 1

    end if

  end subroutine userelmesh

end module stokes_usergauss_m


! the actual program

program xfem14

  use tfem_m
  use stokes_elements_m
  use hsl_ma57_m
  use hsl_ma41_m
  use io_utils_m
  use figplot_m
  use subs14_m
  use functions_xf14_m
  use stokes_usergauss_m
  use timer_m
  use io_utils_m
  use limits_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    nx=10,              & ! number of elements in x
    ny=10,              & ! number of elements in y
    nz=10,              & ! number of elements in z
    numsplit = 2,       & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit
    numsplitmin = 1,    & ! relative size of all subelements will be at least
                          ! as small as 1/2^numsplitmin
    numsplit_plot = 2,  & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit for plotting
    gausse = 3,         & ! Gauss integration for the eltree subelements
    gausssub_small = 12,& ! Gauss integration for the submesh subelements for
                          ! small integration volumes
    gausssub = 4,       & ! Gauss integration for the submesh subelements
    gauss_ie = 6,       & ! Gauss integration for the interface elements
    gauss = 3             ! 3x3 Gauss integration

  real(dp), parameter :: &
    lx = 8._dp,        & ! width of domain
    ly = 8._dp,        & ! depth of domain
    lz = 8._dp,        & ! height of domain
    radius = lx/4,     & ! radius of the sphere
    center(3) = [ 0._dp, 0._dp, 0.75_dp * radius ], & ! center of the sphere
    eta = 1._dp,       & ! viscosity
    kappa = 1._dp,     & ! embedded Dirichlet viscosity parameter
    KN = 25._dp,       & ! Nitsche's factor
    flowrate = 10._dp, & ! flowrate
    split_threshold = 1.e-2_dp, & ! levelset value for being "close" enough
                                  ! to the interface for tree splitting
    epsvol = 1e-6_dp, & ! elements with integration volume smaller are removed
                        ! from the eltree_array. This means that:
                        ! 1) these elements are "fully inside" and no
                        !    "virtual/extended degrees" are generated.
                        ! 2) the interface integration is ignored and therefore
                        !    the interface has a small "hole".
                        ! Note, that epsvol is relative to the volume of an
                        ! element, so epsvol=1 is a full element.
    epsvol_small = 1e-2_dp, & ! elements with integration volume smaller
                              ! are integrated with integration rule
                              ! gausssub_small, otherwise with gaussub
    epsdef = 0.0_dp    ! deformation of the mesh

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t), target :: mesh
  type(mesh_t) :: mesh_plot
  type(input_probdef_t) :: input_probdef, input_probdef_plot
  type(problem_t) :: problem, problem_plot
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, dudy
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(sample_t) :: sample_p, sample_v, sample_dudy
  type(solver_options_ma57_t) :: solver_options_ma57
  type(solver_options_ma41_t) :: solver_options_ma41

  integer :: i, nnodes, plotobj, nbx, nby, nbz

  integer, dimension(:), allocatable :: nodes
  real(dp), dimension(:,:), allocatable :: coor

  integer :: presnod(8) = [1,3,9,7,19,21,27,25]

! array of pointers to an eltree

  type(eltree_p), target, dimension(:,:), allocatable :: eltree_array

  logical, parameter :: Nitsche = .false., symmetric = .true.

  real(dp) :: up(6)

  logical, parameter :: full_figplot_mesh = .true.

! set some parameters in modules

! module limits_m

  !SUBDIVIDE_HEX6 = .true.  ! if SUBDIVIDE_HEX6=.true., the interface is
                           ! continuous, but it is more expensive.



! module timer_m

  timer = .true. ! set to .true. to show cpu time output

! module functions_xf14_m

  rpl = radius ! radius of the circular object
  cpl = center ! initial position of the center of the object

! module stokes_usergauss_m

  intrule_e = gausse  ! Gauss integration for the eltree subelements
  intrule_sub = gausssub  ! Gauss integration for the eltree subelements
  intrule_sub_small = gausssub_small  ! Gauss integration for the eltree
                                      ! subelements (small integration volumes)
  vol_small = epsvol_small  ! elements with integration volume smaller
                            ! are integrated with integration rule
                            ! gausssub_small, otherwise with gaussub

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
  coefficients%i(46) = 1  ! transposed open boundary term:
                          ! 0: absent, -1: Baumann-Oden, 1: Nitsche (symmetric)

  coefficients%r(1) = eta
  coefficients%r(2:) = 0
  coefficients%r(6) = flowrate
  coefficients%r(11) = kappa
  coefficients%r(15) = eta*KN*nx/lx ! Nitsche's factor

  coefficients%r(16:18) = center


  coefficients%set_ninti_user => set_ninti_user
  coefficients%set_Gauss_integration_user => set_Gauss_integration_user

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  call set_mesh_options ( meshgen_options, elshape=14, nx=nx, ny=ny, nz=nz, &
    lx=lx, ly=ly, lz=lz, ox=-lx/2, oy=-ly/2, oz=-lz/2, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

! deform the mesh a little bit to avoid zero or very small integration volumes

  mesh%coor(:,1) = mesh%coor(:,1) * ( 1 - epsdef*(1-2*abs(mesh%coor(:,1))/lx) )
  mesh%coor(:,2) = mesh%coor(:,2) * ( 1 - epsdef*(1-2*abs(mesh%coor(:,2))/lx) )
  mesh%coor(:,3) = mesh%coor(:,3) * ( 1 - epsdef*(1-2*abs(mesh%coor(:,3))/lx) )

! create object for particle boundary (only for sampling)

  nnodes = 8 * nx * radius / lx

  allocate ( coor(nnodes,3) )

  rp = radius ! radius of the circular object
  xpc = center ! initial position of the center of the object

  call objectscoor ( 1, coor )

  call add_to_mesh ( mesh, object='coordinates', coor=coor )

  lmesh => mesh  ! supply mesh to mapcoor for mapping coordinates

! allocate array of pointers to an eltree (eltree for all elements, one group)

  allocate ( eltree_array(1,mesh%grpnumel(1)), vol(mesh%grpnumel(1)) )

! nodeset for Dirichlet on internal nodes

  allocate ( nodes(mesh%nnodes), d(mesh%nnodes) )

  d = levelset ( mesh%coor ) ! set levelset for all nodes

  call create_eltree_in_elements ( numsplit )  ! eltree for elements
                                               ! crossing the interface

  call create_elementset_from_eltree ! elementset for elements
                                     ! crossing the interface

  call add_to_mesh ( mesh, elementset='nodes', elementsetnr=1 )

  call find_internal_zero_nodes

  call add_to_mesh ( mesh, nodeset='nodes', nodes=nodes(1:nnodes) )

! object for plotting the nodeset

  call add_to_mesh ( mesh, object='coordinates', &
    coor=mesh%coor(nodes(1:nnodes),:) )

  nbx = nint(sqrt(real(nx)))
  nby = nint(sqrt(real(ny)))
  nbz = nint(sqrt(real(nz)))

  call add_to_mesh ( mesh, blocks=[nbx,nby,nbz] )

  call fill_mesh_parts ( mesh )

  plot_options%viewpoint = [ 1._dp, 1.7_dp, 0.5_dp ]

  call plot_points_curves ( plot_options, mesh, 'curves_xf14.fig' )

  if ( full_figplot_mesh ) then

    call plot_mesh ( plot_options, mesh, 'mesh_xf14.fig' )
    plot_options%objectpointcolor=4
    plot_options%objectpointsize=0.4
    call plot_objects ( plot_options, mesh, 'mesh_xf14.fig', append=.true., &
      object1=1 )
    plot_options%objectpointcolor=5
    call plot_objects ( plot_options, mesh, 'mesh_xf14.fig', append=.true.,  &
      object1=2 )

  else

    call plot_mesh ( plot_options, mesh, 'mesh_xf14.fig', surfaces=[3,4,6] )

  end if

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(presnod,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, surface1=1, surface2=2, physq=1 )
  call define_essential ( mesh, input_probdef, surface1=4, physq=1 )
  call define_essential ( mesh, input_probdef, surface1=6, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )
  call define_essential ( mesh, input_probdef, nodeset1=1 )

! define constraints on surfaces

  call define_constraint ( mesh, input_probdef, &
    physq=1, surface1=3, nglobalc=1 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, surface1=3, surface2=5, discretization='collocation', &
    excludecurves=[2,6,7,10] )

! define particle using a global constraint on elementset

  call define_constraint ( mesh, input_probdef, elementset1=1, nglobalc=6, &
    diagonal_block=.true. )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

  sol%u = 0

! create the structure oldvectors

  call create ( oldvectors, nsysvec=1, nelta=1 )

  oldvectors%ea(1)%p => eltree_array
  ea => eltree_array ! userelmesh does not have oldvectors as an argument

! create system matrix

  if ( symmetric ) then
    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric=.true. )
  else
    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  end if

  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

  call tic

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, oldvectors=oldvectors, &
    coefficients=coefficients, buildvector=.false. )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=3, elemsub5=stokes_open_boundary_particle_eltree, &
    oldvectors=oldvectors, coefficients=coefficients, &
    addmat=.true. )

  if ( Nitsche ) then

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=3, elemsub4=stokes_Nitsche_particle_eltree, &
      oldvectors=oldvectors, coefficients=coefficients, &
      addmatvec=.true. )

  else

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=3, elemsub4=stokes_embedded_particle_eltree, &
      oldvectors=oldvectors, coefficients=coefficients, &
      addmatvec=.true. )

  end if

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficients  )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call toc ( 'build' )

  if ( symmetric ) then

    solver_options_ma57%real_storage = 1.1

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_ma57 )

  else

    !solver_options_ma41%real_storage = 2.0

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_ma41 )

  end if

  call toc ( 'solve' )

  call get_sysvector_constraint ( mesh, problem, sol, constraint=3, u=up )

  print *, 'up = ', up

  call delete_eltree_in_elements ! remove eltree for all elements

! post-processing

! create eltree for all elements crossing the interface

  call create_eltree_in_elements ( numsplit_plot )

  oldvectors%s(1)%p => sol

! sampling for gnuplot

  open( unit=15, file='sample_p_xf14.out', recl=300 )

! sample pressure on object

  call fill_sample ( mesh, problem, sample_p, ndegfd=1, object=1, &
    elemsub=stokes_sample_pressure, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_p%nnodes
    write( unit=15, fmt=* ) sample_p%coor(i,:), sample_p%u(i,:)
  end do

  close( unit=15 )

  open( unit=15, file='sample_v_xf14.out', recl=300 )

! sample velocity on sphere

  call fill_sample ( mesh, problem, sample_v, ndegfd=3, object=1, &
    elemsub=stokes_sample_velocity, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_v%nnodes
    write( unit=15, fmt=* ) sample_v%coor(i,:), sample_v%u(i,:)
  end do

  close( unit=15 )

  open( unit=15, file='sample_dudy_xf14.out', recl=300 )

! sample dudy on sphere

  coefficients%i(13)=2

  call fill_sample ( mesh, problem, sample_dudy, ndegfd=1, object=1, &
    elemsub=stokes_sample_deriv, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_dudy%nnodes
    write( unit=15, fmt=* ) sample_dudy%coor(i,:), sample_dudy%u(i,:)
  end do

  close( unit=15 )

! post-processing: create mesh_plot for plotting

  call mesh_convert ( mesh, mesh_plot, userelmesh=userelmesh, warn=.false. )

  call fill_mesh_parts ( mesh_plot )

  call plot_mesh ( plot_options, mesh_plot, 'mesh_xf14_plot.fig' )

  call write_mesh_vtk ( mesh_plot, filename='mesh_xf14_plot.vtk' )

  call toc ( 'plot_mesh' )

! problem definition for post_processing

  call create_input_probdef ( mesh_plot, input_probdef_plot, nvec=2 )

  do i = 1, mesh_plot%nelgrp
    input_probdef_plot%elementdof(i)%a(:) = 0
    input_probdef_plot%vec_elementdof(i)%a(:,1) = 3  ! velocity
    input_probdef_plot%vec_elementdof(i)%a(:,2) = 1  ! scalar
  end do

  call problem_definition ( input_probdef_plot, mesh_plot, problem_plot )

  call create_vector ( problem_plot, velocity, vec=1 )
  call create_vector ( problem_plot, pressure, vec=2 )
  call create_vector ( problem_plot, dudy, vec=2 )

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

  call fill_sample ( mesh, problem, sample_v, ndegfd=3, object=plotobj, &
    elemsub=stokes_sample_velocity, coefficients=coefficients, &
    oldvectors=oldvectors )

  velocity%u = reshape ( transpose(sample_v%u), [3*mesh_plot%nnodes] )

! sample dudy

  coefficients%i(13)=2

  call fill_sample ( mesh, problem, sample_dudy, ndegfd=1, object=plotobj, &
    elemsub=stokes_sample_deriv, coefficients=coefficients, &
    oldvectors=oldvectors )

  dudy%u = sample_dudy%u(:,1)

  call toc ( 'fill_sample' )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=pressure, &
    dataname='pressure', filename='sphere_xf14.vtk' )

  call write_vector_vtk ( mesh_plot, problem_plot, &
    filename='sphere_xf14.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true. )

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=dudy, &
    filename='sphere_xf14.vtk', &
    dataname='dudy', append=.true. )

! delete all data including all allocated memory

  call delete_eltree_in_elements

  call delete ( problem, problem_plot )
  call delete ( input_probdef, input_probdef_plot )
  call delete ( mesh, mesh_plot )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, dudy )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( sample_p, sample_dudy, sample_v )

  deallocate ( coor, nodes, d, eltree_array, vol )

contains


! create eltree in the elements which are crossed by the interface

  subroutine create_eltree_in_elements ( numsplit )

    integer, intent(in) :: numsplit

    integer :: nod(mesh%element(1)%numnod), elem
    real(dp) :: coor(2,3)
    real(dp) :: minvalvol


!   loop all elements

    minvalvol = 1

    do elem = 1, mesh%grpnumel(1)

!     nodal points

      nod = mesh%topology(1)%a(:,elem)

      if ( all ( d(nod) > split_threshold ) .or. &
                         all ( d(nod) < -split_threshold ) ) then

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

        call subdivide ( eltree_array(1,elem)%p, levelset=levelset, &
          numsplit=numsplit, split_threshold=split_threshold, &
          numsplitmin=numsplitmin, mapcoor=mapcoor, submesh=.true., &
          intmesh=.true. )

!       check whether interface is in element

        if ( any( number_of_subelements_vector(eltree_array(1,elem)%p,&
                        &lsign=[-1,1]) == 0 ) ) then

!         no interface, remove eltree

          call delete(eltree_array(1,elem)%p)  ! delete actual eltree
          deallocate( eltree_array(1,elem)%p ) ! delete pointer target
                                               ! (pointer becomes disassociated)

        end if

      end if

      if ( associated(eltree_array(1,elem)%p) ) then

!       element crosses the interface

        vol(elem) = volume ( eltree_array(1,elem)%p, lsign=[1] ) / 8

        if ( vol(elem) <= epsvol ) then
!         remove eltree
          call delete(eltree_array(1,elem)%p)  ! delete actual eltree
          deallocate( eltree_array(1,elem)%p ) ! delete pointer target
                                               ! (pointer becomes disassociated)
          print *, 'element at interface removed, elem = ', elem
        end if

        minvalvol = min ( vol(elem), minvalvol )

      end if

    end do

    print *, 'minvalvol = ', minvalvol

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

    call add_to_mesh ( mesh, elementset='elements', &
                       elements=elements(1:numelem), replace=replace )

  end subroutine create_elementset_from_eltree


! delete eltree in the elements

  subroutine delete_eltree_in_elements

    integer :: elem

!   loop all elements

    do elem = 1, size(eltree_array,2)

      if ( associated(eltree_array(1,elem)%p) ) then

        call delete(eltree_array(1,elem)%p) ! delete actual eltree
        deallocate(eltree_array(1,elem)%p) ! delete pointer target
                                           ! (pointer becomes disassociated)

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

end program xfem14
