! Steady scalar diffusion problem on a unit square with Dirichlet boundary
! conditions. There is a circular hole in the domain, where the scalar value
! is imposed using a weak Dirichlet boundary condition.
! XFEM is used to solve the problem.
! Constant or varying alpha coefficient using a function.

module subs5_m

  use math_defs_m

  implicit none

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

end module subs5_m


! functions module for diffusion5

module functions_df5_m

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

end module functions_df5_m

! module for the user gauss routines

module usergauss_m

  use tfem_elem_m
  use functions_df5_m
  use eltree_m

  implicit none

  save

! the levelset function in all the nodes (it is a signed distance function
! in this case of a simple cylinder).
  real(dp), dimension(:), allocatable :: d

! integration area near interface elements
  real(dp), dimension(:), allocatable :: vol

! the standard integration scheme:
!   ninti_s: number of integration points
!   xs, ws: points and weights
! the integration scheme on subelements in the eltree:
!   intrule_e: the integration rule
!   ninti_e: number of integration points
!   xe, we: points and weights
! the integration scheme on subelements in the submesh (triangles):
!   intrule_sub: the integration rule
!   intrule_sub_small: the integration rule for small integration areas
!   ninti_sub: number of integration points
!   ninti_sub_small: number of integration points for small integration are
!   xsub, wsub: points and weights
!   xsub_small, wsub_small: points and weights for small integration areas

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

    use convection_diffusion_globals_m

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
      gausssub%inttype = 3  ! use numerical Gauss values

      call set_ninti ( gausssub, ninti_sub )

      allocate ( wsub(ninti_sub), xsub(ninti_sub,2) )

      call set_Gauss_integration ( gausssub, xsub, wsub )

!     Gauss rule for the subdivided elements in the submesh (small areas)

      gausssub_small%globalshape = 'triangle'
      gausssub_small%intrule = intrule_sub_small
      gausssub_small%inttype = 3  ! use numerical Gauss values

      call set_ninti ( gausssub_small, ninti_sub_small )

      allocate ( wsub_small(ninti_sub_small), xsub_small(ninti_sub_small,2) )

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

!     fully in or removed element due to small integration area

      ninti = 0

    end if

  end subroutine set_ninti_user


! set the composite gauss integration rule with a user subroutine

  subroutine set_Gauss_integration_user ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors )

    use convection_diffusion_globals_m

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

!     fully in or removed element due to small integration area
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

!     fully in or removed element due to small integration area

      indicator = 1

    end if

  end subroutine userelmesh

end module usergauss_m



! the actual program

program diffusion5

  use tfem_m
  use hsl_ma57_m
  use hsl_ma41_m
  use convection_diffusion_elements_m
  use convection_diffusion_functions_m
  use subs5_m
  use poisson_functions_m
  use usergauss_m
  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3 integration of quads
    gaussb = 3,         & ! 3 point integration of boundary elements
    numsplit = 3,       & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit
    numsplitmin = 1,    & ! relative size of all subelements will be at least
                          ! as small as 1/2^numsplitmin
    numsplit_plot = 2,  & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit for plotting
    gausse = 3,         & ! Gauss integration for the eltree subelements
    gausssub_small = 8, & ! Gauss integration for the submesh subelements for
                          ! small integration volumes
    gausssub = 4,       & ! Gauss integration for the submesh subelements
    gauss_ie = 3,       & ! Gauss integration for the interface elements
    nx=50,              & ! number of elements in x
    ny=50,              & ! number of elements in y
    funcnr=8,           & ! function number for the right-hand side
    funcnr_alpha=1,     & ! function number for the alpha coefficient
    coeff_alpha=0         ! how to compute alpha coefficient:
                          ! 0: constant
                          ! 1: given by a function

  real(dp), parameter :: &
    lx = 1._dp,        & ! width of domain
    ly = 1._dp,        & ! height of domain
    alpha = 1._dp,     & ! diffusion coefficient
    KN = 55._dp,       & ! Nitsche coefficient
    radius = lx/4,     & ! radius of the cylinder
    dx = 1e-1, dy = 1e-2, & ! displacement of the cylinder relative to size of
                            ! an element
    center(2) = [ lx/2+dx*lx/nx, ly/2+dy*ly/ny ], & ! center position of
                                                      ! cylinder
    kappa = 1.0_dp,     & ! embedded Dirichlet scalar parameter
    c = 1._dp, & ! imposed scalar (cbar) on the cylinder if
                 ! coefficients%i(42) = 0
    epsvol = 1e-6_dp, & ! elements with integration area smaller are removed
                        ! from the eltree_array. This means that:
                        ! 1) these elements are "fully inside" and no
                        !    "virtual/extended degrees" are generated.
                        ! 2) the interface integration is ignored and therefore
                        !    the interface has a small "hole".
                        ! Note, that epsvol is relative to the area of an
                        ! element, so epsvol=1 is a full element.
    epsvol_small = 1e-2_dp, & ! elements with integration area smaller
                              ! are integrated with integration rule
                              ! gausssub_small, otherwise with gaussub
    epsdef = 0.0_dp,  & ! deformation of the mesh. Note, that epsdef is
                        ! relative to the element size. The deformation is
                        ! maximum near the center of the domain.
    split_threshold = 1.e-2_dp    ! levelset value for being "close" enough
                                  ! to the interface for tree splitting

  logical, parameter :: Nitsche = .false., transposed = .false., &
    symmetric = .false.

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t), target :: mesh
  type(mesh_t) :: mesh_plot
  type(input_probdef_t) :: input_probdef, input_probdef_plot
  type(problem_t) :: problem, problem_plot
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors
  type(sample_t) :: sample_c, sample_gradc
  type(vector_t) :: scalar, gradscalar


  integer :: i, nnodes, plotobj, nbx, nby

  integer, dimension(:), allocatable :: nodes
  real(dp), dimension(:,:), allocatable :: coor

! array of pointers to an eltree

  type(eltree_p), target, dimension(:,:), allocatable :: eltree_array

! set some parameters in modules

! module functions_m

  rpl = radius ! radius of the circular object
  cpl = center ! initial position of the center of the object

! module usergauss_m

  intrule_e = gausse  ! Gauss integration for the eltree subelements
  intrule_sub = gausssub  ! Gauss integration for the eltree subelements
  intrule_sub_small = gausssub_small  ! Gauss integration for the eltree
                                      ! subelements (small integration areas)
  vol_small = epsvol_small  ! elements with integration areas smaller
                            ! are integrated with integration rule
                            ! gausssub_small, otherwise with gaussub

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(2) = coeff_alpha
  coefficients%i(3) = funcnr_alpha
  coefficients%i(4) = 2 ! Gauss points defined by user subroutines
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]
  coefficients%i(17) = -1 ! normal to fluid points inwards to the object
  coefficients%i(41) = gauss_ie ! Gauss points for interface elements
  coefficients%i(42) = 0   ! -1: impose zero ubar, 0: impose non-zero
                           ! ubar c

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0
  coefficients%r(18) = c
  coefficients%r(19) = kappa   ! embedded factor
  coefficients%r(20) = alpha*KN*nx/lx ! Nitsche's factor


  coefficients%func => func

  coefficients%func1(1)%p => alphafunc   ! scalar function for alpha

  coefficients%set_ninti_user => set_ninti_user
  coefficients%set_Gauss_integration_user => set_Gauss_integration_user


! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%lx = lx
  meshgen_options%ly = ly

  call quadrilateral2d ( mesh, meshgen_options )

! deform the mesh a little bit to avoid zero or very small integration areas

  mesh%coor(:,1) = mesh%coor(:,1) * ( 1 - &
                        epsdef*(1-2*abs(mesh%coor(:,1))/lx-0.5_dp) / nx )
  mesh%coor(:,2) = mesh%coor(:,2) * ( 1 - &
                        epsdef*(1-2*abs(mesh%coor(:,2))/ly-0.5_dp) / ny )

! create object for particle boundary (only for plotting)

  nnodes = 8 * nx * radius / lx

  allocate ( coor(nnodes,2) )

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

! eltree for elements crossing the interface

  call create_eltree_in_elements ( numsplit )

! elementset for elements crossing the interface

  call create_elementset_from_eltree

  call find_internal_zero_nodes

  call add_to_mesh ( mesh, nodeset='nodes', nodes=nodes(1:nnodes) )

! object for plotting the nodeset

  call add_to_mesh ( mesh, object='coordinates', &
    coor=mesh%coor(nodes(1:nnodes),:) )

  nbx = nint(sqrt(real(nx)))
  nby = nint(sqrt(real(ny)))

  call add_to_mesh ( mesh, blocks=[nbx,nby] )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves_df5.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh_df5.fig' )
  plot_options%objectpointcolor=4
  plot_options%objectpointsize=0.4
  call plot_objects ( plot_options, mesh, 'mesh_df5.fig', append=.true., &
    object1=1 )
  plot_options%objectpointcolor=5
  call plot_objects ( plot_options, mesh, 'mesh_df5.fig', append=.true.,  &
    object1=2 )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a(:,1) = 1  ! for alpha in the nodes

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4 )
  call define_essential ( mesh, input_probdef, nodeset1=1 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, func=func, funcnr=7 )
  call fill_sysvector ( mesh, problem, sol, &
    nodeset1=1, value=0._dp )

! create the structure oldvectors

  call create ( oldvectors, nsysvec=1, nelta=1 )

  oldvectors%ea(1)%p => eltree_array
  ea => eltree_array ! userelmesh does not have oldvectors as an argument

! create system matrix

  if ( symmetric ) then
    call create_sysmatrix_structure ( sysmatrix, mesh, problem, &
      symmetric=.true. )
  else
    call create_sysmatrix_structure ( sysmatrix, mesh, problem )
  end if

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_diffusion_elem, coefficients=coefficients, &
    oldvectors=oldvectors )

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_diffusion_open_boundary_eltree, oldvectors=oldvectors, &
    coefficients=coefficients, addmatvec=.true., elementset=1 )

  if ( Nitsche ) then

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_diffusion_Nitsche_eltree, oldvectors=oldvectors, &
      coefficients=coefficients, addmatvec=.true., elementset=1 )

  else

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_diffusion_embedded_dirichlet_eltree, &
      oldvectors=oldvectors, coefficients=coefficients, &
      addmatvec=.true., elementset=1 )

  end if

  if ( transposed ) then

    coefficients%i(44) = 1 ! sign of the transposed flux term
    coefficients%i(2) = 0  ! constant alpha value

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_diffusion_open_boundary_transposed_eltree, &
      oldvectors=oldvectors, coefficients=coefficients, &
      addmatvec=.true., elementset=1 )

  end if

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  if ( symmetric ) then
    call solve_system_ma57 ( sysmatrix, rhsd, sol )
  else
    call solve_system_ma41 ( sysmatrix, rhsd, sol )
  end if

! print average sol to standard output

  print *, sum( abs(sol%u) ) / sol%n

  call delete_eltree_in_elements ! remove eltree for all elements

! post-processing

! create eltree for all elements crossing the interface

  call create_eltree_in_elements ( numsplit_plot )

  oldvectors%s(1)%p => sol

! sampling for gnuplot

  open( unit=15, file='sample_c_df5.out', recl=300 )

! sample scalar on object

  call fill_sample ( mesh, problem, sample_c, ndegfd=1, object=1, &
    elemsub=poisson_sample, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_c%nnodes
    write( unit=15, fmt=* ) sample_c%coor(i,:), sample_c%u(i,:)
  end do

  close( unit=15 )

  open( unit=15, file='sample_gradc_df5.out', recl=300 )

! sample gradient on object

  call fill_sample ( mesh, problem, sample_gradc, ndegfd=2, object=1, &
    elemsub=poisson_sample_deriv, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_gradc%nnodes
    write( unit=15, fmt=* ) sample_gradc%coor(i,:), sample_gradc%u(i,:)
  end do

  close( unit=15 )


! post-processing: create mesh_plot for plotting

  call mesh_convert ( mesh, mesh_plot, userelmesh=userelmesh, warn=.false. )

  call fill_mesh_parts ( mesh_plot )

  call plot_mesh ( plot_options, mesh_plot, 'mesh_df5_plot.fig' )

! problem definition for post_processing

  call create_input_probdef ( mesh_plot, input_probdef_plot, nvec=2 )

  do i = 1, mesh_plot%nelgrp
    input_probdef_plot%elementdof(i)%a(:) = 0
    input_probdef_plot%vec_elementdof(i)%a(:,1) = 1  ! scalar
    input_probdef_plot%vec_elementdof(i)%a(:,2) = 2  ! gradient of scalar
  end do

  call problem_definition ( input_probdef_plot, mesh_plot, problem_plot )

  call create_vector ( problem_plot, scalar, vec=1 )
  call create_vector ( problem_plot, gradscalar, vec=2 )

! create an object for interpolation on the original mesh

  warn_add_to_mesh_after_meshgen_parts = .false. ! .false. to suppress warning
  call add_to_mesh ( mesh, object='coordinates', coor=mesh_plot%coor )
  plotobj = mesh%nobjects
  call fill_mesh_parts_objects ( mesh, object1=plotobj )

! sample scalar

  call fill_sample ( mesh, problem, sample_c, ndegfd=1, object=plotobj, &
    elemsub=poisson_sample, coefficients=coefficients, &
    oldvectors=oldvectors )

  scalar%u = sample_c%u(:,1)

! sample gradient of scalar

  call fill_sample ( mesh, problem, sample_gradc, ndegfd=2, object=plotobj, &
    elemsub=poisson_sample_deriv, coefficients=coefficients, &
    oldvectors=oldvectors )

  gradscalar%u = reshape ( transpose(sample_gradc%u), [2*mesh_plot%nnodes] )

  plot_options%plotboundary2=.false.

  call plot_color_contour ( plot_options, mesh_plot, problem_plot, &
    'scalar_plot.fig', vector=scalar )
  call plot_boundary ( plot_options, mesh, 'scalar_plot.fig', append=.true. )

  call plot_color_contour ( plot_options, mesh_plot, problem_plot, &
    'gradxscalar_plot.fig', vector=gradscalar, degfd=1 )
  call plot_boundary ( plot_options, mesh, 'gradxscalar_plot.fig', &
     append=.true. )

  call plot_color_contour ( plot_options, mesh_plot, problem_plot, &
    'gradyscalar_plot.fig', vector=gradscalar, degfd=2 )
  call plot_boundary ( plot_options, mesh, 'gradyscalar_plot.fig', &
     append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( sample_c, sample_gradc )
  call delete ( scalar, gradscalar )
  call delete ( mesh_plot )
  call delete ( input_probdef_plot )
  call delete ( problem_plot )

  deallocate ( coor, nodes, d, eltree_array, vol )

contains


! create eltree in the elements which are crossed by the interface

  subroutine create_eltree_in_elements ( numsplit )

    integer, intent(in) :: numsplit

    integer :: nod(mesh%element(1)%numnod), elem
    real(dp) :: coor(2,2)
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

        vol(elem) = volume ( eltree_array(1,elem)%p, lsign=[1] ) / 4

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
                                           !(pointer becomes disassociated)

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

end program diffusion5

