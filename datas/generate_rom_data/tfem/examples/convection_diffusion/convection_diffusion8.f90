! Steady scalar convection-diffusion problem on a unit square with Dirichlet
! and natural boundary conditions. There is a circular region in the domain,
! where on the boundary a jump in flux is applied (heat supply).
! This is basically an extension of diffusion7 to a convection-diffusion problem
! Weak embedded interface conditions are imposed at the cylindrical interface.
! XFEM is used to solve the problem.
! Constant or varying alpha coefficient using a function.
! Varying velocity coefficient using a function.
! Optional different weighting for both sides of the interface.
! Optional Nitsche's method.
! Optional Baumann-Oden method.
! Optional Lesaint-Raviart upwinding along the interface.


module subs8_m

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

end module subs8_m


! functions module for convection_diffusion8

module functions_cd8_m

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

end module functions_cd8_m

! module for the user gauss routines

module usergauss_m

  use tfem_elem_m
  use functions_cd8_m
  use eltree_m

  implicit none

  save

! integration area near interface elements (on the positive side)
  real(dp), dimension(:), allocatable :: vol

! indicates if an element is outside (.true.) or inside (.false.) for elements
! not at the interface. For elements at the interface it is undefined.
  logical, dimension(:), allocatable :: outside

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
  logical :: printint = .false.

! array of pointers to an eltree (for userelmesh only)

  type(eltree_p), dimension(:,:), pointer :: ea

! treshold for small volumes

  real(dp) :: vol_small = 1e-2_dp

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

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors

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

    if ( associated(eltree1) ) then

!     element crosses the interface

      if ( intregion == 1 ) then

!       outside region

!       composite integration scheme for levelset > 0

        if ( vol(elem) > vol_small ) then
          ninti = number_of_integration_points ( eltree1, lsign=[1],&
            ninti=ninti_e, nintis=ninti_sub, integration_options=intopt )
        else
          ninti = number_of_integration_points ( eltree1, lsign=[1],&
            ninti=ninti_e, nintis=ninti_sub_small, integration_options=intopt )
        end if

      else if ( intregion == 2 ) then

!       inside region

!       composite integration scheme for levelset < 0

        if ( 1 - vol(elem) > vol_small ) then
          ninti = number_of_integration_points ( eltree1, lsign=[-1],&
            ninti=ninti_e, nintis=ninti_sub, integration_options=intopt )
        else
          ninti = number_of_integration_points ( eltree1, lsign=[-1],&
            ninti=ninti_e, nintis=ninti_sub_small, integration_options=intopt )
        end if

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

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors

!   set eltree pointer

    eltree1 => oldvectors%ea(1)%p(elgrp,elem)%p

    if ( associated(eltree1) ) then

!     element crosses the interface

      if ( intregion == 1 ) then

!       outside region

!       composite integration scheme for levelset > 0

        if ( vol(elem) > vol_small ) then
          call integration_points ( eltree1, xig, wg, lsign=[1], &
            ninti=ninti_e, xe=xe, we=we, nintis=ninti_sub, xs=xsub, ws=wsub, &
            integration_options=intopt )
        else
          call integration_points ( eltree1, xig, wg, lsign=[1], &
            ninti=ninti_e, xe=xe, we=we, nintis=ninti_sub_small, &
            xs=xsub_small, ws=wsub_small, integration_options=intopt )
        end if

      else if ( intregion == 2 ) then

!       inside region

!       composite integration scheme for levelset < 0

        if ( 1 - vol(elem) > vol_small ) then
          call integration_points ( eltree1, xig, wg, lsign=[-1], &
            ninti=ninti_e, xe=xe, we=we, nintis=ninti_sub, xs=xsub, ws=wsub, &
            integration_options=intopt )
        else
          call integration_points ( eltree1, xig, wg, lsign=[-1], &
            ninti=ninti_e, xe=xe, we=we, nintis=ninti_sub_small, &
            xs=xsub_small, ws=wsub_small, integration_options=intopt )
        end if

      end if

    else

!     fully in or out: standard integration

      call set_Gauss_integration ( gauss, xig, wg )

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

    type(eltree_t), pointer :: eltree


!   set eltree pointer for convenience

    eltree => ea(elgrp,elem)%p

    if ( associated(eltree) ) then

!     element crosses the interface

      lelgrp = elgrp
      lelem = elem

!     add only elements in the outside or inside region

      if ( meshregion == 1 ) then
!       outside region
        call eltree_to_mesh ( eltree, elmesh, lsign=[1], mapcoor=mapcoor )
      else if ( meshregion == 2 ) then
!       inside region
        call eltree_to_mesh ( eltree, elmesh, lsign=[-1], mapcoor=mapcoor )
      end if

      indicator = 2 ! replace element with submesh

    else if ( outside(elem) ) then

!     fully outside

      if ( meshregion == 1 ) then
!       outside region
        indicator = 0  ! leave element as it is
      else if ( meshregion == 2 ) then
!       inside region
        indicator = 1  ! remove element
      end if

    else if ( .not. outside(elem) ) then

!     fully inside

      if ( meshregion == 1 ) then
!       outside region
        indicator = 1  ! remove element
      else if ( meshregion == 2 ) then
!       inside region
        indicator = 0  ! leave element as it is
      end if

    end if

  end subroutine userelmesh

end module usergauss_m



! the actual program

program convection_diffusion8

  use tfem_m
  use hsl_ma41_m
  use convection_diffusion_elements_m
  use convection_diffusion_functions_m
  use subs8_m
  use poisson_functions_m
  use usergauss_m
  use figplot_m
  use io_utils_m

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
    funcnr_alpha1=1,    & ! function number for the alpha coefficient (outside)
    funcnr_alpha2=1,    & ! function number for the alpha coefficient (inside)
    coeff_alpha=0,      & ! how to compute alpha coefficient:
                          ! 0: constant
                          ! 1: given by a function
    vfuncnr_velo=1,     & ! function number for the velocity coefficient
    coeff_velo=0,       & ! how to compute velocity coefficient:
                          ! 0: constant
                          ! 1: given by a function
                          ! 2: given by nodal values in a vector
    omega_weights=0       ! how to compute the weights (omegaw):
                          ! 0: equally (omegaw=0.5)
                          ! 1: fully on the first layer (omegaw=0)
                          ! 2: fully on the second layer (omegaw=1)
                          ! 3: non-equally


  real(dp), parameter :: &
    lx = 1._dp,        & ! width of domain
    ly = 1._dp,        & ! height of domain
    omega_par = 0.5_dp, & ! weight factor omegaw in case of omega_weights=3
    gamma = 1._dp, &     ! coefficient of the time derivative "rho c_p"
    alpha1 = 1.0_dp,    & ! diffusion coefficient (outside)
    alpha2 = 1.0_dp,    & ! diffusion coefficient (inside)
    velocity(2)=[10._dp,5._dp], & ! constant velocity for coeff_velo=0
    KN = 55._dp,       & ! Nitsche coefficient
    radius = lx/4,     & ! radius of the cylinder
    dx = 1e-1, dy = 1e-2, & ! displacement of the cylinder relative to size of
                            !! an element
    center(2) = [ lx/2+dx*lx/nx, ly/2+dy*ly/ny ], & ! center position of
                                                      ! cylinder
    kappa1 = 1.0_dp,   & ! embedded Dirichlet scalar parameter (outside)
    kappa2 = 1.0_dp,   & ! embedded Dirichlet scalar parameter (inside)
    fjump = 100._dp, & ! imposed scalar flux jump (jbar) on the cylinder if
                       ! coefficients%i(47) = 0
    epsvol = 1e-5_dp, & ! elements with integration area smaller are removed
                        ! from the eltree_array. This means that:
                        ! 1) these elements are "fully inside/outside" and no
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
    embedded = .true., Baumann_Oden = .false., Lesaint_Raviart = .true.

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t), target :: mesh
  type(mesh_t) :: mesh_plot, mesh_plot1, mesh_plot2
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
  type(solver_options_ma41_t) :: solver_options


  integer :: i, nnodes, nnodes1, nnodes2, plotobj1, plotobj2, nbx, nby

  integer, dimension(:), allocatable :: nodes1, nodes2
  real(dp), dimension(:,:), allocatable :: coor

  real(dp) :: omegaw

! plot nodesets
  logical, parameter :: plot_nodesets = .true.

! array of pointers to an eltree

  type(eltree_p), target, dimension(:,:), allocatable :: eltree_array

! the levelset function in all the nodes (it is a signed distance function
! in this case of a simple cylinder).
  real(dp), dimension(:), allocatable :: d


! set omegaw for the weighting between the two sides of the interface

  select case (omega_weights)
  case(0)
    omegaw = 0.5_dp
  case(1)
    omegaw = 0
  case(2)
    omegaw = 1
  case(3)
    omegaw = omega_par
  case default
    write(*,'(/a,i0/)') 'Error: wrong value omega_weights: ', omega_weights
    stop
  end select

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
  coefficients%i(4) = 2 ! Gauss points defined by user subroutines
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]
  coefficients%i(24) = coeff_velo
  coefficients%i(25) = vfuncnr_velo
  coefficients%i(41) = gauss_ie ! Gauss points for interface elements
  coefficients%i(45) = -1  ! -1: impose zero ibar
  coefficients%i(47) = 0   ! -1: impose zero jbar, 0: impose non-zero
                           ! jbar fjump

  coefficients%r = 0

  coefficients%r(4:5) = velocity
  coefficients%r(8) = gamma
  coefficients%r(22) = fjump


  coefficients%func => func
  coefficients%vfunc => vfunc


  coefficients%vfunc1(1)%p => velocityvfunc   ! vector function for velocity
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

  allocate ( eltree_array(1,mesh%grpnumel(1)) )
  allocate ( vol(mesh%grpnumel(1)), outside(mesh%grpnumel(1)) )

! nodes for layers, levelset

  allocate ( nodes1(mesh%nnodes), nodes2(mesh%nnodes), d(mesh%nnodes) )

  d = levelset ( mesh%coor ) ! set levelset for all nodes

! eltree for elements crossing the interface

  call create_eltree_and_nodes ( numsplit )

! elementset for elements crossing the interface

  call create_elementset_from_eltree

! nodesets for the layers

  call add_to_mesh ( mesh, nodeset='nodes', nodes=nodes1(1:nnodes1) )
  call add_to_mesh ( mesh, nodeset='nodes', nodes=nodes2(1:nnodes2) )

! objects for plotting the nodesets

  call add_to_mesh ( mesh, object='coordinates', &
    coor=mesh%coor(nodes1(1:nnodes1),:) )
  call add_to_mesh ( mesh, object='coordinates', &
    coor=mesh%coor(nodes2(1:nnodes2),:) )

  nbx = nint(sqrt(real(nx)))
  nby = nint(sqrt(real(ny)))

  call add_to_mesh ( mesh, blocks=[nbx,nby] )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves_cd8.fig' )

  if ( plot_nodesets ) then

    call plot_mesh_nodeset ( nodeset=1, filename='mesh_outside_cd8.fig' )
    call plot_mesh_nodeset ( nodeset=2, filename='mesh_inside_cd8.fig' )

  else

    call plot_mesh ( plot_options, mesh, 'mesh_cd8.fig' )

    plot_options%objectpointcolor=4
    plot_options%objectpointsize=0.4
    call plot_objects ( plot_options, mesh, 'mesh_cd8.fig', append=.true., &
      object1=1 )

  end if

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1, numlayers=2 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a(:,1) = 1  ! for alpha in the nodes

  input_probdef%layers = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, layer=1 )
  call define_essential ( mesh, input_probdef, curve1=4, layer=1 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, func=func, funcnr=7, layer=1 )

! create the structure oldvectors

  call create ( oldvectors, nsysvec=1, nelta=1 )

  oldvectors%ea(1)%p => eltree_array
  ea => eltree_array ! userelmesh does not have oldvectors as an argument

  if ( coeff_velo == 2 ) then

!   fill velocity in a vector

    write(*,'(a/)') ' Velocity in the nodes not yet available.'
    stop

  end if

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )


! build (assemble) matrix and vector from elements


! Build from the outside (layer 1)
! ================================

  intregion = 1 ! outside

  coefficients%i(2) = coeff_alpha
  coefficients%i(3) = funcnr_alpha1
  coefficients%i(17) = -1 ! normal to fluid points opposite to the interface
  coefficients%i(49) = 1 ! outside (clayer=1)
  coefficients%r(1) = alpha1

! volume integral

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_diffusion_elem, coefficients=coefficients, &
    oldvectors=oldvectors, layer=1 )

! convection term

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_ugradc_elem, coefficients=coefficients, &
    oldvectors=oldvectors, addmatvec=.true., buildvector=.false., layer=1 )

! flux at the interface (NOTE: order='DN' must be specified if no physical
! quantities have been defined).

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_diffusion_flux_interface_eltree, &
    oldvectors=oldvectors, coefficients=coefficients, &
    factormat=1-omegaw, factorvec=omegaw, &
    addmatvec=.true., elementset=1, order='DN' )

! transposed term from the outside

  if ( transposed ) then

    coefficients%i(2) = 0  ! constant alpha value
    if ( Baumann_Oden ) then
      coefficients%i(44) = -1 ! sign of the transposed flux term
    else
      coefficients%i(44) = 1 ! sign of the transposed flux term
    end if

!   NOTE minus sign for rhs (important only if cjump /= 0).

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_diffusion_flux_interface_transposed_eltree, &
      oldvectors=oldvectors, coefficients=coefficients, &
      factormat=1-omegaw, factorvec=1-omegaw, &
      addmatvec=.true., elementset=1, order='DN' )

  end if

! stabilizing term from the outside (layer=1)

  if ( embedded ) then

    coefficients%r(19) = kappa1   ! embedded factor (outside)

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_diffusion_embedded_interface_eltree, &
      oldvectors=oldvectors, coefficients=coefficients, &
      factormat=(1-omegaw)**2, factorvec=(1-omegaw)**2, &
      addmatvec=.true., elementset=1, order='DN' )

  end if


! Build from the inside (layer 2)
! ================================

  intregion = 2 ! inside

  coefficients%i(2) = coeff_alpha
  coefficients%i(3) = funcnr_alpha2
  coefficients%i(17) = 1 ! normal to fluid points like the interface
  coefficients%i(49) = 2 ! inside (clayer=2)
  coefficients%r(1) = alpha2

! volume integral

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_diffusion_elem, coefficients=coefficients, &
    oldvectors=oldvectors, addmatvec=.true., layer=2 )

! convection term

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_ugradc_elem, coefficients=coefficients, &
    oldvectors=oldvectors, addmatvec=.true., buildvector=.false., layer=2 )

! flux at the interface

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_diffusion_flux_interface_eltree, &
    oldvectors=oldvectors, coefficients=coefficients, &
    factormat=omegaw, factorvec=1-omegaw, &
    addmatvec=.true., elementset=1, order='DN' )

! transposed term from the inside

  if ( transposed ) then

    coefficients%i(2) = 0  ! constant alpha value
    if ( Baumann_Oden ) then
      coefficients%i(44) = -1 ! sign of the transposed flux term
    else
      coefficients%i(44) = 1 ! sign of the transposed flux term
    end if

!   NOTE minus sign for rhs (important only if cjump /= 0).

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_diffusion_flux_interface_transposed_eltree, &
      oldvectors=oldvectors, coefficients=coefficients, &
      factormat=omegaw, factorvec=-omegaw, &
      addmatvec=.true., elementset=1, order='DN' )

  end if

! stabilizing term from the inside (layer=2)

  if ( embedded ) then

    coefficients%r(19) = kappa2   ! embedded factor (inside)

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_diffusion_embedded_interface_eltree, &
      oldvectors=oldvectors, coefficients=coefficients, &
      factormat=omegaw**2, factorvec=omegaw**2, &
      addmatvec=.true., elementset=1, order='DN' )

  end if


! two sides at once
! =================

  if ( Lesaint_Raviart ) then

!   Lesaint-Raviart upwinding along interface (no separate layers needed)

    coefficients%i(17) = -1 ! normal from layer 1 to layer 2

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_ugradc_interface_eltree, &
      oldvectors=oldvectors, coefficients=coefficients, &
      addmatvec=.true., elementset=1, order='DN' )

  end if

! optional Nitsche stabilization (no separate layers needed)

  if ( Nitsche ) then

    coefficients%r(20) = max(alpha1,alpha2)*KN*nx/lx ! Nitsche's factor

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_diffusion_Nitsche_interface_eltree, &
      oldvectors=oldvectors, coefficients=coefficients, &
      addmatvec=.true., elementset=1, order='DN' )

  end if

! build missing matrix elements with zero

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    coefficients=coefficients, zeromatvec=.true., &
    addmatvec=.true., exclude_single_layer =.true. )


  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call check ( sysmatrix )


  solver_options%real_storage = 2._dp

  call solve_system_ma41 ( sysmatrix, rhsd, sol, solver_options=solver_options )

! print average sol to standard output

  print *, sum( abs(sol%u) ) / sol%n

  call delete_eltree_in_elements ! remove eltree for all elements

! post-processing

  oldvectors%s(1)%p => sol

! sample scalar on cylinder

  call sample_scalar_on_cylinder ( layer=1, filename='sample_c_out_cd8.out' )
  call sample_scalar_on_cylinder ( layer=2, filename='sample_c_in_cd8.out' )

! sample gradient on cylinder

  call sample_gradient_on_cylinder ( layer=1, filename='sample_gc_out_cd8.out' )
  call sample_gradient_on_cylinder ( layer=2, filename='sample_gc_in_cd8.out' )


! create eltree for all elements crossing the interface

  call create_eltree_and_nodes ( numsplit_plot )


! post-processing: create mesh_plot for plotting

  meshregion = 1
  call mesh_convert ( mesh, mesh_plot1, userelmesh=userelmesh, warn=.false. )
  meshregion = 2
  call mesh_convert ( mesh, mesh_plot2, userelmesh=userelmesh, warn=.false. )

  call mesh_merge ( mesh_plot1, mesh_plot2, mesh_plot )

  call fill_mesh_parts ( mesh_plot )

  call plot_mesh ( plot_options, mesh_plot, 'mesh_plot_cd8.fig' )


! create objects for interpolation on the original mesh

  warn_add_to_mesh_after_meshgen_parts = .false. ! .false. to suppress warning
  call add_to_mesh ( mesh, object='coordinates', coor=mesh_plot1%coor, &
    onlynodeset=1 )  ! use onlynodeset to avoid intersecting the wrong layer.
  plotobj1 = mesh%nobjects
  call add_to_mesh ( mesh, object='coordinates', coor=mesh_plot2%coor, &
    onlynodeset=2 )  ! use onlynodeset to avoid intersecting the wrong layer.
  plotobj2 = mesh%nobjects
  call fill_mesh_parts_objects ( mesh, object1=plotobj1, object2=plotobj2 )

  call delete ( mesh_plot1, mesh_plot2 )

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

! sample and plot scalar on mesh_plot

  call sample_scalar
  call plot_contour_on_mesh_plot ( filename='scalar_cd8.fig', vector=scalar)

! sample and plot gradient of scalar on mesh_plot

  call sample_gradient_scalar
  call plot_contour_on_mesh_plot ( filename='gradc_x.fig', &
    vector=gradscalar, degfd=1 )
  call plot_contour_on_mesh_plot ( filename='gradc_y.fig', &
    vector=gradscalar, degfd=2 )


! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=scalar, &
    dataname='scalar', filename='convection_diffusion8.vtk' )

  call write_vector_vtk ( mesh_plot, problem_plot, &
    filename='convection_diffusion8.vtk', &
    dataname='gradc', vector=gradscalar, append=.true. )


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

  deallocate ( coor, nodes1, nodes2, d, eltree_array, vol, outside )

contains


! Create eltree in the elements which are crossed by the interface
! Also set nodes for layer=1 (outside) and layer=2 (inside)

  subroutine create_eltree_and_nodes ( numsplit )

    integer, intent(in) :: numsplit

    integer :: nod(mesh%element(1)%numnod), numsubel(2), elem, i
    integer :: tmpnodes(mesh%nnodes)
    integer, parameter :: vertQ2(4) = [1,3,5,7]  ! vertices of quad Q2 element
    real(dp) :: coor(2,2)
    real(dp) :: minvalvol


!   loop all elements

    minvalvol = 1
    nodes1 = 0
    nodes2 = 0

    do elem = 1, mesh%grpnumel(1)

!     nodal points

      nod = mesh%topology(1)%a(:,elem)

      if ( all ( d(nod(vertQ2)) > split_threshold ) ) then

!       element not subdivided (no interface)
!       all nodes far from the interface (outside layer)

        nodes1(nod) = 1
        outside(elem) = .true.

        cycle

      else if ( all ( d(nod(vertQ2)) < -split_threshold ) ) then

!       element not subdivided (no interface)
!       all nodes far from the interface (inside layer)

        nodes2(nod) = 1
        outside(elem) = .false.

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

        numsubel = number_of_subelements_vector(eltree_array(1,elem)%p,&
                         &lsign=[-1,1])

        if ( any( numsubel == 0 ) ) then

!         no interface, remove eltree (NEEDS TO BE RECHECKED FOR EDGES!!)

          call delete(eltree_array(1,elem)%p)  ! delete actual eltree
          deallocate( eltree_array(1,elem)%p ) ! delete pointer target
                                               ! (pointer becomes disassociated)

          if ( numsubel(1) == 0 ) then

!           nodes in layer 1 (outside)

            nodes1(nod) = 1
            outside(elem) = .true.

          else

!           nodes in layer 2 (inside)

            nodes2(nod) = 1
            outside(elem) = .false.

          end if

        end if

      end if

      if ( associated(eltree_array(1,elem)%p) ) then

!       element crosses the interface

!       volume on the positive side of the interface (outside)

        vol(elem) = volume ( eltree_array(1,elem)%p, lsign=[1] ) / 4

        if ( vol(elem) <= epsvol .or. vol(elem) >= 1 - epsvol ) then

!         remove eltree

          call delete(eltree_array(1,elem)%p)  ! delete actual eltree
          deallocate( eltree_array(1,elem)%p ) ! delete pointer target
                                               ! (pointer becomes disassociated)
          print *, 'element at interface removed, elem = ', elem

!         set nodes

          if ( vol(elem) <= epsvol ) then
!           element is part of the inside (nodes in layer 2)
            nodes2(nod) = 1
            outside(elem) = .false.
          else
!           element is part of the outside (nodes in layer 1)
            nodes1(nod) = 1
            outside(elem) = .true.
          end if

        else

!         set nodes (all layers present)

          nodes1(nod) = 1
          nodes2(nod) = 1

        end if

        minvalvol = minval ( [ vol(elem), 1._dp - vol(elem), minvalvol ] )

      end if

    end do

    print *, 'minvalvol = ', minvalvol

!   count and collect the nodes

    tmpnodes = nodes1

    nnodes1 = 0
    do i = 1, mesh%nnodes
      if ( tmpnodes(i) == 0 ) cycle
      nnodes1 = nnodes1 + 1
      nodes1(nnodes1) = i
    end do

    tmpnodes = nodes2

    nnodes2 = 0
    do i = 1, mesh%nnodes
      if ( tmpnodes(i) == 0 ) cycle
      nnodes2 = nnodes2 + 1
      nodes2(nnodes2) = i
    end do

  end subroutine create_eltree_and_nodes


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


  subroutine sample_scalar_on_cylinder ( layer, filename )

    integer, intent(in) :: layer
    character(len=*), intent(in) :: filename

    integer :: i

    open( unit=15, file=filename, recl=300 )

    coefficients%i(38) = layer

!   sample on cylinder (errorobject=1 to warn for points not in layer)

    call fill_sample ( mesh, problem, sample_c, ndegfd=1, object=1, &
      layer=layer, elemsub=poisson_sample, coefficients=coefficients, &
      oldvectors=oldvectors, errorobject=1 )

    do i = 1, sample_c%nnodes
      write( unit=15, fmt=* ) sample_c%coor(i,:), sample_c%u(i,:)
    end do

    close( unit=15 )

  end subroutine sample_scalar_on_cylinder

  subroutine sample_gradient_on_cylinder ( layer, filename )

    integer, intent(in) :: layer
    character(len=*), intent(in) :: filename

    integer :: i

    open( unit=15, file=filename, recl=300 )

    coefficients%i(38) = layer

!   sample on cylinder (errorobject=1 to warn for points not in layer)

    call fill_sample ( mesh, problem, sample_gradc, ndegfd=2, object=1, &
      layer=layer, elemsub=poisson_sample_deriv, coefficients=coefficients, &
      oldvectors=oldvectors, errorobject=1 )

    do i = 1, sample_gradc%nnodes
      write( unit=15, fmt=* ) sample_gradc%coor(i,:), sample_gradc%u(i,:)
    end do

    close( unit=15 )

  end subroutine sample_gradient_on_cylinder


! sample scalar on mesh_plot

  subroutine sample_scalar

    integer :: n

!   sample outside (use layer as argument to detect points not in layer)

    coefficients%i(38) = 1

    call fill_sample ( mesh, problem, sample_c, ndegfd=1, object=plotobj1, &
      elemsub=poisson_sample, coefficients=coefficients, &
      oldvectors=oldvectors, layer=1, errorobject=1, allnodes=.true. )

    n = sample_c%nnodes

    scalar%u(1:n) = sample_c%u(1:n,1)

!   sample inside (use layer as argument to detect points not in layer)

    coefficients%i(38) = 2

    call fill_sample ( mesh, problem, sample_c, ndegfd=1, object=plotobj2, &
      elemsub=poisson_sample, coefficients=coefficients, &
      oldvectors=oldvectors, layer=2, errorobject=1, allnodes=.true. )

    scalar%u(n+1:) = sample_c%u(1:sample_c%nnodes,1)

  end subroutine sample_scalar


! sample gradient of scalar on mesh_plot

  subroutine sample_gradient_scalar

    integer :: n

!   sample outside (use layer as argument to detect points not in layer)

    coefficients%i(38) = 1

    call fill_sample ( mesh, problem, sample_gradc, ndegfd=2, object=plotobj1, &
      elemsub=poisson_sample_deriv, coefficients=coefficients, &
      oldvectors=oldvectors, layer=1, errorobject=1, allnodes=.true. )

    n = 2*sample_gradc%nnodes

    gradscalar%u(1:n) = reshape ( transpose(sample_gradc%u), [n] )

!   sample inside (use layer as argument to detect points not in layer)

    coefficients%i(38) = 2

    call fill_sample ( mesh, problem, sample_gradc, ndegfd=2, object=plotobj2, &
      elemsub=poisson_sample_deriv, coefficients=coefficients, &
      oldvectors=oldvectors, layer=2, errorobject=1, allnodes=.true. )

    gradscalar%u(n+1:) = reshape ( transpose(sample_gradc%u), &
                                              [2*sample_gradc%nnodes] )

  end subroutine sample_gradient_scalar


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

end program convection_diffusion8

