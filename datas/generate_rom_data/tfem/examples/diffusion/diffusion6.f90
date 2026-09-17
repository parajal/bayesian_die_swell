! Steady scalar diffusion problem on a unit cube with Dirichlet boundary
! conditions. There is a spherical hole in the domain, where the scalar value
! is imposed using a weak Dirichlet boundary condition.
! XFEM is used to solve the problem.
! Constant or varying alpha coefficient using a function.

module subs6_m

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

end module subs6_m


! functions module for diffusion5

module functions_df6_m

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

end module functions_df6_m

! module for the user gauss routines

module usergauss_m

  use tfem_elem_m
  use functions_df6_m
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
      gausssub%inttype = 3 ! use numerical Gauss values

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

end module usergauss_m



! the actual program

program diffusion6

  use tfem_m
  use hsl_ma57_m
  use hsl_ma41_m
  use convection_diffusion_elements_m
  use convection_diffusion_functions_m
  use subs6_m
  use poisson_functions_m
  use usergauss_m
  use figplot_m
  use io_utils_m
  use limits_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3x3 integration of hexahedrons
    gaussb = 3,         & ! 3x3 point integration of boundary elements (quads)
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
    nx=10,              & ! number of elements in x
    ny=10,              & ! number of elements in y
    nz=10,              & ! number of elements in z
    funcnr=10,          & ! function number for the right-hand side
    funcnr_alpha=2,     & ! function number for the alpha coefficient
    coeff_alpha=0         ! how to compute alpha coefficient:
                          ! 0: constant
                          ! 1: given by a function

  real(dp), parameter :: &
    lx = 1._dp,        & ! width of domain
    ly = 1._dp,        & ! depth of domain
    lz = 1._dp,        & ! height of domain
    alpha = 1._dp,     & ! diffusion coefficient
    KN = 55._dp,       & ! Nitsche coefficient
    radius = lx/4,     & ! radius of the sphere
    dx = 1e-1,         & ! displacement of the sphere relative to size of
    dy = 1e-2,         & ! an element
    dz = 1e-2,         & !
                         ! center position of sphere
    center(3) = [ lx/2+dx*lx/nx, ly/2+dy*ly/ny, lz/2+dz*lz/nz ], &
    kappa = 1.0_dp,     & ! embedded Dirichlet scalar parameter
    c = 1._dp, & ! imposed scalar (cbar) on the sphere if
                 ! coefficients%i(42) = 0
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
    epsdef = 0.0_dp,  & ! deformation of the mesh. Note, that epsdef is
                        ! relative to the element size. The deformation is
                        ! maximum near the center of the domain.
    split_threshold = 1.e-2_dp    ! levelset value for being "close" enough
                                  ! to the interface for tree splitting

  logical, parameter :: Nitsche = .false., transposed = .false., &
    symmetric = .false.

  logical, parameter :: full_figplot_mesh = .false.

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


  integer :: i, nnodes, plotobj, nbx, nby, nbz

  integer, dimension(:), allocatable :: nodes
  real(dp), dimension(:,:), allocatable :: coor

! array of pointers to an eltree

  type(eltree_p), target, dimension(:,:), allocatable :: eltree_array

! set some parameters in modules

! module limits_m

  !SUBDIVIDE_HEX6 = .true.  ! if SUBDIVIDE_HEX6=.true., the interface is
                           ! continuous, but it is more expensive.

! module functions_m

  rpl = radius ! radius of the spherical object
  cpl = center ! initial position of the center of the object

! module usergauss_m

  intrule_e = gausse  ! Gauss integration for the eltree subelements
  intrule_sub = gausssub  ! Gauss integration for the eltree subelements
  intrule_sub_small = gausssub_small  ! Gauss integration for the eltree
                                      ! subelements (small integration volumes)
  vol_small = epsvol_small  ! elements with integration volumes smaller
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
  coefficients%i(40) = 3 ! use numerical Gauss points
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

  call set_mesh_options ( meshgen_options, elshape=14, nx=nx, ny=ny, nz=nz, &
    lx=lx, ly=ly, lz=lz, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

! deform the mesh a little bit to avoid zero or very small integration volumes

  mesh%coor(:,1) = mesh%coor(:,1) * ( 1 - &
                        epsdef*(1-2*abs(mesh%coor(:,1))/lx-0.5_dp) / nx )
  mesh%coor(:,2) = mesh%coor(:,2) * ( 1 - &
                        epsdef*(1-2*abs(mesh%coor(:,2))/ly-0.5_dp) / ny )
  mesh%coor(:,3) = mesh%coor(:,3) * ( 1 - &
                        epsdef*(1-2*abs(mesh%coor(:,3))/lz-0.5_dp) / nz )

! create object for particle boundary (only for sampling)

  nnodes = 8 * nx * radius / lx

  allocate ( coor(nnodes,3) )

  rp = radius ! radius of the spherical object
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
  nbz = nint(sqrt(real(nz)))

  call add_to_mesh ( mesh, blocks=[nbx,nby,nbz] )

  call fill_mesh_parts ( mesh )

  plot_options%viewpoint = [ 1._dp, 1.7_dp, 0.5_dp ]
  call plot_points_curves ( plot_options, mesh, 'curves_df6.fig' )

  if ( full_figplot_mesh ) then

    call plot_mesh ( plot_options, mesh, 'mesh_df6.fig' )
    plot_options%objectpointcolor=4
    plot_options%objectpointsize=0.4
    call plot_objects ( plot_options, mesh, 'mesh_df6.fig', append=.true., &
      object1=1 )
    plot_options%objectpointcolor=5
    call plot_objects ( plot_options, mesh, 'mesh_df6.fig', append=.true.,  &
      object1=2 )

  else

    call plot_mesh ( plot_options, mesh, 'mesh_df6.fig', surfaces=[3,4,6] )

  end if

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a(:,1) = 1  ! for alpha in the nodes

  call define_essential ( mesh, input_probdef, surface1=1, surface2=6 )
  call define_essential ( mesh, input_probdef, nodeset1=1 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    surface1=1, surface2=6, func=func, funcnr=9 )
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

  open( unit=15, file='sample_c_df6.out', recl=300 )

! sample scalar on object

  call fill_sample ( mesh, problem, sample_c, ndegfd=1, object=1, &
    elemsub=poisson_sample, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_c%nnodes
    write( unit=15, fmt=* ) sample_c%coor(i,:), sample_c%u(i,:)
  end do

  close( unit=15 )

  open( unit=15, file='sample_gradc_df6.out', recl=300 )

! sample gradient on object

  call fill_sample ( mesh, problem, sample_gradc, ndegfd=3, object=1, &
    elemsub=poisson_sample_deriv, coefficients=coefficients, &
    oldvectors=oldvectors )

  do i = 1, sample_gradc%nnodes
    write( unit=15, fmt=* ) sample_gradc%coor(i,:), sample_gradc%u(i,:)
  end do

  close( unit=15 )


! post-processing: create mesh_plot for plotting

  call mesh_convert ( mesh, mesh_plot, userelmesh=userelmesh, warn=.false. )

  call fill_mesh_parts ( mesh_plot )

! problem definition for post_processing

  call create_input_probdef ( mesh_plot, input_probdef_plot, nvec=2 )

  do i = 1, mesh_plot%nelgrp
    input_probdef_plot%elementdof(i)%a(:) = 0
    input_probdef_plot%vec_elementdof(i)%a(:,1) = 1  ! scalar
    input_probdef_plot%vec_elementdof(i)%a(:,2) = 3  ! gradient of scalar
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

  call fill_sample ( mesh, problem, sample_gradc, ndegfd=3, object=plotobj, &
    elemsub=poisson_sample_deriv, coefficients=coefficients, &
    oldvectors=oldvectors )

  gradscalar%u = reshape ( transpose(sample_gradc%u), [3*mesh_plot%nnodes] )

  call write_scalar_vtk ( mesh_plot, problem_plot, filename='diffusion6.vtk', &
    dataname='scalar', vector=scalar )

  call write_vector_vtk ( mesh_plot, problem_plot, filename='diffusion6.vtk', &
    dataname='grad_scalar', vector=gradscalar, append=.true. )

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

end program diffusion6

