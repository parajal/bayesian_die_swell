

!  mesh build from 3 submeshes:
!
!                  --------------
!                  |            |
!              n4  |            |
!                  |     3      |
!                S |            |
!                  |            |
!         -----------------------
!         |        |            |
!         |        |            |
!         |        |            |
!     n3  |        |            |
!         |   1    |     2      |
!       R |        |            |
!         |        |            |
!         |        |            |
!         |        |            |
!         |        |            |
!         -----------------------
!        -L1       0            L2
!             n1         n2
!

! NOTE: one single group is created. Use nogroupmerge=.true. in heading of
! mesh_merge to obtain 2 separate element groups.

module mesh_extrudate_swell_xfem

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use meshgen_extra_m
  !use figplot_m
  use io_utils_m

  implicit none


  type(mesh_t), target, save :: mesh

contains

  subroutine mesh_xfem

  type(meshgen_options_t) :: mesh_options
  type(mesh_t) :: mesh1, mesh2, mesh3

  integer, parameter :: elshape=6

!  real(dp), parameter :: L1=5, L2=5
!  real(dp), parameter :: R=1, S=0.5
!  integer, parameter :: n1=8, n2=8, n3=12, n4=5

!  real(dp), parameter :: f1=7, f2=6, f3=8, f4=7

  real(dp), parameter :: L1=2, L2=3
  real(dp), parameter :: R=1, S=0.5
  integer, parameter :: n1=2, n2=3, n3=2, n4=3

  real(dp), parameter :: f1=30, f2=30, f3=30, f4=10

! 1

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=n3, ox=-L1, &
    lx=L1, ly=R, ratio=[3,3,1,1], factor=[f1,f3,f1,f3] )

  call quadrilateral2d ( mesh1, mesh_options )
  !call plot_points_curves ( plot_options, mesh1, 'curves1.fig' )

! 2

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=n3, &
    lx=L2, ly=R, ratio=[1,3,3,1], factor=[f2,f3,f2,f3] )

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, curve1=2, curve2=-4 )

  !plot_options%fontsize=8
  !call plot_points_curves ( plot_options, mesh3, 'curves2.fig' )
  call delete ( mesh1, mesh2 )

! 3

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=n4, &
    lx=L2, ly=S, oy=R, ratio=[1,1,3,3], factor=[f2,f4,f2,f4] )

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh3, mesh1, mesh, curve1=7, curve2=-1 )

  call delete ( mesh1, mesh3 )


! add curves

! outflow
  call add_to_mesh ( mesh, curve=[6,8] ) ! curve 11
! centerline
  call add_to_mesh ( mesh, curve=[1,5] ) ! curve 12
! full boundary
  call add_to_mesh ( mesh, curve=[1,5,6,8,9,10,3,4] ) ! curve 13


  !call fill_mesh_parts ( mesh )

  !call plot_points_curves ( plot_options, mesh, 'curves_esx.fig' )
  !call plot_mesh ( plot_options, mesh, 'mesh_esx.fig' )

!  plot_options%xmin = -1
!  plot_options%xmax = 1
!  call plot_mesh ( plot_options, mesh, 'mesh_esx_zoom.fig' )

  !call write_mesh ( mesh, 'mesh_esx.out' )

  !call printinfo ( mesh, printlevel=2 )

  !call delete ( mesh )

  end subroutine mesh_xfem

end module mesh_extrudate_swell_xfem



! module with additional functions for the extrudate swell problem with xfem

module functions_xfem_m

  use tfem_elem_m

  implicit none
  save

! local pointers to data in main program
  type(mesh_t), pointer :: lmesh_sf_adv => null()
  type(problem_t), pointer :: lproblem_sf_adv => null()
  type(sysvector_t), pointer :: lsol_sf_adv_pred => null()


! initial height
  real(dp) :: h0

! mesh used for mapping of reference element
  type(mesh_t), pointer :: lmesh => null()

  integer :: lelem, lelgrp

contains

! levelset function for the interface

  function levelset ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1)) :: levelset

    integer :: i, elem, nod(2)

    real(dp) :: x0, x1, x2, y0, y1, y2, yp, xtmp(2), u(2)
    logical :: find


    do i = 1, size(x,1)

      x0 = x(i,1)
      y0 = x(i,2)

      if ( x0 <= 0._dp ) then

!       inside the die

        levelset(i) = h0 - y0

      else

!       outside the die

        find  = .false.

        do elem = 1, lmesh_sf_adv%nelem

          nod = lmesh_sf_adv%topology(1)%a(1:2,elem)
          xtmp = lmesh_sf_adv%coor(nod,1)

          x1 = xtmp(1); x2 = xtmp(2)

          if ( (x0 >= x1) .and. (x0 <= x2) ) then

!           get surface height for element elem

            u = lsol_sf_adv_pred%u( &
              lproblem_sf_adv%degfdperm(lproblem_sf_adv%nodnumdegfd(nod)+1,2) )
            y1 = u(1); y2 = u(2)

!           height at position x0
            yp = y1 + (y2-y1)*(x0-x1)/(x2-x1)

            levelset(i) = yp - y0

            find = .true.

            exit

          end if

        end do

        if ( .not. find ) then
          print*, 'Error: Finding levelset failed.'
          print*, i, x0, y0
          print*, 'Program Stop!'
          stop
        end if

      end if

    end do


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


! function for the incremental vertical displacement of the interface
! Used for computing the temporary ALE mesh (md=mesh displacement).

  function func_md ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func_md

    real(dp) :: disp(1), coor(1,size(x))

    if ( x(1) <= 0 ) then
      func_md = 0._dp
    else
      coor(1,:) = x
      disp = levelset( coor )
      func_md = disp(1)
    end if

  end function func_md


end module functions_xfem_m
! module with additional routines for the extrudate swell problem with xfem

module usergauss_xfem_m

  use tfem_elem_m
  use functions_xfem_m

  implicit none
  save

! the levelset function in all the nodes
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
!   ninti_sub_small: number of integration points for small integration areas
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

!     Gauss rule for the subdivided elements in the submesh (small volumes)

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

!     fully inside fluid

      ninti = ninti_s

    else

!     fully outside fluid or removed element due to small integration area

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

!     fully inside fluid: standard integration

      call set_Gauss_integration ( gauss, xig, wg )

    else

!     fully outside fluid or removed element due to small integration area
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

!     fully inside fluid: leave element as it is.

      indicator = 0

    else

!     fully outside fluid or removed element due to small integration area

      indicator = 1

    end if

  end subroutine userelmesh


end module usergauss_xfem_m

module subs_extrudate_swell_m

  use tfem_elem_m

  implicit none

contains


! Special element routine for the surface tension force at the boundary point
! Used together with add_boundary_elements_point

  subroutine surface_tension_boundary_point1 ( mesh, problem, point, matrix, &
    vector, coefficients, oldvectors, elemmat, elemvec )

    use math_defs_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: point
    logical, intent(in) :: matrix, vector
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer, parameter :: ndim = 2
    integer, parameter :: ndf = 3
    integer, parameter :: nodalp = 3
    integer, parameter :: ninti = 1

    real(dp) :: curvel(ninti), normal(ninti,ndim)
    real(dp) :: xig(ninti,1), x(nodalp,ndim)
    real(dp) :: phi(ninti,ndf), dphi(ninti,ndf,1), dxdxi(ninti,ndim)
    real(dp) :: normalc(ndim)

    integer :: curve, elem, coorsys
    real(dp) :: gammac


    if ( matrix ) then
      write(*,'(/3(a/))') 'Error in surface_tension_boundary_point1:', &
        'No matrix to build.', &
        'Call build_system with buildmatrix=.false.'
      stop
      elemmat = 0._dp
    end if

    xig(1,1) = -1._dp  ! assume first point

    ! surface tension
    curve = coefficients%i(501)
    elem = coefficients%i(502)
    coorsys = coefficients%i(23)
    gammac = coefficients%r(19)

    call shape_line_P2 ( xig, phi, dphi )

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )

    call isoparametric_deformation_curve ( x, dphi(:,:,1), dxdxi, curvel, &
      normal )

!   outside normal to the outside "curve"

    normalc(1) =  normal(1,2)
    normalc(2) = -normal(1,1)

    elemvec = normalc * gammac

    if ( coorsys == 1 ) then
!     axisymmetric
      elemvec = 2*pi*x(1,2)*elemvec
    end if

  end subroutine surface_tension_boundary_point1


! Special element routine for the surface tension force at the boundary point
! Used together with buildsystem and object/onobjectnodes

  subroutine surface_tension_boundary_point2 ( mesh, problem, eleminfo, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(eleminfo_t), intent(in) :: eleminfo
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer, parameter :: ndf = 9
    integer, parameter :: ninti = 1
    integer, parameter :: ndim = 2

    integer :: object

    real(dp), dimension(ndim) :: x1, x2, normalc
    real(dp) :: phi(ninti,ndf), xi(1,ndim)
    real(dp) :: d, gammac

!   test

    if ( matrix ) then
      write(*,'(/3(a/))') 'Error in surface_tension_boundary_point2:', &
        'No matrix to build.', &
        'Call build_system with buildmatrix=.false.'
      stop
      elemmat = 0._dp
    end if

!   set object

    object = eleminfo%object

!   set coordinates (input)

    x1 = mesh%objects(object)%coor(1,:)
    x2 = coefficients%r(501:502)

!   shape function of fluid element

    xi(1,:) = mesh%objects(object)%refcoor(1,:)

    call shape_quad_Q2 ( xi, phi )

!   outside normal to the outside "curve"

    d = sqrt ( dot_product ( x1-x2, x1-x2 ) )

    normalc = ( x1 - x2 ) / d

!   surface tension coefficient

    gammac = coefficients%r(19)

!   element vector

    elemvec(    1:ndf  ) = gammac * normalc(1) * phi(1,:)
    elemvec(ndf+1:2*ndf) = gammac * normalc(2) * phi(1,:)

  end subroutine surface_tension_boundary_point2

end module subs_extrudate_swell_m

! Viscoelastic extrudate swell problem
! Planar flow for a UCM/Oldroyd-B model
! DEVSS-G/SUPG
! XFEM approach
! Implicit bilinear terms of CE in momentum balance.
! second-order Gear time integration
! Optionally include surface tension
! NOTE: restart not yet implemented

program extrudate_swell5

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use io_utils_m
  use surface_advection_elements_m
  use usergauss_xfem_m
  use subs_extrudate_swell_m  ! subs for extrudate swell
  !use figplot_m
  !use timer_m
  use mesh_extrudate_swell_xfem

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    gintpl = 4,         & ! Q1 gradients
    cintpl = 4,         & ! Q1 conformation
    ointpl = 6,         & ! P2 shape of object elements
    coorsys = 0,        & ! planar Cartesian (0) coordinate system
    physqgrad = 1,      & ! physical quantity nr of the gradients
    physqvel = 2,       & ! physical quantity nr of the velocities
    physqpress = 3,     & ! physical quantity nr of the pressures
    numsplit = 3,       & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit
    numsplitmin = 2,    & ! relative size of all subelements will be at least
                          ! as small as 1/2^numsplitmin
    !numsplit_plot = 2,  & ! relative size of smallest elements near interface
    !                      ! is 1/2^numsplit for plotting
    gausse = 3,         & ! Gauss integration for the eltree subelements
    gausssub_small = 8, & ! Gauss integration for the submesh subelements for
                          ! small integration volumes
    gausssub = 4,       & ! Gauss integration for the submesh subelements
    gauss_ie = 3,       & ! Gauss integration for the interface elements
    gauss = 3,          & ! 3x3 integration of quads
    ncompc = 3,         & ! number of conformation tensor components
    nmodes = 1,         & ! number of modes
    startm = 503,       & ! start of material model data
    model = 2             ! UCM/Oldroyd-B


  real(dp), parameter :: &
    split_threshold = 1.e-12_dp, & ! levelset value for being "close" enough
                                   ! to the interface for tree splitting
    epsvol = 1e-6_dp, & ! elements with integration area smaller are removed
                        ! from the eltree_array. This means that:
                        ! 1) these elements are "fully inside" and no
                        !    "virtual/extended degrees" are generated.
                        ! 2) the interface integration is ignored and therefore
                        !    the interface has a small "hole".
                        ! Note, that epsvol is relative to the area of an
                        ! element, so epsvol=1 is a full element.
    epsvol_small = 1e-2_dp ! elements with integration area smaller
                           ! are integrated with integration rule
                           ! gausssub_small, otherwise with gaussub

! definitions

  !type(mesh_t), target :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc, input_probdefc_proj, &
                           input_probdef_md
  type(problem_t), target :: problem, problemc, problemc_proj, problem_md
  type(sysmatrix_t) :: sysmatrix, sysmatrixc, sysmatrix_md, sysmatrixc_proj
  type(sysvector_t), target :: sol, solm1, sol_md
  type(sysvector_t) :: rhsd, rhsd_md
  type(oldvectors_t) :: oldvectors_ve, oldvectors_md
  type(coefficients_t) :: coefficients, coefficients_md
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc, solcm1, solc_proj
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc, rhsc_proj
  type(lu_ma41_t) :: luc
  type(solver_options_ma41_t) :: solver_options_u, solver_options_c
  type(mesh_t), target :: mesh_n, mesh_nm1
  type(mesh_t), target :: mesh_ALE_np1, mesh_ALE_n


  integer, dimension(:), allocatable :: nodes
  real(dp), dimension(:,:), allocatable :: coor

! array of pointers to an eltree

  type(eltree_p), target, dimension(:,:), allocatable :: eltree_array


! constants surface advection

  integer, parameter :: &
    hintpl = 2,         & ! interpolation for the height, P1 lines
    ninti_sf_adv = 2,   & ! number of Gauss points
    method = 1            ! discretization method 0: Galerkin, 1: SUPG

  real(dp), parameter :: &
    beta_sf_adv = 1.0_dp        ! beta parameter for SUPG


! 1D height function for surface advection

  type(mesh_t), target :: mesh_sf_adv
  type(problem_t), target :: problem_sf_adv
  type(meshgen_options_t) :: mesh_options
  type(input_probdef_t) :: input_probdef_sf_adv
  type(sysmatrix_t) :: sysmatrix_sf_adv
  type(sysvector_t), target :: sol_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1
  type(sysvector_t), target :: sol_sf_adv_pred
  type(sysvector_t) :: rhsd_sf_adv
  type(oldvectors_t) :: oldvectors_sf_adv
  type(coefficients_t) :: coefficients_sf_adv
  type(vector_t), target :: velocity_sf_adv
  !type(plot_options_t) :: plot_options
  type(subscript_t) :: hgt, hgt_end


! variables

  integer :: &
    timeint1 = 1,         & ! (first-order) Euler time integration (first step)
    timeint2 = 7,         & ! (second-order) semi-implicit Gear time integration
    numtimesteps = 5,   & ! number of time steps
    logc = 1,             & ! standard scheme or log transformation
    step0 = 0,            & ! initial step number
    restart = 0             ! restart=1: restart,
                            ! restart=2: restart+Euler first time step

  real(dp) :: &
    eta_s = 0.0_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    gammac = 0.1_dp,   & ! surface tension coefficient
    lambda = 0.2_dp,   & ! relaxation time
    deltat = 1.e-2_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    time0 = 0._dp,     & ! initial time
    U_avg = 1._dp,     & ! average velocity at the entry
    eps_h0 = 1.e-10_dp,& ! initial height of the free boundary points is given
                         ! by h0-eps_h0 to avoid that interface is exactly
                         ! at element edges.
    htype = 2,         & ! upwind parameter
    Uscaling = 3,      & ! upwind parameter
    rs_gup = 1.5_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.5_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.5_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 1.5_dp       ! integer_storage for conformation LU (HSL)


  integer :: icomp, nnsf, nsf
  integer :: step, i, j, nnodes, nbx, nby, obj_ob, obj_sh, els_srf, obj_bp

  real(dp) :: flowrate, H, coor_bp(1,2)=0
  real(dp) :: alpha, G



  logical :: surface_tension = .true.  ! include surface tension


! namelist for input of variables; read from standard input

!  namelist /comppar/ numtimesteps, plot_mesh_every, restart, ctime, &
!    file_append, logc, eta_s, eta_p, lambda, &
!    deltat, htype, Uscaling, &
!    gammac, rs_gup, is_gup, rs_c, is_c, U_avg, eps_h0, surface_tension

!  read ( unit=*, nml=comppar )


! timer

  !timer = .false.

! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus

! usergauss

  intrule_e = gausse  ! Gauss integration for the eltree subelements
  intrule_sub = gausssub  ! Gauss integration for the eltree subelements
  intrule_sub_small = gausssub_small  ! Gauss integration for the eltree
                                      ! subelements (small integration areas)
  vol_small = epsvol_small  ! elements with integration areas smaller
                            ! are integrated with integration rule
                            ! gausssub_small, otherwise with gaussub


! read mesh

  !call read_mesh ( mesh, filename='mesh_esx.out' )
  call mesh_xfem

! add object for open boundary condition

  call add_to_mesh ( mesh, object='curve', objectcurve=4, topology=.true., &
    intrule=3 )
  obj_ob = mesh%nobjects

  if ( surface_tension ) then

!   add object for assembling boundary point

    call add_to_mesh ( mesh, object='coordinates', coor=coor_bp )
    obj_bp = mesh%nobjects

  end if

  nbx = nint(sqrt(real(mesh%curves(12)%nnodes)))
  nby = nint(sqrt(real(mesh%curves(11)%nnodes)))

  call add_to_mesh ( mesh, blocks=[nbx,nby] )

  call fill_mesh_parts ( mesh )

  lmesh => mesh  ! supply mesh to mapcoor for mapping coordinates


! set flowrate

  h = mesh%coor(mesh%points(3),2) ! height is given by y-coordinate of P3
  flowrate = - h * U_avg

! set h0 in functions_xfem

  h0 = mesh%coor(mesh%points(3),2) ! height is given by y-coordinate of P3


! meshes for temporary ALE

  call mesh_convert_copy ( mesh, mesh_ALE_np1 )
  call add_to_mesh ( mesh_ALE_np1, blocks=[nbx,nby] )
  call fill_mesh_parts ( mesh_ALE_np1 )

  call mesh_convert_copy ( mesh, mesh_ALE_n )
  call add_to_mesh ( mesh_ALE_n, blocks=[nbx,nby] )
  call fill_mesh_parts ( mesh_ALE_n )

! For the current problem, mesh_n and mesh_nm1 are not really necessary
! since the computational mesh does not change. (mesh = mesh_n = mesh_nm1)
! However, if local mesh refinements are used,
! the computational mesh can be different at each time step
! Anyway, we still use mesh_n and mesh_nm1 in temporary ALE modules for a
! later general use

  call mesh_convert_copy ( mesh, mesh_n )
  call add_to_mesh ( mesh_n, blocks=[nbx,nby] )
  call fill_mesh_parts ( mesh_n )

  call mesh_convert_copy ( mesh, mesh_nm1 )
  call add_to_mesh ( mesh_nm1, blocks=[nbx,nby] )
  call fill_mesh_parts ( mesh_nm1 )


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=504 )

  coefficients%i = 0
  coefficients%i(1:23) = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gauss,    cintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      logc,     0,          coorsys &
    ]

  coefficients%i(29) = htype
  coefficients%i(31) = Uscaling
  coefficients%i(35) = ointpl
  coefficients%i(37) = 2 ! Gauss points defined by user subroutines
  coefficients%i(41) = gauss_ie ! Gauss points for interface elements
  coefficients%i(49) = 1 ! projection of exp(s) on c
  coefficients%i(50) = 1 ! Courant number dep. tau

  coefficients%r = 0
  coefficients%r(1:10) = &
    [ eta_s,    0._dp,  0._dp, alpha, 0._dp, &
      flowrate, 0._dp, deltat,  beta, 0._dp  &
    ]
  coefficients%r(19) = gammac
  coefficients%r(503:504) = [ G, lambda ]

  coefficients%set_ninti_user => set_ninti_user
  coefficients%set_Gauss_integration_user => set_Gauss_integration_user


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=3 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [  4,0,4,0,4,0,4,0,0,   &  ! G
                   2,2,2,2,2,2,2,2,2,   &  ! velocity
                   1,0,1,0,1,0,1,0,0,   &  ! pressure
                   1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                   [9,4] )

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1


! Dirichlet boundary conditions

! Add nodeset for essential boundary conditions outside of the fluid.
! Temporarily added nodeset with node '1'. Actual nodeset will be replaced
! in find_outside_zero_nodes
  call add_to_mesh ( mesh, nodeset='nodes', nodes=[1] )

! center line
  call define_essential ( mesh, input_probdef, &
    curve1=12, physq=physqvel, degfd=[0,1] )
! outflow
  call define_essential ( mesh, input_probdef, &
    curve1=11, physq=physqvel, degfd=[0,1] )
! wall
  call define_essential ( mesh, input_probdef, &
    curve1=3, physq=physqvel )
! outside of the fluid
  call define_essential ( mesh, input_probdef, nodeset1=1 )


! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=4, nglobalc=1 )

  call problem_definition ( input_probdef, mesh, problem )


! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, solm1, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0
  solm1%u = 0


! fill coefficients for surface advection problem

  call create_coefficients ( coefficients_sf_adv, ncoefi=100, ncoefr=50 )

  coefficients_sf_adv%i = 0
  coefficients_sf_adv%i(1) = ninti_sf_adv
  coefficients_sf_adv%i(2) = 2 ! velocity given by nodal values
  coefficients_sf_adv%i(4) = method ! 0: Galerkin, 1: SUPG
  coefficients_sf_adv%i(6) = hintpl ! height interpolation
  coefficients_sf_adv%i(9) = 3 ! numerical table for Gauss

  coefficients_sf_adv%r = 0
  coefficients_sf_adv%r(4) = deltat
  coefficients_sf_adv%r(5) = beta_sf_adv


! create mesh for surface advection

  mesh_options%elshape = 1   ! two-node line elements

  mesh_options%nx = mesh%curves(7)%nelem * 2

  call line1d ( mesh_sf_adv, mesh_options )

! set coordinates

  do i = 1, mesh%curves(7)%nnodes
    j = mesh%curves(7)%nnodes - i + 1 ! reverse order to get increasing x
    mesh_sf_adv%coor(i,1) = mesh%coor(mesh%curves(7)%nodes(j),1)
  end do

  call fill_mesh_parts ( mesh_sf_adv )


! add object for sampling velocity in height function advection

  allocate ( coor(mesh_sf_adv%nnodes,2) )
  coor(:,1) = mesh_sf_adv%coor(1:mesh_sf_adv%nnodes,1)
  coor(:,2) = h0

  warn_add_to_mesh_after_meshgen_parts = .false.
! real coordinate should be set in surface_corrector
  call add_to_mesh ( mesh, object='coordinates', coor=coor )
  obj_sh = mesh%nobjects
  call fill_mesh_parts_objects ( mesh, object1=obj_sh )
  warn_add_to_mesh_after_meshgen_parts = .true.

  deallocate ( coor )


! create a mesh for plotting the surface

!  call mesh_skeleton ( mesh1, nnodes=mesh_sf_adv%nnodes, &
!    nelem=mesh_sf_adv%nelem, elshape=1, ndim=2 )
!  mesh1%coor = mesh%objects(obj_sh)%coor
!  mesh1%topology = mesh_sf_adv%topology

! merge normal mesh and surface mesh

!  call mesh_merge ( mesh, mesh1, mesh_surface, warn=.false. )
!  call fill_mesh_parts ( mesh_surface )
!  call delete ( mesh1)

!  if ( restart == 0 ) then

!   plot mesh + initial surface

!    write(filename,'(a,i4.4,a)') 'mesh_swell_', 0, '.fig'
!
!    call plot_mesh ( plot_options, mesh_surface, filename=filename, groups=[1] )
!    plot_options%meshcolor = 4
!    call plot_mesh ( plot_options, mesh_surface, filename=filename, &
!      groups=[2], append=.true. )
!    plot_options%meshcolor = 0
!
!  end if


! create input_probdef

  call create_input_probdef ( mesh_sf_adv, input_probdef_sf_adv, nvec=2, &
    nphysq=1 )

  input_probdef_sf_adv%vec_elementdof(1)%a = &
      reshape ( [ 1,1,    &  ! height function
                  2,2 ],  &  ! velocity
                 [2,2] )

  input_probdef_sf_adv%physq = [1]
  input_probdef_sf_adv%probnr = 2

  call define_essential ( mesh_sf_adv, input_probdef_sf_adv, point=1, physq=1 )


! define problem

  call problem_definition ( input_probdef_sf_adv, mesh_sf_adv, problem_sf_adv )

  call create_sysvector ( problem_sf_adv, sol_sf_adv, rhsd_sf_adv )
  call create_sysvector ( problem_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1 )
  call create_sysvector ( problem_sf_adv, sol_sf_adv_pred )


! set pointers in module functions to data

  lmesh_sf_adv => mesh_sf_adv
  lproblem_sf_adv => problem_sf_adv
  lsol_sf_adv_pred => sol_sf_adv_pred


! fill solution vector with essential boundary conditions

  sol_sf_adv%u = h0-eps_h0
  call fill_sysvector ( mesh_sf_adv, problem_sf_adv, sol_sf_adv, point=1, &
    physq=1, value=h0 )

  sol_sf_adv_n%u = sol_sf_adv%u
  sol_sf_adv_pred%u = sol_sf_adv%u


! create system matrix

  call create_sysmatrix_structure ( sysmatrix_sf_adv, mesh_sf_adv, &
    problem_sf_adv )

  call create_sysmatrix_data ( sysmatrix_sf_adv )


! create vectors

  call create_vector ( problem_sf_adv, velocity_sf_adv, vec=2 )

  call create_oldvectors ( oldvectors_sf_adv, nsysvec=2, nvec=1 )

  oldvectors_sf_adv%s(1)%p => sol_sf_adv_n    ! corrector at n
  oldvectors_sf_adv%s(2)%p => sol_sf_adv_nm1  ! corrector at nm1

  oldvectors_sf_adv%v(1)%p => velocity_sf_adv ! advection velocity at np1


! create subscript for the height values

  call create ( mesh_sf_adv, problem_sf_adv, hgt )
  call create ( mesh_sf_adv, problem_sf_adv, hgt_end, points=[2] )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=2, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a =  &
      reshape ( [ 1,0,1,0,1,0,1,0,0,    &  ! c
                  1,1,1,1,1,1,1,1,1 ],  &  ! scalar for plotting
                  [9,2] )

  input_probdefc%physq = [1]
  input_probdefc%probnr = 3

  call define_essential ( mesh, input_probdefc, nodeset1=1 )

  call problem_definition ( input_probdefc, mesh, problemc )

! create system vectors (solution and right-hand side) for conformation and

  call create ( problemc, solc, solcm1, rhsc )

! initialize vectors with zero stress

  if ( logc == 0 ) then ! standard
    solc(1,1)%u = 1   ! initial cxx
    solc(2,1)%u = 0   ! initial cxy
    solc(3,1)%u = 1   ! initial cyy
  else if ( logc == 1 ) then ! log scheme
    solc(1,1)%u = 0   ! initial sxx
    solc(2,1)%u = 0   ! initial sxy
    solc(3,1)%u = 0   ! initial syy
  end if

  call copy ( solc, solcm1 )


! problem definition for projected "c=exp(s)" in case of log conformation

  call create_input_probdef ( mesh, input_probdefc_proj, nvec=1, nphysq=1 )

  input_probdefc_proj%vec_elementdof(1)%a = &
      reshape ( [ 1,0,1,0,1,0,1,0,0 ], &
                   [9,1] )

  input_probdefc_proj%physq = [1]
  input_probdefc_proj%probnr = 4

  call define_essential ( mesh, input_probdefc_proj, nodeset1=1 )

  call problem_definition ( input_probdefc_proj, mesh, problemc_proj )

  call create ( problemc_proj, solc_proj, rhsc_proj )

  solc_proj(1,1)%u = 1   ! initial cxx
  solc_proj(2,1)%u = 0   ! initial cxy
  solc_proj(3,1)%u = 1   ! initial cyy


! problem definition for mesh displacement used for temporary ALE mesh
! generation

  call create_coefficients ( coefficients_md, ncoefi=100, ncoefr=50 )

  coefficients_md%i = 0
  coefficients_md%i(1) = uintpl
  coefficients_md%i(10) = gauss
  coefficients_md%i(17) = -1 ! normal to fluid points inwards to the object
  coefficients_md%i(41) = gauss_ie ! Gauss points for interface elements
  coefficients_md%i(42) = 1   ! -1: impose zero, 0: impose non-zero value
                              !  1: impose func_md
  coefficients_md%i(43) = 1   ! funcnr, but not used in func_md
  coefficients_md%func1(4)%p => func_md

  coefficients_md%r = 0
  coefficients_md%r(1) = 1._dp
  coefficients_md%r(19) = 1._dp   ! embedded factor


  call create_input_probdef ( mesh, input_probdef_md, nvec=1, nphysq=1 )

  input_probdef_md%vec_elementdof(1)%a = 1

  input_probdef_md%physq = [1]
  input_probdef_md%probnr = 5

  ! wall
  call define_essential ( mesh, input_probdef_md, curve1=3, physq=1 )
  ! inflow
  call define_essential ( mesh, input_probdef_md, curve1=4, physq=1 )
  ! center line
  call define_essential ( mesh, input_probdef_md, curve1=12, physq=1 )
  ! top of swell zone
  call define_essential ( mesh, input_probdef_md, curve1=9, physq=1 )


  call problem_definition ( input_probdef_md, mesh, problem_md )

  call create_sysvector ( problem_md, sol_md, rhsd_md )

  sol_md%u = 0

  call create_sysmatrix_structure ( sysmatrix_md, mesh, problem_md )
  call create_sysmatrix_data ( sysmatrix_md )


! allocate array of pointers to an eltree (eltree for all elements, one group)

  allocate ( eltree_array(1,mesh%grpnumel(1)), vol(mesh%grpnumel(1)) )

  allocate ( nodes(mesh%nnodes), d(mesh%nnodes) )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec2=3, nprob=3, &
    nelta=1, nmesh=4 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s(2)%p => solm1
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%s2(2)%p => solcm1
  oldvectors_ve%s2(3)%p => solc_proj
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc
  oldvectors_ve%p(3)%p => problemc_proj
  oldvectors_ve%ea(1)%p => eltree_array
  oldvectors_ve%m(1)%p => mesh_ALE_np1
  oldvectors_ve%m(2)%p => mesh_n
  oldvectors_ve%m(3)%p => mesh_ALE_n
  oldvectors_ve%m(4)%p => mesh_nm1
  ea => eltree_array ! userelmesh does not have oldvectors as an argument


! oldvectors for mesh displacement

  call create_oldvectors ( oldvectors_md, nsysvec=1, nelta=1 )
  oldvectors_md%s(1)%p => sol_md
  oldvectors_md%ea(1)%p => eltree_array



! time stepping

  coefficients%i(22) = timeint1  ! first step integration scheme
  coefficients_sf_adv%i(5) = 1 ! start with first-order scheme

  do step = 1, numtimesteps

    if ( step >= 2 .or. restart == 1 ) then

      coefficients%i(22) = timeint2
      coefficients_sf_adv%i(5) = 2  ! second-order scheme

!     predict position of the surface
      sol_sf_adv_pred%u = 2._dp*sol_sf_adv_n%u - sol_sf_adv_nm1%u

      call get_temporary_ale_mesh

      call delete_eltree_in_elements


    end if


!   eltree for elements crossing the interface

    d = levelset ( mesh%coor ) ! set levelset for all nodes

    call create_eltree_in_elements ( numsplit )

!   nodeset for Dirichlet on outside nodes

    call find_outside_zero_nodes

!   elementset for elements crossed by the surface

    if ( step == 1 ) then
      call create_elementset_from_eltree
      els_srf = mesh%nelementsets
    else
      call create_elementset_from_eltree ( replace=els_srf )
    end if


!   re-define problems

    call redefine_problems


!   create system matrix of gradient/velocity/pressure problem

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )
    call create_sysmatrix_data ( sysmatrix )

    call create_sysmatrix_structure ( sysmatrixc, mesh, problemc )
    call create_sysmatrix_data ( sysmatrixc )


!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG


!   build implicit terms of CE

    if ( step >= 2 .or. restart == 1 ) then
      coefficients%i(48) = 2 ! tALE
    else
      coefficients%i(48) = 0
    end if

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve, physqrow=[physqvel], physqcol=[physqvel], &
      addmatvec=.true., coefficients=coefficients )

    coefficients%i(48) = 0

!   open boundary

    call build_system(mesh, problem, sysmatrix, rhsd, &
      elemsub1=stokes_open_boundary, oldvectors=oldvectors_ve, &
      physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress], &
      coefficients=coefficients, object=obj_ob, addmatvec=.true. )

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub1=implicit_stress_open_boundary, physqrow=[physqvel], &
      physqcol=[physqvel], object=obj_ob, addmatvec=.true., &
      oldvectors=oldvectors_ve, coefficients=coefficients )

!   imposed flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients, oldvectors=oldvectors_ve )


!   surface tension

    if ( surface_tension ) then

!     build elements in the surface

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=surface_tension_eltree, &
        physqrow=[physqvel], physqcol=[physqvel], &
        addmatvec=.true., buildmatrix=.false., coefficients=coefficients, &
        oldvectors=oldvectors_ve, elementset=els_srf )

!     boundary point x1

      nnsf = mesh_sf_adv%nnodes
      nsf = sol_sf_adv_pred%n
      mesh%objects(obj_bp)%coor(1,1) = mesh_sf_adv%coor(nnsf,1)
      mesh%objects(obj_bp)%coor(1,2) = sol_sf_adv_pred%u(hgt%s(nsf))
      call find_refcoor_objects ( mesh, object1=obj_bp )

      if ( mesh%objects(obj_bp)%grpelm(1,1) == 0 ) then
        print*, 'Error: no reference coordinates in boundary point.'
        stop
      end if

!     second point x2

      coefficients%r(501) = mesh_sf_adv%coor(nnsf-1,1)
      coefficients%r(502) = sol_sf_adv_pred%u(hgt%s(nsf-1))

!     build external force in boundary point

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub1=surface_tension_boundary_point2, &
        object=obj_bp, onobjectnodes=.true., &
        physqrow=[physqvel], physqcol=[physqvel], &
        addmatvec=.true., buildmatrix=.false., coefficients=coefficients, &
        oldvectors=oldvectors_ve )


    end if

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

    call copy ( sol, solm1 )

    call delete ( sysmatrix )



!   build (assemble) matrix and vector for conformation problem

    if ( coefficients%i(22) == timeint1 ) then

      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    else

      coefficients%i(48) = 2 ! tALE

      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem_implicit_2nd_order, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

      coefficients%i(48) = 0 ! unset tALE

    end if

    call check ( sysmatrixc )

    call copy ( solc, solcm1 )



!   solve conformation and keep LU decomposition in the loop over components

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do icomp = 1, ncompc
      call solve_system_ma41 ( sysmatrixc, rhsc(icomp,1), solc(icomp,1), luc, &
        solver_options=solver_options_c  )
    end do

    call delete ( luc )  ! remove LU decomposition and rebuild next time step

    call delete ( sysmatrixc )



!   exps projection on c if log-conformation

    if ( logc == 1 ) call solve_exps_projection


!   solve surface advection (corrector)

    call solve_surface_height_corrector



!    if ( mod(step0+step,plot_mesh_every) == 0 ) then
!
!!     plot mesh + surface
!
!      write(filename,'(a,i4.4,a)') 'mesh_swell_', step0+step, '.fig'
!
!      n2 = mesh_surface%nnodes
!      n1 = n2 - mesh%objects(obj_sh)%nnodes + 1
!      mesh_surface%coor(n1:n2,:) = mesh%objects(obj_sh)%coor
!      call plot_mesh ( plot_options, mesh_surface, filename=filename, &
!        groups=[1] )
!      plot_options%meshcolor = 4  ! red
!      call plot_mesh ( plot_options, mesh_surface, filename=filename, &
!        groups=[2], append=.true. )
!      plot_options%meshcolor = 0  ! reset to black
!
!    end if
!
  end do

  close(unit=11)


! delete all data including all allocated memory

  call delete_eltree_in_elements

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, solm1, rhsd )
  call delete ( oldvectors_ve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( solc, solcm1, rhsc )
  call delete ( coefficients )

  call delete ( input_probdef_md )
  call delete ( problem_md )
  call delete ( oldvectors_md )
  call delete ( sol_md, rhsd_md )
  call delete ( coefficients_md )
  call delete ( sysmatrix_md )

  call delete ( input_probdefc_proj )
  call delete ( problemc_proj )
  call delete ( solc_proj, rhsc_proj )

  call delete ( mesh_n, mesh_nm1, mesh_ALE_np1, mesh_ALE_n )

  call delete ( mesh_sf_adv )
  call delete ( problem_sf_adv )
  call delete ( input_probdef_sf_adv )
  call delete ( sol_sf_adv, rhsd_sf_adv )
  call delete ( sol_sf_adv_pred, sol_sf_adv_n, sol_sf_adv_nm1 )
  call delete ( sysmatrix_sf_adv )
  call delete ( oldvectors_sf_adv )
  call delete ( coefficients_sf_adv )
  call delete ( hgt, hgt_end )

  deallocate ( nodes, d, eltree_array, vol )


contains


! build (assemble) matrix and vector for gradient/velocity/pressure problem

  subroutine build_vpG

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, physqrow=[physqvel,physqpress], &
      physqcol=[physqvel,physqpress], coefficients=coefficients, &
      oldvectors=oldvectors_ve )

!   DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_elem, addmatvec=.true., &
      physqrow=[physqgrad,physqvel], physqcol=[physqgrad,physqvel], &
      coefficients=coefficients, oldvectors=oldvectors_ve )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqgrad], physqcol=[physqpress], &
      zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqpress], physqcol=[physqgrad], &
      zeromatvec=.true. )

  end subroutine build_vpG


  subroutine redefine_problems

    integer :: m, i

!   system numbering -> nodal numbering
    sol%u = sol%u( problem%degfdperm(:,2) )
    solm1%u = solm1%u( problem%degfdperm(:,2) )

    do m = 1, nmodes
      do i = 1, ncompc
        solc(i,m)%u = solc(i,m)%u( problemc%degfdperm(:,2) )
        solcm1(i,m)%u = solcm1(i,m)%u( problemc%degfdperm(:,2) )
      end do
    end do

    call delete ( problem )
    call problem_definition ( input_probdef, mesh, problem )

    call delete ( problemc )
    call problem_definition ( input_probdefc, mesh, problemc )

!   nodal numbering -> system numbering
    sol%u = sol%u( problem%degfdperm(:,1) )
    solm1%u = solm1%u( problem%degfdperm(:,1) )

    do m = 1, nmodes
      do i = 1, ncompc
        solc(i,m)%u = solc(i,m)%u( problemc%degfdperm(:,1) )
        solcm1(i,m)%u = solcm1(i,m)%u( problemc%degfdperm(:,1) )
      end do
    end do

    if ( logc == 1 ) then

!     system numbering -> nodal numbering
      do m = 1, nmodes
        do i = 1, ncompc
          solc_proj(i,m)%u = solc_proj(i,m)%u( problemc_proj%degfdperm(:,2) )
        end do
      end do

      call delete ( problemc_proj )
      call problem_definition ( input_probdefc_proj, mesh, problemc_proj )

      call create_sysmatrix_structure ( sysmatrixc_proj, mesh, problemc_proj, &
        symmetric=.true. )
      call create_sysmatrix_data ( sysmatrixc_proj )

!     nodal numbering -> system numbering
      do m = 1, nmodes
        do i = 1, ncompc
          solc_proj(i,m)%u = solc_proj(i,m)%u( problemc_proj%degfdperm(:,1) )
        end do
      end do

    end if

  end subroutine redefine_problems


! build and solve for exps projection

  subroutine solve_exps_projection

    use hsl_ma57_m

    type(solver_options_ma57_t) :: solver_options_ma57
    type(lu_ma57_t) :: lu_exps_proj

    integer :: i, m

    ! MA57 solver storage
    solver_options_ma57%integer_storage = 1.3
    solver_options_ma57%real_storage    = 1.3


    call build_system ( mesh, problemc_proj, sysmatrixc_proj, &
                        m2sysvector=rhsc_proj, elemsub=exps_projection_elem, &
                        oldvectors=oldvectors_ve, coefficients=coefficients )

    call check ( sysmatrixc_proj )




   do m = 1, nmodes
      do i = 1, ncompc

        call add_effect_of_essential_to_rhs ( problemc_proj, sysmatrixc_proj, &
           solc_proj(i,m), rhsc_proj(i,m) )

        call solve_system_ma57 ( sysmatrixc_proj, rhsc_proj(i,m), &
           solc_proj(i,m), lu_exps_proj, solver_options=solver_options_ma57 )

      end do
    end do

    call delete ( lu_exps_proj )

    call delete ( sysmatrixc_proj )


  end subroutine solve_exps_projection


! compute temporary ALE mesh

 subroutine get_temporary_ale_mesh

    use convection_diffusion_elements_m

    implicit none

    type(solver_options_ma41_t) :: solver_options_ma41

!   fictitious domain approach for mesh displacement

    call build_system ( mesh, problem_md, sysmatrix_md, rhsd_md, &
      elemsub=scalar_diffusion_elem, coefficients=coefficients_md, &
      oldvectors=oldvectors_md )

    call build_system ( mesh, problem_md, sysmatrix_md, rhsd_md, &
      elemsub=scalar_diffusion_embedded_dirichlet_eltree, &
      oldvectors=oldvectors_md, coefficients=coefficients_md, &
      addmatvec=.true., elementset=els_srf )

    call add_effect_of_essential_to_rhs ( problem_md, sysmatrix_md, sol_md, &
      rhsd_md )

    solver_options_ma41%integer_storage = 1.5
    solver_options_ma41%real_storage    = 1.5

    call solve_system_ma41 ( sysmatrix_md, rhsd_md, sol_md, &
      solver_options=solver_options_ma41 )


!   temporary ALE mesh

    mesh_ALE_n%coor = mesh_ALE_np1%coor

    mesh_ALE_np1%coor(:,2) = mesh%coor(:,2) &
                                + sol_md%u( problem_md%degfdperm(:,2) )

    call find_bounds_blocks ( mesh_ALE_np1 )
    call find_bounds_blocks ( mesh_ALE_n )

  end subroutine get_temporary_ale_mesh



! solve convection equation for the surface height (corrector)

  subroutine solve_surface_height_corrector

    use postprocessing_m

    type(sample_t) :: sample
    real(dp) :: max_height, end_height
    type(solver_options_ma41_t) :: solver_options_h

    nnodes = mesh_sf_adv%nnodes

    mesh%objects(obj_sh)%coor(1:nnodes,1) = mesh_sf_adv%coor(1:nnodes,1)
    mesh%objects(obj_sh)%coor(1:nnodes,2) = sol_sf_adv_pred%u(hgt%s)

    call find_refcoor_objects ( mesh, object1=obj_sh )

    if ( any (mesh%objects(obj_sh)%grpelm(:,1) == 0) ) then
      print*, 'Error: no reference coordinates in surface advection corrector.'
      print*, 'Program Stop!'
      stop
    end if

    call fill_sample ( mesh, problem, sample, ndegfd=2, object=obj_sh, &
      elemsub=stokes_sample_velocity, coefficients=coefficients, &
      oldvectors=oldvectors_ve )

    velocity_sf_adv%u = reshape ( transpose(sample%u), [2*nnodes] )


    call build_system ( mesh_sf_adv, problem_sf_adv, sysmatrix_sf_adv, &
      rhsd_sf_adv, elemsub=surface_advection_elem, &
      oldvectors=oldvectors_sf_adv, coefficients=coefficients_sf_adv )


    call check_filled_sysmatrix ( sysmatrix_sf_adv )

    call add_effect_of_essential_to_rhs ( problem_sf_adv, sysmatrix_sf_adv, &
      sol_sf_adv, rhsd_sf_adv )


!   MA41 solver storage

    solver_options_h%integer_storage = 2.0
    solver_options_h%real_storage    = 2.0

!   solve system

    call solve_system_ma41 ( sysmatrix_sf_adv, rhsd_sf_adv, sol_sf_adv, &
                             solver_options=solver_options_h )


    call copy ( sol_sf_adv_n, sol_sf_adv_nm1 )
    call copy ( sol_sf_adv, sol_sf_adv_n )

!   compute the maximum, minimum and end radius

    max_height = maxval ( sol_sf_adv%u )
    end_height = sol_sf_adv%u(hgt_end%s(1))


!   write swell height

    !write(11,'(i6,4es16.8)') step0+step, time0+step*deltat, &
    !                         max_height, end_height
    print '(i6,4es16.8)', step0+step, time0+step*deltat, &
                          max_height, end_height

  end subroutine solve_surface_height_corrector


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

    !print *, 'minvalvol = ', minvalvol

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


! find outside nodes that are not connected to elements at the interface

  subroutine find_outside_zero_nodes

    integer :: elem, i

!   set outside nodes

    where ( d < 0._dp )  ! use < 0 (and not <= 0) to avoid wall nodes to
      nodes = 1          ! end up in the nodeset
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

    call add_to_mesh ( mesh, nodeset='nodes', nodes=nodes(1:nnodes), replace=1 )

  end subroutine find_outside_zero_nodes

end program extrudate_swell5
