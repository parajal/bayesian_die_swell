! 2D viscoelastic problem on a square domain with a disk at the center.
! Shear flow boundary conditions (set in "stokes_functions_m").
! Freely rotating conditions on the disk boundaries.
! Gmsh mesh (in the file "meshP2_2.msh").
! Integration of fields over the domain and boundaries.
! Conformation tensor for the unfilled fluid at inflow sections.
!
! DEVSS-G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! first-order time integration to steady state

module subs_m

use tfem_elem_m

  implicit none

contains

  subroutine elementc ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: ncurveconstr
    real(dp) :: xr(1,2)

    ncurveconstr = problem%constraints(constr)%geometry1

    xr(1,:) = mesh%coor(mesh%curves(ncurveconstr)%nodes(node),:)

!   set shape function in the point

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

!   connection through collocation

    elemmat(1,1) = 1
    elemmat(1,2) = 0

    elemmat(2,1) = 0
    elemmat(2,2) = 1

    elemmatadd(1,:) = [ -xr(1,2) ]
    elemmatadd(2,:) = [ xr(1,1) ]

  end subroutine elementc

end module subs_m

program bulkstress2

  use tfem_m
  use hsl_ma41_m
  use viscoelastic_elements_m
  use figplot_m
  use io_utils_m
  use subs_m
  use visco_pure_m
  use stokes_functions_m

  implicit none

  integer, parameter :: &
    uintpl = 6,     &  ! P2 velocities
    pintpl = 2,     &  ! P1 pressures
    gintpl = 2,     &  ! P1 gradients
    cintpl = 2,     &  ! P1 conformation
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    gauss = 6,      & ! 6-point Gauss integration of triangles
    gaussb = 3,     & ! 3 point integration of boundary elements
    ncompc = 3,     & ! number of conformation tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    coorsys = 0,    & ! coordinate system
    model = 3         ! Giesekus model

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc
  type(problem_t), target :: problem, problemc
  type(sysmatrix_t) :: sysmatrix, sysmatrixc
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve, oldvectors_d
  type(coefficients_t) :: coefficients
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc
  type(lu_ma41_t) :: luc
  type(plot_options_t) :: plot_options
  type(vector_t), target :: pressure, conformation_tensor(nmodes)
  type(solver_options_ma41_t) :: solver_options_u, solver_options_c

! variables

  integer :: &
    timeint = 1,         & ! (first-order) time integration
    numtimesteps = 10,   & ! number of time steps
    logc = 1               ! standard scheme or log transformation

  real(dp) :: &
    eta_s = 0.0_dp,     & ! solvent viscosity
    eta_p = 1.0_dp,     & ! polymer viscosity
    lambda = 0.4_dp,    & ! relaxation time
    mobility = 0.5_dp,  & ! mobility parameter in the Giesekus model (alpha)
    shearrate = 1.0_dp, & ! shear rate for calculation of pure fluid
    deltat = 1.e-2_dp,  & ! time step
    beta = 1.0_dp,      & ! upwinding parameter in the SUPG method
    rs_gup = 1.4_dp,    & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.4_dp,    & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.7_dp,     & ! real_storage for the conformation LU (HSL)
    is_c  = 1.8_dp        ! integer_storage for conformation LU (HSL)

  integer :: icomp, step, i, ip
  real(dp) :: alpha, G, conf_table(3)
  real(dp) :: tau_tensor(3), area_pressure(1), area_domain(1), boundary_tau(4), &
    boundary_pressure(4), disk_length(1)

! set some parameters
  alpha = eta_p         ! DEVSS parameter
  G = eta_p / lambda    ! modulus

! fill coefficients
  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+3*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
       physqvel, physqpress, 0,     physqgrad, gauss,  &
       gaussb,   cintpl,     0,     0,         0,      &
       0,        0,          model, nmodes,    startm, &
       logc,     timeint,    coorsys, ( 0, i = 24, 150 )  &
    ]

  coefficients%r = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
       0._dp, 0._dp,  deltat,   beta,  0._dp, &
       ( 0._dp, i = 11, 500 ), &
       G,     lambda,     mobility   ]

! calculation of conformation tensor for a pure fluid
  call visco_pure_calc (model, nmodes, ncompc, numtimesteps, &
    G, lambda, mobility, shearrate, deltat, coorsys, logc)

  open( unit=87, file='boundary_c' )

! read mesh
  call read_mesh_gmsh ( mesh, filename='meshP2_2.msh', ndim=2, physgeom=.true. )

  call fill_mesh_parts ( mesh )

! plot mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  plot_options%fontsize = 10

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition of gradient/velocity/pressure
  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=3 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [ 4,0,4,0,4,0,    &  ! G
                  2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,    &  ! pressure
                  1,1,1,1,1,1,    &  ! scalar, such as vorticity
                  3,3,3,3,3,3 ], &  ! tensor
                  [6,5] )

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

! define essential boundaries

! cell boundaries
  call define_essential ( mesh, input_probdef, curve1=1, curve2=6, physq=2 )

! disk
  call define_constraint ( mesh, input_probdef, curve1=7, physq=2, &
    discretization='collocation', nodedof=2, naddunknowns=1 )

! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=3 )

  call problem_definition ( input_probdef, mesh, problem )

! problem definition conformation tensor
  call create_input_probdef ( mesh, input_probdefc, nvec=3, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a =   &
      reshape ( [ 1,0,1,0,1,0,    &  ! c
                  1,1,1,1,1,1,    &  ! scalar for plotting
                  3,3,3,3,3,3 ], &  ! tensor
                  [6,3] )

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

! essential bc for conformation at the inflow boundary
  call define_essential ( mesh, input_probdefc, curve1=2 )
  call define_essential ( mesh, input_probdefc, curve1=5 )

  call problem_definition ( input_probdefc, mesh, problemc )

! create system vectors (solution and right-hand side)
  call create_sysvector ( problem, sol, rhsd )

! fill solution vector with essential boundary conditions

! cell boundaries
  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=6, physq=2, degfd=1, func=func, funcnr=1 )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=6, physq=2, degfd=2, func=func, funcnr=2 )

! pressure level
  call fill_sysvector ( mesh, problem, sol, point=1, physq=3, value=0._dp )

! create system matrix of gradient/velocity/pressure problem
  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! create system vectors (solution and right-hand side) for conformation and
! initialize vectors with zero stress
  call create ( problemc, solc, rhsc )

  if ( logc == 0 ) then ! standard
    solc(1,1)%u = 1   ! initial cxx
    solc(2,1)%u = 0   ! initial cxy
    solc(3,1)%u = 1   ! initial cyy
  else if ( logc == 1 ) then ! log scheme
    solc(1,1)%u = 0   ! initial sxx
    solc(2,1)%u = 0   ! initial sxy
    solc(3,1)%u = 0   ! initial syy
  end if

! create system matrix for conformation problem
  call create_sysmatrix_structure ( sysmatrixc, mesh, problemc )
  call create_sysmatrix_data ( sysmatrixc )

! create the structure oldvectors_ve
  call create_oldvectors ( oldvectors_ve, nsysvec=1, nsysvec2=1, nvec1=1, &
    nprob=2 )
  call create_oldvectors ( oldvectors_d, nsysvec=1, nvec=2 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc

  oldvectors_d%s(1)%p => sol

! time stepping

  do step = 1, numtimesteps

!   build (assemble) matrix/vector for gradient/velocity/pressure problem
    call build_vpG

!   build implicit terms
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2],&
      addmatvec=.true., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve gradient/velocity/pressure problem
    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

!   build (assemble) matrix and vector for conformation problem
    call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
      elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
      coefficients=coefficients )

    read(87,'(1X,3F14.6)') conf_table

!   fill solution vector with essential boundary conditions for
!   conformation tensor
    do icomp = 1, ncompc
      call fill_sysvector ( mesh, problemc, solc(icomp,1), &
        curve1=2, value=conf_table(icomp) )
      call fill_sysvector ( mesh, problemc, solc(icomp,1), &
        curve1=5, value=conf_table(icomp) )
    end do

    call check ( sysmatrixc )

    do icomp = 1, ncompc
      call add_effect_of_essential_to_rhs ( problemc, sysmatrixc, &
        solc(icomp,1), rhsc(icomp,1) )
    end do

!   solve conformation and keep LU decomposition in the loop over components
    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do icomp = 1, ncompc
      call solve_system_ma41 ( sysmatrixc, rhsc(icomp,1), solc(icomp,1), luc, &
        solver_options=solver_options_c  )
    end do

    call delete ( luc )  ! remove LU decomposition and rebuild next time step

!*****************************************************************************
!************     Post-processing     ****************************************
!**************************************************************************** *

    print *,'**************************************************************'
    print *,'Step: ', step
    print *,'Time: ', step*deltat

    ip = problem%constraints(1)%addnumdegfd(1)
    print *,'Rotation rate of disk: ', sol%u( problem%degfdperm(ip+1,2) )

!*******************
!**  Integration  **
!*******************

!   Integration of quantities over the domain

!   Get pressure
    call create_vector ( problem, pressure, vec=4 )
    call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

!   Area of the domain
    call integrate ( mesh, problem, area_domain, &
      elemsub=stokes_integrate_volume, &
      coefficients=coefficients, oldvectors=oldvectors_ve)

    print *,'Area of the domain: ', area_domain(1)

!   Pressure
    call integrate ( mesh, problem, area_pressure, &
      elemsub=stokes_integrate_pressure, &
      coefficients=coefficients, oldvectors=oldvectors_ve)

    print *,'Integral of pressure over the domain/area: ', &
      area_pressure(1)/400.0_dp

!   Polymer stress tensor
    call integrate ( mesh, problemc, tau_tensor, &
      elemsub=viscoelastic_integrate_tau_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve)

    print *,'xx-component of polymer tensor/area: ', tau_tensor(1)/400.0_dp
    print *,'xy-component of polymer tensor/area: ', tau_tensor(2)/400.0_dp
    print *,'yx-component of polymer tensor/area: ', tau_tensor(2)/400.0_dp
    print *,'yy-component of polymer tensor/area: ', tau_tensor(3)/400.0_dp

!   Get conformation tensor
    call create ( problemc, conformation_tensor, vec=3 )
    call derive_vector ( mesh, problemc, conformation_tensor(1), &
      elemsub=deriv_conformation_tensor, coefficients=coefficients, &
      oldvectors=oldvectors_ve )

    oldvectors_ve%v1(1)%p => conformation_tensor
    oldvectors_d%v(2)%p => pressure

!   Integration of quantities over the disk boundary

!   Length of disk boundary
    call integrate_boundary_elements ( mesh, problem, disk_length, &
      elemsub=stokes_integrate_area_geometry, curve=7, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    print *,'Length of the disk boundary: ', disk_length(1)

!   Particle pressure contribution
    call integrate_boundary_elements ( mesh, problem, boundary_pressure, &
      elemsub=stokes_integrate_pressure_geometry, curve=7, &
      coefficients=coefficients, oldvectors=oldvectors_d )

    print *,'xx component of (p*I.n).x over the disk boundary/area: ', &
      boundary_pressure(1)/400.0_dp
    print *,'xy component of (p*I.n).x over the disk boundary/area: ', &
      boundary_pressure(2)/400.0_dp
    print *,'yx component of (p*I.n).x over the disk boundary/area: ', &
      boundary_pressure(3)/400.0_dp
    print *,'yy component of (p*I.n).x over the disk boundary/area: ', &
      boundary_pressure(4)/400.0_dp

!   Particle polymer stress contribution
    call integrate_boundary_elements ( mesh, problemc, boundary_tau, &
      elemsub=viscoelastic_integrate_tau_geometry, curve=7, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    print *,'xx component of (tau.n).x over the disk boundary/area: ', &
      boundary_tau(1)/400.0_dp
    print *,'xy component of (tau.n).x over the disk boundary/area: ', &
      boundary_tau(2)/400.0_dp
    print *,'yx component of (tau.n).x over the disk boundary/area: ', &
      boundary_tau(3)/400.0_dp
    print *,'yy component of (tau.n).x over the disk boundary/area: ', &
      boundary_tau(4)/400.0_dp

    print *,'**************************************************************'

    call delete ( pressure )
    call delete ( conformation_tensor )

  end do

! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( solc, rhsc )
  call delete ( coefficients )
  call delete ( oldvectors_ve, oldvectors_d )

contains

  subroutine build_vpG

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, physqrow=[2,3], physqcol=[2,3], &
      coefficients=coefficients )

!   DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_elem, addmatvec=.true., &
      physqrow=[1,2], physqcol=[1,2], coefficients=coefficients )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[1], physqcol=[3], zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[3], physqcol=[1], zeromatvec=.true. )

!   constraints
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      elemsub=elementc, addmatvec=.true. )

  end subroutine build_vpG

end program bulkstress2
