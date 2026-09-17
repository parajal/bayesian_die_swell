! Elastoviscoplastic problem on a unit square.
! Saramito model (2021) with Drucker-Prager plasticity.
! Start-up of plane Poiseuille flow having a constant flow rate.
! Upper half of the domain using symmetry boundary conditions at the bottom.
! DEVSS-G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! L2-projection of c=exp(s) on the finite element approximation space in the
! momentum balance for logc=1, but is recommended to use logc=0 for this
! problem.
! second-order time integration

program poiseuille4

  use tfem_m
  use viscoelastic_elements_m
  use inertia_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities
    pintpl = 4,     & ! Q1 pressures
    gintpl = 4,     & ! Q1 gradients
    cintpl = 4,     & ! Q1 conformation
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    gauss = 3,      & ! 3x3 integration of quads
    ncompc = 4,     & ! number of conformation tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    coorsys = 0       ! planar Cartesian coordinate system


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc, input_probdefc_proj
  type(problem_t), target :: problem, problemc, problemc_proj
  type(sysmatrix_t) :: sysmatrix, sysmatrixc, sysmatrixc_proj
  type(sysvector_t), target :: sol, soln, solnm1, sol_hat
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve, oldvectors
  type(coefficients_t) :: coefficients
  type(subscript_t) :: velx, cval, vel, press
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc, solcm1, solc_proj
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc, rhsc_proj
  type(lu_ma41_t) :: luc
  type(solver_options_ma41_t) :: solver_options_u, solver_options_c
  type(lu_ma57_t) :: lu_exps_proj


! variables

  integer :: &
    nx=20,               & ! number of elements in x
    ny=20,               & ! number of elements in y
    timeint1 = 1,        & ! (first-order) time integration (first step)
    timeint2 = 7,        & ! (second-order) time integration
    numtimesteps = 600,  & ! number of time steps
    inertia = 0,         & ! include instationary inertia term
    logc = 1,            & ! standard scheme or log transformation
    model = 22             ! Saramito model with Drucker-Prager plasticity

  real(dp) :: &
    rho = 1.0_dp,      & ! density
    eta_s = 0.1_dp,    & ! solvent viscosity
    deltat = 5.e-3_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    U_avg = 1._dp,     & ! average velocity at the entry
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.0_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 1.4_dp       ! integer_storage for conformation LU (HSL)

  real(dp), dimension(nmodes) :: &
    G,         & ! modulus
    lambda,    & ! relaxation time
    tau_y,     & ! yield stress (cohesion)
    mu           ! friction coefficient

  integer :: icomp, step, i, j, m
  real(dp) :: alpha, mfac
  real(dp) :: flowrate, H
  real(dp) :: gamma0=1.5_dp, alpha0=2._dp, alpha1=-0.5_dp


! namelist for input of variables; read from standard input

  namelist /comppar/ nx, ny, timeint1, timeint2, numtimesteps, logc, inertia, &
    rho, eta_s, G, lambda, tau_y, mu, deltat, &
    rs_gup, is_gup, rs_c, is_c, U_avg

  read ( unit=*, nml=comppar )


! set some parameters

  alpha = sum(G*lambda)  ! DEVSS parameter


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=250, ncoefr=500*4*nmodes )

  coefficients%i = 0
  coefficients%i(1:23) = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gauss,    cintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      logc,     timeint1,   coorsys                   &
    ]

  coefficients%i(49) = 1  ! exp(s) projection =.true. for logc=1

  coefficients%r = 0
  coefficients%r(1:10) = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
      0._dp, 0._dp,  deltat,   beta,  0._dp  &
    ]
  coefficients%r(151) = rho

  coefficients%r(501:500+4*nmodes) = &
         [ ( G(i), lambda(i), tau_y(i), mu(i), i=1,nmodes ) ]


! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )

! set flowrate

  H = mesh%coor(mesh%points(4),2) ! height is given by y-coordinate of P4
  flowrate = - H * U_avg

  coefficients%r(6) = flowrate


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=3 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [ 4,0,4,0,4,0,4,0,0,    &  ! G
                  2,2,2,2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,1,0,0,    &  ! pressure
                  1,1,1,1,1,1,1,1,1,    &  ! scalar, such as vorticity
                  3,3,3,3,3,3,3,3,3 ], &  ! tensor
                  [9,5] )

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

! Dirichlet boundary conditions

! center line
  call define_essential ( mesh, input_probdef, &
    curve1=1, physq=physqvel, degfd=[0,1] )

! wall
  call define_essential ( mesh, input_probdef, &
    curve1=3, physq=physqvel )

! pressure in point 1
  call define_essential ( mesh, input_probdef, point=1, physq=3 )

! constraints for periodical boundary conditions

  call define_constraint ( mesh, input_probdef, &
    physq=2, curve1=2, curve2=5, discretization='weak', elementdof=[2,0,2] )

! gradients
  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, curve2=5, discretization='collocation' )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=4, nglobalc=1 )

  call problem_definition ( input_probdef, mesh, problem )

! create a vector subscripts for post processing
! horizontal velocity
  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
! velocities
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )
! pressure
  call create_subscript ( mesh, problem, press, physqarr=[physqpress], &
    points=[1,2] )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=3, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a =   &
      reshape ( [ 1,0,1,0,1,0,1,0,0,    &  ! c
                  1,1,1,1,1,1,1,1,1,    &  ! scalar for plotting
                  4,4,4,4,4,4,4,4,4 ], &  ! tensor
                  [9,3] )

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

! constraint for periodical boundary conditions of the conformation
  call define_constraint ( mesh, input_probdefc, curve1=2, curve2=5, &
    discretization='collocation' )

  call problem_definition ( input_probdefc, mesh, problemc )


! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemc, cval, physqarr=[1] )


! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, soln, solnm1, rhsd, sol_hat )


! fill solution vector with essential boundary conditions

  sol%u = 0

! initialize solution at tn to zero

  soln%u = 0

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! fill oldvectors

  oldvectors%s(1)%p => sol_hat  ! alpha0*un+alpha1*un-1


! create system matrix of gradient/velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! create system vectors (solution and right-hand side) for conformation and
! initialize vectors with steady stress and a small random perturbation.

  call create ( problemc, solc, solcm1, rhsc )


  if ( logc == 0 ) then ! standard
    do i = 1, nmodes
      solc(1,i)%u = 1 ! initial cxx
      solc(2,i)%u = 0 ! initial cxy
      solc(3,i)%u = 1 ! initial cyy
      solc(4,i)%u = 1 ! initial czz
    end do
  else if ( logc == 1 ) then ! log scheme
    do i = 1, nmodes
      solc(1,i)%u = 0 ! initial sxx
      solc(2,i)%u = 0 ! initial sxy
      solc(3,i)%u = 0 ! initial syy
      solc(4,i)%u = 0 ! initial szz
    end do
  end if


  if ( logc == 1 ) then

!   problem definition for projected "c=exp(s)" of the log conformation s

    call create_input_probdef ( mesh, input_probdefc_proj, nvec=1, nphysq=1 )

    input_probdefc_proj%vec_elementdof(1)%a = &
        reshape ( [ 1,0,1,0,1,0,1,0,0 ], &
                     [9,1] )

    input_probdefc_proj%physq = [1]
    input_probdefc_proj%probnr = 3

    call problem_definition ( input_probdefc_proj, mesh, problemc_proj )

    call create ( problemc_proj, solc_proj, rhsc_proj )

    do i = 1, nmodes
      solc_proj(1,i)%u = 1   ! initial cxx
      solc_proj(2,i)%u = 0   ! initial cxy
      solc_proj(3,i)%u = 1   ! initial cyy
      solc_proj(4,i)%u = 1   ! initial cyy
    end do

  end if


! create system matrix for conformation problem

  call create_sysmatrix_structure_base ( sysmatrixc, mesh, problemc )
  call create_sysmatrix_structure_constraint ( sysmatrixc, mesh, problemc )
  call finalize_sysmatrix_structure ( sysmatrixc )

  call create_sysmatrix_data ( sysmatrixc )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=1, nsysvec2=3, nprob=3 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%s2(2)%p => solcm1
  oldvectors_ve%s2(3)%p => solc_proj
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc
  oldvectors_ve%p(3)%p => problemc_proj


  if ( logc == 1 ) then

!   create and build system matrix for projection problem
!   NOTE matrix remains constant and needs to be build once.

    call create_sysmatrix_structure ( sysmatrixc_proj, mesh, problemc_proj, &
      symmetric=.true. )
    call create_sysmatrix_data ( sysmatrixc_proj )

    call build_system ( mesh, problemc_proj, sysmatrixc_proj, &
      m2sysvector=rhsc_proj, elemsub=exps_projection_elem, &
      oldvectors=oldvectors_ve, coefficients=coefficients, &
      buildvector=.false. )

    call check ( sysmatrixc_proj )

  end if


! open monitor file

  open ( unit=13, file='out', recl=300 )


! time stepping

  do step = 1, numtimesteps

    if ( step >= 2 ) then
      coefficients%i(22) = timeint2
    end if

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG


!   instationary inertia term (convection not included)

    if ( inertia == 1 ) then

!     instationary term ( rho du/dt )

      if ( step == 1 ) then
        sol_hat%u(vel%s) = soln%u(vel%s)
        mfac = 1
      else
        sol_hat%u(vel%s) = alpha0 * soln%u(vel%s) + alpha1 * solnm1%u(vel%s)
        mfac = gamma0
      end if

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=inertia_elem_dudt, coefficients=coefficients, &
        oldvectors=oldvectors, physqcol=[physqvel], physqrow=[physqvel], &
        factormat=mfac, addmatvec=.true. )

    end if


!   exps projection (for log conformation)

    if ( logc == 1 ) call solve_exps_projection


!   build implicit terms of CE with rhs in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      addmatvec=.true., coefficients=coefficients )


    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )


!   build (assemble) matrix and vector for conformation problem

    if ( coefficients%i(22) == timeint1 ) then

      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    else

      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem_implicit_2nd_order, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    end if


!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problemc, sysmatrixc, &
      m2sysvector=rhsc, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

    call check ( sysmatrixc )

    call copy ( solc, solcm1 )


!   solve conformation and keep LU decomposition in the loop over components

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do m = 1, nmodes
      do icomp = 1, ncompc
        call solve_system_ma41 ( sysmatrixc, rhsc(icomp,m), solc(icomp,m), &
          luc, solver_options=solver_options_c  )
      end do
    end do

    call delete ( luc )  ! remove LU decomposition and rebuild next time step

!   copy data for the next time step

    call copy ( soln, solnm1 )
    call copy ( sol, soln )

!   write monitor data

    write(unit=13,fmt=*) &
      step, step*deltat, &
      sol%u(press%s(1))-sol%u(press%s(2)), &  ! pressure difference
      minval(solc(1,1)%u(cval%s)), maxval(solc(1,1)%u(cval%s)), &
      maxval(sol%u(velx%s))

  end do


! close monitor data file

  close(unit=13)


! write data for post-processing

  call write_mesh ( mesh, filename='mesh.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh, input_probdefc, filename='probdefc.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol%u
  write(10) ((solc(i,j)%u, i=1,ncompc), j=1,nmodes)

  close(unit=10)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, soln, solnm1, rhsd )
  call delete ( sysmatrix )
  call delete ( oldvectors_ve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( solc, solcm1, rhsc )

  if ( logc == 1 ) then
    call delete ( sysmatrixc_proj )
    call delete ( problemc_proj )
    call delete ( input_probdefc_proj )
    call delete ( solc_proj, rhsc_proj )
    call delete ( lu_exps_proj )
  end if

  call delete ( coefficients )
  call delete ( cval )

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

!   periodical condition on velocities and gradients

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_elem_conn, addmatvec=.true., &
      coefficients=coefficients )

!   periodical condition on gradients

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true. )

!   imposed flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=3, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients )

  end subroutine build_vpG


! project c=exp(s) on discrete fem space

  subroutine solve_exps_projection

    type(solver_options_ma57_t) :: solver_options_ma57

    integer :: i, m

!   build vector only (matrix is constant)

    call build_system ( mesh, problemc_proj, sysmatrixc_proj, &
      m2sysvector=rhsc_proj, elemsub=exps_projection_elem, &
      oldvectors=oldvectors_ve, coefficients=coefficients, &
      buildmatrix=.false. )

    ! MA57 solver storage
    solver_options_ma57%integer_storage = 1.3
    solver_options_ma57%real_storage    = 1.3

!   LU decomposition is done in the first call only

    do m = 1, nmodes
      do i = 1, ncompc

        call add_effect_of_essential_to_rhs ( problemc_proj, sysmatrixc_proj, &
           solc_proj(i,m), rhsc_proj(i,m) )

        call solve_system_ma57 ( sysmatrixc_proj, rhsc_proj(i,m), &
           solc_proj(i,m), lu_exps_proj, solver_options=solver_options_ma57 )

      end do
    end do

  end subroutine solve_exps_projection


end program poiseuille4
