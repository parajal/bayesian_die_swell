! Startup of flow around a rotating cylinder in a square box.
! Compressible EGP model with compressible flow.
! 14-mode model with 12 alpha process viscoelastic modes, one elastic mode for
! strain hardening and one mode for the beta process.
! Data for PMMA from van Breemen et al. (2012) DOI: 10.1002/polb.23199
! Pressure (J) dependent relaxation term using adapted lambda.
! Plastic strain softening.
! Contravariant deformation tensor formulation (b-tensor).
! Fully-implicit with Newton-Raphson.
! Compute drag and torque on cylinder using reaction forces.
! Optionally integral pressure constraint imposed.
! Optionally impose pressure in a single point (P1).
! Optionally after the final time step do an iteration to steady-state.
!   This option is not applicable to the model used here because of the
!   elastic mode, which gives solid material behavior.
! Optionally include Jacobian of SUPG test function.
! Optionally discard certain zero blocks from the sparse matrix.

program cylinder_egp7_hJ_gammap

  use tfem_m
  use viscoelastic_elements_m
  use compressible_fluid_elements_m
  use hsl_ma41_m
  use io_utils_m
!  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities
    pintpl = 4,     & ! Q1 pressures
    gintpl = 4,     & ! Q1 gradients
    bintpl = 4,     & ! Q1 b-tensor
    gauss = 3,      & ! 3x3 integration of quads
    gaussb = 3,     & ! 3 point integration of boundary elements
    timeint1 = 8,   & ! BDF1, first time step
    timeint2 = 10,  & ! BDF2, after first time step
    timeint3 = 10,  & ! BDF2, for steady time discretization
    numtimesteps = 1, & ! number of time steps
    maxnumiterations = 20, & ! maximum number of Newton-Raphson iterations
    ncompb = 5,     & ! number of b-tensor components
    nmodes = 4,    & ! number of modes
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    physqb = 4,     & ! physical quantity nr of the b-tensor (mode 1)
    physqgammap = 4+nmodes,& ! physical quantity nr of the equiv. plastic strain
    am = nmodes-2,  & ! 1:am range of modes for alpha process
    em = nmodes-1,  & ! elastic mode
    bm = nmodes,    & ! mode of beta process
    startm = 501,   & ! start of material model data
    model = 24,     & ! EGP compressible
    alamJ_model = 1,& ! adapted lambda model for h(J)
    alam_gammap_model = 1,& ! adapted lambda model for g(gammap)
    alam_model = 4    ! adapted lambda model: 4: Eyring

  real(dp), dimension(nmodes) :: &
    G,       & ! modulus
    lambda,  & ! relaxation time
    tau_ref, &  ! reference von Mises
    Sa,      & ! state parameter Sa
    r0,      & ! fitting parameter r0
    r1,      & ! fitting parameter r1
    r2         ! fitting parameter r2

  real(dp) :: &
    eta_p,    & ! polymer viscosity
    eta_s,    & ! solvent viscosity
    eta_0       ! zero-shear viscosity
!    lambda_0    ! effective relaxation time for defining Wi

  real(dp), parameter :: &
    Kmod = 3.e3_dp,    & ! compression modulus
    p0 = 0._dp,        & ! reference pressure where J=1.
    beta_s = 0.0_dp,   & ! beta viscosity parameter = eta_s/(eta_s+eta_p)
    betaJ = 237._dp,   & ! beta exponent for the h(J) factor = mu*K/tau_ref
!    H = 1.5_dp,        & ! (half-) height of the box (must match mesh)
    R = 1._dp,         & ! radius of the cylinder (must match mesh)
    omega = 1.e-3_dp       ! angular velocity of the cylinder

  real(dp), parameter :: &
    deltat = 0.1_dp,    & ! time step for BDF1/BDF2
    deltat3 = 0.1_dp,   & ! time step for steady iteration scheme
    epsconf = 1e-8_dp,  & ! Newton-Raphson convergence threshold
    beta_supg = 1,      & ! SUPG factor
    rs = 1.2_dp,        & ! real_storage for LU (HSL)
    is = 1.6_dp           ! integer_storage for LU (HSL)

  logical, parameter :: &
    steadylaststep = .false., & ! iterate to steady-state solution in last step
    pressure_constraint = .false., & ! impose pressure constraint
    check_pressure_constraint = .false., & ! check pressure constraint
    pressure_essential = .false., & ! p=p0 in point P1 (lower-left corner)
    physqmask = .false., & ! use physical quantity masking
    physqmask_grad = .true., & ! use physqmask for G-b coupling
    physqmask_press = .true., & ! use physqmask for p-b and b-p coupling
    physqmask_offdiag = .false., & ! use physqmask for b-b off-diagonal coupling
    JacobianSUPG = .true., & ! take Jacobian of SUPG test function into account
    printscreen = .true.  ! print on standard output

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: soln, solnm1, soliter
  type(sysvector_t) :: dsol, rhsd, reacf
!  type(vector_t) :: velocity, pressure, vorticity, Jvol
  type(subscript_t) :: velx, vely, vel, pres, solsc, velxcylinder, &
    velycylinder, gp
  type(subscript_t), dimension(nmodes) :: bval
  type(subscript_t), dimension(ncompb,nmodes) :: bvalcmp
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options

  type(subscriptvec_t) :: cxx, cxy, cyy, czz
  type(vector_t) :: ctensor, btensor!, vonmises

  logical :: steadystep
  integer :: step, i, iter, m, extrastep, mode1, mode2, npar
  integer :: vertices(4) = [1,3,5,7]
  real(dp) :: alpha, epsu, epsp, epsb, epsgp, K_dragx, K_dragy, K_torqz
  real(dp) :: int_press_domain(1), V0(1)

!  timer = .false.

! set material parameters

! alpha process
  G(1:am) = &
    [ 1.08e1_dp, 5.69e2_dp ]
  lambda(1:am) = &
    [ 1.73e4_dp,  1.55e4_dp ]
  tau_ref(1:am) = 2.53_dp ! tau_ref for alpha
  Sa(1:am) = 7.4_dp       ! Sa for alpha
  r0(1:am) = 0.96_dp      ! r0 for alpha
  r1(1:am) = 20.0_dp      ! r1 for alpha
  r2(1:am) = -2.0_dp      ! r2 for alpha

! elastic mode
  G(em) = 2.10e1_dp
  lambda(em) = 1.0_dp ! lambda for determinant stabilization only
  tau_ref(em) = 0._dp ! dummy tau_ref value for elastic mode
  Sa(em) = 0._dp      ! dummy Sa for elastic mode
  r0(em) = 0._dp      ! dummy r0 for elastic mode
  r1(em) = 0._dp      ! dummy r1 for elastic mode
  r2(em) = 0._dp      ! dummy r2 for elastic mode

! beta process
  G(bm) = 1.13e2_dp
  lambda(bm) = 3.65e-3_dp
  tau_ref(bm) = 2.71_dp ! tau_ref for beta
  Sa(bm) = 7.4_dp       ! Sa for alpha
  r0(bm) = 0.96_dp      ! r0 for alpha
  r1(bm) = 20.0_dp      ! r1 for alpha
  r2(bm) = -2.0_dp      ! r2 for alpha

! viscosity parameters
  eta_p = sum( G(1:am)*lambda(1:am) + G(bm)*lambda(bm) ) ! polymer viscosity
  eta_s = beta_s/(1-beta_s)*eta_p       ! solvent viscosity
  eta_0 = eta_s+eta_p    ! zero-shear viscosity

! effective lambda
  !lambda_0 = sum( G(1:am)*lambda(1:am)**2 + G(bm)*lambda(bm)**2 ) / eta_p

!  print *, 'eta_p, eta_s, eta_0, lambda_0'
!  print *, eta_p, eta_s, eta_0, lambda_0

! set other parameters

!  alpha = G * lambda   ! DEVSS parameter
  alpha = sum(G)*deltat   ! DEVSS parameter

! fill coefficients

  if ( alamJ_model == 1 ) then
    npar = 8
  else
    npar = 7
  end if

  call create_coefficients ( coefficients, ncoefi=350, ncoefr=500+npar*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
       physqvel, physqpress, 0,     physqgrad, gauss,  &
       gaussb,   bintpl,     0,     0,         0,      &
       0,        0,          model, nmodes,    startm, &
       0,     timeint1,   ( 0, i = 23, 350 )  &
    ]

  coefficients%i(39) = 1 ! do not build the continuity equation
  coefficients%i(71) = 1  ! CDT for b-formulation
  coefficients%i(80) = alam_model ! adapted lambda
  coefficients%i(83) = physqb ! physical quantity of first mode
  coefficients%i(84) = 3  ! storage of b tensor: all modes in sysvector
  coefficients%i(86) = 2  ! sysvector number for velocity in SUPG = soln
  coefficients%i(87) = 1  ! sysvector number for iteration b
  coefficients%i(88) = 2  ! sysvector number of b at tn
  coefficients%i(90) = 1  ! rotation reinitialization on element level

  if ( JacobianSUPG ) then
    coefficients%i(86) = 1  ! sysvector number for velocity in SUPG = soliter
    coefficients%i(91) = 1  ! Jacobian of SUPG
  end if

  coefficients%i(92) = 1  ! global von Mises
  coefficients%i(93) = em ! elastic mode (norlxmode)
  coefficients%i(96) = alamJ_model ! adapted lambda for J
  coefficients%i(97) = physqgammap ! physical quantity of equiv. plastic strain
  coefficients%i(98) = alam_gammap_model ! adapted lambda for g(gammap)
  coefficients%i(99) = 1  ! mode number used in computing gammap (max lambda)
  coefficients%i(100) = 2 ! storage gammap in sysvector (multiple quantities)
  coefficients%i(101) = 1 ! sysvector number for iteration gammap
  coefficients%i(102) = 2 ! sysvector number of gammap at tn
  coefficients%i(103) = 1 ! Include the Jacobian delta u.grad(gammap)

  coefficients%i(301) = 3 ! start with first order time integration
  coefficients%i(302) = 2 ! position sol_n in oldvectors
  coefficients%i(305) = 1 ! sysvector number for iteration velocity/pressure
  if ( pressure_constraint ) then
    coefficients%i(307) = 1 ! constraint number of the pressure constraint
  end if
  coefficients%i(308) = 1 ! type of pressure constraint (1,2,3 or 4)


  coefficients%r(1:500) = &
    [ eta_s,    0._dp,  0._dp,     alpha, 0._dp, &
         0._dp, 0._dp, deltat, beta_supg, 0._dp, &
       ( 0._dp, i = 11, 500 ) ]

  coefficients%r(29) = Kmod
  coefficients%r(30) = p0

  if ( alamJ_model == 1 ) then
    coefficients%r(501:500+nmodes*8) = &
                        [ ( G(m), lambda(m), tau_ref(m), betaJ, &
                            Sa(m), r0(m), r1(m), r2(m), m=1,nmodes ) ]
  else
    coefficients%r(501:500+nmodes*7) = &
                        [ ( G(m), lambda(m), tau_ref(m), &
                            Sa(m), r0(m), r1(m), r2(m), m=1,nmodes ) ]
  end if

  coefficients%r(10) = omega * R ! global scaling velocity for SUPG, when needed

  coefficients%r(251) = deltat
  coefficients%r(252) = Kmod
  coefficients%r(253) = p0

! read mesh

  call read_mesh ( mesh, filename='mesh_cyl_in_box.out' )

  if ( pressure_constraint ) then

!   add elementset for integral pressure constraint

    call add_to_mesh ( mesh, elementset='elements', &
                                 elements=[(i,i=1,mesh%nelem)] )
    call add_to_mesh ( mesh, elementset='nodes', elementsetnr=1 )

  end if

  call fill_mesh_parts ( mesh )

!  call printinfo ( mesh, printlevel=1 )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=7+nmodes, &
    nphysq=4+nmodes )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 4  ! G
  input_probdef%vec_elementdof(1)%a(:,2) = 2         ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(vertices,4:3+nmodes) = ncompb  ! b
  input_probdef%vec_elementdof(1)%a(vertices,4+nmodes:4+nmodes) = 1  ! gammap
  input_probdef%vec_elementdof(1)%a(:,5+nmodes) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,6+nmodes) = 4  ! tensor for plotting
  input_probdef%vec_elementdof(1)%a(:,7+nmodes) = 5  ! tensor for plotting

  input_probdef%physq = [1,2,3,(3+m,m=1,nmodes),4+nmodes]
  input_probdef%probnr = 1

  if ( physqmask ) then

!   discard certain zero blocks from the sparse matrix

    if ( physqmask_grad ) call fill_grad_physqmask         ! G-b block
    if ( physqmask_press ) call fill_press_physqmask       ! p-b and b-p blocks
    if ( physqmask_offdiag ) call fill_offdiag_physqmask   ! off-diagonal b

  end if

! define essential boundaries

! cylinder
  call define_essential ( mesh, input_probdef, curve1=13, physq=physqvel )
! walls
  call define_essential ( mesh, input_probdef, curves=[4,7,10,12], &
    physq=physqvel )

  if ( pressure_essential ) then

!   set pressure level in point P1

    call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

  end if

  if ( pressure_constraint ) then

!   constraint for integral pressure constraint

    call define_constraint ( mesh, input_probdef, &
      physq=physqpress, elementset1=1, nglobalc=1 )

  end if

  call problem_definition ( input_probdef, mesh, problem )

 ! create vector subscripts for solution (excluding contraint forces)

  call create_subscript ( mesh, problem, solsc, &
    physqarr = &
      [physqgrad,physqvel,physqpress,physqgammap,(physqb+m-1,m=1,nmodes)] )

! create vector subscripts for the velocity

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )
! velocity degrees in x-direction on the cylinder
  call create_subscript ( mesh, problem, velxcylinder, physqarr=[physqvel], &
    degfd=1, curves=[13] )
! velocity degrees in y-direction on the cylinder
  call create_subscript ( mesh, problem, velycylinder, physqarr=[physqvel], &
    degfd=2, curves=[13], fillnodes=.true. )

! create vector subscript for the pressure

  call create_subscript ( mesh, problem, pres, physqarr=[physqpress] )

! create vector subscript for the conformation tensor

  do m = 1, nmodes
    call create_subscript ( mesh, problem, bval(m), physqarr=[physqb+m-1] )
    do i = 1, ncompb
      call create_subscript ( mesh, problem, bvalcmp(i,m), &
        physqarr=[physqb+m-1], degfd=i )
    end do
  end do
  call create_subscript ( mesh, problem, cxx, degfd=1, vec=6+nmodes )
  call create_subscript ( mesh, problem, cxy, degfd=2, vec=6+nmodes )
  call create_subscript ( mesh, problem, cyy, degfd=3, vec=6+nmodes )
  call create_subscript ( mesh, problem, czz, degfd=4, vec=6+nmodes )

! create vector subscript for the equivalent plastic strain

  call create_subscript ( mesh, problem, gp, physqarr=[physqgammap] )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, soln, solnm1, soliter )
  call create_sysvector ( problem, dsol, rhsd, reacf )


! create a vector for conformation
! tensor for post processing

  call create_vector ( problem, ctensor, vec=6+nmodes )
  call create_vector ( problem, btensor, vec=7+nmodes )

! fill solution vector with essential boundary conditions

  soln%u = 0

! initial pressure

  soln%u(pres%s) = p0

! initial solution b tensor

  do m = 1, nmodes
    soln%u(bvalcmp(1,m)%s) = 1 ! initial bxx
    soln%u(bvalcmp(2,m)%s) = 0 ! initial bxy
    soln%u(bvalcmp(3,m)%s) = 0 ! initial byx
    soln%u(bvalcmp(4,m)%s) = 1 ! initial byy
    soln%u(bvalcmp(5,m)%s) = 1 ! initial bzz
  end do

  call copy(soln,solnm1)
  call copy(soln,soliter)

  dsol%u = 0

! set velocity_theta = omega*R on cylinder and pressure level in first step
  call fill_sysvector ( mesh, problem, dsol, &
    curve1=13, physq=physqvel, vfunc=vfunc, vfuncnr=1 )


! create system matrix

  if ( pressure_constraint ) then

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      usephysqmask=physqmask )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

  else

    call create_sysmatrix_structure ( sysmatrix, mesh, problem, &
      usephysqmask=physqmask )

  end if

  call create_sysmatrix_data ( sysmatrix )

  print *, 'nnz = ', sysmatrix%Suu%nnz


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=3, nprob=1 )


! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => soliter
  oldvectors_ve%s(2)%p => soln
  oldvectors_ve%s(3)%p => solnm1
  oldvectors_ve%p(1)%p => problem  ! needed for implicit_ce_supg_elem


! Weisenberg number
  !Wi = lambda_0 * omega * R / (H-R)

!  print *, ' omega * R = ', omega * R, 'Wi = ', Wi
!  print *


!  open ( unit=13, file='iter.out', recl=300 )
!  open ( unit=14, file='Kdrag.out', recl=300 )

  solver_options%real_storage=rs
  solver_options%integer_storage=is

! area of the (original) domain

  call integrate ( mesh, problem, V0, elemsub=stokes_integrate_volume, &
    coefficients=coefficients )

  coefficients%r(254) = V0(1)

! perform a steady step?

  if ( steadylaststep ) then
    extrastep = 1
  else
    extrastep = 0
  end if

!  call tic

! time stepping

  do step = 1, numtimesteps + extrastep  ! extra step for steady-state

    if ( step <= numtimesteps ) then

!     regular time stepping

      steadystep = .false.

!      write(13,*) 'step = ', step
      if ( printscreen) print *, 'step = ', step

      if ( step == 1 ) then
        coefficients%i(78) = 1 ! check for zero velocity in SUPG
      else
        coefficients%i(78) = 0
      end if

      if ( step >= 2 ) then
        coefficients%i(22) = timeint2
        coefficients%i(301) = 4 ! second order time integration
      end if

    else

!     steady-state iteration

      steadystep = .true.

 !     write(13,*) 'step = ', step, 'steady-state iteration'
!      if ( printscreen) print *, 'step = ', step, 'steady-state iteration'

      if ( numtimesteps == 0 ) then
        coefficients%i(78) = 1 ! check for zero velocity in SUPG
      else
        coefficients%i(78) = 0
      end if

      coefficients%i(22) = timeint3
      coefficients%r(8) = deltat3

!     exclude time-derivative in implicit_pressure_elem
      coefficients%i(306) = 1

    end if

    iter = 0

    do

      iter = iter + 1

      if ( steadystep ) then
        if ( numtimesteps == 0 .and. iter == 2 ) then
 !        set dsol = 0 on cylinder
          call fill_sysvector ( mesh, problem, dsol, &
            curve1=13, physq=physqvel, value=0._dp )
        end if
      else
        if ( step == 1 .and. iter == 2 ) then
 !        set dsol = 0 on cylinder
          call fill_sysvector ( mesh, problem, dsol, &
            curve1=13, physq=physqvel, value=0._dp )
        end if
      end if

      if ( iter > maxnumiterations ) then
        write(*,'(3(a,i0/))') &
          ' Maximum number of iterations reached = ', &
          maxnumiterations, ' step = ', step
        stop
      end if

!     build (assemble) matrix/vector for gradient/velocity/pressure part

      call build_vpG

!      call toc ( 'build_vpG' )

!     build (assemble) matrix/vector for conformation part

!     alpha process of modes 1 to am
      mode1 = 1; mode2 = am
      call build_b

!     elastic mode em
      mode1 = em; mode2 = em
      call build_b

!     beta process mode bm
      mode1 = bm; mode2 = bm
      call build_b

      !call toc ( 'build_b' )

      if ( .not. physqmask .or. &
              physqmask .and. .not. physqmask_offdiag ) then

!       set to zero off-diagonal conformation blocks

        call build_offdiag_b

        !call toc ( 'build_offdiag_b' )

      end if

!     build (assemble) matrix/vector for gammap part

!     use alpha process of modes 1 to am for vonmises stress
      mode1 = 1; mode2 = am
      call build_gammap

      !call toc ( 'build_gammap' )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, dsol, rhsd )

      call check ( sysmatrix )

!     solve

      call solve_system_ma41 ( sysmatrix, rhsd, dsol, &
        solver_options=solver_options  )

      !call toc ( 'solve' )

      soliter%u(solsc%s) = soliter%u(solsc%s) + dsol%u(solsc%s)

!     convergence test

      epsu = maxval(abs(dsol%u(vel%s))) / ( omega * R )
      epsp = maxval(abs(dsol%u(pres%s))) * R / ( ( omega * R ) * eta_0 )
      epsb = 0
      do m = 1, nmodes
        epsb = max ( epsb, maxval(abs(dsol%u(bval(m)%s))) )
      end do
      epsgp = maxval(abs(dsol%u(gp%s)))

      if ( printscreen) print *, iter, epsu, epsp, epsb, epsgp
!      write(*,*) iter, epsu, epsp, epsb, epsgp

      if ( maxval([ epsu, epsp, epsb, epsgp ]) < epsconf ) exit

    end do

    !call toc ( 'after iteration loop' )


    if ( pressure_constraint .and. check_pressure_constraint ) then

!     check pressure constraint

      call integrate ( mesh, problem, int_press_domain, &
        elemsub=integrate_pressure_function, &
        coefficients=coefficients, oldvectors=oldvectors_ve )

      if ( printscreen) write(*,*) 'constraint error =', int_press_domain(1)

    end if


!   reaction forces

    call reaction_forces ( problem, sysmatrix, dsol, rhsd, reacf )

    K_dragx = - sum ( reacf%u(velxcylinder%s) ) / ( eta_0*( omega * R ) )
    K_dragy = - sum ( reacf%u(velycylinder%s) ) / ( eta_0*( omega * R ) )
    K_torqz = - sum ( mesh%coor(velycylinder%nodes,1) * &
                            reacf%u(velycylinder%s) - &
                      mesh%coor(velycylinder%nodes,2) * &
                            reacf%u(velxcylinder%s) ) &
                                     / ( eta_0*( omega * R**2 ) )

    if ( printscreen) print *, 'Kdragx = ', K_dragx, 'Kdragy = ', K_dragy, &
      'Ktorqz = ', K_torqz
!    write(14,*) step, step * deltat, K_dragx, K_dragy, K_torqz

!    call toc ( 'after reaction forces' )


!   copy solution to older time step for next time step

    call copy ( soln, solnm1 )
    call copy ( soliter, soln )

!   write max and mean values of conformation and b-tensor tensor to a file

    call derive_vector ( mesh, problem, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

!    if ( step == 1 ) then
!      open(unit=11, recl=600, status='replace', file='cval.out')
!    else
!      open(unit=11, recl=600, position='append', file='cval.out')
!    end if
    write(*, *) step * deltat, sum(ctensor%u(cxx%s))/size(cxx%s), &
                                    sum(ctensor%u(cxy%s))/size(cxy%s), &
                                    sum(ctensor%u(cyy%s))/size(cyy%s), &
                                    sum(ctensor%u(czz%s))/size(czz%s), &
                                    maxval(soln%u(gp%s)), &
                                    minval(soln%u(gp%s))
!    close(unit=11)

!    if ( printscreen ) &
!         print *, 'step = ', step, 'max cxx = ', maxval(ctensor%u(cxx%s))

!    call toc ( 'one step' )

  end do

!  call toc ( 'all steps' )



! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( soln, solnm1, soliter, dsol, rhsd )
  call delete ( reacf )
  call delete ( sysmatrix )

  call delete ( oldvectors_ve )
  do m = 1, nmodes
    do i = 1, ncompb
      call delete ( bvalcmp(i,m), bval(m) )
    end do
  end do
  call delete ( cxx, cxy, cyy, czz )
  call delete ( gp )

contains

  function vfunc ( n, nr, x )
    use kind_defs_m
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    vfunc = omega * [ -x(2), x(1) ] / norm2(x)

  end function vfunc

  subroutine fill_grad_physqmask

    integer :: j

    do j = 1, nmodes
      input_probdef%physqmask(physqgrad,physqb+j-1) = .false.
    end do

  end subroutine fill_grad_physqmask


  subroutine fill_press_physqmask

    integer :: j

    do j = 1, nmodes
      input_probdef%physqmask(physqpress,physqb+j-1) = .false.
    end do

  end subroutine fill_press_physqmask


  subroutine fill_offdiag_physqmask

    integer :: i, j

    do i = em, bm
      do j = 1, nmodes
        if ( i == j ) cycle
          input_probdef%physqmask(physqb+i-1,physqb+j-1) = .false.
          input_probdef%physqmask(physqb+j-1,physqb+i-1) = .false.
      end do
    end do

  end subroutine fill_offdiag_physqmask


  subroutine build_vpG

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients, &
      physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress] )

!   pressure evolution equation
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_pressure_ce_elem, oldvectors=oldvectors_ve, &
      physqrow=[physqpress], physqcol=[physqvel,physqpress], &
      coefficients=coefficients, addmatvec=.true. )

!   DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_elem, coefficients=coefficients, addmatvec=.true., &
      physqrow=[physqgrad,physqvel], physqcol=[physqgrad,physqvel] )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqgrad], physqcol=[physqpress], &
      zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqpress], physqcol=[physqgrad], &
      zeromatvec=.true. )

!   build DEVSS-G right-hand side

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_rhs_div, oldvectors=oldvectors_ve, &
      coefficients=coefficients, addmatvec=.true., &
      buildmatrix=.false., physqrow=[physqvel] )

!   build stokes right-hand side

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_rhs_divsigma, oldvectors=oldvectors_ve, &
      coefficients=coefficients, addmatvec=.true., &
      buildmatrix=.false., physqrow=[physqvel] )

!   build -div(tau) part in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_divtau, oldvectors=oldvectors_ve, &
      physqrow=[physqvel], physqcol=[(physqb+m-1,m=1,nmodes)], &
      addmatvec=.true., coefficients=coefficients )

!   build -div(dtau/dJ) part in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_divtau_J, oldvectors=oldvectors_ve, &
      physqrow=[physqvel], physqcol=[physqpress], &
      addmatvec=.true., coefficients=coefficients, buildvector=.false. )

    if ( pressure_constraint ) then

!     build pressure contraint

!     pressure-pressure part of the Jacobian
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=pp_constraint_Jacobian_elem, coefficients=coefficients, &
        physqrow=[physqpress], physqcol=[physqpress], addmatvec=.true., &
        oldvectors=oldvectors_ve, buildvector=.false. )

!     constraint equation
      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        constraint1=1, elemsub=constraint_pressure_int_elementset, &
        addmatvec=.true., coefficients=coefficients, &
        oldvectors=oldvectors_ve )

    end if

  end subroutine build_vpG


  subroutine build_b

!   set modes of subprocess
    coefficients%i(26:27) = [ mode1, mode2 ]

!   build (assemble) matrix and vector for b-tensor mode m

    if ( steadystep ) then

!     exclude time-derivative in implicit_ce_supg_elem
      coefficients%i(89) = 1

      if ( JacobianSUPG ) then

!       add velocity and diagonal block of steady-state timed derivative
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=steady_ce_timederiv_supg_elem, coefficients=coefficients, &
          physqrow=[(physqb+m-1,m=mode1,mode2)], &
          physqcol=[physqvel,(physqb+m-1,m=mode1,mode2)], &
          oldvectors=oldvectors_ve, addmatvec=.true. )

      else

!       add diagonal block of steady-state timed derivative
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=steady_ce_timederiv_supg_elem, coefficients=coefficients, &
          physqrow=[(physqb+m-1,m=mode1,mode2)], &
          physqcol=[(physqb+m-1,m=mode1,mode2)], &
          oldvectors=oldvectors_ve, addmatvec=.true. )

      end if

    end if

    if ( JacobianSUPG ) then

!     pressure, velocity, b-tensor and gammap
       call build_system ( mesh, problem, sysmatrix, rhsd, &
         elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
         physqrow=[(physqb+m-1,m=mode1,mode2)], &
         physqcol=[physqpress,physqgammap,physqvel,(physqb+m-1,m=mode1,mode2)],&
         oldvectors=oldvectors_ve, addmatvec=.true. )

    else

!     pressure, b-tensor and gammap
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
        physqrow=[(physqb+m-1,m=mode1,mode2)], &
        physqcol=[physqpress,physqgammap,(physqb+m-1,m=mode1,mode2)], &
        oldvectors=oldvectors_ve, addmatvec=.true. )

    end if

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_ce_vel_supg_elem, coefficients=coefficients, &
      physqrow=[(physqb+m-1,m=mode1,mode2)], physqcol=[physqgrad,physqvel], &
      oldvectors=oldvectors_ve, addmatvec=.true., buildvector=.false. )

    if ( .not. physqmask .or. &
            physqmask .and. .not. physqmask_grad ) then
!     set to zero off-diagonal blocks gradient-b
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        physqrow=[physqgrad], physqcol=[(physqb+m-1,m=mode1,mode2)], &
        buildvector=.false., zeromatvec=.true. )
    end if

    if ( .not. physqmask .or. &
            physqmask .and. .not. physqmask_press ) then
!     set to zero off-diagonal block pressure-b
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        physqrow=[physqpress], physqcol=[(physqb+m-1,m=mode1,mode2)], &
        buildvector=.false., zeromatvec=.true. )
    end if

!   reset modes to default
    coefficients%i(26:27) = 0

  end subroutine build_b

  subroutine build_offdiag_b

    integer :: i, j

    do i = em, bm
      do j = 1, nmodes
        if ( i == j ) cycle
          call build_system ( mesh, problem, sysmatrix, rhsd, &
            addmatvec=.true., physqrow=[physqb+i-1], physqcol=[physqb+j-1], &
            buildvector=.false., zeromatvec=.true. )
          call build_system ( mesh, problem, sysmatrix, rhsd, &
            addmatvec=.true., physqrow=[physqb+j-1], physqcol=[physqb+i-1], &
            buildvector=.false., zeromatvec=.true. )
      end do
    end do

  end subroutine build_offdiag_b

  subroutine build_gammap

!   set modes of subprocess
    coefficients%i(26:27) = [ mode1, mode2 ]

!   build (assemble) matrix and vector for gammap

    if ( steadystep ) then

!     exclude time-derivative in implicit_gammap_supg_elem
      coefficients%i(89) = 1

    end if

    if ( JacobianSUPG ) then

!     pressure, velocity, b-tensor and gammap
       call build_system ( mesh, problem, sysmatrix, rhsd, &
         elemsub=implicit_gammap_supg_elem, coefficients=coefficients, &
         physqrow=[physqgammap], &
         physqcol=[physqpress,physqvel,(physqb+m-1,m=mode1,mode2),physqgammap],&
         oldvectors=oldvectors_ve, addmatvec=.true. )

    else

!     pressure, b-tensor and gammap
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=implicit_gammap_supg_elem, coefficients=coefficients, &
        physqrow=[physqgammap], &
        physqcol=[physqpress,(physqb+m-1,m=mode1,mode2),physqgammap], &
        oldvectors=oldvectors_ve, addmatvec=.true. )

    end if

!   set to zero off-diagonal blocks gammap-gradient, zero b
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      physqrow=[physqgammap], physqcol=[physqgrad, &
                  (physqb+m-1,m=1,mode1-1),(physqb+m-1,m=mode2+1,nmodes)], &
      buildvector=.false., zeromatvec=.true. )

!   set to zero off-diagonal blocks gradient,velocity,pressure-gammap
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      physqrow=[physqgrad,physqvel,physqpress], physqcol=[physqgammap], &
      buildvector=.false., zeromatvec=.true. )

!   reset modes to default
    coefficients%i(26:27) = 0

  end subroutine build_gammap

end program cylinder_egp7_hJ_gammap
