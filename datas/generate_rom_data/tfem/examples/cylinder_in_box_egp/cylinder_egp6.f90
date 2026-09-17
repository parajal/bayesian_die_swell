! Sequence of steady states of flow around a rotating cylinder in a square box.
! Compressible EGP model with compressible flow.
! Contravariant deformation tensor formulation (b-tensor).
! Iteration to steady-state with Newton-Raphson.
! Compute drag and torque on cylinder using reaction forces.
! Optionally integral pressure constraint imposed.
! Optionally impose pressure in a single point (P1).
! Optionally include Jacobian of SUPG test function.
! Optionally discard certain zero blocks from the sparse matrix.

program cylinder_egp6

  use tfem_m
  use viscoelastic_elements_m
  use compressible_fluid_elements_m
  use hsl_ma41_m
  use io_utils_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities
    pintpl = 4,     & ! Q1 pressures
    gintpl = 4,     & ! Q1 gradients
    bintpl = 4,     & ! Q1 b-tensor
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    physqb = 4,     & ! physical quantity nr of the b-tensor (mode 1)
    gauss = 3,      & ! 3x3 integration of quads
    gaussb = 3,     & ! 3 point integration of boundary elements
    timeint3 = 10,  & ! BDF2, for steady time discretization
    maxnumiterations = 30, & ! maximum number of Newton-Raphson iterations
    ncompb = 5,     & ! number of b-tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 24,     & ! EGP compressible
    alam_model = 4    ! adapted lambda model: 4: Eyring, 5: Ree-Eyring

  real(dp), parameter :: &
    G     = 4.0_dp,    & ! modulus
    lambda = 1.0_dp,   & ! relaxation time
    Kmod = 1.e2_dp,    & ! compression modulus
    p0 = 0._dp,        & ! reference pressure where J=1.
    eta_p = lambda*G,  & ! polymer viscosity
    beta_s = 0.1_dp,   & ! beta viscosity parameter = eta_s/(eta_s+eta_p)
    eta_s = beta_s/(1-beta_s)*eta_p,    & ! solvent viscosity
    eta_0 = eta_s+eta_p, & ! zero-shear viscosity
    tau_ref = 1.0_dp,    & ! reference von Mises
    tau_ref1 = 1.0_dp,   & ! reference von Mises 1 RE
    tau_ref2 = 4.0_dp,   & ! reference von Mises 2 RE
    f1 = 0.6_dp,         & ! weight factor for first term RE
    H = 1.5_dp,          & ! (half-) height of the box (must match mesh)
    R = 1._dp            ! radius of the cylinder (must match mesh)

  real(dp), parameter :: &
    deltat3 = 0.05_dp,   & ! time step for steady iteration scheme
    beta_supg = 1,       & ! SUPG factor
    rs = 2.0_dp,         & ! real_storage for LU (HSL)
    is = 1.6_dp            ! integer_storage for LU (HSL)

  integer, parameter :: &
    ndU = 1           ! number intervals with different dU increment

  real(dp), parameter :: &
!   increments in imposed angular velocity
    dU(ndU) = [ 0.05_dp ], &
!   end values of the intervals in U
    Ue(ndU) = [ 1.0_dp ], &
!   Newton-Raphson convergence threshold for each interval
    epsconfin(ndU) = [ 1e-8_dp ]

  logical, parameter :: &
    pressure_constraint = .true., & ! impose pressure constraint
    check_pressure_constraint = .false., & ! check pressure constraint
    pressure_essential = .false., & ! p=p0 in point P1 (lower-left corner)
    physqmask = .true., & ! use physical quantity masking
    physqmask_grad = .true., & ! use physqmask for G-b coupling
    physqmask_press = .false., & ! use physqmask for p-b and b-p coupling
    physqmask_offdiag = .true., & ! use physqmask for b-b off-diagonal coupling
    JacobianSUPG = .true., & ! take Jacobian of SUPG test function into account
    printscreen = .true.  ! print on standard output

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: solsupg, soliter
  type(sysvector_t) :: dsol, rhsd, reacf
  type(vector_t) :: velocity, pressure, vorticity, Jvol, divu
  type(subscript_t) :: velx, vely, vel, pres, solsc, velxcylinder, velycylinder
  type(subscript_t), dimension(nmodes) :: bval
  type(subscript_t), dimension(ncompb,nmodes) :: bvalcmp
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options

  type(subscriptvec_t) :: cxx, cxy, cyy, czz
  type(vector_t) :: ctensor, btensor, vonmises

  integer :: step, iter, m, i, nsteps(ndU), numsteps, j, npar
  integer :: vertices(4) = [1,3,5,7]
  real(dp) :: alpha, domega, epsconf, omega, K_dragx, K_dragy, K_torqz
  real(dp) :: epsu, epsp, epsb, Wi
  real(dp), dimension(:), allocatable :: Uspec, epsc
  real(dp) :: int_press_domain(1), V0(1)

  timer = .false.


! number of incremental steps in omega to perform

  nsteps(1) = nint((Ue(1))/dU(1)) ! first interval
  do i = 2, ndU
    nsteps(i) = nint((Ue(i)-Ue(i-1))/dU(i)) ! second and higher intervals
  end do
  numsteps = sum(nsteps)

  if ( printscreen ) print *, 'numsteps = ', numsteps


! omega values specified

  allocate( Uspec(0:numsteps), epsc(numsteps) )

  Uspec = [ 0.0_dp, (i*dU(1), i=1,nsteps(1)), &
                    ((Ue(j-1)+dU(j)*i,i=1,nsteps(j)),j=2,ndU)  ]

  if ( printscreen ) print *, 'Uspec = ', Uspec(1:numsteps)

  epsc = [ ( ( epsconfin(j), i=1,nsteps(j) ), j=1,ndU ) ]


! set some parameters

  alpha = eta_p   ! DEVSS parameter

! fill coefficients

  if ( alam_model == 4 ) then
    npar = 3
  else if ( alam_model == 5 ) then
    npar = 5
  end if

  call create_coefficients ( coefficients, ncoefi=350, ncoefr=500+npar*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gaussb,   bintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      0,     timeint3,   ( 0, i = 23, 350 )  &
    ]

  coefficients%i(39) = 1 ! do not build the continuity equation
  coefficients%i(71) = 1  ! CDT for b-formulation
  coefficients%i(80) = alam_model ! adapted lambda
  coefficients%i(83) = physqb ! physical quantity of first mode
  coefficients%i(84) = 3  ! storage of b tensor: all modes in sysvector
  coefficients%i(86) = 2  ! sysvector number for velocity in SUPG
  coefficients%i(87) = 1  ! sysvector number for iteration
  coefficients%i(89) = 1  ! exclude time-derivative in implicit_ce_supg_elem
  coefficients%i(90) = 1  ! rotation reinitialization on element level

  if ( JacobianSUPG ) coefficients%i(91) = 1  ! Jacobian of SUPG

  coefficients%i(301) = 4 ! second order time integration
  coefficients%i(305) = 1 ! sysvector number for iteration velocity/pressure
  coefficients%i(306) = 1 ! exclude time-derivative in implicit_pressure_ce_elem
  if ( pressure_constraint ) then
    coefficients%i(307) = 1 ! constraint number of the pressure constraint
  end if
  coefficients%i(308) = 1 ! type of pressure constraint (1,2,3 or 4)

  coefficients%r(1:500) = &
    [ eta_s, 0._dp,   0._dp,     alpha, 0._dp, &
      0._dp, 0._dp, deltat3, beta_supg, 0._dp, &
      ( 0._dp, i = 11, 500 ) ]

  coefficients%r(29) = Kmod
  coefficients%r(30) = p0

  if ( alam_model == 4 ) then
    coefficients%r(501:503) = [ ( G/nmodes, lambda, tau_ref, m=1,nmodes ) ]
  else if ( alam_model == 5 ) then
    coefficients%r(501:505) = &
              [ ( G/nmodes, lambda, tau_ref1, f1, tau_ref2, m=1,nmodes ) ]
  end if

  coefficients%r(252) = Kmod
  coefficients%r(253) = p0

! read mesh

  call read_mesh ( mesh, filename='mesh.out' )

  if ( pressure_constraint ) then

!   add elementset for integral pressure constraint

    call add_to_mesh ( mesh, elementset='elements', &
                                 elements=[(i,i=1,mesh%nelem)] )
    call add_to_mesh ( mesh, elementset='nodes', elementsetnr=1 )

  end if

  call fill_mesh_parts ( mesh )

  if ( printscreen ) call printinfo ( mesh, printlevel=1 )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=6+nmodes, &
    nphysq=3+nmodes )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 4  ! G
  input_probdef%vec_elementdof(1)%a(:,2) = 2         ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(vertices,4:3+nmodes) = ncompb  ! b
  input_probdef%vec_elementdof(1)%a(:,4+nmodes) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5+nmodes) = 4  ! tensor for plotting
  input_probdef%vec_elementdof(1)%a(:,6+nmodes) = 5  ! tensor for plotting

  input_probdef%physq = [1,2,3,(3+m,m=1,nmodes)]
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
    physqarr=[physqgrad,physqvel,physqpress,(physqb+m-1,m=1,nmodes)] )

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
  call create_subscript ( mesh, problem, cxx, degfd=1, vec=5+nmodes )
  call create_subscript ( mesh, problem, cxy, degfd=2, vec=5+nmodes )
  call create_subscript ( mesh, problem, cyy, degfd=3, vec=5+nmodes )
  call create_subscript ( mesh, problem, czz, degfd=4, vec=5+nmodes )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, solsupg, soliter )
  call create_sysvector ( problem, dsol, rhsd, reacf )


! create a vector for conformation
! tensor for post processing

  call create_vector ( problem, ctensor, vec=5+nmodes )
  call create_vector ( problem, btensor, vec=6+nmodes )

! fill solution vector with essential boundary conditions

  soliter%u = 0

! initial pressure

  soliter%u(pres%s) = p0

! initial solution b tensor

  do m = 1, nmodes
    soliter%u(bvalcmp(1,m)%s) = 1 ! initial bxx
    soliter%u(bvalcmp(2,m)%s) = 0 ! initial bxy
    soliter%u(bvalcmp(3,m)%s) = 0 ! initial byx
    soliter%u(bvalcmp(4,m)%s) = 1 ! initial byy
    soliter%u(bvalcmp(5,m)%s) = 1 ! initial bzz
  end do

  call copy(soliter,solsupg)

  dsol%u = 0

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

  if ( printscreen ) print *, 'nnz = ', sysmatrix%Suu%nnz


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nprob=1 )


! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => soliter
  oldvectors_ve%s(2)%p => solsupg
  oldvectors_ve%p(1)%p => problem  ! needed for implicit_ce_supg_elem


  open ( unit=13, file='iter.out', recl=300 )
  open ( unit=14, file='Kdrag.out', recl=300 )

  solver_options%real_storage=rs
  solver_options%integer_storage=is

! area of the (original) domain

  call integrate ( mesh, problem, V0, elemsub=stokes_integrate_volume, &
    coefficients=coefficients )

  coefficients%r(254) = V0(1)

  call tic


! stepping through specified U values and continue from previous solution

  do step = 1, numsteps

!   average velocity in this step
    omega = Uspec(step)

!   Weisenberg number
    Wi = lambda * omega / (H-R)

!   increase in omega
    domega = Uspec(step) - Uspec(step-1)

  ! set velocity_theta = domega*R on cylinder
    call fill_sysvector ( mesh, problem, dsol, &
      curve1=13, physq=physqvel, vfunc=vfunc, vfuncnr=1 )

    coefficients%r(10) = omega*R ! global scaling velocity for SUPG, when needed

    write(13,*) 'step = ', step, ' omega = ', omega, 'Wi = ', Wi
    if ( printscreen) print *, 'step = ', step, ' omega = ', omega, 'Wi = ', Wi

    if ( step == 1 ) then
      coefficients%i(78) = 1 ! check for zero velocity in SUPG
    else
      coefficients%i(78) = 0
    end if

!   Newton-Raphson convergence criterium
    epsconf = epsc(step)

    iter = 0

    do

      iter = iter + 1

      if ( JacobianSUPG ) then
!       copy iterative solution velocity for SUPG
        call copy ( soliter, solsupg )
      else if ( iter == 2 ) then
!       copy iterative solution velocity for SUPG for correct velocity bc
!       Also solsupg is constant during the iteration process.
        call copy ( soliter, solsupg )
      end if

      if ( iter == 2 ) then
!       set dsol = 0 on cylinder
        call fill_sysvector ( mesh, problem, dsol, &
          curve1=13, physq=physqvel, value=0._dp )
      end if

      if ( iter > maxnumiterations ) then
        write(*,'(3(a,i0/))') &
          ' Maximum number of iterations reached = ', &
          maxnumiterations, ' step = ', step
        stop
      end if

!     build (assemble) matrix/vector for gradient/velocity/pressure part

      call build_vpG

      call toc ( 'build_vpG' )

!     build (assemble) matrix/vector for conformation part

      do m = 1, nmodes

        coefficients%i(85) = m   ! set mode number

        call build_b

      end do

      call toc ( 'build_b' )

      if ( .not. physqmask .or. &
              physqmask .and. .not. physqmask_offdiag ) then

!       set to zero off-diagonal conformation blocks

        call build_offdiag_b

        call toc ( 'build_offdiag_b' )

      end if

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, dsol, rhsd )

      call check ( sysmatrix )

!     solve

      call solve_system_ma41 ( sysmatrix, rhsd, dsol, &
        solver_options=solver_options  )

      call toc ( 'solve' )

      soliter%u(solsc%s) = soliter%u(solsc%s) + dsol%u(solsc%s)

!     convergence test

      epsu = maxval(abs(dsol%u(vel%s))) / ( omega * R )
      epsp = maxval(abs(dsol%u(pres%s))) * R / ( ( omega * R ) * eta_0 )
      epsb = maxval(abs(dsol%u(bval(1)%s)))

      if ( printscreen) print *, iter, epsu, epsp, epsb
      write(13,*) iter, epsu, epsp, epsb

      if ( maxval([ epsu, epsp, epsb ]) < epsconf ) exit

    end do

    call toc ( 'after iteration loop' )


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
    write(14,*) step, Wi, K_dragx, K_dragy, K_torqz

    call toc ( 'after reaction forces' )


!   copy solution to velocity for supg for next step

    call copy ( soliter, solsupg )

!   write max and mean values of conformation and b-tensor tensor to a file

    call derive_vector ( mesh, problem, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    if ( step == 1 ) then
      open(unit=11, recl=600, status='replace', file='cval.out')
    else
      open(unit=11, recl=600, position='append', file='cval.out')
    end if
    write(11, fmt=*) step, maxval(ctensor%u(cxx%s)), &
                           maxval(ctensor%u(cxy%s)), &
                           maxval(ctensor%u(cyy%s)), &
                           maxval(ctensor%u(czz%s)), &
                           sum(ctensor%u(cxx%s))/size(cxx%s), &
                           sum(ctensor%u(cxy%s))/size(cxy%s), &
                           sum(ctensor%u(cyy%s))/size(cyy%s), &
                           sum(ctensor%u(czz%s))/size(czz%s), &
                           maxval(soliter%u(velx%s)), &
                           maxval(soliter%u(vely%s)), &
                           maxval(soliter%u(pres%s)), &
                           minval(soliter%u(pres%s))

    close(unit=11)

    call toc ( 'one step' )

  end do

  call toc ( 'all steps' )


! post-processing

  call create_vector ( problem, velocity, physq=2 )
  call create_vector ( problem, pressure, vec=4+nmodes )
  call create_vector ( problem, Jvol, vec=4+nmodes )
  call create_vector ( problem, vorticity, vec=4+nmodes )
  call create_vector ( problem, divu, vec=4+nmodes )

  call extract_physvector ( mesh, problem, soliter, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! derive vectors

  oldvectors%s(1)%p => soliter

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

! compute J directly from the pressure
  Jvol%u = exp(-(pressure%u-p0)/Kmod)

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=7

  call derive_vector ( mesh, problem, divu, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='divu', vector=divu, append=.true. )
  call write_vector_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='velocity', vector=velocity )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='pressure', vector=pressure, append=.true. )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='J', vector=Jvol, append=.true. )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='vorticity', vector=vorticity, append=.true. )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='divu', vector=divu, append=.true. )

  call printtofile ( mesh, problem, filename='velocity_cl.out', curve=11, &
    vector=velocity )
  call printtofile ( mesh, problem, filename='pressure_cl.out', curve=12, &
    vector=pressure )

! conformation tensor

  call derive_vector ( mesh, problem, ctensor, &
    elemsub=deriv_conformation_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problem, filename='c.vtk', &
    dataname='conformation_tensor', vector=ctensor, assume33=.true. )

  call delete ( ctensor )

! b tensor

  call derive_vector ( mesh, problem, btensor, &
    elemsub=deriv_conformation_tensor_std, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problem, filename='b.vtk', &
    dataname='b-tensor', vector=btensor, symmetric=.false., assume33=.true. )

  call delete ( btensor )

! von Mises stress

  call create_vector ( problem, vonmises, vec=4+nmodes )

  coefficients%i(13)=1 ! von Mises
  coefficients%i(28)=1 ! mode number

  call derive_vector ( mesh, problem, vonmises, &
    elemsub=deriv_viscoelastic_stress_scalar, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_scalar_vtk ( mesh, problem, filename='c.vtk', &
    dataname='von_mises', vector=vonmises, append=.true. )

  call delete ( vonmises )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( solsupg, soliter, dsol, rhsd, reacf )
  call delete ( velocity, pressure, vorticity, Jvol, divu )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

  call delete ( oldvectors_ve )
  do m = 1, nmodes
    do i = 1, ncompb
      call delete ( bvalcmp(i,m), bval(m) )
    end do
  end do
  call delete ( cxx, cxy, cyy, czz )

contains


  function vfunc ( n, nr, x )
    use kind_defs_m
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    vfunc = domega * [ -x(2), x(1) ] / norm2(x)

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

    if ( nmodes > 1 ) then
      do i = 1, nmodes
        do j = 1, nmodes
          if ( i == j ) cycle
            input_probdef%physqmask(physqb+i-1,physqb+j-1) = .false.
        end do
      end do
    end if

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

!   build (assemble) matrix and vector for b-tensor mode m

    if ( JacobianSUPG ) then

!     add diagonal block of steady-state timed derivative
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=steady_ce_timederiv_supg_elem, coefficients=coefficients, &
        physqrow=[physqb+m-1], physqcol=[physqvel,physqb+m-1], &
        oldvectors=oldvectors_ve, addmatvec=.true. )

!     pressure, velocity and diagonal block
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
        physqrow=[physqb+m-1], physqcol=[physqpress,physqvel,physqb+m-1], &
        oldvectors=oldvectors_ve, addmatvec=.true. )

    else

!      add diagonal block of steady-state timed derivative
       call build_system ( mesh, problem, sysmatrix, rhsd, &
         elemsub=steady_ce_timederiv_supg_elem, coefficients=coefficients, &
         physqrow=[physqb+m-1], physqcol=[physqb+m-1], &
         oldvectors=oldvectors_ve, addmatvec=.true. )

!      pressure and diagonal block
       call build_system ( mesh, problem, sysmatrix, rhsd, &
         elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
         physqrow=[physqb+m-1], physqcol=[physqpress,physqb+m-1], &
         oldvectors=oldvectors_ve, addmatvec=.true. )

    end if

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_ce_vel_supg_elem, coefficients=coefficients, &
      physqrow=[physqb+m-1], physqcol=[physqgrad,physqvel], &
      oldvectors=oldvectors_ve, addmatvec=.true., buildvector=.false. )

    if ( .not. physqmask .or. &
            physqmask .and. .not. physqmask_grad ) then
!     set to zero off-diagonal blocks gradient-b
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        physqrow=[physqgrad], physqcol=[physqb+m-1], &
        buildvector=.false., zeromatvec=.true. )
    end if

    if ( .not. physqmask .or. &
            physqmask .and. .not. physqmask_press ) then
!     set to zero off-diagonal block pressure-b
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        physqrow=[physqpress], physqcol=[physqb+m-1], &
        buildvector=.false., zeromatvec=.true. )
    end if

  end subroutine build_b

  subroutine build_offdiag_b

    integer :: i, j

    if ( nmodes > 1 ) then
      do i = 1, nmodes
        do j = 1, nmodes
          if ( i == j ) cycle
          call build_system ( mesh, problem, sysmatrix, rhsd, &
            addmatvec=.true., physqrow=[physqb+i-1], physqcol=[physqb+j-1], &
            buildvector=.false., zeromatvec=.true. )
        end do
      end do
    end if

  end subroutine build_offdiag_b

end program cylinder_egp6
