#if METIS5

! Sequence of steady states of flow around a cylinder for an Oldroyd-B or
! Giesekus model using the conformation tensor formulation.
! Iteration to steady-state with Newton-Raphson.
! Periodical boundary conditions.
! Compute drag on cylinder using reaction forces.
! Optionally include Jacobian of SUPG test function.
! Optionally discard certain zero blocks from the sparse matrix.
! Optionally restart from a previously saved solution.
! Optionally scale previous solution to match newly imposed flow rate.

program cylinder_c10

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use io_utils_m
  use figplot_m
  use timer_m
  use metis5_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,     & ! P2 velocities
    pintpl = 2,     & ! P1 pressures
    gintpl = 2,     & ! P1 gradients
    cintpl = 2,     & ! P1 conformation
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    physqc = 4,     & ! physical quantity nr of the conformation tensor (mode 1)
    gauss = 6,      & ! 6 point integration of triangles
    gaussb = 3,     & ! 3 point integration of boundary elements
    logc = 0,       & ! standard scheme or log transformation
    ncompc = 3,     & ! number of conformation tensor comp
    nmodes = 1,     & ! number of modes
    startm = 501      ! start of material model data

  real(dp), parameter :: &
    H = 2._dp,           & ! (half-) height of the channel (must match mesh)
    R = 1._dp            ! radius of the cylinder (must match mesh)

  real(dp), parameter :: &
    beta_supg = 1          ! SUPG factor

  integer, parameter :: &
    ndU = 1           ! number intervals with different dU increment

  integer :: &
    timeint3 = 10,     & ! steady time discretization, BDF1: 8, BDF2 10
    maxnumnriter = 30, & ! maximum number of Newton-Raphson iterations
    restart = 0,       & ! restart=0: no restart
                         ! restart=1: normal restart
                         ! restart=2: read restart file + post processing
    restart_every = 1, & ! write restart file every restart_every steps.
                         ! -1: means none
    model = 3            ! model number 2: Oldroyd-B 3: Giesekus

  real(dp) :: &
    G      = 1.0_dp,    & ! modulus
    lambda = 1.0_dp,    & ! relaxation time
    beta_s = 0.0_dp,    & ! beta viscosity parameter = eta_s/(eta_s+eta_p)
    mobility = 0.01_dp, & ! mobility parameter for the Giesekus model
    deltat3 = 0.001_dp, & ! time step for steady iteration scheme
    rs = 1.2_dp,        & ! real_storage for LU (HSL)
    is = 1.6_dp,        & ! integer_storage for LU (HSL)
    dU(ndU),            & ! increments in imposed average velocity
    Ue(ndU),            & ! end values of the intervals in U
    ed(ndU),            & ! Newton-Raphson convergence threshold (difference)
    er(ndU)               ! Newton-Raphson convergence threshold (residual)

  logical :: &
    physqmask = .true., & ! use physical quantity masking
    physqmask_grad = .true., & ! use physmask for G-b coupling
    physqmask_press = .false., & ! use physmask for p-b and b-p coupling
    physqmask_offdiag = .true., & ! use physmask for b-b off-diagonal coupling
    JacobianSUPG = .true., & ! take Jacobian of SUPG test function into account
    scaleU = .false., & ! scale velocity/gradients at the start of a step
    metisnodes = .true., & ! renumber nodes with metis5
    metisdegrees = .true., & ! renumber degrees of freedom with metis5
    postprocess = .false., & ! perform postprocessing after the stepping
    printscreen = .true.  ! print on standard output

  character(len=30) :: meshfile = '', inrestartfile ='', outrestartfile =''

! definitions

  type(mesh_t) :: mesh, mesh1
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: solsupg, soliter
  type(sysvector_t) :: dsol, rhsd, reacf
  type(vector_t) :: cten
  type(subscript_t) :: velx, vely, vel, pres, solsc, velxcylinder, grad
  type(subscript_t) :: res_vel
  type(subscript_t), dimension(nmodes) :: cvalxx, cvalxy, cvalyy, cval, res_cval
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options

  type(subscriptvec_t) :: cxx, cxy, cyy

  integer :: step, iter, m, i, numsteps, ncoefr
  integer :: vertices(3) = [1,3,5]
  real(dp) :: alpha, dflowrate, epsdconf, epsrconf, U, K_drag, U_prev
  real(dp) :: epsu, epsp, epsc, Wi, Ru, Rb
  real(dp) :: eta_p, eta_s, eta_0
  real(dp), dimension(:), allocatable :: Uspec, epsd, epsr

! namelist for input of variables; read from standard input

  namelist /comppar/ timer, meshfile, inrestartfile, outrestartfile, &
    timeint3, restart, restart_every, model, G, lambda, beta_s, mobility, &
    deltat3, rs, is, dU, Ue, ed, er, scaleU, postprocess, &
    maxnumnriter, metisnodes, metisdegrees

  read ( unit=*, nml=comppar )

  if ( restart == 2 ) postprocess = .true.

! set some parameters

  eta_p = lambda*G                ! polymer viscosity
  eta_s = beta_s/(1-beta_s)*eta_p ! solvent viscosity
  eta_0 = eta_s+eta_p             ! zero-shear viscosity
  alpha = eta_p                   ! DEVSS parameter

  if ( printscreen ) then
    print *
    print '(1x,a,i0)', 'model      = ', model
    print *, 'G          =', G
    if ( model == 3 ) print *, 'mobility =', mobility
    print *, 'lambda     =', lambda
    print *, 'beta_s     =', beta_s
    print *, 'eta_p      =', eta_p
    print *, 'eta_s      =', eta_s
    print *, 'eta_0      =', eta_0
  end if

! fill coefficients

  if ( model==2 ) then
    ncoefr = 500 + 2*nmodes
  else if ( model==3 ) then
    ncoefr = 500 + 3*nmodes
  end if

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=ncoefr )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gaussb,   cintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      logc,     timeint3,  ( 0, i = 23, 150 )  &
    ]

  coefficients%i(83) = physqc ! physical quantity of first mode
  coefficients%i(84) = 3  ! storage of c tensor: all modes in sysvector
  coefficients%i(86) = 2  ! sysvector number for velocity in SUPG
  coefficients%i(87) = 1  ! sysvector number for iteration
  coefficients%i(89) = 1  ! exclude time-derivative in implicit_ce_supg_elem

  if ( JacobianSUPG ) coefficients%i(91) = 1  ! Jacobian of SUPG

  coefficients%r(:500) = &
     [ eta_s, 0._dp,   0._dp,     alpha, 0._dp, &
       0._dp, 0._dp, deltat3, beta_supg, 0._dp, &
       ( 0._dp, i = 11, 500 ) ]

  if ( model==2 ) then
    coefficients%r(501:) = [ ( G/nmodes, lambda, m=1,nmodes ) ]
  else if ( model==3 ) then
    coefficients%r(501:) = [ ( G/nmodes, lambda, mobility, m=1,nmodes ) ]
  end if

! read mesh

  call read_mesh_gmsh ( mesh1, filename=meshfile, ndim=2 )
  call mesh_convert ( mesh1, mesh, remove_isolated_nodes=.true. )
  call delete ( mesh1 )

! cylinder
  call add_to_mesh ( mesh, curve=[-2] ) ! curve 9
! top
  call add_to_mesh ( mesh, curve=[5,6,7] ) ! curve 10
! centerline+cylinder
  call add_to_mesh ( mesh, curve=[1,-2,3] ) ! curve 11
! connect curves for periodical bc in DG
  call add_to_mesh ( mesh, curve=[-8] )     ! curve 12

  call fill_mesh_parts ( mesh )

  call tic

  if ( metisnodes ) then

    call renumber_metis ( mesh )

    call toc ('renumber_metis')

  end if

  if ( printscreen ) call printinfo ( mesh, printlevel=1 )

! plot mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  plot_options%fontsize = 10

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=5+nmodes, &
     nphysq=3+nmodes )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 4  ! G
  input_probdef%vec_elementdof(1)%a(:,2) = 2         ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(vertices,4:3+nmodes) = ncompc  ! c
  input_probdef%vec_elementdof(1)%a(:,4+nmodes) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5+nmodes) = 3  ! tensor for plotting

  input_probdef%physq = [1,2,3,(3+m,m=1,nmodes)]
  input_probdef%probnr = 1

  if ( physqmask ) then

!   discard certain zero blocks from the sparse matrix

    if ( physqmask_grad ) call fill_grad_physqmask         ! G-c block
    if ( physqmask_press ) call fill_press_physqmask       ! p-c and c-p blocks
    if ( physqmask_offdiag ) call fill_offdiag_physqmask   ! off-diagonal c

  end if

! define essential boundaries

! cylinder
  call define_essential ( mesh, input_probdef, curve1=9, physq=physqvel )
! top (wall)
  call define_essential ( mesh, input_probdef, curve1=10, physq=physqvel )
! bottom (center line)
  call define_essential ( mesh, input_probdef, curve1=1, &
    physq=physqvel, degfd=[0,1] )
  call define_essential ( mesh, input_probdef, curve1=3, &
    physq=physqvel, degfd=[0,1] )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=12, nglobalc=1 )

! constraints for periodical boundary conditions

! velocities (use weak connection)
  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=4, curve2=12, discretization='weak', &
    elementdof=[2,0,2] )

! gradients (use collocation)
  call define_constraint ( mesh, input_probdef, &
    physq=physqgrad, curve1=4, curve2=12, discretization='collocation' )

! conformation tensor (use collocation)
  do m = 1, nmodes
    call define_constraint ( mesh, input_probdef, &
      physq=physqc+m-1, curve1=4, curve2=12, discretization='collocation' )
  end do

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for solution (excluding contraint forces)

  call create_subscript ( mesh, problem, solsc, &
    physqarr=[physqgrad,physqvel,physqpress,(physqc+m-1,m=1,nmodes)] )

! create vector subscripts for the velocity

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )
! velocity degrees in x-direction on the cylinder
  call create_subscript ( mesh, problem, velxcylinder, physqarr=[physqvel], &
    degfd=1, curves=[9] )
! excluding essential and constraint degrees for computing the residual Ru
  call create_subscript ( mesh, problem, res_vel, physqarr=[physqvel], &
    essentialpart=.false., excludecurves=[4,12] )

! create vector subscript for the gradients

  call create_subscript ( mesh, problem, grad, physqarr=[physqgrad] )

! create vector subscript for the pressure

  call create_subscript ( mesh, problem, pres, physqarr=[physqpress] )

! create vector subscript for the conformation tensor

  do m = 1, nmodes
    call create_subscript ( mesh, problem, cval(m), physqarr=[physqc+m-1] )
!   excluding essential and constraint degrees for computing the residual Rb
    call create_subscript ( mesh, problem, res_cval(m), &
      physqarr=[physqc+m-1], essentialpart=.false., excludecurves=[4,12] )
    call create_subscript ( mesh, problem, cvalxx(m), physqarr=[physqc+m-1], &
      degfd=1 )
    call create_subscript ( mesh, problem, cvalxy(m), physqarr=[physqc+m-1], &
      degfd=2 )
    call create_subscript ( mesh, problem, cvalyy(m), physqarr=[physqc+m-1], &
      degfd=3 )
  end do
  call create_subscript ( mesh, problem, cxx, degfd=1, vec=5+nmodes )
  call create_subscript ( mesh, problem, cxy, degfd=2, vec=5+nmodes )
  call create_subscript ( mesh, problem, cyy, degfd=3, vec=5+nmodes )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, solsupg, soliter )
  call create_sysvector ( problem, dsol, rhsd, reacf )

! create vector cten for output of max min c during time stepping

  call create_vector ( problem, cten, vec=5+nmodes )

! fill solution vector with essential boundary conditions

  soliter%u = 0

! initial solution conformation tensor

  do m = 1, nmodes
    if ( logc == 0 ) then ! standard
      soliter%u(cvalxx(m)%s) = 1 ! initial cxx
      soliter%u(cvalxy(m)%s) = 0 ! initial cxy
      soliter%u(cvalyy(m)%s) = 1 ! initial cyy
    else if ( logc == 1 ) then ! log scheme
      soliter%u(cval(m)%s) = 0 ! initial s
    end if
  end do

  call copy(soliter,solsupg)
  call copy(soliter,dsol)

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    usephysqmask=physqmask )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

  if ( printscreen ) then
!   print size parameters
    print *, 'n   = ', sysmatrix%Suu%n
    print *, 'nu  = ', size(vel%s)
    print *, 'ng  = ', size(grad%s)
    print *, 'np  = ', size(pres%s)
    print *, 'nb  = ', size(cval(1)%s)
    print *, 'nl  = ', soliter%n - size(solsc%s)
    print *, 'nnz = ', sysmatrix%Suu%nnz
    print *
  end if


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nprob=1 )


! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => soliter
  oldvectors_ve%s(2)%p => solsupg
  oldvectors_ve%p(1)%p => problem  ! needed for implicit_ce_supg_elem


  open ( unit=17, file='iter.out', recl=300 )
  open ( unit=14, file='Kdrag.out', recl=300 )

  solver_options%real_storage=rs
  solver_options%integer_storage=is

  if ( metisdegrees ) then

!   use pivot order from the permutation arrays

    solver_options%pivot_order = 1

  end if

! restart

  if ( any ( restart == [1,2] ) ) then
    call input_restart
  else
    U_prev = 0 ! set previous U to zero if not restart
  end if

! set incremental steps

  if ( restart == 2 ) then
    numsteps = 0
  else
    call set_steps
  end if

! stepping through specified U values and continue from previous solution

  call toc ( 'start of stepping' )

  do step = 1, numsteps

!   average velocity in this step
    U = Uspec(step)

!   Weisenberg number
    Wi = lambda * U / R

!   increase in flow rate in (half) the channel
    dflowrate = H*(Uspec(step) - Uspec(step-1))

    if ( scaleU ) then
!     scale velocities/gradients to new flowrate
      if ( step == 1 .and. restart == 0 ) then
!       impose flowrate
        coefficients%r(6) = dflowrate
      else
!       scale and set delta(flowrate) to zero
        soliter%u(vel%s) = Uspec(step) / Uspec(step-1) * soliter%u(vel%s)
        soliter%u(grad%s) = Uspec(step) / Uspec(step-1) * soliter%u(grad%s)
        coefficients%r(6) = 0
      end if
    else
      coefficients%r(6) = dflowrate
    end if

    coefficients%r(10) = U ! global scaling velocity for SUPG, when needed

    write(17,*) 'step = ', step, ' U = ', U, 'Wi = ', Wi
    if ( printscreen) print *, 'step = ', step, ' U = ', U, 'Wi = ', Wi

    if ( step == 1 .and. restart == 0 ) then
      coefficients%i(78) = 1 ! check for zero velocity in SUPG
    else
      coefficients%i(78) = 0
    end if

!   Newton-Raphson convergence criteria for difference and residual
    epsdconf = epsd(step)
    epsrconf = epsr(step)

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

      if ( iter == 2 ) coefficients%r(6) = 0  ! set flowrate in dsol to 0

      if ( iter > maxnumnriter ) then
        write(*,'(3(a,i0/))') &
          ' Maximum number of iterations reached = ', &
          maxnumnriter, ' step = ', step
        stop
      end if

!     build (assemble) matrix/vector for gradient/velocity/pressure part

      call build_vpG

      call toc ( 'build_vpG' )

!     build (assemble) matrix/vector for conformation part

      do m = 1, nmodes

        coefficients%i(85) = m   ! set mode number

        call build_c

      end do

      call toc ( 'build_c' )

      if ( .not. physqmask .or. &
              physqmask .and. .not. physqmask_offdiag ) then

!       set to zero off-diagonal conformation blocks

        call build_offdiag_c

        call toc ( 'build_offdiag_c' )

      end if

      if ( metisdegrees .and. step == 1 .and. iter == 1 ) then

!       generate permutation arrays with metis

        call renumber_metis ( sysmatrix )

        call toc ( 'renumber metis sysmatrix' )

      end if

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, dsol, rhsd )

      call check ( sysmatrix )

!     compute residual norms

      Ru = norm2(rhsd%u(res_vel%s)) &
                 / sqrt(real(size(res_vel%s),dp)) * R**2 / ( U * eta_0 )
      Rb = norm2(rhsd%u(res_cval(1)%s)) &
                / sqrt(real(size(res_cval(1)%s),dp)) * R / U

!     solve

      call solve_system_ma41 ( sysmatrix, rhsd, dsol, &
        solver_options=solver_options  )

      call toc ( 'solve' )

      soliter%u(solsc%s) = soliter%u(solsc%s) + dsol%u(solsc%s)

!     convergence test

      epsu = maxval(abs(dsol%u(vel%s))) / U
      epsp = maxval(abs(dsol%u(pres%s))) * R / ( U * eta_0 )
      epsc = maxval(abs(dsol%u(cval(1)%s)))

      if ( printscreen) print *, iter, epsu, epsp, epsc, Ru, Rb
      write(17,*) iter, epsu, epsp, epsc, Ru, Rb

      if ( maxval([ epsu, epsp, epsc ]) < epsdconf .and. &
           maxval([ Ru, Rb ]) < epsrconf ) exit

    end do

    call toc ( 'after iteration loop' )


!   reaction forces

    call reaction_forces ( problem, sysmatrix, dsol, rhsd, reacf )

    K_drag = - 2 * sum ( reacf%u(velxcylinder%s) ) / ( eta_0*U )

    if ( printscreen) print *, 'Kdrag = ', K_drag
    write(14,*) step, Wi, K_drag

    call toc ( 'after reaction forces' )


!   copy solution to velocity for supg for next step

    call copy ( soliter, solsupg )

!   write max and mean values of conformation tensor to a file

    call derive_vector ( mesh, problem, cten, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    if ( step == 1 ) then
      open(unit=11, recl=300, status='replace', file='cval.out')
    else
      open(unit=11, recl=300, position='append', file='cval.out')
    end if
    write(11, fmt=*) step, maxval(cten%u(cxx%s)), &
                           maxval(cten%u(cxy%s)), &
                           maxval(cten%u(cyy%s)), &
                           sum(cten%u(cxx%s))/size(cxx%s), &
                           sum(cten%u(cxy%s))/size(cxy%s), &
                           sum(cten%u(cyy%s))/size(cyy%s), &
                           maxval(soliter%u(velx%s)), &
                           maxval(soliter%u(vely%s))
    close(unit=11)

    call toc ( 'one step' )

!   write data for restart

    if ( restart_every > 0 ) then
      if ( mod(step,restart_every) == 0 ) then

        call output_restart

      end if
    end if

  end do

  call toc ( 'all steps' )


! postprocessing

  if ( postprocess ) call postprocessing


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( solsupg, soliter, dsol, rhsd, reacf )
  call delete ( cten )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( vel, grad, pres, res_vel )

  call delete ( oldvectors_ve )
  do m = 1, nmodes
    call delete ( cvalxx(m), cvalxy(m), cvalyy(m), cval(m), res_cval(m) )
  end do
  call delete ( solsc )
  call delete ( cxx, cxy, cyy )

contains


  subroutine fill_grad_physqmask

    integer :: j

    do j = 1, nmodes
      input_probdef%physqmask(physqgrad,physqc+j-1) = .false.
    end do

  end subroutine fill_grad_physqmask


  subroutine fill_press_physqmask

    integer :: j

    do j = 1, nmodes
      input_probdef%physqmask(physqpress,physqc+j-1) = .false.
      input_probdef%physqmask(physqc+j-1,physqpress) = .false.
    end do

  end subroutine fill_press_physqmask


  subroutine fill_offdiag_physqmask

    integer :: i, j

    if ( nmodes > 1 ) then
      do i = 1, nmodes
        do j = 1, nmodes
          if ( i == j ) cycle
            input_probdef%physqmask(physqc+i-1,physqc+j-1) = .false.
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

!   stokes add (q,divu) to rhs of mass balance
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_rhs_divu, oldvectors=oldvectors_ve, &
      coefficients=coefficients, addmatvec=.true., &
      physqrow=[physqpress], physqcol=[physqvel], buildmatrix=.false. )

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
      physqrow=[physqvel], physqcol=[(physqc+m-1,m=1,nmodes)], &
      addmatvec=.true., coefficients=coefficients )

!   flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients )

!   periodical condition on velocities

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=stokes_constr_elem_conn, addmatvec=.true., &
      coefficients=coefficients )

!   periodical condition on gradients

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=3, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients )

  end subroutine build_vpG


  subroutine build_c

!   build (assemble) matrix and vector for conformation tensor mode m

    if ( JacobianSUPG ) then

!     diagonal block
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
        physqrow=[physqc+m-1], physqcol=[physqvel,physqc+m-1], &
        oldvectors=oldvectors_ve, addmatvec=.true. )

    else

!     diagonal block
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
        physqrow=[physqc+m-1], physqcol=[physqc+m-1], &
        oldvectors=oldvectors_ve, addmatvec=.true. )

    end if

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_ce_vel_supg_elem, coefficients=coefficients, &
      physqrow=[physqc+m-1], physqcol=[physqgrad,physqvel], &
      oldvectors=oldvectors_ve, addmatvec=.true., buildvector=.false. )

    if ( .not. physqmask .or. &
            physqmask .and. .not. physqmask_grad ) then
!     set to zero off-diagonal blocks gradient-b
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        physqrow=[physqgrad], physqcol=[physqc+m-1], &
        buildvector=.false., zeromatvec=.true. )
    end if

    if ( .not. physqmask .or. &
            physqmask .and. .not. physqmask_press ) then
!     set to zero off-diagonal blocks pressure-b
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        physqrow=[physqpress], physqcol=[physqc+m-1], &
        buildvector=.false., zeromatvec=.true. )
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        physqrow=[physqc+m-1], physqcol=[physqpress], &
        buildvector=.false., zeromatvec=.true. )
    end if

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problem, sysmatrix, &
      sysvector=rhsd, constraint1=3+m, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

  end subroutine build_c

  subroutine build_offdiag_c

    integer :: i, j

    if ( nmodes > 1 ) then
      do i = 1, nmodes
        do j = 1, nmodes
          if ( i == j ) cycle
          call build_system ( mesh, problem, sysmatrix, rhsd, &
            addmatvec=.true., physqrow=[physqc+i-1], physqcol=[physqc+j-1], &
            buildvector=.false., zeromatvec=.true. )
        end do
      end do
    end if

  end subroutine build_offdiag_c

! output restart data

  subroutine output_restart

    write(17,*) 'Output restart data at step ', step, ' U = ', U, 'Wi = ', Wi
    if ( printscreen) &
      print *, 'Output restart data at step ', step, ' U = ', U, 'Wi = ', Wi

    open ( unit=15, form='unformatted', file=outrestartfile )

    write(15) U, 0._dp, soliter%u

    close(unit=15)

  end subroutine output_restart


! open files and input restart data

  subroutine input_restart

    real(dp) :: dummy

    open ( unit=15, form='unformatted', file=inrestartfile )

    read(15) U_prev, dummy, soliter%u

    close(unit=15)

    write(17,*) 'Input restart data, U = ', U_prev
    if ( printscreen) &
      print *, 'Input restart data, U = ', U_prev

  end subroutine input_restart


! set incremental steps in U

  subroutine set_steps

    integer :: i, j, nsteps(ndU)

!   number of incremental steps in U to perform

    if ( restart == 1 .and. abs(dU(1)) <= tiny(1._dp) ) then
      nsteps(1) = 1  ! first interval after restart and no increase in velocity
    else
      nsteps(1) = nint( (Ue(1)-U_prev)/dU(1) ) ! first interval
    end if

    do i = 2, ndU
      nsteps(i) = nint( (Ue(i)-Ue(i-1))/dU(i) ) ! second and higher intervals
    end do
    numsteps = sum(nsteps)

    if ( printscreen ) print *, 'numsteps = ', numsteps

!   U values specified

    allocate( Uspec(0:numsteps), epsd(numsteps), epsr(numsteps) )

    Uspec = [ U_prev, ( U_prev + i*dU(1), i=1,nsteps(1) ), &
                      ( ( Ue(j-1) + i*dU(j), i=1,nsteps(j)), j=2,ndU )  ]

    if ( printscreen ) print *, 'Uspec = ', Uspec(1:numsteps); print *

    epsd = [ ( ( ed(j), i=1,nsteps(j) ), j=1,ndU ) ]
    epsr = [ ( ( er(j), i=1,nsteps(j) ), j=1,ndU ) ]

  end subroutine set_steps


! post-processing

  subroutine postprocessing

    type(vector_t) :: velocity, pressure, vorticity, ctensor

    call create_vector ( problem, velocity, physq=2 )
    call create_vector ( problem, pressure, vec=4+nmodes )
    call create_vector ( problem, vorticity, vec=4+nmodes )

    call extract_physvector ( mesh, problem, soliter, velocity )

!   create the structure oldvectors

    call create_oldvectors ( oldvectors, nsysvec=1 )

!   derive vectors

    oldvectors%s(1)%p => soliter

    call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
      coefficients=coefficients, oldvectors=oldvectors )

    coefficients%i(13)=5

    call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
      coefficients=coefficients, oldvectors=oldvectors )

    call write_vector_vtk ( mesh, problem, filename='sol.vtk', &
      dataname='velocity', vector=velocity )
    call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
      dataname='pressure', vector=pressure, append=.true. )
    call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
      dataname='vorticity', vector=vorticity, append=.true. )

    call printtofile ( mesh, problem, filename='velocity_cl.out', curve=11, &
      vector=velocity )
    call printtofile ( mesh, problem, filename='pressure_cl.out', curve=11, &
      vector=pressure )

    call delete ( velocity, pressure, vorticity )
    call delete ( oldvectors )

!   conformation tensor

    call create_vector ( problem, ctensor, vec=5+nmodes )

    call derive_vector ( mesh, problem, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    call write_tensor_vtk ( mesh, problem, filename='c.vtk', &
      dataname='conformation_tensor', vector=ctensor )

    call printtofile ( mesh, problem, filename='conformation_cl.out', &
      curve=11, vector=ctensor )

    call delete ( ctensor )

  end subroutine postprocessing

end program cylinder_c10

#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on metis5', &
    ' - add metis5 lib for linking', &
    ' - set preprocessing macro METIS5 in Mdefs.mk'
end
#endif
