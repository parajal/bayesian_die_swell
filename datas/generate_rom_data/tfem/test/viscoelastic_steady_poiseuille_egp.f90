! Sequence of steady states of plane Poiseuille flow having a constant flow
! rate for an incompressible EGP model using the conformation tensor
! formulation.
! Iteration to steady-state with Newton-Raphson.
! Optionally include Jacobian of SUPG test function.
! Optionally discard certain zero blocks from the sparse matrix.

program poiseuille6

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use io_utils_m
!  use figplot_m
  use timer_m

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
    physqc = 4,     & ! physical quantity nr of the conformation tensor (mode 1)
    gauss = 3,      & ! 3x3 integration of quads
    gaussb = 3,     & ! 3 point integration of boundary elements
    timeint  = 8,   & ! (first-order) time integration
    maxnumiterations = 30, & ! maximum number of Newton-Raphson iterations
    logc = 0,       & ! standard scheme or log transformation
    ncompc = 4,     & ! number of conformation tensor comp
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 23        ! incompressible EGP

  real(dp), parameter :: &
    G     = 4.0_dp,      & ! modulus
    lambda = 10.0_dp,    & ! relaxation time
    eta_s = 1.0_dp,      & ! solvent viscosity
    eta_p = lambda*G,    & ! polymer viscosity
    eta_0 = eta_s+eta_p, & ! zero-shear viscosity
    alam_model = 4,      & ! adapted lambda model: 4: Eyring, 5: Ree-Eyring
    tau_ref = 1.0_dp,    & ! reference von Mises
    tau_ref1 = 1.0_dp,   & ! reference von Mises 1 RE
    tau_ref2 = 4.0_dp,   & ! reference von Mises 2 RE
    f1 = 0.6_dp            ! weight factor for first term RE

  real(dp), parameter :: &
    deltat = 0.2_dp,     & ! time step
    beta_supg = 1,       & ! SUPG factor
    rs    = 1.2_dp,      & ! real_storage for LU (HSL)
    is    = 1.6_dp         ! integer_storage for LU (HSL)

  integer, parameter :: &
    ndU = 2           ! number intervals with different dU increment

  real(dp), parameter :: &
!   increments in imposed average velocity
    dU(ndU) = [ 0.02_dp, 0.01_dp ], & ! increments in imposed average velocity
!   end values of the intervals in U
    Ue(ndU) = [ 0.04_dp, 0.05_dp ], & ! end values of the intervals in U
!   Newton-Raphson convergence threshold for each interval
    epsconfin(ndU) = [ 1e-10_dp, 1e-10_dp ]

  logical :: &
    physqmask = .true., & ! use physical quantity masking
    physqmask_grad = .true., & ! use physmask for G-b coupling
    physqmask_press = .false., & ! use physmask for p-b and b-p coupling
    physqmask_offdiag = .true., & ! use physmask for b-b off-diagonal coupling
    JacobianSUPG = .true., & ! take Jacobian of SUPG test function into account
    printscreen = .true.  ! print on standard output

! variables

  integer :: &
    nx=3,               & ! number of elements in x
    ny=3                  ! number of elements in y



! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: solsupg, soliter
  type(sysvector_t) :: dsol, rhsd
  type(subscript_t) :: velx, vely, vel, pres, solsc
  type(subscript_t), dimension(nmodes) :: cvalxx, cvalxy, cvalyy, cvalzz, cval
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options

  type(subscriptvec_t) :: cxx, cxy, cyy, czz
  type(vector_t) :: ctensor

  integer :: step, iter, m, i, nsteps(ndU), numsteps, j, nn
  integer :: vertices(4) = [1,3,5,7]
  real(dp) :: alpha, dflowrate, epsconf, U, H
  real(dp) :: epsu, epsp, epsb, Wi
  real(dp), dimension(:), allocatable :: Uspec, epsc

  timer = .false.

! number of incremental steps in U to perform
  nsteps(1) = nint((Ue(1))/dU(1)) ! first interval
  do i = 2, ndU
    nsteps(i) = nint((Ue(i)-Ue(i-1))/dU(i)) ! second and higher intervals
  end do
  numsteps = sum(nsteps)

  if ( printscreen ) print *, 'numsteps = ', numsteps

! U values specified
  allocate( Uspec(0:numsteps), epsc(numsteps) )

  Uspec = [ 0.0_dp, (i*dU(1), i=1,nsteps(1)), &
                    ((Ue(j-1)+dU(j)*i,i=1,nsteps(j)),j=2,ndU)  ]

  if ( printscreen ) print *, 'Uspec = ', Uspec(1:numsteps)

  epsc = [ ( ( epsconfin(j), i=1,nsteps(j) ), j=1,ndU ) ]


! set some parameters

  alpha = eta_p   ! DEVSS parameter
!  alpha = G*deltat   ! DEVSS parameter


! fill coefficients

  if ( alam_model == 4 ) then
    nn = 1
  else if ( alam_model == 5 ) then
    nn = 3
  end if

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+(2+nn)*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
       physqvel, physqpress, 0,     physqgrad, gauss,  &
       gaussb,   cintpl,     0,     0,         0,      &
       0,        0,          model, nmodes,    startm, &
       logc,     timeint,   ( 0, i = 23, 150 )  &
    ]

  coefficients%i(80) = alam_model ! adapted lambda
  coefficients%i(83) = physqc ! physical quantity of first mode
  coefficients%i(84) = 3  ! storage of c tensor: all modes in sysvector
  coefficients%i(86) = 2  ! sysvector number for velocity in SUPG
  coefficients%i(87) = 1  ! sysvector number for iteration
  coefficients%i(89) = 1  ! exclude time-derivative in implicit_ce_supg_elem

  if ( JacobianSUPG ) coefficients%i(91) = 1  ! Jacobian of SUPG

  coefficients%r(1:500) = &
     [ eta_s, 0._dp, 0._dp,  alpha,     0._dp, &
       0._dp, 0._dp, deltat, beta_supg, 0._dp, &
       ( 0._dp, i = 11, 500 ) ]

  if ( alam_model == 4 ) then
    coefficients%r(501:503) = [ ( G/nmodes, lambda, tau_ref, m=1,nmodes ) ]
  else if ( alam_model == 5 ) then
    coefficients%r(501:505) = &
              [ ( G/nmodes, lambda, tau_ref1, f1, tau_ref2, m=1,nmodes ) ]
  end if

! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )

  H = mesh%coor(mesh%points(4),2) ! height is given by y-coordinate of P4

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=6+nmodes, &
     nphysq=3+nmodes )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 4  ! G
  input_probdef%vec_elementdof(1)%a(:,2) = 2         ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(vertices,4:3+nmodes) = ncompc  ! c
  input_probdef%vec_elementdof(1)%a(:,4+nmodes) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5+nmodes) = 4  ! tensor for plotting
  input_probdef%vec_elementdof(1)%a(:,6+nmodes) = 3  ! tensor for plotting 2D

  input_probdef%physq = [1,2,3,(3+m,m=1,nmodes)]
  input_probdef%probnr = 1

  if ( physqmask ) then

!   discard certain zero blocks from the sparse matrix

    if ( physqmask_grad ) call fill_grad_physqmask         ! G-c block
    if ( physqmask_press ) call fill_press_physqmask       ! p-c and c-p blocks
    if ( physqmask_offdiag ) call fill_offdiag_physqmask   ! off-diagonal c

  end if

! define essential boundaries

! center line
  call define_essential ( mesh, input_probdef, &
    curve1=1, physq=physqvel, degfd=[0,1] )

! wall
  call define_essential ( mesh, input_probdef, &
    curve1=3, physq=physqvel )

! pressure in point 1
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=4, nglobalc=1 )

! constraints for periodical boundary conditions

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=2, curve2=5, discretization='weak', &
    elementdof=[2,0,2] )

! gradients
  call define_constraint ( mesh, input_probdef, &
    physq=physqgrad, curve1=2, curve2=5, discretization='collocation' )

! conformation tensor (use collocation)
  do m = 1, nmodes
    call define_constraint ( mesh, input_probdef, &
      physq=physqc+m-1, curve1=2, curve2=5, discretization='collocation' )
  end do

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for solution (excluding contraint forces)

  call create_subscript ( mesh, problem, solsc, &
    physqarr=[physqgrad,physqvel,physqpress,(physqc+m-1,m=1,nmodes)] )

! create vector subscripts for the velocity

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )

! create vector subscript for the pressure

  call create_subscript ( mesh, problem, pres, physqarr=[physqpress] )

! create vector subscript for the conformation tensor

  do m = 1, nmodes
    call create_subscript ( mesh, problem, cval(m), physqarr=[physqc+m-1] )
    call create_subscript ( mesh, problem, cvalxx(m), physqarr=[physqc+m-1], &
      degfd=1 )
    call create_subscript ( mesh, problem, cvalxy(m), physqarr=[physqc+m-1], &
      degfd=2 )
    call create_subscript ( mesh, problem, cvalyy(m), physqarr=[physqc+m-1], &
      degfd=3 )
    call create_subscript ( mesh, problem, cvalzz(m), physqarr=[physqc+m-1], &
      degfd=4 )
  end do
  call create_subscript ( mesh, problem, cxx, degfd=1, vec=5+nmodes )
  call create_subscript ( mesh, problem, cxy, degfd=2, vec=5+nmodes )
  call create_subscript ( mesh, problem, cyy, degfd=3, vec=5+nmodes )
  call create_subscript ( mesh, problem, czz, degfd=4, vec=5+nmodes )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, solsupg, soliter )
  call create_sysvector ( problem, dsol, rhsd )

! create a vector for conformation
! tensor for post processing

  call create_vector ( problem, ctensor, vec=5+nmodes )

! fill solution vector with essential boundary conditions

  soliter%u = 0

! initial solution conformation tensor

  do m = 1, nmodes
    if ( logc == 0 ) then ! standard
      soliter%u(cvalxx(m)%s) = 1 ! initial cxx
      soliter%u(cvalxy(m)%s) = 0 ! initial cxy
      soliter%u(cvalyy(m)%s) = 1 ! initial cyy
      soliter%u(cvalzz(m)%s) = 1 ! initial cyy
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

  if ( printscreen ) print *, 'nnz = ', sysmatrix%Suu%nnz

! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nprob=1 )


! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => soliter
  oldvectors_ve%s(2)%p => solsupg
  oldvectors_ve%p(1)%p => problem  ! needed for implicit_ce_supg_elem


!  open ( unit=13, file='iter.out', recl=300 )
!  open ( unit=14, file='Kdrag.out', recl=300 )

  solver_options%real_storage=rs
  solver_options%integer_storage=is

  call tic


! stepping through specified U values and continue from previous solution

  do step = 1, numsteps

!   average velocity in this step
    U = Uspec(step)

!   Weisenberg number
    Wi = lambda * U / H

!   increase in flow rate in (half) the channel
    dflowrate = H*(Uspec(step) - Uspec(step-1))

    coefficients%r(6) = dflowrate
    coefficients%r(10) = U ! global scaling velocity for SUPG, when needed

    write(*,*) 'step = ', step, ' U = ', U, 'Wi = ', Wi
    if ( printscreen) print *, 'step = ', step, ' U = ', U, 'Wi = ', Wi

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

      if ( iter == 2 ) coefficients%r(6) = 0  ! set flowrate in dsol to 0

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

        call build_c

      end do

      call toc ( 'build_c' )

      if ( .not. physqmask .or. &
              physqmask .and. .not. physqmask_offdiag ) then

!       set to zero off-diagonal conformation blocks

        call build_offdiag_c

        call toc ( 'build_offdiag_c' )

      end if

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, dsol, rhsd )

      call check ( sysmatrix )

!     solve

      call solve_system_ma41 ( sysmatrix, rhsd, dsol, &
        solver_options=solver_options  )

      call toc ( 'solve' )

      soliter%u(solsc%s) = soliter%u(solsc%s) + dsol%u(solsc%s)

!     convergence test

      epsu = maxval(abs(dsol%u(vel%s))) / U
      epsp = maxval(abs(dsol%u(pres%s))) * H / ( U * eta_0 )
      epsb = maxval(abs(dsol%u(cval(1)%s)))

      if ( printscreen) print *, iter, epsu, epsp, epsb

      if ( maxval([ epsu, epsp, epsb ]) < epsconf ) exit

    end do

    call toc ( 'after iteration loop' )


!   copy solution to velocity for supg for next step

    call copy ( soliter, solsupg )

!   write max and mean values of conformation tensor to a file

    call derive_vector ( mesh, problem, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

!    if ( step == 1 ) then
!      open(unit=11, recl=300, status='replace', file='cval.out')
!    else
!      open(unit=11, recl=300, position='append', file='cval.out')
!    end if
    write(*, fmt=*) step, maxval(ctensor%u(cxx%s)), &
                           maxval(ctensor%u(cxy%s)), &
                           maxval(ctensor%u(cyy%s)), &
                           maxval(ctensor%u(czz%s)), &
                           sum(ctensor%u(cxx%s))/size(cxx%s), &
                           sum(ctensor%u(cxy%s))/size(cxy%s), &
                           sum(ctensor%u(cyy%s))/size(cyy%s), &
                           sum(ctensor%u(czz%s))/size(czz%s), &
                           maxval(soliter%u(velx%s)), &
                           maxval(soliter%u(vely%s))
    close(unit=11)

    call toc ( 'one step' )

  end do

  call toc ( 'all steps' )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( solsupg, soliter, dsol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

  call delete ( oldvectors_ve )
  do m = 1, nmodes
    call delete ( cvalxx(m), cvalxy(m), cvalyy(m), cvalzz(m), cval(m) )
  end do
  call delete ( solsc )
  call delete ( cxx, cxy, cyy, czz )

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

end program poiseuille6
