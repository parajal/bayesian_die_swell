! Sequence of steady states of flow around a cylinder for a b-representation
! differential Giesekus model.
! Iteration to steady-state with Newton-Raphson.
! Periodical boundary conditions.
! Compute drag on cylinder using reaction forces.
! Optionally include Jacobian of SUPG test function.
! Optionally discard certain zero blocks from the sparse matrix.

program cylinder8

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use io_utils_m
  use figplot_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,     & ! P2 velocities
    pintpl = 2,     & ! P1 pressures
    gintpl = 2,     & ! P1 gradients
    bintpl = 2,     & ! P1 b-tensor
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    physqb = 4,     & ! physical quantity nr of the b-tensor (mode 1)
    gauss = 6,      & ! 6 point integration of triangles
    gaussb = 3,     & ! 3 point integration of boundary elements
    timeint3 = 10,  & ! BDF2, for steady time discretization
    maxnumiterations = 30, & ! maximum number of Newton-Raphson iterations
    ncompb = 4,     & ! number of b-tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 3         ! model number 3: Giesekus

  real(dp), parameter :: &
    G     = 1.0_dp,    & ! modulus
    lambda = 1.0_dp,   & ! relaxation time
    eta_p = lambda*G,  & ! polymer viscosity
    beta_s = 0.59_dp,  & ! beta viscosity parameter = eta_s/(eta_s+eta_p)
    eta_s = beta_s/(1-beta_s)*eta_p,    & ! solvent viscosity
    eta_0 = eta_s+eta_p, & ! zero-shear viscosity
    mobility = 0.01_dp,& ! mobility parameter
    H = 2._dp,         & ! (half-) height of the channel (must match mesh)
    R = 1._dp            ! radius of the cylinder (must match mesh)

  real(dp), parameter :: &
    deltat3 = 0.05_dp,   & ! time step for steady iteration scheme
    beta_supg = 1,       & ! SUPG factor
    rs = 1.2_dp,         & ! real_storage for LU (HSL)
    is = 1.6_dp            ! integer_storage for LU (HSL)

  integer, parameter :: &
    ndU = 5           ! number intervals with different dU increment

  real(dp), parameter :: &
!   increments in imposed average velocity
    dU(ndU) = [ 0.1_dp, 0.4_dp, 0.5_dp, 0.75_dp, 0.5_dp ], &
!   end values of the intervals in U
    Ue(ndU) = [ 0.1_dp, 0.5_dp, 1.0_dp, 49.0_dp, 100.0_dp ], &
!   Newton-Raphson convergence threshold for each interval
    epsconfin(ndU) = [ 1e-8_dp, 1e-8_dp, 1e-7_dp, 1e-7_dp, 1e-7_dp ]

  logical :: &
    physqmask = .true., & ! use physical quantity masking
    physqmask_grad = .true., & ! use physqmask for G-b coupling
    physqmask_press = .false., & ! use physqmask for p-b and b-p coupling
    physqmask_offdiag = .true., & ! use physqmask for b-b off-diagonal coupling
    JacobianSUPG = .true., & ! take Jacobian of SUPG test function into account
    printscreen = .true.  ! print on standard output

! definitions

  type(mesh_t) :: mesh, mesh1
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: solsupg, soliter
  type(sysvector_t) :: dsol, rhsd, reacf
  type(vector_t) :: velocity, pressure, vorticity
  type(subscript_t) :: velx, vely, vel, pres, solsc, velxcylinder
  type(subscript_t), dimension(nmodes) :: bval
  type(subscript_t), dimension(ncompb,nmodes) :: bvalcmp
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options

  type(subscriptvec_t) :: cxx, cxy, cyy
  type(vector_t) :: ctensor, btensor

  integer :: step, iter, m, i, nsteps(ndU), numsteps, j
  integer :: vertices(3) = [1,3,5]
  real(dp) :: alpha, dflowrate, epsconf, U, K_drag
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

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+3*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gaussb,   bintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      0,     timeint3,   ( 0, i = 23, 150 )  &
    ]

  coefficients%i(71) = 1  ! CDT for b-formulation
  coefficients%i(83) = physqb ! physical quantity of first mode
  coefficients%i(84) = 3  ! storage of b tensor: all modes in sysvector
  coefficients%i(86) = 2  ! sysvector number for velocity in SUPG
  coefficients%i(87) = 1  ! sysvector number for iteration
  coefficients%i(89) = 1  ! exclude time-derivative in implicit_ce_supg_elem
  coefficients%i(90) = 1  ! rotation reinitialization on element level

  if ( JacobianSUPG ) coefficients%i(91) = 1  ! Jacobian of SUPG

  coefficients%r = &
    [ eta_s, 0._dp,   0._dp,     alpha, 0._dp, &
      0._dp, 0._dp, deltat3, beta_supg, 0._dp, &
      ( 0._dp, i = 11, 500 ), &
      ( G/nmodes, lambda, mobility, m=1,nmodes ) ]

! read mesh

  call read_mesh_gmsh ( mesh1, filename='confined_cylinder1.msh', ndim=2 )
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

  if ( printscreen ) call printinfo ( mesh, printlevel=1 )

! plot mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  plot_options%fontsize = 10

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


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
  input_probdef%vec_elementdof(1)%a(:,5+nmodes) = 3  ! tensor for plotting
  input_probdef%vec_elementdof(1)%a(:,6+nmodes) = 4  ! tensor for plotting

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
      physq=physqb+m-1, curve1=4, curve2=12, discretization='collocation' )
  end do

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
    degfd=1, curves=[9] )

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

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, solsupg, soliter )
  call create_sysvector ( problem, dsol, rhsd, reacf )


! create a vector for conformation
! tensor for post processing

  call create_vector ( problem, ctensor, vec=5+nmodes )
  call create_vector ( problem, btensor, vec=6+nmodes )

! fill solution vector with essential boundary conditions

  soliter%u = 0

! initial solution b tensor

  do m = 1, nmodes
    soliter%u(bvalcmp(1,m)%s) = 1 ! initial bxx
    soliter%u(bvalcmp(2,m)%s) = 0 ! initial bxy
    soliter%u(bvalcmp(3,m)%s) = 0 ! initial byx
    soliter%u(bvalcmp(4,m)%s) = 1 ! initial byy
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


  open ( unit=13, file='iter.out', recl=300 )
  open ( unit=14, file='Kdrag.out', recl=300 )

  solver_options%real_storage=rs
  solver_options%integer_storage=is

  call tic


! stepping through specified U values and continue from previous solution

  do step = 1, numsteps

!   average velocity in this step
    U = Uspec(step)

!   Weisenberg number
    Wi = lambda * U / R

!   increase in flow rate in (half) the channel
    dflowrate = H*(Uspec(step) - Uspec(step-1))

    coefficients%r(6) = dflowrate
    coefficients%r(10) = U ! global scaling velocity for SUPG, when needed

    write(13,*) 'step = ', step, ' U = ', U, 'Wi = ', Wi
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

      epsu = maxval(abs(dsol%u(vel%s))) / U
      epsp = maxval(abs(dsol%u(pres%s))) * R / ( U * eta_0 )
      epsb = maxval(abs(dsol%u(bval(1)%s)))

      if ( printscreen) print *, iter, epsu, epsp, epsb
      write(13,*) iter, epsu, epsp, epsb

      if ( maxval([ epsu, epsp, epsb ]) < epsconf ) exit

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

!   write max and mean values of conformation and b-tensor tensor to a file

    call derive_vector ( mesh, problem, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    if ( step == 1 ) then
      open(unit=11, recl=300, status='replace', file='cval.out')
    else
      open(unit=11, recl=300, position='append', file='cval.out')
    end if
    write(11, fmt=*) step, maxval(ctensor%u(cxx%s)), &
                           maxval(ctensor%u(cxy%s)), &
                           maxval(ctensor%u(cyy%s)), &
                           sum(ctensor%u(cxx%s))/size(cxx%s), &
                           sum(ctensor%u(cxy%s))/size(cxy%s), &
                           sum(ctensor%u(cyy%s))/size(cyy%s), &
                           maxval(soliter%u(velx%s)), &
                           maxval(soliter%u(vely%s))
    close(unit=11)

    call toc ( 'one step' )

  end do

  call toc ( 'all steps' )


! post-processing

  call create_vector ( problem, velocity, physq=2 )
  call create_vector ( problem, pressure, vec=4+nmodes )
  call create_vector ( problem, vorticity, vec=4+nmodes )

  call extract_physvector ( mesh, problem, soliter, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! derive vectors

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

! conformation tensor

  call derive_vector ( mesh, problem, ctensor, &
    elemsub=deriv_conformation_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problem, filename='c.vtk', &
    dataname='conformation_tensor', vector=ctensor )

  call printtofile ( mesh, problem, filename='conformation_cl.out', curve=11, &
    vector=ctensor )

  call delete ( ctensor )

! b tensor

  call derive_vector ( mesh, problem, btensor, &
    elemsub=deriv_conformation_tensor_std, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problem, filename='b.vtk', &
    dataname='b-tensor', vector=btensor, symmetric=.false. )

  call printtofile ( mesh, problem, filename='b-tensor_cl.out', &
    curve=11, vector=btensor )

  call delete ( btensor )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( solsupg, soliter, dsol, rhsd, reacf )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

  call delete ( oldvectors_ve )
  do m = 1, nmodes
    do i = 1, ncompb
      call delete ( bvalcmp(i,m), bval(m) )
    end do
  end do
  call delete ( cxx, cxy, cyy )

contains


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
      input_probdef%physqmask(physqb+j-1,physqpress) = .false.
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
      physqrow=[physqvel], physqcol=[(physqb+m-1,m=1,nmodes)], &
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


  subroutine build_b

!   build (assemble) matrix and vector for b-tensor mode m

    if ( JacobianSUPG ) then

!     add diagonal block of steady-state timed derivative
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=steady_ce_timederiv_supg_elem, coefficients=coefficients, &
        physqrow=[physqb+m-1], physqcol=[physqvel,physqb+m-1], &
        oldvectors=oldvectors_ve, addmatvec=.true. )

!     diagonal block
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
        physqrow=[physqb+m-1], physqcol=[physqvel,physqb+m-1], &
        oldvectors=oldvectors_ve, addmatvec=.true. )

    else

!      add diagonal block of steady-state timed derivative
       call build_system ( mesh, problem, sysmatrix, rhsd, &
         elemsub=steady_ce_timederiv_supg_elem, coefficients=coefficients, &
         physqrow=[physqb+m-1], physqcol=[physqb+m-1], &
         oldvectors=oldvectors_ve, addmatvec=.true. )

!      diagonal block
       call build_system ( mesh, problem, sysmatrix, rhsd, &
         elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
         physqrow=[physqb+m-1], physqcol=[physqb+m-1], &
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
!     set to zero off-diagonal blocks pressure-b
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        physqrow=[physqpress], physqcol=[physqb+m-1], &
        buildvector=.false., zeromatvec=.true. )
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        physqrow=[physqb+m-1], physqcol=[physqpress], &
        buildvector=.false., zeromatvec=.true. )
    end if

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problem, sysmatrix, &
      sysvector=rhsd, constraint1=3+m, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

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

end program cylinder8
