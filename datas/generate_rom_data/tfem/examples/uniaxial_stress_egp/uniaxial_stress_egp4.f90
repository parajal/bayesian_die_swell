! Uniaxial compression/extension using a single element in axisymmetric flow.
! Compressible EGP model with compressible flow.
! 14-mode model with 12 alpha process viscoelastic modes, one elastic mode for
! strain hardening and one mode for the beta process.
! Data for PLLA from van Breemen et al. (2012) DOI: 10.1002/polb.23199
! Pressure (J) dependent relaxation term using adapted lambda.
! Plastic strain softening.
! Contravariant deformation tensor formulation (b-tensor).
! Fully-implicit with Newton-Raphson.
! Optionally include Jacobian of SUPG test function.
! Optionally discard certain zero blocks from the sparse matrix.

program uniaxial_stress_egp4

  use tfem_m
  use viscoelastic_elements_m
  use compressible_fluid_elements_m
  use hsl_ma41_m
  use io_utils_m
  use timer_m
  use material_data_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities
    pintpl = 4,     & ! Q1 pressures
    gintpl = 4,     & ! Q1 gradients
    bintpl = 4,     & ! Q1 b-tensor
    gauss = 3,      & ! 3x3 integration of quads
    gaussb = 3,     & ! 3 point integration of boundary elements
    nz=1,           & ! number of elements in z
    nr=1,           & ! number of elements in r
    timeint1 = 8,   & ! BDF1, first time step
    timeint2 = 10,  & ! BDF2, after first time step
    numtimesteps = 120, & ! number of time steps
    maxnumiterations = 30, & ! maximum number of Newton-Raphson iterations
    ncompb = 5,     & ! number of b-tensor components
    nmodes = 3,     & ! number of modes
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    physqb = 4,     & ! physical quantity nr of the b-tensor (mode 1)
    physqgammap = 4+nmodes,& ! physical quantity nr of the equiv. plastic strain
    am = nmodes-2,  & ! 1:am range of modes for alpha process
    em = nmodes-1,  & ! elastic mode
    bm = nmodes,    & ! mode of beta process
    startm = 501,   & ! start of material model data
    endm   = 700,   & ! end of material model data
    model = 24,     & ! EGP compressible
    alamJ_model = 1,& ! adapted lambda model for h(J)
    alam_gammap_model = 1,& ! adapted lambda model for g(gammap)
    alam_model = 4, & ! adapted lambda model: 4: Eyring
    every = 1,    & ! output every every steps
    coorsys = 1       ! axisymmetric coordinate system

  real(dp) :: &
    eta_p,    & ! polymer viscosity
    eta_s,    & ! solvent viscosity
    eta_0,    & ! zero-shear viscosity
    lambda_0    ! effective relaxation time for defining Wi

  real(dp), parameter :: &
    Kmod = 3.5e3_dp,     & ! compression modulus
    p0 = 0._dp,         & ! reference pressure where J=1.
    beta_s = 0.0_dp,    & ! beta viscosity parameter = eta_s/(eta_s+eta_p)
    R = 1._dp,          & ! radius of the cylinder
    L = 1._dp,          & ! length of the cylinder
    epsilondot = -1.e-2_dp ! extension rate of the cylinder

  real(dp), parameter :: &
    deltat = 0.5_dp,    & ! time step for BDF1/BDF2
    epsconf = 1e-8_dp,  & ! Newton-Raphson convergence threshold
    beta_supg = 1,      & ! SUPG factor
    rs = 1.2_dp,        & ! real_storage for LU (HSL)
    is = 1.6_dp           ! integer_storage for LU (HSL)

  logical, parameter :: &
    physqmask = .false., & ! use physical quantity masking
    physqmask_grad = .true., & ! use physqmask for G-b coupling
    physqmask_press = .true., & ! use physqmask for p-b and b-p coupling
    physqmask_offdiag = .false., & ! use physqmask for b-b off-diagonal coupling
    JacobianSUPG = .true., & ! take Jacobian of SUPG test function into account
    printscreen = .true.  ! print on standard output

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: soln, solnm1, soliter
  type(sysvector_t) :: dsol, rhsd
  type(vector_t) :: velocity, pressure, vorticity, Jvol
  type(subscript_t) :: velz, vely, vel, pres, solsc, gp
  type(subscript_t), dimension(nmodes) :: bval
  type(subscript_t), dimension(ncompb,nmodes) :: bvalcmp
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options

  type(subscriptvec_t) :: czz, czr, crr, ctt
  type(vector_t) :: ctensor, btensor, vonmises, gammap, stresstensor

  integer :: step, i, iter, m, mode1, mode2
  integer :: vertices(4) = [1,3,5,7]
  real(dp) :: alpha, epsu, epsp, epsb, epsgp, Wi
  type(vemodel_t) :: vemodel

  timer = .false.


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=350, ncoefr=endm )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gaussb,   bintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      0,     timeint1,  ( 0, i = 23, 350 )  &
    ]

  coefficients%i(23) = coorsys

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

! fill material data

  coefficients%r = 0 ! intialize to zero

  call fill_data ( coefficients, 'PLLA', vemodel )

  alpha = sum(vemodel%modulus)*deltat   ! DEVSS parameter

! viscosity parameters
  eta_p = sum( vemodel%modulus(1:am)*vemodel%lambda(1:am) ) + &
               vemodel%modulus(bm)*vemodel%lambda(bm)    ! polymer viscosity
  eta_s = beta_s/(1-beta_s)*eta_p       ! solvent viscosity
  eta_0 = eta_s+eta_p    ! zero-shear viscosity

  coefficients%r(1:10) = &
    [ eta_s, 0._dp,  0._dp,     alpha, 0._dp, &
      0._dp, 0._dp, deltat, beta_supg, 0._dp ]

  coefficients%r(29) = Kmod
  coefficients%r(30) = p0

! global scaling velocity for SUPG, when needed
  coefficients%r(10) = abs( epsilondot * L )

  coefficients%r(251) = deltat
  coefficients%r(252) = Kmod
  coefficients%r(253) = p0

! effective lambda
  lambda_0 = ( sum( vemodel%modulus(1:am)*vemodel%lambda(1:am)**2 ) + &
                    vemodel%modulus(bm)*vemodel%lambda(bm)**2 ) / eta_p

  print *, 'eta_p, eta_s, eta_0, lambda_0'
  print *, eta_p, eta_s, eta_0, lambda_0


! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nz
  meshgen_options%ny = nr
  meshgen_options%lx = L
  meshgen_options%ly = R

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )

  call printinfo ( mesh, printlevel=1 )


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

! axis of cylinder
  call define_essential ( mesh, input_probdef, curves=[1], degsfd=[2], &
    physq=physqvel )

! cylinder ends
  call define_essential ( mesh, input_probdef, curves=[2,4], degsfd=[1], &
    physq=physqvel )

! constraints for periodical boundary conditions

! gradients
  call define_constraint ( mesh, input_probdef, &
    physq=physqgrad, curve1=2, curve2=5, discretization='collocation' )

! conformation
  do m = 1, nmodes
    call define_constraint ( mesh, input_probdef, &
      physq=physqb-1+m, curve1=2, curve2=5, discretization='collocation' )
  end do

! gammap
  call define_constraint ( mesh, input_probdef, &
    physq=physqgammap, curve1=2, curve2=5, discretization='collocation' )

  call problem_definition ( input_probdef, mesh, problem )

 ! create vector subscripts for solution (excluding contraint forces)

  call create_subscript ( mesh, problem, solsc, &
    physqarr = &
      [physqgrad,physqvel,physqpress,physqgammap,(physqb+m-1,m=1,nmodes)] )

! create vector subscripts for the velocity

  call create_subscript ( mesh, problem, velz, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )

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
  call create_subscript ( mesh, problem, czz, degfd=1, vec=6+nmodes )
  call create_subscript ( mesh, problem, czr, degfd=2, vec=6+nmodes )
  call create_subscript ( mesh, problem, crr, degfd=3, vec=6+nmodes )
  call create_subscript ( mesh, problem, ctt, degfd=4, vec=6+nmodes )

! create vector subscript for the equivalent plastic strain

  call create_subscript ( mesh, problem, gp, physqarr=[physqgammap] )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, soln, solnm1, soliter )
  call create_sysvector ( problem, dsol, rhsd )


! create a vector for conformation
! tensor for post processing

  call create_vector ( problem, stresstensor, vec=6+nmodes )
  call create_vector ( problem, ctensor, vec=6+nmodes )
  call create_vector ( problem, btensor, vec=7+nmodes )

! fill solution vector with essential boundary conditions

  soln%u = 0

! initial pressure

  soln%u(pres%s) = p0

! initial solution b tensor

  do m = 1, nmodes
    soln%u(bvalcmp(1,m)%s) = 1 ! initial bzz
    soln%u(bvalcmp(2,m)%s) = 0 ! initial bzr
    soln%u(bvalcmp(3,m)%s) = 0 ! initial brz
    soln%u(bvalcmp(4,m)%s) = 1 ! initial brr
    soln%u(bvalcmp(5,m)%s) = 1 ! initial btt
  end do

  call copy(soln,solnm1)
  call copy(soln,soliter)

  dsol%u = 0

! set velocity = epsilondot*L on cylinder end
  call fill_sysvector ( mesh, problem, dsol, &
    curve1=2, physq=physqvel, degsfd=[1], value=epsilondot*L )


! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    usephysqmask=physqmask )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

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
  Wi = abs ( lambda_0 * epsilondot )

  print *, ' epsilondot = ', epsilondot, 'Wi = ', Wi
  print *


  open ( unit=13, file='iter.out', recl=300 )
  open ( unit=14, file='stress.out', recl=300 )

  solver_options%real_storage=rs
  solver_options%integer_storage=is

  call tic

! time stepping

  do step = 1, numtimesteps

    write(13,*) 'step = ', step
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

    iter = 0

    do

      iter = iter + 1

      if ( step == 1 .and. iter == 2 ) then
 !      set dsol = 0 on end plate
        call fill_sysvector ( mesh, problem, dsol, &
          curve1=2, physq=physqvel, degsfd=[1], value=0._dp )
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

!     alpha process of modes 1 to am
      mode1 = 1; mode2 = am
      call build_b

!     elastic mode em
      mode1 = em; mode2 = em
      call build_b

!     beta process mode bm
      mode1 = bm; mode2 = bm
      call build_b

      call toc ( 'build_b' )

      if ( .not. physqmask .or. &
              physqmask .and. .not. physqmask_offdiag ) then

!       set to zero off-diagonal conformation blocks

        call build_offdiag_b

        call toc ( 'build_offdiag_b' )

      end if

!     build (assemble) matrix/vector for gammap part

!     use alpha process of modes 1 to am for vonmises stress
      mode1 = 1; mode2 = am
      call build_gammap

      call toc ( 'build_gammap' )

!     periodical condition on gradients, conformation and gammap

      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        elemsub=stokes_constr_node_conn, addmatvec=.true. )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, dsol, rhsd )

      call check ( sysmatrix )

!     solve

      call solve_system_ma41 ( sysmatrix, rhsd, dsol, &
        solver_options=solver_options  )

      call toc ( 'solve' )

      soliter%u(solsc%s) = soliter%u(solsc%s) + dsol%u(solsc%s)

!     convergence test

      epsu = maxval(abs(dsol%u(vel%s))) / abs( epsilondot * L )
      epsp = maxval(abs(dsol%u(pres%s))) / abs( sum(vemodel%modulus) )
      epsb = 0
      do m = 1, nmodes
        epsb = max ( epsb, maxval(abs(dsol%u(bval(m)%s))) )
      end do
      epsgp = maxval(abs(dsol%u(gp%s)))

      if ( printscreen) print *, iter, epsu, epsp, epsb, epsgp
      write(13,*) iter, epsu, epsp, epsb, epsgp

      if ( maxval([ epsu, epsp, epsb, epsgp ]) < epsconf ) exit

    end do

    call toc ( 'after iteration loop' )


!   copy solution to older time step for next time step

    call copy ( soln, solnm1 )
    call copy ( soliter, soln )


!   output data

    if ( mod ( step, every ) == 0 ) call output


    call toc ( 'one step' )

  end do

  call toc ( 'all steps' )


! post-processing

  call create_vector ( problem, velocity, physq=2 )
  call create_vector ( problem, pressure, vec=5+nmodes )
  call create_vector ( problem, Jvol, vec=5+nmodes )
  call create_vector ( problem, vorticity, vec=5+nmodes )

  call extract_physvector ( mesh, problem, soln, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! derive vectors

  oldvectors%s(1)%p => soln

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

! compute J directly from the pressure
  Jvol%u = exp(-(pressure%u-p0)/Kmod)

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call write_vector_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='velocity', vector=velocity )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='pressure', vector=pressure, append=.true. )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='J', vector=Jvol, append=.true. )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='vorticity', vector=vorticity, append=.true. )

  call printtofile ( mesh, problem, filename='velocity_cl.out', curve=11, &
    vector=velocity )
  call printtofile ( mesh, problem, filename='pressure_cl.out', curve=12, &
    vector=pressure )

! stress tensor

  coefficients%i(28) = 0 ! all modes added

  call derive_vector ( mesh, problem, stresstensor, &
    elemsub=deriv_viscoelastic_stress_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problem, filename='stress.vtk', &
    dataname='stress_tensor', vector=stresstensor, assume33=.true. )

  call delete ( stresstensor )

! conformation tensor

  coefficients%i(28) = em ! use elastic mode

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

  call create_vector ( problem, vonmises, vec=5+nmodes )

  coefficients%i(13)=1 ! von Mises
  coefficients%i(28)=0 ! global von Mises

! alpha process of modes 1 to am
  mode1 = 1; mode2 = am
! set modes
  coefficients%i(26:27) = [ mode1, mode2 ]

  call derive_vector ( mesh, problem, vonmises, &
    elemsub=deriv_viscoelastic_stress_scalar, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_scalar_vtk ( mesh, problem, filename='c.vtk', &
    dataname='von_mises_alpha', vector=vonmises, append=.true. )

! beta process of modes bm to bm
  mode1 = 1; mode2 = bm
! set modes
  coefficients%i(26:27) = [ mode1, mode2 ]

  call derive_vector ( mesh, problem, vonmises, &
    elemsub=deriv_viscoelastic_stress_scalar, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_scalar_vtk ( mesh, problem, filename='c.vtk', &
    dataname='von_mises_beta', vector=vonmises, append=.true. )

! reset modes to default
  coefficients%i(26:27) = 0

  call delete ( vonmises )

! gammap

  call create_vector ( problem, gammap, vec=5+nmodes )

  call derive_vector ( mesh, problem, gammap, &
    elemsub=deriv_gammap, coefficients=coefficients, &
    oldvectors=oldvectors_ve )

  call write_scalar_vtk ( mesh, problem, filename='gammap.vtk', &
    dataname='gammap', vector=gammap )

  call delete ( gammap )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( soln, solnm1, soliter, dsol, rhsd )
  call delete ( velocity, pressure, vorticity, Jvol )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

  call delete ( oldvectors_ve )
  do m = 1, nmodes
    do i = 1, ncompb
      call delete ( bvalcmp(i,m), bval(m) )
    end do
  end do
  call delete ( czz, czr, crr, ctt )
  call delete ( gp )

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

  end subroutine build_vpG


  subroutine build_b

!   set modes of subprocess
    coefficients%i(26:27) = [ mode1, mode2 ]

!   build (assemble) matrix and vector for b-tensor mode m

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

  subroutine output

!   write max and mean values of conformation and b-tensor tensor to a file

    call derive_vector ( mesh, problem, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    if ( step == 1 ) then
      open(unit=11, recl=600, status='replace', file='cval.out')
    else
      open(unit=11, recl=600, position='append', file='cval.out')
    end if
    write(11, fmt=*) step * deltat, maxval(ctensor%u(czz%s)), &
                                    maxval(ctensor%u(czr%s)), &
                                    maxval(ctensor%u(crr%s)), &
                                    maxval(ctensor%u(ctt%s)), &
                                    sum(ctensor%u(czz%s))/size(czz%s), &
                                    sum(ctensor%u(czr%s))/size(czr%s), &
                                    sum(ctensor%u(crr%s))/size(crr%s), &
                                    sum(ctensor%u(ctt%s))/size(ctt%s), &
                                    maxval(soln%u(velz%s)), &
                                    maxval(soln%u(vely%s)), &
                                    maxval(soln%u(pres%s)), &
                                    minval(soln%u(pres%s)), &
                                    maxval(soln%u(gp%s)), &
                                    minval(soln%u(gp%s))
    close(unit=11)

    if ( printscreen ) &
         print *, 'step = ', step, 'max czz = ', maxval(ctensor%u(czz%s))

!   stress tensor

    coefficients%i(28) = 0 ! all modes added

    call derive_vector ( mesh, problem, stresstensor, &
      elemsub=deriv_viscoelastic_stress_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    if ( printscreen ) &
         print *, 'step = ', step, 'time = ', step*deltat, &
         'ln strain = ', epsilondot*step*deltat, &
         'strain = ', exp(epsilondot*step*deltat), &
         'tauzz - taurr = ', maxval(stresstensor%u(czz%s) &
                                             - stresstensor%u(ctt%s)), &
          'taurr =', maxval(stresstensor%u(crr%s)), &
          'pressure =', maxval(soln%u(pres%s)), &
          'J =', maxval(exp(-(soln%u(pres%s)-p0)/Kmod))

    write(14,*) step, step*deltat, epsilondot*step*deltat, &
      exp(epsilondot*step*deltat), &
      maxval(stresstensor%u(czz%s) - stresstensor%u(ctt%s)), &
      maxval(stresstensor%u(crr%s)), &
      maxval(soln%u(pres%s)), &
      maxval(exp(-(soln%u(pres%s)-p0)/Kmod))

  end subroutine output

end program uniaxial_stress_egp4
