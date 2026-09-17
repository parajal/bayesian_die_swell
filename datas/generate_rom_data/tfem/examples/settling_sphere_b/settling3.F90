#if METIS5

! Sequence of steady states of settling sphere for a b-representation
! full 3D problem on an angular cutout of the cylinder/sphere
! differential Giesekus model.
! Iteration to steady-state with Newton-Raphson.
! Compute drag on sphere using reaction forces.
! Optionally include Jacobian of SUPG test function.
! Optionally discard certain zero blocks from the sparse matrix.
! This a problem similar to settling2, but now it uses metis5 for renumbering
! the mesh

program settling3

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use io_utils_m
  use timer_m
  use write_vtks_m
  use metis5_m
  use math_defs_m

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
    gauss = 8,      & ! 8 point integration of tets
    gaussb = 3,     & ! 3 point integration of boundary elements
    coorsys = 0,    & ! 3D coordinate system
    timeint3 = 10,  & ! BDF2, for steady time discretization
    maxnumiterations = 30, & ! maximum number of Newton-Raphson iterations
    ncompb = 9,     & ! number of b-tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 3         ! model number 3: Giesekus

  logical :: &
    write_vtk = .true.    ! write vtk files

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
    ndU = 1           ! number intervals with different dU increment

  real(dp), parameter :: &
!   increments in imposed average velocity
    dU(ndU) = [ 0.2_dp ], &
!   end values of the intervals in U
    Ue(ndU) = [ 2._dp ], &
!   Newton-Raphson convergence threshold for each interval
    epsconfin(ndU) = [ 1e-6_dp ]

  logical :: &
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
  type(subscript_t) :: velx, vely, velz, vel, pres, solsc, velzsphere
  type(subscript_t) :: solres, solresu, solresb
  type(subscript_t), dimension(nmodes) :: bval
  type(subscript_t), dimension(ncompb,nmodes) :: bvalcmp
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options_ma41

  type(subscriptvec_t) :: cxx, cxy, cyy, czz
  type(vector_t) :: ctensor, btensor

  integer :: step, iter, m, i, nsteps(ndU), numsteps, j, ipost
  integer :: vertices(4) = [1,3,5,10]
  real(dp) :: alpha, epsconf, U, K_drag, res, resu, resb
  real(dp) :: epsu, epsp, epsb, Wi
  real(dp), dimension(:), allocatable :: Uspec, epsc
  character(len=399) :: filename

  timer = .false.


  !call omp_set_num_threads(2) ! number of threads (overrules OMP_NUM_THREADS)


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
      0,     timeint3,  coorsys, ( 0, i = 24, 150 )  &
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
! Surface 1: inflow surface
! Surface 2: outflow surface
! Surface 3: particle surface
! Surface 4: front surface
! Surface 5: back surface
! Surface 6: cylinder wall

  call read_mesh_gmsh ( mesh, filename='meshes/meshp_3d.msh', ndim=3 )

  call fill_mesh_parts ( mesh )

  call renumber_metis ( mesh )

  call toc ('renumber_metis')

  if ( printscreen ) call printinfo ( mesh, printlevel=1 )

  !call write_geometries_vtk ( mesh, write_normals=.true. )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=6+nmodes, &
    nphysq=3+nmodes )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 9  ! G
  input_probdef%vec_elementdof(1)%a(:,2) = 3         ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(vertices,4:3+nmodes) = ncompb  ! b
  input_probdef%vec_elementdof(1)%a(:,4+nmodes) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5+nmodes) = 6  ! tensor for plotting
  input_probdef%vec_elementdof(1)%a(:,6+nmodes) = 9  ! tensor for plotting

  input_probdef%physq = [1,2,3,(3+m,m=1,nmodes)]
  input_probdef%probnr = 1

  if ( physqmask ) then

!   discard certain zero blocks from the sparse matrix

    if ( physqmask_grad ) call fill_grad_physqmask         ! G-b block
    if ( physqmask_press ) call fill_press_physqmask       ! p-b and b-p blocks
    if ( physqmask_offdiag ) call fill_offdiag_physqmask   ! off-diagonal b

  end if

! define essential boundaries

! inflow boundary
  call define_essential ( mesh, input_probdef, surface1=1, physq=physqvel )
! outflow boundary
  call define_essential ( mesh, input_probdef, surface1=2, physq=physqvel )
! sphere boundary
  call define_essential ( mesh, input_probdef, surface1=3, physq=physqvel )
! bottom boundary (in xz-plane)
  call define_essential ( mesh, input_probdef, surface1=4, physq=physqvel, &
    degfd=[0,1,0] )
! top boundary (in zy-plane)
  call define_essential ( mesh, input_probdef, surface1=5, physq=physqvel, &
    degfd=[1,0,0] )
! cylinder boundary (in zy-plane)
  call define_essential ( mesh, input_probdef, surface1=6, physq=physqvel )

  call define_transformation ( mesh, input_probdef, surface=5, physq=physqvel, &
    normalvector=1, v2=[0._dp,0._dp,1._dp], excludesurfaces=[1,2,3,4,6] )

! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! inflow
  call define_essential ( mesh, input_probdef, surface1=1, physq=physqb )

  call problem_definition ( input_probdef, mesh, problem )

  ! call printinfo_problem ( problem, printlevel=4 )

  call build_transformation_matrix ( mesh, problem )

! create vector subscripts for solution (excluding contraint forces)

  call create_subscript ( mesh, problem, solsc, &
    physqarr=[physqgrad,physqvel,physqpress,(physqb+m-1,m=1,nmodes)] )

! create vector subscripts for the velocity

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, velz, physqarr=[physqvel], degfd=3 )
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )
! velocity degrees in z-direction on the cylinder
  call create_subscript ( mesh, problem, velzsphere, physqarr=[physqvel], &
    degfd=3, surfaces=[3] )

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

  call create_subscript ( mesh, problem, solres, &
    physqarr=[physqgrad,physqvel,physqpress,(physqb+m-1,m=1,nmodes)], &
    essentialpart=.false. )

  call create_subscript ( mesh, problem, solresu, physqarr=[physqvel], &
    essentialpart=.false. )

  call create_subscript ( mesh, problem, solresb, &
    physqarr=[(physqb+m-1,m=1,nmodes)], essentialpart=.false. )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, solsupg, soliter )
  call create_sysvector ( problem, dsol, rhsd, reacf )


! create a vector for conformation
! tensor for post processing

  call create_vector ( problem, ctensor, vec=5+nmodes )
  call create_vector ( problem, btensor, vec=6+nmodes )

! fill solution vector with essential boundary conditions

  soliter%u = 0
  dsol%u = 0

! initial solution b tensor

  do m = 1, nmodes
    soliter%u(bvalcmp(1,m)%s) = 1  ! initial bxx
    soliter%u(bvalcmp(2,m)%s) = 0  ! initial bxy
    soliter%u(bvalcmp(3,m)%s) = 0  ! initial bxz
    soliter%u(bvalcmp(4,m)%s) = 0  ! initial byx
    soliter%u(bvalcmp(5,m)%s) = 1  ! initial byy
    soliter%u(bvalcmp(6,m)%s) = 0  ! initial byz
    soliter%u(bvalcmp(7,m)%s) = 0  ! initial bzx
    soliter%u(bvalcmp(8,m)%s) = 0  ! initial bzy
    soliter%u(bvalcmp(9,m)%s) = 1  ! initial bzz
  end do

  call copy(soliter,solsupg)
  !call copy(soliter,dsol)


! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    usephysqmask=physqmask )
  !call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

  if ( printscreen ) print *, 'nnz = ', sysmatrix%Suu%nnz


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nprob=1 )


! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => soliter
  oldvectors_ve%s(2)%p => solsupg
  oldvectors_ve%p(1)%p => problem  ! needed for implicit_ce_supg_elem


  open ( unit=15, file='iter.out', recl=300 )
  open ( unit=14, file='Kdrag.out', recl=300 )

  solver_options_ma41%real_storage = rs
  solver_options_ma41%integer_storage = is
  solver_options_ma41%scaling = 1

  call tic

! write vtks

  ipost = 0
  if ( write_vtk ) then
    write(filename,'(a,i4.4,a)') 'flow', ipost, '.vtk'
    call write_vtks ( mesh, problem, soliter, coefficients, filename )
    ipost = ipost + 1
  end if

! stepping through specified U values and continue from previous solution

  do step = 1, numsteps

!   average velocity in this step
    U = Uspec(step)

!   Weisenberg number
    Wi = lambda * U / R

    coefficients%r(10) = U ! global scaling velocity for SUPG, when needed

    write(15,*) 'step = ', step, ' U = ', U, 'Wi = ', Wi
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

      if ( iter == 1 ) then
        call fill_sysvector ( mesh, problem, dsol, physq=physqvel, degfd=3, &
          surface1=6, value=(Uspec(step) - Uspec(step-1)) )
        call fill_sysvector ( mesh, problem, dsol, physq=physqvel, degfd=3, &
          surface1=1, surface2=2, value=(Uspec(step) - Uspec(step-1)) )
      else
        dsol%u = 0
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
        solver_options=solver_options_ma41  )

      call toc ( 'solve' )

      call transform_to_global ( problem, dsol )

      soliter%u(solsc%s) = soliter%u(solsc%s) + dsol%u(solsc%s)

      res = sqrt ( ( 1._dp/size(solres%s,1) ) * &
                     dot_product ( rhsd%u(solres%s), rhsd%u(solres%s) ) )

      resu = ( R**2 / ( eta_0 * U ) ) * sqrt ( ( 1._dp/size(solresu%s,1) ) * &
                     dot_product ( rhsd%u(solresu%s), rhsd%u(solresu%s) ) )

      resb = ( R / U ) * sqrt ( ( 1._dp/size(solresb%s,1) ) * &
                     dot_product ( rhsd%u(solresb%s), rhsd%u(solresb%s) ) )

!     convergence test

      epsu = maxval(abs(dsol%u(vel%s))) / U
      epsp = maxval(abs(dsol%u(pres%s))) * R / ( U * eta_0 )
      epsb = maxval(abs(dsol%u(bval(1)%s)))

      if ( printscreen) print *, iter, epsu, epsp, epsb
      write(15,*) iter, epsu, epsp, epsb, res, resu, resb

      if ( maxval([ epsu, epsp, epsb ]) < epsconf ) exit

    end do

    call toc ( 'after iteration loop' )

!   write vtks

    if ( write_vtk ) then
       write(filename,'(a,i4.4,a)') 'flow', ipost, '.vtk'
       call write_vtks ( mesh, problem, soliter, coefficients, filename )
       ipost = ipost + 1
    end if

!   reaction forces

    call reaction_forces ( problem, sysmatrix, dsol, rhsd, reacf )

    call transform_to_global ( problem, reacf )

    K_drag = - 8 * sum ( reacf%u(velzsphere%s) ) / ( 6*pi*eta_0*U*R )

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
      open(unit=11, recl=800, status='replace', file='cval.out')
    else
      open(unit=11, recl=800, position='append', file='cval.out')
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
                           maxval(soliter%u(velz%s))
    close(unit=11)

    call toc ( 'one step' )

  end do

  call toc ( 'all steps' )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( solsupg, soliter, dsol, rhsd, reacf )
  call delete ( sysmatrix )
  call delete ( coefficients )

  call delete ( oldvectors_ve )
  do m = 1, nmodes
    do i = 1, ncompb
      call delete ( bvalcmp(i,m), bval(m) )
    end do
  end do
  call delete ( cxx, cxy, cyy, czz )

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


end program settling3

#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on metis5', &
    ' - add metis5 lib for linking', &
    ' - set preprocessing macro METIS5 in Mdefs.mk'
end
#endif
