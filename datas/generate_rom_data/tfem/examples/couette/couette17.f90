! Viscoelastic problem on a unit square
! linear Couette flow
! DEVSS-G/SUPG
! second-order time integration
! Fully-implicit with Newton-Raphson for the conformation tensor.
! NOTE: Similar to couette2, but now with implicit conformation.

program couette17

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,     &  ! Q2 velocities
    pintpl = 4,     &  ! Q1 pressures
    gintpl = 4,     &  ! Q1 gradients
    cintpl = 4,     &  ! Q1 conformation
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    gauss = 3,      & ! 3x3 integration of quads
    ncompc = 3,     & ! number of conformation tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 2         ! Oldroyd-B model


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc
  type(problem_t), target :: problem, problemc
  type(sysmatrix_t) :: sysmatrix, sysmatrixc
  type(sysvector_t), target :: soln, solnm1, solhat
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vely, cvalxx, cvalxy, cvalyy, cval, cvalc2, cvalc5
  type(sysvector_t), dimension(nmodes), target :: solcn, solcnm1, solciter
  type(sysvector_t) :: dsolc, rhsc
  type(lu_ma57_t) :: lu_u
  type(solver_options_ma57_t) :: solver_options_u
  type(solver_options_ma41_t) :: solver_options_c


! variables

  logical :: &
    printscreen = .true. ! print on standard output

  integer :: &
    nx=20,             & ! number of elements in x
    ny=20,             & ! number of elements in y
    timeint1 = 8,      & ! first-order time integration, first time step
    timeint2 = 10,     & ! second-order time integration after first time step
    maxnumiterations = 20, & ! maximum number of Newton-Raphson iterations
    numtimesteps = 2,  & ! number of time steps
    logc = 0             ! standard scheme or log transformation

  real(dp) :: &
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 1.0_dp,   & ! relaxation time
    deltat = 5.e-3_dp, & ! time step
    thetapar = 0.55_dp, & ! theta parameter in the theta method
    epsconf = 1e-12_dp, & ! Newton-Raphson convergence threshold
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.0_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 1.4_dp       ! integer_storage for conformation LU (HSL)

  integer :: step, i, iter, m
  real(dp) :: alpha, G, csteady(3), lsteady(2), eigv(2,2)


! namelist for input of variables; read from standard input

  namelist /comppar/ printscreen, nx, ny, timeint1, timeint2, &
    numtimesteps, maxnumiterations, logc, thetapar, epsconf, eta_s, &
    eta_p, lambda, deltat, beta, rs_gup, is_gup, rs_c, is_c

  read ( unit=*, nml=comppar )


! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+2*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gauss,    cintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      logc,     timeint1,    ( 0, i = 23, 150 )  &
    ]

  coefficients%i(84) = 2  ! storage of c tensor: all components, separate modes

  coefficients%r = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
      0._dp, 0._dp,  deltat,   beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G,     lambda &
    ]

  coefficients%r(28) = thetapar

! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=3 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [ 4,0,4,0,4,0,4,0,0,    &  ! G
                  2,2,2,2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,1,0,0,    &  ! pressure
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                  [9,4] )

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1


! Dirichlet boundary conditions

! velocity on lower boundary
  call define_essential ( mesh, input_probdef, curve1=1, physq=2 )
! velocity on upper boundary
  call define_essential ( mesh, input_probdef, curve1=3, physq=2 )
! pressure in point 1
  call define_essential ( mesh, input_probdef, point=1, physq=3 )

! constraints for periodical boundary conditions

! velocities
  call define_constraint ( mesh, input_probdef, &
    physq=2, curve1=2, curve2=5, discretization='collocation', exclude=3 )

! gradients
  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, curve2=5, discretization='collocation' )

  call problem_definition ( input_probdef, mesh, problem )


! create a vector subscript for the vertical velocity for post processing

  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=2, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a =   &
      reshape ( [ ncompc,0,ncompc,0,ncompc,0,ncompc,0,0,    &  ! c
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar for plotting
                  [9,2] )

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

! constraint for periodical boundary conditions of the conformation
  call define_constraint ( mesh, input_probdefc, curve1=2, curve2=5, &
    discretization='collocation' )

  call problem_definition ( input_probdefc, mesh, problemc )


! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemc, cval, physqarr=[1] )
  call create_subscript ( mesh, problemc, cvalc2, physqarr=[1], curves=[2] )
  call create_subscript ( mesh, problemc, cvalc5, physqarr=[1], curves=[5] )
  call create_subscript ( mesh, problemc, cvalxx, physqarr=[1], degfd=1 )
  call create_subscript ( mesh, problemc, cvalxy, physqarr=[1], degfd=2 )
  call create_subscript ( mesh, problemc, cvalyy, physqarr=[1], degfd=3 )


! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, soln, solnm1, rhsd, solhat )


! fill solution vector with essential boundary conditions

! set velocity = 0 on lower boundary
  call fill_sysvector ( mesh, problem, soln, &
    curve1=1, physq=2, value=0._dp )
! set velocity = (1,0) on upper boundary
  call fill_sysvector ( mesh, problem, soln, &
    curve1=3, physq=2, degfd=1, value=1._dp )
  call fill_sysvector ( mesh, problem, soln, &
    curve1=3, physq=2, degfd=2, value=0._dp )
! set pressure level = 0 in lower left corner
  call fill_sysvector ( mesh, problem, soln, &
    point=1, physq=3, value=0._dp )


! create system matrix of gradient/velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! build (assemble) matrix and vector for gradient/velocity/pressure problem

! stokes velocity/pressure
  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, physqrow=[2,3], physqcol=[2,3], &
    coefficients=coefficients )

! DEVSS-G
  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=devssg_elem, addmatvec=.true., &
    physqrow=[1,2], physqcol=[1,2], coefficients=coefficients )

! set to zero off-diagonal blocks gradient-pressure
  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., physqrow=[1], physqcol=[3], zeromatvec=.true. )
  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., physqrow=[3], physqcol=[1], zeromatvec=.true. )

! periodical condition on velocities and gradients

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_constr_node_conn, addmatvec=.true. )

  call check ( sysmatrix )


! solve initial gradient/velocity/pressure consistent with initial conditions
! and keep decomposition of the matrix
! This generates a linear profile.

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, soln, rhsd )

  solver_options_u%real_storage=rs_gup
  solver_options_u%integer_storage=is_gup

  call solve_system_ma57 ( sysmatrix, rhsd, soln, lu_u, &
    solver_options=solver_options_u )


! create system vectors (solution and right-hand side) for conformation and
! initialize vectors with steady stress and a small random perturbation.

  call create ( problemc, solcn, solcnm1, solciter )
  call create ( problemc, dsolc, rhsc )

! steady state solution

  csteady(1) = 1 + 2 * lambda ** 2
  csteady(2) = lambda
  csteady(3) = 1

  call random_number (solcn(1)%u)
  if ( logc == 0 ) then ! standard
    solcn(1)%u(cvalxx%s) = csteady(1) + 1e-3*solcn(1)%u(cvalxx%s)! initial cxx
    solcn(1)%u(cvalxy%s) = csteady(2) + 1e-3*solcn(1)%u(cvalxy%s)! initial cxy
    solcn(1)%u(cvalyy%s) = csteady(3) + 1e-3*solcn(1)%u(cvalyy%s)! initial cyy
  else if ( logc == 1 ) then ! log scheme
    call eig2x2 ( csteady, lsteady, eigv )
    lsteady = log ( lsteady )
    call inveig2x2 ( csteady, lsteady, eigv )
    solcn(1)%u(cvalxx%s) = csteady(1) + 1e-3*solcn(1)%u(cvalxx%s)! initial sxx
    solcn(1)%u(cvalxy%s) = csteady(2) + 1e-3*solcn(1)%u(cvalxy%s)! initial sxy
    solcn(1)%u(cvalyy%s) = csteady(3) + 1e-3*solcn(1)%u(cvalyy%s)! initial syy
  end if
  solcn(1)%u(cvalc2%s) = solcn(1)%u(cvalc5%s) ! make disturbance periodic
  call copy(solcn,solciter)

! create system matrix for conformation problem

  call create_sysmatrix_structure_base ( sysmatrixc, mesh, problemc )
  call create_sysmatrix_structure_constraint ( sysmatrixc, mesh, problemc )
  call finalize_sysmatrix_structure ( sysmatrixc )


  call create_sysmatrix_data ( sysmatrixc )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec1=3, nprob=2 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => solhat
  oldvectors_ve%s(2)%p => solnm1
  oldvectors_ve%s1(1)%p => solcn
  oldvectors_ve%s1(2)%p => solcnm1
  oldvectors_ve%s1(3)%p => solciter
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc


! open monitor file

  open ( unit=13, file='out', recl=300 )
  open ( unit=14, file='iter.out', recl=300 )


! time stepping

  do step = 1, numtimesteps

    write(14,*) 'step = ', step
    if ( printscreen) print *, 'step = ', step

    if ( step == 2 ) then
!     change time integration scheme at the second time step
      coefficients%i(22) = timeint2
    end if

!   prediction of velocity and gradient for new time step
    if ( step == 1 ) then
 !    first order prediction
      solhat%u = soln%u
    else if ( step >= 2 .and. any ( timeint2 == [9,10,11] ) ) then
!     second order prediction
      solhat%u = 2*soln%u - solnm1%u
    end if

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do m = 1, nmodes

      coefficients%i(85) = m   ! set mode number

      if ( nmodes > 1 ) then
        write(14,*) 'mode = ', m
        if ( printscreen) print *, 'mode = ', m
      end if

      iter = 0

      do

        iter = iter + 1

        if ( iter > maxnumiterations ) then
          write(*,'(3(a,i0/))') &
            ' Maximum number of iterations reached = ', &
            maxnumiterations, ' step = ', step, 'mode = ', m
          stop
        end if

!       build (assemble) matrix and vector for conformation tensor problem

        call build_system ( mesh, problemc, sysmatrixc, sysvector=rhsc, &
          elemsub=implicit_ce_supg_elem, oldvectors=oldvectors_ve, &
          coefficients=coefficients )

!       periodical condition on conformation tensor
        call build_system_constraint ( mesh, problemc, sysmatrixc, &
          sysvector=rhsc, elemsub=stokes_constr_node_conn, &
          addmatvec=.true. )

        call check ( sysmatrixc )

!       solve c tensor

        call solve_system_ma41 ( sysmatrixc, rhsc, dsolc, &
          solver_options=solver_options_c  )

        solciter(m)%u(cval%s) = solciter(m)%u(cval%s) + dsolc%u(cval%s)

        if ( printscreen) print *, maxval(abs(dsolc%u(cval%s)))
        write(14,*) maxval(abs(dsolc%u(cval%s)))

        if ( maxval(abs(dsolc%u(cval%s))) < epsconf ) exit

      end do

    end do

!   copy previous c-tensor solution to older time step

    call copy ( solcn, solcnm1 )
    call copy ( solciter, solcn )


!   build (assemble) vector for gradient/velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_divtau, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      buildmatrix = .false., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, soln, rhsd )


!   copy previous gradient/velocity/pressure solution to older time step

    call copy ( soln, solnm1 )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma57 ( sysmatrix, rhsd, soln, lu_u, &
      solver_options=solver_options_u )


!   write monitor data

    write(unit=13,fmt=*) &
      step, step*deltat, sqrt( &
      sum( ( solcn(1)%u(cvalxx%s) - csteady(1) )**2 + &
           ( solcn(1)%u(cvalxy%s) - csteady(2) )**2 + &
           ( solcn(1)%u(cvalyy%s) - csteady(3) )**2 ) / 3 / size(cval%s) ), &
      minval(solcn(1)%u(cval%s)), &
      maxval(soln%u(vely%s)), minval(soln%u(vely%s))

  end do


! close monitor data file

  close(unit=13)
  close(unit=14)


! write data for post-processing

  call write_mesh ( mesh, filename='mesh.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh, input_probdefc, filename='probdefc.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) soln%u
  write(10) solcn(1)%u

  close(unit=10)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( soln, solnm1, rhsd, solhat )
  call delete ( sysmatrix )
  call delete ( oldvectors_ve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( lu_u )
  call delete ( solcn, solcnm1, solciter )
  call delete ( dsolc, rhsc )

  call delete ( coefficients )
  call delete ( vely, cvalxx, cvalxy, cvalyy, cval )
  call delete ( cvalc2, cvalc5 )

end program couette17
