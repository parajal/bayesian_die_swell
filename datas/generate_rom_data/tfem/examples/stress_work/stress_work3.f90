! Start-up of an Oldroyd-B fluid flow and temperature problem with stress work
! for the steady state.
! Periodic pipe flow. Constant coefficients.
! Error compared to analytical solution for velocity, conformation and
! temperature.

program stress_work3

  use tfem_m
  use math_defs_m
  use viscoelastic_elements_m
  use generalized_stokes_elements_m, only: fill_viscous_dissipation_gauss
  use convection_diffusion_elements_m
  use stress_work_functions_m
  use hsl_ma57_m
  use hsl_ma41_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,         & ! velocity interpolation
    pintpl = 2,         & ! pressure interpolation
    gintpl = 2,         & ! P1 gradients
    cintpl = 2,         & ! P1 conformation
    tintpl = uintpl,    & ! temperature interpolation
    physqg = 1,         & ! physical quantity nr of the gradients
    physqv = 2,         & ! physical quantity nr of the velocity
    physqp = 3,         & ! physical quantity nr of the pressure
    ncompc = 4,         & ! number of conformation tensor components
    nmodes = 1,         & ! number of modes
    startm = 501,       & ! start of material model data
    gauss = 6,          & ! Gauss integration on volumes
    gaussb = 3,         & ! Gauss integration on curves
    coorsys = 1,        & ! cylindrical coordinate system
    timeint = 1,        & ! (first-order) time integration
    numtimesteps = 2000, & ! number of time steps
    logc = 0,           & ! standard scheme or log transformation
    model = 2,          & ! 2: Oldroyd-B
    poselvec = 4,       & ! elvector position stress work (=3+max(1,numnlpar))
    nz = 2,             & ! number of elements in z-direction
    nr = 10               ! number of elements in r-direction

  real(dp), parameter :: &
    eta_s = 0.5_dp,    & ! solvent viscosity
    eta_p = 0.5_dp,    & ! polymer viscosity
    lambda = 1.0_dp,   & ! relaxation time
    deltat = 1.e-2_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    kappa = 5.0_dp,    & ! thermal conductivity coefficient
    Uavg = 1.0_dp,     & ! imposed averaged velocity
    T0 = 0.5_dp,       & ! imposed wall temperature
    lz = 0.6_dp,       & ! length of the domain in z-direction
    R = 3._dp,         & ! radius of the cylinder
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_t = 1.0_dp,     & ! real_storage temperature LU (HSL)
    is_t = 1.0_dp,     & ! integer_storage temperature LU (HSL)
    rs_c  = 1.5_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 1.5_dp       ! integer_storage for conformation LU (HSL)

! variables

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc, input_probdeft
  type(problem_t), target :: problem, problemc
  type(problem_t) :: problemt
  type(sysmatrix_t) :: sysmatrix, sysmatrixc, sysmatrixt
  type(sysvector_t), target :: sol
  type(sysvector_t) :: sol_exact, rhsd
  type(solver_options_ma57_t) :: solver_options_u, solver_optionst
  type(lu_ma57_t) :: lu_u
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc
  type(sysvector_t) :: solc_exact
  type(vector_t) :: tauviscous, tauviscoelastic
  type(lu_ma41_t) :: luc
  type(sysvector_t), target :: solt
  type(sysvector_t) :: solt_exact, rhsdt
  type(subscript_t) :: svel1, sT, cval, velx, vely
  type(coefficients_t) :: coefficients, coefficientst
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(solver_options_ma41_t) :: solver_options_c
  type(elvector_t), target :: sw_source
  type(elvector_t) :: sw_source_tmp

  integer :: icomp, step, i, grp, elem
  integer :: vertices(3) = [1,3,5]
  real(dp) :: alpha, G, flowrate


! Set function parameters

  Uavg_l = Uavg
  R_l = R
  eta_l = eta_s+eta_p
  lambda_l = lambda
  kappa_l = kappa
  T0_l = T0

! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus
  flowrate = Uavg * pi * R**2

! fill coefficients flow problem

  call create_coefficients ( coefficients, ncoefi=300, ncoefr=502 )

  coefficients%i = 0
  coefficients%i(1:23) = &
    [ uintpl, pintpl,  0,     0,      gintpl, &
      physqv, physqp,  0,     physqg, gauss,  &
      gaussb, cintpl,  0,     0,      0,      &
      0,      0,       model, nmodes, startm, &
      logc,   timeint, coorsys &
    ]

  coefficients%r = 0
  coefficients%r(1:9) = &
    [ eta_s,    0._dp,   0._dp, alpha, 0._dp, &
      flowrate, 0._dp,  deltat, beta &
    ]
  coefficients%r(501) = G
  coefficients%r(502) = lambda

! fill coefficients temperature problem

  call create_coefficients ( coefficientst, ncoefi=100, ncoefr=50 )

  coefficientst%i = 0
  coefficientst%i(1) = tintpl
  coefficientst%i(10:11) = [ gauss, gaussb ]
  coefficientst%i(18) = 3 ! right-hand side given per element in int. points
  coefficientst%i(23)  = coorsys

  coefficientst%r(1) = kappa
  coefficientst%r(2:) = 0


! create mesh

  meshgen_options%elshape = 4 ! 6-node triangles
  meshgen_options%nx = nz
  meshgen_options%ny = nr
  meshgen_options%lx = lz
  meshgen_options%ly = R

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] ) ! curve 5

  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=3 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 4  ! gradients
  input_probdef%vec_elementdof(1)%a(:,2) = 2  ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5) = 4  ! symmetric tensor

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

! center line
  call define_essential ( mesh, input_probdef, curve1=1, physq=physqv, &
    degfd=[0,1], exclude=2 )
! wall
  call define_essential ( mesh, input_probdef, curve1=3, physq=physqv )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=physqp )

! constraint for flow rate
  call define_constraint ( mesh, input_probdef, physq=physqv, curve1=2, &
    nglobalc=1 )

! constraint for periodical boundary conditions velocities
  call define_constraint ( mesh, input_probdef, &
    physq=physqv, curve1=2, curve2=5, discretization='collocation', exclude=2 )

! constraint for periodical boundary conditions gradients
  call define_constraint ( mesh, input_probdef, &
    physq=physqg, curve1=2, curve2=5, discretization='collocation' )

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for the cross velocity for post processing

  call create_subscript ( mesh, problem, velx, physqarr=[physqv], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqv], degfd=2 )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=3, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a(:,1) = 0
  input_probdefc%vec_elementdof(1)%a(vertices,1) = 1  ! c
  input_probdefc%vec_elementdof(1)%a(:,2) = 1         ! scalar for plotting
  input_probdefc%vec_elementdof(1)%a(:,3) = 4         ! symmetric tensor

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

 ! constraint for periodical boundary conditions of the conformation
  call define_constraint ( mesh, input_probdefc, curve1=2, curve2=5, &
    discretization='collocation' )

 call problem_definition ( input_probdefc, mesh, problemc )

! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemc, cval, physqarr=[1] )



! problem definition for the temperature problem

  call create_input_probdef ( mesh, input_probdeft, nvec=2, nphysq=1 )

  input_probdeft%vec_elementdof(1)%a =   &
      reshape ( [1,1,1,1,1,1,   &  ! temperature
                 1,1,1,1,1,1 ], &  ! scalar in in all nodes
                   [6,2] )

  input_probdeft%physq = [1]
  input_probdeft%probnr = 3

  call define_essential ( mesh, input_probdeft, curve1=3 )

  call define_constraint ( mesh, input_probdeft, &
    curve1=2, curve2=5, discretization='collocation', exclude=2 )

  call problem_definition ( input_probdeft, mesh, problemt )


! create system vectors (solution and right-hand side)

  call create ( problem, sol, rhsd )
  call create ( problemt, solt, rhsdt )

! create vector defined per element for scalar values of right-hand

  call create ( mesh, sw_source, nreal1d=1 )

! oldvectors

  call create ( oldvectors, nsysvec=1, nelvec=2 )
  oldvectors%s(1)%p => sol
  oldvectors%e(2)%p => sw_source

! set solution vectors to zero (including ess bc)

  sol%u = 0
  solt%u = 0


! create system matrix flow problem

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

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr_curve, &
    addmatvec=.true., coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, constraint2=3, elemsub=stokes_constr_node_conn, &
    addmatvec=.true., coefficients=coefficients )

  call check ( sysmatrix )

! solve initial gradient/velocity/pressure consistent with initial conditions
! and keep decomposition of the matrix
! This generates a Stokes profile.

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options_u%real_storage=rs_gup
  solver_options_u%integer_storage=is_gup

  call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
    solver_options=solver_options_u  )


! compute viscous dissipation

  call loop_over_elements ( mesh, problem, &
    elemsub=fill_viscous_dissipation_gauss, &
    coefficients=coefficients, oldvectors=oldvectors )

  call copy ( sw_source, sw_source_tmp )


! create system vectors (solution and right-hand side) for conformation and
! initialize vectors to zero stress

  call create ( problemc, solc, rhsc )

  if ( logc == 0 ) then ! standard
    solc(1,1)%u = 1
    solc(2,1)%u = 0
    solc(3,1)%u = 0
    solc(4,1)%u = 1
  else if ( logc == 1 ) then ! log scheme
    do i = 1, ncompc
      solc(i,1)%u = 0
    end do
  end if

! create system matrix for conformation problem

  call create_sysmatrix_structure_base ( sysmatrixc, mesh, problemc )
  call create_sysmatrix_structure_constraint ( sysmatrixc, mesh, problemc )
  call finalize_sysmatrix_structure ( sysmatrixc )

  call create_sysmatrix_data ( sysmatrixc )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=1, nsysvec2=1, nprob=2, &
    nelvec=poselvec )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc
  oldvectors_ve%e(poselvec)%p => sw_source

! open monitor file

  open ( unit=13, file='out', recl=300 )


! time stepping

  do step = 1, numtimesteps


!   build (assemble) matrix and vector for conformation problem

    call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
      elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
      coefficients=coefficients )

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problemc, sysmatrixc, &
      m2sysvector=rhsc, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

    call check ( sysmatrixc )


!   solve conformation and keep LU decomposition in the loop over components

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do icomp = 1, ncompc
      call solve_system_ma41 ( sysmatrixc, rhsc(icomp,1), solc(icomp,1), luc, &
        solver_options=solver_options_c  )
    end do

    call delete ( luc )  ! remove LU decomposition and rebuild next time step


!   build (assemble) vector for gradient/velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_divtau, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      buildmatrix = .false., coefficients=coefficients )

!   imposed flow rate (right-hand side only)

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr_curve, addmatvec=.true., &
      buildmatrix = .false., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
      solver_options=solver_options_u  )


!   write monitor data

    write(unit=13,fmt=*) &
      step, step*deltat, maxval(solc(1,1)%u(cval%s)), &
      minval(solc(2,1)%u(cval%s)), &
      maxval(solc(3,1)%u(cval%s)), &
      maxval(solc(4,1)%u(cval%s)), &
      maxval(abs(sol%u(velx%s))), maxval(abs(sol%u(vely%s)))

  end do


! close monitor data file

  close(unit=13)


! compute viscoelastic stress work

  call loop_over_elements ( mesh, problemc, &
    elemsub=fill_viscoelastic_stress_work_gauss, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

! add the two source terms

  do grp = 1, mesh%nelgrp
    do elem = 1, mesh%grpnumel(grp)
      sw_source%g(grp)%r1(elem,1)%a = &
          sw_source%g(grp)%r1(elem,1)%a + sw_source_tmp%g(grp)%r1(elem,1)%a
    end do
  end do


! create system matrix temperature problem

  call create_sysmatrix_structure_base ( sysmatrixt, mesh, problemt, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrixt, mesh, problemt )
  call finalize_sysmatrix_structure ( sysmatrixt )

  call create_sysmatrix_data ( sysmatrixt )

! solve temperature field

  call fill_sysvector ( mesh, problemt, solt, curve1=3, value=T0 )

  call build_system ( mesh, problemt, sysmatrixt, rhsdt, &
    elemsub=scalar_diffusion_elem, coefficients=coefficientst, &
    oldvectors=oldvectors )

  call build_system_constraint ( mesh, problemt, sysmatrixt, rhsdt, &
    constraint1=1, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficientst  )

  call add_effect_of_essential_to_rhs ( problemt, sysmatrixt, solt, rhsdt )

  solver_optionst%real_storage=rs_t
  solver_optionst%integer_storage=is_t

  call solve_system_ma57 ( sysmatrixt, rhsdt, solt, &
    solver_options=solver_optionst )


! postprocessing


! create exact solution for velocity

  call create_sysvector ( problem, sol_exact )

  call fill_sysvector ( mesh, problem, sol_exact, &
    node1=1, node2=mesh%nnodes, physq=physqv, degfd=1, func=func, funcnr=3 )

  call create_subscript ( mesh, problem, svel1, physqarr=[physqv], degfd=1 )

! print maximum difference of sol-solexact to standard output

  print *, 'max u_z error', maxval( abs(sol%u(svel1%s)-sol_exact%u(svel1%s)) )


! create exact solution for c

  call create_sysvector ( problemc, solc_exact )

! print maximum difference of solc-solc_exact to standard output

  call fill_sysvector ( mesh, problemc, solc_exact, &
    node1=1, node2=mesh%nnodes, physq=1, func=func, funcnr=8 )
  print *, 'max czz error', &
                      maxval( abs(solc(1,1)%u(cval%s)-solc_exact%u(cval%s)) )
  call fill_sysvector ( mesh, problemc, solc_exact, &
    node1=1, node2=mesh%nnodes, physq=1, func=func, funcnr=9 )
  print *, 'max czr error', &
                      maxval( abs(solc(2,1)%u(cval%s)-solc_exact%u(cval%s)) )


! create exact solution for temperature

  call create_sysvector ( problemt, solt_exact )

  call fill_sysvector ( mesh, problemt, solt_exact, &
    node1=1, node2=mesh%nnodes, physq=1, func=func, funcnr=4 )

  call create_subscript ( mesh, problemt, st, physqarr=[1] )

! print maximum difference of solt-solt_exact to standard output

  print *, 'max T error', maxval( abs(solt%u(st%s)-solt_exact%u(st%s)) )


! write vtk

  call create_vector ( problem, tauviscous, vec=5 )
  call create_vector ( problemc, tauviscoelastic, vec=3 )

  call derive_vector ( mesh, problem, tauviscous, elemsub=stokes_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problemc, tauviscoelastic, &
    elemsub=deriv_viscoelastic_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_vector_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='velocity', sysvector=sol, physq=physqv )

  call write_tensor_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='tauviscous', vector=tauviscous, append=.true. )

  call write_tensor_vtk ( mesh, problemc, filename='sol.vtk', &
    dataname='tauviscoelastic', vector=tauviscoelastic, append=.true. )

  call write_scalar_vtk ( mesh, problemt, filename='sol.vtk', &
    dataname='temperature', sysvector=solt, append=.true. )


! delete all data including all allocated memory

  call delete ( problem, problemt )
  call delete ( input_probdef, input_probdeft )
  call delete ( mesh )
  call delete ( sol, solt, sol_exact, solt_exact )
  call delete ( rhsd, rhsdt )
  call delete ( sw_source, sw_source_tmp )
  call delete ( sysmatrix, sysmatrixt )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( svel1, sT, cval, velx, vely )
  call delete ( tauviscous, tauviscoelastic )


end program stress_work3
