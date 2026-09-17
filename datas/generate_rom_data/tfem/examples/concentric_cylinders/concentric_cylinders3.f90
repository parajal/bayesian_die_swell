! Viscoelastic problem on a rectangular (2D) domain with 3D velocities
! (axisymmetric).
! Cylindrical coordinates with circumferential velocity included (swirl).
! Periodical conditions and imposed flow rate in z-direction.
! DEVSS-G/SUPG
! second-order time integration to steady state

program concentric_cylinders3

  use tfem_m
  use math_defs_m
  use hsl_ma41_m
  use hsl_ma57_m
  use viscoelastic_elements_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  logical, parameter :: concentric = .true.  ! F: Single pipe
                                             ! T: Two concentric walls

  integer, parameter :: &
    uintpl = 8,  & ! Q2 velocities
    pintpl = 4,  & ! Q1 pressures
    gintpl = 4,  & ! Q1 gradients
    cintpl = 4,  & ! Q1 conformation
    physqg = 1,  & ! physical quantity nr of the gradients
    physqv = 2,  & ! physical quantity nr of the velocities
    physqp = 3,  & ! physical quantity nr of the pressures
    gauss  = 3,  & ! 3x3 integration of quads
    ncompc = 6,  & ! number of conformation tensor components
    nmodes = 1,  & ! number of modes
    startm = 501   ! start of material model data

  real(dp), parameter :: &
    R1 = 0.5_dp, & ! inner radius (for concentric=.true.)
    R2 = 1._dp,  & ! outer radius
    er = 10._dp    ! ratio of element length z/r.

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc
  type(problem_t), target :: problem, problemc
  type(sysmatrix_t) :: sysmatrix, sysmatrixc
  type(sysvector_t), target :: soln, solnm1
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, gammadot, Dtensor, Ltensor
  type(vector_t) :: tauviscous, tauviscoelastic, ctensor
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: velx, vely, velz, cval
  type(sysvector_t), dimension(ncompc,nmodes), target :: solcn, solcnm1
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc
  type(lu_ma57_t) :: lu_u
  type(lu_ma41_t) :: luc
  type(plot_options_t) :: plot_options
  type(solver_options_ma57_t) :: solver_options_u
  type(solver_options_ma41_t) :: solver_options_c

! variables

  integer :: &
    nz=1,                & ! number of elements in z
    nr=20,               & ! number of elements in r
    timeint1 = 1,        & ! first-order time integration, first time step
    timeint2 = 2,        & ! second-order time integration after first time step
    numtimesteps = 2000, & ! number of time steps
    logc = 0,            & ! standard scheme or log transformation
    model = 3              ! 2: Oldroyd-B 3: Giesekus

  real(dp) :: &
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 1.0_dp,   & ! relaxation time
    alphapar = 0.1_dp, & ! alpha parameter in the Giesekus model
    deltat = 5.e-3_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    U = 1._dp,         & ! imposed average velocity in z direction
    W = 1._dp,         & ! Velocity in theta direction of the (outer) wall.
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.4_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 1.4_dp       ! integer_storage for conformation LU (HSL)

  integer :: icomp, step, i, npar
  integer :: vertices(4) = [1,3,5,7]
  real(dp) :: alpha, G, flowrate

! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus

  if ( model == 2 ) then
    npar = 2 ! Oldroyd-B
  else if ( model == 3 ) then
    npar = 3 ! Giesekus
  end if

  if ( concentric ) then
    flowrate = U * pi * ( R2 ** 2 - R1 ** 2 )
  else
    flowrate = U * pi * R2 ** 2
  end if

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+npar*nmodes )

  coefficients%i = &
    [ uintpl, pintpl, 0, 0,      gintpl, &
      physqv, physqp, 0, physqg, gauss,  &
      gauss,  cintpl, 0, 0,      0,      &
      0,      0,      model, nmodes, startm, &
      logc,   timeint1, ( 0, i = 23, 150 )  &
    ]

  coefficients%i(23) = 1  ! cylindrical coordinate system (axisymmetrical)
  coefficients%i(67) = 1  ! 3D velocity (swirl)

  coefficients%r(1:502) = &
    [ eta_s,    0._dp,   0._dp, alpha, 0._dp, &
      flowrate, 0._dp,  deltat, beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G, lambda  &
    ]

  if ( model == 3 ) then
    coefficients%r(503) = alphapar ! Giesekus
  end if

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nz
  meshgen_options%ny = nr
  if ( concentric ) then
    meshgen_options%oy = R1
    meshgen_options%ly = R2 - R1
  else
    meshgen_options%ly = R2
  end if
  meshgen_options%lx = er * nz * meshgen_options%ly / nr

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=6, nphysq=3 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 6  ! gradients
  input_probdef%vec_elementdof(1)%a(:,2) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5) = 6  ! symmetric tensor
  input_probdef%vec_elementdof(1)%a(:,6) = 9  ! unsymmetric tensor

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

  if ( concentric ) then
!   inner cilinder wall
    call define_essential ( mesh, input_probdef, curve1=1, physq=physqv )
  else
!   center line
    call define_essential ( mesh, input_probdef, curve1=1, physq=physqv, &
      degfd=[0,1,1], exclude=1 )
  end if
! outer wall
  call define_essential ( mesh, input_probdef, curve1=3, physq=physqv )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=physqp )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqv, curve1=2, nglobalc=1 )

! constraint for periodical boundary conditions

! velocities
  if ( concentric ) then
    call define_constraint ( mesh, input_probdef, &
      physq=physqv, curve1=2, curve2=5, discretization='collocation', exclude=3 )
  else
    call define_constraint ( mesh, input_probdef, &
      physq=physqv, curve1=2, curve2=5, discretization='collocation', exclude=2 )
  end if

! gradients
  call define_constraint ( mesh, input_probdef, &
    physq=physqg, curve1=2, curve2=5, discretization='collocation' )

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for the cross velocity for post processing

  call create_subscript ( mesh, problem, velx, physqarr=[physqv], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqv], degfd=2 )
  call create_subscript ( mesh, problem, velz, physqarr=[physqv], degfd=3 )

! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=3, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a(:,1) = 0
  input_probdefc%vec_elementdof(1)%a(vertices,1) = 1  ! c
  input_probdefc%vec_elementdof(1)%a(:,2) = 1         ! scalar for plotting
  input_probdefc%vec_elementdof(1)%a(:,3) = 6         ! symmetric tensor

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

 ! constraint for periodical boundary conditions of the conformation
  call define_constraint ( mesh, input_probdefc, curve1=2, curve2=5, &
    discretization='collocation' )

 call problem_definition ( input_probdefc, mesh, problemc )

! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemc, cval, physqarr=[1] )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, soln, solnm1, rhsd )

! fill solution vector with essential boundary conditions

  soln%u = 0
  call fill_sysvector ( mesh, problem, soln, &
    curve1=3, physq=physqv, degfd=3, value=W )
  call fill_sysvector ( mesh, problem, soln, &
    point=1, physq=physqp, value=0._dp )

! create system matrix

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

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, soln, rhsd )

  solver_options_u%real_storage=rs_gup
  solver_options_u%integer_storage=is_gup

  call solve_system_ma57 ( sysmatrix, rhsd, soln, lu_u, &
    solver_options=solver_options_u  )

! create system vectors (solution and right-hand side) for conformation and
! initialize vectors to zero stress

  call create ( problemc, solcn, solcnm1, rhsc )

  if ( logc == 0 ) then ! standard
    solcn(1,1)%u = 1
    solcn(2,1)%u = 0
    solcn(3,1)%u = 0
    solcn(4,1)%u = 1
    solcn(5,1)%u = 0
    solcn(6,1)%u = 1
  else if ( logc == 1 ) then ! log scheme
    do i = 1, ncompc
      solcn(i,1)%u = 0
    end do
  end if

! create system matrix for conformation problem

  call create_sysmatrix_structure_base ( sysmatrixc, mesh, problemc )
  call create_sysmatrix_structure_constraint ( sysmatrixc, mesh, problemc )
  call finalize_sysmatrix_structure ( sysmatrixc )

  call create_sysmatrix_data ( sysmatrixc )



! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec2=2, nprob=2 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => soln
  oldvectors_ve%s(2)%p => solnm1
  oldvectors_ve%s2(1)%p => solcn
  oldvectors_ve%s2(2)%p => solcnm1
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc

! open monitor file

  open ( unit=13, file='out', recl=300 )


! time stepping

  do step = 1, numtimesteps

    if ( step == 2 ) then
!     change time integration scheme at the second time step
      coefficients%i(22) = timeint2
    end if

!   build (assemble) matrix and vector for conformation problem

    call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
      elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
      coefficients=coefficients )

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problemc, sysmatrixc, &
      m2sysvector=rhsc, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

    call check ( sysmatrixc )

!   copy previous conformation solution to older time step

    call copy ( solcn, solcnm1 )


!   solve conformation and keep LU decomposition in the loop over components

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do icomp = 1, ncompc
      call solve_system_ma41 ( sysmatrixc, rhsc(icomp,1), solcn(icomp,1), luc, &
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

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, soln, rhsd )


!   copy previous gradient/velocity/pressure solution to older time step

    call copy ( soln, solnm1 )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma57 ( sysmatrix, rhsd, soln, lu_u, &
      solver_options=solver_options_u  )


!   write monitor data

    write(unit=13,fmt=*) &
      step, step*deltat, maxval(solcn(1,1)%u(cval%s)), &
      maxval(abs(soln%u(velx%s))), maxval(abs(soln%u(vely%s))), &
      maxval(abs(soln%u(velz%s)))

  end do


! close monitor data file

  close(unit=13)

! post-processing

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  call create_vector ( problem, velocity, physq=physqv )
  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, gammadot, vec=4 )
  call create_vector ( problem, Dtensor, vec=5 )
  call create_vector ( problem, Ltensor, vec=6 )
  call create_vector ( problem, tauviscous, vec=5 )

  call extract_physvector ( mesh, problem, soln, velocity )

! create the structure oldvectors

  oldvectors%s(1)%p => soln

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=11

  call derive_vector ( mesh, problem, gammadot, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Dtensor, elemsub=stokes_D_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Ltensor, elemsub=stokes_gradu_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, tauviscous, elemsub=stokes_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='concyl3.vtk' )

  call write_vector_vtk ( mesh, problem, filename='concyl3.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true., &
    assume3D=.true. )

  call write_vector_vtk ( mesh, problem, filename='concyl3.vtk', &
    dataname='velocity_vector_planar', vector=velocity, append=.true., &
    degfd=[1,2] )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='concyl3.vtk', dataname='gammadot', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='concyl3.vtk', &
    dataname='D', vector=Dtensor, append=.true., assume3D=.true. )

  call write_tensor_vtk ( mesh, problem, filename='concyl3.vtk', &
    dataname='L', vector=Ltensor, append=.true., assume3D=.true., &
    symmetric=.false. )

  call write_tensor_vtk ( mesh, problem, filename='concyl3.vtk', &
    dataname='tauviscous', vector=tauviscous, append=.true., assume3D=.true. )

  call create_vector ( problemc, tauviscoelastic, vec=3 )

  call derive_vector ( mesh, problemc, tauviscoelastic, &
    elemsub=deriv_viscoelastic_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemc, filename='concyl3.vtk', &
    dataname='tauviscoelastic', vector=tauviscoelastic, append=.true., &
    assume3D=.true. )

  call create_vector ( problemc, ctensor, vec=3 )

  call derive_vector ( mesh, problemc, ctensor, &
    elemsub=deriv_conformation_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemc, filename='concyl3.vtk', &
    dataname='c', vector=ctensor, append=.true., &
    assume3D=.true. )

! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( soln, solnm1, rhsd )
  call delete ( sysmatrix )
  call delete ( velocity, pressure, gammadot, Dtensor, Ltensor, tauviscous, &
                tauviscoelastic, ctensor )
  call delete ( oldvectors, oldvectors_ve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( lu_u )
  call delete ( solcn, solcnm1, rhsc )

  call delete ( coefficients )
  call delete ( velx, vely, velz, cval )

end program concentric_cylinders3
