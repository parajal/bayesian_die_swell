! Viscoelastic problem on a rectangular (2D) domain with 3D velocities
! (axisymmetric).
! Cylindrical coordinates with circumferential velocity included (swirl).
! Periodical conditions and imposed flow rate in z-direction.
! DEVSS-G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! L2-projection of c=exp(s) on the finite element approximation space in the
! momentum balance.
! second-order time integration to steady state

program concentric_cylinders4

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
  type(input_probdef_t) :: input_probdef, input_probdefc, input_probdefc_proj
  type(problem_t), target :: problem, problemc, problemc_proj
  type(sysmatrix_t) :: sysmatrix, sysmatrixc, sysmatrixc_proj
  type(sysvector_t), target :: sol, solm1
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, gammadot, Dtensor, Ltensor
  type(vector_t) :: tauviscous, tauviscoelastic, ctensor
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: velx, vely, velz, cval
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc, solcm1, solc_proj
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc, rhsc_proj
  type(lu_ma41_t) :: luc
  type(plot_options_t) :: plot_options
  type(solver_options_ma41_t) :: solver_options_u, solver_options_c
  type(lu_ma57_t) :: lu_exps_proj

! variables

  integer :: &
    nz=1,                & ! number of elements in z
    nr=20,               & ! number of elements in r
    timeint1 = 1,        & ! first-order time integration, first time step
    timeint2 = 7,        & ! second-order time integration after first time step
    numtimesteps = 100,  & ! number of time steps
    logc = 0,            & ! standard scheme or log transformation
    model = 3              ! 2: Oldroyd-B 3: Giesekus

  real(dp) :: &
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 1.0_dp,   & ! relaxation time
    alphapar = 0.1_dp, & ! alpha parameter in the Giesekus model
    deltat = 1.e-1_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    U = 1._dp,         & ! imposed average velocity in z direction
    W = 1._dp,         & ! Velocity in theta direction of the (outer) wall.
    rs_gup = 1.4_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
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

  coefficients%i(23) = 1  ! coorsys = 1, cilindrical coordinates, axisymmetric
  coefficients%i(49) = 1  ! exp(s) projection =.true. for logc=1
  coefficients%i(61) = 0  ! projected G=0, direct velocity gradient=1 in CE
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

  call create_sysvector ( problem, sol, solm1, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=physqv, degfd=3, value=W )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqp, value=0._dp )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! create system vectors (solution and right-hand side) for conformation and
! initialize vectors to zero stress

  call create ( problemc, solc, solcm1, rhsc )

  if ( logc == 0 ) then ! standard
    solc(1,1)%u = 1
    solc(2,1)%u = 0
    solc(3,1)%u = 0
    solc(4,1)%u = 1
    solc(5,1)%u = 0
    solc(6,1)%u = 1
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


! problem definition for projected "c=exp(s)" of the log conformation s

  call create_input_probdef ( mesh, input_probdefc_proj, nvec=1, nphysq=1 )

  input_probdefc_proj%vec_elementdof(1)%a(:,1) = 0
  input_probdefc_proj%vec_elementdof(1)%a(vertices,1) = 1  ! c

  input_probdefc_proj%physq = [1]
  input_probdefc_proj%probnr = 3

  call problem_definition ( input_probdefc_proj, mesh, problemc_proj )

  call create ( problemc_proj, solc_proj, rhsc_proj )

! initialize solc_proj
  solc_proj(1,1)%u = 1
  solc_proj(2,1)%u = 0
  solc_proj(3,1)%u = 0
  solc_proj(4,1)%u = 1
  solc_proj(5,1)%u = 0
  solc_proj(6,1)%u = 1



! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec2=3, nprob=3 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s(2)%p => solm1
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%s2(2)%p => solcm1
  oldvectors_ve%s2(3)%p => solc_proj
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc
  oldvectors_ve%p(3)%p => problemc_proj

! create and build system matrix for projection problem
! NOTE matrix remains constant and needs to be build once.

  call create_sysmatrix_structure ( sysmatrixc_proj, mesh, problemc_proj, &
    symmetric=.true. )
  call create_sysmatrix_data ( sysmatrixc_proj )

  call build_system ( mesh, problemc_proj, sysmatrixc_proj, &
    m2sysvector=rhsc_proj, elemsub=exps_projection_elem, &
    oldvectors=oldvectors_ve, coefficients=coefficients, &
    buildvector=.false. )

  call check ( sysmatrixc_proj )


! open monitor file

  open ( unit=13, file='out', recl=300 )


! time stepping

  do step = 1, numtimesteps

    if ( step == 2 ) then
!     change time integration scheme at the second time step
      coefficients%i(22) = timeint2
    end if

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG


!   exps projection (for log conformation)

    if ( logc == 1 ) call solve_exps_projection


!   build implicit terms of CE with rhs in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      addmatvec=.true., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

    call copy ( sol, solm1 )


!   build (assemble) matrix and vector for conformation problem

    if ( coefficients%i(22) == timeint1 ) then

      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    else

      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem_implicit_2nd_order, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    end if

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problemc, sysmatrixc, &
      m2sysvector=rhsc, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

    call check ( sysmatrixc )

    call copy ( solc, solcm1 )


!   solve conformation and keep LU decomposition in the loop over components

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do icomp = 1, ncompc
      call solve_system_ma41 ( sysmatrixc, rhsc(icomp,1), solc(icomp,1), luc, &
        solver_options=solver_options_c  )
    end do

    call delete ( luc )  ! remove LU decomposition and rebuild next time step


!   write monitor data

    write(unit=13,fmt=*) &
      step, step*deltat, maxval(solc(1,1)%u(cval%s)), &
      maxval(abs(sol%u(velx%s))), maxval(abs(sol%u(vely%s))), &
      maxval(abs(sol%u(velz%s)))

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

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  oldvectors%s(1)%p => sol

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
    filename='concyl4.vtk' )

  call write_vector_vtk ( mesh, problem, filename='concyl4.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true., &
    assume3D=.true. )

  call write_vector_vtk ( mesh, problem, filename='concyl4.vtk', &
    dataname='velocity_vector_planar', vector=velocity, append=.true., &
    degfd=[1,2] )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='concyl4.vtk', dataname='gammadot', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='concyl4.vtk', &
    dataname='D', vector=Dtensor, append=.true., assume3D=.true. )

  call write_tensor_vtk ( mesh, problem, filename='concyl4.vtk', &
    dataname='L', vector=Ltensor, append=.true., assume3D=.true., &
    symmetric=.false. )

  call write_tensor_vtk ( mesh, problem, filename='concyl4.vtk', &
    dataname='tauviscous', vector=tauviscous, append=.true., assume3D=.true. )

  call create_vector ( problemc, tauviscoelastic, vec=3 )

  call derive_vector ( mesh, problemc, tauviscoelastic, &
    elemsub=deriv_viscoelastic_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemc, filename='concyl4.vtk', &
    dataname='tauviscoelastic', vector=tauviscoelastic, append=.true., &
    assume3D=.true. )

  call create_vector ( problemc, ctensor, vec=3 )

  call derive_vector ( mesh, problemc, ctensor, &
    elemsub=deriv_conformation_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemc, filename='concyl4.vtk', &
    dataname='c', vector=ctensor, append=.true., &
    assume3D=.true. )

! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, solm1, rhsd )
  call delete ( sysmatrix )
  call delete ( velocity, pressure, gammadot, Dtensor, Ltensor, tauviscous, &
                tauviscoelastic, ctensor )
  call delete ( oldvectors, oldvectors_ve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( solc, solcm1, rhsc )

  call delete ( problemc_proj )
  call delete ( input_probdefc_proj )
  call delete ( solc_proj, rhsc_proj )
  if ( logc == 1 ) call delete ( lu_exps_proj )

  call delete ( coefficients )
  call delete ( velx, vely, velz, cval )

contains

  subroutine build_vpG

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, physqrow=[2,3], physqcol=[2,3], &
      coefficients=coefficients )

!   DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_elem, addmatvec=.true., &
      physqrow=[1,2], physqcol=[1,2], coefficients=coefficients )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[1], physqcol=[3], zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[3], physqcol=[1], zeromatvec=.true. )

!   imposed flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr_curve, addmatvec=.true., &
      coefficients=coefficients )

!   periodical boundaries

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, constraint2=3, elemsub=stokes_constr_node_conn, &
      addmatvec=.true., coefficients=coefficients )

  end subroutine build_vpG


! project c=exp(s) on discrete fem space

  subroutine solve_exps_projection

    type(solver_options_ma57_t) :: solver_options_ma57

    integer :: i, m

!   build vector only (matrix is constant)

    call build_system ( mesh, problemc_proj, sysmatrixc_proj, &
      m2sysvector=rhsc_proj, elemsub=exps_projection_elem, &
      oldvectors=oldvectors_ve, coefficients=coefficients, &
      buildmatrix=.false. )

    ! MA57 solver storage
    solver_options_ma57%integer_storage = 1.3
    solver_options_ma57%real_storage    = 1.3

!   LU decomposition is done in the first call only

    do m = 1, nmodes
      do i = 1, ncompc

        call add_effect_of_essential_to_rhs ( problemc_proj, sysmatrixc_proj, &
           solc_proj(i,m), rhsc_proj(i,m) )

        call solve_system_ma57 ( sysmatrixc_proj, rhsc_proj(i,m), &
           solc_proj(i,m), lu_exps_proj, solver_options=solver_options_ma57 )

      end do
    end do

  end subroutine solve_exps_projection

end program concentric_cylinders4
