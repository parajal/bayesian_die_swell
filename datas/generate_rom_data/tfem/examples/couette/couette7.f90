! Viscoelastic problem on a unit square
! linear Couette flow
! DEVSS/DG
! first-order time integration

program couette7

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma57_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,     &  ! Q2 velocities
    pintpl = 2,     &  ! P1 pressures
    eintpl = 4,     &  ! Q1 gradients
    cintpl = 4,     &  ! Q1 conformation
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    gauss = 3,      & ! 3x3 integration of quads
    nintbc = 2,     & ! 2 point Gauss integration of DG boundary integrals
    ncompc = 3,     & ! number of conformation tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 2         ! Oldroyd-B


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc
  type(problem_t), target :: problem, problemc
  type(sysmatrix_t) :: sysmatrix, sysmatrixc
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vely
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc
  type(sysvector_t) :: cdot
  type(lu_ma57_t) :: luc, lu_u
  type(solver_options_ma57_t) :: solver_options_u, solver_options_c



! variables

  integer :: &
    nx=20,               & ! number of elements in x
    ny=20,               & ! number of elements in y
    DGtype = 1,          & ! explicit DG
    numtimesteps = 2,    & ! number of time steps
    logc = 0               ! standard scheme or log transformation

  real(dp) :: &
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 1.0_dp,   & ! relaxation time
    deltat = 5.e-3_dp, & ! time step
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.0_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 1.4_dp       ! integer_storage for conformation LU (HSL)

  logical :: buildmatrix
  integer :: icomp, step, i
  real(dp) :: alpha, G, csteady(3), lsteady(2), eigv(2,2)


! namelist for input of variables; read from standard input

  namelist /comppar/ nx, ny, numtimesteps, logc, eta_s, eta_p, &
    lambda, deltat, rs_gup, is_gup, rs_c, is_c

  read ( unit=*, nml=comppar )


! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+2*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     eintpl,         0, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gauss,    cintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      logc,     0,          0,     nintbc,    DGtype, &
      ( 0, i = 26, 150 )  &
    ]

  coefficients%r = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
      0._dp, 0._dp,  deltat,  0._dp,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G,     lambda &
    ]


! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call glue_mesh ( mesh, curve1=2, curve2=5 )

  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=3 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [ 3,0,3,0,3,0,3,0,0,    &  ! E
                  2,2,2,2,2,2,2,2,2,    &  ! velocity
                  0,0,0,0,0,0,0,0,3,    &  ! pressure
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                  [9,4] )

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

! Dirichlet boundary conditions

! velocity on lower boundary
  call define_essential ( mesh, input_probdef, curve1=1, physq=2 )
! velocity on upper boundary
  call define_essential ( mesh, input_probdef, curve1=3, physq=2 )
! pressure in element 1, center node, first degree
  call define_essential ( mesh, input_probdef, element=1, elnode=9, &
    physq=3, degfd=[1] )

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
      reshape ( [ 0,0,0,0,0,0,0,0,4,    &  ! c
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar for plotting
                  [9,2] )

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

  call problem_definition ( input_probdefc, mesh, problemc )


! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )


! fill solution vector with essential boundary conditions

  sol%u = 0

! set velocity_x = 1 on upper boundary
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=2, degfd=1, value=1._dp )


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
    elemsub=devss_elem, addmatvec=.true., &
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

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options_u%real_storage=rs_gup
  solver_options_u%integer_storage=is_gup

  call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
    solver_options=solver_options_u )


! create system vectors (solution and right-hand side) for conformation and
! initialize vectors with steady stress and a small random perturbation.

  call create ( problemc, cdot )
  call create ( problemc, solc, rhsc )

! steady state solution

  csteady(1) = 1 + 2 * lambda ** 2
  csteady(2) = lambda
  csteady(3) = 1

  if ( logc == 0 ) then ! standard
    call random_number (solc(1,1)%u)
    solc(1,1)%u = csteady(1) + 1e-3*solc(1,1)%u  ! initial cxx
    call random_number (solc(2,1)%u)
    solc(2,1)%u = csteady(2) + 1e-3*solc(2,1)%u  ! initial cxy
    call random_number (solc(3,1)%u)
    solc(3,1)%u = csteady(3) + 1e-3*solc(3,1)%u   ! initial cyy
  else if ( logc == 1 ) then ! log scheme
    call eig2x2 ( csteady, lsteady, eigv )
    lsteady = log ( lsteady )
    call inveig2x2 ( csteady, lsteady, eigv )
    call random_number (solc(1,1)%u)
    solc(1,1)%u = csteady(1) + 1e-3*solc(1,1)%u   ! initial sxx
    call random_number (solc(2,1)%u)
    solc(2,1)%u = csteady(2) + 1e-3*solc(2,1)%u   ! initial sxy
    call random_number (solc(3,1)%u)
    solc(3,1)%u = csteady(3) + 1e-3*solc(3,1)%u   ! initial syy
  end if

! create system matrix for conformation problem

  call create_sysmatrix_structure ( sysmatrixc, mesh, problemc, &
    symmetric=.true. )

  call create_sysmatrix_data ( sysmatrixc )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=1, nsysvec2=1, nprob=2 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc


! open monitor file

  open ( unit=13, file='out', recl=300 )

  buildmatrix = .true.

! time stepping

  do step = 1, numtimesteps


!   build (assemble) matrix and vector for conformation problem

    if ( step >= 2 ) buildmatrix = .false.

    call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
      elemsub=ce_dg_elem, oldvectors=oldvectors_ve, &
      coefficients=coefficients, buildmatrix=buildmatrix )


    call check ( sysmatrixc )

!   explicit Euler:
!
!       M ( c^n+1 - c^n ) / deltat = rhsd(c^n,t^n)
!   or
!       c^n+1 = c^n + deltat * M^-1 * rhsd(c^n,t^n)

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do icomp = 1, ncompc

!     solve conformation and keep LU decomposition for all timesteps

      call solve_system_ma57 ( sysmatrixc, rhsc(icomp,1), cdot, luc, &
        solver_options=solver_options_c )

      solc(icomp,1)%u = solc(icomp,1)%u + deltat * cdot%u

    end do


!   build (assemble) vector for gradient/velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_divtau, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      buildmatrix = .false., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
      solver_options=solver_options_u )


!   write monitor data

    write(unit=13,fmt=*) &
      step, step*deltat, sqrt( &
      sum( ( solc(1,1)%u - csteady(1) )**2 + &
           ( solc(2,1)%u - csteady(2) )**2 + &
           ( solc(3,1)%u - csteady(3) )**2 ) / 3 / size(solc(1,1)%u) ), &
      minval(solc(1,1)%u), &
      maxval(sol%u(vely%s)), minval(sol%u(vely%s))

  end do


! close monitor data file

  close(unit=13)


! write data for post-processing

  call write_mesh ( mesh, filename='mesh.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh, input_probdefc, filename='probdefc.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol%u
  write(10) (solc(i,1)%u, i=1,3)

  close(unit=10)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( oldvectors_ve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( lu_u, luc )
  call delete ( cdot )
  call delete ( solc, rhsc )

  call delete ( coefficients )
  call delete ( vely )

end program couette7
