! Axisymmetrical viscoelastic problem: a sphere moving with velocity U on the
! center line of a cylindrical container with a bottom (at infinity).
! The sphere is fixed and the cylinder wall is moved with velocity -U
! (moving frame).
! At the ends fluid is moving with velocity -U as well, representing zero
! flux in the stationary frame (bottom). At the right inflow wall the stress is
! assumed to be zero.
! Dirichlet boundary conditions for the velocities are used.
! DEVSS-G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! second-order time integration of the CE using a Gear/conformation
! time prediction scheme
! L2-projection of c=exp(s) on the finite element approximation space in the
! momentum balance.
! To increase the time step limit, c/lambda has been added to both sides of
! the discretized CE and this term in the LHS is treated implicitly, while
! it is treated explicitly on the RHS
! Giesekus model with four viscoelastic modes:

program sphere24

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m

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
    gauss = 3,      & ! 3x3 integration of quads
    ncompc = 4,     & ! number of conformation tensor components
    nmodes = 4,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 3,      & ! Giesekus
    coorsys = 1       ! axisymmetric coordinate system

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc, input_probdefc_proj
  type(problem_t), target :: problem, problemc, problemc_proj
  type(sysmatrix_t) :: sysmatrix, sysmatrixc, sysmatrixc_proj
  type(sysvector_t), target :: sol_np1
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: velr
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc_np1, solc_n, &
    solc_nm1, solc_proj
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc, rhsc_proj
  type(lu_ma41_t) :: luc
  type(lu_ma57_t) :: lu_exps_proj
  type(solver_options_ma41_t) :: solver_options_u, solver_options_c


! variables

  integer :: &
    timeint1 = 1,        & ! first-order, semi-implicit Euler
    timeint2 = 7,        & ! second order, semi-implicit Gear
    numtimesteps = 200,  & ! number of time steps
    restart = 0,         & ! do a restart
    logc = 1               ! standard scheme or log transformation

  real(dp) :: &
    eta_s = 0.0_dp,    & ! solvent viscosity
    eta_p(nmodes) = 1.0_dp,& ! polymer viscosities
    lambda(nmodes) = 1.0_dp , & ! relaxation times
    alphaG(nmodes) = 0.0_dp , & ! alpha parameter Giesekus
    deltat = 1.e-3_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.0_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 1.5_dp,    & ! integer_storage for conformation LU (HSL)
    U = 1._dp            ! velocity of the sphere

  logical :: printstepvel = .false.
  integer :: icomp, step, i, m
  real(dp) :: alpha, G(nmodes), velr_max_prev


! namelist for input of variables; read from standard input

  namelist /comppar/ timeint1, timeint2, numtimesteps, logc, eta_s, G, &
    lambda, alphaG, deltat, beta, rs_gup, is_gup, rs_c, is_c, U, restart

  read ( unit=*, nml=comppar )


! set some parameters

  eta_p = G * lambda
  alpha = sum(eta_p)  ! DEVSS parameter


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+3*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gauss,    cintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      logc,     timeint1,   coorsys, ( 0, i = 24, 150 )  &
    ]
  coefficients%i(49) = 1  ! exp(s) projection =.true. for logc=1
  coefficients%i(57) = 1  ! add and subtract c/lambda

  coefficients%r = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
      0._dp, 0._dp,  deltat,   beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      (G(i), lambda(i), alphaG(i), i=1,nmodes)  &
    ]


! read mesh

  call read_mesh_gmsh ( mesh, filename='mesh1.msh', ndim=2 )


  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=3 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [ 4,0,4,0,4,0,4,0,0,    &  ! G
                  2,2,2,2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,1,0,0,    &  ! pressure
                  1,1,1,1,1,1,1,1,1,    &  ! scalar, such as vorticity
                  4,4,4,4,4,4,4,4,4 ], &  ! tensor
                  [9,5] )

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1


! define essential boundaries

! center line
  call define_essential ( mesh, input_probdef, curve1=1, curve2=2, physq=2, &
    degfd=[0,1] )
  call define_essential ( mesh, input_probdef, curve1=7, curve2=8, physq=2, &
    degfd=[0,1] )
! sphere
  call define_essential ( mesh, input_probdef, curve1=3, curve2=6, physq=2 )
! ends and cylinder wall
  call define_essential ( mesh, input_probdef, curve1=9, curve2=14, physq=2 )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=3 )

  call problem_definition ( input_probdef, mesh, problem )


! create a vector subscript for the vertical velocity for post processing

  call create_subscript ( mesh, problem, velr, physqarr=[physqvel], degfd=2 )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=3, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a =   &
      reshape ( [ 1,0,1,0,1,0,1,0,0,    &  ! c
                  1,1,1,1,1,1,1,1,1,    &  ! scalar for plotting
                  4,4,4,4,4,4,4,4,4 ], &  ! tensor
                  [9,3] )

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

! essential bc for conformation at the right inflow boundary
  call define_essential ( mesh, input_probdefc, curve1=9 )

  call problem_definition ( input_probdefc, mesh, problemc )


! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol_np1, rhsd )


! fill solution vector with essential boundary conditions

! set velocity = 0 on boundaries
  sol_np1%u = 0
  call fill_sysvector ( mesh, problem, sol_np1, &
    curve1=9, curve2=14, physq=2, degfd=1, value=-U )


! create system matrix of gradient/velocity/pressure problem

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )


! create system vectors (solution and right-hand side) for conformation and
! initialize vectors with zero stress

  call create ( problemc, solc_np1, solc_n, solc_nm1, rhsc )

  do m = 1, nmodes
    if ( logc == 0 ) then ! standard
      solc_np1(1,m)%u = 1   ! initial czz
      solc_np1(2,m)%u = 0   ! initial czr
      solc_np1(3,m)%u = 1   ! initial crr
      solc_np1(4,m)%u = 1   ! initial ctt
    else if ( logc == 1 ) then ! log scheme
      solc_np1(1,m)%u = 0   ! initial szz
      solc_np1(2,m)%u = 0   ! initial szr
      solc_np1(3,m)%u = 0   ! initial srr
      solc_np1(4,m)%u = 0   ! initial ctt
    end if
  end do


! problem definition for projected "c=exp(s)" of the log conformation s

  call create_input_probdef ( mesh, input_probdefc_proj, nvec=1, nphysq=1 )

  input_probdefc_proj%vec_elementdof(1)%a = &
      reshape ( [ 1,0,1,0,1,0,1,0,0 ], &
                   [9,1] )

  input_probdefc_proj%physq = [1]
  input_probdefc_proj%probnr = 3

  call problem_definition ( input_probdefc_proj, mesh, problemc_proj )

  call create ( problemc_proj, solc_proj, rhsc_proj )

  do m = 1, nmodes
    solc_proj(1,m)%u = 1   ! initial czz
    solc_proj(2,m)%u = 0   ! initial czr
    solc_proj(3,m)%u = 1   ! initial crr
    solc_proj(4,m)%u = 1   ! initial ctt
  end do


! create system matrix for conformation problem

  call create_sysmatrix_structure ( sysmatrixc, mesh, problemc )

  call create_sysmatrix_data ( sysmatrixc )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec2=3, nprob=3 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol_np1
  oldvectors_ve%s2(1)%p => solc_n
  oldvectors_ve%s2(2)%p => solc_nm1
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


! restart: read solution from file

  if ( restart == 1 ) then

    open ( unit=10, form='unformatted', file='data.out' )

    read(10) sol_np1%u
    do m = 1, nmodes
      read(10) (solc_np1(i,m)%u, i=1,ncompc)
    end do

    close(unit=10)

  end if


! copy old values

  call copy ( solc_np1, solc_n )


! time stepping

  do step = 1, numtimesteps

    print *,'**************************************'
    print *,'step ',step

    if ( step == 2 ) then
      coefficients%i(22) = timeint2
    end if

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG


!   exps projection (for log conformation)

    if ( logc == 1 ) call solve_exps_projection


!   build implicit terms

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      addmatvec=.true., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol_np1, rhsd )

!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol_np1, &
      solver_options=solver_options_u  )


!   loop over modes

    do m = 1, nmodes


!     pass the current mode to the element routine

      coefficients%i(26) = m  ! mode1
      coefficients%i(27) = m  ! mode2


!     build (assemble) matrix and vector for conformation problem

      if ( step == 1 .or. timeint2 == 1 ) then
        call build_system ( mesh, problemc, sysmatrixc, msysvector=rhsc(:,m), &
          elemsub=ce_supg_elem1, oldvectors=oldvectors_ve, &
          coefficients=coefficients )
      else
        call build_system ( mesh, problemc, sysmatrixc, msysvector=rhsc(:,m), &
          elemsub=ce_supg_elem_implicit_2nd_order, &
          oldvectors=oldvectors_ve, coefficients=coefficients )
      end if

      call check ( sysmatrixc )

      do icomp = 1, ncompc
        call add_effect_of_essential_to_rhs ( problemc, sysmatrixc, &
          solc_np1(icomp,m), rhsc(icomp,m) )
      end do


!     solve conformation and keep LU decomposition in the loop over components

      solver_options_c%real_storage=rs_c
      solver_options_c%integer_storage=is_c

      do icomp = 1, ncompc
        call solve_system_ma41 ( sysmatrixc, rhsc(icomp,m), solc_np1(icomp,m), &
          luc, solver_options=solver_options_c  )
      end do

      call delete ( luc )  ! remove LU decomposition and rebuild next mode

    end do


!   set the modes back to zero (to include all modes in the momentum balance)

    coefficients%i(26) = 0  ! mode1
    coefficients%i(27) = 0  ! mode2


!   print the difference in maximum r-velocity with the previous time step

    if ( printstepvel ) then

      if ( step > 1) &
        print *, 'delta velr = ', maxval(abs(sol_np1%u(velr%s))) - velr_max_prev

      velr_max_prev = maxval(abs(sol_np1%u(velr%s)))

    end if

!   copy old values of the conformation problem

    call copy ( solc_n, solc_nm1 )
    call copy ( solc_np1, solc_n )

  end do

! write data for post-processing

  call write_mesh ( mesh, filename='mesh.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh, input_probdefc, filename='probdefc.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol_np1%u
  do m = 1, nmodes
    write(10) (solc_np1(i,m)%u, i=1,ncompc)
  end do

  close(unit=10)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol_np1, rhsd )
  call delete ( sysmatrix )
  call delete ( oldvectors_ve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( solc_np1, solc_nm1, rhsc )
  call delete ( solc_n )

  call delete ( coefficients )
  call delete ( velr )

  call delete ( sysmatrixc_proj )

  call delete ( problemc_proj )
  call delete ( input_probdefc_proj )
  call delete ( solc_proj, rhsc_proj )
  call delete ( lu_exps_proj )


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


end program sphere24
