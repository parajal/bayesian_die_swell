! Flow around a cylinder confined between two plates for a FENE dumbbell,
! solved using BCF. The entry and exit sides are connected by periodical
! boundary conditions and a constant flow rate is imposed.
! Only the upper-half of the domain is modelled and symmetry conditions are
! assumed (zero tractions).
! Dirichlet boundary conditions for the velocities are used.
! DEVSS/DG
! first-order implicit time integration

program cylinder3

  use tfem_m
  use bcf_elements_m
  use bcf_functions_m
  use hsl_ma57_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,     &  ! Q2 velocities
    pintpl = 2,     &  ! P1 pressures
    eintpl = 4,     &  ! Q1 gradients
    qintpl = 4,     &  ! Q1 Q-vectors
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    gauss = 3,      & ! 3x3 integration of quads
    nintbc = 2,     & ! 2 point Gauss integration of DG boundary integrals
    ncompc = 4,     & ! number of structure tensor/stress tensor components
    ndfq = 4,       & ! nodal/Gauss points in the element
    ncompq = 3,     & ! number of Q-vector components
    startm = 501,   & ! start of material model data
    model = 3         ! FENE dumbbell

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefq
  type(problem_t), target :: problem, problemq
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(sysvector_t), allocatable, dimension(:), target :: solqn, solqnp1
  type(lu_ma57_t) :: lu_u
  type(subscript_t) :: vely
  type(vector_t) :: ctensor
  type(solver_options_ma57_t) :: solver_options


! variables

  integer :: &
    timeint = 2,         & ! first-order implicit
    eqnumtimesteps = 2,  & ! number of time steps for equilibration
    numtimesteps = 2,    & ! number of time steps
    printtimesteps = 5,  & ! print every .. timesteps
    nfield = 1000,       & ! number of fields
    restart = 0,         & ! do a restart of data
    logq = 0               ! standard scheme or log transformation

  real(dp) :: &
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 1.0_dp,   & ! relaxation time
    b = 50.0_dp,       & ! nonlinear parameter for FENE
    deltat = 5.e-3_dp, & ! time step
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    H = 2._dp,         & ! (half-)height of the channel
    U = 1._dp            ! average velocity in the channel


  integer :: step, i, k, j
  real(dp) :: alpha, modulus
  real(dp), allocatable, dimension(:,:) :: qn, qnp1
  real(dp), allocatable, dimension(:,:,:) :: work


! namelist for input of variables; read from standard input


  namelist /comppar/ timeint, eqnumtimesteps, numtimesteps, printtimesteps, &
    nfield, logq, eta_s, eta_p, lambda, deltat, rs_gup, is_gup, U, b, restart

  read ( unit=*, nml=comppar )


! allocate the arrays for Q-vectors (not data!)

  allocate ( solqn(nfield), solqnp1(nfield) )

! allocate the Brownian vector and Q-vectors for equilibration

  allocate ( brown_bcf(nfield,ncompq), qn(nfield,ncompq), qnp1(nfield,ncompq) )



! set some parameters

  if ( model  == 2 ) then   ! FENE-P

    modulus = eta_p / lambda * ( b + 3 ) / b

  else if ( model == 3 ) then  ! FENE

    modulus = eta_p / lambda * ( b + 5 ) / b

  else

    modulus = eta_p / lambda  ! modulus nkT

  end if

  alpha = eta_p  ! DEVSS parameter = polymer viscosity


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=500, ncoefr=500+3 )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     eintpl,         0, &
      physqvel, physqpress, 0,     physqgrad,  gauss, &
      gauss,    0,          0,     0,              0, &
        ( 0, i = 16, 450 ),                           &
      qintpl,   model, startm,     logq,      nfield, &
      timeint,  nintbc,     0,     0,              0, &
        ( 0, i = 461, 500 ) &
    ]

  coefficients%r = &
    [ eta_s,  0._dp,  0._dp,   alpha, 0._dp, &
      U*H,    0._dp,  0._dp,   0._dp, 0._dp, &
      ( 0._dp, i = 11, 400 ), &
      deltat, 0._dp,  0._dp,   0._dp, 0._dp, &
      ( 0._dp, i = 406, 500 ), &
      modulus, lambda, b  &
    ]

  coefficients%qinflow => qinflow
  coefficients%brownian_vector => brownian_vector

  call write_coefficients ( coefficients, filename='coefficients.out' )


! read mesh

  call read_mesh ( mesh, filename='mesh.out' )


! connect curves for periodical bc in DG

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 25

  call glue_mesh ( mesh, curve1=18, curve2=25 )


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


! define essential boundaries

! center line
  call define_essential ( mesh, input_probdef, curve1=23, curve2=24, physq=2, &
    degfd=[0,1] )
! cylinder
  call define_essential ( mesh, input_probdef, curve1=20, physq=2 )
! wall
  call define_essential ( mesh, input_probdef, curve1=21, physq=2 )
! pressure level
  call define_essential ( mesh, input_probdef, element=1, elnode=9, &
    physq=3, degfd=[1] )



! constraints for periodical boundary conditions

! flow rate
  call define_constraint ( mesh, input_probdef, &
    physq=2, curve1=25, nglobalc=1 )

! velocities (use weak connection)
  call define_constraint ( mesh, input_probdef, &
    physq=2, curve1=18, curve2=25, discretization='weak', elementdof=[2,0,2] )

! gradients (use collocation)
  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=18, curve2=25, discretization='collocation' )


  call problem_definition ( input_probdef, mesh, problem )


! create a vector subscript for the vertical velocity for post processing

  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )


! problem definition Q-vectors

  call create_input_probdef ( mesh, input_probdefq, nvec=3, nphysq=1 )

  input_probdefq%vec_elementdof(1)%a(:,1) = [ (0,i=1,8), ndfq*ncompq ] ! Q-vec
  input_probdefq%vec_elementdof(1)%a(:,2) = 1 ! scalar for plotting
  input_probdefq%vec_elementdof(1)%a(:,3) = ncompc ! tensor for plotting

  input_probdefq%physq = [1]
  input_probdefq%probnr = 2

  call problem_definition ( input_probdefq, mesh, problemq )


! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )


! fill solution vector with essential boundary conditions

  sol%u = 0


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

! DEVSS
  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=devss_elem, addmatvec=.true., &
    physqrow=[1,2], physqcol=[1,2], coefficients=coefficients )

! set to zero off-diagonal blocks gradient-pressure
  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., physqrow=[1], physqcol=[3], zeromatvec=.true. )
  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., physqrow=[3], physqcol=[1], zeromatvec=.true. )

! flow rate

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

! periodical condition on velocities and gradients

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_elem_conn, &
    addmatvec=.true., coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=3, elemsub=stokes_constr_node_conn, &
    addmatvec=.true. )

  call check ( sysmatrix )




! solve initial gradient/velocity/pressure consistent with initial conditions
! and keep decomposition of the matrix

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%real_storage=rs_gup
  solver_options%integer_storage=is_gup

  call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
    solver_options=solver_options )


! create system vectors Q
! initialize vectors with steady stress and a small random perturbation.

  call create ( problemq, solqn, solqnp1 )

! create steady state solution

  call equilibrate


! Fill initial Q-vector

  allocate ( work(ndfq,ncompq,mesh%nelem) )

  do i = 1, nfield
    do k = 1, mesh%nelem
      do j = 1, ndfq
         work(j,:,k) = qnp1(i,:)
      end do
    end do
    solqn(i)%u = reshape ( work, [ndfq*ncompq*mesh%nelem] )
  end do

  deallocate ( work )


! vector for structure tensor

  call create_vector ( problemq, ctensor, vec=3 )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=1, nsysvec1=2, nprob=2 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s1(1)%p => solqn
  oldvectors_ve%s1(2)%p => solqnp1
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemq

! restart: read solution from file

  if ( restart == 1 ) then

    open ( unit=10, form='unformatted', file='data_cont.out')

    read(10) sol%u
    read(10) (solqn(i)%u, i=1,nfield)

    close(unit=10)

  end if

! open monitor file

  open ( unit=13, file='out', recl=300 )


! time stepping

  do step = 1, numtimesteps


!   fill Brownian vector with uniform distribution, variance 1

    call random_number ( brown_bcf )

    brown_bcf = ( 2 * brown_bcf - 1 ) * sqrt(3._dp)

!   compute Q-vector at next time step (convection)

    call loop_over_elements ( mesh, problemq, elemsub=bcf_conv_dg_elem, &
      oldvectors=oldvectors_ve, coefficients=coefficients )

!   compute Q-vector at next time step (model)

    call loop_over_elements ( mesh, problemq, elemsub=bcf_model_dg_elem, &
      oldvectors=oldvectors_ve, coefficients=coefficients )


!   build (assemble) vector for gradient/velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=bcf_rhs_divtau, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      buildmatrix=.false., coefficients=coefficients )

!   flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      buildmatrix=.false., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve gradient/velocity/pressure problem

    solver_options%real_storage=rs_gup
    solver_options%integer_storage=is_gup

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
      solver_options=solver_options )


!   copy vector

    call copy ( solqnp1, solqn )

    if ( MOD(step,printtimesteps)==0 ) then

!     derive structure tensor

       call derive_vector ( mesh, problemq, ctensor, &
         elemsub=deriv_structure_tensor, &
         coefficients=coefficients, oldvectors=oldvectors_ve )

!     write monitor data

      write(unit=13,fmt=*) &
        step, step*deltat, maxval(sol%u(vely%s)), maxval(ctensor%u)

    end if

  end do


! close monitor data file

  close(unit=13)

! write data for post-processing

  call write_mesh ( mesh, filename='mesh.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh, input_probdefq, filename='probdefq.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol%u
  write(10) (solqn(i)%u, i=1,nfield)

  close(unit=10)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( oldvectors_ve )
  call delete ( problemq )
  call delete ( input_probdefq )
  call delete ( lu_u )
  call delete ( solqn, solqnp1 )

  call delete ( coefficients )
  call delete ( ctensor )
  call delete ( vely )

  deallocate ( solqn, solqnp1 )
  deallocate ( brown_bcf, qn, qnp1 )

contains

  subroutine equilibrate

    type(stmodel_t) :: stmodel
    type(stnumpar_t) :: stnumpar
    integer :: step
    real(dp) :: gradv(4)

!   define the model

    call create_stochastic_model ( model, stmodel )

!   set material parameters

    stmodel%modulus = modulus
    stmodel%lambda = lambda
    stmodel%nonlin = b

!   initialize velocity gradient

    gradv = 0

!   initialize q

    qn = 0

!   equilibrate

    stnumpar%nfield = nfield
    stnumpar%timeint = timeint
    stnumpar%tstep = deltat

    do step = 1, eqnumtimesteps

!     fill Brownian vector with uniform distribution, variance 1

      call random_number ( brown_bcf )

      brown_bcf = ( 2 * brown_bcf - 1 ) * sqrt(3._dp)

!     step Q

      call ststep1_2D ( stmodel, stnumpar, gradv, qn, qnp1, brown_bcf )

      qn = qnp1

    end do

!   delete the model

    call delete ( stmodel )

  end subroutine equilibrate

end program cylinder3
