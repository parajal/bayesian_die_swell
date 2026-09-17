! diffuse-interface problem on a unit square
! three phase ternary system where the two compositions and chemical potentials
!    are solved  coupled
! linear Couette flow
! constant viscosity
! first-order time integration
! output to Tecplot

program diffuse_interface10

  use tfem_m
  use diffuse_interface_elements_m
  use diffuse_interface_functions_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m
  use tec_utils_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities and (c,mu)
    pintpl = 4,     & ! Q1 pressures
    physqvel = 1,   & ! physical quantity nr of the velocities
    physqpress = 2, & ! physical quantity nr of the pressures
    gauss = 3         ! 3x3 integration of quads


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefdi
  type(problem_t), target :: problem, problemdi
  type(sysmatrix_t) :: sysmatrix, sysmatrixdi
  type(lu_ma57_t) :: lu_u
  type(sysvector_t), target :: sol, rhsd
  type(oldvectors_t) :: oldvectors_di
  type(coefficients_t) :: coefficients
  type(sysvector_t), target :: soldi, rhsdi, soldin, soldi_prev
  type(subscript_t) :: cval, muval
  type(subscript_t) :: c2val, mu2val
  type(solver_options_ma57_t) :: solver_options_up
  type(solver_options_ma41_t) :: solver_options_di
  type(oldvectors_t) :: oldvectors, oldvectors_tec
  type(vector_t), target :: vel_x, vel_y, c, mu
  type(vector_t), target :: velocity, pressure, vorticity
  type(vector_t), target :: concentration
  type(vector_t), target :: concentration2
  type(vector_t), target :: concentration3


! variables

  integer :: &
    nx=20,               & ! number of elements in x
    ny=20,               & ! number of elements in y
    numtimesteps = 2,    & ! number of time steps
    itermax = 10           ! maximum number of Picard iterations

  real(dp) :: &
    eta = 0.1_dp,      & ! viscosity
    deltat = 5.e-3_dp, & ! time step
    rho = 1.0_dp,      & ! density in diffuse-interface method
    alpha = 1.0_dp,    & ! parameter in the diffuse-interface method
    beta = 1.0_dp,     & ! parameter in the diffuse-interface method
    gamm = 1.0_dp,     & ! parameter in the diffuse-interface method
    Mcoef = 1.0_dp,    & ! parameter in the diffuse-interface method
    kappa = 1.0_dp,    & ! parameter in the diffuse-interface method
    U = 1.0_dp,        & ! velocity difference upper-lower wall
    epsdi = 1.e-3_dp,  & ! accuracy in Picard iteration
    rs_up = 1.0_dp,    & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.0_dp,    & ! integer_storage velocity-pressure LU (HSL)
    rs_di  = 1.0_dp,   & ! real_storage for the diffuse-interface LU (HSL)
    is_di  = 1.4_dp      ! integer_storage for diffuse-interface LU (HSL)

  integer ::  step, iter, i
  real(dp) :: cmax, cdiff, mumax, mudiff, c2max, c2diff

  integer,parameter :: num_tec_data = 7  ! number of data sets in tecplot output

  character(len=20), dimension(num_tec_data) :: dataname

  logical :: append


! namelist for input of variables; read from standard input

  namelist /comppar/ nx, ny, numtimesteps, rho, eta, deltat, &
    alpha, beta, gamm, Mcoef, kappa, U, epsdi, itermax, rs_up, is_up, rs_di, is_di

  read ( unit=*, nml=comppar )


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=200, ncoefr=150 )

! Stokes
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

! define type of diffusion
!    0:		Lowengrub type
!			the mobilty is constant
!    1:		Mauri, Molin, Anderson type
!			the mobilty is a function of c
  coefficients%i(151) = 1

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

! diffuse interface
  coefficients%r(101:107) = [ deltat, Mcoef, alpha, beta, kappa, rho, gamm ]


! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5
  call add_to_mesh ( mesh, curve=[-3] )     ! curve 6

  call fill_mesh_parts ( mesh )


! problem definition of velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

  input_probdef%probnr = 1

! Dirichlet boundary conditions
! fix velocity and pressure in one point!

  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )
  call define_essential ( mesh, input_probdef, point=1, physq=physqvel )

! constraints for periodical boundary conditions

! velocities
  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=2, curve2=5, discretization='collocation' )
  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=1, curve2=6, discretization='collocation', &
    exclude=1 )

  call problem_definition ( input_probdef, mesh, problem )


! problem definition diffuse interface

  call create_input_probdef ( mesh, input_probdefdi, nvec=5, nphysq=4 )

  input_probdefdi%vec_elementdof(1)%a =   &
      reshape ( [ 1,1,1,1,1,1,1,1,1,    &  ! c1
                  1,1,1,1,1,1,1,1,1,    &  ! c2
                  1,1,1,1,1,1,1,1,1,    &  ! mu1
                  1,1,1,1,1,1,1,1,1,    &  ! mu2
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar for plotting
                  [9,5] )

  input_probdefdi%physq = [1,2,3,4]
  input_probdefdi%probnr = 2

! constraint for periodical boundary conditions of the (c,mu)
  call define_constraint ( mesh, input_probdefdi, curve1=2, curve2=5, &
    discretization='collocation' )
  call define_constraint ( mesh, input_probdefdi, curve1=1, curve2=6, &
    discretization='collocation', exclude=1 )

  call problem_definition ( input_probdefdi, mesh, problemdi )


! create a vector subscript for the concentration and mu without Lagr. multipl.

  call create_subscript ( mesh, problemdi, cval, physqarr=[1] )
  call create_subscript ( mesh, problemdi, c2val, physqarr=[2] )
  call create_subscript ( mesh, problemdi, muval, physqarr=[3] )
  call create_subscript ( mesh, problemdi, mu2val, physqarr=[4] )


! create system vectors for velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )


! set pressure level = 0 in lower left corner
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqpress, value=0._dp )

! set velocity = 0 in lower left corner
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqvel, value=0._dp )

! create system matrix of velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! build (assemble) matrix and vector for velocity/pressure problem

! stokes velocity/pressure
  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

! periodical condition on velocities

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_constr_node_conn, addmatvec=.true. )

  call check ( sysmatrix )


! solve initial velocity/pressure consistent with initial conditions
! and keep decomposition of the matrix
! This generates a linear profile.

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options_up%real_storage=rs_up
  solver_options_up%integer_storage=is_up

  call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
    solver_options=solver_options_up )


! create system vectors (solution and right-hand side) for diffuse interface and
! initialize vectors plus a small random perturbation.

  call create_sysvector ( problemdi, soldi, rhsdi, soldin, soldi_prev )

  call random_number (soldin%u)

  call fill_sysvector ( mesh, problemdi, soldin, degfd=1, &
       node1=1, node2=mesh%nnodes, func=difunc, funcnr=13)

  soldin%u(c2val%s) = 0.7_dp - soldin%u(cval%s)

  soldin%u(muval%s) = 0.0_dp
  soldin%u(mu2val%s) = 0.0_dp

  soldi%u = soldin%u


  cmax  = maxval(abs(soldi%u(cval%s))) ! c maximum
  c2max  = maxval(abs(soldi%u(c2val%s))) ! c maximum
  print *, 'cmax = ', cmax!, 'mumax = ', mumax
  print *, 'c2max = ', c2max

! create system matrix for diffuse-interface problem

  call create_sysmatrix_structure_base ( sysmatrixdi, mesh, problemdi )
  call create_sysmatrix_structure_constraint ( sysmatrixdi, mesh, problemdi )
  call finalize_sysmatrix_structure ( sysmatrixdi )

  call create_sysmatrix_data ( sysmatrixdi )


! create vectors to extract scalar info

  call create_vector ( problemdi, concentration, physq=1 )
  call create_vector ( problemdi, concentration2, physq=2 )
  call create_vector ( problemdi, concentration3, vec=5 )
  call create_vector ( problemdi, c,     vec=5 )
  call create_vector ( problemdi, mu,     vec=5 )

  call create_vector ( problem, velocity,  physq=1 )
  call create_vector ( problem, pressure,  vec=3 )
  call create_vector ( problem, vorticity, vec=3 )
  call create_vector ( problem, vel_x,     vec=3 )
  call create_vector ( problem, vel_y,     vec=3 )

! create the structure oldvectors_di

  call create_oldvectors ( oldvectors_di, nsysvec=3, nprob=2 )

! create the structure oldvectors, used to obtain pressure

  call create_oldvectors ( oldvectors,     nsysvec=1 )

! create the structure oldvectors_tec for writing Tecplot files

  call create_oldvectors ( oldvectors_tec, nvec=num_tec_data )

! define data set for tecplot output

  dataname(1) = 'u'
  dataname(2) = 'v'
  dataname(3) = 'pressure'
  dataname(4) = 'vorticity'
  dataname(5) = 'c1'
  dataname(6) = 'c2'
  dataname(7) = 'c3'

  oldvectors_tec%v(1)%p => vel_x
  oldvectors_tec%v(2)%p => vel_y
  oldvectors_tec%v(3)%p => pressure
  oldvectors_tec%v(4)%p => vorticity
  oldvectors_tec%v(5)%p => concentration
  oldvectors_tec%v(6)%p => concentration2
  oldvectors_tec%v(7)%p => concentration3

! store solution vectors and problem structures

  sol%u = 0.0_dp
  oldvectors_di%s(1)%p => sol
  oldvectors_di%s(2)%p => soldin
  oldvectors_di%s(3)%p => soldi
  oldvectors_di%p(1)%p => problem
  oldvectors_di%p(2)%p => problemdi

  time = 0.0_dp

! time stepping

  do step = 1, numtimesteps

    time = time + deltat

    print *, ' '
    print *, '----------------------------------------- '
    print *, 'Step: ', step
    print *, 'deltat: ', deltat

!   iteration scheme for (c,mu)

    iter = 0

    picard: do

      print *, '   iter: ', iter

      iter = iter + 1

!     build (assemble) matrix and vector for diffuse-interface problem

      call build_system ( mesh, problemdi, sysmatrixdi, sysvector=rhsdi, &
        elemsub=diffuse_interface_tp_coupled_elem, oldvectors=oldvectors_di, &
        coefficients=coefficients )

!     periodical condition on (c,mu)

      call build_system_constraint ( mesh, problemdi, sysmatrixdi, &
        sysvector=rhsdi, elemsub=stokes_constr_node_conn, &
        addmatvec=.true. )

      call check ( sysmatrixdi )

      soldi_prev%u = soldi%u

!     solve (c1,c2,mu1,mu2)

      solver_options_di%real_storage=rs_di
      solver_options_di%integer_storage=is_di

      call solve_system_ma41 ( sysmatrixdi, rhsdi, soldi, &
        solver_options=solver_options_di )


      cmax  = maxval(abs(soldi%u(cval%s))) ! c maximum
      c2max  = maxval(abs(soldi%u(c2val%s))) ! c2 maximum
      cdiff = maxval(abs(soldi%u(cval%s)-soldi_prev%u(cval%s))) ! c difference
      c2diff = maxval(abs(soldi%u(c2val%s)-soldi_prev%u(c2val%s))) ! c2 difference
      mumax = maxval(abs(soldi%u(muval%s))) ! mu maximum
      mudiff = maxval(abs(soldi%u(muval%s)-soldi_prev%u(muval%s))) ! mu diff

      print *, 'cmax = ', cmax, 'mumax = ', mumax
      print *, 'c2max = ', c2max
      print *, 'cdiff = ', cdiff, 'mudiff = ', mudiff
      print *, 'c2diff = ', c2diff

  !    if ( cdiff < epsdi * cmax .and. mudiff < epsdi * mumax ) exit picard

      if ( cdiff < epsdi * cmax ) exit picard

      if ( iter >= itermax ) then
        write(*,'(a,i0)') ' too many iterations: ', itermax
        stop
      end if

    end do picard

    soldin%u = soldi%u

!   build (assemble) vector for velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_tp_coupled_mugradc, &
      oldvectors=oldvectors_di, physqrow=[physqvel], physqcol=[physqvel], &
      buildmatrix = .false., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve velocity/pressure problem

    solver_options_up%real_storage=rs_up
    solver_options_up%integer_storage=is_up

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
      solver_options=solver_options_up )

!   extract velocity
    call extract_physvector ( mesh, problem, sol, velocity )

    oldvectors%s(1)%p => sol

!   extract pressure and vorticity

    call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
      coefficients=coefficients, oldvectors=oldvectors )

    coefficients%i(13)=5
    call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
      coefficients=coefficients, oldvectors=oldvectors )

!   compute concentration3
    call extract_physvector ( mesh, problemdi, soldi, concentration )
    call extract_physvector ( mesh, problemdi, soldi, concentration2 )

    concentration3%u = 1.0_dp - concentration%u - concentration2%u

!   Tecplot output

    do i = 1, mesh%nnodes
      vel_x%u(i) = velocity%u(problem%vec_nodnumdegfd(i,velocity%vec) + 1)
      vel_y%u(i) = velocity%u(problem%vec_nodnumdegfd(i,velocity%vec) + 2)
    end do

    if (step == 1) then
      append = .false.
    else
      append = .true.
    end if

    if (mod(step,1) == 0) then
      call write_mdata_tecplot (mesh, problem, 'dim_tp.plt', dataname, &
        oldvectors_tec, append=append, time=time )
    end if

  end do


! write data for post-processing

  call write_mesh ( mesh, filename='mesh.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh, input_probdefdi, filename='probdefdi.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol%u
  write(10) soldi%u

  close(unit=10)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( lu_u )
  call delete ( oldvectors_di )
  call delete ( problemdi )
  call delete ( input_probdefdi )
  call delete ( sysmatrixdi )
  call delete ( soldi, rhsdi, soldin, soldi_prev )
  call delete ( cval, muval )
  call delete ( c2val, muval )
  call delete ( coefficients )

end program diffuse_interface10
