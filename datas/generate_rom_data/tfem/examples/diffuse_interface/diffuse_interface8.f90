! diffuse-interface problem on a square
! no external flow, bi-periodic boundary conditions
! & and a time-periodic source term for the momentum equation
! constant viscosity
! first-order time integration
! uses tec_utils for output

program diffuse_interface8

  use tfem_m
  use stokes_functions_m
  use diffuse_interface_functions_m
  use diffuse_interface_elements_m
  use hsl_ma41_m
  use io_utils_m


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
  type(lu_ma41_t) :: lu_u
  type(sysvector_t) , target:: sol, rhsd
  type(oldvectors_t) :: oldvectors_di
  type(coefficients_t) :: coefficients
  type(sysvector_t) , target:: soldi, rhsdi, soldin, soldi_prev
  type(subscript_t) :: cval, muval
  type(vector_t), target :: concentration
  type(oldvectors_t) :: oldvectors, oldvectors_tec
  type(vector_t),target :: vel_x, vel_y, c, mu
  type(vector_t),target :: velocity, pressure, vorticity
  !type(solver_options_ma41_t) :: solver_options_up


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
    Mcoef = 1.0_dp,    & ! parameter in the diffuse-interface method
    kappa = 1.0_dp,    & ! parameter in the diffuse-interface method
    U = 0.0_dp,        & ! velocity difference upper-lower wall
    epsdi = 1.e-3_dp,  & ! accuracy in Picard iteration
    rs_up = 1.0_dp,    & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.0_dp,    & ! integer_storage velocity-pressure LU (HSL)
    rs_di  = 1.0_dp,   & ! real_storage for the diffuse-interface LU (HSL)
    is_di  = 1.4_dp      ! integer_storage for diffuse-interface LU (HSL)

  integer ::  step, iter, istep, i
  real(dp) :: cmax, cdiff, period!, mudiff, mumax

  integer, parameter :: num_tec_data = 5 ! number of data sets in tecplot output

  character(len=20) :: filename
  character(len=20), dimension(num_tec_data) :: dataname

  logical :: append





! namelist for input of variables; read from standard input

  namelist /comppar/ nx, ny, numtimesteps, rho, eta, deltat, &
    alpha, beta, Mcoef, kappa, U, epsdi, itermax, rs_up, is_up, rs_di, is_di

  read ( unit=*, nml=comppar )


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=200, ncoefr=150 )

! Stokes
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0
  coefficients%i(14)  = 0

  coefficients%r(1)  = eta
  coefficients%r(2:) = 0

! diffuse interface
  coefficients%r(101:106) = [ deltat, Mcoef, alpha, beta, kappa, rho ]

  coefficients%vfunc => vfunc

! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny


! meshgen_options%regionshape = 2     ! quadrilateral with straight boundaries
  meshgen_options%x2d = &
    reshape ( [ -1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp,    &
                -1.0_dp,-1.0_dp, 1.0_dp,  1.0_dp ], &
              [4,2] )

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

  call define_essential ( mesh, input_probdef, point=1, physq=physqvel )
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! constraints for periodical boundary conditions

! velocities
  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=2, curve2=5, discretization='collocation' )
  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=1, curve2=6, discretization='collocation', &
    exclude=1 )


  call problem_definition ( input_probdef, mesh, problem )




! problem definition diffuse interface

  call create_input_probdef ( mesh, input_probdefdi, nvec=3, nphysq=2 )

  input_probdefdi%vec_elementdof(1)%a =   &
      reshape ( [ 1,1,1,1,1,1,1,1,1,    &  ! c
                  1,1,1,1,1,1,1,1,1,    &  ! mu
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar for plotting
                  [9,3] )

  input_probdefdi%physq = [1,2]
  input_probdefdi%probnr = 2


! constraint for periodical boundary conditions of the (c,mu)
  call define_constraint ( mesh, input_probdefdi, curve1=2, curve2=5, &
    discretization='collocation' )
  call define_constraint ( mesh, input_probdefdi, curve1=1, curve2=6, &
    discretization='collocation', exclude=1 )


  call problem_definition ( input_probdefdi, mesh, problemdi )


! create a vector subscript for the concentration and mu without Lagr. multipl.

  call create_subscript ( mesh, problemdi, cval, physqarr=[1] )
  call create_subscript ( mesh, problemdi, muval, physqarr=[2] )


! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )


! fill solution vector with essential boundary conditions


! set velocity level = 0 in lower left corner
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqvel, value=0._dp )
! set pressure level = 0 in lower left corner
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqpress, value=0._dp )


! create system matrix of gradient/velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! build (assemble) matrix and vector for gradient/velocity/pressure problem

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

  !solver_options_up%real_storage=rs_up
  !solver_options_up%integer_storage=is_up

  call solve_system_ma41 ( sysmatrix, rhsd, sol, lu_u )


! create system vectors (solution and right-hand side) for diffuse interface and
! initialize vectors plus a small random perturbation.

  call create_sysvector ( problemdi, soldi, rhsdi, soldin, soldi_prev )

! create concentration vector
!   funcnr = 1 : unstable stratified layers
!   funcnr = 2 : retraction of non-circular drop
!   funcnr = 3 : spinodal decomposition
!   funcnr = 4 : a supercritical and a subcritical nucleus

  soldin%u = 0

  call fill_sysvector ( mesh, problemdi, soldin, degfd=1, &
       node1=1, node2=mesh%nnodes, func=difunc, funcnr=3)

  soldi%u = soldin%u

! create system matrix for diffuse-interface problem

  call create_sysmatrix_structure_base ( sysmatrixdi, mesh, problemdi )
  call create_sysmatrix_structure_constraint ( sysmatrixdi, mesh, problemdi )
  call finalize_sysmatrix_structure ( sysmatrixdi )

  call create_sysmatrix_data ( sysmatrixdi )


! create the structure oldvectors_di

  call create_oldvectors ( oldvectors_di, nsysvec=3, nprob=2 )

  call create_oldvectors ( oldvectors,     nsysvec=1 )
  call create_oldvectors ( oldvectors_tec, nvec=num_tec_data )


  call create_vector ( problem, velocity,  physq=1 )
  call create_vector ( problem, pressure,  vec=3 )
  call create_vector ( problem, vorticity, vec=3 )
  call create_vector ( problem, vel_x,     vec=3 )
  call create_vector ( problem, vel_y,     vec=3 )

  call create_vector ( problemdi, concentration, physq=1 )

  call create_vector ( problemdi, c,     vec=3 )
  call create_vector ( problemdi, mu,     vec=3 )

! define data set for tecplot output

  dataname(1) = 'u'
  dataname(2) = 'v'
  dataname(3) = 'pressure'
  dataname(4) = 'vorticity'
  dataname(5) = 'c'

  oldvectors_tec%v(1)%p => vel_x
  oldvectors_tec%v(2)%p => vel_y
  oldvectors_tec%v(3)%p => pressure
  oldvectors_tec%v(4)%p => vorticity
  oldvectors_tec%v(5)%p => concentration


! store solution vectors and problem structures

  oldvectors_di%s(1)%p => sol
  oldvectors_di%s(2)%p => soldin
  oldvectors_di%s(3)%p => soldi
  oldvectors_di%p(1)%p => problem
  oldvectors_di%p(2)%p => problemdi

  time = 0._dp
  period = 0.5_dp

! time stepping

  do step = 1, numtimesteps

    time = time + deltat
!
!   create the time period flow by alternating the source term
!     first half period:  coefficient(14) = 1
!     second half period: coefficient(14) = 2
!
    if  ( mod(time,period) <= period/2._dp)  coefficients%i(14)  =  1
    if  ( mod(time,period) > period/2._dp)  coefficients%i(14)  =  2

    print *, 'time, step = ',time, step

!   iteration scheme for (c,mu)

    iter = 0

    picard: do

      iter = iter + 1

!     build (assemble) matrix and vector for diffuse-interface problem

      call build_system ( mesh, problemdi, sysmatrixdi, sysvector=rhsdi, &
        elemsub=diffuse_interface_elem, oldvectors=oldvectors_di, &
        coefficients=coefficients )

!     periodical condition on (c,mu)
      call build_system_constraint ( mesh, problemdi, sysmatrixdi, &
        sysvector=rhsdi, elemsub=stokes_constr_node_conn, &
        addmatvec=.true. )

      call check ( sysmatrixdi )

      soldi_prev%u = soldi%u

!     solve (c,mu)

      !solver_options_up%real_storage=rs_up
      !solver_options_up%integer_storage=is_up

      call solve_system_ma41 ( sysmatrixdi, rhsdi, soldi )


      cmax  = maxval(abs(soldi%u(cval%s))) ! c maximum
      cdiff = maxval(abs(soldi%u(cval%s)-soldi_prev%u(cval%s))) ! c difference
      !mumax = maxval(abs(soldi%u(muval%s))) ! mu maximum
      !mudiff = maxval(abs(soldi%u(muval%s)-soldi_prev%u(muval%s))) ! mu diff

  !   print *, 'cmax = ', cmax, 'mumax = ', mumax
  !   print *, 'cdiff = ', cdiff, 'mudiff = ', mudiff

  !    if ( cdiff < epsdi * cmax .and. mudiff < epsdi * mumax ) exit picard

      if ( cdiff < epsdi * cmax ) exit picard

      if ( iter >= itermax ) then
        write(*,'(a,i0)') ' too many iterations: ', itermax
        stop
      end if

    end do picard

    soldin%u = soldi%u

!   build (assemble) vector for velocity/pressure problem
    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
      coefficients=coefficients,  buildmatrix = .false. )

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_mugradc, &
      oldvectors=oldvectors_di, physqrow=[physqvel], physqcol=[physqvel], &
      buildmatrix = .false., coefficients=coefficients, addmatvec=.true.)

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve velocity/pressure problem

    !solver_options_up%real_storage=rs_up
    !solver_options_up%integer_storage=is_up

    call solve_system_ma41 ( sysmatrix, rhsd, sol, lu_u )

!   extract velocity
    call extract_physvector ( mesh, problem, sol, velocity )

    oldvectors%s(1)%p => sol

!   extract pressure and vorticity

    call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
      coefficients=coefficients, oldvectors=oldvectors )

    coefficients%i(13)=5
    call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
      coefficients=coefficients, oldvectors=oldvectors )

    call extract_physvector ( mesh, problemdi, soldi, concentration )

    write(filename,'(a,i4.4)') 'conc', step

    open ( unit=10, file=filename )

    do istep=1,mesh%nnodes
      write( unit=10, fmt='(1pE13.5E2,1pE13.5E2,1pE13.5E2)' ) mesh%coor(istep,1), &
           mesh%coor(istep,2), &
           concentration%u(istep)
    end do

    close ( unit=10 )

!   Tecplot output

    do i = 1, mesh%nnodes
      vel_x%u(i) = velocity%u(problem%vec_nodnumdegfd(i,velocity%vec) + 1)
      vel_y%u(i) = velocity%u(problem%vec_nodnumdegfd(i,velocity%vec) + 2)
    end do

    if (step == 10) then
      append = .false.
    else
      append = .true.
    end if

    if (mod(step,10) == 0) then
      call write_mdata_tecplot (mesh, problem, 'chaotic_dim.plt', dataname, &
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

  call delete ( concentration )
  call delete ( coefficients )

end program diffuse_interface8
