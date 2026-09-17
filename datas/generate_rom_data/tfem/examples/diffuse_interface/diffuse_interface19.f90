! periodic condition (left-right)
! zero velocity on lower and upper boundaries, zero pressure on left upper coner
! first-order implicit Euler time integration
! Picard iteration

program diffuse_interface19

  use tfem_m
  use diffuse_interface_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m
  use gmsh_utils_m

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
  type(sysvector_t), target :: sol, soln, soldi, soldin
  type(sysvector_t) :: rhsd, rhsdi, soldi_prev
  type(vector_t) :: concentration1, concentration2, concentration3
  type(oldvectors_t) :: oldvectors_di
  type(coefficients_t) :: coefficients
  type(subscript_t) :: cval, muval, c2val, mu2val
  type(solver_options_ma41_t) :: solver_options_di
  type(solver_options_ma57_t) :: solver_options_up
  type(lu_ma57_t) :: lu_u

! variables

  integer :: &
    nx=20,               & ! number of elements in x
    ny=20,               & ! number of elements in y
    numtimesteps = 2,    & ! number of time steps
    itermax = 10           ! maximum number of Picard iterations

  real(dp) :: &
    eta = 0.1_dp,      & ! viscosity
    deltat = 5.e-3_dp, & ! time step
    rho = 1.0_dp,      & ! density in the diffuse-interface
    alpha = 1.0_dp,    & ! parameter in the diffuse-interface
    beta = 1.0_dp,     & ! parameter in the diffuse-interface
    gamm = 1.0_dp,     & ! parameter in the diffuse-interface
    omega = 1.0_dp,    & ! parameter in the diffuse-interface
    Mcoef = 1.0_dp,    & ! parameter in the diffuse-interface
    kappa = 1.0_dp,    & ! parameter in the diffuse-interface
    kappa2 = 1.0_dp,   & ! parameter in the diffuse-interface
    kappa3 = 1.0_dp,   & ! parameter in the diffuse-interface
    epsdi = 1.e-3_dp,  & ! accuracy in Picard iteration (c,mu)
    rs_up = 1.0_dp,    & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.0_dp,    & ! integer_storage velocity-pressure LU (HSL)
    rs_di  = 1.0_dp,   & ! real_storage for the diffuse-interface LU (HSL)
    is_di  = 1.4_dp      ! integer_storage for diffuse-interface LU (HSL)

  integer ::  step, iter, i

  real(dp) :: cmax, cdiff, mumax!, mudiff
  real(dp) :: c2max, c2diff, mu2max!, mu2diff

  character(len=20) :: filename

! namelist for input of variables; read from standard input

  namelist /comppar/ nx, ny, numtimesteps, rho, eta, deltat, &
    alpha, beta, gamm, omega, Mcoef, kappa, kappa2, kappa3, &
    epsdi, itermax, rs_up, is_up, rs_di, is_di

  read ( unit=*, nml=comppar )

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=200, ncoefr=150 )

! Stokes
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

! diffuse interface
  coefficients%r(101:107) = [ deltat, Mcoef, alpha, beta, kappa, rho, gamm ]
  coefficients%r(109:111) = [ kappa2, kappa3, omega ]

! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  meshgen_options%regionshape = 2     ! quadrilateral with straight boundaries
  meshgen_options%x2d = &
    reshape ( [ -1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp,    &
                -1.0_dp,-1.0_dp, 1.0_dp,  1.0_dp ], &
              [4,2] )

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )

! problem definition of velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar for plot
                 [9,3] )

  input_probdef%physq = [1,2]

  input_probdef%probnr = 1

! Dirichlet boundary conditions

  call define_essential ( mesh, input_probdef, curve1=1, physq=physqvel )
  call define_essential ( mesh, input_probdef, curve1=3, physq=physqvel )
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! constraints for periodical boundary conditions

! velocities

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=2, curve2=5, discretization='collocation', &
    exclude=3 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors for velocity/pressure

  call create_sysvector ( problem, sol, soln, rhsd )

! fill solution vector with essential boundary conditions

! set velocity level = 0 in lower boudary
  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, physq=physqvel, value=0._dp )
! set velocity level = 0 in upper boudary
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=physqvel, value=0._dp )
! set pressure = 0 at left lower coner
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqpress, value=0._dp )

! create system matrix of velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! problem definition diffuse interface

  call create_input_probdef ( mesh, input_probdefdi, nvec=5, nphysq=4 )

  input_probdefdi%vec_elementdof(1)%a =  &
      reshape ( [ 1,1,1,1,1,1,1,1,1,    &  ! c1
                  1,1,1,1,1,1,1,1,1,    &  ! c2
                  1,1,1,1,1,1,1,1,1,    &  ! mu1
                  1,1,1,1,1,1,1,1,1,    &  ! mu2
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar for plot
                  [9,5] )

  input_probdefdi%physq = [1,2,3,4]
  input_probdefdi%probnr = 2

! constraint for periodical boundary conditions of c and mu

  call define_constraint ( mesh, input_probdefdi, curve1=2, curve2=5, &
    discretization='collocation' )

  call problem_definition ( input_probdefdi, mesh, problemdi )


! create a vector subscript for c and mu

  call create_subscript ( mesh, problemdi, cval, physqarr=[1] )
  call create_subscript ( mesh, problemdi, c2val, physqarr=[2] )
  call create_subscript ( mesh, problemdi, muval, physqarr=[3] )
  call create_subscript ( mesh, problemdi, mu2val, physqarr=[4] )

! create system vectors for diffuse interface

  call create_sysvector ( problemdi, soldi, rhsdi, soldin, soldi_prev )

! initialize vectors

  soldin%u = 0

! funcnr = 1 : circular drop (phase 1)
  call fill_sysvector ( mesh, problemdi, soldin, physq=1, &
       node1=1, node2=mesh%nnodes, func=difunc_ternary, funcnr=1 )

! funcnr = 2 : stable layer (phase 2)
  call fill_sysvector ( mesh, problemdi, soldin, physq=2, &
       node1=1, node2=mesh%nnodes, func=difunc_ternary, funcnr=2 )

  soldin%u(c2val%s) = soldin%u(c2val%s) - soldin%u(cval%s)

  do i = 1, mesh%nnodes
    if ( soldin%u(c2val%s(i)) < 0._dp ) &
      soldin%u(c2val%s(i)) = 0._dp
  end do

  soldin%u(muval%s) = 0._dp
  soldin%u(mu2val%s) = 0._dp

  soldi%u = soldin%u

! create system matrix for diffuse-interface problem

  call create_sysmatrix_structure_base ( sysmatrixdi, mesh, problemdi )
  call create_sysmatrix_structure_constraint ( sysmatrixdi, mesh, problemdi )
  call finalize_sysmatrix_structure ( sysmatrixdi )

  call create_sysmatrix_data ( sysmatrixdi )

! create the structure oldvectors

  call create_oldvectors ( oldvectors_di, nsysvec=4, nprob=2 )

! create vectors

  call create_vector ( problemdi, concentration1, physq=1 )
  call create_vector ( problemdi, concentration2, physq=2 )
  call create_vector ( problemdi, concentration3, vec=5 )

! store solution vectors, problem structures

  oldvectors_di%s(1)%p => sol
  oldvectors_di%s(2)%p => soldin
  oldvectors_di%s(3)%p => soldi
  oldvectors_di%s(4)%p => soln

  oldvectors_di%p(1)%p => problem
  oldvectors_di%p(2)%p => problemdi

! build (assemble) matrix and vector for velocity/pressure problem
! stokes velocity/pressure

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

! periodical condition on velocities

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_constr_node_conn, addmatvec=.true. )

  call check ( sysmatrix )

! solve initial velocity/pressure of Stokes flow

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options_up%real_storage=rs_up
  solver_options_up%integer_storage=is_up

  call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
    solver_options=solver_options_up )

  soln%u = sol%u

! time stepping

  do step = 1, numtimesteps

    print *, ' '
    print *, 'step ', step, '------------------------------------------'

!   iteration scheme for (c,mu)

    iter = 0

    picard: do

      iter = iter + 1

!     build (assemble) matrix and vector for diffuse-interface problem

      call build_system ( mesh, problemdi, sysmatrixdi, sysvector=rhsdi, &
        elemsub=diffuse_interface_tp2_coupled_elem, oldvectors=oldvectors_di, &
        coefficients=coefficients )

!     periodical condition on (c,mu)

      call build_system_constraint ( mesh, problemdi, sysmatrixdi, &
        sysvector=rhsdi, elemsub=stokes_constr_node_conn, &
        addmatvec=.true. )

      call check ( sysmatrixdi )

      soldi_prev%u = soldi%u

!     solve (c,mu)

      solver_options_di%real_storage=rs_di
      solver_options_di%integer_storage=is_di

      call solve_system_ma41 ( sysmatrixdi, rhsdi, soldi, &
        solver_options=solver_options_di )

      cmax  = maxval(abs(soldi%u(cval%s)))   ! c1 maximum
      c2max  = maxval(abs(soldi%u(c2val%s))) ! c2 maximum
      cdiff = maxval(abs(soldi%u(cval%s)-soldi_prev%u(cval%s)))    ! c1 difference
      c2diff = maxval(abs(soldi%u(c2val%s)-soldi_prev%u(c2val%s))) ! c2 difference
      mumax = maxval(abs(soldi%u(muval%s)))   ! mu1 maximum
      mu2max = maxval(abs(soldi%u(mu2val%s))) ! mu2 maximum
      !mudiff = maxval(abs(soldi%u(muval%s)-soldi_prev%u(muval%s)))    ! mu1 diff
      !mu2diff = maxval(abs(soldi%u(mu2val%s)-soldi_prev%u(mu2val%s))) ! mu2 diff

      write(*,'(a,es13.3,a,es13.3)') ' cmax  = ', cmax,  ' mumax  = ', mumax
      write(*,'(a,es13.3,a,es13.3)') ' c2max = ', c2max, ' mu2max = ', mu2max

      if ( cdiff < epsdi * cmax .and. c2diff < epsdi * c2max ) then
        write(*,'(a,i10)') ' iteration: ', iter
        write(*,'(a,2es13.3)') ' cdiff, c2diff: ', cdiff, c2diff
        exit picard
      end if

      if ( iter >= itermax ) then
        write(*,'(a,i0)') ' too many iterations: ', itermax
        stop
      end if

    end do picard

    soldin%u = soldi%u

!   build (assemble) vector for velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_tp2_coupled_mugradc, &
      oldvectors=oldvectors_di, physqrow=[physqvel], physqcol=[physqvel], &
      buildmatrix = .false., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve velocity/pressure problem

    solver_options_up%real_storage=rs_up
    solver_options_up%integer_storage=is_up

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
      solver_options=solver_options_up )

    soln%u = sol%u

!   compute c3 vector for post processing

    call extract_physvector ( mesh, problemdi, soldi, concentration1)
    call extract_physvector ( mesh, problemdi, soldi, concentration2)
    concentration3%u = 1._dp - concentration1%u - concentration2%u

!   write files

    if ( step == 1 .or. mod(step,50) == 0 ) then

      write(filename,'(a,i4.4,a)') 'post', step, '.vtk'

      call write_scalar_vtk ( mesh, problemdi, physq=1,&
          filename=filename, sysvector=soldi, dataname='c1', &
          append = .false. )
      call write_scalar_vtk ( mesh, problemdi, physq=2,&
          filename=filename, sysvector=soldi, dataname='c2', &
          append = .true. )
      call write_scalar_vtk ( mesh, problemdi, &
          filename=filename, vector=concentration3, dataname='c3', &
          append = .true. )
      call write_vector_vtk ( mesh, problem, physq=1, &
          filename=filename, sysvector=sol, dataname='vel', &
          append = .true. )

    end if

  end do

  call delete ( problem, problemdi )
  call delete ( sysmatrix, sysmatrixdi )
  call delete ( input_probdef, input_probdefdi )
  call delete ( mesh )
  call delete ( sol, soln, soldi, soldin )
  call delete ( rhsd, rhsdi, soldi_prev )
  call delete ( oldvectors_di )
  call delete ( coefficients )
  call delete ( cval, c2val, muval, mu2val )
  call delete ( concentration1, concentration2, concentration3 )
  call delete ( lu_u )

contains

  function difunc_ternary ( nr, x )

    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: difunc_ternary
    real(dp) :: r, C

    C = 2.e-2_dp

    select case(nr)
      case(1)
!       circular drop
        r = ( - sqrt( (x(1))**2 + (x(2))**2 ) + 0.3_dp ) / ( 2._dp * C )
        difunc_ternary = 0.5_dp * ( tanh(r) + 1._dp )

      case(2)
!       stable layer
        r = ( x(2) ) / ( 2._dp * C )
        difunc_ternary = 0.5_dp * ( tanh(r) + 1._dp )

      case default
        write(*,'(/a,i0/)') 'Error difunc_ternary: wrong function number: ', nr
        stop
    end select

  end function difunc_ternary

end program diffuse_interface19
