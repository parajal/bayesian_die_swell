! one-dimensional diffuse-interface problem
! three phase problem
! no flux boundary conditions
! Note: only Cahn_Hilliard equation (vp element only used as dummy element)

program diffuse_interface9

  use tfem_m
  use diffuse_interface_elements_m
  use diffuse_interface_functions_m
  use hsl_ma41_m
  use hsl_ma57_m
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

  type(meshgen_options_t) :: mesh_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefdi
  type(problem_t),target :: problem, problemdi
  type(sysmatrix_t) :: sysmatrixdi
  type(sysvector_t), target :: sol, rhsd
  type(vector_t) :: concentration1, concentration2
  type(oldvectors_t) :: oldvectors_di
  type(coefficients_t) :: coefficients
  type(sysvector_t), target :: soldi, rhsdi, soldin, soldi_prev
  type(subscript_t) :: cval, muval
  !type(solver_options_ma41_t) :: solver_options_di


! variables

  integer :: &
    nx=20,               & ! number of elements in x
    numtimesteps = 2,    & ! number of time steps
    itermax = 10,        & ! maximum number of Picard iterations
    outputstep = 10        ! Step for writing output files

  real(dp) :: &
    eta = 0.1_dp,      & ! viscosity
    deltat = 5.e-3_dp, & ! time step
    rho = 1.0_dp,      & ! density in diffuse-interface method
    alpha = 1.0_dp,    & ! parameter in the diffuse-interface method
    beta = 1.0_dp,     & ! parameter in the diffuse-interface method
    gamm = 0.0_dp,     & ! parameter in the diffuse-interface method
    Mcoef = 1.0_dp,    & ! parameter in the diffuse-interface method
    kappa = 1.0_dp,    & ! parameter in the diffuse-interface method
    epsdi = 1.e-3_dp,  & ! accuracy in Picard iteration
    rs_di  = 1.0_dp,   & ! real_storage for the diffuse-interface LU (HSL)
    is_di  = 1.4_dp      ! integer_storage for diffuse-interface LU (HSL)

  integer ::  step, iter, istep
  real(dp) :: cmax, cdiff, mumax, mudiff

  character(len=25) :: filename


! namelist for input of variables; read from standard input

  namelist /comppar/ nx, numtimesteps, outputstep, rho, eta, deltat, &
    alpha, beta, gamm, Mcoef, kappa, epsdi, itermax, rs_di, is_di

  read ( unit=*, nml=comppar )


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=200, ncoefr=150 )

! Stokes
! note: just to make a dummy list to use existing CH element (with velocities)

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

! diffuse interface
  coefficients%r(101:107) = [ deltat, Mcoef, alpha, beta, kappa, rho, gamm  ]
  coefficients%i(151) = 1

! create mesh

  mesh_options%elshape = 2 ! two-node line elements

  mesh_options%nx = nx  ! number of elements, equidistant
  mesh_options%lx = 1.0_dp ! length of interval

  call line1d ( mesh, mesh_options )

  call fill_mesh_parts ( mesh )


! problem definition of velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [1,1,1,   &  ! velocity
                 1,0,1,   &  ! pressure
                 1,1,1 ], &  ! scalar, such as vorticity
                 [3,3] )

  input_probdef%physq = [1,2]

  input_probdef%probnr = 1

! Dirichlet boundary conditions

  call problem_definition ( input_probdef, mesh, problem )


! problem definition diffuse interface phase 1

  call create_input_probdef ( mesh, input_probdefdi, nvec=5, nphysq=4 )

  input_probdefdi%vec_elementdof(1)%a =   &
      reshape ( [ 1,1,1,    &  ! c1
                  1,1,1,    &  ! c2
                  1,1,1,    &  ! mu1
                  1,1,1,    &  ! mu2
                  1,1,1  ], &  ! scalar for plotting
                  [3,5] )

  input_probdefdi%physq = [1,2,3,4]
  input_probdefdi%probnr = 2

  call problem_definition ( input_probdefdi, mesh, problemdi )


! create a vector subscript for the concentration and mu without Lagr. multipl.

  call create_subscript ( mesh, problemdi, cval, physqarr=[1] )
  call create_subscript ( mesh, problemdi, muval, physqarr=[3] )

! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )

  sol%u = 0.0_dp

! create system vectors (solution and right-hand side) for diffuse interface and
! initialize vectors plus a small random perturbation.

  call create_sysvector ( problemdi, soldi, rhsdi, soldin, soldi_prev )

  soldin%u = 0

  call fill_sysvector ( mesh, problemdi, soldin, degfd=1, &
       node1=1, node2=mesh%nnodes, func=difunc, funcnr=12)
  call fill_sysvector ( mesh, problemdi, soldin, degfd=2, &
       node1=1, node2=mesh%nnodes, func=difunc, funcnr=12)

  soldin%u(muval%s) = 0

  soldi%u = soldin%u

! create system matrix for diffuse-interface problem

  call create_sysmatrix_structure_base ( sysmatrixdi, mesh, problemdi )
  call finalize_sysmatrix_structure ( sysmatrixdi )

  call create_sysmatrix_data ( sysmatrixdi )

! create the structure oldvectors_di

  call create_oldvectors ( oldvectors_di, nsysvec=3, nprob=2 )

! store solution vectors and problem structures

  oldvectors_di%s(1)%p => sol
  oldvectors_di%s(2)%p => soldin
  oldvectors_di%s(3)%p => soldi
  oldvectors_di%p(1)%p => problem
  oldvectors_di%p(2)%p => problemdi

!
  step = 0

  call create_vector ( problemdi, concentration1, physq=1 )
  call extract_physvector ( mesh, problemdi, soldi, concentration1 )
  call create_vector ( problemdi, concentration2, physq=1 )
  call extract_physvector ( mesh, problemdi, soldi, concentration2 )

  write(filename,'(a,i5.5)') 'conc', step

  open ( unit=10, file=filename )

  do istep=1,mesh%nnodes
      write( unit=10, fmt='(1pE13.5E2,1pE13.5E2,1pE13.5E2)' ) mesh%coor(istep,1), &
           concentration1%u(istep)
  end do
  call delete ( concentration1 )
  call delete ( concentration2 )


! time stepping

  do step = 1, numtimesteps

!   iteration scheme for (c,mu)

    iter = 0

    picard: do

      iter = iter + 1

!     build (assemble) matrix and vector for diffuse-interface problem

      coefficients%r(102) = Mcoef
      coefficients%r(103) = alpha
      coefficients%r(104) = beta
      coefficients%r(105) = kappa

      call build_system ( mesh, problemdi, sysmatrixdi, sysvector=rhsdi, &
        elemsub=diffuse_interface_tp_coupled_elem, oldvectors=oldvectors_di, &
        coefficients=coefficients )

      call check ( sysmatrixdi )

      soldi_prev%u = soldi%u

!     solve (c,mu)

      !solver_options_di%real_storage=rs_di
      !solver_options_di%integer_storage=is_di

      call solve_system_ma41 ( sysmatrixdi, rhsdi, soldi )


      cmax  = maxval(abs(soldi%u(cval%s))) ! c1 maximum
      cdiff = maxval(abs(soldi%u(cval%s)-soldi_prev%u(cval%s))) ! c1 difference
      mumax = maxval(abs(soldi%u(muval%s))) ! mu1 maximum
      mudiff = maxval(abs(soldi%u(muval%s)-soldi_prev%u(muval%s))) ! mu1 diff

      print *, 'cmax = ', cmax, 'mumax = ', mumax
      print *, 'cdiff = ', cdiff, 'mudiff = ', mudiff
      print *, 'Picard iter = ', iter
      print *, 'Step = ', step

!    if ( cdiff < epsdi * cmax .and. mudiff < epsdi * mumax ) exit picard

      if ( cdiff < epsdi * cmax ) exit picard
      if ( iter >= itermax ) then
        write(*,'(a,i0)') ' too many iterations: ', itermax
        stop
      end if

    end do picard

    soldin%u = soldi%u


!   Writing output for the struct_factor program every outputstep time steps

    if ( mod( step,outputstep ) == 0 .or. step==1 ) then

    write(filename,'(a)') 'mesh_nodes.out'
    open ( unit=10, file=filename )

    write( unit=10, fmt=* ) mesh%nnodes

    close ( unit=10 )

    open ( unit=10, file='nstep.data')
    write ( unit=10, fmt=*) step, outputstep
    close ( unit=10 )

    call create_vector ( problemdi, concentration1, physq=1 )
    call create_vector ( problemdi, concentration2, physq=2 )
    call extract_physvector ( mesh, problemdi, soldi, concentration1 )
    call extract_physvector ( mesh, problemdi, soldi, concentration2 )

    write(filename,'(a,i5.5)') 'conc', step

    open ( unit=10, file=filename )

    do istep=1,mesh%nnodes
      write( unit=10, fmt='(1pE13.5E2,1pE13.5E2,1pE13.5E2)' ) mesh%coor(istep,1), &
           concentration1%u(istep)
    end do

    close ( unit=10 )

    write(filename,'(a,i5.5)') 'conc2', step

    open ( unit=10, file=filename )

    do istep=1,mesh%nnodes
      write( unit=10, fmt='(1pE13.5E2,1pE13.5E2,1pE13.5E2)' ) mesh%coor(istep,1), &
           concentration2%u(istep)
    end do

    close ( unit=10 )



    call delete ( concentration1 )
    call delete ( concentration2 )

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

  open ( unit=10, form='formatted', file='vel.out' )
  write(10,*) sol%u
  close(unit=10)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( oldvectors_di )
  call delete ( problemdi )
  call delete ( input_probdefdi )
  call delete ( sysmatrixdi )
  call delete ( soldi, rhsdi, soldin, soldi_prev )
  call delete ( cval, muval )
  call delete ( coefficients )

end program diffuse_interface9
