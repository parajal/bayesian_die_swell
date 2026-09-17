! diffuse-interface problem in a brick
! no external flow, homogeneous essential boundary conditions for velocities
! constant viscosity
! first-order time integration
! output using vtk
! the initial concentration is such for spinodal decomposition

program diffuse_interface4

  use tfem_m
  use diffuse_interface_functions_m
  use diffuse_interface_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m
  use figplot_m

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
  type(sysvector_t), target :: sol, soldi, soldin
  type(sysvector_t) :: rhsd, rhsdi, soldi_prev
  type(oldvectors_t) :: oldvectors_di
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients
  type(subscript_t) :: cval, muval
  type(solver_options_ma57_t) :: solver_options_up
  type(solver_options_ma41_t) :: solver_options_di



! variables


  integer :: presnod(8) = [1,3,9,7,19,21,27,25]


  integer :: &
    nx=10,               & ! number of elements in x
    ny=10,               & ! number of elements in y
    nz=10,               & ! number of elements in y
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
    epsdi = 1.e-3_dp     ! accuracy in Picard iteration
  real(dp) :: &
    rs_up = 1.0_dp,    & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.0_dp,    & ! integer_storage velocity-pressure LU (HSL)
    rs_di  = 1.0_dp,   & ! real_storage for the diffuse-interface LU (HSL)
    is_di  = 1.4_dp,   & ! integer_storage for diffuse-interface LU (HSL)
    lx = 1._dp,        & ! size in x-direction
    ly = 1._dp,        & ! size in y-direction
    lz = 1._dp           ! size in z-direction


  integer ::  step, iter
  real(dp) :: cmax, cdiff, mudiff, mumax

  character(len=20) :: filename



! namelist for input of variables; read from standard input

  namelist /comppar/ nx, ny, nz, numtimesteps, rho, eta, deltat, &
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

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

! diffuse interface
  coefficients%r(101:106) = [ deltat, Mcoef, alpha, beta, kappa, rho ]


! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, &
    lz=lz, elshape=14, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )


! plot curves, surfaces and mesh

  plot_options%viewpoint=[1.,0.8,0.4]
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig', surfaces=[3,4,6] )


! problem definition of velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(presnod,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2]

  input_probdef%probnr = 1

! Dirichlet boundary conditions

  call define_essential ( mesh, input_probdef, surface1=1, surface2=6, physq=physqvel )
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

  call problem_definition ( input_probdef, mesh, problem )


! problem definition diffuse interface

  call create_input_probdef ( mesh, input_probdefdi, nvec=3, nphysq=2 )

  input_probdefdi%vec_elementdof(1)%a(:,1) = 1  ! c
  input_probdefdi%vec_elementdof(1)%a(:,2) = 1  ! mu
  input_probdefdi%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity

  input_probdefdi%physq = [1,2]
  input_probdefdi%probnr = 2


  call problem_definition ( input_probdefdi, mesh, problemdi )


! create a vector subscript for the concentration and mu without Lagr. multipl.

  call create_subscript ( mesh, problemdi, cval, physqarr=[1] )
  call create_subscript ( mesh, problemdi, muval, physqarr=[2] )


! create system vectors for velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )


! fill solution vector with essential boundary conditions


  call fill_sysvector ( mesh, problem, sol, &
    surface1=1, surface2=6, physq=physqvel, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqpress, value=0._dp )


! create system matrix of velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! build (assemble) matrix and vector for velocity/pressure problem

! stokes velocity/pressure
  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

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

! create concentration vector
!   funcnr = 1 : unstable stratified layers
!   funcnr = 2 : retraction of non-circular drop
!   funcnr = 3 : spinodal decomposition
!   funcnr = 4 : a supercritical and a subcritical nucleus
!   funcnr = 5 : spinodal decomposition (3D)

  call fill_sysvector ( mesh, problemdi, soldin, degfd=1, &
       node1=1, node2=mesh%nnodes, func=difunc, funcnr=5)

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


! time stepping

  do step = 1, numtimesteps

    print *, ' '
    print *, 'step ', step, '-------------------------------------------------------'

!   iteration scheme for (c,mu)

    iter = 0

    picard: do

      iter = iter + 1

!     build (assemble) matrix and vector for diffuse-interface problem

      call build_system ( mesh, problemdi, sysmatrixdi, sysvector=rhsdi, &
        elemsub=diffuse_interface_elem, oldvectors=oldvectors_di, &
        coefficients=coefficients )

      call check ( sysmatrixdi )

      soldi_prev%u = soldi%u

!     solve (c,mu)

      solver_options_di%real_storage=rs_di
      solver_options_di%integer_storage=is_di

      call solve_system_ma41 ( sysmatrixdi, rhsdi, soldi, &
        solver_options=solver_options_di )


      cmax  = maxval(abs(soldi%u(cval%s))) ! c maximum
      cdiff = maxval(abs(soldi%u(cval%s)-soldi_prev%u(cval%s))) ! c difference
      mumax = maxval(abs(soldi%u(muval%s))) ! mu maximum
      mudiff = maxval(abs(soldi%u(muval%s)-soldi_prev%u(muval%s))) ! mu diff

      print *, 'cmax = ', cmax, 'mumax = ', mumax
      print *, 'cdiff = ', cdiff, 'mudiff = ', mudiff

  !    if ( cdiff < epsdi * cmax .and. mudiff < epsdi * mumax ) exit picard

      if ( cdiff < epsdi * cmax ) exit picard

      if ( iter >= itermax ) then
        write(*,'(a,i0)') ' too many iterations: ', itermax
        stop
      end if

    end do picard

    soldin%u = soldi%u

!   build (assemble) vector for velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_mugradc, &
      oldvectors=oldvectors_di, physqrow=[physqvel], physqcol=[physqvel], &
      buildmatrix = .false., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve velocity/pressure problem

    solver_options_up%real_storage=rs_up
    solver_options_up%integer_storage=is_up

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
      solver_options=solver_options_up )

    write(filename,'(a,i3.3,a)') 'c_color_fill', step, '.vtk'
    call write_scalar_vtk ( mesh, problemdi,  degfd=1,&
        filename=filename, sysvector=soldi, dataname='concentration', &
        append = .true. )

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

  call delete ( coefficients )

end program diffuse_interface4
