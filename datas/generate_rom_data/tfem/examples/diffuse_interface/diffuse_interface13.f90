! axisymmetric diffuse-interface problem of a droplet on a surface
! Stokes problem with mu grad c term for the surface tension
! constant viscosity
! variable (static) contact angle
! first order time integration
! Newton-Raphson iterations for the non-linear term in the
! Cahn-Hilliard equation

program diffuse_interface13

  use tfem_m
  use math_defs_m
  use functions_m
  use subs_m
  use diffuse_interface_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    coorsys = 1,    & ! axisymmetric
    uintpl = 8,     & ! Q2 velocities and (c,mu)
    pintpl = 4,     & ! Q1 pressures
    physqvel = 1,   & ! physical quantity nr of the velocities
    physqpress = 2, & ! physical quantity nr of the pressures
    gauss = 3         ! 3 point integration


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
  type(coefficients_t) :: coefficients
  type(subscript_t) :: cval, muval
  type(solver_options_ma57_t) :: solver_options_up
  type(solver_options_ma41_t) :: solver_options_di
  type(subscript_t) :: pval, velx, vely


! variables

  integer :: &
    nz=80,                  & ! number of elements in x
    nr=80,                  & ! number of elements in y
    numtimesteps = 100,     & ! number of time steps
    write_vtk_every = 10,   & ! write a vtk every ... timesteps
    itermax = 20              ! maximum number of Newton iterations

  real(dp) :: &
    eta = 1.0_dp,      & ! viscosity
    deltat = 1.e-3_dp, & ! time step
    rho = 1._dp,       & ! density in diffuse-interface method
    alpha = 1.0_dp,    & ! parameter in the diffuse-interface method
    beta = 1.0_dp,     & ! parameter in the diffuse-interface method
    Mcoef = 0.01_dp,   & ! parameter in the diffuse-interface method
    kappa = 4.e-4_dp , & ! parameter in the diffuse-interface method
    theta_c = 70_dp,   & ! contact angle measured through drop (degrees)
    epsdi = 1.e-12_dp, & ! accuracy in Newton iteration
    drop_z0 = 1.0_dp,  & ! initial z-coordinate of center of drop
    drop_r = 0.4_dp,   & ! initial radius of the drop
    time = 0.0_dp,     & ! initial time
    rs_up = 1.4_dp,    & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.4_dp,    & ! integer_storage velocity-pressure LU (HSL)
    rs_di  = 1.4_dp,   & ! real_storage for the diffuse-interface LU (HSL)
    is_di  = 1.4_dp      ! integer_storage for diffuse-interface LU (HSL)

  integer ::  step=0, iter, ipost=0
  real(dp) :: cmax, cdiff, phi_di
  real(dp) :: surf_ten, c_B, x_radius(2), x_height(2)

  character(len=20) :: filename


! pass value to functions_m

  xi = sqrt(kappa/alpha)
  allocate ( c(2) )
  c(1) = drop_z0
  c(2) = 0.0_dp
  r = drop_r


! determine the wetting potential phi

  c_B = sqrt(alpha/beta)
  surf_ten = rho*(2*sqrt(2.0_dp)/3) * (kappa * c_B**2) / xi

  phi_di = (surf_ten/rho) * cos(theta_c*2*pi/360._dp)/(2*(c_B-(c_B**3)/3.0_dp))


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=200, ncoefr=150 )

! Stokes
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0
  coefficients%i(23) = coorsys

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

! diffuse interface
  coefficients%r(101:106) = [ deltat, Mcoef, alpha, beta, kappa, rho ]
  coefficients%r(108) = phi_di
  coefficients%i(152) = 2 ! Newton-Raphson iteration for di equation


! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nz
  meshgen_options%ny = nr

  call quadrilateral2d ( mesh, meshgen_options )

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

  call define_essential ( mesh, input_probdef, curve1=1, degfd=[0,1], &
    physq=physqvel )
  call define_essential ( mesh, input_probdef, curve1=2, curve2=4, &
    physq=physqvel )
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

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

  call problem_definition ( input_probdefdi, mesh, problemdi )


! create a vector subscript for c and mu

  call create_subscript ( mesh, problemdi, cval, physqarr=[1] )
  call create_subscript ( mesh, problemdi, muval, physqarr=[2] )


! create system vectors for velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )

! create subscripts
  call create_subscript ( mesh, problem, velx, physqarr=[1], &
    degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[1], &
    degfd=2 )
  call create_subscript ( mesh, problem, pval, physqarr=[2] )


! fill solution vector with essential boundary conditions

! set velocity = 0
  call fill_sysvector ( mesh, problem, sol, &
    curve1=2, curve2=4, physq=physqvel, degfd=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, physq=physqvel, degfd=2, value=0._dp )
! set pressure level = 0 in lower left corner
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqpress, value=0._dp )



! create system vectors (solution and right-hand side) for diffuse interface and
! initialize vectors plus a small random perturbation.

  call create_sysvector ( problemdi, soldi, rhsdi, soldin, soldi_prev )


! sample the initial c-field and then find the initial mu-field

  call find_mu_init

  soldin%u = soldi%u


! create system matrix for diffuse-interface problem

  call create_sysmatrix_structure_base ( sysmatrixdi, mesh, problemdi )
  call finalize_sysmatrix_structure ( sysmatrixdi )

  call create_sysmatrix_data ( sysmatrixdi )


! create system matrix of velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! create the structure oldvectors_di

  call create_oldvectors ( oldvectors_di, nsysvec=3, nprob=2 )
    oldvectors_di%s(1)%p => sol
    oldvectors_di%s(2)%p => soldin
    oldvectors_di%s(3)%p => soldi
    oldvectors_di%p(1)%p => problem
    oldvectors_di%p(2)%p => problemdi


! build (assemble) matrix and vector for velocity/pressure problem

! stokes velocity/pressure

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

! build (assemble) vector for velocity/pressure problem

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_mugradc, &
    oldvectors=oldvectors_di, physqrow=[physqvel], physqcol=[physqvel], &
    buildmatrix = .false., coefficients=coefficients )

  call check ( sysmatrix )


! solve initial velocity/pressure consistent with initial conditions
! and keep decomposition of the matrix

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options_up%real_storage=rs_up
  solver_options_up%integer_storage=is_up

  call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
    solver_options=solver_options_up )


! postprocessing

  call postprocessing

! find the coordinates of c=0 on the line (z,r) = (1,0) to (1,1)
  call find_c0 ( mesh, problemdi, soldi, x1=[1.0_dp, 0.0_dp], &
    x2=[1.0_dp,1.0_dp], x_c0=x_radius, tol=1.e-12_dp )
! find the coordinates of c=0 on the line (z,r) = (0,0) to (1,0)
  call find_c0 ( mesh, problemdi, soldi, x1=[0.0_dp, 0.0_dp], &
    x2=[1.0_dp,0.0_dp], x_c0=x_height, tol=1.e-12_dp )

  open (unit = 100, file='output_drop.out')
  write(100, '(i5,3es16.8)') step, time, x_radius(2),1.0_dp-x_height(1)
  close(unit=100)

  print *,'Initial drop radius is ',x_radius(2)
  print *,'Initial drop height is ',1.0_dp-x_height(1)


! time stepping

  do step = 1, numtimesteps

    time = time + deltat

    print *,'==============================================================='
    print *,'step = ',step
    print *,'==============================================================='


!   iteration scheme for (c,mu)

    iter = 0

    newton: do

      iter = iter + 1

!     build (assemble) matrix and vector for diffuse-interface problem

      call build_system ( mesh, problemdi, sysmatrixdi,  &
        sysvector=rhsdi, elemsub=diffuse_interface_elem1, &
        oldvectors=oldvectors_di, coefficients=coefficients )

!     add the boundary conditon for the contact angle

      call add_boundary_elements ( mesh, problemdi, rhsdi, curve=2, &
        sysmatrix=sysmatrixdi, elemsub=di_natboun_curve, &
        coefficients=coefficients, oldvectors=oldvectors_di )

      call check ( sysmatrixdi )

      soldi_prev%u = soldi%u


!     solve (c,mu)

      solver_options_di%real_storage=rs_di
      solver_options_di%integer_storage=is_di

      call solve_system_ma41 ( sysmatrixdi, rhsdi, soldi, &
        solver_options=solver_options_di )


      cmax  = maxval(abs(soldi%u(cval%s))) ! c maximum
      cdiff = maxval(abs(soldi%u(cval%s)-soldi_prev%u(cval%s))) ! c difference

      print *,'iter = ',iter, '   cdiff = ', cdiff

      if ( cdiff < epsdi * cmax ) exit newton

      if ( iter >= itermax ) then
        write(*,'(a,i0)') ' too many iterations: ', itermax
        stop
      end if

    end do newton


!   copy di solution vector

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


!   postprocessing

    if ( mod(step,write_vtk_every) == 0 ) call postprocessing


!   find the coordinates of c=0 on the line (z,r) = (1,0) to (1,1)
    call find_c0 ( mesh, problemdi, soldi, x1=[1.0_dp, 0.0_dp], &
      x2=[1.0_dp,1.0_dp], x_c0=x_radius, tol=1.e-12_dp )
!   find the coordinates of c=0 on the line (z,r) = (0,0) to (1,0)
    call find_c0 ( mesh, problemdi, soldi, x1=[0.0_dp, 0.0_dp], &
      x2=[1.0_dp,0.0_dp], x_c0=x_height, tol=1.e-12_dp )

    open (unit = 100, file='output_drop.out', position='APPEND')
    write(100, '(i5,3es16.8)') step, time, x_radius(2),1.0_dp-x_height(1)
    close(unit=100)

    print *,'drop radius is ',x_radius(2)
    print *,'drop height is ',1.0_dp-x_height(1)

  end do


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

  contains


! find the initual mu field

  subroutine find_mu_init

    type(input_probdef_t) :: input_probdefdi_init
    type(problem_t), target :: problemdi_init
    type(sysvector_t), target :: soldi_init, rhsdi_init
    type(subscript_t) :: cval_init, muval_init
    type(sysmatrix_t) :: sysmatrixdi_init
    type(oldvectors_t) :: oldvectors_di_init

!   problem definition for the diffuse-interface problem
    call create_input_probdef ( mesh, input_probdefdi_init, nvec=2, nphysq=2 )

    input_probdefdi_init%vec_elementdof(1)%a =   &
      reshape ( [ 1,1,1,1,1,1,1,1,1,       &  ! c
                   1,1,1,1,1,1,1,1,1 ],    &  ! mu
                   [9,2] )

    input_probdefdi_init%physq = [1,2]
    input_probdefdi_init%probnr = 3

!   c is known at t=0
    call define_essential ( mesh, input_probdefdi_init, elgroup1=1, &
      physq=1 )
    call problem_definition ( input_probdefdi_init, mesh, problemdi_init )

!   create system vectors (solution and right-hand side) for di problem
    call create ( problemdi_init, soldi_init, rhsdi_init )

    call fill_sysvector ( mesh, problemdi_init, soldi_init, physq=1, &
      node1=1, node2=mesh%nnodes, func=difunc, funcnr=1)

!   create subscripts for c and mu
    call create_subscript ( mesh, problemdi_init, cval_init, physqarr=[1] )
    call create_subscript ( mesh, problemdi_init, muval_init, physqarr=[2] )

!   set velocity to zero to find the inital mu-field
    sol%u = 0

!   create oldvectors
    call create_oldvectors ( oldvectors_di_init, nsysvec=5, nprob=2, nvec=1 )
      oldvectors_di_init%s(1)%p => sol
      oldvectors_di_init%s(2)%p => soldi_init
      oldvectors_di_init%s(3)%p => soldi_init
      oldvectors_di_init%p(1)%p => problem
      oldvectors_di_init%p(2)%p => problemdi_init

!   create system matrix for the di problem
    call create_sysmatrix_structure_base ( sysmatrixdi_init, mesh, problemdi_init )
    call finalize_sysmatrix_structure ( sysmatrixdi_init )
    call create_sysmatrix_data ( sysmatrixdi_init )

!   build (assemble) matrix and vector for diffuse-interface problem
    call build_system ( mesh, problemdi_init, sysmatrixdi_init, &
      sysvector=rhsdi_init, elemsub=diffuse_interface_elem1, &
        oldvectors=oldvectors_di_init, coefficients=coefficients )

!   add the boundary conditon for the contact angle
    call add_boundary_elements ( mesh, problemdi_init, sysvector=rhsdi_init, &
      curve=2, sysmatrix=sysmatrixdi_init, elemsub=di_natboun_curve, &
      coefficients=coefficients, oldvectors=oldvectors_di_init )

    call check ( sysmatrixdi_init )

    call add_effect_of_essential_to_rhs ( problemdi_init, sysmatrixdi_init, &
      soldi_init, rhsdi_init )

!   solve (c,mu)
    solver_options_di%real_storage=rs_di
    solver_options_di%integer_storage=is_di

    call solve_system_ma41 ( sysmatrixdi_init, rhsdi_init, soldi_init, &
      solver_options=solver_options_di )

    soldi%u(cval%s) = soldi_init%u(cval_init%s)
    soldi%u(muval%s) = soldi_init%u(muval_init%s)

    call delete ( input_probdefdi_init )
    call delete ( problemdi_init )
    call delete ( soldi_init, rhsdi_init )
    call delete ( cval_init, muval_init )

  end subroutine find_mu_init


! write data to vtks

  subroutine postprocessing

    type(vector_t) :: pressure

    call create_vector ( problem, pressure, vec=3 )

!   derive the pressure in all nodes
    call derive_vector ( mesh, problem, pressure, &
      elemsub=stokes_pressure, coefficients=coefficients, &
      oldvectors=oldvectors_di )

    write(filename,'(a,i4.4,a)') 'flow', ipost, '.vtk'
    call write_scalar_vtk ( mesh, problem, vector=pressure, &
      dataname='modified_pressure',  filename=filename )

    write(filename,'(a,i4.4,a)') 'flow', ipost, '.vtk'
    call write_vector_vtk ( mesh, problem, filename=filename, &
      dataname='velocity_vector', sysvector=sol, physq=physqvel, &
      append=.true. )

    write(filename,'(a,i4.4,a)') 'di', ipost, '.vtk'
    call write_scalar_vtk ( mesh, problemdi, sysvector=soldi, physq=1, &
      dataname='c',  filename=filename )

    call write_scalar_vtk ( mesh, problemdi, sysvector=soldi, physq=2, &
      dataname='mu',  filename=filename, append=.true. )

    call delete ( pressure )

    ipost = ipost+1

  end subroutine postprocessing

end program diffuse_interface13
