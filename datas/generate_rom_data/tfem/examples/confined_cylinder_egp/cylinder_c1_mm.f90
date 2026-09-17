! Startup of flow around a cylinder for an incompressible EGP model using
! the conformation tensor formulation. Explicit div tau term.
! Problem 1: Flow around a confined cylinder.
!   Periodical boundary conditions.
!   Stress-explicit formulation of the momentum balance.
! Problem 2: Conformation tensor
! Multi-mode EGP model with the option to use a global von Mises stress.

program cylinder_c1_mm

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma57_m
  use hsl_ma41_m
  use io_utils_m
  use figplot_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,         & ! P2 velocities
    pintpl = 2,         & ! P1 pressures
    gintpl = 2,         & ! P1 gradients
    cintpl = 2,         & ! P1 conformation
    physqgrad = 1,      & ! physical quantity nr of the gradients
    physqvel = 2,       & ! physical quantity nr of the velocities
    physqpress = 3,     & ! physical quantity nr of the pressures
    gauss = 6,          & ! 6 point integration of triangles
    gaussb = 3,         & ! 3 point integration of boundary elements
    timeint1 = 1,       & ! first-order time integration, first time step
    timeint2 = 5,       & ! second-order time integration after first time step
    numtimesteps = 200, & ! number of time steps
    logc = 1,           & ! standard scheme or log transformation
    ncompc = 4,         & ! number of conformation tensor comp
    nmodes = 2,         & ! number of modes
    startm = 501,       & ! start of material model data
    model = 23,         & ! EGP incompressible
    alam_model = 4        ! adapted lambda model: 4: Eyring, 5: Ree-Eyring

  real(dp), parameter, dimension(nmodes) :: &
!   TEST case: by splitting G of one mode to G/2 per mode for two modes
!   and using a global von Mises gives the same results as the one mode problem.
    G =        [ 2.0_dp, 2.0_dp ], & ! modulus
    lambda =   [ 1.0_dp, 1.0_dp ], & ! relaxation time
    tau_ref =  [ 1.0_dp, 1.0_dp ], & ! reference von Mises
    tau_ref1 = [ 1.0_dp, 1.0_dp ], & ! reference von Mises 1 RE
    tau_ref2 = [ 4.0_dp, 4.0_dp ], & ! reference von Mises 2 RE
    f1 =       [ 0.6_dp, 0.6_dp ]    ! weight factor for first term RE

  real(dp), parameter :: &
    eta_p = sum(lambda*G), & ! polymer viscosity
    beta_s = 0.59_dp,    & ! beta viscosity parameter = eta_s/(eta_s+eta_p)
    eta_s = beta_s/(1-beta_s)*eta_p,    & ! solvent viscosity
!    eta_0 = eta_s+eta_p, & ! zero-shear viscosity
    H = 2._dp,           & ! (half-) height of the channel (must match mesh)
    R = 1._dp,           & ! radius of the cylinder (must match mesh)
    U = 0.5_dp,          & ! average velocity in the channel
    flowrate = H*U         ! flow rate in (half) the channel

  real(dp), parameter :: &
    deltat = 2.e-2_dp,  & ! time step
    beta_supg = 1,      & ! SUPG factor
    rs_gup = 1.2_dp, & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp, & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.0_dp,  & ! real_storage for the conformation tensor LU (HSL)
    is_c  = 1.6_dp     ! integer_storage for conformation tensor LU (HSL)

! definitions

  type(mesh_t) :: mesh, mesh1
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, solm1
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
  type(subscript_t) :: velx, vely
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma57_t) :: solver_options_u
  type(lu_ma57_t) :: lu_u

  type(input_probdef_t) :: input_probdefc
  type(problem_t), target :: problemc
  type(sysmatrix_t) :: sysmatrixc
  type(oldvectors_t) :: oldvectors_ve
  type(subscript_t) :: cval
  type(subscriptvec_t) :: cxx, cxy, cyy, czz
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc, solcm1
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc
  type(vector_t) :: ctensor, vonmises
  type(lu_ma41_t) :: luc
  type(solver_options_ma41_t) :: solver_options_c

  integer :: icomp, step, i, nn, m
  integer :: vertices(3) = [1,3,5]
  real(dp) :: alpha, Wi
  logical :: printscreen = .true.


  timer = .false.

! set some parameters

  alpha = eta_p  ! DEVSS parameter


! fill coefficients

  if ( alam_model == 4 ) then
    nn = 1
  else if ( alam_model == 5 ) then
    nn = 3
  end if

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+(2+nn)*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gaussb,   cintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      logc,  timeint1,   ( 0, i = 23, 150 )  &
    ]

  coefficients%i(80) = alam_model ! adapted lambda
  coefficients%i(92) = 1  ! global von Mises

  coefficients%r(1:500) = &
    [ eta_s,    0._dp, 0._dp,  alpha,     0._dp, &
      flowrate, 0._dp, deltat, beta_supg, 0._dp, &
      ( 0._dp, i = 11, 500 ) &
    ]

  if ( alam_model == 4 ) then
    coefficients%r(501:500+3*nmodes) = &
      [ (G(m), lambda(m), tau_ref(m), m=1,nmodes)]
  else if ( alam_model == 5 ) then
    coefficients%r(501:500+5*nmodes) = &
      [ (G(m), lambda(m), tau_ref1(m), f1(m), tau_ref2(m), m=1,nmodes) ]
  end if

! read mesh

  call read_mesh_gmsh ( mesh1, filename='confined_cylinder1.msh', ndim=2 )
  call mesh_convert ( mesh1, mesh, remove_isolated_nodes=.true. )
  call delete ( mesh1 )

! cylinder
  call add_to_mesh ( mesh, curve=[-2] ) ! curve 9
! top
  call add_to_mesh ( mesh, curve=[5,6,7] ) ! curve 10
! centerline+cylinder
  call add_to_mesh ( mesh, curve=[1,-2,3] ) ! curve 11
! connect curves for periodical bc
  call add_to_mesh ( mesh, curve=[-8] )     ! curve 12

  call fill_mesh_parts ( mesh )

  call printinfo ( mesh, printlevel=1 )

! plot mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  plot_options%fontsize = 10

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=3 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 4  ! G
  input_probdef%vec_elementdof(1)%a(:,2) = 2         ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

! define essential boundaries

! cylinder
  call define_essential ( mesh, input_probdef, curve1=9, physq=physqvel )
! top (wall)
  call define_essential ( mesh, input_probdef, curve1=10, physq=physqvel )
! bottom (center line)
  call define_essential ( mesh, input_probdef, curve1=1, &
    physq=physqvel, degfd=[0,1] )
  call define_essential ( mesh, input_probdef, curve1=3, &
    physq=physqvel, degfd=[0,1] )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=12, nglobalc=1 )

! constraints for periodical boundary conditions

! velocities (use weak connection)
  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=4, curve2=12, discretization='weak', &
    elementdof=[2,0,2] )

! gradients (use collocation)
  call define_constraint ( mesh, input_probdef, &
    physq=physqgrad, curve1=4, curve2=12, discretization='collocation' )

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for the velocity for post processing

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix/vector for gradient/velocity/pressure problem

  call build_vpG

! solve initial gradient/velocity/pressure consistent with initial conditions
! and keep decomposition of the matrix
! This generates a linear profile.

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options_u%real_storage=rs_gup
  solver_options_u%integer_storage=is_gup

  call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
    solver_options=solver_options_u  )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=3, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a(vertices,1) = 1 ! c
  input_probdefc%vec_elementdof(1)%a(:,2) = 1 ! scalar for plotting
  input_probdefc%vec_elementdof(1)%a(:,3) = 4 ! tensor for plotting

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

! constraint for periodical boundary conditions of the b tensor
  call define_constraint ( mesh, input_probdefc, curve1=4, curve2=12, &
    discretization='collocation' )

  call problem_definition ( input_probdefc, mesh, problemc )


! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemc, cval, physqarr=[1] )
  call create_subscript ( mesh, problemc, cxx, degfd=1, vec=3 )
  call create_subscript ( mesh, problemc, cxy, degfd=2, vec=3 )
  call create_subscript ( mesh, problemc, cyy, degfd=3, vec=3 )
  call create_subscript ( mesh, problemc, czz, degfd=4, vec=3 )

! create a vector for conformation tensor for post processing

  call create_vector ( problemc, ctensor, vec=3 )

! create system vectors (solution and right-hand side) for conformation and
! initialize vectors

  call create ( problemc, solc, solcm1, rhsc )

! initial solution

  if ( logc == 0 ) then ! standard
    do m = 1, nmodes
      solc(1,m)%u = 1 ! initial cxx
      solc(2,m)%u = 0 ! initial cxy
      solc(3,m)%u = 1 ! initial cyy
      solc(4,m)%u = 1 ! initial czz
    end do
  else if ( logc == 1 ) then ! log scheme
    do m = 1, nmodes
      solc(1,m)%u = 0 ! initial sxx
      solc(2,m)%u = 0 ! initial sxy
      solc(3,m)%u = 0 ! initial syy
      solc(4,m)%u = 0 ! initial szz
    end do
  end if

! create system matrix for conformation tensor problem

  call create_sysmatrix_structure_base ( sysmatrixc, mesh, problemc )
  call create_sysmatrix_structure_constraint ( sysmatrixc, mesh, problemc )
  call finalize_sysmatrix_structure ( sysmatrixc )

  call create_sysmatrix_data ( sysmatrixc )

! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec2=2, nprob=2 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s(2)%p => solm1
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%s2(2)%p => solcm1
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc


! Weisenberg number
  Wi = sum( G*lambda**2 )/eta_p * U / R

  print *, ' U = ', U, 'Wi = ', Wi
  print *


! time stepping

  call tic

  do step = 1, numtimesteps

    if ( step == 2 ) then
!     change time integration scheme at the second time step
      coefficients%i(22) = timeint2
    end if

!   build (assemble) matrix and vector for conformation tensor problem

    call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
      elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
      coefficients=coefficients )

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problemc, sysmatrixc, &
      m2sysvector=rhsc, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

    call check ( sysmatrixc )


!   copy previous c-tensor solution to older time step

    call copy ( solc, solcm1 )


!   solve c tensor and keep LU decomposition in the loop over modes/components

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do m = 1, nmodes
      do icomp = 1, ncompc
        call solve_system_ma41 ( sysmatrixc, rhsc(icomp,m), solc(icomp,m), &
          luc, solver_options=solver_options_c  )
      end do
    end do

    call delete ( luc )  ! remove LU decomposition and rebuild next time step


!   build (assemble) vector for gradient/velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_divtau, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      buildmatrix=.false., coefficients=coefficients )

!   imposed flow rate (right-hand side only)

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      buildmatrix=.false., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   copy previous gradient/velocity/pressure solution to older time step

    call copy ( sol, solm1 )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
      solver_options=solver_options_u  )


!   write max and mean values of conformation tensor to a file

    coefficients%i(28) = 1 ! mode to use

    call derive_vector ( mesh, problemc, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    if ( step == 1 ) then
      open(unit=11, recl=300, status='replace', file='cval.out')
    else
      open(unit=11, recl=300, position='append', file='cval.out')
    end if
    write(11, fmt=*) step * deltat, maxval(ctensor%u(cxx%s)), &
                                    maxval(ctensor%u(cxy%s)), &
                                    maxval(ctensor%u(cyy%s)), &
                                    maxval(ctensor%u(czz%s)), &
                                    sum(ctensor%u(cxx%s))/size(cxx%s), &
                                    sum(ctensor%u(cxy%s))/size(cxy%s), &
                                    sum(ctensor%u(cyy%s))/size(cyy%s), &
                                    sum(ctensor%u(czz%s))/size(czz%s)
    close(unit=11)

    if ( printscreen ) &
        print *, 'step = ', step, 'max cxx = ', maxval(ctensor%u(cxx%s))

    call toc ( 'one step' )

  end do

  call toc ( 'all steps' )


! post-processing

  call create_vector ( problem, velocity, physq=2 )
  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, vorticity, vec=4 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! derive vectors

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call write_vector_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='velocity', vector=velocity )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='pressure', vector=pressure, append=.true. )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='vorticity', vector=vorticity, append=.true. )

  call printtofile ( mesh, problem, filename='velocity_cl.out', curve=11, &
    vector=velocity )
  call printtofile ( mesh, problem, filename='pressure_cl.out', curve=11, &
    vector=pressure )

! conformation tensor

  coefficients%i(28) = 1 ! mode to use

  call derive_vector ( mesh, problemc, ctensor, &
    elemsub=deriv_conformation_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemc, filename='c.vtk', &
    dataname='conformation_tensor', vector=ctensor, assume33=.true. )

  call printtofile ( mesh, problemc, filename='conformation_cl.out', curve=11, &
    vector=ctensor )

  call delete ( ctensor )

! von Mises stress

  call create_vector ( problemc, vonmises, vec=2 )

  coefficients%i(13)=1 ! von Mises

  if ( coefficients%i(92) == 1 ) then
!   global von Mises
    coefficients%i(28)=0
!    coefficients%i(26) = 1  ! mode1
!    coefficients%i(27) = 2  ! mode2
  else
    coefficients%i(28)=1 ! mode number
  end if

  call derive_vector ( mesh, problemc, vonmises, &
    elemsub=deriv_viscoelastic_stress_scalar, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_scalar_vtk ( mesh, problemc, filename='c.vtk', &
    dataname='von_mises', vector=vonmises, append=.true. )

  call delete ( vonmises )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( lu_u )

  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( solc, solcm1 )
  call delete ( rhsc )
  call delete ( sysmatrixc )
  call delete ( oldvectors_ve )
  call delete ( cval )
  call delete ( cxx, cxy, cyy, czz )

contains

  subroutine build_vpG

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients, &
      physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress] )

!   DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_elem, coefficients=coefficients, addmatvec=.true., &
      physqrow=[physqgrad,physqvel], physqcol=[physqgrad,physqvel] )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqgrad], physqcol=[physqpress], &
      zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqpress], physqcol=[physqgrad], &
      zeromatvec=.true. )

!   flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients )

!   periodical condition on velocities

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=stokes_constr_elem_conn, addmatvec=.true., &
      coefficients=coefficients )

!   periodical condition on gradients

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=3, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients )

  end subroutine build_vpG

end program cylinder_c1_mm
