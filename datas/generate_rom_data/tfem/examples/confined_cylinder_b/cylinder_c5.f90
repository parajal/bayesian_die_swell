! Startup of flow around a cylinder for a Giesekus model using the conformation
! tensor formulation.
! Fully-implicit with Newton-Raphson
! Periodical boundary conditions.
! Compute drag on cylinder using reaction forces.
! Optionally include Jacobian of SUPG test function.

program cylinder_c5

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use io_utils_m
  use figplot_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,     & ! P2 velocities
    pintpl = 2,     & ! P1 pressures
    gintpl = 2,     & ! P1 gradients
    cintpl = 2,     & ! P1 conformation
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    physqc = 4,     & ! physical quantity nr of the conformation tensor (mode 1)
    gauss = 6,      & ! 6 point integration of triangles
    gaussb = 3,     & ! 3 point integration of boundary elements
    timeint1 = 8,   & ! (first-order) time integration (first step)
    timeint2 = 10,  & ! (second-order) time integration
    numtimesteps = 30, & ! number of time steps
    maxnumiterations = 20, & ! maximum number of Newton-Raphson iterations
    logc = 0,       & ! standard scheme or log transformation
    ncompc = 3,     & ! number of conformation tensor comp
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 3         ! Giesekus

  real(dp), parameter :: &
    G     = 1.0_dp,      & ! modulus
    lambda = 1.0_dp,     & ! relaxation time
    eta_p = lambda*G,  & ! polymer viscosity
    beta_s = 0.59_dp,  & ! beta viscosity parameter = eta_s/(eta_s+eta_p)
    eta_s = beta_s/(1-beta_s)*eta_p,    & ! solvent viscosity
    eta_0 = eta_s+eta_p, & ! zero-shear viscosity
    mobility = 0.01_dp,  & ! mobility parameter
    H = 2._dp,           & ! (half-) height of the channel (must match mesh)
    R = 1._dp,         & ! radius of the cylinder (must match mesh)
    U = 0.5_dp,          & ! average velocity in the channel
    flowrate = H*U         ! flow rate in (half) the channel

  real(dp), parameter :: &
    deltat = 0.4_dp,     & ! time step
    thetapar = 0.55_dp,  & ! theta parameter in the theta method
    beta_supg = 1,       & ! SUPG factor
    epsconf = 1e-12_dp,  & ! Newton-Raphson convergence threshold
    rs    = 1.2_dp,      & ! real_storage for LU (HSL)
    is    = 1.6_dp         ! integer_storage for LU (HSL)

  logical, parameter :: &
    JacobianSUPG = .true., & ! take Jacobian of SUPG test function into account
    printscreen = .true.  ! print on standard output

! definitions

  type(mesh_t) :: mesh, mesh1
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: soln, solnm1, soliter
  type(sysvector_t) :: dsol, rhsd, reacf
  type(vector_t) :: velocity, pressure, vorticity
  type(subscript_t) :: velx, vely, vel, pres, solsc, velxcylinder
  type(subscript_t), dimension(nmodes) :: cvalxx, cvalxy, cvalyy, cval
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options

  type(subscriptvec_t) :: cxx, cxy, cyy
  type(vector_t) :: ctensor

  integer :: step, i, iter, m
  integer :: vertices(3) = [1,3,5]
  real(dp) :: alpha, epsu, epsp, epsc, Wi, K_drag

  timer = .false.

! set some parameters

  alpha = eta_p   ! DEVSS parameter

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+3*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gaussb,   cintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      logc,     timeint1,   ( 0, i = 23, 150 )  &
    ]

  coefficients%i(83) = physqc ! physical quantity if first mode
  coefficients%i(84) = 3  ! storage of c tensor: all modes in sysvector
  coefficients%i(86) = 2  ! sysvector number for velocity in SUPG = soln
  coefficients%i(87) = 1  ! sysvector number for iteration conformation
  coefficients%i(88) = 2  ! sysvector number of conformation at tn

  if ( JacobianSUPG ) then
    coefficients%i(86) = 1  ! sysvector number for velocity in SUPG = soliter
    coefficients%i(91) = 1  ! Jacobian of SUPG
  end if

  coefficients%r = &
    [ eta_s,    0._dp,  0._dp,     alpha, 0._dp, &
      flowrate, 0._dp, deltat, beta_supg, 0._dp, &
      ( 0._dp, i = 11, 500 ), &
      ( G/nmodes, lambda, mobility, m=1,nmodes ) ]

  coefficients%r(10) = U ! global scaling velocity for SUPG, when needed
  coefficients%r(28) = thetapar

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
! connect curves for periodical bc in DG
  call add_to_mesh ( mesh, curve=[-8] )     ! curve 12

  call fill_mesh_parts ( mesh )

  call printinfo ( mesh, printlevel=1 )

! plot mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  plot_options%fontsize = 10

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=5+nmodes, &
     nphysq=3+nmodes )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 4  ! G
  input_probdef%vec_elementdof(1)%a(:,2) = 2         ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(vertices,4:3+nmodes) = ncompc  ! c
  input_probdef%vec_elementdof(1)%a(:,4+nmodes) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5+nmodes) = 3  ! tensor for plotting

  input_probdef%physq = [1,2,3,(3+m,m=1,nmodes)]
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

! conformation tensor (use collocation)
  do m = 1, nmodes
    call define_constraint ( mesh, input_probdef, &
      physq=physqc+m-1, curve1=4, curve2=12, discretization='collocation' )
  end do

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for solution (excluding contraint forces)

  call create_subscript ( mesh, problem, solsc, &
    physqarr=[physqgrad,physqvel,physqpress,(physqc+m-1,m=1,nmodes)] )

! create vector subscripts for the velocity

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )
! velocity degrees in x-direction on the cylinder
  call create_subscript ( mesh, problem, velxcylinder, physqarr=[physqvel], &
    degfd=1, curves=[9] )

! create vector subscript for the pressure

  call create_subscript ( mesh, problem, pres, physqarr=[physqpress] )

! create vector subscript for the conformation tensor

  do m = 1, nmodes
    call create_subscript ( mesh, problem, cval(m), physqarr=[physqc+m-1] )
    call create_subscript ( mesh, problem, cvalxx(m), physqarr=[physqc+m-1], &
      degfd=1 )
    call create_subscript ( mesh, problem, cvalxy(m), physqarr=[physqc+m-1], &
      degfd=2 )
    call create_subscript ( mesh, problem, cvalyy(m), physqarr=[physqc+m-1], &
      degfd=3 )
  end do
  call create_subscript ( mesh, problem, cxx, degfd=1, vec=5+nmodes )
  call create_subscript ( mesh, problem, cxy, degfd=2, vec=5+nmodes )
  call create_subscript ( mesh, problem, cyy, degfd=3, vec=5+nmodes )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, soln, solnm1, soliter )
  call create_sysvector ( problem, dsol, rhsd, reacf )

! create a vector for conformation
! tensor for post processing

  call create_vector ( problem, ctensor, vec=5+nmodes )

! fill solution vector with essential boundary conditions

  soln%u = 0

! initial solution conformation tensor

  do m = 1, nmodes
    if ( logc == 0 ) then ! standard
      soln%u(cvalxx(m)%s) = 1 ! initial cxx
      soln%u(cvalxy(m)%s) = 0 ! initial cxy
      soln%u(cvalyy(m)%s) = 1 ! initial cyy
    else if ( logc == 1 ) then ! log scheme
      soln%u(cval(m)%s) = 0 ! initial s
    end if
  end do

  call copy(soln,solnm1)
  call copy(soln,soliter)
  call copy(soln,dsol)

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=3, nprob=1 )


! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => soliter
  oldvectors_ve%s(2)%p => soln
  oldvectors_ve%s(3)%p => solnm1
  oldvectors_ve%p(1)%p => problem  ! needed for implicit_ce_supg_elem


! Weisenberg number
  Wi = lambda * U / R

  print *, ' U = ', U, 'Wi = ', Wi
  print *


  open ( unit=13, file='iter.out', recl=300 )
  open ( unit=14, file='Kdrag.out', recl=300 )

  solver_options%real_storage=rs
  solver_options%integer_storage=is

  call tic

! time stepping

  do step = 1, numtimesteps

    write(13,*) 'step = ', step
    if ( printscreen) print *, 'step = ', step

    if ( step == 1 ) then
      coefficients%i(78) = 1 ! check for zero velocity in SUPG
    else
      coefficients%i(78) = 0
    end if

    if ( step >= 2 ) then
      coefficients%i(22) = timeint2
    end if

    iter = 0

    do

      iter = iter + 1

      if ( step == 1 .and. iter == 2 ) then
        coefficients%r(6) = 0  ! set flowrate in dsol=0
      end if

      if ( iter > maxnumiterations ) then
        write(*,'(3(a,i0/))') &
          ' Maximum number of iterations reached = ', &
          maxnumiterations, ' step = ', step
        stop
      end if

!     build (assemble) matrix/vector for gradient/velocity/pressure part

      call build_vpG

      call toc ( 'build_vpG' )

!     build (assemble) matrix/vector for conformation part

      do m = 1, nmodes

        coefficients%i(85) = m   ! set mode number

        call build_c

      end do

      call toc ( 'build_c' )

!     set to zero off-diagonal conformation blocks

      call build_offdiag_c

      call toc ( 'build_offdiag_c' )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, dsol, rhsd )

      call check ( sysmatrix )

!     solve

      call solve_system_ma41 ( sysmatrix, rhsd, dsol, &
        solver_options=solver_options  )

      call toc ( 'solve' )

      soliter%u(solsc%s) = soliter%u(solsc%s) + dsol%u(solsc%s)

!     convergence test

      epsu = maxval(abs(dsol%u(vel%s))) / U
      epsp = maxval(abs(dsol%u(pres%s))) * R / ( U * eta_0 )
      epsc = maxval(abs(dsol%u(cval(1)%s)))

      if ( printscreen) print *, iter, epsu, epsp, epsc
      write(13,*) iter, epsu, epsp, epsc

      if ( maxval([ epsu, epsp, epsc ]) < epsconf ) exit

    end do

    call toc ( 'after iteration loop' )

!   reaction forces

    call reaction_forces ( problem, sysmatrix, dsol, rhsd, reacf )

    K_drag = - 2 * sum ( reacf%u(velxcylinder%s) ) / ( eta_0*U )

    if ( printscreen) print *, 'Kdrag = ', K_drag
    write(14,*) step, step * deltat, K_drag

    call toc ( 'after reaction forces' )


!   copy solution to older time step for next time step

    call copy ( soln, solnm1 )
    call copy ( soliter, soln )

!   write max and mean values of conformation tensor to a file

    call derive_vector ( mesh, problem, ctensor, &
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
                                    sum(ctensor%u(cxx%s))/size(cxx%s), &
                                    sum(ctensor%u(cxy%s))/size(cxy%s), &
                                    sum(ctensor%u(cyy%s))/size(cyy%s), &
                                    maxval(soln%u(velx%s)), &
                                    maxval(soln%u(vely%s))
    close(unit=11)

    if ( printscreen ) &
         print *, 'step = ', step, 'max cxx = ', maxval(ctensor%u(cxx%s))

    call toc ( 'one step' )

  end do

  call toc ( 'all steps' )


! post-processing

  call create_vector ( problem, velocity, physq=2 )
  call create_vector ( problem, pressure, vec=4+nmodes )
  call create_vector ( problem, vorticity, vec=4+nmodes )

  call extract_physvector ( mesh, problem, soln, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! derive vectors

  oldvectors%s(1)%p => soln

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

  call derive_vector ( mesh, problem, ctensor, &
    elemsub=deriv_conformation_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problem, filename='c.vtk', &
    dataname='conformation_tensor', vector=ctensor )

  call printtofile ( mesh, problem, filename='conformation_cl.out', curve=11, &
    vector=ctensor )

  call delete ( ctensor )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( soln, solnm1, soliter, dsol, rhsd )
  call delete ( reacf )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

  call delete ( oldvectors_ve )
  do m = 1, nmodes
    call delete ( cvalxx(m), cvalxy(m), cvalyy(m), cval(m) )
  end do
  call delete ( solsc )
  call delete ( cxx, cxy, cyy )

contains

  subroutine build_vpG

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients, &
      physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress] )

!   stokes add (q,divu) to rhs of mass balance
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_rhs_divu, oldvectors=oldvectors_ve, &
      coefficients=coefficients, addmatvec=.true., &
      physqrow=[physqpress], physqcol=[physqvel], buildmatrix=.false. )

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

!   build DEVSS-G right-hand side

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_rhs_div, oldvectors=oldvectors_ve, &
      coefficients=coefficients, addmatvec=.true., &
      buildmatrix=.false., physqrow=[physqvel] )

!   build stokes right-hand side

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_rhs_divsigma, oldvectors=oldvectors_ve, &
      coefficients=coefficients, addmatvec=.true., &
      buildmatrix=.false., physqrow=[physqvel] )

!   build -div(tau) part in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_divtau, oldvectors=oldvectors_ve, &
      physqrow=[physqvel], physqcol=[(physqc+m-1,m=1,nmodes)], &
      addmatvec=.true., coefficients=coefficients )

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


  subroutine build_c

!   build (assemble) matrix and vector for conformation tensor mode m

    if ( JacobianSUPG ) then

!     velocity + diagonal block
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
        physqrow=[physqc+m-1], physqcol=[physqvel,physqc+m-1], &
        oldvectors=oldvectors_ve, addmatvec=.true. )

    else

!     diagonal block
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
        physqrow=[physqc+m-1], physqcol=[physqc+m-1], &
        oldvectors=oldvectors_ve, addmatvec=.true. )

    end if

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_ce_vel_supg_elem, coefficients=coefficients, &
      physqrow=[physqc+m-1], physqcol=[physqgrad,physqvel], &
      oldvectors=oldvectors_ve, addmatvec=.true., buildvector=.false. )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      physqrow=[physqgrad,physqpress], physqcol=[physqc+m-1], &
      buildvector=.false., zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      physqrow=[physqc+m-1], physqcol=[physqpress], &
      buildvector=.false., zeromatvec=.true. )

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problem, sysmatrix, &
      sysvector=rhsd, constraint1=3+m, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

  end subroutine build_c

  subroutine build_offdiag_c

    integer :: i, j

    if ( nmodes > 1 ) then
      do i = 1, nmodes
        do j = 1, nmodes
          if ( i == j ) cycle
          call build_system ( mesh, problem, sysmatrix, rhsd, &
            addmatvec=.true., physqrow=[physqc+i-1], physqcol=[physqc+j-1], &
            buildvector=.false., zeromatvec=.true. )
        end do
      end do
    end if

  end subroutine build_offdiag_c

end program cylinder_c5
