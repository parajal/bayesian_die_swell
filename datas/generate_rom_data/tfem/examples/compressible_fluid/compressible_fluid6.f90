! Compressible viscous fluid problem in a (2D) channel or axisymmetric pipe.
! Flow with a parabolic infow profile and either an open boundary for
! entry pressure or and an imposed pressure.
! At the exit the x-direction is free (zero traction) and in the y-direction
! the velocity is set to zero.
! The constitutive equation for the extra stress is a Newtonian fluid.
! Fully-implicit time integration with Newton-Raphson iteration for the
! pressure constitutive equation.
! Optionally after the final time step do an iteration to steady-state.

program compressible_fluid6

  use tfem_m
  use hsl_ma41_m
  use stokes_elements_m
  use compressible_fluid_elements_m
  use io_utils_m
  use figplot_m
  use functions_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    coorsys = 0,        & ! 2D channel (0) or axisymmetric (1)
    funcnr = 3,         & ! function number for the entry profile
    numtimesteps = 10, & ! number of time steps
    maxnumiterations = 20, & ! maximum number of Newton-Raphson iterations
    nx=100,             & ! number of elements in x
    ny=10                 ! number of elements in y

  real(dp), parameter :: &
    Lx = 20._dp, Ly = 1._dp,  & ! size of the domain
    eta = 1._dp,  & ! viscosity
    Kbulk = 200._dp,  & ! bulk modulus
    p0 = 0._dp,        & ! reference pressure where J=1.
    deltat = 1.e0_dp,  & ! time step
    epsconf = 1e-10_dp,  & ! Newton-Raphson convergence threshold
    U = 1._dp, & ! maximum inflow velocity
    Pin = 50.0_dp ! entry pressure

  logical, parameter :: &
    open_boundary_pressure = .true., &  ! use open boundary for entry pressure
    steadylaststep = .true., & ! iterate to steady-state solution in last step
    printscreen = .true.  ! print on standard output

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: soln, solnm1, soliter
  type(sysvector_t) :: dsol, rhsd
  type(vector_t) :: velocity, pressure, vorticity, divu, Jvol
  type(subscript_t) :: velx, vely, vel, pres, solsc
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(sample_t) :: sample
  type(solver_options_ma41_t) :: solver_options

  logical :: steadystep
  integer :: step, obj_point, iter, extrastep
  real(dp) :: xs(1,2), time
  real(dp) :: epsu, epsp


! set variables in module functions_m

  U0 = U
  h0 = Ly
  L0 = Lx

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=350, ncoefr=300 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%i(23) = coorsys ! coordinate system

  coefficients%i(39) = 1 ! do not build the continuity equation

  coefficients%i(301) = 3 ! start with first order time integration
  coefficients%i(302) = 2 ! position sol_n in oldvectors
  coefficients%i(305) = 1 ! sysvector number for iteration velocity/pressure

  coefficients%r = 0
  coefficients%r(1) = eta

  coefficients%r(251) = deltat
  coefficients%r(252) = Kbulk

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%lx = Lx
  meshgen_options%ly = Ly

  call quadrilateral2d ( mesh, meshgen_options )

! one object for sampling in a single point

  xs(1,:) = [ Lx, 0._dp ]

  call add_to_mesh ( mesh, object='coordinates', nnodes=1, coor=xs )

  obj_point = mesh%nobjects

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

! center line
  call define_essential ( mesh, input_probdef, curve1=1, physq=1, degfd=[0,1] )
! exit (on the right)
  call define_essential ( mesh, input_probdef, curve1=2, physq=1, degfd=[0,1] )
! upper wall
  call define_essential ( mesh, input_probdef, curve1=3, physq=1 )
! entry (on the left)
  call define_essential ( mesh, input_probdef, curve1=4, physq=1 )

  if ( .not. open_boundary_pressure ) then
!   specify pressure on entry
    call define_essential ( mesh, input_probdef, curve1=4, physq=2 )
  end if

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for solution (excluding contraint forces)

  call create_subscript ( mesh, problem, solsc, physqarr=[physqvel,physqpress] )

! create vector subscripts for the velocity

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )

! create vector subscript for the pressure

  call create_subscript ( mesh, problem, pres, physqarr=[physqpress] )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, soln, solnm1, soliter )
  call create_sysvector ( problem, dsol, rhsd )

! fill solution vector with essential boundary conditions and initial condition
! for the pressure

  soln%u = 0

! initial pressure

  soln%u(pres%s) = p0

  call copy(soln,solnm1)
  call copy(soln,soliter)

  dsol%u = 0

! impose velocity profile on the left side of the channel in first step
  call fill_sysvector ( mesh, problem, dsol, &
    curve1=4, physq=1, degfd=1, funcnr=funcnr, func=func )

  if ( .not. open_boundary_pressure ) then
!   impose pressure on the left side of the channel in the first step
    call fill_sysvector ( mesh, problem, dsol, &
      curve1=4, physq=2, value=Pin )
  end if

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=3 )

! store solution vectors and problem structures

  oldvectors%s(1)%p => soliter
  oldvectors%s(2)%p => soln
  oldvectors%s(3)%p => solnm1

  open( unit=13, file='sample_point.out', recl=300 )

! perform a steady step?

  if ( steadylaststep ) then
    extrastep = 1
  else
    extrastep = 0
  end if

  time = 0

! time stepping

  do step = 1, numtimesteps + extrastep  ! extra step for steady-state

    if ( step <= numtimesteps ) then

!     regular time stepping

      steadystep = .false.

      time = deltat * step

      if ( printscreen) print *, 'step = ', step, ' time = ', time

      if ( step >= 2 ) then
        coefficients%i(301) = 4 ! second order time integration
      end if

    else

!     steady-state iteration

      steadystep = .true.

      if ( printscreen) print *, 'step = ', step, 'steady-state iteration'

!     exclude time-derivative in implicit_pressure_elem
      coefficients%i(306) = 1

    end if

    iter = 0

    do

      iter = iter + 1

      if ( steadystep .and. ( numtimesteps == 0 .and. iter == 2 ) .or. &
           .not. steadystep .and. ( step == 1 .and. iter == 2 ) ) then

!       impose velocity profile dsol=0 on the left side of the channel
        call fill_sysvector ( mesh, problem, dsol, &
          curve1=4, physq=1, degfd=1, value=0._dp )
        if ( .not. open_boundary_pressure ) then
!         impose pressure dsol=0 on the left side of the channel
          call fill_sysvector ( mesh, problem, dsol, &
            curve1=4, physq=2, value=0._dp )
        end if

      end if

      if ( iter > maxnumiterations ) then
        write(*,'(3(a,i0/))') &
          ' Maximum number of iterations reached = ', &
          maxnumiterations, ' step = ', step
        stop
      end if

!     build (assemble) matrix/vector

      call build_vpG

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, dsol, rhsd )

      call check ( sysmatrix )

!     solve

      call solve_system_ma41 ( sysmatrix, rhsd, dsol, &
        solver_options=solver_options  )

      soliter%u(solsc%s) = soliter%u(solsc%s) + dsol%u(solsc%s)

!     convergence test

      epsu = maxval(abs(dsol%u(vel%s)))
      epsp = maxval(abs(dsol%u(pres%s)))

      if ( printscreen) print *, iter, epsu, epsp

      if ( maxval([ epsu, epsp ]) < epsconf ) exit

    end do

!   sample in one point

    call fill_sample ( mesh, problem, sample, ndegfd=2, &
      object=obj_point, elemsub=stokes_sample_velocity, &
      coefficients=coefficients, oldvectors=oldvectors )

!   write sample to a file

    write( unit=13, fmt=* ) time, sample%u(1,:)

!   copy solution to older time step for next time step

    call copy ( soln, solnm1 )
    call copy ( soliter, soln )

    if ( printscreen ) print *, 'step = ', step
    if ( printscreen ) print *, &
       'max abs velx = ', maxval(abs(soln%u(velx%s))), &
       'max abs vely = ', maxval(abs(soln%u(vely%s))), &
       'min pres = ', minval(soln%u(pres%s)), &
       'max pres = ', maxval(soln%u(pres%s))

  end do

  close ( unit=13 )

! post-processing

  oldvectors%s(1)%p => soln

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, Jvol, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )
  call create_vector ( problem, divu, vec=3 )

  call extract_physvector ( mesh, problem, soln, velocity )

! write profile to a file for plotting with gnuplot

  call printtofile ( mesh, problem, 'vprofile.out', curve=2, &
    vector=velocity )

! derive vectors

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

! compute J directly from the pressure
  Jvol%u = exp(-(pressure%u-p0)/Kbulk)

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=7

  call derive_vector ( mesh, problem, divu, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a fig file for plotting

  plot_options%scalevector=0.05
  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  call plot_color_contour ( plot_options, mesh, problem, &
    'u_contour.fig', vector=velocity, degfd=1 )

  call plot_color_contour ( plot_options, mesh, problem, &
    'v_contour.fig', vector=velocity, degfd=2 )

  call plot_color_contour ( plot_options, mesh, problem, &
    'J_contour.fig', vector=Jvol )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

  call plot_color_contour ( plot_options, mesh, problem, &
    'divu_contour.fig', vector=divu )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( soln, solnm1, soliter, rhsd )
  call delete ( velocity, pressure, vorticity, divu, Jvol )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

contains

  subroutine build_vpG

!   build (assemble) matrix and vector for velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients )

!   pressure evolution equation
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_pressure_ce_elem, oldvectors=oldvectors, &
      physqrow=[physqpress], physqcol=[physqvel,physqpress], &
      coefficients=coefficients, addmatvec=.true. )

!   build stokes right-hand side

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_rhs_divsigma, oldvectors=oldvectors, &
      coefficients=coefficients, addmatvec=.true., &
      buildmatrix=.false., physqrow=[physqvel] )

  end subroutine build_vpG

end program compressible_fluid6
