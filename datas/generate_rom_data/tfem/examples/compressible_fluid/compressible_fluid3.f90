! Compressible viscous fluid problem in a (2D) channel or axisymmetric pipe.
! Flow with a parabolic infow profile and either an open boundary for
! entry pressure or and an imposed pressure.
! At the exit the x-direction is free (zero traction) and in the y-direction
! the velocity is set to zero.
! The constitutive equation for the extra stress is a Newtonian fluid.
! Second-order time integration for the pressure constitutive equation.

program compressible_fluid3

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
    numtimesteps = 300, & ! number of time steps
    nx=100,             & ! number of elements in x
    ny=10                 ! number of elements in y

  real(dp), parameter :: &
    Lx = 20._dp, Ly = 1._dp,  & ! size of the domain
    eta = 1._dp,  & ! viscosity
    Kbulk = 200._dp,  & ! bulk modulus
    deltat = 1.e-1_dp,  & ! time step
    U = 1._dp, & ! maximum inflow velocity
    Pin = 50.0_dp ! entry pressure

  logical, parameter :: &
    open_boundary_pressure = .true.  ! use open boundary for entry pressure

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, sol_n, sol_nm1, solhat
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity, divu
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(sample_t) :: sample

  integer :: step, obj_point
  real(dp) :: xs(1,2), time


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

  coefficients%i(301) = 1 ! first-order time integration
  coefficients%i(302) = 2 ! position sol_n in oldvectors

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

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, sol_n, sol_nm1, solhat )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions and initial condition
! for the pressure (=0)

  sol%u = 0._dp
  sol_n%u = 0._dp
  sol_nm1%u = 0._dp

! impose velocity profile On the left side of the channel
  call fill_sysvector ( mesh, problem, sol, &
    curve1=4, physq=1, degfd=1, funcnr=funcnr, func=func )
! impose pressure On the left side of the channel
  call fill_sysvector ( mesh, problem, sol, &
    curve1=4, physq=2, value=Pin )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=3 )

  oldvectors%s(1)%p => solhat
  oldvectors%s(2)%p => sol_n
  oldvectors%s(3)%p => sol_nm1

  open( unit=13, file='sample_point.out', recl=300 )

  time = 0

! time stepping

  do step = 1, numtimesteps

    if ( step > 1 ) coefficients%i(301) = 2 ! second-order time integration

    time = deltat * step

    call copy ( sol_n, sol_nm1 )
    call copy ( sol, sol_n )

    if ( coefficients%i(301) == 1 ) then
      solhat%u = sol_n%u
    else
      solhat%u = 2._dp*sol_n%u-sol_nm1%u
    end if

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
      coefficients=coefficients )

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=pressure_ce_elem, oldvectors=oldvectors, &
      physqrow=[physqpress], physqcol=[physqvel,physqpress], &
      coefficients=coefficients, addmatvec=.true. )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    call solve_system_ma41 ( sysmatrix, rhsd, sol )

!   sample in one point

    oldvectors%s(1)%p => sol

    call fill_sample ( mesh, problem, sample, ndegfd=2, &
      object=obj_point, elemsub=stokes_sample_velocity, &
      coefficients=coefficients, oldvectors=oldvectors )

!   write sample to a file

    write( unit=13, fmt=* ) time, sample%u(1,:)

    oldvectors%s(1)%p => solhat

  end do

  close ( unit=13 )

! post-processing

  oldvectors%s(1)%p => sol

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )
  call create_vector ( problem, divu, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

! write profile to a file for plotting with gnuplot

  call printtofile ( mesh, problem, 'vprofile.out', curve=2, &
    vector=velocity )

! derive vectors

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

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
    'pressure_contour.fig', vector=pressure )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

  call plot_color_contour ( plot_options, mesh, problem, &
    'divu_contour.fig', vector=divu )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, sol_n, sol_nm1, solhat, rhsd )
  call delete ( velocity, pressure, vorticity, divu )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program compressible_fluid3
