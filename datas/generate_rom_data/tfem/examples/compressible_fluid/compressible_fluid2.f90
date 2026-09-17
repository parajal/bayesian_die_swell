! Compressible viscous fluid problem in a 2D channel or axisymmetrical pipe
! A block infow/outflow profile.
! The mesh is moving with the fluid (ALE) such that the fluid is moving
! between two pistons.
! To avoid a severe pressure singularity, slip is created on the wall over a
! certain length.
! The constitutive equation for the extra stress is a Newtonian fluid.
! First-order time integration for the pressure constitutive equation.

program compressible_fluid2

  use tfem_m
  use hsl_ma41_m
  use stokes_elements_m
  use compressible_fluid_elements_m
  use functions_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    coorsys = 0,        & ! 2D channel (0) or axisymmetric (1)
    funcnr = 2,         & ! function number for the slip velocity
    numtimesteps = 400, & ! number of time steps
    nx=100,             & ! number of elements in x
    ny=10,              & ! number of elements in y
    ncx=2*nx+1,         & ! number of sample points in x
    ncy=2*ny+1            ! number of sample points in y

  real(dp), parameter :: &
    Lx = 40._dp, Ly = 1.0_dp,  & ! size of the domain
    eta = 1._dp,  & ! viscosity
    Kbulk = 200._dp,  & ! bulk modulus
    deltat = 1.e-1_dp,  & ! time step
    sl = 1.0_dp, & ! size of the slip zone
    U = 1._dp ! inflow and outflow velocity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, soln
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity, divu
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(vector_t), target :: meshvel
  type(subscriptvec_t) :: meshvelx, meshvely
  type(sample_t) :: sample

  integer :: step, i, obj_point, obj_line_x, obj_line_y
  real(dp) :: xs(1,2), xc(ncx,2), yc(ncy,2), time, v(2), p

! set variables in module functions_m

  U0 = U
  w0 = sl
  L0 = Lx

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=350, ncoefr=300 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%i(23) = coorsys ! coordinate system

  coefficients%i(48) = 1 ! use mesh velocity for ALE formulation
  coefficients%i(39) = 1 ! do not build the continuity equation
  coefficients%i(301) = 1 ! first-order time integration

  coefficients%r = 0
  coefficients%r(1) = eta

  coefficients%r(251) = deltat
  coefficients%r(252) = Kbulk

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%lx = Lx
  meshgen_options%ly = Ly

  call quadrilateral2d ( mesh, meshgen_options )

! one object for sampling in a single point

  xs(1,:) = [ Lx/2, 0._dp ]

  call add_to_mesh ( mesh, object='coordinates', nnodes=1, coor=xs )

  obj_point = mesh%nobjects

! one object for sampling on a horizontal line

  xc(:,1) = [ (i*Lx/(ncx-1), i=0,ncx-1) ]
  xc(:,2) = 0

  call add_to_mesh ( mesh, object='coordinates', nnodes=ncx, coor=xc )

  obj_line_x = mesh%nobjects

! one object for sampling on a vertical line

  yc(:,1) = Lx/2
  yc(:,2) = [ (i*Ly/(ncy-1), i=0,ncy-1) ]

  call add_to_mesh ( mesh, object='coordinates', nnodes=ncy, coor=yc )

  obj_line_y = mesh%nobjects

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
! left/right and wall
  call define_essential ( mesh, input_probdef, curve1=2, curve2=4, physq=1 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, soln )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions and initial condition
! for the pressure (=0)

  sol%u = 0

! impose slip velocity on wall
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=1, degfd=1, funcnr=funcnr, func=func )
! impose block profile at the left and right side of the channel
  call fill_sysvector ( mesh, problem, sol, &
    curve1=2, physq=1, degfd=1, value=U )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=4, physq=1, degfd=1, value=U )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! create subscript for the mesh velocity vector

  call create_subscript_vector ( mesh, problem, meshvelx, vec=physqvel, &
    degfd=1 )
  call create_subscript_vector ( mesh, problem, meshvely, vec=physqvel, &
    degfd=2 )

! create mesh velocity vector
  call create ( problem, meshvel, vec=physqvel )
  meshvel%u(meshvelx%s) = U
  meshvel%u(meshvely%s) = 0

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1, nvec=1 )

  oldvectors%v(1)%p => meshvel

  open( unit=13, file='sample_point.out', recl=300 )

  time = 0

! time stepping

  do step = 1, numtimesteps

    time = deltat * step

    call copy ( sol, soln )

!   build (assemble) matrix and vector from elements

    oldvectors%s(1)%p => soln

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

    v = sample%u(1,:)

    call fill_sample ( mesh, problem, sample, ndegfd=1, &
      object=obj_point, elemsub=stokes_sample_pressure, &
      coefficients=coefficients, oldvectors=oldvectors )

    p = sample%u(1,1)

!   write samples to a file

    write( unit=13, fmt=* ) time, v, p

  end do

  close ( unit=13 )


! post-processing

  oldvectors%s(1)%p => sol

! sample velocity on a horizontal line

  call fill_sample ( mesh, problem, sample, ndegfd=2, object=obj_line_x, &
    elemsub=stokes_sample_velocity, coefficients=coefficients, &
    oldvectors=oldvectors )

  open( unit=10, file='velocity_line_x.out', recl=300 )
  do i = 1, ncx
    write( unit=10, fmt=* ) sample%coor(i,1:2), sample%u(i,1:2)
  end do
  close( unit=10 )

! sample pressure on a horizontal line

  call fill_sample ( mesh, problem, sample, ndegfd=1, object=obj_line_x, &
    elemsub=stokes_sample_pressure, coefficients=coefficients, &
    oldvectors=oldvectors )

  open( unit=10, file='pressure_line_x.out', recl=300 )
  do i = 1, ncx
    write( unit=10, fmt=* ) sample%coor(i,1:2), sample%u(i,1)
  end do
  close( unit=10 )

! sample velocity on a vertical line

  call fill_sample ( mesh, problem, sample, ndegfd=2, object=obj_line_y, &
    elemsub=stokes_sample_velocity, coefficients=coefficients, &
    oldvectors=oldvectors )

  open( unit=10, file='velocity_line_y.out', recl=300 )
  do i = 1, ncy
    write( unit=10, fmt=* ) sample%coor(i,1:2), sample%u(i,1:2)
  end do
  close( unit=10 )

! sample pressure on a vertical line

  call fill_sample ( mesh, problem, sample, ndegfd=1, object=obj_line_y, &
    elemsub=stokes_sample_pressure, coefficients=coefficients, &
    oldvectors=oldvectors )

  open( unit=10, file='pressure_line_y.out', recl=300 )
  do i = 1, ncy
    write( unit=10, fmt=* ) sample%coor(i,1:2), sample%u(i,1)
  end do
  close( unit=10 )

! sample divergence on a vertical line

  coefficients%i(13)=7

  call fill_sample ( mesh, problem, sample, ndegfd=1, object=obj_line_y, &
    elemsub=stokes_sample_deriv, coefficients=coefficients, &
    oldvectors=oldvectors )

  open( unit=10, file='divu_line_y.out', recl=300 )
  do i = 1, ncy
    write( unit=10, fmt=* ) sample%coor(i,1:2), sample%u(i,1)
  end do
  close( unit=10 )


  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )
  call create_vector ( problem, divu, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

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
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity, divu )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program compressible_fluid2
