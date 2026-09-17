! Navier-Stokes problem on a square cavity with Dirichlet boundary conditions.
! Unsteady flow. Compared to exact (manufactured) solution.
! Crank-Nicolson time discretization with a predictor for the convection
! velocity at tn+1 based on extrapolation from two previous time steps:
!   u_hat=2*un-un-1
! First time step: predictor is un.

program navier_stokes6

  use tfem_m
  use hsl_ma41_m
  use stokes_elements_m
  use stokes_functions_m
  use inertia_elements_m
  use io_utils_m
  use figplot_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=40,              & ! number of elements in x
    ny=40                 ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp           ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, sol_hat, soln, sol_exact
  type(sysvector_t) :: rhsd, solnm1
  type(vector_t) :: velocity, pressure, vorticity, vorticity2, pressure2
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors, oldvectors_hat
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vel, pres

  integer :: step
  integer :: numtimesteps, numeulertimesteps
  real(dp) :: H, U, rho, Re
  real(dp) :: deltat
  real(dp) :: errv, errp, tn, tnp1, maxerrv, maxerrp

! density
  rho = 100._dp
! number time steps:
!  numtimesteps = 100
  numtimesteps = 200
! number Euler time steps (at least one!):
  numeulertimesteps = 1
! time step
  deltat = 1.e-2_dp
!  deltat = 5.e-3_dp

! mesh width and heigth:
  H = 1._dp
! velocity
  U = 1._dp
! Write out the Reynolds number
  Re=rho * H * U / eta
  write(*,'(a,es12.4)') 'Re = ', Re

! fill coefficients
  call create_coefficients ( coefficients, ncoefi=250, ncoefr=200 )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(8) = deltat
  coefficients%r(151) = rho

  coefficients%vfunc => vfunc

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh
  meshgen_options%elshape = 6
  meshgen_options%ox = 0._dp
  meshgen_options%oy = 0._dp
  meshgen_options%lx = H
  meshgen_options%ly = H
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, sol_hat, soln, sol_exact )
  call create_sysvector ( problem, solnm1 )
  call create_sysvector ( problem, rhsd )

! define vector subscripts for direct manipulation of sysvector data

! velocities
  call create_subscript ( mesh, problem, vel, physqarr=[1] )
! pressures
  call create_subscript ( mesh, problem, pres, physqarr=[2] )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )
  call create_oldvectors ( oldvectors_hat, nsysvec=1 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! set intial value of solution to zero

  soln%u = 0

! fill oldvectors

  oldvectors%s(1)%p => soln ! solution at tn
  oldvectors_hat%s(1)%p => sol_hat ! estimated solution at tn+1

! start time integration

  call tic

  maxerrv = 0
  maxerrp = 0

  tn = 0

  integration: do step = 1, numtimesteps

    tnp1 = tn + deltat

    write(*,'(/a,i0,a,es12.4,a/)') ' ** time step = ', step, &
                                   ' time = ', tnp1, ' ** '

!   constants in functions

    eta_m  = coefficients%r(1)
    rho_m  = coefficients%r(151)

!   fill exact solution at time tn+1

    time_m = tnp1

    call fill_sysvector ( mesh, problem, sol_exact, physq=physqvel, degfd=1, &
      elgroup1=1, func=func, funcnr=1 )
    call fill_sysvector ( mesh, problem, sol_exact, physq=physqvel, degfd=2, &
      elgroup1=1, func=func, funcnr=2 )
    call fill_sysvector ( mesh, problem, sol_exact, physq=physqpress, &
      elgroup1=1, func=func, funcnr=3 )

!   viscous term and pressure at tn
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_rhs_divsigma, coefficients=coefficients, &
      oldvectors=oldvectors, physqrow=[physqvel], buildmatrix=.false., &
      factorvec=0.5_dp )

!   rho un grad un term
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=inertia_elem_ungradun, coefficients=coefficients, &
      oldvectors=oldvectors, physqrow=[physqvel], &
      addmatvec=.true., buildmatrix=.false., factorvec=0.5_dp )

!   body force vector of manufactured solution at tn and tn+1

    coefficients%i(14) = 1

    time_m = tn
    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
      coefficients=coefficients, addmatvec=.true., buildmatrix=.false., &
      factorvec=0.5_dp )
    time_m = tnp1
    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
      coefficients=coefficients, addmatvec=.true., buildmatrix=.false., &
      factorvec=0.5_dp )

    coefficients%i(14) = 0

!   build (assemble) viscous/pressure matrix and vector from elements
    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
      coefficients=coefficients, buildvector=.false., factormat=0.5_dp )

!   predictor of convection velocity at tn+1

    if ( step <= numeulertimesteps ) then
!     predictor = velocity at previous time step
      sol_hat%u(vel%s) = soln%u(vel%s)
    else
!     compute predictor by extrapolation
      sol_hat%u(vel%s) = 2 * soln%u(vel%s) - solnm1%u(vel%s)
    end if

!   u^hat*grad u^n+1
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=inertia_elem_conv_picard, coefficients=coefficients, &
      oldvectors=oldvectors_hat, physqcol=[physqvel], &
      physqrow=[physqvel], addmatvec=.true., &
      factormat=0.5_dp, buildvector=.false. )

!   instationary term ( rho du/dt )

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=inertia_elem_dudt, coefficients=coefficients, &
      oldvectors=oldvectors, physqcol=[physqvel], physqrow=[physqvel], &
      addmatvec=.true. )

!   fill solution vector with analytical solution, for the boundary conditions
    call copy ( sol_exact, sol )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    call solve_system_ma41 ( sysmatrix, rhsd, sol )

!   write out the errors between the analytical and the numerical solution
    errv = maxval ( sqrt ( ( sol%u(vel%s) - sol_exact%u(vel%s) )**2 ) )
    errp = maxval ( sqrt ( ( sol%u(pres%s) - sol_exact%u(pres%s) )**2 ) )
    write(*,'(/2(a,es12.4))') 'Error_v: ', errv, ' Error_p: ', errp

    maxerrv = max(errv,maxerrv)
    maxerrp = max(errp,maxerrp)

!   copy the solution to the old step
    call copy ( soln, solnm1 )
    call copy ( sol, soln )
    tn = tnp1

  end do integration

  call toc

  write(*,'(/2(a,es12.4))') 'Maximum error_v: ', maxerrv, &
                           ' Maximum error_p: ', maxerrp

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, pressure2, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )
  call create_vector ( problem, vorticity2, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  oldvectors%s(1)%p => sol_exact

  call derive_vector ( mesh, problem, pressure2, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity2, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a fig file for plotting

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  call plot_color_fill ( plot_options, mesh, problem, 'pressure_color.fig', &
    vector=pressure )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour_exact.fig', vector=pressure2 )

  call plot_color_fill ( plot_options, mesh, problem, 'vorticity_color.fig', &
    vector=vorticity )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour_exact.fig', vector=vorticity2 )

! write binary file for reading by streamfunction

  open(unit=10,form='unformatted',file='velocity_bin.out')

  write(10) velocity%u

  close(unit=10)

! write vtk data
  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='cavity_inertia.vtk' )

  call write_vector_vtk ( mesh, problem, filename='cavity_inertia.vtk', &
    dataname='velocity_vector', sysvector=sol, append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, sol_hat, soln, sol_exact)
  call delete ( solnm1)
  call delete ( rhsd )
  call delete ( velocity, pressure, vorticity, vorticity2, pressure2 )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors, oldvectors_hat )
  call delete ( vel, pres )

end program navier_stokes6
