#if HSL_EXTRA

! Stokes problem on a square cavity with Dirichlet boundary conditions.
! Unsteady flow. Compared to exact (manufactured) solution.
! Semi-implicit Gear/Karniadakis time discretization
! Do a real semi-implicit Euler at the first step.
! Iterative solver with block preconditioning (using both pressure mass and
! diffusion matrix). M_1^{-1}+M_2^{-1} form for approximate inverse
! Schur complement.


program navier_stokes14

  use tfem_m
  use stokes_elements_m
  use stokes_functions_m
  use inertia_elements_m
  use io_utils_m
  use figplot_m
  use timer_m
  use hsl_solve2_mi20_m


  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    itsolver = 5,       & ! iterative solver 5=BiCGSTAB, 8=GMRES
    nx=40,              & ! number of elements in x
    ny=40                 ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp           ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefp
  type(problem_t) :: problem, problemp
  type(sysmatrix_t) :: sysmatrix, sysmatrix_M1, sysmatrix_M2
  type(sysvector_t), target :: sol, soln, sol_exact, sol_hat
  type(sysvector_t) :: solnm1, rhsd, nonlinn, nonlinnm1
  type(vector_t) :: velocity, pressure, vorticity, vorticity2, pressure2
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors, oldvectors_hat
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vel, pres
  type(solver_options2_mi20_t) :: solver_options
  type(subscript_t) :: ssv, ssp
  type(prec2_mi20_t) :: prec


  integer :: step, num_picard, maxnumiterations
  integer :: numtimesteps, numeulertimesteps
  real(dp) :: H, U, rho, Re
  real(dp) :: thresh, deltat
  real(dp) :: errv, errp, tn, tnp1, maxerrv, maxerrp
  real(dp) :: gamma0=1.5_dp, alpha0=2._dp, alpha1=-0.5_dp
  real(dp) :: beta0=2.0_dp, beta1=-1.0_dp

! density
  rho = 100._dp
! number time steps:
  numtimesteps = 50
!  numtimesteps = 50
!  numtimesteps = 100
!  numtimesteps = 200
! number Euler time steps (at least one!):
  numeulertimesteps = 1
! time step
  deltat = 1.e-2_dp
!  deltat = 5.e-3_dp
!  deltat = 2.5e-3_dp

! threshold for the iteration process
  thresh = 1.e-8_dp
! maximum number of iterations allowed to obtain convergence
  maxnumiterations = 100
! number of Picard iterations, followed by Newton iterations
  num_picard = 0

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

  call create_subscript ( mesh, problem, ssv, physqarr=[1], &
    essentialpart=.false. )
  call create_subscript ( mesh, problem, ssp, physqarr=[2], &
    essentialpart=.false. )

! problem definition pressure preconditioning

  call create_input_probdef ( mesh, input_probdefp )

  input_probdefp%elementdof(1)%a(:) = [1,0,1,0,1,0,1,0,0]

  call define_essential ( mesh, input_probdefp, point=1 )

  call problem_definition ( input_probdefp, mesh, problemp )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, soln, sol_exact )
  call create_sysvector ( problem, solnm1, sol_hat )
  call create_sysvector ( problem, rhsd )
  call create_sysvector ( problem, nonlinn, nonlinnm1 )

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

! create system matrix M1 and M2

  call create_sysmatrix_structure ( sysmatrix_M1, mesh, problemp )

  call create_sysmatrix_data ( sysmatrix_M1 )

  call copy ( sysmatrix_M1, sysmatrix_M2 )

! build (assemble) matrix from elements for pressure mass matrix

  call build_system ( mesh, problemp, sysmatrix_M1, &
    elemsub=pressure_mass_elem, factormat=1.0_dp/eta, &
    coefficients=coefficients, buildvector=.false. )

  call build_system ( mesh, problemp, sysmatrix_M2, &
    elemsub=pressure_diffusion_elem, factormat=deltat/rho/gamma0, &
    coefficients=coefficients, buildvector=.false. )

  call set_solver_options ( solver_options, printlevel=2, eps_rel=1e-9_dp, &
     itsolver=itsolver, mgmres=40, cg_fixednumits=8 )

  prec%control_F%c_fail=2
  !prec%control_M2%v_iterations=2
  !prec%control%print_level=2


! set intial value of solution to zero

  soln%u = 0

! fill oldvectors

  oldvectors%s(1)%p => soln ! solution at tn
  oldvectors_hat%s(1)%p => sol_hat  ! alpha0*un+alpha1*un-1

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

!   rho un grad un term
    call build_system ( mesh, problem, sysmatrix, nonlinn, &
      elemsub=inertia_elem_ungradun, coefficients=coefficients, &
      oldvectors=oldvectors, physqrow=[physqvel], &
      buildmatrix=.false. )

!   body force vector of manufactured solution at tn+1

    coefficients%i(14) = 1

    time_m = tnp1
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_Laplace_elem, &
      coefficients=coefficients, buildmatrix=.false. )

    coefficients%i(14) = 0


    if ( step <= numeulertimesteps ) then

!     Euler forward
      rhsd%u(vel%s) = rhsd%u(vel%s) + nonlinn%u(vel%s)

!     build (assemble) viscous/pressure matrix and vector from elements
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=stokes_Laplace_elem, &
        coefficients=coefficients, buildvector=.false. )

!     instationary term ( rho du/dt ), matrix + vector

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=inertia_elem_dudt, coefficients=coefficients, &
        oldvectors=oldvectors, physqcol=[physqvel], physqrow=[physqvel], &
        addmatvec=.true. )

    else

!     combine previous time steps of the non linear term
      rhsd%u(vel%s) = rhsd%u(vel%s) + &
                        beta0 * nonlinn%u(vel%s) + beta1 * nonlinnm1%u(vel%s)

!     combine previous time steps of the old solution
      sol_hat%u(vel%s) = alpha0*soln%u(vel%s) + alpha1*solnm1%u(vel%s)

      if ( step == numeulertimesteps + 1 ) then

!       build (assemble) viscous/pressure matrix and vector from elements
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=stokes_Laplace_elem, coefficients=coefficients, &
          buildvector=.false.)

!       instationary term ( rho du/dt ), matrix + vector

        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_dudt, coefficients=coefficients, &
          oldvectors=oldvectors_hat, physqcol=[physqvel], &
          physqrow=[physqvel], addmatvec=.true., factormat=gamma0 )

      else

!       instationary term ( rho du/dt ), vector only

        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_dudt, coefficients=coefficients, &
          oldvectors=oldvectors_hat, physqcol=[physqvel], &
          physqrow=[physqvel], buildmatrix=.false., addmatvec=.true. )

      end if

    end if

!   fill solution vector with analytical solution, for the boundary conditions
    call copy ( sol_exact, sol )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    !so%real_storage = 1.2_dp

    if ( step <= numeulertimesteps ) then
      call solve_system2_mi20 ( sysmatrix, rhsd, sol, ssv%s, ssp%s, &
        sysmatrix_M1, sysmatrix_M2, solver_options=solver_options )
    else
      call solve_system2_mi20 ( sysmatrix, rhsd, sol, ssv%s, ssp%s, &
        sysmatrix_M1, sysmatrix_M2, prec=prec, solver_options=solver_options )
    end if

!   write out the errors between the analytical and the numerical solution
    errv = maxval ( sqrt ( ( sol%u(vel%s) - sol_exact%u(vel%s) )**2 ) )
    errp = maxval ( sqrt ( ( sol%u(pres%s) - sol_exact%u(pres%s) )**2 ) )
    write(*,'(/2(a,es12.4))') 'Error_v: ', errv, ' Error_p: ', errp

    maxerrv = max(errv,maxerrv)
    maxerrp = max(errp,maxerrp)

!   copy the solution and nonlinear vector to the old step(s)
    call copy ( soln, solnm1 )
    call copy ( sol, soln )
    call copy ( nonlinn, nonlinnm1 )
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

  call delete ( problem, problemp )
  call delete ( input_probdef, input_probdefp )
  call delete ( mesh )
  call delete ( sol, soln, sol_exact, sol_hat )
  call delete ( solnm1, rhsd, nonlinn, nonlinnm1 )
  call delete ( velocity, pressure, vorticity, vorticity2, pressure2 )
  call delete ( sysmatrix, sysmatrix_M1, sysmatrix_M2 )
  call delete ( coefficients )
  call delete ( oldvectors, oldvectors_hat )
  call delete ( vel, pres )
  call delete ( ssv, ssp )
  call delete ( prec )

end program navier_stokes14

#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on hsl_extra', &
    ' - add the libhsl3 library for linking', &
    ' - set preprocessing macro HSL_EXTRA in Mdefs.mk'
end
#endif

