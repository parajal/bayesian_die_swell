! Navier-Stokes problem on a square cavity with Dirichlet boundary conditions.
! Unsteady flow. Compared to exact (manufactured) solution.
! Crank-Nicolson time discretization


module stokes_functions_m

  use math_defs_m

  implicit none

  save

  real(dp) :: time_m, rho_m, eta_m

contains

  function func ( nr, xin )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: xin
    real(dp) :: func

    real(dp) :: x, y, t

    x = xin(1)
    y = xin(2)
    t = time_m

    select case ( nr )

    case (1)

!     velocity solution in x-direction

      func = - sin ( pi * x ) * cos ( pi * y ** 2 ) * y * sin ( pi * t )

    case (2)

!     velocity solution in y-direction

      func = cos( pi * x ) * sin ( pi * y ** 2 ) / 2 * sin( pi * t )

    case (3)

!     the pressure solution

      func = - cos ( pi * x ) * sin ( pi * y / 2 ) * sin( pi * t )

    case default

      write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
      stop

    end select

  end function func

  function vfunc ( n, nr, xin )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: xin
    real(dp), dimension(n) :: vfunc

    real(dp) :: x, y, t, rho, eta

    x = xin(1)
    y = xin(2)
    t = time_m
    rho = rho_m
    eta = eta_m

    select case ( nr )

    case (1)

!     rhs solution

      vfunc(1) = 1._dp/(4._dp ) * ( sin( pi * x ) * &
        ( -4._dp * pi * y * rho * cos ( pi * t ) * cos( pi * y**2 ) + &
        sin ( pi * t ) * ( 4._dp * pi * ( -1._dp *pi * ( y + 4._dp * y**3) * &
        eta * cos( pi * y**2) + sin( pi * y / 2._dp ) - 6._dp * y * eta * &
        sin( pi * y**2 )) - rho * cos( pi * x ) * sin( pi * t) * &
        ( -4._dp * pi * y**2 + sin( 2._dp * pi * y**2)))))

      vfunc(2) = 1._dp/(4._dp ) * ( 2._dp * pi * cos ( pi * x ) * ( &
        -1._dp *( cos( pi * y / 2._dp ) + 2._dp * eta * cos( pi * y**2 ) ) * &
        sin ( pi * t ) + ( rho * cos( pi * t ) + pi * ( 1._dp  + 4._dp * y**2)*&
        eta * sin( pi * t )) * sin( pi * y**2) ) + pi * y * rho * &
        sin( pi * t ) **2 * sin( 2._dp * pi * y**2 ) )

    case default

      write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
      stop

    end select

  end function vfunc

end module stokes_functions_m

program inertia_CN

  use tfem_m
  use hsl_ma41_m
  use stokes_elements_m
  use stokes_functions_m
  use inertia_elements_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=5,              & ! number of elements in x
    ny=5                 ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp           ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol_iter, soln
  type(sysvector_t) :: sol, rhsd, rhsdn, sol_exact
  type(oldvectors_t) :: oldvectors, oldvectors_iter
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vel, pres

  integer :: step, iterstep, num_picard, maxnumiterations
  integer :: numtimesteps
  real(dp) :: H, U, rho, Re
  real(dp) :: diffv, thresh, diffp, deltat
  real(dp) :: errv, errp, tn, tnp1, maxerrv, maxerrp

! density
  rho = 100._dp
! number time steps:
  numtimesteps = 2
!  numtimesteps = 200
! time step
  deltat = 1.e-2_dp

! threshold for the iteration process
  thresh = 1.e-8_dp
! maximum number of iterations allowed to obtain convergence
  maxnumiterations = 100
! number of Picard iterations, followed by Newton iterations
  num_picard = 0

! mesh width and heigth:
  H = 1._dp
! velocity of the lid
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

! create mesh
  meshgen_options%elshape = 6
  meshgen_options%ox = 0._dp
  meshgen_options%oy = 0._dp
  meshgen_options%lx = H
  meshgen_options%ly = H
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

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

  call create_sysvector ( problem, sol, sol_iter, soln, sol_exact )
  call create_sysvector ( problem, rhsd, rhsdn )

! define vector subscripts for direct manipulation of sysvector data

! velocities
  call create_subscript ( mesh, problem, vel, physqarr=[1] )
! pressures
  call create_subscript ( mesh, problem, pres, physqarr=[2] )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )
  call create_oldvectors ( oldvectors_iter, nsysvec=1 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! set intial value of solution to zero

  soln%u = 0
  sol_iter%u = soln%u

! fill oldvectors

  oldvectors%s(1)%p => soln  ! solution at tn
  oldvectors_iter%s(1)%p => sol_iter ! solution at end of previous iteration

! start time integration

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

!   fill right-hand side that is constant during iteration

!   viscous term and pressure at tn
    call build_system ( mesh, problem, sysmatrix, rhsdn, &
      elemsub=stokes_rhs_divsigma, coefficients=coefficients, &
      oldvectors=oldvectors, physqrow=[physqvel], buildmatrix=.false., &
      factorvec=0.5_dp )

!   rho un grad un term
    call build_system ( mesh, problem, sysmatrix, rhsdn, &
      elemsub=inertia_elem_ungradun, coefficients=coefficients, &
      oldvectors=oldvectors, physqrow=[physqvel], &
      addmatvec=.true., buildmatrix=.false., factorvec=0.5_dp )

!   body force vector of manufactured solution at tn and tn+1

    coefficients%i(14) = 1

    time_m = tn
    call build_system ( mesh, problem, sysmatrix, rhsdn, elemsub=stokes_elem, &
      coefficients=coefficients, addmatvec=.true., buildmatrix=.false., &
      factorvec=0.5_dp )
    time_m = tnp1
    call build_system ( mesh, problem, sysmatrix, rhsdn, elemsub=stokes_elem, &
      coefficients=coefficients, addmatvec=.true., buildmatrix=.false., &
      factorvec=0.5_dp )

    coefficients%i(14) = 0

!   start iteration for nonlinear term

    iterstep = 0

!   fill solution vector with analytical solution, for the boundary conditions
    call copy ( sol_exact, sol )

    iterate: do

      iterstep = iterstep + 1

!     copy stored rhsdn
      call copy ( rhsdn, rhsd )

!     build (assemble) viscous/pressure matrix and vector from elements
      call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
        coefficients=coefficients, buildvector=.false., factormat=0.5_dp )

!     first num_picard steps Picard iteration, then Newton iteration for
!     the term un+1.grad un+1
      if ( iterstep <= num_picard ) then
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_picard, coefficients=coefficients, &
          oldvectors=oldvectors_iter, physqcol=[physqvel], &
          physqrow=[physqvel], addmatvec=.true., &
          factormat=0.5_dp, factorvec=0.5_dp )
      else
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_newton, coefficients=coefficients, &
          oldvectors=oldvectors_iter, physqcol=[physqvel], &
          physqrow=[physqvel], addmatvec=.true., &
          factormat=0.5_dp, factorvec=0.5_dp )
      end if

!     instationary term ( rho du/dt )
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=inertia_elem_dudt, coefficients=coefficients, &
        oldvectors=oldvectors, physqcol=[physqvel], physqrow=[physqvel], &
        addmatvec=.true. )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

      call solve_system_ma41 ( sysmatrix, rhsd, sol )

!     compute the difference between iteration steps
      diffv = maxval ( sqrt ( ( sol%u(vel%s) - sol_iter%u(vel%s) )** 2 ) )
      diffp = maxval ( sqrt ( ( sol%u(pres%s) - sol_iter%u(pres%s) )** 2 ) )

      write(*,'(a,i0,2(a,es24.16))') &
        'iterstep = ', iterstep, ' diffv = ', diffv, ' diffp = ', diffp

      call copy ( sol, sol_iter )

      if ( diffv < thresh ) exit iterate

      if ( iterstep == maxnumiterations ) then
        write(*,'(a,i0)') &
          'Maximum number of iterations reached = ', maxnumiterations
        stop
      end if

    end do iterate

!   write out the errors between the analytical and the numerical solution
    errv = maxval ( sqrt ( ( sol%u(vel%s) - sol_exact%u(vel%s) )**2 ) )
    errp = maxval ( sqrt ( ( sol%u(pres%s) - sol_exact%u(pres%s) )**2 ) )
    write(*,'(/2(a,es24.16))') 'Error_v: ', errv, ' Error_p: ', errp

    maxerrv = max(errv,maxerrv)
    maxerrp = max(errp,maxerrp)

!   copy the converged solution to the old step
    call copy ( sol, soln )
    tn = tnp1

  end do integration

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, sol_iter, soln, sol_exact)
  call delete ( rhsd, rhsdn )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors, oldvectors_iter )
  call delete ( vel, pres )

end program inertia_CN
