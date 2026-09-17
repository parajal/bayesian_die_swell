! Startup of flow around a rotating cylinder in a square box.
! Compressible fluid flow.
! Fully-implicit with Newton-Raphson.
! Compute drag and torque on cylinder using reaction forces.
! Optionally integral pressure constraint imposed.
! Optionally impose pressure in a single point (P1).
! Optionally after the final time step do an iteration to steady-state.

program compressible_fluid7

  use tfem_m
  use stokes_elements_m
  use compressible_fluid_elements_m
  use hsl_ma41_m
  use io_utils_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities
    pintpl = 4,     & ! Q1 pressures
    physqvel = 1,   & ! physical quantity nr of the velocities
    physqpress = 2, & ! physical quantity nr of the pressures
    gauss = 3,      & ! 3x3 integration of quads
    gaussb = 3,     & ! 3 point integration of boundary elements
    numtimesteps = 5, & ! number of time steps
    maxnumiterations = 20    ! maximum number of Newton-Raphson iterations

  real(dp), parameter :: &
    Kmod = 1.e2_dp,    & ! compression modulus
    p0 = 0._dp,        & ! reference pressure where J=1.
    eta = 1.0_dp,      & ! fluid viscosity
!    H = 1.5_dp,         & ! (half-) height of the box (must match mesh)
    R = 1._dp,         & ! radius of the cylinder (must match mesh)
    omega = 0.2_dp       ! angular velocity of the cylinder

  real(dp), parameter :: &
    deltat = 0.3_dp,     & ! time step for BDF1/BDF2
    epsconf = 1e-10_dp,  & ! Newton-Raphson convergence threshold
    rs = 1.2_dp,         & ! real_storage for LU (HSL)
    is = 1.6_dp            ! integer_storage for LU (HSL)

  logical, parameter :: &
    steadylaststep = .true., & ! iterate to steady-state solution in last step
    pressure_constraint = .true., & ! impose pressure constraint
    check_pressure_constraint = .false., & ! check pressure constraint
    pressure_essential = .false., & ! p=p0 in point P1 (lower-left corner)
    printscreen = .true.  ! print on standard output

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: soln, solnm1, soliter
  type(sysvector_t) :: dsol, rhsd, reacf
  type(vector_t) :: velocity, pressure, vorticity, Jvol
  type(subscript_t) :: velx, vely, vel, pres, solsc, velxcylinder, velycylinder
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options

  logical :: steadystep
  integer :: step, i, iter, extrastep
  integer :: vertices(4) = [1,3,5,7]
  real(dp) :: epsu, epsp, K_dragx, K_dragy, K_torqz
  real(dp) :: int_press_domain(1), V0(1)

  timer = .false.

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=350, ncoefr=300 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gaussb ]
  coefficients%i(12:) = 0

  coefficients%i(39) = 1 ! do not build the continuity equation

  coefficients%i(301) = 3 ! start with first order time integration
  coefficients%i(302) = 2 ! position sol_n in oldvectors
  coefficients%i(305) = 1 ! sysvector number for iteration velocity/pressure
  if ( pressure_constraint ) then
    coefficients%i(307) = 1 ! constraint number of the pressure constraint
  end if
  coefficients%i(308) = 1 ! type of pressure constraint (1,2,3 or 4)


  coefficients%r = 0
  coefficients%r(1) = eta

  coefficients%r(251) = deltat
  coefficients%r(252) = Kmod
  coefficients%r(253) = p0

! read mesh

  call read_mesh ( mesh, filename='mesh.out' )

  if ( pressure_constraint ) then

!   add elementset for integral pressure constraint

    call add_to_mesh ( mesh, elementset='elements', &
                                 elements=[(i,i=1,mesh%nelem)] )
    call add_to_mesh ( mesh, elementset='nodes', elementsetnr=1 )

  end if

  call fill_mesh_parts ( mesh )

  call printinfo ( mesh, printlevel=1 )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 2         ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(vertices,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1         ! scalar

  input_probdef%physq = [1,2]
  input_probdef%probnr = 1

! define essential boundaries

! cylinder
  call define_essential ( mesh, input_probdef, curve1=13, physq=physqvel )
! walls
  call define_essential ( mesh, input_probdef, curves=[4,7,10,12], &
    physq=physqvel )

  if ( pressure_essential ) then

!   set pressure level in point P1

    call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

  end if

  if ( pressure_constraint ) then

!   constraint for integral pressure constraint

    call define_constraint ( mesh, input_probdef, &
      physq=physqpress, elementset1=1, nglobalc=1 )

  end if

  call problem_definition ( input_probdef, mesh, problem )

 ! create vector subscripts for solution (excluding contraint forces)

  call create_subscript ( mesh, problem, solsc, physqarr=[physqvel,physqpress] )

! create vector subscripts for the velocity

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )
! velocity degrees in x-direction on the cylinder
  call create_subscript ( mesh, problem, velxcylinder, physqarr=[physqvel], &
    degfd=1, curves=[13] )
! velocity degrees in y-direction on the cylinder
  call create_subscript ( mesh, problem, velycylinder, physqarr=[physqvel], &
    degfd=2, curves=[13], fillnodes=.true. )

! create vector subscript for the pressure

  call create_subscript ( mesh, problem, pres, physqarr=[physqpress] )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, soln, solnm1, soliter )
  call create_sysvector ( problem, dsol, rhsd, reacf )


! fill solution vector with essential boundary conditions

  soln%u = 0

! initial pressure

  soln%u(pres%s) = p0

  call copy(soln,solnm1)
  call copy(soln,soliter)

  dsol%u = 0

! set velocity_theta = omega*R on cylinder and pressure level in first step
  call fill_sysvector ( mesh, problem, dsol, &
    curve1=13, physq=physqvel, vfunc=vfunc, vfuncnr=1 )


! create system matrix

  if ( pressure_constraint ) then

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

  else

    call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  end if

  call create_sysmatrix_data ( sysmatrix )

  print *, 'nnz = ', sysmatrix%Suu%nnz


! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=3, nprob=1 )


! store solution vectors and problem structures

  oldvectors%s(1)%p => soliter
  oldvectors%s(2)%p => soln
  oldvectors%s(3)%p => solnm1

  print *, ' omega * R = ', omega * R
  print *


  open ( unit=13, file='iter.out', recl=300 )
  open ( unit=14, file='Kdrag.out', recl=300 )

  solver_options%real_storage=rs
  solver_options%integer_storage=is

! area of the (original) domain

  call integrate ( mesh, problem, V0, elemsub=stokes_integrate_volume, &
    coefficients=coefficients )

  coefficients%r(254) = V0(1)

! perform a steady step?

  if ( steadylaststep ) then
    extrastep = 1
  else
    extrastep = 0
  end if

  call tic

! time stepping

  do step = 1, numtimesteps + extrastep  ! extra step for steady-state

    if ( step <= numtimesteps ) then

!     regular time stepping

      steadystep = .false.

      write(13,*) 'step = ', step
      if ( printscreen) print *, 'step = ', step

      if ( step >= 2 ) then
        coefficients%i(301) = 4 ! second order time integration
      end if

    else

!     steady-state iteration

      steadystep = .true.

      write(13,*) 'step = ', step, 'steady-state iteration'
      if ( printscreen) print *, 'step = ', step, 'steady-state iteration'

!     exclude time-derivative in implicit_pressure_elem
      coefficients%i(306) = 1

    end if

    iter = 0

    do

      iter = iter + 1

      if ( steadystep ) then
        if ( numtimesteps == 0 .and. iter == 2 ) then
 !        set dsol = 0 on cylinder
          call fill_sysvector ( mesh, problem, dsol, &
            curve1=13, physq=physqvel, value=0._dp )
        end if
      else
        if ( step == 1 .and. iter == 2 ) then
 !        set dsol = 0 on cylinder
          call fill_sysvector ( mesh, problem, dsol, &
            curve1=13, physq=physqvel, value=0._dp )
        end if
      end if

      if ( iter > maxnumiterations ) then
        write(*,'(3(a,i0/))') &
          ' Maximum number of iterations reached = ', &
          maxnumiterations, ' step = ', step
        stop
      end if

!     build (assemble) matrix/vector for velocity/pressure part

      call build_vpG

      call toc ( 'build_vpG' )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, dsol, rhsd )

      call check ( sysmatrix )

!     solve

      call solve_system_ma41 ( sysmatrix, rhsd, dsol, &
        solver_options=solver_options  )

      call toc ( 'solve' )

      soliter%u(solsc%s) = soliter%u(solsc%s) + dsol%u(solsc%s)

!     convergence test

      epsu = maxval(abs(dsol%u(vel%s))) / ( omega * R )
      epsp = maxval(abs(dsol%u(pres%s))) * R / ( ( omega * R ) * eta )

      if ( printscreen) print *, iter, epsu, epsp
      write(13,*) iter, epsu, epsp

      if ( maxval([ epsu, epsp ]) < epsconf ) exit

    end do

    call toc ( 'after iteration loop' )


    if ( pressure_constraint .and. check_pressure_constraint ) then

!     check pressure constraint

      call integrate ( mesh, problem, int_press_domain, &
        elemsub=integrate_pressure_function, &
        coefficients=coefficients, oldvectors=oldvectors )

      if ( printscreen) write(*,*) 'constraint error =', int_press_domain(1)

    end if


!   reaction forces

    call reaction_forces ( problem, sysmatrix, dsol, rhsd, reacf )

    K_dragx = - sum ( reacf%u(velxcylinder%s) ) / ( eta*( omega * R ) )
    K_dragy = - sum ( reacf%u(velycylinder%s) ) / ( eta*( omega * R ) )
    K_torqz = - sum ( mesh%coor(velycylinder%nodes,1) * &
                            reacf%u(velycylinder%s) - &
                      mesh%coor(velycylinder%nodes,2) * &
                            reacf%u(velxcylinder%s) ) &
                                     / ( eta*( omega * R**2 ) )

    if ( printscreen) print *, 'Kdragx = ', K_dragx, 'Kdragy = ', K_dragy, &
      'Ktorqz = ', K_torqz
    write(14,*) step, step * deltat, K_dragx, K_dragy, K_torqz

    call toc ( 'after reaction forces' )


!   copy solution to older time step for next time step

    call copy ( soln, solnm1 )
    call copy ( soliter, soln )

    if ( step == 1 ) then
      open(unit=11, recl=600, status='replace', file='cval.out')
    else
      open(unit=11, recl=600, position='append', file='cval.out')
    end if
    write(11, fmt=*) step * deltat, maxval(soln%u(velx%s)), &
                                    maxval(soln%u(vely%s)), &
                                    maxval(soln%u(pres%s)), &
                                    minval(soln%u(pres%s))
    close(unit=11)

    if ( printscreen ) print *, 'step = ', step
    if ( printscreen ) print *, &
       'max abs velx = ', maxval(abs(soln%u(velx%s))), &
       'max abs vely = ', maxval(abs(soln%u(vely%s))), &
       'min pres = ', minval(soln%u(pres%s)), &
       'max pres = ', maxval(soln%u(pres%s))

    call toc ( 'one step' )

  end do

  call toc ( 'all steps' )


! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, Jvol, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )

  call extract_physvector ( mesh, problem, soln, velocity )

! derive vectors

  oldvectors%s(1)%p => soln

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

! compute J directly from the pressure
  Jvol%u = exp(-(pressure%u-p0)/Kmod)

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call write_vector_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='velocity', vector=velocity )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='pressure', vector=pressure, append=.true. )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='J', vector=Jvol, append=.true. )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='vorticity', vector=vorticity, append=.true. )

  call printtofile ( mesh, problem, filename='velocity_cl.out', curve=11, &
    vector=velocity )
  call printtofile ( mesh, problem, filename='pressure_cl.out', curve=12, &
    vector=pressure )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( soln, solnm1, soliter, dsol, rhsd )
  call delete ( reacf )
  call delete ( velocity, pressure, vorticity, Jvol )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

contains

  function vfunc ( n, nr, x )
    use kind_defs_m
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    vfunc = omega * [ -x(2), x(1) ] / norm2(x)

  end function vfunc

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

    if ( pressure_constraint ) then

!     build pressure contraint

!     pressure-pressure part of the Jacobian
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=pp_constraint_Jacobian_elem, coefficients=coefficients, &
        physqrow=[physqpress], physqcol=[physqpress], addmatvec=.true., &
        oldvectors=oldvectors, buildvector=.false. )

!     constraint equation
      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        constraint1=1, elemsub=constraint_pressure_int_elementset, &
        addmatvec=.true., coefficients=coefficients, &
        oldvectors=oldvectors )

    end if

  end subroutine build_vpG

end program compressible_fluid7
