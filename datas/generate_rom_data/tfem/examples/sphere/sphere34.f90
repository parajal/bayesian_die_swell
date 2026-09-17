! Axisymmetrical generalized Newtonian problem: a sphere moving with
! velocity U on the center line of a cylindrical container with a bottom
! The sphere is moving with velocity U. At the ends fluid is at rest
! (zero velocity) representing a bottom.
! Dirichlet boundary conditions for the velocities are used.
! Picard iteration followed by Newton-Raphson iteration.
! All models can be tested, including the yield models.

program sphere34

  use tfem_m
  use hsl_ma57_m
  use generalized_stokes_elements_m
  use figplot_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    gnmodel = 0,        & ! generalized Newtonian model
    yieldmodel = 2,     & ! yield model 1: Bingham, 2: Regularized Bingham.
    nPicard=20,         & ! number of Picard iterations before Newton-Raphson
    itermax=100,        & ! maximum interations
    coorsys = 1           ! axisymmetric coordinate system

  real(dp), parameter :: &
    eta = 1._dp,  & ! viscosity
    rs_gup = 10.5_dp,  & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.4_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    epsvd = 1.e-12,    & ! maximum velocity difference for iteration
    epspd = 1.e-7,     & ! maximum pressure difference for iteration
    U = 1._dp            ! velocity of the sphere

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, soln
  type(sysvector_t) :: rhsd, reacf
  type(vector_t) :: velocity, pressure, vorticity, viscosity, vmstress, &
    gammadot
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma57_t) :: solver_options
  type(subscript_t) :: velzparticle, vel, pres

! variables

  integer :: crv, crvpart, iter
  real(dp) :: vd, pd

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=300, ncoefr=250 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0
  coefficients%i(23) = coorsys

  coefficients%i(254) = gnmodel
  coefficients%i(265) = yieldmodel

  coefficients%r = 0
  coefficients%r(1) = eta

  coefficients%r(201) = 1._dp   ! m (for Power-law only)
  coefficients%r(202) = 0.6_dp  ! n
  coefficients%r(203) = 1._dp   ! eta_0
  coefficients%r(204) = 0.0_dp  ! eta_inf
  coefficients%r(205) = 1._dp   ! lambda
  coefficients%r(206) = 0.7_dp  ! a (for Carreau-Yasuda)
  coefficients%r(207) = 1.0_dp  ! eta_const
  coefficients%r(209) = 10.0_dp ! tau_y
  coefficients%r(210) = 100._dp ! regularizing parameter for yield

  call write_coefficients ( coefficients, filename='coefficients.out' )

! read mesh

  call read_mesh_gmsh ( mesh, filename='mesh1.msh', ndim=2 )

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

! curve for printing data on the centerline+sphere

  call add_to_mesh ( mesh, curve=[1,2,3,4,5,6,7,8], newnr=crv )

! curve for the sphere

  call add_to_mesh ( mesh, curve=[3,4,5,6], newnr=crvpart )

  call fill_mesh_parts ( mesh )

! plot mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  plot_options%fontsize = 10

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a = &
      reshape ( [2,2,2,2,2,2,2,2,2,   &  ! velocity
                 1,0,1,0,1,0,1,0,0,   &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                   [9,3] )

  input_probdef%physq = [1,2]

! define essential boundaries

! center line
  call define_essential ( mesh, input_probdef, curve1=1, curve2=2, &
    physq=physqvel, degfd=[0,1] )
  call define_essential ( mesh, input_probdef, curve1=7, curve2=8, &
    physq=physqvel, degfd=[0,1] )
! sphere
  call define_essential ( mesh, input_probdef, curve1=crvpart, physq=physqvel )
! ends and cylinder wall
  call define_essential ( mesh, input_probdef, curve1=9, curve2=14, &
    physq=physqvel )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, soln, rhsd, reacf )

! define vector subscripts for direct manipulation of sysvector data

! all velocities
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )
! pressures
  call create_subscript ( mesh, problem, pres, physqarr=[physqpress])
! velocity degrees in z-direction on the particle
  call create_subscript ( mesh, problem, velzparticle, physqarr=[physqvel], &
    degfd=1, curves=[crvpart] )


! Create oldvectors

  call create ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => soln

! fill solution vector with essential boundary conditions

  sol%u = 0

  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, curve2=6, physq=physqvel, degfd=1, value=U )


! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )


! starting vector: Stokes flow

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%real_storage=rs_gup
  solver_options%integer_storage=is_gup

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

  call copy ( sol, soln )


! Start Newton-Raphson iteration

  iter = 0

  iterate: do

    iter = iter + 1

    if ( iter <= nPicard ) then
      coefficients%i(253) = 0  ! Picard=0
    else
      coefficients%i(253) = 1  ! Newton-Raphson=1
!     set velocity of increment to zero
      call fill_sysvector ( mesh, problem, sol, &
        curve1=3, curve2=6, physq=physqvel, degfd=1, value=0._dp )
    end if

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=generalized_newtonian_elem, coefficients=coefficients, &
      oldvectors=oldvectors )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    solver_options%real_storage=rs_gup
    solver_options%integer_storage=is_gup

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options )

    call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

    if ( coefficients%i(253) == 1 ) then
!     Newton-Raphson
      vd = maxval(abs(sol%u(vel%s)))               ! velocity difference
      sol%u(vel%s) = sol%u(vel%s) + soln%u(vel%s)  ! velocity update
    else
      vd = maxval(abs(sol%u(vel%s)-soln%u(vel%s)))  ! velocity difference
    end if

    pd = maxval(abs(sol%u(pres%s)-soln%u(pres%s)))  ! pressure difference

    print *, iter, vd , pd

    call copy ( sol, soln )

    if ( vd < epsvd .and. pd < epspd ) exit iterate

    if ( iter >= itermax ) then
      write(*,'(a,i0)') ' too many iterations: ', itermax
      exit iterate
    end if

  end do iterate

  print *,'Fp   = ', sum ( reacf%u(velzparticle%s) )

  call delete ( sysmatrix )

  call create_vector ( problem, velocity, physq=physqvel )
  call extract_physvector ( mesh, problem, sol, velocity )

! post processing

  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )
  call create_vector ( problem, viscosity, vec=3 )
  call create_vector ( problem, vmstress, vec=3 )
  call create_vector ( problem, gammadot, vec=3 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=8

  call derive_vector ( mesh, problem, gammadot, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, viscosity, &
    elemsub=generalized_newtonian_viscosity, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5 ! tau = eta * gammadot = von Mises stress

  call derive_vector ( mesh, problem, vmstress, &
    elemsub=generalized_newtonian_stress, &
    coefficients=coefficients, oldvectors=oldvectors )


! write to vtk

  call write_vector_vtk ( mesh, problem, vector=velocity, &
    dataname='velocity', filename='sphere34.vtk' )

  call write_scalar_vtk ( mesh, problem, vector=pressure, &
    dataname='pressure', filename='sphere34.vtk', append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=vorticity, &
    dataname='vorticity', filename='sphere34.vtk', append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    dataname='gammadot', filename='sphere34.vtk', append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=viscosity, &
    dataname='viscosity', filename='sphere34.vtk', append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=vmstress, &
    dataname='vmstress', filename='sphere34.vtk', append=.true. )


! write binary file for reading by streamfunction

  open(unit=10,form='unformatted',file='velocity_bin.out')

  write(10) velocity%u

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, soln, reacf, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( gammadot, viscosity, vmstress )
  call delete ( coefficients )
  call delete ( oldvectors )

end program sphere34
