! Generalized Newtonian problem on a unit square.
! Periodical boundary conditions.
! Poiseuille flow. Collocation of nodes.
! Similar to generalized_newtonian1 but now Newton-Raphson iteration
! instead of Picard iteration.

program generalized_newtonian_newton_raphson

  use tfem_m
  use hsl_ma57_m
  use generalized_stokes_elements_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    gnmodel = 2,        & ! generalized Newtonian model
    nPicard=2,          & ! number of Picard iterations before Newton-Raphson
    itermax=100,        & ! maximum interations
    nx=2,               & ! number of elements in x
    ny=25                 ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp,       & ! viscosity
    epsvd = 1.e-14,  & ! maximum velocity difference for iteration
    epspd = 1.e-13,  & ! maximum pressure difference for iteration
    flowrate = 2._dp

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, soln
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vel, pres
  type(solver_options_ma57_t) :: solver_options

  integer :: iter
  real(dp) :: vd, pd

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=300, ncoefr=250 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gauss ]
  coefficients%i(12:) = 0

  coefficients%i(254) = gnmodel

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(6) = flowrate

  coefficients%r(201) = 1._dp   ! m (for Power-law only)
  coefficients%r(202) = 0.7_dp  ! n
  coefficients%r(203) = 1._dp   ! eta_0
  coefficients%r(204) = 0.0_dp  ! eta_inf
  coefficients%r(205) = 1._dp   ! lambda
  coefficients%r(206) = 0.7_dp   ! a (for Carreau-Yasuda)

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,1,0,0,    &  ! pressure
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                   [9,3] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=3, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, nglobalc=1 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, curve2=5, discretization='collocation', exclude=3 )

  call problem_definition ( input_probdef, mesh, problem )

! define vector subscripts for direct manipulation of sysvector data

! all velocities
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )
! pressures
  call create_subscript ( mesh, problem, pres, physqarr=[physqpress] )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, soln )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

  soln%u = 0 ! initialize old solution to zero

! Create vector defined per element

  call create ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => soln

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! starting vector: Stokes flow

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, coefficients=coefficients, &
    oldvectors=oldvectors )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficients  )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%real_storage=1.5

  call solve_system_ma57 ( sysmatrix, rhsd, sol, &
    solver_options=solver_options )

  call copy ( sol, soln )



! Start Newton-Raphson iteration

  iter = 0

  iterate: do

    iter = iter + 1

    if ( iter <= nPicard ) then
      coefficients%i(253) = 0  ! Picard=0
    else
      coefficients%i(253) = 1  ! Newton-Raphson=1
      coefficients%r(6) = 0    ! set flowrate of increment to zero
    end if

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=generalized_newtonian_elem, coefficients=coefficients, &
      oldvectors=oldvectors )

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients )

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients  )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    solver_options%real_storage=10.5

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options )

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
!      stop
      exit iterate
    end if

  end do iterate


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, soln )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( vel, pres )

end program generalized_newtonian_newton_raphson
