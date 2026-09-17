! Stokes problem on a unit square. Periodical boundary conditions.
! Poiseuille flow. Collocation of nodes.
! Varying visosity eta using a separate element to fill the Gauss
! point values of eta according to a generalized Newtonian model.
! Picard iteration.

program generalized_stokes_channel

  use tfem_m
  use hsl_ma57_m
  use generalized_stokes_elements_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    gnmodel = 2,        & ! generalized Newtonian model
    itermax=100,        & ! maximum interations
    nx=2,               & ! number of elements in x
    ny=20                 ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp,       & ! viscosity
    epsvd = 1.e-1_dp,  & ! maximum velocity difference for iteration
    epspd = 1.e-1_dp,  & ! maximum pressure difference for iteration
    flowrate = 2._dp

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, soln
  type(sysvector_t) :: rhsd
  type(vector_t) :: stress
  type(elvector_t), target :: elvector
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

  coefficients%i(251) = 3  ! compute eta using Gauss values in separate routine
  coefficients%i(254) = gnmodel

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(6) = flowrate

  coefficients%r(201) = 1._dp   ! m (for Power-law only)
  coefficients%r(202) = 0.3_dp  ! n
  coefficients%r(203) = 1._dp   ! eta_0
  coefficients%r(204) = 0.0_dp  ! eta_inf
  coefficients%r(205) = 1._dp   ! lambda

!  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

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

  call create_elvector ( mesh, elvector, nreal1d=1 )

  call create ( oldvectors, nsysvec=1, nelvec=1 )

  oldvectors%s(1)%p => soln
  oldvectors%e(1)%p => elvector

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! Start Picard iteration

  iter = 0

  iterate: do

    iter = iter + 1

!   build (assemble) matrix and vector from elements

    if ( gnmodel == 1 .and. iter == 1 ) then

!     for power-law model use Stokes for first iteration

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=stokes_elem, coefficients=coefficients, &
        oldvectors=oldvectors )

    else

!     store Gauss values of eta

      call loop_over_elements ( mesh, problem, elemsub=fill_eta_gn, &
        coefficients=coefficients, oldvectors=oldvectors )

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=generalized_stokes_elem, coefficients=coefficients, &
        oldvectors=oldvectors )

    end if

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

    vd = maxval(abs(sol%u(vel%s))-soln%u(vel%s))  ! velocity difference
    pd = maxval(abs(sol%u(pres%s)-soln%u(pres%s)))  ! pressure difference

    print *, iter, vd , pd

    call copy ( sol, soln )

    if ( vd < epsvd .and. pd < epspd ) exit iterate

    if ( iter >= itermax ) then
      write(*,'(a,i0)') ' too many iterations: ', itermax
      stop
    end if

  end do iterate

  call create_vector ( problem, stress, vec=3 )

! derive vectors

  oldvectors%s(1)%p => sol

! store nodal values of eta in a vector defined per element

  coefficients%i(255)=1

  call loop_over_elements ( mesh, problem, elemsub=fill_eta_gn, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=2

  call derive_vector ( mesh, problem, stress, &
    elemsub=generalized_stokes_stress, &
    coefficients=coefficients, oldvectors=oldvectors )

  print *, 'max stress', maxval(stress%u)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, soln )
  call delete ( stress )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( vel, pres )
  call delete ( elvector )

end program generalized_stokes_channel
