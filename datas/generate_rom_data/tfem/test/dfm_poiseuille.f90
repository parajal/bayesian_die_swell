! Deformation fields for a steady Stokes flow.
! Problem 1: Flow based on the stokes2 problem
!   Stokes problem on a unit square. Periodical boundary conditions.
!   Poiseuille flow.
!   Collocation of nodes.
! Problem 2: Deformation fields computed as a function of time, assuming the
!   flow starts instantaneously at t=0 and was at rest for t<0.
! Implicit DG with P1 or P2 interpolation for tau-convection.
! Extrapolated BDF2 time discretization.

program df3

  use tfem_m
  use hsl_ma57_m
  use hsl_ma41_m
  use stokes_elements_m
  use dfm_elements_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=2,               & ! number of elements in x
    ny=2                  ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp,  & ! viscosity
    flowrate = 2._dp

  integer, parameter :: &
    fintpl = 4,         & ! Q1 interpolation in space for F
    fintpltau = 1,      & ! P1 interpolation in age (tau)
!    fintpltau = 2,      & ! P2 interpolation in age (tau)
    nintvaltau = 20,    & ! number of age intervals
    ncompf = 4,         & ! number of components of deformation tensor F
    ninttau = fintpltau+1, & ! number of fields within an interval
    numtimesteps = 3    ! number of time steps

  real(dp), parameter :: &
    deltat = 1.e-2_dp,   & ! time step
    beta = 1,            & ! SUPG factor
    dtau1 = 1e-1_dp,     & ! length of first interval in age (tau) direction
    tauc = 2._dp,        & ! maximum age (cutoff)
    rs_up = 1.1_dp,      & ! real_storage for the velocity-pressure LU (HSL)
    is_up = 1.0_dp,      & ! integer_storage for the velocity-pressure LU (HSL)
    rs_f  = 1.0_dp,      & ! real_storage for the deformation LU (HSL)
    is_f  = 1.4_dp         ! integer_storage for the deformation LU (HSL)

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma57_t) :: solver_options

  type(input_probdef_t) :: input_probdef_f
  type(problem_t), target :: problem_f
  type(sysmatrix_t) :: sysmatrix_f
  type(subscript_t) :: fval
  type(oldvectors_t) :: oldvectors_f

  type(sysvector_t), dimension(ncompf,nintvaltau), target :: &
    sol_f, soln_f, solnm1_f
  type(sysvector_t), dimension(ncompf) :: rhs_f

  type(lu_ma41_t) :: lu_f
  type(solver_options_ma41_t) :: solver_options_f

  integer :: icomp, step, j
  real(dp) :: x(0:nintvaltau)


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=400, ncoefr=350, &
    ncoefra=[nintvaltau+1] )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gauss ]

  coefficients%i(351:353) = [ fintpl, fintpltau, nintvaltau ]
  coefficients%i(354) = 1  ! first-order time integration
  coefficients%i(355) = 1  ! SUPG
  coefficients%i(357) = 1  ! compute velocity gradient directly


  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(6) = flowrate

  coefficients%r(301) = deltat
  coefficients%r(302) = beta

  call distribute_elements ( nintvaltau, x, ratio=7, factor=dtau1/tauc )

  coefficients%ra(1)%a = tauc * x(0:nintvaltau)

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
  input_probdef%probnr = 1

  call define_essential ( mesh, input_probdef, curves=[1,3], physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, nglobalc=1 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, curve2=5, discretization='collocation', exclude=3 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curves=[1,3], physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficients  )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%real_storage=rs_up
  solver_options%integer_storage=is_up

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! derive vectors

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )


! problem definition deformation tensor

  call create_input_probdef ( mesh, input_probdef_f, nvec=2, nphysq=1 )

  input_probdef_f%vec_elementdof(1)%a(:,1) = [1,0,1,0,1,0,1,0,0] * ninttau ! f
  input_probdef_f%vec_elementdof(1)%a(:,2) = 1 ! scalar for plotting

  input_probdef_f%physq = [1]
  input_probdef_f%probnr = 2

! constraint for periodical boundary conditions of the deformation tensor
  call define_constraint ( mesh, input_probdef_f, curve1=2, curve2=5, &
    discretization='collocation' )

  call problem_definition ( input_probdef_f, mesh, problem_f )

! create a vector subscript for the deformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problem_f, fval, physqarr=[1] )

! create system vectors (solution and right-hand side) for F and initialize

  call create ( problem_f, sol_f, soln_f, solnm1_f )
  call create ( problem_f, rhs_f )

! initial solution

  do j = 1, nintvaltau
    soln_f(1,j)%u = 1 ! initial F_xx
    soln_f(2,j)%u = 0 ! initial F_xy
    soln_f(3,j)%u = 0 ! initial F_yx
    soln_f(4,j)%u = 1 ! initial F_yy
  end do


! create system matrix for conformation problem

  call create_sysmatrix_structure_base ( sysmatrix_f, mesh, problem_f )
  call create_sysmatrix_structure_constraint ( sysmatrix_f, mesh, problem_f )
  call finalize_sysmatrix_structure ( sysmatrix_f )

  call create_sysmatrix_data ( sysmatrix_f )


! create the structure oldvectors_f

  call create_oldvectors ( oldvectors_f, nsysvec=1, nsysvec2=3, nprob=2 )

! store solution vectors and problem structures

  oldvectors_f%s(1)%p => sol
  oldvectors_f%s2(1)%p => sol_f
  oldvectors_f%s2(2)%p => soln_f
  oldvectors_f%s2(3)%p => solnm1_f
  oldvectors_f%p(1)%p => problem
  oldvectors_f%p(2)%p => problem_f

! open monitor file


! time stepping

  do step = 1, numtimesteps

    if ( step == 1 ) then
      coefficients%i(354) = 1 ! Euler first time step
    else
      coefficients%i(354) = 2 ! second-order BDF2 extrapolated
    end if

    do j = 1, nintvaltau

!     build (assemble) matrix and vector for the deformation problem

      coefficients%i(362) = j  ! interval number

      call build_system ( mesh, problem_f, sysmatrix_f, msysvector=rhs_f, &
        elemsub=dfm_elem2, oldvectors=oldvectors_f, &
        coefficients=coefficients )

!     periodical condition on conformation tensor
      call build_system_constraint ( mesh, problem_f, sysmatrix_f, &
        msysvector=rhs_f, elemsub=stokes_constr_node_conn, &
        addmatvec=.true. )

      call check ( sysmatrix_f )

!     solve conformation and keep LU decomposition in the loop over components

      solver_options_f%real_storage=rs_f
      solver_options_f%integer_storage=is_f

      do icomp = 1, ncompf
        call solve_system_ma41 ( sysmatrix_f, rhs_f(icomp), &
          sol_f(icomp,j), lu_f, solver_options=solver_options_f  )
      end do

      call delete ( lu_f )  ! remove LU decomposition and rebuild next interval

    end do

    call copy ( soln_f, solnm1_f )
    call copy ( sol_f, soln_f )

!   write monitor data

    j = nintvaltau

    write(*,*) &
      step, step*deltat, &
      maxval(sol_f(1,j)%u(fval%s)), minval(sol_f(1,j)%u(fval%s)), &
      maxval(sol_f(2,j)%u(fval%s)), minval(sol_f(2,j)%u(fval%s)), &
      maxval(sol_f(3,j)%u(fval%s)), minval(sol_f(3,j)%u(fval%s)), &
      maxval(sol_f(4,j)%u(fval%s)), minval(sol_f(4,j)%u(fval%s))

  end do



! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

  call delete ( problem_f )
  call delete ( input_probdef_f )
  call delete ( sol_f, soln_f, solnm1_f )
  call delete ( rhs_f )
  call delete ( sysmatrix_f )
  call delete ( oldvectors_f )
  call delete ( fval )

end program df3
