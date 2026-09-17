! Deformation fields for a steady Stokes flow.
! Problem 1: Flow based on the stokes2 problem
!   Stokes problem on a unit square. Periodical boundary conditions.
!   Poiseuille flow.
!   Collocation of nodes.
! Problem 2: Deformation fields computed as a function of time, assuming the
!   flow starts instantaneously at t=0 and was at rest for t<0.
! Problem 3: Stress tensor projection on the interpolation space for F.
! Implicit DG with P1 or P2 interpolation for tau-convection.
! Extrapolated BDF2 time discretization.
! Compute stress tensor for an integral model (postprocessing only).
! Problem based on df6 with the extension:
! Compute stress tensor for an integral model by L2-projection on the finite
! element approximation space of the deformation gradient tensor
! (postprocessing only).


program df8

  use tfem_m
  use hsl_ma57_m
  use hsl_ma41_m
  use stokes_elements_m
  use dfm_elements_m
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
    nx=20,              & ! number of elements in x
    ny=20                 ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp,  & ! viscosity
    flowrate = 2._dp

  integer, parameter :: &
    type_of_model = 0 , & ! 0: separable Rivlin-Sawyers,
                          ! 1: non-separable RS (not yet available)
    spectrum = 0,       & ! 0: discrete spectrum
                          ! 1: Mittag-Leffler (not yet available)
    nummodes = 1,       & ! number of modes in the discrete spectrum
    dampfunc = 1          ! damping function: 1: Lodge, 2: McKinley, 3: PSM

  real(dp), parameter :: &
    G = 1._dp,            & ! modulus
    lambda = 1._dp,       & ! relaxation time
    a = 0.3_dp,           & ! parameter in McKinley damping function
    alpha_PSM = 2.5e4_dp, & ! alpha parameter in PSM damping function
    beta_PSM = 0.25_dp      ! beta parameter in PSM damping function

  integer, parameter :: &
    fintpl = 4,         & ! Q1 interpolation in space for F
    fintpltau = 1,      & ! P1 interpolation in age (tau)
!    fintpltau = 2,      & ! P2 interpolation in age (tau)
    nintvaltau = 20,    & ! number of age intervals
    ncompf = 4,         & ! number of components of deformation tensor F
    ncompt = 4,         & ! number of stress components
    ninttau = fintpltau+1, & ! number of fields within an interval
    numtimesteps = 100    ! number of time steps

  real(dp), parameter :: &
    deltat = 2.e-1_dp,   & ! time step
    beta = 1,            & ! SUPG factor
    dtau1 = 1.5e-1_dp,   & ! length of first interval in age (tau) direction
    tauc = 15._dp,       & ! maximum age (cutoff)
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
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

  type(input_probdef_t) :: input_probdef_f
  type(problem_t), target :: problem_f
  type(sysmatrix_t) :: sysmatrix_f
  type(subscript_t) :: fval
  type(oldvectors_t) :: oldvectors_f

  type(sysvector_t), dimension(ncompf,nintvaltau), target :: &
    sol_f, soln_f, solnm1_f
  type(sysvector_t), dimension(ncompf) :: rhs_f

  type(vector_t) :: ffield, stress_tensor

  type(lu_ma41_t) :: lu_f
  type(solver_options_ma41_t) :: solver_options_f

  type(input_probdef_t) :: input_probdef_s_proj
  type(problem_t), target :: problem_s_proj
  type(sysmatrix_t) :: sysmatrix_s_proj
  type(sysvector_t), dimension(ncompt), target :: sol_s_proj
  type(sysvector_t), dimension(ncompt) :: rhs_s_proj
  type(lu_ma57_t) :: lu_s_proj

  integer :: icomp, step, i, j
  real(dp) :: x(0:nintvaltau)
  real(dp) :: tau_nodal_points(ninttau,nintvaltau)


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=400, ncoefr=350, &
    ncoefra=[nintvaltau+1,nummodes,nummodes] )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]

  coefficients%i(351:353) = [ fintpl, fintpltau, nintvaltau ]
  coefficients%i(354) = 1  ! first-order time integration
  coefficients%i(355) = 1  ! SUPG
  coefficients%i(357) = 1  ! compute velocity gradient directly

  coefficients%i(363) = type_of_model
  coefficients%i(364) = spectrum
  coefficients%i(365) = nummodes
  coefficients%i(366) = dampfunc


  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(6) = flowrate

  coefficients%r(301) = deltat
  coefficients%r(302) = beta

  coefficients%r(303) = a
  coefficients%r(304) = alpha_PSM
  coefficients%r(305) = beta_PSM

  call distribute_elements ( nintvaltau, x, ratio=7, factor=dtau1/tauc )

  coefficients%ra(1)%a = tauc * x(0:nintvaltau)
  coefficients%ra(2)%a = [ G ]
  coefficients%ra(3)%a = [ lambda ]

  write(*,'(/a,g0.4/)') ' non-equidistant grid; last interval = ', &
     coefficients%ra(1)%a(nintvaltau+1)-coefficients%ra(1)%a(nintvaltau)

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

  call define_essential ( mesh, input_probdef, curve1=1, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=3, physq=1 )
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
    curve1=1, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=1, value=0._dp )
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

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

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

! write to a fig file for plotting

  plot_options%scalevector=0.05
  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  plot_options%numlevels=9
  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )


! problem definition deformation tensor

  call create_input_probdef ( mesh, input_probdef_f, nvec=3, nphysq=1 )

  input_probdef_f%vec_elementdof(1)%a(:,1) = [1,0,1,0,1,0,1,0,0] * ninttau ! f
  input_probdef_f%vec_elementdof(1)%a(:,2) = 1 ! scalar for plotting
  input_probdef_f%vec_elementdof(1)%a(:,3) = 4 ! tensor for plotting

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


! create system matrix for df fields problem

  call create_sysmatrix_structure_base ( sysmatrix_f, mesh, problem_f )
  call create_sysmatrix_structure_constraint ( sysmatrix_f, mesh, problem_f )
  call finalize_sysmatrix_structure ( sysmatrix_f )

  call create_sysmatrix_data ( sysmatrix_f )


! create the structure oldvectors_f

  call create_oldvectors ( oldvectors_f, nsysvec=1, nsysvec1=1, nsysvec2=3, &
    nprob=3 )

! store solution vectors and problem structures

  oldvectors_f%s(1)%p => sol
  oldvectors_f%s2(1)%p => sol_f
  oldvectors_f%s2(2)%p => soln_f
  oldvectors_f%s2(3)%p => solnm1_f
  oldvectors_f%p(1)%p => problem
  oldvectors_f%p(2)%p => problem_f


! problem definition for L2-projection of the stress tensor

  call create_input_probdef ( mesh, input_probdef_s_proj, nvec=1, nphysq=1 )

  input_probdef_s_proj%vec_elementdof(1)%a = &
      reshape ( [ 1,0,1,0,1,0,1,0,0 ], &
                   [9,1] )

  input_probdef_s_proj%physq = [1]
  input_probdef_s_proj%probnr = 3

  call problem_definition ( input_probdef_s_proj, mesh, problem_s_proj )

  call create ( problem_s_proj, sol_s_proj, rhs_s_proj )

! store solution vectors and problem structures

  oldvectors_f%s1(1)%p => sol_s_proj
  oldvectors_f%p(3)%p => problem_s_proj

! create and build system matrix for projection problem
! NOTE matrix remains constant and needs to be build once.

  call create_sysmatrix_structure ( sysmatrix_s_proj, mesh, problem_s_proj, &
    symmetric=.true. )
  call create_sysmatrix_data ( sysmatrix_s_proj )

  call build_system ( mesh, problem_s_proj, sysmatrix_s_proj, &
    msysvector=rhs_s_proj, elemsub=stress_projection_elem, &
    oldvectors=oldvectors_f, coefficients=coefficients, &
    buildvector=.false. )

  call check ( sysmatrix_s_proj )


! open monitor file

  open ( unit=13, file='out', recl=300 )


! time stepping

  call tic

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

!     periodical condition on df fields tensor
      call build_system_constraint ( mesh, problem_f, sysmatrix_f, &
        msysvector=rhs_f, elemsub=stokes_constr_node_conn, &
        addmatvec=.true. )

      call check ( sysmatrix_f )

!     solve df fields and keep LU decomposition in the loop over components

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

    write(unit=13,fmt=*) &
      step, step*deltat, &
      maxval(sol_f(1,j)%u(fval%s)), minval(sol_f(1,j)%u(fval%s)), &
      maxval(sol_f(2,j)%u(fval%s)), minval(sol_f(2,j)%u(fval%s)), &
      maxval(sol_f(3,j)%u(fval%s)), minval(sol_f(3,j)%u(fval%s)), &
      maxval(sol_f(4,j)%u(fval%s)), minval(sol_f(4,j)%u(fval%s))

  end do

  call toc ( 'all steps' )


! plot a field

  call create_vector ( problem_f, ffield, vec=2 )

  coefficients%i(358)=2 ! Fxy
  coefficients%i(359)=nintvaltau ! interval
  coefficients%i(360)=ninttau ! field

  call derive_vector ( mesh, problem_f, ffield, &
    elemsub=dfm_deriv_field, &
    coefficients=coefficients, oldvectors=oldvectors_f )

  call plot_color_contour ( plot_options, mesh, problem_f, &
    'field.fig', vector=ffield )


! output values of a field as a function of tau for a certain point in space

  open ( unit=14, file='out2', recl=300 )

  call dfm_tau_nodal_points ( coefficients, tau_nodal_points )

  call create_subscript ( mesh, problem_f, fval, physqarr=[1], points=[1] )

  do j = 1, nintvaltau
    do i = 1, ninttau
      write(unit=14,fmt=*) tau_nodal_points(i,j), &
        sol_f(1,j)%u(fval%s(i)), sol_f(2,j)%u(fval%s(i)), &
        sol_f(3,j)%u(fval%s(i)), sol_f(4,j)%u(fval%s(i))
    end do
    write(unit=14,fmt=*)
  end do


! post processing of the stress tensor

  call create_vector ( problem_f, stress_tensor, vec=3 )

  call derive_vector ( mesh, problem_f, stress_tensor, &
    elemsub=dfm_deriv_stress_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_f )

  call plot_color_contour ( plot_options, mesh, problem_f, &
    'stress_xx.fig', vector=stress_tensor, degfd=1 )
  call plot_color_contour ( plot_options, mesh, problem_f, &
    'stress_xy.fig', vector=stress_tensor, degfd=2 )
  call plot_color_contour ( plot_options, mesh, problem_f, &
    'stress_yy.fig', vector=stress_tensor, degfd=3 )
  if ( dampfunc /= 1 ) then
    call plot_color_contour ( plot_options, mesh, problem_f, &
      'stress_zz.fig', vector=stress_tensor, degfd=4 )
  end if


! post processing of the projected stress tensor

  call solve_stress_projection

  call derive_vector ( mesh, problem_f, stress_tensor, &
    elemsub=dfm_deriv_stress_tensor_proj, &
    coefficients=coefficients, oldvectors=oldvectors_f )

  call plot_color_contour ( plot_options, mesh, problem_f, &
    'stress_xx_proj.fig', vector=stress_tensor, degfd=1 )
  call plot_color_contour ( plot_options, mesh, problem_f, &
    'stress_xy_proj.fig', vector=stress_tensor, degfd=2 )
  call plot_color_contour ( plot_options, mesh, problem_f, &
    'stress_yy_proj.fig', vector=stress_tensor, degfd=3 )
  if ( dampfunc /= 1 ) then
    call plot_color_contour ( plot_options, mesh, problem_f, &
      'stress_zz_proj.fig', vector=stress_tensor, degfd=4 )
  end if


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
  call delete ( ffield, stress_tensor )
  call delete ( fval )

  call delete ( sysmatrix_s_proj )
  call delete ( problem_s_proj )
  call delete ( input_probdef_s_proj )
  call delete ( sol_s_proj, rhs_s_proj )
  call delete ( lu_s_proj )

contains

  subroutine solve_stress_projection

    type(solver_options_ma57_t) :: solver_options_ma57

    integer :: i

!   build vector only (matrix is constant)

    call build_system ( mesh, problem_s_proj, sysmatrix_s_proj, &
      msysvector=rhs_s_proj, elemsub=stress_projection_elem, &
      oldvectors=oldvectors_f, coefficients=coefficients, &
      buildmatrix=.false. )

    ! MA57 solver storage
    solver_options_ma57%integer_storage = 1.3
    solver_options_ma57%real_storage    = 1.3

!   LU decomposition is done in the first call only

    do i = 1, ncompt

      call add_effect_of_essential_to_rhs ( problem_s_proj, sysmatrix_s_proj, &
         sol_s_proj(i), rhs_s_proj(i) )

      call solve_system_ma57 ( sysmatrix_s_proj, rhs_s_proj(i), &
         sol_s_proj(i), lu_s_proj, solver_options=solver_options_ma57 )

    end do

  end subroutine solve_stress_projection

end program df8
