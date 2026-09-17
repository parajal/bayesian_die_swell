! Deformation fields for a steady Stokes flow.
! Problem 1: Flow based on the stokes2 problem
!   Stokes problem on a unit square. Periodical boundary conditions.
!   Poiseuille flow.
!   Collocation of nodes.
! Problem 2: Deformation fields computed as a function of time, assuming the
!   flow starts instantaneously at t=0 and was at rest for t<0.
! Decoupled fields, explicit DG with P2 interpolation for tau-convection.
!   (similar to the 2001 paper).
! Extrapolated BDF3 time discretization.

program df2

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
    nx=2,               & ! number of elements in x
    ny=2                  ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp,  & ! viscosity
    flowrate = 2._dp

  integer, parameter :: &
    fintpl = 4,         & ! Q1 interpolation in space for F
    fintpltau = 2,      & ! P1 interpolation in age (tau)
    nintvaltau = 20,    & ! number of age intervals
    ncompf = 4,         & ! number of components of deformation tensor F
    ninttau = fintpltau+1, & ! number of fields within an interval
    numtimesteps = 100    ! number of time steps

  real(dp), parameter :: &
    deltat = 8.e-3_dp,  & ! time step
    beta = 1,            & ! SUPG factor
    dtau1 = 1e-1_dp,     & ! length of first interval in age (tau) direction
    tauc = 2._dp,        & ! maximum age (cutoff)
    rs_f  = 1.0_dp,      & ! real_storage for the deformation LU (HSL)
    is_f  = 1.5_dp         ! integer_storage for the deformation LU (HSL)

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

  type(sysvector_t), dimension(ncompf,ninttau,nintvaltau), target :: &
    sol_f, solm1_f, solm2_f
  type(sysvector_t), dimension(ncompf,ninttau,nintvaltau) :: rhs_f

  type(vector_t) :: ffield

  type(lu_ma41_t) :: lu_f
  type(solver_options_ma41_t) :: solver_options_f

  integer :: icomp, step, i, node, j
  real(dp) :: x(0:nintvaltau)
  real(dp) :: tau_nodal_points(ninttau,nintvaltau)


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

  call create_input_probdef ( mesh, input_probdef_f, nvec=2, nphysq=1 )

  input_probdef_f%vec_elementdof(1)%a =   &
      reshape ( [ 1,0,1,0,1,0,1,0,0,    &  ! f
                   1,1,1,1,1,1,1,1,1 ], &  ! scalar for plotting
                   [9,2] )

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

  call create ( problem_f, sol_f, solm1_f, rhs_f )

! initial solution

  do i = 1, ninttau
    do j = 1, nintvaltau
      sol_f(1,i,j)%u = 1 ! initial F_xx
      sol_f(2,i,j)%u = 0 ! initial F_xy
      sol_f(3,i,j)%u = 0 ! initial F_yx
      sol_f(4,i,j)%u = 1 ! initial F_yy
    end do
  end do

! create system matrix for df fields problem

  call create_sysmatrix_structure_base ( sysmatrix_f, mesh, problem_f )
  call create_sysmatrix_structure_constraint ( sysmatrix_f, mesh, problem_f )
  call finalize_sysmatrix_structure ( sysmatrix_f )

  call create_sysmatrix_data ( sysmatrix_f )


! create the structure oldvectors_f

  call create_oldvectors ( oldvectors_f, nsysvec=1, nsysvec3=3, nprob=2 )

! store solution vectors and problem structures

  oldvectors_f%s(1)%p => sol
  oldvectors_f%s3(1)%p => sol_f
  oldvectors_f%s3(2)%p => solm1_f
  oldvectors_f%s3(3)%p => solm2_f
  oldvectors_f%p(1)%p => problem
  oldvectors_f%p(2)%p => problem_f

! open monitor file

  open ( unit=13, file='out', recl=300 )

! time stepping

  call tic

  do step = 1, numtimesteps

    if ( step == 1 ) then
      coefficients%i(354) = 1 ! Euler first time step
    else if ( step == 2 ) then
      coefficients%i(354) = 2 ! 2nd order second time step
    else
      coefficients%i(354) = 0 ! default
    end if

!   build (assemble) matrix and vector for the deformation problem

    call build_system ( mesh, problem_f, sysmatrix_f, m3sysvector=rhs_f, &
      elemsub=dfm_elem1, oldvectors=oldvectors_f, &
      coefficients=coefficients )

!   periodical condition on df fields tensor
    call build_system_constraint ( mesh, problem_f, sysmatrix_f, &
      m3sysvector=rhs_f, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

    call check ( sysmatrix_f )

    call copy ( solm1_f, solm2_f )
    call copy ( sol_f, solm1_f )

!   solve df fields and keep LU decomposition in the loop over components

    solver_options_f%real_storage=rs_f
    solver_options_f%integer_storage=is_f

    do icomp = 1, ncompf
      do i = 1, ninttau
        do j = 1, nintvaltau
          call solve_system_ma41 ( sysmatrix_f, rhs_f(icomp,i,j), &
            sol_f(icomp,i,j), lu_f, solver_options=solver_options_f  )
        end do
      end do
    end do

    call delete ( lu_f )  ! remove LU decomposition and rebuild next time step

!   write monitor data for last field

    i = ninttau
    j = nintvaltau

    write(unit=13,fmt=*) &
      step, step*deltat, &
      maxval(sol_f(1,i,j)%u(fval%s)), minval(sol_f(1,i,j)%u(fval%s)), &
      maxval(sol_f(2,i,j)%u(fval%s)), minval(sol_f(2,i,j)%u(fval%s)), &
      maxval(sol_f(3,i,j)%u(fval%s)), minval(sol_f(3,i,j)%u(fval%s)), &
      maxval(sol_f(4,i,j)%u(fval%s)), minval(sol_f(4,i,j)%u(fval%s))

  end do

  call toc ( 'all steps' )


! post processing of the fields

  call dfm_interpolate_to_tau_nodes ( coefficients, fval, tau_nodal_points, &
    sol_f )

  call create_vector ( problem_f, ffield, vec=2 )

  coefficients%i(358)=2 ! Fxy
  coefficients%i(359)=nintvaltau ! interval
  coefficients%i(360)=ninttau    ! point
  coefficients%i(361)=1 ! location of in oldvectors of the fields is s3

  call derive_vector ( mesh, problem_f, ffield, elemsub=dfm_deriv_field, &
    coefficients=coefficients, oldvectors=oldvectors_f )

  call plot_color_contour ( plot_options, mesh, problem_f, &
    'field.fig', vector=ffield )


! write fields in tau-direction for a single point in space

  open ( unit=14, file='out2', recl=300 )

  node = mesh%points(1)

  do j = 1, nintvaltau
    do i = 1, ninttau
      write(unit=14,fmt=*) tau_nodal_points(i,j), &
        sol_f(1,i,j)%u(fval%s(node)), sol_f(2,i,j)%u(fval%s(node)), &
        sol_f(3,i,j)%u(fval%s(node)), sol_f(4,i,j)%u(fval%s(node))
    end do
    write(unit=14,fmt=*)
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
  call delete ( sol_f, solm1_f, solm2_f, rhs_f )
  call delete ( sysmatrix_f )
  call delete ( oldvectors_f )
  call delete ( ffield )
  call delete ( fval )

end program df2
