! Startup of the flow in a cylinder of a viscoelastic fluid with temperature
! dependent modulus and relaxation time.
! Viscoelastic model: 2-mode UCM/Oldroyd-B fluid.
! Flow computed using a flow rate constraint.
! Temperature distribution computed using the unsteady energy balance equation
! including the stress work as a source term:
!
!   rho cp ( dT/dt + u.nabla T ) - kappa nabla^2 T = tau:D
!
! The initial temperature of the fluid and the temperature of the fluid
! entering the domain is t0.
! Zero-flux (isolated) boundary conditions are assumed on the other boundaries.
! The temperature of the fluid will increase due to the stress work term.

program energy3

  use tfem_m
  use hsl_ma41_m
  use hsl_ma57_m
  use generalized_stokes_elements_m
  use convection_diffusion_supg_elements_m
  use io_utils_m
  use subs_m
  use viscoelastic_elements_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,         & ! velocity interpolation
    pintpl = 2,         & ! pressure interpolation
    tintpl = 2,         & ! temperature interpolation
    gintpl = 2,         & ! gradients interpolation.
    cintpl = 2,         & ! conformation interpolation
    physqvel = 1,       & ! physical quantity nr of the velocity
    physqpress = 2,     & ! physical quantity nr of the pressure
    gauss = 6,          & ! 6 point integration of triangles
    gaussb = 4,         & ! 4 point integration of lines
    ncompc = 4,         & ! number of conformation tensor components
    nmodes = 2,         & ! number of modes
    ncompg = 4,         & ! number of velocity gradient tensor components
    startm = 501,       & ! start of material model data
    model = 2,          & ! UCM/Oldroyd-B model
    poselvec = 4,       & ! elvector position stress work (=3+max(1,numnlpar))
    coorsys = 1,        & ! axisymmetric coordinate system
    nz = 20,            & ! number of elements in z-direction
    nr = 10,            & ! number of elements in r-direction
    ibeta = 2,          & ! method to compute beta in SUPG
    numtimesteps = 200    ! number of time steps

! variables

  integer :: &
    timeint_t1 = 1,     & ! first-order, semi-implicit Euler
    timeint_t2 = 2,     & ! second-order, semi-implicit BDF2
    timeint_ve1 = 1,    & ! first-order, semi-implicit Euler
    timeint_ve2 = 7,    & ! second-order, semi-implicit BDF2 with conformation
                          ! prediction
    logc = 1              ! standard scheme or log transformation

  real(dp), dimension(nmodes) ::  &
    G0 = [10._dp, 2._dp],     & ! maximum polymer modulus for each mode
    lambda0 = [1._dp, 2._dp]    ! maximum polymer relaxation time for each mode

  real(dp) ::  &
    eta_s = 0.1_dp,      & ! (constant) solvent viscosity
    kappa = 0.1_dp,      & ! thermal conductivity coefficient
    rhocp = 1._dp,       & ! density times heat capacity (at constant pressure)
    flowrate = 1.6_dp,   & ! flowrate of the inflow
    t0 = 0._dp,          & ! initial and inflow temperature
    deltat = 0.02_dp,    & ! time step
    oz = 0._dp,          & ! z-coordinate of origin of the mesh
    lz = 2._dp,          & ! length of the domain in z-direction
    R = 1._dp,           & ! radius of the domain
    rs_proj = 1.3_dp,    & ! real_storage for the projection LU (HSL)
    is_proj = 1.3_dp,    & ! integer_storage for projection LU (HSL)
    rs_t = 1.6_dp,       & ! real_storage for the temperature LU (HSL)
    is_t = 1.6_dp,       & ! integer_storage for temperature LU (HSL)
    rs_c = 1.4_dp,       & ! real_storage for the conformation LU (HSL)
    is_c = 1.6_dp,       & ! integer_storage for conformation LU (HSL)
    rs_up = 1.6_dp,      & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.6_dp         ! integer_storage velocity-pressure LU (HSL)

! variables

  integer :: i, step = 0, grp, elem

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdeft
  type(input_probdef_t) :: input_probdefc, input_probdefc_projc
  type(problem_t), target :: problem, problemt
  type(problem_t), target :: problemc, problemc_projc
  type(sysmatrix_t) :: sysmatrix, sysmatrixt, sysmatrixc, sysmatrixc_projc
  type(sysvector_t), target :: sol, sol_n, sol_nm1, solhat, rhsd
  type(sysvector_t), target :: solt, solt_n, solt_nm1, rhsdt, solthat
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc_np1, solc_n, &
    solc_nm1, solchat, rhsc, solc_projcn, solc_projcnm1, rhsc_projc, &
    solchat_proj
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors_ve, oldvectors_en
  type(oldvectors_t) :: oldvectors_sw1, oldvectors_sw2
  type(solver_options_ma41_t) :: solver_options
  type(vector_t) :: modulus_post, lambda_post, temperature_post
  type(vector_t) :: czz_post, crr_post, czr_post, ctt_post
  type(vector_t) :: tauviscous, tauviscoelastic
  type(elvector_t), target :: energy_source1, energy_source2
  type(elvector_t), target :: modulus, lambda, modulus_nodes
  type(lu_ma41_t) :: lu_c

! definitions used in the projection of the gradients

  type(input_probdef_t) :: input_probdef_grad
  type(problem_t) :: problem_grad
  type(sysmatrix_t) :: sysmatrix_grad
  type(sysvector_t) :: sol_grad, rhsd_grad(ncompg)
  type(vector_t), target :: gradients

  character(len=199) :: filename, dataname

  logical :: sw_cproj = .true. ! use c projection in stress work
  integer :: m, icomp
  integer, dimension(3) :: vertices=[1,3,5]

! pass variables to module

  leta0 = eta_s
  lt0   = t0
! the following statement not needed for >=f2003: allocation on assignment
!  allocate ( lG0(nmodes), llam0(nmodes) )
  lG0 = G0
  llam0 = lambda0

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=450, ncoefr=504 )

  coefficients%i = 0

  coefficients%i(1:23) = &
    [ uintpl,   pintpl,      0,     0,      gintpl, &
      physqvel, physqpress,  0,     0,      gauss,  &
      gaussb,   cintpl,      0,     0,      0,      &
      0,        0,           model, nmodes, startm, &
      logc,     timeint_ve1, coorsys ]

  coefficients%i(23) = coorsys
  coefficients%i(49) = 1  ! exp(s) projection =.true. for logc=1
  coefficients%i(60) = 1  ! Separate vector for velocity gradient (No DEVSS).

  coefficients%i(401) = timeint_t1
  coefficients%i(402) = tintpl
  coefficients%i(403) = 1 ! 0: Galerkin, 1: SUPG
  coefficients%i(404) = 1 ! include diffusion term
  coefficients%i(405) = 0 ! no self-dependent source term
  coefficients%i(406) = 4 ! f given in integration points
  coefficients%i(411) = ibeta  ! method to compute beta in SUPG

  coefficients%r = 0
  coefficients%r(1) = eta_s
  coefficients%r(6) = flowrate
  coefficients%r(8) = deltat
  coefficients%r(351:353) = [ deltat, rhocp, kappa ]

! NOTE: the modulus supplied here is only used if coefficients%i(68)==0 and
! the relaxation time supplied here is only used if coefficients%i(69)==0
  coefficients%r(501:504) = [G0(1), lambda0(1), G0(2), lambda0(2)]

  coefficients%i(68) = 3 ! modulus given in integration points
  coefficients%i(69) = 3 ! relaxation time given in integration points
  coefficients%i(73) = 1 ! use varying viscoelastic coefficients
  if ( sw_cproj ) then
    coefficients%i(74) = 1 ! 1: use the least square projection of c=exp(s)
                           ! in the viscoelastic stress work term
  else
    coefficients%i(74) = 0 ! use c=exp(s) directly in the viscoelastic
                           ! stress work term
  end if

! create mesh

  meshgen_options%elshape = 4 ! 6-node triangles
  meshgen_options%nx = nz
  meshgen_options%ny = nr
  meshgen_options%ox = oz
  meshgen_options%lx = lz
  meshgen_options%ly = R

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] ) ! curve 5

  call fill_mesh_parts ( mesh )

! problem definition for the flow problem

  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
     reshape ( [ 2,2,2,2,2,2,   & ! velocity
                 1,0,1,0,1,0,   & ! pressure
                 1,1,1,1,1,1,   & ! scalar, such as vorticity
                 4,0,4,0,4,0,   & ! gradients
                 4,4,4,4,4,4 ], & ! symmetric tensor
                   [6,5] )

  input_probdef%physq = [1,2]
  input_probdef%probnr = 1

  call define_essential ( mesh, input_probdef, curve1=1, curve2=2, &
    physq=physqvel, degfd=[0,1] )
  call define_essential ( mesh, input_probdef, curve1=3, physq=physqvel, &
    degfd=[1,1] )
  call define_essential ( mesh, input_probdef, curve1=4, physq=physqvel, &
    degfd=[0,1] )

  call define_constraint ( mesh, input_probdef, physq=physqvel, curve1=5, &
    nglobalc=1 )

  call problem_definition ( input_probdef, mesh, problem )

! problem definition for the temperature problem

  call create_input_probdef ( mesh, input_probdeft, nvec=2, nphysq=1 )

  input_probdeft%vec_elementdof(1)%a =   &
     reshape ( [ 1,0,1,0,1,0,   &  ! temperature
                 1,1,1,1,1,1 ], &  ! scalar in all nodes
                   [6,2] )

  input_probdeft%physq = [1]
  input_probdef%probnr = 2

  call define_essential ( mesh, input_probdeft, curve1=4 )

  call problem_definition ( input_probdeft, mesh, problemt )

! define viscoelastic problem, essential conditions and constraints

  call create_input_probdef ( mesh, input_probdefc, nvec=3, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a =   &
      reshape ( [ 1,0,1,0,1,0,   &  ! c
                  1,1,1,1,1,1,   &  ! scalar for plotting
                  4,4,4,4,4,4 ], &  ! symmetric tensor
                  [6,3] )

  call define_essential ( mesh, input_probdefc, curve1=4 )

  input_probdefc%physq = [1]
  input_probdefc%probnr = 3

  call problem_definition ( input_probdefc, mesh, problemc )

! problem definition for projected "c=exp(s)" of the log conformation s

  call create_input_probdef ( mesh, input_probdefc_projc, nvec=1, nphysq=1 )

  input_probdefc_projc%vec_elementdof(1)%a = &
          reshape ( [ 1,0,1,0,1,0 ], &
                      [6,1] )

  input_probdefc_projc%physq = [1]
  input_probdefc_projc%probnr = 4

  call problem_definition ( input_probdefc_projc, mesh, problemc_projc )

! create system matrix

  call create_sysmatrix_structure ( sysmatrixc_projc, mesh, problemc_projc, &
    symmetric=.true. )
  call create_sysmatrix_data ( sysmatrixc_projc )

  call create ( problemc_projc, solc_projcn, solc_projcnm1, rhsc_projc )

! problem definition for gradients

  call create_input_probdef ( mesh, input_probdef_grad, nvec=2, nphysq=1 )

  input_probdef_grad%vec_elementdof(1)%a(:,1) = 0
  input_probdef_grad%vec_elementdof(1)%a(vertices,1) = 1  ! gradient component

  input_probdef_grad%physq = [1]
  input_probdef_grad%probnr = 5

  call problem_definition ( input_probdef_grad, mesh, problem_grad )

! create system vectors (solution and right-hand side)

  call create ( problem_grad, sol_grad )
  call create ( problem_grad, rhsd_grad )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix_grad, mesh, problem_grad, &
    symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix_grad )

! create system vectors (solution and right-hand side)

  call create ( problem, sol, sol_n, sol_nm1, solhat, rhsd )
  call create ( problemt, solt, solt_n, solt_nm1, solthat, rhsdt )
  call create ( problemc, solc_np1, solc_n, solc_nm1, rhsc )
  if ( sw_cproj ) then
    call create ( problemc_projc, solchat_proj )
  else
    call create ( problemc, solchat )
  end if

! create gradient vector

  call create ( problem, gradients, vec=4 )

! create vector defined per element for eta

  call create ( mesh, modulus, nreal1d=nmodes )
  call create_elvector ( mesh, lambda, nreal1d=nmodes )

! create vector defined per element for scalar values of right-hand

  call create ( mesh, energy_source1, nreal1d=1 )
  call create ( mesh, energy_source2, nreal1d=1 )

! oldvectors for the flow problem

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec2=3, nprob=3, &
    nvec=3, nelvec=2 )
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemt
  oldvectors_ve%p(3)%p => problemc_projc
  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s(2)%p => solt
  oldvectors_ve%s2(1)%p => solc_n
  oldvectors_ve%s2(2)%p => solc_nm1
  oldvectors_ve%s2(3)%p => solc_projcn
  oldvectors_ve%v(3)%p => gradients
  oldvectors_ve%e(1)%p => modulus
  oldvectors_ve%e(2)%p => lambda

! oldvectors for the viscous dissipation computation

  call create_oldvectors ( oldvectors_sw1, nsysvec=1, nelvec=2 )
  oldvectors_sw1%s(1)%p => solhat
  oldvectors_sw1%e(2)%p => energy_source1

! oldvectors for the viscoelastic stress work computation

  call create_oldvectors ( oldvectors_sw2, nprob=3, nsysvec=2, nsysvec2=4, &
    nelvec=poselvec )
  oldvectors_sw2%p(1)%p => problem
  oldvectors_sw2%p(2)%p => problemt
  oldvectors_sw2%p(3)%p => problemc_projc
  oldvectors_sw2%s(1)%p => solhat
  oldvectors_sw2%s(2)%p => solthat
  if ( sw_cproj ) then
    oldvectors_sw2%s2(2)%p => solchat_proj
  else
    oldvectors_sw2%s2(1)%p => solchat
  end if
  oldvectors_sw2%e(1)%p => modulus
  oldvectors_sw2%e(2)%p => lambda
  oldvectors_sw2%e(poselvec)%p => energy_source2

! oldvectors for the temperature problem

  call create_oldvectors ( oldvectors_en, nsysvec=4, nprob=2, nelvec=2 )
  oldvectors_en%p(1)%p => problem
  oldvectors_en%s(1)%p => solhat
  oldvectors_en%s(2)%p => solt_n
  oldvectors_en%s(3)%p => solt_nm1
  oldvectors_en%e(2)%p => energy_source1


! initialize solution vectors

  sol%u = 0
  solt%u = 0
  do m = 1, nmodes
    solc_projcn(1,m)%u = 1   ! initial czz
    solc_projcn(2,m)%u = 0   ! initial czr
    solc_projcn(3,m)%u = 1   ! initial crr
    solc_projcn(4,m)%u = 1   ! initial ctt
  end do

! fill solution vector with essential boundary conditions for velocity field

  call fill_sysvector ( mesh, problem, sol, curve1=1, curve2=2, &
    physq=physqvel, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, curve1=3, physq=physqvel, &
    value=0._dp )
  call fill_sysvector ( mesh, problem, sol, curve1=4, physq=physqvel, &
    degfd=2, value=0._dp )

! fill initial temperature field

  call fill_sysvector ( mesh, problemt, solt, node1=1, node2=mesh%nnodes, &
    value=t0 )

! initialize vectors with zero stress

  do m = 1, nmodes
    solc_np1(1,m)%u = 0   ! initial szz
    solc_np1(2,m)%u = 0   ! initial szr
    solc_np1(3,m)%u = 0   ! initial srr
    solc_np1(4,m)%u = 0   ! initial ctt
  end do

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! create system matrix

  call create_sysmatrix_structure ( sysmatrixt, mesh, problemt )
  call create_sysmatrix_data ( sysmatrixt )

! create system matrix for the ve problem

  call create_sysmatrix_structure ( sysmatrixc, mesh, problemc )
  call create_sysmatrix_data ( sysmatrixc )

! solve initial velocity field (w/o viscoelastic stress)

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, oldvectors=oldvectors_ve, coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%real_storage=rs_up
  solver_options%integer_storage=is_up

  call solve_system_ma41 ( sysmatrix, rhsd, sol, solver_options=solver_options )

! postprocessing

  call create ( problem, modulus_post, vec=3 )
  call create ( problem, lambda_post, vec=3 )
  call create ( problemt, temperature_post, vec=2 )
  call create ( problemc, czz_post, vec=2 )
  call create ( problemc, czr_post, vec=2 )
  call create ( problemc, crr_post, vec=2 )
  call create ( problemc, ctt_post, vec=2 )
  call create ( mesh, modulus_nodes, nreal1d=nmodes )
  call create_vector ( problem, tauviscous, vec=5 )
  call create_vector ( problemc, tauviscoelastic, vec=3 )

  call postprocessing

! start time stepping

  call copy ( sol, sol_n )
  call copy ( solt, solt_n )
  call copy ( solc_n, solc_nm1 )
  call copy ( solc_np1, solc_n )

  do step = 1, numtimesteps

    if ( step == 2 ) then
      coefficients%i(22) = timeint_ve2
      coefficients%i(401) = timeint_t2
    end if

!   project the conformation tensor from c=exps(s)

    call solve_exps_projection

!   predict several fields that are needed in the stress work term

    if ( step == 1 ) then ! first-order

      solhat%u  = sol_n%u
      solthat%u = solt_n%u

      if ( sw_cproj ) then
        do m = 1, nmodes
          do i = 1, ncompc
            solchat_proj(i,m)%u = solc_projcn(i,m)%u
          end do
        end do
      else
        do m = 1, nmodes
          do i = 1, ncompc
            solchat(i,m)%u = solc_n(i,m)%u
          end do
        end do
      end if

    else ! second-order

      solhat%u  = 2*sol_n%u - sol_nm1%u
      solthat%u = 2*solt_n%u - solt_nm1%u

      if ( sw_cproj ) then
        do m = 1, nmodes
          do i = 1, ncompc
            solchat_proj(i,m)%u = 2*solc_projcn(i,m)%u - solc_projcnm1(i,m)%u
          end do
        end do
      else
        do m = 1, nmodes
          do i = 1, ncompc
            solchat(i,m)%u = 2*solc_n(i,m)%u - solc_nm1(i,m)%u
          end do
        end do
      end if

    end if

!   viscous dissipation of eta_s
    call loop_over_elements ( mesh, problem, &
      elemsub=fill_viscous_dissipation_gauss, &
      coefficients=coefficients, oldvectors=oldvectors_sw1 )

!   modulus and lambda in Gauss points for temperature field solthat
    call loop_over_elements ( mesh, problem, elemsub=fill_modlam_gauss, &
      coefficients=coefficients, oldvectors=oldvectors_sw2 )

!   stress work of viscoelastic stress
    call loop_over_elements ( mesh, problemc, &
      elemsub=fill_viscoelastic_stress_work_gauss, &
      coefficients=coefficients, oldvectors=oldvectors_sw2 )

!   add the two stress work source terms
    do grp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(grp)
        energy_source1%g(grp)%r1(elem,1)%a = &
                                 energy_source1%g(grp)%r1(elem,1)%a &
                                         + energy_source2%g(grp)%r1(elem,1)%a
      end do
    end do

!   build diffusion matrix and right-hand side

    call build_system ( mesh, problemt, sysmatrixt, rhsdt, &
      elemsub=conv_diff_supg_elem, coefficients=coefficients, &
      oldvectors=oldvectors_en )

    call add_effect_of_essential_to_rhs ( problemt, sysmatrixt, solt, rhsdt )

    solver_options%real_storage=rs_t
    solver_options%integer_storage=is_t

    call solve_system_ma41 ( sysmatrixt, rhsdt, solt, &
      solver_options=solver_options )

!   build velocity-pressure matrix and vector

    call loop_over_elements ( mesh, problem, elemsub=fill_modlam_gauss, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, oldvectors=oldvectors_ve, coefficients=coefficients )

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, oldvectors=oldvectors_ve, &
      physqrow=[physqvel], physqcol=[physqvel], addmatvec=.true., &
      coefficients=coefficients )

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    solver_options%real_storage=rs_up
    solver_options%integer_storage=is_up

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options )

!   build and solve velocity projection problem

    call build_and_solve_proj_grad

!   build (assemble) matrix and vector for conformation problem

    if ( step == 1 ) then
      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem1, oldvectors=oldvectors_ve, &
        coefficients=coefficients )
    else
      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem_implicit_2nd_order, &
        oldvectors=oldvectors_ve, coefficients=coefficients )
    end if

    call check ( sysmatrixc )

!   solve conformation and keep LU decomposition in the loop over components

    solver_options%real_storage=rs_c
    solver_options%integer_storage=is_c

    do m = 1, nmodes
      do icomp = 1, ncompc
        call solve_system_ma41 ( sysmatrixc, rhsc(icomp,m), solc_np1(icomp,m), &
          lu_c, solver_options=solver_options  )
      end do
    end do

    call delete ( lu_c )

    write(*,'(a,i0,a,es12.4/2(a,es26.18)/)') &
      'step = ', step, ' time = ', step*deltat, &
      ' maxsolczz = ', maxval(solc_np1(1,1)%u), &
      ' maxsolt = ', maxval(solt%u)

!   postprocessing

    call postprocessing

!   copy old values

    call copy ( sol_n, sol_nm1 )
    call copy ( sol, sol_n )
    call copy ( solc_n, solc_nm1 )
    call copy ( solc_np1, solc_n )
    call copy ( solt_n, solt_nm1 )
    call copy ( solt, solt_n )
    call copy ( solc_projcn, solc_projcnm1 )

  end do


! delete all data including all allocated memory

  call delete ( problem, problemt )
  call delete ( input_probdef, input_probdeft )
  call delete ( mesh )
  call delete ( sol, sol_n, sol_nm1, solhat )
  call delete ( solt, solt_n, solt_nm1 )
  call delete ( solc_n, solc_nm1, solc_np1 )
  call delete ( solc_projcn, solc_projcnm1 )
  if ( sw_cproj ) then
    call delete ( solchat_proj )
  else
    call delete ( solchat )
  end if
  call delete ( rhsd, rhsdt )
  call delete ( sysmatrix, sysmatrixt )
  call delete ( coefficients )
  call delete ( oldvectors_ve, oldvectors_en )
  call delete ( oldvectors_sw1, oldvectors_sw2 )
  call delete ( modulus, lambda, modulus_nodes )
  call delete ( modulus_post, lambda_post, temperature_post )
  call delete ( czz_post, czr_post, crr_post, ctt_post )
  call delete ( tauviscous, tauviscoelastic )

contains


! subroutine to write vtks

  subroutine postprocessing

    integer :: m

    type(oldvectors_t) :: oldvectors_post
    type(vector_t) :: temperature

    call create ( oldvectors_post, nsysvec=1, nsysvec2=1, nelvec=1 )
    oldvectors_post%s(1)%p => solt

!   derive temperature in all nodes

    call create_vector ( problemt, temperature, vec=2 )

    call derive_vector ( mesh, problemt, temperature, &
      elemsub=derive_q, coefficients=coefficients, &
      oldvectors=oldvectors_post )

!   write vtks

    write(filename,'(a,i4.4,a)') 'sol', step, '.vtk'
    call write_scalar_vtk ( mesh, problemt, vector=temperature, &
      filename=filename, dataname='temperature' )

    write(filename,'(a,i4.4,a)') 'velo', step, '.vtk'
    call write_vector_vtk ( mesh, problem, filename=filename, &
      dataname='velocity', sysvector=sol, physq=physqvel )

    do m = 1, nmodes

!     compute local modulus

      modulus_post%u = compute_modulus ( temperature%u, m )

!     compute local modulus

      lambda_post%u = compute_relaxation_time ( temperature%u, m )

      write(dataname,'(a,i1)') 'G', m
      call write_scalar_vtk ( mesh, problem, filename=filename, &
        dataname=trim(dataname), vector=modulus_post, append=.true. )

      write(dataname,'(a,i1)') 'lambda', m
      call write_scalar_vtk ( mesh, problem, filename=filename, &
        dataname=trim(dataname), vector=lambda_post, append=.true. )

    end do

!   derive the conformation tensor

    oldvectors_post%s2(1)%p => solc_np1

    do m = 1, nmodes

      coefficients%i(28) = m

      coefficients%i(13)=1
      call derive_vector ( mesh, problemc, czz_post, &
        elemsub=deriv_conformation, coefficients=coefficients, &
        oldvectors=oldvectors_post )
      coefficients%i(13)=2
      call derive_vector ( mesh, problemc, czr_post, &
        elemsub=deriv_conformation, coefficients=coefficients, &
        oldvectors=oldvectors_post )
      coefficients%i(13)=3
      call derive_vector ( mesh, problemc, crr_post, &
        elemsub=deriv_conformation, coefficients=coefficients, &
        oldvectors=oldvectors_post )
      coefficients%i(13)=4
      call derive_vector ( mesh, problemc, ctt_post, &
        elemsub=deriv_conformation, coefficients=coefficients, &
        oldvectors=oldvectors_post )

      write(filename,'(a,i1,a,i4.4,a)') 'c_mode',m,'_', step, '.vtk'

      call write_scalar_vtk ( mesh, problemc, vector=czz_post, &
        filename=filename, dataname='czz' )

      call write_scalar_vtk ( mesh, problemc, vector=czr_post, &
        filename=filename, dataname='czr', append=.true. )

      call write_scalar_vtk ( mesh, problemc, vector=crr_post, &
        filename=filename, dataname='crr', append=.true. )

      call write_scalar_vtk ( mesh, problemc, vector=ctt_post, &
        filename=filename, dataname='ctt', append=.true. )

    end do

!   derive viscous stress tensor

    oldvectors_post%s(1)%p => sol

    call derive_vector ( mesh, problem, tauviscous, &
      elemsub=stokes_stress_tensor, coefficients=coefficients, &
      oldvectors=oldvectors_post )

    write(filename,'(a,i4.4,a)') 'stress', step, '.vtk'
    call write_tensor_vtk ( mesh, problem, filename=filename, &
      dataname='tauviscous', vector=tauviscous )

!   derive viscoelastic stress tensor

    oldvectors_post%s(1)%p => solt
    oldvectors_post%e(1)%p => modulus_nodes

    call loop_over_elements ( mesh, problemt, elemsub=fill_modulus_nodes, &
      coefficients=coefficients, oldvectors=oldvectors_post )

    oldvectors_post%s2(1)%p => solc_np1

    coefficients%i(28) = 0  ! total stress

    call derive_vector ( mesh, problemc, tauviscoelastic, &
      elemsub=deriv_viscoelastic_stress_tensor,&
      coefficients=coefficients, oldvectors=oldvectors_post )

    write(filename,'(a,i4.4,a)') 'stress', step, '.vtk'
    call write_tensor_vtk ( mesh, problemc, filename=filename, &
      dataname='tauviscoelastic', vector=tauviscoelastic, append=.true. )

  end subroutine postprocessing


! project c=exp(s) on discrete fem space

  subroutine solve_exps_projection

    type(solver_options_ma57_t) :: solver_options_projc
    type(oldvectors_t) :: oldvectors_proj
    type(lu_ma57_t) :: lu_projc

    integer :: i, m

    call create_oldvectors ( oldvectors_proj, nprob=2, nsysvec2=1 )
    oldvectors_proj%s2(1)%p => solc_n
    oldvectors_proj%p(2)%p => problemc

!   build system matrix and vector for projection problem

    call build_system ( mesh, problemc_projc, sysmatrixc_projc, &
      m2sysvector=rhsc_projc, elemsub=exps_projection_elem, &
      oldvectors=oldvectors_proj, coefficients=coefficients )

    call check ( sysmatrixc_projc )

!   MA57 solver storage

    solver_options_projc%integer_storage = is_proj
    solver_options_projc%real_storage    = rs_proj

!   LU decomposition is done in the first loop traversing

    do m = 1, nmodes
      do i = 1, ncompc

        call add_effect_of_essential_to_rhs ( problemc_projc, &
          sysmatrixc_projc, solc_projcn(i,m), rhsc_projc(i,m) )

        call solve_system_ma57 ( sysmatrixc_projc, rhsc_projc(i,m), &
           solc_projcn(i,m), lu_projc, solver_options=solver_options_projc )

      end do
    end do

  end subroutine solve_exps_projection


! project the gradients of the velocity

  subroutine build_and_solve_proj_grad

    integer :: i

    type(solver_options_ma57_t) :: solver_options_grad
    type(oldvectors_t) :: oldvectors_grad
    type(lu_ma57_t) :: lu_projgrad

    call create_oldvectors ( oldvectors_grad, nprob=1, nsysvec=1 )
    oldvectors_grad%p(1)%p => problem
    oldvectors_grad%s(1)%p => sol

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem_grad, sysmatrix_grad, &
      msysvector=rhsd_grad, elemsub=gradient_from_velocity_elem, &
      coefficients=coefficients, oldvectors=oldvectors_grad )

!   solve the projection problem

    solver_options_grad%real_storage=1.6_dp
    solver_options_grad%integer_storage=1.6_dp

    do i = 1, ncompg

      call add_effect_of_essential_to_rhs ( problem_grad, sysmatrix_grad, &
        sol_grad, rhsd_grad(i) )

      call solve_system_ma57 ( sysmatrix_grad, rhsd_grad(i), sol_grad, &
        solver_options=solver_options_grad, lu=lu_projgrad )

      call transfer_data ( mesh, problem_grad, problem, sysvector1=sol_grad, &
        vector2=gradients, degfd2=[i] )

    end do

  end subroutine build_and_solve_proj_grad

end program energy3
