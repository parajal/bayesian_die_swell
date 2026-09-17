! Elastoviscoplastic problem on a unit square.
! Start-up of plane Poiseuille flow having a constant flow rate.
! Upper half of the domain using symmetry boundary conditions at the bottom.
! DEVSS-G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! second-order time integration
! This is similar to poiseulle2, but now using the b-formulation.

program poiseuille3

  use tfem_m
  use viscoelastic_elements_m
  use inertia_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities
    pintpl = 4,     & ! Q1 pressures
    gintpl = 4,     & ! Q1 gradients
    bintpl = 4,     & ! Q1 b-tensor
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    gauss = 3,      & ! 3x3 integration of quads
    ncompb = 4,     & ! number of b-tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    coorsys = 0       ! planar Cartesian coordinate system


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, soln, solnm1, sol_hat
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(subscript_t) :: velx, vel, press
  type(solver_options_ma41_t) :: solver_options_u

  type(input_probdef_t) :: input_probdefb
  type(problem_t), target :: problemb
  type(sysmatrix_t) :: sysmatrixb
  type(oldvectors_t) :: oldvectors_ve
  type(subscript_t) :: bval
  type(subscriptvec_t) :: cxx, cxy, cyy
  type(sysvector_t), dimension(ncompb,nmodes), target :: solb, solbm1
  type(sysvector_t), dimension(ncompb,nmodes) :: rhsb
  type(vector_t) :: ctensor
  type(lu_ma41_t) :: lub
  type(solver_options_ma41_t) :: solver_options_b

  type(input_probdef_t) :: input_probdef_proj
  type(problem_t), target :: problem_proj
  type(sysmatrix_t) :: sysmatrix_proj
  type(sysvector_t), dimension(ncompb-1,nmodes), target :: solc_proj
  type(sysvector_t), dimension(ncompb-1,nmodes) :: rhsc_proj
  type(lu_ma57_t) :: lu_proj




! variables

  logical :: &
    cproj = .true.         ! projection of c=b*b^T for cn in momentum balance

  integer :: &
    nx=20,               & ! number of elements in x
    ny=20,               & ! number of elements in y
    timeint1 = 1,        & ! (first-order) time integration (first step)
    timeint2 = 7,        & ! (second-order) time integration
    numtimesteps = 600,  & ! number of time steps
    inertia = 0,         & ! include instationary inertia term
    model = 2,           & ! 2:Oldroyd-B, 3:Giesekus, 5,6:PTT
    alam_model = 2         ! 0:none, 1:elastic, 2: Saramito1, 3: Saramito2

  real(dp) :: &
    rho = 1.0_dp,      & ! density
    eta_s = 0.1_dp,    & ! solvent viscosity
    deltat = 5.e-3_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    U_avg = 1._dp,     & ! average velocity at the entry
    critval_reinit = 0.0_dp,  & ! criterion for reinitializing b
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_b  = 1.0_dp,    & ! real_storage for the b-tensor LU (HSL)
    is_b  = 1.6_dp       ! integer_storage for b-tensor LU (HSL)

  real(dp), dimension(nmodes) :: &
    eta_p,     & ! polymer viscosity
    lambda,    & ! relaxation time
    G,         & ! modulus = eta_p/lambda
    nlnpar,    & ! nonlinear parameter for models
    tau_y,     & ! yield stress
    K,         & ! power-law coefficient
    n            ! power-law exponent

  integer :: icomp, step, i, j, m, npar
  real(dp) :: alpha, mfac
  real(dp) :: flowrate, H
  real(dp) :: gamma0=1.5_dp, alpha0=2._dp, alpha1=-0.5_dp


! namelist for input of variables; read from standard input

  namelist /comppar/ nx, ny, timeint1, timeint2, numtimesteps, inertia, &
    model, rho, eta_s, eta_p, lambda, nlnpar, alam_model, tau_y, K, n, deltat, &
    rs_gup, is_gup, rs_b, is_b, U_avg

  read ( unit=*, nml=comppar )


! set some parameters

  alpha = sum(eta_p)  ! DEVSS parameter
  G = eta_p / lambda  ! modulus


! fill coefficients

  if ( model == 2 ) then
    select case(alam_model)
    case(0,1)
      npar = 2
    case(2)
      npar = 3
    case(3)
      npar = 5
    case default
      write(*,'(/a,i0/)') 'Error: wrong value alam_model: ', alam_model
      stop
    end select
  else if ( any( model == [3,5,6] ) ) then
    select case(alam_model)
    case(0,1)
      npar = 3
    case(2)
      npar = 4
    case(3)
      npar = 6
    case default
      write(*,'(/a,i0/)') 'Error: wrong value alam_model: ', alam_model
      stop
    end select
  end if

  call create_coefficients ( coefficients, ncoefi=250, ncoefr=500+npar*nmodes )

  coefficients%i = 0
  coefficients%i(1:23) = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gauss,    bintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      0,        timeint1,   coorsys                   &
    ]

  coefficients%i(49) = 1  ! exp(s) projection =.true. for logc=1
  coefficients%i(71) = 1  ! CDT variant for b-formulation
  if ( cproj ) coefficients%i(72) = 1  ! use c projection in momentum balance
  coefficients%i(80) = alam_model ! adapted lambda

  coefficients%r = 0
  coefficients%r(1:10) = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
      0._dp, 0._dp,  deltat,   beta,  0._dp  &
    ]
  coefficients%r(151) = rho

  if ( model == 2 ) then
    select case(alam_model)
    case(0,1)
      coefficients%r(501:500+2*nmodes) = &
         [ ( G(i), lambda(i), i=1,nmodes ) ]
    case(2)
      coefficients%r(501:500+3*nmodes) = &
         [ ( G(i), lambda(i), tau_y(i), i=1,nmodes ) ]
    case(3)
      coefficients%r(501:500+5*nmodes) = &
         [ ( G(i), lambda(i), tau_y(i), K(i), n(i), i=1,nmodes ) ]
    case default
      write(*,'(/a,i0/)') 'Error: wrong value alam_model: ', alam_model
      stop
    end select
  else if ( any( model == [3,5,6] ) ) then
    select case(alam_model)
    case(0,1)
      coefficients%r(501:500+3*nmodes) = &
         [ ( G(i), lambda(i), nlnpar(i), i=1,nmodes ) ]
    case(2)
      coefficients%r(501:500+4*nmodes) = &
         [ ( G(i), lambda(i), nlnpar(i), tau_y(i), i=1,nmodes ) ]
    case(3)
      coefficients%r(501:500+6*nmodes) = &
         [ ( G(i), lambda(i), nlnpar(i), tau_y(i), K(i), n(i), i=1,nmodes ) ]
    case default
      write(*,'(/a,i0/)') 'Error: wrong value alam_model: ', alam_model
      stop
    end select
  end if


! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )

! set flowrate

  H = mesh%coor(mesh%points(4),2) ! height is given by y-coordinate of P4
  flowrate = - H * U_avg

  coefficients%r(6) = flowrate


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=3 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [ 4,0,4,0,4,0,4,0,0,    &  ! G
                  2,2,2,2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,1,0,0,    &  ! pressure
                  1,1,1,1,1,1,1,1,1,    &  ! scalar, such as vorticity
                  3,3,3,3,3,3,3,3,3 ], &  ! tensor
                  [9,5] )

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

! Dirichlet boundary conditions

! center line
  call define_essential ( mesh, input_probdef, &
    curve1=1, physq=physqvel, degfd=[0,1] )

! wall
  call define_essential ( mesh, input_probdef, &
    curve1=3, physq=physqvel )

! pressure in point 1
  call define_essential ( mesh, input_probdef, point=1, physq=3 )

! constraints for periodical boundary conditions

  call define_constraint ( mesh, input_probdef, &
    physq=2, curve1=2, curve2=5, discretization='weak', elementdof=[2,0,2] )

! gradients
  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, curve2=5, discretization='collocation' )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=4, nglobalc=1 )

  call problem_definition ( input_probdef, mesh, problem )

! create a vector subscripts for post processing
! horizontal velocity
  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
! velocities
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )
! pressure
  call create_subscript ( mesh, problem, press, physqarr=[physqpress], &
    points=[1,2] )


! problem definition b-tensor

  call create_input_probdef ( mesh, input_probdefb, nvec=4, nphysq=1 )

  input_probdefb%vec_elementdof(1)%a =   &
      reshape ( [ 1,0,1,0,1,0,1,0,0,    &  ! b
                  1,1,1,1,1,1,1,1,1,    &  ! scalar for plotting
                  3,3,3,3,3,3,3,3,3,    &  ! tensor for plotting
                  4,4,4,4,4,4,4,4,4 ], &  ! tensor for plotting
                  [9,4] )

  input_probdefb%physq = [1]
  input_probdefb%probnr = 2

! constraint for periodical boundary conditions of the b-tensor
  call define_constraint ( mesh, input_probdefb, curve1=2, curve2=5, &
    discretization='collocation' )

  call problem_definition ( input_probdefb, mesh, problemb )


! create a vector subscript for the b-tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemb, bval, physqarr=[1] )
  call create_subscript ( mesh, problemb, cxx, degfd=1, vec=3 )
  call create_subscript ( mesh, problemb, cxy, degfd=2, vec=3 )
  call create_subscript ( mesh, problemb, cyy, degfd=3, vec=3 )

! create a vector for conformation tensor for post processing

  call create_vector ( problemb, ctensor, vec=3 )

! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, soln, solnm1, rhsd, sol_hat )


! fill solution vector with essential boundary conditions

  sol%u = 0

! initialize solution at tn to zero

  soln%u = 0

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! fill oldvectors

  oldvectors%s(1)%p => sol_hat  ! alpha0*un+alpha1*un-1


! create system matrix of gradient/velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! create system vectors (solution and right-hand side) for conformation and
! initialize vectors with steady stress and a small random perturbation.

  call create ( problemb, solb, solbm1, rhsb )

  do i = 1, nmodes
    solb(1,i)%u = 1 ! initial bxx
    solb(2,i)%u = 0 ! initial bxy
    solb(3,i)%u = 0 ! initial byx
    solb(4,i)%u = 1 ! initial byy
  end do


! problem definition for projected "c=b*b^T" of the b-tensor

  call create_input_probdef ( mesh, input_probdef_proj, nvec=1, nphysq=1 )

  input_probdef_proj%vec_elementdof(1)%a = &
      reshape ( [ 1,0,1,0,1,0,1,0,0 ], &
                   [9,1] )

  input_probdef_proj%physq = [1]
  input_probdef_proj%probnr = 3

  call problem_definition ( input_probdef_proj, mesh, problem_proj )

  call create ( problem_proj, solc_proj, rhsc_proj )

  do i = 1, nmodes
    solc_proj(1,i)%u = 1   ! initial cxx
    solc_proj(2,i)%u = 0   ! initial cxy
    solc_proj(3,i)%u = 1   ! initial cyy
  end do


! create system matrix for b-tensor problem

  call create_sysmatrix_structure_base ( sysmatrixb, mesh, problemb )
  call create_sysmatrix_structure_constraint ( sysmatrixb, mesh, problemb )
  call finalize_sysmatrix_structure ( sysmatrixb )

  call create_sysmatrix_data ( sysmatrixb )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=1, nsysvec2=3, nprob=3 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s2(1)%p => solb
  oldvectors_ve%s2(2)%p => solbm1
  oldvectors_ve%s2(3)%p => solc_proj
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemb
  oldvectors_ve%p(3)%p => problem_proj


! create and build system matrix for projection problem
! NOTE matrix remains constant and needs to be build once.

  call create_sysmatrix_structure ( sysmatrix_proj, mesh, problem_proj, &
    symmetric=.true. )
  call create_sysmatrix_data ( sysmatrix_proj )

  call build_system ( mesh, problem_proj, sysmatrix_proj, &
    m2sysvector=rhsc_proj, elemsub=c_projection_elem, &
    oldvectors=oldvectors_ve, coefficients=coefficients, &
    buildvector=.false. )

  call check ( sysmatrix_proj )


! open monitor file

  open ( unit=13, file='out', recl=300 )


! time stepping

  do step = 1, numtimesteps

    if ( step >= 2 ) then
      coefficients%i(22) = timeint2
    end if

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG


!   instationary inertia term (convection not included)

    if ( inertia == 1 ) then

!     instationary term ( rho du/dt )

      if ( step == 1 ) then
        sol_hat%u(vel%s) = soln%u(vel%s)
        mfac = 1
      else
        sol_hat%u(vel%s) = alpha0 * soln%u(vel%s) + alpha1 * solnm1%u(vel%s)
        mfac = gamma0
      end if

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=inertia_elem_dudt, coefficients=coefficients, &
        oldvectors=oldvectors, physqcol=[physqvel], physqrow=[physqvel], &
        factormat=mfac, addmatvec=.true. )

    end if


!   c=b*b^T projection for cn in momentum balance

    if ( cproj ) &
              call solve_projection ( solc_proj, rhsc_proj, c_projection_elem )


!   build implicit terms of CE with rhs in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      addmatvec=.true., coefficients=coefficients )


    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )


!   build (assemble) matrix and vector for b-tensor problem

    if ( coefficients%i(22) == timeint1 ) then

      call build_system ( mesh, problemb, sysmatrixb, m2sysvector=rhsb, &
        elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    else

      call build_system ( mesh, problemb, sysmatrixb, m2sysvector=rhsb, &
        elemsub=ce_supg_elem_implicit_2nd_order, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    end if


!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problemb, sysmatrixb, &
      m2sysvector=rhsb, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

    call check ( sysmatrixb )

    call copy ( solb, solbm1 )


!   solve conformation and keep LU decomposition in the loop over components

    solver_options_b%real_storage=rs_b
    solver_options_b%integer_storage=is_b

    do m = 1, nmodes
      do icomp = 1, ncompb
        call solve_system_ma41 ( sysmatrixb, rhsb(icomp,m), solb(icomp,m), &
          lub, solver_options=solver_options_b  )
      end do
    end do

    call delete ( lub )  ! remove LU decomposition and rebuild next time step


!   reinitialize b to b' = sqrt(c)

    call reinitialize_b


!   copy data for the next time step

    call copy ( soln, solnm1 )
    call copy ( sol, soln )


!   write monitor data

    call derive_vector ( mesh, problemb, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    write(13, fmt=*) step, step * deltat, &
      sol%u(press%s(1))-sol%u(press%s(2)), &  ! pressure difference
      maxval(ctensor%u(cxx%s)), maxval(ctensor%u(cxy%s)), &
      maxval(ctensor%u(cyy%s)), maxval(sol%u(velx%s))

  end do


! close monitor data file

  close(unit=13)


! write data for post-processing

  call write_mesh ( mesh, filename='mesh.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh, input_probdefb, filename='probdefc.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol%u
  write(10) ((solb(i,j)%u, i=1,ncompb), j=1,nmodes)

  close(unit=10)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, soln, solnm1, rhsd )
  call delete ( sysmatrix )
  call delete ( oldvectors_ve )
  call delete ( problemb )
  call delete ( input_probdefb )
  call delete ( sysmatrixb )
  call delete ( solb, solbm1, rhsb )

  if ( cproj ) then
    call delete ( sysmatrix_proj )
    call delete ( problem_proj )
    call delete ( input_probdef_proj )
    call delete ( solc_proj, rhsc_proj )
    call delete ( lu_proj )
  end if

  call delete ( coefficients )
  call delete ( bval )

contains

  subroutine build_vpG

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, physqrow=[2,3], physqcol=[2,3], &
      coefficients=coefficients )

!   DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_elem, addmatvec=.true., &
      physqrow=[1,2], physqcol=[1,2], coefficients=coefficients )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[1], physqcol=[3], zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[3], physqcol=[1], zeromatvec=.true. )

!   periodical condition on velocities and gradients

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_elem_conn, addmatvec=.true., &
      coefficients=coefficients )

!   periodical condition on gradients

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true. )

!   imposed flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=3, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients )

  end subroutine build_vpG



  subroutine solve_projection ( sol_proj, rhs_proj, elemsub_proj )

    type(sysvector_t), dimension(:,:), intent(inout) :: sol_proj

    type(sysvector_t), dimension(:,:), intent(inout) :: rhs_proj

    interface
      subroutine elemsub_proj ( mesh, problem, elgrp, elem, matrix, vector, &
        first, last, coefficients, oldvectors, elemmat, elemvec )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: elgrp, elem
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat
        real(dp), intent(out), dimension(:) :: elemvec
      end subroutine elemsub_proj
    end interface

    type(solver_options_ma57_t) :: solver_options_ma57

    integer :: i, m

!   build vector only (matrix is constant)

    call build_system ( mesh, problem_proj, sysmatrix_proj, &
      m2sysvector=rhs_proj, elemsub=elemsub_proj, &
      oldvectors=oldvectors_ve, coefficients=coefficients, &
      buildmatrix=.false. )

    ! MA57 solver storage
    solver_options_ma57%integer_storage = 1.3
    solver_options_ma57%real_storage    = 1.3

!   LU decomposition is done in the first call only

    do m = 1, nmodes
      do i = 1, ncompb-1

        call add_effect_of_essential_to_rhs ( problem_proj, sysmatrix_proj, &
           sol_proj(i,m), rhs_proj(i,m) )

        call solve_system_ma57 ( sysmatrix_proj, rhs_proj(i,m), &
           sol_proj(i,m), lu_proj, solver_options=solver_options_ma57 )

      end do
    end do

  end subroutine solve_projection


! reinitialize b to b = sqrt(c)

  subroutine reinitialize_b

    real(dp) :: b(size(bval%s),ncompb), bm1(size(bval%s),ncompb)
    real(dp) :: RT(size(bval%s),2,2)
    integer :: i

!   obtain the solution of b at current time step

    do i = 1,ncompb
      b(:,i) = solb(i,1)%u(bval%s)
    end do

!   check skew norm and perform reinitialization if exceeded

    if ( any ( skew_norm_2D_b(b) >= critval_reinit ) ) then

!     determine the square root of b*b^T

      call sqrtc_2D_b ( b, RT=RT )

!     obtain the solution bn-1 at previous time step

      do i = 1,ncompb
        bm1(:,i) = solbm1(i,1)%u(bval%s)
      end do

!     rotate bn-1 according to b

      call rotate_2D_b ( bm1, RT )

!     store b and bn-1 in solution vectors

      do i = 1,ncompb
        solb(i,1)%u(bval%s) = b(:,i)
        solbm1(i,1)%u(bval%s) = bm1(:,i)
      end do

    end if

  end subroutine reinitialize_b

end program poiseuille3
