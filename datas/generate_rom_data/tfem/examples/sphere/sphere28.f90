! Axisymmetrical elastoviscoplastic problem: a sphere moving with velocity U
! on the center line of a cylindrical container with a bottom (at infinity).
! The sphere is fixed and the cylinder wall is moved with velocity -U
! (moving frame).
! At the ends fluid is moving with velocity -U as well, representing zero
! flux in the stationary frame (bottom). At the right inflow wall the stress is
! assumed to be zero.
! Dirichlet boundary conditions for the velocities are used.
! DEVSS-G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! second-order time integration
! NOTE: this is similar to sphere27, but now using the CDT formulation.


program sphere28

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,     &  ! Q2 velocities
    pintpl = 4,     &  ! Q1 pressures
    gintpl = 4,     &  ! Q1 gradients
    bintpl = 4,     &  ! Q1 b-tensor
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    gauss = 3,      & ! 3x3 integration of quads
    ncompb = 5,     & ! number of b-tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    coorsys = 1,    & ! axisymmetric coordinate system
    model = 2,      & ! UCM model/Oldroyd-B
    alam_model = 3    ! 0:none, 1:elastic, 2: Saramito1, 3: Saramito2


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, solm1
  type(sysvector_t) :: rhsd
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options_u

  type(input_probdef_t) :: input_probdefb
  type(problem_t), target :: problemb
  type(sysmatrix_t) :: sysmatrixb
  type(oldvectors_t) :: oldvectors_ve
  type(subscript_t) :: bval
  type(subscriptvec_t) :: cxx, cxy, cyy, ctt
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
    cproj = .true. ! projection of c=b*b^T for cn in momentum balance

 integer :: &
    timeint1 = 1,        & ! (first-order) time integration (first step)
    timeint2 = 6,        & ! (second-order) time integration
    numtimesteps = 2,    & ! number of time steps
    restart = 0            ! do a restart

  real(dp) :: &
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 1.0_dp,   & ! relaxation time
    tau_y = 1.0_dp,    & ! yield stress
    K = 1.0_dp,        & ! power-law coefficient
    n = 0.5_dp,        & ! power-law exponent
    deltat = 5.e-3_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    critval_reinit = 0.0_dp,  & ! criterion for reinitializing b
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_b  = 1.0_dp,    & ! real_storage for the b-tensor LU (HSL)
    is_b  = 1.4_dp,    & ! integer_storage for b-tensor LU (HSL)
    U = 1._dp            ! velocity of the sphere

  integer :: icomp, step, i
  real(dp) :: alpha, G


! namelist for input of variables; read from standard input

  namelist /comppar/ timeint1, timeint2, numtimesteps, eta_s, eta_p, &
    lambda, tau_y, K, n, deltat, beta, rs_gup, is_gup, rs_b, is_b, U, restart

  read ( unit=*, nml=comppar )


! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+5*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gauss,    bintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      0,        timeint1,   coorsys, ( 0, i = 24, 150 )  &
    ]

  coefficients%i(71) = 1  ! CDT variant for b-formulation
  if ( cproj ) coefficients%i(72) = 1  ! use c projection in momentum balance
  coefficients%i(80) = alam_model ! adapted lambda

  coefficients%r = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
      0._dp, 0._dp,  deltat,   beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G, lambda, tau_y, K, n &
    ]


  call write_coefficients ( coefficients, filename='coefficients.out' )


! read mesh

  call read_mesh_gmsh ( mesh, filename='mesh1.msh', ndim=2 )

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=3 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [ 4,0,4,0,4,0,4,0,0,    &  ! G
                  2,2,2,2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,1,0,0,    &  ! pressure
                  1,1,1,1,1,1,1,1,1,    &  ! scalar, such as vorticity
                  4,4,4,4,4,4,4,4,4 ], &  ! tensor
                  [9,5] )

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1


! define essential boundaries

! center line
  call define_essential ( mesh, input_probdef, curve1=1, curve2=2, physq=2, &
    degfd=[0,1] )
  call define_essential ( mesh, input_probdef, curve1=7, curve2=8, physq=2, &
    degfd=[0,1] )
! sphere
  call define_essential ( mesh, input_probdef, curve1=3, curve2=6, physq=2 )
! ends and cylinder wall
  call define_essential ( mesh, input_probdef, curve1=9, curve2=14, physq=2 )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=3 )

  call problem_definition ( input_probdef, mesh, problem )


! problem definition b-tensor

  call create_input_probdef ( mesh, input_probdefb, nvec=4, nphysq=1 )

  input_probdefb%vec_elementdof(1)%a =   &
      reshape ( [ 1,0,1,0,1,0,1,0,0,    &  ! b
                  1,1,1,1,1,1,1,1,1,    &  ! scalar for plotting
                  4,4,4,4,4,4,4,4,4,    &  ! tensor
                  5,5,5,5,5,5,5,5,5 ], &  ! tensor
                  [9,4] )

  input_probdefb%physq = [1]
  input_probdefb%probnr = 2

! essential bc for b-tensor at the right inflow boundary
  call define_essential ( mesh, input_probdefb, curve1=9 )

  call problem_definition ( input_probdefb, mesh, problemb )

! create vector subscripts for the conformation tensor components

  call create_subscript ( mesh, problemb, bval, physqarr=[1] )
  call create_subscript ( mesh, problemb, cxx, degfd=1, vec=3 )
  call create_subscript ( mesh, problemb, cxy, degfd=2, vec=3 )
  call create_subscript ( mesh, problemb, cyy, degfd=3, vec=3 )
  call create_subscript ( mesh, problemb, ctt, degfd=4, vec=3 )

! create a vector for conformation and b-tensor
! tensor for post processing

  call create_vector ( problemb, ctensor, vec=3 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, solm1, rhsd )


! fill solution vector with essential boundary conditions

  sol%u = 0

  call fill_sysvector ( mesh, problem, sol, &
    curve1=9, curve2=14, physq=2, degfd=1, value=-U )


! create system matrix of gradient/velocity/pressure problem

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )


! create system vectors (solution and right-hand side) for b-tensor and
! initialize vectors with zero stress

  call create ( problemb, solb, solbm1, rhsb )

  solb(1,1)%u = 1 ! initial bzz
  solb(2,1)%u = 0 ! initial bzr
  solb(3,1)%u = 0 ! initial brz
  solb(4,1)%u = 1 ! initial brr
  solb(5,1)%u = 1 ! initial btt


! problem definition for projected "c=b*b^T" of the b-tensor

  call create_input_probdef ( mesh, input_probdef_proj, nvec=1, nphysq=1 )

  input_probdef_proj%vec_elementdof(1)%a = &
      reshape ( [ 1,0,1,0,1,0,1,0,0 ], &
                   [9,1] )

  input_probdef_proj%physq = [1]
  input_probdef_proj%probnr = 3

  call problem_definition ( input_probdef_proj, mesh, problem_proj )

  call create ( problem_proj, solc_proj, rhsc_proj )

  solc_proj(1,1)%u = 1   ! initial czz
  solc_proj(2,1)%u = 0   ! initial czr
  solc_proj(3,1)%u = 1   ! initial crr
  solc_proj(4,1)%u = 1   ! initial ctt


! create system matrix for b-tensor problem

  call create_sysmatrix_structure ( sysmatrixb, mesh, problemb )

  call create_sysmatrix_data ( sysmatrixb )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec2=4, nprob=3 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s(2)%p => solm1
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


! restart: read solution from file

  if ( restart == 1 ) then

    open ( unit=10, form='unformatted', file='data.out' )

    read(10) sol%u
    read(10) (solb(i,1)%u, i=1,ncompb)

    close(unit=10)

    call fill_sysvector ( mesh, problem, sol, &
      curve1=9, curve2=14, physq=2, degfd=1, value=-U )

  end if


! open monitor file

  open ( unit=13, file='out', recl=300 )



! time stepping

  do step = 1, numtimesteps

    if ( step >= 2 ) then
      coefficients%i(22) = timeint2
    end if

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG


!   exps projection (for log conformation)

    if ( cproj ) &
             call solve_projection ( solc_proj, rhsc_proj, c_projection_elem )

!   build implicit terms

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2],&
      addmatvec=.true., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    call copy ( sol, solm1 )

!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )


!   build (assemble) matrix and vector for conformation problem

    if ( coefficients%i(22) == timeint1 ) then

      call build_system ( mesh, problemb, sysmatrixb, m2sysvector=rhsb, &
        elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    else

      call build_system ( mesh, problemb, sysmatrixb, m2sysvector=rhsb, &
        elemsub=ce_supg_elem_implicit_2nd_order, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    end if


    call check ( sysmatrixb )

    call copy ( solb, solbm1 )

    do icomp = 1, ncompb
      call add_effect_of_essential_to_rhs ( problemb, sysmatrixb, &
        solb(icomp,1), rhsb(icomp,1) )
    end do


!   solve b-tensor and keep LU decomposition in the loop over components

    solver_options_b%real_storage=rs_b
    solver_options_b%integer_storage=is_b

    do icomp = 1, ncompb
      call solve_system_ma41 ( sysmatrixb, rhsb(icomp,1), solb(icomp,1), lub, &
        solver_options=solver_options_b  )
    end do

    call delete ( lub )  ! remove LU decomposition and rebuild next time step


!   reinitialize b to b' = sqrt(c)

    call reinitialize_b


!   write monitor data

    call derive_vector ( mesh, problemb, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    write(unit=13,fmt=*) &
      step, step*deltat, maxval(ctensor%u(cxx%s)), &
                         maxval(ctensor%u(cxy%s)), &
                         maxval(ctensor%u(cyy%s)), &
                         maxval(ctensor%u(ctt%s))

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
  write(10) (solb(i,1)%u, i=1,ncompb)

  close(unit=10)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, solm1, rhsd )
  call delete ( sysmatrix )
  call delete ( oldvectors_ve )
  call delete ( problemb )
  call delete ( input_probdefb )
  call delete ( sysmatrixb )
  call delete ( solb, solbm1, rhsb )

  call delete ( coefficients )

  if ( cproj ) then
    call delete ( sysmatrix_proj )
    call delete ( problem_proj )
    call delete ( input_probdef_proj )
    call delete ( solc_proj, rhsc_proj )
    call delete ( lu_proj )
  end if


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


end program sphere28
