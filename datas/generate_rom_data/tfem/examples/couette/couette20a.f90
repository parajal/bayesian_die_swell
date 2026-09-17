! Viscoelastic problem on a unit square using b-formulation.
! linear Couette flow
! DEVSS-G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! L2-projection of c=b.b^T on the finite element approximation space in the
! momentum balance.
! second-order time integration
! Fully-implicit with Newton-Raphson for the b-tensor.
! NOTE: To find the reference steady state conformation tensor, the flow is
!       first run starting from the exact steady state without an imposed
!       disturbance until timestep dist_step.
! NOTE: This is similar to couette15, but now with implicit conformation.

program couette20a

  use tfem_m
  use viscoelastic_elements_m
  use math_defs_m, only: chol
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities
    pintpl = 4,     & ! Q1 pressures
    gintpl = 4,     & ! Q1 gradients
    bintpl = 4,     & ! Q1 contravariant deformation
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    gauss = 3,      & ! 3x3 integration of quads
    ncompb = 4,     & ! number of contravariant deformation tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 2         ! Oldroyd-B


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefb, input_probdef_proj
  type(problem_t), target :: problem, problemb, problem_proj
  type(sysmatrix_t) :: sysmatrix, sysmatrixb, sysmatrix_proj
  type(sysvector_t), target :: sol, solm1
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: bvalcmp(ncompb)
  type(subscript_t) :: vely, bval, bvalc2, bvalc5
  type(subscriptvec_t) :: cxx, cxy, cyy
  type(sysvector_t), dimension(nmodes), target :: solbn, solbnm1, solbiter
  type(sysvector_t) :: rhsb, dsolb
  type(sysvector_t), dimension(ncompb-1,nmodes), target :: solc_proj
  type(sysvector_t), dimension(ncompb-1,nmodes) :: rhsc_proj
  type(vector_t) :: ctensor, btensor
  type(solver_options_ma41_t) :: solver_options_u, solver_options_b
  type(lu_ma57_t) :: lu_proj


! variables

  logical :: &
    printscreen = .true. ! print on standard output

  integer :: &
    nx=20,               & ! number of elements in x
    ny=20,               & ! number of elements in y
    timeint1 = 8,        & ! (first-order) time integration (first step)
    timeint2 = 10,        & ! (second-order) time integration
    dist_step = 1,       & ! time step of disturbance
    maxnumiterations = 20, & ! maximum number of Newton-Raphson iterations
    numtimesteps = 2,    & ! number of time steps
    bvariant = 1           ! b-formulation:
                           ! 1: CDT 2: symmetric
                           ! 3: Cholesky 4: Cholesky with log
  real(dp) :: &
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 1.0_dp,   & ! relaxation time
    deltat = 5.e-3_dp, & ! time step
    epsconf = 1e-12_dp, & ! Newton-Raphson convergence threshold
    ampl = 0.0_dp,     & ! amplitude of initial disturbance
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_b  = 1.0_dp,    & ! real_storage for contravariant deformation LU (HSL)
    is_b  = 1.4_dp      ! integer_storage for contravariant deformation LU (HSL)

  integer :: step, i, iter, m
  integer :: vertices(4) = [1,3,5,7]
  real(dp) :: alpha, G, csteady(3), lsteady(2), eigv(2,2), scsteady(3)
  real(dp), dimension(:), allocatable :: tmpar

  logical :: &
    cproj = .true.         ! projection of c=b*b^T for cn in momentum balance



! namelist for input of variables; read from standard input

  namelist /comppar/ printscreen, nx, ny, timeint1, timeint2, &
    numtimesteps, maxnumiterations, dist_step, epsconf, &
    cproj, eta_s, eta_p, lambda, deltat, ampl, beta, &
    rs_gup, is_gup, rs_b, is_b, bvariant

  read ( unit=*, nml=comppar )


! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+2*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gauss,    bintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      0,        timeint1,   ( 0, i = 23, 150 )  &
    ]

  coefficients%i(61) = 0  ! projected G=0, direct velocity gradient=1 in CE
  coefficients%i(71) = bvariant  ! variant for b-formulation
  if ( cproj ) coefficients%i(72) = 1  ! use c projection in momentum balance
  coefficients%i(84) = 2  ! storage of b tensor: all components, separate modes
  coefficients%i(90) = 1  ! rotation reinitialization on element level

  coefficients%r = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
      0._dp, 0._dp,  deltat,   beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G,     lambda &
    ]

  coefficients%r(10) = 1 ! global scaling velocity for SUPG, when needed

! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=3 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [ 4,0,4,0,4,0,4,0,0,    &  ! G
                  2,2,2,2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,1,0,0,    &  ! pressure
                  1,1,1,1,1,1,1,1,1 ],  &  ! scalar, such as vorticity
                  [9,4] )

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

! Dirichlet boundary conditions

! velocity on lower boundary
  call define_essential ( mesh, input_probdef, curve1=1, physq=2 )
! velocity on upper boundary
  call define_essential ( mesh, input_probdef, curve1=3, physq=2 )
! pressure in point 1
  call define_essential ( mesh, input_probdef, point=1, physq=3 )

! constraints for periodical boundary conditions

! velocities
  call define_constraint ( mesh, input_probdef, &
    physq=2, curve1=2, curve2=5, discretization='collocation', exclude=3 )

! gradients
  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, curve2=5, discretization='collocation' )

  call problem_definition ( input_probdef, mesh, problem )


! create a vector subscript for the vertical velocity for post processing

  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )


! problem definition contravariant deformation tensor

  call create_input_probdef ( mesh, input_probdefb, nvec=4, nphysq=1 )

  input_probdefb%vec_elementdof(1)%a(vertices,1) = ncompb ! b
  input_probdefb%vec_elementdof(1)%a(:,2) = 1 ! scalar for plotting
  input_probdefb%vec_elementdof(1)%a(:,3) = 3 ! tensor for plotting
  input_probdefb%vec_elementdof(1)%a(:,4) = 4 ! tensor for plotting

  input_probdefb%physq = [1]
  input_probdefb%probnr = 2

! constraint for periodical boundary conditions of the contravariant deformation
  call define_constraint ( mesh, input_probdefb, curve1=2, curve2=5, &
    discretization='collocation' )

  call problem_definition ( input_probdefb, mesh, problemb )


! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemb, bval, physqarr=[1] )
  call create_subscript ( mesh, problemb, bvalc2, physqarr=[1], curves=[2] )
  call create_subscript ( mesh, problemb, bvalc5, physqarr=[1], curves=[5] )
  do i = 1, ncompb
    call create_subscript ( mesh, problemb, bvalcmp(i), physqarr=[1], degfd=i )
  end do
  call create_subscript ( mesh, problemb, cxx, degfd=1, vec=3 )
  call create_subscript ( mesh, problemb, cxy, degfd=2, vec=3 )
  call create_subscript ( mesh, problemb, cyy, degfd=3, vec=3 )

! create a vector for conformation and contravariant deformation
! tensor for post processing

  call create_vector ( problemb, ctensor, vec=3 )
  call create_vector ( problemb, btensor, vec=4 )

! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, solm1, rhsd )


! fill solution vector with essential boundary conditions

! set velocity = 0 on lower boundary
  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, physq=2, value=0._dp )
! set velocity = (1,0) on upper boundary
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=2, degfd=1, value=1._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=2, degfd=2, value=0._dp )
! set pressure level = 0 in lower left corner
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=3, value=0._dp )
  call copy(sol,solm1)


! create system matrix of gradient/velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! create system vectors (solution and right-hand side) for conformation and
! initialize vectors with steady stress and a small random perturbation.

  call create ( problemb, solbn, solbnm1, solbiter )
  call create ( problemb, rhsb, dsolb )

! exact steady state solution

  csteady(1) = 1 + 2 * lambda ** 2
  csteady(2) = lambda
  csteady(3) = 1

  if ( any ( bvariant == [1,2] ) ) then
    call eig2x2 ( csteady, lsteady, eigv )
    lsteady = sqrt ( lsteady )
    call inveig2x2 ( scsteady, lsteady, eigv )
    solbn(1)%u(bvalcmp(1)%s) = scsteady(1) ! initial bxx
    solbn(1)%u(bvalcmp(2)%s) = scsteady(2) ! initial bxy
    solbn(1)%u(bvalcmp(3)%s) = solbn(1)%u(bvalcmp(2)%s) ! initial byx
    solbn(1)%u(bvalcmp(4)%s) = scsteady(3) ! initial byy
  else if ( bvariant == 3 ) then
    scsteady = chol ( csteady )
    solbn(1)%u(bvalcmp(1)%s) = scsteady(1) ! initial bxx
    solbn(1)%u(bvalcmp(2)%s) = 0           ! initial bxy
    solbn(1)%u(bvalcmp(3)%s) = scsteady(2) ! initial byx
    solbn(1)%u(bvalcmp(4)%s) = scsteady(3) ! initial byy
  else if ( bvariant == 4 ) then
    scsteady = chol ( csteady )
    solbn(1)%u(bvalcmp(1)%s) = log(scsteady(1)) ! initial bxx
    solbn(1)%u(bvalcmp(2)%s) = 0                ! initial bxy
    solbn(1)%u(bvalcmp(3)%s) = scsteady(2)      ! initial byx
    solbn(1)%u(bvalcmp(4)%s) = log(scsteady(3)) ! initial byy
  end if
  solbn(1)%u(bvalc2%s) = solbn(1)%u(bvalc5%s) ! make disturbance periodic
  call copy(solbn,solbiter)

! problem definition for projected "c=b*b^T" of the contravariant deformation b

  call create_input_probdef ( mesh, input_probdef_proj, nvec=1, nphysq=1 )

  input_probdef_proj%vec_elementdof(1)%a = &
      reshape ( [ 1,0,1,0,1,0,1,0,0 ], &
                  [9,1] )

  input_probdef_proj%physq = [1]
  input_probdef_proj%probnr = 3

  call problem_definition ( input_probdef_proj, mesh, problem_proj )

  call create ( problem_proj, solc_proj, rhsc_proj )

  solc_proj(1,1)%u = 1   ! initial cxx
  solc_proj(2,1)%u = 0   ! initial cxy
  solc_proj(3,1)%u = 1   ! initial cyy


! create system matrix for conformation problem

  call create_sysmatrix_structure_base ( sysmatrixb, mesh, problemb )
  call create_sysmatrix_structure_constraint ( sysmatrixb, mesh, problemb )
  call finalize_sysmatrix_structure ( sysmatrixb )

  call create_sysmatrix_data ( sysmatrixb )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec1=3, nsysvec2=3, &
    nprob=3 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s(2)%p => solm1
  oldvectors_ve%s1(1)%p => solbn
  oldvectors_ve%s1(2)%p => solbnm1
  oldvectors_ve%s1(3)%p => solbiter
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
  open ( unit=14, file='iter.out', recl=300 )

  allocate ( tmpar(size(solbn(1)%u)) )

! time stepping

  do step = 1, numtimesteps

    write(14,*) 'step = ', step
    if ( printscreen) print *, 'step = ', step

   if ( step >= 2 ) then
      coefficients%i(22) = timeint2
    end if

!   disturbance

    if ( step == dist_step ) then

      !coefficients%i(61) = 1  ! projected G=0, direct velocity gradient=1 in CE

      csteady(1) = maxval(ctensor%u(cxx%s))
      csteady(2) = maxval(ctensor%u(cxy%s))
      csteady(3) = maxval(ctensor%u(cyy%s))

      call random_number (tmpar)
      solbn(1)%u(bvalcmp(1)%s) = &
          solbn(1)%u(bvalcmp(1)%s) + ampl*tmpar   ! initial bxx
      call random_number (tmpar)
      solbn(1)%u(bvalcmp(2)%s) = &
          solbn(1)%u(bvalcmp(2)%s) + ampl*tmpar   ! initial bxy
      call random_number (tmpar)
      solbn(1)%u(bvalcmp(3)%s) = &
          solbn(1)%u(bvalcmp(3)%s) + ampl*tmpar   ! initial bxy
      call random_number (tmpar)
      solbn(1)%u(bvalcmp(4)%s) = &
          solbn(1)%u(bvalcmp(4)%s) + ampl*tmpar   ! initial byy

      write(unit=13,fmt=*)

    end if


!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG


!   c=b*b^T projection for cn in momentum balance

    if ( cproj ) call solve_projection


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


    solver_options_b%real_storage=rs_b
    solver_options_b%integer_storage=is_b

    do m = 1, nmodes

      coefficients%i(85) = m   ! set mode number

      if ( nmodes > 1 ) then
        write(14,*) 'mode = ', m
        if ( printscreen) print *, 'mode = ', m
      end if

      iter = 0

      do

        iter = iter + 1

        if ( iter > maxnumiterations ) then
          write(*,'(3(a,i0/))') &
            ' Maximum number of iterations reached = ', &
            maxnumiterations, ' step = ', step, ' mode = ', m
          stop
        end if

!       build (assemble) matrix and vector for b tensor problem

        call build_system ( mesh, problemb, sysmatrixb, sysvector=rhsb, &
          elemsub=implicit_ce_supg_elem, oldvectors=oldvectors_ve, &
          coefficients=coefficients )

!       periodical condition on conformation tensor
        call build_system_constraint ( mesh, problemb, sysmatrixb, &
          sysvector=rhsb, elemsub=stokes_constr_node_conn, &
          addmatvec=.true. )

        call check ( sysmatrixb )

!       solve b tensor

        call solve_system_ma41 ( sysmatrixb, rhsb, dsolb, &
          solver_options=solver_options_b  )

        solbiter(m)%u(bval%s) = solbiter(m)%u(bval%s) + dsolb%u(bval%s)

        if ( printscreen) print *, maxval(abs(dsolb%u(bval%s)))
        write(14,*) maxval(abs(dsolb%u(bval%s)))

        if ( maxval(abs(dsolb%u(bval%s))) < epsconf ) exit

      end do

    end do

!   copy b-tensor solution to older time step for next time step

    call copy ( solbn, solbnm1 )
    call copy ( solbiter, solbn )

!   copy velocity solution to older time step for next time step

    call copy ( sol, solm1 )


!   conformation tensor

    call derive_vector ( mesh, problemb, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

!   write monitor data

    write(unit=13,fmt=*) &
      step, (step-dist_step+1)*deltat, sqrt( &
      sum( ( ctensor%u(cxx%s) - csteady(1) )**2 + &
           ( ctensor%u(cxy%s) - csteady(2) )**2 + &
           ( ctensor%u(cyy%s) - csteady(3) )**2 ) / 3 / size(cxx%s) ), &
      maxval(sol%u(vely%s)), minval(sol%u(vely%s)), &
      minval(ctensor%u(cxx%s)), &
      minval(ctensor%u(cxy%s)), &
      minval(ctensor%u(cyy%s))

  end do

  deallocate ( tmpar )

! close monitor data file

  close(unit=13)
  close(unit=14)


! write data for post-processing

  call write_mesh ( mesh, filename='mesh.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh, input_probdefb, filename='probdefc.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol%u
  write(10) solbn(1)%u

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
  call delete ( dsolb, rhsb )
  call delete ( solbn, solbnm1, solbiter )
  call delete ( sysmatrixb )
  call delete ( sysmatrix_proj )

  call delete ( problem_proj )
  call delete ( input_probdef_proj )
  call delete ( solc_proj, rhsc_proj )
  if ( cproj ) call delete ( lu_proj )


  call delete ( coefficients )
  call delete ( vely, bval, bvalc2, bvalc5 )
  do i = 1, ncompb
    call delete ( bvalcmp(i) )
  end do
  call delete ( cxx, cxy, cyy )

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
      elemsub=stokes_constr_node_conn, addmatvec=.true. )

  end subroutine build_vpG


! project c=b.b^T on discrete fem space

  subroutine solve_projection

    type(solver_options_ma57_t) :: solver_options_ma57

    integer :: i, m

!   build vector only (matrix is constant)

    call build_system ( mesh, problem_proj, sysmatrix_proj, &
      m2sysvector=rhsc_proj, elemsub=c_projection_elem, &
      oldvectors=oldvectors_ve, coefficients=coefficients, &
      buildmatrix=.false. )

    ! MA57 solver storage
    solver_options_ma57%integer_storage = 1.3
    solver_options_ma57%real_storage    = 1.3

!   LU decomposition is done in the first call only

    do m = 1, nmodes
      do i = 1, ncompb - 1

        call add_effect_of_essential_to_rhs ( problem_proj, sysmatrix_proj, &
           solc_proj(i,m), rhsc_proj(i,m) )

        call solve_system_ma57 ( sysmatrix_proj, rhsc_proj(i,m), &
           solc_proj(i,m), lu_proj, solver_options=solver_options_ma57 )

      end do
    end do

  end subroutine solve_projection


end program couette20a
