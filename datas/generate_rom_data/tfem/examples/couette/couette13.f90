! Viscoelastic problem on a unit square
! linear Couette flow
! G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! first-order time integration
! NOTE: similar to couette3, but omitting DEVSS.

program couette13

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
    cintpl = 4,     &  ! Q1 conformation
    physqvel = 1,   & ! physical quantity nr of the velocities
    physqpress = 2, & ! physical quantity nr of the pressures
    gauss = 3,      & ! 3x3 integration of quads
    ncompc = 3,     & ! number of conformation tensor components
    ncompg = 4,     & ! number of velocity gradient tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 2         ! UCM/Oldroyd-B


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc
  type(problem_t), target :: problem, problemc
  type(sysmatrix_t) :: sysmatrix, sysmatrixc
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vely, cval, cvalc2, cvalc5
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc
  type(lu_ma41_t) :: luc
  type(solver_options_ma41_t) :: solver_options_u, solver_options_c

  type(input_probdef_t) :: input_probdef_grad
  type(problem_t) :: problem_grad
  type(sysmatrix_t) :: sysmatrix_grad
  type(sysvector_t) :: sol_grad, rhsd_grad(ncompg)
  type(oldvectors_t) :: oldvectors_grad
  type(vector_t), target :: gradients
  type(lu_ma57_t) :: lu_g


! variables

  integer :: &
    nx=20,               & ! number of elements in x
    ny=20,               & ! number of elements in y
    timeint = 1,         & ! (first-order) time integration
    numtimesteps = 2,    & ! number of time steps
    logc = 0               ! standard scheme or log transformation

  real(dp) :: &
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 1.0_dp,   & ! relaxation time
    deltat = 5.e-3_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.0_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 1.4_dp       ! integer_storage for conformation LU (HSL)

  logical :: first
  integer :: icomp, step, i
  integer, dimension(4) :: vertices=[1,3,5,7]
  real(dp) :: alpha, G, csteady(3), lsteady(2), eigv(2,2)


! namelist for input of variables; read from standard input

  namelist /comppar/ nx, ny, timeint, numtimesteps, logc, eta_s, eta_p, &
    lambda, deltat, beta, rs_gup, is_gup, rs_c, is_c

  read ( unit=*, nml=comppar )


! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+2*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     0,         gauss,  &
      gauss,    cintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      logc,     timeint,    ( 0, i = 23, 150 )  &
    ]

  coefficients%i(60) = 1  ! Separate vector for velocity gradient (No DEVSS).

  coefficients%r = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
      0._dp, 0._dp,  deltat,   beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G,     lambda &
    ]


! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [ 2,2,2,2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,1,0,0,    &  ! pressure
                  4,0,4,0,4,0,4,0,0,    &  ! G
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                  [9,4] )

  input_probdef%physq = [1,2]
  input_probdef%probnr = 1

! Dirichlet boundary conditions

! velocity on lower boundary
  call define_essential ( mesh, input_probdef, curve1=1, physq=physqvel )
! velocity on upper boundary
  call define_essential ( mesh, input_probdef, curve1=3, physq=physqvel )
! pressure in point 1
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! constraints for periodical boundary conditions

! velocities
  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=2, curve2=5, discretization='collocation', &
    exclude=3 )

  call problem_definition ( input_probdef, mesh, problem )


! create a vector subscript for the vertical velocity for post processing

  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=2, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a =   &
      reshape ( [ 1,0,1,0,1,0,1,0,0,    &  ! c
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar for plotting
                  [9,2] )

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

! constraint for periodical boundary conditions of the conformation
  call define_constraint ( mesh, input_probdefc, curve1=2, curve2=5, &
    discretization='collocation' )

  call problem_definition ( input_probdefc, mesh, problemc )


! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemc, cval, physqarr=[1] )
  call create_subscript ( mesh, problemc, cvalc2, physqarr=[1], curves=[2] )
  call create_subscript ( mesh, problemc, cvalc5, physqarr=[1], curves=[5] )

! create system vectors for velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )

! create gradient vector

  call create ( problem, gradients, vec=3 )


! fill solution vector with essential boundary conditions

! set velocity = 0 on lower boundary
  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, physq=physqvel, value=0._dp )
! set velocity = (1,0) on upper boundary
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=physqvel, degfd=1, value=1._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=physqvel, degfd=2, value=0._dp )
! set pressure level = 0 in lower left corner
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqpress, value=0._dp )


! create system matrix of velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! create system vectors (solution and right-hand side) for conformation and
! initialize vectors with steady stress and a small random perturbation.

  call create ( problemc, solc, rhsc )

! steady state solution

  csteady(1) = 1 + 2 * lambda ** 2
  csteady(2) = lambda
  csteady(3) = 1

  if ( logc == 0 ) then ! standard
    call random_number (solc(1,1)%u)
    solc(1,1)%u = csteady(1) + 1e-3*solc(1,1)%u  ! initial cxx
    call random_number (solc(2,1)%u)
    solc(2,1)%u = csteady(2) + 1e-3*solc(2,1)%u  ! initial cxy
    call random_number (solc(3,1)%u)
    solc(3,1)%u = csteady(3) + 1e-3*solc(3,1)%u   ! initial cyy
  else if ( logc == 1 ) then ! log scheme
    call eig2x2 ( csteady, lsteady, eigv )
    lsteady = log ( lsteady )
    call inveig2x2 ( csteady, lsteady, eigv )
    call random_number (solc(1,1)%u)
    solc(1,1)%u = csteady(1) + 1e-3*solc(1,1)%u   ! initial sxx
    call random_number (solc(2,1)%u)
    solc(2,1)%u = csteady(2) + 1e-3*solc(2,1)%u   ! initial sxy
    call random_number (solc(3,1)%u)
    solc(3,1)%u = csteady(3) + 1e-3*solc(3,1)%u   ! initial syy
  end if
  do i = 1, 3
    solc(i,1)%u(cvalc2%s) = solc(i,1)%u(cvalc5%s) ! make disturbance periodic
  end do

! create system matrix for conformation problem

  call create_sysmatrix_structure_base ( sysmatrixc, mesh, problemc )
  call create_sysmatrix_structure_constraint ( sysmatrixc, mesh, problemc )
  call finalize_sysmatrix_structure ( sysmatrixc )

  call create_sysmatrix_data ( sysmatrixc )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=1, nsysvec2=1, nprob=2, &
    nvec=3 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc
  oldvectors_ve%v(3)%p => gradients

! define gradient projection problem

  call proj_grad_definition


! open monitor file

  open ( unit=13, file='out', recl=300 )


  first = .true.

! time stepping

  do step = 1, numtimesteps

    if ( step >= 2 ) first = .false.

!   build (assemble) matrix/vector for velocity/pressure problem

    call build_vp

!   build implicit terms of CE with rhs in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve, physqrow=[physqvel], physqcol=[physqvel], &
      addmatvec=.true., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )


!   build and solve velocity projection problem

    call build_and_solve_proj_grad ( first )


!   build (assemble) matrix and vector for conformation problem

    call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
      elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
      coefficients=coefficients )

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problemc, sysmatrixc, &
      m2sysvector=rhsc, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

    call check ( sysmatrixc )


!   solve conformation and keep LU decomposition in the loop over components

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do icomp = 1, ncompc
      call solve_system_ma41 ( sysmatrixc, rhsc(icomp,1), solc(icomp,1), luc, &
        solver_options=solver_options_c  )
    end do

    call delete ( luc )  ! remove LU decomposition and rebuild next time step



!   write monitor data

    write(unit=13,fmt=*) &
      step, step*deltat, sqrt( &
      sum( ( solc(1,1)%u(cval%s) - csteady(1) )**2 + &
           ( solc(2,1)%u(cval%s) - csteady(2) )**2 + &
           ( solc(3,1)%u(cval%s) - csteady(3) )**2 ) / 3 / size(cval%s) ), &
      minval(solc(1,1)%u(cval%s)), &
      maxval(sol%u(vely%s)), minval(sol%u(vely%s))

  end do


! close monitor data file

  close(unit=13)


! write data for post-processing

  call write_mesh ( mesh, filename='mesh.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh, input_probdefc, filename='probdefc.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol%u
  write(10) (solc(i,1)%u, i=1,3)

  close(unit=10)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( oldvectors_ve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( solc, rhsc )

  call delete ( coefficients )
  call delete ( vely, cval, cvalc2, cvalc5 )

  call delete ( input_probdef_grad )
  call delete ( problem_grad )
  call delete ( sysmatrix_grad )
  call delete ( sol_grad )
  call delete ( rhsd_grad )
  call delete ( oldvectors_grad )
  call delete ( gradients )

  call delete ( lu_g )

contains

  subroutine build_vp

!   build (assemble) matrix and vector for velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients )

!   periodical condition on velocities
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_constr_node_conn, addmatvec=.true. )

  end subroutine build_vp

  subroutine proj_grad_definition

!   problem definition for gradients

    call create_input_probdef ( mesh, input_probdef_grad, nvec=2, nphysq=1 )

    input_probdef_grad%vec_elementdof(1)%a(:,1) = 0
    input_probdef_grad%vec_elementdof(1)%a(vertices,1) = 1  ! gradient component

    input_probdef_grad%physq = [1]
    input_probdef_grad%probnr = 3

    call problem_definition ( input_probdef_grad, mesh, problem_grad )

!   create system vectors (solution and right-hand side)

    call create ( problem_grad, sol_grad )
    call create ( problem_grad, rhsd_grad )

!   create system matrix

    call create_sysmatrix_structure ( sysmatrix_grad, mesh, problem_grad, &
      symmetric=.true. )

    call create_sysmatrix_data ( sysmatrix_grad )

!   oldvectors

    call create_oldvectors ( oldvectors_grad, nsysvec=1, nprob=1 )

    oldvectors_grad%s(1)%p => sol
    oldvectors_grad%p(1)%p => problem

  end subroutine proj_grad_definition


  subroutine build_and_solve_proj_grad ( first )

    use set_optional_m

!   optional argument to indicate it is called the first time so the matrix
!   is also built. NOTE: the LU decomposition lu_g is kept, so calling this
!   routine again with first=.true. only makes sense if lu_g is deleted before
!   the call.
!   default=.false.
    logical, optional, intent(in) :: first

    integer :: i

    if ( set_optional ( variable=first, default=.false. ) ) then

!     build (assemble) matrix and vector from elements

      call build_system ( mesh, problem_grad, sysmatrix_grad, &
        msysvector=rhsd_grad, elemsub=gradient_from_velocity_elem, &
        coefficients=coefficients, oldvectors=oldvectors_grad )

    else

!     build (assemble) vector from elements

      call build_system ( mesh, problem_grad, sysmatrix_grad, &
        msysvector=rhsd_grad, elemsub=gradient_from_velocity_elem, &
        coefficients=coefficients, oldvectors=oldvectors_grad, &
        buildmatrix=.false. )

    end if

!   solve

    do i = 1, ncompg

      call add_effect_of_essential_to_rhs ( problem_grad, sysmatrix_grad, &
        sol_grad, rhsd_grad(i) )

      call solve_system_ma57 ( sysmatrix_grad, rhsd_grad(i), sol_grad, lu=lu_g )

      call transfer_data ( mesh, problem_grad, problem, sysvector1=sol_grad, &
        vector2=gradients, degfd2=[i] )

    end do

  end subroutine build_and_solve_proj_grad

end program couette13
