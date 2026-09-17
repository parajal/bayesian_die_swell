! Viscoelastic problem on a unit square (b formulation)
! linear Couette flow
! DEVSS-G/SUPG
! second-order time integration
! Fully-implicit with Newton-Raphson
! Optionally include Jacobian of SUPG test function.
! Optionally use L=G or L=(nabla u)^T in constitutive equation.
! NOTE: To find the reference steady state conformation tensor, the flow is
!       first solved for the exact steady state without an imposed
!       disturbance using Newton-Raphson.

program couette27

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,     &  ! Q2 velocities
    pintpl = 4,     &  ! Q1 pressures
    gintpl = 4,     &  ! Q1 gradients
    bintpl = 4,     &  ! Q1 conformation
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    physqb = 4,     & ! physical quantity nr of the conformation tensor (mode 1)
    gauss = 3,      & ! 3x3 integration of quads
    ncompb = 4,     & ! number of conformation tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 2,      & ! UCM/Oldroyd-B model
    directgrad = 0, & ! direct velocity gradient, 1: yes, 0: no
    bvariant = 1      ! b-formulation:
                      ! 1: CDT 2: symmetric
                      ! 3: Cholesky 4: Cholesky with log

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: soln, solnm1, soliter
  type(sysvector_t) :: dsol, rhsd
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscriptvec_t) :: cxx, cxy, cyy
  type(subscript_t) :: velx, vely, vel, pres, solsc, gxy
  type(subscript_t), dimension(ncompb,nmodes) :: bvalcmp
  type(subscript_t), dimension(nmodes) :: bval, bvalc2, bvalc5

  type(vector_t) :: ctensor
  type(solver_options_ma41_t) :: solver_options


! variables

  logical :: &
    JacobianSUPG = .true., & ! take Jacobian of SUPG test function into account
    printscreen = .true. ! print on standard output

  integer :: &
    nx=20,             & ! number of elements in x
    ny=20,             & ! number of elements in y
    timeint1 = 8,      & ! first-order time integration, first time step
    timeint2 = 10,     & ! second-order time integration after first time step
    dist_step = 1,     & ! time step of disturbance
    maxnumiterations = 20, & ! maximum number of Newton-Raphson iterations
    numtimesteps = 2     ! number of time steps

  real(dp) :: &
    eta_s = 0.1_dp,     & ! solvent viscosity
    eta_p = 1.0_dp,     & ! polymer viscosity
    lambda = 1.0_dp,    & ! relaxation time
    deltat = 0.4_dp,    & ! time step
    thetapar = 0.55_dp, & ! theta parameter in the theta method
    epsconf = 1e-12_dp, & ! Newton-Raphson convergence threshold
    ampl = 0.0_dp,     & ! amplitude of initial disturbance
    beta = 1.0_dp,      & ! upwinding parameter in the SUPG method
    rs    = 1.2_dp,     & ! real_storage for LU (HSL)
    is    = 1.6_dp        ! integer_storage for LU (HSL)

  integer :: step, i, iter, m
  real(dp) :: alpha, G, bsteady(4), csteady(4)
  real(dp), dimension(:), allocatable :: tmpar


! namelist for input of variables; read from standard input

  namelist /comppar/ printscreen, nx, ny, timeint1, timeint2, &
    numtimesteps, maxnumiterations, epsconf, eta_s, &
    eta_p, lambda, deltat, ampl, beta, rs, is

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
      0,        timeint1,    ( 0, i = 23, 150 )  &
    ]

  coefficients%i(61) = directgrad ! direct velocity gradient, 1: yes, 0: no
  coefficients%i(71) = bvariant ! variant for b-formulation
  coefficients%i(83) = physqb ! physical quantity of first mode
  coefficients%i(84) = 3  ! storage of b tensor: all modes in sysvector
  coefficients%i(86) = 2  ! sysvector number for velocity in SUPG = soln
  coefficients%i(87) = 1  ! sysvector number for iteration b tensor
  coefficients%i(88) = 2  ! sysvector number of b tensor at tn
  coefficients%i(90) = 1  ! rotation reinitialization on element level

  if ( JacobianSUPG ) then
    coefficients%i(86) = 1  ! sysvector number for velocity in SUPG = soliter
    coefficients%i(91) = 1  ! Jacobian of SUPG
  end if

  coefficients%r = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
      0._dp, 0._dp,  deltat,   beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G,     lambda &
    ]

  coefficients%r(28) = thetapar

! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=6, nphysq=4 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [ 4,0,4,0,4,0,4,0,0,    &  ! G
                  2,2,2,2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,1,0,0,    &  ! pressure
                  ncompb,0,ncompb,0,ncompb,0,ncompb,0,0, &  ! b
                  1,1,1,1,1,1,1,1,1,    &  ! scalar, such as vorticity
                  3,3,3,3,3,3,3,3,3 ], &  ! tensor
                  [9,6] )

  input_probdef%physq = [1,2,3,(3+m,m=1,nmodes)]
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

! b tensor (use collocation)
  do m = 1, nmodes
    call define_constraint ( mesh, input_probdef, &
      physq=physqb+m-1, curve1=2, curve2=5, discretization='collocation' )
  end do

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for solution (excluding contraint forces)

  call create_subscript ( mesh, problem, solsc, &
    physqarr=[physqgrad,physqvel,physqpress,(physqb+m-1,m=1,nmodes)] )

! create vector subscripts for the gradient xy

  call create_subscript ( mesh, problem, gxy, physqarr=[physqgrad], degfd=2 )

! create vector subscripts for the velocity

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )

! create vector subscript for the pressure

  call create_subscript ( mesh, problem, pres, physqarr=[physqpress] )

! create vector subscript for the b tensor

  do m = 1, nmodes
    call create_subscript ( mesh, problem, bval(m), physqarr=[physqb+m-1] )
    call create_subscript ( mesh, problem, bvalc2(m), physqarr=[physqb+m-1], &
      curves=[2] )
    call create_subscript ( mesh, problem, bvalc5(m), physqarr=[physqb+m-1], &
      curves=[5] )
    do i = 1, ncompb
      call create_subscript ( mesh, problem, bvalcmp(i,m), &
        physqarr=[physqb+m-1], degfd=i )
    end do
  end do
  call create_subscript ( mesh, problem, cxx, degfd=1, vec=6 )
  call create_subscript ( mesh, problem, cxy, degfd=2, vec=6 )
  call create_subscript ( mesh, problem, cyy, degfd=3, vec=6 )

! create a vector for conformation and contravariant deformation
! tensor for post processing

  call create_vector ( problem, ctensor, vec=6 )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, soln, solnm1, soliter )
  call create_sysvector ( problem, dsol, rhsd )


! fill solution vector with essential boundary conditions

  soln%u = 0
! set gxy = 1 on domain
  soln%u(gxy%s) = 1
! set velx = y on domain
  call fill_sysvector ( mesh, problem, soln, &
    node1=1, node2=mesh%nnodes, physq=2, degfd=1, func=func, funcnr=1 )
  dsol%u = 0

! initial solution b tensor

  call steady_state ( shearrate=1._dp, bsteady=bsteady, csteady=csteady )

  solnm1%u(bvalcmp(1,1)%s) = bsteady(1) ! initial bxx
  solnm1%u(bvalcmp(2,1)%s) = bsteady(2) ! initial bxy
  solnm1%u(bvalcmp(3,1)%s) = bsteady(3) ! initial byx
  solnm1%u(bvalcmp(4,1)%s) = bsteady(4) ! initial byy

  allocate ( tmpar(size(bvalcmp(1,1)%s)) )

  call random_number (tmpar)
  soln%u(bvalcmp(1,1)%s) = bsteady(1) + ampl*tmpar   ! initial bxx
  call random_number (tmpar)
  soln%u(bvalcmp(2,1)%s) = bsteady(2) + ampl*tmpar   ! initial bxy
  call random_number (tmpar)
  soln%u(bvalcmp(3,1)%s) = bsteady(3) + ampl*tmpar   ! initial byx
  call random_number (tmpar)
  soln%u(bvalcmp(4,1)%s) = bsteady(4) + ampl*tmpar   ! initial byy
  soln%u(bvalc2(1)%s) = soln%u(bvalc5(1)%s) ! make disturbance periodic

  call copy(soln,soliter)

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=3, nprob=1 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => soliter
  oldvectors_ve%s(2)%p => soln
  oldvectors_ve%s(3)%p => solnm1
  oldvectors_ve%p(1)%p => problem  ! needed for implicit_ce_supg_elem


! open monitor file

  open ( unit=13, file='out', recl=300 )
  open ( unit=14, file='iter.out', recl=300 )

  solver_options%real_storage=rs
  solver_options%integer_storage=is


! time stepping

  do step = 1, numtimesteps

    write(14,*) 'step = ', step
    if ( printscreen) print *, 'step = ', step

    coefficients%i(22) = timeint2

    iter = 0

    do

      iter = iter + 1

      if ( iter > maxnumiterations ) then
        write(*,'(3(a,i0/))') &
          ' Maximum number of iterations reached = ', &
          maxnumiterations, ' step = ', step
        stop
      end if

!     build (assemble) matrix/vector for gradient/velocity/pressure part

      call build_vpG

!     build (assemble) matrix/vector for conformation part

      do m = 1, nmodes

        coefficients%i(85) = m   ! set mode number

        call build_b

      end do

!     set to zero off-diagonal conformation blocks

      call build_offdiag_b

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, dsol, rhsd )

      call check ( sysmatrix )

!     solve

      call solve_system_ma41 ( sysmatrix, rhsd, dsol, &
        solver_options=solver_options  )

      soliter%u(solsc%s) = soliter%u(solsc%s) + dsol%u(solsc%s)

      if ( printscreen) print *, maxval(abs(dsol%u(solsc%s))), &
        maxval(abs(dsol%u(vel%s))), maxval(abs(dsol%u(pres%s))), &
        maxval(abs(dsol%u(bval(1)%s)))
      write(14,*) maxval(abs(dsol%u(solsc%s))), &
        maxval(abs(dsol%u(vel%s))), maxval(abs(dsol%u(pres%s))), &
        maxval(abs(dsol%u(bval(1)%s)))

      if ( maxval(abs(dsol%u(solsc%s))) < epsconf ) exit

    end do

!   copy solution to older time step for next time step

    call copy ( soln, solnm1 )
    call copy ( soliter, soln )

!   conformation tensor

    call derive_vector ( mesh, problem, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

!   write monitor data

    write(unit=13,fmt=*) &
      step, (step-dist_step+1)*deltat, sqrt( &
      sum( ( ctensor%u(cxx%s) - csteady(1) )**2 + &
           ( ctensor%u(cxy%s) - csteady(2) )**2 + &
           ( ctensor%u(cyy%s) - csteady(3) )**2 ) / 3 / size(cxx%s) ), &
      maxval(soln%u(vely%s)), minval(soln%u(vely%s)), &
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

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) soln%u

  close(unit=10)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( soln, solnm1, soliter, dsol, rhsd )
  call delete ( sysmatrix )
  call delete ( oldvectors_ve )

  call delete ( coefficients )
  do m = 1, nmodes
    do i = 1, ncompb
      call delete ( bvalcmp(i,m), bval(m) )
    end do
  end do
  call delete ( solsc, gxy )

contains

  subroutine build_vpG

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients, &
      physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress] )

!   stokes add (q,divu) to rhs of mass balance
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_rhs_divu, oldvectors=oldvectors_ve, &
      coefficients=coefficients, addmatvec=.true., &
      physqrow=[physqpress], physqcol=[physqvel], buildmatrix=.false. )

!   DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_elem, coefficients=coefficients, addmatvec=.true., &
      physqrow=[physqgrad,physqvel], physqcol=[physqgrad,physqvel] )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqgrad], physqcol=[physqpress], &
      zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqpress], physqcol=[physqgrad], &
      zeromatvec=.true. )

!   build DEVSS-G right-hand side

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_rhs_div, oldvectors=oldvectors_ve, &
      coefficients=coefficients, addmatvec=.true., &
      buildmatrix=.false., physqrow=[physqvel] )

!   build stokes right-hand side

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_rhs_divsigma, oldvectors=oldvectors_ve, &
      coefficients=coefficients, addmatvec=.true., &
      buildmatrix=.false., physqrow=[physqvel] )

!   build -div(tau) part in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_divtau, oldvectors=oldvectors_ve, &
      physqrow=[physqvel], physqcol=[(physqb+m-1,m=1,nmodes)], &
      addmatvec=.true., coefficients=coefficients )

!   periodical condition on velocities

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients )

!   periodical condition on gradients

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients )

  end subroutine build_vpG


  subroutine build_b

!   build (assemble) matrix and vector for conformation tensor mode m

    if ( JacobianSUPG ) then

!     velocity + diagonal block
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
        physqrow=[physqb+m-1], physqcol=[physqvel,physqb+m-1], &
        oldvectors=oldvectors_ve, addmatvec=.true. )

    else

!     diagonal block
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
        physqrow=[physqb+m-1], physqcol=[physqb+m-1], &
        oldvectors=oldvectors_ve, addmatvec=.true. )

    end if

    if ( coefficients%i(61) == 1 ) then

!     direct velocity gradient

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=implicit_ce_vel_supg_elem, coefficients=coefficients, &
        physqrow=[physqb+m-1], physqcol=[physqvel], &
        oldvectors=oldvectors_ve, addmatvec=.true., buildvector=.false. )

      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        physqcol=[physqgrad], physqrow=[physqb+m-1], &
        buildvector=.false., zeromatvec=.true. )

    else

!     G- method

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=implicit_ce_vel_supg_elem, coefficients=coefficients, &
        physqrow=[physqb+m-1], physqcol=[physqgrad,physqvel], &
        oldvectors=oldvectors_ve, addmatvec=.true., buildvector=.false. )

    end if

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      physqrow=[physqgrad,physqpress], physqcol=[physqb+m-1], &
      buildvector=.false., zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      physqrow=[physqb+m-1], physqcol=[physqpress], &
      buildvector=.false., zeromatvec=.true. )

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problem, sysmatrix, &
      sysvector=rhsd, constraint1=2+m, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

  end subroutine build_b

  subroutine build_offdiag_b

    integer :: i, j

    if ( nmodes > 1 ) then
      do i = 1, nmodes
        do j = 1, nmodes
          if ( i == j ) cycle
          call build_system ( mesh, problem, sysmatrix, rhsd, &
            addmatvec=.true., physqrow=[physqb+i-1], physqcol=[physqb+j-1], &
            buildvector=.false., zeromatvec=.true. )
        end do
      end do
    end if

  end subroutine build_offdiag_b

  function func ( nr, x )
    use kind_defs_m
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func
    func = x(2)
  end function func

  subroutine steady_state ( shearrate, bsteady, csteady )

    real(dp), intent(in) :: shearrate
    real(dp), dimension(ncompb), intent(out) :: bsteady
    real(dp), dimension(3), intent(out) :: csteady

    integer, parameter :: &
             maxnumiterations = 30     ! maximum number of iterations

    real(dp), parameter :: eps = 1.e-12_dp ! Newton-Raphson accuracy

    real(dp) :: gradv(1,4)
    real(dp) :: b(1,ncompb,1), rhs(1,ncompb,1), c(1,3)
    real(dp) :: Hb(1,ncompb,1), dHb(1,ncompb,ncompb,1)
    real(dp) :: smat(ncompb,ncompb), fvec(ncompb)

    type(vemodel_t) :: vemodel
    type(vemopt_t) :: vemopt
    type(tdpar_t) :: tdpar

    integer :: info, ipiv(ncompb)

!   define the model

    call create_viscoelastic_model ( model, vemodel, nmodes=1, flowtype=0 )

    vemodel%modulus = G
    vemodel%lambda = lambda

    call create_vemopt ( vemopt, compute_drhs=.true., ncomp=ncompb, nmodes=1 )

!   set numerical parameters for time discretization

    tdpar%timestep = deltat
!   method 1: BDF1, 2: BDF2
    select case (timeint2)
    case(8)
      tdpar%method = 1
    case(10)
      tdpar%method = 2
    case default
      write(*,'(/a,i0/)') 'Error: wrong value timeint2: ', timeint2
      stop
    end select

!   velocity gradient

    gradv(1,:) = [ 0._dp, shearrate, 0._dp, 0._dp ]

!   initialize b
    b(1,:,1) = [ 2*(lambda*shearrate)**2, lambda*shearrate, 0._dp, 1._dp ]
    call sqrtc_2D_b( b(:,:,1) )

!   Newton-Raphson

    iter = 0

    do

      iter = iter + 1

      if ( iter > maxnumiterations ) then
        write(*,'(a,i0/)') &
          ' Maximum number of iterations reached = ', &
          maxnumiterations
        stop
      end if

!     right-hand side

      call rhs_viscoelastic_2D_b ( vemodel, gradv, b, rhs, vemopt=vemopt )

      call NRtd_steady_viscoelastic_2D_b ( vemodel, tdpar, b, Hb, dHb )

!     solve J^-1 * residual

      smat = - vemopt%drhs(1,:,:,1) + dHb(1,:,:,1)
      fvec = - rhs(1,:,1) + Hb(1,:,1)

      call dgesv( ncompb, 1, smat, ncompb, ipiv, fvec, ncompb, info )

      if ( info /= 0 ) then
        print *, 'Lapack DGESV: info = ', info
        stop
      end if

!     update for Newton-Raphson: un+1=un-J(un)^-1*residual(un)

      b(1,:,1) = b(1,:,1) - fvec

      if (printscreen) print *, iter, maxval(abs(fvec))

      if ( maxval(abs(fvec)) < eps ) exit

    end do

    bsteady = b(1,:,1)

    call conformation_2D_b ( b(:,:,1), c )

    csteady = c(1,:)

  end subroutine steady_state

end program couette27
