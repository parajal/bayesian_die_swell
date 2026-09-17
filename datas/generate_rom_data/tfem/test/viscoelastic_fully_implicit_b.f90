! Startup of flow around a cylinder for a b-representation differential
! Giesekus model.
! Fully-implicit with Newton-Raphson
! Periodical boundary conditions.

program cylinder5

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use io_utils_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,     & ! P2 velocities
    pintpl = 2,     & ! P1 pressures
    gintpl = 2,     & ! P1 gradients
    bintpl = 2,     & ! P1 b-tensor
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    physqb = 4,     & ! physical quantity nr of the b-tensor (mode 1)
    gauss = 6,      & ! 6 point integration of triangles
    gaussb = 3,     & ! 3 point integration of boundary elements
    timeint1 = 8,   & ! first-order time integration, first time step
    timeint2 = 10,  & ! second-order time integration after first time step
    numtimesteps = 2, & ! number of time steps
    maxnumiterations = 20, & ! maximum number of Newton-Raphson iterations
    ncompb = 4,     & ! number of b-tensor components
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    model = 3,      & ! Giesekus
    bvariant = 1      ! b-formulation:
                      ! 1: CDT 2: symmetric
                      ! 3: Cholesky 4: Cholesky with log

  real(dp), parameter :: &
    G     = 1.0_dp,    & ! modulus
    lambda = 1.0_dp,   & ! relaxation time
    eta_s = 0.59_dp/0.41_dp*lambda*G,    & ! solvent viscosity
    mobility = 0.01_dp,   & ! mobility parameter
    H = 2._dp,         & ! (half-) height of the channel (must match mesh)
    U = 0.5_dp,        & ! average velocity in the channel
    flowrate = H*U       ! flow rate in (half) the channel

  real(dp), parameter :: &
    deltat = 0.4_dp,     & ! time step
    critval_reinit = 0.0_dp,  & ! criterion for reinitializing b
    thetapar = 0.55_dp,  & ! theta parameter in the theta method
    epsconf = 1e-12_dp,  & ! Newton-Raphson convergence threshold
    beta = 1,            & ! SUPG factor
    rs = 1.2_dp,         & ! real_storage for LU (HSL)
    is = 1.6_dp            ! integer_storage for LU (HSL)

!  logical :: &
!    printscreen = .true.  ! print on standard output


! definitions

  type(mesh_t) :: mesh, mesh1
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: soln, solnm1, soliter
  type(sysvector_t) :: dsol, rhsd
!  type(vector_t) :: velocity, pressure, vorticity
  type(subscript_t) :: velx, vely, vel, pres, solsc
  type(subscript_t), dimension(nmodes) :: bval
  type(subscript_t), dimension(ncompb,nmodes) :: bvalcmp
!  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options

  type(subscriptvec_t) :: cxx, cxy, cyy
  type(vector_t) :: ctensor, btensor

  integer :: step, i, iter, m
  integer :: vertices(3) = [1,3,5]
  real(dp) :: alpha

  timer = .false.

! set some parameters

  alpha = G * lambda   ! DEVSS parameter

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+3*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
       physqvel, physqpress, 0,     physqgrad, gauss,  &
       gaussb,   bintpl,     0,     0,         0,      &
       0,        0,          model, nmodes,    startm, &
       0,     timeint1,   ( 0, i = 23, 150 )  &
    ]

  coefficients%i(71) = bvariant  ! variant for b-formulation
  coefficients%i(83) = physqb ! physical quantity if first mode
  coefficients%i(84) = 3  ! storage of b tensor: all modes in sysvector
  coefficients%i(86) = 2  ! sysvector number for velocity in SUPG
  coefficients%i(87) = 1  ! sysvector number for iteration b
  coefficients%i(88) = 2  ! sysvector number of b at tn

  coefficients%r = &
    [ eta_s,    0._dp,   0._dp,   alpha, 0._dp, &
       flowrate, 0._dp,  deltat,   beta,  0._dp, &
       ( 0._dp, i = 11, 500 ), &
       G,     lambda, mobility &
    ]

  coefficients%r(10) = U ! global scaling velocity for SUPG, when needed
  coefficients%r(28) = thetapar

! read mesh

  call read_mesh_gmsh ( mesh1, filename='confined_cylinder1.msh', ndim=2 )
  call mesh_convert ( mesh1, mesh, remove_isolated_nodes=.true. )
  call delete ( mesh1 )

! cylinder
  call add_to_mesh ( mesh, curve=[-2] ) ! curve 9
! top
  call add_to_mesh ( mesh, curve=[5,6,7] ) ! curve 10
! centerline+cylinder
  call add_to_mesh ( mesh, curve=[1,-2,3] ) ! curve 11
! connect curves for periodical bc in DG
  call add_to_mesh ( mesh, curve=[-8] )     ! curve 12

  call fill_mesh_parts ( mesh )

!  call printinfo ( mesh, printlevel=1 )

!! plot mesh
!
!  plot_options%fontsize = 6
!  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
!
!  plot_options%fontsize = 10
!
!  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=7, nphysq=4 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 4  ! G
  input_probdef%vec_elementdof(1)%a(:,2) = 2         ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(vertices,4) = ncompb  ! b
  input_probdef%vec_elementdof(1)%a(:,5) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,6) = 3  ! tensor for plotting
  input_probdef%vec_elementdof(1)%a(:,7) = 4  ! tensor for plotting

  input_probdef%physq = [1,2,3,(3+m,m=1,nmodes)]
  input_probdef%probnr = 1

! define essential boundaries

! cylinder
  call define_essential ( mesh, input_probdef, curve1=9, physq=physqvel )
! top (wall)
  call define_essential ( mesh, input_probdef, curve1=10, physq=physqvel )
! bottom (center line)
  call define_essential ( mesh, input_probdef, curve1=1, &
    physq=physqvel, degfd=[0,1] )
  call define_essential ( mesh, input_probdef, curve1=3, &
    physq=physqvel, degfd=[0,1] )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=12, nglobalc=1 )

! constraints for periodical boundary conditions

! velocities (use weak connection)
  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=4, curve2=12, discretization='weak', &
    elementdof=[2,0,2] )

! gradients (use collocation)
  call define_constraint ( mesh, input_probdef, &
    physq=physqgrad, curve1=4, curve2=12, discretization='collocation' )

! conformation tensor (use collocation)
  do m = 1, nmodes
    call define_constraint ( mesh, input_probdef, &
      physq=physqb+m-1, curve1=4, curve2=12, discretization='collocation' )
  end do

  call problem_definition ( input_probdef, mesh, problem )

 ! create vector subscripts for solution (excluding contraint forces)

  call create_subscript ( mesh, problem, solsc, &
    physqarr=[physqgrad,physqvel,physqpress,(physqb+m-1,m=1,nmodes)] )

! create vector subscripts for the velocity

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )

! create vector subscript for the pressure

  call create_subscript ( mesh, problem, pres, physqarr=[physqpress] )

! create vector subscript for the conformation tensor

  do m = 1, nmodes
    call create_subscript ( mesh, problem, bval(m), physqarr=[physqb+m-1] )
    do i = 1, ncompb
      call create_subscript ( mesh, problem, bvalcmp(i,m), &
        physqarr=[physqb+m-1], degfd=i )
    end do
  end do
  call create_subscript ( mesh, problem, cxx, degfd=1, vec=6 )
  call create_subscript ( mesh, problem, cxy, degfd=2, vec=6 )
  call create_subscript ( mesh, problem, cyy, degfd=3, vec=6 )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, soln, solnm1, soliter )
  call create_sysvector ( problem, dsol, rhsd )


! create a vector for conformation
! tensor for post processing

  call create_vector ( problem, ctensor, vec=6 )
  call create_vector ( problem, btensor, vec=7 )

! fill solution vector with essential boundary conditions

  soln%u = 0

! initial solution b tensor

  do m = 1, nmodes
    soln%u(bvalcmp(1,m)%s) = 1 ! initial bxx
    soln%u(bvalcmp(2,m)%s) = 0 ! initial bxy
    soln%u(bvalcmp(3,m)%s) = 0 ! initial byx
    soln%u(bvalcmp(4,m)%s) = 1 ! initial byy
  end do

  call copy(soln,solnm1)
  call copy(soln,soliter)
  call copy(soln,dsol)


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


!  open ( unit=13, file='iter.out', recl=300 )

  solver_options%real_storage=rs
  solver_options%integer_storage=is

  call tic

! time stepping

  do step = 1, numtimesteps

    write(*,*) 'step = ', step
!    if ( printscreen) print *, 'step = ', step

    !if ( step == 5 ) coefficients%r(8) = 1.0_dp

    if ( step == 1 ) then
      coefficients%i(78) = 1 ! check for zero velocity in SUPG
    else
      coefficients%i(78) = 0
    end if

    if ( step >= 2 ) then
      coefficients%i(22) = timeint2
    end if

    iter = 0

    do

      iter = iter + 1

      if ( step == 1 .and. iter == 2 ) then
        coefficients%r(6) = 0  ! set flowrate in dsol=0
      end if

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

!      if ( printscreen) print *, maxval(abs(dsol%u(solsc%s))), &
!        maxval(abs(dsol%u(vel%s))), maxval(abs(dsol%u(pres%s))), &
!        maxval(abs(dsol%u(bval(1)%s)))
      write(*,*) maxval(abs(dsol%u(solsc%s))), &
        maxval(abs(dsol%u(vel%s))), maxval(abs(dsol%u(pres%s))), &
        maxval(abs(dsol%u(bval(1)%s)))

      if ( maxval(abs(dsol%u(solsc%s))) < epsconf ) exit

    end do

!   copy solution to older time step for next time step

    call copy ( soln, solnm1 )
    call copy ( soliter, soln )

!   reinitialize b to b' = sqrt(c)

    if ( bvariant == 1 ) then
      call reinitialize_b
    end if

!   write max and mean values of conformation and b-tensor tensor to a file

    call derive_vector ( mesh, problem, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

!    if ( step == 1 ) then
!      open(unit=11, recl=300, status='replace', file='cval.out')
!    else
!      open(unit=11, recl=300, position='append', file='cval.out')
!    end if
    write(*, fmt=*) step * deltat, maxval(ctensor%u(cxx%s)), &
                                    maxval(ctensor%u(cxy%s)), &
                                    maxval(ctensor%u(cyy%s)), &
                                    sum(ctensor%u(cxx%s))/size(cxx%s), &
                                    sum(ctensor%u(cxy%s))/size(cxy%s), &
                                    sum(ctensor%u(cyy%s))/size(cyy%s), &
                                    maxval(soln%u(velx%s)), &
                                    maxval(soln%u(vely%s))
 !   close(unit=11)

    call toc ( 'one step' )

  end do

  call toc ( 'all steps' )



! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( soln, solnm1, soliter, dsol, rhsd )
!  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
!  call delete ( oldvectors )

  call delete ( oldvectors_ve )
  do m = 1, nmodes
    do i = 1, ncompb
      call delete ( bvalcmp(i,m), bval(m) )
    end do
  end do
  call delete ( cxx, cxy, cyy )

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

!   flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients )

!   periodical condition on velocities

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=stokes_constr_elem_conn, addmatvec=.true., &
      coefficients=coefficients )

!   periodical condition on gradients

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=3, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients )

  end subroutine build_vpG


  subroutine build_b

!   build (assemble) matrix and vector for b-tensor mode m

!   diagonal block
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_ce_supg_elem, coefficients=coefficients, &
      physqrow=[physqb+m-1], physqcol=[physqb+m-1], &
      oldvectors=oldvectors_ve, addmatvec=.true. )

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=implicit_ce_vel_supg_elem, coefficients=coefficients, &
      physqrow=[physqb+m-1], physqcol=[physqgrad,physqvel], &
      oldvectors=oldvectors_ve, addmatvec=.true., buildvector=.false. )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      physqrow=[physqgrad,physqpress], physqcol=[physqb+m-1], &
      buildvector=.false., zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      physqrow=[physqb+m-1], physqcol=[physqpress], &
      buildvector=.false., zeromatvec=.true. )

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problem, sysmatrix, &
      sysvector=rhsd, constraint1=3+m, elemsub=stokes_constr_node_conn, &
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


! reinitialize b to b = sqrt(c)

  subroutine reinitialize_b

    real(dp) :: b(size(bvalcmp(1,1)%s),ncompb), bm1(size(bvalcmp(1,1)%s),ncompb)
    real(dp) :: RT(size(bvalcmp(1,1)%s),2,2)
    integer :: i, m

    do m = 1, nmodes

!     obtain the solution of b at current time step

      do i = 1,ncompb
        b(:,i) = soln%u(bvalcmp(i,m)%s)
      end do

!     check skew norm and perform reinitialization if exceeded

      if ( any ( skew_norm_2D_b(b) >= critval_reinit ) ) then

!       determine the square root of b*b^T

        call sqrtc_2D_b ( b, RT=RT )

!       obtain the solution bn-1 at previous time step

        do i = 1,ncompb
          bm1(:,i) = solnm1%u(bvalcmp(i,m)%s)
        end do

!       rotate bn-1 according to b

        call rotate_2D_b ( bm1, RT )

!       store b and bn-1 in solution vectors

        do i = 1,ncompb
          soln%u(bvalcmp(i,m)%s) = b(:,i)
          solnm1%u(bvalcmp(i,m)%s) = bm1(:,i)
        end do

      end if

    end do

  end subroutine reinitialize_b

end program cylinder5
