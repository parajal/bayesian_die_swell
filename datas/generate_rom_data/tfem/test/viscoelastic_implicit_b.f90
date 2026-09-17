! Startup of flow around a cylinder for a Giesekus model using the b-tensor
! formulation. Explicit div tau term.
! Fully-implicit with Newton-Raphson for the b-tensor.
! Problem 1: Flow around a confined cylinder.
!   Periodical boundary conditions.
!   Stress-explicit formulation of the momentum balance.
! Problem 2: b tensor

program cylinder3

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma57_m
  use hsl_ma41_m
  use io_utils_m
  use timer_m

  implicit none

! constants

!  logical, parameter :: &
!    printscreen = .true.  ! print on standard output

  integer, parameter :: &
    uintpl = 6,         & ! P2 velocities
    pintpl = 2,         & ! P1 pressures
    gintpl = 2,         & ! P1 gradients
    bintpl = 2,         & ! P1 b-tensor
    physqgrad = 1,      & ! physical quantity nr of the gradients
    physqvel = 2,       & ! physical quantity nr of the velocities
    physqpress = 3,     & ! physical quantity nr of the pressures
    gauss = 6,          & ! 6 point integration of triangles
    gaussb = 3,         & ! 3 point integration of boundary elements
    timeint1 = 8,       & ! first-order time integration, first time step
    timeint2 = 9,      & ! second-order time integration after first time step
    numtimesteps = 2,  & ! number of time steps
    maxnumiterations = 20, & ! maximum number of Newton-Raphson iterations
    ncompb = 4,         & ! number of b-tensor comp
    nmodes = 1,         & ! number of modes
    startm = 501,       & ! start of material model data
    model = 3,          & ! Giesekus
    bvariant = 1          ! b-formulation:
                          ! 1: CDT 2: symmetric
                          ! 3: Cholesky 4: Cholesky with log

  real(dp), parameter :: &
    G = 1.0_dp,          & ! modulus
    lambda = 1.0_dp,     & ! relaxation time
    eta_s = 0.59_dp/0.41_dp*lambda*G, & ! solvent viscosity
    mobility = 0.01_dp,  & ! mobility parameter
    H = 2._dp,           & ! (half-) height of the channel (must match mesh)
    U = 0.5_dp,          & ! average velocity in the channel
    flowrate = H*U         ! flow rate in (half) the channel

  real(dp), parameter :: &
    deltat = 0.2_dp,   & ! time step
    critval_reinit = 0.0_dp,  & ! criterion for reinitializing b
    thetapar = 0.55_dp, & ! theta parameter in the theta method
    epsconf = 1e-12_dp, & ! Newton-Raphson convergence threshold
    beta = 1,        & ! SUPG factor
    rs_gup = 1.2_dp, & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp, & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_b  = 1.0_dp,  & ! real_storage for the b tensor LU (HSL)
    is_b  = 1.6_dp     ! integer_storage for b tensor LU (HSL)

! definitions

  type(mesh_t) :: mesh, mesh1
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, solm1, solhat
  type(sysvector_t) :: rhsd
!  type(vector_t) :: velocity, pressure, vorticity
  type(subscript_t) :: velx, vely
!  type(plot_options_t) :: plot_options
!  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma57_t) :: solver_options_u
  type(lu_ma57_t) :: lu_u

  type(input_probdef_t) :: input_probdefb
  type(problem_t), target :: problemb
  type(sysmatrix_t) :: sysmatrixb
  type(oldvectors_t) :: oldvectors_ve
  type(subscript_t) :: bvalcmp(ncompb), bval
  type(subscriptvec_t) :: cxx, cxy, cyy
  type(sysvector_t), dimension(nmodes), target :: solbn, solbnm1, solbiter
  type(sysvector_t) :: dsolb, rhsb
  type(vector_t) :: ctensor, btensor
  type(solver_options_ma41_t) :: solver_options_b

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
  coefficients%i(84) = 2  ! storage of b tensor: all components, separate modes

  coefficients%r = &
    [ eta_s,    0._dp,   0._dp,   alpha, 0._dp, &
       flowrate, 0._dp,  deltat,   beta,  0._dp, &
       ( 0._dp, i = 11, 500 ), &
       G,     lambda, mobility &
    ]

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

! plot mesh

!  plot_options%fontsize = 6
!  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
!
!  plot_options%fontsize = 10
!
!  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=3 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 4  ! G
  input_probdef%vec_elementdof(1)%a(:,2) = 2         ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2,3]
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

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for the velocity for post processing

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd, solhat )

! fill solution vector with essential boundary conditions

  sol%u = 0

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix/vector for gradient/velocity/pressure problem

  call build_vpG

! solve initial gradient/velocity/pressure consistent with initial conditions
! and keep decomposition of the matrix

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options_u%real_storage=rs_gup
  solver_options_u%integer_storage=is_gup

  call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
    solver_options=solver_options_u  )


! problem definition b tensor

  call create_input_probdef ( mesh, input_probdefb, nvec=4, nphysq=1 )

  input_probdefb%vec_elementdof(1)%a(vertices,1) = ncompb ! b
  input_probdefb%vec_elementdof(1)%a(:,2) = 1 ! scalar for plotting
  input_probdefb%vec_elementdof(1)%a(:,3) = 3 ! tensor for plotting
  input_probdefb%vec_elementdof(1)%a(:,4) = 4 ! tensor for plotting

  input_probdefb%physq = [1]
  input_probdefb%probnr = 2

! constraint for periodical boundary conditions of the b tensor
  call define_constraint ( mesh, input_probdefb, curve1=4, curve2=12, &
    discretization='collocation' )

  call problem_definition ( input_probdefb, mesh, problemb )


! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemb, bval, physqarr=[1] )
  do i = 1, ncompb
    call create_subscript ( mesh, problemb, bvalcmp(i), physqarr=[1], degfd=i )
  end do
  call create_subscript ( mesh, problemb, cxx, degfd=1, vec=3 )
  call create_subscript ( mesh, problemb, cxy, degfd=2, vec=3 )
  call create_subscript ( mesh, problemb, cyy, degfd=3, vec=3 )

! create a vector for conformation tensor for post processing

  call create_vector ( problemb, ctensor, vec=3 )
  call create_vector ( problemb, btensor, vec=4 )

! create system vectors (solution and right-hand side) for b tensor and
! initialize vectors.

  call create ( problemb, solbn, solbnm1, solbiter )
  call create ( problemb, dsolb, rhsb )

! initial solution

  do m = 1, nmodes
    if ( bvariant == 4 ) then
!     log Cholesky
      solbn(m)%u = 0
    else
      solbn(m)%u(bvalcmp(1)%s) = 1 ! initial bxx
      solbn(m)%u(bvalcmp(2)%s) = 0 ! initial bxy
      solbn(m)%u(bvalcmp(3)%s) = 0 ! initial bxy
      solbn(m)%u(bvalcmp(4)%s) = 1 ! initial byy
    end if
  end do

  call copy(solbn,solbiter)


! create system matrix for b tensor problem

  call create_sysmatrix_structure_base ( sysmatrixb, mesh, problemb )
  call create_sysmatrix_structure_constraint ( sysmatrixb, mesh, problemb )
  call finalize_sysmatrix_structure ( sysmatrixb )

  call create_sysmatrix_data ( sysmatrixb )

! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec1=3, nprob=2 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => solhat
  oldvectors_ve%s(2)%p => solm1
  oldvectors_ve%s1(1)%p => solbn
  oldvectors_ve%s1(2)%p => solbnm1
  oldvectors_ve%s1(3)%p => solbiter
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemb


  !open ( unit=13, file='iter.out', recl=300 )

! time stepping

  call tic

  do step = 1, numtimesteps

    write(*,*) 'step = ', step
!    if ( printscreen) print *, 'step = ', step

!   prediction of velocity and gradient for new time step
    if ( step == 1 ) then
 !    first order prediction
      solhat%u = sol%u
    else if ( step >= 2 .and. any ( timeint2 == [9,10,11] ) ) then
!     second order prediction
      solhat%u = 2*sol%u - solm1%u
!      solhat%u = sol%u
    end if

    if ( step == 2 ) then
!     change time integration scheme at the second time step
      coefficients%i(22) = timeint2
    end if

    solver_options_b%real_storage=rs_b
    solver_options_b%integer_storage=is_b

    do m = 1, nmodes

      coefficients%i(85) = m   ! set mode number

      if ( nmodes > 1 ) then
        write(*,*) 'mode = ', m
!        if ( printscreen) print *, 'mode = ', m
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

!       build (assemble) matrix and vector for b-tensor problem

        call build_system ( mesh, problemb, sysmatrixb, sysvector=rhsb, &
          elemsub=implicit_ce_supg_elem, oldvectors=oldvectors_ve, &
          coefficients=coefficients )

!       periodical condition on b-tensor tensor
        call build_system_constraint ( mesh, problemb, sysmatrixb, &
          sysvector=rhsb, elemsub=stokes_constr_node_conn, &
          addmatvec=.true. )

        call check ( sysmatrixb )

!       solve b tensor

        call solve_system_ma41 ( sysmatrixb, rhsb, dsolb, &
          solver_options=solver_options_b  )

        solbiter(m)%u(bval%s) = solbiter(m)%u(bval%s) + dsolb%u(bval%s)

!        if ( printscreen) print *, maxval(abs(dsolb%u(bval%s)))
        write(*,*) maxval(abs(dsolb%u(bval%s)))

        if ( maxval(abs(dsolb%u(bval%s))) < epsconf ) exit

      end do

    end do

!   copy previous b-tensor solution to older time step

    call copy ( solbn, solbnm1 )
    call copy ( solbiter, solbn )

!   reinitialize b to b' = sqrt(c)

    if ( bvariant == 1 ) then
      call reinitialize_b
!      call copy ( solbn, solbiter)
    end if

!   build (assemble) vector for gradient/velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_divtau, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      buildmatrix=.false., coefficients=coefficients )

!   imposed flow rate (right-hand side only)

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      buildmatrix=.false., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   copy previous gradient/velocity/pressure solution to older time step

    call copy ( sol, solm1 )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
      solver_options=solver_options_u  )


!   write max and mean values of conformation tensor to a file

    call derive_vector ( mesh, problemb, ctensor, &
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
                                    maxval(sol%u(velx%s)), &
                                    maxval(sol%u(vely%s))
!    close(unit=11)

    call toc ( 'one step' )

  end do

  call toc ( 'all steps' )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, solm1, rhsd, solhat )
  call delete ( sysmatrix )
  call delete ( coefficients )
!  call delete ( oldvectors )
  call delete ( lu_u )

  call delete ( problemb )
  call delete ( input_probdefb )
  call delete ( solbn, solbnm1, solbiter )
  call delete ( dsolb, rhsb )
  call delete ( sysmatrixb )
  call delete ( oldvectors_ve )
  call delete ( bval )
  do i = 1, ncompb
    call delete ( bvalcmp(i) )
  end do

contains

  subroutine build_vpG

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients, &
      physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress] )

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


! reinitialize b to b = sqrt(c)

  subroutine reinitialize_b

    real(dp) :: b(size(bvalcmp(1)%s),ncompb), bm1(size(bvalcmp(1)%s),ncompb)
    real(dp) :: RT(size(bvalcmp(1)%s),2,2)
    integer :: i, m

    do m = 1, nmodes

!     obtain the solution of b at current time step

      do i = 1,ncompb
        b(:,i) = solbn(m)%u(bvalcmp(i)%s)
      end do

!     check skew norm and perform reinitialization if exceeded

      if ( any ( skew_norm_2D_b(b) >= critval_reinit ) ) then

!       determine the square root of b*b^T

        call sqrtc_2D_b ( b, RT=RT )

!       obtain the solution bn-1 at previous time step

        do i = 1,ncompb
          bm1(:,i) = solbnm1(m)%u(bvalcmp(i)%s)
        end do

!       rotate bn-1 according to b

        call rotate_2D_b ( bm1, RT )

!       store b and bn-1 in solution vectors

        do i = 1,ncompb
          solbn(m)%u(bvalcmp(i)%s) = b(:,i)
          solbnm1(m)%u(bvalcmp(i)%s) = bm1(:,i)
        end do

      end if

    end do

  end subroutine reinitialize_b

end program cylinder3
