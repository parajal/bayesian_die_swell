! Startup of flow around a cylinder for a b-representation differential
! Giesekus model. Semi-implicit div tau term.
! Fully-implicit with Newton-Raphson for the b-tensor.
! Problem 1: Flow around a confined cylinder.
!   Periodical boundary conditions.
!   Stress-implicit formulation of the momentum balance.
! Problem 2: b tensor
! Problem 3: b*b^T projection

program cylinder4

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma57_m
  use hsl_ma41_m
  use io_utils_m
  use figplot_m
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
    gauss = 6,      & ! 6 point integration of triangles
    gaussb = 3,     & ! 3 point integration of boundary elements
    timeint1 = 8,   & ! first-order time integration, first time step
    timeint2 = 10,  & ! second-order time integration after first time step
    numtimesteps = 30, & ! number of time steps
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
    eta_p = lambda*G,  & ! polymer viscosity
    beta_s = 0.59_dp,  & ! beta viscosity parameter = eta_s/(eta_s+eta_p)
    eta_s = beta_s/(1-beta_s)*eta_p,    & ! solvent viscosity
!    eta_0 = eta_s+eta_p, & ! zero-shear viscosity
    mobility = 0.01_dp,   & ! mobility parameter
    H = 2._dp,         & ! (half-) height of the channel (must match mesh)
    R = 1._dp,         & ! radius of the cylinder (must match mesh)
    U = 0.5_dp,        & ! average velocity in the channel
    flowrate = H*U       ! flow rate in (half) the channel

  real(dp), parameter :: &
    deltat = 0.6_dp,     & ! time step
    epsconf = 1e-12_dp,  & ! Newton-Raphson convergence threshold
    beta_supg = 1,       & ! SUPG factor
    rs_gup = 1.2_dp,     & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,     & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_b  = 1.0_dp,      & ! real_storage for the b-tensor LU (HSL)
    is_b  = 1.6_dp         ! integer_storage for b-tensor LU (HSL)

  logical :: &
    cproj = .true.,      & ! projection of c=b*b^T for cn in momentum balance
    printscreen = .true.  ! print on standard output


! definitions

  type(mesh_t) :: mesh, mesh1
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, solm1
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
  type(subscript_t) :: velx, vely
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options_u

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

  type(input_probdef_t) :: input_probdef_proj
  type(problem_t), target :: problem_proj
  type(sysmatrix_t) :: sysmatrix_proj
  type(sysvector_t), dimension(ncompb-1,nmodes), target :: solc_proj
  type(sysvector_t), dimension(ncompb-1,nmodes) :: rhsc_proj
  type(lu_ma57_t) :: lu_proj

  integer :: step, i, iter, m
  integer :: vertices(3) = [1,3,5]
  real(dp) :: alpha, Wi

  timer = .false.

! set some parameters

  alpha = eta_p   ! DEVSS parameter

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
  if ( cproj ) coefficients%i(72) = 1  ! use c projection in momentum balance
  coefficients%i(84) = 2  ! storage of b tensor: all components, separate modes
  coefficients%i(90) = 1  ! rotation reinitialization on element level

  coefficients%r = &
    [ eta_s,    0._dp, 0._dp,  alpha,     0._dp, &
      flowrate, 0._dp, deltat, beta_supg, 0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G, lambda, mobility &
    ]

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

  call printinfo ( mesh, printlevel=1 )

! plot mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  plot_options%fontsize = 10

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


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

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0
  call copy(sol,solm1)


! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! problem definition b-tensor

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


! create a vector subscript for the conformation and b-tensor
! tensor without Lagr. multipl. for post processing the data

  call create_subscript ( mesh, problemb, bval, physqarr=[1] )
  do i = 1, ncompb
    call create_subscript ( mesh, problemb, bvalcmp(i), physqarr=[1], degfd=i )
  end do
  call create_subscript ( mesh, problemb, cxx, degfd=1, vec=3 )
  call create_subscript ( mesh, problemb, cxy, degfd=2, vec=3 )
  call create_subscript ( mesh, problemb, cyy, degfd=3, vec=3 )

! create a vector for conformation and b-tensor
! tensor for post processing

  call create_vector ( problemb, ctensor, vec=3 )
  call create_vector ( problemb, btensor, vec=4 )

! create system vectors (solution and right-hand side) for b-tensor
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


! create system matrix for b-tensor problem

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
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemb


! problem definition for projected "c=b*b^T" of the b-tensor

  call create_input_probdef ( mesh, input_probdef_proj, nvec=1, nphysq=1 )

  input_probdef_proj%vec_elementdof(1)%a = &
      reshape ( [ 1,0,1,0,1,0 ], &
                   [6,1] )

  input_probdef_proj%physq = [1]
  input_probdef_proj%probnr = 3

  call problem_definition ( input_probdef_proj, mesh, problem_proj )

  call create ( problem_proj, solc_proj, rhsc_proj )

  do m = 1, nmodes
    solc_proj(1,m)%u = 1   ! initial cxx
    solc_proj(2,m)%u = 0   ! initial cxy
    solc_proj(3,m)%u = 1   ! initial cyy
  end do

! store solution vectors and problem structures

  oldvectors_ve%s2(3)%p => solc_proj
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


! Weisenberg number
  Wi = lambda * U / R

  print *, ' U = ', U, 'Wi = ', Wi
  print *


  open ( unit=13, file='iter.out', recl=300 )

! time stepping

  call tic

  do step = 1, numtimesteps

    write(13,*) 'step = ', step
    if ( printscreen) print *, 'step = ', step

    if ( step >= 2 ) then
      coefficients%i(22) = timeint2
    end if

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG

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


    solver_options_b%real_storage=rs_b
    solver_options_b%integer_storage=is_b

    do m = 1, nmodes

      coefficients%i(85) = m   ! set mode number

      if ( nmodes > 1 ) then
        write(13,*) 'mode = ', m
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

!       build (assemble) matrix and vector for b-tensor problem

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

        if ( printscreen) print *, iter, maxval(abs(dsolb%u(bval%s)))
        write(13,*) iter, maxval(abs(dsolb%u(bval%s)))

        if ( maxval(abs(dsolb%u(bval%s))) < epsconf ) exit

      end do

    end do

!   copy b-tensor solution to older time step for next time step

    call copy ( solbn, solbnm1 )
    call copy ( solbiter, solbn )

!   copy velocity solution to older time step for next time step

    call copy ( sol, solm1 )


!   write max and mean values of conformation and b-tensor tensor to a file

    call derive_vector ( mesh, problemb, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    if ( step == 1 ) then
      open(unit=11, recl=300, status='replace', file='cval.out')
    else
      open(unit=11, recl=300, position='append', file='cval.out')
    end if
    write(11, fmt=*) step * deltat, maxval(ctensor%u(cxx%s)), &
                                    maxval(ctensor%u(cxy%s)), &
                                    maxval(ctensor%u(cyy%s)), &
                                    sum(ctensor%u(cxx%s))/size(cxx%s), &
                                    sum(ctensor%u(cxy%s))/size(cxy%s), &
                                    sum(ctensor%u(cyy%s))/size(cyy%s), &
                                    maxval(sol%u(velx%s)), &
                                    maxval(sol%u(vely%s))
    close(unit=11)

    if ( printscreen ) &
         print *, 'step = ', step, 'max cxx = ', maxval(ctensor%u(cxx%s))

    call toc ( 'one step' )

  end do

  call toc ( 'all steps' )


! post-processing

  call create_vector ( problem, velocity, physq=2 )
  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, vorticity, vec=4 )

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

  call write_vector_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='velocity', vector=velocity )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='pressure', vector=pressure, append=.true. )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='vorticity', vector=vorticity, append=.true. )

  call printtofile ( mesh, problem, filename='velocity_cl.out', curve=11, &
    vector=velocity )
  call printtofile ( mesh, problem, filename='pressure_cl.out', curve=11, &
    vector=pressure )

! conformation tensor

  call derive_vector ( mesh, problemb, ctensor, &
    elemsub=deriv_conformation_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemb, filename='c.vtk', &
    dataname='conformation_tensor', vector=ctensor )

  call printtofile ( mesh, problemb, filename='conformation_cl.out', curve=11, &
    vector=ctensor )

  call delete ( ctensor )

! b tensor

  call derive_vector ( mesh, problemb, btensor, &
    elemsub=deriv_conformation_tensor_std, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemb, filename='b.vtk', &
    dataname='b-tensor', vector=btensor, symmetric=.false. )

  call printtofile ( mesh, problemb, filename='b-tensor_cl.out', &
    curve=11, vector=btensor )

  call delete ( btensor )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

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
  call delete ( cxx, cxy, cyy )

  call delete ( sysmatrix_proj )
  call delete ( problem_proj )
  call delete ( input_probdef_proj )
  call delete ( solc_proj, rhsc_proj )
  if ( cproj ) call delete ( lu_proj )

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

end program cylinder4
