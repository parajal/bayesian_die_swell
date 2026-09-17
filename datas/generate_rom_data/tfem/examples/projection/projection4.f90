! 2D viscoelastic problem: flow around a cylinder confined between two walls.
! The entry and exit sides are connected by periodical boundary conditions
! and a constant flow rate is imposed.
! Only the upper-half of the domain is modelled and symmetry conditions are
! assumed (zero tractions).
!
! The average velocity is U.
!
! DEVSS-G/SUPG
! first-order time integration to steady state
!
! Gmsh generated mesh.
!
! The final solution of the gradient and conformation is projected
! onto a new mesh.

program projection4

  use tfem_m
  use hsl_ma41_m
  use hsl_ma57_m
  use viscoelastic_elements_m
  use projection_elements_m
  use figplot_m
  use io_utils_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,     & ! P2 velocities
    pintpl = 2,     & ! P1 pressures
    gintpl = 2,     & ! P1 gradients
    cintpl = 2,     & ! P1 conformation
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    gauss = 6,      & ! 6 point integration of triangles
    gaussb = 3,     & ! 3 point integration of boundary elements
    gaussp = 6,     & ! 6 point integration of triangles for projection
    ncompc = 3,     & ! number of conformation tensor components
    nmodes = 1,     & ! number of modes
    ncompv = 7,     & ! number of components input vector for projection
    startm = 501      ! start of material model data

! definitions

  type(mesh_t), target :: mesh
  type(mesh_t) :: mesh1
  type(input_probdef_t) :: input_probdef, input_probdefc
  type(problem_t), target :: problem, problemc
  type(sysmatrix_t) :: sysmatrix, sysmatrixc
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: velx, vely, cval
  type(subscript_t) :: gxxval, gxyval, gyxval, gyyval
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc
  type(vector_t) :: conformation_tensor, gradient_tensor
  !type(vector_t) :: Gxx
  type(vector_t), target :: vec_p
  type(lu_ma57_t) :: lu_u
  type(lu_ma41_t) :: luc
  type(solver_options_ma57_t) :: solver_options_u
  type(solver_options_ma41_t) :: solver_options_c
  type(plot_options_t) :: plot_options

  type(mesh_t) :: mesh2
  type(input_probdef_t) :: input_probdef2
  type(problem_t) :: problem2
  type(sysmatrix_t) :: sysmatrix2
  type(sysvector_t), dimension(ncompv), target :: sol2
  type(sysvector_t), dimension(ncompv) :: rhsd2
  type(oldvectors_t) :: oldvectors2
  type(coefficients_t) :: coefficients2
  type(lu_ma57_t) :: lu2

! variables

  integer :: &
    timeint = 1,         & ! (first-order) time integration
    numtimesteps = 2,    & ! number of time steps
    logc = 0,            & ! standard scheme or log transformation
    model = 2              ! Oldroyd-B

  real(dp) :: &
    ly = 2._dp,        & ! size in y-direction
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 1.0_dp,   & ! relaxation time
    alphapar = 0.1_dp, & ! alpha parameter in the Giesekus model
    deltat = 5.e-3_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    U = 1._dp,         & ! average velocity in the channel
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.0_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 1.4_dp       ! integer_storage for conformation LU (HSL)


  integer :: icomp, step, i, npar, nb
  integer :: vertices(3) = [1,3,5]
  real(dp) :: alpha, G, flowrate


! namelist for input of variables; read from standard input

  namelist /comppar/ ly, timeint, numtimesteps, logc, &
    model, eta_s, eta_p, lambda, alphapar, deltat, beta, U, &
    rs_gup, is_gup, rs_c, is_c

  read ( unit=*, nml=comppar )

! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus

  flowrate = U * ly

  if ( model == 2 ) then
    npar = 2 ! Oldroyd-B
  else if ( model == 3 ) then
    npar = 3 ! Giesekus
  end if


  !call tic

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+npar*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
       physqvel, physqpress, 0,     physqgrad, gauss,  &
       gaussb,   cintpl,     0,     0,         0,      &
       0,        0,          model, nmodes,    startm, &
       logc,     timeint,    ( 0, i = 23, 150 )  &
    ]

  coefficients%r(1:502) = &
    [ eta_s,    0._dp,   0._dp,   alpha, 0._dp, &
      flowrate, 0._dp,  deltat,   beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G,     lambda  &
    ]

  if ( model == 3 ) then
    coefficients%r(503) = alphapar ! Giesekus
  end if

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

  nb = nint( real(mesh%nelem) ** 0.25 ); print *, 'nb=', nb

  call add_to_mesh ( mesh, blocks=[nb,nb] )

  call fill_mesh_parts ( mesh )

  call printinfo ( mesh, printlevel=4 )

! plot mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  plot_options%fontsize = 10

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=6, nphysq=3 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 4  ! gradients
  input_probdef%vec_elementdof(1)%a(:,2) = 2  ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5) = 3  ! symmetric tensor
  input_probdef%vec_elementdof(1)%a(:,6) = 4  ! unsymmetric tensor

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


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=4, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a(:,1) = 0
  input_probdefc%vec_elementdof(1)%a(vertices,1) = 1  ! c
  input_probdefc%vec_elementdof(1)%a(:,2) = 1         ! scalar for plotting
  input_probdefc%vec_elementdof(1)%a(:,3) = 3         ! symmetric tensor
  input_probdefc%vec_elementdof(1)%a(:,4) = 0              ! transfer to
  input_probdefc%vec_elementdof(1)%a(vertices,4) = ncompv  ! projection module


  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

! constraint for periodical boundary conditions of the conformation
  call define_constraint ( mesh, input_probdefc, curve1=4, curve2=12, &
    discretization='collocation' )

  call problem_definition ( input_probdefc, mesh, problemc )

! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemc, cval, physqarr=[1] )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )


! fill solution vector with essential boundary conditions

  sol%u = 0

! create system matrix of gradient/velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector for gradient/velocity/pressure problem

! stokes velocity/pressure
  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, physqrow=[2,3], physqcol=[2,3], &
    coefficients=coefficients )

! DEVSS-G
  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=devssg_elem, addmatvec=.true., &
    physqrow=[1,2], physqcol=[1,2], coefficients=coefficients )

! set to zero off-diagonal blocks gradient-pressure
  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., physqrow=[1], physqcol=[3], zeromatvec=.true. )
  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., physqrow=[3], physqcol=[1], zeromatvec=.true. )

! imposed flow rate

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
    coefficients=coefficients )

! periodical condition on velocities

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_elem_conn, addmatvec=.true., &
    coefficients=coefficients )

! periodical condition on gradients

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=3, elemsub=stokes_constr_node_conn, addmatvec=.true., &
    coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

! solve initial gradient/velocity/pressure consistent with initial conditions
! and keep decomposition of the matrix
! This generates a Stokes profile.

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options_u%real_storage=rs_gup
  solver_options_u%integer_storage=is_gup

  call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
    solver_options=solver_options_u  )


! create system vectors (solution and right-hand side) for conformation and
! initialize vectors to zero stress

  call create ( problemc, solc, rhsc )

  if ( logc == 0 ) then ! standard
    solc(1,1)%u = 1
    solc(2,1)%u = 0
    solc(3,1)%u = 1
  else if ( logc == 1 ) then ! log scheme
    do i = 1, ncompc
      solc(i,1)%u = 0
    end do
  end if

! create system matrix for conformation problem

  call create_sysmatrix_structure_base ( sysmatrixc, mesh, problemc )
  call create_sysmatrix_structure_constraint ( sysmatrixc, mesh, problemc )
  call finalize_sysmatrix_structure ( sysmatrixc )

  call create_sysmatrix_data ( sysmatrixc )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=1, nsysvec2=1, nprob=2 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc


! open monitor file

  open ( unit=13, file='out', recl=300 )

  !call toc ( 'just before stepping' )

! time stepping

  do step = 1, numtimesteps


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


!   build (assemble) vector for gradient/velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_divtau, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      buildmatrix = .false., coefficients=coefficients )

!   imposed flow rate (right-hand side only)

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      buildmatrix = .false., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
      solver_options=solver_options_u  )


!   write monitor data

    write(unit=13,fmt=*) &
      step, step*deltat, maxval(solc(1,1)%u(cval%s)), &
      maxval(abs(sol%u(velx%s))), maxval(abs(sol%u(vely%s)))

    !call toc ( 'one step' )

  end do


! gradient tensor

  call create_vector ( problem, gradient_tensor, vec=6 )

  call derive_vector ( mesh, problem, gradient_tensor, &
    elemsub=deriv_gradient_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problem, filename='cylinder.vtk', &
    dataname='gradient_tensor', vector=gradient_tensor, symmetric=.false. )

  call delete ( gradient_tensor )


! conformation tensor

  call create_vector ( problemc, conformation_tensor, vec=3 )

  call derive_vector ( mesh, problemc, conformation_tensor, &
    elemsub=deriv_conformation_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemc, filename='cylinder.vtk', &
    dataname='conformation_tensor', vector=conformation_tensor, append=.true. )

  call delete ( conformation_tensor )


! create vector and copy data for input of projection problem

  call create ( problemc, vec_p, vec=4 )

  call transfer_data ( mesh, problem, problemc, sysvector1=sol, vector2=vec_p, &
    physq1=[1], degfd1=[1,2,3,4], degfd2=[1,2,3,4] )

  do icomp = 1, ncompc
    call transfer_data ( mesh, problemc, sysvector1=solc(icomp,1), &
      vector2=vec_p, degfd1=[1], degfd2=[icomp+4] )
  end do


!****************************************************************
! PROJECTION PROBLEM
!****************************************************************

! read mesh from gmsh output file

  call read_mesh_gmsh ( mesh1, filename='confined_cylinder2.msh', ndim=2 )
  call mesh_convert ( mesh1, mesh2, remove_isolated_nodes=.true. )
  call delete ( mesh1 )

  call fill_mesh_parts ( mesh2 )

  call printinfo ( mesh2, printlevel=4 )

! problem definition for projection

  call create_input_probdef ( mesh2, input_probdef2, nvec=1 )

  input_probdef2%elementdof(1)%a = 0
  input_probdef2%elementdof(1)%a(vertices) = 1  ! one component! one component
  input_probdef2%vec_elementdof(1)%a(:,1) = 1  ! plotting

  call problem_definition ( input_probdef2, mesh2, problem2 )

! create system vectors (solution and right-hand side)

  call create ( problem2, sol2 )
  call create ( problem2, rhsd2 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix2, mesh2, problem2, &
    symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix2 )

! fill coefficients

  call create_coefficients ( coefficients2, ncoefi=100, ncoefr=50 )

  coefficients2%i = 0
  coefficients2%i(1) = cintpl
  coefficients2%i(2) = ncompv
  coefficients2%i(3) = uintpl
  coefficients2%i(4) = cintpl
  coefficients2%i(10) = gaussp

  coefficients2%r = 0

! oldvectors

  call create_oldvectors ( oldvectors2, nsysvec=1, nvec=1, nprob=1, nmesh=1 )

  oldvectors2%v(1)%p => vec_p
  oldvectors2%p(1)%p => problemc
  oldvectors2%m(1)%p => mesh

! build (assemble) matrix and vector from elements

  call build_system ( mesh2, problem2, sysmatrix2, msysvector=rhsd2, &
    elemsub=projection_elem, coefficients=coefficients2, &
    oldvectors=oldvectors2 )

  do i = 1, ncompv

    call add_effect_of_essential_to_rhs ( problem2, sysmatrix2, sol2(i), &
      rhsd2(i) )

    call solve_system_ma57 ( sysmatrix2, rhsd2(i), sol2(i), lu=lu2 )

  end do

  call delete ( lu2 )  ! remove LU decomposition and rebuild next time step

! plot projection of Gxx
! (this is a way to quickly plot a projection vector, but has been
! commented out in favor of redefining the problem and plotting).
!
!  call create ( problem2, Gxx, vec=1 )
!
!  oldvectors2%s(1)%p => sol2(1)
!
!  call derive_vector ( mesh2, problem2, Gxx, elemsub=deriv_projection, &
!    coefficients=coefficients2, oldvectors=oldvectors2 )
!
!  call write_scalar_vtk ( mesh2, problem2, filename='cylinder2.vtk', &
!    dataname='Gxx', vector=Gxx )
!
!  call delete ( Gxx )


! Redefine the problem on mesh2 for the gradient tensor for plotting

  call delete ( sol )
  call delete ( input_probdef )
  call delete ( problem )

! problem definition of gradient

  call create_input_probdef ( mesh, input_probdef, nvec=2, nphysq=1 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 4  ! gradients
  input_probdef%vec_elementdof(1)%a(:,2) = 4  ! unsymmetric tensor

  input_probdef%physq = [1]
  input_probdef%probnr = 1

  call problem_definition ( input_probdef, mesh2, problem )

  call create_subscript ( mesh2, problem, gxxval, physqarr=[1], degfd=1 )
  call create_subscript ( mesh2, problem, gxyval, physqarr=[1], degfd=2 )
  call create_subscript ( mesh2, problem, gyxval, physqarr=[1], degfd=3 )
  call create_subscript ( mesh2, problem, gyyval, physqarr=[1], degfd=4 )

  call create ( problem, sol )

  sol%u(gxxval%s) = sol2(1)%u
  sol%u(gxyval%s) = sol2(2)%u
  sol%u(gyxval%s) = sol2(3)%u
  sol%u(gyyval%s) = sol2(4)%u

! gradient tensor

  call create_vector ( problem, gradient_tensor, vec=2 )

  call derive_vector ( mesh2, problem, gradient_tensor, &
    elemsub=deriv_gradient_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh2, problem, filename='cylinder2.vtk', &
    dataname='gradient_tensor', vector=gradient_tensor, symmetric=.false. )

  call delete ( gradient_tensor )


! Redefine the problem on mesh2 for the conformation tensor for plotting

  call delete ( solc )
  call delete ( input_probdefc )
  call delete ( problemc )

! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=3, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a(:,1) = 0
  input_probdefc%vec_elementdof(1)%a(vertices,1) = 1  ! c
  input_probdefc%vec_elementdof(1)%a(:,2) = 1         ! scalar for plotting
  input_probdefc%vec_elementdof(1)%a(:,3) = 3         ! symmetric tensor

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

  call problem_definition ( input_probdefc, mesh2, problemc )

  call create_subscript ( mesh2, problemc, cval, physqarr=[1] )

  call create ( problemc, solc )

  solc(1,1)%u(cval%s) = sol2(5)%u
  solc(2,1)%u(cval%s) = sol2(6)%u
  solc(3,1)%u(cval%s) = sol2(7)%u

! conformation tensor

  call create_vector ( problemc, conformation_tensor, vec=3 )

  call derive_vector ( mesh2, problemc, conformation_tensor, &
    elemsub=deriv_conformation_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh2, problemc, filename='cylinder2.vtk', &
    dataname='conformation_tensor', vector=conformation_tensor, append=.true. )

  call delete ( conformation_tensor )


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
  call delete ( lu_u )
  call delete ( solc, rhsc )

  call delete ( coefficients )
  call delete ( velx, vely, cval )
  call delete ( gxxval, gxyval, gyxval, gyyval )

  call delete ( mesh2 )
  call delete ( problem2 )
  call delete ( input_probdef2 )
  call delete ( sol2, rhsd2 )
  call delete ( sysmatrix2 )
  call delete ( oldvectors2 )
  call delete ( coefficients2 )

end program projection4
