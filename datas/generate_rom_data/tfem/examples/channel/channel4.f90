! 3D viscoelastic problem in a channel with a square cross section.
! Computed is the developed flow using only one layer of elements in the flow
! direction (x) and assuming periodical boundary conditions.
!
! DEVSS-G/SUPG
! first-order time integration to steady state
!
! This problem is the same as channel2, but now using coupling of objects
! in realizing the periodical boundary conditions.

module subs_m

  use tfem_elem_m
  implicit none


contains


! Element for constraint on Q2 interpolations (velocities)
! (connection through collocation on objects)

  subroutine constraints_object_conn_Q2 ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


    integer :: object, object2, i
    integer, parameter :: ndf = 27
    real(dp) :: phi(1,ndf), phi2(1,ndf), xr(1,3), xr2(1,3)

!   connection through collocation on objects

    object = problem%constraints(constr)%object
    object2 = problem%constraints(constr)%object2

    xr(1,:) = mesh%objects(object)%refcoor(node,:)
    xr2(1,:) = mesh%objects(object2)%refcoor(node,:)

    call shape_hexa_Q2 ( xr, phi )
    call shape_hexa_Q2 ( xr2, phi2 )

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

      elemmat = 0
      do i = 1, size(elemmat,1)
        elemmat(i,ndf*(i-1)+1:i*ndf) = phi(1,:)
      end do

      elemmat2 = 0
      do i = 1, size(elemmat2,1)
        elemmat2(i,ndf*(i-1)+1:i*ndf) = -phi2(1,:)
      end do
    end if


  end subroutine constraints_object_conn_Q2


! Element for constraint on the Q1 interpolations (gradients, conformation)
! (connection through collocation on objects)

  subroutine constraints_object_conn_Q1 ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


    integer :: object, object2, i
    integer, parameter :: ndf = 8
    real(dp) :: phi(1,ndf), phi2(1,ndf), xr(1,3), xr2(1,3)

!   connection through collocation on objects

    object = problem%constraints(constr)%object
    object2 = problem%constraints(constr)%object2

    xr(1,:) = mesh%objects(object)%refcoor(node,:)
    xr2(1,:) = mesh%objects(object2)%refcoor(node,:)

    call shape_hexa_Q1_reg ( xr, phi )
    call shape_hexa_Q1_reg ( xr2, phi2 )

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

      elemmat = 0
      do i = 1, size(elemmat,1)
        elemmat(i,ndf*(i-1)+1:i*ndf) = phi(1,:)
      end do

      elemmat2 = 0
      do i = 1, size(elemmat2,1)
        elemmat2(i,ndf*(i-1)+1:i*ndf) = -phi2(1,:)
      end do

    end if

  end subroutine constraints_object_conn_Q1

end module subs_m

program channel4

  use tfem_m
  use hsl_ma41_m
  use hsl_ma57_m
  use viscoelastic_elements_m
  use figplot_m
  use io_utils_m
  use subs_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities
    pintpl = 4,     & ! Q1 pressures
    gintpl = 4,     & ! Q1 gradients
    cintpl = 4,     & ! Q1 conformation
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    gauss = 3,      & ! 3x3 integration of quads
    ncompc = 6,     & ! number of conformation tensor components
    nmodes = 1,     & ! number of modes
    startm = 501      ! start of material model data

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh1
  type(input_probdef_t) :: input_probdef, input_probdefc
  type(problem_t), target :: problem, problemc
  type(sysmatrix_t) :: sysmatrix, sysmatrixc
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: velx, vely, velz, cval
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc
  type(lu_ma57_t) :: lu_u
  type(lu_ma41_t) :: luc
  type(plot_options_t) :: plot_options
  type(solver_options_ma57_t) :: solver_options_u
  type(solver_options_ma41_t) :: solver_options_c

! variables

  integer :: &
    nx=1,                & ! number of elements in x
    ny=10,               & ! number of elements in y
    nz=10,               & ! number of elements in z
    timeint = 1,         & ! (first-order) time integration
    numtimesteps = 2,    & ! number of time steps
    logc = 0,            & ! standard scheme or log transformation
    model = 2              ! Oldroyd-B

  real(dp) :: &
    lx = 0.1_dp,       & ! size in x-direction
    ly = 1._dp,        & ! size in y-direction
    lz = 1._dp,        & ! size in z-direction
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 1.0_dp,   & ! relaxation time
    alphapar = 0.1_dp, & ! alpha parameter in the Giesekus model
    deltat = 5.e-3_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    U = 1._dp,         & ! imposed average velocity
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.0_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 1.4_dp       ! integer_storage for conformation LU (HSL)


  integer :: icomp, step, i, npar
  integer :: vertices(8) = [1,3,9,7,19,21,27,25]
  real(dp) :: alpha, G, flowrate

  logical, parameter :: connect_gradients = .false.


! namelist for input of variables; read from standard input

  namelist /comppar/ nx, ny, nz, lx, ly, lz, timeint, numtimesteps, logc, &
    model, eta_s, eta_p, lambda, alphapar, deltat, beta, U, &
    rs_gup, is_gup, rs_c, is_c

  read ( unit=*, nml=comppar )

! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus

  flowrate = U * ly * lz

  if ( model == 2 ) then
    npar = 2 ! Oldroyd-B
  else if ( model == 3 ) then
    npar = 3 ! Giesekus
  end if


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+npar*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gauss,    cintpl,     0,     0,         0,      &
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

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, &
    lz=lz, elshape=14, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )


! additional objects for connecting velocities

  call add_to_mesh ( mesh, object='surface', objectsurface=5, &
    excludecurves=[4,5,12,8] )

  call add_to_mesh ( mesh, object='surface', objectsurface=5, &
    excludecurves=[4,5,12,8] )

! shift object 2 from back to front
  mesh%objects(2)%coor(:,1) = mesh%objects(2)%coor(:,1) + lx


! additional objects for connecting gradients and conformation

! temporary mesh using Q1 elements

  meshgen_options%elshape=13

  call hexahedron ( mesh1, meshgen_options )

  call add_to_mesh ( mesh, object='surface', objectsurface=5, objectmesh=mesh1 )

  call add_to_mesh ( mesh, object='surface', objectsurface=5, objectmesh=mesh1 )

! shift object 4 from back to front
  mesh%objects(4)%coor(:,1) = mesh%objects(4)%coor(:,1) + lx

  call delete ( mesh1 ) ! delete temporary mesh


  call fill_mesh_parts ( mesh )


! plot curves, surfaces and mesh

  plot_options%viewpoint=[1.,0.8,0.4]
  plot_options%fontsize=10
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  plot_options%plotboundary=.false.
  call plot_objects ( plot_options, mesh, 'curves.fig', append=.true. )
  plot_options%plotboundary=.true.
  call plot_mesh ( plot_options, mesh, 'mesh.fig', surfaces=[3,4,6] )


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=3 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 9  ! gradients
  input_probdef%vec_elementdof(1)%a(:,2) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

  call define_essential ( mesh, input_probdef, surface1=1, surface2=2, &
    physq=physqvel )
  call define_essential ( mesh, input_probdef, surface1=4, physq=physqvel )
  call define_essential ( mesh, input_probdef, surface1=6, physq=physqvel )
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, surface1=5, nglobalc=1 )

! constraints for periodical boundary conditions

! velocities
  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, object=1, object2=2, discretization='collocation', &
    nodedof=3 )

  if ( connect_gradients ) then

!   gradients
    call define_constraint ( mesh, input_probdef, &
      physq=physqgrad, object=3, object2=4, discretization='collocation', &
      nodedof=9 )

  end if

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for the cross velocity for post processing

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, velz, physqarr=[physqvel], degfd=3 )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=2, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a(:,1) = 0
  input_probdefc%vec_elementdof(1)%a(vertices,1) = 1  ! c
  input_probdefc%vec_elementdof(1)%a(:,2) = 1         ! scalar for plotting

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

! constraint for periodical boundary conditions of the conformation

  call define_constraint ( mesh, input_probdefc, &
    object=3, object2=4, discretization='collocation', nodedof=1 )

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

! periodical condition on velocities and gradients

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=constraints_object_conn_Q2, &
    addmatvec=.true., coefficients=coefficients )

  if ( connect_gradients ) then

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=3, elemsub=constraints_object_conn_Q1, &
      addmatvec=.true., coefficients=coefficients )

  end if

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
    solc(3,1)%u = 0
    solc(4,1)%u = 1
    solc(5,1)%u = 0
    solc(6,1)%u = 1
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


! time stepping

  do step = 1, numtimesteps


!   build (assemble) matrix and vector for conformation problem

    call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
      elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
      coefficients=coefficients )

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problemc, sysmatrixc, &
      m2sysvector=rhsc, elemsub=constraints_object_conn_Q1, &
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
      maxval(abs(sol%u(velx%s))), maxval(abs(sol%u(vely%s))), &
      maxval(abs(sol%u(velz%s)))

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
  write(10) (solc(i,1)%u, i=1,ncompc)

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
  call delete ( lu_u )
  call delete ( solc, rhsc )

  call delete ( coefficients )
  call delete ( velx, vely, velz, cval )


end program channel4
