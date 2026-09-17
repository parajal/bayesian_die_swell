#if SKIT2 && METIS5

! 3D viscoelastic problem: flow around a sphere confined between four walls.
! The entry and exit sides are connected by periodical boundary conditions
! and a constant flow rate is imposed.
! The average velocity is U.
! Drag computations.
!
! DEVSS-G/SUPG
! first-order time integration to steady state
!
! sparskit solver
! Metis renumbering


module subs_m

  use tfem_elem_m

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


    integer :: object, object2
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

program sphere12

  use tfem_m
  use subs_m
  use sk_solve_m
  use viscoelastic_elements_m
  use figplot_m
  use io_utils_m
  use timer_m
  use metis5_m

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
  type(oldvectors_t) :: oldvectors_ve, oldvectors, oldvectors_dve, oldvectors_d
  type(coefficients_t) :: coefficients
  type(subscript_t) :: velx, vely, velz, cval
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc
  type(vector_t), target :: pressure, stress_tensor, conformation_tensor(nmodes)
  type(vector_t) :: velocity
  type(ilu_sk_t) :: ilu
  type(ilu_sk_t) :: iluc
  type(plot_options_t) :: plot_options
  type(solver_options_sk_t) :: solver_options


! variables

  integer :: &
    nx = 4,              & ! number of element in the x-direction
    ny = 4,              & ! number of element in the y-direction
    restart = 1,         & ! restart
    timeint = 1,         & ! (first-order) time integration
    numtimesteps = 2,    & ! number of time steps
    logc = 0,            & ! standard scheme or log transformation
    htype = 0,           & ! upwind parameter
    Uscaling = 0,        & ! upwind parameter
    model = 2              ! Oldroyd-B

  real(dp) :: &
    ox = -1.5_dp,      & ! origin shift in x-direction
    oy = -1.5_dp,      & ! origin shift in y-direction
    oz = -1.5_dp,      & ! origin shift in z-direction
    lx = 4._dp,        & ! size in x-direction
    ly = 4._dp,        & ! size in y-direction
    lz = 4._dp,        & ! size in z-direction
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

  logical, parameter :: connect_gradients = .false.

  integer :: icomp, step, i, m, npar, ios
  integer :: vertices(8) = [1,3,9,7,19,21,27,25]
  real(dp) :: alpha, G, flowrate
  real(dp) :: dragv(3), dragve(3)



! namelist for input of variables; read from standard input

  namelist /comppar/ restart, nx, ny, timeint, numtimesteps, logc, &
    model, eta_s, eta_p, lambda, alphapar, deltat, &
    beta, htype, Uscaling, U, rs_gup, is_gup, rs_c, is_c

  read ( unit=*, nml=comppar )

! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus

  flowrate = - U * lx * ly

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
      gauss,    cintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      logc,     timeint,    ( 0, i = 23, 150 )  &
    ]

  coefficients%r(1:502) = &
    [ eta_s,    0._dp,   0._dp,   alpha, 0._dp, &
      flowrate, 0._dp,  deltat,   beta,      U, &
      ( 0._dp, i = 11, 500 ), &
      G,     lambda  &
    ]

  if ( model == 3 ) then
    coefficients%r(503) = alphapar ! Giesekus
  end if

  coefficients%i(29) = htype
  coefficients%i(31) = Uscaling


! read mesh

  call read_mesh_sepran ( mesh, filename='sphere7.mout' )

  call add_to_mesh ( mesh, object='surface', objectsurface=11, &
    excludecurves=[13,17,18,21] )

  call add_to_mesh ( mesh, object='surface', objectsurface=11, &
    excludecurves=[13,17,18,21] )


! shift object 2 from bottom to top
  mesh%objects(2)%coor(:,3) = mesh%objects(2)%coor(:,3) + lz


! temporary 3D mesh using Q1 elements

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=1, lx=lx, &
    ly=ly, ox=ox, oy=oy, oz=oz, elshape=13, regionshape=4 )

  call hexahedron ( mesh1, meshgen_options )

  call add_to_mesh ( mesh, object='surface', objectsurface=1, objectmesh=mesh1 )

  call add_to_mesh ( mesh, object='surface', objectsurface=1, objectmesh=mesh1 )

  call delete ( mesh1 )


! shift object 4 from bottom to top
  mesh%objects(4)%coor(:,3) = mesh%objects(4)%coor(:,3) + lz


  call fill_mesh_parts ( mesh )

  call renumber_metis ( mesh )

! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=3, &
    nphysqshifted=1 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 9  ! gradients
  input_probdef%vec_elementdof(1)%a(:,2) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5) = 6  ! symmetric tensor

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1
  input_probdef%physqshifted = [3]

! define essential boundaries

! sphere
  call define_essential ( mesh, input_probdef, surface1=25, physq=physqvel )

! ends and walls
  call define_essential ( mesh, input_probdef, surface1=7, surface2=10, &
    physq=physqvel )

! pressure level
  call define_essential ( mesh, input_probdef, point=11, physq=physqpress )


! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, surface1=11, nglobalc=1 )


! constraints for periodical boundary conditions

! velocities

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, object=1, object2=2, discretization='collocation', &
    nodedof=3 )

! gradients

  if ( connect_gradients ) then

  call define_constraint ( mesh, input_probdef, &
    physq=physqgrad, object=3, object2=4, discretization='collocation', &
    nodedof=9 )

  end if

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for the velocity for post processing

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, velz, physqarr=[physqvel], degfd=3 )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=3, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a(:,1) = 0
  input_probdefc%vec_elementdof(1)%a(vertices,1) = 1  ! c
  input_probdefc%vec_elementdof(1)%a(:,2) = 1         ! scalar for plotting
  input_probdefc%vec_elementdof(1)%a(:,3) = 6         ! symmetric tensor

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

! constraint for periodical boundary conditions of the conformation

  call define_constraint ( mesh, input_probdefc, &
    object=3, object2=4, discretization='collocation' , nodedof=1)

  call problem_definition ( input_probdefc, mesh, problemc )

! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemc, cval, physqarr=[1] )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )


! fill solution vector with essential boundary conditions

  sol%u = 0

  call fill_sysvector ( mesh, problem, sol, &
    surface1=7, surface2=10, physq=physqvel, degfd=3, value=U )

! create system matrix of gradient/velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
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
    constraint1=2, elemsub=constraints_object_conn_Q2, addmatvec=.true., &
    coefficients=coefficients )

  if ( connect_gradients ) then

!   periodical condition on gradients

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=3, elemsub=constraints_object_conn_Q1, addmatvec=.true., &
      coefficients=coefficients )

  end if

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

! solve initial gradient/velocity/pressure consistent with initial conditions
! and keep decomposition of the matrix
! This generates a Stokes profile.

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

! set parameters for iterative solver (for both v/p/G and c system)

  call set_solver_options ( solver_options, printlevel=1, type_prec=2, &
    maxmvm=300, itsolver=8, mgmres=30, prec_store=1._dp, preconditioner=1, &
    droptol=1e-3_dp, fillin=0.5_dp, eps_rel = 1e-7_dp, eps_abs = 1.e-3_dp )

  call solve_system_sk ( sysmatrix, rhsd, sol, ilu, solver_options=solver_options )

  print *, 'solved vp '


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


! Do a restart from previous data?

  if ( restart == 1 ) then

!   read data

    open ( unit = 10, file='data.out', form='unformatted', iostat=ios, &
           status='old' )

    if ( ios /= 0 ) then
      write(*,'(/2a/)') 'Error: cannot open file data.out '
      stop
    end if

    read(10) sol%u
    read(10) (solc(i,1)%u, i=1,ncompc)

    close ( unit=10 )

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

    print *,' step =', step

!   build (assemble) matrix and vector for conformation problem

    call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
      elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
      coefficients=coefficients )

!   periodical condition on conformation tensor
    call build_system_constraint ( mesh, problemc, sysmatrixc, &
      m2sysvector=rhsc, elemsub=constraints_object_conn_Q1, &
      coefficients=coefficients , addmatvec=.true. )

    call check ( sysmatrixc )

!   solve conformation and keep LU decomposition in the loop over components

    print *, 'start solve for solc'

    do icomp = 1, ncompc
      call solve_system_sk ( sysmatrixc, rhsc(icomp,1), solc(icomp,1), iluc, &
        initsol=.true., solver_options=solver_options )
    end do

    print *, 'end solve for solc'

    call delete ( iluc )  ! remove ILU decomposition and rebuild next time step


!   build (assemble) vector for gradient/velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_divtau, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      buildmatrix = .false., coefficients=coefficients )

!   imposed flow rate (right-hand side only)

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      buildmatrix = .false., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


    print *, 'start solve for sol'

!   solve gradient/velocity/pressure problem

    call solve_system_sk ( sysmatrix, rhsd, sol, ilu, initsol=.true., &
      solver_options=solver_options )

    print *, 'end solve for sol'

!   write monitor data

    write(unit=13,fmt=*) &
      step, step*deltat, maxval(solc(1,1)%u(cval%s)), &
      maxval(abs(sol%u(velx%s))), maxval(abs(sol%u(vely%s))), &
      maxval(abs(sol%u(velz%s)))

    !call toc ( 'one step' )

  end do


! close monitor data file

  close(unit=13)

! post processing

! write data for post-processing

  call write_mesh ( mesh, filename='mesh.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh, input_probdefc, filename='probdefc.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol%u
  write(10) (solc(i,1)%u, i=1,ncompc)

  close(unit=10)
!

! post processing

  call create_vector ( problem, velocity, physq=2 )
  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, stress_tensor, vec=5 )
  call extract_physvector ( mesh, problem, sol, velocity )


! fill oldvectors for stokes problem

  call create_oldvectors ( oldvectors, nsysvec=1 )
  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, stress_tensor, &
    elemsub=stokes_stress_tensor, coefficients=coefficients, &
    oldvectors=oldvectors, elsurfaces=[25] ) ! only include surface 25


! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='sphere.vtk' )

  call write_vector_vtk ( mesh, problem, filename='sphere.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true. )
!  call write_vector_vtk ( mesh, problem, filename='sphere.vtk', &
!    dataname='velocity_vector', sysvector=sol, physq=physqvel, append=.true. )

  call delete ( velocity, pressure, stress_tensor )


  print *, ' create output for conformation tensor'

! create oldvectors for viscoelastic problem

  call create_oldvectors ( oldvectors_dve, nsysvec2=1 )

  oldvectors_dve%s2(1)%p => solc

! conformation tensor

  call create ( problemc, conformation_tensor, vec=3 )

  call derive_vector ( mesh, problemc, conformation_tensor(1), &
    elemsub=deriv_conformation_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_dve )

  call write_tensor_vtk ( mesh, problemc, filename='sphere.vtk', &
    dataname='conformation_tensor', vector=conformation_tensor(1), &
    append=.true. )

  call delete ( conformation_tensor )


! fill oldvectors for stokes problem

  call create_oldvectors ( oldvectors_d, nsysvec=1, nvec=2 )

  oldvectors_d%s(1)%p => sol

! pressure

  call create_vector ( problem, pressure, vec=4 )

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors_d )

! viscous stress

  call create_vector ( problem, stress_tensor, vec=5 )

  call derive_vector ( mesh, problem, stress_tensor, &
    elemsub=stokes_stress_tensor, coefficients=coefficients, &
    oldvectors=oldvectors_d, elsurfaces=[25] ) ! only include surface 25

  oldvectors_d%v(1)%p => stress_tensor
  oldvectors_d%v(2)%p => pressure

! integrate drag force

  call integrate_boundary_elements ( mesh, problem, dragv, &
    elemsub=stokes_drag, surface=25, coefficients=coefficients, &
    oldvectors=oldvectors_d )

  !print *, 'viscous drag = ', dragv(1)
  print *, 'viscous drag = ', dragv

  call delete ( pressure, stress_tensor )

! create oldvectors for viscoelastic problem

  call delete ( oldvectors_dve )
  call create_oldvectors ( oldvectors_dve, nsysvec2=1, nvec1=1 )

  oldvectors_dve%s2(1)%p => solc

! conformation tensor in the nodes

  call create ( problemc, conformation_tensor, vec=3 )

  call derive_vector ( mesh, problemc, conformation_tensor(1), &
    elemsub=deriv_conformation_tensor_std, coefficients=coefficients, &
    oldvectors=oldvectors_dve, elsurfaces=[25] ) ! only include surface 25

  oldvectors_dve%v1(1)%p => conformation_tensor

! integrate drag force

  call integrate_boundary_elements ( mesh, problemc, dragve, &
    elemsub=viscoelastic_drag, surface=25, coefficients=coefficients, &
    oldvectors=oldvectors_dve )

  !print *, 'viscoelastic drag = ', dragve(1)
  !print *, 'total drag = ', dragv(1)+dragve(1)
  print *, 'viscoelastic drag = ', dragve
  print *, 'total drag = ', dragv+dragve

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
  call delete ( ilu )
  call delete ( solc, rhsc )

  call delete ( coefficients )
  call delete ( velx, vely, velz, cval )


end program sphere12

#else
  print '(4(a/),a)', &
    'To run this example:', &
    ' - compile add-ons sk_solve and metis5', &
    ' - install and compile sparskit2', &
    ' - add metis5 lib for linking', &
    ' - set preprocessing macro SKIT2 and METIS5 in Mdefs.mk'
end
#endif
