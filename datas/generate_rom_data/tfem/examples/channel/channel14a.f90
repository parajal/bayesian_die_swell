! Viscoelastic problem in a curved channel with a square cross section.
! 3D velocities on a 2D domain (developed flow). Axisymmetric.
! DEVSS-G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! second-order time integration to steady state
! Similar to channel11 but now with a b-tensor formulation.

program channel14a

  use tfem_m
  use hsl_ma41_m
  use hsl_ma57_m
  use viscoelastic_elements_m
  use figplot_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities
    pintpl = 4,     & ! Q1 pressures
    gintpl = 4,     & ! Q1 gradients
    bintpl = 4,     & ! Q1 b-tensor
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    gauss = 3,      & ! 3x3 integration of quads
    ncompc = 6,     & ! number of b-tensor components
    ncompb = 9,     & ! number of b-tensor components
    nmodes = 1,     & ! number of modes
    startm = 501      ! start of material model data

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefb, input_probdef_proj
  type(problem_t), target :: problem, problemb, problem_proj
  type(sysmatrix_t) :: sysmatrix, sysmatrixb, sysmatrix_proj
  type(sysvector_t), target :: sol, solm1
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, gammadot, Dtensor, Ltensor
  type(vector_t) :: tauviscous, tauviscoelastic, ctensor, btensor
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: velx, vely, velz, bval
  type(sysvector_t), dimension(ncompb,nmodes), target :: solb, solbm1
  type(sysvector_t), dimension(ncompb,nmodes) :: rhsb
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc_proj
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc_proj
  type(lu_ma41_t) :: lub
  type(plot_options_t) :: plot_options
  type(solver_options_ma41_t) :: solver_options_u, solver_options_b
  type(lu_ma57_t) :: lu_proj

! variables

  logical :: &
    cproj = .true.         ! projection of c=b*b^T for cn in momentum balance

  integer :: &
    nx=10,               & ! number of elements in z
    ny=10,               & ! number of elements in r
    timeint1 = 1,        & ! first-order time integration, first time step
    timeint2 = 7,        & ! second-order time integration after first time step
    numtimesteps = 2,    & ! number of time steps
    bvariant = 1,        & ! b-formulation:
                           ! 1: CDT 2: symmetric
                           ! 3: Cholesky 4: Cholesky with log
    model = 2              ! Oldroyd-B

  real(dp) :: &
    oy = 1._dp,        & ! left lower corner in r-direction
    lx = 1._dp,        & ! size in z-direction
    ly = 1._dp,        & ! size in r-direction
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 1.0_dp,   & ! relaxation time
    alphapar = 0.1_dp, & ! alpha parameter in the Giesekus model
    deltat = 5.e-3_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    U = 1._dp,         & ! imposed average velocity
    rs_gup = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_b  = 1.0_dp,    & ! real_storage for the conformation LU (HSL)
    is_b  = 1.4_dp       ! integer_storage for conformation LU (HSL)


  integer :: icomp, step, i, npar
  integer :: vertices(4) = [1,3,5,7]
  real(dp) :: alpha, G, flowrate


! namelist for input of variables; read from standard input

  namelist /comppar/ nx, ny, lx, ly, oy, timeint1, timeint2, numtimesteps, &
    bvariant, model, eta_s, eta_p, lambda, alphapar, deltat, beta, U, &
    cproj, rs_gup, is_gup, rs_b, is_b

  read ( unit=*, nml=comppar )

! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus

  flowrate = U * lx * ly

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
      gauss,    bintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      0,        timeint1,   ( 0, i = 23, 150 )  &
    ]

  coefficients%i(23) = 1  ! coorsys = 1, cilindrical coordinates, axisymmetric
  coefficients%i(61) = 0  ! projected G=0, direct velocity gradient=1 in CE
  coefficients%i(67) = 1  ! 3D velocity
  coefficients%i(71) = bvariant  ! variant for b-formulation
  if ( cproj ) coefficients%i(72) = 1  ! use c projection in momentum balance
  coefficients%i(90) = 1  ! rotation reinitialization on element level

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

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, lx=lx, ly=ly, &
    oy=oy, elshape=6 )

  call quadrilateral2d ( mesh, meshgen_options )

! make domain into a surface for the flow rate constraint
  call add_to_mesh ( mesh, surfacefromgroups=[1] )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=6, nphysq=3 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 6  ! gradients
  input_probdef%vec_elementdof(1)%a(:,2) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5) = 6  ! symmetric tensor
  input_probdef%vec_elementdof(1)%a(:,6) = 9  ! unsymmetric tensor

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

  call define_essential ( mesh, input_probdef, &
    curve1=1, curve2=4, physq=physqvel )
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, surface1=1, nglobalc=1 )

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for the cross velocity for post processing

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
  call create_subscript ( mesh, problem, velz, physqarr=[physqvel], degfd=3 )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefb, nvec=4, nphysq=1 )

  input_probdefb%vec_elementdof(1)%a(:,1) = 0
  input_probdefb%vec_elementdof(1)%a(vertices,1) = 1  ! b
  input_probdefb%vec_elementdof(1)%a(:,2) = 1         ! scalar for plotting
  input_probdefb%vec_elementdof(1)%a(:,3) = 6         ! symmetric tensor
  input_probdefb%vec_elementdof(1)%a(:,4) = 9         ! unsymmetric tensor

  input_probdefb%physq = [1]
  input_probdefb%probnr = 2

  call problem_definition ( input_probdefb, mesh, problemb )

! create a vector subscript for the b-tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemb, bval, physqarr=[1] )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, solm1, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0

! create system matrix of gradient/velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! create system vectors (solution and right-hand side) for b-tensor and
! initialize vectors to zero stress

  call create ( problemb, solb, solbm1, rhsb )

  if ( bvariant == 4 ) then
!   log Cholesky
    do i = 1, ncompb
      solb(i,1)%u = 0
    end do
  else
    do i = 1, ncompb
      solb(i,1)%u = 0
    end do
    solb(1,1)%u = 1 ! initial bxx
    solb(5,1)%u = 1 ! initial byy
    solb(9,1)%u = 1 ! initial bzz
  end if

! create system matrix for b-tensor problem

  call create_sysmatrix_structure ( sysmatrixb, mesh, problemb )
  call create_sysmatrix_data ( sysmatrixb )


! problem definition for projected "c=b*b^T" of the b-tensor

  call create_input_probdef ( mesh, input_probdef_proj, nvec=1, nphysq=1 )

  input_probdef_proj%vec_elementdof(1)%a(:,1) = 0
  input_probdef_proj%vec_elementdof(1)%a(vertices,1) = 1  ! c

  input_probdef_proj%physq = [1]
  input_probdef_proj%probnr = 3

  call problem_definition ( input_probdef_proj, mesh, problem_proj )

  call create ( problem_proj, solc_proj, rhsc_proj )

! initialize solc_proj
  solc_proj(1,1)%u = 1
  solc_proj(2,1)%u = 0
  solc_proj(3,1)%u = 0
  solc_proj(4,1)%u = 1
  solc_proj(5,1)%u = 0
  solc_proj(6,1)%u = 1


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec2=3, nprob=3 )

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


! open monitor file

  open ( unit=13, file='out', recl=300 )


! time stepping

  do step = 1, numtimesteps

    if ( step == 1 ) then
      coefficients%i(78) = 1 ! avoid dividing by zero in SUPG in first step
      coefficients%r(10) = U
    else if ( step == 2 ) then
      coefficients%i(78) = 0
    end if

    if ( step == 2 ) then
!     change time integration scheme at the second time step
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

    call copy ( sol, solm1 )


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


!   solve conformation and keep LU decomposition in the loop over components

    solver_options_b%real_storage=rs_b
    solver_options_b%integer_storage=is_b

    do icomp = 1, ncompb
      call solve_system_ma41 ( sysmatrixb, rhsb(icomp,1), solb(icomp,1), lub, &
        solver_options=solver_options_b  )
    end do

    call delete ( lub )  ! remove LU decomposition and rebuild next time step

!   write monitor data

    write(unit=13,fmt=*) &
      step, step*deltat, maxval(solb(1,1)%u(bval%s)), &
      maxval(abs(sol%u(velx%s))), maxval(abs(sol%u(vely%s))), &
      maxval(abs(sol%u(velz%s)))

  end do


! close monitor data file

  close(unit=13)

! post-processing

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  call create_vector ( problem, velocity, physq=2 )
  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, gammadot, vec=4 )
  call create_vector ( problem, Dtensor, vec=5 )
  call create_vector ( problem, Ltensor, vec=6 )
  call create_vector ( problem, tauviscous, vec=5 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=11

  call derive_vector ( mesh, problem, gammadot, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Dtensor, elemsub=stokes_D_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Ltensor, elemsub=stokes_gradu_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, tauviscous, elemsub=stokes_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='channel14.vtk' )

  call write_vector_vtk ( mesh, problem, filename='channel14.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true., &
    assume3D=.true. )

  call write_vector_vtk ( mesh, problem, filename='channel14.vtk', &
    dataname='velocity_vector_planar', vector=velocity, append=.true., &
    degfd=[1,2] )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='channel14.vtk', dataname='gammadot', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='channel14.vtk', &
    dataname='D', vector=Dtensor, append=.true., assume3D=.true. )

  call write_tensor_vtk ( mesh, problem, filename='channel14.vtk', &
    dataname='L', vector=Ltensor, append=.true., assume3D=.true., &
    symmetric=.false. )

  call write_tensor_vtk ( mesh, problem, filename='channel14.vtk', &
    dataname='tauviscous', vector=tauviscous, append=.true., assume3D=.true. )

  call create_vector ( problemb, tauviscoelastic, vec=3 )

  call derive_vector ( mesh, problemb, tauviscoelastic, &
    elemsub=deriv_viscoelastic_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemb, filename='channel14.vtk', &
    dataname='tauviscoelastic', vector=tauviscoelastic, append=.true., &
    assume3D=.true. )

  call create_vector ( problemb, ctensor, vec=3 )

  call derive_vector ( mesh, problemb, ctensor, &
    elemsub=deriv_conformation_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemb, filename='channel14.vtk', &
    dataname='c', vector=ctensor, append=.true., &
    assume3D=.true. )

  call create_vector ( problemb, btensor, vec=4 )

  call derive_vector ( mesh, problemb, btensor, &
    elemsub=deriv_conformation_tensor_std,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemb, filename='channel14.vtk', &
    dataname='b', vector=btensor, symmetric=.false., append=.true., &
    assume3D=.true. )

! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, solm1, rhsd )
  call delete ( sysmatrix )
  call delete ( velocity, pressure, gammadot, Dtensor, Ltensor, tauviscous, &
                tauviscoelastic, ctensor, btensor )
  call delete ( oldvectors, oldvectors_ve )
  call delete ( problemb )
  call delete ( input_probdefb )
  call delete ( sysmatrixb )
  call delete ( solb, solbm1, rhsb )

  call delete ( problem_proj )
  call delete ( input_probdef_proj )
  call delete ( solc_proj, rhsc_proj )
  if ( cproj ) call delete ( lu_proj )

  call delete ( coefficients )
  call delete ( velx, vely, velz, bval )

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

!   imposed flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr_surface, addmatvec=.true., &
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
      do i = 1, ncompc

        call add_effect_of_essential_to_rhs ( problem_proj, sysmatrix_proj, &
           sol_proj(i,m), rhs_proj(i,m) )

        call solve_system_ma57 ( sysmatrix_proj, rhs_proj(i,m), &
           sol_proj(i,m), lu_proj, solver_options=solver_options_ma57 )

      end do
    end do

  end subroutine solve_projection

end program channel14a
