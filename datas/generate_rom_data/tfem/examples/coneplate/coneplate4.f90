! Viscoelastic cone-plate problem with a free surface.
! Compute the torque using both the stress integral and the reaction forces.
! (2D) domain with 3D velocities (axisymmetric).
! Cylindrical coordinates with circumferential velocity included (swirl).
! DEVSS-G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! L2-projection of c=exp(s) on the finite element approximation space in the
! momentum balance.
! second-order time integration

program coneplate4

  use tfem_m
  use math_defs_m
  use hsl_ma41_m
  use hsl_ma57_m
  use viscoelastic_elements_m
  use surface_advection_elements_m
  use io_utils_m
  use torque_elements_m
  use update_mesh_nodes_cp_m
!m  use metis_5_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,  & ! P2 velocities
    pintpl = 2,  & ! P1 pressures
    gintpl = 2,  & ! P1 gradients
    cintpl = 2,  & ! P1 conformation
    hintpl = 6,  & ! P2 height function
    gauss  = 5,  & ! order of Gauss integration of triangles
    gaussb = 5,  & ! number of integration points for boundary elements
    ncompc = 6,  & ! number of conformation tensor components
    ncompg = 6,  & ! number of velocity gradient tensor components swirl
    nmodes = 1,  & ! number of modes
    startm = 501   ! start of material model data

  real(dp), parameter :: &
    Ro = 1.0_dp,      & ! outer radius (of fluid contact point on the plate)
    theta_cp = 0.1_dp,& ! Cone angle
    dx_Ro = 0.01_dp,  & ! Element size at the fluid outside surface
    dx_0 = 0.004_dp     ! Element size at the center (r=0)

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc, input_probdefc_proj
  type(problem_t), target :: problem, problemc, problemc_proj
  type(sysmatrix_t) :: sysmatrix, sysmatrixc, sysmatrixc_proj
  type(sysvector_t), target :: sol, solm1, sol_p
  type(sysvector_t) :: rhsd
  type(sysvector_t) :: reacf
  type(vector_t) :: velocity, pressure, gammadot, Dtensor, Ltensor, &
    velocity_rot
  type(vector_t) :: tauviscous, tauviscoelastic, ctensor, ctensor_rot
  type(vector_t), target :: meshvel
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: velz, velr, velt, cval, velz3, velr3
  type(subscriptvec_t) :: szz, szr, szt, srr, srt, stt, mvelz3, mvelr3
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc, solcm1, &
    solc_proj, solc_p
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc, rhsc_proj
  type(lu_ma41_t) :: luc
  type(solver_options_ma41_t) :: solver_options_u, solver_options_c
  type(lu_ma57_t) :: lu_exps_proj
  type(sample_t) :: sample_p

  real(dp), allocatable, dimension(:,:) :: sample_coor

! variables

  integer :: &
    timeint1 = 1,        & ! first-order time integration, first time step
    timeint2 = 7,        & ! second-order time integration after first time step
    numtimesteps = 500,  & ! number of time steps
    vtkevery = 50,       & ! plot vtk every .... steps
    logc = 1,            & ! standard scheme or log transformation
    model = 3              ! 2: Oldroyd-B 3: Giesekus

  real(dp) :: &
    eta_s = 0.1_dp,      & ! solvent viscosity
    eta_p = 1.0_dp,      & ! polymer viscosity
    lambda = 0.2_dp,     & ! relaxation time
    alphapar = 0.1_dp,   & ! alpha parameter in the Giesekus model
    gamma_fa = 0.02_dp,  & ! surface tension fluid air interface
    dgamma_wa = -0.004_dp,& ! surface tension difference fluid air w.r.t. wall
    deltat = 0.04_dp,    & ! time step
    beta = 1.0_dp,       & ! upwinding parameter in the SUPG method
    betai = 0.5_dp,      & ! upwinding parameter in the SUPG method on the
                           ! free surface
    omega = 0.4_dp,      & ! Angular velocity of the cone
    rs_gup = 3.4_dp,     & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 3.0_dp,     & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 3.4_dp,      & ! real_storage for the conformation LU (HSL)
    is_c  = 3.4_dp         ! integer_storage for conformation LU (HSL)


! definition for height function (surface advection) and mesh motion

  integer, parameter :: &
    ninti_sf_adv = 3,   & ! number of Gauss points
    method = 1            ! discretization method 0: Galerkin, 1: SUPG

  type(mesh_t) :: mesh_sf_adv
  type(input_probdef_t) :: input_probdef_sf_adv
  type(problem_t) :: problem_sf_adv
  type(sysmatrix_t) :: sysmatrix_sf_adv
  type(sysvector_t), target :: sol_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1
  type(sysvector_t), target :: rhsd_sf_adv, sol_sf_adv_pred, sol_sf_adv_pred_n
  type(oldvectors_t) :: oldvectors_sf_adv
  type(subscriptvec_t) :: velr_sf, velt_sf
  type(coefficients_t) :: coefficients_sf_adv
  type(vector_t), target :: velocity_sf_adv
  real(dp), allocatable, dimension(:) :: Yhat, Yhatn, Xhat, Xhatn, displ_plate
  real(dp), allocatable, dimension(:,:) :: meshcoor_initial, meshcoor_n, &
    meshcoor_nm1, coorinitial_cone, coor_cone, &
    displ_cone, displ_cone_n
  type(problem_t), target :: problem_laplx, problem_laply
  type(subscript_t) :: subsh

! type definitions for projection of velocity gradients (if DEVSS is not used)

  type(input_probdef_t) :: input_probdef_grad
  type(problem_t) :: problem_grad
  type(sysmatrix_t) :: sysmatrix_grad
  type(sysvector_t) :: sol_grad, rhsd_grad(ncompg)
  type(oldvectors_t) :: oldvectors_grad
  type(vector_t), target :: gradients

! various locals
  character(len=30) :: filename

  integer :: icomp, step, i, k, m, npar, physqv, physqp, physqg, ipost=0
  integer :: vertices(3) = [1,3,5], obj_p
  real(dp) :: alpha, G, rpos_plate, rpos_cone, radpos, current_r, new_r
  real(dp), dimension(1) :: torque_stress
  real(dp) :: phi
  real(dp) :: Amat(3,3), Qmat(3,3), cmat(3,3), cmat2(3,3)
  logical :: surface_tension=.false., devss=.false.

! set some parameters

  if ( devss ) then
    alpha = eta_p  ! DEVSS parameter
    physqg = 1
    physqv = 2
    physqp = 3
  else
    physqv = 1
    physqp = 2
  end if

  G = eta_p / lambda  ! modulus

  if ( model == 2 ) then
    npar = 2 ! Oldroyd-B
  else
    npar = 3
  end if

  open ( unit=15, file='torque.out', status='unknown')
  open ( unit=16, file='torque_stress.out', status='unknown')

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=600, ncoefr=500+npar*nmodes )

  coefficients%i = &
    [ uintpl, pintpl, 0, 0,      gintpl, &
      physqv, physqp, 0, 0, gauss,  &
      gaussb,  cintpl, 0, 0,      0,      &
      0,      0,      model, nmodes, startm, &
      logc,   timeint1, ( 0, i = 23, 600 )  &
    ]

  coefficients%i(23) = 1  ! coorsys = 1, cilindrical coordinates, axisymmetric
  coefficients%i(40) = 3  ! numerical values for Gauss integration
  coefficients%i(48) = 1  ! use mesh velocity for ALE formulation
  coefficients%i(49) = 1  ! exp(s) projection =.true. for logc=1
  coefficients%i(61) = 0  ! projected G=0, direct velocity gradient=1 in CE
  coefficients%i(67) = 1  ! 3D velocity (swirl)

  coefficients%i(599) = 1 ! use also viscoelastic stress for torque

  coefficients%r(1:502) = &
    [ eta_s,    0._dp,   0._dp, 0._dp, 0._dp, &
      0._dp, 0._dp,  deltat, beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G, lambda  &
    ]

  if ( devss ) then
    coefficients%r(4) = alpha
    coefficients%i(9) = physqg
  else
    coefficients%i(60) = 1 ! Separate vector for velocity gradient
  end if

  if ( surface_tension ) then
    coefficients%r(19) = gamma_fa ! surface tension
  end if

  if ( model == 3 ) then
    coefficients%r(503) = alphapar ! Giesekus
  end if

! create mesh

  call create_mesh

  call fill_mesh_parts ( mesh )

! define some arrays for the ALE mesh position at old times
  allocate ( meshcoor_n(mesh%nnodes,mesh%ndim), &
    meshcoor_nm1(mesh%nnodes,mesh%ndim) )

  allocate ( meshcoor_initial(mesh%nnodes,mesh%ndim) )
    meshcoor_initial = mesh%coor

  meshcoor_n = mesh%coor(:,:)
  meshcoor_nm1 = mesh%coor(:,:)

  allocate ( coorinitial_cone(mesh%curves(2)%nnodes,mesh%ndim), &
             coor_cone(mesh%curves(2)%nnodes, mesh%ndim ) )
  allocate ( displ_cone(mesh%curves(2)%nnodes,mesh%ndim) )
  allocate ( displ_plate(mesh%curves(1)%nnodes) )
  allocate ( displ_cone_n(mesh%curves(2)%nnodes,mesh%ndim) )

  coorinitial_cone = mesh%coor(mesh%curves(2)%nodes,:)

! write mesh for plotting

  call write_mesh_vtk ( mesh, 'mesh.vtk' )

  do i = 1, mesh%ncurves
    write(filename,'(a,i4.4,a)') 'curve_',i,'.vtk'
    call write_geometry_vtk ( mesh, curve=i, filename=filename )
  end do

  call printinfo ( mesh, printlevel=2 )

! problem definition of gradient/velocity/pressure

  if ( devss ) then

    call create_input_probdef ( mesh, input_probdef, nvec=7, nphysq=3 )

    input_probdef%vec_elementdof(1)%a(:,1) = 0
    input_probdef%vec_elementdof(1)%a(vertices,1) = 6  ! gradients
    input_probdef%vec_elementdof(1)%a(:,2) = 3  ! velocity
    input_probdef%vec_elementdof(1)%a(:,3) = 0
    input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
    input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity
    input_probdef%vec_elementdof(1)%a(:,5) = 6  ! symmetric tensor
    input_probdef%vec_elementdof(1)%a(:,6) = 9  ! unsymmetric tensor
    input_probdef%vec_elementdof(1)%a(:,7) = 2  ! meshvel

    input_probdef%physq = [1,2,3]

  else

    call create_input_probdef ( mesh, input_probdef, nvec=7, nphysq=2 )

    input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
    input_probdef%vec_elementdof(1)%a(:,2) = 0
    input_probdef%vec_elementdof(1)%a(vertices,2) = 1  ! pressure
    input_probdef%vec_elementdof(1)%a(:,3) = 0
    input_probdef%vec_elementdof(1)%a(vertices,3) = 6  ! gradients
    input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity
    input_probdef%vec_elementdof(1)%a(:,5) = 6  ! symmetric tensor
    input_probdef%vec_elementdof(1)%a(:,6) = 9  ! unsymmetric tensor
    input_probdef%vec_elementdof(1)%a(:,7) = 2  ! meshvel

    input_probdef%physq = [1,2]

  end if

! Cone
  call define_essential ( mesh, input_probdef, curve1=2, physq=physqv, &
    degfd=[1,0,1] )

! Plate
  call define_essential ( mesh, input_probdef, curve1=1, physq=physqv, &
    degfd=[1,0,1])

  call define_essential ( mesh, input_probdef, point=1, physq=physqv, &
    degfd=[1,1,1])

! Rotation matrix A, for u=Au', where u are the components
! in the cylindrical system (z,r,theta) and u' in the spherical system
! (phi,rho,theta).

  Amat(1,:) = [  cos(theta_cp), sin(theta_cp), 0._dp ]
  Amat(2,:) = [ -sin(theta_cp), cos(theta_cp), 0._dp ]
  Amat(3,:) = [          0._dp,         0._dp, 1._dp ]

! full matrix
  call define_transformation ( mesh, input_probdef, curve=2, physq=physqv, &
    Amatrix=Amat )

  input_probdef%probnr = 1

  call problem_definition ( input_probdef, mesh, problem )

  call build_transformation_matrix ( mesh, problem )

! create vector subscripts for the velocity components for monitoring

  call create_subscript ( mesh, problem, velz, physqarr=[physqv], degfd=1, &
    fillnodes=.true. )
  call create_subscript ( mesh, problem, velr, physqarr=[physqv], degfd=2, &
    fillnodes=.true. )
  call create_subscript ( mesh, problem, velt, physqarr=[physqv], degfd=3, &
    fillnodes=.true. )

  call create_subscript ( mesh, problem, mvelz3, vec=7, degfd=1, &
    curves=[3], fillnodes=.true. )
  call create_subscript ( mesh, problem, mvelr3, vec=7, degfd=2, &
    curves=[3], fillnodes=.true. )

  call create_subscript ( mesh, problem, velz3, physqarr=[physqv], &
           degfd=1, curves=[3], fillnodes=.true. )
  call create_subscript ( mesh, problem, velr3, physqarr=[physqv], &
           degfd=2, curves=[3] )

! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=3, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a(:,1) = 0
  input_probdefc%vec_elementdof(1)%a(vertices,1) = 1  ! c
  input_probdefc%vec_elementdof(1)%a(:,2) = 1         ! scalar for plotting
  input_probdefc%vec_elementdof(1)%a(:,3) = 6         ! symmetric tensor

  input_probdefc%physq = [1]

  input_probdefc%probnr = 2

  call problem_definition ( input_probdefc, mesh, problemc )

! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemc, cval, physqarr=[1], fillnodes=.true. )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, solm1, sol_p, rhsd )

! create system vector ( reaction forces)
  call create_sysvector ( problem, reacf )

  call create_vector ( problem, meshvel, vec=7 )
  meshvel%u = 0._dp ! initialize

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! fill coefficients for surface advection problem

  call create_coefficients ( coefficients_sf_adv, ncoefi=100, ncoefr=50 )

  coefficients_sf_adv%i(1:4) = [ ninti_sf_adv, 2, 0, method ]
  coefficients_sf_adv%i(5:) = 0
  coefficients_sf_adv%i(5) = 1 ! time-integration
  coefficients_sf_adv%i(6) = hintpl ! height interpolation
  coefficients_sf_adv%i(9) = 3 ! numerical table for Gauss

  coefficients_sf_adv%r = 0

  coefficients_sf_adv%r(4) = deltat
  coefficients_sf_adv%r(5) = betai

! create mesh for surface advection

  call mesh_skeleton ( mesh_sf_adv, mesh%curves(3)%nnodes, &
    mesh%curves(3)%nelem, mesh%curves(3)%element%elshape, ndim=1 )

   mesh_sf_adv%topology(1)%a = mesh%curves(3)%topology(:,:,1)

! set coordinates (angles)

  do i = 1, mesh%curves(3)%nnodes
    phi = find_angle2 ( mesh%coor(mesh%curves(3)%nodes(i),1:2) )
    mesh_sf_adv%coor(i,1) = phi
  end do

  call fill_mesh_parts ( mesh_sf_adv )

! problem definition for surface advection

  call create_input_probdef ( mesh_sf_adv, input_probdef_sf_adv, nvec=3, &
    nphysq=1 )

  input_probdef_sf_adv%vec_elementdof(1)%a = &
      reshape ( [ 1,1,1,    &  ! height function
                  2,2,2,    &  ! velocity
                  3,3,3  ], &  ! meshvel
                 [3,3] )

  input_probdef_sf_adv%physq = [1]

! define problem

  input_probdef%probnr = 3

  call problem_definition ( input_probdef_sf_adv, mesh_sf_adv, problem_sf_adv )

  call create_sysvector ( problem_sf_adv, sol_sf_adv, rhsd_sf_adv )
  call create_sysvector ( problem_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1 )
  call create_sysvector ( problem_sf_adv, sol_sf_adv_pred, sol_sf_adv_pred_n )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix_sf_adv, mesh_sf_adv, &
    problem_sf_adv )

  call create_sysmatrix_data ( sysmatrix_sf_adv )

! create vectors

  call create_vector ( problem_sf_adv, velocity_sf_adv, vec=2 )

! subscripts for velocity

  call create_subscript ( mesh_sf_adv, problem_sf_adv, velt_sf, vec=2, degfd=1 )
  call create_subscript ( mesh_sf_adv, problem_sf_adv, velr_sf, vec=2, degfd=2 )
  call create ( mesh_sf_adv, problem_sf_adv, subsh, fillnodes=.true. )

  sol_sf_adv%u(subsh%s) = sqrt( mesh%coor(velz3%nodes,1)**2 + &
                                mesh%coor(velz3%nodes,2)**2 )

! fill solution vector with essential boundary conditions

  sol_sf_adv_n%u = sol_sf_adv%u
  sol_sf_adv_nm1%u =  sol_sf_adv_n%u

! allocate arrays for the ALE displacement problem

  allocate ( Xhat(mesh_sf_adv%nnodes), Xhatn(mesh_sf_adv%nnodes) )
  allocate ( Yhat(mesh_sf_adv%nnodes), Yhatn(mesh_sf_adv%nnodes) )
  Xhat = mesh%coor(mesh%curves(3)%nodes,1)
  Yhat = mesh%coor(mesh%curves(3)%nodes,2)

  call create_oldvectors ( oldvectors_sf_adv, nsysvec=2, nvec=3 )

  oldvectors_sf_adv%s(1)%p => sol_sf_adv_n    ! corrector at n
  oldvectors_sf_adv%s(2)%p => sol_sf_adv_nm1  ! corrector at nm1

  oldvectors_sf_adv%v(1)%p => velocity_sf_adv ! advection velocity at np1

! create system vectors (solution and right-hand side) for conformation and
! initialize vectors to zero stress

  call create ( problemc, solc, solcm1, solc_p, rhsc )

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
  call finalize_sysmatrix_structure ( sysmatrixc )

  call create_sysmatrix_data ( sysmatrixc )

! problem definition for projected "c=exp(s)" of the log conformation s

  call create_input_probdef ( mesh, input_probdefc_proj, nvec=1, nphysq=1 )

  input_probdefc_proj%vec_elementdof(1)%a(:,1) = 0
  input_probdefc_proj%vec_elementdof(1)%a(vertices,1) = 1  ! c

  input_probdefc_proj%physq = [1]

  input_probdef%probnr = 4

  call problem_definition ( input_probdefc_proj, mesh, problemc_proj )

  call create ( problemc_proj, solc_proj, rhsc_proj )

! initialize solc_proj
  solc_proj(1,1)%u = 1
  solc_proj(2,1)%u = 0
  solc_proj(3,1)%u = 0
  solc_proj(4,1)%u = 1
  solc_proj(5,1)%u = 0
  solc_proj(6,1)%u = 1


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec2=3, nprob=3, &
           nvec=3 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s(2)%p => solm1
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%s2(2)%p => solcm1
  oldvectors_ve%s2(3)%p => solc_proj
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc
  oldvectors_ve%p(3)%p => problemc_proj
  oldvectors_ve%v(1)%p => meshvel
  oldvectors_ve%v(3)%p => gradients

  if ( .not. devss ) then

!   create gradient vector
    call create ( problem, gradients, vec=3 )

!   define gradient projection problem
    call proj_grad_definition

  end if

! create system matrix for projection problem

  call create_sysmatrix_structure ( sysmatrixc_proj, mesh, problemc_proj, &
    symmetric=.true. )
  call create_sysmatrix_data ( sysmatrixc_proj )

! time stepping

  do step = 1, numtimesteps

    if ( step == 2 ) then
!     change time integration scheme at the second time step
      coefficients%i(22) = timeint2
      coefficients_sf_adv%i(5) = 2 ! second-order scheme
    end if

    if ( step >= 2 ) then

!     predict position of the surface and adapt mesh accordingly
      Yhatn = mesh%coor(mesh%curves(3)%nodes,2) !Yhat
      Xhatn = mesh%coor(mesh%curves(3)%nodes,1) !Xhat

      meshcoor_nm1 = meshcoor_n
      meshcoor_n = mesh%coor

      Xhat = 2 * meshcoor_n(mesh%curves(3)%nodes,1) &
                     - meshcoor_nm1(mesh%curves(3)%nodes(:),1)
      Yhat = 2 * meshcoor_n(mesh%curves(3)%nodes,2) &
                     - meshcoor_nm1(mesh%curves(3)%nodes(:),2)

      current_r = sqrt( mesh%coor(mesh%points(3),1)**2 + &
                        mesh%coor(mesh%points(3),2)**2 )

      new_r = 2 * current_r - sqrt( meshcoor_nm1(mesh%points(3),1)**2 + &
                                    meshcoor_nm1(mesh%points(3),2)**2 )


      do i = 1,mesh%curves(2)%nnodes

!       r position of current node
        radpos = sqrt( mesh%coor(mesh%curves(2)%nodes(i),1)**2 + &
                       mesh%coor(mesh%curves(2)%nodes(i),2)**2 )

!       scale proportionally to point 3
        rpos_cone = new_r/current_r * radpos

        displ_cone(i,1:2) = rpos_cone * [ sin(theta_cp), cos(theta_cp) ] &
                              - mesh%coor(mesh%curves(2)%nodes(i),1:2)

      end do

      current_r = mesh%coor(mesh%points(2),2)

      new_r = 2 * current_r - meshcoor_nm1(mesh%points(2),2)


      do i = 1,mesh%curves(1)%nnodes

!       r position of current node
        radpos = mesh%coor(mesh%curves(1)%nodes(i),2)

!       scale proportionally to point 3
        rpos_plate = new_r/current_r * radpos

        displ_plate(i) = rpos_plate - mesh%coor(mesh%curves(1)%nodes(i),2)

      end do


      call update_mesh_nodes ( mesh, problem_laplx, problem_laply, &
        disp1=Xhat-Xhatn, disp2=Yhat-Yhatn, disp3=displ_cone(:,1), &
        disp4=displ_cone(:,2), disp5=displ_plate )

      call find_bounds_blocks ( mesh )

!     mesh velocity

      meshvel%u = reshape ( transpose ( &
         (1.5_dp*mesh%coor - 2*meshcoor_n + 0.5_dp*meshcoor_nm1 ) /deltat ), &
                             [2*mesh%nnodes] )

    end if

!   fill solution vector with essential boundary conditions for this step
    sol%u = 0._dp
    call fill_sysvector ( mesh, problem, sol, &
    curve1=2, physq=physqv, degfd=3, func=func, funcnr=1 )

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG

!   exps projection (for log conformation)

    if ( logc == 1 ) call solve_exps_projection

!   build implicit terms of CE with rhs in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve, physqrow=[physqv], physqcol=[physqv], &
      addmatvec=.true., coefficients=coefficients )

!   surface tension

    if ( surface_tension ) then

!     build surface integral on free surface

      coefficients%r(19) = gamma_fa

      call add_boundary_elements ( mesh, problem, rhsd, &
        elemsub=surface_tension_curve, &
        curve=3, coefficients=coefficients, physq=[physqv] )

!     surface tension between fluid and walls

      coefficients%r(19) = dgamma_wa
      call add_boundary_elements ( mesh, problem, rhsd, &
        elemsub=surface_tension_curve, &
        curve=1, coefficients=coefficients, physq=[physqv] )
      call add_boundary_elements ( mesh, problem, rhsd, &
        elemsub=surface_tension_curve, &
        curve=2, coefficients=coefficients, physq=[physqv] )

    end if


    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve gradient/velocity/pressure problem
!m    if (step == 1 ) then
!m       solver_options_u%pivot_order = 1
!m       solver_options_u%scaling = 1
!m       call renumber_metis_sysmatrix(sysmatrix)
!m    end if

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

    call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

    call transform_to_global ( problem, sol, reacf )

    if ( .not. devss ) call build_and_solve_proj_grad

    call copy ( sol, solm1 )

!   build (assemble) matrix and vector for conformation problem

    if ( coefficients%i(22) == timeint1 ) then

      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    else

      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem_implicit_2nd_order, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    end if

    call check ( sysmatrixc )

    call copy ( solc, solcm1 )


!   solve conformation and keep LU decomposition in the loop over components
!m    if (step == 1 ) then
!m       solver_options_c%pivot_order = 1
!m       solver_options_c%scaling = 1
!m       call renumber_metis_sysmatrix(sysmatrixc)
!m    end if

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do icomp = 1, ncompc
      call solve_system_ma41 ( sysmatrixc, rhsc(icomp,1), solc(icomp,1), luc, &
        solver_options=solver_options_c  )
    end do

    call delete ( luc )  ! remove LU decomposition and rebuild next time step

    call solve_surface_height_corrector

!   compute torque

    call calc_torque

    if ( step == numtimesteps ) then

      allocate ( sample_coor ( 51, mesh%ndim) )
      sample_coor(:,2) = [(0._dp+(i-1)*((Ro)/50), &
                      i=1,(50+1))]
      sample_coor(:,1) = 0.0_dp

      warn_add_to_mesh_after_meshgen_parts = .false.
      call add_to_mesh ( mesh, object='coordinates', nnodes=51, &
         coor=sample_coor )
      obj_p = mesh%nobjects
      call fill_mesh_parts_objects ( mesh, object1=obj_p )
      warn_add_to_mesh_after_meshgen_parts = .true.

      call fill_sample(mesh, problem, sample_p, object=obj_p, ndegfd=1, &
             elemsub=stokes_sample_pressure, coefficients=coefficients, &
             oldvectors=oldvectors_ve )
      open ( unit=500, file='convergence_p.out' )
      do i = 1, sample_p%nnodes
         write ( 500, '(2es16.8)') sample_p%u(i,1), sample_p%coor(i,2)
      end do
      close(unit=500)
      call delete ( sample_p )
    end if

  end do

! post-processing

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )
  oldvectors%s(1)%p => sol

  call create_vector ( problem, velocity, physq=physqv )
  call create_vector ( problem, velocity_rot, physq=physqv )
  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, gammadot, vec=4 )
  call create_vector ( problem, Dtensor, vec=5 )
  call create_vector ( problem, Ltensor, vec=6 )
  call create_vector ( problem, tauviscous, vec=5 )

  call extract_physvector ( mesh, problem, sol, velocity )

! derive the components of the velocity in a (phi,rho,theta) spherical
! coordinate system

  do i = 1, size(velz%s)

!   get angle of position vector of current node w.r.t. the y-axis
    phi = find_angle2([mesh%coor(velz%nodes(i),1),mesh%coor(velz%nodes(i),2)])

!   rotation matrix Q, for u'=Qu, where u'are the components
!   in the spherical system (phi,rho,theta) and u in the cylindrical system
!   (z,r,theta).
    Qmat(1,:) = [ cos(phi), -sin(phi), 0._dp ]
    Qmat(2,:) = [ sin(phi),  cos(phi), 0._dp ]
    Qmat(3,:) = [    0._dp,     0._dp, 1._dp ]

!   calculate velocity components in the spherical coordinate system
    sol_p%u([velz%s(i), velr%s(i), velt%s(i)]) = matmul ( Qmat, &
                     [sol%u(velz%s(i)), sol%u(velr%s(i)), sol%u(velt%s(i))] )

  end do

  call extract_physvector ( mesh, problem, sol_p, velocity_rot )

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
    filename='cp4.vtk' )

  call write_vector_vtk ( mesh, problem, filename='cp4.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true., &
    assume3D=.true. )

  call write_vector_vtk ( mesh, problem, filename='cp4.vtk', &
    dataname='velocity_vector_planar', vector=velocity, append=.true., &
    degfd=[1,2] )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='cp4.vtk', dataname='gammadot', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='cp4.vtk', &
    dataname='D', vector=Dtensor, append=.true., assume3D=.true. )

  call write_tensor_vtk ( mesh, problem, filename='cp4.vtk', &
    dataname='L', vector=Ltensor, append=.true., assume3D=.true., &
    symmetric=.false. )

  call write_tensor_vtk ( mesh, problem, filename='cp4.vtk', &
    dataname='tauviscous', vector=tauviscous, append=.true., assume3D=.true. )

  call write_vector_vtk ( mesh, problem, filename='cp4.vtk', &
    dataname='velocity_rot', vector=velocity_rot, append=.true., &
    assume3D=.true. )

  call create_vector ( problemc, tauviscoelastic, vec=3 )

  call derive_vector ( mesh, problemc, tauviscoelastic, &
    elemsub=deriv_viscoelastic_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemc, filename='cp4.vtk', &
    dataname='tauviscoelastic', vector=tauviscoelastic, append=.true., &
    assume3D=.true. )

  call create_vector ( problemc, ctensor, vec=3 )
  call create_vector ( problemc, ctensor_rot, vec=3 )

! loop over the viscoelastic modes
  do m = 1, nmodes

!   derive the components of the conformation in a (phi,rho,theta) spherical
!   coordinate system
    do i = 1, size(cval%s)

!     get angle of position vector of current node w.r.t. the y-axis
      phi = find_angle2([mesh%coor(cval%nodes(i),1),mesh%coor(cval%nodes(i),2)])

!     rotation matrix Q, for u'=Qu, where u' and u are the components
!     in the spherical system (phi,rho,theta) and the cylindrical system
!     (z,r,theta).
      Qmat(1,:) = [ cos(phi), -sin(phi), 0._dp ]
      Qmat(2,:) = [ sin(phi),  cos(phi), 0._dp ]
      Qmat(3,:) = [    0._dp,     0._dp, 1._dp ]

!     get the orginal tensor components
      cmat = reshape ([solc(1,m)%u(cval%s(i)), &
                       solc(2,m)%u(cval%s(i)), &
                       solc(3,m)%u(cval%s(i)), &
                       solc(2,m)%u(cval%s(i)), &
                       solc(4,m)%u(cval%s(i)), &
                       solc(5,m)%u(cval%s(i)), &
                       solc(3,m)%u(cval%s(i)), &
                       solc(5,m)%u(cval%s(i)), &
                       solc(6,m)%u(cval%s(i)) ], [3,3] )

!     calculate tensor components in the spherical coordinate system
      cmat2 = matmul ( Qmat, matmul ( cmat, transpose ( Qmat ) ) )

!     put components in solution vector (symmetry is automatically satisfied!)
      solc_p(1,m)%u(cval%s(i)) = cmat2(1,1)   ! phi,phi
      solc_p(2,m)%u(cval%s(i)) = cmat2(1,2)   ! phi,rho
      solc_p(3,m)%u(cval%s(i)) = cmat2(1,3)   ! phi,theta
      solc_p(4,m)%u(cval%s(i)) = cmat2(2,2)   ! rho,rho
      solc_p(5,m)%u(cval%s(i)) = cmat2(2,3)   ! rho,theta
      solc_p(6,m)%u(cval%s(i)) = cmat2(3,3)   ! theta,theta

    end do

  end do

  call derive_vector ( mesh, problemc, ctensor, &
    elemsub=deriv_conformation_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  oldvectors_ve%s2(1)%p => solc_p

  call derive_vector ( mesh, problemc, ctensor_rot, &
    elemsub=deriv_conformation_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemc, filename='cp4.vtk', &
    dataname='c', vector=ctensor, append=.true., &
    assume3D=.true. )

  call write_tensor_vtk ( mesh, problemc, filename='cp4.vtk', &
    dataname='c_rot', vector=ctensor_rot, append=.true., &
    assume3D=.true. )

  call create ( mesh, problem, szz, vec=5, degfd=1, curves=[1], &
    fillnodes=.true. )
  call create ( mesh, problem, szr, vec=5, degfd=2, curves=[1] )
  call create ( mesh, problem, szt, vec=5, degfd=3, curves=[1] )
  call create ( mesh, problem, srr, vec=5, degfd=4, curves=[1] )
  call create ( mesh, problem, srt, vec=5, degfd=5, curves=[1] )
  call create ( mesh, problem, stt, vec=5, degfd=6, curves=[1] )

  open( unit=84, file='stress_on_plate4.out')

  do i = 1, size(szz%s)
    k = szz%nodes(i)
    write(84,'(15es16.8)') mesh%coor(k,1), mesh%coor(k,2), pressure%u(k), &
      tauviscous%u(szz%s(i)), tauviscous%u(szr%s(i)), &
      tauviscous%u(szt%s(i)), tauviscous%u(srr%s(i)), &
      tauviscous%u(srt%s(i)), tauviscous%u(stt%s(i)), &
      tauviscoelastic%u(szz%s(i)), tauviscoelastic%u(szr%s(i)), &
      tauviscoelastic%u(szt%s(i)), tauviscoelastic%u(srr%s(i)), &
      tauviscoelastic%u(srt%s(i)), tauviscoelastic%u(stt%s(i))
  end do

  close(unit=84)

  call create ( mesh, problem, szz, vec=5, degfd=1, curves=[3], &
    fillnodes=.true. )
  call create ( mesh, problem, szr, vec=5, degfd=2, curves=[3] )
  call create ( mesh, problem, szt, vec=5, degfd=3, curves=[3] )
  call create ( mesh, problem, srr, vec=5, degfd=4, curves=[3] )
  call create ( mesh, problem, srt, vec=5, degfd=5, curves=[3] )
  call create ( mesh, problem, stt, vec=5, degfd=6, curves=[3] )

  open( unit=85, file='stress_on_surface4.out')

  do i = 1, size(szz%s)
     k = szz%nodes(i)
     write(85,'(15es16.8)') mesh%coor(k,1), mesh%coor(k,2), pressure%u(k), &
       tauviscous%u(szz%s(i)), tauviscous%u(szr%s(i)), &
       tauviscous%u(szt%s(i)), tauviscous%u(srr%s(i)), &
       tauviscous%u(srt%s(i)), tauviscous%u(stt%s(i)), &
       tauviscoelastic%u(szz%s(i)), tauviscoelastic%u(szr%s(i)), &
       tauviscoelastic%u(szt%s(i)), tauviscoelastic%u(srr%s(i)), &
       tauviscoelastic%u(srt%s(i)), tauviscoelastic%u(stt%s(i))
  end do

  close(unit=85)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, solm1, rhsd, sol_p )
  call delete ( sysmatrix )
  call delete ( velocity, pressure, gammadot, Dtensor, Ltensor, tauviscous, &
                tauviscoelastic, ctensor, velocity_rot, ctensor_rot )
  call delete ( oldvectors, oldvectors_ve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( solc, solcm1, rhsc, solc_p )
  call delete ( reacf )

  call delete ( problemc_proj )
  call delete ( input_probdefc_proj )
  call delete ( solc_proj, rhsc_proj )

  call delete ( coefficients )
  call delete ( velz, velr, velt, cval )
  call delete ( szz, szr, szt, srr, srt, stt )

  if ( .not. devss ) then
    call delete ( gradients )
    call delete ( input_probdef_grad )
    call delete ( problem_grad )
    call delete ( sysmatrix_grad )
    call delete ( sol_grad )
    call delete ( rhsd_grad )
    call delete ( oldvectors_grad )
  end if


contains

! create the mesh for the cone-plate geometry using gmsh

  subroutine create_mesh

    open ( unit=25, file='ConePlate2d.geo' )

    write ( 25, '(1X,A,F18.14,A)' ) 'z[1] = ', Ro*cos(pi/2-theta_cp), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'r[1] = ', Ro*sin(pi/2-theta_cp), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'theta_cp = ', theta_cp, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_Ro = ', dx_Ro, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_0 = ', dx_0, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'Ro = ', Ro, ';'
    write ( 25, '(/1x,a)' ) 'Include "ConePlate2d.igo";'

    close ( 25 )

    call execute_command_line ( 'gmsh -2 -order 2 &
      & -o ConePlate2d.msh ConePlate2d.geo > outputmesh.out' )

!   read mesh generated by gmsh
    call read_mesh_gmsh ( mesh, filename='ConePlate2d.msh', ndim=2, &
      physgeom=.true., sortphys=.true. )

  end subroutine create_mesh

  subroutine build_vpG

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, physqrow=[physqv,physqp], &
      physqcol=[physqv,physqp], &
      coefficients=coefficients )

    if ( devss ) then

!     DEVSS-G
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=devssg_elem, addmatvec=.true., &
        physqrow=[physqg,physqv], physqcol=[physqg,physqv], &
        coefficients=coefficients )

!     set to zero off-diagonal blocks gradient-pressure
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        buildvector=.false., physqrow=[physqg], physqcol=[physqp], &
        zeromatvec=.true. )
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        buildvector=.false., physqrow=[physqp], physqcol=[physqg], &
        zeromatvec=.true. )

    end if

  end subroutine build_vpG


! project c=exp(s) on discrete fem space

  subroutine solve_exps_projection

    type(solver_options_ma57_t) :: solver_options_ma57

    integer :: i, m

    call build_system ( mesh, problemc_proj, sysmatrixc_proj, &
      m2sysvector=rhsc_proj, elemsub=exps_projection_elem, &
      oldvectors=oldvectors_ve, coefficients=coefficients )

    ! MA57 solver storage
    solver_options_ma57%integer_storage = 1.3
    solver_options_ma57%real_storage    = 1.3

!   LU decomposition is done in the first call only

    do m = 1, nmodes
      do i = 1, ncompc

        call add_effect_of_essential_to_rhs ( problemc_proj, sysmatrixc_proj, &
           solc_proj(i,m), rhsc_proj(i,m) )

        call solve_system_ma57 ( sysmatrixc_proj, rhsc_proj(i,m), &
           solc_proj(i,m), lu_exps_proj, solver_options=solver_options_ma57 )

      end do
    end do

    call delete ( lu_exps_proj )

  end subroutine solve_exps_projection

! function for velocity on wall

  function func ( nr, x )

    integer, intent(in) ::  nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    select case(nr)
      case(1)
        func = omega*x(2)
      case default
        write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
        stop
    end select

  end function func

  subroutine calc_torque

    real(dp) :: torque1, radius
    type(subscript_t) :: s_force1
    type(vector_t) :: tractionr1

    integer :: i

    call create_vector ( problem, tractionr1, vec=4)
    call create_subscript ( mesh, problem, s_force1, physqarr=[physqv], &
                            degfd=3, curves=[2], fillnodes=.true.)

!   initialize
    tractionr1%u = 0._dp

    do i = 1,size(s_force1%s)
      k = s_force1%nodes(i)
      radius = mesh%coor(k,2)
      tractionr1%u(i) = reacf%u(s_force1%s(i))*radius
    end do

    torque1 = sum(tractionr1%u(:))

    call calc_torque_from_stress
    print *, torque1, torque_stress

    write(unit=15,fmt='(4(es19.9))') torque1, torque_stress

    call delete (s_force1 )

  end subroutine calc_torque


  subroutine calc_torque_from_stress

    integer :: i

    type(vector_t), target :: pressure1, stokes_tensor1, radius_vec
    type(vector_t), dimension(nmodes), target :: solc_modes1
    type(oldvectors_t) :: oldvectors_torque_stress1
    real(dp), dimension(1) :: torque_stress1

    do i = 1,nmodes
      call create_vector ( problemc, solc_modes1(i), vec=3 )
    end do

    call create_vector ( problem, pressure1, vec=4 )
    call create_vector ( problem, radius_vec, vec=4 )
    call create_vector ( problem, stokes_tensor1, vec=5 )

    do i =1,mesh%nnodes
      radius_vec%u(i) = mesh%coor(i,2)
    end do

!   create the structure oldvectors

    call create_oldvectors ( oldvectors_torque_stress1, nsysvec=1, &
      nsysvec2=1, nprob=2, nvec=3, nvec1=1 )

    oldvectors_torque_stress1%s(1)%p => sol
    oldvectors_torque_stress1%s2(1)%p => solc
    oldvectors_torque_stress1%p(1)%p => problem
    oldvectors_torque_stress1%p(2)%p => problemc
    oldvectors_torque_stress1%v1(1)%p => solc_modes1
    oldvectors_torque_stress1%v(1)%p => stokes_tensor1
    oldvectors_torque_stress1%v(2)%p => pressure1
    oldvectors_torque_stress1%v(3)%p => radius_vec

    torque_stress1 = 0._dp

    do i = 1,nmodes
      coefficients%i(28) = i
      call derive_vector ( mesh, problemc, solc_modes1(i), &
        elemsub=deriv_conformation_tensor_std, elcurves=[2], &
        coefficients=coefficients, oldvectors=oldvectors_torque_stress1 )
    end do

    coefficients%i(28) = 0

    call derive_vector ( mesh, problem, stokes_tensor1, &
      elemsub=stokes_stress_tensor, elcurves=[2], &
      coefficients=coefficients, oldvectors=oldvectors_torque_stress1 )

    call derive_vector ( mesh, problem, pressure1, elemsub=stokes_pressure, &
      elcurves=[2], coefficients=coefficients, &
      oldvectors=oldvectors_torque_stress1 )

    call integrate_boundary_elements ( mesh, problem, torque_stress1, &
      elemsub=torque_from_stress, curve=2, coefficients=coefficients, &
      oldvectors=oldvectors_torque_stress1 )

    call delete (oldvectors_torque_stress1)

    write(unit=16,fmt='(4(es19.9))') torque_stress1

    torque_stress = torque_stress1

  end subroutine calc_torque_from_stress


! solve convection equation for the surface height (corrector)

  subroutine solve_surface_height_corrector

    use postprocessing_m

    integer :: i
    type(solver_options_ma41_t) :: solver_options_h

!   rotate the (u_z,u_r) velocity vector components to (u_phi,u_rho) components
!   in the polar coordinate system

    do i = 1,size(velz3%s)

      k = velz3%nodes(i)

      phi =  find_angle2(mesh%coor(k,1:2))

      velocity_sf_adv%u(velt_sf%s(i)) = &
            sol%u(velz3%s(i))*cos(phi) - sin(phi)*sol%u(velr3%s(i))
      velocity_sf_adv%u(velr_sf%s(i)) = &
            sol%u(velz3%s(i))*sin(phi) + cos(phi)*sol%u(velr3%s(i))

    end do

    call build_system ( mesh_sf_adv, problem_sf_adv, sysmatrix_sf_adv, &
      rhsd_sf_adv, elemsub=surface_advection_polar_elem, &
      oldvectors=oldvectors_sf_adv, coefficients=coefficients_sf_adv )

    call check_filled_sysmatrix ( sysmatrix_sf_adv )

    call add_effect_of_essential_to_rhs ( problem_sf_adv, sysmatrix_sf_adv, &
      sol_sf_adv, rhsd_sf_adv )

!   MA41 solver storage

    solver_options_h%integer_storage = 5.0
    solver_options_h%real_storage    = 5.0

!   solve system

    call solve_system_ma41 ( sysmatrix_sf_adv, rhsd_sf_adv, sol_sf_adv, &
                             solver_options=solver_options_h )

    do i = 1,size(subsh%s)

      k = subsh%nodes(i)
      mesh%coor(mesh%curves(3)%nodes(i),1) = &
                sol_sf_adv%u(subsh%s(i))*sin(mesh_sf_adv%coor(k,1))
      mesh%coor(mesh%curves(3)%nodes(i),2) = &
                sol_sf_adv%u(subsh%s(i))*cos(mesh_sf_adv%coor(k,1))

    end do


    if ( step == 1 .or. mod(step,vtkevery) == 0 ) then
      write (filename,'(a,i4.4,a)')'meshintime', ipost, '.vtk'
      call write_mesh_vtk ( mesh, filename=filename )
      ipost = ipost + 1
    end if

    call copy ( sol_sf_adv_n, sol_sf_adv_nm1 )
    call copy ( sol_sf_adv, sol_sf_adv_n )

    print *, step, step*deltat, maxval(sol_sf_adv%u)

  end subroutine solve_surface_height_corrector


! find the angle (0 <= angle < 2*pi) of the position vector with positive y-axis
! NOTE: points on the positive y-axis (to precision eps) have angle zero

  function find_angle2 ( x )

    real(dp), parameter :: eps=1e-14_dp

    real(dp), intent(in) :: x(:)
    real(dp) :: find_angle2

    if ( all(abs(x) < eps) ) then
      find_angle2 = 0._dp
    end if

    if ( abs(x(1)) < eps ) then

      find_angle2 = 0._dp

    else

      find_angle2 = atan2 ( x(2), x(1) )

      if ( x(2) < 0._dp ) then ! bottom half-plane
        find_angle2 = pi/2 + abs(find_angle2)
      else
        if ( x(1) > 0._dp ) then ! right half-plane
          find_angle2 = pi/2 - abs(find_angle2)
        else
          find_angle2 = 3*pi/2+(pi-abs(find_angle2))
        end if
      end if

     end if

   end function find_angle2

! project the velocity gradients

  subroutine build_and_solve_proj_grad

    integer :: i
    type(lu_ma41_t) :: lu_g_ma41
    type(solver_options_ma41_t) :: solver_options_ma41

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem_grad, sysmatrix_grad, &
      msysvector=rhsd_grad, elemsub=gradient_from_velocity_elem, &
      coefficients=coefficients, oldvectors=oldvectors_grad )

!m    if (step == 1 .or. remeshing) then
!m      solver_options_ma41%pivot_order = 1
!m      solver_options_ma41%scaling = 1
!m      call renumber_metis_sysmatrix(sysmatrix_grad)
!m    end if

!   MA57 solver storage
    solver_options_ma41%integer_storage = 2.3
    solver_options_ma41%real_storage    = 2.3


!   solve

    do i = 1, ncompg

      call add_effect_of_essential_to_rhs ( problem_grad, sysmatrix_grad, &
        sol_grad, rhsd_grad(i) )

        call solve_system_ma41 ( sysmatrix_grad, rhsd_grad(i), sol_grad, &
          lu=lu_g_ma41, solver_options=solver_options_ma41 )

      call transfer_data ( mesh, problem_grad, problem, &
        sysvector1=sol_grad, vector2=gradients, degfd2=[i] )

    end do

    call delete ( lu_g_ma41 )

  end subroutine build_and_solve_proj_grad

! define problem for the projection of the velocity gradients

  subroutine proj_grad_definition

    integer, dimension(3), parameter :: vertices=[1,3,5]

!   problem definition for gradients

    call create_input_probdef ( mesh, input_probdef_grad, nvec=2, nphysq=1 )

    input_probdef_grad%vec_elementdof(1)%a(:,1) = 0
    input_probdef_grad%vec_elementdof(1)%a(vertices,1) = 1  ! gradient component

    input_probdef_grad%physq = [1]
    input_probdef_grad%probnr = 5

    call problem_definition ( input_probdef_grad, mesh, problem_grad )

!   create system vectors (solution and right-hand side)

    call create ( problem_grad, sol_grad )
    call create ( problem_grad, rhsd_grad )

!   create system matrix

    call create_sysmatrix_structure ( sysmatrix_grad, mesh, problem_grad )

    call create_sysmatrix_data ( sysmatrix_grad )

!   oldvectors

    call create_oldvectors ( oldvectors_grad, nsysvec=1, nprob=1 )

    oldvectors_grad%s(1)%p => sol
    oldvectors_grad%p(1)%p => problem

  end subroutine proj_grad_definition

end program coneplate4
