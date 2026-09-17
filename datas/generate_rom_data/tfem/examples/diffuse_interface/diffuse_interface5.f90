! diffuse-interface problem on a unit square
! linear Couette flow
! constant viscosity
! first-order time integration
! Lees Edwards boundary conditions

module subs5_m

  use tfem_elem_m
  implicit none

  real(dp), save :: Udiff = 0  ! velocity difference between upper and lower

contains


! weak connection element for two objects (and two degrees of freedom)

  subroutine elementc ( mesh, problem, constr, elem, node, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemmat2, elemmatadd, &
    elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer, parameter :: ndf = 9, nodalpb = 3, ndfb=3, ndflb = 2

    integer :: object, object2, i, j
    real(dp) :: phi1(1,ndf), phi2(1,ndf), xr(1,2), psi(1,ndflb), x(nodalpb,2)
    real(dp) :: theta(1,ndfb), dtheta(1,ndfb,1), dxdxi(1,2), curvel(1)

    object = problem%constraints(constr)%object
    object2 = problem%constraints(constr)%object2

!   reference coordinates and shape function in object

    xr(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)
    call shape_quad_Q2 ( xr, phi1 )

!   reference coordinates and shape function in object2

    xr(1,:) = mesh%objects(object2)%refcoor_int(node,:,elem)
    call shape_quad_Q2 ( xr, phi2 )

!   shape function of the Lagrange multiplier

    call shape_line_P1 ( mesh%objects(object)%xig(node:node,1), psi )

!   shape function of the quadratic curve in the integration point

    call shape_line_P2 ( mesh%objects(object)%xig(node:node,1), theta, &
      dtheta(:,:,1) )

!   compute geometry of element deformed

    call get_coordinates_object ( mesh, elem, x, object )

    call isoparametric_deformation_curve ( x, dtheta(:,:,1), dxdxi, curvel )

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

      elemmat = 0
      elemmat2 = 0

!     diagonal block for first dof

      do i = 1, ndflb
        do j = 1, ndf
          elemmat(i,j) = psi(1,i) * phi1(1,j) * &
                                   curvel(1)*mesh%objects(object)%wg(node)
          elemmat2(i,j) = - psi(1,i) * phi2(1,j) * &
                                   curvel(1)*mesh%objects(object)%wg(node)
        end do
      end do

!     fill diagonal block for second dof

      elemmat(ndflb+1:2*ndflb,ndf+1:2*ndf) = elemmat(1:ndflb,1:ndf)
      elemmat2(ndflb+1:2*ndflb,ndf+1:2*ndf) = elemmat2(1:ndflb,1:ndf)

    end if

  end subroutine elementc


! weak connection element for two objects (and two degrees of freedom)
! velocity difference between objects in right-hand side

  subroutine elementc2 ( mesh, problem, constr, elem, node, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemmat2, elemmatadd, &
    elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer, parameter :: ndf = 9, nodalpb = 3, ndfb=3, ndflb = 2

    integer :: object, object2, i, j
    real(dp) :: phi1(1,ndf), phi2(1,ndf), xr(1,2), psi(1,ndflb), x(nodalpb,2)
    real(dp) :: theta(1,ndfb), dtheta(1,ndfb,1), dxdxi(1,2), curvel(1)

    object = problem%constraints(constr)%object
    object2 = problem%constraints(constr)%object2

!   reference coordinates and shape function in object

    xr(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)
    call shape_quad_Q2 ( xr, phi1 )

!   reference coordinates and shape function in object2

    xr(1,:) = mesh%objects(object2)%refcoor_int(node,:,elem)
    call shape_quad_Q2 ( xr, phi2 )

!   shape function of the Lagrange multiplier

    call shape_line_P1 ( mesh%objects(object)%xig(node:node,1), psi )

!   shape function of the quadratic curve in the integration point

    call shape_line_P2 ( mesh%objects(object)%xig(node:node,1), theta, &
      dtheta(:,:,1) )

!   compute geometry of element deformed

    call get_coordinates_object ( mesh, elem, x, object )

    call isoparametric_deformation_curve ( x, dtheta(:,:,1), dxdxi, curvel )

    if ( vector ) then

!     rhs for velocity difference between upper and lower "wall"

      elemvec = 0
      do i = 1, ndflb
        elemvec(i) = - psi(1,i) * Udiff * &
                                 curvel(1)*mesh%objects(object)%wg(node)
      end do

    end if

    if ( matrix ) then

      elemmat = 0
      elemmat2 = 0

!     diagonal block for first dof

      do i = 1, ndflb
        do j = 1, ndf
          elemmat(i,j) = psi(1,i) * phi1(1,j) * &
                                   curvel(1)*mesh%objects(object)%wg(node)
          elemmat2(i,j) = - psi(1,i) * phi2(1,j) * &
                                   curvel(1)*mesh%objects(object)%wg(node)
        end do
      end do

!     fill diagonal block for second dof

      elemmat(ndflb+1:2*ndflb,ndf+1:2*ndf) = elemmat(1:ndflb,1:ndf)
      elemmat2(ndflb+1:2*ndflb,ndf+1:2*ndf) = elemmat2(1:ndflb,1:ndf)

    end if

  end subroutine elementc2

end module subs5_m

program diffuse_interface5

  use tfem_m
  use diffuse_interface_elements_m
  use diffuse_interface_functions_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m
  use subs5_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities and (c,mu)
    pintpl = 4,     & ! Q1 pressures
    physqvel = 1,   & ! physical quantity nr of the velocities
    physqpress = 2, & ! physical quantity nr of the pressures
    gauss = 3         ! 3x3 integration of quads


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefdi
  type(problem_t), target :: problem, problemdi
  type(sysmatrix_t) :: sysmatrix, sysmatrixdi
  type(sysvector_t), target :: sol, soldi, soldin
  type(sysvector_t) :: rhsd, rhsdi, soldi_prev
  type(oldvectors_t) :: oldvectors_di
  type(coefficients_t) :: coefficients
  type(subscript_t) :: cval, muval
  type(solver_options_ma57_t) :: solver_options_up
  type(solver_options_ma41_t) :: solver_options_di



! variables

  integer :: &
    nx=20,               & ! number of elements in x
    ny=20,               & ! number of elements in y
    numtimesteps = 2,    & ! number of time steps
    intrule = 2,         & ! integration rule for objects
    nsubint = 5,         & ! number of subintervals in integration
    itermax = 10           ! maximum number of Picard iterations

  real(dp) :: &
    eta = 0.1_dp,      & ! viscosity
    deltat = 1.e-2_dp, & ! time step
    rho = 1.0_dp,      & ! density in diffuse-interface method
    alpha = 1.0_dp,    & ! parameter in the diffuse-interface method
    beta = 1.0_dp,     & ! parameter in the diffuse-interface method
    Mcoef = 1.0_dp,    & ! parameter in the diffuse-interface method
    kappa = 1.0_dp,    & ! parameter in the diffuse-interface method
    U = 1.0_dp,        & ! velocity difference upper-lower wall
    epsdi = 1.e-3_dp,  & ! accuracy in Picard iteration
    rs_up = 1.0_dp,    & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.0_dp,    & ! integer_storage velocity-pressure LU (HSL)
    rs_di  = 1.0_dp,   & ! real_storage for the diffuse-interface LU (HSL)
    is_di  = 1.4_dp      ! integer_storage for diffuse-interface LU (HSL)

  integer ::  step, iter
  real(dp) :: cmax, cdiff, mumax, mudiff, point(1,2)


! namelist for input of variables; read from standard input

  namelist /comppar/ nx, ny, numtimesteps, rho, eta, deltat, &
    alpha, beta, Mcoef, kappa, U, epsdi, itermax, rs_up, is_up, rs_di, is_di, &
    intrule, nsubint

  read ( unit=*, nml=comppar )


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=200, ncoefr=150 )

! Stokes
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

! diffuse interface
  coefficients%r(101:106) = [ deltat, Mcoef, alpha, beta, kappa, rho ]


! create mesh

  meshgen_options%elshape = 6 ! 9-node quads
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5
  point(1,:) = [0._dp,0.5_dp]
  call add_to_mesh ( mesh, points=point )  ! point 5

! define objects for the periodical boundary conditions

! left boundary
  call add_to_mesh ( mesh, object='curve', objectcurve=4, topology=.true., &
    intrule=gauss )
! right boundary
  call add_to_mesh ( mesh, object='curve', objectcurve=4, topology=.true., &
    intrule=gauss )
  mesh%objects(2)%coor(:,1) = mesh%objects(2)%coor(:,1) + 1 ! shift right

! lower boundary
  call add_to_mesh ( mesh, object='curve', objectcurve=1, topology=.true., &
    intrule=intrule, nsubint=nsubint )
! upper boundary
  call add_to_mesh ( mesh, object='curve', objectcurve=1, topology=.true., &
    intrule=intrule, nsubint=nsubint )
  mesh%objects(4)%coor(:,2) = mesh%objects(4)%coor(:,2) + 1 ! shift up

  call fill_mesh_parts ( mesh )


! problem definition of velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

  input_probdef%probnr = 1

! Dirichlet boundary conditions

! velocity level in x and y
  call define_essential ( mesh, input_probdef, point=5, physq=physqvel )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! constraints for periodical boundary conditions

! velocities x-x
  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, object=1, object2=2, discretization='weak', &
    elementdof=[2,0,2] )
!  call define_constraint ( mesh, input_probdef, &
!    physq=physqvel, curve1=2, curve2=5, discretization='collocation' )!, &
!    exclude=1 )

! velocities vertically (Lees-Edwards)
  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, object=3, object2=4, discretization='weak', &
    elementdof=[2,0,2] )

  call problem_definition ( input_probdef, mesh, problem )


! problem definition diffuse interface

  call create_input_probdef ( mesh, input_probdefdi, nvec=3, nphysq=2 )

  input_probdefdi%vec_elementdof(1)%a =   &
      reshape ( [ 1,1,1,1,1,1,1,1,1,    &  ! c
                  1,1,1,1,1,1,1,1,1,    &  ! mu
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar for plotting
                  [9,3] )

  input_probdefdi%physq = [1,2]
  input_probdefdi%probnr = 2

! constraint for periodical boundary conditions of the (c,mu)
  call define_constraint ( mesh, input_probdefdi, &
    object=1, object2=2, discretization='weak', elementdof=[2,0,2] )
!  call define_constraint ( mesh, input_probdefdi, curve1=2, curve2=5, &
!    discretization='collocation' )

! vertically (Lees-Edwards)
  call define_constraint ( mesh, input_probdefdi, &
    object=3, object2=4, discretization='weak', elementdof=[2,0,2] )

  call problem_definition ( input_probdefdi, mesh, problemdi )


! create a vector subscript for the concentration and mu without Lagr. multipl.

  call create_subscript ( mesh, problemdi, cval, physqarr=[1] )
  call create_subscript ( mesh, problemdi, muval, physqarr=[2] )


! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )


! fill solution vector with essential boundary conditions

  sol%u = 0

! create system matrix of velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! build (assemble) matrix and vector for velocity/pressure problem

! stokes velocity/pressure
  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

! periodical condition on velocities x-x

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=elementc, addmatvec=.true. )
!  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
!    constraint1=1, elemsub=stokes_constr_node_conn, addmatvec=.true. )

! periodical condition on velocities y-y

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=elementc, addmatvec=.true. )

  call check ( sysmatrix )


! solve initial velocity/pressure consistent with initial conditions
! and keep decomposition of the matrix
! This generates a linear profile.

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options_up%real_storage=rs_up
  solver_options_up%integer_storage=is_up

  call solve_system_ma57 ( sysmatrix, rhsd, sol, &
    solver_options=solver_options_up )

  call delete ( sysmatrix )

! create system vectors (solution and right-hand side) for diffuse interface and
! initialize vectors plus a small random perturbation.

  call create_sysvector ( problemdi, soldi, rhsdi, soldin, soldi_prev )

  call random_number (soldin%u)

  soldin%u(cval%s) = 1e-3*soldin%u(cval%s) ! initial c
  soldin%u(muval%s) = 0

  soldi%u = soldin%u


! create the structure oldvectors_di

  call create_oldvectors ( oldvectors_di, nsysvec=3, nprob=2 )

! store solution vectors and problem structures

  oldvectors_di%s(1)%p => sol
  oldvectors_di%s(2)%p => soldin
  oldvectors_di%s(3)%p => soldi
  oldvectors_di%p(1)%p => problem
  oldvectors_di%p(2)%p => problemdi


! Lees-Edwards conditions:

  H = 1; L = 1  ! height and length of the domain
  gammadot = U / H  ! shear rate
  time = 0 ! advancing time for deformation
  Udiff = U  ! velocity difference between upper and lower "wall"


! time stepping

  do step = 1, numtimesteps

    print *, 'step = ', step

!   shift mapped coordinates forward in time for Lees-Edwards conditions

    time = time + deltat

    call find_refcoor_objects ( mesh, object1=4, mapcoor=mapcoor )

!   create system matrix for diffuse-interface problem

    call create_sysmatrix_structure_base ( sysmatrixdi, mesh, problemdi )
    call create_sysmatrix_structure_constraint ( sysmatrixdi, mesh, problemdi )
    call finalize_sysmatrix_structure ( sysmatrixdi )

    call create_sysmatrix_data ( sysmatrixdi )

!   iteration scheme for (c,mu)

    iter = 0

    picard: do

      iter = iter + 1

!     build (assemble) matrix and vector for diffuse-interface problem

      call build_system ( mesh, problemdi, sysmatrixdi, sysvector=rhsdi, &
        elemsub=diffuse_interface_elem, oldvectors=oldvectors_di, &
        coefficients=coefficients )

!     periodical condition on (c,mu) x-x
      call build_system_constraint ( mesh, problemdi, sysmatrixdi, &
        constraint1=1, sysvector=rhsdi, elemsub=elementc, &
        addmatvec=.true. )
!      call build_system_constraint ( mesh, problemdi, sysmatrixdi, &
!        constraint1=1, sysvector=rhsdi, elemsub=stokes_constr_node_conn, &
!        addmatvec=.true. )

!     periodical condition on (c,mu) y-y
      call build_system_constraint ( mesh, problemdi, sysmatrixdi, &
        constraint1=2, sysvector=rhsdi, elemsub=elementc, &
        addmatvec=.true. )

      call check ( sysmatrixdi )

      soldi_prev%u = soldi%u

!     solve (c,mu)

      solver_options_di%real_storage=rs_di
      solver_options_di%integer_storage=is_di

      call solve_system_ma41 ( sysmatrixdi, rhsdi, soldi, &
        solver_options=solver_options_di )


      cmax  = maxval(abs(soldi%u(cval%s))) ! c maximum
      cdiff = maxval(abs(soldi%u(cval%s)-soldi_prev%u(cval%s))) ! c difference
      mumax = maxval(abs(soldi%u(muval%s))) ! mu maximum
      mudiff = maxval(abs(soldi%u(muval%s)-soldi_prev%u(muval%s))) ! mu diff

      print *, 'iter = ', iter
      print *, 'cmax = ', cmax, 'mumax = ', mumax
      print *, 'cdiff = ', cdiff, 'mudiff = ', mudiff

  !    if ( cdiff < epsdi * cmax .and. mudiff < epsdi * mumax ) exit picard

      if ( cdiff < epsdi * cmax ) exit picard

      if ( iter >= itermax ) then
        write(*,'(a,i0)') ' too many iterations: ', itermax
        stop
      end if

    end do picard

    call delete ( sysmatrixdi )

    soldin%u = soldi%u

!   create system matrix of velocity/pressure problem

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric=.true. )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

    call create_sysmatrix_data ( sysmatrix )

!   stokes velocity/pressure

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
      coefficients=coefficients )

!   build (assemble) vector for velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_mugradc, &
      oldvectors=oldvectors_di, physqrow=[physqvel], physqcol=[physqvel], &
      buildmatrix=.false., coefficients=coefficients, addmatvec=.true. )

!   periodical condition on velocities x-x

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=elementc, addmatvec=.true. )
!    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
!      constraint1=1, elemsub=stokes_constr_node_conn, addmatvec=.true. )

!   periodical condition on velocities y-y

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=elementc2, addmatvec=.true. )

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve velocity/pressure problem

    solver_options_up%real_storage=rs_up
    solver_options_up%integer_storage=is_up

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_up )

    call delete ( sysmatrix )

  end do


! write data for post-processing

  call write_mesh ( mesh, filename='mesh.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh, input_probdefdi, filename='probdefdi.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol%u
  write(10) soldi%u

  close(unit=10)


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, rhsd )
  call delete ( oldvectors_di )
  call delete ( problemdi )
  call delete ( input_probdefdi )
  call delete ( soldi, rhsdi, soldin, soldi_prev )
  call delete ( cval, muval )

  call delete ( coefficients )

end program diffuse_interface5
