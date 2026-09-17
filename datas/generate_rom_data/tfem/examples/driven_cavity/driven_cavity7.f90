! Viscoelastic problem on a unit square with Dirichlet boundary conditions.
! DEVSS-G/SUPG (explicit stress formulation)
! Giesekus model with two viscoelastic modes
! Second-order semi-implicit Gear time integration. First step: Euler implicit.
! One freely moving and rotating object.
! Weak constraint on the particle domain.
! Neutrally buoyant particle ( density of particle = density of fluid )

module subs_m

  use kind_defs_m

  implicit none

  real(dp) :: rp = 0, xp(2) = 0   ! radius and position of the particle

contains

! element subroutine for the weak constraints on the objects

  subroutine elementc ( mesh, problem, constr, elem, node, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemmat2, elemmatadd, &
    elemvec, elemvecadd )

    use tfem_elem_m

    implicit none

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer, parameter :: ndf = 9
    integer, parameter :: nodalpb = 9
    integer, parameter :: ndflb = 4  ! number of degrees of freedom of the Lag. Mult.

    integer :: object
    integer :: i, j
    real(dp) :: phi(1,ndf), xr(1,2), r(2), psi(1,ndflb), x(nodalpb,2)
    real(dp) :: theta(1,nodalpb), dtheta(1,nodalpb,2)
    real(dp) :: F(1,2,2), Finv(1,2,2), detF(1)


    object = problem%constraints(constr)%object

!   reference coordinates of the integration point (node) of the element (elem)

    xr(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)
    r = mesh%objects(object)%coor_int(node,:,elem) - xp

!   shape function of the velocity in the integration point

    call shape_quad_Q2 ( xr, phi )

!   shape function of the Langrangian multiplier in the integration point

    call shape_quad_Q1 ( mesh%objects(object)%xig(node:node,1:2), psi )

!   shape function of the quadratic particle element in the integration point

    call shape_quad_Q2 ( mesh%objects(object)%xig(node:node,1:2), theta, dtheta )

!   compute deformed element

    call get_coordinates_object ( mesh, elem, x, object )

    call isoparametric_deformation ( x, dtheta, F, Finv, detF )

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

!   two constraints (vectorial)
!
!      u - up - omega x r = 0
!      -   -      -     -   -
!
!   or in components
!
!      u - up + omega x ry = 0
!      v - vp - omega x rx = 0
!
!   the matrix A therefore becomes
!
!     A =  [ phi    0 ]
!          [   0  phi ]
!
!   and the matrix A_add of the additional unknowns (up,vp,omega)
!
!     A_add =  [ -1  0  ry ]
!              [  0 -1 -rx ]
!
!   In the weak version the equations are multiplied by the test function psi
!   and integrated over the element.
!   NOTE: this is the contribution of a single integration point (node) in
!   a single element (elem).

    if ( matrix ) then

      elemmat = 0

      do i = 1, ndflb
        do j = 1, ndf
          elemmat(i,j) = psi(1,i) * phi(1,j) * &
                                   detF(1) * mesh%objects(object)%wg(node)
        end do
      end do
      elemmat(ndflb+1:2*ndflb,ndf+1:2*ndf) = elemmat(1:ndflb,1:ndf)

      do i = 1, ndflb
        elemmatadd(i,:) = psi(1,i) * [ -1._dp,  0._dp,  r(2) ] * &
                                   detF(1) * mesh%objects(object)%wg(node)
        elemmatadd(i+ndflb,:) = psi(1,i) * [  0._dp, -1._dp, -r(1) ] * &
                                   detF(1) * mesh%objects(object)%wg(node)
      end do

    end if

    elemmat2 = 0._dp


  end subroutine elementc

end module subs_m



program driven_cavity7

  use tfem_m
  use hsl_ma41_m
  use stokes_elements_m
  use inertia_elements_m
  use viscoelastic_elements_m
  use io_utils_m
  use figplot_m
  use subs_m
  use mesh_circular_particle_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    gintpl = 4,         &  ! Q1 gradients
    cintpl = 4,         &  ! Q1 conformation
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    physqgrad = 3,      & ! physical quantity nr of the gradients
    gauss = 3,          & ! 3x3 integration of quads
    ncompc = 3,         & ! number of conformation tensor components
    nmodes = 2,         & ! number of modes
    startm = 501,       & ! start of material model data
    model = 3,          & ! Giesekus model
    nx=50,              & ! number of elements in x
    ny=50,              & ! number of elements in y
    nxp=3,              & ! number of elements in x (particle subdomain)
    nyp=5,              & ! number of elements in y (particle subdomain)
    timeint1 = 1,       & ! first-order implicit Euler time integration, first time step
    timeint2 = 5,       & ! second-order implicit Gear time integration after first time step
    logc = 1,           & ! standard scheme or log transformation
    num_picard = 1,     & ! number of Picard iterations, followed by Newton iterations
    numeulertimesteps = 1 ! number Euler time steps (at least one!)


  real(dp), parameter :: &
    rho = 1._dp,         & ! density of fluid (=density of particle)
    eta_s = 0.1_dp,      & ! solvent viscosity
    eta_p(nmodes) = [ 1.0_dp, 0.1_dp ],  & ! polymer viscosity
    lambda(nmodes) = [ 1.0_dp, 0.1_dp ], & ! relaxation time
    mobility(nmodes) = [ 0.01_dp, 0.01_dp ], & ! mobility parameter
    beta = 1.0_dp          ! upwinding parameter in the SUPG method


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh_particle
  type(input_probdef_t) :: input_probdef, input_probdefc
  type(problem_t), target :: problem, problemc
  type(sysmatrix_t) :: sysmatrix, sysmatrix_stored, sysmatrixc
  type(sysvector_t), target :: sol_iter, sol_hat, soln, solnm1
  type(sysvector_t) :: sol, rhsd, rhsd_stored
  type(sysvector_t), dimension(ncompc,nmodes), target :: solcn, solcnm1
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc
  type(vector_t) :: velocity, pressure, vorticity
  type(vector_t) :: Cxx, Cxy, Cyy
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors, oldvectors_hat, oldvectors_iter
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vel, pres
  type(lu_ma41_t) :: luc
  type(solver_options_ma41_t) :: solver_options

  integer :: step, iterstep, maxnumiterations
  integer :: numtimesteps
  real(dp) :: H, U
  real(dp) :: diffv, thresh, diffp, deltat, mfac
  real(dp) :: tn, tnp1
  real(dp) :: gamma0=1.5_dp, alpha0=2._dp, alpha1=-0.5_dp
  real(dp) :: alpha, G(nmodes)

  integer :: i, m
  real(dp) :: up(3), unm1(3), un(3), pp(3), xp_n(2)


  solver_options%real_storage=1.5
  solver_options%integer_storage=2.0

! number time steps:
  numtimesteps = 10
! time step
  deltat = 5.e-3_dp

! threshold for the iteration process
  thresh = 1.e-8_dp
! maximum number of iterations allowed to obtain convergence
  maxnumiterations = 100

! mesh width and heigth:
  H = 1._dp
! velocity of the lid
  U = 1._dp

! set some parameters

  alpha = sum(eta_p)  ! DEVSS parameter
  G = eta_p / lambda  ! modulus


! fill coefficients
  call create_coefficients ( coefficients, ncoefi=250, ncoefr=500+3*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gauss,    cintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      logc,     timeint1,    ( 0, i = 23, 250 )  &
    ]

  coefficients%r(1:500) = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
      0._dp, 0._dp,  deltat,   beta,  0._dp, &
      ( 0._dp, i = 11, 500 ) &
    ]
  coefficients%r(501:500+3*nmodes) = &
       [  ( G(i), lambda(i), mobility(i), i=1,nmodes )  ]
  coefficients%r(151) = rho

  call write_coefficients ( coefficients, filename='coefficients.out' )


! create mesh for fluid domain
  meshgen_options%elshape = 6
  meshgen_options%ox = 0._dp
  meshgen_options%oy = 0._dp
  meshgen_options%lx = H
  meshgen_options%ly = H
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )


! create mesh for particle domain

  rp = 0.075   ! radius of particle
  xp = [0.5_dp,0.5_dp] ! initial position of the center of the particle
  pp = [ (xp(i),i=1,2), 0._dp ] ! initial position and rotation

  call create_mesh_circular_particle ( mesh_particle, elshape=6, &
    nx=nxp, ny=nyp, rp=rp, xp=xp )

! add particle objects to the mesh

  call add_to_mesh ( mesh, object='mesh', objectmesh=mesh_particle, &
                    topology=.true., intrule=1, nsubint=3 )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )
  plot_options%objectpointcolor=4
  plot_options%objectpointsize=0.1
  call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true. )

  call fill_mesh_parts ( mesh_particle )
  call plot_mesh ( plot_options, mesh_particle, 'mesh_particle.fig' )
  call delete ( mesh_particle )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=3 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 4,0,4,0,4,0,4,0,0,    &  ! G
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,4] )

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=physqvel )
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! define constraints on the object

  call define_constraint ( mesh, input_probdef, object=1, physq=physqvel, &
    discretization='weak', elementdof=[2,0,2,0,2,0,2,0,0], naddunknowns=3 )

  call problem_definition ( input_probdef, mesh, problem )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=2, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a =   &
      reshape ( [ 1,0,1,0,1,0,1,0,0,    &  ! c
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar for plotting
                  [9,2] )

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

  call problem_definition ( input_probdefc, mesh, problemc )

! create system matrix for conformation problem

  call create_sysmatrix_structure ( sysmatrixc, mesh, problemc )

  call create_sysmatrix_data ( sysmatrixc )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, sol_iter, soln )
  call create_sysvector ( problem, solnm1, sol_hat )
  call create_sysvector ( problem, rhsd, rhsd_stored )

  call create ( problemc, solcn, solcnm1, rhsc )

! set initial value of solution to zero

  soln%u = 1.e-15
  sol_iter%u = soln%u
  un = 0._dp

  do m = 1, nmodes
    if ( logc == 0 ) then ! standard
      solcn(1,m)%u = 1   ! initial cxx
      solcn(2,m)%u = 0   ! initial cxy
      solcn(3,m)%u = 1   ! initial cyy
    else if ( logc == 1 ) then ! log scheme
      solcn(1,m)%u = 0   ! initial sxx
      solcn(2,m)%u = 0   ! initial sxy
      solcn(3,m)%u = 0   ! initial syy
    end if
  end do


! define vector subscripts for direct manipulation of sysvector data

  ! velocities
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )
  ! pressures
  call create_subscript ( mesh, problem, pres, physqarr=[physqpress] )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=2, nsysvec2=2, nprob=2 )

  call create_oldvectors ( oldvectors_hat, nsysvec=1 )
  call create_oldvectors ( oldvectors_iter, nsysvec=1 )

! fill oldvectors

  oldvectors%s(1)%p => soln
  oldvectors%s(2)%p => solnm1
  oldvectors%s2(1)%p => solcn
  oldvectors%s2(2)%p => solcnm1
  oldvectors%p(1)%p => problem
  oldvectors%p(2)%p => problemc

  oldvectors_hat%s(1)%p => sol_hat  ! alpha0*un+alpha1*un-1
  oldvectors_iter%s(1)%p => sol_iter ! solution at end of previous iteration

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, physq=physqvel, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=physqvel, degfd=1, exclude=3, value=U )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqpress, value=0._dp )



! start time integration

  call tic

  tn = 0

  integration: do step = 1, numtimesteps

    tnp1 = tn + deltat

    write(*,'(/a,i0,a,es12.4,a/)') ' ** time step = ', step, &
                                   ' time = ', tnp1, ' ** '


    ! advance particle position with 2nd order Adams-Bashforth

    if ( step >= 2 ) then

      unm1 = un
      un = up
      pp = pp + deltat*(3*un/2 - unm1/2)

      xp_n = xp
      xp = pp(1:2)

      ! update coordinates of object
      mesh%objects(1)%coor(:,1) = mesh%objects(1)%coor(:,1) + ( xp(1) - xp_n(1) )
      mesh%objects(1)%coor(:,2) = mesh%objects(1)%coor(:,2) + ( xp(2) - xp_n(2) )

      ! update reference coordinates of object
      call find_refcoor_objects ( mesh )

    end if


!   change time integration scheme at the second time step
    if ( step == 2 ) then
      coefficients%i(22) = timeint2
    end if


!   build (assemble) matrix and vector for conformation problem

    call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
      elemsub=ce_supg_elem, oldvectors=oldvectors, &
      coefficients=coefficients )

    call check ( sysmatrixc )


!   copy previous conformation solution to older time step

    call copy ( solcn, solcnm1 )


!   solve conformation and keep LU decomposition in the loop over components

    do m = 1, nmodes
      do i = 1, ncompc
        call solve_system_ma41 ( sysmatrixc, rhsc(i,m), solcn(i,m), luc, &
          solver_options=solver_options )
      end do
    end do

    call delete ( luc )  ! remove LU decomposition and rebuild next time step


!   create system matrix

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, symmetric=.false. )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )
    call create_sysmatrix_data ( sysmatrix )


!   build (assemble) matrix and vector from elements

    ! Stokes flow : velocity / pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
      physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress], &
      coefficients=coefficients )

    ! DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=devssg_elem, &
      physqrow=[physqgrad,physqvel], physqcol=[physqgrad,physqvel], &
      addmatvec=.true., coefficients=coefficients )

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_divtau, &
      oldvectors=oldvectors, physqrow=[physqvel], physqcol=[physqvel], &
      buildmatrix = .false., addmatvec=.true., coefficients=coefficients )

    ! set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqgrad], physqcol=[physqpress], &
      zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqpress], physqcol=[physqgrad], &
      zeromatvec=.true. )

    ! weak constraints on particle domain
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, elemsub=elementc, &
      addmatvec=.true., coefficients=coefficients )

    ! build instability term ( rho du/dt )
    if ( step <= numeulertimesteps ) then
      ! implicit Euler
      sol_hat%u(vel%s) = soln%u(vel%s)
      mfac = 1._dp
    else
      ! second-order implicit Gear after Euler time steps
      sol_hat%u(vel%s) = alpha0*soln%u(vel%s) + alpha1*solnm1%u(vel%s)
      mfac = gamma0
    end if

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=inertia_elem_dudt, &
      coefficients=coefficients, oldvectors=oldvectors_hat, &
      physqcol=[physqvel], physqrow=[physqvel], &
      factormat=mfac, addmatvec=.true. )


!   store sysmatrix & rhsd for advection iteration

    call create_sysmatrix ( sysmatrix, sysmatrix_stored )
    call copy ( sysmatrix, sysmatrix_stored )
    call copy ( rhsd, rhsd_stored )

    open ( unit=10, file='particle_motion.out' )

!   interation procedure for advection term

    iterstep = 0

    iterate: do

      iterstep = iterstep + 1

      call copy ( sysmatrix_stored, sysmatrix )
      call copy ( rhsd_stored, rhsd )


!     first num_picard steps Picard iteration, then Newton iteration
!     the term un+1.grad un+1
      if ( iterstep <= num_picard ) then
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_picard, coefficients=coefficients, &
          oldvectors=oldvectors_iter, physqcol=[physqvel], &
          physqrow=[physqvel], addmatvec=.true. )
      else
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_newton, coefficients=coefficients, &
          oldvectors=oldvectors_iter, physqcol=[physqvel], &
          physqrow=[physqvel], addmatvec=.true. )
      end if


      call check ( sysmatrix )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

      call solve_system_ma41 ( sysmatrix, rhsd, sol )

!     compute the difference between iteration steps
      diffv = maxval ( sqrt ( ( sol%u(vel%s) - sol_iter%u(vel%s) )** 2 ) )
      diffp = maxval ( sqrt ( ( sol%u(pres%s) - sol_iter%u(pres%s) )** 2 ) )

      write(*,'(a,i0,2(a,es12.4))') &
        '    iterstep = ', iterstep, ' diffv = ', diffv, ' diffp = ', diffp

      call copy ( sol, sol_iter )

      if ( diffv < thresh ) exit iterate

      if ( iterstep == maxnumiterations ) then
        write(*,'(a,i0)') &
          'Maximum number of iterations reached = ', maxnumiterations
        stop
      end if

    end do iterate

    call delete_sysmatrix ( sysmatrix, sysmatrix_stored )

!   copy the converged solution to the old step

    call copy ( soln, solnm1 )
    call copy ( sol, soln )

!   extract info on particle

    call get_sysvector_constraint ( mesh, problem, soln, constraint=1, &
      addunknowns=.true., u=up )

    print *, 'pp = ', pp
    print *, 'up = ', up

    write(10, '(i5,7es16.8)') step, tnp1, pp, up

    tn = tnp1

  end do integration

  call toc

  close(unit=10)

! post-processing

  call create_vector ( problem, velocity, physq=physqvel )
  call create_vector ( problem, pressure, vec=input_probdef%nvec )
  call create_vector ( problem, vorticity, vec=input_probdef%nvec )

  call create_vector ( problemc, Cxx, vec=input_probdefc%nvec )
  call create_vector ( problemc, Cxy, vec=input_probdefc%nvec )
  call create_vector ( problemc, Cyy, vec=input_probdefc%nvec )


  call extract_physvector ( mesh, problem, sol, velocity )


  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=1
  call derive_vector ( mesh, problemc, Cxx, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=2
  call derive_vector ( mesh, problemc, Cxy, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=3
  call derive_vector ( mesh, problemc, Cyy, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a fig file for plotting

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  call plot_color_fill ( plot_options, mesh, problem, 'pressure_color.fig', &
    vector=pressure )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )

  call plot_color_fill ( plot_options, mesh, problem, 'vorticity_color.fig', &
    vector=vorticity )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

  call plot_color_contour ( plot_options, mesh, problemc, 'Cxx_contour.fig', vector=Cxx )

  call plot_color_contour ( plot_options, mesh, problemc, 'Cxy_contour.fig', vector=Cxy )

  call plot_color_contour ( plot_options, mesh, problemc, 'Cyy_contour.fig', vector=Cyy )


! write binary file for reading by streamfunction

  open(unit=10,form='unformatted',file='velocity_bin.out')

  write(10) velocity%u

  close(unit=10)

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

! write vtk data
  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='cavity_inertia.vtk' )

  call write_vector_vtk ( mesh, problem, filename='cavity_inertia.vtk', &
    dataname='velocity_vector', sysvector=sol, append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, sol_iter, soln )
  call delete ( solnm1, sol_hat)
  call delete ( rhsd, rhsd_stored )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( solcn, solcnm1, rhsc )
  call delete ( velocity, pressure, vorticity )
  call delete ( Cxx, Cxy, Cyy )
  call delete ( coefficients )
  call delete ( oldvectors, oldvectors_iter, oldvectors_hat )
  call delete ( vel, pres )

end program driven_cavity7
