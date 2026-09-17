! Navier-Stokes problem on gravitational sedimentation of a particle in a
! rectangular box.
! Second-order implicit Gear discretization. First step: Euler implicit.
! One freely moving and rotating object. Time-stepping.
! Weak constraint on the particle domain.

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



program navier_stokes9

  use tfem_m
  use math_defs_m
  use hsl_ma41_m
  use stokes_elements_m
  use inertia_elements_m
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
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=50,              & ! number of elements in x
    ny=75,              & ! number of elements in y
    nxp=2,              & ! number of elements in x (particle subdomain)
    nyp=3                 ! number of elements in y (particle subdomain)

  real(dp), parameter :: &
    eta = 0.1_dp           ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh_particle
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix, sysmatrix_stored
  type(sysvector_t), target :: sol_iter, soln, sol_hat
  type(sysvector_t) :: sol, solnm1, rhsd, rhsd_stored
  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors, oldvectors_hat, oldvectors_iter
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vel, pres

  integer :: step, iterstep, num_picard, maxnumiterations
  integer :: numtimesteps, numeulertimesteps
  real(dp) :: lx, ly, rho
  real(dp) :: diffv, thresh, diffp, deltat, mfac
  real(dp) :: tn, tnp1
  real(dp) :: gamma0=1.5_dp, alpha0=2._dp, alpha1=-0.5_dp

  integer :: i
  real(dp) :: rho_p, mass, moment, rhofactor
  real(dp) :: up(3), unm1(3), un(3), pp(3), xp_n(2)
  real(dp) :: up_hat(3), matdiag(3), vec(3), gr

! gravitational acceleration (in CGS units)
  gr = 9.81e+2
! density
  rho = 1._dp
  rho_p = 1.25_dp
! number time steps:
  numtimesteps = 170
! number Euler time steps (at least one!):
  numeulertimesteps = 1
! time step
  deltat = 5.e-3_dp

! threshold for the iteration process
  thresh = 1.e-8_dp
! maximum number of iterations allowed to obtain convergence
  maxnumiterations = 10
! number of Picard iterations, followed by Newton iterations
  num_picard = 1

! mesh width and height:
  lx = 2._dp
  ly = 3._dp

! fill coefficients
  call create_coefficients ( coefficients, ncoefi=250, ncoefr=200 )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(8) = deltat
  coefficients%r(151) = rho

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh for fluid domain
  meshgen_options%elshape = 6
  meshgen_options%ox = 0._dp
  meshgen_options%oy = 0._dp
  meshgen_options%lx = lx
  meshgen_options%ly = ly
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )


! create mesh for particle domain

  rp = 0.1   ! radius of particle
  xp = [1.5_dp,2.5_dp] ! initial position of the center of the particle
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

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

! define constraints on the object

  call define_constraint ( mesh, input_probdef, object=1, physq=1, &
    discretization='weak', elementdof=[2,0,2,0,2,0,2,0,0], naddunknowns=3 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, sol_iter, soln )
  call create_sysvector ( problem, solnm1, sol_hat )
  call create_sysvector ( problem, rhsd, rhsd_stored )

! define vector subscripts for direct manipulation of sysvector data

  ! velocities
  call create_subscript ( mesh, problem, vel, physqarr=[1] )
  ! pressures
  call create_subscript ( mesh, problem, pres, physqarr=[2] )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )
  call create_oldvectors ( oldvectors_hat, nsysvec=1 )
  call create_oldvectors ( oldvectors_iter, nsysvec=1 )


! set initial value of solution to zero

  soln%u = 0._dp
  sol_iter%u = soln%u
  un = 0._dp

! fill oldvectors

  oldvectors%s(1)%p => soln
  oldvectors_hat%s(1)%p => sol_hat  ! alpha0*un+alpha1*un-1
  oldvectors_iter%s(1)%p => sol_iter ! solution at end of previous iteration

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )


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


!   create system matrix

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem, &
      diagonal_addunknowns=.true. )
    call finalize_sysmatrix_structure ( sysmatrix )
    call create_sysmatrix_data ( sysmatrix )


!   build (assemble) matrix and vector from elements

    ! Stokes flow : velocity / pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
      coefficients=coefficients )

    ! weak constraints on particle domain
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      elemsub=elementc, addmatvec=.true., coefficients=coefficients )

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

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=inertia_elem_dudt, &
      coefficients=coefficients, oldvectors=oldvectors_hat, &
      physqcol=[physqvel], physqrow=[physqvel], &
      factormat=mfac, addmatvec=.true. )

!   build particle buoyancy due to the density difference

    mass = rho_p * pi * rp**2
    moment = 0.5 * mass * rp**2
    rhofactor = 1.0 - rho/rho_p

    if ( step <= numeulertimesteps ) then
      ! first order Euler
      up_hat = un
      mfac = 1._dp
    else
      ! second-order BDF2 after Euler time steps
      up_hat = alpha0*un + alpha1*unm1
      mfac = gamma0
    end if

    matdiag(1) = rhofactor * mass * mfac / deltat
    matdiag(2) = rhofactor * mass * mfac / deltat
    matdiag(3) = rhofactor * moment * mfac / deltat

    vec(1) = rhofactor * mass * up_hat(1) / deltat
    vec(2) = rhofactor * mass * up_hat(2) / deltat
    vec(3) = rhofactor * moment * up_hat(3) / deltat

    ! gravitational force on particle
    vec(2) = vec(2) - rhofactor * mass * gr

    call build_system_constraint_addunknowns ( mesh, problem, &
      sysmatrix=sysmatrix, sysvector=rhsd, constraint=1, &
      addmatvec=.true., matdiag=matdiag, rhs=vec )

!   store sysmatrix & rhsd for advection iteration

    call create_sysmatrix ( sysmatrix, sysmatrix_stored )
    call copy ( sysmatrix, sysmatrix_stored )
    call copy ( rhsd, rhsd_stored )

    open ( unit=10, file='out9' )

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

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )


  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
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
  call delete ( velocity, pressure, vorticity )
  call delete ( coefficients )
  call delete ( oldvectors, oldvectors_iter, oldvectors_hat )
  call delete ( vel, pres )

end program navier_stokes9


