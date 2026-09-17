#if SKIT2 && METIS5

! Drop consisting of a Newtonian fluid suspended in a Newtonian fluid
! with complex interfacial rheology as a boundary condition.
! Lagrangian description of the interface.
! Stokes flow (no inertia).
! Flow imposed either:
!   1) on the whole boundary using a function, or
!   2) on the upper and lower boundary using a function and
!      periodic boundary conditions in x- and z-direction (periodic=.true.)
!   3) periodic boundary conditions in x-direction and on the other boundaries
!      imposed by a function (periodic_channel=.true.)
!   4) periodic boundary conditions in xyz-directions (triperiodic=.true.)
!      and imposing an applied macroscopic shear rate.
! Optionally apply the method of manufacturing solutions to verify
! correctness of implementation and validate the order of convergence.
! Compute and plot the interfacial stress.
! ALE with interface tracking and optional displacement of the mesh in
! x-direction to keep the center of volume at the same position in the mesh.
! Optionally remeshing is performed when the maximum normalized volume or
! aspect ratio of the mesh elements is higher than a threshold.
! Remeshing is based on position of nodes of the interface mesh, but the
! interface mesh itself is not modified.
! Optionally move the nodes of the interface mesh to highly curved regions.
! (for a Lagrangian description of the interface, this does not make sense!).
! Optionally impose flowrate if periodic=.true.
! Optionally compute the error if the exact solution is known.
! Optionally compute the volume and center of volume.
! Three-dimensional model.
! Optionally Sparskit iterative solvers.
! Metis 5 renumbering.

program drop17

  use tfem_m
  use hsl_ma57_m
  use hsl_ma41_m
  use sk_solve_m
  use stokes_elements_m
  use interfacial_rheology_elements_3D_surface_m
  use interface_tracking_elements_m
  use poisson_gd_ma57_m
  use curvature_gd_ma57_m
  use update_mesh_nodes2_m
  use functions_m
  use io_utils_m
  use figplot_m
  use subs_m
  use metis5_m
  use timer_m
  use eig2D3D_m, only: eig3x3

  implicit none


! constants

  integer, parameter ::  &
    uintpl = 6,          & ! P2 velocity and interface position interpolation
    pintpl = 2,          & ! P1 pressure interpolation
    physqvel = 1,        & ! physical quantity nr of velocity
    physqpress = 2,      & ! physical quantity nr of pressure
    ndim = 3               ! dimension of space


! variables

  logical ::     &
    periodic = .false.,     & ! Use periodical boundaries in x and z-direction
                              ! Only makes sense for vfuncnr=0 (zero velocity)
                              ! or vfuncnr=2 (shear flow).
    triperiodic = .false.,  & ! Use periodical boundaries in xyz directions
    periodic_channel = .false., & ! Use periodical boundaries in x-direction.
                              ! Only makes sense for vfuncnr=0 (zero velocity)
                              ! or vfuncnr=2 (shear flow).
    flowrate = .false.,     & ! Impose flowrate in x-direction (for
                              ! periodic=.true. or periodic_channel=.true. only)
    move_mesh_x = .true.,   & ! Move the mesh in x-direction to keep the center
                              ! of volume of the drop at the same position wrt
                              ! to the mesh
    monitor_func = .false., & ! Use a monitor function to move the surface
                              ! elements tangential to high curvature regions
    ctime = .false.,        & ! continue time in output after restart
    mesh_plot = .false.,    & ! write mesh plot files/info (curves, mesh etc.)
    compute_error = .false.,& ! Compute error. Only makes sense for
                              ! manu_sol=.true. or if eta_e=eta_i
                              ! (no viscosity difference) and gammac=0 (no
                              ! surface tension) or if vfuncnr=0 (no flow)
    compute_volume = .false.,&! Compute volume and center of volume.
    remesh = .true., &        ! perform remeshing
    manu_sol = .true., &      ! Use method of manufactured solutions
    iter_solv = .true., &     ! Use iterative solver
    lagrange = .true.         ! Use full Lagrangian interface tracking

  integer ::  &
    restart=0,           & ! restart=0: no restart
                           ! restart=1: normal restart
                           ! restart=2: restart, but start with Euler step
    step0=0,             & ! initial step number
    post0=0,             & ! initial post number
    n=1,                 & ! number of elements in x, y and z direction of box
    nc=5,                & ! number of elements on the equator of the sphere
    gauss  = 4,          & ! order of integration for tetrahedral elements
    gaussb = 4,          & ! order of integration for triangular elements
    vfuncnr = 2,         & ! function number for the velocity field u
    choice_um = 1,       & ! choice for um in the motion of the interface:
                           !  0: um = 0
                           !  1: um = average velocity of all interface nodes
                           !  2: um = velocity of the center of volume
    vtkevery = 5 ,       & ! vtk file every vtkevery steps. 0: means none
    restart_every = 100, & ! write restart file every restart_every steps.
                           ! 0: means none
    timeint=1,           & ! time integration scheme
    numtimesteps=400,    & ! number of time steps
    imodel = 1,          & ! interface model number:
                           ! 0: surface tension,
                           ! 1: Boussinesq-Scriven,
                           ! 2: Verwijlen,
                           ! 3: Kelvin-Voigt with Verwijlen,
                           ! 4: Huetter-Tervoort,
                           ! 5: Kelvin-Voigt with Huetter-Tervoort.
    surfimpl = 1           ! implicit(=1) or explicit(=0) surface tension

  real(dp) ::  &
    ox = -0.5_dp,         & ! x-coordinate left lower corner
    oy = -0.5_dp,         & ! y-coordinate left lower corner
    oz = -0.5_dp,         & ! z-coordinate left lower corner
    lx = 1.0_dp,          & ! length of the domain
    ly = 1.0_dp,          & ! width of the domain
    lz = 1.0_dp,          & ! height of the domain
    xc(ndim) = [0.0_dp,0.0_dp,0.0_dp],  & ! center of the initial sphere
    rc = 0.1_dp,          & ! radius of the initial sphere
    eta_e = 1.0_dp,       & ! external fluid viscosity
    eta_i = 1.0_dp,       & ! internal fluid viscosity
    muc = 1.0_dp,         & ! surface shear viscosity
    kappac = 1.0_dp,      & ! surface dilatational viscosity
    Gc = 1.0_dp,          & ! surface shear modulus
    Kc = 1.0_dp,          & ! surface dilatational modulus
    gammac = 0.5_dp,      & ! surface tension coefficient
    factor = 100._dp,     & ! factor for the interface grid deformation
    rs_up  = 1.4_dp,      & ! real_storage velocity-pressure LU HSL
    is_up  = 1.5_dp,      & ! integer_storage velocity-pressure LU HSL
    U_avg  = 1.0_dp,      & ! imposed average velocity for flowrate=.true.
    shear_rate  = 1.0_dp, & ! imposed shear rate for triperiodic=.true.
    time0  = 0.0_dp,      & ! initial time
    beta = 0.5_dp,        & ! factor beta in SUPG
    deltat = 7.e-3_dp,    & ! time step
    eps_sk = 1.e-6_dp,    & ! relative tolerance used in sparskit
    fac_min = 0.1_dp,     & ! factor setting dx_min for the monitor function
    fac_max = 10.0_dp,    & ! factor setting dx_max for the monitor function
    volume_threshold = 1.39_dp,      & ! remeshing volume threshold
    aspect_ratio_threshold = 1.39_dp   ! remeshing aspect ratio threshold

! flow problem

  type(mesh_t) :: mesh_ext, mesh_int, mesh, mesh_conform
  type(mesh_t), target :: mesh_0
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients(2)
  type(plot_options_t) :: plot_options
  type(solver_options_ma57_t) :: solver_options_ma57
  type(solver_options_sk_t) :: solver_options
  type(ilu_sk_t) :: ilu
  type(subscript_t) :: ssvel_b, press

! curvature problem

  type(problem_t) :: problem_curv
  type(sysmatrix_t) :: sysmatrix_curv

! interface tangential motion problem

  type(problem_t), target :: problem_pois
  type(sysmatrix_t) :: sysmatrix_pois
  type(sysvector_t), target :: sol_pois

! interface tracking problem

  type(mesh_t), target :: mesh_b, mesh_b0
  type(mesh_t) :: mesh_b_temp
  type(input_probdef_t) :: input_probdef_b
  type(problem_t), target :: problem_b
  type(sysmatrix_t) :: sysmatrix_b
  type(sysvector_t), target :: sol_b
  type(sysvector_t) :: rhsd_b(ndim)
  type(vector_t) :: xnp1, error, stress_scal, stress_ten, tause, tauerr
  type(vector_t) :: lambda(2), eigv(2)
! Interface mesh coordinates at previous times. NOTE: the coordinates are
! defined with respect to the initial coordinate system.
  type(vector_t), target :: xn, xnm1
  type(vector_t), target :: vel_b
  type(oldvectors_t) :: oldvectors_b
  type(coefficients_t) :: coefficients_b, coefficients_i
  type(lu_ma41_t) :: lu_b
  type(solver_options_ma41_t) :: so


! ALE mesh motion problem

  type(problem_t), target ::  problem_update_mesh
  real(dp), allocatable, dimension(:,:) :: disp

! some more misc variables

  logical :: w1, euler_step, euler_step2, maxmvmreached
  character(len=30) :: filename

  integer :: nnodes, nelem, elshape, i, c_cpl, c_perx, c_pery, c_perz, c_flr
  integer :: vertices(4)=[1,3,5,10], post, step, a(1)
  real(dp) :: xp(1,ndim), rp(1), dx_box, dx_part, resultsum(1+ndim)
  real(dp) :: dropvolume(1), L2norm(2), intpress(1), inertia_ten(1+2*ndim)
  real(dp) :: normf, eigval(3), eigvec(3,3), D1, D2, theta, l1, l2, l3

! deformation variables for remeshing
  real(dp), dimension(:), allocatable :: init_volume, volumev
  real(dp), dimension(:), allocatable :: init_aspect_ratio, aspect_ratio
  real(dp) :: norm_volume, norm_aspect_ratio
  real(dp), dimension(:), allocatable :: f_monitor

! Center of volume position of the drop at current and previous times.
! Used for the additional ALE movement of the whole mesh in x-direction to keep
! the x-coordinate of the center of volume of the droplet, with respect to the
! moving mesh, constant in time. Also used for um = velocity of the center of
! volume.
! NOTE: these are defined according to the position of the full mesh, which is
!       based on the PREDICTION of the position of the interface mesh.
! NOTE: the coordinates are defined with respect to the initial coordinate
!       system. Hence, the actual coordinates of the mesh are xm+cv-xc,
!       where xm are the coordinates in mesh%coor and xc are the initial
!       coordinates of the center of volume.
  real(dp), dimension(ndim) :: cv, cvn, cvnm1

! center of volume velocity
  real(dp), dimension(ndim) :: cvv


! namelist for input of variables; read from standard input

  namelist /comppar/ periodic, triperiodic, periodic_channel, flowrate, &
    move_mesh_x, remesh, manu_sol, lagrange, iter_solv, eps_sk, &
    monitor_func, mesh_plot, compute_error, compute_volume, n, nc, &
    gauss, gaussb, choice_um, vfuncnr, vtkevery, &
    ctime, restart, restart_every, timeint, numtimesteps, &
    ox, oy, oz, lx, ly, lz, xc, rc, fac_min, fac_max, eta_e, eta_i, &
    imodel, surfimpl, muc, kappac, Gc, Kc, gammac, factor, U_avg, shear_rate, &
    rs_up, is_up, deltat, beta, volume_threshold, aspect_ratio_threshold

  read ( unit=*, nml=comppar )


! Fill coefficients

! External fluid (element group 1)

  call create ( coefficients, ncoefi=600, ncoefr=600 )

  coefficients(1)%i = 0
  coefficients(1)%i(1:11) = &
  [ uintpl,   pintpl,     0,     0,         0,  &
    physqvel, physqpress, 0,     0,     gauss,  &
    gaussb ]
  coefficients(1)%i(40) = 3 ! standard Gauss-Legendre (numerical tables)

  coefficients(1)%r = 0
  coefficients(1)%r(1) = eta_e
  coefficients(1)%r(6) = - U_avg * ly * lz ! flow rate

! Internal fluid (element group 2)

  coefficients(2)%i = coefficients(1)%i
  coefficients(2)%i(63) = vfuncnr

  coefficients(2)%r = 0
  coefficients(2)%r(1) = eta_i

  coefficients(2)%vfunc => vfunc_3D

! Interface tracking

  call create_coefficients ( coefficients_b, ncoefi=100, ncoefr=50 )

  coefficients_b%i = 0
  coefficients_b%i(1) = uintpl
  coefficients_b%i(3) = gaussb
  coefficients_b%i(5) = 3  ! standard Gauss-Legendre (numerical tables)
  coefficients_b%i(6) = 2  ! velocity u on the interface given by a vector
  if ( .not. lagrange ) &
    coefficients_b%i(10) = 1  ! tangential grid velocity using Poisson problem
  coefficients_b%i(11) = timeint

  coefficients_b%r = 0
  coefficients_b%r(7) = factor
  coefficients_b%r(8) = beta
  coefficients_b%r(9) = deltat

! Interfacial rheology

  call create_coefficients ( coefficients_i, ncoefi=100, ncoefr=50 )

  coefficients_i%i = 0
  coefficients_i%i(1) = uintpl
  coefficients_i%i(3) = gaussb
  coefficients_i%i(5) = 3  ! standard Gauss-Legendre (numerical tables)
  coefficients_i%i(10) = vfuncnr
  coefficients_i%i(11) = imodel
  coefficients_i%i(12) = surfimpl

  coefficients_i%r = 0
  coefficients_i%r(1) = deltat
  coefficients_i%r(2) = gammac
  coefficients_i%r(3) = kappac
  coefficients_i%r(4) = muc
  coefficients_i%r(5) = Kc
  coefficients_i%r(6) = Gc


! create mesh

  dx_box = lx/n
  dx_part = 2*pi*rc/nc

  if ( restart == 0 ) then

!   generate fluid mesh using gmsh

    xp(1,:) = xc
    rp(1) = rc

    call generate_and_read_mesh

  else

!   restart

    call read_mesh ( mesh, 'mesh.out' )

  end if

  call fill_mesh_parts ( mesh )

  if ( mesh_plot ) then

    call plot_points_curves ( plot_options, mesh, filename='curves.fig' )
    call plot_mesh ( plot_options, mesh, filename='mesh.fig' )

    call write_mesh_gmsh ( mesh, filename='mesh.msh' )
    call write_mesh_vtk ( mesh, filename='mesh.vtk' )

    call printinfo ( mesh, printlevel=4 )

  end if


! generate interface mesh from the fluid mesh

  nnodes = mesh%surfaces(7)%nnodes
  nelem = mesh%surfaces(7)%nelem
  elshape = mesh%surfaces(7)%element%elshape

  call mesh_skeleton ( mesh_b, nnodes, nelem, elshape, ndim )

  mesh_b%topology(1)%a=mesh%surfaces(7)%topology(:,:,1)
  mesh_b%coor=mesh%coor(mesh%surfaces(7)%nodes,:)

  call fill_mesh_parts ( mesh_b )

  mesh_b0 = mesh_b

  if ( mesh_plot ) then

    call printinfo ( mesh_b, printlevel=4 )

    call write_mesh_vtk ( mesh_b, filename='mesh_b.vtk' )

  end if


! Velocity/pressure problem definition in external and internal fluid

  if ( iter_solv ) then
    call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2, &
      nphysqshifted=1 )
  else
    call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )
  end if

  input_probdef%vec_elementdof(1)%a(:,1) = ndim      ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0         ! pressure
  input_probdef%vec_elementdof(1)%a(vertices,2) = 1
  input_probdef%vec_elementdof(1)%a(:,3) = 1         ! scalar
  input_probdef%vec_elementdof(1)%a(:,4) = ndim*(ndim+1)/2  ! symmetric tensor

  input_probdef%vec_elementdof(2)%a = input_probdef%vec_elementdof(1)%a

  input_probdef%physq = [physqvel,physqpress]
  if ( iter_solv ) input_probdef%physqshifted = [2]
  input_probdef%probnr = 1

  if ( periodic ) then

!   dirichlet on lower and upper wall

    call define_essential ( mesh, input_probdef, surface1=1, physq=physqvel )
    call define_essential ( mesh, input_probdef, surface1=3, physq=physqvel )

!   velocity periodic in x-direction

    call define_constraint ( mesh, input_probdef, physq=physqvel, &
      surface1=2, surface2=4, discretization='collocation', &
      excludesurfaces=[1,3], num=c_perx )

!   velocity periodic in z-direction

    call define_constraint ( mesh, input_probdef, physq=physqvel, &
      surface1=5, surface2=6, discretization='collocation', &
      excludesurfaces=[1,2,3], num=c_perz )

  else if ( triperiodic  ) then

    call define_essential ( mesh, input_probdef, point=1, physq=physqvel )

!   velocity periodic in x-direction

    call define_constraint ( mesh, input_probdef, physq=physqvel, &
      surface1=2, surface2=4, discretization='collocation', num=c_perx )

!   velocity periodic in z-direction

    call define_constraint ( mesh, input_probdef, physq=physqvel, &
      surface1=5, surface2=6, discretization='collocation', &
      excludesurfaces=[2], num=c_perz )

!   velocity periodic in y-direction

    call define_constraint ( mesh, input_probdef, physq=physqvel, &
      surface1=1, surface2=3, discretization='collocation', &
      excludesurfaces=[2,5], num=c_pery )

  else if ( periodic_channel ) then

!   dirichlet on side walls of the channel

    call define_essential ( mesh, input_probdef, surface1=1, physq=physqvel )
    call define_essential ( mesh, input_probdef, surface1=3, physq=physqvel )
    call define_essential ( mesh, input_probdef, surface1=5, surface2=6, &
      physq=physqvel )

!   velocity periodic in x-direction

    call define_constraint ( mesh, input_probdef, physq=physqvel, &
      surface1=2, surface2=4, discretization='collocation', &
      excludesurfaces=[1,3,5,6], num=c_perx )

  else

!   dirichlet on all flowcell boundaries

    call define_essential ( mesh, input_probdef, surface1=1, surface2=6, &
      physq=physqvel )

  end if

  if ( flowrate .and. ( periodic .or. periodic_channel ) ) then

!   impose flowrate

    call define_constraint ( mesh, input_probdef, physq=physqvel, &
      surface1=2, nglobalc=1, num=c_flr )

  end if

  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! Internal fluid - external fluid coupling through Lagrangian multipliers

  call define_constraint ( mesh, input_probdef, surface1=7, surface2=8, &
    discretization='collocation', physq=physqvel, num=c_cpl )

! Problem definition and create sysmatrix, vectors

  call problem_definition_create_sysmatrix_vectors

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nmesh=2, nprob=1, nvec=1, nsysvec=1 )

! derive vectors

  mesh_0 = mesh
  oldvectors%s(1)%p => sol
  oldvectors%m(1)%p => mesh_b
  oldvectors%m(2)%p => mesh_0
  oldvectors%p(1)%p => problem_b
  oldvectors%v(1)%p => xn

! monitor function

  if ( monitor_func ) allocate ( f_monitor(mesh_b%nnodes) )


! problem definition of the interface

  call create_input_probdef ( mesh_b, input_probdef_b, nvec=3 )

  input_probdef_b%elementdof(1)%a = 1
  input_probdef_b%vec_elementdof(1)%a(:,1) = ndim  ! 2D vector
  input_probdef_b%vec_elementdof(1)%a(:,2) = 1  ! error
  input_probdef_b%vec_elementdof(1)%a(:,3) = ndim*(ndim+1)/2 ! symmetric tensor

  input_probdef_b%probnr = 2

  call problem_definition ( input_probdef_b, mesh_b, problem_b )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem_b, sol_b )
  call create_msysvector ( problem_b, rhsd_b )

! storage for position at tn-1, tn, tnp1 in vector

  call create ( problem_b, xn, vec=1 )
  call create ( problem_b, xnp1, vec=1 )

  if ( timeint == 2 ) then
    call create ( problem_b, xnm1, vec=1 )
  end if

! Initialize xnp1 (which will be the xn of the new time step)

  xnp1%u = reshape ( transpose(mesh_b%coor), [ndim*mesh_b%nnodes] )

! storage for velocity in vector

  call create ( problem_b, vel_b, vec=1 )

  if ( compute_error ) then

!   storage for error

    call create ( problem_b, error, vec=2 )
    call create ( problem_b, tause, vec=3 )
    call create ( problem_b, tauerr, vec=3 )

  end if

  call create ( problem_b, eigv, vec=1 )
  call create ( problem_b, lambda, vec=2 )
  call create ( problem_b, stress_scal, vec=2 )
  call create ( problem_b, stress_ten, vec=3 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix_b, mesh_b, problem_b )

  call create_sysmatrix_data ( sysmatrix_b )

! create oldvectors

  call create ( oldvectors_b, nsysvec=1, nvec=4, nprob=1, nmesh=1 )

  oldvectors_b%v(1)%p => vel_b  ! velocity u on the interface
  oldvectors_b%v(3)%p => xn     ! position at tn
  if ( timeint == 2 ) then
    oldvectors_b%v(4)%p => xnm1  ! position at tn-1
  end if
  oldvectors_b%m(1)%p => mesh_b0


! allocate arrays for the ALE displacement problem

  allocate ( disp(mesh_b%nnodes,ndim) )

! Create aspect ratio arrays

  call create_aspect_ratio

! inititalize cv, ... to xc because they are all written to restart

  cv = xc; cvn = xc; cvnm1 = xc


! open files and restart

  call open_files_and_restart

! imposed shear rate for stokes_constr_node_conn_shear
  lshear_rate = shear_rate
  lly = ly

! time stepping

  post = post0
  euler_step = .true.
  euler_step2 = .true.

  call tic

  do step = 1, numtimesteps

    write(*,'(a,i0)') 'step = ', step0+step

!   compute current aspect ratio
    call compute_element_aspect_ratio ( mesh, volumev, aspect_ratio )

!   compute maximum normalized aspect ratio and volume
!   (with respect to the initial values)
    norm_volume = maxval( abs(log(volumev/init_volume)) )
    norm_aspect_ratio = maxval( abs(log(aspect_ratio/init_aspect_ratio)) )

!   remeshing criterion

    if ( ( norm_volume >= volume_threshold .or. &
           norm_aspect_ratio >= aspect_ratio_threshold ) .and. &
           remesh .and. step /= numtimesteps ) then

      print *,'Doing remeshing...'

      call remeshing

    end if

!   move the interface position and solve flow

    call solve_flow_and_move_interface

!   error of interface position

    if ( compute_error ) then

!     fill error vector

      time = time0 + step * deltat

      call fill_vector ( mesh_b, problem_b, error, &
        node1=1, node2=mesh_b%nnodes, func=func_error_3D, funcnr=vfuncnr )

      write(unit=10,fmt=*) time, maxval(error%u), &
                           sqrt(sum(error%u**2)/mesh_b%nnodes)

    end if

!   volume and center of volume of the interface

    if ( compute_volume ) then

!     compute volume and center of volume

      call integrate ( mesh_b, problem_b, resultsum, elemsub=center_of_volume, &
         coefficients=coefficients_b )

      time = time0 + step * deltat

      write(unit=12,fmt=*) time, resultsum(1), resultsum(2:)/resultsum(1)

!     compute inertia tensor for Taylor deformation and inclination

      mesh_b_temp = mesh_b
      do i = 1, mesh_b%nnodes
        mesh_b_temp%coor(i,:) = mesh_b_temp%coor(i,:) - &
        resultsum(2:)/resultsum(1)
      end do

      call integrate ( mesh_b_temp, problem_b, inertia_ten, &
        elemsub=second_moment_of_area, coefficients=coefficients_b )

      call eig3x3 ( inertia_ten(2:), eigval, eigvec )

      a = maxloc(eigval)
      l1 = sqrt(eigval(a(1)))
      l2 = sqrt(maxval(eigval,[(i, i=1,size(eigval))]/=a(1)))
      l3 = sqrt(minval(eigval))
      D1 = (l1-l3)/(l1+l3)
      D2 = (l1-l2)/(l1+l2)
      theta = atan2(eigvec(2,a(1)),eigvec(1,a(1)))

      write(unit=15,fmt=*) time, D1, D2, theta
      write(unit=16,fmt=*) time, inertia_ten

    end if

!   output L2-norm of pressure and velocity in drop

    if ( manu_sol ) then
      call integrate ( mesh, problem, dropvolume, stokes_integrate_volume, &
        coefficients(2), oldvectors=oldvectors, elgroup1=2 )

      call integrate ( mesh, problem, intpress, stokes_integrate_pressure, &
        coefficients(2), oldvectors=oldvectors, elgroup1=2 )

      coefficients(2)%r(23) = intpress(1) / dropvolume(1)

      call integrate ( mesh, problem, L2norm, stokes_L2_norm, &
        coefficients(2), oldvectors=oldvectors, elgroup1=2 )

      write(unit=18,fmt=*) time0+step*deltat, sqrt(L2norm/dropvolume(1))

!     output minimum, maximum and average value of pressure

      write(unit=14,fmt=*) time0+step*deltat, minval(sol%u(press%s)), &
                                              maxval(sol%u(press%s)), &
                                              intpress/dropvolume
    end if

!   output to vtk

    if ( vtkevery > 0 ) then
      if ( mod(step0+step,vtkevery) == 0 ) then

        post = post + 1

        call output_data

      end if
    end if

!   write data for restart

    if ( restart_every > 0 ) then
      if ( mod(step0+step,restart_every) == 0 ) then

        call output_restart

      end if
    end if

    call toc('end of time step')

  end do

  call toc

  if ( compute_error ) close(unit=10)
  if ( compute_error ) close(unit=17)
  if ( compute_volume ) close(unit=12)
  if ( compute_volume ) close(unit=15)
  if ( compute_volume ) close(unit=16)
  if ( manu_sol ) close(unit=14)
  if ( manu_sol ) close(unit=18)


! delete all data including all allocated memory

! flow problem

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh, mesh_0 )
  call delete ( sol )
  call delete ( rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( ssvel_b )
  if ( manu_sol ) call delete ( press )
  call delete ( oldvectors )
  if ( iter_solv ) call delete ( ilu )

! monitor function

  if ( monitor_func ) then

    deallocate ( f_monitor )

!   curvature problem

    call delete ( problem_curv )
    call delete ( sysmatrix_curv )

  end if

! interface tangential motion problem

  if ( monitor_func .or. .not. lagrange ) then
    call delete ( problem_pois )
    call delete ( sysmatrix_pois )
    call delete ( sol_pois )
  end if

! interface tracking problem

  call delete ( coefficients_b, coefficients_i )
  call delete ( mesh_b, mesh_b0 )
  call delete ( input_probdef_b )
  call delete ( problem_b )
  call delete ( sol_b )
  call delete ( rhsd_b )
  call delete ( sysmatrix_b )
  call delete ( vel_b, xnp1, xn )
  if ( timeint == 2 ) call delete ( xnm1 )
  if ( compute_error ) call delete ( error, tauerr, tause )
  call delete ( stress_ten, stress_scal )
  call delete ( lambda, eigv )
  deallocate ( disp )

! ALE mesh motion problem

  if ( numtimesteps > 1 ) call delete ( problem_update_mesh )


contains


! generate and read mesh

  subroutine generate_and_read_mesh

!   generate external and internal mesh

    call write_gmsh_parameters ( ox, oy, oz, lx, ly, lz, xp, rp, &
      dx_box, dx_part )

    call execute_command_line ( 'gmsh -3 -format msh2 -o mesh_conform.msh &
      &mesh_conform.geo > outputmesh_conform.out' )

!   read mesh generated by gmsh

    call read_mesh_gmsh ( mesh_conform, filename='mesh_conform.msh', &
      physgeom=.true. )

!   split mesh in external and internal mesh

    call mesh_convert ( mesh_conform, mesh_ext, only_groups=[1] )
    call mesh_convert ( mesh_conform, mesh_int, only_groups=[2] )

    call delete ( mesh_conform )

!   merge internal and external meshes

    call mesh_merge ( mesh_ext, mesh_int, mesh, nogroupmerge=.true. )

    call delete ( mesh_ext, mesh_int )

    if ( periodic ) then

!     local node numbering surface 4 to match local node numbering of surface 2

      call add_to_mesh ( mesh, matchingsurface=[4,2], replace=4, &
        displacement=[-lx,0._dp,0._dp] )

!     local node numbering surface 6 to match local node numbering of surface 5

      call add_to_mesh ( mesh, matchingsurface=[6,5], replace=6, &
        displacement=[0._dp,0._dp,-lz] )

    else if ( triperiodic ) then

!     local node numbering surface 4 to match local node numbering of surface 2

      call add_to_mesh ( mesh, matchingsurface=[4,2], replace=4, &
        displacement=[-lx,0._dp,0._dp] )

!     local node numbering surface 6 to match local node numbering of surface 5

      call add_to_mesh ( mesh, matchingsurface=[6,5], replace=6, &
        displacement=[0._dp,0._dp,-lz] )

!     local node numbering surface 3 to match local node numbering of surface 1

      call add_to_mesh ( mesh, matchingsurface=[1,3], replace=1, &
        displacement=[0._dp,-ly,0._dp] )

    else if ( periodic_channel ) then

!     local node numbering surface 4 to match local node numbering of surface 2

      call add_to_mesh ( mesh, matchingsurface=[4,2], replace=4, &
        displacement=[-lx,0._dp,0._dp] )

    end if

!   Change of local node numbering of surface 8 to match local node numbering of
!   surface 7

!    call add_to_mesh ( mesh, matchingsurface=[8,7], replace=8 )

  end subroutine generate_and_read_mesh


  subroutine write_gmsh_parameters ( ox, oy, oz, lx, ly, lz, xp, rp, &
    dx_box, dx_part )

    integer :: i, nobj
    real(dp), intent(in) :: ox, oy, oz, lx, ly, lz, xp(:,:), rp(:), &
      dx_box, dx_part

    character(len=*), parameter :: fmti = '(1x,a,i0,a)'
    character(len=*), parameter :: fmtr = '(1x,a,f18.14,a)'
    character(len=*), parameter :: fmtir = '(1x,a,i0,a,f18.14,a)'

    nobj = size(rp)

    open ( unit=25, file='mesh_conform.geo' )

    write ( 25, fmtr ) 'ox = ', ox, ';'
    write ( 25, fmtr ) 'oy = ', oy, ';'
    write ( 25, fmtr ) 'oz = ', oz, ';'
    write ( 25, fmtr ) 'lx = ', lx, ';'
    write ( 25, fmtr ) 'ly = ', ly, ';'
    write ( 25, fmtr ) 'lz = ', lz, ';'
    write ( 25, fmtr ) 'dx_box = ', dx_box, ';'
    write ( 25, fmti ) 'nobj = ', nobj, ';'

    do i = 1, nobj
      write ( 25, fmtir ) 'xp[', i, '] = ', xp(i,1), ';'
      write ( 25, fmtir ) 'yp[', i, '] = ', xp(i,2), ';'
      write ( 25, fmtir ) 'zp[', i, '] = ', xp(i,3), ';'
      write ( 25, fmtir ) 'rp[', i, '] = ', rp(i), ';'
    end do

    write ( 25, fmtr ) 'dx_part = ', dx_part, ';'

    write ( 25, '(/1x,a)' ) 'Include "mesh_options_3D.igo";'
    write ( 25, '(/1x,a)' ) 'Include "particles_in_a_box_3D_full.igo";'

    close ( 25 )

  end subroutine write_gmsh_parameters


! build and solve for velocity and pressure

  subroutine build_solve_velocity_pressure

    integer :: nu

!   fill boundary conditions

    if ( triperiodic ) then
      call fill_sysvector ( mesh, problem, sol, point=1, physq=physqvel, &
        value=0._dp )
    else
      call fill_sysvector ( mesh, problem, sol, surface1=1, surface2=6, &
        physq=physqvel, vfunc=vfunc_3D, vfuncnr=vfuncnr )
    end if

    call fill_sysvector ( mesh, problem, sol, point=1, physq=physqpress, &
      value=0._dp )

!   Stokes velocity/pressure in external and internal fluid

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, mcoefficients=coefficients )

    if ( periodic ) then

!     periodical condition on velocities in x- and z-direction

      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        constraint1=c_perx, constraint2=c_perz, &
        elemsub=stokes_constr_node_conn, addmatvec=.true. )

    else if ( triperiodic ) then

!     periodical condition on velocities

      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        constraint1=c_perx, constraint2=c_perz, &
        elemsub=stokes_constr_node_conn, addmatvec=.true. )

      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        constraint1=c_pery, elemsub=stokes_constr_node_conn_shear, &
        addmatvec=.true. )

    else if ( periodic_channel ) then

!     periodical condition on velocities in x-direction

      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        constraint1=c_perx, elemsub=stokes_constr_node_conn, addmatvec=.true. )

    end if

    if ( flowrate .and. ( periodic .or. periodic_channel ) ) then

!     impose flowrate

      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        constraint1=c_flr, elemsub=stokes_constr_flowr, addmatvec=.true., &
        coefficients=coefficients(1) )

    end if

!   Coupling of fluids

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=c_cpl, elemsub=stokes_constr_node_conn, addmatvec=.true. )

!   Surface tension and interfacial rheology

    call add_boundary_elements ( mesh, problem, rhsd, &
      elemsub=interfacial_rheology_surface, surface=7, &
      sysmatrix=sysmatrix, coefficients=coefficients_i, &
      oldvectors=oldvectors, physq=[physqvel] )

!   Cancel jump in traction for method of manufactured solutions
    if ( manu_sol ) then
      time = time0 + step * deltat
      call add_boundary_elements ( mesh, problem, rhsd, &
        elemsub=interfacial_rheology_surface_mansol, surface=7, &
        coefficients=coefficients_i, physq=[physqvel] )
    end if

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    if ( iter_solv ) then

!     find the 2-norm of the RHS (which will be used to manually set a relative
!     tolerance)
      nu = problem%numundegfd
      normf = sqrt(dot_product(rhsd%u(:nu),rhsd%u(:nu)))

!     Generate permutation arrays with metis

      if ( .not. sysmatrix%renumber ) call renumber_metis ( sysmatrix )

!     set parameters for iterative solver

      call set_solver_options ( solver_options, printlevel=1, type_prec=2, &
        maxmvm=300, itsolver=8, mgmres=30, prec_store=3._dp, preconditioner=1, &
        droptol=1e-6_dp, fillin=3._dp, eps_rel=0._dp, eps_abs=eps_sk*normf, &
        use_renumber=.true. )

      if ( step == 1 ) then
        call toc('build preconditioner')
        solver_options%stoponmaxmvm = .true.
        call solve_system_sk ( sysmatrix, rhsd, sol, ilu, &
          solver_options=solver_options )
        call toc('done building preconditioner')
      else
        solver_options%stoponmaxmvm = .false.
        call solve_system_sk ( sysmatrix, rhsd, sol, ilu, initsol=.true., &
          solver_options=solver_options, maxmvmreached=maxmvmreached )
        if ( maxmvmreached ) then
          call toc('build preconditioner')
          call delete_ilu_sk ( ilu )
          solver_options%stoponmaxmvm = .true.
          call solve_system_sk ( sysmatrix, rhsd, sol, ilu, initsol=.true., &
            solver_options=solver_options )
          call toc('done building preconditioner')
        end if
      end if

    else

      solver_options_ma57%real_storage=rs_up
      solver_options_ma57%integer_storage=is_up

      call solve_system_ma57 ( sysmatrix, rhsd, sol, &
        solver_options=solver_options_ma57  )

    end if

  end subroutine build_solve_velocity_pressure


! build and solve for xnp1 on interface

  subroutine build_solve_interface

    integer :: i

!   Velocity of the interface nodes (reshaped)
    real(dp), dimension(ndim,mesh_b%nnodes) :: drop_vel

    real(dp) :: dx_min, dx_max

    if ( monitor_func ) then

!     Curvature problem on interface mesh

      call curvature_gd_ma57 ( mesh_b, problem_curv, coefficients_b, &
        sysmatrix_curv, curvature=f_monitor, square=.true. )

      f_monitor = 4 * dx_part ** 2 / rc ** 2 / f_monitor

      dx_min = fac_min * dx_part
      dx_max = fac_max * dx_part

      f_monitor = min ( dx_max**2, max(f_monitor,dx_min**2) )

!     Poisson problem on interface mesh for making nodes "diffuse" and
!     redistribute to highly curved regions

      call poisson_gd_ma57 ( mesh_b, problem_pois, coefficients_b, &
        sysmatrix_pois, sol_pois, f_monitor=f_monitor )

    else if ( .not. lagrange ) then

!     Poisson problem on interface mesh for making nodes "diffuse" and
!     redistribute uniformly

      call poisson_gd_ma57 ( mesh_b, problem_pois, coefficients_b, &
        sysmatrix_pois, sol_pois )

    end if

    oldvectors_b%s(1)%p => sol_pois
    oldvectors_b%p(1)%p => problem_pois

    vel_b%u=sol%u(ssvel_b%s) ! extract velocities from velocity solution

    if ( choice_um == 1 ) then

!     average nodal velocities

      drop_vel = reshape(vel_b%u,[ndim,mesh_b%nnodes]) ! reshape
      do i = 1, ndim
        coefficients_b%r(3+i) = sum(drop_vel(i,:))/mesh_b%nnodes
      end do

    else if ( choice_um == 2 ) then

!     velocity of the center of volume

      coefficients_b%r(4:3+ndim) = cvv

    end if

    if ( lagrange ) then

      if ( timeint == 1 .or. &
               ( euler_step2 .and. timeint ==2 .and. restart /= 1 ) ) then

        xnp1%u = xn%u + vel_b%u * deltat

        if ( timeint == 2 ) euler_step2 = .false.

      else if ( timeint == 2 .or. &
                    ( step == 2 .and. timeint > 2 .and. restart /= 1 ) ) then

        xnp1%u = 4._dp/3._dp * xn%u - 1._dp/3._dp * xnm1%u + &
          2._dp/3._dp * vel_b%u * deltat

      end if

    else

!     build (assemble) matrix and vector from elements

      call build_system ( mesh_b, problem_b, sysmatrix_b, msysvector=rhsd_b, &
        elemsub=interface_tracking_elem_lagrange, coefficients=coefficients_b, &
        oldvectors=oldvectors_b )

!     Generate permutation arrays with metis

      if ( .not. sysmatrix_b%renumber ) call renumber_metis ( sysmatrix_b )

      do i = 1, ndim

        so%integer_storage=1.8_dp
        so%pivot_order=1 ! use pivot order from the permutation arrays

        call solve_system_ma41 ( sysmatrix_b, rhsd_b(i), sol_b, lu=lu_b, &
          solver_options=so )

        call transfer_data ( mesh_b, problem1=problem_b, &
          sysvector1=sol_b, vector2=xnp1, degfd1=[1], degfd2=[i] )

      end do

      call delete(lu_b)

    end if

  end subroutine build_solve_interface


! solve flow and move the interface position (time discretized step)

  subroutine solve_flow_and_move_interface

    if ( timeint == 1 .or. &
               ( euler_step .and. timeint ==2 .and. restart /= 1 ) ) then

      if ( timeint == 2 ) euler_step = .false.

!     Euler

      call copy ( xnp1, xn )

      coefficients_b%i(11) = 1

!     update the nodes of the mesh

      call update_mesh_nodes_from_displacement ( order=1 )

    else if ( timeint == 2 .or. &
                    ( step == 2 .and. timeint > 2 .and. restart /= 1 ) ) then

!     second order Gear with prediction

      call copy ( xn, xnm1 )
      call copy ( xnp1, xn )

!     2nd order prediction of interface coordinates
      mesh_b%coor = &
                 transpose( reshape( 2*xn%u - xnm1%u, [ndim,mesh_b%nnodes] ) )

      coefficients_b%i(11) = 2

!     update the nodes of the mesh

      call update_mesh_nodes_from_displacement ( order=2 )

    else

      write (*,'(/a,i0/)') 'Error: invalid timeint = ', timeint
      stop

    end if

!   solve velocity and pressure

    call build_solve_velocity_pressure

!   solve for interface velocity

    call build_solve_interface

!   Updating of interface mesh coordinates with new xnp1

    mesh_b%coor = transpose( reshape( xnp1%u, [ndim,mesh_b%nnodes] ) )

  end subroutine solve_flow_and_move_interface


! update the nodes of the mesh and compute the mesh velocity

  subroutine update_mesh_nodes_from_displacement ( order )

    integer, intent(in) :: order

!   displacement of the predicted interface with respect to the previous
!   mesh coordinates

    disp = mesh_b%coor - mesh%coor(mesh%surfaces(7)%nodes,:)

    if ( move_mesh_x .or. choice_um == 2 ) then

!     update center of volume

      cvnm1 = cvn
      cvn = cv

      call integrate ( mesh_b, problem_b, resultsum, &
        elemsub=center_of_volume, coefficients=coefficients_b )

      cv = resultsum(2:)/resultsum(1)

      if ( order == 1 ) then
!       BDF1
        cvv = ( cv - cvn  ) / deltat
      else
!       BDF2
        cvv = ( 1.5_dp*cv - 2*cvn + 0.5_dp*cvnm1 ) / deltat
      end if

    end if

    if ( move_mesh_x ) then

!     additional ALE movement of the mesh in x-direction

      disp(:,1) = disp(:,1) - ( cv(1) - xc(1) )

    end if

!   update mesh nodes

    call update_mesh_nodes ( mesh, problem_update_mesh, disp=disp )

  end subroutine update_mesh_nodes_from_displacement


! output various data

  subroutine output_data

    use output_fields1_m

    write(*,'(a,i0,a,es10.3)') &
      'Output data at step ', step0+step, ', time =', time0+step*deltat

    write(filename,'(a,i4.4,a)') 'output', post, '.vtk'

    if ( compute_error ) then
!     error
      call write_scalar_vtk ( mesh_b, problem_b, filename=filename, &
        dataname='error', vector=error )
    end if

!   surface divergence of velocity
    coefficients_i%i(13) = 1
    call derive_vector ( mesh_b, problem_b, stress_scal, &
      elemsub=interfacial_stress_scal, coefficients=coefficients_i, &
      oldvectors=oldvectors_b )
    call write_scalar_vtk ( mesh_b, problem_b, filename=filename, &
      dataname='divsu', vector=stress_scal, append=.false. )

!   rate of strain gammadot = sqrt(2D_s:D_s)
    coefficients_i%i(13) = 2
    call derive_vector ( mesh_b, problem_b, stress_scal, &
      elemsub=interfacial_stress_scal, coefficients=coefficients_i, &
      oldvectors=oldvectors_b )
    call write_scalar_vtk ( mesh_b, problem_b, filename=filename, &
      dataname='gammadot', vector=stress_scal, append=.true. )

!   surface vorticity = omega . normal vector
    coefficients_i%i(13) = 3
    call derive_vector ( mesh_b, problem_b, stress_scal, &
      elemsub=interfacial_stress_scal, coefficients=coefficients_i, &
      oldvectors=oldvectors_b )
    call write_scalar_vtk ( mesh_b, problem_b, filename=filename, &
      dataname='vorticity', vector=stress_scal, append=.true. )

!   surface determinant of F_s
    coefficients_i%i(13) = 4
    call derive_vector ( mesh_b, problem_b, stress_scal, &
      elemsub=interfacial_stress_scal, coefficients=coefficients_i, &
      oldvectors=oldvectors_b )
    call write_scalar_vtk ( mesh_b, problem_b, filename=filename, &
      dataname='detFs', vector=stress_scal, append=.true. )

!   derive interfacial stress tensor
    call derive_vector ( mesh_b, problem_b, stress_ten, &
      elemsub=interfacial_stress_ten, coefficients=coefficients_i, &
      oldvectors=oldvectors_b )

!   interfacial stress error
    if ( compute_error ) then

      call derive_vector ( mesh_b, problem_b, tause, &
        elemsub=interfacial_stress_ten_mansol, coefficients=coefficients_i, &
        oldvectors=oldvectors_b )

      tauerr%u = tause%u - stress_ten%u

      write(unit=17,fmt=*) time0 + step * deltat, &
                           maxval(sqrt(sum(reshape(tauerr%u**2, &
                           [ndim*(ndim+1)/2,mesh_b%nnodes]),1))), &
                           sqrt(sum(tauerr%u**2)/mesh_b%nnodes)

    end if

!   principal interfacial stresses and stress vectors

    call eig ( stress_ten%u, lambda(1)%u, lambda(2)%u, eigv(1)%u, eigv(2)%u )

    call write_scalar_vtk ( mesh_b, problem_b, filename=filename, &
      dataname='sigmas1', vector=lambda(1), append=.true. )
    call write_vector_vtk ( mesh_b, problem_b, filename=filename, &
      dataname='stressvec1', vector=eigv(1), append=.true. )
    call write_scalar_vtk ( mesh_b, problem_b, filename=filename, &
      dataname='sigmas2', vector=lambda(2), append=.true. )
    call write_vector_vtk ( mesh_b, problem_b, filename=filename, &
      dataname='stressvec2', vector=eigv(2), append=.true. )

!   principal interfacial strains and strain directions

    call derive_vector ( mesh_b, problem_b, stress_ten, &
      elemsub=interfacial_Bs_ten, coefficients=coefficients_i, &
      oldvectors=oldvectors_b )

    call eig ( stress_ten%u, lambda(1)%u, lambda(2)%u, eigv(1)%u, eigv(2)%u )
    lambda(1)%u = sqrt(lambda(1)%u) - 1._dp
    lambda(2)%u = sqrt(lambda(2)%u) - 1._dp

    call write_scalar_vtk ( mesh_b, problem_b, filename=filename, &
      dataname='epsilon1', vector=lambda(1), append=.true. )
    call write_vector_vtk ( mesh_b, problem_b, filename=filename, &
      dataname='strainvec1', vector=eigv(1), append=.true. )
    call write_scalar_vtk ( mesh_b, problem_b, filename=filename, &
      dataname='epsilon2', vector=lambda(2), append=.true. )
    call write_vector_vtk ( mesh_b, problem_b, filename=filename, &
      dataname='strainvec2', vector=eigv(2), append=.true. )

!   output coordinates
    write(filename,'(a,i4.4,a)') 'coor', post, '.out'
    open ( unit=11, file=filename, recl=300 )
    write(unit=11,fmt='(i0,3e16.8)') &
                           (i, mesh_b%coor(i,:), i=1,mesh_b%nnodes)
    close ( unit=11 )

!   output fields to vtk files

    call output_fields ( mesh, problem, coefficients, sol, post )

  end subroutine output_data


! output restart data

  subroutine output_restart

    write(*,'(a,i0,a,es10.3)') &
      'Output restart data at step ', step0+step, ', time =', time0+step*deltat

    call write_mesh ( mesh, filename='mesh.out' )
    call write_mesh ( mesh_0, filename='mesh_0.out' )

    open ( unit=11, form='unformatted', file='data.out' )

    write(11) step0 + step, post, time0 + step*deltat
    write(11) mesh_b%coor
    write(11) xn%u, xnp1%u
    write(11) cv, cvn, cvnm1
    w1 = timeint == 2
    write(11) w1
    if ( w1 ) write(11) xnm1%u

    close(unit=11)

  end subroutine output_restart


! open files and input restart data

  subroutine open_files_and_restart

    integer :: endstep, endpost
    real(dp) :: endtime

    if ( restart >= 1 ) then

!   Restart: initial mesh

      call delete ( mesh_0 )
      call read_mesh ( mesh_0, 'mesh_0.out' )

!     restart: read solution from file

      open ( unit=11, form='unformatted', file='data.out' )

      read(11) endstep, endpost, endtime
      read(11) mesh_b%coor
      read(11) xn%u, xnp1%u
      read(11) cv, cvn, cvnm1
      read(11) w1
      if ( w1 ) read(11) xnm1%u

      close(unit=11)

      if ( ctime ) then
        time0 = endtime
        post0 = endpost
        step0 = endstep
      end if

      if ( compute_error ) then

!       open error file

        if ( ctime ) then
          open ( unit=10, file='error.out', recl=300, status='old', &
            position='append' )
          open ( unit=17, file='stress_error.out', recl=300, status='old', &
            position='append' )
        else
          open ( unit=10, file='error.out', recl=300, status='unknown' )
          open ( unit=17, file='stress_error.out', recl=300, status='unknown' )
        end if

      end if

      if ( compute_volume ) then

!       open volume, axes and inertia tensor file

        if ( ctime ) then
          open ( unit=12, file='volume.out', recl=300, status='old', &
            position='append' )
          open ( unit=15, file='axes.out', recl=300, status='old', &
            position='append' )
          open ( unit=16, file='inertia_ten.out', recl=300, status='old', &
            position='append' )
        else
          open ( unit=12, file='volume.out', recl=300, status='unknown' )
          open ( unit=15, file='axes.out', recl=300, status='unknown' )
          open ( unit=16, file='inertia_ten.out', recl=300, status='unknown' )
        end if

      end if

!     open L2-norm of pressure and velocity file

      if ( manu_sol ) then
        if ( ctime ) then
          open ( unit=18, file='L2norm.out', recl=300, status='old', &
            position='append' )
        else
          open ( unit=18, file='L2norm.out', recl=300, status='unknown' )
        end if

!       open min and max pressure file

        if ( ctime ) then
          open ( unit=14, file='minmaxpres.out', recl=300, status='old', &
            position='append' )
        else
          open ( unit=14, file='minmaxpres.out', recl=300, status='unknown' )
        end if
      end if

    else

!     fresh start

      if ( compute_error ) then

!       open error file

        open ( unit=10, file='error.out', recl=300 )
        open ( unit=17, file='stress_error.out', recl=300 )

      end if

      if ( compute_volume ) then

!       open volume, axes and inertia tensor file

        open ( unit=12, file='volume.out', recl=300 )
        open ( unit=15, file='axes.out', recl=300 )
        open ( unit=16, file='inertia_ten.out', recl=300 )

      end if

!     open L2-norm of pressure and velocity file

      if ( manu_sol ) then
        open ( unit=18, file='L2norm.out', recl=300 )

!       open min and max pressure file

        open ( unit=14, file='minmaxpres.out', recl=300 )
      end if

    end if

  end subroutine open_files_and_restart


! Problem definition of velocity/pressure system, sysmatrix and vectors

  subroutine problem_definition_create_sysmatrix_vectors

  call problem_definition ( input_probdef, mesh, problem )

! Creation of system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )

! Creation of system matrix for velocity/pressure problem

  if ( iter_solv ) then
    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  else
    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric=.true. )
  end if
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! subscript for extracting velocities from mesh nodes on interface

  call create_subscript ( mesh, problem, ssvel_b, surfaces=[7], &
    physqarr=[physqvel] )

! subscript for extracting pressures

  if ( manu_sol ) then
    call create_subscript ( mesh, problem, press, physqarr=[physqpress] )
  end if

  end subroutine problem_definition_create_sysmatrix_vectors


! Delete problem of velocity/pressure system, sysmatrix and vectors

  subroutine delete_problem_sysmatrix_vectors

    call delete ( problem )
    call delete ( sol, rhsd )
    call delete ( sysmatrix )
    call delete ( ssvel_b )
    if ( manu_sol ) call delete ( press )

  end subroutine delete_problem_sysmatrix_vectors


! Create and initialize the aspect ratio arrays needed for remesh criterion

  subroutine create_aspect_ratio

!   allocate aspect ratio arrays
    allocate ( init_volume(mesh%nelem), volumev(mesh%nelem) )
    allocate ( init_aspect_ratio(mesh%nelem), aspect_ratio(mesh%nelem) )

!   compute initial element aspect ratio
    call compute_element_aspect_ratio ( mesh, init_volume, init_aspect_ratio )

  end subroutine create_aspect_ratio


! Delete the aspect ratio arrays

  subroutine delete_aspect_ratio

!   deallocate aspect ratio arrays
    deallocate ( init_volume, volumev )
    deallocate ( init_aspect_ratio, aspect_ratio )

  end subroutine delete_aspect_ratio


! compute the aspect ratio and volume for each mesh element

  subroutine compute_element_aspect_ratio ( mesh, volume, asp_ratio )

    type(mesh_t), intent(in) :: mesh
    real(dp), dimension(:), intent(out) :: asp_ratio, volume

    integer :: elem, node(4), grp, totelem
    real(dp) :: la, lb, lc, ld, le, lf
    real(dp) :: vert(4,3)

    node = [1,3,5,10]

    totelem = 0

    do grp = 1, mesh%nelgrp
      do elem = 1,mesh%grpnumel(grp)

        totelem = totelem + 1

!       get coordinates of tetrahedron vertices
        vert = mesh%coor(mesh%topology(grp)%a(node,elem),:)

!       compute side lengths
        la = lv ( vert(1,:) - vert(2,:) )
        lb = lv ( vert(2,:) - vert(3,:) )
        lc = lv ( vert(3,:) - vert(1,:) )
        ld = lv ( vert(1,:) - vert(4,:) )
        le = lv ( vert(2,:) - vert(4,:) )
        lf = lv ( vert(3,:) - vert(4,:) )

!       compute volume of tetrahedron
        volume(totelem) = abs ( dot_product ( vert(1,:) - vert(4,:), &
          cross_product ( vert(2,:) - vert(4,:), vert(3,:) - vert(4,:) ) ) ) / 6

!       compute aspect ratio
        asp_ratio(totelem) = max(la,lb,lc,ld,le,lf)**3/volume(totelem)

      end do
    end do

  end subroutine compute_element_aspect_ratio


! length function

  function lv ( a )

    real(dp), intent(in), dimension(3) :: a
    real(dp) :: lv

    lv = sqrt ( a(1)**2 + a(2)**2 + a(3)**2 )

  end function lv


! Remesh based on the interface mesh and the outside box. The interface mesh
! remains as it is. Only the volume is remeshed using gmsh.

  subroutine remeshing

    real(dp), dimension(:,:), allocatable :: coor_boun

    character(len=*), parameter :: fmti = '(1x,a,i0,a)'
    integer, parameter :: nobj = 1

    call write_geometries_gmsh ( filename = 'mesh_geometries.msh', mesh=mesh, &
      binary=.true. )

    open ( unit=25, file='mesh_conform.geo' )

    write ( 25, fmti ) 'nobj = ', nobj, ';'
    write ( 25, '(/1x,a)' ) 'Include "mesh_options_3D.igo";'
    write ( 25, '(/1x,a)' ) 'Merge "mesh_geometries.msh";'
    write ( 25, '(/1x,a)' ) 'Include "particles_in_a_box_3D_full_2.igo";'

    close ( 25 )

!   save coordinates of the boundary mesh

    allocate ( coor_boun(mesh%surfaces(7)%nnodes,mesh%ndim) )
    coor_boun = mesh%coor(mesh%surfaces(7)%nodes,:)

    call execute_command_line ( 'gmsh -3 -bin -format msh2 -o &
      &mesh_conform.msh mesh_conform.geo > outputmesh_conform.out' )

!   read mesh generated by gmsh

    call read_mesh_gmsh ( mesh_conform, filename='mesh_conform.msh', &
      physgeom=.true. )

!   match numbering/topology of surface=7 with old surface in mesh, since
!   after remeshing with gmsh this has been changed. Use minimum distance
!   algorithm since coordinates have also been changed.

    call add_to_mesh ( mesh_conform, matchingsurface=[7,7], replace=7, &
      matchingmesh=mesh, mindistance=.true. )

!   split mesh in external and internal mesh

    call mesh_convert ( mesh_conform, mesh_ext, only_groups=[1] )
    call mesh_convert ( mesh_conform, mesh_int, only_groups=[2] )

    call delete ( mesh_conform, mesh )

!   merge internal and external meshes

    call mesh_merge ( mesh_ext, mesh_int, mesh, nogroupmerge=.true. )

    call delete ( mesh_ext, mesh_int )

    if ( periodic ) then

!     local node numbering surface 4 to match local node numbering of surface 2

      call add_to_mesh ( mesh, matchingsurface=[4,2], replace=4, &
        displacement=[-lx,0._dp,0._dp] )

!     local node numbering surface 6 to match local node numbering of surface 5

      call add_to_mesh ( mesh, matchingsurface=[6,5], replace=6, &
        displacement=[0._dp,0._dp,-lz] )

    else if ( triperiodic ) then

!     local node numbering surface 4 to match local node numbering of surface 2

      call add_to_mesh ( mesh, matchingsurface=[4,2], replace=4, &
        displacement=[-lx,0._dp,0._dp] )

!     local node numbering surface 6 to match local node numbering of surface 5

      call add_to_mesh ( mesh, matchingsurface=[6,5], replace=6, &
        displacement=[0._dp,0._dp,-lz] )

!     local node numbering surface 3 to match local node numbering of surface 1

      call add_to_mesh ( mesh, matchingsurface=[1,3], replace=1, &
        displacement=[0._dp,-ly,0._dp] )

    else if ( periodic_channel ) then

!     local node numbering surface 4 to match local node numbering of surface 2

      call add_to_mesh ( mesh, matchingsurface=[4,2], replace=4, &
        displacement=[-lx,0._dp,0._dp] )

    end if

!   Change of local node numbering of surface 8 to match local node numbering of
!   surface 7

!    call add_to_mesh ( mesh, matchingsurface=[8,7], replace=8 )

   call fill_mesh_parts ( mesh )

!  delete various structures, arrays

   call delete_problem_sysmatrix_vectors
   call delete_aspect_ratio
   call delete ( problem_update_mesh )

!  recreate the structures, arrays based on the new mesh

   call problem_definition_create_sysmatrix_vectors
   call create_aspect_ratio

!  move boundary nodes to previous position (gmsh does not preserve quadratic
!  element shape).

   mesh%coor(mesh%surfaces(7)%nodes,:) = coor_boun
   mesh%coor(mesh%surfaces(8)%nodes,:) = coor_boun

   euler_step = .true.
   euler_step2 = .true.

   deallocate ( coor_boun )

 end subroutine remeshing


! Boundary element for surface tension and interfacial rheology on a surface.
! Cancel jump in traction for method of manufactured solutions.

  subroutine interfacial_rheology_surface_mansol ( mesh, problem, surface, &
    elem, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemvec )

    use interfacial_rheology_globals_m
    use functions_mansol_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: surface, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, j, vfuncnr, imodel
    real(dp) :: gammac, muc, kappac, Gc, Kc, err


    if ( first ) then

!     first element on this surface

!     set globals

      call set_globals_interfacial_bc ( mesh, coefficients, &
        surface=surface )

!     allocate arrays

      allocate ( wg(ninti), surfl(ninti) )
      allocate ( tmp(ndf,ndim), work(ninti) )
      allocate ( xig(ninti,2), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim), xe(ninti,ndim) )
      allocate ( dphi(ninti,ndf,2), dxdxis(ninti,ndim,2) )
      allocate ( rhsf(ninti,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

!   get position of integration points

    call get_coordinates_geometry ( mesh, elem, x, surface=surface )

    call isoparametric_deformation_surface ( x, dphi, dxdxis, surfl )

    call isoparametric_coordinates ( x, phi, xg )

!   surface tension coefficient

    gammac = coefficients%r(2)  ! surface tension
    kappac = coefficients%r(3)  ! surface dilatational viscosity
    muc = coefficients%r(4)     ! surface shear viscosity
    Kc = coefficients%r(5)      ! surface dilatational elasticity
    Gc = coefficients%r(6)      ! surface shear elasticity

    vfuncnr = coefficients%i(10)
    imodel = coefficients%i(11)

    if ( matrix ) then
      write(*,'(/3(a/))') 'Error in interfacial_rheology_surface_mansol:', &
        ' matrix == .true. Make sure add_boundary_elements is called ', &
        ' with buildmatrix=.false.'
      stop
    end if

    if ( vector ) then

      do j = 1, ninti
!       determine exact position
        err = func_error_mansol_3D ( vfuncnr, xg(j,:) )
        xe(j,:) = xg(j,:) + err * xg(j,:) / &
          sqrt(dot_product(xg(j,:),xg(j,:)))
!       surface tension term
        rhsf(j,:) = gammac * divsIs3D ( xe(j,:) )
      end do

!     constant surface tension
      do j = 1, ndim
        do N = 1, ndf
          work = phi(:,N) * rhsf(:,j)
          tmp(N,j) = sum ( work * surfl * wg )
        end do
      end do

      select case ( imodel )
      case ( 1,3,5 ) ! Boussinesq-Scriven or Kelvin-Voigt
        do j = 1, ninti
          rhsf(j,:) = muc * divsDs3D ( xe(j,:) ) + &
            ( kappac - muc ) * divsdivsuIs3D ( xe(j,:) )
        end do
        do j = 1, ndim
          do N = 1, ndf
            work = phi(:,N) * rhsf(:,j)
            tmp(N,j) = tmp(N,j) + sum ( work * surfl * wg )
          end do
        end do
      case ( 0,2,4 ) ! Constant surface tension or elastic model
      case default
        write(*,'(a,a,i0)') 'Error in interfacial_rheology_surface_mansol: ', &
          'wrong interface model number = ', imodel
        stop
      end select

      select case ( imodel )
      case ( 0,1 ) ! Constant surface tension or Boussinesq-Scriven
      case ( 2:5 ) ! Elastic or Kelvin-Voigt
        do j = 1, ninti
          rhsf(j,:) = Kc * divstauhyd3D ( imodel, xe(j,:) ) + &
                      Gc * divstaudev3D ( imodel, xe(j,:) )
        end do
        do j = 1, ndim
          do N = 1, ndf
            work = phi(:,N) * rhsf(:,j)
            tmp(N,j) = tmp(N,j) + sum ( work * surfl * wg )
          end do
        end do
      case default
        write(*,'(a,a,i0)') 'Error in interfacial_rheology_surface_mansol: ', &
          'wrong interface model number = ', imodel
        stop
      end select

      elemvec = reshape ( tmp, [ndf*ndim] )

    end if

    if ( last ) then

!     last element on this surface

      deallocate ( wg, surfl )
      deallocate ( xig, phi, x )
      deallocate ( tmp, work )
      deallocate ( xg, xe, rhsf )
      deallocate ( dphi, dxdxis )

    end if

  end subroutine interfacial_rheology_surface_mansol


! Compute interfacial stress in all nodes

  subroutine interfacial_stress_ten_mansol ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec, elemwts )

    use interfacial_rheology_globals_m
    use functions_mansol_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    real(dp) :: gammac, kappac, Kc, muc, Gc, err
    integer :: i, imodel, vfuncnr

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_interfacial_stress ( mesh, coefficients, elgrp )

!     allocate arrays

      allocate ( x(nodalp,ndim), work3(nodalp,2*ndim) )

    end if

!   interface model nr, flow nr and material parameters

    vfuncnr = coefficients%i(10)
    imodel = coefficients%i(11)
    gammac = coefficients%r(2)
    kappac = coefficients%r(3)
    muc = coefficients%r(4)
    Kc = coefficients%r(5)
    Gc = coefficients%r(6)

!   position

    call get_coordinates ( mesh, elgrp, elem, x )

!   interfacial tension
    do i = 1, nodalp
      err = func_error_mansol_3D ( vfuncnr, x(i,:) )
      x(i,:) = x(i,:) + err * x(i,:) / &
        sqrt(dot_product(x(i,:),x(i,:)))
      work3(i,:) = gammac * Is3D ( x(i,:) )
    end do

!   viscous stress
    select case ( imodel )
    case ( 0,2,4 )
    case ( 1,3,5 )
      do i = 1, nodalp
        work3(i,:) = work3(i,:) + &
                     (kappac-muc) * divsuIs3D ( x(i,:) ) + &
                            muc  * Ds3D ( x(i,:) )
      end do
    case default
      print *, 'Erorr in interfacial_stress_ten_mansol: ', &
               'interface model imodel = ', imodel, ' not available.'
      stop
    end select

!   elastic stress
    select case ( imodel )
    case ( 0,1 )
    case ( 2:5 )
      do i = 1, nodalp
        work3(i,:) = work3(i,:) + &
                     Kc * tauhyd3D ( imodel, x(i,:) ) + &
                     Gc * taudev3D ( imodel, x(i,:) )
      end do
    case default
      print *, 'Erorr in interfacial_stress_ten_mansol: ', &
               'interface model imodel = ', imodel, ' not available.'
      stop
    end select

!   interfacial stress tensor

    do i = 1, nodalp
      elemvec(i:i+5*nodalp:nodalp) = work3(i,:)
    end do

    elemwts = 1

    if ( last ) then

!     last element in this group

!     deallocate arrays

      deallocate ( x, work3 )

    end if

  end subroutine interfacial_stress_ten_mansol

end program drop17

#else
  print '(4(a/),a)', &
    'To run this example:', &
    ' - compile add-ons sk_solve and metis5', &
    ' - install and compile sparskit2', &
    ' - add metis5 lib for linking', &
    ' - set preprocessing macro SKIT2 and METIS5 in Mdefs.mk'
end
#endif

