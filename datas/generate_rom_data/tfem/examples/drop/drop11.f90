! Drop consisting of a viscoelastic fluid suspended in a Newtonian fluid.
! An elastic soft object as a limit. No inertia.
! Flow imposed either:
!   1) on the whole boundary using a function, or
!   2) on the upper and lower boundary using a function and
!      periodic boundary conditions in the x-direction (periodic=.true.).
! ALE with interface tracking and optional displacement of the mesh in
! x-direction to keep the center of volume at the same position in the mesh.
! Remeshing is performed when the maximum normalized area or aspect ratio
! of the mesh elements is higher than a threshold.
! Remeshing is based on position of nodes of the interface mesh, but the
! interface mesh itself not modified.
! Optionally move the nodes of the interface mesh to highly curved regions.
! Optionally impose flowrate if periodic=.true.
! Optionally compute the error if the exact solution is known.
! Optionally compute the volume and center of volume.
! Two-dimensional "toy" model for testing purposes.

program drop11

  use tfem_m
  use math_defs_m
  use hsl_ma41_m
  use hsl_ma57_m
  use viscoelastic_elements_m
  use interface_tracking_elements_m
  use poisson_gd_ma57_m
  use curvature_gd_ma57_m
  use update_mesh_nodes1_m
  use functions_m
  use io_utils_m
  use figplot_m
  use projection_elements_m

  implicit none


! constants

  integer, parameter ::  &
    gintpl = 2,          & ! P1 gradients
    uintpl = 6,          & ! P2 velocity and interface position interpolation
    pintpl = 2,          & ! P1 pressure interpolation
    cintpl = 2,          & ! P1 conformation
    physqgrad = 1,       & ! physical quantity nr of gradients
    physqvel = 2,        & ! physical quantity nr of velocity
    physqpress = 3,      & ! physical quantity nr of pressure
    ncompc = 3,          & ! number of conformation tensor components
    nmodes = 1,          & ! number of modes
    startm = 501,        & ! start of material model data
    ndim = 2,            & ! dimension of space
    ncomp_P2 = 6,        & ! number of components for P2 projection:
                           ! x,y-coordinates at n and n-1
    ncomp_P1 = 6           ! number of components for P1 projection:
                           ! conformation components at n and n-1

  real(dp), parameter ::  &
    beta = 1.0_dp          ! upwinding parameter in SUPG method

! variables

  logical ::     &
    solid = .true.,         & ! Internal fluid is a solid (lambda = infty)
    periodic = .false.,     & ! Use periodical boundaries in x direction
                              ! Only makes sense for vfuncnr=0 (zero velocity)
                              ! or vfuncnr=2 (shear flow).
    flowrate = .false.,     & ! Impose flowrate in x-direction (for
                              ! periodic=.true. only)
    move_mesh_x = .true.,   & ! Move the mesh in x-direction to keep the center
                              ! of volume of the drop at the same position wrt
                              ! to the mesh
    monitor_func = .false., & ! Use a monitor function to move the surface
                              ! elements tangential to high curvature regions
    ctime = .false.,        & ! continue time in output after restart
    mesh_plot = .false.,    & ! write mesh plot files/info (curves, mesh etc.)
    compute_error = .false.,& ! Compute error. Only makes sense
                              ! if vfuncnr=0 (no flow) or vfuncnr=3 (rotation).
    compute_volume = .false.  ! Compute volume and center of volume.

  integer ::  &
    restart=0,           & ! restart=0: no restart
                           ! restart=1: normal restart
                           ! restart=2: restart, but start with Euler step
    step0=0,             & ! initial step number
    post0=0,             & ! initial post number
    n=5,                 & ! number of elements in x and y direction of box
    nc=40,               & ! number of elements on the circle
    model = 2,           & ! viscoelastic model: UCM
    gauss  = 4,          & ! order of integration for triangular elements
    gaussb = 4,          & ! number of integration points of line elements
    gauss_proj = 6,      & ! integration rule for the projection problem
    inttype_proj = 3,    & ! standard Gauss-Legendre (numerical table)
    vfuncnr=2,           & ! function number for the velocity field u
    choice_um = 1,       & ! choice for um in the motion of the interface:
                           !  0: um = 0
                           !  1: um = average velocity of all interface nodes
                           !  2: um = velocity of the center of volume
    vtkevery = 5 ,       & ! vtk file every vtkevery steps. -1: means none
    restart_every = 100, & ! write restart file every restart_every steps.
                           ! -1: means none
    logc = 1,            & ! log representation of c?
    timeint1 = 1,        & ! (first-order) time integration for viscoelastic
    timeint2 = 7,        & ! (second-order) time integration for viscoelastic
    timeint=2,           & ! time integration scheme
    numtimesteps=400       ! number of time steps

  real(dp) ::  &
    ox = -0.5_dp,         & ! x-coordinate left lower corner
    oy = -0.5_dp,         & ! y-coordinate left lower corner
    lx = 1.0_dp,          & ! length of the domain
    ly = 1.0_dp,          & ! height of the domain
    xc(ndim) = [0.0_dp,0.0_dp],  & ! center of the initial circle
    rc = 0.1_dp,          & ! initial radius of the circle
    eta_e = 1.0_dp,       & ! external fluid viscosity
    eta_i = 0.0_dp,       & ! internal fluid viscosity (solvent)
    eta_p = 1.0_dp,       & ! internal fluid viscosity (polymer)
    lambda_p = 1.0e-1_dp, & ! internal fluid relaxation time
    G_s = 2._dp,          & ! elastic solid modulus
    lambda_s = 1.0e12_dp, & ! "solid" relaxation rime
    gammac = 0.5_dp,      & ! surface tension coefficient
    factor = 100._dp,     & ! factor for the interface grid deformation
    rs_gup  = 1.4_dp,     & ! real_storage gradient-velocity-pressure LU HSL
    is_gup  = 1.5_dp,     & ! integer_storage gradient-velocity-pressure LU HSL
    rs_c = 1.4_dp,        & ! real_storage conformation LU (HSL)
    is_c = 1.6_dp,        & ! integer_storage conformation LU (HSL)
    U_avg  = 1.0_dp,      & ! imposed average velocity for flowrate=.true.
    time0  = 0.0_dp,      & ! initial time
    betai = 0.5_dp,       & ! factor beta in SUPG
    deltat = 7.e-3_dp,    & ! time step
    fac_min = 0.1_dp,     & ! factor setting dx_min for the monitor function
    fac_max = 10.0_dp,    & ! factor setting dx_max for the monitor function
    area_threshold = 1.39_dp,        & ! remeshing area threshold
    aspect_ratio_threshold = 1.39_dp   ! remeshing aspect ratio threshold

! flow problem

  type(mesh_t) :: mesh_ext, mesh_int, mesh, mesh_conform
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, solm1
  type(sysvector_t) :: rhsd
  type(coefficients_t) :: coefficients(2)
  type(plot_options_t) :: plot_options
  type(solver_options_ma41_t) :: solver_options_u
  type(subscript_t) :: ssvel_b
  type(vector_t), target :: meshvel
  type(subscriptvec_t) :: ssmeshvelx

! conformation problem (CE)

  type(input_probdef_t) :: input_probdef_c
  type(problem_t), target :: problem_c
  type(sysmatrix_t) :: sysmatrix_c
  type(sysvector_t), dimension(ncompc, nmodes), target :: sol_c, sol_cm1
  type(sysvector_t), dimension(ncompc, nmodes) :: rhs_c
  type(oldvectors_t) :: oldvectors_c
  type(lu_ma41_t) :: lu_c
  type(solver_options_ma41_t) :: solver_options_c

! conformation projection

  type(input_probdef_t) :: input_probdef_c_proj
  type(problem_t), target :: problem_c_proj
  type(sysmatrix_t) :: sysmatrix_c_proj
  type(sysvector_t), dimension(ncompc,nmodes), target :: sol_c_proj
  type(sysvector_t), dimension(ncompc,nmodes) :: rhs_c_proj
  type(lu_ma57_t) :: lu_exps_proj

! curvature problem

  type(problem_t) :: problem_curv
  type(sysmatrix_t) :: sysmatrix_curv

! interface tangential motion problem

  type(problem_t), target :: problem_pois
  type(sysmatrix_t) :: sysmatrix_pois
  type(sysvector_t), target :: sol_pois

! interface tracking problem

  type(mesh_t) :: mesh_b
  type(input_probdef_t) :: input_probdef_b
  type(problem_t), target :: problem_b
  type(sysmatrix_t) :: sysmatrix_b
  type(sysvector_t), target :: sol_b
  type(sysvector_t) :: rhsd_b(ndim)
  type(vector_t) :: xnp1, error
! Interface mesh coordinates at previous times. NOTE: the coordinates are
! defined with respect to the initial coordinate system.
  type(vector_t), target :: xn, xnm1
  type(vector_t), target :: vel_b
  type(oldvectors_t) :: oldvectors_b
  type(coefficients_t) :: coefficients_b
  type(lu_ma41_t) :: lu_b
  type(solver_options_ma41_t) :: so

! ALE mesh motion problem

  type(problem_t), target ::  problem_update_mesh
! Mesh coordinates at previous times
  real(dp), allocatable, dimension(:,:) :: meshcoor_n, meshcoor_nm1
! Displacement of interface
  real(dp), allocatable, dimension(:,:) :: disp

! Temporary definitions used to store old definitions before projection

  type(mesh_t), target :: mesh_temp
  type(sysvector_t) :: sol_temp
  type(sysvector_t), dimension(ncompc,nmodes) :: sol_c_temp, sol_cm1_temp
  type(vector_t) :: meshcoor_n_temp, meshcoor_nm1_temp
  type(problem_t), target :: problem_temp, problem_c_temp

! some more misc variables

  logical :: w1, euler_step
  character(len=30) :: filename

  integer :: nnodes, nelem, elshape, step, post, c_cpl, c_per, c_flr
  integer :: vertices(3)=[1,3,5]
  real(dp) :: xp(1,ndim), rp(1), dx_box, dx_part, alpha, G, lambda, &
    resultsum(1+ndim)

! deformation variables for remeshing
  real(dp), dimension(:), allocatable :: init_area, areav
  real(dp), dimension(:), allocatable :: init_aspect_ratio, aspect_ratio
  real(dp) :: norm_area, norm_aspect_ratio
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

  namelist /comppar/ solid, periodic, flowrate, move_mesh_x, monitor_func, &
    mesh_plot, compute_error, compute_volume, n, nc, model, gauss, gaussb, &
    gauss_proj, inttype_proj, &
    choice_um, vfuncnr, vtkevery, ctime, restart, restart_every, logc, &
    timeint1, timeint2, timeint, numtimesteps, &
    ox, oy, lx, ly, xc, rc, fac_min, fac_max, eta_e, eta_i, eta_p, lambda_p, &
    G_s, lambda_s, gammac, factor, U_avg, rs_gup, is_gup, rs_c, is_c, deltat, &
    betai, area_threshold, aspect_ratio_threshold

  read ( unit=*, nml=comppar )


! Fill coefficients

! G and DEVSS parameter

  if ( solid ) then
    G = G_s
    alpha = G_s * deltat
    lambda = lambda_s
  else
    G = eta_p / lambda_p
    alpha = eta_p
    lambda = lambda_p
  end if

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
  coefficients(1)%r(6) = U_avg * ly
  coefficients(1)%r(19) = gammac

! Internal fluid (element group 2)

  coefficients(2)%i = 0
  coefficients(2)%i(1:22) = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss, &
      gaussb,   cintpl,     7,     0,         0, &
      0,        0,          model, nmodes,    startm, &
      logc,     timeint1 ]
  coefficients(2)%i(40) = 3 ! standard Gauss-Legendre (numerical tables)
  coefficients(2)%i(48) = 1 ! use mesh velocity for ALE formulation
  coefficients(2)%i(49) = 1 ! exp(s) projection is used for logc=1

  coefficients(2)%r = 0
  coefficients(2)%r(1:9) = &
    [ eta_i, 0._dp, 0._dp,  alpha, 0._dp, &
      0._dp, 0._dp, deltat, beta ]
  coefficients(2)%r(501:502) = [ G, lambda ]

! Interface

  call create_coefficients ( coefficients_b, ncoefi=100, ncoefr=50 )

  coefficients_b%i = 0
  coefficients_b%i(1) = uintpl
  coefficients_b%i(3) = gaussb
  coefficients_b%i(5) = 3  ! standard Gauss-Legendre (numerical tables)
  coefficients_b%i(6) = 2  ! velocity u on the interface given by a vector
  coefficients_b%i(10) = 1  ! tangential grid velocity using Poisson problem
  coefficients_b%i(11) = timeint

  coefficients_b%r = 0
  coefficients_b%r(7) = factor
  coefficients_b%r(8) = betai
  coefficients_b%r(9) = deltat

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
  call find_bounds_blocks ( mesh )

  if ( mesh_plot ) then

    call plot_points_curves ( plot_options, mesh, filename='curves.fig' )
    call plot_mesh ( plot_options, mesh, filename='mesh.fig' )

    call write_mesh_vtk ( mesh, filename='mesh.vtk' )

    call printinfo ( mesh, printlevel=4 )

  end if

! generate interface mesh from the fluid mesh

  nnodes = mesh%curves(5)%nnodes
  nelem = mesh%curves(5)%nelem
  elshape = mesh%curves(5)%element%elshape

  call mesh_skeleton ( mesh_b, nnodes, nelem, elshape, ndim )

  mesh_b%topology(1)%a=mesh%curves(5)%topology(:,:,1)
  mesh_b%coor=mesh%coor(mesh%curves(5)%nodes,:)

  call fill_mesh_parts ( mesh_b )

  if ( mesh_plot ) then

    call plot_mesh ( plot_options, mesh_b, 'mesh_b.fig' )

    call printinfo ( mesh_b, printlevel=4 )

  end if

! Velocity/pressure problem definition in external and internal fluid

  call create_input_probdef ( mesh, input_probdef, nvec=6, nphysq=3 )

! External fluid (group 1)
  input_probdef%vec_elementdof(1)%a(:,1) = 0         ! grad (set to zero)
  input_probdef%vec_elementdof(1)%a(:,2) = ndim      ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0         ! pressure
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1
  input_probdef%vec_elementdof(1)%a(:,4) = 1         ! scalar
  input_probdef%vec_elementdof(1)%a(:,5) = ndim*(ndim+1)/2  ! symmetric tensor
  input_probdef%vec_elementdof(1)%a(:,6) = ncomp_P2  ! projection

! Internal fluid (group 2)
  input_probdef%vec_elementdof(2)%a = input_probdef%vec_elementdof(1)%a
  input_probdef%vec_elementdof(2)%a(vertices,1) = ndim**2 ! grad

  input_probdef%physq = [physqgrad,physqvel,physqpress]
  input_probdef%probnr = 1

  if ( periodic ) then

!   dirichlet on lower and upper wall

    call define_essential ( mesh, input_probdef, curve1=1, physq=physqvel )
    call define_essential ( mesh, input_probdef, curve1=3, physq=physqvel )

!   velocity periodic in x-direction

    call define_constraint ( mesh, input_probdef, physq=physqvel, &
      curve1=2, curve2=4, discretization='collocation', excludecurves=[1,3], &
      num=c_per )

    if ( flowrate ) then

!     impose flowrate

      call define_constraint ( mesh, input_probdef, physq=physqvel, curve1=2, &
        nglobalc=1, num=c_flr )

    end if

  else

!   dirichlet on all flowcell boundaries

    call define_essential ( mesh, input_probdef, curve1=1, curve2=4, &
      physq=physqvel )

  end if

  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! Internal fluid - external fluid coupling through Lagrangian multipliers

  call define_constraint ( mesh, input_probdef, curve1=5, curve2=6, &
    discretization='collocation', physq=physqvel, num=c_cpl )

! Problem definition and create sysmatrix, vectors

  call problem_definition_create_sysmatrix_vectors

! monitor function

  if ( monitor_func ) allocate ( f_monitor(mesh_b%nnodes) )

! conformation problem

  call create_input_probdef ( mesh, input_probdef_c, nvec=4, nphysq=1 )

! External fluid (group 1)
  input_probdef_c%vec_elementdof(1)%a = 0 ! no degrees of freedom in external

! Internal fluid (group 2)
  input_probdef_c%vec_elementdof(2)%a(:,1) = 0               ! c
  input_probdef_c%vec_elementdof(2)%a(vertices,1) = 1
  input_probdef_c%vec_elementdof(2)%a(:,2) = 1               ! scalar
  input_probdef_c%vec_elementdof(2)%a(:,3) = ndim*(ndim+1)/2 ! symmetric tensor
  input_probdef_c%vec_elementdof(2)%a(:,4) = 0               ! projection
  input_probdef_c%vec_elementdof(2)%a(vertices,4) = ncomp_P1 ! projection

  input_probdef_c%physq = [1]
  input_probdef_c%probnr = 2

! Problem definition and create sysmatrix, vectors of conformation problem

  call problem_definition_create_sysmatrix_vectors_c

! conformation projection problem

  call create_input_probdef ( mesh, input_probdef_c_proj, nvec=1, nphysq=1 )

! External fluid (group 1)
  input_probdef_c_proj%vec_elementdof(1)%a = 0 ! no degrees of freedom

! Internal fluid (group 2)
  input_probdef_c_proj%vec_elementdof(2)%a(:,1) = 0               ! c
  input_probdef_c_proj%vec_elementdof(2)%a(vertices,1) = 1

  input_probdef_c_proj%physq = [1]
  input_probdef_c_proj%probnr = 3

  call problem_definition_create_sysmatrix_vectors_c_proj

! create the structure oldvectors_c
  call create_oldvectors ( oldvectors_c, nsysvec=2, nsysvec2=3, nprob=3, &
    nvec=1 )

! store solution vectors and problem structures

  oldvectors_c%s(1)%p => sol
  oldvectors_c%s(2)%p => solm1
  oldvectors_c%s2(1)%p => sol_c
  oldvectors_c%s2(2)%p => sol_cm1
  oldvectors_c%s2(3)%p => sol_c_proj
  oldvectors_c%p(1)%p => problem
  oldvectors_c%p(2)%p => problem_c
  oldvectors_c%p(3)%p => problem_c_proj
  oldvectors_c%v(1)%p => meshvel

! problem definition of the interface

  call create_input_probdef ( mesh_b, input_probdef_b, nvec=2 )

  input_probdef_b%elementdof(1)%a = 1
  input_probdef_b%vec_elementdof(1)%a(:,1) = ndim  ! 2D vector
  input_probdef_b%vec_elementdof(1)%a(:,2) = 1  ! error

  input_probdef_b%probnr = 4

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

  end if

! create system matrix

  call create_sysmatrix_structure ( sysmatrix_b, mesh_b, problem_b )

  call create_sysmatrix_data ( sysmatrix_b )

! create oldvectors

  call create ( oldvectors_b, nsysvec=1, nvec=4, nprob=1 )

  oldvectors_b%v(1)%p => vel_b  ! velocity u on the interface
  oldvectors_b%v(3)%p => xn     ! position at tn
  if ( timeint == 2 ) then
    oldvectors_b%v(4)%p => xnm1  ! position at tn-1
  end if

! allocate arrays for the ALE displacement problem

  allocate ( disp(mesh_b%nnodes,ndim) )

! define some arrays for the ALE mesh position at old times

  allocate ( meshcoor_n(mesh%nnodes,ndim), meshcoor_nm1(mesh%nnodes,ndim) )
  meshcoor_n = mesh%coor
  meshcoor_nm1 = mesh%coor

! Create aspect ratio arrays

  call create_aspect_ratio

! initialize cv, ... to xc because they are all written to restart

  cv = xc; cvn = xc; cvnm1 = xc

! open files and restart

  call open_files_and_restart

! time stepping

  post = post0
  euler_step = .true.

  do step = 1, numtimesteps

    write(*,'(a,i0)') 'step = ', step0+step

!   compute current aspect ratio
    call compute_element_aspect_ratio ( mesh, areav, aspect_ratio )

!   compute maximum normalized aspect ratio and area
!   (with respect to the initial values)
    norm_area = maxval( abs(log(areav/init_area)) )
    norm_aspect_ratio = maxval( abs(log(aspect_ratio/init_aspect_ratio)) )

!   remeshing criterion

    if ( ( norm_area >= area_threshold .or. &
           norm_aspect_ratio >= aspect_ratio_threshold ) .and. &
           step /= numtimesteps ) then

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
        node1=1, node2=mesh_b%nnodes, func=func_error_2D, funcnr=vfuncnr )

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

  end do

  if ( compute_error ) close(unit=10)
  if ( compute_volume ) close(unit=12)


! delete all data including all allocated memory

! flow problem

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, solm1 )
  call delete ( rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( ssvel_b )
  call delete ( meshvel )
  call delete ( ssmeshvelx )

! conformation problem

  call delete ( input_probdef_c )
  call delete ( problem_c )
  call delete ( sol_c, sol_cm1, rhs_c )
  call delete ( oldvectors_c )

! conformation projection

  call delete ( input_probdef_c_proj )
  call delete ( problem_c_proj )
  call delete ( sol_c_proj )

! monitor function

  if ( monitor_func ) then

    deallocate ( f_monitor )

!   curvature problem

    call delete ( problem_curv )
    call delete ( sysmatrix_curv )

  end if

! interface tangential motion problem

  call delete ( problem_pois )
  call delete ( sysmatrix_pois )
  call delete ( sol_pois )

! interface tracking problem

  call delete ( coefficients_b )
  call delete ( mesh_b )
  call delete ( input_probdef_b )
  call delete ( problem_b )
  call delete ( sol_b )
  call delete ( rhsd_b )
  call delete ( sysmatrix_b )
  call delete ( vel_b, xnp1, xn )
  if ( timeint == 2 ) call delete ( xnm1 )
  if ( compute_error ) call delete ( error )

! ALE mesh motion problem

  deallocate ( meshcoor_n, meshcoor_nm1 )
  deallocate ( disp )
  if ( numtimesteps > 1 ) call delete ( problem_update_mesh )

  call delete_aspect_ratio

contains

! generate and read mesh

  subroutine generate_and_read_mesh

!   generate external and internal mesh

    call write_gmsh_parameters ( ox, oy, lx, ly, xp, rp, dx_box, dx_part )

    call execute_command_line ( 'gmsh -2 -o mesh_conform.msh &
      &mesh_conform.geo > outputmesh_conform.out' )

!   read mesh generated by gmsh

    call read_mesh_gmsh ( mesh_conform, filename='mesh_conform.msh', ndim=2, &
      physgeom=.true. )

!   split mesh in external and internal mesh

    call mesh_convert ( mesh_conform, mesh_ext, only_groups=[1] )
    call mesh_convert ( mesh_conform, mesh_int, only_groups=[2] )

    call delete ( mesh_conform )

!   merge internal and external meshes

    call mesh_merge ( mesh_ext, mesh_int, mesh, nogroupmerge=.true. )

    call delete ( mesh_ext, mesh_int )

    if ( periodic ) then

!     local node numbering of curve 4 to match local node numbering of curve 2

      call add_to_mesh ( mesh, matchingcurve=[4,2], replace=4, &
        displacement=[-lx,0._dp] )

    end if

  end subroutine generate_and_read_mesh

  subroutine write_gmsh_parameters ( ox, oy, lx, ly, xp, rp, dx_box, dx_part )

    real(dp), intent(in) :: ox, oy, lx, ly, xp(:,:), rp(:), dx_box, dx_part

    integer :: i, nobj
    character(len=*), parameter :: fmti = '(1x,a,i0,a)'
    character(len=*), parameter :: fmtr = '(1x,a,f18.14,a)'
    character(len=*), parameter :: fmtir = '(1x,a,i0,a,f18.14,a)'

    nobj = size(rp)

    open ( unit=25, file='mesh_conform.geo' )

    write ( 25, fmtr ) 'ox = ', ox, ';'
    write ( 25, fmtr ) 'oy = ', oy, ';'
    write ( 25, fmtr ) 'lx = ', lx, ';'
    write ( 25, fmtr ) 'ly = ', ly, ';'
    write ( 25, fmtr ) 'dx_box = ', dx_box, ';'
    write ( 25, fmti ) 'nobj = ', nobj, ';'

    do i = 1, nobj
      write ( 25, fmtir ) 'xp[', i, '] = ', xp(i,1), ';'
      write ( 25, fmtir ) 'yp[', i, '] = ', xp(i,2), ';'
      write ( 25, fmtir ) 'rp[', i, '] = ', rp(i), ';'
    end do

    write ( 25, fmtr ) 'dx_part = ', dx_part, ';'

    write ( 25, '(/1x,a)' ) 'Include "mesh_options_2D.igo";'
    write ( 25, '(/1x,a)' ) 'Include "particles_in_a_box_2D_full.igo";'

    close ( 25 )

  end subroutine write_gmsh_parameters

! build and solve for velocity and pressure

  subroutine build_solve_velocity_pressure_conformation

    integer :: icomp

!   copy solution
    call copy ( sol, solm1 )

!   fill boundary conditions

    call fill_sysvector ( mesh, problem, sol, curve1=1, curve2=4, &
      physq=physqvel, vfunc=vfunc_2D, vfuncnr=vfuncnr )

    call fill_sysvector ( mesh, problem, sol, point=1, physq=physqpress, &
      value=0._dp )

!   Stokes velocity/pressure in external and internal fluid

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress], &
      elemsub=stokes_elem, mcoefficients=coefficients )

!   DEVSS-G (internal fluid)

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      physqrow=[physqgrad,physqvel], physqcol=[physqgrad,physqvel], &
      elemsub=devssg_elem, coefficients=coefficients(2), elgroup1=2, &
      addmatvec=.true. )

!   Setting to zero of off-diagonal blocks gradient-pressure (internal fluid)

    call build_system ( mesh, problem, sysmatrix, rhsd, buildvector=.false., &
      physqrow=[physqgrad], physqcol=[physqpress], &
      zeromatvec=.true., elgroup1=2, addmatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, buildvector=.false., &
      physqrow=[physqpress], physqcol=[physqgrad], &
      zeromatvec=.true., elgroup1=2, addmatvec=.true. )

    if ( periodic ) then

!     periodical condition on velocities in x-direction

      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        constraint1=c_per, elemsub=stokes_constr_node_conn, addmatvec=.true. )

      if ( flowrate ) then

!       impose flowrate

        call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
          constraint1=c_flr, elemsub=stokes_constr_flowr, addmatvec=.true., &
          coefficients=coefficients(1) )

      end if

    end if

!   Coupling of fluids

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=c_cpl, elemsub=stokes_constr_node_conn, addmatvec=.true. )

!   Surface tension

    call add_boundary_elements ( mesh, problem, rhsd, &
      elemsub=surface_tension_curve, curve=5, &
      coefficients=coefficients(1), physq=[physqvel] )

    call check ( sysmatrix )

!   exp(s) projection (for log conformation)

    if ( logc == 1 ) call solve_exps_projection

!   Building of implicit terms of CE with rhs in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, oldvectors=oldvectors_c, &
      physqrow=[physqvel], physqcol=[physqvel], &
      coefficients=coefficients(2), elgroup1=2, addmatvec=.true. )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

!   Building (assembling) of matrix and vector for conformation problem

    if ( coefficients(2)%i(22) == timeint1 ) then
      call build_system ( mesh, problem_c, sysmatrix_c, m2sysvector=rhs_c, &
        elemsub=ce_supg_elem, oldvectors=oldvectors_c, &
        coefficients=coefficients(2), elgroup1=2 )
    else
      call build_system ( mesh, problem_c, sysmatrix_c, m2sysvector=rhs_c, &
        elemsub=ce_supg_elem_implicit_2nd_order, oldvectors=oldvectors_c, &
        coefficients=coefficients(2), elgroup1=2 )
    end if

    call check ( sysmatrix_c )
    call copy ( sol_c, sol_cm1 )

!   solve conformation, keeping LU decomposition in the loop over components

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do icomp = 1, ncompc

      call solve_system_ma41 ( sysmatrix_c, rhs_c(icomp,1), sol_c(icomp,1), &
        lu_c, solver_options=solver_options_c )

    end do

    call delete ( lu_c ) ! remove LU decomposition for rebuild next time step

  end subroutine build_solve_velocity_pressure_conformation

! project c=exp(s) on discrete fem space

  subroutine solve_exps_projection

    type(solver_options_ma57_t) :: solver_options_ma57

    integer :: i, m

!   Building of system matrix and vector for projection problem

    call build_system ( mesh, problem_c_proj, sysmatrix_c_proj, &
      m2sysvector=rhs_c_proj, elemsub=exps_projection_elem, &
      oldvectors=oldvectors_c, coefficients=coefficients(2), &
      elgroup1=2 )

    call check ( sysmatrix_c_proj )

!   MA57 solver storage

    solver_options_ma57%integer_storage = 1.3
    solver_options_ma57%real_storage = 1.3

!   Solve projection problem, keep LU decomposition in the loop over components

    do m = 1, nmodes
      do i = 1, ncompc

        call add_effect_of_essential_to_rhs ( problem_c_proj, &
          sysmatrix_c_proj, sol_c_proj(i,m), rhs_c_proj(i,m) )

        call solve_system_ma57 ( sysmatrix_c_proj, rhs_c_proj(i,m), &
          sol_c_proj(i,m), lu_exps_proj, solver_options=solver_options_ma57 )

      end do
    end do

!   remove LU decomposition for rebuild next time step

    call delete ( lu_exps_proj )

  end subroutine solve_exps_projection

! build and solve for xnp1 on interface

  subroutine build_solve_interface

!   Velocity of the interface nodes (reshaped)
    real(dp), dimension(ndim,mesh_b%nnodes) :: drop_vel

    integer :: i
    real(dp) :: dx_min, dx_max

    if ( monitor_func ) then

!     Curvature problem on interface mesh

      call curvature_gd_ma57 ( mesh_b, problem_curv, coefficients_b, &
        sysmatrix_curv, curvature=f_monitor )

      f_monitor = dx_part / rc / f_monitor

      dx_min = fac_min * dx_part
      dx_max = fac_max * dx_part

      f_monitor = min ( dx_max, max(f_monitor,dx_min) )

!     Poisson problem on interface mesh for making nodes "diffuse" and
!     redistribute to highly curved regions

      call poisson_gd_ma57 ( mesh_b, problem_pois, coefficients_b, &
        sysmatrix_pois, sol_pois, f_monitor=f_monitor )

    else

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

!   build (assemble) matrix and vector from elements

    call build_system ( mesh_b, problem_b, sysmatrix_b, msysvector=rhsd_b, &
      elemsub=interface_tracking_elem_supg_impl, coefficients=coefficients_b, &
      oldvectors=oldvectors_b )

    do i = 1, ndim

      so%integer_storage=1.8_dp

      call solve_system_ma41 ( sysmatrix_b, rhsd_b(i), sol_b, lu=lu_b, &
        solver_options=so )

      call transfer_data ( mesh_b, problem1=problem_b, &
        sysvector1=sol_b, vector2=xnp1, degfd1=[1], degfd2=[i] )

    end do

    call delete(lu_b)

  end subroutine build_solve_interface

! solve flow and move the interface position (time discretized step)

  subroutine solve_flow_and_move_interface

!   copy mesh coordinates

    meshcoor_nm1 = meshcoor_n
    meshcoor_n = mesh%coor

    if ( timeint == 1 .or. &
               ( euler_step .and. timeint == 2 .and. restart /= 1 ) ) then

      if ( timeint == 2 ) euler_step = .false.

!     Euler

      call copy ( xnp1, xn )

      coefficients_b%i(11) = 1

!     update the nodes of the mesh

      call update_mesh_nodes_and_mesh_velocity ( order=1 )

    else if ( timeint == 2 ) then

!     Change of time integration scheme for internal fluid

      coefficients(2)%i(22) = timeint2

!     second order Gear with prediction

      call copy ( xn, xnm1 )
      call copy ( xnp1, xn )

!     2nd order prediction of interface coordinates
      mesh_b%coor =  &
                transpose( reshape( 2*xn%u - xnm1%u, [ndim,mesh_b%nnodes] ) )

      coefficients_b%i(11) = 2

!     update the nodes of the mesh

      call update_mesh_nodes_and_mesh_velocity ( order=2 )

    else

      write (*,'(/a,i0/)') 'Error: invalid timeint = ', timeint
      stop

    end if

!   solve velocity and pressure

    call build_solve_velocity_pressure_conformation

!   solve for interface velocity

    call build_solve_interface

!   Updating of interface mesh coordinates with new xnp1
    mesh_b%coor = transpose( reshape( xnp1%u, [ndim,mesh_b%nnodes] ) )

  end subroutine solve_flow_and_move_interface

! update the nodes of the mesh and compute the mesh velocity

  subroutine update_mesh_nodes_and_mesh_velocity ( order )

    integer, intent(in) :: order

!   displacement of the predicted interface with respect to the previous
!   mesh coordinates

    disp = mesh_b%coor - mesh%coor(mesh%curves(5)%nodes,:)

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

    if ( order == 1 ) then

!   Mesh velocity (BDF1)

      meshvel%u = reshape ( transpose ( &
        ( mesh%coor - meshcoor_n ) / deltat ), [ndim*mesh%nnodes] )

    else

!   Mesh velocity (BDF2)

      meshvel%u = reshape ( transpose ( &
        ( 1.5_dp*mesh%coor - 2*meshcoor_n + 0.5_dp*meshcoor_nm1 ) / deltat ), &
        [ndim*mesh%nnodes] )

    end if

    if ( move_mesh_x ) then

!   Additional ALE movement of the mesh in x-direction

      meshvel%u(ssmeshvelx%s) = meshvel%u(ssmeshvelx%s) + cvv(1)

    end if

  end subroutine update_mesh_nodes_and_mesh_velocity

! output various data

  subroutine output_data

    use output_fields3_m

    integer :: i

    write(*,'(a,i0,a,es10.3)') &
      'Output data at step ', step0+step, ', time =', time0+step*deltat

    write(filename,'(a,i4.4,a)') 'output', post, '.vtk'

    if ( compute_error ) then
!     error
      call write_scalar_vtk ( mesh_b, problem_b, filename=filename, &
        dataname='error', vector=error )
    end if

!   output coordinates
    write(filename,'(a,i4.4,a)') 'coor', post, '.out'
    open ( unit=11, file=filename, recl=300 )
    write(unit=11,fmt='(i0,2e16.8)') &
                           (i, mesh_b%coor(i,:), i=1,mesh_b%nnodes)
    close ( unit=11 )

!   output fields to vtk files

    call output_fields ( mesh, problem, problem_c, coefficients, sol, sol_c, &
      post )

  end subroutine output_data

! output restart data

  subroutine output_restart

    integer :: i

    write(*,'(a,i0,a,es10.3)') &
      'Output restart data at step ', step0+step, ', time =', time0+step*deltat

    call write_mesh ( mesh, filename='mesh.out' )

    open ( unit=11, form='unformatted', file='data.out' )

    write(11) step0 + step, post, time0 + step*deltat
    write(11) sol%u
    write(11) (sol_c(i,1)%u, i=1,ncompc)
    write(11) solm1%u
    write(11) (sol_cm1(i,1)%u, i=1,ncompc)
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

    integer :: endstep, endpost, i
    real(dp) :: endtime

    if ( restart >= 1 ) then

!     restart: read solution from file

      open ( unit=11, form='unformatted', file='data.out' )

      read(11) endstep, endpost, endtime
      read(11) sol%u
      read(11) (sol_c(i,1)%u, i=1,ncompc)
      read(11) solm1%u
      read(11) (sol_cm1(i,1)%u, i=1,ncompc)
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
        else
          open ( unit=10, file='error.out', recl=300, status='unknown' )
        end if

      end if

      if ( compute_volume ) then

!       open volume file

        if ( ctime ) then
          open ( unit=12, file='volume.out', recl=300, status='old', &
            position='append' )
        else
          open ( unit=12, file='volume.out', recl=300, status='unknown' )
        end if

      end if

    else

!     fresh start

      if ( compute_error ) then

!       open error file

        open ( unit=10, file='error.out', recl=300 )

      end if

      if ( compute_volume ) then

!       open volume file

        open ( unit=12, file='volume.out', recl=300 )

      end if

    end if

  end subroutine open_files_and_restart


! Problem definition of velocity/pressure system, sysmatrix and vectors

  subroutine problem_definition_create_sysmatrix_vectors

    call problem_definition ( input_probdef, mesh, problem )

!   Creation of mesh velocity vector

    call create ( problem, meshvel, physq=physqvel )

!   Creation of system vectors (solution and right-hand side)

    call create_sysvector ( problem, sol, solm1, rhsd )

!   Creation of system matrix for velocity/pressure problem

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

    call create_sysmatrix_data ( sysmatrix )

!   subscript for extracting velocities from mesh nodes on interface

    call create_subscript ( mesh, problem, ssvel_b, curves=[5], &
      physqarr=[physqvel] )

!   Creation of subscript for mesh velocity vector in x-direction

    call create_subscript_vector ( mesh, problem, ssmeshvelx, physq=physqvel, &
      degfd=1 )

  end subroutine problem_definition_create_sysmatrix_vectors

! Problem definition of conformation system, sysmatrix and vectors

  subroutine problem_definition_create_sysmatrix_vectors_c

    call problem_definition ( input_probdef_c, mesh, problem_c )

!   Creation of system vectors (solution and right-hand side) for conformation
!   and initialisation of vectors with zero stress
    call create ( problem_c, sol_c, sol_cm1, rhs_c )

    if ( logc == 0 ) then      ! standard
      sol_c(1,1)%u = 1         ! initial cxx
      sol_c(2,1)%u = 0         ! initial cxy
      sol_c(3,1)%u = 1         ! initial cyy
    else if ( logc == 1 ) then ! log scheme
      sol_c(1,1)%u = 0         ! initial cxx
      sol_c(2,1)%u = 0         ! initial cxy
      sol_c(3,1)%u = 0         ! initial cyy
    end if

!   Creation of system matrix for cnfromation problem
    call create_sysmatrix_structure ( sysmatrix_c, mesh, problem_c )
    call create_sysmatrix_data ( sysmatrix_c )

  end subroutine problem_definition_create_sysmatrix_vectors_c

  subroutine problem_definition_create_sysmatrix_vectors_c_proj

    call problem_definition ( input_probdef_c_proj, mesh, problem_c_proj )

!   Creation and initialisation of system vectors (solution and right-hand side)
!   for exp(s) projection
    call create ( problem_c_proj, sol_c_proj, rhs_c_proj )

    sol_c_proj(1,1)%u = 1         ! initial cxx
    sol_c_proj(2,1)%u = 0         ! initial cxy
    sol_c_proj(3,1)%u = 1         ! initial cyy

!   Creation of system matrix for projection problem
    call create_sysmatrix_structure ( sysmatrix_c_proj, mesh, problem_c_proj, &
      symmetric=.true. )
    call create_sysmatrix_data ( sysmatrix_c_proj )

  end subroutine problem_definition_create_sysmatrix_vectors_c_proj

! Delete problem of velocity/pressure system, sysmatrix and vectors

  subroutine delete_problem_sysmatrix_vectors

    call delete ( problem )
    call delete ( sol, rhsd )
    call delete ( sysmatrix )
    call delete ( ssvel_b )

    call delete ( meshvel )
    call delete ( solm1 )

    deallocate ( meshcoor_n, meshcoor_nm1 )

  end subroutine delete_problem_sysmatrix_vectors

  subroutine delete_problem_sysmatrix_vectors_c

    call delete ( problem_c )
    call delete ( sysmatrix_c )
    call delete ( sol_c, sol_cm1 )
    call delete ( rhs_c )

  end subroutine delete_problem_sysmatrix_vectors_c

  subroutine delete_problem_sysmatrix_vectors_c_proj

    call delete ( problem_c_proj )
    call delete ( sol_c_proj, rhs_c_proj )
    call delete ( sysmatrix_c_proj )

  end subroutine delete_problem_sysmatrix_vectors_c_proj

! Create and initialize the aspect ratio arrays needed for remesh criterion

  subroutine create_aspect_ratio

!   allocate aspect ratio arrays
    allocate ( init_area(mesh%nelem), areav(mesh%nelem) )
    allocate ( init_aspect_ratio(mesh%nelem), aspect_ratio(mesh%nelem) )

!   compute initial element aspect ratio
    call compute_element_aspect_ratio ( mesh, init_area, init_aspect_ratio )

  end subroutine create_aspect_ratio


! Delete the aspect ratio arrays

  subroutine delete_aspect_ratio

!   deallocate aspect ratio arrays
    deallocate ( init_area, areav )
    deallocate ( init_aspect_ratio, aspect_ratio )

  end subroutine delete_aspect_ratio

! compute the aspect ratio and area for each mesh element

  subroutine compute_element_aspect_ratio ( mesh, area, asp_ratio )

    type(mesh_t), intent(in) :: mesh
    real(dp), dimension(:), intent(out) :: asp_ratio, area

    integer :: elem, node(3), grp, totelem
    real(dp) :: la, lb, lc, sper
    real(dp) :: vert(3,2)

    node = [1,3,5]

    totelem = 0

    do grp = 1, mesh%nelgrp
      do elem = 1,mesh%grpnumel(grp)

        totelem = totelem + 1

!       get coordinates of triangle vertices (local node numbers: 1, 3, 5)
        vert = mesh%coor(mesh%topology(grp)%a(node,elem),:)

!       compute side lengths
        la = sqrt( (vert(1,1) - vert(2,1))**2 + (vert(1,2) - vert(2,2))**2 )
        lb = sqrt( (vert(2,1) - vert(3,1))**2 + (vert(2,2) - vert(3,2))**2 )
        lc = sqrt( (vert(3,1) - vert(1,1))**2 + (vert(3,2) - vert(1,2))**2 )

!       compute triangle semiperimeter and area
        sper = (la + lb + lc)/2.0_dp
        area(totelem) = sqrt( sper*(sper - la)*(sper - lb)*(sper - lc) )

!       compute aspect ratio
        asp_ratio(totelem) = max(la,lb,lc)**2/area(totelem)

      end do
    end do

  end subroutine compute_element_aspect_ratio

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
    write ( 25, '(/1x,a)' ) 'Include "mesh_options_2D.igo";'
    write ( 25, '(/1x,a)' ) 'Merge "mesh_geometries.msh";'
    write ( 25, '(/1x,a)' ) 'Include "particles_in_a_box_2D_full_2.igo";'

    close ( 25 )

!   save coordinates of the boundary mesh

    allocate ( coor_boun(mesh%curves(5)%nnodes,mesh%ndim) )
    coor_boun = mesh%coor(mesh%curves(5)%nodes,:)

    call execute_command_line ( 'gmsh -2 -bin -o mesh_conform.msh &
      &mesh_conform.geo > outputmesh_conform.out' )

!   create vectors to store the mesh coordinates, these coordinates will
!   be used in the projection problem to find the new mesh at t_n and t_nm1

    call create ( problem, meshcoor_n_temp, vec=2 )
    call create ( problem, meshcoor_nm1_temp, vec=2 )

    meshcoor_n_temp%u = reshape ( transpose(meshcoor_n), [size(meshcoor_n)] )
    meshcoor_nm1_temp%u = reshape ( transpose(meshcoor_nm1), &
                                    [size(meshcoor_nm1)] )

!   copy current definitions to temporary definitions used for the projection
!   problem

    call copy ( mesh, mesh_temp )
    call fill_mesh_parts ( mesh_temp )
    call find_bounds_blocks ( mesh_temp )
    call copy ( sol, sol_temp )
    call copy ( sol_c, sol_c_temp )
    call copy ( sol_cm1, sol_cm1_temp )
    call problem_definition ( input_probdef, mesh, problem_temp )
    call problem_definition ( input_probdef_c, mesh, problem_c_temp )

!   read mesh generated by gmsh

    call read_mesh_gmsh ( mesh_conform, filename='mesh_conform.msh', ndim=2, &
      physgeom=.true. )

!   split mesh in external and internal mesh

    call mesh_convert ( mesh_conform, mesh_ext, only_groups=[1] )
    call mesh_convert ( mesh_conform, mesh_int, only_groups=[2] )

    call delete ( mesh_conform, mesh )

!   merge internal and external meshes

    call mesh_merge ( mesh_ext, mesh_int, mesh, nogroupmerge=.true. )

    call delete ( mesh_ext, mesh_int )

    if ( periodic ) then

!     local node numbering of curve 4 to match local node numbering of curve 2

      call add_to_mesh ( mesh, matchingcurve=[4,2], replace=4, &
        displacement=[-lx,0._dp] )

    end if

    call fill_mesh_parts ( mesh )
    call find_bounds_blocks ( mesh )

!   delete various structures, arrays

    call delete_problem_sysmatrix_vectors
    call delete_problem_sysmatrix_vectors_c
    call delete_problem_sysmatrix_vectors_c_proj
    call delete_aspect_ratio
    call delete ( problem_update_mesh )

!   recreate the structures, arrays based on the new mesh

    call problem_definition_create_sysmatrix_vectors
    call problem_definition_create_sysmatrix_vectors_c
    call problem_definition_create_sysmatrix_vectors_c_proj
    call create_aspect_ratio

!   Project the old solutions onto the new mesh and fill in the values

    call project_old_solutions_P2
    call project_old_solutions_P1

!   move boundary nodes to previous position (gmsh does not preserve quadratic
!   element shape).

    mesh%coor(mesh%curves(5)%nodes,:) = coor_boun
    mesh%coor(mesh%curves(6)%nodes,:) = coor_boun

    deallocate ( coor_boun )

!   Delete the temporary definitions

    call delete ( meshcoor_n_temp, meshcoor_nm1_temp )
    call delete ( mesh_temp )
    call delete ( sol_temp )
    call delete ( sol_c_temp, sol_cm1_temp )
    call delete ( problem_temp )
    call delete ( problem_c_temp )

  end subroutine remeshing

  subroutine project_old_solutions_P2

    type(input_probdef_t) :: input_probdef_proj
    type(problem_t) :: problem_proj
    type(sysmatrix_t) :: sysmatrix_proj
    type(sysvector_t) :: rhsd_proj(ncomp_P2)
    type(sysvector_t) :: sol_proj(ncomp_P2)
    type(oldvectors_t) :: oldvectors_proj
    type(coefficients_t) :: coefficients_proj
    type(vector_t), target :: vec_proj
    type(subscript_t) :: velx, vely
    type(lu_ma57_t) :: lu_proj

    integer :: i

!   create vector for projection
    call create ( problem_temp, vec_proj, vec=6 )

!   transfer the coordinates at n
    call transfer_data ( mesh_temp, problem_temp, vector1=meshcoor_n_temp, &
      vector2=vec_proj, degfd2=[1,2] )

!   transfer the coordinates at n-1
    call transfer_data ( mesh_temp, problem_temp, vector1=meshcoor_nm1_temp, &
      vector2=vec_proj, degfd2=[3,4] )

!   transfer the velocities at n
    call transfer_data ( mesh_temp, problem_temp, sysvector1=sol_temp, &
      vector2=vec_proj, physq1=[physqvel], degfd1=[1,2], degfd2=[5,6] )

!   problem definition for projection
    call create_input_probdef ( mesh, input_probdef_proj )

    input_probdef_proj%elementdof(1)%a = 1
    input_probdef_proj%elementdof(2)%a = 1
    call problem_definition ( input_probdef_proj, mesh, problem_proj )

!   create system vectors (solution and right-hand side)
    call create ( problem_proj, sol_proj )
    call create ( problem_proj, rhsd_proj )

!   create system matrix
    call create_sysmatrix_structure ( sysmatrix_proj, mesh, problem_proj, &
      symmetric=.true. )
    call create_sysmatrix_data ( sysmatrix_proj )

!   fill coefficients
    call create_coefficients ( coefficients_proj, ncoefi=100, ncoefr=50 )
    coefficients_proj%i = 0
    coefficients_proj%i(1:4) = [ uintpl, ncomp_P2, uintpl, uintpl ]
    coefficients_proj%i(10) = gauss_proj
    coefficients_proj%i(40) = inttype_proj
    coefficients_proj%r = 0

!   oldvectors
    call create_oldvectors ( oldvectors_proj, nvec=1, nprob=1, nmesh=1 )
    oldvectors_proj%v(1)%p => vec_proj
    oldvectors_proj%p(1)%p => problem_temp
    oldvectors_proj%m(1)%p => mesh_temp

!   build (assemble) matrix and vector from elements
    call build_system ( mesh, problem_proj, sysmatrix_proj, &
      msysvector=rhsd_proj, elemsub=projection_elem, &
      coefficients=coefficients_proj, oldvectors=oldvectors_proj )

!   solve the projection problem
    do i = 1, ncomp_P2
      call add_effect_of_essential_to_rhs ( problem_proj, sysmatrix_proj, &
        sol_proj(i), rhsd_proj(i) )
      call solve_system_ma57 ( sysmatrix_proj, rhsd_proj(i), sol_proj(i), &
        lu=lu_proj )
    end do

!   create meshes at t_n and t_nm1 and update the old coordinates
    allocate ( meshcoor_n(mesh%nnodes,ndim), meshcoor_nm1(mesh%nnodes,ndim) )
    meshcoor_n(:,1) = sol_proj(1)%u
    meshcoor_n(:,2) = sol_proj(2)%u
    meshcoor_nm1(:,1) = sol_proj(3)%u
    meshcoor_nm1(:,2) = sol_proj(4)%u
!   update sol at t_n
    call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
    call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
    sol%u(velx%s) = sol_proj(5)%u
    sol%u(vely%s) = sol_proj(6)%u
    call delete ( velx, vely )

!   refill meshvel using the projected coordinates
    if ( timeint == 1 ) then
!     mesh velocity (BDF1)
      meshvel%u = reshape ( transpose ( &
        ( mesh%coor - meshcoor_n ) / deltat ), [ndim*mesh%nnodes] )
    else
!     mesh velocity (BDF2)
      meshvel%u = reshape ( transpose ( &
        ( 1.5_dp*mesh%coor - 2*meshcoor_n + 0.5_dp*meshcoor_nm1 ) / deltat ), &
        [ndim*mesh%nnodes] )
    end if

!   delete definitions used for the projection problem
    call delete ( input_probdef_proj )
    call delete ( problem_proj )
    call delete ( sysmatrix_proj )
    call delete ( rhsd_proj )
    call delete ( sol_proj )
    call delete ( oldvectors_proj )
    call delete ( coefficients_proj )
    call delete ( vec_proj )
    call delete ( lu_proj )

  end subroutine project_old_solutions_P2

  subroutine project_old_solutions_P1

    type(input_probdef_t) :: input_probdef_proj
    type(problem_t) :: problem_proj
    type(sysmatrix_t) :: sysmatrix_proj
    type(sysvector_t) :: rhsd_proj(ncomp_P1)
    type(sysvector_t) :: sol_proj(ncomp_P1)
    type(oldvectors_t) :: oldvectors_proj
    type(coefficients_t) :: coefficients_proj
    type(vector_t), target :: vec_proj
    type(lu_ma57_t) :: lu_proj
    type(subscript_t) :: cc, cc_proj

    integer :: i

!   create vector for projection
    call create ( problem_c_temp, vec_proj, vec=4 )

!   transfer the conformation tensor at n
    call transfer_data ( mesh_temp, problem_c_temp, &
      sysvector1=sol_c_temp(1,1), vector2=vec_proj, physq1=[1], degfd2=[1] )
    call transfer_data ( mesh_temp, problem_c_temp, &
      sysvector1=sol_c_temp(2,1), vector2=vec_proj, physq1=[1], degfd2=[2] )
    call transfer_data ( mesh_temp, problem_c_temp, &
      sysvector1=sol_c_temp(3,1), vector2=vec_proj, physq1=[1], degfd2=[3] )

!   transfer the conformation tensor at n-1
    call transfer_data ( mesh_temp, problem_c_temp, &
      sysvector1=sol_cm1_temp(1,1), vector2=vec_proj, physq1=[1], degfd2=[4] )
    call transfer_data ( mesh_temp, problem_c_temp, &
      sysvector1=sol_cm1_temp(2,1), vector2=vec_proj, physq1=[1], degfd2=[5] )
    call transfer_data ( mesh_temp, problem_c_temp, &
      sysvector1=sol_cm1_temp(3,1), vector2=vec_proj, physq1=[1], degfd2=[6] )

!   problem definition for projection
    call create_input_probdef ( mesh, input_probdef_proj )

    input_probdef_proj%elementdof(1)%a = 0
    input_probdef_proj%elementdof(2)%a = [1,0,1,0,1,0]
    call problem_definition ( input_probdef_proj, mesh, problem_proj )

!   create system vectors (solution and right-hand side)
    call create ( problem_proj, sol_proj )
    call create ( problem_proj, rhsd_proj )

!   create system matrix
    call create_sysmatrix_structure ( sysmatrix_proj, mesh, problem_proj, &
      symmetric=.true. )
    call create_sysmatrix_data ( sysmatrix_proj )

!   fill coefficients
    call create_coefficients ( coefficients_proj, ncoefi=100, ncoefr=50 )
    coefficients_proj%i = 0
    coefficients_proj%i(1:4) = [ gintpl, ncomp_P1, uintpl, gintpl ]
    coefficients_proj%i(10) = gauss_proj
    coefficients_proj%i(40) = inttype_proj
    coefficients_proj%r = 0

!   oldvectors
    call create_oldvectors ( oldvectors_proj, nvec=1, nprob=1, nmesh=1 )
    oldvectors_proj%v(1)%p => vec_proj
    oldvectors_proj%p(1)%p => problem_c_temp
    oldvectors_proj%m(1)%p => mesh_temp

!   build (assemble) matrix and vector from elements
    call build_system ( mesh, problem_proj, sysmatrix_proj, &
      msysvector=rhsd_proj, elemsub=projection_elem, &
      coefficients=coefficients_proj, oldvectors=oldvectors_proj, elgroup1=2 )

!   solve the projection problem
    do i = 1, ncomp_P1
      call add_effect_of_essential_to_rhs ( problem_proj, sysmatrix_proj, &
        sol_proj(i), rhsd_proj(i) )
      call solve_system_ma57 ( sysmatrix_proj, rhsd_proj(i), sol_proj(i), &
        lu=lu_proj )
    end do

!   create meshes at t_n and t_nm1 and update the old coordinates
    call create_subscript ( mesh, problem_proj, cc_proj, groups=[2] )
    call create_subscript ( mesh, problem_c, cc, physqarr=[1] )
    sol_c(1,1)%u(cc%s) = sol_proj(1)%u(cc_proj%s)
    sol_c(2,1)%u(cc%s) = sol_proj(2)%u(cc_proj%s)
    sol_c(3,1)%u(cc%s) = sol_proj(3)%u(cc_proj%s)
    sol_cm1(1,1)%u(cc%s) = sol_proj(4)%u(cc_proj%s)
    sol_cm1(2,1)%u(cc%s) = sol_proj(5)%u(cc_proj%s)
    sol_cm1(3,1)%u(cc%s) = sol_proj(6)%u(cc_proj%s)
    call delete ( cc, cc_proj )

!   delete definitions used for the projection problem
    call delete ( input_probdef_proj )
    call delete ( problem_proj )
    call delete ( sysmatrix_proj )
    call delete ( rhsd_proj )
    call delete ( sol_proj )
    call delete ( oldvectors_proj )
    call delete ( coefficients_proj )
    call delete ( vec_proj )
    call delete ( lu_proj )

  end subroutine project_old_solutions_P1

end program drop11
