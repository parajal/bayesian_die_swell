! Multiple drops consisting of a Newtonian fluid suspended in a Newtonian fluid.
! Stokes flow (no inertia).
! Flow imposed either:
!   1) on the whole boundary using a function, or
!   2) on the upper and lower boundary using a function and
!      periodic boundary conditions in x- and z-direction (periodic=.true.)
!   3) periodic boundary conditions in x-direction and on the other boundaries
!      imposed by a function (periodic_channel=.true.)
! ALE with interface tracking and optional displacement of the mesh in
! x-direction to keep the center of volume at the same position in the mesh.
! Remeshing is performed when the maximum normalized volume or aspect ratio
! of the mesh elements is higher than a threshold.
! Remeshing is based on position of nodes of the interface meshes, but the
! interface mesh itself are not modified.
! Optionally move the nodes of the interfaces mesh to highly curved regions.
! Optionally impose flowrate if periodic=.true.
! Optionally compute the volume and center of volume.
! Three-dimensional model.

program drop8

  use tfem_m
  use math_defs_m
  use hsl_ma57_m
  use hsl_ma41_m
  use stokes_elements_m
  use interface_tracking_elements_m
  use poisson_gd_ma57_m
  use curvature_gd_ma57_m
  use update_mesh_nodes8_m
  use functions_m
  use io_utils_m
  use figplot_m

  implicit none


! constants

  integer, parameter ::  &
    uintpl = 6,          & ! P2 velocity and interface position interpolation
    pintpl = 2,          & ! P1 pressure interpolation
    physqvel = 1,        & ! physical quantity nr of velocity
    physqpress = 2,      & ! physical quantity nr of pressure
    nobj=2,              & ! number of drops
    ndim = 3               ! dimension of space

! variables

  logical ::     &
    periodic = .false.,     & ! Use periodical boundaries in x and z-direction
                              ! Only makes sense for vfuncnr=0 (zero velocity)
                              ! or vfuncnr=2 (shear flow).
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
    compute_volume = .false.  ! Compute volume and center of volume.

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
    vtkevery = 5 ,       & ! vtk file every vtkevery steps. -1: means none
    restart_every = 100, & ! write restart file every restart_every steps.
                           ! -1: means none
    timeint=1,           & ! time integration scheme
    numtimesteps=400       ! number of time steps

  real(dp) ::  &
    ox = -0.5_dp,         & ! x-coordinate left lower corner
    oy = -0.5_dp,         & ! y-coordinate left lower corner
    oz = -0.5_dp,         & ! z-coordinate left lower corner
    lx = 1.0_dp,          & ! length of the domain
    ly = 1.0_dp,          & ! width of the domain
    lz = 1.0_dp,          & ! height of the domain
    rc = 0.1_dp,          & ! radius of the initial sphere
    eta_e = 1.0_dp,       & ! external fluid viscosity
    eta_i = 1.0_dp,       & ! internal fluid viscosity
    gammac = 0.5_dp,      & ! surface tension coefficient
    factor = 100._dp,     & ! factor for the interface grid deformation
    rs_up  = 1.4_dp,      & ! real_storage velocity-pressure LU HSL
    is_up  = 1.5_dp,      & ! integer_storage velocity-pressure LU HSL
    U_avg  = 1.0_dp,      & ! imposed average velocity for flowrate=.true.
    time0  = 0.0_dp,      & ! initial time
    beta = 0.5_dp,        & ! factor beta in SUPG
    deltat = 7.e-3_dp,    & ! time step
    fac_min = 0.1_dp,     & ! factor setting dx_min for the monitor function
    fac_max = 10.0_dp,    & ! factor setting dx_max for the monitor function
    volume_threshold = 1.39_dp,      & ! remeshing volume threshold
    aspect_ratio_threshold = 1.39_dp   ! remeshing aspect ratio threshold

! flow problem

  type(mesh_t) :: mesh_ext, mesh_int, mesh, mesh_conform
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(coefficients_t) :: coefficients(2)
  type(plot_options_t) :: plot_options
  type(solver_options_ma57_t) :: solver_options_ma57
  type(subscript_t) :: ssvel_b(nobj)

! curvature problem

  type(problem_t) :: problem_curv(nobj)
  type(sysmatrix_t) :: sysmatrix_curv(nobj)

! interface tangential motion problem

  type(problem_t), target :: problem_pois(nobj)
  type(sysmatrix_t) :: sysmatrix_pois(nobj)
  type(sysvector_t), target :: sol_pois(nobj)

! interface tracking problem

  type(mesh_t) :: mesh_b(nobj)
  type(input_probdef_t) :: input_probdef_b(nobj)
  type(problem_t), target :: problem_b(nobj)
  type(sysmatrix_t) :: sysmatrix_b(nobj)
  type(sysvector_t), target :: sol_b(nobj)
  type(sysvector_t) :: rhsd_b(nobj,ndim)
   type(vector_t) :: xnp1(nobj)
! Interface mesh coordinates at previous times. NOTE: the coordinates are
! defined with respect to the initial coordinate system.
  type(vector_t), target :: xn(nobj), xnm1(nobj)
  type(vector_t), target :: vel_b(nobj)
  type(oldvectors_t) :: oldvectors_b(nobj)
  type(coefficients_t) :: coefficients_b
  type(lu_ma41_t) :: lu_b
  type(solver_options_ma41_t) :: so


! ALE mesh motion problem

  type(problem_t), target ::  problem_update_mesh
  type(real_array_2d_t), dimension(nobj) :: disp

! some more misc variables

  logical :: w1, euler_step
  character(len=30) :: filename

  integer, parameter :: reclv = (1+(ndim+1)*nobj)*26
  integer :: nnodes, nelem, elshape, step, post, i, c_cpl(nobj), c_perx, &
    c_perz, c_flr, j
  integer :: vertices(4)=[1,3,5,10]
  real(dp) :: xp(nobj,ndim), rp(nobj), dx_box, dx_part, resultsum(1+ndim), &
              res(nobj,1+ndim)

! deformation variables for remeshing
  real(dp), dimension(:), allocatable :: init_volume, volumev
  real(dp), dimension(:), allocatable :: init_aspect_ratio, aspect_ratio
  real(dp) :: norm_volume, norm_aspect_ratio
  type(real_array_1d_t), dimension(nobj) :: f_monitor

! Center of volume position of the drops at current and previous times.
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
  real(dp), dimension(nobj,ndim) :: cv, cvn, cvnm1

! center of volume velocity
  real(dp), dimension(nobj,ndim) :: cvv


! namelist for input of variables; read from standard input

  namelist /comppar/ periodic, periodic_channel, flowrate, move_mesh_x, &
    monitor_func, mesh_plot, compute_volume, n, nc, gauss, gaussb, &
    choice_um, vfuncnr, vtkevery, ctime, restart, restart_every, &
    timeint, numtimesteps, ox, oy, oz, lx, ly, lz, rc, eta_e, eta_i, &
    gammac, factor, U_avg, rs_up, is_up, deltat, beta, fac_min, fac_max, &
    volume_threshold, aspect_ratio_threshold

  read ( unit=*, nml=comppar )
  read ( unit=*, fmt=* )
  read ( unit=*, fmt=* ) ( xp(i,:), i=1,nobj )


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
  coefficients(1)%r(19) = gammac

! Internal fluid (element group 2)

  coefficients(2)%i = coefficients(1)%i

  coefficients(2)%r = 0
  coefficients(2)%r(1) = eta_i

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
  coefficients_b%r(8) = beta
  coefficients_b%r(9) = deltat

! create mesh

  dx_box = lx/n
  dx_part = 2*pi*rc/nc

  if ( restart == 0 ) then

!   generate fluid mesh using gmsh

    rp = rc

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

  do j = 1, nobj

    nnodes = mesh%surfaces(6+j)%nnodes
    nelem = mesh%surfaces(6+j)%nelem
    elshape = mesh%surfaces(6+j)%element%elshape

    call mesh_skeleton ( mesh_b(j), nnodes, nelem, elshape, ndim )

    mesh_b(j)%topology(1)%a=mesh%surfaces(6+j)%topology(:,:,1)
    mesh_b(j)%coor=mesh%coor(mesh%surfaces(6+j)%nodes,:)

    call fill_mesh_parts ( mesh_b(j) )

  end do

! Velocity/pressure problem definition in external and internal fluid

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = ndim      ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0         ! pressure
  input_probdef%vec_elementdof(1)%a(vertices,2) = 1
  input_probdef%vec_elementdof(1)%a(:,3) = 1         ! scalar
  input_probdef%vec_elementdof(1)%a(:,4) = ndim*(ndim+1)/2  ! symmetric tensor

  input_probdef%vec_elementdof(2)%a = input_probdef%vec_elementdof(1)%a

  input_probdef%physq = [physqvel,physqpress]
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

  do j = 1, nobj

    call define_constraint ( mesh, input_probdef, surface1=6+j, &
      surface2=6+nobj+j, discretization='collocation', physq=physqvel, &
      num=c_cpl(j) )

  end do

! Problem definition and create sysmatrix, vectors

  call problem_definition_create_sysmatrix_vectors

! monitor function

  if ( monitor_func ) then

    do j = 1, nobj
      allocate ( f_monitor(j)%a(mesh_b(j)%nnodes) )
    end do

  end if

! problem definition of the interface

  call problem_definition_interface_create_sysmatrix_vectors

! allocate arrays for the ALE displacement problem

  do j = 1, nobj
    allocate ( disp(j)%a(mesh_b(j)%nnodes,ndim) )
  end do

! Create aspect ratio arrays

  call create_aspect_ratio

! inititalize cv, ... to xc because they are all written to restart

  cv = xp; cvn = xp; cvnm1 = xp


! open files and restart

  call open_files_and_restart

! time stepping

  post = post0
  euler_step = .true.

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
           step /= numtimesteps ) then

      print *,'Doing remeshing...'

      call remeshing

    end if

!   move the interface position and solve flow

    call solve_flow_and_move_interface

!   volume and center of volume of the interface

    if ( compute_volume ) then

!     compute volume and center of volume

      do j = 1, nobj

        call integrate ( mesh_b(j), problem_b(j), resultsum, &
          elemsub=center_of_volume, coefficients=coefficients_b )

        res(j,:) = [ resultsum(1), resultsum(2:)/resultsum(1) ]

      end do

      time = time0 + step * deltat

      write(unit=12,fmt=*) time, ( res(j,:), j = 1, nobj )

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

  if ( compute_volume ) close(unit=12)


! delete all data including all allocated memory

! flow problem

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol )
  call delete ( rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )
  do j = 1, nobj
    call delete ( ssvel_b(j) )
  end do

! monitor function

  if ( monitor_func ) then

    do j = 1, nobj

      deallocate ( f_monitor(j)%a )

!     curvature problem

      call delete ( problem_curv(j) )
      call delete ( sysmatrix_curv(j) )

    end do

  end if

! interface tangential motion problem

  do j = 1, nobj
    call delete ( problem_pois(j) )
    call delete ( sysmatrix_pois(j) )
    call delete ( sol_pois(j) )
  end do

! interface tracking problem

  call delete_problem_interface_create_sysmatrix_vectors

! ALE mesh motion problem

  do j = 1, nobj
    deallocate ( disp(j)%a )
  end do

  if ( numtimesteps > 1 ) call delete ( problem_update_mesh )


contains


! generate and read mesh

  subroutine generate_and_read_mesh

!   generate external and internal mesh

    call write_gmsh_parameters ( ox, oy, oz, lx, ly, lz, xp, rp, &
      dx_box, dx_part )

    call execute_command_line ( 'gmsh -3 -o mesh_conform.msh &
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

    else if ( periodic_channel ) then

!     local node numbering surface 4 to match local node numbering of surface 2

      call add_to_mesh ( mesh, matchingsurface=[4,2], replace=4, &
        displacement=[-lx,0._dp,0._dp] )

    end if

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

    integer :: j

!   fill boundary conditions

    call fill_sysvector ( mesh, problem, sol, surface1=1, surface2=6, &
      physq=physqvel, vfunc=vfunc_3D, vfuncnr=vfuncnr )

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
      constraint1=c_cpl(1), constraint2=c_cpl(nobj), &
      elemsub=stokes_constr_node_conn, addmatvec=.true. )

!   Surface tension

    do j = 1, nobj
      call add_boundary_elements ( mesh, problem, rhsd, &
        elemsub=surface_tension_surface, surface=6+j, &
        coefficients=coefficients(1), physq=[physqvel] )
    end do

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    solver_options_ma57%real_storage=rs_up
    solver_options_ma57%integer_storage=is_up

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_ma57  )

  end subroutine build_solve_velocity_pressure


! build and solve for xnp1 on interface

  subroutine build_solve_interface

    integer :: i, j

!   Velocity of the interface nodes (reshaped)
    real(dp), allocatable, dimension(:,:) :: drop_vel

    real(dp) :: dx_min, dx_max

    do j = 1, nobj

      if ( monitor_func ) then

!       Curvature problem on interface mesh

        call curvature_gd_ma57 ( mesh_b(j), problem_curv(j), coefficients_b, &
          sysmatrix_curv(j), curvature=f_monitor(j)%a, square=.true. )

        f_monitor(j)%a = 4 * dx_part ** 2 / rc ** 2 / f_monitor(j)%a

        dx_min = fac_min * dx_part
        dx_max = fac_max * dx_part

        f_monitor(j)%a = min ( dx_max**2, max(f_monitor(j)%a,dx_min**2) )

!       Poisson problem on interface mesh for making nodes "diffuse" and
!       redistribute to highly curved regions

        call poisson_gd_ma57 ( mesh_b(j), problem_pois(j), coefficients_b, &
          sysmatrix_pois(j), sol_pois(j), f_monitor=f_monitor(j)%a )

      else

!       Poisson problem on interface mesh for making nodes "diffuse" and
!       redistribute uniformly

        call poisson_gd_ma57 ( mesh_b(j), problem_pois(j), coefficients_b, &
          sysmatrix_pois(j), sol_pois(j) )

      end if

      oldvectors_b(j)%s(1)%p => sol_pois(j)
      oldvectors_b(j)%p(1)%p => problem_pois(j)

      vel_b(j)%u=sol%u(ssvel_b(j)%s) ! extract velocities from velocity solution

      if ( choice_um == 1 ) then

!       average nodal velocities

        nnodes = mesh_b(j)%nnodes

        allocate(drop_vel(ndim,nnodes))

        drop_vel = reshape(vel_b(j)%u,[ndim,nnodes]) ! reshape
        do i = 1, ndim
          coefficients_b%r(3+i) = sum(drop_vel(i,:))/nnodes
        end do

        deallocate(drop_vel)

      else if ( choice_um == 2 ) then

!       velocity of the center of volume

        coefficients_b%r(4:3+ndim) = cvv(j,:)

      end if

!     build (assemble) matrix and vector from elements

      call build_system ( mesh_b(j), problem_b(j), sysmatrix_b(j), &
        msysvector=rhsd_b(j,:), elemsub=interface_tracking_elem_supg_impl, &
        coefficients=coefficients_b, oldvectors=oldvectors_b(j) )

      do i = 1, ndim

        so%integer_storage=1.8_dp

        call solve_system_ma41 ( sysmatrix_b(j), rhsd_b(j,i), sol_b(j), &
          lu=lu_b, solver_options=so )

        call transfer_data ( mesh_b(j), problem1=problem_b(j), &
          sysvector1=sol_b(j), vector2=xnp1(j), degfd1=[1], degfd2=[i] )

      end do

      call delete(lu_b)

    end do

  end subroutine build_solve_interface


! solve flow and move the interface position (time discretized step)

  subroutine solve_flow_and_move_interface

    integer :: j

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
      do j = 1, nobj
        mesh_b(j)%coor = &
          transpose( reshape( 2*xn(j)%u - xnm1(j)%u, [ndim,mesh_b(j)%nnodes] ) )
      end do

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
    do j = 1, nobj
      mesh_b(j)%coor = &
             transpose( reshape( xnp1(j)%u, [ndim,mesh_b(j)%nnodes] ) )
    end do

  end subroutine solve_flow_and_move_interface


! update the nodes of the mesh and compute the mesh velocity

  subroutine update_mesh_nodes_from_displacement ( order )

    integer, intent(in) :: order

    integer :: j

!   displacement of the predicted interface with respect to the previous
!   mesh coordinates

    do j = 1, nobj
      disp(j)%a = mesh_b(j)%coor - mesh%coor(mesh%surfaces(6+j)%nodes,:)
    end do

    if ( move_mesh_x .or. choice_um == 2 ) then

!     update center of volume

      cvnm1 = cvn
      cvn = cv

      do j = 1, nobj

        call integrate ( mesh_b(j), problem_b(j), resultsum, &
          elemsub=center_of_volume, coefficients=coefficients_b )

        cv(j,:) = resultsum(2:)/resultsum(1)

      end do

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

      do j = 1, nobj
        disp(j)%a(:,1) = disp(j)%a(:,1) - sum ( cv(:,1) - xp(:,1) ) / nobj
      end do

    end if

!   update mesh nodes

    call update_mesh_nodes ( mesh, problem_update_mesh, disp=disp )

  end subroutine update_mesh_nodes_from_displacement


! output various data

  subroutine output_data

    use output_fields1_m

    integer :: i, j

    write(*,'(a,i0,a,es10.3)') &
      'Output data at step ', step0+step, ', time =', time0+step*deltat

    write(filename,'(a,i4.4,a)') 'output', post, '.vtk'

!   output coordinates
    write(filename,'(a,i4.4,a)') 'coor', post, '.out'
    open ( unit=11, file=filename, recl=300 )
    do j = 1, nobj
      write(unit=11,fmt='(i0,3e16.8)') &
                           (i, mesh_b(j)%coor(i,:), i=1,mesh_b(j)%nnodes)
      write(unit=11,fmt=*)
    end do
    close ( unit=11 )

!   output fields to vtk files

    call output_fields ( mesh, problem, coefficients, sol, post )

  end subroutine output_data


! output restart data

  subroutine output_restart

    integer :: j

    write(*,'(a,i0,a,es10.3)') &
      'Output restart data at step ', step0+step, ', time =', time0+step*deltat

    call write_mesh ( mesh, filename='mesh.out' )

    open ( unit=11, form='unformatted', file='data.out' )

    write(11) step0 + step, post, time0 + step*deltat
    write(11) (mesh_b(j)%coor, j=1,nobj)
    write(11) (xn(j)%u, xnp1(j)%u, j=1,nobj)
    write(11) cv, cvn, cvnm1
    w1 = timeint == 2
    write(11) w1
    if ( w1 ) write(11) (xnm1(j)%u, j=1,nobj)

    close(unit=11)

  end subroutine output_restart


! open files and input restart data

  subroutine open_files_and_restart

    integer :: j

    integer :: endstep, endpost
    real(dp) :: endtime

    if ( restart >= 1 ) then

!     restart: read solution from file

      open ( unit=11, form='unformatted', file='data.out' )

      read(11) endstep, endpost, endtime
      read(11) (mesh_b(j)%coor, j=1,nobj)
      read(11) (xn(j)%u, xnp1(j)%u, j=1,nobj)
      read(11) cv, cvn, cvnm1
      read(11) w1
      if ( w1 ) read(11) (xnm1(j)%u, j=1,nobj)

      close(unit=11)

      if ( ctime ) then
        time0 = endtime
        post0 = endpost
        step0 = endstep
      end if

      if ( compute_volume ) then

!       open volume file

        if ( ctime ) then
          open ( unit=12, file='volume.out', recl=reclv, status='old', &
            position='append' )
        else
          open ( unit=12, file='volume.out', recl=reclv, status='unknown' )
        end if

      end if

    else

!     fresh start

      if ( compute_volume ) then

!       open volume file

        open ( unit=12, file='volume.out', recl=reclv )

      end if

    end if

  end subroutine open_files_and_restart


! Problem definition of velocity/pressure system, sysmatrix and vectors

  subroutine problem_definition_create_sysmatrix_vectors

  integer :: j

  call problem_definition ( input_probdef, mesh, problem )

! Creation of system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )

! Creation of system matrix for velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! subscript for extracting velocities from mesh nodes on interface

  do j = 1, nobj
    call create_subscript ( mesh, problem, ssvel_b(j), surfaces=[6+j], &
      physqarr=[physqvel] )
  end do

  end subroutine problem_definition_create_sysmatrix_vectors


! Delete problem of velocity/pressure system, sysmatrix and vectors

  subroutine delete_problem_sysmatrix_vectors

    integer :: j

    call delete ( problem )
    call delete ( sol, rhsd )
    call delete ( sysmatrix )
    do j = 1, nobj
      call delete ( ssvel_b(j) )
    end do

  end subroutine delete_problem_sysmatrix_vectors


! Problem definition of interface system, sysmatrix and vectors

  subroutine problem_definition_interface_create_sysmatrix_vectors

    integer :: j

    do j = 1, nobj

!     problem definition of the interface

      call create_input_probdef ( mesh_b(j), input_probdef_b(j), nvec=2 )

      input_probdef_b(j)%elementdof(1)%a = 1
      input_probdef_b(j)%vec_elementdof(1)%a(:,1) = ndim  ! 3D vector
      input_probdef_b(j)%vec_elementdof(1)%a(:,2) = 1  ! error

      input_probdef_b(j)%probnr = 1+j

      call problem_definition ( input_probdef_b(j), mesh_b(j), problem_b(j) )

!     create system vectors (solution and right-hand side)

      call create_sysvector ( problem_b(j), sol_b(j) )
      call create_msysvector ( problem_b(j), rhsd_b(j,:) )

!     storage for position at tn-1, tn, tnp1 in vector

      call create ( problem_b(j), xn(j), vec=1 )
      call create ( problem_b(j), xnp1(j), vec=1 )

      if ( timeint == 2 ) then
        call create ( problem_b(j), xnm1(j), vec=1 )
      end if

!     Initialize xnp1 (which will be the xn of the new time step)

      xnp1(j)%u = reshape ( transpose(mesh_b(j)%coor), [ndim*mesh_b(j)%nnodes] )

!     storage for velocity in vector

      call create ( problem_b(j), vel_b(j), vec=1 )

!     create system matrix

      call create_sysmatrix_structure ( sysmatrix_b(j), mesh_b(j), &
        problem_b(j) )

      call create_sysmatrix_data ( sysmatrix_b(j) )

!     create oldvectors

      call create ( oldvectors_b(j), nsysvec=1, nvec=4, nprob=1 )

      oldvectors_b(j)%v(1)%p => vel_b(j)  ! velocity u on the interface
      oldvectors_b(j)%v(3)%p => xn(j)     ! position at tn
      if ( timeint == 2 ) then
        oldvectors_b(j)%v(4)%p => xnm1(j)  ! position at tn-1
      end if

    end do

  end subroutine problem_definition_interface_create_sysmatrix_vectors


! Delete problem of interface system, sysmatrix and vectors

  subroutine delete_problem_interface_create_sysmatrix_vectors

    integer :: j

    call delete ( coefficients_b )

    do j = 1, nobj
      call delete ( mesh_b(j) )
      call delete ( input_probdef_b(j) )
      call delete ( problem_b(j) )
      call delete ( sol_b(j) )
      call delete ( rhsd_b(j,:) )
      call delete ( sysmatrix_b(j) )
      call delete ( vel_b(j), xnp1(j), xn(j) )
      if ( timeint == 2 ) call delete ( xnm1(j) )
    end do

  end subroutine delete_problem_interface_create_sysmatrix_vectors


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

    use math_defs_m

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

    integer :: j

    type(real_array_2d_t), dimension(:), allocatable :: coor_boun

    character(len=*), parameter :: fmti = '(1x,a,i0,a)'

    call write_geometries_gmsh ( filename = 'mesh_geometries.msh', mesh=mesh, &
      binary=.true. )

    open ( unit=25, file='mesh_conform.geo' )

    write ( 25, fmti ) 'nobj = ', nobj, ';'
    write ( 25, '(/1x,a)' ) 'Include "mesh_options_3D.igo";'
    write ( 25, '(/1x,a)' ) 'Merge "mesh_geometries.msh";'
    write ( 25, '(/1x,a)' ) 'Include "particles_in_a_box_3D_full_2.igo";'

    close ( 25 )

!   save coordinates of the boundary mesh

    allocate ( coor_boun(nobj) )
    do j = 1, nobj
      allocate ( coor_boun(j)%a(mesh%surfaces(6+j)%nnodes,mesh%ndim) )
      coor_boun(j)%a = mesh%coor(mesh%surfaces(6+j)%nodes,:)
    end do

    call execute_command_line ( 'gmsh -3 -bin -o mesh_conform.msh &
      &mesh_conform.geo > outputmesh_conform.out' )

!   read mesh generated by gmsh

    call read_mesh_gmsh ( mesh_conform, filename='mesh_conform.msh', &
      physgeom=.true. )

!   match numbering/topology of surfaces with old surfaces in mesh, since
!   after remeshing with gmsh this has been changed. Use minimum distance
!   algorithm since coordinates have also been changed.

    do j = 1, nobj
      call add_to_mesh ( mesh_conform, matchingsurface=[6+j,6+j], replace=6+j, &
        matchingmesh=mesh, mindistance=.true. )
    end do

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

    else if ( periodic_channel ) then

!     local node numbering surface 4 to match local node numbering of surface 2

      call add_to_mesh ( mesh, matchingsurface=[4,2], replace=4, &
        displacement=[-lx,0._dp,0._dp] )

    end if

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

   do j = 1, nobj
     mesh%coor(mesh%surfaces(6+j)%nodes,:) = coor_boun(j)%a
     mesh%coor(mesh%surfaces(6+nobj+j)%nodes,:) = coor_boun(j)%a
   end do

   euler_step = .true.

   deallocate ( coor_boun )

 end subroutine remeshing

end program drop8
