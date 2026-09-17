! Drop consisting of a Newtonian fluid suspended in a Newtonian fluid.
! Navier-Stokes flow.
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
!
! Similar to drop5, but extended with inertia.

program drop9

  use tfem_m
  use math_defs_m
  use hsl_ma57_m
  use hsl_ma41_m
  use stokes_elements_m
  use interface_tracking_elements_m
  use poisson_gd_ma57_m
  use curvature_gd_ma57_m
  use update_mesh_nodes1_m
  use functions_m
  use io_utils_m
  use figplot_m
  use inertia_elements_m
  use projection_elements_m

  implicit none


! constants

  integer, parameter ::  &
    uintpl = 6,          & ! P2 velocity and interface position interpolation
    pintpl = 2,          & ! P1 pressure interpolation
    physqvel = 1,        & ! physical quantity nr of velocity
    physqpress = 2,      & ! physical quantity nr of pressure
    ndim = 2,            & ! dimension of space
    gauss_proj = 6,      & ! integration rule for the projection problem
    inttype_proj = 3,    & ! standard Gauss-Legendre (numerical table)
    ncomp_P2 = 8           ! number of components for P2 projection:
                           ! x,y-coordinates and x,y-velocities at t_n and t_nm1

! variables

  logical ::     &
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
    compute_error = .false.,& ! Compute error. Only makes sense for eta_e=eta_i
                              ! (no viscosity difference) and gammac=0 (no
                              ! surface tension) or if vfuncnr=0 (no flow)
                              ! or vfuncnr=3 (rotation).
    compute_volume = .false.  ! Compute volume and center of volume.

  integer ::  &
    restart=0,           & ! restart=0: no restart
                           ! restart=1: normal restart
                           ! restart=2: restart, but start with Euler step
    step0=0,             & ! initial step number
    post0=0,             & ! initial post number
    n=5,                 & ! number of elements in x and y direction of box
    nc=40,               & ! number of elements on the circle
    gauss  = 4,          & ! order of integration for triangular elements
    gaussb = 4,          & ! number of integration points of line elements
    vfuncnr=2,           & ! function number for the velocity field u
    choice_um = 1,       & ! choice for um in the motion of the interface:
                           !  0: um = 0
                           !  1: um = average velocity of all interface nodes
                           !  2: um = velocity of the center of volume
    vtkevery = 5 ,       & ! vtk file every vtkevery steps. -1: means none
    restart_every = 100, & ! write restart file every restart_every steps.
                           ! -1: means none
    timeint=2,           & ! time integration scheme for the interface.
    numtimesteps=400,    & ! number of time steps
    num_picard = 1,      & ! number of Picard iterations, followed by Newton
                           ! iterations
    maxiter = 100,       & ! maximum number of iterations allowed to obtain
                           ! convergence
    num_euler = 1          ! number of Euler time steps (at least one) for
                           ! the inertia terms.

  real(dp) ::  &
    ox = -0.5_dp,         & ! x-coordinate left lower corner
    oy = -0.5_dp,         & ! y-coordinate left lower corner
    lx = 1.0_dp,          & ! length of the domain
    ly = 1.0_dp,          & ! height of the domain
    xc(ndim) = [0.0_dp,0.0_dp],  & ! center of the initial circle
    rc = 0.1_dp,          & ! initial radius of the circle
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
    startup_time = 1._dp, & ! start-up time of the flow
    fac_min = 0.1_dp,     & ! factor setting dx_min for the monitor function
    fac_max = 10.0_dp,    & ! factor setting dx_max for the monitor function
    area_threshold = 1.39_dp,        & ! remeshing area threshold
    aspect_ratio_threshold = 1.39_dp,& ! remeshing aspect ratio threshold
    rho_e = 1.e0_dp,      & ! external fluid density
    rho_i = 1.e0_dp,      & ! internal fluid density
    threshold = 1.e-5_dp, & ! threshold for the iteration process
    gamma0 = 1.5_dp,      & ! coefficient for Gear scheme
    alpha0 = 2._dp,       & ! coefficient for Gear scheme
    alpha1 = -0.5_dp        ! coefficient for Gear scheme

! flow problem

  type(mesh_t) :: mesh_ext, mesh_int, mesh, mesh_conform
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, sol_iter, sol_n, sol_hat
  type(sysvector_t) :: rhsd, rhsd_stored, sol_nm1
  type(coefficients_t) :: coefficients(2)
  type(plot_options_t) :: plot_options
  type(subscript_t) :: ssvel_b, vel, press
  type(oldvectors_t) :: oldvectors_iter, oldvectors_hat
  type(vector_t), target :: meshvel
  type(subscriptvec_t) :: meshvelx, meshvely
  type(subscript_t) :: velx, vely

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

! temp definitions used to store old definitions before projection

  type(mesh_t), target :: mesh_temp
  type(sysvector_t) :: sol_n_temp, sol_nm1_temp
  type(vector_t) :: coords_n_temp, coords_nm1_temp
  type(problem_t), target :: problem_temp

! ALE mesh motion problem

  type(problem_t), target ::  problem_update_mesh
! Mesh coordinates at previous times
  real(dp), allocatable, dimension(:,:) :: meshcoor_n, meshcoor_nm1
! Displacement of interface
  real(dp), allocatable, dimension(:,:) :: disp

! some more misc variables

  logical :: w1
  character(len=30) :: filename

  integer :: nnodes, nelem, elshape, step, post, c_cpl, c_per, c_flr, iter
  integer :: vertices(3)=[1,3,5]
  real(dp) :: xp(1,ndim), rp(1), dx_box, dx_part, resultsum(1+ndim), &
    diffvel, diffpress, mfac

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

  namelist /comppar/ periodic, flowrate, move_mesh_x, monitor_func, mesh_plot, &
    compute_error, compute_volume, n, nc, gauss, gaussb, choice_um, vfuncnr, &
    vtkevery, ctime, restart, restart_every, timeint, numtimesteps, &
    ox, oy, lx, ly, xc, rc, fac_min, fac_max, &
    eta_e, eta_i, gammac, factor, U_avg, startup_time, rs_up, is_up, &
    deltat, beta, area_threshold, aspect_ratio_threshold, rho_e, rho_i

  read ( unit=*, nml=comppar )


! Fill coefficients

  call create ( coefficients, ncoefi=600, ncoefr=600 )

! External fluid (element group 1)

  coefficients(1)%i = 0
  coefficients(1)%i(1:11) = &
  [ uintpl,   pintpl,     0,     0,         0,  &
    physqvel, physqpress, 0,     0,     gauss,  &
    gaussb ]
  coefficients(1)%i(40) = 3 ! standard Gauss-Legendre (numerical tables)
  coefficients(1)%i(48) = 1 ! use mesh velocity for ALE formulation

  coefficients(1)%r = 0
  coefficients(1)%r(1) = eta_e
  coefficients(1)%r(6) = U_avg * ly
  coefficients(1)%r(8) = deltat
  coefficients(1)%r(19) = gammac
  coefficients(1)%r(151) = rho_e

! Internal fluid (element group 2)

  coefficients(2)%i = coefficients(1)%i

  coefficients(2)%r = 0
  coefficients(2)%r(1) = eta_i
  coefficients(2)%r(8) = deltat
  coefficients(2)%r(151) = rho_i

! Interface

  call create_coefficients ( coefficients_b, ncoefi=100, ncoefr=50 )

  coefficients_b%i = 0
  coefficients_b%i(1) = uintpl
  coefficients_b%i(3) = gaussb
  coefficients_b%i(5) = 3  ! standard Gauss-Legendre (numerical tables)
  coefficients_b%i(6) = 2  ! velocity u on the interface given by a vector
  coefficients_b%i(10) = 1 ! tangential grid velocity using Poisson problem
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

! Problem definition and create sysmatrix, vectors

  call problem_definition_create_sysmatrix_vectors

! monitor function

  if ( monitor_func ) allocate ( f_monitor(mesh_b%nnodes) )

! problem definition of the interface

  call create_input_probdef ( mesh_b, input_probdef_b, nvec=2 )

  input_probdef_b%elementdof(1)%a = 1
  input_probdef_b%vec_elementdof(1)%a(:,1) = ndim  ! 2D vector
  input_probdef_b%vec_elementdof(1)%a(:,2) = 1  ! error

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

! Create aspect ratio arrays

  call create_aspect_ratio

! define some arrays for the ALE mesh position at old times

  allocate ( meshcoor_n(mesh%nnodes,ndim), meshcoor_nm1(mesh%nnodes,ndim) )
  meshcoor_n = mesh%coor
  meshcoor_nm1 = mesh%coor

! initialize cv, ... to xc because they are all written to restart

  cv = xc; cvn = xc; cvnm1 = xc

! open files and restart

  call open_files_and_restart

! time stepping

  post = post0
  ts = startup_time

  do step = 1, numtimesteps

    write(*,'(a,i0)') 'step = ', step0+step

    time = time0 + step * deltat

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

!   copy mesh coordinates

    meshcoor_nm1 = meshcoor_n
    meshcoor_n = mesh%coor

  end do

  if ( compute_error ) close(unit=10)
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
  call delete ( ssvel_b )

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

  deallocate ( disp )
  if ( numtimesteps > 1 ) call delete ( problem_update_mesh )
  deallocate ( meshcoor_n, meshcoor_nm1 )

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

  subroutine build_solve_velocity_pressure

!   fill boundary conditions

    call fill_sysvector ( mesh, problem, sol, curve1=1, curve2=4, &
      physq=physqvel, vfunc=vfunc_2D, vfuncnr=vfuncnr )

    call fill_sysvector ( mesh, problem, sol, point=1, physq=physqpress, &
      value=0._dp )

!   Stokes velocity/pressure in external and internal fluid

    call build_system ( mesh, problem, sysmatrix, rhsd_stored, &
      elemsub=stokes_elem, mcoefficients=coefficients, buildmatrix=.false. )

!   Surface tension

    call add_boundary_elements ( mesh, problem, rhsd_stored, &
      elemsub=surface_tension_curve, curve=5, &
      coefficients=coefficients(1), physq=[physqvel] )

    iter = 0

    iterate: do

      iter = iter + 1

      call copy ( rhsd_stored, rhsd )

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=stokes_elem, mcoefficients=coefficients, buildvector=.false. )

      if ( periodic ) then

!       periodical condition on velocities in x-direction

        call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
          constraint1=c_per, elemsub=stokes_constr_node_conn, addmatvec=.true. )

        if ( flowrate ) then

!         impose flowrate

          call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
            constraint1=c_flr, elemsub=stokes_constr_flowr, addmatvec=.true., &
            coefficients=coefficients(1) )

        end if

      end if

!     Coupling of fluids

      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        constraint1=c_cpl, elemsub=stokes_constr_node_conn, addmatvec=.true. )

!     Advection term: first num_picard steps Picard iteration, then Newton
!     iteration

      if ( iter <= num_picard ) then
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_picard, mcoefficients=coefficients, &
          oldvectors=oldvectors_iter, physqcol=[physqvel], &
          physqrow=[physqvel], addmatvec=.true. )
      else
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_newton, mcoefficients=coefficients, &
          oldvectors=oldvectors_iter, physqcol=[physqvel], &
          physqrow=[physqvel], addmatvec=.true. )
      end if

!     du/dt

      if ( step <= num_euler .and. restart /= 1 ) then
!       implicit Euler
        sol_hat%u(vel%s) = sol_n%u(vel%s)
        mfac = 1._dp
      else
!       second-order implicit Gear after Euler time steps
        sol_hat%u(vel%s) = alpha0*sol_n%u(vel%s) + alpha1*sol_nm1%u(vel%s)
        mfac = gamma0
      end if

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=inertia_elem_dudt, mcoefficients=coefficients, &
        oldvectors=oldvectors_hat, physqcol=[physqvel], &
        physqrow=[physqvel], factormat=mfac, addmatvec=.true. )

      call check ( sysmatrix )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

      so%real_storage=rs_up
      so%integer_storage=is_up

      call solve_system_ma41 ( sysmatrix, rhsd, sol )

!     compute the difference between iteration steps

      diffvel = maxval( sqrt ( ( sol%u(vel%s) - sol_iter%u(vel%s) )**2 ) )
      diffpress = maxval( sqrt ( ( sol%u(press%s) - sol_iter%u(press%s) )**2 ) )
      write(*,'(a,i0,2(a,es12.4))')&
        'iter = ', iter, ' diffvel = ', diffvel, ' diffpress = ', diffpress
      call copy ( sol, sol_iter )

      if ( diffvel <= threshold ) exit iterate
      if ( iter == maxiter ) then
        write(*,'(a,i0)')&
          'Maximum number of iterations reached = ', maxiter
        stop
      end if

    end do iterate

    call copy ( sol_n, sol_nm1 )
    call copy ( sol, sol_n )

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

    if ( timeint == 1 .or. &
               ( step == 1 .and. timeint == 2 .and. restart /= 1 ) ) then

!     Euler

      call copy ( xnp1, xn )

      coefficients_b%i(11) = 1

!     update the nodes of the mesh

      call update_mesh_nodes_and_mesh_velocity ( order=1 )

    else if ( timeint == 2 ) then

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

    call build_solve_velocity_pressure

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

    call find_bounds_blocks ( mesh )

    if ( order == 1 ) then

!     Mesh velocity (BDF1)

      meshvel%u = reshape ( transpose ( &
                  ( mesh%coor - meshcoor_n ) / deltat ), [ndim*mesh%nnodes] )

    else

!     mesh velocity (BDF2)

      meshvel%u = reshape ( transpose ( &
        ( 1.5_dp*mesh%coor - 2*meshcoor_n + 0.5_dp*meshcoor_nm1 ) / deltat ), &
             [ndim*mesh%nnodes] )

    end if

    if ( move_mesh_x ) then

!     additional ALE movement of the mesh in x-direction

      meshvel%u(meshvelx%s) = meshvel%u(meshvelx%s) + cvv(1)

    end if

  end subroutine update_mesh_nodes_and_mesh_velocity


! output various data

  subroutine output_data

    use output_fields1_m

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

    call output_fields ( mesh, problem, coefficients, sol, post )

  end subroutine output_data


! output restart data

  subroutine output_restart

    write(*,'(a,i0,a,es10.3)') &
      'Output restart data at step ', step0+step, ', time =', time0+step*deltat

    call write_mesh ( mesh, filename='mesh.out' )

    open ( unit=11, form='unformatted', file='data.out' )

    write(11) step0 + step, post, time0 + step*deltat
    write(11) mesh_b%coor
    write(11) xn%u, xnp1%u
    write(11) cv, cvn, cvnm1
    write(11) meshcoor_n, meshcoor_nm1
    write(11) sol_n%u, sol_nm1%u, sol_iter%u
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

!     restart: read solution from file

      open ( unit=11, form='unformatted', file='data.out' )

      read(11) endstep, endpost, endtime
      read(11) mesh_b%coor
      read(11) xn%u, xnp1%u
      read(11) cv, cvn, cvnm1
      read(11) meshcoor_n, meshcoor_nm1
      read(11) sol_n%u, sol_nm1%u, sol_iter%u
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

!   Velocity/pressure problem definition in external and internal fluid

    call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=2 )

    input_probdef%vec_elementdof(1)%a(:,1) = ndim      ! velocity
    input_probdef%vec_elementdof(1)%a(:,2) = 0         ! pressure
    input_probdef%vec_elementdof(1)%a(vertices,2) = 1
    input_probdef%vec_elementdof(1)%a(:,3) = 1         ! scalar
    input_probdef%vec_elementdof(1)%a(:,4) = ndim*(ndim+1)/2
                                                       ! symmetric tensor
    input_probdef%vec_elementdof(1)%a(:,5) = ncomp_P2  ! projection

    input_probdef%vec_elementdof(2)%a = input_probdef%vec_elementdof(1)%a

    input_probdef%physq = [physqvel,physqpress]
    input_probdef%probnr = 1

    if ( periodic ) then

!     dirichlet on lower and upper wall

      call define_essential ( mesh, input_probdef, curve1=1, physq=physqvel )
      call define_essential ( mesh, input_probdef, curve1=3, physq=physqvel )

!     velocity periodic in x-direction

      call define_constraint ( mesh, input_probdef, physq=physqvel, &
        curve1=2, curve2=4, discretization='collocation', excludecurves=[1,3], &
        num=c_per )

      if ( flowrate ) then

!       impose flowrate

        call define_constraint ( mesh, input_probdef, physq=physqvel, &
          curve1=2, nglobalc=1, num=c_flr )

      end if

    else

!     dirichlet on all flowcell boundaries

      call define_essential ( mesh, input_probdef, curve1=1, curve2=4, &
        physq=physqvel )

    end if

    call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

!   Internal fluid - external fluid coupling through Lagrangian multipliers

    call define_constraint ( mesh, input_probdef, curve1=5, curve2=6, &
      discretization='collocation', physq=physqvel, num=c_cpl )

!   problem definition

    call problem_definition ( input_probdef, mesh, problem )

!   Creation of system vectors (solution and right-hand side)

    call create_sysvector ( problem, sol, sol_iter, sol_n )
    call create_sysvector ( problem, sol_nm1, sol_hat )
    call create_sysvector ( problem, rhsd, rhsd_stored )
    sol_n%u = 0._dp
    sol_iter%u = sol_n%u

!   Create mesh velocity vectors

    call create ( problem, meshvel, vec=1 )
    call create_subscript_vector ( mesh, problem, meshvelx, vec=1, degfd=1 )
    call create_subscript_vector ( mesh, problem, meshvely, vec=1, degfd=2 )
    call create_subscript ( mesh, problem, velx, physqarr=[1], degfd=1 )
    call create_subscript ( mesh, problem, vely, physqarr=[1], degfd=2 )

!   Creation of system matrix for velocity/pressure problem

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

    call create_sysmatrix_data ( sysmatrix )

!   subscript for extracting velocities from mesh nodes on interface

    call create_subscript ( mesh, problem, ssvel_b, curves=[5], &
      physqarr=[physqvel] )
    call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )
    call create_subscript ( mesh, problem, press, physqarr=[physqpress] )

!   Oldvectors

    call create_oldvectors ( oldvectors_hat, nsysvec=1 )
    call create_oldvectors ( oldvectors_iter, nsysvec=1, nprob=1, nvec=1 )
    oldvectors_hat%s(1)%p => sol_hat
    oldvectors_iter%s(1)%p => sol_iter
    oldvectors_iter%p(1)%p => problem
    oldvectors_iter%v(1)%p => meshvel

  end subroutine problem_definition_create_sysmatrix_vectors


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

    call create ( problem, coords_n_temp, vec=1 )
    call create ( problem, coords_nm1_temp, vec=1 )

    coords_n_temp%u = reshape ( transpose ( meshcoor_n ), [size(meshcoor_n)] )
    coords_nm1_temp%u = reshape ( transpose ( meshcoor_nm1 ), &
                                  [size(meshcoor_nm1)] )

!   copy current definitions to temporary definitions used for the projection
!   problem

    call copy ( mesh, mesh_temp )
    call fill_mesh_parts ( mesh_temp )
    call copy ( sol_n, sol_n_temp )
    call copy ( sol_nm1, sol_nm1_temp )
    call problem_definition ( input_probdef, mesh, problem_temp )

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

!   Change of local node numbering of curve 5 to match local node numbering of
!   curve 6

!    call add_to_mesh ( mesh, matchingcurve=[6,5], replace=6 )

   call fill_mesh_parts ( mesh )

!  delete various structures, arrays

   call delete_old_problems
   call delete_aspect_ratio
   call delete ( problem_update_mesh )

!  recreate the structures, arrays based on the new mesh

   call problem_definition_create_sysmatrix_vectors
   call create_aspect_ratio

!  project the old solutions onto the new mesh and fill in the values

   call project_old_solutions_P2

!  move boundary nodes to previous position (gmsh does not preserve quadratic
!  element shape).

   mesh%coor(mesh%curves(5)%nodes,:) = coor_boun
   mesh%coor(mesh%curves(6)%nodes,:) = coor_boun

!  delete the temporary definitions

   call delete ( coords_n_temp, coords_nm1_temp )
   call delete ( mesh_temp )
   call delete ( sol_n_temp )
   call delete ( sol_nm1_temp )
   call delete ( problem_temp )

   deallocate ( coor_boun )

 end subroutine remeshing

 subroutine delete_old_problems

   call delete ( problem )
   call delete ( sysmatrix )
   deallocate ( meshcoor_n, meshcoor_nm1 )
   call delete ( sol, sol_iter, sol_n )
   call delete ( sol_nm1, sol_hat )
   call delete ( rhsd, rhsd_stored )
   call delete ( oldvectors_hat, oldvectors_iter )
   call delete ( meshvel )
   call delete ( meshvelx, meshvely )
   call delete ( velx, vely )
   call delete ( ssvel_b )
   call delete ( vel, press )
   call delete ( input_probdef )

 end subroutine delete_old_problems

 subroutine project_old_solutions_P2

   type(input_probdef_t) :: input_probdef_proj
   type(problem_t) :: problem_proj
   type(sysmatrix_t) :: sysmatrix_proj
   type(sysvector_t) :: rhsd_proj(ncomp_P2)
   type(sysvector_t) :: sol_proj(ncomp_P2)
   type(oldvectors_t) :: oldvectors_proj
   type(coefficients_t) :: coefficients_proj
   type(vector_t), target :: vec_proj
   type(lu_ma57_t) :: lu_proj

   integer :: i

!  create vector for projection
   call create ( problem_temp, vec_proj, vec=5 )

!  transfer the coordinates at t_n
   call transfer_data ( mesh_temp, problem_temp, &
     vector1=coords_n_temp, vector2=vec_proj, degfd2=[1,2] )

!  transfer the velocities at t_n
   call transfer_data ( mesh_temp, problem_temp, &
     sysvector1=sol_n_temp, vector2=vec_proj, physq1=[1], degfd1=[1,2], &
     degfd2=[3,4] )

!  transfer the coordinates at t_nm1
   call transfer_data ( mesh_temp, problem_temp, &
     vector1=coords_nm1_temp, vector2=vec_proj, degfd2=[5,6] )

!  transfer the velocities at t_nm1
   call transfer_data ( mesh_temp, problem_temp, &
     sysvector1=sol_nm1_temp, vector2=vec_proj, physq1=[1], degfd1=[1,2], &
     degfd2=[7,8] )

!  problem definition for projection
   call create_input_probdef ( mesh, input_probdef_proj )

   input_probdef_proj%elementdof(1)%a = 1
   input_probdef_proj%elementdof(2)%a = 1
   call problem_definition ( input_probdef_proj, mesh, problem_proj )

!  create system vectors (solution and right-hand side)
   call create ( problem_proj, sol_proj )
   call create ( problem_proj, rhsd_proj )

!  create system matrix
   call create_sysmatrix_structure ( sysmatrix_proj, mesh, &
     problem_proj, symmetric=.true. )
   call create_sysmatrix_data ( sysmatrix_proj )

!  fill coefficients
   call create_coefficients ( coefficients_proj, ncoefi=100, ncoefr=50 )
   coefficients_proj%i = 0
   coefficients_proj%i(1:4) = [ uintpl, ncomp_P2, uintpl, uintpl ]
   coefficients_proj%i(10) = gauss_proj
   coefficients_proj%i(40) = inttype_proj
   coefficients_proj%r = 0

!  oldvectors
   call create_oldvectors ( oldvectors_proj, nvec=1, nprob=1, nmesh=1 )
   oldvectors_proj%v(1)%p => vec_proj
   oldvectors_proj%p(1)%p => problem_temp
   oldvectors_proj%m(1)%p => mesh_temp

!  build (assemble) matrix and vector from elements
   call build_system ( mesh, problem_proj, sysmatrix_proj, &
     msysvector=rhsd_proj, elemsub=projection_elem, &
     coefficients=coefficients_proj, oldvectors=oldvectors_proj )

!  solve the projection problem
   do i = 1, ncomp_P2
     call add_effect_of_essential_to_rhs ( problem_proj, sysmatrix_proj, &
       sol_proj(i), rhsd_proj(i) )
     call solve_system_ma57 ( sysmatrix_proj, rhsd_proj(i), sol_proj(i), &
       lu=lu_proj )
   end do

!  create meshes and solution vectors at t_n and t_nm1 and update the
!  old coordinates

   allocate ( meshcoor_n(mesh%nnodes,ndim), meshcoor_nm1(mesh%nnodes,ndim) )
   meshcoor_n(:,1) = sol_proj(1)%u
   meshcoor_n(:,2) = sol_proj(2)%u
   meshcoor_nm1(:,1) = sol_proj(5)%u
   meshcoor_nm1(:,2) = sol_proj(6)%u

!  update the old velocities
   sol_n%u(velx%s) = sol_proj(3)%u
   sol_n%u(vely%s) = sol_proj(4)%u
   sol_nm1%u(velx%s) = sol_proj(7)%u
   sol_nm1%u(vely%s) = sol_proj(8)%u

   sol_iter%u = sol_n%u

 end subroutine project_old_solutions_P2

end program drop9
