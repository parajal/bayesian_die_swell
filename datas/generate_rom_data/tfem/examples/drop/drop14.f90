! Axisymmetric drop consisting of a Newtonian fluid suspended in a Newtonian
! fluid. Note in this example (x,y) coordinates should be read as (z,r).
! Stokes flow (no inertia) containing insoluble surfactants on the
! interface.
! Flow imposed either:
!   1) on the whole boundary using a function, or
!   2) on the upper and lower boundary using a function and
!      periodic boundary conditions in the z-direction (periodic=.true.).
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

program drop14

  use tfem_m
  use math_defs_m
  use hsl_ma57_m
  use hsl_ma41_m
  use stokes_elements_m
  use interface_tracking_elements_m
  use surface_convection_diffusion_elements_m
  use poisson_gd_ma57_m
  use curvature_gd_ma57_m
  use update_mesh_nodes13_m
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
    coorsys = 1,         & ! axisymmetric coordinate system
    ndim = 2               ! dimension of space

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
    vfuncnr=4,           & ! function number for the velocity field u
    choice_um = 1,       & ! choice for um in the motion of the interface:
                           !  0: um = 0
                           !  1: um = average velocity of all interface nodes
                           !  2: um = velocity of the center of volume
    vtkevery = 5 ,       & ! vtk file every vtkevery steps. -1: means none
    restart_every = 100, & ! write restart file every restart_every steps.
                           ! -1: means none
    timeint=2,           & ! time integration scheme
    timeint_sc=2,        & ! time integration scheme for surfactant
                           ! concentration
    numtimesteps=400       ! number of time steps

  real(dp) ::  &
    ox = -0.5_dp,         & ! x-coordinate left lower corner
    oy = -0.5_dp,         & ! y-coordinate left lower corner
    lx = 1.0_dp,          & ! length of the domain
    ly = 1.0_dp,          & ! height of the domain
    xc(ndim) = [0.0_dp,0.0_dp],  & ! center of the initial circle
    rc = 0.1_dp,          & ! initial radius of the circle
    eta_e = 1.0_dp,       & ! external fluid viscosity
    eta_i = 1.0_dp,       & ! internal fluid viscosity
    gammac = 0.1_dp,      & ! initial surface tension
    c0 = 0.1_dp,          & ! initial surfactant concentration
    beta_st = 0.1_dp,     & ! surface tension beta coefficient
    Dcoef = 1.0_dp,       & ! surfactant concentration diffusion
    factor = 100._dp,     & ! factor for the interface grid deformation
    rs_up  = 1.4_dp,      & ! real_storage velocity-pressure LU HSL
    is_up  = 1.5_dp,      & ! integer_storage velocity-pressure LU HSL
    U_avg  = 1.0_dp,      & ! imposed average velocity for flowrate=.true.
    time0  = 0.0_dp,      & ! initial time
    beta = 0.5_dp,        & ! factor beta in SUPG
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
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(coefficients_t) :: coefficients(2)
  type(plot_options_t) :: plot_options
  type(solver_options_ma41_t) :: solver_options_ma41
  type(subscript_t) :: ssvel_b

! curvature problem

  type(problem_t) :: problem_curv
  type(sysmatrix_t) :: sysmatrix_curv

! interface tangential motion problem

  type(problem_t), target :: problem_pois
  type(sysmatrix_t) :: sysmatrix_pois
  type(sysvector_t), target :: sol_pois

! interface tracking problem

  type(mesh_t), target :: mesh_b, mesh_c
  type(input_probdef_t) :: input_probdef_b
  type(problem_t), target :: problem_b
  type(sysmatrix_t) :: sysmatrix_b
  type(sysvector_t), target :: sol_b
  type(sysvector_t) :: rhsd_b
  type(vector_t) :: xnp1, error
! Interface mesh coordinates at previous times. NOTE: the coordinates are
! defined with respect to the initial coordinate system.
  type(vector_t), target :: xn, xnm1, um
  type(vector_t), target :: vel_b
  type(oldvectors_t) :: oldvectors_b
  type(coefficients_t) :: coefficients_b
  type(lu_ma41_t) :: lu_b
  type(solver_options_ma41_t) :: so

! surfactant concentration problem

  type(input_probdef_t) :: input_probdef_c
  type(problem_t), target :: problem_c
  type(sysmatrix_t) :: sysmatrix_c
  type(sysvector_t), target :: sol_c, sol_cm1
  type(sysvector_t) :: rhsd_c
  type(oldvectors_t) :: oldvectors_c
  type(coefficients_t) :: coefficients_c

! ALE mesh motion problem

  type(problem_t), target ::  problem_update_mesh
  type(real_array_2d_t) :: disp

! some more misc variables

  logical :: w1, w2, euler_step, euler_step_sc, print_averagec=.true.
  character(len=30) :: filename

  integer :: nnodes, nelem, elshape, step, post, c_cpl, c_per, c_flr
  integer :: vertices(3)=[1,3,5]
  real(dp) :: xp(1,ndim), rp(1), dx_box, dx_part, resultsum(2)
  real(dp) :: resc(1), resarea(1)

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
! NOTE: only the first dimension is relevant since the rest should be
!       always zero for the axisymmetric case.
  real(dp) :: cv, cvn, cvnm1

! center of volume velocity
  real(dp) :: cvv


! namelist for input of variables; read from standard input

  namelist /comppar/ periodic, flowrate, move_mesh_x, monitor_func, mesh_plot, &
    compute_error, compute_volume, n, nc, gauss, gaussb, choice_um, vfuncnr, &
    vtkevery, ctime, restart, restart_every, timeint, timeint_sc, &
    numtimesteps, ox, oy, lx, ly, xc, rc, fac_min, fac_max, &
    eta_e, eta_i, beta_st, gammac, c0, Dcoef, factor, U_avg, rs_up, is_up, &
    deltat, beta, area_threshold, aspect_ratio_threshold

  read ( unit=*, nml=comppar )


! Fill coefficients

  call create ( coefficients, ncoefi=600, ncoefr=600 )

! External fluid (element group 1)

  coefficients(1)%i = 0
  coefficients(1)%i(1:11) = &
  [ uintpl,   pintpl,     0,     0,         0,  &
    physqvel, physqpress, 0,     0,     gauss,  &
    gaussb ]
  coefficients(1)%i(23) = coorsys ! axisymmetric
  coefficients(1)%i(40) = 3 ! standard Gauss-Legendre (numerical tables)

  coefficients(1)%r = 0
  coefficients(1)%r(1) = eta_e
  coefficients(1)%r(6) = U_avg * ly
  coefficients(1)%r(25) = beta_st
  coefficients(1)%r(26) = gammac
  coefficients(1)%r(27) = c0

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
  coefficients_b%i(10) = 1 ! tangential grid velocity using Poisson problem
  coefficients_b%i(11) = timeint
  coefficients_b%i(12) = 1 ! build matrix for both degrees of freedom
  coefficients_b%i(13) = coorsys

  coefficients_b%r = 0
  coefficients_b%r(7) = factor
  coefficients_b%r(8) = beta
  coefficients_b%r(9) = deltat

! Surfactant concentration

  call create_coefficients ( coefficients_c, ncoefi=100, ncoefr=50 )

  coefficients_c%i = 0
  coefficients_c%i(1) = uintpl
  coefficients_c%i(3) = gaussb
  coefficients_c%i(5) = 3  ! standard Gauss-Legendre (numerical tables)
  coefficients_c%i(6) = timeint_sc
  coefficients_c%i(7) = 1 ! coorsys

  coefficients_c%r = 0
  coefficients_c%r(4) = Dcoef
  coefficients_c%r(5) = deltat

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

  nnodes = mesh%curves(2)%nnodes
  nelem = mesh%curves(2)%nelem
  elshape = mesh%curves(2)%element%elshape

  call mesh_skeleton ( mesh_b, nnodes, nelem, elshape, ndim )

  mesh_b%topology(1)%a=mesh%curves(2)%topology(:,:,1)
  mesh_b%coor=mesh%coor(mesh%curves(2)%nodes,:)

  call add_to_mesh ( mesh_b, point=[mesh_b%coor(1,:)] )
  call add_to_mesh ( mesh_b, point=[mesh_b%coor(mesh_b%nnodes,:)] )

  call fill_mesh_parts ( mesh_b )
  call copy ( mesh_b, mesh_c )
  call fill_mesh_parts ( mesh_c )

  if ( mesh_plot ) then

    call plot_mesh ( plot_options, mesh_b, 'mesh_b.fig' )

    call printinfo ( mesh_b, printlevel=4 )

  end if


! Velocity/pressure problem definition in external and internal fluid

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = ndim      ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0         ! pressure
  input_probdef%vec_elementdof(1)%a(vertices,2) = 1
  input_probdef%vec_elementdof(1)%a(:,3) = 1         ! scalar
  input_probdef%vec_elementdof(1)%a(:,4) = ndim*(ndim+1)/2+coorsys !sym tensor

  input_probdef%vec_elementdof(2)%a = input_probdef%vec_elementdof(1)%a

  input_probdef%physq = [physqvel,physqpress]
  input_probdef%probnr = 1

  if ( periodic ) then

!   dirichlet on lower and upper wall

    call define_essential ( mesh, input_probdef, curve1=1, physq=physqvel, &
           degfd = [0,1] )
    call define_essential ( mesh, input_probdef, curve1=3, physq=physqvel, &
           degfd = [0,1] )
    call define_essential ( mesh, input_probdef, curve1=5, physq=physqvel )
    call define_essential ( mesh, input_probdef, curve1=8, physq=physqvel, &
           degfd = [0,1], excludecurves=[7] )

!   velocity periodic in x-direction

    call define_constraint ( mesh, input_probdef, physq=physqvel, &
      curve1=4, curve2=6, discretization='collocation', excludecurves=[1,3,5], &
      num=c_per )

    if ( flowrate ) then

!     impose flowrate

      call define_constraint ( mesh, input_probdef, physq=physqvel, curve1=4, &
        nglobalc=1, num=c_flr )

    end if

  else

!   dirichlet on all flowcell boundaries

    call define_essential ( mesh, input_probdef, curve1=1, physq=physqvel, &
           degfd = [0,1] )
    call define_essential ( mesh, input_probdef, curve1=3, physq=physqvel, &
           degfd = [0,1] )
    call define_essential ( mesh, input_probdef, curve1=4, curve2=6, &
           physq=physqvel )
    call define_essential ( mesh, input_probdef, curve1=8, physq=physqvel, &
           degfd = [0,1], excludecurves=[7] )

  end if

  call define_essential ( mesh, input_probdef, point=3, physq=physqpress )

! Internal fluid - external fluid coupling through Lagrangian multipliers

  call define_constraint ( mesh, input_probdef, curve1=2, curve2=7, &
    discretization='collocation', physq=physqvel, num=c_cpl )

! Problem definition and create sysmatrix, vectors

  call problem_definition_create_sysmatrix_vectors

! monitor function

  if ( monitor_func ) allocate ( f_monitor(mesh_b%nnodes) )

! problem definition of the interface

  call create_input_probdef ( mesh_b, input_probdef_b, nvec=2 )

  input_probdef_b%elementdof(1)%a = ndim ! position
  input_probdef_b%vec_elementdof(1)%a(:,1) = ndim  ! 2D vector
  input_probdef_b%vec_elementdof(1)%a(:,2) = 1  ! error

  input_probdef_b%probnr = 2

! define essential boundary conditions for the r component at
! the symmetry axis for the interface tracking problem
  call define_essential ( mesh_b, input_probdef_b, point=1, &
         degfd = [0,1] )
  call define_essential ( mesh_b, input_probdef_b, point=2, &
         degfd = [0,1] )

  call problem_definition ( input_probdef_b, mesh_b, problem_b )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem_b, sol_b )
  call create_sysvector ( problem_b, rhsd_b )

! storage for position at tn-1, tn, tnp1 in vector

  call create ( problem_b, xn, vec=1 )
  call create ( problem_b, xnp1, vec=1 )
  call create ( problem_b, um, vec=1 )

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

  allocate ( disp%a(mesh_b%nnodes,ndim) )

! Create aspect ratio arrays

  call create_aspect_ratio

! initialize cv, ... to xc because they are all written to restart

  cv = xc(1); cvn = xc(1); cvnm1 = xc(1)

! define surfactant concentration problem

  call define_problem_surfactant

! define initial concentration of surfactants

  sol_c%u = c0

! open files and restart

  call open_files_and_restart

! time stepping

  post = post0
  euler_step = .true.
  euler_step_sc = .true.

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

    call build_solve_surfactant

    if ( print_averagec ) then

!     surfactant concentration on the interface

      call integrate ( mesh_b, problem_c, resc, oldvectors=oldvectors_c, &
        elemsub=integrate_surface_convection_diffusion, &
        coefficients=coefficients_c )

!     interface area

      call integrate ( mesh_b, problem_c, resarea, oldvectors=oldvectors_c, &
        elemsub=integrate_surface_area_elem, coefficients=coefficients_c )

      print *, 'average concentration', resc/resarea

    end if

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

        write(filename,'(a,i4.4,a)') 'surfconc_',post,'.vtk'

        call write_scalar_vtk ( mesh_c, problem_c, sysvector=sol_c, &
          dataname='c', filename=filename )

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
  call delete ( mesh_b, mesh_c )
  call delete ( input_probdef_b )
  call delete ( problem_b )
  call delete ( sol_b )
  call delete ( rhsd_b )
  call delete ( sysmatrix_b )
  call delete ( vel_b, xnp1, xn )
  if ( timeint == 2 ) call delete ( xnm1 )
  if ( compute_error ) call delete ( error )

! ALE mesh motion problem

  deallocate  ( disp%a )
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

!     local node numbering of curve 6 to match local node numbering of curve 4

      call add_to_mesh ( mesh, matchingcurve=[6,4], replace=6, &
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
    write ( 25, '(/1x,a)' ) 'Include "particles_in_a_box_axissym.igo";'

    close ( 25 )

  end subroutine write_gmsh_parameters


! build and solve for velocity and pressure

  subroutine build_solve_velocity_pressure

!   fill boundary conditions

    call fill_sysvector ( mesh, problem, sol, curve1=1, curve2=8, &
      physq=physqvel, vfunc=vfunc_2D, vfuncnr=0 )

    call fill_sysvector ( mesh, problem, sol, curve1=4, curve2=6, &
      physq=physqvel, vfunc=vfunc_2D, vfuncnr=vfuncnr )

    call fill_sysvector ( mesh, problem, sol, point=3, physq=physqpress, &
      value=0._dp )

!   Stokes velocity/pressure in external and internal fluid

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, mcoefficients=coefficients )

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
      elemsub=surface_tension_surfactants_curve, curve=2, &
      coefficients=coefficients(1), physq=[physqvel], &
      oldvectors=oldvectors_c )

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    solver_options_ma41%real_storage=rs_up
    solver_options_ma41%integer_storage=is_up

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_ma41  )

  end subroutine build_solve_velocity_pressure


! build and solve for xnp1 on interface

  subroutine build_solve_interface

!   Velocity of the interface nodes (reshaped)
    real(dp), dimension(ndim,mesh_b%nnodes) :: drop_vel

    real(dp) :: dx_min, dx_max

    call fill_sysvector ( mesh_b, problem_b, sol_b, point=1, &
      value=0._dp, degfd=2 )
    call fill_sysvector ( mesh_b, problem_b, sol_b, point=2, &
      value=0._dp, degfd=2 )

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
      coefficients_b%r(4) = sum(drop_vel(1,:))/mesh_b%nnodes

    else if ( choice_um == 2 ) then

!     velocity of the center of volume

      coefficients_b%r(4) = cvv

    end if

!   build (assemble) matrix and vector from elements

    call build_system ( mesh_b, problem_b, sysmatrix_b, sysvector=rhsd_b, &
!      elemsub=interface_tracking_elem_lagrange, coefficients=coefficients_b, &
      elemsub=interface_tracking_elem_supg_impl, coefficients=coefficients_b, &
      oldvectors=oldvectors_b, order='DN' )


    so%integer_storage=1.8_dp

    call solve_system_ma41 ( sysmatrix_b, rhsd_b, sol_b, lu=lu_b, &
      solver_options=so )

    call transfer_data ( mesh_b, problem1=problem_b, &
      sysvector1=sol_b, vector2=xnp1 )


    call delete(lu_b)

  end subroutine build_solve_interface


! solve flow and move the interface position (time discretized step)

  subroutine solve_flow_and_move_interface

    if ( timeint == 1 .or. &
               ( euler_step .and. timeint == 2 .and. restart /= 1 ) ) then

      if ( timeint == 2 ) euler_step = .false.

!     Euler

      call copy ( xnp1, xn )

      coefficients_b%i(11) = 1

!     update the nodes of the mesh

      call update_mesh_nodes_from_displacement ( order=1 )

    else if ( timeint == 2 ) then

!     second order Gear with prediction

      call copy ( xn, xnm1 )
      call copy ( xnp1, xn )

!     2nd order prediction of interface coordinates
      mesh_b%coor =  &
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

    disp%a = mesh_b%coor - mesh%coor(mesh%curves(2)%nodes,:)

    if ( move_mesh_x .or. choice_um == 2 ) then

!     update center of volume

      cvnm1 = cvn
      cvn = cv

      call integrate ( mesh_b, problem_b, resultsum, &
        elemsub=center_of_volume, coefficients=coefficients_b )

      cv = resultsum(2)/resultsum(1)

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

      disp%a(:,1) = disp%a(:,1) - ( cv - xc(1) )

    end if

!   update mesh nodes

    call update_mesh_nodes ( mesh, problem_update_mesh, disp=disp )

    mesh_c%coor = mesh%coor(mesh%curves(2)%nodes,:)

  end subroutine update_mesh_nodes_from_displacement


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
    write(11) sol_c%u
    write(11) cv, cvn, cvnm1
    w1 = timeint == 2
    w2 = timeint_sc == 2
    write(11) w1, w2
    if ( w1 ) write(11) xnm1%u
    if ( w2 ) write(11) sol_cm1%u

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
      read(11) sol_c%u
      read(11) cv, cvn, cvnm1
      read(11) w1, w2
      if ( w1 ) read(11) xnm1%u
      if ( w2 ) read(11) sol_cm1%u

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

!   Creation of system vectors (solution and right-hand side)

    call create_sysvector ( problem, sol, rhsd )

!   Creation of system matrix for velocity/pressure problem

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric=.false. )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

    call create_sysmatrix_data ( sysmatrix )

!   subscript for extracting velocities from mesh nodes on interface

    call create_subscript ( mesh, problem, ssvel_b, curves=[2], &
      physqarr=[physqvel] )

  end subroutine problem_definition_create_sysmatrix_vectors


! Delete problem of velocity/pressure system, sysmatrix and vectors

  subroutine delete_problem_sysmatrix_vectors

    call delete ( problem )
    call delete ( sol, rhsd )
    call delete ( sysmatrix )
    call delete ( ssvel_b )

  end subroutine delete_problem_sysmatrix_vectors


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
    real(dp), dimension(:), allocatable :: c_temp

    character(len=*), parameter :: fmti = '(1x,a,i0,a)'
    integer, parameter :: nobj = 1

    call write_geometries_gmsh ( filename = 'mesh_geometries.msh', mesh=mesh, &
      binary=.false. )

    open ( unit=25, file='mesh_conform.geo' )

    write ( 25, fmti ) 'nobj = ', nobj, ';'
    write ( 25, '(/1x,a)' ) 'Include "mesh_options_2D.igo";'
    write ( 25, '(/1x,a)' ) 'Merge "mesh_geometries.msh";'
    write ( 25, '(/1x,a)' ) 'Include "particles_in_a_box_axissym_2.igo";'

    close ( 25 )

!   save coordinates of the boundary mesh

    allocate ( coor_boun(mesh%curves(2)%nnodes,mesh%ndim) )
    coor_boun = mesh%coor(mesh%curves(2)%nodes,:)

!   save surfactant concentration

    allocate ( c_temp(mesh%curves(7)%nnodes ) )
    c_temp = sol_c%u

    call execute_command_line ( 'gmsh -2 -bin -o mesh_conform.msh &
      &mesh_conform.geo > outputmesh_conform.out' )

!   read mesh generated by gmsh

    call read_mesh_gmsh ( mesh_conform, filename='mesh_conform.msh', ndim=2, &
      physgeom=.true. )

    call add_to_mesh ( mesh_conform, matchingcurve=[7,7], replace=7, &
      matchingmesh=mesh, mindistance=.true. )

!   split mesh in external and internal mesh

    call mesh_convert ( mesh_conform, mesh_ext, only_groups=[1] )
    call mesh_convert ( mesh_conform, mesh_int, only_groups=[2] )

    call delete ( mesh_conform, mesh )

!   merge internal and external meshes

    call mesh_merge ( mesh_ext, mesh_int, mesh, nogroupmerge=.true. )

    call delete ( mesh_ext, mesh_int )

    if ( periodic ) then

!     local node numbering of curve 6 to match local node numbering of curve 4

      call add_to_mesh ( mesh, matchingcurve=[6,4], replace=6, &
        displacement=[-lx,0._dp] )

    end if

   call fill_mesh_parts ( mesh )

!  delete various structures, arrays

   call delete_problem_sysmatrix_vectors
   call delete_surfactant
   call delete_aspect_ratio
   call delete ( problem_update_mesh )

!  recreate the structures, arrays based on the new mesh

   call problem_definition_create_sysmatrix_vectors
   call create_aspect_ratio
   call define_problem_surfactant

!  move boundary nodes to previous position (gmsh does not preserve quadratic
!  element shape).

   mesh%coor(mesh%curves(2)%nodes,:) = coor_boun
   mesh%coor(mesh%curves(7)%nodes,:) = coor_boun

!  restore the surfactant concetration vector
   sol_c%u = c_temp

   euler_step = .true.
   euler_step_sc = .true.

   deallocate ( coor_boun )

 end subroutine remeshing

 subroutine define_problem_surfactant

!  surfactant concentration problem definition

   call create_input_probdef ( mesh_b, input_probdef_c, nvec=1, nphysq=1 )

   input_probdef_c%vec_elementdof(1)%a(:,1) = 1         ! surfactant
                                                        ! concentration

   input_probdef_c%physq = [1]
   input_probdef_c%probnr = 3

   call problem_definition ( input_probdef_c, mesh_b, problem_c )

!  Creation of system vectors (solution and right-hand side)

   call create_sysvector ( problem_c, sol_c, rhsd_c )
   if ( timeint_sc == 2 ) then
     call create_sysvector ( problem_c, sol_cm1 )
     sol_cm1%u = 0._dp
   end if

!  Creation of system matrix for velocity/pressure problem

   call create_sysmatrix_structure_base ( sysmatrix_c, mesh_b, problem_c, &
     symmetric=.false. )
   call finalize_sysmatrix_structure ( sysmatrix_c )

   call create_sysmatrix_data ( sysmatrix_c )

!  create oldvectors

   call create ( oldvectors_c, nsysvec=2, nvec=2, nprob=2, nmesh=1 )

   oldvectors_c%s(1)%p => sol_c      ! surfactant concentration at tn
   if (timeint_sc == 2 ) &
     oldvectors_c%s(2)%p => sol_cm1  ! surfactant concentration at tn-1
   oldvectors_c%v(1)%p => vel_b      ! velocity u on the interface
   oldvectors_c%v(2)%p => um         ! velocity of the grid on the interface
   oldvectors_c%m(1)%p => mesh_b
   oldvectors_c%p(1)%p => problem_c
   oldvectors_c%p(2)%p => problem_b

 end subroutine define_problem_surfactant

 subroutine build_solve_surfactant

    integer :: i

!   Velocity of the interface nodes (reshaped)
    real(dp), dimension(ndim,mesh_b%nnodes) :: drop_vel

    if ( timeint_sc == 1 .or. &
            ( euler_step_sc .and. timeint_sc == 2 .and. restart/=1 ) ) then

      if ( timeint_sc == 2 ) euler_step_sc = .false.

      coefficients_c%i(6) = 1

!     calculate grid velocity on the interface

      um%u = ( xnp1%u - xn%u ) / deltat

    else if ( timeint_sc == 2 ) then

      coefficients_c%i(6) = 2
      call copy ( sol_c, sol_cm1 )

!     calculate grid velocity on the interface

      um%u = ( 1.5 * xnp1%u - 2 * xn%u + 0.5 * xnm1%u ) / deltat

    else

      write (*,'(/a,i0/)') 'Error: invalid timeint_sc = ', timeint_sc
      stop

    end if

    drop_vel = reshape(vel_b%u,[ndim,mesh_b%nnodes]) ! reshape
    do i = 1, ndim
      coefficients_c%r(i) = sum(drop_vel(i,:))/mesh_b%nnodes
    end do

!   build (assemble) matrix and vector from elements

    call build_system ( mesh_b, problem_c, sysmatrix_c, &
      sysvector=rhsd_c, elemsub=surface_convection_diffusion_elem, &
      coefficients=coefficients_c, oldvectors=oldvectors_c )

    call solve_system_ma41 ( sysmatrix_c, rhsd_c, sol_c, solver_options=so )

 end subroutine build_solve_surfactant

 subroutine delete_surfactant

   call delete ( input_probdef_c )
   call delete ( problem_c )
   call delete ( sol_c )
   call delete ( rhsd_c )
   call delete ( sysmatrix_c )
   if ( timeint_sc == 2 ) call delete ( sol_cm1 )
   call delete ( oldvectors_c )

 end subroutine delete_surfactant

end program drop14
