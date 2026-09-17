! Drop consisting of a Newtonian fluid suspended in a Newtonian fluid.
! Stokes flow (no inertia).
! Flow imposed either:
!   1) on the whole boundary using a function, or
!   2) on the upper and lower boundary using a function and
!      periodic boundary conditions in x- and z-direction (periodic=.true.)
!   3) periodic boundary conditions in x-direction and on the other boundaries
!      imposed by a function (periodic_channel=.true.)
! ALE with interface tracking and optional displacement of the mesh in
! x-direction to keep the center of volume at the same position in the mesh.
! Optionally impose flowrate if periodic=.true.
! Optionally compute the error if the exact solution is known.
! Optionally compute the volume and center of volume.
! Three-dimensional model.
! Similar to drop2, but now with:
!   P1 shape function for the position of the interface.

program drop2b

  use tfem_m
  use math_defs_m
  use hsl_ma57_m
  use hsl_ma41_m
  use stokes_elements_m
  use interface_tracking_elements_m
  use poisson_gd_ma57_m
  use update_mesh_nodes2_m
  use functions_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  integer, parameter ::  &
    ndim = 3               ! dimension of space

! constants flow problem

  integer, parameter ::  &
    uintpl = 6,          & ! P2 velocity interpolation
    pintpl = 2,          & ! P1 pressure interpolation
    physqvel = 1,        & ! physical quantity nr of velocity
    physqpress = 2         ! physical quantity nr of pressure

! constants interface tracking

  integer, parameter :: &
    xintpl = 2,       & ! P1 position
    vintpl = 6          ! P2 velocity


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
    ctime = .false.,        & ! continue time in output after restart
    mesh_plot = .false.,    & ! write mesh plot files/info (curves, mesh etc.)
    compute_error = .false.,& ! Compute error. Only makes sense for eta_e=eta_i
                              ! (no viscosity difference) and gammac=0 (no
                              ! surface tension) or if vfuncnr=0 (no flow)
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
    xc(ndim) = [0.0_dp,0.0_dp,0.0_dp],  & ! center of the initial sphere
    rc = 0.1_dp,          & ! radius of the initial sphere
    eta_e = 1.0_dp,       & ! external fluid viscosity
    eta_i = 1.0_dp,       & ! internal fluid viscosity
    gammac = 0.5_dp,      & ! surface tension coefficient
    factor = 100._dp,     & ! factor for the interface grid deformation
    rs_up  = 1.4_dp,      & ! real_storage velocity-pressure LU HSL
    is_up  = 1.5_dp,      & ! integer_storage velocity-pressure LU HSL
    U_avg  = 1.0_dp,      & ! imposed average velocity for flowrate=.true.
    time0  = 0.0_dp,      & ! initial time
    beta = 1.0_dp,        & ! factor beta in SUPG
    deltat = 7.e-3_dp       ! time step

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
  type(subscript_t) :: ssvel_b

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
  type(vector_t) :: xnp1, error, vec_b
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
  real(dp), allocatable, dimension(:,:) :: disp

! some more misc variables

  logical :: w1
  character(len=30) :: filename

  integer :: nnodes, nelem, elshape, step, post, c_cpl, c_perx, c_perz, c_flr
  integer :: vertices(4)=[1,3,5,10]
  real(dp) :: xp(1,ndim), rp(1), dx_box, dx_part, resultsum(1+ndim)

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

  namelist /comppar/ periodic, periodic_channel, flowrate, move_mesh_x, &
    mesh_plot, compute_error, compute_volume, n, nc, gauss, gaussb, &
    choice_um, vfuncnr, vtkevery, ctime, restart, restart_every, &
    timeint, numtimesteps, ox, oy, oz, lx, ly, lz, xc, rc, eta_e, eta_i, &
    gammac, factor, U_avg, rs_up, is_up, deltat, beta

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
  coefficients(1)%r(19) = gammac

! Internal fluid (element group 2)

  coefficients(2)%i = coefficients(1)%i

  coefficients(2)%r = 0
  coefficients(2)%r(1) = eta_i

! Interface

  call create_coefficients ( coefficients_b, ncoefi=100, ncoefr=50 )

  coefficients_b%i = 0
  coefficients_b%i(1) = xintpl
  coefficients_b%i(3) = gaussb
  coefficients_b%i(5) = 3  ! standard Gauss-Legendre (numerical tables)
  coefficients_b%i(6) = 2  ! velocity u on the interface given by a vector
  coefficients_b%i(10) = 1  ! tangential grid velocity using Poisson problem
  coefficients_b%i(11) = timeint
  coefficients_b%i(14) = vintpl
  coefficients_b%i(16) = 1 ! use shapefunc for velocity for element shape

  coefficients_b%r = 0
  coefficients_b%r(7) = factor
  coefficients_b%r(8) = beta
  coefficients_b%r(9) = deltat


  if ( restart == 0 ) then

!   generate fluid mesh using gmsh

    xp(1,:) = xc
    rp(1) = rc
    dx_box = lx/n
    dx_part = 2*pi*rc/nc

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

  if ( mesh_plot ) then

    call printinfo ( mesh_b, printlevel=4 )

    call write_mesh_vtk ( mesh_b, filename='mesh_b.vtk' )

  end if


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

  call define_constraint ( mesh, input_probdef, surface1=7, surface2=8, &
    discretization='collocation', physq=physqvel, num=c_cpl )

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

  call create_subscript ( mesh, problem, ssvel_b, surfaces=[7], &
    physqarr=[physqvel] )


! problem definition of the interface

  call create_input_probdef ( mesh_b, input_probdef_b, nvec=2 )

  input_probdef_b%elementdof(1)%a = [1,0,1,0,1,0]
  input_probdef_b%vec_elementdof(1)%a(:,1) = ndim  ! 2D vector
  input_probdef_b%vec_elementdof(1)%a(:,2) = 1  ! error

  input_probdef_b%probnr = 2

  call problem_definition ( input_probdef_b, mesh_b, problem_b )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem_b, sol_b )
  call create_msysvector ( problem_b, rhsd_b )

! create vector for storage of position component in all nodes

  call create ( problem_b, vec_b, vec=2 )

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

  call create ( oldvectors_b, nsysvec=2, nvec=4, nprob=1 )

  oldvectors_b%v(1)%p => vel_b  ! velocity u on the interface
  oldvectors_b%v(3)%p => xn     ! position at tn
  if ( timeint == 2 ) then
    oldvectors_b%v(4)%p => xnm1  ! position at tn-1
  end if


! allocate arrays for the ALE displacement problem

  allocate ( disp(mesh_b%nnodes,ndim) )


! inititalize cv, ... to xc because they are all written to restart

  cv = xc; cvn = xc; cvnm1 = xc


! open files and restart

  call open_files_and_restart

! time stepping

  post = post0

  do step = 1, numtimesteps

    write(*,'(a,i0)') 'step = ', step0+step

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
  call delete ( sol )
  call delete ( rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( ssvel_b )

! interface tangential motion problem

  call delete ( problem_pois )
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
  call delete ( vec_b )
  if ( compute_error ) call delete ( error )
  deallocate ( disp )

! ALE mesh motion problem

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
      constraint1=c_cpl, elemsub=stokes_constr_node_conn, addmatvec=.true. )

!   Surface tension

    call add_boundary_elements ( mesh, problem, rhsd, &
      elemsub=surface_tension_surface, surface=7, &
      coefficients=coefficients(1), physq=[physqvel] )

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    solver_options_ma57%real_storage=rs_up
    solver_options_ma57%integer_storage=is_up

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_ma57  )

  end subroutine build_solve_velocity_pressure


! build and solve for xnp1 on interface

  subroutine build_solve_interface

    integer :: i

!   Velocity of the interface nodes (reshaped)
    real(dp), dimension(ndim,mesh_b%nnodes) :: drop_vel

!   Poisson problem on interface mesh for making nodes "diffuse" and
!   redistribute uniformly

    call poisson_gd_ma57 ( mesh_b, problem_pois, coefficients_b, &
      sysmatrix_pois, sol_pois, elementdof=problem_b%elnumdegfd )

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

    oldvectors_b%s(2)%p => sol_b

    do i = 1, ndim

      so%integer_storage=1.8_dp

      call solve_system_ma41 ( sysmatrix_b, rhsd_b(i), sol_b, lu=lu_b, &
        solver_options=so )

!     compute position component in all nodes

      call derive_vector ( mesh_b, problem_b, vec_b, &
        elemsub=position_component_deriv, coefficients=coefficients_b, &
        oldvectors=oldvectors_b )

      call transfer_data ( mesh_b, problem1=problem_b, &
        vector1=vec_b, vector2=xnp1, degfd1=[1], degfd2=[i] )

    end do

    call delete(lu_b)

  end subroutine build_solve_interface


! solve flow and move the interface position (time discretized step)

  subroutine solve_flow_and_move_interface

    if ( timeint == 1 .or. &
               ( step == 1 .and. timeint ==2 .and. restart /= 1 ) ) then

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

end program drop2b
