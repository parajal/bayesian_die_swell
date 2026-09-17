! Stokes problem on a unit square with Dirichlet boundary conditions.
! Lid-driven cavity flow.
! One freely moving and rotating object. Time-stepping.
! Grid deformation.
! This problem is similar to example stokes8, but now with grid deformation.

module subs_m

  use tfem_elem_m
  use math_defs_m

  implicit none

  integer :: ndf = 9
  real(dp) :: rp = 1, xp(2) = 0

contains

! objectscoor defines the coordinates of the objects

  subroutine objectscoor ( objectnr, coor )
    integer, intent(in) :: objectnr
    real(dp), dimension(:,:), intent(inout) :: coor

    integer :: np, i
    real(dp) :: p(size(coor,1))

    np = size(coor,1)
    p = [(2*pi/np*(i-1),i=1,np)]

    coor(:,1) = rp * sin(p) + xp(1)
    coor(:,2) = rp * cos(p) + xp(2)

  end subroutine objectscoor

! elementc is the element subroutine for the constraints on the objects

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

    integer :: object
    real(dp) :: phi(1,ndf), xr(1,2), r(2)

    object = problem%constraints(constr)%object

!   reference coordinates of the collocation point (node)

    xr(1,:) = mesh%objects(object)%refcoor(node,:)
    r = mesh%objects(object)%coor(node,:) - xp

!   shape function of the velocity in the collocation point

    call shape_quad_Q2 ( xr, phi )

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

    if ( matrix ) then

      elemmat(1,1:ndf) = phi(1,:)
      elemmat(1,ndf+1:) = 0
      elemmat(2,1:ndf) = 0
      elemmat(2,ndf+1:) = phi(1,:)
      elemmatadd(1,:) = [ -1._dp,  0._dp,  r(2) ]
      elemmatadd(2,:) = [  0._dp, -1._dp, -r(1) ]

    end if

  end subroutine elementc

end module subs_m

program grid_deformation6

  use tfem_m
  use hsl_ma57_m
  use meshgen_grid_deformation_m
  use deform_mesh_ma57_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m
  use subs_m
  use subs_grid_deformation6_m
  use timer_m


  implicit none

! constants

  integer, parameter :: &
    numsteps = 56,      & ! number of time steps
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=20,              & ! number of elements in x
    ny=20,              & ! number of elements in y
    max_gridsteps = 20    ! number of pseudo time steps in the grid deformation

  real(dp), parameter :: &
!   The (relative) size of the elements near the zero levelset (min_f) and far
!   away from the the zero levelset (max_f).
    min_f = 0.01,  &
    max_f = 0.4,  &
    eta = 1._dp     ! viscosity


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma57_t) :: solver_options

  type(mesh_t) :: mesh_gd
  type(problem_t) :: problem_gd
  type(coefficients_t) :: coefficients_gd
  type(sysmatrix_t) :: sysmatrix_gd
  type(lu_ma57_t) :: lu

  integer :: obj_gd
  real(dp), allocatable, dimension(:) :: f_monitor


! variables

  logical :: buildmatrix

  integer :: i, ip
  real(dp) :: up(3), deltat, unm1(3), un(3), pp(3), time
  character(len=20) :: filename


  timer = .false.

! set values in subs module

  lmin_f = min_f; lmax_f = max_f

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

  call tic

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

! one moving object

  rp = 0.1_dp ! radius of the circular object
  xp = [0.5_dp,0.85_dp] ! initial position of the center of the particle
  pp = [ (xp(i),i=1,2), 0._dp ] ! initial position and rotation

  call add_to_mesh ( mesh, object='coordinates', nnodes=20, &
    objectsub=objectscoor, objectcoornr=1 )

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )

! plot initial mesh

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh_initial.fig' )
  plot_options%objectpointcolor=4
  plot_options%objectpointsize=0.4
  call plot_objects ( plot_options, mesh, 'mesh_initial.fig', append=.true. )

  call toc('mesh')


! copy mesh to grid def. mesh

  call mesh_convert_copy ( mesh, mesh_gd )

! use sblocks to speed up the search for interpolation points

  call add_to_mesh ( mesh_gd, sblocks=[4*nx,4*ny] )

  call fill_mesh_parts ( mesh_gd )

! fill coefficients for the Poisson problem for grid deformation

  call create_coefficients ( coefficients_gd, ncoefi=100, ncoefr=50 )

  coefficients_gd%i = 0
  coefficients_gd%i(1) = uintpl
  coefficients_gd%i(10) = gauss

  coefficients_gd%r = 0

! allocate monitor function f

  allocate ( f_monitor(mesh%nnodes) )



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
    discretization='collocation', nodedof=2, naddunknowns=3 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=1, degfd=1, value=1._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  deltat=1e-1_dp ! time step

  call create_vector ( problem, velocity, physq=1 )

  time = 0
  buildmatrix = .true.
  obj_gd = 0  ! set to non-existing object initially

  do i = 1, numsteps

    if ( i >= 2 ) then

!     advance particle position (based on previous step)

      time = time + deltat

      ip = problem%constraints(1)%addnumdegfd(1)

      if ( i == 2 ) then
!       advance with forward Euler
        un = sol%u( problem%degfdperm(ip+1:ip+3,2) )
        pp = pp + deltat*un
      else
!       advance with 2nd order Adams-Bashforth
        unm1 = un
        un = sol%u( problem%degfdperm(ip+1:ip+3,2) )
        pp = pp + deltat*(3*un/2 - unm1/2)
      end if

      xp = pp(1:2)

!     update coordinates of object
      call objectscoor ( objectnr=1, coor=mesh%objects(1)%coor )

    end if

    print *, '------'
    print *, 'time = ', time
    print *, 'pp = ', pp

!   preform grid deformation

!   radius and center position

    rpl = rp
    cpl = xp

!   set monitor function f

    call set_monitor_function ( mesh_gd, f_monitor )

!   perform the grid deformation by pseudo time stepping

    if ( i >= 2 ) buildmatrix = .false.

    call deform_mesh_ma57 ( mesh_gd, obj_gd, problem_gd, coefficients_gd, &
      sysmatrix_gd, f_monitor, numgridsteps=max_gridsteps, lu=lu, &
      project_boundary=project_boundary, buildmatrix=buildmatrix, &
      plotfigs=.false. )

!   update coordinates of the mesh

    mesh%coor = mesh_gd%objects(obj_gd)%coor

    call find_bounds_blocks ( mesh )

!   update reference coordinates of object
    call find_refcoor_objects ( mesh, object1=1 )


!   create system matrix

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric=.true. )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

    call create_sysmatrix_data ( sysmatrix )

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
      coefficients=coefficients )

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      elemsub=elementc, addmatvec=.true. )

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    solver_options%real_storage = 1.5

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options )

    call delete ( sysmatrix )

!   extract info on particle

    ip = problem%constraints(1)%addnumdegfd(1)

    up = sol%u( problem%degfdperm(ip+1:ip+3,2) ) ! velocity and rotation rate

    print *, 'up = ', up

    write(filename,'(a,i3.3,a)') 'v_contour', i, '.fig'

    plot_options%objectpointcolor=0
    call extract_physvector ( mesh, problem, sol, velocity )
    call plot_color_contour ( plot_options, mesh, problem, filename, &
      vector=velocity, degfd=2 )
    call plot_objects ( plot_options, mesh, filename, append=.true., object1=1 )

    write(filename,'(a,i3.3,a)') 'mesh', i, '.fig'

    call plot_mesh ( plot_options, mesh, filename )
    plot_options%objectpointcolor=4
    plot_options%objectpointsize=0.4
    call plot_objects ( plot_options, mesh, filename, append=.true. )

  end do

  close(unit=10)

  oldvectors%s(1)%p => sol

  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5
  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a fig file for plotting

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )
  call plot_objects ( plot_options, mesh, 'velocity.fig', append=.true. )

  call plot_color_fill ( plot_options, mesh, problem, 'velocity_color.fig', &
    vector=velocity, degfd=1 )
  call plot_objects ( plot_options, mesh, 'velocity_color.fig', append=.true. )

  call plot_color_fill ( plot_options, mesh, problem, 'vorticity_color.fig', &
    vector=vorticity )
  call plot_objects ( plot_options, mesh, 'vorticity_color.fig', append=.true. )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( coefficients )
  call delete ( oldvectors )

  call delete ( problem_gd )
  call delete ( mesh_gd )
  call delete ( coefficients_gd )
  call delete ( sysmatrix_gd )

end program grid_deformation6
