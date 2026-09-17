! Stokes problem on a unit square with Dirichlet boundary conditions.
! One freely moving and rotating object. Time-stepping.
! Tecplot output

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

program stokes9

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use figplot_m
  use tec_utils_m
  use io_utils_m
  use subs_m

  implicit none

! constants

  integer, parameter :: &
    nc = 101,           & ! number of nodes in cross section
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=20,              & ! number of elements in x
    ny=20                 ! number of elements in y

  integer,parameter :: num_tec_data = 4  ! number of data sets in tecplot output

  real(dp), parameter :: &
    eta = 1._dp     ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t), target :: velocity, pressure, vorticity
  type(vector_t), target :: vel_x, vel_y
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors, oldvectors_tec
  type(coefficients_t) :: coefficients

  type(subscript_t) :: velx, vely

! variables

  integer ::  numsteps = 56    ! number of time steps


  real(dp) :: xs(1,2), xc(nc,2)

  integer :: i, istep, ip, k
  real(dp) :: up(3), deltat, unm1(3), un(3), pp(3), r(2), time

  character(len=20), dimension(num_tec_data) :: dataname

  logical :: append


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

! one object for sampling the velocity in a single point

  xs(1,:) = [0.5_dp,0.85_dp]

  call add_to_mesh ( mesh, object='coordinates', nnodes=1, coor=xs )

! one object for sampling the velocity on a line (diagonal)

  xc(:,1) = [ (i*0.01_dp, i=0,nc-1) ]
  xc(:,2) = [ (i*0.01_dp, i=0,nc-1) ]

  call add_to_mesh ( mesh, object='coordinates', nnodes=1000, coor=xc )

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )

! plot mesh and objects

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )
  plot_options%objectpointcolor=4
  plot_options%objectpointsize=0.4
  call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true. )

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

  call create_oldvectors ( oldvectors,     nsysvec=1 )
  call create_oldvectors ( oldvectors_tec, nvec=num_tec_data )

  call create_vector ( problem, velocity,  physq=1 )
  call create_vector ( problem, pressure,  vec=3 )
  call create_vector ( problem, vorticity, vec=3 )
  call create_vector ( problem, vel_x,     vec=3 )
  call create_vector ( problem, vel_y,     vec=3 )

! create subscript for sysvector

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )

! define data set for tecplot output

  dataname(1) = 'u'
  dataname(2) = 'v'
  dataname(3) = 'p'
  dataname(4) = 'vorticity'

  oldvectors_tec%v(1)%p => vel_x
  oldvectors_tec%v(2)%p => vel_y
  oldvectors_tec%v(3)%p => pressure
  oldvectors_tec%v(4)%p => vorticity

  deltat=1e-1_dp ! time step

!  numsteps = numsteps * 10
!  deltat   = deltat / 10.

  time = 0

  do istep = 1, numsteps

    if ( istep >= 2 ) then

!     advance particle position (based on previous step)

      time = time + deltat

      ip = problem%constraints(1)%addnumdegfd(1)

      if ( istep == 2 ) then
!       advance with forward Euler
        un = sol%u( problem%degfdperm(ip+1:ip+3,2) )
        pp = pp + deltat*un
        do k = 1, mesh%objects(1)%nnodes
          mesh%objects(1)%coor(k,:) = mesh%objects(1)%coor(k,:) + deltat*un(1:2)
        end do
      else
!       advance with 2nd order Adams-Bashforth
        unm1 = un
        un = sol%u( problem%degfdperm(ip+1:ip+3,2) )
        pp = pp + deltat*(3*un/2 - unm1/2)
        do k = 1, mesh%objects(1)%nnodes
          mesh%objects(1)%coor(k,:) = mesh%objects(1)%coor(k,:) + &
                                      deltat*(3*un(1:2)/2 - unm1(1:2)/2)
        end do
      end if

      xp = pp(1:2)   ! particle center

!     rotate and update collocation points, rotation angle = pp(3)

      do k = 1, mesh%objects(1)%nnodes
        r  = mesh%objects(1)%coor(k,:) - xp(:)
        mesh%objects(1)%coor(k,1) = r(1)*cos(pp(3)) - r(2)*sin(pp(3)) + xp(1)
        mesh%objects(1)%coor(k,2) = r(1)*sin(pp(3)) + r(2)*cos(pp(3)) + xp(2)
      end do

!     update coordinates of object
!      call objectscoor ( objectnr=1, coor=mesh%objects(1)%coor )

!     update reference coordinates of object

      call find_refcoor_objects ( mesh )

    end if

    print *, '------'
    print *, 'time = ', time
    print *, 'pp = ', pp

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

    call solve_system_ma57 ( sysmatrix, rhsd, sol )

    call delete ( sysmatrix )

!   extract info on particle

    ip = problem%constraints(1)%addnumdegfd(1)

    up = sol%u( problem%degfdperm(ip+1:ip+3,2) ) ! velocity and rotation rate

    print *, 'up = ', up

!   extract velocity

    call extract_physvector ( mesh, problem, sol, velocity )

    oldvectors%s(1)%p => sol

!   extract pressure and vorticity

    call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
      coefficients=coefficients, oldvectors=oldvectors )

    coefficients%i(13)=5
    call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
      coefficients=coefficients, oldvectors=oldvectors )

!   Tecplot output

    do i = 1, mesh%nnodes
      vel_x%u(i) = velocity%u(problem%vec_nodnumdegfd(i,velocity%vec) + 1)
      vel_y%u(i) = velocity%u(problem%vec_nodnumdegfd(i,velocity%vec) + 2)
    end do

    if (istep == 1) then
      append = .false.
    else
      append = .true.
    end if

    call write_mdata_tecplot (mesh, problem, 'cavity.plt', dataname, &
      oldvectors_tec, append=append, time=time )

    call write_objects_tecplot ( mesh, 'object', 'cavity_object.plt', &
      object1=1, append=append, time=time )

  end do

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( vel_x, vel_y )
  call delete ( coefficients )
  call delete ( oldvectors, oldvectors_tec )
  call delete ( velx, vely )

end program stokes9
