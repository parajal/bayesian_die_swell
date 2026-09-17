! 2D Stokes problem on a square domain with a single particle.
! Brownian forces using fluctuating hydrodynamics.
! Weak Lagrange multiplier implementation.
! Movement of the particle.

module subs_m

  use tfem_elem_m
  use math_defs_m

  implicit none

  integer :: ndf = 9
  integer :: ndfb = 3, nodalpb = 3
  integer :: ndflb = 2  ! number of degrees of freedom of the Lag. Mult.
  integer :: intplvel = 8
  real(dp) :: rp = 1, xpc(2) = 0

contains

! objectscoor defines the coordinates of the objects

  subroutine objectscoor ( objectnr, coor )
    integer, intent(in) :: objectnr
    real(dp), dimension(:,:), intent(inout) :: coor

    integer :: np, i
    real(dp) :: p(size(coor,1))

    np = size(coor,1)
    p = [(2*pi/np*(i-1),i=1,np)]

    coor(:,1) = rp * sin(p) + xpc(1)
    coor(:,2) = rp * cos(p) + xpc(2)

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

    integer :: object, i, j
    real(dp) :: phi(1,ndf), xr(1,2), r(2), psi(1,ndflb), x(nodalpb,2)
    real(dp) :: theta(1,ndfb), dtheta(1,ndfb,1), dxdxi(1,2), curvel(1)


    object = problem%constraints(constr)%object

!   reference coordinates of the integration point (node) of the element (elem)

    xr(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)
    r = mesh%objects(object)%coor_int(node,:,elem) - xpc

!   shape function of the velocity in the integration point

    if ( intplvel == 10 ) then
      call shape_quad_Q1isoQ2 ( xr, phi )
    else
      call shape_quad_Q2 ( xr, phi )
    end if

!   shape function of the Langrangian multiplier in the integration point

    call shape_line_P1 ( mesh%objects(object)%xig(node:node,1), psi )

!   shape function of the quadratic curve in the integration point

    call shape_line_P2 ( mesh%objects(object)%xig(node:node,1), theta, &
      dtheta(:,:,1) )

!   compute deformed element

    call get_coordinates_object ( mesh, elem, x, object )

    call isoparametric_deformation_curve ( x, dtheta(:,:,1), dxdxi, curvel )

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
                                   curvel(1) * mesh%objects(object)%wg(node)
        end do
      end do
      elemmat(ndflb+1:2*ndflb,ndf+1:2*ndf) = elemmat(1:ndflb,1:ndf)

      do i = 1, ndflb
        elemmatadd(i,:) = psi(1,i) * [ -1._dp,  0._dp,  r(2) ] * &
                                   curvel(1) * mesh%objects(object)%wg(node)
        elemmatadd(i+ndflb,:) = psi(1,i) * [  0._dp, -1._dp, -r(1) ] * &
                                   curvel(1) * mesh%objects(object)%wg(node)
      end do

    end if

  end subroutine elementc

end module subs_m

program brownian_motion3

  use tfem_m
  use generalized_stokes_elements_m
  use hsl_ma57_m
  use io_utils_m
  use figplot_m
  use subs_m

  implicit none

! constants

  integer(si), parameter :: &
    seed = 30             ! seed for ziggurat

  integer, parameter :: &
!    uintpl = 10,        & ! Q1isoQ2 velocities
    uintpl = 8,        & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    nx=40,              & ! number of elements in x
    ny=40,              & ! number of elements in y
    ran = 3,            & ! ziggurat: 2=uniform, 3=normal distribution
    everystep = 40,     & ! plot every 40 steps
    nsteps=1000

  real(dp), parameter :: &
    eta = 1._dp, &  ! viscosity
    deltat=1e-2_dp   ! time step

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh_particle
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients

  integer :: step, i, elem, gauss, nelem
  real(dp) :: up(3), pp(3), time, pp1(3)
  character(len=20) :: filename

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=300, ncoefr=250 )

  select case ( uintpl )
  case(10)
    gauss = 2  ! 2x2 integration of sub quad elements
  case default
    gauss = 3  ! 3x3 integration of quads
  end select

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  if ( uintpl == 10 ) then
!   divide the Gauss integration into four subquads (volume) and two line
!   intervals (boundary)
    coefficients%i(32:33) = [ 2, 2 ]
  end if

  coefficients%i(256) = ran  ! random number generator

  coefficients%r(1) = eta
  coefficients%r(2:) = 0
  coefficients%r(208) = 1 ! kT

  call write_coefficients ( coefficients, filename='coefficients.out' )

! set interpolation in element

  intplvel = uintpl

! initialize rng

  if ( ran == 2 .or. ran == 3 ) then
!   ziggurat
    call zigset ( seed )
  end if

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%lx = 8
  meshgen_options%ly = 8
  meshgen_options%ox = -4
  meshgen_options%oy = -4

  call quadrilateral2d ( mesh, meshgen_options )

! create mesh for particle boundary

  nelem = nx / 4 * 2

  call mesh_skeleton ( mesh_particle, nnodes=2*nelem, nelem=nelem, elshape=2, &
    ndim=2 )

  rp = 1.0_dp ! radius of the circular object
  xpc = [2.5_dp,0._dp] ! initial position of the center of the particle
  pp = [ (xpc(i),i=1,2), 0._dp ] ! initial position and rotation

  call objectscoor ( 1, mesh_particle%coor )

  do elem = 1, mesh_particle%nelem
    mesh_particle%topology(1)%a(:,elem) = [ 2*elem-1, 2*elem, 2*elem + 1 ]
  end do
  mesh_particle%topology(1)%a(3,mesh_particle%nelem) = 1 ! close circle

! one object

  call add_to_mesh ( mesh, object='mesh', objectmesh=mesh_particle, &
    topology=.true., intrule=1, nsubint=5 )

  if ( uintpl == 10 ) then
!   change elshape to 34 for isoparametric mapping
    mesh%element(:)%elshape = 34
  end if

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh3.fig' )
  plot_options%objectpointcolor=4
  plot_options%objectpointsize=0.4
  call plot_objects ( plot_options, mesh, 'mesh3.fig', append=.true. )

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
    discretization='weak', elementdof=[2,0,2], naddunknowns=3 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

  sol%u = 0

  open ( unit=10, file='out3', recl=300 )

  time = 0; step = 0

  do

!   create system matrix

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric=.true. )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

    call create_sysmatrix_data ( sysmatrix )

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
      coefficients=coefficients, buildvector=.false. )

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      elemsub=elementc, addmat=.true., coefficients=coefficients )

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=generalized_stokes_rhs_fh, physqrow=[1], physqcol=[1], &
      coefficients=coefficients, buildmatrix=.false. )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    call solve_system_ma57 ( sysmatrix, rhsd, sol )

    call delete ( sysmatrix )

    call get_sysvector_constraint ( mesh, problem, sol, constraint=1, &
      addunknowns=.true., u=up )

    pp1 = pp + up * sqrt(deltat)  ! Brownian dynamics

    if ( any ( abs(pp1) >= 3._dp ) ) then
!     reject move
      print *, 'reject step', step
      cycle
    else
      pp = pp1
    end if

    xpc = pp(1:2)

    time = time + deltat
    step = step + 1

    write(10,*) time, pp

!   update coordinates of object
    call objectscoor ( objectnr=1, coor=mesh%objects(1)%coor )

!   update reference coordinates of object
    call find_refcoor_objects ( mesh, object1=1 )

    if ( mod(step,everystep) == 0 ) then

      write(filename,'(a,i3.3,a)') 'particle', step/everystep, '.fig'

      plot_options%objectpointcolor=0
      call plot_objects ( plot_options, mesh, filename, object1=1 )

    end if

    if ( step == nsteps ) exit

  end do

  close(10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh, mesh_particle )
  call delete ( sol, rhsd )
  call delete ( coefficients )

end program brownian_motion3
