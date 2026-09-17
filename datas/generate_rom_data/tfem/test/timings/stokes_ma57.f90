! Stokes problem on a unit square with Dirichlet boundary conditions.
! One freely moving and rotating object. Time-stepping.
! Sample in single point and sample on a line.
!
! This problem is the same as stokes9.f90 but now using the stokes elements
! available in the addon viscoelastic.

module subs_m

  use tfem_elem_m

  implicit none

  integer :: ndf = 9
  real(dp) :: rp = 1, xp(2) = 0

contains

! objectscoor defines the coordinates of the objects

  subroutine objectscoor ( objectnr, coor )
    use math_defs_m, only: pi
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

program stokes

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use subs_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    nc = 101,           & ! number of nodes in cross section
    numsteps = 5,       & ! number of time steps
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=100,             & ! number of elements in x
    ny=100                ! number of elements in y

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

! variables

  real(dp), parameter :: &
    eta = 1._dp     ! viscosity

  real(dp) :: xs(1,2), xc(nc,2)

  integer :: i, step, ip
  real(dp) :: up(3), deltat, unm1(3), un(3), pp(3), time
  character(len=20) :: filename

  integer :: num_build, num_solve
  real :: delta_build_avg, delta_solve_avg, delta

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  write(*,*)
  write(*,*) ' stokes_ma57 '
  write(*,*) ' ----------- '
  write(*,*)

  call tic

  call quadrilateral2d ( mesh, meshgen_options )

! one moving object

  rp = 0.1_dp ! radius of the circular object
  xp = [0.5_dp,0.85_dp] ! initial position of the center of the particle
  pp = [ (xp(i),i=1,2), 0._dp ] ! initial position and rotation

  call add_to_mesh ( mesh, object='coordinates', nnodes=20, &
    objectsub=objectscoor, objectcoornr=1 )

! write mesh (read by streamfunction computation)

  call fill_mesh_parts ( mesh )

  call toc ( 'mesh' )

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

  call toc ( 'problem definition' )

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

  time = 0

  call toc ( 'vectors' )

  do step = 1, numsteps

    if ( step >= 2 ) then

!     advance particle position (based on previous step)

      time = time + deltat

      ip = problem%constraints(1)%addnumdegfd(1)

      if ( step == 2 ) then
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

!     update reference coordinates of object
      call find_refcoor_objects ( mesh, object1=1 )

      call toc ( 'find_refcoor_objects' )

    end if

!   create system matrix

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric=.true. )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

    call create_sysmatrix_data ( sysmatrix )

    call toc ( 'create_sysmatrix' )

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
      coefficients=coefficients )

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      elemsub=elementc, addmatvec=.true. )

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    call toc ( 'build', delta )

    if ( step <= 1 ) then
      num_build = 0
      delta_build_avg = 0
    else
      num_build = num_build + 1
      delta_build_avg = delta_build_avg + delta
    end if

    call solve_system_ma57 ( sysmatrix, rhsd, sol )

    call toc ( 'solve', delta )

    if ( step <= 1 ) then
      num_solve = 0
      delta_solve_avg = 0
    else
      num_solve = num_solve + 1
      delta_solve_avg = delta_solve_avg + delta
    end if

    call delete ( sysmatrix )

  end do

  write ( *, '(/a,f0.5,a/)' ) 'Average CPU time = ', &
    delta_build_avg / num_build, ' seconds in build '
  write ( *, '(a,f0.5,a/)' ) 'Average CPU time = ', &
    delta_solve_avg / num_solve, ' seconds in solve '

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( coefficients )
  call delete ( oldvectors )

end program stokes
