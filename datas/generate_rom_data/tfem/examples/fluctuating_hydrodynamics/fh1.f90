! 2D Stokes problem on a square domain with a single particle.
! Brownian forces using fluctuating hydrodynamics.
! Boundary fitted mesh.
! No movement of the particle in order to estimate the diffusion coefficient
! from the sample and check the Stokes-Einstein-Sutherland relation.

module subs_m

  use tfem_elem_m

  implicit none

  integer :: ndf = 9
  real(dp) :: xpc(2) = 0

contains

! elementc is the element subroutine for the constraints on the curve

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

    integer :: curve, nodenr
    real(dp) :: r(2)

    curve = problem%constraints(constr)%geometry1

!   coordinates of node

    nodenr = mesh%curves(curve)%nodes(node)
    r = mesh%coor(nodenr,:) - xpc

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
!     A =  [ 1  0 ]
!          [ 0  1 ]
!
!   and the matrix A_add of the additional unknowns (up,vp,omega)
!
!     A_add =  [ -1  0  ry ]
!              [  0 -1 -rx ]
!

    if ( matrix ) then

      elemmat(1,1) = 1
      elemmat(1,2) = 0
      elemmat(2,1) = 0
      elemmat(2,2) = 1
      elemmatadd(1,:) = [ -1._dp,  0._dp,  r(2) ]
      elemmatadd(2,:) = [  0._dp, -1._dp, -r(1) ]

    end if

  end subroutine elementc

end module subs_m

program fh1

  use tfem_m
  use generalized_stokes_elements_m
  use hsl_ma57_m
  use io_utils_m
  use figplot_m
  use subs_m

  implicit none

! constants

  integer(si), parameter :: &
    seed = 16             ! seed for ziggurat

  integer, parameter :: &
!    uintpl = 10,        & ! Q1isoQ2 velocities
    uintpl = 8,        & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    ran = 3,            & ! ziggurat: 2=uniform, 3=normal distribution
    nsteps=20000

  real(dp), parameter :: &
    eta = 1._dp      ! viscosity

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: pressure
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(lu_ma57_t) :: lu

  integer :: step, gauss
  real(dp) :: up(3)


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

! initialize rng

  if ( ran == 2 .or. ran == 3 ) then
!   ziggurat
    call zigset ( seed )
  end if

! read mesh

  call read_mesh ( mesh, filename='mesh1.out' )

  if ( uintpl == 10 ) then
!   change elshape to 34 for isoparametric mapping
    mesh%element(:)%elshape = 34
  end if

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh1.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=12, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=10, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=7, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=4, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

! define constraints on curve

  call define_constraint ( mesh, input_probdef, curve1=13, physq=1, &
    discretization='collocation', nodedof=2, naddunknowns=3 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

  sol%u = 0

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients, buildvector=.false. )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    elemsub=elementc, addmat=.true., coefficients=coefficients )

  open ( unit=10, file='out1', recl=300 )

  do step = 1, nsteps

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=generalized_stokes_rhs_fh, physqrow=[1], physqcol=[1], &
      coefficients=coefficients, buildmatrix=.false. )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu )

    call get_sysvector_constraint ( mesh, problem, sol, constraint=1, &
      addunknowns=.true., u=up )

    write(10,*) up

  end do

  close(10)

! post-processing

  call create_vector ( problem, pressure, vec=3 )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  plot_options%scalevector=0.1
  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    sysvector=sol )

  call plot_color_fill ( plot_options, mesh, problem, 'vx_color.fig', &
    sysvector=sol, physq=1, degfd=1 )
  call plot_mesh ( plot_options, mesh, 'vx_color.fig', append=.true. )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vx_contour.fig', sysvector=sol, physq=1, degfd=1 )
  call plot_mesh ( plot_options, mesh, 'vx_contour.fig', append=.true. )

  call plot_color_fill ( plot_options, mesh, problem, 'vy_color.fig', &
    sysvector=sol, physq=1, degfd=2 )
  call plot_mesh ( plot_options, mesh, 'vy_color.fig', append=.true. )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vy_contour.fig', sysvector=sol, physq=1, degfd=2 )
  call plot_mesh ( plot_options, mesh, 'vy_contour.fig', append=.true. )

  call plot_color_fill ( plot_options, mesh, problem, 'pressure_color.fig', &
    vector=pressure )
  call plot_mesh ( plot_options, mesh, 'pressure_color.fig', append=.true. )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )
  call plot_mesh ( plot_options, mesh, 'pressure_contour.fig', append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( pressure )
  call delete ( coefficients )
  call delete ( lu )
  call delete ( oldvectors )

end program fh1
