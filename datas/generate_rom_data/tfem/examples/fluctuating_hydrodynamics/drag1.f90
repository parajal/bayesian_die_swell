! 2D Stokes problem on a square domain with a single particle.
! Boundary fitted mesh.
! Set forces/torque to 1 to find the inverse drag coefficients.

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
      if ( first ) then
        elemvecadd = 1  ! set forces/torque to 1
      end if

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

program drag1

  use tfem_m
  use generalized_stokes_elements_m
  use hsl_ma57_m
  use io_utils_m
  use figplot_m
  use subs_m

  implicit none

! constants

  integer, parameter :: &
!    uintpl = 10,        & ! Q1isoQ2 velocities
    uintpl = 8,        & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2        ! physical quantity nr of the pressures

  real(dp), parameter :: &
    eta = 1._dp        ! viscosity

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(vector_t) :: velocity
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients
  type(lu_ma57_t) :: lu

  integer :: gauss
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

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

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

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol, lu )

  call get_sysvector_constraint ( mesh, problem, sol, constraint=1, &
    addunknowns=.true., u=up )

  print *, 'inverse drag coefficient:'
  print *, up

! post-processing

  plot_options%scalevector=10
  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    sysvector=sol )

! write binary file for reading by streamfunction

  call create_vector ( problem, velocity, physq=1 )
  call extract_physvector ( mesh, problem, sol, velocity )

  open(unit=10,form='unformatted',file='velocity_bin.out')

  write(10) velocity%u

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( lu )

end program drag1
