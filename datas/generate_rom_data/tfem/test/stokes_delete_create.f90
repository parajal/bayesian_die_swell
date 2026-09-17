! Stokes problem on a unit square with Dirichlet boundary conditions.
! Lid-driven cavity flow.
! Two fixed objects.

module subs_m

  use tfem_elem_m
  use math_defs_m

  implicit none

  integer :: ndf = 9
  real(dp) :: rp(2) = 1, xp(2,2) = 0

contains

! objectscoor defines the coordinates of the objects

  subroutine objectscoor ( objectnr, coor )
    integer, intent(in) :: objectnr
    real(dp), dimension(:,:), intent(inout) :: coor

    integer :: np, i
    real(dp) :: p(size(coor,1))

    np = size(coor,1)
    p = [(2*pi/np*(i-1),i=1,np)]

    coor(:,1) = rp(objectnr) * sin(p) + xp(objectnr,1)
    coor(:,2) = rp(objectnr) * cos(p) + xp(objectnr,2)

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
    real(dp) :: phi(1,ndf), xr(1,2)

    object = problem%constraints(constr)%object

!   reference coordinates of the collocation point (node)

    xr(1,:) = mesh%objects(object)%refcoor(node,:)

!   shape function of the velocity in the collocation point

    call shape_quad_Q2 ( xr, phi )

    if ( vector ) then

      elemvec = 0 ! no rhs in constraint

    end if

!   two constraints: u=0, v=0, depending on both u and v
!   the matrix A therefore becomes
!
!     A =  [ phi    0 ]
!          [   0  phi ]

    if ( matrix ) then

      elemmat(1,1:ndf) = phi(1,:)
      elemmat(1,ndf+1:) = 0
      elemmat(2,1:ndf) = 0
      elemmat(2,ndf+1:) = phi(1,:)

    end if

  end subroutine elementc

end module subs_m

program stokes6

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use io_utils_m
!  use figplot_m
  use io_utils_m
  use subs_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=20,              & ! number of elements in x
    ny=20                 ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp     ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh1
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
  !type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients


! fill coefficients

  call fill_coefficients

!  call write_coefficients ( coefficients, filename='coefficients.out' )

  call delete ( coefficients)

  call fill_coefficients

! create mesh

  call fill_mesh

  call delete(mesh)

  call fill_mesh

  call copy(mesh,mesh1)
  call delete(mesh)
  call copy(mesh1,mesh)
  call delete(mesh1)

  call fill_mesh_parts(mesh)

  !stop

! plot mesh and objects

!  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
!  call plot_mesh ( plot_options, mesh, 'mesh.fig' )
!  plot_options%objectpointcolor=4
!  plot_options%objectpointsize=0.4
!  call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true. )

! problem definition

  call fill_problem

  call delete(input_probdef)
  call delete(problem)
  call delete ( sol, rhsd )

  call fill_problem


! create system matrix

  call solve

  call delete ( sysmatrix )

  call solve

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  print *, sum(abs(velocity%u))/velocity%n
  print *, sum(abs(pressure%u))/pressure%n

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

contains

  subroutine fill_coefficients

    call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

    coefficients%i(1:11) = &
      [ uintpl,   pintpl,     0,     0,         0,  &
         physqvel, physqpress, 0,     0,     gauss,  &
         gauss ]
    coefficients%i(12:) = 0

    coefficients%r(1) = eta
    coefficients%r(2:) = 0

  end subroutine fill_coefficients

  subroutine fill_mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

! two fixed circular objects

  rp = [0.1_dp,0.1_dp] ! radii of the two circular objects
  xp(1,:) = [0.5_dp,0.5_dp]  ! centre of circular 1
  xp(2,:) = [0.5_dp,0.82_dp] ! centre of circular 2

  call add_to_mesh ( mesh, object='coordinates', nnodes=20, &
    objectsub=objectscoor, objectcoornr=1 )
  call add_to_mesh ( mesh, object='coordinates', nnodes=20, &
    objectsub=objectscoor, objectcoornr=2 )

! write mesh (read by streamfunction computation)

!  call write_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )

  end subroutine fill_mesh

  subroutine fill_problem

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                  1,0,1,0,1,0,1,0,0,    &  ! pressure
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                   [9,3] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

! define constraints on the two objects

  call define_constraint ( mesh, input_probdef, object=1, physq=1, &
    discretization='collocation', nodedof=2 )
  call define_constraint ( mesh, input_probdef, object=2, physq=1, &
    discretization='collocation', nodedof=2 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curves=[1,2,4], physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=1, vvalue=[1._dp,0._dp] )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

  end subroutine fill_problem

  subroutine solve

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    elemsub=elementc, addmatvec=.true., coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  end subroutine solve

end program stokes6
