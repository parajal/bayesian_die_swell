! Stokes problem on a square domain with a disk at the center.
! Shear flow boundary conditions (set in "stokes_functions_m").
! Freely rotating conditions on the disk boundaries.
! Constraint on the boundary using collocation.
! Integration of fields over the domain and boundaries.

module subs_m

  use tfem_elem_m
  use math_defs_m

  implicit none

  integer :: ndim = 2
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
    integer, parameter :: ndf = 9
    real(dp) :: phi(1,ndf), xr(1,2), r(2)

    object = problem%constraints(constr)%object

    xr(1,:) = mesh%objects(object)%refcoor(node,:)
    r = mesh%objects(object)%coor(node,:) - xp

    call shape_quad_Q2 ( xr, phi )

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

    if ( matrix ) then

      elemmat(1,1:ndf) = phi(1,:)
      elemmat(1,ndf+1:) = 0
      elemmat(2,1:ndf) = 0
      elemmat(2,ndf+1:) = phi(1,:)
      elemmatadd(1,:) = [ -1._dp,  0._dp, -r(2) ]
      elemmatadd(2,:) = [  0._dp, -1._dp,  r(1) ]

    end if

  end subroutine elementc

! Element for computing the traction forces on a particle by Lagrange
! multipliers

  subroutine lm_traction ( mesh, problem, object, oldvectors, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemvec
    integer, intent(in) :: object

    integer :: i

    real(dp) :: u(ndim)

    elemvec(:,:) = 0.0

    do i = 1,size(mesh%objects(object)%coor(:,1))

      call get_sysvector_constraint ( mesh, problem, oldvectors%s(1)%p, &
        object, node=i, u=u )

      if ( ndim == 2 ) then

        elemvec(1,1) = elemvec(1,1) + mesh%objects(object)%coor(i,1)*u(1)
        elemvec(1,2) = elemvec(1,2) + mesh%objects(object)%coor(i,1)*u(2)
        elemvec(2,1) = elemvec(2,1) + mesh%objects(object)%coor(i,2)*u(1)
        elemvec(2,2) = elemvec(2,2) + mesh%objects(object)%coor(i,2)*u(2)

      else if ( ndim == 3 ) then

        elemvec(1,1) = elemvec(1,1) + mesh%objects(object)%coor(i,1)*u(1)
        elemvec(1,2) = elemvec(1,2) + mesh%objects(object)%coor(i,1)*u(2)
        elemvec(1,3) = elemvec(1,3) + mesh%objects(object)%coor(i,1)*u(3)
        elemvec(2,1) = elemvec(2,1) + mesh%objects(object)%coor(i,2)*u(1)
        elemvec(2,2) = elemvec(2,2) + mesh%objects(object)%coor(i,2)*u(2)
        elemvec(2,3) = elemvec(2,3) + mesh%objects(object)%coor(i,2)*u(3)
        elemvec(3,1) = elemvec(3,1) + mesh%objects(object)%coor(i,3)*u(1)
        elemvec(3,2) = elemvec(3,2) + mesh%objects(object)%coor(i,3)*u(2)
        elemvec(3,3) = elemvec(3,3) + mesh%objects(object)%coor(i,3)*u(3)

      end if

    end do

  end subroutine lm_traction

end module subs_m

program bulkstress3

  use tfem_m
  use stokes_elements_m
  use figplot_m
  use subs_m
  use hsl_ma57_m
  use stokes_functions_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    gaussb = 3,         & ! 3 point integration of boundary elements
    nx = 100,           & ! number of elements in x
    ny = 100              ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp           ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

  integer :: ip, objectnr
  real(dp) :: up, vp, omega, total_stress_domain(3), lm_integral(2,2)
  real(dp) :: area_domain(1), ox, oy, lx, ly

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gaussb ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

! create mesh

  ox = -10.0_dp
  oy = -10.0_dp
  lx =  20.0_dp
  ly =  20.0_dp

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%ox = ox
  meshgen_options%oy = oy
  meshgen_options%lx = lx
  meshgen_options%ly = ly

  call quadrilateral2d ( mesh, meshgen_options )

! object

  rp = 1.0 ! radius of the circular object

! set initial positions

  xp = [ 0.0_dp, 0.0_dp ]

  call add_to_mesh ( mesh, object='coordinates', nnodes=28, &
    objectsub=objectscoor, objectcoornr=1 )

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

! define essential boundaries

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, &
    physq=physqvel )
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! define constraints on the objects

  call define_constraint ( mesh, input_probdef, object=1, physq=physqvel, &
    discretization='collocation', nodedof=2, naddunknowns=3 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, curve1=1, curve2=4, &
    physq=physqvel, degfd=1, func=func, funcnr=1 )
  call fill_sysvector ( mesh, problem, sol, curve1=1, curve2=4, &
    physq=physqvel, degfd=2, func=func, funcnr=2 )

  call fill_sysvector ( mesh, problem, sol, point=1, physq=physqpress, &
    value=0._dp )

! create system matrix

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

! post-processing of the velocity and rotation rate of the particle

  ip = problem%constraints(1)%addnumdegfd(1)

  up = sol%u( problem%degfdperm(ip+1,2) )
  vp = sol%u( problem%degfdperm(ip+2,2) )
  omega = sol%u( problem%degfdperm(ip+3,2) )

  write (*,*) ' up = ', up
  write (*,*) ' vp = ', vp
  write (*,*) ' omega = ', omega

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

! total stress integration over the domain

  call integrate ( mesh, problem, total_stress_domain, &
    elemsub=stokes_integrate_total_stress, &
    coefficients=coefficients, oldvectors=oldvectors)

! traction on particle boundary through Lagrange multipliers

  objectnr = 1

  call lm_traction ( mesh, problem, objectnr, oldvectors, lm_integral )

! area of the domain

  call integrate ( mesh, problem, area_domain, &
    elemsub=stokes_integrate_volume, &
    coefficients=coefficients, oldvectors=oldvectors)

  write(*,*) ' xx-component of bulk stress: ', &
    (total_stress_domain(1) + lm_integral(1,1))/area_domain
  write(*,*) ' xy-component of bulk stress: ', &
    (total_stress_domain(2) + lm_integral(1,2))/area_domain
  write(*,*) ' yx-component of bulk stress: ', &
    (total_stress_domain(2) + lm_integral(2,1))/area_domain
  write(*,*) ' yy-component of bulk stress: ', &
    (total_stress_domain(3) + lm_integral(2,2))/area_domain

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( coefficients )
  call delete ( oldvectors )

end program bulkstress3
