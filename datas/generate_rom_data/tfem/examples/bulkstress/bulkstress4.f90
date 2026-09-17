! Stokes problem on a square domain with a disk at the center.
! Shear flow boundary conditions (set in "stokes_functions_m").
! Freely rotating conditions on the disk boundaries.
! Weak constraint on the boundary.
! Integration of fields over the domain and boundaries.

module subs_m

  use tfem_elem_m
  use math_defs_m

  implicit none

  integer :: ndf = 9, ndim = 2
  integer :: ndfb = 3, nodalpb = 3
  integer :: ndflb = 2  ! number of degrees of freedom of the Lag. Mult.
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

    integer :: object, i, j
    real(dp) :: phi(1,ndf), xr(1,2), r(2), psi(1,ndflb), x(nodalpb,2)
    real(dp) :: theta(1,ndfb), dtheta(1,ndfb,1), dxdxi(1,2), curvel(1)

    object = problem%constraints(constr)%object

    xr(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)
    r = mesh%objects(object)%coor_int(node,:,elem) - xp(:)

    call shape_quad_Q2 ( xr, phi )

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
        elemmatadd(i,:) = psi(1,i) * [ -1._dp,  0._dp,  -r(2) ] * &
                                   curvel(1) * mesh%objects(object)%wg(node)
        elemmatadd(i+ndflb,:) = psi(1,i) * [  0._dp, -1._dp, r(1) ] * &
                                   curvel(1) * mesh%objects(object)%wg(node)
      end do

    end if

  end subroutine elementc

! Element for computing the Lagrange multipliers integration on an object

  subroutine lm_integration ( mesh, problem, ei, first, last, coefficients, &
    oldvectors, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(eleminfo_t), intent(in) :: ei
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: object, elemo, nintps, ip, constrnr

    real(dp) :: phib(ei%nintps,ndfb), x(nodalpb,ndim)
    real(dp) :: dphib(ei%nintps,ndfb,1), dxdxi(ei%nintps,ndim)
    real(dp) :: curvel(ei%nintps), normal(ei%nintps,ndim)
    real(dp) :: tmp(ndflb,ndim), lmvector(ei%nintps,ndim), xg(ei%nintps,ndim)
    real(dp) :: u(ndflb*ndim), psi(ei%nintps,ndflb),tauten(ei%nintps,ndim,ndim)

    object = ei%object
    nintps = ei%nintps
    elemo  = ei%elemo

!   number of constraint where Lagrange multipliers are imposed

    constrnr = 1

!   shape function of the Langrangian multiplier in the integration point

    call shape_line_P1 ( mesh%objects(object)%xig(ei%intps,1), psi )

!   shape function of the quadratic curve in the integration points

    call shape_line_P2 ( mesh%objects(object)%xig(ei%intps,1), phib, &
      dphib(:,:,1) )

!   compute deformed element

    call get_coordinates_object ( mesh, elemo, x, object )

    call isoparametric_deformation_curve ( x, dphib(:,:,1), dxdxi, curvel, &
      normal )

    call isoparametric_coordinates ( x, phib, xg )

!   get Lagrange multiplier values

    call get_sysvector_constraint ( mesh, problem, oldvectors%s(1)%p, &
      constraint=constrnr, elem=elemo, u=u )

    tmp = reshape ( u, [ndflb,ndim] )

    lmvector = matmul ( psi, tmp )

    do ip = 1, nintps
      tauten(ip,1,:) = lmvector(ip,1) * xg(ip,:)
      tauten(ip,2,:) = lmvector(ip,2) * xg(ip,:)
    end do

    elemvec(1) = sum ( tauten(:,1,1) * curvel * &
                                mesh%objects(object)%wg(ei%intps) )
    elemvec(2) = sum ( tauten(:,1,2) * curvel * &
                                mesh%objects(object)%wg(ei%intps) )
    elemvec(3) = sum ( tauten(:,2,1) * curvel * &
                                mesh%objects(object)%wg(ei%intps) )
    elemvec(4) = sum ( tauten(:,2,2) * curvel * &
                                mesh%objects(object)%wg(ei%intps) )

  end subroutine lm_integration

end module subs_m

program bulkstress4

  use tfem_m
  use stokes_elements_m
  use figplot_m
  use subs_m
  use hsl_ma57_m
  use stokes_functions_m

  implicit none

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
  type(mesh_t) :: mesh, mesh_particle
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

  integer :: ip, elem
  real(dp) :: up, vp, omega, total_stress_domain(3), lm_integral(4)
  real(dp) :: ox, oy, lx, ly

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

! create mesh for particle boundary

  call mesh_skeleton ( mesh_particle, nnodes=20, nelem=10, elshape=2, ndim=2 )

! object

  rp = 1.0_dp ! radius of the circular object

! set initial positions

  xp = [ 0.0_dp, 0.0_dp ]

  call objectscoor ( 1, mesh_particle%coor )

  do elem = 1, mesh_particle%nelem
    mesh_particle%topology(1)%a(:,elem) = [ 2*elem-1, 2*elem, 2*elem + 1 ]
  end do
  mesh_particle%topology(1)%a(3,mesh_particle%nelem) = 1 ! close circle

! one object

  call add_to_mesh ( mesh, object='mesh', objectmesh=mesh_particle, &
    topology=.true., intrule=1, nsubint=10 )

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
    discretization='weak', elementdof=[2,0,2], naddunknowns=3 )

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
    constraint1=1, elemsub=elementc, addmatvec=.true. )

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

! integration of Lagrange multipliers over the object

  call integrate_object ( mesh, problem, lm_integral, elemsub1=lm_integration, &
    object=1, coefficients=coefficients, oldvectors=oldvectors )

  write(*,*) ' xx-component of bulk stress: ', &
    (total_stress_domain(1) + lm_integral(1))/400.0_dp
  write(*,*) ' xy-component of bulk stress: ', &
    (total_stress_domain(2) + lm_integral(2))/400.0_dp
  write(*,*) ' yx-component of bulk stress: ', &
    (total_stress_domain(2) + lm_integral(3))/400.0_dp
  write(*,*) ' yy-component of bulk stress: ', &
    (total_stress_domain(3) + lm_integral(4))/400.0_dp

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( coefficients )
  call delete ( oldvectors )

end program bulkstress4
