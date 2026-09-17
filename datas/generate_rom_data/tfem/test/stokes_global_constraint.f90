! Stokes problem on a unit square. Periodical boundary conditions
! weak coupling
! Flow rate imposed using a global constraint on an object
! The position of the constraint can be located at any position
! Flow rate computation in any cross-section using integrate_object

module subs_m

  use tfem_elem_m

  implicit none

  integer :: ndf = 9, ndim = 2
  integer :: ndfb = 3, nodalpb = 3
  real(dp) :: xp(2) = 0

contains


! Element for the flowrate constraint on an object

  subroutine elementc ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


    integer :: object, i, j
    real(dp) :: flowrate

    real(dp) :: phi(1,ndf), xr(1,ndim), phib(1,ndfb), x(nodalpb,ndim)
    real(dp) :: dphib(1,ndfb,1), dxdxi(1,ndim), curvel(1), normal(1,ndim)
    real(dp) :: tmp(ndf,ndim)


!   imposed flow rate

    object = problem%constraints(constr)%object

!   reference coordinates of the integration point (node) of the element (elem)

    xr(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)

!   shape function of the velocity in the integration point

    call shape_quad_Q2 ( xr, phi )

!   shape function of the quadratic curve in the integration point

    call shape_line_P2 ( mesh%objects(object)%xig(node:node,1), phib, &
      dphib(:,:,1) )

!   compute deformed element

    call get_coordinates_object ( mesh, elem, x, object )

    call isoparametric_deformation_curve ( x, dphib(:,:,1), dxdxi, curvel, &
      normal )

    if ( vector ) then

!     specify flowrate

      if ( first ) then
!       first integration point: rhs's of equations
!       specify flowrate
        flowrate = coefficients%r(6)
        elemvec(1) = flowrate
      else
        elemvec(1) = 0
      end if

    end if

    if ( matrix ) then

      do i = 1, ndf
        do j = 1, ndim
          tmp(i,j) = &
             phi(1,i) * curvel(1) * mesh%objects(object)%wg(node) * normal(1,j)
        end do
      end do

      elemmat(1,:) = reshape ( tmp, [ndf*ndim] )

    end if

  end subroutine elementc


! Element for computing the flowrate on an object

  subroutine element_flowrate ( mesh, problem, ei, first, last, coefficients, &
    oldvectors, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(eleminfo_t), intent(in) :: ei
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec


    integer, parameter :: physqvel = 1
    integer :: elgrp, elem, object, elemo, nintps, ip

    real(dp) :: phi(ei%nintps,ndf), xr(ei%nintps,ndim)
    real(dp) :: phib(ei%nintps,ndfb), x(nodalpb,ndim)
    real(dp) :: dphib(ei%nintps,ndfb,1), dxdxi(ei%nintps,ndim)
    real(dp) :: curvel(ei%nintps), normal(ei%nintps,ndim)
    real(dp) :: tmp(ndf,ndim), uvector(ei%nintps,ndim)
    real(dp) :: un(ei%nintps), u(ndf*ndim)

    elgrp  = ei%elgrp
    elem   = ei%elem
    object = ei%object
    nintps = ei%nintps
    elemo  = ei%elemo

!   reference coordinates in the element of the mesh

    xr = mesh%objects(object)%refcoor_int(ei%intps,:,elemo)

!   shape function of the velocity

    call shape_quad_Q2 ( xr, phi )

!   reference coordinates in the element of the object

    xr(:,1) = mesh%objects(object)%xig(ei%intps,1)

!   shape function of the quadratic curve in the integration points

    call shape_line_P2 ( xr(:,1), phib, dphib(:,:,1) )

!   compute deformed element

    call get_coordinates_object ( mesh, elemo, x, object )

    call isoparametric_deformation_curve ( x, dphib(:,:,1), dxdxi, curvel, &
      normal )

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel] )

    tmp = reshape ( u, [ndf,ndim] )

    uvector = matmul ( phi, tmp )

!   normal velocity

    do ip = 1, nintps
      un(ip) = dot_product ( uvector(ip,:), normal(ip,:) )
    end do

    elemvec(1) = sum ( un * curvel * mesh%objects(object)%wg(ei%intps) )

  end subroutine element_flowrate

end module subs_m

program stokes18

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use io_utils_m
  use subs_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
!    pintpl = 2,         & ! P1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=10,              & ! number of elements in x
    ny=5                  ! number of elements in y

  real(dp), parameter :: &
    shiftx = 0.48_dp,  & ! shift object to the left by shiftx from initial x=1
    eta = 1._dp,  & ! viscosity
    flowrate = 2._dp

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

  integer :: ip
  real(dp) :: resultsum(1)


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gauss ]
  coefficients%i(12:) = 0

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(6) = flowrate

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

  call add_to_mesh ( mesh, object='curve', objectcurve=2, &
    topology=.true., intrule=2, nsubint=1 )

  mesh%objects(1)%coor(:,1) = mesh%objects(1)%coor(:,1) - shiftx

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 2 ! velocity
  if ( pintpl == 2 ) then
!   P1 discontinuous pressure
    input_probdef%vec_elementdof(1)%a(:,2) = [ 0,0,0,0,0,0,0,0,3 ]
  else if ( pintpl == 4 ) then
!   Q1 continuous pressure
    input_probdef%vec_elementdof(1)%a(:,2) = [ 1,0,1,0,1,0,1,0,0 ]
  end if
  input_probdef%vec_elementdof(1)%a(:,3) = 1 ! scalar, such as vorticity

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=3, physq=1 )
  if ( pintpl == 2 ) then
!   P1 discontinuous pressure
    call define_essential ( mesh, input_probdef, element=1, elnode=9, &
      degfd=[1], physq=2 )
  else if ( pintpl == 4 ) then
!   Q1 continuous pressure
    call define_essential ( mesh, input_probdef, point=1, physq=2 )
  end if

  call define_constraint ( mesh, input_probdef, &
    physq=1, object=1, nglobalc=1, discretization='weak' )

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, curve2=5, discretization='weak', elementdof=[2,0,2] )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=1, value=0._dp )
  if ( pintpl == 2 ) then
!   P1 discontinuous pressure
    call fill_sysvector ( mesh, problem, sol, &
      element=1, elnode=9, degfd=1, physq=2, value=0._dp )
  else if ( pintpl == 4 ) then
!   Q1 continuous pressure
    call fill_sysvector ( mesh, problem, sol, &
      point=1, physq=2, value=0._dp )
  end if

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
    constraint1=1, elemsub=elementc, addmatvec=.true., &
     coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=2, elemsub=stokes_constr_elem_conn, addmatvec=.true., &
     coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  ip = problem%constraints(1)%globnumdegfd(1)

  print *, 'pressure gradient =', sol%u( problem%degfdperm(ip+1,2) )

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

! check flow rate

  call integrate_object ( mesh, problem, resultsum, elemsub1=element_flowrate,&
    object=1, oldvectors=oldvectors )

  print *, 'flowrate in object = ', resultsum

! check flow rate

  call integrate_object ( mesh, problem, resultsum, elemsub1=element_flowrate,&
    object=1, oldvectors=oldvectors )

  mesh%objects(1)%coor(:,1) = 0.45_dp

  call find_refcoor_objects ( mesh, object1=1 )

  call integrate_object ( mesh, problem, resultsum, elemsub1=element_flowrate,&
    object=1, oldvectors=oldvectors )

  print *, 'flowrate in object (shifted) ', resultsum

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program stokes18
