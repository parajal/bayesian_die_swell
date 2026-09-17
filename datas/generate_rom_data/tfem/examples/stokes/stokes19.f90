! Stokes problem on a unit square with Dirichlet boundary conditions.
! Lid-driven cavity flow.
! Non-conforming coupling of two domains.

module subs_m

  use tfem_elem_m

  implicit none

  integer, parameter :: &
!    lintpl = 0, ndflb = 1    ! P0 Lagrange multiplier
    lintpl = 1, ndflb = 2     ! P1 Lagrange multiplier

contains

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

    integer, parameter :: ndf = 9, nodalpb = 3, ndfb = 3

    integer :: object, i, j
    real(dp) :: phi1(1,ndf), phi2(1,ndf), xr(1,2), psi(1,ndflb), x(nodalpb,2)
    real(dp) :: theta(1,ndfb), dtheta(1,ndfb,1), dxdxi(1,2), curvel(1)

    object = problem%constraints(constr)%object

!   reference coordinates and shape function in first intersection

    xr(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)
    call shape_quad_Q2 ( xr, phi1 )

!   reference coordinates and shape function in second intersection

    xr(1,:) = mesh%objects(object)%refcoor2_int(node,:,elem)
    call shape_quad_Q2 ( xr, phi2 )

!   shape function of the Lagrange multiplier

    if ( lintpl == 0 ) then
!     constant Lagrange multiplier
      psi = 1
    else if ( lintpl == 1 ) then
!     linear Lagrange multiplier
      call shape_line_P1 ( mesh%objects(object)%xig(node:node,1), psi )
    end if

!   shape function of the quadratic curve in the integration point

    call shape_line_P2 ( mesh%objects(object)%xig(node:node,1), theta, &
      dtheta(:,:,1) )

!   compute geometry of element deformed

    call get_coordinates_object ( mesh, elem, x, object )

    call isoparametric_deformation_curve ( x, dtheta(:,:,1), dxdxi, curvel )

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     first velocity component (diagonal block)

      do i = 1, ndflb
        do j = 1, ndf
          elemmat(i,j) = psi(1,i) * phi1(1,j) * &
                                   curvel(1)*mesh%objects(object)%wg(node)
          elemmat2(i,j) = - psi(1,i) * phi2(1,j) * &
                                   curvel(1)*mesh%objects(object)%wg(node)
        end do
      end do

!     second velocity component (diagonal block)

      elemmat(ndflb+1:2*ndflb,ndf+1:2*ndf) = elemmat(1:ndflb,1:ndf)
      elemmat2(ndflb+1:2*ndflb,ndf+1:2*ndf) = elemmat2(1:ndflb,1:ndf)

!     zero off-diagonal blocks

      elemmat(1:ndflb,ndf+1:2*ndf) = 0
      elemmat(ndflb+1:2*ndflb,1:ndf) = 0
      elemmat2(1:ndflb,ndf+1:2*ndf) = 0
      elemmat2(ndflb+1:2*ndflb,1:ndf) = 0

    end if

  end subroutine elementc

end module subs_m

program stokes19

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use io_utils_m
  use figplot_m
  use subs_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
!    pintpl = 4,         & ! Q1 pressures
    pintpl = 2,         & ! P1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    lc=1,               & ! define object with lower (1) or upper (2) domain
    nx1=20,             & ! number of elements in x (lower domain)
    ny1=16,             & ! number of elements in y (lower domain)
    e1=2, e2=3            ! for each e1 elements in the lower domain there are
                          ! e2 elements in the upper domain.

  real(dp), parameter :: &
    eta = 1._dp     ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh1, mesh2
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

  integer :: nx2, ny2, elementdof(3)


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

! test number of elements

  if ( mod(nx1,e1) /= 0 ) then
    write(*,'(/a/)') 'Error: nx1 must be a multiple of e1 '
    stop
  else
!   refinement factor is e2/e1
    nx2 = nx1 * e2 / e1
    ny2 = nx2 / 5       ! height is five times smaller
  end if

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%ly = 0.8_dp
  meshgen_options%nx = nx1
  meshgen_options%ny = ny1

  call quadrilateral2d ( mesh1, meshgen_options )

  meshgen_options%elshape = 6
  meshgen_options%ly = 0.2_dp
  meshgen_options%oy = 0.8_dp
  meshgen_options%nx = nx2
  meshgen_options%ny = ny2

  call quadrilateral2d ( mesh2, meshgen_options )

  call mesh_merge ( mesh1, mesh2, mesh, nogroupmerge=.true. )

! create object at the interface between the two domains

  if ( lc == 1 ) then
!   use curve on lower mesh
    call add_to_mesh ( mesh, object='curve', objectcurve=3, typeofobject=2, &
      excludegroups=[2], excludegroups2=[1], topology=.true., &
      intrule=2, nsubint=e2 )  ! number of intervals on e1 elements is e1*e2
  else if ( lc == 2 ) then
!   use curve on upper mesh
    call add_to_mesh ( mesh, object='curve', objectcurve=5, typeofobject=2, &
      excludegroups=[2], excludegroups2=[1], topology=.true., &
      intrule=2, nsubint=e1 )  ! number of intervals on e2 elements is e1*e2
  end if

  call delete ( mesh1, mesh2 )

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

  plot_options%objectpointsize=0.4
  plot_options%objectpointcolor=4
  call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true. )

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

  input_probdef%vec_elementdof(2)%a =  input_probdef%vec_elementdof(1)%a

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, curve2=2, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=4, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=6, curve2=8, physq=1 )
  if ( pintpl == 2 ) then
!   P1 discontinuous pressure
    call define_essential ( mesh, input_probdef, element=1, elnode=9, &
      degfd=[1], physq=2 )
  else if ( pintpl == 4 ) then
!   Q1 continuous pressure
    call define_essential ( mesh, input_probdef, point=1, physq=2 )
  end if

  if ( lintpl == 0 ) then
!   constant Lagrange multiplier
    elementdof = [0,2,0]
  else if ( lintpl == 1 ) then
!   linear Lagrange multiplier
    elementdof = [2,0,2]
  end if

  call define_constraint ( mesh, input_probdef, object=1, &
    discretization='weak', physq=1, elementdof=elementdof )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0

  call fill_sysvector ( mesh, problem, sol, &
    curve1=7, exclude=3, physq=1, degfd=1, value=1._dp )

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
    elemsub=elementc, addmatvec=.true. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%real_storage = 2.0

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

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

! write to a fig file for plotting

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  call plot_color_fill ( plot_options, mesh, problem, 'pressure_color.fig', &
    vector=pressure )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )

  call plot_color_fill ( plot_options, mesh, problem, 'vorticity_color.fig', &
    vector=vorticity )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

! write binary file for reading by streamfunction

  open(unit=10,form='unformatted',file='velocity_bin.out')

  write(10) velocity%u

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program stokes19
