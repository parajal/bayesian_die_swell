! Stokes problem on a unit square with Dirichlet boundary conditions.
! Lid-driven cavity flow.
! Connection of two domains (horizontal velocity, slip layer) using curves.

module subs_m

  use tfem_elem_m
  implicit none

contains

! constraint element for vertical velocity

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

!   connection through collocation

    elemmat(1,1)  = 0
    elemmat(1,2)  = 1
    elemmat2 = - elemmat
    elemvec  = 0

  end subroutine elementc

! connection element for horizontal velocity

  subroutine elementcn ( mesh, problem, conn, elem, node, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat11, elemmat12, &
    elemmat21, elemmat22, elemvec1, elemvec2 )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: conn, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat11, elemmat12, &
      elemmat21, elemmat22
    real(dp), intent(out), dimension(:) :: elemvec1, elemvec2

    integer, parameter :: ndf = 3, nodalp = 3, ninti = 3

    integer :: i, j, curve
    real(dp) :: phi(ninti,ndf), xr(ninti,1), x(nodalp,2), wg(ninti)
    real(dp) :: dphi(ninti,ndf,1)
    real(dp) :: dxdxis(ninti,2,1), curvel(ninti)
    real(dp) :: etaslip

!   first curve for connection

    curve = problem%connections(conn)%geometry1

!   Gauss points

    call Gauss_Legendre_line ( ninti, xr(:,1), wg )

!   shape function of the quadratic curve in the integration points

    call shape_line_P2 ( xr, phi, dphi )

!   compute geometry of element deformed

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )

    call isoparametric_deformation_curved ( x, dphi, dxdxis, curvel )

    if ( vector ) then

      elemvec1 = 0
      elemvec2 = 0

    end if

    if ( matrix ) then

      etaslip = coefficients%r(501)

      elemmat11 = 0
      elemmat12 = 0
      elemmat21 = 0
      elemmat22 = 0

!     first velocity component only

      do i = 1, ndf
        do j = 1, ndf
          elemmat11(i,j) = etaslip * sum( phi(:,i) * phi(:,j) * curvel * wg )
          elemmat12(i,j) = - etaslip * sum( phi(:,i) * phi(:,j) * curvel * wg )
          elemmat22(i,j) = etaslip * sum( phi(:,i) * phi(:,j) * curvel * wg )
        end do
      end do

      elemmat21(1:ndf,1:ndf) = transpose(elemmat12(1:ndf,1:ndf))

    end if

  end subroutine elementcn

end module subs_m

program stokes40

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
!  use io_utils_m
!  use figplot_m
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
    nx=20,              & ! number of elements in x
    ny1=16,             & ! number of elements in y (lower domain)
    ny2=4                 ! number of elements in y (upper domain)

  real(dp), parameter :: &
    etaslip = 1.e0_dp, & ! slip viscosity: tau = etaslip * ( u1 - u2 )
    eta = 1._dp          ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh1, mesh2
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
!  type(vector_t) :: velocity, pressure, vorticity
!  type(plot_options_t) :: plot_options
!  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma57_t) :: solver_options


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=501 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0
  coefficients%r(501) = etaslip

!  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%ly = 0.8_dp
  meshgen_options%nx = nx
  meshgen_options%ny = ny1

  call quadrilateral2d ( mesh1, meshgen_options )

  meshgen_options%elshape = 6
  meshgen_options%ly = 0.2_dp
  meshgen_options%oy = 0.8_dp
  meshgen_options%nx = nx
  meshgen_options%ny = ny2

  call quadrilateral2d ( mesh2, meshgen_options )

  call mesh_merge ( mesh1, mesh2, mesh, nogroupmerge=.true. )

  call add_to_mesh ( mesh, curve=[-3] )     ! curve 9

  call delete ( mesh1, mesh2 )

! write mesh (read by streamfunction computation)

!  call write_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )

!  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
!  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

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

  call define_constraint ( mesh, input_probdef, curve1=9, curve2=5, &
    discretization='collocation', physq=1, nodedof=1, exclude=3 )

  call define_connection ( mesh, input_probdef, curve1=9, curve2=5, &
    discretization='weak', physq=1 )

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
  call create_sysmatrix_structure_connection ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

 call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    elemsub=elementc, addmatvec=.true. )

 call build_system_connection ( mesh, problem, sysmatrix, rhsd, &
    elemsub=elementcn, addmatvec=.true., coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%real_storage = 2.0

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

  print *, sum( sol%u ) / sol%n

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program stokes40
