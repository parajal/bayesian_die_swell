! Poisson problem on a rectangle with Dirichlet boundary conditions
! Weak coupling of two domains using mortars.

module subs_m

  use tfem_elem_m
  implicit none

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

    integer, parameter :: ndf = 9, nodalpb = 3, ndfb = 3, ndflb = 2

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

    call shape_line_P1 ( mesh%objects(object)%xig(node:node,1), psi )

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

      do i = 1, ndflb
        do j = 1, ndf
          elemmat(i,j) = psi(1,i) * phi1(1,j) * &
                                   curvel(1)*mesh%objects(object)%wg(node)
          elemmat2(i,j) = - psi(1,i) * phi2(1,j) * &
                                   curvel(1)*mesh%objects(object)%wg(node)
        end do
      end do

    end if

  end subroutine elementc

end module subs_m

program poisson8

  use tfem_m
  use hsl_ma57_m
  use poisson_elements_m
  use poisson_functions_m
  use figplot_m
  use subs_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3 integration of quads
    gaussb = 3,         & ! 3 point integration of boundary elements
    nx1=5,              & ! number of elements in x (left square)
    ny1=5,              & ! number of elements in y (left square)
    nx2=15,             & ! number of elements in x (right square)
    ny2=15,             & ! number of elements in y (right square)
    lc=1,               & ! define object with left (1) or right (2) domain
    funcnr=4              ! function number for the right-hand side

  real(dp), parameter :: &
    alpha = 1._dp     ! diffusion coefficient


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh1, mesh2
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd, solexact
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients
  type(subscript_t) :: sols


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0

  coefficients%func => func


! create mesh

  meshgen_options%ox = 0.2_dp  ! shift to the right
  meshgen_options%elshape = 6
  meshgen_options%nx = nx1
  meshgen_options%ny = ny1

  call quadrilateral2d ( mesh1, meshgen_options )

  meshgen_options%ox = 1.2_dp  ! shift to the right
  meshgen_options%elshape = 6
  meshgen_options%nx = nx2
  meshgen_options%ny = ny2

  call quadrilateral2d ( mesh2, meshgen_options )

  call mesh_merge ( mesh1, mesh2, mesh, nogroupmerge=.true. )

  call delete ( mesh1, mesh2 )

! create object at the interface between the two domains

  if ( lc == 1 ) then
!   use curve on left mesh
    call add_to_mesh ( mesh, object='curve', objectcurve=2, typeofobject=2, &
      excludegroups=[2], excludegroups2=[1], topology=.true., &
      intrule=2, nsubint=3 )
  else if ( lc == 2 ) then
!   use curve on right mesh
    call add_to_mesh ( mesh, object='curve', objectcurve=8, typeofobject=2, &
      excludegroups=[2], excludegroups2=[1], topology=.true., &
      intrule=2, nsubint=1 )
  end if

  call fill_mesh_parts ( mesh )

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

  plot_options%objectpointsize=0.4
  plot_options%objectpointcolor=4
  call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true. )

! problem definition

  call create_input_probdef ( mesh, input_probdef )

  input_probdef%elementdof(1)%a = 1
  input_probdef%elementdof(2)%a = 1

  call define_essential ( mesh, input_probdef, curves=[1,3,4,5,6,7] )

  call define_constraint ( mesh, input_probdef, object=1, &
    discretization='weak', elementdof=[1,0,1] )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

  call create_subscript ( mesh, problem, sols, degfd=1 )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=7, func=func, funcnr=3 )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=poisson_elem, &
    coefficients=coefficients )

 call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    elemsub=elementc, addmatvec=.true. )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! create exact solution

  call create_sysvector ( problem, solexact )

  call fill_sysvector ( mesh, problem, solexact, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=3 )

! print maximum difference of sol-solexact to standard output

  print *, maxval( abs(sol%u(sols%s) -solexact%u(sols%s)) )

! plot solution

  call plot_color_contour ( plot_options, mesh, problem, filename='sol.fig', &
    sysvector=sol )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, solexact, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program poisson8
