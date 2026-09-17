! Streamfunction problem on a unit square with Neumann boundary conditions.
! Plot with figplot.

! Internal streamfunction element routine using the Poisson Equation:
!
!   - nabla^2 psi = omega
!

! Boundary element for a natural boundary on a curve for the streamfunction
! equation (Poisson equation):
!
!    dpsidn = -v * nx + u * ny   (= -tangential velocity)
!
! since dpsidx = -v, dpsidy = u.
!

! Dirichlet in a single point (P1) psi=0

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

program streamfunction20

  use tfem_m
  use streamfunction_elements_m
  use hsl_ma57_m
  use figplot_m
  use io_utils_m
  use subs_m

  implicit none


! definitions

  type(mesh_t) :: mesh
  type(plot_options_t) :: plot_options
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(vector_t), target :: velocity
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma57_t) :: solver_options

  integer :: crvboun


! print filename in title and date at footer

  plot_options%printfilename = .true.
  plot_options%printdateandtime = .true.


! read mesh

  call read_mesh ( mesh, filename='mesh.out' )

  call add_to_mesh ( mesh, curve=[1,2,3,4,5,6,7,8] )  ! full boundary

  crvboun = mesh%ncurves

  call fill_mesh_parts ( mesh )


! read coefficients

  call read_coefficients ( coefficients, filename='coefficients.out' )


! problem definition stream function

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a = 2
  input_probdef%elementdof(2)%a = 1
  input_probdef%vec_elementdof(2)%a = 2

  call define_essential ( mesh, input_probdef, point=1 ) ! Dirichlet in P1

  call define_constraint ( mesh, input_probdef, object=1, &
    discretization='weak', elementdof=[1,0,1] )

  call problem_definition ( input_probdef, mesh, problem )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )


! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )


! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! read velocity

  call create_vector ( problem, velocity, vec=1 )

  open ( unit=10, form='unformatted', file='velocity_bin.out' )

  read ( unit=10 ) velocity%u

  close ( unit=10 )


! build (assemble) matrix and vector from elements

  call create_oldvectors ( oldvectors, nvec=1 )
  oldvectors%v(1)%p => velocity

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=streamfunction_elem, coefficients=coefficients, &
    oldvectors=oldvectors )

 call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    elemsub=elementc, addmatvec=.true. )

  call add_boundary_elements ( mesh, problem, rhsd, curve=crvboun, &
    elemsub=streamfunction_natboun_curve, coefficients=coefficients, &
    oldvectors=oldvectors )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%integer_storage=1.4

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

  call plot_color_contour ( plot_options, mesh, problem, 'streamfunction.fig', &
    sysvector=sol )


! delete data

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program streamfunction20
