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

program streamfunction1

  use tfem_m
  use streamfunction_elements_m
  use hsl_ma57_m
  use figplot_m
  use io_utils_m

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

  call read_mesh ( mesh, filename='mesh1.out' )

  call add_to_mesh ( mesh, curve=[12,10,7,4,13] )  ! full boundary

  crvboun = mesh%ncurves

  call fill_mesh_parts ( mesh )


! read coefficients

  call read_coefficients ( coefficients, filename='coefficients.out' )


! problem definition stream function

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a = 2

  call define_essential ( mesh, input_probdef, point=1 ) ! Dirichlet in P1

  call problem_definition ( input_probdef, mesh, problem )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )


! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )


! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

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

  call add_boundary_elements ( mesh, problem, rhsd, curve=crvboun, &
    elemsub=streamfunction_natboun_curve, coefficients=coefficients, &
    oldvectors=oldvectors )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%integer_storage=1.4

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

  call plot_color_contour ( plot_options, mesh, problem, 'streamfunction1.fig',&
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

end program streamfunction1
