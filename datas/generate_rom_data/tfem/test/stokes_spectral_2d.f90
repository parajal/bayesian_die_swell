! Stokes problem on a unit square with Dirichlet boundary conditions.
! Lid-driven cavity flow.
! Physical quantities.
! Spectral elements.
! This problem is similar to stokes1, but now with spectral elements.

program stokes34

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use io_utils_m
!  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 13,    & ! Qp velocities
    pintpl = 14,    & ! Pp-1 pressures
    p = 8,          & ! polynomial order
    q = p-1,        & ! polynomial order pressure
    physqvel = 1,   & ! physical quantity nr of the velocities
    physqpress = 2, & ! physical quantity nr of the pressures
    gauss = p+1,    & ! number of Gauss points
    gaussb = p+1,   & ! number of Gauss points on the boundary
    inttype = 1,    & ! Gauss-Legendre-Lobatto for integration
    nx=2,           & ! number of elements in x
    ny=2              ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp     ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh!, mesh_plot
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
!  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gaussb ]
  coefficients%i(12:) = 0
  coefficients%i(40) = inttype
  coefficients%i(51) = p ! velocity polynomial order
  coefficients%i(52) = q ! pressure polynomial order

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

!  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  meshgen_options%elshape = 102
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%p = p

  call quadrilateral2d ( mesh, meshgen_options )

! write mesh (read by streamfunction computation)

!  call write_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )

! Create plot mesh

!  call mesh_convert ( mesh, mesh_plot, elementshapes='spectraltolinear', &
!    warn=.false. )

! write mesh_plot (read by streamfunction computation)

!  call write_mesh ( mesh_plot, filename='mesh_plot.out' )

!  call fill_mesh_parts( mesh_plot )

!  call plot_points_curves ( plot_options, mesh_plot, filename='curves.fig' )
!  call plot_mesh ( plot_options, mesh_plot, filename='mesh_plot.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 2  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0  ! put pressure degrees in one
  input_probdef%vec_elementdof(1)%a(p+3,2) = (q+1)*(q+2)/2 ! internal point
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )
  call define_essential ( mesh, input_probdef, element=1, elnode=p+3, &
    physq=2, degfd=[1] )  ! Set one pressure degree to zero

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=1, degfd=1, value=1._dp )
  call fill_sysvector ( mesh, problem, sol, &
    element=1, elnode=p+3, physq=2, degfd=1, value=0._dp )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

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

  print *, sum(velocity%u)/velocity%n
  print *, sum(pressure%u)/pressure%n
  print *, sum(vorticity%u)/vorticity%n

! write to a fig file for plotting

!  call plot_vector ( plot_options, mesh_plot, problem, 'velocity.fig', &
!    vector=velocity )

!  call plot_color_fill ( plot_options, mesh_plot, problem, &
!    'pressure_color.fig', vector=pressure )

!  call plot_color_contour ( plot_options, mesh_plot, problem, &
!    'pressure_contour.fig', vector=pressure )

!  call plot_color_fill ( plot_options, mesh_plot, problem, &
!    'vorticity_color.fig', vector=vorticity )

!  call plot_color_contour ( plot_options, mesh_plot, problem, &
!    'vorticity_contour.fig', vector=vorticity )

! write binary file for reading by streamfunction

!  open(unit=10,form='unformatted',file='velocity_bin.out')

!  write(10) velocity%u

!  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program stokes34
