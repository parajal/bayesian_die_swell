! Steady scalar convection-diffusion problem on a unit square with
! Dirichlet boundary conditions and natural boundary conditions.
! Varying alpha tensor coefficient using a function.
! Varying velocity coefficient using a function.

program convection_diffusion2

  use tfem_m
  use hsl_ma41_m
  use convection_diffusion_elements_m
  use poisson_functions_m
  use convection_diffusion_functions_m
  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3 integration of quads
    gaussb = 3,         & ! 3 point integration of boundary elements
    nx=100,             & ! number of elements in x
    ny=100,             & ! number of elements in y
    funcnr=8,           & ! function number for the right-hand side
    tfuncnr_alpha=1,    & ! function number for the alpha coefficient
    coeff_alpha=1,      & ! how to compute alpha coefficient:
                          ! 0: constant
                          ! 1: given by a function
                          ! 2: given by nodal values in a vector
    vfuncnr_velo=1,     & ! function number for the velocity coefficient
    coeff_velo=1,       & ! how to compute velocity coefficient:
                          ! 0: constant
                          ! 1: given by a function
                          ! 2: given by nodal values in a vector
    vfuncnr=1             ! function number for the flux vector

  real(dp), parameter :: &
    gamma = 1._dp     ! coefficient of the time derivative "rho c_p"

  real(dp), parameter :: &
    velocity(2)=[100._dp,50._dp], & ! constant velocity for coeff_velo=0
!   constant alpha for coeff_alpha=0 (column wise storage)
    alphaten(4)=[1._dp,-0.1_dp,-0.1_dp,2._dp]

! variables

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(vector_t), target :: alpha_nodes
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(2) = coeff_alpha
  coefficients%i(3) = tfuncnr_alpha
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]
  coefficients%i(16) = vfuncnr
  coefficients%i(24) = coeff_velo
  coefficients%i(25) = vfuncnr_velo

  coefficients%r = 0
  coefficients%r(4:5) = velocity
  coefficients%r(8) = gamma
  coefficients%r(9:12) = alphaten

  coefficients%func => func
  coefficients%vfunc => vfunc

  coefficients%vfunc1(1)%p => velocityvfunc   ! vector function for velocity
  coefficients%tfunc1(1)%p => alphatfunc      ! tensor function for alpha


! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a(:,1) = 4  ! for alpha in the nodes

  call define_essential ( mesh, input_probdef, curve1=1 )
  call define_essential ( mesh, input_probdef, curve1=3, curve2=4 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, func=func, funcnr=7 )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, curve2=4, func=func, funcnr=7 )

! oldvectors

  call create ( oldvectors, nvec=1 )

  if ( coeff_alpha == 2 ) then

!   fill alpha in a vector

    write(*,'(a/)') ' alpha in the nodes not yet available.'
    stop

  end if

  if ( coeff_velo == 2 ) then

!   fill velocity in a vector

    write(*,'(a/)') ' Velocity in the nodes not yet available.'
    stop

  end if

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! build diffusion matrix and right-hand side

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_diffusion2_elem, coefficients=coefficients, &
    oldvectors=oldvectors )

  call add_boundary_elements ( mesh, problem, rhsd, curve=2, &
    elemsub=poisson_natboun, coefficients=coefficients )

! convection term

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_ugradc_elem, coefficients=coefficients, &
    oldvectors=oldvectors, addmatvec=.true., buildvector=.false. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma41 ( sysmatrix, rhsd, sol )

  print *, ' average sol = ', sum( abs(sol%u) ) / sol%n

! plot solution

  call plot_color_contour ( plot_options, mesh, problem, filename='sol.fig', &
    sysvector=sol )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol )
  call delete ( rhsd )
  if ( coeff_alpha == 2 ) call delete ( alpha_nodes )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program convection_diffusion2

