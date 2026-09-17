! Time-dependent scalar convection-diffusion problem on a unit square with
! Dirichlet boundary conditions and natural boundary conditions.
! Varying alpha coefficient using a function or a vector in the nodes.
! Varying beta coefficient using a function or a vector in the nodes.
! Varying velocity coefficient using a function.
! Second-order implicit Gear discretization. First step: Euler implicit.

program convection_diffusion4

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
    numtimesteps=30,    & ! number of time steps
    funcnr=8,           & ! function number for the right-hand side
    funcnr_alpha=1,     & ! function number for the alpha coefficient
    coeff_alpha=1,      & ! how to compute alpha coefficient:
                          ! 0: constant
                          ! 1: given by a function
                          ! 2: given by nodal values in a vector
    funcnr_beta=1,      & ! function number for the beta coefficient
    coeff_beta=0,       & ! how to compute beta coefficient:
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
    alpha = 1._dp, &  ! diffusion coefficient
    beta  = 200._dp, & ! beta (relaxation,reaction) coefficient
    gamma = 1._dp, &  ! coefficient of the time derivative "rho c_p"
    deltat = 1.e-2_dp ! time step

  real(dp), parameter :: &
    gamma0=1.5_dp, alpha0=2._dp, alpha1=-0.5_dp, &  ! Gear parameters
    velocity(2)=[100._dp,10._dp]  ! constant velocity for coeff_velo=0

! variables

  logical :: buildmatrix

  integer :: step

  real(dp) :: mfac

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, soln, solnm1, rhsd
  type(sysvector_t), target :: solhat
  type(vector_t), target :: alpha_nodes
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors
  type(lu_ma41_t) :: lu


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(2) = coeff_alpha
  coefficients%i(3) = funcnr_alpha
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]
  coefficients%i(16) = vfuncnr
  coefficients%i(24) = coeff_velo
  coefficients%i(25) = vfuncnr_velo
  coefficients%i(26) = coeff_beta
  coefficients%i(27) = funcnr_beta

  coefficients%r = 0
  coefficients%r(1) = alpha
  coefficients%r(3) = deltat
  coefficients%r(4:5) = velocity
  coefficients%r(7) = beta
  coefficients%r(8) = gamma

  coefficients%func => func
  coefficients%vfunc => vfunc

  coefficients%func1(1)%p => alphafunc   ! scalar function for alpha
  coefficients%func1(2)%p => betafunc    ! scalar function for beta
  coefficients%vfunc1(1)%p => velocityvfunc   ! vector function for velocity


! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a(:,1) = 1  ! for alpha in the nodes

  call define_essential ( mesh, input_probdef, curve1=1 )
  call define_essential ( mesh, input_probdef, curve1=3, curve2=4 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, soln, solnm1, solhat )
  call create_sysvector ( problem, rhsd )

! fill initial solution vector with zero

  soln%u = 0

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, func=func, funcnr=7 )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, curve2=4, func=func, funcnr=7 )

! oldvectors

  call create ( oldvectors, nsysvec=1, nvec=1 )

  oldvectors%s(1)%p => solhat

  if ( coeff_alpha == 2 ) then

!   fill alpha in a vector

    call create_vector ( problem, alpha_nodes, vec=1 )

!   as an example use a function to fill alpha in the nodes
    call fill_vector ( mesh, problem, alpha_nodes, &
      node1=1, node2=mesh%nnodes, func=alphafunc, funcnr=funcnr_alpha )

    oldvectors%v(1)%p => alpha_nodes

  end if

  if ( coeff_velo == 2 ) then

!   fill velocity in a vector

    write(*,'(a/)') ' Velocity in the nodes not yet available.'
    stop

  end if

  if ( coeff_beta == 2 ) then

!   fill beta in a vector

    write(*,'(a/)') ' beta in the nodes not yet available.'
    stop

  end if

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! start time stepping

  buildmatrix = .true.

  do step = 1, numtimesteps

    if ( step == 3 ) buildmatrix = .false.  ! matrix is constant

!   build diffusion matrix and right-hand side

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_diffusion_elem, coefficients=coefficients, &
      oldvectors=oldvectors, buildmatrix=buildmatrix )

    call add_boundary_elements ( mesh, problem, rhsd, curve=2, &
      elemsub=poisson_natboun, coefficients=coefficients )

!   instationary term ( gamma du/dt )

    if ( step == 1 ) then
!     Euler implicit
      solhat%u = soln%u
      mfac = 1
    else
!     2nd order Gear (implicit)
      solhat%u = alpha0 * soln%u + alpha1 * solnm1%u
      mfac = gamma0
    end if

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_dcdt_elem, coefficients=coefficients, &
      oldvectors=oldvectors, factormat=mfac, addmatvec=.true., &
      buildmatrix=buildmatrix )

!   convection term

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_ugradc_elem, coefficients=coefficients, &
      oldvectors=oldvectors, addmatvec=.true., &
      buildmatrix=buildmatrix, buildvector=.false. )

!   beta * c term (relaxation,reaction)

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_betac_elem, coefficients=coefficients, &
      oldvectors=oldvectors, addmatvec=.true., &
      buildmatrix=buildmatrix, buildvector=.false. )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    call solve_system_ma41 ( sysmatrix, rhsd, sol, lu )

    if ( step == 1 ) call delete(lu) ! perform LU-decomposition at second step

    write(*,'(a,i0,a,es11.4,a,es15.8/)') &
      'step = ', step, ' time = ', step*deltat, &
      ' average sol = ', sum( abs(sol%u) ) / sol%n

    call copy ( soln, solnm1 )
    call copy ( sol, soln )

  end do

! plot solution

  call plot_color_contour ( plot_options, mesh, problem, filename='sol.fig', &
    sysvector=sol )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, soln, solnm1, solhat )
  call delete ( rhsd )
  if ( coeff_alpha == 2 ) call delete ( alpha_nodes )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program convection_diffusion4

