module poisson_functions_m

  use math_defs_m

  implicit none

contains

  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    select case(nr)
      case(1)
        func = 0
      case(2)
        func = 1
      case(3)
        func = 1 + cos(pi*x(1))*cos(pi*x(2))
      case(4)
        func = 2*pi**2*cos(pi*x(1))*cos(pi*x(2))
      case(5)
        func = -pi*sin(pi*x(1))*cos(pi*x(2))
      case(6)
        func = -pi*cos(pi*x(1))*sin(pi*x(2))
      case(7)
        func = cos(pi*x(1))*cos(pi*x(2)) + (x(1)*x(2))**3
      case(8)
        func = 2*pi**2*cos(pi*x(1))*cos(pi*x(2)) &
                    -6*(x(1)*x(2)**3+x(1)**3*x(2))
      case(9)
        func = cos(pi*x(1))*cos(pi*x(2))*cos(pi*x(3)) + (x(1)*x(2)*x(3))**3
      case(10)
        func = 3*pi**2*cos(pi*x(1))*cos(pi*x(2))*cos(pi*x(3)) &
                    -6*(x(1)*x(2)**3*x(3)**3+x(1)**3*x(2)*x(3)**3 + &
                        x(1)**3*x(2)**3*x(3) )
      case(11)
        func = cos(2*pi*x(1)+1._dp)*cos(pi*x(2))*cos(pi*x(3)) + (x(2)*x(3))**3
      case(12)
        func = 6*pi**2*cos(2*pi*x(1)+1._dp)*cos(pi*x(2))*cos(pi*x(3)) &
                    -6*( x(2)*x(3)**3 + x(2)**3*x(3) )
      case(13)
        func = - ( -pi*sin(pi*x(1))*cos(pi*x(2)) + 3*x(1)**2*x(2)**3  &
                    + cos(pi*x(1))*cos(pi*x(2)) + (x(1)*x(2))**3 )
      case default
        write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
        stop
    end select

  end function func

  function vfunc ( n, nr, x )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    select case(nr)
      case(1)
        vfunc = - [ -pi*sin(pi*x(1))*cos(pi*x(2)) + 3*x(1)**2*x(2)**3,  &
                     -pi*cos(pi*x(1))*sin(pi*x(2)) + 3*x(1)**3*x(2)**2 ]
      case(2)
        vfunc = - [ -pi*sin(pi*x(1))*cos(pi*x(2)) ,  &
                     -pi*cos(pi*x(1))*sin(pi*x(2)) ]
      case(3)
        vfunc = - [ -pi*sin(pi*x(1))*cos(pi*x(2))*cos(pi*x(3)) &
                     + 3*x(1)**2*x(2)**3*x(3)**3,  &
                     -pi*cos(pi*x(1))*sin(pi*x(2))*cos(pi*x(3)) &
                     + 3*x(1)**3*x(2)**2*x(3)**3, &
                     -pi*cos(pi*x(1))*cos(pi*x(2))*sin(pi*x(3)) &
                     + 3*x(1)**3*x(2)**3*x(3)**2 ]
      case default
        write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
        stop
    end select

  end function vfunc

end module poisson_functions_m



module convection_diffusion_functions_m

  use math_defs_m

  implicit none

contains

!  function alphafunc ( nr, x )
!    integer, intent(in) :: nr
!    real(dp), intent(in), dimension(:) :: x
!    real(dp) :: alphafunc
!
!    select case(nr)
!      case(1)
!        alphafunc = 1.5_dp + cos(pi*x(1))*cos(pi*x(2))
!      case default
!        write(*,'(/a,i0/)') 'Error alphafunc: wrong function number: ', nr
!        stop
!    end select

!  end function alphafunc

  function alphatfunc ( n, nr, x )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n,n) :: alphatfunc

    real(dp) :: c(2,2)

    select case(nr)
      case(1)
        c(:,1) = [1._dp,-0.1_dp]
        c(:,2) = [-0.1_dp,2._dp]
        alphatfunc = c * ( 1.5_dp + cos(pi*x(1))*cos(pi*x(2)) )
      case(2)
        c(:,1) = [1._dp,0._dp]
        c(:,2) = [0._dp,1._dp]
        alphatfunc = c * ( 1.5_dp + cos(pi*x(1))*cos(pi*x(2)) )
      case default
        write(*,'(/a,i0/)') 'Error alphatfunc: wrong function number: ', nr
        stop
    end select

  end function alphatfunc

  function velocityvfunc ( n, nr, x )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: velocityvfunc

    select case(nr)
      case(1)
        velocityvfunc = [100._dp,50._dp] * &
                             ( 1.5_dp + cos(pi*x(1))*cos(pi*x(2)))
      case default
        write(*,'(/a,i0/)') 'Error velocityvfunc: wrong function number: ', nr
        stop
    end select

  end function velocityvfunc

end module convection_diffusion_functions_m

! Time-dependent scalar convection-diffusion problem on a unit square with
! Dirichlet boundary conditions and natural boundary conditions.
! Varying alpha tensor coefficient using a function.
! Varying beta coefficient using a function or a vector in the nodes.
! Varying velocity coefficient using a function.
! Crank-Nicolson time discretization.

program convection_diffusion6

  use tfem_m
  use hsl_ma41_m
  use convection_diffusion_elements_m
  use convection_diffusion_functions_m
  use poisson_functions_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3 integration of quads
    gaussb = 3,         & ! 3 point integration of boundary elements
    nx=10,             & ! number of elements in x
    ny=10,             & ! number of elements in y
    numtimesteps=2,    & ! number of time steps
    funcnr=8,           & ! function number for the right-hand side
    tfuncnr_alpha=1,    & ! function number for the alpha coefficient
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
    beta  = 200._dp, & ! beta (relaxation,reaction) coefficient
    gamma = 1._dp, &  ! coefficient of the time derivative "rho c_p"
    deltat = 1.e-3_dp ! time step

  real(dp), parameter :: &
    velocity(2)=[100._dp,10._dp], &  ! constant velocity for coeff_velo=0
!   constant alpha for coeff_alpha=0 (column wise storage)
    alphaten(4)=[1._dp,-0.1_dp,-0.1_dp,2._dp]


! variables

  logical :: buildmatrix

  integer :: step

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(sysvector_t), target :: soln
  type(vector_t), target :: alpha_nodes
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors
  type(lu_ma41_t) :: lu


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
  coefficients%i(26) = coeff_beta
  coefficients%i(27) = funcnr_beta

  coefficients%r = 0
  coefficients%r(3) = deltat
  coefficients%r(4:5) = velocity
  coefficients%r(7) = beta
  coefficients%r(8) = gamma
  coefficients%r(9:12) = alphaten

  coefficients%func => func
  coefficients%vfunc => vfunc
  coefficients%vfunc1(1)%p => velocityvfunc
  coefficients%tfunc1(1)%p => alphatfunc

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

  call define_essential ( mesh, input_probdef, curves=[1,3,4] )
!  call define_essential ( mesh, input_probdef, curve1=1 )
!  call define_essential ( mesh, input_probdef, curve1=3, curve2=4 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, soln )
  call create_sysvector ( problem, rhsd )

! fill initial solution vector with zero

  soln%u = 0

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curves=[1,3,4], func=func, funcnr=7 )
!  call fill_sysvector ( mesh, problem, sol, &
!    curve1=1, func=func, funcnr=7 )
!  call fill_sysvector ( mesh, problem, sol, &
!    curve1=3, curve2=4, func=func, funcnr=7 )

! oldvectors

  call create ( oldvectors, nsysvec=1, nvec=1 )

  oldvectors%s(1)%p => soln

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

    if ( step == 2 ) buildmatrix = .false.  ! matrix is constant

!   build diffusion matrix and right-hand side

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_diffusion2_elem, coefficients=coefficients, &
      oldvectors=oldvectors, factormat=0.5_dp, buildmatrix=buildmatrix )

!   build diffusion rhs

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_diffusion2_rhs_elem, coefficients=coefficients, &
      oldvectors=oldvectors, factorvec=0.5_dp, buildmatrix=.false., &
      addmatvec=.true. )

    call add_boundary_elements ( mesh, problem, rhsd, curve=2, &
      elemsub=poisson_natboun, coefficients=coefficients )

!   instationary term ( gamma du/dt )

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_dcdt_elem, coefficients=coefficients, &
      oldvectors=oldvectors, addmatvec=.true., buildmatrix=buildmatrix )

!   convection term

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_ugradc_elem, coefficients=coefficients, &
      oldvectors=oldvectors, addmatvec=.true., factormat=0.5_dp, &
      buildmatrix=buildmatrix, buildvector=.false. )

!   convection term (rhs)

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_ugradcn_elem, coefficients=coefficients, &
      oldvectors=oldvectors, addmatvec=.true., factorvec=0.5_dp, &
      buildmatrix=.false. )

!   beta * c term (relaxation,reaction)

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_betac_elem, coefficients=coefficients, &
      oldvectors=oldvectors, addmatvec=.true., factormat=0.5_dp, &
      buildmatrix=buildmatrix, buildvector=.false. )

!   beta * c term (relaxation,reaction) right-hand side

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_betacn_elem, coefficients=coefficients, &
      oldvectors=oldvectors, addmatvec=.true., factorvec=0.5_dp, &
      buildmatrix=.false. )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    call solve_system_ma41 ( sysmatrix, rhsd, sol, lu )

    write(*,'(a,i0,a,es10.3,a,es23.14/)') &
      'step = ', step, ' time = ', step*deltat, &
      ' average sol = ', sum( abs(sol%u) ) / sol%n

    call copy ( sol, soln )

  end do

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, soln )
  call delete ( rhsd )
  if ( coeff_alpha == 2 ) call delete ( alpha_nodes )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program convection_diffusion6

