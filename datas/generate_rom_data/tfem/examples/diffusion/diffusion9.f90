! Time-dependent temperature problem in a cylindrical rod.
! This problem solves the example in Sec. 3.2.2. of Dantzig & Tucker:
! A solid cylindrical rod of radius R and length L, at an initial
! temperature T1, is quenched in a fluid at a lower temperature T0.
! Second-order implicit Gear discretization. First step: Euler implicit.
! Half the domain, due to symmetry.
! Plots as a function of time.

program diffusion9

  use tfem_m
  use hsl_ma57_m
  use convection_diffusion_elements_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3 integration of quads
    gaussb = 3,         & ! 3 point integration of boundary elements
    nz=50,              & ! number of elements in z-direction
    nr=10,              & ! number of elements in r-direction
    vtkevery = 1,       & ! write vtk file every vtkevery steps. -1: means none
    numtimesteps=80       ! number of time steps

  real(dp), parameter :: &
    L = 5._dp,  & ! length of the rod
    R = 1._dp,  & ! radius of the rod
    k_coef = 1._dp, & ! thermal conductivity
    rho_cp = 1._dp, & ! coefficient of the time derivative "rho c_p"
    T1 = 1._dp, & ! initial temperature of the rod
    T0 = 0._dp, & ! temperature of fluid bath
    deltat = 1.e-2_dp ! time step

  real(dp), parameter :: &
    gamma0=1.5_dp, alpha0=2._dp, alpha1=-0.5_dp  ! Gear parameters

! variables

  logical :: buildmatrix

  integer :: step, ipostv = 0

  real(dp) :: mfac

  character(len=30) :: filename

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: solnm1, rhsd
  type(sysvector_t), target :: sol, solhat, soln
  type(vector_t), target :: grad
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors
  type(lu_ma57_t) :: lu


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(10:11) = [ gauss, gaussb ]
  coefficients%i(23) = 1  ! axi-symmetrical

  coefficients%r = 0
  coefficients%r(1) = k_coef
  coefficients%r(3) = deltat
  coefficients%r(8) = rho_cp

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%lx = L/2
  meshgen_options%ly = R
  meshgen_options%nx = nz/2
  meshgen_options%ny = nr

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a(:,1) = 2  ! for gradT in the nodes

  call define_essential ( mesh, input_probdef, curve1=3, curve2=4 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, soln, solnm1, solhat )
  call create_sysvector ( problem, rhsd )

! fill initial temperature

  soln%u = T1

! create vector for gradient

  call create ( problem, grad, vec=1 )

! oldvectors

  call create ( oldvectors, nsysvec=1, nvec=1 )

! write a vtk file of initial condition

  if ( vtkevery > 0 ) then

    write(filename,'(a,i4.4,a)') 'sol', 0, '.vtk'
    call write_scalar_vtk ( mesh, problem, filename=filename, &
      dataname='temperature', sysvector=soln )

    oldvectors%s(1)%p => soln
    call derive_vector ( mesh, problem, grad, elemsub=poisson_deriv, &
      coefficients=coefficients, oldvectors=oldvectors )
    call write_vector_vtk ( mesh, problem, filename=filename, &
      dataname='gradT', vector=grad, append=.true. )

  end if

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, curve1=2, curve2=4, value=T0 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! start time stepping

  buildmatrix = .true.

  do step = 1, numtimesteps

    if ( step == 3 ) buildmatrix = .false.  ! matrix is constant

!   build diffusion matrix and right-hand side

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_diffusion_elem, coefficients=coefficients, &
      oldvectors=oldvectors, buildmatrix=buildmatrix )

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

    oldvectors%s(1)%p => solhat

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=scalar_dcdt_elem, coefficients=coefficients, &
      oldvectors=oldvectors, factormat=mfac, addmatvec=.true., &
      buildmatrix=buildmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu )

    if ( step == 1 ) call delete(lu) ! perform LU-decomposition at second step

    write(*,'(a,i0,a,es11.4,a,es15.8/)') &
      'step = ', step, ' time = ', step*deltat, &
      ' average sol = ', sum( abs(sol%u) ) / sol%n

!   write a vtk file

    if ( vtkevery > 0 ) then
      if ( mod(step,vtkevery) == 0 ) then
!
        ipostv = ipostv + 1
        write(filename,'(a,i4.4,a)') 'sol', ipostv, '.vtk'
        call write_scalar_vtk ( mesh, problem, filename=filename, &
          dataname='temperature', sysvector=sol )

        oldvectors%s(1)%p => sol
        call derive_vector ( mesh, problem, grad, elemsub=poisson_deriv, &
          coefficients=coefficients, oldvectors=oldvectors )
        call write_vector_vtk ( mesh, problem, filename=filename, &
          dataname='gradT', vector=grad, append=.true. )

      end if
    end if

    call copy ( soln, solnm1 )
    call copy ( sol, soln )

  end do

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, soln, solnm1, solhat )
  call delete ( rhsd )
  call delete ( grad )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program diffusion9

