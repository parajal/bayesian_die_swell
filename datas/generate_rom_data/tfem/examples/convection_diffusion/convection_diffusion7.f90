! Steady scalar convection-diffusion problem in 1D on an interval (0,L):
!
!   u dc/dx - d/dx ( alpha * dc/dx ) = 0 on (0,L)
!
! with c(0)=c0 and c(L)=cL, using a spectral element method with a uniform grid.
! The standard elements are used.
!
program convection_diffusion7

  use tfem_m
  use hsl_ma41_m
  use convection_diffusion_elements_m
  use printtofile_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 13,        & ! scalar interpolation
    p = 9,              & ! polynomial order
    ne = 5 ,            & ! number of elements
    inttype = 1,        & ! Gauss-Legendre-Lobatto for integration
    gauss = p+1           ! number of Gauss points

  real(dp), parameter :: &
    alpha = 2.e-2_dp, &  ! diffusion coefficient
    gamma = 1._dp, &  ! coefficient of the time derivative "rho c_p"
    u = 1._dp,     &  ! convection speed
    c0 = 0._dp,    &  ! value at x=0
    cL = 1._dp,    &  ! value at x=L
    L = 1._dp         ! lengh of the domain

! variables

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(coefficients_t) :: coefficients

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(10) = gauss
  coefficients%i(39) = p
  coefficients%i(40) = inttype

  coefficients%r = 0
  coefficients%r(1) = alpha
  coefficients%r(4) = u
  coefficients%r(8) = gamma

! some info on Peclet number
  print *
  print *, 'Pe= ', abs(u)*L/alpha, 'Peh=', abs(u)*L/alpha/ne/2


! create mesh

  meshgen_options%elshape = 101 ! spectral line elements

  meshgen_options%nx = ne  ! number of elements, equidistant
  meshgen_options%lx = L   ! length of interval
  meshgen_options%p = p    ! polynomial order

  call line1d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )


! problem definition

  call create_input_probdef ( mesh, input_probdef )

  input_probdef%elementdof(1)%a = 1

  call define_essential ( mesh, input_probdef, point=1 )
  call define_essential ( mesh, input_probdef, point=2 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=1, value=c0 )
  call fill_sysvector ( mesh, problem, sol, point=2, value=cL )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! build diffusion matrix and right-hand side

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_diffusion_elem, coefficients=coefficients )

! convection term

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_ugradc_elem, coefficients=coefficients, &
    addmatvec=.true., buildvector=.false. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma41 ( sysmatrix, rhsd, sol )


! write output to file

  call printtofile ( mesh, problem, filename='sol.out', sysvector=sol )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol )
  call delete ( rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program convection_diffusion7

