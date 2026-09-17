
! solve:
!
!  u dc/dx - d/dx ( alpha * dc/dx ) + gamma c = f on (0,L)
!
! using a FEM (Galerkin, SU, SUPG) with a uniform grid
!

program convdiff1

  use tfem_m
  use element_m
  use hsl_ma41_m
  use printtofile_m

  implicit none

  type(meshgen_options_t) :: mesh_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: ipd
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(coefficients_t) :: coefficients
  type(lu_ma41_t) :: lu
  type(solver_options_ma41_t) :: solver_options

! parameters for the element

  integer :: isource = 2, ibeta=5, imeth = 2, ninti = 2
  real(dp) :: u = 1, alpha = 1.e-16_dp, gamma = 0, betac = 1, dcdxL = 1

! parameters for the problem

  integer :: ne = 10, iboun = 1
  real(dp) :: c0 = 0, cL = 1, L = 1

! namelist /elementpar/:
!
! isource   0: source=0
!           1: source=const
!           2: source=const/sqrt(pi)/sigma*exp(-(x-Lh)^2/sigma^2)
!           3: source=0 for x<Lh, const for x>=Lh
!
! ibeta   0: beta=0
!         1: beta=1
!         2: beta=0 for Peh<1
!                 1-1/Peh for Peh>1
!         3: beta=coth(Peh)-1/Peh
!         4: beta=Peh/3 for Peh<3
!                 1 for Peh>3
!         5: beta=betac
!
! imeth  0: Galerkin
!        1: SU
!        2: SUPG
!
! u      convection speed (>0)
! alpha  diffusion coefficient
! gamma  reaction coefficient
! betac  constant beta value for ibeta=5
! ninti  number of Gauss integration points
  namelist /elementpar/ isource, ibeta, imeth, u, alpha, gamma, betac, ninti

! problempar:
!
!  ne  number of elements, equidistant
!  L   length of the interval (interval: [0,L])
!  iboun  1: c(L) = cL (Dirichlet)
!         2: dcdx(L) = dcdxL (nat bc for flux A*dcdx)
!  c0  Dirichlet value at x=0
!  cL  Dirichlet value at x=L
!  dcdxL  Neumann value at x=L
!
  namelist /problempar/ ne, L, iboun, c0, cL, dcdxL


! read data with namelists

  open(unit=9,file='input_convdiff1.txt')

  read(unit=9,nml=sourcepar)
  read(unit=9,nml=elementpar)
  read(unit=9,nml=problempar)

  close(unit=9)


! write data with namelists

  write(*,nml=sourcepar)
  write(*,nml=elementpar)
  write(*,nml=problempar)
  print *
  print *, 'Pe= ', abs(u)*L/alpha, 'Peh=', abs(u)*L/alpha/ne/2


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=50, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1:4) = &
    [ isource, ibeta, imeth, ninti ]

  coefficients%r = 0
  coefficients%r(1:5) = [ u, alpha, gamma, betac, dcdxL ]


! create mesh

  mesh_options%elshape = 1 ! two-node line elements

  mesh_options%nx = ne  ! number of elements, equidistant
  mesh_options%lx = L  ! length of interval

  call line1d ( mesh, mesh_options )

  call fill_mesh_parts ( mesh )


! problem definition

  call create_input_probdef ( mesh, ipd )

  ipd%elementdof(1)%a = 1

  call define_essential ( mesh, ipd, point=1 )
  if ( iboun == 1 ) then
    call define_essential ( mesh, ipd, point=2 )
  end if

  call problem_definition ( ipd, mesh, problem )


! create system vectors

  call create_sysvector ( problem, sol, rhsd )


! fill solution vector with essential boundary conditions

  sol%u = 0

  call fill_sysvector ( mesh, problem, sol, point=1, value=c0 )
  if ( iboun == 1 ) then
    call fill_sysvector ( mesh, problem, sol, point=2, value=cL )
  end if


! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )


! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=element, &
    coefficients=coefficients )

  call check_filled_sysmatrix ( sysmatrix )

  if ( iboun == 2 ) then
    call add_boundary_elements_point ( mesh, problem, rhsd, point=2, &
      elemsub=bounelement, coefficients=coefficients )
  end if

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


! solve system

  solver_options%integer_storage = 2.0

  call solve_system_ma41 ( sysmatrix, rhsd, sol, lu, &
    solver_options=solver_options )


! write output to file

  call printtofile ( mesh, problem, filename='sol.out', sysvector=sol )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( ipd )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( lu )

end program convdiff1

