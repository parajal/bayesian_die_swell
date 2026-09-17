! A beam structure problem of exam Mechanics 8-11-2019. Gravity load.
! Forces in point 2 and 3. Builtin at point 1 and 5
! P1/P3 interpolation. 2D elements.
! Gmsh mesh.

program structures12

  use tfem_m
  use hsl_ma57_m
  use structures_elements_m
  use structures_post_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    ninti = 3 ! number of Gauss points

  real(dp), parameter :: &
    rho = 7850._dp,  & ! density
    g = 9.8_dp,    & ! gravitational acceleration
    Emod = 200e9_dp, &  ! E modulus
    Ac = 1.e-4_dp,   &  ! cross-sectional area
    Iz = 1.e-8_dp/12,   &  ! second moment of area
!   external forces in point 2 and 3
    F = 100._dp, & ! value of F
    F2(2) = [ 3*F, -F ], & ! external force vector in point 2
    F3(2) = [ -2*F, -F ]   ! external force vector in point 3

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: load, rhsd, reacf
  type(coefficients_t) :: coefficients


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(2) = ninti
  coefficients%i(3) = 0 ! P3 Hermite interpolation for v
  coefficients%i(4) = 0 ! constant distributed force
  coefficients%i(12) = 2 ! P1 interpolation for u

  coefficients%r = 0
  coefficients%r(1) = Ac
  coefficients%r(2) = Emod
  coefficients%r(3) = Iz
  coefficients%r(5:6) = [ 0._dp, -rho*g*Ac ]


! read mesh

  call read_mesh_gmsh ( mesh, filename='portal.msh', ndim=2 )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  input_probdef%elementdof(1)%a = 3          ! (u,v,phi) in two end nodes
  input_probdef%vec_elementdof(1)%a = 1      ! scalar quantity
  input_probdef%vec_elementdof(1)%a(:,2) = 2 ! vector quantity

! builtin at points 1 and 5

  call define_essential ( mesh, input_probdef, point=1 )
  call define_essential ( mesh, input_probdef, point=5 )

  call problem_definition ( input_probdef, mesh, problem )

  call printinfo ( mesh, printlevel=4 )
  call printinfo ( problem, printlevel=4 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, load, rhsd, reacf )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=5, value=0._dp )

! fill rhs vector with discrete forces at point 1 and 3

  load%u = 0

  call fill_sysvector ( mesh, problem, load, point=2, degfd=1, value=F2(1) )
  call fill_sysvector ( mesh, problem, load, point=2, degfd=2, value=F2(2) )
  call fill_sysvector ( mesh, problem, load, point=3, degfd=1, value=F3(1) )
  call fill_sysvector ( mesh, problem, load, point=3, degfd=2, value=F3(2) )

  call copy ( load, rhsd )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=beam_elem2, &
    order='ND', coefficients=coefficients, addvec=.true. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

! post processing

  call postprocessing_beam2 ( mesh, problem, coefficients=coefficients, &
    mshfilename='structures12.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures12

