! A two-beam structure problem. Gravity load.
! Forces and moment in point 2 (tip of triangle).
! P1/P3 interpolation. 2D elements.
! Gmsh mesh.

program structures11

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
    rho = 1._dp,  & ! density
    g = 1._dp,    & ! gravitational acceleration
    Emod = 1._dp, &  ! E modulus
    Ac = 1._dp,   &  ! cross-sectional area
    Iz = 1._dp,   &  ! second moment of area
!   external forces/moments in point 2
    F(2) = [ 1._dp,  2._dp ], &  ! external force vector
    M = 3._dp     ! external moment

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

  call read_mesh_gmsh ( mesh, filename='triangle.msh', ndim=2 )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  input_probdef%elementdof(1)%a = 3          ! (u,v,phi) in two end nodes
  input_probdef%vec_elementdof(1)%a = 1      ! scalar quantity
  input_probdef%vec_elementdof(1)%a(:,2) = 2 ! vector quantity

! builtin at points 1 and 3

  call define_essential ( mesh, input_probdef, point=1 )
  call define_essential ( mesh, input_probdef, point=3 )

  call problem_definition ( input_probdef, mesh, problem )

  call printinfo ( mesh, printlevel=4 )
  call printinfo ( problem, printlevel=4 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, load, rhsd, reacf )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=3, value=0._dp )

! fill rhs vector with discrete forces/monents at point 2

  load%u = 0

  call fill_sysvector ( mesh, problem, load, point=2, degfd=1, value=F(1) )
  call fill_sysvector ( mesh, problem, load, point=2, degfd=2, value=F(2) )
  call fill_sysvector ( mesh, problem, load, point=2, degfd=3, value=M )

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
    mshfilename='structures11.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures11

