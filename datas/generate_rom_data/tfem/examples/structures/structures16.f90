! A two-beam structure problem. Gravity load. 3D beam elements.
! Forces and moments in point 2 (tip of triangle).
! Gmsh mesh.
! P1/P3/P1 interpolation for axial displacement, transverse displacement
! and axial rotation, respectively.
! Extension of structures11 to a 3D beam structure.

program structures16

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
    rho = 0._dp,  & ! density
    g = 1._dp,    & ! gravitational acceleration
    Emod = 1._dp, &  ! E modulus
    Gmod = 1._dp, &  ! G modulus
    Ac = 1._dp,   &  ! cross-sectional area
    Izz = 1._dp,  &  ! second moment of area Izz
    Iyy = 2._dp,  &  ! second moment of area Iyy
    Jrr = 1._dp,  &  ! polar moment of area J
!   external forces/moments in point 2
    F(3) = [ 1._dp, 2._dp, 10._dp ], &  ! external force vector
    M(3) = [ 10._dp, 0._dp, 3._dp ]    ! external moment

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
  coefficients%i(13) = 2 ! P1 interpolation for theta

  coefficients%r = 0
  coefficients%r(1) = Ac
  coefficients%r(2) = Emod
  coefficients%r(3) = Izz
  coefficients%r(5:7) = [ 0._dp, -rho*g*Ac, 0._dp ]
  coefficients%r(8:10) = [ 0._dp, 0._dp, 1._dp ] ! local z-axis direction
  coefficients%r(11) = Iyy
  coefficients%r(12) = 0   ! Izy = 0
  coefficients%r(13) = Jrr
  coefficients%r(14) = Gmod


! read mesh

  call read_mesh_gmsh ( mesh, filename='triangle.msh', ndim=3 )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  input_probdef%elementdof(1)%a = 6  ! (u,v,w,phix,phiy,phiz) in two end nodes
  input_probdef%vec_elementdof(1)%a = 1      ! scalar quantity
  input_probdef%vec_elementdof(1)%a(:,2) = 3 ! vector quantity

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
  call fill_sysvector ( mesh, problem, load, point=2, degfd=3, value=F(3) )
  call fill_sysvector ( mesh, problem, load, point=2, degfd=4, value=M(1) )
  call fill_sysvector ( mesh, problem, load, point=2, degfd=5, value=M(2) )
  call fill_sysvector ( mesh, problem, load, point=2, degfd=6, value=M(3) )

  call copy ( load, rhsd )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=beam_elem3, &
    order='ND', coefficients=coefficients, addvec=.true. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

! post processing

  call postprocessing_beam3 ( mesh, problem, coefficients=coefficients, &
    mshfilename='structures16.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures16

