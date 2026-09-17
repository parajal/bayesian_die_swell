! A two-beam structure problem: an L-shaped wrench. 3D beam elements.
! A force down in point 3. Example from the Mechanics course 4RA00.
! Gmsh mesh.
! P1/P3/P1 interpolation for axial displacement, transverse displacement
! and axial rotation, respectively.

program structures18

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
    Emod = 1._dp, &  ! E modulus
    Gmod = 1._dp, &  ! G modulus
    Ac = 1._dp,   &  ! cross-sectional area
    Izz = 1._dp,  &  ! second moment of area Izz
    Iyy = 2._dp,  &  ! second moment of area Iyy
    Jrr = 1._dp,  &  ! polar moment of area J
!   external force/moments in point 3
    F = 1._dp  ! external force vector

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: load, rhsd, reacf
  type(coefficients_t) :: coefficients(2)

  integer :: i


! fill coefficients

  call create_coefficients ( coefficients(1), ncoefi=100, ncoefr=50 )

  coefficients(1)%i = 0
  coefficients(1)%i(2) = ninti
  coefficients(1)%i(3) = 0 ! P3 Hermite interpolation for v
  coefficients(1)%i(4) = 0 ! constant distributed force
  coefficients(1)%i(12) = 2 ! P1 interpolation for u
  coefficients(1)%i(13) = 2 ! P1 interpolation for theta

  coefficients(1)%r = 0
  coefficients(1)%r(1) = Ac
  coefficients(1)%r(2) = Emod
  coefficients(1)%r(3) = Izz
  coefficients(1)%r(8:10) = [ 0._dp, -1._dp, 0._dp ] ! local z-axis direction
  coefficients(1)%r(11) = Iyy
  coefficients(1)%r(12) = 0   ! Izy = 0
  coefficients(1)%r(13) = Jrr
  coefficients(1)%r(14) = Gmod

  coefficients(2) = coefficients(1)
  coefficients(2)%r(8:10) = [ 1._dp, 0._dp, 0._dp ] ! local z-axis direction


! read mesh

  call read_mesh_gmsh ( mesh, filename='ikea.msh', ndim=3 )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  do i = 1, 2
    input_probdef%elementdof(i)%a = 6 ! (u,v,w,phix,phiy,phiz) in two end nodes
    input_probdef%vec_elementdof(i)%a = 1      ! scalar quantity
    input_probdef%vec_elementdof(i)%a(:,2) = 3 ! vector quantity
  end do

! builtin at points 1

  call define_essential ( mesh, input_probdef, point=1 )

  call problem_definition ( input_probdef, mesh, problem )

  call printinfo ( mesh, printlevel=4 )
  call printinfo ( problem, printlevel=4 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, load, rhsd, reacf )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )

! fill rhs vector with discrete forces/monents at point 2

  load%u = 0

  call fill_sysvector ( mesh, problem, load, point=3, degfd=3, value=-F )

  call copy ( load, rhsd )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=beam_elem3, &
    order='ND', mcoefficients=coefficients, addvec=.true. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

! post processing

  call postprocessing_beam3 ( mesh, problem, mcoefficients=coefficients, &
    mshfilename='structures18.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures18

