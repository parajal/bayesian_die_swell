! A simple 1D cantilever beam problem
! P4 interpolation

program structures4

  use tfem_m
  use hsl_ma57_m
  use structures_elements_m
  use structures_post_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    ninti = 3 , & ! number of Gauss points
    ne = 5        ! number of elements

  real(dp), parameter :: &
    Emod = 1_dp, &  ! E modulus
    Iz = 1_dp, &    ! second moment of area
    L =  1_dp, &    ! length of the beam
    F = 1_dp, &     ! Force at x=L
    M = 2_dp, &     ! Moment at x=L
    q = 1_dp        ! distributed load

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd, reacf, load
  type(meshgen_options_t) :: mesh_options
  type(coefficients_t) :: coefficients


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(2) = ninti
  coefficients%i(3) = 1 ! P4 Hermite interpolation
  coefficients%i(4) = 0 ! constant distributed force

  coefficients%r = 0
  coefficients%r(2) = Emod
  coefficients%r(3) = Iz
  coefficients%r(4) = -q


! create mesh

  mesh_options%elshape = 2 ! three-node line elements

  mesh_options%nx = ne  ! number of elements, equidistant
  mesh_options%lx = L  ! length of interval

  call line1d ( mesh, mesh_options )

  call fill_mesh_parts ( mesh )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4 )

! (v,phi) in two end nodes + bubble v in the middle node
  input_probdef%elementdof(1)%a = [ 2, 1, 2 ]
  input_probdef%vec_elementdof(1)%a = 1      ! scalar quantity
  input_probdef%vec_elementdof(1)%a(:,2) = 2 ! vector quantity
  input_probdef%vec_elementdof(1)%a(:,3) = 2 ! all nodes all degrees
  input_probdef%vec_elementdof(1)%a(:,4) = [ 2, 0, 2 ] ! end nodes all degrees

  call define_essential ( mesh, input_probdef, point=1 )

  call problem_definition ( input_probdef, mesh, problem )

  call printinfo ( mesh, printlevel=4 )
  call printinfo ( problem, printlevel=4 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd, reacf, load )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )

! fill rhs vector with discrete forces

  load%u = 0

  call fill_sysvector ( mesh, problem, load, point=2, degfd=1, value=-F )
  call fill_sysvector ( mesh, problem, load, point=2, degfd=2, value=M )

  call copy ( load, rhsd )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=beam_elem1, &
    order='ND', coefficients=coefficients, addvec=.true. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

! post processing

  call postprocessing_beam1 ( mesh, problem, coefficients=coefficients, &
    mshfilename='structures4.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures4

