! A simply supported or built-in 1D beam problem
! P3 interpolation

program structures5

  use tfem_m
  use hsl_ma57_m
  use structures_elements_m
  use structures_post_m
  use io_utils_m

  implicit none


! constants

  logical :: &
    builtin = .true.

  integer, parameter :: &
    ninti = 3,  & ! number of Gauss points
    ne1 = 10,   & ! number of elements left
    ne2 = 10      ! number of elements right

  real(dp), parameter :: &
    Emod = 1._dp, &  ! E modulus
    Iz = 1._dp,   &  ! second moment of area
    L = 1._dp,    &  ! length of the beam
    a = 0.5_dp,   &  ! position of the force F and/or moment M
    F = 1._dp,    &  ! Force at x=a
    M = 2._dp,    &  ! Moment at x=a
    q = 0._dp        ! distributed load

! definitions

  type(mesh_t) :: mesh1, mesh2, mesh
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
  coefficients%i(3) = 0 ! P3 Hermite interpolation
  coefficients%i(4) = 0 ! constant distributed force

  coefficients%r = 0
  coefficients%r(2) = Emod
  coefficients%r(3) = Iz
  coefficients%r(4) = -q


! create mesh

  mesh_options%elshape = 1 ! two-node line elements

  mesh_options%nx = ne1  ! number of elements, equidistant
  mesh_options%lx = a    ! length of interval

  call line1d ( mesh1, mesh_options )

  mesh_options%nx = ne2  ! number of elements, equidistant
  mesh_options%lx = L-a  ! length of interval
  mesh_options%ox = a    ! origin of interval

  call line1d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh, point1=2, point2=1 )

  call fill_mesh_parts ( mesh )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  input_probdef%elementdof(1)%a = 2          ! (v,phi) in two end nodes
  input_probdef%vec_elementdof(1)%a = 1      ! scalar quantity
  input_probdef%vec_elementdof(1)%a(:,2) = 2 ! vector quantity

  if ( builtin ) then
    call define_essential ( mesh, input_probdef, point=1 )
    call define_essential ( mesh, input_probdef, point=3 )
  else
    call define_essential ( mesh, input_probdef, point=1, degfd=[1,0] )
    call define_essential ( mesh, input_probdef, point=3, degfd=[1,0] )
  end if

  call problem_definition ( input_probdef, mesh, problem )

  call printinfo ( mesh, printlevel=4 )
  call printinfo ( problem, printlevel=4 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd, reacf, load )

! fill solution vector with essential boundary conditions

  if ( builtin ) then
    call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )
    call fill_sysvector ( mesh, problem, sol, point=3, value=0._dp )
  else
    call fill_sysvector ( mesh, problem, sol, point=1, degfd=1, value=0._dp )
    call fill_sysvector ( mesh, problem, sol, point=3, degfd=1, value=0._dp )
  end if

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
    mshfilename='structures5.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures5

