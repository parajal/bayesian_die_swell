! A simple 3D truss problem.

program structures2

  use tfem_m
  use hsl_ma57_m
  use structures_elements_m
  use io_utils_m
  use limits_m

  implicit none


! constants

  real(dp), parameter :: &
    Emod = 1.e11_dp, &  ! E modulus
    A = 1.e-4_dp, &     ! A area
    F = -1.e4_dp        ! Force


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd, reacf
  type(vector_t) :: dvec
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors

  WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = .false.

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0

  coefficients%r = 0
  coefficients%r(1) = A
  coefficients%r(2) = Emod


! read mesh

  call read_mesh_gmsh ( mesh, filename='structures_truss_3D.msh' )

  call fill_mesh_parts ( mesh )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 3
  input_probdef%vec_elementdof(1)%a = 1

  call define_essential ( mesh, input_probdef, point=1 )
  call define_essential ( mesh, input_probdef, point=8 )
  call define_essential ( mesh, input_probdef, point=6, degfd=[0,1] )
  call define_essential ( mesh, input_probdef, point=11, degfd=[0,1] )

  call problem_definition ( input_probdef, mesh, problem )

  call printinfo ( mesh, printlevel=4 )
  call printinfo ( problem, printlevel=4 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd, reacf )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=8, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=6, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=11, degfd=2, value=0._dp )

! fill rhs vector with discrete forces

  rhsd%u = 0
  call fill_sysvector ( mesh, problem, rhsd, point=2, degfd=2, value=F )

! write load

  !call write_vector_gmsh ( mesh, problem, filename='structures2.msh', &
    !dataname='load', sysvector=rhsd )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=truss_elem1, &
    order='ND', coefficients=coefficients, buildvector=.false. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

! post processing

  call create ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call create ( problem, dvec, vec=1, elementwise=.true. )

  call derive_write_scalar ( 'strain' )
  call derive_write_scalar ( 'stress' )
  call derive_write_scalar ( 'force' )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf )
  call delete ( dvec )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

contains

! generic routine for post-processing data per element

  subroutine derive_write_scalar ( dataname )

    character (len=*), intent(in) :: dataname

    select case ( dataname )
    case('strain')
      coefficients%i(1) = 1 ! strain
    case('stress')
      coefficients%i(1) = 2 ! stress
    case('force')
      coefficients%i(1) = 3 ! force
      case default
        write(*,'(/a/)') 'Error: wrong function number: '
        stop
    end select

!   derive

    call derive_vector ( mesh, problem, dvec, elemsub=truss_deriv1, &
      coefficients=coefficients, oldvectors=oldvectors )

!   write

!    call write_scalar_gmsh ( mesh, problem, filename='structures2.msh', &
!      dataname=dataname, vector=dvec, append=.true. )

    print *, maxval(abs(dvec%u))

  end subroutine derive_write_scalar

end program structures2

