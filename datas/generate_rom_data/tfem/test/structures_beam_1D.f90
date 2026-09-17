! A simply supported or built-in 1D beam problem
! Different distributed load on parts of the beam using two groups of
! elements with difference coefficients. Alternatively a function can be
! defined for w.
! P4 interpolation

program structures7

  use tfem_m
  use hsl_ma57_m
  use structures_elements_m
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
    M = 0._dp,    &  ! Moment at x=a
    w1 = 1._dp,   &  ! distributed load on the left part
    w2 = 0._dp       ! distributed load on the right part

! definitions

  type(mesh_t) :: mesh1, mesh2, mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd, reacf
  type(vector_t) :: dvec
  type(meshgen_options_t) :: mesh_options
  type(coefficients_t) :: coefficients(2)
  type(oldvectors_t) :: oldvectors

  integer :: i


! fill coefficients

  call create_coefficients ( coefficients(1), ncoefi=100, ncoefr=50 )

  coefficients(1)%i = 0
  coefficients(1)%i(2) = ninti
  coefficients(1)%i(3) = 1 ! P4 Hermite interpolation
  coefficients(1)%i(4) = 0 ! constant distributed force

  coefficients(1)%r = 0
  coefficients(1)%r(2) = Emod
  coefficients(1)%r(3) = Iz
  coefficients(1)%r(4) = -w1

  coefficients(2) = coefficients(1)
  coefficients(2)%r(4) = -w2


! create mesh

  mesh_options%elshape = 2 ! three-node line elements

  mesh_options%nx = ne1  ! number of elements, equidistant
  mesh_options%lx = a    ! length of interval

  call line1d ( mesh1, mesh_options )

  mesh_options%nx = ne2  ! number of elements, equidistant
  mesh_options%lx = L-a  ! length of interval
  mesh_options%ox = a    ! origin of interval

  call line1d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh, point1=2, point2=1, &
    nogroupmerge=.true. )

  call fill_mesh_parts ( mesh )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  do i = 1, 2
    input_probdef%elementdof(i)%a = [ 2, 1, 2 ]
    input_probdef%vec_elementdof(i)%a = 1
  end do

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
  call create_sysvector ( problem, rhsd, reacf )

! fill solution vector with essential boundary conditions

  if ( builtin ) then
    call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )
    call fill_sysvector ( mesh, problem, sol, point=3, value=0._dp )
  else
    call fill_sysvector ( mesh, problem, sol, point=1, degfd=1, value=0._dp )
    call fill_sysvector ( mesh, problem, sol, point=3, degfd=1, value=0._dp )
  end if

! fill rhs vector with discrete forces

  rhsd%u = 0
  call fill_sysvector ( mesh, problem, rhsd, point=2, degfd=1, value=-F )
  call fill_sysvector ( mesh, problem, rhsd, point=2, degfd=2, value=M )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=beam_elem1, &
    order='ND', mcoefficients=coefficients, addvec=.true. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

! post processing

  call create ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call create ( problem, dvec, vec=1 )

  call derive_write_scalar ( 'vertical displacement' )
  call derive_write_scalar ( 'slope', append=.true. )

  call delete ( dvec)
  call create ( problem, dvec, vec=1, elementwise=.true. )

  call derive_write_scalar ( 'second derivative', append=.true. )
  call derive_write_scalar ( 'bending moment', append=.true. )
  call derive_write_scalar ( 'shear force', append=.true. )


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

  subroutine derive_write_scalar ( dataname, append )

    character (len=*), intent(in) :: dataname
    logical, intent(in), optional :: append
    integer :: i


    do i = 1, 2

      select case ( dataname )
      case('vertical displacement')
        coefficients(i)%i(1) = 1 ! vertical displacement
      case('slope')
        coefficients(i)%i(1) = 2 ! slope
      case('second derivative')
        coefficients(i)%i(1) = 3 ! second derivative
      case('bending moment')
        coefficients(i)%i(1) = 4 ! bending moment
      case('shear force')
        coefficients(i)%i(1) = 5 ! shear force
      case default
        write(*,'(/a/)') 'Error: wrong function number: '
        stop
      end select

    end do

!   derive

    call derive_vector ( mesh, problem, dvec, elemsub=beam_deriv1, &
      mcoefficients=coefficients, oldvectors=oldvectors )

!   write

  !  call write_scalar_gmsh ( mesh, problem, filename='structures4.msh', &
  !    dataname=dataname, vector=dvec, append=append )


    print *, maxval(abs(dvec%u))

  end subroutine derive_write_scalar

end program structures7

