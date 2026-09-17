! Poisson problem on a unit cube with Dirichlet boundary conditions and
! periodical boundary conditions in x-direction (weak connection)

module poisson_functions_m

  use math_defs_m

  implicit none

contains

  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    select case(nr)
      case(1)
        func = 0
      case(2)
        func = 1
      case(3)
        func = 1 + cos(pi*x(1))*cos(pi*x(2))
      case(4)
        func = 2*pi**2*cos(pi*x(1))*cos(pi*x(2))
      case(5)
        func = -pi*sin(pi*x(1))*cos(pi*x(2))
      case(6)
        func = -pi*cos(pi*x(1))*sin(pi*x(2))
      case(7)
        func = cos(pi*x(1))*cos(pi*x(2)) + (x(1)*x(2))**3
      case(8)
        func = 2*pi**2*cos(pi*x(1))*cos(pi*x(2)) &
                    -6*(x(1)*x(2)**3+x(1)**3*x(2))
      case(9)
        func = cos(pi*x(1))*cos(pi*x(2))*cos(pi*x(3)) + (x(1)*x(2)*x(3))**3
      case(10)
        func = 3*pi**2*cos(pi*x(1))*cos(pi*x(2))*cos(pi*x(3)) &
                    -6*(x(1)*x(2)**3*x(3)**3+x(1)**3*x(2)*x(3)**3 + &
                        x(1)**3*x(2)**3*x(3) )
      case(11)
        func = cos(2*pi*x(1)+1._dp)*cos(pi*x(2))*cos(pi*x(3)) + (x(2)*x(3))**3
      case(12)
        func = 6*pi**2*cos(2*pi*x(1)+1._dp)*cos(pi*x(2))*cos(pi*x(3)) &
                    -6*( x(2)*x(3)**3 + x(2)**3*x(3) )
      case default
        write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
        stop
    end select

  end function func

  function vfunc ( n, nr, x )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    select case(nr)
      case(1)
        vfunc = [ -pi*sin(pi*x(1))*cos(pi*x(2)) + 3*x(1)**2*x(2)**3,  &
                   -pi*cos(pi*x(1))*sin(pi*x(2)) + 3*x(1)**3*x(2)**2 ]
      case(2)
        vfunc = [ -pi*sin(pi*x(1))*cos(pi*x(2)) ,  &
                   -pi*cos(pi*x(1))*sin(pi*x(2)) ]
      case(3)
        vfunc = [ -pi*sin(pi*x(1))*cos(pi*x(2))*cos(pi*x(3)) &
                   + 3*x(1)**2*x(2)**3*x(3)**3,  &
                   -pi*cos(pi*x(1))*sin(pi*x(2))*cos(pi*x(3)) &
                   + 3*x(1)**3*x(2)**2*x(3)**3, &
                   -pi*cos(pi*x(1))*cos(pi*x(2))*sin(pi*x(3)) &
                   + 3*x(1)**3*x(2)**3*x(3)**2 ]
      case default
        write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
        stop
    end select

  end function vfunc

end module poisson_functions_m



program poisson7

  use tfem_m
  use hsl_ma57_m
  use poisson_elements_m
  use poisson_functions_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3x3 integration of quads
    gaussb = 3,         & ! 3x3 point integration of boundary elements
    nx=3,               & ! number of elements in x
    ny=3,               & ! number of elements in y
    nz=3,               & ! number of elements in z
    funcnr=12             ! function number for the right-hand side

  real(dp), parameter :: &
    alpha = 1._dp     ! diffusion coefficient

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd, solexact
  type(meshgen_options_t) :: meshgen_options
  type(coefficients_t) :: coefficients

  integer :: n


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0

  coefficients%func => func


! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, regionshape=4, &
    elshape=14 )

  call hexahedron ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef )

  input_probdef%elementdof(1)%a = 1

  call define_essential ( mesh, input_probdef, surfaces=[1,2,4,6] )
!  call define_essential ( mesh, input_probdef, surface1=1, surface2=2 )
!  call define_essential ( mesh, input_probdef, surface1=4 )
!  call define_essential ( mesh, input_probdef, surface1=6 )

  call define_constraint ( mesh, input_probdef, &
    surface1=3, surface2=5, discretization='weak', &
    elementdof=[1,0,1,0,1,0,1,0,0] )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    surfaces=[1,2,4,6], func=func, funcnr=11 )
!  call fill_sysvector ( mesh, problem, sol, &
!    surface1=1, surface2=2, func=func, funcnr=11 )
!  call fill_sysvector ( mesh, problem, sol, &
!    surface1=4, func=func, funcnr=11 )
!  call fill_sysvector ( mesh, problem, sol, &
!    surface1=6, func=func, funcnr=11 )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=poisson_elem, &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    elemsub=poisson_constr_elem_conn, addmatvec=.true., &
    coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! create exact solution

  call create_sysvector ( problem, solexact )

  call fill_sysvector ( mesh, problem, solexact, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=11 )

! print maximum difference of sol-solexact to standard output

  n = problem%numnodaldegfd

  print *, maxval( abs(sol%u(problem%degfdperm(:n,2)) - &
                       solexact%u(problem%degfdperm(:n,2)) ) )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, solexact, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program poisson7

