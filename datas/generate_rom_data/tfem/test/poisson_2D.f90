! Poisson problem on a unit square with Dirichlet boundary conditions and
! natural boundary conditions. Curved domain.

module functions_m

  use kind_defs_m

  implicit none

contains

  function cfunc ( nr, xr )
    integer, intent(in) :: nr
    real(dp), intent(in) :: xr
    real(dp), dimension(2) :: cfunc

    select case(nr)
      case(1)
        cfunc = [ 2*xr, 5*xr*(1-xr) ]
      case(2)
        cfunc = [ 2+xr*(1-xr), 2*xr ]
      case(3)
        cfunc = [ -1+3*(1-xr), 2+xr*(1-xr) ]
      case(4)
        cfunc = [ xr-1+xr*(1-xr), 3*(1-xr) ]
      case default
        write(*,'(/a,i0/)') 'Error cfunc: wrong function number: ', nr
        stop
    end select

  end function cfunc

end module functions_m

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
        vfunc = - [ -pi*sin(pi*x(1))*cos(pi*x(2)) ,  &
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

program poisson3

  use tfem_m
  use hsl_ma57_m
  use poisson_elements_m
  use poisson_functions_m
  use functions_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3 integration of quads
    gaussb = 3,         & ! 3 point integration of boundary elements
    nx=20,              & ! number of elements in x
    ny=20,              & ! number of elements in y
    funcnr=4,           & ! function number for the right-hand side
    vfuncnr=2             ! function number for the flux vector

  real(dp), parameter :: &
    alpha = 1._dp     ! diffusion coefficient

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd, solexact
  type(coefficients_t) :: coefficients


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]
  coefficients%i(16) = vfuncnr

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0

  coefficients%func => func
  coefficients%vfunc => vfunc


! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  meshgen_options%regionshape = 3
  meshgen_options%x2d = &
    reshape ( [ 0.0_dp, 2.0_dp, 2.0_dp, -1.0_dp,    &
                 0.0_dp, 0.0_dp, 2.0_dp,  3.0_dp ], &
              [4,2] )
  meshgen_options%curved(1:4) = .true.
  meshgen_options%funcnr(1:4) = [ 1, 2, 3, 4 ]

  call quadrilateral2d ( mesh, meshgen_options, func=cfunc )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef )

  input_probdef%elementdof(1)%a = [1,1,1,1,1,1,1,1,1]

  call define_essential ( mesh, input_probdef, curve1=1 )
  call define_essential ( mesh, input_probdef, curve1=3, curve2=4 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, func=func, funcnr=3 )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, curve2=4, func=func, funcnr=3 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=poisson_elem, &
    coefficients=coefficients )

  call add_boundary_elements ( mesh, problem, rhsd, curve=2, &
    elemsub=poisson_natboun, coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! create exact solution

  call create_sysvector ( problem, solexact )

  call fill_sysvector ( mesh, problem, solexact, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=3 )

! print maximum difference of sol-solexact to standard output

  print *, maxval( abs(sol%u -solexact%u) )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, solexact, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program poisson3

