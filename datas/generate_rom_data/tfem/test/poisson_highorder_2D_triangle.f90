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
      case(13)
        func = - ( -pi*sin(pi*x(1))*cos(pi*x(2)) + 3*x(1)**2*x(2)**3  &
                    + cos(pi*x(1))*cos(pi*x(2)) + (x(1)*x(2))**3 )
      case(14)
        func = 10 - 8*2*log(3/x(2))  ! u1=10, R0=2, R1=3, h0=8
      case(15)
        func = 8
      case(16)
        func = cos(pi*x(1))*cos(pi*x(2))*cos(pi*x(3))
      case(17)
        func = 3*pi**2*cos(pi*x(1))*cos(pi*x(2))*cos(pi*x(3))
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
        vfunc = - [ -pi*sin(pi*x(1))*cos(pi*x(2)) + 3*x(1)**2*x(2)**3,  &
                     -pi*cos(pi*x(1))*sin(pi*x(2)) + 3*x(1)**3*x(2)**2 ]
      case(2)
        vfunc = - [ -pi*sin(pi*x(1))*cos(pi*x(2)) ,  &
                     -pi*cos(pi*x(1))*sin(pi*x(2)) ]
      case(3)
        vfunc = - [ -pi*sin(pi*x(1))*cos(pi*x(2))*cos(pi*x(3)) &
                     + 3*x(1)**2*x(2)**3*x(3)**3,  &
                     -pi*cos(pi*x(1))*sin(pi*x(2))*cos(pi*x(3)) &
                     + 3*x(1)**3*x(2)**2*x(3)**3, &
                     -pi*cos(pi*x(1))*cos(pi*x(2))*sin(pi*x(3)) &
                     + 3*x(1)**3*x(2)**3*x(3)**2 ]
      case(4)
        vfunc = [ -pi*sin(pi*x(1))*cos(pi*x(2))*cos(pi*x(3)) &
                   + 3*x(1)**2*x(2)**3*x(3)**3,  &
                   -pi*cos(pi*x(1))*sin(pi*x(2))*cos(pi*x(3)) &
                   + 3*x(1)**3*x(2)**2*x(3)**3, &
                   -pi*cos(pi*x(1))*cos(pi*x(2))*sin(pi*x(3)) &
                   + 3*x(1)**3*x(2)**3*x(3)**2 ]
      case(5)
        vfunc = - [ -pi*sin(pi*x(1))*cos(pi*x(2))*cos(pi*x(3)), &
                     -pi*cos(pi*x(1))*sin(pi*x(2))*cos(pi*x(3)), &
                     -pi*cos(pi*x(1))*cos(pi*x(2))*sin(pi*x(3)) ]
      case(6)
        vfunc = [ -pi*sin(pi*x(1))*cos(pi*x(2))*cos(pi*x(3)), &
                   -pi*cos(pi*x(1))*sin(pi*x(2))*cos(pi*x(3)), &
                   -pi*cos(pi*x(1))*cos(pi*x(2))*sin(pi*x(3)) ]
      case default
        write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
        stop
    end select

  end function vfunc

end module poisson_functions_m



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




! Poisson problem on a curved domain with Dirichlet boundary conditions and
! natural boundary conditions. Derived quantity.
! High-order triangular elements.
! This problem is similar to poisson4, but now with high-order elements.

program poisson37

  use tfem_m
  use hsl_ma57_m
  use poisson_elements_m
  use poisson_functions_m
  use functions_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 19,  & ! scalar interpolation
    p = 5,        & ! polynomial order
    gauss = 2*p-1,& ! number of Gauss points
    gaussb = p+1, & ! number of Gauss points on the boundary
    inttype = 3,  & ! numerical Gauss points
    nx=5,        & ! number of elements in x
    ny=5,        & ! number of elements in y
    funcnr=4,     & ! function number for the right-hand side
    vfuncnr=2       ! function number for the flux vector

  real(dp), parameter :: &
    alpha = 1._dp     ! diffusion coefficient

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh_plot
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd, solexact
  type(vector_t) :: grad, gradexact
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(6) = 2  ! use a vector function for h_N
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]
  coefficients%i(16) = vfuncnr
  coefficients%i(39) = p
  coefficients%i(40) = inttype

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0

  coefficients%func => func
  coefficients%vfunc => vfunc

! create mesh

  meshgen_options%elshape = 104
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%p = p
  meshgen_options%l = 0 ! equidistant nodal distribution

  meshgen_options%regionshape = 3
  meshgen_options%x2d = &
    reshape ( [ 0.0_dp, 2.0_dp, 2.0_dp, -1.0_dp,    &
                 0.0_dp, 0.0_dp, 2.0_dp,  3.0_dp ], &
              [4,2] )
  meshgen_options%curved(1:4) = [ .false., .true., .true., .true. ]
  meshgen_options%funcnr(1:4) = [ 1, 2, 3, 4 ]

  call quadrilateral2d ( mesh, meshgen_options, func=cfunc )

  call fill_mesh_parts ( mesh )

! Create plot mesh

  call mesh_convert ( mesh, mesh_plot, elementshapes='spectraltolinear', &
    warn=.false. )

  call fill_mesh_parts( mesh_plot )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1

  input_probdef%vec_elementdof(1)%a = 2

  call define_essential ( mesh, input_probdef, curves=[1,3,4] )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curves=[1,3,4], func=func, funcnr=3 )

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

! create grad

  call create ( problem, grad, vec=1 )

! create oldvectors

  call create ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, grad, elemsub=poisson_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! create exact solution of grad

  call create_vector ( problem, gradexact, vec=1 )

  call fill_vector ( mesh, problem, gradexact, degfd=1, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=5 )
  call fill_vector ( mesh, problem, gradexact, degfd=2, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=6 )

! print maximum difference of grad-gradexact to standard output

  print *, maxval( abs(grad%u -gradexact%u) )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh, mesh_plot )
  call delete ( sol, rhsd, solexact )
  call delete ( grad, gradexact )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program poisson37

