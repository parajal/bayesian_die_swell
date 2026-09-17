! Scalar diffusion problem on a unit square with Dirichlet boundary conditions
! and natural boundary conditions.
! Varying gauss integration rule.
! Varying alpha coefficient using a separate element to fill the Gauss
! point values of alpha.

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
        func = 1.5_dp + cos(pi*x(1))*cos(pi*x(2))
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

module subs3_m

  use convection_diffusion_elements_m
  implicit none

contains


! Internal element routine to fill the Gauss rule.
! This element should be used together with the routine loop_over_elements.
! NOTE: this element can, of course, be merged with fill_alpha_gauss.
! However, for the sake of a clarity the two operations are separated here.

  subroutine fill_gauss ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors


    integer :: i1(2)

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_poisson ( mesh, coefficients, elgrp )

      allocate ( x(nodalp,ndim) )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

!   Fill Gauss rule. This is just a simple example.

    if ( x(1,1) < 0.1_dp ) then

!     intrule
      i1(1) = 1
!     nsubint
      i1(2) = 10

    else

!     intrule
      i1(1) = coefficients%i(10)
!     nsubint
      i1(2) = get_coefficient ( coefficients, index=32, default=1 )

    end if

    call put_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, i1=i1 )

    if ( last ) then

!     last element in this group

      deallocate ( x )

    end if

  end subroutine fill_gauss


! Internal element routine to fill the position dependent coefficient alpha
! in the Gauss points. This element should be used together with the routine
! loop_over_elements. Requires the Gauss rule to be stored element wise first.

  subroutine fill_alpha_gauss ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use convection_diffusion_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors


    integer :: ip, funcnr, i1(2)


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_poisson ( mesh, coefficients, elgrp )

      allocate ( x(nodalp,ndim) )

    end if

!   get element values of intrule and nsubint

    call get_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, i1=i1 )

    intrule = i1(1)
    nsubint = i1(2)

    gauss%intrule = intrule
    gauss%nsubint = nsubint

    call set_ninti ( gauss, ninti )

    allocate ( wg(ninti), xig(ninti,ndim), phi(ninti,ndf) )
    allocate ( xg(ninti,ndim), alphag(ninti) )

!   set Gauss integration and shape function

    call set_Gauss_integration ( gauss, xig, wg )

    call set_shape_function ( shapefunc, xig, phi )

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_coordinates ( x, phi, xg )

!   determine coefficient alpha by a function

    funcnr = coefficients%i(3)
    do ip = 1, ninti
      alphag(ip) = coefficients%func1(1)%p ( funcnr, xg(ip,:) )
    end do

    call put_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, r1=alphag )

    if ( last ) then

!     last element in this group

      deallocate ( x )

    end if

    deallocate ( wg, xig, phi )
    deallocate ( xg, alphag )

  end subroutine fill_alpha_gauss

end module subs3_m

program diffusion3

  use tfem_m
  use hsl_ma57_m
  use convection_diffusion_elements_m
  use poisson_functions_m
  use subs3_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3 integration of quads
    gaussb = 3,         & ! 3 point integration of boundary elements
    nx=20,              & ! number of elements in x
    ny=20,              & ! number of elements in y
    funcnr=8,           & ! function number for the right-hand side
    funcnr_alpha=13,    & ! function number for the alpha coefficient
    vfuncnr=1             ! function number for the flux vector

  real(dp), parameter :: &
    alpha = 1._dp     ! diffusion coefficient


! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(elvector_t), target :: elvector
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(2) = 3  ! compute alpha using Gauss values in separate routine
  coefficients%i(3) = funcnr_alpha
  coefficients%i(4) = 1  ! numerical integration different per element

  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]
  coefficients%i(16) = vfuncnr

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0

  coefficients%func => func
  coefficients%vfunc => vfunc
  coefficients%func1(1)%p => func


! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a(:,1) = 1

  call define_essential ( mesh, input_probdef, curves=[1,3,4] )
!  call define_essential ( mesh, input_probdef, curve1=1 )
!  call define_essential ( mesh, input_probdef, curve1=3, curve2=4 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curves=[1,3,4], func=func, funcnr=7 )
!  call fill_sysvector ( mesh, problem, sol, &
!    curve1=1, func=func, funcnr=7 )
!  call fill_sysvector ( mesh, problem, sol, &
!    curve1=3, curve2=4, func=func, funcnr=7 )


! Create vector defined per element and store integration rule and
! Gauss values of alpha

  call create_elvector ( mesh, elvector, nint1d=1, nreal1d=1 )

  call create ( oldvectors, nelvec=1 )

  oldvectors%e(1)%p => elvector

  call loop_over_elements ( mesh, problem, elemsub=fill_gauss, &
    coefficients=coefficients, oldvectors=oldvectors )
  call loop_over_elements ( mesh, problem, elemsub=fill_alpha_gauss, &
    coefficients=coefficients, oldvectors=oldvectors )


! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=scalar_diffusion_elem, coefficients=coefficients, &
    oldvectors=oldvectors )

  call add_boundary_elements ( mesh, problem, rhsd, curve=2, &
    elemsub=poisson_natboun, coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! print average sol to standard output

  print *, sum( abs(sol%u) ) / sol%n

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( elvector )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program diffusion3

