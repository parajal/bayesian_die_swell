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
        func = 1 + sin(pi*x(1))*sin(pi*x(2))
      case(4)
        func = 2*pi**2*sin(pi*x(1))*sin(pi*x(2))
      case(5)
        func = pi*cos(pi*x(1))*sin(pi*x(2))
      case default
        write(*,'(/a/)') 'Error: wrong function number: '
        stop
    end select

  end function func

  function vfunc ( n, nr, x )
    use kind_defs_m
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    select case(nr)
      case(1)
        vfunc = [ pi*cos(pi*x(1))*sin(pi*x(2)),    &
                   pi*sin(pi*x(1))*cos(pi*x(2)) ]
      case default
        write(*,'(/a/)') 'Error: wrong function number: '
        stop
    end select

  end function vfunc

end module poisson_functions_m


! Element routines for the Poisson equation:
!
!    - nabla^2 u = f
!

module poisson_m

  use kind_defs_m
  use mesh_m, only: mesh_t, get_coordinates, get_coordinates_geometry
  use problem_m, only: problem_t
  use element_defs_m, only: coefficients_t, oldvectors_t
  use shapefunc_m
  use gauss_m
  use poisson_functions_m

  implicit none


! global parameters for this element

  integer :: funcnr = 0, vfuncnr = 0

  integer :: ndim = 2, nint1D = 3, ninti = 9, nodalp = 9, ndf = 9
  integer :: nodalpb = 3, ndfb = 3

contains

! Internal element routine for the Poisson Equation

  subroutine poisson_elem ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip
    real(dp), allocatable, dimension(:), save :: wg, fg, detF
    real(dp), allocatable, dimension(:,:), save :: xig, phi, x, tmp, xg
    real(dp), allocatable, dimension(:,:,:), save :: dphi, F, Finv, dphidx


    if ( funcnr <= 0 ) then

      write(*,'(/a/)') &
        'Error in Poisson_elem: funcnr <= 0'
      stop

    end if

    if ( first ) then

!     first element in this group

      allocate ( wg(ninti), fg(ninti), detF(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim), tmp(ndf,ndf) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )

      call Gauss_Legendre_quad ( nint1D, xig, wg )

      call shape_quad_Q2 ( xig, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    call isoparametric_coordinates ( x, phi, xg )

    call shape_derivative ( dphi, Finv, dphidx )

    if ( vector ) then

      do ip = 1, ninti
        fg(ip) = func ( funcnr, xg(ip,:) )
      end do

      do i = 1, ndf
        elemvec(i) = dot_product ( fg * phi(:,i), detF * wg )
      end do

    end if

    if ( matrix ) then

      elemmat = 0
      do ip = 1, ninti
        do i = 1, ndf
          do j = 1, ndf
            tmp(i,j) = dot_product ( dphidx(ip,i,:), dphidx(ip,j,:) )
          end do
        end do
        elemmat = elemmat + tmp * detF(ip) * wg(ip)
      end do

    end if

    if ( last ) then

!     last element in this group

      deallocate ( Finv, dphidx )
      deallocate ( dphi, F )
      deallocate ( xg )
      deallocate ( xig, phi, x, tmp )
      deallocate ( wg, fg, detF )

    end if

  end subroutine poisson_elem


! Boundary element for a natural boundary on a curve for the Poisson equation

  subroutine poisson_natboun ( mesh, problem, curve, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, ip
    real(dp), allocatable, dimension(:), save :: wg, fg, curvel, xig
    real(dp), allocatable, dimension(:,:), save :: phi, x, xg
    real(dp), allocatable, dimension(:,:), save :: dphi, dxdxi, normal


    if ( first ) then

!     first element on this curve

      allocate ( wg(nint1D), fg(nint1D), curvel(nint1D), normal(nint1D,ndim) )
      allocate ( xig(nint1D), phi(nint1D,ndfb), x(nodalpb,ndim) )
      allocate ( xg(nint1D,ndim) )
      allocate ( dphi(nint1D,ndfb), dxdxi(nint1D,ndim) )

      call Gauss_Legendre_line ( nint1D, xig, wg )

      call shape_line_P2 ( xig, phi, dphi )

    end if

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )

    call isoparametric_deformation_curve ( x, dphi, dxdxi, curvel, normal )

    call isoparametric_coordinates ( x, phi, xg )

    if ( funcnr > 0 ) then
      do ip = 1, nint1D
        fg(ip) = func ( funcnr, xg(ip,:) )
      end do
    else if ( vfuncnr > 0 ) then
      do ip = 1, nint1D
        fg(ip) = &
          dot_product ( normal(ip,:), vfunc ( ndim, vfuncnr, xg(ip,:) ) )
      end do
    else
      write(*,'(/2(a/))') 'Error in poisson_natboun:', &
        ' either funcnr > 0 or vfuncnr > 0.'
      stop
    end if

    do i = 1, ndfb
      elemvec(i) = dot_product ( fg * phi(:,i), curvel * wg )
    end do

    if ( last ) then

!     last element on this curve

      deallocate ( wg, fg, curvel, normal )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi, dxdxi )

    end if

  end subroutine poisson_natboun

end module poisson_m
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
        write(*,'(/a/)') 'Error: wrong function number: '
        stop
    end select

  end function cfunc

end module functions_m


program core_poisson

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use problem_m
  use system_m
  use functions_m
  use hsl_ma57_m
  use poisson_m
  use postprocessing_m

  implicit none

  type(meshgen_options_t) :: mesh_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: ipd
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd, solexact
!  type(sample_t) :: sample
  type(solver_options_ma57_t) :: solver_options

!  integer :: i

  mesh_options%elshape = 6
  mesh_options%regionshape = 3
  mesh_options%x2d = &
    reshape ( [ 0.0_dp, 2.0_dp, 2.0_dp, -1.0_dp,    &
                 0.0_dp, 0.0_dp, 2.0_dp,  3.0_dp ], &
              [4,2] )
!  mesh_options%x2d = &
!    reshape ( (/ 0.0_dp, 2.0_dp, 2.0_dp, -1.0_dp,    &
!                 0.0_dp, 0.0_dp, 2.0_dp,  3.0_dp /), &
!              (/4,2/) )
  mesh_options%nx = 20
  mesh_options%ny = 20
  mesh_options%curved(1) = .true.
  mesh_options%funcnr(1) = 1
  mesh_options%curved(2) = .true.
  mesh_options%funcnr(2) = 2
  mesh_options%curved(3) = .true.
  mesh_options%funcnr(3) = 3
  mesh_options%curved(4) = .true.
  mesh_options%funcnr(4) = 4

  call quadrilateral2d ( mesh, mesh_options, func=cfunc )

! other parts of the mesh

  call fill_mesh_parts ( mesh )

  call create_input_probdef ( mesh, ipd )

  ipd%elementdof(1)%a = [1,1,1,1,1,1,1,1,1]

  call define_essential &
    ( mesh, ipd, degfd=[1], curves=[1,3,4] )
!  call define_essential &
!    ( mesh, ipd, degfd=(/1/), curve1=1 )
!
!  call define_essential &
!    ( mesh, ipd, degfd=(/1/), curve1=3, curve2=4 )

  call problem_definition ( ipd, mesh, problem )

  call create_sysvector ( problem, sol )

  sol%u = 0

  call fill_sysvector ( mesh, problem, sol, &
    degfd=1, curves=[1,3,4], func=func, funcnr=3 )
!  call fill_sysvector ( mesh, problem, sol, &
!    degfd=1, curve1=1, func=func, funcnr=3 )
!
!  call fill_sysvector ( mesh, problem, sol, &
!    degfd=1, curve1=3, curve2=4, func=func, funcnr=3 )

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

  call create_sysvector ( problem, rhsd )

  funcnr = 4

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=poisson_elem )

  call check ( sysmatrix )

  funcnr = 0
  vfuncnr = 1

  call add_boundary_elements ( mesh, problem, rhsd, &
    elemsub=poisson_natboun, curve=2 )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%integer_storage = 1.7

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

  call create_sysvector ( problem, solexact )

  call fill_sysvector ( mesh, problem, solexact, &
    degfd=1, node1=1, node2=mesh%nnodes, func=func, funcnr=3 )

  print *, maxval( abs(sol%u -solexact%u) )

!  open(unit=10,file='out')
!
!  do i = 1, mesh%nnodes
!    write(10,*) mesh%coor(i,1), mesh%coor(i,2), sol%u(problem%degfdperm(i,2))
!    if ( mod(i,2*mesh_options%nx+1) == 0 ) write(10,*)
!  end do
!
!  close(unit=10)
!
!  open(unit=10,file='out2')
!
!  call fill_sample ( mesh, problem, sample, ndegfd=1, curve=2, sysvector=sol )
!
!  do i = 1, sample%nnodes
!    write(10,*) sample%coor(i,1), sample%coor(i,2), sample%u(i,1)
!  end do
!
!  close(unit=10)

  call delete ( problem )
  call delete ( ipd )
  call delete ( mesh )
  call delete ( sol, solexact, rhsd )
  call delete ( sysmatrix )
!  call delete ( sample )

end program core_poisson

