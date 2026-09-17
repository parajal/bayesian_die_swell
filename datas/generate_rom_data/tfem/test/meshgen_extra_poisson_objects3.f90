module functions_poisson_m

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

end module functions_poisson_m


! Element routines for the Poisson equation:
!
!    - nabla^2 u = f
!

module poisson_m

  use kind_defs_m
  use mesh_m
  use problem_m, only: problem_t
  use system_m, only: sysvector_t, get_sysvector
  use element_defs_m, only: coefficients_t, oldvectors_t
  use shapefunc_m
  use gauss_m
  use functions_poisson_m

  implicit none


! global parameters for this element

  integer :: funcnr = 0, vfuncnr = 0

  integer :: ndim = 2, nint1D = 3, ninti = 9, nodalp = 9, ndf = 9
  integer :: nodalpb = 3, ndfb = 3
  integer :: ndflb = 2

  real(dp) :: uspecified(2) = [-10,10]

  type(sysvector_t), pointer :: solution

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
        ' either funcnr > 0 or vfuncnr > 0'
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

  subroutine elementc ( mesh, problem, constr, elem, node, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemmat2, elemmatadd, &
    elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


    integer :: object, i, j
    real(dp) :: phi1(1,ndf), phi2(1,ndf), xr(1,2), psi(1,ndflb), x(nodalpb,2)
    real(dp) :: theta(1,ndfb), dtheta(1,ndfb,1), dxdxi(1,2), curvel(1)

    object = problem%constraints(constr)%object

    xr(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)
    call shape_quad_Q2 ( xr, phi1 )
    xr(1,:) = mesh%objects(object)%refcoor2_int(node,:,elem)
    call shape_quad_Q2 ( xr, phi2 )

    call shape_line_P1 ( mesh%objects(object)%xig(node:node,1), psi )
    call shape_line_P2 ( mesh%objects(object)%xig(node:node,1), theta, &
      dtheta(:,:,1) )

    call get_coordinates_object ( mesh, elem, x, object )

    call isoparametric_deformation_curve ( x, dtheta(:,:,1), dxdxi, curvel )

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

      do i = 1, ndflb
        do j = 1, ndf
          elemmat(i,j) = psi(1,i)*phi1(1,j)*curvel(1)*mesh%objects(object)%wg(node)
          elemmat2(i,j) = -psi(1,i)*phi2(1,j)*curvel(1)*mesh%objects(object)%wg(node)
        end do
      end do

    end if

  end subroutine elementc


! sample value in one node of the object

  subroutine sample_node ( mesh, problem, object, nodeobj, coefficients, &
    oldvectors, u )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: object, nodeobj
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: u


    integer :: elgrp, elem
    real(dp) :: phi(1,ndf), xr(1,2), udf(ndf)


    xr(1,:) = mesh%objects(object)%refcoor(nodeobj,:)

    call shape_quad_Q2 ( xr, phi )

    elgrp = mesh%objects(object)%grpelm(nodeobj,1)
    elem = mesh%objects(object)%grpelm(nodeobj,2)

    call get_sysvector ( mesh, problem, solution, elgrp, elem, udf )

    u = dot_product( phi(1,:), udf )

  end subroutine sample_node

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


program meshgen_extra_poisson_objects3

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use problem_m
  use system_m
  use functions_m
  use hsl_ma57_m
  use poisson_m
!  use figplot_m
  use postprocessing_m
  use io_utils_m

  implicit none

  type(mesh_t) :: mesh
  type(input_probdef_t) :: ipd
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(solver_options_ma57_t) :: solver_options

  integer :: i

  call read_mesh ( mesh, filename='mesh2.out' )

! other parts of the mesh

  call fill_mesh_parts ( mesh )

  do i = 1, mesh%nobjects
    print *, 'object = ', i
    print *, 'mesh%objects(:)%refcoor_int'
    print *, mesh%objects(i)%refcoor_int
    print *, 'mesh%objects(:)%grpelm_int'
    print *, mesh%objects(i)%grpelm_int(:,1,:)
    print *, mesh%objects(i)%grpelm_int(:,2,:)
    print *, 'mesh%objects(:)%refcoor2'
    print *, mesh%objects(i)%refcoor2_int
    print *, 'mesh%objects(:)%grpelm2_int'
    print *, mesh%objects(i)%grpelm2_int(:,1,:)
    print *, mesh%objects(i)%grpelm2_int(:,2,:)
  end do

!  call plot_mesh ( plot_options, mesh, 'mesh.fig' )
!  plot_options%objectpointcolor=4
!  plot_options%objectpointsize=0.4
!  call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true. )

  call create_input_probdef ( mesh, ipd )

  ipd%elementdof(1)%a = [1,1,1,1,1,1,1,1,1]
  ipd%elementdof(2)%a = [1,1,1,1,1,1,1,1,1]

  call define_essential &
    ( mesh, ipd, degfd=[1], curve1=1 )

  call define_essential &
    ( mesh, ipd, degfd=[1], curve1=3, curve2=4 )

  call define_essential &
    ( mesh, ipd, degfd=[1], curve1=5, exclude=1 )

  call define_essential &
    ( mesh, ipd, degfd=[1], curve1=7, exclude=2 )

  call define_constraint ( mesh, ipd, object=1, discretization='weak', &
    elementdof=[1,0,1] )

  call problem_definition ( ipd, mesh, problem )

  call create_sysvector ( problem, sol )

  sol%u = 0

  call fill_sysvector ( mesh, problem, sol, &
    degfd=1, curve1=1, func=func, funcnr=3 )

  call fill_sysvector ( mesh, problem, sol, &
    degfd=1, curve1=3, curve2=4, func=func, funcnr=3 )

  call fill_sysvector ( mesh, problem, sol, &
    degfd=1, curve1=5, func=func, funcnr=3 )

  call fill_sysvector ( mesh, problem, sol, &
    degfd=1, curve1=7, func=func, funcnr=3 )

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

  call create_sysvector ( problem, rhsd )

  funcnr = 4

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=poisson_elem )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    elemsub=elementc, addmatvec=.true. )

  call check ( sysmatrix )

  funcnr = 0
  vfuncnr = 1

  call add_boundary_elements ( mesh, problem, rhsd, &
    elemsub=poisson_natboun, curve=6 )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%integer_storage = 1.7
  solver_options%real_storage = 1.7

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

  !plot_options%maxvalue=10
  !plot_options%minvalue=-10
!  call plot_color_contour ( plot_options, mesh, problem, filename='sol.fig', &
!    sysvector=sol )
!  plot_options%objectpointcolor=0
!  plot_options%objectpointsize=0.4
!  call plot_objects ( plot_options, mesh, 'sol.fig', append=.true. )

  print *, maxval(sol%u(problem%degfdperm(1:problem%numnodaldegfd,2)))
  print *, minval(sol%u(problem%degfdperm(1:problem%numnodaldegfd,2)))

  call delete ( problem )
  call delete ( ipd )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )

end program meshgen_extra_poisson_objects3

