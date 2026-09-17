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
  integer :: ndfl = 4, nodalpl = 4

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
    real(dp) :: phi(1,ndf), xr(1,2), psi(1,ndfl), x(nodalpl,2)
    real(dp) :: dpsi(1,ndfl,2), F(1,2,2), Finv(1,2,2), detF(1)

    object = problem%constraints(constr)%object

    xr(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)

    call shape_quad_Q2 ( xr, phi )
    call shape_quad_Q1 ( mesh%objects(object)%xig(node:node,:), psi, dpsi )

    call get_coordinates_object ( mesh, elem, x, object )

    call isoparametric_deformation ( x, dpsi, F, Finv, detF )

    if ( vector ) then

      elemvec = &
        psi(1,:)*uspecified(object)*detF(1)*mesh%objects(object)%wg(node)

    end if

    if ( matrix ) then

      do i = 1, ndfl
        do j = 1, ndf
          elemmat(i,j) = psi(1,i)*phi(1,j)*detF(1)*mesh%objects(object)%wg(node)
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


program core_poisson_objects

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use problem_m
  use system_m
  use functions_m
  use hsl_ma57_m
  use poisson_m
  !use figplot_m
  use postprocessing_m

  implicit none

  type(meshgen_options_t) :: mesh_options
  !type(plot_options_t) :: plot_options
  type(mesh_t) :: mesh, mesh1
  type(input_probdef_t) :: ipd
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(sample_t) :: sample
  type(solver_options_ma57_t) :: solver_options

  integer :: object, nx, ny

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
  nx = 20
  ny = 20
  mesh_options%nx = nx
  mesh_options%ny = ny
  mesh_options%curved(1) = .true.
  mesh_options%funcnr(1) = 1
  mesh_options%curved(2) = .true.
  mesh_options%funcnr(2) = 2
  mesh_options%curved(3) = .true.
  mesh_options%funcnr(3) = 3
  mesh_options%curved(4) = .true.
  mesh_options%funcnr(4) = 4

  call quadrilateral2d ( mesh, mesh_options, func=cfunc )

  call set_mesh_options ( mesh_options, elshape=5, nx=3, ny=3, &
    ox=1.0_dp, oy=1.5_dp, lx=0.3_dp, ly=0.3_dp )
  call quadrilateral2d ( mesh1, mesh_options )
  call add_to_mesh ( mesh, object='mesh', objectmesh=mesh1, topology=.true., &
     intrule=1, nsubint=5 )
  call delete ( mesh1 )

  call set_mesh_options ( mesh_options, elshape=5, nx=3, ny=3, &
    ox=0.0_dp, oy=1.5_dp, lx=0.3_dp, ly=0.3_dp )
  call quadrilateral2d ( mesh1, mesh_options )
  call add_to_mesh ( mesh, object='mesh', objectmesh=mesh1, topology=.true., &
     intrule=1, nsubint=5 )
  call delete ( mesh1 )

! other parts of the mesh

  call fill_mesh_parts ( mesh )

!  call plot_mesh ( plot_options, mesh, 'mesh.fig' )
!  plot_options%objectpointcolor=4
!  plot_options%objectpointsize=0.4
!  call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true. )

  do object = 1, mesh%nobjects

    write(*,'(a,i0,a)') 'mesh%objects(',object,')%coor'
    write(*,*) mesh%objects(object)%coor
    write(*,'(a,i0,a)') 'mesh%objects(',object,')%grpelm'
    write(*,*) mesh%objects(object)%grpelm
    write(*,'(a,i0,a)') 'mesh%objects(',object,')%refcoor'
    write(*,*) mesh%objects(object)%refcoor
    write(*,'(a,i0,a)') 'mesh%objects(',object,')%grpelm_int'
    write(*,*) mesh%objects(object)%grpelm_int
    write(*,'(a,i0,a)') 'mesh%objects(',object,')%refcoor_int'
    write(*,*) mesh%objects(object)%refcoor_int

  end do

  call create_input_probdef ( mesh, ipd )

  ipd%elementdof(1)%a = [1,1,1,1,1,1,1,1,1]

  call define_essential &
    ( mesh, ipd, degfd=[1], curves=[1,3,4] )
!  call define_essential &
!    ( mesh, ipd, degfd=(/1/), curve1=1 )
!
!  call define_essential &
!    ( mesh, ipd, degfd=(/1/), curve1=3, curve2=4 )

  call define_constraint ( mesh, ipd, object=1, discretization='weak', &
    elementdof=[1,1,1,1] )
  call define_constraint ( mesh, ipd, object=2, discretization='weak', &
    elementdof=[1,1,1,1] )

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
    elemsub=poisson_natboun, curve=2 )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options%integer_storage = 1.7
  solver_options%real_storage = 3.0

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=solver_options )

!  plot_options%maxvalue=10
!  plot_options%minvalue=-10
!  call plot_color_contour ( plot_options, mesh, problem, filename='sol.fig', &
!    sysvector=sol )
!  plot_options%objectpointcolor=0
!  plot_options%objectpointsize=0.4
!  call plot_objects ( plot_options, mesh, 'sol.fig', append=.true. )

  print *, maxval(sol%u(problem%degfdperm(1:problem%numnodaldegfd,2)))
  print *, minval(sol%u(problem%degfdperm(1:problem%numnodaldegfd,2)))

!  open(unit=10,file='out')
!
!  do i = 1, mesh%nnodes
!    write(10,*) mesh%coor(i,1), mesh%coor(i,2), sol%u(problem%degfdperm(i,2))
!    if ( mod(i,2*nx+1) == 0 ) write(10,*)
!  end do
!
!  close(unit=10)

  solution => sol

  call fill_sample ( mesh, problem, sample, ndegfd=1, object=1, &
    elemsub=sample_node )

  print *, 'sample%u'
  print *, sample%u
  print *, 'sample%coor'
  print *, sample%coor

  call delete ( sample )
  call create_sample ( mesh, sample, ndegfd=1, object=2 )
  call fill_sample ( mesh, problem, sample, object=2, elemsub=sample_node )

  print *, 'sample%u'
  print *, sample%u
  print *, 'sample%coor'
  print *, sample%coor

  call delete ( sample )

  call delete ( problem )
  call delete ( ipd )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )

end program core_poisson_objects

