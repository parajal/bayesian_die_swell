module functions_convecreac_m

  use math_defs_m

  implicit none

contains

  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    real(dp) :: c, dcdx, dcdy

    c=1+sin(pi*x(1))*sin(pi*x(2))
    dcdx=pi*cos(pi*x(1))*sin(pi*x(2))
    dcdy=pi*sin(pi*x(1))*cos(pi*x(2))

    select case(nr)
      case(1)
        func=dcdx+2*dcdy+c
      case(2)
        func=c
      case default
        write(*,'(/a/)') 'Error: wrong function number: '
        stop
    end select

  end function func

  function vfunc ( n, nr, x )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    select case(nr)
      case(1)
        vfunc = [ 1, 2 ]
      case default
        write(*,'(/a/)') 'Error: wrong function number: '
        stop
    end select

  end function vfunc

end module functions_convecreac_m


! Element routines for the steady convection-reaction equation:
!
!    u . nabla c + c = f
!
! using DG with P_2 interpolation.

module convecreac_m

  use kind_defs_m
  use mesh_m
  use problem_m, only: problem_t
  use system_m, only: sysvector_t, get_sysvector
  use element_defs_m, only: coefficients_t, oldvectors_t
  use shapefunc_m
  use gauss_m
  use vector_m
  use functions_convecreac_m

  implicit none


! global parameters for this element

  integer :: funcnr = 1, vfuncnr = 1, funcnrinflow = 2

  integer :: ndim = 2, nint1D = 3, nintb = 3, nodalp = 9, nodalpb = 3
  integer :: ndf = 9, ndfb = 3, nsides = 4, ninti = 9

  type(sysvector_t), pointer :: solution


contains

! Internal element routine for the steady convection-reaction equation

  subroutine convecreac_elem ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: i, j, ip, side
    real(dp), allocatable, dimension(:), save :: wg, fg, detF
    real(dp), allocatable, dimension(:,:), save :: xig, phi
    real(dp), allocatable, dimension(:), save :: xigb, wgb, cinflow
    real(dp), allocatable, dimension(:,:), save :: phib, dphib, curvel
    real(dp), allocatable, dimension(:,:,:), save :: xigs, phis
    real(dp), allocatable, dimension(:,:), save :: x, xg, ug, ugradphi
    real(dp), allocatable, dimension(:,:), save :: dxdxi, ugn, normal
    real(dp), allocatable, dimension(:,:,:), save :: xs, xgs
    real(dp), allocatable, dimension(:,:,:), save :: dphi, F, Finv, dphidx


    if ( first ) then

!     first element in this group

      allocate ( wg(ninti), fg(ninti), detF(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xigb(nintb), wgb(nintb), cinflow(nintb) )
      allocate ( phib(nintb,ndfb), dphib(nintb,ndfb) )
      allocate ( xigs(nintb,ndim,nsides), phis(nintb,ndf,nsides) )
      allocate ( ugn(nintb,nsides), xgs(nintb,ndim,nsides) )
      allocate ( dxdxi(nintb,ndim), curvel(nintb,nsides) )
      allocate ( normal(nintb,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( xs(nodalpb,ndim,nsides) )
      allocate ( ug(ninti,ndim), ugradphi(ninti,ndf) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )

      call Gauss_Legendre_quad ( nint1D, xig, wg )

      call shape_quad_Q2 ( xig, phi, dphi )

      call Gauss_Legendre_line ( nintb, xigb, wgb )

      call shape_line_P2 ( xigb, phib, dphib )

      call spread_to_sides_2D ( xigb, xigs )

      do side = 1, nsides
        call shape_quad_Q2 ( xigs(:,:,side), phis(:,:,side) )
      end do

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call get_coordinates_sides ( mesh, elgrp, elem, xs )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    call isoparametric_coordinates ( x, phi, xg )

    call shape_derivative ( dphi, Finv, dphidx )

    do side = 1, nsides
      call isoparametric_coordinates ( x, phis(:,:,side), xgs(:,:,side) )
      call isoparametric_deformation_curve ( xs(:,:,side), dphib, dxdxi, &
        curvel(:,side), normal )
      do ip = 1, nintb
        ugn(ip,side) = &
          dot_product( normal(ip,:), vfunc ( ndim, vfuncnr, xgs(ip,:,side) ) )
      end do
    end do

    if ( vector ) then

      do ip = 1, ninti
        fg(ip) = func ( funcnr, xg(ip,:) )
      end do

      do i = 1, ndf
        elemvec(i) = sum ( fg * phi(:,i) * detF * wg )
      end do

      do side = 1, nsides
        if ( mesh%sidelem(elgrp)%a(side,elem,1) == 0 ) then
          do ip = 1, nintb
            if ( ugn(ip,side) < 0 ) then
              cinflow(ip) = func ( funcnrinflow, xgs(ip,:,side) )
            else
              cinflow(ip) = 0
            end if
          end do
          do i = 1, ndf
            elemvec(i) = elemvec(i) &
              - sum ( ugn(:,side) * phis(:,i,side) * cinflow * curvel(:,side) &
                      * wgb, mask = ugn(:,side) < 0 )
          end do
        end if
      end do

    end if

    if ( matrix ) then

      do ip = 1, ninti
        ug(ip,:) = vfunc ( ndim, vfuncnr, xg(ip,:) )
      end do

      do ip = 1, ninti
        do i = 1, ndf
          ugradphi(ip,i) = dot_product( ug(ip,:), dphidx(ip,i,:) )
        end do
      end do

      do i = 1, ndf
        do j = 1, ndf
          elemmat(i,j) = &
             sum ( phi(:,i) * ( ugradphi(:,j) + phi(:,j) ) * detF * wg )
        end do
      end do

      do side = 1, nsides
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = elemmat(i,j) &
              - sum ( ugn(:,side) * phis(:,i,side) * phis(:,j,side) * &
                      curvel(:,side) * wgb, mask = ugn(:,side) < 0 )
          end do
        end do
      end do

    end if

    if ( last ) then

!     last element in this group

      deallocate ( wg, fg, detF )
      deallocate ( xig, phi, x )
      deallocate ( xigb, wgb, cinflow  )
      deallocate ( phib, dphib )
      deallocate ( xigs, phis )
      deallocate ( ugn, xgs )
      deallocate ( dxdxi, curvel )
      deallocate ( normal )
      deallocate ( xg )
      deallocate ( xs )
      deallocate ( ug, ugradphi )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )

    end if

  end subroutine convecreac_elem

! Side element routine for the steady convection-reaction equation

  subroutine convecreac_sidelem ( mesh, problem, elgrp, elem, side, &
    first, last, coefficients, oldvectors, elemmat )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem, side
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat

    integer :: i, j, ip
    integer :: sidenr, elgrpnr
    real(dp), allocatable, dimension(:), save :: ugn
    real(dp), allocatable, dimension(:), save :: xigb, wgb, curvel
    real(dp), allocatable, dimension(:,:), save :: phib, dphib
    real(dp), allocatable, dimension(:,:,:), save :: xigs, phis
    real(dp), allocatable, dimension(:,:), save :: x
    real(dp), allocatable, dimension(:,:), save :: dxdxi, normal
    real(dp), allocatable, dimension(:,:), save :: xs, xgs
    logical, allocatable, dimension(:,:), save :: elemmatzeros


    if ( first ) then

!     first element in this group

      allocate ( x(nodalp,ndim) )
      allocate ( xigb(nintb), wgb(nintb) )
      allocate ( phib(nintb,ndfb), dphib(nintb,ndfb) )
      allocate ( xigs(nintb,ndim,nsides), phis(nintb,ndf,nsides) )
      allocate ( ugn(nintb), xgs(nintb,ndim) )
      allocate ( dxdxi(nintb,ndim), curvel(nintb) )
      allocate ( normal(nintb,ndim) )
      allocate ( xs(nodalpb,ndim) )
      allocate ( elemmatzeros(ndf,ndf) )

      call Gauss_Legendre_line ( nintb, xigb, wgb )

      call shape_line_P2 ( xigb, phib, dphib )

      call spread_to_sides_2D ( xigb, xigs )

      do sidenr = 1, nsides
        call shape_quad_Q2 ( xigs(:,:,sidenr), phis(:,:,sidenr) )
      end do

    end if

    elgrpnr = mesh%sidelem(elgrp)%a(side,elem,1)

    if ( elgrpnr /= 0 ) then

      sidenr  = mesh%sidelem(elgrp)%a(side,elem,3)

      call convecreac_sidelem_zeros ( mesh, problem, elgrp, elem, side, &
        first, last, elemmatzeros )

      call get_coordinates ( mesh, elgrp, elem, x )

      call get_coordinates_side ( mesh, elgrp, elem, side, xs )

      call isoparametric_coordinates ( x, phis(:,:,side), xgs )
      call isoparametric_deformation_curve ( xs, dphib, dxdxi, &
        curvel, normal )
      do ip = 1, nintb
        ugn(ip) = &
          dot_product( normal(ip,:), vfunc ( ndim, vfuncnr, xgs(ip,:) ) )
      end do

      do i = 1, ndf
        do j = 1, ndf
          if ( elemmatzeros(i,j) ) cycle
          elemmat(i,j) = &
              sum ( ugn * phis(:,i,side) * phis(nintb:1:-1,j,sidenr) * &
                    curvel * wgb, mask = ugn < 0 )
        end do
      end do

    end if

    if ( last ) then

!     last element in this group

      deallocate ( elemmatzeros )
      deallocate ( x )
      deallocate ( xigb, wgb )
      deallocate ( phib, dphib )
      deallocate ( xigs, phis )
      deallocate ( ugn, xgs )
      deallocate ( dxdxi, curvel )
      deallocate ( normal )
      deallocate ( xs )

    end if

  end subroutine convecreac_sidelem

! element routine for elemmatzeros

  subroutine convecreac_sidelem_zeros ( mesh, problem, elgrp, elem, side, &
    first, last, elemmatzeros )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem, side
    logical, intent(in) :: first, last
    logical, intent(out), dimension(:,:) :: elemmatzeros


    integer :: sidenr, elgrpnr, node1, node2, n1, n2

    if ( first ) then

!     first element in this group

    end if

    elgrpnr = mesh%sidelem(elgrp)%a(side,elem,1)

    if ( elgrpnr /= 0 ) then

      sidenr  = mesh%sidelem(elgrp)%a(side,elem,3)

      elemmatzeros = .true.

      do node1 = 1, nodalpb
        n1 = mesh%element(elgrp)%sidnod(node1,side)
        do node2 = 1, nodalpb
          n2 = mesh%element(elgrp)%sidnod(node2,sidenr)
          elemmatzeros(n1,n2) = .false.
        end do
      end do

    end if

    if ( last ) then

!     last element in this group

    end if

  end subroutine convecreac_sidelem_zeros

! element routine to find continuous nodal values

  subroutine convecreac_vecelem ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts

    call get_sysvector ( mesh, problem, solution, elgrp, elem, elemvec )

    elemwts = 1

  end subroutine convecreac_vecelem

! element routine to fill u with the exact solution using fill_sysvector_elem

  subroutine convecreac_fillexact ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec )
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: ip
    real(dp), allocatable, dimension(:,:), save :: x


    if ( first ) then

!     first element in this group

      allocate ( x(nodalp,ndim) )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    do ip = 1, ndf
      elemvec(ip) = func ( funcnr, x(ip,:) )
    end do

    if ( last ) then

!     last element in this group

      deallocate ( x )

    end if

  end subroutine convecreac_fillexact

end module convecreac_m



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


program core_DG_2d

  use kind_defs_m
  use mesh_m
  use meshgen_m
  use problem_m
  use system_m
  use functions_m
  use hsl_ma41_m
  use convecreac_m
  use postprocessing_m

  implicit none

  type(meshgen_options_t) :: mesh_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: ipd
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd, solexact
  type(vector_t) :: vector, vectorexact
!  type(sample_t) :: sample
  type(solver_options_ma41_t) :: solver_options

  integer :: i, nx, ny

  mesh_options%regionshape = 3
  mesh_options%elshape = 6
!  mesh_options%x2d = &
!    reshape ( (/ 0.0_dp, 2.0_dp, 2.5_dp, 0.0_dp,    &
!                 0.0_dp, 0.0_dp, 1.5_dp, 1.0_dp /), &
!              (/4,2/) )
  mesh_options%x2d = &
    reshape ( [ 0.0_dp, 2.0_dp, 2.0_dp, -1.0_dp,    &
                 0.0_dp, 0.0_dp, 2.0_dp,  3.0_dp ], &
              [4,2] )

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

! other parts of the mesh

  call fill_mesh_parts ( mesh )

  call create_input_probdef ( mesh, ipd, nvec=1 )

  ipd%elementdof(1)%a = [(0,i=1,8),9]
  ipd%vec_elementdof(1)%a = reshape ( [(1,i=1,9)], [9,1] )

  call problem_definition ( ipd, mesh, problem )

  call create_sysvector ( problem, sol )

  sol%u = 0

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )

!  call create_sysmatrix_structure_dg ( sysmatrix, mesh, problem )

  call create_sysmatrix_structure_dg ( sysmatrix, mesh, problem, &
    elemsubzeros=convecreac_sidelem_zeros )

  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

  call create_sysvector ( problem, rhsd )

  funcnr = 1
  funcnrinflow = 2
  vfuncnr = 1

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=convecreac_elem )

!  call build_system_dg ( mesh, problem, sysmatrix, elemsub=convecreac_sidelem, &
!.    addmatrix=.true. )

  call build_system_dg ( mesh, problem, sysmatrix, elemsub=convecreac_sidelem, &
    addmatrix=.true., elemsubzeros=convecreac_sidelem_zeros )

  call check ( sysmatrix )

  solver_options%integer_storage = 1.7
  !solver_options%printlevel = 1

  call solve_system_ma41 ( sysmatrix, rhsd, sol, solver_options=solver_options )

  call create_sysvector ( problem, solexact )

  funcnr = 2

  call fill_sysvector_elem ( mesh, problem, solexact, &
    elemsub=convecreac_fillexact )

  print *, maxval( abs(sol%u -solexact%u) )

  call create_vector ( problem, vector, vec=1 )

  solution => sol

  call derive_vector ( mesh, problem, vector, elemsub=convecreac_vecelem )

  call create_vector ( problem, vectorexact, vec=1 )
  call fill_vector ( mesh, problem, vectorexact, node1=1, node2=mesh%nnodes, &
    func=func, funcnr=2 )
  print *, maxval( abs(vector%u -vectorexact%u) )

!  open(unit=10,file='out')
!
!  do i = 1, mesh%nnodes
!    write(10,*) mesh%coor(i,1), mesh%coor(i,2), vector%u(i)
!    if ( mod(i,2*mesh_options%nx+1) == 0 ) write(10,*)
!  end do
!
!  close(unit=10)
!
!  open(unit=10,file='out2')
!
!  call fill_sample ( mesh, problem, sample, ndegfd=1, curve=2, vector=vector )
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
  call delete ( vector, vectorexact )
!  call delete ( sample )

end program core_DG_2d

