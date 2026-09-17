
! Element routines for the steady convection-reaction equation:
!
!    u . nabla c + c = f
!
! using DG with P_2 interpolation.

module convecreac_elements_m

  use tfem_elem_m
  use convecreac_functions_m

  implicit none


! global parameters for this element

  integer :: funcnr = 0, vfuncnr = 0, funcnrinflow = 0

  integer :: ndim = 2, nint1D = 3, nintb = 3, nodalp = 9, nodalpb = 3
  integer :: ndf = 9, ndfb = 3, nsides = 4, ninti = 9

  type(sysvector_t), save :: solution


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


    if ( funcnr <= 0 ) then
      write(*,'(/a/)') &
        'Error in convecreac_elem: funcnr <= 0'
      stop
    end if

    if ( vfuncnr <= 0 ) then
      write(*,'(/a/)') &
        'Error in convecreac_elem: vfuncnr <= 0'
      stop
    end if

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

    if ( vfuncnr <= 0 ) then
      write(*,'(/a/)') &
        'Error in convecreac_sideelem: vfuncnr <= 0'
      stop
    end if

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

end module convecreac_elements_m
