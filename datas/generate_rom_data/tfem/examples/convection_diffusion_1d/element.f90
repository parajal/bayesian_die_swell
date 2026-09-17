
! element routines

module element_m

  use tfem_elem_m
  use functions_m

  implicit none

contains


! internal fem element routine (Galerkin, SU, SUPG)

  subroutine element ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use element_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip
    real(dp) :: beta, Peh, h

    integer :: isource, ibeta, imeth
    real(dp) :: u, alpha, gamma, betac, s

    if ( first ) then

!     first element in this group

!     set globals

      ninti = coefficients%i(4)
      nodalp = mesh%element(elgrp)%numnod
      ndf = 2

!     allocate arrays

      allocate ( xig(ninti), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,1) )
      allocate ( x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( F(ninti,ndim,ndim), Finv(ninti,ndim,ndim) )
      allocate ( fg(ninti), detF(ninti) )
      allocate ( dphidx(ninti,ndf,ndim) )

!     Gauss rule and shape function

      call Gauss_Legendre_line ( ninti, xig, wg )

      call shape_line_P1 ( xig, phi, dphi(:,:,1) )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    call isoparametric_coordinates ( x, phi, xg )

    call shape_derivative ( dphi, Finv, dphidx )

!   set coefficients

    isource = coefficients%i(1)
    ibeta   = coefficients%i(2)
    imeth   = coefficients%i(3)
    u       = coefficients%r(1)
    alpha   = coefficients%r(2)
    gamma   = coefficients%r(3)
    betac   = coefficients%r(4)

    s = sign ( 1._dp, u )
    h = x(nodalp,1)-x(1,1)
    Peh = abs(u)*h/2/alpha

    call upwindparameter ( beta, Peh )

    if ( matrix ) then

      if ( imeth == 0 ) then

!       Galerkin
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( &
                            ( alpha * dphidx(:,i,1) * dphidx(:,j,1) +  &
                              u  * phi(:,i) * dphidx(:,j,1) + &
                              gamma * phi(:,i) * phi(:,j) ) &
                                 * detF * wg )
          end do
        end do

      else if ( imeth == 1 ) then

!       SU
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( &
                            ( alpha * dphidx(:,i,1) * dphidx(:,j,1) +  &
         ( phi(:,i) + s * beta * h / 2 * dphidx(:,i,1) ) * u * dphidx(:,j,1) + &
                              gamma * phi(:,i) * phi(:,j) ) &
                                 * detF * wg )
          end do
        end do

      else if ( imeth == 2 ) then

!       SUPG
        do i = 1, ndf
          do j = 1, ndf
            elemmat(i,j) = sum ( &
                            ( alpha * dphidx(:,i,1) * dphidx(:,j,1) +  &
             ( phi(:,i) + s * beta * h / 2 * dphidx(:,i,1) ) * ( &
                      u * dphidx(:,j,1) + gamma * phi(:,j) ) &
                             ) * detF * wg )
          end do
        end do

      end if

    end if

    if ( vector ) then

      do ip = 1, ninti
        fg(ip) = source ( isource, xg(ip,1) )
      end do

      if ( imeth == 2 ) then

!       SUPG
        do i = 1, ndf
          elemvec(i) = &
            sum ( fg * ( phi(:,i) + &
                             s * beta * h / 2 * dphidx(:,i,1) ) * detF * wg )
        end do

      else

        do i = 1, ndf
          elemvec(i) = sum ( fg * phi(:,i) * detF * wg )
        end do

      end if

    end if

    if ( last ) then

!     last element in this group

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( x, xg )
      deallocate ( F, Finv )
      deallocate ( fg, detF )
      deallocate ( dphidx )

    end if

  contains


! inlined routine for upwind parameter

    subroutine upwindparameter ( beta, Peh )

      real(dp), intent(out) :: beta ! upwind parameter beta
      real(dp), intent(in) :: Peh  ! mesh Peclet number, Peh > 0

      ! beta parameter

      if ( ibeta == 0 ) then
         beta = 0
      else if ( ibeta == 1 ) then
         beta = 1
      else if ( ibeta == 2 ) then
         if ( Peh < 1 ) then
            beta = 0
         else
            beta = 1 - 1/Peh
         end if
      else if ( ibeta == 3 ) then
         beta = (exp(Peh)+exp(-Peh))/(exp(Peh)-exp(-Peh)) - 1/Peh
      else if ( ibeta == 4 ) then
         if ( Peh < 3 ) then
            beta = Peh/3
         else
            beta = 1
         end if
      else if ( ibeta == 5 ) then
         beta = betac
      end if

      beta = sign(beta,u)

    end subroutine upwindparameter

  end subroutine element


! internal spectral element routine

  subroutine element2 ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use element_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip

    integer :: isource
    real(dp) :: u, alpha, gamma

    if ( first ) then

!     first element in this group

!     set globals

      ninti = coefficients%i(4)
      nodalp = mesh%element(elgrp)%numnod
      ndf = ninti

!     allocate arrays

      allocate ( xig(ninti), wg(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,1) )
      allocate ( x(nodalp,ndim), xg(ninti,ndim) )
      allocate ( F(ninti,ndim,ndim), Finv(ninti,ndim,ndim) )
      allocate ( fg(ninti), detF(ninti) )
      allocate ( dphidx(ninti,ndf,ndim) )

!     Gauss rule and shape function

      call Gauss_Legendre_Lobatto_line ( ninti, xig, wg )

      call shape_line_at_GLL ( xig, phi, dphi(:,:,1) )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    call isoparametric_coordinates ( x, phi, xg )

    call shape_derivative ( dphi, Finv, dphidx )

!   set coefficients

    isource = coefficients%i(1)
    u       = coefficients%r(1)
    alpha   = coefficients%r(2)
    gamma   = coefficients%r(3)

    if ( matrix ) then

!     Galerkin
      do i = 1, ndf
        do j = 1, ndf
          elemmat(i,j) = sum ( &
                          ( alpha * dphidx(:,i,1) * dphidx(:,j,1) +  &
                            u  * phi(:,i) * dphidx(:,j,1) + &
                            gamma * phi(:,i) * phi(:,j) ) &
                               * detF * wg )
        end do
      end do

    end if

    if ( vector ) then

      do ip = 1, ninti
        fg(ip) = source ( isource, xg(ip,1) )
      end do

      do i = 1, ndf
        elemvec(i) = sum ( fg * phi(:,i) * detF * wg )
      end do

    end if

    if ( last ) then

!     last element in this group

      deallocate ( xig, wg )
      deallocate ( phi, dphi )
      deallocate ( x, xg )
      deallocate ( F, Finv )
      deallocate ( fg, detF )
      deallocate ( dphidx )

    end if

  end subroutine element2


! boundary element for the natural boundary condition

  subroutine bounelement ( mesh, problem, point, matrix, vector, coefficients, &
    oldvectors, elemmat, elemvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: point
    logical, intent(in) :: matrix, vector
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    real(dp) :: alpha, dcdxL

    alpha   = coefficients%r(2)
    dcdxL   = coefficients%r(5)

    elemvec = alpha * dcdxL

  end subroutine bounelement

end module element_m

