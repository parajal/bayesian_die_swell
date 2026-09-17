
! Element routines for the Laplace equation (with alpha variable):
!
!    - nabla( alpha nabla u ) = 0
!
! where alpha is constant per element.
!

module laplace_elements_ale_m

  use tfem_elem_m
  use convection_diffusion_set_globals_m

  implicit none


contains


! Internal element routine for the Laplace Equation
! (with constant coefficient alpha per element = 1/volume of the element).

  subroutine laplace_elem_ale ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use poisson_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip
    real(dp) :: alpha


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_poisson ( mesh, coefficients, elgrp )

      allocate ( wg(ninti), detF(ninti), work(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    call isoparametric_coordinates ( x, phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     1/area or volume of the element
      alpha = 1 / sum ( detF * wg )

      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = alpha * sum ( dphidx(ip,i,:) * dphidx(ip,j,:) )
          end do
          elemmat(i,j) = sum ( work * detF * wg )
          elemmat(j,i) = elemmat(i,j) ! symmetry
        end do
      end do

    end if

    if ( last ) then

!     last element in this group

      deallocate ( Finv, dphidx )
      deallocate ( dphi, F )
      deallocate ( xg )
      deallocate ( xig, phi, x, work )
      deallocate ( wg, detF )

    end if

  end subroutine laplace_elem_ale

end module laplace_elements_ale_m

