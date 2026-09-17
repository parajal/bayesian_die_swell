module subs7_m

  use tfem_elem_m
  use math_defs_m
  use eltree_m
  use convection_diffusion_set_globals_m
  use functions_eltree_m

  implicit none

  save

! array of pointers to an eltree
  type(eltree_p), target, dimension(:,:), allocatable :: eltree_array

! the integration scheme for interface elements in the eltree
!   intrule_ie: the integration rule
!   nintig_ie: number of integration points
!   x_ie, w_ie: points and weights
  integer :: intrule_ie, nintig_ie
  real(dp), allocatable :: w_ie(:), x_ie(:,:)

! the reference coordinates and weights of the Gaussian integration points
  real(dp), allocatable ::xg_ie(:,:), wg_ie(:)

! wn_ie: weight * reference normal in the integration points on the interface
! dan: weight * normal in the integration points on the interface (after mapping)
! da: weight in the integration points on the interface (after mapping)
  real(dp), allocatable :: wn_ie(:,:), dan(:,:), da(:)

contains

! Element routine for use with "integrate" in postprocessing

  subroutine integrate_interface ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors, elemvec )

    use poisson_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    integer :: ip
    type(gauss_t) :: gauss_ie

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_poisson ( mesh, coefficients, elgrp )

      allocate ( x(nodalp,ndim) )

      gauss_ie%globalshape = 'line'
      gauss_ie%intrule = intrule_ie

      call set_ninti ( gauss_ie, nintig_ie )

      allocate ( wg_ie(nintig_ie), xg_ie(nintig_ie,1) )

      call set_Gauss_integration ( gauss_ie, xg_ie, wg_ie )

    end if

    if ( associated(eltree_array(1,elem)%p) ) then

      lelem = elem

!     element crosses the interface

!     allocate arrays

      ninti = number_of_interface_integration_points (eltree_array(1,elem)%p, &
                                                       ninti=nintig_ie )

      allocate ( x_ie(ninti,ndim), w_ie(ninti), wn_ie(ninti,ndim) )

!     composite integration scheme on the interface

      call interface_integration_points ( eltree_array(1,elem)%p, x_ie, w_ie, &
        wn_ie, ninti=nintig_ie, xe=xg_ie, we=wg_ie )

!     allocate arrays

      allocate ( dan(ninti,ndim), da(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim) )
      allocate ( detF(ninti), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( xg(ninti,ndim) )

!     set shape function

      call set_shape_function ( shapefunc, x_ie, phi, dphi )

!     element geometric info

      call get_coordinates ( mesh, elgrp, elem, x )

      call isoparametric_deformation ( x, dphi, F, Finv, detF )

!     find the locations of the integration points in the deformed element

      xg = matmul ( phi, x )

      if ( coorsys == 1 ) then
        detF = 2 * pi * xg(:,2) * detF
      end if

      call shape_derivative ( dphi, Finv, dphidx )

!     transform reference area*normal to actual area*normal: da = J F^{-T} dA

      do ip = 1, ninti
        dan(ip,:) = detF(ip) * matmul ( wn_ie(ip,:), Finv(ip,:,:) )
        da(ip) = sqrt(dot_product(dan(ip,:),dan(ip,:)))
      end do

!     integral over the surface

      elemvec(1) = sum(da)
      elemvec(2) = sum(xg(:,1)**2*da)
      elemvec(3:ndim+2) = sum(dan,dim=1)

!     deallocate arrays

      deallocate ( x_ie, w_ie, wn_ie )
      deallocate ( dan, da )
      deallocate ( phi, dphi )
      deallocate ( detF, F )
      deallocate ( Finv, dphidx )
      deallocate ( xg )

    else

!     normal element inside or outside

      elemvec = 0

    end if

    if ( last ) then

!     last element in this group

      deallocate ( x )
      deallocate ( wg_ie, xg_ie )

    end if

  end subroutine integrate_interface


! define the mesh for each element (interface)

  subroutine userelmesh1 ( mesh, elgrp, elem, indicator, elmesh )

    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: elgrp, elem
    integer, intent(out) :: indicator
    type(mesh_t), intent(inout) :: elmesh

    if ( associated(eltree_array(1,elem)%p) ) then

      lelem = elem

!     add only interface elements
      call eltree_to_interface_mesh( eltree_array(1,elem)%p, elmesh, &
        mapcoor=mapcoor_triangle)

      indicator = 2

    else

!     fully in or out: remove element

      indicator = 1

    end if

  end subroutine userelmesh1

end module subs7_m

