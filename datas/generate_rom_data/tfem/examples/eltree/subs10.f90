module subs10_m

  use tfem_elem_m
  use math_defs_m
  use eltree_m
  use convection_diffusion_set_globals_m
  use functions_eltree_m

  implicit none

  save

! array of pointers to an eltree
  type(eltree_p), target, dimension(:,:), allocatable :: eltree_array

! the integration scheme for volume elements in the eltree
! NOTE: for all elements not containing the interface, subelements of
!       element trees, and submeshelements, the same integration rule is used
!   intrule_ve: the integration rule
!   nintig_ve: number of integration points
!   x_ve, w_ie: points and weights
  integer :: intrule_ve, nintig_ve
  real(dp), allocatable :: w_ve(:), x_ve(:,:)

! the reference coordinates and weights of the Gaussian integration points
  real(dp), allocatable ::xg_ve(:,:), wg_ve(:)

! da: weight in the integration points on the interface (after mapping)
  real(dp), allocatable ::  da(:)

! the levelset function in all the nodes
  real(dp), dimension(:), allocatable :: d

contains

! Element routine for use with "integrate" in postprocessing

  subroutine integrate_volume ( mesh, problem, elgrp, elem, first, &
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
    type(gauss_t) :: gauss_ve

    integer :: nod(mesh%element(1)%numnod)

    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_poisson ( mesh, coefficients, elgrp )

      allocate ( x(nodalp,ndim) )

      gauss_ve%globalshape = 'tetrahedron'
      gauss_ve%intrule = intrule_ve

      call set_ninti ( gauss_ve, nintig_ve )

      allocate ( wg_ve(nintig_ve), xg_ve(nintig_ve,3) )

      call set_Gauss_integration ( gauss_ve, xg_ve, wg_ve )

    end if

!   node numbers for this element
    nod = mesh%topology(1)%a(:,elem)

    if ( associated(eltree_array(1,elem)%p) ) then

      lelem = elem

!     element crosses the interface

!     allocate arrays

      ninti = number_of_integration_points ( eltree_array(1,elem)%p, &
                         lsign=[-1], ninti=nintig_ve, nintis=nintig_ve )

      allocate ( x_ve(ninti,ndim), w_ve(ninti) )

!     composite integration scheme on the interface

      call integration_points_eltree ( eltree_array(1,elem)%p, x_ve, w_ve, &
        lsign=[-1], ninti=nintig_ve, xe=xg_ve, we=wg_ve, nintis=nintig_ve, &
        xs=xg_ve, ws=wg_ve )

!     allocate arrays

      allocate ( da(ninti) )
      allocate ( phi(ninti,ndf), dphi(ninti,ndf,ndim) )
      allocate ( detF(ninti), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( xg(ninti,ndim) )

!     set shape function

      call set_shape_function ( shapefunc, x_ve, phi, dphi )

!     element geometric info

      call get_coordinates ( mesh, elgrp, elem, x )

      call isoparametric_deformation ( x, dphi, F, Finv, detF )

      xg = matmul ( phi, x )

      if ( coorsys == 1 ) then
        detF = 2 * pi * xg(:,2) * detF
      end if

      call shape_derivative ( dphi, Finv, dphidx )

!     transform reference area to actual area
      do ip = 1, ninti
         da(ip) = detF(ip) * w_ve(ip)
      end do

!     integral over the surface
      elemvec(1) = sum(da)
      elemvec(2) = sum(xg(:,1)**2*da)

!     deallocate arrays

      deallocate ( w_ve, x_ve )
      deallocate ( da )
      deallocate ( phi, dphi )
      deallocate ( detF, F )
      deallocate ( Finv, dphidx )
      deallocate ( xg )

    else if ( all ( d(nod) < 0 ) ) then

!     element inside the interface

      allocate ( phi(nintig_ve,ndf), dphi(nintig_ve,ndf,ndim) )
      allocate ( detF(nintig_ve), F(nintig_ve,ndim,ndim) )
      allocate ( Finv(nintig_ve,ndim,ndim), dphidx(nintig_ve,ndf,ndim) )
      allocate ( xg(nintig_ve,ndim) )
      allocate ( da(nintig_ve) )
      allocate ( x_ve(nintig_ve,ndim), w_ve(nintig_ve) )


!     element inside or outside

      call set_shape_function ( shapefunc, xg_ve, phi, dphi )

!     element geometric info

      call get_coordinates ( mesh, elgrp, elem, x )

      call isoparametric_deformation ( x, dphi, F, Finv, detF )

!     find the locations of the integration points in the deformed element

      xg = matmul ( phi, x )

      if ( coorsys == 1 ) then
        detF = 2 * pi * x_ve(:,2) * detF
      end if

      call shape_derivative ( dphi, Finv, dphidx )

!     transform reference area to actual area*normal: da = J F^{-T} dA
      do ip = 1, nintig_ve
         da(ip) = detF(ip) * wg_ve(ip)
      end do

!     integral over the surface
      elemvec(1) = sum(da)
      elemvec(2) = sum(xg(:,1)**2*da)

      deallocate ( da, xg)
      deallocate ( phi, dphi )
      deallocate ( detF, F )
      deallocate ( Finv, dphidx )
      deallocate ( x_ve, w_ve )

    else
      elemvec = 0

    end if

    if ( last ) then

!     last element in this group
      deallocate ( wg_ve, xg_ve )
      deallocate ( x )

    end if

  end subroutine integrate_volume


! define the mesh for each element (interface)

  subroutine userelmesh1 ( mesh, elgrp, elem, indicator, elmesh )

    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: elgrp, elem
    integer, intent(out) :: indicator
    type(mesh_t), intent(inout) :: elmesh

    if ( associated(eltree_array(1,elem)%p) ) then

      lelem = elem

!     add only interface elements
      call eltree_to_interface_mesh ( eltree_array(1,elem)%p, elmesh, &
        mapcoor=mapcoor_tet )

      indicator = 2

    else

!     fully in or out: remove element

      indicator = 1

    end if

  end subroutine userelmesh1

end module subs10_m

