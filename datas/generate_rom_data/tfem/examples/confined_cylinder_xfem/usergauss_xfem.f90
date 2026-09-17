
! module with additional routines for the cylinder problem with xfem

module usergauss_xfem_m

  use tfem_elem_m
  use functions_xfem_m

  implicit none
  save

! the levelset function in all the nodes
  real(dp), dimension(:), allocatable :: d

! integration area near interface elements
  real(dp), dimension(:), allocatable :: vol

! the standard integration scheme:
!   ninti_s: number of integration points
!   xs, ws: points and weights
! the integration scheme on subelements in the eltree:
!   intrule_e: the integration rule
!   ninti_e: number of integration points
!   xe, we: points and weights
! the integration scheme on subelements in the submesh (triangles):
!   intrule_sub: the integration rule
!   intrule_sub_small: the integration rule for small integration areas
!   ninti_sub: number of integration points
!   ninti_sub_small: number of integration points for small integration areas
!   xsub, wsub: points and weights
!   xsub_small, wsub_small: points and weights for small integration areas

  integer :: ninti_s, intrule_e, ninti_e, intrule_sub, intrule_sub_small, &
             ninti_sub, ninti_sub_small

  real(dp), allocatable :: ws(:), xs(:,:), we(:), xe(:,:), wsub(:), xsub(:,:), &
                           wsub_small(:), xsub_small(:,:)

! integration options for setting the composite Gauss integration on the eltree
  type(integration_options_t) :: intopt

! some logicals
  logical :: printnumel = .false., printint = .false.

! array of pointers to an eltree (for userelmesh only)

  type(eltree_p), dimension(:,:), pointer :: ea

! treshold for small volumes

  real(dp) :: vol_small = 1e-2_dp


contains


! user routine for setting the number of integration points

  subroutine set_ninti_user ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors

    integer :: nod(nodalp)
    type(gauss_t) :: gausse, gausssub, gausssub_small

    if ( first ) then

!     first element in this group

!     Gauss rule for the standard case

      ninti_s = ninti

      allocate ( ws(ninti_s), xs(ninti_s,2) )

      call set_Gauss_integration ( gauss, xs, ws )

!     Gauss rule for the subdivided elements in the eltree

      gausse = gauss
      gausse%intrule = intrule_e

      call set_ninti ( gausse, ninti_e )

      allocate ( we(ninti_e), xe(ninti_e,2) )

      call set_Gauss_integration ( gausse, xe, we )

!     Gauss rule for the subdivided elements in the submesh

      gausssub%globalshape = 'triangle'
      gausssub%intrule = intrule_sub
      gausssub%inttype = 3  ! use numerical Gauss values

      call set_ninti ( gausssub, ninti_sub )

      allocate ( wsub(ninti_sub), xsub(ninti_sub,2) )

      call set_Gauss_integration ( gausssub, xsub, wsub )

!     Gauss rule for the subdivided elements in the submesh (small volumes)

      gausssub_small%globalshape = 'triangle'
      gausssub_small%intrule = intrule_sub_small
      gausssub_small%inttype = 3  ! use numerical Gauss values

      call set_ninti ( gausssub_small, ninti_sub_small )

      allocate ( wsub_small(ninti_sub_small), xsub_small(ninti_sub_small,2) )

      call set_Gauss_integration ( gausssub_small, xsub_small, wsub_small )

    end if

!   set eltree pointer (alias)

    eltree1 => oldvectors%ea(1)%p(elgrp,elem)%p

!   Set ninti

    nod = mesh%topology(elgrp)%a(:,elem)

    if ( associated(eltree1) ) then

!     element crosses the interface

      if ( printnumel ) then
        print *, 'number of subelements outside = ', &
          number_of_subelements ( eltree1, lsign=[1] )
      end if

!     composite integration scheme for levelset >= 0

      if ( vol(elem) > vol_small ) then
        ninti = number_of_integration_points ( eltree1, lsign=[1],&
          ninti=ninti_e, nintis=ninti_sub, integration_options=intopt )
      else
        ninti = number_of_integration_points ( eltree1, lsign=[1],&
          ninti=ninti_e, nintis=ninti_sub_small, integration_options=intopt )
      end if

      if ( printint ) then
        print *, 'number of integration points = ', ninti
      end if

    else if ( all ( d(nod) >= 0._dp ) ) then

!     fully out

      ninti = ninti_s

    else

!     fully in or removed element due to small integration area

      ninti = 0

    end if

  end subroutine set_ninti_user


! set the composite gauss integration rule with a user subroutine

  subroutine set_Gauss_integration_user ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors

    integer :: nod(nodalp)

!   set eltree pointer

    eltree1 => oldvectors%ea(1)%p(elgrp,elem)%p

    nod = mesh%topology(elgrp)%a(:,elem)

    if ( associated(eltree1) ) then

!     element crosses the interface

!     composite integration scheme for levelset >= 0

      if ( vol(elem) > vol_small ) then
        call integration_points ( eltree1, xig, wg, lsign=[1], &
          ninti=ninti_e, xe=xe, we=we, nintis=ninti_sub, xs=xsub, ws=wsub, &
          integration_options=intopt )
      else
        call integration_points ( eltree1, xig, wg, lsign=[1], &
          ninti=ninti_e, xe=xe, we=we, nintis=ninti_sub_small, &
          xs=xsub_small, ws=wsub_small, integration_options=intopt )
      end if

    else if ( all ( d(nod) >= 0._dp ) ) then

!     fully out: standard integration

      call set_Gauss_integration ( gauss, xig, wg )

    else

!     fully in or removed element due to small integration area
!     generate fully zero elemmat and elemvec

      xig = 0; wg = 0

    end if

!   delete some allocated memory

    if ( last ) then

      deallocate ( we, xe, ws, xs, wsub, xsub, wsub_small, xsub_small )

    end if

  end subroutine set_Gauss_integration_user


! define the mesh for each element

  subroutine userelmesh ( mesh, elgrp, elem, indicator, elmesh )

    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: elgrp, elem
    integer, intent(out) :: indicator
    type(mesh_t), intent(inout) :: elmesh

    integer :: nod(mesh%element(elgrp)%numnod)
    type(eltree_t), pointer :: eltree


!   set eltree pointer for convenience

    eltree => ea(elgrp,elem)%p

    nod = mesh%topology(elgrp)%a(:,elem)

    if ( associated(eltree) ) then

!     element crosses the interface

      lelgrp = elgrp
      lelem = elem

!     add only elements in the outside region

      call eltree_to_mesh ( eltree, elmesh, lsign=[1], mapcoor=mapcoor )

      indicator = 2 ! replace element with submesh

    else if ( all ( d(nod) >= 0._dp ) ) then

!     fully out: leave element as it is.

      indicator = 0

    else

!     fully in or removed element due to small integration area

      indicator = 1

    end if

  end subroutine userelmesh

end module usergauss_xfem_m
