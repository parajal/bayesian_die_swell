
! Copyright (C) 2012-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the L2-projection of vectors defined on one mesh to
! another mesh
!

module projection_elements_m

  use tfem_elem_m

  implicit none

contains


! Project a vector on one mesh to another using L2-projection

  subroutine projection_elem ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use projection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

    if ( coefficients%i(5) == 0 ) then

!     projection on a different mesh

      call projection_elem1 ( mesh, problem, elgrp, elem, matrix, vector, &
        first, last, coefficients, oldvectors, elemmat, elemvec )

    else if ( coefficients%i(5) == 1 ) then

!     projection on the same mesh

      call projection_elem2 ( mesh, problem, elgrp, elem, matrix, vector, &
        first, last, coefficients, oldvectors, elemmat, elemvec )

    else

      write(*,'(/a/a,i0/)') 'Error in projection_elem:', &
      ' incorrect value of coefficients%i(5) = ', coefficients%i(5)
      stop

    end if

  end subroutine projection_elem


! Project a vector to another using L2-projection
! Different mesh.

  subroutine projection_elem1 ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use projection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip


!   set globals, gauss, shapefunctions, ...

    call set_projection_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors, oldvectors%m(1)%p )

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndfs,:), dphis, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndfs,:), phis, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( vector ) then

!     compute reference coordinates in old mesh
!     note: groups not used

      call find_refcoor_points ( oldvectors%m(1)%p, coor=xg, &
        grpelm=grpelm_n, refcoor=xig_n )

      if ( any ( grpelm_n == 0 ) ) then
        write(*,'(/a/a/)') 'Error projection_elem:', &
          ' No reference coordinates found in old mesh.'
        stop
      end if

      if ( any ( oldvectors%m(1)%p%element(grpelm_n(:,1))%globalshape /= &
                 oldvectors%m(1)%p%element(grpelm_n(:,1))%globalshape ) ) then
        write(*,'(/a/a/)') 'Error velocity_projection_elem:', &
          ' Overlapping elements with different globalshape not implemented. '
        stop
      end if

!     shape functions in old elements

      call set_shape_function ( shapefuncv, xig_n, phi_n )

!     get vector in the Gauss points (integration points can be in
!     different elements).

      do ip = 1, ninti

!       get old vector for integration point ip

        call get_vector ( oldvectors%m(1)%p, oldvectors%p(1)%p, &
          oldvectors%v(1)%p, grpelm_n(ip,1), grpelm_n(ip,2), &
          u=un, layer=layer )

        tmp = reshape ( un, [ndfv,ncompv] )

        ung(ip,:) = matmul ( phi_n(ip,:), tmp )

      end do

      do j = 1, ncompv
        do i = 1, ndf
          work2(i,j) = sum ( phi(:,i) * ung(:,j) * detF * wg )
        end do
      end do

      elemvec = reshape ( work2, [ ndf*ncompv ] )

    end if

    if ( matrix ) then

      elemmat = 0

!     mass matrix

      do i = 1, ndf
        do j = 1, ndf
          elemmat(i,j) = sum ( phi(:,i) * phi(:,j) * detF * wg )
        end do
      end do

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_projection_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if


  contains

    subroutine allocate_arrays

      allocate ( un(ncompv*ndfv), ung(ninti,ncompv) )
      allocate ( grpelm_n(ninti,2), xig_n(ninti,ndim) )
      allocate ( phi_n(ninti,ndfv), work2(ndf,ncompv) )
      allocate ( tmp(ndfv,ncompv) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( un, ung )
      deallocate ( grpelm_n, xig_n )
      deallocate ( phi_n, work2 )
      deallocate ( tmp )

    end subroutine deallocate_arrays

  end subroutine projection_elem1


! Project a vector to another using L2-projection
! Same mesh.

  subroutine projection_elem2 ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use projection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j


!   set globals, gauss, shapefunctions, ...

    call set_projection_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors, mesh )

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndfs,:), dphis, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndfs,:), phis, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( vector ) then

!     shape functions in old elements

      call set_shape_function ( shapefuncv, xig, phi_n )

!     get vector in the Gauss points

      call get_vector ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, &
        elgrp, elem, u=un, layer=layer )

      tmp = reshape ( un, [ndfv,ncompv] )

      ung = matmul ( phi_n, tmp )

      do j = 1, ncompv
        do i = 1, ndf
          work2(i,j) = sum ( phi(:,i) * ung(:,j) * detF * wg )
        end do
      end do

      elemvec = reshape ( work2, [ ndf*ncompv ] )

    end if

    if ( matrix ) then

      elemmat = 0

!     mass matrix

      do i = 1, ndf
        do j = 1, ndf
          elemmat(i,j) = sum ( phi(:,i) * phi(:,j) * detF * wg )
        end do
      end do

    end if

!   unset globals, gauss, shapefunctions, ...

    call unset_projection_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if


  contains

    subroutine allocate_arrays

      allocate ( un(ncompv*ndfv), ung(ninti,ncompv) )
      allocate ( phi_n(ninti,ndfv), work2(ndf,ncompv) )
      allocate ( tmp(ndfv,ncompv) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( un, ung )
      deallocate ( phi_n, work2 )
      deallocate ( tmp )

    end subroutine deallocate_arrays

  end subroutine projection_elem2


! compute the projection solution in the nodes for plotting (with derive_vector)

  subroutine deriv_projection ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use projection_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    if ( first ) then

!     first element in this group

      call set_globals_projection ( mesh, coefficients, elgrp )

!     set Gauss integration and shape function

      allocate ( xrnod(nodalp,ndim), phi(nodalp,ndf) )
      allocate ( u(ndf) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefunc, xrnod, phi )

    end if

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      layer=layer )

    elemvec = matmul ( phi, u )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( xrnod, phi )
      deallocate ( u )

    end if

  end subroutine deriv_projection


! set global parameters (internal element)

  subroutine set_globals_projection ( mesh, coefficients, elgrp, meshv )

    use projection_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp
    type(mesh_t), intent(in), optional :: meshv

    integer :: ndimr

!   check size of coefficients

    call check ( coefficients, 'set_globals_projection', ncoefi=100, &
      ncoefr=50, indexarray=[23,37], minimum=[0,0], maximum=[1,2] )

    layer = coefficients%i(38)

    if ( present(meshv) ) then

!     old mesh

      ndim = meshv%element(elgrp)%ndim
      nodalpv = meshv%element(elgrp)%numnod
      globalshapev = meshv%element(elgrp)%globalshape

!     set number of degrees of freedom for the input vector (old mesh)

      intpolv = coefficients%i(1)

      ndimr = mesh%element(elgrp)%ndimr
      if ( any ( intpolv == [ 13 ] ).and. &
                   any(mesh%element(elgrp)%p(:ndimr,2) /= 1 ) ) then
        write(*,'(/2(a/))') 'Error in set_globals_projection:', &
        ' For spectral elements GLL nodal distribution is required. '
        stop
      else if ( any ( intpolv == [ 19, 20 ] ) .and. &
                            any(meshv%element(elgrp)%p(:ndimr,2) /= 0 ) ) then
        write(*,'(/2(a/))') 'Error in set_globals_projection:', &
        ' For high-order elements equidistant nodal distribution is required. '
        stop
      end if

      shapefuncv%globalshape = globalshapev
      shapefuncv%interpolation = intpolv
      shapefuncv%numbering = 'standard'
      shapefuncv%p = coefficients%i(6)
      shapefuncv%spec_eval = 'gauss'

      call set_ndf ( shapefuncv, 'set_globals_projection', ndf=ndfv )

      ncompv = coefficients%i(2)

      if ( mesh%element(elgrp)%ndim /= ndim ) then
        write(*,'(/3(a/))') 'Error in set_globals_projection:', &
          ' Dimension of old and new mesh is different'
        stop
      end if

    end if

!   new mesh

    nodalp = mesh%element(elgrp)%numnod
    if ( ndim == 2 ) then
      coorsys = coefficients%i(23)
    else
      coorsys = 2
    end if
    globalshape = mesh%element(elgrp)%globalshape


!   set shape of the element (new mesh)

    intpols = coefficients%i(3)

    shapefuncs%globalshape = globalshape
    shapefuncs%interpolation = intpols
    shapefuncs%numbering = 'standard'
    shapefuncs%p = coefficients%i(7)
    shapefuncs%spec_eval = 'gauss'

    call set_ndf ( shapefuncs, 'set_globals_projection', ndf=ndfs )

    if ( coefficients%i(34) == 0 .and. ndfs /= nodalp ) then
      write(*,'(/5(a/))') 'Error in set_globals_projection:', &
        ' Number of degrees of freedom of the shape function (ndfs) ', &
        ' is different from the number of nodal points (nodalp) ', &
        ' If you really want ndfs /= nodalp, set coefficients%i(34) == 1,', &
        ' but you are on your own!! '
      stop
    end if


!   set number of degrees of freedom projection sysvector (new mesh)

    intpol = coefficients%i(4)

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(8)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, 'set_globals_projection', ndf=ndf )

!   set integration

    inttype = coefficients%i(40)

    if ( coefficients%i(41) == 1 .or. SET_GAUSS_BY_ORDER ) then
      intrule = set_intrule ( globalshape, inttype, order=coefficients%i(10) )
      intrule2 = set_intrule2 ( globalshape, inttype, order=coefficients%i(11) )
    else
      intrule = coefficients%i(10)
      intrule2 = coefficients%i(11)
    end if
    nsubint = get_coefficient ( coefficients, index=32, default=1 )

    if ( globalshape == 'prism' .and. intrule2 == 0 ) then
      write(*,'(/a/3a/)') 'Error in set_globals_projection:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape
      stop
    end if

    if ( globalshape == 'pyramid' .and. inttype == 2 .and. intrule2 == 0 ) then
      write(*,'(/a/a/3a,i0/)') 'Error in set_globals_projection:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape, ' and inttype = ', inttype
      stop
    end if

    gauss%globalshape = globalshape
    gauss%intrule = intrule
    gauss%intrule2 = intrule2
    gauss%nsubint = nsubint
    gauss%inttype = inttype

    call set_ninti ( gauss, ninti )

  end subroutine set_globals_projection


! Preamble for the projection element

  subroutine set_projection_elem ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, meshv )

    use projection_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    type(mesh_t), intent(in), optional :: meshv


    integer :: i1(3)


    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

!     set globals

      call set_globals_projection ( mesh, coefficients, elgrp, meshv )

!     allocate arrays

      allocate ( x(nodalp,ndim) )

      allocate ( xig(ninti,ndim), wg(ninti) )
      allocate ( phis(ninti,ndfs), dphis(ninti,ndfs,ndim) )
      allocate ( phi(ninti,ndf) )

      allocate ( detF(ninti), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim) )
      allocate ( xg(ninti,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefuncs, xig, phis, dphis )
      call set_shape_function ( shapefunc, xig, phi )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      if ( first ) then

!       first element in this group

!       set globals

        call set_globals_projection ( mesh, coefficients, elgrp, meshv )

!       allocate arrays

        allocate ( x(nodalp,ndim) )

      end if

!     now choose between elvector and user subroutines

      if ( coefficients%i(37) == 1 ) then

!       get element values of intrule, intrule2 and nsubint

        if ( gauss%globalshape == 'prism' .or. &
             gauss%globalshape == 'pyramid' .and. gauss%inttype == 2 ) then

!         intrule2 required
          call get_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, i1=i1 )

          if ( coefficients%i(41) == 1 .or. SET_GAUSS_BY_ORDER ) then
            intrule2 = set_intrule2 ( gauss%globalshape, gauss%inttype, &
                                      order=i1(3) )
          else
            intrule2 = i1(3)
          end if

          gauss%intrule2 = intrule2

        else

          call get_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, i1=i1(1:2) )

        end if

        if ( coefficients%i(41) == 1 .or. SET_GAUSS_BY_ORDER ) then
          intrule = set_intrule ( gauss%globalshape, gauss%inttype, &
                                  order=i1(1) )
        else
          intrule = i1(1)
        end if
        nsubint = i1(2)

        gauss%intrule = intrule
        gauss%nsubint = nsubint

        call set_ninti ( gauss, ninti )

!       allocate arrays

        allocate ( xig(ninti,ndim), wg(ninti) )

!       set Gauss integration

        call set_Gauss_integration ( gauss, xig, wg )

      else if ( coefficients%i(37) == 2 ) then

!       set ninti

        call coefficients%set_ninti_user ( mesh, problem, elgrp, elem, first, &
          last, coefficients, oldvectors )

!       allocate arrays

        allocate ( xig(ninti,ndim), wg(ninti) )

!       set Gauss integration

        call coefficients%set_Gauss_integration_user ( mesh, problem, &
          elgrp, elem, first, last, coefficients, oldvectors )

      end if

!     allocate arrays

      allocate ( phis(ninti,ndfs), dphis(ninti,ndfs,ndim) )
      allocate ( phi(ninti,ndf) )

      allocate ( detF(ninti), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim) )
      allocate ( xg(ninti,ndim) )

!     set shape function

      call set_shape_function ( shapefuncs, xig, phis, dphis )
      call set_shape_function ( shapefunc, xig, phi )

    end if

  end subroutine set_projection_elem


! Unset the preamble for the projection element
! (deallocate arrays allocated in set_... )

  subroutine unset_projection_elem ( last, coefficients )

    use projection_globals_m

    logical, intent(in) :: last
    type(coefficients_t), intent(in) :: coefficients

!   deallocate memory

    if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      deallocate ( x )
      deallocate ( xig, wg )
      deallocate ( phis, dphis )
      deallocate ( phi )

      deallocate ( detF, F )
      deallocate ( Finv )
      deallocate ( xg )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      if ( last ) then

!       last element in this group

        deallocate ( x )

      end if

      deallocate ( xig, wg )
      deallocate ( phis, dphis )
      deallocate ( phi )

      deallocate ( detF, F )
      deallocate ( Finv )
      deallocate ( xg )

    end if

  end subroutine unset_projection_elem

end module projection_elements_m

