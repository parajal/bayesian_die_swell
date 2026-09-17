
! Copyright (C) 2012-2015 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system (e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Routines for performing grid deformation according to:

! Decheng Wan and Stefan Turek, "Fictitious boundary and moving mesh methods
! for the numerical simulation of rigid particulate flows", J. Comp. Phys.,
! 222 (2007) 28-56

! M. Grajewski, M. Kester, and S. Turek., "Mathematical and numerical analysis
! of a robust and efficient grid deformation method in the finite element
! context", SIAM Journal on Scientific Computing, 31 (2009) 1539-1557


module meshgen_grid_deformation_m

  use tfem_elem_m
  use convection_diffusion_set_globals_m

  implicit none

contains


! Element routine to compute the g area function

  subroutine poisson_area_deriv ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors, elemvec, elemwts )

    use poisson_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec, elemwts


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_poisson ( mesh, coefficients, elgrp )

      allocate ( wg(ninti), detF(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    elemvec = sum ( detF * wg )

    elemwts = 1

    if ( last ) then

!     last element in this group

      deallocate ( wg, detF )
      deallocate ( xig, phi, x )
      deallocate ( dphi, F )
      deallocate ( Finv, dphidx )

    end if

  end subroutine poisson_area_deriv


! Element routine to scale a function 1/f of a function f (given in the nodes):
!
!  c_f \int_Omega 1/f dx = Omega
!
! Computed are the elemvec(1) = integral 1/f, elemvec(2) = Omega
! For use with "integrate" in postprocessing.

  subroutine integrate_inverse_griddef_elem ( mesh, problem, elgrp, elem, &
    first, last, coefficients, oldvectors, elemvec )

    use poisson_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec


    if ( first ) then

!     first element

!     set globals

      call set_globals_poisson ( mesh, coefficients, elgrp )

!     allocate arrays

      allocate ( wg(ninti), detF(ninti) )
      allocate ( work(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim) )
      allocate ( u(ndf) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

!   get f in the nodes

    call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )

!   1/f by interpolation of f from the nodes

    work = 1 / matmul ( phi, u )

!   integration of 1/f over the element

    elemvec(1) = sum ( work * detF * wg )

!   volume

    elemvec(2) = sum ( detF * wg )

    if ( last ) then

!     last element

      deallocate ( wg, detF )
      deallocate ( work )
      deallocate ( xig, phi, x )
      deallocate ( dphi, F )
      deallocate ( Finv )
      deallocate ( u )

    end if

  end subroutine integrate_inverse_griddef_elem


! poisson element with additional elemvec 1/f-1/g

  subroutine poisson_griddef_elem ( mesh, problem, elgrp, elem, matrix, &
    vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use poisson_globals_m

    implicit none

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: i, j, ip


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_poisson ( mesh, coefficients, elgrp )

      allocate ( wg(ninti), fg(ninti), detF(ninti), work(ninti) )
      allocate ( xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( dphi(ninti,ndf,ndim), F(ninti,ndim,ndim) )
      allocate ( Finv(ninti,ndim,ndim), dphidx(ninti,ndf,ndim) )
      allocate ( u(nodalp) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    if ( vector ) then

!     1/f

      call get_vector ( mesh, problem, oldvectors%v(2)%p, elgrp, elem, u )
      work = 1/ matmul ( phi, u )

!     1/g

      call get_vector ( mesh, problem, oldvectors%v(1)%p, elgrp, elem, u )
      work = work - 1 / matmul ( phi, u )

      elemvec = matmul ( work * detF * wg, phi )

    end if

    if ( matrix ) then

      do i = 1, ndf
        do j = i, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,i,:) * dphidx(ip,j,:) )
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
      deallocate ( xig, phi, x, work )
      deallocate ( wg, fg, detF )
      deallocate ( u )

    end if

  end subroutine poisson_griddef_elem


! compute the new grid coordinates by time stepping

  subroutine find_grid_coordinates ( mesh, object, coefficients, &
    grid_velocity, fmon, garea, numgridsteps, move_object_points, &
    project_boundary )

!   Upon entry, this contains the initial mesh. It must contain an object
!   that is not used for something else. The object nr is given by the argument
!   object.
!   Upon exit, the coordinates in the object (mesh%objects(object)coor) contains
!   the coordinates of deformed mesh.
    type(mesh_t), intent(inout) :: mesh

!   The object that is being used for obtaining the deformed coordinates. It
!   needs to be defined, but it can be taken as simple as possible.
    integer, intent(in) :: object

!   Coefficients of the Poisson problem
    type(coefficients_t), intent(in) :: coefficients

!   The grid velocity, the scaled monitor function and the scaled area function.
    type(vector_t), intent(in) :: grid_velocity, fmon, garea

!   Number of pseudo time steps on the time interval [0,1]
    integer, intent(in) :: numgridsteps

!   move object points before intersection with the mesh
!   NOTE: this only affects the coordinates before intersecting. After the time
!   step the coordinates are updated correctly. This is for handling periodical
!   domains.
    optional move_object_points
    interface
      subroutine move_object_points ( mesh, object )
        use mesh_m
        implicit none
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: object
      end subroutine move_object_points
    end interface

!   correct for points moved (numerically) out of the domain at the Neumann
!   boundaries after a time step.
    optional project_boundary
    interface
      subroutine project_boundary ( mesh, object )
        use mesh_m
        implicit none
        type(mesh_t), intent(inout) :: mesh
        integer, intent(in) :: object
      end subroutine project_boundary
    end interface

    integer :: step, i
    real(dp) :: grid_deltat
    real(dp), allocatable, dimension(:) :: f_n, g_n, work
    real(dp), allocatable, dimension(:,:) :: eta_n, eta_nm1, eta_nm2, coorn, v_n
    real(dp) :: tn

    allocate ( f_n(mesh%nnodes), g_n(mesh%nnodes), work(mesh%nnodes) )
    allocate ( eta_n(mesh%nnodes,mesh%ndim), eta_nm1(mesh%nnodes,mesh%ndim) )
    allocate ( eta_nm2(mesh%nnodes,mesh%ndim), coorn(mesh%nnodes,mesh%ndim) )
    allocate ( v_n(mesh%nnodes,mesh%ndim) )

!   Time interval [0,1]

    grid_deltat = 1._dp / numgridsteps

    eta_nm1 = 0
    eta_n = 0

    do step = 1, numgridsteps

      tn = grid_deltat * ( step - 1 )

      coorn = mesh%objects(object)%coor

      if ( step == 1 ) then

        v_n = transpose ( &
                  reshape ( grid_velocity%u, [mesh%ndim,mesh%nnodes] ) )
        f_n = fmon%u
        g_n = garea%u

      else

        if ( present( move_object_points ) ) &
          call move_object_points ( mesh, object )

        call find_refcoor_objects ( mesh, object1=object )

        if ( any( mesh%objects(object)%grpelm(:,1) == 0 ) ) then
          write(*,'(/a/)') 'Error find_grid_coordinates: ', &
            ' no reference coordinates.'
          stop
        end if

!       find values of v, f_tilde and g_tilde in the object points

        call find_vn_fn_gn ( mesh, object, coefficients, grid_velocity, &
          fmon, garea, v_n, f_n, g_n )

      end if

      eta_nm2 = eta_nm1
      eta_nm1 = eta_n

      work = tn/f_n + (1-tn)/g_n

      do i = 1, mesh%ndim
        eta_n(:,i) = v_n(:,i) / work
      end do

      if ( step == 1 ) then
!       Explicit Euler at first step
        mesh%objects(object)%coor = coorn + grid_deltat * eta_n
      else if ( step == 2 ) then
!       AB2 at second step
        mesh%objects(object)%coor = coorn &
                         + grid_deltat * ( 3._dp/2._dp*eta_n - 0.5_dp*eta_nm1 )
      else
!       AB3 after second step
        mesh%objects(object)%coor = coorn &
                 + grid_deltat * ( 23._dp/12._dp*eta_n - 4._dp/3._dp*eta_nm1 &
                                  + 5._dp/12._dp*eta_nm2 )
      end if

!     correct for points moving numerically out of the domain

      if ( present( project_boundary ) ) &
        call project_boundary ( mesh, object )

    end do

    deallocate ( f_n, g_n, work )
    deallocate ( eta_n, eta_nm1 )
    deallocate ( eta_nm2, coorn )
    deallocate ( v_n )

  end subroutine find_grid_coordinates


! find the values of v, f, g in the object points

  subroutine find_vn_fn_gn ( mesh, object, coefficients, grid_velocity, &
    fmon, garea, v_n, f_n, g_n )

    use tfem_elem_m
    use poisson_globals_m

    type(mesh_t), intent(inout) :: mesh

!   The object that is being used for obtaining the deformed coordinates. It
!   needs to be defined, but it can be taken as simple as possible.
    integer, intent(in) :: object

!   Coefficients of the Poisson problem
    type(coefficients_t), intent(in) :: coefficients

!   The grid velocity, the scaled inverse monitor function and the inverse
!   scaled area function.
    type(vector_t), intent(in) :: grid_velocity, fmon, garea

!   The interpolated values of the grid velocity, fmon and garea
    real(dp), dimension(:,:), intent(out) :: v_n
    real(dp), dimension(:), intent(out) :: f_n, g_n

    integer :: node, elem, elgrp, elgrp1
    integer :: i

!   set globals

    elgrp1 = mesh%objects(object)%grpelm(1,1)

    call set_globals_poisson ( mesh, coefficients, elgrp1 )

!   allocate arrays

    allocate ( xig(1,ndim), phi(1,ndf), work(ndf) )


    do node = 1, mesh%objects(object)%nnodes

      elgrp = mesh%objects(object)%grpelm(node,1)

!     check change of element group to account for a change in shape

      if ( elgrp /= elgrp1 ) then

!       reset globals

        call set_globals_poisson ( mesh, coefficients, elgrp )

        elgrp1 = elgrp

!       deallocate arrays

        deallocate ( xig, phi, work )

!       allocate arrays

        allocate ( xig(1,ndim), phi(1,ndf), work(ndf) )

      end if

      elem = mesh%objects(object)%grpelm(node,2)
      xig(1,:) = mesh%objects(object)%refcoor(node,:)

      call set_shape_function ( shapefunc, xig, phi )

      do i = 1, mesh%ndim
        work = &
          grid_velocity%u ( mesh%ndim*(mesh%topology(elgrp)%a(:,elem)-1) + i )
        v_n(node:node,i) = matmul ( phi, work )
      end do

      work = fmon%u ( mesh%topology(elgrp)%a(:,elem) )
      f_n(node:node) = matmul ( phi, work )

      work = garea%u ( mesh%topology(elgrp)%a(:,elem) )
      g_n(node:node) = matmul ( phi, work )

    end do

!   deallocate arrays

    deallocate ( xig, phi, work )

  end subroutine find_vn_fn_gn

end module meshgen_grid_deformation_m
