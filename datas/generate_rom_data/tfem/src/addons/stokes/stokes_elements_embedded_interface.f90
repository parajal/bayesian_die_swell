
! Copyright (C) 2011-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for embedded interface conditions for the
! Stokes equation using an eltree to divide the domain and/or introduce
! an interface.

module stokes_elements_embedded_interface_m

  use tfem_elem_m
  use stokes_elements_embedded_boundary_m, only: set_stokes_eltree, &
      unset_stokes_eltree
  use stokes_set_globals_m
  use stokes_elements_generic_m, only: set_stokes_shape_function_global

  implicit none

contains


! embedded interface element for the tractions using an eltree
! NOTE: order='DN' in build_system, which is the default if physical quantities
! have been defined.

  subroutine stokes_traction_interface_eltree ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   variables

    integer :: i, j, N, M, ip, orientation, clayer, nd
    real(dp) :: eta


!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( work(ninti), work1(ninti) )
    allocate ( pos1(ndim*ndf+ndfp), pos2(ndim*ndf+ndfp) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
    end do

!   set orientation

    call check ( coefficients, 'stokes_traction_interface_eltree', &
      indexarray=[36,47], minimum=[-1,0], maximum=[1,2] )

    clayer = get_coefficient ( coefficients, index=47, default=1 )

    nd = ndim * ndf

    select case ( clayer )
    case(1)
      pos1 = [ pos, posp + nd ]
      pos2 = [ pos + nd, posp + nd + ndfp ]
    case(2)
      pos2 = [ pos, posp + nd ]
      pos1 = [ pos + nd, posp + nd + ndfp ]
    case default
      call errormsg_case_default ( 'stokes_traction_interface_eltree', &
        'clayer', int_value=clayer )
    end select

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   vector

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     viscous traction term

      eta = coefficients%r(1)

!     contribution to the first side

!     diagonal blocks

      do N = 1, ndf
        do M = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,M,:) * dan(ip,:) )
          end do
          do i = 1, ndim
            work1 = phi(:,N) * ( work + dphidx(:,M,i) * dan(:,i) )
            elemmat( pos1(pos(N,i)), pos1(pos(M,i)) ) = - eta * sum ( work1 )
          end do
        end do
      end do

!     off-diagonal blocks

      do N = 1, ndf
        do M = 1, ndf
          do i = 1, ndim
            do j = 1, ndim
              if ( i == j ) cycle
              work = phi(:,N) * dphidx(:,M,i) * dan(:,j)
              elemmat( pos1(pos(N,i)), pos1(pos(M,j)) ) = - eta * sum ( work )
            end do
          end do
        end do
      end do

!     velocity-pressure part

      do N = 1, ndf
        do M = 1, ndfp
          do i = 1, ndim
            elemmat( pos1(pos(N,i)), pos1(posp(M)) ) = &
                        sum ( phi(:,N) * psi(:,M) * dan(:,i) )
          end do
        end do
      end do

!     pressure-pressure part

      elemmat( pos1(posp), pos1(posp) ) = 0

!     pressure-velocity part

      do i = 1, ndim
        elemmat( pos1(posp), pos1(pos(:,i)) ) = 0
      end do

!     contribution to the second side

      elemmat(pos2,pos1) = - elemmat(pos1,pos1)

!     fill zero to the contribution from the other side (needs separate build)

      elemmat(:,pos2) = 0

    end if

!   deallocate arrays

    deallocate ( work, work1 )
    deallocate ( pos1, pos2 )

    call unset_stokes_eltree

  end subroutine stokes_traction_interface_eltree


! (transposed or Baumann-Oden) embedded interface element using an eltree
! NOTE: order='DN' in build_system, which is the default if physical quantities
! have been defined.

  subroutine stokes_traction_interface_transposed_eltree ( mesh, problem, &
    elgrp, elem, matrix, vector, first, last, coefficients, oldvectors, &
    elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   variables

    integer :: i, j, N, M, ip, clayer, nd
    integer :: orientation, signviscousterm, signpressureterm
    real(dp) :: eta


!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( work(ninti), work1(ninti) )
    allocate ( pos1(ndim*ndf+ndfp), pos2(ndim*ndf+ndfp) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

    if ( any( coefficients%i(62) == [1,2] ) ) then
!     global shape function for pressure
      call set_stokes_shape_function_global ( shapefuncp, x, xg, psi, &
        coefficients%i(62)==1 )
    end if

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
    end do

!   set orientation

    call check ( coefficients, 'stokes_traction_interface_transposed_eltree', &
      indexarray=[36,44,45,47], minimum=[-1,-1,-1,0], maximum=[1,1,1,2] )

    clayer = get_coefficient ( coefficients, index=47, default=1 )

    nd = ndim * ndf

    select case ( clayer )
    case(1)
      pos1 = [ pos, posp + nd ]
      pos2 = [ pos + nd, posp + nd + ndfp ]
    case(2)
      pos2 = [ pos, posp + nd ]
      pos1 = [ pos + nd, posp + nd + ndfp ]
    case default
      call errormsg_case_default ( &
        'stokes_traction_interface_transposed_eltree', &
        'clayer', int_value=clayer )
    end select

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   set sign of viscous term (+1=symmetric (Nitsche),-1=Baumann-Oden)

    signviscousterm = get_coefficient ( coefficients, index=44, default=1 )

!   set sign of pressure term (+1,-1=sign, 0=pressure term is absent)

    signpressureterm = get_coefficient ( coefficients, index=45, default=0 )

!   viscosity

    eta = coefficients%r(1)

!   vector

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     viscous traction term

!     contribution of the first side

!     diagonal blocks

      do N = 1, ndf
        do M = 1, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dan(ip,:) )
          end do
          do i = 1, ndim
            work1 = phi(:,M) * ( work + dphidx(:,N,i) * dan(:,i) )
            elemmat( pos1(pos(N,i)), pos1(pos(M,i)) ) = &
                   - eta * sum ( work1 ) * real ( signviscousterm, kind=dp )
          end do
        end do
      end do

!     off-diagonal blocks

      do N = 1, ndf
        do M = 1, ndf
          do i = 1, ndim
            do j = 1, ndim
              if ( i == j ) cycle
              work = phi(:,M) * dphidx(:,N,j) * dan(:,i)
              elemmat( pos1(pos(N,i)), pos1(pos(M,j)) ) = &
                    - eta * sum ( work ) * real ( signviscousterm, kind=dp )
            end do
          end do
        end do
      end do

!     velocity-pressure part

      do i = 1, ndim
        elemmat( pos1(pos(:,i)), pos1(posp) ) = 0
      end do

!     pressure-pressure part

      elemmat( pos1(posp), pos1(posp) ) = 0

      if ( signpressureterm /= 0 ) then

!       pressure-velocity part

        do N = 1, ndfp
          do M = 1, ndf
            do j = 1, ndim
              elemmat( pos1(posp(N)), pos1(pos(M,j)) ) = &
                     sum ( psi(:,N) * phi(:,M) * dan(:,j) ) &
                        * real ( signpressureterm, kind=dp )
            end do
          end do
        end do

      else

        do j = 1, ndim
          elemmat( pos1(posp), pos1(pos(:,j)) ) = 0
        end do

      end if

!     contribution of the second side

      elemmat(pos1,pos2) = - elemmat(pos1,pos1)

!     fill zero the contribution to the other side (needs separate build)

      elemmat(pos2,:) = 0

    end if

!   deallocate arrays

    deallocate ( work, work1 )
    deallocate ( pos1, pos2 )

    call unset_stokes_eltree

  end subroutine stokes_traction_interface_transposed_eltree


! Embedded interface condition for velocity using an eltree
! NOTE: order='DN' in build_system, which is the default if physical quantities
! have been defined.

  subroutine stokes_embedded_interface_eltree ( mesh, problem, elgrp, elem, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec

!   variables

    integer :: i, j, k, N, M, ip, orientation, info, nd
    real(dp) :: kappa, tr


!   first build scalar mass matrix (volume integration)

!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first=.true., &
      last=.true., coefficients=coefficients, oldvectors=oldvectors )

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    xg =  matmul ( phi, x(1:ndf,:) )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

!   allocate mass matrix

    allocate ( work2(ndf,ndf) )

!   fill lower triangle of mass matrix

    do N = 1, ndf
      do M = 1, N
        work2(N,M) = sum ( phi(:,N) * phi(:,M) * wg * detF )
      end do
    end do

!   Cholesky factorization of mass matrix

    call dpotrf ( 'L', ndf, work2, ndf, info )  ! Cholesky factorization

    if ( info /= 0 ) then
      write(*,'(/a/a/a,i0/2(a,i0)/)') &
        'Error in stokes_embedded_interface_eltree: ', &
        ' Cholesky factorization failed', &
        ' info = ', info, &
        ' elgrp = ', elgrp, ' elem = ', elem
      stop
    end if

!   unset volume integration (only deallocates memory)

    call unset_stokes_elem ( last=.true., coefficients=coefficients )


!   boundary integral part

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( work4(ndf,ndim), work5(ndf,ndf,ndim,ndim) )
    allocate ( work6(ndf,ndf,ndim), work7(ndf,ndf,ndim) )
    allocate ( work9(ndf,ndim,ndim) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
    end do

!   set orientation (this does not have an effect)

    call check ( coefficients, 'stokes_embedded_interface_eltree', &
      indexarray=[36], minimum=[-1], maximum=[1] )

    orientation = get_coefficient ( coefficients, index=36, default=1 )

    dan = real ( orientation, kind=dp ) * dan

!   kappa parameter

    kappa = coefficients%r(11)

!   (phi_N,phi_M n)

    do N = 1, ndf
      do M = N, ndf
        do i = 1, ndim
          work6(N,M,i) = sum ( phi(:,N) * phi(:,M) * dan(:,i) )
          work6(M,N,i) = work6(N,M,i) ! symmetry
        end do
      end do
    end do

!   save matrix

    work7 = work6

!   M^-1 * (phi_N,phi_M n)  or [(phi_N,phi_M n) * M^-1]^T  due to symmetry

    do i = 1, ndim

      call dpotrs ( 'L', ndf, ndf, work2, ndf, work6(:,:,i), ndf, info )

      if ( info /= 0 ) then
        write(*,'(/a/4(a,i0/))') 'Error in stokes_embedded_interface_eltree: ',&
          ' Cholesky solve failed. dim = ', i, &
          ' info = ', info, &
          ' elgrp = ', elgrp, &
          ' elem = ', elem
        stop
      end if

    end do

!   vector

    if ( vector ) then

      elemvec = 0

    end if

!   matrix

    if ( matrix ) then

!     (phi_N n,phi_k) * M_km^-1 * (phi_m,phi_M n)

      do i = 1, ndim
        do j = 1, ndim
          work5(:,:,i,j) = matmul ( work7(:,:,i), work6(:,:,j) )
        end do
      end do

!     fill matrix

      do N = 1, ndf
        do M = 1, ndf

!         first part (full matrix part)

          do i = 1, ndim
            do j = 1, ndim
              elemmat( pos(N,i), pos(M,j) ) = work5(N,M,j,i)
            end do
          end do

!         second part (diagonal block matrix part)

          tr = 0
          do k = 1, ndim
            tr = tr + work5(N,M,k,k)
          end do

          do i = 1, ndim
            elemmat( pos(N,i), pos(M,i) ) = elemmat( pos(N,i), pos(M,i) ) + tr
          end do

        end do
      end do

      nd = size(pos)

      elemmat(1:nd,1:nd) = kappa * elemmat(1:nd,1:nd)
      elemmat(nd+1:2*nd,1:nd) = - elemmat(1:nd,1:nd)
      elemmat(:,nd+1:2*nd) = - elemmat(:,1:nd)

    end if

!   deallocate arrays

    deallocate ( work2 )
    deallocate ( work4, work5 )
    deallocate ( work6, work7 )
    deallocate ( work9 )

    call unset_stokes_eltree

  end subroutine stokes_embedded_interface_eltree


! Nitsche term using an eltree
! NOTE: order='DN' in build_system, which is the default if physical quantities
! have been defined.

  subroutine stokes_Nitsche_interface_eltree ( mesh, problem, elgrp, &
    elem, matrix, vector, first, last, coefficients, oldvectors, &
    elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   variables

    integer :: i, N, M, ip, nd
    real(dp) :: factor, he


!   set globals fluid element

    call set_globals_stokes_vp ( mesh, coefficients, elgrp )

!   set globals eltree

    call set_stokes_eltree ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate arrays

    allocate ( da(ninti) )
    allocate ( work2(ndf,ndf) )

!   compute deformed fluid element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    xg =  matmul ( phi, x )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    call shape_derivative ( dphi, Finv, dphidx )

!   transform reference area*normal to actual area*normal: da = J F^{-T} dA

    do ip = 1, ninti
      dan(ip,:) = detF(ip) * matmul ( wng(ip,:), Finv(ip,:,:) )
      da(ip) = sqrt(dot_product(dan(ip,:),dan(ip,:)))
    end do

!   element length scaling

    he = 1 ! no elementscaling yet

!   Nitsche parameter

    factor = coefficients%r(15)

!   vector

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

!     diagonal blocks

      do N = 1, ndf
        do M = N, ndf
          work2(N,M) = sum ( phi(:,N) * phi(:,M) * da )
          work2(M,N) = work2(N,M) ! symmetry
        end do
      end do

      work2 = factor * work2 / he

      elemmat = 0

      do i = 1, ndim
        elemmat( pos(:,i), pos(:,i) ) = work2
      end do

      nd = size(pos)

      elemmat(1:nd,1:nd) = factor * elemmat(1:nd,1:nd)
      elemmat(nd+1:2*nd,1:nd) = - elemmat(1:nd,1:nd)
      elemmat(:,nd+1:2*nd) = - elemmat(:,1:nd)

    end if

!   deallocate arrays

    deallocate ( da )
    deallocate ( work2 )

    call unset_stokes_eltree

  end subroutine stokes_Nitsche_interface_eltree

end module stokes_elements_embedded_interface_m

