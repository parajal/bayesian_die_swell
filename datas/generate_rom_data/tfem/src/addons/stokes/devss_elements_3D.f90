
! Copyright (C) 2005-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the DEVSS and DEVSS-G formulations
!
! DEVSS:
!
!   ... - 2 * alpha * div ( D(u) - E ) = ....
!
!                           D(u) - E   = 0
!
! Note: the last equation is multiplied by 2*alpha in the implementation to
! get a symmetric system matrix.
!
! DEVSS-G:
!
!   ... - alpha * div ( nabla u - G^T ) = ....
!
!                       nabla u - G^T   = 0
!
! Note: by default the last equation is multiplied by alpha in the
! implementation to get a symmetric system matrix. This can be avoided by
! setting a coefficient.
!
! This module contains the 3D elements

module devss_elements_3D_m

  use tfem_elem_m
  use stokes_set_globals_m
  use devss_set_globals_m

  implicit none

contains


! Internal element routine for the DEVSS Equation

  subroutine devss_elem_3D ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: N, M, i1, i2, i3, i4, i5, i6, i7, i8, i9, ip
    real(dp) :: alpha


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_devss_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )


!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of the element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

!   build equations

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

      alpha = coefficients%r(4)

      do N = 1, ndf
        do M = N, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dphidx(ip,M,:) )
          end do
          work1 = work + dphidx(:,N,1) * dphidx(:,M,1)
          Suu(N,M) = alpha * sum ( work1 * detF * wg )
          Suu(M,N) = Suu(N,M)
          work1 = work + dphidx(:,N,2) * dphidx(:,M,2)
          Svv(N,M) = alpha * sum ( work1 * detF * wg )
          Svv(M,N) = Svv(N,M)
          work1 = work + dphidx(:,N,3) * dphidx(:,M,3)
          Sww(N,M) = alpha * sum ( work1 * detF * wg )
          Sww(M,N) = Sww(N,M)
        end do
      end do

      do N = 1, ndf
        do M = 1, ndf
          work = dphidx(:,N,2) * dphidx(:,M,1)
          Suv(N,M) = alpha * sum ( work * detF * wg )
          work = dphidx(:,N,3) * dphidx(:,M,1)
          Suw(N,M) = alpha * sum ( work * detF * wg )
          work = dphidx(:,N,3) * dphidx(:,M,2)
          Svw(N,M) = alpha * sum ( work * detF * wg )
        end do
      end do

      do N = 1, ndf
        do M = 1, ndfe
          Cu11(N,M) = sum ( dphidx(:,N,1) * zeta(:,M) * detF * wg )
          Cu12(N,M) = sum ( dphidx(:,N,2) * zeta(:,M) * detF * wg )
          Cu13(N,M) = sum ( dphidx(:,N,3) * zeta(:,M) * detF * wg )
        end do
      end do

      Cu11 = -2*alpha*Cu11
      Cu12 = -2*alpha*Cu12
      Cu13 = -2*alpha*Cu13

      do N = 1, ndfe
        do M = 1, ndfe
          Dmat(N,M) = sum ( zeta(:,N) * zeta(:,M) * detF * wg )
        end do
      end do

      Dmat = 2*alpha*Dmat

!     pointers in unknowns

      i1 = ndfe        ! E11
      i2 = 2*ndfe      ! E12
      i3 = 3*ndfe      ! E13
      i4 = 4*ndfe      ! E22
      i5 = 5*ndfe      ! E23
      i6 = 6*ndfe      ! E33
      i7 = i6 + ndf    ! u
      i8 = i6 + 2*ndf  ! v
      i9 = i6 + 3*ndf  ! w

!     2*alpha*(F,E)

      elemmat(    1:i6,    1:i6 ) = 0
      elemmat(    1:i1,    1:i1 ) = Dmat
      elemmat( i1+1:i2, i1+1:i2 ) = 2*Dmat
      elemmat( i2+1:i3, i2+1:i3 ) = 2*Dmat
      elemmat( i3+1:i4, i3+1:i4 ) = Dmat
      elemmat( i4+1:i5, i4+1:i5 ) = 2*Dmat
      elemmat( i5+1:i6, i5+1:i6 ) = Dmat

!     -( (nabla v)^T, 2*alpha*E  )

      elemmat( i6+1:i7,    1:i1 ) = Cu11
      elemmat( i6+1:i7, i1+1:i2 ) = Cu12
      elemmat( i6+1:i7, i2+1:i3 ) = Cu13
      elemmat( i6+1:i7, i3+1:i6 ) = 0
      elemmat( i7+1:i8,    1:i1 ) = 0
      elemmat( i7+1:i8, i1+1:i2 ) = Cu11
      elemmat( i7+1:i8, i2+1:i3 ) = 0
      elemmat( i7+1:i8, i3+1:i4 ) = Cu12
      elemmat( i7+1:i8, i4+1:i5 ) = Cu13
      elemmat( i7+1:i8, i5+1:i6 ) = 0
      elemmat( i8+1:i9,    1:i2 ) = 0
      elemmat( i8+1:i9, i2+1:i3 ) = Cu11
      elemmat( i8+1:i9, i3+1:i4 ) = 0
      elemmat( i8+1:i9, i4+1:i5 ) = Cu12
      elemmat( i8+1:i9, i5+1:i6 ) = Cu13

!     transposed of -( (nabla v)^T, 2*alpha*E )

      elemmat(    1:i6, i6+1:i9 ) = transpose( elemmat( i6+1:i9,   1:i6 ) )

!     velocity - velocity part (viscous matrix)
!     ( (nabla v)^T, 2*alpha*D  )

      elemmat( i6+1:i7, i6+1:i7 ) = Suu
      elemmat( i6+1:i7, i7+1:i8 ) = Suv
      elemmat( i6+1:i7, i8+1:i9 ) = Suw
      elemmat( i7+1:i8, i6+1:i7 ) = transpose(Suv)
      elemmat( i7+1:i8, i7+1:i8 ) = Svv
      elemmat( i7+1:i8, i8+1:i9 ) = Svw
      elemmat( i8+1:i9, i6+1:i7 ) = transpose(Suw)
      elemmat( i8+1:i9, i7+1:i8 ) = transpose(Svw)
      elemmat( i8+1:i9, i8+1:i9 ) = Sww

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_devss_elem ( last, coefficients )

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

      allocate ( Cu11(ndf,ndfe), Cu12(ndf,ndfe), Cu13(ndf,ndfe) )
      allocate ( Dmat(ndfe,ndfe) )
      allocate ( Suu(ndf,ndf), Suv(ndf,ndf), Svv(ndf,ndf) )
      allocate ( Suw(ndf,ndf), Svw(ndf,ndf), Sww(ndf,ndf) )
      allocate ( work(ninti), work1(ninti) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( Cu11, Cu12, Cu13 )
      deallocate ( Dmat )
      deallocate ( Suu, Suv, Svv )
      deallocate ( Suw, Svw, Sww )
      deallocate ( work, work1 )

    end subroutine deallocate_arrays

  end subroutine devss_elem_3D


! Internal element routine for the DEVSS-G Equation

  subroutine devssg_elem_3D ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    integer :: k, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, N, M, ip
    real(dp) :: alpha


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_devssg_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )


!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

    end if


!   start of the element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

!   build equations

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix ) then

      alpha = coefficients%r(4)

      do N = 1, ndf
        do M = N, ndf
          do ip = 1, ninti
            work(ip) = sum ( dphidx(ip,N,:) * dphidx(ip,M,:) )
          end do
          Suu(N,M) = alpha * sum ( work * detF * wg )
          Suu(M,N) = Suu(N,M)
        end do
      end do

      do N = 1, ndf
        do M = 1, ndfg
          do k = 1, ndim
            Bmat(N,M,k) = - sum ( dphidx(:,N,k) * zeta(:,M) * detF * wg )
          end do
        end do
      end do

      do N = 1, ndfg
        do M = 1, ndfg
          Dmat(N,M) = sum ( zeta(:,N) * zeta(:,M) * detF * wg )
        end do
      end do

      if ( coefficients%i(30) == 0 ) then
!       symmetric matrix
        Dmat = alpha * Dmat
        Bmat = alpha * Bmat
      end if

!     pointers in unknowns

      i1 = ndfg     ! G11
      i2 = 2*ndfg   ! G12
      i3 = 3*ndfg   ! G13
      i4 = 4*ndfg   ! G21
      i5 = 5*ndfg   ! G22
      i6 = 6*ndfg   ! G23
      i7 = 7*ndfg   ! G31
      i8 = 8*ndfg   ! G32
      i9 = 9*ndfg   ! G33
      i10 = i9 + ndf    ! u
      i11 = i9 + 2*ndf  ! v
      i12 = i9 + 3*ndf  ! w

!     (H,G^T)

      elemmat(    1:i9,    1:i9 ) = 0
      elemmat(    1:i1,    1:i1 ) = Dmat
      elemmat( i1+1:i2, i1+1:i2 ) = Dmat
      elemmat( i2+1:i3, i2+1:i3 ) = Dmat
      elemmat( i3+1:i4, i3+1:i4 ) = Dmat
      elemmat( i4+1:i5, i4+1:i5 ) = Dmat
      elemmat( i5+1:i6, i5+1:i6 ) = Dmat
      elemmat( i6+1:i7, i6+1:i7 ) = Dmat
      elemmat( i7+1:i8, i7+1:i8 ) = Dmat
      elemmat( i8+1:i9, i8+1:i9 ) = Dmat

!     -( (nabla v)^T, G^T  )

      elemmat(  i9+1:i10,    1:i1 ) = Bmat(:,:,1)
      elemmat(  i9+1:i10, i1+1:i2 ) = Bmat(:,:,2)
      elemmat(  i9+1:i10, i2+1:i3 ) = Bmat(:,:,3)
      elemmat(  i9+1:i10, i3+1:i9 ) = 0
      elemmat( i10+1:i11,    1:i3 ) = 0
      elemmat( i10+1:i11, i3+1:i4 ) = Bmat(:,:,1)
      elemmat( i10+1:i11, i4+1:i5 ) = Bmat(:,:,2)
      elemmat( i10+1:i11, i5+1:i6 ) = Bmat(:,:,3)
      elemmat( i10+1:i11, i6+1:i9 ) = 0
      elemmat( i11+1:i12,    1:i6 ) = 0
      elemmat( i11+1:i12, i6+1:i7 ) = Bmat(:,:,1)
      elemmat( i11+1:i12, i7+1:i8 ) = Bmat(:,:,2)
      elemmat( i11+1:i12, i8+1:i9 ) = Bmat(:,:,3)

!     transposed of -( (nabla v)^T, G^T )

      elemmat(    1:i9, i9+1:i12 ) = transpose( elemmat( i9+1:i12,   1:i9 ) )

      if ( coefficients%i(30) == 1 ) then

!       non-symmetric matrix: multiply momentum balance terms with alpha factor
!         -( (nabla v)^T, alpha * G^T  )

        elemmat(  i9+1:i10,    1:i3 ) = alpha * elemmat(  i9+1:i10,    1:i3 )
        elemmat( i10+1:i11, i3+1:i6 ) = alpha * elemmat( i10+1:i11, i3+1:i6 )
        elemmat( i11+1:i12, i6+1:i9 ) = alpha * elemmat( i11+1:i12, i6+1:i9 )

      end if

!     velocity - velocity part (viscous matrix)
!     ( (nabla v)^T, alpha*nabla u  )

      elemmat(   i9+1:i10,  i9+1:i10 ) = Suu
      elemmat(   i9+1:i10, i10+1:i12 ) = 0
      elemmat(  i10+1:i11,  i9+1:i10 ) = 0
      elemmat(  i10+1:i11, i10+1:i11 ) = Suu
      elemmat(  i10+1:i11, i11+1:i12 ) = 0
      elemmat(  i11+1:i12,  i9+1:i11 ) = 0
      elemmat(  i11+1:i12, i11+1:i12 ) = Suu

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_devssg_elem ( last, coefficients )

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

      allocate ( work(ninti) )
      allocate ( Bmat(ndf,ndfg,ndim), Dmat(ndfg,ndfg) )
      allocate ( Suu(ndf,ndf) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work )
      deallocate ( Bmat, Dmat )
      deallocate ( Suu )

    end subroutine deallocate_arrays

  end subroutine devssg_elem_3D

end module devss_elements_3D_m

