! Copyright (C) 2012-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for stabilization of the Stokes equation
!
!    - div( eta(nabla u+nabla u^T) ) + nabla p = f
!      div u = 0
!
! or other incompressible problems. The stabilization is derived from the
! least-squares term
!
!     - ( p - Pp, p - Pp ) / 2
!
! where Pp is a local projection of the pressure on a lower-order space.

! The stabilization has been described in:
!
! C.R. Dohrmann and P.B. Bochev, "A stabilized finite element method for the
! Stokes problem based on polynomial pressure projections", Int. J. Numer.
! Meth. Fluids 46(2004)183-201
!

module stokes_elements_stabilized_m

  use tfem_elem_m
  use stokes_set_globals_m

  implicit none

contains


! Internal element routine for stabilizing Stokes
! (projection on lower-order space)

  subroutine stokes_elem_stab ( mesh, problem, elgrp, elem, matrix, vector, &
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


    integer :: N, M, info
    real(dp) :: eta_stab


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors, maxvel3D=1 )

    if ( first ) then

!     set shape function and number of degrees of freedom pressure
!     for the projection space

      intpolp1 = coefficients%i(53)

      shapefuncp1%globalshape = globalshape
      shapefuncp1%interpolation = intpolp1
      shapefuncp1%numbering = 'regular'
      shapefuncp1%p = coefficients%i(54)
      shapefuncp1%spec_eval = 'gauss'
!     needed for spectral quads/hexahedra only to call Pp_at_GLL routine:
      shapefuncp1%intrule = shapefunc%p + 1

      call set_ndf ( shapefuncp1, 'stokes_elem_stab', ndf=ndfp1 )

    end if

!   allocate more arrays

    if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      call allocate_arrays

!     set shape function

      call set_shape_function ( shapefuncp1, xig, psi1 )

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call allocate_arrays

!     set shape function

      call set_shape_function ( shapefuncp1, xig, psi1 )

    end if


!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call isoparametric_coordinates ( x(1:ndf,:), phi, xg )

    if ( coorsys == 1 ) then
      detF = 2 * pi * xg(:,2) * detF
    end if

    if ( vector ) then

      elemvec = 0

    end if

    if ( matrix .and. ninti > 0 ) then

      eta_stab = coefficients%r(2)

!     full mass matrix

      do N = 1, ndfp
        do M = N, ndfp
          elemmat(N,M) = sum ( psi(:,N) * psi(:,M) * detF * wg )
          elemmat(M,N) = elemmat(N,M) ! symmetry
        end do
      end do

!     fill lower triangle of sub mass matrix

      do N = 1, ndfp1
        do M = 1, N
          work2(N,M) = sum ( psi1(:,N) * psi1(:,M) * wg * detF )
        end do
      end do

!     Cholesky factorization of sub mass matrix

      call dpotrf ( 'L', ndfp1, work2, ndfp1, info )  ! Cholesky factorization

      if ( info /= 0 ) then
        write(*,'(/a/a/a,i0/2(a,i0)/)') &
          'Error in stokes_element_stab: ', &
          ' Cholesky factorization failed', &
          ' info = ', info, &
          ' elgrp = ', elgrp, ' elem = ', elem
        stop
      end if

!     (psi1_N,psi_M)

      do N = 1, ndfp1
        do M = 1, ndfp
          work4(N,M) = sum ( psi1(:,N) * psi(:,M) * wg * detF )
        end do
      end do

!     transpose and save matrix

      work10 = transpose(work4)

!     M^-1 * (psi1_N,psi_M)

      call dpotrs ( 'L', ndfp1, ndfp, work2, ndfp1, work4, ndfp1, info )

      if ( info /= 0 ) then
        write(*,'(/a/a/3(a,i0/))') 'Error in stokes_element_stab: ', &
          ' Cholesky solve failed.', &
          ' info = ', info, &
          ' elgrp = ', elgrp, &
          ' elem = ', elem
        stop
      end if

!     -(psi_N,psi1_k) * M_km^-1 * (psi1_m,psi_M)

      elemmat = elemmat - matmul ( work10, work4 )

      elemmat = - elemmat / eta_stab

    else if ( matrix ) then

      elemmat = 0

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

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

      allocate ( work(ninti), work1(ninti) )

      allocate ( work2(ndfp1,ndfp1) )
      allocate ( work4(ndfp1,ndfp), work10(ndfp,ndfp1) )

      allocate ( psi1(ninti,ndfp1) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( work, work1 )
      deallocate ( work2 )
      deallocate ( work4, work10 )

      deallocate ( psi1 )

    end subroutine deallocate_arrays

  end subroutine stokes_elem_stab

end module stokes_elements_stabilized_m

