
! Copyright (C) 2009-2020 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Element routines for the constitutive equation based on SUPG.
! Fluctuating hydrodynamics.
!

module viscoelastic_elements_fh_m

  use viscoelastic_elements_m
  use generalized_stokes_elements_m

  implicit none

  save

contains


! Internal element routine for the additional term in the constitutive
! equations with fluctuating hydrodynamics. SUPG.

  subroutine ce_supg_fh_elem ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


!   select depending on integration scheme

    select case ( coefficients%i(22) )
      case(1,3) ! first-order
        call ce_supg_fh_elem1 ( mesh, problem, elgrp, elem, matrix, vector, &
          first, last, coefficients, oldvectors, elemmat, elemvec )
      case default
        write(*,'(/a/a,i0/)') 'Error in ce_supg_fh_elem:', &
        ' incorrect time integration scheme = ', coefficients%i(22)
        stop
    end select

  end subroutine ce_supg_fh_elem


! Internal element routine for the additional term in the constitutive
! equations with fluctuating hydrodynamics. SUPG.
! First-order semi-implicit time-integration

  subroutine ce_supg_fh_elem1 ( mesh, problem, elgrp, elem, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemvec )

    use generalized_stokes_globals_m
    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat
    real(dp), intent(out), dimension(:) :: elemvec


    real(dp), parameter :: tol = 1.e-10_dp ! cut-off for equal eigenvalues
    integer :: i, j, ip, m, numstress
    real(dp) :: beta, deltat, kT, esize
    real(dp) :: a1, a2, a3, s1, s2, s3, c1, c2, c3, adiag(3), cdiag(3)


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    call set_globals_generalized_stokes ( mesh, coefficients, elgrp )

    call set_viscoelastic_elem ( mesh, problem, elgrp, elem, first, &
      last, coefficients, oldvectors )

    if ( coorsys == 1 ) then
      write(*,'(a/)') &
        ' Error ce_supg_fh_elem1: axisymmetric not applicable '
      stop
    end if

    if ( ndim == 2 .and. ncompc /=3 ) then
      write(*,'(a/)') &
        ' Error ce_supg_fh_elem1: not applicable for models with more ', &
        ' than three conformation tensor components in 2D.'
      stop
    end if


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

    call shape_derivative ( dtheta, Finv, dthetadx )

!   get velocity vector

    call get_sysvector ( mesh, oldvectors%p(1)%p, oldvectors%s(1)%p, &
      elgrp, elem, u, physq=[physqvel], layer=layer )

    tmp = reshape ( u, [ndf,ndim] )

    uvecn = matmul ( phi, tmp )

!   un.grad operator

    do ip = 1, ninti
      ungradtheta(ip,:) = matmul ( dthetadx(ip,:,:), uvecn(ip,:) )
    end do

!   get conformation tensor at previous time step

    do m = mode1, mode2

      call get_sysvector ( mesh, problem, oldvectors%s2(1)%p(1,m), &
        elgrp, elem, cn(:,1), posu=posc, layer=layer )

      do j = 2, ncompc
        cn(:,j) = oldvectors%s2(1)%p(j,m)%u(posc)
      end do

      cng(:,:,m) = matmul ( theta, cn )

    end do

!   upwinding parameter and other factors

    beta = coefficients%r(9)

    if ( htype == 3 ) then
      esize = sum ( detF * wg )
    else
      esize = 1
    end if

!   compute h/U

    call supg_hU ( ndim, globalshape, hoverU, &
      htype=htype, Uscaling=Uscaling, Uglobal=coefficients%r(10), uvec=uvecn, &
      Finv=Finv, x=x, esize=esize )

    tau = beta * hoverU / 2  ! upwinding parameter

    deltat = coefficients%r(8)

    if ( matrix ) then

      write(*,'(a/a/)') &
        ' Error ce_supg_fh_elem1: element matrix not applicable. ', &
        '   call build_system with buildmatrix=.false. '
      stop

    end if

    if ( vector ) then

!     viscoelastic rhs for fluctuating hydrodynamics

!     fluctuating term

      do m = mode1, mode2

!       generate random numbers for stress mode m

        call generate_brownf ( brownf, coefficients )

        do ip = 1, ninti

!         compute principal values of conformation tensor

          if ( ndim == 2 ) then

!           2D

!           compute principal values/directions

            call eig2x2 ( cng(ip,1:3,m), eigvalue, eigv )

            if ( logc ) then

!             log-conformation formulation

              s1 = eigvalue(1)
              s2 = eigvalue(2)
              s3 = 0

              c1 = exp(s1)
              c2 = exp(s2)
              c3 = 1

            else

!             standard formulation

              c1 = eigvalue(1)
              c2 = eigvalue(2)
              c3 = 1

            end if

          else if ( ndim == 3 ) then

!           3D

!           compute principal values/directions

            call eig3x3 ( cng(ip,1:6,m), eigvalue, eigv )

            if ( logc ) then

!             log-conformation formulation

              s1 = eigvalue(1)
              s2 = eigvalue(2)
              s3 = eigvalue(3)

              c1 = exp(s1)
              c2 = exp(s2)
              c3 = exp(s3)

            else

!             standard formulation

              c1 = eigvalue(1)
              c2 = eigvalue(2)
              c3 = eigvalue(3)

            end if

          end if

          cdiag = [ c1, c2, c3 ]

          call principal_components_of_a ( vemodel, m, cdiag, adiag )

          a1 = adiag(1)
          a2 = adiag(2)
          a3 = adiag(3)

!         fluctuating flux tensor dDi

          if ( ndim == 2 ) then

!           2D

            dDi(1,1) = sqrt( a1 ) * brownf(ip,1)
            dDi(1,2) = sqrt( a1 + a2 ) / 2 * brownf(ip,2)
            dDi(2,1) = dDi(1,2)
            dDi(2,2) = sqrt( a2 ) * brownf(ip,3)

          else if ( ndim == 3 ) then

!           3D

            dDi(1,1) = sqrt( a1 ) * brownf(ip,1)
            dDi(1,2) = sqrt( a1 + a2 ) / 2 * brownf(ip,2)
            dDi(1,3) = sqrt( a1 + a3 ) / 2 * brownf(ip,3)
            dDi(2,1) = dDi(1,2)
            dDi(2,2) = sqrt( a2 ) * brownf(ip,4)
            dDi(2,3) = sqrt( a2 + a3 ) / 2 * brownf(ip,5)
            dDi(3,1) = dDi(1,3)
            dDi(3,2) = dDi(2,3)
            dDi(3,3) = sqrt( a3 ) * brownf(ip,6)

          end if

!         final fluctuating term

          if ( ndim == 2 ) then

!           2D

            if ( logc ) then

!             log-conformation formulation

              work2(1,1) = - 2 * dDi(1,1)
              work2(2,2) = - 2 * dDi(2,2)

              if ( abs(c1-c2) > tol ) then
                work2(1,2) = - (s1 - s2)*(c1+c2)*dDi(1,2)/(c1-c2)
              else
                work2(1,2) = - 2 * dDi(1,2)
              end if
              work2(2,1) = work2(1,2)

            else

!             standard formulation

              work2(1,1) = - 2 * c1 * dDi(1,1)
              work2(1,2) = - ( c1 + c2 ) * dDi(1,2)
              work2(2,1) = work2(1,2)
              work2(2,2) = - 2 * c2 * dDi(2,2)

            end if

!           transform back to global system

            work2 = matmul(eigv,matmul(work2,transpose(eigv)))

            fng(ip,1,m) = work2(1,1)
            fng(ip,2,m) = work2(1,2)
            fng(ip,3,m) = work2(2,2)

          else if ( ndim == 3 ) then

!           3D

            if ( logc ) then

!             log-conformation formulation

              work2(1,1) = - 2 * dDi(1,1)
              work2(2,2) = - 2 * dDi(2,2)
              work2(3,3) = - 2 * dDi(3,3)

              if ( abs(c1-c2) > tol ) then
                work2(1,2) = - (s1 - s2)*(c1+c2)*dDi(1,2)/(c1-c2)
              else
                work2(1,2) = - 2 * dDi(1,2)
              end if
              work2(2,1) = work2(1,2)

              if ( abs(c2-c3) > tol ) then
                work2(2,3) = - (s2 - s3)*(c2+c3)*dDi(1,2)/(c2-c3)
              else
                work2(2,3) = - 2 * dDi(2,3)
              end if
              work2(3,2) = work2(2,3)

              if ( abs(c3-c1) > tol ) then
                work2(3,1) = - (s3 - s1)*(c3+c1)*dDi(1,3)/(c3-c1)
              else
                work2(3,1) = - 2 * dDi(3,1)
              end if
              work2(1,3) = work2(3,1)

            else

!             standard formulation

              work2(1,1) = - 2 * c1 * dDi(1,1)
              work2(1,2) = - ( c1 + c2 ) * dDi(1,2)
              work2(1,3) = - ( c1 + c3 ) * dDi(1,3)
              work2(2,1) = work2(1,2)
              work2(2,2) = - 2 * c2 * dDi(2,2)
              work2(2,3) = - ( c2 + c3 ) * dDi(2,3)
              work2(3,1) = work2(1,3)
              work2(3,2) = work2(2,3)
              work2(3,3) = - 2 * c3 * dDi(3,3)

            end if

!           transform back to global system

            work2 = matmul(eigv,matmul(work2,transpose(eigv)))

            fng(ip,1,m) = work2(1,1)
            fng(ip,2,m) = work2(1,2)
            fng(ip,3,m) = work2(1,3)
            fng(ip,4,m) = work2(2,2)
            fng(ip,5,m) = work2(2,3)
            fng(ip,6,m) = work2(3,3)

          end if

        end do

      end do

      kT = coefficients%r(208)

      work = sqrt ( 2 * kT * detF * wg / deltat )

!     semi-implicit or explicit time integration

      do m = mode1, mode2
        do j = 1, ncompc
          do i = 1, ndfc
            work6(i,j,m) = &
              sum ( ( theta(:,i) + tau * ungradtheta(:,i) ) * &
                                fng(:,j,m) * work )
          end do
        end do
      end do

      elemvec = reshape ( work6, [ ndfc * ncompc * ( mode2 - mode1 + 1 ) ] )

    end if


!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

    call unset_viscoelastic_elem ( last, coefficients )


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

      select case ( coorsys )
        case(0)
          numstress = 3
        case(2)
          numstress = 6
        case default
          call errormsg_case_default ( 'ce_supg_fh_elem1', &
            'coorsys', int_value=coorsys )
      end select

      allocate ( eigv(ndim,ndim), eigvalue(ndim) )
      allocate ( brownf(ninti,numstress) )
      allocate ( dDi(ndim,ndim) )

      allocate ( tmp(ndf,ndim), work2(ndim,ndim) )
      allocate ( work(ninti) )
      allocate ( tau(ninti), hoverU(ninti) )
      allocate ( ungradtheta(ninti,ndfc) )

      allocate ( posc(ndfc), u(ndim*ndf) )
      allocate ( uvecc(ndim), uvecn(ninti,ndim), cn(ndfc,ncompc) )
      allocate ( cng(ninti,ncompc,nmodes) )
      allocate ( fng(ninti,ncompc,nmodes) )
      allocate ( work6(ndfc,ncompc,mode2-mode1+1) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( eigv, eigvalue )
      deallocate ( brownf )
      deallocate ( dDi )

      deallocate ( tmp, work2 )
      deallocate ( work )
      deallocate ( tau, hoverU )
      deallocate ( ungradtheta )

      deallocate ( posc, u )
      deallocate ( uvecc, uvecn, cn )
      deallocate ( cng )
      deallocate ( fng )
      deallocate ( work6 )

    end subroutine deallocate_arrays

  end subroutine ce_supg_fh_elem1


! Principal values of the tensor a for different models

  subroutine principal_components_of_a ( vemodel, mode, cdiag, adiag )

    type(vemodel_t), intent(in) :: vemodel

!   the mode number
    integer, intent(in) :: mode

!   principal components of the conformation tensor c
    real(dp), dimension(:), intent(in) :: cdiag

!   principal components of the tensor a
    real(dp), dimension(:), intent(out) :: adiag

!   The tensor a relates the irreversible rate-of-deformation tensor to the
!   elastic stress tensor tau as follows:
!
!       D_i = a * tau
!

    real(dp) :: G, lambda, fac, alpha


    select case ( vemodel%model )

    case(2)

!     Maxwell/Oldroyd

      G = vemodel%modulus(mode)
      lambda = vemodel%lambda(mode)

      fac = 1 / ( 2 * lambda * G )

      adiag = fac / cdiag

    case(3)

!     Giesekus

      G = vemodel%modulus(mode)
      lambda = vemodel%lambda(mode)
      alpha  = vemodel%nonlin(1,mode)

      fac = 1 / ( 2 * lambda * G )

      adiag = fac * ( 1 + alpha * ( cdiag - 1 ) ) / cdiag

    case default

      write(*,'(a,i0)') &
        'Error in principal_components_of_a: model not available: ', &
         vemodel%model
      stop

    end select

  end subroutine principal_components_of_a

end module viscoelastic_elements_fh_m

