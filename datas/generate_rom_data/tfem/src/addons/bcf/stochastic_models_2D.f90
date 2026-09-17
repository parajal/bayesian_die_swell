
! Copyright (C) 2006-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


!
! stochastic viscoelastic models in 2D and axisymmetrical coordinate systems
! solved using Brownian dynamics.
!

module stochastic_models_2D_m

  use kind_defs_m
  use stochastic_models_defs_m

  implicit none


contains


! step forward for 2D and axisymmetric coordinate systems (multiple points)

  subroutine ststep_2D ( stmodel, stnumpar, gradv, qn, qnp1, brownf, addvec )

    type(stmodel_t), intent(in) :: stmodel
    type(stnumpar_t), intent(in) :: stnumpar

!   velocity gradient L
!   gradv(ip,comp), is component comp in point ip
!     comp = 1 Lyy
!     comp = 2 Lxy
!     comp = 3 Lyx
!     comp = 4 Lyy
!     comp = 5 Ltt (axisymmetric hoop strain rate)
    real(dp), dimension(:,:), intent(in) :: gradv

!   field vector Q at the start and end of the time step
!   q(field,comp,ip) is component comp of the field in point ip
!     comp = 1 Qx
!     comp = 2 Qy
!     comp = 3 Qz when non-zero or
!     comp = 3 Qt axisymmetric hoop component
!   NOTE: if addvec=.true., the qnp1 vector should be intialized to parts
!   of the `right-hand side' already computed before, such as the convection
!   term in BCF.

    real(dp), dimension(:,:,:), intent(in) :: qn
    real(dp), dimension(:,:,:), intent(inout) :: qnp1

!   The Brownian force vector scaled to a variance of 1
    real(dp), dimension(:,:), intent(in) :: brownf

!   if addvec = .true., the vector qnp1 already contains data, such as the
!   convection term in BCF. If addvec = .false. qnp1 is set to zero first.
!   Default addvec = .false.
    logical, intent(in), optional :: addvec


    logical :: laddvec
    integer :: ip


    if ( present(addvec) ) then
      laddvec = addvec
    else
      laddvec = .false.
    end if

!   loop over all points

    do ip = 1, size(gradv,1)

      call ststep1_2D ( stmodel, stnumpar, gradv(ip,:), qn(:,:,ip), &
        qnp1(:,:,ip), brownf, laddvec )

    end do

  end subroutine ststep_2D


! step forward for 2D and axisymmetric coordinate systems (one point only)

  subroutine ststep1_2D ( stmodel, stnumpar, gradv, qn, qnp1, brownf, addvec )

    type(stmodel_t), intent(in) :: stmodel
    type(stnumpar_t), intent(in) :: stnumpar

!   velocity gradient L
!   gradv(comp) is component comp
!     comp = 1 Lyy
!     comp = 2 Lxy
!     comp = 3 Lyx
!     comp = 4 Lyy
!     comp = 5 Ltt (axisymmetric hoop strain rate)
    real(dp), dimension(:), intent(in) :: gradv

!   field vector Q at the start and end of the time step
!   q(field,comp) is component comp of the field
!     comp = 1 Qx
!     comp = 2 Qy
!     comp = 3 Qz when non-zero or
!     comp = 3 Qt axisymmetric hoop component
!   NOTE: if addvec=.true., the qnp1 vector should be intialized to parts
!   of the `right-hand side' already computed before, such as the convection
!   term in BCF.

    real(dp), dimension(:,:), intent(in) :: qn
    real(dp), dimension(:,:), intent(inout) :: qnp1

!   The Brownian force vector scaled to a variance of 1
    real(dp), dimension(:,:), intent(in) :: brownf

!   if addvec = .true., the vector qnp1 already contains data, such as the
!   convection term in BCF. If addvec = .false. qnp1 is set to zero first.
!   Default addvec = .false.
    logical, intent(in), optional :: addvec


    logical :: initzero
    integer :: coorsys, timeint
    real(dp) :: lambda, g2lamb, a11, a12, a21, a22, a33, tstep, fac
    real(dp) :: sq, f, b, sqmean, grdtt
    real(dp) :: a1, a2, a3, RHS, root, temp, g4lamb, gfaclamb


    if ( stmodel%model == 0 ) then

      write(*,'(/3(a/))') &
        'Error in ststep_2D: model not set.', &
        'Call create_stochastic_model and fill material parameters in', &
        'stmodel%modulus, stmodel%lambda, stmodel%nonlin'
      stop

    end if

    if ( stmodel%coorsys /= 0 .and. stmodel%coorsys /= 1 ) then

      write(*,'(/a/2a/)') &
        'Error in ststep_2D: incorrect value of coorsys.', &
        'Call create_stochastic_model with coorsys=0 (2D) or ', &
        'coorsys=1 (axisymmetric) '
      stop

    end if

!   initialize

    if ( present(addvec) ) then
      initzero = .not. addvec
    else
      initzero = .true.
    end if

    if ( initzero ) qnp1 = 0

    coorsys = stmodel%coorsys

    timeint = stnumpar%timeint

    select case ( timeint )

    case (1)

!     explicit

      call explicit_integration

    case (2,3)

!     implicit

      call implicit_integration

    case default

      write(*,'(a,i0)') &
        'Error in ststep_2D: time integration not available: ', timeint
      stop

    end select

  contains


!   explicit integration

    subroutine explicit_integration

      integer :: i

      select case ( stmodel%model )

      case(1)

!       Hookean dumbbell

        lambda = stmodel%lambda(1)
        g2lamb = 0.5_dp / lambda

!       Coefficients

        a11 = gradv(1) - g2lamb
        a12 = gradv(2)
        a21 = gradv(3)
        a22 = gradv(4) - g2lamb
        if ( coorsys == 1 ) then
          a33 = gradv(5) - g2lamb
        end if

!       Add deterministic part

        qnp1(:,1) = qnp1(:,1) + a11 * qn(:,1) + a12 * qn(:,2)
        qnp1(:,2) = qnp1(:,2) + a21 * qn(:,1) + a22 * qn(:,2)

        if ( coorsys == 1 ) then
           qnp1(:,3) = qnp1(:,3) + a33 * qn(:,3)
        end if

!       Euler forward + Brownian force

        tstep = stnumpar%tstep
        fac = sqrt ( tstep / lambda )

        qnp1 = qn + tstep * qnp1 + fac * brownf

      case(2)

!       FENE-P

        lambda = stmodel%lambda(1)
        b      = stmodel%nonlin(1,1)
        g2lamb = 0.5_dp / lambda

!       calculate <Q^2>,

        sqmean = sum ( qn ** 2 ) / stnumpar%nfield

!       force= 1/(1-<Q^2>/b)

        f = 1 / ( 1 - sqmean / b )

!       Coefficients

        a11 = gradv(1) - f * g2lamb
        a12 = gradv(2)
        a21 = gradv(3)
        a22 = gradv(4) - f * g2lamb
        if ( coorsys == 1 ) then
           a33 = gradv(5) - f * g2lamb
        else
           a33 = - f * g2lamb
        end if

        qnp1(:,1) = qnp1(:,1) + a11 * qn(:,1) + a12 * qn(:,2)
        qnp1(:,2) = qnp1(:,2) + a21 * qn(:,1) + a22 * qn(:,2)
        qnp1(:,3) = qnp1(:,3) + a33 * qn(:,3)

!       Euler forward + Brownian force

        tstep = stnumpar%tstep
        fac = sqrt ( tstep / lambda )

        qnp1 = qn + tstep * qnp1 + fac * brownf

      case(3)

!       FENE

        lambda = stmodel%lambda(1)
        b      = stmodel%nonlin(1,1)
        g2lamb = 0.5_dp / lambda

        if ( coorsys == 1 ) then
          grdtt = gradv(5)
        else
          grdtt = 0
        end if

        do i = 1, stnumpar%nfield

!         force= 1/(1-Q^2/b)

          sq  = sum ( qn(i,:) ** 2 )
          f   = 1 / ( 1 - sq / b )

!         Coefficients

          a11 = gradv(1) - f * g2lamb
          a12 = gradv(2)
          a21 = gradv(3)
          a22 = gradv(4) - f * g2lamb
          a33 = grdtt   - f * g2lamb

!         Add deterministic part

          qnp1(i,1) = qnp1(i,1) + a11 * qn(i,1) + a12 * qn(i,2)
          qnp1(i,2) = qnp1(i,2) + a21 * qn(i,1) + a22 * qn(i,2)
          qnp1(i,3) = qnp1(i,3) + a33 * qn(i,3)

        end do

!       Euler forward + Brownian force

        tstep = stnumpar%tstep
        fac = sqrt ( tstep / lambda )

        qnp1 = qn + tstep * qnp1 + fac * brownf

      case default

        write(*,'(a,i0)') &
          'Error in ststep_2D: model not available: ', stmodel%model
        stop

      end select

    end subroutine explicit_integration


!   implicit integration
!   This routine is based on:
!   C. Mangoubi, M.A. Hulsen, R. Kupferman
!   ``Numerical stability of the method of Brownian configuration fields''
!   JNNFM, (2009)

    subroutine implicit_integration

      integer :: i

      select case ( stmodel%model )

      case(1)

!       Hookean dumbbell  - implicit  backward Euler scheme

        lambda = stmodel%lambda(1)
        g2lamb = 0.5_dp / lambda

!       Coefficients

        a11 = gradv(1)
        a12 = gradv(2)
        a21 = gradv(3)
        a22 = gradv(4)
        if ( coorsys == 1 ) then
          a33 = gradv(5)
        end if

!       Add deterministic part

        qnp1(:,1) = qnp1(:,1) + a11 * qn(:,1) + a12 * qn(:,2)
        qnp1(:,2) = qnp1(:,2) + a21 * qn(:,1) + a22 * qn(:,2)

        if ( coorsys == 1 ) then
           qnp1(:,3) = qnp1(:,3) + a33 * qn(:,3)
        end if

!       Euler + Brownian force

        tstep = stnumpar%tstep
        fac = sqrt ( tstep / lambda )

        qnp1 = ( qn + tstep * qnp1 + fac * brownf ) / (1 + tstep * g2lamb)

      case(2)

!       FENE-P

        lambda = stmodel%lambda(1)
        b      = stmodel%nonlin(1,1)
        g2lamb = 0.5_dp / lambda

!       gradv coefficients

        a11 = gradv(1)
        a12 = gradv(2)
        a21 = gradv(3)
        a22 = gradv(4)
        if ( coorsys == 1 ) then
           a33 = gradv(5)
        else
           a33 = 0
        end if

!       compute deterministic part of D_j (the RHS)

        qnp1(:,1) = qnp1(:,1) + a11 * qn(:,1) + a12 * qn(:,2)
        qnp1(:,2) = qnp1(:,2) + a21 * qn(:,1) + a22 * qn(:,2)
        qnp1(:,3) = qnp1(:,3) + a33 * qn(:,3)

!       add Brownian force to deterministic part
!       this gives what I called D_j in my notes

        tstep = stnumpar%tstep
        fac = sqrt ( tstep / lambda )

        qnp1 = qn + tstep * qnp1 + fac * brownf

!       compute < D_j ** D_j >  and extract cubic root
!       there is only one computation because it is the
!       same root for all dumbbells

        RHS =  sum ( qnp1 ** 2 ) / stnumpar%nfield

        a1 = - ( ( 2 + tstep / lambda ) * b + RHS )

        a2 = ( 1 + tstep * g2lamb ) ** 2 * b ** 2 + 2 * b * RHS

        a3 = - b ** 2 * RHS

        call solvefenepcubic ( a1 , a2 , a3 , root )

!       check length of trace < QQ >

        if ( root >= b  ) then

          write(*,'(/2(a/),a,e13.4)') &
            'Error in ststep_2D: ', &
            '  FENE-P model maximum trace length exceeded ', &
            '  trace < Q Q > = ', root
          stop

        end if

!       recover the updated vector

        temp = ( 1 + tstep * g2lamb / ( 1 - root / b ) ) ** ( -1 )

        qnp1 = temp * qnp1


      case(3)

!       FENE implicit backward Euler scheme / Crank-Nicolson

        lambda = stmodel%lambda(1)
        b      = stmodel%nonlin(1,1)
        g2lamb = 0.5_dp / lambda
        g4lamb = 0.25_dp / lambda


        if ( coorsys == 1 ) then
          grdtt = gradv(5)
        else
          grdtt = 0
        end if

        do i = 1, stnumpar%nfield

!         force= 1/(1-Q^2/b)

          sq  = sum ( qn(i,:) ** 2 )
          f   = 1 / ( 1 - sq / b )

!         Coefficients

          if ( timeint == 2 ) then

!           Backward Euler

            a11 = gradv(1)
            a12 = gradv(2)
            a21 = gradv(3)
            a22 = gradv(4)
            a33 = grdtt
            gfaclamb = g2lamb

          else if ( timeint == 3 ) then

!           Crank - Nicolson

            a11 = gradv(1) - g4lamb * f
            a12 = gradv(2)
            a21 = gradv(3)
            a22 = gradv(4) - g4lamb * f
            a33 = grdtt - g4lamb * f
            gfaclamb = g4lamb

          end if

!         Compute RHS; deterministic part

          qnp1(i,1) = qnp1(i,1) + a11 * qn(i,1) + a12 * qn(i,2)
          qnp1(i,2) = qnp1(i,2) + a21 * qn(i,1) + a22 * qn(i,2)
          qnp1(i,3) = qnp1(i,3) + a33 * qn(i,3)

        end do

!       Add Brownian force

        tstep = stnumpar%tstep
        fac = sqrt ( tstep / lambda )

        qnp1 = qn + tstep * qnp1 + fac * brownf ! RHS vector


!       Take norm and extract cubic root

        do i = 1, stnumpar%nfield

          RHS = sqrt ( sum ( qnp1 (i,:) ** 2 ) )

          a1 = -RHS
          a2 = -b * ( 1 + tstep * gfaclamb )
          a3 = b * RHS

          call solvecubic ( a1 , a2 , a3 , root )


!         Recover the updated vector

          temp = 1 - root ** 2 / b

          qnp1(i,:) = ( 1 + tstep * gfaclamb / temp ) ** ( -1 ) * qnp1(i,:)

        end do


      case default

        write(*,'(a,i0)') &
          'Error in ststep_2D: model not available: ', stmodel%model
        stop

      end select

    end subroutine implicit_integration


!   computing the cubic root using trigonometric solutions for FENE

    subroutine solvecubic ( a1 , a2 , a3 , solution)

      use math_defs_m

      real(dp) , intent(in) :: a1 , a2 , a3
      real(dp) , intent(out) :: solution
      real(dp) :: q, r, theta!, dd

      q = ( 3 * a2 - a1 ** 2 ) / 9
      r = ( 9 * a1 * a2 - 27 * a3 - 2 * a1 ** 3 ) / 54
      !dd = q ** 3 + r ** 2
      theta = acos ( r / sqrt ( - q ** 3 ) )
      solution = 2 * sqrt ( - q ) * cos ( ( theta + 4 * pi ) / 3 ) - a1 / 3

    end subroutine solvecubic


!   computing the cubic root using trigonometric solutions for FENE-P

    subroutine solvefenepcubic ( a1 , a2 , a3 , solution)

      use math_defs_m

      real(dp) , intent(in) :: a1 , a2 , a3
      real(dp) , intent(out) :: solution
      real(dp) :: q, r, theta!, dd

      q = ( 3 * a2 - a1 ** 2 ) / 9
      r = ( 9 * a1 * a2 - 27 * a3 - 2 * a1 ** 3 ) / 54
      !dd = q ** 3 + r ** 2
      theta = acos ( r / sqrt ( - q ** 3 ) )
      solution = 2 * sqrt ( - q ) * cos ( ( theta + 2 * pi ) / 3 ) - a1 / 3

    end subroutine solvefenepcubic

  end subroutine ststep1_2D


! stress tensor for 2D and axisymmetric coordinates (multiple points)

  subroutine stress_stochastic_2D ( stmodel, q, tau )

    type(stmodel_t), intent(in) :: stmodel

!   field vector Q
!   q(field,comp,ip) is component comp of the field in point ip
!     comp = 1 Qx
!     comp = 2 Qy
!     comp = 3 Qz when non-zero or
!     comp = 3 Qt axisymmetric hoop component
    real(dp), dimension(:,:,:), intent(in) :: q

!   stress tensor tau of the stochastic model
!   tau(ip,comp) is component comp in point ip
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop strain
    real(dp), dimension(:,:), intent(out) :: tau


    integer :: ip

    do ip = 1, size(tau,1)
      call stress_stochastic1_2D ( stmodel, q(:,:,ip), tau(ip,:) )
    end do

  end subroutine stress_stochastic_2D


! stress tensor for 2D and axisymmetric coordinates (single point)

  subroutine stress_stochastic1_2D ( stmodel, q, tau )

    type(stmodel_t), intent(in) :: stmodel

!   field vector Q
!   q(field,comp) is component comp of the field
!     comp = 1 Qx
!     comp = 2 Qy
!     comp = 3 Qz when non-zero or
!     comp = 3 Qt axisymmetric hoop component
    real(dp), dimension(:,:), intent(in) :: q

!   stress tensor tau of the stochastic model
!   tau(comp) is component comp
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop strain
    real(dp), dimension(:), intent(inout) :: tau ! this should have intent(out)
                                                 ! but due to a bug in lf95v6.0
                                                 ! with checking on I put it to
                                                 ! inout. The bug is gone in 6.2

    integer :: nfield, coorsys, i
    real(dp) :: G, sqmean, b, f
    real(dp), dimension(size(q,1)) :: sq, fa


    nfield = size(q,1)
    coorsys = stmodel%coorsys

    select case ( stmodel%model )

    case(1)

!     Hookean dumbbell

      G = stmodel%modulus(1)

      tau(1) = G * ( sum(q(:,1)**2) / nfield - 1 )
      tau(2) = G *   sum(q(:,1)*q(:,2)) / nfield
      tau(3) = G * ( sum(q(:,2)**2) / nfield - 1 )

      if ( coorsys == 1 ) then
        tau(4) = G * ( sum(q(:,3)**2) / nfield - 1 )
      end if

    case(2)

!     FENE-P

      G = stmodel%modulus(1)
      b = stmodel%nonlin(1,1)

!     calculate <Q^2>,

      sqmean = sum (q**2) / nfield

!     force= 1/(1-<Q^2>/b)

      f = 1 / ( 1 - sqmean / b )

!     stress

      tau(1) = G * ( f * sum(q(:,1)**2) / nfield - 1 )
      tau(2) = G *   f * sum(q(:,1)*q(:,2)) / nfield
      tau(3) = G * ( f * sum(q(:,2)**2) / nfield - 1 )
      tau(4) = G * ( f * sum(q(:,3)**2) / nfield - 1 )

    case(3)

!     FENE

      G = stmodel%modulus(1)
      b = stmodel%nonlin(1,1)

!     test length of vectors

      do i = 1, nfield
        sq(i) = sum(q(i,:)**2)
      end do

      if ( any( sq >= b ) ) then

        write(*,'(/2(a/),a,e13.4)') &
          'Error in stress-stochastic1_2D: ', &
          '  FENE model maximum length exceeded ', &
          '  max Q^2 = ', maxval(sq)
        stop

      end if

      fa = 1 / ( 1 - sq / b )

      tau(1) = G * ( sum(fa*q(:,1)**2) / nfield - 1 )
      tau(2) = G *   sum(fa*q(:,1)*q(:,2)) / nfield
      tau(3) = G * ( sum(fa*q(:,2)**2) / nfield - 1 )
      tau(4) = G * ( sum(fa*q(:,3)**2) / nfield - 1 )

    case default

      write(*,'(a,i0)') &
        'Error in stress_stochastic1_2D: model not available: ', stmodel%model
      stop

    end select

  end subroutine stress_stochastic1_2D


! structure tensor for 2D and axisymmetric coordinates

  subroutine structure_tensor_2D ( stmodel, q, c )

    type(stmodel_t), intent(in) :: stmodel

!   field vector Q
!   q(field,comp) is component comp of the field
!     comp = 1 Qx
!     comp = 2 Qy
!     comp = 3 Qz when non-zero or
!     comp = 3 Qt axisymmetric hoop component
    real(dp), dimension(:,:), intent(in) :: q

!   structure tensor c of the stochastic model
!   c(comp) is component comp
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when computed
!     comp = 4 ctt axisymmetric hoop strain
    real(dp), dimension(:), intent(inout) :: c ! this should have intent(out)
                                               ! but due to a bug in lf95v6.0
                                               ! with checking on I put it to
                                               ! inout. The bug is gone in 6.2


    integer :: nfield, coorsys


    nfield = size(q,1)
    coorsys = stmodel%coorsys

    select case ( stmodel%model )

    case(1)

!     Hookean dumbbell

      c(1) = sum(q(:,1)**2) / nfield
      c(2) = sum(q(:,1)*q(:,2)) / nfield
      c(3) = sum(q(:,2)**2) / nfield

      if ( coorsys == 1 ) then
        c(4) = sum(q(:,3)**2) / nfield
      end if

    case(2,3)

!     FENE-P, FENE

      c(1) = sum(q(:,1)**2) / nfield
      c(2) = sum(q(:,1)*q(:,2)) / nfield
      c(3) = sum(q(:,2)**2) / nfield
      c(4) = sum(q(:,3)**2) / nfield

    case default

      write(*,'(a,i0)') &
        'Error in structure_tensor_2D: model not available: ', stmodel%model
      stop

    end select

  end subroutine structure_tensor_2D

end module stochastic_models_2D_m

