
! Copyright (C) 2005-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


!
! viscoelastic models in 2D and axisymmetrical coordinate systems (log)
!

module viscoelastic_models_2D_log_m

  use glob_defs_m
  use viscoelastic_models_defs_m
  use viscoelastic_models_adapted_lambda_m
  use viscoelastic_models_2D_m, only: stress_viscoelastic_2D_single_mode, &
      vonmises_2D, vm_2D_sm, vm_2D_mm
  use viscoelastic_models_plastic_strain_m
  use eig2D3D_m, only: eig2x2, inveig2x2

  implicit none

contains


! right-hand side for 2D and axisymmetric viscoelastic models.

  subroutine rhs_viscoelastic_2D_log ( vemodel, gradv, s, rhs, mode1, mode2, &
    mvemodel, vemopt )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(s,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   velocity gradient L
!   gradv(ip,comp) is component comp in point ip with
!     comp = 1 Lxx
!     comp = 2 Lxy
!     comp = 3 Lyx
!     comp = 4 Lyy
!     comp = 5 Ltt (axisymmetric hoop strain rate)
    real(dp), dimension(:,:), intent(in) :: gradv

!   (log) conformation tensor s
!   s(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 sxx
!     comp = 2 sxy
!     comp = 3 syy
!     comp = 4 szz when non-zero or
!     comp = 4 stt axisymmetric hoop strain
!     comp = 5 xi = log(lambda) for XPP
    real(dp), dimension(:,:,:), intent(in) :: s

!   right-hand side of constitutive visco-elastic models:
!        .
!        s = rhs(s,L)
!
!   where s = log c and the conformation tensor is given by
!        .
!        c =  L * c + c * L^T - relaxation.
!
!   rhs(ip,comp,mode) is component comp in point ip for mode number mode
!   components similar to the conformation tensor c
    real(dp), dimension(:,:,:), intent(inout) :: rhs

!   if present: limit the range of modes to mode1--mode2
    integer, intent(in), optional :: mode1, mode2

!   if present: additional optional parameters (see type description)
    type(vemopt_t), intent(inout), optional :: vemopt


    logical :: deviatoric
    real(dp), parameter :: tol = 1.e-10_dp ! cut-off for equal eigenvalues
    integer  :: mode, flowtype, ncompc, np, ip, modenr1, modenr2, vm1, vm2
    real(dp) :: L(2,2), Le(2,2), r(2,2), rx(2,2), Lzz
    real(dp) :: eigvalue(2), eigv(2,2), ddss, trL3
    real(dp) :: s1, s2, s3, c1, c2, c3, lambdab, rxi, rzz
    real(dp) :: rlxm(vemodel%ncompc), sm(vemodel%ncompc)
    real(dp) :: rmod(vemodel%ncompc)
    type(vemmod_t) :: vemmod


    if ( vemodel%model == 0 ) then

      write(*,'(/3(a/))') &
        'Error in rhs_viscoelastic_2D_log: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 0 .and. vemodel%flowtype /= 1 ) then

      write(*,'(/a/2a/)') &
        'Error in rhs_viscoelastic_2D_log: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=0 (2D) or ', &
        'flowtype=1 (2D, axisymmetric) '
      stop

    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(s,1) ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D_log:', &
          ' dimension of mvemodel must be identical to size(s,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D_log:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( vemodel%ixsi > 0  ) then

!     models with slip

      write(*,'(/2(a/))') &
        'Error in rhs_viscoelastic_2D_log:', &
        '  models with xsi-parameter not available for log conformation'
      stop

    end if

    if ( vemodel%model == 17 ) then

!     XPP log conformation backbonestretch

      write(*,'(/3(a/))') &
        'Error in rhs_viscoelastic_2D_log:', &
        '  XPP log conformation of backbonestretch not available for', &
        '  full log conformation. Use rhs_viscoelastic_2D instead.'
      stop

    end if

    if ( present(mode1) .and. present(mode2) ) then
      if ( mode1 < 1 .or. mode1 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D_log:', &
          ' mode1 out of range '
        stop
      end if
      if ( mode2 < 1 .or. mode2 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D_log:', &
          ' mode2 out of range '
        stop
      end if
      modenr1 = mode1
      modenr2 = mode2
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    call check ( vemodel, vemopt, 'rhs_viscoelastic_2D_log', &
      components=['dep_J     ', 'dep_gammap', 'drhsdJ    '] )

    if ( present(vemopt) ) then

!     copy logicals from vemopt to vemmod
      vemmod%dep_J = vemopt%dep_J
      vemmod%dep_gammap = vemopt%dep_gammap

      if ( vemodel%vmglobal .and. &
                 any([vemopt%vmmode1,vemopt%vmmode2]==0) ) then
!       set vmmode to modenr
        vemmod%vmmode1 = modenr1
        vemmod%vmmode2 = modenr2
      else
!       set vmmode to values in vemopt
        vemmod%vmmode1 = vemopt%vmmode1
        vemmod%vmmode2 = vemopt%vmmode2
      end if

      if ( vemopt%compute_drhs .or. vemopt%compute_drhs_mm .or. &
           vemopt%compute_drhsdgammap ) then
        write(*,'(/2(a/))') 'Error in rhs_viscoelastic_2D_log:', &
          ' Jacobian for log c formulation not available.'
        stop
      end if

    else if ( vemodel%vmglobal ) then

!     set default for vmglobal=T: modes for von Mises equal to modenr1:modenr2
      vemmod%vmmode1 = modenr1
      vemmod%vmmode2 = modenr2

    end if

!   always allocate arrays in vemmod for the associate contruct

    allocate(vemmod%drlxmod(vemodel%ncompc,vemodel%ncompc))
    allocate(vemmod%drlxmoddJ(vemodel%ncompc))

    if ( vemodel%vmglobal ) then
      allocate ( vemmod%s(size(s,2),size(s,3)) )
      vm1 = vemmod%vmmode1
      vm2 = vemmod%vmmode2
    end if

    ncompc = vemodel%ncompc
    flowtype = vemodel%flowtype
    deviatoric = vemodel%deviatoric
    np = size(s,1)

    do ip = 1, np

!     set J

      if ( vemmod%dep_J ) vemmod%J = vemopt%J(ip)

!     set gammap

      if ( vemmod%dep_gammap ) vemmod%gammap = vemopt%gammap(ip)

!     von Mises for multiple modes

      if ( vemodel%vmglobal ) vemmod%s(:,vm1:vm2) = s(ip,:,vm1:vm2)

      do mode = modenr1, modenr2

!       extract sm

        sm = s(ip,1:ncompc,mode)

!       extract L

        L(1,1) = gradv(ip,1)
        L(1,2) = gradv(ip,2)
        L(2,1) = gradv(ip,3)
        L(2,2) = gradv(ip,4)
        if ( flowtype == 1 ) then
          Lzz  = gradv(ip,5)  ! axi-symmetric
        else
          Lzz = 0
        end if

!       deviatoric

        if ( deviatoric ) then

!         replace L with L^d = L - (tr L)/3 I

          trL3 = ( L(1,1) + L(2,2) + Lzz ) / 3
          L(1,1) = L(1,1) - trL3
          L(2,2) = L(2,2) - trL3
          Lzz = Lzz - trL3

        end if

!       compute principal values/directions

        call eig2x2 ( sm(1:3), eigvalue, eigv )

        s1 = eigvalue(1)
        s2 = eigvalue(2)

        c1 = exp(s1)
        c2 = exp(s2)

        if ( ncompc >= 4 ) then
          s3 = sm(4)
          c3 = exp(s3)
        end if

        if ( vemodel%model == 14 ) then
!         XPP double equation
          lambdab = exp(sm(5))
        end if

!       components of velocity gradient in principal system

        Le = matmul(transpose(eigv),matmul(L,eigv))

!       right-hand side part of bilinear terms

        r(1,1) = 2 * Le(1,1)
        r(2,2) = 2 * Le(2,2)

        if ( abs(c1-c2) > tol ) then
          r(1,2) = (s1 - s2)*(c2*Le(1,2)+c1*Le(2,1))/(c1-c2)
        else
          r(1,2) = Le(1,2)+Le(2,1)
        end if
        r(2,1) = r(1,2)

!       * zz/tt component (axi-symmetric)

        if ( flowtype == 1 .or. deviatoric ) then
           rzz = 2 * Lzz
        else
           rzz = 0
        end if

!       xpp model

        if ( vemodel%model == 14 ) then
          ddss   = Le(1,1) * c1 + Le(2,2) * c2
          r(1,1) = r(1,1) - 2._dp*ddss
          r(2,2) = r(2,2) - 2._dp*ddss
          rzz    = rzz    - 2._dp*ddss
          rxi    = ddss
        end if

        if ( present(mvemodel) ) then
          call rhs_relaxation_2D_log ( mvemodel(ip), mode, c1, c2, c3, &
            lambdab, rlxm, vemmod )
        else
          call rhs_relaxation_2D_log ( vemodel, mode, c1, c2, c3, lambdab, &
            rlxm, vemmod )
        end if

!       Add relaxation term to diagonal

        r(1,1) = r(1,1) - rlxm(1)
        r(2,2) = r(2,2) - rlxm(2)
        if ( ncompc >= 4 ) then
          rzz = rzz - rlxm(3)
        end if
        if ( vemodel%model == 14 ) then
          rxi = rxi - rlxm(4)
        end if

!       * transform back to global system

        rx = matmul(eigv,matmul(r,transpose(eigv)))  ! should be made
                                                     ! symmetric for speed

        rmod(1) = rx(1,1)
        rmod(2) = rx(1,2)
        rmod(3) = rx(2,2)
        if ( ncompc >= 4 ) then
          rmod(4) = rzz
        end if
        if ( vemodel%model == 14 ) then
          rmod(5) = rxi
        end if

!       fill right-hand side

        rhs(ip,:ncompc,mode) = rmod

      end do

    end do

  end subroutine rhs_viscoelastic_2D_log


! relaxation term in the models

  subroutine rhs_relaxation_2D_log ( vemodel, mode, c1, c2, c3, lambdab, &
    rlxm, vemmod )

    use limits_m, only: USE_EGP_STRESS_TENSOR_FORM, &
                        USE_ADAP_LAMBDA_FOR_EGP_DET_STAB

    type(vemodel_t), intent(in) :: vemodel

!   the mode number
    integer, intent(in) :: mode

!   principal components of the conformation tensor c
    real(dp), intent(in) :: c1, c2, c3

!   the backbone stretch in the xpp model
    real(dp), intent(in) ::  lambdab

!   relaxation term for a single mode in a single point for
!        .
!        c =  L * c + c * L^T - relaxation(c).
!
!   rlxm(comp) is component comp=1,2,3 (principal directions) of
!
!      c^{-1} * relaxation(c)

!                                    .
!   which is the relaxation term for s with s=log c.

    real(dp), dimension(:), intent(out) :: rlxm

!   additional optional parameters (see type description)
    type(vemmod_t), intent(inout) :: vemmod


    integer :: flowtype, nr
    real(dp) :: lambda, alpha, fac, eps, yfac, expfac, lambdas, nu
    real(dp) :: yfacgl, f, feq, bpar, trb, rL2, relmod, trc, sqtrc
    real(dp) :: taud, taur, rel1, rel2, beta, delta, trc3
    real(dp) :: c(vemodel%ncompc), eta, tau_y, mu, trtau, tau_e, chpar
    real(dp) :: tau(vemodel%ncompt)


    if ( mode == vemodel%norlxmode ) then

!     no relaxation term: set relaxation term to zero
      rlxm = 0

!     add EGP determinant stabilization even if relaxation is set to zero
      call add_egp_det_stabilization

      return

    end if


    flowtype = vemodel%flowtype

    select case ( vemodel%model )

    case(2)

!     Maxwell/Oldroyd

      lambda = vemodel%lambda(mode)

      rlxm(1) = ( c1 - 1 ) / c1 / lambda
      rlxm(2) = ( c2 - 1 ) / c2 / lambda
      if ( flowtype == 1 ) then
        rlxm(3) = ( c3 - 1 ) / c3 / lambda
      end if

    case(3)

!     Giesekus

      lambda = vemodel%lambda(mode)
      alpha  = vemodel%nonlin(1,mode)

      rlxm(1) = ( c1 - 1 + alpha * ( c1 - 1 ) ** 2 ) / c1 / lambda
      rlxm(2) = ( c2 - 1 + alpha * ( c2 - 1 ) ** 2 ) / c2 / lambda
      if ( flowtype == 1 ) then
        rlxm(3) = ( c3 - 1 + alpha * ( c3 - 1 ) ** 2 ) / c3 / lambda
      end if

    case(5,6)

!     Phan-Thien/Tanner

      lambda = vemodel%lambda(mode)
      eps    = vemodel%nonlin(1,mode)

      if ( vemodel%model == 5 ) then

!        linear factor

         if ( flowtype == 0 ) then
            yfac = 1 + eps * ( c1 + c2 - 2 )
         else if ( flowtype == 1 ) then
            yfac = 1 + eps * ( c1 + c2 + c3 - 3 )
         end if

      else if ( vemodel%model == 6 ) then

!       exponential factor

        if ( flowtype == 0 ) then
           yfac = exp ( eps * ( c1 + c2 - 2 ) )
        else if ( flowtype == 1 ) then
           yfac = exp ( eps * ( c1 + c2 + c3 - 3 ) )
        end if

     end if

     yfacgl = yfac / lambda

     rlxm(1) = yfacgl * ( c1 - 1 ) / c1
     rlxm(2) = yfacgl * ( c2 - 1 ) / c2
     if ( flowtype == 1 ) then
        rlxm(3) = yfacgl * ( c3 - 1 ) / c3
     end if

    case(8)

!     Chilcott-Rallison

      lambda = vemodel%lambda(mode)
      rL2    = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
        trb = c1 + c2 + 1
      else if ( flowtype == 1 ) then
        trb = c1 + c2 + c3
      end if

      f = 1 / ( 1 - trb / rL2 )
      fac = f / lambda
      rlxm(1) = fac * ( c1 - 1 ) / c1
      rlxm(2) = fac * ( c2 - 1 ) / c2
      if ( flowtype == 1 ) then
        rlxm(3) = fac * ( c3 - 1 ) / c3
      end if

    case(9)

!     FENE-P

      lambda = vemodel%lambda(mode)
      bpar   = vemodel%nonlin(1,mode)

      trb = c1 + c2 + c3
      f   = 1 / ( 1 - trb / ( bpar + 3 ) )
      feq = ( bpar + 3 ) / bpar

      rlxm(1) = ( f - feq / c1 ) / lambda
      rlxm(2) = ( f - feq / c2 ) / lambda
      rlxm(3) = ( f - feq / c3 ) / lambda

    case(14)

!     XPP double equation

      lambda  = vemodel%lambda(mode)
      lambdas = vemodel%nonlin(1,mode)
      nu      = vemodel%nonlin(2,mode)

      expfac  = exp ( nu * ( lambdab - 1._dp ) )
      relmod  = 1._dp / ( 3._dp * lambda * lambdab**2 )

      rlxm(1)  = relmod * ( 3._dp - 1._dp/c1 )
      rlxm(2)  = relmod * ( 3._dp - 1._dp/c2 )
      rlxm(3)  = relmod * ( 3._dp - 1._dp/c3 )
      rlxm(4)  = expfac / lambdas * ( 1._dp - 1._dp/(lambdab**2) )

    case(15)

!     Rolie-Poly

      taud    = vemodel%lambda(mode)
      taur    = vemodel%nonlin(1,mode)
      beta    = vemodel%nonlin(2,mode)
      delta   = vemodel%nonlin(3,mode)

      trc     = c1 + c2 + c3
      sqtrc   = sqrt(3._dp/trc)
      rel1    = ( 1._dp - sqtrc ) / taur
      rel2    = beta * ( trc / 3._dp ) ** delta

      rlxm(1) = ( 1._dp - 1._dp / c1  ) / taud + &
                rel1 * ( 1._dp + rel2 * ( 1._dp - 1._dp / c1 ) )
      rlxm(2) = ( 1._dp - 1._dp / c2 ) / taud +  &
                rel1 * ( 1._dp + rel2 * ( 1._dp - 1._dp / c2 ) )
      rlxm(3) = ( 1._dp - 1._dp / c3 ) / taud  + &
                rel1 * ( 1._dp + rel2 * ( 1._dp - 1._dp / c3 ) )

    case(16)

!     XPP single equation

      lambda   = vemodel%lambda(mode)
      lambdas  = vemodel%nonlin(1,mode)
      nu       = vemodel%nonlin(2,mode)

      trc      = c1 + c2 + c3
      trc3     = 3._dp / trc
      expfac   = 2._dp * exp ( nu * ( sqrt(trc/3._dp) - 1._dp ) )
      expfac   = expfac * ( 1._dp - trc3 ) / lambdas

      rlxm(1)  = expfac + ( trc3 - 1._dp / c1 ) / lambda
      rlxm(2)  = expfac + ( trc3 - 1._dp / c2 ) / lambda
      rlxm(3)  = expfac + ( trc3 - 1._dp / c3 ) / lambda

    case(18,19)

!     PTT-XPP single equation
!       model=18: four components in 2D
!       model=19: three components in 2D, czz=1

      lambda   = vemodel%lambda(mode)
      lambdas  = vemodel%nonlin(1,mode)
      nu       = vemodel%nonlin(2,mode)

      if ( vemodel%model == 19 .and. flowtype == 0 ) then
        trc = c1 + c2 + 1._dp
      else
        trc = c1 + c2 + c3
      end if
      trc3     = 3._dp / trc
      expfac   = 2._dp * exp ( nu * ( sqrt(trc/3._dp) - 1._dp ) )
      expfac   = expfac * ( 1._dp - trc3 ) / lambdas

      rlxm(1)  = expfac - expfac/c1 + ( trc3 - trc3 / c1 ) / lambda
      rlxm(2)  = expfac - expfac/c2 + ( trc3 - trc3 / c2 ) / lambda
      if ( vemodel%model == 18 .or. flowtype == 1 ) then
        rlxm(3)  = expfac - expfac/c3 + ( trc3 - trc3 / c3 ) / lambda
      end if

    case(20)

!     FENE-P (see Wapperom & Hulsen 1998, Huetter, Hulsen & Anderson 2018)

      lambda = vemodel%lambda(mode)
      bpar   = vemodel%nonlin(1,mode)

      trb = c1 + c2 + c3
      f   = bpar / ( bpar + 3 - trc )

      rlxm(1) = ( f - 1 / c1 ) / lambda
      rlxm(2) = ( f - 1 / c2 ) / lambda
      rlxm(3) = ( f - 1 / c3 ) / lambda

    case(21)

!     Chilcott-Rallison (alternative formulation)

      lambda = vemodel%lambda(mode)
      rL2    = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
        trb = c1 + c2 + 1
      else if ( flowtype == 1 ) then
        trb = c1 + c2 + c3
      end if

      f = ( rl2 - 3 ) / ( rL2 - trc )
      fac = f / lambda
      rlxm(1) = fac * ( c1 - 1 ) / c1
      rlxm(2) = fac * ( c2 - 1 ) / c2
      if ( flowtype == 1 ) then
        rlxm(3) = fac * ( c3 - 1 ) / c3
      end if

    case(22)

!     Saramito elastoviscoplastic (Drucker-Prager)

      eta = vemodel%modulus(mode) * vemodel%lambda(mode)
      tau_y = vemodel%nonlin(1,mode)
      mu = vemodel%nonlin(2,mode)

      c(1:4) = [ c1, 0._dp, c2, c3 ]

      call stress_viscoelastic_2D_single_mode ( vemodel, c, tau, mode )

      trtau = tau(1) + tau(3) + tau(4)  ! trace tau
      tau_e = vonmises_2D ( tau )       ! von Mises equivalent shear stress

!     find the regime the stress state is in

      chpar = mu * trtau - 3 * tau_y

      if ( chpar <= - 3 * tau_e ) then

!       regime I: sticking

        rlxm = 0

      else if ( chpar >= 2 * mu**2 * tau_e ) then

!       regime III: loosing contact

        fac = 1 / eta
        f = tau_y / mu

        rlxm(1) = fac * ( tau(1) - f ) / c1
        rlxm(2) = fac * ( tau(3) - f ) / c2
        rlxm(3) = fac * ( tau(4) - f ) / c3

      else

!       regime II: sliding

        fac = ( tau_e - tau_y + mu * trtau / 3 ) &
                    / ( ( 1 + 2*mu**2/3 ) * tau_e * eta )
        f = trtau / 3 - 2 * mu * tau_e / 3

        rlxm(1) = fac * ( tau(1) - f ) / c1
        rlxm(2) = fac * ( tau(3) - f ) / c2
        rlxm(3) = fac * ( tau(4) - f ) / c3

      end if

    case(23,24)

!     EGP

      lambda = vemodel%lambda(mode)

      if ( USE_EGP_STRESS_TENSOR_FORM ) then

        eta = vemodel%modulus(mode) * lambda

        c(1:4) = [ c1, 0._dp, c2, c3 ]

        call stress_viscoelastic_2D_single_mode ( vemodel, c, tau, mode, &
          vemmod )

        rlxm(1:3) = tau([1,3,4]) / eta

      else

        trc = c1 + c2 + c3

        rlxm(1:3) = ( [ c1, c2, c3 ] - trc/3 * [ 1, 1, 1 ] ) / lambda

        if ( vemodel%model == 24 ) rlxm(1:3) = rlxm(1:3) / vemmod%J

      end if

      if ( USE_ADAP_LAMBDA_FOR_EGP_DET_STAB ) then
!       Add stabilization for determinant here
        rlxm(1:3) = rlxm(1:3) + ( c1*c2*c3 - 1 ) / 3 * [ c1, c2, c3 ] / lambda
      end if

    case default

      write(*,'(a,i0)') &
        'Error in rhs_relaxation_2D_log: model not available: ', vemodel%model
      stop

    end select


    if ( vemodel%alam_model > 0 ) then

!     adapted lambda

      c(1:3) = [ c1, 0._dp, c2 ]
      if ( vemodel%ncompc >= 4 ) then
        c(4) = c3
        nr = 3
      else
        nr = 2
      end if

      call rhs_adapted_lambda ( vemodel, mode, c, vonmises=vm_2D_sm_log, &
        vonmises_mm=vm_2D_mm_log, rlxmod=rlxm(1:nr), vemmod=vemmod )

    end if

    if ( .not. USE_ADAP_LAMBDA_FOR_EGP_DET_STAB ) then
!     Add stabilization for determinant here with original lambda constant
      call add_egp_det_stabilization
    end if


    if ( subtractlinear ) then

!     Subtract the linear part.

      write(*,'(/a/)') &
        'Error in rhs_relaxation_2D_log: subtractlinear not available: '
      stop

    end if

  contains


!   add EGP determinant stabilization term

    subroutine add_egp_det_stabilization

      real(dp) :: lambda

      if ( all(vemodel%model/=[23,24]) ) return

      lambda = vemodel%lambda(mode)

 !    Add stabilization for determinant

      rlxm(1:3) = rlxm(1:3) + ( c1*c2*c3 - 1 ) / 3 * [ c1, c2, c3 ] / lambda

    end subroutine add_egp_det_stabilization

  end subroutine rhs_relaxation_2D_log


! von Mises equivalent shear stress for a single mode

  subroutine vm_2D_sm_log ( vemodel, c, mode, vemmod, vm, dvm, dvmdJ )

    type(vemodel_t), intent(in) :: vemodel
!   c(:ncompc)
    real(dp), dimension(:), intent(in) :: c
    integer, intent(in) :: mode
    type(vemmod_t), intent(in) :: vemmod
    real(dp), intent(out), optional :: vm
!   dvm(:ncompc)
    real(dp), dimension(:), intent(out), optional :: dvm
    real(dp), intent(out), optional :: dvmdJ

    real(dp) :: tau(vemodel%ncompt)

    call stress_viscoelastic_2D_single_mode ( vemodel, c, tau, mode, vemmod )

    if ( present(vm) ) vm = vonmises_2D ( tau )

    if ( present(dvm) .or. present(dvmdJ) ) then
      write(*,'(a,i0)') &
        'Error in vm_2D_sm_log: derivative not implemented'
      stop
    end if

  end subroutine vm_2D_sm_log


! von Mises equivalent shear stress based on multiple modes

  subroutine vm_2D_mm_log ( vemodel, vemmod, vm, dvm, dvmdJ )

    type(vemodel_t), intent(in) :: vemodel
    type(vemmod_t), intent(in) :: vemmod
    real(dp), intent(out), optional :: vm
!   dvm(:ncompc,:nmodes), only mode range vmmode1:vmmode2 is computed
    real(dp), dimension(:,:), intent(out), optional :: dvm
    real(dp), intent(out), optional :: dvmdJ

    integer :: mode
    real(dp), dimension(vemodel%ncompt) :: tau, tsm

    tau = 0

    do mode = vemmod%vmmode1, vemmod%vmmode2

      call stress_viscoelastic_2D_single_mode_log ( vemodel, &
        vemmod%s(:,mode), tsm, mode, vemmod )

      tau = tau + tsm

    end do

    if ( present(vm) ) vm = vonmises_2D(tau)

    if ( present(dvm) .or. present(dvmdJ) ) then
      write(*,'(a,i0)') &
        'Error in vm_2D_mm_log: derivative not implemented'
      stop
    end if

  end subroutine vm_2D_mm_log


! stress tensor for 2D and axisymmetric viscoelastic models (single mode)

  subroutine stress_viscoelastic_2D_single_mode_log ( vemodel, s, tau, mode, &
    vemmod )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   log conformation tensor s = log c
!   s(comp) is component comp
!     comp = 1 sxx
!     comp = 2 sxy
!     comp = 3 syy
!     comp = 4 szz when non-zero or
!     comp = 4 stt axisymmetric hoop strain
    real(dp), dimension(:), intent(in) :: s

!   stress tensor tau of the viscoelastic model
!   tau(comp) is component comp
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop strain
    real(dp), dimension(:), intent(out) :: tau

!   the mode number
    integer, intent(in) :: mode

!   additional optional parameters (see type description)
!   NOTE: this argument is set to optional for compatibility reasons.
!   It is checked to be present, when needed.
    type(vemmod_t), intent(in), optional :: vemmod


    real(dp) :: cm(vemodel%ncompc)
    real(dp) :: cmdiag(2), eigvalue(2), eigv(2,2)


!   extract single mode value of s

!   compute exp(s)

    call eig2x2 ( s(1:3), eigvalue, eigv )

    cmdiag(1) = exp(eigvalue(1))
    cmdiag(2) = exp(eigvalue(2))

    call inveig2x2 ( cm(1:3), cmdiag, eigv )

    if ( vemodel%ncompc >= 4 ) then
      cm(4) = exp(s(4))
    end if
    if ( vemodel%model == 14 ) then ! XPP
      cm(5) = exp(s(5))
    end if

    call stress_viscoelastic_2D_single_mode ( vemodel, cm, tau, mode, vemmod )

  end subroutine stress_viscoelastic_2D_single_mode_log


! stress tensor for 2D and axisymmetric viscoelastic models.

  subroutine stress_viscoelastic_2D_log ( vemodel, s, tau, mode, mode1, &
    mode2, mvemodel, vemopt )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(s,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   log conformation tensor s = log c
!   s(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 sxx
!     comp = 2 sxy
!     comp = 3 syy
!     comp = 4 szz when non-zero or
!     comp = 4 stt axisymmetric hoop strain
    real(dp), dimension(:,:,:), intent(in) :: s

!   stress tensor tau of the viscoelastic model
!   tau(ip,comp) is component comp in point ip
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop strain
    real(dp), dimension(:,:), intent(out) :: tau

!   when present only a single mode contribution is computed
    integer, intent(in), optional :: mode

!   if present: limit the range of modes to mode1--mode2
    integer, intent(in), optional :: mode1, mode2

!   if present: control the additional optional parameters below
    type(vemopt_t), intent(in), optional :: vemopt


    integer :: modenr1, modenr2, imode, np, ip, ncompc
    real(dp) :: tsm(vemodel%ncompt), sm(vemodel%ncompc)
    type(vemmod_t) :: vemmod


    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(s,1) ) then
        write(*,'(/a/a/)') 'Error in stress_viscoelastic_2D_log:', &
          ' dimension of mvemodel must be identical to size(s,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in stress_viscoelastic_2D_log:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    call check ( vemodel, vemopt, 'stress_viscoelastic_2D_log', &
      components=['dep_J'] )

    if ( present(vemopt) ) then

!     copy logicals from vemopt to vemmod
      vemmod%dep_J = vemopt%dep_J

    end if

    np = size(s,1)
    ncompc = vemodel%ncompc

    if ( present(mode) ) then
      modenr1 = mode
      modenr2 = mode
    else if ( present(mode1) .and. present(mode2) ) then
      modenr1 = mode1
      modenr2 = mode2
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    tau = 0

    do ip = 1, np

!     set J

      if ( vemmod%dep_J ) vemmod%J = vemopt%J(ip)

      do imode = modenr1, modenr2

!       extract single mode value of s

        sm = s(ip,1:ncompc,imode)

        if ( present(mvemodel) ) then
          call stress_viscoelastic_2D_single_mode_log ( mvemodel(ip), &
             sm, tsm, imode, vemmod )
        else
          call stress_viscoelastic_2D_single_mode_log ( vemodel, &
             sm, tsm, imode, vemmod )
        end if

        tau(ip,:) = tau(ip,:) + tsm

      end do

    end do

  end subroutine stress_viscoelastic_2D_log


! conformation tensor for 2D and axisymmetric viscoelastic models.

  subroutine conformation_2D_log ( vemodel, s, c )

    type(vemodel_t), intent(in) :: vemodel

!   log conformation tensor s = log c
!   s(ip,comp) is component comp in point ip
!     comp = 1 sxx
!     comp = 2 sxy
!     comp = 3 syy
!     comp = 4 szz when non-zero or
!     comp = 4 stt axisymmetric hoop strain
!     comp = 5 xi = log(lambda) for XPP
    real(dp), dimension(:,:), intent(in) :: s

!   conformation tensor c of the viscoelastic model
!   c(ip,comp) is component comp in point ip
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambda for XPP
!     NOTE: for the XPP model comp=1:4 represent the orientation tensor S.
    real(dp), dimension(:,:), intent(out) :: c


    integer :: np, ip, ncompc
    real(dp) :: spn(vemodel%ncompc), cpn(vemodel%ncompc)
    real(dp) :: cmdiag(2), eigvalue(2), eigv(2,2)


    np = size(s,1)
    ncompc = vemodel%ncompc

    do ip = 1, np

!     extract value of s in one point

      spn = s(ip,1:ncompc)

!     compute exp(spn)

      call eig2x2 ( spn(1:3), eigvalue, eigv )

      cmdiag(1) = exp(eigvalue(1))
      cmdiag(2) = exp(eigvalue(2))

      call inveig2x2 ( cpn(1:3), cmdiag, eigv )
      if ( ncompc >= 4 ) then
        cpn(4) = exp(spn(4))
      end if
      if ( vemodel%model == 14 ) then ! XPP
        cpn(5) = exp(spn(5))
      end if

      c(ip,:) = cpn

    end do

  end subroutine conformation_2D_log


! Right-hand side for plastic strain evolution of 2D log models.

  subroutine rhs_plastic_strain_2D_log ( vemodel, s, gammap, rhs_gammap, &
    mode1, mode2, mvemodel, vemopt_gammap )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(c,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   (log) conformation tensor s
!   s(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 sxx
!     comp = 2 sxy
!     comp = 3 sxz
!     comp = 4 syy
!     comp = 5 syz
!     comp = 6 szz
!     comp = 7 xi = log(lambda) for XPP
    real(dp), dimension(:,:,:), intent(in) :: s

!   plastic strain gammap
!   gammap(ip) is plastic strain in point ip
    real(dp), dimension(:), intent(in) :: gammap

!   if present: limit the range of modes to mode1--mode2
!   NOTE: this only affects the coupling of modes for von Mises stress
    integer, intent(in), optional :: mode1, mode2

!   right-hand side of the evolution equations of the plastic strain
!   rhs_gammap(ip) in point ip
    real(dp), dimension(:), intent(inout) :: rhs_gammap

!   if present: additional optional parameters (see type description)
    type(vemopt_gammap_t), intent(inout), optional :: vemopt_gammap

    integer :: mode, np, ip, modenr1, modenr2
    real(dp) :: c(1,size(s,2))
    type(vemmod_t) :: vemmod


    if ( vemodel%model == 0 ) then

      write(*,'(/3(a/))') &
        'Error in rhs_plastic_strain_2D_log: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 0 .and. vemodel%flowtype /= 1 ) then

      write(*,'(/a/2a/)') &
        'Error in rhs_plastic_strain_2D_log: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=0 (2D) or ', &
        'flowtype=1 (2D, axisymmetric) '
      stop

    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(s,1) ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_2D_log:', &
          ' dimension of mvemodel must be identical to size(s,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_2D_log:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( present(mode1) .and. present(mode2) ) then
      if ( mode1 < 1 .or. mode1 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_2D_log:', &
          ' mode1 out of range '
        stop
      end if
      if ( mode2 < 1 .or. mode2 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_2D_log:', &
          ' mode2 out of range '
        stop
      end if
      modenr1 = mode1
      modenr2 = mode2
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    vemmod%dep_gammap = .true.

    if ( present(vemopt_gammap) ) then

!     copy logicals from vemopt_gammap to vemmod
      vemmod%dep_J = vemopt_gammap%dep_J

      if ( vemodel%vmglobal .and. &
                 any([vemopt_gammap%vmmode1,vemopt_gammap%vmmode2]==0) ) then
!       set vmmode to modenr
        vemmod%vmmode1 = modenr1
        vemmod%vmmode2 = modenr2
      else
!       set vmmode to values in vemopt_gammap
        vemmod%vmmode1 = vemopt_gammap%vmmode1
        vemmod%vmmode2 = vemopt_gammap%vmmode2
      end if

      if ( vemopt_gammap%compute_drhs .or. vemopt_gammap%compute_drhs_mm .or. &
           vemopt_gammap%compute_drhsdgammap ) then
        write(*,'(/2(a/))') 'Error in rhs_viscoelastic_2D_log:', &
          ' Jacobian for log c formulation not available.'
        stop
      end if

    else if ( vemodel%vmglobal ) then

!     set default for vmglobal=T: modes for von Mises equal to modenr1:modenr2
      vemmod%vmmode1 = modenr1
      vemmod%vmmode2 = modenr2

    end if

!   always allocate arrays in vemmod for the associate contruct

    allocate(vemmod%drlxmod(1,vemodel%ncompc))
    allocate(vemmod%drlxmoddJ(1))

    mode = vemodel%gammap_mode
    np = size(s,1)

    do ip = 1, np

!     set gammap

      vemmod%gammap = gammap(ip)

!     set J

      if ( vemmod%dep_J ) vemmod%J = vemopt_gammap%J(ip)

!     fill right-hand side

      if ( present(mvemodel) ) then
        call conformation_2D_log ( mvemodel(ip), s(ip:ip,:,mode), c )
        call rhs_plastic_strain_generic ( mvemodel(ip), mode, c(1,:), &
          vonmises=vm_2D_sm, vonmises_mm=vm_2D_mm, &
          rhs_gammap=rhs_gammap(ip), vemmod=vemmod )
      else
        call conformation_2D_log ( vemodel, s(ip:ip,:,mode), c )
        call rhs_plastic_strain_generic ( vemodel, mode, c(1,:), &
          vonmises=vm_2D_sm, vonmises_mm=vm_2D_mm, &
          rhs_gammap=rhs_gammap(ip), vemmod=vemmod )
      end if

    end do

  end subroutine rhs_plastic_strain_2D_log

end module viscoelastic_models_2D_log_m

