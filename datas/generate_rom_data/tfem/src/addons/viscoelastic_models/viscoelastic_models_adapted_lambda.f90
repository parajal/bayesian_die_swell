! Copyright (C) 2021-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


!
! Viscoelastic models with adapted lambda
!

module viscoelastic_models_adapted_lambda_m

  use glob_defs_m
  use viscoelastic_models_defs_m

  implicit none


contains


! adapt lambda in relaxation term in the models. Optional Jacobian.

  subroutine rhs_adapted_lambda ( vemodel, mode, c, vonmises, vonmises_mm, &
    rlxmod, vemmod )

    use kind_defs_m

    type(vemodel_t), intent(in) :: vemodel

!   the mode number
    integer, intent(in) :: mode

!   actual argument for c in call to vonmises
!   c(comp) is component comp
    real(dp), dimension(:), intent(in) :: c

!   function routine for the von mises shear stress for a single mode
    interface
      subroutine vonmises ( vemodel, c, mode, vemmod, vm, dvm, dvmdJ )
        use kind_defs_m
        use viscoelastic_models_defs_m, only: vemodel_t, vemmod_t
        implicit none
        type(vemodel_t), intent(in) :: vemodel
        real(dp), dimension(:), intent(in) :: c
        integer, intent(in) :: mode
        type(vemmod_t), intent(in) :: vemmod
        real(dp), intent(out), optional :: vm
        real(dp), dimension(:), intent(out), optional :: dvm
        real(dp), intent(out), optional :: dvmdJ
      end subroutine vonmises
    end interface

!   function routine for the von mises shear stress for multi mode
    interface
      subroutine vonmises_mm ( vemodel, vemmod, vm, dvm, dvmdJ )
        use kind_defs_m
        use viscoelastic_models_defs_m, only: vemodel_t, vemmod_t
        implicit none
        type(vemodel_t), intent(in) :: vemodel
        type(vemmod_t), intent(in) :: vemmod
        real(dp), intent(out), optional :: vm
        real(dp), dimension(:,:), intent(out), optional :: dvm
        real(dp), intent(out), optional :: dvmdJ
      end subroutine vonmises_mm
    end interface

!   relaxation term for a single mode in a single point
    real(dp), dimension(:), intent(inout) :: rlxmod

!   additional optional parameters (see type description)
    type(vemmod_t), intent(inout) :: vemmod


    integer :: m
    real(dp) :: tau_d, tau_y, n, K, G, lambda, f, df, dtau_ddJ
    real(dp) :: tau_ref, x, tau_ref1, f1, tau_ref2, f2, x1, x2, shx, shx1, shx2
    real(dp) :: beta, h, dhdJ
    real(dp) :: Sa, r0, r1, r2, r2m1dr1, fac, y1, R, e, dRdgp, dedgp
    real(dp), dimension(size(c)) :: dtau_d
    real(dp), dimension(size(c),vemodel%nmodes) :: dtau_d_mm


    associate ( drlxmod => vemmod%drlxmod, drlxmoddJ => vemmod%drlxmoddJ )

    select case ( vemodel%alam_model )

    case(1)

!     elastic

      rlxmod = 0
      if ( vemmod%compute_drlxmod ) call drlx_set_to_zero
      if ( vemmod%compute_drhsdJ ) drlxmoddJ = 0

    case(2)

!     elastoviscoplastic Saramito (2007)

      tau_y = vemodel%alam(1,mode)

      call vonmises_stress

      if ( tau_d > tau_y ) then
        f = ( tau_d - tau_y ) / tau_d
        df = tau_y / tau_d**2
      else
        f = 0; df = 0
      end if

      if ( vemmod%compute_drlxmod ) then

!       compute Jacobian wrt to c

        if ( tau_d > tau_y ) then

!         Compute derivatives with von Mises stress based on f and df

          call drlx_with_f_and_df

        else
          call drlx_set_to_zero
        end if

      end if

      if ( vemmod%compute_drhsdJ ) then

!       compute Jacobian wrt to J

        if ( tau_d > tau_y ) then
          drlxmoddJ = f * drlxmoddJ
          drlxmoddJ = drlxmoddJ + df * rlxmod * dtau_ddJ
        else
          drlxmoddJ = 0
        end if

      end if

      if ( tau_d > tau_y ) then
        rlxmod = f * rlxmod
      else
        rlxmod = 0
      end if

    case(3)

!     elastoviscoplastic Saramito (2009)

      G = vemodel%modulus(mode)
      lambda = vemodel%lambda(mode)
      tau_y = vemodel%alam(1,mode)
      K = vemodel%alam(2,mode)
      n = vemodel%alam(3,mode)

      call vonmises_stress

      if ( tau_d > tau_y ) then
        f = G * lambda * ( ( tau_d - tau_y ) / K ) ** (1/n) / tau_d
        df = G * lambda * &
                     ( -f/tau_d + ((tau_d - tau_y)/K)**(1/n-1)/(n*K*tau_d) )
      else
        f = 0; df = 0
      end if

      if ( vemmod%compute_drlxmod ) then

!       compute Jacobian

        if ( tau_d > tau_y ) then

!         Compute derivatives with von Mises stress based on f and df

          call drlx_with_f_and_df

        else
          call drlx_set_to_zero
        end if

      end if

      if ( vemmod%compute_drhsdJ ) then

!       compute Jacobian wrt to J

        if ( tau_d > tau_y ) then
          drlxmoddJ = f * drlxmoddJ
          drlxmoddJ = drlxmoddJ + df * rlxmod * dtau_ddJ
        else
          drlxmoddJ = 0
        end if

      end if

      if ( tau_d > tau_y ) then
        rlxmod = f * rlxmod
      else
        rlxmod = 0
      end if

    case(4)

!     Eyring

      tau_ref = vemodel%alam(1,mode)

      call vonmises_stress

      x = tau_d / tau_ref

      if ( x > tiny(x) ) then
        shx = sinh(x)
        f = shx / x
        df = ( cosh(x)*x - shx ) / x**2 / tau_ref
      end if

      if ( vemmod%compute_drlxmod ) then

!       compute Jacobian

        if ( x > tiny(x) ) then

!         Compute derivatives with von Mises stress based on f and df

          call drlx_with_f_and_df

        else if ( vemodel%vmglobal ) then

!         start from zero (x=0) for multi mode von Mises

          call drlx_start_from_zero_mm

        end if

      end if

      if ( vemmod%compute_drhsdJ ) then

!       compute Jacobian wrt to J

        if ( x > tiny(x) ) then
          drlxmoddJ = f * drlxmoddJ
          drlxmoddJ = drlxmoddJ + df * rlxmod * dtau_ddJ
        else
          drlxmoddJ = 0
        end if

      end if

      if ( x > tiny(x) ) then
        rlxmod = f * rlxmod
      end if

    case(5)

!     Ree-Eyring

      tau_ref1 = vemodel%alam(1,mode)
      f1 = vemodel%alam(2,mode)
      tau_ref2 = vemodel%alam(3,mode)
      f2 = 1 - f1

      call vonmises_stress

      x1 = tau_d / tau_ref1
      x2 = tau_d / tau_ref2

      if ( x1 > tiny(x1) ) then
        shx1 = sinh(x1)
        shx2 = sinh(x2)
        f = 1 / ( f1*x1/shx1 + f2*x2/shx2 )
        df = - f**2 * ( f1*(shx1-x1*cosh(x1))/shx1**2 / tau_ref1 + &
                        f2*(shx2-x2*cosh(x2))/shx2**2 / tau_ref2 )
      end if

      if ( vemmod%compute_drlxmod ) then

!       compute Jacobian

        if ( x1 > tiny(x1) ) then

!         Compute derivatives with von Mises stress based on f and df

          call drlx_with_f_and_df

        else if ( vemodel%vmglobal ) then

!         start from zero (x=0) for multi mode von Mises

          call drlx_start_from_zero_mm

        end if

      end if

      if ( vemmod%compute_drhsdJ ) then

!       compute Jacobian wrt to J

        if ( x1 > tiny(x1) ) then
          drlxmoddJ = f * drlxmoddJ
          drlxmoddJ = drlxmoddJ + df * rlxmod * dtau_ddJ
        else
          drlxmoddJ = 0
        end if

      end if

      if ( x1 > tiny(x1) ) then
        rlxmod = f * rlxmod
      end if

    case default

      write(*,'(a,i0)') &
        'Error in rhs_adapted_lambda: alam_model not available: ', &
        vemodel%alam_model
      stop

    end select


!   dependance on J

    if ( vemodel%alamJ_model > 0 ) then

!     h(J) factor for dependence on J

      select case ( vemodel%alamJ_model )

      case(1)

        beta = vemodel%alamJ(1,mode)

        h = vemmod%J ** beta

        if ( vemmod%compute_drlxmod ) then

          if ( vemodel%vmglobal ) then

!           multi mode von Mises

!           diagonal
            vemmod%drlxmod_mm(:,:,mode) = h * vemmod%drlxmod_mm(:,:,mode)

!           off-diagonal
            do m = vemmod%vmmode1, vemmod%vmmode2
              if ( m == mode ) cycle ! diagonal already done
              vemmod%drlxmod_mm(:,:,m) = h * vemmod%drlxmod_mm(:,:,m)
            end do

          else

!           single mode von Mises

            drlxmod = h * drlxmod

          end if

        end if

        if ( vemmod%compute_drhsdJ ) then

!         compute Jacobian wrt to J

          dhdJ = beta * h / vemmod%J
          drlxmoddJ = h * drlxmoddJ + dhdJ * rlxmod

        end if

        rlxmod = h * rlxmod

      case default

        write(*,'(a,i0)') &
          'Error in rhs_adapted_lambda: alamJ_model not available: ', &
          vemodel%alamJ_model
        stop

      end select

    end if


!   dependance on plastic strain gammap

    if ( vemodel%alam_gammap_model > 0 ) then

!     e(gammap) factor for dependence on gammap

      select case ( vemodel%alam_gammap_model )

      case(1)

!       get parameters of this mode

        Sa = vemodel%alam_gammap(1,mode)
        r0 = vemodel%alam_gammap(2,mode)
        r1 = vemodel%alam_gammap(3,mode)
        r2 = vemodel%alam_gammap(4,mode)
        r2m1dr1 = (r2-1)/r1
        fac = (1+r0**r1)**r2m1dr1

!       e(gammap)

        y1 = (r0*exp(vemmod%gammap))**r1
        R = (1+y1)**r2m1dr1/fac
        e = exp(-Sa*R)

        if ( vemmod%compute_drlxmod ) then

          if ( vemodel%vmglobal ) then

!           multi mode von Mises

!           diagonal
            vemmod%drlxmod_mm(:,:,mode) = e * vemmod%drlxmod_mm(:,:,mode)

!           off-diagonal
            do m = vemmod%vmmode1, vemmod%vmmode2
              if ( m == mode ) cycle ! diagonal already done
              vemmod%drlxmod_mm(:,:,m) = e * vemmod%drlxmod_mm(:,:,m)
            end do

          else

!           single mode von Mises

            drlxmod = e * drlxmod

          end if

        end if

        if ( vemmod%compute_drhsdJ ) then

!         compute Jacobian wrt to J

          drlxmoddJ = e * drlxmoddJ

        end if

        if ( vemmod%compute_drhsdgammap ) then

!         compute Jacobian wrt to gammap

          dRdgp = (r2-1)*y1*(1+y1)**(r2m1dr1-1)/fac
          dedgp = -Sa*e*dRdgp
          vemmod%drlxmoddgammap = dedgp * rlxmod

        end if

        rlxmod = e * rlxmod

      case default

        write(*,'(a,i0)') &
          'Error in rhs_adapted_lambda: alam_gammap_model not available: ', &
          vemodel%alam_gammap_model
        stop

      end select

    end if


    end associate


  contains

!   compute von Mises shear stress and derivatives

    subroutine vonmises_stress

      if ( vemodel%vmglobal ) then

!       von Mises stress based on multiple modes

!       von Mises shear stress

        call vonmises_mm ( vemodel, vemmod, tau_d )

        if ( vemmod%compute_drlxmod ) then

!         compute derivative wrt to c

          call vonmises_mm ( vemodel, vemmod, dvm=dtau_d_mm )

        end if

        if ( vemmod%compute_drhsdJ ) then

!         compute derivative wrt to J

          call vonmises_mm ( vemodel, vemmod, dvmdJ=dtau_ddJ )

        end if

      else

!       von Mises stress based on a single mode

!       von Mises shear stress

        call vonmises ( vemodel, c, mode, vemmod, tau_d )

        if ( vemmod%compute_drlxmod ) then

!         compute derivative wrt to c

          call vonmises ( vemodel, c, mode, vemmod, dvm=dtau_d )

        end if

        if ( vemmod%compute_drhsdJ ) then

!         compute derivative wrt to J

          call vonmises ( vemodel, c, mode, vemmod, dvmdJ=dtau_ddJ )

        end if

      end if

    end subroutine vonmises_stress


!   Compute derivatives with von Mises stress based on f and df

    subroutine drlx_with_f_and_df

      integer :: i, m

      associate ( drlxmod => vemmod%drlxmod )

      if ( vemodel%vmglobal ) then

!       multi mode von Mises

!       diagonal: derivative of tensor part
        vemmod%drlxmod_mm(:,:,mode) = f * drlxmod

        do m = vemmod%vmmode1, vemmod%vmmode2

!         set off-diagonal to zero
          if ( m /= mode ) vemmod%drlxmod_mm(:,:,m) = 0

!         add derivative of adapted lambda part

          do i = 1, size(vemmod%drlxmod_mm,1)
            vemmod%drlxmod_mm(i,:,m) = vemmod%drlxmod_mm(i,:,m) &
                                          + df * rlxmod(i) * dtau_d_mm(:,m)
          end do

        end do

      else

!       single mode von Mises

        drlxmod = f * drlxmod
        do i = 1, size(drlxmod,1)
          drlxmod(i,:) = drlxmod(i,:) + df * rlxmod(i) * dtau_d
        end do

      end if

      end associate

    end subroutine drlx_with_f_and_df


!   start from zero (x=0) for multi mode von Mises

    subroutine drlx_start_from_zero_mm

      integer :: m

      vemmod%drlxmod_mm(:,:,mode) = vemmod%drlxmod

      do m = vemmod%vmmode1, vemmod%vmmode2
        if ( m /= mode ) vemmod%drlxmod_mm(:,:,m) = 0
      end do

    end subroutine drlx_start_from_zero_mm


!   set drlxmod or drlxmod_mm to zero

    subroutine drlx_set_to_zero

      if (vemodel%vmglobal) then
        vemmod%drlxmod_mm = 0
      else
        vemmod%drlxmod = 0
      end if

    end subroutine drlx_set_to_zero

  end subroutine rhs_adapted_lambda


! von Mises equivalent shear stress (2D)

  function vonmises_2D ( tau )

!   stress tensor tau of the viscoelastic model (single mode)
!   tau(comp) is component comp
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop strain
    real(dp), dimension(:), intent(in) :: tau
    real(dp) :: vonmises_2D

    if ( size(tau) == 4 ) then
      vonmises_2D = sqrt ( ( ( tau(1) - tau(3) ) ** 2 + &
                             ( tau(3) - tau(4) ) ** 2 + &
                             ( tau(4) - tau(1) ) ** 2 ) / 6 + &
                             tau(2) ** 2 )
    else
      vonmises_2D = sqrt ( ( ( tau(1) - tau(3) ) ** 2 + &
                             tau(3) ** 2 + tau(1) ** 2 ) / 6 + &
                             tau(2) ** 2 )
    end if

  end function vonmises_2D


! derivative von Mises equivalent shear stress (2D)

  function dvonmises_2D ( tau )

!   stress tensor tau of the viscoelastic model (single mode)
!   tau(comp) is component comp
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop strain
    real(dp), dimension(:), intent(in) :: tau
    real(dp), dimension(size(tau)) :: dvonmises_2D
    real(dp) :: txx, txy, tyy, tzz, tmp

    tmp =  vonmises_2D(tau)

    if ( tmp > tiny(tmp) ) then

      txx = tau(1); txy=tau(2); tyy=tau(3)

      if ( size(tau) == 4 ) then
        tzz=tau(4)
        dvonmises_2D = [ 2*txx-tyy-tzz, txy, 2*tyy-tzz-txx, 2*tzz-txx-tyy ]
        dvonmises_2D([1,3,4]) = dvonmises_2D([1,3,4]) / 6
      else
        dvonmises_2D = [ 2*txx-tyy, txy, 2*tyy-txx ]
        dvonmises_2D([1,3]) = dvonmises_2D([1,3]) / 6
      end if

      dvonmises_2D = dvonmises_2D / tmp

    else

      dvonmises_2D = 0 ! set to zero for tau=0

    end if

  end function dvonmises_2D


! von Mises equivalent shear stress (3D)

  function vonmises_3D ( tau )

!   stress tensor tau of the viscoelastic model (single mode)
!   tau(comp) is component comp
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauxz
!     comp = 4 tauyy
!     comp = 5 tauyz
!     comp = 6 tauzz
    real(dp), dimension(:), intent(in) :: tau
    real(dp) :: vonmises_3D

    vonmises_3D = sqrt ( ( ( tau(1) - tau(4) ) ** 2 + &
                           ( tau(4) - tau(6) ) ** 2 + &
                           ( tau(6) - tau(1) ) ** 2 ) / 6 + &
                           tau(2) ** 2 + tau(3) ** 2 + tau(5) ** 2 )

  end function vonmises_3D


! derivative of von Mises equivalent shear stress (3D)

  function dvonmises_3D ( tau )

!   stress tensor tau of the viscoelastic model (single mode)
!   tau(comp) is component comp
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauxz
!     comp = 4 tauyy
!     comp = 5 tauyz
!     comp = 6 tauzz
    real(dp), dimension(:), intent(in) :: tau
    real(dp), dimension(size(tau)) :: dvonmises_3D

    real(dp) :: txx, txy, txz, tyy, tyz, tzz, tmp

    tmp =  vonmises_3D(tau)

    if ( tmp > tiny(tmp) ) then

      txx = tau(1); txy=tau(2); txz=tau(3); tyy=tau(4); tyz=tau(5); tzz=tau(6)

      dvonmises_3D = &
             [ 2*txx-tyy-tzz, txy, txz, 2*tyy-tzz-txx, tyz, 2*tzz-txx-tyy ]
      dvonmises_3D([1,4,6]) = dvonmises_3D([1,4,6]) / 6
      dvonmises_3D = dvonmises_3D / tmp

    else

      dvonmises_3D = 0 ! set to zero for tau=0

    end if

  end function dvonmises_3D

end module viscoelastic_models_adapted_lambda_m

