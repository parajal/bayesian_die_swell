! Copyright (C) 2025-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


!
! Viscoelastic models with plastic strain
!

module viscoelastic_models_plastic_strain_m

  use glob_defs_m
  use viscoelastic_models_defs_m
  use viscoelastic_models_adapted_lambda_m

  implicit none

contains


! Right-hand side for plastic strain evolution. Optional Jacobian.

  subroutine rhs_plastic_strain_generic ( vemodel, mode, c, vonmises, &
    vonmises_mm, rhs_gammap, vemmod )

    use kind_defs_m

    type(vemodel_t), intent(in) :: vemodel

!   the mode number, on which the plastic strain is based on
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

!   right-hand side of plastic strain evolution in a single point
    real(dp), intent(out) :: rhs_gammap

!   additional optional parameters (see type description)
    type(vemmod_t), intent(inout) :: vemmod


    real(dp) :: tau_d, G, lambda, eta, dtau_ddJ
    real(dp), dimension(size(c)) :: dtau_d
    real(dp), dimension(1) :: rlxmod


!   right-hand side of effective plastic strain rate = tau_d / eta for mode

!   compute tau_d for mode

    call vonmises_stress

!   tau_d / eta = tau_d / ( G * lambda )

    G = vemodel%modulus(mode)
    lambda = vemodel%lambda(mode)
    eta = G * lambda

    rhs_gammap = tau_d / eta

    if ( vemmod%compute_drlxmod ) then

!     compute Jacobian wrt to c

      vemmod%drlxmod(1,:) = dtau_d / eta

    end if

    if ( vemmod%compute_drhsdJ ) then

!     compute Jacobian wrt to J

      vemmod%drlxmoddJ(1) = dtau_ddJ / eta

    end if

!   adapted lambda

    if ( vemodel%alam_model > 0 ) then

!     adapted lambda

      rlxmod = rhs_gammap

      call rhs_adapted_lambda ( vemodel, mode, c, vonmises, vonmises_mm, &
        rlxmod, vemmod )

      rhs_gammap = rlxmod(1)

    end if

  contains

!   compute von Mises shear stress and derivatives

    subroutine vonmises_stress

!     von Mises stress based on a single mode

!     von Mises shear stress

      call vonmises ( vemodel, c, mode, vemmod, tau_d )

      if ( vemmod%compute_drlxmod ) then

!       compute derivative wrt to c

        call vonmises ( vemodel, c, mode, vemmod, dvm=dtau_d )

      end if

      if ( vemmod%compute_drhsdJ ) then

!       compute derivative wrt to J

        call vonmises ( vemodel, c, mode, vemmod, dvmdJ=dtau_ddJ )

      end if

    end subroutine vonmises_stress

  end subroutine rhs_plastic_strain_generic

end module viscoelastic_models_plastic_strain_m

