! Startup from zero stress of flow with a constant velocity gradient
! 2D flow. Combined shear and planar extensional flow.
! Elastoviscoplastic model with Drucker-Prager plasticity (Saramito 2021).
! The integration scheme is second-order Runge Kutta (Heun)

program startup20

  use viscoelastic_models_2D_m

  implicit none

  integer, parameter :: model    = 22, & ! Saramito DP
                        nmodes   = 1, &  ! number of modes
                        flowtype = 0, &  ! 2D
                        ncompc   = 4, &  ! number of components of c
                        ncompt   = 4, &  ! number of components of tau
                        ncompg   = 4, &  ! number of components of gradv
                  numtimesteps   = 1000 ! number of time steps

  real(dp), parameter :: modulus = 10._dp, & ! modulus
                         lambda  = 0.1_dp, & ! relaxation time
                         tau_y   = 0.2_dp, & ! yield stress (cohesion)
                         mu      = 0.2_dp, & ! friction coefficient
                         epsilondot = 0.0_dp, &  ! elongational rate
                         gammadot = 0.1_dp, &  ! shear rate
                         time1 = 200.0_dp, &  ! time for flow stop
                         time2 = 200._dp, &  ! time for flow reversal
                         timestep = 2e-3_dp  ! time step

  integer :: step
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), tau(1,ncompt)
  real(dp) :: k1(1,ncompc,nmodes), gradv(1,ncompg), L(ncompg)
  real(dp) :: trtau, tau_e

  type(vemodel_t) :: vemodel


! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype )

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  vemodel%nonlin(:,1) = [ tau_y, mu ]

! initialize velocity gradient

  L = [ epsilondot, gammadot, 0._dp, -epsilondot ]

  gradv(1,:) = L

! initialize c

  c(1,:,1) = [ 1, 0, 1, 1 ]

! open output file

  open ( unit=13, file='out', recl=600 )

! stepping using RK2

  do step = 1, numtimesteps

    if ( step * timestep > time2 ) then
      gradv(1,:) = - L
    else if ( step * timestep > time1 ) then
      gradv = 0
    end if

!   right-hand side

    call rhs_viscoelastic_2D ( vemodel, gradv, c, rhs )

!   do intermediate step

    k1 = timestep * rhs

!   right-hand side

    call rhs_viscoelastic_2D ( vemodel, gradv, c+k1, rhs )

!   do step

    c = c + ( k1 + timestep * rhs ) / 2

!   compute stress tensor

    call stress_viscoelastic_2D ( vemodel, c, tau )

    trtau = tau(1,1) + tau(1,3) + tau(1,4)
    tau_e = vonmises_2D(tau(1,:))

    write(13,*) step * timestep, tau, tau_e, trtau/3, regime()

  end do

  close(13)

! delete the model

  call delete ( vemodel )

contains

  function regime ()

    integer :: regime

    real(dp) :: chpar

    chpar = mu * trtau - 3 * tau_y

    if ( chpar <= - 3 * tau_e ) then

!     regime I: sticking

      regime = 1

    else if ( chpar >= 2 * mu**2 * tau_e ) then

!     regime III: loosing contact

      regime = 3

    else

!     regime II: sliding

      regime = 2

    end if

  end function regime

end program startup20
