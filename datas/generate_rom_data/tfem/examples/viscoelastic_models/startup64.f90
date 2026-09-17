! Startup from zero stress of flow with a constant velocity gradient
! 2D flow. Combined shear and uniaxial extensional flow.
! The model is an EGP elastoviscoplastic model (single mode)
! Plastic strain softening.
! The integration scheme is second-order Runge Kutta (Heun)
! Similar to startup58, but now with 2D flow

program startup64

  use viscoelastic_models_2D_m
  use limits_m, only: USE_EGP_STRESS_TENSOR_FORM

  implicit none

  integer, parameter :: model    = 23, & ! 23: EGP model incompressible
                                         ! 24: EGP model compressible
                        nmodes   = 1, &  ! number of modes
                        flowtype = 0, &  ! 2D
                        alam_model = 4, & ! adapted lambda model:
                                          ! 4: Eyring, 5: Ree-Eyring
                        alam_gammap_model = 1, & ! adapted lambda gammap model:
                                          ! 1: EGP strain softening
                        ncompc   = 4, &  ! number of components of c
                        ncompt   = 4, &  ! number of components of tau
                        ncompg   = 4, &  ! number of components of gradv
                  numtimesteps   = 1000 ! number of time steps

  real(dp), parameter :: modulus = 4000, & ! modulus
                         lambda  = 10, & ! relaxation time
                         tau_ref = 500, & ! reference von Mises
                         tau_ref1 = 500, & ! reference von Mises 1 RE
                         tau_ref2 = 2000, & ! reference von Mises 2 RE
                         f1 = 0.6_dp, & ! weight factor for first term RE
                         Sa = 7.4_dp,  & ! state parameter
                         r0 = 0.96_dp, & ! fitting parameter r0
                         r1 = 20.0_dp, & ! fitting parameter r1
                         r2 = -2.0_dp, & ! fitting parameter r2
                         epsilondot = 0._dp, &  ! elongational rate
                         gammadot = 0.25_dp, &  ! shear rate
                         Jvol = 0.5_dp, &  ! change in volume for model=24
                         time1 = 1.0e4_dp, &  ! time for flow stop
                         time2 = 1.5e4_dp, &  ! time for flow reversal
                         timestep = 1.e-2_dp  ! time step

  logical, parameter :: dep_J = model==24

  integer :: step
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), tau(1,ncompt)
  real(dp) :: k1(1,ncompc,nmodes), gradv(1,ncompg), L(ncompg)
  real(dp) :: gammap(1), rhs_gammap(1), k1gp(1)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemopt
  type(vemopt_gammap_t) :: vemopt_gammap

  USE_EGP_STRESS_TENSOR_FORM = .false.

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype, &
    alam_model=alam_model, alam_gammap_model=alam_gammap_model, &
    gammap_mode=1 )

! set optional input parameters

  call create_vemopt ( vemopt, dep_J=dep_J, dep_gammap=.true. )
  call create_vemopt_gammap ( vemopt_gammap, dep_J=dep_J )

  if ( dep_J ) then
    vemopt%J = Jvol
    vemopt_gammap%J = Jvol
  end if

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  if ( alam_model == 4 ) then
    vemodel%alam = tau_ref
  else if ( alam_model == 5 ) then
    vemodel%alam(1,1) = tau_ref1
    vemodel%alam(2,1) = f1
    vemodel%alam(3,1) = tau_ref2
  end if
  vemodel%alam_gammap(1,1) = Sa
  vemodel%alam_gammap(2,1) = r0
  vemodel%alam_gammap(3,1) = r1
  vemodel%alam_gammap(4,1) = r2

! initialize velocity gradient

  L = [ epsilondot, gammadot, 0._dp, -epsilondot ]

  gradv(1,:) = L

! initialize c

  c(1,:,1) = [ 1, 0, 1, 1 ]

! initialize gammap

  gammap(1) = 0
  vemopt%gammap = gammap

! open output file

  open ( unit=13, file='out', recl=300 )

! stepping using RK2

  do step = 1, numtimesteps

    if ( step * timestep > time2 ) then
      gradv(1,:) = - L
    else if ( step * timestep > time1 ) then
      gradv = 0
    end if

!   right-hand side

    call rhs_viscoelastic_2D ( vemodel, gradv, c, rhs, vemopt=vemopt )
    call rhs_plastic_strain_2D ( vemodel, c, gammap, rhs_gammap, &
      vemopt_gammap=vemopt_gammap )

!   do intermediate step

    k1 = timestep * rhs
    k1gp = timestep * rhs_gammap

!   right-hand side

    vemopt%gammap = gammap + k1gp
    call rhs_viscoelastic_2D ( vemodel, gradv, c+k1, rhs, vemopt=vemopt )
    call rhs_plastic_strain_2D ( vemodel, c, gammap + k1gp, rhs_gammap, &
      vemopt_gammap=vemopt_gammap )

!   do step

    c = c + ( k1 + timestep * rhs ) / 2
    gammap = gammap + ( k1gp + timestep * rhs_gammap ) / 2
    vemopt%gammap = gammap

!   compute stress tensor

    call stress_viscoelastic_2D ( vemodel, c, tau, vemopt=vemopt )

    write(13,*) step * timestep, tau, vonmises_2D(tau(1,:)), &
      detc_2D(c(1,:,1)), gammap(1)

  end do

  close(13)

! delete the model

  call delete ( vemodel )
  call delete ( vemopt )
  call delete ( vemopt_gammap )

end program startup64
