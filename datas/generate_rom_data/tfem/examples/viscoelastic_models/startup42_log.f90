! Startup from zero stress of flow with a constant velocity gradient
! 3D flow. Combined shear and uniaxial extensional flow.
! The model is an EGP elastoviscoplastic model.
! The integration scheme is second-order Runge Kutta (Heun)

program startup42_log

  use viscoelastic_models_3D_log_m
  use viscoelastic_models_3D_m, only: detc_3D
  use limits_m, only: USE_EGP_STRESS_TENSOR_FORM

  implicit none

  integer, parameter :: model    = 24, & ! 23: EGP model incompressible
                                         ! 24: EGP model compressible
                        nmodes   = 1, &  ! number of modes
                        flowtype = 2, &  ! 3D
                        alam_model = 4, & ! adapted lambda model:
                                          ! 4: Eyring, 5: Ree-Eyring
                        ncompc   = 6, &  ! number of components of c
                        ncompt   = 6, &  ! number of components of tau
                        ndim     = 3, &  ! number of coordinate directions
                  numtimesteps   = 10000 ! number of time steps

  real(dp), parameter :: modulus = 4000, & ! modulus
                         lambda  = 10, & ! relaxation time
                         tau_ref = 500, & ! reference von Mises
                         tau_ref1 = 500, & ! reference von Mises 1 RE
                         tau_ref2 = 2000, & ! reference von Mises 2 RE
                         f1 = 0.6_dp, & ! weight factor for first term RE
                         epsilondot = 0._dp, &  ! elongational rate
                         gammadot = 0.25_dp, &  ! shear rate
                         Lxx_add = 0.01_dp, & ! added value for Lxx (compr.)
                         Jvol = 0.5_dp, &  ! change in volume for model=24
                         time1 = 1.0e4_dp, &  ! time for flow stop
                         time2 = 1.5e4_dp, &  ! time for flow reversal
                         timestep = 1.e-2_dp  ! time step

  integer :: step
  real(dp) :: s(1,ncompc,nmodes), rhs(1,ncompc,nmodes), tau(1,ncompt)
  real(dp) :: k1(1,ncompc,nmodes), gradv(1,ndim,ndim), L(ndim,ndim)
  real(dp) :: c(1,ncompc)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemopt

  USE_EGP_STRESS_TENSOR_FORM = .false.

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype, &
    alam_model=alam_model )

  if ( model == 24 ) then
!   EGP compressible: input J
    call create_vemopt ( vemopt, dep_J=.true. )
    vemopt%J = Jvol
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

! initialize velocity gradient

  L(1,:) = [ epsilondot,      gammadot,         0._dp ]
  L(2,:) = [      0._dp, -epsilondot/2,         0._dp ]
  L(3,:) = [      0._dp,         0._dp, -epsilondot/2 ]

  if ( model == 24 ) L(1,1) = L(1,1) + Lxx_add  ! EGP compressible trL/=0

  gradv(1,:,:) = L

! initialize s

  s(1,:,1) = 0

! open output file

  open ( unit=13, file='out', recl=300 )

! stepping using RK2

  do step = 1, numtimesteps

    if ( step * timestep > time2 ) then
      gradv(1,:,:) = - L
    else if ( step * timestep > time1 ) then
      gradv = 0
    end if

!   right-hand side

    call rhs_viscoelastic_3D_log ( vemodel, gradv, s, rhs, vemopt=vemopt )

!   do intermediate step

    k1 = timestep * rhs

!   right-hand side

    call rhs_viscoelastic_3D_log ( vemodel, gradv, s+k1, rhs, vemopt=vemopt )

!   do step

    s = s + ( k1 + timestep * rhs ) / 2

!   compute stress tensor

    call stress_viscoelastic_3D_log ( vemodel, s, tau, vemopt=vemopt )

    call conformation_3D_log ( vemodel, s(:,:,1), c )

    write(13,*) step * timestep, tau, vonmises_3D(tau(1,:)), detc_3D(c(1,:))

  end do

  close(13)

! delete the model

  call delete ( vemodel )

end program startup42_log
