! Startup from zero stress of flow with a constant velocity gradient
! 3D flow. Combined shear and uniaxial extensional flow.
! The model is an elastoviscoplastic model (Saramito 2009).
! The integration scheme is second-order Runge Kutta (Heun)

program startup18

  use viscoelastic_models_3D_m

  implicit none

  integer, parameter :: model    = 2, &  ! 2: Oldroyd model 3: Giesekus
                        nmodes   = 1, &  ! number of modes
                        flowtype = 2, &  ! 3D
                        alam_model = 3, &  ! adapted lambda model
                        ncompc   = 6, &  ! number of components of c
                        ncompt   = 6, &  ! number of components of tau
                        ndim     = 3, &  ! number of coordinate directions
                  numtimesteps   = 1000 ! number of time steps

  real(dp), parameter :: modulus = 1000, & ! modulus
                         lambda  = 0.1, & ! relaxation time
                         mobility = 0.3_dp, & ! mobility parameter
                         tau_y   = 2000, & ! yield stress
                         Kfac    = 100, & ! K viscosity factor
                         nexp    = 0.5, & ! n power-law exponent
                         epsilondot = 0._dp, &  ! elongational rate
                         gammadot = 15._dp, &  ! shear rate
                         time1 = 1.0_dp, &  ! time for flow stop
                         time2 = 1.5_dp, &  ! time for flow reversal
                         timestep = 2.5e-3_dp  ! time step

  integer :: step
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), tau(1,ncompt)
  real(dp) :: k1(1,ncompc,nmodes), gradv(1,ndim,ndim), L(ndim,ndim)

  type(vemodel_t) :: vemodel


! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype, &
    alam_model=alam_model )

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  if ( model == 3 ) vemodel%nonlin(1,:) = mobility
  vemodel%alam(1:3,1) = [ tau_y, Kfac, nexp ]

! initialize velocity gradient

  L(1,:) = [ epsilondot,      gammadot,         0._dp ]
  L(2,:) = [      0._dp, -epsilondot/2,         0._dp ]
  L(3,:) = [      0._dp,         0._dp, -epsilondot/2 ]

  gradv(1,:,:) = L

! initialize c

  c(1,:,1) = [ 1, 0, 0, 1, 0, 1 ]

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

    call rhs_viscoelastic_3D ( vemodel, gradv, c, rhs )

!   do intermediate step

    k1 = timestep * rhs

!   right-hand side

    call rhs_viscoelastic_3D ( vemodel, gradv, c+k1, rhs )

!   do step

    c = c + ( k1 + timestep * rhs ) / 2

!   compute stress tensor

    call stress_viscoelastic_3D ( vemodel, c, tau )

    write(13,*) step * timestep, tau, vonmises_3D(tau(1,:))

  end do

  close(13)

! delete the model

  call delete ( vemodel )

end program startup18
