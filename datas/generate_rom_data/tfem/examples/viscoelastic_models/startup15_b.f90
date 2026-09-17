! Startup from zero stress of flow with a constant velocity gradient
! 2D flow. Combined shear and planar extensional flow.
! The model is an elastoviscoplastic model (Saramito 2009).
! The integration scheme is second-order Runge Kutta (Heun)

program startup15_b

  use viscoelastic_models_2D_b_m

  implicit none

  integer, parameter :: model    = 2, &  ! 2: Oldroyd model 3: Giesekus
                        nmodes   = 1, &  ! number of modes
                        flowtype = 0, &  ! 2D
                        alam_model = 3, &  ! adapted lambda model
                        ncompb   = 4, &  ! number of components of b
                        ncompt   = 3, &  ! number of components of tau
                        ncompg   = 4, &  ! number of components of gradv
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
  real(dp) :: b(1,ncompb,nmodes), rhs(1,ncompb,nmodes), tau(1,ncompt)
  real(dp) :: k1(1,ncompb,nmodes), gradv(1,ncompg), L(ncompg)

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

  L = [ epsilondot, gammadot, 0._dp, -epsilondot ]

  gradv(1,:) = L

! initialize b

  b(1,:,1) = [ 1, 0, 0, 1 ]

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

    call rhs_viscoelastic_2D_b ( vemodel, gradv, b, rhs )

!   do intermediate step

    k1 = timestep * rhs

!   right-hand side

    call rhs_viscoelastic_2D_b ( vemodel, gradv, b+k1, rhs )

!   do step

    b = b + ( k1 + timestep * rhs ) / 2

!   compute stress tensor

    call stress_viscoelastic_2D_b ( vemodel, b, tau )

    write(13,*) step * timestep, tau, vonmises_2D(tau(1,:))

  end do

  close(13)

! delete the model

  call delete ( vemodel )

end program startup15_b
