! Startup from zero stress of flow with a constant velocity gradient
! Planar extensional flow
! The model is a Giesekus model
! The integration scheme is second-order Runge Kutta (Heun)

program startup2

  use viscoelastic_models_2D_m

  implicit none

  integer, parameter :: model    = 3, &  ! Giesekus model
                        nmodes   = 1, &  ! number of modes
                        flowtype = 0, &  ! 2D
                        ncompc   = 3, &  ! number of components of c
                        ncompt   = 3, &  ! number of components of tau
                        ncompg   = 4, &  ! number of components of gradv
                  numtimesteps   = 1000 ! number of time steps

  real(dp), parameter :: modulus = 2, & ! modulus
                         lambda  = 3, & ! relaxation time
                         mobility = 0.1_dp, &  ! mobility parameter
                         epsilondot = 2.5_dp, &  ! elongational rate
                         timestep = 2.5e-2_dp  ! time step

  integer :: step
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), tau(1,ncompt)
  real(dp) :: k1(1,ncompc,nmodes), gradv(1,ncompg)

  type(vemodel_t) :: vemodel


! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype )

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  vemodel%nonlin = mobility

! initialize velocity gradient

  gradv(1,:) = [ epsilondot, 0._dp, 0._dp, -epsilondot ]

! initialize c

  c(1,:,1) = [ 1, 0, 1 ]

! open output file

  open ( unit=13, file='out', recl=300 )

! stepping using RK2

  do step = 1, numtimesteps

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

    write(13,*) step * timestep, tau

  end do

  close(13)

! delete the model

  call delete ( vemodel )

end program startup2
