! Startup from zero stress of flow with a constant velocity gradient
! 3D flow. Combined shear and uniaxial extensional flow.
! The model is a Giesekus model, two-modes
! The integration scheme is second-order Runge Kutta (Heun)

program startup4

  use viscoelastic_models_3D_m

  implicit none

  integer, parameter :: model    = 3, &  ! Giesekus model
                        nmodes   = 2, &  ! number of modes
                        flowtype = 2, &  ! 3D
                        ncompc   = 6, &  ! number of components of c
                        ncompt   = 6, &  ! number of components of tau
                        ndim     = 3, &  ! number of coordinate directions
                  numtimesteps   = 1000 ! number of time steps

  real(dp), parameter :: modulus(nmodes) = [2,1], & ! modulus
                         lambda(nmodes)  = [3._dp,0.3_dp], & ! relaxation time
                         mobility(nmodes) = 0.1_dp, &  ! mobility parameter
                         epsilondot = 2.5_dp, &  ! elongational rate
                         gammadot = 1.5_dp, &  ! shear rate
                         timestep = 2.5e-2_dp  ! time step

  integer :: mode, step
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), tau(1,ncompt)
  real(dp) :: k1(1,ncompc,nmodes), gradv(1,ndim,ndim)

  type(vemodel_t) :: vemodel


! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype )

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  vemodel%nonlin(1,:) = mobility

! initialize velocity gradient

  gradv(1,1,:) = [ epsilondot,      gammadot,         0._dp ]
  gradv(1,2,:) = [      0._dp, -epsilondot/2,         0._dp ]
  gradv(1,3,:) = [      0._dp,         0._dp, -epsilondot/2 ]

! initialize c

  do mode = 1, nmodes
    c(1,:,mode) = [ 1, 0, 0, 1, 0, 1 ]
  end do

! open output file

  open ( unit=13, file='out', recl=300 )

! stepping using RK2

  do step = 1, numtimesteps

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

    write(13,*) step * timestep, tau

  end do

  close(13)

! delete the model

  call delete ( vemodel )

end program startup4
