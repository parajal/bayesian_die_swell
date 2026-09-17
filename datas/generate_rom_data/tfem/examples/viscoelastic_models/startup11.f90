! Startup from zero stress of flow with a constant velocity gradient
! Planar flow
! The model is a FENE-P model
! The integration scheme is second-order Runge Kutta (Heun)
! Uses the contravariant deformation tensor (b).
! Choice for all four variants of the b-formulation.

program startup11

  use viscoelastic_models_2D_b_m

  implicit none

  integer, parameter :: model    = 20, & ! FENE-P model (choose model=9 or 20)
                        nmodes   = 1,  & ! number of modes
                        flowtype = 0,  & ! 2D
                        bvariant = 1,  & ! b-formulation:
                                         ! 1: CDT 2: symmetric
                                         ! 3: Cholesky 4: Cholesky with log
                        ncompb   = 5,  & ! number of components of b
                        ncompt   = 4,  & ! number of components of tau
                        ncompg   = 4,  & ! number of components of gradv
                  numtimesteps   = 1000 ! number of time steps

  real(dp), parameter :: modulus = 2, & ! modulus
                         lambda  = 3, & ! relaxation time
                         bpar = 50._dp, &  ! b-parameter
                         gammadot = 2.5_dp, &  ! shear rate
                         epsilondot = 0.15_dp, &  ! elongational rate
                         timestep = 2.5e-2_dp  ! time step

  integer :: step
  real(dp) :: b(1,ncompb,nmodes), rhs(1,ncompb,nmodes), tau(1,ncompt)
  real(dp) :: k1(1,ncompb,nmodes), gradv(1,ncompg)

  type(vemodel_t) :: vemodel


! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype, bvariant )

! set material parameters

  vemodel%modulus = modulus
  if ( model == 9 ) then
    vemodel%lambda = lambda*(bpar+3)/bpar     ! lambda_H=lambda*(b+3)/b
  else
    vemodel%lambda = lambda
  end if
  vemodel%nonlin = bpar

! initialize velocity gradient

  gradv(1,:) = [ epsilondot, gammadot, 0._dp, -epsilondot ]

! initialize b

  if ( bvariant == 4 ) then
    b = 0
  else
    b(1,:,1) = [ 1, 0, 0, 1, 1 ]
  end if

! open output file

  open ( unit=13, file='out', recl=300 )
  open ( unit=14, file='outb', recl=300 )

! stepping using RK2

  do step = 1, numtimesteps

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

    write(13,*) step * timestep, tau

    if ( bvariant == 4 ) then
      write(14,*) step * timestep, &
                  exp(b(1,1,1)), b(1,2:3,1), exp(b(1,4,1)), b(1,5,1)
    else
      write(14,*) step * timestep, b
    end if

  end do

  close(13)
  close(14)

! delete the model

  call delete ( vemodel )

end program startup11
