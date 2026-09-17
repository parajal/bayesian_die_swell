! Startup from zero stress of flow with a constant velocity gradient
! 3D flow.
! The model is a Giesekus model
! The integration scheme is second-order Runge Kutta (Heun)
! Using the contravariant deformation tensor (b).
! Choice for all four variants of the b-formulation.

program startup12

  use viscoelastic_models_3D_b_m

  implicit none

  integer, parameter :: model    = 3, & ! Giesekus model
                        nmodes   = 1, & ! number of modes
                        flowtype = 2, & ! 3D
                        bvariant = 1, & ! b-formulation:
                                        ! 1: CDT 2: symmetric
                                        ! 3: Cholesky 4: Cholesky with log
                        ncompb   = 9, & ! number of components of b
                        ncompt   = 6, & ! number of components of tau
                        ndim     = 3, & ! number of coordinate directions
                  numtimesteps   = 20 ! number of time steps

  real(dp), parameter :: modulus = 2, & ! modulus
                         lambda  = 3, & ! relaxation time
                         mobility = 0.1_dp, &  ! mobility parameter
                         gammadot = 2.5_dp, &  ! shear rate
                         epsilondot = 0.15_dp, &  ! elongational rate
                         timestep = 2.5e-2_dp  ! time step

  integer :: step
  real(dp) :: b(1,ncompb,nmodes), rhs(1,ncompb,nmodes), tau(1,ncompt)
  real(dp) :: k1(1,ncompb,nmodes), gradv(1,ndim,ndim)

  type(vemodel_t) :: vemodel


! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype, bvariant )

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  vemodel%nonlin = mobility

! initialize velocity gradient

  gradv(1,1,:) = [ epsilondot,      gammadot,         0._dp ]
  gradv(1,2,:) = [      0._dp, -epsilondot/2,         0._dp ]
  gradv(1,3,:) = [      0._dp,         0._dp, -epsilondot/2 ]

! initialize b

  if ( bvariant == 4 ) then
    b = 0
  else
    b(1,:,1) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ]
  end if

! open output file

!  open ( unit=13, file='out', recl=300 )
!  open ( unit=14, file='outb', recl=300 )

! stepping using RK2

  do step = 1, numtimesteps

!   right-hand side

    call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs )

!   do intermediate step

    k1 = timestep * rhs

!   right-hand side

    call rhs_viscoelastic_3D_b ( vemodel, gradv, b+k1, rhs )

!   do step

    b = b + ( k1 + timestep * rhs ) / 2

!   compute stress tensor

    call stress_viscoelastic_3D_b ( vemodel, b, tau )

    write(*,*) step * timestep, tau

    if ( bvariant == 4 ) then
      write(*,*) step * timestep, exp(b(1,1,1)), b(1,2:4,1), exp(b(1,5,1)), &
        b(1,6:8,1), exp(b(1,9,1))
    else
      write(*,*) step * timestep, b
    end if

  end do

!  close(13)
!  close(14)

! delete the model

  call delete ( vemodel )

end program startup12
