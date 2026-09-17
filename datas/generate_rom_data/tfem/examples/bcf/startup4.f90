! Startup from zero stress of flow with a constant velocity gradient
! Simple shear flow
! The model is a Hookean dumbbell (stochastic model)
! The integration scheme is implicit Euler

program startup4

  use stochastic_models_2D_m

  implicit none

  integer, parameter :: model   = 1, &  ! Hookean dumbbell
                        coorsys = 0, &  ! 2D Cartesian
                        ncompq  = 2, &  ! number of components of Q
                        ncompt  = 3, &  ! number of components of tau
                        ncompg  = 4, &  ! number of components of gradv
                        nfield  = 2000, &  ! number of Q fields
                        timeint = 2, &  ! implicit Euler
                eqnumtimesteps  = 2000, & ! number of time steps to equilibrate
                  numtimesteps  = 2500 ! number of time steps

  real(dp), parameter :: modulus = 1, & ! modulus
                         lambda  = 1, & ! relaxation time
                         shearrate = 1.0_dp, &  ! shearrate
                         timestep = 1.e-2_dp  ! time step

  integer :: step
  real(dp) :: qn(nfield,ncompq), qnp1(nfield,ncompq), brownf(nfield,ncompq)
  real(dp) :: tau(ncompt)
  real(dp) :: eqgradv(ncompg), gradv(ncompg)

  type(stmodel_t) :: stmodel
  type(stnumpar_t) :: stnumpar


! define the model

  call create_stochastic_model ( model, stmodel, coorsys )

! set material parameters

  stmodel%modulus = modulus
  stmodel%lambda = lambda

! initialize velocity gradient

  gradv(:) = [ 0._dp, shearrate, 0._dp, 0._dp ]

! initialize q

  qn = 0

! open output file

  open ( unit=13, file='out', recl=300 )

! stepping

  stnumpar%nfield = nfield
  stnumpar%timeint = timeint
  stnumpar%tstep = timestep

! equilibrate

  eqgradv = 0

  do step = 1, eqnumtimesteps

!   fill Brownian vector with uniform distribution, variance 1

    call random_number ( brownf )

    brownf = ( 2 * brownf - 1 ) * sqrt(3._dp)

!   step Q

    call ststep1_2D ( stmodel, stnumpar, eqgradv, qn, qnp1, brownf )

    qn = qnp1

  end do

! real stepping

  do step = 1, numtimesteps

!   fill Brownian vector with uniform distribution, variance 1

    call random_number ( brownf )

    brownf = ( 2 * brownf - 1 ) * sqrt(3._dp)

!   step Q

    call ststep1_2D ( stmodel, stnumpar, gradv, qn, qnp1, brownf )

!   compute stress tensor

    call stress_stochastic1_2D ( stmodel, qnp1, tau )
    !call structure_tensor_2D ( stmodel, qnp1, c )

    qn = qnp1

    write(13,*) step * timestep, tau
    !write(13,*) step * timestep, c

  end do

  close(13)

! delete the model

  call delete ( stmodel )

end program startup4
