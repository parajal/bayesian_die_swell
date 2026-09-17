! Startup from zero stress of flow with a constant velocity gradient
! The model is a multi-mode PTT-XPP model. 2D model.
! Output is the viscoelastic (polymer) extra-stress.
! The integration scheme is second-order Runge Kutta (Heun)

program startup8_log

  use viscoelastic_models_2D_log_m

  implicit none

  integer, parameter :: nmodes   = 2,  & ! number of modes
                        model    = 18, & ! PTT-XPP model
                        flowtype = 0,  & ! 0: 2D planar 1: 2D axisym.
                        ncompc   = 4,  & ! number of components of c
                        ncompt   = 4,  & ! number of components of tau
                        ncompg   = 4+flowtype ! number of components of gradv

  integer :: flow = 3, &  ! type of flow:
                          ! 1: shear
                          ! 2: planar elongation
                          ! 3: uniaxial elongation
             numtimesteps  = 1000 ! number of time steps

  real(dp)  :: modulus(nmodes) = [2._dp,1._dp], & ! modulus
               lambda(nmodes)  = [3._dp,0.3_dp], & ! relaxation time
               lambda_s(nmodes) = [1._dp,0.1_dp], & ! lambda_s
               nu(nmodes) = [0.2_dp,0.2_dp]    ! nu

  real(dp) :: shearrate = 2.5_dp, &   ! shear rate
              epsilondot = 2.5_dp, &  ! elongational rate
              timestep = 2.5e-2_dp    ! time step

  integer :: mode, step

  real(dp) :: tau(1,ncompt), gradv(1,ncompg)
  real(dp) :: s(1,ncompc,nmodes), rhs(1,ncompc,nmodes), k1(1,ncompc,nmodes)

  type(vemodel_t) :: vemodel


! namelist for input of variables; read from standard input

  namelist /comppar/ flow, modulus, lambda, lambda_s, nu, &
                     shearrate, epsilondot, timestep, numtimesteps

  read ( unit=*, nml=comppar )


! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype )

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  vemodel%nonlin(1,:) = lambda_s
  vemodel%nonlin(2,:) = nu

! initialize velocity gradient

  select case ( flow )
  case(1)
    gradv(1,1:4) = [ 0._dp, shearrate, 0._dp, 0._dp ]
    if ( flowtype == 1 ) print *, 'error flowtype = 1 for flow=1'
  case(2)
    gradv(1,1:4) = [ epsilondot, 0._dp, 0._dp, -epsilondot ]
    if ( flowtype == 1 ) print *, 'error flowtype = 1 for flow=2'
  case(3)
    gradv(1,1:4) = [ epsilondot, 0._dp, 0._dp, -epsilondot/2 ]
    if ( flowtype == 0 ) print *, 'error flowtype = 0 for flow=3'
    gradv(1,ncompg) = -epsilondot/2
  case default
    write(*,'(/a,i0/)') 'Error: wrong value flow: ', flow
    stop
  end select

! initialize s

  do mode = 1, nmodes
    s(1,:,mode) = 0
  end do

! open output file

  open ( unit=13, file='out', recl=300 )

! stepping using RK2

  do step = 1, numtimesteps

!   right-hand side

    call rhs_viscoelastic_2D_log ( vemodel, gradv, s, rhs )

!   do intermediate step

    k1 = timestep * rhs

!   right-hand side

    call rhs_viscoelastic_2D_log ( vemodel, gradv, s+k1, rhs )

!   do step

    s = s + ( k1 + timestep * rhs ) / 2

!   compute stress tensor

    call stress_viscoelastic_2D_log ( vemodel, s, tau )

    write(13,*) step * timestep, tau

  end do

  close(13)

! delete the model

  call delete ( vemodel )

end program startup8_log
