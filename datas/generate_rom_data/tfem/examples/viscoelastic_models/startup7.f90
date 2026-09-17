! Startup from zero stress of flow with a constant velocity gradient
! The model is a multi-mode PTT-XPP model.
! Output is the viscoelastic (polymer) extra-stress.
! The integration scheme is second-order Runge Kutta (Heun)

program startup7

  use viscoelastic_models_3D_m

  implicit none

  integer, parameter :: nmodes   = 2,  & ! number of modes
                        model    = 18, & ! PTT-XPP model
                        flowtype = 2,  & ! 3D
                        ncompc   = 6,  & ! number of components of c
                        ncompt   = 6,  & ! number of components of tau
                        ndim     = 3     ! number of coordinate directions

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

  real(dp) :: tau(1,ncompt), gradv(1,ndim,ndim)
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), k1(1,ncompc,nmodes)

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
    gradv(1,1,:) = [ 0._dp, shearrate, 0._dp ]
    gradv(1,2,:) = [ 0._dp,     0._dp, 0._dp ]
    gradv(1,3,:) = [ 0._dp,     0._dp, 0._dp ]
  case(2)
    gradv(1,1,:) = [ epsilondot,       0._dp,         0._dp ]
    gradv(1,2,:) = [      0._dp, -epsilondot,         0._dp ]
    gradv(1,3,:) = [      0._dp,       0._dp,         0._dp ]
  case(3)
    gradv(1,1,:) = [ epsilondot,         0._dp,         0._dp ]
    gradv(1,2,:) = [      0._dp, -epsilondot/2,         0._dp ]
    gradv(1,3,:) = [      0._dp,         0._dp, -epsilondot/2 ]
  case default
    write(*,'(/a,i0/)') 'Error: wrong value flow: ', flow
    stop
  end select

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

end program startup7
