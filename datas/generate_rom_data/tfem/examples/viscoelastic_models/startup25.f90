! Startup from zero stress of a flow with a constant velocity gradient.
! Multiple continuation with a different velocity gradient to obtain
! the steady state values.
! The model is an evp model with Drucker-Prager plasticity (Saramito 2021).
! Output is the viscoelastic (polymer) extra-stress.
! Both the time-dependent and the final steady data points are output.
! The integration scheme is second-order Runge Kutta (Heun)

program startup25

  use meshgen_basic_common_m
  use viscoelastic_models_3D_m

  implicit none

  integer, parameter :: nmodes  = 1,  & ! number of modes
                        model   = 22, & ! Saramito DP
                        coorsys = 2,  & ! 3D Cartesian
                        ncompc  = 6,  & ! number of components of c
                        ncompt  = 6,  & ! number of components of tau
                        ndim    = 3     ! number of coordinate directions

  integer :: flowtype = 3, &  ! flow type:
                              ! 1: shear
                              ! 2: planar extension
                              ! 3: uniaxial extension
             numdatapoints = 10, & ! number of steady data points
             numtimesteps  = 1000  ! number of time steps in a time interval

  real(dp) :: timestep = 2.5e-2_dp,   &  ! time step
              start_rate = 2.5e-2_dp, &  ! start data point
              end_rate = 2.5e-2_dp,   &  ! end data point
              factor = 2._dp             ! factor last/first data point distance

  real(dp), dimension(nmodes) :: &
               G,        & ! modulus
               lambda,   & ! relaxation time
               tau_y,    & ! yield stress (cohesion)
               mu          ! friction coefficient

  real(dp), allocatable, dimension(:) :: rate ! shear rate (for flowtype=1)
                                              ! strain rate (for flowtype=2,3)

  integer :: mode, step, i

  real(dp) :: time0, tau_e, trtau
  real(dp) :: tau(1,ncompt), gradv(1,ndim,ndim)
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), k1(1,ncompc,nmodes)

  type(vemodel_t) :: vemodel


! namelist for input of variables; read from standard input

  namelist /comppar/ flowtype, G, lambda, tau_y, mu, &
                     numdatapoints, start_rate, end_rate, factor, &
                     timestep, numtimesteps

  read ( unit=*, nml=comppar )

! set data points

  allocate ( rate(numdatapoints) )

  call distribute_elements ( numdatapoints - 1, rate, ratio=1, factor=factor )

  rate = start_rate + rate * (end_rate - start_rate)

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, coorsys )

! set material parameters

  vemodel%modulus = G
  vemodel%lambda = lambda
  vemodel%nonlin(1,:) = tau_y
  vemodel%nonlin(2,:) = mu

! initialize c

  do mode = 1, nmodes
    c(1,:,mode) = [ 1, 0, 0, 1, 0, 1 ]
  end do

! open output files

  open ( unit=13, file='startup.out', recl=300 )
  open ( unit=14, file='steady25.out', recl=300 )

! time stepping using Heun

  time0 = 0

  do i = 1, numdatapoints

    do step = 1, numtimesteps

!     initialize velocity gradient

      select case ( flowtype )
      case(1)
        gradv(1,1,:) = [ 0._dp, rate(i), 0._dp ]
        gradv(1,2,:) = [ 0._dp,   0._dp, 0._dp ]
        gradv(1,3,:) = [ 0._dp,   0._dp, 0._dp ]
      case(2)
        gradv(1,1,:) = [ rate(i),    0._dp, 0._dp ]
        gradv(1,2,:) = [   0._dp, -rate(i), 0._dp ]
        gradv(1,3,:) = [   0._dp,    0._dp, 0._dp ]
      case(3)
        gradv(1,1,:) = [ rate(i),      0._dp,      0._dp ]
        gradv(1,2,:) = [   0._dp, -rate(i)/2,      0._dp ]
        gradv(1,3,:) = [   0._dp,      0._dp, -rate(i)/2 ]
      case default
        write(*,'(/a,i0/)') 'Error: wrong value flowtype: ', flowtype
        stop
      end select

!     right-hand side

      call rhs_viscoelastic_3D ( vemodel, gradv, c, rhs )

!     do intermediate step

      k1 = timestep * rhs

!     right-hand side

      call rhs_viscoelastic_3D ( vemodel, gradv, c+k1, rhs )

!     do step

      c = c + ( k1 + timestep * rhs ) / 2

!     compute stress tensor

      call stress_viscoelastic_3D ( vemodel, c, tau )

      if ( nmodes == 1 ) then

        trtau = tau(1,1) + tau(1,4) + tau(1,6)
        tau_e = vonmises_3D(tau(1,:))

        write(13,'(9es16.8,1x,i0)') time0 + step * timestep, tau, &
          tau_e, trtau/3, regime()

      else

        write(13,'(7es16.8)') time0 + step * timestep, tau

      end if

    end do

    time0 = time0 + numtimesteps * timestep

    if ( nmodes == 1 ) then
      trtau = tau(1,1) + tau(1,4) + tau(1,6)
      tau_e = vonmises_3D(tau(1,:))
      write(14,'(9es16.8,1x,i0)') rate(i), tau, tau_e, trtau/3, regime()
    else
      write(14,'(7es16.8)') rate(i), tau
    end if

  end do

  close(13)
  close(14)

! delete the model

  call delete ( vemodel )

  deallocate ( rate )

contains

  function regime ()

    integer :: regime

    real(dp) :: chpar

    chpar = mu(1) * trtau - 3 * tau_y(1)

    if ( chpar <= - 3 * tau_e ) then

!     regime I: sticking

      regime = 1

    else if ( chpar >= 2 * mu(1)**2 * tau_e ) then

!     regime III: loosing contact

      regime = 3

    else

!     regime II: sliding

      regime = 2

    end if

  end function regime

end program startup25

