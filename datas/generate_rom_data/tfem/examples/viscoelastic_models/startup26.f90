! Startup from zero stress of flow with an imposed shear stress function.
! Instantaneous elasticity has been handled with a finite but fast ramp
! towards the constant stress value.
! The model is a Giesekus model with multiple modes. 2D flow.
! The integration scheme is second-order Runge Kutta (Heun)

program startup26

  use viscoelastic_models_2D_m

  implicit none

  integer, parameter :: model    = 3, &  ! 2: Oldroyd model 3: Giesekus
                        nmodes   = 2, &  ! number of modes
                        flowtype = 0, &  ! 2D
                        ncompc   = 3, &  ! number of components of c
                        ncompt   = 3, &  ! number of components of tau
                        ncompg   = 4, &  ! number of components of L
                  numtimesteps1  = 10, & ! number of time steps in ramp
                  numtimesteps2  = 1000  ! number of time steps after ramp

  real(dp), parameter :: modulus(nmodes) = [2,1], & ! modulus
                         lambda(nmodes)  = [3._dp,0.3_dp], & ! relaxation time
                         mobility(nmodes) = 0.1_dp, &  ! mobility parameter
                         time1 = 1.0e-3_dp, &  ! time shear stress ramp
                         tauxy_imp = 2, &  ! shear stress imposed
                         timestep1 = time1/numtimesteps1, &! time step in ramp
                         timestep2 = 1.0e-3_dp   ! time step after ramp

  integer :: step, i, zone, mode
  real(dp) :: gdot1, gdot2, shearstrain, timestep, tn, tnp1
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), tau(1,ncompt)
  real(dp) :: k1(1,ncompc,nmodes), gradv(1,ncompg), L(ncompg)

  type(vemodel_t) :: vemodel


! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype )

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  if ( model == 3 ) vemodel%nonlin(1,:) = mobility

! initialize shearstrain

  shearstrain = 0

! initialize c

  do mode = 1, nmodes
    c(1,:,mode) = [ 1, 0, 1 ]
  end do

! open output file

  open ( unit=13, file='out', recl=300 )

! stepping using Heun

  zone = 1
  timestep = timestep1
  tn = 0
  write(13,*) (0._dp,i=1,6)

  do step = 1, numtimesteps1 + numtimesteps2

    if ( step == numtimesteps1 + 1 ) then
      timestep = timestep2
      zone = 2
    end if

    tnp1 = tn + timestep

!   initial step

    gdot1 = gammadot ( zone, tn, c )

!   velocity gradient

    L = [ 0._dp, gdot1, 0._dp, 0._dp ]

    gradv(1,:) = L

!   right-hand side

    call rhs_viscoelastic_2D ( vemodel, gradv, c, rhs )

!   do intermediate step

    k1 = timestep * rhs

    gdot2 = gammadot ( zone, tnp1, c+k1 )

!   velocity gradient

    L = [ 0._dp, gdot2, 0._dp, 0._dp ]

    gradv(1,:) = L

!   right-hand side

    call rhs_viscoelastic_2D ( vemodel, gradv, c+k1, rhs )

!   do step

    shearstrain = shearstrain + timestep * ( gdot1 + gdot2 ) / 2
    c = c + ( k1 + timestep * rhs ) / 2

!   compute stress tensor

    call stress_viscoelastic_2D ( vemodel, c, tau )

    write(13,*) tnp1, shearstrain, tau, abs(tau(1,2)-tauxy(zone,tnp1))

    tn = tnp1

  end do

  close(13)

! delete the model

  call delete ( vemodel )

contains


! shear rate function

  function gammadot ( zone, t, c )

    integer, intent(in) :: zone
    real(dp), intent(in) :: t
    real(dp), intent(in) :: c(:,:,:)

    real(dp) :: gammadot
    real(dp) :: sumfxy, sumcyy

    call rhs_viscoelastic_relax_2D ( vemodel, c, rhs )

    sumfxy = - sum( modulus * rhs(1,2,:) )
    sumcyy = sum( modulus * c(1,3,:) )

    gammadot = ( tauxydot(zone,t) + sumfxy ) / sumcyy

  end function gammadot


! stress rate as a function of time (imposed)

  function tauxydot ( zone, t )

    integer, intent(in) :: zone
    real(dp), intent(in) :: t

    real(dp) :: tauxydot

    select case (zone)
    case(1)
      tauxydot = tauxy_imp / time1
    case(2)
      tauxydot = 0
    case default
      stop 'tauxydot: zone not defined'
    end select

  end function tauxydot


! stress as a function of time (imposed, but only for computing actual error)

  function tauxy ( zone, t )

    integer, intent(in) :: zone
    real(dp), intent(in) :: t

    real(dp) :: tauxy

    select case (zone)
    case(1)
      tauxy = tauxy_imp / time1 * t
    case(2)
      tauxy = tauxy_imp
    case default
      stop 'tauxy: zone not defined'
    end select

  end function tauxy

end program startup26
