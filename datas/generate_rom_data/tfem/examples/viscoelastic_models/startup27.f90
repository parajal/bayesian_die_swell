! Startup from zero stress of flow with an imposed shear stress function.
! Instantaneous elasticity has been handled with a finite but fast ramp
! towards the constant stress value.
! The model is an elastoviscoplastic model: SRM1 (Saramito 2007) or SRM2
! (Saramito 2009). Single mode. 2D flow.
! The integration scheme is second-order Runge Kutta (Heun)

program startup27

  use viscoelastic_models_2D_m

  implicit none

  integer, parameter :: model    = 2, &  ! 2: Oldroyd model 3: Giesekus
                        nmodes   = 1, &  ! number of modes
                        flowtype = 0, &  ! 2D
                        alam_model = 3, &  ! adapted lambda model
!                         0: constant lambda (not adapted)
!                         1  elastic model 1/lambda=0 (no aditional parameters)
!                         2  elastoviscoplastic SRM1, parameter: tau_y
!                         3  elastoviscoplastic SRM2, parameters: tau_y, K, n
                        ncompc   = 3, &  ! number of components of c
                        ncompt   = 3, &  ! number of components of tau
                        ncompg   = 4, &  ! number of components of L
                  numtimesteps1  = 10, & ! number of time steps in ramp
                  numtimesteps2  = 1000  ! number of time steps after ramp

  real(dp), parameter :: modulus = 10000, & ! modulus
                         lambda  = 0.1, & ! relaxation time
                         mobility = 0.3_dp, & ! mobility parameter
                         tau_y   = 2000, & ! yield stress
                         Kfac    = 100, & ! K viscosity factor
                         nexp    = 0.5, & ! n power-law exponent
                         time1 = 1.0e-3_dp, &  ! time shear stress ramp
                         tauxy_imp = 2100, &  ! shear stress imposed
                         timestep1 = time1/numtimesteps1, &! time step in ramp
                         timestep2 = 1.0e-3_dp   ! time step after ramp

  integer :: step, i, zone
  real(dp) :: gdot1, gdot2, shearstrain, timestep, tn, tnp1
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), tau(1,ncompt)
  real(dp) :: k1(1,ncompc,nmodes), gradv(1,ncompg), L(ncompg)

  type(vemodel_t) :: vemodel


! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype, &
    alam_model=alam_model )

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  if ( model == 3 ) vemodel%nonlin(1,:) = mobility
  select case ( alam_model )
    case(2)
      vemodel%alam(1:1,1) = [ tau_y ]
    case(3)
      vemodel%alam(1:3,1) = [ tau_y, Kfac, nexp ]
    case default
      write(*,'(/a,i0/)') 'Error: wrong value alam_model: ', alam_model
      stop
  end select

! initialize shearstrain

  shearstrain = 0

! initialize c

  c(1,:,1) = [ 1, 0, 1 ]

! open output file

  open ( unit=13, file='out', recl=300 )

! stepping using Heun

  zone = 1
  timestep = timestep1
  tn = 0
  write(13,*) (0._dp,i=1,7)

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

    write(13,*) tnp1, shearstrain, tau, vonmises_2D(tau(1,:)), &
                abs(tau(1,2)-tauxy(zone,tnp1))

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

end program startup27
