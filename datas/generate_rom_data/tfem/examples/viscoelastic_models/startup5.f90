! Startup from zero stress of simple shear flow with the xpp model.

program startup5

  use viscoelastic_models_2D_m

  implicit none

  integer, parameter :: model    = 14, & !XPP
                        nmodes   = 1,  & !one mode
                        flowtype = 0,  & !2D flow
                        ncompc   = 5,  & !5 comp conformation tensor
                        ncompt   = 4,  & !4 comp stress tensor
                        ncompg   = 4,  & !4 comp gradv
                        numtimesteps = 10000  ! number of timesteps

  real(dp), parameter :: modulus   = 1._dp, & !modulus
                         lambda    = 1._dp, & !relaxation time
                         !ratio between backbone and arm relaxation
                         rratio    = 150._dp/912._dp, &
                         q         = 5._dp, & !number of arms
                         shearrate = 1.e3_dp,& !shearrate
                         timestep  = 1.e-4_dp ! timestep
  integer :: step
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes),tau(1,ncompt), &
              rhsold(1,ncompc,nmodes)
  real(dp) :: gradv(1,ncompg)
  real(dp) :: lambdas, nu
  type(vemodel_t) :: vemodel

! model definition

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype )

! calculate parameters

  lambdas = rratio*lambda
  nu      = 2._dp/q

! set parameters

  vemodel%modulus = modulus
  vemodel%lambda  = lambda
  vemodel%nonlin(1,1) = lambdas
  vemodel%nonlin(2,1) = nu

! set gradv

  gradv(1,:) = [ 0._dp, shearrate, 0._dp, 0._dp ]

! set initial c

  c(1,:,1) = 1._dp/3._dp*[ 1._dp, 0._dp, 1._dp, 1._dp, 3._dp ]

! open output file

  open (unit=13,file='out',recl=300 )

! integrate using Adams-Bashforth 2nd order with a single Euler-forward
! starting step

  do step =1, numtimesteps

    if ( step >= 2 ) rhsold =rhs

    call rhs_viscoelastic_2D ( vemodel, gradv, c, rhs )

    if( step == 1 ) then
      c = c + timestep * rhs
    else
      c = c + timestep * ( 3._dp/2._dp*rhs - 1._dp/2._dp*rhsold )
    end if

    call stress_viscoelastic_2D ( vemodel, c, tau )

    write(13,*) step*timestep, tau, c(1,1:4,1), c(1,5,1), &
      c(1,1,1)+c(1,3,1)+c(1,4,1)

  end do

  close(13)

  call delete (vemodel)

end program startup5
