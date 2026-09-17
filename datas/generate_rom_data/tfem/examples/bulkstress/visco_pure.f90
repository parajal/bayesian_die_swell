module visco_pure_m

  use viscoelastic_models_2D_m
  use viscoelastic_models_2D_log_m

  implicit none

contains

  subroutine visco_pure_calc (model, nmodes, ncompc, numtimesteps, &
    modulus, lambda, mobility, shearrate, timestep, flowtype, logc)

  integer, intent(in) :: model, nmodes, ncompc, numtimesteps, flowtype, logc
  real(dp), intent(in) :: modulus, lambda, mobility, shearrate, timestep

  integer, parameter :: ncompg = 4

  integer :: step
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), tau(1,ncompc)
  real(dp) :: gradv(1,ncompg)

  type(vemodel_t) :: vemodel

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype )

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  vemodel%nonlin = mobility

! initialize velocity gradient

  gradv(1,:) = [ 0._dp, shearrate, 0._dp, 0._dp ]

! initialize c

  if ( logc == 1 ) then
    c(1,:,1) = [ 0, 0, 0 ]
  else
    c(1,:,1) = [ 1, 0, 1 ]
  end if

! open output file

  open ( unit=13, file='boundary_c', recl=300 )
  open ( unit=14, file='boundary_tau', recl=300 )

! stepping using explicit Euler

  do step = 1, numtimesteps

!   right-hand side

    if ( logc == 1 ) then
      call rhs_viscoelastic_2D_log ( vemodel, gradv, c, rhs )
    else
      call rhs_viscoelastic_2D ( vemodel, gradv, c, rhs )
    end if

!   do step

    c = c + timestep * rhs

!   compute stress tensor

    if ( logc == 1 ) then
      call stress_viscoelastic_2D_log ( vemodel, c, tau )
    else
      call stress_viscoelastic_2D ( vemodel, c, tau )
    end if

    write(14,*) step * timestep, tau

    write(13,'(1X,3F14.6)') c

  end do

  close(13)

! delete the model

  call delete ( vemodel )

  end subroutine visco_pure_calc

end module visco_pure_m
