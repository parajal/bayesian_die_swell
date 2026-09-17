! Startup from zero stress of flow with a constant velocity gradient
! 3D flow. Combined shear and uniaxial extensional flow.
! The model is an EGP elastoviscoplastic model.
! The integration scheme is second-order Runge Kutta (Heun)
! Using the contravariant deformation tensor (b).
! 13-mode model with 12 viscoelastic modes and one elastic mode for
! strain hardening.
! Data for PMMA from van Breemen et al. (2012) DOI: 10.1002/polb.23199
! Similar to startup50, but now using the b-formulation

program startup51

  use viscoelastic_models_3D_b_m
  use limits_m, only: USE_CONFORMATION_MODEL, USE_EGP_STRESS_TENSOR_FORM
  use material_data_m

  implicit none

  logical, parameter :: reinitialize = .false.

  integer, parameter :: model    = 24, & ! 23: EGP model incompressible
                                         ! 24: EGP model compressible
                        nmodes   = 13, &  ! number of modes
                        am = nmodes-1, & ! 1:am range of modes for alpha process
                        em = nmodes,  &  ! elastic mode
                        flowtype = 2, &  ! 3D
                        alam_model = 4, & ! adapted lambda model:
                                          ! 4: Eyring
                        ncompb   = 9, &  ! number of components of b
                        ncompt   = 6, &  ! number of components of tau
                        ndim     = 3, &  ! number of coordinate directions
                  numtimesteps   = 1000000 ! number of time steps

  real(dp), parameter :: epsilondot = 0._dp, &  ! elongational rate
                         gammadot = 1e-2_dp, &  ! shear rate
                         Jvol = 1.0_dp, &  ! change in volume for model=24
                         time1 = 1.0e4_dp, &  ! time for flow stop
                         time2 = 1.5e4_dp, &  ! time for flow reversal
                         timestep = 1.e-5_dp  ! time step

  integer :: step, m

  real(dp) :: b(1,ncompb,nmodes), rhs(1,ncompb,nmodes), tau(1,ncompt)
  real(dp) :: k1(1,ncompb,nmodes), gradv(1,ndim,ndim), L(ndim,ndim)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemopt

  USE_CONFORMATION_MODEL = .false.
  USE_EGP_STRESS_TENSOR_FORM = .false.

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype, &
    alam_model=alam_model, vmglobal=.true., norlxmode=em )

! fill material parameters

  call fill_data ( vemodel, 'PMMA' )

  if ( model == 24 ) then
!   EGP compressible: input J
    call create_vemopt ( vemopt, dep_J=.true. )
    vemopt%J = Jvol
  else
    call create_vemopt ( vemopt )
  end if

! initialize velocity gradient

  L(1,:) = [ epsilondot,      gammadot,         0._dp ]
  L(2,:) = [      0._dp, -epsilondot/2,         0._dp ]
  L(3,:) = [      0._dp,         0._dp, -epsilondot/2 ]

  gradv(1,:,:) = L

! initialize b

  do m = 1, nmodes
    b(1,:,m) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ]
  end do

! open output file

  open ( unit=13, file='out', recl=300 )

! stepping using RK2

  do step = 1, numtimesteps

    if ( step * timestep > time2 ) then
      gradv(1,:,:) = - L
    else if ( step * timestep > time1 ) then
      gradv = 0
    end if

!   right-hand side

    call rhs_3D_b ( b )

!   do intermediate step

    k1 = timestep * rhs

!   right-hand side

    call rhs_3D_b ( b+k1 )

!   do step

    b = b + ( k1 + timestep * rhs ) / 2

!   compute stress tensor

    call stress_viscoelastic_3D_b ( vemodel, b, tau, vemopt=vemopt )

    write(13,*) step * timestep, tau, vonmises_3D(tau(1,:))
!    write(13,*) step * timestep, tau, vonmises_3D(tau(1,:)), &
!      (detb_3D(b(1,:,m)), m=1,nmodes)

    if ( reinitialize ) then
      do m = 1, nmodes
        call sqrtc_3D_b ( b(:,:,m) )
      end do
    end if

  end do

  close(13)

! delete the model

  call delete ( vemodel )

contains

! build rhs for given b

  subroutine rhs_3D_b ( b )

    real(dp), dimension(:,:,:), intent(in) :: b

!   alpha process
    call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs, vemopt=vemopt, &
      mode1=1, mode2=am )

!   elastic mode
    call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs, vemopt=vemopt, &
      mode1=em, mode2=em )

  end subroutine rhs_3D_b

end program startup51
