! Startup from zero stress of flow with a constant velocity gradient
! 3D flow. Combined shear and uniaxial extensional flow.
! The model is an EGP elastoviscoplastic model
! Implicit trapezoidal time integration
! Using the contravariant deformation tensor (b).
! Similar to startup43, but now using the b-formulation

program startup45

  use viscoelastic_models_3D_b_m
  use limits_m, only: USE_CONFORMATION_MODEL, USE_EGP_STRESS_TENSOR_FORM

  implicit none

  logical, parameter :: reinitialize = .false.

  integer, parameter :: model    = 24, & ! 23: EGP model incompressible
                                         ! 24: EGP model compressible
                        nmodes   = 1, &  ! number of modes
                        flowtype = 2, &  ! 3D
                        alam_model = 4, & ! adapted lambda model:
                                          ! 4: Eyring, 5: Ree-Eyring
                        ncompb   = 9, &  ! number of components of b
                        ncompt   = 6, &  ! number of components of tau
                        ndim     = 3, &  ! number of coordinate directions
                maxnumiterations = 100, &! maximum number of iterations
                  numtimesteps   = 1000 ! number of time steps

  real(dp), parameter :: modulus = 4000, & ! modulus
                         lambda  = 10, & ! relaxation time
                         tau_ref = 500, & ! reference von Mises
                         tau_ref1 = 500, & ! reference von Mises 1 RE
                         tau_ref2 = 2000, & ! reference von Mises 2 RE
                         f1 = 0.6_dp, & ! weight factor for first term RE
                         epsilondot = 0.0_dp, &  ! elongational rate
                         gammadot = 0.25_dp, &   ! shear rate
                         Lxx_add = 0.01_dp, & ! added value for Lxx (compr.))
                         Jvol = 0.5_dp, &  ! change in volume for model=24
                         timestep = 0.1_dp,& ! time step
                              eps = 1.e-10_dp   ! Newton-Raphson accuracy


  integer :: mode, step, iter, i, info, ipiv(ncompb)
  real(dp) :: b(1,ncompb,nmodes), rhs(1,ncompb,nmodes), tau(1,ncompt)
  real(dp) :: gradv(1,ndim,ndim)
  real(dp) :: smat(ncompb,ncompb), fvec(ncompb), fvecn(1,ncompb,nmodes)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemoptn, vemopt

  USE_CONFORMATION_MODEL = .false.
  USE_EGP_STRESS_TENSOR_FORM = .false.

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype, &
    alam_model=alam_model )

  if ( model == 24 ) then
!   EGP compressible: input J
    call create_vemopt ( vemoptn, dep_J=.true. )
    vemoptn%J = Jvol
    call create_vemopt ( vemopt, dep_J=.true., compute_drhs=.true., &
      ncomp=ncompb, nmodes=nmodes )
    vemopt%J = Jvol
  else
    call create_vemopt ( vemoptn )
    call create_vemopt ( vemopt, compute_drhs=.true., ncomp=ncompb, &
      nmodes=nmodes )
  end if

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  if ( alam_model == 4 ) then
    vemodel%alam = tau_ref
  else if ( alam_model == 5 ) then
    vemodel%alam(1,1) = tau_ref1
    vemodel%alam(2,1) = f1
    vemodel%alam(3,1) = tau_ref2
  end if

! initialize velocity gradient

  gradv(1,1,:) = [ epsilondot,      gammadot,         0._dp ]
  gradv(1,2,:) = [      0._dp, -epsilondot/2,         0._dp ]
  gradv(1,3,:) = [      0._dp,         0._dp, -epsilondot/2 ]

! EGP compressible trL/=0
  if ( model == 24 ) gradv(1,1,1) = gradv(1,1,1) + Lxx_add

! initialize b

  do mode = 1, nmodes
    b(1,:,mode) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ]
  end do

! open output file

  open ( unit=13, file='iter.out', recl=300 )
  open ( unit=14, file='out', recl=300 )

! Solve dc/dt = rhs(c) using trapezoidal method:
!
!  cn+1 = cn + deltat * (rhs(cn+1) + rhs(cn)) / 2
!
! or find cn+1 by solving f(c) = 0 using Newton-Raphson, where
!
!  f(c) = c - cn - deltat * ( rhs(c) + rhs(cn) ) / 2
!
! and Jacobian
!
!  df/dc = I - deltat / 2 * drhs(c)/dc

  do step = 1, numtimesteps

!   constant part of right-hand side

    call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs, vemopt=vemoptn )

    fvecn = -b - timestep * rhs / 2

!   loop over all modes

    do mode = 1, nmodes

!     Newton-Raphson for each mode

      iter = 0

      do

        iter = iter + 1

        !print *, iter

        if ( iter > maxnumiterations ) then
          write(*,'(3(a,i0/))') &
            ' Maximum number of iterations reached = ', &
            maxnumiterations, ' step = ', step, ' mode = ', mode
          stop
        end if

!       right-hand side for mode

        call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs, mode, mode, &
          vemopt=vemopt )

!       solve J^-1 * residual

        smat = - timestep * vemopt%drhs(1,:,:,mode) / 2
        do i = 1, ncompb
          smat(i,i) = smat(i,i) + 1
        end do
        fvec = b(1,:,mode) - timestep * rhs(1,:,mode) / 2 + fvecn(1,:,mode)

        call dgesv( ncompb, 1, smat, ncompb, ipiv, fvec, ncompb, info )

        if ( info /= 0 ) then
          print *, 'Lapack DGESV: info = ', info
          stop
        end if

!       update for Newton-Raphson: un+1=un-J(un)^-1*residual(un)

        b(1,:,mode) = b(1,:,mode) - fvec

        !print *, maxval(abs(fvec))
        write(13,*) maxval(abs(fvec))

        if ( maxval(abs(fvec)) < eps ) exit

      end do

    end do

!   compute stress tensor

    call stress_viscoelastic_3D_b ( vemodel, b, tau, vemopt=vemopt )

    write(14,*) step * timestep, tau, vonmises_3D(tau(1,:)), detb_3D(b(1,:,1))

    if ( reinitialize ) call sqrtc_3D_b ( b(:,:,1) )

  end do

  close(13)
  close(14)

! delete the model

  call delete ( vemodel )
  call delete ( vemoptn, vemopt )

end program startup45
