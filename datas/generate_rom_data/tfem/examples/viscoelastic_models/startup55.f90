! Startup from zero stress of flow with a constant velocity gradient
! 3D flow. Combined shear and uniaxial extensional flow.
! The model is an EGP elastoviscoplastic model
! Implicit BDF2 time integration
! Using the contravariant deformation tensor (b).
! 13-mode model with 12 viscoelastic modes and one elastic mode for
! strain hardening.
! Data for PMMA from van Breemen et al. (2012) DOI: 10.1002/polb.23199
! Similar to startup54, but now using the b-formulation

program startup55

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
                maxnumiterations = 100, &! maximum number of iterations
                  numtimesteps   = 100 ! number of time steps

  integer, parameter :: n = ncompb*nmodes ! number of equations to solve

  real(dp), parameter :: epsilondot = 0.0_dp, &  ! elongational rate
                         gammadot = 1.e-2_dp, &   ! shear rate
                         Jvol = 1.0_dp, &  ! change in volume for model=24
                         timestep = 0.1_dp, & ! time step
                         eps = 1.e-10_dp   ! Newton-Raphson accuracy


  integer :: m, step, iter, i, info, ipiv(n)
  real(dp) :: b(1,ncompb,nmodes), rhs(1,ncompb,nmodes), tau(1,ncompt)
  real(dp) :: bnm1(1,ncompb,nmodes)
  real(dp) :: gradv(1,ndim,ndim)
  real(dp) :: smat(n,n), fvec(n), fvecn(n)

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
    call create_vemopt ( vemopt, dep_J=.true., compute_drhs_mm=.true., &
      ncomp=ncompb, nmodes=nmodes )
    vemopt%J = Jvol
  else
    call create_vemopt ( vemopt, compute_drhs_mm=.true., ncomp=ncompb, &
      nmodes=nmodes )
  end if

! initialize velocity gradient

  gradv(1,1,:) = [ epsilondot,      gammadot,         0._dp ]
  gradv(1,2,:) = [      0._dp, -epsilondot/2,         0._dp ]
  gradv(1,3,:) = [      0._dp,         0._dp, -epsilondot/2 ]

! initialize c

  do m = 1, nmodes
    b(1,:,m) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ]
  end do

! open output file

  open ( unit=13, file='iter.out', recl=300 )
  open ( unit=14, file='out', recl=300 )

! Solve dc/dt = rhs(c) using BDF2 method:
!
!  step 1:
!    cn+1 = cn + deltat * rhs(cn+1)
!  step >1:
!    1.5 * cn+1 = 2 * cn - 0.5 * cn-1 + deltat * rhs(cn+1)
!
! or find cn+1 by solving f(c) = 0 using Newton-Raphson, where
!
!  step 1:
!    f(c) = c - cn - deltat * rhs(c)
!  step >1:
!    f(c) = 1.5 * c - 2 * cn + 0.5 * cn-1 - deltat * rhs(c)
!
! and Jacobian
!
!  step 1:
!    df/dc = I - deltat * drhs(c)/dc
!  step >1:
!    df/dc = 1.5 * I - deltat * drhs(c)/dc

  do step = 1, numtimesteps

!   constant part of right-hand side

    if ( step == 1 ) then
      fvecn = - reshape ( b, [n] )
    else
      fvecn = - reshape ( 2 * b - 0.5_dp * bnm1, [n] )
    end if

!   copy bn to bn-1

    bnm1 = b

!   Newton-Raphson for all modes at once

    iter = 0

    do

      iter = iter + 1

      !print *, iter

      if ( iter > maxnumiterations ) then
        write(*,'(2(a,i0/))') &
          ' Maximum number of iterations reached = ', &
          maxnumiterations, ' step = ', step
        stop
      end if

!     right-hand side for mode

      vemopt%drhs_mm = 0 ! initialize all to zero to fill off-diagonal terms

      call rhs_3D_b

!     solve J^-1 * residual

      if ( step == 1 ) then

        smat = - reshape ( timestep * vemopt%drhs_mm, [n,n] )
        do i = 1, n
          smat(i,i) = smat(i,i) + 1
        end do
        fvec = reshape ( b - timestep * rhs, [n] ) + fvecn

      else

        smat = - reshape ( timestep * vemopt%drhs_mm, [n,n] )
        do i = 1, n
          smat(i,i) = smat(i,i) + 1.5_dp
        end do
        fvec = reshape ( 1.5_dp * b - timestep * rhs, [n] ) + fvecn

      end if

      call dgesv( n, 1, smat, n, ipiv, fvec, n, info )

      if ( info /= 0 ) then
        print *, 'Lapack DGESV: info = ', info
        stop
      end if

!     update for Newton-Raphson: un+1=un-J(un)^-1*residual(un)

      b = b - reshape ( fvec, [1,ncompb,nmodes] )

      !print *, maxval(abs(fvec))
      write(13,*) maxval(abs(fvec))

      if ( maxval(abs(fvec)) < eps ) exit

    end do

!   compute stress tensor

    call stress_viscoelastic_3D_b ( vemodel, b, tau, vemopt=vemopt )

    write(14,*) step * timestep, tau, vonmises_3D(tau(1,:)), detb_3D(b(1,:,1))

    if ( reinitialize ) then
      do m = 1, nmodes
        call sqrtc_3D_b ( b(:,:,m) )
      end do
    end if

  end do

  close(13)
  close(14)

! delete the model

  call delete ( vemodel )
  call delete ( vemopt )

contains

! build rhs and matrix

  subroutine rhs_3D_b

!   alpha process
    call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs, vemopt=vemopt, &
      mode1=1, mode2=am )

!   elastic mode
    call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs, vemopt=vemopt, &
      mode1=em, mode2=em )

  end subroutine rhs_3D_b

end program startup55
