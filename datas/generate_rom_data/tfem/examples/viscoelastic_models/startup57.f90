! Startup from zero stress of flow with a constant velocity gradient
! 3D flow. Combined shear and uniaxial extensional flow.
! The model is an EGP elastoviscoplastic model
! Implicit trapezoidal time integration
! Using the contravariant deformation tensor (b).
! 14-mode model with 12 alpha process viscoelastic modes, one elastic mode for
! strain hardening and one mode for the beta process.
! Data for PMMA from van Breemen et al. (2012) DOI: 10.1002/polb.23199

program startup57

  use viscoelastic_models_3D_b_m
  use limits_m, only: USE_CONFORMATION_MODEL, USE_EGP_STRESS_TENSOR_FORM
  use material_data_m

  implicit none

  integer, parameter :: model    = 24, & ! 23: EGP model incompressible
                                         ! 24: EGP model compressible
                        nmodes   = 14, &  ! number of modes
                        am = nmodes-2, & ! 1:am range of modes for alpha process
                        em = nmodes-1, &  ! elastic mode
                        bm = nmodes,   &  ! mode of beta process
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
  real(dp) :: gradv(1,ndim,ndim)
  real(dp) :: smat(n,n), fvec(n), fvecn(n)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemoptn, vemopt

  USE_CONFORMATION_MODEL = .false.
  USE_EGP_STRESS_TENSOR_FORM = .false.

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype, &
    alam_model=alam_model, vmglobal=.true., norlxmode=em )

! fill material parameters

  call fill_data ( vemodel, 'PMMA' )

  if ( model == 24 ) then
!   EGP compressible: input J
    call create_vemopt ( vemoptn, dep_J=.true. )
    vemoptn%J = Jvol
    call create_vemopt ( vemopt, dep_J=.true., compute_drhs_mm=.true., &
      ncomp=ncompb, nmodes=nmodes )
    vemopt%J = Jvol
  else
    call create_vemopt ( vemoptn )
    call create_vemopt ( vemopt, compute_drhs_mm=.true., ncomp=ncompb, &
      nmodes=nmodes )
  end if

! initialize velocity gradient

  gradv(1,1,:) = [ epsilondot,      gammadot,         0._dp ]
  gradv(1,2,:) = [      0._dp, -epsilondot/2,         0._dp ]
  gradv(1,3,:) = [      0._dp,         0._dp, -epsilondot/2 ]

! initialize b

  do m = 1, nmodes
    b(1,:,m) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ]
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

    call rhs_3D_b ( vemoptn )

    fvecn = - reshape ( b + timestep * rhs / 2, [n] )

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

!     right-hand side and Jacobian for all modes at once

      vemopt%drhs_mm = 0 ! initialize all to zero to fill off-diagonal terms

      call rhs_3D_b ( vemopt )

!     solve J^-1 * residual

      smat = - reshape ( timestep * vemopt%drhs_mm / 2, [n,n] )
      do i = 1, n
        smat(i,i) = smat(i,i) + 1
      end do
      fvec = reshape ( b - timestep * rhs / 2, [n] ) + fvecn

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

  end do

  close(13)
  close(14)

! delete the model

  call delete ( vemodel )
  call delete ( vemoptn, vemopt )

contains

! build rhs for given vemopt

  subroutine rhs_3D_b ( vemopt )

    type(vemopt_t), intent(inout) :: vemopt

!   alpha process
    call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs, mode1=1, &
        mode2=am, vemopt=vemopt )

!   elastic mode
    call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs, mode1=em, &
      mode2=em, vemopt=vemopt )

!   beta process
    call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs, mode1=bm, &
        mode2=bm, vemopt=vemopt )

  end subroutine rhs_3D_b

end program startup57
