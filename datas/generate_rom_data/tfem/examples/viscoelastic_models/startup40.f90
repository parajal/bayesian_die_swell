! Startup from zero stress of flow with a constant velocity gradient
! 3D flow. Combined shear and uniaxial extensional flow.
! The model is a Maxwell, Giesekus of FENE-P model, possible multiple modes.
! The time is split into three intervals. The integration scheme in the
! first interval is either the implicit Euler or the trapezoidal method.
! In the second and third interval there is an additional choice for BDF2.
! In this way it is possible to start with small time steps and go to large
! ones in the third interval to obtain a steady state quickly.
! Using the contravariant deformation tensor (b).
!
! IMPLICIT EULER:
!
! Solve db/dt = rhs(b):
!
!  bn+1 = bn + deltat * rhs(bn+1)
!
! or find bn+1 by solving f(b) = 0 using Newton-Raphson, where
!
!  f(b) = b - bn - deltat * rhs(b)
!
! and Jacobian
!
!  df/db = I - deltat * drhs(b)/db
!
! TRAPEZOIDAL:
!
! Solve db/dt = rhs(b):
!
!  bn+1 = bn + deltat * (rhs(bn+1) + rhs(bn)) / 2
!
! or find bn+1 by solving f(b) = 0 using Newton-Raphson, where
!
!  f(b) = b - bn - deltat * ( rhs(b) + rhs(bn) ) / 2
!
! and Jacobian
!
!  df/db = I - deltat / 2 * drhs(b)/db
!
! BDF2:
!
! Solve db/dt = rhs(b):
!
!  gamma0*bn+1 - alpah0*bn - alpha1*bnm1 = deltat * rhs(bn+1)
!
! or find bn+1 by solving f(b) = 0 using Newton-Raphson, where
!
!  f(b) = gamma0*b - alpha0*bn - alpha1*bnm1 - deltat * rhs(b)
!
! and Jacobian
!
!  df/db = gamma0*I - deltat * drhs(b)/db


program startup40

  use viscoelastic_models_3D_b_m
!  use limits_m, only: USE_CONFORMATION_MODEL

  implicit none

  integer, parameter :: model    = 3, & ! 2: Maxwell 3: Giesekus 20: FENE-P
                        nmodes   = 1,  & ! number of modes
                        flowtype = 2,  & ! 3D
                        ncompb   = 9,  & ! number of components of b
                        ncompt   = 6,  & ! number of components of tau
                        ndim     = 3,  & ! number of coordinate directions
                maxnumiterations = 100, &! maximum number of iterations
                        method1  = 2,  & ! 1: Euler 2: trapezoidal
                        method2  = 3,  & ! 1: Euler 2: trapezoidal 3: BDF2
                        method3  = 3,  & ! 1: Euler 2: trapezoidal 3: BDF2
                  numtimesteps1  = 1,  & ! number of time steps first interval
                  numtimesteps2  = 12,  & ! number of time steps second interval
                  numtimesteps3  = 80  ! number of time steps third interval

  real(dp), parameter :: modulus(nmodes) = [1], & ! modulus
                         lambda(nmodes)  = [1._dp], & ! relaxation time
                         mobility(nmodes)  = 0.1_dp, & ! alpha
                         bpar(nmodes)  = 100.0_dp, & ! bparameter
                         epsilondot = 0.0_dp, &  ! elongational rate
                         gammadot = 10.0_dp, &   ! shear rate
                         timestep1 = 1.0e-1_dp,& ! time step interval 1
                         timestep2 = 2.0e-1_dp,& ! time step interval 2
                         timestep3 = 2.0e-1_dp,& ! time step interval 3
                              eps = 1.e-10_dp   ! Newton-Raphson accuracy


  integer :: mode, step, iter, i, info, ipiv(ncompb), method
  real(dp) :: timestep, time, r, gamma0, alpha0, alpha1
  real(dp) :: b(1,ncompb,nmodes), rhs(1,ncompb,nmodes), tau(1,ncompt)
  real(dp) :: bn(1,ncompb,nmodes), bnm1(1,ncompb,nmodes)
  real(dp) :: gradv(1,ndim,ndim)
  real(dp) :: smat(ncompb,ncompb), fvec(ncompb), fvecn(1,ncompb,nmodes)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemopt

!  USE_CONFORMATION_MODEL=.true.

! create vemopt structure

  call create ( vemopt, compute_drhs=.true., ncomp=ncompb, nmodes=nmodes )

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype )

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  if ( model == 3 ) then
    vemodel%nonlin(1,:) = mobility
  else if ( model == 20 ) then
    vemodel%nonlin(1,:) = bpar
  end if

! initialize velocity gradient

  gradv(1,1,:) = [ epsilondot,      gammadot,         0._dp ]
  gradv(1,2,:) = [      0._dp, -epsilondot/2,         0._dp ]
  gradv(1,3,:) = [      0._dp,         0._dp, -epsilondot/2 ]

! initialize b

  do mode = 1, nmodes
    b(1,:,mode) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ]
  end do

  bn = b
  bnm1 = b

! open output file

  open ( unit=13, file='iter.out', recl=300 )
  open ( unit=14, file='out', recl=300 )
  open ( unit=15, file='outb', recl=300 )

! start stepping

  time = 0

  do step = 1, numtimesteps1 + numtimesteps2 + numtimesteps3

!   update b in previous steps for new time step

    bnm1 = bn
    bn = b

    if ( step <= numtimesteps1 ) then
      method = method1
      timestep = timestep1
    else if ( step <= numtimesteps1+numtimesteps2 ) then
      method = method2
      r = timestep2/timestep
      timestep = timestep2
    else
      method = method3
      r = timestep3/timestep
      timestep = timestep3
    end if

!   constant part of right-hand side

    select case ( method )
    case(1) ! Implicit Euler
      fvecn = - bn
    case(2) ! Trapezoidal
      call rhs_viscoelastic_3D_b ( vemodel, gradv, bn, rhs )
      fvecn = - bn - timestep * rhs / 2
    case(3) ! BDF2
      gamma0 = (1+2*r)/(1+r)
      alpha0 = 1+r
      alpha1 = -r**2/(1+r)
      fvecn = - alpha0 * bn - alpha1 * bnm1
    case default
      write(*,'(/a,i0/)') 'Error: wrong value method: ', method
      stop
    end select

    time = time + timestep

!   loop over all modes

    do mode = 1, nmodes

!     Newton-Raphson for each mode

      iter = 0

      do

        iter = iter + 1

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

        select case ( method )
        case(1) ! Implicit Euler
          smat = - timestep * vemopt%drhs(1,:,:,mode)
          do i = 1, ncompb
            smat(i,i) = smat(i,i) + 1
          end do
          fvec = b(1,:,mode) - timestep * rhs(1,:,mode) + fvecn(1,:,mode)
        case(2) ! Trapezoidal
          smat = - timestep * vemopt%drhs(1,:,:,mode) / 2
          do i = 1, ncompb
            smat(i,i) = smat(i,i) + 1
          end do
          fvec = b(1,:,mode) - timestep * rhs(1,:,mode) / 2 + fvecn(1,:,mode)
        case(3) ! BDF2
          smat = - timestep * vemopt%drhs(1,:,:,mode)
          do i = 1, ncompb
            smat(i,i) = smat(i,i) + gamma0
          end do
          fvec = gamma0*b(1,:,mode) - timestep * rhs(1,:,mode) + fvecn(1,:,mode)
        case default
          write(*,'(/a,i0/)') 'Error: wrong value method: ', method
          stop
        end select

        call dgesv( ncompb, 1, smat, ncompb, ipiv, fvec, ncompb, info )

        if ( info /= 0 ) then
          print *, 'Lapack DGESV: info = ', info
          stop
        end if

!       update for Newton-Raphson: un+1=un-J(un)^-1*residual(un)

        b(1,:,mode) = b(1,:,mode) - fvec

!        print *, maxval(abs(fvec))
        write(13,*) maxval(abs(fvec))

        if ( maxval(abs(fvec)) < eps ) exit

      end do

    end do

!   compute stress tensor

    call stress_viscoelastic_3D_b ( vemodel, b, tau )

    write(14,*) time, tau
    write(15,*) time, b

  end do

  close(13)
  close(14)
  close(15)

! delete the model

  call delete ( vemodel )
  call delete ( vemopt )

end program startup40
