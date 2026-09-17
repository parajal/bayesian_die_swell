! Startup from zero stress of flow with a constant velocity gradient
! 3D flow. Combined shear and uniaxial extensional flow.
! The model is a Maxwell, Giesekus of FENE-P model, one-modes
! The integration scheme is the trapezoidal method.
! Using the contravariant deformation tensor (b).
! Output of eigenvalues and eigenvectors of Jacobian of rhs.
! Similar to startup37, but now just one mode and computation of eigen problem.

program startup41

  use viscoelastic_models_3D_b_m
!  use limits_m, only: USE_CONFORMATION_MODEL

  implicit none

  logical, parameter :: reinitialize = .true.

  integer, parameter :: model    = 3, & ! 2: Maxwell 3: Giesekus 20: FENE-P
                        nmodes   = 1,  & ! number of modes
                        flowtype = 2,  & ! 3D
                        ncompb   = 9,  & ! number of components of b
                        ncompt   = 6,  & ! number of components of tau
                        ndim     = 3,  & ! number of coordinate directions
                maxnumiterations = 100, &! maximum number of iterations
                  numtimesteps   = 1000 ! number of time steps

  real(dp), parameter :: modulus   = 2, & ! modulus
                         lambda    = 10._dp, & ! relaxation time
                         mobility  = 0.1_dp, & ! alpha
                         bpar      = 100.0_dp, & ! bparameter
                         epsilondot = 0.0_dp, &  ! elongational rate
                         gammadot = 1.0_dp, &   ! shear rate
                         timestep = 1.e-1_dp,& ! time step
                              eps = 1.e-10_dp   ! Newton-Raphson accuracy


  integer :: mode, step, iter, i, info, ipiv(ncompb)
  real(dp) :: b(1,ncompb,nmodes), rhs(1,ncompb,nmodes), tau(1,ncompt)
  real(dp) :: gradv(1,ndim,ndim)
  real(dp) :: smat(ncompb,ncompb), fvec(ncompb), fvecn(1,ncompb,nmodes)

  integer, parameter :: lwork = 780
  real(dp) :: wr(ncompb), wi(ncompb), vl(ncompb,ncompb), vr(ncompb,ncompb)
  real(dp) :: work(lwork)

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

! open output file

  open ( unit=13, file='iter.out', recl=300 )
  open ( unit=14, file='out', recl=300 )
  open ( unit=15, file='eigen41.out', recl=300 )

! Solve db/dt = rhs(b) using trapezoidal method:
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

  do step = 1, numtimesteps

!   constant part of right-hand side

    call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs )

    fvecn = -b - timestep * rhs / 2

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

!        print *, maxval(abs(fvec))
        write(13,*) maxval(abs(fvec))

        if ( maxval(abs(fvec)) < eps ) exit

      end do

    end do

    if ( reinitialize ) call sqrtc_3D_b ( b(:,:,1) )

!   compute stress tensor

    call stress_viscoelastic_3D_b ( vemodel, b, tau )

    write(14,*) step * timestep, tau
!    write(14,*) step * timestep, b(:,:,1)

    call output_eigen

  end do

  close(13)
  close(14)
  close(15)

! delete the model

  call delete ( vemodel )
  call delete ( vemopt )

contains

! output eigenvalues and eigenvector of the Jacobian of the rhs

  subroutine output_eigen

    integer :: j

!   right-hand side for mode=1

    call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs, mode1=1, mode2=1, &
      vemopt=vemopt )

!   eigenvalues

    call dgeev( 'N', 'V', ncompb, vemopt%drhs(1,:,:,1), ncompb, wr, wi, vl, &
      ncompb, vr, ncompb, work, lwork, info )

    write(15,'(a,g0)') 'gammadot = ', gammadot
    write(15,'(a,g0)') 'epsilondot = ', epsilondot
    write(15,'(a,g0)') 'time = ', step * timestep
    write(15,'(a,i0)') 'info = ', info
    if ( info == 0 ) &
              write(15,'(a,i0)') 'optimal lwork: work(1) = ', nint(work(1))
    write(15,'(a)') ' Eigenvalues: j, wr(j), wi(j) '
    do j = 1, ncompb
      write(15,'(i3,2es16.8)') j, wr(j), wi(j)
    end do
    write(15,'(a)') ' Eigenvectors: j, v(j) '
    do j = 1, ncompb
      write(15,'(i3,9es16.8)') j, vr(:,j)
    end do

  end subroutine output_eigen

end program startup41
