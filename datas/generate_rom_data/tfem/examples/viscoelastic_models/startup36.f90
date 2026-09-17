! Startup from zero stress of flow with a constant velocity gradient
! 3D flow. Combined shear and uniaxial extensional flow.
! The model is a Giesekus model, two-modes
! The integration scheme is the trapezoidal method.
! Similar to startup4, but now using implicit trapezoidal time integration.

program startup36

  use viscoelastic_models_3D_m

  implicit none

  integer, parameter :: model    = 3, &  ! Giesekus model
                        nmodes   = 2, &  ! number of modes
                        flowtype = 2, &  ! 3D
                        ncompc   = 6, &  ! number of components of c
                        ncompt   = 6, &  ! number of components of tau
                        ndim     = 3, &  ! number of coordinate directions
                maxnumiterations = 100, &! maximum number of iterations
                  numtimesteps   = 250 ! number of time steps

  real(dp), parameter :: modulus(nmodes) = [2,1], & ! modulus
                         lambda(nmodes)  = [3._dp,0.3_dp], & ! relaxation time
                         mobility(nmodes) = 0.1_dp, &  ! mobility parameter
                         epsilondot = 2.5_dp, &  ! elongational rate
                         gammadot = 1.5_dp, &   ! shear rate
                         timestep = 1.e-1_dp,& ! time step
                              eps = 1.e-10_dp   ! Newton-Raphson accuracy


  integer :: mode, step, iter, i, info, ipiv(ncompc)
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), tau(1,ncompt)
  real(dp) :: gradv(1,ndim,ndim)
  real(dp) :: smat(ncompc,ncompc), fvec(ncompc), fvecn(1,ncompc,nmodes)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemopt

! create vemopt structure

  call create ( vemopt, compute_drhs=.true., ncomp=ncompc, nmodes=nmodes )

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype )

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
  vemodel%nonlin(1,:) = mobility

! initialize velocity gradient

  gradv(1,1,:) = [ epsilondot,      gammadot,         0._dp ]
  gradv(1,2,:) = [      0._dp, -epsilondot/2,         0._dp ]
  gradv(1,3,:) = [      0._dp,         0._dp, -epsilondot/2 ]

! initialize c

  do mode = 1, nmodes
    c(1,:,mode) = [ 1, 0, 0, 1, 0, 1 ]
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

    call rhs_viscoelastic_3D ( vemodel, gradv, c, rhs )

    fvecn = -c - timestep * rhs / 2

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

        call rhs_viscoelastic_3D ( vemodel, gradv, c, rhs, mode, mode, &
          vemopt=vemopt )

!       solve J^-1 * residual

        smat = - timestep * vemopt%drhs(1,:,:,mode) / 2
        do i = 1, ncompc
          smat(i,i) = smat(i,i) + 1
        end do
        fvec = c(1,:,mode) - timestep * rhs(1,:,mode) / 2 + fvecn(1,:,mode)

        call dgesv( ncompc, 1, smat, ncompc, ipiv, fvec, ncompc, info )

        if ( info /= 0 ) then
          print *, 'Lapack DGESV: info = ', info
          stop
        end if

!       update for Newton-Raphson: un+1=un-J(un)^-1*residual(un)

        c(1,:,mode) = c(1,:,mode) - fvec

!        print *, maxval(abs(fvec))
        write(13,*) maxval(abs(fvec))

        if ( maxval(abs(fvec)) < eps ) exit

      end do

    end do

!   compute stress tensor

    call stress_viscoelastic_3D ( vemodel, c, tau )

    write(14,*) step * timestep, tau

  end do

  close(13)
  close(14)

! delete the model

  call delete ( vemodel )
  call delete ( vemopt )

end program startup36
