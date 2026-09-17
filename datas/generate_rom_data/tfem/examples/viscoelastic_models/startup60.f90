! Startup from zero stress of flow with a constant velocity gradient
! 3D flow. Combined shear and uniaxial extensional flow.
! The model is an EGP elastoviscoplastic model
! Plastic strain softening.
! Implicit trapezoidal time integration
! Using the standard conformation tensor (c).
! 14-mode model with 12 alpha process viscoelastic modes, one elastic mode for
! strain hardening and one mode for the beta process.
! Data for PMMA from van Breemen et al. (2012) DOI: 10.1002/polb.23199
! Similar to startup56, but now with plastic strain softening.

program startup60

  use viscoelastic_models_3D_m
  use limits_m, only: USE_EGP_STRESS_TENSOR_FORM
  use material_data_m

  implicit none

  integer, parameter :: model    = 23, & ! 23: EGP model incompressible
                                         ! 24: EGP model compressible
                        nmodes   = 14, &  ! number of modes
                        am = nmodes-2, & ! 1:am range of modes for alpha process
                        em = nmodes-1, &  ! elastic mode
                        bm = nmodes,   &  ! mode of beta process
                        flowtype = 2, &  ! 3D
                        alam_model = 4, & ! adapted lambda model:
                                          ! 4: Eyring
                        alam_gammap_model = 1, & ! adapted lambda gammap model:
                                          ! 1: EGP strain softening
                        ncompc   = 6, &  ! number of components of c
                        ncompt   = 6, &  ! number of components of tau
                        ndim     = 3, &  ! number of coordinate directions
                maxnumiterations = 100, &! maximum number of iterations
                  numtimesteps   = 100 ! number of time steps

  integer, parameter :: nc = ncompc*nmodes, & ! number of conformation unknows
                        n = nc + 1    ! number of equations to solve

  real(dp), parameter :: epsilondot = 0.e-2_dp, &  ! elongational rate
                         gammadot = 1.e-2_dp, &   ! shear rate
                         Jvol = 1.0_dp, &  ! change in volume for model=24
                         timestep = 1.0_dp, & ! time step
                         eps = 1.e-10_dp   ! Newton-Raphson accuracy

  logical, parameter :: dep_J = model==24

  integer :: m, step, iter, i, info, ipiv(n)
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), tau(1,ncompt)
  real(dp) :: gradv(1,ndim,ndim)
  real(dp) :: smat(n,n), fvec(n), fvecn(1,ncompc,nmodes), fvec_gpn(1)
  real(dp) :: gammap(1), rhs_gammap(1)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemoptn, vemopt
  type(vemopt_gammap_t) :: vemopt_gammapn, vemopt_gammap

  USE_EGP_STRESS_TENSOR_FORM = .false.

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype, &
    alam_model=alam_model, vmglobal=.true., norlxmode=em, &
    alam_gammap_model=alam_gammap_model, gammap_mode=1 )

! fill material parameters

  call fill_data ( vemodel, 'PMMA' )

! set optional input parameters

  call create_vemopt ( vemoptn, dep_J=dep_J, dep_gammap=.true. )
  call create_vemopt ( vemopt, dep_J=dep_J, compute_drhs_mm=.true., &
    ncomp=ncompc, nmodes=nmodes, dep_gammap=.true., &
    compute_drhsdgammap=.true. )
  call create_vemopt_gammap ( vemopt_gammapn, dep_J=dep_J )
  call create_vemopt_gammap ( vemopt_gammap, dep_J=dep_J, &
    compute_drhs_mm=.true., ncomp=ncompc, nmodes=nmodes, &
    compute_drhsdgammap=.true. )

  if ( dep_J ) then
    vemoptn%J = Jvol
    vemopt_gammapn%J = Jvol
    vemopt%J = Jvol
    vemopt_gammap%J = Jvol
  end if

! initialize velocity gradient

  gradv(1,1,:) = [ epsilondot,      gammadot,         0._dp ]
  gradv(1,2,:) = [      0._dp, -epsilondot/2,         0._dp ]
  gradv(1,3,:) = [      0._dp,         0._dp, -epsilondot/2 ]

! initialize c

  do m = 1, nmodes
    c(1,:,m) = [ 1, 0, 0, 1, 0, 1 ]
  end do

! initialize gammap

  gammap(1) = 0
  vemoptn%gammap = gammap
  vemopt%gammap = gammap

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

    call rhs_3D ( vemoptn, vemopt_gammapn )

    fvecn = -c - timestep * rhs / 2
    fvec_gpn = -gammap - timestep * rhs_gammap / 2

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
      vemopt_gammap%drhs_mm = 0 ! initialize all to zero

      call rhs_3D ( vemopt, vemopt_gammap )

!     solve J^-1 * residual

      smat(:nc,:nc) = - timestep * reshape ( vemopt%drhs_mm, [nc,nc] ) / 2
      smat(:nc,n) = - timestep * reshape ( vemopt%drhsdgammap, [nc] ) / 2
      smat(n,:nc) = - timestep * reshape ( vemopt_gammap%drhs_mm, [nc] ) / 2
      smat(n,n) = - timestep * vemopt_gammap%drhsdgammap(1) / 2

      do i = 1, n
        smat(i,i) = smat(i,i) + 1
      end do

      fvec(:nc) = reshape ( c - timestep * rhs / 2 + fvecn, [nc] )
      fvec(n:) = gammap - timestep * rhs_gammap / 2 + fvec_gpn

      call dgesv( n, 1, smat, n, ipiv, fvec, n, info )

      if ( info /= 0 ) then
        print *, 'Lapack DGESV: info = ', info
        stop
      end if

!     update for Newton-Raphson: un+1=un-J(un)^-1*residual(un)

      c = c - reshape ( fvec(:nc), [1,ncompc,nmodes] )
      gammap = gammap - fvec(n)
      vemopt%gammap = gammap

      !print *, maxval(abs(fvec)), abs(fvec(n))
      write(13,*) maxval(abs(fvec)), abs(fvec(n))

      if ( maxval(abs(fvec)) < eps ) exit

    end do

!   update vemoptn for next time step

    vemoptn%gammap = gammap

!   compute stress tensor

    call stress_viscoelastic_3D ( vemodel, c, tau, vemopt=vemopt )

    write(14,*) step * timestep, tau, vonmises_3D(tau(1,:)), &
      detc_3D(c(1,:,1)), gammap(1)

  end do

  close(13)
  close(14)

! delete the model

  call delete ( vemodel )
  call delete ( vemoptn, vemopt )
  call delete ( vemopt_gammapn, vemopt_gammap )

contains

! build rhs for given vemopt

  subroutine rhs_3D ( vemopt, vemopt_gammap )

    type(vemopt_t), intent(inout) :: vemopt
    type(vemopt_gammap_t), intent(inout) :: vemopt_gammap

!   alpha process
    call rhs_viscoelastic_3D ( vemodel, gradv, c, rhs, mode1=1, &
        mode2=am, vemopt=vemopt )

!   elastic mode
    call rhs_viscoelastic_3D ( vemodel, gradv, c, rhs, mode1=em, &
      mode2=em, vemopt=vemopt )

!   beta process
    call rhs_viscoelastic_3D ( vemodel, gradv, c, rhs, mode1=bm, &
        mode2=bm, vemopt=vemopt )

!   plastic strain
    call rhs_plastic_strain_3D ( vemodel, c, gammap, rhs_gammap, &
      mode1=1, mode2=am, vemopt_gammap=vemopt_gammap )

  end subroutine rhs_3D

end program startup60
