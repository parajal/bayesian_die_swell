! Startup from zero stress of flow with a constant velocity gradient
! 2D flow. Combined shear and uniaxial extensional flow.
! The model is an EGP elastoviscoplastic model
! Plastic strain softening.
! Implicit trapezoidal time integration
! Using the contravariant deformation tensor (b).
! 14-mode model with 12 alpha process viscoelastic modes, one elastic mode for
! strain hardening and one mode for the beta process.
! Data for PMMA from van Breemen et al. (2012) DOI: 10.1002/polb.23199
! Similar to startup61, but now with 2D flow.

program startup63

  use viscoelastic_models_2D_b_m
  use limits_m, only: USE_CONFORMATION_MODEL, USE_EGP_STRESS_TENSOR_FORM
  use material_data_m

  implicit none

  logical, parameter :: reinitialize = .false.

  integer, parameter :: model = 23,  & ! 23: EGP model incompressible
                                       ! 24: EGP model compressible
                        nmodes = 14, & ! number of modes
                        am = 12,     & ! 1:am range of modes for alpha process
                        em = 13,     & ! elastic mode
                        bm = 14,     & ! mode of beta process
                        flowtype = 0, &  ! 2D
                        alam_model = 4, & ! adapted lambda model:
                                          ! 4: Eyring
                        alam_gammap_model = 1, & ! adapted lambda gammap model:
                                          ! 1: EGP strain softening
                        ncompb   = 5, &  ! number of components of b
                        ncompt   = 4, &  ! number of components of tau
                        ncompg   = 4, &  ! number of components of gradv
                maxnumiterations = 100, &! maximum number of iterations
                  numtimesteps   = 100 ! number of time steps

  integer, parameter :: nb = ncompb*nmodes, & ! number of conformation unknows
                        n = nb + 1    ! number of equations to solve

  real(dp), parameter :: epsilondot = 0.e-2_dp, &  ! elongational rate
                         gammadot = 1.e-2_dp, &   ! shear rate
                         Jvol = 1.0_dp, &  ! change in volume for model=24
                         timestep = 1.0_dp, & ! time step
                         eps = 1.e-10_dp   ! Newton-Raphson accuracy

  logical, parameter :: dep_J = model==24

  integer :: m, step, iter, i, info, ipiv(n)
  real(dp) :: b(1,ncompb,nmodes), rhs(1,ncompb,nmodes), tau(1,ncompt)
  real(dp) :: gradv(1,ncompg), L(ncompg)
  real(dp) :: smat(n,n), fvec(n), fvecn(1,ncompb,nmodes), fvec_gpn(1)
  real(dp) :: gammap(1), rhs_gammap(1)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemoptn, vemopt
  type(vemopt_gammap_t) :: vemopt_gammapn, vemopt_gammap

  USE_CONFORMATION_MODEL = .false.
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
    ncomp=ncompb, nmodes=nmodes, dep_gammap=.true., &
    compute_drhsdgammap=.true. )
  call create_vemopt_gammap ( vemopt_gammapn, dep_J=dep_J )
  call create_vemopt_gammap ( vemopt_gammap, dep_J=dep_J, &
    compute_drhs_mm=.true., ncomp=ncompb, nmodes=nmodes, &
    compute_drhsdgammap=.true. )

  if ( dep_J ) then
    vemoptn%J = Jvol
    vemopt_gammapn%J = Jvol
    vemopt%J = Jvol
    vemopt_gammap%J = Jvol
  end if

! initialize velocity gradient

  L = [ epsilondot, gammadot, 0._dp, -epsilondot ]

  gradv(1,:) = L

! initialize b

  do m = 1, nmodes
    b(1,:,m) = [ 1, 0, 0, 1, 1 ]
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

    call rhs_2D_b ( vemoptn, vemopt_gammapn )

    fvecn = -b - timestep * rhs / 2
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

      call rhs_2D_b ( vemopt, vemopt_gammap )

!     solve J^-1 * residual

      smat(:nb,:nb) = - timestep * reshape ( vemopt%drhs_mm, [nb,nb] ) / 2
      smat(:nb,n) = - timestep * reshape ( vemopt%drhsdgammap, [nb] ) / 2
      smat(n,:nb) = - timestep * reshape ( vemopt_gammap%drhs_mm, [nb] ) / 2
      smat(n,n) = - timestep * vemopt_gammap%drhsdgammap(1) / 2

      do i = 1, n
        smat(i,i) = smat(i,i) + 1
      end do

      fvec(:nb) = reshape ( b - timestep * rhs / 2 + fvecn, [nb] )
      fvec(n:) = gammap - timestep * rhs_gammap / 2 + fvec_gpn

      call dgesv( n, 1, smat, n, ipiv, fvec, n, info )

      if ( info /= 0 ) then
        print *, 'Lapack DGESV: info = ', info
        stop
      end if

!     update for Newton-Raphson: un+1=un-J(un)^-1*residual(un)

      b = b - reshape ( fvec(:nb), [1,ncompb,nmodes] )
      gammap = gammap - fvec(n)
      vemopt%gammap = gammap

      !print *, maxval(abs(fvec)), abs(fvec(n))
      write(13,*) maxval(abs(fvec)), abs(fvec(n))

      if ( maxval(abs(fvec)) < eps ) exit

    end do

!   update vemoptn for next time step

    vemoptn%gammap = gammap

!   compute stress tensor

    call stress_viscoelastic_2D_b ( vemodel, b, tau, vemopt=vemopt )

    write(14,*) step * timestep, tau, vonmises_2D(tau(1,:)), &
      detb_2D(b(1,:,1)), gammap(1)

    if ( reinitialize ) then
      do m = 1, nmodes
        call sqrtc_2D_b ( b(:,:,m) )
      end do
    end if

  end do

  close(13)
  close(14)

! delete the model

  call delete ( vemodel )
  call delete ( vemoptn, vemopt )
  call delete ( vemopt_gammapn, vemopt_gammap )

contains

! build rhs for given vemopt

  subroutine rhs_2D_b ( vemopt, vemopt_gammap )

    type(vemopt_t), intent(inout) :: vemopt
    type(vemopt_gammap_t), intent(inout) :: vemopt_gammap

!   alpha process
    call rhs_viscoelastic_2D_b ( vemodel, gradv, b, rhs, mode1=1, &
        mode2=am, vemopt=vemopt )

!   elastic mode
    call rhs_viscoelastic_2D_b ( vemodel, gradv, b, rhs, mode1=em, &
      mode2=em, vemopt=vemopt )

!   beta process
    call rhs_viscoelastic_2D_b ( vemodel, gradv, b, rhs, mode1=bm, &
        mode2=bm, vemopt=vemopt )

!   plastic strain
    call rhs_plastic_strain_2D_b ( vemodel, b, gammap, rhs_gammap, &
      mode1=1, mode2=am, vemopt_gammap=vemopt_gammap )

  end subroutine rhs_2D_b

end program startup63
