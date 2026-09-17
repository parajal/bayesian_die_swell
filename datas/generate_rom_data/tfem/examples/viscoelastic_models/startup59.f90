! Startup from zero stress of flow with a constant velocity gradient
! 3D flow. Combined shear and uniaxial extensional flow.
! The model is an EGP elastoviscoplastic model (single mode)
! Plastic strain softening.
! Implicit trapezoidal time integration
! Similar to startup43, but now with plastic strain softening. Single mode only.
! Similar to startup57, but now using implicit trapezoidal time integration.

program startup59

  use viscoelastic_models_3D_m
  use limits_m, only: USE_EGP_STRESS_TENSOR_FORM

  implicit none

  integer, parameter :: model    = 23, & ! 23: EGP model incompressible
                                         ! 24: EGP model compressible
                        flowtype = 2, &  ! 3D
                        alam_model = 4, & ! adapted lambda model:
                                          ! 4: Eyring, 5: Ree-Eyring
                        alam_gammap_model = 1, & ! adapted lambda gammap model:
                                          ! 1: EGP strain softening
                        ncompc   = 6, &  ! number of components of c
                        ncompt   = 6, &  ! number of components of tau
                        ndim     = 3, &  ! number of coordinate directions
                maxnumiterations = 100, &! maximum number of iterations
                  numtimesteps   = 100 ! number of time steps

  integer, parameter :: n = ncompc + 1 ! number of equations to solve

  real(dp), parameter :: modulus = 4000, & ! modulus
                         lambda  = 10, & ! relaxation time
                         tau_ref = 500, & ! reference von Mises
                         tau_ref1 = 500, & ! reference von Mises 1 RE
                         tau_ref2 = 2000, & ! reference von Mises 2 RE
                         f1 = 0.6_dp, & ! weight factor for first term RE
                         Sa = 7.4_dp,  & ! state parameter
                         r0 = 0.96_dp, & ! fitting parameter r0
                         r1 = 20.0_dp, & ! fitting parameter r1
                         r2 = -2.0_dp, & ! fitting parameter r2
                         epsilondot = 0.0_dp, &  ! elongational rate
                         gammadot = 0.25_dp, &   ! shear rate
                         Lxx_add = 0.01_dp, & ! added value for Lxx (compr.)
                         Jvol = 0.5_dp, &  ! change in volume for model=24
                         timestep = 0.1_dp,& ! time step
                         eps = 1.e-10_dp   ! Newton-Raphson accuracy

  logical, parameter :: dep_J = model==24

  integer :: step, iter, i, info, ipiv(n)
  real(dp) :: c(1,ncompc,1), rhs(1,ncompc,1), tau(1,ncompt)
  real(dp) :: gradv(1,ndim,ndim)
  real(dp) :: smat(n,n), fvec(n), fvecn(1,ncompc,1), fvec_gpn(1)
  real(dp) :: gammap(1), rhs_gammap(1)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemoptn, vemopt
  type(vemopt_gammap_t) :: vemopt_gammapn, vemopt_gammap

  USE_EGP_STRESS_TENSOR_FORM = .false.

! define the model

  call create_viscoelastic_model ( model, vemodel, flowtype=flowtype, &
    alam_model=alam_model, alam_gammap_model=alam_gammap_model, &
    gammap_mode=1 )

! set optional input parameters

  call create_vemopt ( vemoptn, dep_J=dep_J, dep_gammap=.true. )
  call create_vemopt ( vemopt, dep_J=dep_J, compute_drhs=.true., &
    ncomp=ncompc, nmodes=1, dep_gammap=.true., &
    compute_drhsdgammap=.true. )
  call create_vemopt_gammap ( vemopt_gammapn, dep_J=dep_J )
  call create_vemopt_gammap ( vemopt_gammap, dep_J=dep_J, &
    compute_drhs=.true., ncomp=ncompc, compute_drhsdgammap=.true. )

  if ( dep_J ) then
    vemoptn%J = Jvol
    vemopt_gammapn%J = Jvol
    vemopt%J = Jvol
    vemopt_gammap%J = Jvol
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
  if ( alam_gammap_model == 1 ) then
    vemodel%alam_gammap(1,1) = Sa
    vemodel%alam_gammap(2,1) = r0
    vemodel%alam_gammap(3,1) = r1
    vemodel%alam_gammap(4,1) = r2
  end if

! initialize velocity gradient

  gradv(1,1,:) = [ epsilondot,      gammadot,         0._dp ]
  gradv(1,2,:) = [      0._dp, -epsilondot/2,         0._dp ]
  gradv(1,3,:) = [      0._dp,         0._dp, -epsilondot/2 ]

! EGP compressible trL/=0
  if ( model == 24 ) gradv(1,1,1) = gradv(1,1,1) + Lxx_add

! initialize c

  c(1,:,1) = [ 1, 0, 0, 1, 0, 1 ]

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

    call rhs_viscoelastic_3D ( vemodel, gradv, c, rhs, vemopt=vemoptn )
    call rhs_plastic_strain_3D ( vemodel, c, gammap, rhs_gammap, &
      vemopt_gammap=vemopt_gammapn )

    fvecn = -c - timestep * rhs / 2
    fvec_gpn = -gammap - timestep * rhs_gammap / 2

!   Newton-Raphson for single mode

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

!     right-hand side

      call rhs_viscoelastic_3D ( vemodel, gradv, c, rhs, vemopt=vemopt )
      call rhs_plastic_strain_3D ( vemodel, c, gammap, rhs_gammap, &
        vemopt_gammap=vemopt_gammap )

!     solve J^-1 * residual

      smat(:ncompc,:ncompc) = - timestep * vemopt%drhs(1,:,:,1) / 2
      smat(:ncompc,n) = - timestep * vemopt%drhsdgammap(1,:,1) / 2
      smat(n,:ncompc) = - timestep * vemopt_gammap%drhs(1,:) / 2
      smat(n,n) = - timestep * vemopt_gammap%drhsdgammap(1) / 2

      do i = 1, n
        smat(i,i) = smat(i,i) + 1
      end do

      fvec(:ncompc) = reshape( c - timestep * rhs / 2 + fvecn, [ncompc] )
      fvec(n:) = gammap - timestep * rhs_gammap / 2 + fvec_gpn

      call dgesv( ncompc, 1, smat, n, ipiv, fvec, n, info )

      if ( info /= 0 ) then
        print *, 'Lapack DGESV: info = ', info
        stop
      end if

!     update for Newton-Raphson: un+1=un-J(un)^-1*residual(un)

      c = c - reshape ( fvec(:ncompc), [1,ncompc,1] )
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

end program startup59
