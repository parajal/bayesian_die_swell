! Steady flow with a constant velocity gradient.
! Newton-Raphson iteration.
! The model is an evp model with Drucker-Prager plasticity (Saramito 2021).
! Output is the viscoelastic (polymer) extra-stress.

program steady1

  use meshgen_basic_common_m
  use viscoelastic_models_3D_m

  implicit none

  integer, parameter :: nmodes  = 1,  & ! number of modes
                        model   = 22, & ! Saramito DP
                        coorsys = 2,  & ! 3D Cartesian
                        ncompc  = 6,  & ! number of components of c
                        ncompt  = 6,  & ! number of components of tau
                        ndim    = 3     ! number of coordinate directions

  logical :: includesolvent = .false. ! solvent included only if set to .true.

  integer :: flowtype = 1, &  ! flow type:
                              ! 1: shear
                              ! 2: planar extension
                              ! 3: uniaxial extension
             numdatapoints = 10, & ! number of steady data points
             maxnumiterations = 100     ! maximum number of iterations

  real(dp) :: eta_s = 0._dp,          &  ! solvent viscosity
              start_rate = 2.5e-2_dp, &  ! start data point
              end_rate = 2.5e-2_dp,   &  ! end data point
              factor = 2._dp,         &  ! factor last/first data point distance
              eps = 1.e-10_dp            ! Newton-Raphson accuracy

  real(dp), dimension(nmodes) :: &
               G,        & ! modulus
               lambda,   & ! relaxation time
               tau_y,    & ! yield stress (cohesion)
               mu          ! friction coefficient

  real(dp), allocatable, dimension(:) :: rate ! shear rate (for flowtype=1)
                                              ! strain rate (for flowtype=2,3)

  integer :: mode, iter, i, info, ipiv(ncompc)

  real(dp) :: tau(1,ncompt), gradv(1,ndim,ndim)
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes)
  real(dp) :: smat(ncompc,ncompc), taumod(ncompt), fvec(ncompc)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemopt

  G      = 10            ! modulus G
  lambda = 0.1           ! relaxation time
  tau_y  = 0.2           ! tau_y
  mu     = 0.2           ! mu
  eps = 1.e-10           ! Newton-Raphson accuracy
  maxnumiterations = 100 ! maximum number of iterations
  start_rate = 0.1     ! start data point
  end_rate = 10          ! end data point
  factor = 100         ! data point distance: last=factor * first


! namelist for input of variables; read from standard input

!  namelist /comppar/ flowtype, includesolvent, eta_s, G, lambda, tau_y, mu, &
!                     numdatapoints, start_rate, end_rate, factor, eps, &
!                     maxnumiterations

!  read ( unit=*, nml=comppar )

! create vemopt structure

  call create ( vemopt, compute_drhs=.true., ncomp=ncompc, nmodes=nmodes )

! set data points

  allocate ( rate(numdatapoints) )

  call distribute_elements ( numdatapoints - 1, rate, ratio=1, factor=factor )

  rate = start_rate + rate * (end_rate - start_rate)

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, coorsys )

! set material parameters

  vemodel%modulus = G
  vemodel%lambda = lambda
  vemodel%nonlin(1,:) = tau_y
  vemodel%nonlin(2,:) = mu

! open output files

!  open ( unit=13, file='iter.out', recl=300 )
!  open ( unit=14, file='steady1.out', recl=300 )

! initialize c, depending on flow type

  select case ( flowtype )
  case(1) ! shear
    do mode = 1, nmodes
      c(1,:,mode) = [ 1, 0, 0, 1, 0, 1 ]*(1+tau_y(mode)/mu(mode)/G(mode)/3) + &
                    [ 0, 1, 0, 0, 0, 0 ]*tau_y(mode)/G(mode)
    end do
  case(2) ! planar extension
    do mode = 1, nmodes
      c(1,:,mode) = [ 1, 0, 0, 1, 0, 1 ] + &
                    [ 1, 0, 0, 0, 0, 0 ]*tau_y(mode)/mu(mode)/G(mode)
    end do
  case(3) ! uniaxial extension
    do mode = 1, nmodes
      c(1,:,mode) = [ 1, 0, 0, 1, 0, 1 ] + &
                    [ 1, 0, 0, 0, 0, 0 ]*tau_y(mode)/mu(mode)/G(mode)
    end do
  case default
    write(*,'(/a/)') 'Error: wrong function number: '
    stop
  end select

! loop over data points

  do i = 1, numdatapoints

!   initialize velocity gradient

    select case ( flowtype )
    case(1)
      gradv(1,1,:) = [ 0._dp, rate(i), 0._dp ]
      gradv(1,2,:) = [ 0._dp,   0._dp, 0._dp ]
      gradv(1,3,:) = [ 0._dp,   0._dp, 0._dp ]
    case(2)
      gradv(1,1,:) = [ rate(i),    0._dp, 0._dp ]
      gradv(1,2,:) = [   0._dp, -rate(i), 0._dp ]
      gradv(1,3,:) = [   0._dp,    0._dp, 0._dp ]
    case(3)
      gradv(1,1,:) = [ rate(i),      0._dp,      0._dp ]
      gradv(1,2,:) = [   0._dp, -rate(i)/2,      0._dp ]
      gradv(1,3,:) = [   0._dp,      0._dp, -rate(i)/2 ]
    case default
      write(*,'(/a/)') 'Error: wrong function number: '
      stop
    end select

!   loop over all modes

    do mode = 1, nmodes

!     Newton-Raphson for each mode

      iter = 0

      do

        iter = iter + 1

        if ( iter > maxnumiterations ) then
          write(*,'(2(a,i0/),a,g0/)') &
            ' Maximum number of iterations reached = ', &
            maxnumiterations, ' datapoint = ', i, ' rate = ', rate(i)
          stop
        end if

!       right-hand side for mode

        call rhs_viscoelastic_3D ( vemodel, gradv, c, rhs, mode, mode, &
          vemopt=vemopt )

!       stress tensor, trace tau and von Mises for mode

        call stress_viscoelastic_3D_single_mode ( vemodel, c(1,:,mode), &
          taumod, mode )

!       solve J^-1 * residual

        smat = - vemopt%drhs(1,:,:,1)
        fvec = - rhs(1,:,mode)

        call dgesv( ncompc, 1, smat, ncompc, ipiv, fvec, ncompc, info )

        if ( info /= 0 ) then
          print *, 'Lapack DGESV: info = ', info
          stop
        end if

!       update for Newton-Raphson: un+1=un-J(un)^-1*residual(un)

        c(1,:,mode) = c(1,:,mode) - fvec

!        print *, maxval(abs(fvec))
        write(*,*) maxval(abs(fvec))

        if ( maxval(abs(fvec)) < eps ) exit

      end do

    end do

!   Output

    call stress_viscoelastic_3D ( vemodel, c, tau )

    call output ( rate(i), L=gradv(1,:,:), tau=tau(1,:) )

  end do

  close(13)
  close(14)

! delete the model

  call delete ( vemodel )
  call delete ( vemopt )

  deallocate ( rate )

contains


! Regime for a single mode of the evp model of Saramito (2021) based on
! Drucker-Prager plasticity.

  function regime ( tau_y, mu, trtau, tau_e )

    real(dp), intent(in) :: tau_y, mu, trtau, tau_e

    integer :: regime

    real(dp) :: chpar

    chpar = mu * trtau - 3 * tau_y

    if ( chpar <= - 3 * tau_e ) then

!     regime I: sticking

      regime = 1

    else if ( chpar >= 2 * mu**2 * tau_e ) then

!     regime III: loosing contact

      regime = 3

    else

!     regime II: sliding

      regime = 2

    end if

  end function regime


! solvent stress
! NOTE: the following variables are inherited from the host:
!  ndim, ncompt, eta_s

  function tau_s ( L )

!  velocity gradient
   real(dp), dimension(ndim,ndim), intent(in) :: L
!  stress tensor
   real(dp), dimension(ncompt) :: tau_s

   tau_s(1:3) = eta_s * ( L(1,:) + L(:,1) )
   tau_s(4:5) = eta_s * ( L(2,2:3) + L(2:3,2) )
   tau_s(6)   = 2 * eta_s * L(2,2)

  end function tau_s


! output
! NOTE: the following variables are inherited from the host:
!  nmodes, includesolvent, tau_y, mu

  subroutine output ( rate, L, tau )

    real(dp), intent(in) :: rate
!   velocity gradient
    real(dp), dimension(ndim,ndim), intent(in) :: L
!   stress tensor
    real(dp), dimension(ncompt), intent(in) :: tau

    real(dp) :: trtau, tau_e

    if ( nmodes == 1 ) then

      trtau = tau(1) + tau(4) + tau(6)
      tau_e = vonmises_3D(tau)

      if ( includesolvent) then
        write(*,'(9es16.8,1x,i0)') rate, tau + tau_s(L), tau_e, trtau/3, &
          regime(tau_y(1), mu(1), trtau, tau_e )
      else
        write(*,'(9es16.8,1x,i0)') rate, tau, tau_e, trtau/3, &
          regime(tau_y(1), mu(1), trtau, tau_e )
      end if

    else

      if ( includesolvent) then
        write(*,'(7es16.8)') rate, tau + tau_s(L)
      else
        write(*,'(7es16.8)') rate, tau
      end if

    end if

  end subroutine output

end program steady1

