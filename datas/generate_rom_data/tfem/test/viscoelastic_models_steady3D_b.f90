! Steady flow with a constant velocity gradient.
! Newton-Raphson iteration.
! The model is an evp model (Saramito 2007, SRM1) with Oldroyd-B or Giesekus
! relaxation expression.
! Output is the viscoelastic (polymer) extra-stress.
! Contravariant deformation formulation. 3D flow.

program steady9

  use meshgen_basic_common_m
  use viscoelastic_models_3D_b_m
  use limits_m

  implicit none

  integer, parameter :: nmodes  = 1,  & ! number of modes
                        coorsys = 2,  & ! 3D Cartesian
                        ncompb  = 9,  & ! number of components of b
                        ncompt  = 6,  & ! number of components of tau
                        ndim    = 3     ! number of coordinate directions

  logical :: includesolvent = .false.    ! solvent included only if set .true.

  integer :: model   = 2,  & ! 2: Oldroyd model 3: Giesekus
             alam_model = 2, &  ! adapted lambda model SRM1
             flowtype = 1, &  ! flow type:
                              ! 1: shear
                              ! 2: planar extension
                              ! 3: uniaxial extension
             method = 1,   &  ! method 1: BDF1, 2: BDF2
             numdatapoints = 10, & ! number of steady data points
             maxnumiterations = 100     ! maximum number of iterations

  real(dp) :: eta_s = 0._dp,          &  ! solvent viscosity
              start_rate = 2.5e-2_dp, &  ! start data point
              end_rate = 2.5e-1_dp,   &  ! end data point
              timestep = 1.e-3_dp,    &  ! time step for steady
              factor = 2._dp,         &  ! factor last/first data point distance
              eps = 1.e-10_dp            ! Newton-Raphson accuracy

  real(dp), dimension(nmodes) :: &
               G,        & ! modulus
               lambda,   & ! relaxation time
               mobility = 0.3_dp, & ! mobility parameter
               tau_y       ! yield stress

  real(dp), allocatable, dimension(:) :: rate ! shear rate (for flowtype=1)
                                              ! strain rate (for flowtype=2,3)

  integer :: mode, iter, i, info, ipiv(ncompb)

  real(dp) :: tau(1,ncompt), gradv(1,ndim,ndim)
  real(dp) :: b(1,ncompb,nmodes), rhs(1,ncompb,nmodes)
  real(dp) :: Hb(1,ncompb,nmodes), dHb(1,ncompb,ncompb,nmodes)
  real(dp) :: smat(ncompb,ncompb), taumod(ncompt), fvec(ncompb)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemopt

  type(tdpar_t) :: tdpar

  G      = 10            ! modulus G
  lambda = 0.1           ! relaxation time
  tau_y  = 0.2           ! tau_y
  eps = 1.e-10           ! Newton-Raphson accuracy
  maxnumiterations = 100 ! maximum number of iterations
  start_rate = 0.1     ! start data point
  end_rate = 10          ! end data point
  factor = 100         ! data point distance: last=factor * first


! namelist for input of variables; read from standard input

!  namelist /comppar/ flowtype, includesolvent, eta_s, model, G, lambda, &
!                     mobility, tau_y, numdatapoints, start_rate, end_rate, &
!                     factor, eps, maxnumiterations, method, timestep, &
!                     alam_model

!  read ( unit=*, nml=comppar )

  USE_CONFORMATION_MODEL = .false.


! create vemopt structure

  call create ( vemopt, compute_drhs=.true., ncomp=ncompb, nmodes=nmodes )

! set data points

  allocate ( rate(numdatapoints) )

  call distribute_elements ( numdatapoints - 1, rate, ratio=1, factor=factor )

  rate = start_rate + rate * (end_rate - start_rate)

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, coorsys, &
    alam_model=alam_model  )

! set material parameters

  vemodel%modulus = G
  vemodel%lambda = lambda
  if ( model == 3 ) vemodel%nonlin(1,:) = mobility
  if ( alam_model == 2 ) vemodel%alam(1,:) = tau_y

! set numerical parameters for time discretization

  tdpar%timestep = timestep
  tdpar%method = method

! open output files

!  open ( unit=13, file='iter.out', recl=300 )
!  open ( unit=14, file='steady9.out', recl=300 )

! initialize c, depending on flow type

  if ( alam_model == 2 ) then

    select case ( flowtype )
    case(1) ! shear
      do mode = 1, nmodes
        b(1,:,mode) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ] + &
                      [ 0, 1, 0, 0, 0, 0, 0, 0, 0 ]*tau_y(mode)/G(mode)
        call sqrtc_3D_b( b(:,:,mode) )
      end do
    case(2) ! planar extension
      do mode = 1, nmodes
        b(1,:,mode) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ] + &
                      [ 1, 0, 0, 0, 0, 0, 0, 0, 0 ]*tau_y(mode)/G(mode)
        call sqrtc_3D_b( b(:,:,mode) )
      end do
    case(3) ! uniaxial extension
      do mode = 1, nmodes
        b(1,:,mode) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ] + &
                      [ 1, 0, 0, 0, 0, 0, 0, 0, 0 ]*tau_y(mode)/G(mode)
        call sqrtc_3D_b( b(:,:,mode) )
      end do
    case default
      write(*,'(/a/)') 'Error: wrong function number: '
      stop
    end select

  else

    do mode = 1, nmodes
      b(1,:,mode) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ]
    end do

  end if


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

        call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs, mode, mode, &
          vemopt=vemopt )

        call NRtd_steady_viscoelastic_3D_b ( vemodel, tdpar, b, Hb, dHb, &
          mode, mode )

!       stress tensor, trace tau and von Mises for mode

        call stress_viscoelastic_3D_single_mode_b ( vemodel, b(1,:,mode), &
          taumod, mode )

!       solve J^-1 * residual

        smat = - vemopt%drhs(1,:,:,mode) + dHb(1,:,:,mode)
        fvec = - rhs(1,:,mode) + Hb(1,:,mode)

        call dgesv( ncompb, 1, smat, ncompb, ipiv, fvec, ncompb, info )

        if ( info /= 0 ) then
          print *, 'Lapack DGESV: info = ', info
          stop
        end if

!       update for Newton-Raphson: un+1=un-J(un)^-1*residual(un)

        b(1,:,mode) = b(1,:,mode) - fvec

!        print *, maxval(abs(fvec))
        write(*,*) iter, maxval(abs(fvec))

        if ( maxval(abs(fvec)) < eps ) exit

      end do

!      print *, detb_3D(b(1,:,mode)) - 1._dp

    end do

!   Output

    call stress_viscoelastic_3D_b ( vemodel, b, tau )

    call output ( rate(i), L=gradv(1,:,:), tau=tau(1,:) )

  end do

!  close(13)
!  close(14)

! delete the model

  call delete ( vemodel )
  call delete ( vemopt )

  deallocate ( rate )

contains


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
!  nmodes, includesolvent, tau_y

  subroutine output ( rate, L, tau )

    real(dp), intent(in) :: rate
!   velocity gradient
    real(dp), dimension(ndim,ndim), intent(in) :: L
!   stress tensor
    real(dp), dimension(ncompt), intent(in) :: tau

    if ( nmodes == 1 ) then

      if ( includesolvent) then
        write(*,'(9es16.8,1x,i0)') rate, tau + tau_s(L), vonmises_3D(tau)
      else
        write(*,'(9es16.8,1x,i0)') rate, tau, vonmises_3D(tau)
      end if

    else

      if ( includesolvent) then
        write(*,'(7es16.8)') rate, tau + tau_s(L)
      else
        write(*,'(7es16.8)') rate, tau
      end if

    end if

  end subroutine output

end program steady9

