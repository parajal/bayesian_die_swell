! Steady flow with a constant velocity gradient. 2D flow.
! Newton-Raphson iteration.
! The model is a EGP elastoviscoplastic model for incompressible flow
! Output is the viscoelastic (polymer) extra-stress.
! Optional output of eigenvalues and eigenvectors of Jacobian of rhs.
! Contravariant deformation formulation.

program steady13

  use meshgen_basic_common_m
  use viscoelastic_models_2D_b_m
  use limits_m, only: USE_CONFORMATION_MODEL, USE_EGP_STRESS_TENSOR_FORM

  implicit none

  integer, parameter :: nmodes  = 1,  & ! number of modes
                        coorsys = 0,  & ! 3D Cartesian
                        ncompb  = 5,  & ! number of components of b
                        ncompt  = 4,  & ! number of components of tau
                        ncompg  = 4     ! number of components of gradv

  logical :: includesolvent = .false.    ! solvent included only if set .true.

  integer :: model = 24,      & ! 23: EGP model incompressible
                                ! 24: EGP model compressible
             alam_model = 4,  & ! adapted lambda model:
                                ! 4: Eyring, 5: Ree-Eyring
             flowtype = 1, &  ! flow type:
                              ! 1: shear
                              ! 2: planar extension
                              ! 3: combined
             method = 1,   &  ! method 1: BDF1, 2: BDF2
             numdatapoints = 10, & ! number of steady data points
             maxnumiterations = 100     ! maximum number of iterations

  real(dp) :: eta_s = 0._dp,          & ! solvent viscosity
              start_rate = 2.5e-2_dp, &  ! start data point
              end_rate = 2.5e-1_dp,   &  ! end data point
              factor = 2._dp,         &  ! factor last/first data point distance
              timestep = 1.e-3_dp,    &  ! time step for steady
              Lxx_add = 0.01_dp,      &  ! added value for Lxx (compressible))
              Jvol = 0.5_dp,          &  ! change in volume for model=24
              eps = 1.e-10_dp            ! Newton-Raphson accuracy

  real(dp), dimension(nmodes) :: &
               G,        & ! modulus
               lambda,   & ! relaxation time
               tau_ref,  & ! reference von Mises
               tau_ref1, & ! reference von Mises 1 RE
               tau_ref2, & ! reference von Mises 2 RE
               f1          ! weight factor for first term RE

  real(dp), allocatable, dimension(:) :: rate ! shear rate (for flowtype=1)
                                              ! strain rate (for flowtype=2,3)

  integer :: mode, iter, i, info, ipiv(ncompb)

  real(dp) :: tau(1,ncompt), gradv(1,ncompg)
  real(dp) :: b(1,ncompb,nmodes), rhs(1,ncompb,nmodes)
  real(dp) :: Hb(1,ncompb,nmodes), dHb(1,ncompb,ncompb,nmodes)
  real(dp) :: smat(ncompb,ncompb), fvec(ncompb)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemopt

  type(tdpar_t) :: tdpar

! namelist for input of variables; read from standard input

  namelist /comppar/ flowtype, includesolvent, eta_s, G, lambda, alam_model, &
                     tau_ref, tau_ref1, tau_ref2, f1, numdatapoints, &
                     start_rate, end_rate, factor, eps, maxnumiterations, &
                     method, timestep

  read ( unit=*, nml=comppar )

  USE_CONFORMATION_MODEL = .false.
  USE_EGP_STRESS_TENSOR_FORM = .false.

! set numerical parameters for time discretization

  tdpar%timestep = timestep
  tdpar%method = method

! set data points

  allocate ( rate(numdatapoints) )

  call distribute_elements ( numdatapoints - 1, rate, ratio=1, factor=factor )

  rate = start_rate + rate * (end_rate - start_rate)

! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, coorsys, &
    alam_model=alam_model )

  if ( model == 24 ) then
!   EGP compressible: input J
    call create_vemopt ( vemopt, dep_J=.true., compute_drhs=.true., &
      ncomp=ncompb, nmodes=nmodes )
    vemopt%J = Jvol
  else
    call create_vemopt ( vemopt, compute_drhs=.true., ncomp=ncompb, &
      nmodes=nmodes )
  end if

! set material parameters

  vemodel%modulus = G
  vemodel%lambda = lambda
  if ( alam_model == 4 ) then
    vemodel%alam(1,:) = tau_ref
  else if ( alam_model == 5 ) then
    vemodel%alam(1,:) = tau_ref1
    vemodel%alam(2,:) = f1
    vemodel%alam(3,:) = tau_ref2
  end if

! open output files

  open ( unit=13, file='iter.out', recl=300 )
  open ( unit=14, file='steady13.out', recl=300 )
  open ( unit=16, file='detb13.out', recl=300 )

! initialize b

  do mode = 1, nmodes
    b(1,:,mode) = [ 1, 0, 0, 1, 1 ]
  end do

! loop over data points

  do i = 1, numdatapoints

!   initialize velocity gradient

    select case ( flowtype )
    case(1)
      gradv(1,:) = [ 0._dp, rate(i), 0._dp, 0._dp ]
      if ( model == 24 ) gradv(1,1) = gradv(1,1) + Lxx_add ! EGP compressible trL/=0
    case(2)
      gradv(1,:) = [ rate(i), 0._dp, 0._dp, -rate(i) ]
    case(3)
      gradv(1,:) = [ rate(i), rate(i), 0._dp, -rate(i) ]
    case default
      write(*,'(/a,i0/)') 'Error: wrong value flowtype: ', flowtype
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

        call rhs_viscoelastic_2D_b ( vemodel, gradv, b, rhs, mode, mode, &
          vemopt=vemopt )

        call NRtd_steady_viscoelastic_2D_b ( vemodel, tdpar, b, Hb, dHb, &
          mode, mode )

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
        write(13,*) maxval(abs(fvec))

        if ( maxval(abs(fvec)) < eps ) exit

      end do

    end do

!   Output

    call stress_viscoelastic_2D_b ( vemodel, b, tau, vemopt=vemopt )

    call output ( rate(i), L=gradv(1,:), tau=tau(1,:) )

    write(16,*) rate(i), detb_2D(b(1,:,1))

  end do

  close(13)
  close(14)
  close(16)

! delete the model

  call delete ( vemodel )
  call delete ( vemopt )

  deallocate ( rate )

contains


! solvent stress
! NOTE: the following variables are inherited from the host:
!  ncompg, ncompt, eta_s

  function tau_s ( L )

!  velocity gradient
   real(dp), dimension(ncompg), intent(in) :: L
!  stress tensor
   real(dp), dimension(ncompt) :: tau_s

   tau_s(1) = 2 * eta_s * L(1)
   tau_s(2) = eta_s * ( L(2) + L(3) )
   tau_s(3) = 2 * eta_s * L(4)
   tau_s(4) = 0

  end function tau_s


! output
! NOTE: the following variables are inherited from the host:
!  nmodes, includesolvent

  subroutine output ( rate, L, tau )

    real(dp), intent(in) :: rate
!   velocity gradient
    real(dp), dimension(ncompg), intent(in) :: L
!   stress tensor
    real(dp), dimension(ncompt), intent(in) :: tau

    if ( nmodes == 1 ) then

      if ( includesolvent) then
        write(14,'(9es16.8,1x,i0)') rate, tau + tau_s(L), vonmises_2D(tau)
      else
        write(14,'(9es16.8,1x,i0)') rate, tau, vonmises_2D(tau)
      end if

    else

      if ( includesolvent) then
        write(14,'(7es16.8)') rate, tau + tau_s(L)
      else
        write(14,'(7es16.8)') rate, tau
      end if

    end if

  end subroutine output

end program steady13

