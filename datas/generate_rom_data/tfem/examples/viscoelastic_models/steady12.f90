! Steady flow with a constant velocity gradient. 2D flow.
! Newton-Raphson iteration.
! The model is a EGP elastoviscoplastic model for incompressible flow
! Output is the viscoelastic (polymer) extra-stress.
! Optional output of eigenvalues and eigenvectors of Jacobian of rhs.

program steady12

  use meshgen_basic_common_m
  use viscoelastic_models_2D_m
  use limits_m, only: USE_EGP_STRESS_TENSOR_FORM

  implicit none

  integer, parameter :: nmodes  = 1,  & ! number of modes
                        coorsys = 0,  & ! 3D Cartesian
                        ncompc  = 4,  & ! number of components of c
                        ncompt  = 4,  & ! number of components of tau
                        ncompg  = 4     ! number of components of gradv

  logical :: includesolvent = .false., & ! solvent included only if set .true.
             eigen = .false. ! compute eigenvalues/vectors of Jacobian of rhs

  integer :: model = 24,      & ! 23: EGP model incompressible
                                ! 24: EGP model compressible
             alam_model = 4,  & ! adapted lambda model:
                                ! 4: Eyring, 5: Ree-Eyring
             flowtype = 1, &  ! flow type:
                              ! 1: shear
                              ! 2: planar extension
                              ! 3: combined
             numdatapoints = 10, & ! number of steady data points
             maxnumiterations = 100     ! maximum number of iterations

  real(dp) :: eta_s = 0._dp,          & ! solvent viscosity
              start_rate = 2.5e-2_dp, &  ! start data point
              end_rate = 2.5e-1_dp,   &  ! end data point
              factor = 2._dp,         &  ! factor last/first data point distance
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

  integer :: mode, iter, i, info, ipiv(ncompc)

  real(dp) :: tau(1,ncompt), gradv(1,ncompg)
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes)
  real(dp) :: smat(ncompc,ncompc), fvec(ncompc)

  integer, parameter :: lwork = 780
  real(dp) :: wr(ncompc), wi(ncompc), vl(ncompc,ncompc), vr(ncompc,ncompc)
  real(dp) :: work(lwork)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemopt

! namelist for input of variables; read from standard input

  namelist /comppar/ flowtype, includesolvent, eta_s, G, lambda, alam_model, &
                     tau_ref, tau_ref1, tau_ref2, f1, numdatapoints, &
                     start_rate, end_rate, factor, eps, maxnumiterations, eigen

  read ( unit=*, nml=comppar )

  USE_EGP_STRESS_TENSOR_FORM = .false.

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
      ncomp=ncompc, nmodes=nmodes )
    vemopt%J = Jvol
  else
    call create_vemopt ( vemopt, compute_drhs=.true., ncomp=ncompc, &
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
  open ( unit=14, file='steady12.out', recl=300 )
  if ( eigen ) open ( unit=15, file='eigen12.out', recl=300 )
  open ( unit=16, file='detc12.out', recl=300 )

! initialize c

  do mode = 1, nmodes
    c(1,:,mode) = [ 1, 0, 1, 1 ]
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

        call rhs_viscoelastic_2D ( vemodel, gradv, c, rhs, mode, mode, &
          vemopt=vemopt )

!       solve J^-1 * residual

        smat = - vemopt%drhs(1,:,:,mode)
        fvec = - rhs(1,:,mode)

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

!   Output

    call stress_viscoelastic_2D ( vemodel, c, tau, vemopt=vemopt )

    call output ( rate(i), L=gradv(1,:), tau=tau(1,:) )

    if ( eigen ) call output_eigen ( rate(i) )

    write(16,*) rate(i), detc_2D(c(1,:,1))

  end do

  close(13)
  close(14)
  if ( eigen ) close(15)
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

! output eigenvalues and eigenvector of the Jacobian of the rhs

  subroutine output_eigen ( rate )

    real(dp), intent(in) :: rate
    integer :: j

!   solve (right) eigen problem

    call dgeev( 'N', 'V', ncompc, -smat, ncompc, wr, wi, vl, ncompc, vr, &
      ncompc, work, lwork, info )

    write(15,'(a,g0)') 'rate = ', rate
    write(15,'(a,i0)') 'info = ', info
    if ( info == 0 ) &
              write(15,'(a,i0)') 'optimal lwork: work(1) = ', nint(work(1))
    write(15,'(a)') ' Eigenvalues: j, wr(j), wi(j) '
    do j = 1, ncompc
      write(15,'(i3,2es16.8)') j, wr(j), wi(j)
    end do
    write(15,'(a)') ' Eigenvectors: j, v(j) '
    do j = 1, ncompc
      write(15,'(i3,6es16.8)') j, vr(:,j)
    end do

  end subroutine output_eigen

end program steady12

