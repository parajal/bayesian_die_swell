! Steady flow with a constant velocity gradient.
! Newton-Raphson iteration.
! The model is an evp model (Saramito 2009, SRM2) with Phan-Thien Tanner
! relaxation expression.
! Output is the viscoelastic (polymer) extra-stress.
! Optional output of eigenvalues and eigenvectors of Jacobian of rhs.

program steady4

  use meshgen_basic_common_m
  use viscoelastic_models_3D_m

  implicit none

  integer, parameter :: nmodes  = 1,  & ! number of modes
                        coorsys = 2,  & ! 3D Cartesian
                        alam_model = 3, &  ! adapted lambda model SRM2
                        ncompc  = 6,  & ! number of components of c
                        ncompt  = 6,  & ! number of components of tau
                        ndim    = 3     ! number of coordinate directions

  logical :: includesolvent = .false., & ! solvent included only if set .true.
             eigen = .false. ! compute eigenvalues/vectors of Jacobian of rhs

  integer :: model = 5,    &  ! 5: PTT linear model 6: PTT exponential
             flowtype = 1, &  ! flow type:
                              ! 1: shear
                              ! 2: planar extension
                              ! 3: uniaxial extension
             numdatapoints = 10, & ! number of steady data points
             maxnumiterations = 100     ! maximum number of iterations

  real(dp) :: eta_s = 0._dp,          &  ! solvent viscosity
              start_rate = 2.5e-2_dp, &  ! start data point
              end_rate = 2.5e-1_dp,   &  ! end data point
              factor = 2._dp,         &  ! factor last/first data point distance
              eps = 1.e-10_dp            ! Newton-Raphson accuracy

  real(dp), dimension(nmodes) :: &
               G,        & ! modulus
               lambda,   & ! relaxation time
               epsptt = 0.3_dp, & ! non-linear PTT parameter
               tau_y, &    ! yield stress
               Kfac    = 1, & ! K viscosity factor
               nexp    = 0.5_dp ! n power-law exponent

  real(dp), allocatable, dimension(:) :: rate ! shear rate (for flowtype=1)
                                              ! strain rate (for flowtype=2,3)

  integer :: mode, iter, i, info, ipiv(ncompc)

  real(dp) :: tau(1,ncompt), gradv(1,ndim,ndim)
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes)
  real(dp) :: smat(ncompc,ncompc), taumod(ncompt), fvec(ncompc)

  integer, parameter :: lwork = 780
  real(dp) :: wr(ncompc), wi(ncompc), vl(ncompc,ncompc), vr(ncompc,ncompc)
  real(dp) :: work(lwork)

  type(vemodel_t) :: vemodel
  type(vemopt_t) :: vemopt


! namelist for input of variables; read from standard input

  namelist /comppar/ flowtype, includesolvent, eta_s, model, G, lambda, &
                     epsptt, tau_y, Kfac, nexp, numdatapoints, &
                     start_rate, end_rate, factor, eps, maxnumiterations, eigen

  read ( unit=*, nml=comppar )

! create vemopt structure

  call create ( vemopt, compute_drhs=.true., ncomp=ncompc, nmodes=nmodes )

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
  vemodel%nonlin(1,:) = epsptt
  vemodel%alam(1:3,1) = [ tau_y, Kfac, nexp ]

! open output files

  open ( unit=13, file='iter.out', recl=300 )
  open ( unit=14, file='steady4.out', recl=300 )
  if ( eigen ) open ( unit=15, file='eigen4.out', recl=300 )

! initialize c, depending on flow type

  select case ( flowtype )
  case(1) ! shear
    do mode = 1, nmodes
      c(1,:,mode) = [ 1, 0, 0, 1, 0, 1 ] + &
                    [ 0, 2, 0, 0, 0, 0 ]*tau_y(mode)/G(mode)
    end do
  case(2) ! planar extension
    do mode = 1, nmodes
      c(1,:,mode) = [ 1, 0, 0, 1, 0, 1 ] + &
                    [ 2, 0, 0, 0, 0, 0 ]*tau_y(mode)/G(mode)
    end do
  case(3) ! uniaxial extension
    do mode = 1, nmodes
      c(1,:,mode) = [ 1, 0, 0, 1, 0, 1 ] + &
                    [ 3, 0, 0,-1, 0,-1 ]*tau_y(mode)/G(mode)
    end do
  case default
    write(*,'(/a,i0/)') 'Error: wrong value flowtype: ', flowtype
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

        call rhs_viscoelastic_3D ( vemodel, gradv, c, rhs, mode, mode, &
          vemopt=vemopt )


!       stress tensor, trace tau and von Mises for mode

        call stress_viscoelastic_3D_single_mode ( vemodel, c(1,:,mode), &
          taumod, mode )

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

    call stress_viscoelastic_3D ( vemodel, c, tau )

    call output ( rate(i), L=gradv(1,:,:), tau=tau(1,:) )

    if ( eigen ) call output_eigen ( rate(i) )

  end do

  close(13)
  close(14)

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
        write(14,'(9es16.8,1x,i0)') rate, tau + tau_s(L), vonmises_3D(tau)
      else
        write(14,'(9es16.8,1x,i0)') rate, tau, vonmises_3D(tau)
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

end program steady4
