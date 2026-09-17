! Test for the derivative of the stress expression using the contravariant
! deformation tensor (b). 3D flow.
! The model is a Maxwell, Giesekus or FENE-P model.

program stresstest1

  use viscoelastic_models_3D_b_m

  implicit none

  integer, parameter :: model    = 3, &  ! 2: Maxwell 3: Giesekus 20: FENE-P
                        nmodes   = 1, &  ! number of modes
                        flowtype = 2, &  ! 3D
                        ncompb   = 9, & ! number of components of b
                        ncompt   = 6, & ! number of components of tau
                maxnumiterations = 100   ! maximum number of iterations

  real(dp), parameter :: modulus(nmodes) = [2._dp], & ! modulus
                         lambda(nmodes)  = [3._dp], & ! relaxation time
                         mobility(nmodes)  = 0.1_dp, & ! alpha
                         bpar(nmodes)  = 100.0_dp, & ! bparameter
                         c = 1._dp, & ! constant
                         eps = 1.e-10_dp   ! Newton-Raphson accuracy


  integer :: iter, i, info, ipiv(ncompb)
  real(dp) :: b(ncompb), tau(ncompt)
  real(dp) :: dtau(ncompt,ncompb)
  real(dp) :: smat(ncompb,ncompb), fvec(ncompb), f0(ncompb)

  type(vemodel_t) :: vemodel

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

! set f0

  f0 = [ 2, 1, 0, 1, 2, -1, 0, -1, 2 ]

! initialize b

  b = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ]

! open output file

  open ( unit=13, file='iter.out', recl=300 )
  open ( unit=14, file='out', recl=300 )

! Solve f(b) = 0 with
!
!    f(b) = tau(b) + c b - f0
!
!  using Newton-Raphson, with Jacobian
!
!    df/db = dtaudb + c I

  iter = 0

  do

    iter = iter + 1

    if ( iter > maxnumiterations ) then
      write(*,'(a,i0/)') &
        ' Maximum number of iterations reached = ', maxnumiterations
      stop
    end if

!   stress

    call stress_viscoelastic_3D_single_mode_b ( vemodel, b, tau, mode=1 )
    call dstress_viscoelastic_3D_single_mode_b ( vemodel, b, dtau, mode=1 )

!   solve J^-1 * residual

    fvec([1,2,3,5,6,9]) = tau
    fvec([4,7,8]) = tau([2,3,5])
    fvec = fvec + c * b - f0
    smat([1,2,3,5,6,9],:) = dtau
    smat([4,7,8],:) = dtau([2,3,5],:)
    do i = 1, ncompb
      smat(i,i) = smat(i,i) + c
    end do

    call dgesv( ncompb, 1, smat, ncompb, ipiv, fvec, ncompb, info )

    if ( info /= 0 ) then
      print *, 'Lapack DGESV: info = ', info
      stop
    end if

!   update for Newton-Raphson: un+1=un-J(un)^-1*residual(un)

    b = b - fvec

!    print *, maxval(abs(fvec))
    write(13,*) maxval(abs(fvec))
    write(14,*) iter, b

    if ( maxval(abs(fvec)) < eps ) exit

  end do

  close(13)
  close(14)

! delete the model

  call delete ( vemodel )

end program stresstest1
