! Equilibration or start-up of flow with thermal
! fluctuations with a constant velocity gradient.
! 3D single-point calculations
! Equilibrium, shear or uniaxial extensional flow
! The model is a UCM, FENE-P or Giesekus model
! The integration scheme is forward Euler
! Computes average stress, conformation,
! contravariant deformation, variance and
! autocorrelation functions over independent
! realizations.

program startup14

  use viscoelastic_models_3D_b_m
  use generalized_stokes_elements_m

  implicit none

  integer(si), parameter :: &
    seed = 10             ! seed for starting random-number
                          ! generator

  integer, parameter :: model   = 2, &     ! 2: UCM model
                                           ! 3: Giesekus model
                                           ! 20: FENE-P model
                        nmodes  = 1, &     ! number of modes
                        flowtype = 2, &    ! 3D Cartesian
                        flow    = 1, &     ! flow:
                                           ! 1: simple shear
                                           ! 2: uniaxial elongation
                        ncompb  = 9, &     ! number of components of b
                        ncompt  = 6, &     ! number of components of tau
                        ndim    = 3, &     ! number of coordinate directions
                        nsim = 1e3, &      ! number of realizations
                        nsteps = 1e5, &    ! number of time steps
                        correvery = 10     ! correlation function every
                                           ! correvery steps

  real(dp), parameter :: modulus = 1, &      ! modulus
                         lambda  = 1, &      ! relaxation time
!                         nonlinpar = 0.1, &  ! non-linear material parameter
!                                             ! (if model=3,20)
                         rate = 0._dp, &     ! strain rate
                         timestep = 1.e-2_dp, &  ! time step
                         theta = 0.1_dp      ! measure for the importance of
                                             ! fluctuations

  integer :: i, j, sim, step, ncorr
  real(dp) :: b(1,ncompb,nmodes), c(1,ncompt)!, c0
  real(dp) :: rhs(1,ncompb,nmodes), tau(1,ncompt), gradv(1,ndim,ndim)
  real(dp) :: dW1(1,ncompb), dW2(1,ncompb)
  real(dp) :: rhs_fluct(1,ncompb,nmodes), rhs_stoch(1,ncompb,nmodes)
  real(dp) :: bt(nsteps,3), ct(nsteps,4), taut(nsteps,4)
  real(dp) :: bmean(nsteps,3), cmean(nsteps,4), taumean(nsteps,4)
  real(dp) :: varc(nsteps,2), tauend(nsim,3)
  real(dp), dimension(:,:), allocatable :: acb, acc
  real(dp), dimension(:), allocatable :: ac

  type(vemodel_t) :: vemodel


! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype )

! set material parameters

  vemodel%modulus = modulus
  vemodel%lambda = lambda
!  vemodel%nonlin = nonlinpar

! initialize velocity gradient

  select case ( flow )
    case ( 1 )
      gradv(1,1,:) = [ 0._dp,  rate, 0._dp ]
      gradv(1,2,:) = [ 0._dp, 0._dp, 0._dp ]
      gradv(1,3,:) = [ 0._dp, 0._dp, 0._dp ]
    case ( 2 )
      gradv(1,1,:) = [  rate,   0._dp,   0._dp ]
      gradv(1,2,:) = [ 0._dp, -rate/2,   0._dp ]
      gradv(1,3,:) = [ 0._dp,   0._dp, -rate/2 ]
    case default
      write(*,'(/a,i0/)') 'Error: wrong value flow: ', flow
      stop
  end select

! initialize random-number generator

  call zigset ( seed )

! initialize time stepping

  ncorr = nsteps/correvery
  allocate( acb(ncorr,2), acc(ncorr,2), ac(ncorr) )

!  c0 = 1 + (ndim+1)*theta

  bmean = 0._dp
  cmean = 0._dp
  taumean = 0._dp
  acb = 0._dp
  acc = 0._dp
  varc = 0._dp


! perform nsim realizations

  do sim = 1, nsim

!   initialize b

    b(1,:,1) = [ 1, 0, 0, 0, 1, 0, 0, 0, 1 ]
!    b = sqrt(c0)*b


!   stepping using explicit Euler

    do step = 1, nsteps

!     create random numbers with normal distribution, variance 1
!     use ziggurat

      do i = 1, size(dW1,1)
        do j = 1, size(dW1,2)
          dW1(i,j) = rnor()
          if ( model == 3 ) dW2(i,j) = rnor()
        end do
      end do

!     right-hand side

      call rhs_viscoelastic_3D_b ( vemodel, gradv, b, rhs )

      if ( model == 3 ) then
        call rhs_viscoelastic_fluct_3D_b ( vemodel, b, rhs_fluct, &
          rhs_stoch, dW1, dW2 )
      else
        call rhs_viscoelastic_fluct_3D_b ( vemodel, b, rhs_fluct, &
          rhs_stoch, dW1 )
      end if

!     do step

      b = b + ( rhs + theta * rhs_fluct ) * timestep

      b = b + sqrt ( theta ) * rhs_stoch * sqrt ( timestep )

!     compute conformation and stress

      call conformation_3D_b ( b(:,:,1), c )

      call stress_viscoelastic_3D_b ( vemodel, b, tau )

!     store b, c and tau in matrices for further processing

      bt(step,:) = [ step*timestep, b(1,1,1), b(1,2,1) ]
      ct(step,:) = [ step*timestep, c(1,1), c(1,2), c(1,4) ]
      taut(step,:) = [ step*timestep, tau(1,1), tau(1,2), tau(1,4) ]

    end do

!   compute averages and variance

    bmean = bmean + bt/nsim
    cmean = cmean + ct/nsim
    varc = varc + ct(:,2:3)**2/nsim
    taumean = taumean + taut/nsim
    tauend(sim,:) = taut(size(taut,1),2:size(taut,2))

!   compute autocorrelation functions for b and c

    call correlation ( bt(1:nsteps:correvery,2), ac )
    acb(:,1) = acb(:,1) + ac/nsim
    call correlation ( bt(1:nsteps:correvery,3), ac )
    acb(:,2) = acb(:,2) + ac/nsim
    call correlation ( ct(1:nsteps:correvery,2), ac )
    acc(:,1) = acc(:,1) + ac/nsim
    call correlation ( ct(1:nsteps:correvery,3), ac )
    acc(:,2) = acc(:,2) + ac/nsim

  end do

! compute variance of c

  varc = varc - cmean(:,2:3)**2

! write averages to files

  open(unit=11, recl=300, file='b.out', &
     status='replace')
  open(unit=13, recl=300, file='acb.out', &
     status='replace')
  open(unit=12, recl=300, file='c.out', &
     status='replace')
  open(unit=14, recl=300, file='acc.out', &
     status='replace')
  open(unit=15, recl=300, file='varc.out', &
     status='replace')
  open(unit=16, recl=300, file='tau.out', &
     status='replace')
  open(unit=17, recl=300, file='tauend.out', &
     status='replace')
  do i=1,nsteps
    write(11, fmt=*) bmean(i,:)
    write(12, fmt=*) cmean(i,:)
    write(15, fmt=*) cmean(i,1), varc(i,:)
    write(16, fmt=*) taumean(i,:)
  end do
  do i=1,nsim
    write(17, fmt=*) tauend(i,:)
  end do
  do i=1,ncorr
    write(13, fmt=*) bmean((i-1)*correvery+1,1)-timestep, acb(i,:)
    write(14, fmt=*) cmean((i-1)*correvery+1,1)-timestep, acc(i,:)
  end do
  close(unit=11)
  close(unit=13)
  close(unit=12)
  close(unit=14)
  close(unit=15)
  close(unit=16)
  close(unit=17)

! delete the model

  call delete ( vemodel )


contains

  subroutine correlation ( A, C )

    real(dp), dimension(:), intent(in) :: A
    real(dp), dimension(:), intent(out) :: C
    real(dp) :: C0
    integer :: i, n

    n = size(A)
    C0 = sum(A**2)/n - (sum(A)/n)**2

    do i = 1, size(C)

      C(i) = sum(A(i:n)*A(1:n-i+1))/(n-i+1) - sum(A(i:n))/(n-i+1) * &
        sum(A(1:n-i+1))/(n-i+1)

    end do

    C = C/C0

  end subroutine correlation

end program startup14
