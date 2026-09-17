! Steady stress tensor of an integral model for a range of rates in shear and
! extension.

program df_stress1

  use tfem_m
  use dfm_elements_m

  implicit none

  integer, parameter :: &
    type_of_model = 0 , & ! 0: separable Rivlin-Sawyers,
                          ! 1: non-separable RS (not yet available)
    spectrum = 0,       & ! 0: discrete spectrum
                          ! 1: Mittag-Leffler (not yet available)
                          ! 2: discrete+Carreau
    numdismodes = 0,    & ! number of discrete modes (for spectrum=2 only).
    nummodes = 2,       & ! number of modes in the spectrum
    dampfunc = 2          ! damping function: 1: Lodge, 2: McKinley, 3: PSM

  integer, parameter :: &
    flowtype = 1, &  ! flow type:
                     ! 1: shear
                     ! 2: planar extension
                     ! 3: uniaxial extension
    numdatapoints = 120   ! number of steady data points

  real(dp) :: &
    start_rate = 1.e-1_dp, & ! start data point
    end_rate = 1.e3_dp,   &  ! end data point
    factor = 1.e4_dp         ! factor last/first data point distance

  real(dp), parameter :: &
    eta_s = 0.0_dp,  & ! solvent viscosity
    G(nummodes) =      [  5._dp, 1._dp ],  & ! modulus
    lambda(nummodes) = [ 0.2_dp, 1._dp ],  & ! relaxation time
    nexpo(nummodes) =  [  2._dp, 2._dp ],  & ! power-law exponent
                                             ! (for spectrum=2)
    a = 0.0001_dp,        & ! parameter in McKinley damping function
    alpha_PSM = 2.5e4_dp, & ! alpha parameter in PSM damping function
    beta_PSM = 0.25_dp      ! beta parameter in PSM damping function

  integer, parameter :: &
    fintpl = 4,         & ! Q1 interpolation in space for F (fake)
    fintpltau = 1,      & ! P1 interpolation in age (tau)
!    fintpltau = 2,      & ! P2 interpolation in age (tau)
    nintvaltau = 200,   & ! number of age intervals
    ncompf = 9,         & ! number of components of deformation tensor F
    ninttau = fintpltau+1    ! number of fields within an interval

  real(dp), parameter :: &
    dtau1 = 3*minval(lambda)/nintvaltau, & ! length of first interval in
                                           ! age (tau) direction
    tauc = 15._dp*maxval(lambda)     ! maximum age (cutoff)

  type(coefficients_t) :: coefficients

  integer :: i, j, k
  real(dp) :: x(0:nintvaltau), stress_ten(3,3), val
  real(dp) :: tau_nodal_points(ninttau,nintvaltau)
  real(dp) :: F_vec(ncompf,ninttau,nintvaltau)
  real(dp), dimension(numdatapoints) :: rate ! shear rate (for flowtype=1)
                                             ! strain rate (for flowtype=2,3)


! set data points

  call distribute_elements ( numdatapoints - 1, rate, ratio=1, factor=factor )

  rate = start_rate + rate * (end_rate - start_rate)

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=400, ncoefr=350, &
    ncoefra=[nintvaltau+1,nummodes,nummodes,nummodes] )

  coefficients%i = 0

  coefficients%i(351:353) = [ fintpl, fintpltau, nintvaltau ]

  coefficients%i(363) = type_of_model
  coefficients%i(364) = spectrum
  coefficients%i(365) = nummodes
  coefficients%i(366) = dampfunc
  coefficients%i(367) = numdismodes

  coefficients%r = 0

  coefficients%r(303) = a
  coefficients%r(304) = alpha_PSM
  coefficients%r(305) = beta_PSM

  call distribute_elements ( nintvaltau, x, ratio=7, factor=dtau1/tauc )

  coefficients%ra(1)%a = tauc * x(0:nintvaltau)
  coefficients%ra(2)%a = G
  coefficients%ra(3)%a = lambda
  coefficients%ra(4)%a = nexpo

  write(*,'(/a,g0.4/)') ' non-equidistant grid; last interval = ', &
     coefficients%ra(1)%a(nintvaltau+1)-coefficients%ra(1)%a(nintvaltau)

  call dfm_tau_nodal_points ( coefficients, tau_nodal_points )

! open output file

  open ( unit=14, file='steady.out', recl=300 )

  do k = 1, numdatapoints

    select case ( flowtype )
    case(1)
      F_vec = 0
      F_vec(1,:,:) = 1 ! steady F_xx
      F_vec(5,:,:) = 1 ! steady F_yy
      F_vec(9,:,:) = 1 ! steady F_zz
      do j = 1, nintvaltau
        do i = 1, ninttau
          F_vec(2,i,j) = rate(k) * tau_nodal_points(i,j) ! steady F_xy
        end do
      end do
    case(2)
      F_vec = 0
      F_vec(9,:,:) = 1 ! steady F_zz
      do j = 1, nintvaltau
        do i = 1, ninttau
          val = min( 100._dp, rate(k) * tau_nodal_points(i,j) ) ! avoid overflow
          F_vec(1,i,j) = exp(val)  ! steady F_xx
          F_vec(5,i,j) = exp(-val) ! steady F_yy
        end do
      end do
    case(3)
      F_vec = 0
      do j = 1, nintvaltau
        do i = 1, ninttau
          val = min( 150._dp, rate(k) * tau_nodal_points(i,j) ) ! avoid overflow
          F_vec(1,i,j) = exp(val)    ! steady F_xx
          F_vec(5,i,j) = exp(-val/2) ! steady F_yy
          F_vec(9,i,j) = exp(-val/2) ! steady F_zz
        end do
      end do
    case default
      write(*,'(/a,i0/)') 'Error: wrong value flowtype: ', flowtype
      stop
    end select

    call dfm_stress_tensor_point ( coefficients, F_vec, stress_ten )

    select case ( flowtype )
    case(1)
      stress_ten(1,2) = stress_ten(1,2) + eta_s * rate(k)
      stress_ten(2,1) = stress_ten(2,1) + eta_s * rate(k)
    case(2)
      stress_ten(1,1) = stress_ten(1,1) + 2 * eta_s * rate(k)
      stress_ten(2,2) = stress_ten(2,2) - 2 * eta_s * rate(k)
    case(3)
      stress_ten(1,1) = stress_ten(1,1) + 2 * eta_s * rate(k)
      stress_ten(2,2) = stress_ten(2,2) - 2 * eta_s * rate(k) / 2
      stress_ten(3,3) = stress_ten(2,2) - 2 * eta_s * rate(k) / 2
    case default
      write(*,'(/a,i0/)') 'Error: wrong value flowtype: ', flowtype
      stop
    end select

    write(14,'(10es16.8)') rate(k), stress_ten

  end do

end program df_stress1
