! Viscoelastic data for various materials

module material_data_m

  use viscoelastic_models_defs_m
  use element_defs_m, only: coefficients_t, check

  implicit none

! interface for generic fill subroutine

  interface fill_data
    module procedure fill_data_v, fill_data_c
  end interface fill_data

contains


  subroutine fill_data_v ( vemodel, data )

    type(vemodel_t), intent(inout) :: vemodel
    character(len=*), intent(in) :: data

    select case ( data )
    case('PMMA')
      call fill_PMMA_data_v ( vemodel )
    case('iPP')
      call fill_iPP_data_v ( vemodel )
    case('PS')
      call fill_PS_data_v ( vemodel )
    case('PLLA')
      call fill_PLLA_data_v ( vemodel )
    case default
      write(*,'(/a/2a/)') &
        'Error in fill_data_v:', &
        '  unknown data: ', data
      stop
    end select

  end subroutine fill_data_v


! fill data for coefficients

  subroutine fill_data_c ( coefficients, data, vemodel )

    type(coefficients_t), intent(inout) :: coefficients
    character(len=*), intent(in) :: data

!   If present: output of a structure of type vemodel_t that contains
!   the filled data in a more accessible format.
!   WARNING: It should only be used to access the filled data! All other
!   components might not be correct, such as flowtype, bvariant, vmglobal
!   and norlxmode.
    type(vemodel_t), optional, intent(out) :: vemodel

    type(vemodel_t) :: lvemodel

    call create_vemodel_without_data ( coefficients, lvemodel )
    call fill_data_v ( lvemodel, data )
    call fill_data_in_coefficients ( lvemodel, coefficients )

    if ( present(vemodel) ) vemodel = lvemodel

  end subroutine fill_data_c


! Data for PMMA from van Breemen et al. (2012) DOI: 10.1002/polb.23199

  subroutine fill_PMMA_data_v ( vemodel )

    type(vemodel_t), intent(inout) :: vemodel

    integer, parameter :: am = 12, & ! 1:am range of modes for alpha process
                          em = 13, & ! elastic mode
                          bm = 14    ! mode of beta process

    real(dp) :: T, kB, Vstar, eta0ref, tau_ref, Sa, r0, r1, r2, Gr, Gtot, &
      mu, betaJ, kappa

!   check vemodel

    if ( .not. vemodel%created ) then
      write(*,'(/2(a/))') &
        'Error in fill_PMMA_data_v:', &
        '  vemodel has not been created.'
      stop
    end if

!   check number of modes

    if ( .not. any( vemodel%nmodes == [em,bm] ) ) then
      write(*,'(/a/,2(a,i0)/a,i0/)') &
        'Error in fill_PMMA_data_v:', &
        '  number of modes must be ', em, ' or ', bm, &
        '  nmodes = ', vemodel%nmodes
      stop
    end if

!   set material parameters

    T = 295.15_dp  ! K  (22 degr C)
    kB = 1.380649e-23_dp  ! m^2 kg s^-2 K^-1
    kappa = 3.e3_dp ! Mpa

!   alpha process

    vemodel%modulus(1:am) = &
      [ 1.08e1_dp, 5.69e2_dp, 1.13e2_dp, 2.87e1_dp, 4.86e1_dp, &
        5.53e1_dp, 3.10e1_dp, 2.78e1_dp, 2.78e1_dp, 4.00e1_dp, &
        2.38e1_dp, 3.71e1_dp ] ! MPa

    vemodel%lambda(1:am) = &
      [ 1.73e4_dp,  1.55e4_dp, 7.50e2_dp, 1.48e2_dp, 1.19e2_dp, &
        2.91e1_dp,  9.87e0_dp, 4.64e0_dp, 2.17e0_dp, 7.33e-1_dp, &
        1.96e-1_dp, 8.41e-2_dp ] ! s

    Vstar = 1.60e-27_dp  ! m^3
    Sa = 7.4_dp  ! [-]
    r0 = 0.96_dp ! [-]
    r1 = 20.0_dp ! [-]
    r2 = -2.0_dp ! [-]
    mu = 0.2_dp ! [-]

    tau_ref = kB * T / Vstar * 1e-6_dp ! MPa
    betaJ = mu * kappa / tau_ref ! [-]

    vemodel%alam(1,1:am) = tau_ref

    if ( vemodel%alam_gammap_model == 1 ) then
      vemodel%alam_gammap(1,1:am) = Sa ! state parameter Sa
      vemodel%alam_gammap(2,1:am) = r0 ! fitting parameter r0
      vemodel%alam_gammap(3,1:am) = r1 ! fitting parameter r1
      vemodel%alam_gammap(4,1:am) = r2 ! fitting parameter r2
    end if

    if ( vemodel%alamJ_model == 1 ) then
      vemodel%alamJ(1,1:am) = betaJ ! exponent pressure dependence
    end if

!   elastic mode

    Gr = 21.0_dp ! MPa

    vemodel%modulus(em) = Gr
    vemodel%lambda(em) = 1.0_dp ! lambda for determinant stabilization only
    vemodel%alam(1,em) = 0._dp ! dummy tau_ref value for elastic mode

    if ( vemodel%alam_gammap_model == 1 ) then
      vemodel%alam_gammap(:,em) = 0._dp ! dummy softening parameters
    end if

    if ( vemodel%alamJ_model == 1 ) then
      vemodel%alamJ(1,em) = 0._dp ! dummy exponent pressure dependence
    end if

!   beta process

    if ( vemodel%nmodes == bm ) then

      Gtot = 1.13e2_dp ! MPa
      eta0ref = 4.13e-1_dp ! MPa.s
      Vstar = 1.49e-27_dp  ! m^3
      Sa = 7.4_dp  ! [-]
      r0 = 0.96_dp ! [-]
      r1 = 20.0_dp ! [-]
      r2 = -2.0_dp ! [-]
      mu = 0.2_dp ! [-]

      tau_ref = kB * T / Vstar * 1e-6_dp ! MPa
      betaJ = mu * kappa / tau_ref ! [-]

      vemodel%modulus(bm) = Gtot
      vemodel%lambda(bm) = eta0ref / Gtot ! s
      vemodel%alam(1,bm) = tau_ref

      if ( vemodel%alam_gammap_model == 1 ) then
        vemodel%alam_gammap(1,bm) = Sa ! state parameter Sa
        vemodel%alam_gammap(2,bm) = r0 ! fitting parameter r0
        vemodel%alam_gammap(3,bm) = r1 ! fitting parameter r1
        vemodel%alam_gammap(4,bm) = r2 ! fitting parameter r2
      end if

      if ( vemodel%alamJ_model == 1 ) then
        vemodel%alamJ(1,bm) = betaJ ! exponent pressure dependence
      end if

    end if

  end subroutine fill_PMMA_data_v


! Data for iPP from van Breemen et al. (2012) DOI: 10.1002/polb.23199

  subroutine fill_iPP_data_v ( vemodel )

    type(vemodel_t), intent(inout) :: vemodel

    integer, parameter :: am = 12, & ! 1:am range of modes for alpha process
                          em = 13, & ! elastic mode
                          bm = 14    ! mode of beta process

    real(dp) :: T, kB, Vstar, eta0ref, tau_ref, Sa, r0, r1, r2, Gr, Gtot, &
      mu, betaJ, kappa

!   check vemodel

    if ( .not. vemodel%created ) then
      write(*,'(/2(a/))') &
        'Error in fill_iPP_data_v:', &
        '  vemodel has not been created.'
      stop
    end if

!   check number of modes

    if ( .not. any( vemodel%nmodes == [em,bm] ) ) then
      write(*,'(/a/,2(a,i0)/a,i0/)') &
        'Error in fill_iPP_data_v:', &
        '  number of modes must be ', em, ' or ', bm, &
        '  nmodes = ', vemodel%nmodes
      stop
    end if

!   set material parameters

    T = 293.15_dp  ! K  (20 degr C)
    kB = 1.380649e-23_dp  ! m^2 kg s^-2 K^-1
    kappa = 1.66e3_dp ! Mpa

!   alpha process

    vemodel%modulus(1:am) = &
      [ 1.71e1_dp, 1.45e2_dp, 3.77e1_dp, 2.41e1_dp, 5.80e0_dp, &
        1.31e1_dp, 1.26e1_dp, 1.22e1_dp, 4.94e0_dp, 1.17e1_dp, &
        1.57e1_dp, 2.35e1_dp ] ! MPa
    vemodel%lambda(1:am) = &
      [ 2.51e8_dp,  2.20e8_dp, 1.67e7_dp, 2.85e6_dp, 7.82e5_dp, &
        5.49e5_dp,  1.69e5_dp, 5.87e4_dp, 2.28e4_dp, 1.27e4_dp, &
        3.09e3_dp,  4.70e2_dp ] ! s

    Vstar = 3.35e-27_dp  ! m^3
    Sa = 4.9_dp  ! [-]
    r0 = 0.955_dp ! [-]
    r1 = 2.0_dp ! [-]
    r2 = -1.0_dp ! [-]
    mu = 0.0_dp ! [-]

    tau_ref = kB * T / Vstar * 1e-6_dp ! MPa
    betaJ = mu * kappa / tau_ref ! [-]

    vemodel%alam(1,1:am) = tau_ref

    if ( vemodel%alam_gammap_model == 1 ) then
      vemodel%alam_gammap(1,1:am) = Sa ! state parameter Sa
      vemodel%alam_gammap(2,1:am) = r0 ! fitting parameter r0
      vemodel%alam_gammap(3,1:am) = r1 ! fitting parameter r1
      vemodel%alam_gammap(4,1:am) = r2 ! fitting parameter r2
    end if

    if ( vemodel%alamJ_model == 1 ) then
      vemodel%alamJ(1,1:am) = betaJ ! exponent pressure dependence
    end if

!   elastic mode

    Gr = 1.2_dp ! MPa

    vemodel%modulus(em) = Gr
    vemodel%lambda(em) = 1.0_dp ! lambda for determinant stabilization only
    vemodel%alam(1,em) = 0._dp ! dummy tau_ref value for elastic mode

    if ( vemodel%alam_gammap_model == 1 ) then
      vemodel%alam_gammap(:,em) = 0._dp ! dummy softening parameters
    end if

    if ( vemodel%alamJ_model == 1 ) then
      vemodel%alamJ(1,em) = 0._dp ! dummy exponent pressure dependence
    end if

!   beta process

    if ( vemodel%nmodes == bm ) then

      Gtot = 3.24e2_dp ! MPa
      eta0ref = 1.91e-5_dp ! MPa.s
      Vstar = 7.44e-27_dp  ! m^3
      Sa = 21._dp  ! [-]
      r0 = 1.0_dp ! [-]
      r1 = 0.2_dp ! [-]
      r2 = -0.3_dp ! [-]
      mu = 0.0_dp ! [-]

      tau_ref = kB * T / Vstar * 1e-6_dp ! MPa
      betaJ = mu * kappa / tau_ref ! [-]

      vemodel%modulus(bm) = Gtot
      vemodel%lambda(bm) = eta0ref / Gtot ! s
      vemodel%alam(1,bm) = tau_ref

      if ( vemodel%alam_gammap_model == 1 ) then
        vemodel%alam_gammap(1,bm) = Sa ! state parameter Sa
        vemodel%alam_gammap(2,bm) = r0 ! fitting parameter r0
        vemodel%alam_gammap(3,bm) = r1 ! fitting parameter r1
        vemodel%alam_gammap(4,bm) = r2 ! fitting parameter r2
      end if

      if ( vemodel%alamJ_model == 1 ) then
        vemodel%alamJ(1,bm) = betaJ ! exponent pressure dependence
      end if

    end if

  end subroutine fill_iPP_data_v


! Data for PS from van Breemen et al. (2012) DOI: 10.1002/polb.23199

  subroutine fill_PS_data_v ( vemodel )

    type(vemodel_t), intent(inout) :: vemodel

    integer, parameter :: am = 1, & ! 1:am range of modes for alpha process
                          em = 2, & ! elastic mode
                          bm = 3    ! mode of beta process

    real(dp) :: T, kB, Vstar, eta0ref, tau_ref, Sa, r0, r1, r2, Gr, Gtot, &
      mu, betaJ, kappa

!   check vemodel

    if ( .not. vemodel%created ) then
      write(*,'(/2(a/))') &
        'Error in fill_PS_data_v:', &
        '  vemodel has not been created.'
      stop
    end if

!   check number of modes

    if ( .not. any( vemodel%nmodes == [em,bm] ) ) then
      write(*,'(/a/,2(a,i0)/a,i0/)') &
        'Error in fill_PS_data_v:', &
        '  number of modes must be ', em, ' or ', bm, &
        '  nmodes = ', vemodel%nmodes
      stop
    end if

!   set material parameters

    T = 296.15_dp  ! K  (23 degr C)
    kB = 1.380649e-23_dp  ! m^2 kg s^-2 K^-1
    kappa = 3.5e3_dp ! Mpa

!   alpha process

    Gtot = 5.50e2_dp ! MPa
    eta0ref = 2.7e12_dp ! MPa.s

    Vstar = 4.17e-27_dp  ! m^3
    Sa = 14.0_dp  ! [-]
    r0 = 0.99_dp ! [-]
    r1 = 50.0_dp ! [-]
    r2 = -3.0_dp ! [-]
    mu = 0.14_dp ! [-]

    tau_ref = kB * T / Vstar * 1e-6_dp ! MPa
    betaJ = mu * kappa / tau_ref ! [-]

    vemodel%modulus(1:am) = [ Gtot ]
    vemodel%lambda(1:am) = [ eta0ref / Gtot ] ! s
    vemodel%alam(1,1:am) = tau_ref

    if ( vemodel%alam_gammap_model == 1 ) then
      vemodel%alam_gammap(1,1:am) = Sa ! state parameter Sa
      vemodel%alam_gammap(2,1:am) = r0 ! fitting parameter r0
      vemodel%alam_gammap(3,1:am) = r1 ! fitting parameter r1
      vemodel%alam_gammap(4,1:am) = r2 ! fitting parameter r2
    end if

    if ( vemodel%alamJ_model == 1 ) then
      vemodel%alamJ(1,1:am) = betaJ ! exponent pressure dependence
    end if

!   elastic mode

    Gr = 8.0_dp ! MPa

    vemodel%modulus(em) = Gr
    vemodel%lambda(em) = 1.0_dp ! lambda for determinant stabilization only
    vemodel%alam(1,em) = 0._dp ! dummy tau_ref value for elastic mode

    if ( vemodel%alam_gammap_model == 1 ) then
      vemodel%alam_gammap(:,em) = 0._dp ! dummy softening parameters
    end if

    if ( vemodel%alamJ_model == 1 ) then
      vemodel%alamJ(1,em) = 0._dp ! dummy exponent pressure dependence
    end if

!   beta process

    if ( vemodel%nmodes == bm ) then

      Gtot = 5.50e2_dp ! MPa
      eta0ref = 1.0e-1_dp ! MPa.s
      Vstar = 2.65e-27_dp  ! m^3
      Sa = 14._dp  ! [-]
      r0 = 0.99_dp ! [-]
      r1 = 50.0_dp ! [-]
      r2 = -4.0_dp ! [-]
      mu = 0.14_dp ! [-]

      tau_ref = kB * T / Vstar * 1e-6_dp ! MPa
      betaJ = mu * kappa / tau_ref ! [-]

      vemodel%modulus(bm) = Gtot
      vemodel%lambda(bm) = eta0ref / Gtot ! s
      vemodel%alam(1,bm) = tau_ref

      if ( vemodel%alam_gammap_model == 1 ) then
        vemodel%alam_gammap(1,bm) = Sa ! state parameter Sa
        vemodel%alam_gammap(2,bm) = r0 ! fitting parameter r0
        vemodel%alam_gammap(3,bm) = r1 ! fitting parameter r1
        vemodel%alam_gammap(4,bm) = r2 ! fitting parameter r2
      end if

      if ( vemodel%alamJ_model == 1 ) then
        vemodel%alamJ(1,bm) = betaJ ! exponent pressure dependence
      end if

    end if

  end subroutine fill_PS_data_v


! Data for PLLA from van Breemen et al. (2012) DOI: 10.1002/polb.23199

  subroutine fill_PLLA_data_v ( vemodel )

    type(vemodel_t), intent(inout) :: vemodel

    integer, parameter :: am = 1, & ! 1:am range of modes for alpha process
                          em = 2, & ! elastic mode
                          bm = 3    ! mode of beta process

    real(dp) :: T, kB, Vstar, eta0ref, tau_ref, Sa, r0, r1, r2, Gr, Gtot, &
      mu, betaJ, kappa

!   check vemodel

    if ( .not. vemodel%created ) then
      write(*,'(/2(a/))') &
        'Error in fill_PLLA_data_v:', &
        '  vemodel has not been created.'
      stop
    end if

!   check number of modes

    if ( .not. any( vemodel%nmodes == [em,bm] ) ) then
      write(*,'(/a/,2(a,i0)/a,i0/)') &
        'Error in fill_PLLA_data_v:', &
        '  number of modes must be ', em, ' or ', bm, &
        '  nmodes = ', vemodel%nmodes
      stop
    end if

!   set material parameters

    T = 295.15_dp  ! K  (22 degr C)
    kB = 1.380649e-23_dp  ! m^2 kg s^-2 K^-1
    kappa = 3.5e3_dp ! Mpa

!   alpha process

    Gtot = 5.50e2_dp ! MPa
    eta0ref = 1.0e15_dp ! MPa.s

    Vstar = 4.75e-27_dp  ! m^3
    Sa = 12.25_dp  ! [-]
    r0 = 0.99_dp ! [-]
    r1 = 20.0_dp ! [-]
    r2 = -5.0_dp ! [-]
    mu = 0.0_dp ! [-]

    tau_ref = kB * T / Vstar * 1e-6_dp ! MPa
    betaJ = mu * kappa / tau_ref ! [-]

    vemodel%modulus(1:am) = [ Gtot ]
    vemodel%lambda(1:am) = [ eta0ref / Gtot ] ! s
    vemodel%alam(1,1:am) = tau_ref

    if ( vemodel%alam_gammap_model == 1 ) then
      vemodel%alam_gammap(1,1:am) = Sa ! state parameter Sa
      vemodel%alam_gammap(2,1:am) = r0 ! fitting parameter r0
      vemodel%alam_gammap(3,1:am) = r1 ! fitting parameter r1
      vemodel%alam_gammap(4,1:am) = r2 ! fitting parameter r2
    end if

    if ( vemodel%alamJ_model == 1 ) then
      vemodel%alamJ(1,1:am) = betaJ ! exponent pressure dependence
    end if

!   elastic mode

    Gr = 3.45_dp ! MPa

    vemodel%modulus(em) = Gr
    vemodel%lambda(em) = 1.0_dp ! lambda for determinant stabilization only
    vemodel%alam(1,em) = 0._dp ! dummy tau_ref value for elastic mode

    if ( vemodel%alam_gammap_model == 1 ) then
      vemodel%alam_gammap(:,em) = 0._dp ! dummy softening parameters
    end if

    if ( vemodel%alamJ_model == 1 ) then
      vemodel%alamJ(1,em) = 0._dp ! dummy exponent pressure dependence
    end if

!   beta process

    if ( vemodel%nmodes == bm ) then

      Gtot = 5.50e2_dp ! MPa
      eta0ref = 3.5e1_dp ! MPa.s
      Vstar = 1.75e-27_dp  ! m^3
      Sa = 12.25_dp  ! [-]
      r0 = 0.99_dp ! [-]
      r1 = 100.0_dp ! [-]
      r2 = -12.5_dp ! [-]
      mu = 0.0 ! [-]

      tau_ref = kB * T / Vstar * 1e-6_dp ! MPa
      betaJ = mu * kappa / tau_ref ! [-]

      vemodel%modulus(bm) = Gtot
      vemodel%lambda(bm) = eta0ref / Gtot ! s
      vemodel%alam(1,bm) = tau_ref

      if ( vemodel%alam_gammap_model == 1 ) then
        vemodel%alam_gammap(1,bm) = Sa ! state parameter Sa
        vemodel%alam_gammap(2,bm) = r0 ! fitting parameter r0
        vemodel%alam_gammap(3,bm) = r1 ! fitting parameter r1
        vemodel%alam_gammap(4,bm) = r2 ! fitting parameter r2
      end if

      if ( vemodel%alamJ_model == 1 ) then
        vemodel%alamJ(1,bm) = betaJ ! exponent pressure dependence
      end if

    end if

  end subroutine fill_PLLA_data_v


! create vemodel from coefficients without any material data yet

  subroutine create_vemodel_without_data ( coefficients, vemodel )

    type(coefficients_t), intent(in) :: coefficients
    type(vemodel_t), intent(inout) :: vemodel

    call create_viscoelastic_model ( model=coefficients%i(18), &
      vemodel=vemodel, nmodes=coefficients%i(19), &
      alam_model=coefficients%i(80), alamJ_model=coefficients%i(96), &
      alam_gammap_model=coefficients%i(98), gammap_mode=coefficients%i(99) )

  end subroutine create_vemodel_without_data


! fill material data in coefficients from vemodel

  subroutine fill_data_in_coefficients ( vemodel, coefficients )

    type(vemodel_t), intent(in) :: vemodel
    type(coefficients_t), intent(inout) :: coefficients

    integer :: smp, emp, np_mode, nq_mode, nr_mode, npar_mode, npar_tot
    integer :: ns_mode

!   fill material parameters in coefficients

    np_mode = size(vemodel%nonlin,1)
    nq_mode = size(vemodel%alam,1)
    nr_mode = size(vemodel%alamJ,1)
    ns_mode = size(vemodel%alam_gammap,1)
    npar_mode = 2 + np_mode + nq_mode + nr_mode + ns_mode
    npar_tot  = vemodel%nmodes * npar_mode

!   material parameters

    if ( coefficients%i(20) > 0 ) then

!     material parameters in real array r

      call parameters_in_r

    else if ( coefficients%i(20) < 0 ) then

!     material parameters in 2D real array ra2

      call parameters_in_ra2

    else

      write(*,'(/3(a/))') &
        'Error fill_data_in_coefficients:', &
        ' The start index of the viscoelastic material parameters ', &
        ' coefficients%i(20) has not been set.'
      stop

    end if

  contains

!   material parameters in real array r

    subroutine parameters_in_r

      integer :: k

!     start pointer

      smp = coefficients%i(20)   ! start of material parameters

      if ( smp < 501 .or. smp > 1000 ) then
        write(*,'(/a/a,i0/3(a/))') &
          'Error fill_data_in_coefficients:', &
          ' The start index of the viscoelastic material parameters = ', smp, &
          ' This is outside the free range of 501:1000.', &
          ' Move the parameters either to the free range, ', &
          ' or use a 2D array in coefficients%ra2 (use coefficients%i(20)<0).'
        stop
      end if

      emp = smp + npar_tot - 1   ! end of material parameters

!     check size of real array

      call check ( coefficients, 'fill_data_in_coefficients', ncoefr=emp )

      if ( emp > 1000 ) then
        write(*,'(/a/a,i0/3(a/))') &
          'Error fill_data_in_coefficients:', &
          ' The end index of the viscoelastic material parameters = ', smp, &
          ' This is beyond the free range of 501:1000.', &
          ' The parameters must be either fully within the free range, ', &
          ' or use a 2D array in coefficients%ra2 (use coefficients%i(20)<0).'
        stop
      end if

      coefficients%r(smp:emp:npar_mode) = vemodel%modulus
      coefficients%r(smp+1:emp:npar_mode) = vemodel%lambda
      do k = 1, np_mode
        coefficients%r(smp+1+k:emp:npar_mode) = vemodel%nonlin(k,:)
      end do
      do k = 1, nq_mode
        coefficients%r(smp+1+np_mode+k:emp:npar_mode) = vemodel%alam(k,:)
      end do
      do k = 1, nr_mode
        coefficients%r(smp+1+np_mode+nq_mode+k:emp:npar_mode) = &
                                                     vemodel%alamJ(k,:)
      end do
      do k = 1, ns_mode
        coefficients%r(smp+1+np_mode+nq_mode+nr_mode+k:emp:npar_mode) = &
                                                     vemodel%alam_gammap(k,:)
      end do

    end subroutine parameters_in_r


!   material parameters in 2D real array ra2

    subroutine parameters_in_ra2

      integer :: m, ib, ie, sh(2)

!     start pointer

      m = - coefficients%i(20)   ! start of material parameters

      if ( m > size(coefficients%ra2) ) then
        write(*,'(/a/a,i0/a,i0/)') &
          'Error fill_data_in_coefficients:', &
          ' The index for the viscoelastic material parameters = ', m, &
          ' which is larger than the size of coefficients%ra2 = ', &
          size(coefficients%ra2)
        stop
      else
        sh = shape(coefficients%ra2(m)%a)
        if ( any( sh /= [npar_mode,vemodel%nmodes] ) ) then
          write(*,'(/a/a,2(a,i0,a,i0,a/))') &
            'Error fill_data_in_coefficients:', &
            ' The shape of ra2 for the viscoelastic material parameters is', &
            ' (', sh(1), ',', sh(2), '),', &
            ' whereas it must be (', npar_mode, ',', vemodel%nmodes, ').'
          stop
        end if
      end if

      coefficients%ra2(m)%a(1,:) = vemodel%modulus
      coefficients%ra2(m)%a(2,:) = vemodel%lambda
      ib = 3; ie = 2 + np_mode
      coefficients%ra2(m)%a(ib:ie,:) = vemodel%nonlin
      ib = ie+1; ie = ie + nq_mode
      coefficients%ra2(m)%a(ib:ie,:) = vemodel%alam
      ib = ie+1; ie = ie + nr_mode
      coefficients%ra2(m)%a(ib:ie,:) = vemodel%alamJ
      ib = ie+1; ie = ie + ns_mode
      coefficients%ra2(m)%a(ib:ie,:) = vemodel%alam_gammap

    end subroutine parameters_in_ra2

  end subroutine fill_data_in_coefficients

end module material_data_m
