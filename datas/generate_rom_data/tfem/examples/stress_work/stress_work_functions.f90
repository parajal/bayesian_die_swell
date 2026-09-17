module stress_work_functions_m

  use kind_defs_m

  implicit none

  real(dp) :: eta_l = 1.0_dp, kappa_l = 1.0_dp, T0_l = 1.0_dp
  real(dp) :: Uavg_l = 1.0_dp, R_l = 1.0_dp
  real(dp) :: R1_l = 1.0_dp, R2_l = 1.0_dp, Uth1_l = 1.0_dp, Uth2_l = 1.0_dp
  real(dp) :: lambda_l = 1.0_dp

contains

  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    real(dp) :: c1, c2, e1, fac, ratio, gd

    select case(nr)
      case(1)
        func = 0
      case(2)
        func = 1
      case(3)
!       full pipe developed flow (note431)
        func = 2 * Uavg_l * ( 1 - x(2)**2 / R_l**2 )
      case(4)
!       full pipe developed temperature (note431)
        func = eta_l * Uavg_l**2 / kappa_l * ( 1 - x(2)**4/R_l**4 ) + T0_l
      case(5)
!       two concentric pipes developed flow (Dantzig & Tucker p.145; note432).
        ratio = R1_l/R2_l
        fac = 1 + ratio**2 - (1-ratio**2)/log(1/ratio)
        func = 2 * Uavg_l / fac * ( 1 - x(2)**2 / R2_l**2 &
                    + (1-ratio**2)/log(1/ratio) * log(x(2)/R2_l) )
      case(6)
!       two concentric pipes swirl developed flow (note432).
        c1 = ( Uth2_l*R2_l - Uth1_l*R1_l ) / ( R2_l**2 - R1_l**2 )
        c2 = R1_l*R2_l * ( Uth1_l*R2_l - Uth2_l*R1_l ) / ( R2_l**2 - R1_l**2 )
        func = c1 * x(2) + c2 / x(2)
      case(7)
!       two concentric pipes swirl only developed temperature (note432).
        c2 = R1_l*R2_l * ( Uth1_l*R2_l - Uth2_l*R1_l ) / ( R2_l**2 - R1_l**2 )
        e1 = (1/R1_l**2-1/R2_l**2)/log(R2_l/R1_l)
        func = T0_l+eta_l*c2**2/kappa_l*(1/R2_l**2-1/x(2)**2-e1*log(x(2)/R2_l))
      case(8)
!       full pipe developed czz for Oldroyd-B
        gd = - 4 * Uavg_l *x(2) / R_l**2
        func = 1 + 2*(lambda_l*gd)**2
      case(9)
!       full pipe developed czr for Oldroyd-B
        gd = - 4 * Uavg_l *x(2) / R_l**2
        func = lambda_l*gd
      case default
        write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
        stop
    end select

  end function func

end module stress_work_functions_m



