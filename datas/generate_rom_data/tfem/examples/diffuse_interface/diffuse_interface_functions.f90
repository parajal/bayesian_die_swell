module diffuse_interface_functions_m

  use math_defs_m

  implicit none

  real(dp), save :: gammadot=0, H=1, time=0, L=1

contains

  function mapcoor ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(size(x)) :: mapcoor

!   Lees-Edwards mapping
    mapcoor(1) = mod ( x(1) + gammadot * H * time, L )
    mapcoor(2) = x(2)

  end function mapcoor

  function difunc ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: difunc

    real(dp) :: r, r01, r1, r2, x1, x2, y1, y2, C

    C = 0.04_dp

    select case(nr)
      case(1)
!       unstably stratified layers
        y1 = (x(2)+0.15_dp)/(sqrt(2._dp)*C)
        y2 = (x(2)-0.15_dp)/(sqrt(2._dp)*C)
        y1 = y1 - 0.05_dp*cos(pi*x(1))/(sqrt(2._dp)*C)
        y2 = y2 + 0.05_dp*cos(pi*x(1))/(sqrt(2._dp)*C)
        difunc = -(exp(y1)-exp(-y1))/(exp(y1)+exp(-y1))  +  &
                  (exp(y2)-exp(-y2))/(exp(y2)+exp(-y2)) + 1._dp
      case(2)
!       non-circular drop
        r = (sqrt(3._dp*x(1)*x(1) + x(2)*x(2)) - .5_dp ) / (0.02_dp*sqrt (2._dp))
        difunc = (exp(r)-exp(-r))/(exp(r)+exp(-r))
      case(3)
!       spinodal decomposition
        r = 100._dp + 100._dp*(x(1)+10._dp)*(x(2)+50._dp)
        call random_number (r)
        difunc = (r - 0.5_dp) * 1.e-5_dp
      case(4)
!       a supercritical and a subcritical nucleus
        r01 = 0.15_dp
        x1 = 0.3_dp
        y1 = 0.3_dp
        x2 = -0.3_dp
        y2 = -0.3_dp
        r1 = sqrt((x(1)-x1)**2._dp+(x(2)-y1)**2._dp)
        r2 = sqrt((x(1)-x2)**2._dp+(x(2)-y2)**2._dp)
        difunc = -0.7_dp + 0.7_dp*exp(-(r1/r01)**2._dp)  &
                 + 0.35_dp*exp(-(r2/r01)**2._dp)
      case(5)
!       spinodal decomposition in 3D
        r = 1._dp + 1000._dp*(x(1)+1._dp)*(x(2)+10._dp)*(x(3)+100._dp)
        call random_number (r)
        difunc = (r - 0.5_dp) * 1.e-4_dp
      case(6)
!       non-circular axisymmetrix drop
        r = (sqrt(3._dp*x(1)*x(1) + 3._dp*x(2)*x(2)) - 1.25_dp ) / (0.02_dp*sqrt (2._dp))
        difunc = (exp(r)-exp(-r))/(exp(r)+exp(-r))
      case(7)
!      spinodal decomposition in 1D
        r = 1._dp + 1000._dp*(x(1)+1._dp)
        call random_number (r)
        difunc = 0.3_dp*(0.5_dp-r)
      case(8)
!       spinodal decomposition in 1D
        r = 100._dp + 100._dp*(x(1)+10._dp)*50._dp
        call random_number (r)
        difunc = (r - 0.5_dp) * 1.e-5_dp
      case(10)
        difunc = -0.5_dp * x(2)
      case(11)
        difunc = x(1)
      case(12)
        difunc = 0.25_dp + 0.01_dp * cos(3.0_dp * pi * x(1)) + &
                           0.04_dp * cos(5.0_dp * pi * x(1))
      case(13)
        difunc = 0.1_dp
        r01 = 0.25_dp
        x1 = 0.5_dp
        y1 = 0.45_dp
        r1 = sqrt(3.0_dp*(x(1)-x1)**2._dp+(x(2)-y1)**2._dp)
        if (r1 < r01) difunc = 0.30_dp
      case default
        write(*,'(/a,i0/)') 'Error difunc: wrong function number: ', nr
        stop
    end select

  end function difunc

end module diffuse_interface_functions_m
