
module stokes_functions_m

  use math_defs_m

  implicit none

  save

  real(dp) :: time_m, rho_m, eta_m

contains

  function func ( nr, xin )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: xin
    real(dp) :: func

    real(dp) :: x, y, t

    x = xin(1)
    y = xin(2)
    t = time_m

    select case ( nr )

    case (1)

!     velocity solution in x-direction

      func = - sin ( pi * x ) * cos ( pi * y ** 2 ) * y * sin ( pi * t )

    case (2)

!     velocity solution in y-direction

      func = cos( pi * x ) * sin ( pi * y ** 2 ) / 2 * sin( pi * t )

    case (3)

!     the pressure solution

      func = - cos ( pi * x ) * sin ( pi * y / 2 ) * sin( pi * t )

    case default

      write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
      stop

    end select

  end function func

  function vfunc ( n, nr, xin )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: xin
    real(dp), dimension(n) :: vfunc

    real(dp) :: x, y, t, rho, eta

    x = xin(1)
    y = xin(2)
    t = time_m
    rho = rho_m
    eta = eta_m

    select case ( nr )

    case (1)

!     rhs solution

      vfunc(1) = 1._dp/(4._dp ) * ( sin( pi * x ) * &
        ( -4._dp * pi * y * rho * cos ( pi * t ) * cos( pi * y**2 ) + &
        sin ( pi * t ) * ( 4._dp * pi * ( -1._dp *pi * ( y + 4._dp * y**3) * &
        eta * cos( pi * y**2) + sin( pi * y / 2._dp ) - 6._dp * y * eta * &
        sin( pi * y**2 )) - rho * cos( pi * x ) * sin( pi * t) * &
        ( -4._dp * pi * y**2 + sin( 2._dp * pi * y**2)))))

      vfunc(2) = 1._dp/(4._dp ) * ( 2._dp * pi * cos ( pi * x ) * ( &
        -1._dp *( cos( pi * y / 2._dp ) + 2._dp * eta * cos( pi * y**2 ) ) * &
        sin ( pi * t ) + ( rho * cos( pi * t ) + pi * ( 1._dp  + 4._dp * y**2)*&
        eta * sin( pi * t )) * sin( pi * y**2) ) + pi * y * rho * &
        sin( pi * t ) **2 * sin( 2._dp * pi * y**2 ) )

    case default

      write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
      stop

    end select

  end function vfunc

end module stokes_functions_m
