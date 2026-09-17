module velocity_functions_m

  use tfem_m
  implicit none

contains

  function vel_func ( nr, x )
    use kind_defs_m
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: vel_func

    if ( nr == 1 ) then

      vel_func = 1._dp

    else if ( nr == 2 ) then

      vel_func = 16._dp * x(1)**2 * ( 1._dp - x(1) )**2

    else if ( nr == 3 ) then

      if (x(2) > 0.5) then
        vel_func = 10._dp*(x(2)-1._dp)*(0.5_dp - x(2))
      else
        vel_func = 0._dp
      end if

    else if ( nr == 4 ) then

      vel_func = -10._dp*(x(2)-1._dp)*x(2)

    else if ( nr == 5 ) then

      if (x(2) > 0.5_dp) then
        vel_func = 100._dp*(x(2)-0.5_dp)*(x(2)-1.0_dp) * x(3)*(x(3)-1.0_dp)
      else
        vel_func = 100._dp*(x(2)-0.5_dp)*x(2)          * x(3)*(x(3)-1.0_dp)
      end if

    else

      print*, ' wrong number in velocity function ...'
      stop

    end if

  end function vel_func



end module velocity_functions_m

