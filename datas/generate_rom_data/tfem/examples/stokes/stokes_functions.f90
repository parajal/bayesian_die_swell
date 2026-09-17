module stokes_functions_m

  use kind_defs_m

  implicit none

contains

  function wall_velocity ( n, nr, x )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: wall_velocity

    select case ( nr )
      case(1)
        wall_velocity = [ 0.9_dp, 0.0_dp ]
      case default
        write(*,'(/a,i0/)') 'Error wall_velocity: wrong function number: ', nr
        stop
    end select

  end function wall_velocity

end module stokes_functions_m
