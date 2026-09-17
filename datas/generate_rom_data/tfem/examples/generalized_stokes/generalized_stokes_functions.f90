module generalized_stokes_functions_m

  use kind_defs_m

  implicit none

contains

  function etafunc ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: etafunc

    select case(nr)
      case(1)
        etafunc = 1
      case(2)
        if ( x(2) < 0.2_dp .or. x(2) > 0.8_dp ) then
          etafunc = 0.2_dp
        else
          etafunc = 1_dp
        end if
      case default
        write(*,'(/a,i0/)') 'Error etafunc: wrong function number: ', nr
        stop
    end select

  end function etafunc

end module generalized_stokes_functions_m
