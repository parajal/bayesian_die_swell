module functions_m

  use kind_defs_m

  implicit none

  save

  real(dp) :: R_w = 1, U_max = 1

contains

  function uprofile ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: uprofile

    select case (nr)
    case(1)
      uprofile = U_max * ( 1 - ( x(2) / R_w ) ** 2 )
    case default
      write(*,'(/a,i0/)') 'Error uprofile: wrong function number: ', nr
      stop
    end select

  end function uprofile

end module functions_m
