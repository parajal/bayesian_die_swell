module functions_m

  use kind_defs_m

  implicit none

  real(dp) :: U0 = 1._dp, w0 = 0.1_dp, L0 = 3._dp
  real(dp) :: h0 = 1._dp

contains

  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    select case(nr)
      case(1)
        if ( x(1) < w0 ) then
          func = U0 * ( 1 - x(1)/w0 )
        else if ( x(1) > L0 - w0 ) then
          func = U0 * ( 1 - (L0 - x(1))/w0 )
        else
          func = 0
        end if
      case(2)
        if ( x(1) < w0 ) then
          func = U0 * ( 1 - x(1)/w0 ) ** 2
        else if ( x(1) > L0 - w0 ) then
          func = U0 * ( 1 - (L0 - x(1))/w0 ) ** 2
        else
          func = 0
        end if
      case(3)
        func = U0 * ( 1 - (x(2)/h0)**2 )
      case default
        write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
        stop
    end select

  end function func

end module functions_m

