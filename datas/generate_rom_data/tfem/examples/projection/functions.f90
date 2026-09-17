module functions_m

  use math_defs_m

  implicit none

contains

  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    select case(nr)
      case(1)
        func = x(2)
      case(2)
        func = 0
      case(3)
        func = 1 + cos(pi*x(1))*cos(pi*x(2))
      case(4)
        func = 2*pi**2*cos(pi*x(1))*cos(pi*x(2))
      case default
        write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
        stop
    end select

  end function func

!  function vfunc ( n, nr, x )
!
!    integer, intent(in) :: n, nr
!    real(dp), intent(in), dimension(:) :: x
!    real(dp), dimension(n) :: vfunc
!
!    select case(nr)
!      case default
!        write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
!        stop
!    end select
!
!  end function vfunc

end module functions_m
