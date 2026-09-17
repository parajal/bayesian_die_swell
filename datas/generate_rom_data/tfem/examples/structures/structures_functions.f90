module structures_functions_m

  use math_defs_m

  implicit none

  real(dp) :: Lb = 1._dp, q0b = 1._dp

contains

  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    select case(nr)
      case(1)
        func = 0
      case(2)
        func = 1
      case(3)
        func = 1 + cos(pi*x(1))
      case(4)
        func = - q0b * ( 1 - x(1) / Lb )
      case(5)
        func = - q0b * ( 1 - x(1) / Lb ) ** 2
      case default
        write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
        stop
    end select

  end function func

  function vfunc ( n, nr, x )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    select case(nr)
      case(1)
        vfunc = - [ -pi*sin(pi*x(1))*cos(pi*x(2)) ,  &
                    -pi*cos(pi*x(1))*sin(pi*x(2)) ]
      case default
        write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
        stop
    end select

  end function vfunc

end module structures_functions_m
