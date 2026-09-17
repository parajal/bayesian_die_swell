module convecreac_functions_m

  use math_defs_m

  implicit none

contains

  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    real(dp) :: c, dcdx, dcdy

    c=1+sin(pi*x(1))*sin(pi*x(2))
    dcdx=pi*cos(pi*x(1))*sin(pi*x(2))
    dcdy=pi*sin(pi*x(1))*cos(pi*x(2))

    select case(nr)
      case(1)
        func=dcdx+2*dcdy+c
      case(2)
        func=c
      case default
        write(*,*) 'invalid function number in func'
    end select

  end function func

  function vfunc ( n, nr, x )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    select case(nr)
      case(1)
        vfunc = [ 1, 2 ]
      case default
        write(*,*) 'invalid function number in vfunc'
    end select

  end function vfunc

end module convecreac_functions_m
