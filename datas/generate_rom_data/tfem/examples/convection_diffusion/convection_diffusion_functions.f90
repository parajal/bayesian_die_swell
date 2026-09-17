module convection_diffusion_functions_m

  use math_defs_m

  implicit none

contains

  function alphafunc ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: alphafunc

    select case(nr)
      case(1)
        alphafunc = 1.5_dp + cos(pi*x(1))*cos(pi*x(2))
      case default
        write(*,'(/a,i0/)') 'Error alphafunc: wrong function number: ', nr
        stop
    end select

  end function alphafunc

  function alphatfunc ( n, nr, x )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n,n) :: alphatfunc

    real(dp) :: c(2,2)

    select case(nr)
      case(1)
        c(:,1) = [1._dp,-0.1_dp]
        c(:,2) = [-0.1_dp,2._dp]
        alphatfunc = c * ( 1.5_dp + cos(pi*x(1))*cos(pi*x(2)) )
      case(2)
        c(:,1) = [1._dp,0._dp]
        c(:,2) = [0._dp,1._dp]
        alphatfunc = c * ( 1.5_dp + cos(pi*x(1))*cos(pi*x(2)) )
      case default
        write(*,'(/a,i0/)') 'Error alphatfunc: wrong function number: ', nr
        stop
    end select

  end function alphatfunc

  function gammafunc ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: gammafunc

    write(*,*) 'function gammafunc has not been defined'
    gammafunc = 0
    stop

  end function gammafunc

  function betafunc ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: betafunc

    write(*,*) 'function betafunc has not been defined'
    betafunc = 0
    stop

  end function betafunc

  function velocityvfunc ( n, nr, x )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: velocityvfunc

    select case(nr)
      case(1)
        velocityvfunc = [100._dp,50._dp] * &
                             ( 1.5_dp + cos(pi*x(1))*cos(pi*x(2)))
      case default
        write(*,'(/a,i0/)') 'Error velocityvfunc: wrong function number: ', nr
        stop
    end select

  end function velocityvfunc

end module convection_diffusion_functions_m
