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
      case(2)
        alphafunc = 1.5_dp + cos(pi*x(1))*cos(pi*x(2))**cos(pi*x(3))
      case default
        write(*,'(/a,i0/)') 'Error alphafunc: wrong function number: ', nr
        stop
    end select

  end function alphafunc

end module convection_diffusion_functions_m
