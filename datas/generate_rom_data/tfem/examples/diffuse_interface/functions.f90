module functions_m

  use kind_defs_m
  use tfem_elem_m

  implicit none

  real(dp) :: xi
  real(dp) :: r
  real(dp), allocatable :: c(:)

contains

! functions for the initial c-field

  function difunc ( nr, x )

    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: difunc, dist

!   function to create a circular drop
    select case(nr)
      case(1)
        dist = sqrt( sum((x-c)**2) ) - r

        difunc = tanh( dist / ( sqrt(2.0_dp) * xi))

      case default
        write(*,'(/a,i0/)') 'Error difunc: wrong function number: ', nr
        stop
    end select

  end function difunc

end module functions_m
