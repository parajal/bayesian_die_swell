module functions_m

  use kind_defs_m
  use tfem_elem_m

  implicit none

  real(dp) :: xi
  real(dp), allocatable, dimension(:) :: rdrop
  real(dp), allocatable, dimension(:,:) :: cdrop

contains

! functions for the initial c-field

  function difunc ( nr, x )

    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: difunc

    integer :: i
    real(dp), allocatable :: dist(:)

    allocate(dist(size(cdrop,1)))

!   function to create multiple circular droplets
    select case(nr)
      case(1)

        do i = 1, size(cdrop,1)
          dist(i) = sqrt( sum((x-cdrop(i,:))**2) ) - rdrop(i)
        end do

        difunc = tanh( minval(dist) / ( sqrt(2.0_dp) * xi))

      case default
        write(*,'(/a,i0/)') 'Error difunc: wrong function number: ', nr
        stop
    end select

    deallocate(dist)

  end function difunc

end module functions_m
