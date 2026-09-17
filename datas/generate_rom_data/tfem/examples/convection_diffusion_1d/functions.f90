module functions_m

  use math_defs_m

  implicit none

  real(dp) :: const = 1, sigma = 0.05, Lh=0.5_dp

! namelist /sourcepar/:
!
! alpha, sigma, Lh: parameters in the source functions:
!
! nr   0: source=0
!      1: source=const
!      2: source=const/sqrt(pi)/sigma*exp(-(x-Lh)^2/sigma^2)
!      3: source=0 for x<Lh, const for x>=Lh
  namelist /sourcepar/ const, sigma, Lh

contains

  function source ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in) :: x
    real(dp) :: source

    ! source term

    select case(nr)
      case(0)
        source = 0
      case(1)
        source = const
      case(2)
        source = const / sqrt(pi) / sigma * exp(-(x-Lh)**2/sigma**2)
      case(3)
        if ( x < Lh ) then
          source = 0
        else
          source = const
        end if
      case default
        write(*,*) 'invalid function number in source function'
    end select

  end function source

end module functions_m
