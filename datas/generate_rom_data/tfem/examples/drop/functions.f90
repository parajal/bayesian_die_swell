module functions_m

  use math_defs_m

  implicit none

  real(dp), save :: epsilondot = 1._dp, gammadot = 1._dp, omega = 1._dp, &
    time = 1._dp, R0 = 0.1_dp, ts = -1._dp

contains

  function vfunc_2D ( n, nr, x )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc_2D

    select case(nr)
      case(0) ! zero
        vfunc_2D = 0
      case(1) ! planar elongation
        vfunc_2D(1) = x(1)
        vfunc_2D(2) = -x(2)
        vfunc_2D = epsilondot * vfunc_2D
      case(2) ! shear
        vfunc_2D(1) = gammadot * x(2)
        vfunc_2D(2) = 0
      case(3) ! rotation
        vfunc_2D(1) = x(2)
        vfunc_2D(2) = -x(1)
        vfunc_2D = omega * vfunc_2D
      case(4) ! uniaxial elongation
        vfunc_2D(1) = x(1)
        vfunc_2D(2) = -x(2)/2
        vfunc_2D = epsilondot * vfunc_2D
      case default
        write(*,'(/a,i0/)') 'Error vfunc_2D: wrong function number: ', nr
        stop
    end select

    if ( ts > 0 ) then

!     start-up for Navier-Stokes

      if ( time <= ts ) then
        vfunc_2D = vfunc_2D * ( 1 - cos(pi*time/ts) ) / 2
      end if

    end if

  end function vfunc_2D

  function func_error_2D ( nr, x )

    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func_error_2D
    real(dp) :: a, b, lambda, x0, y0

    select case(nr)
      case(0) ! zero
        func_error_2D = sqrt( x(1)**2 + x(2) ** 2 ) - R0
      case(1) ! planar elongation
        a = exp(epsilondot*time)
        b = exp(-epsilondot*time)
        lambda = R0 / sqrt( (x(1)/a)**2 + (x(2)/b)**2 )
        func_error_2D = abs(lambda-1)* sqrt( x(1)**2 + x(2) ** 2 )
      case(2) ! shear
        x0 = x(1) - gammadot * x(2) * time
        y0 = x(2)
        lambda = R0 / sqrt( x0**2 + y0**2 )
        func_error_2D = abs(lambda-1)* sqrt( x(1)**2 + x(2) ** 2 )
      case(3) ! rotation
        func_error_2D = sqrt( x(1)**2 + x(2) ** 2 ) - R0
      case(4) ! uniaxial elongation
        a = exp(epsilondot*time)
        b = exp(-epsilondot*time/2)
        lambda = R0 / sqrt( (x(1)/a)**2 + (x(2)/b)**2 )
        func_error_2D = abs(lambda-1)* sqrt( x(1)**2 + x(2) ** 2 )
      case default
        write(*,'(/a,i0/)') 'Error func_error_2D: wrong function number: ', nr
        stop
    end select

  end function func_error_2D

  function vfunc_3D ( n, nr, x )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc_3D

    select case(nr)
      case(0) ! zero
        vfunc_3D = 0
      case(1) ! planar elongation
        vfunc_3D(1) = x(1)
        vfunc_3D(2) = -x(2)
        vfunc_3D(3) = 0
        vfunc_3D = epsilondot * vfunc_3D
      case(2) ! shear
        vfunc_3D(1) = gammadot * x(2)
        vfunc_3D(2) = 0
        vfunc_3D(3) = 0
      case(3) ! rotation
        vfunc_3D(1) = x(2)
        vfunc_3D(2) = -x(1)
        vfunc_3D(3) = 0
        vfunc_3D = omega * vfunc_3D
      case(4) ! uniaxial elongation
        vfunc_3D(1) = x(1)
        vfunc_3D(2) = -x(2)/2
        vfunc_3D(3) = -x(3)/2
        vfunc_3D = epsilondot * vfunc_3D
      case default
        write(*,'(/a,i0/)') 'Error vfunc_3D: wrong function number: ', nr
        stop
    end select

    if ( ts > 0 ) then

!     start-up for Navier-Stokes

      if ( time <= ts ) then
        vfunc_3D = vfunc_3D * ( 1 - cos(pi*time/ts) ) / 2
      end if

    end if

  end function vfunc_3D

  function func_error_3D ( nr, x )

    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func_error_3D
    real(dp) :: a, b, lambda, x0, y0, z0

    select case(nr)
      case(0) ! zero
        func_error_3D = sqrt( x(1)**2 + x(2) ** 2 + x(3) ** 2 ) - R0
      case(1) ! planar elongation
        a = R0 * exp(epsilondot*time)
        b = R0 * exp(-epsilondot*time)
        lambda = 1 / sqrt( (x(1)/a)**2 + (x(2)/b)**2 + (x(3)/R0)**2 )
        func_error_3D = abs(lambda-1)* sqrt( x(1)**2 + x(2) ** 2 + x(3) ** 2 )
      case(2) ! shear
        x0 = x(1) - gammadot * x(2) * time
        y0 = x(2)
        z0 = x(3)
        lambda = R0 / sqrt( x0**2 + y0**2 + z0**2 )
        func_error_3D = abs(lambda-1)* sqrt( x(1)**2 + x(2) ** 2 + x(3) ** 2 )
      case(3) ! rotation
        func_error_3D = sqrt( x(1)**2 + x(2) ** 2 + x(3) ** 2 ) - R0
      case(4) ! uniaxial elongation
        a = exp(epsilondot*time)
        b = exp(-epsilondot*time/2)
        lambda = R0 / sqrt( (x(1)/a)**2 + (x(2)/b)**2 + (x(3)/b)**2 )
        func_error_3D = abs(lambda-1)* sqrt( x(1)**2 + x(2) ** 2 + x(3) ** 2 )
      case default
        write(*,'(/a,i0/)') 'Error func_error_3D: wrong function number: ', nr
        stop
    end select

  end function func_error_3D

end module functions_m
