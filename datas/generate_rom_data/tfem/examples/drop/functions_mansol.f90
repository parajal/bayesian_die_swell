module functions_mansol_m

! Functions for (the surface divergence of) the interfacial stress tensor,
! and the interfacial position error.
! Used for the method of manufactured solutions.
! 2D functions: uniaxial extensional flow
! 3D functions: simple shear flow

  use functions_m

  implicit none

contains


  function divsIs2D ( x )

    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(2) :: divsIs2D
    real(dp) :: t, eps

    t = time
    eps = epsilondot

    divsIs2D(1) = (exp(3*eps*t)*x(1)*(2*x(1)**2 + (exp(3*eps*t) + &
      exp(6*eps*t))*x(2)**2))/(x(1)**2 + exp(6*eps*t)*x(2)**2)**2
    divsIs2D(2) = (exp(6*eps*t)*x(2)*(2*x(1)**2 + (exp(3*eps*t) + &
      exp(6*eps*t))*x(2)**2))/(x(1)**2 + exp(6*eps*t)*x(2)**2)**2

  end function divsIs2D

  function divsDs2D ( x )

    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(2) :: divsDs2D
    real(dp) :: t, eps

    t = time
    eps = epsilondot

    divsDs2D(1) = (exp(3*eps*t)*eps*x(1)*(-2*x(1)**4 + exp(3*eps*t)*(-1 + &
      9*exp(3*eps*t))*x(1)**2*x(2)**2 + 2*exp(9*eps*t)*(4 + exp(3*eps*t))* &
      x(2)**4))/(x(1)**2 + exp(6*eps*t)*x(2)**2)**3
    divsDs2D(2) = -((exp(6*eps*t)*eps*x(2)*(11*x(1)**4 + exp(3*eps*t)*(7 + &
      3*exp(3*eps*t))*x(1)**2*x(2)**2 + exp(9*eps*t)*(-2 + exp(3*eps*t))* &
      x(2)**4))/(x(1)**2 + exp(6*eps*t)*x(2)**2)**3)

  end function divsDs2D

  function divsdivsuIs2D ( x )

    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(2) :: divsdivsuIs2D
    real(dp) :: t, eps

    t = time
    eps = epsilondot

    divsdivsuIs2D(1) = (exp(3*eps*t)*eps*x(1)*(-4*x(1)**4 + 2*exp(3*eps*t)* &
      (-1 + 3*exp(3*eps*t))*x(1)**2*x(2)**2 + exp(9*eps*t)*(7 + exp(3*eps*t))* &
      x(2)**4))/(2*(x(1)**2 + exp(6*eps*t)*x(2)**2)**3)
    divsdivsuIs2D(2) = (exp(6*eps*t)*eps*x(2)*(-10*x(1)**4 - 8*exp(3*eps*t)* &
      x(1)**2*x(2)**2 + exp(9*eps*t)*(1 + exp(3*eps*t))*x(2)**4))/(2*(x(1)**2 +&
      exp(6*eps*t)*x(2)**2)**3)

  end function divsdivsuIs2D

  function divstauhyd2D ( imodel, x )

    integer, intent(in) :: imodel
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(2) :: divstauhyd2D
    real(dp) :: eps, t

    eps = epsilondot
    t = time

    select case ( imodel )
    case ( 2:3 )
    divstauhyd2D(1) = (exp(3*eps*t)*x(1)*(exp(3*eps*t)*(-1 + exp(3*eps*t))* &
      x(2)**2 + (Log((x(1)**2 + exp(6*eps*t)*x(2)**2)/(exp(2*eps*t)* &
      (x(1)**2 + exp(3*eps*t)*x(2)**2)))*(2*x(1)**2 + (exp(3*eps*t) + &
      exp(6*eps*t))*x(2)**2))/2))/(x(1)**2 + exp(6*eps*t)*x(2)**2)**2
    divstauhyd2D(2) = (exp(3*eps*t)*x(2)*(-((-1 + exp(3*eps*t))*x(1)**2) + &
      (exp(3*eps*t)*Log((x(1)**2 + exp(6*eps*t)*x(2)**2)/(exp(2*eps*t)* &
      (x(1)**2 + exp(3*eps*t)*x(2)**2)))*(2*x(1)**2 + (exp(3*eps*t) + &
      exp(6*eps*t))*x(2)**2))/2))/(x(1)**2 + exp(6*eps*t)*x(2)**2)**2
    case ( 4:5 )
    divstauhyd2D(1) = (exp(3*eps*t)*Log(Sqrt((x(1)**2 + exp(6*eps*t)*x(2)**2)/ &
      (exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)*x(2)**2))))*x(1))/((x(1)**2 + &
      exp(6*eps*t)*x(2)**2)*Sqrt((x(1)**2 + exp(6*eps*t)*x(2)**2)/ &
      (exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)*x(2)**2)))) - (exp(3*eps*t)*x(1)* &
      (-(exp(3*eps*t)*(-1 + exp(3*eps*t))*x(2)**2) - (Log((x(1)**2 + &
      exp(6*eps*t)*x(2)**2)/(exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)*x(2)**2)))* &
      (x(1)**2 - exp(3*eps*t)*(-2 + exp(3*eps*t))*x(2)**2))/2))/((x(1)**2 + &
      exp(6*eps*t)*x(2)**2)**2*Sqrt((x(1)**2 + exp(6*eps*t)*x(2)**2)/ &
      (exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)*x(2)**2))))
    divstauhyd2D(2) = (exp(6*eps*t)*Log(Sqrt((x(1)**2 + exp(6*eps*t)*x(2)**2)/ &
      (exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)*x(2)**2))))*x(2))/((x(1)**2 + &
      exp(6*eps*t)*x(2)**2)*Sqrt((x(1)**2 + exp(6*eps*t)*x(2)**2)/ &
      (exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)*x(2)**2)))) + (exp(3*eps*t)*x(2)* &
      (-2*(-1 + exp(3*eps*t))*x(1)**2 + Log((x(1)**2 + exp(6*eps*t)*x(2)**2)/ &
      (exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)*x(2)**2)))*((-1 + 2*exp(3*eps*t))* &
      x(1)**2 + exp(6*eps*t)*x(2)**2)))/(2*(x(1)**2 + exp(6*eps*t)* &
      x(2)**2)**2*Sqrt((x(1)**2 + exp(6*eps*t)*x(2)**2)/(exp(2*eps*t)* &
      (x(1)**2 + exp(3*eps*t)*x(2)**2))))
    case default
      write(*,'(/a,i0/)') 'Error divstauhyd2D: wrong value imodel: ', imodel
      stop
    end select

  end function divstauhyd2D

  function divstaudev2D ( imodel, x )

    integer, intent(in) :: imodel
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(2) :: divstaudev2D
    real(dp) :: eps, t

    eps = epsilondot
    t = time

    select case ( imodel )
    case ( 2:3 )
    divstaudev2D(1) = -((exp(4*eps*t)*x(1)*Sqrt((x(1)**2 + exp(6*eps*t)* &
      x(2)**2)/(exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)*x(2)**2)))*(2*x(1)**4* &
      (-1 + exp(eps*t)*Sqrt((x(1)**2 + exp(6*eps*t)*x(2)**2)/(exp(2*eps*t)* &
      (x(1)**2 + exp(3*eps*t)*x(2)**2)))) + exp(7*eps*t)*x(2)**4*(-2* &
      exp(5*eps*t) + Sqrt((x(1)**2 + exp(6*eps*t)*x(2)**2)/(exp(2*eps*t)* &
      (x(1)**2 + exp(3*eps*t)*x(2)**2))) + exp(3*eps*t)*Sqrt((x(1)**2 + &
      exp(6*eps*t)*x(2)**2)/(exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)* &
      x(2)**2)))) + exp(4*eps*t)*x(1)**2*x(2)**2*(-4*exp(2*eps*t) + 3* &
      Sqrt((x(1)**2 + exp(6*eps*t)*x(2)**2)/(exp(2*eps*t)*(x(1)**2 + &
      exp(3*eps*t)*x(2)**2))) + exp(3*eps*t)*Sqrt((x(1)**2 + exp(6*eps*t)* &
      x(2)**2)/(exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)*x(2)**2))))))/ &
      (x(1)**2 + exp(6*eps*t)*x(2)**2)**3)
    divstaudev2D(2) = -((exp(4*eps*t)*x(2)*Sqrt((x(1)**2 + exp(6*eps*t)* &
      x(2)**2)/(exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)*x(2)**2)))* &
      (exp(10*eps*t)*x(2)**4*(-2*exp(2*eps*t) + Sqrt((x(1)**2 + exp(6*eps*t)* &
      x(2)**2)/(exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)*x(2)**2))) + &
      exp(3*eps*t)*Sqrt((x(1)**2 + exp(6*eps*t)*x(2)**2)/(exp(2*eps*t)* &
      (x(1)**2 + exp(3*eps*t)*x(2)**2)))) + 2*x(1)**4*(-1 + exp(4*eps*t)* &
      Sqrt((x(1)**2 + exp(6*eps*t)*x(2)**2)/(exp(2*eps*t)*(x(1)**2 + &
      exp(3*eps*t)*x(2)**2)))) + exp(6*eps*t)*x(1)**2*x(2)**2*(-4 + 3* &
      exp(eps*t)*Sqrt((x(1)**2 + exp(6*eps*t)*x(2)**2)/(exp(2*eps*t)* &
      (x(1)**2 + exp(3*eps*t)*x(2)**2))) + exp(4*eps*t)*Sqrt((x(1)**2 + &
      exp(6*eps*t)*x(2)**2)/(exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)* &
      x(2)**2))))))/(x(1)**2 + exp(6*eps*t)*x(2)**2)**3)
    case ( 4:5 )
    divstaudev2D(1) = (3*exp(7*eps*t)*(-1 + exp(3*eps*t))*x(1)*x(2)**2* &
      (x(1)**2 + exp(3*eps*t)*x(2)**2))/(2*(x(1)**2 + exp(6*eps*t)* &
      x(2)**2)**3) + (exp(4*eps*t)*x(1))/(x(1)**2 + exp(6*eps*t)* &
      x(2)**2) - (exp(5*eps*t)*x(1)*(x(1)**2 + exp(3*eps*t)*x(2)**2)* &
      (exp(-(eps*t)) + x(1)**2/(exp(eps*t)*(x(1)**2 + exp(3*eps*t)* &
      x(2)**2)) + (exp(5*eps*t)*x(2)**2)/(x(1)**2 + exp(3*eps*t)* &
      x(2)**2)))/(2*(x(1)**2 + exp(6*eps*t)*x(2)**2)**2)
    divstaudev2D(2) = (exp(4*eps*t)*(-1 + exp(3*eps*t))*x(2)* &
      (-2*x(1)**4 + exp(3*eps*t)*(-2 + exp(3*eps*t))*x(1)**2* &
      x(2)**2 + exp(9*eps*t)*x(2)**4))/(2*(x(1)**2 + exp(6*eps*t)* &
      x(2)**2)**3) - (exp(2*eps*t)*(x(1)**2 + exp(3*eps*t)*x(2)**2)* &
      (-exp(-(eps*t)) + x(1)**2/(exp(eps*t)*(x(1)**2 + exp(3*eps*t)* &
      x(2)**2))))/(x(2)*(x(1)**2 + exp(6*eps*t)*x(2)**2)) - &
      (exp(8*eps*t)*x(2)*(x(1)**2 + exp(3*eps*t)*x(2)**2)* &
      (exp(-(eps*t)) + x(1)**2/(exp(eps*t)*(x(1)**2 + exp(3*eps*t)* &
      x(2)**2)) + (exp(5*eps*t)*x(2)**2)/(x(1)**2 + exp(3*eps*t)* &
      x(2)**2)))/(2*(x(1)**2 + exp(6*eps*t)*x(2)**2)**2)
    case default
      write(*,'(/a,i0/)') 'Error divstaudev2D: wrong value imodel: ', imodel
      stop
    end select

  end function divstaudev2D

  function func_error_mansol_2D ( nr, x )

    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func_error_mansol_2D
    real(dp) :: a, b, lambda

    select case(nr)
      case(4) ! uniaxial elongation
        a = exp(epsilondot*time)
        b = exp(-epsilondot*time/2)
        lambda = R0 / sqrt( (x(1)/a)**2 + (x(2)/b)**2 )
        func_error_mansol_2D = (lambda-1)* sqrt( x(1)**2 + x(2) ** 2 )
      case default
        write(*,'(/a,a,i0/)') 'Error func_error_mansol_2D: wrong function ', &
          'number: ', nr
        stop
    end select

  end function func_error_mansol_2D

  function divsIs3D ( x )

    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(3) :: divsIs3D
    real(dp) :: t, rate

    t = time
    rate = gammadot

    divsIs3D(1) = ((x(1) - rate*t*x(2))*((2 + rate**2*t**2)*x(1)**2 - 2*rate* &
      t*(3 + rate**2*t**2)*x(1)*x(2) + (2 + 4*rate**2*t**2 + rate**4*t**4)* &
      x(2)**2 + (2 + rate**2*t**2)*x(3)**2))/((1 + rate**2*t**2)*x(1)**2 - &
      2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + &
      rate**4*t**4)*x(2)**2 + x(3)**2)**2
    divsIs3D(2) = ((-(rate*t*x(1)) + x(2) + rate**2*t**2*x(2))*((2 + rate**2* &
      t**2)*x(1)**2 - 2*rate*t*(3 + rate**2*t**2)*x(1)*x(2) + (2 + &
      4*rate**2*t**2 + rate**4*t**4)*x(2)**2 + (2 + rate**2*t**2)*x(3)**2))/ &
      ((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + &
      (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**2
    divsIs3D(3) = (x(3)*((2 + rate**2*t**2)*x(1)**2 - 2*rate*t*(3 + rate**2* &
      t**2)*x(1)*x(2) + (2 + 4*rate**2*t**2 + rate**4*t**4)*x(2)**2 + (2 + &
      rate**2*t**2)*x(3)**2))/((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + &
      rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)* &
      x(2)**2 + x(3)**2)**2

  end function divsIs3D

  function divsDs3D ( x )

    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(3) :: divsDs3D
    real(dp) :: t, rate

    t = time
    rate = gammadot

    divsDs3D(1) = -((rate*(5*rate**9*t**9*x(1)*x(2)**4 - rate**10*t**10* &
      x(2)**5 - rate**8*t**8*x(2)**3*(10*x(1)**2 + 7*x(2)**2 + x(3)**2) + &
      rate**7*t**7*x(1)*x(2)**2*(10*x(1)**2 + 30*x(2)**2 + 3*x(3)**2) - &
      rate**2*t**2*x(2)*(-19*x(1)**4 + 9*x(2)**4 + 16*x(2)**2*x(3)**2 + &
      7*x(3)**4 - 8*x(1)**2*(x(2)**2 - 2*x(3)**2)) + rate*t*x(1)*(-5* &
      x(1)**4 + 7*x(2)**4 + 13*x(2)**2*x(3)**2 + 6*x(3)**4 + x(1)**2*(-18* &
      x(2)**2 + x(3)**2)) + rate**3*t**3*x(1)*(2*x(1)**4 + 18*x(2)**4 + 32* &
      x(2)**2*x(3)**2 + x(3)**4 + 2*x(1)**2*(-7*x(2)**2 + x(3)**2)) + x(2)* &
      (5*x(1)**4 + 2*x(1)**2*(x(2)**2 + x(3)**2) - 3*(x(2)**2 + &
      x(3)**2)**2) + rate**5*t**5*x(1)*(x(1)**4 + 35*x(2)**4 + 12*x(2)**2* &
      x(3)**2 + x(1)**2*(40*x(2)**2 + x(3)**2)) - rate**6*t**6*x(2)*(5* &
      x(1)**4 + 13*x(2)**4 + 5*x(2)**2*x(3)**2 + x(1)**2*(50*x(2)**2 + &
      3*x(3)**2)) - rate**4*t**4*x(2)*(15*x(1)**4 + 13*x(2)**4 + 17* &
      x(2)**2*x(3)**2 + x(3)**4 + x(1)**2*(22*x(2)**2 + 9*x(3)**2))))/ &
      ((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + &
      (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**3)
    divsDs3D(2) = -((rate*((-3 + 4*rate**2*t**2 + rate**4*t**4)*x(1)**5 + &
      rate*t*(7 - 24*rate**2*t**2 - 5*rate**4*t**4)*x(1)**4*x(2) + &
      x(1)**3*(2*(1 + 6*rate**2*t**2 + 28*rate**4*t**4 + 5*rate**6*t**6)* &
      x(2)**2 + (-6 - rate**2*t**2 + rate**4*t**4)*x(3)**2) - rate*t* &
      x(1)**2*x(2)*(2*(9 + 24*rate**2*t**2 + 32*rate**4*t**4 + 5*rate**6* &
      t**6)*x(2)**2 + (-13 - 4*rate**2*t**2 + 3*rate**4*t**4)*x(3)**2) + &
      x(1)*((5 + 30*rate**2*t**2 + 47*rate**4*t**4 + 36*rate**6*t**6 + 5* &
      rate**8*t**8)*x(2)**4 + (2 - 13*rate**2*t**2 - 5*rate**4*t**4 + 3* &
      rate**6*t**6)*x(2)**2*x(3)**2 - (3 + 5*rate**2*t**2)*x(3)**4) - &
      rate*t*x(2)*((5 + 14*rate**2*t**2 + 15*rate**4*t**4 + 8*rate**6*t**6 + &
      rate**8*t**8)*x(2)**4 + (-1 - 6*rate**2*t**2 - 2*rate**4*t**4 + &
      rate**6*t**6)*x(2)**2*x(3)**2 - (6 + 5*rate**2*t**2)*x(3)**4)))/ &
      ((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + &
      (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**3)
    divsDs3D(3) = (rate*x(3)*(-15*rate**6*t**6*x(1)*x(2)**3 + &
      4*rate**7*t**7*x(2)**4 - 8*x(1)*x(2)*(x(1)**2 + x(2)**2 + x(3)**2) + &
      rate**5*t**5*x(2)**2*(21*x(1)**2 + 21*x(2)**2 + 4*x(3)**2) - &
      rate**4*t**4*x(1)*x(2)*(13*x(1)**2 + 69*x(2)**2 + 8*x(3)**2) - &
      rate**2*t**2*x(1)*x(2)*(45*x(1)**2 + 55*x(2)**2 + 17*x(3)**2) + &
      rate**3*t**3*(3*x(1)**4 + 4*x(1)**2*(21*x(2)**2 + x(3)**2) + &
      5*x(2)**2*(5*x(2)**2 + 2*x(3)**2)) + rate*t*(9*x(1)**4 + 9*x(2)**4 + &
      7*x(2)**2*x(3)**2 - 2*x(3)**4 + x(1)**2*(38*x(2)**2 + 7*x(3)**2))))/ &
      ((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + &
      (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**3

  end function divsDs3D

  function divsdivsuIs3D ( x )

    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(3) :: divsdivsuIs3D
    real(dp) :: t, rate

    t = time
    rate = gammadot

    divsdivsuIs3D(1) = -((rate*(rate**6*t**6*x(2)**3*(2*x(2)**2 + x(3)**2) - &
      3*rate**5*t**5*x(1)*x(2)**2*(4*x(2)**2 + x(3)**2) - x(2)*(-3*x(1)**4 - &
      2*x(1)**2*(x(2)**2 + x(3)**2) + (x(2)**2 + x(3)**2)**2) + rate**4*t**4* &
      x(2)*(2*x(2)**4 - x(2)**2*x(3)**2 + 3*x(1)**2*(9*x(2)**2 + x(3)**2)) - &
      rate*t*x(1)*(3*x(1)**4 + x(2)**4 - x(2)**2*x(3)**2 - 2*x(3)**4 + x(1)**2*&
      (14*x(2)**2 + x(3)**2)) + rate**2*t**2*x(2)*(15*x(1)**4 - x(2)**4 - &
      3*x(2)**2*x(3)**2 - 2*x(3)**4 + x(1)**2*(21*x(2)**2 + x(3)**2)) - &
      rate**3*t**3*x(1)*(12*x(2)**4 - x(2)**2*x(3)**2 + x(1)**2*(29*x(2)**2 + &
      x(3)**2))))/((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)* &
      x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**3)
    divsdivsuIs3D(2) = (rate*(-((-1 + 3*rate**2*t**2 + rate**4*t**4)*x(1)**5) +&
      rate*t*(1 + 19*rate**2*t**2 + 5*rate**4*t**4)*x(1)**4*x(2) - x(1)**3* &
      (2*(1 + 11*rate**2*t**2 + 23*rate**4*t**4 + 5*rate**6*t**6)*x(2)**2 + &
      (-2 + rate**2*t**2 + rate**4*t**4)*x(3)**2) + rate*t*x(1)**2*x(2)*(2* &
      (7 + 25*rate**2*t**2 + 27*rate**4*t**4 + 5*rate**6*t**6)*x(2)**2 + (-1 + &
      4*rate**2*t**2 + 3*rate**4*t**4)*x(3)**2) + rate*t*(1 + rate**2*t**2)* &
      x(2)*((3 + 7*rate**2*t**2 + 6*rate**4*t**4 + rate**6*t**6)*x(2)**4 + (1 +&
      rate**2*t**2 + rate**4*t**4)*x(2)**2*x(3)**2 - 2*x(3)**4) - x(1)*((3 + &
      22*rate**2*t**2 + 43*rate**4*t**4 + 31*rate**6*t**6 + 5*rate**8*t**8)* &
      x(2)**4 + (2 + 3*rate**2*t**2 + 5*rate**4*t**4 + 3*rate**6*t**6)*x(2)**2*&
      x(3)**2 - (1 + 2*rate**2*t**2)*x(3)**4)))/((1 + rate**2*t**2)*x(1)**2 - &
      2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4* &
      t**4)*x(2)**2 + x(3)**2)**3
    divsdivsuIs3D(3) = (rate*x(3)*(-11*rate**6*t**6*x(1)*x(2)**3 + 3*rate**7* &
      t**7*x(2)**4 - 4*x(1)*x(2)*(x(1)**2 + x(2)**2 + x(3)**2) + rate**5*t**5* &
      x(2)**2*(15*x(1)**2 + 14*x(2)**2 + 3*x(3)**2) - rate**4*t**4*x(1)*x(2)* &
      (9*x(1)**2 + 44*x(2)**2 + 6*x(3)**2) - rate**2*t**2*x(1)*x(2)*(26* &
      x(1)**2 + 31*x(2)**2 + 14*x(3)**2) + 5*rate*t*(x(1)**4 + x(2)**2* &
      (x(2)**2 + x(3)**2) + x(1)**2*(4*x(2)**2 + x(3)**2)) + rate**3*t**3* &
      (2*x(1)**4 + 15*x(2)**4 + 9*x(2)**2*x(3)**2 + 3*x(1)**2*(17*x(2)**2 + &
      x(3)**2))))/((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)* &
      x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**3

  end function divsdivsuIs3D

  function divstauhyd3D ( imodel, x )

    integer, intent(in) :: imodel
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(3) :: divstauhyd3D
    real(dp) :: rate, t

    rate = gammadot
    t = time

    select case ( imodel )
    case ( 2:3 )
    divstauhyd3D(1) = (rate*t*(-(rate*t*x(1)) + x(2) + rate**2*t**2*x(2))* &
      (-x(1)**2 + rate*t*x(1)*x(2) + x(2)**2 + x(3)**2) + Log(Sqrt(((x(1) - &
      rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)/ &
      (x(2)**2 + (x(1) - rate*t*x(2))**2 + x(3)**2)))*(x(1) - rate*t*x(2))* &
      ((2 + rate**2*t**2)*x(1)**2 - 2*rate*t*(3 + rate**2*t**2)*x(1)*x(2) + &
      (2 + 4*rate**2*t**2 + rate**4*t**4)*x(2)**2 + (2 + rate**2*t**2)* &
      x(3)**2))/((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)* &
      x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**2
    divstauhyd3D(2) = -((-(Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 +&
      x(3)**2)))*(-(rate*t*x(1)) + x(2) + rate**2*t**2*x(2))*((2 + rate**2* &
      t**2)*x(1)**2 - 2*rate*t*(3 + rate**2*t**2)*x(1)*x(2) + (2 + 4*rate**2* &
      t**2 + rate**4*t**4)*x(2)**2 + (2 + rate**2*t**2)*x(3)**2)) + rate*t* &
      (-x(1)**3 + 2*rate*t*x(1)**2*x(2) - x(1)*((-1 + rate**2*t**2)*x(2)**2 + &
      (1 + rate**2*t**2)*x(3)**2) + rate*t*x(2)*(-x(2)**2 + (2 + rate**2*t**2)* &
      x(3)**2)))/((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*&
      x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**2)
    divstauhyd3D(3) = (x(3)*(rate*t*(-2*x(1)*x(2) - 6*rate**2*t**2*x(1)*x(2) - &
      2*rate**4*t**4*x(1)*x(2) + rate**5*t**5*x(2)**2 + rate*t*(2*x(1)**2 + 3* &
      x(2)**2) + rate**3*t**3*(x(1)**2 + 4*x(2)**2)) + Log(Sqrt(((x(1) - rate*t*&
      x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + &
      (x(1) - rate*t*x(2))**2 + x(3)**2)))*((2 + rate**2*t**2)*x(1)**2 - 2*rate*&
      t*(3 + rate**2*t**2)*x(1)*x(2) + (2 + 4*rate**2*t**2 + rate**4*t**4)* &
      x(2)**2 + (2 + rate**2*t**2)*x(3)**2)))/((1 + rate**2*t**2)*x(1)**2 - 2* &
      rate*t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*&
      x(2)**2 + x(3)**2)**2
    case ( 4:5 )
    divstauhyd3D(1) = -((-(rate*t*(-(rate*t*x(1)) + x(2) + &
      rate**2*t**2*x(2))*(-x(1)**2 + rate*t*x(1)*x(2) + x(2)**2 + x(3)**2)) + &
      Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2)))*(-2*x(1)**3 + rate*t*(7 + rate**2*t**2)*x(1)**2*x(2) - &
      2*x(1)*((1 + 5*rate**2*t**2 + rate**4*t**4)*x(2)**2 + (1 + &
      rate**2*t**2)*x(3)**2) + rate*t*x(2)*((3 + 5*rate**2*t**2 + &
      rate**4*t**4)*x(2)**2 + (3 + 2*rate**2*t**2)*x(3)**2)))/ &
      (((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + &
      (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**2* &
      Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2))))
    divstauhyd3D(2) = -((rate*t*(-x(1)**3 + 2*rate*t*x(1)**2*x(2) - &
      x(1)*((-1 + rate**2*t**2)*x(2)**2 + (1 + rate**2*t**2)*x(3)**2) + &
      rate*t*x(2)*(-x(2)**2 + (2 + rate**2*t**2)*x(3)**2)) - &
      Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2)))*(-3*rate**5*t**5*x(1)*x(2)**2 + rate**6*t**6*x(2)**3 + &
      2*x(2)*(x(1)**2 + x(2)**2 + x(3)**2) + rate**4*t**4*x(2)*(3*x(1)**2 + &
      5*x(2)**2 + 2*x(3)**2) - rate**3*t**3*x(1)*(x(1)**2 + 13*x(2)**2 + &
      2*x(3)**2) - rate*t*x(1)*(3*x(1)**2 + 7*x(2)**2 + 3*x(3)**2) + &
      rate**2*t**2*x(2)*(11*x(1)**2 + 5*(x(2)**2 + x(3)**2))))/ &
      (((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + &
      (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**2* &
      Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2))))
    divstauhyd3D(3) = (x(3)*(rate*t*(-2*x(1)*x(2) - 6*rate**2*t**2*x(1)* &
      x(2) - 2*rate**4*t**4*x(1)*x(2) + rate**5*t**5*x(2)**2 + &
      rate*t*(2*x(1)**2 + 3*x(2)**2) + rate**3*t**3*(x(1)**2 + 4*x(2)**2)) + &
      Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2)))*(-((-2 + rate**2*t**2 + rate**4*t**4)*x(1)**2) + &
      2*rate*t*(-2 + 2*rate**2*t**2 + rate**4*t**4)*x(1)*x(2) + &
      (2 + rate**2*t**2 - 3*rate**4*t**4 - rate**6*t**6)*x(2)**2 + (2 + &
      rate**2*t**2)*x(3)**2)))/(((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + &
      rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + &
      x(3)**2)**2*Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2)))
    case default
      write(*,'(/a,i0/)') 'Error divstauhyd3D: wrong value imodel: ', imodel
      stop
    end select

  end function divstauhyd3D

  function divstaudev3D ( imodel, x )

    integer, intent(in) :: imodel
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(3) :: divstaudev3D
    real(dp) :: rate, t

    rate = gammadot
    t = time

    select case ( imodel )
    case ( 2:3 )
    divstaudev3D(1) = -((x(1)**3 - 3*rate*t*x(1)**2*x(2) + x(1)*x(2)**2 + 3* &
      rate**2*t**2*x(1)*x(2)**2 - rate*t*x(2)**3 - rate**3*t**3*x(2)**3 + x(1)*&
      x(3)**2 + rate**2*t**2*x(1)*x(3)**2 - rate*t*x(2)*x(3)**2 - rate**3*t**3*&
      x(2)*x(3)**2 + (x(1) - rate*t*x(2))*((1 + rate**2*t**2)*x(1)**2 - 2*rate*&
      t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)* &
      x(2)**2 + x(3)**2) - 2*x(1)*(x(1)**2 - 2*rate*t*x(1)*x(2) + (1 + rate**2*&
      t**2)*x(2)**2 + x(3)**2)*(((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + &
      rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + &
      x(3)**2)/(x(1)**2 - 2*rate*t*x(1)*x(2) + (1 + rate**2*t**2)*x(2)**2 + &
      x(3)**2))**1.5_dp)/((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2* &
      t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + &
      x(3)**2)**2)
    divstaudev3D(2) = -((((2*x(2) + rate*t*(1 + rate**2*t**2)*(-((2 + rate**2* &
      t**2)*x(1)) + rate*t*(3 + rate**2*t**2)*x(2)))*((1 + rate**2*t**2)* &
      x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + &
      rate**4*t**4)*x(2)**2 + x(3)**2))/(1 + rate**2*t**2)**2 - 2*x(2)* &
      (x(1)**2 - 2*rate*t*x(1)*x(2) + (1 + rate**2*t**2)*x(2)**2 + x(3)**2)* &
      (((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + &
      (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)/(x(1)**2 - 2*rate*&
      t*x(1)*x(2) + (1 + rate**2*t**2)*x(2)**2 + x(3)**2))**1.5_dp + (rate*t* &
      (x(2)**2*(-((-2 + 3*rate**2*t**2 + rate**4*t**4)*x(1)) + rate*t*(1 + 4* &
      rate**2*t**2 + rate**4*t**4)*x(2)) + rate*t*(-(rate*t*(2 + 3*rate**2* &
      t**2 + rate**4*t**4)*x(1)) + (4 + 5*rate**2*t**2 + 4*rate**4*t**4 + &
      rate**6*t**6)*x(2))*x(3)**2))/(1 + rate**2*t**2)**2)/((1 + rate**2*t**2)* &
      x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + &
      rate**4*t**4)*x(2)**2 + x(3)**2)**2)
    divstaudev3D(3) = -((x(3)*(-2*rate**2*t**2*x(1)**2 - rate**4*t**4*x(1)**2 + &
      2*rate*t*x(1)*x(2) + 6*rate**3*t**3*x(1)*x(2) + 2*rate**5*t**5*x(1)*x(2) -&
      3*rate**2*t**2*x(2)**2 - 4*rate**4*t**4*x(2)**2 - rate**6*t**6*x(2)**2 + &
      (2 + rate**2*t**2)*((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2* &
      t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2) -&
      2*(x(1)**2 - 2*rate*t*x(1)*x(2) + (1 + rate**2*t**2)*x(2)**2 + x(3)**2)* &
      (((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + &
      (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)/(x(1)**2 - 2*rate*&
      t*x(1)*x(2) + (1 + rate**2*t**2)*x(2)**2 + x(3)**2))**1.5_dp))/((1 + &
      rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3* &
      rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**2)
    case ( 4:5 )
    divstaudev3D(1) = (rate*t*(-3*rate**7*t**7*x(1)*x(2)**4 + &
      rate**8*t**8*x(2)**5 + 2*rate**5*t**5*x(1)*x(2)**2*(x(1)**2 - &
      6*x(2)**2 - 5*x(3)**2) + 2*rate**6*t**6*x(2)**3*(x(1)**2 + 3*x(2)**2 + &
      2*x(3)**2) + rate**4*t**4*x(2)*(-3*x(1)**4 + 14*x(2)**4 + &
      20*x(2)**2*x(3)**2 + 3*x(3)**4 - 4*x(1)**2*(x(2)**2 - 2*x(3)**2)) + &
      rate**3*t**3*x(1)*(x(1)**4 - 24*x(2)**4 - 40*x(2)**2*x(3)**2 - &
      3*x(3)**4 + x(1)**2*(24*x(2)**2 - 2*x(3)**2)) + 4*rate*t*x(1)*(x(1)**4 - &
      3*x(2)**4 - 5*x(2)**2*x(3)**2 - 2*x(3)**4 + x(1)**2*(3*x(2)**2 - &
      x(3)**2)) + 4*x(2)*(-x(1)**4 + (x(2)**2 + x(3)**2)**2) + &
      2*rate**2*t**2*x(2)*(-9*x(1)**4 + 6*x(2)**4 + 11*x(2)**2*x(3)**2 + &
      5*x(3)**4 + x(1)**2*(x(2)**2 + 12*x(3)**2))))/(2*((1 + rate**2*t**2)* &
      x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + &
      rate**4*t**4)*x(2)**2 + x(3)**2)**3)
    divstaudev3D(2) = -(rate*t*(-((4 + 2*rate**2*t**2 + rate**4*t**4)* &
      x(1)**5) + rate*t*(20 + 13*rate**2*t**2 + 5*rate**4*t**4)*x(1)**4* &
      x(2) + 2*rate*t*x(1)**2*x(2)*((-2 + 23*rate**2*t**2 + 19*rate**4*t**4 + &
      5*rate**6*t**6)*x(2)**2 + (16 + 23*rate**2*t**2 + 6*rate**4*t**4)* &
      x(3)**2) - x(1)*((-4 - 8*rate**2*t**2 + 26*rate**4*t**4 + &
      22*rate**6*t**6 + 5*rate**8*t**8)*x(2)**4 + 2*rate**2*t**2*(27 + &
      28*rate**2*t**2 + 6*rate**4*t**4)*x(2)**2*x(3)**2 + (4 + &
      10*rate**2*t**2 + 3*rate**4*t**4)*x(3)**4) + rate*t*x(2)*((-4 - &
      4*rate**2*t**2 + 6*rate**4*t**4 + 5*rate**6*t**6 + rate**8*t**8)* &
      x(2)**4 + 2*(4 + 15*rate**2*t**2 + 11*rate**4*t**4 + 2*rate**6*t**6)* &
      x(2)**2*x(3)**2 + (12 + 13*rate**2*t**2 + 3*rate**4*t**4)*x(3)**4) - &
      2*x(1)**3*(5*rate**6*t**6*x(2)**2 + 4*x(3)**2 + 2*rate**4*t**4*(8* &
      x(2)**2 + x(3)**2) + 3*rate**2*t**2*(7*x(2)**2 + 2*x(3)**2))))/(2*((1 + &
      rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3* &
      rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**3)
    divstaudev3D(3) = (rate*t*x(3)*(-8*rate**8*t**8*x(1)*x(2)**3 + 2* &
      rate**9*t**9*x(2)**4 - 8*x(1)*x(2)*(x(1)**2 + x(2)**2 + x(3)**2) - &
      4*rate**6*t**6*x(1)*x(2)*(2*x(1)**2 + 13*x(2)**2 + x(3)**2) + rate**7* &
      t**7*x(2)**2*(12*x(1)**2 + 15*x(2)**2 + 2*x(3)**2) - 4*rate**2*t**2* &
      x(1)*x(2)*(12*x(1)**2 + 17*x(2)**2 + 4*x(3)**2) - 4*rate**4*t**4*x(1)* &
      x(2)*(9*x(1)**2 + 27*x(2)**2 + 4*x(3)**2) + 4*rate*t*(2*x(1)**4 + 3* &
      x(2)**4 + 2*x(2)**2*x(3)**2 - x(3)**4 + x(1)**2*(10*x(2)**2 + &
      x(3)**2)) + 2*rate**5*t**5*(x(1)**4 + 19*x(2)**4 + 5*x(2)**2*x(3)**2 + &
      x(1)**2*(33*x(2)**2 + x(3)**2)) + rate**3*t**3*(7*x(1)**4 + 36*x(2)**4 + &
      14*x(2)**2*x(3)**2 - x(3)**4 + 2*x(1)**2*(55*x(2)**2 + 3*x(3)**2))))/ &
      (2*((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + &
      (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**3)
    case default
      write(*,'(/a,i0/)') 'Error divstaudev3D: wrong value imodel: ', imodel
      stop
    end select

  end function divstaudev3D

  function Is3D ( x )

    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(6) :: Is3D
    real(dp) :: rate, t

    rate = gammadot
    t = time

    Is3D(1) = 1 - (x(1) - rate*t*x(2))**2/((x(1) - rate*t*x(2))**2 + (x(2) + &
      rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)
    Is3D(2) = -(((x(1) - rate*t*x(2))*(-(rate*t*x(1)) + x(2) + rate**2*t**2* &
      x(2)))/((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t* &
      x(2)))**2 + x(3)**2))
    Is3D(3) = -(((x(1) - rate*t*x(2))*x(3))/((x(1) - rate*t*x(2))**2 + &
      (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2))
    Is3D(4) = 1 - (-(rate*t*x(1)) + x(2) + rate**2*t**2*x(2))**2/((x(1) - &
      rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)
    Is3D(5) = -(((-(rate*t*x(1)) + x(2) + rate**2*t**2*x(2))*x(3))/((x(1) - &
      rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2))
    Is3D(6) = 1 - x(3)**2/((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)

  end function Is3D

  function divsuIs3D ( x )

    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(6) :: divsuIs3D
    real(dp) :: rate, t

    rate = gammadot
    t = time

    divsuIs3D(1) = -((rate*(x(1) - rate*t*x(2))*(-(rate*t*x(1)) + x(2) + &
      rate**2*t**2*x(2))*(1 - (x(1) - rate*t*x(2))**2/((x(1) - rate*t* &
      x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)))/ &
      ((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + &
      x(3)**2))
    divsuIs3D(2) = (rate*(x(1) - rate*t*x(2))**2*(-(rate*t*x(1)) + x(2) + &
      rate**2*t**2*x(2))**2)/((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2)**2
    divsuIs3D(3) = (rate*(x(1) - rate*t*x(2))**2*(-(rate*t*x(1)) + x(2) + &
      rate**2*t**2*x(2))*x(3))/((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2)**2
    divsuIs3D(4) = (rate*(-x(1) + rate*t*x(2))*(-(rate*t*x(1)) + x(2) + &
      rate**2*t**2*x(2))*(x(1)**2 - 2*rate*t*x(1)*x(2) + rate**2*t**2*x(2)**2 +&
      x(3)**2))/((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*&
      x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**2
    divsuIs3D(5) = (rate*(x(1) - rate*t*x(2))*(-(rate*t*x(1)) + x(2) + &
      rate**2*t**2*x(2))**2*x(3))/((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2)**2
    divsuIs3D(6) = -((rate*(x(1) - rate*t*x(2))*(-(rate*t*x(1)) + x(2) + &
      rate**2*t**2*x(2))*(1 - x(3)**2/((x(1) - rate*t*x(2))**2 + (x(2) + rate* &
      t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)))/((x(1) - rate*t*x(2))**2 + &
      (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2))

  end function divsuIs3D

  function Ds3D ( x )

    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(6) :: Ds3D
    real(dp) :: rate, t

    rate = gammadot
    t = time

    Ds3D(1) = (-2*rate*(x(1) - rate*t*x(2))*(-(rate*t*x(1)) + x(2) + rate**2* &
      t**2*x(2))*(1 - (x(1) - rate*t*x(2))**2/((x(1) - rate*t*x(2))**2 + &
      (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)))/((x(1) - rate*t* &
      x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)
    Ds3D(2) = (rate*(x(1) - rate*t*x(2))**2*(-(rate*t*x(1)) + x(2) + rate**2* &
      t**2*x(2))**2 + rate*(x(1)**2 - 2*rate*t*x(1)*x(2) + rate**2*t**2* &
      x(2)**2 + x(3)**2)*((x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + &
      x(3)**2))/((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t* &
      x(2)))**2 + x(3)**2)**2
    Ds3D(3) = -((rate*(-(rate*t*x(1)) + x(2) + rate**2*t**2*x(2))*x(3)*((-1 + &
      rate**2*t**2)*x(1)**2 - 2*rate**3*t**3*x(1)*x(2) + (1 + rate**2*t**2 + &
      rate**4*t**4)*x(2)**2 + x(3)**2))/((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*&
      (2 + rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)* &
      x(2)**2 + x(3)**2)**2)
    Ds3D(4) = (2*rate*(-x(1) + rate*t*x(2))*(-(rate*t*x(1)) + x(2) + rate**2* &
      t**2*x(2))*(x(1)**2 - 2*rate*t*x(1)*x(2) + rate**2*t**2*x(2)**2 + &
      x(3)**2))/((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)* &
      x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**2
    Ds3D(5) = -((rate*(-x(1) + rate*t*x(2))*x(3)*((-1 + rate**2*t**2)* &
      x(1)**2 - 2*rate**3*t**3*x(1)*x(2) + (1 + rate**2*t**2 + rate**4*t**4)* &
      x(2)**2 - x(3)**2))/((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*&
      t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + &
      x(3)**2)**2)
    Ds3D(6) = (2*rate*(x(1) - rate*t*x(2))*(-(rate*t*x(1)) + x(2) + rate**2* &
      t**2*x(2))*x(3)**2)/((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)**2

  end function Ds3D

  function tauhyd3D ( imodel, x )

    integer, intent(in) :: imodel
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(6) :: tauhyd3D
    real(dp) :: rate, t

    rate = gammadot
    t = time

    select case ( imodel )
    case ( 2:3 )
    tauhyd3D(1) = Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2)))*(1 - (x(1) - rate*t*x(2))**2/((x(1) - rate*t*x(2))**2 + &
      (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2))
    tauhyd3D(2) = -((Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 +&
      x(3)**2)))*(x(1) - rate*t*x(2))*(-(rate*t*x(1)) + x(2) + rate**2*t**2* &
      x(2)))/((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t* &
      x(2)))**2 + x(3)**2))
    tauhyd3D(3) = -((Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 +&
      x(3)**2)))*(x(1) - rate*t*x(2))*x(3))/((x(1) - rate*t*x(2))**2 + (x(2) + &
      rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2))
    tauhyd3D(4) = Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2)))*(1 - (-(rate*t*x(1)) + x(2) + rate**2*t**2*x(2))**2/((x(1) - &
      rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2))
    tauhyd3D(5) = -((Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 +&
      x(3)**2)))*(-(rate*t*x(1)) + x(2) + rate**2*t**2*x(2))*x(3))/((x(1) - &
      rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2))
    tauhyd3D(6) = Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2)))*(1 - x(3)**2/((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2))
    case ( 4:5 )
    tauhyd3D(1) = (Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2)))*(1 - (x(1) - rate*t*x(2))**2/((x(1) - rate*t*x(2))**2 + &
      (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)))/Sqrt(((x(1) - rate*&
      t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)/ &
      (x(2)**2 + (x(1) - rate*t*x(2))**2 + x(3)**2))
    tauhyd3D(2) = -((Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 +&
      x(3)**2)))*(x(1) - rate*t*x(2))*(-(rate*t*x(1)) + x(2) + rate**2*t**2* &
      x(2)))/((x(2)**2 + (x(1) - rate*t*x(2))**2 + x(3)**2)*(((x(1) - rate*t* &
      x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 +&
     (x(1) - rate*t*x(2))**2 + x(3)**2))**1.5_dp))
    tauhyd3D(3) = -((Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 +&
      x(3)**2)))*(x(1) - rate*t*x(2))*x(3))/((x(2)**2 + (x(1) - rate*t* &
      x(2))**2 + x(3)**2)*(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2))**1.5_dp))
    tauhyd3D(4) = (Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2)))*(1 - (-(rate*t*x(1)) + x(2) + rate**2*t**2*x(2))**2/((x(1) - &
      rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)))/ &
      Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t* &
      x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + x(3)**2))
    tauhyd3D(5) = -((Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 +&
      x(3)**2)))*(-(rate*t*x(1)) + x(2) + rate**2*t**2*x(2))*x(3))/((x(2)**2 + &
      (x(1) - rate*t*x(2))**2 + x(3)**2)*(((x(1) - rate*t*x(2))**2 + (x(2) + &
      rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t* &
      x(2))**2 + x(3)**2))**1.5_dp))
    tauhyd3D(6) = (Log(Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + &
      rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + &
      x(3)**2)))*(1 - x(3)**2/((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2)))/Sqrt(((x(1) - rate*t*x(2))**2 + &
      (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - &
      rate*t*x(2))**2 + x(3)**2))
    case default
      write(*,'(/a,i0/)') 'Error tauhyd3D: wrong value imodel: ', imodel
      stop
    end select

  end function tauhyd3D

  function taudev3D ( imodel, x )

    integer, intent(in) :: imodel
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(6) :: taudev3D
    real(dp) :: rate, t

    rate = gammadot
    t = time

    select case ( imodel )
    case ( 2:3 )
    taudev3D(1) = -1 + (x(1) - rate*t*x(2))**2/((x(1) - rate*t*x(2))**2 + &
      (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2) + (-2*rate*t* &
      x(1)*x(2) - 2*rate**3*t**3*x(1)*x(2) + x(2)**2 + rate**4*t**4*x(2)**2 + &
      x(3)**2 + rate**2*t**2*(x(1)**2 + 2*x(2)**2 + x(3)**2))/((x(1)**2 - 2* &
      rate*t*x(1)*x(2) + (1 + rate**2*t**2)*x(2)**2 + x(3)**2)*Sqrt(((x(1) - &
      rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)/ &
      (x(2)**2 + (x(1) - rate*t*x(2))**2 + x(3)**2)))
    taudev3D(2) = ((x(1) - rate*t*x(2))*(-(rate*t*x(1)) + x(2) + rate**2* &
      t**2*x(2)))/((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t* &
      x(2)))**2 + x(3)**2) + (-(x(1)*x(2)) - 2*rate**2*t**2*x(1)*x(2) + &
      rate**3*t**3*x(2)**2 + rate*t*(x(1)**2 + x(2)**2 + x(3)**2))/((x(1)**2 - &
      2*rate*t*x(1)*x(2) + (1 + rate**2*t**2)*x(2)**2 + x(3)**2)*Sqrt(((x(1) - &
      rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)/ &
      (x(2)**2 + (x(1) - rate*t*x(2))**2 + x(3)**2)))
    taudev3D(3) = x(3)*((x(1) - rate*t*x(2))/((x(1) - rate*t*x(2))**2 + (x(2) +&
      rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2) - x(1)/((x(1)**2 - 2*rate*t* &
      x(1)*x(2) + (1 + rate**2*t**2)*x(2)**2 + x(3)**2)*Sqrt(((x(1) - rate*t* &
      x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 +&
      (x(1) - rate*t*x(2))**2 + x(3)**2))))
    taudev3D(4) = -1 + (-(rate*t*x(1)) + x(2) + rate**2*t**2*x(2))**2/((x(1) - &
      rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2) + &
      (x(1)**2 - 2*rate*t*x(1)*x(2) + rate**2*t**2*x(2)**2 + x(3)**2)/ &
      ((x(1)**2 - 2*rate*t*x(1)*x(2) + (1 + rate**2*t**2)*x(2)**2 + x(3)**2)* &
      Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t* &
      x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + x(3)**2)))
    taudev3D(5) = x(3)*((-(rate*t*x(1)) + x(2) + rate**2*t**2*x(2))/((x(1) - &
      rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2) - &
      x(2)/((x(1)**2 - 2*rate*t*x(1)*x(2) + (1 + rate**2*t**2)*x(2)**2 + &
      x(3)**2)*Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + rate*t*(-x(1) + rate*t* &
      x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t*x(2))**2 + x(3)**2))))
    taudev3D(6) = -1 + x(3)**2/((x(1) - rate*t*x(2))**2 + (x(2) + rate*t* &
      (-x(1) + rate*t*x(2)))**2 + x(3)**2) + (x(1)**2 - 2*rate*t*x(1)*x(2) + &
      (1 + rate**2*t**2)*x(2)**2)/((x(1)**2 - 2*rate*t*x(1)*x(2) + (1 + &
      rate**2*t**2)*x(2)**2 + x(3)**2)*Sqrt(((x(1) - rate*t*x(2))**2 + (x(2) + &
      rate*t*(-x(1) + rate*t*x(2)))**2 + x(3)**2)/(x(2)**2 + (x(1) - rate*t* &
      x(2))**2 + x(3)**2)))
    case ( 4:5 )
    taudev3D(1) = (rate*t*(-4*rate**6*t**6*x(1)*x(2)**3 + rate**7*t**7* &
      x(2)**4 - 2*x(1)*x(2)*(x(2)**2 + x(3)**2) + rate**5*t**5*x(2)**2*(6* &
      x(1)**2 + 4*x(2)**2 + x(3)**2) - 2*rate**4*t**4*x(1)*x(2)*(2*x(1)**2 + &
      6*x(2)**2 + x(3)**2) - 2*rate**2*t**2*x(1)*x(2)*(2*x(1)**2 + 5*x(2)**2 + &
      4*x(3)**2) + rate**3*t**3*(x(1)**4 + 5*x(2)**2*(x(2)**2 + x(3)**2) + &
      x(1)**2*(12*x(2)**2 + x(3)**2)) + rate*t*(2*x(2)**4 + 3*x(2)**2*x(3)**2 +&
      x(3)**4 + x(1)**2*(5*x(2)**2 + 3*x(3)**2))))/(2*((1 + rate**2*t**2)* &
      x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + &
      rate**4*t**4)*x(2)**2 + x(3)**2)**2)
    taudev3D(2) = (rate*t*(-4*rate**5*t**5*x(1)*x(2)**3 + rate**6*t**6* &
      x(2)**4 + 2*(x(1)**2 + x(3)**2)*(x(2)**2 + x(3)**2) + rate**4*t**4* &
      x(2)**2*(6*x(1)**2 + 3*x(2)**2 + x(3)**2) - rate**3*t**3*x(1)*x(2)*(4* &
      x(1)**2 + 9*x(2)**2 + 2*x(3)**2) - rate*t*x(1)*x(2)*(3*x(1)**2 + 4* &
      x(2)**2 + 7*x(3)**2) + rate**2*t**2*(x(1)**4 + 2*x(2)**4 + 5*x(2)**2* &
      x(3)**2 + x(1)**2*(9*x(2)**2 + x(3)**2))))/(2*((1 + rate**2*t**2)* &
      x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2* &
      t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**2)
    taudev3D(3) = -(rate*t*x(3)*(-(rate**3*t**3*x(1)*x(2)**2) + rate**4*t**4*&
      x(2)**3 + rate*t*x(1)*(x(1)**2 - 4*x(2)**2 - x(3)**2) + 2*x(2)* &
      (x(2)**2 + x(3)**2) + rate**2*t**2*x(2)*(-x(1)**2 + 4*x(2)**2 + &
      x(3)**2)))/(2*((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2* &
      t**2)*x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + &
      x(3)**2)**2)
    taudev3D(4) = (rate*t*(x(1)**2 - 2*rate*t*x(1)*x(2) + rate**2*t**2* &
      x(2)**2 + x(3)**2)*(-2*x(1)*x(2) - 2*rate**2*t**2*x(1)*x(2) + rate**3* &
      t**3*x(2)**2 + rate*t*(x(1)**2 + 2*x(2)**2 - x(3)**2)))/(2*((1 + rate**2*&
      t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2* &
      t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**2)
    taudev3D(5) = (rate*t*(-((2 + rate**2*t**2)*x(1)) + rate*t*(3 + rate**2* &
      t**2)*x(2))*x(3)*(x(1)**2 - 2*rate*t*x(1)*x(2) + rate**2*t**2*x(2)**2 + &
      x(3)**2))/(2*((1 + rate**2*t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)* &
      x(1)*x(2) + (1 + 3*rate**2*t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**2)
    taudev3D(6) = -(rate*t*(-4*rate**6*t**6*x(1)*x(2)**3 + rate**7*t**7* &
      x(2)**4 - 8*rate**2*t**2*x(1)*x(2)*(x(1)**2 + 2*x(2)**2 + x(3)**2) + &
      rate**5*t**5*x(2)**2*(6*x(1)**2 + 5*x(2)**2 + x(3)**2) - 2*rate**4*t**4*&
      x(1)*x(2)*(2*x(1)**2 + 8*x(2)**2 + x(3)**2) - 2*x(1)*x(2)*(x(1)**2 + &
      x(2)**2 + 2*x(3)**2) + rate**3*t**3*(x(1)**4 + 7*x(2)**4 + 5*x(2)**2* &
      x(3)**2 + x(1)**2*(18*x(2)**2 + x(3)**2)) + rate*t*(x(1)**4 + 2*x(2)**4 +&
      5*x(2)**2*x(3)**2 + x(1)**2*(11*x(2)**2 + 3*x(3)**2))))/(2*((1 + rate**2*&
      t**2)*x(1)**2 - 2*rate*t*(2 + rate**2*t**2)*x(1)*x(2) + (1 + 3*rate**2* &
      t**2 + rate**4*t**4)*x(2)**2 + x(3)**2)**2)
    case default
      write(*,'(/a,i0/)') 'Error taudev3D: wrong value imodel: ', imodel
      stop
    end select

  end function taudev3D

  function func_error_mansol_3D ( nr, x )

    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func_error_mansol_3D
    real(dp) :: lambda, x0, y0, z0

    select case(nr)
      case(2) ! shear
        x0 = x(1) - gammadot * x(2) * time
        y0 = x(2)
        z0 = x(3)
        lambda = R0 / sqrt( x0**2 + y0**2 + z0**2 )
        func_error_mansol_3D = (lambda-1)* sqrt( x(1)**2 + x(2) ** 2 + x(3) ** 2 )
      case default
        write(*,'(/a,i0/)') 'Error func_error_mansol_3D: wrong function number: ', &
          nr
        stop
    end select

  end function func_error_mansol_3D

end module functions_mansol_m
