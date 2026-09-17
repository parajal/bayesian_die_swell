module subs_m

  use tfem_elem_m

  implicit none

  real(dp) :: ltime

contains


! functions for the source term

  function sourcefunc ( nr, x )
    use math_defs_m
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: sourcefunc, sigma

    select case(nr)
      case(1)
        sigma = 0.05_dp
        sourcefunc = exp(-(x(1)+x(2)-1)**2/(2*sigma**2))/(sqrt(pi)*sigma)
        sourcefunc = maxval ( [sourcefunc, 1e-12_dp] )
      case(2)
        sourcefunc = sin(2*x(1)+x(2)-ltime) ! exact q1
      case(3)
        sourcefunc = 2*cos(2*x(1)+x(2)-ltime) ! f q1
      case(4)
        sourcefunc = cos(x(1)+2*x(2)-ltime) ! exact q2
      case(5)
        sourcefunc = -2*sin(x(1)+2*x(2)-ltime) - sin(2*x(1)+x(2)-ltime) ! f q2
    case default
      write(*,'(/a,i0/)') 'Error sourcefunc: wrong function number: ', nr
      stop
    end select

  end function sourcefunc


! Internal element routine to fill the scalar right-hand side terms of the
! energy equation in the Gauss points. This element should be used together
! with the routine loop_over_elements.

  subroutine fill_f_elem ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors )

    use stokes_elements_m
    use stokes_globals_m
    use convection_diffusion_supg_elements_m
    use convection_diffusion_supg_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors

    integer :: ip

    if ( first ) then

!     first element
!     set globals, gauss, shapefunctions, ...

      call set_stokes_elem ( mesh, problem, elgrp, elem, &
        first, last, coefficients, oldvectors )

      call set_conv_diff_supg_elem ( mesh, problem, elgrp, elem, &
        first, last, coefficients, oldvectors )

      allocate ( q_n(ndf), q_ng(ninti) )
      allocate ( f_g(ninti) )

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_coordinates ( x, phi, xg )

!   get value of q

    call get_sysvector ( mesh, problem, oldvectors%s(4)%p, elgrp, &
      elem, q_n, physq=[1] )

    q_ng = matmul ( chi, q_n )

!   add the manufactured solution f2
    do ip = 1, ninti
      f_g(ip) = sourcefunc ( nr=5, x=xg(ip,:) ) + q_ng(ip)
    end do

    call put_elvector ( mesh, oldvectors%e(2)%p, elgrp, elem, r1=f_g )

    if ( last ) then

!     last element
!     unset globals, gauss, shapefunctions, ...

      call unset_stokes_elem ( last, coefficients )

      call unset_conv_diff_supg_elem ( last, coefficients )

      deallocate ( q_n, q_ng )
      deallocate ( f_g )

    end if

  end subroutine fill_f_elem


! Internal element routine to fill multiple scalar right-hand sides terms of the
! energy equation in the Gauss points. This element should be used together
! with the routine loop_over_elements.

  subroutine fill_f_elem_multi ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors )

    use stokes_elements_m
    use stokes_globals_m
    use convection_diffusion_supg_elements_m
    use convection_diffusion_supg_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors

    integer :: ip

    if ( first ) then

!     first element
!     set globals, gauss, shapefunctions, ...

      call set_stokes_elem ( mesh, problem, elgrp, elem, &
        first, last, coefficients, oldvectors )

      call set_conv_diff_supg_elem ( mesh, problem, elgrp, elem, &
        first, last, coefficients, oldvectors )

      allocate ( q_n(ndf), q_ng(ninti) )
      allocate ( f_g(ninti) )

    end if

!   start of element

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_coordinates ( x, phi, xg )

!   RHS of the first equation

    do ip = 1, ninti
      f_g(ip) = sourcefunc ( nr=3, x=xg(ip,:) )
    end do

    call put_elvector ( mesh, oldvectors%e(2)%p, elgrp, elem, r1=f_g, nr=1 )

!   RHS of the second equation

!   get value of q1hat

    call get_sysvector ( mesh, problem, oldvectors%s1(3)%p(1), elgrp, &
      elem, q_n, physq=[1] )

    q_ng = matmul ( chi, q_n )

!   add the manufactured solution f2
    do ip = 1, ninti
      f_g(ip) = sourcefunc ( nr=5, x=xg(ip,:) ) + q_ng(ip)
    end do

    call put_elvector ( mesh, oldvectors%e(2)%p, elgrp, elem, r1=f_g, nr=2 )

    if ( last ) then

!     last element
!     unset globals, gauss, shapefunctions, ...

      call unset_stokes_elem ( last, coefficients )

      call unset_conv_diff_supg_elem ( last, coefficients )

      deallocate ( q_n, q_ng )
      deallocate ( f_g )

    end if

  end subroutine fill_f_elem_multi


end module subs_m
