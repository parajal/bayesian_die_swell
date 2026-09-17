
! Copyright (C) 2009-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Utility routines for the element subroutines

module element_utils_m

  use element_defs_m
  use vector_m, only: get_vector
  use elvector_m, only: get_elvector

  implicit none


  interface evaluate_coefficient
    module procedure evaluate_scalar_coefficient_simple, &
      evaluate_scalar_coefficient, evaluate_vector_coefficient, &
      evaluate_tensor_coefficient
  end interface evaluate_coefficient

contains


! evaluate a scalar coefficient as function of position (simple)

  subroutine evaluate_scalar_coefficient_simple ( choice, value, func, funcnr, &
    x, coef )

!   The choice parameter for evaluating the scalar coefficient coef:
!      0: constant, given by value
!      1: given by a function (func) with specified funcnr
    integer, intent(in) :: choice

!   constant value for coef (choice=0)
    real(dp), intent(in) :: value

!   specification of the function (choice=1)
    interface
      function func ( nr, x ) result(f)
        use kind_defs_m
        implicit none
        integer, intent(in) :: nr
        real(dp), intent(in), dimension(:) :: x
        real(dp) :: f
      end function func
    end interface

!   function number (choice=1)
    integer, intent(in) :: funcnr

!   coordinates used for the function func: x(1:np,1:ndim), where np is the
!   number of points (choice=1)
    real(dp), dimension(:,:), intent(in) :: x

!   the scalar computed in np points: coef(1:np)
    real(dp), dimension(:), intent(out) :: coef


    integer :: ip


    select case ( choice )

      case(0) ! coef is constant

        coef = value

      case(1) ! coef given by function

        do ip = 1, size(coef)
          coef(ip) = func ( funcnr, x(ip,:) )
        end do

      case default

        write(*,'(/a/a,i0/)') 'Error in evaluate_scalar_coefficient_simple:', &
        ' incorrect value of choice = ', choice
        stop

    end select

  end subroutine evaluate_scalar_coefficient_simple


! evaluate a scalar coefficient as function of position

  subroutine evaluate_scalar_coefficient ( mesh, problem, oldvectors, elgrp, &
    elem, choice, value, func, funcnr, x, indx_v, indx_v1, indx_p1, layer, &
    phi, indx_e, coef )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(oldvectors_t), intent(in) :: oldvectors

!   element group and element number
    integer, intent(in) :: elgrp, elem

!   The choice parameter for evaluating the scalar coefficient coef:
!      0: constant, given by value
!      1: given by a function (func) with specified funcnr
!      2: given by nodal values in the specified vector
!      3: given by the vector per element elvectors
    integer, intent(in) :: choice

!   constant value for coef (choice=0)
    real(dp), intent(in) :: value

!   specification of the function (choice=1)
    interface
      function func ( nr, x )
        use kind_defs_m
        implicit none
        integer, intent(in) :: nr
        real(dp), intent(in), dimension(:) :: x
        real(dp) :: func
      end function func
    end interface

!   function number (choice=1)
    integer, intent(in) :: funcnr

!   coordinates used for the function func: x(1:np,1:ndim), where np is the
!   number of points (choice=1)
    real(dp), dimension(:,:), intent(in) :: x

!   index for the vector for evaluating the coefficient for choice=2:
!   vector=oldvectors%v(indx_v)%p
!     or
!   vector=oldvectors%v1(indx_v1)%p(indx_p1)
    integer, intent(in), optional :: indx_v, indx_v1, indx_p1

!   if layer > 0 the element degrees of freedom in vector are restricted
!   to the specified layer (for choice=2).
!   For layer=0 this parameter is ignored.
    integer, intent(in) :: layer

!   shape function used for choice=2
    real(dp), dimension(:,:), intent(in) :: phi

!   index for the elvector for evaluating the coefficient for choice=3:
!   elvector=oldvectors%e(indx_e)%p
    integer, intent(in) :: indx_e

!   the scalar computed in np points: coef(1:np)
    real(dp), dimension(:), intent(out) :: coef


    integer :: ip
    real(dp), dimension(size(phi,2)) :: u


    select case ( choice )

      case(0) ! coef is constant

        coef = value

      case(1) ! coef given by function

        do ip = 1, size(coef)
          coef(ip) = func ( funcnr, x(ip,:) )
        end do

      case(2) ! coef given by nodal point values

        if ( present(indx_v) ) then

          call get_vector ( mesh, problem, oldvectors%v(indx_v)%p, &
            elgrp, elem, u, layer=layer )

        else if ( present(indx_v1) ) then

          call get_vector ( mesh, problem, oldvectors%v1(indx_v1)%p(indx_p1), &
            elgrp, elem, u, layer=layer )

        end if

        coef = matmul ( phi, u )

      case(3) ! coef given in points per element

        call get_elvector ( mesh, oldvectors%e(indx_e)%p, elgrp, elem, r1=coef )

      case default

        write(*,'(/a/a,i0/)') 'Error in evaluate_scalar_coefficient:', &
        ' incorrect value of choice = ', choice
        stop

    end select

  end subroutine evaluate_scalar_coefficient


! evaluate a vector coefficient as function of position (simple)

  subroutine evaluate_vector_coefficient_simple ( choice, value, vfunc, &
    vfuncnr, x, coef )

!   The choice parameter for evaluating the vector coefficient coef:
!      0: constant, given by value
!      1: given by a function (vfunc) with specified funcnr
    integer, intent(in) :: choice

!   constant value for coef (choice=0)
    real(dp), dimension(:), intent(in) :: value

!   specification of the function (choice=1)
    interface
      function vfunc ( n, nr, x ) result(v)
        use kind_defs_m
        implicit none
        integer, intent(in) :: n, nr
        real(dp), intent(in), dimension(:) :: x
        real(dp), dimension(n) :: v
      end function vfunc
    end interface

!   function number (choice=1)
    integer, intent(in) :: vfuncnr

!   coordinates used for the function func: x(1:np,1:ndim), where np is the
!   number of points (choice=1)
    real(dp), dimension(:,:), intent(in) :: x

!   the vector computed in np points: coef(1:np,:)
    real(dp), dimension(:,:), intent(out) :: coef


    integer :: ip, np, nd


    np = size(coef,1)
    nd = size(coef,2)

    select case ( choice )

      case(0) ! coef is constant

        do ip = 1, np
          coef(ip,:) = value
        end do

      case(1) ! coef given by function

        do ip = 1, np
          coef(ip,:) = vfunc ( nd, vfuncnr, x(ip,:) )
        end do

      case default

        write(*,'(/a/a,i0/)') 'Error in evaluate_vector_coefficient_simple:', &
        ' incorrect value of choice = ', choice
        stop

    end select

  end subroutine evaluate_vector_coefficient_simple


! evaluate a vector coefficient as function of position

  subroutine evaluate_vector_coefficient ( mesh, problem, oldvectors, elgrp, &
    elem, choice, value, vfunc, vfuncnr, x, indx_v, indx_v1, indx_p1, layer, &
    phi, indx_e, coef )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(oldvectors_t), intent(in) :: oldvectors

!   element group and element number
    integer, intent(in) :: elgrp, elem

!   The choice parameter for evaluating the vector coefficient coef:
!      0: constant, given by value
!      1: given by a function (vfunc) with specified funcnr
!      2: given by nodal values in the specified vector
!      3: given by the vector per element elvectors
    integer, intent(in) :: choice

!   constant value for coef (choice=0)
    real(dp), dimension(:), intent(in) :: value

!   specification of the function (choice=1)
    interface
      function vfunc ( n, nr, x )
        use kind_defs_m
        implicit none
        integer, intent(in) :: n, nr
        real(dp), intent(in), dimension(:) :: x
        real(dp), dimension(n) :: vfunc
      end function vfunc
    end interface

!   function number (choice=1)
    integer, intent(in) :: vfuncnr

!   coordinates used for the function func: x(1:np,1:ndim), where np is the
!   number of points (choice=1)
    real(dp), dimension(:,:), intent(in) :: x

!   index for the vector for evaluating the coefficient for choice=2:
!   vector=oldvectors%v(indx_v)%p
!     or
!   vector=oldvectors%v1(indx_v1)%p(indx_p1)
    integer, intent(in), optional :: indx_v, indx_v1, indx_p1

!   if layer > 0 the element degrees of freedom in vector are restricted
!   to the specified layer (for choice=2).
!   For layer=0 this parameter is ignored.
    integer, intent(in) :: layer

!   shape function used for choice=2
    real(dp), dimension(:,:), intent(in) :: phi

!   index for the elvector for evaluating the coefficient for choice=3:
!   elvector=oldvectors%e(indx_e)%p
    integer, intent(in) :: indx_e

!   the vector computed in np points: coef(1:np,:)
    real(dp), dimension(:,:), intent(out) :: coef


    integer :: ip, np, nd
    real(dp), dimension(size(phi,2)*size(coef,2)) :: u


    np = size(coef,1)
    nd = size(coef,2)

    select case ( choice )

      case(0) ! coef is constant

        do ip = 1, np
          coef(ip,:) = value
        end do

      case(1) ! coef given by function

        do ip = 1, np
          coef(ip,:) = vfunc ( nd, vfuncnr, x(ip,:) )
        end do

      case(2) ! coef given by nodal point values

        if ( present(indx_v) ) then

          call get_vector ( mesh, problem, oldvectors%v(indx_v)%p, &
            elgrp, elem, u, layer=layer )

        else if ( present(indx_v1) ) then

          call get_vector ( mesh, problem, oldvectors%v1(indx_v1)%p(indx_p1), &
            elgrp, elem, u, layer=layer )

        end if

        coef = matmul ( phi, reshape ( u, [size(phi,2),nd] ) )

      case(3) ! coef given in points per element

        call get_elvector ( mesh, oldvectors%e(indx_e)%p, elgrp, elem, r2=coef )

      case default

        write(*,'(/a/a,i0/)') 'Error in evaluate_vector_coefficient:', &
        ' incorrect value of choice = ', choice
        stop

    end select

  end subroutine evaluate_vector_coefficient


! evaluate a tensor coefficient as function of position

  subroutine evaluate_tensor_coefficient ( mesh, problem, oldvectors, elgrp, &
    elem, choice, value, tfunc, tfuncnr, x, indx_v, indx_v1, indx_p1, layer, &
    phi, indx_e, coef )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(oldvectors_t), intent(in) :: oldvectors

!   element group and element number
    integer, intent(in) :: elgrp, elem

!   The choice parameter for evaluating the tensor coefficient coef:
!      0: constant, given by value
!      1: given by a function (vfunc) with specified funcnr
!      2: given by nodal values in the specified vector
!      3: given by the vector per element elvectors
    integer, intent(in) :: choice

!   constant value for coef (choice=0)
    real(dp), dimension(:,:), intent(in) :: value

!   specification of the function (choice=1)
    interface
      function tfunc ( n, nr, x )
        use kind_defs_m
        implicit none
        integer, intent(in) :: n, nr
        real(dp), intent(in), dimension(:) :: x
        real(dp), dimension(n,n) :: tfunc
      end function tfunc
    end interface

!   function number (choice=1)
    integer, intent(in) :: tfuncnr

!   coordinates used for the function func: x(1:np,1:ndim), where np is the
!   number of points (choice=1)
    real(dp), dimension(:,:), intent(in) :: x

!   index for the vector for evaluating the coefficient for choice=2:
!   vector=oldvectors%v(indx_v)%p
!     or
!   vector=oldvectors%v1(indx_v1)%p(indx_p1)
    integer, intent(in), optional :: indx_v, indx_v1, indx_p1

!   if layer > 0 the element degrees of freedom in vector are restricted
!   to the specified layer (for choice=2).
!   For layer=0 this parameter is ignored.
    integer, intent(in) :: layer

!   shape function used for choice=2
    real(dp), dimension(:,:), intent(in) :: phi

!   index for the elvector for evaluating the coefficient for choice=3:
!   elvector=oldvectors%e(indx_e)%p
    integer, intent(in) :: indx_e

!   the tensor computed in np points: coef(1:np,:,:)
    real(dp), dimension(:,:,:), intent(out) :: coef


    integer :: ip, nd, np, id
    real(dp), dimension(size(phi,2)*size(coef,2)**2) :: u
    real(dp), dimension(size(phi,2),size(coef,2),size(coef,2)) :: tmp

    np = size(coef,1)
    nd = size(coef,2)

    select case ( choice )

      case(0) ! coef is constant

        do ip = 1, np
          coef(ip,:,:) = value
        end do

      case(1) ! coef given by function

        do ip = 1, np
          coef(ip,:,:) = tfunc ( nd, tfuncnr, x(ip,:) )
        end do

      case(2) ! coef given by nodal point values

        if ( present(indx_v) ) then

          call get_vector ( mesh, problem, oldvectors%v(indx_v)%p, &
            elgrp, elem, u, layer=layer )

        else if ( present(indx_v1) ) then

          call get_vector ( mesh, problem, oldvectors%v1(indx_v1)%p(indx_p1), &
            elgrp, elem, u, layer=layer )

        end if

        tmp = reshape ( u, [size(phi,2),nd,nd] )

        do id = 1, nd
          coef(:,:,id) = matmul ( phi, tmp(:,:,id) )
        end do

      case(3) ! coef given in points per element

        call get_elvector ( mesh, oldvectors%e(indx_e)%p, elgrp, elem, r3=coef )

      case default

        write(*,'(/a/a,i0/)') 'Error in evaluate_tensor_coefficient:', &
        ' incorrect value of choice = ', choice
        stop

    end select

  end subroutine evaluate_tensor_coefficient

end module element_utils_m
