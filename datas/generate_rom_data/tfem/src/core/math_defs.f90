
! Copyright (C) 2019-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Define common math constants and functions in one place.
! Includes kind_defs_m module

module math_defs_m

  use kind_defs_m

  implicit none

  real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

! interfaces to support generic interfaces

  interface norm
    module procedure norm, norm_mat_dim, norm_mat
  end interface norm

  interface sumsq
    module procedure sumsq, sumsq_mat_dim, sumsq_mat
  end interface sumsq

contains

! cross product function

  function cross_product ( a, b )

    real(dp), intent(in), dimension(3) :: a, b
    real(dp), dimension(3) :: cross_product

    cross_product = [ a(2) * b(3) - a(3) * b(2), &
                      a(3) * b(1) - a(1) * b(3), &
                      a(1) * b(2) - a(2) * b(1) ]

  end function cross_product


! norm of a vector

  function norm ( uvector )

    real(dp), intent(in), dimension(:) :: uvector
    real(dp) :: norm

    norm = norm2 ( uvector )

  end function norm


! norm of columns (dim=1) or rows (dim=2) of a matrix

  function norm_mat_dim ( umatrix, dim )

    real(dp), intent(in), dimension(:,:) :: umatrix
    integer, intent(in) :: dim
    real(dp), dimension(:), allocatable :: norm_mat_dim

    norm_mat_dim = norm2 ( umatrix, dim=dim )

  end function norm_mat_dim


! norm of a matrix

  function norm_mat ( umatrix )

    real(dp), intent(in), dimension(:,:) :: umatrix
    real(dp) :: norm_mat

    norm_mat = norm2 ( umatrix )

  end function norm_mat


! sum of squares of a vector

  function sumsq ( uvector )

    real(dp), intent(in), dimension(:) :: uvector
    real(dp) :: sumsq

    sumsq = sum ( uvector**2 )

  end function sumsq


! sum of squares of columns (dim=1) or rows (dim=2) of a matrix

  function sumsq_mat_dim ( umatrix, dim )

    real(dp), intent(in), dimension(:,:) :: umatrix
    integer, intent(in) :: dim
    real(dp), dimension(:), allocatable :: sumsq_mat_dim

    sumsq_mat_dim = sum ( umatrix**2, dim=dim )

  end function sumsq_mat_dim


! sum of squares of a matrix

  function sumsq_mat ( umatrix )

    real(dp), intent(in), dimension(:,:) :: umatrix
    real(dp) :: sumsq_mat

    sumsq_mat = sum ( umatrix**2 )

  end function sumsq_mat


! Performs a Cholesky decomposition of a symmetric positive definite matrix
! Only 2x2 or 3x3 matrices available

  function chol ( a )

!   matrix input, lower triangle column wise
    real(dp), dimension(:), intent(in) :: a

!   Cholesky lower triangular column wise output
    real(dp), dimension(size(a)) :: chol

    real(dp) :: tmp, b(size(a))

    select case ( size(a) )
    case(3) ! 2x2 matrix

      if ( a(1) < 0 ) then
        write(*,'(/a/)') 'Error chol: a11 < 0'
        stop
      end if
      b(1) = sqrt(a(1))
      b(2) = a(2) / b(1)
      tmp = a(3) - b(2)**2
      if ( tmp < 0 ) then
        write(*,'(/a/)') 'Error chol: a22 - l21^2 < 0'
        stop
      end if
      b(3) = sqrt(tmp)

    case(6) ! 3x3 matrix

      if ( a(1) < 0 ) then
        write(*,'(/a/)') 'Error chol: a11 < 0'
        stop
      end if
      b(1) = sqrt(a(1))
      b(2) = a(2) / b(1)
      b(3) = a(3) / b(1)
      tmp = a(4) - b(2)**2
      if ( tmp < 0 ) then
        write(*,'(/a/)') 'Error chol: a22 - l21^2 < 0'
        stop
      end if
      b(4) = sqrt(tmp)
      b(5) = (a(5)-b(3)*b(2))/b(4)
      tmp = a(6)-b(3)**2-b(5)**2
      if ( tmp < 0 ) then
        write(*,'(/a/)') 'Error chol: a33 - l31^2 - l32^2 < 0'
        stop
      end if
      b(6) = sqrt(tmp)

    case default

      write(*,'(/a,i0/)') 'Error chol: invalid size(a) = ', size(a)

    end select

    chol = b

  end function chol

end module math_defs_m
