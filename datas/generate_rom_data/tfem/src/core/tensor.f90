
! Copyright (C) 2023-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, routines, for tensors

module tensor_m

  use kind_defs_m

  implicit none

contains


! Convert tensor2 to vector

  function tensor2_to_vector ( ndim, A ) result(B)

!   dimension of space
    integer, intent(in) :: ndim

!   tensor2 T_ij with i,j=1,...,ndim
    real(dp), dimension(:,:), intent(in) :: A

!   vector with ordering according to a row-major tensor2, i.e.
!   ndim=2: 11, 12, 21, 22
!   ndim=3: 11, 12, 13, 21, 22, 23, 31, 32, 33
    real(dp), dimension(ndim**2) :: B

    if ( any( ndim /= shape(A) ) ) then
      write(*,'(/2(a/))') ' Error in tensor2_to_vector: ', &
        ' Extend of the dimensions of A not equal to specified ndim.'
      stop
    end if

    B = reshape ( transpose(A), [ndim**2] )

  end function tensor2_to_vector


! Convert vector to tensor2

  function vector_to_tensor2 ( ndim, A ) result(B)

!   dimension of space
    integer, intent(in) :: ndim

!   vector with ordering according to a row-major tensor2, i.e.
!   ndim=2: 11, 12, 21, 22
!   ndim=3: 11, 12, 13, 21, 22, 23, 31, 32, 33
    real(dp), dimension(:), intent(in) :: A

!   tensor2 T_ij with i,j=1,...,ndim
    real(dp), dimension(ndim,ndim) :: B

    if ( size(A) /= ndim**2 ) then
      write(*,'(/2(a/))') ' Error in vector_to_tensor2: ', &
        ' Dimension of A not equal to specified ndim.'
      stop
    end if

    B = reshape ( A, [ndim,ndim], order=[2,1] )

  end function vector_to_tensor2


! Convert tensor2 to vector, with tensor assumed symmetric

  function tensor2_to_vector_symmetric ( ndim, A ) result(B)

!   dimension of space
    integer, intent(in) :: ndim

!   tensor2 T_ij with i,j=1,..,ndim
!   NOTE: it is assumed that the 2nd-order tensor is symmetric, i.e.
!         T_ij=T_jl and only the upper-diagonal of the tensor is
!         used in the vector. The symmetry is not checked.
    real(dp), dimension(:,:), intent(in) :: A

!   vector with ordering according to a symmetric row-major tensor2, i.e.
!   ndim=2: 11, 12, 22
!   ndim=3: 11, 12, 13, 22, 23, 33
!   NOTE: Only components corresponding the upper-triangle of the tensor
!         are retained and the "lower-triangle" of the tensor is
!         discarded. This means that the resulting matrix does not represent
!         the full tensorial properties since the lower-diagonal is missing.
    real(dp), dimension(ndim*(ndim+1)/2) :: B

    integer :: ii, i, j

    if ( any( ndim /= shape(A) ) ) then
      write(*,'(/2(a/))') ' Error in tensor2_to_vector_symmetric: ', &
        ' Extend of the dimensions of tensor A not equal to specified ndim.'
      stop
    end if

    ii = 0
    do i = 1, ndim
      do j = i, ndim
        ii = ii + 1
        B(ii) = A(i,j)
      end do
    end do

  end function tensor2_to_vector_symmetric


! Convert vector to tensor2, with tensor assumed symmetric

  function vector_to_tensor2_symmetric ( ndim, A ) result(B)

!   dimension of space
    integer, intent(in) :: ndim

!   vector with ordering according to a symmetric row-major tensor2, i.e.
!   ndim=2: 11, 12, 22
!   ndim=3: 11, 12, 13, 22, 23, 33
    real(dp), dimension(:), intent(in) :: A

!   tensor2 T_ij with i,j=1,..,ndim
!   NOTE: it is assumed that the 2nd-order tensor is symmetric, i.e.
!         T_ij=T_jl and only the upper-diagonal of the tensor is taken from
!         the vector. The lower-diagonal is filled using symmetry of the tensor
    real(dp), dimension(ndim,ndim) :: B

    integer :: ii, i, j

    if ( size(A) /= ndim*(ndim+1)/2 ) then
      write(*,'(/2(a/))') ' Error in vector_to_tensor2_symmetric: ', &
        ' Dimension of A not consistent with to specified ndim.'
      stop
    end if

!   upper diagonal

    ii = 0
    do i = 1, ndim
      do j = i, ndim
        ii = ii + 1
        B(i,j) = A(ii)
      end do
    end do

!   assume symmetry

    do i = 1, ndim
      do j = i+1, ndim
        B(j,i) = B(i,j)
      end do
    end do

  end function vector_to_tensor2_symmetric


! Convert tensor4 to matrix

  function tensor4_to_matrix ( ndim, A ) result(B)

!   dimension of space
    integer, intent(in) :: ndim

!   tensor4 T_ijkl with i,j,k,l=1,...,ndim
    real(dp), dimension(:,:,:,:), intent(in) :: A

!   matrix with rows/columns ordering according to a row-major tensor2, i.e.
!   ndim=2: 11, 12, 21, 22
!   ndim=3: 11, 12, 13, 21, 22, 23, 31, 32, 33
    real(dp), dimension(ndim**2,ndim**2) :: B

    integer :: ii, j, l
    integer, dimension(:,:), allocatable :: indx

    if ( any( ndim /= shape(A) ) ) then
      write(*,'(/2(a/))') ' Error in tensor4_to_matrix: ', &
        ' Extend of the dimensions of A not equal to specified ndim.'
      stop
    end if

    indx = reshape ( [(ii,ii=1,ndim**2)], [ndim,ndim], order=[2,1] )

    do l = 1, ndim
      do j = 1, ndim
        B(indx(:,j),indx(:,l)) = A(:,j,:,l)
      end do
    end do

  end function tensor4_to_matrix


! Convert matrix to tensor4

  function matrix_to_tensor4 ( ndim, A ) result(B)

!   dimension of space
    integer, intent(in) :: ndim

!   matrix with rows/columns ordering according to a row-major tensor2, i.e.
!   ndim=2: 11, 12, 21, 22
!   ndim=3: 11, 12, 13, 21, 22, 23, 31, 32, 33
    real(dp), dimension(:,:), intent(in) :: A

!   tensor4 T_ijkl with i,j,k,l=1,...,ndim
    real(dp), dimension(ndim,ndim,ndim,ndim) :: B

    integer :: ii, j, l
    integer, dimension(:,:), allocatable :: indx

    if ( any( ndim**2 /= shape(A) ) ) then
      write(*,'(/2(a/))') ' Error in matrix_to_tensor4: ', &
        ' Extend of the dimensions of A not equal to specified ndim**2.'
      stop
    end if

    indx = reshape ( [(ii,ii=1,ndim**2)], [ndim,ndim], order=[2,1] )

    do l = 1, ndim
      do j = 1, ndim
         B(:,j,:,l) = A(indx(:,j),indx(:,l))
      end do
    end do

  end function matrix_to_tensor4


! Convert tensor4 to matrix, with "left" tensor assumed symmetric

  function tensor4_to_matrix_left_symmetric ( ndim, A ) result(B)

!   dimension of space
    integer, intent(in) :: ndim

!   tensor4 T_ijkl with i,j,k,l=1,...,ndim
!   NOTE: it is assumed that the 4th-order tensor is left-symmetric, i.e.
!         T_ijkl=T_jikl and only the upper-diagonal of the "left" tensor is
!         used in the matrix. The left symmetry is not checked.
    real(dp), dimension(:,:,:,:), intent(in) :: A

!   matrix with rows ordering according to a symmetric row-major tensor2, i.e.
!   ndim=2: 11, 12, 22
!   ndim=3: 11, 12, 13, 22, 23, 33
!   and with columns ordering according to a row-major tensor2, i.e.
!   ndim=2: 11, 12, 21, 22
!   ndim=3: 11, 12, 13, 21, 22, 23, 31, 32, 33
!   NOTE: Only the rows corresponding the upper-triangle of the "left" tensor
!         are retained and the "lower-triangle" of the "left" tensor is
!         discarded. This means that the resulting matrix does not represent
!         the full tensorial properties since the "lower-diagonal" is missing.
    real(dp), dimension(ndim*(ndim+1)/2,ndim**2) :: B

    integer :: ii, i, j, l, is(ndim*(ndim+1)/2)
    integer, dimension(:,:), allocatable :: indx
    real(dp), dimension(ndim**2,ndim**2) :: Blocal

    if ( any( ndim /= shape(A) ) ) then
      write(*,'(/2(a/))') ' Error in tensor4_to_matrix_left_symmetric: ', &
        ' Extend of the dimensions of tensor A not equal to specified ndim.'
      stop
    end if

    indx = reshape ( [(ii,ii=1,ndim**2)], [ndim,ndim], order=[2,1] )

    ii = 0
    do i = 1, ndim
      do j = i, ndim
        ii = ii + 1
        is(ii) = indx(i,j)
      end do
    end do

    do l = 1, ndim
      do j = 1, ndim
        Blocal(indx(:,j),indx(:,l)) = A(:,j,:,l)
      end do
    end do

    B = Blocal(is,:)

  end function tensor4_to_matrix_left_symmetric


! Convert matrix to tensor4, with "left" tensor assumed symmetric

  function matrix_to_tensor4_left_symmetric ( ndim, A ) result(B)

!   dimension of space
    integer, intent(in) :: ndim

!   matrix with rows ordering according to a symmetric row-major tensor2, i.e.
!   ndim=2: 11, 12, 22
!   ndim=3: 11, 12, 13, 22, 23, 33
!   and with columns ordering according to a row-major tensor2, i.e.
!   ndim=2: 11, 12, 21, 22
!   ndim=3: 11, 12, 13, 21, 22, 23, 31, 32, 33
!   NOTE: Only the rows corresponding the upper-triangle of the "left" tensor
!         are assumed to be in the matrix and the "lower-triangle" of the
!         "left" tensor can be found by symmetry.
    real(dp), dimension(:,:), intent(in) :: A

!   tensor4 T_ijkl with i,j,k,l=1,...,ndim
!   NOTE: it is assumed that the 4th-order tensor is left-symmetric, i.e.
!         T_ijkl=T_jikl and only the upper-diagonal of the "left" tensor is
!         used in the matrix. The lower-diagonal is filled by assuming
!         left-symmetry.
    real(dp), dimension(ndim,ndim,ndim,ndim) :: B

    integer :: ii, i, j, k, l
    integer, dimension(:,:), allocatable :: indx

    if ( any( [ndim*(ndim+1)/2,ndim**2] /= shape(A) ) ) then
      write(*,'(/2(a/))') ' Error in matrix_to_tensor4_left_symmetric: ', &
        ' Shape of matrix A not consistent with specified ndim.'
      stop
    end if

    indx = reshape ( [(ii,ii=1,ndim**2)], [ndim,ndim], order=[2,1] )

!   upper-diagonal part

    ii = 0
    do i = 1, ndim
      do j = i, ndim
        ii = ii + 1
        do k = 1, ndim
          do l = 1, ndim
             B(i,j,k,l) = A(ii,indx(k,l))
          end do
        end do
      end do
    end do

!   assume left symmetric

    do i = 1, ndim
      do j = i+1, ndim
        B(j,i,:,:) = B(i,j,:,:)
      end do
    end do

  end function matrix_to_tensor4_left_symmetric

end module tensor_m
