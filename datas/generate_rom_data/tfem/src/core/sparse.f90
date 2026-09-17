
! Copyright (C) 2004-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, routines, for sparse matrices in CSR format

module sparse_m

  use kind_defs_m
  use misc_m

  implicit none


! Type definition of sparse matrix (Compressed Sparse Row (CSR) format)

  type sparsematrix_t

    integer :: n    = 0        ! row dimension of the matrix
    integer :: m    = 0        ! column dimension of the matrix
    integer :: nnz  = 0        ! number of non-zeros

!   Array containing the nnz non-zeros
    real(dp), allocatable, dimension(:) :: a

!   Array containing the nnz column numbers of the non-zeros
    integer, allocatable, dimension(:) :: ja

!   Array of length n+1 containing the pointers to beginning of each row in the
!   arrays a and ja. Thus ia(i) is position in arrays a and ja where row i
!   starts. Last pointer: ia(n+1) = ia(1) + nnz, which points to a fictitious
!   row n+1.
    integer, allocatable, dimension(:) :: ia

  end type sparsematrix_t

! Interface for generic delete subroutine

  interface delete
    module procedure delete_sparsematrix
  end interface delete

! Interface for generic diagonal subroutine

  interface diagonal
    module procedure diagonal_sparsematrix
  end interface diagonal

! Interface for generic extract_submat subroutine

  interface extract_submat
    module procedure extract_submat_sparsematrix
  end interface extract_submat

contains


! Delete matrix

  subroutine delete_sparsematrix ( matrix )

    type(sparsematrix_t), intent(out) :: matrix

  end subroutine delete_sparsematrix


! Transpose of sparse matrix A

  function stranspose ( A ) result(B)

    type ( sparsematrix_t ), intent(in) :: A

    type ( sparsematrix_t ) :: B

    integer :: i, j, k, pos

!   set dimensions and allocate sparse matrix B
    B%nnz = A%nnz
    B%n = A%m
    B%m = A%n
    allocate ( B%ia(B%n+1), B%ja(B%nnz), B%a(B%nnz) )

!   find number of non-zero columns in each row of B
    B%ia = 0
    do i = 1, A%n
      do k = A%ia(i), A%ia(i+1) - 1
        j = A%ja(k) + 1
        B%ia(j) = B%ia(j) + 1
      end do
    end do

!   accumulate B%ia
    B%ia(1) = 1
    do i = 2, B%n+1
      B%ia(i) = B%ia(i) + B%ia(i-1)
    end do

!   copy the data and increment B%ia
    do i = 1, A%n
      do k = A%ia(i), A%ia(i+1) - 1
        j = A%ja(k)
        pos = B%ia(j)
        B%ja(pos) = i
        B%a(pos) = A%a(k)
        B%ia(j) = pos + 1
      end do
    end do

!   shift B%ia one up
    do i = B%n, 1, -1
      B%ia(i+1) = B%ia(i)
    end do
    B%ia(1) = 1

  end function stranspose


! Sparse matrix-vector multiplication.
! NOTE: assumes all non-zeros are present!

  function smatvec ( A, x ) result(b)

    type ( sparsematrix_t ), intent(in) :: A
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(A%n) :: b

    integer :: i, j
    real(dp) :: tmp

    do i = 1, A%n
      tmp = 0
      do j = A%ia(i), A%ia(i+1) - 1
        tmp = tmp + A%a(j) * x(A%ja(j))
      end do
      b(i) = tmp
    end do

  end function smatvec


! Sparse matrix-vector multiplication A^T x.
! NOTE: assumes all non-zeros are present!

  function stmatvec ( A, x ) result(b)

    type ( sparsematrix_t ), intent(in) :: A
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(A%m) :: b

    integer :: i, j, k

    b = 0

    do i = 1, A%n
      do k = A%ia(i), A%ia(i+1) - 1
        j = A%ja(k)
        b(j) = b(j) + A%a(k) * x(i)
      end do
    end do

  end function stmatvec


! Multiplication of sparse matrix A and B
! The result is a sparse matrix C = A B
! NOTE: one third of the CPU time is due to the "dry run" to obtain the
! size of the matrix C.

  function smatmul ( A, B ) result(C)

    type ( sparsematrix_t ), intent(in) :: A, B

    type ( sparsematrix_t ) :: C

    integer :: i, j, k, m, p, nc, jj, pold
!   if workc(m)/=0 it stores the position of the non zero element for row i
!   and column m (it is restarted every new row i).
    integer, allocatable, dimension(:) :: workc
!   wcol(1:nc) stores the column numbers of row i that have a non zero element
!   (it is restarted every new row i).
    integer, allocatable, dimension(:) :: wcol

    allocate ( workc(B%m), wcol(B%m) )

    workc = 0

!   dry run to obtain the number of non-zero elements

    p = 0
    nc = 0

    do i = 1, A%n  ! row i
      do j = A%ia(i), A%ia(i+1) - 1
        k = A%ja(j) ! column number k of matrix A
        do jj = B%ia(k), B%ia(k+1) - 1  ! the k^th row of matrix B
          m = B%ja(jj) ! column number m of matrix B
          if ( workc(m) == 0 ) then
!           new element found
            p = p + 1
            workc(m) = p
            nc = nc + 1
            wcol(nc) = m
          end if
        end do
      end do
!     set workc back to zero for row i
      workc(wcol(1:nc)) = 0
      nc = 0
    end do

!   set number of non zeros and initialize matrix C

    C%nnz = p
    C%n = A%n
    C%m = B%m
    allocate ( C%ia(C%n+1), C%ja(C%nnz), C%a(C%nnz) )

!   fill matrix C = A B

    p = 0

    C%ia(1) = 1

    do i = 1, A%n  ! row i
      do j = A%ia(i), A%ia(i+1) - 1
        k = A%ja(j) ! column number k of matrix A
        do jj = B%ia(k), B%ia(k+1) - 1 ! the k^th row of matrix B
          m = B%ja(jj)  ! column number m of matrix B
          pold = workc(m)
          if ( pold == 0 ) then
!           new element found
            p = p + 1
            workc(m) = p
            C%ja(p) = m
            C%a(p) = A%a(j) * B%a(jj)
          else
            C%a(pold) = C%a(pold) + A%a(j) * B%a(jj)
          end if
        end do
      end do
      C%ia(i+1) = p + 1
!     set workc back to zero for row i
      workc(C%ja(C%ia(i):p)) = 0
    end do

    deallocate ( workc, wcol )

  end function smatmul


! Multiplication of sparse matrix A and full matrix B
! The result is a full matrix C = A B

  function sfmatmul ( A, B ) result(C)

    type ( sparsematrix_t ), intent(in) :: A
    real(dp), dimension(:,:), intent(in) :: B

    real(dp), dimension(A%n,size(B,2)) :: C

    integer :: i, j, k, m

!   fill matrix C = A B

    C = 0

    do i = 1, A%n  ! row i of matrix A
      do m = A%ia(i), A%ia(i+1) - 1
        k = A%ja(m) ! column number k of matrix A
        do j = 1, size(B,2)
          C(i,j) = C(i,j) + A%a(m) * B(k,j)
        end do
      end do
    end do

  end function sfmatmul


! Multiplication of full matrix A and sparse matrix B
! The result is a full matrix C = A B

  function fsmatmul ( A, B ) result(C)

    real(dp), dimension(:,:), intent(in) :: A
    type ( sparsematrix_t ), intent(in) :: B

    real(dp), dimension(size(A,1),B%m) :: C

    integer :: i, j, k, m

!   fill matrix C = A B

    C = 0

    do k = 1, B%n  ! row k of matrix B
      do m = B%ia(k), B%ia(k+1) - 1
        j = B%ja(m) ! column number j of matrix B
        do i = 1, size(A,1)
          C(i,j) = C(i,j) + A(i,k) * B%a(m)
        end do
      end do
    end do

  end function fsmatmul


! Multiplication of the transpose of sparse matrix A and full matrix B
! The result is a full matrix C = A^T B

  function stfmatmul ( A, B ) result(C)

    type ( sparsematrix_t ), intent(in) :: A
    real(dp), dimension(:,:), intent(in) :: B

    real(dp), dimension(A%m,size(B,2)) :: C

    integer :: i, j, k, m

!   fill matrix C = A^T B

    C = 0

    do k = 1, A%n  ! column k of matrix A^T
      do m = A%ia(k), A%ia(k+1) - 1
        i = A%ja(m) ! row number i of matrix A^T
        do j = 1, size(B,2)
          C(i,j) = C(i,j) + A%a(m) * B(k,j)
        end do
      end do
    end do

  end function stfmatmul


! Multiplication of full matrix A and the transpose of sparse matrix B
! The result is a full matrix C = A B^T

  function fstmatmul ( A, B ) result(C)

    real(dp), dimension(:,:), intent(in) :: A
    type ( sparsematrix_t ), intent(in) :: B

    real(dp), dimension(size(A,1),B%n) :: C

    integer :: i, j, k, m

!   fill matrix C = A B

    C = 0

    do j = 1, B%n  ! column j of matrix B^T
      do m = B%ia(j), B%ia(j+1) - 1
        k = B%ja(m) ! row number k of matrix B
        do i = 1, size(A,1)
          C(i,j) = C(i,j) + A(i,k) * B%a(m)
        end do
      end do
    end do

  end function fstmatmul


! Sort rows of matrix for increasing column number

  subroutine sort_sparsematrix ( A )

    type(sparsematrix_t), intent(inout) :: A

    integer :: i, j1, j2
    integer, allocatable, dimension(:) :: order

    allocate ( order(A%n) )

    do i = 1, A%n
      j1 = A%ia(i)
      j2 = A%ia(i+1) - 1
      call sort ( A%ja(j1:j2), order )
      A%a(j1:j2) = A%a(j1-1+order(1:j2-j1+1))
    end do

    deallocate ( order )

  end subroutine sort_sparsematrix


! Get (main) diagonal of matrix

  function diagonal_sparsematrix ( A, error ) result(d)

    type ( sparsematrix_t ), intent(in) :: A

!   error=.true. if not all diagonal elements are present in the matrix.
!                The corresponding diagonal element is set to zero.
    logical, intent(out) :: error

    real(dp), dimension(A%n) :: d

    integer :: i, j

    error = .false.

    row: do i = 1, A%n
      do j = A%ia(i), A%ia(i+1) - 1
        if ( A%ja(j) == i ) then
          d(i) = A%a(j)
          cycle row
        end if
      end do
      d(i) = 0
      error = .true.
    end do row

  end function diagonal_sparsematrix


! Extract sub-matrix from sparse matrix
! The result is a sparse matrix S

  function extract_submat_sparsematrix ( A, rs, cs ) result(S)

    type ( sparsematrix_t ), intent(in) :: A

!   index (vector subscript) arrays for rows and columns
!      S = A( rs, cs )
!   NOTE: duplicate indices are not allowed.
    integer, dimension(:), intent(in) :: rs, cs

    type ( sparsematrix_t ) :: S

    integer :: i, j, k, m, nr, nc
    logical, allocatable, dimension(:) :: workr
    integer, allocatable, dimension(:) :: workc

    allocate ( workr(A%n), workc(A%m) )

!   test

    if ( any( rs > A%n ) .or. any( rs < 1 ) .or. &
         any( cs > A%m ) .or. any( cs < 1 ) ) then
      write(*,'(/2(a/))') ' Error in extract_submat_sparsematrix: ', &
        ' Index arrays rs and/or cs are out of range '
      stop
    end if

!   use work array workr to find nr

    workr = .false.
    workr(rs) = .true.
    nr = count(workr)

!   set work array with new column numbers for S

    workc = 0
    nc = 0
    do k = 1, size(cs)
      j = cs(k)
      if ( workc(j) > 0 ) cycle ! already counted
      nc = nc + 1
      workc(j) = nc
    end do

!   test duplicate entries

    if ( nr /= size(rs) .or. nc /= size(cs) ) then
      write(*,'(/2(a/))') ' Error in extract_submat_sparsematrix: ', &
        ' Index arrays rs and/or cs have duplicate entries '
      stop
    end if

!   Initialize matrix S

    S%n = nr
    S%m = nc
    S%nnz = 0
    allocate(S%ia(nr+1))
    S%ia = 0

!   count matrix entries

    do k = 1, nr
      i = rs(k)
      do j = A%ia(i), A%ia(i+1) - 1
        if ( workc(A%ja(j)) > 0 ) then
          S%ia(k+1) = S%ia(k+1) + 1  ! column found
        end if
      end do
    end do

!   accumulate ia

    S%ia(1) = 1
    do k = 1, nr
      S%ia(k+1) = S%ia(k+1) + S%ia(k)
    end do
    S%nnz = S%ia(nr+1) - 1

!   allocate data of matrix S

    allocate( S%a(S%nnz), S%ja(S%nnz) )

!   fill matrix entries

    m = 0
    do k = 1, nr
      i = rs(k)
      do j = A%ia(i), A%ia(i+1) - 1
        if ( workc(A%ja(j)) > 0 ) then ! column found
          m = m + 1
          S%ja(m) = workc(A%ja(j))
          S%a(m) = A%a(j)
        end if
      end do
    end do

    deallocate ( workr, workc )

  end function extract_submat_sparsematrix


! clear rows of matrix

  subroutine clear_rows_sparsematrix ( rows, A )

!   index (vector subscript) arrays for rows
    integer, dimension(:), intent(in) :: rows
    type(sparsematrix_t), intent(inout) :: A

    integer :: i, k, j1, j2

    do i = 1, size(rows)
      k = rows(i)
      if ( k < 1 .or. k > A%n ) cycle
      j1 = A%ia(k)
      j2 = A%ia(k+1) - 1
      A%a(j1:j2) = 0
    end do

  end subroutine clear_rows_sparsematrix


! reorder (permute) rows of matrix

  subroutine reorder_rows_sparsematrix ( perm, A, B )

!   permutation array to reorder the rows of the sparse matrix A
!   perm(i) is the new row number of row i
    integer, dimension(:), intent(in) :: perm

    type(sparsematrix_t), intent(inout) :: A

!   if present the resulting matrix is written to B and A is kept.
!   Otherwise matrix A will be overwritten.
    type(sparsematrix_t), optional, intent(out) :: B

    integer, allocatable, dimension(:) :: wia, wja
    real(dp), allocatable, dimension(:) :: wa

    integer :: i, j, k, m

!   allocate work storage

    allocate( wia(A%n+1), wja(A%nnz), wa(A%nnz) )

!   collect number of non-zeros in each new row

    do i = 1, A%n
      k = perm(i) ! new row
      wia(k+1) = A%ia(i+1) - A%ia(i)
    end do

!   accumulate

    wia(1) = 1
    do k = 1, A%n
      wia(k+1) = wia(k+1) + wia(k)
    end do

!   copy column numbers and matrix elements

    do i = 1, A%n
      m = wia(perm(i))
      do j = A%ia(i), A%ia(i+1) - 1
        wja(m) = A%ja(j)
        wa(m) = A%a(j)
        m = m + 1
      end do
    end do

    if ( present(B) ) then

      B%n = A%n
      B%m = A%n
      B%nnz = A%nnz

      call move_alloc ( from=wia, to=B%ia )
      call move_alloc ( from=wja, to=B%ja )
      call move_alloc ( from=wa, to=B%a )

    else

      call move_alloc ( from=wia, to=A%ia )
      call move_alloc ( from=wja, to=A%ja )
      call move_alloc ( from=wa, to=A%a )

    end if

  end subroutine reorder_rows_sparsematrix


! reorder (permute) columns of matrix

  subroutine reorder_columns_sparsematrix ( perm, A, B )

!   permutation array to reorder the columns of the sparse matrix A
!   perm(j) is the new column number of column j
    integer, dimension(:), intent(in) :: perm

    type(sparsematrix_t), intent(inout) :: A

!   if present the resulting matrix is written to B and A is kept.
!   Otherwise matrix A will be overwritten.
    type(sparsematrix_t), optional, intent(out) :: B

!   new column numbers

    if ( present(B) ) then

      B%n = A%n
      B%m = A%n
      B%nnz = A%nnz
      B%ia = A%ia
      B%a = A%a

      B%ja = perm(A%ja)

    else

      A%ja = perm(A%ja)

    end if

  end subroutine reorder_columns_sparsematrix


! Add sparse matrix, A + B -> A in place
! Row and column dimensions need to be the same.
! Matrix A needs to have already sufficient space reserved for adding B,
! denoted by column nrs set to 0 in A%ja at the end of a row.

  subroutine add_sparsematrix ( A, B, symmetric )

    type(sparsematrix_t), intent(inout) :: A
    type(sparsematrix_t), intent(in) :: B

!   if present and .true. the lower triangle is ignored
    logical, optional, intent(in) :: symmetric

    logical :: sym
    integer :: i, pf, pa, pb, cola, colb
    real(dp), allocatable, dimension(:) :: w

    if ( any ( [ A%n, A%m ] /= [ B%n, B%m ] ) ) then
      write(*,'(/2(a/))') ' Error in add_sparsematrix: ', &
        '    Matrices A and B do not have the same dimensions.'
      stop
    end if

    if ( present(symmetric) ) then
      sym = symmetric
    else
      sym = .false.
    end if

    if ( sym .and. A%n /= A%m ) then
      write(*,'(/2(a/))') ' Error in add_sparsematrix: ', &
        '    For symmetric=.true. matrices A and B need to square.'
      stop
    end if

    allocate ( w(A%m) )  ! work array

    w = 0

    do i = 1, A%n

!     fill w for non-zeros already in matrix A for row i

!     pf is position in row i of A for new data; start with none available
      pf = A%ia(i+1)
      do pa = A%ia(i), A%ia(i+1) - 1
        cola = A%ja(pa)
        if ( cola == 0 ) then
!         matrix row not completely filled
          pf = pa
          exit
        end if
        w(cola) = pa ! put in w position of non-zero entry
      end do

!     add row i of B
      do pb = B%ia(i), B%ia(i+1) - 1
        colb = B%ja(pb)
        if ( sym .and. i > colb ) cycle ! ignore lower triangle for symmetric
        pa = w(colb)
        if ( pa == 0 .and. pf /= A%ia(i+1) ) then
!         put new non-zero in matrix
          A%a(pf) = B%a(pb)
          A%ja(pf) = colb
          pf = pf + 1
        else if ( pa /= 0 ) then
!         add to existing non-zero
          A%a(pa) = A%a(pa) + B%a(pb)
        else if ( pa == 0 .and. pf == A%ia(i+1) ) then
!         new non-zero but no free space
          write(*,'(/a,i0/a/)') &
          'Error add_sparsematrix: no more space in row ', i, &
          '   of matrix A to add row of matrix B'
          stop
        end if
      end do

!     put w back to zero

      do pa = A%ia(i), A%ia(i+1) - 1
        cola = A%ja(pa)
        if ( cola == 0 ) exit
        w(cola) = 0
      end do

    end do

  end subroutine add_sparsematrix

end module sparse_m
