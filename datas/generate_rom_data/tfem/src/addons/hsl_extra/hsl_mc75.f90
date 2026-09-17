
! Copyright (C) 2010-2010 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Routines, for computing the condition number of the system of equations
! using MC75 from the HSL library (http://www.numerical.rl.ac.uk/hsl)

module hsl_mc75_m

  use glob_defs_m
  use sparse_m, only: sparsematrix_t
  use system_defs_m
  use set_optional_m

  implicit none


contains


! Compute the condition number of the matrix of the system of equations
! using MC75.

  subroutine condition_number_system_mc75 ( sysmatrix, cond, mem )

    type(sysmatrix_t), intent(in) :: sysmatrix

!   Condition number output. From the manual of mc75:
!   On return COND(1) holds the value of the classical condition number
!   in the infinity–norm, and COND(2) holds the value of Skeel’s condition
!   number (but see Section 2.5 for the case when A is singular). In the
!   case of an error return, COND(1) and COND(2) are set to zero.
    real(dp), dimension(2), intent(out) :: cond

!   If present: the value of la (length of A array) is mem * nz, where
!   nz is the number of non-zeros in the matrix.
!   Default = 2log n, with n=order of the matrix (number of unknowns). This
!   value is the maximum value according to the manual. So in practice mem
!   can be given a lower value to reduce memory requirements.
    real(dp), intent(in), optional :: mem


!   some testing first

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(/a/)') &
        'Error in condition_number_system_mc75: no system matrix structure.'
      stop
    end if

    if ( .not. sysmatrix%finalized ) then
      write(*,'(/a/a/)') &
        'Error in condition_number_system_mc75:', &
        '  system matrix has not been finalized.'
      stop
    end if

    if ( .not. sysmatrix%allocated_data ) then
      write(*,'(/a/a/)') &
        'Error in condition_number_system_mc75:', &
        '  data in system matrix not allocated.'
      stop
    end if

    if ( sysmatrix%symmetric ) then
      write(*,'(/a/)') &
        'Error in condition_number_system_mc75: symmetric matrix not allowed.'
      stop
    end if

    call condition_number_mc75 ( sysmatrix%Suu, cond, mem )

  end subroutine condition_number_system_mc75


! Compute the condition number of the matrix A using MC75.

  subroutine condition_number_mc75 ( A, cond, mem )

!   The matrix A in CSR format
    type(sparsematrix_t), intent(in) :: A

!   Condition number output. From the manual of mc75:
!   On return COND(1) holds the value of the classical condition number
!   in the infinity–norm, and COND(2) holds the value of Skeel’s condition
!   number (but see Section 2.5 for the case when A is singular). In the
!   case of an error return, COND(1) and COND(2) are set to zero.
    real(dp), dimension(2), intent(out) :: cond

!   If present: the value of la (length of A array) is mem * nz, where
!   nz is the number of non-zeros in the matrix.
!   Default = 2log n, with n=order of the matrix (number of unknowns). This
!   value is the maximum value according to the manual. So in practice mem
!   can be given a lower value to reduce memory requirements.
    real(dp), intent(in), optional :: mem


    integer :: n, nz, la, row, liw, lw, i
    integer, dimension(5) :: icntl, info
    real(dp) :: lmem
    integer, dimension(:), allocatable :: irn, jcn, iw
    real(dp), dimension(:), allocatable :: Amat, w



!   initialize local parameters

    icntl = 0
    info = 0

!   set dimensions of problem

    n  = A%n
    nz = A%nnz

!   set parameters

    lmem = set_optional ( variable=mem, default=log(real(n,dp))/log(2._dp) )
    la = max ( 4*nz, nint(lmem*nz) )
    liw = 19*n+7
    lw = 8*n

!   allocate memory

    allocate ( irn(la), jcn(la), iw(liw), Amat(la), w(lw) )

!   initialize

    call mc75id ( icntl )

!   fill row number array to transform CSR -> coordinate format

    do row = 1, n
      irn( A%ia(row):A%ia(row+1)-1 ) = row
    end do

!   copy arrays

    jcn(1:nz) = A%ja
    Amat(1:nz) = A%a

!   compute condition numbers

    call mc75ad ( n, nz, la, amat, irn, jcn, cond, liw, iw, lw, w, icntl, info )

    if ( info(1) < 0) then

!     Error

      write(*,'(/a,i0/)') 'Error in mc75ad, INFO(1) = ', info(1)

      write(*,'(/a)') 'INFO:'

      do i = 1, 5
         write(*,'(i5,1x,i0)') i, info(i)
      end do

      stop

    else if ( info(1) > 0) then

!     Warning

      write(*,'(/a,i0/)') 'Warning in mc75ad, INFO(1) = ', info(1)

      write(*,'(/a)') 'INFO:'

      do i = 1, 5
         write(*,'(i5,1x,i0)') i, info(i)
      end do

    end if

!   deallocate memory

    deallocate ( irn, jcn, iw, Amat, w )

  end subroutine condition_number_mc75

end module hsl_mc75_m
