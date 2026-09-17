
! Copyright (C) 2022-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

!
! Interfaces to the UMFPACK C Libary.
!

module umfpack_interfaces_m

  use iso_c_binding
  use iso_fortran_env

  implicit none

  interface

!   standard integer interfaces (32bit)

    subroutine umfpack_di_defaults ( Control ) &
      bind(c, name='umfpack_di_defaults')
      import :: c_double
      real(c_double), dimension(*) :: Control
    end subroutine umfpack_di_defaults

    integer(c_int) function umfpack_di_symbolic ( n_row, n_col, Ap, Ai, Ax, &
      Symbolic, Control, Info ) bind(c, name='umfpack_di_symbolic')
      import :: c_int, c_ptr, c_double
      integer(c_int), value :: n_row, n_col
      integer(c_int), dimension(*), intent(in) :: Ap, Ai
      real(c_double), dimension(*), intent(in) :: Ax
      type(c_ptr) :: Symbolic
      real(c_double), dimension(*) :: Control, Info
    end function umfpack_di_symbolic

    integer(c_int) function umfpack_di_numeric ( Ap, Ai, Ax, Symbolic, &
      Numeric, Control, Info ) bind(c, name='umfpack_di_numeric')
      import :: c_int, c_ptr, c_double
      integer(c_int), dimension(*), intent(in) :: Ap, Ai
      real(c_double), dimension(*), intent(in) :: Ax
      type(c_ptr), intent(in), value :: Symbolic
      type(c_ptr) :: Numeric
      real(c_double), dimension(*) :: Control, Info
    end function umfpack_di_numeric

    integer(c_int) function umfpack_di_solve ( sys, Ap, Ai, Ax, X, B, &
      Numeric, Control, Info ) bind(c, name='umfpack_di_solve')
      import :: c_int, c_ptr, c_double
      integer(c_int), value :: sys
      integer(c_int), dimension(*), intent(in) :: Ap, Ai
      real(c_double), dimension(*), intent(in) :: Ax
      real(c_double), dimension(*) :: X
      real(c_double), dimension(*), intent(in) :: B
      type(c_ptr), intent(in), value :: Numeric
      real(c_double), dimension(*) :: Control, Info
    end function umfpack_di_solve

    subroutine umfpack_di_free_symbolic ( Symbolic ) &
      bind(c, name='umfpack_di_free_symbolic')
      import :: c_ptr
      type(c_ptr) :: Symbolic
    end subroutine umfpack_di_free_symbolic

    subroutine umfpack_di_free_numeric ( Numeric ) &
      bind(c, name='umfpack_di_free_numeric')
      import :: c_ptr
      type(c_ptr) :: Numeric
    end subroutine umfpack_di_free_numeric

    subroutine umfpack_di_report_status ( Control, status ) &
      bind(c, name='umfpack_di_report_status')
      import c_int, c_double
      real(c_double), dimension(*) :: Control
      integer(c_int), value :: status
    end subroutine umfpack_di_report_status

    subroutine umfpack_di_report_control ( Control ) &
      bind(c, name='umfpack_di_report_control')
      import c_double
      real(c_double), dimension(*) :: Control
    end subroutine umfpack_di_report_control

    subroutine umfpack_di_report_info ( Control, Info ) &
      bind(c, name='umfpack_di_report_info')
      import c_double
      real(c_double), dimension(*) :: Control, Info
    end subroutine umfpack_di_report_info

!   long integer interfaces (64bit)

    subroutine umfpack_dl_defaults ( Control ) &
      bind(c, name='umfpack_dl_defaults')
      import :: c_double
      real(c_double), dimension(*) :: Control
    end subroutine umfpack_dl_defaults

    integer(c_long) function umfpack_dl_symbolic ( n_row, n_col, Ap, Ai, Ax, &
      Symbolic, Control, Info ) bind(c, name='umfpack_dl_symbolic')
      import :: c_long, c_ptr, c_double
      integer(c_long), value :: n_row, n_col
      integer(c_long), dimension(*), intent(in) :: Ap, Ai
      real(c_double), dimension(*), intent(in) :: Ax
      type(c_ptr) :: Symbolic
      real(c_double), dimension(*) :: Control, Info
    end function umfpack_dl_symbolic

    integer(c_long) function umfpack_dl_numeric ( Ap, Ai, Ax, Symbolic, &
      Numeric, Control, Info ) bind(c, name='umfpack_dl_numeric')
      import :: c_long, c_ptr, c_double
      integer(c_long), dimension(*), intent(in) :: Ap, Ai
      real(c_double), dimension(*), intent(in) :: Ax
      type(c_ptr), intent(in), value :: Symbolic
      type(c_ptr) :: Numeric
      real(c_double), dimension(*) :: Control, Info
    end function umfpack_dl_numeric

    integer(c_long) function umfpack_dl_solve ( sys, Ap, Ai, Ax, X, B, &
      Numeric, Control, Info ) bind(c, name='umfpack_dl_solve')
      import :: c_long, c_ptr, c_double
      integer(c_long), value :: sys
      integer(c_long), dimension(*), intent(in) :: Ap, Ai
      real(c_double), dimension(*), intent(in) :: Ax
      real(c_double), dimension(*) :: X
      real(c_double), dimension(*), intent(in) :: B
      type(c_ptr), intent(in), value :: Numeric
      real(c_double), dimension(*) :: Control, Info
    end function umfpack_dl_solve

    subroutine umfpack_dl_free_symbolic ( Symbolic ) &
      bind(c, name='umfpack_dl_free_symbolic')
      import :: c_ptr
      type(c_ptr) :: Symbolic
    end subroutine umfpack_dl_free_symbolic

    subroutine umfpack_dl_free_numeric ( Numeric ) &
      bind(c, name='umfpack_dl_free_numeric')
      import :: c_ptr
      type(c_ptr) :: Numeric
    end subroutine umfpack_dl_free_numeric

    subroutine umfpack_dl_report_status ( Control, status ) &
      bind(c, name='umfpack_dl_report_status')
      import c_long, c_double
      real(c_double), dimension(*) :: Control
      integer(c_long), value :: status
    end subroutine umfpack_dl_report_status

    subroutine umfpack_dl_report_control ( Control ) &
      bind(c, name='umfpack_dl_report_control')
      import c_double
      real(c_double), dimension(*) :: Control
    end subroutine umfpack_dl_report_control

    subroutine umfpack_dl_report_info ( Control, Info ) &
      bind(c, name='umfpack_dl_report_info')
      import c_double
      real(c_double), dimension(*) :: Control, Info
    end subroutine umfpack_dl_report_info

  end interface

end module umfpack_interfaces_m
