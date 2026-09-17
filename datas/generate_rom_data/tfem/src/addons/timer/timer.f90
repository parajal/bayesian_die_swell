
! Copyright (C) 2005-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! timer routines:
!   tic and toc for measuring the cpu time
!   tic_w and toc_w for measuring the wall clock time

module timer_m

  use, intrinsic :: iso_fortran_env, only: int64

  implicit none

! logical to set timer on and off
  logical :: timer = .true.

  real, private :: t_start, t_finish, t_before
  integer(int64), private :: clock_start, clock_end, clock_before, clock_rate

  save

contains


! start the timer

  subroutine tic

    if ( .not. timer ) return

    call cpu_time ( t_start )

    t_before = t_start

  end subroutine tic


! print the elapsed time

  subroutine toc ( string, delta )

!   if string is present the incremental time since the last call of toc (or
!   after a restart using tic) is printed along with the string
    character (len=*), optional, intent(in) :: string

!   if delta is present: the incremental time
    real, optional, intent(out) :: delta

    if ( .not. timer ) return

    call cpu_time ( t_finish )

    if ( present(string) ) then
      write ( *, '(2(a,f0.5),2a)' ) 'Elapsed CPU time = ', &
        t_finish - t_before, ' incremental and ', t_finish - t_start, &
        ' total seconds in ', string
    else
      write ( *, '(a,f0.5,a)' ) 'Elapsed CPU time = ', t_finish - t_start, &
        ' seconds'
    end if

    if ( present (delta) ) delta = t_finish - t_before

    t_before = t_finish

  end subroutine toc


! start the wall-clock timer

  subroutine tic_w

    if ( .not. timer ) return

    call system_clock ( count_rate = clock_rate ) ! get rate (ticks/sec)
    call system_clock ( count = clock_start )     ! start ticking

    clock_before = clock_start

  end subroutine tic_w


! print the elapsed wall-clock time

  subroutine toc_w ( string, delta )

!   if string is present the incremental time since the last call of toc_w (or
!   after a restart using tic_w) is printed along with the string
    character (len=*), optional, intent(in) :: string

!   if delta is present: the incremental time
    real, optional, intent(out) :: delta

    real :: elapsed_time_total, elapsed_time_incr

    if ( .not. timer ) return

    call system_clock ( count = clock_end )       ! get current ticks

    elapsed_time_incr = real(clock_end-clock_before)/real(clock_rate)
    elapsed_time_total  = real(clock_end-clock_start)/real(clock_rate)

    if ( present(string) ) then
      write ( *, '(2(a,f0.5),2a)' ) 'Elapsed wall-clock time = ', &
        elapsed_time_incr, ' incremental and ', elapsed_time_total, &
        ' total seconds in ', string
    else
      write ( *, '(a,f0.5,a)' ) 'Elapsed wall-clock time = ', &
                elapsed_time_total, ' seconds'
    end if

    if ( present (delta) ) delta = elapsed_time_incr

    call system_clock ( count = clock_before )    ! get current ticks

  end subroutine toc_w

end module timer_m
