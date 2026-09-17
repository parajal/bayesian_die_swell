
! Copyright (C) 2005-2006 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


module stokes_functions_m

  use math_defs_m

  implicit none

contains

  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    write(*,*) 'function func has not been defined'
    func = 0
    stop

  end function func

  function vfunc ( n, nr, x )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc
    real(dp) :: alpha, period

    vfunc = 0
    alpha = .4_dp
    period = .5_dp
    if (nr == 1) then
      vfunc(1) = - 8._dp*pi**3*alpha*sin(2*pi*x(2))/period
      vfunc(2) =   0._dp
    end if
    if (nr == 2) then
      vfunc(1) =   0._dp
      vfunc(2) =   8*pi**3*alpha*sin(2*pi*x(1))/period
    end if
!   write(*,*) 'function vfunc is defined', n,nr, x, vfunc

  end function vfunc

end module stokes_functions_m
