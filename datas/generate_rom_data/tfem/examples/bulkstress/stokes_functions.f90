
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

  use kind_defs_m

  implicit none

contains

  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    select case(nr)
    case(1)
      func = x(2)
    case(2)
      func = 0.0
    case default
      write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
      stop
   end select

  end function func

end module stokes_functions_m
