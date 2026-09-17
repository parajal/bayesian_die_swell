
! Copyright (C) 2007-2014 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


!
! set value routines for optional parameters
!

module set_optional_m

  use kind_defs_m

  implicit none

  save

  interface set_optional
    module procedure set_optional_l, set_optional_c, &
                     set_optional_i, set_optional_r
  end interface set_optional

contains


! set value for logical optional parameter

  function set_optional_l ( default, variable )

    logical, intent(in) :: default
    logical, intent(in), optional :: variable
    logical :: set_optional_l

    if ( present(variable) ) then
      set_optional_l = variable
    else
      set_optional_l = default
    end if

  end function set_optional_l


! set value for character optional parameter

  function set_optional_c ( default, variable )

    character (len=*), intent(in) :: default
    character (len=*), intent(in), optional :: variable
    character (len=:), allocatable :: set_optional_c

    if ( present(variable) ) then
      set_optional_c = variable
    else
      set_optional_c = default
    end if

  end function set_optional_c


! set value for integer optional parameter

  function set_optional_i ( default, variable )

    integer, intent(in) :: default
    integer, intent(in), optional :: variable
    integer :: set_optional_i

    if ( present(variable) ) then
      set_optional_i = variable
    else
      set_optional_i = default
    end if

  end function set_optional_i


! set value for real optional parameter

  function set_optional_r ( default, variable )

    real(dp), intent(in) :: default
    real(dp), intent(in), optional :: variable
    real(dp) :: set_optional_r

    if ( present(variable) ) then
      set_optional_r = variable
    else
      set_optional_r = default
    end if

  end function set_optional_r


end module set_optional_m

