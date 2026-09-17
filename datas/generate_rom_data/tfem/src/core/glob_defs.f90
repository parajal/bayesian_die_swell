
! Copyright (C) 2026-2026 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Define global parameters and routines. Includes kind_defs module.

module glob_defs_m

  use kind_defs_m

  implicit none

contains

! Give error message for case default in select case statement

  subroutine errormsg_case_default ( name_of_routine, name_of_select_var, &
    int_value, char_value )

!   name of calling routine and select variable
    character (len=*), intent(in) :: name_of_routine, name_of_select_var

!   integer selector
    integer, intent(in), optional :: int_value

!   character selector
    character (len=*), intent(in), optional :: char_value

    write(*,'(2(/2a))') 'Error in ', trim(name_of_routine)//':', &
       '  No case available for case selector variable ', &
       trim(name_of_select_var)
    if ( present(int_value) ) write(*,'(a,i0/)') '  Value = ', int_value
    if ( present(char_value) ) write(*,'(2a/)') '  Value = ', trim(char_value)
    stop

  end subroutine errormsg_case_default

end module glob_defs_m
