
! Copyright (C) 2012-2012 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Contains various routines for stream I/O

module streamio_m

  implicit none

contains


! read line from file using stream I/O

  subroutine read_line ( unit, line, nc )

    integer, intent(in) :: unit
    character(len=*), intent(out) :: line
    integer, intent(out), optional :: nc ! number of characters read

    character(len=1), parameter ::  NL=char(10) !  Unix end of line
    character(len=1) :: c
    integer :: i, ios, l

    l = len(line)

    i = 0

    do

      read ( unit=unit, iostat=ios ) c

      if ( ios /= 0 ) then
        write(*,'(/a/)') 'Error in read_line: cannot find EOL.'
        stop
      end if

      if ( c == NL ) exit

      i = i + 1

      if ( i > l ) then
        write(*,'(/a,i0/)') 'Error in read_line: line length beyond ', l
        stop
      end if

      line(i:i) = c

    end do

    line(i+1:) = ''

    if ( present(nc) ) nc = i

  end subroutine read_line


! write line to file using stream I/O

  subroutine write_line ( unit, line )

    integer, intent(in) :: unit
    character(len=*), intent(in) :: line

    character(len=1), parameter ::  NL=char(10) !  Unix end of line

    write ( unit=unit ) line, NL

  end subroutine write_line


end module streamio_m
