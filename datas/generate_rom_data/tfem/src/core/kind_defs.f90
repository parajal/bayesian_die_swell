
! Copyright (C) 2004-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Define kinds for integers and reals in one place

module kind_defs_m

  implicit none

  integer, parameter :: si=selected_int_kind(r=9)  ! integer*4
  integer, parameter :: di=selected_int_kind(r=18) ! integer*8

  integer, parameter :: sp=selected_real_kind(P=6)  ! single precision
  integer, parameter :: dp=selected_real_kind(P=15) ! double precision

end module kind_defs_m
