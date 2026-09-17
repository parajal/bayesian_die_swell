
! Copyright (C) 2004-2010 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines arrays encapsulated in a type

module array_defs_m

  use kind_defs_m

  implicit none

  type logical_array_1d_t
    logical, allocatable, dimension(:) :: a
  end type logical_array_1d_t

  type logical_array_2d_t
    logical, allocatable, dimension(:,:) :: a
  end type logical_array_2d_t

  type logical_array_3d_t
    logical, allocatable, dimension(:,:,:) :: a
  end type logical_array_3d_t

  type int_array_1d_t
    integer, allocatable, dimension(:) :: a
  end type int_array_1d_t

  type int_array_2d_t
    integer, allocatable, dimension(:,:) :: a
  end type int_array_2d_t

  type int_array_3d_t
    integer, allocatable, dimension(:,:,:) :: a
  end type int_array_3d_t

  type real_array_1d_t
    real(dp), allocatable, dimension(:) :: a
  end type real_array_1d_t

  type real_array_2d_t
    real(dp), allocatable, dimension(:,:) :: a
  end type real_array_2d_t

  type real_array_3d_t
    real(dp), allocatable, dimension(:,:,:) :: a
  end type real_array_3d_t

end module array_defs_m

