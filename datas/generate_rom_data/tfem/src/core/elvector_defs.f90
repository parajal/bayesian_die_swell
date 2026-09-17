
! Copyright (C) 2007-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, for the elvector (vector defined per element)

module elvector_defs_m

  use kind_defs_m
  use array_defs_m

  implicit none


! type definition of elvector for a single element group

  type elvector_group_t

!   one-dimensional integer arrays (for each element)
!   i1(elem,i)%a gives array number i for element elem
    type(int_array_1d_t), allocatable, dimension(:,:) :: i1

!   one-dimensional real arrays (for each element)
!   r1(elem,i)%a gives array number i for element elem
    type(real_array_1d_t), allocatable, dimension(:,:) :: r1

!   two-dimensional real arrays (for each element)
!   r2(elem,i)%a gives array number i for element elem
    type(real_array_2d_t), allocatable, dimension(:,:) :: r2

!   three-dimensional real arrays (for each element)
!   r3(elem,i)%a gives array number i for element elem
    type(real_array_3d_t), allocatable, dimension(:,:) :: r3

  end type elvector_group_t


! type definition of elvector

  type elvector_t

!   has the vector been created?
!   Note, that created=.true. only means that the arrays in the current
!   structure and elvector_group_t (i1,r1,r2,r3) have been allocated.
!   The allocation/deallocation of data arrays in the structures itself
!   is the responsibility of the user.
    logical :: created = .false.

!   g(elgrp) gives the elvector data for element group elgrp
    type(elvector_group_t), allocatable, dimension(:) :: g

  end type elvector_t

end module elvector_defs_m
