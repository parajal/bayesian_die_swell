
! Copyright (C) 2004-2008 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, for the vectors

module vector_defs_m

  use kind_defs_m

  implicit none


! type definition of vector

  type vector_t

    integer :: probnr = 0  ! problem number. If 0 no problem number is defined.
    integer :: n = 0    ! size of the vector = number of components
    integer :: vec = 0  ! vector number
    logical :: elementwise = .false.  ! store elementwise
    logical :: created = .false.  ! has the vector been created

!   the components
    real(dp), allocatable, dimension(:) :: u

!   the weighting factors
    real(dp), allocatable, dimension(:) :: w

  end type vector_t


! type definition of a vector subscript in a vector

  type subscriptvec_t

!   the subscript. Data can be accessed with vector%u(subscript%s).
    integer, allocatable, dimension(:) :: s

!   the nodes involved in the subscript that have degrees of freedom in the
!   subscript. This is only filled when fillnodes=.true. in create_subscript.
    integer, allocatable, dimension(:) :: nodes

  end type subscriptvec_t


end module vector_defs_m
