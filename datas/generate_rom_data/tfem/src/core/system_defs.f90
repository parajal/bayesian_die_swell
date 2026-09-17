
! Copyright (C) 2004-2016 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, for the system of equations

module system_defs_m

  use kind_defs_m
  use sparse_m, only: sparsematrix_t

  implicit none


! type definition of system matrix
!
!  [  Suu   Sup  ]
!  [  Spu   Spp  ]
!
  type sysmatrix_t

!   the structure of sysmatrix (the arrays ia) has been initialized
    logical :: initialized_structure = .false.

!   the structure of sysmatrix (the arrays ia) has been finalized
    logical :: finalized = .false.

!   the data in the sysmatrix has been allocated
    logical :: allocated_data = .false.

!   symmetry of the matrix
    logical :: symmetric = .false.

!   renumbering of unknown degrees for minimal fill-in (e.g. metis).
!   this logical only indicates whether the arrays permu and ipermu are present.
    logical :: renumber = .false.

!   permutation arrays for minimal fill-in
!     permu(i) is the old number of new degree i
!     ipermu(i) is the new number of old degree i
    integer, allocatable, dimension(:) :: permu, ipermu

!   the square system matrix of the unknown degrees of freedom
    type(sparsematrix_t) :: Suu

!   the non-square system matrix coupling the unknown degrees of freedom with
!   the prescribed degrees of freedom.
    type(sparsematrix_t) :: Sup

!   the non-square system matrix coupling the prescribed degrees of freedom with
!   the unknown degrees of freedom.
    type(sparsematrix_t) :: Spu

!   the square system matrix of the prescribed degrees of freedom
    type(sparsematrix_t) :: Spp

  end type sysmatrix_t


! type definition of the system vector (solution, rhs)

  type sysvector_t

    integer :: probnr = 0  ! problem number. If 0 no problem number is defined.
    integer :: n = 0  ! size of the vector = number of degrees of freedom
    logical :: created = .false.  ! has the sysvector been created

!   the components sequentially
!   a partioning like the system matrix is implicitly assumed:
!
!    [ u_u ]
!    [ u_p ]
!
    real(dp), allocatable, dimension(:) :: u

  end type sysvector_t


! type definition of a vector subscript in the system vector

  type subscript_t

!   the subscript. Data can be accessed with sysvector%u(subscript%s).
    integer, allocatable, dimension(:) :: s

!   the nodes involved in the subscript that have degrees of freedom in the
!   subscript. This is only filled when fillnodes=.true. in create_subscript.
    integer, allocatable, dimension(:) :: nodes

  end type subscript_t

end module system_defs_m
