
! Copyright (C) 2005-2006 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Global variables for the viscoelastic elements

module bcf_globals_m

  use kind_defs_m
  use shapefunc_m, only: shapefunc_t
  use viscoelastic_globals_m
  use stochastic_models_defs_m

  implicit none

  save

! use log transformation? (not yet used)
  logical :: logq = .false.

! interpolation number (shape function) for the Q-vector
  integer :: intpolq = 0

! number of Q-vector degrees of freedom with respect to one component only
  integer :: ndfq = 0

! number of degrees of freedom of Q for one mode (one point only)
  integer :: ncompq = 0

! number of Q fields
  integer :: nfield = 0


! various allocatable arrays:

  integer, allocatable, dimension(:) :: posq, posqs

  real(dp), allocatable, dimension(:) :: qv

  real(dp), allocatable, dimension(:,:) :: brownf, tauvecq

  real(dp), allocatable, dimension(:,:,:) :: qn, qnside, ungradqn, qnp1, qg, &
    qnp1tr

  real(dp), allocatable, dimension(:,:,:,:) :: qnjump


  type(stmodel_t) :: stmodel
  type(stnumpar_t) :: stnumpar

! shapefunc types

  type(shapefunc_t) :: shapefuncq

end module bcf_globals_m
