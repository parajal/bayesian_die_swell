
! Copyright (C) 2007-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Global variables for the generalized Stokes elements

module generalized_stokes_globals_m

  use stokes_globals_m

  implicit none

  save


! number of points
  integer :: nump = 0


! various allocatable arrays:

  real(dp), allocatable, dimension(:) :: etaconst, etag, detag, eta0, etainf, &
    lambda, tauyield, tauvm

  real(dp), allocatable, dimension(:,:) :: brownf


end module generalized_stokes_globals_m
