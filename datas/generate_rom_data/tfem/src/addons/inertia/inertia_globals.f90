
! Copyright (C) 2007-2011 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Global variables for the inertia elements

module inertia_globals_m

  use kind_defs_m
  use stokes_globals_m

  implicit none

  save


! various allocatable arrays:

  real(dp), allocatable, dimension(:) :: uvec

  real(dp), allocatable, dimension(:,:) :: ungradphi, uvecn, uvecmeshn, &
    uhatgradphi, uvecnm1, uvecmeshnm1, ungradun

  real(dp), allocatable, dimension(:,:,:) :: dudx

end module inertia_globals_m
