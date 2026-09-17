
! Copyright (C) 2006-2012 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Global variables for the compressible fluid elements

module compressible_fluid_globals_m

  use kind_defs_m
  use stokes_globals_m

  implicit none

  save


! various allocatable arrays:

  real(dp), allocatable, dimension(:) :: p, png, pnm1g, png_iter, &
    unp1gradpiter, divgiter, lagmul, pg_iter

  real(dp), allocatable, dimension(:,:) :: ungradpsi, uvecn, uvecnm1, &
    uvecmeshn, uvecnp1, uvecmeshnp1, uvechat, uvecmeshhat, uvec_iter, &
    unp1gradpsi, gradpiter

end module compressible_fluid_globals_m
