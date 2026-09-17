
! Copyright (C) 2016-2016 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Global variables for the convection_diffusion_supg elements

module convection_diffusion_supg_globals_m

  use kind_defs_m
  use stokes_globals_m

  implicit none

  save

! number of c degrees of freedom
  integer :: ndfq = 0

! shapefunc types
  type(shapefunc_t) :: shapefuncq

! SUPG parameters
  integer :: htypeq, Uscalingq

! various allocatable arrays:
  real(dp), allocatable, dimension(:) :: htimesU, q_ng, q_nm1g, q_n, &
    q_nm1, Peh, f_g, tauq, hoverU

  real(dp), allocatable, dimension(:,:) :: Smat, Mmat, convmat, u_n, &
    ungradchi, chi, mq_n, mq_nm1, mq_ng, mq_nm1g, mf_g

  real(dp), allocatable, dimension(:,:,:) :: dchi, dchidx

end module convection_diffusion_supg_globals_m
