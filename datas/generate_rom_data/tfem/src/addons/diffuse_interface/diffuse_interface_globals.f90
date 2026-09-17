
! Copyright (C) 2006-2014 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Global variables for the diffuse-interface elements

module diffuse_interface_globals_m

  use kind_defs_m
  use stokes_globals_m

  implicit none

  save


! various allocatable arrays:

  real(dp), allocatable, dimension(:) :: citer, citer_g, c_ng, c_nm1g, c_n, &
    c_nm1, mu, mu_g, c2iter, c2iter_g, c2_ng, c2_n, mu2, mu2_g, evec, &
    c_i, ungradc, gradc_sq

  real(dp), allocatable, dimension(:,:) :: S, S1, S2, S3, Mmat, LLbm, LUbm, &
    LLbm1, LLBm2, LLbm3, ungradphi, u_n, u_nm1, gradc, gradc2, ung, &
    emat, gradmu, gradcdotgradphi, iso_mat, gradcgradc_mat, u_iter

  real(dp), allocatable, dimension(:,:,:) :: gradcgradc, tau_c

  real(dp), allocatable, dimension(:,:,:,:) :: gradcgradphi

end module diffuse_interface_globals_m
