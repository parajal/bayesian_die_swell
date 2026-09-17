
! Copyright (C) 2012-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Global variables for the convection-diffusion elements

module structures_globals_m


  use math_defs_m

  implicit none

  save

! global element shape

  character (len=13) :: globalshape = ''

! space dimension
  integer :: ndim = 0

! number of degrees of freedom with respect to one direction only (truss)
  integer :: ndf = 0

! interpolation number (shape function) for the axial displacement
  integer :: intpolu = 0

! number of axial degrees of freedom
  integer :: ndfu = 0

! interpolation number (shape function) for the lateral displacement
  integer :: intpolv = 0

! number of lateral degrees of freedom with respect to one direction only
  integer :: ndfv = 0

! interpolation number (shape function) for the axial (torsion) rotation
  integer :: intpolt = 0

! number of axial rotation degrees of freedom
  integer :: ndft = 0

! number of nodal points of the interior element
  integer :: nodalp = 0

! number of integration points for the interior of an element
  integer :: ninti = 0


  integer, allocatable, dimension(:) :: pos1, pos2, por1, por2, pos, &
    posu, posv, posw, post

  real(dp), allocatable, dimension(:) :: e, n, u, v, xig, wg, qg, Iz, dIz, Ac, &
    fu, fv, fw, qu, qv, qw, sol, fe, qnod, Izz, Iyz, Iyy, nz, ny, Jrr, &
    dIzz, dIyz, dIyy, w, t

  real(dp), allocatable, dimension(:,:) :: x, emat, uvec, psi_r, d2psi_r, xg, &
    psi, d2psi, xrnod, dpsi_r, d3psi_r, dpsi, d3psi, phi, dphi_r, dphi, &
    zeta, dzeta_r, dzeta, qvecg, Suu, Svv, Svw, Sww, Stt, Q, qvecnod, Qr

! work arrays

  real(dp), allocatable, dimension(:,:) :: work

end module structures_globals_m
