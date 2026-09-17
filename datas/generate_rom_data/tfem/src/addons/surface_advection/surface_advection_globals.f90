
! Copyright (C) 2011-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Global variables for the surface advection elements

module surface_advection_globals_m

  use tfem_elem_m

  implicit none

  save

! global element shape

  character (len=13) :: globalshape = ''

! space dimension
  integer :: ndim = 0

! interpolation number (shape function) for the position
  integer :: intpol = 0

! interpolation number (shape function) for the velocity
  integer :: intpolv = 0

! shape function used for the "isoparametric" deformation of the element
! 0=position, 1=velocity
  integer :: isoshape = 0

! integration rule for the interior of an element
  integer :: intrule = 0

! integration type
  integer :: inttype = 0

! number of integration points
  integer :: ninti = 0

! number of subdomains for the integration (generalized integration rules)
  integer :: nsubint = 0

! number of nodal points
  integer :: nodalp = 0

! number of degrees of freedom of the height function
  integer :: ndf = 0

! number of degrees of freedom of the velocity
  integer :: ndfv = 0

! number of velocity components
  integer :: ncompv = 0

! various allocatable arrays

  real(dp), allocatable, dimension(:) :: wg, fg, detF, xig1, Hn, Hng, fac, &
    Hnm1, Hnm1g, tau, hoverU, Hf

  real(dp), allocatable, dimension(:,:) :: phi, x, xg, uvecg, xig, ugradphi, &
    xrnod, phiv, Hngl, Hnm1gl, work

  real(dp), allocatable, dimension(:,:,:) :: dphi, dphiv, F, Finv, dphidx

! gauss types

  type(gauss_t) :: gauss

! shapefunc types

  type(shapefunc_t) :: shapefunc, shapefuncv

end module surface_advection_globals_m
