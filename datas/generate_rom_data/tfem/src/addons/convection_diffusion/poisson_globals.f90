
! Copyright (C) 2007-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Global variables for the Poisson elements

module poisson_globals_m

  use math_defs_m
  use gauss_defs_m
  use shapefunc_m, only: shapefunc_t

  implicit none

  save

! global element shape

  character (len=13) :: globalshape = ''

! space dimension
  integer :: ndim = 0

! coordinate system: 0=Cartesian (x,y), 2D
!                    1=Cylindrical (z,r), axisymmetrical
!                    2=Cartesian (x,y,z), 3D
  integer :: coorsys = 0

! interpolation number (shape function)
  integer :: intpol = 0

! number of degrees of freedom with respect to one direction only
  integer :: ndf = 0

! interpolation number (shape function) for the Lagrange multiplier
  integer :: intpoll = 0

! number of Lagrange multiplier degrees of freedom
  integer :: ndfl = 0

! integration rule for the interior of an element
  integer :: intrule = 0

! secondary integration rule for the interior of an element
  integer :: intrule2 = 0

! integration type
  integer :: inttype = 0

! number of integration points for the interior of an element
  integer :: ninti = 0

! number of subdomains for the integration (generalized integration rules)
  integer :: nsubint = 0

! number of nodal points of the interior element
  integer :: nodalp = 0

! number of nodal points on a single side of the boundary of an interior element
  integer :: nodalpb = 0

! number of sides of an interior element
  integer :: nsides = 0

! layer number
  integer :: layer = 0


! various allocatable arrays:

  real(dp), allocatable, dimension(:) :: wg, fg, detF, curvel, u, surfl, sigma

  real(dp), allocatable, dimension(:,:) :: xig, phi, x, xg, xrnod, dxdxi, &
    normal, psi

  real(dp), allocatable, dimension(:,:,:) :: dphi, F, Finv, dphidx, dxdxis


! various work arrays for general use:

  real(dp), allocatable, dimension(:) :: work

  real(dp), allocatable, dimension(:,:) :: tmp


! gauss types

  type(gauss_t) :: gauss


! shapefunc types

  type(shapefunc_t) :: shapefunc, shapefuncl


end module poisson_globals_m
