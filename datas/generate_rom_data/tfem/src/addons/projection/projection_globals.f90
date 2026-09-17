
! Copyright (C) 2012-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Global variables for the projections elements

module projection_globals_m

  use math_defs_m
  use gauss_defs_m
  use shapefunc_m, only: shapefunc_t

  implicit none

  save


! OLD MESH

! global element shape old mesh

  character (len=13) :: globalshapev = ''

! space dimension
  integer :: ndim = 0

! interpolation number (shape function) for the input vector
! (vector of type vector_t)
  integer :: intpolv = 0

! number of degrees of freedom of the input vector (one direction only)
  integer :: ndfv = 0

! number of degrees of freedom of the input vector (one point only)
  integer :: ncompv = 0

! number of nodal points of the interior element (old mesh)
  integer :: nodalpv = 0

! layer number
  integer :: layer = 0


! NEW MESH

! global element shape

  character (len=13) :: globalshape = ''

! coordinate system: 0=Cartesian (x,y), 2D
!                    1=Cylindrical (z,r), axisymmetrical
!                    2=Cartesian (x,y,z), 3D
  integer :: coorsys = 0

! interpolation number (shape function) for the shape of the element
  integer :: intpols = 0

! number of shape degrees of freedom with respect to one direction only
  integer :: ndfs = 0

! interpolation number (shape function) for the output vector
! (solution of type sysvector_t)
  integer :: intpol = 0

! number of degrees of freedom with respect to one direction only
! (solution of type sysvector_t)
  integer :: ndf = 0

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

! number of nodal points of the interior element (new mesh)
  integer :: nodalp = 0


! various allocatable arrays:

  integer, allocatable, dimension(:,:) :: grpelm_n, grpelm_nm1

  real(dp), allocatable, dimension(:) :: wg, detF, u, un

  real(dp), allocatable, dimension(:,:) :: xig, phis, x, xg, xig_n, phi_n, &
    ung, phi, xrnod

  real(dp), allocatable, dimension(:,:,:) :: dphis, F, Finv



! various work arrays for general use:

  real(dp), allocatable, dimension(:) :: work, work1, work8

  real(dp), allocatable, dimension(:,:) :: tmp, work2, work4, work10, work11

  real(dp), allocatable, dimension(:,:,:) :: work6, work7, work9

  real(dp), allocatable, dimension(:,:,:,:) :: work5

  real(dp), allocatable, dimension(:,:,:,:,:) :: work3


! gauss types

  type(gauss_t) :: gauss


! shapefunc types

  type(shapefunc_t) :: shapefunc, shapefuncv, shapefuncs


end module projection_globals_m
