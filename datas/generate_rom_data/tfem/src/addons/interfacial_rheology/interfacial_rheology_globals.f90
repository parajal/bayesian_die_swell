!
! Copyright (C) 2004-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system (e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Global variables for the interfacial rheology elements

module interfacial_rheology_globals_m

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

! include third velocity component in a 2D space (coorsys=0,1)
  integer :: vel3D = 0

! interpolation number (shape function) for the velocity
  integer :: intpol = 0

! number of velocity degrees of freedom with respect to one direction only
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

! number of nodal points of the interior element
  integer :: nodalp = 0

! layer number
  integer :: layer = 0


! various allocatable arrays:

  real(dp), allocatable, dimension(:) :: wg, curvel, surfl, &
    surfl0, Jinv, Jdet, xar, u

  real(dp), allocatable, dimension(:,:) :: xig, phi, x, xg, &
    Suu, Suv, Svv, Suw, Svw, Sww, normal, rhsf, xe, xn, xng, &
    x0, xg0, xrnod, uvec

  real(dp), allocatable, dimension(:,:,:) :: dphi, dxdxis, &
    gi_up, gij_up, gi_up0, dxdxis0, Fs, Bs, dxndxis


! various work arrays for general use:

  real(dp), allocatable, dimension(:) :: work, work_I, work_B, &
    work_gijup, work_axi, trBs, work1

  real(dp), allocatable, dimension(:,:) :: tmp, work2, work3, &
    work4, work5, work6, work7, work8

  real(dp), allocatable, dimension(:,:,:) :: work_upup, &
    work_upup2, work_updown, work9, work10, work11, work12

! gauss types

  type(gauss_t) :: gauss


! shapefunc types

  type(shapefunc_t) :: shapefunc


end module interfacial_rheology_globals_m
