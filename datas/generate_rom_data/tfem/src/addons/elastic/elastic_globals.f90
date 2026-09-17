
! Copyright (C) 2007-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Global variables for the Mooney-Rivlin incompressible elastic solid

module elastic_globals_m

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

! interpolation number (shape function) for the displacement
  integer :: intpol = 0

! number of displacement degrees of freedom with respect to one direction only
  integer :: ndf = 0

! interpolation number (shape function) for the pressure
  integer :: intpolp = 0

! number of pressure degrees of freedom
  integer :: ndfp = 0

! interpolation number (shape function) for the Lagrange multiplier
  integer :: intpoll = 0

! secondary integration rule for the interior of an element
  integer :: intrule2 = 0

! number of Lagrange multiplier degrees of freedom
  integer :: ndfl = 0

! integration rule for the interior of an element
  integer :: intrule = 0

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

! physical quantity number of the displacement and pressure in the system vector
  integer :: physqdisp = 0, physqpress = 0

! layer number
  integer :: layer = 0

! Elastic model: 1=neo-Hookean, 2=Mooney-Rivlin
  integer :: model = 1

! the material parameters for the elastic model
  real(dp) :: C1, C2


! various allocatable arrays:

  real(dp), allocatable, dimension(:) :: wg, detF0, detFc, Jv, u, p, &
    pg, fpg, Jwg, tr, curvel, surfl
  real(dp), allocatable, dimension(:,:) :: xig, phi, x, xg, fg, psi, &
    Lu, Lv, Lw, x0, xrnod, xg0
  real(dp), allocatable, dimension(:,:,:) :: dphi, F, Finv, dphidx, F0, &
    F0inv, Fc, Fcinv, B, Binv, tau
  real(dp), allocatable, dimension(:,:,:,:) :: Smat


! various work arrays for general use:

  real(dp), allocatable, dimension(:) :: work4

  real(dp), allocatable, dimension(:,:) :: work1, work2, work5, tmp

  real(dp), allocatable, dimension(:,:,:) :: work6

  real(dp), allocatable, dimension(:,:,:,:) :: work, work3


! gauss types

  type(gauss_t) :: gauss


! shapefunc types

  type(shapefunc_t) :: shapefunc, shapefuncp, shapefuncl

end module elastic_globals_m
