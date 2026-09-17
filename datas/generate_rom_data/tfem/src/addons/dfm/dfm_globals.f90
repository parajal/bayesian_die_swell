
! Copyright (C) 2016-2017 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Global variables for the viscoelastic elements

module dfm_globals_m

  use kind_defs_m
  use gauss_defs_m
  use shapefunc_m, only: shapefunc_t
  use stokes_globals_m

  implicit none

  save

! logicals to save CPU time to avoid compute terms depend on I1, I2 or B^-1
  logical :: I1dep =.true., I2dep = .true., Binv = .true.

! interpolation number (shape function) for the deformation field components
  integer :: intpolf = 0

! number of deformation field degrees of freedom with respect to one component
! in one field only
  integer :: ndff = 0

! number of degrees of freedom of F for one field (one point only)
  integer :: ncompf = 0

! deformation tensor interpolation for DG (in age coordinate tau)
  integer :: intpoltau = 0

! number of intervals N in the age coordinate tau
  integer :: nintvaltau = 0

! number of fields within an interval in the age coordinate tau (=number of
! Gauss points within an interval)
  integer :: ninttau = 0

! number of deformation fields
  integer :: nfields = 0

! time integration
  integer :: timeint = 0

! number of components of the total stress (one point only)
  integer :: ncompt = 0

! integral model type
  integer :: intmodel = 0

! integral model spectrum
  integer :: spectrum = 0

! number of modes in a spectrum
  integer :: nummodes = 0

! number of discrete modes in a combined spectrum
  integer :: numdismodes = 0

! the damping function
  integer :: dampingf = 0

! real values

  real(dp) :: mf0h0 ! coefficient of the unity tensor for dottau

! various allocatable arrays:

  integer, allocatable, dimension(:) :: posf, pos_s

  real(dp), allocatable, dimension(:) :: taug, wtaug, g, tauSUPG, hoverU, &
    fac1, fu, taug2, wtaug2, Gmod, lam

  real(dp), allocatable, dimension(:,:) :: phitau, dphitau, phitaulr, &
    uvec, gvec, ugradtheta, uvecmesh, fn, fnm1, fnm2, fmat, phitaunod, &
    mmatrix, nmatrix, Fjump2, mf, dmf, tau_gp, tau_gw, stressvec, s, &
    uvecmeshnp1, rhs, taun, relaxstressvec

  real(dp), allocatable, dimension(:,:,:) :: Fjump, fn2, fn2m1, fn2m2, fnp1, &
    fn2g, fn2m1g, fn2m2g, fhat2g, fnp1g, rhs1g, rhs2g, ivr, dpf, gradtaun, &
    taunten, rhsten, Fvec1

  real(dp), allocatable, dimension(:,:,:,:) :: fng, fnm1g, fnm2g, fhatg, rhsg, &
    rhs3g, Fvec, Ften, Bten, gradtauten, ddpf

  real(dp), allocatable, dimension(:,:,:,:,:) :: doublestressten

! shapefunc types

  type(shapefunc_t) :: shapefuncf


end module dfm_globals_m
