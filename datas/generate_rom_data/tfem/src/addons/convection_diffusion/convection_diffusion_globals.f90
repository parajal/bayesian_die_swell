
! Copyright (C) 2007-2011 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Global variables for the convection-diffusion elements

module convection_diffusion_globals_m

  use poisson_globals_m
  use eltree_basic_m, only: eltree_t

  implicit none

  save


! the integration scheme on interface elements in the eltree:
!   intrule_ie: the integration rule
!   ninti_ie: number of integration points
!   x_ie, w_ie: points and weights
  integer :: intrule_ie, ninti_ie
  real(dp), allocatable :: w_ie(:), x_ie(:,:)

! wng: weight * reference normal in the integration points on the interface
! dan: weight * normal in the integration points on the interface
!              (after mapping)
! da: weight in the integration points on the interface (after mapping)
  real(dp), allocatable :: wng(:,:), dan(:,:), da(:)


! various allocatable arrays:

  integer, allocatable, dimension(:) :: pos1, pos2

  real(dp), allocatable, dimension(:) :: alphag, gammag, cn, cng, betag, &
    ugradcn, cbar, jbar, ibar, ung

  real(dp), allocatable, dimension(:,:) :: uvecg, ugradphi, gradcn

  real(dp), allocatable, dimension(:,:,:) :: alphatg


! various work arrays for general use:

  real(dp), allocatable, dimension(:,:) :: work2, work4

  real(dp), allocatable, dimension(:,:,:) :: work6, work7


! gauss types

  type(gauss_t) :: gauss_ie

! eltree type (pointer alias)

  type(eltree_t), pointer :: eltree1


end module convection_diffusion_globals_m
