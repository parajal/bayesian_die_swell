
! Copyright (C) 2020-2020 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Global variables for the linear elastic elements

module linear_elastic_globals_m

  use math_defs_m
  use gauss_defs_m
  use shapefunc_m, only: shapefunc_t
  use eltree_basic_m, only: eltree_t
  use linear_elastic_material_m, only: lemodel_t

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

! include third displacement component in a 2D space (coorsys=0,1)
  integer :: disp3D = 0

! interpolation number (shape function) for the displacement
  integer :: intpol = 0

! number of displacement degrees of freedom with respect to one direction only
  integer :: ndf = 0

! number of displacement degrees of freedom (one point only)
  integer :: ncompu = 0

! number of displacement degrees of freedom on the boundary of an interior
! element
  integer :: ndfb = 0

! interpolation number (shape function) for the pressure
  integer :: intpolp = 0

! number of pressure degrees of freedom
  integer :: ndfp = 0

! interpolation number (shape function) for the pressure projection space
  integer :: intpolp1 = 0

! number of pressure degrees of freedom of the projection space
  integer :: ndfp1 = 0

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

! number of nodal points on of an element of an object
  integer :: nodalpo = 0

! number of sides of an interior element
  integer :: nsides = 0

! physical quantity number of the displacement and pressure in the system vector
  integer :: physqdisp = 0, physqpress = 0

! layer number
  integer :: layer = 0


! various allocatable arrays:

  integer, allocatable, dimension(:) :: posp, pos1, pos2

  integer, allocatable, dimension(:,:) :: pos, grpelm_n, grpelm_nm1

  real(dp), allocatable, dimension(:) :: wg, detF, curvel, u, surfl, dudy, &
    dvdx, st, pr, sigmatt, press, w, detF_n, p_exact, trac, csg, car, gammacg, &
    prs, dTemp, fp

  real(dp), allocatable, dimension(:,:) :: xig, phi, x, xg, xrnod, fg, &
    psi, Suu, Suv, Svv, Lu, Lv, theta, Bu11, Bu12, Amat, zeta, &
    Cu11, Cu12, Cu13, Dmat, normal, gradpsi, Suw, Svw, Sww, Lw, uvector, &
    dxdxi, xigo, xel, wgradpsi, wvec, ubar, theta1, xg_n, xg_nm1, &
    xig_n, phi_n, xig_nm1, phi_nm1, g1_up, psi1, u_exact, ugvector, &
    Cmat, fuv

  real(dp), allocatable, dimension(:,:,:) :: dphi, F, Finv, dphidx, Bmat, &
    dxdxis, gradu, tauten, sigmaten, dtheta, dthetadx, dpsi, dpsidx, &
    dtheta1, dtheta1dx, dphi_n, F_n, Finv_n, gi_up, Dten


! various work arrays for general use:

  real(dp), allocatable, dimension(:) :: work, work1, work8

  real(dp), allocatable, dimension(:,:) :: tmp, work2, work4, work10, work11

  real(dp), allocatable, dimension(:,:,:) :: work6, work7, work9

  real(dp), allocatable, dimension(:,:,:,:) :: work5

  real(dp), allocatable, dimension(:,:,:,:,:) :: work3


! gauss types

  type(gauss_t) :: gauss


! shapefunc types

  type(shapefunc_t) :: shapefunc, shapefuncp, shapefuncl, shapefunco
  type(shapefunc_t) :: shapefuncp1

! model parameters

  type(lemodel_t) :: lemodel

contains

  subroutine errormsg_notdisp3D ( name_of_routine )

!   name of calling routine
    character (len=*), intent(in) :: name_of_routine

    write(*,'(2a/)') 'Error: disp3D not set for routine: ', &
       trim(name_of_routine)
    stop

  end subroutine errormsg_notdisp3D

  subroutine errormsg_notmixed_incompressible ( name_of_routine )

!   name of calling routine
    character (len=*), intent(in) :: name_of_routine

    write(*,'(/a/2a/)') &
      'Error: mixed (u,p) elements needed for incompressible materials.', &
      ' Calling routine: ', trim(name_of_routine)
    stop

  end subroutine errormsg_notmixed_incompressible

end module linear_elastic_globals_m
