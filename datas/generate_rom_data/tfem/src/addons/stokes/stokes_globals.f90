
! Copyright (C) 2005-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Global variables for the Stokes elements

module stokes_globals_m

  use math_defs_m
  use gauss_defs_m
  use shapefunc_m, only: shapefunc_t
  use eltree_basic_m, only: eltree_t

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

! number of velocity degrees of freedom (one point only)
  integer :: ncompu = 0

! number of velocity degrees of freedom on the boundary of an interior element
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

! interpolation number (shape function) for the stress in mixed stokes
  integer :: intpolt = 0

! number of stress degrees of freedom with respect to one direction only
! (mixed stokes element only)
  integer :: ndft = 0

! interpolation number (shape function) for the gradient E in DEVSS
  integer :: intpole = 0

! number of degrees of freedom of E in DEVSS (one direction only)
  integer :: ndfe = 0

! number of degrees of freedom of E in DEVSS (one point only)
  integer :: ncompe = 0

! interpolation number (shape function) for the gradient G in DEVSS-G
  integer :: intpolg = 0

! number of degrees of freedom of G in DEVSS-G (one direction only)
  integer :: ndfg = 0

! number of degrees of freedom of G in DEVSS-G (one point only)
  integer :: ncompg = 0

! physical quantity number of the velocity and pressure in the system vector
  integer :: physqvel = 0, physqpress = 0

! physical quantity number of E or G  in the sysvector (DEVSS and DEVSS-G)
  integer :: physqgrad = 0

! physical quantity number of the stress in the sysvector (mixed stokes only)
  integer :: physqstress = 0

! layer number
  integer :: layer = 0

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

  integer, allocatable, dimension(:) :: posp, pos1, pos2

  integer, allocatable, dimension(:,:) :: pos, grpelm_n, grpelm_nm1

  real(dp), allocatable, dimension(:) :: wg, detF, curvel, u, surfl, dudy, &
    dvdx, st, pr, sigmatt, press, w, detF_n, p_exact, trac, csg, car, gammacg, &
    um, gtmp

  real(dp), allocatable, dimension(:,:) :: xig, phi, x, xg, xrnod, fg, &
    psi, Suu, Suv, Svv, Lu, Lv, theta, Bu11, Bu12, Amat, zeta, &
    Cu11, Cu12, Cu13, Dmat, normal, gradpsi, Suw, Svw, Sww, Lw, uvector, &
    dxdxi, xigo, xel, wgradpsi, wvec, ubar, theta1, xg_n, xg_nm1, &
    xig_n, phi_n, xig_nm1, phi_nm1, g1_up, psi1, u_exact, ugvector, &
    eigval, umesh, u_mesh

  real(dp), allocatable, dimension(:,:,:) :: dphi, F, Finv, dphidx, Bmat, &
    dxdxis, gradu, tauten, sigmaten, dtheta, dthetadx, dpsi, dpsidx, &
    dtheta1, dtheta1dx, dphi_n, F_n, Finv_n, gi_up, Dten, &
    eigvec, ggrad


! various work arrays for general use:

  real(dp), allocatable, dimension(:) :: work, work1, work8

  real(dp), allocatable, dimension(:,:) :: tmp, work2, work4, work10, work11, &
    work13

  real(dp), allocatable, dimension(:,:,:) :: work6, work7, work9

  real(dp), allocatable, dimension(:,:,:,:) :: work5, work12

  real(dp), allocatable, dimension(:,:,:,:,:) :: work3, work15

  real(dp), allocatable, dimension(:,:,:,:,:,:) :: work14


! gauss types

  type(gauss_t) :: gauss, gauss_ie


! shapefunc types

  type(shapefunc_t) :: shapefunc, shapefuncp, shapefuncl, shapefunco
  type(shapefunc_t) :: shapefunce, shapefuncg
  type(shapefunc_t) :: shapefunct
  type(shapefunc_t) :: shapefuncp1

! eltree type (pointer alias)

  type(eltree_t), pointer :: eltree1

contains

  subroutine errormsg_notvel3D ( name_of_routine )

!   name of calling routine
    character (len=*), intent(in) :: name_of_routine

    write(*,'(2a/)') 'vel3D not set for routine: ', &
       trim(name_of_routine)
    stop

  end subroutine errormsg_notvel3D

end module stokes_globals_m
