
! Copyright (C) 2005-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Global variables for the viscoelastic elements

module viscoelastic_globals_m

  use kind_defs_m
  use gauss_defs_m
  use shapefunc_m, only: shapefunc_t
  use stokes_globals_m
  use viscoelastic_models_defs_m

  implicit none

  save

! use log conformation transformation?
  logical :: logc = .false.

! use projection of exp(s) in stress-implicit methods
  logical :: exps = .false.

! use b-tensor formulation?
  logical :: bten = .false.

! use projection of c=b.b^T in stress-implicit methods
  logical :: cproj = .false.

! compute c in integration points directly from c=exp(s) or c=b.b^T in
! stress-implicit methods
  logical :: c_direct = .false.

! add c/lambda to both sides of the constitutive equation
  logical :: coverl = .false.

! variable model parameters
  logical :: varpar = .false.

! model is dependent on J
  logical :: vemcompressible = .false.

! model is dependent on the equivalent plastic strain gammap
  logical :: vemgammap = .false.

! interpolation number (shape function) for the tensor components
  integer :: intpolc = 0

! number of conformation degrees of freedom with respect to one component only
  integer :: ndfc = 0

! number of degrees of freedom of conformation tensor for one mode
! (one point only)
  integer :: ncompc = 0

! number of degrees of freedom of tensor that is convected (b or c) for one mode
! (one point only).
  integer :: ncomp = 0

! number of modes
  integer :: nmodes = 0

! mode numbers mode1--mode2 that need to be computed
  integer :: mode1 = 1, mode2 = 0

! number of sub mode numbers = max(mode2-mode1+1,0)
  integer :: nsubm = 0

! start of the physical quantities for the conformation (only when part of
! the system vector).
  integer :: physqc = 0

! physical quantity number for the equivalent plastic strain
  integer :: physqgammap = 0

! number of components of the total stress (one point only)
  integer :: ncompt = 0

! number of integration points along the boundary of an element for DG
  integer :: nintbc = 0

! parameters for the h/U factor in SUPG
  integer :: htype, Uscaling

! storage type of the conformation
  integer :: cstorage = 0

! storage type of the equivalent plastic strain
  integer :: gpstorage = 0

! various allocatable arrays:

  integer, allocatable, dimension(:) :: posc, posu, posg

  real(dp), allocatable, dimension(:) :: g, uvecc, Gmod, tau, hoverU, &
    eigvalue, Gmodl, lambda_small, facv, cm, &
    gpn, gpng, gpnm1, gpnm1g, gpiter, gpiterg, unp1gradgpiter, ungradgpn, &
    fgpiterg, fgpng, gprhsmodel, gpnod, unm1gradgpnm1, fgpnm1g, gphatg, &
    fgphatg

  real(dp), allocatable, dimension(:,:) :: uvecn, uvecnm1, &
    cn, cnm1, gvecn, gvecnm1, c, tauvec, s, snod, cnod, ungradtheta, &
    uhatgradtheta, unm1gradtheta, uhat, ugsn, curvels, us, cnside, dDi, eigv, &
    chat, uvecnp1, gvecnp1, unp1gradtheta, uvecmeshnp1, uvecmeshn, &
    uvecmeshnm1, theta_n, theta_nm1, b, bnod, sg, citer, uvec_supg, &
    ugradtheta_supg, dhoverU, dtauvecdJ, gradgpiter

  real(dp), allocatable, dimension(:,:,:) :: cng, cnm1g, &
    fng, fnm1g, cg, cten, hten, rhsd, ungradcn, ungradcnm1, &
    unm1gradcnm1, xigs, phis, xs, xgs, thetas, normals, &
    chatg, fhatg, rhsmodel, dtheta_n, dthetadx_n, dtauvec, &
    citerg, fiterg, unp1gradciter, RT, dtautendJ

  real(dp), allocatable, dimension(:,:,:,:) :: gradcn, gradcten, cnjump, &
    dfiterg, gradciter, drhsdLg

! arrays needed for local material parameters

  real(dp), allocatable, dimension(:) :: lambda_g_mode1, fac2_g

  real(dp), allocatable, dimension(:,:) :: Gmod_g, lambda_g, lambda_small_g, &
    Gmodl_g

  real(dp), allocatable, dimension(:,:,:) :: nonlin_g


! special allocatables for boundary integrals in DG (needs generalization):
  real(dp), allocatable, dimension(:) :: wgb, xigb
  real(dp), allocatable, dimension(:,:) :: phib, dphib



! viscoelastic model type

  type(vemodel_t) :: vemodel
  type(vemodel_t), allocatable, dimension(:) :: mvemodel


! shapefunc types

  type(shapefunc_t) :: shapefuncc


! vemopt structure for optional arguments hided in a structure

  type(vemopt_t) :: vemopt, vemoptn, vemoptnm1, vemoptnp1, vemoptiter
  type(vemopt_gammap_t) :: vemopt_gammap, vemoptn_gammap, vemoptnm1_gammap, &
    vemoptnp1_gammap, vemoptiter_gammap

contains

  subroutine errormsg_notvel3D_coorsys ( name_of_routine )

!   name of calling routine
    character (len=*), intent(in) :: name_of_routine

    write(*,'(/a,2(a/))') 'Error in ', trim(name_of_routine)//':', &
       ' coorsys == 1 .and. vel3D == 1 not yet available.'
    stop

  end subroutine errormsg_notvel3D_coorsys

  subroutine warningmsg_vel3D_coorsys ( name_of_routine )

!   name of calling routine
    character (len=*), intent(in) :: name_of_routine

    write(*,'(/a,5(a/))') 'Warning in ', trim(name_of_routine)//':', &
       ' vel3D == 1 not tested (for both coorsys=0 and 1).', &
       ' Use at your own risk!', &
       ' If the routine works as expected and you think this warning ', &
       ' message is not needed anymore, notify me (m.a.hulsen@tue.nl).'

  end subroutine warningmsg_vel3D_coorsys

end module viscoelastic_globals_m
