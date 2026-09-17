
! Copyright (C) 2005-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

!
! definitions for viscoelastic models
!

module viscoelastic_models_defs_m

  use glob_defs_m
  use set_optional_m

  implicit none


! type definition of a visco elastic model

  type vemodel_t

    logical :: created = .false. ! structure created?
    logical :: deviatoric = .false. ! replace L with L^d = L - (tr L)/3 I
    logical :: compressible = .false. ! model depends on the input of J
    logical :: vmglobal = .false. ! adapted lambda uses global vonmises stress
    logical :: gammap = .false. ! model depends on plastic strain
    integer :: ncompb = 0  ! number of components of the b-tensor
    integer :: ncompc = 0  ! number of components of the conformation tensor
    integer :: ncompt = 0  ! number of components of the stress tensor
    integer :: nmodes = 0  ! number of modes
    integer :: flowtype = 0  ! the flow type:
                           !     0=2D
                           !     1=2D, axisymmetrical
                           !     2=3D
    integer :: bvariant = 0  ! Variant of the b-equation (for the b-tensor):
!                          !    1=CDT (Carrozza et al.)
!                          !    2=symmetric (Balci et al.)
!                          !    3=Cholesky (Vaithianathan & Collins)
!                          !    4=Cholesky with log (Vaithianathan & Collins)
    integer :: ixsi = 0    ! Position of the slip parameter xsi in the
                           ! nonlinear parameters stored in nonlin.
                           ! ixsi=0 means no slip.
    integer :: norlxmode = 0  ! mode number where the relaxation term will be
                              ! set to zero.
                              ! norlxmode = 0 means none
    integer :: gammap_mode = 0  ! mode number for the computation
!                               ! of the plastic strain gamma_p for
!                               ! alam_gammap_model == 1
    integer :: model = 0   ! the model number:
!
!           model     name                    nonlinear parameters
!
!              1      Leonov
!              2      UCM/Oldroyd-B
!              3      Giesekus                alpha
!              4      Larson                  gamma
!              5      Phan-Thien/Tanner
!                     linear factor           epsilon
!              6      Phan-Thien/Tanner
!                     expon. factor           epsilon
!              7      Modified Leonov         alpha    beta
!              8      Chilcott-Rallison       L^2
!              9      FENE-P (macro)          b
!             10      Johnson-Segalman        xsi
!             11      Phan-Thien/Tanner
!                     linear factor
!                     with slip               epsilon  xsi
!             12      Phan-Thien/Tanner
!                     expon. factor
!                     with slip               epsilon  xsi
!             13      Marrucci                beta
!             14      XPP double equation     lambdas  nu
!             15      Rolie-Poly              lambdar  beta  delta
!             16      XPP single equation     lambdas  nu
!             17      XPP log-conformation    lambdas  nu
!                     backbonestretch
!             18      PTT-XPP single equation lambdas  nu
!                     (XPP compatible: four
!                     components in 2D)
!             19      PTT-XPP single equation lambdas  nu
!                     (three components in 2D)
!             20      FENE-P (macro)          b
!                     (consistent with Wapperom & Hulsen 1998)
!             21      Chilcott-Rallison       L^2
!                     (G, lambda are of linear Maxwell model)
!             22      Saramito (2021) Drucker-Prager plasticity
!                     elastoviscoplastic      tau_y    mu
!             23      EGP for incompressible
!                     flow
!             24      EGP for compressible
!                     flow

    integer :: alam_model = 0
!   The adapted lambda model number.
!   For alam_model > 0 the factor 1/lambda is replaced by:
!
!           1   0    elastic model (no parameters)
!
!           2     1           tau_d-tau_y       elastoviscoplastic
!              -------- max(0,-----------)      Saramito (2007)
!               lambda           tau_d          parameter: tau_y
!
!           3     G         tau_d-tau_y         elastoviscoplastic
!               ----- max(0,-----------)^1/n    Saramito (2009)
!               tau_d            K              parameters: tau_y, K, n
!
!           4     1     sinh(tau_d/tau_ref)
!              -------- -------------------           Eyring
!               lambda     tau_d/tau_ref              parameter: tau_ref
!
!           5     1             tau_d/tau_ref1
!              -------- [ f1 --------------------     Ree-Eyring
!               lambda       sinh(tau_d/tau_ref1)     parameters: tau_ref1, f1
!                                                                 tau_ref2
!                               tau_d/tau_ref2     -1 note: f2=1-f1
!                       + f2 -------------------- ]
!                            sinh(tau_d/tau_ref2)

!   Here, tau_d is the von-Mises equivalent shear stress:
!
!      (tau_d)^2 = [(tau_xx-tau_yy)^2+(tau_yy-tau_zz)^2+(tau_zz-tau_xx)^2]/6+
!                        (tau_xy)^2+(tau_yz)^2+(tau_zx)^2

    integer :: alamJ_model = 0
!   The adapted lambda model number for the dependence on J.
!   For alam_model > 1 (see above) there is an additional factor h(J) for
!   the replacement of 1/lambda. Thus if the factor f(tau_d)
!   represents the adapted lambda for alam_model > 1, the replacement for
!   alamJ_model > 0 becomes h(J)*f(tau_d).
!
!   For alamJ_model > 0 the factor h(J) is given by:
!
!           1   J^beta       EGP viscosity pressure dependence (see userguide)
!                            parameter: beta

    integer :: alam_gammap_model = 0
!   The adapted lambda model number for the dependence on the plastic strain
!   gamma_p. For alam_model > 1 (see above) there is an additional factor
!   e(gamma_p) for the replacement of 1/lambda. Thus if the factor f(tau_d)
!   represents the adapted lambda for alam_model > 1, the replacement for
!   alam_gammap_model > 0 becomes e(gamma_p)*f(tau_d).
!
!   For alam_gammap_model > 0 the factor e(gamma_p) is given by:
!
!           1   exp(-Sa R)   EGP viscosity plastic strain dependence
!                            (see userguide)
!                            where R is function of gamma_p:
!
!               R = (1+(r0 exp(gamma_p))^r1)^(r2-1)/r1
!                   ----------------------------------
!                         (1+r0^r1)^(r2-1)/r1
!
!                            parameters: Sa, r0, r1, r2


!   material parameters determining the linear spectrum:
!     modulus(mode):  G modulus parameter for mode
!     lambda(mode):  relaxation time parameter for mode
    real(dp), dimension(:), allocatable :: modulus, lambda

!   nonlinear material parameters
!     nonlin(i,mode) the ith nonlinear parameter for mode
    real(dp), dimension(:,:), allocatable :: nonlin

!   material parameters for the adapted lambda models
!     alam(i,mode) the ith adapted lambda parameter for mode
    real(dp), dimension(:,:), allocatable :: alam

!   material parameters for the adapted lambda dependent on J models
!     alamJ(i,mode) the ith adapted lambda parameter for mode
    real(dp), dimension(:,:), allocatable :: alamJ

!   material parameters for the adapted lambda dependent on gamma_p models
!     alam_gammap(i,mode) the ith adapted lambda parameter for mode
    real(dp), dimension(:,:), allocatable :: alam_gammap

  end type vemodel_t


! type definition for the time discretization in steady routines for the
! b-formulation

  type tdpar_t

!   time discretization method
!   1: BDF1 (implicit Euler)
!   2: BDF2 (2nd order implicit Gear)
    integer :: method = 1

!   discrete time step
    real(dp) :: timestep = 1.e-3_dp

  end type tdpar_t


! type definition for using optional input/output of data to the main
! routines (multi-point, multi-mode) of the constitutive viscoelastic model.
! This can be easily extended without adding more arguments to the routines.

  type vemopt_t

    logical :: created = .false. ! structure created?

!   take dependence of the model on relative change in volume J into account.
    logical :: dep_J = .false.

!   take dependence of the model on the plastic strain gammap into account.
    logical :: dep_gammap = .false.

!   compute the derivative of the right-hand side with respect to the
!   conformation tensor without mode coupling.
    logical :: compute_drhs = .false.

!   compute the derivative of the right-hand side with respect to the
!   conformation tensor with mode coupling.
!   NOTE: if the model does not have mode coupling, drhs_mm wil contain a
!         "diagonal" matrix version of drhs.
    logical :: compute_drhs_mm = .false.

!   compute the derivative of the right-hand side wrt to J
    logical :: compute_drhsdJ = .false.

!   compute the derivative of the right-hand side wrt to gammap
    logical :: compute_drhsdgammap = .false.

!   actual arguments in the call to the vonmises_mm routine in the
!   rhs_adapted_lambda subroutine. In this way the modes used for the
!   computation of the von Mises stress can be imposed by the user.
!   NOTE: vmmode1=0, vmmode2=0 is the default, which means that the actual
!         values used for vmmode1 and vmmode2 are taken from the actual
!         mode1, mode2 values in the routine that uses them.
!   NOTE: setting vmmode1/=mode1, vmmode2/=mode2 is possible and creates
!         special models, where vmmode1/vmmode2 can be even be outside
!         mode1:mode2 range. However, this will not be supported by the
!         flow modules for now.
!   intent(in)
    integer :: vmmode1 = 0, vmmode2 = 0

!   relative change in volume J.
!   J(ip) is the value of J in point ip
!   intent(in)
    real(dp), dimension(:), allocatable :: J

!   plastic strain gammap.
!   gammap(ip) is the value of the plastic strain in point ip
!   intent(in)
    real(dp), dimension(:), allocatable :: gammap

!   Jacobian of the right-hand side of constitutive visco-elastic
!   models without mode coupling (see rhs in the description of
!   the rhs_viscoelastic routines).
!   drhs(ip,i,j,mode) is derivative of component i with respect to component j
!   in point ip for mode number mode
!   components similar to the conformation tensor c or contravariant
!   deformation tensor b.
!   intent(out)
    real(dp), dimension(:,:,:,:), allocatable :: drhs

!   Jacobian of the right-hand side of constitutive visco-elastic
!   models with mode coupling (see rhs in the description of
!   the rhs_viscoelastic routines).
!   drhs_mm(ip,i,mi,j,mj) is derivative of component i of mode mi with respect
!   to component j of mode mj in point ip
!   components similar to the conformation tensor c or contravariant
!   deformation tensor b.
!   intent(out)
    real(dp), dimension(:,:,:,:,:), allocatable :: drhs_mm

!   derivative of the right-hand side with respect to J
!   drhsdJ(ip,comp,mode) is component comp in point ip for mode number mode
!   components similar to the conformation tensor c or contravariant
!   deformation tensor b.
!   intent(out)
    real(dp), dimension(:,:,:), allocatable :: drhsdJ

!   derivative of the relaxation terms of right-hand side with respect to gammap
!   drhsdgammap(ip,comp,mode) is component comp in point ip for mode number mode
!   components similar to the conformation tensor c or contravariant
!   deformation tensor b.
!   intent(out)
    real(dp), dimension(:,:,:), allocatable :: drhsdgammap

  end type vemopt_t


! type definition for using optional input/output of data to the lower level
! routines (single point, single mode) of the constitutive viscoelastic model.
! This can be easily extended without adding more arguments to the routines.

  type vemmod_t

!   compute drlxmod or drlxmod_mm
!   NOTE: Either drlxmod or drlxmod_mm is computed, not both. Which one is
!         computed depends on the coupling between modes. Also drlxmod is
!         used as a temporary storage for computing drlxmod_mm, but is not
!         a final usable result in that case.
!   intent(in)
    logical :: compute_drlxmod = .false.

!   compute ddefmod
!   intent(in)
    logical :: compute_ddefmod = .false.

!   take dependence of the model on relative change in volume J into account.
!   intent(in)
    logical :: dep_J = .false.

!   compute the derivative of the right-hand side wrt to J
!   intent(in)
    logical :: compute_drhsdJ = .false.

!   take dependence of the model on the plastic strain gammap into account.
!   intent(in)
    logical :: dep_gammap = .false.

!   compute the derivative of the right-hand side wrt to gammap
!   intent(in)
    logical :: compute_drhsdgammap = .false.

!   actual arguments in the call to the vonmises_mm routine in the
!   rhs_adapted_lambda subroutine. The values of vmmode1/vmmode2 in vemopt
!   determine the values here. See the desciption in vemopt_t.
    integer :: vmmode1, vmmode2

!   Jacobian of deformation term for a single mode in a single point
!        .
!        c = deformation(L,c) - relaxation
!
!   or
!        .
!        b = deformation(L,b) - relaxation
!
!   ddefmod(i,j) is derivative of component i with respect to component j.
!   component sequence similar to the conformation tensor c or contravariant
!   deformation tensor b, respectively
!   intent(out)
    real(dp), dimension(:,:), allocatable :: ddefmod

!   Jacobian of relaxation term for a single mode in a single point
!        .
!        c =  L * c + c * L^T - relaxation
!
!   or
!        .
!        b =  L * b - relaxation
!
!   drlxmod(i,j) is derivative of component i with respect to component j.
!   component sequence similar to the conformation tensor c or contravariant
!   deformation tensor b, respectively
!   intent(out)
    real(dp), dimension(:,:), allocatable :: drlxmod

!   Jacobian of relaxation term for multi mode in a single point
!        .
!        c =  L * c + c * L^T - relaxation
!
!   or
!        .
!        b =  L * b - relaxation
!
!   drlxmod_mm(i,j,k) is derivative of component i with respect to component j
!   of mode k.
!   component sequence similar to the conformation tensor c or contravariant
!   deformation tensor b, respectively
!   intent(out)
    real(dp), dimension(:,:,:), allocatable :: drlxmod_mm

!   relative change in volume J.
!   intent(in)
    real(dp) :: J = 1._dp

!   derivative of the relaxation terms of right-hand side with respect to J
!   drlxmoddJ(comp) is component comp
!   components similar to the conformation tensor c or contravariant
!   deformation tensor b.
!   intent(out)
    real(dp), dimension(:), allocatable :: drlxmoddJ

!   plastic strain gammap.
!   intent(in)
    real(dp) :: gammap = 0._dp

!   derivative of the relaxation terms of right-hand side with respect to gammap
!   drlxmoddgammap(comp) is component comp
!   components similar to the conformation tensor c or contravariant
!   deformation tensor b.
!   intent(out)
    real(dp), dimension(:), allocatable :: drlxmoddgammap

!   conformation tensor c, contravariant deformation tensor b or s=logc for
!   use in the actual call to the vonmises_mm routine in the
!   rhs_adapted_lambda subroutine.
!   c(comp,mode), b(comp,mode) is component comp for mode number mode
!   intent(in)
    real(dp), dimension(:,:), allocatable :: c, b, s

  end type vemmod_t


! type definition for using optional input/output of data to the routines
! (multi-point) of the evolution of plastic strain for the constitutive
! viscoelastic model.
! This can be easily extended without adding more arguments to the routines.

  type vemopt_gammap_t

    logical :: created = .false. ! structure created?

!   take dependence of the model on relative change in volume J into account.
    logical :: dep_J = .false.

!   compute the derivative of the right-hand side with respect to the
!   conformation tensor without mode coupling.
    logical :: compute_drhs = .false.

!   compute the derivative of the right-hand side with respect to the
!   conformation tensor with mode coupling.
!   NOTE: if the model does not have mode coupling, drhs_mm wil contain
!         only the single mode part of drhs.
    logical :: compute_drhs_mm = .false.

!   compute the derivative of the right-hand side wrt to J
    logical :: compute_drhsdJ = .false.

!   compute the derivative of the right-hand side wrt to gammap
    logical :: compute_drhsdgammap = .false.

!   actual arguments in the call to the vonmises_mm routine in the
!   rhs_adapted_lambda subroutine. In this way the modes used for the
!   computation of the von Mises stress can be imposed by the user.
!   NOTE: vmmode1=0, vmmode2=0 is the default, which means that the actual
!         values used for vmmode1 and vmmode2 are taken from the actual
!         mode1, mode2 values in the routine that uses them.
!   NOTE: setting vmmode1/=mode1, vmmode2/=mode2 is possible and creates
!         special models, where vmmode1/vmmode2 can be even be outside
!         mode1:mode2 range. However, this will not be supported by the
!         flow modules for now.
!   intent(in)
    integer :: vmmode1 = 0, vmmode2 = 0

!   relative change in volume J.
!   J(ip) is the value of J in point ip
!   intent(in)
    real(dp), dimension(:), allocatable :: J

!   Jacobian of the right-hand side of the plastic strain equation with
!   with respect to the conformation tensor without mode coupling.
!   drhs(ip,j) is derivative with respect to component j of specified mode
!   in point ip
!   components similar to the conformation tensor c or contravariant
!   deformation tensor b.
!   intent(out)
    real(dp), dimension(:,:), allocatable :: drhs

!   Jacobian of the right-hand side of the plastic strain equation with
!   with respect to the conformation tensor with mode coupling.
!   drhs_mm(ip,j,mj) is derivative with respect to component j of
!   mode mj in point ip
!   components similar to the conformation tensor c or contravariant
!   deformation tensor b.
!   intent(out)
    real(dp), dimension(:,:,:), allocatable :: drhs_mm

!   derivative of the right-hand side with respect to J
!   drhsdJ(ip) in point ip
!   intent(out)
    real(dp), dimension(:), allocatable :: drhsdJ

!   derivative of the right-hand side with respect to gammap
!   drhsdgammap(ip) in point ip
!   deformation tensor b.
!   intent(out)
    real(dp), dimension(:), allocatable :: drhsdgammap

  end type vemopt_gammap_t


! global parameters

! if subtractlinear = .true. the linear relaxation term is extracted from the
! the right-hand side of the constitutive equation. This term should appear in
! the implicit terms (left-hand side).
  logical :: subtractlinear = .false.


! interface for generic create subroutine

  interface create
    module procedure create_viscoelastic_model, create_vemopt, &
      create_vemopt_gammap
  end interface create


! interface for generic delete subroutine

  interface delete
    module procedure delete_viscoelastic_model, delete_vemopt, &
      delete_vemopt_gammap
  end interface delete

! interface for generic check subroutine

  interface check
    module procedure check_vemopt
  end interface check


contains


! Create viscoelastic model in the structure vemodel

  subroutine create_viscoelastic_model ( model, vemodel, nmodes, flowtype, &
    bvariant, alam_model, vmglobal, alamJ_model, norlxmode, alam_gammap_model, &
    gammap_mode )

!   the model number: see vemodel_t for the possibilities
    integer, intent(in) :: model

!   struct defining the model parameters
    type(vemodel_t), intent(inout) :: vemodel

!   number of modes
!   default: 1
    integer, intent(in), optional :: nmodes

!   flow type: 0=2D
!              1=2D, Axisymmetrical
!              2=3D
!   default: 0
    integer, intent(in), optional :: flowtype

!   Variant of the b-equation (for b-formulation):
!        0=set to the default value (CDT)
!        1=CDT (Carrozza et al.)
!        2=symmetric (Balci et al.)
!        3=Cholesky (Vaithianathan & Collins)
!        4=Cholesky with log (Vaithianathan & Collins)
!   default=1
    integer, intent(in), optional :: bvariant

!   The adapted lambda model number: see vemodel_t for the possibilities
!   default = 0 (no adapted lambda model)
    integer, intent(in), optional :: alam_model

!   Use global von Mises equivalent stress in the adapted lambda model for
!   a multi-mode model
!   default = .false.
    logical, intent(in), optional :: vmglobal

!   The adapted lambda model number for the dependence on J: see vemodel_t
!   for the possibilities.
!   NOTE: This only applies if alam_model>1.
!   default = 0 (no adapted lambda model for J)
    integer, intent(in), optional :: alamJ_model

!   The mode selected will have the relaxation term set to zero.
!   default = 0 (no modes will have zero relaxation term)
    integer, intent(in), optional :: norlxmode

!   The adapted lambda model number for the dependence on gamma_p: see vemodel_t
!   for the possibilities.
!   NOTE: This only applies if alam_model>1.
!   default = 0 (no adapted lambda model for gamma_p)
    integer, intent(in), optional :: alam_gammap_model

!   The mode selected determines the plastic strain gamma_p.
!   Required if alam_gammap_model == 1.
!   default = 0 (means no mode set)
    integer, intent(in), optional :: gammap_mode

    integer, parameter :: MAXMODEL = 24, MAXALAMMODEL = 5, MAXALAMJMODEL = 1, &
      MAXALAMGAMMAPMODEL = 1

!   arrays that must be filled as follows:
!
!     ncompbpar(:,1) ncompb for each model when coorsys = 0
!     ncompbpar(:,2) ncompb for each model when coorsys = 1
!     ncompbpar(:,3) ncompb for each model when coorsys = 2
!
!   where ncompb is the number of components of the contravariant deformation
!   tensor.
!
!     ncompcpar(:,1) ncompc for each model when flowtype = 0
!     ncompcpar(:,2) ncompc for each model when flowtype = 1
!     ncompcpar(:,3) ncompc for each model when flowtype = 2
!
!   where ncompc is the number of components of the conformation tensor
!
!     numnlnpar(:) number of nonlinear parameters for each model
!
!     ixsipar(:) ixsi for each model
!
!     ncomptpar(:,1) ncompt for each model when flowtype = 0
!     ncomptpar(:,2) ncompt for each model when flowtype = 1
!     ncomptpar(:,3) ncompt for each model when flowtype = 2
!
!   where ncompt is the number of stress components
!
!     numalampar(:) number of adapted lambda parameters for each alam_model
!     numalamJpar(:) number of adapted lambda parameters for each alamJ_model
!     numalamgammappar(:) number of adapted lambda parameters for each
!                       alam_gammap_model
!
    integer :: ncompbpar(MAXMODEL,3), ncompcpar(MAXMODEL,3)
    integer :: ncomptpar(MAXMODEL,3)
    integer :: numnlnpar(MAXMODEL), ixsipar(MAXMODEL)
    integer :: numalampar(0:MAXALAMMODEL), numalamJpar(0:MAXALAMJMODEL), &
               numalamgammappar(0:MAXALAMGAMMAPMODEL)


    if ( vemodel%created ) then
      write(*,'(/a/)') &
        'Error in create_viscoelastic_model: vemodel has already been created '
      stop
    end if

!   fill par arrays

    ncompbpar(:,1) = [ 4,4,4,4,4, 4,4,4,5,4, 4,4,4, 6,5, 5, 6,5,4,5, 4,5,5,5 ]
    ncompbpar(:,2) = [ 5,5,5,5,5, 5,5,5,5,5, 5,5,5, 6,5, 5, 6,5,5,5, 5,5,5,5 ]
    ncompbpar(:,3) = [ 9,9,9,9,9, 9,9,9,9,9, 9,9,9,10,9, 9,10,9,9,9, 9,9,9,9 ]
    ncompcpar(:,1) = [ 3,3,3,3,3, 3,3,3,4,3, 3,3,3, 5,4, 4, 5,4,3,4, 3,4,4,4 ]
    ncompcpar(:,2) = [ 4,4,4,4,4, 4,4,4,4,4, 4,4,4, 5,4, 4, 5,4,4,4, 4,4,4,4 ]
    ncompcpar(:,3) = [ 6,6,6,6,6, 6,6,6,6,6, 6,6,6, 7,6, 6, 7,6,6,6, 6,6,6,6 ]
    numnlnpar(:)   = [ 0,0,1,1,1, 1,2,1,1,1, 2,2,1, 2,3, 2, 2,2,2,1, 1,2,0,0 ]
    ixsipar(:)     = [ 0,0,0,0,0, 0,0,0,0,1, 2,2,0, 0,0, 0, 0,0,0,0, 0,0,0,0 ]
    numalampar(0)  = 0
    numalampar(1:) = [ 0,1,3,1,3 ]
    numalamJpar(0)  = 0
    numalamJpar(1:) = [ 1 ]
    numalamgammappar(0)  = 0
    numalamgammappar(1:) = [ 4 ]

    ncomptpar = ncompcpar
!   XPP  (models 14 & 17 only ) has different number of stress components
!   than for tensor c (one less):
    ncomptpar(14,:) = [ 4,4,6 ]
    ncomptpar(17,:) = [ 4,4,6 ]

!   set model

    if ( model < 1 .or. model > MAXMODEL ) then
      write(*,'(/a,i0/)') &
        'Error in create_viscoelastic_model: invalid model = ', model
      stop
    else
      vemodel%model = model
    end if

!   set nmodes

    if ( present(nmodes) ) then
      vemodel%nmodes = nmodes
    else
      vemodel%nmodes = 1
    end if

!   set flowtype

    if ( present(flowtype) ) then
      if ( flowtype < 0 .or. flowtype > 2 ) then
        write(*,'(/a,i0/)') &
          'Error in create_viscoelastic_model: invalid flowtype = ', flowtype
        stop
      end if
      vemodel%flowtype = flowtype
    else
      vemodel%flowtype = 0
    end if

!   set bvariant

    if ( present(bvariant) ) then
      if ( bvariant == 0 ) then
        vemodel%bvariant = 1  ! set to CDT
      else if ( bvariant < 1 .or. bvariant > 4 ) then
        write(*,'(/a,i0/)') &
          'Error in create_viscoelastic_model: invalid bvariant = ', bvariant
        stop
      end if
      vemodel%bvariant = bvariant
    else
      vemodel%bvariant = 1  ! CDT is default
    end if

!   set norlxmode

    if ( present(norlxmode) ) then
      if ( norlxmode < 0 .or. norlxmode > vemodel%nmodes ) then
        write(*,'(/a,i0/)') &
          'Error in create_viscoelastic_model: invalid norlxmode = ', norlxmode
        stop
      end if
      vemodel%norlxmode = norlxmode
    else
      vemodel%norlxmode = 0  ! none is default
    end if

!   set gammap_mode

    if ( present(gammap_mode) ) then
      if ( gammap_mode < 0 .or. &
           gammap_mode > vemodel%nmodes ) then
        write(*,'(/2a,i0/)') &
          'Error in create_viscoelastic_model:', &
          ' invalid gammap_mode = ', gammap_mode
        stop
      end if
      vemodel%gammap_mode = gammap_mode
    else
      vemodel%gammap_mode = 0  ! no mode is set
    end if

!   set alam_model

    if ( present(alam_model) ) then
      if ( alam_model == 0 ) then
        vemodel%alam_model = 0  ! no adapted lambda model
        vemodel%vmglobal = .false.
      else if ( any ( vemodel%model == [ 9,14,15,16,17,18,19,22 ] ) ) then
        write(*,'(/a/a,i0,a/)') &
          'Error in create_viscoelastic_model:', '  model = ', &
          vemodel%model, ' is incompatible with adapted lambda models.'
        stop
      else if ( alam_model < 1 .or. alam_model > MAXALAMMODEL ) then
        write(*,'(/a,i0/)') &
          'Error in create_viscoelastic_model: invalid alam_model = ', &
          alam_model
        stop
      else
        vemodel%alam_model = alam_model
        vemodel%vmglobal = set_optional ( variable=vmglobal, default=.false.)
        if ( vemodel%vmglobal .and. &
               all ( alam_model /= [ 2,3,4,5 ] ) ) then
          write(*,'(/a/a,i0,a/)') &
            'Error in create_viscoelastic_model:', ' alam_model = ', &
             alam_model, ' is incompatible with vmglobal = .true. '
          stop
        end if
      end if
    else
      vemodel%alam_model = 0  ! default: no adapted lambda model
      vemodel%vmglobal = .false.
    end if

!   set alamJ_model

    if ( present(alamJ_model) ) then
      if ( alamJ_model == 0 ) then
        vemodel%alamJ_model = 0  ! no adapted lambda model for J
      else if ( alamJ_model < 1 .or. alamJ_model > MAXALAMJMODEL ) then
        write(*,'(/a,i0/)') &
          'Error in create_viscoelastic_model: invalid alamJ_model = ', &
          alamJ_model
        stop
      else if ( vemodel%alam_model <= 1 ) then
        write(*,'(/a,i0/a,i0/)') &
          'Error in create_viscoelastic_model: alam_model = ', &
          vemodel%alam_model, &
          ' must be > 1 for alamJ_model = ', alamJ_model
        stop
      else
        vemodel%alamJ_model = alamJ_model
      end if
    else
      vemodel%alamJ_model = 0  ! default: no adapted lambda model for J
    end if

!   set alam_gammap_model

    if ( present(alam_gammap_model) ) then
      if ( alam_gammap_model == 0 ) then
        vemodel%alam_gammap_model = 0  ! no model set for gamma_p
      else if ( all ( vemodel%model /= [ 23,24 ] ) ) then
        write(*,'(/a/a,i0,a/)') &
          'Error in create_viscoelastic_model:', '  model = ', &
          vemodel%model, ' is incompatible with plastic strain softening'
        stop
      else if ( alam_gammap_model < 1 .or. &
                alam_gammap_model > MAXALAMGAMMAPMODEL ) then
        write(*,'(/a,i0/)') &
          'Error in create_viscoelastic_model: invalid alam_gammap_model = ', &
          alam_gammap_model
        stop
      else
        vemodel%alam_gammap_model = alam_gammap_model
        if ( vemodel%gammap_mode == 0 .and. &
                 any(vemodel%alam_gammap_model==[1]) ) then
          write(*,'(/a/a,i0,a/)') &
            'Error in create_viscoelastic_model:', ' gammap_mode = ', &
             gammap_mode, ' not set '
          stop
        end if
      end if
    else
      vemodel%alam_gammap_model = 0  ! default: no model set for gamma_p
    end if

!   set deviatoric

    vemodel%deviatoric = any ( vemodel%model == [ 24 ] )

!   set compressible

    vemodel%compressible = any ( vemodel%model == [ 24 ] ) &
                                 .or. vemodel%alamJ_model > 0

!   set plastic_strain

    vemodel%gammap = vemodel%alam_gammap_model > 0

!   set some model parameters

    vemodel%ncompb = ncompbpar(model,vemodel%flowtype+1)
    vemodel%ncompc = ncompcpar(model,vemodel%flowtype+1)
    vemodel%ncompt = ncomptpar(model,vemodel%flowtype+1)
    vemodel%ixsi   = ixsipar(model)

!   allocate material parameters

    allocate(vemodel%modulus(vemodel%nmodes))
    allocate(vemodel%lambda(vemodel%nmodes))
    allocate(vemodel%nonlin(numnlnpar(model),vemodel%nmodes))
    allocate(vemodel%alam(numalampar(vemodel%alam_model),vemodel%nmodes))
    allocate(vemodel%alamJ(numalamJpar(vemodel%alamJ_model),vemodel%nmodes))
    allocate(vemodel%alam_gammap(numalamgammappar(vemodel%alam_gammap_model),&
                                                         &vemodel%nmodes))

!   set material parameters to zero (must be filled by the user)

    vemodel%modulus = 0
    vemodel%lambda = 0
    vemodel%nonlin = 0
    vemodel%alam = 0
    vemodel%alamJ = 0
    vemodel%alam_gammap = 0

    vemodel%created = .true.

  end subroutine create_viscoelastic_model


! Delete viscoelastic model in the structure vemodel

  subroutine delete_viscoelastic_model ( vemodel )

    type(vemodel_t), intent(inout) :: vemodel


    if ( .not. vemodel%created ) then
      write(*,'(/2(a/))') &
        'Error in delete_viscoelastic_model:', &
        ' vemodel has not been created and cannot be deleted '
      stop
    end if

    vemodel%deviatoric = .false.
    vemodel%compressible = .false.
    vemodel%vmglobal = .false.
    vemodel%gammap = .false.
    vemodel%model = 0
    vemodel%nmodes = 0
    vemodel%flowtype = 0
    vemodel%bvariant = 0
    vemodel%ncompb = 0
    vemodel%ncompc = 0
    vemodel%ncompt = 0
    vemodel%ixsi   = 0
    vemodel%norlxmode   = 0
    vemodel%alam_model = 0
    vemodel%alamJ_model = 0
    vemodel%alam_gammap_model = 0

!   deallocate material parameters

    deallocate(vemodel%modulus)
    deallocate(vemodel%lambda)
    deallocate(vemodel%nonlin)
    deallocate(vemodel%alam)
    deallocate(vemodel%alamJ)
    deallocate(vemodel%alam_gammap)

    vemodel%created = .false.

  end subroutine delete_viscoelastic_model


! Check the structure vemopt

  subroutine check_vemopt ( vemodel, vemopt, name_of_routine, components )

    type(vemodel_t), intent(in) :: vemodel
    type(vemopt_t), intent(in), optional :: vemopt
    character(len=*), intent(in) :: name_of_routine

!   which component to check?
!   possibilities: 'dep_J', 'dep_gammap', 'drhsdJ', 'nodrhsdJ'
!   NOTE: when specifying with an array constructor, all lengths must the same
    character(len=*), dimension(:), intent(in) :: components


    if ( present(vemopt) ) then

!     check vemopt

      if ( vemodel%vmglobal ) then
        if ( any([vemopt%vmmode1,vemopt%vmmode2]<0) .or. &
           any([vemopt%vmmode1,vemopt%vmmode2]>vemodel%nmodes) ) then
          write(*,'(/3a/a/)') 'Error in ', name_of_routine,':', &
            ' vmmode1 and/or vmmode2 component in vemopt out of range'
          stop
        end if
      end if

      if ( any(components=='dep_J') ) then

!       check dep_J

        if ( vemopt%dep_J ) then

          if ( .not. vemodel%compressible ) then
            write(*,'(/3a/2(a,i0)/)') 'Error in ', name_of_routine,':', &
              ' vemopt%depJ=.true. not applicable to model = ', &
              vemodel%model, ' with alamJ_model = ', vemodel%alamJ_model
            stop
          endif

        else

          if ( vemodel%compressible ) then
            write(*,'(/3a/2(a,i0)/)') 'Error in ', name_of_routine,':', &
              ' vemopt%dep_J=.false. whereas vemodel%compressible=.true.', &
              ' model = ', vemodel%model, &
              ' with alamJ_model = ', vemodel%alamJ_model
            stop
          end if

        end if

      end if

      if ( any(components=='dep_gammap') ) then

!       check dep_gammap

        if ( vemopt%dep_gammap ) then

          if ( .not. vemodel%gammap ) then
            write(*,'(/3a/2(a,i0)/)') 'Error in ', name_of_routine,':', &
              ' vemopt%dep_gammap=.true. not applicable to model = ', &
              vemodel%model, &
              ' with alam_gammap_model = ', vemodel%alam_gammap_model
            stop
          endif

        else

          if ( vemodel%gammap ) then
            write(*,'(/3a/2(a,i0)/)') 'Error in ', name_of_routine,':', &
              ' vemopt%dep_gammap=.false. whereas vemodel%gammap=.true.', &
              ' model = ', vemodel%model, &
              ' with alam_gammap_model = ', vemodel%alam_gammap_model
            stop
          end if

        end if

      end if

!     check compute_drhsd

      if ( any(components=='drhsdJ') ) then

        if ( vemopt%compute_drhsdJ ) then

          if ( .not. vemodel%compressible ) then
            write(*,'(/3a/2(a,i0)/)') 'Error in ', name_of_routine,':', &
              ' vemopt%compute_drhsdJ=.true. not applicable to model = ', &
              vemodel%model, ' with alamJ_model = ', vemodel%alamJ_model
            stop
          endif

        end if

      end if

!     check compute_drhsd if it should not be false

      if ( any(components=='nodrhsdJ') ) then

!       check compute_drhsdJ

        if ( vemopt%compute_drhsdJ ) then
          write(*,'(/3a/a/)') 'Error in ', name_of_routine,':', &
            ' vemopt%compute_drhsdJ=.true. not applicable.'
          stop
        end if

      end if

    else if ( vemodel%compressible .or. vemodel%gammap ) then

      write(*,'(/3a/a,i0/)') 'Error in ', name_of_routine,':', &
        ' vemopt must be present if vemodel%compressible=.true.'
      stop

    end if

  end subroutine check_vemopt


! create structure vemopt.

  subroutine create_vemopt ( vemopt, dep_J, dep_gammap, compute_drhs, &
    compute_drhs_mm, compute_drhsdJ, compute_drhsdgammap, &
    np, ncomp, nmodes, vmmode1, vmmode2 )

    type(vemopt_t), intent(inout) :: vemopt

!   take dependence of the model on relative change in volume J into account.
!   default=.false.
    logical, intent(in), optional :: dep_J

!   take dependence of the model on plastic strain softeming into account.
!   default=.false.
    logical, intent(in), optional :: dep_gammap

!   compute the single mode derivative of the right-hand side wrt to
!   conformation tensor
!   default=.false.
    logical, intent(in), optional :: compute_drhs

!   compute the multi-mode derivative of the right-hand side wrt to
!   conformation tensor
!   default=.false.
    logical, intent(in), optional :: compute_drhs_mm

!   compute the derivative of the right-hand side wrt to J
!   default=.false.
    logical, intent(in), optional :: compute_drhsdJ

!   compute the derivative of the right-hand side wrt to gammap
    logical, intent(in), optional :: compute_drhsdgammap

!   number of modes
!   required if compute_drhs=.true. or compute_drhs_mm=.true.
!   or compute_drhsdJ=.true.
    integer, intent(in), optional :: nmodes

!   number of points
!   default: 1
    integer, intent(in), optional :: np

!   number of components of the rhs of the CE
!   required if compute_drhs=.true. or compute_drhs_mm=.true.
!   or compute_drhsdJ=.true.
    integer, intent(in), optional :: ncomp

!   actual arguments in the call to the vonmises_mm routine in the
!   rhs_adapted_lambda subroutine. In this way the modes used for the
!   computation if the von Mises stress can be imposed by the user.
!   See desciption in the type vemopt_t.
    integer, intent(in), optional :: vmmode1, vmmode2

    integer :: lnp

    if ( vemopt%created ) then
      write(*,'(/a/)') &
        'Error in create_vemopt: vemopt has already been created '
      stop
    end if

    vemopt%dep_J = set_optional ( variable=dep_J, default=.false. )
    vemopt%dep_gammap = set_optional ( variable=dep_gammap, default=.false. )
    vemopt%compute_drhs = &
                 set_optional ( variable=compute_drhs, default=.false. )
    vemopt%compute_drhs_mm = &
                 set_optional ( variable=compute_drhs_mm, default=.false. )
    vemopt%compute_drhsdJ = &
                 set_optional ( variable=compute_drhsdJ, default=.false. )
    vemopt%compute_drhsdgammap = &
                 set_optional ( variable=compute_drhsdgammap, default=.false. )
    lnp = set_optional ( variable=np, default=1 )

    if ( vemopt%dep_J ) allocate(vemopt%J(lnp))

    if ( vemopt%compute_drhs .or. vemopt%compute_drhs_mm .or. &
         vemopt%compute_drhsdJ .or. vemopt%compute_drhsdgammap ) then
      if ( .not. present(ncomp) ) then
        write(*,'(/a/)') &
          'Error in create_vemopt: ncomp argument missing.'
        stop
      end if
      if ( .not. present(nmodes) ) then
        write(*,'(/a/)') &
          'Error in create_vemopt: nmodes argument missing.'
        stop
      end if
    end if

    if ( vemopt%compute_drhs ) allocate(vemopt%drhs(lnp,ncomp,ncomp,nmodes))

    if ( vemopt%compute_drhs_mm ) &
                       allocate(vemopt%drhs_mm(lnp,ncomp,nmodes,ncomp,nmodes))

    if ( vemopt%compute_drhsdJ ) allocate(vemopt%drhsdJ(lnp,ncomp,nmodes))
    if ( vemopt%compute_drhsdgammap ) &
                       allocate(vemopt%drhsdgammap(lnp,ncomp,nmodes))

    if ( present(vmmode1) .and. present(vmmode2) ) then
      vemopt%vmmode1 = vmmode1
      vemopt%vmmode2 = vmmode2
    else if ( present(vmmode1) .or. present(vmmode2) ) then
       write(*,'(/a/)') &
        'Error in create_vemopt: vmmode1 or vmmode2 argument missing.'
      stop
    else
      vemopt%vmmode1 = 0
      vemopt%vmmode2 = 0
    end if

    vemopt%created = .true.

  end subroutine create_vemopt


! Delete single vemopt

  subroutine delete_single_vemopt ( vemopt, nr )

    type(vemopt_t), intent(inout) :: vemopt
    integer, intent(in) :: nr

    if ( .not. vemopt%created ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error in delete_vemopt:', &
        ' vemopt has not been created and cannot be deleted ', &
        ' vemopt number in heading = ', nr
      stop
    end if

!   deallocate all allocatables and intialize to defaults
    call deall ( vemopt )

  contains

    subroutine deall ( vemopt )
      type(vemopt_t), intent(out) :: vemopt
    end subroutine deall

  end subroutine delete_single_vemopt


! Delete system vector

  subroutine delete_vemopt ( vemopt1, vemopt2, vemopt3, &
    vemopt4, vemopt5 )

    type(vemopt_t), intent(inout) :: vemopt1
    type(vemopt_t), intent(inout), optional :: vemopt2, vemopt3, &
      vemopt4, vemopt5

    call delete_single_vemopt(vemopt1,1)
    if ( present(vemopt2) ) call delete_single_vemopt(vemopt2,2)
    if ( present(vemopt3) ) call delete_single_vemopt(vemopt3,3)
    if ( present(vemopt4) ) call delete_single_vemopt(vemopt4,4)
    if ( present(vemopt5) ) call delete_single_vemopt(vemopt5,5)

  end subroutine delete_vemopt


! create structure vemopt_gammap.

  subroutine create_vemopt_gammap ( vemopt_gammap, dep_J, compute_drhs, &
    compute_drhs_mm, compute_drhsdJ, compute_drhsdgammap, &
    np, ncomp, nmodes, vmmode1, vmmode2 )

    type(vemopt_gammap_t), intent(inout) :: vemopt_gammap

!   take dependence of the model on relative change in volume J into account.
!   default=.false.
    logical, intent(in), optional :: dep_J

!   compute the single mode derivative of the right-hand side wrt to
!   conformation tensor
!   default=.false.
    logical, intent(in), optional :: compute_drhs

!   compute the multi-mode derivative of the right-hand side wrt to
!   conformation tensor
!   default=.false.
    logical, intent(in), optional :: compute_drhs_mm

!   compute the derivative of the right-hand side wrt to J
!   default=.false.
    logical, intent(in), optional :: compute_drhsdJ

!   compute the derivative of the right-hand side wrt to gammap
    logical, intent(in), optional :: compute_drhsdgammap

!   number of modes
!   required if compute_drhs=.true.
    integer, intent(in), optional :: nmodes

!   number of points
!   default: 1
    integer, intent(in), optional :: np

!   number of components of the conformation tensor
!   required if compute_drhs=.true.
    integer, intent(in), optional :: ncomp

!   actual arguments in the call to the vonmises_mm routine in the
!   rhs_adapted_lambda subroutine. In this way the modes used for the
!   computation if the von Mises stress can be imposed by the user.
!   See desciption in the type vemopt_gammap_t.
    integer, intent(in), optional :: vmmode1, vmmode2

    integer :: lnp

    if ( vemopt_gammap%created ) then
      write(*,'(/a/)') &
        'Error in create_vemopt_gammap: vemopt_gammap has already been created '
      stop
    end if

    vemopt_gammap%dep_J = set_optional ( variable=dep_J, default=.false. )
    vemopt_gammap%compute_drhs = &
                 set_optional ( variable=compute_drhs, default=.false. )
    vemopt_gammap%compute_drhs_mm = &
                 set_optional ( variable=compute_drhs_mm, default=.false. )
    vemopt_gammap%compute_drhsdJ = &
                 set_optional ( variable=compute_drhsdJ, default=.false. )
    vemopt_gammap%compute_drhsdgammap = &
                 set_optional ( variable=compute_drhsdgammap, default=.false. )
    lnp = set_optional ( variable=np, default=1 )

    if ( vemopt_gammap%dep_J ) allocate(vemopt_gammap%J(lnp))

    if ( vemopt_gammap%compute_drhs .or. vemopt_gammap%compute_drhs_mm ) then
      if ( .not. present(ncomp) ) then
        write(*,'(/a/)') &
          'Error in create_vemopt_gammap: ncomp argument missing.'
        stop
      end if
    end if

    if ( vemopt_gammap%compute_drhs_mm ) then
      if ( .not. present(nmodes) ) then
        write(*,'(/a/)') &
          'Error in create_vemopt_gammap: nmodes argument missing.'
        stop
      end if
    end if

    if ( vemopt_gammap%compute_drhs ) &
                       allocate(vemopt_gammap%drhs(lnp,ncomp))
    if ( vemopt_gammap%compute_drhs_mm ) &
                       allocate(vemopt_gammap%drhs_mm(lnp,ncomp,nmodes))

    if ( vemopt_gammap%compute_drhsdJ ) allocate(vemopt_gammap%drhsdJ(lnp))
    if ( vemopt_gammap%compute_drhsdgammap ) &
                       allocate(vemopt_gammap%drhsdgammap(lnp))

    if ( present(vmmode1) .and. present(vmmode2) ) then
      vemopt_gammap%vmmode1 = vmmode1
      vemopt_gammap%vmmode2 = vmmode2
    else if ( present(vmmode1) .or. present(vmmode2) ) then
       write(*,'(/a/)') &
        'Error in create_vemopt_gammap: vmmode1 or vmmode2 argument missing.'
      stop
    else
      vemopt_gammap%vmmode1 = 0
      vemopt_gammap%vmmode2 = 0
    end if

    vemopt_gammap%created = .true.

  end subroutine create_vemopt_gammap


! Delete single vemopt_gammap

  subroutine delete_single_vemopt_gammap ( vemopt_gammap, nr )

    type(vemopt_gammap_t), intent(inout) :: vemopt_gammap
    integer, intent(in) :: nr

    if ( .not. vemopt_gammap%created ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error in delete_vemopt_gammap:', &
        ' vemopt_gammap has not been created and cannot be deleted ', &
        ' vemopt_gammap number in heading = ', nr
      stop
    end if

!   deallocate all allocatables and intialize to defaults
    call deall ( vemopt_gammap )

  contains

    subroutine deall ( vemopt_gammap )
      type(vemopt_gammap_t), intent(out) :: vemopt_gammap
    end subroutine deall

  end subroutine delete_single_vemopt_gammap


! Delete system vector

  subroutine delete_vemopt_gammap ( vemopt_gammap1, vemopt_gammap2, &
    vemopt_gammap3, vemopt_gammap4, vemopt_gammap5 )

    type(vemopt_gammap_t), intent(inout) :: vemopt_gammap1
    type(vemopt_gammap_t), intent(inout), optional :: vemopt_gammap2, &
      vemopt_gammap3, vemopt_gammap4, vemopt_gammap5

    call delete_single_vemopt_gammap(vemopt_gammap1,1)
    if ( present(vemopt_gammap2) ) &
                           call delete_single_vemopt_gammap(vemopt_gammap2,2)
    if ( present(vemopt_gammap3) ) &
                           call delete_single_vemopt_gammap(vemopt_gammap3,3)
    if ( present(vemopt_gammap4) ) &
                           call delete_single_vemopt_gammap(vemopt_gammap4,4)
    if ( present(vemopt_gammap5) ) &
                           call delete_single_vemopt_gammap(vemopt_gammap5,5)

  end subroutine delete_vemopt_gammap

end module viscoelastic_models_defs_m
