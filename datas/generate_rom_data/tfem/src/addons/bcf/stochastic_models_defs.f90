
! Copyright (C) 2006-2006 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


!
! definition for stochastic models
!

module stochastic_models_defs_m

  use kind_defs_m

  implicit none


! type definition of a stochastic model

  type stmodel_t

    logical :: created = .false. ! structure created?
    integer :: ncompq = 0  ! number of components of the Q vector
    integer :: ncompt = 0  ! number of components of the stress tensor
    integer :: coorsys = 0 ! coordinate system: 0=Cartesian planar
                           !                    1=Axisymmetrical
                           !                    2=Cartesian, 3D
    integer :: model = 0   ! the model number:
!
!           modelnr   name                      nonlinear parameters
!
!              1      Hookean dumbbell (2D,axi)
!              2      FENE-P                    b
!              3      FENE                      b

!   material parameters determining the linear spectrum:
!     modulus(mode):  G modulus parameter for mode
!     lambda(mode):  relaxation time parameter for mode
    real(dp), dimension(:), allocatable :: modulus, lambda

!   nonlinear material parameters
!     nonlin(i,mode) the ith nonlinear parameter for mode
    real(dp), dimension(:,:), allocatable :: nonlin

  end type stmodel_t


! type definition for the numerical parameters of a stochastic model

  type stnumpar_t

    integer :: nfield = 0  ! number of fields

!   Time integration scheme for the Q-equation:
!     1: explicit Euler
!     2: backward Euler (implicit) for the force part
!     3: Crank-Nicolson (implicit) for the force part (only available for FENE)
    integer :: timeint = 0

    real(dp) :: tstep = 0  ! time step

  end type stnumpar_t


! interface for generic create subroutine

  interface create
    module procedure create_stochastic_model
  end interface create


! interface for generic delete subroutine

  interface delete
    module procedure delete_stochastic_model
  end interface delete

contains


! Create stochastic model in the structure stmodel

  subroutine create_stochastic_model ( model, stmodel, coorsys )

!   the model number: see stmodel_t for the possibilities
    integer, intent(in) :: model

    type(stmodel_t), intent(inout) :: stmodel

!   coordinate system: 0=Cartesian planar
!                      1=Axisymmetrical
!                      2=Cartesian, 3D
!   default: 0
    integer, intent(in), optional :: coorsys


    integer, parameter :: MAXMODEL = 3

!   arrays that must be filled as follows:
!
!     ncompqpar(:,1) ncompq for each model when coorsys = 0
!     ncompqpar(:,2) ncompq for each model when coorsys = 1
!     ncompqpar(:,3) ncompq for each model when coorsys = 2
!
!   where ncompq is the number of components of the Q-vector
!
!     numnlnpar(:) number of nonlinear parameters for each model
!
!     ncomptpar(:,1) ncompt for each model when coorsys = 0
!     ncomptpar(:,2) ncompt for each model when coorsys = 1
!     ncomptpar(:,3) ncompt for each model when coorsys = 2
!
!   where ncompt is the number of stress components
    integer :: ncompqpar(MAXMODEL,3), ncomptpar(MAXMODEL,3)
    integer :: numnlnpar(MAXMODEL)


    if ( stmodel%created ) then
      write(*,'(/a/)') &
        'Error in create_stochastic_model: stmodel has already been created '
      stop
    end if

!   fill par arrays

    ncompqpar(:,1) = [ 2,3,3 ]
    ncompqpar(:,2) = [ 3,3,3 ]
    ncompqpar(:,3) = [ 3,3,3 ]
    numnlnpar(:)   = [ 0,1,1 ]
    ncomptpar(:,1) = [ 3,4,4 ]
    ncomptpar(:,2) = [ 4,4,4 ]
    ncomptpar(:,3) = [ 6,6,6 ]

!   set model

    if ( model < 1 .or. model > MAXMODEL ) then
      write(*,'(/a,i0/)') &
        'Error in create_stochastic_model: invalid model = ', model
      stop
    else
      stmodel%model = model
    end if

!   set coorsys

    if ( present(coorsys) ) then
      if ( coorsys < 0 .or. coorsys > 2 ) then
        write(*,'(/a,i0/)') &
          'Error in create_stochastic_model: invalid coorsys = ', coorsys
        stop
      end if
      stmodel%coorsys = coorsys
    else
      stmodel%coorsys = 0
    end if

!   set some model parameters

    stmodel%ncompq = ncompqpar(model,stmodel%coorsys+1)
    stmodel%ncompt = ncomptpar(model,stmodel%coorsys+1)

!   allocate material parameters (single mode for now)

    allocate(stmodel%modulus(1))
    allocate(stmodel%lambda(1))
    allocate(stmodel%nonlin(numnlnpar(model),1))

!   set material parameters to zero (must be filled by the user)

    stmodel%modulus = 0
    stmodel%lambda = 0
    stmodel%nonlin = 0

    stmodel%created = .true.

  end subroutine create_stochastic_model


! Delete stochastic model in the structure vemodel

  subroutine delete_stochastic_model ( stmodel )

    type(stmodel_t), intent(inout) :: stmodel


    if ( .not. stmodel%created ) then
      write(*,'(/2(a/))') &
        'Error in delete_stochastic_model:', &
        ' stmodel has not been created and cannot be deleted '
      stop
    end if

    stmodel%model = 0
    stmodel%coorsys = 0
    stmodel%ncompq = 0
    stmodel%ncompt = 0

!   deallocate material parameters

    deallocate(stmodel%modulus)
    deallocate(stmodel%lambda)
    deallocate(stmodel%nonlin)

    stmodel%created = .false.

  end subroutine delete_stochastic_model



end module stochastic_models_defs_m

