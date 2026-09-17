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
! Viscoelastic models in 2D and axisymmetrical coordinate systems
! Standard case
! Optional Jacobian matrix.
!

module viscoelastic_models_2D_m

  use glob_defs_m
  use viscoelastic_models_defs_m
  use viscoelastic_models_adapted_lambda_m
  use tensor_m
  use viscoelastic_models_plastic_strain_m

  implicit none


contains


! right-hand side for 2D and axisymmetric viscoelastic models.

  subroutine rhs_viscoelastic_2D ( vemodel, gradv, c, rhs, mode1, mode2, &
    mvemodel, drhs, vemopt )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(c,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   velocity gradient L
!   gradv(ip,comp) is component comp in point ip with
!     comp = 1 Lxx
!     comp = 2 Lxy
!     comp = 3 Lyx
!     comp = 4 Lyy
!     comp = 5 Ltt (axisymmetric hoop strain rate)
    real(dp), dimension(:,:), intent(in) :: gradv

!   conformation tensor c
!   c(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambdab for double equation XPP model
    real(dp), dimension(:,:,:), intent(in) :: c

!   right-hand side of constitutive visco-elastic models:
!        .
!        c =  L * c + c * L^T - relaxation.
!
!   rhs(ip,comp,mode) is component comp in point ip for mode number mode
!   components sequence similar to the conformation tensor c
    real(dp), dimension(:,:,:), intent(inout) :: rhs

!   Obsolete and gives an error message if present.
!   Argument will be removed in a future version.
    real(dp), dimension(:,:,:,:), intent(out), optional :: drhs

!   if present: limit the range of modes to mode1--mode2
    integer, intent(in), optional :: mode1, mode2

!   if present: additional optional parameters (see type description)
    type(vemopt_t), intent(inout), optional :: vemopt


    logical :: compute_drhs, compute_drhs_mm
    integer :: mode, ncompc, np, ip, modenr1, modenr2, m, vm1, vm2
    real(dp), dimension(vemodel%ncompc) :: defmod, rlxmod
    type(vemmod_t) :: vemmod


    if ( vemodel%model == 0 ) then

      write(*,'(/3(a/))') &
        'Error in rhs_viscoelastic_2D: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 0 .and. vemodel%flowtype /= 1 ) then

      write(*,'(/a/2a/)') &
        'Error in rhs_viscoelastic_2D: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=0 (2D) or ', &
        'flowtype=1 (2D, axisymmetric) '
      stop

    end if

    if ( present(mode1) .and. present(mode2) ) then
      if ( mode1 < 1 .or. mode1 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D:', &
          ' mode1 out of range '
        stop
      end if
      if ( mode2 < 1 .or. mode2 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D:', &
          ' mode2 out of range '
        stop
      end if
      modenr1 = mode1
      modenr2 = mode2
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(c,1) ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( present(drhs) ) then
      write(*,'(/4(a/))') 'Error in rhs_viscoelastic_2D:', &
        ' Argument drhs is obsolete and will be removed in the future.', &
        ' The output array drhs has been moved to vemopt. Use the routine ', &
        ' create_vemopt with argument compute_drhs=.true. instead. '
      stop
    end if

    call check ( vemodel, vemopt, 'rhs_viscoelastic_2D', &
      components=['dep_J     ', 'dep_gammap', 'drhsdJ    '] )

    compute_drhs = .false.
    compute_drhs_mm = .false.

    if ( present(vemopt) ) then

!     copy logicals from vemopt to vemmod
      vemmod%dep_J = vemopt%dep_J
      vemmod%compute_drhsdJ = vemopt%compute_drhsdJ
      vemmod%dep_gammap = vemopt%dep_gammap
      vemmod%compute_drhsdgammap = vemopt%compute_drhsdgammap

      if ( vemodel%vmglobal .and. &
                 any([vemopt%vmmode1,vemopt%vmmode2]==0) ) then
!       set vmmode to modenr
        vemmod%vmmode1 = modenr1
        vemmod%vmmode2 = modenr2
      else
!       set vmmode to values in vemopt
        vemmod%vmmode1 = vemopt%vmmode1
        vemmod%vmmode2 = vemopt%vmmode2
      end if

      if ( vemopt%compute_drhs ) then

        if ( vemodel%vmglobal ) then
          write(*,'(/a/a,i0/)') 'Error in rhs_viscoelastic_2D:', &
            ' single mode Jacobian matrix not available for multi-mode', &
            ' von-Mises stress models '
          stop
        endif
        if ( any( vemodel%model == [14,17] ) ) then
          write(*,'(/a/a,i0/)') 'Error in rhs_viscoelastic_2D:', &
            ' Jacobian matrix not available for model = ', vemodel%model
          stop
        endif
        if ( subtractlinear ) then
          write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D:', &
            ' Jacobian matrix not available for subtractlinear=.true.'
          stop
        endif

        vemmod%compute_ddefmod = .true.
        vemmod%compute_drlxmod = .true.

        compute_drhs = .true.

      end if

      if ( vemopt%compute_drhs_mm ) then

        if ( any( vemodel%model == [14,17] ) ) then
          write(*,'(/a/a,i0/)') 'Error in rhs_viscoelastic_2D:', &
            ' Jacobian matrix not available for model = ', vemodel%model
          stop
        endif
        if ( subtractlinear ) then
          write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D:', &
            ' Jacobian matrix not available for subtractlinear=.true.'
          stop
        endif

        vemmod%compute_ddefmod = .true.
        vemmod%compute_drlxmod = .true.

        compute_drhs_mm = .true.

        allocate(&
          vemmod%drlxmod_mm(vemodel%ncompc,vemodel%ncompc,vemodel%nmodes) )

      end if

      if ( vemopt%compute_drhsdgammap ) then

        allocate(vemmod%drlxmoddgammap(vemodel%ncompc) )

      end if

    else if ( vemodel%vmglobal ) then

!     set default for vmglobal=T: modes for von Mises equal to modenr1:modenr2
      vemmod%vmmode1 = modenr1
      vemmod%vmmode2 = modenr2

    end if

!   always allocate arrays in vemmod for the associate contruct

    allocate(vemmod%ddefmod(vemodel%ncompc,vemodel%ncompc))
    allocate(vemmod%drlxmod(vemodel%ncompc,vemodel%ncompc))
    allocate(vemmod%drlxmoddJ(vemodel%ncompc))

    if ( vemodel%vmglobal ) then
      allocate ( vemmod%c(size(c,2),size(c,3)) )
      vm1 = vemmod%vmmode1
      vm2 = vemmod%vmmode2
    end if

    ncompc = vemodel%ncompc
    np = size(c,1)

    do ip = 1, np

!     set J

      if ( vemmod%dep_J ) vemmod%J = vemopt%J(ip)

!     set gammap

      if ( vemmod%dep_gammap ) vemmod%gammap = vemopt%gammap(ip)

!     von Mises for multiple modes

      if ( vemodel%vmglobal ) vemmod%c(:,vm1:vm2) = c(ip,:,vm1:vm2)

      do mode = modenr1, modenr2

!       deformation term

        if ( present(mvemodel) ) then
          call rhs_deformation_2D ( mvemodel(ip), mode, gradv(ip,:), &
            c(ip,:,mode), defmod, vemmod )
        else
          call rhs_deformation_2D ( vemodel, mode, gradv(ip,:), &
            c(ip,:,mode), defmod, vemmod )
        end if

!       relaxation term

        if ( present(mvemodel) ) then
          call rhs_relaxation_2D ( mvemodel(ip), mode, c(ip,:,mode), rlxmod, &
            vemmod )
        else
          call rhs_relaxation_2D ( vemodel, mode, c(ip,:,mode), rlxmod, vemmod )
        end if

!       fill right-hand side

        rhs(ip,:ncompc,mode) = defmod - rlxmod

!       fill derivative (single mode) of right-hand side

        if ( compute_drhs ) then
          vemopt%drhs(ip,:ncompc,:ncompc,mode) = &
                                            vemmod%ddefmod - vemmod%drlxmod
        end if

!       fill derivative (multi mode) of right-hand side

        if ( compute_drhs_mm ) then

          if ( vemodel%vmglobal ) then

!           mode coupling due to global von Mises stress

!           diagonal
            vemopt%drhs_mm(ip,:ncompc,mode,:ncompc,mode) = &
                                  vemmod%ddefmod - vemmod%drlxmod_mm(:,:,mode)

!           off-diagonal
            do m = vemmod%vmmode1, vemmod%vmmode2
              if ( m == mode ) cycle ! diagonal already accounted for
              vemopt%drhs_mm(ip,:ncompc,mode,:ncompc,m) = &
                                                 - vemmod%drlxmod_mm(:,:,m)
            end do

          else

!           no mode coupling: diagonal drhs

            vemopt%drhs_mm(ip,:ncompc,mode,:ncompc,modenr1:mode-1) = 0
            vemopt%drhs_mm(ip,:ncompc,mode,:ncompc,mode) = &
                                            vemmod%ddefmod - vemmod%drlxmod
            vemopt%drhs_mm(ip,:ncompc,mode,:ncompc,mode+1:modenr2) = 0

          end if

        end if

!       fill derivative of right-hand side with respect to J

        if ( vemmod%compute_drhsdJ ) then
          vemopt%drhsdJ(ip,:ncompc,mode) = - vemmod%drlxmoddJ
        end if

!       fill derivative of right-hand side with respect to gammap

        if ( vemmod%compute_drhsdgammap ) then
          vemopt%drhsdgammap(ip,:ncompc,mode) = - vemmod%drlxmoddgammap
        end if

      end do

    end do

  end subroutine rhs_viscoelastic_2D


! right-hand side for 2D and axisymmetric viscoelastic models
! (only relaxation terms)

  subroutine rhs_viscoelastic_relax_2D ( vemodel, c, rhs, mode1, mode2, &
    mvemodel, drhs, vemopt )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(c,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   conformation tensor c
!   c(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambdab for double equation XPP model
    real(dp), dimension(:,:,:), intent(in) :: c

!   relaxation terms only of the right-hand side of constitutive
!   visco-elastic models:
!        .
!        c =  L * c + c * L^T - relaxation.
!
!   rhs(ip,comp,mode) is component comp in point ip for mode number mode
!   components similar to the conformation tensor c
    real(dp), dimension(:,:,:), intent(inout) :: rhs

!   Obsolete and gives an error message if present.
!   Argument will be removed in a future version.
    real(dp), dimension(:,:,:,:), intent(out), optional :: drhs

!   if present: limit the range of modes to mode1--mode2
    integer, intent(in), optional :: mode1, mode2

!   if present: additional optional parameters (see type description)
    type(vemopt_t), intent(inout), optional :: vemopt


    logical :: compute_drhs, compute_drhs_mm
    integer :: modenr1, modenr2, mode, ncompc, np, ip, m, vm1, vm2
    real(dp), dimension(vemodel%ncompc) :: rlxmod
    type(vemmod_t) :: vemmod


    if ( vemodel%model == 0 ) then

      write(*,'(/3(a/))') &
        'Error in rhs_viscoelastic_relax_2D: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 0 .and. vemodel%flowtype /= 1 ) then

      write(*,'(/a/2a/)') &
        'Error in rhs_viscoelastic_relax_2D: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=0 (2D) or ', &
        'flowtype=1 (2D, axisymmetric) '
      stop

    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(c,1) ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_relax_2D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_relax_2D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( present(drhs) ) then
      write(*,'(/4(a/))') 'Error in rhs_viscoelastic_relax_2D:', &
        ' Argument drhs is obsolete and will be removed in the future.', &
        ' The output array drhs has been moved to vemopt. Use the routine ', &
        ' create_vemopt with argument compute_drhs=.true. instead. '
      stop
    end if

    call check ( vemodel, vemopt, 'rhs_viscoelastic_relax_2D', &
      components=['dep_J ', 'drhsdJ'] )

    if ( present(mode1) .and. present(mode2) ) then
      modenr1 = mode1
      modenr2 = mode2
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    compute_drhs = .false.
    compute_drhs_mm = .false.

    if ( present(vemopt) ) then

!     copy logicals from vemopt to vemmod
      vemmod%dep_J = vemopt%dep_J
      vemmod%compute_drhsdJ = vemopt%compute_drhsdJ
      vemmod%dep_gammap = vemopt%dep_gammap
      vemmod%compute_drhsdgammap = vemopt%compute_drhsdgammap

      if ( vemodel%vmglobal .and. &
                 any([vemopt%vmmode1,vemopt%vmmode2]==0) ) then
!       set vmmode to modenr
        vemmod%vmmode1 = modenr1
        vemmod%vmmode2 = modenr2
      else
!       set vmmode to values in vemopt
        vemmod%vmmode1 = vemopt%vmmode1
        vemmod%vmmode2 = vemopt%vmmode2
      end if

      if ( vemopt%compute_drhs ) then
        vemmod%compute_drlxmod = .true.
        compute_drhs = .true.
      end if

      if ( vemopt%compute_drhs_mm ) then

        vemmod%compute_drlxmod = .true.
        compute_drhs_mm = .true.

        allocate(&
          vemmod%drlxmod_mm(vemodel%ncompc,vemodel%ncompc,vemodel%nmodes) )

      end if

      if ( vemopt%compute_drhsdgammap ) then

        allocate(vemmod%drlxmoddgammap(vemodel%ncompc) )

      end if

    else if ( vemodel%vmglobal ) then

!     set default for vmglobal=T: modes for von Mises equal to modenr1:modenr2
      vemmod%vmmode1 = modenr1
      vemmod%vmmode2 = modenr2

    end if

!   always allocate arrays in vemmod for the associate contruct

    allocate(vemmod%drlxmod(vemodel%ncompc,vemodel%ncompc))
    allocate(vemmod%drlxmoddJ(vemodel%ncompc))

    if ( vemodel%vmglobal ) then
      allocate ( vemmod%c(size(c,2),size(c,3)) )
      vm1 = vemmod%vmmode1
      vm2 = vemmod%vmmode2
    end if

    ncompc = vemodel%ncompc
    np = size(c,1)

    do ip = 1, np

!     set J

      if ( vemmod%dep_J ) vemmod%J = vemopt%J(ip)

!     set gammap

      if ( vemmod%dep_gammap ) vemmod%gammap = vemopt%gammap(ip)

!     von Mises for multiple modes

      if ( vemodel%vmglobal ) vemmod%c(:,vm1:vm2) = c(ip,:,vm1:vm2)

      do mode = modenr1, modenr2

!       relaxation terms

        if ( present(mvemodel) ) then
          call rhs_relaxation_2D ( mvemodel(ip), mode, c(ip,:,mode), rlxmod, &
            vemmod )
        else
          call rhs_relaxation_2D ( vemodel, mode, c(ip,:,mode), rlxmod, vemmod )
        end if

!       fill right-hand side

        rhs(ip,:ncompc,mode) = - rlxmod

!       fill derivative (single mode) of right-hand side

        if ( compute_drhs ) then
          vemopt%drhs(ip,:ncompc,:ncompc,mode) = - vemmod%drlxmod
        end if

!       fill derivative (multi mode) of right-hand side

        if ( compute_drhs_mm ) then

          if ( vemodel%vmglobal ) then

!           mode coupling due to global von Mises stress

!           diagonal
            vemopt%drhs_mm(ip,:ncompc,mode,:ncompc,mode) = &
                                                 - vemmod%drlxmod_mm(:,:,mode)

!           off-diagonal
            do m = vemmod%vmmode1, vemmod%vmmode2
              if ( m == mode ) cycle ! diagonal already accounted for
              vemopt%drhs_mm(ip,:ncompc,mode,:ncompc,m) = &
                                                 - vemmod%drlxmod_mm(:,:,m)
            end do

          else

!           no mode coupling: diagonal drhs

            vemopt%drhs_mm(ip,:ncompc,mode,:ncompc,modenr1:mode-1) = 0
            vemopt%drhs_mm(ip,:ncompc,mode,:ncompc,mode) = - vemmod%drlxmod
            vemopt%drhs_mm(ip,:ncompc,mode,:ncompc,mode+1:modenr2) = 0

          end if

        end if

!       fill derivative of right-hand side with respect to J

        if ( vemmod%compute_drhsdJ ) then
          vemopt%drhsdJ(ip,:ncompc,mode) = - vemmod%drlxmoddJ
        end if

!       fill derivative of right-hand side with respect to gammap

        if ( vemmod%compute_drhsdgammap ) then
          vemopt%drhsdgammap(ip,:ncompc,mode) = - vemmod%drlxmoddgammap
        end if

      end do

    end do

  end subroutine rhs_viscoelastic_relax_2D


! deformation terms for 3D viscoelastic models.

  subroutine rhs_deformation_2D ( vemodel, mode, gradv, c, defmod, vemmod )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   the mode number
    integer, intent(in) :: mode

!   velocity gradient L
!   gradv(comp) is component comp
!     comp = 1 Lxx
!     comp = 2 Lxy
!     comp = 3 Lyx
!     comp = 4 Lyy
!     comp = 5 Ltt (axisymmetric hoop strain rate)
    real(dp), dimension(:), intent(in) :: gradv

!   conformation tensor c
!   c(comp) is component comp
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambdab for double equation XPP
    real(dp), dimension(:), intent(in) :: c

!   deformation term for a single mode in a single point
!        .
!        c =  deformation(c,L) - relaxation.
!
!   defmod(comp) is component comp
!   components similar to the conformation tensor c
    real(dp), dimension(:), intent(out) :: defmod

!   additional optional parameters (see type description)
    type(vemmod_t), intent(inout) :: vemmod

    logical :: deviatoric
    integer :: flowtype, ncompc
    real(dp) :: cxx, cxy, cyy, czz, Lxx, Lxy, Lyx, Lyy, Lzz, xsi
    real(dp) :: Lsxx, Lsxy, Lsyx, Lsyy, Lszz, ddss
    real(dp) :: trL3
    real(dp), parameter :: z=0._dp


    associate ( ddefmod => vemmod%ddefmod )


    ncompc = vemodel%ncompc
    flowtype = vemodel%flowtype
    deviatoric = vemodel%deviatoric

!   extract c

    cxx = c(1)
    cxy = c(2)
    cyy = c(3)

    if ( flowtype == 1 .or. deviatoric &
                       .or. any(vemodel%model==[14,17]) ) then
      czz = c(4)
    end if

!   extract L

    Lxx = gradv(1)
    Lxy = gradv(2)
    Lyx = gradv(3)
    Lyy = gradv(4)
    if ( flowtype == 1 ) Lzz  = gradv(5)  ! axi-symmetric

!   slip term

    if ( vemodel%ixsi > 0  ) then

!     models with slip

      xsi = vemodel%nonlin(vemodel%ixsi,mode)

      Lsxx = Lxx - xsi * Lxx
      Lsxy = Lxy - xsi * ( Lxy + Lyx ) / 2
      Lsyx = Lyx - xsi * ( Lxy + Lyx ) / 2
      Lsyy = Lyy - xsi * Lyy
      if ( flowtype == 1 ) then
        Lszz = Lzz - xsi * Lzz  ! axi-symmetric
      else
        Lszz = 0
      end if

    else

      Lsxx = Lxx
      Lsxy = Lxy
      Lsyx = Lyx
      Lsyy = Lyy
      if ( flowtype == 1 ) then
        Lszz = Lzz ! axi-symmetric
      else
        Lszz = 0
      end if

    end if

!   deviatoric

    if ( deviatoric ) then

!     replace L with L^d = L - (tr L)/3 I

      trL3 = ( Lsxx + Lsyy + Lszz ) / 3
      Lsxx = Lsxx - trL3
      Lsyy = Lsyy - trL3
      Lszz = Lszz - trL3

    end if

!   Quasi-linear terms

    defmod(1) = 2._dp * ( cxx * Lsxx + cxy * Lsxy )
    defmod(2) = cxx * Lsyx + cxy * ( Lsxx + Lsyy ) + cyy * Lsxy
    defmod(3) = 2._dp * ( cxy * Lsyx + cyy * Lsyy )

    if ( flowtype == 1 .or. deviatoric ) then
      defmod(4) = 2._dp * czz * Lszz
      defmod(5:ncompc) = 0._dp
    else
      defmod(4:ncompc) = 0._dp
    end if

    if ( vemmod%compute_ddefmod ) then

!     Jacobian of quasi-linear terms

      ddefmod(1,1:3) = 2 * [ Lsxx,      Lsxy,    z ]
      ddefmod(2,1:3) =     [ Lsyx, Lsxx+Lsyy, Lsxy ]
      ddefmod(3,1:3) = 2 * [    z,      Lsyx, Lsyy ]

      if ( flowtype == 1 .or. deviatoric ) then
        ddefmod(1:3,4) = 0
        ddefmod(4,:) = [ z, z, z, 2 * Lszz ]
        ddefmod(5:ncompc,:) = 0
        ddefmod(:,5:ncompc) = 0
      else
        ddefmod(4:ncompc,:) = 0
        ddefmod(:,4:ncompc) = 0
      end if

    end if

    if ( vemodel%model == 14 .or. vemodel%model == 17 ) then

!     XPP double equation model

      ddss  = cxx*Lsxx + cyy*Lsyy + cxy*(Lsxy+Lsyx)

      defmod(1) = defmod(1) - 2._dp*ddss*cxx
      defmod(2) = defmod(2) - 2._dp*ddss*cxy
      defmod(3) = defmod(3) - 2._dp*ddss*cyy
      defmod(4) = defmod(4) - 2._dp*ddss*czz

!     Check if backbonestretch is in logconformation (model=17)

      if ( vemodel%model == 14 ) then
        defmod(5) = defmod(5) + ddss*c(5)
      else if ( vemodel%model == 17 ) then
        defmod(5) = defmod(5) + ddss
      end if

    end if


    end associate

  end subroutine rhs_deformation_2D


! relaxation term in the models

  subroutine rhs_relaxation_2D ( vemodel, mode, c, rlxmod, vemmod )

    use math_defs_m
    use limits_m, only: USE_EGP_STRESS_TENSOR_FORM, &
                        USE_ADAP_LAMBDA_FOR_EGP_DET_STAB


    type(vemodel_t), intent(in) :: vemodel

!   the mode number
    integer, intent(in) :: mode

!   conformation tensor c
!   c(comp) is component comp
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambdab for double equation XPP

    real(dp), dimension(:), intent(in) :: c

!   relaxation term for a single mode in a single point
!        .
!        c =  L * c + c * L^T - relaxation.
!
!   rlxmod(comp) is component comp
!   components similar to the conformation tensor c
    real(dp), dimension(:), intent(out) :: rlxmod

!   additional optional parameters (see type description)
    type(vemmod_t), intent(inout) :: vemmod


    integer :: flowtype

    real(dp) :: lambda, alpha, i1, i2, i1mi2, lambd2, fac, beta, eps, yfac
    real(dp) :: fi1, gi1pi2, yfacgl, f, feq, bpar, trc, xsi, epsxsi, gamma
    real(dp) :: rL2, lambdas, nu, relmod, expfac, delta, fretr, fccr, trc3
    real(dp) :: cm(2,2), cm2(2,2), sqtrc, taur, taud, lambdab, eta, tau_y, mu
    real(dp) :: trtau, tau_e, chpar, fm2(2,2)
    real(dp) :: tau(vemodel%ncompt)


    if ( mode == vemodel%norlxmode ) then

!     no relaxation term: set relaxation term to zero
      call set_rlx_to_zero

!     add EGP determinant stabilization even if relaxation is set to zero
      call add_egp_det_stabilization

      return

    end if


    flowtype = vemodel%flowtype

    select case ( vemodel%model )

    case(1)

!     Leonov

      lambda = vemodel%lambda(mode)
      rlxmod(1) = c(1) ** 2 + c(2) ** 2 - 1
      rlxmod(3) = c(3) ** 2 + c(2) ** 2 - 1
      if ( flowtype == 1 ) then
        i1 = c(1) + c(3) + c(4)
        i2 = c(1) * c(3) + c(3) * c(4) + c(4) * c(1) - c(2) ** 2
        i1mi2 = ( i1 - i2 ) / 3
        rlxmod(1) = rlxmod(1) - i1mi2 * c(1)
        rlxmod(3) = rlxmod(3) - i1mi2 * c(3)
      end if
      lambd2 = 2 * lambda
      rlxmod(1) = rlxmod(1) / lambd2
      rlxmod(2) = ( c(1) + c(3) )  * c(2)
      rlxmod(3) = rlxmod(3) / lambd2
      if ( flowtype == 1 ) then
        rlxmod(2) = rlxmod(2) - i1mi2 * c(2)
        rlxmod(4) = ( c(4) ** 2 - 1 - i1mi2 * c(4) ) / lambd2
      end if
      rlxmod(2) = rlxmod(2) / lambd2

    case(2)

!     Maxwell/Oldroyd

      lambda = vemodel%lambda(mode)
      rlxmod(1) = ( c(1) - 1 ) / lambda
      rlxmod(2) = c(2) / lambda
      rlxmod(3) = ( c(3) - 1 ) / lambda
      if ( flowtype == 1 ) then
        rlxmod(4) = ( c(4) - 1 ) / lambda
      end if

    case(3)

!     Giesekus

      lambda = vemodel%lambda(mode)
      alpha  = vemodel%nonlin(1,mode)
      rlxmod(1) = ( c(1) - 1 + alpha * ( ( c(1) - 1 ) ** 2 +  &
                  c(2) ** 2 ) ) / lambda
      rlxmod(2) = ( c(2) + alpha * ( c(1) + c(3) - 2 ) * c(2) ) / lambda
      rlxmod(3) = ( c(3) - 1 + alpha * ( ( c(3) - 1 ) ** 2 +  &
                  c(2) ** 2 ) ) / lambda
      if ( flowtype == 1 ) then
        rlxmod(4) = ( c(4) - 1 + alpha * ( c(4) - 1 ) ** 2 ) / lambda
      end if

    case(4)

!     Larson

      lambda = vemodel%lambda(mode)
      gamma  = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
        fac = 1 + gamma * ( c(1) + c(3) - 2 )
      else if ( flowtype == 1 ) then
        fac = 1 + gamma * ( c(1) + c(3) + c(4) - 3 )
      end if
      fac = fac / lambda

      rlxmod(1) = fac * ( c(1) - 1 )
      rlxmod(2) = fac * c(2)
      rlxmod(3) = fac * ( c(3) - 1 )
      if ( flowtype == 1 ) then
         rlxmod(4) = fac * ( c(4) - 1 )
      end if

    case(5,6)

!     Phan-Thien/Tanner

      lambda = vemodel%lambda(mode)
      eps    = vemodel%nonlin(1,mode)

      if ( vemodel%model == 5 ) then

!       linear factor

        if ( flowtype == 0 ) then
           yfac = 1 + eps * ( c(1) + c(3) - 2 )
        else if ( flowtype == 1 ) then
           yfac = 1 + eps * ( c(1) + c(3) + c(4) - 3 )
        end if

      else if ( vemodel%model == 6 ) then

!       exponential factor

        if ( flowtype == 0 ) then
          yfac = exp ( eps * ( c(1) + c(3) - 2 ) )
        else if ( flowtype == 1 ) then
          yfac = exp ( eps * ( c(1) + c(3) + c(4) - 3 ) )
        end if

      end if

      yfacgl = yfac / lambda

      rlxmod(1) = yfacgl * ( c(1) - 1 )
      rlxmod(2) = yfacgl * c(2)
      rlxmod(3) = yfacgl * ( c(3) - 1 )
      if ( flowtype == 1 ) then
         rlxmod(4) = yfacgl * ( c(4) - 1 )
      end if

    case(7)

!     Extended Leonov

      lambda = vemodel%lambda(mode)
      alpha  = vemodel%nonlin(1,mode)
      beta   = vemodel%nonlin(2,mode)

      if ( flowtype == 0 ) then

!       planar

        i1   = c(1) + c(3)
        fi1  = 1._dp / ( 1 + 2 * alpha / pi * &
                        atan( beta / 2 * ( i1 - 2 ) )  )
        lambd2 = 2 * lambda
        fac    = fi1 / lambd2

        rlxmod(1) = fac * ( c(1) ** 2 + c(2) ** 2 - 1 )
        rlxmod(2) = fac * ( c(1) + c(3) )  * c(2)
        rlxmod(3) = fac * ( c(3) ** 2 + c(2) ** 2 - 1 )

      else if ( flowtype == 1 ) then

!       axisymmetrical

        i1 = c(1) + c(3) + c(4)
        i2 = c(1) * c(3) + c(3) * c(4) + c(4) * c(1) - c(2) ** 2
        gi1pi2  = 1._dp / (  1 + 2 * alpha / pi *  &
                        atan( beta / 4 * ( i1 + i2 - 6 ) )  )
        i1mi2  = ( i1 - i2 ) / 3
        lambd2 = 2 * lambda
        fac    = gi1pi2 / lambd2

        rlxmod(1) = fac * ( c(1) ** 2 + c(2) ** 2 - 1 - i1mi2 * c(1) )
        rlxmod(2) = fac * ( ( c(1) + c(3) ) * c(2) - i1mi2 * c(2) )
        rlxmod(3) = fac * ( c(3) ** 2 + c(2) ** 2 - 1 - i1mi2 * c(3) )
        rlxmod(4) = fac * ( c(4) ** 2 - 1 - i1mi2 * c(4) )

      end if

    case(8)

!     Chilcott-Rallison

      lambda = vemodel%lambda(mode)
      rL2    = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
         trc = c(1) + c(3) + 1
      else if ( flowtype == 1 ) then
         trc = c(1) + c(3) + c(4)
      end if
      f = 1 / ( 1 - trc / rL2 )
      fac = f / lambda

      rlxmod(1) = fac * ( c(1) - 1 )
      rlxmod(2) = fac * c(2)
      rlxmod(3) = fac * ( c(3) - 1 )
      if ( flowtype == 1 ) then
         rlxmod(4) = fac * ( c(4) - 1 )
      end if

    case(9)

!     FENE-P

      lambda = vemodel%lambda(mode)
      bpar   = vemodel%nonlin(1,mode)

      trc = c(1) + c(3) + c(4)
      f   = 1 / ( 1 - trc / ( bpar + 3 ) )
      feq = ( bpar + 3 ) / bpar

      rlxmod(1) = ( f * c(1) - feq ) / lambda
      rlxmod(2) = f * c(2) / lambda
      rlxmod(3) = ( f * c(3) - feq )  / lambda
      rlxmod(4) = ( f * c(4) - feq ) / lambda

    case(10)

!     Johnson-Segalman

      lambda = vemodel%lambda(mode)
      rlxmod(1) = ( c(1) - 1 ) / lambda
      rlxmod(2) = c(2) / lambda
      rlxmod(3) = ( c(3) - 1 ) / lambda
      if ( flowtype == 1 ) then
         rlxmod(4) = ( c(4) - 1 ) / lambda
      end if

    case(11,12)

!     Phan-Thien/Tanner with slip

      lambda = vemodel%lambda(mode)
      eps    = vemodel%nonlin(1,mode)
      xsi    = vemodel%nonlin(2,mode)
      epsxsi = eps / ( 1 - xsi )

      if ( vemodel%model == 11 ) then

!       linear factor

        if ( flowtype == 0 ) then
           yfac = 1 + epsxsi * ( c(1) + c(3) - 2 )
        else if ( flowtype == 1 ) then
           yfac = 1 + epsxsi * ( c(1) + c(3) + c(4) - 3 )
        end if

      else if ( vemodel%model == 12 ) then

!       exponential factor

        if ( flowtype == 0 ) then
           yfac = exp ( epsxsi * ( c(1) + c(3) - 2 ) )
        else if ( flowtype == 1 ) then
           yfac = exp ( epsxsi * ( c(1) + c(3) + c(4) - 3 ) )
        end if

      end if

      yfacgl = yfac / lambda

      rlxmod(1) = yfacgl * ( c(1) - 1 )
      rlxmod(2) = yfacgl * c(2)
      rlxmod(3) = yfacgl * ( c(3) - 1 )
      if ( flowtype == 1 ) then
         rlxmod(4) = yfacgl * ( c(4) - 1 )
      end if

    case(13)

!     Marrucci

      lambda = vemodel%lambda(mode)
      beta   = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
         yfac = 1._dp / ( 1 - beta * ( c(1) + c(3) - 2 ) )
      else if ( flowtype == 1 ) then
         yfac = 1._dp / ( 1 - beta * ( c(1) + c(3) + c(4) - 3 ) )
      end if

      yfacgl = yfac / lambda

      rlxmod(1) = yfacgl * ( c(1) - 1 )
      rlxmod(2) = yfacgl * c(2)
      rlxmod(3) = yfacgl * ( c(3) - 1 )
      if ( flowtype == 1 ) then
         rlxmod(4) = yfacgl * ( c(4) - 1 )
      end if

    case(14,17)

!     XPP double equation
!     14 -> c(5) = lambdab
!     17 -> exp(c(5)) = lambdab

      lambda  = vemodel%lambda(mode)
      lambdas = vemodel%nonlin(1,mode)
      nu      = vemodel%nonlin(2,mode)

      if ( vemodel%model == 14 ) then
        lambdab  = c(5)
      else if ( vemodel%model == 17 ) then
!       backbonestretch is in log conformation
        lambdab  = exp(c(5))
      end if

      relmod  = 1._dp / ( 3._dp * lambda * lambdab**2 )
      expfac  = exp ( nu * ( lambdab - 1._dp ) )

      rlxmod(1) = relmod * ( 3._dp * c(1) - 1._dp )
      rlxmod(2) = relmod * ( 3._dp * c(2) )
      rlxmod(3) = relmod * ( 3._dp * c(3) - 1._dp )
      rlxmod(4) = relmod * ( 3._dp * c(4) - 1._dp )

      if ( vemodel%model == 14 ) then
        rlxmod(5) = expfac / lambdas * ( lambdab - 1._dp/lambdab )
      else if ( vemodel%model == 17 ) then
!       if backbonestretch is in log conformation
        rlxmod(5) = expfac / lambdas * ( 1._dp - 1._dp/lambdab**2 )
      end if

    case(15)

!     Rolie-Poly

      taud     = vemodel%lambda(mode)
      taur     = vemodel%nonlin(1,mode)
      beta     = vemodel%nonlin(2,mode)
      delta    = vemodel%nonlin(3,mode)

      trc      = c(1) + c(3) + c(4)
      sqtrc    = sqrt ( 3._dp / trc )
      fretr    = (1._dp - sqtrc) / taur
      fccr     = beta * ( (trc/3._dp) ** delta )

      rlxmod(1) = (c(1) - 1._dp)/taud + fretr * (c(1) + fccr * ( c(1) - 1._dp ))
      rlxmod(2) = (c(2)        )/taud + fretr * (c(2) + fccr *   c(2))
      rlxmod(3) = (c(3) - 1._dp)/taud + fretr * (c(3) + fccr * ( c(3) - 1._dp ))
      rlxmod(4) = (c(4) - 1._dp)/taud + fretr * (c(4) + fccr * ( c(4) - 1._dp ))

    case(16)

!     XPP single equation

      lambda    = vemodel%lambda(mode)
      lambdas   = vemodel%nonlin(1,mode)
      nu        = vemodel%nonlin(2,mode)

      trc       = c(1) + c(3) + c(4)
      trc3      = 3._dp / trc
      expfac    = 2._dp * exp ( nu * ( sqrt(trc/3._dp) - 1._dp ) )
      expfac    = expfac * (1._dp - trc3) / lambdas

      rlxmod(1) = expfac * c(1) + ( trc3 * c(1) - 1._dp ) / lambda
      rlxmod(2) = expfac * c(2) + ( trc3 * c(2) ) / lambda
      rlxmod(3) = expfac * c(3) + ( trc3 * c(3) - 1._dp ) / lambda
      rlxmod(4) = expfac * c(4) + ( trc3 * c(4) - 1._dp ) / lambda

    case(18,19)

!     PTT-XPP single equation
!       model=18: four components in 2D
!       model=19: three components in 2D, czz=1

      lambda    = vemodel%lambda(mode)
      lambdas   = vemodel%nonlin(1,mode)
      nu        = vemodel%nonlin(2,mode)

      if ( vemodel%model == 19 .and. flowtype == 0 ) then
        trc = c(1) + c(3) + 1._dp
      else
        trc = c(1) + c(3) + c(4)
      end if
      trc3   = 3._dp / trc
      expfac = 2._dp * exp ( nu * ( sqrt(trc/3._dp) - 1._dp ) )
      expfac = expfac * (1._dp - trc3) / lambdas

      rlxmod(1) = expfac * c(1) - expfac + ( trc3 * c(1) - trc3 ) / lambda
      rlxmod(2) = expfac * c(2) + ( trc3 * c(2) ) / lambda
      rlxmod(3) = expfac * c(3) - expfac + ( trc3 * c(3) - trc3 ) / lambda
      if ( vemodel%model == 18 .or. flowtype == 1 ) then
        rlxmod(4) = expfac * c(4) - expfac + ( trc3 * c(4) - trc3 ) / lambda
      end if

    case(20)

!     FENE-P (see Wapperom & Hulsen 1998, Huetter, Hulsen & Anderson 2018)

      lambda = vemodel%lambda(mode)
      bpar   = vemodel%nonlin(1,mode)

      trc = c(1) + c(3) + c(4)
      f   = bpar / ( bpar + 3 - trc )

      rlxmod(1) = ( f * c(1) - 1 ) / lambda
      rlxmod(2) = f * c(2) / lambda
      rlxmod(3) = ( f * c(3) - 1 )  / lambda
      rlxmod(4) = ( f * c(4) - 1 ) / lambda

    case(21)

!     Chilcott-Rallison (alternative formulation)

      lambda = vemodel%lambda(mode)
      rL2    = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
         trc = c(1) + c(3) + 1
      else if ( flowtype == 1 ) then
         trc = c(1) + c(3) + c(4)
      end if
      f = ( rl2 - 3 ) / ( rL2 - trc )
      fac = f / lambda

      rlxmod(1) = fac * ( c(1) - 1 )
      rlxmod(2) = fac * c(2)
      rlxmod(3) = fac * ( c(3) - 1 )
      if ( flowtype == 1 ) then
         rlxmod(4) = fac * ( c(4) - 1 )
      end if

    case(22)

!     Saramito elastoviscoplastic (Drucker-Prager)

      eta = vemodel%modulus(mode) * vemodel%lambda(mode)
      tau_y = vemodel%nonlin(1,mode)
      mu = vemodel%nonlin(2,mode)

      call stress_viscoelastic_2D_single_mode ( vemodel, c, tau, mode )

      trtau = tau(1) + tau(3) + tau(4)  ! trace tau
      tau_e = vonmises_2D ( tau )       ! von Mises equivalent shear stress

!     find the regime the stress state is in

      chpar = mu * trtau - 3 * tau_y

      if ( chpar <= - 3 * tau_e ) then

!       regime I: sticking

        rlxmod = 0

      else if ( chpar >= 2 * mu**2 * tau_e ) then

!       regime III: loosing contact

        fac = 1 / eta
        f = tau_y / mu

        rlxmod(1) = fac * ( tau(1) - f )
        rlxmod(2) = fac * tau(2)
        rlxmod(3) = fac * ( tau(3) - f )
        rlxmod(4) = fac * ( tau(4) - f )

      else

!       regime II: sliding

        fac = ( tau_e - tau_y + mu * trtau / 3 ) &
                    / ( ( 1 + 2*mu**2/3 ) * tau_e * eta )
        f = trtau / 3 - 2 * mu * tau_e / 3

        rlxmod(1) = fac * ( tau(1) - f )
        rlxmod(2) = fac * tau(2)
        rlxmod(3) = fac * ( tau(3) - f )
        rlxmod(4) = fac * ( tau(4) - f )

      end if

    case(23,24)

!     EGP

      lambda = vemodel%lambda(mode)

      if ( USE_EGP_STRESS_TENSOR_FORM ) then

        eta = vemodel%modulus(mode) * lambda

        call stress_viscoelastic_2D_single_mode ( vemodel, c, tau, mode, &
          vemmod )

        fm2 = matmul ( vector_to_tensor2_symmetric (ndim=2,A=c(1:3)), &
                       vector_to_tensor2_symmetric (ndim=2,A=tau(1:3)) ) / eta

        rlxmod(1:3) = tensor2_to_vector_symmetric ( ndim=2, A=fm2 )
        rlxmod(4) = c(4) * tau(4) / eta

      else

        cm = vector_to_tensor2_symmetric (ndim=2,A=c(1:3))
        cm2 = matmul ( cm, cm )
        trc = c(1) + c(3) + c(4)
        rlxmod(1:3) = &
                 tensor2_to_vector_symmetric (ndim=2,A=cm2-trc/3*cm ) / lambda
        rlxmod(4) = ( c(4)**2 - trc/3*c(4) ) / lambda

       if ( vemodel%model == 24 ) rlxmod = rlxmod / vemmod%J

      end if

!     Add stabilization for determinant

      if ( USE_ADAP_LAMBDA_FOR_EGP_DET_STAB ) then
!       Add stabilization for determinant here
        rlxmod = rlxmod + ( detc_2D(c) - 1 ) / 3 * c / lambda
      end if

    case default

      write(*,'(a,i0)') &
        'Error in rhs_relaxation_2D: model not available: ', vemodel%model
      stop

    end select

    if ( vemmod%compute_drlxmod ) then

!     Jacobian of relaxation terms

      call drhs_relaxation_2D ( vemodel, mode, c, vemmod )

    end if

    if ( vemmod%compute_drhsdJ ) then

!     Derivative of relaxation terms with respect to J

      if ( vemodel%model == 24 ) then
        vemmod%drlxmoddJ = - rlxmod / vemmod%J
      else
        vemmod%drlxmoddJ = 0
      end if

    end if

    if ( vemodel%alam_model > 0 ) then

!     adapted lambda

      call rhs_adapted_lambda ( vemodel, mode, c, vonmises=vm_2D_sm, &
        vonmises_mm=vm_2D_mm, rlxmod=rlxmod, vemmod=vemmod )

    end if

    if ( .not. USE_ADAP_LAMBDA_FOR_EGP_DET_STAB ) then
!     Add stabilization for determinant here with original lambda constant
      call add_egp_det_stabilization
    end if

    if ( subtractlinear ) then

!     Subtract the linear part.

      lambda = vemodel%lambda(mode)
      rlxmod(1) = rlxmod(1) - ( c(1) - 1 ) / lambda
      rlxmod(2) = rlxmod(2) - c(2) / lambda
      rlxmod(3) = rlxmod(3) - ( c(3) - 1 ) / lambda
      if ( vemodel%ncompc == 4 ) then
         rlxmod(4) = rlxmod(4) - ( c(4) - 1 ) / lambda
      end if

    end if

  contains


!   no relaxation term: set relaxation term to zero

    subroutine set_rlx_to_zero

      rlxmod = 0
      if ( vemmod%compute_drlxmod ) then
        vemmod%drlxmod = 0
        if ( vemodel%vmglobal ) vemmod%drlxmod_mm = 0
      end if
      if ( vemmod%compute_drhsdJ ) vemmod%drlxmoddJ = 0
      if ( vemmod%compute_drhsdgammap ) vemmod%drlxmoddgammap = 0

    end subroutine set_rlx_to_zero


!   add EGP determinant stabilization term

    subroutine add_egp_det_stabilization

      real(dp) :: lambda
      integer :: i

      if ( all(vemodel%model/=[23,24]) ) return

      lambda = vemodel%lambda(mode)

 !    Add stabilization for determinant

      rlxmod = rlxmod + ( detc_2D(c) - 1 ) / 3 * c / lambda

      associate ( drlxmod => vemmod%drlxmod )

      if ( vemmod%compute_drlxmod ) then

!       Jacobian of stabilization term for determinant

        do i = 1, 4
          drlxmod(i,:) = drlxmod(i,:) + ddetc_2D(c) / 3 * c(i) / lambda
          drlxmod(i,i) = drlxmod(i,i) + ( detc_2D(c) - 1 ) / 3 / lambda
        end do

      end if

      end associate

    end subroutine add_egp_det_stabilization

  end subroutine rhs_relaxation_2D


! Jacobian of relaxation term in the models

  subroutine drhs_relaxation_2D ( vemodel, mode, c, vemmod )

    use math_defs_m
    use limits_m, only: USE_EGP_STRESS_TENSOR_FORM, &
                        USE_ADAP_LAMBDA_FOR_EGP_DET_STAB


    type(vemodel_t), intent(in) :: vemodel

!   the mode number
    integer, intent(in) :: mode

!   conformation tensor c
!   c(comp) is component comp
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambdab for double equation XPP
    real(dp), dimension(:), intent(in) :: c

!   additional optional parameters (see type description)
    type(vemmod_t), intent(inout) :: vemmod


    integer :: flowtype, i

    real(dp) :: lambda, alpha, fac, eps, yfac, beta
    real(dp) :: yfacgl, f, bpar, trc, xsi, epsxsi, gamma, rL2, eta, tau_y, mu
    real(dp) :: trtau, tau_e, chpar
    real(dp) :: tau(vemodel%ncompt), dtau(vemodel%ncompt,vemodel%ncompc)
    real(dp), dimension(vemodel%ncompc) :: dtrtau, dtau_e, s, dyfac, dyfacgl, &
      dfac, df
    real(dp) :: gxx, gxy, gyy, cminI(4)
    real(dp), parameter :: z=0._dp
    real(dp) :: cxx, cxy, cyy, czz, i1, i2, i1mi2, di1mi2(4), lambda2
    real(dp) :: taud, taur, sqtrc, rel1, rel2, rel3, rel4, rel5, lambdas, facc2
    real(dp) :: delta, dtrc(4), facdtrc(4), facc, lambdab, nu, trc3, expfac
    real(dp) :: txx, txy, tyy, tzz, w2(4,4)

    flowtype = vemodel%flowtype


    associate ( drlxmod => vemmod%drlxmod )


    select case ( vemodel%model )

    case(1)

!     Leonov

      lambda = vemodel%lambda(mode)

      cxx=c(1); cxy=c(2); cyy=c(3)

      if ( flowtype == 1 ) then

        czz=c(4)

!       c^2

        drlxmod(1,:) = 2 * [ cxx,     cxy,    z,    z ]
        drlxmod(2,:) =     [ cxy, cxx+cyy,  cxy,    z ]
        drlxmod(3,:) = 2 * [   z,     cxy,  cyy,    z ]
        drlxmod(4,:) = 2 * [   z,       z,    z,  czz ]

!       -(i1-i2)c/3

        i1 = cxx + cyy + czz
        i2 = cxx*cyy + cyy*czz + czz*cxx - cxy**2
        i1mi2 = ( i1 - i2 ) / 3
        di1mi2 = [1-cyy-czz,2*cxy,1-cxx-czz,1-cxx-cyy]/3

        do i = 1, 4
          drlxmod(i,:) = drlxmod(i,:) - di1mi2 * c(i)
          drlxmod(i,i) = drlxmod(i,i) - i1mi2
        end do

      else

!       c^2

        drlxmod(1,:) = 2 * [ cxx,     cxy,    z ]
        drlxmod(2,:) =     [ cxy, cxx+cyy,  cxy ]
        drlxmod(3,:) = 2 * [   z,     cxy,  cyy ]

      end if

      lambda2 = 2 * lambda

      drlxmod = drlxmod / lambda2

    case(2)

!     Maxwell/Oldroyd

      lambda = vemodel%lambda(mode)

      drlxmod = 0
      do i = 1, 3
        drlxmod(i,i) = 1 / lambda
      end do
      if ( flowtype == 1 ) then
        drlxmod(4,4) = 1 / lambda
      end if

    case(3)

!     Giesekus

      lambda = vemodel%lambda(mode)
      alpha  = vemodel%nonlin(1,mode)

!     s=c-I
      s([1,3]) = c([1,3]) - 1
      s(2) = c(2)
      if ( flowtype == 1 ) then
        s(4) = c(4) - 1
      end if

      drlxmod = 0
      do i = 1, 3
        drlxmod(i,i) = 1
      end do
      if ( flowtype == 1 ) then
        drlxmod(4,4) = 1
      end if

      gxx=s(1); gxy=s(2); gyy=s(3)

      drlxmod(1,1:3) = drlxmod(1,1:3) + 2 * alpha * [ gxx,     gxy,   z ]
      drlxmod(2,1:3) = drlxmod(2,1:3) +     alpha * [ gxy, gxx+gyy, gxy ]
      drlxmod(3,1:3) = drlxmod(3,1:3) + 2 * alpha * [   z,     gxy, gyy ]
      if ( flowtype == 1 ) then
        drlxmod(4,4) = drlxmod(4,4) + 2 * alpha * s(4)
      end if

      drlxmod = drlxmod / lambda

    case(4)

!     Larson

      lambda = vemodel%lambda(mode)
      gamma  = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
        fac = 1 + gamma * ( c(1) + c(3) - 2 )
        dfac = gamma * [1,0,1] / lambda
      else if ( flowtype == 1 ) then
        fac = 1 + gamma * ( c(1) + c(3) + c(4) - 3 )
        dfac = gamma * [1,0,1,1] / lambda
      end if
      fac = fac / lambda

      drlxmod(1,:) = dfac * ( c(1) - 1 )
      drlxmod(2,:) = dfac * c(2)
      drlxmod(3,:) = dfac * ( c(3) - 1 )
      do i = 1, 3
        drlxmod(i,i) = drlxmod(i,i) + fac
      end do
      if ( flowtype == 1 ) then
        drlxmod(4,:) = dfac * ( c(4) - 1 )
        drlxmod(4,4) = drlxmod(4,4) + fac
      end if

    case(5,6)

!     Phan-Thien/Tanner

      lambda = vemodel%lambda(mode)
      eps    = vemodel%nonlin(1,mode)

      if ( vemodel%model == 5 ) then

!       linear factor

        if ( flowtype == 0 ) then
          yfac = 1 + eps * ( c(1) + c(3) - 2 )
          dyfac = eps * [1,0,1]
        else if ( flowtype == 1 ) then
          yfac = 1 + eps * ( c(1) + c(3) + c(4) - 3 )
          dyfac = eps * [1,0,1,1]
        end if

      else if ( vemodel%model == 6 ) then

!       exponential factor

        if ( flowtype == 0 ) then
          yfac = exp ( eps * ( c(1) + c(3) - 2 ) )
          dyfac = eps * yfac * [1,0,1]
        else if ( flowtype == 1 ) then
          yfac = exp ( eps * ( c(1) + c(3) + c(4) - 3 ) )
          dyfac = eps * yfac * [1,0,1,1]
        end if

      end if

      yfacgl = yfac / lambda
      dyfacgl = dyfac / lambda

      drlxmod(1,:) = dyfacgl * ( c(1) - 1 )
      drlxmod(2,:) = dyfacgl * c(2)
      drlxmod(3,:) = dyfacgl * ( c(3) - 1 )
      do i = 1, 3
        drlxmod(i,i) = drlxmod(i,i) + yfacgl
      end do
      if ( flowtype == 1 ) then
        drlxmod(4,:) = dyfacgl * ( c(4) - 1 )
        drlxmod(4,4) = drlxmod(4,4) + yfacgl
      end if

    case(8)

!     Chilcott-Rallison

      lambda = vemodel%lambda(mode)
      rL2    = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
         trc = c(1) + c(3) + 1
      else if ( flowtype == 1 ) then
         trc = c(1) + c(3) + c(4)
      end if
      f = 1 / ( 1 - trc / rL2 )
      fac = f / lambda
      if ( flowtype == 0 ) then
        dfac = f**2 * [1,0,1] / ( rL2 * lambda )
      else if ( flowtype == 1 ) then
        dfac = f**2 * [1,0,1,1] / ( rL2 * lambda )
      end if

      drlxmod(1,:) = dfac * ( c(1) - 1 )
      drlxmod(2,:) = dfac * c(2)
      drlxmod(3,:) = dfac * ( c(3) - 1 )
      do i = 1, 3
        drlxmod(i,i) = drlxmod(i,i) + fac
      end do
      if ( flowtype == 1 ) then
        drlxmod(4,:) = dfac * ( c(4) - 1 )
        drlxmod(4,4) = drlxmod(4,4) + fac
      end if

    case(9)

!     FENE-P

      lambda = vemodel%lambda(mode)
      bpar   = vemodel%nonlin(1,mode)

      trc = c(1) + c(3) + c(4)
      f   = 1 / ( 1 - trc / ( bpar + 3 ) )
      df  = f**2 * [1,0,1,1] / ( bpar + 3 )

      drlxmod(1,:) = df * c(1) / lambda
      drlxmod(2,:) = df * c(2) / lambda
      drlxmod(3,:) = df * c(3) / lambda
      drlxmod(4,:) = df * c(4) / lambda
      do i = 1, 4
        drlxmod(i,i) = drlxmod(i,i) + f / lambda
      end do

    case(10)

!     Johnson-Segalman

      lambda = vemodel%lambda(mode)

      drlxmod = 0
      do i = 1, 3
        drlxmod(i,i) = 1 / lambda
      end do
      if ( flowtype == 1 ) then
        drlxmod(4,4) = 1 / lambda
      end if

    case(11,12)

!     Phan-Thien/Tanner with slip

      lambda = vemodel%lambda(mode)
      eps    = vemodel%nonlin(1,mode)
      xsi    = vemodel%nonlin(2,mode)
      epsxsi = eps / ( 1 - xsi )

      if ( vemodel%model == 11 ) then

!       linear factor

        if ( flowtype == 0 ) then
          yfac = 1 + epsxsi * ( c(1) + c(3) - 2 )
          dyfac = epsxsi * [1,0,1]
        else if ( flowtype == 1 ) then
          yfac = 1 + epsxsi * ( c(1) + c(3) + c(4) - 3 )
          dyfac = epsxsi * [1,0,1,1]
        end if

      else if ( vemodel%model == 12 ) then

!       exponential factor

        if ( flowtype == 0 ) then
          yfac = exp ( epsxsi * ( c(1) + c(3) - 2 ) )
          dyfac = epsxsi * yfac * [1,0,1]
        else if ( flowtype == 1 ) then
          yfac = exp ( epsxsi * ( c(1) + c(3) + c(4) - 3 ) )
          dyfac = epsxsi * yfac * [1,0,1,1]
        end if

      end if

      yfacgl = yfac / lambda
      dyfacgl = dyfac / lambda

      drlxmod(1,:) = dyfacgl * ( c(1) - 1 )
      drlxmod(2,:) = dyfacgl * c(2)
      drlxmod(3,:) = dyfacgl * ( c(3) - 1 )
      do i = 1, 3
        drlxmod(i,i) = drlxmod(i,i) + yfacgl
      end do
      if ( flowtype == 1 ) then
        drlxmod(4,:) = dyfacgl * ( c(4) - 1 )
        drlxmod(4,4) = drlxmod(4,4) + yfacgl
      end if
!
    case(13)

!     Marrucci

      lambda = vemodel%lambda(mode)
      beta   = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
        yfac = 1._dp / ( 1 - beta * ( c(1) + c(3) - 2 ) )
        dyfac = yfac**2 * beta * [1,0,1]
      else if ( flowtype == 1 ) then
        yfac = 1._dp / ( 1 - beta * ( c(1) + c(3) + c(4) - 3 ) )
        dyfac = yfac**2 * beta * [1,0,1,1]
      end if

      yfacgl = yfac / lambda
      dyfacgl = dyfac / lambda

      drlxmod(1,:) = dyfacgl * ( c(1) - 1 )
      drlxmod(2,:) = dyfacgl * c(2)
      drlxmod(3,:) = dyfacgl * ( c(3) - 1 )
      do i = 1, 3
        drlxmod(i,i) = drlxmod(i,i) + yfacgl
      end do
      if ( flowtype == 1 ) then
        drlxmod(4,:) = dyfacgl * ( c(4) - 1 )
        drlxmod(4,4) = drlxmod(4,4) + yfacgl
      end if

    case(15)

!     Rolie-Poly

      taud     = vemodel%lambda(mode)
      taur     = vemodel%nonlin(1,mode)
      beta     = vemodel%nonlin(2,mode)
      delta    = vemodel%nonlin(3,mode)

      trc   = c(1) + c(3) + c(4)
      sqtrc = sqrt( 3._dp / trc )
      rel1  = ( 1._dp - sqtrc ) / taur
      rel2  = beta * ( trc / 3._dp ) ** delta
      rel3  = sqrt( 3._dp ) / ( 2 * taur ) / trc ** 1.5_dp
      rel4  = beta * delta / 3 * ( trc / 3._dp ) ** ( delta - 1 )
      rel5  = rel1 * rel4

      facc = 1._dp/taud + rel1 * ( 1 + rel2 )
      facdtrc(1) = rel3 * ( c(1) + rel2 * (c(1)-1) ) + rel5 * (c(1)-1)
      facdtrc(2) = rel3 * ( c(2) + rel2 * c(2) ) + rel5 * c(2)
      facdtrc(3) = rel3 * ( c(3) + rel2 * (c(3)-1) ) + rel5 * (c(3)-1)
      facdtrc(4) = rel3 * ( c(4) + rel2 * (c(4)-1) ) + rel5 * (c(4)-1)
      dtrc = [1,0,1,1]

      drlxmod(1,:) = dtrc * facdtrc(1)
      drlxmod(2,:) = dtrc * facdtrc(2)
      drlxmod(3,:) = dtrc * facdtrc(3)
      drlxmod(4,:) = dtrc * facdtrc(4)
      do i = 1, 4
        drlxmod(i,i) = drlxmod(i,i) + facc
      end do

    case(16)

!     XPP single equation

      lambdab   = vemodel%lambda(mode)
      lambdas   = vemodel%nonlin(1,mode)
      nu        = vemodel%nonlin(2,mode)

      trc     = c(1) + c(3) + c(4)
      trc3    = 3._dp / trc
      expfac  = 2._dp * exp ( nu * ( sqrt(trc/3._dp) - 1._dp ) ) / lambdas
      facc    = expfac * (1._dp - trc3) + trc3 / lambdab
      facc2   = expfac * ( nu * (1._dp - trc3) / ( 2 * sqrt(3*trc ) ) &
                  + 3/trc**2 ) - 3/trc**2/lambdab

      dtrc = [1,0,1,1] * facc2

      do i = 1, 4
        drlxmod(i,:) = dtrc * c(i)
        drlxmod(i,i) = drlxmod(i,i) + facc
      end do

    case(18,19)

!     PTT-XPP single equation
!       model=18: four components in 2D
!       model=19: three components in 2D, czz=1

      lambdab   = vemodel%lambda(mode)
      lambdas   = vemodel%nonlin(1,mode)
      nu        = vemodel%nonlin(2,mode)

      if ( vemodel%model == 19 .and. flowtype == 0 ) then
        trc = c(1) + c(3) + 1._dp
      else
        trc = c(1) + c(3) + c(4)
      end if
      trc3   = 3._dp / trc
      expfac  = 2._dp * exp ( nu * ( sqrt(trc/3._dp) - 1._dp ) ) / lambdas
      facc    = expfac * (1._dp - trc3) + trc3 / lambdab
      facc2   = expfac * ( nu * (1._dp - trc3) / ( 2 * sqrt(3*trc ) ) &
                  + 3/trc**2 ) - 3/trc**2/lambdab

      dtrc = [1,0,1,1] * facc2

      cminI(2) = c(2)
      cminI([1,3]) = c([1,3]) - 1

      if ( vemodel%model == 18 .or. flowtype == 1 ) then

        cminI(4) = c(4) - 1

        do i = 1, 4
          drlxmod(i,:) = dtrc * cminI(i)
          drlxmod(i,i) = drlxmod(i,i) + facc
        end do

      else

        do i = 1, 3
          drlxmod(i,:) = dtrc(1:3) * cminI(i)
          drlxmod(i,i) = drlxmod(i,i) + facc
        end do

      end if

    case(20)

!     FENE-P (see Wapperom & Hulsen 1998, Huetter, Hulsen & Anderson 2018)

      lambda = vemodel%lambda(mode)
      bpar   = vemodel%nonlin(1,mode)

      trc = c(1) + c(3) + c(4)
      f   = bpar / ( bpar + 3 - trc )
      df  = bpar / ( bpar + 3 - trc )**2 * [1,0,1,1]

      drlxmod(1,:) = df * c(1) / lambda
      drlxmod(2,:) = df * c(2) / lambda
      drlxmod(3,:) = df * c(3) / lambda
      drlxmod(4,:) = df * c(4) / lambda
      do i = 1, 4
        drlxmod(i,i) = drlxmod(i,i) + f / lambda
      end do

    case(21)

!     Chilcott-Rallison (alternative formulation)

      lambda = vemodel%lambda(mode)
      rL2    = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
         trc = c(1) + c(3) + 1
      else if ( flowtype == 1 ) then
         trc = c(1) + c(3) + c(4)
      end if
      f = ( rl2 - 3 ) / ( rL2 - trc )
      if ( flowtype == 0 ) then
        df = ( rl2 - 3 ) / ( rL2 - trc )**2 * [1,0,1]
      else if ( flowtype == 1 ) then
        df = ( rl2 - 3 ) / ( rL2 - trc )**2 * [1,0,1,1]
      end if
      fac = f / lambda
      dfac = df / lambda

      drlxmod(1,:) = dfac * ( c(1) - 1 )
      drlxmod(2,:) = dfac * c(2)
      drlxmod(3,:) = dfac * ( c(3) - 1 )
      do i = 1, 3
        drlxmod(i,i) = drlxmod(i,i) + fac
      end do

      if ( flowtype == 1 ) then
        drlxmod(4,:) = dfac * ( c(4) - 1 )
        drlxmod(4,4) = drlxmod(4,4) + fac
      end if

    case(22)

!     Saramito elastoviscoplastic (Drucker-Prager)

      eta = vemodel%modulus(mode) * vemodel%lambda(mode)
      tau_y = vemodel%nonlin(1,mode)
      mu = vemodel%nonlin(2,mode)

      call stress_viscoelastic_2D_single_mode ( vemodel, c, tau, mode )

      trtau = tau(1) + tau(3) + tau(4)  ! trace tau
      tau_e = vonmises_2D ( tau )       ! von Mises equivalent shear stress

!     find the regime the stress state is in

      chpar = mu * trtau - 3 * tau_y
!
      if ( chpar <= - 3 * tau_e ) then

!       regime I: sticking

        drlxmod = 0

      else if ( chpar >= 2 * mu**2 * tau_e ) then

!       regime III: loosing contact

        call dstress_viscoelastic_2D_single_mode ( vemodel, c, dtau, mode )

        drlxmod = dtau / eta

      else

!       regime II: sliding

        call dstress_viscoelastic_2D_single_mode ( vemodel, c, dtau, mode )

        dtrtau = matmul ( [1,0,1,1], dtau )
        dtau_e = matmul ( dvonmises_2D ( tau ), dtau )

        fac = ( tau_e - tau_y + mu * trtau / 3 ) &
                    / ( ( 1 + 2*mu**2/3 ) * tau_e * eta )
        dfac = ( mu/3 * tau_e * dtrtau + &
                       ( tau_y - mu/3 * trtau ) * dtau_e ) &
                    / ( ( 1 + 2*mu**2/3 ) * tau_e**2 * eta )

        f = trtau / 3 - 2 * mu * tau_e / 3
        df = dtrtau / 3 - 2 * mu * dtau_e / 3

        drlxmod = fac * dtau

        drlxmod(1,:) = drlxmod(1,:) - dfac * f - fac * df
        drlxmod(3,:) = drlxmod(3,:) - dfac * f - fac * df
        drlxmod(4,:) = drlxmod(4,:) - dfac * f - fac * df

        do i = 1, 4
          drlxmod(i,:) = drlxmod(i,:) + dfac * tau(i)
        end do

      end if

    case(23,24)

!     EGP

      lambda = vemodel%lambda(mode)

      if ( USE_EGP_STRESS_TENSOR_FORM ) then

        eta = vemodel%modulus(mode) * lambda

        call stress_viscoelastic_2D_single_mode ( vemodel, c, tau, mode, &
          vemmod )
        call dstress_viscoelastic_2D_single_mode ( vemodel, c, dtau, mode, &
          vemmod )

!       delta c.tau

        txx=tau(1); txy=tau(2); tyy=tau(3); tzz=tau(4)

        drlxmod(1,:) = [ txx, txy,    z,    z ]
        drlxmod(2,:) = [ txy, tyy,    z,    z ]
        drlxmod(3,:) = [   z, txy,  tyy,    z ]
        drlxmod(4,:) = [   z,   z,    z,  tzz ]

!       c.delta tau

        cxx=c(1); cxy=c(2); cyy=c(3); czz=c(4)

        w2(1,:) = [ cxx, cxy,    z,    z ]
        w2(2,:) = [   z, cxx,  cxy,    z ]
        w2(3,:) = [   z, cxy,  cyy,    z ]
        w2(4,:) = [   z,   z,    z,  czz ]

        drlxmod = drlxmod + matmul( w2, dtau )

        drlxmod = drlxmod / eta

      else

!       c^2

        cxx=c(1); cxy=c(2); cyy=c(3); czz=c(4)

        drlxmod(1,:) = 2 * [ cxx,     cxy,    z,    z ]
        drlxmod(2,:) =     [ cxy, cxx+cyy,  cxy,    z ]
        drlxmod(3,:) = 2 * [   z,     cxy,  cyy,    z ]
        drlxmod(4,:) = 2 * [   z,       z,    z,  czz ]

!       -trc / 3 * c

        trc = c(1) + c(3) + c(4)
        dtrc = [1,0,1,1]

        do i = 1, 4
          drlxmod(i,:) = drlxmod(i,:) - dtrc / 3 * c(i)
          drlxmod(i,i) = drlxmod(i,i) - trc / 3
        end do

        drlxmod = drlxmod / lambda

        if ( vemodel%model == 24 ) drlxmod = drlxmod / vemmod%J

      end if

!     Jacobian of stabilization term for determinant

      if ( USE_ADAP_LAMBDA_FOR_EGP_DET_STAB ) then
!       Add Jacobian of stabilization term for determinant here
        do i = 1, 4
          drlxmod(i,:) = drlxmod(i,:) + ddetc_2D(c) / 3 * c(i) / lambda
          drlxmod(i,i) = drlxmod(i,i) + ( detc_2D(c) - 1 ) / 3 / lambda
        end do
      end if

    case default

      write(*,'(a,i0)') &
        'Error in drhs_relaxation_2D: model not available: ', vemodel%model
      stop

    end select


    end associate


  end subroutine drhs_relaxation_2D


! derivative of right-hand side wrt L for 2D and axisymm. viscoelastic models.

  subroutine drhsdL_viscoelastic_2D ( vemodel, c, drhsdL, mode1, mode2, &
    mvemodel )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(c,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   conformation tensor c
!   c(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambdab for double equation XPP model
    real(dp), dimension(:,:,:), intent(in) :: c

!   derivative of the right-hand side wrt L of visco-elastic models:
!        .
!        c =  L * c + c * L^T - relaxation.
!
!   drhsdL(ip,i,j,mode) is derivative of component i of rhs with respect
!   to component j of L in point ip for mode number mode
!   column component sequence similar to the conformation tensor c
!   row sequence ordering according to row-major storage of L:
!      2D:           [ Lxx, Lxy, Lyx, Lyy ]
!      axisymmetric: [ Lxx, Lxy, Lyx, Lyy, Ltt ]
!   where Ltt is the hoop strain rate.
    real(dp), dimension(:,:,:,:), intent(out) :: drhsdL

!   if present: limit the range of modes to mode1--mode2
    integer, intent(in), optional :: mode1, mode2


    integer :: mode, flowtype, ncompc, np, ip, modenr1, modenr2
    real(dp) :: cxx, cxy, cyy, czz, xsi
    real(dp) :: fac
    real(dp), parameter :: z=0._dp
    real(dp), dimension(3,4) :: dfL, dfL1

    if ( vemodel%model == 0 ) then

      write(*,'(/3(a/))') &
        'Error in drhsdL_viscoelastic_2D: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 0 .and. vemodel%flowtype /= 1 ) then

      write(*,'(/a/2a/)') &
        'Error in drhsdL_viscoelastic_2D: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=0 (2D) or ', &
        'flowtype=1 (2D, axisymmetric) '
      stop

    end if

    if ( present(mode1) .and. present(mode2) ) then
      if ( mode1 < 1 .or. mode1 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_2D:', &
          ' mode1 out of range '
        stop
      end if
      if ( mode2 < 1 .or. mode2 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_2D:', &
          ' mode2 out of range '
        stop
      end if
      modenr1 = mode1
      modenr2 = mode2
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(c,1) ) then
        write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_2D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_2D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( any( vemodel%model == [14,17] ) ) then
      write(*,'(/a/a,i0/)') 'Error in drhsdL_viscoelastic_2D:', &
        ' derivative matrix not yet available for model = ', vemodel%model
      stop
    endif

    ncompc = vemodel%ncompc
    flowtype = vemodel%flowtype
    np = size(c,1)

    do mode = modenr1, modenr2

      do ip = 1, np

!       extract c

        cxx = c(ip,1,mode); cxy = c(ip,2,mode); cyy = c(ip,3,mode)

        if ( flowtype == 1 ) czz = c(ip,4,mode) ! axisymmetric

!       upper-convected terms

        dfL(1,:) = 2 * [ cxx, cxy, z,   z   ]
        dfL(2,:) =     [ cxy, cyy, cxx, cxy ]
        dfL(3,:) = 2 * [   z,   z, cxy, cyy ]

!       slip term

        if ( vemodel%ixsi > 0  ) then

!         models with slip

          if ( present(mvemodel) ) then
            xsi = mvemodel(ip)%nonlin(vemodel%ixsi,mode)
          else
            xsi = vemodel%nonlin(vemodel%ixsi,mode)
          end if

!         lower-convected terms

          dfL1(1,1:4) = 2 * [ cxx,   z,   cxy,   z ]
          dfL1(2,1:4) =     [ cxy, cxx,   cyy, cxy ]
          dfL1(3,1:4) = 2 * [   z, cxy,   z,   cyy ]

!         Gordon-Schowalter

          dfL = ( 1 - xsi/2 ) * dfL - xsi/2 * dfL1

        end if

        drhsdL(ip,1:3,1:4,mode) = dfL

        drhsdL(ip,4:ncompc,1:4,mode) = 0

        if ( flowtype == 1 ) then

          drhsdL(ip,1:3,5,mode) = 0
          if ( vemodel%ixsi > 0  ) then
            drhsdL(ip,4,5,mode) = 2 * ( 1 - xsi ) * czz
          else
            drhsdL(ip,4,5,mode) = 2 * czz
          end if
          drhsdL(ip,5:ncompc,5,mode) = 0

        end if

!       replace L with L^d = L - (tr L)/3 I

        if ( vemodel%deviatoric ) then

          if ( vemodel%ixsi > 0  ) then
            fac = 2 * ( 1 - xsi ) / 3
          else
            fac = 2._dp / 3
          end if

          drhsdL(ip,1:4,1,mode) = drhsdL(ip,1:4,1,mode) - fac * c(ip,1:4,mode)
          drhsdL(ip,1:4,4,mode) = drhsdL(ip,1:4,4,mode) - fac * c(ip,1:4,mode)
          if ( flowtype == 1 ) then
            drhsdL(ip,1:4,5,mode) = drhsdL(ip,1:4,5,mode) - fac * c(ip,1:4,mode)
          end if

        end if

      end do

    end do

  end subroutine drhsdL_viscoelastic_2D


! von Mises equivalent shear stress for a single mode

  subroutine vm_2D_sm ( vemodel, c, mode, vemmod, vm, dvm, dvmdJ )

    type(vemodel_t), intent(in) :: vemodel
!   c(:ncompc)
    real(dp), dimension(:), intent(in) :: c
    integer, intent(in) :: mode
    type(vemmod_t), intent(in) :: vemmod
    real(dp), intent(out), optional :: vm
!   dvm(:ncompc)
    real(dp), dimension(:), intent(out), optional :: dvm
    real(dp), intent(out), optional :: dvmdJ

    real(dp) :: tau(vemodel%ncompt)
    real(dp) :: dtau(vemodel%ncompt,vemodel%ncompc), dtaudJ(vemodel%ncompt)

    call stress_viscoelastic_2D_single_mode ( vemodel, c, tau, mode, vemmod )

    if ( present(vm) ) vm = vonmises_2D ( tau )

    if ( present(dvm) ) then

      call dstress_viscoelastic_2D_single_mode ( vemodel, c, dtau, mode, &
        vemmod )

      dvm = matmul ( dvonmises_2D(tau), dtau )

    end if

    if ( present(dvmdJ) ) then

      call dstressdJ_viscoelastic_2D_single_mode ( vemodel, c, dtaudJ, mode, &
        vemmod )

      dvmdJ = dot_product ( dvonmises_2D(tau), dtaudJ )

    end if

  end subroutine vm_2D_sm


! von Mises equivalent shear stress based on multiple modes

  subroutine vm_2D_mm ( vemodel, vemmod, vm, dvm, dvmdJ )

    type(vemodel_t), intent(in) :: vemodel
    type(vemmod_t), intent(in) :: vemmod
    real(dp), intent(out), optional :: vm
!   dvm(:ncompc,:nmodes), only mode range vmmode1:vmmode2 is computed
    real(dp), dimension(:,:), intent(out), optional :: dvm
    real(dp), intent(out), optional :: dvmdJ

    integer :: mode, vmmode1, vmmode2
    real(dp), dimension(vemodel%ncompt) :: tau, tsm, dvm2Dtau
    real(dp) :: dtau(vemodel%ncompt,vemodel%ncompc)
    real(dp) :: dtaudJsm(vemodel%ncompt), dtaudJ(vemodel%ncompt)

    vmmode1 = vemmod%vmmode1
    vmmode2 = vemmod%vmmode2

    tau = 0

    do mode = vmmode1, vmmode2

      call stress_viscoelastic_2D_single_mode ( vemodel, vemmod%c(:,mode), &
        tsm, mode, vemmod )

      tau = tau + tsm

    end do

    if ( present(vm) ) vm = vonmises_2D(tau)

    if ( present(dvm) ) then

      dvm2Dtau = dvonmises_2D(tau)

      do mode = vmmode1, vmmode2

        call dstress_viscoelastic_2D_single_mode ( vemodel, vemmod%c(:,mode), &
          dtau, mode, vemmod )

        dvm(:,mode) = matmul ( dvm2Dtau, dtau )

      end do

    end if

    if ( present(dvmdJ) ) then

      dtaudJ = 0

      do mode = vmmode1, vmmode2

        call dstressdJ_viscoelastic_2D_single_mode ( vemodel, &
          vemmod%c(:,mode), dtaudJsm, mode, vemmod )

        dtaudJ = dtaudJ + dtaudJsm

      end do

      dvmdJ = dot_product ( dvonmises_2D(tau), dtaudJ )

    end if

  end subroutine vm_2D_mm


! stress tensor for 2D and axisymmetric viscoelastic models.

  subroutine stress_viscoelastic_2D ( vemodel, c, tau, mode, mode1, mode2, &
    mvemodel, vemopt )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(c,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   conformation tensor c
!   c(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambdab for double equation XPP model

    real(dp), dimension(:,:,:), intent(in) :: c

!   stress tensor tau of the viscoelastic model
!   tau(ip,comp) is component comp in point ip
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop strain
    real(dp), dimension(:,:), intent(out) :: tau

!   if present only a single mode contribution is computed
    integer, intent(in), optional :: mode

!   if present: limit the range of modes to mode1--mode2
    integer, intent(in), optional :: mode1, mode2

!   if present: control the additional optional parameters below
    type(vemopt_t), intent(in), optional :: vemopt


    integer :: modenr1, modenr2, imode, np, ip
    real(dp) :: tsm(vemodel%ncompt)
    type(vemmod_t) :: vemmod


    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(c,1) ) then
        write(*,'(/a/a/)') 'Error in stress_viscoelastic_2D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in stress_viscoelastic_2D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    call check ( vemodel, vemopt, 'stress_viscoelastic_2D', &
      components=['dep_J'] )

    if ( present(vemopt) ) then

!     copy logicals from vemopt to vemmod
      vemmod%dep_J = vemopt%dep_J

    end if

    np = size(c,1)

    if ( present(mode) ) then
      modenr1 = mode
      modenr2 = mode
    else if ( present(mode1) .and. present(mode2) ) then
      modenr1 = mode1
      modenr2 = mode2
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    tau = 0

    do ip = 1, np

!     set J

      if ( vemmod%dep_J ) vemmod%J = vemopt%J(ip)

      do imode = modenr1, modenr2

        if ( present(mvemodel) ) then
          call stress_viscoelastic_2D_single_mode ( mvemodel(ip), &
            c(ip,:,imode), tsm, imode, vemmod )
        else
          call stress_viscoelastic_2D_single_mode ( vemodel, c(ip,:,imode), &
            tsm, imode, vemmod )
        end if

        tau(ip,:) = tau(ip,:) + tsm

      end do

    end do

  end subroutine stress_viscoelastic_2D


! stress tensor for 2D and axisymmetric viscoelastic models (single mode)

  subroutine stress_viscoelastic_2D_single_mode ( vemodel, c, tau, mode, &
    vemmod )

    type(vemodel_t), intent(in) :: vemodel

!   conformation tensor c (single mode)
!   c(comp) is component comp
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambdab for XPP double equation model

    real(dp), dimension(:), intent(in) :: c

!   stress tensor tau of the viscoelastic model (single mode)
!   tau(comp) is component comp
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop stress
    real(dp), dimension(:), intent(out) :: tau

!   the mode number
    integer, intent(in) :: mode

!   additional optional parameters (see type description)
!   NOTE: this argument is set to optional for compatibility reasons.
!   It is checked to be present, when needed.
    type(vemmod_t), intent(in), optional :: vemmod


    integer :: flowtype
    real(dp) :: fac, f, feq, bpar, trc, gamma, xsi, rL2, G
    real(dp) :: cd(4), trc3


    flowtype = vemodel%flowtype

    select case ( vemodel%model )

    case(1,2,3,5,6,7,13,19)

!     Leonov, Maxwell/Oldroyd, Giesekus, Phan-Thien/Tanner,
!     Extended Leonov, Marucci, PTT-XPP (19)

      G = vemodel%modulus(mode)

      tau(1) = G * ( c(1) - 1 )
      tau(2) = G * c(2)
      tau(3) = G * ( c(3) - 1 )
      if ( flowtype == 1 ) then
        tau(4) = G * ( c(4) - 1 )
      end if

    case(4)

!     Larson

      G = vemodel%modulus(mode)
      gamma  = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
        fac = 1 + gamma * ( c(1) + c(3) - 2 )
      else if ( flowtype == 1 ) then
        fac = 1 + gamma * ( c(1) + c(3) + c(4) - 3 )
      end if

      tau(1) = G * ( c(1) / fac  - 1 )
      tau(2) = G * c(2) / fac
      tau(3) = G * ( c(3) / fac  - 1 )
      if ( flowtype == 1 ) then
        tau(4) = G * ( c(4) / fac  - 1 )
      end if

    case(8)

!     Chilcott-Rallison

      G   = vemodel%modulus(mode)
      rL2 = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
         trc = c(1) + c(3) + 1
      else if ( flowtype == 1 ) then
         trc = c(1) + c(3) + c(4)
      end if
      f   = 1 / ( 1 - trc / rL2 )
      feq = rL2 / ( rL2 - 3 )

      tau(1) = G * ( f * c(1) - feq )
      tau(2) = G * f * c(2)
      tau(3) = G * ( f * c(3) - feq )
      if ( flowtype == 1 ) then
         tau(4) = G * ( f * c(4) - feq )
      end if

    case(9)

!     FENE-P

      G    = vemodel%modulus(mode)
      bpar = vemodel%nonlin(1,mode)

      trc = c(1) + c(3) + c(4)
      if ( trc >= bpar + 3 ) then
        write(*,'(/a/a,es16.6/)') &
          'Error in stress_viscoelastic_2D_single_mode:', &
          ' trc > b+3 for FENE-P model. trc=', trc
        stop
      end if
      f   = 1 / ( 1 - trc / ( bpar + 3 ) )
      fac = bpar / ( bpar + 3 ) * f

      tau(1) = G * ( fac * c(1) - 1 )
      tau(2) = G * fac * c(2)
      tau(3) = G * ( fac * c(3) - 1 )
      tau(4) = G * ( fac * c(4) - 1 )

    case(10,11,12)

!     Johnson-Segalman, Phan-Thien/Tanner with slip

      G   = vemodel%modulus(mode)
      xsi = vemodel%nonlin(vemodel%ixsi,mode)

      fac  = G / ( 1 - xsi )

      tau(1) = fac * ( c(1) - 1 )
      tau(2) = fac * c(2)
      tau(3) = fac * ( c(3) - 1 )
      if ( flowtype == 1 ) then
        tau(4) = fac * ( c(4) - 1 )
      end if

    case(14,17)

!     XPP double equation

      G = vemodel%modulus(mode)

!     if backbonestretch is log conformation

      if ( vemodel%model == 14 ) then
        fac = c(5)**2  ! stretch factor
      else if ( vemodel%model == 17 ) then
        fac = exp ( 2._dp * c(5) ) ! stretch factor
      end if

      tau(1) = G * ( 3._dp * fac * c(1) - 1._dp )
      tau(2) = G * ( 3._dp * fac * c(2) )
      tau(3) = G * ( 3._dp * fac * c(3) - 1._dp )
      tau(4) = G * ( 3._dp * fac * c(4) - 1._dp )

    case(15,16,18,22)

!     Rolie-Poly, XPP single eq., PTT-XPP single eq. (18), Saramito DP

      G = vemodel%modulus(mode)

      tau(1) = G * ( c(1) - 1 )
      tau(2) = G * c(2)
      tau(3) = G * ( c(3) - 1 )
      tau(4) = G * ( c(4) - 1 )

    case(20)

!     FENE-P (see Wapperom & Hulsen 1998, Huetter, Hulsen & Anderson 2018)

      G    = vemodel%modulus(mode)
      bpar = vemodel%nonlin(1,mode)

      trc = c(1) + c(3) + c(4)
      if ( trc >= bpar + 3 ) then
        write(*,'(/a/a,es16.6/)') &
          'Error in stress_viscoelastic_2D_single_mode:', &
          ' trc > b+3 for FENE-P model. trc=', trc
        stop
      end if
      f = bpar / ( bpar + 3 - trc )

      tau(1) = G * ( f * c(1) - 1 )
      tau(2) = G * f * c(2)
      tau(3) = G * ( f * c(3) - 1 )
      tau(4) = G * ( f * c(4) - 1 )

    case(21)

!     Chilcott-Rallison (alternative formulation)

      G   = vemodel%modulus(mode)
      rL2 = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
         trc = c(1) + c(3) + 1
      else if ( flowtype == 1 ) then
         trc = c(1) + c(3) + c(4)
      end if
      f = ( rl2 - 3 ) / ( rL2 - trc )

      tau(1) = G * ( f * c(1) - 1 )
      tau(2) = G * f * c(2)
      tau(3) = G * ( f * c(3) - 1 )
      if ( flowtype == 1 ) then
         tau(4) = G * ( f * c(4) - 1 )
      end if

    case(23,24)

!     EGP

      G = vemodel%modulus(mode)

      trc3 = (c(1) + c(3) + c(4))/3

      cd(1) = c(1) - trc3
      cd(2) = c(2)
      cd(3) = c(3) - trc3
      cd(4) = c(4) - trc3

      tau = G * cd

      if ( vemodel%model == 24 ) tau = tau / vemmod%J

    case default

      write(*,'(a,i0)') &
        'Error in stress_relaxation_2D_single_mode: model not available: ', &
        vemodel%model
      stop

    end select

  end subroutine stress_viscoelastic_2D_single_mode


! Jacobian of stress tensor for 2D and axisymmetric viscoelastic models
! (single mode) in multiple points with optional variable coefficients

  subroutine dstress_viscoelastic_2D ( vemodel, c, dtau, mode, mvemodel, &
    vemopt )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(c,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   conformation tensor c
!   c(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambdab for double equation XPP model

    real(dp), dimension(:,:,:), intent(in) :: c

!   derivative of stress tensor tau of the viscoelastic model (single mode)
!   with respect to the conformation tensor
!   dtau(i,j,k) is the derivative in point i of stress component j wrt to
!   component k of c
!     j = 1 tauxx
!     j = 2 tauxy
!     j = 3 tauyy
!     j = 4 tauzz when non-zero or
!     j = 4 tautt axisymmetric hoop stress
    real(dp), dimension(:,:,:), intent(out) :: dtau

!   mode computed
    integer, intent(in) :: mode

!   if present: additional optional parameters (see type description)
    type(vemopt_t), intent(in), optional :: vemopt


    integer :: np, ip
    real(dp) :: dtsm(vemodel%ncompt,vemodel%ncompc)
    type(vemmod_t) :: vemmod

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(c,1) ) then
        write(*,'(/a/a/)') 'Error in dstress_viscoelastic_2D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in dstress_viscoelastic_2D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    call check ( vemodel, vemopt, 'dstress_viscoelastic_2D', &
      components=['dep_J'] )

    if ( present(vemopt) ) then

!     copy logicals from vemopt to vemmod
      vemmod%dep_J = vemopt%dep_J

    end if

    np = size(c,1)

    do ip = 1, np

!     set J

      if ( vemmod%dep_J ) vemmod%J = vemopt%J(ip)

      if ( present(mvemodel) ) then
        call dstress_viscoelastic_2D_single_mode ( mvemodel(ip), &
          c(ip,:,mode), dtsm, mode, vemmod )
      else
        call dstress_viscoelastic_2D_single_mode ( vemodel, c(ip,:,mode), &
          dtsm, mode, vemmod )
      end if

      dtau(ip,:,:) = dtsm

    end do

  end subroutine dstress_viscoelastic_2D


! Jacobian of stress tensor for 2D and axisymmetric viscoelastic models
! (single mode)

  subroutine dstress_viscoelastic_2D_single_mode ( vemodel, c, dtau, mode, &
    vemmod )

    type(vemodel_t), intent(in) :: vemodel

!   conformation tensor c (single mode)
!   c(comp) is component comp
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambdab for XPP double equation model

    real(dp), dimension(:), intent(in) :: c

!   derivative of stress tensor tau of the viscoelastic model (single mode)
!   with respect to the conformation tensor
!   dtau(i,j) is the derivative of stress component i wrt to component j of c
!   tau(comp) is component comp
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop strain
    real(dp), dimension(:,:), intent(out) :: dtau

!   the mode number
    integer, intent(in) :: mode

!   additional optional parameters (see type description)
!   NOTE: this argument is set to optional for compatibility reasons.
!   It is checked to be present, when needed.
    type(vemmod_t), intent(in), optional :: vemmod


    integer :: flowtype, i
    real(dp) :: fac, f, bpar, trc, gamma, xsi, rL2, G
    real(dp), dimension(vemodel%ncompc) :: dfac, df
    real(dp) :: dcd(4,4), dtrc3(4)


    flowtype = vemodel%flowtype

    select case ( vemodel%model )

    case(1,2,3,5,6,7,13,19)

!     Leonov, Maxwell/Oldroyd, Giesekus, Phan-Thien/Tanner,
!     Extended Leonov, Marucci, PTT-XPP (19)

      G = vemodel%modulus(mode)

      dtau = 0
      do i = 1, 3
        dtau(i,i) = G
      end do
      if ( flowtype == 1 ) then
        dtau(4,4) = G
      end if

    case(4)

!     Larson

      G = vemodel%modulus(mode)
      gamma  = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
        fac = 1 + gamma * ( c(1) + c(3) - 2 )
        dfac = gamma * [1,0,1]
      else if ( flowtype == 1 ) then
        fac = 1 + gamma * ( c(1) + c(3) + c(4) - 3 )
        dfac = gamma * [1,0,1,1]
      end if

      dtau(1,:) = - G * dfac * c(1) / fac**2
      dtau(2,:) = - G * dfac * c(2) / fac**2
      dtau(3,:) = - G * dfac * c(3) / fac**2
      do i = 1, 3
        dtau(i,i) = dtau(i,i) + G / fac
      end do
      if ( flowtype == 1 ) then
        dtau(4,:) = - G * dfac * c(4) / fac**2
        dtau(4,4) = dtau(4,4) + G / fac
      end if

    case(8)

!     Chilcott-Rallison

      G   = vemodel%modulus(mode)
      rL2 = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
        trc = c(1) + c(3) + 1
      else if ( flowtype == 1 ) then
        trc = c(1) + c(3) + c(4)
      end if
      f = 1 / ( 1 - trc / rL2 )
      if ( flowtype == 0 ) then
        df = f**2 * [1,0,1] / rL2
      else if ( flowtype == 1 ) then
        df = f**2 * [1,0,1,1] / rL2
      end if

      dtau(1,:) = G * df * c(1)
      dtau(2,:) = G * df * c(2)
      dtau(3,:) = G * df * c(3)
      do i = 1, 3
        dtau(i,i) = dtau(i,i) + G * f
      end do
      if ( flowtype == 1 ) then
        dtau(4,:) = G * df * c(4)
        dtau(4,4) = dtau(4,4) + G * f
      end if

    case(9)

!     FENE-P

      G    = vemodel%modulus(mode)
      bpar = vemodel%nonlin(1,mode)

      trc = c(1) + c(3) + c(4)
      if ( trc >= bpar + 3 ) then
        write(*,'(/a/a,es16.6/)') &
          'Error in stress_viscoelastic_2D_single_mode:', &
          ' trc > b+3 for FENE-P model. trc=', trc
        stop
      end if
      f   = 1 / ( 1 - trc / ( bpar + 3 ) )
      fac = bpar / ( bpar + 3 ) * f
      df  = f**2 * [1,0,1,1] / ( bpar + 3 )
      dfac = bpar / ( bpar + 3 ) * df

      dtau(1,:) = G * dfac * c(1)
      dtau(2,:) = G * dfac * c(2)
      dtau(3,:) = G * dfac * c(3)
      dtau(4,:) = G * dfac * c(4)
      do i = 1, 4
        dtau(i,i) = dtau(i,i) + G * fac
      end do

    case(10,11,12)

!     Johnson-Segalman, Phan-Thien/Tanner with slip

      G   = vemodel%modulus(mode)
      xsi = vemodel%nonlin(vemodel%ixsi,mode)

      fac  = G / ( 1 - xsi )

      dtau = 0
      do i = 1, 3
        dtau(i,i) = fac
      end do
      if ( flowtype == 1 ) then
        dtau(4,4) = fac
      end if

    case(20)

!     FENE-P (see Wapperom & Hulsen 1998, Huetter, Hulsen & Anderson 2018)

      G    = vemodel%modulus(mode)
      bpar = vemodel%nonlin(1,mode)

      trc = c(1) + c(3) + c(4)
      if ( trc >= bpar + 3 ) then
        write(*,'(/a/a,es16.6/)') &
          'Error in stress_viscoelastic_2D_single_mode:', &
          ' trc > b+3 for FENE-P model. trc=', trc
        stop
      end if
      f = bpar / ( bpar + 3 - trc )
      df = bpar / ( bpar + 3 - trc )**2 * [1,0,1,1]

      dtau(1,:) = G * df * c(1)
      dtau(2,:) = G * df * c(2)
      dtau(3,:) = G * df * c(3)
      dtau(4,:) = G * df * c(4)
      do i = 1, 4
        dtau(i,i) = dtau(i,i) + G * f
      end do

    case(21)

!     Chilcott-Rallison (alternative formulation)

      G   = vemodel%modulus(mode)
      rL2 = vemodel%nonlin(1,mode)

      if ( flowtype == 0 ) then
         trc = c(1) + c(3) + 1
      else if ( flowtype == 1 ) then
         trc = c(1) + c(3) + c(4)
      end if
      f = ( rl2 - 3 ) / ( rL2 - trc )
      if ( flowtype == 0 ) then
        df = ( rl2 - 3 ) / ( rL2 - trc )**2 * [1,0,1]
      else if ( flowtype == 1 ) then
        df = ( rl2 - 3 ) / ( rL2 - trc )**2 * [1,0,1,1]
      end if

      dtau(1,:) = G * df * c(1)
      dtau(2,:) = G * df * c(2)
      dtau(3,:) = G * df * c(3)
      do i = 1, 3
        dtau(i,i) = dtau(i,i) + G * f
      end do
      if ( flowtype == 1 ) then
        dtau(4,:) = G * df * c(4)
        dtau(4,4) = dtau(4,4) + G * f
      end if

    case(22)

!     Saramito DP

      G = vemodel%modulus(mode)

      dtau = 0
      do i = 1, 4
        dtau(i,i) = G
      end do

    case(23,24)

!     EGP

      G = vemodel%modulus(mode)

      dtrc3 = [1,0,1,1] / 3._dp

      dcd = 0
      do i = 1, 4
        dcd(i,i) = 1
      end do
      dcd(1,:) = dcd(1,:) - dtrc3
      dcd(3,:) = dcd(3,:) - dtrc3
      dcd(4,:) = dcd(4,:) - dtrc3

      dtau = G * dcd

      if ( vemodel%model == 24 ) dtau = dtau / vemmod%J

    case default

      write(*,'(a,i0)') &
        'Error in dstress_relaxation_2D_single_mode: model not available: ', &
        vemodel%model
      stop

    end select

  end subroutine dstress_viscoelastic_2D_single_mode


! Jacobian with respect to J of stress tensor for 2D viscoelastic models
! (single mode) in multiple points with optional variable coefficients

  subroutine dstressdJ_viscoelastic_2D ( vemodel, c, dtaudJ, mode, mode1, &
    mode2, mvemodel, vemopt )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(c,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   conformation tensor c
!   c(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambdab for double equation XPP model
    real(dp), dimension(:,:,:), intent(in) :: c

!   stress tensor tau of the viscoelastic model
!   tau(ip,comp) is component comp in point ip
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop strain
    real(dp), dimension(:,:), intent(out) :: dtaudJ

!   if present only a single mode contribution is computed
    integer, intent(in), optional :: mode

!   if present: limit the range of modes to mode1--mode2
    integer, intent(in), optional :: mode1, mode2

!   if present: additional optional parameters (see type description)
    type(vemopt_t), intent(in), optional :: vemopt


    integer :: modenr1, modenr2, imode, np, ip
    real(dp) :: dtsm(vemodel%ncompt)
    type(vemmod_t) :: vemmod


    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(c,1) ) then
        write(*,'(/a/a/)') 'Error in dstressdJ_viscoelastic_2D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in dstressdJ_viscoelastic_2D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    call check ( vemodel, vemopt, 'dstressdJ_viscoelastic_2D', &
      components=['dep_J'] )

    if ( present(vemopt) ) then

!     copy logicals from vemopt to vemmod
      vemmod%dep_J = vemopt%dep_J

    end if

    np = size(c,1)

    if ( present(mode) ) then
      modenr1 = mode
      modenr2 = mode
    else if ( present(mode1) .and. present(mode2) ) then
      modenr1 = mode1
      modenr2 = mode2
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    dtaudJ = 0

    do ip = 1, np

!     set J

      if ( vemmod%dep_J ) vemmod%J = vemopt%J(ip)

      do imode = modenr1, modenr2

        if ( present(mvemodel) ) then
          call dstressdJ_viscoelastic_2D_single_mode ( mvemodel(ip), &
            c(ip,:,imode), dtsm, imode, vemmod )
        else
          call dstressdJ_viscoelastic_2D_single_mode ( vemodel, c(ip,:,imode), &
            dtsm, imode, vemmod )
        end if

        dtaudJ(ip,:) = dtaudJ(ip,:) + dtsm

      end do

    end do

  end subroutine dstressdJ_viscoelastic_2D


! Jacobian of stress tensor wrt J for 2D viscoelastic models (single mode)

  subroutine dstressdJ_viscoelastic_2D_single_mode ( vemodel, c, dtaudJ, mode, &
    vemmod )

    type(vemodel_t), intent(in) :: vemodel

!   conformation tensor c (single mode)
!   c(comp) is component comp
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
!     comp = 7 lambdab for XPP double equation model
    real(dp), dimension(:), intent(in) :: c

!   derivative of stress tensor tau of the viscoelastic model (single mode)
!   with respect to J
!   dtau(i)is the derivative of stress component i wrt to J
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauxz
!     comp = 4 tauyy
!     comp = 5 tauyz
!     comp = 6 tauzz
    real(dp), dimension(:), intent(out) :: dtaudJ

!   the mode number
    integer, intent(in) :: mode

!   additional optional parameters (see type description)
    type(vemmod_t), intent(in) :: vemmod



    real(dp) :: tau(vemodel%ncompt)


    select case ( vemodel%model )

    case(24)

!     EGP compressible

      call stress_viscoelastic_2D_single_mode ( vemodel, c, tau, mode, &
        vemmod )
      dtaudJ = - tau / vemmod%J

    case default

      dtaudJ = 0

    end select

  end subroutine dstressdJ_viscoelastic_2D_single_mode


! square root of conformation
! for 2D and axisymmetric viscoelastic models.

  subroutine sqrtc_2D ( c, sqrtc )

    use eig2D3D_m, only: eig2x2, inveig2x2

!   conformation tensor c
!   c(ip,comp) is component comp in point ip
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
    real(dp), dimension(:,:), intent(in) :: c

!   square root of conformation tensor sqrt(c)
!   sqrtc(ip,comp) is component comp in point ip
!   components similar to the conformation tensor c
    real(dp), dimension(:,:), intent(out) :: sqrtc


    integer :: np, ip
    real(dp) :: cpn(size(c,2)), sqrtcpn(size(c,2))
    real(dp) :: sqrtcmdiag(2), eigvalue(2), eigv(2,2)


    np = size(c,1)

    do ip = 1, np

!     extract value of c in one point

      cpn = c(ip,1:size(c,2))

!     compute sqrt(cpn)

      call eig2x2 ( cpn(1:3), eigvalue, eigv )

      sqrtcmdiag = sqrt(eigvalue)

      call inveig2x2 ( sqrtcpn(1:3), sqrtcmdiag, eigv )

      sqrtc(ip,1:3) = sqrtcpn(1:3)
      if ( size(c,2) == 4 ) then
        sqrtc(ip,4:size(c,2)) = sqrt(cpn(4:size(c,2)))
      end if

    end do

  end subroutine sqrtc_2D


! determinant of conformation tensor c

  function detc_2D ( c )

!   conformation tensor c
!   c(comp) is component comp
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
    real(dp), dimension(:), intent(in) :: c
    real(dp) :: detc_2D

    detc_2D = c(1)*c(3) - c(2)**2
    if ( size(c) == 4 ) detc_2D = detc_2D * c(4)

  end function detc_2D


! Jacobian determinant of conformation tensor c

  function ddetc_2D ( c )

!   conformation tensor c
!   c(comp) is component comp
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
    real(dp), dimension(:), intent(in) :: c
    real(dp), dimension(size(c)) :: ddetc_2D

    ddetc_2D(1:3) = [ c(3), -2*c(2), c(1) ]
    if ( size(c) == 4 ) then
      ddetc_2D(1:3) = ddetc_2D(1:3) * c(4)
      ddetc_2D(4) = c(1)*c(3) - c(2)**2
    end if

  end function ddetc_2D


! Right-hand side for plastic strain evolution of 2D models. Optional Jacobian.

  subroutine rhs_plastic_strain_2D ( vemodel, c, gammap, rhs_gammap, &
    mode1, mode2, mvemodel, vemopt_gammap )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(c,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   conformation tensor c
!   c(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
!     comp = 5 lambdab for double equation XPP model
    real(dp), dimension(:,:,:), intent(in) :: c

!   plastic strain gammap
!   gammap(ip) is plastic strain in point ip
    real(dp), dimension(:), intent(in) :: gammap

!   if present: limit the range of modes to mode1--mode2
!   NOTE: this only affects the coupling of modes for von Mises stress
    integer, intent(in), optional :: mode1, mode2

!   right-hand side of the evolution equations of the plastic strain
!   rhs_gammap(ip) in point ip
    real(dp), dimension(:), intent(inout) :: rhs_gammap

!   if present: additional optional parameters (see type description)
    type(vemopt_gammap_t), intent(inout), optional :: vemopt_gammap

    logical :: compute_drhs, compute_drhs_mm
    integer :: mode, ncompc, np, ip, modenr1, modenr2, m, vm1, vm2
    type(vemmod_t) :: vemmod


    if ( vemodel%model == 0 ) then

      write(*,'(/3(a/))') &
        'Error in rhs_plastic_strain_2D: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 0 .and. vemodel%flowtype /= 1 ) then

      write(*,'(/a/2a/)') &
        'Error in rhs_plastic_strain_2D: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=0 (2D) or ', &
        'flowtype=1 (2D, axisymmetric) '
      stop

    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(c,1) ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_2D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_2D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( present(mode1) .and. present(mode2) ) then
      if ( mode1 < 1 .or. mode1 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_2D:', &
          ' mode1 out of range '
        stop
      end if
      if ( mode2 < 1 .or. mode2 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_2D:', &
          ' mode2 out of range '
        stop
      end if
      modenr1 = mode1
      modenr2 = mode2
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    compute_drhs = .false.
    compute_drhs_mm = .false.

    vemmod%dep_gammap = .true.

    if ( present(vemopt_gammap) ) then

!     copy logicals from vemopt_gammap to vemmod
      vemmod%dep_J = vemopt_gammap%dep_J
      vemmod%compute_drhsdJ = vemopt_gammap%compute_drhsdJ
      vemmod%compute_drhsdgammap = vemopt_gammap%compute_drhsdgammap

      if ( vemodel%vmglobal .and. &
                 any([vemopt_gammap%vmmode1,vemopt_gammap%vmmode2]==0) ) then
!       set vmmode to modenr
        vemmod%vmmode1 = modenr1
        vemmod%vmmode2 = modenr2
      else
!       set vmmode to values in vemopt_gammap
        vemmod%vmmode1 = vemopt_gammap%vmmode1
        vemmod%vmmode2 = vemopt_gammap%vmmode2
      end if

      if ( vemopt_gammap%compute_drhs ) then

        if ( vemodel%vmglobal ) then
          write(*,'(/a/a,i0/)') 'Error in rhs_plastic_strain_2D:', &
            ' single mode Jacobian matrix not available for multi-mode', &
            ' von-Mises stress models '
          stop
        endif

        vemmod%compute_drlxmod = .true.

        compute_drhs = .true.

      end if

      if ( vemopt_gammap%compute_drhs_mm ) then

        vemmod%compute_drlxmod = .true.

        allocate( vemmod%drlxmod_mm(1,vemodel%ncompc,vemodel%nmodes) )

        compute_drhs_mm = .true.

      end if

      if ( vemopt_gammap%compute_drhsdgammap ) then
        allocate( vemmod%drlxmoddgammap(1) )
      end if

    else if ( vemodel%vmglobal ) then

!     set default for vmglobal=T: modes for von Mises equal to modenr1:modenr2
      vemmod%vmmode1 = modenr1
      vemmod%vmmode2 = modenr2

    end if

!   always allocate arrays in vemmod for the associate contruct

    allocate(vemmod%drlxmod(1,vemodel%ncompc))
    allocate(vemmod%drlxmoddJ(1))

    if ( vemodel%vmglobal ) then
      allocate ( vemmod%c(size(c,2),size(c,3)) )
      vm1 = vemmod%vmmode1
      vm2 = vemmod%vmmode2
    end if

    mode = vemodel%gammap_mode
    ncompc = vemodel%ncompc
    np = size(c,1)

    do ip = 1, np

!     set gammap

      vemmod%gammap = gammap(ip)

!     set J

      if ( vemmod%dep_J ) vemmod%J = vemopt_gammap%J(ip)

!     von Mises for multiple modes

      if ( vemodel%vmglobal ) vemmod%c(:,vm1:vm2) = c(ip,:,vm1:vm2)

!     fill right-hand side

      if ( present(mvemodel) ) then
        call rhs_plastic_strain_generic ( mvemodel(ip), mode, c(ip,:,mode), &
          vonmises=vm_2D_sm, vonmises_mm=vm_2D_mm, rhs_gammap=rhs_gammap(ip), &
          vemmod=vemmod )
      else
        call rhs_plastic_strain_generic ( vemodel, mode, c(ip,:,mode), &
          vonmises=vm_2D_sm, vonmises_mm=vm_2D_mm, rhs_gammap=rhs_gammap(ip), &
          vemmod=vemmod )
      end if

!     fill derivative (single mode) of right-hand side

      if ( compute_drhs ) then
        vemopt_gammap%drhs(ip,:ncompc) = vemmod%drlxmod(1,:)
      end if

!     fill derivative (multi-mode) with respect to conformation tensor

      if ( compute_drhs_mm ) then

        if ( vemodel%vmglobal ) then

!         mode coupling due to global von Mises stress

!         diagonal
          vemopt_gammap%drhs_mm(ip,:ncompc,mode) = vemmod%drlxmod_mm(1,:,mode)

!         off-diagonal,
          do m = vemmod%vmmode1, vemmod%vmmode2
            if ( m == mode ) cycle ! diagonal already accounted for
            vemopt_gammap%drhs_mm(ip,:ncompc,m) = vemmod%drlxmod_mm(1,:,m)
          end do

        else

!         no mode coupling: only mode given by mode

          vemopt_gammap%drhs_mm(ip,:ncompc,modenr1:mode-1) = 0
          vemopt_gammap%drhs_mm(ip,:ncompc,mode) = vemmod%drlxmod(1,:)
          vemopt_gammap%drhs_mm(ip,:ncompc,mode+1:modenr2) = 0

        end if

      end if

!     fill derivative of right-hand side with respect to J

      if ( vemmod%compute_drhsdJ ) then
        vemopt_gammap%drhsdJ(ip) = vemmod%drlxmoddJ(1)
      end if

!     fill derivative of right-hand side with respect to gammap

      if ( vemmod%compute_drhsdgammap ) then
        vemopt_gammap%drhsdgammap(ip) = vemmod%drlxmoddgammap(1)
      end if

    end do

  end subroutine rhs_plastic_strain_2D

end module viscoelastic_models_2D_m

