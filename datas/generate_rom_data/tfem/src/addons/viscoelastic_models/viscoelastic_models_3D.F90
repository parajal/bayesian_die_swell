
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
! Viscoelastic models in a 3D coordinate system
! Standard case
! Optional Jacobian matrix.
!

module viscoelastic_models_3D_m

  use kind_defs_m
  use viscoelastic_models_defs_m
  use viscoelastic_models_adapted_lambda_m
  use tensor_m
  use viscoelastic_models_plastic_strain_m


  implicit none


contains


! right-hand side for 3D viscoelastic models.

  subroutine rhs_viscoelastic_3D ( vemodel, gradv, c, rhs, mode1, mode2, &
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
!   gradv(ip,i,j) is component L_ij in point ip with
    real(dp), dimension(:,:,:), intent(in) :: gradv

!   conformation tensor c
!   c(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
!     comp = 7 lambda for the XPP double equation
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
        'Error in rhs_viscoelastic_3D: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 2 ) then

      write(*,'(/2(a/))') &
        'Error in rhs_viscoelastic_3D: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=2 '
      stop

    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(c,1) ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_3D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_3D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( present(mode1) .and. present(mode2) ) then
      if ( mode1 < 1 .or. mode1 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_3D:', &
          ' mode1 out of range '
        stop
      end if
      if ( mode2 < 1 .or. mode2 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_3D:', &
          ' mode2 out of range '
        stop
      end if
      modenr1 = mode1
      modenr2 = mode2
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    if ( present(drhs) ) then
      write(*,'(/4(a/))') 'Error in rhs_viscoelastic_3D:', &
        ' Argument drhs is obsolete and will be removed in the future.', &
        ' The output array drhs has been moved to vemopt. Use the routine ', &
        ' create_vemopt with argument compute_drhs=.true. instead. '
      stop
    end if

    call check ( vemodel, vemopt, 'rhs_viscoelastic_3D', &
      components=['dep_J     ', 'dep_gammap', 'drhsdJ    '] )

    compute_drhs = .false.
    compute_drhs_mm = .false.

    if ( present(vemopt) ) then

!     copy logicals from vemopt to vemmod
      vemmod%dep_J = vemopt%dep_J
      vemmod%dep_gammap = vemopt%dep_gammap
      vemmod%compute_drhsdJ = vemopt%compute_drhsdJ
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
          write(*,'(/a/a,i0/)') 'Error in rhs_viscoelastic_3D:', &
            ' single mode Jacobian matrix not available for multi-mode', &
            ' von-Mises stress models '
          stop
        endif
        if ( any( vemodel%model == [14,17] ) ) then
          write(*,'(/a/a,i0/)') 'Error in rhs_viscoelastic_3D:', &
            ' Jacobian matrix not available for model = ', vemodel%model
          stop
        endif
        if ( subtractlinear ) then
          write(*,'(/a/a/)') 'Error in rhs_viscoelastic_3D:', &
            ' Jacobian matrix not available for subtractlinear=.true.'
          stop
        endif

        vemmod%compute_ddefmod = .true.
        vemmod%compute_drlxmod = .true.

        compute_drhs = .true.

      end if

      if ( vemopt%compute_drhs_mm ) then

        if ( any( vemodel%model == [14,17] ) ) then
          write(*,'(/a/a,i0/)') 'Error in rhs_viscoelastic_3D:', &
            ' Jacobian matrix not available for model = ', vemodel%model
          stop
        endif
        if ( subtractlinear ) then
          write(*,'(/a/a/)') 'Error in rhs_viscoelastic_3D:', &
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
          call rhs_deformation_3D ( mvemodel(ip), mode, gradv(ip,:,:), &
            c(ip,:,mode), defmod, vemmod )
        else
          call rhs_deformation_3D ( vemodel, mode, gradv(ip,:,:), &
            c(ip,:,mode), defmod, vemmod )
        end if

!       relaxation term

        if ( present(mvemodel) ) then
          call rhs_relaxation_3D ( mvemodel(ip), mode, c(ip,:,mode), rlxmod, &
            vemmod )
        else
          call rhs_relaxation_3D ( vemodel, mode, c(ip,:,mode), rlxmod, vemmod )
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

  end subroutine rhs_viscoelastic_3D


! right-hand side for 3D viscoelastic models.
! (only relaxation terms)
! Optional Jacobian

  subroutine rhs_viscoelastic_relax_3D ( vemodel, c, rhs, mode1, mode2, &
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
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
!     comp = 7 lambda for the XPP double equation
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
        'Error in rhs_viscoelastic_relax_3D: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 2 ) then

      write(*,'(/2(a/))') &
        'Error in rhs_viscoelastic_relax_3D: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=2 '
      stop

    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(c,1) ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_relax_3D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_relax_3D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( present(drhs) ) then
      write(*,'(/4(a/))') 'Error in rhs_viscoelastic_relax_3D:', &
        ' Argument drhs is obsolete and will be removed in the future.', &
        ' The output array drhs has been moved to vemopt. Use the routine ', &
        ' create_vemopt with argument compute_drhs=.true. instead. '
      stop
    end if

    call check ( vemodel, vemopt, 'rhs_viscoelastic_relax_3D', &
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
      vemmod%dep_gammap = vemopt%dep_gammap
      vemmod%compute_drhsdJ = vemopt%compute_drhsdJ
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
          call rhs_relaxation_3D ( mvemodel(ip), mode, c(ip,:,mode), rlxmod, &
            vemmod )
        else
          call rhs_relaxation_3D ( vemodel, mode, c(ip,:,mode), rlxmod, vemmod )
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

  end subroutine rhs_viscoelastic_relax_3D


! deformation terms for 3D viscoelastic models.

  subroutine rhs_deformation_3D ( vemodel, mode, gradv, c, defmod, vemmod )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   the mode number
    integer, intent(in) :: mode

!   velocity gradient L
!   gradv(i,j) is component L_ij
    real(dp), dimension(:,:), intent(in) :: gradv

!   conformation tensor c
!   c(comp) is component comp
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
!     comp = 7 lambdab for XPP model
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


    real(dp) :: xsi, ddss
    real(dp) :: Lxx, Lxy, Lxz, Lyx, Lyy, Lyz, Lzx, Lzy, Lzz
    real(dp), parameter :: z=0._dp
    real(dp) :: cm(3,3), L(3,3), fm(3,3), trL3


    associate ( ddefmod => vemmod%ddefmod )


!   extract cm

    cm(1,:) = c(1:3)
    cm(2,2:3) = c(4:5)
    cm(3,3) = c(6)
    cm(2,1) = cm(1,2)
    cm(3,1) = cm(1,3)
    cm(3,2) = cm(2,3)

!   extract L

    L = gradv

!   slip term

    if ( vemodel%ixsi > 0  ) then

!     models with slip

      xsi = vemodel%nonlin(vemodel%ixsi,mode)

      L = ( 1 - xsi/2 ) * L - xsi/2 * transpose(L)

    end if

!   deviatoric

    if ( vemodel%deviatoric ) then

!     replace L with L^d = L - (tr L)/3 I

      trL3 = ( L(1,1) + L(2,2) + L(3,3) ) / 3
      L(1,1) = L(1,1) - trL3
      L(2,2) = L(2,2) - trL3
      L(3,3) = L(3,3) - trL3

    end if

!   Quasi-linear terms

    fm = matmul ( L, cm )

    defmod(1) = 2 * fm(1,1)
    defmod(2) = ( fm(1,2) + fm(2,1) )
    defmod(3) = ( fm(1,3) + fm(3,1) )
    defmod(4) = 2 * fm(2,2)
    defmod(5) = ( fm(2,3) + fm(3,2) )
    defmod(6) = 2 * fm(3,3)

    if ( vemmod%compute_ddefmod ) then

!     Jacobian of quasi-linear terms

      Lxx = L(1,1); Lxy=L(1,2); Lxz=L(1,3)
      Lyx = L(2,1); Lyy=L(2,2); Lyz=L(2,3)
      Lzx = L(3,1); Lzy=L(3,2); Lzz=L(3,3)

      ddefmod(1,:) = 2 * [ Lxx,     Lxy,     Lxz,   z,       z,   z ]
      ddefmod(2,:) =     [ Lyx, Lxx+Lyy,     Lyz, Lxy,     Lxz,   z ]
      ddefmod(3,:) =     [ Lzx,     Lzy, Lxx+Lzz,   z,     Lxy, Lxz ]
      ddefmod(4,:) = 2 * [   z,     Lyx,       z, Lyy,     Lyz,   z ]
      ddefmod(5,:) =     [   z,     Lzx,     Lyx, Lzy, Lyy+Lzz, Lyz ]
      ddefmod(6,:) = 2 * [   z,       z,     Lzx,   z,     Lzy, Lzz ]

    end if

!   For XPP double equation add D:S
!   Since all components of L are potentially full this is:

    if ( vemodel%model == 14 .or. vemodel%model == 17 ) then

      ddss = cm(1,1) * L(1,1)  + cm(2,2) * L(2,2) + cm(3,3) * L(3,3) + &
             cm(1,2) * (L(1,2) + L(2,1)) + &
             cm(1,3) * (L(1,3) + L(3,1)) + &
             cm(2,3) * (L(2,3) + L(3,2))

      defmod(1) = defmod(1) - 2._dp * ddss * cm(1,1)
      defmod(2) = defmod(2) - 2._dp * ddss * cm(1,2)
      defmod(3) = defmod(3) - 2._dp * ddss * cm(1,3)
      defmod(4) = defmod(4) - 2._dp * ddss * cm(2,2)
      defmod(5) = defmod(5) - 2._dp * ddss * cm(2,3)
      defmod(6) = defmod(6) - 2._dp * ddss * cm(3,3)

      if ( vemodel%model == 14 ) then
        defmod(7) = ddss * c(7)
      else if ( vemodel%model == 17 ) then
        defmod(7) = ddss
      end if

    end if


    end associate

  end subroutine rhs_deformation_3D


! relaxation term in the models

  subroutine rhs_relaxation_3D ( vemodel, mode, c, rlxmod, vemmod )

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
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
!     comp = 7 lambdab for XPP model
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


    real(dp) :: lambda, alpha, i1, i2, i1mi2, lambd2, fac, beta, eps, yfac
    real(dp) :: gi1pi2, yfacgl, f, feq, bpar, trc, xsi, epsxsi, gamma
    real(dp) :: rL2, expfac, lambdas, nu, relmod, taud, taur, sqtrc, rel1, rel2
    real(dp) :: cm(3,3), cm2(3,3), delta, trc3, lambdab, eta, tau_y, mu
    real(dp) :: trtau, tau_e, chpar, fm2(3,3)
    real(dp) :: tau(vemodel%ncompt)


    if ( mode == vemodel%norlxmode ) then

!     no relaxation term: set relaxation term to zero
      call set_rlx_to_zero

!     add EGP determinant stabilization even if relaxation is set to zero
      call add_egp_det_stabilization

      return

    end if


    select case ( vemodel%model )

    case(1)

!     Leonov

      lambda = vemodel%lambda(mode)

      cm(1,:) = c(1:3)
      cm(2,2:3) = c(4:5)
      cm(3,3) = c(6)
      cm(2,1) = cm(1,2)
      cm(3,1) = cm(1,3)
      cm(3,2) = cm(2,3)

      cm2 = matmul ( cm, cm )

      i1 = cm(1,1) + cm(2,2) + cm(3,3)
      i2 = ( i1 ** 2 - ( cm2(1,1) + cm2(2,2) + cm2(3,3) ) ) / 2
      i1mi2 = ( i1 - i2 ) / 3
      lambd2 = 2 * lambda

      rlxmod(1) = ( cm2(1,1) - i1mi2 * cm(1,1) - 1 ) / lambd2
      rlxmod(2) = ( cm2(1,2) - i1mi2 * cm(1,2) ) / lambd2
      rlxmod(3) = ( cm2(1,3) - i1mi2 * cm(1,3) ) / lambd2
      rlxmod(4) = ( cm2(2,2) - i1mi2 * cm(2,2) - 1 ) / lambd2
      rlxmod(5) = ( cm2(2,3) - i1mi2 * cm(2,3) ) / lambd2
      rlxmod(6) = ( cm2(3,3) - i1mi2 * cm(3,3) - 1 ) / lambd2

    case(2)

!     Maxwell/Oldroyd

      lambda = vemodel%lambda(mode)

      rlxmod(1) = ( c(1) - 1 ) / lambda
      rlxmod(2) = c(2) / lambda
      rlxmod(3) = c(3) / lambda
      rlxmod(4) = ( c(4) - 1 ) / lambda
      rlxmod(5) = c(5) / lambda
      rlxmod(6) = ( c(6) - 1 ) / lambda

    case(3)

!     Giesekus

      lambda = vemodel%lambda(mode)
      alpha  = vemodel%nonlin(1,mode)

      cm(1,:) = c(1:3)
      cm(2,2:3) = c(4:5)
      cm(3,3) = c(6)
      cm(2,1) = cm(1,2)
      cm(3,1) = cm(1,3)
      cm(3,2) = cm(2,3)

      cm(1,1) = cm(1,1) - 1
      cm(2,2) = cm(2,2) - 1
      cm(3,3) = cm(3,3) - 1

      cm2 = matmul ( cm, cm )

      rlxmod(1) = ( cm(1,1) + alpha * cm2(1,1) ) / lambda
      rlxmod(2) = ( cm(1,2) + alpha * cm2(1,2) ) / lambda
      rlxmod(3) = ( cm(1,3) + alpha * cm2(1,3) ) / lambda
      rlxmod(4) = ( cm(2,2) + alpha * cm2(2,2) ) / lambda
      rlxmod(5) = ( cm(2,3) + alpha * cm2(2,3) ) / lambda
      rlxmod(6) = ( cm(3,3) + alpha * cm2(3,3) ) / lambda

    case(4)

!     Larson

      lambda = vemodel%lambda(mode)
      gamma  = vemodel%nonlin(1,mode)

      fac = 1 + gamma * ( c(1) + c(4) + c(6) - 3 )
      fac = fac / lambda

      rlxmod(1) = fac * ( c(1) - 1 )
      rlxmod(2) = fac * c(2)
      rlxmod(3) = fac * c(3)
      rlxmod(4) = fac * ( c(4) - 1 )
      rlxmod(5) = fac * c(5)
      rlxmod(6) = fac * ( c(6) - 1 )

    case(5,6)

!     Phan-Thien/Tanner

      lambda = vemodel%lambda(mode)
      eps    = vemodel%nonlin(1,mode)

      if ( vemodel%model == 5 ) then

!       linear factor

        yfac = 1 + eps * ( c(1) + c(4) + c(6) - 3 )

      else if ( vemodel%model == 6 ) then

!       exponential factor

        yfac = exp ( eps * ( c(1) + c(4) + c(6) - 3 ) )

      end if

      yfacgl = yfac / lambda

      rlxmod(1) = yfacgl * ( c(1) - 1 )
      rlxmod(2) = yfacgl * c(2)
      rlxmod(3) = yfacgl * c(3)
      rlxmod(4) = yfacgl * ( c(4) - 1 )
      rlxmod(5) = yfacgl * c(5)
      rlxmod(6) = yfacgl * ( c(6) - 1 )

    case(7)

!     Extended Leonov

      lambda = vemodel%lambda(mode)
      alpha  = vemodel%nonlin(1,mode)
      beta   = vemodel%nonlin(2,mode)

      cm(1,:) = c(1:3)
      cm(2,2:3) = c(4:5)
      cm(3,3) = c(6)
      cm(2,1) = cm(1,2)
      cm(3,1) = cm(1,3)
      cm(3,2) = cm(2,3)

      cm2 = matmul ( cm, cm )

      i1 = cm(1,1) + cm(2,2) + cm(3,3)
      i2 = ( i1 ** 2 - ( cm2(1,1) + cm2(2,2) + cm2(3,3) ) ) / 2
      i1mi2 = ( i1 - i2 ) / 3
      gi1pi2  = 1 / (  1 + 2 * alpha / pi * &
                        atan( beta / 4 * ( i1 + i2 - 6 ) )  )
      i1mi2  = ( i1 - i2 ) / 3
      lambd2 = 2 * lambda
      fac    = gi1pi2 / lambd2

      rlxmod(1) = fac * ( cm2(1,1) - i1mi2 * cm(1,1) - 1 )
      rlxmod(2) = fac * ( cm2(1,2) - i1mi2 * cm(1,2) )
      rlxmod(3) = fac * ( cm2(1,3) - i1mi2 * cm(1,3) )
      rlxmod(4) = fac * ( cm2(2,2) - i1mi2 * cm(2,2) - 1 )
      rlxmod(5) = fac * ( cm2(2,3) - i1mi2 * cm(2,3) )
      rlxmod(6) = fac * ( cm2(3,3) - i1mi2 * cm(3,3) - 1 )

    case(8)

!     Chilcott-Rallison

      lambda = vemodel%lambda(mode)
      rL2    = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      f = 1 / ( 1 - trc / rL2 )
      fac = f / lambda

      rlxmod(1) = fac * ( c(1) - 1 )
      rlxmod(2) = fac * c(2)
      rlxmod(3) = fac * c(3)
      rlxmod(4) = fac * ( c(4) - 1 )
      rlxmod(5) = fac * c(5)
      rlxmod(6) = fac * ( c(6) - 1 )

    case(9)

!     FENE-P

      lambda = vemodel%lambda(mode)
      bpar   = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      f   = 1 / ( 1 - trc / ( bpar + 3 ) )
      feq = ( bpar + 3 ) / bpar

      rlxmod(1) = ( f * c(1) - feq ) / lambda
      rlxmod(2) = f * c(2) / lambda
      rlxmod(3) = f * c(3) / lambda
      rlxmod(4) = ( f * c(4) - feq )  / lambda
      rlxmod(5) = f * c(5) / lambda
      rlxmod(6) = ( f * c(6) - feq ) / lambda

    case(10)

!     Johnson-Segalman

      lambda = vemodel%lambda(mode)

      rlxmod(1) = ( c(1) - 1 ) / lambda
      rlxmod(2) = c(2) / lambda
      rlxmod(3) = c(3) / lambda
      rlxmod(4) = ( c(4) - 1 ) / lambda
      rlxmod(5) = c(5) / lambda
      rlxmod(6) = ( c(6) - 1 ) / lambda

    case(11,12)

!     Phan-Thien/Tanner with slip

      lambda = vemodel%lambda(mode)
      eps    = vemodel%nonlin(1,mode)
      xsi    = vemodel%nonlin(2,mode)
      epsxsi = eps / ( 1 - xsi )

      if ( vemodel%model == 11 ) then

!       linear factor

        yfac = 1 + epsxsi * ( c(1) + c(4) + c(6) - 3 )

      else if ( vemodel%model == 12 ) then

!       exponential factor

        yfac = exp ( epsxsi * ( c(1) + c(4) + c(6) - 3 ) )

      end if

      yfacgl = yfac / lambda

      rlxmod(1) = yfacgl * ( c(1) - 1 )
      rlxmod(2) = yfacgl * c(2)
      rlxmod(3) = yfacgl * c(3)
      rlxmod(4) = yfacgl * ( c(4) - 1 )
      rlxmod(5) = yfacgl * c(5)
      rlxmod(6) = yfacgl * ( c(6) - 1 )

    case(13)

!     Marrucci

      lambda = vemodel%lambda(mode)
      beta   = vemodel%nonlin(1,mode)

      yfac = 1._dp / ( 1 - beta * ( c(1) + c(4) + c(6) - 3 ) )

      yfacgl = yfac / lambda

      rlxmod(1) = yfacgl * ( c(1) - 1 )
      rlxmod(2) = yfacgl * c(2)
      rlxmod(3) = yfacgl * c(3)
      rlxmod(4) = yfacgl * ( c(4) - 1 )
      rlxmod(5) = yfacgl * c(5)
      rlxmod(6) = yfacgl * ( c(6) - 1 )

    case(14, 17)

!     XPP double equation

      lambda  = vemodel%lambda(mode)
      lambdas = vemodel%nonlin(1,mode)
      nu      = vemodel%nonlin(2,mode)

      if ( vemodel%model == 14 ) then
        lambdab = c(7)
      else if ( vemodel%model == 17 ) then
!       log conformation of backbonestretch
        lambdab = exp( c(7) )
      end if

      expfac  = exp( nu * ( lambdab - 1._dp ) )
      relmod  = 1._dp / ( lambda * lambdab**2 )

      rlxmod(1) = relmod * ( c(1) - 1._dp/3._dp )
      rlxmod(2) = relmod *   c(2)
      rlxmod(3) = relmod *   c(3)
      rlxmod(4) = relmod * ( c(4) - 1._dp/3._dp )
      rlxmod(5) = relmod *   c(5)
      rlxmod(6) = relmod * ( c(6) - 1._dp/3._dp )

      if ( vemodel%model == 14 ) then
        rlxmod(7) = expfac / lambdas * ( lambdab - 1._dp/lambdab )
      else if ( vemodel%model == 17 ) then
        rlxmod(7) = expfac / lambdas * ( 1._dp  - 1._dp/lambdab**2 )
      end if

    case(15)

!     Rolie-Poly

      taud  = vemodel%lambda(mode)
      taur  = vemodel%nonlin(1,mode)
      beta  = vemodel%nonlin(2,mode)
      delta = vemodel%nonlin(3,mode)

      trc   = c(1) + c(4) + c(6)
      sqtrc = sqrt( 3._dp / trc )
      rel1  = ( 1._dp - sqtrc ) / taur
      rel2  = beta * ( trc / 3._dp ) ** delta

      rlxmod(1) = (c(1) - 1._dp)/taud + rel1 * (c(1) + rel2 * ( c(1) - 1._dp ))
      rlxmod(2) = c(2)/taud + rel1 * ( c(2) + rel2 * c(2) )
      rlxmod(3) = c(3)/taud + rel1 * ( c(3) + rel2 * c(3) )
      rlxmod(4) = (c(4) - 1._dp)/taud + rel1 * (c(4) + rel2 * ( c(4) - 1._dp ))
      rlxmod(5) = c(5)/taud + rel1 * ( c(5) + rel2 * c(5) )
      rlxmod(6) = (c(6) - 1._dp)/taud + rel1 * (c(6) + rel2 * ( c(6) - 1._dp ))

    case(16)

!     XPP single equation

      lambdab = vemodel%lambda(mode)
      lambdas = vemodel%nonlin(1,mode)
      nu      = vemodel%nonlin(2,mode)

      trc     = c(1) + c(4) + c(6)
      trc3    = 3._dp / trc
      expfac  = 2._dp * exp ( nu * ( sqrt(trc/3._dp) - 1._dp ) )
      expfac  = expfac * (1._dp - trc3) / lambdas

      rlxmod(1) = expfac * c(1) + ( trc3 * c(1)  - 1._dp ) / lambdab
      rlxmod(2) = expfac * c(2) + ( trc3 * c(2) ) / lambdab
      rlxmod(3) = expfac * c(3) + ( trc3 * c(3) ) / lambdab
      rlxmod(4) = expfac * c(4) + ( trc3 * c(4)  - 1._dp ) / lambdab
      rlxmod(5) = expfac * c(5) + ( trc3 * c(5) ) / lambdab
      rlxmod(6) = expfac * c(6) + ( trc3 * c(6)  - 1._dp ) / lambdab

    case(18,19)

!     PTT-XPP single equation

      lambdab = vemodel%lambda(mode)
      lambdas = vemodel%nonlin(1,mode)
      nu      = vemodel%nonlin(2,mode)

      trc     = c(1) + c(4) + c(6)
      trc3    = 3._dp / trc
      expfac  = 2._dp * exp ( nu * ( sqrt(trc/3._dp) - 1._dp ) )
      expfac  = expfac * (1._dp - trc3) / lambdas

      rlxmod(1) = expfac * c(1) - expfac + ( trc3 * c(1)  - trc3 ) / lambdab
      rlxmod(2) = expfac * c(2) + ( trc3 * c(2) ) / lambdab
      rlxmod(3) = expfac * c(3) + ( trc3 * c(3) ) / lambdab
      rlxmod(4) = expfac * c(4) - expfac + ( trc3 * c(4)  - trc3 ) / lambdab
      rlxmod(5) = expfac * c(5) + ( trc3 * c(5) ) / lambdab
      rlxmod(6) = expfac * c(6) - expfac + ( trc3 * c(6)  - trc3 ) / lambdab

    case(20)

!     FENE-P (see Wapperom & Hulsen 1998, Huetter, Hulsen & Anderson 2018)

      lambda = vemodel%lambda(mode)
      bpar   = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      f   = bpar / ( bpar + 3 - trc )

      rlxmod(1) = ( f * c(1) - 1 ) / lambda
      rlxmod(2) = f * c(2) / lambda
      rlxmod(3) = f * c(3) / lambda
      rlxmod(4) = ( f * c(4) - 1 )  / lambda
      rlxmod(5) = f * c(5) / lambda
      rlxmod(6) = ( f * c(6) - 1 ) / lambda

    case(21)

!     Chilcott-Rallison (alternative formulation)

      lambda = vemodel%lambda(mode)
      rL2    = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      f = ( rl2 - 3 ) / ( rL2 - trc )
      fac = f / lambda

      rlxmod(1) = fac * ( c(1) - 1 )
      rlxmod(2) = fac * c(2)
      rlxmod(3) = fac * c(3)
      rlxmod(4) = fac * ( c(4) - 1 )
      rlxmod(5) = fac * c(5)
      rlxmod(6) = fac * ( c(6) - 1 )

    case(22)

!     Saramito elastoviscoplastic (Drucker-Prager)

      eta = vemodel%modulus(mode) * vemodel%lambda(mode)
      tau_y = vemodel%nonlin(1,mode)
      mu = vemodel%nonlin(2,mode)

      call stress_viscoelastic_3D_single_mode ( vemodel, c, tau, mode )

      trtau = tau(1) + tau(4) + tau(6)  ! trace tau
      tau_e = vonmises_3D ( tau )       ! von Mises equivalent shear stress

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
        rlxmod(3) = fac * tau(3)
        rlxmod(4) = fac * ( tau(4) - f )
        rlxmod(5) = fac * tau(5)
        rlxmod(6) = fac * ( tau(6) - f )

      else

!       regime II: sliding

        fac = ( tau_e - tau_y + mu * trtau / 3 ) &
                    / ( ( 1 + 2*mu**2/3 ) * tau_e * eta )
        f = trtau / 3 - 2 * mu * tau_e / 3

        rlxmod(1) = fac * ( tau(1) - f )
        rlxmod(2) = fac * tau(2)
        rlxmod(3) = fac * tau(3)
        rlxmod(4) = fac * ( tau(4) - f )
        rlxmod(5) = fac * tau(5)
        rlxmod(6) = fac * ( tau(6) - f )

      end if

    case(23,24)

!     EGP

      lambda = vemodel%lambda(mode)

      if ( USE_EGP_STRESS_TENSOR_FORM ) then

        eta = vemodel%modulus(mode) * lambda

        call stress_viscoelastic_3D_single_mode ( vemodel, c, tau, mode, &
          vemmod )

        fm2 = matmul ( vector_to_tensor2_symmetric (ndim=3,A=c), &
                       vector_to_tensor2_symmetric (ndim=3,A=tau) ) / eta

        rlxmod = tensor2_to_vector_symmetric ( ndim=3, A=fm2 )

      else

        cm = vector_to_tensor2_symmetric (ndim=3,A=c)
        cm2 = matmul ( cm, cm )
        trc = c(1) + c(4) + c(6)
        rlxmod = tensor2_to_vector_symmetric (ndim=3,A=cm2-trc/3*cm ) / lambda

        if ( vemodel%model == 24 ) rlxmod = rlxmod / vemmod%J

      end if

      if ( USE_ADAP_LAMBDA_FOR_EGP_DET_STAB ) then
!       Add stabilization for determinant here
        rlxmod = rlxmod + ( detc_3D(c) - 1 ) / 3 * c / lambda
      end if

    case default

      write(*,'(a,i0)') &
        'Error in rhs_relaxation_3D: model not available: ', vemodel%model
      stop

    end select

    if ( vemmod%compute_drlxmod ) then

!     Jacobian of relaxation terms

      call drhs_relaxation_3D ( vemodel, mode, c, vemmod )

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

      call rhs_adapted_lambda ( vemodel, mode, c, vonmises=vm_3D_sm, &
        vonmises_mm=vm_3D_mm, rlxmod=rlxmod, vemmod=vemmod )

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
      rlxmod(3) = rlxmod(3) - c(3) / lambda
      rlxmod(4) = rlxmod(4) - ( c(4) - 1 ) / lambda
      rlxmod(5) = rlxmod(5) - c(5) / lambda
      rlxmod(6) = rlxmod(6) - ( c(6) - 1 ) / lambda

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

      rlxmod = rlxmod + ( detc_3D(c) - 1 ) / 3 * c / lambda

      associate ( drlxmod => vemmod%drlxmod )

      if ( vemmod%compute_drlxmod ) then

!       Jacobian of stabilization term for determinant

        do i = 1, 6
          drlxmod(i,:) = drlxmod(i,:) + ddetc_3D(c) / 3 * c(i) / lambda
          drlxmod(i,i) = drlxmod(i,i) + ( detc_3D(c) - 1 ) / 3 / lambda
        end do

      end if

      end associate

    end subroutine add_egp_det_stabilization

  end subroutine rhs_relaxation_3D


! Jacobian of relaxation term in the models

  subroutine drhs_relaxation_3D ( vemodel, mode, c, vemmod )

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
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
!     comp = 7 lambdab for XPP model
    real(dp), dimension(:), intent(in) :: c

!   additional optional parameters (see type description)
    type(vemmod_t), intent(inout) :: vemmod

    integer :: i
    real(dp) :: lambda, alpha, fac, beta, eps, yfac
    real(dp) :: yfacgl, f, bpar, trc, xsi, epsxsi, gamma
    real(dp) :: rL2, eta, tau_y, mu, cminI(6)
    real(dp) :: trtau, tau_e, chpar, s(6), dyfac(6), dyfacgl(6), dfac(6), df(6)
    real(dp) :: tau(vemodel%ncompt), dtau(vemodel%ncompt,vemodel%ncompc)
    real(dp), dimension(vemodel%ncompc) :: dtrtau, dtau_e
    real(dp) :: gxx, gxy, gxz, gyy, gyz, gzz
    real(dp), parameter :: z=0._dp
    real(dp) :: cxx, cxy, cxz, cyy, cyz, czz, i1, i2, i1mi2, di1mi2(6), lambda2
    real(dp) :: taud, taur, sqtrc, rel1, rel2, rel3, rel4, rel5, lambdas, facc2
    real(dp) :: delta, dtrc(6), facdtrc(6), facc, lambdab, nu, trc3, expfac
    real(dp) :: txx, txy, txz, tyy, tyz, tzz, w2(6,6)


    associate ( drlxmod => vemmod%drlxmod )


    select case ( vemodel%model )

    case(1)

!     Leonov

      lambda = vemodel%lambda(mode)

      cxx = c(1); cxy=c(2); cxz=c(3); cyy=c(4); cyz=c(5); czz=c(6)

!     c^2

      drlxmod(1,:) = 2 * [ cxx,     cxy,     cxz,   z,       z,   z ]
      drlxmod(2,:) =     [ cxy, cxx+cyy,     cyz, cxy,     cxz,   z ]
      drlxmod(3,:) =     [ cxz,     cyz, cxx+czz,   z,     cxy, cxz ]
      drlxmod(4,:) = 2 * [   z,     cxy,       z, cyy,     cyz,   z ]
      drlxmod(5,:) =     [   z,     cxz,     cxy, cyz, cyy+czz, cyz ]
      drlxmod(6,:) = 2 * [   z,       z,     cxz,   z,     cyz, czz ]

!     -(i1-i2)c/3

      i1 = cxx + cyy + czz
      i2 = cxx*cyy + cyy*czz + czz*cxx - cxy**2 - cyz**2 - cxz**2
      i1mi2 = ( i1 - i2 ) / 3
      di1mi2 = [1-cyy-czz,2*cxy,2*cxz,1-cxx-czz,2*cyz,1-cxx-cyy]/3

      do i = 1, 6
        drlxmod(i,:) = drlxmod(i,:) - di1mi2 * c(i)
        drlxmod(i,i) = drlxmod(i,i) - i1mi2
      end do

      lambda2 = 2 * lambda

      drlxmod = drlxmod / lambda2

    case(2)

!     Maxwell/Oldroyd

      lambda = vemodel%lambda(mode)

      drlxmod = 0
      do i = 1, 6
        drlxmod(i,i) = 1 / lambda
      end do

    case(3)

!     Giesekus

      lambda = vemodel%lambda(mode)
      alpha  = vemodel%nonlin(1,mode)

!     s=c-I
      s([2,3,5]) = c([2,3,5])
      s([1,4,6]) = c([1,4,6]) - 1

      drlxmod = 0
      do i = 1, 6
        drlxmod(i,i) = 1
      end do

      gxx=s(1); gxy=s(2); gxz=s(3); gyy=s(4); gyz=s(5); gzz=s(6)

      drlxmod(1,:) = drlxmod(1,:) + &
                       2 * alpha * [ gxx,     gxy,     gxz,   z,       z,   z ]
      drlxmod(2,:) = drlxmod(2,:) + &
                           alpha * [ gxy, gxx+gyy,     gyz, gxy,     gxz,   z ]
      drlxmod(3,:) = drlxmod(3,:) + &
                           alpha * [ gxz,     gyz, gxx+gzz,   z,     gxy, gxz ]
      drlxmod(4,:) = drlxmod(4,:) + &
                       2 * alpha * [   z,     gxy,       z, gyy,     gyz,   z ]
      drlxmod(5,:) = drlxmod(5,:) + &
                           alpha * [   z,     gxz,     gxy, gyz, gyy+gzz, gyz ]
      drlxmod(6,:) = drlxmod(6,:) + &
                       2 * alpha * [   z,       z,     gxz,   z,     gyz, gzz ]

      drlxmod = drlxmod / lambda

    case(4)

!     Larson

      lambda = vemodel%lambda(mode)
      gamma  = vemodel%nonlin(1,mode)

      fac = 1 + gamma * ( c(1) + c(4) + c(6) - 3 )
      fac = fac / lambda
      dfac = gamma * [1,0,0,1,0,1] / lambda

      drlxmod(1,:) = dfac * ( c(1) - 1 )
      drlxmod(2,:) = dfac * c(2)
      drlxmod(3,:) = dfac * c(3)
      drlxmod(4,:) = dfac * ( c(4) - 1 )
      drlxmod(5,:) = dfac * c(5)
      drlxmod(6,:) = dfac * ( c(6) - 1 )
      do i = 1, 6
        drlxmod(i,i) = drlxmod(i,i) + fac
      end do

    case(5,6)

!     Phan-Thien/Tanner

      lambda = vemodel%lambda(mode)
      eps    = vemodel%nonlin(1,mode)

      if ( vemodel%model == 5 ) then

!       linear factor

        yfac = 1 + eps * ( c(1) + c(4) + c(6) - 3 )
        dyfac = eps * [1,0,0,1,0,1]

      else if ( vemodel%model == 6 ) then

!       exponential factor

        yfac = exp ( eps * ( c(1) + c(4) + c(6) - 3 ) )
        dyfac = eps * yfac * [1,0,0,1,0,1]

      end if

      yfacgl = yfac / lambda
      dyfacgl = dyfac / lambda

      drlxmod(1,:) = dyfacgl * ( c(1) - 1 )
      drlxmod(2,:) = dyfacgl * c(2)
      drlxmod(3,:) = dyfacgl * c(3)
      drlxmod(4,:) = dyfacgl * ( c(4) - 1 )
      drlxmod(5,:) = dyfacgl * c(5)
      drlxmod(6,:) = dyfacgl * ( c(6) - 1 )
      do i = 1, 6
        drlxmod(i,i) = drlxmod(i,i) + yfacgl
      end do

    case(8)

!     Chilcott-Rallison

      lambda = vemodel%lambda(mode)
      rL2    = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      f = 1 / ( 1 - trc / rL2 )
      fac = f / lambda
      dfac = f**2 * [1,0,0,1,0,1] / ( rL2 * lambda )

      drlxmod(1,:) = dfac * ( c(1) - 1 )
      drlxmod(2,:) = dfac * c(2)
      drlxmod(3,:) = dfac * c(3)
      drlxmod(4,:) = dfac * ( c(4) - 1 )
      drlxmod(5,:) = dfac * c(5)
      drlxmod(6,:) = dfac * ( c(6) - 1 )
      do i = 1, 6
        drlxmod(i,i) = drlxmod(i,i) + fac
      end do

    case(9)

!     FENE-P

      lambda = vemodel%lambda(mode)
      bpar   = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      f   = 1 / ( 1 - trc / ( bpar + 3 ) )
      df  = f**2 * [1,0,0,1,0,1] / ( bpar + 3 )

      drlxmod(1,:) = df * c(1) / lambda
      drlxmod(2,:) = df * c(2) / lambda
      drlxmod(3,:) = df * c(3) / lambda
      drlxmod(4,:) = df * c(4) / lambda
      drlxmod(5,:) = df * c(5) / lambda
      drlxmod(6,:) = df * c(6) / lambda
      do i = 1, 6
        drlxmod(i,i) = drlxmod(i,i) + f / lambda
      end do

    case(10)

!     Johnson-Segalman

      lambda = vemodel%lambda(mode)

      drlxmod = 0
      do i = 1, 6
        drlxmod(i,i) = 1 / lambda
      end do

    case(11,12)

!     Phan-Thien/Tanner with slip

      lambda = vemodel%lambda(mode)
      eps    = vemodel%nonlin(1,mode)
      xsi    = vemodel%nonlin(2,mode)
      epsxsi = eps / ( 1 - xsi )

      if ( vemodel%model == 11 ) then

!       linear factor

        yfac = 1 + epsxsi * ( c(1) + c(4) + c(6) - 3 )
        dyfac = epsxsi * [1,0,0,1,0,1]

      else if ( vemodel%model == 12 ) then

!       exponential factor

        yfac = exp ( epsxsi * ( c(1) + c(4) + c(6) - 3 ) )
        dyfac = epsxsi * yfac * [1,0,0,1,0,1]

      end if

      yfacgl = yfac / lambda
      dyfacgl = dyfac / lambda

      drlxmod(1,:) = dyfacgl * ( c(1) - 1 )
      drlxmod(2,:) = dyfacgl * c(2)
      drlxmod(3,:) = dyfacgl * c(3)
      drlxmod(4,:) = dyfacgl * ( c(4) - 1 )
      drlxmod(5,:) = dyfacgl * c(5)
      drlxmod(6,:) = dyfacgl * ( c(6) - 1 )
      do i = 1, 6
        drlxmod(i,i) = drlxmod(i,i) + yfacgl
      end do

    case(13)

!     Marrucci

      lambda = vemodel%lambda(mode)
      beta   = vemodel%nonlin(1,mode)

      yfac = 1._dp / ( 1 - beta * ( c(1) + c(4) + c(6) - 3 ) )
      dyfac = yfac**2 * beta * [1,0,0,1,0,1]

      yfacgl = yfac / lambda
      dyfacgl = dyfac / lambda

      drlxmod(1,:) = dyfacgl * ( c(1) - 1 )
      drlxmod(2,:) = dyfacgl * c(2)
      drlxmod(3,:) = dyfacgl * c(3)
      drlxmod(4,:) = dyfacgl * ( c(4) - 1 )
      drlxmod(5,:) = dyfacgl * c(5)
      drlxmod(6,:) = dyfacgl * ( c(6) - 1 )
      do i = 1, 6
        drlxmod(i,i) = drlxmod(i,i) + yfacgl
      end do

    case(15)

!     Rolie-Poly

      taud  = vemodel%lambda(mode)
      taur  = vemodel%nonlin(1,mode)
      beta  = vemodel%nonlin(2,mode)
      delta = vemodel%nonlin(3,mode)

      trc   = c(1) + c(4) + c(6)
      sqtrc = sqrt( 3._dp / trc )
      rel1  = ( 1._dp - sqtrc ) / taur
      rel2  = beta * ( trc / 3._dp ) ** delta
      rel3  = sqrt( 3._dp ) / ( 2 * taur ) / trc ** 1.5_dp
      rel4  = beta * delta / 3 * ( trc / 3._dp ) ** ( delta - 1 )
      rel5  = rel1 * rel4

      facc = 1._dp/taud + rel1 * ( 1 + rel2 )
      facdtrc(1) = rel3 * ( c(1) + rel2 * (c(1)-1) ) + rel5 * (c(1)-1)
      facdtrc(2) = rel3 * ( c(2) + rel2 * c(2) ) + rel5 * c(2)
      facdtrc(3) = rel3 * ( c(3) + rel2 * c(3) ) + rel5 * c(3)
      facdtrc(4) = rel3 * ( c(4) + rel2 * (c(4)-1) ) + rel5 * (c(4)-1)
      facdtrc(5) = rel3 * ( c(5) + rel2 * c(5) ) + rel5 * c(5)
      facdtrc(6) = rel3 * ( c(6) + rel2 * (c(6)-1) ) + rel5 * (c(6)-1)
      dtrc = [1,0,0,1,0,1]

      drlxmod(1,:) = dtrc * facdtrc(1)
      drlxmod(2,:) = dtrc * facdtrc(2)
      drlxmod(3,:) = dtrc * facdtrc(3)
      drlxmod(4,:) = dtrc * facdtrc(4)
      drlxmod(5,:) = dtrc * facdtrc(5)
      drlxmod(6,:) = dtrc * facdtrc(6)
      do i = 1, 6
        drlxmod(i,i) = drlxmod(i,i) + facc
      end do

    case(16)

!     XPP single equation

      lambdab = vemodel%lambda(mode)
      lambdas = vemodel%nonlin(1,mode)
      nu      = vemodel%nonlin(2,mode)

      trc     = c(1) + c(4) + c(6)
      trc3    = 3._dp / trc
      expfac  = 2._dp * exp ( nu * ( sqrt(trc/3._dp) - 1._dp ) ) / lambdas
      facc    = expfac * (1._dp - trc3) + trc3 / lambdab
      facc2   = expfac * ( nu * (1._dp - trc3) / ( 2 * sqrt(3*trc ) ) &
                  + 3/trc**2 ) - 3/trc**2/lambdab

      dtrc = [1,0,0,1,0,1] * facc2

      do i = 1, 6
        drlxmod(i,:) = dtrc * c(i)
        drlxmod(i,i) = drlxmod(i,i) + facc
      end do

    case(18,19)

!     PTT-XPP single equation

      lambdab = vemodel%lambda(mode)
      lambdas = vemodel%nonlin(1,mode)
      nu      = vemodel%nonlin(2,mode)

      trc     = c(1) + c(4) + c(6)
      trc3    = 3._dp / trc
      expfac  = 2._dp * exp ( nu * ( sqrt(trc/3._dp) - 1._dp ) ) / lambdas
      facc    = expfac * (1._dp - trc3) + trc3 / lambdab
      facc2   = expfac * ( nu * (1._dp - trc3) / ( 2 * sqrt(3*trc ) ) &
                  + 3/trc**2 ) - 3/trc**2/lambdab

      dtrc = [1,0,0,1,0,1] * facc2

      cminI([2,3,5]) = c([2,3,5])
      cminI([1,4,6]) = c([1,4,6]) - 1

      do i = 1, 6
        drlxmod(i,:) = dtrc * cminI(i)
        drlxmod(i,i) = drlxmod(i,i) + facc
      end do

    case(20)

!     FENE-P (see Wapperom & Hulsen 1998, Huetter, Hulsen & Anderson 2018)

      lambda = vemodel%lambda(mode)
      bpar   = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      f   = bpar / ( bpar + 3 - trc )
      df  = bpar / ( bpar + 3 - trc )**2 * [1,0,0,1,0,1]

      drlxmod(1,:) = df * c(1) / lambda
      drlxmod(2,:) = df * c(2) / lambda
      drlxmod(3,:) = df * c(3) / lambda
      drlxmod(4,:) = df * c(4) / lambda
      drlxmod(5,:) = df * c(5) / lambda
      drlxmod(6,:) = df * c(6) / lambda
      do i = 1, 6
        drlxmod(i,i) = drlxmod(i,i) + f / lambda
      end do

    case(21)

!     Chilcott-Rallison (alternative formulation)

      lambda = vemodel%lambda(mode)
      rL2    = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      f = ( rl2 - 3 ) / ( rL2 - trc )
      df = ( rl2 - 3 ) / ( rL2 - trc )**2 * [1,0,0,1,0,1]
      fac = f / lambda
      dfac = df / lambda

      drlxmod(1,:) = dfac * ( c(1) - 1 )
      drlxmod(2,:) = dfac * c(2)
      drlxmod(3,:) = dfac * c(3)
      drlxmod(4,:) = dfac * ( c(4) - 1 )
      drlxmod(5,:) = dfac * c(5)
      drlxmod(6,:) = dfac * ( c(6) - 1 )
      do i = 1, 6
        drlxmod(i,i) = drlxmod(i,i) + fac
      end do

    case(22)

!     Saramito elastoviscoplastic (Drucker-Prager)

      eta = vemodel%modulus(mode) * vemodel%lambda(mode)
      tau_y = vemodel%nonlin(1,mode)
      mu = vemodel%nonlin(2,mode)

      call stress_viscoelastic_3D_single_mode ( vemodel, c, tau, mode )

      trtau = tau(1) + tau(4) + tau(6)  ! trace tau
      tau_e = vonmises_3D ( tau )       ! von Mises equivalent shear stress

!     find the regime the stress state is in

      chpar = mu * trtau - 3 * tau_y

      if ( chpar <= - 3 * tau_e ) then

!       regime I: sticking

        drlxmod = 0

      else if ( chpar >= 2 * mu**2 * tau_e ) then

!       regime III: loosing contact

        call dstress_viscoelastic_3D_single_mode ( vemodel, c, dtau, mode )

        drlxmod = dtau / eta

      else

!       regime II: sliding

        call dstress_viscoelastic_3D_single_mode ( vemodel, c, dtau, mode )

        dtrtau = matmul ( [1,0,0,1,0,1], dtau )
        dtau_e = matmul ( dvonmises_3D ( tau ), dtau )

        fac = ( tau_e - tau_y + mu * trtau / 3 ) &
                    / ( ( 1 + 2*mu**2/3 ) * tau_e * eta )
        dfac = ( mu/3 * tau_e * dtrtau + &
                       ( tau_y - mu/3 * trtau ) * dtau_e ) &
                    / ( ( 1 + 2*mu**2/3 ) * tau_e**2 * eta )

        f = trtau / 3 - 2 * mu * tau_e / 3
        df = dtrtau / 3 - 2 * mu * dtau_e / 3

        drlxmod = fac * dtau

        drlxmod(1,:) = drlxmod(1,:) - dfac * f - fac * df
        drlxmod(4,:) = drlxmod(4,:) - dfac * f - fac * df
        drlxmod(6,:) = drlxmod(6,:) - dfac * f - fac * df

        do i = 1, 6
          drlxmod(i,:) = drlxmod(i,:) + dfac * tau(i)
        end do

      end if

    case(23,24)

!     EGP

      lambda = vemodel%lambda(mode)

      if ( USE_EGP_STRESS_TENSOR_FORM ) then

        eta = vemodel%modulus(mode) * lambda

        call stress_viscoelastic_3D_single_mode ( vemodel, c, tau, mode, &
          vemmod )
        call dstress_viscoelastic_3D_single_mode ( vemodel, c, dtau, mode, &
          vemmod )

!       delta c.tau

        txx=tau(1); txy=tau(2); txz=tau(3); tyy=tau(4); tyz=tau(5); tzz=tau(6)

        drlxmod(1,:) = [ txx, txy, txz,   z,   z,   z ]
        drlxmod(2,:) = [ txy, tyy, tyz,   z,   z,   z ]
        drlxmod(3,:) = [ txz, tyz, tzz,   z,   z,   z ]
        drlxmod(4,:) = [   z, txy,   z, tyy, tyz,   z ]
        drlxmod(5,:) = [   z, txz,   z, tyz, tzz,   z ]
        drlxmod(6,:) = [   z,   z, txz,   z, tyz, tzz ]

!       c.delta tau

        cxx=c(1); cxy=c(2); cxz=c(3); cyy=c(4); cyz=c(5); czz=c(6)

        w2(1,:) = [ cxx, cxy, cxz,   z,   z,   z ]
        w2(2,:) = [   z, cxx,   z, cxy, cxz,   z ]
        w2(3,:) = [   z,   z, cxx,   z, cxy, cxz ]
        w2(4,:) = [   z, cxy,   z, cyy, cyz,   z ]
        w2(5,:) = [   z,   z, cxy,   z, cyy, cyz ]
        w2(6,:) = [   z,   z, cxz,   z, cyz, czz ]

        drlxmod = drlxmod + matmul( w2, dtau )

        drlxmod = drlxmod / eta

      else

!       c^2

        cxx=c(1); cxy=c(2); cxz=c(3); cyy=c(4); cyz=c(5); czz=c(6)

        drlxmod(1,:) = 2 * [ cxx,     cxy,     cxz,   z,       z,   z ]
        drlxmod(2,:) =     [ cxy, cxx+cyy,     cyz, cxy,     cxz,   z ]
        drlxmod(3,:) =     [ cxz,     cyz, cxx+czz,   z,     cxy, cxz ]
        drlxmod(4,:) = 2 * [   z,     cxy,       z, cyy,     cyz,   z ]
        drlxmod(5,:) =     [   z,     cxz,     cxy, cyz, cyy+czz, cyz ]
        drlxmod(6,:) = 2 * [   z,       z,     cxz,   z,     cyz, czz ]

!       -trc / 3 * c

        trc = c(1) + c(4) + c(6)
        dtrc = [1,0,0,1,0,1]

        do i = 1, 6
          drlxmod(i,:) = drlxmod(i,:) - dtrc / 3 * c(i)
          drlxmod(i,i) = drlxmod(i,i) - trc / 3
        end do

        drlxmod = drlxmod / lambda

        if ( vemodel%model == 24 ) drlxmod = drlxmod / vemmod%J

      end if

      if ( USE_ADAP_LAMBDA_FOR_EGP_DET_STAB ) then
!       Add Jacobian of stabilization term for determinant here
        do i = 1, 6
          drlxmod(i,:) = drlxmod(i,:) + ddetc_3D(c) / 3 * c(i) / lambda
          drlxmod(i,i) = drlxmod(i,i) + ( detc_3D(c) - 1 ) / 3 / lambda
        end do
      end if

    case default

      write(*,'(a,i0)') &
        'Error in drhs_relaxation_3D: model not available: ', vemodel%model
      stop

    end select


    end associate


  end subroutine drhs_relaxation_3D


! derivative of right-hand side wrt L for 3D viscoelastic models.

  subroutine drhsdL_viscoelastic_3D ( vemodel, c, drhsdL, mode1, mode2, &
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
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
!     comp = 7 lambda for the XPP double equation
    real(dp), dimension(:,:,:), intent(in) :: c

!   derivative of the right-hand side wrt L of visco-elastic models:
!        .
!        c =  L * c + c * L^T - relaxation.
!
!   drhsdL(ip,i,j,mode) is derivative of component i of rhs with respect
!   to component j of L in point ip for mode number mode
!   column component sequence similar to the conformation tensor c
!   row sequence ordering according to row-major storage of L:
!       [ Lxx, Lxy, Lxz, Lyx, Lyy, Lyz, Lzx, Lzy, Lzz ]
    real(dp), dimension(:,:,:,:), intent(out) :: drhsdL

!   if present: limit the range of modes to mode1--mode2
    integer, intent(in), optional :: mode1, mode2


    integer :: mode, np, ip, modenr1, modenr2, ncompc
    real(dp) :: xsi, fac
    real(dp), dimension(6,9) :: dfL, dfL1
    real(dp) :: cxx, cxy, cxz, cyy, cyz, czz
    real(dp), parameter :: z=0._dp


    if ( vemodel%model == 0 ) then

      write(*,'(/3(a/))') &
        'Error in drhsdL_viscoelastic_3D: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 2 ) then

      write(*,'(/2(a/))') &
        'Error in drhsdL_viscoelastic_3D: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=2 '
      stop

    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(c,1) ) then
        write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_3D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_3D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( present(mode1) .and. present(mode2) ) then
      if ( mode1 < 1 .or. mode1 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_3D:', &
          ' mode1 out of range '
        stop
      end if
      if ( mode2 < 1 .or. mode2 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_3D:', &
          ' mode2 out of range '
        stop
      end if
      modenr1 = mode1
      modenr2 = mode2
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    if ( any( vemodel%model == [14,17] ) ) then
      write(*,'(/a/a,i0/)') 'Error in drhsdL_viscoelastic_3D:', &
        ' derivative matrix not yet available for model = ', vemodel%model
      stop
    endif

    ncompc = vemodel%ncompc
    np = size(c,1)

    do mode = modenr1, modenr2

      do ip = 1, np

!       extract c

        cxx = c(ip,1,mode); cxy = c(ip,2,mode); cxz = c(ip,3,mode)
        cyy = c(ip,4,mode); cyz = c(ip,5,mode); czz = c(ip,6,mode)

!       upper-convected terms

        dfL(1,:) = 2 * [ cxx, cxy, cxz,   z,   z,   z,   z,   z,   z ]
        dfL(2,:) =     [ cxy, cyy, cyz, cxx, cxy, cxz,   z,   z,   z ]
        dfL(3,:) =     [ cxz, cyz, czz,   z,   z,   z, cxx, cxy, cxz ]
        dfL(4,:) = 2 * [   z,   z,   z, cxy, cyy, cyz,   z,   z,   z ]
        dfL(5,:) =     [   z,   z,   z, cxz, cyz, czz, cxy, cyy, cyz ]
        dfL(6,:) = 2 * [   z,   z,   z,   z,   z,   z, cxz, cyz, czz ]

!       slip term

        if ( vemodel%ixsi > 0  ) then

          if ( present(mvemodel) ) then
            xsi = mvemodel(ip)%nonlin(vemodel%ixsi,mode)
          else
            xsi = vemodel%nonlin(vemodel%ixsi,mode)
          end if

!         lower-convected terms

          dfL1(1,:) = 2 * [ cxx,   z,   z, cxy,   z,   z, cxz,   z,   z ]
          dfL1(2,:) =     [ cxy, cxx,   z, cyy, cxy,   z, cyz, cxz,   z ]
          dfL1(3,:) =     [ cxz,   z, cxx, cyz,   z, cxy, czz,   z, cxz ]
          dfL1(4,:) = 2 * [   z, cxy,   z,   z, cyy,   z,   z, cyz,   z ]
          dfL1(5,:) =     [   z, cxz, cxy,   z, cyz, cyy,   z, czz, cyz ]
          dfL1(6,:) = 2 * [   z,   z, cxz,   z,   z, cyz,   z,   z, czz ]

!         Gordon-Schowalter

          dfL = ( 1 - xsi/2 ) * dfL - xsi/2 * dfL1

        end if

!       replace L with L^d = L - (tr L)/3 I

        if ( vemodel%deviatoric ) then

          if ( vemodel%ixsi > 0  ) then
            fac = 2 * ( 1 - xsi ) / 3
          else
            fac = 2._dp / 3
          end if

          dfL(:,1) = dfL(:,1) - fac * c(ip,1:6,mode)
          dfL(:,5) = dfL(:,5) - fac * c(ip,1:6,mode)
          dfL(:,9) = dfL(:,9) - fac * c(ip,1:6,mode)

        end if

        drhsdL(ip,1:6,:,mode) = dfL
        drhsdL(ip,7:ncompc,:,mode) = 0

      end do

    end do

  end subroutine drhsdL_viscoelastic_3D


! von Mises equivalent shear stress for a single mode

  subroutine vm_3D_sm ( vemodel, c, mode, vemmod, vm, dvm, dvmdJ )

    type(vemodel_t), intent(in) :: vemodel
!   c(:ncompc)
    real(dp), dimension(:), intent(in) :: c
    integer, intent(in) :: mode
    type(vemmod_t), intent(in) :: vemmod
    real(dp), intent(out), optional :: vm
!   dvm(:ncompc)
    real(dp), dimension(:), intent(out), optional :: dvm
    real(dp), intent(out), optional :: dvmdJ

    real(dp), dimension(vemodel%ncompt) :: tau
    real(dp) :: dtau(vemodel%ncompt,vemodel%ncompc), dtaudJ(vemodel%ncompt)

    call stress_viscoelastic_3D_single_mode ( vemodel, c, tau, mode, vemmod )

    if ( present(vm) ) vm = vonmises_3D(tau)

    if ( present(dvm) ) then

      call dstress_viscoelastic_3D_single_mode ( vemodel, c, dtau, mode, &
        vemmod )

      dvm = matmul ( dvonmises_3D(tau), dtau )

    end if

    if ( present(dvmdJ) ) then

      call dstressdJ_viscoelastic_3D_single_mode ( vemodel, c, dtaudJ, mode, &
        vemmod )

      dvmdJ = dot_product ( dvonmises_3D(tau), dtaudJ )

    end if

  end subroutine vm_3D_sm


! von Mises equivalent shear stress based on multiple modes

  subroutine vm_3D_mm ( vemodel, vemmod, vm, dvm, dvmdJ )

    type(vemodel_t), intent(in) :: vemodel
    type(vemmod_t), intent(in) :: vemmod
    real(dp), intent(out), optional :: vm
!   dvm(:ncompc,:nmodes), only mode range vmmode1:vmmode2 is computed
    real(dp), dimension(:,:), intent(out), optional :: dvm
    real(dp), intent(out), optional :: dvmdJ

    integer :: mode, vmmode1, vmmode2
    real(dp), dimension(vemodel%ncompt) :: tau, tsm, dvm3Dtau
    real(dp) :: dtau(vemodel%ncompt,vemodel%ncompc)
    real(dp) :: dtaudJsm(vemodel%ncompt), dtaudJ(vemodel%ncompt)

    vmmode1 = vemmod%vmmode1
    vmmode2 = vemmod%vmmode2

    tau = 0

    do mode = vmmode1, vmmode2

      call stress_viscoelastic_3D_single_mode ( vemodel, vemmod%c(:,mode), &
        tsm, mode, vemmod )

      tau = tau + tsm

    end do

    if ( present(vm) ) vm = vonmises_3D(tau)

    if ( present(dvm) ) then

      dvm3Dtau = dvonmises_3D(tau)

      do mode = vmmode1, vmmode2

        call dstress_viscoelastic_3D_single_mode ( vemodel, vemmod%c(:,mode), &
          dtau, mode, vemmod )

        dvm(:,mode) = matmul ( dvm3Dtau, dtau )

      end do

    end if

    if ( present(dvmdJ) ) then

      dtaudJ = 0

      do mode = vmmode1, vmmode2

        call dstressdJ_viscoelastic_3D_single_mode ( vemodel, &
          vemmod%c(:,mode), dtaudJsm, mode, vemmod )

        dtaudJ = dtaudJ + dtaudJsm

      end do

      dvmdJ = dot_product ( dvonmises_3D(tau), dtaudJ )

    end if

  end subroutine vm_3D_mm


! stress tensor for 3D viscoelastic models in multiple points with optional
! variable coefficients

  subroutine stress_viscoelastic_3D ( vemodel, c, tau, mode, mode1, mode2, &
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
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
!     comp = 7 lambdab for XPP double equation model
    real(dp), dimension(:,:,:), intent(in) :: c

!   stress tensor tau of the viscoelastic model
!   tau(ip,comp) is component comp in point ip
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauxz
!     comp = 4 tauyy
!     comp = 5 tauyz
!     comp = 6 tauzz
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
        write(*,'(/a/a/)') 'Error in stress_viscoelastic_3D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in stress_viscoelastic_3D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    call check ( vemodel, vemopt, 'stress_viscoelastic_3D', &
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
          call stress_viscoelastic_3D_single_mode ( mvemodel(ip), &
            c(ip,:,imode), tsm, imode, vemmod )
        else
          call stress_viscoelastic_3D_single_mode ( vemodel, c(ip,:,imode), &
            tsm, imode, vemmod )
        end if

        tau(ip,:) = tau(ip,:) + tsm

      end do

    end do

  end subroutine stress_viscoelastic_3D


! stress tensor for 3D viscoelastic models (single mode)

  subroutine stress_viscoelastic_3D_single_mode ( vemodel, c, tau, mode, &
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

!   stress tensor tau of the viscoelastic model (single mode)
!   tau(comp) is component comp
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauxz
!     comp = 4 tauyy
!     comp = 5 tauyz
!     comp = 6 tauzz
    real(dp), dimension(:), intent(out) :: tau

!   the mode number
    integer, intent(in) :: mode

!   additional optional parameters (see type description)
!   NOTE: this argument is set to optional for compatibility reasons.
!   It is checked to be present, when needed.
    type(vemmod_t), intent(in), optional :: vemmod


    real(dp) :: fac, f, feq, bpar, trc, gamma, xsi, rL2, G
    real(dp) :: cd(6), trc3


    select case ( vemodel%model )

    case(1,2,3,5,6,7,13,15,16,18,19,22)

!     Leonov, Maxwell/Oldroyd, Giesekus, Phan-Thien/Tanner,
!     Extended Leonov, Marucci, Rolie-Poly, XPP single eq.,
!     PTT-XPP single eq., Saramito DP

      G = vemodel%modulus(mode)

      tau(1) = G * ( c(1) - 1 )
      tau(2) = G * c(2)
      tau(3) = G * c(3)
      tau(4) = G * ( c(4) - 1 )
      tau(5) = G * c(5)
      tau(6) = G * ( c(6) - 1 )

    case(4)

!     Larson

      G = vemodel%modulus(mode)
      gamma  = vemodel%nonlin(1,mode)

      fac = 1 + gamma * ( c(1) + c(4) + c(6) - 3 )

      tau(1) = G * ( c(1) / fac - 1 )
      tau(2) = G * c(2) / fac
      tau(3) = G * c(3) / fac
      tau(4) = G * ( c(4) / fac - 1 )
      tau(5) = G * c(5) / fac
      tau(6) = G * ( c(6) / fac - 1 )

    case(8)

!     Chilcott-Rallison

      G   = vemodel%modulus(mode)
      rL2 = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      f   = 1 / ( 1 - trc / rL2 )
      feq = rL2 / ( rL2 - 3 )

      tau(1) = G * ( f * c(1) - feq )
      tau(2) = G * f * c(2)
      tau(3) = G * f * c(3)
      tau(4) = G * ( f * c(4) - feq )
      tau(5) = G * f * c(5)
      tau(6) = G * ( f * c(6) - feq )

    case(9)

!     FENE-P

      G    = vemodel%modulus(mode)
      bpar = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      if ( trc >= bpar + 3 ) then
        write(*,'(/a/a,es16.6/)') &
          'Error in stress_viscoelastic_3D_single_mode:', &
          ' trc > b+3 for FENE-P model. trc=', trc
        stop
      end if
      f   = 1 / ( 1 - trc / ( bpar + 3 ) )
      fac = bpar / ( bpar + 3 ) * f

      tau(1) = G * ( fac * c(1) - 1 )
      tau(2) = G * fac * c(2)
      tau(3) = G * fac * c(3)
      tau(4) = G * ( fac * c(4) - 1 )
      tau(5) = G * fac * c(5)
      tau(6) = G * ( fac * c(6) - 1 )

    case(10,11,12)

!     Johnson-Segalman, Phan-Thien/Tanner with slip

      G   = vemodel%modulus(mode)
      xsi = vemodel%nonlin(vemodel%ixsi,mode)

      fac  = G / ( 1 - xsi )

      tau(1) = fac * ( c(1) - 1 )
      tau(2) = fac * c(2)
      tau(3) = fac * c(3)
      tau(4) = fac * ( c(4) - 1 )
      tau(5) = fac * c(5)
      tau(6) = fac * ( c(6) - 1 )

    case(14,17)

!     XPP

      G  = vemodel%modulus(mode)

      if ( vemodel%model == 14 ) then
        fac = c(7)**2  ! stretch factor
      else if ( vemodel%model == 17 ) then
        fac = exp ( 2._dp * c(7) )
      end if

      tau(1) = G * ( 3._dp * fac * c(1) - 1._dp )
      tau(2) = G * ( 3._dp * fac * c(2) )
      tau(3) = G * ( 3._dp * fac * c(3) )
      tau(4) = G * ( 3._dp * fac * c(4) - 1._dp )
      tau(5) = G * ( 3._dp * fac * c(5) )
      tau(6) = G * ( 3._dp * fac * c(6) - 1._dp )

    case(20)

!     FENE-P (see Wapperom & Hulsen 1998, Huetter, Hulsen & Anderson 2018)

      G    = vemodel%modulus(mode)
      bpar = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      if ( trc >= bpar + 3 ) then
        write(*,'(/a/a,es16.6/)') &
          'Error in stress_viscoelastic_3D_single_mode:', &
          ' trc > b+3 for FENE-P model. trc=', trc
        stop
      end if
      f = bpar / ( bpar + 3 - trc )

      tau(1) = G * ( f * c(1) - 1 )
      tau(2) = G * f * c(2)
      tau(3) = G * f * c(3)
      tau(4) = G * ( f * c(4) - 1 )
      tau(5) = G * f * c(5)
      tau(6) = G * ( f * c(6) - 1 )

    case(21)

!     Chilcott-Rallison (alternative formulation)

      G   = vemodel%modulus(mode)
      rL2 = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      f = ( rl2 - 3 ) / ( rL2 - trc )

      tau(1) = G * ( f * c(1) - 1 )
      tau(2) = G * f * c(2)
      tau(3) = G * f * c(3)
      tau(4) = G * ( f * c(4) - 1 )
      tau(5) = G * f * c(5)
      tau(6) = G * ( f * c(6) - 1 )

    case(23,24)

!     EGP

      G = vemodel%modulus(mode)

      trc3 = (c(1) + c(4) + c(6))/3

      cd(1) = c(1) - trc3
      cd(2) = c(2)
      cd(3) = c(3)
      cd(4) = c(4) - trc3
      cd(5) = c(5)
      cd(6) = c(6) - trc3

      tau = G * cd

      if ( vemodel%model == 24 ) tau = tau / vemmod%J

    case default

      write(*,'(a,i0)') &
        'Error in stress_relaxation_3D_single_mode: model not available: ', &
        vemodel%model
      stop

    end select

  end subroutine stress_viscoelastic_3D_single_mode


! Jacobian of stress tensor for 3D viscoelastic models (single mode) in
! multiple points with optional variable coefficients

  subroutine dstress_viscoelastic_3D ( vemodel, c, dtau, mode, mvemodel, &
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
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
!     comp = 7 lambdab for XPP double equation model
    real(dp), dimension(:,:,:), intent(in) :: c

!   derivative of stress tensor tau of the viscoelastic model (single mode)
!   with respect to the conformation tensor
!   dtau(i,j,k) is the derivative in point i of stress component j wrt to
!   component k of c
!     j = 1 tauxx
!     j = 2 tauxy
!     j = 3 tauxz
!     j = 4 tauyy
!     j = 5 tauyz
!     j = 6 tauzz
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
        write(*,'(/a/a/)') 'Error in dstress_viscoelastic_3D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in dstress_viscoelastic_3D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    call check ( vemodel, vemopt, 'dstress_viscoelastic_3D', &
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
        call dstress_viscoelastic_3D_single_mode ( mvemodel(ip), &
          c(ip,:,mode), dtsm, mode, vemmod )
      else
        call dstress_viscoelastic_3D_single_mode ( vemodel, c(ip,:,mode), &
          dtsm, mode, vemmod )
      end if

      dtau(ip,:,:) = dtsm

    end do

  end subroutine dstress_viscoelastic_3D


! Jacobian of stress tensor for 3D viscoelastic models (single mode)

  subroutine dstress_viscoelastic_3D_single_mode ( vemodel, c, dtau, mode, &
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
!   with respect to the conformation tensor
!   dtau(i,j) is the derivative of stress component i wrt to component j of c
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauxz
!     comp = 4 tauyy
!     comp = 5 tauyz
!     comp = 6 tauzz
    real(dp), dimension(:,:), intent(out) :: dtau

!   the mode number
    integer, intent(in) :: mode

!   additional optional parameters (see type description)
!   NOTE: this argument is set to optional for compatibility reasons.
!   It is checked to be present, when needed.
    type(vemmod_t), intent(in), optional :: vemmod


    integer :: i
    real(dp) :: fac, f, bpar, trc, gamma, xsi, rL2, G, dfac(6), df(6)
    real(dp) :: dcd(6,6), dtrc3(6)


    select case ( vemodel%model )

    case(1,2,3,5,6,7,13,15,16,18,19,22)

!     Leonov, Maxwell/Oldroyd, Giesekus, Phan-Thien/Tanner,
!     Extended Leonov, Marucci, Rolie-Poly, XPP single eq.,
!     PTT-XPP single eq., Saramito DP

      G = vemodel%modulus(mode)

      dtau = 0
      do i = 1, 6
        dtau(i,i) = G
      end do

    case(4)

!     Larson

      G = vemodel%modulus(mode)
      gamma  = vemodel%nonlin(1,mode)

      fac = 1 + gamma * ( c(1) + c(4) + c(6) - 3 )
      dfac = gamma * [1,0,0,1,0,1]

      dtau(1,:) = - G * dfac * c(1) / fac**2
      dtau(2,:) = - G * dfac * c(2) / fac**2
      dtau(3,:) = - G * dfac * c(3) / fac**2
      dtau(4,:) = - G * dfac * c(4) / fac**2
      dtau(5,:) = - G * dfac * c(5) / fac**2
      dtau(6,:) = - G * dfac * c(6) / fac**2
      do i = 1, 6
        dtau(i,i) = dtau(i,i) + G / fac
      end do

    case(8)

!     Chilcott-Rallison

      G   = vemodel%modulus(mode)
      rL2 = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      f   = 1 / ( 1 - trc / rL2 )
      df  = f**2 * [1,0,0,1,0,1] / rL2

      dtau(1,:) = G * df * c(1)
      dtau(2,:) = G * df * c(2)
      dtau(3,:) = G * df * c(3)
      dtau(4,:) = G * df * c(4)
      dtau(5,:) = G * df * c(5)
      dtau(6,:) = G * df * c(6)
      do i = 1, 6
        dtau(i,i) = dtau(i,i) + G * f
      end do

    case(9)

!     FENE-P

      G    = vemodel%modulus(mode)
      bpar = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      if ( trc >= bpar + 3 ) then
        write(*,'(/a/a,es16.6/)') &
          'Error in dstress_viscoelastic_3D_single_mode:', &
          ' trc > b+3 for FENE-P model. trc=', trc
        stop
      end if
      f   = 1 / ( 1 - trc / ( bpar + 3 ) )
      fac = bpar / ( bpar + 3 ) * f
      df  = f**2 * [1,0,0,1,0,1] / ( bpar + 3 )
      dfac = bpar / ( bpar + 3 ) * df

      dtau(1,:) = G * dfac * c(1)
      dtau(2,:) = G * dfac * c(2)
      dtau(3,:) = G * dfac * c(3)
      dtau(4,:) = G * dfac * c(4)
      dtau(5,:) = G * dfac * c(5)
      dtau(6,:) = G * dfac * c(6)
      do i = 1, 6
        dtau(i,i) = dtau(i,i) + G * fac
      end do

    case(10,11,12)

!     Johnson-Segalman, Phan-Thien/Tanner with slip

      G   = vemodel%modulus(mode)
      xsi = vemodel%nonlin(vemodel%ixsi,mode)

      fac  = G / ( 1 - xsi )

      dtau = 0
      do i = 1, 6
        dtau(i,i) = fac
      end do

    case(20)

!     FENE-P (see Wapperom & Hulsen 1998, Huetter, Hulsen & Anderson 2018)

      G    = vemodel%modulus(mode)
      bpar = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      if ( trc >= bpar + 3 ) then
        write(*,'(/a/a,es16.6/)') &
          'Error in stress_viscoelastic_3D_single_mode:', &
          ' trc > b+3 for FENE-P model. trc=', trc
        stop
      end if
      f = bpar / ( bpar + 3 - trc )
      df = bpar / ( bpar + 3 - trc )**2 * [1,0,0,1,0,1]

      dtau(1,:) = G * df * c(1)
      dtau(2,:) = G * df * c(2)
      dtau(3,:) = G * df * c(3)
      dtau(4,:) = G * df * c(4)
      dtau(5,:) = G * df * c(5)
      dtau(6,:) = G * df * c(6)
      do i = 1, 6
        dtau(i,i) = dtau(i,i) + G * f
      end do

    case(21)

!     Chilcott-Rallison (alternative formulation)

      G   = vemodel%modulus(mode)
      rL2 = vemodel%nonlin(1,mode)

      trc = c(1) + c(4) + c(6)
      f = ( rl2 - 3 ) / ( rL2 - trc )
      df = ( rl2 - 3 ) / ( rL2 - trc )**2 * [1,0,0,1,0,1]

      dtau(1,:) = G * df * c(1)
      dtau(2,:) = G * df * c(2)
      dtau(3,:) = G * df * c(3)
      dtau(4,:) = G * df * c(4)
      dtau(5,:) = G * df * c(5)
      dtau(6,:) = G * df * c(6)
      do i = 1, 6
        dtau(i,i) = dtau(i,i) + G * f
      end do

    case(23,24)

!     EGP

      G = vemodel%modulus(mode)

      dtrc3 = [1,0,0,1,0,1] / 3._dp

      dcd = 0
      do i = 1, 6
        dcd(i,i) = 1
      end do
      dcd(1,:) = dcd(1,:) - dtrc3
      dcd(4,:) = dcd(4,:) - dtrc3
      dcd(6,:) = dcd(6,:) - dtrc3

      dtau = G * dcd

      if ( vemodel%model == 24 ) dtau = dtau / vemmod%J

    case default

      write(*,'(a,i0)') &
        'Error in dstress_relaxation_3D_single_mode: model not available: ', &
        vemodel%model
      stop

    end select

  end subroutine dstress_viscoelastic_3D_single_mode


! Jacobian with respect to J of stress tensor for 3D viscoelastic models
! (single mode) in multiple points with optional variable coefficients

  subroutine dstressdJ_viscoelastic_3D ( vemodel, c, dtaudJ, mode, mode1, &
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
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
!     comp = 7 lambdab for XPP double equation model
    real(dp), dimension(:,:,:), intent(in) :: c

!   derivative of stress tensor tau of the viscoelastic model (single mode)
!   with respect to J.
!   dtaudJ(i,j) is the derivative in point i of stress component j wrt to J
!     j = 1 tauxx
!     j = 2 tauxy
!     j = 3 tauxz
!     j = 4 tauyy
!     j = 5 tauyz
!     j = 6 tauzz
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
        write(*,'(/a/a/)') 'Error in dstressdJ_viscoelastic_3D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in dstressdJ_viscoelastic_3D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    call check ( vemodel, vemopt, 'dstressdJ_viscoelastic_3D', &
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
          call dstressdJ_viscoelastic_3D_single_mode ( mvemodel(ip), &
            c(ip,:,imode), dtsm, imode, vemmod )
        else
          call dstressdJ_viscoelastic_3D_single_mode ( vemodel, c(ip,:,imode), &
            dtsm, imode, vemmod )
        end if

        dtaudJ(ip,:) = dtaudJ(ip,:) + dtsm

      end do

    end do

  end subroutine dstressdJ_viscoelastic_3D


! Jacobian of stress tensor wrt J for 3D viscoelastic models (single mode)

  subroutine dstressdJ_viscoelastic_3D_single_mode ( vemodel, c, dtaudJ, mode, &
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


    real(dp) :: tau(6)


    select case ( vemodel%model )

    case(24)

!     EGP compressible

      call stress_viscoelastic_3D_single_mode ( vemodel, c, tau, mode, &
        vemmod )
      dtaudJ = - tau / vemmod%J

    case default

      dtaudJ = 0

    end select

  end subroutine dstressdJ_viscoelastic_3D_single_mode


! square root of conformation for 3D viscoelastic models.

  subroutine sqrtc_3D ( c, sqrtc )

#if NO_LIBHSL
    use eig2D3D_m, only: eig3x3, inveig3x3
#else
    use eig2D3D_m, only: inveig3x3
    use hsl_ea23_m, only: eig3x3
#endif

!   conformation tensor c
!   c(ip,comp) is component comp in point ip
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
    real(dp), dimension(:,:), intent(in) :: c

!   square root of conformation tensor sqrt(c)
!   sqrtc(ip,comp) is component comp in point ip
!   components similar to the conformation tensor c
    real(dp), dimension(:,:), intent(out) :: sqrtc


    integer :: np, ip
    real(dp) :: cpn(size(c,2)), sqrtcpn(size(c,2))
    real(dp) :: sqrtcmdiag(3), eigvalue(3), eigv(3,3)


    np = size(c,1)

    do ip = 1, np

!     extract value of c in one point

      cpn = c(ip,1:size(c,2))

!     compute sqrt(cpn)

      call eig3x3 ( cpn(1:6), eigvalue, eigv )

      sqrtcmdiag = sqrt(eigvalue)

      call inveig3x3 ( sqrtcpn(1:6), sqrtcmdiag, eigv )

      sqrtc(ip,1:6) = sqrtcpn(1:6)

    end do

  end subroutine sqrtc_3D


! determinant of conformation tensor c

  function detc_3D ( c )

!   conformation tensor c
!   c(comp) is component comp
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
    real(dp), dimension(:), intent(in) :: c
    real(dp) :: detc_3D

    detc_3D = c(1)*c(4)*c(6) + 2*c(2)*c(3)*c(5) &
              - c(1)*c(5)**2 - c(4)*c(3)**2 - c(6)*c(2)**2

  end function detc_3D


! Jacobian of determinant of conformation tensor c

  function ddetc_3D ( c )

!   conformation tensor c
!   c(comp) is component comp
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
    real(dp), dimension(:), intent(in) :: c
    real(dp), dimension(size(c)) :: ddetc_3D

    ddetc_3D = [ c(4)*c(6) - c(5)**2, 2*c(3)*c(5) - 2*c(6)*c(2), &
                 2*c(2)*c(5) - 2*c(4)*c(3), c(1)*c(6) - c(3)**2, &
                 2*c(2)*c(3) - 2*c(1)*c(5), c(1)*c(4) - c(2)**2 ]

  end function ddetc_3D


! Right-hand side for plastic strain evolution of 3D models. Optional Jacobian.

  subroutine rhs_plastic_strain_3D ( vemodel, c, gammap, rhs_gammap, &
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
!     comp = 3 cxz
!     comp = 4 cyy
!     comp = 5 cyz
!     comp = 6 czz
!     comp = 7 lambda for the XPP double equation
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
        'Error in rhs_plastic_strain_3D: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 2 ) then

      write(*,'(/2(a/))') &
        'Error in rhs_plastic_strain_3D: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=2 '
      stop

    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(c,1) ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_3D:', &
          ' dimension of mvemodel must be identical to size(c,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_3D:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( present(mode1) .and. present(mode2) ) then
      if ( mode1 < 1 .or. mode1 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_3D:', &
          ' mode1 out of range '
        stop
      end if
      if ( mode2 < 1 .or. mode2 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_3D:', &
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
          write(*,'(/a/a,i0/)') 'Error in rhs_plastic_strain_3D:', &
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
          vonmises=vm_3D_sm, vonmises_mm=vm_3D_mm, rhs_gammap=rhs_gammap(ip), &
          vemmod=vemmod )
      else
        call rhs_plastic_strain_generic ( vemodel, mode, c(ip,:,mode), &
          vonmises=vm_3D_sm, vonmises_mm=vm_3D_mm, rhs_gammap=rhs_gammap(ip), &
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

  end subroutine rhs_plastic_strain_3D

end module viscoelastic_models_3D_m

