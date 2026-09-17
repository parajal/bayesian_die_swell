! Copyright (C) 2018-2025 Martien A. Hulsen
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
! Contravariant deformation tensor formulation
! Optional Jacobian matrix.
!

module viscoelastic_models_2D_b_m

  use glob_defs_m
  use tensor_m
  use set_optional_m
  use viscoelastic_models_defs_m
  use viscoelastic_models_adapted_lambda_m
  use viscoelastic_models_2D_m, only: stress_viscoelastic_2D_single_mode, &
      rhs_relaxation_2D, vonmises_2D, dstress_viscoelastic_2D_single_mode
  use viscoelastic_models_plastic_strain_m

  use eig2D3D_m, only: eig2x2, inveig2x2

  implicit none

contains


! right-hand side for 2D and axisymmetric viscoelastic models.

  subroutine rhs_viscoelastic_2D_b ( vemodel, gradv, b, rhs, mode1, mode2, &
    mvemodel, drhs, vemopt )

    use limits_m, only: USE_CONFORMATION_MODEL

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(b,1) must be valid.
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

!   contravariant deformation tensor b
!   b(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:,:,:), intent(in) :: b

!   right-hand side of constitutive visco-elastic models:
!        .
!        b =  L * b - relaxation.
!
!   rhs(ip,comp,mode) is component comp in point ip for mode number mode
!   components sequence similar to the contravariant deformation tensor b
    real(dp), dimension(:,:,:), intent(inout) :: rhs

!   Obsolete and gives an error message if present.
!   Argument will be removed in a future version.
    real(dp), dimension(:,:,:,:), intent(out), optional :: drhs

!   if present: limit the range of modes to mode1--mode2
    integer, intent(in), optional :: mode1, mode2

!   if present: additional optional parameters (see type description)
    type(vemopt_t), intent(inout), optional :: vemopt


    logical :: compute_drhs, compute_drhs_mm
    integer :: mode, ncompb, np, ip, modenr1, modenr2, bvariant, m, vm1, vm2
    real(dp), dimension(vemodel%ncompb) :: defmod, rlxmod, bvec
    type(vemmod_t) :: vemmod


    if ( USE_CONFORMATION_MODEL ) then

      write(*,'(/2(a/))') &
        'Error in rhs_viscoelastic_2D_b: USE_CONFORMATION_MODEL=.true. ', &
        ' not available. Needs to be fixed.'
      stop

    end if

    if ( vemodel%model == 0 ) then

      write(*,'(/3(a/))') &
        'Error in rhs_viscoelastic_2D_b: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 0 .and. vemodel%flowtype /= 1 ) then

      write(*,'(/a/2a/)') &
        'Error in rhs_viscoelastic_2D_b: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=0 (2D) or ', &
        'flowtype=1 (axisymmetric) '
      stop

    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(b,1) ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D_b:', &
          ' dimension of mvemodel must be identical to size(b,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) .or. &
           any( mvemodel(:)%bvariant /= vemodel%bvariant ) ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D_b:', &
          ' model, flowtype, bvariant in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( present(mode1) .and. present(mode2) ) then
      if ( mode1 < 1 .or. mode1 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D_b:', &
          ' mode1 out of range '
        stop
      end if
      if ( mode2 < 1 .or. mode2 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D_b:', &
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
      write(*,'(/4(a/))') 'Error in rhs_viscoelastic_2D_b:', &
        ' Argument drhs is obsolete and will be removed in the future.', &
        ' The output array drhs has been moved to vemopt. Use the routine ', &
        ' create_vemopt with argument compute_drhs=.true. instead. '
      stop
    end if

    call check ( vemodel, vemopt, 'rhs_viscoelastic_2D_b', &
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
          write(*,'(/a/a,i0/)') 'Error in rhs_viscoelastic_2D_b:', &
            ' single mode Jacobian matrix not available for multi-mode', &
            ' von-Mises stress models '
          stop
        endif
        if ( vemodel%bvariant >= 2 ) then
          write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D_b:', &
            ' Jacobian matrix not available for bvariant >= 2.'
          stop
        endif
        if ( subtractlinear ) then
          write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D_b:', &
            ' Jacobian matrix not available for subtractlinear=.true.'
          stop
        endif

        vemmod%compute_ddefmod = .true.
        vemmod%compute_drlxmod = .true.

        compute_drhs = .true.

      end if

      if ( vemopt%compute_drhs_mm ) then

        if ( vemodel%bvariant >= 2 ) then
          write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D_b:', &
            ' Jacobian matrix not available for bvariant >= 2.'
          stop
        endif
        if ( subtractlinear ) then
          write(*,'(/a/a/)') 'Error in rhs_viscoelastic_2D_b:', &
            ' Jacobian matrix not available for subtractlinear=.true.'
          stop
        endif

        vemmod%compute_ddefmod = .true.
        vemmod%compute_drlxmod = .true.

        compute_drhs_mm = .true.

        allocate(&
          vemmod%drlxmod_mm(vemodel%ncompb,vemodel%ncompb,vemodel%nmodes) )

      end if

      if ( vemopt%compute_drhsdgammap ) then

        allocate(vemmod%drlxmoddgammap(vemodel%ncompb) )

      end if

    else if ( vemodel%vmglobal ) then

!     set default for vmglobal=T: modes for von Mises equal to modenr1:modenr2
      vemmod%vmmode1 = modenr1
      vemmod%vmmode2 = modenr2

    end if

!   always allocate arrays in vemmod for the associate contruct

    allocate(vemmod%ddefmod(vemodel%ncompb,vemodel%ncompb))
    allocate(vemmod%drlxmod(vemodel%ncompb,vemodel%ncompb))
    allocate(vemmod%drlxmoddJ(vemodel%ncompb))

    if ( vemodel%vmglobal ) then
      allocate ( vemmod%b(size(b,2),size(b,3)) )
      vm1 = vemmod%vmmode1
      vm2 = vemmod%vmmode2
    end if

    ncompb = vemodel%ncompb
    bvariant = vemodel%bvariant
    np = size(b,1)

    do ip = 1, np

!     set J

      if ( vemmod%dep_J ) vemmod%J = vemopt%J(ip)

!     set gammap

      if ( vemmod%dep_gammap ) vemmod%gammap = vemopt%gammap(ip)

!     von Mises for multiple modes

      if ( vemodel%vmglobal ) vemmod%b(:,vm1:vm2) = b(ip,:,vm1:vm2)

      do mode = modenr1, modenr2

!       extract b

        bvec = b(ip,:ncompb,mode)

        if ( bvariant == 4 ) then
!         Cholesky log
          bvec(1) = exp(bvec(1))
          bvec(4) = exp(bvec(4))
          if ( ncompb == 5 ) bvec(5) = exp(bvec(5))
        end if

!       deformation term

        if ( present(mvemodel) ) then
          call rhs_deformation_2D_b ( mvemodel(ip), mode, gradv(ip,:), bvec, &
            defmod, vemmod )
        else
          call rhs_deformation_2D_b ( vemodel, mode, gradv(ip,:), bvec, &
            defmod, vemmod )
        end if

!       relaxation term

        if ( USE_CONFORMATION_MODEL ) then

          if ( present(mvemodel) ) then
            call rhs_relaxation_2D_cb ( mvemodel(ip), mode, bvec, rlxmod, &
              vemmod )
          else
            call rhs_relaxation_2D_cb ( vemodel, mode, bvec, rlxmod, vemmod )
          end if

        else

          if ( present(mvemodel) ) then
            call rhs_relaxation_2D_b ( mvemodel(ip), mode, bvec, rlxmod, &
              vemmod )
          else
            call rhs_relaxation_2D_b ( vemodel, mode, bvec, rlxmod, vemmod )
          end if

        end if

!       fill right-hand side

        rhs(ip,:ncompb,mode) = defmod - rlxmod

!       modify rhs for b:A variants

        if ( bvariant >= 2 ) then

          call rhs_bvariants

        end if

!       fill derivative (single mode) of right-hand side

        if ( compute_drhs ) then
          vemopt%drhs(ip,:ncompb,:ncompb,mode) = &
                                          vemmod%ddefmod - vemmod%drlxmod
        end if

!       fill derivative (multi mode) of right-hand side

        if ( compute_drhs_mm ) then

          if ( vemodel%vmglobal ) then

!           mode coupling due to global von Mises stress

!           diagonal
            vemopt%drhs_mm(ip,:ncompb,mode,:ncompb,mode) = &
                                  vemmod%ddefmod - vemmod%drlxmod_mm(:,:,mode)

!           off-diagonal
            do m = vemmod%vmmode1, vemmod%vmmode2
              if ( m == mode ) cycle ! diagonal already accounted for
              vemopt%drhs_mm(ip,:ncompb,mode,:ncompb,m) = &
                                                 - vemmod%drlxmod_mm(:,:,m)
            end do

          else

!           no mode coupling: diagonal drhs

            vemopt%drhs_mm(ip,:ncompb,mode,:ncompb,modenr1:mode-1) = 0
            vemopt%drhs_mm(ip,:ncompb,mode,:ncompb,mode) = &
                                            vemmod%ddefmod - vemmod%drlxmod
            vemopt%drhs_mm(ip,:ncompb,mode,:ncompb,mode+1:modenr2) = 0

          end if

        end if

!       fill derivative of right-hand side with respect to J

        if ( vemmod%compute_drhsdJ ) then
          vemopt%drhsdJ(ip,:ncompb,mode) = - vemmod%drlxmoddJ
        end if

!       fill derivative of right-hand side with respect to gammap

        if ( vemmod%compute_drhsdgammap ) then
          vemopt%drhsdgammap(ip,:ncompb,mode) = - vemmod%drlxmoddgammap
        end if

     end do

    end do

  contains

!   modify rhs for b:A variants

    subroutine rhs_bvariants

      real(dp) :: a3

      select case(bvariant)

      case(2)

!       symmetric variant

        a3 = ( defmod(2) - defmod(3) ) / ( bvec(1) + bvec(4) )

      case(3,4)

!       Cholesky variant

        a3 = rhs(ip,2,mode) / bvec(1)

      case default

        call errormsg_case_default ( 'rhs_bvariants', &
          'bvariant', int_value=bvariant )

      end select

!     add term b.A

      rhs(ip,1,mode) = rhs(ip,1,mode) + a3 * bvec(2)
      rhs(ip,2,mode) = rhs(ip,2,mode) - a3 * bvec(1)
      rhs(ip,3,mode) = rhs(ip,3,mode) + a3 * bvec(4)
      rhs(ip,4,mode) = rhs(ip,4,mode) - a3 * bvec(3)

      if ( bvariant == 4 ) then

!       Cholesky log

        rhs(ip,1,mode) = rhs(ip,1,mode) / bvec(1)
        rhs(ip,4,mode) = rhs(ip,4,mode) / bvec(4)
        if ( ncompb == 5 ) rhs(ip,5,mode) = rhs(ip,5,mode) / bvec(5)

      end if

    end subroutine rhs_bvariants

  end subroutine rhs_viscoelastic_2D_b


! deformation terms for 2D viscoelastic models.

  subroutine rhs_deformation_2D_b ( vemodel, mode, gradv, b, defmod, vemmod )

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

!   contravariant deformation tensor b
!   b(comp) is component comp
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:), intent(in) :: b

!   relaxation term for a single mode in a single point
!        .
!        b =  deformation(b,L) - relaxation.
!
!   defmod(comp) is component comp
!   components similar to the contravariant deformation tensor b
    real(dp), dimension(:), intent(out) :: defmod

!   additional optional parameters (see type description)
    type(vemmod_t), intent(inout) :: vemmod


    logical :: deviatoric
    integer :: flowtype, ncompb, i, j
    real(dp) :: Lxx, Lxy, Lyx, Lyy, Lzz, xsi
    real(dp) :: bxx, bxy, byx, byy, bzz
    real(dp) :: Lsxx, Lsxy, Lsyx, Lsyy, Lszz
    real(dp) :: trL3
    real(dp), dimension(2,2,2,2) :: w4


    associate ( ddefmod => vemmod%ddefmod )


    ncompb = vemodel%ncompb
    flowtype = vemodel%flowtype
    deviatoric = vemodel%deviatoric

!   extract b

    bxx = b(1)
    bxy = b(2)
    byx = b(3)
    byy = b(4)

    if ( flowtype == 1 .or. deviatoric ) bzz = b(5)

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
        Lszz  =  Lzz - xsi * Lzz  ! axi-symmetric
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

    defmod(1) = bxx * Lsxx + byx * Lsxy
    defmod(2) = bxy * Lsxx + byy * Lsxy
    defmod(3) = bxx * Lsyx + byx * Lsyy
    defmod(4) = bxy * Lsyx + byy * Lsyy

    if ( flowtype == 1 .or. deviatoric ) then
      defmod(5) = bzz * Lszz
    else
      defmod(5:ncompb) = 0._dp
    end if

    if ( vemmod%compute_ddefmod ) then

!     Jacobian of quasi-linear terms

      do i = 1, 2
        w4(:,i,:,i) = reshape ( [Lsxx,Lsxy,Lsyx,Lsyy], [2,2], order=[2,1] )
        do j = i+1, 2
          w4(:,i,:,j) = 0
          w4(:,j,:,i) = 0
        end do
      end do

      ddefmod(1:4,1:4) = tensor4_to_matrix ( ndim=2, A=w4 )

      if ( flowtype == 1 .or. deviatoric ) then
        ddefmod(5,1:4) = 0
        ddefmod(1:4,5) = 0
        ddefmod(5,5) = Lszz
      else
        ddefmod(5:ncompb,:) = 0
        ddefmod(:,5:ncompb) = 0
      end if

    end if


    end associate


  end subroutine rhs_deformation_2D_b


! relaxation term in the models

  subroutine rhs_relaxation_2D_b ( vemodel, mode, b, rlxmod, vemmod )

    use limits_m, only: USE_EGP_STRESS_TENSOR_FORM, &
                        USE_ADAP_LAMBDA_FOR_EGP_DET_STAB

    type(vemodel_t), intent(in) :: vemodel

!   the mode number
    integer, intent(in) :: mode

!   contravariant deformation tensor b
!   b(comp) is component comp
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain

    real(dp), dimension(:), intent(in) :: b

!   relaxation term for a single mode in a single point
!        .
!        b =  L * b - relaxation.
!
!   rlxmod(comp) is component comp
!   components similar to the contravariant deformation tensor b
    real(dp), dimension(:), intent(out) :: rlxmod

!   additional optional parameters (see type description)
    type(vemmod_t), intent(inout) :: vemmod


    integer :: flowtype

    real(dp) :: lambda, alpha, f, bpar, trc
    real(dp) :: bm(2,2), bminv(2,2), cm(2,2), cm2(2,2), rlxmodm(2,2)
    real(dp) :: tau(vemodel%ncompt), fm2(2,2), eta


    if ( mode == vemodel%norlxmode ) then

!     no relaxation term: set relaxation term to zero
      call set_rlx_to_zero

!     add EGP determinant stabilization even if relaxation is set to zero
      call add_egp_det_stabilization

      return

    end if


    flowtype = vemodel%flowtype

    select case ( vemodel%model )

    case(2)

!     Maxwell/Oldroyd

      lambda = vemodel%lambda(mode)
      bm = reshape ( b(1:4), [2,2], order=[2,1] )
      call matinv2 ( bm, bminv )
      rlxmod(1:4) = reshape ( transpose ( bm - transpose ( bminv ) ), [4] )
      rlxmod(1:4) = rlxmod(1:4) / ( 2 * lambda )
      if ( flowtype == 1 ) then
        rlxmod(5) = ( b(5) - 1 / b(5) ) / ( 2 * lambda )
      end if

    case(3)

!     Giesekus

      lambda = vemodel%lambda(mode)
      alpha  = vemodel%nonlin(1,mode)
      bm = reshape ( b(1:4), [2,2], order=[2,1] )
      call matinv2 ( bm, bminv )
      cm = matmul ( bm, transpose ( bm ) )
      cm2 = matmul ( cm, cm )
      rlxmodm = - ( 1 - 2 * alpha ) * cm - alpha * cm2
      rlxmodm(1,1) = rlxmodm(1,1) + ( 1 - alpha )
      rlxmodm(2,2) = rlxmodm(2,2) + ( 1 - alpha )
      rlxmodm = matmul ( rlxmodm, transpose ( bminv ) )
      rlxmodm = - rlxmodm / ( 2 * lambda )
      rlxmod(1:4) = reshape ( transpose ( rlxmodm ), [4] )
      if ( flowtype == 1 ) then
        rlxmod(5) = - ( ( 1 - alpha ) - ( 1 - 2 * alpha ) * b(5)**2 &
                        - alpha * b(5)**4 ) / ( 2 * lambda * b(5) )
      end if

    case(20)

!     FENE-P

      lambda = vemodel%lambda(mode)
      bpar   = vemodel%nonlin(1,mode)

      bm = reshape ( b(1:4), [2,2], order=[2,1] )
      call matinv2 ( bm, bminv )

      cm = matmul ( bm, transpose ( bm ) )
      trc = cm(1,1) + cm(2,2) + b(5)**2
      f   = bpar / ( bpar + 3 - trc )

      rlxmod(1:4) = reshape ( transpose ( f * bm - transpose ( bminv ) ), [4] )
      rlxmod(5) = f * b(5) - 1 / b(5)
      rlxmod(1:5) = rlxmod(1:5) / ( 2 * lambda )

    case(23,24)

!     EGP

      lambda = vemodel%lambda(mode)

      if ( USE_EGP_STRESS_TENSOR_FORM ) then

        eta = vemodel%modulus(mode) * lambda

        call stress_viscoelastic_2D_single_mode_b ( vemodel, b, tau, mode, &
          vemmod )

        fm2 = matmul ( vector_to_tensor2_symmetric (ndim=2,A=tau(1:3)), &
                       vector_to_tensor2 (ndim=2,A=b(1:4)) ) / (2*eta)

        rlxmod(1:4) = tensor2_to_vector ( ndim=2, A=fm2 )
        rlxmod(5) = tau(4) * b(5) / (2*eta)

      else

        bm = reshape ( b(1:4), [2,2], order=[2,1] )
        cm = matmul ( bm, transpose ( bm ) )
        trc = cm(1,1) + cm(2,2) + b(5)**2
        rlxmod(1:4) = &
               tensor2_to_vector (ndim=2,A=matmul(cm,bm)-trc/3*bm)/(2*lambda)
        rlxmod(5) = ( b(5)**2 * b(5) - trc/3*b(5) ) / (2*lambda)

        if ( vemodel%model == 24 ) rlxmod = rlxmod / vemmod%J

      end if

      if ( USE_ADAP_LAMBDA_FOR_EGP_DET_STAB ) then
!       Add stabilization for determinant here
        rlxmod = rlxmod + ( detb_2D(b) - 1 ) * b / (3*lambda)
!        rlxmod = rlxmod + ( detb_2D(b)**2 - 1 ) * b / (6*lambda)
      end if

    case default

      write(*,'(a,i0,3(/a)/)') &
        'Error in rhs_relaxation_2D_b: model not available: ', vemodel%model, &
        '  Set USE_CONFORMATION_MODEL = .true. in the limits_m model ', &
        '  to use a transformation of the conformation model formulation ', &
        '  to the contravariant deformation formulation.'
      stop

    end select

    if ( vemmod%compute_drlxmod ) then

!     Jacobian of relaxation terms

      call drhs_relaxation_2D_b ( vemodel, mode, b, vemmod )

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

      call rhs_adapted_lambda ( vemodel, mode, b, vonmises=vm_2D_sm_b, &
        vonmises_mm=vm_2D_mm_b, rlxmod=rlxmod, vemmod=vemmod )

    end if

    if ( .not. USE_ADAP_LAMBDA_FOR_EGP_DET_STAB ) then
!     Add stabilization for determinant here with original lambda constant
      call add_egp_det_stabilization
    end if


    if ( subtractlinear ) then

!     Subtract the linear part.

      write(*,'(/2(a/))') &
        'Error in rhs_relaxation_2D_b: subtractlinear', &
        'not yet implemented for the b-formulation'
      stop

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

      rlxmod = rlxmod + ( detb_2D(b) - 1 ) * b / (3*lambda)
!      rlxmod = rlxmod + ( detb_2D(b)**2 - 1 ) * b / (6*lambda)

      associate ( drlxmod => vemmod%drlxmod )

      if ( vemmod%compute_drlxmod ) then

!       Jacobian of stabilization term for determinant

        do i = 1, 5
          drlxmod(i,:) = drlxmod(i,:) + ddetb_2D(b) * b(i) / ( 3 * lambda )
          drlxmod(i,i) = drlxmod(i,i) + ( detb_2D(b)- 1 ) / ( 3 * lambda )
!          drlxmod(i,:) = drlxmod(i,:) &
!                         + 2 * detb_2D(b) * ddetb_2D(b) * b(i) / ( 6 * lambda )
!          drlxmod(i,i) = drlxmod(i,i) + ( detb_2D(b)**2 - 1 ) / ( 6 * lambda )
        end do

      end if

      end associate

    end subroutine add_egp_det_stabilization

  end subroutine rhs_relaxation_2D_b


! relaxation term in the models derived from the conformation models

  subroutine rhs_relaxation_2D_cb ( vemodel, mode, b, rlxmod, vemmod )

    type(vemodel_t), intent(in) :: vemodel

!   the mode number
    integer, intent(in) :: mode

!   contravariant deformation tensor b
!   b(comp) is component comp
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain

    real(dp), dimension(:), intent(in) :: b

!   relaxation term for a single mode in a single point
!        .
!        b =  L * b - relaxation.
!
!   rlxmod(comp) is component comp
!   components similar to the contravariant deformation tensor b
    real(dp), dimension(:), intent(out) :: rlxmod

!   additional optional parameters (see type description)
    type(vemmod_t), intent(inout) :: vemmod


    real(dp), dimension(vemodel%ncompc) :: c, rlxmodc
    real(dp), dimension(2,2) :: bm, bminv, cm, rlxmodm
    type(vemmod_t) :: vemmodc

!   fill vemmodc

    vemmodc%compute_drlxmod = vemmod%compute_drlxmod
    allocate(vemmodc%drlxmod(vemodel%ncompc,vemodel%ncompc))
    vemmodc%dep_J = vemmod%dep_J
    if ( vemmodc%dep_J ) vemmodc%J = vemmod%J
    vemmodc%compute_drhsdJ = vemmod%compute_drhsdJ
    allocate(vemmodc%drlxmoddJ(vemodel%ncompc))

!   b-tensor and inverse

    bm = reshape ( b(1:4), [2,2], order=[2,1] )
    call matinv2 ( bm, bminv )

!   conformation tensor

    cm = matmul ( bm, transpose ( bm ) )
    c(1:3) = [ cm(1:2,1), cm(2,2) ]
    if ( vemodel%ncompb == 5 ) then
      c(4) = b(5)**2
    end if

!   rhs for the conformation tensor

    call rhs_relaxation_2D ( vemodel, mode, c, rlxmodc, vemmodc )

!   transform to rhs of the b-tensor: rhs(c).b^-T / 2

    rlxmodm = matmul ( vector_to_tensor2_symmetric ( ndim=2, A=rlxmodc(1:3) ), &
                         transpose(bminv) ) / 2
    rlxmod(1:4) = tensor2_to_vector ( ndim=2, A = rlxmodm )

    if ( vemodel%ncompb == 5 ) then
      rlxmod(5) = rlxmodc(4) / b(5) / 2
    end if

    if ( vemmod%compute_drlxmod ) then

!     Jacobian of relaxation terms

      if ( vemodel%ncompb == 5 ) then
        call compute_drlxmod_2D_with_zz
      else
        call compute_drlxmod_2D
      end if

    end if

!   TODO: vemmod%drlxmoddJ from vemmodc%drlmoddcdJ

    if ( vemodel%alam_model > 0 ) then

!     adapted lambda

      call rhs_adapted_lambda ( vemodel, mode, b, vonmises=vm_2D_sm_b, &
        vonmises_mm=vm_2D_mm_b, rlxmod=rlxmod, vemmod=vemmod )

    end if

    if ( subtractlinear ) then

!     Subtract the linear part.

      write(*,'(/2(a/))') &
        'Error in rhs_relaxation_2D_cb: subtractlinear', &
        'not yet implemented for the b-formulation'
      stop

    end if

  contains


!   internal subroutine for computing drlxmod (2D flow)

    subroutine compute_drlxmod_2D

      integer :: i, j, k, l
      real(dp) :: drlxten4(2,2,2,2), ten4(2,2,2,2)

      do j = 1, 2
        do k = 1, 2
          drlxten4(:,j,k,:) = - rlxmodm * bminv(j,k)
        end do
      end do

!     delta_ik b_jl

      do i = 1, 2
        ten4(i,:,i,:) = bm
        do k = i+1, 2
          ten4(i,:,k,:) = 0
          ten4(k,:,i,:) = 0
        end do
      end do

!     b_il delta_jk

      do k = 1, 2
        ten4(:,k,k,:) = ten4(:,k,k,:) + bm
      end do

!     dh(c)/dc * dc/db / 2

      ten4 = matrix_to_tensor4_left_symmetric ( ndim=2, &
                  A = matmul ( vemmodc%drlxmod, &
                  tensor4_to_matrix_left_symmetric ( ndim=2, A=ten4 ) ) ) / 2

!     ten4(i,m,k,l) * b^{-1}_jm

      do k = 1, 2
        do l = 1, 2
          drlxten4(:,:,k,l) = drlxten4(:,:,k,l) + &
                                matmul ( ten4(:,:,k,l), transpose(bminv) )
        end do
      end do

!     transform to matrix format

      vemmod%drlxmod = tensor4_to_matrix ( ndim=2, A=drlxten4 )

    end subroutine compute_drlxmod_2D


!   internal subroutine for computing drlxmod (2D flow + zz component)

    subroutine compute_drlxmod_2D_with_zz

      integer :: i, j, k, l
      integer, parameter :: i1(4)=[1,2,4,6], i2(2)=[3,5], i3(5)=[1,2,4,5,9]
      real(dp) :: bm_zz(3,3), bminv_zz(3,3), drlxmodc_zz(6,6)
      real(dp) :: rlxmodm_zz(3,3), drlxmod_zz(9,9)
      real(dp) :: drlxten4(3,3,3,3), ten4(3,3,3,3)

!     prepare 3D tensors/matrices from 2D with zz ones

      bm_zz(1:2,1:2) = bm
      bm_zz(1:2,3) = 0
      bm_zz(3,1:2) = 0
      bm_zz(3,3) = b(5)

      bminv_zz(1:2,1:2) = bminv
      bminv_zz(1:2,3) = 0
      bminv_zz(3,1:2) = 0
      bminv_zz(3,3) = 1/b(5)

      drlxmodc_zz(i1,i1) = vemmodc%drlxmod
      drlxmodc_zz(i2,i1) = 0
      drlxmodc_zz(i1,i2) = 0
      drlxmodc_zz(i2,i2) = 0

      rlxmodm_zz(1:2,1:2) = rlxmodm
      rlxmodm_zz(1:2,3) = 0
      rlxmodm_zz(3,1:2) = 0
      rlxmodm_zz(3,3) = rlxmod(5)

!     all 3D now

      do j = 1, 3
        do k = 1, 3
          drlxten4(:,j,k,:) = 0 - rlxmodm_zz * bminv_zz(j,k)
        end do
      end do

!     delta_ik b_jl

      do i = 1, 3
        ten4(i,:,i,:) = bm_zz
        do k = i+1, 3
          ten4(i,:,k,:) = 0
          ten4(k,:,i,:) = 0
        end do
      end do

!     b_il delta_jk

      do k = 1, 3
        ten4(:,k,k,:) = ten4(:,k,k,:) + bm_zz
      end do

!     dh(c)/dc * dc/db / 2

      ten4 = matrix_to_tensor4_left_symmetric ( ndim=3, &
                A = matmul ( drlxmodc_zz, &
                  tensor4_to_matrix_left_symmetric ( ndim=3, A=ten4 ) ) ) / 2

!     ten4(i,m,k,l) * b^{-1}_jm

      do k = 1, 3
        do l = 1, 3
          drlxten4(:,:,k,l) = drlxten4(:,:,k,l) + &
                                matmul ( ten4(:,:,k,l), transpose(bminv_zz) )
        end do
      end do

!     transform to matrix format

      drlxmod_zz = tensor4_to_matrix ( ndim=3, A=drlxten4 )
      vemmod%drlxmod = drlxmod_zz(i3,i3)

    end subroutine compute_drlxmod_2D_with_zz

  end subroutine rhs_relaxation_2D_cb


! Jacobian of relaxation term in the models

  subroutine drhs_relaxation_2D_b ( vemodel, mode, b, vemmod )

    use limits_m, only: USE_EGP_STRESS_TENSOR_FORM, &
                        USE_ADAP_LAMBDA_FOR_EGP_DET_STAB

    type(vemodel_t), intent(in) :: vemodel

!   the mode number
    integer, intent(in) :: mode

!   contravariant deformation tensor b
!   b(comp) is component comp
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:), intent(in) :: b

!   additional optional parameters (see type description)
    type(vemmod_t), intent(inout) :: vemmod


    integer :: i, j, k, l, flowtype
    real(dp) :: lambda, alpha, f, df, bpar, trc, c5, eta, tau(4), dtau(4,5)
    real(dp), dimension(2,2) :: bm, bminv, cm, cm1, bm1, taum, dtaum5
    real(dp), dimension(2,2,2,2) :: ten4, v4, w4, dtaum

    flowtype = vemodel%flowtype


    associate ( drlxmod => vemmod%drlxmod )


    select case ( vemodel%model )

    case(2)

!     Maxwell/Oldroyd

      lambda = vemodel%lambda(mode)

      bm = reshape ( b(1:4), [2,2], order=[2,1] )
      call matinv2 ( bm, bminv )

      do l = 1, 2
        do k = 1, 2
          do j = 1, 2
            do i = 1, 2
              ten4(i,j,k,l) = - bminv(l,i) * bminv(j,k)
            end do
          end do
        end do
      end do

      drlxmod(1:4,1:4) = - tensor4_to_matrix ( ndim=2, A=ten4 )
      do i = 1, 4
        drlxmod(i,i) = drlxmod(i,i) + 1
      end do

      if ( flowtype == 1 ) then
        drlxmod(1:4,5) = 0
        drlxmod(5,1:4) = 0
        drlxmod(5,5) = 1 + 1 / b(5)**2
      end if

      drlxmod = drlxmod / ( 2 * lambda )

    case(3)

!     Giesekus

      lambda = vemodel%lambda(mode)
      alpha  = vemodel%nonlin(1,mode)

      bm = reshape ( b(1:4), [2,2], order=[2,1] )
      call matinv2 ( bm, bminv )

!     bm1 = b - b^-T
      bm1 = bm - transpose ( bminv )

!     cm1 = (1-alpha) I + alpha c
      cm1 = alpha * matmul ( bm, transpose ( bm ) )
      cm1(1,1) = 1 - alpha + cm1(1,1)
      cm1(2,2) = 1 - alpha + cm1(2,2)

!     delta_ik b_jl

      do i = 1, 2
        v4(i,:,i,:) = bm
        do k = i+1, 2
          v4(i,:,k,:) = 0
          v4(k,:,i,:) = 0
        end do
      end do

!     b_il delta_jk

      do k = 1, 2
        v4(:,k,k,:) = v4(:,k,k,:) + bm
      end do

      v4 = alpha * v4

!     b^{-1}_li b^{-1}_jk

      do l = 1, 2
        do k = 1, 2
          do j = 1, 2
            do i = 1, 2
              w4(i,j,k,l) = bminv(l,i) * bminv(j,k)
            end do
          end do
        end do
      end do

!     delta_ik delta_jl

      do i = 1, 2
        do j = 1, 2
          w4(i,j,i,j) = w4(i,j,i,j) + 1
        end do
      end do

      do k = 1, 2
        do l = 1, 2
          ten4(:,:,k,l) = &
                    matmul( v4(:,:,k,l), bm1 ) + matmul( cm1, w4(:,:,k,l) )
        end do
      end do

      drlxmod(1:4,1:4) = tensor4_to_matrix ( ndim=2, A=ten4 )

      if ( flowtype == 1 ) then
        c5 = b(5)**2
        drlxmod(5,5) = ( 1 + c5 + alpha * (c5-1)*(3*c5+1) ) / c5
        drlxmod(1:4,5) = 0
        drlxmod(5,1:4) = 0
      end if

      drlxmod = drlxmod / ( 2 * lambda )

    case(20)

!     FENE-P

      lambda = vemodel%lambda(mode)
      bpar   = vemodel%nonlin(1,mode)

      bm = reshape ( b(1:4), [2,2], order=[2,1] )
      call matinv2 ( bm, bminv )

      cm = matmul ( bm, transpose ( bm ) )
      trc = cm(1,1) + cm(2,2) + b(5)**2
      f   = bpar / ( bpar + 3 - trc )
      df  = bpar / ( bpar + 3 - trc )**2

!     b^{-1}_li b^{-1}_jk

      do l = 1, 2
        do k = 1, 2
          do j = 1, 2
            do i = 1, 2
              w4(i,j,k,l) = bminv(l,i) * bminv(j,k)
            end do
          end do
        end do
      end do

!     f(trc) delta_ik delta_jl

      do i = 1, 2
        do j = 1, 2
          w4(i,j,i,j) = w4(i,j,i,j) + f
        end do
      end do

!     f'(trc) 2 bb

      do l = 1, 2
        do k = 1, 2
          w4(:,:,k,l) = w4(:,:,k,l) + 2 * df * bm * bm(k,l)
        end do
      end do

      drlxmod(1:4,1:4) = tensor4_to_matrix ( ndim=2, A=w4 )

      drlxmod(5,5) = 2 * df * b(5)**2 + f + 1 / b(5)**2
      drlxmod(5,1:4) = 2 * df * b(5) * b(1:4)
      drlxmod(1:4,5) = 2 * df * b(1:4) * b(5)

      drlxmod = drlxmod / ( 2 * lambda )

    case(23,24)

!     EGP

      lambda = vemodel%lambda(mode)

      if ( USE_EGP_STRESS_TENSOR_FORM ) then

        eta = vemodel%modulus(mode) * lambda

        call stress_viscoelastic_2D_single_mode_b ( vemodel, b, tau, mode, &
          vemmod )
        call dstress_viscoelastic_2D_single_mode_b ( vemodel, b, dtau, mode, &
          vemmod )

        taum = vector_to_tensor2_symmetric ( ndim=2, A=tau(1:3) )
        dtaum = matrix_to_tensor4_left_symmetric ( ndim=2, A=dtau(1:3,1:4) )
        dtaum5 = vector_to_tensor2_symmetric ( ndim=2, A=dtau(1:3,5) )

        bm = reshape ( b(1:4), [2,2], order=[2,1] )

        do k = 1, 2
          do l = 1, 2
            ten4(:,:,k,l) = matmul( dtaum(:,:,k,l), bm )
          end do
        end do

!       tau_ik delta_jl

        do j = 1, 2
          ten4(:,j,:,j) = ten4(:,j,:,j) + taum
        end do

        drlxmod(1:4,1:4) = tensor4_to_matrix ( ndim=2, A=ten4 )
        drlxmod(1:4,5) = tensor2_to_vector ( ndim=2, A=matmul( dtaum5, bm ) )
        drlxmod(5,1:4) = dtau(4,1:4) * b(5)
        drlxmod(5,5) = dtau(4,5) * b(5) + tau(4)
        drlxmod = drlxmod / ( 2 * eta )

      else

        bm = reshape ( b(1:4), [2,2], order=[2,1] )

!       cm1 = c - trc/3 I

        cm1 = matmul ( bm, transpose ( bm ) )
        trc = cm1(1,1) + cm1(2,2) + b(5)**2
        do i = 1, 2
          cm1(i,i) = cm1(i,i) - trc/3
        end do

!       delta_ik b_jl

        do i = 1, 2
          v4(i,:,i,:) = bm
          do k = i+1, 2
            v4(i,:,k,:) = 0
            v4(k,:,i,:) = 0
          end do
        end do

!       b_il delta_jk

        do k = 1, 2
          v4(:,k,k,:) = v4(:,k,k,:) + bm
        end do

!       - 2 * delta_ij b_kl / 3

        do i = 1, 2
          v4(i,i,:,:) = v4(i,i,:,:) - 2 * bm / 3
        end do

        do k = 1, 2
          do l = 1, 2
            ten4(:,:,k,l) = matmul( v4(:,:,k,l), bm )
          end do
        end do

!       cm1_ik delta_jl

        do j = 1, 2
          ten4(:,j,:,j) = ten4(:,j,:,j) + cm1
        end do

        drlxmod(1:4,1:4) = tensor4_to_matrix ( ndim=2, A=ten4 )
        drlxmod(1:4,5) = -2*b(5)/3*b(1:4)
        drlxmod(5,1:4) = -2*b(1:4)/3*b(5)
        drlxmod(5,5) = (7*b(5)**2 - trc)/3
!        drlxmod(5,5) = 3*b(5)**2 - trc/3 - 2*b(5)/3 * b(5)
        drlxmod = drlxmod / ( 2 * lambda )

        if ( vemodel%model == 24 ) drlxmod = drlxmod / vemmod%J

      end if

      if ( USE_ADAP_LAMBDA_FOR_EGP_DET_STAB ) then
!       Add Jacobian of stabilization term for determinant here
        do i = 1, 5
          drlxmod(i,:) = drlxmod(i,:) + ddetb_2D(b) * b(i) / ( 3 * lambda )
          drlxmod(i,i) = drlxmod(i,i) + ( detb_2D(b)- 1 ) / ( 3 * lambda )
!          drlxmod(i,:) = drlxmod(i,:) &
!                         + 2 * detb_2D(b) * ddetb_2D(b) * b(i) / ( 6 * lambda )
!          drlxmod(i,i) = drlxmod(i,i) + ( detb_2D(b)**2 - 1 ) / ( 6 * lambda )
        end do
      end if

    case default

      write(*,'(a,i0,3(/a)/)') &
        'Error in drhs_relaxation_2D_b: model not available: ', vemodel%model,&
        '  Set USE_CONFORMATION_MODEL = .true. in the limits_m model ', &
        '  to use a transformation of the conformation model formulation ', &
        '  to the contravariant deformation formulation.'
      stop

    end select


    end associate


  end subroutine drhs_relaxation_2D_b


! derivative of right-hand side wrt L for 2D and axisym. viscoelastic models.

  subroutine drhsdL_viscoelastic_2D_b ( vemodel, b, drhsdL, mode1, mode2, &
    mvemodel )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(b,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   contravariant deformation tensor b
!   b(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:,:,:), intent(in) :: b

!   derivative of the right-hand side wrt L of visco-elastic models:
!        .
!        b =  L * b - relaxation.
!
!   drhsdL(ip,i,j,mode) is derivative of component i of rhs with respect
!   to component j of L in point ip for mode number mode
!   column component sequence similar to the b tensor
!   row sequence ordering according to row-major storage of L:
!      2D:           [ Lxx, Lxy, Lyx, Lyy ]
!      axisymmetric: [ Lxx, Lxy, Lyx, Lyy, Ltt ]
!   where Ltt is the hoop strain rate.
    real(dp), dimension(:,:,:,:), intent(out) :: drhsdL

!   if present: limit the range of modes to mode1--mode2
    integer, intent(in), optional :: mode1, mode2


    integer :: mode, flowtype, ncompb, np, ip, modenr1, modenr2, i, j
    real(dp) :: bzz, xsi, fac
    real(dp), dimension(vemodel%ncompb) :: bvec
    real(dp), dimension(2,2) :: bm
    real(dp), dimension(4,4) :: dfL, dfL1
    real(dp), dimension(2,2,2,2) :: w4


    if ( vemodel%model == 0 ) then

      write(*,'(/3(a/))') &
        'Error in drhsdL_viscoelastic_2D_b: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 0 .and. vemodel%flowtype /= 1 ) then

      write(*,'(/a/2a/)') &
        'Error in drhsdL_viscoelastic_2D_b: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=0 (2D) or ', &
        'flowtype=1 (axisymmetric) '
      stop

    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(b,1) ) then
        write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_2D_b:', &
          ' dimension of mvemodel must be identical to size(b,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) .or. &
           any( mvemodel(:)%bvariant /= vemodel%bvariant ) ) then
        write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_2D_b:', &
          ' model, flowtype, bvariant in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( present(mode1) .and. present(mode2) ) then
      if ( mode1 < 1 .or. mode1 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_2D_b:', &
          ' mode1 out of range '
        stop
      end if
      if ( mode2 < 1 .or. mode2 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_2D_b:', &
          ' mode2 out of range '
        stop
      end if
      modenr1 = mode1
      modenr2 = mode2
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    if ( vemodel%bvariant >= 2 ) then
       write(*,'(/a/a/)') 'Error in drhsdL_viscoelastic_2D_b:', &
         ' derivative matrix not available for bvariant >= 2.'
      stop
    endif

    if ( any( vemodel%model == [14,17] ) ) then
      write(*,'(/a/a,i0/)') 'Error in drhsdL_viscoelastic_2D_b:', &
        ' derivative matrix not yet available for model = ', vemodel%model
      stop
    endif

    ncompb = vemodel%ncompb
    flowtype = vemodel%flowtype
    np = size(b,1)

    do mode = modenr1, modenr2

      do ip = 1, np

!       extract b

        bvec = b(ip,:ncompb,mode)

        bm = reshape ( bvec(1:4), [2,2], order=[2,1] )

        if ( flowtype == 1 ) bzz = bvec(5) ! axisymmetric

        do i = 1, 2
          w4(i,:,i,:) = transpose(bm)
          do j = i+1, 2
            w4(i,:,j,:) = 0
            w4(j,:,i,:) = 0
          end do
        end do

        dfL = tensor4_to_matrix ( ndim=2, A=w4 )

!       slip term

        if ( vemodel%ixsi > 0  ) then

!         models with slip

          if ( present(mvemodel) ) then
            xsi = mvemodel(ip)%nonlin(vemodel%ixsi,mode)
          else
            xsi = vemodel%nonlin(vemodel%ixsi,mode)
          end if

          do i = 1, 2
            w4(i,:,:,i) = transpose(bm)
            do j = i+1, 2
              w4(i,:,:,j) = 0
              w4(j,:,:,i) = 0
            end do
          end do

          dfL1 = tensor4_to_matrix ( ndim=2, A=w4 )

!         Gordon-Schowalter

          dfL = ( 1 - xsi/2 ) * dfL - xsi/2 * dfL1

        end if

        drhsdL(ip,1:4,1:4,mode) = dfL

        drhsdL(ip,5:ncompb,1:4,mode) = 0

        if ( flowtype == 1 ) then

          drhsdL(ip,1:4,5,mode) = 0
          if ( vemodel%ixsi > 0  ) then
            drhsdL(ip,5,5,mode) = ( 1 - xsi ) * bzz
          else
            drhsdL(ip,5,5,mode) = bzz
          end if
          drhsdL(ip,5,1:4,mode) = 0

        end if

!       replace L with L^d = L - (tr L)/3 I

        if ( vemodel%deviatoric ) then

          if ( vemodel%ixsi > 0  ) then
            fac = ( 1 - xsi ) / 3
          else
            fac = 1._dp / 3
          end if

          drhsdL(ip,1:5,1,mode) = drhsdL(ip,1:5,1,mode) - fac * b(ip,1:5,mode)
          drhsdL(ip,1:5,4,mode) = drhsdL(ip,1:5,4,mode) - fac * b(ip,1:5,mode)
          if ( flowtype == 1 ) then
            drhsdL(ip,1:5,5,mode) = drhsdL(ip,1:5,5,mode) - fac * b(ip,1:5,mode)
          end if

        end if

      end do

    end do

  end subroutine drhsdL_viscoelastic_2D_b


! Newton-Raphson iteration terms for steady state time-discretization for
! 2D and axisymmetric viscoelastic models.

  subroutine NRtd_steady_viscoelastic_2D_b ( vemodel, tdpar, b, Hb, dHb, &
    mode1, mode2 )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   Parameters regarding time discretization
    type(tdpar_t), intent(in) :: tdpar

!   contravariant deformation tensor b
!   b(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:,:,:), intent(in) :: b

!   steady state time-discretization term (see note 519):
!
!        H(b) ... = 0
!
!   where H(b) depending on the BDF discretization scheme:
!
!     BDF1:  H = (b-b.RT)/deltat = (b-V)/deltat
!
!     BDF2:  H = (3b-4b.RT+b.RT^2)/2/deltat = (3b-4V+V.RT)/2/deltat
!
!   where the polar decomposition of is written as b=V.R and RT = transpose(R).
!
!   Hb(ip,comp,mode) is component comp in point ip for mode number mode
!   components sequence similar to the contravariant deformation tensor b
    real(dp), dimension(:,:,:), intent(out) :: Hb

!   Jacobian of the H(b) term (see note 519)
!   dHb(ip,i,j,mode) is derivative of component i with respect to component j
!   in point ip for mode number mode
!   component sequence similar to the contravariant deformation tensor b
    real(dp), dimension(:,:,:,:), intent(out) :: dHb

!   if present: limit the range of modes to mode1--mode2
    integer, intent(in), optional :: mode1, mode2

    integer :: mode, flowtype, ncompb, np, ip, modenr1, modenr2, i, j, k, l, &
      alpha, beta
    real(dp) :: fac
    real(dp), dimension(2,2,2,2) :: w4, ften4, dVdb, diad4, dHbten4, dKten4

    real(dp), dimension(1,vemodel%ncompb) :: bvec
    real(dp), dimension(1,vemodel%ncompc) :: Vvec
    real(dp), dimension(1,2) :: Veigval1
    real(dp), dimension(1,2,2) :: RT1, Veigvec1
    real(dp), dimension(2) :: Veigval
    real(dp), dimension(2,2) :: bten, RT, Vten, Hbten, Veigvec, bteninv, Vbinv


    if ( vemodel%model == 0 ) then

      write(*,'(/3(a/))') &
        'Error in NRtd_steady_viscoelastic_2D_b: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 0 .and. vemodel%flowtype /= 1 ) then

      write(*,'(/a/2a/)') &
        'Error in NRtd_steady_viscoelastic_2D_b: incorrect value of flowtype.',&
        'Call create_viscoelastic_model with flowtype=0 (2D) or ', &
        'flowtype=1 (axisymmetric) '
      stop

    end if

    if ( present(mode1) .and. present(mode2) ) then
      if ( mode1 < 1 .or. mode1 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in NRtd_steady_viscoelastic_2D_b:', &
          ' mode1 out of range '
        stop
      end if
      if ( mode2 < 1 .or. mode2 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in NRtd_steady_viscoelastic_2D_b:', &
          ' mode2 out of range '
        stop
      end if
      modenr1 = mode1
      modenr2 = mode2
    else if ( present(mode1) .or. present(mode2) ) then
      write(*,'(/a/a/)') 'Error in NRtd_steady_viscoelastic_2D_b:', &
        ' mode1 and mode2 need to be both present'
      stop
    else
      modenr1 = 1
      modenr2 = vemodel%nmodes
    end if

    if ( vemodel%bvariant /= 1 ) then
      write(*,'(/a/a,i0/)') 'Error in NRtd_steady_viscoelastic_2D_b:', &
        ' Only bvariant = 1 can be used. bvariant = ', vemodel%bvariant
      stop
    endif

    ncompb = vemodel%ncompb
    flowtype = vemodel%flowtype
    np = size(b,1)

    do mode = modenr1, modenr2

      do ip = 1, np

!       extract b

        bvec(1,:) = b(ip,:ncompb,mode)

        bten = vector_to_tensor2 ( ndim=2, A=bvec(1,1:4) )

!       compute V and get eigenvectors/values

        call sqrtc_2D_b ( bvec, sqrtc=Vvec, RT=RT1, Veigval=Veigval1, &
          Veigvec=Veigvec1 )

        RT=RT1(1,:,:); Veigval=Veigval1(1,:); Veigvec=Veigvec1(1,:,:)

        Vten = vector_to_tensor2_symmetric ( ndim=2, A=Vvec(1,1:3) )

!       compute H(b)

        select case ( tdpar%method )
        case (1) ! BDF1
          Hbten = (bten - Vten)/tdpar%timestep
        case (2) ! BDF2
          Hbten = (3*bten - 4*Vten + matmul(Vten,RT))/(2*tdpar%timestep)
        case default
          call errormsg_case_default ( 'NRtd_steady_viscoelastic_2D_b', &
            'tdpar%method', int_value=tdpar%method )
        end select

        Hb(ip,1:4,mode) = tensor2_to_vector ( ndim=2, A=Hbten )

        if ( flowtype == 1 ) then
          Hb(ip,5,mode) = 0
        else
          Hb(ip,5:ncompb,mode) = 0
        end if

!       compute dc/db

!       delta_ik b_jl

        do i = 1, 2
          w4(i,:,i,:) = bten
          do k = i+1, 2
            w4(i,:,k,:) = 0
            w4(k,:,i,:) = 0
          end do
        end do

!       b_il delta_jk

        do k = 1, 2
          w4(:,k,k,:) = w4(:,k,k,:) + bten
        end do

!       compute f^(kl)_alpha,beta = n_alpha.dc/db_kl.n_beta

        do k = 1, 2
          do l = 1, 2
            ften4(:,:,k,l) = &
                 matmul ( transpose(Veigvec), matmul(w4(:,:,k,l),Veigvec) )
          end do
        end do

!       compute n_alpha n_beta diadics

        do alpha = 1, 2
          do beta = 1, 2
            do i = 1, 2
              do j = 1, 2
                diad4(alpha,beta,i,j) = Veigvec(i,alpha) * Veigvec(j,beta)
              end do
            end do
          end do
        end do

!       compute dVdb

        dVdb = 0

        do alpha = 1, 2
          do beta = 1, 2
            fac = 1 / ( Veigval(alpha) + Veigval(beta) )
            do k = 1, 2
              do l = 1, 2
                dVdb(:,:,k,l) = dVdb(:,:,k,l) + &
                   ften4(alpha,beta,k,l) * fac * diad4(alpha,beta,:,:)
              end do
            end do
          end do
        end do

!       compute dH(b)

        select case ( tdpar%method )

        case (1) ! BDF1

          dHbten4 = - dVdb

!         delta_ik delta_jl

          do i = 1, 2
            do j = 1, 2
              dHbten4(i,j,i,j) = dHbten4(i,j,i,j) + 1
            end do
          end do
          dHbten4 = dHbten4 / tdpar%timestep

        case (2) ! BDF2

!         Compute dK/db with K=V.R^T=V.b^-1.V

!         term I: dV/db_kl.R^T

          do k = 1, 2
            do l = 1, 2
              dKten4(:,:,k,l) = matmul( dVdb(:,:,k,l), RT )
            end do
          end do

!         term II: (V.b^-1)_ik R^T_lj

          call matinv2 ( bten, bteninv )
          Vbinv = matmul ( Vten, bteninv )

          do i = 1, 2
            do j = 1, 2
              do k = 1, 2
                do l = 1, 2
                  dKten4(i,j,k,l) = dKten4(i,j,k,l) - Vbinv(i,k)*RT(l,j)
                end do
              end do
            end do
          end do

!         term III: V.b^-1.dV/db_kl

          do k = 1, 2
            do l = 1, 2
              dKten4(:,:,k,l) = dKten4(:,:,k,l) + matmul( Vbinv, dVdb(:,:,k,l) )
            end do
          end do

!         dHb

          dHbten4 = - 4 * dVdb + dKten4

!         3 * delta_ik delta_jl

          do i = 1, 2
            do j = 1, 2
              dHbten4(i,j,i,j) = dHbten4(i,j,i,j) + 3
            end do
          end do

          dHbten4 = dHbten4 / ( 2 * tdpar%timestep )

        case default

          call errormsg_case_default ( 'NRtd_steady_viscoelastic_2D_b', &
            'tdpar%method', int_value=tdpar%method )

        end select

        dHb(ip,1:4,1:4,mode) = tensor4_to_matrix ( ndim=2, A=dHbten4 )

        if ( flowtype == 1 ) then
          dHb(ip,5,1:5,mode) = 0
          dHb(ip,1:4,5,mode) = 0
        else
          dHb(ip,5:ncompb,:,mode) = 0
          dHb(ip,:,5:ncompb,mode) = 0
        end if

      end do

    end do

  end subroutine NRtd_steady_viscoelastic_2D_b


! von Mises equivalent shear stress for a single mode

  subroutine vm_2D_sm_b ( vemodel, b, mode, vemmod, vm, dvm, dvmdJ )

    type(vemodel_t), intent(in) :: vemodel
!   b(:ncompb)
    real(dp), dimension(:), intent(in) :: b
    integer, intent(in) :: mode
    type(vemmod_t), intent(in) :: vemmod
    real(dp), intent(out), optional :: vm
!   dvm(:ncompb)
    real(dp), dimension(:), intent(out), optional :: dvm
    real(dp), intent(out), optional :: dvmdJ

    real(dp) :: tau(vemodel%ncompt)
    real(dp) :: dtau(vemodel%ncompt,vemodel%ncompb), dtaudJ(vemodel%ncompt)

    call stress_viscoelastic_2D_single_mode_b ( vemodel, b, tau, mode, vemmod )

    if ( present(vm) ) vm = vonmises_2D ( tau )

    if ( present(dvm) ) then

      call dstress_viscoelastic_2D_single_mode_b ( vemodel, b, dtau, mode, &
        vemmod )

      dvm = matmul ( dvonmises_2D(tau), dtau )

    end if

    if ( present(dvmdJ) ) then

      call dstressdJ_viscoelastic_2D_single_mode_b ( vemodel, b, dtaudJ, mode, &
        vemmod )

      dvmdJ = dot_product ( dvonmises_2D(tau), dtaudJ )

    end if

  end subroutine vm_2D_sm_b


! von Mises equivalent shear stress based on multiple modes

  subroutine vm_2D_mm_b ( vemodel, vemmod, vm, dvm, dvmdJ )

    type(vemodel_t), intent(in) :: vemodel
    type(vemmod_t), intent(in) :: vemmod
    real(dp), intent(out), optional :: vm
!   dvm(:ncompb,:nmodes), only mode range vmmode1:vmmode2 is computed
    real(dp), dimension(:,:), intent(out), optional :: dvm
    real(dp), intent(out), optional :: dvmdJ

    integer :: mode, vmmode1, vmmode2
    real(dp), dimension(vemodel%ncompt) :: tau, tsm, dvm2Dtau
    real(dp) :: dtau(vemodel%ncompt,vemodel%ncompb)
    real(dp) :: dtaudJsm(vemodel%ncompt), dtaudJ(vemodel%ncompt)

    vmmode1 = vemmod%vmmode1
    vmmode2 = vemmod%vmmode2

    tau = 0

    do mode = vmmode1, vmmode2

      call stress_viscoelastic_2D_single_mode_b ( vemodel, vemmod%b(:,mode), &
        tsm, mode, vemmod )

      tau = tau + tsm

    end do

    if ( present(vm) ) vm = vonmises_2D(tau)

    if ( present(dvm) ) then

      dvm2Dtau = dvonmises_2D(tau)

      do mode = vmmode1, vmmode2

        call dstress_viscoelastic_2D_single_mode_b ( vemodel, &
          vemmod%b(:,mode), dtau, mode, vemmod )

        dvm(:,mode) = matmul ( dvm2Dtau, dtau )

      end do

    end if

    if ( present(dvmdJ) ) then

      dtaudJ = 0

      do mode = vmmode1, vmmode2

        call dstressdJ_viscoelastic_2D_single_mode_b ( vemodel, &
          vemmod%b(:,mode), dtaudJsm, mode, vemmod )

        dtaudJ = dtaudJ + dtaudJsm

      end do

      dvmdJ = dot_product ( dvonmises_2D(tau), dtaudJ )

    end if

  end subroutine vm_2D_mm_b


! stress tensor for 2D and axisymmetric viscoelastic models (single mode)

  subroutine stress_viscoelastic_2D_single_mode_b ( vemodel, b, tau, mode, &
    vemmod )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   contravariant deformation tensor b
!   b(comp) is component comp
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain

    real(dp), dimension(:), intent(in) :: b

!   stress tensor tau of the viscoelastic model
!   tau(comp) is component comp
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop strain
    real(dp), dimension(:), intent(out) :: tau

!   the mode number
    integer, intent(in), optional :: mode

!   additional optional parameters (see type description)
!   NOTE: this argument is set to optional for compatibility reasons.
!   It is checked to be present, when needed.
    type(vemmod_t), intent(in), optional :: vemmod


    integer :: ncompb
    real(dp) :: bm(2,2), cm(2,2), c(vemodel%ncompc), bm5


    ncompb = vemodel%ncompb

    bm = reshape ( b(1:4), [2,2], order=[2,1] )
    if ( ncompb == 5 ) bm5 = b(5)

    if ( vemodel%bvariant == 4 ) then
!     Cholesky log
      bm(1,1) = exp(bm(1,1))
      bm(2,2) = exp(bm(2,2))
      if ( ncompb == 5 ) bm5 = exp(bm5)
    end if

    cm = matmul ( bm, transpose ( bm ) )
    c(1:3) = [ cm(1:2,1), cm(2,2) ]
    if ( ncompb == 5 ) then
      c(4) = bm5**2
    end if

    call stress_viscoelastic_2D_single_mode ( vemodel, c, tau, mode, vemmod )

  end subroutine stress_viscoelastic_2D_single_mode_b


! stress tensor for 2D and axisymmetric viscoelastic models.

  subroutine stress_viscoelastic_2D_b ( vemodel, b, tau, mode, mode1, mode2, &
    mvemodel, vemopt )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(b,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   contravariant deformation tensor b
!   b(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain

    real(dp), dimension(:,:,:), intent(in) :: b

!   stress tensor tau of the viscoelastic model
!   tau(ip,comp) is component comp in point ip
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop stress
    real(dp), dimension(:,:), intent(out) :: tau

!   when present only a single mode contribution is computed
    integer, intent(in), optional :: mode

!   if present: limit the range of modes to mode1--mode2
    integer, intent(in), optional :: mode1, mode2

!   if present: control the additional optional parameters below
    type(vemopt_t), intent(in), optional :: vemopt


    integer :: modenr1, modenr2, imode, np, ip, ncompb
    real(dp) :: tsm(vemodel%ncompt), bm(vemodel%ncompb)
    type(vemmod_t) :: vemmod


    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(b,1) ) then
        write(*,'(/a/a/)') 'Error in stress_viscoelastic_2D_b:', &
          ' dimension of mvemodel must be identical to size(b,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) .or. &
           any( mvemodel(:)%bvariant /= vemodel%bvariant ) ) then
        write(*,'(/a/a/)') 'Error in stress_viscoelastic_2D_b:', &
          ' model, flowtype, bvariant in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    call check ( vemodel, vemopt, 'stress_viscoelastic_2D_b', &
      components=['dep_J'] )

    if ( present(vemopt) ) then

!     copy logicals from vemopt to vemmod
      vemmod%dep_J = vemopt%dep_J

    end if

    np = size(b,1)
    ncompb = vemodel%ncompb

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

        bm = b(ip,1:ncompb,imode)

        if ( present(mvemodel) ) then
          call stress_viscoelastic_2D_single_mode_b ( mvemodel(ip), &
            bm, tsm, imode, vemmod )
        else
          call stress_viscoelastic_2D_single_mode_b ( vemodel, &
            bm, tsm, imode, vemmod )
        end if

        tau(ip,:) = tau(ip,:) + tsm

      end do

    end do

  end subroutine stress_viscoelastic_2D_b


! Jacobian of stress tensor for 2D viscoelastic models (single mode)

  subroutine dstress_viscoelastic_2D_single_mode_b ( vemodel, b, dtau, mode, &
    vemmod )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   contravariant deformation tensor b
!   b(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:), intent(in) :: b

!   derivative of stress tensor tau of the viscoelastic model with respect
!   to the contravariant deformation tensor b
!   dtau(i,j) is the derivative of stress component i wrt to component j of b
!     comp = 1 tauxx
!     comp = 2 tauxy
!     comp = 3 tauyy
!     comp = 4 tauzz when non-zero or
!     comp = 4 tautt axisymmetric hoop strain
    real(dp), dimension(:,:), intent(out) :: dtau

!   the mode number
    integer, intent(in), optional :: mode

!   additional optional parameters (see type description)
!   NOTE: this argument is set to optional for compatibility reasons.
!   It is checked to be present, when needed.
    type(vemmod_t), intent(in), optional :: vemmod


    integer :: i, k, ncompb
    real(dp) :: bm(2,2), cm(2,2), c(vemodel%ncompc), ten4(2,2,2,2)
    real(dp) :: dtauc(size(dtau,1),vemodel%ncompc)
    real(dp) :: wrk(vemodel%ncompc,vemodel%ncompb)


    if ( vemodel%bvariant >= 2 ) then
      write(*,'(/a/a/)') 'Error in dstress_viscoelastic_2D_single_mode_b:', &
        ' Jacobian matrix not available for bvariant >= 2.'
      stop
    endif

    ncompb = vemodel%ncompb

    bm = reshape ( b(1:4), [2,2], order=[2,1] )

    cm = matmul ( bm, transpose ( bm ) )
    c(1:3) = [ cm(1:2,1), cm(2,2) ]
    if ( ncompb == 5 ) then
      c(4) = b(5)**2
    end if

    call dstress_viscoelastic_2D_single_mode ( vemodel, c, dtauc, mode, vemmod )

!   delta_ik b_jl

    do i = 1, 2
      ten4(i,:,i,:) = bm
      do k = i+1, 2
        ten4(i,:,k,:) = 0
        ten4(k,:,i,:) = 0
      end do
    end do

!   b_il delta_jk

    do k = 1, 2
      ten4(:,k,k,:) = ten4(:,k,k,:) + bm
    end do

    if ( ncompb == 5 ) then
      wrk(1:3,1:4) = tensor4_to_matrix_left_symmetric(ndim=2,A=ten4)
      wrk(4,1:4) = 0
      wrk(1:3,5) = 0
      wrk(4,5) = 2 * b(5)
      dtau = matmul ( dtauc, wrk )
    else
      dtau = matmul ( dtauc, tensor4_to_matrix_left_symmetric(ndim=2,A=ten4) )
    end if

  end subroutine dstress_viscoelastic_2D_single_mode_b


! Jacobian of stress tensor for 2D and axisymmetric viscoelastic models
! (single mode) in multiple points with optional variable coefficients

  subroutine dstress_viscoelastic_2D_b ( vemodel, b, dtau, mode, mvemodel, &
    vemopt )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(b,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   contravariant deformation tensor b
!   b(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:,:,:), intent(in) :: b

!   derivative of the stress tensor tau of the viscoelastic model
!   dtau(i,j,k) is the derivative in point i of stress component j wrt to
!   component k of b
!     j = 1 tauxx
!     j = 2 tauxy
!     j = 3 tauyy
!     j = 4 tauzz when non-zero or
!     j = 4 tautt axisymmetric hoop stress
    real(dp), dimension(:,:,:), intent(out) :: dtau

!   mode computed
    integer, intent(in) :: mode

!   if present: control the additional optional parameters below
    type(vemopt_t), intent(in), optional :: vemopt


    integer :: np, ip, ncompb
    real(dp) :: dtsm(vemodel%ncompt,vemodel%ncompb), bm(vemodel%ncompb)
    type(vemmod_t) :: vemmod


    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(b,1) ) then
        write(*,'(/a/a/)') 'Error in dstress_viscoelastic_2D_b:', &
          ' dimension of mvemodel must be identical to size(b,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) .or. &
           any( mvemodel(:)%bvariant /= vemodel%bvariant ) ) then
        write(*,'(/a/a/)') 'Error in dstress_viscoelastic_2D_b:', &
          ' model, flowtype, bvariant in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    call check ( vemodel, vemopt, 'dstress_viscoelastic_2D_b', &
      components=['dep_J'] )

    if ( present(vemopt) ) then

!     copy logicals from vemopt to vemmod
      vemmod%dep_J = vemopt%dep_J

    end if

    np = size(b,1)
    ncompb = vemodel%ncompb

    do ip = 1, np

!     set J

      if ( vemmod%dep_J ) vemmod%J = vemopt%J(ip)

      bm = b(ip,1:ncompb,mode)

      if ( present(mvemodel) ) then
        call dstress_viscoelastic_2D_single_mode_b ( mvemodel(ip), &
          bm, dtsm, mode, vemmod )
      else
        call dstress_viscoelastic_2D_single_mode_b ( vemodel, &
          bm, dtsm, mode, vemmod )
      end if

      dtau(ip,:,:) = dtsm

    end do

  end subroutine dstress_viscoelastic_2D_b


! Jacobian with respect to J of stress tensor for 2D viscoelastic models
! (single mode) in multiple points with optional variable coefficients

  subroutine dstressdJ_viscoelastic_2D_b ( vemodel, b, dtaudJ, mode, mode1, &
    mode2, mvemodel, vemopt )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(b,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   contravariant deformation tensor b
!   b(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:,:,:), intent(in) :: b

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
      if ( size(mvemodel,1) /= size(b,1) ) then
        write(*,'(/a/a/)') 'Error in dstressdJ_viscoelastic_2D_b:', &
          ' dimension of mvemodel must be identical to size(b,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in dstressdJ_viscoelastic_2D_b:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    call check ( vemodel, vemopt, 'dstressdJ_viscoelastic_2D_b', &
      components=['dep_J'] )

    if ( present(vemopt) ) then

!     copy logicals from vemopt to vemmod
      vemmod%dep_J = vemopt%dep_J

    end if

    np = size(b,1)

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
          call dstressdJ_viscoelastic_2D_single_mode_b ( mvemodel(ip), &
            b(ip,:,imode), dtsm, imode, vemmod )
        else
          call dstressdJ_viscoelastic_2D_single_mode_b ( vemodel, &
            b(ip,:,imode), dtsm, imode, vemmod )
        end if

        dtaudJ(ip,:) = dtaudJ(ip,:) + dtsm

      end do

    end do

  end subroutine dstressdJ_viscoelastic_2D_b


! Jacobian of stress tensor wrt J for 2D viscoelastic models (single mode)

  subroutine dstressdJ_viscoelastic_2D_single_mode_b ( vemodel, b, dtaudJ, &
    mode, vemmod )

    type(vemodel_t), intent(in) :: vemodel

!   contravariant deformation tensor b
!   b(comp) is component comp
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:), intent(in) :: b

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

      call stress_viscoelastic_2D_single_mode_b ( vemodel, b, tau, mode, &
        vemmod )
      dtaudJ = - tau / vemmod%J

    case default

      dtaudJ = 0

    end select

  end subroutine dstressdJ_viscoelastic_2D_single_mode_b


! conformation tensor for 2D and axisymmetric viscoelastic models.

  subroutine conformation_2D_b ( b, c, bvariant )

!   contravariant deformation tensor b in c = b*b^T
!   b(ip,comp) is component comp in point ip
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:,:), intent(in) :: b

!   conformation tensor c of the viscoelastic model
!   c(ip,comp) is component comp in point ip
!     comp = 1 cxx
!     comp = 2 cxy
!     comp = 3 cyy
!     comp = 4 czz when non-zero or
!     comp = 4 ctt axisymmetric hoop strain
    real(dp), dimension(:,:), intent(out) :: c

!   Optionally the b variant can be given.
!   Only needed for bvariant=4 (Cholesky-log)
!   default=1
    integer, intent(in), optional :: bvariant

    integer :: np, ip, ncompb, lbvariant
    real(dp) :: cm(2,2), bm(2,2), bm5


    lbvariant = set_optional ( variable=bvariant, default=1 )

    np = size(b,1)
    ncompb = size(b,2)

    do ip = 1, np

!     extract value of b in one point and store in matrix

      bm = reshape ( b(ip,1:4), [2,2], order=[2,1] )
      if ( ncompb == 5 ) bm5 = b(ip,5)

      if ( lbvariant == 4 ) then
!       Cholesky log
        bm(1,1) = exp(bm(1,1))
        bm(2,2) = exp(bm(2,2))
        if ( ncompb == 5 ) bm5 = exp(bm5)
      end if

!     compute c = b*b^T

      cm = matmul ( bm, transpose ( bm ) )

      c(ip,1:3) = [ cm(1:2,1), cm(2,2) ]
      if ( ncompb == 5 ) then
        c(ip,4) = bm5**2
      end if

    end do

  end subroutine conformation_2D_b


! square root of conformation from the b tensor for 2D and axisymmetric
! viscoelastic models. Optionally compute R from the polar decomposition
! b = V.R, with V=sqrt(b.b^T)=sqrt(c).

  subroutine sqrtc_2D_b ( b, sqrtc, RT, bvariant, Veigval, Veigvec )

!   INPUT:
!     contravariant deformation tensor b in c = b*b^T
!       b(ip,comp) is component comp in point ip
!       comp = 1 bxx
!       comp = 2 bxy
!       comp = 3 byx
!       comp = 4 byy
!       comp = 5 bzz when non-zero or
!       comp = 5 btt axisymmetric hoop strain
!   OUTPUT:
!     if present(sqrtc)
!       b is not modified
!     otherwise it contains the V-tensor in unsymmetric storage:
!       comp = 1 Vxx
!       comp = 2 Vxy
!       comp = 3 Vyx=Vxy
!       comp = 4 Vyy
!       comp = 5 Vzz when non-zero or
!       comp = 5 Vtt axisymmetric hoop strain
    real(dp), dimension(:,:), intent(inout) :: b

!   If present: the symmetric V tensor = sqrt(b.b^T) = sqrt(c)
!   sqrtc(ip,comp) is component comp in point ip
!     comp = 1 Vxx
!     comp = 2 Vxy
!     comp = 3 Vyy
!     comp = 4 Vzz when non-zero or
!     comp = 4 Vtt axisymmetric hoop strain
!   otherwise the V tensor is output via the argument b.
    real(dp), dimension(:,:), intent(out), optional :: sqrtc

!   optional (transpose of) the rotation tensor from polar decomposition
    real(dp), dimension(size(b,1),2,2), intent(out), optional :: RT

!   optional eigenvalues and eigenvectors of the symmetric tensor V
    real(dp), dimension(size(b,1),2), intent(out), optional :: Veigval
    real(dp), dimension(size(b,1),2,2), intent(out), optional :: Veigvec

!   Optionally the b variant can be given.
!   Only needed for bvariant=4 (Cholesky-log)
!   default=1
    integer, intent(in), optional :: bvariant


    integer  :: np, ip, ncompb, lbvariant
    real(dp) :: bpn(size(b,2)), cpn(3)
    real(dp) :: Vdiag(2), eigvalue(2), eigv(2,2)
    real(dp) :: bpnm(2,2), bpnminv(2,2)
    real(dp) :: V(3), Vm(2,2)

    lbvariant = set_optional ( variable=bvariant, default=1 )

    np = size(b,1)
    ncompb = size(b,2)

    do ip = 1, np

!     extract value of b and compute c in one point

      bpn = b(ip,:)

      if ( lbvariant == 4 ) then
!       Cholesky log
        bpn(1) = exp(bpn(1))
        bpn(4) = exp(bpn(4))
        if ( ncompb == 5 ) bpn(5) = exp(bpn(5))
      end if

      cpn(1) = bpn(1)**2 + bpn(2)**2
      cpn(2) = bpn(1)*bpn(3) + bpn(2)*bpn(4)
      cpn(3) = bpn(3)**2 + bpn(4)**2

!     compute sqrt(cpn)

      call eig2x2 ( cpn, eigvalue, eigv )

      if ( present(Veigvec) ) Veigvec(ip,:,:) = eigv

      if ( all( eigvalue >= 0 ) ) then
        Vdiag = sqrt(eigvalue)
      else
        write(*,*) ' Warning sqrtc_2D_b: negative eigvalue = ',  eigvalue
        where ( eigvalue >= 0 )
          Vdiag = sqrt(eigvalue)
        else where
          Vdiag = 10*tiny(1._dp)  ! set to very small positive value
        end where
      end if

      if ( present(Veigval) ) Veigval(ip,:) = Vdiag

      call inveig2x2 ( V, Vdiag, eigv )

      Vm = reshape ( [ V(1:2), V(2), V(3) ], [2,2], order=[2,1] )

      if ( present(sqrtc) ) then
!       fill sqrtc with V
        sqrtc(ip,1:3) = V
        if ( ncompb == 5 ) sqrtc(ip,4) = bpn(5)
      else
!       fill b with V
        b(ip,1:4) = reshape ( Vm, [4] )
        if ( ncompb == 5 ) then
          b(ip,5) = bpn(5)
        end if
      end if

      if ( present(RT) ) then

!       determine (transpose of) rotation tensor R

        bpnm = reshape ( bpn(1:4), [2,2], order=[2,1] )
        call matinv2 ( bpnm, bpnminv )
        RT(ip,:,:) = matmul ( bpnminv, Vm ) ! R^T = b^{-1}.V

      end if

    end do

  end subroutine sqrtc_2D_b


! rotate contravariant deformation for 2D and axisymmetric
! viscoelastic models.

  subroutine rotate_2D_b ( b, RT )

!   INPUT
!     contravariant deformation tensor b in c = b*b^T
!     b(ip,comp) is component comp in point ip
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
!   OUTPUT
!     components of b.R^T
    real(dp), dimension(:,:), intent(inout) :: b

!   (transpose of) rotation tensor R
    real(dp), dimension(size(b,1),2,2), intent(in) :: RT

    integer :: np, ip
    real(dp) :: bpn(2,2)

    np = size(b,1)

    do ip = 1, np

!     extract value of b in one point

      bpn = reshape ( b(ip,1:4), [2,2], order=[2,1] )

!     rotate b according to b' = b*R^T

      b(ip,1:4) = reshape ( transpose ( matmul ( bpn, RT(ip,:,:) ) ), [4] )

    end do

  end subroutine rotate_2D_b


! performs a direct calculation of the inverse of a 2x2 matrix

  subroutine matinv2 ( a, b )

    real(dp), intent(in) :: a(2,2)   ! Matrix
    real(dp), intent(out) :: b(2,2)  ! Inverse matrix
    real(dp) :: detinv               ! Inverse determinant

!   calculate the inverse determinant of the matrix
    detinv = 1/(a(1,1)*a(2,2) - a(1,2)*a(2,1))

!   calculate the inverse of the matrix
    b(1,1) =  detinv * a(2,2)
    b(2,1) = -detinv * a(2,1)
    b(1,2) = -detinv * a(1,2)
    b(2,2) =  detinv * a(1,1)

  end subroutine matinv2


! skew norm contravariant deformation tensor for 2D and axisymmetric
! viscoelastic models. Not applicable to Cholesly-log.

  function skew_norm_2D_b ( b )

    use math_defs_m, only: sumsq

!   contravariant deformation tensor b in c = b*b^T
!   b(ip,comp) is component comp in point ip
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:,:), intent(in) :: b

    real(dp), dimension(size(b,1)) :: skew_norm_2D_b

    skew_norm_2D_b = ( b(:,2) - b(:,3) )**2 / sumsq( b(:,1:4), dim=2 ) / 2

  end function skew_norm_2D_b


! determinant of contravariant deformation tensor b

  function detb_2D ( b )

!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:), intent(in) :: b
    real(dp) :: detb_2D

    detb_2D = b(1)*b(4) - b(2)*b(3)
    if ( size(b) == 5 ) detb_2D = detb_2D * b(5)

  end function detb_2D


! Jacobian of determinant of contravariant deformation tensor b

  function ddetb_2D ( b )

!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 byx
!     comp = 4 byy
!     comp = 5 bzz when non-zero or
!     comp = 5 btt axisymmetric hoop strain
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(size(b)) :: ddetb_2D

    ddetb_2D(1:4) = [ b(4), -b(3), -b(2), b(1) ]
    if ( size(b) == 5 ) then
      ddetb_2D(1:4) = ddetb_2D(1:4) * b(5)
      ddetb_2D(5) = b(1)*b(4) - b(2)*b(3)
    end if

  end function ddetb_2D


! Right-hand side for plastic strain evolution of 2D models. Optional Jacobian.

  subroutine rhs_plastic_strain_2D_b ( vemodel, b, gammap, rhs_gammap, &
    mode1, mode2, mvemodel, vemopt_gammap )

!   The definition of the viscoelastic model
    type(vemodel_t), intent(in) :: vemodel

!   If present, mvemodel(ip) provides alternative data (modulus,lambda,nonlin)
!   for point ip. Note, that
!    - size(mvemodel,1)==size(b,1) must be valid.
!    - all mvemodel(:) must only have alternative data (modulus,lambda,nonlin),
!      all other components (model, .. ) must be identical to vemodel!
    type(vemodel_t), dimension(:), intent(in), optional :: mvemodel

!   contravariant deformation tensor b
!   b(ip,comp,mode) is component comp in point ip for mode number mode
!     comp = 1 bxx
!     comp = 2 bxy
!     comp = 3 bxz
!     comp = 4 byx
!     comp = 5 byy
!     comp = 6 byz
!     comp = 7 bzx
!     comp = 8 bzy
!     comp = 9 bzz
    real(dp), dimension(:,:,:), intent(in) :: b

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
    integer :: mode, ncompb, np, ip, modenr1, modenr2, m, vm1, vm2
    type(vemmod_t) :: vemmod


    if ( vemodel%model == 0 ) then

      write(*,'(/3(a/))') &
        'Error in rhs_plastic_strain_2D_b: model not set.', &
        'Call create_viscoelastic_model and fill material parameters in', &
        'vemodel%modulus, vemodel%lambda, vemodel%nonlin'
      stop

    end if

    if ( vemodel%flowtype /= 0 .and. vemodel%flowtype /= 1 ) then

      write(*,'(/a/2a/)') &
        'Error in rhs_plastic_strain_2D_b: incorrect value of flowtype.', &
        'Call create_viscoelastic_model with flowtype=0 (2D) or ', &
        'flowtype=1 (2D, axisymmetric) '
      stop

    end if

    if ( present(mvemodel) ) then
      if ( size(mvemodel,1) /= size(b,1) ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_2D_b:', &
          ' dimension of mvemodel must be identical to size(b,1)'
        stop
      end if
      if ( any( mvemodel(:)%model /= vemodel%model ) .or. &
           any( mvemodel(:)%flowtype /= vemodel%flowtype ) ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_2D_b:', &
          ' model and/or flowtype in mvemodel inconsistent with vemodel'
        stop
      end if
    end if

    if ( present(mode1) .and. present(mode2) ) then
      if ( mode1 < 1 .or. mode1 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_2D_b:', &
          ' mode1 out of range '
        stop
      end if
      if ( mode2 < 1 .or. mode2 > vemodel%nmodes ) then
        write(*,'(/a/a/)') 'Error in rhs_plastic_strain_2D_b:', &
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
          write(*,'(/a/a,i0/)') 'Error in rhs_plastic_strain_2D_b:', &
            ' single mode Jacobian matrix not available for multi-mode', &
            ' von-Mises stress models '
          stop
        endif

        vemmod%compute_drlxmod = .true.

        compute_drhs = .true.

      end if

      if ( vemopt_gammap%compute_drhs_mm ) then

        vemmod%compute_drlxmod = .true.

        allocate( vemmod%drlxmod_mm(1,vemodel%ncompb,vemodel%nmodes) )

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

    allocate(vemmod%drlxmod(1,vemodel%ncompb))
    allocate(vemmod%drlxmoddJ(1))

    if ( vemodel%vmglobal ) then
      allocate ( vemmod%b(size(b,2),size(b,3)) )
      vm1 = vemmod%vmmode1
      vm2 = vemmod%vmmode2
    end if

    mode = vemodel%gammap_mode
    ncompb = vemodel%ncompb
    np = size(b,1)

    do ip = 1, np

!     set gammap

      vemmod%gammap = gammap(ip)

!     set J

      if ( vemmod%dep_J ) vemmod%J = vemopt_gammap%J(ip)

!     von Mises for multiple modes

      if ( vemodel%vmglobal ) vemmod%b(:,vm1:vm2) = b(ip,:,vm1:vm2)

!     fill right-hand side

      if ( present(mvemodel) ) then
        call rhs_plastic_strain_generic ( mvemodel(ip), mode, b(ip,:,mode), &
          vonmises=vm_2D_sm_b, vonmises_mm=vm_2D_mm_b, &
          rhs_gammap=rhs_gammap(ip), vemmod=vemmod )
      else
        call rhs_plastic_strain_generic ( vemodel, mode, b(ip,:,mode), &
          vonmises=vm_2D_sm_b, vonmises_mm=vm_2D_mm_b, &
          rhs_gammap=rhs_gammap(ip), vemmod=vemmod )
      end if

!     fill derivative (single mode) of right-hand side

      if ( compute_drhs ) then
        vemopt_gammap%drhs(ip,:ncompb) = vemmod%drlxmod(1,:)
      end if

!     fill derivative (multi-mode) with respect to conformation tensor

      if ( compute_drhs_mm ) then

        if ( vemodel%vmglobal ) then

!         mode coupling due to global von Mises stress

!         diagonal
          vemopt_gammap%drhs_mm(ip,:ncompb,mode) = vemmod%drlxmod_mm(1,:,mode)

!         off-diagonal,
          do m = vemmod%vmmode1, vemmod%vmmode2
            if ( m == mode ) cycle ! diagonal already accounted for
            vemopt_gammap%drhs_mm(ip,:ncompb,m) = vemmod%drlxmod_mm(1,:,m)
          end do

        else

!         no mode coupling: only mode given by mode

          vemopt_gammap%drhs_mm(ip,:ncompb,modenr1:mode-1) = 0
          vemopt_gammap%drhs_mm(ip,:ncompb,mode) = vemmod%drlxmod(1,:)
          vemopt_gammap%drhs_mm(ip,:ncompb,mode+1:modenr2) = 0

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

  end subroutine rhs_plastic_strain_2D_b

end module viscoelastic_models_2D_b_m

