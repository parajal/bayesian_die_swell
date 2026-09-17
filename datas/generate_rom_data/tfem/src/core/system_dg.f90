
! Copyright (C) 2004-2008 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Routines, for the system of equations

module system_dg_m

  use kind_defs_m
  use mesh_m
  use problem_defs_m
  use sparse_m
  use system_defs_m
  use system_matrix_m
  use element_defs_m

  implicit none


contains


! Create system matrix DG structure (ia not yet accumulated)

  subroutine create_sysmatrix_structure_dg ( sysmatrix, mesh, problem, physq, &
    blocksize, elemsubzeros )

    type(sysmatrix_t), intent(inout) :: sysmatrix
    type(mesh_t),      intent(in)  :: mesh
    type(problem_t),   intent(in)  :: problem

!   if present, physq specifies the physical quantities involved in the DG
!   system matrix part. For example physq=(/2,3)/, means physical quanties 2
!   and 3 are involved and connected. It is possible to call this routines
!   more than once using different physical quantities. Note, that calling this
!   routine two times: once with physq=(/2/) and once with physq=(/3)/ is
!   different from one time with physq=(/2,3)/, because of the assumed
!   connections.
!   If not present all degrees of freedom in the internal nodes are involved
!   and connected.
    integer, intent(in), dimension(:), optional :: physq

!   if present, blocksize specifies that the connections between internal
!   variables is blockwise only. The size of the array blockwise must be equal
!   to the number of element groups.
!   To give an example: with triangles using linear interpolation DG,
!   there are 3 unknowns per element. If the computations
!   involve two components (a vector) all 6 unknowns (3 for each component)
!   will be connected. With respect to the DG part the components
!   might act independently, for example independent convection. In that case,
!   ordering the unknowns like 3+3 and using blocksize=(/3/) will make only
!   connections _within_ each component.
    integer, intent(in), dimension(:), optional :: blocksize

!   If present this is an element subroutine that must be supplied by the
!   calling routine. If not present it is assumed that elemmatzeros=.false.
!   (full matrix).
!   The matrix elemmatzeros has the same shape as elemmat of the element
!   routine elemsub, but is a logical. It gives the entries in the matrix that
!   are zero (.true.) or .non-zero (.false.).
!   The row unknowns are with respect to the current element elgrp/elem and
!   the column unknowns are with respect to the adjacent element at side given
!   by the side heading parameter. The row and column unknowns are affected by
!   the presence of physq.
!   Note that the element matrix elemmatzeros is adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
!   NOTE: elemsubzeros cannot be used together with blocksize
    optional :: elemsubzeros
    interface
      subroutine elemsubzeros ( mesh, problem, elgrp, elem, side, first, &
        last, elemmatzeros )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: elgrp, elem, side
        logical, intent(in) :: first, last
        logical, intent(out), dimension(:,:) :: elemmatzeros
      end subroutine elemsubzeros
    end interface


!   This routine fills the arrays ia with the number of non-zeros in each row:
!
!       ia(i+1) = number of non-zeros in row i
!
!   ia is not yet accumulated and may still be adapted (add more non-zeros)
!   Here the discontinuous Galerkin structure (connections between internal
!   nodes of adjacent elements) is created.

    integer :: numess, numund, elgrp, elem, ndofm, dof_i, side
    integer :: ext_elgrp, ext_elem, dof_e, deg_r, deg_c, row, col, rowp
    integer :: block, deg_c1, deg_c2, nblocks
    integer, dimension(mesh%nelgrp) :: ndof
    integer, dimension(:), allocatable :: pos_i, pos_e
    logical, allocatable, dimension(:,:) :: elemmatzeros
    logical :: first, last


    call check ( mesh, 'create_sysmatrix_structure_dg' )
    call check ( problem, 'create_sysmatrix_structure_dg', mesh )

    if ( sysmatrix%finalized ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysmatrix_structure_dg: ', &
        ' sysmatrix has already been finalized '
      stop
    end if

    if ( sysmatrix%symmetric ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysmatrix_structure_dg: ', &
        ' a symmetric sysmatrix is incompatible with the DG structure'
      stop
    end if

!   do some testing on optional parameters

    if ( present(physq) ) then
      if ( problem%nphysq == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in create_sysmatrix_structure_dg: ',&
          ' specification of physical quantities (physq) ', &
          ' is only available if physical quantities have been defined.'
        stop
      end if
      if ( any( physq < 1 ) .or.  any( physq > problem%nphysq ) ) then
        write(*,'(2(/a)/)') &
          'Error in create_sysmatrix_structure_dg: ', &
          ' physical quantities in physq are out of range.'
        stop
      end if
    end if

    if ( present(blocksize) ) then
      if ( present(elemsubzeros) ) then
        write(*,'(2(/a)/)') &
          'Error in create_sysmatrix_structure_dg: ', &
          ' blocksize cannot be used with elemsubzeros present'
        stop
      end if
      if ( size(blocksize) == mesh%nelgrp ) then
        write(*,'(2(/a)/)') &
          'Error in create_sysmatrix_structure_dg: ', &
          ' size of blocksize must be equal to the number of element groups'
        stop
      end if
      if ( any( blocksize < 1 ) ) then
        write(*,'(2(/a)/)') &
          'Error in create_sysmatrix_structure_dg: ', &
          ' blocksize must be larger than 0'
        stop
      end if
    end if

    numess = problem%numessdegfd
    numund = problem%numundegfd

    if ( .not. sysmatrix%initialized_structure ) then

!     initialize system_matrix

      call initialize_sysmatrix_structure ( sysmatrix, numund, numess )

    end if

!   find maximum number of degrees of freedom in the internal nodes

    do elgrp = 1, mesh%nelgrp
      if ( present(physq) ) then
        ndof(elgrp) = sum ( problem%vec_elnumdegfd(elgrp)&
                            &%a(mesh%element(elgrp)%internnod,physq) )
      else
        ndof(elgrp) = sum ( problem%elnumdegfd(elgrp)&
                            &%a(mesh%element(elgrp)%internnod) )
      end if
    end do
    ndofm = maxval(ndof)

!   reserve memory for positions and elemmatzeros

    allocate ( pos_i(ndofm), pos_e(ndofm) )
    if ( present(elemsubzeros) ) then
      allocate ( elemmatzeros(ndofm,ndofm) )
    end if

!   scan all internal nodes for degrees of freedom and look for connections

    do elgrp = 1, problem%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)

!       get internal degrees of freedom
        if ( present(physq) ) then
          call pos_array_dg ( mesh, problem, elgrp, elem, dof_i, pos_i, physq )
        else
          call pos_array_dg ( mesh, problem, elgrp, elem, dof_i, pos_i )
        end if

!       loop over all sides
        do side = 1, mesh%element(elgrp)%numsides

          ext_elgrp = mesh%sidelem(elgrp)%a(side,elem,1)
          ext_elem  = mesh%sidelem(elgrp)%a(side,elem,2)

          if (  ext_elgrp /= 0 ) then

!           get external degrees of freedom
            if ( present(physq) ) then
              call pos_array_dg ( mesh, problem, ext_elgrp, ext_elem, dof_e, &
                                  pos_e, physq )
            else
              call pos_array_dg ( mesh, problem, ext_elgrp, ext_elem, dof_e, &
                                  pos_e )
            end if

          else

!          no element connected, but element routine is always entered to be
!          consistent and let elemsubzeros possibly deallocate memory

           dof_e = 0

          end if

!         zero elements?

          if ( present(elemsubzeros) ) then

            first = elem == 1 .and. side == 1
            last  = elem == mesh%grpnumel(elgrp) .and. &
                      side == mesh%element(elgrp)%numsides

            call elemsubzeros ( mesh, problem, elgrp, elem, side, first, &
              last, elemmatzeros(1:dof_i,1:dof_e) )

          end if

          if ( ext_elgrp == 0 ) cycle

          if ( present(blocksize) ) then

            nblocks = dof_i / blocksize(elgrp)

            if ( nblocks * blocksize(elgrp) /= dof_i .or. &
                 nblocks * blocksize(ext_elgrp) /= dof_e ) then
              write(*,'(/2a//2(3(a,i0)/))') &
                'Error in create_sysmatrix_structure_dg: ', &
                ' blocksize incompatible with degrees of freedom:', &
                'elgrp = ', elgrp, ', dof_i = ', dof_i, &
                ', blocksize = ', blocksize(elgrp), &
                'ext_elgrp = ', ext_elgrp, ', dof_e = ', dof_e, &
                ', blocksize = ', blocksize(ext_elgrp)
                stop
            end if

          end if

!         loop over internal degrees
          do deg_r = 1, dof_i

            row = pos_i(deg_r)

            if ( present(blocksize) ) then
              block = (deg_r-1) / blocksize(elgrp)
              deg_c1 = block * blocksize(ext_elgrp) + 1
              deg_c2 = min(deg_c1+blocksize(ext_elgrp)-1,dof_e)
            else
              deg_c1 = 1
              deg_c2 = dof_e
            end if

            do deg_c = deg_c1, deg_c2

              col = pos_e(deg_c)

              if ( present(elemsubzeros) ) then
                if ( elemmatzeros(deg_r,deg_c) ) cycle
              end if

              if ( row <= numund ) then

!               unknown degree of freedom (row)

                if ( col <= numund ) then

!                 unknown degree of freedom (column)
                  sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1

                else if ( col > numund ) then

!                 prescribed degree of freedom (column)
                  sysmatrix%Sup%ia(row+1) = sysmatrix%Sup%ia(row+1) + 1

                end if

              else if ( row > numund ) then

!               prescribed degree of freedom (row)

                rowp = row - numund

                if ( col <= numund ) then

!                 unknown degree of freedom (column)
                  sysmatrix%Spu%ia(rowp+1) = sysmatrix%Spu%ia(rowp+1) + 1

                else if ( col > numund ) then

!                 prescribed degree of freedom (column)
                  sysmatrix%Spp%ia(rowp+1) = sysmatrix%Spp%ia(rowp+1) + 1

                end if

              end if

            end do

          end do

        end do

      end do
    end do

    deallocate ( pos_i, pos_e )

  end subroutine create_sysmatrix_structure_dg


! Position array DG

  subroutine pos_array_dg ( mesh, problem, elgrp, elem, dof, pos, physqarr )

    type(mesh_t), intent(in)  :: mesh
    type(problem_t), intent(in)  :: problem
    integer, intent(in) :: elgrp, elem
    integer, intent(out) :: dof
    integer, dimension(:), intent(out) :: pos
    integer, dimension(:), intent(in), optional :: physqarr

!   this routine computes the positions of the degrees of freedom of the
!   internal degrees of freedom. The ordering is such that the inner loop is
!   always the nodes (relevant only if there is more than 1) and the outer
!   loop are the physical quantities (if physqarr is present).

    integer :: nodenr, node, physq, bp, nndof, deg, maxdeg, intnode

    if ( present(physqarr) ) then

      dof = 0

      do physq = 1, size(physqarr)

        maxdeg = maxval( problem%vec_elnumdegfd(elgrp)&
          &%a(mesh%element(elgrp)%internnod,problem%physq(physqarr(physq))) )

        do deg = 1, maxdeg

          do node = 1, mesh%element(elgrp)%numinternnod

            intnode = mesh%element(elgrp)%internnod(node)
            nodenr = mesh%topology(elgrp)%a(intnode,elem)

!           positions

            bp = problem%nodnumdegfd(nodenr) + &
                   sum( problem%vec_elnumdegfd(elgrp)%a(intnode,&
                           &problem%physq(1:physqarr(physq)-1) ) )
            nndof = problem%vec_elnumdegfd(elgrp)&
                         &%a(intnode,problem%physq(physqarr(physq)))

!           positions in sysvector and sysmatrix (renumbered)
            if ( deg <= nndof ) then
              pos(dof+1) = problem%degfdperm(bp+deg,2)
              dof = dof + 1
            end if

          end do

        end do

      end do

    else

!     all degrees of freedom

      dof = 0

      maxdeg = maxval(problem%elnumdegfd(elgrp)&
                       &%a(mesh%element(elgrp)%internnod) )

      do deg = 1, maxdeg

        do node = 1, mesh%element(elgrp)%numinternnod

          intnode = mesh%element(elgrp)%internnod(node)
          nodenr = mesh%topology(elgrp)%a(intnode,elem)

!         positions

          bp = problem%nodnumdegfd(nodenr)
          nndof = problem%elnumdegfd(elgrp)%a(intnode)

!         positions in sysvector and sysmatrix (renumbered)
          if ( deg <= nndof ) then
            pos(dof+1) = problem%degfdperm(bp+deg,2)
            dof = dof + 1
          end if

        end do

      end do

    end if

  end subroutine pos_array_dg


! Assemble system matrix for DG

  subroutine build_system_dg ( mesh, problem, sysmatrix, elemsub, &
    coefficients, mcoefficients, oldvectors, physq, blocksize, elemsubzeros, &
    elgroup1, elgroup2, groups, addmatrix )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the system to be assembled
    type(sysmatrix_t), intent(inout) :: sysmatrix

!   this is the element subroutine that must be supplied by the calling routine
!   the row unknowns are with respect to the current element elgrp/elem and
!   the column unknowns are with respect to the adjacent element at side given
!   by the side heading parameter. The row and column unknowns are affected by
!   the presence of physq but _not_ by the presence of blocksize.
!   Note that the element matrix elemmat is adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
    interface
      subroutine elemsub ( mesh, problem, elgrp, elem, side, first, last, &
        coefficients, oldvectors, elemmat )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: elgrp, elem, side
        logical, intent(in) :: first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat
      end subroutine elemsub
    end interface

!   if coefficients is present, it is passed through to the element subroutine
    type(coefficients_t), intent(in), optional :: coefficients

!   array of coefficients (length=number of element groups)
!   If present only one of them is passed through to the element subroutine
!   based on the element group number
    type(coefficients_t), dimension(:), intent(in), optional :: mcoefficients

!   if oldvectors is present, it is passed through to the element subroutine.
!   This structure can be used to supply old sysvectors/vectors from possibly
!   other problems to the element subroutines. Note that the problem structures
!   of the other problems must be included in oldvectors.
    type(oldvectors_t), intent(in), optional :: oldvectors

!   if present, physq specifies the physical quantities involved in the DG
!   system matrix part. For example physq=(/2,3)/, means physical quanties 2
!   and 3 are involved and connected. It is possible to call this routines
!   more than once using different physical quantities. Note, that calling this
!   routine two times: once with physq=(/2/) and once with physq=(/3)/ is
!   different from one time with physq=(/2,3)/, because of the assumed
!   connections.
!   Note, that the element matrix is adjusted to the physical quantities
!   involved.
!   If not present all degrees of freedom in the internal nodes are involved
!   and connected.
    integer, intent(in), dimension(:), optional :: physq

!   if present, blocksize specifies that the connections between internal
!   variables is blockwise only. The size of the array blockwise must be equal
!   to the number of element groups.
!   To give an example: with triangles using linear interpolation DG,
!   there are 3 unknowns per element. If the computations
!   involve two components (a vector) all 6 unknowns (3 for each component)
!   will be connected. With respect to the DG part the components
!   might act independently, for example independent convection. In that case,
!   ordering the unknowns like 3+3 and using blocksize=(/3/) will make only
!   connections _within_ each component.
!   Note, that the element matrix is _not_ adjusted by blocksize. In the
!   example, the shape is (6,6), although only the two (3,3) blocks on the
!   diagonal will be used in the assembly.
!   NOTE: blocksize cannot be used together with elemsubzeros
    integer, intent(in), dimension(:), optional :: blocksize

!   If present this is an element subroutine that must be supplied by the
!   calling routine. If not present it is assumed that elemmatzeros=.false.
!   (full matrix).
!   The matrix elemmatzeros has the same shape as elemmat of the element
!   routine elemsub, but is a logical. It gives the entries in the matrix that
!   are zero (.true.) or .non-zero (.false.).
!   The row unknowns are with respect to the current element elgrp/elem and
!   the column unknowns are with respect to the adjacent element at side given
!   by the side heading parameter. The row and column unknowns are affected by
!   the presence of physq.
!   Note that the element matrix elemmatzeros is adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
!   NOTE: elemsubzeros cannot be used together with blocksize
    optional :: elemsubzeros
    interface
      subroutine elemsubzeros ( mesh, problem, elgrp, elem, side, first, &
        last, elemmatzeros )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: elgrp, elem, side
        logical, intent(in) :: first, last
        logical, intent(out), dimension(:,:) :: elemmatzeros
      end subroutine elemsubzeros
    end interface

!   if these are present the assembling takes place for element groups
!   elgroup1,...,elgroup2 only. If only elgroup1 is present one group is
!   assembled only. Note, that adjacent elements can still be in another group,
!   so that the usefullness of these parameters can be questioned.
    integer, intent(in), optional :: elgroup1, elgroup2

!   if present: the element groups to be assembled
!   For example groups=(/2,4/) will assemble for element groups 2 and 4.
!   Default: all groups
    integer, dimension(:), intent(in), optional :: groups

!   if addmatrix is set to .true. the matrix is not cleared before
!   the assembling and thus the element matrices are added to an
!   existing system. Note, that this parameter is required if this is not the
!   first build stage, such when a call of build_system for a standard FEM
!   build of the matrix has already been preformed.
!   The default is .false. (clearing)
    logical, intent(in), optional :: addmatrix


!   This routine performs the assembly proces and nothing more. This means that
!   the routine can be called as many times as needed with addmatrix=.true.
!   Other parameters, such as material parameters, function numbers,
!   old vectors etc., can be supplied to the element routine by either
!   modules variables (global module variables) or via the heading parameters
!   coefficients and oldvectors.


    logical :: addmat, first, last
    integer :: numess, numund, elgrp, elem, ndofm, dof_i, side, grp
    integer :: nblocks, block
    integer :: ext_elgrp, ext_elem, dof_e
    integer :: deg_r1, deg_r2, deg_c1, deg_c2
    integer :: lgroups(mesh%nelgrp), lnelgrp, i
    integer, dimension(mesh%nelgrp) :: ndof
    integer, dimension(:), allocatable :: pos_i, pos_e, w1, w2
    real(dp), allocatable, dimension(:,:) :: elemmat
    logical, allocatable, dimension(:,:) :: elemmatzeros
    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl

    call check ( mesh, 'build_system_dg' )
    call check ( problem, 'build_system_dg', mesh )

!   test inactive groups

    if ( problem%numinactivegroups > 0 ) then
      write(*,'(2(/a)/)') &
        'Error in build_system_dg: ',&
        ' inactive groups not allowed in the DG system build '
      stop
    end if

!   check physq and blocksize

    if ( present(physq) ) then
      if ( problem%nphysq == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in build_system_dg: ',&
          ' specification of physical quantities (physq) ', &
          ' is only available if physical quantities have been defined.'
        stop
      end if
      if ( any( physq < 1 ) .or.  any( physq > problem%nphysq ) ) then
        write(*,'(2(/a)/)') &
          'Error in build_system_dg: ', &
          ' physical quantities in physq are out of range.'
        stop
      end if
    end if

    if ( present(blocksize) ) then
      if ( present(elemsubzeros) ) then
        write(*,'(2(/a)/)') &
          'Error in build_system_dg: ', &
          ' blocksize cannot be used with elemsubzeros present'
        stop
      end if
      if ( size(blocksize) /= mesh%nelgrp ) then
        write(*,'(2(/a)/)') &
          'Error in build_system_dg: ', &
          ' size of blocksize must be equal to the number of element groups'
        stop
      end if
      if ( any( blocksize < 1 ) ) then
        write(*,'(2(/a)/)') &
          'Error build_system_dg: ', &
          ' blocksize must be larger than 0'
        stop
      end if
    end if

    if ( present(mcoefficients) ) then
      if ( size(mcoefficients) /= mesh%nelgrp ) then
        write(*,'(/2(a/))') &
          'Error in build_system_dg: ', &
          '   size of mcoefficients /= number of element groups'
        stop
      end if
    end if

!   initialization

    if ( present(oldvectors) ) oldvl = oldvectors
    if ( present(coefficients) ) coeffl = coefficients

!   which element groups?

    if ( present(elgroup1) .and. present(elgroup2) ) then
!     specified range of groups only
      lnelgrp = elgroup2 - elgroup1 + 1
      lgroups(1:lnelgrp) = [ (i,i=elgroup1,elgroup2) ]
    else if ( present(elgroup1) ) then
!     one group only
      lnelgrp = 1
      lgroups(1) = elgroup1
    else if ( present(groups) ) then
      lnelgrp = size(groups)
      lgroups(1:lnelgrp) = groups
    else
!     all groups
      lnelgrp = mesh%nelgrp
      lgroups = [ (i,i=1,lnelgrp) ]
    end if

    if ( any ( lgroups(1:lnelgrp) < 1 ) .or. &
         any ( lgroups(1:lnelgrp) > mesh%nelgrp ) ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error: elgroup1, elgroup2, or groups in the heading of ', &
        ' build_system_dg is out of range: ', &
        ' some groups are < 1 or larger than the number of groups ', mesh%nelgrp
      stop
    end if

    do grp = 1, lnelgrp

      elgrp = lgroups(grp)

      if ( all( problem%elnumdegfd(elgrp)&
                   %a(mesh%element(elgrp)%internnod) == 0 ) ) then

        write(*,'(/a/a,i0/)') &
          'Error: element group in the heading of build_system_dg has ', &
          'no degrees of freedom. elgrp is ', elgrp
        stop

      end if

    end do

!   test presence of matrix?

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(/a/)') &
        'Error in build_system_dg: no system matrix structure.'
      stop
    end if

    if ( .not. sysmatrix%finalized ) then
      write(*,'(/a/)') &
        'Error in build_system_dg: system matrix has not been finalized.'
      stop
    end if

    if ( .not. sysmatrix%allocated_data ) then
      write(*,'(/a/)') &
        'Error in build_system_dg: data in system matrix not allocated.'
      stop
    end if

!   clear matrix?

    if ( present(addmatrix) ) then
      addmat = addmatrix
    else
      addmat = .false.
    end if

    if ( .not. addmat ) then

!     clear matrix
      call clear_sysmatrix ( sysmatrix )

    end if

!   temporary work arrays to store information on current large matrix row

    numess = problem%numessdegfd
    numund = problem%numundegfd

    allocate( w1(numund), w2(numess) )

    w1 = 0
    w2 = 0

!   find maximum number of degrees of freedom in the internal nodes

    do elgrp = 1, mesh%nelgrp
      if ( present(physq) ) then
        ndof(elgrp) = sum ( problem%vec_elnumdegfd(elgrp)&
                            &%a(mesh%element(elgrp)%internnod,physq) )
      else
        ndof(elgrp) = sum ( problem%elnumdegfd(elgrp)&
                            &%a(mesh%element(elgrp)%internnod) )
      end if
    end do
    ndofm = maxval(ndof)

!   reserve memory for positions and element matrix

    allocate ( pos_i(ndofm), pos_e(ndofm) )
    allocate ( elemmat(ndofm,ndofm) )
    if ( present(elemsubzeros) ) then
      allocate ( elemmatzeros(ndofm,ndofm) )
    end if

!   start loop over all elements

    do grp = 1, lnelgrp

      elgrp = lgroups(grp)

      if ( present(mcoefficients) ) coeffl = mcoefficients(elgrp)

      do elem = 1, mesh%grpnumel(elgrp)

!       get internal degrees of freedom
        if ( present(physq) ) then
          call pos_array_dg ( mesh, problem, elgrp, elem, dof_i, pos_i, physq )
        else
          call pos_array_dg ( mesh, problem, elgrp, elem, dof_i, pos_i )
        end if

!       loop over all sides
        do side = 1, mesh%element(elgrp)%numsides

          ext_elgrp = mesh%sidelem(elgrp)%a(side,elem,1)
          ext_elem  = mesh%sidelem(elgrp)%a(side,elem,2)

          if ( ext_elgrp /= 0 ) then

!           get external degrees of freedom
            if ( present(physq) ) then
              call pos_array_dg ( mesh, problem, ext_elgrp, ext_elem, dof_e, &
                                  pos_e, physq )
            else
              call pos_array_dg ( mesh, problem, ext_elgrp, ext_elem, dof_e, &
                                  pos_e )
            end if

          else

!          no element connected, but element routine is always entered to be
!          consistent and let elemsub possibly deallocate memory

           dof_e = 0

          end if

          first = elem == 1 .and. side == 1
          last  = elem == mesh%grpnumel(elgrp) .and. &
                    side == mesh%element(elgrp)%numsides

!         compute element matrix

          call elemsub ( mesh, problem, elgrp, elem, side, first, last, &
            coeffl, oldvl, elemmat(1:dof_i,1:dof_e) )

!         zero elements?

          if ( present(elemsubzeros) ) then
            call elemsubzeros ( mesh, problem, elgrp, elem, side, first, &
              last, elemmatzeros(1:dof_i,1:dof_e) )
          end if

          if ( ext_elgrp == 0 ) cycle

!         add element matrix row by row

          if ( present(blocksize) ) then

!           in blocks

            nblocks = dof_i / blocksize(elgrp)

            if ( nblocks * blocksize(elgrp) /= dof_i .or. &
                 nblocks * blocksize(ext_elgrp) /= dof_e ) then
              write(*,'(/2a//2(3(a,i0)/))') &
                'Error in build_system_dg: blocksize incompatible ', &
                'with element matrix:', &
                'elgrp = ', elgrp, ', dof_i = ', dof_i, &
                ', blocksize = ', blocksize(elgrp), &
                'ext_elgrp = ', ext_elgrp, ', dof_e = ', dof_e, &
                ', blocksize = ', blocksize(ext_elgrp)
                stop
            end if

            do block = 0, nblocks - 1
              deg_r1 = block * blocksize(elgrp) + 1
              deg_r2 = deg_r1+blocksize(elgrp)-1
              deg_c1 = block * blocksize(ext_elgrp) + 1
              deg_c2 = deg_c1+blocksize(ext_elgrp)-1
              call add_elemmat_to_sysmatrix &
                ( sysmatrix, elemmat(deg_r1:deg_r2,deg_c1:deg_c2), &
                  pos_i(deg_r1:deg_r2), pos_e(deg_c1:deg_c2), w1, w2 )
            end do

          else

!           full matrix

            if ( present(elemsubzeros) ) then
              call add_elemmat_to_sysmatrix &
                ( sysmatrix, elemmat(1:dof_i,1:dof_e), pos_i(1:dof_i), &
                  pos_e(1:dof_e), w1, w2, elemmatzeros(1:dof_i,1:dof_e) )
            else
              call add_elemmat_to_sysmatrix &
                ( sysmatrix, elemmat(1:dof_i,1:dof_e), pos_i(1:dof_i), &
                  pos_e(1:dof_e), w1, w2 )
            end if

          end if

        end do

      end do
    end do

    if ( present(elemsubzeros) ) then
      deallocate ( elemmatzeros )
    end if
    deallocate ( elemmat )
    deallocate ( pos_i, pos_e )
    deallocate ( w1, w2 )

  end subroutine build_system_dg

end module system_dg_m
