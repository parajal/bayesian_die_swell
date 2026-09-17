
! Copyright (C) 2004-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Routines, for the system of equations (matrix)

module system_matrix_m

  use glob_defs_m
  use mesh_m
  use problem_defs_m
  use pos_array_m
  use sparse_m
  use system_defs_m
  use system_vector_m
  use element_defs_m

  implicit none

! type for integer array using a pointer

  type int_array_1d_p
    integer, pointer, dimension(:) :: a => null()
  end type int_array_1d_p

! interface for generic delete subroutine

  interface delete
    module procedure delete_sysmatrix
  end interface delete

! interface for generic check subroutine

  interface check
    module procedure check_filled_sysmatrix
  end interface check

! interface for generic create subroutine

  interface create
    module procedure create_sysmatrix
  end interface create

! interface for generic copy subroutine

  interface copy
    module procedure copy_sysmatrix
  end interface copy

! interface for generic clear_rows subroutine

  interface clear_rows
    module procedure clear_rows_sysmatrix
  end interface clear_rows

contains


! Create system matrix basic structure.

  subroutine create_sysmatrix_structure ( sysmatrix, mesh, problem, &
    symmetric, usephysqmask )

    type(sysmatrix_t), intent(inout) :: sysmatrix
    type(mesh_t),      intent(in)  :: mesh
    type(problem_t),   intent(in)  :: problem

!   if symmetric is present and has the values .true. the matrix is assumed to
!   be symmetric and only the upper-triangle (including the diagonal) of the
!   matrix of unknowns (Suu) is stored.
!   NOTE: this has no effect if sysmatrix has been initialized already.
!   default: symmetric=.false.
    logical, intent(in), optional :: symmetric

!   if present, do not reserve space for the partitions denoted by .false.
!   in the matrix problem%physqmask.
!   default=.false.
    logical, intent(in), optional :: usephysqmask

!   only initializes the basic row structure (ia) and the size (n,m,nnz) of all
!   submatrices of the system matrix.


!   create system matrix basic structure (ia not yet accumulated)

    call check ( mesh, 'create_sysmatrix_structure' )
    call check ( problem, 'create_sysmatrix_structure', mesh )

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric, usephysqmask )

!   accumulate the ia vector and compute nnz

    call finalize_sysmatrix_structure ( sysmatrix )

  end subroutine create_sysmatrix_structure


! Finalize the structure of the sysmatrix, i.e. accumulate the ia vector
! and compute nnz

  subroutine finalize_sysmatrix_structure ( sysmatrix )

    type(sysmatrix_t), intent(inout) :: sysmatrix

!   test

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(2(/a)/)') &
        'Error in finalize_sysmatrix_structure: ', &
        ' sysmatrix structure has not been created '
      stop
    end if

    if ( sysmatrix%finalized ) then
      write(*,'(2(/a)/)') &
        'Error in finalize_sysmatrix_structure: ', &
        ' sysmatrix has already been finalized '
      stop
    end if

!   accumulate the ia vector and compute nnz

    call accumulate_ia_nnz ( sysmatrix%Suu )
    call accumulate_ia_nnz ( sysmatrix%Sup )
    call accumulate_ia_nnz ( sysmatrix%Spu )
    call accumulate_ia_nnz ( sysmatrix%Spp )

    sysmatrix%finalized = .true.

  contains

!   helper routine

    subroutine accumulate_ia_nnz ( matrix )

      type(sparsematrix_t), intent(inout) :: matrix

      integer :: i

      if ( matrix%n == 0 .or. matrix%m == 0 ) return

      matrix%ia(1) = 1
      do i = 1, matrix%n
        matrix%ia(i+1) = matrix%ia(i+1) + matrix%ia(i)
      end do
      matrix%nnz = matrix%ia(matrix%n+1) - 1

    end subroutine accumulate_ia_nnz

  end subroutine finalize_sysmatrix_structure


! Create system matrix basic structure (ia not yet accumulated)

  subroutine create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric, usephysqmask )

    type(sysmatrix_t), intent(inout) :: sysmatrix
    type(mesh_t),      intent(in)  :: mesh
    type(problem_t),   intent(in)  :: problem

!   if symmetric is present and has the values .true. the matrix is assumed to
!   be symmetric and only the non-zero elements of the upper-triangle
!   (including the diagonal) of the matrix of unknowns (Suu) are stored.
!   NOTE: this has no effect if sysmatrix has been initialized already.
!   default: symmetric=.false.
    logical, intent(in), optional :: symmetric

!   if present, do not reserve space for the partitions denoted by .false.
!   in the matrix problem%physqmask.
!   default=.false.
    logical, intent(in), optional :: usephysqmask

!   This routine fills the arrays ia with the number of non-zeros in each row:
!
!       ia(i+1) = number of non-zeros in row i
!
!   ia is not yet accumulated and may still be adapted (add more non-zeros)
!   Here the basic FEM structure is created.

    integer :: numund, numess, nodenr, deg_r1, deg_r2, n1, n2, numnod, deg_r
    integer :: row, i, node, deg_c1, deg_c2, deg_c, rowp, col
    integer, dimension(mesh%maxnodnumnod+1) :: nodes
    integer, allocatable, dimension(:) :: utype
!   array of length numnodaldegfd specifying the physical quantity of each
!   nodal degree of freedom
    integer, allocatable, dimension(:) :: physq_of_degfd
    logical :: lusephysqmask

    allocate ( utype(problem%numdegfd) )

    call check ( mesh, 'create_sysmatrix_structure_base' )
    call check ( problem, 'create_sysmatrix_structure_base', mesh )

    if ( sysmatrix%finalized ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysmatrix_structure_base: ', &
        ' sysmatrix has already been finalized '
      stop
    end if

    lusephysqmask = set_optional ( variable=usephysqmask, default=.false. )

    if ( lusephysqmask .and. problem%nphysq == 0 ) then
      write(*,'(/a/a/a/)') &
        'Error in create_sysmatrix_structure_base: ', &
        '  masking physical quantity blocks ', &
        '  cannot be used since physical quantities have not been defined.'
      stop
    end if

!   fill temporary utype array that gives the type of the degree of freedom
!      1: unknown  3: prescribed

    numess = problem%numessdegfd
    numund = problem%numundegfd

    where ( problem%degfdperm(:,2) > numund )
      utype = 3
    else where
      utype = 1
    end where

    if ( .not. sysmatrix%initialized_structure ) then

!     initialize system_matrix

      call initialize_sysmatrix_structure ( sysmatrix, numund, numess )

!     set symmetry of the matrix

      sysmatrix%symmetric = set_optional ( variable=symmetric, default=.false. )

    end if

    if ( lusephysqmask ) then

!     take problem%physqmask array into account to discard zero blocks

      if ( sysmatrix%symmetric .and. &
         any( transpose(problem%physqmask) .neqv. problem%physqmask ) ) then
        write(*,'(2(/a)/)') &
          'Error in create_sysmatrix_structure_base: ', &
          ' a symmetric sysmatrix requires problem%usephysqmask to be symmetric'
        stop
      end if

      allocate ( physq_of_degfd(problem%numnodaldegfd) )

      call fill_physq_of_degfd

    end if

!   scan all nodes for degrees of freedom and look for connections

    do nodenr = 1, mesh%nnodes

      deg_r1 = problem%nodnumdegfd ( nodenr ) + 1
      deg_r2 = problem%nodnumdegfd ( nodenr + 1 )

!     nodes are current node + all connections

      nodes(1) = nodenr
      n1 = mesh%nodnumnod(nodenr)
      n2 = mesh%nodnumnod(nodenr+1)
      numnod = n2 - n1
      nodes(2:numnod+1) = mesh%nodnod(n1+1:n2)

      do deg_r = deg_r1, deg_r2

        row = problem%degfdperm(deg_r,2)

        do i = 1, numnod + 1

          node = nodes(i)

          deg_c1 = problem%nodnumdegfd ( node ) + 1
          deg_c2 = problem%nodnumdegfd ( node + 1 )

          do deg_c = deg_c1, deg_c2

            if ( lusephysqmask ) then
!             in mask block?
              if ( .not. problem%physqmask(physq_of_degfd(deg_r), &
                                           physq_of_degfd(deg_c) ) ) &
                cycle
            end if

            if ( utype(deg_r) == 1 ) then

!             unknown degree of freedom (row)

              if ( utype(deg_c) == 1 ) then

!               unknown degree of freedom (column)

                if ( sysmatrix%symmetric ) then
!                 symmetric matrix
                  col = problem%degfdperm(deg_c,2)
                  if ( row > col ) cycle  ! element in lower triangle
                end if

                sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1

              else if ( utype(deg_c) == 3 ) then

!               prescribed degree of freedom (column)
                sysmatrix%Sup%ia(row+1) = sysmatrix%Sup%ia(row+1) + 1

              end if

            else if ( utype(deg_r) == 3 ) then

!             prescribed degree of freedom (row)

              rowp = row - numund

              if ( utype(deg_c) == 1 ) then

!               unknown degree of freedom (column)
                sysmatrix%Spu%ia(rowp+1) = sysmatrix%Spu%ia(rowp+1) + 1

              else if ( utype(deg_c) == 3 ) then

!               prescribed degree of freedom (column)
                sysmatrix%Spp%ia(rowp+1) = sysmatrix%Spp%ia(rowp+1) + 1

              end if

            end if

          end do

        end do

      end do

    end do

    deallocate ( utype )
    if ( lusephysqmask ) deallocate ( physq_of_degfd )

  contains


!   fill physq_of_degfd array: physical quantity for each nodal degree

    subroutine fill_physq_of_degfd

      integer :: dof, physq, nndof, nodenr

      dof = 0

!     loop over nodes

      do nodenr = 1, mesh%nnodes

!       positions in the node

        do physq = 1, problem%nphysq

          nndof = problem%vec_nodnumdegfd(nodenr+1,problem%physq(physq)) &
                   - problem%vec_nodnumdegfd(nodenr,problem%physq(physq))

!         fill local positions in node with physq
          physq_of_degfd(dof+1:dof+nndof) = physq

          dof = dof + nndof

        end do

      end do

      if ( dof /= problem%numnodaldegfd ) stop ' stop internal error dof '

    end subroutine fill_physq_of_degfd

  end subroutine create_sysmatrix_structure_base


! Allocate and initialize sysmatrix structure (the array ia)

  subroutine initialize_sysmatrix_structure ( sysmatrix, numund, numess )

    type(sysmatrix_t), intent(inout) :: sysmatrix
    integer, intent(in) :: numund, numess

!   This is a helper routine. Not intended for direct user call.

    sysmatrix%Suu%n = numund
    sysmatrix%Suu%m = numund
    sysmatrix%Suu%nnz = 0
    allocate(sysmatrix%Suu%ia(numund+1))
    sysmatrix%Suu%ia = 0

    sysmatrix%Sup%n = numund
    sysmatrix%Sup%m = numess
    sysmatrix%Sup%nnz = 0
    allocate(sysmatrix%Sup%ia(numund+1))
    sysmatrix%Sup%ia = 0

    sysmatrix%Spu%n = numess
    sysmatrix%Spu%m = numund
    sysmatrix%Spu%nnz = 0
    allocate(sysmatrix%Spu%ia(numess+1))
    sysmatrix%Spu%ia = 0

    sysmatrix%Spp%n = numess
    sysmatrix%Spp%m = numess
    sysmatrix%Spp%nnz = 0
    allocate(sysmatrix%Spp%ia(numess+1))
    sysmatrix%Spp%ia = 0

    sysmatrix%initialized_structure = .true.

  end subroutine initialize_sysmatrix_structure


! Allocate sysmatrix (the matrices data (a) + column numbers (ja))

  subroutine create_sysmatrix_data ( sysmatrix )

    type(sysmatrix_t), intent(inout) :: sysmatrix

!   test

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysmatrix_data: ', &
        ' sysmatrix structure has not been created '
      stop
    end if

    if ( .not. sysmatrix%finalized ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysmatrix_data: ', &
        ' sysmatrix has not been finalized '
      stop
    end if

    if ( sysmatrix%allocated_data ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysmatrix_data: ', &
        ' sysmatrix data has already been allocated '
      stop
    end if

    allocate( sysmatrix%Suu%a(sysmatrix%Suu%nnz) )
    allocate( sysmatrix%Suu%ja(sysmatrix%Suu%nnz) )
    allocate( sysmatrix%Sup%a(sysmatrix%Sup%nnz) )
    allocate( sysmatrix%Sup%ja(sysmatrix%Sup%nnz) )
    allocate( sysmatrix%Spu%a(sysmatrix%Spu%nnz) )
    allocate( sysmatrix%Spu%ja(sysmatrix%Spu%nnz) )
    allocate( sysmatrix%Spp%a(sysmatrix%Spp%nnz) )
    allocate( sysmatrix%Spp%ja(sysmatrix%Spp%nnz) )

    allocate( sysmatrix%permu(0), sysmatrix%ipermu(0) )

    sysmatrix%allocated_data = .true.
    sysmatrix%renumber = .false.

  end subroutine create_sysmatrix_data


! Create a single new system matrix from an existing one.
! NOTE: the new sysmatrix contains the structure and all data has been
! allocated. However the data has NOT been copied from the existing one.
! Use copy_sysmatrix for that.

  subroutine create_single_sysmatrix ( sysmatrix1, sysmatrix2, nr )

    type(sysmatrix_t), intent(inout) :: sysmatrix1, sysmatrix2
    integer, intent(in) :: nr

    if ( .not. sysmatrix1%initialized_structure ) then
      write(*,'(/2a/)') &
        'Error in create_sysmatrix: no system matrix structure in ', &
        'existing sysmatrix'
      stop
    end if

    if ( .not. sysmatrix1%finalized ) then
      write(*,'(/2a/)') &
        'Error in create_sysmatrix: system matrix has not been finalized in ', &
        'existing sysmatrix'
      stop
    end if

    if ( .not. sysmatrix1%allocated_data ) then
      write(*,'(/2a/)') &
        'Error in create_sysmatrix: data in system matrix not allocated in ', &
        'existing sysmatrix'
      stop
    end if

    if ( sysmatrix2%initialized_structure .or. sysmatrix2%finalized .or. &
         sysmatrix2%allocated_data ) then
      write(*,'(/2a/a,i0/)') &
        'Error in create_sysmatrix: new sysmatrix has already been (partly) ', &
        'created.', &
        ' Matrix number in heading = ', nr
      stop
    end if

!   copy scalars

    sysmatrix2%initialized_structure = sysmatrix1%initialized_structure
    sysmatrix2%finalized = sysmatrix1%finalized
    sysmatrix2%allocated_data = sysmatrix1%allocated_data
    sysmatrix2%symmetric = sysmatrix1%symmetric

    call copy_scalars_sparsematrix ( sysmatrix1%Suu, sysmatrix2%Suu )
    call copy_scalars_sparsematrix ( sysmatrix1%Sup, sysmatrix2%Sup )
    call copy_scalars_sparsematrix ( sysmatrix1%Spu, sysmatrix2%Spu )
    call copy_scalars_sparsematrix ( sysmatrix1%Spp, sysmatrix2%Spp )

!   allocate arrays

    call allocate_sparsematrix ( sysmatrix2%Suu )
    call allocate_sparsematrix ( sysmatrix2%Sup )
    call allocate_sparsematrix ( sysmatrix2%Spu )
    call allocate_sparsematrix ( sysmatrix2%Spp )

!   copy structure

    sysmatrix2%Suu%ia = sysmatrix1%Suu%ia
    sysmatrix2%Sup%ia = sysmatrix1%Sup%ia
    sysmatrix2%Spu%ia = sysmatrix1%Spu%ia
    sysmatrix2%Spp%ia = sysmatrix1%Spp%ia

  contains

!   copy scalars

    subroutine copy_scalars_sparsematrix ( matrix1, matrix2 )

      type(sparsematrix_t), intent(inout) :: matrix1, matrix2

      matrix2%n = matrix1%n
      matrix2%m = matrix1%m
      matrix2%nnz = matrix1%nnz

    end subroutine copy_scalars_sparsematrix

!   Allocate matrix

    subroutine allocate_sparsematrix ( matrix )

      type(sparsematrix_t), intent(inout) :: matrix

      allocate( matrix%a(matrix%nnz), matrix%ja(matrix%nnz) )
      allocate( matrix%ia(matrix%n+1) )

    end subroutine allocate_sparsematrix

  end subroutine create_single_sysmatrix


! Create new system matrix from existing one
! NOTE: the new sysmatrix contains the structure and all data has been
! allocated. However the data has NOT been copied from the existing one.
! Use copy_sysmatrix for that.

  subroutine create_sysmatrix ( sysmatrix, sysmatrix1, sysmatrix2, sysmatrix3, &
    sysmatrix4, sysmatrix5 )

    type(sysmatrix_t), intent(inout) :: sysmatrix, sysmatrix1
    type(sysmatrix_t), intent(inout), optional :: sysmatrix2, sysmatrix3, &
      sysmatrix4, sysmatrix5

    call create_single_sysmatrix( sysmatrix, sysmatrix1, 1 )
    if ( present(sysmatrix2) ) &
      call create_single_sysmatrix( sysmatrix, sysmatrix2, 2 )
    if ( present(sysmatrix3) ) &
      call create_single_sysmatrix( sysmatrix, sysmatrix3, 3 )
    if ( present(sysmatrix4) ) &
      call create_single_sysmatrix( sysmatrix, sysmatrix4, 4 )
    if ( present(sysmatrix5) ) &
      call create_single_sysmatrix( sysmatrix, sysmatrix5, 5 )

  end subroutine create_sysmatrix


! Copy system matrix to another system matrix

  subroutine copy_sysmatrix ( sysmatrix1, sysmatrix2 )

    type(sysmatrix_t), intent(inout) :: sysmatrix1, sysmatrix2

!   This routine performs a real data copy of one matrix to the other
!   by the assignment statement (with allocatables)
!
!      sysmatrix2 = sysmatrix1
!

    if ( .not. sysmatrix1%allocated_data ) then
      write(*,'(/a/)') &
        'Error in copy_sysmatrix: data not allocated in sysmatrix1'
      stop
    end if

    sysmatrix2 = sysmatrix1

  end subroutine copy_sysmatrix


! Assemble system matrix and system vector

  subroutine build_system ( mesh, problem, sysmatrix, sysvector, msysvector, &
    m2sysvector, m3sysvector, elemsub, elemsub1, coefficients, mcoefficients, &
    oldvectors, physqrow, physqcol, order, elgroup1, elgroup2, groups, layer, &
    layer_in_all_nodes, exclude_single_layer, buildmatrix, buildvector, &
    addmatvec, addmat, addvec, zeromatvec, factormat, factorvec, object, &
    onobjectnodes, skipelementfunc, elementset, transform, usephysqmask )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the system to be assembled
    type(sysmatrix_t), intent(inout) :: sysmatrix

!   the right-hand side vector to be assembled
    type(sysvector_t), intent(inout), optional :: sysvector

!   multiple right-hand side vectors to be assembled
    type(sysvector_t), dimension(:), intent(inout), optional :: msysvector
    type(sysvector_t), dimension(:,:), intent(inout), optional :: m2sysvector
    type(sysvector_t), dimension(:,:,:), intent(inout), optional :: m3sysvector

!   this is the element subroutine that must be supplied by the calling routine
!   Note that the element matrix elemmat and elemvec are adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
    optional :: elemsub
    interface
      subroutine elemsub ( mesh, problem, elgrp, elem, matrix, vector, &
        first, last, coefficients, oldvectors, elemmat, elemvec )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: elgrp, elem
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat
        real(dp), intent(out), dimension(:) :: elemvec
      end subroutine elemsub
    end interface

!   this is the element subroutine that must be supplied by the calling routine
!   Note that the element matrix elemmat and elemvec are adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
!   The difference with elemsub is the use of eleminfo_t type.
!   Required for assembling on objects.
    optional :: elemsub1
    interface
      subroutine elemsub1 ( mesh, problem, eleminfo, matrix, vector, &
        first, last, coefficients, oldvectors, elemmat, elemvec )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t, eleminfo_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        type(eleminfo_t), intent(in) :: eleminfo
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat
        real(dp), intent(out), dimension(:) :: elemvec
      end subroutine elemsub1
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

!   if these are present the assembling takes place for physical quantities
!   in physqrow and physqcol only. The size of the arrays determine the
!   number of physical quantities involved. The element matrix is defined
!   accordingly, omitting all partitions of physical unknowns not specified.
!   The element matrix can be non-square if physqrow /= physqcol.
!   For example:
!     physqrow=(/2,1/)
!     physqcol=(/3/)
!   means that the element matrix A consists of the two partitions
!      [ A_23 ]
!      [ A_13 ]
!   and the element vector v consists of the two partitions
!      ( v_2 )
!      ( v_1 )
!   The element subroutine must be defined accordingly.
!   The arrays physqrow and physqcol work independently and each can be omitted.
!   Note, that physqrow affects both the element matrix and vector
!   whereas physqcol only affects the (columns of the) element matrix.
!   The arrays physqrow and physqcol affect the local element matrix and vector
!   only and _do not affect_ the global numbering of the unknowns. However the
!   (column) layout of the assembled sparse matrix can be different due to the
!   different sequence of unknowns in the assembling process.
    integer, intent(in), dimension(:), optional :: physqrow, physqcol

!   The parameter order determines the sequence of the degrees of freedom on
!   elementlevel. There are only two possibilities and affect only the two loops
!   over nodal points and degrees of freedom (within a physical quantity if
!   problem%nphysq>0):
!     order = 'ND' : the most inner loop is over the degrees of freedom
!     order = 'DN' : the most inner loop is over the nodal points
!   Specifically we have:   nphysq = 0   nphysq > 0
!         order = 'ND' :       ND           PND
!         order = 'DN' :       DN           PDN
!   NOTE: if nphys >0 the main ordering (outside loop) of the element
!   degrees of freedom are the physical quantities.
!   If problem%numlayers > 0 the degrees are stacked into layers and we have
!                           nphysq = 0   nphysq > 0
!         order = 'ND' :      NLD           PNLD
!         order = 'DN' :      LDN           PLDN
!   See also the userguide for further explanation.
!
!   For example consider an element with two nodes and two physical quanties:
!   one with two degrees (u,v) and one with a single degree p. Then we have:
!     order = 'ND' : u1, v1, u2, v2, p1, p2
!     order = 'DN' : u1, u2, v1, v2, p1, p2
!   The parameter order generates a permutation of the degrees of freedom on
!   element level, making the life of an element programmer easier,
!   but _does not affect_ the global numbering of the unknowns.
!   The (column) layout of the assembled sparse matrix can be different due to
!   the different sequence of unknowns in the assembling process.
!
!   The default of order is:
!     order = 'ND' : if no physical quantities have been defined (nphysq=0)
!     order = 'DN' : if physical quantities have been defined (nphysq>0)
    character(len=*), intent(in), optional :: order

!   if these are present the assembling takes place for element groups
!   elgroup1,...,elgroup2 only. If only elgroup1 is present one group is
!   assembled only.
    integer, intent(in), optional :: elgroup1, elgroup2

!   if present: the element groups to be assembled
!   For example groups=(/2,4/) will assemble for element groups 2 and 4.
!   Default: all groups
    integer, dimension(:), intent(in), optional :: groups

!   if present the assembling takes place within the single layer of degrees of
!   freedom only.
!   If layer_in_all_nodes=.true. (default) only elements where all nodes have
!   degrees in the specified layer are assembled.
!   If layer_in_all_nodes=.false. also the elements where only part of the
!   nodes have degrees in the specified layer are assembled. Be careful here:
!   the number of degrees of freedom in a node not in the layer will be zero
!   and the size of the elementmatrix/elementvector will reflect that.
    integer, intent(in), optional :: layer
    logical, intent(in), optional :: layer_in_all_nodes

!   If exclude_single_layer=.true. elements where all nodes have degrees in
!   a single layer only are excluded from the assembly.
!   If exclude_single_layer=.false. (default) all elements are assembled.
    logical, intent(in), optional :: exclude_single_layer

!   if assigned the value .false. the assembling of the matrix or vector
!   is not done.
    logical, intent(in), optional :: buildmatrix, buildvector

!   if addmatvec is set to .true. the matrix and vector are not cleared before
!   the assembling and thus the element matrices and vectors are added to an
!   existing system. The default is .false. (clearing)
    logical, intent(in), optional :: addmatvec

!   if addmat is set to .true. the matrix is not cleared before
!   the assembling and thus the element matrices and vectors are added to an
!   existing system. The default is .false. (clearing)
!   NOTE: if addmatvec is present, the value of addmat will be ignored.
    logical, intent(in), optional :: addmat

!   if addvec is set to .true. the vector is not cleared before
!   the assembling and thus the element matrices and vectors are added to an
!   existing system. The default is .false. (clearing)
!   NOTE: if addmatvec is present, the value of addvec will be ignored.
    logical, intent(in), optional :: addvec

!   if zeromatvec is set to .true. the element matrices and element vectors
!   are assumed to be zero and elemsub is not called. This parameter is for
!   filling parts of the system matrix with zeros combined with physqrow and
!   physqcol.
    logical, intent(in), optional :: zeromatvec

!   if present the element matrix and/or vector will be multiplied by the
!   specified factor before assembly into the system matrix/vector.
    real(dp), intent(in), optional :: factormat, factorvec

!   if present the assembling takes place on the part of the domain that
!   is occupied by the object.
!   The object must contain integration points (for objectnodes=.false.).
!   This parameter requires elemsub1 to be defined.
    integer, intent(in), optional :: object

!   if onobjectnodes is set to .true. the assembling takes place for the
!   elements intersected by the nodes of an object. The assembling is node for
!   node, i.e. if an element contains more than one node, the elementmatrix is
!   assembled more than once.
!   default=.false.
!   This is a special case for special applications.
    logical, intent(in), optional :: onobjectnodes

!   this element function must be supplied by the calling routine if
!   only part of the elements must be assembled and must be fully
!   under user control.
!   A value of .true. will skip the element given by (elgrp,elem).
!   NOTE: this parameter is basically superseded by the more flexible
!   build on elementset.
    optional :: skipelementfunc
    interface
      function skipelementfunc ( mesh, problem, elgrp, elem, coefficients, &
        oldvectors )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: elgrp, elem
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        logical :: skipelementfunc
      end function skipelementfunc
    end interface

!   if present the assembling takes place only on the elements of the domain
!   that is defined by the elementset.
!   Note, that this argument can be combined with other arguments such as
!   object, layer etc. It is just a restriction on the elements that will be
!   assembled.
    integer, intent(in), optional :: elementset

!   By setting transform=.false. the transformation matrix is not applied
!   to the element matrix and vector. This means, the element degrees of
!   freedom are defined in the local (transformed) system.
!   default = .true.
    logical, intent(in), optional :: transform

!   if present, do not assemble the physical quantity partitions denoted
!   by .false. in the matrix problem%physqmask.
!   default=.false.
    logical, intent(in), optional :: usephysqmask

!   This routine performs the assembly proces and nothing more. This means that
!   after the first call (with addmatvec=.false.), the routine can be called
!   as many times as needed with addmatvec=.true. Combined with the parameters
!   elemsub, physqrow, physqcol, order, elgroup1, elgroup2, buildmatrix,
!   buildvector very flexible structure is created for building the system.
!   Other parameters, such as material parameters, function numbers,
!   old vectors etc., can be supplied to the element routine by either
!   modules variables (global module variables) or via the heading parameters
!   coefficients and oldvectors.


    logical :: matrix, vector, laddmat, laddvec, first, last, zeros, ltransform
    logical :: lonobjectnodes, llayer_in_all_nodes, lexclude_single_layer
    logical :: constant_element_size ! element size constant within a group
    logical, allocatable, dimension(:) :: wkl
    integer :: elem, elgrp, ndofr, ndofc, grp, intp, nelem, ninti, elemo
    integer :: numund, numess, rowg
    integer :: dofr, dofc, nrhsd, rhsd
    integer :: lgroups(mesh%nelgrp), lnelgrp, i, ogroups(mesh%nelgrp)
    integer, allocatable, dimension(:) :: posr, posc, w1, w2, nlw
    integer, allocatable, dimension(:,:) :: wk1, wk2, nwk
    integer, allocatable, dimension(:) :: pqr, pqc
    logical, allocatable, dimension(:,:) :: elemmatzeros
!   work1: array for denoting transformed degrees of freedom
    logical, allocatable, dimension(:) :: work1
    integer :: nwkelem, node, si, sj, j, k
    type(int_array_1d_p), allocatable, dimension(:) :: intps
    real(dp), allocatable, dimension(:,:) :: elemmat
    real(dp), allocatable, dimension(:) :: elemvec
    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl
    type(eleminfo_t) :: eleminfo
    type(logical_array_1d_t), allocatable, dimension(:) :: buildelements
    logical :: lusephysqmask

    allocate ( wkl(mesh%nnodes) )
    allocate ( wk1(mesh%nelgrp,maxval(mesh%grpnumel)) )
    allocate ( wk2(mesh%nelgrp,maxval(mesh%grpnumel)) )
    allocate ( nwk(mesh%nelem,2) )

    call check ( mesh, 'build_system' )
    call check ( problem, 'build_system', mesh )

    if ( present(mcoefficients) ) then
      if ( size(mcoefficients) /= mesh%nelgrp ) then
        write(*,'(/2(a/))') &
          'Error in build_system: ', &
          '   size of mcoefficients /= number of element groups'
        stop
      end if
    end if

    lonobjectnodes = set_optional ( variable=onobjectnodes, default = .false. )

    if ( present(object) ) then

!     check object

      if ( object <=0 .or. object > mesh%nobjects ) then
        write(*,'(/a/a,i0/)') &
          'Error: object in heading of build_system has wrong value:', &
          'object <=0 or object > number objects = ', mesh%nobjects
        stop
      end if

    end if

    if ( present(object) .and. .not. lonobjectnodes ) then

      if ( .not. mesh%objects(object)%intpoints ) then
        write(*,'(/a/a,i0/a/)') &
          'Error build_system: ', &
          ' building system on object ', object, &
          ' however it does not contain integration points.'
        stop
      end if

      if ( any ( mesh%objects(object)%grpelm_int(:,1,:) == 0 ) ) then
!       not all elements connected
        write(*,'(2(/a)/a,i0/)') &
        'Error in build_system: ', &
        ' reference coordinates in object are missing,', &
        ' object = ', object
        stop
      end if

!     fill ogroups

      ogroups = 0

      nelem = mesh%objects(object)%nelem
      ninti = mesh%objects(object)%ninti

      do elemo = 1, nelem
        do intp = 1, ninti
          elgrp = mesh%objects(object)%grpelm_int(intp,1,elemo)
          ogroups(elgrp) = 1
        end do
      end do

    else if ( present(object) .and. lonobjectnodes ) then

      if ( any ( mesh%objects(object)%grpelm(:,1) == 0 ) ) then
!       not all nodes connected
        write(*,'(2(/a)/a,i0/)') &
        'Error in build_system: ', &
        ' reference coordinates in object are missing,', &
        ' object = ', object
        stop
      end if

!     fill ogroups

      ogroups = 0

      do node = 1, mesh%objects(object)%nnodes
        elgrp = mesh%objects(object)%grpelm(node,1)
        ogroups(elgrp) = 1
      end do

    end if

!   initialize local parameters

    nrhsd = 1

    if ( present(zeromatvec) ) then
      zeros = zeromatvec
    else
      if ( .not. present(elemsub) .and. .not. present(elemsub1) ) then
        write(*,'(2(/a)/)') &
          'Error in build_system: ', &
          ' no elemsub or elemsub1 present in heading '
        stop
      end if
      if ( present(elemsub) .and. present(elemsub1) ) then
        write(*,'(2(/a)/)') &
          'Error in build_system: ', &
          ' elemsub and elemsub1 cannot be both present in heading '
        stop
      end if
      zeros = .false.
    end if

    if ( present(oldvectors) ) oldvl = oldvectors
    if ( present(coefficients) ) coeffl = coefficients

!   test order

    if ( present(order) ) then
      if ( all ( order /= [ 'ND', 'DN' ] ) ) then
        write(*,'(2(/a)/)') &
          'Error in build_system: ', &
          ' heading parameter order must be either ''ND'', ''DN''.'
        stop
      end if
    end if

!   check physqrow, physqcol and order

    if ( present(physqrow) .or. present(physqcol) ) then

      if ( problem%nphysq == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in build_system: ',&
          ' partioning of the element matrix using physqrow or physqcol ', &
          ' is only available if physical quantities have been defined.'
        stop
      end if

    end if

    if ( present(physqrow) ) then
      if ( any( physqrow < 1 ) .or.  any( physqrow > problem%nphysq ) ) then
        write(*,'(2(/a)/)') &
          'Error in build_system: ', &
          ' physical quantities in physqrow are out of range.'
        stop
      end if
    end if

    if ( present(physqcol) ) then
      if ( any( physqcol < 1 ) .or. any( physqcol > problem%nphysq ) ) then
        write(*,'(2(/a)/)') &
          'Error in build_system: ', &
          ' physical quantities in physqcol are out of range.'
        stop
      end if
    end if

    lusephysqmask = set_optional ( variable=usephysqmask, default=.false. )

!   layers

    if ( present(layer) ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in build_system: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in build_system: layer out of range.'
        stop
      end if
      wkl = btest(problem%nodlayers,layer-1)
    end if

    llayer_in_all_nodes = &
               set_optional ( variable=layer_in_all_nodes, default = .true. )

    lexclude_single_layer = &
               set_optional ( variable=exclude_single_layer, default = .false. )

    if ( lexclude_single_layer ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in build_system: ',&
          ' exclude_single_layer=.true. ', &
          ' is only available if layers have been defined.'
        stop
      end if
    end if

    if ( present(layer) .and. lexclude_single_layer ) then
      write(*,'(2(/a)/)') &
        'Error in build_system: ', &
        ' exclude_single_layer=.true. is incompatible with layer present'
      stop
    end if

    if ( present(layer) ) then
      constant_element_size = llayer_in_all_nodes
    else
      constant_element_size = problem%numlayers == 0
    end if

    if ( .not. constant_element_size .or. lexclude_single_layer ) then
      allocate ( nlw(mesh%nnodes) )
      nlw = count_layers ( problem%nodlayers, problem%numlayers )
    end if

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
        ' build_system is out of range: ', &
        ' some groups are < 1 or larger than the number of groups ', mesh%nelgrp
      stop
    end if

    do grp = 1, lnelgrp

      elgrp = lgroups(grp)

!     skip inactive groups
      if ( any( elgrp == problem%inactivegroups ) ) cycle

      if ( present(object) ) then
        if ( ogroups(elgrp) == 0 ) cycle  ! group not in object
      end if

      if ( all( problem%elnumdegfd(elgrp)%a == 0 ) ) then

        write(*,'(/a/a,i0/)') &
          'Error: element group in the heading of build_system has ', &
          'no degrees of freedom. elgrp is ', elgrp
        stop

      end if

    end do

!   build on elementset

    if ( present(elementset) ) then

!     check elementset

      if ( elementset <=0 .or. elementset > mesh%nelementsets ) then
        write(*,'(/a/a,i0/)') &
          'Error: elementset in heading of build_system has wrong value:', &
          'elementset <=0 or elementset > number elementsets = ', &
          mesh%nelementsets
        stop
      end if

!     create logical array for elements

      allocate ( buildelements(mesh%nelgrp) )

      do elgrp = 1, mesh%nelgrp

        allocate ( buildelements(elgrp)%a(mesh%grpnumel(elgrp)) )

        buildelements(elgrp)%a = .false.

        buildelements(elgrp)%a(&
           mesh%elementsets(elementset)%elements(elgrp)%a ) = .true.

      end do

    end if

!   compute matrix and vector?

    if ( present(buildmatrix) ) then
      matrix = buildmatrix
    else
      matrix = .true.
    end if

    if ( present(buildvector) ) then
      vector = buildvector
    else
      vector = .true.
    end if

    if ( .not. ( matrix .or. vector ) ) return

!   clear matrix and vector?

    if ( matrix ) then

      if ( .not. sysmatrix%initialized_structure ) then
        write(*,'(/a/)') &
          'Error in build_system: no system matrix structure.'
        stop
      end if

      if ( .not. sysmatrix%finalized ) then
        write(*,'(/a/)') &
          'Error in build_system: system matrix has not been finalized.'
        stop
      end if

      if ( .not. sysmatrix%allocated_data ) then
        write(*,'(/a/)') &
          'Error in build_system: data in system matrix not allocated.'
        stop
      end if

      if ( sysmatrix%Suu%m + sysmatrix%Sup%m /= problem%numdegfd ) then
        write(*,'(/a/a/)') &
          'Error in build_system: number of degrees of freedom of the ', &
          'system matrix is different from the number in problem.'
        stop
      end if

      if ( present(addmatvec) ) then
        laddmat = addmatvec
      else if ( present(addmat) ) then
        laddmat = addmat
      else
        laddmat = .false.
      end if

      if ( .not. laddmat ) then

!       clear matrix
        call clear_sysmatrix ( sysmatrix )

      end if

    end if

    if ( vector ) then

      if ( present(addmatvec) ) then
        laddvec = addmatvec
      else if ( present(addvec) ) then
        laddvec = addvec
      else
        laddvec = .false.
      end if

      if ( present(sysvector) ) then

!       single right-hand side

        if ( .not. sysvector%created ) then
          write(*,'(/a/)') &
            'Error in build_system: sysvector not created.'
          stop
        end if

        if ( sysvector%n /= problem%numdegfd ) then
          write(*,'(/a/a/)') &
            'Error in build_system: number of degrees of freedom of the ', &
            'sysvector vector different from the number in problem.'
          stop
        end if

        if ( .not. laddvec ) then

!         clear vector
          sysvector%u = 0

        end if

      else if ( present(msysvector) ) then

!       multiple right-hand sides

        nrhsd = size(msysvector)

        if ( .not. all(msysvector%created) ) then
          write(*,'(/a/)') &
            'Error in build_system: msysvector not created.'
          stop
        end if

        if ( any( msysvector%n /= problem%numdegfd ) ) then
          write(*,'(/a/a/)') &
            'Error in build_system: number of degrees of freedom of the ', &
            'msysvector vector different from the number in problem.'
          stop
        end if

        if ( .not. laddvec ) then

!         clear vector
          do rhsd = 1, nrhsd
            msysvector(rhsd)%u = 0
          end do

        end if

      else if ( present(m2sysvector) ) then

!       multiple right-hand sides (2d array)

        nrhsd = size(m2sysvector)

        if ( .not. all(m2sysvector%created) ) then
          write(*,'(/a/)') &
            'Error in build_system: m2sysvector not created.'
          stop
        end if

        if ( any( m2sysvector%n /= problem%numdegfd ) ) then
          write(*,'(/a/a/)') &
            'Error in build_system: number of degrees of freedom of the ', &
            'm2sysvector vector different from the number in problem.'
          stop
        end if

        if ( .not. laddvec ) then

!         clear vectors
          do i = 1, size(m2sysvector,1)
            do j = 1, size(m2sysvector,2)
              m2sysvector(i,j)%u = 0
            end do
          end do

        end if

      else if ( present(m3sysvector) ) then

!       multiple right-hand sides (3d array)

        nrhsd = size(m3sysvector)

        if ( .not. all(m3sysvector%created) ) then
          write(*,'(/a/)') &
            'Error in build_system: m3sysvector not created.'
          stop
        end if

        if ( any( m3sysvector%n /= problem%numdegfd ) ) then
          write(*,'(/a/a/)') &
            'Error in build_system: number of degrees of freedom of the ', &
            'm3sysvector vector different from the number in problem.'
          stop
        end if

        if ( .not. laddvec ) then

!         clear vectors
          do i = 1, size(m3sysvector,1)
            do j = 1, size(m3sysvector,2)
              do k = 1, size(m3sysvector,3)
                m3sysvector(i,j,k)%u = 0
              end do
            end do
          end do

        end if

      else

        write(*,'(/2a/)') &
          'Error in build_system: ', &
          ' no sysvector, msysvector, m2sysvector or m3sysvector present'
        stop

      end if

    end if

!   transformations

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in build_system: ', &
        ' transformation matrix has not been build'
      stop
    end if

    ltransform = set_optional ( variable=transform, default=.true. )

    if ( problem%numtransdegfd > 0 .and. ltransform ) then

!     work array to indicate transformed degrees of freedom
      allocate ( work1(problem%numdegfd) )
      work1 = .false.
      work1(problem%degfdtrans) = .true.

    end if

!   temporary work arrays to store information on current large matrix row

    numess = problem%numessdegfd
    numund = problem%numundegfd

    allocate( w1(numund), w2(numess) )

    w1 = 0
    w2 = 0

    if ( present(object) .and. lonobjectnodes ) then

!     special case: loop over nodes of object

      do node = 1, mesh%objects(object)%nnodes

        elgrp = mesh%objects(object)%grpelm(node,1)
        elem  = mesh%objects(object)%grpelm(node,2)

!       skip inactive groups
        if ( any( elgrp == problem%inactivegroups ) ) cycle

!       skip unspecified groups
        if ( .not. any( elgrp == lgroups(1:lnelgrp) ) ) cycle

!       skip elements not in layer
        if ( present(layer) ) then
          if ( .not. element_in_layer() ) cycle
        end if

!       skip elements not in a single layer
        if ( lexclude_single_layer ) then
          if ( element_in_single_layer() ) cycle
        end if

        if ( present(mcoefficients) ) coeffl = mcoefficients(elgrp)

        first = node == 1
        last  = node == mesh%objects(object)%nnodes

!       skip user specified elements
        if ( present(skipelementfunc) ) then
          if ( skipelementfunc ( mesh, problem, elgrp, elem, &
                                 coeffl, oldvl ) ) cycle
        end if

!       only elements on elementset
        if ( present(elementset) ) then
          if ( .not. buildelements(elgrp)%a(elem) ) cycle
        end if

        if ( constant_element_size ) then
          call allocate_elemmatvec
        else
          call allocate_elemmatvec_layers
        end if

        eleminfo%elgrp  = elgrp
        eleminfo%elem   = elem
        eleminfo%object = object
        eleminfo%node   = node

        call compute_assemble_elemmatvec

        deallocate ( elemmat, elemvec, posr, posc )
        if ( lusephysqmask ) deallocate ( elemmatzeros, pqr, pqc )

      end do

    else if ( present(object) ) then

!     loop over elements of object

!     set to zero work arrays for storing counted elements

      wk1 = 0
      wk2 = 0

      do elemo = 1, nelem

!       find mesh elements and integration points

        call find_elements_integrationpoints ( mesh, object, elemo, nwkelem, &
          wk1, wk2, nwk, intps )

        do i = 1, nwkelem

          elgrp = nwk(i,1)
          elem  = nwk(i,2)

!         skip inactive groups
          if ( any( elgrp == problem%inactivegroups ) ) cycle

!         skip unspecified groups
          if ( .not. any( elgrp == lgroups(1:lnelgrp) ) ) cycle

!         skip elements not in layer
          if ( present(layer) ) then
            if ( .not. element_in_layer() ) cycle
          end if

!         skip elements not in a single layer
          if ( lexclude_single_layer ) then
            if ( element_in_single_layer() ) cycle
          end if

          if ( present(mcoefficients) ) coeffl = mcoefficients(elgrp)

          first = i == 1
          last  = i == nwkelem

!         skip user specified elements
          if ( present(skipelementfunc) ) then
            if ( skipelementfunc ( mesh, problem, elgrp, elem, &
                                   coeffl, oldvl ) ) cycle
          end if

!         only elements on elementset
          if ( present(elementset) ) then
            if ( .not. buildelements(elgrp)%a(elem) ) cycle
          end if

          if ( constant_element_size ) then
            call allocate_elemmatvec
          else
            call allocate_elemmatvec_layers
          end if

          eleminfo%elgrp  = elgrp
          eleminfo%elem   = elem
          eleminfo%object = object
          eleminfo%elemo  = elemo
          eleminfo%nintps = size(intps(i)%a)
          eleminfo%intps  => intps(i)%a

          call compute_assemble_elemmatvec

          deallocate ( elemmat, elemvec, posr, posc )
          if ( lusephysqmask ) deallocate ( elemmatzeros, pqr, pqc )

        end do

        call reset_elements_integrationpoints ( nwkelem, wk1, wk2, nwk, intps )

      end do

    else

!     standard loop of groups, elements

!     start loop over groups

      do grp = 1, lnelgrp

        elgrp = lgroups(grp)

!       skip inactive groups
        if ( any( elgrp == problem%inactivegroups ) ) cycle

        if ( present(mcoefficients) ) coeffl = mcoefficients(elgrp)

        if ( constant_element_size ) call allocate_elemmatvec

!       loop over elements in this group

        do elem = 1, mesh%grpnumel(elgrp)

          first = elem == 1
          last  = elem == mesh%grpnumel(elgrp)

!         skip user specified elements
          if ( present(skipelementfunc) ) then
            if ( skipelementfunc ( mesh, problem, elgrp, elem, &
                                   coeffl, oldvl ) ) cycle
!           avoid leaving allocated memory in elements
            first = .true.; last = .true.
          end if

!         skip elements not in layer
          if ( present(layer) ) then
            if ( .not. element_in_layer() ) cycle
!           avoid leaving allocated memory in elements
            first = .true.; last = .true.
          end if

!         skip elements not in a single layer
          if ( lexclude_single_layer ) then
            if ( element_in_single_layer() ) cycle
!           avoid leaving allocated memory in elements
            first = .true.; last = .true.
          end if

!         only elements on elementset
          if ( present(elementset) ) then
            if ( .not. buildelements(elgrp)%a(elem) ) cycle
!           avoid leaving allocated memory in elements
            first = .true.; last = .true.
          end if

          if ( .not. constant_element_size ) call allocate_elemmatvec_layers

          eleminfo%elgrp  = elgrp
          eleminfo%elem   = elem

          call compute_assemble_elemmatvec

          if ( .not. constant_element_size ) then
            deallocate ( elemmat, elemvec, posr, posc )
            if ( lusephysqmask ) deallocate ( elemmatzeros, pqr, pqc )
          end if

        end do

        if ( constant_element_size ) then
          deallocate ( elemmat, elemvec, posr, posc )
          if ( lusephysqmask ) deallocate ( elemmatzeros, pqr, pqc )
        end if

      end do

    end if

    if ( problem%numtransdegfd > 0. .and. ltransform ) then
      deallocate ( work1 )
    end if

    deallocate( w1, w2 )

    if ( .not. constant_element_size .or. lexclude_single_layer ) then
      deallocate ( nlw )
    end if

    if ( present(elementset) ) then
      do elgrp = 1, mesh%nelgrp
        deallocate ( buildelements(elgrp)%a )
      end do
      deallocate ( buildelements )
    end if

    deallocate ( wkl, wk1, wk2, nwk )

  contains


!   test whether element (elgrp,elem) is in layer

    function element_in_layer ()
      logical :: element_in_layer
      element_in_layer = &
          all ( wkl(mesh%topology(elgrp)%a(:,elem)) ) .or. &
          .not. llayer_in_all_nodes .and. &
               any ( wkl(mesh%topology(elgrp)%a(:,elem)) )
    end function element_in_layer

!   test whether element (elgrp,elem) is in a single layer

    function element_in_single_layer ()
      logical :: element_in_single_layer
      element_in_single_layer = &
          all ( nlw(mesh%topology(elgrp)%a(:,elem)) == 1 ) .and. &
          all ( problem%nodlayers(mesh%topology(elgrp)%a(:,elem)) == &
                problem%nodlayers(mesh%topology(elgrp)%a(1,elem)) )
    end function element_in_single_layer


!   allocate the element matrix and vector (standard: based on element groups)

    subroutine allocate_elemmatvec

!     find number of degrees of freedom in an element of this group

      if ( present(physqrow) ) then
        ndofr = &
             sum ( problem%vec_elnumdegfd(elgrp)%a(:,problem%physq(physqrow)) )
      else
        ndofr = sum ( problem%elnumdegfd(elgrp)%a )
      end if

      if ( present(physqcol) ) then
        ndofc = &
             sum ( problem%vec_elnumdegfd(elgrp)%a(:,problem%physq(physqcol)) )
      else
        ndofc = sum ( problem%elnumdegfd(elgrp)%a )
      end if

!     reserve memory for element matrix and vector and positions

      allocate ( elemmat(ndofr,ndofc), elemvec(ndofr*nrhsd) )
      allocate ( posr(ndofr), posc(ndofc) )

      if ( lusephysqmask ) then
        allocate ( elemmatzeros(ndofr,ndofc) )
        allocate ( pqr(ndofr), pqc(ndofc) )
      end if

    end subroutine allocate_elemmatvec


!   allocate the element matrix and vector (variable: based on layers)

    subroutine allocate_elemmatvec_layers

!     find number of degrees of freedom in an element

      integer :: i
      integer, dimension(mesh%elnumnod(elgrp)) :: nl

      if ( problem%numlayers == 0 ) stop 'internal error: numlayers=0'

      if ( present(layer) ) then
!       single layer
        where ( wkl(mesh%topology(elgrp)%a(:,elem)) )
          nl = 1
        else where
          nl = 0
        end where
      else
!       multiple layers
        nl = nlw(mesh%topology(elgrp)%a(:,elem))
      end if

      if ( present(physqrow) ) then
        ndofr = 0
        do i = 1, size(physqrow)
          ndofr = ndofr + sum ( &
            problem%vec_elnumdegfd(elgrp)%a(:,problem%physq(physqrow(i))) * nl )
        end do
      else
        ndofr = sum ( problem%elnumdegfd(elgrp)%a * nl )
      end if

      if ( present(physqcol) ) then
        ndofc = 0
        do i = 1, size(physqcol)
          ndofc = ndofc + sum ( &
            problem%vec_elnumdegfd(elgrp)%a(:,problem%physq(physqcol(i))) * nl )
        end do
      else
        ndofc = sum ( problem%elnumdegfd(elgrp)%a * nl )
      end if

!     reserve memory for element matrix and vector and positions

      allocate ( elemmat(ndofr,ndofc), elemvec(ndofr*nrhsd) )
      allocate ( posr(ndofr), posc(ndofc) )

      if ( lusephysqmask ) then
        allocate ( elemmatzeros(ndofr,ndofc) )
        allocate ( pqr(ndofr), pqc(ndofc) )
      end if

    end subroutine allocate_elemmatvec_layers


!   compute and assemble the element matrix and vector

    subroutine compute_assemble_elemmatvec

      integer :: k, rhsd, row, i, j
      logical :: transformrows, transformcolumns
      type(sparsematrix_t) :: A

!     compute element matrix and vector

      if ( zeros ) then

        elemmat = 0
        elemvec = 0

      else if ( present(elemsub) ) then

!       simple interface

        call elemsub ( mesh, problem, elgrp, elem, matrix, vector, &
          first, last, coeffl, oldvl, elemmat, elemvec )

      else if ( present(elemsub1) ) then

!       extended interface using a type for element info

        call elemsub1 ( mesh, problem, eleminfo, matrix, vector, &
          first, last, coeffl, oldvl, elemmat, elemvec )

      end if

!     compute positions in large matrix/vector of element degrees of freedom

      if ( lusephysqmask ) then

!       row
        call pos_array ( mesh, problem, elgrp, elem, dofr, posr, physqrow, &
          order, layer, pq=pqr )

!       column
        call pos_array ( mesh, problem, elgrp, elem, dofc, posc, physqcol, &
          order, layer, pq=pqc )

!       fill elemmatzeros

        elemmatzeros = .not. problem%physqmask(pqr,pqc)

      else

!       row
        call pos_array ( mesh, problem, elgrp, elem, dofr, posr, physqrow, &
          order, layer )

!       column
        call pos_array ( mesh, problem, elgrp, elem, dofc, posc, physqcol, &
          order, layer )

      end if

      if ( dofr /= ndofr .or. dofc /= ndofc ) then
        write(*,'(/a/a/)') &
          'Internal error in build_system: ', &
          ' dofr /= ndofr .or. dofc /= ndofc'
        stop
      end if

!     transformations

      if ( problem%numtransdegfd > 0 .and. ltransform ) then

        transformrows = any ( work1(posr) )
        transformcolumns = any ( work1(posc) ) .and. matrix

        if ( transformcolumns ) then

!         extract transformation matrix

          A = extract_submat ( problem%Amat, posc, posc )

!         transform matrix elementmatrix * A

          if ( matrix ) elemmat = fsmatmul ( elemmat, A )

          call delete(A)

        end if

        if ( transformrows ) then

!         extract transformation matrix

          A = extract_submat ( problem%Amat, posr, posr )

!         transform matrix A^T * elementmatrix

          if ( matrix ) elemmat = stfmatmul ( A, elemmat )

!         transform vector A^T * vector

          if ( vector ) then

!           possibly multiple right-hand side

            do rhsd = 1, nrhsd
              k = dofr*(rhsd-1)
              elemvec(k+1:k+dofr) = stmatvec ( A, elemvec(k+1:k+dofr) )
            end do

          end if

          call delete(A)

        end if

      end if

!     add matrix

      if ( matrix ) then

!       add element matrix row by row

        if ( present(factormat) ) elemmat = factormat * elemmat

        if ( lusephysqmask ) then
          call add_elemmat_to_sysmatrix ( sysmatrix, elemmat, posr, posc, &
            w1, w2, elemmatzeros )
        else
          call add_elemmat_to_sysmatrix ( sysmatrix, elemmat, posr, posc, &
            w1, w2 )
        end if

      end if

!     add vector

      if ( vector ) then

!       add element vector to large vector

        if ( present(factorvec) ) elemvec = factorvec * elemvec

        if ( present(msysvector) ) then

!         multiple right-hand side

          do row = 1, dofr

            rowg = posr(row) ! global row number

            do rhsd = 1, nrhsd
              msysvector(rhsd)%u(rowg) = msysvector(rhsd)%u(rowg) &
                                          + elemvec( row + dofr*(rhsd-1) )
            end do

          end do

        else if ( present(m2sysvector) ) then

!         multiple right-hand side (matrix)

          si = size(m2sysvector,1)

          do row = 1, dofr

            rowg = posr(row) ! global row number

            do j = 1, size(m2sysvector,2)
              do i = 1, si
                rhsd = i + si*(j-1)
                m2sysvector(i,j)%u(rowg) = m2sysvector(i,j)%u(rowg) &
                                              + elemvec( row + dofr*(rhsd-1) )
              end do
            end do

          end do

        else if ( present(m3sysvector) ) then

!         multiple right-hand side (3D array)

          si = size(m3sysvector,1)
          sj = size(m3sysvector,2)

          do row = 1, dofr

            rowg = posr(row) ! global row number

            do k = 1, size(m3sysvector,3)
              do j = 1, sj
                do i = 1, si
                  rhsd = i + si*(j-1) + si*sj*(k-1)
                  m3sysvector(i,j,k)%u(rowg) = m3sysvector(i,j,k)%u(rowg) &
                                              + elemvec( row + dofr*(rhsd-1) )
                end do
              end do
            end do

          end do

        else

!         single right-hand side

          do row = 1, dofr

            rowg = posr(row) ! global row number

            sysvector%u(rowg) = sysvector%u(rowg) + elemvec(row)

          end do

        end if

      end if

    end subroutine compute_assemble_elemmatvec

  end subroutine build_system


! find mesh elements and integration points in an element of an object

  subroutine find_elements_integrationpoints ( mesh, object, elemo, nwkelem, &
    wk1, wk2, nwk, intps )

    type(mesh_t), intent(in) :: mesh

!   object number and element number in object
    integer, intent(in) :: object, elemo

!   number of mesh elements found that need to be integrated/assembled
    integer, intent(out) :: nwkelem

!   group and element numbers found
!   nwk(i,:) = /(elgrp,elem)/, i=1,nwkelem
    integer, dimension(:,:), intent(out) :: nwk

!   work storage for mesh elements
    integer, dimension(:,:), intent(inout) :: wk1, wk2

!   storage for integration point number in the elements
!   intps(i)%a, i=1,nwkelem gives the integration point numbers for element i
    type(int_array_1d_p), allocatable, dimension(:), intent(out) :: intps


    integer, allocatable, dimension(:) :: numi
    integer :: intp, elgrp, elem, posel, i, ninti


    ninti = mesh%objects(object)%ninti

!   find mesh elements and number of integration points first (first pass)

    nwkelem = 0

    do intp = 1, ninti

      elgrp = mesh%objects(object)%grpelm_int(intp,1,elemo)
      elem  = mesh%objects(object)%grpelm_int(intp,2,elemo)

      if ( wk1(elgrp,elem) == 0 ) then
!       new element found
        nwkelem = nwkelem + 1
        nwk(nwkelem,:) = [elgrp,elem]
        wk2(elgrp,elem) = nwkelem  ! store number
      end if
      wk1(elgrp,elem) = wk1(elgrp,elem) + 1  ! increment number of intps

    end do

    allocate ( intps(nwkelem), numi(nwkelem) )
    do i = 1, nwkelem
      allocate ( intps(i)%a(wk1(nwk(i,1),nwk(i,2))) )
    end do
    numi = 0

!   find integration points of each element (second pass)

    do intp = 1, ninti

      elgrp = mesh%objects(object)%grpelm_int(intp,1,elemo)
      elem  = mesh%objects(object)%grpelm_int(intp,2,elemo)

      posel = wk2(elgrp,elem)
      numi(posel) = numi(posel) + 1
      intps(posel)%a(numi(posel)) = intp

    end do

    deallocate(numi)

  end subroutine find_elements_integrationpoints


! reset mesh elements and integration points in the mesh elements

  subroutine reset_elements_integrationpoints ( nwkelem, wk1, wk2, nwk, intps )

    integer, intent(inout) :: nwkelem

    integer, dimension(:,:), intent(inout) :: wk1, wk2, nwk

    type(int_array_1d_p), allocatable, dimension(:), intent(inout) :: intps

    integer :: i

    do i = 1, nwkelem
      deallocate ( intps(i)%a )
    end do
    deallocate ( intps )

!   set back elements counted in objects

    do i = 1, nwkelem
      wk1(nwk(i,1),nwk(i,2)) = 0
      wk2(nwk(i,1),nwk(i,2)) = 0
      nwk(i,:) = 0
    end do
    nwkelem = 0

  end subroutine reset_elements_integrationpoints


! Add element matrix to system matrix

  subroutine add_elemmat_to_sysmatrix ( sysmatrix, elemmat, posr, posc, &
    w1, w2, elemmatzeros )

!   the system matrix
    type(sysmatrix_t), intent(inout) :: sysmatrix

!   the element matrix
    real(dp), intent(in), dimension(:,:) :: elemmat

!   the positions of the global degrees of freedom for row and column
    integer, intent(in), dimension(:) :: posr, posc

!   the work arrays of minimal size numund and numess, respectively
    integer, intent(inout), dimension(:) :: w1, w2

!   The optional matrix elemmatzeros has the same shape as elemmat
!   but is a logical. It gives the entries in the matrix that
!   are zero (.true.) or .non-zero (.false.). If an entry = .true. the
!   corresponding matrix entry is not added to the system matrix.
    logical, intent(in), dimension(:,:), optional :: elemmatzeros

    integer :: numund, row, col, rowg, colg, rowgr
    integer :: jcol, pa, puu, pup, ppu, ppp


    if ( size(posr) == 0 .or. size(posc) == 0 ) return

    numund = sysmatrix%Suu%m

!   add element matrix row by row

    do row = 1, size(posr)

      rowg = posr(row) ! global row number

      if ( rowg <= numund ) then

!       row of unknown

!       fill w1 for non-zeros already in matrix Suu

!       puu is position in the row of Suu for new data
        puu = sysmatrix%Suu%ia(rowg+1)

        do pa = sysmatrix%Suu%ia(rowg), sysmatrix%Suu%ia(rowg+1) - 1

          jcol = sysmatrix%Suu%ja(pa)

          if ( jcol == 0 ) then
!           matrix row not completely filled
            puu = pa
            exit
          end if

          w1(jcol) = pa  ! put in w1: position in a of non-zero entry

        end do

!       fill w2 for non-zeros already in matrix Sup

!       pup is position in the row of Sup for new data
        pup = sysmatrix%Sup%ia(rowg+1)

        do pa = sysmatrix%Sup%ia(rowg), sysmatrix%Sup%ia(rowg+1) - 1

          jcol = sysmatrix%Sup%ja(pa)

          if ( jcol == 0 ) then
!           matrix row not completely filled
            pup = pa
            exit
          end if

          w2(jcol) = pa  ! put in w2: position in a of non-zero entry

        end do

!       add element row (row) to global matrix row (rowg)

        do col = 1, size(posc)

          if ( present(elemmatzeros) ) then
            if ( elemmatzeros(row,col) ) cycle
          end if

          colg = posc(col)

          if ( colg <= numund ) then

!           column of unknown

            if ( sysmatrix%symmetric ) then
!             symmetric matrix
              if ( rowg > colg ) cycle  ! ignore lower triangle
            end if

            pa = w1(colg)

            if ( pa == 0 .and. puu /= sysmatrix%Suu%ia(rowg+1) ) then

!             put new non-zero in matrix

              sysmatrix%Suu%a(puu)  = elemmat(row,col)
              sysmatrix%Suu%ja(puu) = colg

              puu = puu + 1

            else if ( pa /= 0 ) then

!             add to existing non-zero

              sysmatrix%Suu%a(pa) = sysmatrix%Suu%a(pa) + elemmat(row,col)

            else if ( pa == 0 .and. puu == sysmatrix%Suu%ia(rowg+1) ) then

!             new non-zero but no free space

              write(*,'(/a,i0/a/)') &
                'Internal Error: no more space in row ', rowg, &
                'of system matrix Suu to add element matrix.'
              stop

            end if

          else if ( colg > numund ) then

!           column of prescribed degree of freedom

            pa = w2(colg-numund)

            if ( pa == 0 .and. pup /= sysmatrix%Sup%ia(rowg+1) ) then

!             put new non-zero in matrix

              sysmatrix%Sup%a(pup)  = elemmat(row,col)
              sysmatrix%Sup%ja(pup) = colg - numund

              pup = pup + 1

            else if ( pa /= 0 ) then

!             add to existing non-zero

              sysmatrix%Sup%a(pa) = sysmatrix%Sup%a(pa) + elemmat(row,col)

            else if ( pa == 0 .and. pup == sysmatrix%Sup%ia(rowg+1) ) then

!             new non-zero but no free space

              write(*,'(/a,i0/a/)') &
                'Internal Error: no more space in row ', rowg, &
                'of system matrix Sup to add element matrix.'
              stop

            end if

          end if

        end do

!       put w1 and w2 back to zero

        do pa = sysmatrix%Suu%ia(rowg), sysmatrix%Suu%ia(rowg+1) - 1

          jcol = sysmatrix%Suu%ja(pa)

          if ( jcol == 0 ) exit

          w1(jcol) = 0

        end do

        do pa = sysmatrix%Sup%ia(rowg), sysmatrix%Sup%ia(rowg+1) - 1

          jcol = sysmatrix%Sup%ja(pa)

          if ( jcol == 0 ) exit

          w2(jcol) = 0  ! put in w2: position in a of non-zero entry

        end do

      else if ( rowg > numund ) then

!       row of prescribed degree of freedom

        rowgr = rowg - numund

!       fill w1 for non-zeros already in matrix Spu

!       ppu is position in the row of Spu for new data
        ppu = sysmatrix%Spu%ia(rowgr+1)

        do pa = sysmatrix%Spu%ia(rowgr), sysmatrix%Spu%ia(rowgr+1) - 1

          jcol = sysmatrix%Spu%ja(pa)

          if ( jcol == 0 ) then
!           matrix row not completely filled
            ppu = pa
            exit
          end if

          w1(jcol) = pa  ! put in w1: position in a of non-zero entry

        end do

!       fill w2 for non-zeros already in matrix Spp

!       ppp is position in the row of Spp for new data
        ppp = sysmatrix%Spp%ia(rowgr+1)

        do pa = sysmatrix%Spp%ia(rowgr), sysmatrix%Spp%ia(rowgr+1) - 1

          jcol = sysmatrix%Spp%ja(pa)

          if ( jcol == 0 ) then
!           matrix row not completely filled
            ppp = pa
            exit
          end if

          w2(jcol) = pa  ! put in w2: position in a of non-zero entry

        end do

!       add element row (row) to global matrix row (rowg)

        do col = 1, size(posc)

          if ( present(elemmatzeros) ) then
            if ( elemmatzeros(row,col) ) cycle
          end if

          colg = posc(col)

          if ( colg <= numund ) then

!           column of unknown

            pa = w1(colg)

            if ( pa == 0 .and. ppu /= sysmatrix%Spu%ia(rowgr+1) ) then

!             put new non-zero in matrix

              sysmatrix%Spu%a(ppu)  = elemmat(row,col)
              sysmatrix%Spu%ja(ppu) = colg

              ppu = ppu + 1

            else if ( pa /= 0 ) then

!             add to existing non-zero

              sysmatrix%Spu%a(pa) = sysmatrix%Spu%a(pa) + elemmat(row,col)

            else if ( pa == 0 .and. &
                      ppu == sysmatrix%Spu%ia(rowgr+1) ) then

!             new non-zero but no free space

              write(*,'(/a,i0/a/)') &
                'Internal Error: no more space in row ', rowg, &
                'of system matrix Spu to add element matrix.'
              stop

            end if

          else if ( colg > numund ) then

!           column of prescribed degree of freedom

            pa = w2(colg-numund)

            if ( pa == 0 .and. ppp /= sysmatrix%Spp%ia(rowgr+1) ) then

!             put new non-zero in matrix

              sysmatrix%Spp%a(ppp)  = elemmat(row,col)
              sysmatrix%Spp%ja(ppp) = colg - numund

              ppp = ppp + 1

            else if ( pa /= 0 ) then

!             add to existing non-zero

              sysmatrix%Spp%a(pa) = sysmatrix%Spp%a(pa) + elemmat(row,col)

            else if ( pa == 0 .and. &
                      ppp == sysmatrix%Spp%ia(rowgr+1) ) then

!             new non-zero but no free space

              write(*,'(/a,i0/a/)') &
                'Internal Error: no more space in row ', rowg, &
                'of system matrix Spp to add element matrix.'
              stop

            end if

          end if

        end do

!       put w1 and w2 back to zero

        do pa = sysmatrix%Spu%ia(rowgr), sysmatrix%Spu%ia(rowgr+1) - 1

          jcol = sysmatrix%Spu%ja(pa)

          if ( jcol == 0 ) exit

          w1(jcol) = 0

        end do

        do pa = sysmatrix%Spp%ia(rowgr), sysmatrix%Spp%ia(rowgr+1) - 1

          jcol = sysmatrix%Spp%ja(pa)

          if ( jcol == 0 ) exit

          w2(jcol) = 0  ! put in w2: position in a of non-zero entry

        end do

      end if

    end do

  end subroutine add_elemmat_to_sysmatrix


! Add a one-dimensional vector to the system matrix, for example a diagonal.

  subroutine add_elemvec_to_sysmatrix ( sysmatrix, elemvec, posr, posc, &
    w1, w2 )

!   the system matrix
    type(sysmatrix_t), intent(inout) :: sysmatrix

!   the element matrix
    real(dp), intent(in), dimension(:) :: elemvec

!   the positions of the global degrees of freedom for row and column in
!   the element vector elemvec. Thus: elemvec(i) will be added to the
!   the matrix entry A(posr(i),posc(i)).
    integer, intent(in), dimension(:) :: posr, posc

!   the work arrays of minimal size numund and numess, respectively
    integer, intent(inout), dimension(:) :: w1, w2


!   Add a one-dimensional vector to the system matrix, for example a diagonal.
!   This routine uses the trick of filling the diagonal of an element matrix
!   and combine that with filling of elemmatzeros=.true. for the off-diagonal
!   terms. Therefore this routine is not suitable for large vectors since this
!   would require a full matrix of the (size vector)**2.


    integer :: row
    real(dp) :: elemmat(size(posr),size(posr))
    logical :: elemmatzeros(size(posr),size(posr))


    if ( size(posr) /= size(posc) ) then
      write(*,'(/a/)') 'Internal error: size(posr) =/ size(posc) '
      stop
    end if

    if ( size(posr) == 0 ) return

!   fill elemmat and elemmatzeros

    elemmatzeros = .true.

    do row = 1, size(posr)

!      put vector on diagonal of elemmat

       elemmat(row,row) = elemvec(row)
       elemmatzeros(row,row) = .false.

    end do

!   put diagonal of elemmat in matrix

    call add_elemmat_to_sysmatrix ( sysmatrix, elemmat, posr, posc, w1, w2, &
      elemmatzeros )

  end subroutine add_elemvec_to_sysmatrix


! Check whether sysmatrix is fully filled

  subroutine check_filled_sysmatrix ( sysmatrix )

    type(sysmatrix_t), intent(in) :: sysmatrix

    if ( any(sysmatrix%Suu%ja == 0) .or.  any(sysmatrix%Sup%ja == 0) .or. &
         any(sysmatrix%Spu%ja == 0) .or.  any(sysmatrix%Spp%ja == 0) ) then
      write(*,'(/a/)') &
        'Error check_filled_sysmatrix: sysmatrix not fully filled.'
      stop
    end if

  end subroutine check_filled_sysmatrix


! Add boundary elements to sysvector and (optionally) sysmatrix (in a point)

  subroutine add_boundary_elements_point ( mesh, problem, sysvector, point, &
    elemsub, sysmatrix, coefficients, oldvectors, physq, factorvec, factormat, &
    layer, buildvector, buildmatrix, transform )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the system vector to be assembled
!   NOTE: element vectors are added to the existing vector!
    type(sysvector_t), intent(inout) :: sysvector

!   the assembling of the boundary element is for this point only
    integer, intent(in) :: point

!   this is the element subroutine that must be supplied by the calling routine
    interface
      subroutine elemsub ( mesh, problem, point, matrix, vector, &
        coefficients, oldvectors, elemmat, elemvec )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: point
        logical, intent(in) :: matrix, vector
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat
        real(dp), intent(out), dimension(:) :: elemvec
      end subroutine elemsub
    end interface

!   if present, the system matrix to be assembled
!   NOTE: element matrices are added to the existing matrix!
    type(sysmatrix_t), intent(inout), optional :: sysmatrix

!   if coefficients is present, it is passed through to the element subroutine
    type(coefficients_t), intent(in), optional :: coefficients

!   if oldvectors is present, it is passed through to the element subroutine.
!   This structure can be used to supply old sysvectors/vectors from possibly
!   other problems to the element subroutines. Note that the problem structures
!   of the other problems must be included in oldvectors.
    type(oldvectors_t), intent(in), optional :: oldvectors

!   if present, physq specifies the physical quantities involved in the
!   boundary element. For example physq=(/2,3)/, means physical quanties 2
!   and 3 are involved. If not present all degrees of freedom are involved.
    integer, intent(in), dimension(:), optional :: physq

!   if present the element matrix and/or vector will be multiplied by the
!   specified factor before assembly into the system matrix/vector.
    real(dp), intent(in), optional :: factormat, factorvec

!   if present the assembling takes place within the single layer of degrees of
!   freedom only.
!   Only a point where the node has degrees in the specified layer is
!   assembled.
    integer, intent(in), optional :: layer

!   if assigned the value .false. the assembling of the matrix or vector
!   is not done. Note, that if sysmatrix is not present, buildmatrix is ignored.
    logical, intent(in), optional :: buildmatrix, buildvector

!   By setting transform=.false. the transformation matrix is not applied
!   to the element matrix and vector. This means, the element degrees of
!   freedom are defined in the local (transformed) system.
!   default = .true.
    logical, intent(in), optional :: transform


!   More points can be added by calling this routine for another point.
!   This routine performs the assembly proces for a single points only and
!   nothing more. This means that this routine must
!   be called for each point separately.
!   Other parameters, such as material parameters, function numbers,
!   old vectors etc., can be supplied to the element routine by either
!   modules variables (global module variables) or via the heading parameters
!   coefficients and oldvectors.

!
!   This routine is meant for implementing standard Neumann and Robin boundary
!   conditions. If more complicated boundary conditions, such as using
!   derivatives, are needed, the point should be converted to an
!   object and the system needs to be build on the object using build_system.


    logical :: matrix, vector, donode, ltransform
    integer :: ndof, dof, nodenr
    integer :: row, rowg, numess, numund
    integer, allocatable, dimension(:) :: pos, w1, w2
!   work1: array for denoting transformed degrees of freedom
    logical, allocatable, dimension(:) :: work1
    real(dp), allocatable, dimension(:,:) :: elemmat
    real(dp), allocatable, dimension(:) :: elemvec
    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl
    type(sparsematrix_t) :: A

    call check ( mesh, 'add_boundary_elements_point' )
    call check ( problem, 'add_boundary_elements_point', mesh )

!   test valid point

    if ( point < 1 .or. point > mesh%npoints ) then

      write(*,'(/a/a,i0,/a,i0/)') &
        'Error: point in the heading of add_boundary_elements_point is ', &
        'out of range. point is ', point, &
        'whereas the number of points is ', mesh%npoints
      stop

    end if

!   compute matrix and vector?

    if ( present(buildmatrix) ) then
      matrix = buildmatrix .and. present(sysmatrix)
    else
      matrix = present(sysmatrix)
    end if

    if ( present(buildvector) ) then
      vector = buildvector
    else
      vector = .true.
    end if

    if ( .not. ( matrix .or. vector ) ) return

    if ( matrix ) then

      if ( .not. sysmatrix%initialized_structure ) then
        write(*,'(/a/)') &
          'Error in add_boundary_elements_point: no system matrix structure.'
        stop
      end if

      if ( .not. sysmatrix%finalized ) then
        write(*,'(/2(a/))') &
          'Error in add_boundary_elements_point:', &
          ' system matrix has not been finalized.'
        stop
      end if

      if ( .not. sysmatrix%allocated_data ) then
        write(*,'(/2(a/))') &
          'Error in add_boundary_elements_point:', &
          ' data in system matrix not allocated.'
        stop
      end if

      if ( sysmatrix%Suu%m + sysmatrix%Sup%m /= problem%numdegfd ) then
        write(*,'(/a/a/a/)') &
          'Error in add_boundary_elements_point:', &
          ' number of degrees of freedom of ', &
          ' the system matrix is different from the number in problem.'
        stop
      end if

    end if

    if ( vector ) then

!     test sysvector

      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in add_boundary_elements_point: sysvector vector not created.'
        stop
      end if

      if ( sysvector%n /= problem%numdegfd ) then
        write(*,'(/3(a/))') &
          'Error in add_boundary_elements_point: ', &
          ' number of degrees of freedom ', &
          'of the sysvector vector different from the number in problem.'
        stop
      end if

    end if

!   check physq

    if ( present(physq) ) then

      if ( problem%nphysq == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in add_boundary_elements_point: ',&
          ' physq is only availaba if physical quantities have been defined.'
        stop
      end if

    end if

    if ( present(physq) ) then

      if ( any( physq < 1 ) .or.  any( physq > problem%nphysq ) ) then
        write(*,'(2(/a)/)') &
          'Error in add_boundary_elements_point: ', &
          ' physical quantities in physq are out of range.'
        stop
      end if

    end if

!   layers

    if ( present(layer) ) then
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in add_boundary_elements_point: layer out of range.'
        stop
      end if
    else if ( problem%numlayers > 0 ) then
      write(*,'(3(/a)/)') &
        'Error in add_boundary_elements_point: ', &
        ' layers have been defined and the layer keyword is not present', &
        ' layers not implemented for this routine without the layer keyword'
      stop
    end if

!   transformations

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in add_boundary_elements_point: ', &
        ' transformation matrix has not been build'
      stop
    end if

    ltransform = set_optional ( variable=transform, default=.true. )

    if ( problem%numtransdegfd > 0 .and. ltransform ) then

!     work array to indicate transformed degrees of freedom
      allocate ( work1(problem%numdegfd) )
      work1 = .false.
      work1(problem%degfdtrans) = .true.

    end if

!   initialize

    if ( present(oldvectors) ) oldvl = oldvectors
    if ( present(coefficients) ) coeffl = coefficients

!   temporary work arrays to store information on current large matrix row

    numess = problem%numessdegfd
    numund = problem%numundegfd

    allocate( w1(numund), w2(numess) )

    w1 = 0
    w2 = 0

!   find number of degrees of freedom in this point

    nodenr = mesh%points(point)

    if ( present(physq) ) then
      ndof = sum ( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physq)) &
                 - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physq)) )
    else
      ndof = problem%nodnumdegfd(nodenr+1) - problem%nodnumdegfd(nodenr)
    end if

!   reserve memory for element vector and matrix and positions

    allocate ( elemmat(ndof,ndof), elemvec(ndof), pos(ndof) )

    if ( present(layer) ) then
      donode = btest(problem%nodlayers(nodenr),layer-1) ! layer present
    else
      donode = .true.
    end if

    if ( donode ) then

!     compute element vector

      call elemsub ( mesh, problem, point, matrix, vector, coeffl, oldvl, &
        elemmat, elemvec )

!     compute positions in large vector of element degrees of freedom

      call pos_array_node ( problem, nodenr, dof, pos, physq, layer )

      if ( dof /= ndof ) then
        write(*,'(/a/a/)') &
          'Internal error in add_boundary_elements_point ', &
          ' dof /= ndof '
        stop
      end if

!     transformations

      if ( problem%numtransdegfd > 0 .and. ltransform ) then

        if ( any ( work1(pos) ) ) then

!         extract transformation matrix

          A = extract_submat ( problem%Amat, pos, pos )

!         transform matrix: A^T * elementmatrix * A

          if ( matrix ) then
            elemmat = fsmatmul ( elemmat, A )
            elemmat = stfmatmul ( A, elemmat )
          end if

!         transform vector: A^T * vector

          if ( vector ) then
            elemvec = stmatvec ( A, elemvec )
          end if

          call delete(A)

        end if

      end if

!     add matrix

      if ( matrix ) then

!       add element matrix row by row

        if ( present(factormat) ) elemmat = factormat * elemmat

        call add_elemmat_to_sysmatrix ( sysmatrix, elemmat, pos, pos, w1, w2 )

      end if

!     add vector

      if ( vector ) then

!       add element vector to large vector

        if ( present(factorvec) ) elemvec = factorvec * elemvec

        do row = 1, dof

          rowg = pos(row) ! global row number

          sysvector%u(rowg) = sysvector%u(rowg) + elemvec(row)

        end do

      end if

    end if

    deallocate ( elemmat, elemvec, pos )

    deallocate( w1, w2 )

    if ( problem%numtransdegfd > 0. .and. ltransform ) then
      deallocate ( work1 )
    end if

  end subroutine add_boundary_elements_point


! Add boundary elements (on a geometry) to sysvector and (optionally) to the
! system matrix.

  subroutine add_boundary_elements ( mesh, problem, sysvector, elemsub, &
    sysmatrix, curve, surface, coefficients, oldvectors, physq, order, &
    layer, factorvec, factormat, buildvector, buildmatrix, transform )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the system vector to be assembled
!   NOTE: element vectors are added to the existing vector!
    type(sysvector_t), intent(inout) :: sysvector

!   this is the element subroutine that must be supplied by the calling routine
    interface
      subroutine elemsub ( mesh, problem, geom, elem, matrix, vector, &
        first, last, coefficients, oldvectors, elemmat, elemvec )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: geom, elem
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat
        real(dp), intent(out), dimension(:) :: elemvec
      end subroutine elemsub
    end interface

!   if present, the system matrix to be assembled
!   NOTE: element matrices are added to the existing matrix!
    type(sysmatrix_t), intent(inout), optional :: sysmatrix

!   the assembling of the boundary elements are for this curve only.
    integer, intent(in), optional :: curve

!   the assembling of the boundary elements are for this surface only.
    integer, intent(in), optional :: surface

!   if coefficients is present, it is passed through to the element subroutine
    type(coefficients_t), intent(in), optional :: coefficients

!   if oldvectors is present, it is passed through to the element subroutine.
!   This structure can be used to supply old sysvectors/vectors from possibly
!   other problems to the element subroutines. Note that the problem structures
!   of the other problems must be included in oldvectors.
    type(oldvectors_t), intent(in), optional :: oldvectors

!   if present, physq specifies the physical quantities involved in the
!   boundary element. For example physq=(/2,3)/, means physical quanties 2
!   and 3 are involved. If not present all degrees of freedom are involved.
    integer, intent(in), dimension(:), optional :: physq

!   The parameter order determines the sequence of the degrees of freedom on
!   elementlevel. There are only two possibilities and affect only the two loops
!   over nodal points and degrees of freedom (within a physical quantity if
!   problem%nphysq>0):
!     order = 'ND' : the most inner loop is over the degrees of freedom
!     order = 'DN' : the most inner loop is over the nodal points
!   Specifically we have:   nphysq = 0   nphysq > 0
!         order = 'ND' :       ND           PND
!         order = 'DN' :       DN           PDN
!   NOTE: if nphys >0 the main ordering (outside loop) of the element
!   degrees of freedom are the physical quantities.
!   If problem%numlayers > 0 the degrees are stacked into layers and we have
!                           nphysq = 0   nphysq > 0
!         order = 'ND' :      NLD           PNLD
!         order = 'DN' :      LDN           PLDN
!   See also the userguide for further explanation.
!
!   For example consider an element with two nodes and two physical quanties:
!   one with two degrees (u,v) and one with a single degree p. Then we have:
!     order = 'ND' : u1, v1, u2, v2, p1, p2
!     order = 'DN' : u1, u2, v1, v2, p1, p2
!   The parameter order generates a permutation of the degrees of freedom on
!   element level, making the life of an element programmer easier,
!   but _does not affect_ the global numbering of the unknowns.
!   The (column) layout of the assembled sparse matrix can be different due to
!   the different sequence of unknowns in the assembling process.
!
!   The default of order is:
!     order = 'ND' : if no physical quantities have been defined (nphysq=0)
!     order = 'DN' : if physical quantities have been defined (nphysq>0)
    character(len=*), intent(in), optional :: order

!   if present the assembling takes place within the single layer of degrees of
!   freedom only.
!   Only elements where all nodes have degrees in the specified layer are
!   assembled.
    integer, intent(in), optional :: layer

!   if present the element matrix and/or vector will be multiplied by the
!   specified factor before assembly into the system matrix/vector.
    real(dp), intent(in), optional :: factormat, factorvec

!   if assigned the value .false. the assembling of the matrix or vector
!   is not done. Note, that if sysmatrix is not present, buildmatrix is ignored.
    logical, intent(in), optional :: buildmatrix, buildvector

!   By setting transform=.false. the transformation matrix is not applied
!   to the element matrix and vector. This means, the element degrees of
!   freedom are defined in the local (transformed) system.
!   default = .true.
    logical, intent(in), optional :: transform

!   More curves/surfaces can be added by calling this routine for another
!   curve/surface.
!   This routine performs the assembly proces of the boundary elements
!   for a single curve/surface only and nothing more. This means that this
!   routine must be called for each curve/surface separately.
!   Other parameters, such as material parameters, function numbers,
!   old vectors etc., can be supplied to the element routine by either
!   modules variables (global module variables) or via the heading parameters
!   coefficients and oldvectors.
!
!   This routine is meant for implementing standard Neumann and Robin boundary
!   conditions. If more complicated boundary conditions, such as using
!   derivatives, are needed, the curve/surface should be converted to an
!   object and the system needs to be build on the object using build_system.


    logical :: first, last, matrix, vector, ltransform
    logical, allocatable, dimension(:) :: wkl
    integer :: elem, ndof, dof, geom
    integer :: row, rowg, numess, numund
    integer, allocatable, dimension(:) :: pos, nodes, w1, w2
!   work1: array for denoting transformed degrees of freedom
    logical, allocatable, dimension(:) :: work1
    real(dp), allocatable, dimension(:,:) :: elemmat
    real(dp), allocatable, dimension(:) :: elemvec
    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl
    type(geometry_t) :: geometry
    type(sparsematrix_t) :: A

    allocate ( wkl(mesh%nnodes) )

    call check ( mesh, 'add_boundary_elements' )
    call check ( problem, 'add_boundary_elements', mesh )

!   test valid geometry

    if ( present(curve) .and. present(surface) ) then
      write(*,'(/a/a/)') &
        'Error: curve and surface cannot be both present in the heading', &
        ' of add_boundary_elements'
      stop
    end if

    if ( present(curve) ) then

      if ( curve < 1 .or. curve > mesh%ncurves ) then

        write(*,'(/a/a,i0,/a,i0/)') &
          'Error: curve in the heading of add_boundary_elements is ', &
          'out of range. curve is ', curve, &
          'whereas the number of curves is ', mesh%ncurves
        stop

      end if

      allocate(nodes(mesh%curves(curve)%elnumnod))

    else if ( present(surface) ) then

      if ( surface < 1 .or. surface > mesh%nsurfaces ) then

        write(*,'(/a/a,i0,/a,i0/)') &
          'Error: surface in the heading of add_boundary_elements is ', &
          'out of range. surface is ', surface, &
          'whereas the number of surfaces is ', mesh%nsurfaces
        stop

      end if

      allocate(nodes(mesh%surfaces(surface)%elnumnod))

    else

      write(*,'(/a/a/)') &
        'Error: either curve or surface must be present in the heading', &
        ' of add_boundary_elements'
      stop

    end if

!   compute matrix and vector?

    if ( present(buildmatrix) ) then
      matrix = buildmatrix .and. present(sysmatrix)
    else
      matrix = present(sysmatrix)
    end if

    if ( present(buildvector) ) then
      vector = buildvector
    else
      vector = .true.
    end if

    if ( .not. ( matrix .or. vector ) ) return

    if ( matrix ) then

      if ( .not. sysmatrix%initialized_structure ) then
        write(*,'(/a/)') &
          'Error in add_boundary_elements: no system matrix structure.'
        stop
      end if

      if ( .not. sysmatrix%finalized ) then
        write(*,'(/2a/)') &
          'Error in add_boundary_elements: ', &
          'system matrix has not been finalized.'
        stop
      end if

      if ( .not. sysmatrix%allocated_data ) then
        write(*,'(/a/)') &
          'Error in add_boundary_elements: data in system matrix not allocated.'
        stop
      end if

      if ( sysmatrix%Suu%m + sysmatrix%Sup%m /= problem%numdegfd ) then
        write(*,'(/a/a/)') &
          'Error in add_boundary_elements: number of degrees of freedom of ', &
          ' the system matrix is different from the number in problem.'
        stop
      end if

    end if

    if ( vector ) then

!     test sysvector

      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in add_boundary_elements: sysvector vector not created.'
        stop
      end if

      if ( sysvector%n /= problem%numdegfd ) then
        write(*,'(/a/a/)') &
          'Error in add_boundary_elements: number of degrees of freedom ', &
          'of the sysvector vector different from the number in problem.'
        stop
      end if

    end if

!   ordering; set defaults

    if ( present(order) ) then
      if ( all ( order /= [ 'ND', 'DN' ] ) ) then
        write(*,'(2(/a)/)') &
          'Error in add_boundary_elements: ', &
          ' heading parameter order must be either ''ND'' or ''DN''.'
        stop
      end if
    end if

!   check physq

    if ( present(physq) ) then

      if ( problem%nphysq == 0 ) then
        write(*,'(4(/a)/)') &
          'Error in add_boundary_elements: ',&
          ' using physq or ordering of the degrees of freedom using order', &
          ' is only available if physical quantities have been defined.'
        stop
      end if

    end if

    if ( present(physq) ) then
      if ( any( physq < 1 ) .or.  any( physq > problem%nphysq ) ) then
        write(*,'(2(/a)/)') &
          'Error in add_boundary_elements: ', &
          ' physical quantities in physq are out of range.'
        stop
      end if
    end if

!   layers

    if ( present(layer) ) then
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in add_boundary_elements: layer out of range.'
        stop
      end if
      wkl = btest(problem%nodlayers,layer-1)
    else if ( problem%numlayers > 0 ) then
      write(*,'(3(/a)/)') &
        'Error in add_boundary_elements: ', &
        ' layers have been defined and the layer keyword is not present', &
        ' layers not implemented for this routine without the layer keyword'
      stop
    end if

!   transformations

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in add_boundary_elements: ', &
        ' transformation matrix has not been build'
      stop
    end if

    ltransform = set_optional ( variable=transform, default=.true. )

    if ( problem%numtransdegfd > 0 .and. ltransform ) then

!     work array to indicate transformed degrees of freedom
      allocate ( work1(problem%numdegfd) )
      work1 = .false.
      work1(problem%degfdtrans) = .true.

    end if

!   initialize

    if ( present(oldvectors) ) oldvl = oldvectors
    if ( present(coefficients) ) coeffl = coefficients

!   temporary work arrays to store information on current large matrix row

    numess = problem%numessdegfd
    numund = problem%numundegfd

    allocate( w1(numund), w2(numess) )

    w1 = 0
    w2 = 0

!   find geometry and number of degr of freedom in an element of this geometry

    if ( present(curve) ) then
      geom = curve
      nodes = mesh%curves(curve)%topology(:,1,2)
      geometry = mesh%curves(curve)
    else if ( present(surface) ) then
      geom = surface
      nodes = mesh%surfaces(surface)%topology(:,1,2)
      geometry = mesh%surfaces(surface)
    end if

    if ( present(physq) ) then
      ndof = sum ( problem%vec_nodnumdegfd(nodes+1,&
                                      &problem%physq(physq)) &
                       - problem%vec_nodnumdegfd(nodes,&
                                      &problem%physq(physq)) )
    else
      ndof = sum ( problem%nodnumdegfd(nodes+1) - problem%nodnumdegfd(nodes) )
    end if

!   reserve memory for element vector and matrix and positions

    allocate ( elemmat(ndof,ndof), elemvec(ndof), pos(ndof) )

!   loop over elements in this geometry

    do elem = 1, geometry%nelem

      first = elem == 1
      last  = elem == geometry%nelem

!     skip elements not in layer
      if ( present(layer) ) then
        if ( .not. all( wkl(geometry%topology(:,elem,2)) ) ) cycle
!       avoid leaving allocated memory in elements
        first = .true.; last = .true.
      end if

!     compute element vector

      call elemsub ( mesh, problem, geom, elem, matrix, vector, first, last, &
        coeffl, oldvl, elemmat, elemvec )

!     compute positions in large vector of element degrees of freedom

      if ( present(physq) ) then
!       physical quantities
        call pos_array_geometry ( problem, geometry, elem, dof, pos, &
          physq, order, layer )
      else
        call pos_array_geometry ( problem, geometry, elem, dof, pos, &
          order=order, layer=layer )
      end if

      if ( dof /= ndof ) then
        write(*,'(/a/a/)') &
          'Internal error in add_boundary_elements ', &
          ' dof /= ndof '
        stop
      end if

!     transformations

      if ( problem%numtransdegfd > 0 .and. ltransform ) then

        if ( any ( work1(pos) ) ) then

!         extract transformation matrix

          A = extract_submat ( problem%Amat, pos, pos )

!         transform matrix: A^T * elementmatrix * A

          if ( matrix ) then
            elemmat = fsmatmul ( elemmat, A )
            elemmat = stfmatmul ( A, elemmat )
          end if

!         transform vector: A^T * vector

          if ( vector ) then
            elemvec = stmatvec ( A, elemvec )
          end if

          call delete(A)

        end if

      end if

!     add matrix

      if ( matrix ) then

!       add element matrix row by row

        if ( present(factormat) ) elemmat = factormat * elemmat

        call add_elemmat_to_sysmatrix ( sysmatrix, elemmat, pos, pos, w1, w2 )

      end if

!     add vector

      if ( vector ) then

!       add element vector to large vector

        if ( present(factorvec) ) elemvec = factorvec * elemvec

        do row = 1, dof

          rowg = pos(row) ! global row number

          sysvector%u(rowg) = sysvector%u(rowg) + elemvec(row)

        end do

      end if

    end do

    deallocate ( elemmat, elemvec, pos, nodes )

    deallocate( w1, w2 )

    if ( problem%numtransdegfd > 0. .and. ltransform ) then
      deallocate ( work1 )
    end if

    deallocate ( wkl )

  end subroutine add_boundary_elements


! Add effect of essential boundary conditions to right-hand side

  subroutine add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd, &
    essentialpart )

    type(problem_t), intent(in) :: problem
    type(sysmatrix_t), intent(in) :: sysmatrix

!   the right-hand side.
!   unknown part of rhs is modified:
!     rhs_u = rhs_u - Sup * usol_p
!   the essential part of rhs is modified if essentialpart=.true.:
!     rhs_p = rhs_p - Spp * usol_p
    type(sysvector_t), intent(inout) :: rhsd

!   the vector sol must be filled with the essential boundary conditions
    type(sysvector_t), intent(in) :: sol

!   if essentialpart=.true. the essential part of rhsd is also modified
!   default=.false.
    logical, intent(in), optional :: essentialpart


!   This routine adds the effect of essential boundary conditions to the
!   right-hand side vector. This routine should be called when the assembly
!   process is fully completed. Note that rhsd is modified.


    integer :: numdegfd, numund
    logical :: esspart

    call check ( problem, 'add_effect_of_essential_to_rhs' )

    esspart = set_optional ( variable=essentialpart, default=.false. )

    numdegfd = problem%numdegfd
    numund = problem%numundegfd

!   unknown part of rhs:  rhs_u = rhs_u - Sup * usol_p

    rhsd%u(1:numund) = rhsd%u(1:numund)  &
        - smatvec ( sysmatrix%Sup, sol%u(numund+1:numdegfd) )

    if ( esspart ) then

!     essential part of rhs:  rhs_p = rhs_p - Spp * usol_p

      rhsd%u(numund+1:numdegfd) = rhsd%u(numund+1:numdegfd)  &
          - smatvec ( sysmatrix%Spp, sol%u(numund+1:numdegfd) )

    end if

  end subroutine add_effect_of_essential_to_rhs


! Compute reaction_forces

  subroutine reaction_forces ( problem, sysmatrix, sol, rhsd, reacf, &
    unknownpart )

    type(problem_t), intent(in) :: problem

!   the system matrix
    type(sysmatrix_t), intent(in) :: sysmatrix

!   the vector sol containing the solution
    type(sysvector_t), intent(in) :: sol

!   The right-hand side.
    type(sysvector_t), intent(in) :: rhsd

!   If unknownpart=.false. (default), the unknown part of reacf is set to zero:
!       reacf_u=0
!   If unknownpart=.true. the unknown part of reacf is set to the residual:
!       reacf_u = Suu * usol_u + Sup * rhs_p - rhs_u = 0
!   The result should be approximately zero. Note, that the last two terms are
!   already contained in the rhsd vector.
!   The essential part of reacf is set to the reaction forces:
!       reacf_p = Spu * usol_p + Spp * usol_p - rhs_p
!   WARNING: do not use essentialpart=.true. in the call of
!   add_effect_of_essential_to_rhs for reacf_p to represent the
!   reaction_forces.
    type(sysvector_t), intent(inout) :: reacf

!   if unknownpart=.true. the unknown part of reacf is set to the residual.
!   Otherwise it is set to zero (see above)
!   default=.false.
    logical, intent(in), optional :: unknownpart

!   The system is given by
!
!    [ Suu Sup ] [ usol_u ] = [ rhs_u ]
!    [ Spu Spp ] [ usol_p ] = [ rhs_p ]
!
!   The following system of equations is solved
!
!     Suu * usol_u = rhs_u - Sup * rhs_p
!
!   where the last term is already added using add_effect_of_essential_to_rhs.
!   The residual should be close to zero (except for some small numerical
!   errors):
!
!     reacf_u = Suu * usol_u + Sup * rhs_p - rhs_u = 0
!
!   The reaction forces in Dirichlet nodes are computed from the residual of the
!   essential part:
!
!     reacf_p = Spu * usol_u + Spp * usol_p - rhs_p
!

    integer :: numdegfd, numund
    logical :: unkpart

    call check ( problem, 'reaction_forces' )

    unkpart = set_optional ( variable=unknownpart, default=.false. )

    if ( unkpart .and. sysmatrix%symmetric ) then

      write(*,'(3(a/))') &
        'Error in residual:', &
        ' unknownpart=.true. not allowed for a ', &
        ' symmetric system matrix structure.'
      stop

    end if

    numdegfd = problem%numdegfd
    numund = problem%numundegfd

    if ( unkpart ) then

!     unknown part of reacf:  reacf_u = Suu * usol_u - rhs_u

      reacf%u(1:numund) = smatvec ( sysmatrix%Suu, sol%u(1:numund) ) &
                                - rhsd%u(1:numund)

    else

      reacf%u(1:numund) = 0

    end if

!   essential part of reacf:
!     reacf_p = Spu * usol_p + Spp * usol_p - rhs_p

    reacf%u(numund+1:numdegfd) = &
                       smatvec ( sysmatrix%Spu, sol%u(1:numund) ) &
                     + smatvec ( sysmatrix%Spp, sol%u(numund+1:numdegfd) ) &
                     - rhsd%u(numund+1:numdegfd)

  end subroutine reaction_forces


! Build transformation matrix.

  subroutine build_transformation_matrix ( mesh, problem, transformation1, &
    transformation2, coefficients, oldvectors, nodesub )

    use math_defs_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(inout) :: problem

!   if these are present the building takes place for transformations
!   transformation1,...,transformation2 only. If only transformation1 is
!   present one transformation is build only.
    integer, intent(in), optional :: transformation1, transformation2

!   if coefficients is present, it is passed through to the element subroutine
    type(coefficients_t), intent(in), optional :: coefficients

!   if oldvectors is present, it is passed through to the element subroutine.
!   This structure can be used to supply old sysvectors/vectors from possibly
!   other problems to the element subroutines. Note that the problem structures
!   of the other problems must be included in oldvectors.
    type(oldvectors_t), intent(in), optional :: oldvectors

!   This is the nodal subroutine that must be supplied by the calling routine
!   for building the matrix node for node.
!   Only needed if typetransformation=1 (distributed transformation) and
!   normalvector=0.
    optional nodesub
    interface
      subroutine nodesub ( mesh, problem, geom, node, first, last, &
        coefficients, oldvectors, nodemat )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: geom, node
        logical, intent(in) :: first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: nodemat
      end subroutine nodesub
    end interface


    logical :: first, last, allnodes
    integer :: trans1, trans2, physq, layer, geom, node
    integer :: i, j, trans, ndof, nnodes, nodenr, k, step, exclude, numec
    integer :: crvnum, srfnum, numes, dof, fac
!   work: indicate nodes that need to be excluded from the transformation
    logical, dimension(mesh%nnodes) :: work
!   work1: temporary to indicate degrees that are transformed
    logical, dimension(:), allocatable :: work1
    real(dp) :: lv
    real(dp), allocatable, dimension(:) :: u1, u2, u3
!   nodes: array with nodes to be considered for transformation
    integer, allocatable, dimension(:) :: nodes, pos
    real(dp), allocatable, dimension(:,:) :: nodemat, normal_in_nodes
    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl

    call check ( mesh, 'build_transformation_matrix' )
    call check ( problem, 'build_transformation_matrix', mesh )

!   initialize

    if ( present(oldvectors) ) oldvl = oldvectors
    if ( present(coefficients) ) coeffl = coefficients

!   which transformations?

    if ( present(transformation1) .and. present(transformation2) ) then
!     specified range of transformations only
      trans1 = transformation1
      trans2 = transformation2
    else if ( present(transformation1) ) then
!     one transformations only
      trans1 = transformation1
      trans2 = transformation1
    else
!     all transformations
      trans1 = 1
      trans2 = problem%numtransformations
    end if

    if ( trans1 < 1 .or. trans1 > problem%numtransformations .or.  &
         trans2 < 1 .or. trans2 > problem%numtransformations ) then

      write(*,'(/2(a/),2(a,i0),/a,i0/)') &
        'Error: transformation1 and/or transformation2 in the heading of', &
        ' build_transformation_matrix is out of range ', &
        ' transformation1 is ', trans1, '; transformation2 is ', trans2, &
        ' whereas the number of transformations is ', problem%numtransformations
      stop

    end if

    if ( any( problem%transformations(trans1:trans2)%typetransformation == 1 &
           .and. problem%transformations(trans1:trans2)%normalvector == 0 ) &
           .and. .not. present(nodesub) ) then
      write(*,'(4(/a)/)') &
        'Error in build_transformation_matrix: ', &
        ' nodesub must be present in the heading, since at least one of ', &
        ' the transformations has typetransformation = 1 (distributed) ', &
        ' and normalvector = 0 (user transformation matrix per node). '
      stop
    end if

!   work array to store transformed degrees of freedom

    allocate ( work1(problem%numdegfd) )

    work1 = .false.
    work1(problem%degfdtrans) = .true.

!   Fill diagonal with 1 on rows of non-transformed degrees

    do i = 1, problem%numdegfd
      if ( work1(i) ) cycle
!     degree not transformed (value of 1 on the diagonal)
      j = problem%Amat%ia(i)
      problem%Amat%a(j) = 1
      problem%Amat%ja(j) = i
    end do

    deallocate(work1)

    work = .false.

!   start loop over all transformations

    do trans = trans1, trans2

!     inactive transformation?
      if ( problem%transformations(trans)%typetransformation == 0 ) cycle

      physq = problem%transformations(trans)%physq
      layer = problem%transformations(trans)%layer
      geom = problem%transformations(trans)%geometry

!     exclude nodes

      call exclude_nodes ( .true. )

!     get nodes

      call get_nodes

      allnodes = all ( .not. work(nodes) )

!     max number of degrees in the nodes

      call ndof_nodes

      allocate ( pos(ndof) )

!     normal vector

      if ( any ( problem%transformations(trans)%normalvector == [-1,1] ) ) then
        call compute_normalvector_in_nodes
      end if

!     loop over nodes

      do node = 1, nnodes

        if ( allnodes ) then
          first = node == 1
          last = node == nnodes
        else
          first = .true.
          last = .true.
        end if

        nodenr = nodes(node)

        if ( work(nodenr) ) cycle  ! node excluded

!       get degrees positions

        if ( physq > 0 ) then
!         physical quantity specified
          call pos_array_node ( problem, nodenr, dof, pos, &
            physqarr=[physq], layer=layer )
        else
          call pos_array_node ( problem, nodenr, dof, pos, layer=layer )
        end if

        allocate ( nodemat(dof,dof) )

        select case ( problem%transformations(trans)%typetransformation )
        case(1) ! distributed
          if ( problem%transformations(trans)%normalvector == 0 ) then
!           nodal subroutine
            call nodesub ( mesh, problem, geom, node, first, last, &
              coeffl, oldvl, nodemat )
          else
!           normal direction is the new x'-direction
            fac = sign ( 1, problem%transformations(trans)%normalvector )
            select case ( mesh%ndim )
            case(2) ! 2D
              allocate ( u1(2), u2(2) )
              lv = sqrt ( dot_product ( normal_in_nodes(node,:), &
                                       &normal_in_nodes(node,:) ) )
              u1 = fac * normal_in_nodes(node,:) / lv
              u2 = [ -u1(2), u1(1) ]
              nodemat(:,1) = u1
              nodemat(:,2) = u2
              deallocate ( u1, u2 )
            case(3) ! 3D
              allocate ( u1(3), u2(3), u3(3) )
              lv = sqrt ( dot_product ( normal_in_nodes(node,:), &
                                       &normal_in_nodes(node,:) ) )
              u1 = fac * normal_in_nodes(node,:) / lv
              u2 = cross_product ( u1, problem%transformations(trans)%v2 )
              lv = sqrt ( dot_product ( u2, u2 ) )
              u2 = u2 / lv
              u3 = cross_product ( u1, u2 )
              nodemat(:,1) = u1
              nodemat(:,2) = u2
              nodemat(:,3) = u3
              deallocate ( u1, u2, u3 )
            case default
              call errormsg_case_default ( 'build_transformation_matrix', &
                'problem%transformations(trans)%typetransformation', &
                int_value=problem%transformations(trans)%typetransformation )
            end select
          end if
        case(2) ! global
          nodemat = problem%transformations(trans)%Amat_global
        case default
          call errormsg_case_default ( 'build_transformation_matrix', &
            'mesh%ndim', int_value=mesh%ndim )
        end select

!       insert matrix in Amat

        call insert_nodemat_in_Amat

        deallocate ( nodemat )

        problem%transformations(trans)%build = .true.

      end do

!     set back exclude nodes

      call exclude_nodes ( .false. )

      deallocate ( nodes, pos )

      if ( any ( problem%transformations(trans)%normalvector == [-1,1] ) ) then
        deallocate ( normal_in_nodes )
      end if

    end do

!   set buildAmat

    if ( all( problem%transformations(1:problem%numtransformations)%build &
              .or. problem%transformations(1:problem%numtransformations)%&
               &typetransformation == 0 ) ) then
      problem%buildAmat = .true.
    else
      problem%buildAmat = .false.
    end if

  contains


!   insert nodemat into the sparse matrix Amat

    subroutine insert_nodemat_in_Amat

      integer :: r, c

      do r = 1, dof
        i = pos(r)  ! row number
        k = problem%Amat%ia(i) - 1
        do c = 1, dof
          j = pos(c)  ! column number
          k = k + 1  ! position in ja and a array
          problem%Amat%ja(k) = j
          problem%Amat%a(k) = nodemat(r,c)
        end do
      end do

    end subroutine insert_nodemat_in_Amat

!  set or unset nodes that need to be excluded from transformations

    subroutine exclude_nodes ( val )

      logical, intent(in) :: val

      integer :: curve, surface

      if ( problem%transformations(trans)%typegeometry == 2 ) then
!       curve
        nnodes = mesh%curves(geom)%nnodes
        step = problem%transformations(trans)%step
        if ( step < 0 ) then
!         exclude nodes with increment -step
          work(mesh%curves(geom)%nodes(1:nnodes:-step)) = val
        end if
        exclude = problem%transformations(trans)%exclude
        if ( exclude == 1 .or. exclude == 3 ) then
!         exclude first node
          work(mesh%curves(geom)%nodes(1)) = val
        end if
        if ( exclude == 2 .or. exclude == 3 ) then
!         exclude last node
          work(mesh%curves(geom)%nodes(nnodes)) = val
        end if
      end if
      work ( &
        mesh%points ( problem%transformations(trans)%excludepoints ) ) = val
      numec = size ( problem%transformations(trans)%excludecurves )
      do curve = 1, numec
        crvnum = problem%transformations(trans)%excludecurves(curve)
        work ( mesh%curves(crvnum)%nodes ) = val
      end do
      numes = size ( problem%transformations(trans)%excludesurfaces )
!     set work
      do surface = 1, numes
        srfnum = problem%transformations(trans)%excludesurfaces(surface)
        work ( mesh%surfaces(srfnum)%nodes ) = val
      end do

    end subroutine exclude_nodes


!   determine number of degrees in the nodes (max)

    subroutine ndof_nodes

      integer :: node, nodenr

      ndof = 0

      if ( physq > 0 ) then
        do node = 1, nnodes
          nodenr = nodes(node)
          ndof = max( ndof, problem%vec_nodnumdegfd(nodenr+1,&
                                    &problem%physq(physq)) &
                      - problem%vec_nodnumdegfd(nodenr,&
                                    &problem%physq(physq)) )
        end do
      else
        do node = 1, nnodes
          nodenr = nodes(node)
          ndof = max( ndof, problem%nodnumdegfd(nodenr+1) &
                      - problem%nodnumdegfd(nodenr) )
        end do
      end if

    end subroutine ndof_nodes


!   get the nodes for the transformation

    subroutine get_nodes

      select case ( problem%transformations(trans)%typegeometry )
      case(1) ! transformation in point
        nnodes = 1
        allocate(nodes(1))
        nodes = [ mesh%points(geom) ]
      case(2)
!       transformation on curve
        nnodes = mesh%curves(geom)%nnodes
        allocate(nodes(nnodes))
        nodes = mesh%curves(geom)%nodes
      case(3)
!       transformation on surface
        nnodes = mesh%surfaces(geom)%nnodes
        allocate(nodes(nnodes))
        nodes = mesh%surfaces(geom)%nodes
      case(4)
!       transformation on volume
        nnodes = mesh%volumes(geom)%nnodes
        allocate(nodes(nnodes))
        nodes = mesh%volumes(geom)%nodes
      case(5)
!       transformation on nodeset
        nnodes = size(mesh%nodesets(geom)%a)
        allocate(nodes(nnodes))
        nodes = mesh%nodesets(geom)%a
      case default
        call errormsg_case_default ( 'get_nodes', &
          'problem%transformations(trans)%typegeometry', &
          int_value=problem%transformations(trans)%typegeometry )
      end select

    end subroutine get_nodes


!   Compute the normalvector in the nodes of the geometry
!   Common nodes are added, basically averaging the contribution from elements

    subroutine compute_normalvector_in_nodes

      use shapefunc_m

      integer :: elem, ndim, nodalp, ndimr
      real(dp), dimension(:), allocatable :: surfl
      real(dp), dimension(:,:), allocatable :: xig, x, phi, normal, xrnod
      real(dp), dimension(:,:,:), allocatable :: dphi, dxdxis
      type(geometry_t) :: geometry


      select case ( problem%transformations(trans)%typegeometry )
      case(2)
!       normal vector on curve
        geometry = mesh%curves(geom)
        ndimr = 1
        ndim = mesh%curves(geom)%ndim
        nodalp = mesh%curves(geom)%elnumnod
      case(3)
!       normal vector on surface
        geometry = mesh%surfaces(geom)
        ndimr = 2
        ndim = mesh%surfaces(geom)%ndim
        nodalp = mesh%surfaces(geom)%elnumnod
      case default
        call errormsg_case_default ( 'compute_normalvector_in_nodes', &
          'problem%transformations(trans)%typegeometry', &
          int_value=problem%transformations(trans)%typegeometry )
      end select

      allocate ( normal_in_nodes(geometry%nnodes,ndim) )

!     allocate local arrays

      allocate ( surfl(nodalp), xig(nodalp,ndimr), x(nodalp,ndim) )
      allocate ( phi(nodalp,nodalp), dphi(nodalp,nodalp,ndimr) )
      allocate ( dxdxis(nodalp,ndim,ndimr) )
      allocate ( normal(nodalp,ndim) )
      allocate ( xrnod(nodalp,ndimr) )

!     set shape function of elements

      select case ( geometry%element%elshape )
        case(1) ! two-node line
          xrnod = reshape ( [ -1, 1 ], [2,1] )
          call shape_line_P1_2 ( xrnod, phi, dphi )
        case(2) ! three-node line
          xrnod = reshape ( [ -1, 0, 1 ], [3,1] )
          call shape_line_P2_2 ( xrnod, phi, dphi )
        case(3) ! three-node triangle
          xrnod = reshape ( [ 0, 1, 0, &
                              0, 0, 1 ], [3,2] )
          call shape_triangle_P1 ( xrnod, phi, dphi )
        case(4) ! six-node triangle
          xrnod = &
          reshape ( [ 0._dp, 0.5_dp, 1._dp, 0.5_dp, 0._dp, 0._dp, &
                      0._dp, 0._dp,  0._dp, 0.5_dp, 1._dp, 0.5_dp ], [6,2] )
          call shape_triangle_P2 ( xrnod, phi, dphi )
        case(5) ! four-node quadrilateral
          xrnod = reshape ( [ -1,  1, 1, -1, &
                              -1, -1, 1,  1 ], [4,2] )
          call shape_quad_Q1 ( xrnod, phi, dphi )
        case(6) ! nine-node quadrilateral
          xrnod = reshape ( [ -1,  0,  1, 1, 1, 0, -1, -1, 0, &
                              -1, -1, -1, 0, 1, 1,  1,  0, 0 ], [9,2] )
          call shape_quad_Q2 ( xrnod, phi, dphi )
        case(7) ! seven-node triangle
          xrnod = &
          reshape ( [ 0._dp, 0.5_dp, 1._dp, 0.5_dp, 0._dp,  0._dp, 1._dp/3, &
                      0._dp,  0._dp, 0._dp, 0.5_dp, 1._dp, 0.5_dp, 1._dp/3 ],&
                      [7,2] )
          call shape_triangle_P2plus ( xrnod, phi, dphi )
        case(9) ! five-node quadrilateral
          xrnod = reshape ( [ -1,  1, 1, -1, 0, &
                              -1, -1, 1,  1, 0 ], [5,2] )
          call shape_quad_Q1plus ( xrnod, phi, dphi )
        case(10) ! four-node triangle
          xrnod = reshape ( [ 0._dp,  1._dp, 0._dp, 1._dp/3, &
                              0._dp,  0._dp, 1._dp, 1._dp/3 ], [4,2] )
          call shape_triangle_P1plus ( xrnod, phi, dphi )
        case(30) ! 8-node quadrilateral
          xrnod = reshape ( [ -1,  0,  1, 1, 1, 0, -1, -1, &
                              -1, -1, -1, 0, 1, 1,  1,  0 ], [8,2] )
          call shape_quad_serendipity2 ( xrnod, phi, dphi )
        case default
          write(*,'(/2a,i0/)') 'Error compute_normal_vector_in_nodes: ', &
            ' invalid element shape: ', geometry%element%elshape
          stop
      end select


      normal_in_nodes = 0

!     loop over elements

      do elem = 1, geometry%nelem

        call get_coordinates_geometry ( mesh, elem, x, ndimr=ndimr, &
          geometry=geom )

        call isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, normal )

!       add element normal vector

        normal_in_nodes(geometry%topology(:,elem,1),:) = &
                 normal_in_nodes(geometry%topology(:,elem,1),:) + normal

      end do

!     deallocate local arrays

      deallocate ( surfl, xig, x )
      deallocate ( phi, dphi )
      deallocate ( dxdxis )
      deallocate ( normal )
      deallocate ( xrnod )

    end subroutine compute_normalvector_in_nodes

  end subroutine build_transformation_matrix


! clear some rows in sysmatrix en entries in rhsd

  subroutine clear_rows_sysmatrix ( sysmatrix, rows, essentialpart )

!   the system matrix
!   on output the specfied rows are set to zero
    type(sysmatrix_t), intent(inout) :: sysmatrix

!   array of the degrees of freedom numbers for rows to be set to zero
    integer, dimension(:), intent(in) :: rows

!   if essentialpart=.true. the rows of the sysmatrix with respect to the
!   prescribed part is set to zero (cleared).
!   default=.false.
    logical, intent(in), optional :: essentialpart

    logical :: lessentialpart

    lessentialpart = set_optional ( variable=essentialpart, default=.false. )

    if ( .not. sysmatrix%initialized_structure ) then
      write(*,'(/a/)') &
        'Error in clear_rows_sysmatrix: no system matrix structure.'
      stop
    end if

    if ( .not. sysmatrix%finalized ) then
      write(*,'(/2(a/))') &
        'Error in clear_rows_sysmatrix:, ', &
        '  system matrix has not been finalized.'
      stop
    end if

    if ( .not. sysmatrix%allocated_data ) then
      write(*,'(/a/)') &
        'Error in clear_rows_sysmatrix:, ', &
        '  data in system matrix not allocated.'
      stop
    end if

    if ( sysmatrix%symmetric ) then

      write(*,'(3(a/))') &
        'Error in clear_rows_sysmatrix:', &
        ' clear of rows not available for a ', &
        ' symmetric system matrix structure.'
      stop

    end if

!   clear rows in matrix (unknown part only)
!   note: rows outside the range are ignored

    call clear_rows_sparsematrix ( rows, sysmatrix%Suu )
    call clear_rows_sparsematrix ( rows, sysmatrix%Sup )

    if ( lessentialpart ) then

!     clear rows in matrix (essential part only)

      call clear_rows_sparsematrix ( rows, sysmatrix%Spu )
      call clear_rows_sparsematrix ( rows, sysmatrix%Spp )

    end if

  end subroutine clear_rows_sysmatrix


! Add sparse matrix to one of the sparse matrices of the system matrix

  subroutine add_smat_to_sysmatrix ( sysmatrix, A, part )

    type ( sysmatrix_t ), intent(inout) :: sysmatrix

    type ( sparsematrix_t ), intent(in) :: A

!   indicate which part of the system matrix: uu, up, pu, pp
    character(len=*), intent(in) :: part

    select case ( part )
    case ( 'uu' )
      call smat_add_sys ( sysmatrix%Suu, A, sysmatrix%symmetric )
    case ( 'up' )
      call smat_add_sys ( sysmatrix%Sup, A )
    case ( 'pu' )
      call smat_add_sys ( sysmatrix%Spu, A )
    case ( 'pp' )
      call smat_add_sys ( sysmatrix%Spp, A )
    case default
      write(*,'(/2a/)') &
        'Error add_smat_to_sysmatrix: part must be one of ', &
        '''uu'', ''up'', ''pu'', ''pp''.'
      stop
    end select

  end subroutine add_smat_to_sysmatrix


! Add sparse matrix B to a sparse matrix A (part of the system matrix)

  subroutine smat_add_sys ( A, B, symmetric )

    type ( sparsematrix_t ), intent(inout) :: A

    type ( sparsematrix_t ), intent(in) :: B

!   the matrix A is symmetric and stored symmetric (only upper triangle)
!   default = .false.
    logical, intent(in), optional :: symmetric


    logical :: lsymmetric
    integer :: p, pa, pb, row, col, jcol
    integer, dimension(:), allocatable :: w


    lsymmetric = set_optional ( variable=symmetric, default=.false. )


    if ( any ( [ A%n, A%m ] /= [ B%n, B% m ] ) ) then

      write(*,'(/a/a/)') &
        'Error smat_add_sys: Sparse matrices A and B have different ', &
        '  row and/or column dimension.'
      stop

    end if

    allocate ( w(A%m) )

    w = 0

!   add sparse matrix B to A row by row

    do row = 1, B%n

!     fill w for non-zeros already in matrix A

!     p is position in the row of A for new data

      p = A%ia(row+1)

      do pa = A%ia(row), A%ia(row+1) - 1

        jcol = A%ja(pa)

        if ( jcol == 0 ) then
!         matrix row not completely filled
          p = pa
          exit
        end if

        w(jcol) = pa  ! put in w: position in a of non-zero entry of A

      end do

!     add element row of B to A

      do pb = B%ia(row), B%ia(row+1) - 1

        col = B%ja(pb)

        if ( lsymmetric ) then
!         symmetric matrix
          if ( row > col ) cycle  ! ignore lower triangle
        end if

        pa = w(col)

        if ( pa == 0 .and. p /= A%ia(row+1) ) then

!         put new non-zero in matrix

          A%a(p) = B%a(pb)
          A%ja(p) = col

          p = p + 1

        else if ( pa /= 0 ) then

!         add to existing non-zero

          A%a(pa) = A%a(pa) + B%a(pb)

        else if ( pa == 0 .and. p == A%ia(row+1) ) then

!         new non-zero but no free space

          write(*,'(/a,i0/a/)') &
            'Internal error smat_add_sys: no more space in row ', row, &
            ' of sparse matrix A to add sparse matrix B.'
          stop

        end if

      end do

!     put w back to zero

      do pa = A%ia(row), A%ia(row+1) - 1

        jcol = A%ja(pa)

        if ( jcol == 0 ) exit

        w(jcol) = 0

      end do

    end do

    deallocate ( w )

  end subroutine smat_add_sys


! Clear sysmatrix

  subroutine clear_sysmatrix ( sysmatrix )

    type(sysmatrix_t), intent(inout) :: sysmatrix

    sysmatrix%Suu%a = 0
    sysmatrix%Suu%ja = 0
    sysmatrix%Sup%a = 0
    sysmatrix%Sup%ja = 0
    sysmatrix%Spu%a = 0
    sysmatrix%Spu%ja = 0
    sysmatrix%Spp%a = 0
    sysmatrix%Spp%ja = 0

  end subroutine clear_sysmatrix


! Delete single sysmatrix

  subroutine delete_single_sysmatrix ( sysmatrix, nr )

    type(sysmatrix_t), intent(inout) :: sysmatrix
    integer, intent(in) :: nr

    if ( .not. sysmatrix%initialized_structure .or. &
         .not. sysmatrix%allocated_data ) then
      write(*,'(3(/a)/a,i0/)') &
        'Error in delete_sysmatrix: incomplete matrix ', &
        ' system matrix structure and/or data has not been allocated.', &
        ' sysmatrix cannot be deleted with this routine.', &
        ' sysmatrix number in heading = ', nr
      stop
    end if

!   deallocate all allocatables and intialize to defaults
    call deall ( sysmatrix )

  contains

    subroutine deall ( sysmatrix )
      type(sysmatrix_t), intent(out) :: sysmatrix
    end subroutine deall

  end subroutine delete_single_sysmatrix


! Delete system matrix

  subroutine delete_sysmatrix ( sysmatrix1, sysmatrix2, sysmatrix3, &
    sysmatrix4, sysmatrix5 )

    type(sysmatrix_t), intent(inout) :: sysmatrix1
    type(sysmatrix_t), intent(inout), optional :: sysmatrix2, sysmatrix3, &
      sysmatrix4, sysmatrix5

    call delete_single_sysmatrix(sysmatrix1,1)
    if ( present(sysmatrix2) ) call delete_single_sysmatrix(sysmatrix2,2)
    if ( present(sysmatrix3) ) call delete_single_sysmatrix(sysmatrix3,3)
    if ( present(sysmatrix4) ) call delete_single_sysmatrix(sysmatrix4,4)
    if ( present(sysmatrix5) ) call delete_single_sysmatrix(sysmatrix5,5)

  end subroutine delete_sysmatrix

end module system_matrix_m
