
! Copyright (C) 2023-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Routines, for the system of equations involving dependencies

module system_dependency_m

  use glob_defs_m
  use mesh_m
  use problem_defs_m
  use sparse_m
  use system_defs_m
  use system_matrix_m
  use element_defs_m
  use set_optional_m

  implicit none


! Definition of the matrices Eu, Ep and vector g in
!
!   u_p = [ Eu, Ep ] ( u_u ) + g
!                    ( u_p )
!
! for updating the system matrix to take dependencies into account
!
  type sysdep_t

!   the structure of the matrices (the arrays ia) have been initialized
    logical :: initialized_structure = .false.

!   the structure of the matrices (the arrays ia) have been finalized
    logical :: finalized = .false.

!   the data in the matrices and vector has been allocated
    logical :: allocated_data = .false.

!   dependency sparse matrices, size: Eu(np,nu), Ep(np,np)
    type(sparsematrix_t) :: Eu, Ep

!   dependency vector g, size: g(np)
    real(dp), allocatable, dimension(:) :: g

  end type sysdep_t


! type definition of a vector subscript in the dependency part of a sysvector
! additional unknowns only

  type subscriptdep_t

!   the subscript. Data can be accessed with sysvector%u(subscript%s).
    integer, allocatable, dimension(:) :: s

  end type subscriptdep_t


! interface for generic create subroutine

  interface create
    module procedure create_subscript_dependency
  end interface create

  interface create_subscript
    module procedure create_subscript_dependency
  end interface create_subscript


! interface for generic delete subroutine

  interface delete
    module procedure delete_subscript_dependency
  end interface delete

contains


! Create system matrix structure of dependencies (ia not yet accumulated)

  subroutine create_sysmatrix_structure_dependency ( sysmatrix, sysdep, &
    mesh, problem, symmetric )

    type(sysmatrix_t), intent(inout) :: sysmatrix
    type(sysdep_t),    intent(out) :: sysdep
    type(mesh_t),      intent(in)  :: mesh
    type(problem_t),   intent(in)  :: problem

!   if symmetric is present and has the values .true. the matrix is assumed to
!   be symmetric and only the non-zero elements of the upper-triangle
!   (including the diagonal) of the matrix of unknowns (Suu) are stored.
!   NOTE: this has no effect if sysmatrix has been initialized already.
!   default: symmetric=.false.
    logical, intent(in), optional :: symmetric

!   This routine fills the arrays ia with the number of non-zeros in each row:
!
!       ia(i+1) = number of non-zeros in row i
!
!   ia is not yet accumulated and may still be adapted (add more non-zeros)
!
!   Here space is created in the system matrix for the adding the dependencies
!   to the matrix. Furthermore, space is created in the matrices Eu, Ep
!   contained in sysdep for storing the dependency relations.

    integer :: numund, numess, dep, nnodes1, nnodes2, geom
    integer :: elemset, exclude, step, node, crvnum, srfnum
    integer :: nodenr1, nodenr2, physq, layer, dofr, dof, n1, n2
    integer :: numnod1, numnod2, deg_r1, deg_r2, nadd
    integer, dimension(mesh%maxnodnumnod+1) :: nodes_g1, nodes_g2
    integer, dimension(problem%maxnoddegfd) :: pos, posr, posc
    integer, allocatable, dimension(:) :: utype
    integer, dimension(:), allocatable :: nodes1, nodes2

    allocate ( utype(problem%numdegfd) )

    call check ( mesh, 'create_sysmatrix_structure_dependency' )
    call check ( problem, 'create_sysmatrix_structure_dependency', mesh )

    if ( sysmatrix%finalized ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysmatrix_structure_dependency: ', &
        ' sysmatrix has already been finalized '
      stop
    end if

    if ( any( &
      problem%dependencies(1:problem%numdependencies)%typedependency /= 0 ) &
       ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysmatrix_structure_dependency: ', &
        ' typedependency /= 0 not yet supported. '
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

!   initialize?

    if ( .not. sysmatrix%initialized_structure ) then

!     initialize system_matrix

      call initialize_sysmatrix_structure ( sysmatrix, numund, numess )

!     set symmetry of the matrix

      sysmatrix%symmetric = set_optional ( variable=symmetric, default=.false. )

    end if

!   initialize sysdep

    call initialize_sysdep_structure ( sysdep, numund, numess )

!   loop over all dependencies

    do dep = 1, problem%numdependencies

!     dependency on what geometrical object?

      if ( problem%dependencies(dep)%geometry1 > 0 ) then
!       dependency on geometry
        geom = problem%dependencies(dep)%geometry1
        if ( problem%dependencies(dep)%typegeometry1 == 1 ) then
!         dependency in point
          nnodes1 = 1
          nodes1 = [mesh%points(geom)]
        else if ( problem%dependencies(dep)%typegeometry1 == 2 ) then
!         dependency on curve
          nnodes1 = mesh%curves(geom)%nnodes
          nodes1 = mesh%curves(geom)%nodes
        else if ( problem%dependencies(dep)%typegeometry1 == 3 ) then
!         dependency on surface
          nnodes1 = mesh%surfaces(geom)%nnodes
          nodes1 = mesh%surfaces(geom)%nodes
        else if ( problem%dependencies(dep)%typegeometry1 == 4 ) then
!         dependency on volume
          nnodes1 = mesh%volumes(geom)%nnodes
          nodes1 = mesh%volumes(geom)%nodes
        else if ( problem%dependencies(dep)%typegeometry1 == 5 ) then
!         dependency on nodeset
          nodes1 = size(mesh%nodesets(geom)%a)
          nodes1 = mesh%nodesets(geom)%a
        end if
      else if ( problem%dependencies(dep)%elementset1 > 0 ) then
!       dependency on elementset
        elemset = problem%dependencies(dep)%elementset1
        nnodes1 = mesh%elementsets(elemset)%nnodes
        nodes1 = mesh%elementsets(elemset)%nodes
      end if

!     second geometrical object

      nnodes2 = 0

      if ( problem%dependencies(dep)%geometry2 > 0 ) then
!       dependency on geometry
        geom = problem%dependencies(dep)%geometry2
        if ( problem%dependencies(dep)%typegeometry2 == 1 ) then
!         dependency in point
          nnodes2 = 1
          nodes2 = [mesh%points(geom)]
        else if ( problem%dependencies(dep)%typegeometry2 == 2 ) then
!         dependency on curve
          nnodes2 = mesh%curves(geom)%nnodes
          nodes2 = mesh%curves(geom)%nodes
        else if ( problem%dependencies(dep)%typegeometry2 == 3 ) then
!         dependency on surface
          nnodes2 = mesh%surfaces(geom)%nnodes
          nodes2 = mesh%surfaces(geom)%nodes
        else if ( problem%dependencies(dep)%typegeometry2 == 4 ) then
!         dependency on volume
          nnodes2 = mesh%volumes(geom)%nnodes
          nodes2 = mesh%volumes(geom)%nodes
        else if ( problem%dependencies(dep)%typegeometry2 == 5 ) then
!         dependency on nodeset
          nodes2 = size(mesh%nodesets(geom)%a)
          nodes2 = mesh%nodesets(geom)%a
        end if
      else if ( problem%dependencies(dep)%elementset2 > 0 ) then
!       dependency on elementset
        elemset = problem%dependencies(dep)%elementset2
        nnodes2 = mesh%elementsets(elemset)%nnodes
        nodes2 = mesh%elementsets(elemset)%nodes
      end if

      if ( nnodes2 > 0 ) then
        if ( nnodes1 /= nnodes2 ) then
          write(*,'(2(/a)/)') &
            'Error in create_sysmatrix_structure_dependency: ', &
            ' number of nodes of second part /= dependent (first) part.'
          stop
        end if
      end if

!     exclude nodes

      call exclude_nodes

!     range of additional unknowns degrees for current value of dep

      deg_r1 = problem%dependencies(dep)%addnumdegfd(1) + 1
      deg_r2 = problem%dependencies(dep)%addnumdegfd(2)
      nadd = problem%dependencies(dep)%naddunknowns

!     scan all nodes for degrees of freedom and look for connections

      do node = 1, nnodes1

        if ( nodes1(node) == 0 ) cycle ! node has been excluded

!       nodes on part 1 (only connected nodes)
        nodenr1 = nodes1(node)
        n1 = mesh%nodnumnod(nodenr1)
        n2 = mesh%nodnumnod(nodenr1+1)
        numnod1 = n2 - n1
        nodes_g1(1:numnod1) = mesh%nodnod(n1+1:n2)

!       count matrix entries of diagonal block of additional unknowns

        call count_matrix_entries_addunknowns_diagonal

!       count matrix entries of additional unknowns connecting nodes from part 1

        call count_matrix_entries_addunknowns ( part=1 )

!       count E matrix entries of additional unknowns

        call count_Ematrix_entries_addunknowns

!       part 2 present?

        if ( nnodes2 > 0 ) then

!         part 2 is present

!         nodes on part 2 (current node + connections)
          nodenr2 = nodes2(node)
          nodes_g2(1) = nodenr2
          n1 = mesh%nodnumnod(nodenr2)
          n2 = mesh%nodnumnod(nodenr2+1)
          numnod2 = n2 - n1 + 1
          nodes_g2(2:numnod2) = mesh%nodnod(n1+1:n2)

!         count matrix entries of connecting fem nodes from part 1

          call count_matrix_entries_fem_nodes

!         count matrix entries of connecting nodes from part 2 (addunknowns)

          call count_matrix_entries_addunknowns ( part=2 )

!         count E matrix entries of remaining degrees on part 2

          call count_Ematrix_entries_degrees_from_part2

        end if

      end do

    end do

!   finalized sysdep

    call finalize_sysdep_structure ( sysdep )

  contains


!   exclude nodes from the dependency

    subroutine exclude_nodes

      integer :: node, curve, surface
      logical, allocatable, dimension(:) :: work, work1

!     exclude nodes

      allocate( work(nnodes1) )

      if ( problem%dependencies(dep)%step < 0 ) then
!       exclude nodes with increment -step
        work = .true.
        step = - problem%dependencies(dep)%step
        do node = 1, nnodes1, step
          work(node) = .false.
        end do
      else if ( problem%dependencies(dep)%step > 0 ) then
!       include nodes with increment step
        work = .false.
        step = problem%dependencies(dep)%step
        do node = 1, nnodes1, step
          work(node) = .true.
        end do
      else
!       all nodes
        work = .true.
      end if

      exclude = problem%dependencies(dep)%exclude
      if ( exclude == 1 .or. exclude == 3 ) then
!       exclude first node
        work(1) = .false.
      end if
      if ( exclude == 2 .or. exclude == 3 ) then
!       exclude last node
        work(nnodes1) = .false.
      end if

      if ( size(problem%dependencies(dep)%excludepoints) > 0 .or. &
           size(problem%dependencies(dep)%excludecurves) > 0 .or. &
           size(problem%dependencies(dep)%excludesurfaces) > 0 ) then
!       exclude nodes in points, on curves and/or surfaces
        allocate ( work1(mesh%nnodes) )
        work1 = .false.
!       set work1
        work1 ( &
            mesh%points ( problem%dependencies(dep)%excludepoints ) ) = .true.
!       set work1
        do curve = 1, size ( problem%dependencies(dep)%excludecurves )
          crvnum = problem%dependencies(dep)%excludecurves(curve)
          work1 ( mesh%curves(crvnum)%nodes ) = .true.
        end do
!       set work1
        do surface = 1, size ( problem%dependencies(dep)%excludesurfaces )
          srfnum = problem%dependencies(dep)%excludesurfaces(surface)
          work1 ( mesh%surfaces(srfnum)%nodes ) = .true.
        end do
!       delete nodes
        do node = 1, nnodes1
          if ( work1( nodes1(node) ) ) then
!           exclude node
            work(node) = .false.
          end if
        end do
        deallocate ( work1 )
      end if

!     set node number to 0 for excluded nodes

      where ( .not. work )
        nodes1 = 0
        nodes2 = 0
      end where

      deallocate(work)

    end subroutine exclude_nodes


!   count matrix entries of diagonal block of additonal unknowns

    subroutine count_matrix_entries_addunknowns_diagonal

      integer :: deg_r, deg_c, row, col

      do deg_r = deg_r1, deg_r2

        row = problem%degfdperm(deg_r,2)

        if ( sysmatrix%symmetric ) then

!         symmetric matrix

          do deg_c = deg_r1, deg_r2

            col = problem%degfdperm(deg_c,2)

            if ( col < row ) cycle  ! element in lower diagonal

            sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1

          end do

        else

          sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + nadd

        end if

      end do

    end subroutine count_matrix_entries_addunknowns_diagonal


!   count matrix entries of node connections of nodes on part 1

    subroutine count_matrix_entries_fem_nodes

      integer :: dr, dc, i, nd, deg_r, deg_c, row, col, colp, rowp, dofc

!     row degrees from connecting node on part 2
      physq = problem%dependencies(dep)%physq2
      layer = problem%dependencies(dep)%layer2

      if ( physq > 0 ) then
        call pos_array_local_node ( problem, nodenr2, dof, pos, [physq], &
          layer )
      else
        call pos_array_local_node ( problem, nodenr2, dof, pos, layer=layer )
      end if

!     store and change to absolute positions

      dofr = dof
      posr(1:dof) = problem%nodnumdegfd(nodenr2) + pos(1:dof)

      do dr = 1, dofr

        deg_r = posr(dr)

        row = problem%degfdperm(deg_r,2)

!       count matrix entries of node connections on part 1

        do i = 1, numnod1

          nd = nodes_g1(i)

          physq = problem%dependencies(dep)%physq1
          layer = problem%dependencies(dep)%layer1

          if ( physq > 0 ) then
            call pos_array_local_node ( problem, nd, dof, pos, [physq], &
              layer )
          else
            call pos_array_local_node ( problem, nd, dof, pos, layer=layer )
          end if

!         column degrees

          dofc = dof
          posc(1:dof) = problem%nodnumdegfd(nd) + pos(1:dof)

!         start loop over degrees

          do dc = 1, dofc

            deg_c = posc(dc)

            col = problem%degfdperm(deg_c,2)

            if ( utype(deg_r) == 1 ) then

!             unknown degree of freedom (row)

              if ( utype(deg_c) == 1 ) then

!               unknown degree of freedom (column)

                if ( sysmatrix%symmetric ) then
!                 symmetric matrix
                  if ( col > row ) then
!                   element in upper triangle
                    sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1
                  else
!                   element in lower triangle, add transposed element
                    sysmatrix%Suu%ia(col+1) = sysmatrix%Suu%ia(col+1) + 1
                  end if
                else
!                 unsymmetric matrix
                  sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1
                  sysmatrix%Suu%ia(col+1) = sysmatrix%Suu%ia(col+1) + 1
                end if

              else if ( utype(deg_c) == 3 ) then

!               prescribed degree of freedom (column)
                sysmatrix%Sup%ia(row+1) = sysmatrix%Sup%ia(row+1) + 1
                colp = col - numund
                sysmatrix%Spu%ia(colp+1) = sysmatrix%Spu%ia(colp+1) + 1

              end if

            else if ( utype(deg_r) == 3 ) then

!             prescribed degree of freedom (row)

              rowp = row - numund

              if ( utype(deg_c) == 1 ) then

!               unknown degree of freedom (column)
                sysmatrix%Sup%ia(col+1) = sysmatrix%Sup%ia(col+1) + 1
                sysmatrix%Spu%ia(rowp+1) = sysmatrix%Spu%ia(rowp+1) + 1

              else if ( utype(deg_c) == 3 ) then

!               prescribed degree of freedom (column)
                sysmatrix%Spp%ia(rowp+1) = sysmatrix%Spp%ia(rowp+1) + 1
                colp = col - numund
                sysmatrix%Spp%ia(colp+1) = sysmatrix%Spp%ia(colp+1) + 1

              end if

            end if

          end do

        end do

      end do

    end subroutine count_matrix_entries_fem_nodes


!   count matrix entries of node connections of additional unknowns

    subroutine count_matrix_entries_addunknowns ( part )

      integer, intent(in) :: part

      integer :: dc, i, nd, deg_r, deg_c, row, col, colp, dofc, numnod
      integer, dimension(:), allocatable :: nodes_g

!     choose part 1 or 2

      select case ( part )
      case(1)
        numnod = numnod1
        physq = problem%dependencies(dep)%physq1
        layer = problem%dependencies(dep)%layer1
        nodes_g = nodes_g1(1:numnod1)
      case(2)
        numnod = numnod2
        physq = problem%dependencies(dep)%physq2
        layer = problem%dependencies(dep)%layer2
        nodes_g = nodes_g2(1:numnod2)
      case default
        call errormsg_case_default ( 'count_matrix_entries_addunknowns', &
          'part', int_value=part )
      end select

!     loop over row degrees

      do deg_r = deg_r1, deg_r2

        row = problem%degfdperm(deg_r,2)

!       count matrix entries of node connections

        do i = 1, numnod

          nd = nodes_g(i)

          if ( physq > 0 ) then
            call pos_array_local_node ( problem, nd, dof, pos, [physq], &
              layer )
          else
            call pos_array_local_node ( problem, nd, dof, pos, layer=layer )
          end if

!         column degrees

          dofc = dof
          posc(1:dof) = problem%nodnumdegfd(nd) + pos(1:dof)

!         start loop over degrees

          do dc = 1, dofc

            deg_c = posc(dc)

            col = problem%degfdperm(deg_c,2)

            if ( utype(deg_r) == 1 ) then

!             unknown degree of freedom (row)

              if ( utype(deg_c) == 1 ) then

!               unknown degree of freedom (column)

                if ( sysmatrix%symmetric ) then
!                 symmetric matrix
                  if ( col > row ) then
!                   element in upper triangle
                    sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1
                  else
!                   element in lower triangle, add transposed element
                    sysmatrix%Suu%ia(col+1) = sysmatrix%Suu%ia(col+1) + 1
                  end if
                else
!                 unsymmetric matrix
                  sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1
                  sysmatrix%Suu%ia(col+1) = sysmatrix%Suu%ia(col+1) + 1
                end if

              else if ( utype(deg_c) == 3 ) then

!               prescribed degree of freedom (column)
                sysmatrix%Sup%ia(row+1) = sysmatrix%Sup%ia(row+1) + 1
                colp = col - numund
                sysmatrix%Spu%ia(colp+1) = sysmatrix%Spu%ia(colp+1) + 1

              end if

            else if ( utype(deg_r) == 3 ) then

!             prescribed degree of freedom (row)

              stop 'internal error count_matrix_entries_addunknowns'

            end if

          end do

        end do

      end do

    end subroutine count_matrix_entries_addunknowns


!   count E matrix entries of additional unknowns

    subroutine count_Ematrix_entries_addunknowns

      integer :: dr, deg_r, row, rowp

!     row degrees from node on part 1
      physq = problem%dependencies(dep)%physq1
      layer = problem%dependencies(dep)%layer1

      if ( physq > 0 ) then
        call pos_array_local_node ( problem, nodenr1, dof, pos, [physq], &
          layer )
      else
        call pos_array_local_node ( problem, nodenr1, dof, pos, layer=layer )
      end if

!     store and change to absolute positions

      dofr = dof
      posr(1:dof) = problem%nodnumdegfd(nodenr1) + pos(1:dof)

      if ( any( problem%degfdperm(posr(1:dofr),2) <= numund ) ) then
        write(*,'(/a/a,i0,a/)') &
          'Error in create_sysmatrix_structure_dependency: ', &
          ' dependent degrees in node ', nodenr1, ' not defined essential.'
        stop
      end if

      do dr = 1, dofr

        deg_r = posr(dr)

        row = problem%degfdperm(deg_r,2)
        rowp = row - numund

        sysdep%Eu%ia(rowp+1) = sysdep%Eu%ia(rowp+1) + nadd

      end do

    end subroutine count_Ematrix_entries_addunknowns


!   count E matrix entries of remaining degrees on part 2

    subroutine count_Ematrix_entries_degrees_from_part2

      integer :: dr, dc, deg_r, deg_c, row, rowp, dofc

!     row degrees from node on part 1
      physq = problem%dependencies(dep)%physq1
      layer = problem%dependencies(dep)%layer1

      if ( physq > 0 ) then
        call pos_array_local_node ( problem, nodenr1, dof, pos, [physq], &
          layer )
      else
        call pos_array_local_node ( problem, nodenr1, dof, pos, layer=layer )
      end if

!     store and change to absolute positions

      dofr = dof
      posr(1:dof) = problem%nodnumdegfd(nodenr1) + pos(1:dof)

!     col degrees from node on part 2
      physq = problem%dependencies(dep)%physq2
      layer = problem%dependencies(dep)%layer2

      if ( physq > 0 ) then
        call pos_array_local_node ( problem, nodenr2, dof, pos, [physq], &
          layer )
      else
        call pos_array_local_node ( problem, nodenr2, dof, pos, layer=layer )
      end if

      dofc = dof
      posc(1:dof) = problem%nodnumdegfd(nodenr2) + pos(1:dof)

      do dr = 1, dofr

        deg_r = posr(dr)

        row = problem%degfdperm(deg_r,2)
        rowp = row - numund

        do dc = 1, dofc

          deg_c = posc(dc)

          if ( utype(deg_c) == 1 ) then

!           unknown degree of freedom (column)
            sysdep%Eu%ia(rowp+1) = sysdep%Eu%ia(rowp+1) + 1

          else if ( utype(deg_c) == 3 ) then

!           prescribed degree of freedom (column)
            sysdep%Ep%ia(rowp+1) = sysdep%Ep%ia(rowp+1) + 1

          end if

        end do

      end do

    end subroutine count_Ematrix_entries_degrees_from_part2

  end subroutine create_sysmatrix_structure_dependency


! Finalize the structure of the sysdep, i.e. accumulate the ia vector
! and compute nnz

  subroutine finalize_sysdep_structure ( sysdep )

    type(sysdep_t), intent(inout) :: sysdep

!   test

    if ( .not. sysdep%initialized_structure ) then
      write(*,'(2(/a)/)') &
        'Error in finalize_sysdep_structure: ', &
        ' sysdep structure has not been created '
      stop
    end if

    if ( sysdep%finalized ) then
      write(*,'(2(/a)/)') &
        'Error in finalize_sysdep_structure: ', &
        ' sysdep has already been finalized '
      stop
    end if

!   accumulate the ia vector and compute nnz

    call accumulate_ia_nnz ( sysdep%Eu )
    call accumulate_ia_nnz ( sysdep%Ep )

    sysdep%finalized = .true.

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

  end subroutine finalize_sysdep_structure


! Allocate and initialize sysdep structure (the array ia)

  subroutine initialize_sysdep_structure ( sysdep, numund, numess )

    type(sysdep_t), intent(inout) :: sysdep
    integer, intent(in) :: numund, numess

!   This is a helper routine. Not intended for direct user call.

    sysdep%Eu%n = numess
    sysdep%Eu%m = numund
    sysdep%Eu%nnz = 0
    allocate(sysdep%Eu%ia(numess+1))
    sysdep%Eu%ia = 0

    sysdep%Ep%n = numess
    sysdep%Ep%m = numess
    sysdep%Ep%nnz = 0
    allocate(sysdep%Ep%ia(numess+1))
    sysdep%Ep%ia = 0

    sysdep%initialized_structure = .true.

  end subroutine initialize_sysdep_structure


! Allocate sysdep (the matrices data (a) + column numbers (ja) and vector g)

  subroutine create_sysdep_data ( sysdep )

    type(sysdep_t), intent(inout) :: sysdep

!   test

    if ( .not. sysdep%initialized_structure ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysdep_data: ', &
        ' sysdep structure has not been created '
      stop
    end if

    if ( .not. sysdep%finalized ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysdep_data: ', &
        ' sysdep has not been finalized '
      stop
    end if

    if ( sysdep%allocated_data ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysdep_data: ', &
        ' sysdep data has already been allocated '
      stop
    end if

    allocate( sysdep%Eu%a(sysdep%Eu%nnz) )
    allocate( sysdep%Eu%ja(sysdep%Eu%nnz) )
    allocate( sysdep%Ep%a(sysdep%Ep%nnz) )
    allocate( sysdep%Ep%ja(sysdep%Ep%nnz) )

    allocate( sysdep%g(sysdep%Eu%n) )

    sysdep%allocated_data = .true.

  end subroutine create_sysdep_data


! Create a single new sysdep from an existing one.
! NOTE: the new sysdep contains the structure and all data has been
! allocated. However the data has NOT been copied from the existing one.
! Use copy_sysdep for that.

  subroutine create_single_sysdep ( sysdep1, sysdep2, nr )

    type(sysdep_t), intent(inout) :: sysdep1, sysdep2
    integer, intent(in) :: nr

    if ( .not. sysdep1%initialized_structure ) then
      write(*,'(/2a/)') &
        'Error in create_sysdep: no matrix structure in existing sysdep'
      stop
    end if

    if ( .not. sysdep1%finalized ) then
      write(*,'(/2a/)') &
        'Error in create_sysdep: matrix has not been finalized in ', &
        'existing sysdep'
      stop
    end if

    if ( .not. sysdep1%allocated_data ) then
      write(*,'(/2a/)') &
        'Error in create_sysdep: data not allocated in existing sysdep'
      stop
    end if

    if ( sysdep2%initialized_structure .or. sysdep2%finalized .or. &
         sysdep2%allocated_data ) then
      write(*,'(/2a/a,i0/)') &
        'Error in create_sysdep: new sysdep has already been (partly) ', &
        'created.', &
        ' sysdep number in heading = ', nr
      stop
    end if

!   copy scalars

    sysdep2%initialized_structure = sysdep1%initialized_structure
    sysdep2%finalized = sysdep1%finalized
    sysdep2%allocated_data = sysdep1%allocated_data

    call copy_scalars_sparsematrix ( sysdep1%Eu, sysdep2%Eu )
    call copy_scalars_sparsematrix ( sysdep1%Ep, sysdep2%Ep )

!   allocate arrays

    call allocate_sparsematrix ( sysdep2%Eu )
    call allocate_sparsematrix ( sysdep2%Ep )
    allocate( sysdep2%g(sysdep2%Eu%n) )

!   copy structure

    sysdep2%Eu%ia = sysdep1%Eu%ia
    sysdep2%Ep%ia = sysdep1%Ep%ia

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

  end subroutine create_single_sysdep


! Create new system matrix from existing one
! NOTE: the new sysdep contains the structure and all data has been
! allocated. However the data has NOT been copied from the existing one.
! Use copy_sysdep for that.

  subroutine create_sysdep ( sysdep, sysdep1, sysdep2, sysdep3, &
    sysdep4, sysdep5 )

    type(sysdep_t), intent(inout) :: sysdep, sysdep1
    type(sysdep_t), intent(inout), optional :: sysdep2, sysdep3, &
      sysdep4, sysdep5

    call create_single_sysdep( sysdep, sysdep1, 1 )
    if ( present(sysdep2) ) &
      call create_single_sysdep( sysdep, sysdep2, 2 )
    if ( present(sysdep3) ) &
      call create_single_sysdep( sysdep, sysdep3, 3 )
    if ( present(sysdep4) ) &
      call create_single_sysdep( sysdep, sysdep4, 4 )
    if ( present(sysdep5) ) &
      call create_single_sysdep( sysdep, sysdep5, 5 )

  end subroutine create_sysdep


! Copy system matrix to another system matrix

  subroutine copy_sysdep ( sysdep1, sysdep2 )

    type(sysdep_t), intent(inout) :: sysdep1, sysdep2

!   This routine performs a real data copy of one sysdep to the other
!   by the assignment statement (with allocatables)
!
!      sysdep2 = sysdep1
!

    if ( .not. sysdep1%allocated_data ) then
      write(*,'(/a/)') &
        'Error in copy_sysdep: data not allocated in sysdep1'
      stop
    end if

    sysdep2 = sysdep1

  end subroutine copy_sysdep


! Get additional unknowns in dependency + index from the sysvector.

  subroutine get_sysvector_dependency ( mesh, problem, sysvector, dependency, &
    u, posu, ndofu )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(sysvector_t), intent(in) :: sysvector

!   the dependency number
    integer, intent(in) :: dependency

!   the degrees of freedom, if present
    real(dp), intent(out), dimension(:), optional :: u

!   if present it contains the positions of the components of u in sysvector
    integer, dimension(:), intent(out), optional :: posu

!   if present it gives the number of degrees of freedom stored in u
    integer, intent(out), optional :: ndofu


    integer :: ndof, lpos
    integer, dimension(:), allocatable :: pos

    call check ( problem, 'get_sysvector_dependency', mesh )

!   Some testing

    if ( .not. sysvector%created ) then
      write(*,'(/a/)') &
        'Error in get_sysvector_dependency: sysvector not created.'
      stop
    end if

    if ( problem%probnr /= sysvector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in get_sysvector_dependency: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in sysvector = ',  sysvector%probnr
      stop
    end if

    if ( dependency < 1 .or. dependency > problem%numdependencies ) then
      write(*,'(/2(a/),2(a,i0),/a,i0/)') &
        'Error: dependency in the heading of', &
        ' get_sysvector_dependency is out of range ', &
        ' dependency is ', dependency, &
        ' whereas the number of dependencies is ', problem%numdependencies
      stop
    end if

    lpos = problem%dependencies(dependency)%naddunknowns

    allocate(pos(lpos))

    call pos_array_dependency ( problem, dependency, ndof, pos )

    if ( present(ndofu) ) ndofu = ndof

!   check u array

    if ( present(u) ) then
      if ( size(u) < ndof ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_sysvector_dependency: ', &
          ' u array is too small: size(u) = ', size(u), &
          ' whereas the number of degrees of freedom = ', ndof
        stop
      end if
    end if

!   check posu array

    if ( present(posu) ) then
      if ( size(posu) < ndof ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_sysvector_dependency: ', &
          ' posu array is too small: size(posu) = ', size(posu), &
          ' whereas the number of degrees of freedom = ', ndof
        stop
      end if
    end if

    if ( lpos < ndof ) then
      write(*,'(/a/a/)') &
        'Internal error in get_sysvector_dependency: ', &
        ' lpos < ndof '
      stop
    end if

!   compute u

    if ( present(u) ) u(1:ndof) = sysvector%u(pos(1:ndof))

!   posu array

    if ( present(posu) ) posu(1:ndof) = pos(1:ndof)

    deallocate(pos)

  end subroutine get_sysvector_dependency


! create vector subscript array in sysvector for dependencies (additional
! unknowns only)

  subroutine create_subscript_dependency ( mesh, problem, subscript, &
    dependencies, degfd )

    type(mesh_t), intent(in)  :: mesh
    type(problem_t), intent(in)  :: problem

!   the vector subscript of the degrees of freedom in a sysvector
    type(subscriptdep_t), intent(inout) :: subscript

!   dependencies included for defining the subscript.
!   For example dependencies=(/2,4/), includes dependencies 2 and 4.
!   Default is to include all dependencies.
    integer, dimension(:), intent(in), optional :: dependencies

!   if degfd is present only the degfd^th degree of freedom will be included.
    integer, intent(in), optional :: degfd


    integer :: dof, dep, bp, nndof, ip
    logical :: ldep(problem%numdependencies)
    integer, allocatable, dimension(:) :: work


    call check ( mesh, 'create_subscript_dependency' )
    call check ( problem, 'create_subscript_dependency' )

!   subscript already allocated?
    if (allocated(subscript%s)) deallocate(subscript%s)

!   check dependencies
    if ( present(dependencies) ) then
!     specified dependencies
      if ( any ( dependencies < 1 ) .or. &
           any ( dependencies > problem%numdependencies ) ) then
        write(*,'(/a/a/a,i0/)') &
          'Error: argument dependencies in the heading of', &
          ' create_subscript_dependency is ', &
          ' out of range: some dependencies are < 1 or larger than ', &
           problem%numdependencies
        stop
      end if
      ldep = .false.
      ldep(dependencies) = .true.
    else
!     all dependencies
      ldep = .true.
    end if

    allocate ( work(problem%numdependegfd) )

    dof = 0

    do dep = 1, problem%numdependencies

      if ( .not. ldep(dep) ) cycle  ! do not include dep

!     get positions in array of  additional unknowns

      bp = problem%dependencies(dep)%addnumdegfd(1)
      nndof = problem%dependencies(dep)%addnumdegfd(2) &
                 - problem%dependencies(dep)%addnumdegfd(1)

!     positions in sysvector and sysmatrix (renumbered)

      call fill_work

    end do

    if ( dof == 0 ) then
      write(*,'(/a/)') &
        'Warning in create_subscript_dependency: subscript is empty '
    end if

    subscript%s = work(1:dof)

    deallocate ( work )

  contains

    subroutine fill_work

      integer :: deg

      if ( present(degfd) ) then

!       one degree

        if ( degfd >=1 .and. degfd <= nndof ) then

!         degfd in valid range

          ip =  problem%degfdperm(bp+degfd,2)

!         position in sysvector and sysmatrix (renumbered)
          work(dof+1) = ip

          dof = dof + 1

        end if

      else

!       all degrees

        do deg = 1, nndof

          ip =  problem%degfdperm(bp+deg,2)

!         position in sysvector and sysmatrix (renumbered)
          work(dof+1) = ip

          dof = dof + 1

        end do

      end if

    end subroutine fill_work

  end subroutine create_subscript_dependency


! delete vector subscript array in sysvector for dependencies

  subroutine delete_single_subscript_dependency ( subscript )

!   the vector subscript of the degrees of freedom in a sysvector
    type(subscriptdep_t), intent(out) :: subscript

  end subroutine delete_single_subscript_dependency


! delete vector subscript array in sysvector for dependencies

  subroutine delete_subscript_dependency ( subscript1, subscript2, subscript3, &
    subscript4, subscript5 )

!   the vector subscript of the degrees of freedom in a sysvector
    type(subscriptdep_t), intent(inout) :: subscript1
    type(subscriptdep_t), intent(inout), optional :: subscript2, subscript3, &
      subscript4, subscript5

    call delete_single_subscript_dependency(subscript1)
    if ( present(subscript2) ) &
                  call delete_single_subscript_dependency(subscript2)
    if ( present(subscript3) ) &
                  call delete_single_subscript_dependency(subscript3)
    if ( present(subscript4) ) &
                  call delete_single_subscript_dependency(subscript4)
    if ( present(subscript5) ) &
                  call delete_single_subscript_dependency(subscript5)

  end subroutine delete_subscript_dependency

end module system_dependency_m
