
! Copyright (C) 2004-2022 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Routines, for the system of equations involving constraints

module system_constraint_m

  use kind_defs_m
  use mesh_m
  use problem_defs_m
  use sparse_m
  use system_defs_m
  use system_matrix_m
  use element_defs_m
  use set_optional_m

  implicit none


! If set to .true. a warning is given if there are nodes in constraints on
! objects that do not intersect with the mesh (no reference coordinates).
! NOTE: nodes having no reference coordinates are ignored in the constraint!
  logical, save :: warn_objects_missing_nodes = .true.


! type definition of a vector subscript in the constraint part of a sysvector

  type subscriptcon_t

!   the subscript. Data can be accessed with sysvector%u(subscript%s).
    integer, allocatable, dimension(:) :: s

!   the nodes involved in the subscript that have a non-zero number of
!   degrees of freedom in the subscript. This is only filled when
!   fillnodes=.true. in create_subscript_constraint.
!   NOTE: this applies only to distributed degrees of freedom.
!   NOTE: the values are local node numbers within the geometry or object.
    integer, allocatable, dimension(:) :: nodes

  end type subscriptcon_t


! interface for generic create subroutine

  interface create
    module procedure create_subscript_constraint
  end interface create

  interface create_subscript
    module procedure create_subscript_constraint
  end interface create_subscript


! interface for generic delete subroutine

  interface delete
    module procedure delete_subscript_constraint
  end interface delete


contains


! Create system matrix structure of constraints (ia not yet accumulated)

  subroutine create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem, &
    diagonal, diagonal_block, diagonal_addunknowns, diagonal_block_addunknowns,&
    symmetric )

    type(sysmatrix_t), intent(inout) :: sysmatrix
    type(mesh_t),      intent(in)  :: mesh
    type(problem_t),   intent(in)  :: problem

!   if present and .true. the diagonal of the matrix is included in the
!   system matrix structure for the constraint part, including the additional
!   unknowns part.
!   NOTE: this affects _all_ constraints.
!   default=.false.
    logical, intent(in), optional :: diagonal

!   if present and .true. the diagonal block of the matrix is included in the
!   system matrix structure for the constraint part, including the additional
!   unknowns part.
!   Note, that it is not the full block, but the matrix elements that would
!   naturally be filled with zero when building the contraint point for point
!   (collocation) or element by element (weak).
!   NOTE: this affects _all_ constraints. Setting this for each constraint
!   individually can be done when defining the constraint using
!   define_constraint.
!   default=.false.
    logical, intent(in), optional :: diagonal_block

!   if present and .true. the diagonal of the matrix regarding the
!   the additional unknowns is included in the system matrix structure.
!   NOTE: this affects _all_ constraints.
!   default=.false.
    logical, intent(in), optional :: diagonal_addunknowns

!   if present and .true. the diagonal block of the matrix regarding the
!   the additional unknowns is included in the system matrix structure.
!   NOTE: this affects _all_ constraints. Setting this for each constraint
!   individually can be done when defining the constraint using
!   define_constraint.
!   default=.false.
    logical, intent(in), optional :: diagonal_block_addunknowns

!   if symmetric is present and has the values .true. the matrix is assumed to
!   be symmetric and only the non-zero elements of the upper-triangle
!   (including the diagonal) of the matrix of unknowns (Suu) are stored.
!   default: symmetric=.false.
!   NOTE: this has no effect if sysmatrix has been initialized already.
    logical, intent(in), optional :: symmetric

!   This routine fills the arrays ia with the number of non-zeros in each row:
!
!       ia(i+1) = number of non-zeros in row i
!
!   ia is not yet accumulated and may still be adapted (add more non-zeros)
!   Here the structure of the constraints is created.
!
!   The original equations are
!
!      S u = f
!
!   The constraints are
!
!      A u + B a = g
!
!   where a are the additional unknowns.
!
!   The structure of the matrix system now becomes:
!
!      [  S  A^T  0 ] [ u ] = [ f ]
!      [  A   0   B ] [ l ] = [ g ]
!      [  0  B^T  0 ] [ a ] = [ h ]
!
!   where l are the Lagrange multipliers and h is a right-hand side for the
!   additional unknowns to introduce forcing on the additional unknowns.

    logical :: ldiagonal, ldiagonal_block
    logical :: ldiagonal_addunknowns, ldiagonal_block_addunknowns
    integer :: numund, numess, nodenr, deg_r1, deg_r2, n1, n2, numnod, deg_r
    integer :: row, node, deg_c1, deg_c2, deg_c, rowp
    integer :: constr, nnodes, sizenodes, col, nnodes2
    integer :: object, elgrp, elem, physq, elgrp2, elem2, object2
    integer :: nnodes1, oelem, ipe, numel, nadd
    integer, allocatable, dimension(:) :: utype, work
    integer, allocatable, dimension(:,:) :: wk1, wk2, nwk1, nwk2
    integer, dimension(problem%maxnoddegfd) :: pos
    integer :: nwkelem1, nwkelem2, sizelnodes, numlnod
    integer :: dof, layer, elementset1, elementset2
    integer :: point, nodeset
!   the absolute value of the array nodes contains the nodes that are connected
!   to the Lagrange multipliers in the currently evaluated node. The sign is
!   an indicator whether the unknowns are given by physq1 (>0) or physq2 (<0)
!   in the derive type constraint_t.
    integer, dimension(:), allocatable :: nodes, workc1, workc2
!   the array lnodes contains the local nodes within the geometry or object
!   that are connected to the Lagrange multipliers in the currently
!   evaluated node. Note, that the connection is not just geometrical, but also
!   depends on the discretization (weak or collocation).
!   For global constraints there is only one fake node of value 1.
    integer, dimension(:), allocatable :: lnodes
    type(geometry_t) :: geometry, geometry2

    allocate ( utype(problem%numdegfd), work(problem%numdegfd) )
    allocate ( wk1(mesh%nelgrp,maxval(mesh%grpnumel)) )
    allocate ( wk2(mesh%nelgrp,maxval(mesh%grpnumel)) )
    allocate ( nwk1(mesh%nelem,2), nwk2(mesh%nelem,2) )

    call check ( mesh, 'create_sysmatrix_structure_constraint' )
    call check ( problem, 'create_sysmatrix_structure_constraint', mesh )

    if ( sysmatrix%finalized ) then
      write(*,'(2(/a)/)') &
        'Error in create_sysmatrix_structure_constraint: ', &
        ' sysmatrix has already been finalized '
      stop
    end if

!   test constraints for limit in implementation
!   (limitation: degrees involved in the constraint need to be consecutive, i.e.
!    have "stride=1". This is violated for nphysq>0, physq=0, layer>0)

    if ( problem%nphysq > 0 ) then
!     limit only applies if physical quantities have been defined
      do constr = 1, problem%numconstraints
        if ( problem%constraints(constr)%layer1 > 0 .and. &
             problem%constraints(constr)%physq1 == 0 .or. &
             problem%constraints(constr)%layer2 > 0 .and. &
             problem%constraints(constr)%physq2 == 0 ) then
          write(*,'(3(/a)/a,i0/)') &
            'Error in create_sysmatrix_structure_constraint: ', &
            ' It is not yet possible to specify a layer without specifying ', &
            ' a physical quantity as well.', &
            ' constraint = ', constr
          stop
        end if
      end do
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

!   set work array (degrees present) to zero

    work = 0

!   set to zero work arrays for storing counted elements in weak objects

    nwkelem1 = 0
    wk1 = 0
    nwkelem2 = 0
    wk2 = 0

!   diagonal

    ldiagonal = set_optional ( variable=diagonal, default=.false. )
    ldiagonal_block = &
       set_optional ( variable=diagonal_block, default=.false. )
    ldiagonal_addunknowns = &
       set_optional ( variable=diagonal_addunknowns, default=.false. )
    ldiagonal_block_addunknowns = &
       set_optional ( variable=diagonal_block_addunknowns, default=.false. )

!   initialize?

    if ( .not. sysmatrix%initialized_structure ) then

!     initialize system_matrix

      call initialize_sysmatrix_structure ( sysmatrix, numund, numess )

!     set symmetry of the matrix

      sysmatrix%symmetric = set_optional ( variable=symmetric, default=.false. )

    end if

!   loop over all constraints

    do constr = 1, problem%numconstraints

      if ( problem%constraints(constr)%typeconstraint == 0 ) cycle

      sizenodes = 0
      sizelnodes = 0

!     constraint on what geometrical object?

      if ( problem%constraints(constr)%geometry1 > 0 ) then
!       constraint on geometry
        if ( problem%constraints(constr)%typegeometry == 1 ) then
!         constraint in point; fill fake geometry
          point = problem%constraints(constr)%geometry1
          geometry%nnodes = 1
          geometry%nodes = [mesh%points(point)]
        else if ( problem%constraints(constr)%typegeometry == 2 ) then
!         constraint on curve
          call copy_g ( mesh%curves(problem%constraints(constr)%geometry1), &
                        geometry )
        else if ( problem%constraints(constr)%typegeometry == 3 ) then
!         constraint on surface
          call copy_g ( mesh%surfaces(problem%constraints(constr)%geometry1), &
                        geometry )
        else if ( problem%constraints(constr)%typegeometry == 4 ) then
!         constraint on volume
          call copy_g ( mesh%volumes(problem%constraints(constr)%geometry1), &
                        geometry )
        else if ( problem%constraints(constr)%typegeometry == 5 ) then
!         constraint on nodeset; fill fake geometry
          nodeset = problem%constraints(constr)%geometry1
          geometry%nnodes = size(mesh%nodesets(nodeset)%a)
          geometry%nodes = mesh%nodesets(nodeset)%a
        end if
!       constraint on curve
        nnodes = geometry%nnodes
!       determine size of array nodes and lnodes
        if ( problem%constraints(constr)%typeconstraint == 1 ) then
!         distributed constraint
          sizenodes = geometry%nc%maxnodnumnod+1
          sizelnodes = geometry%nc%maxnodnumnod+1
        else if ( problem%constraints(constr)%typeconstraint == 2 ) then
!         global constraint
          sizenodes = geometry%nnodes
          sizelnodes = 1
        end if
        if ( problem%constraints(constr)%geometry2 > 0 ) then
!         connecting geometry
          if ( problem%constraints(constr)%typegeometry == 1 ) then
  !         constraint in point; fill fake geometry
            point = problem%constraints(constr)%geometry2
            geometry2%nnodes = 1
            geometry2%nodes = [mesh%points(point)]
          else if ( problem%constraints(constr)%typegeometry == 2 ) then
!           constraint on curve
            call copy_g ( mesh%curves(problem%constraints(constr)%geometry2), &
                          geometry2 )
          else if ( problem%constraints(constr)%typegeometry == 3 ) then
!           constraint on surface
            call copy_g ( mesh%surfaces(problem%constraints(constr)%geometry2),&
                          geometry2 )
          else if ( problem%constraints(constr)%typegeometry == 5 ) then
!           constraint on nodeset; fill fake geometry
            nodeset = problem%constraints(constr)%geometry2
            geometry2%nnodes = size(mesh%nodesets(nodeset)%a)
            geometry2%nodes = mesh%nodesets(nodeset)%a
          end if
          if ( problem%constraints(constr)%typeconstraint == 1 ) then
!           distributed constraint
            if ( problem%constraints(constr)%full2 ) then
!             connect to full curve
              sizenodes = sizenodes + geometry2%nnodes
            else
              sizenodes = sizenodes + geometry2%nc%maxnodnumnod + 1
            end if
          else if ( problem%constraints(constr)%typeconstraint == 2 ) then
!           global constraint
            sizenodes = sizenodes + geometry2%nnodes
          end if
        end if
      else if ( problem%constraints(constr)%elementset1 > 0 ) then
!       constraint on elementset
        elementset1 = problem%constraints(constr)%elementset1
        nnodes = mesh%elementsets(elementset1)%nnodes
!       determine size of array nodes and lnodes
        if ( problem%constraints(constr)%typeconstraint == 1 ) then
!         distributed constraint (not yet available)
          stop 'distributed constraint (not yet available) for elementsets'
        else if ( problem%constraints(constr)%typeconstraint == 2 ) then
!         global constraint
          sizenodes = mesh%elementsets(elementset1)%nnodes
          sizelnodes = 1
        end if
        if ( problem%constraints(constr)%elementset2 > 0 ) then
!         connecting elementset
          elementset2 = problem%constraints(constr)%elementset2
          if ( problem%constraints(constr)%typeconstraint == 1 ) then
!           distributed constraint
            stop 'distributed constraint (not yet available) for elementsets'
          else if ( problem%constraints(constr)%typeconstraint == 2 ) then
!           global constraint
            sizenodes = sizenodes + mesh%elementsets(elementset2)%nnodes
          end if
        end if
      else if ( problem%constraints(constr)%object > 0 ) then
!       constraint on object
        object = problem%constraints(constr)%object
        nnodes = mesh%objects(object)%nnodes
!       determine size of array nodes
        if ( problem%constraints(constr)%typeconstraint == 1 ) then
!         distributed constraint
          if ( problem%constraints(constr)%discretization == 0 ) then
!           weak
            sizenodes = mesh%objects(object)%nc%maxnodnumel * &
                        mesh%objects(object)%ninti * &
                        maxval(mesh%elnumnod)
            sizelnodes = mesh%objects(object)%nc%maxnodnumnod + 1
          else if ( problem%constraints(constr)%discretization == 1 ) then
!           collocation
            sizenodes = maxval(mesh%elnumnod)
            sizelnodes = 1
          end if
        else if ( problem%constraints(constr)%typeconstraint == 2 ) then
!         global constraint
          sizenodes = mesh%nelem * maxval(mesh%elnumnod)
          sizelnodes = 1
        end if
        if ( mesh%objects(object)%typeofobject == 2 ) then
!         two-sided connections
          sizenodes = 2*sizenodes
        else if ( problem%constraints(constr)%object2 > 0 ) then
!         connect to other object
          sizenodes = 2*sizenodes
        end if
      end if

      allocate ( nodes(sizenodes), workc1(sizenodes), workc2(sizenodes) )
      allocate ( lnodes(sizelnodes) )

      if ( problem%constraints(constr)%typeconstraint == 1 ) then

!       distributed constraint

!       scan all nodes for degrees of freedom and look for connections

        do nodenr = 1, nnodes

          deg_r1 = problem%constraints(constr)%nodnumdegfd ( nodenr ) + 1
          deg_r2 = problem%constraints(constr)%nodnumdegfd ( nodenr + 1 )

          if ( any(utype(deg_r1:deg_r2) /= 1 ) ) stop 'degr/=1'

!         find all connections to current Lag. multpl.

          if ( problem%constraints(constr)%geometry1 > 0 ) then

!           constraint on geometry

!           nodes are current node + all connections (Lagr mult. connections)
            lnodes(1) = nodenr
            numlnod = 1

            if ( problem%constraints(constr)%discretization == 0 ) then
!             element wise connections (weak form)
              n1 = geometry%nc%nodnumnod(nodenr)
              n2 = geometry%nc%nodnumnod(nodenr+1)
              numlnod = numlnod + n2 - n1
              lnodes(2:numlnod) = geometry%nc%nodnod(n1+1:n2)
            else
!             collocation (no extra nodes)
            end if

!           nodes are current node + all connections
            nodes(1) = geometry%nodes(nodenr)
            numnod = 1

            if ( problem%constraints(constr)%discretization == 0 ) then
!             element wise connections (weak form)
              n1 = geometry%nc%nodnumnod(nodenr)
              n2 = geometry%nc%nodnumnod(nodenr+1)
              numnod = numnod + n2 - n1
              nodes(2:numnod) = geometry%nodes( geometry%nc%nodnod(n1+1:n2) )
            else
  !           collocation (no extra nodes)
            end if

            if ( problem%constraints(constr)%geometry2 > 0 ) then
  !           connection to a second curve
              if ( problem%constraints(constr)%full2 ) then
  !             connect to full curve
                nnodes2 = geometry2%nnodes
                nodes(numnod+1:numnod+nnodes2) = - geometry2%nodes
                numnod = numnod + nnodes2
              else if ( problem%constraints(constr)%discretization == 0 ) then
  !             element wise connections (weak form)
                nodes(numnod+1) = - geometry2%nodes(nodenr)
                numnod = numnod + 1
                n1 = geometry2%nc%nodnumnod(nodenr)
                n2 = geometry2%nc%nodnumnod(nodenr+1)
                nodes(numnod+1:numnod+n2-n1) = - geometry2%nodes( &
                                     geometry2%nc%nodnod(n1+1:n2) )
                numnod = numnod + n2 - n1
              else if ( problem%constraints(constr)%discretization == 1 ) then
  !             collocation
                nodes(numnod+1) = - geometry2%nodes(nodenr)
                numnod = numnod + 1
              end if

            end if

          else if ( problem%constraints(constr)%object > 0 ) then

!           constraints on object

!           nodes are current node + all connections (Lagr mult. connections)
            lnodes(1) = nodenr
            numlnod = 1

            if ( problem%constraints(constr)%discretization == 0 ) then
!             element wise connections (weak form)
              n1 = mesh%objects(object)%nc%nodnumnod(nodenr)
              n2 = mesh%objects(object)%nc%nodnumnod(nodenr+1)
              numlnod = numlnod + n2 - n1
              lnodes(2:numlnod) = mesh%objects(object)%nc%nodnod(n1+1:n2)
            else
!             collocation (no extra nodes)
            end if

!           nodes: all connections in intersecting elements

            numnod = 0

            call nodes_connect_single_node_object

!           set back elements counted in weak constraint on objects

            call reset_wk

          end if

          call count_matrix_entries ( problem%constraints(constr)%nodnumdegfd )

        end do

      else if ( problem%constraints(constr)%typeconstraint == 2 ) then

!       global constraint

        deg_r1 = problem%constraints(constr)%globnumdegfd (1) + 1
        deg_r2 = problem%constraints(constr)%globnumdegfd (2)

        if ( any(utype(deg_r1:deg_r2) /= 1 ) ) stop 'degr/=1'

!       nodes for Lagrange multiplier connections (fake node)
        numlnod = 1
        lnodes(1) = 1

!       nodes are all nodes on geometries, elementsets or objects

        if ( problem%constraints(constr)%geometry1 > 0 ) then

!         constraint on curve or surface

          numnod = nnodes
          nodes(1:numnod) = geometry%nodes

          if ( problem%constraints(constr)%geometry2 > 0 ) then
  !         connection to a second curve
            nnodes2 = geometry2%nnodes
            nodes(nnodes+1:nnodes+nnodes2) = - geometry2%nodes
            numnod = numnod + nnodes2
          end if

        else if ( problem%constraints(constr)%elementset1 > 0 ) then

!         constraint on elementset

          numnod = nnodes
          nodes(1:numnod) = mesh%elementsets(elementset1)%nodes

          if ( problem%constraints(constr)%geometry2 > 0 ) then
  !         connection to a second elementset
            nnodes2 = mesh%elementsets(elementset2)%nnodes
            nodes(nnodes+1:nnodes+nnodes2) = &
                                      - mesh%elementsets(elementset2)%nodes
            numnod = numnod + nnodes2
          end if

        else if ( problem%constraints(constr)%object > 0 ) then

!         constraint on object

!         loop over all nodes of the object

          numnod = 0

          do nodenr = 1, nnodes

            call nodes_connect_single_node_object

          end do

!         set back elements counted in weak constraint on objects

          call reset_wk

        end if

        call count_matrix_entries ( problem%constraints(constr)%globnumdegfd )

      end if

      deallocate ( nodes, workc1, workc2 )
      deallocate ( lnodes )

!     add diagonal in additional unknowns equation

      if ( ldiagonal_block .or. ldiagonal_block_addunknowns .or. &
           problem%constraints(constr)%diagonal_block .or. &
           problem%constraints(constr)%diagonal_block_addunknowns ) then

!       add diagonal_block

        deg_r1 = problem%constraints(constr)%addnumdegfd(1) + 1
        deg_r2 = problem%constraints(constr)%addnumdegfd(2)
        nadd = problem%constraints(constr)%naddunknowns

        do deg_r = deg_r1, deg_r2

          row = problem%degfdperm(deg_r,2)

          if ( sysmatrix%symmetric ) then

!           symmetric matrix

            do deg_c = deg_r1, deg_r2

              col = problem%degfdperm(deg_c,2)

              if ( col < row ) cycle  ! element in lower diagonal

              sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1

            end do

          else

            sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + nadd

          end if

        end do

      else if ( ldiagonal .or. ldiagonal_addunknowns ) then

!       add diagonal

        deg_r1 = problem%constraints(constr)%addnumdegfd(1) + 1
        deg_r2 = problem%constraints(constr)%addnumdegfd(2)

        do deg_r = deg_r1, deg_r2

          row = problem%degfdperm(deg_r,2)

          sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1

        end do

      end if

!     remove memory from geometry (only necessary for typegeometry= 1 or 5)

      if ( allocated(geometry%nodes) ) deallocate(geometry%nodes)
      if ( allocated(geometry2%nodes) ) deallocate(geometry2%nodes)

    end do

    deallocate ( utype, work )
    deallocate ( wk1, wk2 )
    deallocate ( nwk1, nwk2 )

  contains


!   find nodes for connection for a single node (nodenr) of the object

    subroutine nodes_connect_single_node_object

      integer :: el, oip

!     nodes: all connections in intersecting elements

      if ( problem%constraints(constr)%discretization == 0 ) then

!       weak or "element by element" global constraint

!       loop over elements (connected to the node) and integration points

        ipe = mesh%objects(object)%nc%nodnumel(nodenr)
        numel = mesh%objects(object)%nc%nodnumel(nodenr+1) - ipe

        do el = 1, numel

          oelem = mesh%objects(object)%nc%nodelem( ipe + el )

          do oip = 1, mesh%objects(object)%ninti

            elgrp = mesh%objects(object)%grpelm_int(oip,1,oelem)
            elem  = mesh%objects(object)%grpelm_int(oip,2,oelem)

            if ( elgrp == 0 ) then  ! no nodes connected
              if ( warn_objects_missing_nodes ) then
                write(*,'(3(/a)/3(a,i0)/)') &
                  'Warning in create_sysmatrix_structure_constraint: ',&
                  ' reference coordinates in object are missing,', &
                  ' skipping constraint in this node.', &
                  ' object = ', object, ' element = ', oelem, &
                  ' integration point = ', oip
              end if
              cycle
            end if

            if ( wk1(elgrp,elem) == 0 ) then

!             new element found

              wk1(elgrp,elem) = 1
              nwkelem1 = nwkelem1 + 1
              nwk1(nwkelem1,:) = [elgrp,elem]

              nnodes1 = mesh%elnumnod(elgrp)
              nodes(numnod+1:numnod+nnodes1) = &
                                        mesh%topology(elgrp)%a(:,elem)
              numnod = numnod + nnodes1

            end if

            if ( mesh%objects(object)%typeofobject == 2 ) then

!             two-sided connections

              elgrp2 = mesh%objects(object)%grpelm2_int(oip,1,oelem)
              elem2  = mesh%objects(object)%grpelm2_int(oip,2,oelem)

              if ( elgrp2 == 0 ) then   ! no nodes connected
                write(*,'(2(/a)/3(a,i0)/)') &
                  'Error in create_sysmatrix_structure_constraint: ', &
                  ' connecting reference coordinates in two-sided ', &
                  ' object are missing. object = ', object , &
                  ' element = ', oelem, ' integration point = ', oip
                stop
              end if

              if ( wk2(elgrp2,elem2) == 0 ) then

!               new element found

                wk2(elgrp2,elem2) = 1
                nwkelem2 = nwkelem2 + 1
                nwk2(nwkelem2,:) = [elgrp2,elem2]

                nnodes2 = mesh%elnumnod(elgrp2)
                nodes(numnod+1:numnod+nnodes2) = &
                                        - mesh%topology(elgrp2)%a(:,elem2)
                numnod = numnod + nnodes2

              end if

            else if ( problem%constraints(constr)%object2 > 0 ) then

!             connect to other object

              object2 = problem%constraints(constr)%object2
              elgrp2 = mesh%objects(object2)%grpelm_int(oip,1,oelem)
              elem2  = mesh%objects(object2)%grpelm_int(oip,2,oelem)

              if ( elgrp2 == 0 ) then   ! no nodes connected
                write(*,'(2(/a)/3(a,i0)/)') &
                  'Error in create_sysmatrix_structure_constraint: ', &
                  ' connecting reference coordinates in second object',&
                  ' are missing. object = ', object2, &
                  ' element = ', oelem, ' integration point = ', oip
                stop
              end if

              if ( wk2(elgrp2,elem2) == 0 ) then

!               new element found

                wk2(elgrp2,elem2) = 1
                nwkelem2 = nwkelem2 + 1
                nwk2(nwkelem2,:) = [elgrp2,elem2]

                nnodes2 = mesh%elnumnod(elgrp2)
                nodes(numnod+1:numnod+nnodes2) = &
                                        - mesh%topology(elgrp2)%a(:,elem2)
                numnod = numnod + nnodes2

              end if

            end if

          end do
        end do

      else if ( problem%constraints(constr)%discretization == 1 ) then

!       collocation of "node for node" global constraint

        elgrp = mesh%objects(object)%grpelm(nodenr,1)
        elem  = mesh%objects(object)%grpelm(nodenr,2)

        if ( elgrp == 0 ) then  ! no nodes connected
          if ( warn_objects_missing_nodes ) then
            write(*,'(3(/a)/a,i0/)') &
              'Warning in create_sysmatrix_structure_constraint: ', &
              ' reference coordinates in object are missing,', &
              ' skipping constraint in this node.', &
              ' object nodenr = ', nodenr
          end if
          return
        end if

        nnodes1 = mesh%elnumnod(elgrp)
        nodes(numnod+1:numnod+nnodes1) = mesh%topology(elgrp)%a(:,elem)
        numnod = numnod + nnodes1

        if ( mesh%objects(object)%typeofobject == 2 ) then

!         two-sided connections

          elgrp2 = mesh%objects(object)%grpelm2(nodenr,1)
          elem2  = mesh%objects(object)%grpelm2(nodenr,2)

          if ( elgrp2 == 0 ) then   ! no nodes connected
            write(*,'(2(/a)/a,i0/)') &
              'Error in create_sysmatrix_structure_constraint: ', &
              ' connecting reference coordinates in two-sided object', &
              ' are missing. object nodenr = ', nodenr
            stop
          end if

          nnodes2 = mesh%elnumnod(elgrp2)
          nodes(numnod+1:numnod+nnodes2) = &
                                  - mesh%topology(elgrp2)%a(:,elem2)
          numnod = numnod + nnodes2

        else if ( problem%constraints(constr)%object2 > 0 ) then

!         connect to other object

          object2 = problem%constraints(constr)%object2
          elgrp2 = mesh%objects(object2)%grpelm(nodenr,1)
          elem2  = mesh%objects(object2)%grpelm(nodenr,2)

          if ( elgrp2 == 0 ) then   ! no nodes connected
            write(*,'(2(/a)/a,i0/)') &
              'Error in create_sysmatrix_structure_constraint: ', &
              ' connecting reference coordinates in second object', &
              ' are missing. object nodenr = ', nodenr
            stop
          end if

          nnodes2 = mesh%elnumnod(elgrp2)
          nodes(numnod+1:numnod+nnodes2) = &
                                  - mesh%topology(elgrp2)%a(:,elem2)
          numnod = numnod + nnodes2

        end if

      end if

    end subroutine nodes_connect_single_node_object


!   set back elements counted in weak constraint on objects

    subroutine reset_wk

      integer :: i

      do i = 1, nwkelem1
        wk1(nwk1(i,1),nwk1(i,2)) = 0
        nwk1(i,:) = 0
      end do
      nwkelem1 = 0
      do i = 1, nwkelem2
        wk2(nwk2(i,1),nwk2(i,2)) = 0
        nwk2(i,:) = 0
      end do
      nwkelem2 = 0

    end subroutine reset_wk


!   count matrix entries of node connections

    subroutine count_matrix_entries ( nodnumdegfd )

      integer, dimension(:), intent(in) :: nodnumdegfd

      integer :: is1, is2, deg_r, i, deg_c

      do deg_r = deg_r1, deg_r2

        row = problem%degfdperm(deg_r,2)

        if ( ldiagonal_block .or. &
                  problem%constraints(constr)%diagonal_block ) then

!         count matrix entries for the diagonal block of the Lagrange m.

          do i = 1, numlnod

            node = lnodes(i)

            deg_c1 = nodnumdegfd ( node ) + 1
            deg_c2 = nodnumdegfd ( node + 1 )

            do deg_c = deg_c1, deg_c2

              if ( sysmatrix%symmetric ) then
!               symmetric matrix
                col = problem%degfdperm(deg_c,2)
                if ( row > col ) cycle  ! element in lower triangle
              end if

              sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1

            end do

          end do

        else if ( ldiagonal ) then

!         add diagonal in constraint equation

          sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1

        end if

!       count matrix entries of node connections

        do i = 1, numnod

          node = nodes(i)

          if ( node > 0 ) then
            physq = problem%constraints(constr)%physq1
            layer = problem%constraints(constr)%layer1
          else
            physq = problem%constraints(constr)%physq2
            layer = problem%constraints(constr)%layer2
          end if

          node = abs ( node )

          if ( physq > 0 ) then
            call pos_array_local_node ( problem, node, dof, pos, [physq], &
              layer )
          else
            call pos_array_local_node ( problem, node, dof, pos, layer=layer )
          end if

          if ( dof == 0 ) then
!           no degrees, skip node
            is1 = 1
            is2 = 0
          else
            is1 = pos(1)
            is2 = pos(dof)
          end if

!         column degrees: assume degrees are consecutive!
          deg_c1 = problem%nodnumdegfd(node) + is1
          deg_c2 = problem%nodnumdegfd(node) + is2

!         store deg_c1 and deg_c2 for later use

          workc1(i) = deg_c1
          workc2(i) = deg_c2

!         start loop over degrees

          do deg_c = deg_c1, deg_c2

            if ( work(deg_c) /= 0 ) cycle ! degree already counted

            col = problem%degfdperm(deg_c,2)

            if ( utype(deg_c) == 1 ) then

!             unknown degree of freedom (column)

              if ( sysmatrix%symmetric ) then
!               symmetric matrix
                if ( col > row ) then
                  write(*,'(/2(a/))') &
                    'internal error create_sysmatrix_structure_constraint:',&
                    ' col > row, not possible with lagr.multl at vector end'
                  stop
                end if
                sysmatrix%Suu%ia(col+1) = sysmatrix%Suu%ia(col+1) + 1  ! A^T
              else
                sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1  ! A
                sysmatrix%Suu%ia(col+1) = sysmatrix%Suu%ia(col+1) + 1  ! A^T
              end if

            else if ( utype(deg_c) == 3 ) then

!             prescribed degree of freedom (column)

              sysmatrix%Sup%ia(row+1) = sysmatrix%Sup%ia(row+1) + 1  ! A
              rowp = col - numund
              sysmatrix%Spu%ia(rowp+1) = sysmatrix%Spu%ia(rowp+1) + 1  ! A^T

            end if

            work(deg_c) = 1 ! degree counted

          end do

        end do

!       set work (degrees found) back to zero

        do i = 1, numnod
          work( workc1(i):workc2(i) ) = 0
        end do

!       count matrix entries of additional unknowns

        deg_c1 = problem%constraints(constr)%addnumdegfd(1) + 1
        deg_c2 = problem%constraints(constr)%addnumdegfd(2)

        if ( any(utype(deg_c1:deg_c2) /= 1 ) ) stop 'degc/=1'

        do deg_c = deg_c1, deg_c2

          col = problem%degfdperm(deg_c,2)

          if ( sysmatrix%symmetric ) then
!           symmetric matrix
            if ( row > col ) then
              write(*,'(/2(a/))') &
                'internal error create_sysmatrix_structure_constraint:',&
                ' row > col, not possible when add unknowns > lagr mult'
              stop
            end if
            sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1  ! B
          else
            sysmatrix%Suu%ia(row+1) = sysmatrix%Suu%ia(row+1) + 1  ! B
            sysmatrix%Suu%ia(col+1) = sysmatrix%Suu%ia(col+1) + 1  ! B^T
          end if

        end do

      end do

    end subroutine count_matrix_entries

  end subroutine create_sysmatrix_structure_constraint


! Assemble system matrix and system vector for constraints

  subroutine build_system_constraint ( mesh, problem, sysmatrix, sysvector, &
    msysvector, m2sysvector, m3sysvector, elemsub, elemsub1, elemsub2, &
    elemsub3, elemsub4, elemsub5, coefficients, oldvectors, order, &
    constraint1, constraint2, buildmatrix, buildvector, addmatvec, addmat, &
    addvec, diagonal, diagonal_block, diagonal_addunknowns, &
    diagonal_block_addunknowns, transform )

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

!   This is the element subroutine that must be supplied by the calling routine
!   Note that the element matrices and vectors are adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
    optional :: elemsub
    interface
      subroutine elemsub ( mesh, problem, constr, elem, node, matrix, vector, &
        first, last, coefficients, oldvectors, elemmat, elemmat2, elemmatadd, &
        elemvec, elemvecadd )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: constr, elem, node
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
        real(dp), intent(out), dimension(:) :: elemvec, elemvecadd
      end subroutine elemsub
    end interface

!   This is the element subroutine that must be supplied by the calling routine
!   Note that the element matrices and vectors are adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
!   Additional possibilities compared to elemsub are:
!     a. the transposed matrices A^T and B^T must be filled independently on
!        elementlevel. In this way the user is free to devise an unsymmetric
!        matrix for the constraint.
!     b. a contribution to the vector f for the constraint is added.
    optional :: elemsub1
    interface
      subroutine elemsub1 ( mesh, problem, constr, elem, node, matrix, vector, &
        first, last, coefficients, oldvectors, elemmat, elemmatt, elemmat2, &
        elemmat2t, elemmatadd, elemmataddt, elemvecf, elemvecf2, elemvec, &
        elemvecadd )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: constr, elem, node
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat, elemmatt
        real(dp), intent(out), dimension(:,:) :: elemmat2, elemmat2t
        real(dp), intent(out), dimension(:,:) :: elemmatadd, elemmataddt
        real(dp), intent(out), dimension(:) :: elemvecf, elemvecf2
        real(dp), intent(out), dimension(:) :: elemvec, elemvecadd
      end subroutine elemsub1
    end interface

!   This is the element subroutine that must be supplied by the calling routine
!   Note that the element matrices and vectors are adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
!   Additional possibilities compared to elemsub are:
!     a. the matrices C and D on the diagonal:
!              [  S  A^T  0 ] [ u ] = [ f ]
!              [  A   C   B ] [ l ] = [ g ]
!              [  0  B^T  D ] [ a ] = [ h ]
!        can be filled on elementlevel. In elemsub it is assumed that C=0 and
!        D=0.
!        NOTE, that either diagonal_block=.true. or
!              diagonal_block_addunknowns=.true., otherwise the
!              elementmatrices for C and D are ignored and not taken into
!              account.
    optional :: elemsub2
    interface
      subroutine elemsub2 ( mesh, problem, constr, elem, node, matrix, vector, &
        first, last, coefficients, oldvectors, elemmat, elemmat2, elemmatdiag, &
        elemmatadd, elemmatdiagadd, elemvec, elemvecadd )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: constr, elem, node
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
        real(dp), intent(out), dimension(:,:) :: elemmatdiag, elemmatdiagadd
        real(dp), intent(out), dimension(:) :: elemvec, elemvecadd
      end subroutine elemsub2
    end interface

!   This is the element subroutine that must be supplied by the calling routine
!   Note that the element matrices and vectors are adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
!   Additional possibilities compared to elemsub are:
!     a. the transposed matrices A^T and B^T must be filled independently on
!        elementlevel. In this way the user is free to devise an unsymmetric
!        matrix for the constraint.
!     b. a contribution to the vector f for the constraint is added.
!     c. the matrices C and D on the diagonal:
!              [  S  A^T  0 ] [ u ] = [ f ]
!              [  A   C   B ] [ l ] = [ g ]
!              [  0  B^T  D ] [ a ] = [ h ]
!        can be filled on elementlevel. In elemsub it is assumed that C=0 and
!        D=0.
!        NOTE, that either diagonal_block=.true. or
!              diagonal_block_addunknowns=.true., otherwise the
!              elementmatrices for C and D are ignored and not taken into
!              account.
    optional :: elemsub3
    interface
      subroutine elemsub3 ( mesh, problem, constr, elem, node, matrix, vector, &
        first, last, coefficients, oldvectors, elemmat, elemmatt, elemmat2, &
        elemmat2t, elemmatdiag, elemmatadd, elemmataddt, elemmatdiagadd, &
        elemvecf, elemvecf2, elemvec, elemvecadd )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: constr, elem, node
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat, elemmatt
        real(dp), intent(out), dimension(:,:) :: elemmat2, elemmat2t
        real(dp), intent(out), dimension(:,:) :: elemmatadd, elemmataddt
        real(dp), intent(out), dimension(:,:) :: elemmatdiag, elemmatdiagadd
        real(dp), intent(out), dimension(:) :: elemvecf, elemvecf2
        real(dp), intent(out), dimension(:) :: elemvec, elemvecadd
      end subroutine elemsub3
    end interface

!   This is the element subroutine that must be supplied by the calling routine
!   Note that the element matrices and vectors are adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
!   Additional possibilities compared to elemsub are:
!     a. the matrices S, C and D on the diagonal:
!              [  S  A^T  0 ] [ u ] = [ f ]
!              [  A   C   B ] [ l ] = [ g ]
!              [  0  B^T  D ] [ a ] = [ h ]
!        can be filled on elementlevel. In elemsub it is assumed that C=0 and
!        D=0, and S is filled only by the standard elements.
!        NOTE, filling of the diagonal block S has been limited to the diagonal
!              block concerning degrees on the same "side" of the constraint,
!              i.e. the degrees concerning the matrix A1 and A2 have a separate
!              diagonal block S1 and S2 that can be filled. Also, no additional
!              storage with respect to the standard fem storage is reserved.
!              Therefore, constraint%full2=.true. is not allowed.
!        NOTE, that either diagonal_block=.true. or
!              diagonal_block_addunknowns=.true., otherwise the
!              elementmatrices for C and D are ignored and not taken into
!              account.
    optional :: elemsub4
    interface
      subroutine elemsub4 ( mesh, problem, constr, elem, node, matrix, vector, &
        first, last, coefficients, oldvectors, elemmatsdiag, elemmatsdiag2, &
        elemmat, elemmat2, elemmatdiag, elemmatadd, elemmatdiagadd, elemvec, &
        elemvecadd )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: constr, elem, node
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmatsdiag, elemmatsdiag2
        real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
        real(dp), intent(out), dimension(:,:) :: elemmatdiag, elemmatdiagadd
        real(dp), intent(out), dimension(:) :: elemvec, elemvecadd
      end subroutine elemsub4
    end interface

!   This is the element subroutine that must be supplied by the calling routine
!   Note that the element matrices and vectors are adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
!   Additional possibilities compared to elemsub are:
!     a. the transposed matrices A^T and B^T must be filled independently on
!        elementlevel. In this way the user is free to devise an unsymmetric
!        matrix for the constraint.
!     b. a contribution to the vector f for the constraint is added.
!     c. the matrices S, C and D on the diagonal:
!              [  S  A^T  0 ] [ u ] = [ f ]
!              [  A   C   B ] [ l ] = [ g ]
!              [  0  B^T  D ] [ a ] = [ h ]
!        can be filled on elementlevel. In elemsub it is assumed that C=0 and
!        D=0, and S is filled only by the standard elements.
!        NOTE, filling of the diagonal block has S been limited to the diagonal
!              block concerning degrees on the same "side" of the constraint,
!              i.e. the degrees concerning the matrix A1 and A2 have a separate
!              diagonal block S1 and S2 that can be filled. Also, no additional
!              storage with respect to the standard fem storage is reserved.
!              Therefore, constraint%full2=.true. is not allowed.
!        NOTE, that either diagonal_block=.true. or
!              diagonal_block_addunknowns=.true., otherwise the
!              elementmatrices for C and D are ignored and not taken into
!              account.
    optional :: elemsub5
    interface
      subroutine elemsub5 ( mesh, problem, constr, elem, node, matrix, vector, &
        first, last, coefficients, oldvectors, elemmatsdiag, elemmatsdiag2, &
        elemmat, elemmatt, elemmat2, elemmat2t, elemmatdiag, elemmatadd, &
        elemmataddt, elemmatdiagadd, elemvecf, elemvecf2, elemvec, elemvecadd )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: constr, elem, node
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmatsdiag, elemmatsdiag2
        real(dp), intent(out), dimension(:,:) :: elemmat, elemmatt
        real(dp), intent(out), dimension(:,:) :: elemmat2, elemmat2t
        real(dp), intent(out), dimension(:,:) :: elemmatadd, elemmataddt
        real(dp), intent(out), dimension(:,:) :: elemmatdiag, elemmatdiagadd
        real(dp), intent(out), dimension(:) :: elemvecf, elemvecf2
        real(dp), intent(out), dimension(:) :: elemvec, elemvecadd
      end subroutine elemsub5
    end interface

!   if coefficients is present, it is passed through to the element subroutine
    type(coefficients_t), intent(in), optional :: coefficients

!   if oldvectors is present, it is passed through to the element subroutine.
!   This structure can be used to supply old sysvectors/vectors from possibly
!   other problems to the element subroutines. Note that the problem structures
!   of the other problems must be included in oldvectors.
    type(oldvectors_t), intent(in), optional :: oldvectors

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
!
!   Note that in case of a weak constraint the sequence of the Lagrange
!   multipliers will also be affected for order = 'ND' or 'DN'. Although for
!   this case the outside loop of physical quantities is not applicable
!   the loops over degrees and nodes are affected.
    character(len=*), intent(in), optional :: order

!   if these are present the assembling takes place for constraints
!   constraint1,...,constraint2 only. If only constraint1 is present
!   one constraint is assembled only.
    integer, intent(in), optional :: constraint1, constraint2

!   if one of these is assigned the value .false. the assembling of
!   the matrix or vector is not done.
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

!   if present and .true. the diagonal of the matrix is included in the
!   system matrix for the constraint and additional unknown part.
!   NOTE: the subroutine create_sysmatrix_structure_constraint must have been
!   called with diagonal=.true. also.
!   NOTE: the values are set to zero.
!   default=.false.
    logical, intent(in), optional :: diagonal

!   if present and .true. the diagonal block of the matrix is included in the
!   system matrix for the constraint and additional unknown part.
!   NOTE: the subroutine create_sysmatrix_structure_constraint must have been
!   called with diagonal_block=.true. also.
!   NOTE: the values are set to zero when using elemsub or elemsub1.
!   When using elemsub2--elemsub5 the values can be set on elementlevel.
!   NOTE: this affects _all_ constraints. Setting this for each constraint
!   individually can be done when defining the constraint using
!   define_constraint.
!   default=.false.
    logical, intent(in), optional :: diagonal_block

!   if present and .true. the diagonal of the matrix regarding the
!   the additional unknowns is included in the system matrix structure.
!   NOTE: the subroutine create_sysmatrix_structure_constraint must have been
!   called with diagonal_addunknowns=.true. also.
!   NOTE: the values are set to zero.
!   default=.false.
    logical, intent(in), optional :: diagonal_addunknowns

!   if present and .true. the diagonal block of the matrix regarding the
!   the additional unknowns is included in the system matrix structure.
!   NOTE: the subroutine create_sysmatrix_structure_constraint must have been
!   called with diagonal_block_addunknowns=.true. also.
!   NOTE: the values are set to zero when using elemsub or elemsub1.
!   When using elemsub2--elemsub5 the values can be set on elementlevel.
!   NOTE: this affects _all_ constraints. Setting this for each constraint
!   individually can be done when defining the constraint using
!   define_constraint.
!   default=.false.
    logical, intent(in), optional :: diagonal_block_addunknowns

!   By setting transform=.false. the transformation matrix is not applied
!   to the element matrix and vector. This means, the element degrees of
!   freedom are defined in the local (transformed) system.
!   default = .true.
    logical, intent(in), optional :: transform


!   The original equations are
!
!      S u = f
!
!   The constraints are
!
!      A u + B a = g
!
!   where a are the additional unknowns.
!
!   The structure of the matrix system including the constraints now becomes:
!
!      [  S  A^T  0 ] [ u ] = [ f ]
!      [  A   0   B ] [ l ] = [ g ]
!      [  0  B^T  0 ] [ a ] = [ h ]
!
!   where l are the Lagrange multipliers and h is a right-hand side for the
!   additional unknowns to introduce forcing on the additional unknowns.
!   Using elemsub/elemsub2/elemsub4 the matrices A and B and the vectors g and h
!   can be supplied on elementlevel.
!   Using elemsub1/elemsub3/elemsub5 the matrices A, A^T, B and B^T and the
!   vectors f (constraint contribution only), g and h can be supplied on
!   elementlevel.
!
!   For some applications the diagonal "zero blocks" need to be filled with
!   non-zero values. The structure of the matrix system now becomes:
!
!      [  S  A^T  0 ] [ u ] = [ f ]
!      [  A   C   B ] [ l ] = [ g ]
!      [  0  B^T  D ] [ a ] = [ h ]
!
!   Using elemsub2 or elemsub3 the matrices C and D can be supplied on
!   elementlevel.
!   Using elemsub4 or elemsub5 the matrices S, C and D can be supplied on
!   elementlevel.

    logical :: matrix, vector, laddmat, laddvec, first, last
    logical :: ldiagonal, ldiagonal_block
    logical :: ldiagonal_addunknowns, ldiagonal_block_addunknowns
    logical :: transposedmatf, smat, diagmat, ltransform
    logical, dimension(:), allocatable :: lp
    integer :: constr1, constr2, elem, constr, ndofr, ndofc
    integer :: ndofc2, ndofcadd, nodenr, nodenr2
    integer :: numund, numess, row, rowg, physq1, physq2, object, elgrp, elgrp2
    integer :: dofr, dofc, dofc2, dofcadd, nloop, iloop, node, elemc, elemc2
    integer :: nrhsd, rhsd, object2, i, j, si, sj
    integer :: layer1, layer2, elementset1, elementset2, el, grp, nel
    integer :: point, nodeset, k
!   work5: array for denoting transformed degrees of freedom
    logical, allocatable, dimension(:) :: work5
    integer, allocatable, dimension(:) :: posr, posc, posc2, poscadd, w1, w2
    integer, allocatable, dimension(:) :: work1, work2, work3, work4
    real(dp), allocatable, dimension(:,:) :: elemmatsdiag, elemmatsdiag2
    real(dp), allocatable, dimension(:,:) :: elemmat, elemmat2, elemmatadd, &
      elemmatdiag, elemmatdiagadd
    real(dp), allocatable, dimension(:,:) :: elemmatt, elemmat2t, elemmataddt
    real(dp), allocatable, dimension(:) :: elemvecf, elemvecf2
    real(dp), allocatable, dimension(:) :: elemvec, elemvecadd
    real(dp), allocatable, dimension(:) :: diagl, diaga
    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl
    type(geometry_t) :: geometry, geometry2
    type(sparsematrix_t) :: A


    if ( problem%numconstraints == 0 ) return

    call check ( mesh, 'build_system_constraint' )
    call check ( problem, 'build_system_constraint', mesh )

!   initialize local parameters

    nrhsd = 1

!   diagonal

    ldiagonal = set_optional ( variable=diagonal, default=.false. )
    ldiagonal_block = &
       set_optional ( variable=diagonal_block, default=.false. )
    ldiagonal_addunknowns = &
       set_optional ( variable=diagonal_addunknowns, default=.false. )
    ldiagonal_block_addunknowns = &
       set_optional ( variable=diagonal_block_addunknowns, default=.false. )

!   additional possibilities for elemsubx, with x=1,2,3,4,5

!   build "transposed" matrices A^T, B^T and f by user on element level
    transposedmatf = present(elemsub1) .or. present(elemsub3) .or. &
                     present(elemsub5)

!   build diagonal block matrices C and D
    diagmat = present(elemsub2) .or. present(elemsub3) .or. &
              present(elemsub4) .or. present(elemsub5)

!   build diagonal block S
    smat = present(elemsub4) .or. present(elemsub5)

!   ordering

    if ( present(order) ) then
      if ( all ( order /= [ 'ND', 'DN' ] ) ) then
        write(*,'(2(/a)/)') &
          'Error in build_system_constraint: ', &
          ' heading parameter order must be either ''ND'' or ''DN''.'
        stop
      end if
    end if

!   initialize

    if ( present(oldvectors) ) oldvl = oldvectors
    if ( present(coefficients) ) coeffl = coefficients

!   which constraints?

    if ( present(constraint1) .and. present(constraint2) ) then
!     specified range of constraints only
      constr1 = constraint1
      constr2 = constraint2
    else if ( present(constraint1) ) then
!     one constraints only
      constr1 = constraint1
      constr2 = constraint1
    else
!     all constraints
      constr1 = 1
      constr2 = problem%numconstraints
    end if

    if ( constr1 < 1 .or. constr1 > problem%numconstraints .or.  &
         constr2 < 1 .or. constr2 > problem%numconstraints ) then

      write(*,'(/2(a/),2(a,i0),/a,i0/)') &
        'Error: constraint1 and/or constraint2 in the heading of', &
        ' build_system_constraint is out of range ', &
        ' constraint1 is ', constr1, '; constraint2 is ', constr2, &
        ' whereas the number of constraints is ', problem%numconstraints
      stop

    end if

!   test elemsub

    if ( count( [ present(elemsub), present(elemsub1), &
                  present(elemsub2), present(elemsub3), &
                  present(elemsub4), present(elemsub5) ] ) /= 1 ) then
      write(*,'(3(/a)/)') &
        'Error in build_system_constraint: ', &
        ' one of elemsub, elemsub1, elemsub2, elemsub3, elemsub4, elemsub5 ', &
        ' must be present in the heading '
      stop
    end if

    if ( smat .and. any(problem%constraints(constr1:constr2)%full2) ) then
      write(*,'(3(/a)/)') &
        'Error in build_system_constraint: ', &
        ' building the matrix S cannot be combined with ', &
        ' problem%constraints(constr1:constr2)%full2 = .true. '
      stop
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
          'Error in build_system_constraint: no system matrix structure.'
        stop
      end if

      if ( .not. sysmatrix%finalized ) then
        write(*,'(/2a/)') &
          'Error in build_system_constraint:', &
          ' system matrix has not been finalized.'
        stop
      end if

      if ( .not. sysmatrix%allocated_data ) then
        write(*,'(/2a/)') &
          'Error in build_system_constraint:', &
          ' data in system matrix not allocated.'
        stop
      end if

      if ( sysmatrix%Suu%m + sysmatrix%Sup%m /= problem%numdegfd ) then
        write(*,'(/a/a/)') &
          'Error in build_system_constraint: number of degrees of freedom', &
          ' of the system matrix is different from the number in problem.'
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
            'Error in build_system_constraint: sysvector not created.'
          stop
        end if

        if ( sysvector%n /= problem%numdegfd ) then
          write(*,'(/a/a/)') &
            'Error in build_system_constraint: number of degrees of freedom', &
            ' of the sysvector vector different from the number in problem.'
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
            'Error in build_system_constraint: msysvector not created.'
          stop
        end if

        if ( any( msysvector%n /= problem%numdegfd ) ) then
          write(*,'(/2a/a/)') &
            'Error in build_system_constraint: number of degrees of ', &
            ' freedom of the ', &
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
            'Error in build_system_constraint: m2sysvector not created.'
          stop
        end if

        if ( any( m2sysvector%n /= problem%numdegfd ) ) then
          write(*,'(/a/a/)') &
            'Error in build_system_constraint: number of degrees of freedom', &
            ' in m2sysvector vector different from the number in problem.'
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
            'Error in build_system_constraint: m3sysvector not created.'
          stop
        end if

        if ( any( m3sysvector%n /= problem%numdegfd ) ) then
          write(*,'(/a/a/)') &
            'Error in build_system_constraint: number of degrees of freedom', &
            ' in m3sysvector vector different from the number in problem.'
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
          'Error in build_system_constraint: ', &
          ' no sysvector, msysvector, m2sysvector or m3sysvector present'
        stop

      end if

    end if

!   transformations

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in build_system_constraint: ', &
        ' transformation matrix has not been build'
      stop
    end if

    ltransform = set_optional ( variable=transform, default=.true. )

    if ( problem%numtransdegfd > 0 .and. ltransform ) then

!     work array to indicate transformed degrees of freedom
      allocate ( work5(problem%numdegfd) )
      work5 = .false.
      work5(problem%degfdtrans) = .true.

    end if

!   temporary work arrays to store information on current large matrix row

    numess = problem%numessdegfd
    numund = problem%numundegfd

    allocate( w1(numund), w2(numess) )

    w1 = 0
    w2 = 0

    if ( problem%numlayers > 0 ) then
      allocate ( lp( max(mesh%maxelnumnod, &
                     maxval(mesh%curves(1:mesh%ncurves)%nnodes), &
                     maxval(mesh%surfaces(1:mesh%nsurfaces)%nnodes), &
                     maxval(mesh%volumes(1:mesh%nvolumes)%nnodes) ) ) )
    else
      allocate ( lp(0) )
    end if

!   start loop over all constraints

    do constr = constr1, constr2

      if ( problem%constraints(constr)%typeconstraint == 0 ) cycle

      physq1 = problem%constraints(constr)%physq1
      physq2 = problem%constraints(constr)%physq2
      layer1 = problem%constraints(constr)%layer1
      layer2 = problem%constraints(constr)%layer2

!     find (maximum) number of degrees of freedom (Lagrange multp)
!     in an element of this constraint

      if ( problem%constraints(constr)%typeconstraint == 1 ) then
!       distributed constraint
        if ( problem%constraints(constr)%discretization == 0 ) then
!         weak constraint
          ndofr = sum ( problem%constraints(constr)%elnumdegfd )
        else if ( problem%constraints(constr)%discretization == 1 ) then
!         collocation
          ndofr = problem%constraints(constr)%maxnumdegfd
        end if
      else if ( problem%constraints(constr)%typeconstraint == 2 ) then
!       global constraint
        ndofr = problem%constraints(constr)%nglobalc
      end if

      if ( problem%constraints(constr)%geometry1 > 0 ) then

!       constraint on geometry

        if ( problem%constraints(constr)%typegeometry == 1 ) then
!         constraint in point; fill fake geometry
          point = problem%constraints(constr)%geometry1
          geometry%nnodes = 1
          geometry%nodes = [mesh%points(point)]
        else if ( problem%constraints(constr)%typegeometry == 2 ) then
!         constraint on curve
          call copy_g ( mesh%curves(problem%constraints(constr)%geometry1), &
                        geometry )
        else if ( problem%constraints(constr)%typegeometry == 3 ) then
!         constraint on surface
          call copy_g ( mesh%surfaces(problem%constraints(constr)%geometry1), &
                        geometry )
        else if ( problem%constraints(constr)%typegeometry == 4 ) then
!         constraint on volume
          call copy_g ( mesh%volumes(problem%constraints(constr)%geometry1), &
                        geometry )
        else if ( problem%constraints(constr)%typegeometry == 5 ) then
!         constraint on nodeset; fill fake geometry
          nodeset = problem%constraints(constr)%geometry1
          geometry%nnodes = size(mesh%nodesets(nodeset)%a)
          geometry%nodes = mesh%nodesets(nodeset)%a
        end if

        if ( problem%constraints(constr)%discretization == 0 ) then
!         element wise connections (weak form)
          call ndof_geometry_elem ( ndofc, geometry, physq1 )
        else
!         collocation
          call ndof_geometry_node ( ndofc, geometry, physq1 )
        end if

        if ( problem%constraints(constr)%geometry2 > 0 ) then
!         connection to a second geometry
          if ( problem%constraints(constr)%typegeometry == 1 ) then
!           constraint in point; fill fake geometry
            point = problem%constraints(constr)%geometry2
            geometry2%nnodes = 1
            geometry2%nodes = [mesh%points(point)]
          else if ( problem%constraints(constr)%typegeometry == 2 ) then
!           constraint on curve
            call copy_g ( mesh%curves(problem%constraints(constr)%geometry2), &
                          geometry2 )
          else if ( problem%constraints(constr)%typegeometry == 3 ) then
!           constraint on surface
            call copy_g ( mesh%surfaces(problem%constraints(constr)%geometry2),&
                          geometry2 )
          else if ( problem%constraints(constr)%typegeometry == 5 ) then
!           constraint on nodeset; fill fake geometry
            nodeset = problem%constraints(constr)%geometry2
            geometry2%nnodes = size(mesh%nodesets(nodeset)%a)
            geometry2%nodes = mesh%nodesets(nodeset)%a
          end if
          if ( problem%constraints(constr)%full2 ) then
!           connect to full curve
            call ndof_fullgeometry ( ndofc2, geometry2, physq2 )
          else if ( problem%constraints(constr)%discretization == 0 ) then
!           element wise connections (weak form)
            call ndof_geometry_elem ( ndofc2, geometry2, physq2 )
          else if ( problem%constraints(constr)%discretization == 1 ) then
!           collocation
            call ndof_geometry_node ( ndofc2, geometry2, physq2 )
          end if
        else
          ndofc2 = 0
        end if

        ndofcadd = problem%constraints(constr)%naddunknowns

      else if ( problem%constraints(constr)%elementset1 > 0 ) then

!       constraint on elementset

        elementset1 = problem%constraints(constr)%elementset1

        if ( problem%constraints(constr)%discretization == 0 ) then

!         element wise (weak form)

          call ndof_elem ( ndofc, physq1 )

        else

!         collocation

          call ndof_elementset_node ( ndofc, mesh%elementsets(elementset1), &
            physq1 )

        end if

        if ( problem%constraints(constr)%elementset2 > 0 ) then

!         connect to second elementset

          elementset2 = problem%constraints(constr)%elementset2

          if ( problem%constraints(constr)%discretization == 0 ) then

!           element wise (weak form)

            call ndof_elem ( ndofc2, physq2 )

          else

!           collocation

            call ndof_elementset_node ( ndofc2, mesh%elementsets(elementset2), &
              physq2 )

          end if

        else

          ndofc2 = 0

        end if

        ndofcadd = problem%constraints(constr)%naddunknowns

      else if ( problem%constraints(constr)%object > 0 ) then

!       constraint on object

        object = problem%constraints(constr)%object

!       weak and collocation is basically the same (single point matrix)

        call ndof_elem ( ndofc, physq1 )

        if ( mesh%objects(object)%typeofobject == 2 .or. &
             problem%constraints(constr)%object2 > 0 ) then

!         two-sided connections or connect to other object

          call ndof_elem ( ndofc2, physq2 )

        else

          ndofc2 = 0

        end if

        ndofcadd = problem%constraints(constr)%naddunknowns

      end if

!     reserve memory for element matrix and vector and positions


      allocate ( elemmat(ndofr,ndofc), elemmat2(ndofr,ndofc2) )
      allocate ( elemmatadd(ndofr,ndofcadd) )
      allocate ( elemvec(ndofr*nrhsd), elemvecadd(ndofcadd*nrhsd) )
      allocate ( posr(ndofr), posc(ndofc), posc2(ndofc2), poscadd(ndofcadd) )

      if ( smat ) then
        allocate ( elemmatsdiag(ndofc,ndofc), elemmatsdiag2(ndofc2,ndofc2) )
      end if

      if ( ldiagonal .or. ldiagonal_addunknowns ) then
        allocate ( diagl(ndofr), diaga(ndofcadd) )
      end if

      if ( ldiagonal_block .or. ldiagonal_block_addunknowns .or. diagmat .or. &
           problem%constraints(constr)%diagonal_block .or. &
           problem%constraints(constr)%diagonal_block_addunknowns ) then
        allocate ( elemmatdiag(ndofr,ndofr), elemmatdiagadd(ndofcadd,ndofcadd) )
      end if

      if ( transposedmatf ) then
        allocate ( elemmatt(ndofc,ndofr), elemmat2t(ndofc2,ndofr) )
        allocate ( elemmataddt(ndofcadd,ndofr) )
        allocate ( elemvecf(ndofc*nrhsd) )
        allocate ( elemvecf2(ndofc2*nrhsd) )
      end if

!     loop over elements or nodes in this constraint

      if ( problem%constraints(constr)%discretization == 0 ) then
!       weak constraint
        if ( problem%constraints(constr)%geometry1 > 0 ) then
!         constraint on geometry
          nloop = geometry%nelem
        else if ( problem%constraints(constr)%elementset1 > 0 ) then
!         constraint on elementset
          nloop = mesh%elementsets(elementset1)%nelem
        else if ( problem%constraints(constr)%object > 0 ) then
!         constraint on object
          nloop = mesh%objects(object)%nelem * mesh%objects(object)%ninti
        end if
      else if ( problem%constraints(constr)%discretization == 1 ) then
!       collocation
        if ( problem%constraints(constr)%geometry1 > 0 ) then
!         constraint on curve
          nloop = geometry%nnodes
        else if ( problem%constraints(constr)%elementset1 > 0 ) then
!         constraint on elementset
          nloop = mesh%elementsets(elementset1)%nnodes
        else if ( problem%constraints(constr)%object > 0 ) then
!         constraint on object
          nloop = mesh%objects(object)%nnodes
        end if
      end if

!     initialize

      if ( problem%constraints(constr)%discretization == 0 ) then
!       weak constraint
        if ( problem%constraints(constr)%elementset1 > 0 ) then
!         weak constraint on elementset
          allocate ( work1(nloop), work2(nloop) )
          el = 0
          do grp = 1, mesh%nelgrp
            nel = mesh%elementsets(elementset1)%grpnumel(grp)
            work1(el+1:el+nel) = [(i,i=1,nel)]  ! element numbers
            work2(el+1:el+nel) = grp      ! group numbers
            el = el + nel
          end do
          if ( problem%constraints(constr)%elementset2 > 0 ) then
!           second connected elementset
            allocate ( work3(nloop), work4(nloop) )
            el = 0
            do grp = 1, mesh%nelgrp
              nel = mesh%elementsets(elementset2)%grpnumel(grp)
              work3(el+1:el+nel) = [(i,i=1,nel)]  ! element numbers
              work4(el+1:el+nel) = grp      ! group numbers
              el = el + nel
            end do
          end if
        else if ( problem%constraints(constr)%object > 0 ) then
!         weak constraint on object
          elem = 1; node = 0
        end if
      end if

!     start loop

      do iloop = 1, nloop

        first = iloop == 1
        last  = iloop == nloop

!       set elem and node Lagrange multipliers

        if ( problem%constraints(constr)%discretization == 0 ) then
!         weak constraint
          if ( problem%constraints(constr)%geometry1 > 0 ) then
!           weak constraint on geometry
            elem = iloop
            node = 0
          else if ( problem%constraints(constr)%elementset1 > 0 ) then
!           weak constraint on elementset
            elem = work1(iloop)
            node = work2(iloop)  ! node is now the group number !!
            first = elem == 1
            last  = elem == mesh%elementsets(elementset1)%grpnumel(node)
          else if ( problem%constraints(constr)%object > 0 ) then
!           weak constraint on object
            node = node + 1
            if ( node > mesh%objects(object)%ninti ) then
!             enter new element
              node = 1
              elem = elem + 1
            end if
          end if
        else if ( problem%constraints(constr)%discretization == 1 ) then
!         collocation
          elem = 0
          node = iloop
        end if

!       compute positions in large matrix/vector of element degrees of freedom
!       for the Lagrange multipliers

        if ( problem%constraints(constr)%typeconstraint == 1 ) then
!         distributed constraint
          if ( problem%constraints(constr)%discretization == 0 ) then
!           weak constraint
            if ( problem%constraints(constr)%geometry1 > 0 ) then
!             weak constraint on geometry
              call pos_array_constraint ( problem, constr, dofr, posr, &
                geometry=geometry, elem=elem, order=order )
            else if ( problem%constraints(constr)%object > 0 ) then
!             weak constraint on object
              call pos_array_constraint ( problem, constr, dofr, posr, &
                object=mesh%objects(object), elem=elem, order=order )
            end if
          else if ( problem%constraints(constr)%discretization == 1 ) then
!           collocation
            call pos_array_constraint ( problem, constr, dofr, posr, &
              node=node )
          end if
        else if ( problem%constraints(constr)%typeconstraint == 2 ) then
!         global constraint
          call pos_array_constraint ( problem, constr, dofr, posr, &
            globalc=.true. )
        end if

!       if dofr == 0 there are no constraints to fill here
        if ( dofr == 0 ) cycle

!       compute positions in large matrix/vector of element degrees of freedom
!       for the columns (the unknowns in the problem)

        if ( problem%constraints(constr)%discretization == 0 ) then

!         weak constraint

          if ( problem%constraints(constr)%geometry1 > 0 ) then

!           constraint on geometry

            call pos_array_geometryc ( geometry, dofc, posc, physq1, layer1 )

            if ( problem%constraints(constr)%geometry2 > 0 ) then

!             connection to a second curve

              if ( problem%constraints(constr)%full2 ) then

!               connect to full curve

                call pos_array_fullgeometryc ( geometry2, dofc2, posc2, &
                  physq2, layer2 )

              else

!               connect to element of second curve

                call pos_array_geometryc ( geometry2, dofc2, posc2, physq2, &
                  layer2 )

              end if

            else

              dofc2 = 0

            end if

          else if ( problem%constraints(constr)%elementset1 > 0 ) then

!           constraint on elementset

            elgrp = node
            elemc = mesh%elementsets(elementset1)%elements(elgrp)%a(elem)

            call pos_array_elemc ( elgrp, elemc, dofc, posc, physq1, layer1 )

            if ( problem%constraints(constr)%elementset2 > 0 ) then

!             connected elementset

              first = .true.; last = .true.
              elgrp2 = work3(iloop)
              elemc2 = &
                mesh%elementsets(elementset2)%elements(elgrp2)%a(work4(iloop))

              call pos_array_elemc ( elgrp2, elemc2, dofc2, posc2, physq2, &
                layer2 )

            else

              dofc2 = 0

            end if

          else if ( problem%constraints(constr)%object > 0 ) then

!           constraint on object

            elgrp = mesh%objects(object)%grpelm_int(node,1,elem)
            elemc = mesh%objects(object)%grpelm_int(node,2,elem)

            if ( elgrp == 0 ) cycle   ! no degrees connected

            call pos_array_elemc ( elgrp, elemc, dofc, posc, physq1, layer1 )

            if ( mesh%objects(object)%typeofobject == 2 .or. &
                 problem%constraints(constr)%object2 > 0 ) then

!             two-sided constraint or connect to other object

              if ( mesh%objects(object)%typeofobject == 2 ) then
!               two-sided constraint
                elgrp2 = mesh%objects(object)%grpelm2_int(node,1,elem)
                elemc2 = mesh%objects(object)%grpelm2_int(node,2,elem)
              else if ( problem%constraints(constr)%object2 > 0 ) then
!               connect to other object
                object2 = problem%constraints(constr)%object2
                elgrp2 = mesh%objects(object2)%grpelm_int(node,1,elem)
                elemc2 = mesh%objects(object2)%grpelm_int(node,2,elem)
              end if

              if ( elgrp2 == 0 ) then   ! no degrees connected
                write(*,'(2(/a)/a,i0/)') &
                  'Error in build_system_constraint: ', &
                  ' connecting node in two-sided object is missing.', &
                  ' are missing. object = ', object2, &
                  ' element = ', elemc2, ' integration point = ', node
                stop
              end if

              call pos_array_elemc ( elgrp2, elemc2, dofc2, posc2, physq2, &
                layer2 )

            else

              dofc2 = 0

            end if

          end if

        else if ( problem%constraints(constr)%discretization == 1 ) then

!         collocation

          if ( problem%constraints(constr)%geometry1 > 0 ) then

!           constraint on geometry

            nodenr = geometry%nodes(node)

            call pos_array_nodenr ( nodenr, dofc, posc, physq1, layer1 )

            if ( problem%constraints(constr)%geometry2 > 0 ) then

!             connection to a second curve

              if ( problem%constraints(constr)%full2 ) then

!               connect to full curve

                call pos_array_fullgeometryc ( geometry2, dofc2, posc2, &
                  physq2, layer2 )

              else

                nodenr2 = geometry2%nodes(node)

!               connect to node of second geometry

                call pos_array_nodenr ( nodenr2, dofc2, posc2, physq2, layer2 )

              end if

            else

              dofc2 = 0

            end if

          else if ( problem%constraints(constr)%elementset1 > 0 ) then

!           constraint on elementset

            nodenr = mesh%elementsets(elementset1)%nodes(node)

            call pos_array_nodenr ( nodenr, dofc, posc, physq1, layer1 )

            if ( problem%constraints(constr)%elementset2 > 0 ) then

!             connection to a second elementset

              nodenr2 = mesh%elementsets(elementset2)%nodes(node)

              call pos_array_nodenr ( nodenr2, dofc2, posc2, physq2, layer2 )

            else

              dofc2 = 0

            end if

          else if ( problem%constraints(constr)%object > 0 ) then

!           constraint on object

            elgrp = mesh%objects(object)%grpelm(node,1)
            elemc = mesh%objects(object)%grpelm(node,2)

            if ( elgrp == 0 ) cycle   ! no degrees connected

            call pos_array_elemc ( elgrp, elemc, dofc, posc, physq1, layer1 )

            if ( mesh%objects(object)%typeofobject == 2 .or. &
                 problem%constraints(constr)%object2 > 0 ) then

!             two-sided constraint or connect to other object

              if ( mesh%objects(object)%typeofobject == 2 ) then
!               two-sided constraint
                elgrp2 = mesh%objects(object)%grpelm2(node,1)
                elemc2 = mesh%objects(object)%grpelm2(node,2)
              else if ( problem%constraints(constr)%object2 > 0 ) then
!               connect to other object
                object2 = problem%constraints(constr)%object2
                elgrp2 = mesh%objects(object2)%grpelm(node,1)
                elemc2 = mesh%objects(object2)%grpelm(node,2)
              end if

              if ( elgrp2 == 0 ) then   ! no degrees connected
                write(*,'(2(/a)/a,i0/)') &
                  'Error in build_system_constraint: ', &
                  ' connecting node in two-sided object is missing.', &
                  ' object nodenr = ', node
                stop
              end if

              call pos_array_elemc ( elgrp2, elemc2, dofc2, posc2, physq2, &
                layer2 )

            else

              dofc2 = 0

            end if

          end if

        end if

!       additional degrees of freedom

        call pos_array_constraint ( problem, constr, dofcadd, poscadd, &
          addunknowns=.true. )

        if ( dofr > ndofr .or. dofc > ndofc .or.  &
             dofc2 > ndofc2 .or. dofcadd > ndofcadd ) then
          write(*,'(/3(a/))') &
            'Internal error in build_system_constraint: ', &
            ' dofr > ndofr .or. dofc > ndofc', &
            ' dofc2 > ndofc2 .or. dofcadd > ndofcadd'
          stop
        end if

!       compute element matrix and vector

        if ( present(elemsub) ) then

          call elemsub ( mesh, problem, constr, elem, node, matrix, vector, &
            first, last, coeffl, oldvl, &
            elemmat(:dofr,:dofc),  elemmat2(:dofr,:dofc2), &
            elemmatadd(:dofr,:dofcadd), elemvec, elemvecadd )

        else if ( present(elemsub1) ) then

          call elemsub1 ( mesh, problem, constr, elem, node, matrix, vector, &
            first, last, coeffl, oldvl, &
            elemmat(:dofr,:dofc), elemmatt(:dofc,:dofr), &
            elemmat2(:dofr,:dofc2), elemmat2t(:dofc2,:dofr), &
            elemmatadd(:dofr,:dofcadd), elemmataddt(:dofcadd,:dofr), &
            elemvecf, elemvecf2, elemvec, elemvecadd )

        else if ( present(elemsub2) ) then

          call elemsub2 ( mesh, problem, constr, elem, node, matrix, vector, &
            first, last, coeffl, oldvl, &
            elemmat(:dofr,:dofc),  elemmat2(:dofr,:dofc2), &
            elemmatdiag(:dofr,:dofr), elemmatadd(:dofr,:dofcadd), &
            elemmatdiagadd(:dofcadd,:dofcadd), elemvec, elemvecadd )

        else if ( present(elemsub3) ) then

          call elemsub3 ( mesh, problem, constr, elem, node, matrix, vector, &
            first, last, coeffl, oldvl, &
            elemmat(:dofr,:dofc), elemmatt(:dofc,:dofr), &
            elemmat2(:dofr,:dofc2), elemmat2t(:dofc2,:dofr), &
            elemmatdiag(:dofr,:dofr), elemmatadd(:dofr,:dofcadd), &
            elemmataddt(:dofcadd,:dofr), elemmatdiagadd(:dofcadd,:dofcadd), &
            elemvecf, elemvecf2, elemvec, elemvecadd )

        else if ( present(elemsub4) ) then

          call elemsub4 ( mesh, problem, constr, elem, node, matrix, vector, &
            first, last, coeffl, oldvl, &
            elemmatsdiag(:dofc,:dofc), elemmatsdiag2(:dofc2,:dofc2), &
            elemmat(:dofr,:dofc),  elemmat2(:dofr,:dofc2), &
            elemmatdiag(:dofr,:dofr), elemmatadd(:dofr,:dofcadd), &
            elemmatdiagadd(:dofcadd,:dofcadd), elemvec, elemvecadd )

        else if ( present(elemsub5) ) then

          call elemsub5 ( mesh, problem, constr, elem, node, matrix, vector, &
            first, last, coeffl, oldvl, &
            elemmatsdiag(:dofc,:dofc), elemmatsdiag2(:dofc2,:dofc2), &
            elemmat(:dofr,:dofc), elemmatt(:dofc,:dofr), &
            elemmat2(:dofr,:dofc2), elemmat2t(:dofc2,:dofr), &
            elemmatdiag(:dofr,:dofr), elemmatadd(:dofr,:dofcadd), &
            elemmataddt(:dofcadd,:dofr), elemmatdiagadd(:dofcadd,:dofcadd), &
            elemvecf, elemvecf2, elemvec, elemvecadd )

        end if

!       transformations

        if ( problem%numtransdegfd > 0 .and. ltransform ) then

          if ( any ( work5(posc) ) ) then

!           extract transformation matrix

            A = extract_submat ( problem%Amat, posc, posc )

!           transform matrix

            if ( matrix ) then

!             elementmatrix * A
              elemmat = fsmatmul ( elemmat, A )

              if ( smat ) then
!               A^T * elementmatrix * A
                elemmatsdiag = fsmatmul ( elemmatsdiag, A )
                elemmatsdiag = stfmatmul ( A, elemmatsdiag )
              end if

              if ( transposedmatf ) then
!               A^T * elementmatrix
                elemmatt = stfmatmul ( A, elemmatt )
              end if

            end if

!           transform vector

            if ( vector .and. transposedmatf ) then
!             A^T * elementvector (possibly multiple right-hand sides)
              do rhsd = 1, nrhsd
                k = dofc*(rhsd-1)
                elemvecf(k+1:k+dofc) = stmatvec ( A, elemvecf(k+1:k+dofc) )
              end do
            end if

            call delete(A)

          end if

          if ( any ( work5(posc2) ) ) then

!           extract transformation matrix

            A = extract_submat ( problem%Amat, posc2, posc2 )

!           transform matrix

            if ( matrix ) then

!             elementmatrix * A
              elemmat2 = fsmatmul ( elemmat2, A )

              if ( smat ) then
!               A^T * elementmatrix * A
                elemmatsdiag2 = fsmatmul ( elemmatsdiag2, A )
                elemmatsdiag2 = stfmatmul ( A, elemmatsdiag2 )
              end if

              if ( transposedmatf ) then
!               A^T * elementmatrix
                elemmat2t = stfmatmul ( A, elemmat2t )
              end if

            end if

!           transform vector

            if ( vector .and. transposedmatf ) then
!             A^T * elementvector (possibly multiple right-hand sides)
              do rhsd = 1, nrhsd
                k = dofc*(rhsd-1)
                elemvecf2(k+1:k+ndofc) = stmatvec ( A, elemvecf2(k+1:k+dofc) )
              end do
            end if

            call delete(A)

          end if

        end if

!       add matrix

        if ( matrix ) then

!         add element matrix row by row

!         A
          call add_elemmat_to_sysmatrix ( sysmatrix, &
            elemmat(:dofr,:dofc), posr(:dofr), posc(:dofc), w1, w2 )
!         A2
          call add_elemmat_to_sysmatrix ( sysmatrix, &
            elemmat2(:dofr,:dofc2), posr(:dofr), posc2(:dofc2), w1, w2 )
!         B
          call add_elemmat_to_sysmatrix ( sysmatrix, &
            elemmatadd(:dofr,:dofcadd), posr(:dofr), poscadd(:dofcadd), w1, w2 )

          if ( transposedmatf ) then

!           "transpose" filled by user

!           A^T
            call add_elemmat_to_sysmatrix ( sysmatrix, &
              elemmatt(:dofc,:dofr), posc(:dofc), posr(:dofr), w1, w2 )
!           A2^T
            call add_elemmat_to_sysmatrix ( sysmatrix, &
              elemmat2t(:dofc2,:dofr), posc2(:dofc2), posr(:dofr), w1, w2 )
!           B^T
            call add_elemmat_to_sysmatrix ( sysmatrix, &
              elemmataddt(:dofcadd,:dofr), poscadd(:dofcadd), &
              posr(:dofr), w1, w2 )

          else

!           add transpose

!           A^T
            call add_elemmat_to_sysmatrix ( sysmatrix, &
              transpose(elemmat(:dofr,:dofc)), posc(:dofc), posr(:dofr), &
              w1, w2 )
!           A2^T
            call add_elemmat_to_sysmatrix ( sysmatrix, &
              transpose(elemmat2(:dofr,:dofc2)), posc2(:dofc2), posr(:dofr), &
              w1, w2 )
!           B^T
            call add_elemmat_to_sysmatrix ( sysmatrix, &
              transpose(elemmatadd(:dofr,:dofcadd)), poscadd(:dofcadd), &
              posr(:dofr), w1, w2 )

          end if

!         add diagonal

          if ( smat ) then

!           diagonal block S1

            call add_elemmat_to_sysmatrix ( sysmatrix, &
              elemmatsdiag(:dofc,:dofc), posc(:dofc), posc(:dofc), w1, w2 )

!           diagonal block S2

            call add_elemmat_to_sysmatrix ( sysmatrix, &
              elemmatsdiag2(:dofc2,:dofc2), posc(:dofc2), posc(:dofc2), w1, w2 )

          end if

          if ( ldiagonal ) then

!           diagonal constraint equation

            diagl(:dofr) = 0
            call add_elemvec_to_sysmatrix ( sysmatrix, &
              diagl(:dofr), posr(:dofr), posr(:dofr), w1, w2 )

          end if

          if ( ldiagonal .or. ldiagonal_addunknowns ) then

!           diagonal additional unknown equation

            diaga(:dofcadd) = 0
            call add_elemvec_to_sysmatrix ( sysmatrix, &
              diaga(:dofcadd), poscadd(:dofcadd), poscadd(:dofcadd), w1, w2 )

          end if

          if ( ldiagonal_block .or. &
                      problem%constraints(constr)%diagonal_block ) then

!           diagonal block constraint equation

            if ( .not. diagmat ) then
              elemmatdiag(:dofr,:dofr) = 0  ! set to zero if not build
            end if
            call add_elemmat_to_sysmatrix ( sysmatrix, &
              elemmatdiag(:dofr,:dofr), posr(:dofr), posr(:dofr), w1, w2 )

          end if

          if ( ldiagonal_block .or. ldiagonal_block_addunknowns .or. &
               problem%constraints(constr)%diagonal_block .or. &
               problem%constraints(constr)%diagonal_block_addunknowns ) then

!           diagonal block additional unknown equation

            if ( .not. diagmat ) then
              elemmatdiagadd(:dofcadd,:dofcadd) = 0  ! set to zero if not build
            end if
            call add_elemmat_to_sysmatrix ( sysmatrix, &
              elemmatdiagadd(:dofcadd,:dofcadd), poscadd(:dofcadd), &
              poscadd(:dofcadd), w1, w2 )

          end if

        end if

!       add vector

        if ( vector ) then

!         add element vector to large vector

          if ( present(msysvector) ) then

!           multiple right-hand side

            do row = 1, dofr

              rowg = posr(row) ! global row number

              do rhsd = 1, nrhsd
                msysvector(rhsd)%u(rowg) = msysvector(rhsd)%u(rowg) &
                                            + elemvec( row + dofr*(rhsd-1) )
              end do

            end do

            do row = 1, dofcadd

              rowg = poscadd(row) ! global row number

              do rhsd = 1, nrhsd
                msysvector(rhsd)%u(rowg) = msysvector(rhsd)%u(rowg) &
                                     + elemvecadd( row + dofcadd*(rhsd-1) )
              end do

            end do

            if ( transposedmatf ) then

              do row = 1, dofc

                rowg = posc(row) ! global row number

                do rhsd = 1, nrhsd
                  msysvector(rhsd)%u(rowg) = msysvector(rhsd)%u(rowg) &
                                              + elemvecf( row + dofc*(rhsd-1) )
                end do

              end do

              do row = 1, dofc2

                rowg = posc2(row) ! global row number

                do rhsd = 1, nrhsd
                  msysvector(rhsd)%u(rowg) = msysvector(rhsd)%u(rowg) &
                                             + elemvecf2( row + dofc2*(rhsd-1) )
                end do

              end do

            end if

          else if ( present(m2sysvector) ) then

!           multiple right-hand side (matrix)

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

            do row = 1, dofcadd

              rowg = poscadd(row) ! global row number

              do j = 1, size(m2sysvector,2)
                do i = 1, si
                  rhsd = i + si*(j-1)
                  m2sysvector(i,j)%u(rowg) = m2sysvector(i,j)%u(rowg) &
                                       + elemvecadd( row + dofcadd*(rhsd-1) )
                end do
              end do

            end do

            if ( transposedmatf ) then

              do row = 1, dofc

                rowg = posc(row) ! global row number

                do j = 1, size(m2sysvector,2)
                  do i = 1, si
                    rhsd = i + si*(j-1)
                    m2sysvector(i,j)%u(rowg) = m2sysvector(i,j)%u(rowg) &
                                               + elemvecf( row + dofc*(rhsd-1) )
                  end do
                end do

              end do

              do row = 1, dofc2

                rowg = posc2(row) ! global row number

                do j = 1, size(m2sysvector,2)
                  do i = 1, si
                    rhsd = i + si*(j-1)
                    m2sysvector(i,j)%u(rowg) = m2sysvector(i,j)%u(rowg) &
                                             + elemvecf2( row + dofc2*(rhsd-1) )
                  end do
                end do

              end do

            end if

          else if ( present(m3sysvector) ) then

!           multiple right-hand side (3d array)

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

            do row = 1, dofcadd

              rowg = poscadd(row) ! global row number

              do k = 1, size(m3sysvector,3)
                do j = 1, sj
                  do i = 1, si
                    rhsd = i + si*(j-1) + si*sj*(k-1)
                    m3sysvector(i,j,k)%u(rowg) = m3sysvector(i,j,k)%u(rowg) &
                                       + elemvecadd( row + dofcadd*(rhsd-1) )
                  end do
                end do
              end do

            end do

            if ( transposedmatf ) then

              do row = 1, dofc

                rowg = posc(row) ! global row number

                do k = 1, size(m3sysvector,3)
                  do j = 1, sj
                    do i = 1, si
                      rhsd = i + si*(j-1) + si*sj*(k-1)
                      m3sysvector(i,j,k)%u(rowg) = m3sysvector(i,j,k)%u(rowg) &
                                               + elemvecf( row + dofc*(rhsd-1) )
                    end do
                  end do
                end do

              end do

              do row = 1, dofc2

                rowg = posc2(row) ! global row number

                do k = 1, size(m3sysvector,3)
                  do j = 1, sj
                    do i = 1, si
                      rhsd = i + si*(j-1) + si*sj*(k-1)
                      m3sysvector(i,j,k)%u(rowg) = m3sysvector(i,j,k)%u(rowg) &
                                             + elemvecf2( row + dofc2*(rhsd-1) )
                    end do
                  end do
                end do

              end do

            end if

          else

!           single right-hand side

            do row = 1, dofr

              rowg = posr(row) ! global row number

              sysvector%u(rowg) = sysvector%u(rowg) + elemvec(row)

            end do

            do row = 1, dofcadd

              rowg = poscadd(row) ! global row number

              sysvector%u(rowg) = sysvector%u(rowg) + elemvecadd(row)

            end do

            if ( transposedmatf ) then

              do row = 1, dofc

                rowg = posc(row) ! global row number

                sysvector%u(rowg) = sysvector%u(rowg) + elemvecf(row)

              end do

              do row = 1, dofc2

                rowg = posc2(row) ! global row number

                sysvector%u(rowg) = sysvector%u(rowg) + elemvecf2(row)

              end do

            end if

          end if

        end if

      end do

      deallocate ( elemmat, elemmat2, elemmatadd, elemvec, elemvecadd )
      deallocate ( posr, posc, posc2, poscadd )

      if ( smat ) then
        deallocate ( elemmatsdiag, elemmatsdiag2 )
      end if

      if ( ldiagonal .or. ldiagonal_addunknowns ) deallocate ( diagl, diaga )

      if ( ldiagonal_block .or. ldiagonal_block_addunknowns .or. diagmat .or. &
           problem%constraints(constr)%diagonal_block .or. &
           problem%constraints(constr)%diagonal_block_addunknowns ) then
        deallocate ( elemmatdiag, elemmatdiagadd )
      end if

      if ( transposedmatf ) then
        deallocate ( elemmatt, elemmat2t )
        deallocate ( elemmataddt )
        deallocate ( elemvecf, elemvecf2 )
      end if

      if ( problem%constraints(constr)%discretization == 0 ) then
!       weak constraint
        if ( problem%constraints(constr)%elementset1 > 0 ) then
!         weak constraint on elementset
          deallocate ( work1, work2 )
          if ( problem%constraints(constr)%elementset2 > 0 ) then
!           connected elementset
            deallocate ( work3, work4 )
          end if
        end if
      end if

!     remove memory from geometry (only necessary for typegeometry= 1 or 5)

      if ( allocated(geometry%nodes) ) deallocate(geometry%nodes)
      if ( allocated(geometry2%nodes) ) deallocate(geometry2%nodes)

    end do

    deallocate( w1, w2 )

    deallocate ( lp )

    if ( problem%numtransdegfd > 0. .and. ltransform ) then
      deallocate ( work5 )
    end if

  contains


!   determine maximum number of degrees on the geometry (one element)

    subroutine ndof_geometry_elem ( ndof, geometry, physq )

      integer, intent(out) :: ndof
      integer, intent(in)  :: physq
      type(geometry_t), intent(in) :: geometry

      integer :: elem, node, nodenr, dofe

      ndof = 0

      if ( physq > 0 ) then
        do elem = 1, geometry%nelem
          dofe = 0
          do node = 1, geometry%elnumnod
            nodenr = geometry%topology(node,elem,2)
            dofe = dofe + problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physq)) &
                        - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physq))
          end do
          ndof = max( ndof, dofe )
        end do
      else
        do elem = 1, geometry%nelem
          dofe = 0
          do node = 1, geometry%elnumnod
            nodenr = geometry%topology(node,elem,2)
            dofe = dofe + problem%nodnumdegfd(nodenr+1) &
                        - problem%nodnumdegfd(nodenr)
          end do
          ndof = max( ndof, dofe )
        end do
      end if

    end subroutine ndof_geometry_elem


!   determine number of degrees on the curve (in nodes)

    subroutine ndof_geometry_node ( ndof, geometry, physq )

      integer, intent(out) :: ndof
      integer, intent(in)  :: physq
      type(geometry_t), intent(in) :: geometry

      integer :: node, nodenr

      ndof = 0

      if ( physq > 0 ) then
        do node = 1, geometry%nnodes
          nodenr = geometry%nodes(node)
          ndof = max( ndof, problem%vec_nodnumdegfd(nodenr+1,&
                                    &problem%physq(physq)) &
                      - problem%vec_nodnumdegfd(nodenr,&
                                    &problem%physq(physq)) )
        end do
      else
        do node = 1, geometry%nnodes
          nodenr = geometry%nodes(node)
          ndof = max( ndof, problem%nodnumdegfd(nodenr+1) &
                      - problem%nodnumdegfd(nodenr) )
        end do
      end if

    end subroutine ndof_geometry_node


!   determine maximum number of degrees in an element

    subroutine ndof_elem ( ndof, physq )

      integer, intent(out) :: ndof
      integer, intent(in)  :: physq

      integer :: elgrp

      ndof = 0

      do elgrp = 1, mesh%nelgrp

        if ( physq > 0 ) then
          ndof = max( ndof, &
                       sum( problem%vec_elnumdegfd(elgrp)%a(:,physq) ) )
        else
          ndof = max( ndof, sum( problem%elnumdegfd(elgrp)%a ) )
        end if

      end do

      if ( problem%numlayers > 0 ) ndof = ndof * problem%maxnodnumlayers

    end subroutine ndof_elem


!   determine number of degrees on the curve (full curve)

    subroutine ndof_fullgeometry ( ndof, geometry, physq )

      integer, intent(out) :: ndof
      integer, intent(in)  :: physq
      type(geometry_t), intent(in) :: geometry

      integer :: node, nodenr

      ndof = 0

      if ( physq > 0 ) then
        do node = 1, geometry%nnodes
          nodenr = geometry%nodes(node)
          ndof = ndof + problem%vec_nodnumdegfd(nodenr+1,&
                                    &problem%physq(physq)) &
                      - problem%vec_nodnumdegfd(nodenr,&
                                    &problem%physq(physq))
        end do
      else
        do node = 1, geometry%nnodes
          nodenr = geometry%nodes(node)
          ndof = ndof + problem%nodnumdegfd(nodenr+1) &
                      - problem%nodnumdegfd(nodenr)
        end do
      end if

    end subroutine ndof_fullgeometry


!   determine number of degrees on the elementset (in nodes)

    subroutine ndof_elementset_node ( ndof, elementset, physq )

      integer, intent(out) :: ndof
      integer, intent(in)  :: physq
      type(elementset_t), intent(in) :: elementset

      integer :: node, nodenr

      ndof = 0

      if ( physq > 0 ) then
        do node = 1, elementset%nnodes
          nodenr = elementset%nodes(node)
          ndof = max( ndof, problem%vec_nodnumdegfd(nodenr+1,&
                                    &problem%physq(physq)) &
                      - problem%vec_nodnumdegfd(nodenr,&
                                    &problem%physq(physq)) )
        end do
      else
        do node = 1, elementset%nnodes
          nodenr = elementset%nodes(node)
          ndof = max( ndof, problem%nodnumdegfd(nodenr+1) &
                      - problem%nodnumdegfd(nodenr) )
        end do
      end if

    end subroutine ndof_elementset_node


!   check whether all nodes in layer

    subroutine check_layer ( layer, lp )

      integer, intent(in) :: layer
      logical, dimension(:), intent(in) :: lp

      if ( .not. all( lp ) ) then
        if ( problem%constraints(constr)%errorlayer == 0 ) then
          write(*,'(3(/a)/2(a,i0)/)') &
            'Error in build_system_constraint: ', &
            ' Constraint not fully within layer.', &
            ' Stop because errorlayer = 0 ', &
            ' constraint = ', constr, ' layer = ', layer
          stop
        else if ( problem%constraints(constr)%errorlayer == 1 ) then
          write(*,'(3(/a)/2(a,i0)/)') &
            'Warning in build_system_constraint: ', &
            ' Constraint not fully within layer.', &
            ' Warning because errorlayer = 1 ', &
            ' constraint = ', constr, ' layer = ', layer
        end if
      end if

    end subroutine check_layer


!   determine position array for element elemc

    subroutine pos_array_elemc ( elgrp, elemc, dofc, posc, physq, layer )

      integer, intent(in) :: elgrp, elemc, physq, layer
      integer, intent(out) :: dofc
      integer, dimension(:), intent(out) :: posc

      if ( physq > 0 ) then
!       physical quantity specified
        call pos_array ( mesh, problem, elgrp, elemc, dofc, posc, &
          [ physq ], order, layer, lp )
      else
        call pos_array ( mesh, problem, elgrp, elemc, dofc, posc, &
          order=order, layer=layer, lp=lp )
      end if

!     check full layer
      if ( layer > 0 ) then
        call check_layer ( layer, lp(1:mesh%elnumnod(elgrp)) )
      end if

    end subroutine pos_array_elemc


!   determine position array for node nodenr

    subroutine pos_array_nodenr ( nodenr, dofc, posc, physq, layer )

      integer, intent(in) :: nodenr, physq, layer
      integer, intent(out) :: dofc
      integer, dimension(:), intent(out) :: posc

      if ( physq > 0 ) then
!       physical quantity specified
        call pos_array_node ( problem, nodenr, dofc, posc, [ physq ], &
          layer )
      else
        call pos_array_node ( problem, nodenr, dofc, posc, layer=layer )
      end if

!     check layer
      if ( layer > 0 ) then
        lp(1) = btest(problem%nodlayers(nodenr),layer-1)
        call check_layer ( layer, lp(1:1) )
      end if

    end subroutine pos_array_nodenr


!   determine position array for geometry

    subroutine pos_array_geometryc ( geometry, dofc, posc, physq, layer )

      type(geometry_t), intent(in) :: geometry
      integer, intent(in) :: physq, layer
      integer, intent(out) :: dofc
      integer, dimension(:), intent(out) :: posc

      if ( physq > 0 ) then
!       physical quantity specified
        call pos_array_geometry ( problem, geometry, elem, dofc, &
          posc, [ physq ], order, layer, lp )
      else
        call pos_array_geometry ( problem, geometry, elem, dofc, posc, &
          order=order, layer=layer, lp=lp )
      end if

!     check full layer
      if ( layer > 0 ) then
          call check_layer ( layer, lp(1:geometry%elnumnod) )
      end if

    end subroutine pos_array_geometryc


!   determine position array for full geometry

    subroutine pos_array_fullgeometryc ( geometry, dofc, posc, physq, layer )

      type(geometry_t), intent(in) :: geometry
      integer, intent(in) :: physq, layer
      integer, intent(out) :: dofc
      integer, dimension(:), intent(out) :: posc

      if ( physq > 0 ) then
!       physical quantity specified
        call pos_array_fullgeometry ( problem, geometry, dofc, &
          posc, [ physq ], order, layer, lp )
      else
        call pos_array_fullgeometry ( problem, geometry, dofc, &
          posc, order=order, layer=layer, lp=lp )
      end if
!     check full layer
      if ( layer > 0 ) then
        call check_layer ( layer, lp(1:geometry%nnodes) )
      end if

    end subroutine pos_array_fullgeometryc

  end subroutine build_system_constraint


! Assemble system matrix and system vector for additional unknowns in
! a single constraint

  subroutine build_system_constraint_addunknowns ( mesh, problem, sysmatrix, &
    sysvector, msysvector, constraint, addmatvec, addmat, addvec, matdiag, &
    matdiagblock, rhs )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the system to be assembled
    type(sysmatrix_t), intent(inout) :: sysmatrix

!   the right-hand side vector to be assembled
    type(sysvector_t), intent(inout), optional :: sysvector

!   multiple right-hand side vectors to be assembled
    type(sysvector_t), dimension(:), intent(inout), optional :: msysvector

!   the constraint to assemble the additional unknown equation for
    integer, intent(in) :: constraint

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

!   if present: the diagonal of the matrix
    real(dp), dimension(:), intent(in), optional :: matdiag

!   if present: the full diagonal block matrix
    real(dp), dimension(:,:), intent(in), optional :: matdiagblock

!   if present: the right-hand side
    real(dp), dimension(:), intent(in), optional :: rhs

!   The original equations are
!
!      S u = f
!
!   The constraints are
!
!      A u + B a = g
!
!   where a are the additional unknowns.
!
!   The structure of the matrix system now becomes:
!
!      [  S  A^T  0 ] [ u ] = [ f ]
!      [  A   0   B ] [ l ] = [ g ]
!      [  0  B^T  0 ] [ a ] = [ h ]
!
!   where l are the Lagrange multipliers and h is a right-hand side for the
!   additional unknowns to introduce forcing on the additional unknowns.
!   This system can be build using build_system_constraint.
!
!   In the current routine (build_system_constraint_addunknowns), the lower
!   right diagonal zero block regarding the additional unknowns can be filled:
!
!      [  S  A^T  0 ] [ u ] = [ f ]
!      [  A   0   B ] [ l ] = [ g ]
!      [  0  B^T  C ] [ a ] = [ h + r ]
!
!   where C and r are assembled in this routine. The matrix C can be either
!   diagonal (matdiag) or full for a single constraint (matdiagblock). The
!   additional right-hand side r (rhs) can be used to put additional forces
!   on the additional unknowns.
!
!   NOTE: create_sysmatrix_structure_constraint needs to be called with
!         diagonal_addunknowns=.true. if matdiag is present, except if
!         diagonal=.true. is already present.
!
!   NOTE: create_sysmatrix_structure_constraint needs to be called with
!         diagonal_block_addunknowns=.true. if matblockdiag is present.
!

    logical :: matrix, vector, laddmat, laddvec
    integer :: ndofadd, numund, numess, row, rowg
    integer :: dofadd, nrhsd, rhsd
    integer, allocatable, dimension(:) :: posadd, w1, w2


    if ( problem%numconstraints == 0 ) return

    call check ( mesh, 'build_system_constraint_addunknowns' )
    call check ( problem, 'build_system_constraint_addunknowns', mesh )

!   initialize local parameters

    nrhsd = 1

    if ( constraint < 1 .or. constraint > problem%numconstraints ) then

      write(*,'(/2(a/),2(a,i0/))') &
        'Error: constraint in the heading of', &
        ' build_system_constraint_addunknowns is out of range ', &
        ' constraint is ', constraint, &
        ' whereas the number of constraints is ', problem%numconstraints
      stop

    end if

!   compute matrix and vector?

    matrix = present(matdiag) .or. present(matdiagblock)

    vector = present(rhs)

    if ( .not. ( matrix .or. vector ) ) return

!   clear matrix and vector?

    if ( matrix ) then

      if ( .not. sysmatrix%initialized_structure ) then
        write(*,'(/2(a/))') &
          'Error in build_system_constraint_addunknowns:', &
          ' no system matrix structure.'
        stop
      end if

      if ( .not. sysmatrix%finalized ) then
        write(*,'(/2a/)') &
          'Error in build_system_constraint_addunknowns:', &
          ' system matrix has not been finalized.'
        stop
      end if

      if ( .not. sysmatrix%allocated_data ) then
        write(*,'(/2a/)') &
          'Error in build_system_constraint_addunknowns:', &
          ' data in system matrix not allocated.'
        stop
      end if

      if ( sysmatrix%Suu%m + sysmatrix%Sup%m /= problem%numdegfd ) then
        write(*,'(/3(a/))') &
          'Error in build_system_constraint_addunknowns: ', &
          ' number of degrees of freedom of the system matrix is', &
          ' is different from the number in problem.'
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

      if ( present(sysvector) .and. present(msysvector) ) then

        write(*,'(/2a/)') &
          'Error in build_system_constraint_addunknowns:', &
          ' sysvector and msysvector both present'
        stop

      else if ( present(sysvector) ) then

!       single right-hand side

        if ( .not. sysvector%created ) then
          write(*,'(/2(a/))') &
            'Error in build_system_constraint_addunknowns: ', &
            ' sysvector not created.'
          stop
        end if

        if ( sysvector%n /= problem%numdegfd ) then
          write(*,'(/3(a/))') &
            'Error in build_system_constraint_addunknowns: ', &
            ' number of degrees of freedom of the sysvector vector ', &
            'different from the number in problem.'
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
          write(*,'(/2(a/))') &
            'Error in build_system_constraint_addunknowns: ', &
            ' msysvector not created.'
          stop
        end if

        if ( any( msysvector%n /= problem%numdegfd ) ) then
          write(*,'(/3(a/))') &
            'Error in build_system_constraint_addunknowns: ', &
            ' number of degrees of  freedom of the ', &
            'msysvector vector different from the number in problem.'
          stop
        end if

        if ( .not. laddvec ) then

!         clear vector
          do rhsd = 1, nrhsd
            msysvector(rhsd)%u = 0
          end do

        end if

      else

        write(*,'(/2(a/))') &
          'Error in build_system_constraint_addunknowns:', &
          ' no sysvector or msysvector present'
        stop

      end if

    end if

    ndofadd = problem%constraints(constraint)%naddunknowns

!   test arrays

    if ( present(matdiag) ) then
      if ( size(matdiag) /= ndofadd ) then
        write(*,'(/2(a/))') &
          'Error in build_system_constraint_addunknowns:', &
          ' size of matdiag /= number of additional unknowns '
        stop
      end if
    end if

    if ( present(matdiagblock) ) then
      if ( any( shape(matdiagblock) /= [ ndofadd, ndofadd ] ) ) then
        write(*,'(/2(a/))') &
          'Error in build_system_constraint_addunknowns:', &
          ' shape of matdiagblock is /= (/ ndofadd, ndofadd /) '
        stop
      end if
    end if

    if ( present(rhs) ) then
      if ( size(rhs) /= ndofadd * nrhsd ) then
        write(*,'(/2(a/))') &
          'Error in build_system_constraint_addunknowns:', &
          ' size of rhs /= number of additional unknowns * nrhs'
        stop
      end if
    end if

!   temporary work arrays to store information on current large matrix row

    numess = problem%numessdegfd
    numund = problem%numundegfd

    allocate( w1(numund), w2(numess) )

    w1 = 0
    w2 = 0

!   reserve memory positions

    allocate ( posadd(ndofadd) )

!   additional degrees of freedom

    call pos_array_constraint ( problem, constraint, dofadd, posadd, &
      addunknowns=.true. )

    if ( dofadd /= ndofadd ) then
      write(*,'(/3(a/))') &
        'Internal error in build_system_constraint_addunknowns: ', &
        ' dofadd > ndofadd'
      stop
    end if

!   add diagonal block matrix row by row

    if ( present(matdiagblock) ) then

      call add_elemmat_to_sysmatrix ( sysmatrix, &
        matdiagblock, posadd, posadd, w1, w2 )

    end if

!   add diagonal matrix

    if ( present(matdiag) ) then

      call add_elemvec_to_sysmatrix ( sysmatrix, &
        matdiag, posadd, posadd, w1, w2 )

    end if

    if ( vector ) then

!     add element vector to large vector

      if ( present(msysvector) ) then

!       multiple right-hand side

        do row = 1, ndofadd

          rowg = posadd(row) ! global row number

          do rhsd = 1, nrhsd
            msysvector(rhsd)%u(rowg) = msysvector(rhsd)%u(rowg) &
                                 + rhs( row + ndofadd*(rhsd-1) )
          end do

        end do

      else

!       single right-hand side

        do row = 1, ndofadd

          rowg = posadd(row) ! global row number

          sysvector%u(rowg) = sysvector%u(rowg) + rhs(row)

        end do

      end if

      deallocate ( posadd )

    end if

    deallocate( w1, w2 )

  end subroutine build_system_constraint_addunknowns


! Get element/nodal Lagrange multipliers or additional unknowns + index from
! the sysvector.

  subroutine get_sysvector_constraint ( mesh, problem, sysvector, constraint, &
    elem, node, addunknowns, u, order, posu, ndofu )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(sysvector_t), intent(in) :: sysvector

!   the constraint number
    integer, intent(in) :: constraint

!   if present the element number in the geometry or object
    integer, intent(in), optional :: elem

!   if present the node number in the geometry or object
    integer, intent(in), optional :: node

!   if .true. the degrees of freedom are the additional unknowns
    logical, intent(in), optional :: addunknowns

!   the element/nodal degrees of freedom, if present
    real(dp), intent(out), dimension(:), optional :: u

!   The parameter order determines the sequence of the degrees of freedom on
!   elementlevel. There are two possibilities:
!     order = 'ND' : the most inner loop is over the degrees of freedom
!     order = 'DN' : the most inner loop is over the nodal points
!   this only makes sense for the elem parameter
!   default = 'DN'
    character(len=*), intent(in), optional :: order

!   if present it contains the positions of the components of u in sysvector
    integer, dimension(:), intent(out), optional :: posu

!   if present it gives the number of degrees of freedom stored in u
    integer, intent(out), optional :: ndofu


    integer :: ndof, lpos, object, nodeset, point
    integer, dimension(:), allocatable :: pos
    logical :: laddunknowns
    type(geometry_t) :: geometry


    call check ( problem, 'get_sysvector_constraint', mesh )

!   Some testing

    if ( .not. sysvector%created ) then
      write(*,'(/a/)') &
        'Error in get_sysvector_constraint: sysvector not created.'
      stop
    end if

    if ( problem%probnr /= sysvector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in get_sysvector_constraint: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in sysvector = ',  sysvector%probnr
      stop
    end if

!   ordering; set defaults

    if ( present(order) ) then
      if ( all ( order /= [ 'ND', 'DN' ] ) ) then
        write(*,'(2(/a)/)') &
          'Error in get_sysvector_constraint: ', &
          ' heading parameter order must be either ''ND'' or ''DN''.'
        stop
      end if
    end if

    if ( constraint < 1 .or. constraint > problem%numconstraints ) then

      write(*,'(/2(a/),2(a,i0),/a,i0/)') &
        'Error: constraint in the heading of', &
        ' get_sysvector_constraint is out of range ', &
        ' constraint is ', constraint, &
        ' whereas the number of constraints is ', problem%numconstraints
      stop

    end if

    if ( present(elem) ) then
      if ( present(node) .or. present(addunknowns) ) then
        write(*,'(/2(a/))') &
          'Error get_sysvector_constraint: ', &
          ' if elem present, node or addunknowns can not be present'
        stop
      end if
      if ( problem%constraints(constraint)%discretization == 1 ) then
        write(*,'(/2(a/),a,i0/)') &
          'Error get_sysvector_constraint: ', &
          ' elem in heading, whereas the constraint is using collocation', &
          ' constraint = ', constraint
        stop
      end if
      object = problem%constraints(constraint)%object
      if ( object > 0 ) then
        if ( .not. mesh%objects(object)%topol ) then
          write(*,'(/a/,a,i0/)') &
            'Error get_sysvector_constraint: object has no topology ', &
            ' constraint = ', constraint
          stop
        end if
      end if
      if ( problem%constraints(constraint)%typeconstraint == 2 ) then
        write(*,'(/2(a/),a,i0/)') &
          'Error get_sysvector_constraint: ', &
          ' elem in heading, whereas the constraint is global', &
          ' constraint = ', constraint
        stop
      end if
    else if ( present(node) ) then
      if ( present(addunknowns) ) then
        write(*,'(/2(a/))') &
          'Error get_sysvector_constraint: ', &
          ' if node in heading, addunknowns can not be present'
        stop
      end if
      if ( problem%constraints(constraint)%typeconstraint == 2 ) then
        write(*,'(/2(a/),a,i0/)') &
          'Error get_sysvector_constraint: ', &
          ' node in heading, whereas the constraint is global', &
          ' constraint = ', constraint
        stop
      end if
    end if

    laddunknowns = set_optional ( variable=addunknowns, default=.false. )

    lpos = max(problem%constraints(constraint)%maxnumdegfd,&
               problem%constraints(constraint)%nglobalc,&
               problem%constraints(constraint)%naddunknowns)
    if ( problem%constraints(constraint)%typeconstraint == 1 .and. &
         problem%constraints(constraint)%discretization == 0 ) then
!     weak distributed
      lpos = max(lpos,sum(problem%constraints(constraint)%elnumdegfd))
    end if

    allocate(pos(lpos))

!   get positions

    if ( problem%constraints(constraint)%typeconstraint == 0 ) then

!     not active

      ndof = 0

      write(*,'(/2(a/),a,i0/)') &
        'Warning: constraint in the heading of', &
        ' get_sysvector_constraint is not active ', &
        ' constraint number = ', constraint

    else if ( problem%constraints(constraint)%typeconstraint == 1 ) then

!     distributed

      if ( problem%constraints(constraint)%geometry1 > 0 ) then
!       on geometry
        if ( problem%constraints(constraint)%typegeometry == 1 ) then
!         constraint in point; fill fake geometry
          point = problem%constraints(constraint)%geometry1
          geometry%nnodes = 1
          geometry%nodes = [mesh%points(point)]
        else if ( problem%constraints(constraint)%typegeometry == 2 ) then
!         curve
          geometry = mesh%curves(problem%constraints(constraint)%geometry1)
        else if ( problem%constraints(constraint)%typegeometry == 3 ) then
!         surface
          geometry = mesh%surfaces(problem%constraints(constraint)%geometry1)
        else if ( problem%constraints(constraint)%typegeometry == 4 ) then
!         volume
          geometry = mesh%volumes(problem%constraints(constraint)%geometry1)
        else if ( problem%constraints(constraint)%typegeometry == 5 ) then
!         constraint on nodeset; fill fake geometry
          nodeset = problem%constraints(constraint)%geometry1
          geometry%nnodes = size(mesh%nodesets(nodeset)%a)
          geometry%nodes = mesh%nodesets(nodeset)%a
        end if
        call pos_array_constraint ( problem, constraint, ndof, pos, geometry, &
          elem=elem, node=node, addunknowns=addunknowns, order=order )
      else if ( problem%constraints(constraint)%object > 0 ) then
!       on object
        call pos_array_constraint ( problem, constraint, ndof, pos, &
          object=mesh%objects(problem%constraints(constraint)%object), &
          elem=elem, node=node, addunknowns=addunknowns, order=order )
      end if

    else if ( problem%constraints(constraint)%typeconstraint == 2 ) then

!     global

      if ( laddunknowns ) then
        call pos_array_constraint ( problem, constraint, ndof, pos, &
          addunknowns=.true. )
      else
        call pos_array_constraint ( problem, constraint, ndof, pos, &
          globalc=.true. )
      end if

    end if

    if ( present(ndofu) ) ndofu = ndof

!   check u array

    if ( present(u) ) then
      if ( size(u) < ndof ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_sysvector_constraint: ', &
          ' u array is too small: size(u) = ', size(u), &
          ' whereas the number of degrees of freedom = ', ndof
        stop
      end if
    end if

!   check posu array

    if ( present(posu) ) then
      if ( size(posu) < ndof ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_sysvector_constraint: ', &
          ' posu array is too small: size(posu) = ', size(posu), &
          ' whereas the number of degrees of freedom = ', ndof
        stop
      end if
    end if

    if ( lpos < ndof ) then
      write(*,'(/a/a/)') &
        'Internal error in get_sysvector_constraint: ', &
        ' lpos < ndof '
      stop
    end if

!   compute u

    if ( present(u) ) u(1:ndof) = sysvector%u(pos(1:ndof))

!   posu array

    if ( present(posu) ) posu(1:ndof) = pos(1:ndof)

    deallocate(pos)

  end subroutine get_sysvector_constraint


! create vector subscript array in sysvector for constraints

  subroutine create_subscript_constraint ( mesh, problem, subscript, &
    constraints, lagrangian_multipliers, addunknowns, degfd, unknownpart, &
    essentialpart, fillnodes )

    type(mesh_t), intent(in)  :: mesh
    type(problem_t), intent(in)  :: problem

!   the vector subscript of the degrees of freedom in a sysvector (constraints)
    type(subscriptcon_t), intent(inout) :: subscript

!   constraints included for defining the subscript.
!   For example constraints=(/2,4/), includes constraints 2 and 4.
!   Default is to include all constraints.
    integer, dimension(:), intent(in), optional :: constraints

!   include the Lagrange multipliers in the subscript
!   default = .true.
    logical, intent(in), optional :: lagrangian_multipliers

!   include the additional unknowns in the subscript
!   default = .true.
    logical, intent(in), optional :: addunknowns

!   if degfd is present only the degfd^th degree of freedom will be included.
!   NOTE: this applies to distributed and global constraints as well as the
!   the additional unknowns!
    integer, intent(in), optional :: degfd

!   logicals to indicate inclusion of the unknown degrees and the essential
!   degrees of freedom. Defaults: unknownpart=.true., essentialpart=.true.
    logical, intent(in), optional :: unknownpart, essentialpart

!   logical to indicate whether the nodes involved in the subscript should
!   put into the subscript%nodes.
!   default = .false.
!   NOTE: only nodes that have an non-zero number of degrees of freedom in
!   the subcript are stored.
!   NOTE: this applies only to distributed degrees of freedom.
    logical, intent(in), optional :: fillnodes


    integer :: dof, constr, geomnr, nnodes, nodenr, bp, nndof, ip
    integer :: nsnodes
    logical :: unpart, esspart, lams, adduns, lconstr(problem%numconstraints)
    logical :: node_with_degrees, lfillnodes
    integer, allocatable, dimension(:) :: work, snodes


    call check ( mesh, 'create_subscript_constraint' )
    call check ( problem, 'create_subscript_constraint' )

!   subscript already allocated?
    if (allocated(subscript%s)) deallocate(subscript%s)
    if (allocated(subscript%nodes)) deallocate(subscript%nodes)

!   include unknown part?
    unpart = set_optional ( variable=unknownpart, default=.true. )

!   include essential part?
    esspart = set_optional ( variable=essentialpart, default=.true. )

!   include lagrangian multipliers?
    lams = set_optional ( variable=lagrangian_multipliers, default=.true. )

!   include additional unknowns?
    adduns = set_optional ( variable=addunknowns, default=.true. )

!   include nodes?
    lfillnodes = set_optional ( variable=fillnodes, default=.false. )


!   check constraints
    if ( present(constraints) ) then
!     specified constraints
      if ( any ( constraints < 1 ) .or. &
           any ( constraints > problem%numconstraints ) ) then
        write(*,'(/a/a/a,i0/)') &
          'Error: parameter constraints in the heading of', &
          ' create_subscript_constraint is ', &
          ' out of range: some constraints are < 1 or larger than ', &
           problem%numconstraints
        stop
      end if
      lconstr = .false.
      lconstr(constraints) = .true.
    else
!     all constraints
      lconstr = .true.
    end if

    allocate ( work(problem%numconstrdegfd), snodes(problem%numconstrdegfd) )

    dof = 0

    nsnodes = 0

    do constr = 1, problem%numconstraints

      if ( .not. lconstr(constr) ) cycle  ! do not include constr

!     get positions in array

      if ( problem%constraints(constr)%typeconstraint == 0 ) then

!       not active

        write(*,'(/2(a/),a,i0/)') &
          'Warning create_subscript_constraint: ', &
          ' constraint is not active ', &
          ' constraint number = ', constr

        cycle

      end if

      if ( lams ) then

!       include lagrangian multipliers

        if ( problem%constraints(constr)%typeconstraint == 1 ) then

!         distributed

          if ( problem%constraints(constr)%geometry1 > 0 ) then
!           on geometry
            geomnr = problem%constraints(constr)%geometry1
            if ( problem%constraints(constr)%typegeometry == 2 ) then
!             curve
              nnodes = mesh%curves(geomnr)%nnodes
            else if ( problem%constraints(constr)%typegeometry == 3 ) then
!             surface
              nnodes = mesh%surfaces(geomnr)%nnodes
            else if ( problem%constraints(constr)%typegeometry == 4 ) then
!             volume
              nnodes = mesh%volumes(geomnr)%nnodes
            end if
          else if ( problem%constraints(constr)%object > 0 ) then
!           on object
            nnodes = mesh%objects(problem%constraints(constr)%object)%nnodes
          end if

!         scan nodes for degrees of freedom

          do nodenr = 1, nnodes

!           positions in nodes

            bp = problem%constraints(constr)%nodnumdegfd(nodenr)
            nndof = problem%constraints(constr)%nodnumdegfd(nodenr+1) &
                       - problem%constraints(constr)%nodnumdegfd(nodenr)

!           positions in sysvector and sysmatrix (renumbered)

            call fill_work

            if ( lfillnodes .and. node_with_degrees ) then
              nsnodes = nsnodes + 1
              snodes(nsnodes) = nodenr
            end if

          end do

        else if ( problem%constraints(constr)%typeconstraint == 2 ) then

!         global

          bp = problem%constraints(constr)%globnumdegfd(1)
          nndof = problem%constraints(constr)%globnumdegfd(2) &
                     - problem%constraints(constr)%globnumdegfd(1)

!         positions in sysvector and sysmatrix (renumbered)

          call fill_work

        end if

      end if

      if ( adduns ) then

!       include additional unknowns

        bp = problem%constraints(constr)%addnumdegfd(1)
        nndof = problem%constraints(constr)%addnumdegfd(2) &
                   - problem%constraints(constr)%addnumdegfd(1)

!       positions in sysvector and sysmatrix (renumbered)

        call fill_work

      end if

    end do

    if ( dof == 0 ) then
      write(*,'(/a/)') &
        'Warning in create_subscript_constraint: subscript is empty '
    end if

    subscript%s = work(1:dof)

    if ( lfillnodes ) then
      subscript%nodes = snodes(1:nsnodes)
    end if

    deallocate ( work, snodes )

  contains

    subroutine fill_work

      integer :: deg

      node_with_degrees = .false.

      if ( present(degfd) ) then

!       one degree

        if ( degfd >=1 .and. degfd <= nndof ) then

!         degfd in valid range

          ip =  problem%degfdperm(bp+degfd,2)

          if ( ip <= problem%numundegfd .and. unpart .or. &
               ip >  problem%numundegfd .and. esspart ) then

!           position in sysvector and sysmatrix (renumbered)
            work(dof+1) = ip

            dof = dof + 1

            node_with_degrees = .true.

          end if

        end if

      else

!       all degrees

        do deg = 1, nndof

          ip =  problem%degfdperm(bp+deg,2)

          if ( ip <= problem%numundegfd .and. unpart .or. &
               ip >  problem%numundegfd .and. esspart ) then

!           position in sysvector and sysmatrix (renumbered)
            work(dof+1) = ip

            dof = dof + 1

            node_with_degrees = .true.

          end if

        end do

      end if

    end subroutine fill_work

  end subroutine create_subscript_constraint


! delete vector subscript array in sysvector for constraints

  subroutine delete_single_subscript_constraint ( subscript )

!   the vector subscript of the degrees of freedom in a sysvector
    type(subscriptcon_t), intent(out) :: subscript

  end subroutine delete_single_subscript_constraint


! delete vector subscript array in sysvector for constraints

  subroutine delete_subscript_constraint ( subscript1, subscript2, subscript3, &
    subscript4, subscript5 )

!   the vector subscript of the degrees of freedom in a sysvector
    type(subscriptcon_t), intent(inout) :: subscript1
    type(subscriptcon_t), intent(inout), optional :: subscript2, subscript3, &
      subscript4, subscript5

    call delete_single_subscript_constraint(subscript1)
    if ( present(subscript2) ) &
                  call delete_single_subscript_constraint(subscript2)
    if ( present(subscript3) ) &
                  call delete_single_subscript_constraint(subscript3)
    if ( present(subscript4) ) &
                  call delete_single_subscript_constraint(subscript4)
    if ( present(subscript5) ) &
                  call delete_single_subscript_constraint(subscript5)

  end subroutine delete_subscript_constraint

end module system_constraint_m
