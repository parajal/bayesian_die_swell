
! Copyright (C) 2009-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Routines, for the position array of degrees of freedom:
!   local position in the node (pos_array_local)
!   position in the system vector (pos_array)
!   position in the system vector regarding constraints (pos_array_constraint)
!   position in the system vector regarding dependencies (pos_array_dependency)
!   position in a vector (pos_array_vec).

module pos_array_m

  use kind_defs_m
  use mesh_m
  use problem_defs_m

  implicit none


contains


! Local position array for nodes

  subroutine pos_array_local_node ( problem, nodenr, dof, pos, physqarr, layer )

    type(problem_t), intent(in)  :: problem
    integer, intent(in) :: nodenr

!   dof is the number of degrees of freedom involved
    integer, intent(out) :: dof

!   the position array
    integer, dimension(:), intent(out) :: pos

!   if physarr is present only physical quantities in this array are included.
!   For example:
!     physarr=(/2,1/)
!   means that pos will contain the positions of the physical quantities 2 and
!   1, in that order.
    integer, dimension(:), intent(in), optional :: physqarr

!   if layer is present and layer > 0 the element degrees of freedom are
!   restricted to the specified layer.
!   layer=0 is identical to layer not present
    integer, intent(in), optional :: layer


!   this routine computes the local positions of degrees of freedom of a node

    integer :: llayer
    integer :: physq, bp, nndof, k


    llayer = set_optional ( variable=layer, default=0 )

    if ( present(physqarr) ) then

!     physical quantities specified (PD or PLD)

      dof = 0

!     positions

      do physq = 1, size(physqarr)

        bp = sum( problem%vec_nodnumdegfd(nodenr+1,&
                                  &problem%physq(1:physqarr(physq)-1)) &
                   - problem%vec_nodnumdegfd(nodenr,&
                                  &problem%physq(1:physqarr(physq)-1)) )
        nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                  &problem%physq(physqarr(physq))) &
                   - problem%vec_nodnumdegfd(nodenr,&
                                  &problem%physq(physqarr(physq)))

        if ( llayer > 0 ) call adjust_for_layer

!       local positions in node
        pos(dof+1:dof+nndof) = [(k,k=bp+1,bp+nndof)]

        dof = dof + nndof

      end do

    else if ( problem%nphysq > 0 ) then

!     all physical quantities (PD or PLD)

      dof = 0

!     positions

      do physq = 1, problem%nphysq

        bp = sum( problem%vec_nodnumdegfd(nodenr+1,problem%physq(1:physq-1)) &
                   - problem%vec_nodnumdegfd(nodenr,problem%physq(1:physq-1)) )
        nndof = problem%vec_nodnumdegfd(nodenr+1,problem%physq(physq)) &
                   - problem%vec_nodnumdegfd(nodenr,problem%physq(physq))

        if ( llayer > 0 ) call adjust_for_layer

!       local positions in node
        pos(dof+1:dof+nndof) = [(k,k=bp+1,bp+nndof)]

        dof = dof + nndof

      end do

    else

!     no physical quantities (D or LD)

      bp    = 0
      nndof = problem%nodnumdegfd(nodenr+1) - problem%nodnumdegfd(nodenr)

      if ( llayer > 0 ) call adjust_for_layer

      dof = nndof

!     positions in sysvector and sysmatrix (renumbered)
      pos(:dof) = [(k,k=bp+1,bp+nndof)]

    end if


  contains


!   modify start pointer and number of degrees for layer

    subroutine adjust_for_layer

      integer, dimension(1) :: numl, pnuml

      call find_layer_in_nodes ( problem, llayer, nodes=[nodenr], &
          numl=numl, pnuml=pnuml )

      if ( numl(1) > 1 ) then
        nndof = nndof / numl(1)
        bp = bp + nndof * pnuml(1)
      else if ( numl(1) == 0 ) then
        nndof = 0
      end if

    end subroutine adjust_for_layer

  end subroutine pos_array_local_node


! Position array in sysvector

  subroutine pos_array ( mesh, problem, elgrp, elem, dof, pos, physqarr, &
    order, layer, lp, pq )

    type(mesh_t), intent(in)  :: mesh
    type(problem_t), intent(in)  :: problem
    integer, intent(in) :: elgrp, elem

!   dof is the number of degrees of freedom involved
    integer, intent(out) :: dof

!   the position array
    integer, dimension(:), intent(out) :: pos

!   if physarr is present only physical quantities in this array are included.
!   For example:
!     physarr=(/2,1/)
!   means that pos will contain the positions of the physical quantities 2 and
!   1, in that order.
    integer, dimension(:), intent(in), optional :: physqarr

!   order gives the sequence order of the degrees of freedom in the nodes of
!   the element. There are only two possibilities and affect only the two loops
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
!   The default of order is:
!     order = 'ND' : if no physical quantities have been defined (nphysq=0)
!     order = 'DN' : if physical quantities have been defined (nphysq>0)
    character(len=*), intent(in), optional :: order

!   if layer is present (in heading) and layer > 0 the element degrees
!   of freedom are restricted to the specified layer.
!   layer=0 is identical to layer not present
    integer, intent(in), optional :: layer

!   if lp is present and layer is present (in heading) and layer > 0, then
!   lp(i) = .true. if layer is present in i^th node else .false.
    logical, dimension(:), intent(out), optional :: lp

!   if present, of all the degrees denoted in pos the physical quantity number
!   is given. Only filled if physical quantities have been defined
    integer, dimension(:), intent(out), optional :: pq

!   this routine computes the positions of the degrees of freedom of an element.

    integer :: llayer
    integer :: nodenr, node, physq, bp, nndof, deg, maxdeg
    integer, dimension(mesh%elnumnod(elgrp)) :: nodes, numl, pnuml, maxdegfac
    character(len=2) :: lorder

    llayer = set_optional ( variable=layer, default=0 )

!   find layer in the nodes of the element

    if ( llayer > 0 ) then

      call find_layer_in_nodes ( problem, llayer, &
        nodes=mesh%topology(elgrp)%a(:,elem), numl=numl, pnuml=pnuml, lp=lp )

      maxdegfac = max ( numl, 1 )

    else

      maxdegfac = 1

    end if


!   find positions

    if ( present(physqarr) ) then

!     only physical quantities specified

      lorder = set_optional ( variable=order, default='DN' )

      if ( lorder == 'ND' ) then

!       ordering of degrees: physical quantity outside loop,
!       then nodal points and degrees inside the loop (PND of PNLD).

        dof = 0

        do physq = 1, size(physqarr)

          do node = 1, mesh%elnumnod(elgrp)

            nodenr = mesh%topology(elgrp)%a(node,elem)

!           positions

            bp = problem%nodnumdegfd(nodenr) + &
                  sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physqarr(physq)-1)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physqarr(physq)-1)) )
            nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physqarr(physq))) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physqarr(physq)))

            if ( llayer > 0 ) call adjust_for_layer

!           positions in sysvector and sysmatrix (renumbered)
            pos(dof+1:dof+nndof) = problem%degfdperm(bp+1:bp+nndof,2)

            if ( present(pq) ) pq(dof+1:dof+nndof) = physqarr(physq)

            dof = dof + nndof

          end do

        end do

      else if ( lorder == 'DN' ) then

!       ordering of degrees: physical quantity outside loop,
!       then degrees and nodes inside the loop (PDN or PLDN).

        dof = 0

        do physq = 1, size(physqarr)

          nodes = mesh%topology(elgrp)%a(:,elem)
          maxdeg = maxval( ( problem%vec_nodnumdegfd(nodes+1,&
                                      &problem%physq(physqarr(physq))) &
                           - problem%vec_nodnumdegfd(nodes,&
                                      &problem%physq(physqarr(physq))) &
                            ) / maxdegfac )

          do deg = 1, maxdeg

            do node = 1, mesh%elnumnod(elgrp)

              nodenr = mesh%topology(elgrp)%a(node,elem)

!             positions

              bp = problem%nodnumdegfd(nodenr) + &
                  sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physqarr(physq)-1)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physqarr(physq)-1)) )
                nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                    &problem%physq(physqarr(physq))) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physqarr(physq)))

              if ( llayer > 0 ) call adjust_for_layer

!             positions in sysvector and sysmatrix (renumbered)
              if ( deg <= nndof ) then
                pos(dof+1) = problem%degfdperm(bp+deg,2)
                if ( present(pq) ) pq(dof+1) = physqarr(physq)
                dof = dof + 1
              end if

            end do

          end do

        end do

      end if

    else if ( problem%nphysq > 0 ) then

!     all physical quantities

      lorder = set_optional ( variable=order, default='DN' )

      if ( lorder == 'ND' ) then

!       ordering of degrees: physical quantity outside loop,
!       then nodal points and degrees inside the loop (PND of PNLD).

        dof = 0

        do physq = 1, problem%nphysq

          do node = 1, mesh%elnumnod(elgrp)

            nodenr = mesh%topology(elgrp)%a(node,elem)

!           positions

            bp = problem%nodnumdegfd(nodenr) + &
                  sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physq-1)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physq-1)) )
            nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physq)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physq))

            if ( llayer > 0 ) call adjust_for_layer

!           positions in sysvector and sysmatrix (renumbered)
            pos(dof+1:dof+nndof) = problem%degfdperm(bp+1:bp+nndof,2)

            if ( present(pq) ) pq(dof+1:dof+nndof) = physq

            dof = dof + nndof

          end do

        end do

      else if ( lorder == 'DN' ) then

!       ordering of degrees: physical quantity outside loop,
!       then degrees and nodes inside the loop (PDN or PLDN).

        dof = 0

        do physq = 1, problem%nphysq

          nodes = mesh%topology(elgrp)%a(:,elem)
          maxdeg = maxval( ( problem%vec_nodnumdegfd(nodes+1,&
                                      &problem%physq(physq)) &
                           - problem%vec_nodnumdegfd(nodes,  &
                                      &problem%physq(physq)) &
                            ) / maxdegfac )

          do deg = 1, maxdeg

            do node = 1, mesh%elnumnod(elgrp)

              nodenr = mesh%topology(elgrp)%a(node,elem)

!             positions

              bp = problem%nodnumdegfd(nodenr) + &
                  sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physq-1)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physq-1)) )
                nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                    &problem%physq(physq)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physq))

              if ( llayer > 0 ) call adjust_for_layer

!             positions in sysvector and sysmatrix (renumbered)
              if ( deg <= nndof ) then
                pos(dof+1) = problem%degfdperm(bp+deg,2)
                if ( present(pq) ) pq(dof+1) = physq
                dof = dof + 1
              end if

            end do

          end do

        end do

      end if

    else

!     no physical quantities defined

      lorder = set_optional ( variable=order, default='ND' )

      if ( lorder == 'ND' ) then

!       nodal points outside and degrees inside the loop (ND of NLD).

        dof = 0

        do node = 1, mesh%elnumnod(elgrp)

          nodenr = mesh%topology(elgrp)%a(node,elem)
          bp     = problem%nodnumdegfd(nodenr)
          nndof  = problem%nodnumdegfd(nodenr+1) - problem%nodnumdegfd(nodenr)

          if ( llayer > 0 ) call adjust_for_layer

!         positions in sysvector and sysmatrix (renumbered)
          pos(dof+1:dof+nndof) = problem%degfdperm(bp+1:bp+nndof,2)

          dof = dof + nndof

        end do

      else if ( lorder == 'DN' ) then

!       degrees outside and nodes inside the loop (DN or LDN).

        dof = 0

        nodes = mesh%topology(elgrp)%a(:,elem)
        maxdeg = maxval( ( problem%nodnumdegfd(nodes+1) &
                           - problem%nodnumdegfd(nodes) &
                         ) / maxdegfac )

        do deg = 1, maxdeg

          do node = 1, mesh%elnumnod(elgrp)

            nodenr = mesh%topology(elgrp)%a(node,elem)

!           positions

            bp = problem%nodnumdegfd(nodenr)
            nndof  = problem%nodnumdegfd(nodenr+1) - problem%nodnumdegfd(nodenr)

            if ( llayer > 0 ) call adjust_for_layer

!           positions in sysvector and sysmatrix (renumbered)
            if ( deg <= nndof ) then
              pos(dof+1) = problem%degfdperm(bp+deg,2)
              dof = dof + 1
            end if

          end do

        end do

      end if

    end if

  contains

    subroutine adjust_for_layer

      if ( numl(node) > 1 ) then
        nndof = nndof / numl(node)
        bp = bp + nndof * pnuml(node)
      else if ( numl(node) == 0 ) then
        nndof = 0
      end if

    end subroutine adjust_for_layer

  end subroutine pos_array


! Position array in sysvector for nodes

  subroutine pos_array_node ( problem, nodenr, dof, pos, physqarr, layer )

    type(problem_t), intent(in)  :: problem
    integer, intent(in) :: nodenr

!   dof is the number of degrees of freedom involved
    integer, intent(out) :: dof

!   the position array
    integer, dimension(:), intent(out) :: pos

!   if physarr is present only physical quantities in this array are included.
!   For example:
!     physarr=(/2,1/)
!   means that pos will contain the positions of the physical quantities 2 and
!   1, in that order.
    integer, dimension(:), intent(in), optional :: physqarr

!   if layer is present and layer > 0 the element degrees of freedom are
!   restricted to the specified layer.
!   layer=0 is identical to layer not present
    integer, intent(in), optional :: layer


!   this routine computes the positions of the degrees of freedom of a node

    integer :: llayer
    integer :: physq, bp, nndof


    llayer = set_optional ( variable=layer, default=0 )

    if ( present(physqarr) ) then

!     physical quantities specified

      dof = 0

!     positions

      do physq = 1, size(physqarr)

        bp = problem%nodnumdegfd(nodenr) + &
              sum( problem%vec_nodnumdegfd(nodenr+1,&
                                  &problem%physq(1:physqarr(physq)-1)) &
                   - problem%vec_nodnumdegfd(nodenr,&
                                  &problem%physq(1:physqarr(physq)-1)) )
        nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                  &problem%physq(physqarr(physq))) &
                   - problem%vec_nodnumdegfd(nodenr,&
                                  &problem%physq(physqarr(physq)))

        if ( llayer > 0 ) call adjust_for_layer

!       positions in sysvector and sysmatrix (renumbered)
        pos(dof+1:dof+nndof) = problem%degfdperm(bp+1:bp+nndof,2)

        dof = dof + nndof

      end do

    else if ( problem%nphysq > 0 ) then

!     all physical quantities

      dof = 0

!     positions

      do physq = 1, problem%nphysq

        bp = problem%nodnumdegfd(nodenr) + &
              sum( problem%vec_nodnumdegfd(nodenr+1,&
                                  &problem%physq(1:physq-1)) &
                   - problem%vec_nodnumdegfd(nodenr,&
                                  &problem%physq(1:physq-1)) )
        nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                  &problem%physq(physq)) &
                   - problem%vec_nodnumdegfd(nodenr,&
                                  &problem%physq(physq))

        if ( llayer > 0 ) call adjust_for_layer

!       positions in sysvector and sysmatrix (renumbered)
        pos(dof+1:dof+nndof) = problem%degfdperm(bp+1:bp+nndof,2)

        dof = dof + nndof

      end do

    else

!     no physical quantities

      bp    = problem%nodnumdegfd(nodenr)
      nndof = problem%nodnumdegfd(nodenr+1) - problem%nodnumdegfd(nodenr)

      if ( llayer > 0 ) call adjust_for_layer

      dof = nndof

!     positions in sysvector and sysmatrix (renumbered)
      pos(:dof) = problem%degfdperm(bp+1:bp+dof,2)

    end if


  contains


!   modify start pointer and number of degrees for layer

    subroutine adjust_for_layer

      integer, dimension(1) :: numl, pnuml

      call find_layer_in_nodes ( problem, llayer, nodes=[nodenr], &
          numl=numl, pnuml=pnuml )

      if ( numl(1) > 1 ) then
        nndof = nndof / numl(1)
        bp = bp + nndof * pnuml(1)
      else if ( numl(1) == 0 ) then
        nndof = 0
      end if

    end subroutine adjust_for_layer

  end subroutine pos_array_node


! Position array in sysvector for geometries (element)

  subroutine pos_array_geometry ( problem, geometry, elem, dof, pos, &
    physqarr, order, layer, lp )

    type(problem_t), intent(in)  :: problem
    type(geometry_t), intent(in) :: geometry
    integer, intent(in) :: elem

!   dof is the number of degrees of freedom involved
    integer, intent(out) :: dof

!   the position array
    integer, dimension(:), intent(out) :: pos

!   if physarr is present only physical quantities in this array are included.
!   For example:
!     physarr=(/2,1/)
!   means that pos will contain the positions of the physical quantities 2 and
!   1, in that order.
    integer, dimension(:), intent(in), optional :: physqarr

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

!   if layer is present (in heading) and layer > 0 the element degrees
!   restricted to the specified layer.
!   layer=0 is identical to layer not present
    integer, intent(in), optional :: layer

!   if lp is present and layer is present (in heading) and layer > 0, then
!   lp(i) = .true. if layer is present in i^th node else .false.
    logical, dimension(:), intent(out), optional :: lp


!   this routine computes the positions of the degrees of freedom of an element
!   on a geometry

    integer :: llayer
    integer :: nodenr, node, physq, bp, nndof, deg, maxdeg
    integer, dimension(geometry%elnumnod) :: nodes, numl, pnuml, maxdegfac
    character(len=2) :: lorder


    llayer = set_optional ( variable=layer, default=0 )

!   find layer in the nodes of the element

    if ( llayer > 0  ) then

      call find_layer_in_nodes ( problem, llayer, &
        nodes=geometry%topology(:,elem,2), numl=numl, pnuml=pnuml, lp=lp )

      maxdegfac = max ( numl, 1 )

    else

      maxdegfac = 1

    end if


!   find positions

    if ( present(physqarr) ) then

!     only physical quantities specified

      lorder = set_optional ( variable=order, default='DN' )

      if ( lorder == 'ND' ) then

!       ordering of degrees: physical quantity outside loop,
!       then nodal points and degrees inside the loop (PND of PNLD).

        dof = 0

        do physq = 1, size(physqarr)

          do node = 1, geometry%elnumnod

            nodenr = geometry%topology(node,elem,2)

!           positions

            bp = problem%nodnumdegfd(nodenr) + &
                  sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physqarr(physq)-1)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physqarr(physq)-1)) )
            nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physqarr(physq))) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physqarr(physq)))

            if ( llayer > 0 ) call adjust_for_layer

!           positions in sysvector and sysmatrix (renumbered)
            pos(dof+1:dof+nndof) = problem%degfdperm(bp+1:bp+nndof,2)

            dof = dof + nndof

          end do

        end do

      else if ( lorder == 'DN' ) then

!       ordering of degrees: physical quantity outside loop,
!       then degrees and nodes inside the loop (PDN or PLDN).

        dof = 0

        do physq = 1, size(physqarr)

          nodes = geometry%topology(:,elem,2)
          maxdeg = maxval( ( problem%vec_nodnumdegfd(nodes+1,&
                                      &problem%physq(physqarr(physq))) &
                           - problem%vec_nodnumdegfd(nodes,&
                                      &problem%physq(physqarr(physq))) &
                           ) / maxdegfac )

          do deg = 1, maxdeg

            do node = 1, geometry%elnumnod

              nodenr = geometry%topology(node,elem,2)

!             positions

              bp = problem%nodnumdegfd(nodenr) + &
                  sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physqarr(physq)-1)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physqarr(physq)-1)) )
              nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physqarr(physq))) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physqarr(physq)))

              if ( llayer > 0 ) call adjust_for_layer

!             positions in sysvector and sysmatrix (renumbered)
              if ( deg <= nndof ) then
                pos(dof+1) = problem%degfdperm(bp+deg,2)
                dof = dof + 1
              end if

            end do

          end do

        end do

      end if

    else if ( problem%nphysq > 0 ) then

!     all physical quantities

      lorder = set_optional ( variable=order, default='DN' )

      if ( lorder == 'ND' ) then

!       ordering of degrees: physical quantity outside loop,
!       then nodal points and degrees inside the loop (PND of PNLD).

        dof = 0

        do physq = 1,  problem%nphysq

          do node = 1, geometry%elnumnod

            nodenr = geometry%topology(node,elem,2)

!           positions

            bp = problem%nodnumdegfd(nodenr) + &
                  sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physq-1)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physq-1)) )
            nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physq)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physq))

            if ( llayer > 0 ) call adjust_for_layer

!           positions in sysvector and sysmatrix (renumbered)
            pos(dof+1:dof+nndof) = problem%degfdperm(bp+1:bp+nndof,2)

            dof = dof + nndof

          end do

        end do

      else if ( lorder == 'DN' ) then

!       ordering of degrees: physical quantity outside loop,
!       then degrees and nodes inside the loop (PDN or PLDN).

        dof = 0

        do physq = 1, problem%nphysq

          nodes = geometry%topology(:,elem,2)
          maxdeg = maxval( ( problem%vec_nodnumdegfd(nodes+1,&
                                      &problem%physq(physq)) &
                           - problem%vec_nodnumdegfd(nodes,&
                                      &problem%physq(physq)) &
                           ) / maxdegfac )

          do deg = 1, maxdeg

            do node = 1, geometry%elnumnod

              nodenr = geometry%topology(node,elem,2)

!             positions

              bp = problem%nodnumdegfd(nodenr) + &
                  sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physq-1)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physq-1)) )
              nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physq)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physq))

              if ( llayer > 0 ) call adjust_for_layer

!             positions in sysvector and sysmatrix (renumbered)
              if ( deg <= nndof ) then
                pos(dof+1) = problem%degfdperm(bp+deg,2)
                dof = dof + 1
              end if

            end do

          end do

        end do

      end if

    else

!     no physical quantities defined

      lorder = set_optional ( variable=order, default='ND' )

      if ( lorder == 'ND' ) then

!       nodal points outside and degrees inside the loop (ND of NLD).

        dof = 0

        do node = 1, geometry%elnumnod

          nodenr = geometry%topology(node,elem,2)
          bp     = problem%nodnumdegfd(nodenr)
          nndof  = problem%nodnumdegfd(nodenr+1) - problem%nodnumdegfd(nodenr)

          if ( llayer > 0 ) call adjust_for_layer

!         positions in sysvector and sysmatrix (renumbered)
          pos(dof+1:dof+nndof) = problem%degfdperm(bp+1:bp+nndof,2)

          dof = dof + nndof

        end do

      else  if ( lorder == 'DN' ) then

!       degrees outside and nodes inside the loop (DN or LDN).

        dof = 0

        nodes = geometry%topology(:,elem,2)
        maxdeg = maxval( ( problem%nodnumdegfd(nodes+1) &
                           - problem%nodnumdegfd(nodes) &
                         ) / maxdegfac )

        do deg = 1, maxdeg

          do node = 1, geometry%elnumnod

            nodenr = geometry%topology(node,elem,2)

!           positions

            bp = problem%nodnumdegfd(nodenr)
            nndof  = problem%nodnumdegfd(nodenr+1) - problem%nodnumdegfd(nodenr)

            if ( llayer > 0 ) call adjust_for_layer

!           positions in sysvector and sysmatrix (renumbered)
            if ( deg <= nndof ) then
              pos(dof+1) = problem%degfdperm(bp+deg,2)
              dof = dof + 1
            end if

          end do

        end do

      end if

    end if

  contains

    subroutine adjust_for_layer

      if ( numl(node) > 1 ) then
        nndof = nndof / numl(node)
        bp = bp + nndof * pnuml(node)
      else if ( numl(node) == 0 ) then
        nndof = 0
      end if

    end subroutine adjust_for_layer

  end subroutine pos_array_geometry


! Position array in sysvector for geometries (full geometry)

  subroutine pos_array_fullgeometry ( problem, geometry, dof, pos, &
    physqarr, order, layer, lp )

    type(problem_t), intent(in)  :: problem
    type(geometry_t), intent(in) :: geometry

!   dof is the number of degrees of freedom involved
    integer, intent(out) :: dof

!   the position array
    integer, dimension(:), intent(out) :: pos

!   if physarr is present only physical quantities in this array are included.
!   For example:
!     physarr=(/2,1/)
!   means that pos will contain the positions of the physical quantities 2 and
!   1, in that order.
    integer, dimension(:), intent(in), optional :: physqarr

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

!   if layer is present (in heading) and layer > 0 the element degrees
!   restricted to the specified layer.
!   layer=0 is identical to layer not present
    integer, intent(in), optional :: layer

!   if lp is present and layer is present (in heading) and layer > 0, then
!   lp(i) = .true. if layer is present in i^th node else .false.
    logical, dimension(:), intent(out), optional :: lp


!   this routine computes the positions of the degrees of freedom of an element
!   on a geometry

    integer :: llayer
    integer :: nodenr, node, physq, bp, nndof, deg, maxdeg
    integer, allocatable, dimension(:) :: nodes, numl, pnuml, maxdegfac
    character(len=2) :: lorder

    allocate ( nodes(geometry%nnodes), numl(geometry%nnodes), &
               pnuml(geometry%nnodes), maxdegfac(geometry%nnodes) )

    llayer = set_optional ( variable=layer, default=0 )

!   find layer in the nodes of the element

    if ( llayer > 0 ) then

      call find_layer_in_nodes ( problem, llayer, &
        nodes=geometry%nodes, numl=numl, pnuml=pnuml, lp=lp )

      maxdegfac = max ( numl, 1 )

    else

      maxdegfac = 1

    end if

!   find positions

    if ( present(physqarr) ) then

!     only physical quantities specified

      lorder = set_optional ( variable=order, default='DN' )

      if ( lorder == 'ND' ) then

!       ordering of degrees: physical quantity outside loop,
!       then nodal points and degrees inside the loop (PND of PNLD).

        dof = 0

        do physq = 1, size(physqarr)

          do node = 1, geometry%nnodes

            nodenr = geometry%nodes(node)

!           positions

            bp = problem%nodnumdegfd(nodenr) + &
                  sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physqarr(physq)-1)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physqarr(physq)-1)) )
            nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physqarr(physq))) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physqarr(physq)))

            if ( llayer > 0 ) call adjust_for_layer

!           positions in sysvector and sysmatrix (renumbered)
            pos(dof+1:dof+nndof) = problem%degfdperm(bp+1:bp+nndof,2)

            dof = dof + nndof

          end do

        end do

      else if ( lorder == 'DN' ) then

!       ordering of degrees: physical quantity outside loop,
!       then degrees and nodes inside the loop (PDN or PLDN).

        dof = 0

        do physq = 1, size(physqarr)

          nodes = geometry%nodes
          maxdeg = maxval( ( problem%vec_nodnumdegfd(nodes+1,&
                                      &problem%physq(physqarr(physq))) &
                           - problem%vec_nodnumdegfd(nodes,&
                                      &problem%physq(physqarr(physq))) &
                           ) / maxdegfac )

          do deg = 1, maxdeg

            do node = 1, geometry%nnodes

              nodenr = geometry%nodes(node)

!             positions

              bp = problem%nodnumdegfd(nodenr) + &
                  sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physqarr(physq)-1)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physqarr(physq)-1)) )
              nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physqarr(physq))) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physqarr(physq)))

              if ( llayer > 0 ) call adjust_for_layer

!             positions in sysvector and sysmatrix (renumbered)
              if ( deg <= nndof ) then
                pos(dof+1) = problem%degfdperm(bp+deg,2)
                dof = dof + 1
              end if

            end do

          end do

        end do

      end if

    else if ( problem%nphysq > 0 ) then

!     all physical quantities

      lorder = set_optional ( variable=order, default='DN' )

      if ( lorder == 'ND' ) then

!       ordering of degrees: physical quantity outside loop,
!       then nodal points and degrees inside the loop (PND of PNLD).

        dof = 0

        do physq = 1, problem%nphysq

          do node = 1, geometry%nnodes

            nodenr = geometry%nodes(node)

!           positions

            bp = problem%nodnumdegfd(nodenr) + &
                  sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physq-1)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physq-1)) )
            nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physq)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physq))

            if ( llayer > 0 ) call adjust_for_layer

!           positions in sysvector and sysmatrix (renumbered)
            pos(dof+1:dof+nndof) = problem%degfdperm(bp+1:bp+nndof,2)

            dof = dof + nndof

          end do

        end do

      else if ( lorder == 'DN' ) then

!       ordering of degrees: physical quantity outside loop,
!       then degrees and nodes inside the loop (PDN or PLDN).

        dof = 0

        do physq = 1, problem%nphysq

          nodes = geometry%nodes
          maxdeg = maxval( ( problem%vec_nodnumdegfd(nodes+1,&
                                      &problem%physq(physq)) &
                           - problem%vec_nodnumdegfd(nodes,&
                                      &problem%physq(physq)) &
                           ) / maxdegfac )

          do deg = 1, maxdeg

            do node = 1, geometry%nnodes

              nodenr = geometry%nodes(node)

!             positions

              bp = problem%nodnumdegfd(nodenr) + &
                  sum( problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(1:physq-1)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(1:physq-1)) )
              nndof = problem%vec_nodnumdegfd(nodenr+1,&
                                      &problem%physq(physq)) &
                       - problem%vec_nodnumdegfd(nodenr,&
                                      &problem%physq(physq))

              if ( llayer > 0 ) call adjust_for_layer

!             positions in sysvector and sysmatrix (renumbered)
              if ( deg <= nndof ) then
                pos(dof+1) = problem%degfdperm(bp+deg,2)
                dof = dof + 1
              end if

            end do

          end do

        end do

      end if

    else

!     no physical quantities defined

      lorder = set_optional ( variable=order, default='ND' )

      if ( lorder == 'ND' ) then

        dof = 0

        do node = 1, geometry%nnodes

          nodenr = geometry%nodes(node)
          bp     = problem%nodnumdegfd(nodenr)
          nndof  = problem%nodnumdegfd(nodenr+1) - problem%nodnumdegfd(nodenr)

          if ( llayer > 0 ) call adjust_for_layer

!         positions in sysvector and sysmatrix (renumbered)
          pos(dof+1:dof+nndof) = problem%degfdperm(bp+1:bp+nndof,2)

          dof = dof + nndof

        end do

      else if ( lorder == 'DN' ) then

!       degrees outside and nodes inside the loop (DN or LDN).

        dof = 0

        nodes = geometry%nodes
        maxdeg = maxval( ( problem%nodnumdegfd(nodes+1) &
                           - problem%nodnumdegfd(nodes) &
                         ) / maxdegfac )

        do deg = 1, maxdeg

          do node = 1, geometry%nnodes

            nodenr = geometry%nodes(node)

!           positions

            bp = problem%nodnumdegfd(nodenr)
            nndof  = problem%nodnumdegfd(nodenr+1) - problem%nodnumdegfd(nodenr)

            if ( llayer > 0 ) call adjust_for_layer

!           positions in sysvector and sysmatrix (renumbered)
            if ( deg <= nndof ) then
              pos(dof+1) = problem%degfdperm(bp+deg,2)
              dof = dof + 1
            end if

          end do

        end do

      end if

    end if

    deallocate ( nodes, numl, pnuml, maxdegfac )

  contains

    subroutine adjust_for_layer

      if ( numl(node) > 1 ) then
        nndof = nndof / numl(node)
        bp = bp + nndof * pnuml(node)
      else if ( numl(node) == 0 ) then
        nndof = 0
      end if

    end subroutine adjust_for_layer

  end subroutine pos_array_fullgeometry


! Position array in sysvector for constraints

  subroutine pos_array_constraint ( problem, constr, dof, pos, geometry, &
    object, elem, node, globalc, addunknowns, order )

    type(problem_t), intent(in)  :: problem

!   constr is the constraint number
    integer, intent(in) :: constr

!   dof is the number of degrees of freedom involved
    integer, intent(out) :: dof

!   the position array
    integer, dimension(:), intent(out) :: pos

!   NOTE: only one of geometry/object can be present:

!   the geometry when elem is present
    type(geometry_t), intent(in), optional :: geometry

!   the object when elem is present
    type(object_t), intent(in), optional :: object

!   these parameters specify which positions are needed
!     elem: the element number
!     node: the node number
!     globalc: global constraints
!     addunknowns: additional unknowns
    integer, intent(in), optional :: elem, node
    logical, intent(in), optional :: globalc, addunknowns

!   this routine computes the positions of the degrees of freedom of a
!   constraint

!   ordering. For the two loops over nodal points and degrees of freedom we have
!   two possibilities:
!     order = 'ND' : the most inner loop is over the degrees of freedom
!     order = 'DN' : the most inner loop is over the nodal points
!   this only makes sense for the elem parameter.
!   If not present the value order='DN' will be used (default)
    character(len=*), intent(in), optional :: order

    integer :: nodenr, bp, nndof, nnodes, nod, deg, maxdeg
    character(len=2) :: lorder


    if ( present(elem) ) then

!     constraints in element elem

      lorder = set_optional ( variable=order, default='DN' )

      if ( present(geometry) ) then
        nnodes = geometry%elnumnod
      else if ( present(object) ) then
        nnodes = object%elnumnod
      end if

      if ( lorder == 'ND' ) then

!       ordering of degrees: degrees in inside loop

        dof = 0

        do nod = 1, nnodes

          if ( present(geometry) ) then
            nodenr = geometry%topology(nod,elem,1)
          else if ( present(object) ) then
            nodenr = object%topology(nod,elem)
          end if

!         positions

          bp = problem%constraints(constr)%nodnumdegfd(nodenr)
          nndof = problem%constraints(constr)%nodnumdegfd(nodenr+1) &
                     - problem%constraints(constr)%nodnumdegfd(nodenr)

!         positions in sysvector and sysmatrix (renumbered)
          pos(dof+1:dof+nndof) = problem%degfdperm(bp+1:bp+nndof,2)

          dof = dof + nndof

        end do

      else if ( lorder == 'DN' ) then

!       ordering of degrees: nodes in inside loop

        dof = 0

        maxdeg = maxval( problem%constraints(constr)%elnumdegfd )

        do deg = 1, maxdeg

          do nod = 1, nnodes

            if ( present(geometry) ) then
              nodenr = geometry%topology(nod,elem,1)
            else if ( present(object) ) then
              nodenr = object%topology(nod,elem)
            end if

!           positions

            bp = problem%constraints(constr)%nodnumdegfd(nodenr)
            nndof = problem%constraints(constr)%nodnumdegfd(nodenr+1) &
                       - problem%constraints(constr)%nodnumdegfd(nodenr)

!           positions in sysvector and sysmatrix (renumbered)
            if ( deg <= nndof ) then
              pos(dof+1) = problem%degfdperm(bp+deg,2)
              dof = dof + 1
            end if

          end do

        end do

      end if

    else if ( present(node) ) then

!     constraints in node

      bp = problem%constraints(constr)%nodnumdegfd(node)
      nndof = problem%constraints(constr)%nodnumdegfd(node+1) &
                 - problem%constraints(constr)%nodnumdegfd(node)

!     positions in sysvector and sysmatrix (renumbered)
      pos(1:nndof) = problem%degfdperm(bp+1:bp+nndof,2)

      dof = nndof

    else if ( present(globalc) ) then

!     global constraints

      bp = problem%constraints(constr)%globnumdegfd(1)
      nndof = problem%constraints(constr)%globnumdegfd(2) &
                 - problem%constraints(constr)%globnumdegfd(1)

!     positions in sysvector and sysmatrix (renumbered)
      pos(1:nndof) = problem%degfdperm(bp+1:bp+nndof,2)

      dof = nndof

    else if ( present(addunknowns) ) then

!     additional unknowns

      bp = problem%constraints(constr)%addnumdegfd(1)
      nndof = problem%constraints(constr)%addnumdegfd(2) &
                 - problem%constraints(constr)%addnumdegfd(1)

!     positions in sysvector and sysmatrix (renumbered)
      pos(1:nndof) = problem%degfdperm(bp+1:bp+nndof,2)

      dof = nndof

    else

      write(*,'(/2a/)') &
        'Error in pos_array_constraint: ', &
        ' either elem, node, globalc or addunknowns must be present'
      stop

    end if

  end subroutine pos_array_constraint


! Position array in sysvector for dependencies (additional unknowns only)

  subroutine pos_array_dependency ( problem, dep, dof, pos )

    type(problem_t), intent(in)  :: problem

!   constr is the dependency number
    integer, intent(in) :: dep

!   dof is the number of degrees of freedom involved
    integer, intent(out) :: dof

!   the position array
    integer, dimension(:), intent(out) :: pos

!   this routine computes the positions of the degrees of freedom of a
!   dependency (additonal unknowns only)

    integer :: bp, nndof

!   additional unknowns

    bp = problem%dependencies(dep)%addnumdegfd(1)
    nndof = problem%dependencies(dep)%addnumdegfd(2) &
                 - problem%dependencies(dep)%addnumdegfd(1)

!   positions in sysvector and sysmatrix (renumbered)
    pos(1:nndof) = problem%degfdperm(bp+1:bp+nndof,2)

    dof = nndof

  end subroutine pos_array_dependency


! Position array in vector

  subroutine pos_array_vec ( mesh, problem, elgrp, elem, dof, pos, vec, &
    elementwise, order, layer, lp )

    type(mesh_t), intent(in)  :: mesh
    type(problem_t), intent(in)  :: problem
    integer, intent(in) :: elgrp, elem

!   dof is the number of degrees of freedom involved
    integer, intent(out) :: dof

!   the position array
    integer, dimension(:), intent(out) :: pos

!   the vector number
    integer, intent(in) :: vec

!   indicates whether data are stored elementwise (.true.) or not (.false.).
    logical, intent(in) :: elementwise

!   The parameter order determines the sequence of the degrees of freedom on
!   elementlevel. There are only two possibilities and affect only the two loops
!   over nodal points and degrees of freedom
!     order = 'ND' : the most inner loop is over the degrees of freedom
!     order = 'DN' : the most inner loop is over the nodal points
!   If problem%numlayers > 0 the degrees are stacked into layers and we have
!         order = 'ND' :      NLD
!         order = 'DN' :      LDN
!   See also the userguide for further explanation.
!   The default of order is order = 'DN'
    character(len=*), intent(in), optional :: order

!   if layer is present and layer > 0 the element degrees of freedom are
!   restricted to the specified layer.
!   layer=0 is identical to layer not present
!   NOTE: if elementwise=.true., layer has no effect.
    integer, intent(in), optional :: layer

!   if lp is present and layer is present (in heading) and layer > 0, then
!   lp(i) = .true. if layer is present in i^th node else .false.
    logical, dimension(:), intent(out), optional :: lp


!   this routine computes the positions of the components of a vector inside
!   an element.


    integer :: nodenr, node, sp, nndof, deg, maxdeg, i, spe, elg, llayer
    integer, dimension(elgrp) :: w
    integer, dimension(mesh%elnumnod(elgrp)) :: nodes, numl, pnuml, maxdegfac
    character(len=2) :: lorder

    llayer = set_optional ( variable=layer, default=0 )
    lorder = set_optional ( variable=order, default='DN' )

    if ( elementwise ) then

      if ( llayer > 0 ) then
        write(*,'(/a/2a/a/)') &
          'Error in pos_array_vec: ', &
          '  layer not available for elementwise=.true.'
        stop
      end if

!     compute starting position of element data

      do elg = 1, elgrp - 1
        w(elg) = mesh%grpnumel(elg) * &
                   sum( problem%vec_elnumdegfd(elg)%a(:,vec) )
      end do
      w(elgrp) = ( elem - 1) * &
                   sum( problem%vec_elnumdegfd(elgrp)%a(:,vec) )
      spe = sum(w)

    end if

!   find layer in the nodes of the element

    if ( llayer > 0 .and. .not. elementwise ) then

      call find_layer_in_nodes ( problem, llayer, &
        nodes=mesh%topology(elgrp)%a(:,elem), numl=numl, pnuml=pnuml, lp=lp )

      maxdegfac = max ( numl, 1 )

    else

      maxdegfac = 1

    end if

!   compute positions

    if ( lorder == 'ND' ) then

!     nodal points and then degrees in inner loop

      if ( elementwise ) then

        dof = sum ( problem%vec_elnumdegfd(elgrp)%a(:,vec) )

!       positions in vector
        pos(1:dof) = [ (spe+i, i =1,dof) ]

      else

        dof = 0

        do node = 1, mesh%elnumnod(elgrp)

          nodenr = mesh%topology(elgrp)%a(node,elem)
          sp     = problem%vec_nodnumdegfd(nodenr,vec)
          nndof  = problem%vec_nodnumdegfd(nodenr+1,vec) - sp

          if ( llayer > 0 ) call adjust_for_layer

!         positions in vector
          pos(dof+1:dof+nndof) = [ (sp+i, i =1,nndof) ]

          dof = dof + nndof

        end do

      end if

    else if ( lorder == 'DN' ) then

!     degrees and then nodes in inner the loop

      if ( elementwise ) then

!       data stored elementwise

        dof = 0

        maxdeg = maxval(problem%vec_elnumdegfd(elgrp)%a(:,vec))

        do deg = 1, maxdeg

          do node = 1, mesh%elnumnod(elgrp)

!           positions

            sp    = spe + sum ( problem%vec_elnumdegfd(elgrp)%a(1:node-1,vec) )
            nndof = problem%vec_elnumdegfd(elgrp)%a(node,vec)

!           positions in vector
            if ( deg <= nndof ) then
              pos(dof+1) = sp+deg
              dof = dof + 1
            end if

          end do

        end do

      else

!       data stored in nodes

        dof = 0

        nodes = mesh%topology(elgrp)%a(:,elem)
        maxdeg = maxval( ( problem%vec_nodnumdegfd(nodes+1,vec) &
                               - problem%vec_nodnumdegfd(nodes,vec) &
                         ) / maxdegfac )

        do deg = 1, maxdeg

          do node = 1, mesh%elnumnod(elgrp)

            nodenr = mesh%topology(elgrp)%a(node,elem)

!           positions

            sp = problem%vec_nodnumdegfd(nodenr,vec)
            nndof = problem%vec_nodnumdegfd(nodenr+1,vec) &
                      - problem%vec_nodnumdegfd(nodenr,vec)

            if ( llayer > 0 ) call adjust_for_layer

!           positions in vector
            if ( deg <= nndof ) then
              pos(dof+1) = sp+deg
              dof = dof + 1
            end if

          end do

        end do

      end if

    end if

 contains

    subroutine adjust_for_layer

      if ( numl(node) > 1 ) then
        nndof = nndof / numl(node)
        sp = sp + nndof * pnuml(node)
      else if ( numl(node) == 0 ) then
        nndof = 0
      end if

    end subroutine adjust_for_layer

  end subroutine pos_array_vec


! Position array in vector for nodes

  subroutine pos_array_vec_node ( problem, nodenr, dof, pos, vec, layer )

    type(problem_t), intent(in)  :: problem
    integer, intent(in) :: nodenr

!   dof is the number of degrees of freedom involved
    integer, intent(out) :: dof

!   the position array
    integer, dimension(:), intent(out) :: pos

!   the vector number
    integer, intent(in) :: vec

!   if layer is present and layer > 0 the element degrees of freedom are
!   restricted to the specified layer.
!   layer=0 is identical to layer not present
    integer, intent(in), optional :: layer


!   this routine computes the positions of the degrees of freedom of a node

    integer :: llayer
    integer :: bp, nndof, i


    llayer = set_optional ( variable=layer, default=0 )

!   all degrees

    bp    = problem%vec_nodnumdegfd(nodenr,vec)
    nndof = problem%vec_nodnumdegfd(nodenr+1,vec) &
                - problem%vec_nodnumdegfd(nodenr,vec)

    if ( llayer > 0 ) call adjust_for_layer

    dof = nndof

!   positions in vector
    pos(:dof) = [(i,i=bp+1,bp+dof)]

  contains


!   modify start pointer and number of degrees for layer

    subroutine adjust_for_layer

      integer, dimension(1) :: numl, pnuml

      call find_layer_in_nodes ( problem, llayer, nodes=[nodenr], &
          numl=numl, pnuml=pnuml )

      if ( numl(1) > 1 ) then
        nndof = nndof / numl(1)
        bp = bp + nndof * pnuml(1)
      else if ( numl(1) == 0 ) then
        nndof = 0
      end if

    end subroutine adjust_for_layer

  end subroutine pos_array_vec_node


! Position array in vector for geometries

  subroutine pos_array_vec_geometry ( problem, geometry, elem, dof, pos, &
    vec, order, layer, lp )

    type(problem_t), intent(in)  :: problem
    type(geometry_t), intent(in) :: geometry
    integer, intent(in) :: elem

!   dof is the number of degrees of freedom involved
    integer, intent(out) :: dof

!   the position array
    integer, dimension(:), intent(out) :: pos

!   the vector number
    integer, intent(in) :: vec

!   The parameter order determines the sequence of the degrees of freedom on
!   elementlevel. There are only two possibilities and affect only the two loops
!   over nodal points and degrees of freedom
!     order = 'ND' : the most inner loop is over the degrees of freedom
!     order = 'DN' : the most inner loop is over the nodal points
!   If problem%numlayers > 0 the degrees are stacked into layers and we have
!         order = 'ND' :      NLD
!         order = 'DN' :      LDN
!   See also the userguide for further explanation.
!   The default of order is order = 'DN'
    character(len=*), intent(in), optional :: order

!   if layer is present and layer > 0 the element degrees of freedom are
!   restricted to the specified layer.
!   layer=0 is identical to layer not present
    integer, intent(in), optional :: layer

!   if lp is present and layer is present (in heading) and layer > 0, then
!   lp(i) = .true. if layer is present in i^th node else .false.
    logical, dimension(:), intent(out), optional :: lp

!   this routine computes the positions of the components of a vector inside
!   an element on a geometry


    integer :: nodenr, node, sp, nndof, deg, maxdeg, i, llayer
    integer, dimension(geometry%elnumnod) :: nodes, numl, pnuml, maxdegfac
    character(len=2) :: lorder

    llayer = set_optional ( variable=layer, default=0 )
    lorder = set_optional ( variable=order, default='DN' )

!   find layer in the nodes of the element

    if ( llayer > 0  ) then

      call find_layer_in_nodes ( problem, llayer, &
        nodes=geometry%topology(:,elem,2), numl=numl, pnuml=pnuml, lp=lp )

      maxdegfac = max ( numl, 1 )

    else

      maxdegfac = 1

    end if

!   compute positions

    if ( lorder == 'ND' ) then

!     nodal points and then degrees in inner loop

      dof = 0

      do node = 1, geometry%elnumnod

        nodenr = geometry%topology(node,elem,2)
        sp     = problem%vec_nodnumdegfd(nodenr,vec)
        nndof  = problem%vec_nodnumdegfd(nodenr+1,vec) - sp

        if ( llayer > 0 ) call adjust_for_layer

!       positions in vector
        pos(dof+1:dof+nndof) = [ (sp+i, i =1,nndof) ]

        dof = dof + nndof

      end do

    else if ( lorder == 'DN' ) then

!     degrees and then nodes in inner the loop

      dof = 0

      nodes = geometry%topology(:,elem,2)
      maxdeg = maxval( ( problem%vec_nodnumdegfd(nodes+1,vec) &
                             - problem%vec_nodnumdegfd(nodes,vec) &
                       ) / maxdegfac )

      do deg = 1, maxdeg

        do node = 1, geometry%elnumnod

          nodenr = geometry%topology(node,elem,2)

!         positions

          sp = problem%vec_nodnumdegfd(nodenr,vec)
          nndof = problem%vec_nodnumdegfd(nodenr+1,vec) &
                      - problem%vec_nodnumdegfd(nodenr,vec)

          if ( llayer > 0 ) call adjust_for_layer

!         positions in vector
          if ( deg <= nndof ) then
            pos(dof+1) = sp+deg
            dof = dof + 1
          end if

        end do

      end do

    end if

  contains

    subroutine adjust_for_layer

      if ( numl(node) > 1 ) then
        nndof = nndof / numl(node)
        sp = sp + nndof * pnuml(node)
      else if ( numl(node) == 0 ) then
        nndof = 0
      end if

    end subroutine adjust_for_layer

  end subroutine pos_array_vec_geometry


end module pos_array_m
