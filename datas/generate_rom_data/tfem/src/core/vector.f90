
! Copyright (C) 2004-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Routines, for the vectors

module vector_m

  use kind_defs_m
  use mesh_m
  use problem_defs_m
  use pos_array_m
  use system_defs_m
  use vector_defs_m
  use element_defs_m
  use set_optional_m

  implicit none


  logical, save :: firstwarning = .true.


! interface for generic create subroutine

  interface create
    module procedure create_vector, create_mvector, create_subscript_vector
  end interface create

  interface create_subscript
    module procedure create_subscript_vector
  end interface create_subscript


! interface for generic delete subroutine

  interface delete
    module procedure delete_vector, delete_mvector, delete_subscript_vector
  end interface delete


! interface for generic copy subroutine

  interface copy
    module procedure copy_vector, copy_mvector
  end interface copy

contains


! Create vector

  subroutine create_vector ( problem, vector, vec, physq, elementwise )

    type(problem_t), intent(in) :: problem
    type(vector_t), intent(inout) :: vector

!   if vec is present a vector is created with vector nr vec
!   NOTE: one of vec or physq must be present
    integer, intent(in), optional :: vec

!   if physq is present a vector is created with vector nr defining the
!   physical quantity physq (=problem%phys(physq)).
!   NOTE: one of vec or physq must be present
    integer, intent(in), optional :: physq

!   if elementwise is present and elementwise=.true. the vector will be
!   stored elementwise: for each element the values in the nodes of the
!   element are stored separately element by element in the vector.
!   If such a vector is used in the routine derive_vector the elementvectors
!   are stored separately and are not averaged over the nodes.
!   In this way discontinuous properties, such as gradients can be plotted
!   or sampled without preaveraging.
!   NOTE: layers are not available for this structure.
    logical, intent(in), optional :: elementwise


    call check ( problem, 'create_vector' )

    if ( vector%created ) then
      write(*,'(/a/)') &
        'Error in create_vector: vector has already been created '
      stop
    end if

    if ( present(vec) .and. present(physq) ) then
      write(*,'(/a/)') &
        'Error in create_vector: both vec and physq present'
      stop
    else if ( .not. present(vec) .and. .not. present(physq) ) then
      write(*,'(/a/)') &
        'Error in create_vector: either vec or physq must be present'
      stop
    end if

    if ( present(vec) ) then
      if ( vec < 1 .or. vec > problem%nvec ) then
        write(*,'(/a/)') &
          'Error in create_vector: vec out of range'
        stop
      end if
    end if

    if ( present(physq) ) then
      if ( problem%nphysq == 0 ) then
        write(*,'(/a/a/)') &
          'Error in create_vector: physq specified, however there are no ', &
          ' physical quantities defined'
        stop
      end if
      if ( physq < 1 .or. physq > problem%nphysq ) then
        write(*,'(/a/)') &
          'Error in create_vector: physq out of range'
        stop
      end if
    end if

    vector%probnr = problem%probnr

    if ( present(vec) ) then
      vector%vec = vec
    else if ( present(physq) ) then
      vector%vec = problem%physq(physq)
    end if

!   elementwise ?

    if ( present(elementwise) ) then
      if ( elementwise ) then
        vector%n = problem%vec_numdegfd(vector%vec,2)
      else
        vector%n = problem%vec_numdegfd(vector%vec,1)
      end if
      vector%elementwise = elementwise
    else
      vector%n = problem%vec_numdegfd(vector%vec,1)
      vector%elementwise = .false.
    end if

    allocate(vector%u(vector%n))

    vector%created = .true.

  end subroutine create_vector


! Create multiple vector

  subroutine create_mvector ( problem, mvector, vec, physq, elementwise )

    type(problem_t),   intent(in)  :: problem
    type(vector_t), dimension(:), intent(inout) :: mvector

!   if vec is present a vector is created with vector nr vec
!   NOTE: one of vec or physq must be present
    integer, intent(in), optional :: vec

!   if physq is present a vector is created with vector nr defining the
!   physical quantity physq (=problem%phys(physq)).
!   NOTE: one of vec or physq must be present
    integer, intent(in), optional :: physq

!   if elementwise is present and elementwise=.true. the vector will be
!   stored elementwise: for each element the values in the nodes of the
!   element are stored separately element by element in the vector.
!   If such a vector is used in the routine derive_vector the elementvectors
!   are stored separately and are not averaged over the nodes.
!   In this way discontinuous properties, such as gradients can be plotted
!   or sampled without preaveraging.
    logical, intent(in), optional :: elementwise

    integer :: ivec

    if ( any(mvector%created) ) then
      write(*,'(/3(a/)/)') &
        'Error in create_mvector:', &
        ' some vectors in mvector have already been created ', &
        ' mvector cannot be created '
      stop
    end if

    do ivec = 1, size(mvector)
      call create_vector ( problem, mvector(ivec), vec, physq, elementwise )
    end do

  end subroutine create_mvector


! Delete single vector

  subroutine delete_single_vector ( vector, nr )

    type(vector_t), intent(inout) :: vector
    integer, intent(in) :: nr

    if ( .not. vector%created ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error in delete_vector:', &
        ' Vector has not been created and cannot be deleted ', &
        ' Vector number in heading = ', nr
      stop
    end if

    if ( allocated(vector%w) ) then
      write(*,'(/4(a/),a,i0/)') &
        'Warning in delete_vector:', &
        ' the vector with weights w is still there, whereas it should ', &
        ' have been removed with the finalize step in derive_vector.', &
        ' Deallocating it anyway. ', &
        ' Vector number in heading = ', nr
    end if

!   deallocate all allocatables and intialize to defaults
    call deall ( vector )

  contains

    subroutine deall ( vector )
      type(vector_t), intent(out) :: vector
    end subroutine deall

  end subroutine delete_single_vector


! Delete vector

  subroutine delete_vector ( vector1, vector2, vector3, vector4, vector5, &
                             vector6, vector7, vector8, vector9, vector10 )

    type(vector_t), intent(inout) :: vector1
    type(vector_t), intent(inout), optional :: vector2, vector3, vector4, &
      vector5, vector6, vector7, vector8, vector9, vector10

    call delete_single_vector(vector1,1)
    if ( present(vector2) ) call delete_single_vector(vector2,2)
    if ( present(vector3) ) call delete_single_vector(vector3,3)
    if ( present(vector4) ) call delete_single_vector(vector4,4)
    if ( present(vector5) ) call delete_single_vector(vector5,5)
    if ( present(vector6) ) call delete_single_vector(vector6,6)
    if ( present(vector7) ) call delete_single_vector(vector7,7)
    if ( present(vector8) ) call delete_single_vector(vector8,8)
    if ( present(vector9) ) call delete_single_vector(vector9,9)
    if ( present(vector10) ) call delete_single_vector(vector10,10)

  end subroutine delete_vector


! Delete single mvector

  subroutine delete_single_mvector ( mvector, nr )

    type(vector_t), dimension(:), intent(inout) :: mvector
    integer, intent(in) :: nr

!   deallocate all allocatables and intialize to defaults
    call deall ( mvector )

  contains

    subroutine deall ( mvector )
      type(vector_t), dimension(:), intent(out) :: mvector
    end subroutine deall

  end subroutine delete_single_mvector


! Delete mvector

  subroutine delete_mvector ( mvector1, mvector2, mvector3, mvector4, mvector5 )

    type(vector_t), dimension(:), intent(inout) :: mvector1
    type(vector_t), dimension(:), intent(inout), optional :: &
      mvector2, mvector3, mvector4, mvector5

    call delete_single_mvector(mvector1,1)
    if ( present(mvector2) ) call delete_single_mvector(mvector2,2)
    if ( present(mvector3) ) call delete_single_mvector(mvector3,3)
    if ( present(mvector4) ) call delete_single_mvector(mvector4,4)
    if ( present(mvector5) ) call delete_single_mvector(mvector5,5)

  end subroutine delete_mvector


! Copy vector to another vector

  subroutine copy_vector ( vector1, vector2 )

    type(vector_t), intent(inout) :: vector1, vector2

!   This routine performs a real data copy of one vector to the other using
!   the assignment statement:
!
!      vector2 = vector1
!

    if ( .not. vector1%created ) then
      write(*,'(/a/)') &
        'Error in copy_vector:  vector1 has not been created '
      stop
    end if

    vector2 = vector1

  end subroutine copy_vector


! Copy mvector to another mvector

  subroutine copy_mvector ( mvector1, mvector2 )

    type(vector_t), dimension(:), intent(inout) :: mvector1, mvector2

!   This routine performs a real data copy of one vector to the other using
!   the assignment statement:
!
!      vector2 = vector1
!

    integer :: ivec

    if ( size(mvector1) /= size(mvector2) ) then
      write(*,'(/2a/)') &
        'Error in copy_mvector: ', &
        ' dimension of mvector1 and mvector2 is different '
      stop
    end if

    do ivec = 1, size(mvector1)
      call copy_vector ( mvector1(ivec), mvector2(ivec) )
    end do


  end subroutine copy_mvector


! Get element values/index of vector for the element given by (elgrp,elem)

  subroutine get_vector ( mesh, problem, vector, elgrp, elem, u, order, layer, &
    posu, ndofu, sloppy )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(vector_t), intent(in) :: vector
    integer, intent(in) :: elgrp, elem

!   the element degrees of freedom, if present
    real(dp), intent(out), dimension(:), optional :: u

!   The parameter order determines the sequence of the degrees of freedom on
!   elementlevel
!   There are two possibibilities:
!     order = 'ND' : outside loop: nodal points
!                    inner loop: the degrees of freedom
!     order = 'DN' : outside loop: degrees of freedom
!                    inner loop: nodal points
!   The default is order = 'DN'
    character(len=*), intent(in), optional :: order

!   if layer is present and layer > 0 the element degrees of freedom are
!   restricted to the specified layer.
!   layer=0 is identical to layer not present
    integer, intent(in), optional :: layer

!   if present it contains the positions of the components of u in vector
    integer, dimension(:), intent(out), optional :: posu

!   if present it gives the number of element degrees of freedom stored in u
    integer, intent(out), optional :: ndofu

!   sloppy? Don't check nelgrp
    logical, intent(in), optional :: sloppy

    logical, dimension(mesh%elnumnod(elgrp)) :: lp
    integer :: dof, ndof, llayer
    integer, dimension(:), allocatable :: pos
    integer, dimension(mesh%elnumnod(elgrp)) :: nodes

!   testing

    call check ( problem, 'get_vector', mesh, sloppy )

    if ( .not. vector%created ) then
      write(*,'(/a/)') &
        'Error in get_vector: vector not created.'
      stop
    end if

    if ( problem%probnr /= vector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in get_vector: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in vector = ', vector%probnr
      stop
    end if

!   ordering; set defaults

    if ( present(order) ) then
      if ( all ( order /= [ 'ND', 'DN' ] ) ) then
        write(*,'(2(/a)/)') &
          'Error in get_vector: ', &
          ' heading parameter order must be either ''ND'' or ''DN'''
        stop
      end if
    end if

!   layers

    llayer = set_optional ( variable=layer, default=0 )

    if ( present(layer) ) then
      if ( layer < 0 .or. layer > problem%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in get_vector: layer out of range.'
        stop
      end if
    end if

!   number of degrees of freedom in element

    nodes = mesh%topology(elgrp)%a(:,elem)

    if ( llayer > 0 ) then

!     single layer

      lp = btest(problem%nodlayers(nodes),llayer-1)

      ndof = sum ( problem%vec_elnumdegfd(elgrp)%a(:,vector%vec), mask=lp )

    else

!     standard situation

      ndof = sum ( problem%vec_nodnumdegfd( nodes+1, vector%vec ) &
                          - problem%vec_nodnumdegfd( nodes, vector%vec ) )

    end if

    if ( present(ndofu) ) ndofu = ndof

!   check u array

    if ( present(u) ) then
      if ( size(u) < ndof ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_vector: ', &
          ' u array is too small: size(u) = ', size(u), &
          ' whereas the number of degrees of freedom = ', ndof
        stop
      end if
    end if

!   check posu array

    if ( present(posu) ) then
      if ( size(posu) < ndof ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_vector: ', &
          ' posu array is too small: size(posu) = ', size(posu), &
          ' whereas the number of degrees of freedom = ', ndof
        stop
      end if
    end if

    allocate(pos(ndof))

!   compute positions in vector

    call pos_array_vec ( mesh, problem, elgrp, elem, dof, pos, vector%vec, &
      vector%elementwise, order, layer )

    if ( dof /= ndof ) then
      write(*,'(/a/a/)') &
        'Internal error in get_vector: ', &
        ' dof /= ndof '
      stop
    end if

!   compute u

    if ( present(u) ) u(1:ndof) = vector%u(pos)

!   posu array

    if ( present(posu) ) posu(1:ndof) = pos

    deallocate(pos)

  end subroutine get_vector


! Get element values/index of vector for the element on a geometry given by
! either (curve,elem), (surface,elem), (volume,elem) or (ndimr,geometry,elem).

  subroutine get_vector_geometry ( mesh, problem, vector, elem, u, curve, &
    surface, volume, ndimr, geometry, order, layer, posu )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(vector_t), intent(in) :: vector
    integer, intent(in) :: elem

!   the element degrees of freedom, if present
    real(dp), intent(out), dimension(:), optional :: u

!   the geometry number: one of curve, surface, volume must be present (legacy)
    integer, intent(in), optional :: curve, surface, volume

!   dimension of reference space (ndimr) and geometry number
    integer, intent(in), optional :: ndimr, geometry

!   The parameter order determines the sequence of the degrees of freedom on
!   elementlevel
!   There are two possibibilities:
!     order = 'ND' : outside loop: nodal points
!                    inner loop: the degrees of freedom
!     order = 'DN' : outside loop: degrees of freedom
!                    inner loop: nodal points
!   The default is order = 'DN'
    character(len=*), intent(in), optional :: order

!   if layer is present and layer > 0 the element degrees of freedom are
!   restricted to the specified layer.
!   layer=0 is identical to layer not present
    integer, intent(in), optional :: layer

!   if present it contains the positions of the components of u in vector
    integer, dimension(:), intent(out), optional :: posu


    integer :: dof, ndof, llayer
    integer, dimension(:), allocatable :: pos, nodes, nl
    type(geometry_t) :: geom

    integer :: lcurve, lsurface, lvolume

!   traditional interface (legacy)

    lcurve = set_optional ( variable=curve, default=0 )
    lsurface = set_optional ( variable=surface, default=0 )
    lvolume = set_optional ( variable=volume, default=0 )

!   interface with dimension of reference space of geometry

    if ( present(ndimr) .and. present(geometry) ) then
      select case (ndimr)
        case(1); lcurve = geometry
        case(2); lsurface = geometry
        case(3); lvolume = geometry
        case default
          call errormsg_case_default ( 'get_vector_geometry', 'ndimr', &
            int_value=ndimr )
      end select
    end if

    llayer = set_optional ( variable=layer, default=0 )

!   testing

    call check ( problem, 'get_vector_geometry', mesh )

!   ordering; set defaults

    if ( present(order) ) then
      if ( all ( order /= [ 'ND', 'DN' ] ) ) then
        write(*,'(2(/a)/)') &
          'Error in get_vector_geometry: ', &
          ' heading parameter order must be either ''ND'' or ''DN'''
        stop
      end if
    end if

    if ( lcurve /= 0 ) then

      if ( lcurve < 1 .or. lcurve > mesh%ncurves ) then

        write(*,'(/a/a,i0,/a,i0/)') &
          'Error: curve in the heading of get_vector_geometry is ', &
          'out of range. curve is ', lcurve, &
          'whereas the number of curves is ', mesh%ncurves
        stop

      end if

      allocate(nodes(mesh%curves(lcurve)%elnumnod))

    else if ( lsurface /= 0 ) then

      if ( lsurface < 1 .or. lsurface > mesh%nsurfaces ) then

        write(*,'(/a/a,i0,/a,i0/)') &
          'Error: surface in the heading of get_vector_geometry is ', &
          'out of range. surface is ', lsurface, &
          'whereas the number of surfaces is ', mesh%nsurfaces
        stop

      end if

      allocate(nodes(mesh%surfaces(lsurface)%elnumnod))

    else if ( lvolume /= 0 ) then

      if ( lvolume < 1 .or. lvolume > mesh%nvolumes ) then

        write(*,'(/a/a,i0,/a,i0/)') &
          'Error: volume in the heading of get_vector_geometry is ', &
          'out of range. volume is ', lvolume, &
          'whereas the number of volumes is ', mesh%nvolumes
        stop

      end if

      allocate(nodes(mesh%volumes(lvolume)%elnumnod))

    else

      write(*,'(/a/a/)') &
        'Error: either curve, surface, volume or (ndimr,geometry)', &
        ' must be present in the heading of get_vector_geometry'
      stop

    end if

    if ( .not. vector%created ) then
      write(*,'(/a/)') &
        'Error in get_vector_geometry: vector not created.'
      stop
    end if

    if ( problem%probnr /= vector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in get_vector_geometry: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in vector = ', vector%probnr
      stop
    end if

!   check vector

    if ( vector%elementwise ) then
      write(*,'(/2(a/))') &
        'Error in get_vector_geometry: a vector with elementwise ', &
        ' stored data has not been implemented. '
      stop
    end if

!   layers

    if ( present(layer) ) then
      if ( layer < 0 .or. layer > problem%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in get_vector_geometry: layer out of range.'
        stop
      end if
    end if

!   find geometry

    if ( lcurve > 0 ) then
      nodes = mesh%curves(lcurve)%topology(:,elem,2)
      geom = mesh%curves(lcurve)
    else if ( lsurface > 0 ) then
      nodes = mesh%surfaces(lsurface)%topology(:,elem,2)
      geom = mesh%surfaces(lsurface)
    else if ( lvolume > 0 ) then
      nodes = mesh%volumes(lvolume)%topology(:,elem,2)
      geom = mesh%volumes(lvolume)
    end if


!   number of degrees of freedom in element

    if ( llayer > 0 ) then

!     single layer

      allocate ( nl(size(nodes)) )

      where ( btest(problem%nodlayers(nodes),llayer-1) )
        nl = count_layers ( problem%nodlayers(nodes), problem%numlayers )
      else where
        nl = 0
      end where

      ndof = sum ( ( problem%vec_nodnumdegfd( nodes+1, vector%vec ) &
                     - problem%vec_nodnumdegfd( nodes, vector%vec ) ) / nl, &
                     mask = nl > 0 )

    else

!     standard situation

      ndof = sum ( problem%vec_nodnumdegfd( nodes+1, vector%vec ) &
                   - problem%vec_nodnumdegfd( nodes, vector%vec ) )

    end if

!   check u array

    if ( present(u) ) then
      if ( size(u) < ndof ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_vector_geometry: ', &
          ' u array is too small: size(u) = ', size(u), &
          ' whereas the number of degrees of freedom = ', ndof
        stop
      end if
    end if

!   check posu array

    if ( present(posu) ) then
      if ( size(posu) < ndof ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_vector_geometry: ', &
          ' posu array is too small: size(posu) = ', size(posu), &
          ' whereas the number of degrees of freedom = ', ndof
        stop
      end if
    end if

    allocate(pos(ndof))

!   compute positions in vector

    call pos_array_vec_geometry ( problem, geom, elem, dof, pos, &
      vector%vec, order, layer )

    if ( dof /= ndof ) then
      write(*,'(/a/a/)') &
        'Internal error in get_vector_geometry: ', &
        ' dof /= ndof '
      stop
    end if

!   compute u

    if ( present(u) ) u(1:ndof) = vector%u(pos)

!   posu array

    if ( present(posu) ) posu(1:ndof) = pos

    deallocate( pos, nodes )

  end subroutine get_vector_geometry


! Derive vector

  subroutine derive_vector ( mesh, problem, vector, elemsub, coefficients, &
    mcoefficients, oldvectors, order, elgroup1, elgroup2, groups, elpoints, &
    elcurves, elsurfaces, layer, addvec, averaging, finalize )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the vector to be derived
    type(vector_t), intent(inout) :: vector

!   this is the element subroutine that must be supplied by the calling routine
    interface
      subroutine elemsub ( mesh, problem, elgrp, elem, first, last, &
        coefficients, oldvectors, elemvec, elemwts )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: elgrp, elem
        logical, intent(in) :: first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:) :: elemvec, elemwts
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

!   The parameter order determines the sequence of the degrees of freedom on
!   elementlevel
!   There are two possibibilities:
!     order = 'ND' : outside loop: nodal points
!                    inner loop: the degrees of freedom
!     order = 'DN' : outside loop: degrees of freedom
!                    inner loop: nodal points
!   The default is order = 'DN'
    character(len=*), intent(in), optional :: order

!   if these are present the assembling takes place for element groups
!   elgroup1,...,elgroup2 only. If only elgroup1 is present one group is
!   assembled only.
    integer, intent(in), optional :: elgroup1, elgroup2

!   if present: the element groups to be assembled
!   For example groups=(/2,4/) will assemble for element groups 2 and 4.
!   Default: all groups
    integer, dimension(:), intent(in), optional :: groups

!   if one or more of these are present the assembling takes place only
!   for elements CONNECTED to the specified points, curves and/or surfaces.
!   In this way the computations can be much faster if, for example, only
!   the data on a curve is needed for drag computations.
!   NOTE: the values in the nodes of the geometrical object
!   (point,curve,surface) will be the same as when all elements contributions
!   are added and averaged. For the other nodes in the elements connected to
!   the geometrical object this might not be true!
!   NOTE: the first and last element in a group are _always_ included.
    integer, dimension(:), intent(in), optional :: elpoints, elcurves, &
      elsurfaces

!   if present the assembling takes place within the single layer of degrees of
!   freedom only. Only elements where all nodes have degrees in the specified
!   layer are assembled.
    integer, intent(in), optional :: layer

!   if addvec is set to .true. the vectors are not cleared before
!   the assembling and thus the element vectors are added to an
!   existing vector. The default is .false. (clearing)
    logical, intent(in), optional :: addvec

!   if averaging is set to .true. the weight factors for the averaging
!   are assembled and element vectors are added. The actual averaging is only
!   performed if finalize=.true. The weighting factors are stored in vector.
!   When averaging=.false., weighting factors
!   are not assembled and element vectors are just substituted in the big
!   vector, possibly destroying previous values from other elements.
!   When averaging is absent .true. is assumed.
!   NOTE: if vector%elementwise=.true. averaging is not performed.
    logical, intent(in), optional :: averaging

!   if finalize is set to .true. the averaging step takes place and weighting
!   factors in vector are destroyed.
!   Default is .true.
!   NOTE: if vector%elementwise=.true. finalizing is not performed.
    logical, intent(in), optional :: finalize


!   This routine performs the assembly proces and nothing more. This means that
!   after the first call (with addvec=.false., finalize=.false.),
!   the routine can be called as many times as needed with addvec=.true. The
!   final call needs finalize=.true. to compute the averaged quantities.


    logical :: addvector, first, last, aver, fina, allgroups, allelements
    logical, allocatable, dimension(:) :: wkl
    integer :: elem, elgrp, ndof, dof, grp
    integer :: crv, curve, srf, surface
    integer :: row, rowg
    integer :: lgroups(mesh%nelgrp), lnelgrp, i

    integer, allocatable, dimension(:) :: pos
    integer, allocatable, dimension(:,:) :: elements
    real(dp), allocatable, dimension(:) :: elemvec, elemwts
    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl

    allocate ( wkl(mesh%nnodes) )


!   testing

    call check ( mesh, 'derive_vector' )
    call check ( problem, 'derive_vector', mesh )

    if ( present(mcoefficients) ) then
      if ( size(mcoefficients) /= mesh%nelgrp ) then
        write(*,'(/2(a/))') &
          'Error in derive_vector: ', &
          '   size of mcoefficients /= number of element groups'
        stop
      end if
    end if

    if ( .not. vector%created ) then
      write(*,'(/a/)') &
        'Error in derive_vector: vector not created.'
      stop
    end if

    if ( problem%probnr /= vector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in derive_vector: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in vector = ', vector%probnr
      stop
    end if

!   ordering; set defaults

    if ( present(order) ) then
      if ( all ( order /= [ 'ND', 'DN' ] ) ) then
        write(*,'(2(/a)/)') &
          'Error in derive_vector: ', &
          ' heading parameter order must be either ''ND'' or ''DN'''
        stop
      end if
    end if

!   layers

    if ( present(layer) ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in derive_vector: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in derive_vector: layer out of range.'
        stop
      end if
      if ( vector%elementwise ) then
        write(*,'(/2(a/))') &
          'Error in derive_vector:', &
          '  vector%elementwise=.true. is incompatible with layers'
      end if
    else if ( problem%numlayers > 0 ) then
      write(*,'(3(/a)/)') &
        'Error in derive_vector: ', &
        ' layers have been defined and the layer keyword is not present', &
        ' layers not implemented for this routine without the layer keyword'
      stop
    end if

!   initialization

    if ( present(oldvectors) ) oldvl = oldvectors
    if ( present(coefficients) ) coeffl = coefficients

!   which element groups?

    if ( present(elgroup1) .and. present(elgroup2) ) then
!     specified range of groups only
      lnelgrp = elgroup2 - elgroup1 + 1
      lgroups(1:lnelgrp) = [ (i,i=elgroup1,elgroup2) ]
      allgroups = .false.
    else if ( present(elgroup1) ) then
!     one group only
      lnelgrp = 1
      lgroups(1) = elgroup1
      allgroups = .false.
    else if ( present(groups) ) then
      lnelgrp = size(groups)
      lgroups(1:lnelgrp) = groups
      allgroups = .false.
    else
!     all groups
      lnelgrp = mesh%nelgrp
      lgroups = [ (i,i=1,lnelgrp) ]
      allgroups = .true.
    end if

    if ( any ( lgroups(1:lnelgrp) < 1 ) .or. &
         any ( lgroups(1:lnelgrp) > mesh%nelgrp ) ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error: elgroup1, elgroup2, or groups in the heading of ', &
        ' derive_vector is out of range: ', &
        ' some groups are < 1 or larger than the number of groups ', mesh%nelgrp
      stop
    end if

    do grp = 1, lnelgrp

      elgrp = lgroups(grp)

!     skip inactive groups
      if ( any( elgrp == problem%inactivegroups ) ) cycle

      if ( all( problem%vec_elnumdegfd(elgrp)%a(:,vector%vec) == 0 ) ) then

        write(*,'(/a/a,i0/)') &
          'Error: element group in the heading of derive_vector has ', &
          'no degrees of freedom. elgrp is ', elgrp
        stop

      end if

    end do

!   which elements?

    if ( present(elpoints) .or. present(elcurves) &
         .or. present(elsurfaces) .or. present(layer) ) then

!     specified elements only

      allocate(elements(maxval(mesh%grpnumel),mesh%nelgrp))
      elements(:,lgroups(1:lnelgrp)) = 0

!     points

      if ( present(elpoints) ) then

        if ( any ( elpoints <=0 ) .or. any ( elpoints > mesh%npoints ) ) then
          write(*,'(/2(a/),a,i0/)') &
            'Error: elpoints in the heading of derive_vector is ', &
            ' out of range. The point numbers must not be smaller than 1 or ', &
            ' larger than the number of points in the mesh = ', mesh%npoints
          stop
        end if

!       set elements to 1 that need to be included
        call set_elements ( mesh%points(elpoints) )

      end if

!     curves

      if ( present(elcurves) ) then

        if ( any ( elcurves <=0 ) .or. any ( elcurves > mesh%ncurves ) ) then
          write(*,'(/2(a/),a,i0/)') &
            'Error: elcurves in the heading of derive_vector is ', &
            ' out of range. The curve numbers must not be smaller than 1 or ', &
            ' larger than the number of curves in the mesh = ', mesh%ncurves
          stop
        end if

!       set elements to 1 that need to be included
        do crv = 1, size(elcurves)
          curve = elcurves(crv)
          call set_elements ( mesh%curves(curve)%nodes )
        end do

      end if

!     surfaces

      if ( present(elsurfaces) ) then
        if ( any ( elsurfaces <=0 ) .or. &
          any ( elsurfaces > mesh%nsurfaces ) ) then
          write(*,'(/2(a/),a,i0/)') &
            'Error: elsurfaces in the heading of derive_vector is out of ', &
            ' range. The surface numbers must not be smaller than 1 or ', &
            ' larger than the number of surfaces in the mesh = ', mesh%nsurfaces
          stop
        end if

!       set elements to 1 that need to be included
        do srf = 1, size(elsurfaces)
          surface = elsurfaces(srf)
          call set_elements ( mesh%surfaces(surface)%nodes )
        end do

      end if

!     layer

      if ( present(layer) ) then

        wkl = btest(problem%nodlayers,layer-1)

        do grp = 1, lnelgrp

          elgrp = lgroups(grp)

!         skip inactive groups
          if ( any( elgrp == problem%inactivegroups ) ) cycle

          do elem = 1, mesh%grpnumel(elgrp)
            if ( all( wkl(mesh%topology(elgrp)%a(:,elem)) ) ) then
!             element full within layer
              elements(elem,elgrp) = 1
            end if
          end do

        end do

      end if

      allelements = .false.

    else

!     all elements
      allelements = .true.

    end if

    if ( present(addvec) ) then
      addvector = addvec
    else
      addvector = .false.
    end if

    if ( vector%elementwise ) then
      aver = .false.
      fina = .false.
    else
      if ( present(averaging) ) then
        aver = averaging
      else
        aver = .true.
      end if
      if ( present(finalize) ) then
        fina = finalize
      else
        fina = .true.
      end if
    end if

    if ( .not. addvector ) then

!     clear vector
      vector%u = 0

      if (aver) then

        if ( allocated(vector%w) ) then
          write(*,'(/4(a/))') &
            'Warning in derive_vector:', &
            ' the vector with weights w is still there, whereas it should ', &
            ' have been removed with the finalize step in derive_vector.', &
            ' Continuing anyway. '
        else
          allocate(vector%w(vector%n))
        end if

!       clear weights
        vector%w = 0

      end if

    else if ( aver .and. .not. allocated(vector%w) ) then

      write(*,'(/4(a/))') &
        'Error in derive_vector:', &
        ' the vector with weights w is not there, whereas it should ', &
        ' have been created in a previous step with addvec=.false. and', &
        ' finalize=.false.'
      stop

    end if

!   start loop over all elements

    do grp = 1, lnelgrp

      elgrp = lgroups(grp)

!     skip inactive groups
      if ( any( elgrp == problem%inactivegroups ) ) cycle

      if ( present(mcoefficients) ) coeffl = mcoefficients(elgrp)

!     find number of degrees of freedom in an element of this group

      ndof = sum ( problem%vec_elnumdegfd(elgrp)%a(:,vector%vec) )

!     reserve memory for element vector, element weights and positions

      allocate ( elemvec(ndof), elemwts(ndof) )
      allocate ( pos(ndof) )

!     loop over elements in this group

      do elem = 1, mesh%grpnumel(elgrp)

!       skip element?

        if ( .not. allelements ) then
!         avoid leaving allocated memory in elements
          first = .true.; last = .true.
          if ( elements(elem,elgrp) == 0 ) cycle
        else
          first = elem == 1
          last  = elem == mesh%grpnumel(elgrp)
        end if

!       compute element matrix and vector

        call elemsub ( mesh, problem, elgrp, elem, first, last, coeffl, oldvl, &
          elemvec, elemwts )

!       compute positions in vector of element degrees of freedom

        call pos_array_vec ( mesh, problem, elgrp, elem, dof, pos, &
          vector%vec, vector%elementwise, order, layer )

        if ( dof /= ndof ) then
          write(*,'(/a/a/)') &
            'Internal error in derive_vector: ', &
            ' dof /= ndof '
          stop
        end if

!       add element vector to large vector

        do row = 1, dof

          rowg = pos(row) ! global row number

          if ( aver ) then
            vector%u(rowg) = vector%u(rowg) + elemvec(row)
            vector%w(rowg) = vector%w(rowg) + elemwts(row)
          else
            vector%u(rowg) = elemvec(row)
          end if

        end do

      end do

      deallocate ( elemvec, elemwts, pos )

    end do

!   Finalize?

    if ( fina ) then

      if ( allgroups .and. allelements ) then
!       all groups and elements have been processed
        vector%u = vector%u / vector%w
      else
!       not all groups or elements have been processed: zero weights possible
        where ( vector%w /= 0 )
          vector%u = vector%u / vector%w
        end where
      end if

      deallocate(vector%w)

    end if

!   remove memory

    if ( .not. allelements ) deallocate ( elements )

    deallocate ( wkl )

  contains

!   set elements array to 1 for elements connected to the nodes in nodes

    subroutine set_elements ( nodes )

      integer, dimension(:), intent(in) :: nodes

      integer :: node, nodenr, start, numel, elm, elem, elgrp

      do node = 1, size(nodes)
        nodenr = nodes(node)
        start = mesh%nodnumel(nodenr)
        numel = mesh%nodnumel(nodenr+1) - start
        do elm = 1, numel
          elgrp = mesh%nodelem( start + elm, 1 )
          elem  = mesh%nodelem( start + elm, 2 )
        end do
        elements(elem,elgrp) = 1
      end do

    end subroutine set_elements

  end subroutine derive_vector


! extract a physical quantity from the system vector

  subroutine extract_physvector ( mesh, problem, sysvector, physvector )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the sysvector from which physvector must be extracted
    type(sysvector_t), intent(inout) :: sysvector

!   the vector to be extracted
    type(vector_t), intent(inout) :: physvector


    integer :: iq, physq, sp1, sp2, nodenr, elgrp, elem, node, ndof


    call check ( mesh, 'extract_physvector' )
    call check ( problem, 'extract_physvector', mesh )

    if ( .not. sysvector%created ) then
      write(*,'(/a/)') &
        'Error in extract_physvector: sysvector has not been created'
      stop
    end if
    if ( .not. physvector%created ) then
      write(*,'(/a/)') &
        'Error in extract_physvector: physvector has not been created'
      stop
    end if
    if ( problem%probnr /= sysvector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in extract_physvector: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in sysvector = ', sysvector%probnr
      stop
    end if
    if ( problem%probnr /= physvector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in extract_physvector: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in physvector = ', physvector%probnr
      stop
    end if
    if ( all ( problem%physq /= physvector%vec ) ) then
      write(*,'(/a/a/)') &
        'Error in extract_physvector: physvector in heading is not a vector', &
        ' of a physical quantity.'
      stop
    end if

!   which physical quantity?

    do iq = 1, problem%nphysq
      if ( problem%physq(iq) == physvector%vec ) then
        physq = iq
        exit
      end if
    end do

!   extract vector

    if ( physvector%elementwise ) then

!     store elementwise

      sp2 = 0

      do elgrp = 1, problem%nelgrp
        do elem = 1, mesh%grpnumel(elgrp)
          do node = 1, mesh%elnumnod(elgrp)

            nodenr = mesh%topology(elgrp)%a(node,elem)

            sp1 = problem%nodnumdegfd(nodenr) + &
                     sum( problem%vec_elnumdegfd(elgrp)%a(node,&
                                    &problem%physq(1:physq-1)) )

            ndof = problem%vec_elnumdegfd(elgrp)%a(node,physvector%vec)

            physvector%u(sp2+1:sp2+ndof) = &
                    sysvector%u(problem%degfdperm(sp1+1:sp1+ndof,2))

            sp2 = sp2 + ndof

          end do
        end do
      end do

    else

!     store in nodes

      do nodenr = 1, mesh%nnodes

        sp1 = problem%nodnumdegfd(nodenr) + &
                        sum( problem%vec_nodnumdegfd(nodenr+1,&
                                    &problem%physq(1:physq-1)) - &
                             problem%vec_nodnumdegfd(nodenr,&
                                    &problem%physq(1:physq-1)) )

        sp2 = problem%vec_nodnumdegfd(nodenr,physvector%vec)

        ndof = problem%vec_nodnumdegfd(nodenr+1,physvector%vec) &
                 - problem%vec_nodnumdegfd(nodenr,physvector%vec)

        physvector%u(sp2+1:sp2+ndof) = &
                sysvector%u(problem%degfdperm(sp1+1:sp1+ndof,2))


      end do

    end if

  end subroutine extract_physvector


! store a physical quantity in the system vector

  subroutine store_physvector ( mesh, problem, sysvector, physvector )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the sysvector in which physvector must be stored
    type(sysvector_t), intent(inout) :: sysvector

!   the vector to be put in sysvector
    type(vector_t), intent(inout) :: physvector


    integer :: iq, physq, sp1, sp2, nodenr, elgrp, elem, node, ndof


    call check ( mesh, 'store_physvector' )
    call check ( problem, 'store_physvector', mesh )

    if ( .not. sysvector%created ) then
      write(*,'(/a/)') &
        'Error in store_physvector: sysvector has not been created'
      stop
    end if
    if ( .not. physvector%created ) then
      write(*,'(/a/)') &
        'Error in store_physvector: physvector has not been created'
      stop
    end if
    if ( problem%probnr /= sysvector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in store_physvector: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in sysvector = ', sysvector%probnr
      stop
    end if
    if ( problem%probnr /= physvector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in store_physvector: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in physvector = ', physvector%probnr
      stop
    end if
    if ( all ( problem%physq /= physvector%vec ) ) then
      write(*,'(/a/a/)') &
        'Error in store_physvector: physvector in heading is not a vector', &
        ' of a physical quantity.'
      stop
    end if

!   which physical quantity?

    do iq = 1, problem%nphysq
      if ( problem%physq(iq) == physvector%vec ) then
        physq = iq
        exit
      end if
    end do

!   store vector in sysvector

    if ( physvector%elementwise ) then

!     vector stored elementwise

      if ( firstwarning ) then

        write(*,'(/4(a/))') &
        'Warning in store_physvector: physvector in heading has been stored ',&
        ' elementwise. This means that multiple data values in common nodes ',&
        ' are lost in the transfer to the sysvector.', &
        ' The data in the element with the largest number will be stored'

        firstwarning = .false.

      end if

      sp2 = 0

      do elgrp = 1, problem%nelgrp
        do elem = 1, mesh%grpnumel(elgrp)
          do node = 1, mesh%elnumnod(elgrp)

            nodenr = mesh%topology(elgrp)%a(node,elem)

            sp1 = problem%nodnumdegfd(nodenr) + &
                     sum( problem%vec_elnumdegfd(elgrp)%a(node,&
                                    &problem%physq(1:physq-1)) )

            ndof = problem%vec_elnumdegfd(elgrp)%a(node,physvector%vec)

            sysvector%u(problem%degfdperm(sp1+1:sp1+ndof,2)) = &
                                           physvector%u(sp2+1:sp2+ndof)

            sp2 = sp2 + ndof

          end do
        end do
      end do

    else

!     store in nodes

      do nodenr = 1, mesh%nnodes

        sp1 = problem%nodnumdegfd(nodenr) + &
                        sum( problem%vec_nodnumdegfd(nodenr+1,&
                                    &problem%physq(1:physq-1)) - &
                             problem%vec_nodnumdegfd(nodenr,&
                                    &problem%physq(1:physq-1)) )

        sp2 = problem%vec_nodnumdegfd(nodenr,physvector%vec)

        ndof = problem%vec_nodnumdegfd(nodenr+1,physvector%vec) &
                 - problem%vec_nodnumdegfd(nodenr,physvector%vec)

        sysvector%u(problem%degfdperm(sp1+1:sp1+ndof,2)) = &
                                    physvector%u(sp2+1:sp2+ndof)


      end do

    end if

  end subroutine store_physvector


! extract a vector from another vector

  subroutine extract_vector ( mesh, problem, invector, outvector, indegfd, &
    outdegfd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the vector from which outvector must be extracted
    type(vector_t), intent(inout) :: invector

!   the vector to be extracted
    type(vector_t), intent(inout) :: outvector

!   degrees of freedom to be extracted from invector
    integer, dimension(:), intent(in) :: indegfd

!   degrees of freedom used for storing in outvector
!   size of outdegfd must be identical to indegfd
!   default = (/1,2,...,size(indegfd)/)
    integer, dimension(:), intent(in), optional :: outdegfd


!   NOTES:
!    - nodes in invector and/or outvector that have zero degrees of freedom
!      are simply ignored.
!    - all degrees of freedom in invector and outvector must exist, otherwise
!      an error message is produced


    integer :: sp1, sp2, nodenr, elgrp, elem, node, ndof1, ndof2, i
    integer, dimension(size(indegfd)) :: deg1, deg2


    call check ( mesh, 'extract_vector' )
    call check ( problem, 'extract_vector', mesh )

    if ( .not. invector%created ) then
      write(*,'(/a/)') &
        'Error in extract_vector: invector has not been created'
      stop
    end if
    if ( .not. outvector%created ) then
      write(*,'(/a/)') &
        'Error in extract_vector: outvector has not been created'
      stop
    end if
    if ( problem%probnr /= invector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in extract_vector: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in invector = ', invector%probnr
      stop
    end if
    if ( problem%probnr /= outvector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in extract_vector: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in outvector = ', outvector%probnr
      stop
    end if
    if ( invector%elementwise .and. .not. outvector%elementwise  ) then
      write(*,'(/3(a/))') &
        'Error in extract_vector: ',&
        ' it is not possible to extract a vector stored in nodes ', &
        ' from a vector stored elementwise.', &
        ' Use derive_vector instead.'
      stop
    end if

!   degrees of freedom

    deg1 = indegfd

    if ( any( indegfd <= 0 ) ) then
      write(*,'(/2(a/))') &
        'Error in extract_vector: ', &
        ' values in indegfd must be larger than zero '
      stop
    end if
    do elgrp = 1, problem%nelgrp
      if ( maxval(indegfd) > &
           minval(problem%vec_elnumdegfd(elgrp)%a(:,invector%vec), mask = &
                  problem%vec_elnumdegfd(elgrp)%a(:,invector%vec) > 0 ) ) then
         write(*,'(/4(a/))') &
           'Error in extract_vector: ', &
           ' indegfd has values larger than the number of degrees in ', &
           ' some nodes of invector (except possibly for nodes having zero ', &
           ' number of degrees, which are ignored).'
         stop
      end if
    end do

    if ( present(outdegfd) ) then
      if ( size(indegfd) /= size(outdegfd) ) then
        write(*,'(/2(a/))') &
          'Error in extract_vector: ', &
          ' size(outegfd) must be equal to size(indegfd) '
        stop
      end if
      if ( any( outdegfd <= 0 ) ) then
        write(*,'(/2(a/))') &
          'Error in extract_vector: ', &
          ' values in outdegfd must be larger than zero '
        stop
      end if
      do elgrp = 1, problem%nelgrp
        if ( maxval(outdegfd) > &
             minval(problem%vec_elnumdegfd(elgrp)%a(:,outvector%vec), mask = &
                 problem%vec_elnumdegfd(elgrp)%a(:,outvector%vec) > 0 ) ) then
           write(*,'(/4(a/))') &
             'Error in extract_vector: ', &
             ' outdegfd has values larger than the number of degrees in ', &
             ' some nodes of outvector (except possibly for nodes having ', &
             ' zero number of degrees, which are ignored).'
           stop
        end if
      end do
      deg2 = outdegfd
    else
      deg2 = [ (i,i=1,size(deg2)) ]
      do elgrp = 1, problem%nelgrp
        if ( maxval(deg2) > &
             minval(problem%vec_elnumdegfd(elgrp)%a(:,outvector%vec), mask = &
                 problem%vec_elnumdegfd(elgrp)%a(:,outvector%vec) > 0 ) ) then
           write(*,'(/5(a/))') &
             'Error in extract_vector: ', &
             ' unable to store the values in outvector: ', &
             ' degrees has values larger than the number of degrees in ', &
             ' some nodes of outvector (except possibly for nodes having ', &
             ' zero number of degrees, which are ignored).'
           stop
        end if
      end do
    end if

!   extract vector

    if ( invector%elementwise .and. outvector%elementwise ) then

!     read and store elementwise

      sp1 = 0
      sp2 = 0

      do elgrp = 1, problem%nelgrp
        do elem = 1, mesh%grpnumel(elgrp)
          do node = 1, mesh%elnumnod(elgrp)

            ndof1 = problem%vec_elnumdegfd(elgrp)%a(node,invector%vec)
            ndof2 = problem%vec_elnumdegfd(elgrp)%a(node,outvector%vec)

            if ( ndof1 > 0 .and. ndof2 > 0 ) then
              outvector%u(sp2+deg2) = invector%u(sp1+deg1)
            end if

            sp1 = sp1 + ndof1
            sp2 = sp2 + ndof2

          end do
        end do
      end do

    else if ( outvector%elementwise ) then

!     read nodes and store elementwise

      sp2 = 0

      do elgrp = 1, problem%nelgrp
        do elem = 1, mesh%grpnumel(elgrp)
          do node = 1, mesh%elnumnod(elgrp)

            nodenr = mesh%topology(elgrp)%a(node,elem)

            sp1 = problem%vec_nodnumdegfd(nodenr,invector%vec)

            ndof1 = problem%vec_nodnumdegfd(nodenr+1,invector%vec) &
                     - problem%vec_nodnumdegfd(nodenr,invector%vec)

            ndof2 = problem%vec_elnumdegfd(elgrp)%a(node,outvector%vec)

            if ( ndof1 > 0 .and. ndof2 > 0 ) then
              outvector%u(sp2+deg2) = invector%u(sp1+deg1)
            end if

            sp2 = sp2 + ndof2

          end do
        end do
      end do

    else

!     read and store in nodes

      do nodenr = 1, mesh%nnodes

        sp1 = problem%vec_nodnumdegfd(nodenr,invector%vec)

        ndof1 = problem%vec_nodnumdegfd(nodenr+1,invector%vec) &
                 - problem%vec_nodnumdegfd(nodenr,invector%vec)

        sp2 = problem%vec_nodnumdegfd(nodenr,outvector%vec)

        ndof2 = problem%vec_nodnumdegfd(nodenr+1,outvector%vec) &
                 - problem%vec_nodnumdegfd(nodenr,outvector%vec)

        if ( ndof1 > 0 .and. ndof2 > 0 ) then
          outvector%u(sp2+deg2) = invector%u(sp1+deg1)
        end if

      end do

    end if

  end subroutine extract_vector


! Transfer data (node for node) from a sysvector/vector to another
! sysvector/vector. The transfer is from 1 -> 2.

  subroutine transfer_data ( mesh, problem1, problem2, sysvector1, sysvector2, &
    vector1, vector2, degfd1, degfd2, physq1, physq2, layer1, layer2 )

    type(mesh_t), intent(in) :: mesh

!   the problem structure for the input vector
    type(problem_t), intent(in) :: problem1

!   if present the problem structure for the output vector, otherwise problem1
!   is used for the output vector as well
    type(problem_t), intent(in), optional :: problem2

!   the sysvectors
    type(sysvector_t), intent(in), optional :: sysvector1
    type(sysvector_t), intent(inout), optional :: sysvector2

!   the vectors
    type(vector_t), intent(in), optional :: vector1
    type(vector_t), intent(inout), optional :: vector2

!   degrees of freedom arrays
!   if both are present the size must be identical
    integer, dimension(:), intent(in), optional :: degfd1, degfd2

!   physical quantity arrays
    integer, dimension(:), intent(in), optional :: physq1, physq2

!   physical quantity arrays
    integer, intent(in), optional :: layer1, layer2


!   NOTES:
!    - The routine is rather forgiving: if a transfer is not possible, because
!      a degree in either the input or output is missing, the transfer is
!      simply ignored.
!    - In the output vector only copied degrees are modified. Data that is
!      undefined in the output vector (on input) remains undefined on output
!      if no data is copied to the degrees.


    integer :: ndegmax1, ndegmax2, deg, dg1, dg2, dg
    integer :: nodenr, ndof1, ndof2, ndeg, ndeg1, ndeg2, i
    integer, dimension(:), allocatable :: deg1, deg2, pos1, pos2, indx1, indx2


!   testing

    call check ( mesh, 'transfer_data' )
    call check ( problem1, 'transfer_data', mesh )
    if ( present(problem2) ) call check ( problem2, 'transfer_data', mesh )


!   input vectors testing

    if ( present(sysvector1) ) then

      if ( present(vector1) ) then
        write(*,'(/a/a/)') &
          'Error in transfer_data:', &
          ' sysvector1 and vector1 cannot be both present in heading'
        stop
      end if
      if ( .not. sysvector1%created ) then
        write(*,'(/a/)') &
          'Error in transfer_data: sysvector1 has not been created'
        stop
      end if
      if ( problem1%probnr /= sysvector1%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in transfer_data: ',&
          ' problem number in problem1 = ', problem1%probnr, &
          ' whereas problem number in sysvector1 = ', sysvector1%probnr
        stop
      end if

    else if ( present(vector1) ) then

      if ( .not. vector1%created ) then
        write(*,'(/a/)') &
          'Error in transfer_data: vector1 has not been created'
        stop
      end if
      if ( problem1%probnr /= vector1%probnr ) then
        write(*,'(/a/2(a,i0/))') &
          'Error in transfer_data: ',&
          ' problem number in problem1 = ', problem1%probnr, &
          ' whereas problem number in vector1 = ', vector1%probnr
        stop
      end if
      if ( vector1%elementwise ) then
        write(*,'(/2(a/))') &
          'Error in transfer_data: ',&
          ' vector1 is a vector stored elementwise. Not implemented.'
        stop
      end if
      if ( present(physq1) ) then
        write(*,'(/2(a/))') &
          'Warning in transfer_data: ',&
          ' physq1 not applicable for vector1, Ignored '
      end if

    else

      write(*,'(/2(a/))') &
        'Error in transfer_data: ',&
        ' neither sysvector1 or vector1 present '
      stop

    end if

!   output vectors testing

    if ( present(sysvector2) ) then

      if ( present(vector2) ) then
        write(*,'(/a/a/)') &
          'Error in transfer_data:', &
          ' sysvector2 and vector2 cannot be both present in heading'
        stop
      end if
      if ( .not. sysvector2%created ) then
        write(*,'(/a/)') &
          'Error in transfer_data: sysvector2 has not been created'
        stop
      end if
      if ( present(problem2) ) then
        if ( problem2%probnr /= sysvector2%probnr ) then
          write(*,'(/a/2(a,i0/))') &
            'Error in transfer_data: ',&
            ' problem number in problem2 = ', problem2%probnr, &
            ' whereas problem number in sysvector2 = ', sysvector2%probnr
          stop
        end if
      else
        if ( problem1%probnr /= sysvector2%probnr ) then
          write(*,'(/a/2(a,i0/))') &
            'Error in transfer_data: ',&
            ' problem number in problem1 = ', problem1%probnr, &
            ' whereas problem number in sysvector2 = ', sysvector2%probnr
          stop
        end if
      end if

    else if ( present(vector2) ) then

      if ( .not. vector2%created ) then
        write(*,'(/a/)') &
          'Error in transfer_data: vector2 has not been created'
        stop
      end if
      if ( vector2%elementwise ) then
        write(*,'(/2(a/))') &
          'Error in transfer_data: ',&
          ' vector2 is a vector stored elementwise. Not implemented.'
        stop
      end if
      if ( present(problem2) ) then
        if ( problem2%probnr /= vector2%probnr ) then
          write(*,'(/a/2(a,i0/))') &
            'Error in transfer_data: ',&
            ' problem number in problem2 = ', problem2%probnr, &
            ' whereas problem number in vector2 = ', vector2%probnr
          stop
        end if
      else
        if ( problem1%probnr /= vector2%probnr ) then
          write(*,'(/a/2(a,i0/))') &
            'Error in transfer_data: ',&
            ' problem number in problem1 = ', problem1%probnr, &
            ' whereas problem number in vector2 = ', vector2%probnr
          stop
        end if
      end if
      if ( present(physq2) ) then
        write(*,'(/2(a/))') &
          'Warning in transfer_data: ',&
          ' physq2 not applicable for vector2, Ignored '
      end if

    else

      write(*,'(/2(a/))') &
        'Error in transfer_data: ',&
        ' neither sysvector2 or vector2 present '
      stop

    end if

!   degrees of freedom

    if ( present(degfd1) ) then
      if ( any( degfd1 <= 0 ) ) then
        write(*,'(/2(a/))') &
          'Error in transfer_data: ', &
          ' values in degfd1 must be larger than zero '
        stop
      end if
    end if

    if ( present(degfd2) ) then
      if ( any( degfd2 <= 0 ) ) then
        write(*,'(/2(a/))') &
          'Error in transfer_data: ', &
          ' values in degfd2 must be larger than zero '
        stop
      end if
    end if

    if ( present(degfd1) .and. present(degfd2) ) then
      if ( size(degfd1) /= size(degfd2) ) then
        write(*,'(/3(a/))') &
          'Warning in transfer_data: ', &
          ' size of degfd1 and degfd2 different ', &
          ' ignoring unmatched degrees '
      end if
    end if

!   physical quantities

    if ( present(physq1) ) then

      if ( problem1%nphysq == 0 ) then
        write(*,'(2(/a)/)') &
          'Error in transfer_data: ',&
          ' no physical quantities defined in problem1 '
        stop
      end if

      if ( any( physq1 < 1 ) .or. any( physq1 > problem1%nphysq ) ) then
        write(*,'(2(/a)/)') &
          'Error in transfer_data: ', &
          ' physical quantities in physq1 are out of range.'
        stop
      end if

    end if

    if ( present(physq2) ) then

      if ( present(problem2) ) then

        if ( problem2%nphysq == 0 ) then
          write(*,'(2(/a)/)') &
            'Error in transfer_data: ',&
            ' no physical quantities defined in problem2 '
          stop
        end if

        if ( any( physq2 < 1 ) .or. any( physq2 > problem2%nphysq ) ) then
          write(*,'(2(/a)/)') &
            'Error in transfer_data: ', &
            ' physical quantities in physq2 are out of range.'
          stop
        end if

      else

        if ( problem1%nphysq == 0 ) then
          write(*,'(2(/a)/)') &
            'Error in transfer_data: ',&
            ' no physical quantities defined in problem1 '
          stop
        end if

        if ( any( physq2 < 1 ) .or. any( physq2 > problem1%nphysq ) ) then
          write(*,'(2(/a)/)') &
            'Error in transfer_data: ', &
            ' physical quantities in physq2 are out of range.'
          stop
        end if

      end if

    end if

!   layers

    if ( present(layer1) ) then
      if ( problem1%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in transfer_data: ',&
          ' specification of layer1 not allowed because', &
          ' layers have not been defined.'
        stop
      end if
      if ( layer1 < 1 .or. layer1 > problem1%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in build_system: layer1 out of range.'
        stop
      end if
    end if

    if ( present(layer2) ) then
      if ( present(problem2) ) then
        if ( problem2%numlayers == 0 ) then
          write(*,'(3(/a)/)') &
            'Error in transfer_data: ',&
            ' specification of layer2 not allowed because', &
            ' layers have not been defined.'
          stop
        end if
        if ( layer2 < 1 .or. layer2 > problem2%numlayers ) then
          write(*,'(2(/a)/)') &
            'Error in build_system: layer2 out of range.'
          stop
        end if
      else
        if ( problem1%numlayers == 0 ) then
          write(*,'(3(/a)/)') &
            'Error in transfer_data: ',&
            ' specification of layer2 ', &
            ' is only available if layers have been defined.'
          stop
        end if
        if ( layer2 < 1 .or. layer2 > problem1%numlayers ) then
          write(*,'(2(/a)/)') &
            'Error in build_system: layer2 out of range.'
          stop
        end if
      end if
    end if


!   reserve some memory

    ndegmax1 = max ( problem1%maxnoddegfd, problem1%maxvecnoddegfd )
    if ( present(problem2) ) then
      ndegmax2 = max ( problem2%maxnoddegfd, problem2%maxvecnoddegfd )
    else
      ndegmax2 = ndegmax1
    end if
    allocate ( deg1(ndegmax1), deg2(ndegmax2) )
    allocate ( pos1(ndegmax1), pos2(ndegmax2) )

    ndegmax1 = max ( problem1%numdegfd, maxval(problem1%vec_numdegfd(:,1)) )
    if ( present(problem2) ) then
      ndegmax2 = max ( problem2%numdegfd, maxval(problem2%vec_numdegfd(:,1)) )
    else
      ndegmax2 = ndegmax1
    end if
    allocate ( indx1(ndegmax1), indx2(ndegmax2) )

!   loop over nodes to build the index arrays

    deg = 0

    do nodenr = 1, mesh%nnodes

!     find position arrays of data in the nodes

!     input vector

      if ( present(sysvector1) ) then
        call pos_array_node ( problem1, nodenr, ndof1, pos1, physq1, layer1 )
      else if ( present(vector1) ) then
        call pos_array_vec_node ( problem1, nodenr, ndof1, pos1, vector1%vec, &
          layer1 )
      end if

!     output vector

      if ( present(sysvector2) ) then
        if ( present(problem2) ) then
          call pos_array_node ( problem2, nodenr, ndof2, pos2, physq2, layer2 )
        else
          call pos_array_node ( problem1, nodenr, ndof2, pos2, physq2, layer2 )
        end if
      else if ( present(vector2) ) then
        if ( present(problem2) ) then
          call pos_array_vec_node ( problem2, nodenr, ndof2, pos2, &
            vector2%vec, layer2 )
        else
          call pos_array_vec_node ( problem1, nodenr, ndof2, pos2, &
            vector2%vec, layer2 )
        end if
      end if

!     degrees of freedom

      if ( present(degfd1) ) then
        ndeg1 = size(degfd1)
        deg1(1:ndeg1) = degfd1
      else
        ndeg1 = ndof1
        deg1(1:ndeg1) = [(i,i=1,ndof1)]
      end if

      if ( present(degfd2) ) then
        ndeg2 = size(degfd2)
        deg2(1:ndeg2) = degfd2
      else
        ndeg2 = ndof2
        deg2(1:ndeg2) = [(i,i=1,ndof2)]
      end if

      ndeg = min(ndeg1,ndeg2)

!     loop over degrees of freedom: fill index arrays

      do dg = 1, ndeg
        dg1 = deg1(dg)
        dg2 = deg2(dg)
        if ( dg1 > ndof1 .or. dg2 > ndof2 ) cycle
        deg = deg + 1
        indx1(deg) = pos1(dg1)
        indx2(deg) = pos2(dg2)
      end do

    end do

!   transfer actual data

    if ( present(sysvector1) ) then

      if ( present(sysvector2) ) then

        sysvector2%u(indx2(1:deg)) = sysvector1%u(indx1(1:deg))

      else if ( present(vector2) ) then

        vector2%u(indx2(1:deg)) = sysvector1%u(indx1(1:deg))

      end if

    else if ( present(vector1) ) then

      if ( present(sysvector2) ) then

        sysvector2%u(indx2(1:deg)) = vector1%u(indx1(1:deg))

      else if ( present(vector2) ) then

        vector2%u(indx2(1:deg)) = vector1%u(indx1(1:deg))

      end if

    end if

    deallocate ( deg1, deg2, pos1, pos2, indx1, indx2 )

  end subroutine transfer_data


! Fill vector with values

  subroutine fill_vector ( mesh, problem, vector, layer, degfd, degsfd, &
    point, points, curve1, curve2, step, exclude, curves, node1, node2, nodes, &
    surface1, surface2, surfaces, elgroup1, elgroup2, elgroups, nodeset1, &
    nodeset2, nodesets, excludepoints, excludecurves, &
    excludesurfaces, value, func, funcnr, vvalue, vfunc, vfuncnr, &
    ignore_nodegfd, blend  )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(vector_t), intent(inout) :: vector

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered.
    integer, intent(in), optional :: layer

!   if degfd is present
!     the degree of freedom given by degfd is given a value
!   else
!     all degrees are given a value
!   NOTE: degfd cannot be combined with degsfd
    integer, intent(in), optional :: degfd

!   else if degsfd is present
!     the degrees of freedom given by degsfd are given a value
!   else
!     all degrees are given a value
!   NOTE: degsfd cannot be combined with degfd
    integer, intent(in), dimension(:), optional :: degsfd

!   if present and .true.: ignore nodes that do not have a degree of freedom
!   given by degfd or degsfd.
!   default = .false.
    logical, intent(in), optional :: ignore_nodegfd

!   only the nodes of mesh blend are considered to fill the vector
!     blend = -1 all nodes are considered
!     blend = 0  only nodes of the main mesh with element shapes defined by
!                mesh%element are considered
!     blend > 0  only nodes of the blend mesh with element shapes defined by
!                mesh%element_blend are considered.
!   default=0
    integer, intent(in), optional :: blend

!   if point is present: set value in a point
    integer, intent(in), optional :: point

!   if points is present: set value in the point given
!   NOTE: if point is present, points is not used.
    integer, intent(in), dimension(:), optional :: points

!   if curve1 is present
!     if curve2 is present
!       set value along curves curve1,...,curve2
!     otherwise
!       set value along curve curve1
!   step and exclude can be used together with curve1.
!   See the NOTE with step, exclude below for the restrictions.
    integer, intent(in), optional :: curve1, curve2

!   if curves is present set value along the curves given
!   step and exclude can be used together with curves.
!   See the NOTE with step, exclude below for the restrictions.
!   NOTE: if curve1 is present, curves is not used.
    integer, intent(in), dimension(:), optional :: curves

!   if node1 is present
!     if node2 is present
!       set value along for node1,...,node2
!     otherwise
!       set value for node1
!   step can be used together with node1 and node2
    integer, intent(in), optional :: node1, node2

!   if nodes is present set value for the nodes given.
!   step can be used together with nodes.
!   NOTE: if node1 is present, nodes is not used.
    integer, intent(in), dimension(:), optional :: nodes

!   step, exclude:
!     if curve1 is present
!       step is given by
!         step=0 all nodes (default)
!         step>0 nodes 1, 1+step, 1+2*step, ...
!         step<0 except nodes 1, 1-step, 1-2*step, ...
!       note that this is applied for all curves _separately_ and
!       curve1 to curve2 is not treated as a single curve.
!       if exclude is present it is given by
!         exclude=0 no excludes, use all nodes
!         exclude=1 exclude first point of first curve
!         exclude=2 exclude last point of last curve
!         exclude=3 exclude first and last point
!     else if node1 is present
!       if step present, the nodes node1, node+step, node1+2*step, ... are set
!
!   NOTE: step, exclude only work for curves where the local node numbering is
!   in a natural sequence along the curve. This might not be the case for
!   curves constructed from several other curves or curves read from external
!   mesh generators.
    integer, intent(in), optional :: step, exclude

!   if surface1 is present
!     if surface2 is present
!       set value along surfaces surface1,...,surface2
!     otherwise
!       set value along surface surface1
    integer, intent(in), optional :: surface1, surface2

!   if surfaces is present values are set along the surfaces given.
!   NOTE: if surface1 is present, surfaces is not used.
    integer, intent(in), dimension(:), optional :: surfaces

!   if elgroup1 is present
!     if elgroup2 is present
!       set value along elgroups elgroup1,...,elgroup2
!     otherwise
!       set value along elgroup elgroup1
    integer, intent(in), optional :: elgroup1, elgroup2

!   if elgroups is present values are set in the element groups given.
!   NOTE: if elgroup1 is present, elgroups is not used.
    integer, intent(in), dimension(:), optional :: elgroups

!   if nodeset1 is present
!     if nodeset2 is present
!       set value along nodesets nodeset1,...,nodeset2
!     otherwise
!       set value along nodeset nodeset1
    integer, intent(in), optional :: nodeset1, nodeset2

!   if nodesets is present values are set along the nodesets given.
!   NOTE: if nodeset1 is present, nodesets is not used.
    integer, intent(in), dimension(:), optional :: nodesets

!   exclude points:
!   for example excludepoints=(/1,3/) excludes nodes in points P1 and P3
!   to be set. Can be used together with the parameters curves1, curve2,
!   surface1, surface2, nodeset1, nodeset2, elgroup1, elgroup2.
!   Note that only the nodes that are actually on the curves, surfaces,
!   nodesets or groups are excluded.
    integer, intent(in), dimension(:), optional :: excludepoints

!   exclude curves:
!   for example excludecurves=(/1,3/) excludes nodes on curves C1 and C3
!   to be set. Can be used together with the parameters curves1, curve2,
!   surface1, surface2, nodeset1, nodeset2, elgroup1, elgroup2.
!   Note that only the nodes that are actually on the curves, surfaces,
!   nodesets or groups are excluded.
    integer, intent(in), dimension(:), optional :: excludecurves

!   exclude surfaces:
!   for example excludesurfaces=(/1,3/) excludes nodes on surfaces S1 and S3
!   to be set. Can be used together with the parameters curve1, curve2,
!   surface1, surface2, nodeset1, nodeset2, elgroup1 elgroup2.
!   Note that only the nodes that are actually on the curves, surfaces, nodesets
!   or groups are excluded.
    integer, intent(in), dimension(:), optional :: excludesurfaces

!   if value is present
!     value is the assigned value
!   else if func present
!     function routine func is used to assign a value. funcnr is assigned to nr
!     and the nodal coordinates are assigned to x.
!   else if vvalue is present
!     vvalue is the assigned vector value.
!     NOTE: the size of the vector needs to be at least equal to the maximum
!     number of degrees in a node that needs assignment.
!   else if vfunc present
!     function routine vfunc is used to assign a value. The size is assigned to
!     n, vfuncnr is assigned to nr and the nodal coordinates are assigned to x.
!     NOTE, BE CAREFUL: in general the size of the vector depends on the
!     number of degrees in a node that needs to be assigned a value. The number
!     of degrees is given by n and might vary from node to node.
    real(dp), intent(in), optional :: value
    optional :: func
    interface
      function func ( nr, x )
        use kind_defs_m
        implicit none
        integer, intent(in) :: nr
        real(dp), intent(in), dimension(:) :: x
        real(dp) :: func
      end function func
    end interface
    integer, intent(in), optional :: funcnr
    real(dp), dimension(:), intent(in), optional :: vvalue
    optional :: vfunc
    interface
      function vfunc ( n, nr, x )
        use kind_defs_m
        implicit none
        integer, intent(in) :: n, nr
        real(dp), intent(in), dimension(:) :: x
        real(dp), dimension(n) :: vfunc
      end function vfunc
    end interface
    integer, intent(in), optional :: vfuncnr

!   This routine fills a vector with values such as initial values. It is a
!   clone of the routine fill_sysvector. Most options are less useful here
!   though, but are retained anyway. A vector can also be filled using
!   derive_vector+an element subroutine. For simple filling, for example with
!   node1=1, node2=mesh%nnodes, this routine is easier. Note that this routine
!   can only be used for vectors with vector%elementwise=.false. (nodal
!   storage)
!

    logical :: excludenodes, lignore_nodegfd
    integer :: worklength, curve, nn, sz, stp, node, nodenr
    integer :: nodeset, nval
    integer :: surface, numec, numes, enode, elem
    integer :: ndegfd, elgrp
    integer, allocatable, dimension(:) :: wnodes, work, ldegsfd, &
      lcurves, lsurfaces, lnodes, lelgroups, lnodesets
    integer :: pos(problem%maxvecnoddegfd), i, lblend

    lignore_nodegfd  = &
                set_optional ( variable=ignore_nodegfd, default = .false. )

    if ( present(blend) ) then
      if ( blend < -1 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error fill_vector: invalid argument ', &
          '  blend < -1 or larger than ', mesh%nblend
        stop
      end if
      lblend = blend
    else
      lblend = 0 ! default
    end if

!   testing

    call check ( mesh, 'fill_vector' )
    call check ( problem, 'fill_vector', mesh )

    if ( .not. vector%created ) then
      write(*,'(/a/)') &
        'Error in fill_vector: vector not created.'
      stop
    end if

    if ( problem%probnr /= vector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in fill_vector: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in vector = ', vector%probnr
      stop
    end if

    if ( vector%elementwise ) then
      write(*,'(/3(a/))') &
        'Error in fill_vector: filling a vector with elementwise stored ', &
        ' data is not implemented. Use the routine derive_vector with an ', &
        ' element subroutine instead. '
      stop
    end if

    if ( present(layer) ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in fill_vector: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer <=0 .or. layer > problem%numlayers ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_vector: ',&
          ' layer <=0 or layer > number of layers = ', &
          problem%numlayers
        stop
      end if
    end if

    if ( present(excludepoints) ) then
      if ( any(excludepoints <=0) .or. any(excludepoints > mesh%npoints) ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_vector: ',&
          ' excludepoints <=0 or excludepoints > number of points = ', &
           mesh%npoints
        stop
      end if
    end if

    if ( present(excludecurves) ) then
      if ( any(excludecurves <=0) .or. any(excludecurves > mesh%ncurves) ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_vector: ',&
          ' excludecurves <=0 or excludecurves > number of curves = ', &
           mesh%ncurves
        stop
      end if
    end if

    if ( present(excludesurfaces) ) then
      if ( any(excludesurfaces <=0) .or. &
           any(excludesurfaces > mesh%nsurfaces) ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_vector: ',&
          ' excludesurfaces <=0 or excludesurfaces > number of surfaces = ', &
           mesh%nsurfaces
        stop
      end if
    end if

    if ( present(degfd) .and. present(degsfd) ) then
      write(*,'(2(/a)/)') &
        'Error in fill_vector: ',&
        '  Both degfd and degsfd present in heading.'
      stop
    end if

!   set work for excludes

    if ( present(excludepoints) .or. present(excludecurves) .or. &
         present(excludesurfaces) ) then
      excludenodes = .true.
!     create workspace for nodes on points, curves or surfaces
      allocate ( work(mesh%nnodes) )
      work = 0
      if ( present(excludepoints) ) then
!       set work
        work ( mesh%points(excludepoints) ) = 1
      end if
      if ( present(excludecurves) ) then
        numec = size ( excludecurves )
!       set work
        do curve = 1, numec
          work ( mesh%curves(excludecurves(curve))%nodes ) = 1
        end do
      end if
      if ( present(excludesurfaces) ) then
        numes = size ( excludesurfaces )
!       set work
        do surface = 1, numes
          work ( mesh%surfaces(excludesurfaces(surface))%nodes ) = 1
        end do
      end if
    else
      excludenodes = .false.
    end if


!   determine type of call

    if ( present(point) ) then

!     fill for a point

      if ( point <=0 .or. point > mesh%npoints ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_vector: ',&
          ' point <=0 or point > number of points = ', &
           mesh%npoints
        stop
      end if

      wnodes = [mesh%points(point)]

      nn = 1

    else if ( present(points) ) then

      if ( any( points <=0 .or. points > mesh%npoints ) ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_vector: ',&
          ' points <=0 or points > number of points = ', &
           mesh%npoints
        stop
      end if

      wnodes = [mesh%points(points)]

      nn = size(points)

    else if ( present(curve1) .or. present(curves) ) then

!     fill for curves

      if ( present(curve1) ) then

        if ( curve1 <=0 .or. curve1 > mesh%ncurves ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_vector: ',&
            ' curve1 <=0 or curve1 > number of curves = ', &
             mesh%ncurves
          stop
        end if

        if ( present(curve2) ) then

          if ( curve2 <=0 .or. curve2 > mesh%ncurves ) then
            write(*,'(2(/a),i0/)') &
              'Error in fill_vector: ',&
              ' curve2 <=0 or curve2 > number of curves = ', &
               mesh%ncurves
            stop
          end if

          lcurves = [ (i, i=curve1,curve2) ]

        else

          lcurves = [ curve1 ]

        end if

      else if ( present(curves) ) then

        if ( any ( curves <=0 .or. curves > mesh%ncurves ) ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_vector: ',&
            ' curves <=0 or curves > number of curves = ', &
             mesh%ncurves
          stop
        end if

        lcurves = curves

      end if

!     determine worklength for wnodes

      worklength = 0
      do i = 1, size(lcurves)
        worklength = worklength + size(mesh%curves(lcurves(i))%nodes)
      end do

      allocate ( wnodes(worklength) )

!     determine nodes for filling with a value

      nn = 0

      do i = 1, size(lcurves)

        curve = lcurves(i)
        sz = size(mesh%curves(curve)%nodes)

        if ( present(step) ) then
          stp = step
        else
          stp = 0
        end if

        if ( stp == 0 ) then

          if ( excludenodes ) then

!           exclude nodes

            do node = 1, sz
              nodenr = mesh%curves(curve)%nodes(node)
              if ( work(nodenr) == 0 ) then
!               include node
                wnodes(nn+node) = nodenr
              else
!               exclude node
                wnodes(nn+node) = 0
              end if
            end do

          else

!           all nodes
            wnodes(nn+1:nn+sz) = mesh%curves(curve)%nodes

         end if

        else if ( stp > 0 ) then

!         nodes according to step
          wnodes(nn+1:nn+sz) = 0

          if ( excludenodes ) then

!           exclude nodes

            do node = 1, sz, stp
              nodenr = mesh%curves(curve)%nodes(node)
              if ( work(nodenr) == 0 ) then
!               include node
                wnodes(nn+node) = nodenr
              end if
            end do

          else

            wnodes(nn+1:nn+sz:stp) = mesh%curves(curve)%nodes(1:sz:stp)

          end if

        else if ( stp < 0 ) then

!         exclude nodes according to step
          wnodes(nn+1:nn+sz) = mesh%curves(curve)%nodes

          if ( excludenodes ) then

!           exclude nodes

            do node = 1, sz, -stp
              nodenr = mesh%curves(curve)%nodes(node)
              if ( work(nodenr) /= 0 ) then
!               exclude node
                wnodes(nn+node) = 0
              end if
            end do

          else

            wnodes(nn+1:nn+sz:-stp) = 0

          end if

        end if

        nn = nn + sz

      end do

      if ( present(exclude) ) then

        if ( exclude < 0 .or. exclude > 4 ) then

          write(*,'(2(/a)/)') &
            'Error: exclude in the heading of fill_vector can only ', &
            'have the values 0, 1, 2 or 3.'
          stop

        end if

        if ( exclude == 1 .or. exclude == 3 ) then

!         exclude first point
          wnodes(1) = 0

        end if

        if ( exclude == 2 .or. exclude == 3 ) then

!         exclude last point
          wnodes(nn) = 0

        end if

      end if

    else if ( present(surface1) .or. present(surfaces) ) then

!     fill for surfaces

      if ( present(surface1) ) then

        if ( surface1 <=0 .or. surface1 > mesh%nsurfaces ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_vector: ',&
            ' surface1 <=0 or surface1 > number of surfaces = ', &
             mesh%nsurfaces
          stop
        end if

        if ( present(surface2) ) then

          if ( surface2 <=0 .or. surface2 > mesh%nsurfaces ) then
            write(*,'(2(/a),i0/)') &
              'Error in fill_vector: ',&
              ' surface2 <=0 or surface2 > number of surfaces = ', &
               mesh%nsurfaces
            stop
          end if

          lsurfaces = [ (i, i=surface1,surface2) ]

        else

          lsurfaces = [ surface1 ]

        end if

      else if ( present(surfaces) ) then

        if ( any ( surfaces <=0 .or. surfaces > mesh%nsurfaces ) ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_vector: ',&
            ' surfaces <=0 or surfaces > number of surfaces = ', &
             mesh%nsurfaces
          stop
        end if

        lsurfaces = surfaces

      end if

!     determine worklength for wnodes

      worklength = 0
      do i = 1, size(lsurfaces)
        worklength = worklength + size(mesh%surfaces(lsurfaces(i))%nodes)
      end do

      allocate ( wnodes(worklength) )

!     determine nodes for filling with a value

      nn = 0

      do i = 1, size(lsurfaces)

        surface = lsurfaces(i)
        sz = size(mesh%surfaces(surface)%nodes)

        if ( excludenodes ) then

!         exclude nodes on curves

          do node = 1, sz
            nodenr = mesh%surfaces(surface)%nodes(node)
            if ( work(nodenr) == 0 ) then
!             include node
              wnodes(nn+node) = nodenr
            else
!             exclude node
              wnodes(nn+node) = 0
            end if
          end do

        else

!         all nodes
          wnodes(nn+1:nn+sz) = mesh%surfaces(surface)%nodes

        end if

        nn = nn + sz

      end do

    else if ( present(node1) .or. present(nodes) ) then

!     fill for nodes

      if ( present(node1) ) then

        if ( node1 <=0 .or. node1 > mesh%nnodes ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_vector: ',&
            ' node1 <=0 or node1 > number of nodes = ', &
             mesh%nnodes
          stop
        end if

        if ( present(node2) ) then
          if ( node2 <=0 .or. node2 > mesh%nnodes ) then
            write(*,'(2(/a),i0/)') &
              'Error in fill_vector: ',&
              ' node2 <=0 or node1 > number of nodes = ', &
               mesh%nnodes
            stop
          end if
          lnodes = [(i,i=node1,node2)]
        else
          lnodes = [node1]
        end if

      else if ( present(nodes) ) then

        if ( any( nodes <=0 .or. nodes > mesh%nnodes ) ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_vector: ',&
            ' nodes <=0 or nodes > number of nodes = ', &
             mesh%nnodes
          stop
        end if

        lnodes = nodes

      end if

!     determine worklength for wnodes

      worklength = size(lnodes)

      allocate ( wnodes(worklength) )

!     determine nodes for filling with a value

      nn = worklength

      if ( present(step) ) then
!       nodes according to step
        wnodes(1:nn) = 0
        wnodes(1:nn:step) = lnodes(1:nn:step)
      else
!       all nodes (or step=1)
        wnodes(1:nn) = lnodes(1:nn)
      end if

    else if ( present(nodeset1) .or. present(nodesets) ) then

!     fill for nodesets

      if ( present(nodeset1) ) then

        if ( nodeset1 <=0 .or. nodeset1 > mesh%nnodesets ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_sysvector: ',&
            ' nodeset1 <=0 or nodeset1 > number of nodesets = ', &
             mesh%nnodesets
          stop
        end if

        if ( present(nodeset2) ) then

          if ( nodeset2 <=0 .or. nodeset2 > mesh%nnodesets ) then
            write(*,'(2(/a),i0/)') &
              'Error in fill_sysvector: ',&
              ' nodeset2 <=0 or nodeset2 > number of nodesets = ', &
               mesh%nnodesets
            stop
          end if

          lnodesets = [ (i,i=nodeset1,nodeset2) ]

        else

          lnodesets = [ nodeset1 ]

        end if

      else if ( present(nodesets) ) then

        if ( any( nodesets <=0 .or. nodesets > mesh%nnodesets ) ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_sysvector: ',&
            ' nodesets <=0 or nodesets > number of nodesets = ', &
             mesh%nnodesets
          stop
        end if

        lnodesets = nodesets

      end if

!     determine worklength for wnodes

      worklength = 0
      do i = 1, size(lnodesets)
        worklength = worklength + size(mesh%nodesets(lnodesets(i))%a)
      end do

      allocate ( wnodes(worklength) )

!     determine nodes for filling with a value

      nn = 0

      do i = 1, size(lnodesets)

        nodeset = lnodesets(i)
        sz = size(mesh%nodesets(nodeset)%a)

        if ( excludenodes ) then

!         exclude nodes

          do node = 1, sz
            nodenr = mesh%nodesets(nodeset)%a(node)
            if ( work(nodenr) == 0 ) then
!             include node
              wnodes(nn+node) = nodenr
            else
!             exclude node
              wnodes(nn+node) = 0
            end if
          end do

        else

!         all nodes
          wnodes(nn+1:nn+sz) = mesh%nodesets(nodeset)%a

        end if

        nn = nn + sz

      end do

   else if ( present(elgroup1) .or. present(elgroups) ) then

!     fill for groups

      if ( present(elgroup1) ) then

        if ( elgroup1 <=0 .or. elgroup1 > mesh%nelgrp ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_vector: ',&
            ' elgroup1 <=0 or elgroup1 > number of groups = ', &
             mesh%nelgrp
          stop
        end if

        if ( present(elgroup2) ) then

          if ( elgroup2 <=0 .or. elgroup2 > mesh%nelgrp ) then
            write(*,'(2(/a),i0/)') &
              'Error in fill_vector: ',&
              ' elgroup2 <=0 or elgroup2 > number of groups = ', &
               mesh%nelgrp
            stop
          end if

          lelgroups = [ (i,i=elgroup1,elgroup2) ]

        else

          lelgroups = [ elgroup1 ]

        end if

      else if ( present(elgroups) ) then

        if ( any( elgroups <=0 .or. elgroups > mesh%nelgrp ) ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_vector: ',&
            ' elgroups <=0 or elgroups > number of elgroups = ', &
             mesh%nelgrp
          stop
        end if

        lelgroups = elgroups

      end if

!     determine worklength for wnodes

      worklength = 0
      do i = 1, size(lelgroups)
        elgrp = lelgroups(i)
        worklength = worklength + mesh%grpnumel(elgrp) * mesh%elnumnod(elgrp)
      end do

      allocate ( wnodes(worklength) )

!     determine nodes for filling with a value

      nn = 0

      do i = 1, size(lelgroups)

        elgrp = lelgroups(i)
        sz = mesh%grpnumel(elgrp)*mesh%elnumnod(elgrp)

        if ( excludenodes ) then

!         exclude nodes

          node = 0
          do elem = 1, mesh%grpnumel(elgrp)
            do enode = 1, mesh%elnumnod(elgrp)
              node = node + 1
              nodenr = mesh%topology(elgrp)%a(enode,elem)
              if ( work(nodenr) == 0 ) then
!               include node
                wnodes(nn+node) = nodenr
              else
!               exclude node
                wnodes(nn+node) = 0
              end if
            end do
          end do

        else

!         all nodes
          wnodes(nn+1:nn+sz) = reshape ( mesh%topology(elgrp)%a, [sz] )

        end if

        nn = nn + sz

      end do

    else

      write(*,'(20(/a)/)') &
        'Error: could not determine type of fill from the heading of ', &
        'fill_vector. Non of the keywords:',      &
        '  point or points ',     &
        '  curve1 or curves ',    &
        '  surface1 or surfaces ',  &
        '  node1 or nodes  ',    &
        '  elgroup1 or elgroups ',  &
        '  nodeset1 or nodesets ',  &
        'is present.'
      stop

    end if

!   remove work for excludes

    if ( present(excludepoints) .or. present(excludecurves) .or. &
         present(excludesurfaces) ) then
      deallocate ( work )
    end if

!   Fill degree of freedom in nodes

    do node = 1, nn

      nodenr = wnodes(node)

      if ( node < nn ) then
!        node present more than once?
        if ( wnodes(node+1) == nodenr ) cycle
      end if

      if ( lblend >= 0 ) then
!       only consider nodes of the blend mesh lblend
        if ( nodenr <= mesh%nnodes_blend(lblend+1) .or. &
                  nodenr > mesh%nnodes_blend(lblend+2) ) cycle
      end if

      if ( nodenr > 0 ) then

!       node found

        if ( present(layer) ) then
!         cycle if layer absent in node
          if ( .not. btest(problem%nodlayers(nodenr),layer-1) ) cycle
        end if

!       get positions of degrees

        call pos_array_vec_node ( problem, nodenr, ndegfd, pos, vector%vec, &
          layer )

        if ( present(degfd) ) then

!         single degree

          if ( degfd < 1 .or. degfd > ndegfd ) then

            if ( lignore_nodegfd ) cycle

            write(*,'(/a/2(a,i0),/a,i0/)') &
              'Error: degfd in the heading of fill_vector is outside', &
              'range in node ', nodenr, '; degfd is ', degfd, &
              'whereas the number of degrees of freedom is ', ndegfd
            stop

          end if

          ldegsfd = [degfd]

        else if ( present(degsfd) ) then

!         multiple degrees

          if ( any ( degsfd < 1 .or. degsfd > ndegfd ) ) then

            if ( lignore_nodegfd ) cycle

            write(*,'(/a/2(a,i0),/a,20i0/)') &
              'Error: degsfd in the heading of fill_vector is outside', &
              'range in node ', nodenr, &
              'The number of degrees of freedom is ', ndegfd, &
              'whereas degsfd is ', degsfd
            stop

          end if

          ldegsfd = degsfd

        else

!         default: all degrees

          ldegsfd = [(i,i=1,ndegfd)]

        end if

        if ( present(value) ) then

!         value given; takes precedence over function
          vector%u(pos(ldegsfd)) = value

        else if ( present(func) ) then

!         function given

          if ( present(funcnr) ) then

            vector%u(pos(ldegsfd)) = func ( funcnr, mesh%coor(nodenr,:) )

          else

            write(*,'(2(/a)/)') &
              'Error: if func is present in the heading of ', &
              'fill_vector, funcnr must be present as well.'
            stop

          end if

        else if ( present(vvalue) ) then

!         vvalue given; takes precedence over vector function

          nval = size(ldegsfd)

          if ( nval <= size(vvalue) ) then

            vector%u(pos(ldegsfd)) = vvalue(1:nval)

          else

            write(*,'(/a/2(a,i0)/)') &
              'Error in fill_vector: size of vvalue is too small:', &
              '  size = ', size(vvalue), ', required = ', nval
            stop

          end if

        else if ( present(vfunc) ) then

!         vector function given

          if ( present(vfuncnr) ) then

            vector%u(pos(ldegsfd)) = &
                        vfunc ( size(ldegsfd), vfuncnr, mesh%coor(nodenr,:) )

          else

            write(*,'(2(/a)/)') &
              'Error: if vfunc is present in the heading of ', &
              'fill_vector, vfuncnr must be present as well.'
            stop

          end if

        else

          write(*,'(3(/a)/)') &
            'Error: either value, func, vvalue or vfunc must be present in ', &
            ' the heading of fill_vector in order to fill the vector ', &
            ' with values.'
          stop

        end if

      end if

    end do

    deallocate ( wnodes )

  end subroutine fill_vector


! create vector subscript array in vector

  subroutine create_subscript_vector ( mesh, problem, subscript, degfd, &
    degsfd, layer, vec, physq, points, curves, surfaces, nodesets, &
    elementsets, groups, excludepoints, excludecurves, excludesurfaces, &
    excludenodesets, excludeelementsets, xmin, xmax, elementwise, fillnodes, &
    includenodesonce )

    type(mesh_t), intent(in)  :: mesh
    type(problem_t), intent(in)  :: problem

!   the vector subscript of the degrees of freedom in a vector
    type(subscriptvec_t), intent(inout) :: subscript

!   a vector subscript is made for vector nr vec
!   note: if vec is present, the physq cannot be present
    integer, intent(in), optional :: vec

!   if physq is present a vector subscript is made for the vector defining the
!   physical quantity physq (=problem%phys(physq)).
!   note: if physq is present, the vec cannot be present
    integer, intent(in), optional :: physq

!   if degfd is present only the degfd^th degree of freedom will be included,
!   otherwise all degrees will be included.
!   Note, degfd and degsfd cannot be present at the same time.
    integer, intent(in), optional :: degfd

!   if degsfd is present only the degrees of freedom in degsfd will be included,
!   otherwise all degrees will be included.
!   Note, degfd and degsfd cannot be present at the same time.
    integer, intent(in), dimension(:), optional :: degsfd

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered.
    integer, intent(in), optional :: layer

!   element groups included for defining the subscript.
!   For example groups=(/2,4/), includes nodes in element groups 2 and 4.
!   Default is to include all element groups.
!   Note that nodes are included if at least one of the connected elements
!   belongs to one of groups in the array groups.
    integer, dimension(:), intent(in), optional :: groups

!   degrees of freedom in points, on curves and surfaces and in nodesets and
!   elementsets to include in the subscript.
!   For example curves=(/1,4/) includes all degrees on curves 1 and 4.
!   Default is to include all degrees in all nodes.
!   NOTE: the sequence of the degrees is following the nodal numbering of
!         of the points, curves, surfaces and/or nodesets and NOT the global
!         numbering.
!   NOTE: nodes that appear more than once in the list (for example in
!         connecting surfaces or points+curves) will get subscripted only
!         once. See optional argument includenodesonce how to change this.
!   NOTE: elementsets are required to have nodes added using add_to_mesh.
    integer, intent(in), dimension(:), optional :: points, curves, surfaces
    integer, intent(in), dimension(:), optional :: nodesets, elementsets

!   exclude degrees of freedom in points.
!   For example excludepoints=(/1,3/) excludes all degrees of
!   freedom in points 1 and 3.
    integer, intent(in), dimension(:), optional :: excludepoints

!   exclude degrees of freedom on curves.
!   For example excludecurves=(/1,3/) excludes all degrees on curves 1 and 3.
    integer, intent(in), dimension(:), optional :: excludecurves

!   exclude degrees of freedom on surfaces
!   for example excludesurfaces=(/5/) excludes all degrees on surface 5.
    integer, intent(in), dimension(:), optional :: excludesurfaces

!   exclude degrees of freedom in nodes of elementsets
!   for example excludeelementsets=(/2/) excludes all degrees in nodes
!   of elementset 2.
!   NOTE: elementsets are required to have nodes added using add_to_mesh.
    integer, intent(in), dimension(:), optional :: excludeelementsets

!   exclude degrees of freedom in nodes of nodesets
!   for example excludenodesets=(/2/) excludes all degrees in nodes
!   of nodeset 2.
    integer, intent(in), dimension(:), optional :: excludenodesets

!   exclude degrees of freedom in nodes if any of the coordinates is less
!   than the ones given by xmin or any of the coordinates is larger than
!   the ones given by xmax.
!   For example xmin=(/-1._dp,-2._dp/) excludes all nodes that have an
!   x-coordinate less than -1 and/or a y-coordinate less than -2.
    real(dp), intent(in), dimension(:), optional :: xmin, xmax

!   if elementwise is present the subscript is made according to a vector that
!   is stored elementwise. default=.false.
    logical, intent(in), optional :: elementwise

!   logical to indicate whether the nodes involved in the subscript should
!   put into the subscript%nodes.
!   default = .false.
!   NOTE: only nodes that have an non-zero number of degrees of freedom in
!   the subcript are stored.
!   NOTE: nodes that appear more than once in defining the subscript
!   (for example in connecting surfaces or points+curves) will have a single
!   entry in subscript%nodes. See optional argument includenodesonce how to
!   change this.
    logical, intent(in), optional :: fillnodes

!   logical to indicate that nodes that appear more than once in defining
!   the subscript (for example in connecting surfaces or points+curves) will
!   be included only once if set to .true.
!   default = .true.
!   NOTE: this optional argument is to create compatibility with previous
!   versions by setting includenodesonce=.false.
    logical, intent(in), optional :: includenodesonce


!   this routine computes the positions (vector subscript) of degrees of
!   freedom in the vector given by vec or physq.
!   NOTE: if neither degfd of groups is present, all components are included.


    logical :: lelementwise, allnodes, lfillnodes, lincludenodesonce
    logical :: node_with_degrees
    logical, allocatable, dimension(:) :: excludenodes, nodesdone
    integer :: nodenr, bp, nndof, dof, lvec, elgrp, elem, node
    integer :: nnodes, crv, srf, i, nn, nsnodes
    integer, allocatable, dimension(:) :: work, nodes, snodes
    integer :: pos(problem%maxvecnoddegfd), set

    allocate ( excludenodes(mesh%nnodes), nodesdone(mesh%nnodes) )


    call check ( mesh, 'create_subscript_vector' )
    call check ( problem, 'create_subscript_vector', mesh )

    if ( present(vec) .and. present(physq) ) then
      write(*,'(/a/)') &
        'Error in create_subscript_vector: both vec and physq present'
      stop
    else if ( .not. present(vec) .and. .not. present(physq) ) then
      write(*,'(/a/)') &
        'Error in create_subscript_vector: either vec or physq must be present'
      stop
    end if

    if ( present(physq) ) then
      if ( problem%nphysq == 0 ) then
        write(*,'(/a/a/)') &
          'Error in create_subscript_vector: physq specified, however ', &
          ' there are no physical quantities defined'
        stop
      end if
      if ( physq < 1 .or. physq > problem%nphysq ) then
        write(*,'(/a/)') &
          'Error in create_subscript_vector: physq out of range'
        stop
      end if
      lvec = problem%physq(physq)
    else if ( present(vec) ) then
      if ( vec < 1 .or. vec > problem%nvec ) then
        write(*,'(/a/)') &
          'Error in create_subscript_vector: vec out of range'
        stop
      end if
      lvec = vec
    end if

!   elementwise vector?
    lelementwise = set_optional ( variable=elementwise, default=.false. )

    if ( lelementwise .and. ( present(points) .or. present(curves) .or. &
                     present(surfaces) .or. present(excludepoints) .or. &
                     present(excludecurves) .or. present(excludesurfaces) .or. &
                     present(xmin) .or. present(xmax) .or. &
                     present(fillnodes) ) ) then
      write(*,'(4(/a)/)') &
        'Error in create_subscript_vector: for vectors defined elementwise', &
        ' it not possible to specify also one or more of the parameters:', &
        '   points, curves, surfaces, excludepoints, excludecurves', &
        '   excludesurfaces, xmin, xmax, fillnodes'
      stop
    end if

!   subscript already allocated?
    if (allocated(subscript%s)) deallocate(subscript%s)
    if (allocated(subscript%nodes)) deallocate(subscript%nodes)

!   include nodes?
    lfillnodes = set_optional ( variable=fillnodes, default=.false. )

!   include nodes once?
    lincludenodesonce = set_optional ( variable=includenodesonce, &
      default=.true. )

!   check groups
    if ( present(groups) ) then
      if ( any ( groups < 1 ) .or. any ( groups > mesh%nelgrp ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: parameter groups in the heading of create_subscript_vector', &
          ' is out of range: some groups are < 1 or larger than ', mesh%nelgrp
        stop
      end if
    end if

!   check degfd and degsfd
    if ( present(degfd) .and. present(degsfd) ) then
      write(*,'(2(/a)/)') &
        'Error in create_subscript_vector: ',&
        '  Both degfd and degsfd present in heading.'
      stop
    end if

!   determine nnodes

    nnodes = 0

    if ( present(points) ) then
      if ( any ( points < 1 ) .or. any ( points > mesh%npoints ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: parameter points in the heading of create_subscript_vector', &
          ' is out of range: some points are < 1 or larger than ', mesh%npoints
        stop
      end if
      nnodes = nnodes + size(points)
    end if

    if ( present(curves) ) then
      if ( any ( curves < 1 ) .or. any ( curves > mesh%ncurves ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: parameter curves in the heading of create_subscript_vector', &
          ' is out of range: some curves are < 1 or larger than ', mesh%ncurves
        stop
      end if
      nnodes = nnodes + sum(mesh%curves(curves)%nnodes)
    end if

    if ( present(surfaces) ) then
      if ( any ( surfaces < 1 ) .or. any ( surfaces > mesh%nsurfaces ) ) then
        write(*,'(/a/a/a,i0/)') &
          'Error: parameter surfaces in the heading of ', &
          ' create_subscript_vector is ', &
          ' out of range: some surfaces are < 1 or larger than ', mesh%nsurfaces
        stop
      end if
      nnodes = nnodes + sum(mesh%surfaces(surfaces)%nnodes)
    end if

    if ( present(nodesets) ) then
      if ( any ( nodesets < 1 ) .or. any ( nodesets > mesh%nnodesets ) ) then
        write(*,'(/a/a/a,i0/)') &
          'Error: argument nodesets in the heading of', &
          ' create_subscript_vector is ', &
          ' out of range: some nodesets are < 1 or larger than ', mesh%nnodesets
        stop
      end if
      do i = 1, size(nodesets)
        set = nodesets(i)
        nnodes = nnodes + size(mesh%nodesets(set)%a)
      end do
    end if

    if ( present(elementsets) ) then
      if ( any ( elementsets < 1 ) .or. &
           any ( elementsets > mesh%nelementsets ) ) then
        write(*,'(/a/a/a,i0/)') &
          'Error: argument elementsets in the heading of', &
          ' create_subscript_vector is ', &
          ' out of range: some elementsets are < 1 or larger than ', &
          mesh%nelementsets
        stop
      end if
      if ( any ( .not. mesh%elementsets(elementsets)%nodes_created ) ) then
        write(*,'(/4(a/))') &
          'Error: argument elementsets in the heading of', &
          ' create_subscript_vector is ', &
          ' invalid: for some elementsets the array nodes is not filled.', &
          ' Add nodes to the invalid elementsets with add_to_mesh '
        stop
      end if
      nnodes = nnodes + sum(mesh%elementsets(elementsets)%nnodes)
    end if

    allnodes = nnodes == 0

    if ( allnodes ) nnodes = mesh%nnodes

    allocate ( nodes(nnodes) )
    if ( lfillnodes ) allocate ( snodes(nnodes) )


!   determine nodes

    if ( allnodes ) then

!     all nodes

      nodes = [ (i,i=1,mesh%nnodes) ]

    else

!     nodes given by points, curves and/or surfaces

      nn = 0

      if ( present(points) ) then
        nodes(nn+1:nn+size(points)) = mesh%points(points)
        nn = nn + size(points)
      end if

      if ( present(curves) ) then
        do i = 1, size(curves)
          crv = curves(i)
          nodes(nn+1:nn+mesh%curves(crv)%nnodes) = mesh%curves(crv)%nodes
          nn = nn + mesh%curves(crv)%nnodes
        end do
      end if

      if ( present(surfaces) ) then
        do i = 1, size(surfaces)
          srf = surfaces(i)
          nodes(nn+1:nn+mesh%surfaces(srf)%nnodes) = &
                                               mesh%surfaces(srf)%nodes
          nn = nn + mesh%surfaces(srf)%nnodes
        end do
      end if

      if ( present(nodesets) ) then
        do i = 1, size(nodesets)
          set = nodesets(i)
          nodes(nn+1:nn+size(mesh%nodesets(set)%a)) = mesh%nodesets(set)%a
          nn = nn + size(mesh%nodesets(set)%a)
        end do
      end if

      if ( present(elementsets) ) then
        do i = 1, size(elementsets)
          set = elementsets(i)
          nodes(nn+1:nn+mesh%elementsets(set)%nnodes) = &
                                               mesh%elementsets(set)%nodes
          nn = nn + mesh%elementsets(set)%nnodes
        end do
      end if

      if ( nn /= nnodes ) stop 'internal error nn /= nnodes'

    end if

!   excludenodes

    call set_excludenodes


!   fill subscript for node for node

    allocate ( work(problem%maxvecnoddegfd*mesh%maxelnumnod*mesh%nelem) )

    dof = 0
    nsnodes = 0
    nodesdone = .false.

    if ( lelementwise ) then

!     vector stored elementwise

      bp = 0

      do elgrp = 1, problem%nelgrp

        if ( present(groups) ) then
!         check whether group is in included groups
          if ( .not. any( groups == elgrp ) ) cycle
        end if

        do elem = 1, mesh%grpnumel(elgrp)
          do node = 1, mesh%elnumnod(elgrp)

            nndof = problem%vec_elnumdegfd(elgrp)%a(node,lvec)

            call fill_work

            bp = bp + nndof

          end do
        end do

      end do

    else

!     vector stored nodalwise

      do node = 1, nnodes

        nodenr = nodes(node)

        if ( present(groups) ) then
!         check whether nodal point is in included groups
          if ( .not. node_in_groups ( mesh, groups, nodenr ) ) cycle
        end if

        if ( present(layer) ) then
!         cycle if layer absent in node
          if ( .not. btest(problem%nodlayers(nodenr),layer-1) ) cycle
        end if

        if ( excludenodes(nodenr) ) cycle

        if ( lincludenodesonce .and. nodesdone(nodenr) ) cycle

!       get positions of degrees

        call pos_array_vec_node ( problem, nodenr, nndof, pos, lvec, layer )

        if ( nndof > 0 ) bp = pos(1) - 1

        call fill_work

        if ( lfillnodes .and. node_with_degrees ) then
          nsnodes = nsnodes + 1
          snodes(nsnodes) = nodenr
        end if

        nodesdone(nodenr) = .true.

      end do

    end if

    if ( dof == 0 ) then
      write(*,'(/a/)') &
        'Warning in create_subscript_vector: subscript is empty '
    end if

    allocate(subscript%s(dof))

    subscript%s = work(1:dof)

    deallocate ( nodes, work )

    if ( lfillnodes ) then
      allocate(subscript%nodes(nsnodes))
      subscript%nodes = snodes(1:nsnodes)
      deallocate ( snodes )
    end if

    deallocate ( excludenodes, nodesdone )

  contains

    subroutine fill_work

      integer :: deg

      node_with_degrees = .false.

      if ( present(degfd) ) then

!       one degree

        if ( degfd >=1 .and. degfd <= nndof ) then

!         degfd in valid range

!         position in vector
          work(dof+1) = bp + degfd

          dof = dof + 1

          node_with_degrees = .true.

        end if

      else if ( present(degsfd) ) then

!       multiple degrees

        do deg = 1, size(degsfd)

          if ( degsfd(deg) >=1 .and. degsfd(deg) <= nndof ) then

!           degfd in valid range

!           position in vector
            work(dof+1) = bp + degsfd(deg)

            dof = dof + 1

            node_with_degrees = .true.

          end if

        end do

      else

!       all degrees

        do deg = 1, nndof

!         position in vector
          work(dof+1) = bp + deg

          dof = dof + 1

          node_with_degrees = .true.

        end do

      end if

    end subroutine fill_work


!   set excludenodes from points, curves and surfaces

    subroutine set_excludenodes

      integer :: crv, srf, node, set

      excludenodes = .false.

      if ( present(excludepoints) ) then
        if ( any(excludepoints <=0) .or. &
             any(excludepoints > mesh%npoints) ) then
          write(*,'(2(/a),i0/)') &
            'Error in create_subscript_vector: ',&
            ' excludepoints <=0 or excludepoints > number of points = ', &
             mesh%npoints
          stop
        end if
        excludenodes ( mesh%points(excludepoints) ) = .true.
      end if

      if ( present(excludecurves) ) then
        if ( any(excludecurves <=0) .or. &
             any(excludecurves > mesh%ncurves) ) then
          write(*,'(2(/a),i0/)') &
            'Error in create_subscript_vector: ',&
            ' excludecurves <=0 or excludecurves > number of curves = ', &
             mesh%ncurves
          stop
        end if
        do crv = 1, size(excludecurves)
          excludenodes ( mesh%curves(excludecurves(crv))%nodes ) = .true.
        end do
      end if

      if ( present(excludesurfaces) ) then
        if ( any(excludesurfaces <=0) .or. &
             any(excludesurfaces > mesh%nsurfaces) ) then
          write(*,'(2(/a),i0/)') &
            'Error in create_subscript_vector: ',&
            ' excludesurfaces <=0 or excludesurfaces > number of surfaces = ', &
             mesh%nsurfaces
          stop
        end if
        do srf = 1, size(excludesurfaces)
          excludenodes ( mesh%surfaces(excludesurfaces(srf))%nodes ) = .true.
        end do
      end if

      if ( present(excludenodesets) ) then
        if ( any(excludenodesets <=0) .or. &
             any(excludenodesets > mesh%nnodesets) ) then
          write(*,'(2(/a),i0/)') &
            'Error in create_subscript_vector: ',&
            ' excludenodesets <=0 or excludenodesets > number of nodesets = ', &
             mesh%nnodesets
          stop
        end if
        do set = 1, size(excludenodesets)
          excludenodes ( mesh%nodesets(excludenodesets(set))%a ) = .true.
        end do
      end if

      if ( present(excludeelementsets) ) then
        if ( any(excludeelementsets <=0) .or. &
             any(excludeelementsets > mesh%nelementsets) ) then
          write(*,'(2(/2a),i0/)') &
            'Error in create_subscript_vector: ',&
            ' excludeelementsets <=0 or excludeelementsets > ', &
            'number of elementsets = ', mesh%nelementsets
          stop
        end if
        if ( any ( &
              .not. mesh%elementsets(excludeelementsets)%nodes_created ) ) then
          write(*,'(/4(a/))') &
            'Error: argument excludeelementsets in the heading of', &
            ' create_subscript_vector is ', &
            ' invalid: for some elementsets the array nodes is not filled.', &
            ' Add nodes to the invalid elementsets with add_to_mesh '
          stop
        end if
        do set = 1, size(excludeelementsets)
          excludenodes ( &
                     mesh%elementsets(excludeelementsets(set))%nodes ) = .true.
        end do
      end if

      if ( present(xmin) ) then
        if ( size(xmin) /= mesh%ndim ) then
          write(*,'(2(/a),i0/)') &
            'Error in create_subscript_vector: ',&
            ' dimension of xmin must be equal to the space dimension = ', &
             mesh%ndim
          stop
        end if
        do node = 1, mesh%nnodes
          if ( all ( mesh%coor(node,:) > xmin ) ) cycle  ! in region
          excludenodes ( node ) = .true.
        end do
      end if

      if ( present(xmax) ) then
        if ( size(xmax) /= mesh%ndim ) then
          write(*,'(2(/a),i0/)') &
            'Error in create_subscript_vector: ',&
            ' dimension of xmax must be equal to the space dimension = ', &
             mesh%ndim
          stop
        end if
        do node = 1, mesh%nnodes
          if ( all ( mesh%coor(node,:) < xmax ) ) cycle  ! in region
          excludenodes ( node ) = .true.
        end do
      end if

    end subroutine set_excludenodes


  end subroutine create_subscript_vector


! delete vector subscript array in vector

  subroutine delete_single_subscript_vector ( subscript )

!   the vector subscript of the degrees of freedom in a sysvector
    type(subscriptvec_t), intent(out) :: subscript

  end subroutine delete_single_subscript_vector


! delete vector subscript array in vector

  subroutine delete_subscript_vector ( subscript1, subscript2, subscript3, &
    subscript4, subscript5, subscript6, subscript7, subscript8, subscript9, &
    subscript10 )

!   the vector subscript of the degrees of freedom in a sysvector
    type(subscriptvec_t), intent(inout) :: subscript1
    type(subscriptvec_t), intent(inout), optional :: subscript2, subscript3, &
      subscript4, subscript5, subscript6, subscript7, subscript8, subscript9, &
      subscript10


    call delete_single_subscript_vector(subscript1)
    if ( present(subscript2) ) call delete_single_subscript_vector(subscript2)
    if ( present(subscript3) ) call delete_single_subscript_vector(subscript3)
    if ( present(subscript4) ) call delete_single_subscript_vector(subscript4)
    if ( present(subscript5) ) call delete_single_subscript_vector(subscript5)
    if ( present(subscript6) ) call delete_single_subscript_vector(subscript6)
    if ( present(subscript7) ) call delete_single_subscript_vector(subscript7)
    if ( present(subscript8) ) call delete_single_subscript_vector(subscript8)
    if ( present(subscript9) ) call delete_single_subscript_vector(subscript9)
    if ( present(subscript10) ) call delete_single_subscript_vector(subscript10)

  end subroutine delete_subscript_vector


end module vector_m
