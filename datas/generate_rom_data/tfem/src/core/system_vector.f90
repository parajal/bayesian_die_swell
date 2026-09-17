
! Copyright (C) 2004-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Routines, for the system of equations (vectors)

module system_vector_m

  use kind_defs_m
  use mesh_m
  use problem_defs_m
  use pos_array_m
  use sparse_m
  use system_defs_m
  use element_defs_m
  use set_optional_m
  use misc_m, only: sort

  implicit none


! interface for generic create subroutine

  interface create
    module procedure create_sysvector, create_msysvector, &
      create_subscript_sysvector, create_m2sysvector, create_m3sysvector
  end interface create

  interface create_subscript
    module procedure create_subscript_sysvector
  end interface create_subscript


! interface for generic delete subroutine

  interface delete
    module procedure delete_sysvector, delete_msysvector, delete_subscript, &
      delete_m2sysvector, delete_m3sysvector
  end interface delete


! interface for generic copy subroutine

  interface copy
    module procedure copy_sysvector, copy_msysvector, copy_m2sysvector, &
      copy_m3sysvector
  end interface copy


! interface for generic clear_rows subroutine

  interface clear_rows
    module procedure clear_rows_sysvector, clear_rows_msysvector, &
      clear_rows_m2sysvector, clear_rows_m3sysvector
  end interface clear_rows


! interface for generic transform_to_global subroutine

  interface transform_to_global
    module procedure transform_sysvector_to_global, &
      transform_msysvector_to_global, transform_m2sysvector_to_global, &
      transform_m3sysvector_to_global
  end interface transform_to_global


! interface for generic transform_to_local subroutine

  interface transform_to_local
    module procedure transform_sysvector_to_local, &
      transform_msysvector_to_local, transform_m2sysvector_to_local, &
      transform_m3sysvector_to_local
  end interface transform_to_local

contains


! Create single system vector

  subroutine create_single_sysvector ( problem, sysvector, nr )

    type(problem_t),   intent(in)  :: problem
    type(sysvector_t), intent(inout) :: sysvector
    integer, intent(in) :: nr

    if ( sysvector%created ) then
      write(*,'(/a/a,i0/)') &
        'Error in create_sysvector: sysvector has already been created ', &
        ' Vector number in heading = ', nr
      stop
    end if

    sysvector%probnr = problem%probnr
    sysvector%n = problem%numdegfd

    allocate(sysvector%u(sysvector%n))

    sysvector%created = .true.

  end subroutine create_single_sysvector


! Create system vector

  subroutine create_sysvector ( problem, sysvector1, sysvector2, sysvector3, &
    sysvector4, sysvector5 )

    type(problem_t),   intent(in)  :: problem
    type(sysvector_t), intent(inout) :: sysvector1
    type(sysvector_t), intent(inout), optional :: sysvector2, sysvector3, &
      sysvector4, sysvector5

    call check ( problem, 'create_sysvector' )

    call create_single_sysvector( problem, sysvector1, 1 )
    if ( present(sysvector2) ) &
      call create_single_sysvector( problem, sysvector2, 2 )
    if ( present(sysvector3) ) &
      call create_single_sysvector( problem, sysvector3, 3 )
    if ( present(sysvector4) ) &
      call create_single_sysvector( problem, sysvector4, 4 )
    if ( present(sysvector5) ) &
      call create_single_sysvector( problem, sysvector5, 5 )

  end subroutine create_sysvector



! Create (single) multiple system vector

  subroutine create_single_msysvector ( problem, msysvector, nr )

    type(problem_t),   intent(in)  :: problem
    type(sysvector_t), dimension(:), intent(inout) :: msysvector
    integer, intent(in) :: nr

    integer :: isys

    if ( any(msysvector%created) ) then
      write(*,'(/3(a/),a,i0/)') &
        'Error in create_msysvector:', &
        ' some sysvectors in msysvector have already been created ', &
        ' msysvector cannot be created ', &
        ' msysvector number in heading = ', nr
      stop
    end if

    msysvector%probnr = problem%probnr
    msysvector%n = problem%numdegfd

    do isys = 1, size(msysvector)
      allocate(msysvector(isys)%u(problem%numdegfd))
    end do

    msysvector%created = .true.

  end subroutine create_single_msysvector


! Create multiple system vector

  subroutine create_msysvector ( problem, msysvector1, msysvector2, &
    msysvector3, msysvector4, msysvector5 )

    type(problem_t),   intent(in)  :: problem
    type(sysvector_t), dimension(:), intent(inout) :: msysvector1
    type(sysvector_t), dimension(:), intent(inout), optional :: msysvector2, &
      msysvector3, msysvector4, msysvector5

    call check ( problem, 'create_msysvector' )

    call create_single_msysvector( problem, msysvector1, 1 )
    if ( present(msysvector2) ) &
      call create_single_msysvector( problem, msysvector2, 2 )
    if ( present(msysvector3) ) &
      call create_single_msysvector( problem, msysvector3, 3 )
    if ( present(msysvector4) ) &
      call create_single_msysvector( problem, msysvector4, 4 )
    if ( present(msysvector5) ) &
      call create_single_msysvector( problem, msysvector5, 5 )

  end subroutine create_msysvector


! Create (single) matrix of system vectors

  subroutine create_single_m2sysvector ( problem, m2sysvector, nr )

    type(problem_t),   intent(in)  :: problem
    type(sysvector_t), dimension(:,:), intent(inout) :: m2sysvector
    integer, intent(in) :: nr

    integer :: isys1, isys2

    if ( any(m2sysvector%created) ) then
      write(*,'(/3(a/),a,i0/)') &
        'Error in create_m2sysvector:', &
        ' some sysvectors in m2sysvector have already been created ', &
        ' m2sysvector cannot be created ', &
        ' m2sysvector number in heading = ', nr
      stop
    end if

    m2sysvector%probnr = problem%probnr
    m2sysvector%n = problem%numdegfd

    do isys1 = 1, size(m2sysvector,1)
      do isys2 = 1, size(m2sysvector,2)
        allocate(m2sysvector(isys1,isys2)%u(problem%numdegfd))
      end do
    end do

    m2sysvector%created = .true.

  end subroutine create_single_m2sysvector


! Create matrix of system vectors

  subroutine create_m2sysvector ( problem, m2sysvector1, m2sysvector2, &
    m2sysvector3, m2sysvector4, m2sysvector5 )

    type(problem_t),   intent(in)  :: problem
    type(sysvector_t), dimension(:,:), intent(inout) :: m2sysvector1
    type(sysvector_t), dimension(:,:), intent(inout), optional :: &
      m2sysvector2, m2sysvector3, m2sysvector4, m2sysvector5

    call check ( problem, 'create_m2sysvector' )

    call create_single_m2sysvector( problem, m2sysvector1, 1 )
    if ( present(m2sysvector2) ) &
      call create_single_m2sysvector( problem, m2sysvector2, 2 )
    if ( present(m2sysvector3) ) &
      call create_single_m2sysvector( problem, m2sysvector3, 3 )
    if ( present(m2sysvector4) ) &
      call create_single_m2sysvector( problem, m2sysvector4, 4 )
    if ( present(m2sysvector5) ) &
      call create_single_m2sysvector( problem, m2sysvector5, 5 )

  end subroutine create_m2sysvector


! Create (single) 3D array of system vectors

  subroutine create_single_m3sysvector ( problem, m3sysvector, nr )

    type(problem_t),   intent(in)  :: problem
    type(sysvector_t), dimension(:,:,:), intent(inout) :: m3sysvector
    integer, intent(in) :: nr

    integer :: isys1, isys2, isys3

    if ( any(m3sysvector%created) ) then
      write(*,'(/3(a/),a,i0/)') &
        'Error in create_m3sysvector:', &
        ' some sysvectors in m3sysvector have already been created ', &
        ' m3sysvector cannot be created ', &
        ' m3sysvector number in heading = ', nr
      stop
    end if

    m3sysvector%probnr = problem%probnr
    m3sysvector%n = problem%numdegfd

    do isys1 = 1, size(m3sysvector,1)
      do isys2 = 1, size(m3sysvector,2)
        do isys3 = 1, size(m3sysvector,3)
          allocate(m3sysvector(isys1,isys2,isys3)%u(problem%numdegfd))
        end do
      end do
    end do

    m3sysvector%created = .true.

  end subroutine create_single_m3sysvector


! Create 3D array of system vectors

  subroutine create_m3sysvector ( problem, m3sysvector1, m3sysvector2, &
    m3sysvector3, m3sysvector4, m3sysvector5 )

    type(problem_t),   intent(in)  :: problem
    type(sysvector_t), dimension(:,:,:), intent(inout) :: m3sysvector1
    type(sysvector_t), dimension(:,:,:), intent(inout), optional :: &
      m3sysvector2, m3sysvector3, m3sysvector4, m3sysvector5

    call check ( problem, 'create_m3sysvector' )

    call create_single_m3sysvector( problem, m3sysvector1, 1 )
    if ( present(m3sysvector2) ) &
      call create_single_m3sysvector( problem, m3sysvector2, 2 )
    if ( present(m3sysvector3) ) &
      call create_single_m3sysvector( problem, m3sysvector3, 3 )
    if ( present(m3sysvector4) ) &
      call create_single_m3sysvector( problem, m3sysvector4, 4 )
    if ( present(m3sysvector5) ) &
      call create_single_m3sysvector( problem, m3sysvector5, 5 )

  end subroutine create_m3sysvector


! Delete single system vector

  subroutine delete_single_sysvector ( sysvector, nr )

    type(sysvector_t), intent(inout) :: sysvector
    integer, intent(in) :: nr

    if ( .not. sysvector%created ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error in delete_sysvector:', &
        ' sysvector has not been created and cannot be deleted ', &
        ' sysvector number in heading = ', nr
      stop
    end if

!   deallocate all allocatables and intialize to defaults
    call deall ( sysvector )

  contains

    subroutine deall ( sysvector )
      type(sysvector_t), intent(out) :: sysvector
    end subroutine deall

  end subroutine delete_single_sysvector


! Delete system vector

  subroutine delete_sysvector ( sysvector1, sysvector2, sysvector3, &
    sysvector4, sysvector5 )

    type(sysvector_t), intent(inout) :: sysvector1
    type(sysvector_t), intent(inout), optional :: sysvector2, sysvector3, &
      sysvector4, sysvector5

    call delete_single_sysvector(sysvector1,1)
    if ( present(sysvector2) ) call delete_single_sysvector(sysvector2,2)
    if ( present(sysvector3) ) call delete_single_sysvector(sysvector3,3)
    if ( present(sysvector4) ) call delete_single_sysvector(sysvector4,4)
    if ( present(sysvector5) ) call delete_single_sysvector(sysvector5,5)

  end subroutine delete_sysvector


! Delete (single) multiple system vector

  subroutine delete_single_msysvector ( msysvector, nr )

    type(sysvector_t), dimension(:), intent(inout) :: msysvector
    integer, intent(in) :: nr

    if ( .not. all(msysvector%created) ) then
      write(*,'(/3(a/),a,i0/)') &
        'Error in delete_msysvector:', &
        ' not all sysvectors in msysvector have been created ', &
        ' msysvector cannot be deleted ', &
        ' msysvector number in heading = ', nr
      stop
    end if

!   deallocate all allocatables and intialize to defaults
    call deall ( msysvector )

  contains

    subroutine deall ( msysvector )
      type(sysvector_t), dimension(:), intent(out) :: msysvector
    end subroutine deall

  end subroutine delete_single_msysvector


! Delete multiple system vector

  subroutine delete_msysvector ( msysvector1, msysvector2, msysvector3, &
    msysvector4, msysvector5 )

    type(sysvector_t), dimension(:), intent(inout) :: msysvector1
    type(sysvector_t), dimension(:), intent(inout), optional :: msysvector2, &
      msysvector3, msysvector4, msysvector5

    call delete_single_msysvector(msysvector1,1)
    if ( present(msysvector2) ) call delete_single_msysvector(msysvector2,2)
    if ( present(msysvector3) ) call delete_single_msysvector(msysvector3,3)
    if ( present(msysvector4) ) call delete_single_msysvector(msysvector4,4)
    if ( present(msysvector5) ) call delete_single_msysvector(msysvector5,5)

  end subroutine delete_msysvector


! Delete (single) matrix of system vectors

  subroutine delete_single_m2sysvector ( m2sysvector, nr )

    type(sysvector_t), dimension(:,:), intent(inout) :: m2sysvector
    integer, intent(in) :: nr

    if ( .not. all(m2sysvector%created) ) then
      write(*,'(/3(a/),a,i0/)') &
        'Error in delete_m2sysvector:', &
        ' not all sysvectors in m2sysvector have been created ', &
        ' m2sysvector cannot be deleted ', &
        ' m2sysvector number in heading = ', nr
      stop
    end if

!   deallocate all allocatables and intialize to defaults
    call deall ( m2sysvector )

  contains

    subroutine deall ( m2sysvector )
      type(sysvector_t), dimension(:,:), intent(out) :: m2sysvector
    end subroutine deall

  end subroutine delete_single_m2sysvector


! Delete matrix of system vectors

  subroutine delete_m2sysvector ( m2sysvector1, m2sysvector2, m2sysvector3, &
    m2sysvector4, m2sysvector5 )

    type(sysvector_t), dimension(:,:), intent(inout) :: m2sysvector1
    type(sysvector_t), dimension(:,:), intent(inout), optional :: &
      m2sysvector2, m2sysvector3, m2sysvector4, m2sysvector5

    call delete_single_m2sysvector(m2sysvector1,1)
    if ( present(m2sysvector2) ) call delete_single_m2sysvector(m2sysvector2,2)
    if ( present(m2sysvector3) ) call delete_single_m2sysvector(m2sysvector3,3)
    if ( present(m2sysvector4) ) call delete_single_m2sysvector(m2sysvector4,4)
    if ( present(m2sysvector5) ) call delete_single_m2sysvector(m2sysvector5,5)

  end subroutine delete_m2sysvector


! Delete (single) 3D array of system vectors

  subroutine delete_single_m3sysvector ( m3sysvector, nr )

    type(sysvector_t), dimension(:,:,:), intent(inout) :: m3sysvector
    integer, intent(in) :: nr

    if ( .not. all(m3sysvector%created) ) then
      write(*,'(/3(a/),a,i0/)') &
        'Error in delete_m3sysvector:', &
        ' not all sysvectors in m3sysvector have been created ', &
        ' m3sysvector cannot be deleted ', &
        ' m3sysvector number in heading = ', nr
      stop
    end if

!   deallocate all allocatables and intialize to defaults
    call deall ( m3sysvector )

  contains

    subroutine deall ( m3sysvector )
      type(sysvector_t), dimension(:,:,:), intent(out) :: m3sysvector
    end subroutine deall

  end subroutine delete_single_m3sysvector


! Delete 3D array of system vectors

  subroutine delete_m3sysvector ( m3sysvector1, m3sysvector2, m3sysvector3, &
    m3sysvector4, m3sysvector5 )

    type(sysvector_t), dimension(:,:,:), intent(inout) :: m3sysvector1
    type(sysvector_t), dimension(:,:,:), intent(inout), optional :: &
      m3sysvector2, m3sysvector3, m3sysvector4, m3sysvector5

    call delete_single_m3sysvector(m3sysvector1,1)
    if ( present(m3sysvector2) ) call delete_single_m3sysvector(m3sysvector2,2)
    if ( present(m3sysvector3) ) call delete_single_m3sysvector(m3sysvector3,3)
    if ( present(m3sysvector4) ) call delete_single_m3sysvector(m3sysvector4,4)
    if ( present(m3sysvector5) ) call delete_single_m3sysvector(m3sysvector5,5)

  end subroutine delete_m3sysvector


! Copy system vector to another system vector

  subroutine copy_sysvector ( sysvector1, sysvector2 )

    type(sysvector_t), intent(inout) :: sysvector1, sysvector2

!   This routine performs a real data copy of one vector to the other using
!   the assignment statement:
!
!      sysvector2 = sysvector1
!

    if ( .not. sysvector1%created ) then
      write(*,'(/a/)') &
        'Error in copy_sysvector:  sysvector1 has not been created '
      stop
    end if

    sysvector2 = sysvector1

  end subroutine copy_sysvector


! Copy multiple system vectors to another system vector

  subroutine copy_msysvector ( msysvector1, msysvector2 )

    type(sysvector_t), dimension(:), intent(inout) :: msysvector1, msysvector2

    integer :: i

    if ( size(msysvector1) /= size(msysvector2) ) then
      write(*,'(/2a/)') &
        'Error in copy_msysvector: ', &
        ' dimension of msysvector1 and msysvector2 is different '
      stop
    end if

    do i = 1, size(msysvector1)
      call copy_sysvector ( msysvector1(i), msysvector2(i) )
    end do

  end subroutine copy_msysvector


! Copy multiple (matrix) system vectors to another multiple system vector

  subroutine copy_m2sysvector ( m2sysvector1, m2sysvector2 )

    type(sysvector_t), dimension(:,:), intent(inout) :: m2sysvector1, &
      m2sysvector2

    integer :: i, j

    if ( size(m2sysvector1,1) /= size(m2sysvector2,1) .or. &
         size(m2sysvector1,2) /= size(m2sysvector2,2) ) then
      write(*,'(/2a/)') &
        'Error in copy_m2sysvector: ', &
        ' dimension of m2sysvector1 and m2sysvector2 is different '
      stop
    end if

    do i = 1, size(m2sysvector1,1)
      do j = 1, size(m2sysvector1,2)
        call copy_sysvector ( m2sysvector1(i,j), m2sysvector2(i,j) )
      end do
    end do

  end subroutine copy_m2sysvector


! Copy multiple (3D array) system vectors to another multiple system vector

  subroutine copy_m3sysvector ( m3sysvector1, m3sysvector2 )

    type(sysvector_t), dimension(:,:,:), intent(inout) :: m3sysvector1, &
      m3sysvector2

    integer :: i, j, k

    if ( size(m3sysvector1,1) /= size(m3sysvector2,1) .or. &
         size(m3sysvector1,2) /= size(m3sysvector2,2) .or. &
         size(m3sysvector1,3) /= size(m3sysvector2,3) ) then
      write(*,'(/2a/)') &
        'Error in copy_m3sysvector: ', &
        ' dimension of m3sysvector1 and m3sysvector2 is different '
      stop
    end if

    do i = 1, size(m3sysvector1,1)
      do j = 1, size(m3sysvector1,2)
        do k = 1, size(m3sysvector1,3)
          call copy_sysvector ( m3sysvector1(i,j,k), m3sysvector2(i,j,k) )
        end do
      end do
    end do

  end subroutine copy_m3sysvector


! Get element values/index of sysvector for the element given by (elgrp,elem)

  subroutine get_sysvector ( mesh, problem, sysvector, elgrp, elem, u, physq, &
    order, layer, posu, ndofu, sloppy )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(sysvector_t), intent(in) :: sysvector
    integer, intent(in) :: elgrp, elem

!   the element degrees of freedom, if present
    real(dp), intent(out), dimension(:), optional :: u

!   if physq is present only physical quantities in this array are included.
!   For example:
!     physq=(/2,1/)
!   means that u will contain the physical quantities 2 and 1, in that order.
    integer, dimension(:), intent(in), optional :: physq

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

!   if layer is present and layer > 0 the element degrees of freedom are
!   restricted to the specified layer.
!   layer=0 is identical to layer not present
    integer, intent(in), optional :: layer

!   if present it contains the positions of the components of u in sysvector
    integer, dimension(:), intent(out), optional :: posu

!   if present it gives the number of element degrees of freedom stored in u
    integer, intent(out), optional :: ndofu

!   sloppy? Don't check nelgrp
    logical, intent(in), optional :: sloppy


    logical, dimension(mesh%elnumnod(elgrp)) :: lp
    integer :: dof, ndof, llayer, i
    integer, dimension(:), allocatable :: pos
    integer, dimension(mesh%elnumnod(elgrp)) :: nodes


    call check ( problem, 'get_sysvector', mesh, sloppy )

    llayer = set_optional ( variable=layer, default=0 )

!   Some testing

    if ( .not. sysvector%created ) then
      write(*,'(/a/)') &
        'Error in get_sysvector: sysvector not created.'
      stop
    end if

    if ( problem%probnr /= sysvector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in get_sysvector: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in sysvector = ',  sysvector%probnr
      stop
    end if

    if ( present(physq) ) then
      if ( problem%nphysq == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in get_sysvector: ',&
          ' specification of physical quantities (physq) ', &
          ' is only available if physical quantities have been defined.'
        stop
      end if
      if ( any( physq < 1 ) .or.  any( physq > problem%nphysq ) ) then
        write(*,'(2(/a)/)') &
          'Error in get_sysvector: ', &
          ' physical quantities in physq are out of range.'
        stop
      end if
    end if

!   ordering; set defaults

    if ( present(order) ) then
      if ( all ( order /= [ 'ND', 'DN' ] ) ) then
        write(*,'(2(/a)/)') &
          'Error in get_sysvector: ', &
          ' heading parameter order must be either ''ND'' or ''DN''.'
        stop
      end if
    end if

!   layers

    if ( present(layer) ) then
      if ( layer < 0 .or. layer > problem%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in get_sysvector: layer out of range.'
        stop
      end if
    end if

!   number of degrees of freedom in element

    nodes = mesh%topology(elgrp)%a(:,elem)

    if ( llayer > 0 ) then

!     single layer

      lp = btest(problem%nodlayers(nodes),llayer-1)

      if ( present(physq) ) then
        ndof = 0
        do i = 1, size(physq)
          ndof = ndof &
                 + sum ( problem%vec_elnumdegfd(elgrp)%a(:,physq(i)), mask=lp  )
        end do
      else
        ndof = sum ( problem%elnumdegfd(elgrp)%a, mask=lp )
      end if

    else

!     standard situation

      if ( present(physq) ) then
        ndof = sum ( problem%vec_nodnumdegfd(nodes+1,&
                                  &problem%physq(physq)) &
                   - problem%vec_nodnumdegfd(nodes,&
                                  &problem%physq(physq)) )
      else
        ndof = sum ( problem%nodnumdegfd(nodes+1) - problem%nodnumdegfd(nodes) )
      end if

    end if

    if ( present(ndofu) ) ndofu = ndof

!   check u array

    if ( present(u) ) then
      if ( size(u) < ndof ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_sysvector: ', &
          ' u array is too small: size(u) = ', size(u), &
          ' whereas the number of degrees of freedom = ', ndof
        stop
      end if
    end if

!   check posu array

    if ( present(posu) ) then
      if ( size(posu) < ndof ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_sysvector: ', &
          ' posu array is too small: size(posu) = ', size(posu), &
          ' whereas the number of degrees of freedom = ', ndof
        stop
      end if
    end if

    allocate(pos(ndof))

!   compute positions in sysvector

    if ( present(physq) ) then
!     physical quantities
      call pos_array ( mesh, problem, elgrp, elem, dof, pos, physq, &
        order, layer )
    else
      call pos_array ( mesh, problem, elgrp, elem, dof, pos, order=order, &
        layer=layer )
    end if

    if ( dof /= ndof ) then
      write(*,'(/a/a/)') &
        'Internal error in get_sysvector: ', &
        ' dof /= ndof '
      stop
    end if

!   compute u

    if ( present(u) ) u(1:ndof) = sysvector%u(pos)

!   posu array

    if ( present(posu) ) posu(1:ndof) = pos

    deallocate(pos)

  end subroutine get_sysvector


! Get element values/index of sysvector for the element on a geometry given by
! (curve, elem), (surface,elem), (volume,elem) or (ndimr,geometry,elem).

  subroutine get_sysvector_geometry ( mesh, problem, sysvector, elem, u, &
    curve, surface, volume, ndimr, geometry, physq, order, layer, posu )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(sysvector_t), intent(in) :: sysvector
    integer, intent(in) :: elem

!   the element degrees of freedom, if present
    real(dp), intent(out), dimension(:), optional :: u

!   the geometry number: one of curve, surface, volume must be present (legacy)
    integer, intent(in), optional :: curve, surface, volume

!   dimension of reference space (ndimr) and geometry number
    integer, intent(in), optional :: ndimr, geometry

!   if physq is present only physical quantities in this array are included.
!   For example:
!     physq=(/2,1/)
!   means that u will contain the physical quantities 2 and 1, in that order.
    integer, dimension(:), intent(in), optional :: physq

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

!   if layer is present and layer > 0 the element degrees of freedom are
!   restricted to the specified layer.
!   layer=0 is identical to layer not present
    integer, intent(in), optional :: layer

!   if present it contains the positions of the components of u in sysvector
    integer, dimension(:), intent(out), optional :: posu


    integer :: dof, ndof, llayer, i
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
          call errormsg_case_default ( 'get_sysvector_geometry', 'ndimr', &
            int_value=ndimr )
      end select
    end if

    llayer = set_optional ( variable=layer, default=0 )

!   Some testing

    call check ( problem, 'get_sysvector_geometry' )

!   test valid geometry

    if ( lcurve /= 0 ) then

      if ( lcurve < 1 .or. lcurve > mesh%ncurves ) then

        write(*,'(/a/a,i0,/a,i0/)') &
          'Error: curve in the heading of get_sysvector_geometry is ', &
          'out of range. curve is ', lcurve, &
          'whereas the number of curves is ', mesh%ncurves
        stop

      end if

      allocate(nodes(mesh%curves(lcurve)%elnumnod))

    else if ( lsurface /= 0 ) then

      if ( lsurface < 1 .or. lsurface > mesh%nsurfaces ) then

        write(*,'(/a/a,i0,/a,i0/)') &
          'Error: surface in the heading of get_sysvector_geometry is ', &
          'out of range. surface is ', lsurface, &
          'whereas the number of surfaces is ', mesh%nsurfaces
        stop

      end if

      allocate(nodes(mesh%surfaces(lsurface)%elnumnod))

    else if ( lvolume /= 0 ) then

      if ( lvolume < 1 .or. lvolume > mesh%nvolumes ) then

        write(*,'(/a/a,i0,/a,i0/)') &
          'Error: volume in the heading of get_sysvector_geometry is ', &
          'out of range. volume is ', lvolume, &
          'whereas the number of volumes is ', mesh%nvolumes
        stop

      end if

      allocate(nodes(mesh%volumes(lvolume)%elnumnod))

    else

      write(*,'(/a/a/)') &
        'Error: either curve, surface, volume or (ndimr,geometry) must be', &
        ' present in the heading of get_sysvector_geometry'
      stop

    end if

    if ( .not. sysvector%created ) then
      write(*,'(/a/)') &
        'Error in get_sysvector_geometry: sysvector not created.'
      stop
    end if

    if ( problem%probnr /= sysvector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in get_sysvector_geometry: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in sysvector = ',  sysvector%probnr
      stop
    end if

!   ordering; set defaults

    if ( present(order) ) then
      if ( all ( order /= [ 'ND', 'DN' ] ) ) then
        write(*,'(2(/a)/)') &
          'Error in get_sysvector_geometry: ', &
          ' heading parameter order must be either ''ND'' and ''DN''.'
        stop
      end if
    end if

    if ( present(physq) ) then

      if ( problem%nphysq == 0 ) then
        write(*,'(4(/a)/)') &
          'Error in get_sysvector_geometry: ',&
          ' retrieving physical quantities physq', &
          ' is only available if physical quantities have been defined.'
        stop
      end if

    end if

    if ( present(physq) ) then
      if ( problem%nphysq == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in get_sysvector_geometry: ',&
          ' specification of physical quantities (physq) ', &
          ' is only available if physical quantities have been defined.'
        stop
      end if
      if ( any( physq < 1 ) .or.  any( physq > problem%nphysq ) ) then
        write(*,'(2(/a)/)') &
          'Error in get_sysvector_geometry: ', &
          ' physical quantities in physq are out of range.'
        stop
      end if
    end if

!   layers

    if ( present(layer) ) then
      if ( layer < 0 .or. layer > problem%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in get_sysvector_geometry: layer out of range.'
        stop
      end if
    end if

!   find geometry and number of degr of freedom in an element of this geometry

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

    if ( llayer > 0 ) then

!     single layer

      allocate ( nl(size(nodes)) )

      where ( btest(problem%nodlayers(nodes),llayer-1) )
        nl = count_layers ( problem%nodlayers(nodes), problem%numlayers )
      else where
        nl = 0
      end where

      if ( present(physq) ) then
        ndof = 0
        do i = 1, size(physq)
          ndof = ndof + &
                   sum ( ( problem%vec_nodnumdegfd(nodes+1,&
                                        &problem%physq(physq(i))) &
                         - problem%vec_nodnumdegfd(nodes,&
                                        &problem%physq(physq(i))) ) / nl, &
                            mask = nl > 0 )
        end do
      else
        ndof = sum ( ( problem%nodnumdegfd(nodes+1) &
                        - problem%nodnumdegfd(nodes) ) / nl, &
                        mask = nl > 0 )
      end if

      deallocate ( nl )

    else

!     standard situation

      if ( present(physq) ) then
        ndof = sum ( problem%vec_nodnumdegfd(nodes+1,&
                                        &problem%physq(physq)) &
                         - problem%vec_nodnumdegfd(nodes,&
                                        &problem%physq(physq)) )
      else
        ndof = sum ( problem%nodnumdegfd(nodes+1) - problem%nodnumdegfd(nodes) )
      end if

    end if

!   check u array

    if ( present(u) ) then
      if ( size(u) < ndof ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_sysvector_geometry: ', &
          ' u array is too small: size(u) = ', size(u), &
          ' whereas the number of degrees of freedom = ', ndof
        stop
      end if
    end if

!   check posu array

    if ( present(posu) ) then
      if ( size(posu) < ndof ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_sysvector_geometry: ', &
          ' posu array is too small: size(posu) = ', size(posu), &
          ' whereas the number of degrees of freedom = ', ndof
        stop
      end if
    end if

    allocate(pos(ndof))

!   compute positions in large vector of element degrees of freedom

    call pos_array_geometry ( problem, geom, elem, dof, pos, &
      physq, order, layer )

!   compute u

    if ( present(u) ) u(1:ndof) = sysvector%u(pos)

!   posu array

    if ( present(posu) ) posu(1:ndof) = pos

    deallocate ( pos, nodes )

  end subroutine get_sysvector_geometry


! Fill system vector with values

  subroutine fill_sysvector ( mesh, problem, sysvector, physq, layer, degfd, &
    degsfd, point, points, curve1, curve2, step, exclude, curves, &
    node1, node2, nodes, group, element, elnode, surface1, surface2, &
    surfaces, elgroup1, elgroup2, elgroups, nodeset1, nodeset2, nodesets, &
    excludepoints, excludecurves, excludesurfaces, value, func, funcnr, &
    vvalue, vfunc, vfuncnr, ignore_nodegfd, blend )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(sysvector_t), intent(inout) :: sysvector

!   if physq is present
!     only physical quantity physq is involved.
!     degfd is with respect to the physical quantity physq
!   else
!     all degrees of freedom are involved
!     degfd is with respect to the all degrees of freedom
    integer, intent(in), optional :: physq

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

!   If element is present, values are set in a element
!   If also group is present, values in the specified element in the
!   specified group is set. Default=1
!   if also elnode is present, values in the local node elnode of the element
!   is set. Default=1.
    integer, intent(in), optional :: element, group, elnode

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
!         exclude=0 no excludes, all nodes
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



!   This routine fills a sysvector with values such as for essential boundary
!   conditions and initial values.


    logical :: excludenodes, lignore_nodegfd
    integer :: worklength, curve, nn, sz, stp, node, nodenr
    integer :: nodeset, nval
    integer :: surface, numec, numes, enode, elem
    integer :: ndegfd, elnod, elgrp
    integer, allocatable, dimension(:) :: wnodes, work, ldegsfd, &
      lcurves, lsurfaces, lnodes, lelgroups, lnodesets
    integer :: pos(problem%maxnoddegfd), i, lblend


    lignore_nodegfd  = &
                set_optional ( variable=ignore_nodegfd, default = .false. )

    if ( present(blend) ) then
      if ( blend < -1 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error fill_sysvector: invalid argument ', &
          '  blend < -1 or larger than ', mesh%nblend
        stop
      end if
      lblend = blend
    else
      lblend = 0 ! default
    end if

!   testing

    call check ( mesh, 'fill_sysvector' )
    call check ( problem, 'fill_sysvector' )

    if ( .not. sysvector%created ) then
      write(*,'(/a/)') &
        'Error in fill_sysvector: sysvector not created.'
      stop
    end if

    if ( problem%probnr /= sysvector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in fill_sysvector: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in sysvector = ',  sysvector%probnr
      stop
    end if

    if ( present(physq) ) then
      if ( problem%nphysq == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in fill_sysvector: ',&
          ' specification of physical quantities (physq) ', &
          ' is only available if physical quantities have been defined.'
        stop
      end if
      if ( physq <=0 .or. physq > problem%nphysq ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_sysvector: ',&
          ' physq <=0 or physq > number of physical quantities = ', &
          problem%nphysq
        stop
      end if
    end if

    if ( present(layer) ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in fill_sysvector: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_sysvector: ',&
          ' layer < 1 or layer > number of layers = ', &
          problem%numlayers
        stop
      end if
    end if

    if ( present(excludepoints) ) then
      if ( any(excludepoints <=0) .or. any(excludepoints > mesh%npoints) ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_sysvector: ',&
          ' excludepoints <=0 or excludepoints > number of points = ', &
           mesh%npoints
        stop
      end if
    end if

    if ( present(excludecurves) ) then
      if ( any(excludecurves <=0) .or. any(excludecurves > mesh%ncurves) ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_sysvector: ',&
          ' excludecurves <=0 or excludecurves > number of curves = ', &
           mesh%ncurves
        stop
      end if
    end if

    if ( present(excludesurfaces) ) then
      if ( any(excludesurfaces <=0) .or. &
           any(excludesurfaces > mesh%nsurfaces) ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_sysvector: ',&
          ' excludesurfaces <=0 or excludesurfaces > number of surfaces = ', &
           mesh%nsurfaces
        stop
      end if
    end if

    if ( present(degfd) .and. present(degsfd) ) then
      write(*,'(2(/a)/)') &
        'Error in fill_sysvector: ',&
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

!     fill for a single point

      if ( point <=0 .or. point > mesh%npoints ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_sysvector: ',&
          ' point <=0 or point > number of points = ', &
           mesh%npoints
        stop
      end if

      wnodes = [mesh%points(point)]

      nn = 1

    else if ( present(points) ) then

      if ( any( points <=0 .or. points > mesh%npoints ) ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_sysvector: ',&
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
            'Error in fill_sysvector: ',&
            ' curve1 <=0 or curve1 > number of curves = ', &
             mesh%ncurves
          stop
        end if

        if ( present(curve2) ) then

          if ( curve2 <=0 .or. curve2 > mesh%ncurves ) then
            write(*,'(2(/a),i0/)') &
              'Error in fill_sysvector: ',&
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
            'Error in fill_sysvector: ',&
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
            'Error: exclude in the heading of fill_sysvector can only ', &
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
            'Error in fill_sysvector: ',&
            ' surface1 <=0 or surface1 > number of surfaces = ', &
             mesh%nsurfaces
          stop
        end if

        if ( present(surface2) ) then

          if ( surface2 <=0 .or. surface2 > mesh%nsurfaces ) then
            write(*,'(2(/a),i0/)') &
              'Error in fill_sysvector: ',&
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
            'Error in fill_sysvector: ',&
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

!         exclude nodes

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
            'Error in fill_sysvector: ',&
            ' node1 <=0 or node1 > number of nodes = ', &
             mesh%nnodes
          stop
        end if

        if ( present(node2) ) then
          if ( node2 <=0 .or. node2 > mesh%nnodes ) then
            write(*,'(2(/a),i0/)') &
              'Error in fill_sysvector: ',&
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
            'Error in fill_sysvector: ',&
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

    else if ( present(element) ) then

!     fill for an element

      elgrp = 1
      if ( present(group) ) then
        if ( group <=0 .or. group > mesh%nelgrp ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_sysvector: ',&
            ' group <=0 or group > number of groups = ', &
             mesh%nelgrp
          stop
        end if
        elgrp = group
      end if

      elnod = 1
      if ( present(elnode) ) then
        if ( elnode <=0 .or. elnode > mesh%elnumnod(elgrp) ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_sysvector: ',&
            ' elnode <=0 or elnode > number of nodes in element = ', &
             mesh%elnumnod(elgrp)
          stop
        end if
        elnod = elnode
      end if

      if ( element <=0 .or. element > mesh%grpnumel(elgrp) ) then
        write(*,'(/a/a,i0,a,i0/)') &
          'Error in fill_sysvector: ',&
          ' element <=0 or element > number of element in group ', elgrp, &
          ' = ', mesh%grpnumel(elgrp)
        stop
      end if

      allocate ( wnodes(1) )

!     determine nodal point

      wnodes(1)= mesh%topology(elgrp)%a(elnod,element)

      nn = 1

    else if ( present(elgroup1) .or. present(elgroups) ) then

!     fill for groups

      if ( present(elgroup1) ) then

        if ( elgroup1 <=0 .or. elgroup1 > mesh%nelgrp ) then
          write(*,'(2(/a),i0/)') &
            'Error in fill_sysvector: ',&
            ' elgroup1 <=0 or elgroup1 > number of groups = ', &
             mesh%nelgrp
          stop
        end if

        if ( present(elgroup2) ) then

          if ( elgroup2 <=0 .or. elgroup2 > mesh%nelgrp ) then
            write(*,'(2(/a),i0/)') &
              'Error in fill_sysvector: ',&
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
            'Error in fill_sysvector: ',&
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
        'fill_sysvector. Non of the keywords:',   &
        '  point or points ',     &
        '  curve1 or curves ',    &
        '  surface1 or surfaces ',  &
        '  node1 or nodes  ',    &
        '  element ',   &
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

!   sort nodes

    call sort ( wnodes(1:nn) )

!   Fill degrees of freedom in wnodes

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

        if ( present(physq) ) then

          call pos_array_node ( problem, nodenr, ndegfd, pos, [physq], layer )

        else

          call pos_array_node ( problem, nodenr, ndegfd, pos, layer=layer )

        end if

!       find degrees

        if ( present(degfd) ) then

!         single degree

          if ( degfd < 1 .or. degfd > ndegfd ) then

            if ( lignore_nodegfd ) cycle

            write(*,'(/a/2(a,i0),/a,i0/)') &
              'Error: degfd in the heading of fill_sysvector is outside', &
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
              'Error: degsfd in the heading of fill_sysvector is outside', &
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

!       specify value

        if ( present(value) ) then

!         value given; takes precedence over function
          sysvector%u(pos(ldegsfd)) = value

        else if ( present(func) ) then

!         function given

          if ( present(funcnr) ) then

            sysvector%u(pos(ldegsfd)) = &
                                  func ( funcnr, mesh%coor(nodenr,:) )

          else

            write(*,'(2(/a)/)') &
              'Error: if func is present in the heading of ', &
              'fill_sysvector, funcnr must be present as well.'
            stop

          end if

        else if ( present(vvalue) ) then

!         vvalue given; takes precedence over vector function

          nval = size(ldegsfd)

          if ( nval <= size(vvalue) ) then

            sysvector%u(pos(ldegsfd)) = vvalue(1:nval)

          else

            write(*,'(/a/2(a,i0)/)') &
              'Error in fill_sysvector: size of vvalue is too small:', &
              '  size = ', size(vvalue), ', required = ', nval
            stop

          end if

        else if ( present(vfunc) ) then

!         vector function given

          if ( present(vfuncnr) ) then

            sysvector%u(pos(ldegsfd)) = &
                        vfunc ( size(ldegsfd), vfuncnr, mesh%coor(nodenr,:) )

          else

            write(*,'(2(/a)/)') &
              'Error: if vfunc is present in the heading of ', &
              'fill_sysvector, vfuncnr must be present as well.'
            stop

          end if

        else

          write(*,'(3(/a)/)') &
            'Error: either value, func, vvalue or vfunc must be present in ', &
            ' the heading of fill_sysvector in order to fill the vector ', &
            ' with values.'
          stop

        end if

      end if

    end do

    deallocate ( wnodes )

  end subroutine fill_sysvector


! Fill system vector using elements

  subroutine fill_sysvector_elem ( mesh, problem, sysvector, elemsub, &
    coefficients, mcoefficients, oldvectors, physq, order, elgroup1, &
    elgroup2, groups, layer, addvec )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the sysvector vector to be filled
    type(sysvector_t), intent(inout) :: sysvector

!   this is the element subroutine that must be supplied by the calling routine
!   Note that the element matrix elemmat and elemvec are adjusted to the actual
!   size of these matrices. Therefore the shape can be used on element level.
    interface
      subroutine elemsub ( mesh, problem, elgrp, elem, first, last, &
        coefficients, oldvectors, elemvec )
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
        real(dp), intent(out), dimension(:) :: elemvec
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

!   if physq is present the filling takes place for physical quantities
!   in physq only. The size of the arrays determine the
!   number of physical quantities involved. The element vector is defined
!   accordingly.
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

!   if these are present the filling takes place for element groups
!   elgroup1,...,elgroup2 only. If only elgroup1 is present one group is
!   filled only.
    integer, intent(in), optional :: elgroup1, elgroup2

!   if present: the element groups to be filled
!   For example groups=(/2,4/) will fill element groups 2 and 4.
!   Default: all groups
    integer, dimension(:), intent(in), optional :: groups

!   if present the filling takes place within the single layer of degrees of
!   freedom only.
!   Only elements where all nodes have degrees in the specified layer are
!   assembled.
    integer, intent(in), optional :: layer


!   if addvec is set to .true. the vector is cleared before
!   the filling and thus the element vectors are filled in an
!   existing system. The default is .false. (clearing)
    logical, intent(in), optional :: addvec


!   This routine performs the filling proces and nothing more. This means that
!   after the first call (with addvec=.false.), the routine can be called
!   as many times as needed with addvec=.true.


    logical :: addvector, first, last
    logical, allocatable, dimension(:) :: wkl
    integer :: elem, elgrp, ndof, dof, grp
    integer :: row, rowg
    integer :: lgroups(mesh%nelgrp), lnelgrp, i
    integer, allocatable, dimension(:) :: pos
    real(dp), allocatable, dimension(:) :: elemvec
    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl

    allocate ( wkl(mesh%nnodes) )

    call check ( mesh, 'fill_sysvector_elem' )
    call check ( problem, 'fill_sysvector_elem' )

    if ( present(mcoefficients) ) then
      if ( size(mcoefficients) /= mesh%nelgrp ) then
        write(*,'(/2(a/))') &
          'Error in fill_sysvector_elem: ', &
          '   size of mcoefficients /= number of element groups'
        stop
      end if
    end if

!   sysvector created?

    if ( .not. sysvector%created ) then
      write(*,'(/a/)') &
        'Error in fill_sysvector_elem: sysvector not created.'
      stop
    end if

!   ordering; set defaults

    if ( present(order) ) then
      if ( all ( order /= [ 'ND', 'DN' ] ) ) then
        write(*,'(2(/a)/)') &
          'Error in fill_sysvector_elem: ', &
          ' heading parameter order must be either ''ND'' or ''DN''.'
        stop
      end if
    end if

!   check physq

    if ( problem%probnr /= sysvector%probnr ) then
      write(*,'(/a/2(a,i0/))') &
        'Error in fill_sysvector_elem: ',&
        ' problem number in problem = ', problem%probnr, &
        ' whereas problem number in sysvector = ',  sysvector%probnr
      stop
    end if

    if ( present(physq) ) then
      if ( any( physq < 1 ) .or. any( physq > problem%nphysq ) ) then
        write(*,'(2(/a)/)') &
          'Error in fill_sysvector_elem: ', &
          ' physical quantities in physq are out of range.'
        stop
      end if
    end if

!   layers

    if ( present(layer) ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in fill_sysvector_elem: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in fill_sysvector_elem: layer out of range.'
        stop
      end if
      wkl = btest(problem%nodlayers,layer-1)
    else if ( problem%numlayers > 0 ) then
      write(*,'(3(/a)/)') &
        'Error in fill_sysvector_elem: ', &
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
        ' fill_sysvector_elem is out of range: ', &
        ' some groups are < 1 or larger than the number of groups ', mesh%nelgrp
      stop
    end if

    do grp = 1, lnelgrp

      elgrp = lgroups(grp)

!     skip inactive groups
      if ( any( elgrp == problem%inactivegroups ) ) cycle

      if ( all( problem%elnumdegfd(elgrp)%a == 0 ) ) then

        write(*,'(/a/a,i0/)') &
          'Error: element group in the heading of fill_sysvector_elem has ', &
          'no degrees of freedom. elgrp is ', elgrp
        stop

      end if

    end do

    if ( present(addvec) ) then
      addvector = addvec
    else
      addvector = .false.
    end if

    if ( .not. addvector ) then

!     clear vector
      sysvector%u = 0

    end if

!   start loop over all elements

    do grp = 1, lnelgrp

      elgrp = lgroups(grp)

!     skip inactive groups
      if ( any( elgrp == problem%inactivegroups ) ) cycle

      if ( present(mcoefficients) ) coeffl = mcoefficients(elgrp)

!     find number of degrees of freedom in an element of this group

      if ( present(physq) ) then
        ndof = sum ( problem%vec_elnumdegfd(elgrp)%a(:,physq) )
      else
        ndof = sum ( problem%elnumdegfd(elgrp)%a )
      end if

!     reserve memory for element vector and positions

      allocate ( elemvec(ndof) )
      allocate ( pos(ndof) )

!     loop over elements in this group

      do elem = 1, mesh%grpnumel(elgrp)

        first = elem == 1
        last  = elem == mesh%grpnumel(elgrp)

!       skip elements not in layer
        if ( present(layer) ) then
          if ( .not. all ( wkl(mesh%topology(elgrp)%a(:,elem)) ) ) cycle
!         avoid leaving allocated memory in elements
          first = .true.; last = .true.
        end if

!       compute element matrix and vector

        call elemsub ( mesh, problem, elgrp, elem, first, last, coeffl, oldvl, &
          elemvec )

!       compute positions in large matrix/vector of element degrees of freedom

        if ( present(physq) ) then
!         physical quantities
          call pos_array ( mesh, problem, elgrp, elem, dof, pos, physq, &
            order, layer )
        else
          call pos_array ( mesh, problem, elgrp, elem, dof, pos, order=order, &
            layer=layer )
        end if

        if ( dof /= ndof ) then
          write(*,'(/a/a/)') &
            'Internal error in fill_sysvector_elem: ', &
            ' dof /= ndof'
          stop
        end if

!       fill element vector in large vector

        do row = 1, dof

          rowg = pos(row) ! global row number

          sysvector%u(rowg) = elemvec(row)

        end do

      end do

      deallocate ( elemvec, pos )

    end do

    deallocate ( wkl )

  end subroutine fill_sysvector_elem


! Bare bones loop over elements

  subroutine loop_over_elements ( mesh, problem, elemsub, coefficients, &
    mcoefficients, oldvectors, elgroup1, elgroup2, groups, layer, &
    layer_in_all_nodes, exclude_single_layer, skipelementfunc )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   this is the element subroutine that must be supplied by the calling routine
    interface
      subroutine elemsub ( mesh, problem, elgrp, elem, first, last, &
        coefficients, oldvectors )
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
        type(oldvectors_t), intent(inout) :: oldvectors
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
!   Contrary to other build-type routines oldvectors is now an inout variable
!   so the the user can change the data in oldvectors on element level.
    type(oldvectors_t), intent(inout), optional :: oldvectors

!   if these are present the loop takes place for element groups
!   elgroup1,...,elgroup2 only. If only elgroup1 is present one group is
!   used only.
    integer, intent(in), optional :: elgroup1, elgroup2

!   if present: the element groups to be used
!   For example groups=(/2,4/) will loop over element groups 2 and 4.
!   Default: all groups
    integer, dimension(:), intent(in), optional :: groups

!   if present the loop takes place within the single layer of degrees of
!   freedom only.
!   If layer_in_all_nodes=.true. (default) only elements where all nodes have
!   degrees in the specified layer are looped.
!   If layer_in_all_nodes=.false. also the elements where only part of the
!   nodes have degrees in the specified layer are looped.
    integer, intent(in), optional :: layer
    logical, intent(in), optional :: layer_in_all_nodes

!   If exclude_single_layer=.true. elements where all nodes have degrees in
!   a single layer only are excluded from the assembly.
!   If exclude_single_layer=.false. (default) all elements are assembled.
    logical, intent(in), optional :: exclude_single_layer


!   this element function must be supplied by the calling routine if
!   only part of the elements must be looped and must be fully
!   under user control.
!   A value of .true. will skip the element given by (elgrp,elem).
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




!   This is a bare bones loop over elements doing nothing at all. Everything
!   must be done by the writer of the element routine.


    logical :: first, last
    logical :: llayer_in_all_nodes, lexclude_single_layer
    logical, allocatable, dimension(:) :: wkl
    integer :: elem, elgrp, grp
    integer :: lgroups(mesh%nelgrp), lnelgrp, i
    integer, allocatable, dimension(:) :: nlw

    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl

    allocate ( wkl(mesh%nnodes) )

    call check ( mesh, 'loop_over_elements' )
    call check ( problem, 'loop_over_elements' )

    if ( present(mcoefficients) ) then
      if ( size(mcoefficients) /= mesh%nelgrp ) then
        write(*,'(/2(a/))') &
          'Error in loop_over_elements: ', &
          '   size of mcoefficients /= number of element groups'
        stop
      end if
    end if

!   layers

    if ( present(layer) ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in loop_over_elements: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in loop_over_elements: layer out of range.'
        stop
      end if
      wkl = btest(problem%nodlayers,layer-1)
    end if

    llayer_in_all_nodes = &
               set_optional ( variable=layer_in_all_nodes, default = .true. )

    lexclude_single_layer = &
               set_optional ( variable=exclude_single_layer, default = .false. )

    if ( present(layer) .and. lexclude_single_layer ) then
      write(*,'(2(/a)/)') &
        'Error in loop_over_elements: ', &
        ' exclude_single_layer=.true. is incompatible with layer present'
      stop
    end if

    if ( lexclude_single_layer ) then
      allocate ( nlw(mesh%nnodes) )
      nlw = count_layers ( problem%nodlayers, problem%numlayers )
    end if

    if ( lexclude_single_layer ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in loop_over_elements: ',&
          ' exclude_single_layer=.true. ', &
          ' is only available if layers have been defined.'
        stop
      end if
    end if

!   initialize local parameters

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
        ' loop_over_elements is out of range: ', &
        ' some groups are < 1 or larger than the number of groups ', mesh%nelgrp
      stop
    end if

    do grp = 1, lnelgrp

      elgrp = lgroups(grp)

!     skip inactive groups
      if ( any( elgrp == problem%inactivegroups ) ) cycle

      if ( all( problem%elnumdegfd(elgrp)%a == 0 ) ) then

        write(*,'(/a/a,i0/)') &
          'Error: element group in the heading of loop_over_elements has ', &
          'no degrees of freedom. elgrp is ', elgrp
        stop

      end if

    end do

!   start loop over all elements

    do grp = 1, lnelgrp

      elgrp = lgroups(grp)

!     skip inactive groups
      if ( any( elgrp == problem%inactivegroups ) ) cycle

      if ( present(mcoefficients) ) coeffl = mcoefficients(elgrp)

!     loop over elements in this group

      do elem = 1, mesh%grpnumel(elgrp)

        first = elem == 1
        last  = elem == mesh%grpnumel(elgrp)

!       skip user specified elements
        if ( present(skipelementfunc) ) then
          if ( skipelementfunc ( mesh, problem, elgrp, elem, &
                                 coeffl, oldvl ) ) cycle
!         avoid leaving allocated memory in elements
          first = .true.; last = .true.
        end if

!       skip elements not in layer
        if ( present(layer) ) then
          if ( .not. element_in_layer() ) cycle
!         avoid leaving allocated memory in elements
          first = .true.; last = .true.
        end if

!       skip elements not in a single layer
        if ( lexclude_single_layer ) then
          if ( element_in_single_layer() ) cycle
!         avoid leaving allocated memory in elements
          first = .true.; last = .true.
        end if

!       call element routine

        call elemsub ( mesh, problem, elgrp, elem, first, last, coeffl, oldvl )

      end do

    end do

    if ( lexclude_single_layer ) then
      deallocate ( nlw )
    end if

    deallocate ( wkl )

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

  end subroutine loop_over_elements


! create vector subscript array in sysvector

  subroutine create_subscript_sysvector ( mesh, problem, subscript, physqarr, &
    layer, degfd, degsfd, unknownpart, essentialpart, points, curves, &
    surfaces, nodesets, elementsets, groups, excludepoints, excludecurves, &
    excludesurfaces, excludenodesets, excludeelementsets, xmin, xmax, &
    fillnodes, includenodesonce )

    type(mesh_t), intent(in)  :: mesh
    type(problem_t), intent(in)  :: problem

!   the vector subscript of the degrees of freedom in a sysvector
    type(subscript_t), intent(inout) :: subscript

!   if physarr is present only physical quantities in this array are included.
!   For example:
!     physarr=(/1,2/)
!   means that subscript will contain the positions of the physical
!   quantities 1 and 2.
    integer, dimension(:), intent(in), optional :: physqarr

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered.
    integer, intent(in), optional :: layer

!   if degfd is present only the degfd^th degree of freedom will be included.
!   Note that if more than one physical quantity is given in physqarr the
!   degfd^th degree of freedom of the physical quantities combined will be
!   included.
!   Example:
!     physarr not present
!     degfd=2
!   means the second degree of freedom will be included.
!   Example:
!     physarr=[2,1]
!     degfd=2
!   means the second degree of freedom of the physical quantities 2 and 1
!   combined (in that order) will be included.
    integer, intent(in), optional :: degfd

!   if degsfd is present only the degrees of freedom in degsfd will be included.
!   Note that if more than one physical quantity is given in physqarr the
!   degrees of freedom of the physical quantities combined will be
!   included.
!   Example:
!     physarr not present
!     degsfd=[2,3]
!   means the second and third degree of freedom will be included.
!   Example:
!     physarr=[2,1]
!     degsfd=[2,3]
!   means the second and third degree of freedom of the physical quantities
!   2 and 1 combined (in that order) will be included.
    integer, intent(in), dimension(:), optional :: degsfd

!   logicals to indicate inclusion of the unknown degrees and the essential
!   degrees of freedom. Defaults: unknownpart=.true., essentialpart=.true.
    logical, intent(in), optional :: unknownpart, essentialpart

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

!   exclude degrees of freedom in nodes of nodesets
!   for example excludenodesets=(/2/) excludes all degrees in nodes
!   of nodeset 2.
    integer, intent(in), dimension(:), optional :: excludenodesets

!   exclude degrees of freedom in nodes of elementsets
!   for example excludeelementsets=(/2/) excludes all degrees in nodes
!   of elementset 2.
!   NOTE: elementsets are required to have nodes added using add_to_mesh.
    integer, intent(in), dimension(:), optional :: excludeelementsets

!   exclude degrees of freedom in nodes if any of the coordinates is less
!   than the ones given by xmin or any of the coordinates is larger than
!   the ones given by xmax.
!   For example xmin=(/-1._dp,-2._dp/) excludes all nodes that have an
!   x-coordinate less than -1 and/or a y-coordinate less than -2.
    real(dp), intent(in), dimension(:), optional :: xmin, xmax

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
!   freedom in the sysvector.

    logical :: unpart, esspart, allnodes, lfillnodes, lincludenodesonce
    logical :: node_with_degrees
    logical, allocatable, dimension(:) :: excludenodes, nodesdone
    integer :: nodenr, nndof, dof, nnodes, node, crv, srf
    integer :: i, nn, nsnodes
    integer, allocatable, dimension(:) :: work, nodes, snodes
    integer :: pos(problem%maxnoddegfd), set

    allocate ( excludenodes(mesh%nnodes), nodesdone(mesh%nnodes) )

    call check ( mesh, 'create_subscript' )
    call check ( problem, 'create_subscript' )

!   subscript already allocated?
    if (allocated(subscript%s)) deallocate(subscript%s)
    if (allocated(subscript%nodes)) deallocate(subscript%nodes)

!   include unknown part?
    unpart = set_optional ( variable=unknownpart, default=.true. )

!   include essential part?
    esspart = set_optional ( variable=essentialpart, default=.true. )

!   include nodes?
    lfillnodes = set_optional ( variable=fillnodes, default=.false. )

!   include nodes once?
    lincludenodesonce = set_optional ( variable=includenodesonce, &
      default=.true. )

!   check groups
    if ( present(groups) ) then
      if ( any ( groups < 1 ) .or. any ( groups > mesh%nelgrp ) ) then
        write(*,'(/a/a/a,i0/)') &
          'Error: parameter groups in the heading of', &
          ' create_subscript_sysvector is ', &
          ' out of range: some groups are < 1 or larger than ', mesh%nelgrp
        stop
      end if
    end if

!   check degfd and degsfd
    if ( present(degfd) .and. present(degsfd) ) then
      write(*,'(2(/a)/)') &
        'Error in create_subscript_sysvector: ',&
        '  Both degfd and degsfd present in heading.'
      stop
    end if

!   determine nnodes

    nnodes = 0

    if ( present(points) ) then
      if ( any ( points < 1 ) .or. any ( points > mesh%npoints ) ) then
        write(*,'(/a/a/a,i0/)') &
          'Error: parameter points in the heading of', &
          ' create_subscript_sysvector is ', &
          ' out of range: some points are < 1 or larger than ', mesh%npoints
        stop
      end if
      nnodes = nnodes + size(points)
    end if

    if ( present(curves) ) then
      if ( any ( curves < 1 ) .or. any ( curves > mesh%ncurves ) ) then
        write(*,'(/a/a/a,i0/)') &
          'Error: parameter curves in the heading of', &
          ' create_subscript_sysvector is ', &
          ' out of range: some curves are < 1 or larger than ', mesh%ncurves
        stop
      end if
      nnodes = nnodes + sum(mesh%curves(curves)%nnodes)
    end if

    if ( present(surfaces) ) then
      if ( any ( surfaces < 1 ) .or. any ( surfaces > mesh%nsurfaces ) ) then
        write(*,'(/a/a/a,i0/)') &
          'Error: parameter surfaces in the heading of', &
          ' create_subscript_sysvector is ', &
          ' out of range: some surfaces are < 1 or larger than ', mesh%nsurfaces
        stop
      end if
      nnodes = nnodes + sum(mesh%surfaces(surfaces)%nnodes)
    end if

    if ( present(nodesets) ) then
      if ( any ( nodesets < 1 ) .or. any ( nodesets > mesh%nnodesets ) ) then
        write(*,'(/a/a/a,i0/)') &
          'Error: argument nodesets in the heading of', &
          ' create_subscript_sysvector is ', &
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
          ' create_subscript_sysvector is ', &
          ' out of range: some elementsets are < 1 or larger than ', &
          mesh%nelementsets
        stop
      end if
      if ( any ( .not. mesh%elementsets(elementsets)%nodes_created ) ) then
        write(*,'(/4(a/))') &
          'Error: argument elementsets in the heading of', &
          ' create_subscript_sysvector is ', &
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

!   create subscript node for node

    allocate ( work(problem%numdegfd) )

    nsnodes = 0
    dof = 0
    nodesdone = .false.

    do node = 1, nnodes

      nodenr = nodes(node)

!     positions

      if ( present(groups) ) then
!       check whether nodal point is in included groups
        if ( .not. node_in_groups ( mesh, groups, nodenr ) ) cycle
      end if

      if ( present(layer) ) then
!       cycle if layer absent in node
        if ( .not. btest(problem%nodlayers(nodenr),layer-1) ) cycle
      end if

      if ( excludenodes(nodenr) ) cycle

      if ( lincludenodesonce .and. nodesdone(nodenr) ) cycle

!     get positions of degrees

      call pos_array_node ( problem, nodenr, nndof, pos, physqarr, layer )

      call fill_work

      if ( lfillnodes .and. node_with_degrees ) then
        nsnodes = nsnodes + 1
        snodes(nsnodes) = nodenr
      end if

      nodesdone(nodenr) = .true.

    end do

    if ( dof == 0 ) then
      write(*,'(/a/)') &
        'Warning in create_subscript_sysvector: subscript is empty '
    end if

    allocate(subscript%s(dof))

    subscript%s = work(1:dof)

    deallocate(nodes,work)

    if ( lfillnodes ) then
      allocate(subscript%nodes(nsnodes))
      subscript%nodes = snodes(1:nsnodes)
      deallocate ( snodes )
    end if

    deallocate ( excludenodes, nodesdone )

  contains

    subroutine fill_work

      integer :: deg, ip

      node_with_degrees = .false.

      if ( present(degfd) ) then

!       one degree

        if ( degfd >=1 .and. degfd <= nndof ) then

!         degfd in valid range

          ip = pos(degfd)

          if ( ip <= problem%numundegfd .and. unpart .or. &
               ip >  problem%numundegfd .and. esspart ) then

!           position in sysvector and sysmatrix (renumbered)
            work(dof+1) = ip

            dof = dof + 1

            node_with_degrees = .true.

          end if

        end if

      else if ( present(degsfd) ) then

!       multiple degrees

        do deg = 1, size(degsfd)

          ip = pos(degsfd(deg))

          if ( ip <= problem%numundegfd .and. unpart .or. &
               ip >  problem%numundegfd .and. esspart ) then

!           position in sysvector and sysmatrix (renumbered)
            work(dof+1) = ip

            dof = dof + 1

            node_with_degrees = .true.

          end if

        end do

      else

!       all degrees

        do deg = 1, nndof

          ip = pos(deg)

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


!   set excludenodes from points, curves and surfaces

    subroutine set_excludenodes

      integer :: crv, srf, node, set

      excludenodes = .false.

      if ( present(excludepoints) ) then
        if ( any(excludepoints <=0) .or. &
             any(excludepoints > mesh%npoints) ) then
          write(*,'(2(/a),i0/)') &
            'Error in create_subscript_sysvector: ',&
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
            'Error in create_subscript_sysvector: ',&
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
            'Error in create_subscript_sysvector: ',&
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
            'Error in create_subscript_sysvector: ',&
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
          write(*,'(2(/a),a,i0/)') &
            'Error in create_subscript_sysvector: ',&
            ' excludeelementsets <=0 or excludeelementsets > ', &
            'number of elementsets = ', &
             mesh%nelementsets
          stop
        end if
        if ( any ( &
             .not. mesh%elementsets(excludeelementsets)%nodes_created ) ) then
          write(*,'(/4(a/))') &
            'Error: argument excludeelementsets in the heading of', &
            ' create_subscript_sysvector is ', &
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
            'Error in create_subscript_sysvector: ',&
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
            'Error in create_subscript_sysvector: ',&
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

  end subroutine create_subscript_sysvector


! delete vector subscript array in sysvector

  subroutine delete_single_subscript ( subscript )

!   the vector subscript of the degrees of freedom in a sysvector
    type(subscript_t), intent(out) :: subscript

  end subroutine delete_single_subscript


! delete vector subscript array in sysvector

  subroutine delete_subscript ( subscript1, subscript2, subscript3, &
    subscript4, subscript5 )

!   the vector subscript of the degrees of freedom in a sysvector
    type(subscript_t), intent(inout) :: subscript1
    type(subscript_t), intent(inout), optional :: subscript2, subscript3, &
      subscript4, subscript5

    call delete_single_subscript(subscript1)
    if ( present(subscript2) ) call delete_single_subscript(subscript2)
    if ( present(subscript3) ) call delete_single_subscript(subscript3)
    if ( present(subscript4) ) call delete_single_subscript(subscript4)
    if ( present(subscript5) ) call delete_single_subscript(subscript5)

  end subroutine delete_subscript


! clear some rows in sysvector

  subroutine clear_rows_sysvector ( sysvector, rows )

!   The system vector
!   On output the specfied rows (entries) are set to zero
    type(sysvector_t), intent(inout) :: sysvector

!   array of the degrees of freedom numbers for rows to be set to zero
    integer, dimension(:), intent(in) :: rows

    if ( .not. sysvector%created ) then
      write(*,'(/a/)') &
        'Error in clear_rows_sysvector:  sysvector has not been created '
      stop
    end if

    sysvector%u(rows) = 0

  end subroutine clear_rows_sysvector


! clear some rows in multiple (array) sysvector

  subroutine clear_rows_msysvector ( msysvector, rows )

!   The array of system vectors
!   On output the specfied rows (entries) are set to zero
    type(sysvector_t), dimension(:), intent(inout) :: msysvector

!   array of the degrees of freedom numbers for rows to be set to zero
    integer, dimension(:), intent(in) :: rows

    integer :: i

    do i = 1, size(msysvector)
      call clear_rows_sysvector ( msysvector(i), rows )
    end do

  end subroutine clear_rows_msysvector


! clear some rows in multiple (matrix) sysvector

  subroutine clear_rows_m2sysvector ( m2sysvector, rows )

!   The matrix of system vectors
!   On output the specfied rows (entries) are set to zero
    type(sysvector_t), dimension(:,:), intent(inout) :: m2sysvector

!   array of the degrees of freedom numbers for rows to be set to zero
    integer, dimension(:), intent(in) :: rows

    integer :: i, j

    do i = 1, size(m2sysvector,1)
      do j = 1, size(m2sysvector,2)
        call clear_rows_sysvector ( m2sysvector(i,j), rows )
      end do
    end do

  end subroutine clear_rows_m2sysvector


! clear some rows in multiple (3D array) sysvector

  subroutine clear_rows_m3sysvector ( m3sysvector, rows )

!   The matrix of system vectors
!   On output the specfied rows (entries) are set to zero
    type(sysvector_t), dimension(:,:,:), intent(inout) :: m3sysvector

!   array of the degrees of freedom numbers for rows to be set to zero
    integer, dimension(:), intent(in) :: rows

    integer :: i, j, k

    do i = 1, size(m3sysvector,1)
      do j = 1, size(m3sysvector,2)
        do k = 1, size(m3sysvector,3)
          call clear_rows_sysvector ( m3sysvector(i,j,k), rows )
        end do
      end do
    end do

  end subroutine clear_rows_m3sysvector


! transform system vectors from transformed degrees to global degrees

  subroutine transform_sysvector_to_global ( problem, sysvector, &
    sysvector2, sysvector3, sysvector4, sysvector5 )

    type(problem_t), intent(in) :: problem
    type(sysvector_t), intent(inout) :: sysvector
    type(sysvector_t), intent(inout), optional :: sysvector2, sysvector3, &
      sysvector4, sysvector5

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in transform_sysvector_to_global:', &
        ' transformation matrix has not been build'
      stop
    end if

    if ( problem%numtransdegfd > 0 ) then
      call do_one ( sysvector )
      call do_one ( sysvector2 )
      call do_one ( sysvector3 )
      call do_one ( sysvector4 )
      call do_one ( sysvector5 )
    end if

  contains

    subroutine do_one ( sysvector )

      type(sysvector_t), intent(inout), optional :: sysvector

      if ( present(sysvector) ) then
        sysvector%u = smatvec ( problem%Amat, sysvector%u )
      end if

    end subroutine do_one

  end subroutine transform_sysvector_to_global


! transform vector of sysvectors from transformed degrees to the global system

  subroutine transform_msysvector_to_global ( problem, msysvector, &
    msysvector2, msysvector3, msysvector4, msysvector5 )

    type(problem_t), intent(in) :: problem
    type(sysvector_t), dimension(:), intent(inout) :: msysvector
    type(sysvector_t), dimension(:), intent(inout), optional :: msysvector2, &
      msysvector3, msysvector4, msysvector5

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in transform_msysvector_to_global:', &
        ' transformation matrix has not been build'
      stop
    end if

    if ( problem%numtransdegfd > 0 ) then
      call do_one ( msysvector )
      call do_one ( msysvector2 )
      call do_one ( msysvector3 )
      call do_one ( msysvector4 )
      call do_one ( msysvector5 )
    end if

  contains

    subroutine do_one ( msysvector )

      type(sysvector_t), dimension(:), intent(inout), optional :: msysvector

      integer :: i

      if ( present(msysvector) ) then
        do i = 1, size(msysvector)
          msysvector(i)%u = smatvec ( problem%Amat, msysvector(i)%u )
        end do
      end if

    end subroutine do_one

  end subroutine transform_msysvector_to_global


! transform matrix of sysvectors from transformed degrees to the global system

  subroutine transform_m2sysvector_to_global ( problem, m2sysvector, &
    m2sysvector2, m2sysvector3, m2sysvector4, m2sysvector5  )

    type(problem_t), intent(in) :: problem
    type(sysvector_t), dimension(:,:), intent(inout) :: m2sysvector
    type(sysvector_t), dimension(:,:), intent(inout), optional :: &
      m2sysvector2, m2sysvector3, m2sysvector4, m2sysvector5

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in transform_m2sysvector_to_global:', &
        ' transformation matrix has not been build'
      stop
    end if

    if ( problem%numtransdegfd > 0 ) then
      call do_one ( m2sysvector )
      call do_one ( m2sysvector2 )
      call do_one ( m2sysvector3 )
      call do_one ( m2sysvector4 )
      call do_one ( m2sysvector5 )
    end if

  contains

    subroutine do_one ( m2sysvector )

      type(sysvector_t), dimension(:,:), intent(inout), optional :: m2sysvector

      integer :: i, j

      if ( present(m2sysvector) ) then
        do i = 1, size(m2sysvector,1)
          do j = 1, size(m2sysvector,2)
            m2sysvector(i,j)%u = smatvec ( problem%Amat, m2sysvector(i,j)%u )
          end do
        end do
      end if

    end subroutine do_one

  end subroutine transform_m2sysvector_to_global


! transform 3D array of sysvectors from transformed degrees to the global system

  subroutine transform_m3sysvector_to_global ( problem, m3sysvector, &
    m3sysvector2, m3sysvector3, m3sysvector4, m3sysvector5 )

    type(problem_t), intent(in) :: problem
    type(sysvector_t), dimension(:,:,:), intent(inout) :: m3sysvector
    type(sysvector_t), dimension(:,:,:), intent(inout), optional :: &
      m3sysvector2, m3sysvector3, m3sysvector4, m3sysvector5

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in transform_m3sysvector_to_global:', &
        ' transformation matrix has not been build'
      stop
    end if

    if ( problem%numtransdegfd > 0 ) then
      call do_one ( m3sysvector )
      call do_one ( m3sysvector2 )
      call do_one ( m3sysvector3 )
      call do_one ( m3sysvector4 )
      call do_one ( m3sysvector5 )
    end if

  contains

    subroutine do_one ( m3sysvector )

      type(sysvector_t), dimension(:,:,:), intent(inout), optional :: &
                                                                   m3sysvector
      integer :: i, j, k

      if ( present(m3sysvector) ) then
        do i = 1, size(m3sysvector,1)
          do j = 1, size(m3sysvector,2)
            do k = 1, size(m3sysvector,3)
              m3sysvector(i,j,k)%u = &
                                smatvec ( problem%Amat, m3sysvector(i,j,k)%u )
            end do
          end do
        end do
     end if

    end subroutine do_one

  end subroutine transform_m3sysvector_to_global


! transform the system vector from the global system to the transformed degrees

  subroutine transform_sysvector_to_local ( problem, sysvector, &
    sysvector2, sysvector3, sysvector4, sysvector5 )

    type(problem_t), intent(in) :: problem
    type(sysvector_t), intent(inout) :: sysvector
    type(sysvector_t), intent(inout), optional :: sysvector2, sysvector3, &
      sysvector4, sysvector5

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in transform_sysvector_to_local:', &
        ' transformation matrix has not been build'
      stop
    end if

    if ( .not. problem%orthogonal ) then
      write(*,'(/2(a/))') &
        'Error in transform_sysvector_to_local: ', &
        ' non-orthogonal transformation has not been implemented.'
      stop
    end if

    if ( problem%numtransdegfd > 0 ) then
      call do_one ( sysvector )
      call do_one ( sysvector2 )
      call do_one ( sysvector3 )
      call do_one ( sysvector4 )
      call do_one ( sysvector5 )
    end if

  contains

    subroutine do_one ( sysvector )

      type(sysvector_t), intent(inout), optional :: sysvector

      if ( present(sysvector) ) then
        sysvector%u = stmatvec ( problem%Amat, sysvector%u )
      end if

    end subroutine do_one

  end subroutine transform_sysvector_to_local


! transform vector of sysvectors from the global to the transformed degrees

  subroutine transform_msysvector_to_local ( problem, msysvector, &
    msysvector2, msysvector3, msysvector4, msysvector5 )

    type(problem_t), intent(in) :: problem
    type(sysvector_t), dimension(:), intent(inout) :: msysvector
    type(sysvector_t), dimension(:), intent(inout), optional :: msysvector2, &
      msysvector3, msysvector4, msysvector5

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in transform_msysvector_to_local:', &
        ' transformation matrix has not been build'
      stop
    end if

    if ( .not. problem%orthogonal ) then
      write(*,'(/2(a/))') &
        'Error in transform_msysvector_to_local: ', &
        ' non-orthogonal transformation has not been implemented.'
      stop
    end if

    if ( problem%numtransdegfd > 0 ) then
      call do_one ( msysvector )
      call do_one ( msysvector2 )
      call do_one ( msysvector3 )
      call do_one ( msysvector4 )
      call do_one ( msysvector5 )
    end if

  contains

    subroutine do_one ( msysvector )

      type(sysvector_t), dimension(:), intent(inout), optional :: msysvector

      integer :: i

      if ( present(msysvector) ) then
        do i = 1, size(msysvector)
          msysvector(i)%u = stmatvec ( problem%Amat, msysvector(i)%u )
        end do
      end if

    end subroutine do_one

  end subroutine transform_msysvector_to_local


! transform matrix of sysvectors from the global to the transformed degrees

  subroutine transform_m2sysvector_to_local ( problem, m2sysvector, &
    m2sysvector2, m2sysvector3, m2sysvector4, m2sysvector5  )

    type(problem_t), intent(in) :: problem
    type(sysvector_t), dimension(:,:), intent(inout) :: m2sysvector
    type(sysvector_t), dimension(:,:), intent(inout), optional :: &
      m2sysvector2, m2sysvector3, m2sysvector4, m2sysvector5

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in transform_m2sysvector_to_local:', &
        ' transformation matrix has not been build'
      stop
    end if

    if ( .not. problem%orthogonal ) then
      write(*,'(/2(a/))') &
        'Error in transform_m2sysvector_to_local: ', &
        ' non-orthogonal transformation has not been implemented.'
      stop
    end if

    if ( problem%numtransdegfd > 0 ) then
      call do_one ( m2sysvector )
      call do_one ( m2sysvector2 )
      call do_one ( m2sysvector3 )
      call do_one ( m2sysvector4 )
      call do_one ( m2sysvector5 )
    end if

  contains

    subroutine do_one ( m2sysvector )

      type(sysvector_t), dimension(:,:), intent(inout), optional :: m2sysvector

      integer :: i, j

      if ( present(m2sysvector) ) then
        do i = 1, size(m2sysvector,1)
          do j = 1, size(m2sysvector,2)
            m2sysvector(i,j)%u = stmatvec ( problem%Amat, m2sysvector(i,j)%u )
          end do
        end do
      end if

    end subroutine do_one

  end subroutine transform_m2sysvector_to_local


! transform 3D array of sysvectors from the global to the transformed degrees

  subroutine transform_m3sysvector_to_local ( problem, m3sysvector, &
    m3sysvector2, m3sysvector3, m3sysvector4, m3sysvector5 )

    type(problem_t), intent(in) :: problem
    type(sysvector_t), dimension(:,:,:), intent(inout) :: m3sysvector
    type(sysvector_t), dimension(:,:,:), intent(inout), optional :: &
      m3sysvector2, m3sysvector3, m3sysvector4, m3sysvector5

    if ( problem%numtransformations > 0 .and. .not. problem%buildAmat ) then
      write(*,'(/2a/)') &
        'Error in transform_m3sysvector_to_local:', &
        ' transformation matrix has not been build'
      stop
    end if

    if ( .not. problem%orthogonal ) then
      write(*,'(/2(a/))') &
        'Error in transform_m3sysvector_to_local: ', &
        ' non-orthogonal transformation has not been implemented.'
      stop
    end if

    if ( problem%numtransdegfd > 0 ) then
      call do_one ( m3sysvector )
      call do_one ( m3sysvector2 )
      call do_one ( m3sysvector3 )
      call do_one ( m3sysvector4 )
      call do_one ( m3sysvector5 )
    end if

  contains

    subroutine do_one ( m3sysvector )

      type(sysvector_t), dimension(:,:,:), intent(inout), optional :: &
                                                                   m3sysvector
      integer :: i, j, k

      if ( present(m3sysvector) ) then
        do i = 1, size(m3sysvector,1)
          do j = 1, size(m3sysvector,2)
            do k = 1, size(m3sysvector,3)
              m3sysvector(i,j,k)%u = &
                                stmatvec ( problem%Amat, m3sysvector(i,j,k)%u )
            end do
          end do
        end do
     end if

    end subroutine do_one

  end subroutine transform_m3sysvector_to_local

end module system_vector_m
