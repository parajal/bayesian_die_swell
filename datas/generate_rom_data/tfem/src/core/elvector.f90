
! Copyright (C) 2007-2019 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Routines, for the elvectors (vector defined per element)

module elvector_m

  use kind_defs_m
  use mesh_m
  use problem_defs_m
  use system_defs_m
  use vector_defs_m
  use elvector_defs_m
  use element_defs_m
  use set_optional_m

  implicit none


! interface to support legacy interfaces for creating elvectors

  interface create_elvector
    module procedure create_elvector, create_elvector_with_problem
  end interface create_elvector

! interface for generic create subroutine (+ support for legacy interface)

  interface create
    module procedure create_elvector, create_elvector_with_problem
  end interface create

! interface for generic delete subroutine

  interface delete
    module procedure delete_elvector
  end interface delete

! interface to support legacy interfaces for putting elvectors

  interface put_elvector
    module procedure put_elvector, put_elvector_with_problem
  end interface put_elvector

! interface to support legacy interfaces for getting elvectors

  interface get_elvector
    module procedure get_elvector, get_elvector_with_problem
  end interface get_elvector

! interface for generic copy subroutine

  interface copy
    module procedure copy_elvector
  end interface copy

contains


! Create vector

  subroutine create_elvector ( mesh, elvector, nint1d, nreal1d, nreal2d, &
    nreal3d )

    type(mesh_t), intent(in) :: mesh
    type(elvector_t), intent(inout) :: elvector

!   elvector will contain nint1d one-dimensional integer vectors per element
!   NOTE: this number is used for all element groups
    integer, intent(in), optional :: nint1d

!   elvector will contain nreal1d one-dimensional real vectors per element
!   NOTE: this number is used for all element groups
    integer, intent(in), optional :: nreal1d

!   elvector will contain nreal2d two-dimensional real vectors per element
!   NOTE: this number is used for all element groups
    integer, intent(in), optional :: nreal2d

!   elvector will contain nreal3d three-dimensional real vectors per element
!   NOTE: this number is used for all element groups
    integer, intent(in), optional :: nreal3d


    integer :: lnint1d, lnreal1d, lnreal2d, lnreal3d
    integer :: elgrp


!   test

    call check ( mesh, 'create_elvector' )

    if ( elvector%created ) then
      write(*,'(/a/)') &
        'Error in create_elvector: elvector has already been created '
      stop
    end if

    lnint1d = set_optional ( variable=nint1d, default=0 )
    lnreal1d = set_optional ( variable=nreal1d, default=0 )
    lnreal2d = set_optional ( variable=nreal2d, default=0 )
    lnreal3d = set_optional ( variable=nreal3d, default=0 )

    if ( lnint1d < 0 ) then
      write(*,'(/a/)') &
        'Error in create_elvector: nint1d is negative'
      stop
    end if
    if ( lnreal1d < 0 ) then
      write(*,'(/a/)') &
        'Error in create_elvector: nreal1d is negative'
      stop
    end if
    if ( lnreal2d < 0 ) then
      write(*,'(/a/)') &
        'Error in create_elvector: nreal2d is negative'
      stop
    end if
    if ( lnreal3d < 0 ) then
      write(*,'(/a/)') &
        'Error in create_elvector: nreal3d is negative'
      stop
    end if

    if ( lnint1d == 0 .and. lnreal1d == 0 .and. lnreal2d == 0 &
      .and. lnreal3d == 0 ) then
      write(*,'(/a/)') &
        'Warning in create_elvector: empty vector'
    end if

    allocate(elvector%g(mesh%nelgrp))
    do elgrp = 1, mesh%nelgrp
      allocate(elvector%g(elgrp)%i1(mesh%grpnumel(elgrp),lnint1d))
      allocate(elvector%g(elgrp)%r1(mesh%grpnumel(elgrp),lnreal1d))
      allocate(elvector%g(elgrp)%r2(mesh%grpnumel(elgrp),lnreal2d))
      allocate(elvector%g(elgrp)%r3(mesh%grpnumel(elgrp),lnreal3d))
    end do

    elvector%created = .true.

  end subroutine create_elvector


! Create elvector with problem argument (legacy interface)

  subroutine create_elvector_with_problem ( mesh, problem, elvector, nint1d, &
    nreal1d, nreal2d, nreal3d )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(elvector_t), intent(inout) :: elvector
    integer, intent(in), optional :: nint1d
    integer, intent(in), optional :: nreal1d
    integer, intent(in), optional :: nreal2d
    integer, intent(in), optional :: nreal3d

    logical, save :: firstcall = .true.

    if ( firstcall ) write(*,'(/3(a/)/)') &
      'Warning in create_elvector_with_problem:', &
      ' Elvector has been made problem agnostic. Remove problem from the ', &
      ' arguments of create_elvector to get rid of this message.'

    firstcall = .false.

    call create_elvector ( mesh, elvector, nint1d, nreal1d, nreal2d, nreal3d )

  end subroutine create_elvector_with_problem


! Copy element vector to another element vector

  subroutine copy_elvector ( elvector1, elvector2 )

    type(elvector_t), intent(inout) :: elvector1, elvector2

!   This routine performs a real data copy of one elvector to the other using
!   the assignment statement:
!
!      elvector2 = elvector1
!

    if ( .not. elvector1%created ) then
      write(*,'(/a/)') &
        'Error in copy_elvector: elvector1 has not been created '
      stop
    end if

    elvector2 = elvector1

  end subroutine copy_elvector


! Delete single elvector

  subroutine delete_single_elvector ( elvector, nr, deletestructure )

    type(elvector_t), intent(inout) :: elvector
    integer, intent(in) :: nr
    logical, intent(in), optional :: deletestructure

    integer :: elgrp, elem, i
    logical :: ldelstructure

    if ( .not. elvector%created ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error in delete_elvector:', &
        ' Elvector has not been created and cannot be deleted ', &
        ' Vector number in heading = ', nr
      stop
    end if

    ldelstructure = set_optional ( variable=deletestructure, default=.true. )

    do elgrp = 1, size(elvector%g)

!     delete data if allocated

      do elem = 1, size(elvector%g(elgrp)%i1,1)
        do i = 1, size(elvector%g(elgrp)%i1,2)
          if ( allocated(elvector%g(elgrp)%i1(elem,i)%a) ) then
            deallocate ( elvector%g(elgrp)%i1(elem,i)%a )
          end if
        end do
        do i = 1, size(elvector%g(elgrp)%r1,2)
          if ( allocated(elvector%g(elgrp)%r1(elem,i)%a) ) then
            deallocate ( elvector%g(elgrp)%r1(elem,i)%a )
          end if
        end do
        do i = 1, size(elvector%g(elgrp)%r2,2)
          if ( allocated(elvector%g(elgrp)%r2(elem,i)%a) ) then
            deallocate ( elvector%g(elgrp)%r2(elem,i)%a )
          end if
        end do
        do i = 1, size(elvector%g(elgrp)%r3,2)
          if ( allocated(elvector%g(elgrp)%r3(elem,i)%a) ) then
            deallocate ( elvector%g(elgrp)%r3(elem,i)%a )
          end if
        end do
      end do

      if ( ldelstructure ) then

!       delete structure

        deallocate(elvector%g(elgrp)%i1)
        deallocate(elvector%g(elgrp)%r1)
        deallocate(elvector%g(elgrp)%r2)
        deallocate(elvector%g(elgrp)%r3)

      end if

    end do

    if ( ldelstructure ) then

!     delete structure

      deallocate(elvector%g)

      elvector%created = .false.

    end if

  end subroutine delete_single_elvector


! Delete vector

  subroutine delete_elvector ( vector1, vector2, vector3, vector4, vector5, &
    deletestructure )

    type(elvector_t), intent(inout) :: vector1
    type(elvector_t), intent(inout), optional :: vector2, vector3, vector4, &
      vector5

!   logical to indicate if the complete structure needs to be deleted
!   default=.true.
!   if deletestructure=.false. data will be removed and the structure of
!   the vector will be left intact. This has the same effect as a
!   delete + a create.
    logical, intent(in), optional :: deletestructure

    call delete_single_elvector(vector1,1,deletestructure)
    if ( present(vector2) ) &
      call delete_single_elvector(vector2,2,deletestructure)
    if ( present(vector3) ) &
      call delete_single_elvector(vector3,3,deletestructure)
    if ( present(vector4) ) &
      call delete_single_elvector(vector4,4,deletestructure)
    if ( present(vector5) ) &
      call delete_single_elvector(vector5,5,deletestructure)

  end subroutine delete_elvector


! Get element values of elvector for the element given by (elgrp,elem) and array
! number nr.

  subroutine get_elvector ( mesh, elvector, elgrp, elem, nr, i1, r1, r2, r3 )

    type(mesh_t), intent(in) :: mesh
    type(elvector_t), intent(in) :: elvector

!   element group, element number (within group)
    integer, intent(in) :: elgrp, elem

!   array number
!   default=1
    integer, intent(in), optional :: nr

!   the element degrees of freedom (integer 1D array), if present
    integer, intent(out), dimension(:), optional :: i1

!   the element degrees of freedom (real 1D array), if present
    real(dp), intent(out), dimension(:), optional :: r1

!   the element degrees of freedom (real 2D array), if present
    real(dp), intent(out), dimension(:,:), optional :: r2

!   the element degrees of freedom (real 3D array), if present
    real(dp), intent(out), dimension(:,:,:), optional :: r3


    integer :: num, lnr


!   testing

    if ( .not. elvector%created ) then
      write(*,'(/a/)') &
        'Error in get_elvector: elvector not created.'
      stop
    end if

    num = 0
    if ( present(i1) ) num = num + 1
    if ( present(r1) ) num = num + 1
    if ( present(r2) ) num = num + 1
    if ( present(r3) ) num = num + 1

    if ( num > 1 ) then
      write(*,'(/a/a/)') &
        'Error in get_elvector: ', &
          ' Only one of i1, r1, r2, r3 can be present'
      stop
    end if

    lnr = set_optional ( variable=nr, default=1 )

!   i1 array

    if ( present(i1) ) then

      if ( lnr < 1 .or. lnr > size(elvector%g(elgrp)%i1,2) ) then
        write(*,'(/a/)') &
          'Error in get_elvector: nr out of range for array i1 '
        stop
      end if

      if ( .not. allocated(elvector%g(elgrp)%i1(elem,lnr)%a) ) then
        write(*,'(/a/)') &
          'Error in get_elvector: i1 array in elvector not allocated '
        stop
      end if

      if ( any ( shape(i1) /= &
                 shape(elvector%g(elgrp)%i1(elem,lnr)%a) ) ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_elvector: ', &
          ' i1 array in heading has wrong dimensions = ', shape(i1), &
          ' whereas in elvector the dimensions are = ', &
                     shape(elvector%g(elgrp)%i1(elem,lnr)%a)
        stop
      end if

      i1 = elvector%g(elgrp)%i1(elem,lnr)%a

    end if

!   r1 array

    if ( present(r1) ) then

      if ( lnr < 1 .or. lnr > size(elvector%g(elgrp)%r1,2) ) then
        write(*,'(/a/)') &
          'Error in get_elvector: nr out of range for array r1 '
        stop
      end if

      if ( .not. allocated(elvector%g(elgrp)%r1(elem,lnr)%a) ) then
        write(*,'(/a/)') &
          'Error in get_elvector: r1 array in elvector not allocated '
        stop
      end if

      if ( any ( shape(r1) /= &
                 shape(elvector%g(elgrp)%r1(elem,lnr)%a) ) ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_elvector: ', &
          ' r1 array in heading has wrong dimensions = ', shape(r1), &
          ' whereas in elvector the dimensions are = ', &
                     shape(elvector%g(elgrp)%r1(elem,lnr)%a)
        stop
      end if

      r1 = elvector%g(elgrp)%r1(elem,lnr)%a

    end if

!   r2 array

    if ( present(r2) ) then

      if ( lnr < 1 .or. lnr > size(elvector%g(elgrp)%r2,2) ) then
        write(*,'(/a/)') &
          'Error in get_elvector: nr out of range for array r2 '
        stop
      end if

      if ( .not. allocated(elvector%g(elgrp)%r2(elem,lnr)%a) ) then
        write(*,'(/a/)') &
          'Error in get_elvector: r2 array in elvector not allocated '
        stop
      end if

      if ( any ( shape(r2) /= &
                 shape(elvector%g(elgrp)%r2(elem,lnr)%a) ) ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_elvector: ', &
          ' r1 array in heading has wrong dimensions = ', shape(r2), &
          ' whereas in elvector the dimensions are = ', &
                     shape(elvector%g(elgrp)%r2(elem,lnr)%a)
        stop
      end if

      r2 = elvector%g(elgrp)%r2(elem,lnr)%a

    end if

!   r3 array

    if ( present(r3) ) then

      if ( lnr < 1 .or. lnr > size(elvector%g(elgrp)%r3,2) ) then
        write(*,'(/a/)') &
          'Error in get_elvector: nr out of range for array r3 '
        stop
      end if

      if ( .not. allocated(elvector%g(elgrp)%r3(elem,lnr)%a) ) then
        write(*,'(/a/)') &
          'Error in get_elvector: r3 array in elvector not allocated '
        stop
      end if

      if ( any ( shape(r3) /= &
                 shape(elvector%g(elgrp)%r3(elem,lnr)%a) ) ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error in get_elvector: ', &
          ' r1 array in heading has wrong dimensions = ', shape(r3), &
          ' whereas in elvector the dimensions are = ', &
                     shape(elvector%g(elgrp)%r3(elem,lnr)%a)
        stop
      end if

      r3 = elvector%g(elgrp)%r3(elem,lnr)%a

    end if

  end subroutine get_elvector


! Legacy interface to get_elvector

  subroutine get_elvector_with_problem ( mesh, problem, elvector, elgrp, elem, &
    nr, i1, r1, r2, r3 )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(elvector_t), intent(in) :: elvector
    integer, intent(in) :: elgrp, elem
    integer, intent(in), optional :: nr
    integer, intent(out), dimension(:), optional :: i1
    real(dp), intent(out), dimension(:), optional :: r1
    real(dp), intent(out), dimension(:,:), optional :: r2
    real(dp), intent(out), dimension(:,:,:), optional :: r3

    logical, save :: firstcall = .true.

    if ( firstcall ) write(*,'(/3(a/)/)') &
      'Warning in get_elvector_with_problem:', &
      ' Elvector has been made problem agnostic. Remove problem from the ', &
      ' arguments of get_elvector to get rid of this message.'

    firstcall = .false.

    call get_elvector ( mesh, elvector, elgrp, elem, nr, i1, r1, r2, r3 )

  end subroutine get_elvector_with_problem


! Put element values of elvector for the element given by (elgrp,elem) and array
! number nr.

  subroutine put_elvector ( mesh, elvector, elgrp, elem, nr, i1, r1, r2, r3 )

    type(mesh_t), intent(in) :: mesh
    type(elvector_t), intent(inout) :: elvector

!   element group, element number (within group)
    integer, intent(in) :: elgrp, elem

!   array number
!   default=1
    integer, intent(in), optional :: nr

!   the element degrees of freedom (integer 1D array), if present
    integer, intent(in), dimension(:), optional :: i1

!   the element degrees of freedom (real 1D array), if present
    real(dp), intent(in), dimension(:), optional :: r1

!   the element degrees of freedom (real 2D array), if present
    real(dp), intent(in), dimension(:,:), optional :: r2

!   the element degrees of freedom (real 3D array), if present
    real(dp), intent(in), dimension(:,:,:), optional :: r3


    integer :: num, lnr


!   testing

    if ( .not. elvector%created ) then
      write(*,'(/a/)') &
        'Error in put_elvector: elvector not created.'
      stop
    end if

    num = 0
    if ( present(i1) ) num = num + 1
    if ( present(r1) ) num = num + 1
    if ( present(r2) ) num = num + 1
    if ( present(r3) ) num = num + 1

    if ( num > 1 ) then
      write(*,'(/a/a/)') &
        'Error in put_elvector: ', &
          ' Only one of i1, r1, r2, r3 can be present'
      stop
    end if

    lnr = set_optional ( variable=nr, default=1 )

!   i1 array

    if ( present(i1) ) then

      if ( lnr < 1 .or. lnr > size(elvector%g(elgrp)%i1,2) ) then
        write(*,'(/a/)') &
          'Error in put_elvector: nr out of range for array i1 '
        stop
      end if

      if ( allocated(elvector%g(elgrp)%i1(elem,lnr)%a) ) then
        if ( any ( shape(i1) /= &
                   shape(elvector%g(elgrp)%i1(elem,lnr)%a) ) ) then
!         reallocate memory: shapes differ
          deallocate(elvector%g(elgrp)%i1(elem,lnr)%a)
          allocate(elvector%g(elgrp)%i1(elem,lnr)%a(size(i1)))
        end if
      else
        allocate(elvector%g(elgrp)%i1(elem,lnr)%a(size(i1)))
      end if

      elvector%g(elgrp)%i1(elem,lnr)%a = i1

    end if

!   r1 array

    if ( present(r1) ) then

      if ( lnr < 1 .or. lnr > size(elvector%g(elgrp)%r1,2) ) then
        write(*,'(/a/)') &
          'Error in put_elvector: nr out of range for array r1 '
        stop
      end if

      if ( allocated(elvector%g(elgrp)%r1(elem,lnr)%a) ) then
        if ( any ( shape(r1) /= &
                   shape(elvector%g(elgrp)%r1(elem,lnr)%a) ) ) then
!         reallocate memory: shapes differ
          deallocate(elvector%g(elgrp)%r1(elem,lnr)%a)
          allocate(elvector%g(elgrp)%r1(elem,lnr)%a(size(r1)))
        end if
      else
        allocate(elvector%g(elgrp)%r1(elem,lnr)%a(size(r1)))
      end if

      elvector%g(elgrp)%r1(elem,lnr)%a = r1

    end if

!   r2 array

    if ( present(r2) ) then

      if ( lnr < 1 .or. lnr > size(elvector%g(elgrp)%r2,2) ) then
        write(*,'(/a/)') &
          'Error in put_elvector: nr out of range for array r2 '
        stop
      end if

      if ( allocated(elvector%g(elgrp)%r2(elem,lnr)%a) ) then
        if ( any ( shape(r2) /= &
                   shape(elvector%g(elgrp)%r2(elem,lnr)%a) ) ) then
!         reallocate memory: shapes differ
          deallocate(elvector%g(elgrp)%r2(elem,lnr)%a)
          allocate(elvector%g(elgrp)%r2(elem,lnr)%a(size(r2,1),size(r2,2)))
        end if
      else
        allocate(elvector%g(elgrp)%r2(elem,lnr)%a(size(r2,1),size(r2,2)))
      end if

      elvector%g(elgrp)%r2(elem,lnr)%a = r2

    end if

!   r3 array

    if ( present(r3) ) then

      if ( lnr < 1 .or. lnr > size(elvector%g(elgrp)%r3,2) ) then
        write(*,'(/a/)') &
          'Error in put_elvector: nr out of range for array r3 '
        stop
      end if

      if ( allocated(elvector%g(elgrp)%r3(elem,lnr)%a) ) then
        if ( any ( shape(r3) /= &
                   shape(elvector%g(elgrp)%r3(elem,lnr)%a) ) ) then
!         reallocate memory: shapes differ
          deallocate(elvector%g(elgrp)%r3(elem,lnr)%a)
          allocate(elvector%g(elgrp)%r3(elem,lnr)&
                        &%a(size(r3,1),size(r3,2),size(r3,3)))
        end if
      else
        allocate(elvector%g(elgrp)%r3(elem,lnr)&
                      &%a(size(r3,1),size(r3,2),size(r3,3)))
      end if

      elvector%g(elgrp)%r3(elem,lnr)%a = r3

    end if

  end subroutine put_elvector


! Legacy interface for put_elvector

  subroutine put_elvector_with_problem ( mesh, problem, elvector, elgrp, elem,&
    nr, i1, r1, r2, r3 )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    type(elvector_t), intent(inout) :: elvector
    integer, intent(in) :: elgrp, elem
    integer, intent(in), optional :: nr
    integer, intent(in), dimension(:), optional :: i1
    real(dp), intent(in), dimension(:), optional :: r1
    real(dp), intent(in), dimension(:,:), optional :: r2
    real(dp), intent(in), dimension(:,:,:), optional :: r3

    logical, save :: firstcall = .true.

    if ( firstcall ) write(*,'(/3(a/)/)') &
      'Warning in put_elvector_with_problem:', &
      ' Elvector has been made problem agnostic. Remove problem from the ', &
      ' arguments of put_elvector to get rid of this message.'

    call put_elvector ( mesh, elvector, elgrp, elem, nr, i1, r1, r2, r3 )

    firstcall = .false.

  end subroutine put_elvector_with_problem

end module elvector_m
