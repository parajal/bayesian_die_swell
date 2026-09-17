
! Copyright (C) 2005-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types and associated routines for the element subroutines

module element_defs_m

  use kind_defs_m
  use array_defs_m
  use mesh_m, only: mesh_t, me_array_1d_t
  use problem_defs_m, only: problem_t
  use system_defs_m, only: sysvector_t
  use vector_defs_m, only: vector_t
  use elvector_defs_m, only: elvector_t
  use eltree_basic_m, only: eltree_p
  use set_optional_m

  implicit none


! encapsulated function pointers for use in coefficients

  type func_p
    procedure (dum_func), pointer, nopass :: p => null()
  end type func_p

  type vfunc_p
    procedure (dum_vfunc), pointer, nopass :: p => null()
  end type vfunc_p

  type tfunc_p
    procedure (dum_tfunc), pointer, nopass :: p => null()
  end type tfunc_p


! type definition of coefficients
!
  type coefficients_t

!   has the structure been created?
    logical :: created = .false.

!   integer coefficients array
    integer, allocatable, dimension(:) :: i

!   real coefficients array
    real(dp), allocatable, dimension(:) :: r

!   array of integer coefficient arrays
    type(int_array_1d_t), allocatable, dimension(:) :: ia

!   array of real coefficient arrays
    type(real_array_1d_t), allocatable, dimension(:) :: ra

!   array of 2D integer coefficient arrays
    type(int_array_2d_t), allocatable, dimension(:) :: ia2

!   array of 2D real coefficient arrays
    type(real_array_2d_t), allocatable, dimension(:) :: ra2

!   procedure pointer for func
    procedure (dum_func), pointer, nopass :: func => null()

!   procedure pointer for vfunc
    procedure (dum_vfunc), pointer, nopass :: vfunc => null()

!   procedure pointer for tfunc
    procedure (dum_tfunc), pointer, nopass :: tfunc => null()

!   array of func procedure pointers
    type(func_p), allocatable, dimension(:) :: func1

!   array of vfunc procedure pointers
    type(vfunc_p), allocatable, dimension(:) :: vfunc1

!   array of tfunc procedure pointers
    type(tfunc_p), allocatable, dimension(:) :: tfunc1

!   procedure pointers for variable Gauss points
    procedure (dum_set_ninti_user), pointer, nopass :: set_ninti_user  => null()
    procedure (dum_set_Gauss_integration_user), pointer, nopass :: &
      set_Gauss_integration_user => null()

!   procedure pointers for DG
    procedure (dum_cinflow), pointer, nopass :: cinflow  => null()
    procedure (dum_qinflow), pointer, nopass :: qinflow  => null()
    procedure (dum_brownian_vector), pointer, nopass :: &
      brownian_vector  => null()

!   workaround for bug in ifort version 12.0.1.107 (and below).
    integer :: unused_integer = 0

  end type coefficients_t


! pointers for use in oldvectors

  type sysvector_p
    type(sysvector_t), pointer :: p => null()
  end type sysvector_p

  type sysvector_p1
    type(sysvector_t), pointer, dimension(:) :: p => null()
  end type sysvector_p1

  type sysvector_p2
    type(sysvector_t), pointer, dimension(:,:) :: p => null()
  end type sysvector_p2

  type sysvector_p3
    type(sysvector_t), pointer, dimension(:,:,:) :: p => null()
  end type sysvector_p3

  type vector_p
    type(vector_t), pointer :: p => null()
  end type vector_p

  type vector_p1
    type(vector_t), pointer, dimension(:) :: p => null()
  end type vector_p1

  type elvector_p
    type(elvector_t), pointer :: p => null()
  end type elvector_p

  type problem_p
    type(problem_t), pointer :: p => null()
  end type problem_p

  type mesh_p
    type(mesh_t), pointer :: p => null()
  end type mesh_p

  type eltree_array_p
    type(eltree_p), dimension(:,:), pointer :: p => null()
  end type eltree_array_p

  type melements_p
    type(me_array_1d_t), pointer :: p => null()
  end type melements_p

! type definition of oldvectors
!
  type oldvectors_t

!   has the structure been created?
    logical :: created = .false.

!   system vectors
    type(sysvector_p), allocatable, dimension(:) :: s ! sysvectors

    type(sysvector_p1), allocatable, dimension(:) :: s1 ! arrays of sysvectors

    type(sysvector_p2), allocatable, dimension(:) :: s2 !2d arrays of sysvectors

    type(sysvector_p3), allocatable, dimension(:) :: s3 !3d arrays of sysvectors

    type(vector_p), allocatable, dimension(:) :: v ! vectors

    type(vector_p1), allocatable, dimension(:) :: v1 ! arrays of vectors

    type(elvector_p), allocatable, dimension(:) :: e ! elvectors

    type(problem_p), allocatable, dimension(:) :: p  ! problems

    type(mesh_p), allocatable, dimension(:) :: m ! meshes

    type(eltree_array_p), allocatable, dimension(:) :: ea ! eltree arrays

    type(melements_p) :: me ! multilevel elements

  end type oldvectors_t


! type definition of eleminfo
!
  type eleminfo_t

    integer :: elgrp = 0    ! group number
    integer :: elem = 0     ! element number
    integer :: node = 0     ! nodal number
    integer :: point = 0    ! point number
    integer :: curve = 0    ! curve number
    integer :: surface = 0  ! surface number
    integer :: constr = 0   ! constraint number
    integer :: object = 0   ! object number
    integer :: elemo = 0    ! element number in the object
    integer :: nintps = 0   ! number of integration points

!   integration point numbers in object element elemo
    integer, pointer, dimension(:) :: intps => null()

  end type eleminfo_t


! interface for generic create subroutine

  interface create
    module procedure create_coefficients, create_mcoefficients, &
      create_oldvectors
  end interface create


! interface for generic delete subroutine

  interface delete
    module procedure delete_coefficients, delete_mcoefficients, &
                     delete_oldvectors
  end interface delete

! interface for generic check subroutine

  interface check
    module procedure check_coefficients
  end interface check

contains


! construct a structure of type coefficients_t

  subroutine create_coefficients ( coefficients, ncoefi, ncoefr, ncoefia, &
    ncoefra, ncoefia2, ncoefra2 )

    use limits_m, only: MAXFUNCTIONS

    type(coefficients_t), intent(inout) :: coefficients

!   number of coefficients for integer and real arrays
    integer, intent(in) :: ncoefi, ncoefr

!   number of coefficients for arrays of integer and real arrays
    integer, intent(in), dimension(:), optional :: ncoefia, ncoefra

!   number of coefficients for 2D arrays of integer and real arrays
    integer, intent(in), dimension(:), optional :: ncoefia2, ncoefra2

    integer :: i, ni, nr, ni2, nr2

    if ( coefficients%created ) then
      write(*,'(/a/)') &
        'Error in create_coefficients: coefficients has already been created '
      stop
    end if

    allocate(coefficients%i(ncoefi),coefficients%r(ncoefr))

    if ( present(ncoefia) ) then
      ni = size(ncoefia)
      allocate(coefficients%ia(ni))
      do i = 1, ni
        allocate(coefficients%ia(i)%a(ncoefia(i)))
      end do
    else
      ni = 0
      allocate(coefficients%ia(0))
    end if

    if ( present(ncoefra) ) then
      nr = size(ncoefra)
      allocate(coefficients%ra(nr))
      do i = 1, nr
        allocate(coefficients%ra(i)%a(ncoefra(i)))
      end do
    else
      nr = 0
      allocate(coefficients%ra(0))
    end if

    if ( present(ncoefia2) ) then
      ni2 = size(ncoefia2)/2
      allocate(coefficients%ia2(ni2))
      do i = 1, ni2
        allocate(coefficients%ia2(i)%a(ncoefia2(2*i-1),ncoefia2(2*i)))
      end do
    else
      ni2 = 0
      allocate(coefficients%ia2(0))
    end if

    if ( present(ncoefra2) ) then
      nr2 = size(ncoefra2)/2
      allocate(coefficients%ra2(nr2))
      do i = 1, nr2
        allocate(coefficients%ra2(i)%a(ncoefra2(2*i-1),ncoefra2(2*i)))
      end do
    else
      nr2 = 0
      allocate(coefficients%ra2(0))
    end if

    allocate(coefficients%func1(MAXFUNCTIONS))
    allocate(coefficients%vfunc1(MAXFUNCTIONS))
    allocate(coefficients%tfunc1(MAXFUNCTIONS))

!   initialize to zero
    coefficients%i = 0
    coefficients%r = 0
    do i = 1, ni
      coefficients%ia(i)%a = 0
    end do
    do i = 1, nr
      coefficients%ra(i)%a = 0
    end do
    do i = 1, ni2
      coefficients%ia2(i)%a = 0
    end do
    do i = 1, nr2
      coefficients%ra2(i)%a = 0
    end do

!   initialize procedure pointers to dummy routines

    coefficients%func => dum_func
    coefficients%vfunc => dum_vfunc
    coefficients%tfunc => dum_tfunc

    do i = 1, MAXFUNCTIONS
      coefficients%func1(i)%p => dum_func
      coefficients%vfunc1(i)%p => dum_vfunc
      coefficients%tfunc1(i)%p => dum_tfunc
    end do

    coefficients%set_ninti_user => dum_set_ninti_user
    coefficients%set_Gauss_integration_user => dum_set_Gauss_integration_user

    coefficients%cinflow => dum_cinflow
    coefficients%qinflow => dum_qinflow
    coefficients%brownian_vector => dum_brownian_vector

    coefficients%created = .true.

  end subroutine create_coefficients


! construct a array of structures of type coefficients_t

  subroutine create_mcoefficients ( mcoefficients, ncoefi, ncoefr, ncoefia, &
    ncoefra )

    type(coefficients_t), dimension(:), intent(inout) :: mcoefficients

!   number of coefficients for integer and real arrays
    integer, intent(in) :: ncoefi, ncoefr

!   number of coefficients for arrays of integer and real arrays
    integer, intent(in), dimension(:), optional :: ncoefia, ncoefra

    integer :: i

    do i = 1, size(mcoefficients)
      call create_coefficients ( mcoefficients(i), ncoefi, ncoefr, ncoefia, &
        ncoefra )
    end do

  end subroutine create_mcoefficients


! check structure coefficients

  subroutine check_coefficients ( coefficients, name_of_routine, ncoefi, &
    ncoefr, indexarray, minimum, maximum, wncoefi, wncoefr )

    type(coefficients_t), intent(in) :: coefficients
    character(len=*), intent(in) :: name_of_routine

!   the required number of coefficients for integer and real
!   In case wncoefi and/or wncoefr are also present this should be
!   the values before the increase.
    integer, optional, intent(in) :: ncoefi, ncoefr

!   the required minimum and/or maximum for the integer coefficients supplied
!   in indexarray
    integer, dimension(:), optional, intent(in) :: indexarray, minimum, maximum

!   the required number of coefficients for integer and real after they have
!   been increased. This gives an error message if the required number is
!   too low, but also additional useful info for changing the user code. In
!   a future TFEM release the optional argument can be removed again and
!   ncoefi, ncoefr changed to the actual required values.
    integer, optional, intent(in) :: wncoefi, wncoefr


    logical :: error
    integer :: i, j

    error = .false.

    if ( .not. coefficients%created ) then
      write(*,'(3a/)') &
        'Error in ', name_of_routine, ': coefficients have not been created '
      stop
    end if

    if ( present(ncoefi) ) then
      if ( size( coefficients%i ) < ncoefi ) then
        write(*,'(/3a/,2(a,i0/))') &
          'Error ', name_of_routine, ':', &
          ' The number of integer coefficients must be at least ', ncoefi, &
          ' whereas size(coefficients%i) = ', size( coefficients%i )
        error = .true.
      end if
      if ( present(wncoefi) ) then
        if ( size( coefficients%i ) < wncoefi ) then
          write(*,'(/3a/,2(a,i0/),2a,i0,/a/)') &
            'Error ', name_of_routine, ':', &
            ' The number of integer coefficients must be at least ', wncoefi, &
            ' whereas size(coefficients%i) = ', size( coefficients%i ), &
            ' The previous required minimum number of integer coefficients', &
            ' was ', ncoefi, &
            ' Please increase the size of coefficients%i'
          error = .true.
        end if
      end if
    end if

    if ( present(ncoefr) ) then
      if ( size( coefficients%r ) < ncoefr ) then
        write(*,'(/3a/,2(a,i0/))') &
          'Error ', name_of_routine, ':', &
          ' The number of real coefficients must be at least ', ncoefr, &
          ' whereas size(coefficients%r) = ', size( coefficients%r )
        error = .true.
      end if
      if ( present(wncoefr) ) then
        if ( size( coefficients%r ) < wncoefr ) then
          write(*,'(/3a/,2(a,i0/),2a,i0,/a/)') &
            'Error ', name_of_routine, ':', &
            ' The number of real coefficients must be at least ', wncoefr, &
            ' whereas size(coefficients%r) = ', size( coefficients%r ), &
            ' The previous required minimum number of real coefficients', &
            ' was ', ncoefr, &
            ' Please increase the size of coefficients%r'
          error = .true.
        end if
      end if
    end if

    if ( error ) stop

    if ( present(indexarray) ) then

      if ( present(minimum) ) then
        do i = 1, size(indexarray)
          j = indexarray(i)
          if ( coefficients%i(j) < minimum(i) ) then
            write(*,'(/3a/,2(a,i0/))') &
              'Error ', name_of_routine, ':', &
              ' The value of integer coefficient ', j, &
              ' should be at at least ', minimum(i)
            error = .true.
          end if
        end do
      end if
      if ( present(maximum) ) then
        do i = 1, size(indexarray)
          j = indexarray(i)
          if ( coefficients%i(j) > maximum(i) ) then
            write(*,'(/3a/,2(a,i0/))') &
              'Error ', name_of_routine, ':', &
              ' The value of integer coefficient ', j, &
              ' should not be larger than ', maximum(i)
            error = .true.
          end if
        end do
      end if

    end if

    if ( error ) stop

  end subroutine check_coefficients


! get integer coefficient

  function get_coefficient ( coefficients, index, default )

    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: index, default
    integer :: get_coefficient

    if ( coefficients%i(index) == 0 ) then
      get_coefficient = default
    else
      get_coefficient = coefficients%i(index)
    end if

  end function get_coefficient


! delete a (single) structure of type coefficients_t

  subroutine delete_single_coefficients ( coefficients, nr )

    type(coefficients_t), intent(inout) :: coefficients
    integer, intent(in) :: nr

    if ( .not. coefficients%created ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error in delete_coefficients:', &
        ' coefficients has not been created and cannot be deleted ', &
        ' coefficients number in heading = ', nr
      stop
    end if

!   deallocate all allocatables and intialize to defaults
    call deall ( coefficients )

  contains

    subroutine deall ( coefficients )
      type(coefficients_t), intent(out) :: coefficients
    end subroutine deall

  end subroutine delete_single_coefficients


! Delete a structure of type coefficients_t

  subroutine delete_coefficients ( coefficients1, coefficients2, &
    coefficients3, coefficients4, coefficients5 )

    type(coefficients_t), intent(inout) :: coefficients1
    type(coefficients_t), intent(inout), optional :: coefficients2, &
      coefficients3, coefficients4, coefficients5

    call delete_single_coefficients(coefficients1,1)
    if ( present(coefficients2) ) &
      call delete_single_coefficients(coefficients2,2)
    if ( present(coefficients3) ) &
      call delete_single_coefficients(coefficients3,3)
    if ( present(coefficients4) ) &
      call delete_single_coefficients(coefficients4,4)
    if ( present(coefficients5) ) &
      call delete_single_coefficients(coefficients5,5)

  end subroutine delete_coefficients


! delete a (single) array of type coefficients_t

  subroutine delete_single_mcoefficients ( mcoefficients, nr )

    type(coefficients_t), dimension(:), intent(inout) :: mcoefficients
    integer, intent(in) :: nr

    integer :: ic

    if ( .not. all(mcoefficients(:)%created) ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error in delete_mcoefficients:', &
        ' mcoefficients has not been created and cannot be deleted ', &
        ' mcoefficients number in heading = ', nr
      stop
    end if

    do ic = 1, size(mcoefficients)
      call delete_single_coefficients(mcoefficients(ic),nr)
    end do

  end subroutine delete_single_mcoefficients


! Delete a array of type coefficients_t

  subroutine delete_mcoefficients ( mcoefficients1, mcoefficients2, &
    mcoefficients3, mcoefficients4, mcoefficients5 )

    type(coefficients_t), dimension(:), intent(inout) :: mcoefficients1
    type(coefficients_t), dimension(:), intent(inout), optional :: &
      mcoefficients2, mcoefficients3, mcoefficients4, mcoefficients5

    call delete_single_mcoefficients(mcoefficients1,1)
    if ( present(mcoefficients2) ) &
      call delete_single_mcoefficients(mcoefficients2,2)
    if ( present(mcoefficients3) ) &
      call delete_single_mcoefficients(mcoefficients3,3)
    if ( present(mcoefficients4) ) &
      call delete_single_mcoefficients(mcoefficients4,4)
    if ( present(mcoefficients5) ) &
      call delete_single_mcoefficients(mcoefficients5,5)

  end subroutine delete_mcoefficients


! construct a structure of type oldvectors_t

  subroutine create_oldvectors ( oldvectors, nsysvec, nvec, nelvec, nprob, &
    nmesh, nvec1, nsysvec1, nsysvec2, nsysvec3, nelta )

    type(oldvectors_t), intent(inout) :: oldvectors

!   number of sysvectors
    integer, intent(in), optional :: nsysvec

!   number of arrays of sysvectors (one-, two- and three dimensional arrays)
    integer, intent(in), optional :: nsysvec1, nsysvec2, nsysvec3

!   number of vectors
    integer, intent(in), optional :: nvec

!   number of arrays of vectors (one-dimensional arrays)
    integer, intent(in), optional :: nvec1

!   number of elvectors
    integer, intent(in), optional :: nelvec

!   number of problems contained in the structure oldvectors_t
    integer, intent(in), optional :: nprob

!   number of meshes contained in the structure oldvectors_t
    integer, intent(in), optional :: nmesh

!   number of eltree arrays contained in the structure oldvectors_t
    integer, intent(in), optional :: nelta


    integer :: ns, ns1, ns2, ns3, nv, nv1, np, ne, nm, nea


    if ( oldvectors%created ) then
      write(*,'(/a/)') &
        'Error in create_oldvectors: oldvectors has already been created '
      stop
    end if

    ns  = set_optional ( variable=nsysvec, default=0 )
    ns1 = set_optional ( variable=nsysvec1, default=0 )
    ns2 = set_optional ( variable=nsysvec2, default=0 )
    ns3 = set_optional ( variable=nsysvec3, default=0 )
    nv  = set_optional ( variable=nvec, default=0 )
    nv1 = set_optional ( variable=nvec1, default=0 )
    ne  = set_optional ( variable=nelvec, default=0 )
    np  = set_optional ( variable=nprob, default=0 )
    nm  = set_optional ( variable=nmesh, default=0 )
    nea = set_optional ( variable=nelta, default=0 )

    allocate ( oldvectors%s(ns) )
    allocate ( oldvectors%s1(ns1) )
    allocate ( oldvectors%s2(ns2) )
    allocate ( oldvectors%s3(ns3) )
    allocate ( oldvectors%v(nv) )
    allocate ( oldvectors%v1(nv1) )
    allocate ( oldvectors%e(ne) )
    allocate ( oldvectors%p(np) )
    allocate ( oldvectors%m(nm) )
    allocate ( oldvectors%ea(nea) )

    oldvectors%created = .true.

  end subroutine create_oldvectors


! delete a (single) structure of type oldvectors_t

  subroutine delete_single_oldvectors ( oldvectors, nr )

    type(oldvectors_t), intent(inout) :: oldvectors
    integer, intent(in) :: nr

    if ( .not. oldvectors%created ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error in delete_oldvectors:', &
        ' oldvectors has not been created and cannot be deleted ', &
        ' oldvectors number in heading = ', nr
      stop
    end if

!   deallocate all allocatables and intialize to defaults
    call deall ( oldvectors )

  contains

    subroutine deall ( oldvectors )
      type(oldvectors_t), intent(out) :: oldvectors
    end subroutine deall

  end subroutine delete_single_oldvectors


! Delete a structure of type oldvectors_t

  subroutine delete_oldvectors ( oldvectors1, oldvectors2, &
    oldvectors3, oldvectors4, oldvectors5 )

    type(oldvectors_t), intent(inout) :: oldvectors1
    type(oldvectors_t), intent(inout), optional :: oldvectors2, &
      oldvectors3, oldvectors4, oldvectors5

    call delete_single_oldvectors(oldvectors1,1)
    if ( present(oldvectors2) ) call delete_single_oldvectors(oldvectors2,2)
    if ( present(oldvectors3) ) call delete_single_oldvectors(oldvectors3,3)
    if ( present(oldvectors4) ) call delete_single_oldvectors(oldvectors4,4)
    if ( present(oldvectors5) ) call delete_single_oldvectors(oldvectors5,5)

  end subroutine delete_oldvectors


! Now follow all dummy functions used in tfem.


! dummy scalar function

  function dum_func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: dum_func

    write(*,'(/a/)') 'Scalar function (func interface) has not been defined'
    dum_func = 0
    stop

  end function dum_func


! dummy vector function

  function dum_vfunc ( n, nr, x )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: dum_vfunc

    write(*,'(/a/)') 'Vector function (vfunc interface) has not been defined'
    dum_vfunc = 0
    stop

  end function dum_vfunc


! dummy tensor function

  function dum_tfunc ( n, nr, x )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n,n) :: dum_tfunc

    write(*,'(/a/)') 'Tensor function (tfunc interface) has not been defined'
    dum_tfunc = 0
    stop

  end function dum_tfunc


! dummy routine for setting ninti on element level

  subroutine dum_set_ninti_user ( mesh, problem, elgrp, elem, first, &
    last, coefficients, oldvectors )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors

    write(*,'(/a/)') 'Subroutine set_ninti_user has not been defined'
    stop

  end subroutine dum_set_ninti_user


! dummy routine for setting Gauss rule on element level

  subroutine dum_set_Gauss_integration_user ( mesh, problem, elgrp, &
    elem, first, last, coefficients, oldvectors )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors

    write(*,'(/a/)') &
       'Subroutine set_Gauss_integration_user has not been defined'
    stop

  end subroutine dum_set_Gauss_integration_user


! dummy routine for inflow boundary conditions for the c-tensor in DG/BCF

  function dum_cinflow ( n, x, un )

    integer, intent(in) :: n
    real(dp), intent(in), dimension(:) :: x
    real(dp), intent(in) :: un
    real(dp), dimension(n) :: dum_cinflow

    write(*,'(/a/)') 'Function cinflow has not been defined'
    stop
    dum_cinflow = 0
    stop

  end function dum_cinflow


! dummy routine for inflow boundary conditions for the Q-vector in DG/BCF

  function dum_qinflow ( nfield, ncompq, x, un )

    integer, intent(in) :: nfield, ncompq
    real(dp), intent(in), dimension(:) :: x
    real(dp), intent(in) :: un
    real(dp), dimension(nfield,ncompq) :: dum_qinflow

    write(*,'(/a/)') 'Function qinflow has not been defined'
    stop
    dum_qinflow = 0
    stop

  end function dum_qinflow


! dummy routine for the Brownian vector in DG/BCF

  function dum_brownian_vector ( nfield, ncompq )

    integer, intent(in) :: nfield, ncompq
    real(dp), dimension(nfield,ncompq) :: dum_brownian_vector

    write(*,*) 'Function brownian_vector has not been defined'
    dum_brownian_vector = 0
    stop

  end function dum_brownian_vector

end module element_defs_m
