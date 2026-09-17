
! Copyright (C) 2008-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines routines associated to the eltree type with respect to the "volume"

module eltree_volume_m

  use eltree_basic_m
  use meshgen_m

  implicit none

! options for the integration scheme

  type integration_options_t

!   if .true. the weights of elements with eltree%lsign==0 will be divided by
!   two. NOTE: not relevant if there is a submesh at the end of a branch.
    logical :: half = .false.

!   if .true. the elements with eltree%lsign==0 will be filled with a midpoint
!   rule. NOTE: not relevant if there is a submesh at the end of a branch.
    logical :: midp = .false.

!   subelements in smesh (submesh) having a jacobian smaller than
!   epsjac will be omitted.
!   NOTE: only relevant if there is a submesh at the end of a branch.
    real(dp) :: epsjac = 0._dp   ! default=0 (no elements removed).

  end type integration_options_t


! interface for generic volume subroutine

  interface volume
    module procedure volume_eltree
  end interface volume

! interface for generic number_of_subelements

  interface number_of_subelements
    module procedure number_of_subelements_eltree
  end interface number_of_subelements

! interface for generic number_of_subelements (vector)

  interface number_of_subelements_vector
    module procedure number_of_subelements_vector_eltree
  end interface number_of_subelements_vector

! interface for generic number_of_integration_points subroutine

  interface number_of_integration_points
    module procedure number_of_integration_points_eltree
  end interface number_of_integration_points

! interface for generic integration_points subroutine

  interface integration_points
    module procedure integration_points_eltree
  end interface integration_points

contains


! count number of subelements in a tree

  recursive function number_of_subelements_eltree ( eltree, lsign, &
    includeeltree, includesubmesh ) result(numsub)

    type(eltree_t), intent(in) :: eltree

!   if present include elements that have an eltree%lsign or smesh%lsign value
!   given by lsign. Example lsign=[0,1] will only include elements that have
!   an eltree%lsign or smesh%lsign value of 0 or 1.
!   default = all possible values [-1,0,1]
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign

!   include the basic elements in the tree in the counting
!   default = .true.
    logical, intent(in), optional :: includeeltree

!   include the elements in a submesh in the counting
!   default = .true.
    logical, intent(in), optional :: includesubmesh

    integer :: numsub


    integer :: num, inum, ndim
    logical :: lincludeeltree, lincludesubmesh


    numsub = 0

    if ( associated(eltree%subel) ) then

!     continue searching

      ndim = size(eltree%coor,2)
      num = 2**ndim

      do inum = 1, num
        numsub = numsub + number_of_subelements_eltree ( eltree%subel(inum), &
                            lsign, includeeltree, includesubmesh )
      end do

    else

!     end of the branch

      if ( associated(eltree%smesh) ) then

!       further subdivision using submesh

        lincludesubmesh = set_optional ( variable=includesubmesh, &
                                         default=.true. )

        if ( lincludesubmesh ) then
          numsub = number_of_subelements_smesh ( eltree%smesh, lsign )
        else
          numsub = 0
        end if

      else

!       standard case

        lincludeeltree = set_optional ( variable=includeeltree, default=.true. )

        if ( lincludeeltree ) then

          if ( present(lsign) ) then
            if ( any( lsign == eltree%lsign ) ) then
              numsub = 1
            else
              numsub = 0
            end if
          else
            numsub = 1
          end if

        else
          numsub = 0
        end if

      end if

    end if

  end function number_of_subelements_eltree


! count number of subelements in a tree (vector)

  recursive function number_of_subelements_vector_eltree ( eltree, lsign, &
    includeeltree, includesubmesh ) result(numsub)

    type(eltree_t), intent(in) :: eltree

!   count elements that have an eltree%lsign value given by lsign (vector value)
!   Example lsign=[0,1] will count elements that have an eltree%lsign
!   value of 0 or 1. The result is a vector of the same size as lsign.
!   default = all possible values [-1,0,1]
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:) :: lsign

!   include the basic elements in the tree in the counting
!   default = .true.
    logical, intent(in), optional :: includeeltree

!   include the elements in a submesh in the counting
!   default = .true.
    logical, intent(in), optional :: includesubmesh

    integer, dimension(size(lsign)) :: numsub


    integer :: num, inum, ndim
    logical :: lincludeeltree, lincludesubmesh


    numsub = 0

    if ( associated(eltree%subel) ) then

!     continue searching

      ndim = size(eltree%coor,2)
      num = 2**ndim

      do inum = 1, num
        numsub = numsub + number_of_subelements_vector_eltree ( &
                     eltree%subel(inum), lsign, includeeltree, includesubmesh )
      end do

    else

!     end of the branch

      if ( associated(eltree%smesh) ) then

!       further subdivision using submesh

        lincludesubmesh = set_optional ( variable=includesubmesh, &
                                         default=.true. )

        if ( lincludesubmesh ) then
          numsub = number_of_subelements_vector_smesh ( eltree%smesh, lsign )
        else
          numsub = 0
        end if

      else

!       standard case

        lincludeeltree = set_optional ( variable=includeeltree, default=.true. )

        if ( lincludeeltree ) then

          where ( lsign == eltree%lsign )
            numsub = 1
          elsewhere
            numsub = 0
          end where

        else
          numsub = 0
        end if

      end if

    end if

  end function number_of_subelements_vector_eltree


! count number of submesh nodes in a tree

  recursive function number_of_submesh_nodes ( eltree, lsign ) result(numnodes)

    type(eltree_t), intent(inout) :: eltree

!   if present include elements that have an eltree%lsign or smesh%lsign value
!   given by lsign.
!   Example lsign=[0,1] will only include elements that have an eltree%lsign
!   or smesh%lsign value of 0 or 1.
!   default = all possible values [-1,0,1]
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign

    integer :: numnodes

    integer :: num, inum, ndim


    numnodes = 0

    if ( associated(eltree%subel) ) then

!     continue searching

      ndim = size(eltree%coor,2)
      num = 2**ndim

      do inum = 1, num
        numnodes = numnodes + number_of_submesh_nodes ( eltree%subel(inum), &
                                                        lsign )
      end do

    else

!     end of the branch

      if ( associated(eltree%smesh) ) then

!       further subdivision using submesh

        numnodes = number_of_nodes_smesh ( eltree%smesh, lsign )

      else

        numnodes = 0

      end if

    end if

  end function number_of_submesh_nodes


! compute length (1D), area (2D) or volume (3D)
! NOTE: this is the reference volume!

  recursive function volume_eltree ( eltree, lsign ) result(esum)

    type(eltree_t), intent(inout) :: eltree

!   if present include elements that have an eltree%lsign or smesh%lsign value
!   given by lsign.
!   Example lsign=[0,1] will only include elements that have an eltree%lsign
!   or smesh%lsign value of 0 or 1.
!   default = all possible values [-1,0,1]
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign

    integer :: eltype

    real(dp) :: esum


    integer :: num, inum, ndim

    eltype = eltree%eltype


    ndim = size(eltree%coor,2)

    esum = 0

    if ( associated(eltree%subel) ) then

!     continue searching

      num = 2**ndim

      do inum = 1, num
        esum = esum + volume_eltree ( eltree%subel(inum), lsign )
      end do

    else

!     end of the branch

      if ( associated(eltree%smesh) ) then

!       further subdivision using submesh

        esum = volume_smesh ( eltree%smesh, lsign )

      else

!       standard case

        if ( present(lsign) ) then
          if ( .not. any ( lsign == eltree%lsign ) ) return
        end if

!       determine length, area or volume

        select case ( ndim )

          case(1) ! line

            esum = eltree%coor(2,1)-eltree%coor(1,1)

          case(2)

            if ( eltype == 1 ) then ! triangle
              esum = eltree%coor(2,1) - eltree%coor(1,1)
              esum = 0.5_dp * esum * ( eltree%coor(2,2)-eltree%coor(1,2) )
            else if ( eltype == 2 ) then ! quad
              esum = eltree%coor(2,1) - eltree%coor(1,1)
              esum = esum * ( eltree%coor(2,2)-eltree%coor(1,2) )
            end if

          case(3)

            if (eltype == 1 ) then ! tet
              esum = eltree%coor(2,1) - eltree%coor(1,1)
              esum = esum * ( eltree%coor(2,2) - eltree%coor(1,2) ) / 2.0_dp
              esum = esum * ( eltree%coor(2,3) - eltree%coor(1,3) ) / 3.0_dp
            else if ( ndim == 3 .and.  eltype == 2 ) then ! hexa
              esum = eltree%coor(2,1)-eltree%coor(1,1)
              esum = esum * ( eltree%coor(2,2)-eltree%coor(1,2) )
              esum = esum * ( eltree%coor(2,3)-eltree%coor(1,3) )
            end if

          case default

            call errormsg_case_default ( 'volume_eltree', 'ndim', &
              int_value=ndim )

        end select

      end if

    end if

  end function volume_eltree


! count number of integration points defined by a tree

  recursive function number_of_integration_points_eltree ( eltree, lsign, &
    ninti, nintis, integration_options ) result(numint)

    type(eltree_t), intent(inout) :: eltree

!   if present include elements that have an eltree%lsign or smesh%lsign value
!   given by lsign.
!   Example lsign=[0,1] will only include elements that have an eltree%lsign
!   or smesh%lsign value of 0 or 1.
!   default = all possible values [-1,0,1]
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign

!   number of integration points on an eltree subelement (line,quad,hexa).
!   default=1 (mid-point)
    integer, intent(in), optional :: ninti

!   number of integration points on a subelement in smesh (submesh with
!   lines,triangles,tets).
!   NOTE: if subdivide has been called with submesh=.true. nintis must be
!   present.
    integer, intent(in), optional :: nintis

!   The options/parameters for the integration points
!   See the type definition for the defaults.
    type(integration_options_t), intent(in), optional :: integration_options


    logical :: lmidp
    integer :: numint
    type(integration_options_t) :: lintopt ! local options


    integer :: num, inum, ndim


!   initialize local options if in heading
    if ( present(integration_options) ) lintopt = integration_options

    numint = 0

    ndim = size(eltree%coor,2)

    if ( associated(eltree%subel) ) then

!     continue searching

      num = 2**ndim

      do inum = 1, num
        numint = numint + number_of_integration_points_eltree ( &
                                eltree%subel(inum), &
                                lsign, ninti, nintis, integration_options )
      end do

    else

!     end of the branch

      if ( associated(eltree%smesh) ) then

!       further subdivision using submesh

        if ( .not. present(nintis) ) then
          write(*,'(2(/a)/)') &
            'Error in number_of_integration_points_eltree : ', &
            ' nintis not present, whereas smesh has been defined '
          stop
        end if

        numint = number_of_integration_points_smesh ( eltree%smesh, nintis, &
          lsign, lintopt%epsjac )

      else

!       standard case

        lmidp = lintopt%midp

        if ( present(lsign) ) then

          if ( any( lsign == eltree%lsign ) ) then

            if ( .not. present(ninti) .or. &
                 ( lmidp .and. eltree%lsign == 0 ) ) then
              numint = 1 ! mid-point
            else
              numint = ninti
            end if

          else

            numint = 0

          end if

        else

          if ( .not. present(ninti) .or. &
               ( lmidp .and. eltree%lsign == 0 ) ) then
            numint = 1 ! mid-point
          else
            numint = ninti
          end if

        end if

      end if

    end if

  end function number_of_integration_points_eltree


! integration points for a node (recursively)

  recursive subroutine node_integration_points ( eltree, ipnr, x, w, lsign, &
    ninti, xe, we, nintis, xs, ws, integration_options )

    type(eltree_t), intent(inout) :: eltree

!   start pointer for each new subelement
!   should be set to zero on first entry
    integer, intent(inout) :: ipnr

!   the reference coordinates and weights of the integration points
!   in the eltree elements (line,quad,hexa).
!   x(:,ndim), w(:)
    real(dp), intent(inout), dimension(:,:) :: x
    real(dp), intent(inout), dimension(:) :: w

!   if present include elements that have an eltree%lsign or smesh%lsign value
!   given by lsign.
!   Example lsign=[0,1] will only include elements that have an eltree%lsign
!   or smesh%lsign value of 0 or 1.
!   default = all possible values [-1,0,1]
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign

!   number of integration points on a subelement
!   default=1 (mid-point)
    integer, intent(in), optional :: ninti

!   the reference coordinates and weights of the integration points
!   for a single element in the eltree.
!   xe(ninti,ndim), we(ninti)
!   must be present if ninti is present
    real(dp), intent(inout), dimension(:,:), optional :: xe
    real(dp), intent(inout), dimension(:), optional :: we

!   number of integration points on a subelement of smesh
!   NOTE: if subdivide has been called with submesh=.true. nintis must be
!   present.
    integer, intent(in), optional :: nintis

!   the reference coordinates and weights of the integration points
!   in the submesh (smesh)
!   xs(:,ndim), ws(:)
!   NOTE: if subdivide has been called with submesh=.true. xs and ws must be
!   present.
    real(dp), intent(inout), dimension(:,:), optional :: xs
    real(dp), intent(inout), dimension(:), optional :: ws

!   The options/parameters for the integration points
!   See the type definition for the defaults.
    type(integration_options_t), intent(in), optional :: integration_options

    integer :: eltype


    logical :: lhalf, lmidp
    integer :: num, inum, ndim, numint
    real(dp) :: href, h(size(eltree%coor,2))
    type(integration_options_t) :: lintopt ! local options


    eltype = eltree%eltype

!   initialize local options if in heading
    if ( present(integration_options) ) lintopt = integration_options

    ndim = size(eltree%coor,2)

    if ( associated(eltree%subel) ) then

!     continue searching

      num = 2**ndim

      do inum = 1, num
        call node_integration_points ( eltree%subel(inum), ipnr, x, w, lsign, &
                         ninti, xe, we, nintis, xs, ws, integration_options )
      end do

    else

!     end of the branch

      if ( associated(eltree%smesh) ) then

!       further subdivision using submesh

        if ( .not. present(nintis) .or. .not. present(xs) .or. &
             .not. present(ws) ) then
          write(*,'(2(/a)/)') &
            'Error in number_of_integration_points_eltree : ', &
            ' nintis, xs, ws not present, whereas smesh has been defined '
          stop
        end if

!       integration points for smesh

        call integration_points_smesh ( eltree%smesh, x(ipnr+1:,:), &
          w(ipnr+1:), nintis, xs, ws, lsign, lintopt%epsjac, numint )

        ipnr = ipnr + numint

      else

!       standard case

        lmidp = lintopt%midp
        lhalf = lintopt%half

!       size of the reference element

        select case ( ndim )

          case(1) ! line

            href = 2

          case(2)

            if ( eltype == 1  ) then ! triangle
              href = 1
            else if ( eltype == 2 ) then ! quad
              href = 2
            end if

          case(3)

            if ( eltype == 1 ) then ! tetrahedron
              href = 1
            else if ( eltype == 2  ) then ! hexahedron
              href = 2
            end if

          case default

            call errormsg_case_default ( 'node_integration_points', 'ndim', &
              int_value=ndim )

        end select

        h = ( eltree%coor(2,:) - eltree%coor(1,:) ) ! size vector of the element

        if ( present(lsign) ) then

          if ( any( lsign == eltree%lsign ) ) then

            call fill_int

          end if

        else

          call fill_int

        end if

      end if

    end if

  contains

! fill element with integration points and weights

  subroutine fill_int

    integer :: i, j
    real(dp), allocatable :: xnod(:,:)

      if ( .not. present(ninti) .or. ( lmidp .and. eltree%lsign == 0 ) ) then

!       fill with midpoint rule

        ipnr = ipnr + 1

!       find the midpoints

        select case ( ndim )

        case(1) ! line

          x(ipnr,:) = ( eltree%coor(2,:) + eltree%coor(1,:) ) / 2

        case(2)

          if ( eltype == 1  ) then ! triangle

            select case ( eltree%conf )
              case(1)
                x(ipnr,:) = ( eltree%coor(1,:) + &
                              [eltree%coor(2,1), eltree%coor(1,2)] + &
                              [eltree%coor(1,1), eltree%coor(2,2)] ) / 3
              case(2)
                x(ipnr,:) = ( eltree%coor(2,:) + &
                              [eltree%coor(2,1), eltree%coor(1,2)] + &
                              [eltree%coor(1,1), eltree%coor(2,2)] ) / 3
            case default
              call errormsg_case_default ( 'fill_int', 'eltree%conf', &
                int_value=eltree%conf )
            end select

          else if ( eltype == 2 ) then ! quad

            x(ipnr,:) = ( eltree%coor(2,:) + eltree%coor(1,:) ) / 2

          end if

        case(3)

          if ( eltype == 1 ) then ! tetrahedron

            select case ( eltree%conf )
              case(1)
                x(ipnr,:)  = ( mapref ( eltree%coor(1,:) ) + &
                               mapref ( [ eltree%coor(2,1), eltree%coor(1,2), &
                                          eltree%coor(1,3) ] ) + &
                               mapref ( [ eltree%coor(2,1), eltree%coor(2,2), &
                                          eltree%coor(1,3) ] ) + &
                               mapref ( eltree%coor(2,:) ) ) / 4
              case(2)
                x(ipnr,:)  = ( mapref ( eltree%coor(1,:) ) + &
                               mapref ( [ eltree%coor(1,1), eltree%coor(2,2), &
                                          eltree%coor(1,3) ] ) + &
                               mapref ( [ eltree%coor(1,1), eltree%coor(2,2), &
                                          eltree%coor(2,3) ] ) + &
                               mapref ( eltree%coor(2,:) ) ) / 4
              case(3)
                x(ipnr,:)  = ( mapref ( eltree%coor(1,:) ) + &
                               mapref ( [ eltree%coor(1,1), eltree%coor(2,2), &
                                          eltree%coor(1,3) ] ) + &
                               mapref ( [ eltree%coor(2,1), eltree%coor(2,2), &
                                          eltree%coor(1,3) ] ) + &
                               mapref ( eltree%coor(2,:) ) ) / 4
              case(4)
                x(ipnr,:)  = ( mapref ( eltree%coor(1,:) ) + &
                               mapref ( [ eltree%coor(1,1), eltree%coor(1,2), &
                                          eltree%coor(2,3) ] ) + &
                               mapref ( [ eltree%coor(2,1), eltree%coor(1,2), &
                                          eltree%coor(2,3) ] ) + &
                               mapref ( eltree%coor(2,:) ) ) / 4
              case(5)
                x(ipnr,:)  = ( mapref ( eltree%coor(1,:) ) + &
                               mapref ( [ eltree%coor(2,1), eltree%coor(1,2), &
                                          eltree%coor(1,3) ] ) + &
                               mapref ( [ eltree%coor(2,1), eltree%coor(1,2), &
                                          eltree%coor(2,3) ] ) + &
                               mapref ( eltree%coor(2,:) ) ) / 4
              case(6)
                x(ipnr,:)  = ( mapref ( eltree%coor(1,:) ) + &
                               mapref ( [ eltree%coor(1,1), eltree%coor(1,2), &
                                          eltree%coor(2,3) ] ) + &
                               mapref ( [ eltree%coor(1,1), eltree%coor(2,2), &
                                          eltree%coor(2,3) ] ) + &
                               mapref ( eltree%coor(2,:) ) ) / 4
             case default
                call errormsg_case_default ( 'fill_int', 'eltree%conf', &
                  int_value=eltree%conf )
            end select

          else if ( eltype == 2  ) then ! hexahedron

            x(ipnr,:) = ( eltree%coor(2,:) + eltree%coor(1,:) ) / 2

          end if

        case default

          call errormsg_case_default ( 'fill_int', 'ndim', int_value=ndim )

        end select

!       the weights are scaled according to the size of the subelement

        w(ipnr) = h(1)

        do i = 2, ndim
          w(ipnr) = w(ipnr) * h(i)
        end do

!       for triangles, the weight is divided by 2

        if ( eltype == 1 .and. ndim == 2)  w(ipnr) = w(ipnr) / 2

!       for tetrahedrons, the weight is divided by 6

        if ( eltype == 1 .and. ndim == 3)  w(ipnr) = w(ipnr) / 6

        if ( lhalf .and. eltree%lsign == 0 ) w(ipnr) = w(ipnr) / 2

      else

!       fill with supplied Gauss points

        do j = 1, ninti

          select case ( ndim )

            case(1) ! line

              x(ipnr+j,:) = ( eltree%coor(2,:) + eltree%coor(1,:) ) / 2 + &
                              h * xe(j,:) / 2
            case(2)

              if (eltype == 1 ) then ! triangle

!               the location of the vertices

                allocate(xnod(3,2))

                select case ( eltree%conf ) ! find the configuration
                  case(1)
                    xnod(1,:) = eltree%coor(1,:)
                    xnod(2,:) = [ eltree%coor(2,1), eltree%coor(1,2) ]
                    xnod(3,:) = [ eltree%coor(1,1), eltree%coor(2,2) ]
                  case(2)
                    xnod(1,:) = [ eltree%coor(2,1), eltree%coor(1,2) ]
                    xnod(2,:) = eltree%coor(2,:)
                    xnod(3,:) = [ eltree%coor(1,1), eltree%coor(2,2) ]
                  case default
                    call errormsg_case_default ( 'fill_int', 'eltree%conf', &
                      int_value=eltree%conf )
                end select

!               the weights are moved to the correct position inside
!               the subelement using isoparametric elements

                call mapref2 ( xe(j,:), xnod, x(ipnr+j,:) )

                deallocate(xnod)

              else if (eltype == 2 ) then ! quad

                x(ipnr+j,:) = ( eltree%coor(2,:) + eltree%coor(1,:) ) / 2 + &
                                h * xe(j,:) / 2

              end if

          case(3)

            if (eltype == 1 ) then ! tetrahedron

!             the location of the vertices
!             NOTE: mapref is called again to map the reference region for the
!             eltree to the reference region for isoparametric elements.

              allocate(xnod(4,3))

              select case ( eltree%conf ) ! find the configuration
                case(1)
                  xnod(1,:) = mapref ( eltree%coor(1,:) )
                  xnod(2,:) = mapref ( [ eltree%coor(2,1), eltree%coor(1,2), &
                                         eltree%coor(1,3) ] )
                  xnod(3,:) = mapref ( [ eltree%coor(2,1), eltree%coor(2,2), &
                                         eltree%coor(1,3) ] )
                  xnod(4,:) = mapref ( eltree%coor(2,:) )
                case(2)
                  xnod(1,:) = mapref ( eltree%coor(1,:) )
                  xnod(2,:) = mapref ( [ eltree%coor(1,1), eltree%coor(2,2), &
                                         eltree%coor(1,3) ] )
                  xnod(3,:) = mapref ( [ eltree%coor(1,1), eltree%coor(2,2), &
                                         eltree%coor(2,3) ] )
                  xnod(4,:) = mapref ( eltree%coor(2,:) )
                case(3)
                  xnod(1,:) = mapref ( eltree%coor(1,:) )
                  xnod(2,:) = mapref ( [ eltree%coor(1,1), eltree%coor(2,2), &
                                         eltree%coor(1,3) ] )
                  xnod(3,:) = mapref ( [ eltree%coor(2,1), eltree%coor(2,2), &
                                         eltree%coor(1,3) ] )
                  xnod(4,:) = mapref ( eltree%coor(2,:) )
                case(4)
                  xnod(1,:) = mapref ( eltree%coor(1,:) )
                  xnod(2,:) = mapref ( [ eltree%coor(1,1), eltree%coor(1,2), &
                                         eltree%coor(2,3) ] )
                  xnod(3,:) = mapref ( [ eltree%coor(2,1), eltree%coor(1,2), &
                                         eltree%coor(2,3) ] )
                  xnod(4,:) = mapref ( eltree%coor(2,:) )
                case(5)
                  xnod(1,:) = mapref ( eltree%coor(1,:) )
                  xnod(2,:) = mapref ( [ eltree%coor(2,1), eltree%coor(1,2), &
                                         eltree%coor(1,3) ] )
                  xnod(3,:) = mapref ( [ eltree%coor(2,1), eltree%coor(1,2), &
                                         eltree%coor(2,3) ] )
                  xnod(4,:) = mapref ( eltree%coor(2,:) )
                case(6)
                  xnod(1,:) = mapref ( eltree%coor(1,:) )
                  xnod(2,:) = mapref ( [ eltree%coor(1,1), eltree%coor(1,2), &
                                         eltree%coor(2,3) ] )
                  xnod(3,:) = mapref ( [ eltree%coor(1,1), eltree%coor(2,2), &
                                         eltree%coor(2,3) ] )
                  xnod(4,:) = mapref ( eltree%coor(2,:) )
                case default
                  call errormsg_case_default ( 'fill_int', 'eltree%conf', &
                    int_value=eltree%conf )
              end select

!             the weight are moved to the correct position inside the
!             subelement using isoparametric elements

              call mapref2 ( xe(j,:), xnod, x(ipnr+j,:) )

              deallocate(xnod)

            else if (eltype == 2 ) then ! hex

              x(ipnr+j,:) = ( eltree%coor(2,:) + eltree%coor(1,:) ) / 2 + &
                              h * xe(j,:) / 2

            end if

          case default

            call errormsg_case_default ( 'fill_int', 'ndim', int_value=ndim )

          end select

!         the weights are scaled according to the size of the subelement

          w(ipnr+j) = we(j)

          do i = 1, ndim
            w(ipnr+j) = w(ipnr+j) * h(i) / href
          end do

        end do

        if ( lhalf .and. eltree%lsign == 0 ) then
          w(ipnr+1:ipnr+ninti) = w(ipnr+1:ipnr+ninti) / 2
        end if

        ipnr = ipnr + ninti

      end if

    end subroutine fill_int


!   mapref2 is used to find the locations of the integration points (given on
!   the reference element) inside the subeltree. This is using isoparametric
!   elements.
!   NOTE: this is only necessary for triangles and tets because for lines,
!         squares and cubes it is easy to find the locations since the shapes
!         for every subelement and the reference element are equal.

    subroutine mapref2 ( xi, xnod, x )

!     reference coordinates xi(ndim)
      real(dp), dimension(:), intent(in) :: xi

!     nodal coordinates x(nvertices,ndim)
      real(dp), dimension(:,:), intent(in) :: xnod

!     real coordinates x(ndim)
      real(dp), dimension(:), intent(out) :: x

!     shape functions phi(nvertices)
      real(dp), dimension(size(xi,1)+1) :: phi

      integer :: ndim

      ndim = size(xi,1)

      select case(ndim)
        case(2) ! triangle
          phi(1) = 1 - xi(1) - xi(2)
          phi(2) = xi(1)
          phi(3) = xi(2)
        case(3) ! tetrahedron
          phi(1) = 1 - xi(1) - xi(2) - xi(3)
          phi(2) = xi(1)
          phi(3) = xi(2)
          phi(4) = xi(3)
      case default
        call errormsg_case_default ( 'mapref2', 'ndim', int_value=ndim )
      end select

      x = matmul ( phi, xnod )

    end subroutine mapref2

  end subroutine node_integration_points


! integration points for an eltree

  subroutine integration_points_eltree ( eltree, x, w, lsign, ninti, xe, we, &
    nintis, xs, ws, integration_options )

    type(eltree_t), intent(inout) :: eltree

!   the reference coordinates and weights of the integration points
!   x(:,ndim), w(:)
    real(dp), intent(inout), dimension(:,:) :: x
    real(dp), intent(inout), dimension(:) :: w

!   if present include elements that have an eltree%lsign or smesh%lsign value
!   given by lsign.
!   Example lsign=[0,1] will only include elements that have an eltree%lsign
!   or smesh%lsign value of 0 or 1.
!   default = all possible values [-1,0,1]
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign

!   number of integration points on a subelement
!   default=1 (mid-point)
    integer, intent(in), optional :: ninti

!   the reference coordinates and weights of the integration points
!   for a single element
!   xe(ninti,ndim), we(ninti)
!   must be present if ninti is present
    real(dp), intent(inout), dimension(:,:), optional :: xe
    real(dp), intent(inout), dimension(:), optional :: we

!   number of integration points on a subelement of smesh
!   NOTE: if subdivide has been called with submesh=.true. nintis must be
!   present.
    integer, intent(in), optional :: nintis

!   the reference coordinates and weights of the integration points
!   in the submesh (smesh)
!   xs(:,ndim), ws(:)
!   NOTE: if subdivide has been called with submesh=.true. xs and ws must be
!   present.
    real(dp), intent(inout), dimension(:,:), optional :: xs
    real(dp), intent(inout), dimension(:), optional :: ws

!   The options/parameters for the integration points
!   See the type definition for the defaults.
    type(integration_options_t), intent(in), optional :: integration_options

    integer :: ipnr

    ipnr = 0

    call node_integration_points ( eltree, ipnr, x, w, lsign, ninti, xe, we, &
      nintis, xs, ws, integration_options )

  end subroutine integration_points_eltree


! create a mesh from an eltree

  subroutine eltree_to_mesh ( eltree, mesh, lsign, mapcoor )

    type(eltree_t), intent(inout) :: eltree

!   the basic mesh created
!   NOTE: if mapcoor is not present the mesh is defined on the reference element
    type(mesh_t), intent(inout) :: mesh

!   if present include elements that have an eltree%lsign or smesh%lsign value
!   given by lsign.
!   Example lsign=[0,1] will only include elements that have an eltree%lsign
!   or smesh%lsign value of 0 or 1.
!   default = all possible values [-1,0,1]
!   NOTE: for submesh elements, the value of 0 might mean the levelset
!   function is not well defined.
    integer, intent(in), dimension(:), optional :: lsign

    integer :: eltype

!   if present map reference coordinates to real coordinates.
!   NOTE: multiple points need to be handled at the same time,
!   with x(nnodes,ndim)
    optional :: mapcoor
    interface
      function mapcoor ( x ) result(m)
        use kind_defs_m
        implicit none
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1),size(x,2)) :: m
      end function mapcoor
    end interface

    integer :: nelem1, nelem2, nnodes, elshape_mesh, node, elem, ndim
    type(mesh_t) :: mesh1, mesh2

    eltype = eltree%eltype

    ndim = size(eltree%coor,2)

!   Basic eltree mesh

    nelem1 = number_of_subelements ( eltree, lsign, includesubmesh=.false. )

    select case ( ndim )

      case(1)

        nnodes = 2 * nelem1
        elshape_mesh = 1  ! 2-node line element

      case(2)

        if ( eltype == 1) then ! triangle
          nnodes = 3 * nelem1
          elshape_mesh = 3   ! 3-node triangle element
        else if ( eltype == 2 ) then ! quad
          nnodes = 4 * nelem1
          elshape_mesh = 5  ! 4-node quad element
        end if

      case(3)

        if ( eltype == 1 ) then ! tet
          nnodes = 4 * nelem1
          elshape_mesh = 11 ! 4-node tetrahedron element
        else if ( eltype == 2 ) then ! hexa
          nnodes = 8 * nelem1
          elshape_mesh = 13  ! 8-node hexahedron element
        end if

      case default

        call errormsg_case_default ( 'eltree_to_mesh', 'ndim', int_value=ndim )

    end select

    call mesh_skeleton ( mesh1, nnodes, nelem1, elshape_mesh, ndim )

    node = 0
    elem = 0

    call fill_coor_topology_mesh1 ( eltree, mesh1, node, elem, lsign, mapcoor )

!   submesh mesh

    nelem2 = number_of_subelements ( eltree, lsign, includeeltree=.false. )

    if ( nelem2 > 0 ) then

      nnodes = number_of_submesh_nodes ( eltree, lsign )

      select case ( ndim )
        case(1)
          elshape_mesh = 1  ! 2-node line element
        case(2)
          elshape_mesh = 3  ! 3-node triangular element
        case(3)
          elshape_mesh = 11  ! 4-node tetrahedral element
        case default
          call errormsg_case_default ( 'eltree_to_mesh', 'ndim', &
            int_value=ndim )
      end select

      call mesh_skeleton ( mesh2, nnodes, nelem2, elshape_mesh, ndim )

      node = 0
      elem = 0

      call fill_coor_topology_mesh2 ( eltree, mesh2, node, elem, lsign, &
        mapcoor )

    end if

!   copy or merge

    if ( nelem2 > 0 ) then
      call mesh_merge ( mesh1, mesh2, mesh, nogroupmerge=.true. )
      call delete ( mesh1, mesh2 )
    else
      call mesh_convert_copy ( mesh1, mesh )
      call delete ( mesh1 )
    end if

  end subroutine eltree_to_mesh


! fill mesh with coordinates and topology to basic eltree

  recursive subroutine fill_coor_topology_mesh1 ( eltree, mesh, node, elem, &
    lsign, mapcoor )

    type(eltree_t), intent(inout) :: eltree
    type(mesh_t), intent(inout) :: mesh
    integer, intent(inout) :: node, elem
    integer, intent(in), dimension(:), optional :: lsign
    integer :: eltype

!   if present map reference coordinates to real coordinates
!   NOTE: multiple points need to be handled at the same time,
!   with x(nnodes,ndim)
    optional :: mapcoor
    interface
      function mapcoor ( x ) result(m)
        use kind_defs_m
        implicit none
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1),size(x,2)) :: m
      end function mapcoor
    end interface


    integer :: num, inum, ndim, i
    logical :: add

    eltype = eltree%eltype

    ndim = size(eltree%coor,2)

    if ( associated(eltree%subel) ) then

!     continue searching

      num = 2**ndim

      do inum = 1, num
        call fill_coor_topology_mesh1 ( eltree%subel(inum), mesh, node, elem, &
         lsign, mapcoor )
      end do

    else

!     end of the branch

      if ( associated(eltree%smesh) ) return ! submesh elements not included

      if ( present(lsign) ) then
        add = any ( lsign == eltree%lsign )
      else
        add = .true.
      end if

      if ( add ) then

        elem = elem + 1

!       number of vertices

        select case ( ndim )

          case(1)

            num = 2

          case(2)

            if ( eltype == 1 ) then ! triangle
              num = 3
            else if ( eltype == 2 ) then ! quad
              num = 4
            end if

          case(3)

            if ( eltype == 1 ) then ! tet
              num = 4
            else if ( eltype == 2 ) then ! hexa
              num = 8
            end if

          case default

            call errormsg_case_default ( 'fill_coor_topology_mesh1', 'ndim', &
              int_value=ndim )

        end select

        mesh%topology(1)%a(:,elem) = [(node+i,i=1,num)]

!       find the location of the vertices

        select case ( ndim )

          case(1) ! line

            mesh%coor(node+1,:) = eltree%coor(1,:)
            mesh%coor(node+2,:) = eltree%coor(2,:)

          case(2)

            if ( eltype == 1 ) then ! triangle

              select case ( eltree%conf )
                case(1)
                  mesh%coor(node+1,:) = eltree%coor(1,:)
                  mesh%coor(node+2,:) = [ eltree%coor(2,1), eltree%coor(1,2) ]
                  mesh%coor(node+3,:) = [ eltree%coor(1,1), eltree%coor(2,2) ]
                case(2)
                  mesh%coor(node+1,:) = [ eltree%coor(2,1), eltree%coor(1,2) ]
                  mesh%coor(node+2,:) = eltree%coor(2,:)
                  mesh%coor(node+3,:) = [ eltree%coor(1,1), eltree%coor(2,2) ]
                case default
                  call errormsg_case_default ( 'fill_coor_topology_mesh1', &
                    'eltree%conf', int_value=eltree%conf )
              end select

            else if ( eltype == 2 ) then ! quad

              mesh%coor(node+1,:) = eltree%coor(1,:)
              mesh%coor(node+2,:) = [ eltree%coor(2,1), eltree%coor(1,2) ]
              mesh%coor(node+3,:) = eltree%coor(2,:)
              mesh%coor(node+4,:) = [ eltree%coor(1,1), eltree%coor(2,2) ]

            end if

          case(3)

            if ( eltype == 1 ) then ! tet

              select case ( eltree%conf )
                case(1)
                   mesh%coor(node+1,:) = mapref ( eltree%coor(1,:) )
                   mesh%coor(node+2,:) = mapref ( [ eltree%coor(2,1), &
                                  eltree%coor(1,2), eltree%coor(1,3) ] )
                   mesh%coor(node+3,:) = mapref ( [ eltree%coor(2,1), &
                                  eltree%coor(2,2), eltree%coor(1,3) ] )
                   mesh%coor(node+4,:) = mapref ( eltree%coor(2,:) )
                case(2)
                   mesh%coor(node+1,:) = mapref ( eltree%coor(1,:) )
                   mesh%coor(node+2,:) = mapref ( [ eltree%coor(1,1), &
                                  eltree%coor(2,2), eltree%coor(1,3) ] )
                   mesh%coor(node+3,:) = mapref ( [ eltree%coor(1,1), &
                                  eltree%coor(2,2), eltree%coor(2,3) ] )
                   mesh%coor(node+4,:) = mapref ( eltree%coor(2,:) )
                case(3)
                   mesh%coor(node+1,:) = mapref ( eltree%coor(1,:) )
                   mesh%coor(node+2,:) = mapref ( [ eltree%coor(1,1), &
                                  eltree%coor(2,2), eltree%coor(1,3) ] )
                   mesh%coor(node+3,:) = mapref ( [ eltree%coor(2,1), &
                                  eltree%coor(2,2), eltree%coor(1,3) ] )
                   mesh%coor(node+4,:) = mapref ( eltree%coor(2,:) )
                case(4)
                   mesh%coor(node+1,:) = mapref ( eltree%coor(1,:) )
                   mesh%coor(node+2,:) = mapref ( [ eltree%coor(1,1), &
                                  eltree%coor(1,2), eltree%coor(2,3) ] )
                   mesh%coor(node+3,:) = mapref ( [ eltree%coor(2,1), &
                                  eltree%coor(1,2), eltree%coor(2,3) ] )
                   mesh%coor(node+4,:) = mapref ( eltree%coor(2,:) )
                case(5)
                   mesh%coor(node+1,:) = mapref ( eltree%coor(1,:) )
                   mesh%coor(node+2,:) = mapref ( [ eltree%coor(2,1), &
                                  eltree%coor(1,2), eltree%coor(1,3) ] )
                   mesh%coor(node+3,:) = mapref ( [ eltree%coor(2,1), &
                                  eltree%coor(1,2), eltree%coor(2,3) ] )
                   mesh%coor(node+4,:) = mapref ( eltree%coor(2,:) )
                case(6)
                   mesh%coor(node+1,:) = mapref ( eltree%coor(1,:) )
                   mesh%coor(node+2,:) = mapref ( [ eltree%coor(1,1), &
                                  eltree%coor(1,2), eltree%coor(2,3) ] )
                   mesh%coor(node+3,:) = mapref ( [ eltree%coor(1,1), &
                                  eltree%coor(2,2), eltree%coor(2,3) ] )
                   mesh%coor(node+4,:) = mapref ( eltree%coor(2,:) )
                case default
                  call errormsg_case_default ( 'fill_coor_topology_mesh1', &
                    'eltree%conf', int_value=eltree%conf )
              end select

            else if ( eltype == 2 ) then ! hexa

              mesh%coor(node+1,:) = eltree%coor(1,:)
              mesh%coor(node+2,:) = [ eltree%coor(2,1), eltree%coor(1,2), &
                                      eltree%coor(1,3) ]
              mesh%coor(node+3,:) = [ eltree%coor(2,1), eltree%coor(2,2), &
                                      eltree%coor(1,3) ]
              mesh%coor(node+4,:) = [ eltree%coor(1,1), eltree%coor(2,2), &
                                      eltree%coor(1,3) ]
              mesh%coor(node+5,:) = [ eltree%coor(1,1), eltree%coor(1,2), &
                                      eltree%coor(2,3) ]
              mesh%coor(node+6,:) = [ eltree%coor(2,1), eltree%coor(1,2), &
                                      eltree%coor(2,3) ]
              mesh%coor(node+7,:) = eltree%coor(2,:)
              mesh%coor(node+8,:) = [ eltree%coor(1,1), eltree%coor(2,2), &
                                      eltree%coor(2,3) ]

            end if

          case default

            call errormsg_case_default ( 'fill_coor_topology_mesh1', 'ndim', &
              int_value=ndim )

        end select

        if ( present(mapcoor) ) then
          mesh%coor(node+1:node+num,:) = mapcoor( mesh%coor(node+1:node+num,:) )
        end if

        node = node + num

      end if

    end if

  end subroutine fill_coor_topology_mesh1


! fill mesh with coordinates and topology according to submesh

  recursive subroutine fill_coor_topology_mesh2 ( eltree, mesh, node, elem, &
    lsign, mapcoor )

    type(eltree_t), intent(inout) :: eltree
    type(mesh_t), intent(inout) :: mesh
    integer, intent(inout) :: node, elem
    integer, intent(in), dimension(:), optional :: lsign

!   if present map reference coordinates to real coordinates
!   NOTE: multiple points need to be handled at the same time, &
!   with x(nnodes,ndim)
    optional :: mapcoor
    interface
      function mapcoor ( x )
        use kind_defs_m
        implicit none
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(size(x,1),size(x,2)) :: mapcoor
      end function mapcoor
    end interface

    integer :: num, inum, ndim
    type(smesh_t) :: smesh1

    ndim = size(eltree%coor,2)

    if ( associated(eltree%subel) ) then

!     continue searching

      num = 2**ndim

      do inum = 1, num
        call fill_coor_topology_mesh2 ( eltree%subel(inum), mesh, node, elem, &
          lsign, mapcoor )
      end do

    else

!     end of the branch

      if ( .not. associated(eltree%smesh) ) return ! only submesh elements
                                                   ! are included

!     limit elements to lsign

      call smesh_to_smesh ( eltree%smesh, smesh1, lsign )

      if ( smesh1%nelem > 0 ) then

!       add elements

        mesh%topology(1)%a(:,elem+1:elem+smesh1%nelem) = smesh1%topology + node

        num = smesh1%nnodes

        if ( present(mapcoor) ) then
          mesh%coor(node+1:node+num,:) = mapcoor ( smesh1%coor )
        else
          mesh%coor(node+1:node+num,:) = smesh1%coor
        end if

        node = node + smesh1%nnodes
        elem = elem + smesh1%nelem

      end if

      call delete ( smesh1 )

    end if

  end subroutine fill_coor_topology_mesh2

end module eltree_volume_m
