
! Copyright (C) 2008-2016 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines type for the particles to be tracked and associated routines

module particle_tracking_m

  use kind_defs_m

  implicit none


! type for parameters/options for the tracking

  type tracking_options_t

!   accuracy of odeint

    real(dp) :: eps = 1e-4_dp   ! accuracy parameter eps
    real(dp) :: dtmin = 0.0_dp  ! minimum time step

!   splitting of elements

    real(dp) :: hmax = 0.01_dp   ! maximum element length
    real(dp) :: hcmax = 0.05_dp  ! maximum element length when angles too large
    real(dp) :: alphac = 0.5_dp  ! critical angle between elements

!   choice to re-order elements after splitting
    logical :: reorder = .false.

  end type tracking_options_t


! type definition of a curve defined by a number of points and linear elements
! in between the points

  type pcurve_t

    logical :: closed = .false.  ! curve is closed?
    integer :: nnodes = 0        ! number of nodal points
    integer :: nelem  = 0        ! number of elements

    real(dp) :: time = 0._dp     ! current time
    real(dp) :: prevtime = 0._dp ! previous time

    integer :: startp = 0        ! start point (for open curve)
    integer :: endp = 0          ! end point (for open curve)

!   connections of elements to nodal points
!   topology(node,elem) gives the global node number in element elem for local
!   node (node=1 or 2).
    integer, allocatable, dimension(:,:) :: topology

!   elements connected to nodes
!   nodelem(elem,node), elem=1,2 gives the element numbers connected to
!   the global node.
    integer, allocatable, dimension(:,:) :: nodelem

!   coordinates coor(nnodes,2) corresponding to time
    real(dp), allocatable, dimension(:,:) :: coor

!   coordinates coor(nnodes,2) corresponding to prevtime
    real(dp), allocatable, dimension(:,:) :: prevcoor

  end type pcurve_t

! interface for generic length subroutine

  interface length
    module procedure length_pcurve
  end interface length

! interface for generic area subroutine

  interface area
    module procedure area_pcurve
  end interface area

! interface for generic angle subroutine

  interface angle
    module procedure angle_threep
  end interface angle

! interface for generic create subroutine

  interface create
    module procedure create_pcurve
  end interface create


! interface for generic delete subroutine

  interface delete
    module procedure delete_pcurve
  end interface delete

contains


! create initial pcurve

  subroutine create_pcurve ( pcurve, max_nnodes, closed, coor )

    type(pcurve_t), intent(inout) :: pcurve

!   maximum number of nodes
    integer, intent(in) :: max_nnodes

!   is the curve closed or not?
    logical, intent(in) :: closed

!   coordinates of initial points coor(nnodes,2)
    real(dp), dimension(:,:), intent(in) :: coor


    integer :: nnodes, nelem, max_nelem, elem, node


    if ( max_nnodes < size(coor,1) ) then
      write(*,'(/a/)') 'Error create_pcurve: max_nnodes < size(coor,1)'
      stop
    end if

!   some parameters

    nnodes = size(coor,1)
    if ( closed ) then
      max_nelem = max_nnodes
      nelem = nnodes
    else
      max_nelem = max_nnodes - 1
      nelem = nnodes - 1
    end if

!   fill structure

    pcurve%nnodes = nnodes
    pcurve%nelem = nelem

    allocate ( pcurve%topology(2,max_nelem), pcurve%nodelem(2,max_nnodes) )
    allocate ( pcurve%coor(max_nnodes,2), pcurve%prevcoor(max_nnodes,2) )

    pcurve%coor(1:nnodes,:) = coor
    pcurve%prevcoor(1:nnodes,:) = coor

    do elem = 1, nelem
      pcurve%topology(:,elem)= [ elem, elem+1 ]
    end do

    do node = 1, nnodes
      pcurve%nodelem(:,node)= [ node-1, node ]
    end do

    if ( closed ) then
      pcurve%closed = .true.
      pcurve%topology(2,nelem) = 1
      pcurve%nodelem(1,1) = nelem
    else
      pcurve%closed = .false.
      pcurve%nodelem(1,1) = 0
      pcurve%nodelem(2,nnodes) = 0
    end if

!   add start and end point
    pcurve%startp = 1
    pcurve%endp   = pcurve%nnodes

  end subroutine create_pcurve


! delete pcurve

  subroutine delete_pcurve ( pcurve )

    type(pcurve_t), intent(inout) :: pcurve

    deallocate ( pcurve%topology, pcurve%nodelem )
    deallocate ( pcurve%coor, pcurve%prevcoor )

    pcurve%closed = .false.
    pcurve%nnodes = 0
    pcurve%nelem = 0
    pcurve%time = 0
    pcurve%prevtime = 0
    pcurve%startp = 0
    pcurve%endp = 0

  end subroutine delete_pcurve

! step points

  subroutine step_points ( tracking_options, pcurve, time1, time2, deltat1, &
    derivs, point1, point2 )

    use nr

    type(tracking_options_t), intent(in) :: tracking_options
    type(pcurve_t), intent(inout) :: pcurve

!   starting time end time and guessed first time step (h1 in odeint)
    real(dp), intent(in) :: time1, time2, deltat1

    interface
      subroutine derivs(t,y,dydt)
        use kind_defs_m
        implicit none
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(:), intent(out) :: dydt
      end subroutine derivs
    end interface

    integer, intent(in), optional :: point1, point2


    integer :: ip, p1, p2
    real(dp) :: ystart(2)


!   which points?

    if ( present(point1) .and. present(point2) ) then
!     specified range of objects only
      p1 = point1
      p2 = point2
    else if ( present(point1) ) then
!     one point only
      p1 = point1
      p2 = point1
    else
!     all points
      p1 = 1
      p2 = pcurve%nnodes
    end if

!   do one step with odeint

    do ip = p1, p2
      ystart = pcurve%coor(ip,:)
      call odeint ( ystart, time1, time2, tracking_options%eps, deltat1, &
        tracking_options%dtmin, derivs, rkqs )
      pcurve%coor(ip,:) = ystart
    end do

    pcurve%time = time2

  end subroutine step_points


! split elements

  subroutine split_elements ( tracking_options, pcurve, deltat1, derivs )

    type(tracking_options_t), intent(in) :: tracking_options
    type(pcurve_t), intent(inout) :: pcurve

!   guessed first time step (h1 in odeint)
    real(dp), intent(in) :: deltat1

    interface
      subroutine derivs(t,y,dydt)
        use kind_defs_m
        implicit none
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(in) :: y
        real(dp), dimension(:), intent(out) :: dydt
      end subroutine derivs
    end interface


    logical :: markelements(pcurve%nelem), keepl
    integer :: ip, elem1, elem2, p(2), elem
    integer :: nnodesnew, nnodesold, nelemnew, nelemold
    real(dp) :: angles(pcurve%nnodes)
    real(dp) :: h, x1(2), x2(2), x3(2), elmangles(2), v(2)


!   compute angles in nodal points

    if ( pcurve%closed ) then

!     closed curve

      do ip = 1, pcurve%nnodes

        elem1 = pcurve%nodelem(1,ip)
        elem2 = pcurve%nodelem(2,ip)
        x1 = pcurve%coor(pcurve%topology(1,elem1),:)
        x2 = pcurve%coor(pcurve%topology(2,elem1),:)
        x3 = pcurve%coor(pcurve%topology(2,elem2),:)
        angles(ip) = abs ( angle ( x1, x2, x3 ) )

      end do

    else

!     open curve

      do ip = 1, pcurve%nnodes

!       skip start and end point
        if ( ip == pcurve%startp .or. ip == pcurve%endp ) cycle

        elem1 = pcurve%nodelem(1,ip)
        elem2 = pcurve%nodelem(2,ip)
        x1 = pcurve%coor(pcurve%topology(1,elem1),:)
        x2 = pcurve%coor(pcurve%topology(2,elem1),:)
        x3 = pcurve%coor(pcurve%topology(2,elem2),:)
        angles(ip) = abs ( angle ( x1, x2, x3 ) )

      end do

!     zero angle at the start and end points
      angles(pcurve%startp) = 0
      angles(pcurve%endp)   = 0

    end if

!   mark elements for splitting

    markelements = .false.

    do elem = 1, pcurve%nelem

      p = pcurve%topology(:,elem)

      v = pcurve%coor(p(2),:) - pcurve%coor(p(1),:)
      h = sqrt( v(1)**2 + v(2)**2 )

      elmangles = abs(angles(p))

      markelements(elem) = h > tracking_options%hmax .or. &
        ( h > tracking_options%hcmax .and. &
             any ( elmangles > tracking_options%alphac ) )

    end do

!   split elements

    if ( pcurve%nnodes + count(markelements) > size(pcurve%coor,1) ) then
      write(*,'(/a/)') 'Error split_element: not enough space in pcurve'
      stop
    end if

    nnodesnew = pcurve%nnodes
    nelemnew = pcurve%nelem

    do elem = 1, pcurve%nelem

      if ( .not. markelements(elem) ) cycle

      nnodesnew = nnodesnew + 1
      nelemnew  = nelemnew + 1

      p = pcurve%topology(:,elem)

      pcurve%coor(nnodesnew,:) = &
                 ( pcurve%prevcoor(p(2),:) + pcurve%prevcoor(p(1),:) ) / 2

      pcurve%topology(2,elem) = nnodesnew
      pcurve%topology(1,nelemnew) = nnodesnew
      pcurve%topology(2,nelemnew) = p(2)

      pcurve%nodelem(1,p(2)) = nelemnew
      pcurve%nodelem(1,nnodesnew) = elem
      pcurve%nodelem(2,nnodesnew) = nelemnew

    end do

!   update rest of pcurve

    nnodesold = pcurve%nnodes
    nelemold = pcurve%nelem
    pcurve%nnodes = nnodesnew
    pcurve%nelem = nelemnew

!   update coordinates of new points

    call step_points ( tracking_options, pcurve, pcurve%prevtime, pcurve%time, &
      deltat1, derivs, point1=nnodesold+1, point2=nnodesnew )

    if ( tracking_options%reorder ) then
      call reorder_pcurve ( pcurve )
    end if

  end subroutine split_elements


!   Utility function to reorder a pcurve such that all nodes are order from low
!   to high. By default new nodes inserted into the pcurve after calling
!   split_elements are not ordered, which may be impractical for postprocessing
!   purposes.

  subroutine reorder_pcurve ( pcurve1, pcurve2 )

    type(pcurve_t), intent(inout) :: pcurve1
    type(pcurve_t), intent(inout), optional :: pcurve2

    real(dp), dimension(:,:), allocatable :: pcoor, prevpcoor
    real(dp) :: ptime, prevptime
    integer :: npointsmax
    integer :: ielem, inode, ipoint

!   maximum number of points from pcurve coor size
    npointsmax = size(pcurve1%coor,1)

    ielem = 1
    ipoint = 1

    allocate ( pcoor(pcurve1%nnodes,2), prevpcoor(pcurve1%nnodes,2))
    ptime = pcurve1%time
    prevptime = pcurve1%prevtime

    pcoor(1,1) = pcurve1%coor(1,1)
    pcoor(1,2) = pcurve1%coor(1,2)
    prevpcoor(1,1) = pcurve1%prevcoor(1,1)
    prevpcoor(1,2) = pcurve1%prevcoor(1,2)
    inode = pcurve1%topology(2,ielem)

!   Loop over the number of nodes, jumping from one element to the next using
!   the topology and nodelem.

    do while ( inode /= 1 .and. ielem /= 0 )
      ipoint = ipoint + 1
      pcoor(ipoint,1) = pcurve1%coor(inode,1)
      pcoor(ipoint,2) = pcurve1%coor(inode,2)
      prevpcoor(ipoint,1) = pcurve1%prevcoor(inode,1)
      prevpcoor(ipoint,2) = pcurve1%prevcoor(inode,2)
      if ( pcurve1%nodelem(2,inode) == ielem ) then
        ielem = pcurve1%nodelem(1,inode)
      else
        ielem = pcurve1%nodelem(2,inode)
      end if
      inode = pcurve1%topology(2,ielem)
    end do

!   Build the new pcurve from the ordered coordinates, or overwrite the old one

    if ( .not. present(pcurve2) ) then
      call delete ( pcurve1 )
      call create_pcurve ( pcurve1, max_nnodes=npointsmax, closed=.true., &
        coor=pcoor )
      pcurve1%time = ptime
      pcurve1%prevtime = prevptime
      pcurve1%prevcoor = prevpcoor
    else
      call create_pcurve ( pcurve2, max_nnodes=npointsmax, closed=.true., &
        coor=pcoor )
      pcurve2%time = ptime
      pcurve2%prevtime = prevptime
      pcurve2%prevcoor = prevpcoor
    end if

  end subroutine reorder_pcurve


! compute length of curve

  function length_pcurve ( pcurve )

    type(pcurve_t), intent(in) :: pcurve

    real(dp) :: length_pcurve


    integer :: p(2), elem
    real(dp) :: v(2)

    length_pcurve = 0

    do elem = 1, pcurve%nelem
      p = pcurve%topology(:,elem)
      v = pcurve%coor(p(2),:) - pcurve%coor(p(1),:)
      length_pcurve = length_pcurve + sqrt( v(1)**2 + v(2)**2 )
    end do

  end function length_pcurve


! compute area of closed curve

  function area_pcurve ( pcurve )

    type(pcurve_t), intent(in) :: pcurve

    real(dp) :: area_pcurve


    integer :: p(2), ielem
    real(dp) :: v(2), s(2)


    if ( .not. pcurve%closed ) then
      write(*,'(/a/)') 'Error: area computation not valid for open curve'
      stop
    end if

    area_pcurve = 0

    do ielem = 1, pcurve%nelem
      p = pcurve%topology(:,ielem)
      v = pcurve%coor(p(2),:) - pcurve%coor(p(1),:)
      s = pcurve%coor(p(2),:) + pcurve%coor(p(1),:)
      area_pcurve = area_pcurve + (v(1)*s(2) - v(2)*s(1))/4
    end do

  end function area_pcurve


! compute angle between vectors x(2)-x(1) and x(3)-x(2)

  function angle_threep ( x1, x2, x3 )

    real(dp), dimension(2), intent(in) :: x1, x2, x3
    real(dp) :: angle_threep

    real(dp) :: v1(2), v2(2), lv1, lv2, valcos

    v1 = x2 - x1
    v2 = x3 - x2
    lv1 = sqrt(sum(v1**2))
    lv2 = sqrt(sum(v2**2))

    valcos = dot_product(v1,v2)/lv1/lv2

    valcos = min(1.0_dp,valcos)
    valcos = max(-1.0_dp,valcos)

    angle_threep = acos ( valcos )

  end function angle_threep

end module particle_tracking_m
