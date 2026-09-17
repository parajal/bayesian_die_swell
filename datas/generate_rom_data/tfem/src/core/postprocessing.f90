
! Copyright (C) 2005-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Routines, for post-processing

module postprocessing_m

  use kind_defs_m
  use mesh_m
  use problem_defs_m
  use pos_array_m, only: pos_array_node, pos_array_vec_node
  use system_defs_m, only: sysvector_t
  use system_matrix_m, only: find_elements_integrationpoints, &
                             reset_elements_integrationpoints, &
                             int_array_1d_p
  use vector_defs_m, only: vector_t
  use element_defs_m
  use set_optional_m

  implicit none


! type definition of sample

  type sample_t

    logical :: created = .false. ! has the sample been created?
    integer :: nnodes = 0  ! number of actual sampling points
    integer :: ndim = 0    ! dimension of space of sample
    integer :: ndegfd = 0  ! number of degrees of freedom in each sample point

!   the data in the format u(1:nnodes,1:ndegfd)
!   the actual value of size(u,1) can be larger than nnodes
!   the actual value of size(u,2) can be larger than ndegfd
    real(dp), allocatable, dimension(:,:) :: u

!   the coordinates of the data sample in the format coor(1:nnodes,1:ndim)
!   the actual value of size(coor,1) can be larger than nnodes
!   the actual value of size(coor,2) can be larger than ndim
    real(dp), allocatable, dimension(:,:) :: coor

!   valid(1:nnodes) indicates if the corresponding sample value in u is valid.
!   The size of this array is only defined to be larger than zero if
!   create_sample has been called with allnodes=.true.
!   the actual value of size(valid) can be larger than nnodes
    logical, allocatable, dimension(:) :: valid

  end type sample_t


! interface for generic create subroutine

  interface create
    module procedure create_sample
  end interface create


! interface for generic delete subroutine

  interface delete
    module procedure delete_sample
  end interface delete


contains


! create sample (reserve memory only). Delete old sample if necessary.

  subroutine create_sample ( mesh, sample, ndegfd, object, nodes, points, &
    curve, surface, nnodes, allnodes )

    type(mesh_t), intent(in) :: mesh
    type(sample_t), intent(inout) :: sample

!   maximum number of degrees of freedom in each sample point. Note, that
!   ndegfd must be at least equal to the number of degrees of freedom that
!   will actually be filled in the routine fill_sample.
    integer, intent(in) :: ndegfd

!   if present then sample is on objectnr=object
    integer, intent(in), optional :: object

!   if present, the sample is filled from the value in the given nodal points
    integer, dimension(:), intent(in), optional :: nodes

!   if present then sample is in points given by pointnr=points(:)
    integer, dimension(:), intent(in), optional :: points

!   if present then sample is on curvenr=curve
    integer, intent(in), optional :: curve

!   if present then sample is on surfacenr=surface
    integer, intent(in), optional :: surface

!   if present then it gives the maximum number of sampling points
    integer, intent(in), optional :: nnodes

!   if present and .true. samples of "all nodes" are defined (no nodes are
!   skipped). An additional logical array valid in the sample structure is
!   defined that indicates whether the corresponding sample value is valid
!   or not.
!   default=.false.
    logical, intent(in), optional :: allnodes


!   Create the structure sample with enough memory to store all the sample data.
!   Basically the two parameters ndegfd and nnodes (together with mesh%ndim)
!   are sufficient.
!   The parameters object, curve and surface are there to help to
!   define the (maximum) number of sampling points in a more easy way.
!   The value of the object, points, curve or surface is not stored.


    logical :: lallnodes
    integer :: lnnodes


    call check ( mesh, 'create_sample' )
    lallnodes = set_optional ( variable=allnodes, default=.false. )

!   choose keyword

    if ( present(object) ) then

      if ( object < 1 .or. object > mesh%nobjects ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error: object in the heading of create_sample is ', &
          'out of range. object is ', object, &
          'whereas the number of objects is ', mesh%nobjects
        stop
      end if

      lnnodes = mesh%objects(object)%nnodes

    else if ( present(curve) ) then

      if ( curve < 1 .or. curve > mesh%ncurves ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error: curve in the heading of create_sample is ', &
          'out of range. curve is ', curve, &
          'whereas the number of curves is ', mesh%ncurves
        stop
      end if

      lnnodes = mesh%curves(curve)%nnodes

    else if ( present(surface) ) then

      if ( surface < 1 .or. surface > mesh%nsurfaces ) then
        write(*,'(/a/2(a,i0)/)') &
          'Error: surface in the heading of create_sample is ', &
          'out of range. surface is ', surface, &
          'whereas the number of surfaces is ', mesh%nsurfaces
        stop
      end if

      lnnodes = mesh%surfaces(surface)%nnodes

    else if ( present(nodes) ) then

      lnnodes = size(nodes)

    else if ( present(points) ) then

      lnnodes = size(points)

    else if ( present(nnodes) ) then

      if ( nnodes < 1 ) then
        write(*,'(/a,a/)') &
          'Error: nnodes in the heading of create_sample ', &
          'must be positive.'
        stop
      end if

      lnnodes = nnodes

    else

      write(*,'(2a/6(a/))') &
        ' Error create_sample, one of the following keywords must be ', &
        ' present: ', &
        '    object ', &
        '    nodes  ', &
        '    points  ', &
        '    curve  ', &
        '    surface  ', &
        '    nnodes '
      stop

    end if

    if ( sample%created ) then

!     test size

      if ( lnnodes > size(sample%u,1) .or. ndegfd > size(sample%u,2) .or. &
           mesh%ndim > size(sample%coor,2) ) then

!       release memory because size is too small to fit in the sample
        deallocate(sample%u)
        deallocate(sample%coor)

        allocate(sample%u(lnnodes,ndegfd))
        allocate(sample%coor(lnnodes,mesh%ndim))

      end if

      if ( lallnodes ) then

!       all nodes sampled

        if ( lnnodes > size(sample%valid) ) then

!         release memory because size is too small to fit in the sample
          deallocate(sample%valid)
          allocate(sample%valid(lnnodes))

        end if

      else

!       not all nodes sampled (set valid to zero size)

        deallocate(sample%valid)
        allocate(sample%valid(0))

      end if

    else

      allocate(sample%u(lnnodes,ndegfd))
      allocate(sample%coor(lnnodes,mesh%ndim))
      if ( lallnodes ) then
        allocate(sample%valid(lnnodes))
      else
        allocate(sample%valid(0))
      end if

      sample%created = .true.

    end if

!   create empty sample

    sample%nnodes = 0
    sample%ndim = 0
    sample%ndegfd = 0

  end subroutine create_sample


! Delete sample

  subroutine delete_sample ( sample1, sample2, sample3, sample4, sample5 )

    type(sample_t), intent(inout) :: sample1
    type(sample_t), optional, intent(inout) :: sample2, sample3, sample4, &
      sample5

    call delete_single_sample(sample1,1)
    if ( present(sample2) ) call delete_single_sample(sample2,2)
    if ( present(sample3) ) call delete_single_sample(sample3,3)
    if ( present(sample4) ) call delete_single_sample(sample4,4)
    if ( present(sample5) ) call delete_single_sample(sample5,5)

  end subroutine delete_sample


! delete single sample

  subroutine delete_single_sample ( sample, nr )

    type(sample_t), intent(inout) :: sample
    integer, intent(in) :: nr

    if ( .not. sample%created ) then
      write(*,'(/2(a/),a,i0/)') &
        'Error in delete_sample:', &
        ' sample has not been created and cannot be deleted ', &
        ' sample number in heading = ', nr
      stop
    end if

    deallocate(sample%u)
    deallocate(sample%coor)
    deallocate(sample%valid)

    sample%nnodes = 0
    sample%ndim = 0
    sample%ndegfd = 0

    sample%created = .false.

  end subroutine delete_single_sample


! sample nodes

  subroutine sample_nodes ( mesh, problem, sample, nodes, physq, layer, degfd, &
    sysvector, vector, allnodes )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

    type(sample_t), intent(inout) :: sample

    integer, dimension(:), intent(in) :: nodes

!   the physical quantity that is being sampled in sysvector
    integer, intent(in), optional :: physq

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered.
    integer, intent(in), optional :: layer

!   the degrees of freedom that are being sampled. If physq is
!   present, degfd means the degrees of freedom within physq.
!   Example: degfd=(/1,0,1/), means the first and third degree of freedom.
!   If a degree of freedom is included that does not exist, it is simply
!   ignored.
!   In a node where the valid number of degrees of freedom is not sample%ndegfd
!   the node is not included in the sample.
    integer, dimension(:), intent(in), optional :: degfd

!   The vectors that are being sampled. Note, that either sysvector or vector
!   must be present.
    type(sysvector_t), intent(in), optional :: sysvector
    type(vector_t), intent(in), optional :: vector

!   if present and .true. samples of "all nodes" are defined (no nodes are
!   skipped). An additional logical array valid in the sample structure is
!   filled that indicates whether the corresponding sample value is valid
!   or not.
!   default=.false.
    logical, intent(in), optional :: allnodes

!   sample nodes


    logical :: lallnodes
    integer :: ndim, nnodes, ndegfd

    lallnodes = set_optional ( variable=allnodes, default=.false. )

    if ( any( nodes < 1 ) .or. any( nodes > mesh%nnodes ) ) then
      write(*,'(/a/a/)') &
        'Error: nodes in the heading of sample_nodes is ', &
        'out of range.'
      stop
    end if

!   take maximum possible first: all nodes
    nnodes = size(nodes)

    if ( nnodes > size(sample%u,1) .or. nnodes > size(sample%coor,1) ) then
      write(*,'(/a/3(a,i0/)/)') &
        'Error: sample in the heading of sample_nodes is ', &
        'too small to store the data: number of nodes needed = ', nnodes, &
        'whereas size(sample%u,1) = ', size(sample%u,1), &
        'and size(sample%coor,1) = ', size(sample%coor,1)
      stop
    end if

    if ( lallnodes .and. nnodes > size(sample%valid) ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: sample in the heading of sample_nodes is ', &
        'too small to store the data: number of nodes needed = ', nnodes, &
        'whereas size(sample%valid) = ', size(sample%valid)
      stop
    end if

    ndegfd = sample%ndegfd

    if ( ndegfd > size(sample%u,2) ) then
      write(*,'(/a/2(a,i0/)/)') &
        'Error: sample in the heading of sample_nodes is ', &
        'too small to store the data: number of degrees of freedom needed = ', &
        ndegfd, &
        'whereas size(sample%u,2) = ', size(sample%u,2)
      stop
    end if

    ndim = mesh%ndim

    if ( ndim > size(sample%coor,2) ) then
      write(*,'(/a/2(a,i0/)/)') &
        'Error: sample in the heading of sample_nodes is ', &
        'too small to store the data: space dimension in mesh is = ', ndim, &
        'and size(sample%coor,2) = ', size(sample%coor,2)
      stop
    end if

!   take sample

    sample%ndim = ndim

    call sample_in_nodes ( mesh, problem, sample, nodes, physq, layer, degfd, &
      sysvector, vector, allnodes )

  end subroutine sample_nodes


! sample points

  subroutine sample_points ( mesh, problem, sample, points, physq, layer, &
    degfd, sysvector, vector, allnodes )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

    type(sample_t), intent(inout) :: sample

    integer, dimension(:), intent(in) :: points

!   the physical quantity that is being sampled in sysvector
    integer, intent(in), optional :: physq

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered.
    integer, intent(in), optional :: layer

!   the degrees of freedom that are being sampled. If physq is
!   present, degfd means the degrees of freedom within physq.
!   Example: degfd=(/1,0,1/), means the first and third degree of freedom.
!   If a degree of freedom is included that does not exist, it is simply
!   ignored.
!   In a node where the valid number of degrees of freedom is not sample%ndegfd
!   the node is not included in the sample.
    integer, dimension(:), intent(in), optional :: degfd

!   The vectors that are being sampled. Note, that either sysvector or vector
!   must be present.
    type(sysvector_t), intent(in), optional :: sysvector
    type(vector_t), intent(in), optional :: vector

!   if present and .true. samples of "all nodes" are defined (no nodes are
!   skipped). An additional logical array valid in the sample structure is
!   filled that indicates whether the corresponding sample value is valid
!   or not.
!   default=.false.
    logical, intent(in), optional :: allnodes

!   sample points


    logical :: lallnodes
    integer :: ndim, nnodes, ndegfd

    lallnodes = set_optional ( variable=allnodes, default=.false. )

    if ( any( points < 1 ) .or. any( points > mesh%npoints ) ) then
      write(*,'(/a/a/)') &
        'Error: points in the heading of sample_points is ', &
        'out of range.'
      stop
    end if

!   take maximum possible first: all nodes
    nnodes = size(points)

    if ( nnodes > size(sample%u,1) .or. nnodes > size(sample%coor,1) ) then
      write(*,'(/a/3(a,i0/)/)') &
        'Error: sample in the heading of sample_points is ', &
        'too small to store the data: number of nodes needed = ', nnodes, &
        'whereas size(sample%u,1) = ', size(sample%u,1), &
        'and size(sample%coor,1) = ', size(sample%coor,1)
      stop
    end if

    if ( lallnodes .and. nnodes > size(sample%valid) ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: sample in the heading of sample_points is ', &
        'too small to store the data: number of nodes needed = ', nnodes, &
        'whereas size(sample%valid) = ', size(sample%valid)
      stop
    end if

    ndegfd = sample%ndegfd

    if ( ndegfd > size(sample%u,2) ) then
      write(*,'(/a/2(a,i0/)/)') &
        'Error: sample in the heading of sample_points is ', &
        'too small to store the data: number of degrees of freedom needed = ', &
        ndegfd, &
        'whereas size(sample%u,2) = ', size(sample%u,2)
      stop
    end if

    ndim = mesh%ndim

    if ( ndim > size(sample%coor,2) ) then
      write(*,'(/a/2(a,i0/)/)') &
        'Error: sample in the heading of sample_points is ', &
        'too small to store the data: space dimension in mesh is = ', ndim, &
        'and size(sample%coor,2) = ', size(sample%coor,2)
      stop
    end if

!   take sample

    sample%ndim = ndim

    call sample_in_nodes ( mesh, problem, sample, mesh%points(points), physq, &
      layer, degfd, sysvector, vector, allnodes )

  end subroutine sample_points


! sample a curve

  subroutine sample_curve ( mesh, problem, sample, curve, physq, layer, degfd, &
    sysvector, vector, allnodes )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

    type(sample_t), intent(inout) :: sample

    integer, intent(in) :: curve

!   the physical quantity that is being sampled in sysvector
    integer, intent(in), optional :: physq

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered.
    integer, intent(in), optional :: layer

!   the degrees of freedom that are being sampled. If physq is
!   present, degfd means the degrees of freedom within physq.
!   Example: degfd=(/1,0,1/), means the first and third degree of freedom.
!   If a degree of freedom is included that does not exist, it is simply
!   ignored.
!   In a node where the valid number of degrees of freedom is not sample%ndegfd
!   the node is not included in the sample.
    integer, dimension(:), intent(in), optional :: degfd

!   The vectors that are being sampled. Note, that either sysvector or vector
!   must be present.
    type(sysvector_t), intent(in), optional :: sysvector
    type(vector_t), intent(in), optional :: vector

!   if present and .true. samples of "all nodes" are defined (no nodes are
!   skipped). An additional logical array valid in the sample structure is
!   filled that indicates whether the corresponding sample value is valid
!   or not.
!   default=.false.
    logical, intent(in), optional :: allnodes

!   sample a curve


    logical :: lallnodes
    integer :: ndim, nnodes, ndegfd

    lallnodes = set_optional ( variable=allnodes, default=.false. )

    if ( curve < 1 .or. curve > mesh%ncurves ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: curve in the heading of sample_curve is ', &
        'out of range. curve is ', curve, &
        'whereas the number of curves is ', mesh%ncurves
      stop
    end if

!   take maximum possible first: all nodes
    nnodes = mesh%curves(curve)%nnodes

    if ( nnodes > size(sample%u,1) .or. nnodes > size(sample%coor,1) ) then
      write(*,'(/a/3(a,i0)/)') &
        'Error: sample in the heading of sample_curve is ', &
        'too small to store the data: number of nodes needed = ', nnodes, &
        'whereas size(sample%u,1) = ', size(sample%u,1), &
        'and size(sample%coor,1) = ', size(sample%coor,1)
      stop
    end if

    if ( lallnodes .and. nnodes > size(sample%valid) ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: sample in the heading of sample_curve is ', &
        'too small to store the data: number of nodes needed = ', nnodes, &
        'whereas size(sample%valid) = ', size(sample%valid)
      stop
    end if

    ndegfd = sample%ndegfd

    if ( ndegfd > size(sample%u,2) ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: sample in the heading of sample_curve is ', &
        'too small to store the data: number of degrees of freedom needed = ', &
        ndegfd, &
        'whereas size(sample%u,2) = ', size(sample%u,2)
      stop
    end if

    ndim = mesh%curves(curve)%ndim

    if ( ndim > size(sample%coor,2) ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: sample in the heading of sample_curve is ', &
        'too small to store the data: space dimension in curve is = ', ndim, &
        'and size(sample%coor,2) = ', size(sample%coor,2)
      stop
    end if

!   take sample

    sample%ndim = ndim

!     sample sysvector

    call sample_in_nodes ( mesh, problem, sample, mesh%curves(curve)%nodes, &
      physq, layer, degfd, sysvector, vector, allnodes )

  end subroutine sample_curve


! sample a surface

  subroutine sample_surface ( mesh, problem, sample, surface, physq, layer, &
    degfd, sysvector, vector, allnodes )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

    type(sample_t), intent(inout) :: sample

    integer, intent(in) :: surface

!   the physical quantity that is being sampled in sysvector
    integer, intent(in), optional :: physq

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered.
    integer, intent(in), optional :: layer

!   the degrees of freedom that are being sampled. If physq is
!   present, degfd means the degrees of freedom within physq.
!   Example: degfd=(/1,0,1/), means the first and third degree of freedom.
!   If a degree of freedom is included that does not exist, it is simply
!   ignored.
!   In a node where the valid number of degrees of freedom is not sample%ndegfd
!   the node is not included in the sample.
    integer, dimension(:), intent(in), optional :: degfd

!   The vectors that are being sampled. Note, that either sysvector or vector
!   must be present.
    type(sysvector_t), intent(in), optional :: sysvector
    type(vector_t), intent(in), optional :: vector

!   if present and .true. samples of "all nodes" are defined (no nodes are
!   skipped). An additional logical array valid in the sample structure is
!   filled that indicates whether the corresponding sample value is valid
!   or not.
!   default=.false.
    logical, intent(in), optional :: allnodes

!   sample a surface


    logical :: lallnodes
    integer :: ndim, nnodes, ndegfd

    lallnodes = set_optional ( variable=allnodes, default=.false. )

    if ( surface < 1 .or. surface > mesh%nsurfaces ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: surface in the heading of sample_surface is ', &
        'out of range. surface is ', surface, &
        'whereas the number of surface is ', mesh%nsurfaces
      stop
    end if

!   take maximum possible first: all nodes
    nnodes = mesh%surfaces(surface)%nnodes

    if ( nnodes > size(sample%u,1) .or. nnodes > size(sample%coor,1) ) then
      write(*,'(/a/3(a,i0)/)') &
        'Error: sample in the heading of sample_surface is ', &
        'too small to store the data: number of nodes needed = ', nnodes, &
        'whereas size(sample%u,1) = ', size(sample%u,1), &
        'and size(sample%coor,1) = ', size(sample%coor,1)
      stop
    end if

    if ( lallnodes .and. nnodes > size(sample%valid) ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: sample in the heading of sample_surface is ', &
        'too small to store the data: number of nodes needed = ', nnodes, &
        'whereas size(sample%valid) = ', size(sample%valid)
      stop
    end if

    ndegfd = sample%ndegfd

    if ( ndegfd > size(sample%u,2) ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: sample in the heading of sample_surface is ', &
        'too small to store the data: number of degrees of freedom needed = ', &
        ndegfd, &
        'whereas size(sample%u,2) = ', size(sample%u,2)
      stop
    end if

    ndim = mesh%surfaces(surface)%ndim

    if ( ndim > size(sample%coor,2) ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: sample in the heading of sample_surface is ', &
        'too small to store the data: space dimension in surface is = ', ndim, &
        'and size(sample%coor,2) = ', size(sample%coor,2)
      stop
    end if

!   take sample

    sample%ndim = ndim

    call sample_in_nodes ( mesh, problem, sample, &
      mesh%surfaces(surface)%nodes, physq, layer, degfd, sysvector, vector, &
      allnodes )

  end subroutine sample_surface


! sample in actual nodes (help routine for nodes given by nodes, points, curves
! and surfaces).

  subroutine sample_in_nodes ( mesh, problem, sample, nodes, physq, layer, &
    degfd, sysvector, vector, allnodes )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

    type(sample_t), intent(inout) :: sample

    integer, dimension(:), intent(in) :: nodes

!   the physical quantity that is being sampled in sysvector
    integer, intent(in), optional :: physq

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered.
    integer, intent(in), optional :: layer

!   the degrees of freedom that are being sampled. If physq is
!   present, degfd means the degrees of freedom within physq.
!   Example: degfd=(/1,0,1/), means the first and third degree of freedom.
!   If a degree of freedom is included that does not exist, it is simply
!   ignored.
!   In a node where the valid number of degrees of freedom is not sample%ndegfd
!   the node is not included in the sample.
    integer, dimension(:), intent(in), optional :: degfd

!   The vectors that are being sampled. Note, that either sysvector or vector
!   must be present.
    type(sysvector_t), intent(in), optional :: sysvector
    type(vector_t), intent(in), optional :: vector

!   if present and .true. samples of "all nodes" are defined (no nodes are
!   skipped). An additional logical array valid in the sample structure is
!   filled that indicates whether the corresponding sample value is valid
!   or not.
!   default=.false.
    logical, intent(in), optional :: allnodes


!   sample nodes


    logical :: lallnodes, lvalid
    integer :: ndim, nnodes, ndegfd, node, nodep, nodenr
    integer :: dof, ndeg, ideg, i
    integer, dimension(sample%ndegfd) :: pos, ind

    lallnodes = set_optional ( variable=allnodes, default=.false. )

!   take maximum possible first: all nodes

    nnodes = size(nodes)
    ndegfd = sample%ndegfd
    ndim = sample%ndim

!   take sample

    if ( present(sysvector) ) then

!     sample sysvector

      node = 0

      do nodep = 1, nnodes

        lvalid = .true.

!       get positions of degrees of freedom

        nodenr = nodes(nodep)

        if ( present(physq) ) then
          call pos_array_node ( problem, nodenr, dof, pos, [physq], layer )
        else
          call pos_array_node ( problem, nodenr, dof, pos, layer=layer )
        end if

!       fill array ind, pointing to the degrees included

        if ( present(degfd) ) then
          ndeg = min( dof, size(degfd) )
          if ( count( degfd(1:ndeg) > 0 ) /= ndegfd ) then
            if ( lallnodes ) then
              lvalid = .false.  ! invalid sample
            else
              cycle ! ignore this node
            end if
          end if
          ideg = 0
          do i = 1, ndeg
            if ( degfd(i) < 1 ) cycle  ! do not include this degree
            ideg = ideg + 1
            ind(ideg) = i
          end do
        else
          if ( dof < ndegfd ) then
            if ( lallnodes ) then
              lvalid = .false.  ! invalid sample
            else
              cycle ! ignore this node
            end if
          end if
          ind = [ ( i, i=1,ndegfd ) ]
        end if

        node = node + 1

!       fill node

        if ( lvalid ) then

!         valid sample

          sample%u(node,1:ndegfd) = sysvector%u(pos(ind))

        else

!         invalid sample

          sample%u(node,1:ndegfd) = 0 ! set to zero

        end if

        sample%coor(node,1:ndim) = mesh%coor(nodenr,1:ndim)

        if ( lallnodes ) sample%valid(node) = lvalid

      end do

      sample%nnodes = node

    else if ( present(vector) ) then

!     sample vector

      node = 0

      do nodep = 1, nnodes

        lvalid = .true.

!       get positions of degrees of freedom

        nodenr = nodes(nodep)

        call pos_array_vec_node ( problem, nodenr, dof, pos, vector%vec, layer )

!       fill array ind, pointing to the degrees included

        if ( present(degfd) ) then
          ndeg = min( dof, size(degfd) )
          if ( count( degfd(1:ndeg) > 0 ) /= ndegfd ) then
            if ( lallnodes ) then
              lvalid = .false.  ! invalid sample
            else
              cycle ! ignore this node
            end if
          end if
          ideg = 0
          do i = 1, ndeg
            if ( degfd(i) < 1 ) cycle  ! do not include this degree
            ideg = ideg + 1
            ind(ideg) = i
          end do
        else
          if ( dof < ndegfd ) then
            if ( lallnodes ) then
              lvalid = .false.  ! invalid sample
            else
              cycle ! ignore this node
            end if
          end if
          ind = [ ( i, i=1,ndegfd ) ]
        end if

        node = node + 1

!       fill node
        if ( lvalid ) then

!         valid sample

          sample%u(node,1:ndegfd) = vector%u(pos(ind))

        else

!         invalid sample

          sample%u(node,1:ndegfd) = 0 ! set to zero

        end if

        sample%coor(node,1:ndim) = mesh%coor(nodenr,1:ndim)

        if ( lallnodes ) sample%valid(node) = lvalid

      end do

      sample%nnodes = node

    end if

  end subroutine sample_in_nodes


! sample an object

  subroutine sample_object ( mesh, problem, sample, object, elemsub, &
    oldvectors, coefficients, mcoefficients, layer, errorobject, allnodes )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

    integer, intent(in) :: object

    type(sample_t), intent(inout) :: sample

    interface
      subroutine elemsub ( mesh, problem, object, nodeobj, coefficients, &
        oldvectors, u )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_m, only: problem_t
        use element_defs_m, only: oldvectors_t, coefficients_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: object, nodeobj
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:) :: u
      end subroutine elemsub
    end interface

!   coefficients is passed through to the element subroutine
    type(coefficients_t), intent(in), optional :: coefficients

!   array of coefficients (length=number of element groups)
!   If present only one of them is passed through to the element subroutine
!   based on the element group number
    type(coefficients_t), dimension(:), intent(in), optional :: mcoefficients

!   This structure can be used to supply old sysvectors/vectors from possibly
!   other problems to the element subroutines. Note that the problem structures
!   of the other problems must be included in oldvectors.
    type(oldvectors_t), intent(in) :: oldvectors

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered.
!   Only elements where all nodes have degrees in the specified layer are
!   taken into account for the sampling.
    integer, intent(in), optional :: layer

!   error checking indicator for objects:
!     0 : ignore sample points that do not intersect with the mesh/layer
!     1 : warn for sample points that do not intersect with the mesh/layer
!     2 : stop if sample points do not intersect with the mesh/layer
!   NOTE: if allnodes=.true., for both errorobject=0 or 1 the array
!         sample%valid indicates whether a sample is valid or not.
!   default=0 if layer not present
!   default=2 if layer is present
    integer, intent(in), optional :: errorobject

!   if present and .true. samples of "all nodes" are defined (no nodes are
!   skipped). An additional logical array valid in the sample structure is
!   filled that indicates whether the corresponding sample value is valid
!   or not.
!   default=.false.
    logical, intent(in), optional :: allnodes


!   sample an object


    logical :: lallnodes, lvalid
    logical, allocatable, dimension(:) :: wkl
    integer :: ndim, nnodes, ndegfd, node, nodeobj, elgrp, lerror, elem
    integer :: lerrordef
    real(dp) :: u(sample%ndegfd)
    type(coefficients_t) :: coeffl

    allocate ( wkl(mesh%nnodes) )

    lallnodes = set_optional ( variable=allnodes, default=.false. )

    if ( object < 1 .or. object > mesh%nobjects ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: object in the heading of sample_object is ', &
        'out of range. object is ', object, &
        'whereas the number of objects is ', mesh%nobjects
      stop
    end if

    if ( lallnodes ) then
      nnodes = mesh%objects(object)%nnodes
    else
      nnodes = count( mesh%objects(object)%grpelm(:,1) > 0 )
    end if

    if ( nnodes > size(sample%u,1) .or. nnodes > size(sample%coor,1) ) then
      write(*,'(/a/3(a,i0)/)') &
        'Error: sample in the heading of sample_object is ', &
        'too small to store the data: number of nodes needed = ', nnodes, &
        'whereas size(sample%u,1) = ', size(sample%u,1), &
        'and size(sample%coor,1) = ', size(sample%coor,1)
      stop
    end if

    if ( lallnodes .and. nnodes > size(sample%valid) ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: sample in the heading of sample_object is ', &
        'too small to store the data: number of nodes needed = ', nnodes, &
        'whereas size(sample%valid) = ', size(sample%valid)
      stop
    end if

    ndegfd = sample%ndegfd

    if ( ndegfd > size(sample%u,2) ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: sample in the heading of sample_object is ', &
        'too small to store the data: number of degrees of freedom needed = ', &
        ndegfd, &
        'whereas size(sample%u,2) = ', size(sample%u,2)
      stop
    end if

    ndim = mesh%objects(object)%ndim

    if ( ndim > size(sample%coor,2) ) then
      write(*,'(/a/2(a,i0)/)') &
        'Error: sample in the heading of sample_object is ', &
        'too small to store the data: space dimension in object is = ', ndim, &
        'and size(sample%coor,2) = ', size(sample%coor,2)
      stop
    end if

    if ( present(layer) ) then
      wkl = btest(problem%nodlayers,layer-1)
    end if

!   coefficients

    if ( present(coefficients) ) coeffl = coefficients

!   error

    if ( present(layer) ) then
      lerrordef = 2  ! stop on sample points outside mesh or layer
    else
      lerrordef = 0  ! ignore sample points outside mesh
    end if

    lerror = set_optional ( variable=errorobject, default=lerrordef )

!   take sample

    node = 0

    do nodeobj = 1, mesh%objects(object)%nnodes

      lvalid = .true.

      elgrp = mesh%objects(object)%grpelm(nodeobj,1)

!     intersection with mesh?
      if ( elgrp == 0 ) then
        if ( lerror == 0 ) then
          if ( lallnodes ) then
            lvalid = .false.  ! invalid sample
          else
            cycle
          end if
        else if ( lerror == 1 ) then
          write(*,'(/2(a/),2(a,i0)/a,3es16.8/)') &
            'Warning in sample_object: ', &
            ' object node does not intersect with the mesh ', &
            ' object = ', object, ' node = ', nodeobj, &
            ' coordinates = ', mesh%objects(object)%coor(nodeobj,:)
          if ( lallnodes ) then
            lvalid = .false.  ! invalid sample
            write(*,'(/a/)') ' node set to invalid (valid=.false.) '
          else
            cycle
            write(*,'(/a/)') &
              ' node will be skipped and number of sample points reduced '
          end if
        else if ( lerror == 2 ) then
          write(*,'(/2(a/),2(a,i0)/a,3es16.8/)') &
            'Error in sample_object: ', &
            ' object node does not intersect with the mesh ', &
            ' object = ', object, ' node = ', nodeobj, &
            ' coordinates = ', mesh%objects(object)%coor(nodeobj,:)
          stop
        end if
      end if

!     intersection with layer?
      if ( present(layer) ) then
        elem = mesh%objects(object)%grpelm(nodeobj,2)
        if ( .not. all ( wkl(mesh%topology(elgrp)%a(:,elem)) ) ) then
!         not all nodes in layer
          if ( lerror == 0 ) then
            if ( lallnodes ) then
              lvalid = .false.  ! invalid sample
            else
              cycle
            end if
          else if ( lerror == 1 ) then
            write(*,'(/2(a/),3(a,i0)/a,3es16.8/)') &
              'Warning in sample_object: ', &
              ' object node does not intersect with the elements in layer ', &
              ' layer = ', layer, '  object = ', object, ' node = ', nodeobj, &
              ' coordinates = ', mesh%objects(object)%coor(nodeobj,:)
            if ( lallnodes ) then
              lvalid = .false.  ! invalid sample
              write(*,'(/a/)') ' node set to invalid (valid=.false.) '
            else
              write(*,'(/a/)') &
                ' node will be skipped and number of sample points reduced '
              cycle
            end if
          else if ( lerror == 2 ) then
            write(*,'(/2(a/),3(a,i0)/a,3es16.8/)') &
              'Error in sample_object: ', &
              ' object node does not intersect with the elements in layer ', &
              ' layer = ', layer, '  object = ', object, ' node = ', nodeobj, &
              ' coordinates = ', mesh%objects(object)%coor(nodeobj,:)
            stop
          end if
        end if
      end if

      if ( lvalid ) then

!       valid sample

        if ( present(mcoefficients) ) coeffl = mcoefficients(elgrp)

!       node subroutine
        call elemsub ( mesh, problem, object, nodeobj, coeffl, oldvectors, u )

      else

!       invalid sample

        u = 0  ! set to zero

      end if

      node = node + 1

!     fill node
      sample%u(node,1:ndegfd) = u
      sample%coor(node,1:ndim) = mesh%objects(object)%coor(nodeobj,:)
      if ( lallnodes ) sample%valid(node) = lvalid

    end do

    if ( .not. present(layer) .and. nnodes /= node ) then
      stop 'internal error sample_object'
    end if

    sample%nnodes = node
    sample%ndim = ndim

    deallocate ( wkl )

  end subroutine sample_object


! fill a sample on an object, curve ...

  subroutine fill_sample ( mesh, problem, sample, object, elemsub, ndegfd, &
    nodes, points, curve, surface, physq, layer, degfd, sysvector, &
    vector, coefficients, mcoefficients, oldvectors, errorobject, allnodes )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

    type(sample_t), intent(inout) :: sample

!   if present, the sample is filled from the object
    integer, intent(in), optional :: object

!   this is the user routine for filling the sample value at a single node
    optional :: elemsub
    interface
      subroutine elemsub ( mesh, problem, object, nodeobj, coefficients, &
        oldvectors, u )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_m, only: problem_t
        use element_defs_m, only: oldvectors_t, coefficients_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: object, nodeobj
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:) :: u
      end subroutine elemsub
    end interface

!   if present: the number of degrees of freedom in each point of the sample.
!   if not present: the number of degrees of freedom in each point of the
!                   sample is given by the maximum possible value as given
!                   by size(sample%u,2). NOTE: this is only possible if
!                   sample has already been created.
    integer, intent(in), optional :: ndegfd

!   if present, the sample is filled from the nodal points
    integer, dimension(:), intent(in), optional :: nodes

!   if present, the sample is filled from the points
    integer, dimension(:), intent(in), optional :: points

!   if present, the sample is filled from the curve
    integer, intent(in), optional :: curve

!   if present, the sample is filled from the surface
    integer, intent(in), optional :: surface

!   if present, the physical quantity that is being sampled in sysvector
    integer, intent(in), optional :: physq

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered.
!   NOTE: For objects this means that only elements where all nodes have
!   degrees in the specified layer are taken into account for the sampling.
    integer, intent(in), optional :: layer

!   if present, the degrees of freedom that are being sampled on a curve.
!   If physq is present, degfd means the degrees of freedom within physq.
!   Example: degfd=(/1,0,1/), means the first and third degree of freedom.
!   If a degree of freedom is included that does not exist, it is simply
!   ignored.
!   In a node where the valid number of degrees of freedom is not sample%ndegfd
!   the node is not included in the sample.
    integer, dimension(:), intent(in), optional :: degfd

!   The vectors that are being sampled on a curve.
!   Note, that either sysvector or vector must be present.
    type(sysvector_t), intent(in), optional :: sysvector
    type(vector_t), intent(in), optional :: vector

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

!   error checking indicator for objects:
!     0 : ignore sample points that are do not intersect with the mesh/layer
!     1 : warn for sample points that are do not intersect with the mesh/layer
!     2 : stop if sample points do not intersect with the mesh/layer
!   NOTE: if allnodes=.true., for both errorobject=0 or 1 the array
!         sample%valid indicates whether a sample is valid or not.
!   default=0 if layer not present
!   default=2 if layer is present
    integer, intent(in), optional :: errorobject

!   if present and .true. samples of "all nodes" are defined (no nodes are
!   skipped). An additional logical array valid in the sample structure is
!   defined/filled that indicates whether the corresponding sample value is
!   valid or not.
!   default=.false.
    logical, intent(in), optional :: allnodes


!   fill a sample on an object, curve ...


    type(oldvectors_t) :: oldvl
    integer :: lndegfd


    call check ( mesh, 'sample_points' )
    call check ( problem, 'sample_points', mesh )


    if ( .not. sample%created ) then
      if ( present(ndegfd) ) then
        lndegfd = ndegfd
      else
        write(*,'(/a/a/)') &
          'Error in fill_sample: if sample is not yet created the heading', &
          ' parameter ndegfd needs to be specified '
        stop
      end if
    else
      lndegfd = set_optional ( variable=ndegfd, default=size(sample%u,2) )
    end if

    if ( present(physq) .and. present(sysvector) ) then
      if ( physq <=0 .or. physq > problem%nphysq ) then
        write(*,'(/a/a,i0/)') &
          'Error: physq in heading of fill_sample has wrong value:', &
          'physq <=0 or physq > number of physical quantities = ', &
          problem%nphysq
        stop
      end if
    end if

    if ( present(mcoefficients) ) then
      if ( size(mcoefficients) /= mesh%nelgrp ) then
        write(*,'(/2(a/))') &
          'Error in fill_sample: ', &
          '   size of mcoefficients /= number of element groups'
        stop
      end if
    end if

    if ( present(layer) ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in fill_sample: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a),i0/)') &
          'Error in fill_sample: ',&
          ' layer < 1 or layer > number of layers = ', &
          problem%numlayers
        stop
      end if
    end if


!   create sample (or use old one if large enough)

    call create_sample ( mesh, sample, lndegfd, object, nodes, points, &
      curve, surface, allnodes=allnodes )

!   ndegfd in sample

    sample%ndegfd = lndegfd

!   choose geometrical entity

    if ( present(object) ) then

!     object

!     initialize

      if ( present(oldvectors) ) oldvl = oldvectors

      if ( .not. present(elemsub) ) then
        write(*,'(2a)') &
          'Error fill_sample: when sampling an object, elemsub must be ', &
          'specified as well.'
        stop
      end if

      if ( present(errorobject) ) then
        if ( errorobject < 0 .or. errorobject > 2 ) then
          write(*,'(/a/a/)') &
            'Error: errorobject in heading of fill_sample has wrong value:', &
            'errorobject < 0 or errorobject > 2 '
          stop
        end if
      end if

      call sample_object ( mesh, problem, sample, object, elemsub, &
        oldvl, coefficients, mcoefficients, layer, errorobject, allnodes )

    else if ( present(nodes) .or. present(points) .or. present(curve) &
             .or. present(surface) ) then

!     nodes, points, curve or surface

      if ( present(sysvector) .and. present(vector) ) then

        write(*,'(/a/)') &
          'Error in fill_sample: both sysvector and vector present '
        stop

      else if ( present(sysvector) ) then

!       sysvector

        if ( .not. sysvector%created ) then
          write(*,'(/a/)') &
            'Error in fill_sample: sysvector has not been created '
          stop
        end if

      else if ( present(vector) ) then

!       vector

        if ( .not. vector%created ) then
          write(*,'(/a/)') &
            'Error in fill_sample: vector has not been created '
          stop
        end if

      else

        write(*,'(a)') &
          'Error fill_sample: no sysvector or vector present'
        stop

      end if

      if ( present(nodes) ) then
        call sample_nodes ( mesh, problem, sample, nodes, physq, layer, &
          degfd, sysvector, vector, allnodes )
      else if ( present(points) ) then
        call sample_points ( mesh, problem, sample, points, physq, layer, &
          degfd, sysvector, vector, allnodes )
      else if ( present(curve) ) then
        call sample_curve ( mesh, problem, sample, curve, physq, layer, &
          degfd, sysvector, vector, allnodes )
      else if ( present(surface) ) then
        call sample_surface ( mesh, problem, sample, surface, physq, layer, &
          degfd, sysvector, vector, allnodes )
      end if

    else

      write(*,'(8(/a)/)') &
        'Error: could not determine type of sample from the heading', &
        'of fill_sample. None of the keywords:',      &
        '  object ', &
        '  nodes ', &
        '  points ', &
        '  curve ', &
        '  surface ', &
        'is present.'
      stop

    end if

  end subroutine fill_sample


! Integrate over mesh

  subroutine integrate ( mesh, problem, resultsum, elemsub, coefficients, &
    mcoefficients, oldvectors, elgroup1, elgroup2, groups, layer, elementset, &
    add )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the result of the integration
    real(dp), intent(inout), dimension(:) :: resultsum

!   this is the element subroutine that must be supplied by the calling routine
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

!   if these are present the integration takes place for element groups
!   elgroup1,...,elgroup2 only. If only elgroup1 is present one group is
!   assembled only.
    integer, intent(in), optional :: elgroup1, elgroup2

!   if present: the element groups to be used
!   For example groups=(/2,4/) will integrate on element groups 2 and 4.
!   Default: all groups
    integer, dimension(:), intent(in), optional :: groups

!   if layer is present only elements where all nodes have degrees in the
!   specified layer are taken into account for the integration.
    integer, intent(in), optional :: layer

!   if present the integration takes place only on the elements of the domain
!   that is defined by the elementset.
!   Note, that this argument can be combined with other arguments such as
!   groups, layer etc. It is just a restriction on the elements that will be
!   assembled.
    integer, intent(in), optional :: elementset

!   if add is set to .true. resultsum is not cleared before
!   the integration and thus the element integration results are added to the
!   input value. The default is .false. (clearing)
    logical, intent(in), optional :: add


    logical :: clearres, first, last
    logical, allocatable, dimension(:) :: wkl
    integer :: elem, elgrp, grp
    integer :: lgroups(mesh%nelgrp), lnelgrp, i
    real(dp) :: elemvec(size(resultsum))
    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl
    type(logical_array_1d_t), allocatable, dimension(:) :: integrateelements

    allocate ( wkl(mesh%nnodes) )

    call check ( mesh, 'integrate' )
    call check ( problem, 'integrate', mesh )

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
        ' integrate is out of range: ', &
        ' some groups are < 1 or larger than the number of groups ', mesh%nelgrp
      stop
    end if

!   layers

    if ( present(layer) ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in integrate: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a),i0/)') &
          'Error in integrate: ',&
          ' layer < 1 or layer > number of layers = ', &
          problem%numlayers
        stop
      end if
      wkl = btest(problem%nodlayers,layer-1)
    end if

!   test groups

    do grp = 1, lnelgrp

      elgrp = lgroups(grp)

!     skip inactive groups
      if ( any( elgrp == problem%inactivegroups ) ) cycle

      if ( all( problem%elnumdegfd(elgrp)%a == 0 ) ) then

        write(*,'(/a/a,i0/)') &
          'Error: element group in the heading of integrate has ', &
          'no degrees of freedom. elgrp is ', elgrp
        stop

      end if

    end do

!   integrate on elementset

    if ( present(elementset) ) then

!     check elementset

      if ( elementset <=0 .or. elementset > mesh%nelementsets ) then
        write(*,'(/a/a,i0/)') &
          'Error: elementset in heading of integrate has wrong value:', &
          'elementset <=0 or elementset > number elementsets = ', &
          mesh%nelementsets
        stop
      end if

!     create logical array for elements

      allocate ( integrateelements(mesh%nelgrp) )

      do elgrp = 1, mesh%nelgrp

        allocate ( integrateelements(elgrp)%a(mesh%grpnumel(elgrp)) )

        integrateelements(elgrp)%a = .false.

        integrateelements(elgrp)%a(&
           mesh%elementsets(elementset)%elements(elgrp)%a ) = .true.

      end do

    end if


    if ( present(add) ) then
      clearres = .not. add
    else
      clearres = .true.
    end if

    if ( clearres ) resultsum = 0


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

!       skip elements not in layer
        if ( present(layer) ) then
          if ( .not. all ( wkl(mesh%topology(elgrp)%a(:,elem)) ) ) cycle
!         avoid leaving allocated memory in elements
          first = .true.; last = .true.
        end if

!       only elements on elementset
        if ( present(elementset) ) then
          if ( .not. integrateelements(elgrp)%a(elem) ) cycle
!         avoid leaving allocated memory in elements
          first = .true.; last = .true.
        end if

!       compute element vector

        call elemsub ( mesh, problem, elgrp, elem, first, last, coeffl, &
          oldvl, elemvec )

!       add element vector to resultsum

        resultsum = resultsum + elemvec

      end do

    end do

    deallocate ( wkl )

  end subroutine integrate


! Integrate boundary elements (on a geometry)

  subroutine integrate_boundary_elements ( mesh, problem, resultsum, elemsub, &
    curve, surface, coefficients, oldvectors, layer, add )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the result of the integration
    real(dp), intent(inout), dimension(:) :: resultsum

!   this is the element subroutine that must be supplied by the calling routine
    interface
      subroutine elemsub ( mesh, problem, geom, elem, first, last, &
        coefficients, oldvectors, elemvec )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: geom, elem
        logical, intent(in) :: first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:) :: elemvec
      end subroutine elemsub
    end interface

!   the integration on the boundary elements are for this curve only.
    integer, intent(in), optional :: curve

!   the integration of the boundary elements are for this surface only.
    integer, intent(in), optional :: surface

!   if coefficients is present, it is passed through to the element subroutine
    type(coefficients_t), intent(in), optional :: coefficients

!   if oldvectors is present, it is passed through to the element subroutine.
!   This structure can be used to supply old sysvectors/vectors from possibly
!   other problems to the element subroutines. Note that the problem structures
!   of the other problems must be included in oldvectors.
    type(oldvectors_t), intent(in), optional :: oldvectors

!   if layer is present only elements where all nodes have degrees in the
!   specified layer are taken into account for the integration.
    integer, intent(in), optional :: layer

!   if add is set to .true. resultsum is not cleared before
!   the integration and thus the element integration results are added to the
!   input value. The default is .false. (clearing)
    logical, intent(in), optional :: add



    logical :: clearres, first, last
    logical, allocatable, dimension(:) :: wkl
    integer :: elem, geom, nelem
    real(dp), dimension(size(resultsum)) :: elemvec
    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl
    type(geometry_t) :: geometry

    allocate ( wkl(mesh%nnodes) )

    call check ( mesh, 'integrate_boundary_elements' )
    call check ( problem, 'integrate_boundary_elements', mesh )

!   test valid geometry

    if ( present(curve) .and. present(surface) ) then
      write(*,'(/a/a/)') &
        'Error: curve and surface cannot be both present in the heading', &
        ' of integrate_boundary_elements'
      stop
    end if

    if ( present(curve) ) then

      if ( curve < 1 .or. curve > mesh%ncurves ) then

        write(*,'(/a/a,i0,/a,i0/)') &
          'Error: curve in the heading of integrate_boundary_elements is ', &
          'out of range. curve is ', curve, &
          'whereas the number of curves is ', mesh%ncurves
        stop

      end if

    else if ( present(surface) ) then

      if ( surface < 1 .or. surface > mesh%nsurfaces ) then

        write(*,'(/a/a,i0,/a,i0/)') &
          'Error: surface in the heading of integrate_boundary_elements is ', &
          'out of range. surface is ', surface, &
          'whereas the number of surfaces is ', mesh%nsurfaces
        stop

      end if

    else

      write(*,'(/a/a/)') &
        'Error: either curve or surface must be present in the heading', &
        ' of integrate_boundary_elements'
      stop

    end if

!   layers

    if ( present(layer) ) then
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in integrate_boundary_elements: layer out of range.'
        stop
      end if
      wkl = btest(problem%nodlayers,layer-1)
    end if

    if ( present(add) ) then
      clearres = .not. add
    else
      clearres = .true.
    end if

    if ( clearres ) resultsum = 0

!   initialize

    if ( present(oldvectors) ) oldvl = oldvectors
    if ( present(coefficients) ) coeffl = coefficients

!   find geometry and number of degr of freedom in an element of this geometry

    if ( present(curve) ) then
      geom = curve
      nelem = mesh%curves(curve)%nelem
      geometry = mesh%curves(curve)
    else if ( present(surface) ) then
      geom = surface
      nelem = mesh%surfaces(surface)%nelem
      geometry = mesh%surfaces(surface)
    end if

!   loop over elements in this geometry

    do elem = 1, nelem

      first = elem == 1
      last  = elem == nelem

!     skip elements not in layer
      if ( present(layer) ) then
        if ( .not. all( wkl(geometry%topology(:,elem,2)) ) ) cycle
!       avoid leaving allocated memory in elements
        first = .true.; last = .true.
      end if

!     compute element vector

      call elemsub ( mesh, problem, geom, elem, first, last, coeffl, oldvl, &
        elemvec )

!     add element vector to integration result

      resultsum = resultsum + elemvec

    end do

    deallocate ( wkl )

  end subroutine integrate_boundary_elements


! Integrate elements on an object

  subroutine integrate_object ( mesh, problem, resultsum, elemsub1, object, &
    coefficients, mcoefficients, oldvectors, layer, add )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the result of the integration
    real(dp), intent(inout), dimension(:) :: resultsum

!   this is the element subroutine that must be supplied by the calling routine
    interface
      subroutine elemsub1 ( mesh, problem, eleminfo, first, last, &
        coefficients, oldvectors, elemvec )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        type(eleminfo_t), intent(in) :: eleminfo
        logical, intent(in) :: first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:) :: elemvec
      end subroutine elemsub1
    end interface

!   integration on the elements for this object
    integer, intent(in) :: object

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

!   if layer is present only integration points in elements where all
!   nodes have degrees in the specified layer are taken into account for the
!   integration.
    integer, intent(in), optional :: layer

!   if add is set to .true. resultsum is not cleared before
!   the integration and thus the element integration results are added to the
!   input value. The default is .false. (clearing)
    logical, intent(in), optional :: add



    logical :: clearres, first, last
    logical, dimension(mesh%nnodes) :: wkl
    integer :: elem, nelem, elemo, i, elgrp
    real(dp), dimension(size(resultsum)) :: elemvec
    integer, dimension(mesh%nelgrp,maxval(mesh%grpnumel)) :: wk1, wk2
    integer, dimension(mesh%nelem,2) :: nwk
    integer :: nwkelem
    type(int_array_1d_p), allocatable, dimension(:) :: intps
    type(oldvectors_t) :: oldvl
    type(coefficients_t) :: coeffl
    type(eleminfo_t) :: eleminfo


    call check ( mesh, 'integrate_object' )
    call check ( problem, 'integrate_object', mesh )

!   test valid object

    if ( object <=0 .or. object > mesh%nobjects ) then
      write(*,'(/a/a,i0/)') &
        'Error: object in heading of integrate_object has wrong value:', &
        'object <=0 or object > number objects = ', mesh%nobjects
      stop
    end if

    if ( .not. mesh%objects(object)%intpoints ) then
      write(*,'(/a/a,i0,a/)') &
        'Error in integrate_object:', &
        '  object ', object, ' has no integration points '
      stop
    end if

    if ( any ( mesh%objects(object)%grpelm_int(:,1,:) == 0 ) ) then
!     not all elements connected
      write(*,'(2(/a)/a,i0/)') &
      'Error in integrate_object: ', &
      ' reference coordinates in object are missing,', &
      ' object = ', object
      stop
    end if

    if ( mesh%objects(object)%typeofobject == 2 ) then
      if ( any ( mesh%objects(object)%grpelm2_int(:,1,:) == 0 ) ) then
!       not all element connected (second intersection)
        write(*,'(2(/a)/a,i0/)') &
        'Error in integrate_object: ', &
        ' reference coordinates in second intersection of object are missing,',&
        ' object = ', object
        stop
      end if
    end if

!   layers

    if ( present(layer) ) then
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a)/)') &
          'Error in integrate_object: layer out of range.'
        stop
      end if
      wkl = btest(problem%nodlayers,layer-1)
    end if

    if ( present(add) ) then
      clearres = .not. add
    else
      clearres = .true.
    end if

    if ( clearres ) resultsum = 0

!   initialize

    if ( present(oldvectors) ) oldvl = oldvectors
    if ( present(coefficients) ) coeffl = coefficients

    nelem = mesh%objects(object)%nelem

!   set to zero work arrays for storing counted elements

    wk1 = 0
    wk2 = 0

!   loop over elements of this object

    do elemo = 1, nelem

!     find mesh elements and integration points

      call find_elements_integrationpoints ( mesh, object, elemo, nwkelem, &
        wk1, wk2, nwk, intps )

      do i = 1, nwkelem

        elgrp = nwk(i,1)
        elem  = nwk(i,2)

!       skip elements not in layer
        if ( present(layer) ) then
          if ( .not. all( wkl(mesh%topology(elgrp)%a(:,elem)) ) ) cycle
        end if

        if ( present(mcoefficients) ) coeffl = mcoefficients(elgrp)

        first = i == 1
        last  = i == nwkelem

        eleminfo%elgrp  = elgrp
        eleminfo%elem   = elem
        eleminfo%object = object
        eleminfo%elemo  = elemo
        eleminfo%nintps = size(intps(i)%a)
        eleminfo%intps  => intps(i)%a

!       compute element vector

        call elemsub1 ( mesh, problem, eleminfo, first, last, &
          coeffl, oldvl, elemvec )

!       add element vector to integration result

        resultsum = resultsum + elemvec

      end do

      call reset_elements_integrationpoints ( nwkelem, wk1, wk2, nwk, intps )

    end do

  end subroutine integrate_object

end module postprocessing_m
