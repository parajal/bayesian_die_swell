! Example program for using the eltree module for tetrahedra.
! Use eltree to subdivide the region and integrate over the volume of a
! sphere.

program eltree10

  use math_defs_m
  use io_utils_m
  use meshgen_m
  use eltree_m
  use functions_eltree_m
  use gauss_m
  use subs10_m
  use postprocessing_m

  implicit none

! constants

  integer, parameter :: &
    numsplit    = 4, &  ! relative size of all subelements will be at least
                        ! as small as 1/2^numsplitmin
    numsplitmin = 1, &  ! relative size of all subelements will be at least
                        ! as small as 1/2^numsplitmin
    uintpl      = 2, &  ! P1 interpolation
    gauss_ve    = 8     ! Gauss integration rule for the volume elements

  real(dp), parameter :: &
    split_threshold = 1.e-3, &  ! levelset value for being "close" enough
                                ! to the interface for tree splitting
    epsvol = 0.0_dp,  & ! elements with (reference) integration area
                        ! smaller than epsvol are removed from the eltree_array
    radius = 1.0_dp     ! the radius of the sphere/circle

! definitions

  type(mesh_t), target :: mesh
  type(problem_t) :: problem
  type(input_probdef_t) :: input_probdef
  type(coefficients_t) :: coefficients

! integration area near interface elements

  real(dp), dimension(:), allocatable :: vol

  integer :: i, j, nr_sub_el
  real(dp) :: resultsum(2), areae, intge

! pass the radius on to the levelset function in functions_eltree_m

  rpl = radius

! read mesh

  call read_mesh_gmsh ( mesh, filename='mesh10.msh' )

  call fill_mesh_parts ( mesh )

  call write_mesh_vtk ( mesh, filename='mesh10.vtk' )

  call printinfo ( mesh )

  lmesh => mesh  ! supply mesh to mapcoor for mapping coordinates
  lelem = 1
  lelgrp = 1

  allocate ( eltree_array(1,mesh%grpnumel(1)) )
  allocate ( vol(mesh%grpnumel(1)) )
  allocate ( d(mesh%nnodes) )

  d = levelset ( mesh%coor ) ! set levelset for all nodes

  call create_eltree_in_elements ( mesh, numsplit, numsplitmin, &
    split_threshold, d, vol, epsvol)

! find information about the element tree

  nr_sub_el = 0
  j = 0
  do i = 1, mesh%grpnumel(1)
    if ( associated( eltree_array(1,i)%p ) ) then
      nr_sub_el = nr_sub_el + number_of_subelements(eltree_array(1,i)%p)
      j = j + 1
    end if
  end do

  print *, 'number of elements containing the interface = ', j
  print *, 'number of subelements = ', nr_sub_el
  print *, ''

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl

  coefficients%r = 0

! problem definition

  call create_input_probdef ( mesh, input_probdef )

  input_probdef%elementdof(1)%a = 1

  call problem_definition ( input_probdef, mesh, problem )

! integrate on volume of the sphere

  intrule_ve = gauss_ve

  call integrate ( mesh, problem, resultsum, elemsub=integrate_volume, &
    coefficients=coefficients )

! analytical solutions

  areae = (4.0/3.0) * pi * radius**3
  intge = (16.0/12.0)*pi*(radius**4)/5

  print *, 'computed volume = ', resultsum(1)
  print *,  'relative error = ', (resultsum(1)-areae)/areae
  print *,  ''
  print *, 'computed integral of x^2 on whole domain = ', resultsum(2)
  print *,  'relative error = ', (resultsum(2)-intge)/intge

! post processing

  call post_processing

! delete the definitions

  call delete_eltree_in_elements
  deallocate( vol, d, eltree_array )

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( coefficients )

contains

   subroutine post_processing

    integer :: i

    type(mesh_t) :: mesh1, mesh2, mesh_plot

!   plot the interface mesh (double nodes removed)
    call mesh_convert ( mesh, mesh1, userelmesh=userelmesh1, warn=.false. )
    call fill_mesh_parts ( mesh1 )

    call add_to_mesh ( mesh1, nodblocks=[10,10,10] )
    call mesh_convert ( mesh1, mesh2, remove_double_nodes=.true., &
      warn=.false. )
    call fill_mesh_parts ( mesh2 )
    call write_mesh_vtk ( mesh2, filename='imesh10.vtk' )
    call delete ( mesh1, mesh2 )

    call copy( mesh, mesh_plot )

!   plot the eltrees
    do i = 1, mesh%grpnumel(1)
      if ( associated( eltree_array(1,i)%p )) then
        lelem=i

        call eltree_to_mesh ( eltree_array(1,i)%p, mesh1, lsign=[-1,0,1], &
          mapcoor=mapcoor_tet )

        call mesh_merge ( mesh_plot, mesh1, mesh2, warn=.false. )
        call delete ( mesh_plot )
        call copy ( mesh2, mesh_plot )
        call delete (mesh1, mesh2)

      end if
    end do

    call fill_mesh_parts ( mesh_plot )
    call write_mesh_vtk ( mesh_plot, filename='mesh_eltree10.vtk' )
    call delete(mesh_plot)

  end subroutine post_processing


! create eltree in the elements which are crossed by the interface

  subroutine create_eltree_in_elements ( mesh, numsplit, numsplitmin, &
    split_threshold, d, vol, epsvol )

    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: numsplit, numsplitmin
    real(dp), intent(in) :: split_threshold, epsvol
    real(dp), intent(in), dimension(:) :: d
    real(dp), dimension(:), intent(inout) :: vol

    real(dp), dimension(:,:), allocatable :: coor
    integer :: elem, conf, eltype, ndim, nod(mesh%element(1)%numnod)

    ndim = mesh%ndim

!   loop all elements

    do elem = 1, mesh%grpnumel(1)

!     nodal points

      nod = mesh%topology(1)%a(:,elem)

      if ( all ( d(nod) > split_threshold ) .or. &
                        all ( d(nod) < -split_threshold ) ) then

!       all nodes far from the interface
        cycle

      else

!       element possibly contains an interface, start subdivide

        allocate ( eltree_array(1,elem)%p )

        lelgrp = 1
        lelem = elem

!       fill root node of the eltree
!       NOTE: for tets, eltype = 1 and a configuration must be given
!             for the root element: conf = 2 and the coordinates:
!             coor(1) = (0,0) and coor(2) = (1,1)
!             for more info see eltree_t!
        allocate( coor(2,ndim) )
        coor(1,:) = [  0.0, 0.0, 0.0 ]
        coor(2,:) = [  1.0, 1.0, 1.0 ]
        conf = 2
        eltype = 1

        call fill_node_eltree ( eltree_array(1,elem)%p, coor, conf, &
          eltype=eltype )

        call subdivide ( eltree_array(1,elem)%p, levelset=levelset, &
          numsplit=numsplit, split_threshold=split_threshold, &
          numsplitmin=numsplitmin, mapcoor=mapcoor_tet, &
          submesh=.true., intmesh=.true. )

!       check whether interface is in element
        if ( any( number_of_subelements_vector(eltree_array(1,elem)%p,&
                         &lsign=[-1,1]) == 0 ) ) then

!         no interface, remove eltree
          call delete(eltree_array(1,elem)%p)  ! delete actual eltree
          deallocate( eltree_array(1,elem)%p ) ! delete pointer target
                                               ! pointer becomes disassociated
          print *, 'eltree removed (no interface in element), elem = ', elem

        end if

      end if

!     check whether the volume of the eltree is large enough
      if ( associated(eltree_array(1,elem)%p) ) then

!       element crosses the interface (tet with ref vol 1/6).
        vol(elem) = volume ( eltree_array(1,elem)%p, lsign=[1] ) * 6

        if ( vol(elem) <= epsvol ) then

!         remove eltree
          call delete(eltree_array(1,elem)%p)  ! delete actual eltree
          deallocate( eltree_array(1,elem)%p ) ! delete pointer target

          print *, 'element at interface removed, elem = ', elem
          print *, 'vol(elem) =',vol(elem)

        end if

      end if

    deallocate ( coor )

    end do

  end subroutine create_eltree_in_elements

! delete eltree in the elements

  subroutine delete_eltree_in_elements

  integer :: elem

!   loop all elements
    do elem = 1, size(eltree_array,2)

      if ( associated(eltree_array(1,elem)%p) ) then

        call delete(eltree_array(1,elem)%p) ! delete actual eltree
        deallocate(eltree_array(1,elem)%p)  ! delete pointer target
                                            ! (pointer becomes disassociated)
      end if

    end do

  end subroutine delete_eltree_in_elements

end program eltree10
