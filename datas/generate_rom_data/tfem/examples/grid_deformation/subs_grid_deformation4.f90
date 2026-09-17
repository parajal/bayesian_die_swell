module subs_grid_deformation4_m

  use tfem_elem_m

  implicit none
  save

! radius and center of the sphere
  real(dp) :: rpl = 1._dp, cpl(3) = [0._dp,0._dp,0._dp]

! the levelset function in all the nodes
  real(dp), dimension(:), allocatable :: d

! the domain
  real(dp) :: lox, loy, loz, llx, lly, llz

! the minimum and maximum value of the monitor function f
  real(dp) :: lmin_f = 0.1, lmax_f = 1.0

contains

  function levelset ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1)) :: levelset

    real(dp) :: c(size(x,1),size(x,2))

    c(:,1) = cpl(1)
    c(:,2) = cpl(2)
    c(:,3) = cpl(3)

!   sphere with center at cpl and radius rpl
    levelset = sqrt(sum((x-c)**2,dim=2)) - rpl

    c(:,1) = cpl(1) + llx

    levelset = min ( levelset, sqrt(sum((x-c)**2,dim=2)) - rpl )

    c(:,1) = cpl(1) - llx

    levelset = min ( levelset, sqrt(sum((x-c)**2,dim=2)) - rpl )

  end function levelset


! monitor function

  subroutine set_monitor_function ( mesh, f_monitor )

    type(mesh_t), intent(in) :: mesh
    real(dp), dimension(:), intent(out) :: f_monitor

    real(dp) :: d(mesh%nnodes)

    d = levelset ( mesh%coor )
    d = abs( d )

    f_monitor = min( lmax_f, max(d,lmin_f) )

  end subroutine set_monitor_function


! Define constraints for the Poisson problem.

  subroutine constraint_definition ( mesh, input_probdef )

    type(mesh_t), intent(inout) :: mesh
    type(input_probdef_t), intent(inout) :: input_probdef

    call define_constraint ( mesh, input_probdef, &
      surface1=5, surface2=3, discretization='collocation' )

  end subroutine constraint_definition


! move object points before intersection with the mesh
! NOTE: this only affects the coordinates before intersecting. After the time
! step the coordinates are updated correctly. This is for handling periodical
! domains.

  subroutine move_object_points ( mesh, object )

    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: object

    where ( mesh%objects(object)%coor(:,1) < lox )
      mesh%objects(object)%coor(:,1) = mesh%objects(object)%coor(:,1) + llx
    else where ( mesh%objects(object)%coor(:,1) > lox + llx )
      mesh%objects(object)%coor(:,1) = mesh%objects(object)%coor(:,1) - llx
    end where

  end subroutine move_object_points


! correct for points moved (numerically) out of the domain at the Neumann
! boundaries after a time step.

  subroutine project_boundary ( mesh, object )

    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: object

    integer :: surface, surface2

    surface = 1
    mesh%objects(object)%coor(mesh%surfaces(surface)%nodes,3) = loz

    surface = 2
    mesh%objects(object)%coor(mesh%surfaces(surface)%nodes,2) = loy

    surface = 4
    mesh%objects(object)%coor(mesh%surfaces(surface)%nodes,2) = loy + lly

    surface = 6
    mesh%objects(object)%coor(mesh%surfaces(surface)%nodes,3) = loz + llz

!   exact periodic node positions
    surface = 3
    surface2 = 5
    mesh%objects(object)%coor(mesh%surfaces(surface)%nodes,1) = &
           mesh%objects(object)%coor(mesh%surfaces(surface2)%nodes,1) + llx
    mesh%objects(object)%coor(mesh%surfaces(surface)%nodes,2:3) = &
           mesh%objects(object)%coor(mesh%surfaces(surface2)%nodes,2:3)


  end subroutine project_boundary

end module subs_grid_deformation4_m
