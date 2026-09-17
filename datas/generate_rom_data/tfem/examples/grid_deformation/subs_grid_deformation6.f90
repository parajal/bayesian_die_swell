module subs_grid_deformation6_m

  use tfem_elem_m

  implicit none
  save

! radius and center of the circle
  real(dp) :: rpl = 1._dp, cpl(2) = [0._dp,0._dp]

! the levelset function in all the nodes
  real(dp), dimension(:), allocatable :: d

! the minimum and maximum value of the monitor function f
  real(dp) :: lmin_f = 0.1, lmax_f = 1.0


contains

  function levelset ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1)) :: levelset

    real(dp) :: c(size(x,1),size(x,2))

    c(:,1) = cpl(1)
    c(:,2) = cpl(2)

!   circle with center at cpl and radius rpl
    levelset = sqrt(sum((x-c)**2,dim=2)) - rpl

  end function levelset


! monitor function

  subroutine set_monitor_function ( mesh, f_monitor )

    type(mesh_t), intent(in) :: mesh
    real(dp), dimension(:), intent(out) :: f_monitor

    real(dp) :: d(mesh%nnodes)

    d = levelset ( mesh%coor )
    d = (1-mesh%coor(:,2))* abs( d )

    f_monitor = min( lmax_f, max(d,lmin_f) )

  end subroutine set_monitor_function


! correct for points moved (numerically) out of the domain at the Neumann
! boundaries after a time step.

  subroutine project_boundary ( mesh, object )

    type(mesh_t), intent(inout) :: mesh
    integer, intent(in) :: object

    integer :: curve

    curve = 1
    mesh%objects(object)%coor(mesh%curves(curve)%nodes,2) = 0

    curve = 2
    mesh%objects(object)%coor(mesh%curves(curve)%nodes,1) = 1

    curve = 3
    mesh%objects(object)%coor(mesh%curves(curve)%nodes,2) = 1

    curve = 4
    mesh%objects(object)%coor(mesh%curves(curve)%nodes,1) = 0

  end subroutine project_boundary

end module subs_grid_deformation6_m
