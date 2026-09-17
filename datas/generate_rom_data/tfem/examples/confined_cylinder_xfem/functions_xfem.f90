! functions for the cylinder problem

module functions_xfem_m

  use tfem_elem_m
  use math_defs_m

  implicit none
  save

! radius and center of the cylinder
  real(dp) :: rpl = 1._dp, cpl(2) = [0._dp,0._dp]

! mesh used for mapping of reference element
  type(mesh_t), pointer :: lmesh => null()

  integer :: lelem, lelgrp

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


! function for mapping reference coordinates to the real coordinates
! in building the eltree within an element.

  function mapcoor ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1),size(x,2)) :: mapcoor

    real(dp) :: phi(size(x,1),9), xnod(9,2)

    call shape_quad_Q2 ( x, phi )

    call get_coordinates ( lmesh, lelgrp, lelem, xnod )

    mapcoor = matmul ( phi, xnod ) ! isoparametric

  end function mapcoor


! objectscoor defines the coordinates of the objects

  subroutine objectscoor ( objectnr, coor )
    integer, intent(in) :: objectnr
    real(dp), dimension(:,:), intent(inout) :: coor

    integer :: np, i
    real(dp) :: p(size(coor,1))

    np = size(coor,1)
    p = [(2*pi/np*(i-1),i=1,np)]

    coor(:,1) = rpl * sin(p) + cpl(1)
    coor(:,2) = rpl * cos(p) + cpl(2)

  end subroutine objectscoor

end module functions_xfem_m
