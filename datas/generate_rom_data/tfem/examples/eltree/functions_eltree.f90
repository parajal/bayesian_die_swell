! functions to be used with the eltree examples using triangles and tetrahedra.

module functions_eltree_m

  use math_defs_m
  use tfem_elem_m

  implicit none
  save

! mesh used for mapping of reference element
  type(mesh_t), pointer :: lmesh => null()

  integer :: lelem, lelgrp

  real(dp) :: rpl = 1.0_dp

contains

! levelset function for a circle / sphere

  function levelset ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1)) :: levelset

    levelset = sqrt(sum(x**2,dim=2)) - rpl

  end function levelset

! functions for mapping reference coordinates to the real coordinates
! in building the eltree within an element.

! line element
  function mapcoor_line ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1),size(x,2)) :: mapcoor_line

    real(dp) :: phi(size(x,1),2), xnod(2,1)

!   get the values of the shape functions at the reference coordinates
    call shape_line_P1 ( x, phi )

!   get the global coordinates of an element in the vertices
    call get_coordinates ( lmesh, lelgrp, lelem, xnod )

!   find the global coordinates associated with the reference coordinates
!   by multiplying the shapefunctions with the global coordinates of the
!   vertices (isoparametric)
    mapcoor_line = matmul ( phi, xnod )

  end function mapcoor_line

! triangular element
  function mapcoor_triangle ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1),size(x,2)) :: mapcoor_triangle

    real(dp) :: phi(size(x,1),3), xnod(3,2)

    call shape_triangle_P1 ( x, phi )

    call get_coordinates ( lmesh, lelgrp, lelem, xnod )

    mapcoor_triangle = matmul ( phi, xnod )

  end function mapcoor_triangle

! tetrahedral element
  function mapcoor_tet ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1),size(x,2)) :: mapcoor_tet

    real(dp) :: phi(size(x,1),4), xnod(4,3)

    call shape_tetra_P1 ( x, phi )

    call get_coordinates ( lmesh, lelgrp, lelem, xnod )

    mapcoor_tet = matmul ( phi, xnod )

  end function mapcoor_tet

! quadrilateral element
  function mapcoor_quad ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1),size(x,2)) :: mapcoor_quad

    real(dp) :: phi(size(x,1),4), xnod(4,2)

    call shape_quad_Q1 ( x, phi )

    call get_coordinates ( lmesh, lelgrp, lelem, xnod )

    mapcoor_quad = matmul ( phi, xnod )

  end function mapcoor_quad

! hexahedral element
  function mapcoor_hexa ( x )

    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(size(x,1),size(x,2)) :: mapcoor_hexa

    real(dp) :: phi(size(x,1),8), xnod(8,3)

    call shape_hexa_Q1 ( x, phi )

    call get_coordinates ( lmesh, lelgrp, lelem, xnod )

    mapcoor_hexa = matmul ( phi, xnod )

  end function mapcoor_hexa

end module functions_eltree_m
