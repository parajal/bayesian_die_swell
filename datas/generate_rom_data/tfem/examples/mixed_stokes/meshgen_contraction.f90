module meshgen_contraction_m

  use meshgen_m

  implicit none

  type ic_t
    integer :: elshape=6, p=0, l=1
    integer :: n1=5, n2=5, n3=5, n4=3, n5=10, m1=5, m2=7
    real(dp) :: refine=1.0_dp
    real(dp) :: factor1=0.5_dp, factor2=0.2_dp
    real(dp) :: L1=10, L2=1, L3=1, L4=2, L5=10
    real(dp) :: H1=1.5_dp, H2=2.5_dp, H3=1.0_dp
  end type ic_t

contains

! 2D contraction

  subroutine contraction2D ( ic, mesh )

    type(ic_t), intent(in) :: ic
    type(mesh_t), intent(out) :: mesh

    type(meshgen_options_t) :: mesh_options
    type(mesh_t) :: mesh1, mesh2, mesh3

    real(dp) :: x2d(4,2)

    integer :: elshape, n1, n2, n3, n4, n5, m1, m2
    real(dp) :: factor1, factor2, L1, L2, L3, L4, L5, H1, H2, H3
!
!  contraction build from 7 submeshes:
!
!
!         ------------------
!         |        |       |
!  m2  H2 |   2    |   4   |
!         |        |       |
!         |--------|\____  |
!         |        |     \_|---------------------------
!  m1  H1 |   1    |   3   |    |    |                |  H3
!         |        |       | 5  | 6  |      7         |
!       0 ---------------------------------------------
!        -L1      -L2      0    L3   L4               L5
!             n1      n2     n3   n4        n5
!
! NOTE: one single group is created. Use nogroupmerge=.false. in heading of
! mesh_merge to obtain 7 separate element groups.


    elshape = ic%elshape
    n1 = nint( ic%n1 * ic%refine )
    n2 = nint( ic%n2 * ic%refine )
    n3 = nint( ic%n3 * ic%refine )
    n4 = nint( ic%n4 * ic%refine )
    n5 = nint( ic%n5 * ic%refine )
    m1 = nint( ic%m1 * ic%refine )
    m2 = nint( ic%m2 * ic%refine )
    factor1 = ic%factor1
    factor2 = ic%factor2
    L1 = ic%L1
    L2 = ic%L2
    L3 = ic%L3
    L4 = ic%L4
    L5 = ic%L5
    H1 = ic%H1
    H2 = ic%H2
    H3 = ic%H3

! 1

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=m1, ox=-L1, &
    lx=L1-L2, ly=H1, p=ic%p, l=ic%l )

  call quadrilateral2d ( mesh1, mesh_options )

! 2

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=m2, ox=-L1, &
    oy=H1, lx=L1-L2, ly=H2, p=ic%p, l=ic%l )

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, curve1=3, curve2=-1 )
  call delete ( mesh1, mesh2 )

! 3

  x2d(1,:) = [-L2,0._dp]
  x2d(2,:) = [0,0]
  x2d(3,:) = [0._dp,H3]
  x2d(4,:) = [-L2,H1]

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=m1, &
    x2d=x2d, regionshape=2, ratio=[6,6,5,0], p=ic%p, l=ic%l )

  mesh_options%factor = factor1

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh3, mesh2, mesh1, curve1=2, curve2=-4 )
  call delete ( mesh3, mesh2 )

! 4

  x2d(1,:) = [-L2,H1]
  x2d(2,:) = [0._dp,H3]
  x2d(3,:) = [0._dp,H1+H2]
  x2d(4,:) = [-L2,H1+H2]

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=m2, &
    x2d=x2d, regionshape=2, ratio=[6,5,5,0], p=ic%p, l=ic%l )

  mesh_options%factor = factor1
  mesh_options%factor(2) = factor2

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, curves1=[10,5], curves2=[-1,-4] )
  call delete ( mesh1, mesh2 )

! 5

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n3, ny=m1, &
    lx=L3, ly=H3, ratio=[5,6,6,5], p=ic%p, l=ic%l )

  mesh_options%factor = factor1

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh3, mesh1, mesh2, curve1=9, curve2=-4 )
  call delete ( mesh3, mesh1 )

! 6

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n4, ny=m1, ox=L3, &
    lx=L4-L3, ly=H3, ratio=[0,0,0,5], p=ic%p, l=ic%l )

  mesh_options%factor = factor1

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh2, mesh1, mesh3, curve1=14, curve2=-4 )
  call delete ( mesh2, mesh1 )

! 7

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n5, ny=m1, ox=L4, &
    lx=L5-L4, ly=H3, p=ic%p, l=ic%l )

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh3, mesh1, mesh, curve1=17, curve2=-4 )
  call delete ( mesh3, mesh1 )

  end subroutine contraction2D

end module meshgen_contraction_m

