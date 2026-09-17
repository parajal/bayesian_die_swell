! Global variables for the elements needed by the code

module interface_tracking_globals_m

  use tfem_elem_m
  use math_defs_m

  implicit none

  save

! global element shape

  character (len=13) :: globalshape = ''

! space dimension
  integer :: ndim = 0

! coordinate system: 0=Cartesian (x,y), 2D
!                    1=Cylindrical (z,r), axisymmetrical
!                    2=Cartesian (x,y,z), 3D
  integer :: coorsys = 0

! interpolation number (shape function) for the position
  integer :: intpol = 0

! number of position degrees of freedom with respect to one direction only
  integer :: ndf = 0

! interpolation number (shape function) for the velocity
  integer :: intpolv = 0

! number of velocity degrees of freedom
  integer :: ndfv = 0

! shape function used for the "isoparametric" deformation of the element
! 0=position, 1=velocity
  integer :: isoshape = 0

! integration rule for the interior of an element
  integer :: intrule = 0

! integration type
  integer :: inttype = 0

! number of integration points for the interior of an element
  integer :: ninti = 0

! number of subdomains for the integration (generalized integration rules)
  integer :: nsubint = 0

! number of nodal points of the interior element
  integer :: nodalp = 0

! various allocatable arrays
  real(dp), allocatable, dimension(:) :: wg, fg, detF, work, work1, work2, &
    un, surfl, c, curvel, u, xn, x2, Uc, facs, xar
  real(dp), allocatable, dimension(:,:) :: phi, phiv, x, xig, xg, normal, &
    uvecg, uvecmg, uvectg, dxdxi, g1_up, work3, ugradphi, xng, xnm1g, &
    uvecg_conv, uvecg_rhs, tmp, xrnod

  real(dp), allocatable, dimension(:,:,:) :: dphi, dphiv, F, Finv, dphidx, &
    dxdxis, gi_up, gij_up, gij_down

! gauss types

  type(gauss_t) :: gauss

! shapefunc types

  type(shapefunc_t) :: shapefunc, shapefuncv

end module interface_tracking_globals_m
