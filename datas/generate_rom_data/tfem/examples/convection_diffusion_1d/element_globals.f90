
! global variables for the element

module element_globals_m

  use tfem_elem_m

  implicit none

  save

! space dimension
  integer :: ndim = 1

! number of integration points
  integer :: ninti = 0

! number of nodal points
  integer :: nodalp = 0

! number of degrees of freedom
  integer :: ndf = 0

! various allocatable arrays
  real(dp), allocatable, dimension(:) :: wg, fg, detF, work
  real(dp), allocatable, dimension(:) :: xig
  real(dp), allocatable, dimension(:,:) :: phi, x, xg
  real(dp), allocatable, dimension(:,:,:) :: dphi, F, Finv, dphidx

end module element_globals_m
