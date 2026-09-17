! Compute area, centroid and second moment of area tensor for a cross section

program second_moment1

  use tfem_m
  use interface_tracking_elements_m
  use eig2D3d_m
  use io_utils_m
  use limits_m

  implicit none

! constants

  integer, parameter ::  &
    ndim = 2               ! dimension of space

! variables

  integer ::  &
    uintpl = 2, &          ! P1 interface position interpolation
!    uintpl = 6, &         ! P2 interface position interpolation
    gaussb = 4             ! number of integration points of line elements

! use interface tracking problem

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(coefficients_t) :: coefficients

  character (len=20) :: meshfile
  real(dp) :: resultsum(1+ndim*(ndim+1)/2), A, Q(ndim), x_c(ndim), &
    Jxx, Jxy, Jyy, Jxx_c, Jxy_c, Jyy_c, &
    Ixx, Ixy, Iyy, Ixx_c, Ixy_c, Iyy_c, Iten_c(ndim*(ndim+1)/2), &
    lambda(ndim), eigv(ndim,ndim)

! namelist for input of variables; read from standard input

  namelist /comppar/ uintpl, gaussb, meshfile

  read ( unit=*, nml=comppar )

  WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = .false.

! Define coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(3) = gaussb
  coefficients%i(5) = 3  ! standard Gauss-Legendre (numerical tables)

  coefficients%r = 0

! read mesh of cross section

  print *
  write(*,'(2a)') 'Meshfile = ', meshfile
  print *

  call read_mesh_gmsh ( mesh, meshfile, ndim=ndim )

  call fill_mesh_parts ( mesh )

  call printinfo ( mesh, printlevel=1 )


! problem definition of the boundary

  call create_input_probdef ( mesh, input_probdef )

  input_probdef%elementdof(1)%a = 1

  call problem_definition ( input_probdef, mesh, problem )


! compute area and centroid

  call integrate ( mesh, problem, resultsum(:ndim+1), &
    elemsub=center_of_volume, coefficients=coefficients )

  A = resultsum(1)
  Q = resultsum(2:ndim+1)
  x_c = Q / A

! compute J =\int x x dA tensor
!         =       - -

  call integrate ( mesh, problem, resultsum, elemsub=second_moment_of_area,&
    coefficients=coefficients )

  Jxx = resultsum(2)
  Jxy = resultsum(3)
  Jyy = resultsum(4)

! Use Steiner shift towards centroid position:

  Jxx_c = Jxx - A * x_c(1)**2
  Jxy_c = Jxy - A * x_c(1)*x_c(2)
  Jyy_c = Jyy - A * x_c(2)**2

! Compute
!    Ixx=\int y^2 dA, Ixy=\int xy dA, Iyy=\int x^2 dA
! in the global system and the local (centroid) system. Note, that
! Ixx, -Ixy and Iyy are the tensor components of the "moment of inertia" tensor
!
!  I = tr(J) 1 - J
!  =      =  =   =

  Ixx = Jyy
  Ixy = Jxy
  Iyy = Jxx
  Ixx_c = Jyy_c
  Ixy_c = Jxy_c
  Iyy_c = Jxx_c

! Output data

  print *
  write(*,'(a)') 'Note:'
  write(*,'(a)') '  Ixx=\int y^2 dA '
  write(*,'(a)') '  Ixy=\int xy dA '
  write(*,'(a)') '  Iyy=\int x^2 dA '
  print *
  write(*,'(a)') 'Values in global xy-system:'
  print *
  print *, 'area  = ', A
  print *, 'x_c   = ', x_c
  print *, 'Ixx, Ixy, Iyy = '
  print *,  Ixx, Ixy, Iyy
  print *
  write(*,'(a)') 'Values in centroid xy-system:'
  print *
  print *, 'Ixx_c, Ixy_c, Iyy_c = '
  print *,  Ixx_c, Ixy_c, Iyy_c

! principal values

  Iten_c = [ Ixx_c, -Ixy_c, Iyy_c ]
  call eig2x2 ( Iten_c, lambda, eigv )

  print *
  write(*,'(a)') 'Principal values and directions of'
  write(*,'(a)') '          I_tensor_c = [  Ixx_c -Ixy_c ] '
  write(*,'(a)') '                       [ -Ixy_c  Iyy_c ]:'
  print *
  print *, 'lambda = ', lambda
  print *, 'eigv = '
  print *, eigv(1,:)
  print *, eigv(2,:)
  print *

! delete all data including all allocated memory

  call delete ( coefficients )
  call delete ( mesh )
  call delete ( input_probdef )
  call delete ( problem )

end program second_moment1
