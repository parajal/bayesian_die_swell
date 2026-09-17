! Compute volume, centroid and second moment of volume tensor for a body

program second_moment2

  use tfem_m
  use interface_tracking_elements_m
  use eig2D3d_m
  use io_utils_m
  use limits_m

  implicit none

! constants

  integer, parameter ::  &
    ndim = 3               ! dimension of space

! variables

  integer ::  &
    uintpl = 2, &          ! P1 interface position interpolation
!    uintpl = 6, &         ! P2 interface position interpolation
    gaussb = 7             ! order of integration points of triangular elements

! use interface tracking problem

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(coefficients_t) :: coefficients

  character (len=20) :: meshfile
  real(dp) :: resultsum(1+ndim*(ndim+1)/2), V, Q(ndim), x_c(ndim), &
    Jxx, Jxy, Jxz, Jyy, Jyz, Jzz, trJ, &
    Jxx_c, Jxy_c, Jxz_c, Jyy_c, Jyz_c, Jzz_c, trJ_c, &
    Ixx, Ixy, Ixz, Iyy, Iyz, Izz, &
    Ixx_c, Ixy_c, Ixz_c, Iyy_c, Iyz_c, Izz_c, &
    Iten_c(ndim*(ndim+1)/2), lambda(ndim), eigv(ndim,ndim)

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

  V = resultsum(1)
  Q = resultsum(2:ndim+1)
  x_c = Q / V

! compute J =\int x x dA tensor
!         =       - -

  call integrate ( mesh, problem, resultsum, elemsub=second_moment_of_area,&
    coefficients=coefficients )

  Jxx = resultsum(2)
  Jxy = resultsum(3)
  Jxz = resultsum(4)
  Jyy = resultsum(5)
  Jyz = resultsum(6)
  Jzz = resultsum(7)

! Use Steiner shift towards centroid position:

  Jxx_c = Jxx - V * x_c(1)**2
  Jxy_c = Jxy - V * x_c(1)*x_c(2)
  Jxz_c = Jxz - V * x_c(1)*x_c(3)
  Jyy_c = Jyy - V * x_c(2)**2
  Jyz_c = Jyz - V * x_c(2)*x_c(3)
  Jzz_c = Jzz - V * x_c(3)**2

! Compute the tensor components of the "moment of inertia" tensor
!
!  I = tr(J) 1 - J
!  =      =  =   =
! in the global system and the local (centroid) system.

  trJ = Jxx + Jyy + Jzz
  Ixx = trJ - Jxx
  Ixy = -Jxy
  Ixz = -Jxz
  Iyy = trJ - Jyy
  Iyz = -Jyz
  Izz = trJ - Jzz

  trJ_c = Jxx_c + Jyy_c + Jzz_c
  Ixx_c = trJ_c - Jxx_c
  Ixy_c = -Jxy_c
  Ixz_c = -Jxz_c
  Iyy_c = trJ_c - Jyy_c
  Iyz_c = -Jyz_c
  Izz_c = trJ_c - Jzz_c

! Output data

  print *
  write(*,'(a)') 'Note:'
  write(*,'(a)') ' second moment:'
  write(*,'(a)') '  J =\int x x dA tensor'
  write(*,'(a)') '  =       - -'
  write(*,'(a)') ' the "moment of inertia" tensor:'
  write(*,'(a)') '  I = tr(J) 1 - J'
  write(*,'(a)') '  =      =  =   = '
  print *
  write(*,'(a)') 'Values in global xyz-system:'
  print *
  print *, 'volume  = ', V
  print *, 'x_c   = ', x_c
  print *, 'Ixx, Ixy, Ixz = '
  print *,  Ixx, Ixy, Ixz
  print *, 'Iyy, Iyz, Izz = '
  print *,  Iyy, Iyz, Izz
  print *
  write(*,'(a)') 'Values in centroid xyz-system:'
  print *
  print *, 'Ixx_c, Ixy_c, Ixz_c = '
  print *,  Ixx_c, Ixy_c, Ixz_c
  print *, 'Iyy_c, Iyz_c, Izz_c = '
  print *,  Iyy_c, Iyz_c, Izz_c

! principal values

  Iten_c = [ Ixx_c, Ixy_c, Ixz_c, Iyy_c, Iyz_c, Izz_c ]
  call eig3x3 ( Iten_c, lambda, eigv )

  print *
  write(*,'(a)') 'Principal values and directions of'
  write(*,'(a)') '          I_tensor_c = [ Ixx_c Ixy_c Ixz_c ] '
  write(*,'(a)') '                       [ Ixy_c Iyy_c Iyz_c]:'
  write(*,'(a)') '                       [ Ixz_c Iyz_c Izz_c]:'
  print *
  print *, 'lambda = ', lambda
  print *, 'eigv = '
  print *, eigv(1,:)
  print *, eigv(2,:)
  print *, eigv(3,:)
  print *

! delete all data including all allocated memory

  call delete ( coefficients )
  call delete ( mesh )
  call delete ( input_probdef )
  call delete ( problem )

end program second_moment2
