! Tracking of an initially circular interface in a flow given by a function.
! Interface moves in normal direction.
! tfem mesh made using mesh_skeleton.
! Up to fourth-order Gear with prediction.
! Error output.

module functions_m

  use kind_defs_m

  implicit none

  real(dp), save :: epsilondot = 1._dp, gammadot = 1._dp, omega = 1._dp, &
    time = 1._dp, R0 = 1._dp

contains

  function vfunc_2D ( n, nr, x )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc_2D

    select case(nr)
      case(1) ! elongation
        vfunc_2D(1) = x(1)
        vfunc_2D(2) = -x(2)
        vfunc_2D = epsilondot * vfunc_2D
      case(2) ! shear
        vfunc_2D(1) = gammadot * x(2)
        vfunc_2D(2) = 0
      case(3) ! rotation
        vfunc_2D(1) = x(2)
        vfunc_2D(2) = -x(1)
        vfunc_2D = omega * vfunc_2D
      case default
        write(*,'(/a,i0/)') 'Error vfunc_2D: wrong function number: ', nr
        stop
    end select

  end function vfunc_2D

  function func_error_2D ( nr, x )

    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func_error_2D
    real(dp) :: a, b, lambda, x0, y0

    select case(nr)
      case(1) ! elongation
        a = exp(epsilondot*time)
        b = exp(-epsilondot*time)
        lambda = R0 / sqrt( (x(1)/a)**2 + (x(2)/b)**2 )
        func_error_2D = abs(lambda-1)* sqrt( x(1)**2 + x(2) ** 2 )
      case(2) ! shear
        x0 = x(1) - gammadot * x(2) * time
        y0 = x(2)
        lambda = R0 / sqrt( x0**2 + y0**2 )
        func_error_2D = abs(lambda-1)* sqrt( x(1)**2 + x(2) ** 2 )
      case(3) ! rotation
        func_error_2D = sqrt( x(1)**2 + x(2) ** 2 ) - R0
      case default
        write(*,'(/a,i0/)') 'Error func_error_2D: wrong function number: ', nr
        stop
    end select

  end function func_error_2D

  function vfunc_3D ( n, nr, x )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc_3D

    select case(nr)
      case(1) ! elongation
        vfunc_3D(1) = x(1)
        vfunc_3D(2) = -x(2)/2
        vfunc_3D(3) = -x(3)/2
        vfunc_3D = epsilondot * vfunc_3D
      case(2) ! shear
        vfunc_3D(1) = gammadot * x(2)
        vfunc_3D(2) = 0
        vfunc_3D(3) = 0
      case(3) ! rotation
        vfunc_3D(1) = x(2)
        vfunc_3D(2) = -x(1)
        vfunc_3D(3) = 0
        vfunc_3D = omega * vfunc_3D
      case(4) ! plane elongation
        vfunc_3D(1) = x(1)
        vfunc_3D(2) = -x(2)
        vfunc_3D(3) = 0
        vfunc_3D = epsilondot * vfunc_3D
      case default
        write(*,'(/a,i0/)') 'Error vfunc_3D: wrong function number: ', nr
        stop
    end select

  end function vfunc_3D

  function func_error_3D ( nr, x )

    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func_error_3D
    real(dp) :: a, b, lambda, x0, y0, z0

    select case(nr)
      case(1) ! elongation
        a = exp(epsilondot*time)
        b = exp(-epsilondot*time/2)
        lambda = R0 / sqrt( (x(1)/a)**2 + (x(2)/b)**2 + (x(3)/b)**2 )
        func_error_3D = abs(lambda-1)* sqrt( x(1)**2 + x(2) ** 2 + x(3) ** 2 )
      case(2) ! shear
        x0 = x(1) - gammadot * x(2) * time
        y0 = x(2)
        z0 = x(3)
        lambda = R0 / sqrt( x0**2 + y0**2 + z0**2 )
        func_error_3D = abs(lambda-1)* sqrt( x(1)**2 + x(2) ** 2 + x(3) ** 2 )
      case(3) ! rotation
        func_error_3D = sqrt( x(1)**2 + x(2) ** 2 + x(3) ** 2 ) - R0
      case(4) ! plane elongation
        a = R0 * exp(epsilondot*time)
        b = R0 * exp(-epsilondot*time)
        lambda = 1 / sqrt( (x(1)/a)**2 + (x(2)/b)**2 + (x(3)/R0)**2 )
        func_error_3D = abs(lambda-1)* sqrt( x(1)**2 + x(2) ** 2 + x(3) ** 2 )
      case default
        write(*,'(/a,i0/)') 'Error func_error_3D: wrong function number: ', nr
        stop
    end select

  end function func_error_3D

end module functions_m
program interface_tracking1

  use tfem_m
  use math_defs_m
  use hsl_ma57_m
  use interface_tracking_elements_m
  use functions_m
  use io_utils_m
  use limits_m

  implicit none

! constants

  integer, parameter ::  &
    uintpl = 6,          & ! interpolation of position on the interface
    gauss = 4,           & ! number of integration points of line elements
    ndim = 2,            & ! space dimension
    nelem = 100,         & ! number of elements
    vfuncnr=1,           & ! function number for the velocity field u
    !vtkevery = 10,       & ! vtk file every vtkevery steps. 0: means none
    timeint=3,           & ! time integration scheme
    numtimesteps=10       ! number of time steps

  real(dp), parameter :: &
    Um(ndim) = [ 0._dp, 0._dp ], &  ! arbitrary velocity u_m (constant)
    deltat = 5.e-3_dp   ! time step

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd(ndim)
  type(vector_t) :: xdot, error
  type(coefficients_t) :: coefficients
  type(lu_ma57_t) :: lu

! Interface mesh coordinates at previous times
  real(dp), allocatable, dimension(:,:) :: mesh_coor_n, mesh_coor_nm1, &
    mesh_coor_nm2, mesh_coor_nm3


  integer :: nnodes, step, elem, node

  WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = .false.

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(3) = gauss
  coefficients%i(5) = 3  ! standard Gauss-Legendre (numerical tables)
  coefficients%i(6) = 1  ! velocity given by a function
  coefficients%i(7) = vfuncnr
  coefficients%i(8) = 0  ! arbitrary field given by a constant


  coefficients%r = 0
  coefficients%r(4:5) = Um

  coefficients%vfunc1(1)%p => vfunc_2D

! create mesh of a circular curve discretized with 3-node line elements

  call mesh_skeleton ( mesh, nnodes=2*nelem, nelem=nelem, elshape=2, ndim=ndim )

  do elem = 1, mesh%nelem
    mesh%topology(1)%a(:,elem) = [ 2*elem-1, 2*elem, 2*elem+1 ]
  end do

  mesh%topology(1)%a(3,nelem) = 1

  do node = 1, mesh%nnodes
    mesh%coor(node,1) = R0 * cos((node-1)*2*pi/mesh%nnodes)
    mesh%coor(node,2) = R0 * sin((node-1)*2*pi/mesh%nnodes)
  end do

  call fill_mesh_parts ( mesh )

  call printinfo ( mesh, printlevel=4 )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a(:,1) = ndim  ! 2D vector
  input_probdef%vec_elementdof(1)%a(:,2) = 1  ! error

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_msysvector ( problem, rhsd )

! storage for velocity in vector

  call create ( problem, xdot, vec=1 )

! storage for error

  call create ( problem, error, vec=2 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, &
    symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )


! arrays for storing coordinates of interface at old times,

  nnodes = mesh%nnodes
  allocate ( mesh_coor_n(nnodes,ndim) )
  allocate ( mesh_coor_nm1(nnodes,ndim) )
  allocate ( mesh_coor_nm2(nnodes,ndim) )
  allocate ( mesh_coor_nm3(nnodes,ndim) )

! time stepping

  do step = 1, numtimesteps

!   move the interface position

    call move_interface

!   fill error vector

    time = step * deltat

    call fill_vector ( mesh, problem, error, &
      node1=1, node2=mesh%nnodes, func=func_error_2D, funcnr=vfuncnr )

    write(*,*) time, maxval(error%u)

  end do

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol )
  call delete ( rhsd )
  call delete ( xdot )
  call delete ( sysmatrix )
  call delete ( coefficients )

  deallocate ( mesh_coor_n )
  deallocate ( mesh_coor_nm1 )
  deallocate ( mesh_coor_nm2 )
  deallocate ( mesh_coor_nm3 )

contains

! build and solve for xdot on interface

  subroutine build_solve_interface

    integer :: i

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, msysvector=rhsd, &
      elemsub=interface_tracking_elem, coefficients=coefficients )

    do i = 1, ndim

      call solve_system_ma57 ( sysmatrix, rhsd(i), sol, lu=lu )

      call transfer_data ( mesh, problem1=problem, &
        sysvector1=sol, vector2=xdot, degfd1=[1], degfd2=[i] )

    end do

    call delete(lu)

  end subroutine build_solve_interface


! move the interface position (time discretized step)

  subroutine move_interface

    if ( timeint == 1 .or. ( step == 1 .and. timeint > 1 ) ) then

!     Euler

      mesh_coor_n = mesh%coor

      call build_solve_interface

!     Explicit Euler updating of interface mesh coordinates
      mesh%coor = mesh_coor_n + &
                     deltat * transpose( reshape((xdot%u),[ndim,nnodes]) )

    else if ( timeint == 2 .or. ( step == 2 .and. timeint > 2 ) ) then

!     second order Gear with prediction

      mesh_coor_nm1 = mesh_coor_n
      mesh_coor_n = mesh%coor

!     2nd order prediction of interface coordinates
      mesh%coor = 2*mesh_coor_n - mesh_coor_nm1

      call build_solve_interface

!     2nd order Gear updating of interface mesh coordinates
      mesh%coor = 2 * ( 2*mesh_coor_n - mesh_coor_nm1/2 &
                + deltat * transpose( reshape((xdot%u),[ndim,nnodes]) ) ) / 3

    else if ( timeint == 3 .or. ( step == 3 .and. timeint > 3 ) ) then

!     third order Gear with prediction

      mesh_coor_nm2 = mesh_coor_nm1
      mesh_coor_nm1 = mesh_coor_n
      mesh_coor_n = mesh%coor

!     3rd order prediction of interface coordinates
      mesh%coor = 3*mesh_coor_n - 3*mesh_coor_nm1 + mesh_coor_nm2

      call build_solve_interface

!     3rd order Gear updating of interface mesh coordinates
      mesh%coor = 6 * ( 3*mesh_coor_n - 3*mesh_coor_nm1/2 + mesh_coor_nm2/3 &
               + deltat * transpose( reshape((xdot%u),[ndim,nnodes]) ) ) / 11

    else if ( timeint == 4 .or. ( step == 4 .and. timeint > 4 ) ) then

!     fourth order Gear with prediction

      mesh_coor_nm3 = mesh_coor_nm2
      mesh_coor_nm2 = mesh_coor_nm1
      mesh_coor_nm1 = mesh_coor_n
      mesh_coor_n = mesh%coor

!     4th order prediction of interface coordinates
      mesh%coor = 4*mesh_coor_n - 6*mesh_coor_nm1 + 4*mesh_coor_nm2 - &
                  mesh_coor_nm3

      call build_solve_interface

!     4th order Gear updating of interface mesh coordinates
      mesh%coor = 12 * ( 4*mesh_coor_n - 3*mesh_coor_nm1 + 4*mesh_coor_nm2/3 &
                                       - mesh_coor_nm3/4 &
              + deltat * transpose( reshape((xdot%u),[ndim,nnodes]) ) ) / 25

    end if

  end subroutine move_interface

end program interface_tracking1
