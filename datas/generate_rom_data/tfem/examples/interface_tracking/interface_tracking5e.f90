! Tracking of an initially circular interface in a flow given by a function.
! Interface movement is Lagrange based.
! Tangential grid deformation.
! Up to fourth-order explicit time integration with prediction.
! tfem mesh made using mesh_skeleton.
! Error output

program interface_tracking5e

  use tfem_m
  use math_defs_m
  use hsl_ma57_m
  use interface_tracking_elements_m
  use poisson_gd_ma57_m
  use functions_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  integer, parameter ::  &
    uintpl = 6,          & ! interpolation of position on the interface
    gauss = 4,           & ! number of integration points of line elements
    ndim = 2,            & ! space dimension
    nelem = 100,         & ! number of elements
    vfuncnr=1,           & ! function number for the velocity field u
    vtkevery = 10,       & ! vtk file every vtkevery steps. -1: means none
    timeint=3,           & ! time integration scheme
    numtimesteps=100       ! number of time steps

  real(dp), parameter :: &
    factor = 100._dp, &  ! factor for tangential grid deformation
    deltat = 5.e-3_dp   ! time step

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(problem_t), target :: problem_pois
  type(sysmatrix_t) :: sysmatrix, sysmatrix_pois
  type(sysvector_t) :: sol, rhsd(ndim)
  type(sysvector_t), target :: sol_pois
  type(vector_t) :: xdot, error
  type(coefficients_t) :: coefficients
  type(plot_options_t) :: plot_options
  type(lu_ma57_t) :: lu
  type(oldvectors_t) :: oldvectors

! Interface mesh coordinates at previous times
  real(dp), allocatable, dimension(:,:) :: mesh_coor_n, mesh_coor_nm1, &
    mesh_coor_nm2, mesh_coor_nm3


  character(len=30) :: filename

  integer :: i, nnodes, step, post, elem, node


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(3) = gauss
  coefficients%i(5) = 3  ! standard Gauss-Legendre (numerical tables)
  coefficients%i(6) = 1  ! velocity given by a function
  coefficients%i(7) = vfuncnr
  coefficients%i(8) = 0  ! arbitrary field given by a constant
  coefficients%i(10) = 1  ! tangential grid velocity using Poisson problem

  coefficients%r = 0
  coefficients%r(7) = factor

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

  call plot_points_curves ( plot_options, mesh, filename='curves.fig' )

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )
  call write_mesh_vtk ( mesh, 'mesh.vtk' )

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

! create oldvectors

  call create ( oldvectors, nsysvec=1, nprob=1 )


! open monitor file

  open ( unit=10, file='error5.out', recl=300 )


! time stepping

  post = 0

  do step = 1, numtimesteps

!   move the interface position

    call move_interface

!   fill error vector

    time = step * deltat

    call fill_vector ( mesh, problem, error, &
      node1=1, node2=mesh%nnodes, func=func_error_2D, funcnr=vfuncnr )

    write(unit=10,fmt=*) time, maxval(error%u), &
                         sqrt(sum(error%u**2)/mesh%nnodes)

!   output to vtk

    if ( vtkevery > 0 ) then
      if ( mod(step,vtkevery) == 0 ) then

        post = post + 1
        write(filename,'(a,i4.4,a)') 'output', post, '.vtk'

!       xdot
        call write_vector_vtk ( mesh, problem, filename=filename, &
          dataname='xdot', vector=xdot )

!       error
        call write_scalar_vtk ( mesh, problem, filename=filename, &
          dataname='error', vector=error, append=.true. )

!       output coordinates
        write(filename,'(a,i4.4,a)') 'coor', post, '.out'
        open ( unit=11, file=filename, recl=300 )
        write(unit=11,fmt='(i0,2e16.8)') (i, mesh%coor(i,:), i=1,mesh%nnodes)
        close ( unit=11 )

      end if
    end if

  end do

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol )
  call delete ( rhsd )
  call delete ( xdot, error )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( problem_pois )
  call delete ( sol_pois )

  deallocate ( mesh_coor_n )
  deallocate ( mesh_coor_nm1 )
  deallocate ( mesh_coor_nm2 )
  deallocate ( mesh_coor_nm3 )

contains

! build and solve for xdot on interface

  subroutine build_solve_interface

    integer :: i

!   Poisson problem on interface mesh for making nodes "diffuse" and
!   redistribute uniformly

    call poisson_gd_ma57 ( mesh, problem_pois, coefficients, sysmatrix_pois, &
      sol_pois )

    oldvectors%s(1)%p => sol_pois
    oldvectors%p(1)%p => problem_pois

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, msysvector=rhsd, &
      elemsub=interface_tracking_elem_xdot_lagrange, &
      coefficients=coefficients, oldvectors=oldvectors )

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

end program interface_tracking5e
