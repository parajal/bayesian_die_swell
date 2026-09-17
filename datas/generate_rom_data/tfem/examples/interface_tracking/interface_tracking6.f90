! Tracking of an initially spherical interface in a flow given by a function.
! Interface movement is Lagrange based.
! Tangential grid deformation.
! Second-order time integration (performed on element level).
! Gmsh mesh with six-node triangles.
! Error output

program interface_tracking6

  use tfem_m
  use hsl_ma57_m
  use interface_tracking_elements_m
  use poisson_gd_ma57_m
  use functions_m
  use io_utils_m

  implicit none

! constants

  integer, parameter ::  &
!    uintpl = 2,          & ! interpolation of position on the interface
    uintpl = 6,          & ! interpolation of position on the interface
    gauss = 6,           & ! order of integration of triangles
    ndim = 3,            & ! space dimension
    vfuncnr=1,           & ! function number for the velocity field u
    vtkevery = 1,       & ! vtk file every vtkevery steps. -1: means none
    timeint=2,           & ! time integration scheme
    numtimesteps=20       ! number of time steps

  real(dp), parameter :: &
    factor = 3._dp, &  ! factor for tangential grid deformation
    deltat = 5.e-2_dp   ! time step

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(problem_t), target :: problem_pois
  type(sysmatrix_t) :: sysmatrix, sysmatrix_pois
  type(sysvector_t) :: sol, rhsd(ndim)
  type(sysvector_t), target :: sol_pois
  type(vector_t) :: xnp1, error
  type(vector_t), target :: xn, xnm1
  type(coefficients_t) :: coefficients
  type(lu_ma57_t) :: lu
  type(oldvectors_t) :: oldvectors
  type(solver_options_ma57_t) :: so


  character(len=30) :: filename

  integer :: step, post


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
  coefficients%i(11) = timeint

  coefficients%r = 0
  coefficients%r(7) = factor
  coefficients%r(9) = deltat

  coefficients%vfunc1(1)%p => vfunc_3D

  if ( uintpl == 2 ) then
!   read gmsh file of a sphere surface (3-node triangular elements)
    call read_mesh_gmsh ( mesh, 'mesh2l.msh' )
  else if ( uintpl == 6 ) then
!   read gmsh file of a sphere surface (6-node triangular elements)
    call read_mesh_gmsh ( mesh, 'mesh2q.msh' )
  end if

  call fill_mesh_parts ( mesh )

  call write_mesh_vtk ( mesh, 'mesh.vtk' )

  call printinfo ( mesh, printlevel=4 )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a(:,1) = ndim  ! 3D vector
  input_probdef%vec_elementdof(1)%a(:,2) = 1  ! error

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_msysvector ( problem, rhsd )

! storage for position at tn, tn+1 in vector

  call create ( problem, xn, vec=1 )
  call create ( problem, xnp1, vec=1 )

  if ( timeint == 2 ) then

!   storage for position at tn-1 in vector

    call create ( problem, xnm1, vec=1 )

  end if

! Initialize xnp1 (which will be the xn of the new time step)

  xnp1%u = reshape ( transpose(mesh%coor), [ndim*mesh%nnodes] )

! storage for error

  call create ( problem, error, vec=2 )

! oldvectors

  if ( timeint == 1 ) then
    call create_oldvectors ( oldvectors, nsysvec=1, nprob=1, nvec=3 )
    oldvectors%v(3)%p => xn
  else if ( timeint == 2 ) then
    call create_oldvectors ( oldvectors, nsysvec=1, nprob=1, nvec=4 )
    oldvectors%v(3)%p => xn
    oldvectors%v(4)%p => xnm1
  end if

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! open monitor file

  open ( unit=10, file='error4.out', recl=300 )


! time stepping

  post = 0

  do step = 1, numtimesteps

!   move the interface position

    call move_interface

!   fill error vector

    time = step * deltat

    call fill_vector ( mesh, problem, error, &
      node1=1, node2=mesh%nnodes, func=func_error_3D, funcnr=vfuncnr )

    write(unit=10,fmt=*) time, maxval(error%u), &
                         sqrt(sum(error%u**2)/mesh%nnodes)

!   output to vtk

    if ( vtkevery > 0 ) then
      if ( mod(step,vtkevery) == 0 ) then

        post = post + 1
        write(filename,'(a,i4.4,a)') 'output', post, '.vtk'

!       error
        call write_scalar_vtk ( mesh, problem, filename=filename, &
          dataname='error', vector=error )

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
  call delete ( xn, xnp1 )
  if ( timeint == 2 ) call delete ( xnm1 )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( problem_pois )
  call delete ( sol_pois )

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
      elemsub=interface_tracking_elem_lagrange, coefficients=coefficients, &
      oldvectors=oldvectors )

    do i = 1, ndim

      so%integer_storage=1.8_dp

      call solve_system_ma57 ( sysmatrix, rhsd(i), sol, lu=lu, &
        solver_options=so )

      call transfer_data ( mesh, problem1=problem, &
        sysvector1=sol, vector2=xnp1, degfd1=[1], degfd2=[i] )

    end do

    call delete(lu)

  end subroutine build_solve_interface


! move the interface position (time discretized step)

  subroutine move_interface

    if ( timeint == 1 .or. ( step == 1 .and. timeint > 1 ) ) then

!     Euler

      call copy ( xnp1, xn )

      coefficients%i(11) = 1

    else if ( timeint == 2  ) then

!     second order Gear with prediction

      call copy ( xn, xnm1 )
      call copy ( xnp1, xn )

!     2nd order prediction of interface coordinates
      mesh%coor = transpose( reshape( 2*xn%u - xnm1%u, [ndim,mesh%nnodes] ) )

      coefficients%i(11) = 2

    end if

    call build_solve_interface

!   Updating of interface mesh coordinates with new xnp1
    mesh%coor = transpose( reshape( xnp1%u, [ndim,mesh%nnodes] ) )

  end subroutine move_interface

end program interface_tracking6
