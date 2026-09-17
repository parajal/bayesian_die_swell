! Tracking of an initially circular interface in a flow given by a function.
! Interface movement is Lagrangian.
! No tangential grid deformation (u_t=0).
! Only the upper half is solved using the full system with symmetry conditions.
! Therefore only elongational flow (vfuncnr=1) is applicable.
! tfem mesh made using mesh_skeleton.
! Error output.

program interface_tracking13

  use tfem_m
  use math_defs_m
  use hsl_ma57_m
  use interface_tracking_elements_m
  use functions_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  integer, parameter ::  &
!    uintpl = 2,          & ! interpolation of position on the interface
    uintpl = 6,          & ! interpolation of position on the interface
    gauss = 4,           & ! number of integration points of line elements
    ndim = 2,            & ! space dimension
!    nelem = 100,         & ! number of elements
    nelem = 50,          & ! number of elements
    vfuncnr=1,           & ! function number for the velocity field u
    vtkevery = 1,        & ! vtk file every vtkevery steps. -1: means none
    timeint=2,           & ! time integration scheme
    numtimesteps=20       ! number of time steps

  real(dp), parameter :: &
    deltat = 5.e-2_dp   ! time step

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(vector_t) :: xnp1, error
  type(vector_t), target :: xn, xnm1
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors
  type(plot_options_t) :: plot_options
  type(solver_options_ma57_t) :: so


  character(len=30) :: filename

  integer :: i, step, post, elem, node
  real(dp) :: points(2,2)


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(3) = gauss
  coefficients%i(5) = 3  ! standard Gauss-Legendre (numerical tables)
  coefficients%i(6) = 1  ! velocity given by a function
  coefficients%i(7) = vfuncnr
  coefficients%i(11) = timeint
  coefficients%i(12) = 1 ! full vector system

  coefficients%r = 0
  coefficients%r(9) = deltat

  coefficients%vfunc1(1)%p => vfunc_2D

  if ( uintpl == 2 ) then

!   create mesh of a circular curve discretized with 2-node line elements

    call mesh_skeleton ( mesh, nnodes=nelem+1, nelem=nelem, elshape=1, &
      ndim=ndim )

    do elem = 1, mesh%nelem
      mesh%topology(1)%a(:,elem) = [ elem, elem+1 ]
    end do

  else if ( uintpl == 6 ) then

!   create mesh of a circular curve discretized with 3-node line elements

    call mesh_skeleton ( mesh, nnodes=2*nelem+1, nelem=nelem, elshape=2, &
      ndim=ndim )

    do elem = 1, mesh%nelem
      mesh%topology(1)%a(:,elem) = [ 2*elem-1, 2*elem, 2*elem+1 ]
    end do

  end if

  do node = 1, mesh%nnodes
    mesh%coor(node,1) = R0 * cos((node-1)*pi/(mesh%nnodes-1))
    mesh%coor(node,2) = R0 * sin((node-1)*pi/(mesh%nnodes-1))
  end do

  points(1,:) = [ R0, 0._dp ]
  points(2,:) = [ -R0, 0._dp ]

  call add_to_mesh ( mesh, points=points )

  call fill_mesh_parts ( mesh )

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )
  call write_mesh_vtk ( mesh, 'mesh.vtk' )

  call printinfo ( mesh, printlevel=4 )

! output initial coordinates

  write(filename,'(a,i4.4,a)') 'coor', 0, '.out'
  open ( unit=11, file=filename, recl=300 )
  write(unit=11,fmt='(i0,2e16.8)') (i, mesh%coor(i,:), i=1,mesh%nnodes)
  close ( unit=11 )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  input_probdef%elementdof(1)%a = ndim
  input_probdef%vec_elementdof(1)%a(:,1) = ndim  ! 2D vector
  input_probdef%vec_elementdof(1)%a(:,2) = 1  ! error

! symmetry conditions in P1 and P2

  call define_essential ( mesh, input_probdef, point=1, degfd=[0,1] )
  call define_essential ( mesh, input_probdef, point=2, degfd=[0,1] )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

  call fill_sysvector ( mesh, problem, sol, point=1, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=2, degfd=2, value=0._dp )

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
    call create_oldvectors ( oldvectors, nvec=3 )
    oldvectors%v(3)%p => xn
  else if ( timeint == 2 ) then
    call create_oldvectors ( oldvectors, nvec=4 )
    oldvectors%v(3)%p => xn
    oldvectors%v(4)%p => xnm1
  end if

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )


! open monitor files

  open ( unit=10, file='error1.out', recl=300 )


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

!       error
        call write_scalar_vtk ( mesh, problem, filename=filename, &
          dataname='error', vector=error )

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
  call delete ( xn, xnp1 )
  if ( timeint == 2 ) call delete ( xnm1 )
  call delete ( sysmatrix )
  call delete ( coefficients )

contains

! build and solve for xnp1 on interface

  subroutine build_solve_interface

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, sysvector=rhsd, &
      elemsub=interface_tracking_elem_lagrange, coefficients=coefficients, &
      oldvectors=oldvectors, order='DN' )

    so%integer_storage=1.8_dp

    call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=so )

    call transfer_data ( mesh, problem1=problem, sysvector1=sol, vector2=xnp1 )

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

end program interface_tracking13
