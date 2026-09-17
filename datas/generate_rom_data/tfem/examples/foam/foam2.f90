! Stokes problem on a cubic domain with a spherical bubble (void) at the center.
! There is surface tension on the bubble boundary.
! Inside the bubble a constant pressure is imposed and the domain boundaries
! remain straight, but can move according to a zero imposed force.
! A time stepping is performed where the nodes of the mesh are updated
! according to a Lagrangian scheme (forward Euler).
! Gmsh mesh.

program foam2

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m
  use subs1_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,         & ! P2 velocities
    pintpl = 2,         & ! P1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 15,         & ! 15-point integration of tetrahedra
    gaussb = 6            ! 6-point integration of triangles on the boundary

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: pressure
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors
  type(subscript_t) :: vel
  type(solver_options_ma57_t) :: so

! variables

  integer :: &
    vtkevery = 1,         & ! vtk file every vtkevery steps. -1: means none
    numtimesteps = 150      ! number of time steps

  real(dp) ::  &
    eta = 1._dp,          & ! fluid viscosity
    gammac = 0.5_dp,       & ! surface tension coefficient
    bubble_pressure = 5._dp, &  ! pressure inside the bubble
    deltat = 5.e-3_dp       ! time step

  integer :: presnod(4) = [1,3,5,10], step, post, c1, c2, c3


! fill coefficients
  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gaussb ]
  coefficients%i(12:) = 0
  coefficients%r(1:) = 0
  coefficients%r(1) = eta
  coefficients%r(19) = gammac
  coefficients%r(24) = bubble_pressure

! read mesh from gmsh output file
  call read_mesh_gmsh ( mesh, filename='mesh2.msh', ndim=3, physgeom=.true. )

  call fill_mesh_parts ( mesh )

! plot mesh

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(presnod,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar
  input_probdef%vec_elementdof(1)%a(:,4) = 6  ! tensor

  input_probdef%physq = [1,2]
  input_probdef%probnr = 1

! define essential boundaries

! cell boundaries
  call define_essential ( mesh, input_probdef, surface1=4, degfd=[1,0,0], &
    physq=physqvel )
  call define_essential ( mesh, input_probdef, surface1=1, degfd=[0,1,0], &
    physq=physqvel )
  call define_essential ( mesh, input_probdef, surface1=6, degfd=[0,0,1], &
    physq=physqvel )

! force free contraints
  call define_constraint ( mesh, input_probdef, num=c1, &
    physq=physqvel, surface1=2, nodedof=1, discretization='collocation', &
    naddunknowns=1 )
  call define_constraint ( mesh, input_probdef, num=c2, &
    physq=physqvel, surface1=3, nodedof=1, discretization='collocation', &
    naddunknowns=1 )
  call define_constraint ( mesh, input_probdef, num=c3, &
    physq=physqvel, surface1=5, nodedof=1, discretization='collocation', &
    naddunknowns=1 )

  call problem_definition ( input_probdef, mesh, problem )

! create subscript for velocity

  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )

! create system vectors (solution and right-hand side)
  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill boundary conditions

! cell boundaries
  call fill_sysvector ( mesh, problem, sol, &
    surface1=4, physq=physqvel, degfd=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=1, physq=physqvel, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=6, physq=physqvel, degfd=3, value=0._dp )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! time stepping

  open ( unit=10, file='size.out' ) ! output file for size of the foam

  post = 0

  do step = 1, numtimesteps

    print *, 'step = ', step

!   build (assemble) matrix and vector from elements

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
      coefficients=coefficients )

!   pressure on the inside of the bubble

    call add_boundary_elements ( mesh, problem, rhsd, surface=7, &
      physq=[physqvel], elemsub=stokes_natboun_normal, &
      coefficients=coefficients )

!   surface tension on the bubble boundary

    call add_boundary_elements ( mesh, problem, rhsd, &
      elemsub=surface_tension_surface, surface=7, &
      coefficients=coefficients, physq=[physqvel] )

!   force free boundaries
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=c1, elemsub=elementc1, addmatvec=.true. )
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=c2, elemsub=elementc2, addmatvec=.true. )
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=c3, elemsub=elementc3, addmatvec=.true. )

    call check_filled_sysmatrix ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve gradient/velocity/pressure problem

    so%real_storage = 1.5_dp

    call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=so )

!   update mesh coordinates (forward Euler, Lagrangian)

    mesh%coor = mesh%coor + &
           transpose ( reshape ( sol%u(vel%s), [3,mesh%nnodes] ) ) * deltat

!   output length of the domain

    write(10,*) step*deltat, mesh%coor(mesh%points(2),1)
    write(*,'(2(a,g0.4))') &
              'time = ', step*deltat, ' size = ', mesh%coor(mesh%points(2),1)

!   output to vtk

    if ( vtkevery > 0 ) then
      if ( mod(step,vtkevery) == 0 ) then

        print *, 'writing fields to vtk'

        post = post + 1

        call output_fields ( mesh, problem, coefficients, sol, post )

      end if
    end if


  end do

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( pressure )
  call delete ( oldvectors )
  call delete ( coefficients )

end program foam2
