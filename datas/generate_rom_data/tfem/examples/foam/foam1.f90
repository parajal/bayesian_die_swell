! Stokes problem on a square domain with a circular bubble (void) at the center.
! There is surface tension on the bubble boundary.
! Inside the bubble a constant pressure is imposed and the domain boundaries
! remain straight, but can move according to a zero imposed force.
! A time stepping is performed where the nodes of the mesh are updated
! according to a Lagrangian scheme (forward Euler).
! Gmsh mesh.

program foam1

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
    gauss = 6,          & ! 6-point Gauss integration of triangles
    gaussb = 3            ! 3 point integration of boundary elements

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

! variables

  integer :: &
    vtkevery = 1,         & ! vtk file every vtkevery steps. -1: means none
    numtimesteps = 220      ! number of time steps

  real(dp) ::  &
    eta = 1._dp,          & ! fluid viscosity
    gammac = 1._dp,       & ! surface tension coefficient
    bubble_pressure = 2._dp, &  ! pressure inside the bubble
    deltat = 5.e-3_dp       ! time step

  integer :: step, post, c1, c2


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
  call read_mesh_gmsh ( mesh, filename='mesh1.msh', ndim=2, physgeom=.true. )

  call fill_mesh_parts ( mesh )

! plot mesh

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  plot_options%fontsize = 10
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [ 2,2,2,2,2,2,   &  ! velocity
                  1,0,1,0,1,0,   &  ! pressure
                  1,1,1,1,1,1,   &  ! scalar
                  3,3,3,3,3,3 ], &  ! tensor
                    [6,4] )

  input_probdef%physq = [1,2]
  input_probdef%probnr = 1

! define essential boundaries

! cell boundaries
  call define_essential ( mesh, input_probdef, curve1=1, degfd=[0,1], &
    physq=physqvel )
  call define_essential ( mesh, input_probdef, curve1=4, degfd=[1,0], &
    physq=physqvel )

! force free contraints
  call define_constraint ( mesh, input_probdef, num=c1, &
    physq=physqvel, curve1=2, nodedof=1, discretization='collocation', &
    naddunknowns=1 )
  call define_constraint ( mesh, input_probdef, num=c2, &
    physq=physqvel, curve1=3, nodedof=1, discretization='collocation', &
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
    curve1=1, physq=physqvel, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=4, physq=physqvel, degfd=1, value=0._dp )

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

    call add_boundary_elements ( mesh, problem, rhsd, curve=5, &
      physq=[physqvel], elemsub=stokes_natboun_normal, &
      coefficients=coefficients )

!   surface tension on the bubble boundary

    call add_boundary_elements ( mesh, problem, rhsd, &
      elemsub=surface_tension_curve, curve=5, &
      coefficients=coefficients, physq=[physqvel] )

!   force free boundaries
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=c1, elemsub=elementc1, addmatvec=.true. )
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=c2, elemsub=elementc2, addmatvec=.true. )

    call check_filled_sysmatrix ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve gradient/velocity/pressure problem

    call solve_system_ma57 ( sysmatrix, rhsd, sol )

!   update mesh coordinates (forward Euler, Lagrangian)

    mesh%coor = mesh%coor + &
           transpose ( reshape ( sol%u(vel%s), [2,mesh%nnodes] ) ) * deltat

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

! deformed mesh

  plot_options%fontsize = 10
  call plot_mesh ( plot_options, mesh, 'mesh_deformed.fig' )

! plot velocity field

  call plot_color_contour ( plot_options, mesh, problem, 'u.fig', &
    physq=1, degfd=1, sysvector=sol )
  call plot_color_contour ( plot_options, mesh, problem, 'v.fig', &
    physq=1, degfd=2, sysvector=sol )

! plot pressure field

  call create_vector ( problem, pressure, vec=3 )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  call plot_color_contour ( plot_options, mesh, problem, 'p.fig', &
    vector=pressure )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( pressure )
  call delete ( oldvectors )
  call delete ( coefficients )

end program foam1
