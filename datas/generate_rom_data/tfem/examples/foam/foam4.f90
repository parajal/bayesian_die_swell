! Stokes problem on a cubic domain with a spherical bubble (void).
! There is surface tension on the bubble boundary.
! Inside the bubble a constant pressure is imposed.
! Triperiodicity is imposed and the domain can freely stretch (without force)
! in x-, y- and z-direction.
! A time stepping is performed where the nodes of the mesh are updated
! according to a Lagrangian scheme (forward Euler).
! Gmsh mesh.

program foam4

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m
  use subs3_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,         & ! P2 velocities
    pintpl = 2,         & ! P1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 15,         & ! 15-point integration of tetrahedra
    gaussb = 6            ! 6-point integration of triangles on the boundary

  real(dp), parameter :: &
    lx = 4._dp,         & ! length in x of the domain (should be = mesh size)
    ly = 4._dp,         & ! length in y of the domain (should be = mesh size)
    lz = 4._dp            ! length in z of the domain (should be = mesh size)

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
  call read_mesh_gmsh ( mesh, filename='mesh4.msh', ndim=3, physgeom=.true. )

! local node numbering of surface 3 to match local node numbering of surface 1

  call add_to_mesh ( mesh, matchingsurface=[3,1], replace=3, &
    displacement=[0._dp,ly,0._dp] )

! local node numbering of surface 4 to match local node numbering of surface 2

  call add_to_mesh ( mesh, matchingsurface=[4,2], replace=4, &
    displacement=[-lx,0._dp,0._dp] )

! local node numbering of surface 6 to match local node numbering of surface 5

  call add_to_mesh ( mesh, matchingsurface=[6,5], replace=6, &
    displacement=[0._dp,0._dp,-lz] )

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

! define essential bc in one point for the velocity

  call define_essential ( mesh, input_probdef, point=1, physq=physqvel )

! force free contraints

  call define_constraint ( mesh, input_probdef, num=c1, &
    physq=physqvel, surface1=2, surface2=4, discretization='collocation', &
    naddunknowns=1 )

  call define_constraint ( mesh, input_probdef, num=c2, &
    physq=physqvel, surface1=3, surface2=1, discretization='collocation', &
    naddunknowns=1, excludesurfaces=[2] )

  call define_constraint ( mesh, input_probdef, num=c3, &
    physq=physqvel, surface1=5, surface2=6, discretization='collocation', &
    naddunknowns=1, excludesurfaces=[2,3] )

  call problem_definition ( input_probdef, mesh, problem )

! create subscript for velocity

  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )

! create system vectors (solution and right-hand side)
  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill boundary conditions

! set velocity in single point
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqvel, value=0._dp )

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

!   force free periodic constraints
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=c1, elemsub=elementc_x, addmatvec=.true. )
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=c2, elemsub=elementc_y, addmatvec=.true. )
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=c3, elemsub=elementc_z, addmatvec=.true. )

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

end program foam4
