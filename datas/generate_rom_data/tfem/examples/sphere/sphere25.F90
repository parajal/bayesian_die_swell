#if SKIT2 && METIS5

! 3D Stokes problem: a sphere moving with velocity -U on the center
! line of a container with a square cross section and a bottom (at infinity).
! The sphere is fixed and the container wall is moved in z-direction with
! velocity U (moving frame).
! At the ends fluid is moving with velocity U as well, representing zero
! flux in the stationary frame (bottom).
! Dirichlet boundary conditions for the velocities are used.
! Drag computations.
! Sparskit solver
! Metis renumbering
!

program sphere25

  use tfem_m
  use sk_solve_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m
  use metis5_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,         & ! P2 velocities
    pintpl = 2,         & ! P1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 15,         &  ! 15-point integration of tetrahedra
    gaussb = 6             ! 6-point integration of triangles

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t), target :: pressure, stress_tensor
  type(vector_t) :: velocity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_sk_t) :: solver_options

! variables

  real(dp), parameter :: &
    eta = 1._dp,  & ! viscosity
    U = 1._dp       ! velocity of the sphere

  real(dp) :: drag(3)

  integer :: presnod(4) = [1,3,5,10]


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gaussb ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

! read mesh

  call read_mesh_gmsh ( mesh, filename='mesh25.msh', physgeom=.true. )

  call fill_mesh_parts ( mesh )

  call renumber_metis ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2, &
    nphysqshifted=1 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(presnod,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar
  input_probdef%vec_elementdof(1)%a(:,4) = 6  ! tensor

  input_probdef%physq = [1,2]
  input_probdef%physqshifted = [2]

! define essential boundaries

! sphere
  call define_essential ( mesh, input_probdef, surface1=7, physq=1 )
! ends and walls
  call define_essential ( mesh, input_probdef, surface1=1, surface2=6, &
    physq=1 )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0

  call fill_sysvector ( mesh, problem, sol, &
    surface1=1, surface2=6, physq=1, degfd=3, value=U )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1, nvec=2 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

! set parameters for iterative solver

!  ILUk
  call set_solver_options ( solver_options, printlevel=2, type_prec=2, &
    maxmvm=300, itsolver=8, mgmres=30, prec_store=4._dp, preconditioner=5, &
    droptol=1e-4_dp, fillin=1._dp, eps_rel = 1e-7_dp, eps_abs = 1.e-3_dp )

! ILUt
!  call set_solver_options ( solver_options, printlevel=2, type_prec=2, &
!    maxmvm=300, itsolver=8, mgmres=30, prec_store=4._dp, preconditioner=1, &
!    droptol=1e-3_dp, fillin=1.5_dp, eps_rel = 1e-6_dp, eps_abs = 1.e-6_dp )

  call solve_system_sk ( sysmatrix, rhsd, sol, solver_options=solver_options )

  call delete ( sysmatrix )

  call create_vector ( problem, velocity, physq=1 )
  call extract_physvector ( mesh, problem, sol, velocity )

! post processing

  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, stress_tensor, vec=4 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, stress_tensor, &
    elemsub=stokes_stress_tensor, coefficients=coefficients, &
    oldvectors=oldvectors, elsurfaces=[7] ) ! only include surface 7

  oldvectors%v(1)%p => stress_tensor
  oldvectors%v(2)%p => pressure

! integrate drag force

  call integrate_boundary_elements ( mesh, problem, drag, &
    elemsub=stokes_drag, surface=7, coefficients=coefficients, &
    oldvectors=oldvectors )

  print *, 'drag force =', drag


! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='sphere.vtk' )

  call write_vector_vtk ( mesh, problem, filename='sphere.vtk', &
    dataname='velocity_vector', sysvector=sol, append=.true. )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, stress_tensor )
  call delete ( coefficients )
  call delete ( oldvectors )

end program sphere25

#else
  print '(4(a/),a)', &
    'To run this example:', &
    ' - compile add-ons sk_solve and metis5', &
    ' - install and compile sparskit2', &
    ' - add metis5 lib for linking', &
    ' - set preprocessing macro SKIT2 and METIS5 in Mdefs.mk'
end
#endif
