! Mixed Stokes problem for a contraction flow
! Q2/Q1 Taylor-Hood for velocity-pressure, using main mesh for interpolation.
! Qp for stresses using the blend mesh 1 for interpolation.
! Blended mesh:
!  main mesh: 9-node quadrilaterals
!  blend mesh 1: (p+1)^2-node high-order quadrilaterals, equidistant nodes

program mixed_stokes1

  use tfem_m
  use stokes_elements_m
  use mixed_stokes_elements_2D_m
  use io_utils_m
  use figplot_m
  use meshgen_contraction_m
  use functions_m
  use hsl_ma57_m

  implicit none

! constants

  logical, parameter :: &
    CR = .false.           ! use mesh in lecture computational rheology

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    tintpl = 20,        & ! Qp stresses
    p = 4,              & ! polynomial order stresses
    physqvel = 2,       & ! physical quantity nr of the velocities
    physqpress = 3,     & ! physical quantity nr of the pressures
    physqstress = 1,    & ! physical quantity nr of the stress
    gauss = max(p+1,3), & ! p+1 x p+1 integration of quads
    gaussb = max(p+1,3),& ! number of Gauss points on the boundary
    inttype = 3,        & ! numerical rules for integration
    nn0 = 9,            & ! number of nodes in main mesh element
    nn1 = (p+1)**2        ! number of nodes in blend mesh element

  integer, parameter :: &
    funcnr = 3,         & ! function number of traction natural boundary
    direction = 1         ! direction of the traction component

  real(dp), parameter :: &
    refine = 1.0_dp, & ! refinement factor for mesh
    eta = 1._dp       ! viscosity

! definitions

  type(mesh_t) :: meshQ2, meshQp, mesh, mesh_plot_Qp
  type(coefficients_t) :: coefficients
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: pressure, vorticity, divergence
  type(oldvectors_t) :: oldvectors
  type(plot_options_t) :: po
  type(solver_options_ma57_t) :: so
  type(ic_t) :: ic


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,          tintpl, 0,     0,  &
      physqvel, physqpress, physqstress, 0, gauss,  &
      gaussb ]
  coefficients%i(12:) = 0

  coefficients%i(16) = funcnr
  coefficients%i(17) = direction
  coefficients%i(40) = inttype
  coefficients%i(82) = p ! stress polynomial order

  coefficients%r = 0
  coefficients%r(1) = 0    ! zero direct viscosity
  coefficients%r(5) = eta  ! viscosity in mixed element

  coefficients%func => func


! Create blended mesh

  if ( CR ) then
    ic = ic_t ( elshape=6, H1=1._dp, H2=3.0_dp, refine=refine )
  else
    ic = ic_t ( elshape=6, refine=refine )
  end if

  call generate_mesh ( meshQ2 ) ! velocity-pressure mesh

  if ( CR ) then
    ic = ic_t ( elshape=102, H1=1._dp, H2=3.0_dp, p=p, l=0, refine=refine )
  else
    ic = ic_t ( elshape=102, p=p, l=0, refine=refine )
  end if

  call generate_mesh ( meshQp ) ! stress tensor mesh

  call mesh_convert ( meshQ2, blendmesh=meshQp, mesh=mesh ) ! blended mesh

  call fill_mesh_parts ( mesh )

! Create plot mesh for stresses

  call mesh_convert ( meshQp, mesh_plot_Qp, warn=.false., &
    elementshapes='spectraltolinear' )

  call fill_mesh_parts ( mesh_plot_Qp )

! output mesh (print and plot)

  call printinfo ( mesh, printlevel=6 )

  po = plot_options_t ( fontsize=5 )
  call plot_points_curves ( po, mesh, 'curvesQ2.fig', blend=0 )
  call plot_points_curves ( po, mesh, 'curvesQp.fig', blend=1 )

  call write_mesh_gmsh ( mesh, filename='meshQ2.msh', blend=0 )
  call write_mesh_gmsh ( mesh, filename='meshQp.msh', blend=1 )

  call delete ( meshQ2, meshQp )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=3 )

  input_probdef%vec_elementdof(1)%a = 0           ! initialize to zero

  input_probdef%vec_elementdof(1)%a(1:nn0,1) = 2         ! velocity (Q2 mesh)
  input_probdef%vec_elementdof(1)%a([1,3,5,7],2) = 1     ! pressure (Q2 mesh)
  input_probdef%vec_elementdof(1)%a(nn0+1:nn0+nn1,3) = 3 ! stress (Qp mesh)
  input_probdef%vec_elementdof(1)%a(1:nn0,4) = 1         ! scalar (Q2 mesh)

  input_probdef%physq = [3,1,2]

! inflow, outflow, center line
  call define_essential ( mesh, input_probdef, curves=[20,22,24], &
    physq=physqvel, degsfd=[2] )

! wall
  call define_essential ( mesh, input_probdef, curves=[23], physq=physqvel )

  call problem_definition ( input_probdef, mesh, problem )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0


! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )


! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress], &
    coefficients=coefficients )

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=mixed_stokes_elem, addmatvec=.true., coefficients=coefficients, &
    physqrow=[physqstress,physqvel], physqcol=[physqstress,physqvel] )

  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., zeromatvec=.true., coefficients=coefficients, &
    physqrow=[physqstress], physqcol=[physqpress] )

  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., zeromatvec=.true., coefficients=coefficients, &
    physqrow=[physqpress], physqcol=[physqstress] )

  call check (sysmatrix)

  call add_boundary_elements ( mesh, problem, rhsd, curve=24, &
    physq=[physqvel], elemsub=stokes_natboun_curve, coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=so )


! post processing solution

  call write_vector_vtk ( mesh, problem, sysvector=sol, &
    dataname='velocity', filename='mixed_stokes1.vtk', physq=physqvel )

  call write_tensor_vtk ( mesh_plot_Qp, problem, sysvector=sol, &
    dataname='stress', filename='mixed_stokes1_stress.vtk', &
    physq=physqstress, offset=mesh%nnodes_blend(2) )

  if ( CR ) then

    call print_solution_on_curve25

  end if


! post processing using derive

  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, vorticity, vec=4 )
  call create_vector ( problem, divergence, vec=4 )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5
  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=7
  call derive_vector ( mesh, problem, divergence, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! write derivatives to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, &
    dataname='pressure', filename='mixed_stokes1.vtk', append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=vorticity, &
    dataname='vorticity', filename='mixed_stokes1.vtk', append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=divergence, &
    dataname='divergence', filename='mixed_stokes1.vtk', append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh, mesh_plot_Qp )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( pressure, vorticity, divergence )
  call delete ( coefficients )
  call delete ( oldvectors )

contains

  subroutine generate_mesh ( mesh )

    type(mesh_t), intent(out) :: mesh

    call contraction2D ( ic, mesh )

    call add_to_mesh ( mesh, curve=[1,8,13,16,19] )      ! c22 center line
    call add_to_mesh ( mesh, curve=[21,18,15,11,12,6] )  ! c23 wall
    call add_to_mesh ( mesh, curve=[7,4] )               ! c24 inflow boundary
    call add_to_mesh ( mesh, curve=[-10,-15] )           ! c25 curve for plot

  end subroutine generate_mesh

  subroutine print_solution_on_curve25

    type(vector_t) :: velocity, pressure_sol, stress

    call create_vector ( problem, velocity, physq=physqvel)
    call extract_physvector ( mesh, problem, sol, velocity )
    call printtofile ( mesh, problem, filename='velocity_on_curve25.out', &
      curve=25, vector=velocity )
    call delete ( velocity )

    call create_vector ( problem, pressure_sol, physq=physqpress )
    call extract_physvector ( mesh, problem, sol, pressure_sol )
    call printtofile ( mesh, problem, filename='pressure_on_curve25.out', &
      curve=25, vector=pressure_sol )
    call delete ( pressure_sol )

    call create_vector ( problem, stress, physq=physqstress )
    call extract_physvector ( mesh, problem, sol, stress )
    call printtofile ( mesh, problem, filename='stress_on_curve25.out', &
      curve=25, vector=stress )
    call delete ( stress )

  end subroutine print_solution_on_curve25

end program mixed_stokes1

