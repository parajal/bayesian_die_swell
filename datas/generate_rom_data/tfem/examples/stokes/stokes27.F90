#if HSL_EXTRA

! Stokes problem on a 3D square cavity. Only half of the domain is solved
! due to symmetry conditions.
! Iterative solver using block preconditioning (using pressure mass matrix).

program stokes27

  use tfem_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m
  use hsl_solve1_mi20_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    nx = 20,            & ! number of elements in x-direction
    ny = 10,            & ! number of elements in y-direction
    nz = 20,            & ! number of elements in z-direction
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    itsolver = 5,       & ! iterative solver 5=BiCGSTAB, 8=GMRES
    gauss = 3             ! 3x3x3 integration of hexahedra

  real(dp), parameter :: &
    lx = 1._dp,          & ! size in x-direction
    ly = 0.5_dp,         & ! size in y-direction
    lz = 1._dp,          & ! size in z-direction
    eta = 1._dp            ! viscosity

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefp
  type(problem_t) :: problem, problemp
  type(sysmatrix_t) :: sysmatrix, sysmatrix_S
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, dudz
  type(plot_options_t) :: plot_options
  type(meshgen_options_t) :: meshgen_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options1_mi20_t) :: solver_options
  type(subscript_t) :: ssv, ssp
  type(prec1_mi20_t) :: prec

! variables

  integer :: presnod(8) = [1,3,9,7,19,21,27,25]

! additional factor for the pressure mass matrix
  real(dp) :: factor = 1.0_dp

! set pressure level?
  logical :: pressure_level = .false.

! do a second solve and keep the preconditioner data structure (no new setup)?
  logical :: second_solve = .false.


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, &
    lz=lz, elshape=14, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! plot curves, surfaces and mesh

  plot_options%viewpoint=[1.,0.8,0.4]
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig', surfaces=[3,4,6] )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(presnod,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, surface1=1, surface2=3, physq=1 )
  call define_essential ( mesh, input_probdef, surface1=4, physq=1, &
    degfd=[0,1,0] )
!  call define_essential ( mesh, input_probdef, surface1=4, physq=1 )
  call define_essential ( mesh, input_probdef, surface1=5, surface2=6, physq=1 )
  if ( pressure_level ) then
    call define_essential ( mesh, input_probdef, point=1, physq=2 )
  end if

  call problem_definition ( input_probdef, mesh, problem )

! subscript for partitioning the matrix into a velocity and pressure part
  call create_subscript ( mesh, problem, ssv, physqarr=[1], &
    essentialpart=.false. )
  call create_subscript ( mesh, problem, ssp, physqarr=[2], &
    essentialpart=.false. )

! problem definition pressure preconditioning

  call create_input_probdef ( mesh, input_probdefp )

  input_probdefp%elementdof(1)%a(:) = 0
  input_probdefp%elementdof(1)%a(presnod) = 1

  if ( pressure_level ) then
    call define_essential ( mesh, input_probdefp, point=1 )
  end if

  call problem_definition ( input_probdefp, mesh, problemp )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    surface1=1, surface2=3, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=4, physq=1, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=5, surface2=6, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=6, physq=1, degfd=1, value=1._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! create system matrix S

  call create_sysmatrix_structure ( sysmatrix_S, mesh, problemp )

  call create_sysmatrix_data ( sysmatrix_S )


 call tic


! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_Laplace_elem, coefficients=coefficients )


  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

! build (assemble) matrix from elements for pressure mass matrix

  call build_system ( mesh, problemp, sysmatrix_S, elemsub=pressure_mass_elem, &
    coefficients=coefficients, factormat=factor, buildvector=.false. )

  call set_solver_options ( solver_options, printlevel=2, eps_rel=1e-7_dp, &
    itsolver=itsolver, mgmres=40, cg_fixednumits=8 )

  !prec%control%v_iterations=2
  !prec%control%print_level=2

  call toc ( 'build' )

  call solve_system1_mi20 ( sysmatrix, rhsd, sol, ssv%s, ssp%s, sysmatrix_S, &
    prec=prec, solver_options=solver_options )

  call toc ('solve')

  if ( second_solve ) then

    call solve_system1_mi20 ( sysmatrix, rhsd, sol, ssv%s, ssp%s, sysmatrix_S, &
      prec=prec, solver_options=solver_options )

    call toc ('solve2')

  end if

  call delete ( prec )

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, dudz, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=3
  call derive_vector ( mesh, problem, dudz, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )


! write to a fig file for plotting

  plot_options%printlabels=.false.

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity, surfaces=[4,6] )
  call plot_points_curves ( plot_options, mesh, 'velocity.fig', append=.true. )

  call plot_color_fill ( plot_options, mesh, problem, 'velocity_color.fig', &
    vector=velocity, degfd=1, surfaces=[4] )
  call plot_points_curves ( plot_options, mesh, 'velocity_color.fig',  &
   append=.true. )

  call plot_color_contour ( plot_options, mesh, problem, 'dudz_contour.fig', &
    vector=dudz, surfaces=[3,4,6] )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure, surfaces=[3,4,6] )


! delete all data including all allocated memory

  call delete ( problem, problemp )
  call delete ( input_probdef, input_probdefp )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix, sysmatrix_S )
  call delete ( velocity, pressure, dudz )
  call delete ( coefficients )
  call delete ( oldvectors )

end program stokes27

#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on hsl_extra', &
    ' - add the libhsl3 library for linking', &
    ' - set preprocessing macro HSL_EXTRA in Mdefs.mk'
end
#endif

