#if HSL_EXTRA

! Stokes+mass problem on a unit square with Dirichlet boundary conditions.
! Lid-driven cavity flow.
! Physical quantities.
! Iterative solver with block preconditioning (using both pressure mass and
! diffusion matrix). M_1^{-1}+M_2^{-1} form for inverse Schur complement.

program stokes31

  use tfem_m
  use hsl_solve2_mi20_m
  use stokes_elements_m
  use io_utils_m
  use figplot_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    itsolver = 5,       & ! iterative solver 5=BiCGSTAB, 8=GMRES
    nx=50,              & ! number of elements in x
    ny=50                 ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp     ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefp
  type(problem_t) :: problem, problemp
  type(sysmatrix_t) :: sysmatrix, sysmatrix_M1, sysmatrix_M2
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options2_mi20_t) :: solver_options
  type(subscript_t) :: ssv, ssp
  type(prec2_mi20_t) :: prec

! set pressure level?
  logical :: pressure_level = .false.

! factor a of the mass matrix: F_v = a * M_v + eta * A_v
  real(dp) :: factor_a = 1.e4_dp


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

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )
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

  input_probdefp%elementdof(1)%a(:) = [1,0,1,0,1,0,1,0,0]

  if ( pressure_level ) then
    call define_essential ( mesh, input_probdefp, point=1 )
  end if

  call problem_definition ( input_probdefp, mesh, problemp )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=1, degfd=1, value=1._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! create system matrix M1 and M2

  call create_sysmatrix_structure ( sysmatrix_M1, mesh, problemp )

  call create_sysmatrix_data ( sysmatrix_M1 )

  call copy ( sysmatrix_M1, sysmatrix_M2 )

! build (assemble) matrix and vector from elements

  call tic

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_Laplace_elem, &
    coefficients=coefficients )

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_mass_elem, coefficients=coefficients, &
    factormat=factor_a, addmatvec=.true., buildvector=.false. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

! build (assemble) matrix from elements for pressure mass and diffusion matrix

  call build_system ( mesh, problemp, sysmatrix_M1, &
    elemsub=pressure_mass_elem, factormat=1._dp/eta, &
    coefficients=coefficients, buildvector=.false. )

  call build_system ( mesh, problemp, sysmatrix_M2, &
    elemsub=pressure_diffusion_elem, factormat=1._dp/factor_a, &
    coefficients=coefficients, buildvector=.false. )


  call set_solver_options ( solver_options, printlevel=2, eps_rel=1e-7_dp, &
    itsolver=itsolver, mgmres=40, cg_fixednumits=8 )

  prec%control_F%c_fail=2
  !prec%control_M2%v_iterations=2
  !prec%control%print_level=2

  call toc ( 'build' )

  call solve_system2_mi20 ( sysmatrix, rhsd, sol, ssv%s, ssp%s, sysmatrix_M1, &
    sysmatrix_M2, prec=prec, solver_options=solver_options )

  call toc ('solve')

  call delete ( prec )

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, vorticity, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a fig file for plotting

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  call plot_color_fill ( plot_options, mesh, problem, 'pressure_color.fig', &
    vector=pressure )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )

  call plot_color_fill ( plot_options, mesh, problem, 'vorticity_color.fig', &
    vector=vorticity )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

! write binary file for reading by streamfunction

  open(unit=10,form='unformatted',file='velocity_bin.out')

  write(10) velocity%u

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem, problemp )
  call delete ( input_probdef, input_probdefp )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity )
  call delete ( sysmatrix, sysmatrix_M1, sysmatrix_M2 )
  call delete ( coefficients )
  call delete ( oldvectors )

end program stokes31

#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on hsl_extra', &
    ' - add the libhsl3 library for linking', &
    ' - set preprocessing macro HSL_EXTRA in Mdefs.mk'
end
#endif

