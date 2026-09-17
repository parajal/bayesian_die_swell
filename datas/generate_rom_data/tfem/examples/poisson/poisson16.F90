#if HSL_EXTRA

! Poisson problem on a unit cube with Dirichlet boundary conditions and
! natural boundary conditions
! Iterative solver CG with algebraic multigrid preconditioner.
! Change the default parameters of the AMG preconditioner.
! Do a second solve with the same matrix keeping the preconditioner setup.

program poisson16

  use tfem_m
  use hsl_solve_mi20_m
  use poisson_elements_m
  use poisson_functions_m
  use figplot_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,         & ! scalar interpolation
    gauss = 3,          & ! 3x3x3 integration of quads
    gaussb = 3,         & ! 3x3 point integration of boundary elements
    nx=10,              & ! number of elements in x
    ny=10,              & ! number of elements in y
    nz=10,              & ! number of elements in z
    funcnr=10,          & ! function number for the right-hand side
    vfuncnr=3             ! function number for the flux vector

  real(dp), parameter :: &
    alpha = 1._dp     ! diffusion coefficient

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd, solexact
  type(plot_options_t) :: plot_options
  type(meshgen_options_t) :: meshgen_options
  type(coefficients_t) :: coefficients
  type(solver_options_mi20_t) :: solver_options
  type(prec_mi20_t) :: prec


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(10:12) = [ gauss, gaussb, funcnr ]
  coefficients%i(16) = vfuncnr

  coefficients%r(1) = alpha
  coefficients%r(2:) = 0

  coefficients%func => func
  coefficients%vfunc => vfunc

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, regionshape=4, &
    elshape=14 )

  call hexahedron ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! plot curves, surfaces and mesh

  plot_options%viewpoint=[1.,0.8,0.4]
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig', surfaces=[3,4,6] )

! problem definition

  call create_input_probdef ( mesh, input_probdef )

  input_probdef%elementdof(1)%a = 1

  call define_essential ( mesh, input_probdef, surfaces=[1,2,4,5,6] )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    surfaces=[1,2,4,5,6], func=func, funcnr=9 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=poisson_elem, &
    coefficients=coefficients )

  call add_boundary_elements ( mesh, problem, rhsd, surface=3, &
    elemsub=poisson_natboun, coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call set_solver_options ( solver_options, printlevel=2, move_matrix=.true., &
    eps_rel=1e-7_dp, itsolver=1 )

! set control options for hsl_mi20 (AMG preconditioner)

  prec%control%one_pass_coarsen=.true.
  prec%control%c_fail=2

  call solve_system_mi20 ( sysmatrix, rhsd, sol, prec=prec, &
    solver_options=solver_options )

! Do another solve with the same preconditioning matrix setup and
! using the old solution as the initial estimate of the solution.
! In this case you get a more accurate solution, but it is also possible
! to change the rhsd (but the matrix stays the same).

  call solve_system_mi20 ( sysmatrix, rhsd, sol, prec=prec, &
    solver_options=solver_options, initsol=.true. )

  call delete ( prec, sysmatrix=sysmatrix )

! create exact solution

  call create_sysvector ( problem, solexact )

  call fill_sysvector ( mesh, problem, solexact, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=9 )

! print maximum difference of sol-solexact to standard output

  print *, maxval( abs(sol%u -solexact%u) )

! plot solution on surface 3

  plot_options%printlabels=.false.

  call plot_color_fill ( plot_options, mesh, problem, filename='color.fig', &
    surfaces=[3], sysvector=sol )
  call plot_points_curves ( plot_options, mesh, 'color.fig', append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, solexact, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program poisson16

#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on hsl_extra', &
    ' - add the libhsl3 library for linking', &
    ' - set preprocessing macro HSL_EXTRA in Mdefs.mk'
end
#endif

