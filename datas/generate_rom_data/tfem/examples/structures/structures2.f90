! A simple 3D truss problem.

program structures2

  use tfem_m
  use hsl_ma57_m
  use structures_elements_m
  use structures_post_m
  use io_utils_m
  use figplot_m

  implicit none


! constants

  real(dp), parameter :: &
    Emod = 1.e11_dp, &  ! E modulus
    A = 1.e-4_dp, &     ! A area
    F = -1.e4_dp        ! Force


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd, reacf, load
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0

  coefficients%r = 0
  coefficients%r(1) = A
  coefficients%r(2) = Emod


! read mesh

  call read_mesh_gmsh ( mesh, filename='mesh2.msh' )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, filename='curves.fig' )
  call plot_mesh ( plot_options, mesh, filename='mesh.fig' )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  input_probdef%elementdof(1)%a = 3       ! displacement vector
  input_probdef%vec_elementdof(1)%a = 1   ! scalar quantity
  input_probdef%vec_elementdof(1)%a = 3   ! vector quantity

  call define_essential ( mesh, input_probdef, point=1 )
  call define_essential ( mesh, input_probdef, point=8 )
  call define_essential ( mesh, input_probdef, point=6, degfd=[0,1] )
  call define_essential ( mesh, input_probdef, point=11, degfd=[0,1] )

  call problem_definition ( input_probdef, mesh, problem )

  call printinfo ( mesh, printlevel=4 )
  call printinfo ( problem, printlevel=4 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd, reacf, load )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=8, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=6, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=11, degfd=2, value=0._dp )

! fill rhs vector with discrete forces

  load%u = 0

  call fill_sysvector ( mesh, problem, load, point=2, degfd=2, value=F )

  call copy ( load, rhsd )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=truss_elem1, &
    order='ND', coefficients=coefficients, buildvector=.false. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

! post processing

  call postprocessing_truss ( mesh, problem, coefficients=coefficients, &
    mshfilename='structures2.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures2

