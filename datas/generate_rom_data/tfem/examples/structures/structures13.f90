! A simple 2D truss problem
! (question 1 of exam on 8-11-2019 of course Mechanics 4RA00).

program structures13

  use tfem_m
  use hsl_ma57_m
  use structures_elements_m
  use structures_post_m
  use io_utils_m
  use figplot_m

  implicit none


! constants

  real(dp), parameter :: &
    Emod = 1._dp, &  ! E modulus
    A = 1._dp, &     ! A area
    F = 1_dp         ! Force


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: load, rhsd, reacf
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0

  coefficients%r = 0
  coefficients%r(1) = A
  coefficients%r(2) = Emod


! read mesh

  call read_mesh_gmsh ( mesh, filename='bridge.msh', ndim=2 )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, filename='curves.fig' )
  call plot_mesh ( plot_options, mesh, filename='mesh.fig' )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  input_probdef%elementdof(1)%a = 2       ! displacement vector
  input_probdef%vec_elementdof(1)%a = 1   ! scalar quantity
  input_probdef%vec_elementdof(1)%a = 2   ! vector quantity

  call define_essential ( mesh, input_probdef, point=2, degfd=[0,1] )
  call define_essential ( mesh, input_probdef, point=7, degfd=[1,0] )
  call define_essential ( mesh, input_probdef, point=13 )

  call problem_definition ( input_probdef, mesh, problem )

  call printinfo ( mesh, printlevel=4 )
  call printinfo ( problem, printlevel=4 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd, reacf, load )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=2, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=7, degfd=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=13, value=0._dp )

! fill rhs vector with discrete forces

  load%u = 0

  call fill_sysvector ( mesh, problem, load, point=7, degfd=2, value=-F )
  call fill_sysvector ( mesh, problem, load, point=8, degfd=2, value=-F )
  call fill_sysvector ( mesh, problem, load, point=9, degfd=2, value=-F )
  call fill_sysvector ( mesh, problem, load, point=10, degfd=2, value=-F )
  call fill_sysvector ( mesh, problem, load, point=11, degfd=2, value=-F )
  call fill_sysvector ( mesh, problem, load, point=12, degfd=2, value=-F )

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
    mshfilename='structures13.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures13

