! steady convection-reaction equation using DGFEM.
! periodic boundaries

program convectionreaction2

  use tfem_m
  use hsl_ma41_m
  use convecreac_elements_m
  use functions_m
  use figplot_m

  implicit none

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd, solexact
  type(vector_t) :: vector
  type(plot_options_t) :: plot_options

  integer :: i, nx, ny

! mesh generation

  meshgen_options%regionshape = 1
  meshgen_options%lx = 2.0_dp
  meshgen_options%ly = 2.0_dp
  meshgen_options%ox = 0.3_dp
  meshgen_options%oy = 0.3_dp

  nx = 20
  ny = 20

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-3] )  ! curve 5
  call add_to_mesh ( mesh, curve=[-4] )  ! curve 6

  call glue_mesh ( mesh, curve1=1, curve2=5 )
  call glue_mesh ( mesh, curve1=2, curve2=6 )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = [(0,i=1,8),9]
  input_probdef%vec_elementdof(1)%a = reshape ( [(1,i=1,9)], [9,1] )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors sol and rhsd

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )

  call create_sysmatrix_structure_dg ( sysmatrix, mesh, problem, &
    elemsubzeros=convecreac_sidelem_zeros )

  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build system

  funcnr = 1
  funcnrinflow = -1  ! <0: will stop in function func when there is inflow
  vfuncnr = 1

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=convecreac_elem )

  call build_system_dg ( mesh, problem, sysmatrix, elemsub=convecreac_sidelem, &
    addmatrix=.true., elemsubzeros=convecreac_sidelem_zeros )

  call check ( sysmatrix )

! solve system

  call solve_system_ma41 ( sysmatrix, rhsd, sol )

! create exact solution

  call create_sysvector ( problem, solexact )

  funcnr = 2

  call fill_sysvector_elem ( mesh, problem, solexact, &
    elemsub=convecreac_fillexact )

! print maximum difference of sol-solexact to standard output

  print *, maxval( abs(sol%u -solexact%u) )

! create vector for plotting

  call create_vector ( problem, vector, vec=1, elementwise=.true. )

  solution = sol

  call derive_vector ( mesh, problem, vector, elemsub=convecreac_vecelem )

! plot solution using figplot

  call plot_color_fill ( plot_options, mesh, problem, filename='sol_f.fig', &
    vector=vector )
  call plot_color_contour ( plot_options, mesh, problem, filename='sol_c.fig', &
   vector=vector )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, solexact, rhsd )
  call delete ( sysmatrix )
  call delete ( vector )

end program convectionreaction2

