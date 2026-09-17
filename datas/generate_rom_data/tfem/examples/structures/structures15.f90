! A simple 2D truss problem with two groups having a different cross section.
! Similar to structures14, but now includes gravity load (lumped to the nodes).

program structures15

  use tfem_m
  use math_defs_m
  use hsl_ma57_m
  use structures_elements_m
  use structures_post_m
  use io_utils_m
  use figplot_m

  implicit none


! constants

  real(dp), parameter :: &
    rho = 7850._dp,  & ! density
    g = 9.8_dp, & ! gravitational constant
    Emod = 200.e9_dp, &  ! E modulus
    A1 = pi/4*(0.2_dp**2-0.08_dp**2), &  ! A area outside
    A2 = pi/4*(0.2_dp**2-0.14_dp**2), &  ! A area inside
    F = 10000._dp * g  ! Force


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: load, rhsd, reacf
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients(2)

  integer :: i

! fill coefficients

  call create_coefficients ( coefficients(1), ncoefi=100, ncoefr=50 )

  coefficients(1)%i = 0
  coefficients(1)%i(10) = 0  ! constant distributed load (lumped to end nodes)


  coefficients(1)%r = 0
  coefficients(1)%r(1) = A1
  coefficients(1)%r(2) = Emod
  coefficients(1)%r(5:6) = [ 0._dp, -rho*g*A1 ]

  coefficients(2) = coefficients(1)
  coefficients(2)%r(1) = A2
  coefficients(2)%r(5:6) = [ 0._dp, -rho*g*A2 ]


! read mesh

  call read_mesh_gmsh ( mesh, filename='constructie14.msh', ndim=2 )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, filename='curves.fig' )
  call plot_mesh ( plot_options, mesh, filename='mesh.fig' )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  do i = 1, 2
    input_probdef%elementdof(i)%a = 2       ! displacement vector
    input_probdef%vec_elementdof(i)%a = 1   ! scalar quantity
    input_probdef%vec_elementdof(i)%a = 2   ! vector quantity
  end do

  call define_essential ( mesh, input_probdef, point=1 )
  call define_essential ( mesh, input_probdef, point=3, degfd=[0,1] )

  call problem_definition ( input_probdef, mesh, problem )

  call printinfo ( mesh, printlevel=4 )
  call printinfo ( problem, printlevel=4 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd, reacf, load )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=3, degfd=2, value=0._dp )

! fill rhs vector with discrete forces

  load%u = 0

  call fill_sysvector ( mesh, problem, load, point=9, degfd=1, value=F )

  call copy ( load, rhsd )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=truss_elem1, &
    order='ND', mcoefficients=coefficients, addvec=.true. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

! post processing

  call postprocessing_truss ( mesh, problem, mcoefficients=coefficients, &
    mshfilename='structures15.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures15

