! diffuse-interface problem on a unit square
! constant viscosity
! Post-processing

program  diffuse_interface_post

  use tfem_m
  use diffuse_interface_elements_m
  use figplot_m
  use io_utils_m

  implicit none


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefdi
  type(problem_t) :: problem, problemdi
  type(sysvector_t), target :: sol
  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors_d
  type(coefficients_t) :: coefficients

  type(sysvector_t) :: soldi

  integer :: ios


! print filename in title and date at footer

  plot_options%printfilename = .true.
  plot_options%printdateandtime = .true.


! read coefficients

  call read_coefficients ( coefficients, filename="coefficients.out" )


! read mesh

  call read_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )


! problem definition of velocity/pressure

  call read_input_probdef ( mesh, input_probdef, filename='probdef.out' )

  call problem_definition ( input_probdef, mesh, problem )


! problem definition conformation tensor

  call read_input_probdef ( mesh, input_probdefdi, filename='probdefdi.out' )

  call problem_definition ( input_probdefdi, mesh, problemdi )


! create system vectors for velocity/pressure (solution)

  call create_sysvector ( problem, sol )

! create system vectors (solution ) for conformation and

  call create_sysvector ( problemdi, soldi )


! read data for post-processing

  open ( unit = 10, file='data.out', form='unformatted', iostat=ios, &
         status='old' )

  if ( ios /= 0 ) then
    write(*,'(/2a/)') 'Error: cannot open file data.out '
    stop
  end if

  read(10) sol%u
  read(10) soldi%u

  close ( unit=10 )


! post-processing


! mesh

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! velocity vector

  call create_vector ( problem, velocity, physq=1 )
  call extract_physvector ( mesh, problem, sol, velocity )

  plot_options%scalevector=1e-1

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  call plot_color_fill ( plot_options, mesh, problem, 'velocity_color.fig', &
    vector=velocity, degfd=1 )

! write binary file for reading by streamfunction

  open ( unit=10, form='unformatted', file='velocity.out' )

  write( unit=10 ) velocity%u

  close ( unit=10 )

  call delete ( velocity )


! fill oldvectors for stokes problem

  call create_oldvectors ( oldvectors_d, nsysvec=1 )

  oldvectors_d%s(1)%p => sol


! pressure

  call create_vector ( problem, pressure, vec=3 )

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors_d )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )

  call delete ( pressure )


! vorticity

  call create_vector ( problem, vorticity, vec=3 )

  coefficients%i(13)=5
  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors_d )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

  call delete ( vorticity )


! Plotting c and mu

  call plot_color_contour ( plot_options, mesh, problemdi, &
    'c_contour.fig', sysvector=soldi, physq=1 )
  call plot_color_fill ( plot_options, mesh, problemdi, &
    'c_color_fill.fig', sysvector=soldi, physq=1 )
  call plot_color_contour ( plot_options, mesh, problemdi, &
    'mu_contour.fig', sysvector=soldi, physq=2 )

! print data on curve to a file

  call printtofile ( mesh, problemdi, 'c_mu_curve1.out', curve=1, &
    sysvector=soldi )
  call printtofile ( mesh, problemdi, 'c_mu_curve2.out', curve=2, &
    sysvector=soldi )

! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol )
  call delete ( oldvectors_d )
  call delete ( problemdi )
  call delete ( input_probdefdi )
  call delete ( soldi )

  call delete ( coefficients )

end program diffuse_interface_post
