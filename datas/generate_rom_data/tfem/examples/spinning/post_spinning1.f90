
! Post-processing for extrudate_swell1 (planar or axisymmetric Newtonian)

program post_swell1

  use tfem_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m

  implicit none


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysvector_t), target :: sol
  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors_d
  type(coefficients_t) :: coefficients

  integer :: ios


! print filename in title and date at footer

  plot_options%printfilename = .true.
  plot_options%printdateandtime = .true.


! read coefficients

  call read_coefficients ( coefficients, filename="coefficients.out" )


! read mesh

  call read_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call read_input_probdef ( mesh, input_probdef, filename='probdef.out' )

  call problem_definition ( input_probdef, mesh, problem )


! create system vectors for gradient/velocity/pressure (solution)

  call create_sysvector ( problem, sol )


! read data for post-processing

  open ( unit = 10, file='data.out', form='unformatted', iostat=ios, &
         status='old' )

  if ( ios /= 0 ) then
    write(*,'(/2a/)') 'Error: cannot open file data.out '
    stop
  end if

  read(10) sol%u

  close ( unit=10 )


! post-processing


! mesh

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! velocity vector

  call create_vector ( problem, velocity, physq=1 )
  call extract_physvector ( mesh, problem, sol, velocity )

  plot_options%scalevector=0.3

  plot_options%fontsize = 10

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  plot_options%rotatecolorbar=.true.
  plot_options%shiftbary=120

  call plot_color_fill ( plot_options, mesh, problem, 'velocity_color.fig', &
    vector=velocity, degfd=1 )

! write binary file for reading by streamfunction

  open ( unit=10, form='unformatted', file='velocity_bin.out' )

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

  call printtofile ( mesh, problem, 'curve=3.out', curve=3, sysvector=sol )

  call delete ( pressure )


! vorticity

  call create_vector ( problem, vorticity, vec=3 )

  coefficients%i(13)=5
  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors_d )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

  call delete ( vorticity )


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol )
  call delete ( oldvectors_d )

  call delete ( coefficients )

end program post_swell1
