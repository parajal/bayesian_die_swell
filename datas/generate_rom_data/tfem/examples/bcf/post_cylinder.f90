! Flow around a cylinder confined between two plates for a Hookean dumbbell.

! Post-processing

program cylinder_post

  use tfem_m
  use bcf_elements_m
  use figplot_m
  use io_utils_m

  implicit none


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefq
  type(problem_t) :: problem, problemq
  type(sysvector_t), target :: sol
  type(oldvectors_t) :: oldvectors_d, oldvectors_dve
  type(coefficients_t) :: coefficients
  type(sysvector_t), allocatable, dimension(:), target :: solq

  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options

  type(vector_t) :: ctensor

  integer :: i, ios, nfield


! print filename in title and date at footer

  plot_options%printfilename = .true.
  plot_options%printdateandtime = .true.


! read coefficients

  call read_coefficients ( coefficients, filename="coefficients.out" )

! allocate the arrays for Q-vectors (not data!)

  nfield = coefficients%i(455)

  allocate ( solq(nfield) )

! read mesh

  call read_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call read_input_probdef ( mesh, input_probdef, filename='probdef.out' )

  call problem_definition ( input_probdef, mesh, problem )


! problem definition Q vectors

  call read_input_probdef ( mesh, input_probdefq, filename='probdefq.out' )

  call problem_definition ( input_probdefq, mesh, problemq )


! create system vectors for gradient/velocity/pressure (solution)

  call create_sysvector ( problem, sol )


! create system vectors (solution ) for Q-vectors

  call create ( problemq, solq )


! read data for post-processing

  open ( unit = 10, file='data.out', form='unformatted', iostat=ios, &
         status='old' )

  if ( ios /= 0 ) then
    write(*,'(/2a/)') 'Error: cannot open file data.out '
    stop
  end if

  read(10) sol%u
  read(10) (solq(i)%u, i=1,nfield)

  close ( unit=10 )


! post-processing


! mesh

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! velocity vector

  call create_vector ( problem, velocity, physq=2 )
  call extract_physvector ( mesh, problem, sol, velocity )

  plot_options%fontsize = 10

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  plot_options%rotatecolorbar=.true.
  plot_options%shiftbary=100

  call plot_color_fill ( plot_options, mesh, problem, 'velocity_color.fig', &
    vector=velocity, degfd=1 )

  plot_options%shiftbary = 30

  call printtofile ( mesh, problem, filename='vx_on_centerline.out', &
    curve=22, vector=velocity )

  call printtofile ( mesh, problem, filename='vx_cross_section.out', &
    curve=18, vector=velocity )

! write binary file for reading by streamfunction

  open ( unit=10, form='unformatted', file='velocity_bin.out' )

  write( unit=10 ) velocity%u

  close ( unit=10 )

  call delete ( velocity )


! fill oldvectors for stokes problem

  call create_oldvectors ( oldvectors_d, nsysvec=1 )

  oldvectors_d%s(1)%p => sol

! pressure

  call create_vector ( problem, pressure, vec=4 )

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors_d )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )

  call delete ( pressure )


! vorticity

  call create_vector ( problem, vorticity, vec=4 )

  coefficients%i(13)=5
  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors_d )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

  call delete ( vorticity )


! create oldvectors for viscoelastic problem

  call create_oldvectors ( oldvectors_dve, nsysvec1=1 )

  oldvectors_dve%s1(1)%p => solq


! structure tensor

  call create_vector ( problemq, ctensor, vec=3 )

  call derive_vector ( mesh, problemq, ctensor, &
    elemsub=deriv_structure_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_dve )

  call plot_color_contour ( plot_options, mesh, problemq, &
    'cxx_contour.fig', vector=ctensor, degfd=1 )
  call plot_color_contour ( plot_options, mesh, problemq, &
    'cxy_contour.fig', vector=ctensor, degfd=2 )
  call plot_color_contour ( plot_options, mesh, problemq, &
    'cyy_contour.fig', vector=ctensor, degfd=3 )

  call printtofile ( mesh, problemq, filename='c_on_centerline.out', &
    curve=22, vector=ctensor )

  call delete ( ctensor )


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol )
  call delete ( oldvectors_d, oldvectors_dve )
  call delete ( problemq )
  call delete ( input_probdefq )
  call delete ( solq )

  call delete ( coefficients )

  deallocate ( solq )

end program cylinder_post
