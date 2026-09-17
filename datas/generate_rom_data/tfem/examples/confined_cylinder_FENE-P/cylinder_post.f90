! 2D viscoelastic problem: flow around a cylinder confined between two walls.
! Post-processing
! FENE-P model

program cylinder_post

  use tfem_m
  use viscoelastic_elements_m
  use figplot_m
  use io_utils_m

  implicit none


! definitions

  integer, parameter :: ncompc = 4

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc
  type(problem_t) :: problem, problemc
  type(sysvector_t), target :: sol
  type(vector_t) :: velocity, pressure, ctensor, vorticity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors_d, oldvectors_dve
  type(coefficients_t) :: coefficients

  type(sysvector_t), dimension(ncompc,1), target :: solc

  integer :: i, ios


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


! problem definition conformation tensor

  call read_input_probdef ( mesh, input_probdefc, filename='probdefc.out' )

  call problem_definition ( input_probdefc, mesh, problemc )


! create system vectors for gradient/velocity/pressure (solution)

  call create_sysvector ( problem, sol )

! create system vectors (solution ) for conformation and

  call create ( problemc, solc )


! read data for post-processing

  open ( unit = 10, file='data.out', form='unformatted', iostat=ios, &
         status='old' )

  if ( ios /= 0 ) then
    write(*,'(/2a/)') 'Error: cannot open file data.out '
    stop
  end if

  read(10) sol%u
  read(10) (solc(i,1)%u, i=1,ncompc)

  close ( unit=10 )


! post-processing

! mesh


! plot curves, surfaces

  plot_options%fontsize=5
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  plot_options%fontsize=10

  plot_options%printlabels=.false.


! mesh


  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! post processing

! velocity vector

  call create_vector ( problem, velocity, physq=2 )
  call extract_physvector ( mesh, problem, sol, velocity )

  plot_options%fontsize = 10

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

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


  plot_options%shiftbary = 30


! pressure

  call create_vector ( problem, pressure, vec=4 )

  oldvectors_d%s(1)%p => sol

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

  call create_oldvectors ( oldvectors_dve, nsysvec2=1 )

  oldvectors_dve%s2(1)%p => solc


! conformation tensor

  call create_vector ( problemc, ctensor, vec=3 )

  call derive_vector ( mesh, problemc, ctensor, &
    elemsub=deriv_conformation_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_dve )

  call plot_color_contour ( plot_options, mesh, problemc, &
    'cxx_contour.fig', vector=ctensor, degfd=1 )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'cxy_contour.fig', vector=ctensor, degfd=2 )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'cyy_contour.fig', vector=ctensor, degfd=3 )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'czz_contour.fig', vector=ctensor, degfd=4 )

  call printtofile ( mesh, problemc, filename='c_on_centerline.out', &
    curve=22, vector=ctensor )

  call delete ( ctensor )

! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol )
  call delete ( oldvectors_d, oldvectors_dve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( solc )

  call delete ( coefficients )

end program cylinder_post
