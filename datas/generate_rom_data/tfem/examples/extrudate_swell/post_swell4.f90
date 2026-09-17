
! Post-processing for extrudate_swell4 (planar or axisymmetric Newtonian)

program post_swell4

  use tfem_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m

  implicit none


! definitions

  type(mesh_t) :: mesh, mesh_plot
  type(input_probdef_t) :: input_probdef, input_probdef_plot
  type(problem_t) :: problem, problem_plot
  type(sysvector_t), target :: sol
  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(sample_t) :: sample_p, sample_v, sample_vort

  integer :: obj_plot, ios

! print filename in title and date at footer

  plot_options%printfilename = .true.
  plot_options%printdateandtime = .true.


! read coefficients

  call read_coefficients ( coefficients, filename="coefficients.out" )


! read original mesh

  call read_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )

! read plot mesh

  call read_mesh ( mesh_plot, filename='mesh_plot.out' )

  call fill_mesh_parts ( mesh_plot )


! problem definition of velocity/pressure

  call read_input_probdef ( mesh, input_probdef, filename='probdef.out' )

  call problem_definition ( input_probdef, mesh, problem )

! problem definition of plotting

  call read_input_probdef ( mesh_plot, input_probdef_plot, &
    filename='probdef_plot.out' )

  call problem_definition ( input_probdef_plot, mesh_plot, problem_plot )


! create system vectors for velocity/pressure (solution)

  call create_sysvector ( problem, sol )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! store solution vector

  oldvectors%s(1)%p => sol



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

  call plot_mesh ( plot_options, mesh_plot, 'mesh_xf_plot.fig' )

  call create_vector ( problem_plot, velocity, vec=1 )
  call create_vector ( problem_plot, pressure, vec=2 )
  call create_vector ( problem_plot, vorticity, vec=2 )

! create an object for interpolation on the original mesh

  warn_add_to_mesh_after_meshgen_parts = .false. ! .false. to suppress warning
  call add_to_mesh ( mesh, object='coordinates', coor=mesh_plot%coor )
  obj_plot = mesh%nobjects
  call fill_mesh_parts_objects ( mesh, object1=obj_plot )


! sample pressure

  call fill_sample ( mesh, problem, sample_p, ndegfd=1, object=obj_plot, &
    elemsub=stokes_sample_pressure, coefficients=coefficients, &
    oldvectors=oldvectors )

  pressure%u = sample_p%u(:,1)

! sample velocity

  call fill_sample ( mesh, problem, sample_v, ndegfd=2, object=obj_plot, &
    elemsub=stokes_sample_velocity, coefficients=coefficients, &
    oldvectors=oldvectors )

  velocity%u = reshape ( transpose(sample_v%u), [2*mesh_plot%nnodes] )

! sample vorticity

  coefficients%i(13)=5

  call fill_sample ( mesh, problem, sample_vort, ndegfd=1, object=obj_plot, &
    elemsub=stokes_sample_deriv, coefficients=coefficients, &
    oldvectors=oldvectors )

  vorticity%u = sample_vort%u(:,1)


  plot_options%plotboundary2=.false.

  call plot_color_fill ( plot_options, mesh_plot, problem_plot, &
    'velocity_u_xf.fig', vector=velocity, degfd=1 )

  call plot_color_fill ( plot_options, mesh_plot, problem_plot, &
    'velocity_v_xf.fig', vector=velocity, degfd=2 )

  call plot_color_fill ( plot_options, mesh_plot, problem_plot, &
    'pressure_xf.fig', vector=pressure )

  call plot_color_fill ( plot_options, mesh_plot, problem_plot, &
    'vorticity_xf.fig', vector=vorticity )


  plot_options%scalevector=0.3

  plot_options%fontsize = 10

  call plot_vector ( plot_options, mesh_plot, problem_plot, 'velocity.fig', &
    vector=velocity )


! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=pressure, &
    dataname='pressure', filename='extrudate_swell4.vtk' )

  call write_vector_vtk ( mesh_plot, problem_plot, &
    filename='extrudate_swell4.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true. )

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=vorticity, &
    filename='extrudate_swell4.vtk', &
    dataname='vorticity', append=.true. )


! delete all data including all allocated memory

  call delete ( mesh, mesh_plot )
  call delete ( problem, problem_plot )
  call delete ( input_probdef, input_probdef_plot )
  call delete ( sol )
  call delete ( oldvectors )
  call delete ( coefficients )

end program post_swell4
