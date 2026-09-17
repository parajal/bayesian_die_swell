
! Post-processing for extrudate_swell5

program post_swell5

  use tfem_m
  use viscoelastic_elements_m
  use figplot_m
  use io_utils_m

  implicit none


! definitions

  type(mesh_t) :: mesh, mesh_plot
  type(input_probdef_t) :: input_probdef, input_probdef_plot, input_probdefc
  type(problem_t) :: problem, problem_plot, problemc
  type(sysvector_t), target :: sol
  type(vector_t) :: velocity, pressure, vorticity
  type(sysvector_t), dimension(3,1), target :: solc
  type(vector_t) :: cxx, cxy, cyy
  type(sample_t) :: sample_p, sample_v, sample_vort, sample_c
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(coefficients_t) :: coefficients

  integer :: obj_plot, ios, i

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


! problem definition of gradient/velocity/pressure

  call read_input_probdef ( mesh, input_probdef, filename='probdef.out' )

  call problem_definition ( input_probdef, mesh, problem )

! problem definition conformation tensor

  call read_input_probdef ( mesh, input_probdefc, filename='probdefc.out' )

  call problem_definition ( input_probdefc, mesh, problemc )

! problem definition of plotting

  call read_input_probdef ( mesh_plot, input_probdef_plot, &
    filename='probdef_plot.out' )

  call problem_definition ( input_probdef_plot, mesh_plot, problem_plot )


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
  read(10) (solc(i,1)%u, i=1,3)

  close ( unit=10 )


! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! store solution vector

  oldvectors%s(1)%p => sol


! create oldvectors for viscoelastic problem

  call create_oldvectors ( oldvectors_ve, nsysvec2=1 )

  oldvectors_ve%s2(1)%p => solc


! post-processing

  call plot_mesh ( plot_options, mesh_plot, 'mesh_xf_plot.fig' )

  call create_vector ( problem_plot, velocity, vec=1 )
  call create_vector ( problem_plot, pressure, vec=2 )
  call create_vector ( problem_plot, vorticity, vec=2 )

  call create_vector ( problem_plot, cxx, vec=2 )
  call create_vector ( problem_plot, cxy, vec=2 )
  call create_vector ( problem_plot, cyy, vec=2 )

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

! sample conformation tensor

  coefficients%i(28) = 1   ! mode number to be sampled

  call fill_sample ( mesh, problemc, sample_c, ndegfd=3, object=obj_plot, &
    elemsub=sample_conformation_tensor, coefficients=coefficients, &
    oldvectors=oldvectors_ve )

  cxx%u = sample_c%u(:,1)
  cxy%u = sample_c%u(:,2)
  cyy%u = sample_c%u(:,3)


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

  call plot_color_fill ( plot_options, mesh_plot, problem_plot, &
    'cxx_xf.fig', vector=cxx )
  call plot_color_fill ( plot_options, mesh_plot, problem_plot, &
    'cxy_xf.fig', vector=cxy )
  call plot_color_fill ( plot_options, mesh_plot, problem_plot, &
    'cyy_xf.fig', vector=cyy )


! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=pressure, &
    dataname='pressure', filename='extrudate_swell5.vtk' )

  call write_vector_vtk ( mesh_plot, problem_plot, &
    filename='extrudate_swell5.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true. )

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=vorticity, &
    filename='extrudate_swell5.vtk', &
    dataname='vorticity', append=.true. )

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=cxx, &
    filename='extrudate_swell5.vtk', dataname='cxx', append=.true. )

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=cxy, &
    filename='extrudate_swell5.vtk', dataname='cxy', append=.true. )

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=cyy, &
    filename='extrudate_swell5.vtk', dataname='cyy', append=.true. )


! delete all data including all allocated memory

  call delete ( mesh, mesh_plot )
  call delete ( problem, problem_plot, problemc )
  call delete ( input_probdef, input_probdefc, input_probdef_plot )
  call delete ( sol )
  call delete ( oldvectors, oldvectors_ve )
  call delete ( coefficients )

end program post_swell5
