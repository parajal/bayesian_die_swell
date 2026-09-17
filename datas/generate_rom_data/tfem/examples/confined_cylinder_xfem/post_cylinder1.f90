
! Post-processing for cylinder1

program post_cylinder1

  use tfem_m
  use viscoelastic_elements_m
  use figplot_m
  use io_utils_m

  implicit none


! constants
  integer, parameter :: &
    nmodes = 1,         & ! number of modes
    ncompc = 3            ! number of conformation tensor components


! definitions

  type(mesh_t) :: mesh, mesh_plot
  type(input_probdef_t) :: input_probdef, input_probdefc, input_probdef_plot
  type(problem_t), target :: problem, problemc
  type(problem_t) :: problem_plot
  type(sysvector_t), target :: sol
  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

  type(sysvector_t), dimension(ncompc,nmodes), target :: solc
  type(vector_t) :: c_xx, c_xy, c_yy
  type(sample_t) :: sample_p, sample_v, sample_vort, sample_c

  integer :: i, ios, obj_plot


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


! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1, nsysvec2=1, nprob=2 )


! store solution vectors and problem structures

  oldvectors%s(1)%p => sol
  oldvectors%s2(1)%p => solc
  oldvectors%p(1)%p => problem
  oldvectors%p(2)%p => problemc


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

  call plot_mesh ( plot_options, mesh_plot, 'mesh_xf_plot.fig' )


! create vectors

  call create_vector ( problem_plot, velocity, vec=1 )
  call create_vector ( problem_plot, pressure, vec=2 )
  call create_vector ( problem_plot, vorticity, vec=2 )

  call create_vector ( problem_plot, c_xx, vec=2 )
  call create_vector ( problem_plot, c_xy, vec=2 )
  call create_vector ( problem_plot, c_yy, vec=2 )


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

! sample conformation

  coefficients%i(28) = 1   ! mode number to be sampled

  call fill_sample ( mesh, problemc, sample_c, ndegfd=3, object=obj_plot, &
    elemsub=sample_conformation_tensor, coefficients=coefficients, &
    oldvectors=oldvectors )

  c_xx%u = sample_c%u(:,1)
  c_xy%u = sample_c%u(:,2)
  c_yy%u = sample_c%u(:,3)


  plot_options%plotboundary2=.false.

  call plot_color_fill ( plot_options, mesh_plot, problem_plot, &
    'velocity_u_xf.fig', vector=velocity, degfd=1 )

  call plot_color_fill ( plot_options, mesh_plot, problem_plot, &
    'velocity_v_xf.fig', vector=velocity, degfd=2 )

  call plot_color_fill ( plot_options, mesh_plot, problem_plot, &
    'pressure_xf.fig', vector=pressure )

  call plot_color_fill ( plot_options, mesh_plot, problem_plot, &
    'vorticity_xf.fig', vector=vorticity )

  call plot_color_fill ( plot_options, mesh_plot, problem_plot, 'c_xx.fig', &
    vector=c_xx )
  call plot_color_fill ( plot_options, mesh_plot, problem_plot, 'c_xy.fig', &
    vector=c_xy )
  call plot_color_fill ( plot_options, mesh_plot, problem_plot, 'c_yy.fig', &
    vector=c_yy )


! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=pressure, &
    dataname='pressure', filename='cylinder1.vtk' )

  call write_vector_vtk ( mesh_plot, problem_plot, &
    filename='cylinder1.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true. )

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=vorticity, &
    filename='cylinder1.vtk', &
    dataname='vorticity', append=.true. )

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=c_xx, &
    filename='cylinder1.vtk', &
    dataname='c_xx', append=.true. )

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=c_xy, &
    filename='cylinder1.vtk', &
    dataname='c_xy', append=.true. )

  call write_scalar_vtk ( mesh_plot, problem_plot, vector=c_yy, &
    filename='cylinder1.vtk', &
    dataname='c_yy', append=.true. )


! delete all data including all allocated memory

  call delete ( mesh, mesh_plot )
  call delete ( problem, problem_plot )
  call delete ( input_probdef, input_probdef_plot )
  call delete ( sol )
  call delete ( oldvectors )
  call delete ( coefficients )

  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( solc )

  call delete ( pressure, velocity, vorticity )
  call delete ( c_xx, c_xy, c_yy )

end program post_cylinder1
