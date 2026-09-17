! Axisymmetrical stick-slip problem

! Post-processing

program post_stick_slip2

  use tfem_m
  use viscoelastic_elements_m
  use figplot_m
  use io_utils_m

  implicit none


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc
  type(problem_t) :: problem, problemc
  type(sysvector_t), target :: sol
  type(vector_t) :: velocity, pressure, vorticity
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors_d, oldvectors_dve
  type(coefficients_t) :: coefficients

  type(sysvector_t), dimension(4,1), target :: solc
  type(vector_t) :: czz, czr, crr, ctt

  integer :: i, ios

  integer, parameter :: wall_curve=12, center_curve=13


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
  read(10) (solc(i,1)%u, i=1,4)

  close ( unit=10 )


! post-processing


! mesh

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! velocity vector

  call create_vector ( problem, velocity, physq=2 )
  call extract_physvector ( mesh, problem, sol, velocity )

  plot_options%scalevector=0.3

  plot_options%fontsize = 10

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  plot_options%rotatecolorbar=.true.
  plot_options%shiftbary=100

  call plot_color_fill ( plot_options, mesh, problem, 'velocity_color.fig', &
    vector=velocity, degfd=1 )

  plot_options%shiftbary = 30

  call printtofile ( mesh, problem, filename='vz_on_wall.out', &
    curve=wall_curve, vector=velocity )

  call printtofile ( mesh, problem, filename='vz_on_centerline.out', &
    curve=center_curve, vector=velocity )

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

  call printtofile ( mesh, problem, filename='p_on_wall.out', &
    curve=wall_curve, vector=pressure )

  call printtofile ( mesh, problem, filename='p_on_centerline.out', &
    curve=center_curve, vector=pressure )

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

  call create_vector ( problemc, czz, vec=2 )
  call create_vector ( problemc, czr, vec=2 )
  call create_vector ( problemc, crr, vec=2 )
  call create_vector ( problemc, ctt, vec=2 )

  coefficients%i(13)=1
  call derive_vector ( mesh, problemc, czz, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )
  coefficients%i(13)=2
  call derive_vector ( mesh, problemc, czr, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )
  coefficients%i(13)=3
  call derive_vector ( mesh, problemc, crr, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )
  coefficients%i(13)=4
  call derive_vector ( mesh, problemc, ctt, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )

  call plot_color_contour ( plot_options, mesh, problemc, &
    'czz_contour.fig', vector=czz )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'czr_contour.fig', vector=czr )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'crr_contour.fig', vector=crr )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'ctt_contour.fig', vector=ctt )

  call printtofile ( mesh, problemc, filename='czz_on_wall.out', &
    curve=wall_curve, vector=czz )

  call printtofile ( mesh, problemc, filename='czr_on_wall.out', &
    curve=wall_curve, vector=czr )

  call printtofile ( mesh, problemc, filename='crr_on_wall.out', &
    curve=wall_curve, vector=crr )

  call printtofile ( mesh, problemc, filename='ctt_on_wall.out', &
    curve=wall_curve, vector=ctt )

  call delete ( czz, czr, crr, ctt )


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

end program post_stick_slip2
