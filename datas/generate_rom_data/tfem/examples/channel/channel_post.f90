! 3D viscoelastic problem in a channel with a square cross section.
! Computed is the developed flow using only one layer of elements in the flow
! direction (x) and assuming periodical boundary conditions.
! DEVSS-G/SUPG
! Post-processing

program channel_post

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
  type(vector_t) :: velocity, pressure
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors_d, oldvectors_dve
  type(coefficients_t) :: coefficients

  type(sysvector_t), dimension(6,1), target :: solc
  type(vector_t) :: cxx, cxy, cxz, cyy, cyz, czz

  integer :: i, ios
  integer, parameter :: nc = 101
  real(dp) :: xc(nc,3)


! print filename in title and date at footer

  plot_options%printfilename = .true.
  plot_options%printdateandtime = .true.


! read coefficients

  call read_coefficients ( coefficients, filename="coefficients.out" )


! read mesh

  call read_mesh ( mesh, filename='mesh.out' )

! generate a cross section at x=0.5, y=0.5 using an object with nc points

  xc(:,1) = 0.5_dp
  xc(:,2) = 0.5_dp
  xc(:,3) = [ (i*1._dp/(nc-1), i=0,nc-1) ]

  call add_to_mesh ( mesh, object='coordinates', coor=xc )

  call fill_mesh_parts ( mesh )


! problem definition of gradient/velocity/pressure

  call read_input_probdef ( mesh, input_probdef, filename='probdef.out' )

  call problem_definition ( input_probdef, mesh, problem )


! problem definition conformation tensor

  call read_input_probdef ( mesh, input_probdefc, filename='probdefc.out' )

  call problem_definition ( input_probdefc, mesh, problemc )


! create system vectors for gradient/velocity/pressure (solution)

  call create ( problem, sol )

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
  read(10) (solc(i,1)%u, i=1,6)

  close ( unit=10 )


! post-processing


! mesh


! plot curves, surfaces and mesh

  plot_options%viewpoint=[1.,0.8,0.4]
  plot_options%fontsize=10
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig', surfaces=[3,4,6] )

  plot_options%printlabels=.false.

! velocity vector

  call create_vector ( problem, velocity, physq=2 )
  call extract_physvector ( mesh, problem, sol, velocity )

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity, surfaces=[3] )
  call plot_points_curves ( plot_options, mesh, 'velocity.fig', append=.true. )

  call plot_color_fill ( plot_options, mesh, problem, 'vx_color.fig', &
    vector=velocity, degfd=1, surfaces=[3] )
  call plot_points_curves ( plot_options, mesh, 'vx_color.fig',  &
    append=.true. )

  call plot_color_fill ( plot_options, mesh, problem, 'vy_color.fig', &
    vector=velocity, degfd=2, surfaces=[3] )
  call plot_points_curves ( plot_options, mesh, 'vy_color.fig',  &
   append=.true. )

  plot_options%scalevector=100

  call plot_vector ( plot_options, mesh, problem, 'crossv.fig', &
    vector=velocity, degfd=[0,2,3], surfaces=[3] )
  call plot_points_curves ( plot_options, mesh, 'crossv.fig', append=.true. )


  call delete ( velocity )


! fill oldvectors for stokes problem

  call create_oldvectors ( oldvectors_d, nsysvec=1 )

  oldvectors_d%s(1)%p => sol


! pressure

  call create_vector ( problem, pressure, vec=4 )

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors_d )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure, surfaces=[3,4,6] )
  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour_S3.fig', vector=pressure, surfaces=[3] )
  call plot_points_curves ( plot_options, mesh, 'pressure_contour_S3.fig', &
     append=.true. )

  call delete ( pressure )


! create oldvectors for viscoelastic problem

  call create_oldvectors ( oldvectors_dve, nsysvec2=1 )

  oldvectors_dve%s2(1)%p => solc


! conformation tensor

  call create_vector ( problemc, cxx, vec=2 )
  call create_vector ( problemc, cxy, vec=2 )
  call create_vector ( problemc, cxz, vec=2 )
  call create_vector ( problemc, cyy, vec=2 )
  call create_vector ( problemc, cyz, vec=2 )
  call create_vector ( problemc, czz, vec=2 )

  coefficients%i(13)=1
  call derive_vector ( mesh, problemc, cxx, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )
  coefficients%i(13)=2
  call derive_vector ( mesh, problemc, cxy, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )
  coefficients%i(13)=3
  call derive_vector ( mesh, problemc, cxz, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )
  coefficients%i(13)=4
  call derive_vector ( mesh, problemc, cyy, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )
  coefficients%i(13)=5
  call derive_vector ( mesh, problemc, cyz, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )
  coefficients%i(13)=6
  call derive_vector ( mesh, problemc, czz, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )

  call plot_color_contour ( plot_options, mesh, problemc, &
    'cxx_contour.fig', vector=cxx, surfaces=[3] )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'cxy_contour.fig', vector=cxy, surfaces=[3] )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'cxz_contour.fig', vector=cxz, surfaces=[3] )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'cyy_contour.fig', vector=cyy, surfaces=[3] )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'cyz_contour.fig', vector=cyz, surfaces=[3] )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'czz_contour.fig', vector=czz, surfaces=[3] )

! print data on curve to a file

  call printtofile ( mesh, problemc, 'cxx_curve2.out', curve=2, vector=cxx )

  call delete ( cxx, cxy, cxz )
  call delete ( cyy, cyz, czz )

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

end program channel_post
