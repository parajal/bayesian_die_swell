! Viscoelastic problem on a unit square
! linear Couette flow
! DEVSS-G/SUPG
! Post-processing

program couette_post

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

  type(sysvector_t), dimension(3,1), target :: solc
  type(vector_t) :: cxx, cxy, cyy
  type(sample_t) :: sample

  integer :: i, ios
  integer, parameter :: nc = 101
  real(dp) :: xc(nc,2)


! print filename in title and date at footer

  plot_options%printfilename = .true.
  plot_options%printdateandtime = .true.


! read coefficients

  call read_coefficients ( coefficients, filename="coefficients.out" )


! read mesh

  call read_mesh ( mesh, filename='mesh.out' )

! generate a cross section at x=0.9 using an object with nc points

  xc(:,1) = 0.9_dp
  xc(:,2) = [ (i*1._dp/(nc-1), i=0,nc-1) ]

  call add_to_mesh ( mesh, object='coordinates', coor=xc )

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
  read(10) (solc(i,1)%u, i=1,3)

  close ( unit=10 )


! post-processing


! mesh

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! velocity vector

  call create_vector ( problem, velocity, physq=coefficients%i(6) )
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

  call create_oldvectors ( oldvectors_dve, nsysvec2=1 )

  oldvectors_dve%s2(1)%p => solc


! conformation tensor

  call create_vector ( problemc, cxx, vec=2 )
  call create_vector ( problemc, cxy, vec=2 )
  call create_vector ( problemc, cyy, vec=2 )

  coefficients%i(13)=1
  call derive_vector ( mesh, problemc, cxx, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )
  coefficients%i(13)=2
  call derive_vector ( mesh, problemc, cxy, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )
  coefficients%i(13)=3
  call derive_vector ( mesh, problemc, cyy, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors_dve )

  call plot_color_contour ( plot_options, mesh, problemc, &
    'cxx_contour.fig', vector=cxx )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'cxy_contour.fig', vector=cxy )
  call plot_color_contour ( plot_options, mesh, problemc, &
    'cyy_contour.fig', vector=cyy )

! print data on curve to a file

  call printtofile ( mesh, problemc, 'cxx_curve1.out', curve=1, vector=cxx )
  call printtofile ( mesh, problemc, 'cxx_curve2.out', curve=2, vector=cxx )
  call printtofile ( mesh, problemc, 'cxy_curve1.out', curve=1, vector=cxy )
  call printtofile ( mesh, problemc, 'cxy_curve2.out', curve=2, vector=cxy )
  call printtofile ( mesh, problemc, 'cyy_curve1.out', curve=1, vector=cyy )
  call printtofile ( mesh, problemc, 'cyy_curve2.out', curve=2, vector=cyy )

  call delete ( cxx, cxy, cyy )

! print data on cross-section to a file

  call fill_sample ( mesh, problemc, sample, ndegfd=3, object=1, &
    elemsub=sample_conformation_tensor, oldvectors=oldvectors_dve, &
    coefficients=coefficients )

  open( unit=10, file='c_on_x=0.9.out', recl=300 )
  do i = 1, nc
    write( unit=10, fmt=* ) sample%coor(i,1), sample%coor(i,2), &
                            sample%u(i,1:3)
  end do
  close( unit=10 )

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
  call delete ( sample )

end program couette_post
