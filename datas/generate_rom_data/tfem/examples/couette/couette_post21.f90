! Viscoelastic problem on a unit square
! linear Couette flow
! DEVSS-G/SUPG
! Post-processing

program couette_post21

  use tfem_m
  use viscoelastic_elements_m
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
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

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


! create system vector

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

  call create_vector ( problem, velocity, physq=coefficients%i(6) )
  call extract_physvector ( mesh, problem, sol, velocity )

  plot_options%scalevector=1e-1

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  call plot_color_fill ( plot_options, mesh, problem, 'velocity_color.fig', &
    vector=velocity, degfd=1 )

  call write_vector_vtk ( mesh, problem, filename='vel.vtk', &
    dataname='velocity', vector=velocity )

! write binary file for reading by streamfunction

  open ( unit=10, form='unformatted', file='velocity.out' )

  write( unit=10 ) velocity%u

  close ( unit=10 )

  call delete ( velocity )


! fill oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol


! pressure

  call create_vector ( problem, pressure, vec=5 )

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )

  call write_scalar_vtk ( mesh, problem, filename='vel.vtk', &
    dataname='pressure', vector=pressure, append=.true. )

  call delete ( pressure )


! vorticity

  call create_vector ( problem, vorticity, vec=5 )

  coefficients%i(13)=5
  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

  call write_scalar_vtk ( mesh, problem, filename='vel.vtk', &
    dataname='vorticity', vector=vorticity, append=.true. )

  call delete ( vorticity )


! conformation tensor

  call create_vector ( problem, cxx, vec=5 )
  call create_vector ( problem, cxy, vec=5 )
  call create_vector ( problem, cyy, vec=5 )

  coefficients%i(13)=1
  call derive_vector ( mesh, problem, cxx, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors )
  coefficients%i(13)=2
  call derive_vector ( mesh, problem, cxy, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors )
  coefficients%i(13)=3
  call derive_vector ( mesh, problem, cyy, elemsub=deriv_conformation, &
    coefficients=coefficients, oldvectors=oldvectors )

  call plot_color_contour ( plot_options, mesh, problem, &
    'cxx_contour.fig', vector=cxx )
  call plot_color_contour ( plot_options, mesh, problem, &
    'cxy_contour.fig', vector=cxy )
  call plot_color_contour ( plot_options, mesh, problem, &
    'cyy_contour.fig', vector=cyy )

  call write_scalar_vtk ( mesh, problem, filename='c.vtk', &
    dataname='cxx', vector=cxx )
  call write_scalar_vtk ( mesh, problem, filename='c.vtk', &
    dataname='cxy', vector=cxy, append=.true. )
  call write_scalar_vtk ( mesh, problem, filename='c.vtk', &
    dataname='cyy', vector=cyy, append=.true. )

! print data on curve to a file

  call printtofile ( mesh, problem, 'cxx_curve1.out', curve=1, vector=cxx )
  call printtofile ( mesh, problem, 'cxx_curve2.out', curve=2, vector=cxx )
  call printtofile ( mesh, problem, 'cxy_curve1.out', curve=1, vector=cxy )
  call printtofile ( mesh, problem, 'cxy_curve2.out', curve=2, vector=cxy )
  call printtofile ( mesh, problem, 'cyy_curve1.out', curve=1, vector=cyy )
  call printtofile ( mesh, problem, 'cyy_curve2.out', curve=2, vector=cyy )

  call delete ( cxx, cxy, cyy )

! print data on cross-section to a file

  call fill_sample ( mesh, problem, sample, ndegfd=3, object=1, &
    elemsub=sample_conformation_tensor, oldvectors=oldvectors, &
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
  call delete ( oldvectors )

  call delete ( coefficients )
  call delete ( sample )

end program couette_post21
