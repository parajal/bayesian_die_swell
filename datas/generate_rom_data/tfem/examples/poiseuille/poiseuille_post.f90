! Viscoelastic problem on a unit square
! planar Poiseuille flow
! Post-processing

program poiseuille_post

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

  type(sysvector_t), dimension(:,:), allocatable, target :: solc
  type(vector_t) :: cxx, cxy, cyy, viscous_stress, viscoelastic_stress, &
    tau_tensor, workvec1, workvec2, vonmises, trace, regime
  type(subscriptvec_t) :: svec

  logical :: viscous, assume33
  integer :: i, j, ios, nmodes, ncomp
  real(dp):: x, y, p
  real(dp), dimension(:), allocatable :: tau, sigma


! print filename in title and date at footer

  plot_options%printfilename = .true.
  plot_options%printdateandtime = .true.


! read coefficients

  call read_coefficients ( coefficients, filename="coefficients.out" )

  nmodes = coefficients%i(19)

  if ( coefficients%i(71) == 0 .and. coefficients%i(18) /= 22 ) then
    ncomp = 3
  else
    ncomp = 4
  end if

  allocate ( solc(ncomp,nmodes) )


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
  read(10) ((solc(i,j)%u, i=1,ncomp), j=1,nmodes)

  close ( unit=10 )


! post-processing


! mesh

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! velocity vector

  call create_vector ( problem, velocity, physq=2 )
  call extract_physvector ( mesh, problem, sol, velocity )

  plot_options%scalevector=1e-1

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  call plot_color_fill ( plot_options, mesh, problem, 'velocity_color.fig', &
    vector=velocity, degfd=1 )

  call write_vector_vtk ( mesh, problem, filename='poiseuille.vtk', &
    dataname='velocity', vector=velocity )

! write binary file for reading by streamfunction

  open ( unit=10, form='unformatted', file='velocity.out' )

  write( unit=10 ) velocity%u

  close ( unit=10 )

  call printtofile ( mesh, problem, 'velocity_curve2.out', curve=2, &
    vector=velocity )

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

  call printtofile ( mesh, problem, 'pressure_curve2.out', curve=2, &
    vector=pressure )

  call write_scalar_vtk ( mesh, problem, filename='poiseuille.vtk', &
    dataname='pressure', vector=pressure, append=.true. )

  call delete ( pressure )


! vorticity

  call create_vector ( problem, vorticity, vec=4 )

  coefficients%i(13)=5
  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors_d )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

  call printtofile ( mesh, problem, 'vorticity_curve2.out', curve=2, &
    vector=vorticity )

  call write_scalar_vtk ( mesh, problem, filename='poiseuille.vtk', &
    dataname='vorticity', vector=vorticity, append=.true. )

  call delete ( vorticity )


  viscous = coefficients%r(1) >= tiny(1._dp)

  if ( viscous ) then

!   viscous stress

    call create_vector ( problem, viscous_stress, vec=5 )

    call derive_vector ( mesh, problem, viscous_stress, &
      elemsub=stokes_stress_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_d )

    call plot_color_contour ( plot_options, mesh, problem, &
      'taus_xy.fig', vector=viscous_stress, degfd=2 )

    call write_tensor_vtk ( mesh, problem, filename='poiseuille.vtk', &
      dataname='viscous_stress', vector=viscous_stress, append=.true. )

  end if


! create oldvectors for viscoelastic problem

  call create_oldvectors ( oldvectors_dve, nsysvec2=1 )

  oldvectors_dve%s2(1)%p => solc


! viscoelastic stress

  call create_vector ( problemc, viscoelastic_stress, vec=3 )
  call create_vector ( problemc, tau_tensor, vec=3 )
  call create_vector ( problemc, workvec1, vec=2 )
  call create_vector ( problemc, workvec2, vec=2 )

  call derive_vector ( mesh, problemc, viscoelastic_stress, &
    elemsub=deriv_viscoelastic_stress_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_dve )

  call plot_color_contour ( plot_options, mesh, problemc, &
    'tauve_xy.fig', vector=viscoelastic_stress, degfd=2 )

  assume33 = problemc%vec_elnumdegfd(1)%a(1,3) == 4

  call write_tensor_vtk ( mesh, problemc, filename='poiseuille.vtk', &
    dataname='viscoelastic_stress', vector=viscoelastic_stress, &
    append=.true., assume33=assume33 )

! tauve_xx
  call extract_vector ( mesh, problemc, viscoelastic_stress, workvec1, &
    indegfd=[1] )
! tauve_yy
  call extract_vector ( mesh, problemc, viscoelastic_stress, workvec2, &
    indegfd=[3] )

  workvec1%u = workvec1%u - workvec2%u  ! N1

  call plot_color_contour ( plot_options, mesh, problemc, &
    'N1.fig', vector=workvec1 )

  call write_scalar_vtk ( mesh, problemc, filename='poiseuille.vtk', &
    dataname='N1', vector=workvec1, append=.true. )

  tau_tensor%u = viscoelastic_stress%u

  if ( viscous ) then
    call create ( mesh, problemc, svec, degsfd=[1,2,3], vec=3 )
    tau_tensor%u(svec%s) = tau_tensor%u(svec%s) + viscous_stress%u
  end if

  call plot_color_contour ( plot_options, mesh, problemc, &
    'tau_xy.fig', vector=tau_tensor, degfd=2 )

  call printtofile ( mesh, problemc, 'tau_curve2.out', curve=2, &
    vector=tau_tensor )

  allocate ( tau(problemc%vec_elnumdegfd(1)%a(1,3)) )

  open ( unit=10, file='pressure_curve2.out', status='old' )
  open ( unit=11, file='tau_curve2.out', status='old' )
  open ( unit=12, file='sigma_curve2.out', status='replace', recl=300 )

  do

    read (10,*, end=10 ) x, y, p
    read (11,*) x, y, tau
    sigma = tau
    if ( size(tau) == 4 ) then
      sigma([1,3,4]) = tau([1,3,4]) - p
    else
      sigma([1,3]) = tau([1,3]) - p
    end if
    write (12,*) x, y, sigma

  end do

  10 close ( unit=10 ); close ( unit=11 ); close ( unit=12 )

  call delete ( viscoelastic_stress, workvec1, workvec2 )

  call write_tensor_vtk ( mesh, problemc, filename='poiseuille.vtk', &
    dataname='tau_tensor', vector=tau_tensor, &
    append=.true., assume33=assume33 )

  if ( viscous ) call delete ( viscous_stress )


! conformation tensor (mode 1 only)

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

  if ( coefficients%i(80) >= 2 .or. coefficients%i(18) == 22 ) then

!   evp model: von mises shear stress

    call create_vector ( problemc, vonmises, vec=2 )

    coefficients%i(13)=1 ! component
    coefficients%i(28)=1 ! mode

    call derive_vector ( mesh, problemc, vonmises, &
      elemsub=deriv_viscoelastic_stress_scalar, &
      coefficients=coefficients, oldvectors=oldvectors_dve )

    call plot_color_contour ( plot_options, mesh, problemc, &
      'vonmises.fig', vector=vonmises )

    call printtofile ( mesh, problemc, 'vonmises_curve2.out', curve=2, &
      vector=vonmises )

    call write_scalar_vtk ( mesh, problemc, filename='poiseuille.vtk', &
      dataname='von_Mises', vector=vonmises, append=.true. )

  end if

  if ( coefficients%i(18) == 22 ) then

!   Saramito DP model: trace and regime

    call create_vector ( problemc, trace, vec=2 )

    coefficients%i(13)=2 ! component
    coefficients%i(28)=1 ! mode

    call derive_vector ( mesh, problemc, trace, &
      elemsub=deriv_viscoelastic_stress_scalar, &
      coefficients=coefficients, oldvectors=oldvectors_dve )

    call plot_color_contour ( plot_options, mesh, problemc, &
      'trace.fig', vector=trace )

    call printtofile ( mesh, problemc, 'trace_curve2.out', curve=2, &
      vector=trace )

    call write_scalar_vtk ( mesh, problemc, filename='poiseuille.vtk', &
      dataname='trace_tau', vector=trace, append=.true. )

    call create_vector ( problemc, regime, vec=2 )

    coefficients%i(13)=3 ! component
    coefficients%i(28)=1 ! mode

    call derive_vector ( mesh, problemc, regime, &
      elemsub=deriv_viscoelastic_stress_scalar, &
      coefficients=coefficients, oldvectors=oldvectors_dve )

    call plot_color_contour ( plot_options, mesh, problemc, &
      'regime.fig', vector=regime )

    call printtofile ( mesh, problemc, 'regime_curve2.out', curve=2, &
      vector=regime )

    call write_scalar_vtk ( mesh, problemc, filename='poiseuille.vtk', &
      dataname='regime', vector=regime, append=.true. )

  end if

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

end program poiseuille_post
