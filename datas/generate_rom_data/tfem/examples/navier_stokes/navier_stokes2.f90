! Navier-Stokes problem on a 3D cavity. Two moving boundaries,
! and Re=500. Either Picard or Newton or a combination of both can be used for
! the iteration process

program navier_stokes2

  use tfem_m
  use hsl_ma41_m
  use stokes_elements_m
  use inertia_elements_m
  use figplot_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    nx = 10,            & ! number of elements in x-direction
    ny = 5,             & ! number of elements in y-direction
    nz = 10,            & ! number of elements in z-direction
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3             ! 3x3x3 integration of hexahedra

  real(dp), parameter :: &
    U  = 1._dp,          & ! velocity of both the lids
    lx = 1._dp,          & ! size in x-direction
    ly = 0.5_dp,         & ! size in y-direction
    lz = 1._dp,          & ! size in z-direction
    eta = 1._dp            ! viscosity

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, soln
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, dudz
  type(plot_options_t) :: plot_options
  type(meshgen_options_t) :: meshgen_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vel, pres

! variables

  integer :: presnod(8) = [1,3,9,7,19,21,27,25]

  integer :: step, nsteps, iterstep, nitermax, num_picard, maxnumiterations

  real(dp) :: rho, difference, thresh, diffp

! A stepping method for the Reynolds number is used, with an initial solution
! of the previously computed Reynolds number, starting with the Stokes solution
! Re=0
! starting and stepping density
! rho_{i+1} = rho_{i} + rho
  rho = 100._dp
! number of steps:
  nsteps = 5
! final Reynolds number is: rho*U*lz*nsteps/eta, since UH=1

! threshold for the iteration process
  thresh = 1.e-5_dp
! maximum number of iterations allowed to obtain convergence
  maxnumiterations = 100
! number of Picard iterations, followed by Newton iterations
  num_picard = 0

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=250, ncoefr=200 )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(151) = rho

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, &
    lz=lz, elshape=14, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! plot curves, surfaces and mesh

  plot_options%viewpoint=[1.,0.8,0.4]
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig', surfaces=[3,4,6] )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(presnod,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, surface1=1, surface2=3, physq=1 )
  call define_essential ( mesh, input_probdef, surface1=4, physq=1, &
    degfd=[0,1,0] )
  call define_essential ( mesh, input_probdef, surface1=5, surface2=6, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, soln )
  call create_sysvector ( problem, rhsd )

! define vector subscripts for direct manipulation of sysvector data

! velocities
  call create_subscript ( mesh, problem, vel, physqarr=[1] )
! pressures
  call create_subscript ( mesh, problem, pres, physqarr=[2] )


! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    surface1=1, surface2=3, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=4, physq=1, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=5, surface2=6, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=6, physq=1, degfd=1, value=U )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=2, physq=1, degfd=3, value=U )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=2 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements
  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma41 ( sysmatrix, rhsd, sol )

  call copy ( sol, soln )

  oldvectors%s(1)%p => soln

  nitermax = 0

! start increments

  increm: do step = 1, nsteps

!   set the density

    coefficients%r(151) = step * rho

    write(*,*) 'Re = ', coefficients%r(151)*U*lz/eta

    iterstep = 0

!   start iteration within increments

    iterate: do

      iterstep = iterstep + 1

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=stokes_elem, coefficients=coefficients )

!     first num_picard steps Picard iteration, then Newton iteration
      if ( iterstep <= num_picard ) then
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_picard, coefficients=coefficients, &
          oldvectors=oldvectors, physqcol=[physqvel],physqrow=[physqvel],&
          addmatvec=.true. )
      else
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_newton, coefficients=coefficients, &
          oldvectors=oldvectors, physqcol=[physqvel],physqrow=[physqvel],&
          addmatvec=.true. )
      end if

      call check ( sysmatrix )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

      call solve_system_ma41 ( sysmatrix, rhsd, sol )

!     compute the difference between iteration steps
      difference = maxval ( sqrt ( ( sol%u(vel%s) - soln%u(vel%s) )** 2 ) )
      diffp = maxval ( sqrt ( ( sol%u(pres%s) - soln%u(pres%s) )** 2 ) )

      write(*,'(a,i0,2(a,es12.4))') &
        'iterstep = ', iterstep, ' diffv = ', difference, ' diffp = ', diffp

      call copy ( sol, soln )

      if ( difference < thresh ) exit iterate

      if ( iterstep == maxnumiterations ) then
        write(*,'(a,i0)') &
          'Maximum number of iterations reached = ', maxnumiterations
        stop
      end if

    end do iterate

    nitermax = max ( iterstep, nitermax )

  end do increm

  write(*,*) 'Max number of iterations: ', nitermax

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, dudz, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=3
  call derive_vector ( mesh, problem, dudz, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )


! write to a fig file for plotting

  plot_options%printlabels=.false.

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity, surfaces=[4,6] )
  call plot_points_curves ( plot_options, mesh, 'velocity.fig', append=.true. )

  call plot_points_curves ( plot_options, mesh, 'velocity_color.fig',  &
   append=.true. )

  call plot_color_contour ( plot_options, mesh, problem, 'dudz_contour.fig', &
    vector=dudz, surfaces=[3,4,6] )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure, surfaces=[3,4,6] )

! write vtk data
  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='cavity_inertia_3D.vtk' )

  call write_vector_vtk ( mesh, problem, filename='cavity_inertia_3D.vtk', &
    dataname='velocity_vector', sysvector=sol, append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, soln )
  call delete ( sysmatrix )
  call delete ( velocity, pressure, dudz )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( vel, pres )

end program navier_stokes2
