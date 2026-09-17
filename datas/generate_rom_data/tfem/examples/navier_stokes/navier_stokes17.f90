! Navier-Stokes problem on a square cavity with 3D velocities (curved channel).
! Steady flow, gradual increase of the Reynolds number.
! Imposed flow rate in theta-direction. Axisymmetric.

program navier_stokes17

  use tfem_m
  use hsl_ma41_m
  use stokes_elements_m
  use inertia_elements_m
  use io_utils_m
  use figplot_m
  use timer_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=40,              & ! number of elements in z
    ny=40                 ! number of elements in r

  real(dp), parameter :: &
    oy = 1._dp,   &  ! left lower corner in r-direction
    H   = 1._dp,  &  ! height of the cavity
    eta = 1._dp,  &  ! viscosity
    flowrate = 1._dp ! flow rate in third direction

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, soln
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, gammadot, Dtensor, Ltensor
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vel, pres

  integer :: step, nsteps, iterstep, nitermax, num_picard, maxnumiterations
  real(dp) :: rho, difference, thresh, diffp

! A stepping method for the Reynolds number is used, with an initial solution
! of the previously computed Reynolds number, starting with the Stokes solution
! Re=0
! starting and stepping density
! rho_{i+1} = rho_{i} + rho
  rho = 20._dp
! number of incremental steps:
  nsteps = 5
! final Reynolds number is: rho*nsteps*U*H/eta

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

  coefficients%i(23) = 1  ! coorsys = 1, cilindrical coordinates, axisymmetric
  coefficients%i(67) = 1  ! 3D velocity

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(6) = flowrate
  coefficients%r(151) = rho

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%oy = oy
  meshgen_options%lx = H
  meshgen_options%ly = H
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

! make domain into a surface for the flow rate constraint
  call add_to_mesh ( mesh, surfacefromgroups=[1] )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [3,3,3,3,3,3,3,3,3,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1,    &  ! scalar
                 6,6,6,6,6,6,6,6,6,    &  ! symmetric tensor
                 9,9,9,9,9,9,9,9,9 ], &  ! non-symmetric tensor
                 [9,5] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=1, surface1=1, nglobalc=1 )

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
    curve1=1, curve2=4, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr_surface, &
    addmatvec=.true., coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma41 ( sysmatrix, rhsd, sol )
  call copy ( sol, soln )

  oldvectors%s(1)%p => soln

  nitermax = 0

  call tic

! start increments

  increm: do step = 1, nsteps

!   set the density

    coefficients%r(151) = step * rho

    write(*,*) 'Re = ', coefficients%r(151)*flowrate/H/eta

    iterstep = 0

!   start iteration within increments

    iterate: do

      iterstep = iterstep + 1

!     build (assemble) matrix and vector from elements
      call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
        coefficients=coefficients )

!     first num_picard steps Picard iteration, then Newton iteration
      if ( iterstep <= num_picard ) then
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_picard, coefficients=coefficients, &
          oldvectors=oldvectors, physqcol=[physqvel],physqrow=[physqvel], &
          addmatvec=.true. )
      else
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_newton, coefficients=coefficients, &
          oldvectors=oldvectors, physqcol=[physqvel],physqrow=[physqvel], &
          addmatvec=.true. )
      end if

      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        constraint1=1, elemsub=stokes_constr_flowr_surface, &
        addmatvec=.true., coefficients=coefficients )

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

  call toc

  write(*,*) 'Max number of iterations: ', nitermax

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, gammadot, vec=3 )
  call create_vector ( problem, Dtensor, vec=4 )
  call create_vector ( problem, Ltensor, vec=5 )

  call extract_physvector ( mesh, problem, sol, velocity )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=11

  call derive_vector ( mesh, problem, gammadot, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Dtensor, elemsub=stokes_D_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Ltensor, elemsub=stokes_gradu_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='ns17.vtk' )

  call write_vector_vtk ( mesh, problem, filename='ns17.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true., &
    assume3D=.true. )

  call write_vector_vtk ( mesh, problem, filename='ns17.vtk', &
    dataname='velocity_vector_planar', vector=velocity, append=.true., &
    degfd=[1,2] )

  call write_vector_vtk ( mesh, problem, filename='ns17.vtk', &
    dataname='velocity_vector_z', vector=velocity, append=.true., &
    degfd=[0,0,3] )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='ns17.vtk', dataname='gammadot', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='ns17.vtk', &
    dataname='D', vector=Dtensor, append=.true., assume3D=.true. )

  call write_tensor_vtk ( mesh, problem, filename='ns17.vtk', &
    dataname='L', vector=Ltensor, append=.true., assume3D=.true., &
    symmetric=.false. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, soln, rhsd )
  call delete ( velocity, pressure, gammadot, Dtensor, Ltensor )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( vel, pres )

end program navier_stokes17
