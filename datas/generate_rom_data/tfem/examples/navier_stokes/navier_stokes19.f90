! Navier-Stokes problem on a square cavity with 3D velocities (curved channel).
! Unsteady flow.
! Second-order implicit Gear discretization. First step: Euler implicit.

program navier_stokes19

  use tfem_m
  use hsl_ma41_m
  use stokes_elements_m
  use stokes_functions_m
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
    oy = 1._dp,   &       ! left lower corner in r-direction
    flowrate = 1._dp, &   ! flow rate in third direction
    eta = 1._dp           ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, sol_iter, sol_hat
  type(sysvector_t) :: soln, solnm1, rhsd
  type(vector_t) :: velocity, pressure, gammadot, Dtensor, Ltensor
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors, oldvectors_iter
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vel, pres

  integer :: step, iterstep, num_picard, maxnumiterations
  integer :: numtimesteps, numeulertimesteps
  real(dp) :: H, rho, Re
  real(dp) :: diffv, thresh, diffp, deltat, mfac
  real(dp) :: tn, tnp1
  real(dp) :: gamma0=1.5_dp, alpha0=2._dp, alpha1=-0.5_dp

! density
  rho = 100._dp
! number time steps:
  numtimesteps = 200
! number Euler time steps (at least one!):
  numeulertimesteps = 1
! time step
  deltat = 1.e-2_dp

! threshold for the iteration process
  thresh = 1.e-8_dp
! maximum number of iterations allowed to obtain convergence
  maxnumiterations = 100
! number of Picard iterations, followed by Newton iterations
  num_picard = 0

! mesh width and heigth:
  H = 1._dp
! Write out the Reynolds number
  Re=rho * flowrate / H / eta
  write(*,'(a,es12.4)') 'Re = ', Re

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
  coefficients%r(8) = deltat
  coefficients%r(151) = rho

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh
  meshgen_options%elshape = 6
  meshgen_options%ox = 0._dp
  meshgen_options%oy = oy
  meshgen_options%lx = H
  meshgen_options%ly = H
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

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

  call create_sysvector ( problem, sol, sol_iter, soln )
  call create_sysvector ( problem, solnm1, sol_hat )
  call create_sysvector ( problem, rhsd )

! define vector subscripts for direct manipulation of sysvector data

! velocities
  call create_subscript ( mesh, problem, vel, physqarr=[1] )
! pressures
  call create_subscript ( mesh, problem, pres, physqarr=[2] )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )
  call create_oldvectors ( oldvectors_iter, nsysvec=1 )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! set initial value of solution to zero

  soln%u = 0
  sol_iter%u = soln%u

! fill oldvectors

  oldvectors%s(1)%p => sol_hat  ! alpha0*un+alpha1*un-1
  oldvectors_iter%s(1)%p => sol_iter ! solution at end of previous iteration

! start time integration

  call tic

  tn = 0

  integration: do step = 1, numtimesteps

    tnp1 = tn + deltat

    write(*,'(/a,i0,a,es12.4,a/)') ' ** time step = ', step, &
                                   ' time = ', tnp1, ' ** '
!   set flowrate at tn+1

    coefficients%r(6) = flowrate * ( 1 - exp(-10*tnp1) )

!   start iteration for nonlinear term

    iterstep = 0

    iterate: do

      iterstep = iterstep + 1

!     build (assemble) viscous/pressure matrix and vector from elements
      call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
        coefficients=coefficients )

!     first num_picard steps Picard iteration, then Newton iteration
!     the term un+1.grad un+1
      if ( iterstep <= num_picard ) then
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_picard, coefficients=coefficients, &
          oldvectors=oldvectors_iter, physqcol=[physqvel], &
          physqrow=[physqvel], addmatvec=.true. )
      else
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_newton, coefficients=coefficients, &
          oldvectors=oldvectors_iter, physqcol=[physqvel], &
          physqrow=[physqvel], addmatvec=.true. )
      end if

!     instationary term ( rho du/dt )

      if ( step <= numeulertimesteps ) then
        sol_hat%u(vel%s) = soln%u(vel%s)
        mfac = 1
      else
        sol_hat%u(vel%s) = alpha0 * soln%u(vel%s) + alpha1 * solnm1%u(vel%s)
        mfac = gamma0
      end if

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=inertia_elem_dudt, coefficients=coefficients, &
        oldvectors=oldvectors, physqcol=[physqvel], physqrow=[physqvel], &
        factormat=mfac, addmatvec=.true. )

      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        constraint1=1, elemsub=stokes_constr_flowr_surface, &
        addmatvec=.true., coefficients=coefficients )

!     fill solution vector with essential boundary conditions

      call fill_sysvector ( mesh, problem, sol, &
        curve1=1, curve2=4, physq=1, value=0._dp )
      call fill_sysvector ( mesh, problem, sol, &
        point=1, physq=2, value=0._dp )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

      call solve_system_ma41 ( sysmatrix, rhsd, sol )

!     compute the difference between iteration steps
      diffv = maxval ( sqrt ( ( sol%u(vel%s) - sol_iter%u(vel%s) )** 2 ) )
      diffp = maxval ( sqrt ( ( sol%u(pres%s) - sol_iter%u(pres%s) )** 2 ) )

      write(*,'(a,i0,2(a,es12.4))') &
        'iterstep = ', iterstep, ' diffv = ', diffv, ' diffp = ', diffp

      call copy ( sol, sol_iter )

      if ( diffv < thresh ) exit iterate

      if ( iterstep == maxnumiterations ) then
        write(*,'(a,i0)') &
          'Maximum number of iterations reached = ', maxnumiterations
        stop
      end if

    end do iterate

!   copy the converged solution to the old step
    call copy ( soln, solnm1 )
    call copy ( sol, soln )
    tn = tnp1

  end do integration

  call toc

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
    filename='ns19.vtk' )

  call write_vector_vtk ( mesh, problem, filename='ns19.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true., &
    assume3D=.true. )

  call write_vector_vtk ( mesh, problem, filename='ns19.vtk', &
    dataname='velocity_vector_planar', vector=velocity, append=.true., &
    degfd=[1,2] )

  call write_vector_vtk ( mesh, problem, filename='ns19.vtk', &
    dataname='velocity_vector_z', vector=velocity, append=.true., &
    degfd=[0,0,3] )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='ns19.vtk', dataname='gammadot', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='ns19.vtk', &
    dataname='D', vector=Dtensor, append=.true., assume3D=.true. )

  call write_tensor_vtk ( mesh, problem, filename='ns19.vtk', &
    dataname='L', vector=Ltensor, append=.true., assume3D=.true., &
    symmetric=.false. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, sol_iter, soln )
  call delete ( solnm1, sol_hat )
  call delete ( rhsd )
  call delete ( velocity, pressure, gammadot, Dtensor, Ltensor )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( vel, pres )

end program navier_stokes19
