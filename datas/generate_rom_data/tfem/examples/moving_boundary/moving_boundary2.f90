! Moving boundary due to inertia and gravity forces.
! 2D domain, axisymmetric with 3D velocities.
! Wall(s) and bottom rotate with an angular velocities omega.
! Perfect slip (in plane) on the walls.
! Surface tension on fluid-air interface and fluid-wall interface.
! Newtonian fluid.

program moving_boundary2

  use tfem_m
  use stokes_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m
  use inertia_elements_m
  use surface_advection_elements_m
  use update_mesh_nodes1_m
  use figplot_m

  implicit none

! constants

! T: Single wall, symmetry line (planar) or center line (axisymmetric)
!    Note, that for the axisymmetric case R1=0 in the mesh generation!
! F: Two walls
  logical, parameter :: symmetric = .true.

! constants flow problem

  integer, parameter :: &
    uintpl = 8,  & ! Q2 velocities
    pintpl = 4,  & ! Q1 pressures
    physqv = 1,  & ! physical quantity nr of the velocities
    physqp = 2,  & ! physical quantity nr of the pressures
    gauss = 3      ! 3x3 integration of quads

! constants surface advection

  integer, parameter :: &
    hintpl = 6,       & ! P2 height function
    ninti_sf_adv = 3, & ! number of Gauss points
    method = 1          ! discretization method 0: Galerkin, 1: SUPG


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix, sysmatrix_stored
  type(sysvector_t), target :: soln, sol_iter, sol_hat
  type(sysvector_t) :: sol, solnm1, rhsd, rhsd_stored
  type(coefficients_t) :: coefficients
  type(vector_t) :: velocity, pressure, gammadot
  type(solver_options_ma41_t) :: solver_options_u
  type(oldvectors_t) :: oldvectors, oldvectors_hat, oldvectors_iter
  type(subscript_t) :: vel8, vel, pres

! 1D height function for surface advection
  type(mesh_t) :: mesh_sf_adv
  type(input_probdef_t) :: input_probdef_sf_adv
  type(problem_t) :: problem_sf_adv
  type(sysmatrix_t) :: sysmatrix_sf_adv
  type(sysvector_t), target :: sol_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1
  type(sysvector_t) :: rhsd_sf_adv, sol_sf_adv_pred
  type(oldvectors_t) :: oldvectors_sf_adv
  type(coefficients_t) :: coefficients_sf_adv
  type(vector_t), target :: velocity_sf_adv
  type(subscript_t) :: hgt, hgt_end

  type(vector_t), target :: meshvel
  real(dp), allocatable, dimension(:,:) :: meshcoor_initial

! ALE mesh motion problem

  type(problem_t), target ::  problem_update_mesh
  real(dp), allocatable, dimension(:) :: disp

! variables

  integer :: &
    coorsys = 1, & ! planar Cartesian (0) or axisymmetric (1) coor. system
    vfuncnr = 1, & ! function nr for gravity 0=none 1=rho g (-1,0)
    numtimesteps = 250,   & ! number of time steps
    numeulertimesteps = 1, & ! number Euler time steps (at least one!):
    maxnumiterations = 100, & ! maximum number of iterations allowed to obtain
                              ! convergence
    num_picard = 1, &         ! number of Picard iterations, followed by Newton
                              ! iterations
    plot_mesh_every = 10    ! plot mesh every .. steps

  real(dp) :: &
    rho = 1.0_dp,      & ! density
    grv = 1.0_dp,      & ! gravity constant g
    eta = 1.0_dp,      & ! viscosity
    gamma_fa = 0.1_dp, & ! surface tension coefficient fluid-air interface
    dgamma_w = 0.01_dp,& ! surface tension difference fluid-air wrt the wall
    omega = 1.0_dp,    & ! angular rate of the wall(s) and bottom
    deltat = 5.e-2_dp, & ! time step
    thresh = 1.e-8_dp, & ! threshold for the iteration process
    betah = 0.5_dp,    & ! beta parameter for SUPG of the height function
    rs_up = 1.5_dp,    & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.0_dp       ! integer_storage velocity-pressure LU (HSL)

  integer :: vertices(4) = [1,3,5,7]
  integer :: step, nnodes, curvew, iterstep
  real(dp) :: diffv, diffp, mfac
  real(dp) :: tn, tnp1
  real(dp) :: gamma0=1.5_dp, alpha0=2._dp, alpha1=-0.5_dp

  real(dp), allocatable, dimension(:,:) :: meshcoor_n, meshcoor_nm1

  character(len=20) :: filename

  logical :: surface_tension = .true.  ! include surface tension


! namelist for input of variables; read from standard input

  namelist /comppar/ numtimesteps, plot_mesh_every, &
    rho, grv, eta, gamma_fa, dgamma_w, omega, deltat, rs_up, is_up, &
    surface_tension

  read ( unit=*, nml=comppar )


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=600, ncoefr=300 )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl, pintpl, 0, 0,     0, &
      physqv, physqp, 0, 0, gauss, &
      gauss ]
  coefficients%i(14) = vfuncnr
  coefficients%i(23) = coorsys
  coefficients%i(48) = 1 ! use mesh velocity for ALE formulation
  coefficients%i(67) = 1 ! 3D velocity

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(8) = deltat
  coefficients%r(151) = rho

  coefficients%vfunc => vfunc

! create mesh

  call read_mesh ( mesh, filename='mesh2b.out' )

  call fill_mesh_parts ( mesh )

  if ( coorsys==1 ) then
    if ( symmetric .and. abs(mesh%coor(mesh%points(1),2)) > tiny(0._dp) ) then
      write(*,'(a/)') &
        'Error: r coordinates of the center line non-zero. Change mesh.'
      stop
    end if
    if ( mesh%coor(mesh%points(1),2) < 0._dp ) then
      write(*,'(a/)') 'Error: negative r coordinates. Change mesh.'
      stop
    end if
  end if

! problem definition of velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(vertices,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,4) = 2  ! vector for mesh velocity

  input_probdef%physq = [1,2]
  input_probdef%probnr = 1

! Dirichlet boundary conditions

! bottom
  call define_essential ( mesh, input_probdef, &
    curve1=9, physq=physqv )
! walls (perfect slip in z-direction)
  call define_essential ( mesh, input_probdef, &
    curve1=10, physq=physqv, degfd=[0,1,1] )

  call problem_definition ( input_probdef, mesh, problem )


! create system vectors for velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, sol_iter, soln )
  call create_sysvector ( problem, solnm1, sol_hat )
  call create_sysvector ( problem, rhsd )

! define vector subscripts for direct manipulation of sysvector data

! velocities
  call create_subscript ( mesh, problem, vel, physqarr=[1] )
! pressures
  call create_subscript ( mesh, problem, pres, physqarr=[2] )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )
  call create_sysmatrix_data ( sysmatrix )
  call create_sysmatrix ( sysmatrix, sysmatrix_stored )

! set initial value of solution to zero

  soln%u = 0
  sol_iter%u = soln%u

! define some arrays and vector for the ALE mesh movement

  allocate ( meshcoor_initial(mesh%nnodes,2) )
  meshcoor_initial = mesh%coor

  allocate ( meshcoor_n(mesh%nnodes,2), meshcoor_nm1(mesh%nnodes,2) )
  meshcoor_n = mesh%coor
  meshcoor_nm1 = mesh%coor

  call create_vector ( problem, meshvel, vec=4 )
  meshvel%u = 0._dp  ! initialize

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )
  call create_oldvectors ( oldvectors_hat, nsysvec=1 )
  call create_oldvectors ( oldvectors_iter, nsysvec=1, nprob=1, nvec=1 )

! fill oldvectors

  oldvectors%s(1)%p => soln
  oldvectors_hat%s(1)%p => sol_hat  ! alpha0*un+alpha1*un-1
  oldvectors_iter%s(1)%p => sol_iter ! solution at end of previous iteration
  oldvectors_iter%p(1)%p => problem
  oldvectors_iter%v(1)%p => meshvel


! create subscript for velocity sampling on curve 8

  call create_subscript ( mesh, problem, vel8, physqarr=[physqv], curves=[8] )


! fill coefficients for surface advection problem

  call create_coefficients ( coefficients_sf_adv, ncoefi=100, ncoefr=50 )

  coefficients_sf_adv%i = 0
  coefficients_sf_adv%i(1:6) = [ ninti_sf_adv, 2, 0, method, 0, hintpl ]

  coefficients_sf_adv%r = 0

  coefficients_sf_adv%r(4) = deltat
  coefficients_sf_adv%r(5) = betah

! create mesh for surface advection

  call curve_to_1Dmesh ( mesh, mesh%curves(8), dim1D=2, pnts=[2,5], &
    mesh1D=mesh_sf_adv )

  call fill_mesh_parts ( mesh_sf_adv )


! problem definition for surface advection

  call create_input_probdef ( mesh_sf_adv, input_probdef_sf_adv, nvec=2, &
    nphysq=1 )

  input_probdef_sf_adv%vec_elementdof(1)%a(:,1) = 1 ! height function
  input_probdef_sf_adv%vec_elementdof(1)%a(:,2) = 2 ! velocity

  input_probdef_sf_adv%physq = [1]
  input_probdef_sf_adv%probnr = 2

! define problem

  call problem_definition ( input_probdef_sf_adv, mesh_sf_adv, problem_sf_adv )

  call create_sysvector ( problem_sf_adv, sol_sf_adv, rhsd_sf_adv )
  call create_sysvector ( problem_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1 )
  call create_sysvector ( problem_sf_adv, sol_sf_adv_pred )

! fill solution vector with essential boundary conditions

  sol_sf_adv%u = mesh%coor(mesh%points(2),1)
  sol_sf_adv_n%u = sol_sf_adv%u
  sol_sf_adv_pred%u = sol_sf_adv%u

! create system matrix

  call create_sysmatrix_structure ( sysmatrix_sf_adv, mesh_sf_adv, &
    problem_sf_adv )

  call create_sysmatrix_data ( sysmatrix_sf_adv )


  ! create vectors

  call create_vector ( problem_sf_adv, velocity_sf_adv, vec=2 )

  call create_oldvectors ( oldvectors_sf_adv, nsysvec=2, nvec=1 )

  oldvectors_sf_adv%s(1)%p => sol_sf_adv_n    ! corrector at n
  oldvectors_sf_adv%s(2)%p => sol_sf_adv_nm1  ! corrector at nm1

  oldvectors_sf_adv%v(1)%p => velocity_sf_adv ! advection velocity at np1


! create subscript for the height values

  call create ( mesh_sf_adv, problem_sf_adv, hgt )
  call create ( mesh_sf_adv, problem_sf_adv, hgt_end, points=[1,2] )


! allocate arrays for the ALE displacement problem

  allocate ( disp(mesh_sf_adv%nnodes) )


! open files

  open ( unit=11, file='height.out', status='unknown' )

! write initial mesh

  write(filename,'(a,i4.4,a)') 'mesh', 0, '.vtk'
  call write_mesh_vtk ( mesh, filename )


! time stepping

  coefficients_sf_adv%i(5) = 1 ! start with first-order scheme
  tn = 0

  do step = 1, numtimesteps

    tnp1 = tn + deltat

    write(*,'(/a,i0,a,es12.4,a/)') ' ** time step = ', step, &
                                   ' time = ', tnp1, ' ** '
    if ( step >= 2 ) then

      coefficients_sf_adv%i(5) = 2 ! second-order scheme

!     predict position of the surface and adapt mesh accordingly

      call update_mesh_surface_predictor

    end if

!   build (assemble) matrix/vector for velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients )

!   surface tension

    if ( surface_tension ) then

!     build surface integral on free surface

      coefficients%r(19) = gamma_fa

      call add_boundary_elements ( mesh, problem, rhsd, &
        elemsub=surface_tension_curve, &
        curve=8, coefficients=coefficients, physq=[physqv] )

!     build surface integral on walls

      if ( symmetric ) then
        curvew = 6  ! only upper wall
      else
        curvew = 10 ! both walls
      end if

      coefficients%r(19) = dgamma_w

      call add_boundary_elements ( mesh, problem, rhsd, &
        elemsub=surface_tension_curve, &
        curve=curvew, coefficients=coefficients, physq=[physqv] )

    end if

!   build instationary term ( rho du/dt )
    if ( step <= numeulertimesteps ) then
!     implicit Euler
      sol_hat%u(vel%s) = soln%u(vel%s)
      mfac = 1._dp
    else
!     second-order implicit Gear after Euler time steps
      sol_hat%u(vel%s) = alpha0*soln%u(vel%s) + alpha1*solnm1%u(vel%s)
      mfac = gamma0
    end if

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=inertia_elem_dudt, &
      coefficients=coefficients, oldvectors=oldvectors_hat, &
      physqcol=[physqv], physqrow=[physqv], &
      factormat=mfac, addmatvec=.true. )

!   store sysmatrix & rhsd for advection iteration

    call copy ( sysmatrix, sysmatrix_stored )
    call copy ( rhsd, rhsd_stored )

!   fill solution vector with essential boundary conditions

    sol%u = 0
!   wall
    call fill_sysvector ( mesh, problem, sol, curve1=6, physq=physqv, degfd=3, &
      func=vel_theta, funcnr=1 )
!   bottom
    call fill_sysvector ( mesh, problem, sol, curve1=9, physq=physqv, degfd=3, &
      func=vel_theta, funcnr=1 )

!   interation procedure for advection term

    iterstep = 0

    iterate: do

      iterstep = iterstep + 1

      call copy ( sysmatrix_stored, sysmatrix )
      call copy ( rhsd_stored, rhsd )

!     first num_picard steps Picard iteration, then Newton iteration
!     the term un+1.grad un+1
      if ( iterstep <= num_picard ) then
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_picard, coefficients=coefficients, &
          oldvectors=oldvectors_iter, physqcol=[physqv], &
          physqrow=[physqv], addmatvec=.true. )
      else
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=inertia_elem_conv_newton, coefficients=coefficients, &
          oldvectors=oldvectors_iter, physqcol=[physqv], &
          physqrow=[physqv], addmatvec=.true. )
      end if

      call check ( sysmatrix )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

      solver_options_u%real_storage=rs_up
      solver_options_u%integer_storage=is_up

      call solve_system_ma41 ( sysmatrix, rhsd, sol, &
        solver_options=solver_options_u  )

!     compute the difference between iteration steps
      diffv = maxval ( sqrt ( ( sol%u(vel%s) - sol_iter%u(vel%s) )** 2 ) )
      diffp = maxval ( sqrt ( ( sol%u(pres%s) - sol_iter%u(pres%s) )** 2 ) )

      write(*,'(a,i0,2(a,es12.4))') &
        '    iterstep = ', iterstep, ' diffv = ', diffv, ' diffp = ', diffp

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

!   solve surface advection (corrector)

    call solve_surface_height_corrector

    tn = tnp1

    if ( mod(step,plot_mesh_every) == 0 ) then
      write(filename,'(a,i4.4,a)') 'mesh', step, '.vtk'
      call write_mesh_vtk ( mesh, filename )
    end if

  end do

  close(unit=11)

! post-processing

  call create_vector ( problem, velocity, physq=physqv )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, gammadot, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=11

  call derive_vector ( mesh, problem, gammadot, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='moving_boundary2.vtk' )

  call write_vector_vtk ( mesh, problem, filename='moving_boundary2.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true., &
    assume3D=.true. )

  call write_vector_vtk ( mesh, problem, filename='moving_boundary2.vtk', &
    dataname='velocity_vector_planar', vector=velocity, append=.true., &
    degfd=[1,2] )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='moving_boundary2.vtk', dataname='gammadot', append=.true. )


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, rhsd )
  call delete ( sysmatrix, sysmatrix_stored )
  call delete ( coefficients )
  call delete ( vel8 )
  call delete ( pressure, velocity, gammadot )

  deallocate ( meshcoor_initial, meshcoor_n, meshcoor_nm1 )

  call delete ( mesh_sf_adv )
  call delete ( problem_sf_adv )
  call delete ( input_probdef_sf_adv )
  call delete ( sol_sf_adv, rhsd_sf_adv )
  call delete ( sol_sf_adv_pred, sol_sf_adv_n, sol_sf_adv_nm1 )
  call delete ( sysmatrix_sf_adv )
  call delete ( oldvectors_sf_adv )
  call delete ( coefficients_sf_adv )
  call delete ( hgt, hgt_end )

! ALE mesh motion problem

  deallocate ( disp )
  if ( numtimesteps > 1 ) then
    call delete ( problem_update_mesh )
  end if

contains


! Convert curve to a 1D mesh with only a single group.

  subroutine curve_to_1Dmesh ( mesh, geometry, dim1D, pnts, mesh1D )

!   the mesh that contains the geometry
    type(mesh_t), intent(in) :: mesh

!   the geometry of the curve
    type(geometry_t), intent(in) :: geometry

!   the coordinate direction that form the 1D curve coordinates
    integer, intent(in) :: dim1D

!   create points in the new mesh from the points in mesh
    integer, dimension(:), intent(in) :: pnts

!   the new mesh created
    type(mesh_t), intent(inout) :: mesh1D

    real(dp), dimension(size(pnts),1) :: points

!   make skeleton mesh
    call mesh_skeleton ( mesh=mesh1D, nnodes=geometry%nnodes, &
      nelem=geometry%nelem, elshape=geometry%element%elshape, &
      ndim=1, callname='curve_to_1Dmesh' )

!   fill coordinates
    mesh1D%coor(:,1) = mesh%coor(geometry%nodes,dim1D)

!   fill topology
    mesh1D%topology(1)%a = geometry%topology(:,:,1)

!   add points

    points(:,1) = mesh%coor(mesh%points(pnts),dim1D)

    call add_to_mesh ( mesh1D, points=points )

  end subroutine curve_to_1Dmesh

! velocity theta function

  function vel_theta ( nr, xin )

    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: xin
    real(dp) :: vel_theta

    select case ( nr )

    case (1)

!     rotation with omega

      vel_theta  = omega * xin(2)

    case default

      write(*,'(/a,i0/)') 'Error vel_theta: wrong function number: ', nr
      stop

    end select

  end function vel_theta

! body force function

  function vfunc ( n, nr, xin )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: xin
    real(dp), dimension(n) :: vfunc

    select case ( nr )

    case (1)

!     gravity

      vfunc(1)  = - rho * grv
      vfunc(2:) = 0

    case default

      write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
      stop

    end select

  end function vfunc


! update the (ALE) mesh based on a predictor for the surface position

  subroutine update_mesh_surface_predictor

!   predictor for surface height
    sol_sf_adv_pred%u = 2._dp*sol_sf_adv_n%u - sol_sf_adv_nm1%u

    meshcoor_nm1 = meshcoor_n
    meshcoor_n = mesh%coor

    disp = sol_sf_adv_pred%u - mesh%coor(mesh%curves(8)%nodes,1)

!   update nodes of the mesh (only in the first direction)

    call update_mesh_nodes ( mesh, problem_update_mesh, disp, &
      el=mesh_sf_adv%points([1,2]), eg=mesh%points([2,5]), cf=8 )

    call find_bounds_blocks ( mesh )

!   mesh velocity

    meshvel%u = reshape ( transpose ( &
        ( 1.5_dp*mesh%coor - 2*meshcoor_n + 0.5_dp*meshcoor_nm1 ) / deltat ), &
                         [2*mesh%nnodes] )

  end subroutine update_mesh_surface_predictor


! solve convection equation for the surface height (corrector)

  subroutine solve_surface_height_corrector

    use postprocessing_m

    real(dp) :: max_height, min_height, end_height(2)
    type(solver_options_ma41_t) :: solver_options_h
    real(dp), dimension(:,:), allocatable :: tmp

!   get velocities on the moving curve

    nnodes = mesh_sf_adv%nnodes
    tmp = reshape ( sol%u(vel8%s), [3,nnodes] )
    velocity_sf_adv%u = reshape ( tmp([2,1],:), [2*nnodes] )

    call build_system ( mesh_sf_adv, problem_sf_adv, sysmatrix_sf_adv, &
      rhsd_sf_adv, elemsub=surface_advection_elem, &
      oldvectors=oldvectors_sf_adv, coefficients=coefficients_sf_adv )

    call check_filled_sysmatrix ( sysmatrix_sf_adv )

    call add_effect_of_essential_to_rhs ( problem_sf_adv, sysmatrix_sf_adv, &
      sol_sf_adv, rhsd_sf_adv )


!   MA41 solver storage

    solver_options_h%integer_storage = 2.0
    solver_options_h%real_storage    = 2.0

!   solve system

    call solve_system_ma41 ( sysmatrix_sf_adv, rhsd_sf_adv, sol_sf_adv, &
                             solver_options=solver_options_h )


    call copy ( sol_sf_adv_n, sol_sf_adv_nm1 )
    call copy ( sol_sf_adv, sol_sf_adv_n )

!   compute the maximum and end height

    max_height = maxval ( sol_sf_adv%u )
    min_height = minval ( sol_sf_adv%u )
    end_height = sol_sf_adv%u(hgt_end%s(1:2))

!   write swell height

    write(11,'(i6,6es16.8)') step, step*deltat, &
                             max_height, min_height, end_height
    print '(i6,6es16.8)', step, step*deltat, &
                          max_height, min_height, end_height

  end subroutine solve_surface_height_corrector

end program moving_boundary2
