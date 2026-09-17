! Moving boundary due to capillary forces. Perfect slip on the walls.
! Surface tension on fluid-air interface and fluid-wall interface.
! Newtonian fluid (Stokes).

program moving_boundary1

  use tfem_m
  use stokes_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m
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
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(coefficients_t) :: coefficients
  type(vector_t) :: velocity, pressure, gammadot
  type(solver_options_ma57_t) :: solver_options_u
  type(oldvectors_t) :: oldvectors
  type(subscript_t) :: vel8

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

  type(vector_t), target :: meshvel  ! computed, but not used in this program
  real(dp), allocatable, dimension(:,:) :: meshcoor_initial

! ALE mesh motion problem

  type(problem_t), target ::  problem_update_mesh
  real(dp), allocatable, dimension(:) :: disp

! variables

  integer :: &
    coorsys = 0, & ! planar Cartesian (0) or axisymmetric (1) coor. system
    vfuncnr = 1, & ! function nr for gravity 0=none 1=rho g (-1,0)
    numtimesteps = 250,   & ! number of time steps
    plot_mesh_every = 10    ! plot mesh every .. steps

  real(dp) :: &
    rho = 1.0_dp,      & ! density
    grv = 1.0_dp,      & ! gravity constant g
    eta = 1.0_dp,      & ! viscosity
    gamma_fa = 0.1_dp, & ! surface tension coefficient fluid-air interface
    dgamma_w = 0.01_dp,& ! surface tension difference fluid-air wrt the wall
    deltat = 5.e-2_dp, & ! time step
    betah = 0.5_dp,    & ! beta parameter for SUPG of the height function
    rs_up = 1.5_dp,    & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.0_dp       ! integer_storage velocity-pressure LU (HSL)

  integer :: vertices(4) = [1,3,5,7]
  integer :: step, nnodes, curvew

  real(dp), allocatable, dimension(:,:) :: meshcoor_n, meshcoor_nm1

  character(len=20) :: filename

  logical :: surface_tension = .true.  ! include surface tension


! namelist for input of variables; read from standard input

  namelist /comppar/ numtimesteps, plot_mesh_every, coorsys, &
    rho, grv, eta, gamma_fa, dgamma_w, deltat, rs_up, is_up, surface_tension

  read ( unit=*, nml=comppar )


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=600, ncoefr=100 )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl, pintpl, 0, 0,     0, &
      physqv, physqp, 0, 0, gauss, &
      gauss ]
  coefficients%i(14) = vfuncnr
  coefficients%i(23) = coorsys

  coefficients%r = 0
  coefficients%r(1) = eta

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

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 2  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(vertices,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2]
  input_probdef%probnr = 1

! Dirichlet boundary conditions

! bottom
  call define_essential ( mesh, input_probdef, &
    curve1=9, physq=physqv )
! walls (perfect slip)
  call define_essential ( mesh, input_probdef, &
    curve1=10, physq=physqv, degfd=[0,1] )

  call problem_definition ( input_probdef, mesh, problem )


! create system vectors for velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )


! fill solution vector with essential boundary conditions

  sol%u = 0


! create oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )
  oldvectors%s(1)%p => sol

! create system matrix of velocity/pressure problem

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! define some arrays and vector for the ALE mesh movement

  allocate ( meshcoor_initial(mesh%nnodes,2) )
  meshcoor_initial = mesh%coor

  allocate ( meshcoor_n(mesh%nnodes,2), meshcoor_nm1(mesh%nnodes,2) )
  meshcoor_n = mesh%coor
  meshcoor_nm1 = mesh%coor

  call create_vector ( problem, meshvel, physq=physqv )
  meshvel%u = 0._dp  ! initialize

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

  do step = 1, numtimesteps

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

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve velocity/pressure problem

    solver_options_u%real_storage=rs_up
    solver_options_u%integer_storage=is_up

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

!   solve surface advection (corrector)

    call solve_surface_height_corrector

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

  coefficients%i(13)=8

  call derive_vector ( mesh, problem, gammadot, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='moving_boundary1.vtk' )

  call write_vector_vtk ( mesh, problem, filename='moving_boundary1.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='moving_boundary1.vtk', dataname='gammadot', append=.true. )


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( oldvectors )
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
    tmp = reshape ( sol%u(vel8%s), [2,nnodes] )
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

end program moving_boundary1
