! Extrudate swell problem of a Newtonian fluid (Stokes)
! Planar or axisymmetrical flow

program extrudate_swell1

  use tfem_m
  use stokes_elements_m
  use hsl_ma41_m
  use surface_advection_elements_m

  implicit none


! constants flow problem

  integer, parameter :: &
    uintpl = 8,     & ! Q2 velocities
    pintpl = 4,     & ! Q1 pressures
    ointpl = 6,     & ! P2 shape of object elements
    coorsys = 0,    & ! planar Cartesian (0) or axisymmetric (1) coor. system
    physqvel = 1,   & ! physical quantity nr of the velocities
    physqpress = 2, & ! physical quantity nr of the pressures
    gauss = 3         ! 3x3 integration of quads

! constants surface advection

  integer, parameter :: &
    hintpl = 2,         & ! interpolation for the height, P1 lines
    ninti_sf_adv = 2,   & ! number of Gauss points
    method = 1            ! discretization method 0: Galerkin, 1: SUPG

   real(dp), parameter :: &
     beta_sf_adv = 1.0_dp        ! beta parameter for SUP

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options_u
  type(oldvectors_t) :: oldvectors

  ! 1D height function for surface advection
  type(meshgen_options_t) :: mesh_options
  type(mesh_t) :: mesh_sf_adv
  type(input_probdef_t) :: input_probdef_sf_adv
  type(problem_t) :: problem_sf_adv
  type(sysmatrix_t) :: sysmatrix_sf_adv
  type(sysvector_t), target :: sol_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1
  type(sysvector_t) :: rhsd_sf_adv
  type(oldvectors_t) :: oldvectors_sf_adv
  type(coefficients_t) :: coefficients_sf_adv
  type(vector_t), target :: velocity_sf_adv
  type(subscript_t) :: hgt

  type(vector_t), target :: meshvel  ! computed, but not used in this program
  real(dp), allocatable, dimension(:,:) :: meshcoor_initial


! variables

  integer :: &
    numtimesteps = 2

  real(dp) :: &
    eta = 1.0_dp,    & ! viscosity
    deltat = 5.e-2_dp, & ! time step
    rs_gup = 1.5_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.0_dp      ! integer_storage gradient-velocity-pressure LU (HSL)

  integer :: step, i, j, nnodes
  real(dp) :: flowrate = -1._dp

  real(dp), allocatable, dimension(:,:) :: coor
  real(dp) :: initial_h

  real(dp), allocatable, dimension(:,:) :: meshcoor_n, meshcoor_nm1

  type(meshgen_options_t) :: mesh_options1
  type(mesh_t) :: mesh1, mesh2

  integer :: elshape=6
  real(dp) :: L1=5, L2=5
  real(dp) :: R=1
  integer :: n1=8, n2=8, n3=12

  real(dp) :: f1=7, f2=6, f3=8


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,     0, &
      physqvel, physqpress, 0,     0, gauss, &
      gauss ]
  coefficients%i(23) = coorsys

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(6) = flowrate

  coefficients%i(35) = ointpl

! mesh

! 1

  call set_mesh_options ( mesh_options1, elshape=elshape, nx=n1, ny=n3, ox=-L1, &
    lx=L1, ly=R, ratio=[3,3,1,1], factor=[f1,f3,f1,f3] )

  call quadrilateral2d ( mesh1, mesh_options1 )
  !call plot_points_curves ( plot_options, mesh1, 'curves1.fig' )

! 2

  call set_mesh_options ( mesh_options1, elshape=elshape, nx=n2, ny=n3, &
    lx=L2, ly=R, ratio=[1,3,3,1], factor=[f2,f3,f2,f3] )

  call quadrilateral2d ( mesh2, mesh_options1 )

  call mesh_merge ( mesh1, mesh2, mesh, curve1=2, curve2=-4 )

  !call plot_points_curves ( plot_options, mesh, 'curves2.fig' )
  call delete ( mesh1, mesh2 )

! add curves

! wall+slip surface
  call add_to_mesh ( mesh, curve=[-3,-7] ) ! curve 8
! centerline
  call add_to_mesh ( mesh, curve=[1,5] ) ! curve 9
! full boundary
  call add_to_mesh ( mesh, curve=[1,5,6,7,3,4] ) ! curve 10

  call add_to_mesh ( mesh, object='curve', objectcurve=4, topology=.true., &
    intrule=3 )

  call fill_mesh_parts ( mesh )


! problem definition of velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =  &
      reshape ( [ 2,2,2,2,2,2,2,2,2,   &  ! velocity
                  1,0,1,0,1,0,1,0,0,   &  ! pressure
                  1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                [9,3] )

  input_probdef%physq = [1,2]
  input_probdef%probnr = 1

! Dirichlet boundary conditions

! center line
  call define_essential ( mesh, input_probdef, &
    curve1=9, physq=physqvel, degfd=[0,1] )
! outflow
  call define_essential ( mesh, input_probdef, &
    curve1=6, physq=physqvel, degfd=[0,1] )
! wall
  call define_essential ( mesh, input_probdef, &
    curve1=3, physq=physqvel )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=4, nglobalc=1 )

  call problem_definition ( input_probdef, mesh, problem )


! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )


! fill solution vector with essential boundary conditions

  sol%u = 0

! create oldvectors for sampling

  call create_oldvectors ( oldvectors, nsysvec=1 )
  oldvectors%s(1)%p => sol

! create system matrix of gradient/velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! define some arrays and vector for the ALE mesh movement

  allocate ( meshcoor_initial(mesh%nnodes,2) )
  meshcoor_initial = mesh%coor

  allocate ( meshcoor_n(mesh%nnodes,2), meshcoor_nm1(mesh%nnodes,2) )
  meshcoor_n = mesh%coor
  meshcoor_nm1 = mesh%coor

  call create_vector ( problem, meshvel, physq=physqvel )
  meshvel%u = 0._dp  ! initialize


! fill coefficients for surface advection problem

  call create_coefficients ( coefficients_sf_adv, ncoefi=100, ncoefr=50 )

  coefficients_sf_adv%i = 0
  coefficients_sf_adv%i(1) = ninti_sf_adv
  coefficients_sf_adv%i(2) = 2 ! velocity given by nodal values
  coefficients_sf_adv%i(4) = method ! 0: Galerkin, 1: SUPG
  coefficients_sf_adv%i(5) = 1  ! start with first-order time discretization
  coefficients_sf_adv%i(6) = hintpl ! height interpolation
  coefficients_sf_adv%i(9) = 3 ! numerical table for Gauss

  coefficients_sf_adv%r = 0
  coefficients_sf_adv%r(4) = deltat
  coefficients_sf_adv%r(5) = beta_sf_adv

! create mesh for surface advection

  mesh_options%elshape = 1   ! two-node line elements

  mesh_options%nx = mesh%curves(7)%nelem * 2

  call line1d ( mesh_sf_adv, mesh_options )

! set coordinates

  do i = 1, mesh%curves(7)%nnodes
    j = mesh%curves(7)%nnodes - i + 1 ! reverse order to get increasing x
    mesh_sf_adv%coor(i,1) = mesh%coor(mesh%curves(7)%nodes(j),1)
  end do

  call fill_mesh_parts ( mesh_sf_adv )

! add object for sampling velocity in height function advection

  initial_h = mesh%coor(mesh%curves(7)%nodes(1),2)

  allocate ( coor(mesh_sf_adv%nnodes,2) )

  coor(:,1) = mesh_sf_adv%coor(1:mesh_sf_adv%nnodes,1)
  coor(:,2) = initial_h

  warn_add_to_mesh_after_meshgen_parts = .false.
  call add_to_mesh ( mesh, object='coordinates', coor=coor )
  call fill_mesh_parts_objects ( mesh, object1=mesh%nobjects )
  warn_add_to_mesh_after_meshgen_parts = .true.

  deallocate ( coor )


! problem definition for surface advection

  call create_input_probdef ( mesh_sf_adv, input_probdef_sf_adv, nvec=2, &
    nphysq=1 )

  input_probdef_sf_adv%vec_elementdof(1)%a = &
      reshape ( [ 1,1,    &  ! height function
                  2,2 ],  &  ! velocity
                [2,2] )

  input_probdef_sf_adv%physq = [1]
  input_probdef_sf_adv%probnr = 2

  call define_essential ( mesh_sf_adv, input_probdef_sf_adv, point=1, physq=1 )

! define problem

  call problem_definition ( input_probdef_sf_adv, mesh_sf_adv, problem_sf_adv )

  call create_sysvector ( problem_sf_adv, sol_sf_adv, rhsd_sf_adv )
  call create_sysvector ( problem_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1 )

! fill solution vector with essential boundary conditions

  sol_sf_adv%u = initial_h
  sol_sf_adv_n%u = sol_sf_adv%u

  call fill_sysvector ( mesh_sf_adv, problem_sf_adv, sol_sf_adv, point=1, &
    physq=1, value=initial_h )

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


! time stepping

  do step = 1, numtimesteps

    if ( step >= 2 ) then

      coefficients_sf_adv%i(5) = 2

!     predict position of the surface and adapt mesh accordingly

      call update_mesh_surface_predictor

    end if

!   build (assemble) matrix/vector for velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients )

!   open boundary

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub1=stokes_open_boundary, &
      coefficients=coefficients, object=1, addmatvec=.true. )

!   imposed flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients )


    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

!   solve surface advection (corrector)

    call solve_surface_height_corrector

  end do

! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( oldvectors )
  call delete ( hgt )

  call delete ( coefficients )

contains


! update the (ALE) mesh based on a predictor for the surface position

  subroutine update_mesh_surface_predictor

    integer :: elem, nod(2), i

    real(dp) :: x1, x2, y1, y2, yp, x(2), u(2)
    real(dp) :: displ, factor
    logical :: find


    ! predictor for surface height
    sol_sf_adv%u = 2._dp*sol_sf_adv_n%u - sol_sf_adv_nm1%u


    meshcoor_nm1 = meshcoor_n
    meshcoor_n = mesh%coor


    do i = 1, mesh%nnodes

      if ( mesh%coor(i,1) <= 0 ) cycle  ! inside die

!     search for the element

      find  = .false.

      do elem = 1, mesh_sf_adv%nelem

        nod = mesh_sf_adv%topology(1)%a(1:2,elem)
        x = mesh_sf_adv%coor(nod,1)

        x1 = x(1); x2 = x(2)

        if ( (mesh%coor(i,1) >= x1) .and. (mesh%coor(i,1) <= x2) ) then

!         get surface height for element elem

          u = sol_sf_adv%u( &
                 problem_sf_adv%degfdperm(problem_sf_adv%nodnumdegfd(nod)+1,2) )
          y1 = u(1); y2 = u(2)

!         height at position mesh%coor(i,1)
          yp = y1 + (y2-y1)*(mesh%coor(i,1)-x1)/(x2-x1)

          displ = yp - initial_h
          factor = meshcoor_initial(i,2) / initial_h

          mesh%coor(i,2) = meshcoor_initial(i,2) + displ*factor

          find = .true.

          exit

        end if

      end do

      if ( .not. find ) then
        print*, 'Error update_mesh_surface_predictor:'
        print*, 'Finding position in surface failed for node:'
        print*, i, mesh%coor(i,1), mesh%coor(i,2)
        print*, 'Program Stop!'
        stop
      end if

    end do

    call find_bounds_blocks ( mesh )

!   mesh velocity

    meshvel%u = reshape ( transpose ( &
        ( 1.5_dp*mesh%coor - 2*meshcoor_n + 0.5_dp*meshcoor_nm1 ) / deltat ), &
                         [2*mesh%nnodes] )

  end subroutine update_mesh_surface_predictor


! solve convection equation for the surface height (corrector)

  subroutine solve_surface_height_corrector

    use postprocessing_m

    type(sample_t) :: sample
    real(dp) :: max_height
    type(solver_options_ma41_t) :: solver_options_h

    nnodes = mesh_sf_adv%nnodes

    mesh%objects(2)%coor(1:nnodes,1) = mesh_sf_adv%coor(1:nnodes,1)
    mesh%objects(2)%coor(1:nnodes,2) = sol_sf_adv%u(hgt%s)

    call find_refcoor_objects ( mesh, object1=2 )

    if ( any (mesh%objects(2)%grpelm(:,1) == 0) ) then
      print*, 'Error: no reference coordinates in surface advection corrector.'
      print*, 'Program Stop!'
      stop
    end if

    call fill_sample ( mesh, problem, sample, ndegfd=2, object=2, &
      elemsub=stokes_sample_velocity, coefficients=coefficients, &
      oldvectors=oldvectors )

    velocity_sf_adv%u = reshape ( transpose(sample%u), [2*nnodes] )


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

!   compute the maximum height

    max_height = maxval ( sol_sf_adv%u )

    print '(i6,3es25.13)', step, step*deltat, max_height

  end subroutine solve_surface_height_corrector


end program extrudate_swell1
