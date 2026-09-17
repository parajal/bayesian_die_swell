! Extrudate swell problem of a Newtonian fluid (Stokes)
! Optionally include surface tension (using an object).
! Planar or axisymmetrical flow

program extrudate_swell1

  use tfem_m
  use math_defs_m
  use stokes_elements_m
  use hsl_ma41_m
  use io_utils_m
  use surface_advection_elements_m
  use subs_extrudate_swell_m  ! subs for extrudate swell
  use figplot_m

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
    beta_sf_adv = 1.0_dp        ! beta parameter for SUPG


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
  type(sysvector_t) :: rhsd_sf_adv, sol_sf_adv_pred
  type(oldvectors_t) :: oldvectors_sf_adv
  type(coefficients_t) :: coefficients_sf_adv
  type(vector_t), target :: velocity_sf_adv
  type(plot_options_t) :: plot_options
  type(subscript_t) :: hgt, hgt_end

  type(vector_t), target :: meshvel  ! computed, but not used in this program
  real(dp), allocatable, dimension(:,:) :: meshcoor_initial


! variables

  integer :: &
    numtimesteps = 250,   & ! number of time steps
    plot_mesh_every = 10, & ! plot mesh every .. steps
    step0 = 0,            & ! initial step number
    restart = 0             ! restart=1: restart,
                            ! restart=2: restart+Euler first time step

  real(dp) :: &
    eta = 1.0_dp,      & ! viscosity
    gammac = 0.1_dp,   & ! surface tension coefficient
    deltat = 5.e-2_dp, & ! time step
    time0 = 0._dp,     & ! initial time
    U_avg = 1._dp,     & ! average velocity at the entry
    rs_up = 1.5_dp,    & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.0_dp       ! integer_storage velocity-pressure LU (HSL)

  integer :: step, i, j, nnodes, obj_ob, obj_sh, obj_st
  real(dp) :: flowrate, R, H

  real(dp), allocatable, dimension(:,:) :: coor
  real(dp) :: initial_h

  real(dp), allocatable, dimension(:,:) :: meshcoor_n, meshcoor_nm1

  character(len=20) :: filename

  logical :: surface_tension = .false.  ! include surface tension

  logical :: ctime = .false. ! continue time in output after restart

  logical :: file_append = .false.  ! append radius output to existing
                                    ! file after restart

  logical, parameter :: mesh_plot_surface = .false.


! namelist for input of variables; read from standard input

  namelist /comppar/ numtimesteps, plot_mesh_every, restart, ctime, &
    file_append, eta, gammac, deltat, rs_up, is_up, U_avg, surface_tension

  read ( unit=*, nml=comppar )


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=600, ncoefr=100 )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,     0, &
      physqvel, physqpress, 0,     0, gauss, &
      gauss ]
  coefficients%i(23) = coorsys

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(19) = gammac

  coefficients%i(35) = ointpl

! read mesh

  call read_mesh ( mesh, filename='mesh2b.out' )

! add object for open boundary condition

  call add_to_mesh ( mesh, object='curve', objectcurve=4, topology=.true., &
    intrule=3 )
  obj_ob = mesh%nobjects

  if ( surface_tension ) then

!   add object for surface tension

    call add_to_mesh ( mesh, object='curve', objectcurve=7, topology=.true., &
      intrule=3 )
    obj_st = mesh%nobjects

  end if

  call fill_mesh_parts ( mesh )

! set flowrate

  if ( coorsys == 0 ) then
    h = mesh%coor(mesh%points(4),2) ! height is given by y-coordinate of P4
    flowrate = - h * U_avg
  else if ( coorsys == 1 ) then
    R = mesh%coor(mesh%points(4),2) ! radius is given by r-coordinate of P4
    flowrate = - pi * R**2 * U_avg
  end if

  coefficients%r(6) = flowrate


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
! outflow: ur=0
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
  obj_sh = mesh%nobjects
  call fill_mesh_parts_objects ( mesh, object1=obj_sh )

  deallocate ( coor )

! plot mesh with object for surface

  if ( mesh_plot_surface ) then

   call plot_mesh ( plot_options, mesh, 'mesh.fig' )

   plot_options%objectpointcolor = 4   ! red color
   plot_options%objectpointsize = 0.2

   call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true., &
     object1=obj_sh )

  end if


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
  call create_sysvector ( problem_sf_adv, sol_sf_adv_pred )

! fill solution vector with essential boundary conditions

  sol_sf_adv%u = initial_h
  sol_sf_adv_n%u = sol_sf_adv%u
  sol_sf_adv_pred%u = sol_sf_adv%u

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
  call create ( mesh_sf_adv, problem_sf_adv, hgt_end, points=[2] )


! open files and restart

  call open_files_and_restart


! time stepping

  coefficients_sf_adv%i(5) = 1 ! start with first-order scheme

  do step = 1, numtimesteps

    if ( step >= 2 .or. restart == 1 ) then

      coefficients_sf_adv%i(5) = 2 ! second-order scheme

!     predict position of the surface and adapt mesh accordingly

      call update_mesh_surface_predictor

    end if

!   build (assemble) matrix/vector for velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients )

!   open boundary

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub1=stokes_open_boundary, &
      coefficients=coefficients, object=obj_ob, addmatvec=.true. )

!   imposed flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients )

!   surface tension

    if ( surface_tension ) then

!     update coordinates of object obj_st

      mesh%objects(obj_st)%coor = mesh%coor(mesh%curves(7)%nodes,:)
      call find_refcoor_objects ( mesh, object1=obj_st )

!     build surface integral

      call build_system ( mesh, problem, sysmatrix, rhsd,&
        elemsub1=surface_tension_object, &
        physqrow=[physqvel], physqcol=[physqvel], object=obj_st,&
        coefficients=coefficients, addmatvec=.true., buildmatrix=.false. )

!     build boundary force term

      coefficients%i(501) = 7  ! free surface = curve 7
      coefficients%i(502) = 1  ! element number on the curve

      call add_boundary_elements_point ( mesh, problem, rhsd, point=6, &
        elemsub=surface_tension_boundary_point1, physq=[physqvel], &
        coefficients=coefficients, buildmatrix=.false. )

    end if

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve velocity/pressure problem

    solver_options_u%real_storage=rs_up
    solver_options_u%integer_storage=is_up

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

!   solve surface advection (corrector)

    call solve_surface_height_corrector

    if ( mod(step0+step,plot_mesh_every) == 0 ) then
      write(filename,'(a,i4.4,a)') 'mesh_swell_', step0+step, '.fig'
      call plot_mesh ( plot_options, mesh, filename )
    end if

  end do

  close(unit=11)

! write data for post-processing and restart

  call write_mesh ( mesh, filename='mesh.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol%u
  write(10) sol_sf_adv_pred%u, sol_sf_adv_n%u, sol_sf_adv_nm1%u
  write(10) mesh%coor, meshcoor_n
  write(10) step0 + numtimesteps, time0 + numtimesteps*deltat

  close(unit=10)

! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( oldvectors )
  call delete ( coefficients )

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

contains


! update the (ALE) mesh based on a predictor for the surface position

  subroutine update_mesh_surface_predictor

    integer :: elem, nod(2), i

    real(dp) :: x1, x2, y1, y2, yp, x(2), u(2)
    real(dp) :: displ, factor
    logical :: find


    ! predictor for surface height
    sol_sf_adv_pred%u = 2._dp*sol_sf_adv_n%u - sol_sf_adv_nm1%u


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

          u = sol_sf_adv_pred%u( &
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
    real(dp) :: max_height, end_height
    type(solver_options_ma41_t) :: solver_options_h

    nnodes = mesh_sf_adv%nnodes

    mesh%objects(obj_sh)%coor(1:nnodes,1) = mesh_sf_adv%coor(1:nnodes,1)
    mesh%objects(obj_sh)%coor(1:nnodes,2) = sol_sf_adv_pred%u(hgt%s)

    call find_refcoor_objects ( mesh, object1=obj_sh )

    if ( any (mesh%objects(obj_sh)%grpelm(:,1) == 0) ) then
      print*, 'Error: no reference coordinates in surface advection corrector.'
      print*, 'Program Stop!'
      stop
    end if

    call fill_sample ( mesh, problem, sample, ndegfd=2, object=obj_sh, &
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

!   compute the maximum and end height

    max_height = maxval ( sol_sf_adv%u )
    end_height = sol_sf_adv%u(hgt_end%s(1))

!   write swell height

    write(11,'(i6,4es16.8)') step0+step, time0+step*deltat, &
                             max_height, end_height
    print '(i6,4es16.8)', step0+step, time0+step*deltat, &
                          max_height, end_height

  end subroutine solve_surface_height_corrector


  subroutine open_files_and_restart

    integer :: endstep
    real(dp) :: endtime

    if ( restart >= 1 ) then

!     restart: read solution from file

      open ( unit=10, form='unformatted', file='data.out' )

      read(10) sol%u
      read(10) sol_sf_adv_pred%u, sol_sf_adv_n%u, sol_sf_adv_nm1%u
      read(10) mesh%coor, meshcoor_n
      read(10) endstep, endtime

      close(unit=10)

!     coordinates of mesh have been changed
      call find_bounds_blocks ( mesh )
      call find_refcoor_objects ( mesh )

      if ( ctime ) then
        time0 = endtime
        step0 = endstep
      end if

!     open file for writing the radius as a function of time

      if ( file_append ) then
        open ( unit=11, file='swell_height.out', status='old', &
               position='append' )
      else
        open ( unit=11, file='swell_height.out', status='unknown' )
      end if

    else

!     fresh start

      open ( unit=11, file='swell_height.out', status='unknown' )

    end if

  end subroutine open_files_and_restart

end program extrudate_swell1
