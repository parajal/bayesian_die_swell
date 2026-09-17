! Extrudate swell problem of a Newtonian fluid (Stokes)
! Planar or axisymmetrical flow
! XFEM approach
! Optionally include surface tension

program extrudate_swell4

  use tfem_m
  use math_defs_m
  use stokes_elements_m
  use hsl_ma41_m
  use io_utils_m
  use surface_advection_elements_m
  use usergauss_xfem_m
  use subs_extrudate_swell_m  ! subs for extrudate swell
  use figplot_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    ointpl = 6,         & ! P2 shape of object elements
    coorsys = 0,        & ! planar Cartesian (0) or axisymm. (1) coor. system
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    numsplit = 3,       & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit
    numsplitmin = 2,    & ! relative size of all subelements will be at least
                          ! as small as 1/2^numsplitmin
    numsplit_plot = 2,  & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit for plotting
    gausse = 3,         & ! Gauss integration for the eltree subelements
    gausssub_small = 8, & ! Gauss integration for the submesh subelements for
                          ! small integration volumes
    gausssub = 4,       & ! Gauss integration for the submesh subelements
    gauss_ie = 3,       & ! Gauss integration for the interface elements
    gauss = 3             ! 3x3 integration of quads

  real(dp), parameter :: &
    split_threshold = 1.e-12_dp, & ! levelset value for being "close" enough
                                   ! to the interface for tree splitting
    epsvol = 1e-6_dp, & ! elements with integration area smaller are removed
                        ! from the eltree_array. This means that:
                        ! 1) these elements are "fully inside" and no
                        !    "virtual/extended degrees" are generated.
                        ! 2) the interface integration is ignored and therefore
                        !    the interface has a small "hole".
                        ! Note, that epsvol is relative to the area of an
                        ! element, so epsvol=1 is a full element.
    epsvol_small = 1e-2_dp ! elements with integration area smaller
                           ! are integrated with integration rule
                           ! gausssub_small, otherwise with gaussub

! definitions

  type(mesh_t), target :: mesh
  type(mesh_t) :: mesh_plot, mesh1, mesh_surface
  type(input_probdef_t) :: input_probdef, input_probdef_plot
  type(problem_t), target :: problem, problem_plot
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options_u

  integer, dimension(:), allocatable :: nodes
  real(dp), dimension(:,:), allocatable :: coor

! array of pointers to an eltree

  type(eltree_p), target, dimension(:,:), allocatable :: eltree_array


! constants surface advection

  integer, parameter :: &
    hintpl = 2,         & ! interpolation for the height, P1 lines
    ninti_sf_adv = 2,   & ! number of Gauss points
    method = 1            ! discretization method 0: Galerkin, 1: SUPG

  real(dp), parameter :: &
    beta_sf_adv = 1.0_dp        ! beta parameter for SUPG

! 1D height function for surface advection

  type(mesh_t), target :: mesh_sf_adv
  type(problem_t), target :: problem_sf_adv
  type(meshgen_options_t) :: mesh_options
  type(input_probdef_t) :: input_probdef_sf_adv
  type(sysmatrix_t) :: sysmatrix_sf_adv
  type(sysvector_t), target :: sol_sf_adv, sol_sf_adv_n, sol_sf_adv_nm1
  type(sysvector_t), target :: sol_sf_adv_pred
  type(sysvector_t) :: rhsd_sf_adv
  type(oldvectors_t) :: oldvectors_sf_adv
  type(coefficients_t) :: coefficients_sf_adv
  type(vector_t), target :: velocity_sf_adv
  type(plot_options_t) :: plot_options
  type(subscript_t) :: hgt, hgt_end


! variables

  integer :: &
    numtimesteps = 300,   & ! number of time steps
    plot_mesh_every = 10, & ! plot mesh every .. steps
    step0 = 0,            & ! initial step number
    restart = 0             ! restart=1: restart,
                            ! restart=2: restart+Euler first time step

  real(dp) :: &
    eta = 1.0_dp,      & ! viscosity
    gammac = 0.1_dp,   & ! surface tension coefficient
    deltat = 2.e-2_dp, & ! time step
    time0 = 0._dp,     & ! initial time
    U_avg = 1._dp,     & ! average velocity at the entry
    rs_up = 1.5_dp,    & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.0_dp,    & ! integer_storage velocity-pressure LU (HSL)
    eps_h0 = 1.e-10_dp   ! initial height of the free boundary points is given
                         ! by h0-eps_h0 to avoid that interface is exactly
                         ! at element edges.


  integer :: step, i, j, nnodes, nbx, nby, obj_ob, obj_sh, els_srf, obj_bp
  integer :: n1, n2, nnsf, nsf
  real(dp) :: flowrate, R, H, coor_bp(1,2)=0


  character(len=20) :: filename

  logical :: ctime = .false. ! continue time in output after restart

  logical :: file_append = .false.  ! append radius output to existing
                                    ! file after restart

  logical :: surface_tension = .false.  ! include surface tension


! namelist for input of variables; read from standard input

  namelist /comppar/ numtimesteps, plot_mesh_every, restart, ctime, &
    file_append, eta, gammac, deltat, rs_up, is_up, U_avg, eps_h0, &
    surface_tension

  read ( unit=*, nml=comppar )


! usergauss

  intrule_e = gausse  ! Gauss integration for the eltree subelements
  intrule_sub = gausssub  ! Gauss integration for the eltree subelements
  intrule_sub_small = gausssub_small  ! Gauss integration for the eltree
                                      ! subelements (small integration areas)
  vol_small = epsvol_small  ! elements with integration areas smaller
                            ! are integrated with integration rule
                            ! gausssub_small, otherwise with gaussub


! read mesh

  call read_mesh ( mesh, filename='mesh_esx.out' )

! add object for open boundary condition

  call add_to_mesh ( mesh, object='curve', objectcurve=4, topology=.true., &
    intrule=3 )
  obj_ob = mesh%nobjects

  if ( surface_tension ) then

!   add object for assembling boundary point

    call add_to_mesh ( mesh, object='coordinates', coor=coor_bp )
    obj_bp = mesh%nobjects

  end if

  nbx = nint(sqrt(real(mesh%curves(12)%nnodes)))
  nby = nint(sqrt(real(mesh%curves(11)%nnodes)))

  call add_to_mesh ( mesh, blocks=[nbx,nby] )

  call fill_mesh_parts ( mesh )

  lmesh => mesh  ! supply mesh to mapcoor for mapping coordinates


! set flowrate

  if ( coorsys == 0 ) then
    h = mesh%coor(mesh%points(3),2) ! height is given by y-coordinate of P3
    flowrate = - h * U_avg
  else if ( coorsys == 1 ) then
    R = mesh%coor(mesh%points(3),2) ! radius is given by r-coordinate of P3
    flowrate = - pi * R**2 * U_avg
  end if

! set h0 in functions_xfem

  h0 = mesh%coor(mesh%points(3),2) ! height is given by y-coordinate of P3


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=502 )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,     0,  &
      physqvel, physqpress, 0,     0, gauss,  &
      gauss ]
  coefficients%i(23) = coorsys
  coefficients%i(35) = ointpl
  coefficients%i(37) = 2 ! Gauss points defined by user subroutines
  coefficients%i(41) = gauss_ie ! Gauss points for interface elements

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(6) = flowrate
  coefficients%r(8) = deltat
  coefficients%r(19) = gammac

  coefficients%set_ninti_user => set_ninti_user
  coefficients%set_Gauss_integration_user => set_Gauss_integration_user


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

! Add nodeset for essential boundary conditions outside of the fluid.
! Temporarily added nodeset with node '1'. Actual nodeset will be replaced
! in find_outside_zero_nodes
  call add_to_mesh ( mesh, nodeset='nodes', nodes=[1] )

! center line
  call define_essential ( mesh, input_probdef, &
    curve1=12, physq=physqvel, degfd=[0,1] )
! outflow
  call define_essential ( mesh, input_probdef, &
    curve1=11, physq=physqvel, degfd=[0,1] )
! wall
  call define_essential ( mesh, input_probdef, &
    curve1=3, physq=physqvel )
! outside of the fluid
  call define_essential ( mesh, input_probdef, nodeset1=1 )


! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=4, nglobalc=1 )

  call problem_definition ( input_probdef, mesh, problem )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0


! allocate array of pointers to an eltree (eltree for all elements, one group)

  allocate ( eltree_array(1,mesh%grpnumel(1)), vol(mesh%grpnumel(1)) )

  allocate ( nodes(mesh%nnodes), d(mesh%nnodes) )


! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1, nprob=1, nelta=1 )


! store solution vectors and problem structures

  oldvectors%s(1)%p => sol
  oldvectors%p(1)%p => problem
  oldvectors%ea(1)%p => eltree_array
  ea => eltree_array ! userelmesh does not have oldvectors as an argument


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

  allocate ( coor(mesh_sf_adv%nnodes,2) )
  coor(:,1) = mesh_sf_adv%coor(1:mesh_sf_adv%nnodes,1)
  coor(:,2) = h0

  warn_add_to_mesh_after_meshgen_parts = .false.
! real coordinate should be set in surface_corrector
  call add_to_mesh ( mesh, object='coordinates', coor=coor )
  obj_sh = mesh%nobjects
  call fill_mesh_parts_objects ( mesh, object1=obj_sh )
  warn_add_to_mesh_after_meshgen_parts = .true.

  deallocate ( coor )


! create a mesh for plotting the surface

  call mesh_skeleton ( mesh1, nnodes=mesh_sf_adv%nnodes, &
    nelem=mesh_sf_adv%nelem, elshape=1, ndim=2 )
  mesh1%coor = mesh%objects(obj_sh)%coor
  mesh1%topology = mesh_sf_adv%topology

! merge normal mesh and surface mesh

  call mesh_merge ( mesh, mesh1, mesh_surface, warn=.false. )
  call fill_mesh_parts ( mesh_surface )
  call delete ( mesh1)

  if ( restart == 0 ) then

!   plot mesh + initial surface

    write(filename,'(a,i4.4,a)') 'mesh_swell_', 0, '.fig'

    call plot_mesh ( plot_options, mesh_surface, filename=filename, groups=[1] )
    plot_options%meshcolor = 4
    call plot_mesh ( plot_options, mesh_surface, filename=filename, &
      groups=[2], append=.true. )
    plot_options%meshcolor = 0

  end if


! create input_probdef

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


! set pointers in module functions to data

  lmesh_sf_adv => mesh_sf_adv
  lproblem_sf_adv => problem_sf_adv
  lsol_sf_adv_pred => sol_sf_adv_pred


! fill solution vector with essential boundary conditions

  sol_sf_adv%u = h0-eps_h0
  call fill_sysvector ( mesh_sf_adv, problem_sf_adv, sol_sf_adv, point=1, &
    physq=1, value=h0 )

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
  call create ( mesh_sf_adv, problem_sf_adv, hgt_end, points=[2] )


! open files and restart

  call open_files_and_restart


! time stepping

  coefficients_sf_adv%i(5) = 1 ! start with first-order scheme

  do step = 1, numtimesteps

    if ( step >= 2 .or. restart == 1 ) then

      coefficients_sf_adv%i(5) = 2  ! second-order scheme

!     predict position of the surface
      sol_sf_adv_pred%u = 2._dp*sol_sf_adv_n%u - sol_sf_adv_nm1%u

    end if

!   eltree for elements crossing the interface

    d = levelset ( mesh%coor ) ! set levelset for all nodes

    call create_eltree_in_elements ( numsplit )

!   nodeset for Dirichlet on outside nodes

    call find_outside_zero_nodes

!   re-define problem

    call delete ( problem )
    call problem_definition ( input_probdef, mesh, problem )

!   create system matrix of velocity/pressure problem

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

    call create_sysmatrix_data ( sysmatrix )

!   build (assemble) matrix/vector

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, oldvectors=oldvectors, &
      coefficients=coefficients )

!   open boundary

    call build_system( mesh, problem, sysmatrix, rhsd, &
      elemsub1=stokes_open_boundary, coefficients=coefficients, &
      object=obj_ob, addmatvec=.true., oldvectors=oldvectors )

!   imposed flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients, oldvectors=oldvectors )

!   surface tension

    if ( surface_tension ) then

!     elementset for elements crossed by the surface

      if ( step == 1 ) then
        call create_elementset_from_eltree
        els_srf = mesh%nelementsets
      else
        call create_elementset_from_eltree ( replace=els_srf )
      end if

!     build elements in the surface

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=surface_tension_eltree, &
        physqrow=[physqvel], physqcol=[physqvel], &
        addmatvec=.true., buildmatrix=.false., coefficients=coefficients, &
        oldvectors=oldvectors, elementset=els_srf )

!     boundary point x1

      nnsf = mesh_sf_adv%nnodes
      nsf = sol_sf_adv_pred%n
      mesh%objects(obj_bp)%coor(1,1) = mesh_sf_adv%coor(nnsf,1)
      mesh%objects(obj_bp)%coor(1,2) = sol_sf_adv_pred%u(hgt%s(nsf))
      call find_refcoor_objects ( mesh, object1=obj_bp )

      if ( mesh%objects(obj_bp)%grpelm(1,1) == 0 ) then
        print*, 'Error: no reference coordinates in boundary point.'
        stop
      end if

!     second point x2

      coefficients%r(501) = mesh_sf_adv%coor(nnsf-1,1)
      coefficients%r(502) = sol_sf_adv_pred%u(hgt%s(nsf-1))

!     build external force in boundary point

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub1=surface_tension_boundary_point2, &
        object=obj_bp, onobjectnodes=.true., &
        physqrow=[physqvel], physqcol=[physqvel], &
        addmatvec=.true., buildmatrix=.false., coefficients=coefficients, &
        oldvectors=oldvectors )

    end if

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve velocity/pressure problem

    solver_options_u%real_storage=rs_up
    solver_options_u%integer_storage=is_up

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

    call delete ( sysmatrix )


!   solve surface advection (corrector)

    call solve_surface_height_corrector


    if ( mod(step0+step,plot_mesh_every) == 0 ) then

!     plot mesh + surface

      write(filename,'(a,i4.4,a)') 'mesh_swell_', step0+step, '.fig'

      n2 = mesh_surface%nnodes
      n1 = n2 - mesh%objects(obj_sh)%nnodes + 1
      mesh_surface%coor(n1:n2,:) = mesh%objects(obj_sh)%coor
      call plot_mesh ( plot_options, mesh_surface, filename=filename, &
        groups=[1] )
      plot_options%meshcolor = 4  ! red
      call plot_mesh ( plot_options, mesh_surface, filename=filename, &
        groups=[2], append=.true. )
      plot_options%meshcolor = 0  ! reset to black

    end if

    call delete_eltree_in_elements

  end do

  close(unit=11)


! create data structures for post-processing

! create eltree for all elements crossing the interface

  call create_eltree_in_elements ( numsplit_plot )

! post-processing: create mesh_plot for plotting

  call mesh_convert ( mesh, mesh_plot, userelmesh=userelmesh, warn=.false. )

  call fill_mesh_parts ( mesh_plot )

! problem definition for post_processing

  call create_input_probdef ( mesh_plot, input_probdef_plot, nvec=2 )

  do i = 1, mesh_plot%nelgrp
    input_probdef_plot%elementdof(i)%a(:) = 0
    input_probdef_plot%vec_elementdof(i)%a(:,1) = 2  ! velocity
    input_probdef_plot%vec_elementdof(i)%a(:,2) = 1  ! scalar
  end do

  call problem_definition ( input_probdef_plot, mesh_plot, problem_plot )


! write data for post-processing and restart

  call write_mesh ( mesh, filename='mesh.out' )
  call write_mesh ( mesh_plot, filename='mesh_plot.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh_plot, input_probdef_plot, &
    filename='probdef_plot.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol%u
  write(10) sol_sf_adv_pred%u, sol_sf_adv_n%u, sol_sf_adv_nm1%u
  write(10) step0 + numtimesteps, time0 + numtimesteps*deltat

  close(unit=10)


! delete all data including all allocated memory

  call delete_eltree_in_elements

  call delete ( mesh, mesh_plot, mesh_surface )
  call delete ( problem, problem_plot )
  call delete ( input_probdef, input_probdef_plot )
  call delete ( sol, rhsd )
  call delete ( oldvectors )
  call delete ( coefficients )

  call delete ( mesh_sf_adv )
  call delete ( problem_sf_adv )
  call delete ( input_probdef_sf_adv )
  call delete ( sol_sf_adv, rhsd_sf_adv )
  call delete ( sol_sf_adv_pred, sol_sf_adv_n, sol_sf_adv_nm1 )
  call delete ( sysmatrix_sf_adv )
  call delete ( oldvectors_sf_adv )
  call delete ( coefficients_sf_adv )
  call delete ( hgt, hgt_end )

  deallocate ( nodes, d, eltree_array, vol )


contains


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

!   compute the maximum, minimum and end radius

    max_height = maxval ( sol_sf_adv%u )
    end_height = sol_sf_adv%u(hgt_end%s(1))


!   write swell height

    write(11,'(i6,4es16.8)') step0+step, time0+step*deltat, &
                             max_height, end_height
    print '(i6,4es16.8)', step0+step, time0+step*deltat, &
                          max_height, end_height

  end subroutine solve_surface_height_corrector


! create eltree in the elements which are crossed by the interface

  subroutine create_eltree_in_elements ( numsplit )

    integer, intent(in) :: numsplit

    integer :: nod(mesh%element(1)%numnod), elem
    real(dp) :: coor(2,2)
    real(dp) :: minvalvol

!   loop all elements

    minvalvol = 1

    do elem = 1, mesh%grpnumel(1)

!     nodal points

      nod = mesh%topology(1)%a(:,elem)

      if ( all ( d(nod) > split_threshold ) .or. &
                        all ( d(nod) < -split_threshold ) ) then

!       all nodes far from the interface

        cycle

      else

!       element possibly contains an interface, start subdivide

        allocate ( eltree_array(1,elem)%p )

!       fill root node of the eltree

        coor(1,:) = -1._dp  ! lower corner of the reference region (-1,-1,-1)
        coor(2,:) =  1._dp  ! upper corner of the reference region (1,1,1)

        call fill_node_eltree ( eltree_array(1,elem)%p, coor )

!       divide root reference domain into subdomains

        lelgrp = 1
        lelem = elem

        call subdivide ( eltree_array(1,elem)%p, levelset=levelset, &
          numsplit=numsplit, split_threshold=split_threshold, &
          numsplitmin=numsplitmin, mapcoor=mapcoor, submesh=.true., &
          intmesh=.true. )

!       check whether interface is in element

        if ( any( number_of_subelements_vector(eltree_array(1,elem)%p,&
                         &lsign=[-1,1]) == 0 ) ) then

!         no interface, remove eltree

          call delete(eltree_array(1,elem)%p)  ! delete actual eltree
          deallocate( eltree_array(1,elem)%p ) ! delete pointer target
                                               ! (pointer becomes disassociated)

        end if

      end if

      if ( associated(eltree_array(1,elem)%p) ) then

!       element crosses the interface

        vol(elem) = volume ( eltree_array(1,elem)%p, lsign=[1] ) / 4

        if ( vol(elem) <= epsvol ) then
!         remove eltree
          call delete(eltree_array(1,elem)%p)  ! delete actual eltree
          deallocate( eltree_array(1,elem)%p ) ! delete pointer target
                                               ! (pointer becomes disassociated)
          print *, 'element at interface removed, elem = ', elem
        end if

        minvalvol = min ( vol(elem), minvalvol )

      end if

    end do

    !print *, 'minvalvol = ', minvalvol

  end subroutine create_eltree_in_elements


! create elementset of elements which are crossed by the interface

  subroutine create_elementset_from_eltree ( replace )

    integer, intent(in), optional :: replace

    integer :: elements(mesh%grpnumel(1)), numelem, elem


!   loop all elements

    numelem = 0

    do elem = 1, mesh%grpnumel(1)

      if ( associated(eltree_array(1,elem)%p) ) then

!       element crosses the interface

        numelem = numelem + 1

        elements(numelem) = elem

      end if

    end do

!   add elementset

    call add_to_mesh ( mesh, elementset='elements', &
                       elements=elements(1:numelem), replace=replace )

  end subroutine create_elementset_from_eltree


! delete eltree in the elements

  subroutine delete_eltree_in_elements

    integer :: elem

!   loop all elements

    do elem = 1, size(eltree_array,2)

      if ( associated(eltree_array(1,elem)%p) ) then

        call delete(eltree_array(1,elem)%p) ! delete actual eltree
        deallocate(eltree_array(1,elem)%p) ! delete pointer target
                                           !(pointer becomes disassociated)

      end if

    end do

  end subroutine delete_eltree_in_elements


! find outside nodes that are not connected to elements at the interface

  subroutine find_outside_zero_nodes

    integer :: elem, i

!   set outside nodes

    where ( d < 0._dp )  ! use < 0 (and not <= 0) to avoid wall nodes to
      nodes = 1          ! end up in the nodeset
    else where
      nodes = 0
    end where

!   remove nodes connected to elements at the interface

    do elem = 1, mesh%grpnumel(1)
      if ( associated(eltree_array(1,elem)%p) ) then
         nodes(mesh%topology(1)%a(:,elem)) = 0
       end if
    end do

!   count and collect the nodes (overwrite array nodes along the way)

    nnodes = 0
    do i = 1, mesh%nnodes
      if ( nodes(i) == 0 ) cycle
      nnodes = nnodes + 1
      nodes(nnodes) = i
    end do

    call add_to_mesh ( mesh, nodeset='nodes', nodes=nodes(1:nnodes), replace=1 )

  end subroutine find_outside_zero_nodes


  subroutine open_files_and_restart

    integer :: endstep
    real(dp) :: endtime

    if ( restart >= 1 ) then

!     restart: read solution from file

      open ( unit=10, form='unformatted', file='data.out' )

      read(10) sol%u
      read(10) sol_sf_adv_pred%u, sol_sf_adv_n%u, sol_sf_adv_nm1%u
      read(10) endstep, endtime

      close(unit=10)

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

end program extrudate_swell4
