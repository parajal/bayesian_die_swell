! 3D viscous extrudate swell problem
! Square die
! second-order Gear time integration

program extrudate_swell7

  use tfem_m
  use kind_defs_m
  use mesh_m
  use meshgen_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m
  use surface_advection_elements_m
  use subs_extrudate_swell_m
  use update_mesh_nodes_surfaces_m
  use figplot_m
  use timer_m
  use projection_elements_m
  use limits_m

  implicit none


! constants flow problem

  integer, parameter :: &
    uintpl = 6,       & ! P2 velocities
    pintpl = 2,       & ! P1 pressures
    ointpl = 6,       & ! P2 shape of object elements
    hintpl = 6,       & ! Interpolation for the height
    hintpl_mat = 6,   & ! Interpolation for the material lines
    physqvel = 1,     & ! physical quantity nr of the velocities
    physqpress = 2,   & ! physical quantity nr of the pressures
    gauss = 5,        & ! 6-point Gauss integration of tets
    gaussb = 5,       & ! 3-point integration of boundary elements
    coorsys = 0         ! planar Cartesian coordinate system

! constants surface advection

  integer, parameter :: &
    ninti_mat = 6,      & ! number of Gauss points
    ninti_sf_adv = 6,   & ! number of Gauss points
    method = 1            ! discretization method 0: Galerkin, 1: SUPG


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, sol_n, sol_nm1
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(lu_ma41_t) :: luc
  type(solver_options_ma41_t) :: solver_options_u

! 1D height function for surface advection
  type(meshgen_options_t) :: mesh_options
  type(plot_options_t) :: plot_options
  type(subscript_t) :: hgt, hgt_end
  type(subscript_t) :: subsh1, subsh2, subsh3, subsh4
  type(subscript_t) :: velx_inlet, vely_inlet, velz_inlet
  type(subscript_t) ::  velxline1, velyline1, velzline1, &
                        velxline2, velyline2, velzline2, velxline3, &
                        velyline3, velzline3, velxline4, velyline4, &
                        velzline4
  type(subscriptvec_t) :: velx_line1, vely_line1, velz_line1, &
                          velx_line2, vely_line2, velz_line2
  type(subscriptvec_t) :: velx_line3, vely_line3, velz_line3, velx_line4, &
                          vely_line4, velz_line4

! 1D material lines
  type(mesh_t) :: mesh_line1, mesh_line2, mesh_line3, mesh_line4
  type(input_probdef_t) :: input_probdef_line1, input_probdef_line2
  type(input_probdef_t) :: input_probdef_line3, input_probdef_line4
  type(problem_t) :: problem_line1, problem_line2, problem_line3, problem_line4
  type(sysmatrix_t) :: sysmatrix_line1, sysmatrix_line2, sysmatrix_line3, &
                       sysmatrix_line4
  type(sysvector_t), dimension(2), target :: sol_line1, sol_line1_n, &
                                             sol_line1_nm1
  type(sysvector_t), dimension(2), target :: sol_line2, sol_line2_n, &
                                             sol_line2_nm1
  type(sysvector_t), dimension(2), target :: sol_line3, sol_line3_n, &
                                             sol_line3_nm1
  type(sysvector_t), dimension(2), target :: sol_line4, sol_line4_n, &
                                             sol_line4_nm1
  type(sysvector_t), dimension(2) :: rhsd_line1, rhsd_line2, rhsd_line3, &
                                     rhsd_line4
  type(oldvectors_t) :: oldvectors_line1, oldvectors_line2, &
                        oldvectors_line3, oldvectors_line4
  type(coefficients_t) ::coefficients_mat
  type(vector_t), target :: velocity_line1, velocity_line2, velocity_line3, &
                            velocity_line4

! 2D height function for surface advection and subscripts of lines and surfs
  type(mesh_t) :: mesh_surf1, mesh_surf2, mesh_surf3, mesh_surf4
  type(input_probdef_t) :: input_probdef_surf1, input_probdef_surf2, &
                           input_probdef_surf3, input_probdef_surf4
  type(problem_t) :: problem_surf1, problem_surf2, problem_surf3, problem_surf4
  type(sysmatrix_t) :: sysmatrix_surf1, sysmatrix_surf2, sysmatrix_surf3, &
                       sysmatrix_surf4
  type(sysvector_t), target :: sol_surf1, sol_surf1_n, sol_surf1_nm1
  type(sysvector_t), target :: sol_surf2, sol_surf2_n, sol_surf2_nm1
  type(sysvector_t), target :: sol_surf3, sol_surf3_n, sol_surf3_nm1
  type(sysvector_t), target :: sol_surf4, sol_surf4_n, sol_surf4_nm1
  type(sysvector_t) :: rhsd_surf1, rhsd_surf3
  type(sysvector_t) :: rhsd_surf2, rhsd_surf4
  type(oldvectors_t) :: oldvectors_surf1, oldvectors_surf2, oldvectors_surf3, &
                        oldvectors_surf4
  type(coefficients_t) :: coefficients_sf_adv
  type(vector_t), target :: velocity_surf1, velocity_surf2, velocity_surf3, &
                            velocity_surf4
  type(subscript_t) ::  velxsurf1, velysurf1, velzsurf1, &
                        velxsurf2, velysurf2, velzsurf2, &
                        velxsurf3, velysurf3, velzsurf3, &
                        velxsurf4, velysurf4, velzsurf4
  type(subscriptvec_t) :: velx_surf1, vely_surf1, velz_surf1, &
                       velx_surf2, vely_surf2, velz_surf2, &
                       velx_surf3, vely_surf3, velz_surf3, &
                       velx_surf4, vely_surf4, velz_surf4
  type(subscriptvec_t) :: meshvelx_surf1, meshvely_surf1, meshvelx_surf2, &
                          meshvely_surf2, meshvelx_surf3, meshvely_surf3, &
                          meshvelx_surf4, meshvely_surf4
  type(subscript_t) :: subs_surf1, subs_surf2, subs_surf3, subs_surf4
  type(subscript_t) :: subs_curve1, subs_curve2, subs_curve3, subs_curve4
   type(vector_t), target :: meshvel, meshvel_surf1, meshvel_surf2, &
                             meshvel_surf3, meshvel_surf4
  type(problem_t) :: problem_lapl, problem_laply
  type(problem_t) :: problem_lapl_surf, problem_laply_surf, problem_laplz_surf
  real(dp), allocatable, dimension(:,:) :: meshcoor_initial


! variables

  integer :: &
    timeint1 = 1,        & ! (first-order) Euler time integration (first step)
    numtimesteps = 1500, & ! number of time steps
    step0 = 0,           & ! initial step number
    vtkevery = 1           ! write vtk every .. steps

  real(dp) :: &
    eta_s = 1.0_dp,       & ! solvent viscosity
    deltat = 1.e-2_dp,    & ! time step
    time0 = 0._dp,        & ! initial time
    U_avg = 0.5_dp,       & ! average velocity at the entry
    beta_sf_adv = 0.5_dp, & ! upwinding parameter in the SUPG method (surface)
    beta_line = 0.5_dp,   & ! upwinding parameter in the SUPG method (lines)
    Ho = 0.5_dp,          & ! Initial height domain
    dx_box = 0.4_dp,      & ! element spacing on the external boundaries
    dx_wall = 0.1_dp,    & ! element spacing on upper boundary and die exit
    rs_up = 3.5_dp,      & ! real_storage gradient-velocity-pressure LU (HSL)
    is_up = 3.5_dp         ! integer_storage gradient-velocity-pressure LU


  integer :: step=1, ipost=0
  integer :: ipost_in=0
  integer :: curve2, nodenr2
  real(dp) :: time

  real(dp), allocatable, dimension(:) :: initial_line1y, initial_line2y, &
                                         initial_line3y, initial_line4y, &
                                         initial_line1z, initial_line2z, &
                                         initial_line3z, initial_line4z
  real(dp), allocatable, dimension(:) :: initial_hsurf1, initial_hsurf2, &
    initial_hsurf3, initial_hsurf4, initial_hsurf1_line,  &
    initial_hsurf2_line, initial_hsurf3_line, initial_hsurf4_line
  real(dp), allocatable, dimension(:,:) :: bd
  real(dp) :: ox=3.0_dp, L2=3.0_dp

!  real(dp), dimension(:,:), allocatable :: coor_boun
  real(dp), allocatable, dimension(:,:) :: meshcoor_n, meshcoor_nm1
  real(dp), allocatable, dimension(:,:) :: meshcoorsurf1_n, &
            meshcoorsurf1_nm1, meshcoorsurf2_n, meshcoorsurf2_nm1, &
            meshcoorsurf3_n, meshcoorsurf3_nm1, meshcoorsurf4_n, &
            meshcoorsurf4_nm1
  real(dp), allocatable, dimension(:) :: displine1, displine2, displine3, &
                                         displine4
  real(dp), allocatable, dimension(:,:) :: surf1hat, surf1hatn, surf2hat, &
                                         surf2hatn, surf3hat, surf3hatn, &
                                         surf4hat, surf4hatn

  character(len=20) :: filename

! All parameters for inlet problem

! definitions

  type(mesh_t) :: mesh_inlet
  type(input_probdef_t) :: input_probdef_inlet
  type(problem_t), target :: problem_inlet
  type(sysmatrix_t) :: sysmatrix_inlet
  type(sysvector_t), target :: sol_inlet, solm1_inlet
  type(sysvector_t) :: rhsd_inlet
  type(oldvectors_t) :: oldvectors_inlet
  type(coefficients_t) :: coefficients_inlet
  type(subscript_t) :: velx_in, vely_in, velz_in
  type(solver_options_ma41_t) :: solver_options_u_in


  ! variables

  real(dp) :: &
    rs_up_in = 2.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_up_in = 2.0_dp      ! integer_storage gradient-velocity-pressure LU

  real(dp) :: flowrate_in

  SEARCH_BLOCKS_PASSTWO = .true.
  EPSELEMENT2 = 3e-2

! timer

  timer = .false.

  flowrate_in = U_avg*Ho**2

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=600, ncoefr=100 )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,     0, &
      physqvel, physqpress, 0,     0, gauss, &
      gaussb ]
  coefficients%i(23) = 0

  coefficients%r = 0
  coefficients%r(1) = eta_s

  coefficients%i(35) = ointpl
  coefficients%i(40) = 3  ! Standard Gauss-Legendre (numerical)

 call generate_read_mesh

! fill coefficients inlet

  call create_coefficients ( coefficients_inlet, ncoefi=150, ncoefr=100 )

  coefficients_inlet%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,     0,  &
       physqvel, physqpress, 0,    0 , gaussb,     &
       gaussb ]


  coefficients_inlet%i(23) = coorsys  ! coorsys = 1, axisymmetric
  coefficients_inlet%i(40) = 3  ! Standard Gauss-Legendre (numerical)
  coefficients_inlet%i(67) = 1  ! 3D velocity (swirl)

  coefficients_inlet%r = 0
  coefficients_inlet%r(1) = eta_s
  coefficients_inlet%r(6) = flowrate_in

  call define_inlet_problem


  call define_problems_create_vectors

! fill coefficients for surface advection problem
  call create_coefficients ( coefficients_mat, ncoefi=100, ncoefr=50 )

  coefficients_mat%i = 0
  coefficients_mat%i(1) = ninti_mat
  coefficients_mat%i(2) = 2 ! velocity given by nodal values
  coefficients_mat%i(4) = method ! 0: Galerkin, 1: SUPG
  coefficients_mat%i(5) = 1 ! time-integration
  coefficients_mat%i(6) = hintpl_mat ! height interpolation
  coefficients_mat%i(9) = 3 ! numerical table for Gauss

  coefficients_mat%r = 0
  coefficients_mat%r(4) = deltat
  coefficients_mat%r(5) = beta_line

! fill coefficients for surface advection problem

  call create_coefficients ( coefficients_sf_adv, ncoefi=100, ncoefr=50 )

  coefficients_sf_adv%i = 0
  coefficients_sf_adv%i(1) = ninti_sf_adv
  coefficients_sf_adv%i(2) = 2 ! velocity given by nodal values
  coefficients_sf_adv%i(4) = method ! 0: Galerkin, 1: SUPG
  coefficients_sf_adv%i(5) = 1 ! time-integration
  coefficients_sf_adv%i(6) = hintpl ! height interpolation
  coefficients_sf_adv%i(9) = 3 ! numerical table for Gauss

  coefficients_sf_adv%r = 0
  coefficients_sf_adv%r(4) = deltat
  coefficients_sf_adv%r(5) = beta_sf_adv

  call define_line_problems
  call define_surface_problems

  call tic

! time stepping

  time = 0

  coefficients%i(22) = timeint1  ! first step integration scheme

  do step = 1, numtimesteps

    time = time + deltat

    call solve_inlet_problem

    sol%u(velx_inlet%s) = sol_inlet%u(velz_in%s)
    sol%u(vely_inlet%s) = sol_inlet%u(velx_in%s)
    sol%u(velz_inlet%s) = sol_inlet%u(vely_in%s)

    if ( step >= 2 ) then

      coefficients_mat%i(5) = 2 ! second order scheme
      coefficients_sf_adv%i(5) = 2  ! second-order scheme

      meshcoor_nm1 = meshcoor_n
      meshcoor_n = mesh%coor

      surf1hatn(:,1) = mesh%coor(mesh%surfaces(8)%nodes(:),2)
      surf1hatn(:,2) = mesh%coor(mesh%surfaces(8)%nodes(:),3)
      surf1hat(:,1) = 2._dp*meshcoor_n(mesh%surfaces(8)%nodes(:),2)- &
                        meshcoor_nm1(mesh%surfaces(8)%nodes(:),2)
      surf1hat(:,2) = 2._dp*meshcoor_n(mesh%surfaces(8)%nodes(:),3)- &
                        meshcoor_nm1(mesh%surfaces(8)%nodes(:),3)

      surf2hatn(:,1) = mesh%coor(mesh%surfaces(9)%nodes(:),3)
      surf2hatn(:,2) = mesh%coor(mesh%surfaces(9)%nodes(:),2)
      surf2hat(:,1) = 2._dp*meshcoor_n(mesh%surfaces(9)%nodes(:),3)- &
                     meshcoor_nm1(mesh%surfaces(9)%nodes(:),3)
      surf2hat(:,2) = 2._dp*meshcoor_n(mesh%surfaces(9)%nodes(:),2)- &
                     meshcoor_nm1(mesh%surfaces(9)%nodes(:),2)

      surf3hatn(:,1) = mesh%coor(mesh%surfaces(10)%nodes(:),2)
      surf3hatn(:,2) = mesh%coor(mesh%surfaces(10)%nodes(:),3)
      surf3hat(:,1) = 2._dp*meshcoor_n(mesh%surfaces(10)%nodes(:),2)- &
                        meshcoor_nm1(mesh%surfaces(10)%nodes(:),2)
      surf3hat(:,2) = 2._dp*meshcoor_n(mesh%surfaces(10)%nodes(:),3)- &
                        meshcoor_nm1(mesh%surfaces(10)%nodes(:),3)

      surf4hatn(:,1) = mesh%coor(mesh%surfaces(11)%nodes(:),3)
      surf4hatn(:,2) = mesh%coor(mesh%surfaces(11)%nodes(:),2)
      surf4hat(:,1) = 2._dp*meshcoor_n(mesh%surfaces(11)%nodes(:),3)- &
                     meshcoor_nm1(mesh%surfaces(11)%nodes(:),3)
      surf4hat(:,2) = 2._dp*meshcoor_n(mesh%surfaces(11)%nodes(:),2)- &
                     meshcoor_nm1(mesh%surfaces(11)%nodes(:),2)

      call update_mesh_nodes_surfaces ( mesh, problem_lapl_surf, &
           problem_laply_surf, problem_laplz_surf, disp1=surf1hat-surf1hatn, &
           disp2=surf2hat-surf2hatn, disp3=surf3hat-surf3hatn, &
           disp4=surf4hat-surf4hatn )

      call find_bounds_blocks ( mesh )

      call write_mesh_vtk ( mesh, 'meshmove.vtk')

      call toc ( 'update mesh' )

      meshvel%u = reshape ( transpose ( &
            ( 1.5_dp*mesh%coor - 2*meshcoor_n + 0.5_dp*meshcoor_nm1 ) /&
              deltat ), [mesh%ndim*mesh%nnodes] )
    end if

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_solve_vpG

!   solve material lines
    call solve_material_line_corrector

    call expand_surfaces

    call solve_4surfaces_height_corrector

    print '(i6,4es16.8)', step0+step, time0+step*deltat

    if (step ==1 ) then
      call postprocessing ( mesh, sol, ipost=ipost )
      ipost = ipost + 1
    end if

!   write VTK files
    if ( vtkevery > 0 .and. step>1 ) then
      if ( mod(step,vtkevery) == 0 ) then
        call postprocessing ( mesh, sol, ipost )
        ipost = ipost + 1
      end if
    end if

    call toc ( 'solve_surface_system' )

  end do

  call delete_old_problems

contains

! generate and read mesh

  subroutine generate_read_mesh

    type (mesh_t) :: mesh1, mesh2

    integer :: i

    call write_gmsh_parameters ( ox, L2, dx_box, dx_wall)

    call execute_command_line ( 'gmsh -3 -order 2 -o mesh.msh &
                  &mesh.geo > outputmesh.out' )

!   read mesh generated by gmsh
    call read_mesh_gmsh ( mesh, filename='mesh.msh', ndim=3, &
      physgeom=.true. )

    call mesh_convert ( mesh, mesh1, only_groups=[1] )
    call mesh_convert ( mesh, mesh2, only_groups=[2] )
    call delete ( mesh )

    call mesh_merge ( mesh1, mesh2, mesh, surface1=2, surface2=1 )

    call add_to_mesh ( mesh, object='surface', objectsurface=1, &
         topology=.true., intrule=3 )

! add volume made from mesh elements for integral p dV =0

    call add_to_mesh ( mesh, volumefromgroups=[1] )

    call fill_mesh_parts ( mesh )

    do i = 1, mesh%ncurves
      write(filename,'(a,i4.4,a)') 'curve_',i,'.vtk'
      call write_geometry_vtk ( mesh, curve=i, filename=filename )
    end do

    do i = 1, mesh%nsurfaces
       write(filename,'(a,i4.4,a)') 'surface_',i,'.vtk'
       call write_geometry_vtk ( mesh, surface=i, filename=filename )
    end do

    do i = 1, mesh%npoints
       write(filename,'(a,i4.4,a)') 'point_',i,'.vtk'
       call write_geometry_vtk ( mesh, point=i, filename=filename )
    end do

    call plot_points_curves ( plot_options, mesh, 'points.fig' )

    call printinfo ( mesh, printlevel=2 )

!   define some arrays for the ALE mesh position at old times
    allocate ( meshcoor_n(mesh%nnodes,mesh%ndim), &
      meshcoor_nm1(mesh%nnodes,mesh%ndim) )

    allocate ( meshcoor_initial(mesh%nnodes,2) )
      meshcoor_initial = mesh%coor

    meshcoor_n = mesh%coor
    meshcoor_nm1 = mesh%coor

    call write_mesh_vtk ( mesh, filename='mesh.vtk' )

  end subroutine generate_read_mesh

! write the parameters in gmsh format

  subroutine write_gmsh_parameters ( ox, L2, dx_box, dx_wall )

    real(dp), intent(in) :: ox, L2, dx_box, dx_wall

     open ( unit=25, file='mesh.geo' )

     write ( 25, '(1X,A,F18.14,A)' ) 'ox = ', ox, ';'
     write ( 25, '(1X,A,F18.14,A)' ) 'L2 = ', L2, ';'
     write ( 25, '(1X,A,F18.14,A)' ) 'dx_box = ', dx_box, ';'
     write ( 25, '(1X,A,F18.14,A)' ) 'Ho = ', Ho, ';'

     write ( 25, '(1X,A,F18.14,A)' ) 'dx_wall = ', dx_wall, ';'

     write ( 25, '(/1x,a)' ) 'Include "mesh3D.igo";'

     close ( 25 )

  end subroutine write_gmsh_parameters

  subroutine define_problems_create_vectors

! problem definition of gradient/velocity/pressure

!   problem definition of gradient/velocity/pressure
    call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

    input_probdef%vec_elementdof(1)%a =   &
        reshape ( [ 3,3,3,3,3,3,3,3,3,3,    &  ! velocity
                     1,0,1,0,1,0,0,0,0,1,    &  ! pressure
                     1,1,1,1,1,1,1,1,1,1,    &  ! scalar, such as vorticity
                     1,1,1,1,1,1,1,1,1,1 ], &  ! effective shear rate
                     [10,4] )

    input_probdef%physq = [1,2]
    input_probdef%probnr = 1


! Dirichlet boundary conditions

! outflow
   call define_essential ( mesh, input_probdef, &
      surface1=7, physq=physqvel, degfd=[0,1,1] )

! inflow
    call define_essential ( mesh, input_probdef, &
      surface1=1, physq=physqvel )

! wall
    call define_essential ( mesh, input_probdef, &
      surface1=5, physq=physqvel, excludesurfaces=[1] )
    call define_essential ( mesh, input_probdef, &
      surface1=6, physq=physqvel, excludesurfaces=[1] )
    call define_essential ( mesh, input_probdef, &
      surface1=3, physq=physqvel, excludesurfaces=[1])
    call define_essential ( mesh, input_probdef, &
      surface1=4, physq=physqvel, excludesurfaces=[1])

    call problem_definition ( input_probdef, mesh, problem )

! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

    call create_sysvector ( problem, sol, sol_n, sol_nm1, rhsd )

! fill solution vector with essential boundary conditions

    sol%u = 0

! create subscripts for the meshvelocity

    call create_subscript ( mesh, problem, velx_inlet, physqarr=[physqvel], &
      degfd=1, surfaces=[1], fillnodes=.true. )
    call create_subscript ( mesh, problem, vely_inlet, physqarr=[physqvel], &
      degfd=2, surfaces=[1], fillnodes=.true. )
    call create_subscript ( mesh, problem, velz_inlet, physqarr=[physqvel], &
      degfd=3, surfaces=[1], fillnodes=.true. )

    call create_subscript ( mesh, problem, velxline1, physqarr=[physqvel], &
      degfd=1, curves=[19], fillnodes=.true. )
    call create_subscript ( mesh, problem, velyline1, physqarr=[physqvel], &
      degfd=2, curves=[19], fillnodes=.true. )
    call create_subscript ( mesh, problem, velzline1, physqarr=[physqvel], &
      degfd=3, curves=[19], fillnodes=.true. )

    call create_subscript ( mesh, problem, velxline2, physqarr=[physqvel], &
      degfd=1, curves=[20], fillnodes=.true. )
    call create_subscript ( mesh, problem, velyline2, physqarr=[physqvel], &
      degfd=2, curves=[20], fillnodes=.true. )
    call create_subscript ( mesh, problem, velzline2, physqarr=[physqvel], &
      degfd=3, curves=[20], fillnodes=.true. )

    call create_subscript ( mesh, problem, velxline3, physqarr=[physqvel], &
      degfd=1, curves=[17], fillnodes=.true. )
    call create_subscript ( mesh, problem, velyline3, physqarr=[physqvel], &
      degfd=2, curves=[17], fillnodes=.true. )
    call create_subscript ( mesh, problem, velzline3, physqarr=[physqvel], &
      degfd=3, curves=[17], fillnodes=.true. )

    call create_subscript ( mesh, problem, velxline4, physqarr=[physqvel], &
      degfd=1, curves=[18], fillnodes=.true. )
    call create_subscript ( mesh, problem, velyline4, physqarr=[physqvel], &
      degfd=2, curves=[18], fillnodes=.true. )
    call create_subscript ( mesh, problem, velzline4, physqarr=[physqvel], &
      degfd=3, curves=[18], fillnodes=.true. )

    call create_subscript ( mesh, problem, velxsurf1, physqarr=[physqvel], &
      degfd=1, surfaces=[8], fillnodes=.true. )
    call create_subscript ( mesh, problem, velysurf1, physqarr=[physqvel], &
      degfd=2, surfaces=[8], fillnodes=.true. )
    call create_subscript ( mesh, problem, velzsurf1, physqarr=[physqvel], &
      degfd=3, surfaces=[8], fillnodes=.true. )

    call create_subscript ( mesh, problem, velxsurf2, physqarr=[physqvel], &
      degfd=1, surfaces=[9], fillnodes=.true. )
    call create_subscript ( mesh, problem, velysurf2, physqarr=[physqvel], &
      degfd=2, surfaces=[9], fillnodes=.true. )
    call create_subscript ( mesh, problem, velzsurf2, physqarr=[physqvel], &
      degfd=3, surfaces=[9], fillnodes=.true. )

    call create_subscript ( mesh, problem, velxsurf3, physqarr=[physqvel], &
      degfd=1, surfaces=[10], fillnodes=.true. )
    call create_subscript ( mesh, problem, velysurf3, physqarr=[physqvel], &
      degfd=2, surfaces=[10], fillnodes=.true. )
    call create_subscript ( mesh, problem, velzsurf3, physqarr=[physqvel], &
      degfd=3, surfaces=[10], fillnodes=.true. )

    call create_subscript ( mesh, problem, velxsurf4, physqarr=[physqvel], &
      degfd=1, surfaces=[11], fillnodes=.true. )
    call create_subscript ( mesh, problem, velysurf4, physqarr=[physqvel], &
      degfd=2, surfaces=[11], fillnodes=.true. )
    call create_subscript ( mesh, problem, velzsurf4, physqarr=[physqvel], &
      degfd=3, surfaces=[11], fillnodes=.true. )

! create system matrix of gradient/velocity/pressure problem

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

    call create_sysmatrix_data ( sysmatrix )

    call create_vector ( problem, meshvel, physq=physqvel )

!   initialize all vectors to zero
    sol_nm1%u = 0.0_dp
    sol_n%u = 0.0_dp
    meshvel%u = 0._dp

! create the structure oldvectors

    call create_oldvectors ( oldvectors, nsysvec=2, nprob=1, &
      nvec=2 )

! store solution vectors and problem structures

    oldvectors%s(1)%p => sol
    oldvectors%s(2)%p => sol_n
    oldvectors%p(1)%p => problem
    oldvectors%v(1)%p => meshvel

  end subroutine define_problems_create_vectors

  subroutine build_solve_vpG

!   copy old solution
    call copy ( sol, sol_n )
    call copy ( sol_n, sol_nm1 )

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, oldvectors=oldvectors, coefficients=coefficients )

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    call toc ( 'build_flow_system' )

!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_up
    solver_options_u%integer_storage=is_up

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

    call toc ( 'solve_flow_system' )

  end subroutine build_solve_vpG

! write the data to .vtk files

  subroutine postprocessing ( mesh, sol, ipost  )

    type(mesh_t), intent(inout) :: mesh
    type(sysvector_t), intent(in), target :: sol
    integer, intent(in), optional :: ipost

    type(oldvectors_t) :: oldvectors_d
    type(vector_t) :: pressure

    if ( .not. mesh%meshparts) call fill_mesh_parts ( mesh )

    write (filename,'(a,i4.4,a)')'meshtime', ipost, '.vtk'
    call write_mesh_vtk ( mesh, filename=filename )

    call create_oldvectors ( oldvectors_d, nsysvec=1, nsysvec2=1 )
    oldvectors_d%s(1)%p => sol

    call create_vector ( problem, pressure, vec=3 )

!   derive the pressure in all nodes
    call derive_vector ( mesh, problem, pressure, &
      elemsub=stokes_pressure, coefficients=coefficients, &
      oldvectors=oldvectors_d )

    if ( present(ipost) ) write(filename,'(a,i4.4,a)') 'flow', ipost, '.vtk'
    call write_scalar_vtk ( mesh, problem, vector=pressure, &
      dataname='pressure',  filename=filename )

    call write_vector_vtk ( mesh, problem, filename=filename, &
      dataname='velocity', sysvector=sol, physq=physqvel, &
      append=.true. )

    call delete(pressure)

    call delete(oldvectors_d)

  end subroutine postprocessing

  subroutine delete_old_problems


!   delete all data including all allocated memory

    deallocate ( meshcoor_initial )

    call delete ( problem_lapl_surf, problem_laply_surf, problem_laplz_surf )
    call delete ( sysmatrix)
    call delete ( input_probdef)
    call delete ( mesh )
    call delete ( sol, sol_n, sol_nm1 )
    call delete ( rhsd )
    call delete ( oldvectors )
    call delete ( meshvel )

    deallocate ( meshcoor_n, meshcoor_nm1 )

    ! delete all inlet data including all allocated memory

     call delete ( coefficients_mat)
     call delete ( coefficients_sf_adv )
     call delete ( problem_surf1, problem_surf2, problem_surf3, problem_surf4)
     call delete ( sol_surf1, sol_surf2, sol_surf3, sol_surf4 )
     call delete ( sysmatrix_surf1, sysmatrix_surf2, sysmatrix_surf3, &
                   sysmatrix_surf4 )
     call delete ( input_probdef_surf1, input_probdef_surf2, &
                   input_probdef_surf3, input_probdef_surf4)
     call delete ( problem_line1, problem_line2, problem_line3, problem_line4)
     call delete ( sol_line1, sol_line2, sol_line3, sol_line4 )
     call delete ( sysmatrix_line1, sysmatrix_line2, sysmatrix_line3, &
                   sysmatrix_line4)
     call delete ( input_probdef_line1, input_probdef_line2, &
                   input_probdef_line3, input_probdef_line4)

     call delete ( mesh_inlet )
     call delete ( problem_inlet )
     call delete ( input_probdef_inlet )
     call delete ( sol_inlet, rhsd_inlet )
     call delete ( sysmatrix_inlet )
     call delete ( oldvectors_inlet )

     call delete ( coefficients_inlet )
     call delete ( velx_in, vely_in, velz_in )

  end subroutine delete_old_problems

  subroutine define_inlet_problem

!   create mesh

    call generate_inlet_mesh

!   problem definition of gradient/velocity/pressure for the channel

    call create_input_probdef ( mesh_inlet, input_probdef_inlet, nvec=3, &
                                nphysq=2 )

    input_probdef_inlet%vec_elementdof(1)%a =   &
        reshape ( [ 3,3,3,3,3,3,    &  ! velocity
                     1,0,1,0,1,0,    &  ! pressure
                     1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                     [6,3] )

    input_probdef_inlet%physq = [1,2]
    input_probdef_inlet%probnr = 2

!   outer wall
    call define_essential ( mesh_inlet, input_probdef_inlet, curve1=1, &
                            curve2=4, physq=physqvel )
!   pressure level
    call define_essential ( mesh_inlet, input_probdef_inlet, point=1, &
                            physq=physqpress )

!   constraint for flow rate

    call define_constraint ( mesh_inlet, input_probdef_inlet, &
      physq=physqvel, surface1=1, nglobalc=1 )

    call problem_definition ( input_probdef_inlet, mesh_inlet, problem_inlet )

!   create system vectors (solution and right-hand side)

    call create_sysvector ( problem_inlet, sol_inlet, solm1_inlet, rhsd_inlet )


    call create_subscript ( mesh_inlet, problem_inlet, velx_in, &
         physqarr=[physqvel],  degfd=1, surfaces=[1], fillnodes=.true. )
    call create_subscript ( mesh_inlet, problem_inlet, vely_in, &
         physqarr=[physqvel],  degfd=2, surfaces=[1], fillnodes=.true. )
    call create_subscript ( mesh_inlet, problem_inlet, velz_in, &
         physqarr=[physqvel], degfd=3, surfaces=[1], fillnodes=.true. )

!   fill solution vector with essential boundary conditions

    sol_inlet%u = 0
    call fill_sysvector ( mesh_inlet, problem_inlet, sol_inlet, &
      curve1=3, physq=physqvel, degfd=1, value=1.0e-20_dp )
    call fill_sysvector ( mesh_inlet, problem_inlet, sol_inlet, &
      point=1, physq=physqpress, value=0._dp )

!   create system matrix

    call create_sysmatrix_structure_base ( sysmatrix_inlet, mesh_inlet, &
         problem_inlet )
    call create_sysmatrix_structure_constraint ( sysmatrix_inlet, mesh_inlet, &
         problem_inlet )
    call finalize_sysmatrix_structure ( sysmatrix_inlet )

    call create_sysmatrix_data ( sysmatrix_inlet )

!   create the structure oldvectors

    call create_oldvectors ( oldvectors_inlet, nsysvec=2, nprob=1 )

!   store solution vectors and problem structures

    oldvectors_inlet%s(1)%p => sol_inlet
    oldvectors_inlet%s(2)%p => solm1_inlet
    oldvectors_inlet%p(1)%p => problem_inlet

  end subroutine define_inlet_problem

  subroutine solve_inlet_problem

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG_inlet

    call add_effect_of_essential_to_rhs ( problem_inlet, sysmatrix_inlet, &
         sol_inlet, rhsd_inlet )

!   solve gradient/velocity/pressure problem

    solver_options_u_in%real_storage=rs_up_in
    solver_options_u_in%integer_storage=is_up_in

    call solve_system_ma41 ( sysmatrix_inlet, rhsd_inlet, sol_inlet, &
      solver_options=solver_options_u_in  )

    call copy ( sol_inlet, solm1_inlet )

    if (step ==1 ) then
      call postprocessing_inlet ( mesh_inlet, sol_inlet, ipost=ipost_in )
      ipost_in = ipost_in + 1
    end if

!   write VTK files
    if ( vtkevery > 0 ) then
      if ( mod(step,vtkevery) == 0 ) then
         call postprocessing_inlet ( mesh_inlet, sol_inlet, ipost=ipost_in )
         ipost_in = ipost_in + 1
      end if
    end if

  end subroutine solve_inlet_problem

  subroutine generate_inlet_mesh

    integer :: i, curve, node1, node2, node

    call mesh_skeleton ( mesh_inlet, mesh%surfaces(1)%nnodes, &
      mesh%surfaces(1)%nelem, mesh%surfaces(1)%element%elshape, ndim=2 )

    mesh_inlet%coor = mesh%coor(mesh%surfaces(1)%nodes,[2,3])
    mesh_inlet%topology(1)%a = mesh%surfaces(1)%topology(:,:,1)

    allocate(bd(mesh_inlet%ndim,2))

    do i=1,mesh_inlet%ndim
      bd(i,1)=minval(mesh_inlet%coor(:,i))-1.e-10_dp
      bd(i,2)=maxval(mesh_inlet%coor(:,i))+1.e-10_dp
    end do

    call add_to_mesh ( mesh_inlet, blocks=[(5,i=1,mesh_inlet%ndim)], &
         blocksdomain=bd )

    deallocate(bd)

!   add curves to the 2D mesh

    deallocate(mesh_inlet%curves)
    allocate(mesh_inlet%curves(4))
    mesh_inlet%ncurves = 4

    do curve = 1, 4

      if (curve == 1) curve2=1
      if (curve == 2) curve2=2
      if (curve == 3) curve2=3
      if (curve == 4) curve2=4

      mesh_inlet%curves(curve)%nnodes = mesh%curves(curve2)%nnodes
      mesh_inlet%curves(curve)%nelem = mesh%curves(curve2)%nelem
      mesh_inlet%curves(curve)%elnumnod = mesh%curves(curve2)%elnumnod
      allocate (mesh_inlet%curves(curve)%nodes(mesh_inlet%curves(curve)%nnodes))
      allocate (mesh_inlet%curves(curve)%topology( &
          mesh%curves(curve2)%element%numnod, &
          mesh_inlet%curves(curve)%nelem, 2))
      mesh_inlet%curves(curve)%nodes = mesh%curves(curve2)%nodes
      mesh_inlet%curves(curve)%topology = mesh%curves(curve2)%topology
      mesh_inlet%curves(curve)%element%elshape = &
                               mesh%curves(curve2)%element%elshape
      mesh_inlet%curves(curve)%ndim = 2

      do node2 = 1, mesh%curves(curve2)%nnodes

        nodenr2 = mesh%curves(curve2)%nodes(node2)

        do node1 = 1, mesh_inlet%nnodes

          if ( abs(mesh%coor(nodenr2,2)-mesh_inlet%coor(node1,1))<1e-12_dp &
            .and. abs(mesh%coor(nodenr2,3)-mesh_inlet%coor(node1,2))<1e-12_dp )&
             then

!           matching node
            mesh_inlet%curves(curve)%nodes(node2) = node1

            exit ! leave loop

          end if

        end do

      end do

!     update the global node numbers in the toplogy
      do node = 1, mesh_inlet%curves(curve)%nnodes
        do i = 1, mesh_inlet%curves(curve)%nelem
         mesh_inlet%curves(curve)%topology(:,i,2) = &
         mesh_inlet%curves(curve)%nodes(&
                                  mesh_inlet%curves(curve)%topology(:,i,1))
        end do
      end do

      if ( any( mesh_inlet%curves(curve)%nodes == 0 ) ) then
        print *,'Error: no corresponding node found!'
        stop
      end if

      allocate ( mesh_inlet%curves(curve)%element_blend(0) )
      allocate ( mesh_inlet%curves(curve)%nnodes_blend(2) )
      mesh_inlet%curves(curve)%nnodes_blend = [0, &
            mesh_inlet%curves(curve)%nnodes ]

    end do

    WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = .false.

!   make domain into a surface for the flow rate constraint
    call add_to_mesh ( mesh_inlet, surfacefromgroups=[1] )
    call add_to_mesh ( mesh_inlet, point=[Ho,Ho])

    call find_bounds_blocks ( mesh_inlet )

    call fill_mesh_parts ( mesh_inlet )

    call write_mesh_vtk ( mesh_inlet, 'mesh_inlet.vtk' )

    do i = 1, mesh_inlet%ncurves
      write(filename,'(a,i4.4,a)') 'curvein_',i,'.vtk'
      call write_geometry_vtk ( mesh_inlet, curve=i, filename=filename )
    end do

    call write_mesh_vtk ( mesh_inlet, filename='mesh_inlet.vtk' )

  end subroutine generate_inlet_mesh

  subroutine build_vpG_inlet

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh_inlet, problem_inlet, sysmatrix_inlet, &
      rhsd_inlet, elemsub=stokes_elem, oldvectors=oldvectors_inlet, &
      coefficients=coefficients_inlet )

!   imposed flow rate

    call build_system_constraint ( mesh_inlet, problem_inlet, &
      sysmatrix_inlet, rhsd_inlet, constraint1=1, &
      elemsub=stokes_constr_flowr_surface, addmatvec=.true., &
      coefficients=coefficients_inlet )

  end subroutine build_vpG_inlet

  subroutine postprocessing_inlet ( mesh_in, sol_in, ipost  )

    type(vector_t) :: velocity_in, pressure
    type(oldvectors_t) :: oldvectors_in

    type(mesh_t), intent(inout) :: mesh_in
    type(sysvector_t), intent(in), target :: sol_in
    integer, intent(in), optional :: ipost

!   create the structure oldvectors

    call create_oldvectors ( oldvectors_in, nsysvec=1)

    call create_vector ( problem_inlet, velocity_in, physq=physqvel )
    call create_vector ( problem_inlet, pressure, vec=3 )

    call extract_physvector ( mesh_in, problem_inlet, sol_inlet, velocity_in )

    oldvectors_in%s(1)%p => sol_in

    call derive_vector ( mesh_in, problem_inlet, pressure, &
      elemsub=stokes_pressure, coefficients=coefficients_inlet, &
      oldvectors=oldvectors_in )

    coefficients_inlet%i(13)=11

!   write to a vtk file for further post-processing

    if ( present(ipost) ) write(filename,'(a,i4.4,a)') 'inletproblem', &
         ipost, '.vtk'
    call write_scalar_vtk ( mesh_in, problem_inlet, vector=pressure, &
    dataname='pressure', filename=filename )

    if ( present(ipost) ) write(filename,'(a,i4.4,a)') 'inletproblem', &
         ipost, '.vtk'
    call write_vector_vtk ( mesh_in, problem_inlet, filename=filename, &
      dataname='velocity_vector', vector=velocity_in, append=.true., &
      assume3D=.true. )

    call delete(oldvectors_in)
    call delete(velocity_in)
    call delete(pressure)

  end subroutine postprocessing_inlet

  subroutine define_line_problems

    integer :: i

!   create mesh for material lines

    mesh_options%elshape = 2   ! three-node line elements

    mesh_options%nx = mesh%curves(19)%nelem
    call line1d ( mesh_line1, mesh_options )

    mesh_options%nx = mesh%curves(20)%nelem
    call line1d ( mesh_line2, mesh_options )

    mesh_options%nx = mesh%curves(17)%nelem
    call line1d ( mesh_line3, mesh_options )

    mesh_options%nx = mesh%curves(18)%nelem
    call line1d ( mesh_line4, mesh_options )

!   set coordinates
    do i = 1, mesh%curves(19)%nnodes
      mesh_line1%coor(i,1) = mesh%coor(mesh%curves(19)%nodes(i),1)
    end do

    call fill_mesh_parts ( mesh_line1 )

    call write_mesh_vtk ( mesh_line1, filename='mesh_line1.vtk' )

    do i = 1, mesh%curves(20)%nnodes
      mesh_line2%coor(i,1) = mesh%coor(mesh%curves(20)%nodes(i),1)
    end do

    call fill_mesh_parts ( mesh_line2 )
    call write_mesh_vtk ( mesh_line2, filename='mesh_line2.vtk' )

    do i = 1, mesh%curves(17)%nnodes
      mesh_line3%coor(i,1) = mesh%coor(mesh%curves(17)%nodes(i),1)
    end do

    call fill_mesh_parts ( mesh_line3 )
    call write_mesh_vtk ( mesh_line3, filename='mesh_line3.vtk' )

    do i = 1, mesh%curves(18)%nnodes
      mesh_line4%coor(i,1) = mesh%coor(mesh%curves(18)%nodes(i),1)
    end do

    call fill_mesh_parts ( mesh_line4 )
    call write_mesh_vtk ( mesh_line4, filename=' mesh_line4.vtk' )

    allocate(initial_line1y(mesh_line1%nnodes))
    allocate(initial_line2y(mesh_line2%nnodes))
    allocate(initial_line3y(mesh_line3%nnodes))
    allocate(initial_line4y(mesh_line4%nnodes))
    allocate(initial_line1z(mesh_line1%nnodes))
    allocate(initial_line2z(mesh_line2%nnodes))
    allocate(initial_line3z(mesh_line3%nnodes))
    allocate(initial_line4z(mesh_line4%nnodes))

    do i = 1,mesh%curves(19)%nnodes
      initial_line1y(i) = mesh%coor(mesh%curves(19)%nodes(i),2)
      initial_line1z(i) = mesh%coor(mesh%curves(19)%nodes(i),3)
    end do

    do i = 1, mesh%curves(20)%nnodes
      initial_line2y(i) = mesh%coor(mesh%curves(20)%nodes(i),2)
      initial_line2z(i) = mesh%coor(mesh%curves(20)%nodes(i),3)
    end do

    do i = 1, mesh%curves(17)%nnodes
      initial_line3y(i) = mesh%coor(mesh%curves(17)%nodes(i),2)
      initial_line3z(i) = mesh%coor(mesh%curves(17)%nodes(i),3)
    end do

    do i = 1, mesh%curves(18)%nnodes
      initial_line4y(i) = mesh%coor(mesh%curves(18)%nodes(i),2)
      initial_line4z(i) = mesh%coor(mesh%curves(18)%nodes(i),3)
    end do

  ! problem definition for surface advection

    call create_input_probdef ( mesh_line1, input_probdef_line1, nvec=2, &
      nphysq=1 )

    input_probdef_line1%vec_elementdof(1)%a = &
        reshape ( [ 1,1,1,    &  ! height function
                    3,3,3 ],  &  ! velocity
                  [3,2] )

    input_probdef_line1%physq = [1]
    input_probdef_line1%probnr = 6

    call create_input_probdef ( mesh_line2, input_probdef_line2, nvec=2, &
      nphysq=1 )

    input_probdef_line2%vec_elementdof(1)%a = &
        reshape ( [ 1,1,1,    &  ! height function
                    3,3,3 ],  &  ! velocity
                  [3,2] )

    input_probdef_line2%physq = [1]
    input_probdef_line2%probnr = 7

    call create_input_probdef ( mesh_line3, input_probdef_line3, nvec=2, &
      nphysq=1 )

    input_probdef_line3%vec_elementdof(1)%a = &
        reshape ( [ 1,1,1,    &  ! height function
                    3,3,3 ],  &  ! velocity
                  [3,2] )

    input_probdef_line3%physq = [1]
    input_probdef_line3%probnr = 8

    call create_input_probdef ( mesh_line4, input_probdef_line4, nvec=2, &
      nphysq=1 )

    input_probdef_line4%vec_elementdof(1)%a = &
        reshape ( [ 1,1,1,    &  ! height function
                    3,3,3 ],  &  ! velocity
                  [3,2] )

    input_probdef_line4%physq = [1]
    input_probdef_line4%probnr = 9

    call define_essential ( mesh_line1, input_probdef_line1, point=1, physq=1 )
    call define_essential ( mesh_line2, input_probdef_line2, point=1, physq=1 )
    call define_essential ( mesh_line3, input_probdef_line3, point=1, physq=1 )
    call define_essential ( mesh_line4, input_probdef_line4, point=1, physq=1 )

!   define problem
    call problem_definition ( input_probdef_line1, mesh_line1, problem_line1 )
    call create_subscript ( mesh_line1, problem_line1, subsh1, physqarr=[1] )

    call create ( problem_line1, sol_line1, rhsd_line1 )
    call create ( problem_line1, sol_line1_n, sol_line1_nm1 )

    call problem_definition ( input_probdef_line2, mesh_line2, problem_line2 )
    call create_subscript ( mesh_line2, problem_line2, subsh2, physqarr=[1] )

    call create ( problem_line2, sol_line2, rhsd_line2 )
    call create ( problem_line2, sol_line2_n, sol_line2_nm1 )

    call problem_definition ( input_probdef_line3, mesh_line3, problem_line3 )
    call create_subscript ( mesh_line3, problem_line3, subsh3, physqarr=[1] )

    call create ( problem_line3, sol_line3, rhsd_line3 )
    call create ( problem_line3, sol_line3_n, sol_line3_nm1 )

    call problem_definition ( input_probdef_line4, mesh_line4, problem_line4 )
    call create_subscript ( mesh_line4, problem_line4, subsh4, physqarr=[1] )

    call create ( problem_line4, sol_line4, rhsd_line4 )
    call create ( problem_line4, sol_line4_n, sol_line4_nm1 )

!   fill solution vector with essential boundary conditions

    do i = 1,mesh%curves(19)%nnodes
      sol_line1(1)%u(subsh1%s(i)) = initial_line1y(i)
      sol_line1(2)%u(subsh1%s(i)) = initial_line1z(i)
    end do
    sol_line1_n(1)%u = sol_line1(1)%u
    sol_line1_n(2)%u = sol_line1(2)%u

    do i = 1,mesh%curves(20)%nnodes
      sol_line2(1)%u(subsh2%s(i)) = initial_line2y(i)
      sol_line2(2)%u(subsh2%s(i)) = initial_line2z(i)
    end do
    sol_line2_n(1)%u = sol_line2(1)%u
    sol_line2_n(2)%u = sol_line2(2)%u

    do i = 1,mesh%curves(17)%nnodes
      sol_line3(1)%u(subsh3%s(i)) = initial_line3y(i)
      sol_line3(2)%u(subsh3%s(i)) = initial_line3z(i)
    end do
    sol_line3_n(1)%u = sol_line3(1)%u
    sol_line3_n(2)%u = sol_line3(2)%u

    do i = 1,mesh%curves(18)%nnodes
      sol_line4(1)%u(subsh4%s(i)) = initial_line4y(i)
      sol_line4(2)%u(subsh4%s(i)) = initial_line4z(i)
    end do
    sol_line4_n(1)%u = sol_line4(1)%u
    sol_line4_n(2)%u = sol_line4(2)%u

!   create system matrix

    call create_sysmatrix_structure ( sysmatrix_line1, mesh_line1, &
      problem_line1 )
    call create_sysmatrix_data ( sysmatrix_line1 )

    call create_sysmatrix_structure ( sysmatrix_line2, mesh_line2, &
      problem_line2 )
    call create_sysmatrix_data ( sysmatrix_line2 )

    call create_sysmatrix_structure ( sysmatrix_line3, mesh_line3, &
      problem_line3 )
    call create_sysmatrix_data ( sysmatrix_line3 )

    call create_sysmatrix_structure ( sysmatrix_line4, mesh_line4, &
      problem_line4 )
    call create_sysmatrix_data ( sysmatrix_line4 )

!   create vectors

    call create_vector ( problem_line1, velocity_line1, vec=2 )
    velocity_line1%u = 0 ! initialize

    call create_vector ( problem_line2, velocity_line2, vec=2 )
    velocity_line2%u = 0 ! initialize

    call create_vector ( problem_line3, velocity_line3, vec=2 )
    velocity_line3%u = 0 ! initialize

    call create_vector ( problem_line4, velocity_line4, vec=2 )
    velocity_line4%u = 0 ! initialize

    call create_subscript ( mesh_line1, problem_line1, velx_line1, vec=2, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_line1, problem_line1, vely_line1, vec=2, &
        degfd=2, fillnodes=.true. )
    call create_subscript ( mesh_line1, problem_line1, velz_line1, vec=2, &
        degfd=3, fillnodes=.true. )

    call create_subscript ( mesh_line2, problem_line2, velx_line2, vec=2, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_line2, problem_line2, vely_line2, vec=2, &
        degfd=2, fillnodes=.true. )
    call create_subscript ( mesh_line2, problem_line2, velz_line2, vec=2, &
        degfd=3, fillnodes=.true. )

    call create_subscript ( mesh_line3, problem_line3, velx_line3, vec=2, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_line3, problem_line3, vely_line3, vec=2, &
        degfd=2, fillnodes=.true. )
    call create_subscript ( mesh_line3, problem_line3, velz_line3, vec=2, &
        degfd=3, fillnodes=.true. )

    call create_subscript ( mesh_line4, problem_line4, velx_line4, vec=2, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_line4, problem_line4, vely_line4, vec=2, &
        degfd=2, fillnodes=.true. )
    call create_subscript ( mesh_line4, problem_line4, velz_line4, vec=2, &
        degfd=3, fillnodes=.true. )

    call create_oldvectors ( oldvectors_line1, nsysvec1=2, nvec=1 )
    oldvectors_line1%s1(1)%p => sol_line1_n    ! corrector at n
    oldvectors_line1%s1(2)%p => sol_line1_nm1  ! corrector at nm1

    oldvectors_line1%v(1)%p => velocity_line1 ! advection velocity at np1

    call create_oldvectors ( oldvectors_line2, nsysvec1=2, nvec=1 )
    oldvectors_line2%s1(1)%p => sol_line2_n    ! corrector at n
    oldvectors_line2%s1(2)%p => sol_line2_nm1  ! corrector at nm1

    oldvectors_line2%v(1)%p => velocity_line2 ! advection velocity at np1

    call create_oldvectors ( oldvectors_line3, nsysvec1=2, nvec=1 )
    oldvectors_line3%s1(1)%p => sol_line3_n    ! corrector at n
    oldvectors_line3%s1(2)%p => sol_line3_nm1  ! corrector at nm1

    oldvectors_line3%v(1)%p => velocity_line3 ! advection velocity at np1

    call create_oldvectors ( oldvectors_line4, nsysvec1=2, nvec=1 )
    oldvectors_line4%s1(1)%p => sol_line4_n    ! corrector at n
    oldvectors_line4%s1(2)%p => sol_line4_nm1  ! corrector at nm1

    oldvectors_line4%v(1)%p => velocity_line4 ! advection velocity at np1

  end subroutine define_line_problems

! solve material line equations (corrector)

  subroutine solve_material_line_corrector

    integer :: i

    type(solver_options_ma41_t) :: solver_options_h

    velocity_line1%u(velx_line1%s) = sol%u(velxline1%s)
    velocity_line1%u(vely_line1%s) = sol%u(velyline1%s)
    velocity_line1%u(velz_line1%s) = sol%u(velzline1%s)

    velocity_line2%u(velx_line2%s) = sol%u(velxline2%s)
    velocity_line2%u(vely_line2%s) = sol%u(velyline2%s)
    velocity_line2%u(velz_line2%s) = sol%u(velzline2%s)

    velocity_line3%u(velx_line3%s) = sol%u(velxline3%s)
    velocity_line3%u(vely_line3%s) = sol%u(velyline3%s)
    velocity_line3%u(velz_line3%s) = sol%u(velzline3%s)

    velocity_line4%u(velx_line4%s) = sol%u(velxline4%s)
    velocity_line4%u(vely_line4%s) = sol%u(velyline4%s)
    velocity_line4%u(velz_line4%s) = sol%u(velzline4%s)

    call build_system ( mesh_line1, problem_line1, sysmatrix_line1, &
      msysvector=rhsd_line1, elemsub=surface_material_line_elem, &
      oldvectors=oldvectors_line1, coefficients=coefficients_mat )

    call check_filled_sysmatrix ( sysmatrix_line1 )

!   MA41 solver storage

    solver_options_h%integer_storage = 2.0
    solver_options_h%real_storage    = 2.0

!   solve system
    do i = 1,2
      call add_effect_of_essential_to_rhs ( problem_line1, sysmatrix_line1, &
                              sol_line1(i), rhsd_line1(i) )
      call solve_system_ma41 ( sysmatrix_line1, rhsd_line1(i), sol_line1(i), &
                              luc,  solver_options=solver_options_h )
    end do

    call delete (luc)

    call build_system ( mesh_line2, problem_line2, sysmatrix_line2, &
      msysvector=rhsd_line2, elemsub=surface_material_line_elem, &
      oldvectors=oldvectors_line2, coefficients=coefficients_mat )

    call check_filled_sysmatrix ( sysmatrix_line2 )

    do i = 1,2
      call add_effect_of_essential_to_rhs ( problem_line2, sysmatrix_line2, &
                             sol_line2(i), rhsd_line2(i) )
      call solve_system_ma41 ( sysmatrix_line2, rhsd_line2(i), sol_line2(i), &
                             luc, solver_options=solver_options_h )
    end do

    call delete (luc)

    call build_system ( mesh_line3, problem_line3, sysmatrix_line3, &
      msysvector=rhsd_line3, elemsub=surface_material_line_elem, &
      oldvectors=oldvectors_line3, coefficients=coefficients_mat )

    call check_filled_sysmatrix ( sysmatrix_line3 )

    do i = 1,2
      call add_effect_of_essential_to_rhs ( problem_line3, sysmatrix_line3, &
                              sol_line3(i), rhsd_line3(i) )
      call solve_system_ma41 ( sysmatrix_line3, rhsd_line3(i), sol_line3(i), &
                              luc, solver_options=solver_options_h )
    end do

    call delete (luc)

    call build_system ( mesh_line4, problem_line4, sysmatrix_line4, &
      msysvector=rhsd_line4, elemsub=surface_material_line_elem, &
      oldvectors=oldvectors_line4, coefficients=coefficients_mat )

    call check_filled_sysmatrix ( sysmatrix_line4 )

    do i = 1,2
      call add_effect_of_essential_to_rhs ( problem_line4, sysmatrix_line4, &
                             sol_line4(i), rhsd_line4(i) )
      call solve_system_ma41 ( sysmatrix_line4, rhsd_line4(i), sol_line4(i), &
                             luc, solver_options=solver_options_h )
    end do

    call delete (luc)

  end subroutine solve_material_line_corrector

  subroutine define_surface_problems

    integer :: i

    call create_surface_mesh

!   define some arrays for the ALE mesh position at old times
    allocate ( meshcoorsurf1_n(mesh_surf1%nnodes,mesh_surf1%ndim), &
        meshcoorsurf1_nm1(mesh_surf1%nnodes,mesh_surf1%ndim) )
    allocate ( meshcoorsurf2_n(mesh_surf2%nnodes,mesh_surf2%ndim), &
        meshcoorsurf2_nm1(mesh_surf2%nnodes,mesh_surf2%ndim) )
    allocate ( meshcoorsurf3_n(mesh_surf3%nnodes,mesh_surf3%ndim), &
        meshcoorsurf3_nm1(mesh_surf3%nnodes,mesh_surf3%ndim) )
    allocate ( meshcoorsurf4_n(mesh_surf4%nnodes,mesh_surf4%ndim), &
        meshcoorsurf4_nm1(mesh_surf4%nnodes,mesh_surf4%ndim) )

    allocate(initial_hsurf1(mesh_surf1%nnodes))
    allocate(initial_hsurf2(mesh_surf2%nnodes))
    allocate(initial_hsurf3(mesh_surf3%nnodes))
    allocate(initial_hsurf4(mesh_surf4%nnodes))

    initial_hsurf1 = mesh%coor(mesh%surfaces(8)%nodes,3)
    initial_hsurf2 = mesh%coor(mesh%surfaces(9)%nodes,2)
    initial_hsurf3 = mesh%coor(mesh%surfaces(10)%nodes,3)
    initial_hsurf4 = mesh%coor(mesh%surfaces(11)%nodes,2)

    allocate(initial_hsurf1_line(mesh%curves(7)%nnodes))
    allocate(initial_hsurf2_line(mesh%curves(8)%nnodes))
    allocate(initial_hsurf3_line(mesh%curves(5)%nnodes))
    allocate(initial_hsurf4_line(mesh%curves(6)%nnodes))

    do i = 1,mesh%curves(7)%nnodes
      initial_hsurf1_line(i) = mesh%coor(mesh%curves(7)%nodes(i),3)
    end do

    do i = 1, mesh%curves(8)%nnodes
      initial_hsurf2_line(i) = mesh%coor(mesh%curves(8)%nodes(i),2)
    end do

    do i = 1, mesh%curves(5)%nnodes
      initial_hsurf3_line(i) = mesh%coor(mesh%curves(5)%nodes(i),3)
    end do

    do i = 1, mesh%curves(6)%nnodes
      initial_hsurf4_line(i) = mesh%coor(mesh%curves(6)%nodes(i),2)
    end do

    allocate(surf1hat(mesh_surf1%nnodes,2),surf1hatn(mesh_surf1%nnodes,2))
    allocate(surf2hat(mesh_surf2%nnodes,2),surf2hatn(mesh_surf2%nnodes,2))
    allocate(surf3hat(mesh_surf3%nnodes,2),surf3hatn(mesh_surf3%nnodes,2))
    allocate(surf4hat(mesh_surf4%nnodes,2),surf4hatn(mesh_surf4%nnodes,2))

    surf1hat(:,2) = initial_hsurf1
    surf1hat(:,1) = mesh%coor(mesh%surfaces(8)%nodes(:),2)
    surf2hat(:,2) = initial_hsurf2
    surf2hat(:,1) = mesh%coor(mesh%surfaces(9)%nodes(:),3)
    surf3hat(:,2) = initial_hsurf3
    surf3hat(:,1) = mesh%coor(mesh%surfaces(10)%nodes(:),2)
    surf4hat(:,2) = initial_hsurf4
    surf4hat(:,1) = mesh%coor(mesh%surfaces(11)%nodes(:),3)

    surf1hatn = surf1hat
    surf2hatn = surf2hat
    surf3hatn = surf3hat
    surf4hatn = surf4hat

!   problem definition for surface advection

    call create_input_probdef ( mesh_surf1, input_probdef_surf1, nvec=3, &
      nphysq=1 )

    input_probdef_surf1%vec_elementdof(1)%a = &
        reshape ( [ 1,1,1,1,1,1,    &  ! height function
                    3,3,3,3,3,3,    &  ! velocity
                    2,2,2,2,2,2],   &  ! mesh velocity
                  [6,3] )

    input_probdef_surf1%physq = [1]
    input_probdef_surf1%probnr = 10

    call create_input_probdef ( mesh_surf2, input_probdef_surf2, nvec=3, &
      nphysq=1 )

    input_probdef_surf2%vec_elementdof(1)%a = &
        reshape ( [ 1,1,1,1,1,1,    &  ! height function
                    3,3,3,3,3,3,    &  ! velocity
                    2,2,2,2,2,2],   &  ! mesh velocity
                  [6,3] )

    input_probdef_surf2%physq = [1]
    input_probdef_surf2%probnr = 11

    call create_input_probdef ( mesh_surf3, input_probdef_surf3, nvec=3, &
      nphysq=1 )

    input_probdef_surf3%vec_elementdof(1)%a = &
        reshape ( [ 1,1,1,1,1,1,    &  ! height function
                    3,3,3,3,3,3,    &  ! velocity
                    2,2,2,2,2,2],   &  ! mesh velocity
                  [6,3] )

    input_probdef_surf3%physq = [1]
    input_probdef_surf3%probnr = 12

    call create_input_probdef ( mesh_surf4, input_probdef_surf4, nvec=3, &
      nphysq=1 )

    input_probdef_surf4%vec_elementdof(1)%a = &
        reshape ( [ 1,1,1,1,1,1,    &  ! height function
                    3,3,3,3,3,3,    &  ! velocity
                    2,2,2,2,2,2],   &  ! mesh velocity
                  [6,3] )

    input_probdef_surf4%physq = [1]
    input_probdef_surf4%probnr = 13

    call define_essential ( mesh_surf1, input_probdef_surf1, curve1=1, &
         physq=1 )
    call define_essential ( mesh_surf2, input_probdef_surf2, curve1=1, &
         physq=1 )
    call define_essential ( mesh_surf3, input_probdef_surf3, curve1=1, &
         physq=1 )
    call define_essential ( mesh_surf4, input_probdef_surf4, curve1=1, &
         physq=1 )

!   define problem
    call problem_definition ( input_probdef_surf1, mesh_surf1, problem_surf1 )

    call create_subscript ( mesh_surf1, problem_surf1, subs_surf1, &
         physqarr=[1] )
    call create_subscript ( mesh_surf1, problem_surf1, subs_curve1, &
         physqarr=[1], curves=[1])

    call create ( problem_surf1, sol_surf1, rhsd_surf1 )
    call create ( problem_surf1, sol_surf1_n, sol_surf1_nm1 )

    call problem_definition ( input_probdef_surf2, mesh_surf2, problem_surf2 )

    call create_subscript ( mesh_surf2, problem_surf2, subs_surf2, &
         physqarr=[1] )
    call create_subscript ( mesh_surf2, problem_surf2, subs_curve2, &
         physqarr=[1], curves=[1])

    call create ( problem_surf2, sol_surf2, rhsd_surf2 )
    call create ( problem_surf2, sol_surf2_n, sol_surf2_nm1 )

    call problem_definition ( input_probdef_surf3, mesh_surf3, problem_surf3 )

    call create_subscript ( mesh_surf3, problem_surf3, subs_surf3, &
         physqarr=[1] )
    call create_subscript ( mesh_surf3, problem_surf3, subs_curve3, &
         physqarr=[1], curves=[1])

    call create ( problem_surf3, sol_surf3, rhsd_surf3 )
    call create ( problem_surf3, sol_surf3_n, sol_surf3_nm1 )

    call problem_definition ( input_probdef_surf4, mesh_surf4, problem_surf4 )

    call create_subscript ( mesh_surf4, problem_surf4, subs_surf4, &
         physqarr=[1] )
    call create_subscript ( mesh_surf4, problem_surf4, subs_curve4, &
         physqarr=[1], curves=[1])

    call create ( problem_surf4, sol_surf4, rhsd_surf4 )
    call create ( problem_surf4, sol_surf4_n, sol_surf4_nm1 )

!   fill solution vector with essential boundary conditions

    do i = 1, mesh_surf1%nnodes
      sol_surf1%u(subs_surf1%s(i)) = initial_hsurf1(i)
    end do
    sol_surf1_n%u = sol_surf1%u

    do i = 1, mesh_surf2%nnodes
      sol_surf2%u(subs_surf2%s(i)) = initial_hsurf2(i)
    end do
    sol_surf2_n%u = sol_surf2%u

    do i = 1, mesh_surf3%nnodes
      sol_surf3%u(subs_surf3%s(i)) = initial_hsurf3(i)
    end do
    sol_surf3_n%u = sol_surf3%u

    do i = 1, mesh_surf4%nnodes
      sol_surf4%u(subs_surf4%s(i)) = initial_hsurf4(i)
    end do
    sol_surf4_n%u = sol_surf4%u

!   fix z/y position of points right after die

    sol_surf1%u(subs_curve1%s) = initial_hsurf1_line
    sol_surf2%u(subs_curve2%s) = initial_hsurf2_line
    sol_surf3%u(subs_curve3%s) = initial_hsurf3_line
    sol_surf4%u(subs_curve4%s) = initial_hsurf4_line

!   create system matrix

    call create_sysmatrix_structure ( sysmatrix_surf1, mesh_surf1, &
      problem_surf1 )
    call create_sysmatrix_data ( sysmatrix_surf1 )

    call create_sysmatrix_structure ( sysmatrix_surf2, mesh_surf2, &
      problem_surf2 )
    call create_sysmatrix_data ( sysmatrix_surf2 )

    call create_sysmatrix_structure ( sysmatrix_surf3, mesh_surf3, &
      problem_surf3 )
    call create_sysmatrix_data ( sysmatrix_surf3 )

    call create_sysmatrix_structure ( sysmatrix_surf4, mesh_surf4, &
      problem_surf4 )
    call create_sysmatrix_data ( sysmatrix_surf4 )

!   create vectors

    call create_vector ( problem_surf1, velocity_surf1, vec=2 )
    velocity_surf1%u = 0 ! initialize
    call create_vector ( problem_surf1, meshvel_surf1, vec=3 )
    meshvel_surf1%u = 0

    call create_vector ( problem_surf2, velocity_surf2, vec=2 )
    velocity_surf2%u = 0 ! initialize
    call create_vector ( problem_surf2, meshvel_surf2, vec=3 )
    meshvel_surf2%u = 0

    call create_vector ( problem_surf3, velocity_surf3, vec=2 )
    velocity_surf3%u = 0 ! initialize
    call create_vector ( problem_surf3, meshvel_surf3, vec=3 )
    meshvel_surf3%u = 0

    call create_vector ( problem_surf4, velocity_surf4, vec=2 )
    velocity_surf4%u = 0 ! initialize
    call create_vector ( problem_surf4, meshvel_surf4, vec=3 )
    meshvel_surf4%u = 0

    call create_subscript ( mesh_surf1, problem_surf1, velx_surf1, vec=2, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_surf1, problem_surf1, vely_surf1, vec=2, &
        degfd=2, fillnodes=.true. )
    call create_subscript ( mesh_surf1, problem_surf1, velz_surf1, vec=2, &
        degfd=3, fillnodes=.true. )

    call create_subscript ( mesh_surf2, problem_surf2, velx_surf2, vec=2, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_surf2, problem_surf2, vely_surf2, vec=2, &
        degfd=2, fillnodes=.true. )
    call create_subscript ( mesh_surf2, problem_surf2, velz_surf2, vec=2, &
        degfd=3, fillnodes=.true. )

    call create_subscript ( mesh_surf3, problem_surf3, velx_surf3, vec=2, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_surf3, problem_surf3, vely_surf3, vec=2, &
        degfd=2, fillnodes=.true. )
    call create_subscript ( mesh_surf3, problem_surf3, velz_surf3, vec=2, &
        degfd=3, fillnodes=.true. )

    call create_subscript ( mesh_surf4, problem_surf4, velx_surf4, vec=2, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_surf4, problem_surf4, vely_surf4, vec=2, &
        degfd=2, fillnodes=.true. )
    call create_subscript ( mesh_surf4, problem_surf4, velz_surf4, vec=2, &
        degfd=3, fillnodes=.true. )

    call create_subscript ( mesh_surf1, problem_surf1, meshvelx_surf1, vec=3, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_surf1, problem_surf1, meshvely_surf1, vec=3, &
        degfd=2, fillnodes=.true. )

    call create_subscript ( mesh_surf2, problem_surf2, meshvelx_surf2, vec=3, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_surf2, problem_surf2, meshvely_surf2, vec=3, &
        degfd=2, fillnodes=.true. )

    call create_subscript ( mesh_surf3, problem_surf3, meshvelx_surf3, vec=3, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_surf3, problem_surf3, meshvely_surf3, vec=3, &
        degfd=2, fillnodes=.true. )

    call create_subscript ( mesh_surf4, problem_surf4, meshvelx_surf4, vec=3, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_surf4, problem_surf4, meshvely_surf4, vec=3, &
        degfd=2, fillnodes=.true. )

!   create subscript for the height values

    call create ( mesh_surf2, problem_surf2, hgt )
    call create ( mesh_surf2, problem_surf2, hgt_end, curves=[2] )

    call create_oldvectors ( oldvectors_surf1, nsysvec=2, nvec=1 )

    oldvectors_surf1%s(1)%p => sol_surf1_n    ! corrector at n
    oldvectors_surf1%s(2)%p => sol_surf1_nm1  ! corrector at nm1

    oldvectors_surf1%v(1)%p => velocity_surf1 ! advection velocity at np1

    call create_oldvectors ( oldvectors_surf2, nsysvec=2, nvec=1 )
    oldvectors_surf2%s(1)%p => sol_surf2_n    ! corrector at n
    oldvectors_surf2%s(2)%p => sol_surf2_nm1  ! corrector at nm1

    oldvectors_surf2%v(1)%p => velocity_surf2 ! advection velocity at np1

    call create_oldvectors ( oldvectors_surf3, nsysvec=2, nvec=1 )
    oldvectors_surf3%s(1)%p => sol_surf3_n    ! corrector at n
    oldvectors_surf3%s(2)%p => sol_surf3_nm1  ! corrector at nm1

    oldvectors_surf3%v(1)%p => velocity_surf3 ! advection velocity at np1

    call create_oldvectors ( oldvectors_surf4, nsysvec=2, nvec=1 )
    oldvectors_surf4%s(1)%p => sol_surf4_n    ! corrector at n
    oldvectors_surf4%s(2)%p => sol_surf4_nm1  ! corrector at nm1

    oldvectors_surf4%v(1)%p => velocity_surf4 ! advection velocity at np1

  end subroutine define_surface_problems

  subroutine create_surface_mesh

    integer :: i, curve, node1, node2, node

!   create mesh for surface advection

    call mesh_skeleton ( mesh_surf1, mesh%surfaces(8)%nnodes, &
      mesh%surfaces(8)%nelem, mesh%surfaces(8)%element%elshape, ndim=2 )

    mesh_surf1%coor = mesh%coor(mesh%surfaces(8)%nodes,[1,2])
    mesh_surf1%topology(1)%a = mesh%surfaces(8)%topology(:,:,1)

    allocate(bd(mesh_surf1%ndim,2))

    do i=1,mesh_surf1%ndim
      bd(i,1)=minval(mesh_surf1%coor(:,i))-1.e-10_dp
      bd(i,2)=maxval(mesh_surf1%coor(:,i))+1.e-10_dp
    end do

    call add_to_mesh ( mesh_surf1, blocks=[(3,i=1,mesh_surf1%ndim)], &
         blocksdomain=bd )

    deallocate(bd)

!   add curves to the 2D mesh

    deallocate(mesh_surf1%curves)
    allocate(mesh_surf1%curves(4))
    mesh_surf1%ncurves = 4

    do curve = 1, 4

      if (curve == 1) curve2=7
      if (curve == 2) curve2=15
      if (curve == 3) curve2=20
      if (curve == 4) curve2=19

      mesh_surf1%curves(curve)%nnodes = mesh%curves(curve2)%nnodes
      mesh_surf1%curves(curve)%nelem = mesh%curves(curve2)%nelem
      mesh_surf1%curves(curve)%elnumnod = mesh%curves(curve2)%elnumnod
      allocate ( mesh_surf1%curves(curve)%nodes(mesh_surf1%curves(&
        curve)%nnodes))
      allocate ( mesh_surf1%curves(curve)%topology( &
        mesh%curves(curve2)%element%numnod, &
        mesh_surf1%curves(curve)%nelem, 2))
        mesh_surf1%curves(curve)%nodes = mesh%curves(curve2)%nodes
        mesh_surf1%curves(curve)%topology = mesh%curves(curve2)%topology
        mesh_surf1%curves(curve)%element%elshape = &
                              mesh%curves(curve2)%element%elshape
        mesh_surf1%curves(curve)%ndim = 2

      do node2 = 1, mesh%curves(curve2)%nnodes

        nodenr2 = mesh%curves(curve2)%nodes(node2)

        do node1 = 1, mesh_surf1%nnodes

          if ( abs(mesh%coor(nodenr2,1)-mesh_surf1%coor(node1,1))<1e-12_dp &
               .and. abs(mesh%coor(nodenr2,2)-mesh_surf1%coor(node1,2))< &
               1e-12_dp ) then

!           matching node
            mesh_surf1%curves(curve)%nodes(node2) = node1

            exit ! leave loop

          end if

        end do

      end do

!     update the global node numbers in the toplogy
      do node = 1, mesh_surf1%curves(curve)%nnodes
        do i = 1, mesh_surf1%curves(curve)%nelem
          mesh_surf1%curves(curve)%topology(:,i,2) = &
            mesh_surf1%curves(curve)%nodes(mesh_surf1%curves(&
            curve)%topology(:,i,1))
        end do
      end do

      if ( any( mesh_surf1%curves(curve)%nodes == 0 ) ) then
        print *,'Error: no corresponding node found!'
        stop
      end if

      allocate ( mesh_surf1%curves(curve)%element_blend(0) )
      allocate ( mesh_surf1%curves(curve)%nnodes_blend(2) )
      mesh_surf1%curves(curve)%nnodes_blend = [0, &
            mesh_surf1%curves(curve)%nnodes ]

   end do

   WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = .false.

   call find_bounds_blocks ( mesh_surf1 )

   call fill_mesh_parts ( mesh_surf1 )

   call write_mesh_vtk ( mesh_surf1, 'mesh_surf1.vtk' )

   call mesh_skeleton ( mesh_surf2, mesh%surfaces(9)%nnodes, &
     mesh%surfaces(9)%nelem, mesh%surfaces(9)%element%elshape, ndim=2 )

    mesh_surf2%coor = mesh%coor(mesh%surfaces(9)%nodes,[1,3])
    mesh_surf2%topology(1)%a = mesh%surfaces(9)%topology(:,:,1)

    allocate(bd(mesh_surf2%ndim,2))

    do i=1,mesh_surf2%ndim
      bd(i,1)=minval(mesh_surf2%coor(:,i))-1.e-10_dp
      bd(i,2)=maxval(mesh_surf2%coor(:,i))+1.e-10_dp
    end do

    call add_to_mesh ( mesh_surf2, blocks=[(3,i=1,mesh_surf2%ndim)], &
         blocksdomain=bd )

    deallocate(bd)

!   add curves to the 2D mesh

    deallocate(mesh_surf2%curves)
    allocate(mesh_surf2%curves(4))
    mesh_surf2%ncurves = 4

    do curve = 1, 4

      if (curve == 1) curve2=8
      if (curve == 2) curve2=16
      if (curve == 3) curve2=20
      if (curve == 4) curve2=17

      mesh_surf2%curves(curve)%nnodes = mesh%curves(curve2)%nnodes
      mesh_surf2%curves(curve)%nelem = mesh%curves(curve2)%nelem
      mesh_surf2%curves(curve)%elnumnod = mesh%curves(curve2)%elnumnod
      allocate ( mesh_surf2%curves(curve)%nodes(mesh_surf2%curves(&
          curve)%nnodes))
      allocate ( mesh_surf2%curves(curve)%topology( &
          mesh%curves(curve2)%element%numnod, &
          mesh_surf2%curves(curve)%nelem, 2))
      mesh_surf2%curves(curve)%nodes = mesh%curves(curve2)%nodes
      mesh_surf2%curves(curve)%topology = mesh%curves(curve2)%topology
      mesh_surf2%curves(curve)%element%elshape = &
                               mesh%curves(curve2)%element%elshape
      mesh_surf2%curves(curve)%ndim = 2

      do node2 = 1, mesh%curves(curve2)%nnodes

        nodenr2 = mesh%curves(curve2)%nodes(node2)

        do node1 = 1, mesh_surf2%nnodes

          if ( abs(mesh%coor(nodenr2,1)-mesh_surf2%coor(node1,1))<1e-12_dp &
               .and. abs(mesh%coor(nodenr2,3)-mesh_surf2%coor(node1,2))< &
               1e-12_dp ) then

!           matching node
            mesh_surf2%curves(curve)%nodes(node2) = node1

            exit ! leave loop

          end if

        end do

      end do

!     update the global node numbers in the toplogy

      do node = 1, mesh_surf2%curves(curve)%nnodes
        do i = 1, mesh_surf2%curves(curve)%nelem
          mesh_surf2%curves(curve)%topology(:,i,2) = &
            mesh_surf2%curves(curve)%nodes(mesh_surf2%curves(&
             curve)%topology(:,i,1))
        end do
      end do

      if ( any( mesh_surf2%curves(curve)%nodes == 0 ) ) then
        print *,'Error: no corresponding node found!'
        stop
      end if

      allocate ( mesh_surf2%curves(curve)%element_blend(0) )
      allocate ( mesh_surf2%curves(curve)%nnodes_blend(2) )
      mesh_surf2%curves(curve)%nnodes_blend = [0, &
            mesh_surf2%curves(curve)%nnodes ]

    end do

    WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = .false.

    call find_bounds_blocks ( mesh_surf2 )

    call fill_mesh_parts ( mesh_surf2 )

    call write_mesh_vtk ( mesh_surf2, 'mesh_surf2.vtk' )

    call mesh_skeleton ( mesh_surf3, mesh%surfaces(10)%nnodes, &
      mesh%surfaces(10)%nelem, mesh%surfaces(10)%element%elshape, ndim=2 )

    mesh_surf3%coor = mesh%coor(mesh%surfaces(10)%nodes,[1,2])
    mesh_surf3%topology(1)%a = mesh%surfaces(10)%topology(:,:,1)

    allocate(bd(mesh_surf3%ndim,2))

    do i=1,mesh_surf3%ndim
      bd(i,1)=minval(mesh_surf3%coor(:,i))-1.e-10_dp
      bd(i,2)=maxval(mesh_surf3%coor(:,i))+1.e-10_dp
    end do

    call add_to_mesh ( mesh_surf3, blocks=[(3,i=1,mesh_surf3%ndim)], &
         blocksdomain=bd )

    deallocate(bd)

!   add curves to the 2D mesh

    deallocate(mesh_surf3%curves)
    allocate(mesh_surf3%curves(4))
    mesh_surf3%ncurves = 4

    do curve = 1, 4

      if (curve == 1) curve2=5
      if (curve == 2) curve2=13
      if (curve == 3) curve2=17
      if (curve == 4) curve2=18

      mesh_surf3%curves(curve)%nnodes = mesh%curves(curve2)%nnodes
      mesh_surf3%curves(curve)%nelem = mesh%curves(curve2)%nelem
      mesh_surf3%curves(curve)%elnumnod = mesh%curves(curve2)%elnumnod
      allocate ( mesh_surf3%curves(curve)%nodes(mesh_surf3%curves(&
          curve)%nnodes))
      allocate ( mesh_surf3%curves(curve)%topology( &
         mesh%curves(curve2)%element%numnod, &
         mesh_surf3%curves(curve)%nelem, 2))
      mesh_surf3%curves(curve)%nodes = mesh%curves(curve2)%nodes
      mesh_surf3%curves(curve)%topology = mesh%curves(curve2)%topology
      mesh_surf3%curves(curve)%element%elshape = &
                               mesh%curves(curve2)%element%elshape
      mesh_surf3%curves(curve)%ndim = 2

      do node2 = 1, mesh%curves(curve2)%nnodes

        nodenr2 = mesh%curves(curve2)%nodes(node2)

        do node1 = 1, mesh_surf3%nnodes

          if ( abs(mesh%coor(nodenr2,1)-mesh_surf3%coor(node1,1))<1e-12_dp &
               .and. abs(mesh%coor(nodenr2,2)-mesh_surf3%coor(node1,2))< &
               1e-12_dp ) then

!           matching node
            mesh_surf3%curves(curve)%nodes(node2) = node1

            exit ! leave loop

          end if

        end do

      end do

!     update the global node numbers in the toplogy
      do node = 1, mesh_surf3%curves(curve)%nnodes
        do i = 1, mesh_surf3%curves(curve)%nelem
          mesh_surf3%curves(curve)%topology(:,i,2) = &
            mesh_surf3%curves(curve)%nodes(mesh_surf3%curves(&
            curve)%topology(:,i,1))
        end do
      end do

      if ( any( mesh_surf3%curves(curve)%nodes == 0 ) ) then
        print *,'Error: no corresponding node found!'
        stop
      end if

      allocate ( mesh_surf3%curves(curve)%element_blend(0) )
      allocate ( mesh_surf3%curves(curve)%nnodes_blend(2) )
      mesh_surf3%curves(curve)%nnodes_blend = [0, &
            mesh_surf3%curves(curve)%nnodes ]

    end do

    WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = .false.

    call find_bounds_blocks ( mesh_surf3 )

    call fill_mesh_parts ( mesh_surf3 )

    call write_mesh_vtk ( mesh_surf3, 'mesh_surf3.vtk' )

    call mesh_skeleton ( mesh_surf4, mesh%surfaces(11)%nnodes, &
      mesh%surfaces(11)%nelem, mesh%surfaces(11)%element%elshape, ndim=2 )

    mesh_surf4%coor = mesh%coor(mesh%surfaces(11)%nodes,[1,3])
    mesh_surf4%topology(1)%a = mesh%surfaces(11)%topology(:,:,1)

    allocate(bd(mesh_surf4%ndim,2))

    do i=1,mesh_surf4%ndim
      bd(i,1)=minval(mesh_surf4%coor(:,i))-1.e-10_dp
      bd(i,2)=maxval(mesh_surf4%coor(:,i))+1.e-10_dp
    end do

    call add_to_mesh ( mesh_surf4, blocks=[(3,i=1,mesh_surf4%ndim)], &
         blocksdomain=bd )

    deallocate(bd)

!   add curves to the 2D mesh

    deallocate(mesh_surf4%curves)
    allocate(mesh_surf4%curves(4))
    mesh_surf4%ncurves = 4

    do curve = 1, 4

       if (curve == 1) curve2=6
       if (curve == 2) curve2=14
       if (curve == 3) curve2=19
       if (curve == 4) curve2=18

      mesh_surf4%curves(curve)%nnodes = mesh%curves(curve2)%nnodes
      mesh_surf4%curves(curve)%nelem = mesh%curves(curve2)%nelem
      mesh_surf4%curves(curve)%elnumnod = mesh%curves(curve2)%elnumnod
      allocate ( mesh_surf4%curves(curve)%nodes(mesh_surf4%curves(&
          curve)%nnodes))
      allocate ( mesh_surf4%curves(curve)%topology( &
          mesh%curves(curve2)%element%numnod, &
          mesh_surf4%curves(curve)%nelem, 2))
      mesh_surf4%curves(curve)%nodes = mesh%curves(curve2)%nodes
      mesh_surf4%curves(curve)%topology = mesh%curves(curve2)%topology
      mesh_surf4%curves(curve)%element%elshape = &
                               mesh%curves(curve2)%element%elshape
      mesh_surf4%curves(curve)%ndim = 2

      do node2 = 1, mesh%curves(curve2)%nnodes

        nodenr2 = mesh%curves(curve2)%nodes(node2)

        do node1 = 1, mesh_surf4%nnodes

          if ( abs(mesh%coor(nodenr2,1)-mesh_surf4%coor(node1,1))<1e-12_dp &
               .and. abs(mesh%coor(nodenr2,3)-mesh_surf4%coor(node1,2))< &
               1e-12_dp ) then

!           matching node
            mesh_surf4%curves(curve)%nodes(node2) = node1

            exit ! leave loop

          end if

        end do

      end do

!     update the global node numbers in the toplogy
      do node = 1, mesh_surf4%curves(curve)%nnodes
        do i = 1, mesh_surf4%curves(curve)%nelem
          mesh_surf4%curves(curve)%topology(:,i,2) = &
            mesh_surf4%curves(curve)%nodes(mesh_surf4%curves(&
            curve)%topology(:,i,1))
        end do
      end do

      if ( any( mesh_surf4%curves(curve)%nodes == 0 ) ) then
        print *,'Error: no corresponding node found!'
        stop
      end if

      allocate ( mesh_surf4%curves(curve)%element_blend(0) )
      allocate ( mesh_surf4%curves(curve)%nnodes_blend(2) )
      mesh_surf4%curves(curve)%nnodes_blend = [0, &
            mesh_surf4%curves(curve)%nnodes ]

    end do

    WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = .false.

    call find_bounds_blocks ( mesh_surf4 )

    call fill_mesh_parts ( mesh_surf4 )

    call write_mesh_vtk ( mesh_surf4, 'mesh_surf4.vtk' )

  end subroutine create_surface_mesh

  subroutine expand_surfaces

    use update_mesh_nodes_surfaces_m

    if (step == 1) then
      meshcoorsurf1_n = mesh_surf1%coor
    else
      meshcoorsurf1_nm1 = meshcoorsurf1_n
      meshcoorsurf1_n = mesh_surf1%coor
     end if

    allocate(displine1(mesh_surf1%curves(3)%nnodes))
    allocate(displine2(mesh_surf1%curves(4)%nnodes))

    displine1 = sol_line2(1)%u(subsh2%s) - sol_line2_n(1)%u(subsh2%s)
    displine2 = sol_line1(1)%u(subsh1%s) - sol_line1_n(1)%u(subsh1%s)

    call update_mesh_nodes_expand ( mesh_surf1, problem_lapl, &
                     problem_laply, disp1=displine1, disp2=displine2 )

    deallocate(displine1)
    deallocate(displine2)
    allocate(displine2(mesh_surf2%curves(3)%nnodes))
    allocate(displine3(mesh_surf2%curves(4)%nnodes))

    if (step == 1) then
      meshcoorsurf2_n = mesh_surf2%coor
    else
      meshcoorsurf2_nm1 = meshcoorsurf2_n
      meshcoorsurf2_n = mesh_surf2%coor
    end if

    displine2 = sol_line2(2)%u(subsh2%s) - sol_line2_n(2)%u(subsh2%s)
    displine3 = sol_line3(2)%u(subsh3%s) - sol_line3_n(2)%u(subsh3%s)

    call update_mesh_nodes_expand ( mesh_surf2, problem_lapl, &
                     problem_laply, disp1=displine2, disp2=displine3 )

    call find_bounds_blocks ( mesh_surf2 )

    deallocate(displine2)
    deallocate(displine3)
    allocate(displine3(mesh_surf3%curves(3)%nnodes))
    allocate(displine4(mesh_surf3%curves(4)%nnodes))

    if (step == 1) then
      meshcoorsurf3_n = mesh_surf3%coor
    else
      meshcoorsurf3_nm1 = meshcoorsurf3_n
      meshcoorsurf3_n = mesh_surf3%coor
    end if

    displine3 = sol_line3(1)%u(subsh3%s) - sol_line3_n(1)%u(subsh3%s)
    displine4 = sol_line4(1)%u(subsh4%s) - sol_line4_n(1)%u(subsh4%s)

    call update_mesh_nodes_expand ( mesh_surf3, problem_lapl, &
                     problem_laply, disp1=displine3, disp2=displine4 )

    call find_bounds_blocks ( mesh_surf3 )

    deallocate(displine3)
    deallocate(displine4)
    allocate(displine1(mesh_surf4%curves(3)%nnodes))
    allocate(displine4(mesh_surf4%curves(4)%nnodes))

    if (step == 1) then
      meshcoorsurf4_n = mesh_surf4%coor
    else
      meshcoorsurf4_nm1 = meshcoorsurf4_n
      meshcoorsurf4_n = mesh_surf4%coor
    end if

    displine1 = sol_line1(2)%u(subsh1%s) - sol_line1_n(2)%u(subsh1%s)
    displine4 = sol_line4(2)%u(subsh4%s) - sol_line4_n(2)%u(subsh4%s)

    call update_mesh_nodes_expand ( mesh_surf4, problem_lapl, &
                     problem_laply, disp1=displine1, disp2=displine4 )

    call find_bounds_blocks ( mesh_surf4 )

    deallocate(displine1)
    deallocate(displine4)

!   mesh velocity

    if ( step == 1  ) then
      meshvel_surf1%u = reshape ( transpose ( &
                  ( mesh_surf1%coor - meshcoorsurf1_n ) / deltat ), &
                                     [mesh_surf1%ndim*mesh_surf1%nnodes] )
      meshvel_surf2%u = reshape ( transpose ( &
                  ( mesh_surf2%coor - meshcoorsurf2_n ) / deltat ), &
                                     [mesh_surf2%ndim*mesh_surf2%nnodes] )
      meshvel_surf3%u = reshape ( transpose ( &
                  ( mesh_surf3%coor - meshcoorsurf3_n ) / deltat ), &
                                     [mesh_surf3%ndim*mesh_surf3%nnodes] )
      meshvel_surf4%u = reshape ( transpose ( &
                  ( mesh_surf4%coor - meshcoorsurf4_n ) / deltat ), &
                                     [mesh_surf4%ndim*mesh_surf4%nnodes] )
    else
      meshvel_surf1%u = reshape ( transpose ( &
          ( 1.5_dp*mesh_surf1%coor - 2*meshcoorsurf1_n + &
            0.5_dp*meshcoorsurf1_nm1 ) / deltat ), &
                           [mesh_surf1%ndim*mesh_surf1%nnodes] )
      meshvel_surf2%u = reshape ( transpose ( &
          ( 1.5_dp*mesh_surf2%coor - 2*meshcoorsurf2_n + &
            0.5_dp*meshcoorsurf2_nm1 ) / deltat ), &
                           [mesh_surf2%ndim*mesh_surf2%nnodes] )
      meshvel_surf3%u = reshape ( transpose ( &
        ( 1.5_dp*mesh_surf3%coor - 2*meshcoorsurf3_n + &
          0.5_dp*meshcoorsurf3_nm1 ) / deltat ), &
                           [mesh_surf3%ndim*mesh_surf3%nnodes] )
      meshvel_surf4%u = reshape ( transpose ( &
          ( 1.5_dp*mesh_surf4%coor - 2*meshcoorsurf4_n + &
            0.5_dp*meshcoorsurf4_nm1 ) / deltat ), &
                           [mesh_surf4%ndim*mesh_surf4%nnodes] )
    end if

    call copy ( sol_line1_n, sol_line1_nm1 )
    call copy ( sol_line1, sol_line1_n )
    call copy ( sol_line2_n, sol_line2_nm1 )
    call copy ( sol_line2, sol_line2_n )
    call copy ( sol_line3_n, sol_line3_nm1 )
    call copy ( sol_line3, sol_line3_n )
    call copy ( sol_line4_n, sol_line4_nm1 )
    call copy ( sol_line4, sol_line4_n )

  end subroutine expand_surfaces

! solve convection equation for the surface height (corrector)

  subroutine solve_4surfaces_height_corrector

    use postprocessing_m

    type(solver_options_ma41_t) :: solver_options_h

    velocity_surf1%u(velx_surf1%s) = sol%u(velxsurf1%s) - &
                         meshvel_surf1%u(meshvelx_surf1%s)
    velocity_surf1%u(vely_surf1%s) = sol%u(velysurf1%s) - &
                         meshvel_surf1%u(meshvely_surf1%s)
    velocity_surf1%u(velz_surf1%s) = sol%u(velzsurf1%s)

    velocity_surf2%u(velx_surf2%s) = sol%u(velxsurf2%s) - &
                         meshvel_surf2%u(meshvelx_surf2%s)
    velocity_surf2%u(vely_surf2%s) = sol%u(velzsurf2%s) - &
                         meshvel_surf2%u(meshvely_surf2%s)
    velocity_surf2%u(velz_surf2%s) = sol%u(velysurf2%s)

    velocity_surf3%u(velx_surf3%s) = sol%u(velxsurf3%s) - &
                         meshvel_surf3%u(meshvelx_surf3%s)
    velocity_surf3%u(vely_surf3%s) = sol%u(velysurf3%s) - &
                         meshvel_surf3%u(meshvely_surf3%s)
    velocity_surf3%u(velz_surf3%s) = sol%u(velzsurf3%s)

    velocity_surf4%u(velx_surf4%s) = sol%u(velxsurf4%s) - &
                         meshvel_surf4%u(meshvelx_surf4%s)
    velocity_surf4%u(vely_surf4%s) = sol%u(velzsurf4%s) - &
                         meshvel_surf4%u(meshvely_surf4%s)
    velocity_surf4%u(velz_surf4%s) = sol%u(velysurf4%s)

    call build_system ( mesh_surf1, problem_surf1, sysmatrix_surf1, &
      rhsd_surf1, elemsub=surface_advection_elem, &
      oldvectors=oldvectors_surf1, coefficients=coefficients_sf_adv )

    call check_filled_sysmatrix ( sysmatrix_surf1 )

    call build_system ( mesh_surf2, problem_surf2, sysmatrix_surf2, &
      rhsd_surf2, elemsub=surface_advection_elem, &
      oldvectors=oldvectors_surf2, coefficients=coefficients_sf_adv )

    call check_filled_sysmatrix ( sysmatrix_surf2 )

    call build_system ( mesh_surf3, problem_surf3, sysmatrix_surf3, &
      rhsd_surf3, elemsub=surface_advection_elem, &
      oldvectors=oldvectors_surf3, coefficients=coefficients_sf_adv )

    call check_filled_sysmatrix ( sysmatrix_surf3 )

    call build_system ( mesh_surf4, problem_surf4, sysmatrix_surf4, &
      rhsd_surf4, elemsub=surface_advection_elem, &
      oldvectors=oldvectors_surf4, coefficients=coefficients_sf_adv )

    call check_filled_sysmatrix ( sysmatrix_surf4 )

    call add_effect_of_essential_to_rhs ( problem_surf1, sysmatrix_surf1, &
      sol_surf1, rhsd_surf1 )

    call add_effect_of_essential_to_rhs ( problem_surf2, sysmatrix_surf2, &
      sol_surf2, rhsd_surf2 )

    call add_effect_of_essential_to_rhs ( problem_surf3, sysmatrix_surf3, &
      sol_surf3, rhsd_surf3 )

    call add_effect_of_essential_to_rhs ( problem_surf4, sysmatrix_surf4, &
      sol_surf4, rhsd_surf4 )

!   MA41 solver storage

    solver_options_h%integer_storage = 2.0
    solver_options_h%real_storage    = 2.0

!   solve system

    call solve_system_ma41 ( sysmatrix_surf1, rhsd_surf1, sol_surf1, &
                             solver_options=solver_options_h )

!   update nodes in surfaces
    mesh%coor(mesh%surfaces(8)%nodes(:),2) = mesh_surf1%coor(:,2)
    mesh%coor(mesh%surfaces(8)%nodes(:),3) = sol_surf1%u(subs_surf1%s)

    call copy ( sol_surf1_n, sol_surf1_nm1 )
    call copy ( sol_surf1, sol_surf1_n )

    call solve_system_ma41 ( sysmatrix_surf2, rhsd_surf2, sol_surf2, &
                             solver_options=solver_options_h )

!   update nodes in surfaces
    mesh%coor(mesh%surfaces(9)%nodes(:),3) = mesh_surf2%coor(:,2)
    mesh%coor(mesh%surfaces(9)%nodes(:),2) = sol_surf2%u(subs_surf2%s)

    call copy ( sol_surf2_n, sol_surf2_nm1 )
    call copy ( sol_surf2, sol_surf2_n )

    call solve_system_ma41 ( sysmatrix_surf3, rhsd_surf3, sol_surf3, &
                             solver_options=solver_options_h )

!   update nodes in surfaces
    mesh%coor(mesh%surfaces(10)%nodes(:),2) = mesh_surf3%coor(:,2)
    mesh%coor(mesh%surfaces(10)%nodes(:),3) = sol_surf3%u(subs_surf3%s)

    call copy ( sol_surf3_n, sol_surf3_nm1 )
    call copy ( sol_surf3, sol_surf3_n )

    call solve_system_ma41 ( sysmatrix_surf4, rhsd_surf4, sol_surf4, &
                             solver_options=solver_options_h )

!   update nodes in surfaces
    mesh%coor(mesh%surfaces(11)%nodes(:),3) = mesh_surf4%coor(:,2)
    mesh%coor(mesh%surfaces(11)%nodes(:),2) = sol_surf4%u(subs_surf4%s)

    call copy ( sol_surf4_n, sol_surf4_nm1 )
    call copy ( sol_surf4, sol_surf4_n )

  end subroutine solve_4surfaces_height_corrector

end program extrudate_swell7
