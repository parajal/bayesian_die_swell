! Fountain flow problem (2D). Newtonian fluid.
! Free surface position determined using a height function approach.
! Walls move with velocity -U in x-direction (no slip).
! At the left boundary a parabolic profile is assumed with Q=0.
! Optional surface tension on moving surface.
! Similar to fountain_flow1, but now:
!   Q1isoQ2 for velocities
!   P1isoP2 line elements for height function.
!   Correction for velocity (Dirichlet) to make Q=0 for the approx. profile

program fountain_flow1b

  use tfem_m
  use stokes_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m
  use surface_advection_elements_m
  use update_mesh_nodes1_m
  use stream_function1_m
  use figplot_m

  implicit none

! constants flow problem

  integer, parameter :: &
!    uintpl = 8,  &  ! Q2 velocities
    uintpl = 10,  &  ! Q1isoQ2 velocities
    pintpl = 4,  &  ! Q1 pressures
    physqv = 1,  &  ! physical quantity nr of the velocities
    physqp = 2,  &  ! physical quantity nr of the pressures
!    nsubint = 0, &  ! single subintegration domain
    nsubint = 2,  & ! number of subintegration domains
    gauss = 3,   &  ! 3x3 integration of quads
!    nsubintb = 0, & ! single subintegration domain
    nsubintb = 2, & ! number of subintegration domains
    gaussb = 3      ! 3-point integration of boundary elements

! constants surface advection

  integer, parameter :: &
!    hintpl = 6,       & ! P2 height function
    hintpl = 9,       & ! P1isoP2 height function
!    nsubint_adv = 0,  & ! single subintegration domain
    nsubint_adv = 2,  & ! number of subintegration domains
!    vintpl = 0,       & ! velocity the same as height function
    vintpl = 9,       & ! P1isoP2 velocity
!    vintpl = 6,       & ! P2 velocity
    ninti_sf_adv = 3, & ! number of Gauss points
    ndisp = 1,        & ! number of coordinate directions for mesh update
    method = 1          ! discretization method 0: Galerkin, 1: SUPG


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, gammadot, tauviscous, Dtensor, Ltensor
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma57_t) :: solver_options_u
  type(subscript_t) :: vel14, velx, vely, vel

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

  type(problem_t), target :: problem_update_mesh
  real(dp), allocatable, dimension(:,:) :: disp

! stream function problem
  type(problem_t) :: problem_strm
  type(sysvector_t) :: sol_strm

! variables

  integer :: &
    numtimesteps = 125,  & ! number of time steps
    coorsys = 0, & ! planar Cartesian (0) or axisymmetric (1) coor. system
    plot_mesh_every = 10    ! plot mesh every .. steps

  real(dp) :: &
    eta = 1.0_dp,      & ! viscosity
    gamma_fa = 0.1_dp, & ! surface tension coefficient fluid-air interface
    U = 1.0_dp,        & ! -U is the velocity of the walls
    deltat = 5.e-2_dp, & ! time step
    betah = 0.5_dp,    & ! beta parameter for SUPG of the height function
!    betah = 1.0_dp,    & ! beta parameter for SUPG of the height function
    rs_up = 1.5_dp,    & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.0_dp       ! integer_storage velocity-pressure LU (HSL)

  integer :: vertices(4) = [1,3,5,7]
  integer :: step, nnodes
  real(dp) :: H, resultsum(1)
  real(dp) :: Ucorrect ! U correction for nonzero flowrate on entry

  real(dp), allocatable, dimension(:,:) :: meshcoor_n, meshcoor_nm1

  character(len=20) :: filename

  logical :: surface_tension = .false.  ! include surface tension


! namelist for input of variables; read from standard input

  namelist /comppar/ numtimesteps, plot_mesh_every, &
    eta, gamma_fa, U, deltat, rs_up, is_up, surface_tension

  read ( unit=*, nml=comppar )

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl, pintpl, 0, 0,     0, &
      physqv, physqp, 0, 0, gauss, &
      gaussb ]

  coefficients%i(23) = coorsys
  coefficients%i(32) = nsubint
  coefficients%i(33) = nsubintb

  coefficients%r = 0
  coefficients%r(1) = eta

! create mesh

  call read_mesh ( mesh, filename='mesh4b.bin' )

  if ( uintpl == 10 ) then
!   change elshape to 34 for isoparametric mapping
!   Not really needed here. Only needed for finding refcoors in object sampling
    mesh%element(:)%elshape = 34
  end if

  call fill_mesh_parts ( mesh )

  H = mesh%coor(mesh%points(9),2) ! half height of the domain


! problem definition of velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 2  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(vertices,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,4) = 3  ! symmetric tensor
  input_probdef%vec_elementdof(1)%a(:,5) = 4  ! unsymmetric tensor

  input_probdef%physq = [1,2]
  input_probdef%probnr = 1

! Dirichlet boundary conditions

! left boundary
  call define_essential ( mesh, input_probdef, &
    curve1=15, physq=physqv )
! walls
  call define_essential ( mesh, input_probdef, &
    curve1=16, physq=physqv )

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for the velocity components for monitoring

  call create_subscript ( mesh, problem, velx, physqarr=[physqv], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqv], degfd=2 )

! create vector subscript for the velocity for use with the stream function

  call create_subscript ( mesh, problem, vel, physqarr=[physqv] )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )
  oldvectors%s(1)%p => sol

! fill solution vector with essential boundary conditions

  sol%u = 0

  Ucorrect = 0

! entry profile without correction
  call fill_sysvector ( mesh, problem, sol, curve1=15, physq=physqv, degfd=1, &
    func=vel_entry, funcnr=1 )

! flow rate of imposed profile
  call integrate_boundary_elements ( mesh, problem, resultsum, &
    elemsub=stokes_flowrate, curve=15, coefficients=coefficients, &
    oldvectors=oldvectors )

  Ucorrect = resultsum(1) / ( 2 * H )

! entry profile with correction to set Q=0 exactly.
  call fill_sysvector ( mesh, problem, sol, curve1=15, physq=physqv, degfd=1, &
    func=vel_entry, funcnr=1 )

! walls
  call fill_sysvector ( mesh, problem, sol, curve1=16, physq=physqv, degfd=1, &
    value=-U+Ucorrect )

! create system matrix

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


! create subscript for velocity sampling on curve 14

  call create_subscript ( mesh, problem, vel14, physqarr=[physqv], curves=[14] )

! fill coefficients for surface advection problem

  call create_coefficients ( coefficients_sf_adv, ncoefi=100, ncoefr=50 )

  coefficients_sf_adv%i = 0
  coefficients_sf_adv%i(1:6) = [ ninti_sf_adv, 2, 0, method, 0, hintpl ]

  coefficients_sf_adv%i(8) = nsubint_adv
  coefficients_sf_adv%i(10) = vintpl

  coefficients_sf_adv%r = 0

  coefficients_sf_adv%r(4) = deltat
  coefficients_sf_adv%r(5) = betah

! create mesh for surface advection

  call curve_to_1Dmesh ( mesh, mesh%curves(14), dim1D=2, pnts=[2,9], &
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

  allocate ( disp(mesh_sf_adv%nnodes,ndisp) )


! open files

  open ( unit=11, file='height.out', status='unknown' )
  open ( unit=23, file='out', recl=300 )

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
        curve=14, coefficients=coefficients, physq=[physqv] )

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

!   write monitor data

    write(unit=23,fmt=*) &
      step, step*deltat, maxval(abs(sol%u(velx%s))), maxval(abs(sol%u(vely%s)))

  end do

  close(unit=11)
  close(unit=23)

! post-processing

  call create_vector ( problem, velocity, physq=physqv )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, gammadot, vec=3 )
  call create_vector ( problem, Dtensor, vec=4 )
  call create_vector ( problem, Ltensor, vec=5 )
  call create_vector ( problem, tauviscous, vec=4 )

  call extract_physvector ( mesh, problem, sol, velocity )

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=8

  call derive_vector ( mesh, problem, gammadot, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Dtensor, elemsub=stokes_D_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Ltensor, elemsub=stokes_gradu_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, tauviscous, elemsub=stokes_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors )

  call stream_function ( mesh, problem_strm, coefficients, sol, vel, sol_strm )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='fountain_flow1b.vtk' )

  call write_vector_vtk ( mesh, problem, filename='fountain_flow1b.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='fountain_flow1b.vtk', dataname='gammadot', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='fountain_flow1b.vtk', &
    dataname='D', vector=Dtensor, append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='fountain_flow1b.vtk', &
    dataname='L', vector=Ltensor, append=.true., symmetric=.false. )

  call write_tensor_vtk ( mesh, problem, filename='fountain_flow1b.vtk', &
    dataname='tauviscous', vector=tauviscous, append=.true. )

  call write_scalar_vtk ( mesh, problem_strm, sysvector=sol_strm, &
    dataname='stream_function', filename='fountain_flow1b.vtk', append=.true. )


! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( vel14 )
  call delete ( velocity, pressure, gammadot, tauviscous, Dtensor, Ltensor )
  call delete ( oldvectors )

  deallocate ( meshcoor_initial, meshcoor_n, meshcoor_nm1 )

  call delete ( coefficients )
  call delete ( velx, vely )

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

! stream function problem

  call delete ( problem_strm )
  call delete ( sol_strm )
  call delete ( vel )

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

! velocity function

  function vel_entry ( nr, xin )

    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: xin
    real(dp) :: vel_entry

    select case ( nr )

    case (1)

!     entry velocity profile (zero flow rate)

      vel_entry = U * ( 1 - 3 * xin(2)**2 / H**2 ) / 2
      vel_entry = vel_entry + Ucorrect

    case default

      write(*,'(/a,i0/)') 'Error vel_entry: wrong function number: ', nr
      stop

    end select

  end function vel_entry


! update the (ALE) mesh based on a predictor for the surface position

  subroutine update_mesh_surface_predictor

!   predictor for surface height
    sol_sf_adv_pred%u = 2._dp*sol_sf_adv_n%u - sol_sf_adv_nm1%u

    meshcoor_nm1 = meshcoor_n
    meshcoor_n = mesh%coor

    disp(:,1) = sol_sf_adv_pred%u - mesh%coor(mesh%curves(14)%nodes,1)
    disp(:,2:ndisp) = 0

!   update nodes of the mesh (only in the first direction)

    call update_mesh_nodes ( mesh, problem_update_mesh, disp(:,1:ndisp), &
      el=mesh_sf_adv%points([1,2]), eg=mesh%points([2,9]) )

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
    tmp = reshape ( sol%u(vel14%s), [2,nnodes] )
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

    write(11,'(i6,es12.4,5es24.16)') step, step*deltat, &
                             max_height, min_height, end_height
    print '(i6,es12.4,5es24.16)', step, step*deltat, &
                          max_height, min_height, end_height

  end subroutine solve_surface_height_corrector

end program fountain_flow1b
