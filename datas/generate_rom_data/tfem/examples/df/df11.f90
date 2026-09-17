! Startup of flow around a cylinder for an integral model.
! Problem 1: Flow around a confined cylinder.
!   Periodical boundary conditions.
!   Stress-implicit formulation of the momentum balance.
! Problem 2: Deformation fields computed as a function of time, assuming the
!   flow starts instantaneously at t=0 and was at rest for t<0.
! Problem 3: Stress tensor projection on the interpolation space for F.
! Implicit DG with P1 or P2 interpolation for tau-convection.
! Extrapolated BDF2 time discretization.
! Some CPU savings for equidistant grids in tau, since matrix is the same for
! all intervals.

program df11

  use tfem_m
  use hsl_ma57_m
  use hsl_ma41_m
  use stokes_elements_m
  use devss_elements_m
  use dfm_elements_m
  use io_utils_m
  use figplot_m
  use timer_m

  implicit none

! constants

  logical, parameter :: &
    equidistant = .false. ! Use equidistant grid in age (tau).

  integer, parameter :: &
    uintpl = 6,         & ! P2 velocities
    pintpl = 2,         & ! P1 pressures
    gintpl = 2,         & ! P1 gradients
    physqgrad = 1,      & ! physical quantity nr of the gradients
    physqvel = 2,       & ! physical quantity nr of the velocities
    physqpress = 3,     & ! physical quantity nr of the pressures
    gauss = 6,          & ! 6 point integration of triangles
    gaussb = 3            ! 3 point integration of boundary elements

  real(dp), parameter :: &
    eta_s = 0._dp,  & ! solvent viscosity
    H = 2._dp,      & ! (half-) height of the channel (must match mesh)
    U = 0.5_dp,     & ! average velocity in the channel
    flowrate = H*U  ! flow rate in (half) the channel

  integer, parameter :: &
    type_of_model = 0 , & ! 0: separable Rivlin-Sawyers,
                          ! 1: non-separable RS (not yet available)
    spectrum = 0,       & ! 0: discrete spectrum
                          ! 1: Mittag-Leffler (not yet available)
    nummodes = 1,       & ! number of modes in the discrete spectrum
    dampfunc = 1          ! damping function: 1: Lodge, 2: McKinley, 3: PSM

  real(dp), parameter :: &
    G = 1._dp,            & ! modulus
    lambda = 1._dp,       & ! relaxation time
    a = 0.01_dp,          & ! parameter in McKinley damping function
    alpha_PSM = 2.5e4_dp, & ! alpha parameter in PSM damping function
    beta_PSM = 0.25_dp      ! beta parameter in PSM damping function

  integer, parameter :: &
    fintpl = 2,         & ! P1 interpolation in space for F
    fintpltau = 1,      & ! P1 interpolation in age (tau)
!    fintpltau = 2,      & ! P2 interpolation in age (tau)
    nintvaltau = 50,    & ! number of age intervals
    ncompf = 4,         & ! number of components of deformation tensor F
    ncompt = 4,         & ! number of stress components
    ninttau = fintpltau+1, & ! number of fields within an interval
    numtimesteps = 1, & ! number of time steps
    ns = 3                ! number of samples

  real(dp), parameter :: &
    deltat = min(0.1_dp,0.25_dp/U),  & ! time step
    dtau1 = 6e-2_dp,     & ! length of first interval in age (tau) direction:
    tauc = 15.0_dp,      & ! maximum age (cutoff) (only for equidistant=.false.)
    beta = 1,            & ! SUPG factor
    rs_gup = 1.0_dp,     & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.4_dp,     & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_f = 1.0_dp,       & ! real_storage for the deformation LU (HSL)
    is_f = 1.7_dp          ! integer_storage for the deformation LU (HSL)

! definitions

  type(mesh_t) :: mesh, mesh1
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(subscript_t) :: velx, vely
  type(vector_t) :: velocity, pressure, vorticity, gammadot

  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options_u

  type(input_probdef_t) :: input_probdef_f
  type(problem_t), target :: problem_f
  type(sysmatrix_t) :: sysmatrix_f
  type(subscript_t) :: fval
  type(oldvectors_t) :: oldvectors_f

  type(sysvector_t), dimension(ncompf,nintvaltau), target :: &
    sol_f, soln_f, solnm1_f
  type(sysvector_t), dimension(ncompf) :: rhs_f

  type(vector_t) :: Ffield, Ffieldtensor, Bfieldtensor, stress_tensor

  type(lu_ma41_t) :: lu_f
  type(solver_options_ma41_t) :: solver_options_f

  type(input_probdef_t) :: input_probdef_s_proj
  type(problem_t), target :: problem_s_proj
  type(sysmatrix_t) :: sysmatrix_s_proj
  type(sysvector_t), dimension(ncompt), target :: sol_s_proj
  type(sysvector_t), dimension(ncompt) :: rhs_s_proj
  type(lu_ma57_t) :: lu_s_proj

  type(sample_t) :: sample

  logical :: buildmatrix
  character(len=20) :: fname
  integer :: icomp, step, i, j, is
  integer :: vertices(3) = [1,3,5]
  real(dp) :: x(0:nintvaltau), xs(ns,2), alpha
  real(dp) :: tau_nodal_points(ninttau,nintvaltau)
  real(dp) :: fields(ncompf,ninttau,nintvaltau)


! set some parameters

  alpha = G * lambda   ! DEVSS parameter

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=400, ncoefr=350, &
    ncoefra=[nintvaltau+1,nummodes,nummodes] )

  coefficients%i = 0
  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gaussb ]

  coefficients%i(351:353) = [ fintpl, fintpltau, nintvaltau ]
  coefficients%i(354) = 1  ! first-order time integration
  coefficients%i(355) = 1  ! SUPG
  coefficients%i(356) = 0  ! G from DEVSS
  coefficients%i(357) = 0  ! velocity gradient = G

  coefficients%i(363) = type_of_model
  coefficients%i(364) = spectrum
  coefficients%i(365) = nummodes
  coefficients%i(366) = dampfunc

  coefficients%r = 0
  coefficients%r(1) = eta_s
  coefficients%r(4) = alpha
  coefficients%r(6) = flowrate

  coefficients%r(301) = deltat
  coefficients%r(302) = beta

  coefficients%r(303) = a
  coefficients%r(304) = alpha_PSM
  coefficients%r(305) = beta_PSM

  if ( equidistant ) then
!   equidistant tau grid
    coefficients%ra(1)%a = [ (i,i=0,nintvaltau) ] * dtau1
    write(*,'(/a,g0.4/)') ' equidistant grid; cutoff time = ', nintvaltau*dtau1
  else
!   non-equidistant tau grid
    call distribute_elements ( nintvaltau, x, ratio=7, factor=dtau1/tauc )
    coefficients%ra(1)%a = tauc * x(0:nintvaltau)
    write(*,'(/a,g0.4/)') ' non-equidistant grid; last interval = ', &
       coefficients%ra(1)%a(nintvaltau+1)-coefficients%ra(1)%a(nintvaltau)
  end if

  coefficients%ra(2)%a = [ G ]
  coefficients%ra(3)%a = [ lambda ]


! read mesh

  call read_mesh_gmsh ( mesh1, filename='confined_cylinder1.msh', ndim=2 )
  call mesh_convert ( mesh1, mesh, remove_isolated_nodes=.true. )
  call delete ( mesh1 )

! cylinder
  call add_to_mesh ( mesh, curve=[-2] ) ! curve 9
! top
  call add_to_mesh ( mesh, curve=[5,6,7] ) ! curve 10
! centerline+cylinder
  call add_to_mesh ( mesh, curve=[1,-2,3] ) ! curve 11
! connect curves for periodical bc in DG
  call add_to_mesh ( mesh, curve=[-8] )     ! curve 12

! one object for sampling the fields in a points in space

  xs(1,:) =  [ 1.01_dp, 0.0_dp ]
  xs(2,:) =  [ 1.1_dp, 0.0_dp ]
  xs(3,:) =  [ 2.0_dp, 0.0_dp ]

  call add_to_mesh ( mesh, object='coordinates', coor=xs )

  call fill_mesh_parts ( mesh )

  call printinfo ( mesh, printlevel=1 )

! plot mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  plot_options%fontsize = 10

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=3 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 4  ! G
  input_probdef%vec_elementdof(1)%a(:,2) = 2         ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

! define essential boundaries

! cylinder
  call define_essential ( mesh, input_probdef, curve1=9, physq=physqvel )
! top (wall)
  call define_essential ( mesh, input_probdef, curve1=10, physq=physqvel )
! bottom (center line)
  call define_essential ( mesh, input_probdef, curve1=1, &
    physq=physqvel, degfd=[0,1] )
  call define_essential ( mesh, input_probdef, curve1=3, &
    physq=physqvel, degfd=[0,1] )
! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=12, nglobalc=1 )

! constraints for periodical boundary conditions

! velocities (use weak connection)
  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=4, curve2=12, discretization='weak', &
    elementdof=[2,0,2] )

! gradients (use collocation)
  call define_constraint ( mesh, input_probdef, &
    physq=physqgrad, curve1=4, curve2=12, discretization='collocation' )

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for the velocity for post processing

  call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
  call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! problem definition deformation tensor

  call create_input_probdef ( mesh, input_probdef_f, nvec=3, nphysq=1 )

  input_probdef_f%vec_elementdof(1)%a(vertices,1) = ninttau ! f
  input_probdef_f%vec_elementdof(1)%a(:,2) = 1 ! scalar for plotting
  input_probdef_f%vec_elementdof(1)%a(:,3) = 4 ! tensor for plotting

  input_probdef_f%physq = [1]
  input_probdef_f%probnr = 2

! constraint for periodical boundary conditions of the deformation tensor
  call define_constraint ( mesh, input_probdef_f, curve1=4, curve2=12, &
    discretization='collocation' )

  call problem_definition ( input_probdef_f, mesh, problem_f )

! create a vector subscript for the deformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problem_f, fval, physqarr=[1] )

! create system vectors (solution and right-hand side) for F and initialize

  call create ( problem_f, sol_f, soln_f, solnm1_f )
  call create ( problem_f, rhs_f )

! initial solution

  do j = 1, nintvaltau
    soln_f(1,j)%u = 1 ! initial F_xx
    soln_f(2,j)%u = 0 ! initial F_xy
    soln_f(3,j)%u = 0 ! initial F_yx
    soln_f(4,j)%u = 1 ! initial F_yy
  end do


! create system matrix for df fields problem

  call create_sysmatrix_structure_base ( sysmatrix_f, mesh, problem_f )
  call create_sysmatrix_structure_constraint ( sysmatrix_f, mesh, problem_f )
  call finalize_sysmatrix_structure ( sysmatrix_f )

  call create_sysmatrix_data ( sysmatrix_f )


! create the structure oldvectors_f

  call create_oldvectors ( oldvectors_f, nsysvec=1, nsysvec1=1, nsysvec2=3, &
    nprob=3 )

! store solution vectors and problem structures

  oldvectors_f%s(1)%p => sol
  oldvectors_f%s2(1)%p => sol_f
  oldvectors_f%s2(2)%p => soln_f
  oldvectors_f%s2(3)%p => solnm1_f
  oldvectors_f%p(1)%p => problem
  oldvectors_f%p(2)%p => problem_f


! problem definition for L2-projection of the stress tensor

  call create_input_probdef ( mesh, input_probdef_s_proj, nvec=1, nphysq=1 )

  input_probdef_s_proj%vec_elementdof(1)%a = &
      reshape ( [ 1,0,1,0,1,0 ], &
                   [6,1] )

  input_probdef_s_proj%physq = [1]
  input_probdef_s_proj%probnr = 3

  call problem_definition ( input_probdef_s_proj, mesh, problem_s_proj )

  call create ( problem_s_proj, sol_s_proj, rhs_s_proj )

! store solution vectors and problem structures

  oldvectors_f%s1(1)%p => sol_s_proj
  oldvectors_f%p(3)%p => problem_s_proj

! create and build system matrix for projection problem
! NOTE matrix remains constant and needs to be build once.

  call create_sysmatrix_structure ( sysmatrix_s_proj, mesh, problem_s_proj, &
    symmetric=.true. )
  call create_sysmatrix_data ( sysmatrix_s_proj )

  call build_system ( mesh, problem_s_proj, sysmatrix_s_proj, &
    msysvector=rhs_s_proj, elemsub=stress_projection_elem, &
    oldvectors=oldvectors_f, coefficients=coefficients, &
    buildvector=.false. )

  call check ( sysmatrix_s_proj )


! open monitor files

  open ( unit=13, file='out', recl=300 )
  open ( unit=14, file='out_stress', recl=300 )


! time stepping

  call tic

  do step = 1, numtimesteps

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG

!   compute stress projection

    call solve_stress_projection

!   build implicit terms of CE with rhs in momentum balance

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_intmodel_elem, &
      oldvectors=oldvectors_f, physqrow=[physqvel], physqcol=[physqvel], &
      addmatvec=.true., coefficients=coefficients )

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )


!   deformation fields

    if ( step == 1 ) then
      coefficients%i(354) = 1 ! Euler first time step
    else
      coefficients%i(354) = 2 ! second-order BDF2 extrapolated
    end if

    do j = 1, nintvaltau

      if ( equidistant .and. j > 1 ) then
        buildmatrix = .false.
      else
        buildmatrix = .true.
      end if

!     build (assemble) matrix and vector for the deformation problem

      coefficients%i(362) = j  ! interval number

      call build_system ( mesh, problem_f, sysmatrix_f, msysvector=rhs_f, &
        elemsub=dfm_elem2, oldvectors=oldvectors_f, &
        coefficients=coefficients, buildmatrix=buildmatrix )

!     periodical condition on df fields problem
      call build_system_constraint ( mesh, problem_f, sysmatrix_f, &
        msysvector=rhs_f, elemsub=stokes_constr_node_conn, &
        addmatvec=.true., buildmatrix=buildmatrix )

      call check ( sysmatrix_f )


!     solve fields and keep LU decomposition in the loop over components

      solver_options_f%real_storage=rs_f
      solver_options_f%integer_storage=is_f

      do icomp = 1, ncompf
        call solve_system_ma41 ( sysmatrix_f, rhs_f(icomp), &
          sol_f(icomp,j), lu_f, solver_options=solver_options_f  )
      end do

      if ( .not. equidistant ) then
!       remove LU decomposition and rebuild next interval
        call delete ( lu_f )
      end if

    end do

    if ( equidistant ) then
!     remove LU decomposition and rebuild next time step
      call delete ( lu_f )
    end if

    call copy ( soln_f, solnm1_f )
    call copy ( sol_f, soln_f )

!   write monitor data

    j = nintvaltau

    write(unit=13,fmt=*) &
      step, step*deltat, &
      maxval(sol_f(1,j)%u(fval%s)), minval(sol_f(1,j)%u(fval%s)), &
      maxval(sol_f(2,j)%u(fval%s)), minval(sol_f(2,j)%u(fval%s)), &
      maxval(sol_f(3,j)%u(fval%s)), minval(sol_f(3,j)%u(fval%s)), &
      maxval(sol_f(4,j)%u(fval%s)), minval(sol_f(4,j)%u(fval%s))

    write(unit=14,fmt=*) &
      step, step*deltat, &
      maxval(sol_s_proj(1)%u), minval(sol_s_proj(1)%u), &
      maxval(sol_s_proj(2)%u), minval(sol_s_proj(2)%u), &
      maxval(sol_s_proj(3)%u), minval(sol_s_proj(3)%u), &
      maxval(sol_s_proj(4)%u), minval(sol_s_proj(4)%u), &
      maxval(abs(sol%u(velx%s))), maxval(abs(sol%u(vely%s)))

    print *, 'step = ', step, ' time = ', step*deltat
    call toc ( 'one step' )

  end do

  call toc ( 'all steps' )


! post-processing

  call create_vector ( problem, velocity, physq=2 )
  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, vorticity, vec=4 )
  call create_vector ( problem, gammadot, vec=4 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! derive vectors

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=8

  call derive_vector ( mesh, problem, gammadot, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )


! write to a fig file for plotting

  plot_options%scalevector=0.05
  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  plot_options%numlevels=9
  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

  call plot_color_contour ( plot_options, mesh, problem, &
    'gammadot_contour.fig', vector=gammadot )

  call write_vector_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='velocity', vector=velocity )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='pressure', vector=pressure, append=.true. )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='vorticity', vector=vorticity, append=.true. )
  call write_scalar_vtk ( mesh, problem, filename='sol.vtk', &
    dataname='gammadot', vector=gammadot, append=.true. )

  call printtofile ( mesh, problem, filename='velocity_cl.out', curve=11, &
    vector=velocity )
  call printtofile ( mesh, problem, filename='pressure_cl.out', curve=11, &
    vector=pressure )


! plot fields

  call create_vector ( problem_f, Ffield, vec=2 )
  call create_vector ( problem_f, Ffieldtensor, vec=3 )
  call create_vector ( problem_f, Bfieldtensor, vec=3 )

  coefficients%i(359)=nintvaltau ! interval
  coefficients%i(360)=ninttau ! field

  call derive_vector ( mesh, problem_f, Ffieldtensor, &
    elemsub=dfm_deriv_field_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_f )

  coefficients%i(358)=1 ! detF

  call derive_vector ( mesh, problem_f, Ffield, &
    elemsub=dfm_deriv, &
    coefficients=coefficients, oldvectors=oldvectors_f )

  coefficients%i(358)=2 ! B=F.F^T

  call derive_vector ( mesh, problem_f, Bfieldtensor, &
    elemsub=dfm_deriv, &
    coefficients=coefficients, oldvectors=oldvectors_f )

  call plot_color_contour ( plot_options, mesh, problem_f, &
    'Fxx.fig', vector=Ffieldtensor, degfd=1 )
  call plot_color_contour ( plot_options, mesh, problem_f, &
    'Fxy.fig', vector=Ffieldtensor, degfd=2 )
  call plot_color_contour ( plot_options, mesh, problem_f, &
    'Fyx.fig', vector=Ffieldtensor, degfd=3 )
  call plot_color_contour ( plot_options, mesh, problem_f, &
    'Fyy.fig', vector=Ffieldtensor, degfd=4 )
  call plot_color_contour ( plot_options, mesh, problem_f, &
    'detF.fig', vector=Ffield )

  call write_tensor_vtk ( mesh, problem_f, filename='f.vtk', &
    dataname='F-field', vector=Ffieldtensor, symmetric=.false. )
  call write_scalar_vtk ( mesh, problem_f, filename='f.vtk', &
    dataname='detF', vector=Ffield, append=.true. )
  call write_tensor_vtk ( mesh, problem_f, filename='f.vtk', &
    dataname='B-field', vector=Bfieldtensor, symmetric=.false., append=.true. )


! output values of field as a function of tau for points in space

  call dfm_tau_nodal_points ( coefficients, tau_nodal_points )

  call fill_sample ( mesh, problem_f, sample, &
    ndegfd=ncompf*ninttau*nintvaltau, object=1, elemsub=dfm_sample, &
    coefficients=coefficients, oldvectors=oldvectors_f )

  do is = 1, ns

    write(fname,'(a,i0)') 'outs', is
    open ( unit=15, file=fname, recl=300 )

    fields = reshape ( sample%u(is,:), [ ncompf,ninttau,nintvaltau ] )

    do j = 1, nintvaltau
      do i = 1, ninttau
        write(unit=15,fmt=*) tau_nodal_points(i,j), &
          fields(1,i,j), fields(2,i,j), fields(3,i,j), fields(4,i,j)
      end do
      write(unit=15,fmt=*)
    end do

    close(15)

  end do


! post processing of the stress tensor

  call create_vector ( problem_f, stress_tensor, vec=3 )

  call derive_vector ( mesh, problem_f, stress_tensor, &
    elemsub=dfm_deriv_stress_tensor, &
    coefficients=coefficients, oldvectors=oldvectors_f )

  call write_tensor_vtk ( mesh, problem_f, filename='stress.vtk', &
    dataname='stress', vector=stress_tensor, symmetric=.true. )

  call plot_color_contour ( plot_options, mesh, problem_f, &
    'stress_xx.fig', vector=stress_tensor, degfd=1 )
  call plot_color_contour ( plot_options, mesh, problem_f, &
    'stress_xy.fig', vector=stress_tensor, degfd=2 )
  call plot_color_contour ( plot_options, mesh, problem_f, &
    'stress_yy.fig', vector=stress_tensor, degfd=3 )

  if ( dampfunc /= 1 ) then
    call write_scalar_vtk ( mesh, problem_f, filename='stress.vtk', &
      dataname='stress_zz', vector=stress_tensor, degfd=4, append=.true. )
    call plot_color_contour ( plot_options, mesh, problem_f, &
      'stress_zz.fig', vector=stress_tensor, degfd=4 )
  end if

  call printtofile ( mesh, problem_f, filename='stress_cl.out', curve=11, &
    vector=stress_tensor )


! post processing of the projected stress tensor

  call solve_stress_projection

  call derive_vector ( mesh, problem_f, stress_tensor, &
    elemsub=dfm_deriv_stress_tensor_proj, &
    coefficients=coefficients, oldvectors=oldvectors_f )

  call write_tensor_vtk ( mesh, problem_f, filename='stress_proj.vtk', &
    dataname='stress_proj', vector=stress_tensor, symmetric=.true. )

  call plot_color_contour ( plot_options, mesh, problem_f, &
    'stress_xx_proj.fig', vector=stress_tensor, degfd=1 )
  call plot_color_contour ( plot_options, mesh, problem_f, &
    'stress_xy_proj.fig', vector=stress_tensor, degfd=2 )
  call plot_color_contour ( plot_options, mesh, problem_f, &
    'stress_yy_proj.fig', vector=stress_tensor, degfd=3 )
  if ( dampfunc /= 1 ) then
    call write_scalar_vtk ( mesh, problem_f, filename='stress_proj.vtk', &
      dataname='stress_zz_proj', vector=stress_tensor, degfd=4, append=.true. )
    call plot_color_contour ( plot_options, mesh, problem_f, &
      'stress_zz_proj.fig', vector=stress_tensor, degfd=4 )
  end if

  call printtofile ( mesh, problem_f, filename='stress_proj_cl.out', curve=11, &
    vector=stress_tensor )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity, gammadot )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

  call delete ( problem_f )
  call delete ( input_probdef_f )
  call delete ( sol_f, soln_f, solnm1_f )
  call delete ( rhs_f )
  call delete ( sysmatrix_f )
  call delete ( oldvectors_f )
  call delete ( Ffield, Ffieldtensor, Bfieldtensor, stress_tensor )
  call delete ( fval )

  call delete ( sysmatrix_s_proj )
  call delete ( problem_s_proj )
  call delete ( input_probdef_s_proj )
  call delete ( sol_s_proj, rhs_s_proj )
  call delete ( lu_s_proj )

contains

  subroutine build_vpG

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients, &
      physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress] )

!   DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_elem, coefficients=coefficients, addmatvec=.true., &
      physqrow=[physqgrad,physqvel], physqcol=[physqgrad,physqvel] )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqgrad], physqcol=[physqpress], &
      zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqpress], physqcol=[physqgrad], &
      zeromatvec=.true. )

!   flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients )

!   periodical condition on velocities

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=stokes_constr_elem_conn, addmatvec=.true., &
      coefficients=coefficients )

!   periodical condition on gradients

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=3, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients )

  end subroutine build_vpG

  subroutine solve_stress_projection

    type(solver_options_ma57_t) :: solver_options_ma57

    integer :: i

!   build vector only (matrix is constant)

    call build_system ( mesh, problem_s_proj, sysmatrix_s_proj, &
      msysvector=rhs_s_proj, elemsub=stress_projection_elem, &
      oldvectors=oldvectors_f, coefficients=coefficients, &
      buildmatrix=.false. )

    ! MA57 solver storage
    solver_options_ma57%integer_storage = 1.3
    solver_options_ma57%real_storage    = 1.3

!   LU decomposition is done in the first call only

    do i = 1, ncompt

      call add_effect_of_essential_to_rhs ( problem_s_proj, sysmatrix_s_proj, &
         sol_s_proj(i), rhs_s_proj(i) )

      call solve_system_ma57 ( sysmatrix_s_proj, rhs_s_proj(i), &
         sol_s_proj(i), lu_s_proj, solver_options=solver_options_ma57 )

    end do

  end subroutine solve_stress_projection

end program df11
