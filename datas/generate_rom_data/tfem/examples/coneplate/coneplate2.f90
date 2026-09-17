! Viscoelastic cone-plate problem
! (2D) domain with 3D velocities (axisymmetric).
! Cylindrical coordinates with circumferential velocity included (swirl).
! DEVSS-G/SUPG
! first-order time integration to steady state

program coneplate2

  use tfem_m
  use hsl_ma41_m
  use hsl_ma57_m
  use viscoelastic_elements_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,     & ! P2 velocities
    pintpl = 2,     & ! P1 pressures
    gintpl = 2,     & ! P1 gradients
    cintpl = 2,     & ! P1 conformation
    physqg = 1,  & ! physical quantity nr of the gradients
    physqv = 2,  & ! physical quantity nr of the velocities
    physqp = 3,  & ! physical quantity nr of the pressures
    gauss  = 6,  & ! 6-point Gauss integration of triangles
    gaussb = 3,  & ! 3-point integration of boundary elements
    ncompc = 6,  & ! number of conformation tensor components
    nmodes = 1,  & ! number of modes
    startm = 501   ! start of material model data

  real(dp), parameter :: &
    Ro = 1.0_dp,      & ! outer radius (of fluid contact point on the plate)
    theta_cp = 0.1_dp,& ! Cone angle
    dx_Ro = 0.01_dp,  & ! Element size at the fluid outside surface
    dx_0 = 0.001_dp     ! Element size at the center (r=0)

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefc
  type(problem_t), target :: problem, problemc
  type(sysmatrix_t) :: sysmatrix, sysmatrixc
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, gammadot, Dtensor, Ltensor
  type(vector_t) :: tauviscous, tauviscoelastic, ctensor
  type(oldvectors_t) :: oldvectors, oldvectors_ve
  type(coefficients_t) :: coefficients
  type(subscript_t) :: velz, velr, velt, cval
  type(subscriptvec_t) :: szz, szr, szt, srr, srt, stt
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc
  type(lu_ma57_t) :: lu_u
  type(lu_ma41_t) :: luc
  type(solver_options_ma57_t) :: solver_options_u
  type(solver_options_ma41_t) :: solver_options_c

! variables

  integer :: &
    timeint = 1,         & ! (first-order) time integration
    numtimesteps = 2500, & ! number of time steps
    logc = 0,            & ! standard scheme or log transformation
    model = 3              ! 2: Oldroyd-B 3: Giesekus

  real(dp) :: &
    eta_s = 0.1_dp,    & ! solvent viscosity
    eta_p = 1.0_dp,    & ! polymer viscosity
    lambda = 2.0_dp,   & ! relaxation time
    alphapar = 0.1_dp, & ! alpha parameter in the Giesekus model
    deltat = 5.e-3_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    omega = 0.1_dp,    & ! Angular velocity of the cone
    rs_gup = 3.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 3.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 3.4_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 3.4_dp       ! integer_storage for conformation LU (HSL)

  integer :: icomp, step, i, k, npar
  integer :: vertices(3) = [1,3,5]
  real(dp) :: alpha, G

  character(len=30) :: filename

! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus

  if ( model == 2 ) then
    npar = 2 ! Oldroyd-B
  else if ( model == 3 ) then
    npar = 3 ! Giesekus
  end if

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+npar*nmodes )

  coefficients%i = &
    [ uintpl, pintpl, 0, 0,      gintpl, &
      physqv, physqp, 0, physqg, gauss,  &
      gaussb,  cintpl, 0, 0,      0,      &
      0,      0,      model, nmodes, startm, &
      logc,   timeint, ( 0, i = 23, 150 )  &
    ]

  coefficients%i(23) = 1  ! cylindrical coordinate system (axisymmetrical)
  coefficients%i(67) = 1  ! 3D velocity (swirl)

  coefficients%r(1:502) = &
    [ eta_s,    0._dp,   0._dp, alpha, 0._dp, &
      0.0_dp, 0._dp,  deltat, beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G, lambda  &
    ]

  if ( model == 3 ) then
    coefficients%r(503) = alphapar ! Giesekus
  end if

! create mesh

  call create_mesh

  call fill_mesh_parts ( mesh )

! write mesh for plotting

  call write_mesh_vtk ( mesh, 'mesh.vtk' )

  do i = 1, mesh%ncurves
    write(filename,'(a,i4.4,a)') 'curve_',i,'.vtk'
    call write_geometry_vtk ( mesh, curve=i, filename=filename )
  end do


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=6, nphysq=3 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 6  ! gradients
  input_probdef%vec_elementdof(1)%a(:,2) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5) = 6  ! symmetric tensor
  input_probdef%vec_elementdof(1)%a(:,6) = 9  ! unsymmetric tensor

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

! Cone
  call define_essential ( mesh, input_probdef, curve1=2, physq=physqv )

! Plate
  call define_essential ( mesh, input_probdef, curve1=1, physq=physqv )

  call problem_definition ( input_probdef, mesh, problem )

! create vector subscripts for the velocity components for monitoring

  call create_subscript ( mesh, problem, velz, physqarr=[physqv], degfd=1 )
  call create_subscript ( mesh, problem, velr, physqarr=[physqv], degfd=2 )
  call create_subscript ( mesh, problem, velt, physqarr=[physqv], degfd=3 )

! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=3, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a(:,1) = 0
  input_probdefc%vec_elementdof(1)%a(vertices,1) = 1  ! c
  input_probdefc%vec_elementdof(1)%a(:,2) = 1         ! scalar for plotting
  input_probdefc%vec_elementdof(1)%a(:,3) = 6         ! symmetric tensor

  input_probdefc%physq = [1]
  input_probdefc%probnr = 2

  call problem_definition ( input_probdefc, mesh, problemc )

! create a vector subscript for the conformation tensor without Lagr. multipl.
! for post processing the data

  call create_subscript ( mesh, problemc, cval, physqarr=[1] )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0
  call fill_sysvector ( mesh, problem, sol, &
    curve1=2, physq=physqv, degfd=3, func=func, funcnr=1 )
  call fill_sysvector ( mesh, problem, sol, curve1=1, &
    physq=physqv, degfd=1, value=1e-20_dp )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector for gradient/velocity/pressure problem

! stokes velocity/pressure
  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, physqrow=[2,3], physqcol=[2,3], &
    coefficients=coefficients )

! DEVSS-G
  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=devssg_elem, addmatvec=.true., &
    physqrow=[1,2], physqcol=[1,2], coefficients=coefficients )

! set to zero off-diagonal blocks gradient-pressure
  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., physqrow=[1], physqcol=[3], zeromatvec=.true. )
  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., physqrow=[3], physqcol=[1], zeromatvec=.true. )

  call check ( sysmatrix )

! solve initial gradient/velocity/pressure consistent with initial conditions
! and keep decomposition of the matrix
! This generates a Stokes profile.

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  solver_options_u%real_storage=rs_gup
  solver_options_u%integer_storage=is_gup

  call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
    solver_options=solver_options_u  )

! create system vectors (solution and right-hand side) for conformation and
! initialize vectors to zero stress

  call create ( problemc, solc, rhsc )

  if ( logc == 0 ) then ! standard
    solc(1,1)%u = 1
    solc(2,1)%u = 0
    solc(3,1)%u = 0
    solc(4,1)%u = 1
    solc(5,1)%u = 0
    solc(6,1)%u = 1
  else if ( logc == 1 ) then ! log scheme
    do i = 1, ncompc
      solc(i,1)%u = 0
    end do
  end if

! create system matrix for conformation problem

  call create_sysmatrix_structure_base ( sysmatrixc, mesh, problemc )
  call finalize_sysmatrix_structure ( sysmatrixc )

  call create_sysmatrix_data ( sysmatrixc )

! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=1, nsysvec2=1, nprob=2 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc

! open monitor file

  open ( unit=13, file='out', recl=300 )


! time stepping

  do step = 1, numtimesteps


!   build (assemble) matrix and vector for conformation problem

    call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
      elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
      coefficients=coefficients )

    call check ( sysmatrixc )


!   solve conformation and keep LU decomposition in the loop over components

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do icomp = 1, ncompc
      call solve_system_ma41 ( sysmatrixc, rhsc(icomp,1), solc(icomp,1), luc, &
        solver_options=solver_options_c  )
    end do

    call delete ( luc )  ! remove LU decomposition and rebuild next time step


!   build (assemble) vector for gradient/velocity/pressure problem

    call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=rhs_divtau, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      buildmatrix = .false., coefficients=coefficients )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma57 ( sysmatrix, rhsd, sol, lu_u, &
      solver_options=solver_options_u  )


!   write monitor data

    write(unit=13,fmt=*) &
      step, step*deltat, maxval(solc(1,1)%u(cval%s)), &
      maxval(abs(sol%u(velz%s))), maxval(abs(sol%u(velr%s))), &
      maxval(abs(sol%u(velt%s)))

  end do


! close monitor data file

  close(unit=13)

! post-processing

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  call create_vector ( problem, velocity, physq=physqv )
  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, gammadot, vec=4 )
  call create_vector ( problem, Dtensor, vec=5 )
  call create_vector ( problem, Ltensor, vec=6 )
  call create_vector ( problem, tauviscous, vec=5 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

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

  call derive_vector ( mesh, problem, tauviscous, elemsub=stokes_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='cp2.vtk' )

  call write_vector_vtk ( mesh, problem, filename='cp2.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true., &
    assume3D=.true. )

  call write_vector_vtk ( mesh, problem, filename='cp2.vtk', &
    dataname='velocity_vector_planar', vector=velocity, append=.true., &
    degfd=[1,2] )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='cp2.vtk', dataname='gammadot', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='cp2.vtk', &
    dataname='D', vector=Dtensor, append=.true., assume3D=.true. )

  call write_tensor_vtk ( mesh, problem, filename='cp2.vtk', &
    dataname='L', vector=Ltensor, append=.true., assume3D=.true., &
    symmetric=.false. )

  call write_tensor_vtk ( mesh, problem, filename='cp2.vtk', &
    dataname='tauviscous', vector=tauviscous, append=.true., assume3D=.true. )

  call create_vector ( problemc, tauviscoelastic, vec=3 )

  call derive_vector ( mesh, problemc, tauviscoelastic, &
    elemsub=deriv_viscoelastic_stress_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemc, filename='cp2.vtk', &
    dataname='tauviscoelastic', vector=tauviscoelastic, append=.true., &
    assume3D=.true. )

  call create_vector ( problemc, ctensor, vec=3 )

  call derive_vector ( mesh, problemc, ctensor, &
    elemsub=deriv_conformation_tensor,&
    coefficients=coefficients, oldvectors=oldvectors_ve )

  call write_tensor_vtk ( mesh, problemc, filename='cp2.vtk', &
    dataname='c', vector=ctensor, append=.true., &
    assume3D=.true. )

  call create ( mesh, problem, szz, vec=5, degfd=1, curves=[1], &
    fillnodes=.true. )
  call create ( mesh, problem, szr, vec=5, degfd=2, curves=[1] )
  call create ( mesh, problem, szt, vec=5, degfd=3, curves=[1] )
  call create ( mesh, problem, srr, vec=5, degfd=4, curves=[1] )
  call create ( mesh, problem, srt, vec=5, degfd=5, curves=[1] )
  call create ( mesh, problem, stt, vec=5, degfd=6, curves=[1] )

  open( unit=84, file='stress_on_plate2.out')

  do i = 1, size(szz%s)
     k = szz%nodes(i)
     write(84,'(15es16.8)') mesh%coor(k,1), mesh%coor(k,2), pressure%u(k), &
       tauviscous%u(szz%s(i)), tauviscous%u(szr%s(i)), &
       tauviscous%u(szt%s(i)), tauviscous%u(srr%s(i)), &
       tauviscous%u(srt%s(i)), tauviscous%u(stt%s(i)), &
       tauviscoelastic%u(szz%s(i)), tauviscoelastic%u(szr%s(i)), &
       tauviscoelastic%u(szt%s(i)), tauviscoelastic%u(srr%s(i)), &
       tauviscoelastic%u(srt%s(i)), tauviscoelastic%u(stt%s(i))
  end do

  close(unit=84)

  call create ( mesh, problem, szz, vec=5, degfd=1, curves=[3], &
    fillnodes=.true. )
  call create ( mesh, problem, szr, vec=5, degfd=2, curves=[3] )
  call create ( mesh, problem, szt, vec=5, degfd=3, curves=[3] )
  call create ( mesh, problem, srr, vec=5, degfd=4, curves=[3] )
  call create ( mesh, problem, srt, vec=5, degfd=5, curves=[3] )
  call create ( mesh, problem, stt, vec=5, degfd=6, curves=[3] )

  open( unit=85, file='stress_on_surface2.out')

  do i = 1, size(szz%s)
     k = szz%nodes(i)
     write(85,'(15es16.8)') mesh%coor(k,1), mesh%coor(k,2), pressure%u(k), &
       tauviscous%u(szz%s(i)), tauviscous%u(szr%s(i)), &
       tauviscous%u(szt%s(i)), tauviscous%u(srr%s(i)), &
       tauviscous%u(srt%s(i)), tauviscous%u(stt%s(i)), &
       tauviscoelastic%u(szz%s(i)), tauviscoelastic%u(szr%s(i)), &
       tauviscoelastic%u(szt%s(i)), tauviscoelastic%u(srr%s(i)), &
       tauviscoelastic%u(srt%s(i)), tauviscoelastic%u(stt%s(i))
  end do

  close(unit=85)

! delete all data including all allocated memory

  call delete ( mesh )
  call delete ( problem )
  call delete ( input_probdef )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( velocity, pressure, gammadot, Dtensor, Ltensor, tauviscous, &
                tauviscoelastic, ctensor )
  call delete ( oldvectors, oldvectors_ve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( lu_u )
  call delete ( solc, rhsc )

  call delete ( coefficients )
  call delete ( velz, velr, velt, cval )
  call delete ( szz, szr, szt, srr, srt, stt )

contains


! create the mesh for the cone-plate geometry using gmsh

  subroutine create_mesh

    use math_defs_m

    open ( unit=25, file='ConePlate2d.geo' )

    write ( 25, '(1X,A,F18.14,A)' ) 'z[1] = ', Ro*cos(pi/2-theta_cp), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'r[1] = ', Ro*sin(pi/2-theta_cp), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'theta_cp = ', theta_cp, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_Ro = ', dx_Ro, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_0 = ', dx_0, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'Ro = ', Ro, ';'
    write ( 25, '(/1x,a)' ) 'Include "ConePlate2d.igo";'

    close ( 25 )

    call execute_command_line ( 'gmsh -2 -order 2 -algo front2d &
       & -o ConePlate2d.msh ConePlate2d.geo > outputmesh.out' )

!   read mesh generated by gmsh
    call read_mesh_gmsh ( mesh, filename='ConePlate2d.msh', ndim=2, &
      physgeom=.true., sortphys=.true. )

  end subroutine create_mesh


! function for velocity on wall

  function func ( nr, x )

    integer, intent(in) ::  nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    select case(nr)
      case(1)
        func = omega*x(2)
      case default
        write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
        stop
    end select

  end function func

end program coneplate2
