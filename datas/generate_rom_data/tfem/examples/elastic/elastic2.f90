
! Elastic bar in a Stokes fluid. Mooney-Rivlin solid.

module subs_m

  use kind_defs_m
  use tfem_elem_m

  implicit none

! time step
  real(dp) :: deltat = 1.e-2_dp   ! time step

! displacement at tn (full solution vector)
  type(sysvector_t), save :: dispn

! displacement during the iteration within a time step (full solution vector)
  type(sysvector_t), pointer :: dispiter => null()

  integer :: physqdsp = 1, physqvelo = 1, ndf = 9, ndim = 2


contains


! element routine for the constraint between the solid and the fluid

  subroutine elementc ( mesh, problem, constr, elem, node, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemmat2, elemmatadd, &
    elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


    integer :: object, elgrp, elemnr
    real(dp) :: phi(1,ndf), phi2(1,ndf), xr(1,ndim), xr2(1,ndim), u(ndf*ndim)
    real(dp) :: dn(ndim), di(ndim), tmp(ndf,ndim)

    object = problem%constraints(constr)%object

!   reference coordinates in fluid
    xr(1,:) = mesh%objects(object)%refcoor(node,:)
!   reference coordinates in solid
    xr2(1,:) = mesh%objects(object)%refcoor2(node,:)

    call shape_quad_Q2 ( xr, phi )
    call shape_quad_Q2 ( xr2, phi2 )

!   group and element number of the solid
    elgrp = mesh%objects(object)%grpelm2(node,1)
    elemnr = mesh%objects(object)%grpelm2(node,2)

!   get displacement at time tn

    call get_sysvector ( mesh, problem, dispn, elgrp, elemnr, u, &
      physq=[physqdsp] )

    tmp = reshape ( u, [ndf,ndim] )

    dn = matmul(phi2(1,:),tmp)

!   get displacement during the iteration within a time step

    call get_sysvector ( mesh, problem, dispiter, elgrp, elemnr, u, &
      physq=[physqdsp] )

    tmp = reshape ( u, [ndf,ndim] )

    di = matmul(phi2(1,:),tmp)

!   build the constraint:
!          .
!      v - d = 0
!
!   Using first-order backward Euler this becomes:
!
!      v    - ( d -  d ) / delta = 0
!       n+1      n+1  n

    if ( vector ) then

      elemvec = ( di - dn ) / deltat

    end if

    if ( matrix ) then

      elemmat = 0
      elemmat(1,1:9) = phi(1,:)
      elemmat(2,10:18) = phi(1,:)
      elemmat2 = 0
      elemmat2(1,1:9) = -phi2(1,:) / deltat
      elemmat2(2,10:18) = -phi2(1,:) / deltat

    end if

  end subroutine elementc


! element routine for the constraint between the solid and the fluid
! Include term with velocity gradient.

  subroutine elementc2 ( mesh, problem, constr, elem, node, matrix, vector, &
    first, last, coefficients, oldvectors, elemmat, elemmat2, elemmatadd, &
    elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd


    integer :: object, elgrp, elemnr, i
    real(dp) :: phi(1,ndf), phi2(1,ndf), xr(1,ndim), xr2(1,ndim), u(ndf*ndim)
    real(dp) :: dn(ndim), di(ndim), tmp(ndf,ndim)
    real(dp) :: dphi(1,ndf,ndim), dphidx(1,ndf,ndim)
    real(dp) :: x(ndf,ndim), dvdx(ndim,ndim)
    real(dp) :: F(1,ndim,ndim), Finv(1,ndim,ndim), detF(1)

    object = problem%constraints(constr)%object

!   reference coordinates in fluid
    xr(1,:) = mesh%objects(object)%refcoor(node,:)

    call shape_quad_Q2 ( xr, phi, dphi )

!   group and element number of the fluid
    elgrp = mesh%objects(object)%grpelm(node,1)
    elemnr = mesh%objects(object)%grpelm(node,2)

    call get_coordinates ( mesh, elgrp, elemnr, x )

    call isoparametric_deformation ( x, dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

!   get velocity during the iteration within a time step

    call get_sysvector ( mesh, problem, dispiter, elgrp, elemnr, u, &
      physq=[physqvelo] )

    tmp = reshape ( u, [ndf,ndim] )

    dvdx = matmul(transpose(tmp),dphidx(1,:,:))

!   reference coordinates in solid
    xr2(1,:) = mesh%objects(object)%refcoor2(node,:)

    call shape_quad_Q2 ( xr2, phi2 )

!   group and element number of the solid
    elgrp = mesh%objects(object)%grpelm2(node,1)
    elemnr = mesh%objects(object)%grpelm2(node,2)

!   get displacement at time tn

    call get_sysvector ( mesh, problem, dispn, elgrp, elemnr, u, &
      physq=[physqdsp] )

    tmp = reshape ( u, [ndf,ndim] )

    dn = matmul(phi2(1,:),tmp)

!   get displacement during the iteration within a time step

    call get_sysvector ( mesh, problem, dispiter, elgrp, elemnr, u, &
      physq=[physqdsp] )

    tmp = reshape ( u, [ndf,ndim] )

    di = matmul(phi2(1,:),tmp)

!   build the constraint:
!          .
!      v - d = 0
!
!   Using first-order backward Euler this becomes:
!
!      v    - ( d -  d ) / delta = 0
!       n+1      n+1  n

    if ( vector ) then

      elemvec = ( di - dn ) / deltat

    end if

    if ( matrix ) then

      elemmat = 0
      elemmat(1,1:9) = phi(1,:)
      elemmat(2,10:18) = phi(1,:)
      elemmat2 = 0
      elemmat2(1,1:9) = -phi2(1,:) / deltat
      elemmat2(2,10:18) = -phi2(1,:) / deltat
      do i = 1, ndim
        elemmat2(i,1:9)   = elemmat2(i,1:9)   + dvdx(i,1)*phi2(1,:)
        elemmat2(i,10:18) = elemmat2(i,10:18) + dvdx(i,2)*phi2(1,:)
      end do

    end if

  end subroutine elementc2

end module subs_m

module functions_m

  use kind_defs_m

  implicit none

contains

  function flow_profile ( nr, x )

    use kind_defs_m
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: flow_profile

    flow_profile = 1-4*(x(2)-0.5_dp)**2

  end function flow_profile

end module functions_m

program elastic2

  use tfem_m
  use hsl_ma41_m
  use elastic_elements_m
  use stokes_elements_m
  use figplot_m
  use io_utils_m
  use functions_m
  use subs_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities/displacements
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqdisp = 1,      & ! physical quantity nr of the displacements
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3             ! 3x3 integration of quads

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh, mesh1, mesh2, updatemesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(sysvector_t), target :: soln
  type(vector_t) :: displacement, pressure, velocity
  type(plot_options_t) :: plot_options
  type(subscript_t) :: velo, disp, pres, dispobject
  type(subscriptvec_t) :: disp_vec
  type(coefficients_t) :: coefficients_f, coefficients_s
  type(oldvectors_t) :: oldvectors
  type(solver_options_ma41_t) :: solver_options

! variables

  integer :: &
    numtimesteps=10,    & ! number of time steps
    itermax=3,          & ! maximum number of iterations in a time step
    flowreverse=1000,   & ! flow reverse every ... steps
    figoutput=1,        & ! figure output every ... steps
    nx_f = 40,          & ! number of elements of fluid domain x-direction
    ny_f = 10,          & ! number of elements of fluid domain y-direction
    nx_s = 8,           & ! number of elements of solid domain x-direction
    ny_s = 4,           & ! number of elements of solid domain y-direction
    model = 2             ! solid: model=1 (neo-Hookean), 2 (Mooney-Rivlin)

  real(dp) :: &
    ! stop criterium for iterative displacement differences
    epsdd=1e-3_dp, &
    ! stop criterium for iterative pressure differences
    epspd=1e-2_dp, &
    lx_f = 4._dp,     & ! width of fluid domain
    ly_f = 1._dp,     & ! height of fluid domain
    ox_s = 1.95_dp,   & ! position of solid domain
    lx_s = 0.1_dp,    & ! width of solid domain
    ly_s = 0.7_dp,    & ! height of solid domain
    eta = 1.0_dp,     & ! viscosity of the fluid
    rs_up = 1.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_up = 1.0_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    C1 = 1.0e4_dp,    & ! C1 parameter Mooney-Rivlin
    C2 = 0.2e4_dp       ! C2 parameter Mooney-Rivlin

  real(dp), allocatable, dimension(:,:) :: coor_object_orig

  character(len=20) :: filename

! logical for building an object from a curve or a full mesh
  logical :: full_mesh = .false.

  integer :: iter, nno, nnd, step
  real(dp) :: dd, pd


! namelist for input of variables; read from standard input

  namelist /comppar/ numtimesteps, itermax, flowreverse, figoutput, &
    nx_f, ny_f, nx_s, ny_s, model, epsdd, epspd, lx_f, ly_f, ox_s, lx_s, ly_s,&
    eta, C1, C2, deltat, rs_up, is_up

  read ( unit=*, nml=comppar )



! fill fluid coefficients

  call create_coefficients ( coefficients_f, ncoefi=100, ncoefr=50 )

  coefficients_f%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients_f%i(12:) = 0

  coefficients_f%r(1) = eta
  coefficients_f%r(2:) = 0

! fill solid coefficients

  call create_coefficients ( coefficients_s, ncoefi=100, ncoefr=50 )

  coefficients_s%i(1:11) = &
    [ uintpl,   pintpl,      model, 0,     0,  &
      physqdisp, physqpress,     0, 0, gauss,  &
      gauss ]
  coefficients_s%i(12:) = 0

  coefficients_s%r(1:2) = [ C1, C2 ]
  coefficients_s%r(3:) = 0


! create fluid mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx_f
  meshgen_options%ny = ny_f
  meshgen_options%lx = lx_f
  meshgen_options%ly = ly_f

  call quadrilateral2d ( mesh1, meshgen_options )

! create solid mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx_s
  meshgen_options%ny = ny_s
  meshgen_options%ox = ox_s
  meshgen_options%lx = lx_s
  meshgen_options%ly = ly_s

  call quadrilateral2d ( mesh2, meshgen_options )

! merge fluid and solid mesh

  call mesh_merge ( mesh1, mesh2, mesh, nogroupmerge=.true. )

! boundary curve of the solid
  call add_to_mesh ( mesh, curve=[6,7,8] )


! make object for solid-fluid interaction

  if ( full_mesh ) then

!   use full mesh

    call add_to_mesh ( mesh, object='mesh', objectmesh=mesh2, typeofobject=2, &
      excludegroups=[2], excludegroups2=[1], excludecurves=[1] )

  else

!   use boundary curve

    call add_to_mesh ( mesh, object='curve', objectcurve=9, typeofobject=2, &
      excludegroups=[2], excludegroups2=[1], excludecurves=[5] )

  end if

! save original coordinates of the object

  nno = mesh%objects(1)%nnodes
  allocate ( coor_object_orig(nno,2) )
  coor_object_orig = mesh%objects(1)%coor

  call fill_mesh_parts ( mesh )

! do some mesh plotting

  plot_options%fontsize=8
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )
  plot_options%objectpointsize=0.4
  plot_options%objectpointcolor=2
  call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true. )

  plot_options%plotboundary=.false.
  call plot_mesh ( plot_options, mesh, 'mesh1.fig', groups=[1] )
  call plot_mesh ( plot_options, mesh, 'mesh2.fig', groups=[2] )
  plot_options%plotboundary=.true.

  call write_mesh ( mesh, filename='mesh.out' )

! create updatemesh from real mesh:
!   mesh: the original mesh. For the elastic solid the reference configuration.
!   updatemesh: copy of the original mesh, except for the
!               coordinates of the nodal points, which are redefined and
!               updated after each time step.
  call copy ( mesh, updatemesh )
  call fill_mesh_parts ( updatemesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

! fluid group of elements
  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

! solid group of elements
  input_probdef%vec_elementdof(2)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! displacement
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar
                 [9,3] )

  input_probdef%physq = [1,2]

! the fluid boundary
  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )
! the attachment of the solid to the wall
  call define_essential ( mesh, input_probdef, curve1=5, physq=1 )
! the pressure level in one point
  call define_essential ( mesh, input_probdef, point=2, physq=2 )

! constraint between the fluid and the solid

  call define_constraint ( mesh, input_probdef, object=1, &
    discretization='collocation', nodedof=2, physq=1 )

  call problem_definition ( input_probdef, mesh, problem )


! define vector subscripts for direct manipulation of sysvector and vector data

! all velocities of the fluid
  call create_subscript ( mesh, problem, velo, physqarr=[1], groups=[1] )

! all displacements of the solid
  call create_subscript ( mesh, problem, disp, physqarr=[1], groups=[2], &
    fillnodes=.true. )
  nnd = size(disp%nodes) ! number of nodes in displacements
  call create_subscript ( mesh, problem, disp_vec, physq=1, groups=[2] )

! displacements of the object
  if ( full_mesh ) then
!   full mesh
    call create_subscript ( mesh, problem, dispobject, physqarr=[1], &
      groups=[2] )
  else
!   boundary curve
    call create_subscript ( mesh, problem, dispobject, physqarr=[1], &
                            curves=[9], groups=[2], excludecurves=[5] )
  end if

! pressures
  call create_subscript ( mesh, problem, pres, physqarr=[2] )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )  ! Newton-Raphson iteration increment
  call create_sysvector ( problem, rhsd ) ! Right-hand side
  call create_sysvector ( problem, soln ) ! current solution
  call create_sysvector ( problem, dispn ) ! solution at tn used in constraint

  sol%u = 0
  soln%u = 0
  dispn%u = 0

  call create ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => soln ! used in solid element module (pointer copy)

  dispiter => soln ! used in constraint (pointer copy)


! fill velocity solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=2, physq=1, degfd=1, func=flow_profile, funcnr=1 )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=4, physq=1, degfd=1, func=flow_profile, funcnr=1 )


! start time stepping

  stepping: do step = 1, numtimesteps

    write(*,'(/a,i0,a,f10.4/)') 'step = ', step, ' time = ', step * deltat

    write(*,'(a,6x,a,10x,a/)') 'iter', 'dd', 'dp'

!   start Newton-Raphson/Picard iteration within increments

    iter = 0

    iterate: do

      iter = iter + 1

!     build (assemble) matrix and vector from elements

!     revers the flow?

      if ( mod(step+flowreverse/2,flowreverse) == 0 ) then
        sol%u = -sol%u
      end if

!     create system matrix

      call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
      call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
      call finalize_sysmatrix_structure ( sysmatrix )

      call create_sysmatrix_data ( sysmatrix )

!     fluid
      call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
        coefficients=coefficients_f, elgroup1=1 )

!     solid
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=elastic_elem, coefficients=coefficients_s, &
        oldvectors=oldvectors, elgroup1=2, addmatvec=.true., &
        factormat=1._dp/deltat, factorvec=1._dp/deltat )

!     coupling
      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        elemsub=elementc, addmatvec=.true. )


      call check ( sysmatrix )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


      solver_options%real_storage = rs_up
      solver_options%integer_storage = is_up

      call solve_system_ma41 ( sysmatrix, rhsd, sol, &
        solver_options=solver_options )

      call delete ( sysmatrix )

      dd = maxval(abs(sol%u(disp%s))) ! displacement difference
      pd = maxval(abs(sol%u(pres%s)-soln%u(pres%s)))  ! pressure difference

      write(*,'(i4,2es12.4)') iter, dd, pd

      soln%u(velo%s) = sol%u(velo%s) ! copy velocities
      soln%u(disp%s) = soln%u(disp%s) + sol%u(disp%s)  ! update solution
      soln%u(pres%s) = sol%u(pres%s) ! copy pressure

!     update coordinates of the fluid-solid object

      mesh%objects(1)%coor = coor_object_orig + &
         transpose ( reshape ( soln%u(dispobject%s), [2,nno] ) )

!     intersect

      call find_refcoor_objects ( mesh, object1=1, updaterefcoor2=.false. )

!     check convergence

      if ( dd < epsdd .and. pd < epspd ) exit iterate

      if ( iter >= itermax ) then
        write(*,'(/a,i0,a)') 'Maximum number of iterations obtained: ', &
          itermax, '. Going to next step'
        exit iterate
      end if

    end do iterate

!   update coordinates of mesh

    updatemesh%coor(disp%nodes,:) = mesh%coor(disp%nodes,:) + &
           transpose ( reshape (soln%u(disp%s),[2,nnd]) )

    if ( mod(step,figoutput) == 0 ) then

!     plot velocity vector

      write(filename,'(a,i4.4,a)') 'velocity', step, '.fig'

      call plot_vector ( plot_options, updatemesh, problem, filename=filename, &
        physq=1, groups=[1], sysvector=soln )

    end if

    call copy ( soln, dispn )

  end do stepping

! Do some more plotting

  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, displacement, physq=1 )
  call create_vector ( problem, velocity, physq=1 )

! pressure everywhere

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients_f, oldvectors=oldvectors, groups=[1], &
    finalize=.false. )
  call derive_vector ( mesh, problem, pressure, elemsub=elastic_pressure, &
    coefficients=coefficients_s, oldvectors=oldvectors, groups=[2], &
    addvec=.true. )

  call plot_color_contour ( plot_options, updatemesh, problem, &
    'pressure1_contour.fig', groups=[1], vector=pressure )
  call plot_color_contour ( plot_options, updatemesh, problem, &
    'pressure2_contour.fig', groups=[2], vector=pressure )

! displacement vector
  call extract_physvector ( mesh, problem, soln, velocity )
  displacement%u = 0
  displacement%u(disp_vec%s) = velocity%u(disp_vec%s)
  plot_options%scalevector=1
  call plot_vector ( plot_options, mesh, problem, &
    filename='displacement.fig', physq=1, groups=[2], vector=displacement )

! velocity
  plot_options%scalevector=0.1
  call plot_vector ( plot_options, updatemesh, problem, &
    filename='velocity.fig', physq=1, groups=[1], vector=velocity )
  call plot_objects ( plot_options, updatemesh, 'velocity.fig', append=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh, updatemesh )
  call delete ( sol, soln, rhsd, dispn )
  call delete ( displacement, pressure, velocity )
  call delete ( disp, pres, dispobject )
  call delete ( coefficients_f, coefficients_s )
  call delete ( oldvectors )

end program elastic2
