!--------------------------------------------------------------------------
! Example elastic.f90
! <Physical Problem>:
! Elastic bar in a Stokes fluid. Mooney-Rivlin solid.
! <Numerical Approximation>:
! Quadrilateral Taylor-Hood elements for both phases.
! Langrange multipliers for interfacial coupling of the two domains.
! Implicit Euler method is used with constant time step for
! time integration.
!--------------------------------------------------------------------------


!--------------------------------------------------------------------------
!  Module common_definitions
!  Contains all the common parameters for this problem
!--------------------------------------------------------------------------

   module common_definitions

    use kind_defs_m

    implicit none


!   1. constants for the problem definition
    integer, parameter :: &
    uintpl = 8,           & ! Q2 velocities/displacements (Taylor-Hood elements)
    pintpl = 4,           & ! Q1 pressures
    physqvel   = 1,       & ! physical quantity nr of the velocities
    physqdisp  = 1,       & ! physical quantity nr of the displacements
    physqpress = 2,       & ! physical quantity nr of the pressures
    gauss = 3               ! 3x3 integration of quads


!   2. constants for the mesh generation & geometry definition
    integer, parameter :: &
    nx_f  = 100,          & ! number of elements of fluid domain x-direction
    ny_f  = 40,           & ! number of elements of fluid domain y-direction
    nx_s  = 10,           & ! number of elements of solid domain x-direction
    ny_s  = 20              ! number of elements of solid domain y-direction

    real(dp), parameter :: &
    ox_f  = 0.00_dp,      & ! x-position of 1st node of fluid domain
    oy_f  = 0.00_dp,      & ! y-position of 1st node of fluid domain
    lx_f  = 4.00_dp,      & ! width  of fluid domain
    ly_f  = 1.00_dp,      & ! height of fluid domain
    ox_s  = 1.95_dp,      & ! x-position of 1st node of solid domain
    oy_s  = 0.00_dp,      & ! y-position of 1st node of solid domain
    lx_s  = 0.10_dp,      & ! width of solid domain
    ly_s  = 0.70_dp         ! height of solid domain


!   3. properties of both phases
    integer, parameter :: &
    model        = 2     ! solid: model=1 (neo-Hookean), 2 (Mooney-Rivlin)

    real(dp), parameter :: &
    eta = 1.0_dp,          & ! viscosity of the fluid
    rho = 0.0_dp,          & ! density  of the fluid
    C1 = 1.0e2_dp,         & ! C1 parameter Mooney-Rivlin, (=G/2 neo-Hookean)
    C2 = 0.2e2_dp            ! C2 parameter Mooney-Rivlin


!   4. constants for time integration
    integer, parameter :: &
    numtimesteps=200,     & ! number of time steps
    flowreverse =1000       ! flow reverse every ... steps


    real(dp), parameter :: &
    deltat = 1.e-2_dp    ! time step



!   5. constants for Newton-Raphson method
    integer, parameter :: &
    itermax=50           ! maximum number of NR iterations in a time step

    real(dp), parameter :: &
    epsdd = 1e-5_dp,       & ! stop criterium for iterative displacement differences
    epspd = 1e-4_dp          ! stop criterium for iterative pressure differences
    real(dp) :: dd, pd



!   6. constants for the HSL solver
    real(dp), parameter :: &
    rs_up = 1.0_dp,        & ! real_storage gradient-velocity-pressure LU (HSL)
    is_up = 1.0_dp           ! integer_storage gradient-velocity-pressure LU (HSL)


!   7. constants for plotting
    integer, parameter :: &
    figoutput=5            ! figure output every ... steps


!   8. constants for Langrange multiplier (LM)
    character, parameter :: clc_or_wk = 'w'     !  'w' or 'W' for weak imposition & 'c' or 'C' for collocation


!   >>>uncommand one of the following lines for choosing an approximation for the LM.
!      in case of 'collocation' your choice
!    integer :: lintpl_v = 0,  ndflb_v = 1, elementdof_v(3) = (/0,2,0/)    ! constant LM at the middle of element
    integer :: lintpl_v = 1,  ndflb_v = 2, elementdof_v(3) = [2,0,2]    ! piecewise linear continuous    LM
!    integer :: lintpl_v = 2,  ndflb_v = 2, elementdof_v(3) = (/0,4,0/)    ! piecewise linear discontinuous LM



   end module common_definitions


!--------------------------------------------------------------------------
!  Module subs_m
!  Contains subroutines for the fluid-structure constraint
!--------------------------------------------------------------------------

   module subs_m

   use kind_defs_m
   use tfem_elem_m

   use common_definitions, only: deltat, lintpl_v, ndflb_v

   implicit none


!  full solution vector at tn
   type(sysvector_t), pointer :: soln1_loc => null()

!  full solution vector during the iteration within a time step
   type(sysvector_t), pointer :: soln_loc => null()

   integer :: physqdsp = 1, physqvelo = 1, ndf = 9, ndim = 2
   integer :: nodalpb  = 3



   contains


!  element routine for the constraint between the solid and the fluid using galerkin
!  method
   subroutine elementc_weak &
   ( mesh, problem, constr, elem, node, matrix, vector, first, last, coefficients, &
     oldvectors, elemmat, elemmat2, elemmatadd, elemvec, elemvecadd )

!   arguments
    type(mesh_t),          intent(in)  :: mesh
    type(problem_t),       intent(in)  :: problem
    integer,               intent(in)  :: constr, elem, node
    logical,               intent(in)  :: matrix, vector, first, last
    type(coefficients_t),  intent(in)  :: coefficients
    type(oldvectors_t),    intent(in)  :: oldvectors

    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:)   :: elemvec, elemvecadd

!   local variables
    integer, parameter :: ndfb_v  = 3

    integer  :: object, elgrp, elemnr, i, j
    real(dp) :: phi(1,ndf), phi2(1,ndf), xr(1,ndim), xr2(1,ndim), u(ndf*ndim)
    real(dp) :: dn(ndim), di(ndim), tmp(ndf,ndim)

    real(dp) :: psi(1,ndflb_v), x(nodalpb,2)
    real(dp) :: theta_q(1,ndfb_v), dtheta_q(1,ndfb_v,1), dxdxi(1,2), curvel(1)


!   assign to "object" variable the object number
    object = problem%constraints(constr)%object


!   reference coordinates in fluid
!   reference coordinates (-1<ksi<1, -1<eta<1) in quadratic element
!   and shape function in first intersection (object)
    xr(1,:) = mesh%objects(object)%refcoor_int(node,:,elem)
    call shape_quad_Q2 ( xr,  phi  )


!   reference coordinates in solid
!   reference coordinates (-1<ksi<1, -1<eta<1) in quadratic element
!   and shape function in first intersection (object)
    xr2(1,:) = mesh%objects(object)%refcoor2_int(node,:,elem)
    call shape_quad_Q2 ( xr2, phi2 )


!   group and element number of the solid
    elgrp  = mesh%objects(object)%grpelm2_int(node,1,elem)
    elemnr = mesh%objects(object)%grpelm2_int(node,2,elem)

!   get displacement at time tn
    call get_sysvector &
    ( mesh, problem, soln1_loc, elgrp, elemnr, u, physq=[physqdsp] )

    tmp = reshape( u, [ndf,ndim] )
    dn  = matmul(phi2(1,:),tmp)

!   get displacement during the iteration within a time step
    call get_sysvector&
    ( mesh, problem, soln_loc, elgrp, elemnr, u, physq=[physqdsp] )

    tmp = reshape ( u, [ndf,ndim] )
    di  = matmul(phi2(1,:),tmp)


!---------------------------------------------------------------------------------
!   1-D shape function of the Lagrange multiplier

    if ( lintpl_v == 0 ) then
!     constant Lagrange multiplier
      psi = 1
    else if ( lintpl_v == 1 ) then
!     linear Lagrange multiplier
      call shape_line_P1 ( mesh%objects(object)%xig(node:node,1), psi )
    else if ( lintpl_v == 2 ) then
!     discontinuous linear Lagrange multiplier
      psi(1,1) = 1.e0_dp
      psi(1,2) = mesh%objects(object)%xig(node,1)
    end if

!   shape function of the quadratic curve at the integration point
    call shape_line_P2 &
    ( mesh%objects(object)%xig(node:node,1), theta_q, dtheta_q(:,:,1) )


!   compute geometry of element deformed
    call get_coordinates_object ( mesh, elem, x, object )

    call isoparametric_deformation_curve &
    ( x, dtheta_q(:,:,1), dxdxi, curvel )



!--------------------------------------------------------------------------
!   build the constraint:
!          .
!      v - d = 0
!
!   Using first-order backward Euler this becomes:
!
!      v    - ( d -  d ) / delta = 0
!       n+1      n+1  n
!--------------------------------------------------------------------------


    if ( vector ) then
      do i = 1, ndflb_v
      elemvec(i)         =  &
      psi(1,i) *( di(1) - dn(1) )/deltat *curvel(1)* mesh%objects(object)%wg(node)
      elemvec(i+ndflb_v) =  &
      psi(1,i) *( di(2) - dn(2) )/deltat *curvel(1)* mesh%objects(object)%wg(node)
      end do
    end if




    if ( matrix ) then
      elemmat = 0
      do i = 1, ndflb_v
      do j = 1, ndf
      elemmat(i,j)             =  &
      psi(1,i) *phi(1,j) *curvel(1)*mesh%objects(object)%wg(node)
      elemmat(ndflb_v+i,ndf+j) =  &
      psi(1,i) *phi(1,j) *curvel(1)*mesh%objects(object)%wg(node)
      end do
      end do

      elemmat2 = 0
      do i = 1, ndflb_v
      do j = 1, ndf
      elemmat2(i,j)             = &
      - psi(1,i) *phi2(1,j)/deltat *curvel(1)*mesh%objects(object)%wg(node)
      elemmat2(ndflb_v+i,ndf+j) = &
      - psi(1,i) *phi2(1,j)/deltat *curvel(1)*mesh%objects(object)%wg(node)
      end do
      end do
    end if



  end subroutine elementc_weak




!  element routine for the constraint between the solid and the fluid using collocation
!  technique
   subroutine elementc_collocation&
   ( mesh, problem, constr, elem, node, matrix, vector, first, last, coefficients, &
     oldvectors, elemmat, elemmat2, elemmatadd, elemvec, elemvecadd )

!   arguments
    type(mesh_t),          intent(in)     :: mesh
    type(problem_t),       intent(in)     :: problem
    integer,               intent(in)     :: constr, elem, node
    logical,               intent(in)     :: matrix, vector, first, last
    type(coefficients_t),  intent(in)     :: coefficients
    type(oldvectors_t),    intent(in)     :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:)   :: elemvec, elemvecadd

!   local variables
    integer  :: object, elgrp, elemnr
    real(dp) :: phi(1,ndf), phi2(1,ndf), xr(1,ndim), xr2(1,ndim), u(ndf*ndim)
    real(dp) :: dn(ndim), di(ndim), tmp(ndf,ndim)


!   assign to "object" variable the object number
    object = problem%constraints(constr)%object

!   reference coordinates in fluid
    xr(1,:) = mesh%objects(object)%refcoor(node,:)
    call shape_quad_Q2( xr,  phi  )

!   reference coordinates in solid
    xr2(1,:) = mesh%objects(object)%refcoor2(node,:)
    call shape_quad_Q2( xr2, phi2 )

!   group and element number of the solid
    elgrp  = mesh%objects(object)%grpelm2(node,1)
    elemnr = mesh%objects(object)%grpelm2(node,2)

!   get displacement at time tn
    call get_sysvector &
    ( mesh, problem, soln1_loc, elgrp, elemnr, u, physq=[physqdsp] )

    tmp = reshape( u, [ndf,ndim] )
    dn  = matmul(phi2(1,:),tmp)

!   get displacement during the iteration within a time step
    call get_sysvector&
    ( mesh, problem, soln_loc, elgrp, elemnr, u, physq=[physqdsp] )

    tmp = reshape ( u, [ndf,ndim] )
    di  = matmul(phi2(1,:),tmp)


!--------------------------------------------------------------------------
!   build the constraint:
!          .
!      v - d = 0
!
!   Using first-order backward Euler this becomes:
!
!      v    - ( d -  d ) / delta = 0
!       n+1      n+1  n
!--------------------------------------------------------------------------


    if ( vector ) then
      elemvec(:) = ( di(:) - dn(:) ) / deltat
    end if


    if ( matrix ) then
      elemmat(:,:)  = 0.0_dp
      elemmat(1,1:9)    = +phi(1,:)
      elemmat(2,10:18)  = +phi(1,:)
      elemmat2(:,:) = 0.0_dp
      elemmat2(1,1:9)   = -phi2(1,:) / deltat
      elemmat2(2,10:18) = -phi2(1,:) / deltat
    end if


  end subroutine elementc_collocation





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





!--------------------------------------------------------------------------
!--------------------------------------------------------------------------







!--------------------------------------------------------------------------
! Elastic bar in a Stokes fluid. Mooney-Rivlin solid.
!--------------------------------------------------------------------------

    program elastic4

!   core modules
    use tfem_m

!   modules for building equations
    use elastic_elements_m
    use stokes_elements_m

!   modules for solver
    use hsl_ma41_m

!   modules for ploting & io
    use figplot_m
    use io_utils_m

!   user defined modules
    use common_definitions
    use functions_m
    use subs_m


    implicit none


!   definitions based on new types
    type(meshgen_options_t)     :: meshgen_options_f, meshgen_options_s
    type(mesh_t)                :: mesh, mesh_f, mesh_s, updatemesh
    type(input_probdef_t)       :: input_probdef
    type(problem_t)             :: problem
    type(sysmatrix_t)           :: sysmatrix
    type(sysvector_t)           :: rhsd, sol, soltmp
    type(sysvector_t), target   :: soln, soln1
    type(vector_t)              :: displacement, pressure, velocity
    type(plot_options_t)        :: plot_options
    type(subscript_t)           :: velo, disp, pres, dispobject
    type(subscriptvec_t)        :: disp_vec
    type(coefficients_t)        :: coefficients_f, coefficients_s
    type(oldvectors_t)          :: oldvectors, oldvectors_plt
    type(solver_options_ma41_t) :: solver_options



!   definitions based on common types
    real(dp), allocatable, dimension(:,:) :: coor_object_orig
    character(len=20) :: filename

!   local variables
    integer :: time_step
    integer  :: iter, nno, nnd
    real(dp) :: mfac_s



!--------------------------------------------------------------------------
!   fill fluid coefficients
!--------------------------------------------------------------------------

    call create_coefficients &
    ( coefficients_f, ncoefi=600, ncoefr=600 )

    coefficients_f%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
    coefficients_f%i(12:) = 0

    coefficients_f%r(1) = eta
    coefficients_f%r(2:) = 0.0_dp




!--------------------------------------------------------------------------
!   fill solid coefficients
!--------------------------------------------------------------------------

    call create_coefficients&
    ( coefficients_s, ncoefi=100, ncoefr=50 )

    coefficients_s%i(1:11) = &
    [ uintpl,   pintpl,      model, 0,     0,  &
      physqdisp, physqpress,     0, 0, gauss,  &
      gauss ]
    coefficients_s%i(12:) = 0


    coefficients_s%r(1:2) = [ C1, C2 ]
    coefficients_s%r(3:)  = 0.0_dp



!--------------------------------------------------------------------------
!   create fluid mesh
!--------------------------------------------------------------------------

    meshgen_options_f%elshape = 6
    meshgen_options_f%nx = nx_f
    meshgen_options_f%ny = ny_f
    meshgen_options_f%ox = ox_f
    meshgen_options_f%oy = oy_f
    meshgen_options_f%lx = lx_f
    meshgen_options_f%ly = ly_f

    call quadrilateral2d ( mesh_f, meshgen_options_f )  !boundary curves 1-4


!--------------------------------------------------------------------------
!   create solid mesh
!--------------------------------------------------------------------------

    meshgen_options_s%elshape = 6
    meshgen_options_s%nx = nx_s
    meshgen_options_s%ny = ny_s
    meshgen_options_s%ox = ox_s
    meshgen_options_s%oy = oy_s
    meshgen_options_s%lx = lx_s
    meshgen_options_s%ly = ly_s

    call quadrilateral2d ( mesh_s, meshgen_options_s ) !boundary curves 5-8


!--------------------------------------------------------------------------
!   merge fluid and solid mesh
!--------------------------------------------------------------------------

    call mesh_merge&
    ( mesh_f, mesh_s, mesh, nogroupmerge=.true. )





!--------------------------------------------------------------------------
!   boundary curve of the solid
!--------------------------------------------------------------------------

    call add_to_mesh&
    ( mesh, curve=[6,7,8] )             !this is the 9th curve of the whole domain



!--------------------------------------------------------------------------
!   make object for solid-fluid interaction
!--------------------------------------------------------------------------

!     use boundary curve
      select case( clc_or_wk )

      case('c','C')
      call add_to_mesh &
      ( mesh, object='curve', objectcurve=9, typeofobject=2, excludegroups=[2],   &
       excludegroups2=[1], excludecurves=[5] )

      case('w','W')
      call add_to_mesh &
      ( mesh, object='curve', objectcurve=9, typeofobject=2, excludegroups=[2],   &
        excludegroups2=[1], topology=.true., intrule=1, nsubint=10 )

      case default
      write(*,*)'message from main program: incorrect value of "clc_or_wk" variable: ', &
                clc_or_wk
      write(*,*)'program stop'
      stop

      end select




!--------------------------------------------------------------------------
!   save original coordinates of the object
!--------------------------------------------------------------------------

    nno = mesh%objects(1)%nnodes
    allocate ( coor_object_orig(nno,2) )
    coor_object_orig = mesh%objects(1)%coor

    call fill_mesh_parts ( mesh )



!--------------------------------------------------------------------------
!   Plot all meshes & write mesh data in a file
!--------------------------------------------------------------------------

    plot_options%fontsize=8
    call plot_points_curves ( plot_options, mesh, 'curves.fig' )

    call plot_mesh ( plot_options, mesh, 'mesh.fig' )
    plot_options%objectpointsize=0.4
    plot_options%objectpointcolor=2
    call plot_objects ( plot_options, mesh, 'mesh.fig', append=.true. )

    plot_options%plotboundary=.false.
    call plot_mesh ( plot_options, mesh, 'mesh_fluid.fig', groups=[1] )
    call plot_mesh ( plot_options, mesh, 'mesh_solid.fig', groups=[2] )
    plot_options%plotboundary=.true.

    call write_mesh ( mesh, filename='mesh.out' )



!--------------------------------------------------------------------------
!   create updatemesh from real mesh:
!   mesh: the original mesh. For the elastic solid the reference configuration.
!   updatemesh: copy of the original mesh, except for the
!               coordinates of the nodal points, which are redefined and
!               updated after each time step.
!--------------------------------------------------------------------------

    call copy ( mesh, updatemesh )
    call fill_mesh_parts ( updatemesh )



!--------------------------------------------------------------------------
!   problem definition
!--------------------------------------------------------------------------

    call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

!   fluid group of elements
    input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

!   solid group of elements
    input_probdef%vec_elementdof(2)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! displacement
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar
                 [9,3] )

    input_probdef%physq = [1,2]





!--------------------------------------------------------------------------
!   Boundaries
!--------------------------------------------------------------------------

!   the fluid boundary
    call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )

!   the attachment of the solid to the wall
    call define_essential ( mesh, input_probdef, curve1=5, physq=1 )

!   the pressure level in one point
    call define_essential ( mesh, input_probdef, point=2, physq=2 )





!--------------------------------------------------------------------------
!   constraint between the fluid and the solid
!--------------------------------------------------------------------------

    select case( clc_or_wk )

    case('c','C')
    call define_constraint &
    ( mesh, input_probdef, object=1, discretization='collocation', &
      nodedof=2, physq=1 )

    case('w','W')
    call define_constraint &
    ( mesh, input_probdef, object=1, discretization='weak', physq=1, &
      elementdof=elementdof_v )

    case default
      write(*,'(/a,a/)') 'Error: wrong value clc_or_wk: ', clc_or_wk
      stop
    end select


    call problem_definition &
    ( input_probdef, mesh, problem )



!--------------------------------------------------------------------------
!   define vector subscripts for direct manipulation of sysvector and vector data
!--------------------------------------------------------------------------

!   all velocities of the fluid
    call create_subscript &
    ( mesh, problem, velo, physqarr=[1], groups=[1] )

!   all displacements of the solid
    call create_subscript &
    ( mesh, problem, disp, physqarr=[1], groups=[2], fillnodes=.true. )

    nnd = size(disp%nodes) ! number of nodes in displacements
    call create_subscript &
    ( mesh, problem, disp_vec, physq=1, groups=[2] )



!    interfacial boundary curve
     select case( clc_or_wk )
     case('c','C')
     call create_subscript &
     ( mesh, problem, dispobject, physqarr=[1], curves=[9], &
       groups=[2], excludecurves=[5] )
     case('w','W')
     call create_subscript &
     ( mesh, problem, dispobject, physqarr=[1], curves=[9], &
       groups=[2] )
     case default
       write(*,'(/a,a/)') 'Error: wrong value clc_or_wk: ', clc_or_wk
       stop
     end select



!   pressures
    call create_subscript &
    ( mesh, problem, pres, physqarr=[2] )



!--------------------------------------------------------------------------
!   create system vectors (solution and right-hand side)
!--------------------------------------------------------------------------

    call create_sysvector( problem, rhsd )       ! Right-hand side of Ax=b system
    call create_sysvector( problem, sol  )       ! Newton-Raphson iteration increment (mixed dx/x)

    call create_sysvector( problem, soln  )      ! solution vector (x) at each Newton-Rapshon cycle
    call create_sysvector( problem, soln1 )      ! solution vector (x) at previous time step
    call create_sysvector( problem, soltmp)      ! working array



!   initialize system vectors
    sol%u         = 0.0_dp
    soln%u        = 0.0_dp
    soln1%u       = 0.0_dp


!   create oldvectors
    call create ( oldvectors, nsysvec=1 )
    oldvectors%s(1)%p     => soln                ! used in building subroutines

    call create ( oldvectors_plt, nsysvec=1 )
    oldvectors_plt%s(1)%p => soln                ! used in ploting subroutines



!   pointer copy for contraint subroutines
    soln_loc  => soln
    soln1_loc => soln1



!--------------------------------------------------------------------------
!   fill velocity solution vector with essential boundary conditions
!--------------------------------------------------------------------------

    call fill_sysvector&
    ( mesh, problem, sol, curve1=2, physq=1, degfd=1, func=flow_profile, funcnr=1 )

    call fill_sysvector&
    ( mesh, problem, sol, curve1=4, physq=1, degfd=1, func=flow_profile, funcnr=1 )

    call fill_sysvector&
    ( mesh, problem, sol, point=1, physq=2, value=0._dp )


!--------------------------------------------------------------------------
!    preprocessing for plotting
!--------------------------------------------------------------------------

    call create_vector &
    ( problem, pressure, vec=3 )

    call create_vector &
    ( problem, displacement, physq=1 )

    call create_vector &
    ( problem, velocity, physq=1 )



!--------------------------------------------------------------------------
!   start time stepping
!--------------------------------------------------------------------------

    time_integration: do time_step = 1, numtimesteps

    write(*,'(/a,i0,a,f10.4/)') 'step = ', time_step, ' time = ', time_step * deltat
    write(*,'(a,6x,a,10x,a/)') 'iter', 'dd', 'dp'



!   start Newton-Raphson/Picard iteration within increments
    iter = 0


    nr_iterations: do

      iter = iter + 1


!     build (assemble) matrix and vector from elements

!     revers the flow?
      if ( mod(time_step+flowreverse/2,flowreverse) == 0 ) then
        write(*,*)'the flow direction is reversing'
        write(*,*)'newton raphson may not converge at this time step '
        sol%u = -sol%u
      end if


!     create system matrix
      call create_sysmatrix_structure_base       ( sysmatrix, mesh, problem )
      call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
      call finalize_sysmatrix_structure          ( sysmatrix )
      call create_sysmatrix_data                 ( sysmatrix )



!     fluid domain
      call build_system &
      ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, coefficients=coefficients_f,&
        elgroup1=1 )



!     assign value to scaling factor
      mfac_s = 1._dp/deltat       ! it comes from the kinematic constraint

!     solid domain
      call build_system &
      ( mesh, problem, sysmatrix, rhsd, elemsub=elastic_elem, coefficients=coefficients_s, &
        oldvectors=oldvectors, elgroup1=2, addmatvec=.true., factormat=mfac_s, factorvec=mfac_s )


!     coupling
      select case( clc_or_wk )
      case('c','C')
       call build_system_constraint&
       ( mesh, problem, sysmatrix, rhsd, elemsub=elementc_collocation, addmatvec=.true. )
      case('w','W')
       call build_system_constraint&
       ( mesh, problem, sysmatrix, rhsd, elemsub=elementc_weak, addmatvec=.true. )
      case default
        write(*,'(/a,a/)') 'Error: wrong value clc_or_wk: ', clc_or_wk
        stop
      end select


!     check the consistency of matrix
      call check ( sysmatrix )


!     eliminate prescribed degrees of freedom
      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )


!     call direct solver
      solver_options%real_storage    = rs_up
      solver_options%integer_storage = is_up

      call solve_system_ma41 &
      ( sysmatrix, rhsd, sol, solver_options=solver_options )


!     delete the structure
      call delete ( sysmatrix )


!     calculate solution increment during NR cycles
      dd = maxval(abs(sol%u(disp%s)))                 ! displacement difference
      pd = maxval(abs(sol%u(pres%s)-soln%u(pres%s)))  ! pressure difference

      write(*, '(i4,2es12.4)') iter, dd, pd


!     update solution
      soltmp%u(:)      = 0.0_dp                          ! clean working array
      soltmp%u(disp%s) = soln%u(disp%s) + sol%u(disp%s)  ! update position vector (x+dx)
      call copy(sol, soln)                               ! copy whole solution vector
      soln%u(disp%s) = soltmp%u(disp%s)                  ! finally position vector



!     update coordinates of the fluid-solid object
      mesh%objects(1)%coor = coor_object_orig + &
         transpose ( reshape ( soln%u(dispobject%s), [2,nno] ) )


!     find coordinates of intersection points
      call find_refcoor_objects &
      ( mesh, object1=1, updaterefcoor2=.false. )


!     check convergence
      if ( dd < epsdd .and. pd < epspd ) exit nr_iterations

      if ( iter >= itermax ) then
      write(*,'(/a,i0,a)') 'Maximum number of iterations obtained: ', itermax, '. Going to next step'
      exit nr_iterations
      end if

    end do nr_iterations




!--------------------------------------------------------------------------
!   update coordinates of mesh
!--------------------------------------------------------------------------

    updatemesh%coor(disp%nodes,:) = &
    mesh%coor(disp%nodes,:) + transpose ( reshape (soln%u(disp%s),[2,nnd]) )


!--------------------------------------------------------------------------
!   create figures & plots
!--------------------------------------------------------------------------

    if ( mod(time_step,figoutput) == 0 ) then

!     plot velocity vector
      write(filename,'(a,i4.4,a)') 'velocity', time_step, '.fig'
      call plot_vector&
      ( plot_options, updatemesh, problem, filename=filename, physq=1, &
        groups=[1], sysvector=soln )

!     calculate pressure in fluid & solid domains & plot it
      call derive_vector &
      ( mesh, problem, pressure, elemsub=stokes_pressure, coefficients=coefficients_f, &
        oldvectors=oldvectors_plt, groups=[1], finalize=.false. )

      call derive_vector &
      ( mesh, problem, pressure, elemsub=elastic_pressure, coefficients=coefficients_s, &
        oldvectors=oldvectors_plt, groups=[2], addvec=.true. )


      write(filename,'(a,i4.4,a)') 'pressure_fd', time_step, '.fig'
      call plot_color_contour &
      ( plot_options, updatemesh, problem, filename=filename, groups=[1], &
        vector=pressure )


      write(filename,'(a,i4.4,a)') 'pressure_sd', time_step, '.fig'
      call plot_color_contour &
      ( plot_options, updatemesh, problem, filename=filename, groups=[2], &
        vector=pressure )

    end if



!--------------------------------------------------------------------------
!   update solutions vectors
!--------------------------------------------------------------------------

    call copy( soln,  soln1 )


  end do time_integration






!--------------------------------------------------------------------------
!   Do some more plotting
!--------------------------------------------------------------------------

!   Displacement vector
    call extract_physvector ( mesh, problem, soln, velocity )

    displacement%u = 0
    displacement%u(disp_vec%s) = velocity%u(disp_vec%s)
    plot_options%scalevector=1
    call plot_vector &
    ( plot_options, mesh, problem, filename='displacement.fig', &
      physq=1, groups=[2], vector=displacement )





!--------------------------------------------------------------------------
!   Delete all data including all allocated memory
!--------------------------------------------------------------------------

    call delete ( problem )
    call delete ( input_probdef )
    call delete ( mesh, updatemesh )
    call delete ( sol, soln, soln1, rhsd )
    call delete ( soltmp )
    call delete ( displacement, pressure, velocity )
    call delete ( disp, pres, dispobject )
    call delete ( coefficients_f, coefficients_s )
    call delete ( oldvectors, oldvectors_plt )


end program elastic4
