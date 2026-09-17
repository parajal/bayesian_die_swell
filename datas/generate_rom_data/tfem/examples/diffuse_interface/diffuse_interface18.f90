! axisymmetric diffuse-interface problem of a droplet in a hyperbolic
! flow field
! Stokes problem using the stress form or one of the potential forms for
! the surface tension
! variable viscosity
! second order time integration
! adaptive meshing using Gmsh (triangular elements)
! Newton-Raphson iterations for the non-linear term in the
! Cahn-Hilliard equation
! This example is similar to diffuse_interface17, but now the flow and
! diffuse-interface equations are solved fully-implicit (in one system)

program diffuse_interface18

  use tfem_m
  use math_defs_m
  use hsl_ma41_m
  use hsl_ma57_m
  use generalized_stokes_elements_m
  use diffuse_interface_elements_generic_m
  use diffuse_interface_elements_m
  use functions_m
  use subs_m
  use io_utils_m
  use projection_elements_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 6,         & ! P2 velocities
    pintpl = 2,         & ! P1 pressures
    physqc = 1,         & ! physical quantity nr of c
    physqmu = 2,        & ! physical quantity nr of mu
    physqvel = 3,       & ! physical quantity nr of the velocities
    physqpress = 4,     & ! physical quantity nr of the pressures
    gauss = 6,          & ! 6-point Gauss integration of triangles
    gaussb = 3,         & ! 3-point integration of boundary elements
    gauss_proj = 6,     & ! integration rule for the projection problem
    inttype_proj = 3,   & ! standard Gauss-Legendre (numerical table)
    coorsys = 1,        & ! coordinate system
    vtkevery = 1,       & ! vtk file every vtkevery steps. -1: means none
    ncompv_P2 = 5,      & ! number of components for P2 projection: c at n
                          ! and n-1, and velocities and mu at n
    numtimesteps = 100, & ! number of time steps
    capillary_form = 1, & ! 1: potential form mugradc
                          ! 2: potential form -cgradmu
                          ! 3: stress form
    itermax = 50,       & ! maximum number of Newton iterations
    n_init = 1000         ! initial number of refinement points on DI


    logical :: implicit_viscosity = .false. ! use viscosity in iteration
                                            ! (only Picard)


! variables

  integer :: &
    timeint1 = 1,       & ! first-order Euler
    timeint2 = 2          ! second-order Gear

  real(dp) :: &
    eta1 = 1.0_dp,        & ! viscosity of the ambient fluid (c=1)
    eta2 = 2.0_dp,        & ! viscosity of the droplet (c=-1)
    rho = 100.0_dp,       & ! density in diffuse-interface method
    alpha = 1.0_dp,       & ! parameter in the diffuse-interface method
    beta = 1.0_dp,        & ! parameter in the diffuse-interface method
    Mcoef = 1.e-2_dp,     & ! parameter in the diffuse-interface method
    kappa = 1.e-4_dp,     & ! parameter in the diffuse-interface method
    epsdi = 1.e-6_dp,     & ! accuracy in Newton iteration
    elong_rate = 1._dp,   & ! elongation rate
    drop_z0 = 0.5_dp,     & ! initial z-coordinate of center of drop
    drop_r = 0.2_dp,      & ! initial radius of the drop
    time = 0.0_dp,        & ! initial time
    h_coarse = 0.1_dp,    & ! element spacing on the external boundaries
    h_fine = 0.01_dp,     & ! element spacing around the interface
    d_min = 0.02_dp,      & ! minimal width of refined area
    d_max = 0.2_dp,       & ! maximal width of refined area
    oz = 0.0_dp,          & ! z-coordinate of origin of mesh
    or = 0.0_dp,          & ! r-coordinate of origin of mesh
    lz = 1.0_dp,          & ! length of mesh in z-direction
    lr = 1.0_dp,          & ! length of mesh in r-direction
    deltat = 0.04_dp,     & ! time step
    dist_threshold = 0.34_dp,&! elements for which any (abs(c)<dist_threshold)
                              ! are allowed to contain the interface
    rs_up_di = 1.4_dp,    & ! real_storage (HSL)
    is_up_di = 1.6_dp       ! integer_storage LU (HSL)


! definitions for DI and up problems

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol_np1, sol_n, sol_nm1, sol_i
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(solver_options_ma41_t) :: solver_options_ma41
  type(subscript_t) :: velx, vely, pval
  type(subscript_t) :: cval, muval
  type(vector_t), target :: composition
  type(elvector_t), target :: viscosity

! temp definitions used to store old definitions before projection

  type(mesh_t), target :: mesh_temp
  type(sysvector_t) :: sol_n_temp, sol_nm1_temp
  type(problem_t), target :: problem_temp

! type definition for the refinement fields

  type(refinement_fields_t) :: refinement_fields

! variables

  integer :: step=0, ipost=0, iter

  real(dp) :: cmax, cdiff, c_B, surf_ten
  real(dp) :: velxmax, velymax, velxdiff, velydiff, mumax, mudiff
  real(dp), allocatable, dimension(:,:) :: refinement_coor
  logical, allocatable, dimension(:) :: elems_allowed, elems_refine

  character(len=30) :: filename
  logical :: remeshing=.false.


! pass value to functions_m

  xi = sqrt(kappa/alpha)
  allocate(c(2))
  c(1) = drop_z0
  c(2) = 0.0_dp
  r = drop_r
  eta1l = eta1
  eta2l = eta2


! determine the capillary number

  c_B = sqrt(alpha/beta)
  surf_ten = rho*(2*sqrt(2.0_dp)/3) * (kappa * c_B**2) / xi

  print *,'Capillary number = ',eta2*elong_rate*drop_r / surf_ten


! for the initial mesh, the refinement points have to be given manually

  allocate ( refinement_coor(n_init,2) )

  call coors_circle ( refinement_coor )

  call add_refinement_field ( refinement_fields, coor=refinement_coor, &
    distmin=d_min, distmax=d_max, dx_fine=h_fine, dx_coarse=h_coarse )

  deallocate ( refinement_coor )


! fill coefficients for the flow problem

  call create_coefficients ( coefficients, ncoefi=300, ncoefr=250 )

  coefficients%i(:)  = 0
  coefficients%i(1:15)  = [ uintpl,       pintpl,  0,  0,     0, &
                            physqvel, physqpress,  0,  0, gauss, &
                            gaussb,            0,  0,  0,     0 ]
  coefficients%i(23) = coorsys

  coefficients%r = 0
  coefficients%i(251) = 3 ! eta is function of composition

! fill coefficients for the diffuse interface problem

  coefficients%r(101:106)  = [ deltat,  Mcoef,  alpha,  beta,  kappa, rho ]
  coefficients%i(151) = timeint1 ! Euler integration of the di equations
  coefficients%i(152) = 2        ! Newton-Raphson iteration for di equation
  coefficients%i(153) = 2        ! use the iteration velocity u_i in the
                                 ! DI-convection term, instead of the
                                 ! prediction uhat =2*u_nm1 - un
                                 ! (diffuse_interface_elem1)


! generate and read the refined mesh using gmsh

  call generate_read_mesh_refined


! define the up and di problem, create (sys)vectors, subscripts, oldvectors
! and initialize sysvectors

  call define_problems_create_vectors


! sample the initial c-field and then find the initial velocity-field

  call find_velocity_init


! find the elements that are allowed to contain the interface

  call find_elems_allowed


! write VTK files

  call postprocessing ( mesh, sol_np1, ipost=ipost )
  ipost = ipost + 1


! copy old values

  call copy ( sol_np1, sol_n )


! start time integration

  do step = 1, numtimesteps

    time = time + deltat

    write ( *, '(1X,A)' ) '****************************************'
    write ( *, '(1X,A,I0)' ) 'step = ', step


!   change time integration to Gear after the first step

    if ( step > 1 ) coefficients%i(151) = timeint2


!   remeshing criterion

    remeshing = .false.


!   check if the interface is in allowed elements

    call check_elems_allowed


    if ( remeshing ) then


      print *,'Remeshing and projection...'


!     find refinement points

      call find_refinement_points


!     copy current definitions to temporary definitions used for the projection
!     problem

      call copy ( mesh, mesh_temp )
      call fill_mesh_parts ( mesh_temp )
      call copy ( sol_n, sol_n_temp )
      call copy ( sol_nm1, sol_nm1_temp )
      call problem_definition ( input_probdef, mesh, problem_temp )


!     delete the old problems

      call delete_old_problems


!     generate new mesh

      call generate_read_mesh_refined


!     define the up and di problem, create (sys)vectors, subscripts, oldvectors
!     and initialize sysvectors

      call define_problems_create_vectors


!     project the old solutions onto the new mesh and fill in the values

      call project_old_solutions_P2


!     delete the temporary definitions

      call delete ( mesh_temp )
      call delete ( sol_n_temp )
      call delete ( sol_nm1_temp )
      call delete ( problem_temp )


    end if


!   build and solve the di-up system

    call build_solve_up_di


!   find the elements that are allowed to contain the interface

    if ( remeshing ) call find_elems_allowed


!   write VTK files

    if ( vtkevery > 0 ) then
      if ( mod(step,vtkevery) == 0 ) then
        call postprocessing ( mesh, sol_np1, ipost=ipost )
        ipost = ipost + 1
      end if
    end if


!   copy old values

    call copy ( sol_n, sol_nm1 )
    call copy ( sol_np1, sol_n )


  end do


! delete all data including all allocated memory

  call delete_old_problems


contains


! generate and read mesh

  subroutine generate_read_mesh_refined

!   write the mesh parameters
    open ( unit=25, file='mesh.geo' )
    write ( 25, '(1X,A,F18.14,A)' ) 'ox = ', oz, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'oy = ', or, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'lx = ', lz, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'ly = ', lr, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'h_coarse = ', h_coarse, ';'

!   write the refinement fields (the file needs to be open for writing)
    call write_refinement_fields ( refinement_fields, 'mesh.geo' )

    write ( 25, '(/1x,a)' ) 'Include "mesh_quad_refined.igo";'
    close ( 25 )

!   refinement points served their purpose
    call delete_refinement_fields ( refinement_fields )

!   generate the mesh using gmsh
    call execute_command_line ( 'gmsh -2 -order 2 -algo del2d -o mesh.msh &
        &mesh.geo > outputmesh.out' )

!   read mesh generated by gmsh
    call read_mesh_gmsh ( mesh, filename='mesh.msh', ndim=2, &
      physgeom=.true. )

    call fill_mesh_parts ( mesh )

!   allocate some arrays
    allocate ( elems_allowed(mesh%nelem) )
    allocate ( elems_refine(mesh%nelem) )

  end subroutine generate_read_mesh_refined


! find the elements that are allowed to contain the interface

  subroutine find_elems_allowed

    type(vector_t) :: c
    integer :: nod(mesh%element(1)%numnod), elem

    call create_vector ( problem, c, physq=physqc )
    call extract_physvector (mesh, problem, sol_np1, c )

    elems_allowed = .false.

    do elem = 1, mesh%nelem

!     nodal points
      nod = mesh%topology(1)%a(:,elem)

      if ( any ( abs(c%u(nod)) < dist_threshold ) ) then
!       close enough to interface
        elems_allowed(elem) = .true.
      else
!       too far from interface
        cycle
      end if

    end do

    call delete(c)

  end subroutine find_elems_allowed


! check if the interface is in allowed elements

  subroutine check_elems_allowed

    type(vector_t) :: c

    integer :: nod(mesh%element(1)%numnod), elem

    call create_vector ( problem, c, physq=physqc )
    call extract_physvector (mesh, problem, sol_np1, c )

    do elem = 1, mesh%nelem

!     nodal points
      nod = mesh%topology(1)%a(:,elem)

      if ( all ( c%u(nod) > 0.0 ) .or. all ( c%u(nod) < 0.0 ) ) then
!       all nodes far from the interface
        cycle
      else ! element contains interface
        if ( elems_allowed(elem) .eqv. .false. ) then
          remeshing = .true.
          print *,' *******Interface in non-allowed elements: Remeshing*******'
          print *,' '
          exit
        end if
      end if

    end do

    call delete(c)

  end subroutine check_elems_allowed


! find the refinement points of the mesh

  subroutine find_refinement_points

    type(vector_t) :: c
    integer :: nod(mesh%element(1)%numnod), elem, ii
    real(dp) :: c_elem(mesh%element(1)%numnod)

    call create_vector ( problem, c, physq=physqc )
    call extract_physvector (mesh, problem, sol_np1, c )

!   first traversing loop: find elements containing interface
    do elem = 1, mesh%nelem

      nod = mesh%topology(1)%a(:,elem)  ! nodal points

!     only add refinement points for elements containing the interface
      elems_refine(elem) = &
        .not. ( all ( c%u(nod) > 0.0 ) .or. all ( c%u(nod) < 0.0 ) )

     end do

    allocate ( refinement_coor(count(elems_refine),2) )

!   second traversing loop: find the refinement points
    ii = 1
    do elem = 1, mesh%nelem

      if ( elems_refine(elem) ) then

        nod = mesh%topology(1)%a(:,elem)  ! nodal points

!       add a refinement point for the node closes to the interface
        c_elem = abs(c%u(nod))
        refinement_coor(ii,:) = mesh%coor(nod(minloc(c_elem,1)),:)
        ii = ii + 1

      end if

    end do

!   add the refinement points to a refinement field
    call add_refinement_field ( refinement_fields, coor=refinement_coor, &
      distmin=d_min, distmax=d_max, dx_fine=h_fine, dx_coarse=h_coarse )

    call delete(c)

    deallocate ( refinement_coor )

  end subroutine find_refinement_points


! start the problems

  subroutine define_problems_create_vectors

!   problem definition for the stokes problem
    call create_input_probdef ( mesh, input_probdef, nvec=7, nphysq=4 )

    input_probdef%vec_elementdof(1)%a =   &
        reshape ( [1,1,1,1,1,1,    &  ! c
                   1,1,1,1,1,1,    &  ! mu
                   2,2,2,2,2,2,    &  ! velocity
                   1,0,1,0,1,0,    &  ! pressure
                   1,1,1,1,1,1,    &  ! scalar
                   2,2,2,2,2,2,    &  ! coordinates
                   ncompv_P2,ncompv_P2,ncompv_P2,   &
                   ncompv_P2,ncompv_P2,ncompv_P2 ], &
                         [6,7] )

    input_probdef%physq = [1,2,3,4]
    input_probdef%probnr = 1

!   Dirichlet boundary conditions

    call define_essential ( mesh, input_probdef, curve1=1, degfd=[0,1], &
      physq=physqvel )
    call define_essential ( mesh, input_probdef, curve1=2, curve2=4, &
      physq=physqvel )
    call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

    call problem_definition ( input_probdef, mesh, problem )

!   create system vectors for the up problem
    call create_sysvector ( problem, sol_np1, sol_n, sol_nm1, rhsd )
    call create_sysvector ( problem, sol_i )

!   create vector to store the viscosity
    call create_elvector ( mesh, viscosity, nreal1d=1 ) ! eta data

!   create vector for the composition (used to determine the local viscosity)
    call create_vector ( problem, composition, vec=1 )

!   create subscripts for the velocity and pressure
    call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1 )
    call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2 )
    call create_subscript ( mesh, problem, pval, physqarr=[physqpress] )

!   create subscripts for c and mu
    call create_subscript ( mesh, problem, cval, physqarr=[physqc] )
    call create_subscript ( mesh, problem, muval, physqarr=[physqmu] )

!   create system matrix for the up problem
    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric=.false. )
    call finalize_sysmatrix_structure ( sysmatrix )
    call create_sysmatrix_data ( sysmatrix )

!   create oldvectors
!   NOTE: two problems are supplied for backwards compatibility with the
!   diffuse_interface_elem1 routine
    call create_oldvectors ( oldvectors, nsysvec=4, nprob=2, nelvec=1, &
      nvec=1 )
      oldvectors%s(1)%p => null()
      oldvectors%s(2)%p => sol_n
      oldvectors%s(3)%p => sol_np1
      oldvectors%s(4)%p => sol_nm1
      oldvectors%p(1)%p => problem
      oldvectors%p(2)%p => problem
      oldvectors%e(1)%p => viscosity
      oldvectors%v(1)%p => composition

!   initialize all vectors to zero
    sol_nm1%u = 0.0_dp
    sol_n%u = 0.0_dp
    sol_np1%u = 0.0_dp

  end subroutine define_problems_create_vectors


! delete old problem

  subroutine delete_old_problems

    call delete ( problem )
    call delete ( sysmatrix )
    call delete ( input_probdef )
    call delete ( mesh )
    call delete ( sol_np1, sol_n, sol_nm1, sol_i )
    call delete ( rhsd )
    call delete ( oldvectors )
    call delete ( velx, vely, pval )
    call delete ( cval, muval )
    call delete ( viscosity )
    call delete ( composition )

    deallocate ( elems_allowed )
    deallocate ( elems_refine )

  end subroutine delete_old_problems



! build and solve the velocity-pressure-c-mu problem

  subroutine build_solve_up_di

    print *,'Starting iterations'

!   fill the viscosity in the gauss points ( using a prediction for the
!   composition: chat = 2*c_n - c_nm1 )
    if ( step == 1 ) then
      composition%u = sol_n%u(cval%s)
    else
      composition%u = 2._dp*sol_n%u(cval%s) - sol_nm1%u(cval%s)
    end if
    call loop_over_elements ( mesh, problem, elemsub=fill_eta_gauss_c, &
      coefficients=coefficients, oldvectors=oldvectors )


!   initial guess is solution from previous step
    sol_np1%u(cval%s) = sol_n%u(cval%s)
    sol_np1%u(muval%s) = sol_n%u(muval%s)
    sol_np1%u(velx%s) = sol_n%u(velx%s)
    sol_np1%u(vely%s) = sol_n%u(vely%s)

!   fill solution vector with essential boundary conditions
    call fill_sysvector ( mesh, problem, sol_np1, curve1=1, physq=physqvel, &
      degfd=2, value=0._dp )
    call fill_sysvector ( mesh, problem, sol_np1, curve1=2, curve2=4, &
      physq=physqvel, vfunc=velfunc, vfuncnr=1 )
    call fill_sysvector ( mesh, problem, sol_np1, &
      point=1, physq=physqpress, value=0._dp )

    iter = 0

    newton: do

      iter = iter + 1


!     fill the viscosity in the gauss points (using sol_np1)
      if ( implicit_viscosity ) then
        composition%u = sol_np1%u(cval%s)
        call loop_over_elements ( mesh, problem, elemsub=fill_eta_gauss_c, &
          coefficients=coefficients, oldvectors=oldvectors )
      end if


!     Stokes flow : velocity / pressure
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=generalized_stokes_elem, oldvectors=oldvectors, &
         coefficients=coefficients, physqrow=[physqvel,physqpress], &
         physqcol=[physqvel,physqpress] )


      select case ( capillary_form )

      case(1) ! potential form 1

!       add the mugradc term to the matrix/vector
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=mugradc_implicit_elem, oldvectors=oldvectors, &
          physqrow=[physqvel], physqcol=[physqc,physqmu], addmatvec=.true., &
          coefficients=coefficients )

      case(2) ! potential form 2

!       add the -cgradmu term to the matrix/vector
        call build_system ( mesh, problem, sysmatrix, rhsd, &
          elemsub=mincgradmu_implicit_elem, oldvectors=oldvectors, &
          physqrow=[physqvel], physqcol=[physqc,physqmu], addmatvec=.true., &
          coefficients=coefficients )

      case(3) ! stress form

!       add the div tau_c to the matrix/vector
        call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
          elemsub=divstressdi_implicit_elem, oldvectors=oldvectors, &
          physqrow=[physqvel], physqcol=[physqc], coefficients=coefficients )

!       set additional block to zero if stress form is used
        call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
          buildvector=.false., physqrow=[physqvel], &
          physqcol=[physqmu], zeromatvec=.true. )

      case default

        write(*,'(/a,i0/)') 'Error: wrong value capillary_form: ', capillary_form
        stop

      end select

!     build matrix and vector for time integration and explicit convection term
      call build_system ( mesh, problem, sysmatrix, addmatvec=.true., &
        sysvector=rhsd, elemsub=diffuse_interface_elem1, &
        oldvectors=oldvectors, coefficients=coefficients, &
        physqrow=[physqc,physqmu], physqcol=[physqc,physqmu] )


!     build matrix and vector for implicit convection terms
      call build_system ( mesh, problem, sysmatrix, addmatvec=.true., &
        sysvector=rhsd, elemsub=di_convection_implicit_newton_elem, &
        oldvectors=oldvectors, coefficients=coefficients, &
        physqrow=[physqc], physqcol=[physqvel] )


!     set remaining blocks to zero
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        buildvector=.false., physqrow=[physqc], &
        physqcol=[physqpress], zeromatvec=.true. )
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        buildvector=.false., physqrow=[physqmu], &
        physqcol=[physqvel,physqpress], zeromatvec=.true. )
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        buildvector=.false., physqrow=[physqpress], &
        physqcol=[physqc,physqmu], zeromatvec=.true. )


      call check ( sysmatrix )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol_np1, rhsd )

      sol_i%u = sol_np1%u

!     solve (c,mu)
      solver_options_ma41%real_storage=rs_up_di
      solver_options_ma41%integer_storage=is_up_di
      solver_options_ma41%scaling=1

      call solve_system_ma41 ( sysmatrix, rhsd, sol_np1, &
        solver_options=solver_options_ma41 )

      cmax  = maxval(abs(sol_np1%u(cval%s))) ! c maximum
      cdiff = maxval(abs(sol_np1%u(cval%s)-sol_i%u(cval%s))) ! c diff

      mumax  = maxval(abs(sol_np1%u(muval%s))) ! mu maximum
      mudiff = maxval(abs(sol_np1%u(muval%s)-sol_i%u(muval%s))) ! mu diff

      velxmax  = maxval(abs(sol_np1%u(velx%s))) ! velx maximum
      velxdiff = maxval(abs(sol_np1%u(velx%s)-sol_i%u(velx%s))) ! velx diff

      velymax  = maxval(abs(sol_np1%u(vely%s))) ! vely maximum
      velydiff = maxval(abs(sol_np1%u(vely%s)-sol_i%u(vely%s))) ! vely diff

      print *,'iter = ',iter
!      print *,'        cdiff    = ',cdiff/cmax
!      print *,'        mudiff   = ',mudiff/mumax
!      print *,'        velxdiff = ',velxdiff/velxmax
!      print *,'        velydiff = ',velydiff/velymax
      print *,'        maxdiff  = ',maxval([cdiff/cmax, mudiff/mumax, &
                                            velxdiff/velxmax, velydiff/velymax])

      if ( cdiff < epsdi * cmax .and. velxdiff < epsdi * velxmax .and. &
        velydiff < epsdi * velymax .and. mudiff < epsdi * mumax ) exit newton


      if ( iter >= itermax ) then
        write(*,'(a,i0)') ' too many iterations: ', itermax
        stop
      end if

    end do newton

    print *,'Diffuse-interface problem solved!'

  end subroutine build_solve_up_di


! write the data to .vtk files

  subroutine postprocessing ( mesh, sol, ipost )

    type(mesh_t), intent(inout) :: mesh
    type(sysvector_t), intent(in), target :: sol
    integer, intent(in) :: ipost

    type(oldvectors_t) :: oldvectors_post
    type(vector_t) :: pressure

    if ( .not. mesh%meshparts) call fill_mesh_parts ( mesh )

    call create_oldvectors ( oldvectors_post, nsysvec=1, nsysvec2=1 )
    oldvectors_post%s(1)%p => sol

!   temporarily use the vector for c to derive the pressure
    call create_vector ( problem, pressure, vec=1 )

!   derive the pressure in all nodes
    call derive_vector ( mesh, problem, pressure, &
      elemsub=stokes_pressure, coefficients=coefficients, &
      oldvectors=oldvectors_post )

    write(filename,'(a,i4.4,a)') 'di_flow', ipost, '.vtk'
    call write_scalar_vtk ( mesh, problem, vector=pressure, &
      dataname='pressure',  filename=filename )

    call write_vector_vtk ( mesh, problem, filename=filename, &
      dataname='velocity', sysvector=sol, physq=physqvel, &
      append=.true. )

    call write_scalar_vtk ( mesh, problem, sysvector=sol, physq=physqc, &
      dataname='c',  filename=filename, append=.true. )

    call write_scalar_vtk ( mesh, problem, sysvector=sol, physq=physqmu, &
      dataname='mu',  filename=filename, append=.true. )

   call delete(pressure)
   call delete(oldvectors_post)

  end subroutine postprocessing


! P2 project the velocities and coordinates on the new mesh

  subroutine project_old_solutions_P2

    type(input_probdef_t) :: input_probdef_proj
    type(problem_t) :: problem_proj
    type(sysmatrix_t) :: sysmatrix_proj
    type(sysvector_t) :: rhsd_proj(ncompv_P2)
    type(sysvector_t) :: sol_proj(ncompv_P2)
    type(oldvectors_t) :: oldvectors_proj
    type(coefficients_t) :: coefficients_proj
    type(vector_t), target :: vec_p
    type(lu_ma57_t) :: lu2

    integer :: i

!   create vector for projection
    call create ( problem_temp, vec_p, vec=7 )

!   transfer the velocities at n
    call transfer_data ( mesh_temp, problem_temp, &
      sysvector1=sol_n_temp, vector2=vec_p, &
      physq1=[physqvel], degfd1=[1,2], degfd2=[1,2] )

!   transfer c at time step n
    call transfer_data ( mesh_temp, problem_temp, &
      sysvector1=sol_n_temp, vector2=vec_p, physq1=[physqc], degfd2=[3] )

!   transfer c at time step n-1
    call transfer_data ( mesh_temp, problem_temp, &
      sysvector1=sol_nm1_temp, vector2=vec_p, physq1=[physqc], degfd2=[4] )

!   transfer mu at time step n
    call transfer_data ( mesh_temp, problem_temp, &
      sysvector1=sol_n_temp, vector2=vec_p, physq1=[physqmu], degfd2=[5] )

!   problem definition for projection
    call create_input_probdef ( mesh, input_probdef_proj )

    input_probdef_proj%elementdof(1)%a = 1
    call problem_definition ( input_probdef_proj, mesh, problem_proj )

!   create system vectors (solution and right-hand side)
    call create ( problem_proj, sol_proj )
    call create ( problem_proj, rhsd_proj )

!   create system matrix
    call create_sysmatrix_structure ( sysmatrix_proj, mesh, &
      problem_proj, symmetric=.true. )
    call create_sysmatrix_data ( sysmatrix_proj )

!   fill coefficients
    call create_coefficients ( coefficients_proj, ncoefi=100, ncoefr=50 )
    coefficients_proj%i = 0
    coefficients_proj%i(1:4) = [ uintpl, ncompv_P2, uintpl, uintpl ]
    coefficients_proj%i(10) = gauss_proj
    coefficients_proj%i(23) = coorsys
    coefficients_proj%i(40) = inttype_proj
    coefficients_proj%r = 0

!   oldvectors
    call create_oldvectors ( oldvectors_proj, nvec=1, nprob=1, nmesh=1 )
    oldvectors_proj%v(1)%p => vec_p
    oldvectors_proj%p(1)%p => problem_temp
    oldvectors_proj%m(1)%p => mesh_temp

!   build (assemble) matrix and vector from elements
    call build_system ( mesh, problem_proj, sysmatrix_proj, &
      msysvector=rhsd_proj, elemsub=projection_elem, &
      coefficients=coefficients_proj, oldvectors=oldvectors_proj )

!   solve the projection problem
    do i = 1, ncompv_P2
      call add_effect_of_essential_to_rhs ( problem_proj, sysmatrix_proj, &
       sol_proj(i), rhsd_proj(i) )
      call solve_system_ma57 ( sysmatrix_proj, rhsd_proj(i), sol_proj(i), &
        lu=lu2 )
    end do

!   update the old velocites
!   note: in the projection, essential BCs could be used on the walls, since
!   these are known. However, to do all projections with the same matrix, this
!   was not done, and the sol_proj vectors retain nodal numbering.
!   To transfer the values, the 'target' vector should also be given in
!   nodal numbering. The easiest way to do this is by using subscripts (which
!   do not only select the correct dof's, but also give them in nodal
!   numbering).
    sol_n%u(velx%s) = sol_proj(1)%u
    sol_n%u(vely%s) = sol_proj(2)%u

!   update old composition
    sol_n%u(cval%s) = sol_proj(3)%u
    sol_nm1%u(cval%s) = sol_proj(4)%u
    sol_n%u(muval%s) = sol_proj(5)%u

!   delete definitions used for the projection problem
    call delete ( lu2 )
    call delete( input_probdef_proj )
    call delete( problem_proj )
    call delete( sysmatrix_proj )
    call delete( rhsd_proj )
    call delete( sol_proj )
    call delete( oldvectors_proj )
    call delete( coefficients_proj )
    call delete( vec_p )

  end subroutine project_old_solutions_P2


! subroutine to give coordinates of quarter of a circle

  subroutine coors_circle ( coor )

    real(dp), dimension(:,:), intent(inout) :: coor

    integer :: np, i
    real(dp) :: p(size(coor,1))

    np = size(coor,1)
    p = [(2*pi*(i-1)/np,i=1,np)]

    coor(:,1) = drop_r * cos(p) + drop_z0
    coor(:,2) = drop_r * sin(p) + 0.0_dp

  end subroutine coors_circle


! find the initual velocity field

  subroutine find_velocity_init

    type(input_probdef_t) :: input_probdef_init
    type(problem_t), target :: problem_init
    type(sysvector_t), target :: sol_init, rhs_init
    type(subscript_t) :: cval_init, muval_init
    type(sysmatrix_t) :: sysmatrix_init
    type(oldvectors_t) :: oldvectors_init
    type(vector_t), target :: composition_init
    type(elvector_t), target :: viscosity_init
    type(subscript_t) :: velx_init, vely_init, pval_init

!   problem definition for the diffuse-interface problem
    call create_input_probdef ( mesh, input_probdef_init, nvec=4, nphysq=4 )

    input_probdef_init%vec_elementdof(1)%a =   &
      reshape ( [ 1,1,1,1,1,1,       &  ! c
                   1,1,1,1,1,1,       &  ! mu
                   2,2,2,2,2,2,       &  ! velocity
                   1,0,1,0,1,0 ],    &  ! pressure
                   [6,4] )

    input_probdef_init%physq = [1,2,3,4]
    input_probdef_init%probnr = 3

!   Dirichlet boundary conditions

    call define_essential ( mesh, input_probdef_init, curve1=1, degfd=[0,1], &
      physq=physqvel )
    call define_essential ( mesh, input_probdef_init, curve1=2, curve2=4, &
      physq=physqvel )
    call define_essential ( mesh, input_probdef_init, point=1, &
      physq=physqpress )

!   c is known at t=0
    call define_essential ( mesh, input_probdef_init, elgroup1=1, &
      physq=physqc )

    call problem_definition ( input_probdef_init, mesh, problem_init )

!   create system vectors (solution and right-hand side) for di problem
    call create ( problem_init, sol_init, rhs_init )
    call create_elvector ( mesh, viscosity_init, nreal1d=1 )

!   create vector for the composition (used to determine the local viscosity)
    call create_vector ( problem_init, composition_init, vec=1 )

!   create subscripts for c and mu
    call create_subscript ( mesh, problem_init, cval_init, physqarr=[physqc] )
    call create_subscript ( mesh, problem_init, muval_init, physqarr=[physqmu] )
    call create_subscript ( mesh, problem_init, velx_init, &
      physqarr=[physqvel], degfd=1 )
    call create_subscript ( mesh, problem_init, vely_init, &
      physqarr=[physqvel], degfd=2 )
    call create_subscript ( mesh, problem_init, pval_init, &
      physqarr=[physqpress] )

!   initialize vecotr to zero
    sol_init%u = 0

!   fill the initial concentration field
    call fill_sysvector ( mesh, problem_init, sol_init, physq=physqc, &
      node1=1, node2=mesh%nnodes, func=difunc, funcnr=1)

!   create oldvectors
    call create_oldvectors ( oldvectors_init, nsysvec=5, nprob=2, nvec=1, &
      nelvec=1 )
      oldvectors_init%s(1)%p => sol_init
      oldvectors_init%s(2)%p => sol_init
      oldvectors_init%s(3)%p => sol_init
      oldvectors_init%p(1)%p => problem_init
      oldvectors_init%p(2)%p => problem_init
      oldvectors_init%e(1)%p => viscosity_init
      oldvectors_init%v(1)%p => composition_init


!   create system matrix for the di problem
    call create_sysmatrix_structure_base ( sysmatrix_init, mesh, problem_init, &
      symmetric=.false. )
    call finalize_sysmatrix_structure ( sysmatrix_init )
    call create_sysmatrix_data ( sysmatrix_init )


!   fill the viscosity in the gauss points
    composition_init%u = sol_init%u(cval_init%s)
    call loop_over_elements ( mesh, problem_init, elemsub=fill_eta_gauss_c, &
      coefficients=coefficients, oldvectors=oldvectors_init )

!   fill solution vector with essential boundary conditions
    call fill_sysvector ( mesh, problem_init, sol_init, curve1=1, &
      physq=physqvel, degfd=2, value=0._dp )
    call fill_sysvector ( mesh, problem_init, sol_init, curve1=2, &
      curve2=4, physq=physqvel, vfunc=velfunc, vfuncnr=1 )
    call fill_sysvector ( mesh, problem_init, sol_init, &
      point=1, physq=physqpress, value=0._dp )


!   Stokes flow : velocity / pressure
    call build_system ( mesh, problem_init, sysmatrix_init, rhs_init, &
      elemsub=generalized_stokes_elem, oldvectors=oldvectors_init, &
       coefficients=coefficients, physqrow=[physqvel,physqpress], &
       physqcol=[physqvel,physqpress] )


    select case ( capillary_form )

    case(1) ! potential form 1

!     add the mugradc term to the matrix/vector
      call build_system ( mesh, problem_init, sysmatrix_init, rhs_init, &
        elemsub=mugradc_implicit_elem, oldvectors=oldvectors_init, &
        physqrow=[physqvel], physqcol=[physqc,physqmu], addmatvec=.true., &
        coefficients=coefficients )

    case(2) ! potential form 2

!     add the -cgradmu term to the matrix/vector
      call build_system ( mesh, problem_init, sysmatrix_init, rhs_init, &
        elemsub=mincgradmu_implicit_elem, oldvectors=oldvectors_init, &
        physqrow=[physqvel], physqcol=[physqc,physqmu], addmatvec=.true., &
        coefficients=coefficients )

    case(3) ! stress form

!     add the div tau_c to the matrix/vector
      call build_system ( mesh, problem_init, sysmatrix_init, rhs_init, &
        addmatvec=.true., &
        elemsub=divstressdi_implicit_elem, oldvectors=oldvectors_init, &
        physqrow=[physqvel], physqcol=[physqc], coefficients=coefficients )

!     set additional block to zero if stress form is used
      call build_system ( mesh, problem_init, sysmatrix_init, rhs_init, &
        addmatvec=.true., &
        buildvector=.false., physqrow=[physqvel], &
        physqcol=[physqmu], zeromatvec=.true. )

    case default

      write(*,'(/a,i0/)') 'Error: wrong value capillary_form: ', capillary_form
      stop

    end select

!   build matrix and vector for time integration and explicit convection term
    call build_system ( mesh, problem_init, sysmatrix_init, addmatvec=.true., &
      sysvector=rhs_init, elemsub=diffuse_interface_elem1, &
      oldvectors=oldvectors_init, coefficients=coefficients, &
      physqrow=[physqc,physqmu], physqcol=[physqc,physqmu] )


!   build matrix and vector for implicit convection terms
    call build_system ( mesh, problem_init, sysmatrix_init, addmatvec=.true., &
      sysvector=rhs_init, elemsub=di_convection_implicit_newton_elem, &
      oldvectors=oldvectors_init, coefficients=coefficients, &
      physqrow=[physqc], physqcol=[physqvel] )


!   set remaining blocks to zero
    call build_system ( mesh, problem_init, sysmatrix_init, rhs_init, &
      addmatvec=.true., buildvector=.false., physqrow=[physqc], &
      physqcol=[physqpress], zeromatvec=.true. )
    call build_system ( mesh, problem_init, sysmatrix_init, rhs_init, &
      addmatvec=.true., buildvector=.false., physqrow=[physqmu], &
      physqcol=[physqvel,physqpress], zeromatvec=.true. )
    call build_system ( mesh, problem_init, sysmatrix_init, rhs_init, &
      addmatvec=.true., buildvector=.false., physqrow=[physqpress], &
      physqcol=[physqc,physqmu], zeromatvec=.true. )

    call check ( sysmatrix_init )

    call add_effect_of_essential_to_rhs ( problem_init, sysmatrix_init, &
      sol_init, rhs_init )

!   solve (c,mu)
    solver_options_ma41%real_storage=rs_up_di
    solver_options_ma41%integer_storage=is_up_di
    solver_options_ma41%scaling=1

    call solve_system_ma41 ( sysmatrix_init, rhs_init, sol_init, &
      solver_options=solver_options_ma41 )

    sol_np1%u(cval%s) = sol_init%u(cval_init%s)
    sol_np1%u(muval%s)= sol_init%u(muval_init%s)
    sol_np1%u(velx%s) = sol_init%u(velx_init%s)
    sol_np1%u(vely%s) = sol_init%u(vely_init%s)
    sol_np1%u(pval%s) = sol_init%u(pval_init%s)

  end subroutine find_velocity_init


! velocity function to impose a hyperbolic velocity field (the velocity field
! is chosen such that inflow=outflow)

  function velfunc ( n, nr, x )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: velfunc

    select case(nr)

      case(1)

         if ( coorsys == 1 ) then
           velfunc = elong_rate * 2._dp * [ -(x(1)-lz/2._dp), x(2)/2._dp ] &
                       / lr
         else
           velfunc = elong_rate * [ -(x(1)-lz/2._dp), x(2) ] / lr
         end if

      case default
        write(*,'(/a,i0/)') 'Error velfunc: wrong function number: ', nr
        stop

    end select

  end function velfunc



end program diffuse_interface18
