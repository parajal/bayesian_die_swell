! axisymmetric diffuse-interface problem of a droplet on a surface
! Stokes problem with mu grad c term for the surface tension
! constant viscosity
! variable (static) contact angle
! second order time integration
! Newton-Raphson iterations for the non-linear term in the
! Cahn-Hilliard equation
! This problem is similar to diffuse_interface13, but now with
! adaptive meshing using Gmsh (triangular elements)
! and second-order time-integration

program diffuse_interface14

  use tfem_m
  use math_defs_m
  use hsl_ma41_m
  use hsl_ma57_m
  use generalized_stokes_elements_m
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
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 6,          & ! 6-point Gauss integration of triangles
    gaussb = 3,         & ! 3-point integration of boundary elements
    gauss_proj = 6,     & ! integration rule for the projection problem
    inttype_proj = 3,   & ! standard Gauss-Legendre (numerical table)
    coorsys =1,         & ! coordinate system
    vtkevery = 10,      & ! vtk file every vtkevery steps. -1: means none
    ncompv_P2 = 6,      & ! number of components for P2 projection: c,
                          ! and x,y-velocities at n and n-1
    numtimesteps = 100,&  ! number of time steps
    itermax = 50,       & ! maximum number of Newton iterations
    n_init = 1000         ! initial number of refinement points on DI


! variables

  integer :: &
    timeint1 = 1,       & ! first-order, semi-implicit Euler
    timeint2 = 2          ! second-order Gear with velocity prediction

  real(dp) :: &
    eta = 0.1_dp,         & ! viscosity of fluid 1 (c=1)
    rho = 100.0_dp,       & ! density in diffuse-interface method
    alpha = 1.0_dp,       & ! parameter in the diffuse-interface method
    beta = 1.0_dp,        & ! parameter in the diffuse-interface method
    Mcoef = 1.e-2_dp,     & ! parameter in the diffuse-interface method
    kappa = 1.e-4_dp,     & ! parameter in the diffuse-interface method
    theta_c = 70_dp,      & ! contact angle measured through drop (degrees)
    epsdi = 1.e-9_dp,     & ! accuracy in Newton iteration
    drop_z0 = 1.0_dp,     & ! initial z-coordinate of center of drop
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
    deltat = 0.001_dp,    & ! time step
    dist_threshold = 0.34_dp,&! elements for which any (abs(c)<dist_threshold)
                              ! are allowed to contain the interface
    rs_up = 1.4_dp,       &! real_storage velocity-pressure LU (HSL)
    is_up = 1.5_dp,       & ! integer_storage velocity-pressure LU (HSL)
    rs_di  = 1.4_dp,      & ! real_storage for the diffuse interface (HSL)
    is_di  = 1.6_dp         ! integer_storage for diffuse interface LU (HSL)


! definitions for DI and up problems

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefdi
  type(problem_t), target :: problem, problemdi
  type(sysmatrix_t) :: sysmatrix, sysmatrixdi
  type(sysvector_t), target :: sol_np1, sol_n, sol_nm1
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_di
  type(coefficients_t) :: coefficients
  type(sysvector_t), target :: soldi_np1, soldi_n, soldi_i
  type(sysvector_t), target :: soldi_nm1, rhsdi
  type(solver_options_ma41_t) :: solver_options_ma41
  type(solver_options_ma57_t) :: solver_options_ma57
  type(subscript_t) :: velx, vely
  type(subscript_t) :: cval, muval
  type(lu_ma57_t) :: lu_u

! temp definitions used to store old definitions before projection

  type(mesh_t), target :: mesh_temp
  type(sysvector_t) :: soldi_n_temp, soldi_nm1_temp
  type(sysvector_t) :: sol_n_temp, sol_nm1_temp
  type(problem_t), target :: problem_temp, problemdi_temp

! type definition for the refinement fields

  type(refinement_fields_t) :: refinement_fields

! variables

  integer :: step=0, ipost=0, iter

  real(dp) :: cmax, cdiff, phi_di, c_B, surf_ten, x_radius(2), x_height(2)
  real(dp), allocatable, dimension(:,:) :: refinement_coor
  logical, allocatable, dimension(:) :: elems_allowed, elems_refine

  character(len=30) :: filename
  logical :: remeshing=.false.


! pass value to functions_m

  xi = sqrt(kappa/alpha)
  allocate ( c(2) )
  c(1) = drop_z0
  c(2) = 0.0_dp
  r = drop_r


! determine the wetting potential phi

  c_B = sqrt(alpha/beta)
  surf_ten = rho*(2*sqrt(2.0_dp)/3) * (kappa * c_B**2) / xi

  phi_di = (surf_ten/rho) * cos(theta_c*2*pi/360._dp)/(2*(c_B-(c_B**3)/3.0_dp))


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
  coefficients%r(1) = eta


! fill coefficients for the diffuse interface problem

  coefficients%r(101:106)  = [ deltat,  Mcoef,  alpha,  beta,  kappa, rho ]
  coefficients%i(151) = timeint1 ! Euler integration of the di equations
  coefficients%i(152) = 2        ! Newton-Raphson iteration for di equation
  coefficients%r(108) = phi_di   ! wetting potential


! generate and read the refined mesh using gmsh

  call generate_read_mesh_refined


! define the up and di problem, create (sys)vectors, subscripts, oldvectors
! and initialize sysvectors

  call define_problems_create_vectors


! sample the initial c-field and then find the initial mu-field

  call find_mu_init


! copy old value (already here, so soldi_np1 is used in the up problem!)

  call copy ( soldi_np1, soldi_n )


! find the elements that are allowed to contain the interface

  call find_elems_allowed


! build and solve the up system to find initial velocity

  call build_solve_up


! write VTK files

  call postprocessing ( mesh, sol_np1, soldi_np1, ipost=ipost )
  ipost = ipost + 1

! find the coordinates of c=0 on the line (z,r) = (1,0) to (1,1)
  call find_c0 ( mesh, problemdi, soldi_np1, x1=[1.0_dp, 0.0_dp], &
    x2=[1.0_dp,1.0_dp], x_c0=x_radius, tol=1.e-12_dp )
! find the coordinates of c=0 on the line (z,r) = (0,0) to (1,0)
  call find_c0 ( mesh, problemdi, soldi_np1, x1=[0.0_dp, 0.0_dp], &
    x2=[1.0_dp,0.0_dp], x_c0=x_height, tol=1.e-12_dp )

  open (unit = 100, file='output_drop.out')
  write(100, '(i5,3es16.8)') step, time, x_radius(2),1.0_dp-x_height(1)
  close(unit=100)

  print *,'Initial drop radius is ',x_radius(2)
  print *,'Initlal drop height is ',1.0_dp-x_height(1)


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
      call copy ( soldi_n, soldi_n_temp )
      call copy ( sol_n, sol_n_temp )
      call copy ( soldi_nm1, soldi_nm1_temp )
      call copy ( sol_nm1, sol_nm1_temp )
      call problem_definition ( input_probdef, mesh, problem_temp )
      call problem_definition ( input_probdefdi, mesh, problemdi_temp )


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
      call delete ( soldi_n_temp )
      call delete ( sol_n_temp )
      call delete ( soldi_nm1_temp )
      call delete ( sol_nm1_temp )
      call delete ( problemdi_temp )
      call delete ( problem_temp )


    end if


!   build and solve the diffuse interface problem

    call build_solve_di


!   find the elements that are allowed to contain the interface

    if ( remeshing ) call find_elems_allowed


!   copy old value (already here, so soldi_np1 is used in the up problem!)

    call copy ( soldi_n, soldi_nm1 )
    call copy ( soldi_np1, soldi_n )


!   build and solve the up system

    call build_solve_up


!   write VTK files

    if ( vtkevery > 0 ) then
      if ( mod(step,vtkevery) == 0 ) then
        call postprocessing ( mesh, sol_np1, soldi_np1, ipost=ipost )
        ipost = ipost + 1
      end if
    end if

!   find the coordinates of c=0 on the line (z,r) = (1,0) to (1,1)
    call find_c0 ( mesh, problemdi, soldi_np1, x1=[1.0_dp, 0.0_dp], &
      x2=[1.0_dp,1.0_dp], x_c0=x_radius, tol=1.e-12_dp )
!   find the coordinates of c=0 on the line (z,r) = (0,0) to (1,0)
    call find_c0 ( mesh, problemdi, soldi_np1, x1=[0.0_dp, 0.0_dp], &
      x2=[1.0_dp,0.0_dp], x_c0=x_height, tol=1.e-12_dp )

    open (unit = 100, file='output_drop.out', position='APPEND')
    write(100, '(i5,3es16.8)') step, time, x_radius(2),1.0_dp-x_height(1)
    close(unit=100)

    print *,'drop radius is ',x_radius(2)
    print *,'drop height is ',1.0_dp-x_height(1)


!   copy old values
!   note: soldi is already copied in the time integration loop)

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

    call create_vector ( problemdi, c, physq=1 )
    call extract_physvector (mesh, problemdi, soldi_np1, c )

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

    call create_vector ( problemdi, c, physq=1 )
    call extract_physvector (mesh, problemdi, soldi_np1, c )

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

    call create_vector ( problemdi, c, physq=1 )
    call extract_physvector (mesh, problemdi, soldi_np1, c )

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
    call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

    input_probdef%vec_elementdof(1)%a =   &
        reshape ( [2,2,2,2,2,2,     &  ! velocity
                    1,0,1,0,1,0,    &  ! pressure
                    1,1,1,1,1,1,    &  ! scalar
                    3,3,3,3,3,3],   &  ! symmetric tensor
                         [6,4] )

    input_probdef%physq = [1,2]
    input_probdef%probnr = 1

!   Dirichlet boundary conditions

    call define_essential ( mesh, input_probdef, curve1=1, degfd=[0,1], &
      physq=physqvel )
    call define_essential ( mesh, input_probdef, curve1=2, curve2=4, &
      physq=physqvel )
    call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

    call problem_definition ( input_probdef, mesh, problem )

!   problem definition for the diffuse-interface problem
    call create_input_probdef ( mesh, input_probdefdi, nvec=4, nphysq=2 )

    input_probdefdi%vec_elementdof(1)%a =   &
        reshape ( [1,1,1,1,1,1,   &  ! c
                   1,1,1,1,1,1,   &  ! mu
                   2,2,2,2,2,2,   &  ! coordinates
                   ncompv_P2,ncompv_P2,ncompv_P2,   &
                   ncompv_P2,ncompv_P2,ncompv_P2 ], &
                    [6,4] )

    input_probdefdi%physq = [1,2]
    input_probdefdi%probnr = 2

    call problem_definition ( input_probdefdi, mesh, problemdi )

!   create system vectors for the up problem
    call create_sysvector ( problem, sol_np1, sol_n, sol_nm1, rhsd )

!   create subscripts for the mesh velocity
    call create_subscript ( mesh, problem, velx, physqarr=[1], degfd=1 )
    call create_subscript ( mesh, problem, vely, physqarr=[1], degfd=2 )

!   create system vectors (solution and right-hand side) for di problem
    call create ( problemdi, soldi_np1, soldi_n, soldi_nm1, soldi_i, rhsdi )

!   create subscripts for c and mu
    call create_subscript ( mesh, problemdi, cval, physqarr=[1] )
    call create_subscript ( mesh, problemdi, muval, physqarr=[2] )

!   create system matrix for the up problem
    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric=.true. )
    call finalize_sysmatrix_structure ( sysmatrix )
    call create_sysmatrix_data ( sysmatrix )

!   create system matrix for the di problem
    call create_sysmatrix_structure_base ( sysmatrixdi, mesh, problemdi )
    call finalize_sysmatrix_structure ( sysmatrixdi )
    call create_sysmatrix_data ( sysmatrixdi )

!   create oldvectors
    call create_oldvectors ( oldvectors_di, nsysvec=5, nprob=2 )
      oldvectors_di%s(1)%p => sol_n
      oldvectors_di%s(2)%p => soldi_n
      oldvectors_di%s(3)%p => soldi_np1
      oldvectors_di%s(4)%p => soldi_nm1
      oldvectors_di%s(5)%p => sol_nm1
      oldvectors_di%p(1)%p => problem
      oldvectors_di%p(2)%p => problemdi

!   initialize all vectors to zero
    sol_nm1%u = 0.0_dp
    soldi_nm1%u = 0.0_dp
    sol_n%u = 0.0_dp
    soldi_n%u = 0.0_dp
    sol_np1%u = 0.0_dp
    soldi_np1%u = 0.0_dp
    soldi_i%u = 0.0_dp

  end subroutine define_problems_create_vectors


! delete old problem

  subroutine delete_old_problems

    call delete ( problem, problemdi )
    call delete ( sysmatrix, sysmatrixdi )
    call delete ( input_probdef, input_probdefdi )
    call delete ( mesh )
    call delete ( sol_np1, sol_n, sol_nm1 )
    call delete ( soldi_np1, soldi_n, soldi_nm1, soldi_i )
    call delete ( rhsdi )
    call delete ( rhsd )
    call delete ( oldvectors_di )
    call delete ( velx, vely )
    call delete ( cval, muval )
    call delete ( lu_u )

    deallocate ( elems_allowed )
    deallocate ( elems_refine )

  end subroutine delete_old_problems


! build and solve the up system

  subroutine build_solve_up

!   fill solution vector with essential boundary conditions
    sol_np1%u = 0.0_dp

!   set pressure level = 0 in lower left corner
    call fill_sysvector ( mesh, problem, sol_np1, &
      point=1, physq=physqpress, value=0._dp )

!   Stokes flow : velocity / pressure
    if ( remeshing .or. step == 0 ) then
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=generalized_stokes_elem, oldvectors=oldvectors_di, &
        coefficients=coefficients )
    end if

!   add the mugradc term to the rhs
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=rhs_mugradc, oldvectors=oldvectors_di, physqrow=[physqvel], &
      physqcol=[physqvel], buildmatrix = .false., &
      coefficients=coefficients )

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol_np1, rhsd )

!   solve up system
    solver_options_ma57%real_storage=rs_up
    solver_options_ma57%integer_storage=is_up

    call solve_system_ma57 ( sysmatrix, rhsd, sol_np1, lu=lu_u, &
      solver_options=solver_options_ma57 )

  end subroutine build_solve_up


! build and solve the diffuse interface problem

  subroutine build_solve_di

    print *,'Starting iterations'

    soldi_np1%u = soldi_n%u ! initial guess is solution from previous step

    iter = 0

    newton: do

      iter = iter + 1

!     build (assemble) matrix and vector for diffuse-interface problem

      call build_system ( mesh, problemdi, sysmatrixdi, &
        sysvector=rhsdi, elemsub=diffuse_interface_elem1, &
        oldvectors=oldvectors_di, coefficients=coefficients )

!     add the boundary conditon for the contact angle

      call add_boundary_elements ( mesh, problemdi, rhsdi, curve=2, &
        sysmatrix=sysmatrixdi, elemsub=di_natboun_curve, &
        coefficients=coefficients, oldvectors=oldvectors_di )

      call check ( sysmatrixdi )

      call add_effect_of_essential_to_rhs ( problemdi, sysmatrixdi, soldi_np1, rhsdi )

      soldi_i%u = soldi_np1%u

!     solve (c,mu)
      solver_options_ma41%real_storage=rs_di
      solver_options_ma41%integer_storage=is_di

      call solve_system_ma41 ( sysmatrixdi, rhsdi, soldi_np1, &
        solver_options=solver_options_ma41 )

      cmax  = maxval(abs(soldi_np1%u(cval%s))) ! c maximum
      cdiff = maxval(abs(soldi_np1%u(cval%s)-soldi_i%u(cval%s))) ! c diff

      print *,'iter = ',iter,'    cdiff = ',cdiff

      if ( cdiff < epsdi * cmax ) exit newton

      if ( iter >= itermax ) then
        write(*,'(a,i0)') ' too many iterations: ', itermax
        stop
      end if

    end do newton

    print *,'Diffuse-interface problem solved!'

  end subroutine build_solve_di


! find the initual mu field

  subroutine find_mu_init

    type(input_probdef_t) :: input_probdefdi_init
    type(problem_t), target :: problemdi_init
    type(sysvector_t), target :: soldi_init, rhsdi_init
    type(subscript_t) :: cval_init, muval_init
    type(sysmatrix_t) :: sysmatrixdi_init
    type(oldvectors_t) :: oldvectors_di_init

!   problem definition for the diffuse-interface problem
    call create_input_probdef ( mesh, input_probdefdi_init, nvec=2, nphysq=2 )

    input_probdefdi_init%vec_elementdof(1)%a =   &
      reshape ( [ 1,1,1,1,1,1,       &  ! c
                   1,1,1,1,1,1 ],    &  ! mu
                   [6,2] )

    input_probdefdi_init%physq = [1,2]
    input_probdefdi_init%probnr = 3

!   c is known at t=0
    call define_essential ( mesh, input_probdefdi_init, elgroup1=1, &
      physq=1 )
    call problem_definition ( input_probdefdi_init, mesh, problemdi_init )

!   create system vectors (solution and right-hand side) for di problem
    call create ( problemdi_init, soldi_init, rhsdi_init )

    call fill_sysvector ( mesh, problemdi_init, soldi_init, physq=1, &
      node1=1, node2=mesh%nnodes, func=difunc, funcnr=1)

!   create subscripts for c and mu
    call create_subscript ( mesh, problemdi_init, cval_init, physqarr=[1] )
    call create_subscript ( mesh, problemdi_init, muval_init, physqarr=[2] )

!   set velocity to zero to find the inital mu-field
    sol_np1%u = 0

!   create oldvectors
    call create_oldvectors ( oldvectors_di_init, nsysvec=5, nprob=2, nvec=1 )
      oldvectors_di_init%s(1)%p => sol_np1
      oldvectors_di_init%s(2)%p => soldi_init
      oldvectors_di_init%s(3)%p => soldi_init
      oldvectors_di_init%p(1)%p => problem
      oldvectors_di_init%p(2)%p => problemdi_init

!   create system matrix for the di problem
    call create_sysmatrix_structure_base ( sysmatrixdi_init, mesh, problemdi_init )
    call finalize_sysmatrix_structure ( sysmatrixdi_init )
    call create_sysmatrix_data ( sysmatrixdi_init )

!   build (assemble) matrix and vector for diffuse-interface problem
    call build_system ( mesh, problemdi_init, sysmatrixdi_init, &
      sysvector=rhsdi_init, elemsub=diffuse_interface_elem1, &
        oldvectors=oldvectors_di_init, coefficients=coefficients )

!   add the boundary conditon for the contact angle
    call add_boundary_elements ( mesh, problemdi_init, sysvector=rhsdi_init, &
      curve=2, sysmatrix=sysmatrixdi_init, elemsub=di_natboun_curve, &
      coefficients=coefficients, oldvectors=oldvectors_di_init )

    call check ( sysmatrixdi_init )

    call add_effect_of_essential_to_rhs ( problemdi_init, sysmatrixdi_init, &
      soldi_init, rhsdi_init )

!   solve (c,mu)
    solver_options_ma41%real_storage=rs_di
    solver_options_ma41%integer_storage=is_di

    call solve_system_ma41 ( sysmatrixdi_init, rhsdi_init, soldi_init, &
      solver_options=solver_options_ma41 )

    soldi_np1%u(cval%s) = soldi_init%u(cval_init%s)
    soldi_np1%u(muval%s) = soldi_init%u(muval_init%s)

    call delete ( input_probdefdi_init )
    call delete ( problemdi_init )
    call delete ( soldi_init, rhsdi_init )
    call delete ( cval_init, muval_init )

  end subroutine find_mu_init


! write the data to .vtk files

  subroutine postprocessing ( mesh, sol, soldi, ipost )

    type(mesh_t), intent(inout) :: mesh
    type(sysvector_t), intent(in), target :: sol
    type(sysvector_t), intent(in), target :: soldi
    integer, intent(in) :: ipost

    type(oldvectors_t) :: oldvectors_post
    type(vector_t) :: pressure

    if ( .not. mesh%meshparts) call fill_mesh_parts ( mesh )

    call create_oldvectors ( oldvectors_post, nsysvec=1, nsysvec2=1 )
    oldvectors_post%s(1)%p => sol

    call create_vector ( problem, pressure, vec=3 )

!   derive the pressure in all nodes
    call derive_vector ( mesh, problem, pressure, &
      elemsub=stokes_pressure, coefficients=coefficients, &
      oldvectors=oldvectors_post )

    write(filename,'(a,i4.4,a)') 'flow', ipost, '.vtk'
    call write_scalar_vtk ( mesh, problem, vector=pressure, &
      dataname='pressure',  filename=filename )

    call write_vector_vtk ( mesh, problem, filename=filename, &
      dataname='velocity', sysvector=sol, physq=physqvel, &
      append=.true. )

    write(filename,'(a,i4.4,a)') 'di', ipost, '.vtk'
    call write_scalar_vtk ( mesh, problemdi, sysvector=soldi, physq=1, &
      dataname='c',  filename=filename )

    call write_scalar_vtk ( mesh, problemdi, sysvector=soldi, physq=2, &
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
    call create ( problemdi_temp, vec_p, vec=4 )

!   transfer the velocities at n
    call transfer_data ( mesh_temp, problem1=problem_temp, &
      problem2=problemdi_temp,sysvector1=sol_n_temp, vector2=vec_p, &
      physq1=[1], degfd1=[1,2], degfd2=[1,2] )

!   transfer the velocities at n-1
    call transfer_data ( mesh_temp, problem1=problem_temp, &
      problem2=problemdi_temp,sysvector1=sol_nm1_temp, vector2=vec_p, &
      physq1=[1], degfd1=[1,2], degfd2=[3,4] )

!   transfer c at time step n
    call transfer_data ( mesh_temp, problemdi_temp, &
      sysvector1=soldi_n_temp, vector2=vec_p, physq1=[1], degfd2=[5] )

!   transfer c at time step n-1
    call transfer_data ( mesh_temp, problemdi_temp, &
      sysvector1=soldi_nm1_temp, vector2=vec_p, physq1=[1], degfd2=[6] )

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
    oldvectors_proj%p(1)%p => problemdi_temp
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
    sol_nm1%u(velx%s) = sol_proj(3)%u
    sol_nm1%u(vely%s) = sol_proj(4)%u

!   update old concentration
    soldi_n%u(cval%s) = sol_proj(5)%u
    soldi_nm1%u(cval%s) = sol_proj(6)%u

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


end program diffuse_interface14
