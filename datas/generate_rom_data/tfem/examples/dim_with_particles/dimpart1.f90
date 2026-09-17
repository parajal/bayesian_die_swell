!  Particle-covered droplet sheared between two plates.
!  Diffuse-interface model for the fluid-fluid interface.
!  Sharp-interface model for the particle-fluid interface.
!  Boundary-fitted meshes for the particles using an ALE formulation.
!  Surface tension added using the stress form of the Cahn-Hilliard model.
!  Newtonian fluids with different viscosities.
!  Variable (static) contact angle fluid-fluid interface with particle boundary.
!  Second order time integration.
!  Adaptive meshing using Gmsh (triangular elements).
!  Newton-Raphson iterations for the non-linear term in the
!  Cahn-Hilliard equation.
!  Flow and diffuse-interface equations are solved fully-implicit (in one
!  system).

program dimpart1

  use tfem_m
  use math_defs_m
  use hsl_ma41_m
  use hsl_ma57_m
  use generalized_stokes_elements_m
  use diffuse_interface_elements_m
  use io_utils_m
  use subs_m
  use update_mesh_nodes_multi_m
  use projection_elements_m
  use functions_m

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
    vtkevery = 1,       & ! vtk file every vtkevery steps. -1: means none
    ncompv_P2 = 9,      & ! number of components for P2 projection: c at n
                          ! and n-1, and velocities and mu at n,  coords n+nm1
    numtimesteps = 200, & ! number of time steps
    ndrops = 1,         & ! number of drops
    nparticles = 4,     & ! number of particles
    itermax = 50,       & ! maximum number of Newton iterations
    n_init = 100          ! initial number of refinement points on DI

! variables

  integer :: &
    timeint1 = 1,         & ! first-order Euler
    timeint2 = 2            ! second-order Gear

  real(dp) :: &
    eta1 = 1._dp,         & ! viscosity of the ambient fluid (c=1)
    eta2 = 0.2_dp,        & ! viscosity of the droplet (c=-1)
    shear_rate = 1._dp,   & ! imposed shear rate
    radiusdrop = 1._dp,   & ! droplet radius
    radiuspart = 0.2_dp,  & ! particle radius
    rho = 100._dp,        & ! density in diffuse-interface method
    alpha = 1.0_dp,       & ! parameter in the diffuse-interface method
    beta = 1.0_dp,        & ! parameter in the diffuse-interface method
    Mcoef = 0.01_dp,      & ! parameter in the diffuse-interface method
    kappa = 1.e-2_dp ,    & ! parameter in the diffuse-interface method
    theta_c = 120._dp,    & ! contact angle measured through droplet
    epsdi = 1.e-6_dp,     & ! accuracy in Newton iteration
    min_elem = 2._dp,     & ! minimum number of elements between particles
    nelem_part = 16._dp,  & ! number of elements on the particles
    nelem_int = 1._dp,    & ! number of elements across the interface (=xi)
    lx = 5.0_dp,          & ! length of the domain in x direction
    ly = 5.0_dp,          & ! length of the domain in y direction
    deltat = 0.04_dp,     & ! time step
    dist_threshold = 0.34_dp,&! elements for which any (abs(c)<dist_threshold)
                              ! are allowed to contain the interface
    rs_up_di = 2.4_dp,    & ! real_storage (HSL)
    is_up_di = 2.6_dp,    & ! integer_storage LU (HSL)
    area_threshold = 1.39_dp,        & ! remeshing area threshold
    aspect_ratio_threshold = 1.39_dp   ! remeshing aspect ratio threshold

! definitions for DI and up problems

  type(mesh_t) :: mesh, meshup, meshdown
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
  type(vector_t), target :: meshvel
  type(problem_t) :: problem_lapl

! temp definitions used to store old definitions before projection

  type(mesh_t), target :: mesh_temp
  type(sysvector_t) :: sol_n_temp, sol_nm1_temp
  type(problem_t), target :: problem_temp
  type(vector_t) :: coords_n_temp, coords_nm1_temp

! type definition for the refinement fields

  type(refinement_fields_t) :: refinement_fields

! variables

  integer :: step=0, ipost=0, iter, i

  real(dp) :: cmax, cdiff, c_B, surf_ten, phi_di
  real(dp) :: velxmax, velymax, velxdiff, velydiff, mumax, mudiff
  real(dp), allocatable, dimension(:,:) :: refinement_coor
  logical, allocatable, dimension(:) :: elems_allowed, elems_refine
  real(dp), allocatable, dimension(:,:) :: meshcoor_n, meshcoor_nm1
  real(dp), dimension(:), allocatable :: init_area, areav
  real(dp), dimension(:), allocatable :: init_aspect_ratio, aspect_ratio
  real(dp) :: norm_area, norm_aspect_ratio, time

  real(dp) :: up(nparticles,3), unm1(nparticles,3), un(nparticles,3)

  real(dp) :: dx_box, dx_part, dx_di, distmin_betw
  real(dp) :: distmin_int, distmax_int, distmin_part, distmax_part, distmax_betw

  character(len=199) :: filename
  logical :: remeshing=.false.

! allocate memory for the particles / droplets

  allocate(cdrop(ndrops,2),rdrop(ndrops))
  allocate(xp(nparticles,2),xpn(nparticles,2))
  allocate(rp(nparticles))

! pass some variables to the modulus

  rdrop = radiusdrop
  rp = radiuspart
  xi = sqrt(kappa/alpha)

! element sizes

  dx_box = 1._dp*rdrop(1)         ! element spacing on the external boundaries
  dx_part = pi*2*rp(1)/nelem_part ! element spacing around the particle
  dx_di = xi/nelem_int            ! element spacing around the interface

! sizes of the refined areas

  distmin_int = 2*xi
  distmax_int = 20*xi
  distmin_part = rp(1)+rp(1)/2
  distmax_part = rp(1)+10*rp(1)
  distmin_betw = rp(1)/4
  distmax_betw = 3*rp(1)

! pass value to functions_m

  cdrop(1,:) = 0
  eta1l = eta1
  eta2l = eta2

! spread the particle over the surfaces of the drops

  xp=0; xpn=0; up=0; un=0; unm1=0 ! initialize to zero

  allocate ( refinement_coor(nparticles,2) )  ! temporary use
  call coors_circle(refinement_coor,1)
  xp = refinement_coor
  deallocate ( refinement_coor )

! determine the wetting potential phi

  c_B = sqrt(alpha/beta)
  surf_ten = rho*(2*sqrt(2.0_dp)/3) * (kappa * c_B**2) / xi

  phi_di = (surf_ten/rho) * cos(theta_c*2*pi/360._dp)/(2*(c_B-(c_B**3)/3.0_dp))


! determine the capillary number

  write ( *, '(/1X,A,F9.3/)' ) 'Ca = ',shear_rate*eta1/(surf_ten/rdrop(1))

! for the initial mesh, the refinement points have to be given manually

  allocate ( refinement_coor(n_init,2) )

  do i = 1, size(cdrop,1)

    call coors_circle ( refinement_coor, i )

    call add_refinement_field ( refinement_fields, coor=refinement_coor, &
      distmin=distmin_int, distmax=distmax_int, dx_fine=dx_di, dx_coarse=dx_box)

!   add the periodic refinement fields
    refinement_coor(:,1) = refinement_coor(:,1) + lx
    call add_refinement_field ( refinement_fields, coor=refinement_coor, &
      distmin=distmin_int, distmax=distmax_int, dx_fine=dx_di, dx_coarse=dx_box)

    refinement_coor(:,1) = refinement_coor(:,1) - 2*lx
    call add_refinement_field ( refinement_fields, coor=refinement_coor, &
      distmin=distmin_int, distmax=distmax_int, dx_fine=dx_di, dx_coarse=dx_box)

  end do

  deallocate ( refinement_coor )

! fill coefficients for the flow problem

  call create_coefficients ( coefficients, ncoefi=300, ncoefr=250 )

  coefficients%i(:)  = 0
  coefficients%i(1:15)  = [ uintpl,       pintpl,  0,  0,     0, &
                            physqvel, physqpress,  0,  0, gauss, &
                            gaussb,            0,  0,  0,     0 ]

  coefficients%i(48) = 1 ! use mesh velocity for ALE formulation

  coefficients%r = 0
  coefficients%i(251) = 3 ! eta is function of composition

! fill coefficients for the diffuse interface problem

  coefficients%r(101:106)  = [ deltat,  Mcoef,  alpha,  beta,  kappa, rho ]
  coefficients%r(108) = phi_di   ! wetting potential
  coefficients%i(151) = timeint1 ! Euler integration of the di equations
  coefficients%i(152) = 2        ! Newton-Raphson iteration for di equation
  coefficients%i(153) = 2        ! use the iteration velocity u_i in the
                                 ! DI-convection term, instead of the
                                 ! prediction uhat =2*u_nm1 - un
                                 ! (diffuse_interface_elem1)
  call generate_read_mesh

! define the up and di problem, create (sys)vectors, subscripts, oldvectors
! and initialize sysvectors

  call define_problems_create_vectors ( init=.true. )

! fill the initial composition field (copied to sol_np1 in build_solve_up_di)

  call fill_sysvector ( mesh, problem, sol_n, physq=1, &
    node1=1, node2=mesh%nnodes, func=difunc, funcnr=1 )

! build and solve the di-up system

  call build_solve_up_di

! extract info on particle

  do i = 1, nparticles
    call get_sysvector_constraint ( mesh, problem, sol_np1, constraint=i, &
      addunknowns=.true., u=up(i,:) )
  end do

  time  = 0 ! initial time

  write ( *, '(1X,A)' ) '****************************************'
  write ( *, '(1X,A,I0)' ) 'step = ', step
  write ( *, '(1X,A,F10.6)' ) 'time = ', time
  write ( *,* )
  !do i = 1, nparticles
  !  write ( *, '(1X,A,I0)' ) 'particle ',i
  !  write ( *, '(1X,A,6F10.6)' ) 'xp = ', xp(i,:)
  !  write ( *, '(1X,A,6F10.6)' ) 'up = ', up(i,:)
  !end do

! write VTK files

  call postprocessing ( mesh, sol_np1, ipost=ipost )
  ipost = ipost + 1

! delete the old problems

  call delete_old_problems ( keepmesh=.true. )

! define the up and di problem, create (sys)vectors, subscripts, oldvectors
! and initialize sysvectors

  call define_problems_create_vectors ( init=.false. )

! fill the initial composition field again

  call fill_sysvector ( mesh, problem, sol_np1, physq=1, &
    node1=1, node2=mesh%nnodes, func=difunc, funcnr=1 )

! find the elements that are allowed to contain the interface

  call find_elems_allowed

! copy old values

  call copy ( sol_np1, sol_n )
  meshcoor_n = mesh%coor

! start time integration

  do step = 1, numtimesteps

    time = time + deltat

    write ( *, '(1X,A)' ) '****************************************'
    write ( *, '(1X,A,I0)' ) 'step = ', step
    write ( *, '(1X,A,F10.6)' ) 'time = ', time

!   change time integration to Gear after the first step

    if ( step > 1 ) coefficients%i(151) = timeint2


    if ( step == 1 ) then

!     advance particles with forward Euler
      do i = 1, nparticles
        un(i,:) = up(i,:)
        xpn(i,:) = xp(i,:)
        xp(i,:) = xp(i,:) + deltat * un(i,1:2)
      end do

    else

!     advance particle positions with 2nd order Adams-Bashforth
      do i = 1, nparticles
        unm1(i,:) = un(i,:)
        un(i,:) = up(i,:)
        xpn(i,:) = xp(i,:)
        xp(i,:) = xp(i,:) + deltat*(3*un(i,1:2)/2 - unm1(i,1:2)/2)
      end do

    end if

!   update mesh nodes: solve a Laplace's equation

    call update_mesh_nodes_multi ( mesh, problem_lapl, xp(:,1:2)-xpn(:,1:2) )

    call find_bounds_blocks ( mesh )

!   compute current aspect ratio

    call compute_element_aspect_ratio ( mesh, areav, aspect_ratio )

!   compute maximum normalized aspect ratio and area
!   (with respect to the initial values)

    norm_area = maxval( abs(log(areav/init_area)) )
    norm_aspect_ratio = maxval( abs(log(aspect_ratio/init_aspect_ratio)) )

!   remeshing criterion

    remeshing = .false.

    if ( ( norm_area >= area_threshold .or. &
           norm_aspect_ratio >= aspect_ratio_threshold ) ) then
      remeshing = .true.
      print *,'Doing remeshing...'
    end if

!   check if the interface is in allowed elements

    call check_elems_allowed

    if ( remeshing ) then

      print *,'Remeshing and projection...'

!     find refinement points

      call find_refinement_points

!     create vectors to store the mesh coordinates, these coordinates will
!     be used in the projection problem to find the new mesh at t_n and t_nm1

      call create ( problem, coords_n_temp, vec=physqvel )
      call create ( problem, coords_nm1_temp, vec=physqvel )

      coords_n_temp%u = reshape ( transpose ( meshcoor_n ), [size(meshcoor_n)] )
      coords_nm1_temp%u = reshape ( transpose ( meshcoor_nm1 ), &
                                    [size(meshcoor_nm1)] )

!     copy current definitions to temporary definitions used for the projection
!     problem

      call copy ( mesh, mesh_temp )
      call fill_mesh_parts ( mesh_temp )
      call copy ( sol_n, sol_n_temp )
      call copy ( sol_nm1, sol_nm1_temp )
      call problem_definition ( input_probdef, mesh, problem_temp )

!     delete the old problems

      call delete_old_problems ( .false. )

!     generate new mesh

      call generate_read_mesh

!     define the up and di problem, create (sys)vectors, subscripts, oldvectors
!     and initialize sysvectors

      call define_problems_create_vectors ( init=.false. )

!     project the old solutions onto the new mesh and fill in the values

      call project_old_solutions_P2

!     delete the temporary definitions

      call delete ( mesh_temp )
      call delete ( sol_n_temp )
      call delete ( sol_nm1_temp )
      call delete ( problem_temp )
      call delete ( coords_n_temp, coords_nm1_temp )

    end if

!   find the mesh velocity using a backwards differencing scheme

    if ( step == 1  ) then
      meshvel%u = reshape ( transpose ( &
                ( mesh%coor - meshcoor_n ) / deltat ), [mesh%ndim*mesh%nnodes] )
    else
      meshvel%u = reshape ( transpose ( &
        ( 1.5_dp*mesh%coor - 2*meshcoor_n + 0.5_dp*meshcoor_nm1 ) / deltat ), &
             [mesh%ndim*mesh%nnodes] )
    end if

!   build and solve the di-up system

    call build_solve_up_di

!   find the elements that are allowed to contain the interface

    if ( remeshing ) call find_elems_allowed

!   extract info on particle

    do i=1,nparticles
      call get_sysvector_constraint ( mesh, problem, sol_np1, constraint=i, &
        addunknowns=.true., u=up(i,:) )
    end do

    !do i = 1, nparticles
    !  write ( *, '(1X,A,I0)' ) 'particle ',i
    !  write ( *, '(1X,A,6F10.6)' ) 'xp = ', xp(i,:)
    !  write ( *, '(1X,A,6F10.6)' ) 'up = ', up(i,:)
    !end do

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
    meshcoor_nm1 = meshcoor_n
    meshcoor_n = mesh%coor

  end do

! delete all data including all allocated memory

  call delete_old_problems ( .false. )

contains

! generate and read mesh

  subroutine generate_read_mesh

    integer :: i

!   add refinements on and between particles

    call add_refinement_on_and_between_particles

    open ( unit=25, file='mesh.geo' )

    write ( 25, '(1X,A,F18.14,A)' ) 'ox = ', -lx/2, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'oy = ', -ly/2, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'lx = ', lx, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'ly = ', ly, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_box = ', dx_box, ';'
    write ( 25, '(1X,A,I0,A)' ) 'nobj = ', nparticles, ';'

    do i = 1, nparticles
      write ( 25, '(1X,A,I0,A,F18.14,A)' ) 'xp[', i, '] = ', xp(i,1), ';'
      write ( 25, '(1X,A,I0,A,F18.14,A)' ) 'yp[', i, '] = ', xp(i,2), ';'
      write ( 25, '(1X,A,I0,A,F18.14,A)' ) 'rp[', i, '] = ', rp(i), ';'
    end do

    write ( 25, '(1X,A,F18.14,A)' ) 'dx_part = ', dx_part, ';'

    call write_refinement_fields ( refinement_fields, 'mesh.geo' )

    write ( 25, '(/1x,a)' ) 'Include "particles_in_a_box_2D.igo";'

    write ( 25, '(/1x,a)' ) 'Include "refinement.igo";'

    close ( 25 )

    call execute_command_line ( 'gmsh -2 -order 2 -algo front2d -o mesh.msh &
                  &mesh.geo > outputmesh.out' )

!   refinement points served their purpose
    call delete_refinement_fields ( refinement_fields )

!   read mesh generated by gmsh
    call read_mesh_gmsh ( mesh, filename='mesh.msh', ndim=2, &
      physgeom=.true. )

    call add_to_mesh ( mesh, matchingcurve=[4,2], replace=4, &
      displacement=[-lx,0._dp] )

    call fill_mesh_parts ( mesh )

!   allocate some arrays
    allocate ( elems_allowed(mesh%nelem) )
    allocate ( elems_refine(mesh%nelem) )

!   allocate aspect ratio arrays
    allocate ( init_area(mesh%nelem), areav(mesh%nelem) )
    allocate ( init_aspect_ratio(mesh%nelem), aspect_ratio(mesh%nelem) )

!   define some arrays for the ALE mesh position at old times
    allocate ( meshcoor_n(mesh%nnodes,mesh%ndim), &
      meshcoor_nm1(mesh%nnodes,mesh%ndim) )

!   compute initial element aspect ratio
    call compute_element_aspect_ratio ( mesh, init_area, init_aspect_ratio )

    meshcoor_n=0; meshcoor_nm1=0 ! initialize to zero
    norm_area=0; norm_aspect_ratio = 0

    call write_mesh_vtk ( mesh, filename='mesh.vtk' )

  end subroutine generate_read_mesh

! compute the aspect ratio and area for each mesh element

  subroutine compute_element_aspect_ratio ( mesh, area, asp_ratio )

    type(mesh_t), intent(in) :: mesh
    real(dp), dimension(:), intent(out) :: asp_ratio, area

    integer, parameter :: elgrp=1
    integer :: elem, node(3)
    real(dp) :: la, lb, lc, sper
    real(dp) :: vert(3,2)

    do elem = 1,mesh%nelem

!     get coordinates of triangle vertices (local node numbers: 1, 3, 5)
      node = [1,3,5]
      vert = mesh%coor(mesh%topology(elgrp)%a(node,elem),:)

!     compute side lengths
      la = sqrt( (vert(1,1) - vert(2,1))**2 + (vert(1,2) - vert(2,2))**2 )
      lb = sqrt( (vert(2,1) - vert(3,1))**2 + (vert(2,2) - vert(3,2))**2 )
      lc = sqrt( (vert(3,1) - vert(1,1))**2 + (vert(3,2) - vert(1,2))**2 )

!     compute triangle semiperimeter and area
      sper = (la + lb + lc)/2.0_dp
      area(elem) = sqrt( sper*(sper - la)*(sper - lb)*(sper - lc) )

!     compute aspect ratio
      asp_ratio(elem) = max(la,lb,lc)**2/area(elem)

    end do

  end subroutine compute_element_aspect_ratio

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
      distmin=distmin_int, distmax=distmax_int, dx_fine=dx_di, dx_coarse=dx_box)

!   add the periodic refinement fields
    refinement_coor(:,1) = refinement_coor(:,1) + lx
    call add_refinement_field ( refinement_fields, coor=refinement_coor, &
      distmin=distmin_int, distmax=distmax_int, dx_fine=dx_di, dx_coarse=dx_box)

    refinement_coor(:,1) = refinement_coor(:,1) - 2*lx
    call add_refinement_field ( refinement_fields, coor=refinement_coor, &
      distmin=distmin_int, distmax=distmax_int, dx_fine=dx_di, dx_coarse=dx_box)

    call delete(c)

    deallocate ( refinement_coor )

  end subroutine find_refinement_points

! start the problems

  subroutine define_problems_create_vectors ( init )

    logical, intent(in) :: init

    integer :: i

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

!   essential boundary conditions for velocity/pressure
    call define_essential ( mesh, input_probdef, curve1=1, physq=physqvel )
    call define_essential ( mesh, input_probdef, curve1=3, physq=physqvel )
    call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

!   constraint on particle boundaries (freely floating)
    do i = 1, nparticles
      call define_constraint ( mesh, input_probdef, curve1=i+4, &
        physq=physqvel, discretization='collocation', nodedof=2, &
        naddunknowns=3 )
    end do

!   periodic velocities
    call define_constraint ( mesh, input_probdef, &
      physq=physqvel, curve1=2, curve2=4, discretization='collocation', &
      nodedof=2, excludecurves=[1,3] )

!   periodic mu
    call define_constraint ( mesh, input_probdef, &
      physq=physqmu, curve1=2, curve2=4, discretization='collocation' )

!   Initial timestep: c is known; remaining time steps: c periodic
    if ( init ) then
      call define_essential ( mesh, input_probdef, elgroup1=1, physq=physqc )
    else
      call define_constraint ( mesh, input_probdef, &
        physq=physqc, curve1=2, curve2=4, discretization='collocation' )
    end if

    call problem_definition ( input_probdef, mesh, problem )

!   create system vectors for the up problem
    call create_sysvector ( problem, sol_np1, sol_n, sol_nm1, rhsd )
    call create_sysvector ( problem, sol_i )

!   create mesh velocity vectors and vector to store mesh coordinates
    call create ( problem, meshvel, vec=physqvel )

    meshvel%u = 0 ! initialize to zero

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

!   create system matrix for the upg problem
    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric=.false. )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )
    call create_sysmatrix_data ( sysmatrix )

!   create oldvectors
!   NOTE: two problems are supplied for backwards compatibility with the
!   diffuse_interface_elem1 routine
    call create_oldvectors ( oldvectors, nsysvec=4, nprob=2, nelvec=1, &
      nvec=2 )
      oldvectors%s(1)%p => null()
      oldvectors%s(2)%p => sol_n
      oldvectors%s(3)%p => sol_np1
      oldvectors%s(4)%p => sol_nm1
      oldvectors%p(1)%p => problem
      oldvectors%p(2)%p => problem
      oldvectors%e(1)%p => viscosity
      oldvectors%v(1)%p => meshvel
      oldvectors%v(2)%p => composition

!   initialize all vectors to zero
    sol_nm1%u = 0.0_dp
    sol_n%u = 0.0_dp
    sol_np1%u = 0.0_dp

  end subroutine define_problems_create_vectors

! delete old problem

  subroutine delete_old_problems ( keepmesh )

    logical, intent(in) :: keepmesh

    call delete ( problem )
    call delete ( sysmatrix )
    call delete ( input_probdef )
    call delete ( sol_np1, sol_n, sol_nm1, sol_i )
    call delete ( rhsd )
    call delete ( oldvectors )
    call delete ( velx, vely, pval )
    call delete ( cval, muval )
    call delete ( viscosity )
    call delete ( composition )
    call delete ( meshvel )

    if ( problem_lapl%created ) call delete ( problem_lapl )

    if ( .not. keepmesh ) then
      call delete ( mesh )
      deallocate ( elems_allowed )
      deallocate ( elems_refine )
      deallocate ( meshcoor_n, meshcoor_nm1 )
      deallocate ( init_aspect_ratio, aspect_ratio )
      deallocate ( init_area, areav )
      if ( meshup%meshgen ) call delete ( meshup )
      if ( meshdown%meshgen ) call delete ( meshdown )
    end if

  end subroutine delete_old_problems

! build and solve the velocity-pressure-c-mu problem

  subroutine build_solve_up_di

    integer :: i

    print *,'Starting iterations'

!   fill the viscosity in the gauss points ( using a prediction for the
!   composition: chat = 2*c_n - c_nm1 )
    if ( step == 0 .or. step == 1 ) then
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
      degfd=1, value=-ly*shear_rate/2 )
    call fill_sysvector ( mesh, problem, sol_np1, curve1=1, physq=physqvel, &
      degfd=2, value=0.0_dp )
    call fill_sysvector ( mesh, problem, sol_np1, curve1=3, physq=physqvel, &
      degfd=1, value=ly*shear_rate/2 )
    call fill_sysvector ( mesh, problem, sol_np1, curve1=3, physq=physqvel, &
      degfd=2, value=0.0_dp )

!   presssure
    call fill_sysvector ( mesh, problem, sol_np1, &
      point=1, physq=physqpress, value=0._dp )

    iter = 0

    newton: do

      iter = iter + 1

!     Stokes flow : velocity / pressure
      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=generalized_stokes_elem, oldvectors=oldvectors, &
         coefficients=coefficients, physqrow=[physqvel,physqpress], &
         physqcol=[physqvel,physqpress] )

!     constraint on particle boundary (freely floating)

      do i = 1, nparticles
        ipart=i
        call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
          constraint1=i, elemsub=elementc, addmatvec=.true., &
          coefficients=coefficients )
      end do

!     periodical condition on velocities
      call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
        constraint1=nparticles+1, elemsub=stokes_constr_node_conn, &
        addmatvec=.true., coefficients=coefficients )

!     add the div tau_c to the matrix/vector
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        elemsub=divstressdi_implicit_elem, oldvectors=oldvectors, &
        physqrow=[physqvel], physqcol=[physqc], coefficients=coefficients )

!     set additional block to zero if stress form is used
      call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
        buildvector=.false., physqrow=[physqvel], &
        physqcol=[physqmu], zeromatvec=.true. )

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

!     periodical condition on c/mu
      if ( step == 0 ) then
        call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
          constraint1=nparticles+2, elemsub=stokes_constr_node_conn, &
          addmatvec=.true., coefficients=coefficients )
      else
        call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
          constraint1=nparticles+2, constraint2=nparticles+3, &
          elemsub=stokes_constr_node_conn, addmatvec=.true., &
          coefficients=coefficients )
      end if

!     add the boundary conditon for the contact angle
      do i = 1, nparticles
        call add_boundary_elements ( mesh, problem, rhsd, curve=4+i, &
          sysmatrix=sysmatrix, elemsub=di_natboun_curve, &
          physq=[physqc,physqmu], coefficients=coefficients, &
          oldvectors=oldvectors )
      end do

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

      if ( step == 0 ) exit newton

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

!   derive the pressure in all nodes
    call create_oldvectors ( oldvectors_post, nsysvec=1 )
    oldvectors_post%s(1)%p => sol

    call create_vector ( problem, pressure, vec=1 )

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

!   transfer the coordinates at n
    call transfer_data ( mesh_temp, problem_temp, &
      vector1=coords_n_temp, vector2=vec_p, degfd2=[6,7] )

!   transfer the coordinates at n-1
    call transfer_data ( mesh_temp, problem_temp, &
      vector1=coords_nm1_temp, vector2=vec_p, degfd2=[8,9] )

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
    sol_n%u(velx%s) = sol_proj(1)%u
    sol_n%u(vely%s) = sol_proj(2)%u

!   update old composition
    sol_n%u(cval%s) = sol_proj(3)%u
    sol_nm1%u(cval%s) = sol_proj(4)%u
    sol_n%u(muval%s) = sol_proj(5)%u

!   update the old coordinates
    meshcoor_n(:,1) = sol_proj(6)%u
    meshcoor_n(:,2) = sol_proj(7)%u
    meshcoor_nm1(:,1) = sol_proj(8)%u
    meshcoor_nm1(:,2) = sol_proj(9)%u

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

! subroutine to give coordinates on a circle

  subroutine coors_circle ( coor, idrop )

    real(dp), dimension(:,:), intent(inout) :: coor
    integer, intent(in) :: idrop

    integer :: np, i
    real(dp) :: p(size(coor,1))

    np = size(coor,1)
    p = [(2*pi*(i-1)/np,i=1,np)]

    coor(:,1) = rdrop(idrop)* cos(p) + cdrop(idrop,1)
    coor(:,2) = rdrop(idrop)* sin(p) + cdrop(idrop,2)

  end subroutine coors_circle

! add refinements between close particles

  subroutine add_refinement_on_and_between_particles

    real(dp) :: part_dist, coor(1,2), vec(2)
    integer :: i, j

    real(dp), allocatable :: lxp(:,:), lrp(:)

!   add the periodic particles (to ensure refinement is periodic)

    allocate(lxp(nparticles*3,2))
    allocate(lrp(nparticles*3))

    lxp(:nparticles,:) = xp
    lxp(nparticles+1:2*nparticles,1) = xp(:,1) - lx
    lxp(nparticles+1:2*nparticles,2) = xp(:,2)
    lxp(2*nparticles+1:,1) = xp(:,1) + lx
    lxp(2*nparticles+1:,2) = xp(:,2)

    lrp(:nparticles) = rp
    lrp(nparticles+1:2*nparticles) = rp
    lrp(2*nparticles+1:) = rp

!   add refinement for the particles

    call add_refinement_field ( refinement_fields, &
      coor=lxp, distmin=distmin_part, distmax=distmax_part, dx_fine=dx_part,&
      dx_coarse=dx_box )

!   find the minimal distance between every particle (periodic particles
!   included)

    do i=1,size(lxp,1)-1

      do j=i+1,size(lxp,1)
         part_dist = &
           sqrt(dot_product(lxp(i,:)-lxp(j,:),lxp(i,:)-lxp(j,:)))-lrp(i)-lrp(j)

!       check the number of elements between the particles
        if ( part_dist/dx_part < min_elem ) then

!         place coordinate exactly between the particle boundaries
          vec = (lxp(j,:)-lxp(i,:)) ! connecting vector
          vec = vec / sqrt ( vec(1)**2 + vec(2)**2 ) ! normalized
          coor(1,:) = lxp(i,:) + vec*lrp(i) + vec*part_dist/2

          call add_refinement_field ( refinement_fields, coor=coor, &
            distmin=distmin_betw, distmax=distmax_betw, &
            dx_fine=part_dist/min_elem, dx_coarse=dx_box )

        end if

      end do

    end do

    deallocate(lxp,lrp)

  end subroutine add_refinement_on_and_between_particles

end program dimpart1
