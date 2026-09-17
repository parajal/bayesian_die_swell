! Viscoelastic extrudate swell problem
! 3D flow for a Giesekus model using a b-tensor formulation
! DEVSS-G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! second-order Gear time integration

program extrudate_swell3d_ba

  use tfem_m
  use kind_defs_m
  use mesh_m
  use meshgen_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use hsl_ma57_m
  use io_utils_m
  use surface_advection_elements_m
  use update_mesh_nodes_bc_m
  use figplot_m
  use projection_elements_m
  use limits_m
!  use metis5_m
  use math_defs_m

  implicit none


! constants flow problem

  integer, parameter :: &
    uintpl = 6,     & ! P2 velocities
    pintpl = 2,     & ! P1 pressures
    gintpl = 2,     & ! P1 gradients
    bintpl = 2,     & ! P1 conformation
    ointpl = 6,     & ! P2 shape of object elements
    hintpl = 6,     & ! Interpolation for the height
    physqgrad = 1,  & ! physical quantity nr of the gradients
    physqvel = 2,   & ! physical quantity nr of the velocities
    physqpress = 3, & ! physical quantity nr of the pressures
    gauss = 5,      & ! 6-point Gauss integration of tets
    gaussb = 5,     & ! 3-point integration of boundary elements
    ncompb = 9,     & ! number of b-tensor components
    ncompc = 6,     & ! number of conformation tensor components
    bvariant = 1,   & ! b-formulation: 1:CDT 2:symmetric 4:Cholesky with log
    nmodes = 1,     & ! number of modes
    startm = 501,   & ! start of material model data
    coorsys = 0,    & ! planar Cartesian coordinate system
    ncompv_P2 = 6,  & ! number of components for P2 projection:
                      ! x,y-coordinates at n and n-1
    ncompv_P1 = 12, &   ! number of components for P1 projection:
                      ! conformation components at n and n-1
    model = 3         ! 2: UCM/Oldroyd-B, 3: Giesekus, 5: linear PTT

! constants surface advection

  integer, parameter :: &
    ninti_sf_adv = 6,   & ! number of Gauss points
    method = 1            ! discretization method 0: Galerkin, 1: SUPG


! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefb, input_probdefc_projc
  type(problem_t), target :: problem, problemb, problemc_projc
  type(sysmatrix_t) :: sysmatrix, sysmatrixb, sysmatrixc_projc
  type(sysvector_t), target :: sol, sol_n, sol_nm1
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(sysvector_t), dimension(ncompb,nmodes), target :: solb, solb_n, solb_nm1
  type(sysvector_t), dimension(ncompb,nmodes) :: rhsb
  type(sysvector_t), dimension(ncompc,nmodes) :: solc
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc_projc, rhsc_projc
  type(lu_ma41_t) :: lub
  type(solver_options_ma41_t) :: solver_options_u, solver_options_b

  type(refinement_fields_t) :: refinement_fields
  type(vector_t) :: ctensor, btensor, ctensor_in, btensor_in
  type(subscript_t) :: velx, vely, velz, bval, bval_in
  type(subscript_t) :: hgt, hgt_end
  type(subscript_t) :: velx_inlet, vely_inlet, velz_inlet
  type(subscript_t) :: bb_in(ncompb), bb_inlet(ncompb)
  real(dp) :: phi
  type(sysvector_t), target :: sol_p

! 2D height function for surface advection
  type(mesh_t) :: mesh_surf1
  type(input_probdef_t) :: input_probdef_surf1
  type(problem_t) :: problem_surf1
  type(sysmatrix_t) :: sysmatrix_surf1
  type(sysvector_t), target :: sol_surf1, sol_surf1_n, sol_surf1_nm1
  type(sysvector_t) :: rhsd_surf1
  type(oldvectors_t) :: oldvectors_surf1
  type(coefficients_t) :: coefficients_sf_adv
  type(vector_t), target :: velocity_surf1
  type(subscriptvec_t) :: mvelx_surf1, mvely_surf1, mvelz_surf1
  type(subscript_t) ::  velxsurf1, velysurf1, velzsurf1
  type(subscriptvec_t) :: velx_surf1, vely_surf1, velz_surf1
  type(subscriptvec_t) :: meshvelx_surf1, meshvely_surf1
  type(subscript_t) :: subs_surf1, pos1_surf1, pos2_surf1
  type(subscript_t) :: subs_curve1

  type(vector_t), target :: meshvel, meshvel_surf1
  type(problem_t) :: problem_lapl_surf, problem_laply_surf, problem_laplz_surf
  real(dp), allocatable, dimension(:,:) :: meshcoor_initial

! variables

  integer :: &
    timeint1 = 1,         & ! (first-order) Euler time integration (first step)
    timeint2 = 7,         & ! (second-order) semi-implicit Gear time integration
    numtimesteps = 2500,  & ! number of time steps
    step0 = 0,            & ! initial step number
    htype = 2,            & ! upwind parameter
    Uscaling = 3,         & ! upwind parameter
    vtkevery = 1            ! write vtk every .. steps

  real(dp), parameter :: &
    eta_p = 1.0_dp,     & ! polymer viscosity
    G = 1.0_dp,         & ! polymer modulus
    betav = 0.59_dp,    & ! viscosity ratio
    mobility = 0.01_dp, & ! mobility parameter in the Giesekus model (alpha)
    lambda = 1._dp        ! maximum polymer relaxation time for each mode

 real(dp), parameter :: &
    eta_s = betav/(1-betav)*lambda*G,    & ! solvent viscosity
    deltat = 5.e-3_dp,      & ! time step
    time0 = 0._dp,          & ! initial time
    U_avg = 1.0_dp,         & ! average velocity at the entry
    Ro = 1.0_dp               ! radius

  real(dp) :: &
    beta = 1.0_dp,         & ! upwinding parameter in the SUPG method
    beta_sf_adv = 0.5_dp,  & ! upwinding parameter in the SUPG method (interface)
    rmin = 0.1_dp,         & ! refinement parameter
    theta = 45._dp*2*pi/360, & !angle
    dx_box = 0.2_dp,         & ! element spacing on the external boundaries
    dx_wall = 0.1_dp,       & ! element spacing on upper boundary and die exit
    dx_inlet = 0.2_dp,       & ! element spacing for inlet problem
    rs_gup = 8.5_dp,         & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 8.5_dp,         & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 8.5_dp,          & ! real_storage for the conformation LU (HSL)
    is_c  = 8.6_dp             ! integer_storage for conformation LU (HSL)

  integer :: step=1, i, m, ipost=0
  integer :: ipost_in=0
  real(dp) :: alpha, time, r
  real(dp) :: flowrate

  real(dp), allocatable, dimension(:) :: initial_hsurf1
  real(dp), allocatable, dimension(:,:) :: bd
  real(dp) :: ox=3._dp, L2=5._dp, oy=0.125_dp, oz=0.375_dp

  real(dp), allocatable, dimension(:,:) :: meshcoor_n, meshcoor_nm1
  real(dp), allocatable, dimension(:,:) :: meshcoorsurf1_n, meshcoorsurf1_nm1

  real(dp), allocatable, dimension(:,:) :: surf1hat, surf1hatn

  character(len=300) :: filename

  logical :: cproj=.true.

! All parameters for inlet problem

! definitions

  type(mesh_t) :: mesh_inlet
  type(input_probdef_t) :: input_probdef_inlet, input_probdefb_inlet, &
                           input_probdefc_proj_inlet
  type(problem_t), target :: problem_inlet, problemb_inlet, problemc_proj_inlet
  type(sysmatrix_t) :: sysmatrix_inlet, sysmatrixb_inlet, sysmatrixc_proj_inlet
  type(sysvector_t), target :: sol_inlet, solm1_inlet
  type(sysvector_t) :: rhsd_inlet
  type(oldvectors_t) :: oldvectors_ve_inlet
  type(coefficients_t) :: coefficients_inlet
  type(subscript_t) :: velx_in, vely_in, velz_in, cval_in
  type(subscriptvec_t) :: cxx, cxy, cyz, cyy, czz, cxz, &
          cxx13, cxy13, cxz13, cyy13, cyz13, czz13
  type(sysvector_t), dimension(ncompb,nmodes), target :: solb_inlet, solbm1_inlet
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc_inlet
  real(dp) :: cc_chol(ncompc), low_chol(ncompc)
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc_proj_inlet
  type(sysvector_t), dimension(ncompb,nmodes) :: rhsb_inlet
  type(sysvector_t), dimension(ncompc,nmodes) ::  rhsc_proj_inlet
  type(lu_ma41_t) :: lub_inlet
  type(lu_ma57_t) :: luc_proj_inlet
  type(solver_options_ma41_t) :: solver_options_u_in, solver_options_b_in

  ! variables

  real(dp) :: &
    rs_gup_in = 2.0_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup_in = 2.0_dp,   &   ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c_in  = 2.0_dp,    & ! rbal_storage for the conformation LU (HSL)
    is_c_in  = 2.4_dp       ! integer_storage for conformation LU (HSL)


  integer :: npar
  integer :: vertices_in(3) = [1,3,5]

  call execute_command_line ( 'rm *.vtk' )

  open (unit = 405, file = 'conf_line13.out')

! set some parameters

  alpha = eta_p  ! DEVSS parameter

  flowrate = 0.25_dp*U_avg*pi*Ro**2

  if ( model == 2 ) then
    npar = 2 ! Oldroyd-B
  else if ( any(model == [3,5,6]) ) then
    npar = 3 ! Giesekus, PTT
  end if

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=600, ncoefr=500+npar*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gaussb,   bintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      0,     0,          coorsys,  ( 0, i = 24, 600 )  &
    ]

  coefficients%i(29:31) = [ htype, 0, Uscaling ]
  coefficients%i(35) = ointpl
  coefficients%i(48) = 1 ! 1: use mesh velocity for ALE formulation
  coefficients%i(40) = 3

  coefficients%i(71) = bvariant ! variant for b-formulation
  if ( cproj ) coefficients%i(72) = 1 ! use c projection in momentum balance
  coefficients%i(90) = 1  ! rotation reinitialization on element level

  coefficients%r = 0
  coefficients%r(1:10) = &
    [ eta_s, 0._dp,  0._dp, alpha, 0._dp, &
      0._dp, 0._dp, deltat,  beta, 0._dp  &
    ]

  coefficients%r(501:500+2*nmodes) = [G, lambda]

  coefficients%r(503:502+nmodes) = [mobility] ! Giesekus, PTT

  call generate_read_mesh

! fill coefficients inlet

  call create_coefficients ( coefficients_inlet, ncoefi=150, ncoefr=500+npar*nmodes )

  coefficients_inlet%i = &
    [ uintpl, pintpl, 0, 0,      gintpl, &
      physqvel, physqpress, 0, physqgrad, gauss,  &
      gaussb,  bintpl, 0, 0,      0,      &
      0,      0,      model, nmodes, startm, &
      0,   timeint1, ( 0, i = 23, 150 )  &
    ]

  coefficients_inlet%i(23) = 0  ! coorsys = 1, cilindrical coordinates, axisymmetric
  coefficients_inlet%i(61) = 0  ! projected G=0, direct velocity gradient=1 in CE
  coefficients_inlet%i(67) = 1  ! 3D velocity (swirl)
  coefficients_inlet%i(40) = 3

  coefficients_inlet%i(71) = bvariant ! variant for b-formulation
  if ( cproj ) coefficients_inlet%i(72) = 1 ! use c projection in momentum balance
  coefficients_inlet%i(90) = 1  ! rotation reinitialization on element level


  coefficients_inlet%r(1:503) = &
    [ eta_s,    0._dp,   0._dp, alpha, 0._dp, &
      flowrate, 0._dp,  deltat, beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G, lambda, mobility  &
    ]

  call define_inlet_problem

  call define_problems_create_vectors

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

  call define_surface_problems

! time stepping

  time = 0

  coefficients%i(22) = timeint1  ! first step integration scheme

  do step = 1, numtimesteps

    time = time + deltat

    if (step==2) then

      coefficients%i(22) = timeint2
      coefficients_sf_adv%i(5) = 2  ! second-order scheme

    end if

    call solve_inlet_problem

    sol%u(velx_inlet%s) = sol_inlet%u(velz_in%s)
    sol%u(vely_inlet%s) = sol_inlet%u(velx_in%s)
    sol%u(velz_inlet%s) = sol_inlet%u(vely_in%s)

    if ( bvariant == 4 ) then

      do m = 1, nmodes

        do i = 1, size( bb_in(1)%s(:) )

          solc(1,m)%u(bb_inlet(1)%s(i)) = solc_inlet(6,m)%u(bb_in(1)%s(i))
          solc(2,m)%u(bb_inlet(1)%s(i)) = solc_inlet(3,m)%u(bb_in(1)%s(i))
          solc(3,m)%u(bb_inlet(1)%s(i)) = solc_inlet(5,m)%u(bb_in(1)%s(i))
          solc(4,m)%u(bb_inlet(1)%s(i)) = solc_inlet(1,m)%u(bb_in(1)%s(i))
          solc(5,m)%u(bb_inlet(1)%s(i)) = solc_inlet(2,m)%u(bb_in(1)%s(i))
          solc(6,m)%u(bb_inlet(1)%s(i)) = solc_inlet(4,m)%u(bb_in(1)%s(i))

          cc_chol(1) = solc(1,m)%u(bb_inlet(1)%s(i))
          cc_chol(2) = solc(2,m)%u(bb_inlet(1)%s(i))
          cc_chol(3) = solc(3,m)%u(bb_inlet(1)%s(i))
          cc_chol(4) = solc(4,m)%u(bb_inlet(1)%s(i))
          cc_chol(5) = solc(5,m)%u(bb_inlet(1)%s(i))
          cc_chol(6) = solc(6,m)%u(bb_inlet(1)%s(i))

          low_chol = chol(cc_chol)

          solb(1,m)%u(bb_inlet(1)%s(i)) = log(low_chol(1))
          solb(2,m)%u(bb_inlet(1)%s(i)) = 0._dp
          solb(3,m)%u(bb_inlet(1)%s(i)) = 0._dp
          solb(4,m)%u(bb_inlet(1)%s(i)) = low_chol(2)
          solb(5,m)%u(bb_inlet(1)%s(i)) = log(low_chol(4))
          solb(6,m)%u(bb_inlet(1)%s(i)) = 0._dp
          solb(7,m)%u(bb_inlet(1)%s(i)) = low_chol(3)
          solb(8,m)%u(bb_inlet(1)%s(i)) = low_chol(5)
          solb(9,m)%u(bb_inlet(1)%s(i)) = log(low_chol(6))

        end do

      end do

    else

      do m = 1, nmodes

        do i = 1, size( bb_in(1)%s(:) )

          solb(1,m)%u(bb_inlet(1)%s(i)) = solb_inlet(9,m)%u(bb_in(1)%s(i))
          solb(2,m)%u(bb_inlet(1)%s(i)) = solb_inlet(7,m)%u(bb_in(1)%s(i))
          solb(3,m)%u(bb_inlet(1)%s(i)) = solb_inlet(8,m)%u(bb_in(1)%s(i))
          solb(4,m)%u(bb_inlet(1)%s(i)) = solb_inlet(3,m)%u(bb_in(1)%s(i))
          solb(5,m)%u(bb_inlet(1)%s(i)) = solb_inlet(1,m)%u(bb_in(1)%s(i))
          solb(6,m)%u(bb_inlet(1)%s(i)) = solb_inlet(2,m)%u(bb_in(1)%s(i))
          solb(7,m)%u(bb_inlet(1)%s(i)) = solb_inlet(6,m)%u(bb_in(1)%s(i))
          solb(8,m)%u(bb_inlet(1)%s(i)) = solb_inlet(4,m)%u(bb_in(1)%s(i))
          solb(9,m)%u(bb_inlet(1)%s(i)) = solb_inlet(5,m)%u(bb_in(1)%s(i))

         end do

       end do

    end if

    if ( step >= 2 ) then


      meshcoor_nm1 = meshcoor_n
      meshcoor_n = mesh%coor

      surf1hatn(:,1) = mesh%coor(mesh%surfaces(4)%nodes(:),3)
      surf1hatn(:,2) = mesh%coor(mesh%surfaces(4)%nodes(:),2)
      surf1hat(:,1) = 2._dp*meshcoor_n(mesh%surfaces(4)%nodes(:),3)- &
                        meshcoor_nm1(mesh%surfaces(4)%nodes(:),3)
      surf1hat(:,2) = 2._dp*meshcoor_n(mesh%surfaces(4)%nodes(:),2)- &
                        meshcoor_nm1(mesh%surfaces(4)%nodes(:),2)

      call update_mesh_nodes_surfaces_3D ( mesh, problem_lapl_surf, problem_laply_surf, &
           problem_laplz_surf, disp1=surf1hat-surf1hatn )

      call find_bounds_blocks ( mesh )

      write (filename,'(a,i4.4,a)')'meshtime', ipost, '.vtk'
      call write_mesh_vtk ( mesh, filename=filename )

      meshvel%u = reshape ( transpose ( &
            ( 1.5_dp*mesh%coor - 2*meshcoor_n + 0.5_dp*meshcoor_nm1 ) / deltat ), &
                             [mesh%ndim*mesh%nnodes] )
    end if

!   L2-projection c=b*b^T
    if ( cproj ) then
      call solve_b_projection
    end if

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_solve_vpG

    call build_solve_ve

    call solve_4surfaces_height_corrector

!   write max and mean values of conformation and b-tensor tensor to a file

!   create a vector for conformation and b-tensor
!   tensor for post processing

    call create_vector ( problemb, ctensor, vec=3 )
    call create_vector ( problemb, btensor, vec=4 )

    call derive_vector ( mesh, problemb, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    if ( step == 1 ) then
      open(unit=12, recl=300, status='replace', file='cval.out')
    else
      open(unit=12, recl=300, position='append', file='cval.out')
    end if
    write(12, fmt=*) step * deltat, maxval(ctensor%u(cxx%s)), &
                                    maxval(ctensor%u(cxy%s)), &
                                    maxval(ctensor%u(cxz%s)), &
                                    maxval(ctensor%u(cyy%s)), &
                                    maxval(ctensor%u(czz%s)), &
                                    sum(ctensor%u(cxx%s))/size(cxx%s), &
                                    sum(ctensor%u(cxy%s))/size(cxy%s), &
                                    sum(ctensor%u(cxz%s))/size(cxz%s), &
                                    sum(ctensor%u(cyy%s))/size(cyy%s), &
                                    sum(ctensor%u(czz%s))/size(czz%s)
    close(unit=12)

    call delete ( ctensor, btensor )

    print '(i6,4es16.8)', step0+step, time0+step*deltat, maxval(sol_surf1%u)

    if (step ==1 ) then
      call postprocessing ( mesh, sol, ipost=ipost )
      ipost = ipost + 1
    end if

!   write VTK files
    if ( vtkevery > 0 ) then
      if ( mod(step,vtkevery) == 0 ) then
        call postprocessing ( mesh, sol, ipost )
        ipost = ipost + 1
      end if
    end if

  end do

  close(unit=11)

  call delete_old_problems

contains

! generate and read mesh

  subroutine generate_read_mesh

    integer :: i, node

    integer :: surface
    real(dp) :: dr2(2)

    call write_gmsh_parameters ( ox, oy, oz, L2, dx_box, dx_wall, dx_inlet)

    call execute_command_line ( 'gmsh -3 -order 2 -o mesh.msh &
                  &mesh.geo > outputmesh.out' )

!   read mesh generated by gmsh
    call read_mesh_gmsh ( mesh, filename='mesh.msh', ndim=3, &
      physgeom=.true. )

!   move quadratic nodes to surface of cylinder
    do surface = 3,4
      do i = 1, mesh%surfaces(surface)%nnodes
        node = mesh%surfaces(surface)%nodes(i)
        r = sqrt(dot_product(mesh%coor(node,2:3),mesh%coor(node,2:3)))
        dr2 = mesh%coor(node,2:3)        ! vector from origin to current node
        dr2 = dr2/sqrt(dot_product(dr2,dr2)) ! normalize
        dr2 = (Ro - r) * dr2                 ! add wanted displacement
        mesh%coor(node,2:3) = mesh%coor(node,2:3) + dr2
      end do
    end do

!   add volume made from mesh elements for integral p dV =0

    call add_to_mesh ( mesh, volumefromgroups=[1] )
    call add_to_mesh ( mesh, curve=[11,12] ) ! curve13

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

  subroutine write_gmsh_parameters ( ox, oy, oz, L2, dx_box, dx_wall, dx_inlet )

    real(dp), intent(in) :: ox, oy, oz, L2, dx_box, dx_wall, dx_inlet

    call add_refinement_field ( refinement_fields, coor=reshape([0._dp,Ro*cos(0.25*theta),Ro*sin(0.25*theta)],[1,3]), &
     distmin=4*rmin, distmax=8*rmin, dx_fine=dx_wall, dx_coarse=dx_box )
    call add_refinement_field ( refinement_fields, coor=reshape([0._dp,Ro*cos(0.5*theta),Ro*sin(0.5*theta)],[1,3]), &
     distmin=4*rmin, distmax=8*rmin, dx_fine=dx_wall, dx_coarse=dx_box )
    call add_refinement_field ( refinement_fields, coor=reshape([0._dp,Ro*cos(0.75*theta),Ro*sin(0.75*theta)],[1,3]), &
     distmin=4*rmin, distmax=8*rmin, dx_fine=dx_wall, dx_coarse=dx_box )
    call add_refinement_field ( refinement_fields, coor=reshape([0._dp,Ro*cos(theta),Ro*sin(theta)],[1,3]), &
     distmin=4*rmin, distmax=8*rmin, dx_fine=dx_wall, dx_coarse=dx_box )
    call add_refinement_field ( refinement_fields, coor=reshape([0._dp,Ro*cos(1.25*theta),Ro*sin(1.25*theta)],[1,3]), &
     distmin=4*rmin, distmax=8*rmin, dx_fine=dx_wall, dx_coarse=dx_box )
    call add_refinement_field ( refinement_fields, coor=reshape([0._dp,Ro*cos(1.5*theta),Ro*sin(1.5*theta)],[1,3]), &
     distmin=4*rmin, distmax=8*rmin, dx_fine=dx_wall, dx_coarse=dx_box )
    call add_refinement_field ( refinement_fields, coor=reshape([0._dp,Ro*cos(1.75*theta),Ro*sin(1.75*theta)],[1,3]), &
     distmin=4*rmin, distmax=8*rmin, dx_fine=dx_wall, dx_coarse=dx_box )
    call add_refinement_field ( refinement_fields, coor=reshape([0._dp,Ro,0._dp],[1,3]), &
     distmin=4*rmin, distmax=8*rmin, dx_fine=dx_wall, dx_coarse=dx_box )
    call add_refinement_field ( refinement_fields, coor=reshape([0._dp,0._dp,Ro],[1,3]), &
     distmin=4*rmin, distmax=8*rmin, dx_fine=dx_wall, dx_coarse=dx_box )

    open ( unit=25, file='mesh.geo' )

    write ( 25, '(1X,A,F18.14,A)' ) 'ox = ', ox, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'oy = ', oy, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'oz = ', oz, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'L2 = ', L2, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_box = ', dx_box, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'y[1] = ', Ro*cos(theta), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'z[1] = ', Ro*sin(theta), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'y[2] = ', Ro, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'z[2] = ', Ro, ';'

    write ( 25, '(1X,A,F18.14,A)' ) 'dx_wall = ', dx_wall, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_inlet = ', dx_inlet, ';'

!   write the refinement fields (the file needs to be open for writing)
    call write_refinement_fields ( refinement_fields, 'mesh.geo' )

    write ( 25, '(/1x,a)' ) 'Include "mesh_3D_bc.igo";'

    close ( 25 )

!   refinement points served their purpose
    call delete_refinement_fields ( refinement_fields )

  end subroutine write_gmsh_parameters

  subroutine define_problems_create_vectors

! problem definition of gradient/velocity/pressure

    integer :: i, m

!   problem definition of gradient/velocity/pressure
    call create_input_probdef ( mesh, input_probdef, nvec=7, nphysq=3 )

    input_probdef%vec_elementdof(1)%a =   &
        reshape ( [ 9,0,9,0,9,0,0,0,0,9,    &  ! G
                     3,3,3,3,3,3,3,3,3,3,    &  ! velocity
                     1,0,1,0,1,0,0,0,0,1,    &  ! pressure
                     1,1,1,1,1,1,1,1,1,1,    &  ! scalar, such as vorticity
                     6,6,6,6,6,6,6,6,6,6,    &  ! tensor
                     [(ncompv_P2,i=1,10)],   &  ! projection
                     1,1,1,1,1,1,1,1,1,1 ], &  ! effective shear rate
                     [10,7] )

    input_probdef%physq = [1,2,3]
    input_probdef%probnr = 1

! Dirichlet boundary conditions

! outflow
   call define_essential ( mesh, input_probdef, &
      surface1=2, physq=physqvel, degfd=[0,1,1] )

! inflow
    call define_essential ( mesh, input_probdef, &
      surface1=1, physq=physqvel )

! wall
    call define_essential ( mesh, input_probdef, &
      surface1=3, physq=physqvel, excludesurfaces=[1] )
    call define_essential ( mesh, input_probdef, &
      surface1=6, physq=physqvel, degfd=[0,1,0], excludesurfaces=[1] ) !symmetry
    call define_essential ( mesh, input_probdef, &
      surface1=5, physq=physqvel, degfd=[0,0,1], excludesurfaces=[1]) !symmetry

    call problem_definition ( input_probdef, mesh, problem )

! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

    call create_sysvector ( problem, sol, sol_n, sol_nm1, sol_p, rhsd )

! fill solution vector with essential boundary conditions

    sol%u = 0

! create subscripts for the meshvelocity

    call create_subscript ( mesh, problem, velx, physqarr=[physqvel], degfd=1, &
                            fillnodes=.true.)
    call create_subscript ( mesh, problem, vely, physqarr=[physqvel], degfd=2, &
                            fillnodes=.true. )
    call create_subscript ( mesh, problem, velz, physqarr=[physqvel], degfd=3, &
                            fillnodes=.true. )

    call create_subscript ( mesh, problem, velx_inlet, physqarr=[physqvel], &
      degfd=1, surfaces=[1], fillnodes=.true. )
    call create_subscript ( mesh, problem, vely_inlet, physqarr=[physqvel], &
      degfd=2, surfaces=[1], fillnodes=.true. )
    call create_subscript ( mesh, problem, velz_inlet, physqarr=[physqvel], &
      degfd=3, surfaces=[1], fillnodes=.true. )

    call create_subscript ( mesh, problem, mvelx_surf1, vec=physqvel, degfd=1, &
       surfaces=[4] )
    call create_subscript ( mesh, problem, mvely_surf1, vec=physqvel, degfd=2, &
       surfaces=[4] )
    call create_subscript ( mesh, problem, mvelz_surf1, vec=physqvel, degfd=3, &
       surfaces=[4] )

    call create_subscript ( mesh, problem, velxsurf1, physqarr=[physqvel], &
      degfd=1, surfaces=[4], fillnodes=.true. )
    call create_subscript ( mesh, problem, velysurf1, physqarr=[physqvel], &
      degfd=2, surfaces=[4], fillnodes=.true. )
    call create_subscript ( mesh, problem, velzsurf1, physqarr=[physqvel], &
      degfd=3, surfaces=[4], fillnodes=.true. )

! create system matrix of gradient/velocity/pressure problem

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )

    call create_sysmatrix_data ( sysmatrix )

    call create_vector ( problem, meshvel, physq=physqvel )

!   define viscoelastic problem, essential conditions and constraints
    call create_input_probdef ( mesh, input_probdefb, nvec=5, nphysq=1 )

    input_probdefb%vec_elementdof(1)%a =   &
        reshape ( [ 1,0,1,0,1,0,0,0,0,1,    &  ! c
                     1,1,1,1,1,1,1,1,1,1,    &  ! scalar for plotting
                     6,6,6,6,6,6,6,6,6,6,    &  ! tensor
                     9,9,9,9,9,9,9,9,9,9,    &  ! tensor
                     ncompv_P1,0,ncompv_P1,0,ncompv_P1,    & ! projection
                             0,0,        0,0,ncompv_p1 ], &
                     [10,5] )

    input_probdefb%physq = [1]
    input_probdefb%probnr = 2


  call define_essential ( mesh, input_probdefb, surface1=1, physq=1 )

  call problem_definition ( input_probdefb, mesh, problemb )

! create system vectors (solution and right-hand side) for conformation and
! initialize vectors with steady stress and a small random perturbation.

  call create ( problemb, solc, solb, solb_n, solb_nm1, rhsb )

! initialize vectors with zero stress
  if ( bvariant == 4 ) then ! log Cholesky
    do m = 1, nmodes
      solb(1,m)%u = 0 ! initial bxx
      solb(2,m)%u = 0 ! initial bxy
      solb(3,m)%u = 0 ! initial bxz
      solb(4,m)%u = 0 ! initial byx
      solb(5,m)%u = 0 ! initial byy
      solb(6,m)%u = 0 ! initial byz
      solb(7,m)%u = 0 ! initial bzx
      solb(8,m)%u = 0 ! initial bzy
      solb(9,m)%u = 0 ! initial bzz
   end do
  else
    do m = 1, nmodes
      solb(1,m)%u = 1 ! initial bxx
      solb(2,m)%u = 0 ! initial bxy
      solb(3,m)%u = 0 ! initial bxz
      solb(4,m)%u = 0 ! initial byx
      solb(5,m)%u = 1 ! initial byy
      solb(6,m)%u = 0 ! initial byz
      solb(7,m)%u = 0 ! initial bzx
      solb(8,m)%u = 0 ! initial bzy
      solb(9,m)%u = 1 ! initial bzz
    end do
  end if

  do i = 1, ncompb
    call create_subscript ( mesh, problemb, bb_inlet(i), surfaces=[1], &
      physqarr=[1], fillnodes=.true. )
  end do

  call create_subscript ( mesh, problemb, bval, physqarr=[1] )
  call create_subscript ( mesh, problemb, cxx, degfd=1, vec=3 )
  call create_subscript ( mesh, problemb, cxy, degfd=2, vec=3 )
  call create_subscript ( mesh, problemb, cxz, degfd=3, vec=3 )
  call create_subscript ( mesh, problemb, cyy, degfd=4, vec=3 )
  call create_subscript ( mesh, problemb, cyz, degfd=5, vec=3 )
  call create_subscript ( mesh, problemb, czz, degfd=6, vec=3 )

! create system matrix for conformation problem

  call create_sysmatrix_structure ( sysmatrixb, mesh, problemb )

  call create_sysmatrix_data ( sysmatrixb )

!   initialize all vectors to zero
    sol_nm1%u = 0.0_dp
    sol_n%u = 0.0_dp
    do i = 1, ncompb
      do m = 1, nmodes
        solb_n(i,m)%u = 0.0_dp
        solb_nm1(i,m)%u = 0.0_dp
      end do
    end do
    meshvel%u = 0._dp

!   problem definition for projected "c=b*b^T"
    call create_input_probdef ( mesh, input_probdefc_projc, nvec=1, &
      nphysq=1 )

      input_probdefc_projc%vec_elementdof(1)%a = &
          reshape ( [ 1,0,1,0,1,0,0,0,0,1 ], &
                       [10,1] )

    input_probdefc_projc%physq = [1]
    input_probdefc_projc%probnr = 3

    call problem_definition ( input_probdefc_projc, mesh, problemc_projc )

    call create ( problemc_projc, solc_projc, rhsc_projc )

    solc_projc(1,1)%u = 1 ! initial cxx
    solc_projc(2,1)%u = 0 ! initial cxy
    solc_projc(3,1)%u = 0 ! initial cxz
    solc_projc(4,1)%u = 1 ! initial cyy
    solc_projc(5,1)%u = 0 ! initial cyz
    solc_projc(6,1)%u = 1 ! initial czz

!   create system matrix for projection problem
    call create_sysmatrix_structure ( sysmatrixc_projc, mesh, &
      problemc_projc)!, symmetric=.true. )
    call create_sysmatrix_data ( sysmatrixc_projc )

! create the structure oldvectors_ve

    call create_oldvectors ( oldvectors_ve, nsysvec=2,  nsysvec2=3, nprob=3, &
      nvec=2 )

! store solution vectors and problem structures

    oldvectors_ve%s(1)%p => sol
    oldvectors_ve%s(2)%p => sol_n
    oldvectors_ve%s2(1)%p => solb
    oldvectors_ve%s2(2)%p => solb_n
    oldvectors_ve%s2(3)%p => solc_projc
    oldvectors_ve%p(1)%p => problem
    oldvectors_ve%p(2)%p => problemb
    oldvectors_ve%p(3)%p => problemc_projc
    oldvectors_ve%v(1)%p => meshvel

  end subroutine define_problems_create_vectors

  subroutine build_solve_vpG

!   copy old solution
    call copy ( sol, sol_n )
    call copy ( sol_n, sol_nm1 )

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, physqrow=[2,3], physqcol=[2,3], &
      coefficients=coefficients )

!   DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_elem, addmatvec=.true., &
      physqrow=[1,2], physqcol=[1,2], coefficients=coefficients )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[1], physqcol=[3], zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[3], physqcol=[1], zeromatvec=.true. )

!   build implicit terms of CE

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve, physqrow=[2], physqcol=[2], &
      addmatvec=.true., coefficients=coefficients )

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve gradient/velocity/pressure problem

!    if ( step == 1 ) then
!       solver_options_u%pivot_order = 1
!       solver_options_u%scaling = 1
!       call renumber_metis_sysmatrix(sysmatrix)
!    end if

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

  end subroutine build_solve_vpG

  subroutine build_solve_ve

!   build (assemble) matrix and vector for conformation problem

    integer :: i, m

    if ( coefficients%i(22) == timeint1 ) then

      call build_system ( mesh, problemb, sysmatrixb, m2sysvector=rhsb, &
        elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    else

      call build_system ( mesh, problemb, sysmatrixb, m2sysvector=rhsb, &
        elemsub=ce_supg_elem_implicit_2nd_order, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    end if


    call check ( sysmatrixb )
!
    call copy ( solb, solb_n )
    call copy ( solb_n, solb_nm1 )

!   solve conformation and keep LU decomposition in the loop over components

!    if ( step == 1 ) then
!       solver_options_b%pivot_order = 1
!       solver_options_b%scaling = 1
!       call renumber_metis_sysmatrix(sysmatrixb)
!    end if

    solver_options_b%real_storage=rs_c
    solver_options_b%integer_storage=is_c

    do m = 1, nmodes
      do i = 1, ncompb
        call add_effect_of_essential_to_rhs ( problemb, sysmatrixb, &
            solb(i,m), rhsb(i,m) )
        call solve_system_ma41 ( sysmatrixb, rhsb(i,m), &
           solb(i,m), lub, solver_options=solver_options_b )

      end do
    end do


    call delete ( lub )  ! remove LU decomposition and rebuild next time step

  end subroutine build_solve_ve

! write the data to .vtk files

  subroutine postprocessing ( mesh, sol, ipost  )

    type(mesh_t), intent(inout) :: mesh
    type(sysvector_t), intent(in), target :: sol
    integer, intent(in), optional :: ipost

    integer :: i, k

    type(oldvectors_t) :: oldvectors_dve
    type(vector_t) :: pressure, eff_shear

    if ( .not. mesh%meshparts) call fill_mesh_parts ( mesh )

    write (filename,'(a,i4.4,a)')'meshtime', ipost, '.vtk'
    call write_mesh_vtk ( mesh, filename=filename )

    call create_oldvectors ( oldvectors_dve, nsysvec=1, nsysvec2=1 )
    oldvectors_dve%s(1)%p => sol

    call create_vector ( problem, pressure, vec=4 )

!   derive the pressure in all nodes
    call derive_vector ( mesh, problem, pressure, &
      elemsub=stokes_pressure, coefficients=coefficients, &
      oldvectors=oldvectors_dve )

    call create_vector ( problem, eff_shear, vec=7 )

    coefficients%i(13) = 8

!   derive the effective shear rate in all nodes
    call derive_vector ( mesh, problem, eff_shear, &
      elemsub=stokes_deriv, coefficients=coefficients, &
      oldvectors=oldvectors_dve )

    if ( present(ipost) ) write(filename,'(a,i4.4,a)') 'flow', ipost, '.vtk'
    call write_scalar_vtk ( mesh, problem, vector=pressure, &
      dataname='pressure',  filename=filename )

    call write_vector_vtk ( mesh, problem, filename=filename, &
      dataname='velocity', sysvector=sol, physq=physqvel, &
      append=.true. )

    call delete(pressure)

    call delete(oldvectors_dve)

!   conformation tensor
!   create a vector for conformation and b-tensor
!   tensor for post processing

    call create_vector ( problemb, ctensor, vec=3 )
    call create_vector ( problemb, btensor, vec=4 )

    call derive_vector ( mesh, problemb, ctensor, &
      elemsub=deriv_conformation_tensor, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    if ( present(ipost) ) write(filename,'(a,i4.4,a)') 'c_', ipost, '.vtk'
    call write_tensor_vtk ( mesh, problemb, filename=filename, &
      dataname='conformation_tensor', vector=ctensor, append=.true. )

!   b tensor

    call derive_vector ( mesh, problemb, btensor, &
      elemsub=deriv_conformation_tensor_std, &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    if ( present(ipost) ) write(filename,'(a,i4.4,a)') 'b_', ipost, '.vtk'
    call write_tensor_vtk ( mesh, problemb, filename=filename, &
      dataname='b-tensor', vector=btensor, symmetric=.false., append=.true. )

    call delete ( ctensor, btensor )

    if ( step == numtimesteps ) then
       call create ( mesh, problem, cxx13, vec=5, degfd=1, curves=[13], &
         fillnodes=.true. )
       call create ( mesh, problem, cxy13, vec=5, degfd=2, curves=[13] )
       call create ( mesh, problem, cxz13, vec=5, degfd=3, curves=[13] )
       call create ( mesh, problem, cyy13, vec=5, degfd=4, curves=[13] )
       call create ( mesh, problem, cyz13, vec=5, degfd=5, curves=[13] )
       call create ( mesh, problem, czz13, vec=5, degfd=6, curves=[13] )

       call create_vector ( problemb, ctensor, vec=3 )

       call derive_vector ( mesh, problemb, ctensor, &
          elemsub=deriv_conformation_tensor,&
          coefficients=coefficients, oldvectors=oldvectors_ve )

       do i = 1, size(cxx13%s)
          k = cxx13%nodes(i)
          write(405,'(15es16.8)') mesh%coor(k,1), mesh%coor(k,2), &
          ctensor%u(cxx13%s(i)), ctensor%u(cxy13%s(i)), ctensor%u(cxz13%s(i)), &
          ctensor%u(cyy13%s(i)), ctensor%u(cyz13%s(i)), ctensor%u(czz13%s(i))
       end do

      call delete ( ctensor )
      call delete ( cxx13, cxy13, cxz13, cyy13, cyz13, czz13 )

      close(unit=405)

     end if

  end subroutine postprocessing

  subroutine delete_old_problems

!   delete all data including all allocated memory

    deallocate ( meshcoor_initial )

    call delete ( problem, problemb )
    call delete ( sysmatrix, sysmatrixb )
    call delete ( input_probdef, input_probdefb )
    call delete ( mesh )
    call delete ( sol, sol_n, sol_nm1 )
    call delete ( solb, solb_n, solb_nm1 )
    call delete ( rhsb )
    call delete ( rhsd )
    call delete ( oldvectors_ve )
    call delete ( meshvel )
    call delete ( problemc_projc )
    call delete ( input_probdefc_projc )
    call delete ( solc_projc, rhsc_projc )
    call delete ( sysmatrixc_projc )

!   deallocate area and aspect ratio arrays
    deallocate ( meshcoor_n, meshcoor_nm1 )

    ! delete all inlet data including all allocated memory

    call delete ( coefficients_sf_adv )
    call delete ( problem_surf1)
    call delete ( sol_surf1 )
    call delete ( sysmatrix_surf1 )
    call delete ( input_probdef_surf1)

    call delete ( mesh_inlet )
    call delete ( problem_inlet )
    call delete ( input_probdef_inlet )
    call delete ( sol_inlet, rhsd_inlet )
    call delete ( sysmatrix_inlet )
    call delete ( oldvectors_ve_inlet )
    call delete ( problemb_inlet )
    call delete ( input_probdefb_inlet )
    call delete ( sysmatrixb_inlet )
    call delete ( solb_inlet, rhsb_inlet )

    call delete ( coefficients_inlet )
    call delete ( velx_in, vely_in, velz_in, cval_in )

  end subroutine delete_old_problems

  subroutine define_inlet_problem

    integer :: i

!   create mesh

    call generate_inlet_mesh

!   problem definition of gradient/velocity/pressure for the channel

    call create_input_probdef ( mesh_inlet, input_probdef_inlet, nvec=6, nphysq=3 )

    input_probdef_inlet%vec_elementdof(1)%a(:,1) = 0
    input_probdef_inlet%vec_elementdof(1)%a(vertices_in,1) = 6  ! gradients
    input_probdef_inlet%vec_elementdof(1)%a(:,2) = 3  ! velocity
    input_probdef_inlet%vec_elementdof(1)%a(:,3) = 0
    input_probdef_inlet%vec_elementdof(1)%a(vertices_in,3) = 1  ! pressure
    input_probdef_inlet%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity
    input_probdef_inlet%vec_elementdof(1)%a(:,5) = 6  ! symmetric tensor
    input_probdef_inlet%vec_elementdof(1)%a(:,6) = 9  ! unsymmetric tensor

    input_probdef_inlet%physq = [1,2,3]
    input_probdef_inlet%probnr = 1

!   outer wall
    call define_essential ( mesh_inlet, input_probdef_inlet, curve1=1, physq=physqvel )
    call define_essential ( mesh_inlet, input_probdef_inlet, curve1=3, physq=physqvel, degfd=[1,0,0] )
    call define_essential ( mesh_inlet, input_probdef_inlet, curve1=2, physq=physqvel, degfd=[0,1,0] )
!   pressure level
    call define_essential ( mesh_inlet, input_probdef_inlet, point=1, physq=physqpress )

!   constraint for flow rate

    call define_constraint ( mesh_inlet, input_probdef_inlet, &
      physq=physqvel, surface1=1, nglobalc=1 )

    call problem_definition ( input_probdef_inlet, mesh_inlet, problem_inlet )

!   problem definition conformation tensor

    call create_input_probdef ( mesh_inlet, input_probdefb_inlet, nvec=4, nphysq=1 )

    input_probdefb_inlet%vec_elementdof(1)%a(:,1) = 0
    input_probdefb_inlet%vec_elementdof(1)%a(vertices_in,1) = 1  ! b
    input_probdefb_inlet%vec_elementdof(1)%a(:,2) = 1         ! scalar for plotting
    input_probdefb_inlet%vec_elementdof(1)%a(:,3) = 6         ! symmetric tensor
    input_probdefb_inlet%vec_elementdof(1)%a(:,4) = 9         ! unsymmetric tensor

    input_probdefb_inlet%physq = [1]
    input_probdefb_inlet%probnr = 2

    call problem_definition ( input_probdefb_inlet, mesh_inlet, problemb_inlet )

!   create system vectors (solution and right-hand side)

    call create_sysvector ( problem_inlet, sol_inlet, solm1_inlet, rhsd_inlet )

    call create_subscript ( mesh_inlet, problem_inlet, velx_in, physqarr=[physqvel], &
        degfd=1, surfaces=[1], fillnodes=.true. )
    call create_subscript ( mesh_inlet, problem_inlet, vely_in, physqarr=[physqvel], &
        degfd=2, surfaces=[1], fillnodes=.true. )
    call create_subscript ( mesh_inlet, problem_inlet, velz_in, physqarr=[physqvel], &
        degfd=3, surfaces=[1], fillnodes=.true. )

!   fill solution vector with essential boundary conditions

    sol_inlet%u = 0
    call fill_sysvector ( mesh_inlet, problem_inlet, sol_inlet, &
      curve1=1, physq=physqvel, degfd=1, value=0._dp )
    call fill_sysvector ( mesh_inlet, problem_inlet, sol_inlet, &
      point=1, physq=physqpress, value=0._dp )

!   create system matrix

    call create_sysmatrix_structure_base ( sysmatrix_inlet, mesh_inlet, problem_inlet )
    call create_sysmatrix_structure_constraint ( sysmatrix_inlet, mesh_inlet, problem_inlet )
    call finalize_sysmatrix_structure ( sysmatrix_inlet )

    call create_sysmatrix_data ( sysmatrix_inlet )

!   create system vectors (solution and right-hand side) for conformation and
!   initialize vectors to zero stress

    call create ( problemb_inlet, solc_inlet, solb_inlet, solbm1_inlet, rhsb_inlet )

    if ( bvariant == 4 ) then ! log Cholesky
      do i = 1, ncompb
        solb_inlet(i,1)%u = 0
      end do
      do i = 1,ncompc
        solc_inlet(i,1)%u = 0
      end do
    else
      solb_inlet(1,1)%u = 1
      solb_inlet(2,1)%u = 0
      solb_inlet(3,1)%u = 0
      solb_inlet(4,1)%u = 0
      solb_inlet(5,1)%u = 1
      solb_inlet(6,1)%u = 0
      solb_inlet(7,1)%u = 0
      solb_inlet(8,1)%u = 0
      solb_inlet(9,1)%u = 1
    end if

    do i = 1, ncompb
      call create_subscript ( mesh_inlet, problemb_inlet, bb_in(i), surfaces=[1], &
        physqarr=[1], fillnodes=.true. )
    end do

    call create_subscript ( mesh_inlet, problemb_inlet, bval_in, physqarr=[1] )

!   create system matrix for conformation problem

    call create_sysmatrix_structure_base ( sysmatrixb_inlet, mesh_inlet, problemb_inlet )
    call create_sysmatrix_structure_constraint ( sysmatrixb_inlet, mesh_inlet, problemb_inlet )
    call finalize_sysmatrix_structure ( sysmatrixb_inlet )

    call create_sysmatrix_data ( sysmatrixb_inlet )

!   problem definition for projected "c=bb^T"

    call create_input_probdef ( mesh_inlet, input_probdefc_proj_inlet, nvec=1, nphysq=1 )

    input_probdefc_proj_inlet%vec_elementdof(1)%a(:,1) = 0
    input_probdefc_proj_inlet%vec_elementdof(1)%a(vertices_in,1) = 1  ! c

    input_probdefc_proj_inlet%physq = [1]
    input_probdefc_proj_inlet%probnr = 3

    call problem_definition ( input_probdefc_proj_inlet, mesh_inlet, problemc_proj_inlet )

    call create ( problemc_proj_inlet, solc_proj_inlet, rhsc_proj_inlet )

!   initialize solc_proj
    solc_proj_inlet(1,1)%u = 1
    solc_proj_inlet(2,1)%u = 0
    solc_proj_inlet(3,1)%u = 0
    solc_proj_inlet(4,1)%u = 1
    solc_proj_inlet(5,1)%u = 0
    solc_proj_inlet(6,1)%u = 1

!   create the structure oldvectors_ve

    call create_oldvectors ( oldvectors_ve_inlet, nsysvec=2, nsysvec2=3, nprob=3 )

!   store solution vectors and problem structures

    oldvectors_ve_inlet%s(1)%p => sol_inlet
    oldvectors_ve_inlet%s(2)%p => solm1_inlet
    oldvectors_ve_inlet%s2(1)%p => solb_inlet
    oldvectors_ve_inlet%s2(2)%p => solbm1_inlet
    oldvectors_ve_inlet%s2(3)%p => solc_proj_inlet
    oldvectors_ve_inlet%p(1)%p => problem_inlet
    oldvectors_ve_inlet%p(2)%p => problemb_inlet
    oldvectors_ve_inlet%p(3)%p => problemc_proj_inlet

!   create and build system matrix for projection problem
!   NOTE matrix remains constant and needs to be build once.

    call create_sysmatrix_structure ( sysmatrixc_proj_inlet, mesh_inlet, problemc_proj_inlet, &
      symmetric=.true. )
    call create_sysmatrix_data ( sysmatrixc_proj_inlet )

    call build_system ( mesh_inlet, problemc_proj_inlet, sysmatrixc_proj_inlet, &
      m2sysvector=rhsc_proj_inlet, elemsub=exps_projection_elem, &
      oldvectors=oldvectors_ve_inlet, coefficients=coefficients_inlet, &
      buildvector=.false. )

    call check ( sysmatrixc_proj_inlet )

  end subroutine define_inlet_problem

  subroutine solve_inlet_problem

    integer :: icomp

    if ( step == 2 ) then
!     change time integration scheme at the second time step
      coefficients_inlet%i(22) = timeint2
    end if

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG_inlet

    if ( cproj ) &
             call solve_projection ( solc_proj_inlet, rhsc_proj_inlet, &
                   c_projection_elem )

!   build implicit terms of CE with rhs in momentum balance

    call build_system ( mesh_inlet, problem_inlet, sysmatrix_inlet, rhsd_inlet, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve_inlet, physqrow=[2], physqcol=[2], &
      addmatvec=.true., coefficients=coefficients_inlet )

    call add_effect_of_essential_to_rhs ( problem_inlet, sysmatrix_inlet, sol_inlet, rhsd_inlet )

!   solve gradient/velocity/pressure problem

    solver_options_u_in%real_storage=rs_gup_in
    solver_options_u_in%integer_storage=is_gup_in

    call solve_system_ma41 ( sysmatrix_inlet, rhsd_inlet, sol_inlet, &
      solver_options=solver_options_u_in  )

    call copy ( sol_inlet, solm1_inlet )

!   build (assemble) matrix and vector for conformation problem

    if ( coefficients_inlet%i(22) == timeint1 ) then

      call build_system ( mesh_inlet, problemb_inlet, sysmatrixb_inlet, m2sysvector=rhsb_inlet, &
        elemsub=ce_supg_elem, oldvectors=oldvectors_ve_inlet, &
        coefficients=coefficients_inlet )

    else

      call build_system ( mesh_inlet, problemb_inlet, sysmatrixb_inlet, m2sysvector=rhsb_inlet, &
        elemsub=ce_supg_elem_implicit_2nd_order, oldvectors=oldvectors_ve_inlet, &
        coefficients=coefficients_inlet )

    end if


    call check ( sysmatrixb_inlet )

    call copy ( solb_inlet, solbm1_inlet )

!   solve conformation and keep LU decomposition in the loop over components

    solver_options_b_in%real_storage=rs_c_in
    solver_options_b_in%integer_storage=is_c_in

    do icomp = 1, ncompb
      call solve_system_ma41 ( sysmatrixb_inlet, rhsb_inlet(icomp,1), solb_inlet(icomp,1), lub_inlet, &
        solver_options=solver_options_b_in  )
    end do

    call delete ( lub_inlet )  ! remove LU decomposition and rebuild next time step

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

    if ( bvariant == 4 ) call calculate_c_inlet

  end subroutine solve_inlet_problem

  subroutine generate_inlet_mesh

    integer :: i, node, node1, node2, nodenr2, curve, curve2

    call mesh_skeleton ( mesh_inlet, mesh%surfaces(1)%nnodes, &
      mesh%surfaces(1)%nelem, mesh%surfaces(1)%element%elshape, ndim=2 )

     mesh_inlet%coor = mesh%coor(mesh%surfaces(1)%nodes,[2,3])
     mesh_inlet%topology(1)%a = mesh%surfaces(1)%topology(:,:,1)

     allocate(bd(mesh_inlet%ndim,2))

     do i=1,mesh_inlet%ndim
       bd(i,1)=minval(mesh_inlet%coor(:,i))-1.e-10_dp
       bd(i,2)=maxval(mesh_inlet%coor(:,i))+1.e-10_dp
     end do

     call add_to_mesh ( mesh_inlet, blocks=[(5,i=1,mesh_inlet%ndim)], blocksdomain=bd )

     deallocate(bd)

!   add curves to the 2D mesh

    deallocate(mesh_inlet%curves)
    allocate(mesh_inlet%curves(3))
    mesh_inlet%ncurves = 3

    do curve = 1, 3

      if (curve == 1) curve2=1
      if (curve == 2) curve2=2
      if (curve == 3) curve2=3

      mesh_inlet%curves(curve)%nnodes = mesh%curves(curve2)%nnodes
      mesh_inlet%curves(curve)%nelem = mesh%curves(curve2)%nelem
      mesh_inlet%curves(curve)%elnumnod = mesh%curves(curve2)%elnumnod
      allocate ( mesh_inlet%curves(curve)%nodes(mesh_inlet%curves(curve)%nnodes))
      allocate ( mesh_inlet%curves(curve)%topology( &
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

          if ( abs(mesh%coor(nodenr2,2)-mesh_inlet%coor(node1,1))<1e-12_dp .and. &
               abs(mesh%coor(nodenr2,3)-mesh_inlet%coor(node1,2))<1e-12_dp ) then

!           matching node
            mesh_inlet%curves(curve)%nodes(node2) = node1

            exit ! leave loop

          end if

        end do

      end do
      !   update the global node numbers in the toplogy
      do node = 1, mesh_inlet%curves(curve)%nnodes
        do i = 1, mesh_inlet%curves(curve)%nelem
          mesh_inlet%curves(curve)%topology(:,i,2) = &
            mesh_inlet%curves(curve)%nodes(mesh_inlet%curves(curve)%topology(:,i,1))
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
    call add_to_mesh ( mesh_inlet, point=[0._dp,Ro])

    call find_bounds_blocks ( mesh_inlet )

    call fill_mesh_parts ( mesh_inlet )

    call write_mesh_vtk ( mesh_inlet, 'mesh_inlet.vtk' )

    do i = 1, mesh_inlet%ncurves
      write(filename,'(a,i4.4,a)') 'curvein_',i,'.vtk'
      call write_geometry_vtk ( mesh_inlet, curve=i, filename=filename )
    end do

    call write_mesh_vtk ( mesh_inlet, filename='mesh_inlet.vtk' )

  end subroutine generate_inlet_mesh

! write the parameters in gmsh format

  subroutine write_gmsh_inlet (ox, dx_box )

    real(dp), intent(in) ::ox, dx_box

    open ( unit=25, file='mesh_inlet.geo' )

    write ( 25, '(1X,A,F18.14,A)' ) 'ox = ', ox, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_box = ', dx_box, ';'

    write ( 25, '(1X,A,F18.14,A)' ) 'y[1] = ', Ro*cos(theta), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'z[1] = ', Ro*sin(theta), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'y[2] = ', Ro*cos(3*theta), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'z[2] = ', Ro*sin(3*theta), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'y[3] = ', Ro*cos(5*theta), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'z[3] = ', Ro*sin(5*theta), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'y[4] = ', Ro*cos(7*theta), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'z[4] = ', Ro*sin(7*theta), ';'
    write ( 25, '(/1x,a)' ) 'Include "mesh_inlet.igo";'

    close ( 25 )

  end subroutine write_gmsh_inlet

  subroutine build_vpG_inlet

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh_inlet, problem_inlet, sysmatrix_inlet, rhsd_inlet, &
      elemsub=stokes_elem, physqrow=[2,3], physqcol=[2,3], &
      coefficients=coefficients_inlet )

!   DEVSS-G
    call build_system ( mesh_inlet, problem_inlet, sysmatrix_inlet, rhsd_inlet, &
      elemsub=devssg_elem, addmatvec=.true., &
      physqrow=[1,2], physqcol=[1,2], coefficients=coefficients_inlet )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh_inlet, problem_inlet, sysmatrix_inlet, rhsd_inlet, addmatvec=.true., &
      buildvector=.false., physqrow=[1], physqcol=[3], zeromatvec=.true. )
    call build_system ( mesh_inlet, problem_inlet, sysmatrix_inlet, rhsd_inlet, addmatvec=.true., &
      buildvector=.false., physqrow=[3], physqcol=[1], zeromatvec=.true. )

!   imposed flow rate

    call build_system_constraint ( mesh_inlet, problem_inlet, sysmatrix_inlet, rhsd_inlet, &
      constraint1=1, elemsub=stokes_constr_flowr_surface, addmatvec=.true., &
      coefficients=coefficients_inlet )

  end subroutine build_vpG_inlet

  subroutine postprocessing_inlet ( mesh_in, sol_in, ipost  )

    type(vector_t) :: velocity_in, pressure
    type(oldvectors_t) :: oldvectors_in

    type(mesh_t), intent(inout) :: mesh_in
    type(sysvector_t), intent(in), target :: sol_in
    integer, intent(in), optional :: ipost

!   create the structure oldvectors

    call create_oldvectors ( oldvectors_in, nsysvec=1, nsysvec2=1 )

    call create_vector ( problem_inlet, velocity_in, physq=physqvel )
    call create_vector ( problem_inlet, pressure, vec=4 )

    call extract_physvector ( mesh_in, problem_inlet, sol_inlet, velocity_in )

    oldvectors_in%s(1)%p => sol_in

    call derive_vector ( mesh_in, problem_inlet, pressure, elemsub=stokes_pressure, &
      coefficients=coefficients_inlet, oldvectors=oldvectors_in )

    coefficients_inlet%i(13)=11

!   write to a vtk file for further post-processing

    if ( present(ipost) ) write(filename,'(a,i4.4,a)') 'inletproblem', ipost, '.vtk'
    call write_scalar_vtk ( mesh_in, problem_inlet, vector=pressure, dataname='pressure', &
      filename=filename )

    if ( present(ipost) ) write(filename,'(a,i4.4,a)') 'inletproblem', ipost, '.vtk'
    call write_vector_vtk ( mesh_in, problem_inlet, filename=filename, &
      dataname='velocity_vector', vector=velocity_in, append=.true., &
      assume3D=.true. )

    oldvectors_in%s2(1)%p => solb_inlet
    call create_vector ( problemb_inlet, ctensor_in, vec=3 )
    call create_vector ( problemb_inlet, btensor_in, vec=4 )

    call derive_vector ( mesh_in, problemb_inlet, ctensor_in, &
      elemsub=deriv_conformation_tensor,&
      coefficients=coefficients_inlet, oldvectors=oldvectors_ve_inlet )

    if ( present(ipost) ) write(filename,'(a,i4.4,a)') 'inletproblem', ipost, '.vtk'
    call write_tensor_vtk ( mesh_in, problemb_inlet, filename=filename, &
      dataname='c', vector=ctensor_in, append=.true., &
      assume3D=.true. )

!   b tensor

    call derive_vector ( mesh_inlet, problemb_inlet, btensor_in, &
      elemsub=deriv_conformation_tensor_std, &
      coefficients=coefficients_inlet, oldvectors=oldvectors_ve_inlet )

    if ( present(ipost) ) write(filename,'(a,i4.4,a)') 'bin_', ipost, '.vtk'
    call write_tensor_vtk ( mesh_inlet, problemb_inlet, filename=filename, &
      dataname='b-tensor', vector=btensor_in, symmetric=.false., append=.true., &
      assume3D=.true. )

    call delete(oldvectors_in)
    call delete(velocity_in)
    call delete(pressure)
    call delete(ctensor_in, btensor_in)

  end subroutine postprocessing_inlet

  subroutine define_surface_problems

    integer :: i

    call create_surface_mesh

!   define some arrays for the ALE mesh position at old times
    allocate ( meshcoorsurf1_n(mesh_surf1%nnodes,mesh_surf1%ndim), &
        meshcoorsurf1_nm1(mesh_surf1%nnodes,mesh_surf1%ndim) )

    allocate(initial_hsurf1(mesh_surf1%nnodes))

    initial_hsurf1 = sqrt( mesh%coor(mesh%surfaces(4)%nodes,2)**2 + &
                      mesh%coor(mesh%surfaces(4)%nodes,3)**2 )

    allocate(surf1hat(mesh_surf1%nnodes,2),surf1hatn(mesh_surf1%nnodes,2))

    surf1hat(:,2) = mesh%coor(mesh%surfaces(4)%nodes(:),2)
    surf1hat(:,1) = mesh%coor(mesh%surfaces(4)%nodes(:),3)

    surf1hatn = surf1hat

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

    call define_essential ( mesh_surf1, input_probdef_surf1, curve1=1, physq=1 )

!   define problem
    call problem_definition ( input_probdef_surf1, mesh_surf1, problem_surf1 )

    call create_subscript ( mesh_surf1, problem_surf1, subs_surf1, physqarr=[1], fillnodes=.true. )
    call create_subscript ( mesh_surf1, problem_surf1, pos1_surf1, physqarr=[1], curves=[4])
    call create_subscript ( mesh_surf1, problem_surf1, pos2_surf1, physqarr=[1], curves=[3])
    call create_subscript ( mesh_surf1, problem_surf1, subs_curve1, physqarr=[1], curves=[1])

    call create ( problem_surf1, sol_surf1, rhsd_surf1 )
    call create ( problem_surf1, sol_surf1_n, sol_surf1_nm1 )

!   fill solution vector with essential boundary conditions

    do i = 1, mesh_surf1%nnodes
      sol_surf1%u(subs_surf1%s(i)) = initial_hsurf1(i)
    end do
    sol_surf1_n%u = sol_surf1%u
    sol_surf1_nm1%u = sol_surf1%u

!   create system matrix

    call create_sysmatrix_structure ( sysmatrix_surf1, mesh_surf1, &
      problem_surf1 )
    call create_sysmatrix_data ( sysmatrix_surf1 )

!   create vectors

    call create_vector ( problem_surf1, velocity_surf1, vec=2 )
    velocity_surf1%u = 0 ! initialize
    call create_vector ( problem_surf1, meshvel_surf1, vec=3 )
    meshvel_surf1%u = 0

    call create_subscript ( mesh_surf1, problem_surf1, velx_surf1, vec=2, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_surf1, problem_surf1, vely_surf1, vec=2, &
        degfd=2, fillnodes=.true. )
    call create_subscript ( mesh_surf1, problem_surf1, velz_surf1, vec=2, &
        degfd=3, fillnodes=.true. )

    call create_subscript ( mesh_surf1, problem_surf1, meshvelx_surf1, vec=3, &
        degfd=1, fillnodes=.true. )
    call create_subscript ( mesh_surf1, problem_surf1, meshvely_surf1, vec=3, &
        degfd=2, fillnodes=.true. )

!   create subscript for the height values

    call create ( mesh_surf1, problem_surf1, hgt )
    call create ( mesh_surf1, problem_surf1, hgt_end, curves=[2] )

    call create_oldvectors ( oldvectors_surf1, nsysvec=2, nvec=1 )

    oldvectors_surf1%s(1)%p => sol_surf1_n    ! corrector at n
    oldvectors_surf1%s(2)%p => sol_surf1_nm1  ! corrector at nm1

    oldvectors_surf1%v(1)%p => velocity_surf1 ! advection velocity at np1

  end subroutine define_surface_problems

  subroutine create_surface_mesh

    integer :: i, node, nodenr2, curve, node1, node2, curve2
!   create mesh for surface advection

    call mesh_skeleton ( mesh_surf1, mesh%surfaces(4)%nnodes, &
      mesh%surfaces(4)%nelem, mesh%surfaces(4)%element%elshape, ndim=2 )

     mesh_surf1%coor(:,1) = mesh%coor(mesh%surfaces(4)%nodes,1)
!    set coordinates (angles)

     do i = 1, mesh%surfaces(4)%nnodes
        phi = find_angle2 ( mesh%coor(mesh%surfaces(4)%nodes(i),2:3) )
        mesh_surf1%coor(i,2) = phi
     end do

     mesh_surf1%topology(1)%a = mesh%surfaces(4)%topology(:,:,1)

     allocate(bd(mesh_surf1%ndim,2))

     do i=1,mesh_surf1%ndim
       bd(i,1)=minval(mesh_surf1%coor(:,i))-1.e-10_dp
       bd(i,2)=maxval(mesh_surf1%coor(:,i))+1.e-10_dp
     end do

     call add_to_mesh ( mesh_surf1, blocks=[(3,i=1,mesh_surf1%ndim)], blocksdomain=bd )

     deallocate(bd)

!   add curves to the 2D mesh

    deallocate(mesh_surf1%curves)
    allocate(mesh_surf1%curves(4))
    mesh_surf1%ncurves = 4

    do curve = 1, 4

       if (curve == 1) curve2=4
       if (curve == 2) curve2=12
       if (curve == 3) curve2=10
       if (curve == 4) curve2=5

      mesh_surf1%curves(curve)%nnodes = mesh%curves(curve2)%nnodes
      mesh_surf1%curves(curve)%nelem = mesh%curves(curve2)%nelem
      mesh_surf1%curves(curve)%elnumnod = mesh%curves(curve2)%elnumnod
      allocate ( mesh_surf1%curves(curve)%nodes(mesh_surf1%curves(curve)%nnodes))
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

          phi = find_angle2 ( mesh%coor(nodenr2,2:3) )
          if ( abs(phi - mesh_surf1%coor(node1,2))<1e-12_dp .and. &
            abs(mesh%coor(nodenr2,1)-mesh_surf1%coor(node1,1))<1e-12_dp ) then

  !          matching node
             mesh_surf1%curves(curve)%nodes(node2) = node1

             exit ! leave loop

          end if

        end do

      end do
!     update the global node numbers in the toplogy
      do node = 1, mesh_surf1%curves(curve)%nnodes
        do i = 1, mesh_surf1%curves(curve)%nelem
          mesh_surf1%curves(curve)%topology(:,i,2) = &
            mesh_surf1%curves(curve)%nodes(mesh_surf1%curves(curve)%topology(:,i,1))
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

  end subroutine create_surface_mesh

! solve convection equation for the surface height (corrector)

  subroutine solve_4surfaces_height_corrector

    use postprocessing_m

    integer :: i, k

    real(dp) :: max_height, end_height
    type(solver_options_ma41_t) :: solver_options_h

    velocity_surf1%u(velx_surf1%s) = sol%u(velxsurf1%s)
    do i = 1, mesh_surf1%nnodes

      k = subs_surf1%nodes(i)
      velocity_surf1%u(vely_surf1%s(i)) = sol%u(velzsurf1%s(i)) * & ! u_theta
        cos(mesh_surf1%coor(k,2)) - sol%u(velysurf1%s(i)) * &
        sin(mesh_surf1%coor(k,2))
      velocity_surf1%u(velz_surf1%s(i)) = sol%u(velzsurf1%s(i)) * & ! u_r
        sin(mesh_surf1%coor(k,2)) + sol%u(velysurf1%s(i)) * &
        cos(mesh_surf1%coor(k,2))
    end do

    call build_system ( mesh_surf1, problem_surf1, sysmatrix_surf1, &
      rhsd_surf1, elemsub=surface_advection_elem, &
      oldvectors=oldvectors_surf1, coefficients=coefficients_sf_adv )

    call check_filled_sysmatrix ( sysmatrix_surf1 )

    call add_effect_of_essential_to_rhs ( problem_surf1, sysmatrix_surf1, &
      sol_surf1, rhsd_surf1 )

!    if (step == 1 ) then
!       solver_options_h%pivot_order = 1
!       solver_options_h%scaling = 1
!       call renumber_metis_sysmatrix(sysmatrix_surf1)
!    end if

!   MA41 solver storage

    solver_options_h%integer_storage = 2.0
    solver_options_h%real_storage    = 2.0

!   solve system

    call solve_system_ma41 ( sysmatrix_surf1, rhsd_surf1, sol_surf1, &
                             solver_options=solver_options_h )

    call copy ( sol_surf1_n, sol_surf1_nm1 )
    call copy ( sol_surf1, sol_surf1_n )

    do i = 1,size(subs_surf1%s)

      k = subs_surf1%nodes(i)
      mesh%coor(mesh%surfaces(4)%nodes(i),2) = &
                sol_surf1%u(subs_surf1%s(i))*cos(mesh_surf1%coor(k,2))
      mesh%coor(mesh%surfaces(4)%nodes(i),3) = &
                sol_surf1%u(subs_surf1%s(i))*sin(mesh_surf1%coor(k,2))

    end do

!   compute the maximum, minimum and end radius

    max_height = maxval ( sol_surf1%u )
    end_height = sol_surf1%u(hgt_end%s(1))

!   write fiber radius

    write(11,'(i6,4es16.8)') step0+step, time0+step*deltat, &
                             max_height, end_height

  end subroutine solve_4surfaces_height_corrector

! find the angle (0 <= angle < 2*pi) of the position vector with positive y-axis
! NOTE: points on the positive y-axis (to precision eps) have angle zero

  function find_angle2 ( x )

    real(dp), parameter :: eps=1e-14_dp

    real(dp), intent(in) :: x(:)
    real(dp) :: find_angle2

    if ( all(abs(x) < eps) ) then
      find_angle2 = 0._dp
    end if

    if ( abs(x(2)) < eps ) then

      find_angle2 = 0._dp

    else

      find_angle2 = atan2 ( x(1), x(2) )

      if ( x(1) < 0._dp ) then ! bottom half-plane
        find_angle2 = pi/2 + abs(find_angle2)
      else
        if ( x(2) > 0._dp ) then ! right half-plane
          find_angle2 = pi/2 - abs(find_angle2)
        else
          find_angle2 = 3*pi/2+(pi-abs(find_angle2))
        end if
      end if

     end if

   end function find_angle2

  subroutine solve_b_projection

    type(solver_options_ma41_t) :: solver_options_proj
    type(oldvectors_t) :: oldvectors_proj
    type(lu_ma41_t) :: lu_projc

    integer :: i, m

    call create_oldvectors ( oldvectors_proj, nprob=2, nsysvec2=1 )
    oldvectors_proj%s2(1)%p => solb_n
    oldvectors_proj%p(2)%p => problemb

!   build system matrix and vector for projection problem

    call build_system ( mesh, problemc_projc, sysmatrixc_projc, &
      m2sysvector=rhsc_projc, elemsub=c_projection_elem, &
      oldvectors=oldvectors_proj, coefficients=coefficients )

    call check ( sysmatrixc_projc )

!   MA41 solver storage

    solver_options_proj%integer_storage = 2.3
    solver_options_proj%real_storage    = 2.3

!   LU decomposition is done in the first call only

    do m = 1, nmodes
      do i = 1, ncompc

        call add_effect_of_essential_to_rhs ( problemc_projc, sysmatrixc_projc, &
           solc_projc(i,m), rhsc_projc(i,m) )

        call solve_system_ma41 ( sysmatrixc_projc, rhsc_projc(i,m), &
           solc_projc(i,m), lu_projc, solver_options=solver_options_proj )

      end do
    end do

   call delete ( lu_projc )

  end subroutine solve_b_projection

  subroutine solve_projection ( solc_proj_inlet, rhsc_proj_inlet, elemsub_proj )

    type(sysvector_t), dimension(:,:), intent(inout) :: solc_proj_inlet

    type(sysvector_t), dimension(:,:), intent(inout) :: rhsc_proj_inlet

    interface
      subroutine elemsub_proj ( mesh, problem, elgrp, elem, matrix, vector, &
        first, last, coefficients, oldvectors, elemmat, elemvec )
        use kind_defs_m
        use mesh_m, only: mesh_t
        use problem_defs_m, only: problem_t
        use element_defs_m, only: coefficients_t, oldvectors_t
        implicit none
        type(mesh_t), intent(in) :: mesh
        type(problem_t), intent(in) :: problem
        integer, intent(in) :: elgrp, elem
        logical, intent(in) :: matrix, vector, first, last
        type(coefficients_t), intent(in) :: coefficients
        type(oldvectors_t), intent(in) :: oldvectors
        real(dp), intent(out), dimension(:,:) :: elemmat
        real(dp), intent(out), dimension(:) :: elemvec
      end subroutine elemsub_proj
    end interface

    type(solver_options_ma57_t) :: solver_options_ma57

    integer :: i, m

!   build vector only (matrix is constant)

    call build_system ( mesh_inlet, problemc_proj_inlet, sysmatrixc_proj_inlet, &
      m2sysvector=rhsc_proj_inlet, elemsub=elemsub_proj, &
      oldvectors=oldvectors_ve_inlet, coefficients=coefficients_inlet, &
      buildmatrix=.false. )

    ! MA57 solver storage
    solver_options_ma57%integer_storage = 1.3
    solver_options_ma57%real_storage    = 1.3

!   LU decomposition is done in the first call only

    do m = 1, nmodes
      do i = 1, ncompb-3

        call add_effect_of_essential_to_rhs ( problemc_proj_inlet, sysmatrixc_proj_inlet, &
           solc_proj_inlet(i,m), rhsc_proj_inlet(i,m) )

        call solve_system_ma57 ( sysmatrixc_proj_inlet, rhsc_proj_inlet(i,m), &
           solc_proj_inlet(i,m), luc_proj_inlet, solver_options=solver_options_ma57 )

      end do
    end do

    call delete ( luc_proj_inlet )

  end subroutine solve_projection

! calculate c from b

  subroutine calculate_c_inlet

    real(dp) :: bb_inlet(size(bval_in%s),ncompb), cc_inlet(size(bval_in%s),ncompc)
    integer :: i

!   obtain the solution of b at current time step

    do i = 1,ncompb
      bb_inlet(:,i) = solb_inlet(i,1)%u(bval_in%s)
    end do

!   calculate c from b

    call conformation_3D_b (bb_inlet, cc_inlet, bvariant)

!   store c in solution vectors

    do i = 1,ncompc
      solc_inlet(i,1)%u(bval_in%s) = cc_inlet(:,i)
    end do

  end subroutine calculate_c_inlet

end program extrudate_swell3d_ba
