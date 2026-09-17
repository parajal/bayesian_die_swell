! Viscoelastic fluid flow around a stationary cylinder
! DEVSS-G/SUPG
! Implicit bilinear terms of CE in momentum balance.
! second-order Gear time integration
! XFEM approach

program cylinder1

  use tfem_m
  use viscoelastic_elements_m
  use hsl_ma41_m
  use io_utils_m
  use usergauss_xfem_m
  use figplot_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    gintpl = 4,         & ! Q1 gradients
    cintpl = 4,         & ! Q1 conformation
    physqgrad = 1,      & ! physical quantity nr of the gradients
    physqvel = 2,       & ! physical quantity nr of the velocities
    physqpress = 3,     & ! physical quantity nr of the pressures
    numsplit = 3,       & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit
    numsplitmin = 2,    & ! relative size of all subelements will be at least
                          ! as small as 1/2^numsplitmin
    numsplit_plot = 2,  & ! relative size of smallest elements near interface
                          ! is 1/2^numsplit for plotting
    gausse = 3,         & ! Gauss integration for the eltree subelements
    gausssub_small = 8, & ! Gauss integration for the submesh subelements for
                          ! small integration volumes
    gausssub = 4,       & ! Gauss integration for the submesh subelements
    gauss_ie = 3,       & ! Gauss integration for the interface elements
    gauss = 3,          & ! 3x3 integration of quads
    ncompc = 3,         & ! number of conformation tensor components
    nmodes = 1,         & ! number of modes
    startm = 501,       & ! start of material model data
    model = 2             ! UCM/Oldroyd-B


! definitions

  type(mesh_t), target :: mesh
  type(mesh_t) :: mesh_plot
  type(input_probdef_t) :: input_probdef, input_probdefc, input_probdef_plot
  type(problem_t), target :: problem, problemc, problem_plot
  type(sysmatrix_t) :: sysmatrix, sysmatrixc
  type(sysvector_t), target :: sol, solm1
  type(sysvector_t) :: rhsd
  type(oldvectors_t) :: oldvectors_ve
  type(coefficients_t) :: coefficients
  type(sysvector_t), dimension(ncompc,nmodes), target :: solc, solcm1
  type(sysvector_t), dimension(ncompc,nmodes) :: rhsc
  type(lu_ma41_t) :: luc
  type(solver_options_ma41_t) :: solver_options_u, solver_options_c
  type(plot_options_t) :: plot_options

  integer, dimension(:), allocatable :: nodes
  real(dp), dimension(:,:), allocatable :: coor

! array of pointers to an eltree

  type(eltree_p), target, dimension(:,:), allocatable :: eltree_array


! variables

  integer :: &
    timeint1 = 1,        & ! (first-order) Euler time integration (first step)
    timeint2 = 7,        & ! (second-order) semi-implicit Gear time integration
    numtimesteps = 300,  & ! number of time steps
    logc = 1               ! standard scheme or log transformation

  real(dp) :: &
    eta_s = 0.59_dp,   & ! solvent viscosity
    eta_p = 0.41_dp,   & ! polymer viscosity
    lambda = 0.3_dp,   & ! relaxation time
    deltat = 1.e-2_dp, & ! time step
    beta = 1.0_dp,     & ! upwinding parameter in the SUPG method
    rs_gup = 1.5_dp,   & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 1.5_dp,   & ! integer_storage gradient-velocity-pressure LU (HSL)
    rs_c  = 1.5_dp,    & ! real_storage for the conformation LU (HSL)
    is_c  = 1.5_dp       ! integer_storage for conformation LU (HSL)

  real(dp), parameter :: &
    split_threshold = 1.e-12_dp, & ! levelset value for being "close" enough
                                   ! to the interface for tree splitting
    epsvol = 1e-6_dp, & ! elements with integration area smaller are removed
                        ! from the eltree_array. This means that:
                        ! 1) these elements are "fully inside" and no
                        !    "virtual/extended degrees" are generated.
                        ! 2) the interface integration is ignored and therefore
                        !    the interface has a small "hole".
                        ! Note, that epsvol is relative to the area of an
                        ! element, so epsvol=1 is a full element.
    epsvol_small = 1e-2_dp ! elements with integration area smaller
                           ! are integrated with integration rule
                           ! gausssub_small, otherwise with gaussub


! (negative) drag values
  real(dp), dimension(2) :: int_pressure,            & ! -pI
                            int_viscous_stress,      & ! 2*eta*D
                            int_viscoelastic_stress, & ! tau_p
                            int_total_stress           ! -pI + 2*eta*D + tau_p

  integer :: icomp, step, i, nnodes, nbx, nby
  real(dp) :: alpha, G
  real(dp) :: kappa ! embedded Dirichlet viscosity parameter

  real(dp) :: flowrate, h
  real(dp) :: radius = 1._dp, U_avg = 1._dp, center(2) = 0._dp


! namelist for input of variables; read from standard input

  namelist /comppar/ numtimesteps, logc, eta_s, eta_p, lambda, deltat, &
    beta, U_avg, radius, rs_gup, is_gup, rs_c, is_c

  read ( unit=*, nml=comppar )


! set some parameters

  alpha = eta_p  ! DEVSS parameter
  G = eta_p / lambda  ! modulus

! stokes usergauss

  intrule_e = gausse  ! Gauss integration for the eltree subelements
  intrule_sub = gausssub  ! Gauss integration for the eltree subelements
  intrule_sub_small = gausssub_small  ! Gauss integration for the eltree
                                      ! subelements (small integration areas)
  vol_small = epsvol_small  ! elements with integration areas smaller
                            ! are integrated with integration rule
                            ! gausssub_small, otherwise with gaussub


! radius and center position of the cylinder
  rpl = radius ! radius of the circular object
  cpl = center ! initial position of the center of the object


! read mesh

  call read_mesh ( mesh, filename='mesh_channel_3p.out' )

  nbx = nint(sqrt(real(mesh%curves(12)%nnodes)*0.5))
  nby = nint(sqrt(real(mesh%curves(9)%nnodes)))

  call add_to_mesh ( mesh, blocks=[nbx,nby] )

  lmesh => mesh  ! supply mesh to mapcoor for mapping coordinates

! create object for the cylinder boundary (only for plotting)

  nnodes = 50

  allocate ( coor(nnodes,2) )

  call objectscoor ( 1, coor )

  call add_to_mesh ( mesh, object='coordinates', coor=coor )

  deallocate ( coor )

! allocate array of pointers to an eltree (eltree for all elements, one group)

  allocate ( eltree_array(1,mesh%grpnumel(1)), vol(mesh%grpnumel(1)) )

! eltree for elements crossing the interface

  allocate ( nodes(mesh%nnodes), d(mesh%nnodes) )

  d = levelset ( mesh%coor ) ! set levelset for all nodes

  call create_eltree_in_elements ( numsplit )

! elementset for elements crossing the interface

  call create_elementset_from_eltree

! nodeset for Dirichlet on internal nodes

  call find_internal_zero_nodes

  call add_to_mesh ( mesh, nodeset='nodes', nodes=nodes(1:nnodes) )

! object for plotting the nodeset

  call add_to_mesh ( mesh, object='coordinates', &
    coor=mesh%coor(nodes(1:nnodes),:) )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves_xf.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh_xf.fig' )
  plot_options%objectpointcolor=4
  plot_options%objectpointsize=0.4
  call plot_objects ( plot_options, mesh, 'mesh_xf.fig', append=.true., &
    object1=1 )
  plot_options%objectpointcolor=5
  call plot_objects ( plot_options, mesh, 'mesh_xf.fig', append=.true.,  &
    object1=2 )


! set flowrate

  h = 2* mesh%coor(mesh%points(3),2) ! height is given by 2*y-coordinate of P3
  flowrate = h * U_avg


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=500+2*nmodes )

  coefficients%i = &
    [ uintpl,   pintpl,     0,     0,         gintpl, &
      physqvel, physqpress, 0,     physqgrad, gauss,  &
      gauss,    cintpl,     0,     0,         0,      &
      0,        0,          model, nmodes,    startm, &
      logc,     timeint1,   ( 0, i = 23, 150 )  &
    ]

  coefficients%r = &
    [ eta_s, 0._dp,   0._dp,   alpha, 0._dp, &
      flowrate, 0._dp,  deltat,   beta,  0._dp, &
      ( 0._dp, i = 11, 500 ), &
      G,     lambda &
    ]

  coefficients%i(36) = -1 ! normal to fluid points inwards to the cylinder
  coefficients%i(37) = 2 ! Gauss points defined by user subroutines
  coefficients%i(41) = gauss_ie ! Gauss points for interface elements
  coefficients%i(42) = -1  ! -1: impose zero velocity, 0: impose non-zero
                           ! velocity U

  coefficients%i(29) = 2       ! htype
  coefficients%i(31) = 4       ! Uscaling
  coefficients%r(10) = 1.0_dp  ! Uglobal

  kappa = eta_s + eta_p
  coefficients%r(11) = kappa

  coefficients%set_ninti_user => set_ninti_user
  coefficients%set_Gauss_integration_user => set_Gauss_integration_user


! problem definition of gradient/velocity/pressure

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=3 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [  4,0,4,0,4,0,4,0,0,   &  ! G
                   2,2,2,2,2,2,2,2,2,   &  ! velocity
                   1,0,1,0,1,0,1,0,0,   &  ! pressure
                   1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                   [9,4] )

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

! Dirichlet boundary conditions

! walls
  call define_essential ( mesh, input_probdef, curve1=12, physq=physqvel )
! outside of the fluid
  call define_essential ( mesh, input_probdef, nodeset1=1 )
! pressure
  call define_essential ( mesh, input_probdef, point=5, physq=physqpress )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=9, nglobalc=1 )

! constraint for periodic b.c

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, curve1=9, curve2=11, discretization='collocation', &
    exclude=3 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors for gradient/velocity/pressure
! (solution and right-hand side)

  call create_sysvector ( problem, sol, solm1, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0

! create system matrix of gradient/velocity/pressure problem

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )


! problem definition conformation tensor

  call create_input_probdef ( mesh, input_probdefc, nvec=2, nphysq=1 )

  input_probdefc%vec_elementdof(1)%a =  &
      reshape ( [ 1,0,1,0,1,0,1,0,0,    &  ! c
                  1,1,1,1,1,1,1,1,1 ],  &  ! scalar for plotting
                  [9,2] )

  input_probdefc%physq = [1]
  input_probdefc%probnr = 3

  call define_essential ( mesh, input_probdefc, nodeset1=1 )

! constraint for periodic b.c

  call define_constraint ( mesh, input_probdefc, curve1=9, curve2=11, &
    physq=1, discretization='collocation' )

  call problem_definition ( input_probdefc, mesh, problemc )

! create system vectors (solution and right-hand side) for conformation and

  call create ( problemc, solc, solcm1, rhsc )

! initialize vectors with zero stress

  if ( logc == 0 ) then ! standard
    solc(1,1)%u = 1   ! initial cxx
    solc(2,1)%u = 0   ! initial cxy
    solc(3,1)%u = 1   ! initial cyy
  else if ( logc == 1 ) then ! log scheme
    solc(1,1)%u = 0   ! initial sxx
    solc(2,1)%u = 0   ! initial sxy
    solc(3,1)%u = 0   ! initial syy
  end if

! create system matrix for conformation problem

  call create_sysmatrix_structure_base ( sysmatrixc, mesh, problemc )
  call create_sysmatrix_structure_constraint ( sysmatrixc, mesh, problemc )
  call finalize_sysmatrix_structure ( sysmatrixc )

  call create_sysmatrix_data ( sysmatrixc )


! create the structure oldvectors_ve

  call create_oldvectors ( oldvectors_ve, nsysvec=2, nsysvec2=2, nprob=2, &
    nelta=1 )

! store solution vectors and problem structures

  oldvectors_ve%s(1)%p => sol
  oldvectors_ve%s(2)%p => solm1
  oldvectors_ve%s2(1)%p => solc
  oldvectors_ve%s2(2)%p => solcm1
  oldvectors_ve%p(1)%p => problem
  oldvectors_ve%p(2)%p => problemc
  oldvectors_ve%ea(1)%p => eltree_array
  ea => eltree_array ! userelmesh does not have oldvectors as an argument


! open file for writing the drag

  open ( unit=20, file='drag_on_cylinder.out' )


! time stepping

  do step = 1, numtimesteps

!   build (assemble) matrix/vector for gradient/velocity/pressure problem

    call build_vpG

!   build implicit terms of CE

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=divtau_implicit_ce_elem_c, &
      oldvectors=oldvectors_ve, physqrow=[physqvel], physqcol=[physqvel], &
      addmatvec=.true., coefficients=coefficients )

!   embedded Dirichlet b.c on the cylinder

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_open_boundary_eltree, oldvectors=oldvectors_ve, &
      physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress], &
      coefficients=coefficients, addmatvec=.true., elementset=1 )

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      physqrow=[physqvel], physqcol=[physqvel], &
      elemsub=stokes_embedded_dirichlet_eltree, oldvectors=oldvectors_ve, &
      coefficients=coefficients, addmatvec=.true., elementset=1 )

    call build_system ( mesh, problem, sysmatrix, rhsd, &
       physqrow=[physqvel], physqcol=[physqvel], &
       elemsub=implicit_stress_open_boundary_eltree, &
       oldvectors=oldvectors_ve, coefficients=coefficients, &
       addmatvec=.true., elementset=1 )

!   imposed flow rate

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients, oldvectors=oldvectors_ve )

!   periodic boundary condition on velocities

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients, oldvectors=oldvectors_ve )

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve gradient/velocity/pressure problem

    solver_options_u%real_storage=rs_gup
    solver_options_u%integer_storage=is_gup

    call solve_system_ma41 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_u  )

    call copy ( sol, solm1 )


!   build (assemble) matrix and vector for conformation problem

    if ( step == 1 ) then

      coefficients%i(22) = timeint1

      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    else

      coefficients%i(22) = timeint2

      call build_system ( mesh, problemc, sysmatrixc, m2sysvector=rhsc, &
        elemsub=ce_supg_elem_implicit_2nd_order, oldvectors=oldvectors_ve, &
        coefficients=coefficients )

    end if

!   periodical condition on conformation tensor

    call build_system_constraint ( mesh, problemc, sysmatrixc, &
      m2sysvector=rhsc, elemsub=stokes_constr_node_conn, &
      addmatvec=.true. )

    call check ( sysmatrixc )

    call copy ( solc, solcm1 )


!   solve conformation and keep LU decomposition in the loop over components

    solver_options_c%real_storage=rs_c
    solver_options_c%integer_storage=is_c

    do icomp = 1, ncompc
      call solve_system_ma41 ( sysmatrixc, rhsc(icomp,1), solc(icomp,1), luc, &
        solver_options=solver_options_c  )
    end do

    call delete ( luc )  ! remove LU decomposition and rebuild next time step


!   calculate hydrodynamic force on the cylinder

!   NOTE: due to the normal pointing inward, the integrated quantities
!   int_pressure, int_viscous_stress and int_viscoelastic_stress,
!   represent the force on the fluid.

    ! integrate -p*I.n = -p*n
    call integrate ( mesh, problem, int_pressure, &
      elemsub=integration_pressure_drag_eltree, coefficients=coefficients, &
      oldvectors=oldvectors_ve, elementset=1 )

    ! integrate 2*eta*D.n
    call integrate ( mesh, problem, int_viscous_stress, &
      elemsub=integration_viscous_drag_eltree, coefficients=coefficients, &
      oldvectors=oldvectors_ve, elementset=1 )

    ! integrate tau_p.n
    call integrate ( mesh, problemc, int_viscoelastic_stress, &
      elemsub=integration_viscoelastic_drag_eltree, coefficients=coefficients, &
      oldvectors=oldvectors_ve, elementset=1 )

    int_total_stress = &
                  int_pressure + int_viscous_stress + int_viscoelastic_stress

    print '(i6,5es16.8)', step, step*deltat, -int_total_stress
    write(20, '(i6,11es16.8)') step, step*deltat, -int_pressure, &
      -int_viscous_stress, -int_viscoelastic_stress, -int_total_stress

  end do

  call delete_eltree_in_elements

  close(unit=20)


! create data structures for post-processing

! create eltree for all elements crossing the interface

  call create_eltree_in_elements ( numsplit_plot )

! post-processing: create mesh_plot for plotting

  call mesh_convert ( mesh, mesh_plot, userelmesh=userelmesh, warn=.false. )

  call fill_mesh_parts ( mesh_plot )

! problem definition for post_processing

  call create_input_probdef ( mesh_plot, input_probdef_plot, nvec=2 )

  do i = 1, mesh_plot%nelgrp
    input_probdef_plot%elementdof(i)%a(:) = 0
    input_probdef_plot%vec_elementdof(i)%a(:,1) = 2  ! velocity
    input_probdef_plot%vec_elementdof(i)%a(:,2) = 1  ! scalar
  end do

  call problem_definition ( input_probdef_plot, mesh_plot, problem_plot )


! write data for post-processing

  call write_mesh ( mesh, filename='mesh.out' )
  call write_mesh ( mesh_plot, filename='mesh_plot.out' )

  call write_input_probdef ( mesh, input_probdef, filename='probdef.out' )
  call write_input_probdef ( mesh, input_probdefc, filename='probdefc.out' )
  call write_input_probdef ( mesh_plot, input_probdef_plot, &
    filename='probdef_plot.out' )

  call write_coefficients ( coefficients, filename='coefficients.out' )

  open ( unit=10, form='unformatted', file='data.out' )

  write(10) sol%u
  write(10) (solc(i,1)%u, i=1,ncompc)

  close(unit=10)


! delete all data including all allocated memory

  call delete_eltree_in_elements

  call delete ( mesh, mesh_plot )
  call delete ( problem, problem_plot )
  call delete ( input_probdef, input_probdef_plot )
  call delete ( sol, solm1, rhsd )
  call delete ( oldvectors_ve )
  call delete ( problemc )
  call delete ( input_probdefc )
  call delete ( sysmatrixc )
  call delete ( solc, solcm1, rhsc )
  call delete ( coefficients )

  deallocate ( nodes, d, eltree_array, vol )


contains

  subroutine build_vpG

!   build (assemble) matrix and vector for gradient/velocity/pressure problem

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, physqrow=[physqvel,physqpress], &
      physqcol=[physqvel,physqpress], &
      coefficients=coefficients, oldvectors=oldvectors_ve )

!   DEVSS-G
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=devssg_elem, addmatvec=.true., &
      physqrow=[physqgrad,physqvel], physqcol=[physqgrad,physqvel], &
      coefficients=coefficients, oldvectors=oldvectors_ve )

!   set to zero off-diagonal blocks gradient-pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqgrad], physqcol=[physqpress], &
      zeromatvec=.true. )
    call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
      buildvector=.false., physqrow=[physqpress], physqcol=[physqgrad], &
      zeromatvec=.true. )

  end subroutine build_vpG


! create eltree in the elements which are crossed by the interface

  subroutine create_eltree_in_elements ( numsplit )

    integer, intent(in) :: numsplit

    integer :: nod(mesh%element(1)%numnod), elem
    real(dp) :: coor(2,2)
    real(dp) :: minvalvol


!   loop all elements

    minvalvol = 1

    do elem = 1, mesh%grpnumel(1)

!     nodal points

      nod = mesh%topology(1)%a(:,elem)

      if ( all ( d(nod) > split_threshold ) .or. &
                        all ( d(nod) < -split_threshold ) ) then

!       all nodes far from the interface

        cycle

      else

!       element possibly contains an interface, start subdivide

        allocate ( eltree_array(1,elem)%p )

!       fill root node of the eltree

        coor(1,:) = -1._dp  ! lower corner of the reference region (-1,-1,-1)
        coor(2,:) =  1._dp  ! upper corner of the reference region (1,1,1)

        call fill_node_eltree ( eltree_array(1,elem)%p, coor )

!       divide root reference domain into subdomains

        lelgrp = 1
        lelem = elem

        call subdivide ( eltree_array(1,elem)%p, levelset=levelset, &
          numsplit=numsplit, split_threshold=split_threshold, &
          numsplitmin=numsplitmin, mapcoor=mapcoor, submesh=.true., &
          intmesh=.true. )

!       check whether interface is in element

        if ( any( number_of_subelements_vector(eltree_array(1,elem)%p,&
                         &lsign=[-1,1]) == 0 ) ) then

!         no interface, remove eltree

          call delete(eltree_array(1,elem)%p)  ! delete actual eltree
          deallocate( eltree_array(1,elem)%p ) ! delete pointer target
                                               ! (pointer becomes disassociated)

        end if

      end if

      if ( associated(eltree_array(1,elem)%p) ) then

!       element crosses the interface

        vol(elem) = volume ( eltree_array(1,elem)%p, lsign=[1] ) / 4

        if ( vol(elem) <= epsvol ) then
!         remove eltree
          call delete(eltree_array(1,elem)%p)  ! delete actual eltree
          deallocate( eltree_array(1,elem)%p ) ! delete pointer target
                                               ! (pointer becomes disassociated)
          print *, 'element at interface removed, elem = ', elem
        end if

        minvalvol = min ( vol(elem), minvalvol )

      end if

    end do

    !print *, 'minvalvol = ', minvalvol

  end subroutine create_eltree_in_elements


! create elementset of elements which are crossed by the interface

  subroutine create_elementset_from_eltree ( replace )

    integer, intent(in), optional :: replace

    integer :: elements(mesh%grpnumel(1)), numelem, elem


!   loop all elements

    numelem = 0

    do elem = 1, mesh%grpnumel(1)

      if ( associated(eltree_array(1,elem)%p) ) then

!       element crosses the interface

        numelem = numelem + 1

        elements(numelem) = elem

      end if

    end do

!   add elementset

    call add_to_mesh ( mesh, elementset='elements', &
                        elements=elements(1:numelem), replace=replace )

  end subroutine create_elementset_from_eltree


! delete eltree in the elements

  subroutine delete_eltree_in_elements

    integer :: elem

!   loop all elements

    do elem = 1, size(eltree_array,2)

      if ( associated(eltree_array(1,elem)%p) ) then

        call delete(eltree_array(1,elem)%p) ! delete actual eltree
        deallocate(eltree_array(1,elem)%p) ! delete pointer target
                                           !(pointer becomes disassociated)

      end if

    end do

  end subroutine delete_eltree_in_elements


! find internal nodes that are not connected to elements at the interface

  subroutine find_internal_zero_nodes

    integer :: elem, i

!   set internal nodes

    where ( d <= 0._dp )
      nodes = 1
    else where
      nodes = 0
    end where

!   remove nodes connected to elements at the interface

    do elem = 1, mesh%grpnumel(1)
      if ( associated(eltree_array(1,elem)%p) ) then
         nodes(mesh%topology(1)%a(:,elem)) = 0
       end if
    end do

!   count and collect the nodes (overwrite array nodes along the way)

    nnodes = 0
    do i = 1, mesh%nnodes
      if ( nodes(i) == 0 ) cycle
      nnodes = nnodes + 1
      nodes(nnodes) = i
    end do

  end subroutine find_internal_zero_nodes


end program cylinder1
