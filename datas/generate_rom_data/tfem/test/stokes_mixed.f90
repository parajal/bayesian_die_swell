module functions_m

  use kind_defs_m

  implicit none

contains

  function func ( nr, x )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    select case(nr)
      case(1)
        func = 0
      case(2)
        func = 1
      case(3)
        func = 100
      case default
        write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
        stop
    end select

  end function func

  function vfunc ( n, nr, x )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    write(*,*) 'function vfunc has not been defined'
    vfunc = 0
    stop

  end function vfunc

end module functions_m
module meshgen_contraction_m

  use meshgen_m

  implicit none

  type ic_t
    integer :: elshape=6, p=0, l=1
    integer :: n1=5, n2=5, n3=5, n4=3, n5=10, m1=5, m2=7
    real(dp) :: refine=1.0_dp
    real(dp) :: factor1=0.5_dp, factor2=0.2_dp
    real(dp) :: L1=10, L2=1, L3=1, L4=2, L5=10
    real(dp) :: H1=1.5_dp, H2=2.5_dp, H3=1.0_dp
  end type ic_t

contains

! 2D contraction

  subroutine contraction2D ( ic, mesh )

    type(ic_t), intent(in) :: ic
    type(mesh_t), intent(out) :: mesh

    type(meshgen_options_t) :: mesh_options
    type(mesh_t) :: mesh1, mesh2, mesh3

    real(dp) :: x2d(4,2)

    integer :: elshape, n1, n2, n3, n4, n5, m1, m2
    real(dp) :: factor1, factor2, L1, L2, L3, L4, L5, H1, H2, H3
!
!  contraction build from 7 submeshes:
!
!
!         ------------------
!         |        |       |
!  m2  H2 |   2    |   4   |
!         |        |       |
!         |--------|\____  |
!         |        |     \_|---------------------------
!  m1  H1 |   1    |   3   |    |    |                |  H3
!         |        |       | 5  | 6  |      7         |
!       0 ---------------------------------------------
!        -L1      -L2      0    L3   L4               L5
!             n1      n2     n3   n4        n5
!
! NOTE: one single group is created. Use nogroupmerge=.false. in heading of
! mesh_merge to obtain 7 separate element groups.


    elshape = ic%elshape
    n1 = nint( ic%n1 * ic%refine )
    n2 = nint( ic%n2 * ic%refine )
    n3 = nint( ic%n3 * ic%refine )
    n4 = nint( ic%n4 * ic%refine )
    n5 = nint( ic%n5 * ic%refine )
    m1 = nint( ic%m1 * ic%refine )
    m2 = nint( ic%m2 * ic%refine )
    factor1 = ic%factor1
    factor2 = ic%factor2
    L1 = ic%L1
    L2 = ic%L2
    L3 = ic%L3
    L4 = ic%L4
    L5 = ic%L5
    H1 = ic%H1
    H2 = ic%H2
    H3 = ic%H3

! 1

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=m1, ox=-L1, &
    lx=L1-L2, ly=H1, p=ic%p, l=ic%l )

  call quadrilateral2d ( mesh1, mesh_options )

! 2

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n1, ny=m2, ox=-L1, &
    oy=H1, lx=L1-L2, ly=H2, p=ic%p, l=ic%l )

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, curve1=3, curve2=-1 )
  call delete ( mesh1, mesh2 )

! 3

  x2d(1,:) = [-L2,0._dp]
  x2d(2,:) = [0,0]
  x2d(3,:) = [0._dp,H3]
  x2d(4,:) = [-L2,H1]

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=m1, &
    x2d=x2d, regionshape=2, ratio=[6,6,5,0], p=ic%p, l=ic%l )

  mesh_options%factor = factor1

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh3, mesh2, mesh1, curve1=2, curve2=-4 )
  call delete ( mesh3, mesh2 )

! 4

  x2d(1,:) = [-L2,H1]
  x2d(2,:) = [0._dp,H3]
  x2d(3,:) = [0._dp,H1+H2]
  x2d(4,:) = [-L2,H1+H2]

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n2, ny=m2, &
    x2d=x2d, regionshape=2, ratio=[6,5,5,0], p=ic%p, l=ic%l )

  mesh_options%factor = factor1
  mesh_options%factor(2) = factor2

  call quadrilateral2d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, curves1=[10,5], curves2=[-1,-4] )
  call delete ( mesh1, mesh2 )

! 5

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n3, ny=m1, &
    lx=L3, ly=H3, ratio=[5,6,6,5], p=ic%p, l=ic%l )

  mesh_options%factor = factor1

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh3, mesh1, mesh2, curve1=9, curve2=-4 )
  call delete ( mesh3, mesh1 )

! 6

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n4, ny=m1, ox=L3, &
    lx=L4-L3, ly=H3, ratio=[0,0,0,5], p=ic%p, l=ic%l )

  mesh_options%factor = factor1

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh2, mesh1, mesh3, curve1=14, curve2=-4 )
  call delete ( mesh2, mesh1 )

! 7

  call set_mesh_options ( mesh_options, elshape=elshape, nx=n5, ny=m1, ox=L4, &
    lx=L5-L4, ly=H3, p=ic%p, l=ic%l )

  call quadrilateral2d ( mesh1, mesh_options )

  call mesh_merge ( mesh3, mesh1, mesh, curve1=17, curve2=-4 )
  call delete ( mesh3, mesh1 )

  end subroutine contraction2D

end module meshgen_contraction_m

! Mixed Stokes problem for a contraction flow
! Q2/Q1 Taylor-Hood for velocity-pressure, using main mesh for interpolation.
! Qp for stresses using the blend mesh 1 for interpolation.
! Blended mesh:
!  main mesh: 9-node quadrilaterals
!  blend mesh 1: (p+1)^2-node high-order quadrilaterals, equidistant nodes

program mixed_stokes1

  use tfem_m
  use stokes_elements_m
  use mixed_stokes_elements_2D_m
  use meshgen_contraction_m
  use io_utils_m
  use functions_m
  use hsl_ma57_m

  implicit none

! constants

  logical, parameter :: &
    CR = .false.           ! use mesh in lecture computational rheology

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    tintpl = 20,        & ! Qp stresses
    p = 4,              & ! polynomial order stresses
    physqvel = 2,       & ! physical quantity nr of the velocities
    physqpress = 3,     & ! physical quantity nr of the pressures
    physqstress = 1,    & ! physical quantity nr of the stress
    gauss = max(p+1,3), & ! p+1 x p+1 integration of quads
    gaussb = max(p+1,3),& ! number of Gauss points on the boundary
    inttype = 3,        & ! numerical rules for integration
    nn0 = 9,            & ! number of nodes in main mesh element
    nn1 = (p+1)**2        ! number of nodes in blend mesh element

  integer, parameter :: &
    funcnr = 3,         & ! function number of traction natural boundary
    direction = 1         ! direction of the traction component

  real(dp), parameter :: &
    refine = 0.5_dp, & ! refinement factor for mesh
    eta = 1._dp       ! viscosity

! definitions

  type(mesh_t) :: meshQ2, meshQp, mesh, mesh_plot_Qp
  type(coefficients_t) :: coefficients
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: pressure, vorticity, divergence
  type(oldvectors_t) :: oldvectors
  type(solver_options_ma57_t) :: so
  type(ic_t) :: ic


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,          tintpl, 0,     0,  &
      physqvel, physqpress, physqstress, 0, gauss,  &
      gaussb ]
  coefficients%i(12:) = 0

  coefficients%i(16) = funcnr
  coefficients%i(17) = direction
  coefficients%i(40) = inttype
  coefficients%i(82) = p ! stress polynomial order

  coefficients%r = 0
  coefficients%r(1) = 0    ! zero direct viscosity
  coefficients%r(5) = eta  ! viscosity in mixed element

  coefficients%func => func


! Create blended mesh

  if ( CR ) then
    ic = ic_t ( elshape=6, H1=1._dp, H2=3.0_dp, refine=refine )
  else
    ic = ic_t ( elshape=6, refine=refine )
  end if

  call generate_mesh ( meshQ2 ) ! velocity-pressure mesh

  if ( CR ) then
    ic = ic_t ( elshape=102, H1=1._dp, H2=3.0_dp, p=p, l=0, refine=refine )
  else
    ic = ic_t ( elshape=102, p=p, l=0, refine=refine )
  end if

  call generate_mesh ( meshQp ) ! stress tensor mesh

  call mesh_convert ( meshQ2, blendmesh=meshQp, mesh=mesh ) ! blended mesh

  call fill_mesh_parts ( mesh )

! Create plot mesh for stresses

  call mesh_convert ( meshQp, mesh_plot_Qp, warn=.false., &
    elementshapes='spectraltolinear' )

  call fill_mesh_parts ( mesh_plot_Qp )

! output mesh (print and plot)

  call printinfo ( mesh, printlevel=6 )

  call delete ( meshQ2, meshQp )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=3 )

  input_probdef%vec_elementdof(1)%a = 0           ! initialize to zero

  input_probdef%vec_elementdof(1)%a(1:nn0,1) = 2         ! velocity (Q2 mesh)
  input_probdef%vec_elementdof(1)%a([1,3,5,7],2) = 1     ! pressure (Q2 mesh)
  input_probdef%vec_elementdof(1)%a(nn0+1:nn0+nn1,3) = 3 ! stress (Qp mesh)
  input_probdef%vec_elementdof(1)%a(1:nn0,4) = 1         ! scalar (Q2 mesh)

  input_probdef%physq = [3,1,2]

! inflow, outflow, center line
  call define_essential ( mesh, input_probdef, curves=[20,22,24], &
    physq=physqvel, degsfd=[2] )

! wall
  call define_essential ( mesh, input_probdef, curves=[23], physq=physqvel )

  call problem_definition ( input_probdef, mesh, problem )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0


! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )


! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    physqrow=[physqvel,physqpress], physqcol=[physqvel,physqpress], &
    coefficients=coefficients )

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=mixed_stokes_elem, addmatvec=.true., coefficients=coefficients, &
    physqrow=[physqstress,physqvel], physqcol=[physqstress,physqvel] )

  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., zeromatvec=.true., coefficients=coefficients, &
    physqrow=[physqstress], physqcol=[physqpress] )

  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., zeromatvec=.true., coefficients=coefficients, &
    physqrow=[physqpress], physqcol=[physqstress] )

  call check (sysmatrix)

  call add_boundary_elements ( mesh, problem, rhsd, curve=24, &
    physq=[physqvel], elemsub=stokes_natboun_curve, coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol, solver_options=so )

  print *, 'sqrt(sum(sol%u**2)/sol%n)'
  print *, sqrt(sum(sol%u**2)/sol%n)

! post processing using derive

  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, vorticity, vec=4 )
  call create_vector ( problem, divergence, vec=4 )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  print *, 'sqrt(sum(pressure%u**2)/pressure%n)'
  print *, sqrt(sum(pressure%u**2)/pressure%n)

  coefficients%i(13)=5
  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  print *, 'sqrt(sum(vorticity%u**2)/vorticity%n)'
  print *, sqrt(sum(vorticity%u**2)/vorticity%n)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh, mesh_plot_Qp )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( pressure, vorticity, divergence )
  call delete ( coefficients )
  call delete ( oldvectors )

contains

  subroutine generate_mesh ( mesh )

    type(mesh_t), intent(out) :: mesh

    call contraction2D ( ic, mesh )

    call add_to_mesh ( mesh, curve=[1,8,13,16,19] )      ! c22 center line
    call add_to_mesh ( mesh, curve=[21,18,15,11,12,6] )  ! c23 wall
    call add_to_mesh ( mesh, curve=[7,4] )               ! c24 inflow boundary
    call add_to_mesh ( mesh, curve=[-10,-15] )           ! c25 curve for plot

  end subroutine generate_mesh

end program mixed_stokes1

