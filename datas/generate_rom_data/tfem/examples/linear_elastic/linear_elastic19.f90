! Linear elastic problem on a square. Manufactured solution.
! Mixed displacement-pressure formulation. Axisymmetric with torsion.

program linear_elastic19

  use tfem_m
  use hsl_ma57_m
  use linear_elastic_elements_m
  use linear_elastic_post_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 displacements
    pintpl = 4,         & ! Q1 pressures
    physqdisp = 1,      & ! physical quantity nr of the displacements
    physqpress = 2,     & ! physical quantity nr of the pressure
    gauss = 3,          & ! 3x3 integration of quads
    nz=20,              & ! number of elements in z
    nr=20                 ! number of elements in r

  real(dp), parameter :: &
    lz = 1._dp,   &     ! length of domain in z-direction
    lr = 1._dp,   &     ! length of domain in r-direction
    or = 0._dp,   &     ! length of domain in r-direction
    Emod = 2._dp, &     ! Young's modulus
    nu = 0.3_dp         ! Poisson's ratio

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(plot_options_t) :: plot_options
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors
  type(lemodel_t) :: lemodel

  integer :: vertices(4) = [1,3,5,7]

  real(dp) :: lambda, mu
  real(dp) :: area_domain(1), L2_error_norm(2)

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(2) = pintpl
  coefficients%i(4) = 0   !  0: compressible, 1: incompressible
  coefficients%i(6) = physqdisp
  coefficients%i(7) = physqpress
  coefficients%i(10) = gauss
  coefficients%i(11) = gauss
  coefficients%i(14) = 2 ! vfuncnr for exact body force
  coefficients%i(23) = 1 ! axisymmetric, cylindrical coordinates
  coefficients%i(67) = 1 ! 3D displacement (torsion)
  coefficients%i(63) = 1 ! vfuncnr for exact displacement
  coefficients%i(64) = 1 ! funcnr for pressure

  coefficients%r = 0
  coefficients%r(1) = Emod
  coefficients%r(3) = nu

  coefficients%func => func
  coefficients%vfunc => vfunc

! set material parameters

  call set_linear_elastic_material ( coefficients, coorsys=1, lemodel=lemodel )
  lambda = lemodel%lambda
  mu = lemodel%mu

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nz
  meshgen_options%ny = nr
  meshgen_options%lx = lz
  meshgen_options%ly = lr
  meshgen_options%oy = or

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! displacement
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(vertices,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar
  input_probdef%vec_elementdof(1)%a(:,4) = 6  ! symmetric tensor

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, physq=1, vfunc=vfunc, vfuncnr=1 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=linear_elastic_elem, coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call integrate ( mesh, problem, area_domain, &
    elemsub=linear_elastic_integrate_volume, &
    coefficients=coefficients, oldvectors=oldvectors )

  call integrate ( mesh, problem, L2_error_norm, &
    elemsub=linear_elastic_L2_norm, &
    coefficients=coefficients, oldvectors=oldvectors )

  print *, ' ||u-u_exact||_2 ', sqrt(L2_error_norm(1)/area_domain(1))
  print *, ' ||p-p_exact||_2 ', sqrt(L2_error_norm(2)/area_domain(1))

! post-processing

  call postprocessing_linear_elastic ( mesh, problem, 'le19.vtk', sol, &
    vec_scalar=3, vec_tensor=4, coefficients=coefficients, assume3D=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

contains

  function func ( nr, xin )
    integer, intent(in) :: nr
    real(dp), intent(in), dimension(:) :: xin
    real(dp) :: func

    real(dp) :: x, y

    x = xin(1); y = xin(2)

    select case ( nr )
      case(1)
!       exact pressure
        func = -lambda*(x**4-6*y**4/5+x**2*y**3+x**2+3*x**2*y/2-1)
      case(2)
!       divergence
        func = x**4-6*y**4/5+x**2*y**3+x**2+3*x**2*y/2-1
      case default
        write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
        stop
    end select

  end function func

  function vfunc ( n, nr, xin )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: xin
    real(dp), dimension(n) :: vfunc

    real(dp) :: x, y

    x = xin(1); y = xin(2)

    select case ( nr )
      case(1)
        vfunc = [ x + x**5/5 + y**2/2 + x**3*y**3/3, &
                 -y -y**5/5 + x**2/2*y + x**2*y**2/2, &
                  y + y**3/3 + x**2*y**2/2 ]  ! exact displacement
      case(2)
        vfunc = -(lambda+mu)* [ &
                    4*x**3 + 2*x*y**3 + 2*x + 3*x*y, &
                    -24*y**3/5+3*x**2*y**2+3*x**2/2, 0._dp ] - &
                mu * [ 4*x**3+2*x*y**3+2+ 3*x**3*y, &
                               y+y**2-24*y**3/5+3*x**2/2, &
                       y**2 + 8*y/3 + 3*x**2/2 ] ! rhs
      case default
        write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
        stop
    end select

  end function vfunc

end program linear_elastic19
