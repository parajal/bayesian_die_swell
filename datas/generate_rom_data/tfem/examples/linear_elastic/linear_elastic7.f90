! Linear elastic problem on a rectangle. Stretching and rotating a bar.
! Standard displacement formulation. Axisymmetric with torsion.

program linear_elastic7

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
    physqdisp = 1,      & ! physical quantity nr of the displacements
    gauss = 3,          & ! 3x3 integration of quads
    nz=20,              & ! number of elements in z
    nr=20                 ! number of elements in r

  real(dp), parameter :: &
    lz = 10._dp,  &     ! length of domain in z-direction
    lr = 2._dp,  &      ! length of domain in r-direction
    Emod = 2._dp, &     ! Young's modulus
    nu = 0.3_dp, &      ! Poisson's ratio
    thetal = 2e-2_dp, & ! rotation at z = L
    ul = 2e-2_dp        ! displacement at z = L

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

  logical :: uniform = .false. ! uniform deformation, no edge effects


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(6) = physqdisp
  coefficients%i(10) = gauss
  coefficients%i(11) = gauss
  coefficients%i(23) = 1 ! axisymmetric, cylindrical coordinates
  coefficients%i(67) = 1 ! 3D displacement (torsion)


  coefficients%r = 0
  coefficients%r(1) = Emod
  coefficients%r(3) = nu

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nz
  meshgen_options%ny = nr
  meshgen_options%lx = lz
  meshgen_options%ly = lr

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=1 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! displacement
  input_probdef%vec_elementdof(1)%a(:,2) = 1  ! scalar
  input_probdef%vec_elementdof(1)%a(:,3) = 6  ! symmetric tensor

  input_probdef%physq = [1]

  call define_essential ( mesh, input_probdef, curve1=1, physq=1, &
    degfd=[0,1,1] )

  if ( uniform ) then
    call define_essential ( mesh, input_probdef, curves=[2,4], physq=1, &
      degfd=[1,0,1] )
  else
    call define_essential ( mesh, input_probdef, curves=[2,4], physq=1 )
  end if

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, physq=1, degsfd=[2,3], value=0._dp )

  if ( uniform ) then
!   uniform deformation, no edge effects
    call fill_sysvector ( mesh, problem, sol, &
      curve1=2, physq=1, degsfd=[1,3], vfunc=displacement_theta, vfuncnr=1 )
    call fill_sysvector ( mesh, problem, sol, &
      curve1=4, physq=1, degsfd=[1,3], value=0._dp )
  else
!   edge effects
    call fill_sysvector ( mesh, problem, sol, &
      curve1=2, physq=1, degsfd=[1,3], vfunc=displacement_theta, vfuncnr=1 )
    call fill_sysvector ( mesh, problem, sol, &
      curve1=2, physq=1, degfd=2, value=0._dp )
    call fill_sysvector ( mesh, problem, sol, &
      curve1=4, physq=1, value=0._dp )
  end if

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=linear_elastic_elem, coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! post processing

  call postprocessing_linear_elastic ( mesh, problem, 'le7.vtk', sol, &
    vec_scalar=2, vec_tensor=3, coefficients=coefficients, assume3D=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )

contains

  function displacement_theta ( n, nr, x )
    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: displacement_theta

    select case ( nr )
      case(1)
        displacement_theta = [ ul, thetal * x(2) ]
      case default
        write(*,'(/a,i0/)') &
            'Error displacement_theta: wrong function number: ', nr
        stop
    end select

  end function displacement_theta

end program linear_elastic7
