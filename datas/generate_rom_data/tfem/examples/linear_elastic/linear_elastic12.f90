! Linear elastic problem on a square. Manufactured solution.
! Standard displacement formulation. 2D, plane strain or plane stress.

program linear_elastic12

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
    nx=20,              & ! number of elements in x
    ny=40                 ! number of elements in y

  real(dp), parameter :: &
    lx = 1.0_dp,  &  ! length of domain in x-direction
    ly = 1.0_dp,  &  ! length of domain in y-direction
    Emod = 2._dp, &  ! Young's modulus
    nu = 0.3_dp      ! Poisson's ratio

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd, sol_exact
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(lemodel_t) :: lemodel
  type(vector_t) :: scalar, scalar_exact

  real(dp) :: lambda, mu
  real(dp) :: area_domain(1), L2_error_norm(1)

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(1) = uintpl
  coefficients%i(3) = 1   !  0: plane strain, 1: plane stress
  coefficients%i(6) = physqdisp
  coefficients%i(10) = gauss
  coefficients%i(11) = gauss
  coefficients%i(14) = 2  ! vfuncnr for exact body force
  coefficients%i(63) = 1  ! vfuncnr for exact displacement

  coefficients%r = 0
  coefficients%r(1) = Emod
  coefficients%r(3) = nu

  coefficients%func => func
  coefficients%vfunc => vfunc

! set material parameters

  call set_linear_elastic_material ( coefficients, coorsys=0, lemodel=lemodel )
  if ( coefficients%i(3) == 1 ) then
!   plane stress
    lambda = lemodel%lambda_bar
  else
    lambda = lemodel%lambda
  end if
  mu = lemodel%mu

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%lx = lx
  meshgen_options%ly = ly

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=1 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! displacement
                 1,1,1,1,1,1,1,1,1,    &  ! scalar
                 4,4,4,4,4,4,4,4,4 ], &  ! tensor with 4 components
                 [9,3] )

  input_probdef%physq = [1]

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, sol_exact )
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

! post processing

  call postprocessing_linear_elastic ( mesh, problem, 'le12.vtk', sol, &
    vec_scalar=2, vec_tensor=3, coefficients=coefficients, assume33=.true. )

! post processing

  call fill_sysvector ( mesh, problem, sol_exact, &
    node1=1, node2=mesh%nnodes, vfunc=vfunc, vfuncnr=1 )

  call postprocessing_linear_elastic ( mesh, problem, 'le12e.vtk', sol_exact, &
    vec_scalar=2, vec_tensor=3, coefficients=coefficients, assume33=.true. )

  call create_vector ( problem, scalar, vec=2 )
  call create_vector ( problem, scalar_exact, vec=2 )

  call fill_vector ( mesh, problem, scalar_exact, &
    node1=1, node2=mesh%nnodes, func=func, funcnr=2 )

  coefficients%i(13)=1

  call derive_vector ( mesh, problem, scalar, &
    elemsub=linear_elastic_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  print *, ' ||divu-divu_exact||_max ', maxval(abs(scalar%u-scalar_exact%u))

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

    real(dp) :: x, y, tmp

    x = xin(1); y = xin(2)

    select case ( nr )
      case(1)
        func = -lambda*(x**4-y**4+x**2*y**3+x**2*y) ! exact pressure
      case(2)
        tmp = x**4-y**4+x**2*y**3+x**2*y
        if (lemodel%plane_stress) then
          func = (1-2*lemodel%nu)/(1-lemodel%nu)*tmp  ! divergence
        else
          func = tmp ! divergence
        end if
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
                  -y -y**5/5 + x**2/2 + x**2*y**2/2 ]  ! exact displacement
      case(2)
        vfunc = -(lambda+mu)* [ &
                    4*x**3 + 2*x*y**3 + 2*x*y, &
                    -4*y**3+3*x**2*y**2+x**2 ] - &
                mu * [ 4*x**3+2*x*y**3+1+ 2*x**3*y, &
                               1+y**2-4*y**3+x**2 ] ! rhs
      case default
        write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
        stop
    end select

  end function vfunc

end program linear_elastic12
