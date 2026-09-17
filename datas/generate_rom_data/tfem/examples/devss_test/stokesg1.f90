! Stokes problem on a 3D square cavity. Only half of the domain is solved
! due to symmetry conditions.
! Test problem for DEVSS-G
!
! Plot G from the solution vector and as determined from G=grad u^T

program stokesg1

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use devss_elements_m
  use figplot_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    nx = 5,             & ! number of elements in x-direction
    ny = 5,             & ! number of elements in y-direction
    nz = 5,             & ! number of elements in z-direction
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    gintpl = 4,         & ! Q1 gradients
    physqgrad = 1,      & ! physical quantity nr of the gradients
    physqvel = 2,       & ! physical quantity nr of the velocities
    physqpress = 3,     & ! physical quantity nr of the pressures
    gauss = 3             ! 3x3x3 integration of hexahedra

  real(dp), parameter :: &
    lx = 1._dp,          & ! size in x-direction
    ly = 0.5_dp,         & ! size in y-direction
    lz = 1._dp,          & ! size in z-direction
    eta = 1._dp            ! viscosity

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, dudz, gradient_tensor
  type(plot_options_t) :: plot_options
  type(meshgen_options_t) :: meshgen_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

  type(input_probdef_t) :: input_probdef_grad
  type(problem_t) :: problem_grad
  type(sysmatrix_t) :: sysmatrix_grad
  type(sysvector_t) :: sol_grad, rhsd_grad(9)
  type(oldvectors_t) :: oldvectors_grad
  type(coefficients_t) :: coefficients_grad


! variables

  real(dp) :: alpha
  integer :: vertices(8) = [1,3,9,7,19,21,27,25], i


! set alpha

  alpha = eta

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,         0, gintpl,  &
      physqvel, physqpress, 0, physqgrad, gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1:5) = &
    [ eta, 0._dp,   0._dp,   alpha, 0._dp ]
  coefficients%r(6:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, &
    lz=lz, elshape=14, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! plot curves, surfaces and mesh

  plot_options%viewpoint=[1.,0.8,0.4]
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig', surfaces=[3,4,6] )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=3 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 9  ! gradients
  input_probdef%vec_elementdof(1)%a(:,2) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5) = 9  ! unsymmetric tensor

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

  call define_essential ( mesh, input_probdef, surface1=1, surface2=3, &
    physq=physqvel )
  call define_essential ( mesh, input_probdef, surface1=4, physq=physqvel,  &
    degfd=[0,1,0] )
  call define_essential ( mesh, input_probdef, surface1=5, surface2=6, &
    physq=physqvel )
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    surface1=1, surface2=3, physq=physqvel, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=4, physq=physqvel, degfd=2, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=5, surface2=6, physq=physqvel, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    surface1=6, physq=physqvel, degfd=1, value=1._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqpress, value=0._dp )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector for gradient/velocity/pressure problem

! stokes velocity/pressure
  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, physqrow=[2,3], physqcol=[2,3], &
    coefficients=coefficients )

! DEVSS-G
  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=devssg_elem, addmatvec=.true., &
    physqrow=[1,2], physqcol=[1,2], coefficients=coefficients )

! set to zero off-diagonal blocks gradient-pressure
  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., physqrow=[1], physqcol=[3], zeromatvec=.true. )
  call build_system ( mesh, problem, sysmatrix, rhsd, addmatvec=.true., &
    buildvector=.false., physqrow=[3], physqcol=[1], zeromatvec=.true. )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! post-processing

  call create_vector ( problem, velocity, physq=physqvel )
  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, dudz, vec=4 )

  call extract_physvector ( mesh, problem, sol, velocity )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=3
  call derive_vector ( mesh, problem, dudz, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )


! write to a fig file for plotting

  plot_options%printlabels=.false.

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity, surfaces=[4,6] )
  call plot_points_curves ( plot_options, mesh, 'velocity.fig', append=.true. )

  call plot_color_fill ( plot_options, mesh, problem, 'velocity_color.fig', &
    vector=velocity, degfd=1, surfaces=[4] )
  call plot_points_curves ( plot_options, mesh, 'velocity_color.fig',  &
   append=.true. )

  call plot_color_contour ( plot_options, mesh, problem, 'dudz_contour.fig', &
    vector=dudz, surfaces=[3,4,6] )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure, surfaces=[3,4,6] )


! gradient tensor

  call create_vector ( problem, gradient_tensor, vec=5 )

  call derive_vector ( mesh, problem, gradient_tensor, &
    elemsub=deriv_gradient_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call write_tensor_vtk ( mesh, problem, filename='stokesg.vtk', &
    dataname='gradient_tensor', vector=gradient_tensor, symmetric=.false. )

  call delete ( gradient_tensor )


!*******************
! Gradient problem *
!*******************

! problem definition for gradients

  call create_input_probdef ( mesh, input_probdef_grad, nvec=2, nphysq=1 )

  input_probdef_grad%vec_elementdof(1)%a(:,1) = 0
  input_probdef_grad%vec_elementdof(1)%a(vertices,1) = 1  ! gradient component

  input_probdef_grad%physq = [1]
  input_probdef_grad%probnr = 2

  call problem_definition ( input_probdef_grad, mesh, problem_grad )

! create system vectors (solution and right-hand side)

  call create ( problem_grad, sol_grad )
  call create ( problem_grad, rhsd_grad )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix_grad, mesh, problem_grad, &
    symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix_grad )

! fill coefficients
  call create_coefficients ( coefficients_grad, ncoefi=160, ncoefr=100 )

  coefficients_grad%i = coefficients%i

  coefficients_grad%r = 0


! oldvectors

  call create_oldvectors ( oldvectors_grad, nsysvec=1, nprob=1 )

  oldvectors_grad%s(1)%p => sol
  oldvectors_grad%p(1)%p => problem

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem_grad, sysmatrix_grad, msysvector=rhsd_grad,&
    elemsub=gradient_from_velocity_elem, coefficients=coefficients_grad, &
    oldvectors=oldvectors_grad )

  do i = 1, 9

    call add_effect_of_essential_to_rhs ( problem_grad, sysmatrix_grad, &
      sol_grad, rhsd_grad(i) )

    call solve_system_ma57 ( sysmatrix_grad, rhsd_grad(i), sol_grad )

    call transfer_data ( mesh, problem1=problem_grad, &
      sysvector1=sol_grad, problem2=problem, sysvector2=sol, &
      degfd1=[1], degfd2=[i] )

  end do

! gradient tensor

  call create_vector ( problem, gradient_tensor, vec=5 )

  call derive_vector ( mesh, problem, gradient_tensor, &
    elemsub=deriv_gradient_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call write_tensor_vtk ( mesh, problem, filename='stokesg1.vtk', &
    dataname='gradient_tensor', vector=gradient_tensor, symmetric=.false. )

  call delete ( gradient_tensor )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( sysmatrix )
  call delete ( velocity, pressure, dudz )
  call delete ( coefficients )
  call delete ( oldvectors )

  call delete ( problem_grad )
  call delete ( input_probdef_grad )
  call delete ( sol_grad )
  call delete ( rhsd_grad )
  call delete ( sysmatrix_grad )
  call delete ( oldvectors_grad )

end program stokesg1
