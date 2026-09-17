! Stokes problem on a unit square with 3D velocities (developed flow).
! Lid-driven cavity flow in (x,y).
! Imposed flow rate in z-direction.
! Test problem for DEVSS-G
! Plot G from the solution vector and as determined from G=grad u^T


program stokesg2

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use devss_elements_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    gintpl = 4,         & ! Q1 gradients
    physqgrad = 1,      & ! physical quantity nr of the gradients
    physqvel = 2,       & ! physical quantity nr of the velocities
    physqpress = 3,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=20,              & ! number of elements in x
    ny=20                 ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp,     & ! viscosity
    U = 1._dp,       & ! velocity of the "lid"
    flowrate = 1._dp   ! flow rate in third direction

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, dwdx, gammadot, Dtensor, Ltensor
  type(vector_t) :: tauxz, tautensor, Gtensor
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

  type(input_probdef_t) :: input_probdef_grad
  type(problem_t) :: problem_grad
  type(sysmatrix_t) :: sysmatrix_grad
  type(sysvector_t) :: sol_grad, rhsd_grad(6)
  type(oldvectors_t) :: oldvectors_grad
  type(coefficients_t) :: coefficients_grad

! variables

  real(dp) :: alpha
  integer :: vertices(4) = [1,3,5,7], i


! set alpha

  alpha = eta

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,         0, gintpl,  &
      physqvel, physqpress, 0, physqgrad, gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%i(67) = 1  ! 3D velocity

  coefficients%r = 0
  coefficients%r(1) = eta
  coefficients%r(4) = alpha
  coefficients%r(6) = flowrate

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

! make domain into a surface for the flow rate constraint
  call add_to_mesh ( mesh, surfacefromgroups=[1] )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=6, nphysq=3 )

  input_probdef%vec_elementdof(1)%a(:,1) = 0
  input_probdef%vec_elementdof(1)%a(vertices,1) = 6  ! gradients
  input_probdef%vec_elementdof(1)%a(:,2) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,3) = 0
  input_probdef%vec_elementdof(1)%a(vertices,3) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,4) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,5) = 6  ! symmetric tensor
  input_probdef%vec_elementdof(1)%a(:,6) = 9  ! unsymmetric tensor

  input_probdef%physq = [1,2,3]
  input_probdef%probnr = 1

  call define_essential ( mesh, input_probdef, &
    curve1=1, curve2=4, physq=physqvel )
  call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

! constraint for flow rate

  call define_constraint ( mesh, input_probdef, &
    physq=physqvel, surface1=1, nglobalc=1 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, physq=physqvel, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=physqvel, degfd=1, value=U )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=physqpress, value=0._dp )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

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

 call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
    constraint1=1, elemsub=stokes_constr_flowr_surface, &
    addmatvec=.true., coefficients=coefficients )

  call check ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! post-processing

  call create_vector ( problem, velocity, physq=2 )
  call create_vector ( problem, pressure, vec=4 )
  call create_vector ( problem, dwdx, vec=4 )
  call create_vector ( problem, gammadot, vec=4 )
  call create_vector ( problem, Dtensor, vec=5 )
  call create_vector ( problem, Ltensor, vec=6 )
  call create_vector ( problem, tauxz, vec=4 )
  call create_vector ( problem, tautensor, vec=5 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=7

  call derive_vector ( mesh, problem, dwdx, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=11

  call derive_vector ( mesh, problem, gammadot, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Dtensor, elemsub=stokes_D_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Ltensor, elemsub=stokes_gradu_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=3

  call derive_vector ( mesh, problem, tauxz, elemsub=stokes_stress, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, tautensor, elemsub=stokes_stress_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='stokesg2.vtk' )

  call write_vector_vtk ( mesh, problem, filename='stokesg2.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true., &
    assume3D=.true. )

  call write_scalar_vtk ( mesh, problem, vector=dwdx, filename='stokesg2.vtk', &
    dataname='dwdx', append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='stokesg2.vtk', dataname='gammadot', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='stokesg2.vtk', &
    dataname='D', vector=Dtensor, append=.true., assume3D=.true. )

  call write_tensor_vtk ( mesh, problem, filename='stokesg2.vtk', &
    dataname='L', vector=Ltensor, append=.true., assume3D=.true., &
    symmetric=.false. )

  call write_scalar_vtk ( mesh, problem, vector=tauxz, &
    filename='stokesg2.vtk', dataname='tauxz', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='stokesg2.vtk', &
    dataname='tau', vector=tautensor, append=.true., assume3D=.true. )

! gradient tensor

  call create_vector ( problem, Gtensor, vec=6 )

  call derive_vector ( mesh, problem, Gtensor, &
    elemsub=deriv_gradient_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call write_tensor_vtk ( mesh, problem, filename='stokesg2.vtk', &
    dataname='Gtensor', vector=Gtensor, symmetric=.false., assume3D=.true., &
    append=.true. )

  call delete ( Gtensor )


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
  call create_coefficients ( coefficients_grad, ncoefi=150, ncoefr=100 )

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

  do i = 1, 6

    call add_effect_of_essential_to_rhs ( problem_grad, sysmatrix_grad, &
      sol_grad, rhsd_grad(i) )

    call solve_system_ma57 ( sysmatrix_grad, rhsd_grad(i), sol_grad )

    call transfer_data ( mesh, problem1=problem_grad, &
      sysvector1=sol_grad, problem2=problem, sysvector2=sol, &
      degfd1=[1], degfd2=[i] )

  end do

! gradient tensor

  call create_vector ( problem, Gtensor, vec=6 )

  call derive_vector ( mesh, problem, Gtensor, &
    elemsub=deriv_gradient_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call write_tensor_vtk ( mesh, problem, filename='stokesg2.vtk', &
    dataname='Gtensor2', vector=Gtensor, symmetric=.false., assume3D=.true., &
    append=.true. )

  call delete ( Gtensor )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, dwdx, gammadot, Dtensor, Ltensor )
  call delete ( tauxz, tautensor )!, Gtensor )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

  call delete ( problem_grad )
  call delete ( input_probdef_grad )
  call delete ( sol_grad )
  call delete ( rhsd_grad )
  call delete ( sysmatrix_grad )
  call delete ( oldvectors_grad )

end program stokesg2
