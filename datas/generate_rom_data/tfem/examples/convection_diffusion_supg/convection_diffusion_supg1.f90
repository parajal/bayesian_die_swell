! unsteady convection problem with a source term
! constant velocity field
! SUPG for stabilisation
! second-order time integration
! linear interpolation of q

program convection_diffusion_supg1

  use tfem_m
  use hsl_ma41_m
  use convection_diffusion_supg_elements_m
  use io_utils_m
  use subs_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! velocity interpolation
    pintpl = 4,         & ! pressure interpolation
    qintpl = 4,         & ! q interpolation
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    timeint1 = 1,       & ! first-order implicit Euler
    timeint2 = 2,       & ! second-order implicit Gear
    gauss = 5,          & ! order of integration of quads
    coorsys = 0,        & ! planar coordinate system
    nx = 15,            & ! number of elements in x-direction
    ny = 15,            & ! number of elements in y-direction
    ibeta = 1,          & ! method to compute beta in SUPG
    numtimesteps=400      ! number of time steps

  real(dp), parameter :: &
    gamma = 1._dp,       & ! parameter in the convection-diffusion equation
    deltat = 1.e-2_dp,   & ! time step
    ox = 0._dp,          & ! x-coordinate of origin of the mesh
    oy = 0._dp,          & ! y-coordinate of origin of the mesh
    lx = 1._dp,          & ! length of the domain in x-direction
    ly = 1._dp,          & ! length of the domain in y-direction
    rs = 1.4_dp,         & ! real_storage velocity-pressure LU (HSL)
    is = 1.5_dp            ! integer_storage velocity-pressure LU (HSL)

! variables

  logical :: buildmatrix

  integer :: step

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefq
  type(problem_t), target :: problem, problemq
  type(sysmatrix_t) :: sysmatrixq
  type(sysvector_t), target :: sol, solq, solq_n, solq_nm1, rhsdq
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors
  type(lu_ma41_t) :: lu
  type(solver_options_ma41_t) :: solver_options

  character(len=199) :: filename

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=450, ncoefr=400 )

  coefficients%i = 0
  coefficients%i(1:10) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss ]
  coefficients%i(23)  = coorsys
  coefficients%i(40)  = 3 ! numerical table for Gauss
  coefficients%i(401) = timeint1
  coefficients%i(402) = qintpl
  coefficients%i(403) = 1 ! 0: Galerkin, 1: SUPG
  coefficients%i(406) = 2 ! f-term given by function
  coefficients%i(407) = 1 ! funcnr for f-term
  coefficients%i(411) = ibeta  ! method to compute beta in SUPG
!  coefficients%i(412) = 1 ! Courant-dependent tau in SUPG

  coefficients%r = 0
  coefficients%r(351) = deltat
  coefficients%r(352) = gamma

  coefficients%func1(2)%p => sourcefunc

! create mesh

  meshgen_options%elshape = 6 ! 9-node quad
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%ox = ox
  meshgen_options%oy = oy
  meshgen_options%lx = lx
  meshgen_options%ly = ly

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! problem definition for the flow problem

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

  call problem_definition ( input_probdef, mesh, problem )

! problem definition for the convection-diffusion problem

  call create_input_probdef ( mesh, input_probdefq, nvec=2, nphysq=1 )

  input_probdefq%vec_elementdof(1)%a =   &
      reshape ( [1,0,1,0,1,0,1,0,0,    &  ! q
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar
                 [9,2] )

  input_probdefq%physq = [1]

  call define_essential ( mesh, input_probdefq, curve1=4 )
  call define_essential ( mesh, input_probdefq, curve1=1 )

  call problem_definition ( input_probdefq, mesh, problemq )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problemq, solq, solq_n, solq_nm1 )
  call create_sysvector ( problemq, rhsdq )

! oldvectors

  call create ( oldvectors, nprob=2, nsysvec=3 )
  oldvectors%p(1)%p => problem
  oldvectors%p(2)%p => problemq
  oldvectors%s(1)%p => sol
  oldvectors%s(2)%p => solq_n
  oldvectors%s(3)%p => solq_nm1

! initialize solution vectors

  sol%u = 0
  solq%u = 0

! fill velocity field

  call fill_sysvector ( mesh, problem, sol, node1=1, node2=mesh%nnodes, &
    physq=physqvel, value=1/sqrt(2._dp) )

! fill initial q-field

  call fill_sysvector ( mesh, problemq, solq, node1=1, node2=mesh%nnodes, &
    value=0._dp )

! fill q BCs

  call fill_sysvector ( mesh, problemq, solq, curve1=1, value=0._dp )
  call fill_sysvector ( mesh, problemq, solq, curve1=4, value=0._dp )

! create system matrix

  call create_sysmatrix_structure ( sysmatrixq, mesh, problemq )

  call create_sysmatrix_data ( sysmatrixq )

! postprocessing

  step = 0

  call postprocessing

  call copy ( solq, solq_n )

! start time stepping

  buildmatrix = .true.

  do step = 1, numtimesteps

    if ( step == 2 ) coefficients%i(401) = timeint2 ! second-order
    if ( step == 3 ) buildmatrix = .false.  ! matrix is constant

!   build diffusion matrix and right-hand side

    call build_system ( mesh, problemq, sysmatrixq, rhsdq, &
      elemsub=conv_diff_supg_elem, coefficients=coefficients, &
      oldvectors=oldvectors, buildmatrix=buildmatrix )

    call add_effect_of_essential_to_rhs ( problemq, sysmatrixq, solq, rhsdq )

    solver_options%real_storage=rs
    solver_options%integer_storage=is

    call solve_system_ma41 ( sysmatrixq, rhsdq, solq, lu, solver_options )

    if ( step == 1 ) call delete(lu) ! perform LU-decomposition at second step

    write(*,'(a,i0,a,es11.4,a,es15.8/)') &
      'step = ', step, ' time = ', step*deltat, &
      ' average sol = ', sum( abs(solq%u) ) / solq%n

!   postprocessing

    call postprocessing

    call copy ( solq_n, solq_nm1 )
    call copy ( solq, solq_n )

  end do


! delete all data including all allocated memory

  call delete ( problem, problemq )
  call delete ( input_probdef, input_probdefq )
  call delete ( mesh )
  call delete ( sol, solq, solq_n, solq_nm1 )
  call delete ( rhsdq )
  call delete ( sysmatrixq )
  call delete ( coefficients )
  call delete ( oldvectors )

contains

! subroutine to write vtks

  subroutine postprocessing

    type(oldvectors_t) :: oldvectors_post
    type(vector_t) :: q

    call create ( oldvectors_post, nsysvec=1 )
    oldvectors_post%s(1)%p => solq

!   derive q in all nodes

    call create_vector ( problemq, q, vec=2 )

    call derive_vector ( mesh, problemq, q, &
      elemsub=derive_q, coefficients=coefficients, &
      oldvectors=oldvectors_post )

    write(filename,'(a,i4.4,a)') 'solq', step, '.vtk'
    call write_scalar_vtk ( mesh, problemq, vector=q, &
      filename=filename, dataname='t' )

    write(filename,'(a,i4.4,a)') 'velo', step, '.vtk'
    call write_vector_vtk ( mesh, problem, filename=filename, &
      dataname='velocity', sysvector=sol, physq=physqvel )

    call delete ( oldvectors_post )
    call delete ( q )

  end subroutine postprocessing

end program convection_diffusion_supg1
