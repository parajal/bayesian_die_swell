! unsteady convection problem with two coupled equations
! solution of the first equation is used in the RHS of the second
! manufactured solutions for both equations
! constant velocity field
! SUPG for stabilisation
! second-order time integration
! quadratic interpolation of q
!
! NOTE: this example uses SUPG combined with quadratic interpolation of q, which
! appears to work well for pure convection problems. If the diffusion term is
! non-zero, it is advised to use linear interpolation of q when applying SUPG.

program convection_diffusion_supg2

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
    qintpl = 8,         & ! q interpolation
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    timeint1 = 1,       & ! first-order implicit Euler
    timeint2 = 2,       & ! second-order implicit Gear
    gauss = 5,          & ! order of integration of quads
    coorsys = 0,        & ! planar coordinate system
    neq = 2,            & ! number of convection-diffusion equations
    nx = 4,             & ! number of elements in x-direction
    ny = 4,             & ! number of elements in y-direction
    ibeta = 1,          & ! method to compute beta in SUPG
    numtimesteps=100      ! number of time steps

  real(dp), parameter :: &
    gamma = 1._dp,       & ! parameter in the convection-diffusion equation
    deltat = 1e-2_dp,    & ! time step
    ox = 0._dp,          & ! x-coordinate of origin of the mesh
    oy = 0._dp,          & ! y-coordinate of origin of the mesh
    lx = 1._dp,          & ! length of the domain in x-direction
    ly = 1._dp,          & ! length of the domain in y-direction
    rs = 1.8_dp,         & ! real_storage velocity-pressure LU (HSL)
    is = 1.8_dp            ! integer_storage velocity-pressure LU (HSL)

! variables

  logical :: buildmatrix

  integer :: step, i

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef, input_probdefq
  type(problem_t), target :: problem, problemq
  type(sysmatrix_t) :: sysmatrixq
  type(sysvector_t), target :: sol
  type(sysvector_t), target, dimension(neq) :: solq, solq_n, solq_nm1, rhsdq
  type(sysvector_t), target, dimension(neq) :: solexact
  type(elvector_t), target, dimension(neq) :: f
  type(coefficients_t) :: coefficients
  type(oldvectors_t), dimension(neq) :: oldvectors
  type(lu_ma41_t) :: lu
  type(solver_options_ma41_t) :: solver_options

  character(len=199) :: filename, dataname

  real(dp) :: time

  time = 0 ! initial time
  ltime = time ! pass to module

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
  coefficients%i(407) = 3 ! funcnr for f-term
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
      reshape ( [1,1,1,1,1,1,1,1,1,    &  ! q
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar
                 [9,2] )

  input_probdefq%physq = [1]

  call define_essential ( mesh, input_probdefq, curve1=4 )
  call define_essential ( mesh, input_probdefq, curve1=1 )

  call problem_definition ( input_probdefq, mesh, problemq )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_msysvector ( problemq, solq, solq_n, solq_nm1, solexact )
  call create_msysvector ( problemq, rhsdq )

! create vectors / oldvectors for each convection-diffusion equation

  do i = 1, neq

!   create vector defined per element for scalar values of right-hand
    call create_elvector ( mesh, f(i), nreal1d=1 )

!   oldvectors
    call create ( oldvectors(i), nprob=2, nsysvec=4, nelvec=2 )
    oldvectors(i)%p(1)%p => problem
    oldvectors(i)%p(2)%p => problemq
    oldvectors(i)%s(1)%p => sol
    oldvectors(i)%s(2)%p => solq_n(i)
    oldvectors(i)%s(3)%p => solq_nm1(i)
    oldvectors(i)%e(2)%p => f(i)

!   fill initial q-field
    call fill_sysvector ( mesh, problemq, solq(i), node1=1, node2=mesh%nnodes, &
      func=sourcefunc, funcnr=2*i )

  end do

  oldvectors(2)%s(4)%p => solq(1) ! use the updated q1 solution in q2 problem
                                  ! (applies to the fill_f_elem routine)

! fill velocity field

  call fill_sysvector ( mesh, problem, sol, node1=1, node2=mesh%nnodes, &
    physq=physqvel, value=1._dp )

! create system matrix

  call create_sysmatrix_structure ( sysmatrixq, mesh, problemq )

  call create_sysmatrix_data ( sysmatrixq )

! postprocessing

  step = 0

  call postprocessing

  call copy ( solq, solq_n )

! start time stepping

  do step = 1, numtimesteps

    if ( step == 2 ) coefficients%i(401) = timeint2 ! second-order

    time = time + deltat
    ltime = time

!   loop over the equations

    if ( step <= 2 ) buildmatrix = .true. ! after timestep 2 matrix is constant

    do i = 1, neq

!     fill q BCs (change in time!)
      call fill_sysvector ( mesh, problemq, solq(i), curve1=1, &
        func=sourcefunc, funcnr=i*2 )
      call fill_sysvector ( mesh, problemq, solq(i), curve1=4, &
        func=sourcefunc, funcnr=i*2 )

!     set correct variables for the f-term in the current equation
      if ( i == 1 ) then
        coefficients%i(406) = 2 ! f-term given by function
        coefficients%i(407) = 3 ! funcnr for f-term
      else
        coefficients%i(406) = 4 ! f-term given in the integration points

!       fill the term in the integration points (based on the updated q1!)
        call loop_over_elements ( mesh, problemq, elemsub=fill_f_elem, &
          coefficients=coefficients, oldvectors=oldvectors(i) )

      end if

!     build the system and RHS (constant flow: after step 2 system is constant)
      call build_system ( mesh, problemq, sysmatrixq, rhsdq(i), &
        elemsub=conv_diff_supg_elem, coefficients=coefficients, &
        oldvectors=oldvectors(i), buildmatrix=buildmatrix )

      call add_effect_of_essential_to_rhs ( problemq, sysmatrixq, solq(i), &
        rhsdq(i) )

      solver_options%real_storage=rs
      solver_options%integer_storage=is

      call solve_system_ma41 ( sysmatrixq, rhsdq(i), solq(i), lu, &
        solver_options )

      buildmatrix=.false. ! after first pass matrix is constant

    end do

    if ( step == 1 ) call delete(lu) ! perform LU-decomposition at second step

!   postprocessing

    call postprocessing

    write(*,'(/a,i0,a,es11.4,a,es15.8,a,es15.8)') &
      'step = ', step, ' time = ', step*deltat

    do i = 1, neq

      call fill_sysvector ( mesh, problemq, solexact(i), node1=1, &
        node2=mesh%nnodes, func=sourcefunc, funcnr=2*i )

      write(*,'(a,i1,a,es15.8)') &
      ' error q' , i, ' = ', maxval(abs(solq(i)%u-solexact(i)%u))

    end do

!   copy old solutions

    call copy ( solq_n, solq_nm1 )
    call copy ( solq, solq_n )

  end do


! delete all data including all allocated memory

  call delete ( problem, problemq )
  call delete ( input_probdef, input_probdefq )
  call delete ( mesh )
  call delete ( sol )
  call delete ( solq, solq_n, solq_nm1 )
  call delete ( rhsdq )
  call delete ( sysmatrixq )
  call delete ( coefficients )
  do i = 1, neq
    call delete ( oldvectors(i) )
  end do

contains

! subroutine to write vtks

  subroutine postprocessing

    logical :: append
    integer :: i

    write(filename,'(a,i4.4,a)') 'solq', step, '.vtk'
    append=.false.
    do i = 1, neq
      write(dataname,'(a,i1)') 'q', i
      call write_scalar_vtk ( mesh, problemq, sysvector=solq(i), &
        filename=filename, dataname=dataname, append=append )
      append=.true.
    end do

    write(filename,'(a,i4.4,a)') 'velo', step, '.vtk'
    call write_vector_vtk ( mesh, problem, filename=filename, &
      dataname='velocity', sysvector=sol, physq=physqvel )

  end subroutine postprocessing

end program convection_diffusion_supg2
