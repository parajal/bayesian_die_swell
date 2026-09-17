! 3D Stokes problem with an open boundary condition.
! zero integral pressure imposed on the obc

program open_boundary3

  use tfem_m
  use hsl_ma41_m
  use stokes_elements_m
  use velocity_functions_m
  use io_utils_m
  use figplot_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
!    pintpl = 2,         & ! P1 pressures
    ointpl = 8,         & ! Q2 shape of object elements
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=10,              & ! number of elements in x
    ny=10,              & ! number of elements in y
    nz=10                 ! number of elements in z

  real(dp), parameter :: &
    eta = 1._dp     ! viscosity

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, dudy, divergence
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

  integer  :: presnod(8) = [1,3,9,7,19,21,27,25]


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=300, ncoefr=250 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gauss ]
  coefficients%i(12:) = 0

  coefficients%i(35) = ointpl

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  meshgen_options%elshape = 14
  meshgen_options%regionshape = 4
  meshgen_options%nx = nx
  meshgen_options%ny = ny
  meshgen_options%nz = nz

  call hexahedron ( mesh, meshgen_options )

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

  call add_to_mesh ( mesh, object='surface', objectmesh=mesh, &
    objectsurface=3, topology=.true., intrule=3, nsubint=1 )

  call add_to_mesh ( mesh, point=[1.0_dp,0.5_dp,0.5_dp] )

  call fill_mesh_parts ( mesh )

  call plot_points_curves ( plot_options, mesh, 'curves.fig' )
  call plot_mesh ( plot_options, mesh, 'mesh.fig', surfaces=[3,4,6] )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = mesh%ndim ! velocity

  if ( pintpl == 2 ) then                            ! P1 pressure
    input_probdef%vec_elementdof(1)%a(:,2) = 0
    input_probdef%vec_elementdof(1)%a(14,2) = mesh%ndim+1
  else if ( pintpl == 4 )  then                      ! Q1 pressure
    input_probdef%vec_elementdof(1)%a(:,2) = 0
    input_probdef%vec_elementdof(1)%a(presnod,2) = 1
  end if

  input_probdef%vec_elementdof(1)%a(:,3) = 1

  input_probdef%physq = [1,2]

! essential BCs

  call define_essential ( mesh, input_probdef, surface1=1, surface2=2, physq=1 )
  call define_essential ( mesh, input_probdef, surface1=4, surface2=6, physq=1 )

!  if ( pintpl == 2) then
!    call define_essential ( mesh, input_probdef, element=nx, elnode=14, physq=2 )
!  else if ( pintpl == 4) then
!    call define_essential ( mesh, input_probdef, point=2, physq=2 )
!  end if

! constraint for integral of the pressure on the obc

  if ( pintpl == 2) then                     ! P1 pressure
    call define_constraint ( mesh, input_probdef, &
      physq=2, object=1, nglobalc=1 )
  else if ( pintpl == 4 ) then               ! Q1 pressure
    call define_constraint ( mesh, input_probdef, &
      physq=2, surface1=3, nglobalc=1 )
  end if

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0._dp

  call fill_sysvector ( mesh, problem, sol, surface1=5, degfd=1, &
    physq=1, func=vel_func, funcnr=5 )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call build_system(mesh, problem, sysmatrix, rhsd, &
    elemsub1=stokes_open_boundary, coefficients=coefficients, &
    oldvectors=oldvectors, object=1, addmatvec=.true. )

  if ( pintpl == 2) then                     ! P1 pressure
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_zero_pressure_int_object, &
      addmatvec=.true., coefficients=coefficients )
  else if ( pintpl == 4 ) then               ! Q1 pressure
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_zero_pressure_int_surface, &
      addmatvec=.true., coefficients=coefficients )
  end if

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma41 ( sysmatrix, rhsd, sol )

  !call get_sysvector_constraint ( mesh, problem, sysvector=sol, constraint=1, &
    !u=u )

  !print *, 'u=', u

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, dudy, vec=3 )
  call create_vector ( problem, divergence, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=2

  call derive_vector ( mesh, problem, dudy, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=10

  call derive_vector ( mesh, problem, divergence, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a fig file for plotting

  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity, surfaces=[3,5] )

  call plot_color_fill ( plot_options, mesh, problem, 'velocity_color.fig', &
    vector=velocity, degfd=1, surfaces=[3,5]  )

  call plot_color_fill ( plot_options, mesh, problem, 'pressure_color.fig', &
    vector=pressure, surfaces=[3,4,6]  )

  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure, surfaces=[3]  )

  call plot_color_fill ( plot_options, mesh, problem, 'dudy_color.fig', &
    vector=dudy,  surfaces=[3,4,6] )

  call plot_color_contour ( plot_options, mesh, problem, &
    'dudy_contour.fig', vector=dudy,  surfaces=[3] )

  call plot_color_fill ( plot_options, mesh, problem, 'divergence_color.fig', &
    vector=divergence,  surfaces=[3,4,6] )

  call plot_color_contour ( plot_options, mesh, problem, &
    'divergence_contour.fig', vector=divergence, surfaces=[3] )

! write binary file for reading by streamfunction

  open(unit=10,form='unformatted',file='velocity_bin.out')

  write(10) velocity%u

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, dudy, divergence )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )

end program open_boundary3
