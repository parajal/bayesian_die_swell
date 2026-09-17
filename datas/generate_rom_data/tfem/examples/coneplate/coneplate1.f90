! Cone-plate problem using the Stokes equations.
! (2D) domain with 3D velocities (axisymmetric).
! Cylindrical coordinates with circumferential velocity included (swirl).

program coneplate1

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 6,   & ! P2 velocities
    pintpl = 2,   & ! P1 pressures
    physqv = 1,   & ! physical quantity nr of the velocities
    physqp = 2,   & ! physical quantity nr of the pressures
    gauss = 6,    & ! 6-point Gauss integration of triangles
    gaussb = 3      ! 3-point integration of boundary elements

  real(dp), parameter :: &
    Ro = 1.0_dp,        & ! outer radius (of fluid contact point on the plate)
    eta = 1._dp,        & ! viscosity
    theta_cp = 0.1_dp,  & ! Cone angle
    dx_Ro = 0.01_dp,    & ! Element size at the fluid outside surface
    dx_0 = 0.001_dp,    & ! Element size at the center (r=0)
    omega = 0.1_dp        ! Angular velocity of the cone

  integer :: vertices(3) = [1,3,5]

  integer :: i

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, duthetadr, gammadot, Dtensor, Ltensor
  type(vector_t) :: tau_rtheta, tautensor
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

  character(len=30) :: filename


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl, pintpl, 0, 0, 0,      &
      physqv, physqp, 0, 0, gauss,  &
      gaussb ]
  coefficients%i(12:) = 0

  coefficients%i(23) = 1  ! cylindrical coordinate system (axisymmetrical)
  coefficients%i(67) = 1  ! 3D velocity (swirl)

  coefficients%r = 0
  coefficients%r(1) = eta

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  call create_mesh

  call fill_mesh_parts ( mesh )

! write mesh for plotting

  call write_mesh_vtk ( mesh, 'mesh.vtk' )

  do i = 1, mesh%ncurves
    write(filename,'(a,i4.4,a)') 'curve_',i,'.vtk'
    call write_geometry_vtk ( mesh, curve=i, filename=filename )
  end do

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=5, nphysq=2 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 0
  input_probdef%vec_elementdof(1)%a(vertices,2) = 1  ! pressure
  input_probdef%vec_elementdof(1)%a(:,3) = 1  ! scalar, such as vorticity
  input_probdef%vec_elementdof(1)%a(:,4) = 6  ! symmetric tensor
  input_probdef%vec_elementdof(1)%a(:,5) = 9  ! unsymmetric tensor

  input_probdef%physq = [1,2]

! Cone
  call define_essential ( mesh, input_probdef, curve1=2, physq=physqv )

! Plate
  call define_essential ( mesh, input_probdef, curve1=1, physq=physqv )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, rhsd )

! fill solution vector with essential boundary conditions

  sol%u = 0
  call fill_sysvector ( mesh, problem, sol, &
    curve1=2, physq=physqv, degfd=3, func=func, funcnr=1 )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, &
    elemsub=stokes_elem, coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, duthetadr, vec=3 )
  call create_vector ( problem, gammadot, vec=3 )
  call create_vector ( problem, Dtensor, vec=4 )
  call create_vector ( problem, Ltensor, vec=5 )
  call create_vector ( problem, tau_rtheta, vec=3 )
  call create_vector ( problem, tautensor, vec=4 )

  call extract_physvector ( mesh, problem, sol, velocity )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=8

  call derive_vector ( mesh, problem, duthetadr, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=11

  call derive_vector ( mesh, problem, gammadot, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Dtensor, elemsub=stokes_D_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, Ltensor, elemsub=stokes_gradu_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, tau_rtheta, elemsub=stokes_stress, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, tautensor, elemsub=stokes_stress_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a vtk file for further post-processing

  call write_scalar_vtk ( mesh, problem, vector=pressure, dataname='pressure', &
    filename='cp1.vtk' )

  call write_vector_vtk ( mesh, problem, filename='cp1.vtk', &
    dataname='velocity_vector', vector=velocity, append=.true., &
    assume3D=.true. )

  call write_vector_vtk ( mesh, problem, filename='cp1.vtk', &
    dataname='velocity_vector_planar', vector=velocity, append=.true., &
    degfd=[1,2] )

  call write_scalar_vtk ( mesh, problem, vector=velocity, degfd=3, &
    filename='cp1.vtk', dataname='velocity_theta', append=.true. )

  call write_scalar_vtk ( mesh, problem, vector=gammadot, &
    filename='cp1.vtk', dataname='gammadot', append=.true. )

  call write_tensor_vtk ( mesh, problem, filename='cp1.vtk', &
    dataname='D', vector=Dtensor, append=.true., assume3D=.true. )

  call write_tensor_vtk ( mesh, problem, filename='cp1.vtk', &
    dataname='L', vector=Ltensor, append=.true., assume3D=.true., &
    symmetric=.false. )

  call write_tensor_vtk ( mesh, problem, filename='cp1.vtk', &
    dataname='tau', vector=tautensor, append=.true., assume3D=.true. )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, duthetadr, gammadot, Dtensor, Ltensor )
  call delete ( tau_rtheta, tautensor )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )


contains


! create the mesh for the cone-plate geometry using gmsh

  subroutine create_mesh

    use math_defs_m

    open ( unit=25, file='ConePlate2d.geo' )

    write ( 25, '(1X,A,F18.14,A)' ) 'z[1] = ', Ro*cos(pi/2-theta_cp), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'r[1] = ', Ro*sin(pi/2-theta_cp), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'theta_cp = ', theta_cp, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_Ro = ', dx_Ro, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_0 = ', dx_0, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'Ro = ', Ro, ';'
    write ( 25, '(/1x,a)' ) 'Include "ConePlate2d.igo";'

    close ( 25 )

    call execute_command_line ( 'gmsh -2 -order 2 -algo front2d &
       & -o ConePlate2d.msh ConePlate2d.geo > outputmesh.out' )

!   read mesh generated by gmsh
    call read_mesh_gmsh ( mesh, filename='ConePlate2d.msh', ndim=2, &
      physgeom=.true., sortphys=.true. )

  end subroutine create_mesh


! function for velocity on cone

  function func ( nr, x )

    integer, intent(in) ::  nr
    real(dp), intent(in), dimension(:) :: x
    real(dp) :: func

    select case(nr)
      case(1)
        func = omega*x(2)
      case default
        write(*,'(/a,i0/)') 'Error func: wrong function number: ', nr
        stop
    end select

  end function func

end program coneplate1
