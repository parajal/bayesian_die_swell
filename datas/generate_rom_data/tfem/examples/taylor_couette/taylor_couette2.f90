! 2D steady Taylor-Couette flow simulation program
! for verificaiton of the routines computing persistence-of-straining

program taylor_couette2

  use tfem_m
  use gmsh_utils_m
  use math_defs_m
  use hsl_ma57_m
  use stokes_elements_m
  use io_utils_m

  implicit none

! constant parameters

  integer, parameter :: &
    uintpl = 8,       & ! Q2 velocities
    pintpl = 4,       & ! Q1 pressures
    gauss = 3,        & ! 3x3 integration of quad
    physqvel = 1,     & ! physical quantity nr of the velocities
    physqpress = 2      ! physical quantity nr of the pressures

! variables

  real(dp) :: &
    eta = 1._dp,      & ! viscosity
    Uo = 1._dp,       & ! angular velocity of the outer wall
    Ui = -0.5_dp,     & ! angular velocity of the inner wall
    rs_up = 1.0_dp,   & ! real_storage velocity-pressure LU (HSL)
    is_up = 1.0_dp      ! integer_storage velocity-pressure LU (HSL)

  real(dp) :: &
    theta = 0.0_dp     ! angle of Couette slice

  character(50) :: &
    filename = 'post.vtk', & ! output file name
    outdir = 'taylor_couette2_files'  ! output directory name

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(coefficients_t) :: coefficients
  type(oldvectors_t) :: oldvectors
  type(vector_t), target :: strain
  type(vector_t) :: persist, fclass
  type(solver_options_ma57_t) :: solver_options

! namelist for input of variables; read from standard input

  namelist /comppar/ &
    eta, Uo, Ui

! write initial message

  write(*,'(/a/)') 'Example program begins'

! read input file

  read ( unit=*, nml=comppar )

! generate mesh2.msh from mesh2.geo

  write(*, '(a/)') 'Execute gmsh to generate mesh'

  call execute_command_line ( 'gmsh mesh2.geo -save' )

! read mesh2.msh

  call read_mesh_gmsh ( mesh, filename='mesh2.msh', ndim=2, &
    physgeom=.true., sortphys=.true. )

! find angle between two curves

  call angle_curves ( mesh=mesh, curve1=1, curve2=2, theta=theta )

  write(*,'(/a,es14.4,a/)') 'Angle of the slice is ', theta, ' rad'

! match the numbering of the rotated curves

  call add_to_mesh ( mesh, matchingcurve=[2,1], replace=2, &
    rotation_angle=[theta], excludepoints=[1,2] )

  call fill_mesh_parts ( mesh )

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=450, ncoefr=400 )

  coefficients%i = 0
  coefficients%i(1:11) = &
     [ uintpl,   pintpl,     0,     0,        0, &
       physqvel, physqpress, 0,     0,    gauss, &
       gauss ]

  coefficients%r = 0
  coefficients%r(1) = eta

! define stokes problem, create sysvectors, subscripts, oldvectors...

  call define_problems_create_vectors

! solve stokes problem

  call build_solve_stokes

! compute kinematic parameters

  call derive_vector ( mesh, problem, strain, &
    elemsub=stokes_D_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, persist, &
    elemsub=stokes_P_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, fclass, &
    elemsub=stokes_P_scalar, &
    coefficients=coefficients, oldvectors=oldvectors )

! create a folder for output files

  call execute_command_line ( 'mkdir -p '//trim(outdir) )

! write files

  write(filename,'(a,a,a)') './', trim(outdir), '/post.vtk'

  call write_vector_vtk ( mesh, problem, physq=1, &
      filename=filename, sysvector=sol, dataname='vel', &
      append = .false. )

  call write_tensor_vtk ( mesh, problem, &
      filename=filename, vector=strain, dataname='D', &
      append = .true. )

  call write_scalar_vtk ( mesh, problem, &
      filename=filename, vector=persist, dataname='P', &
      append = .true. )

  call write_scalar_vtk ( mesh, problem, &
      filename=filename, vector=fclass, dataname='R', &
      append = .true. )

! delete all data including all allocated memory

  call delete_old_problems


contains

  subroutine define_problems_create_vectors

!   problem definition of velocity/pressure

    call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

    input_probdef%vec_elementdof(1)%a = &
          reshape ( [ 2,2,2,2,2,2,2,2,2,  & ! velocity
                      1,0,1,0,1,0,1,0,0,  & ! pressure
                      1,1,1,1,1,1,1,1,1,  & ! P matrix or scalar
                      3,3,3,3,3,3,3,3,3], & ! D matrix
                      [9,4] )

    input_probdef%physq = [physqvel,physqpress]

    input_probdef%probnr = 1

!   essential boundary condition

    call define_essential ( mesh, input_probdef, curve1=3, curve2=4, &
      physq=physqvel )

    call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

!   periodic velocities

    call define_constraint ( mesh, input_probdef, &
      physq=physqvel, curve1=1, curve2=2, discretization='collocation', &
      excludecurves=[3,4] )

!   define transformation

    call define_transformation ( mesh, input_probdef, curve=1, &
      physq=physqvel, normalvector=1, excludecurves=[3,4] )

    call define_transformation ( mesh, input_probdef, curve=2, &
      physq=physqvel, normalvector=1, excludecurves=[3,4] )

!   define problem

    call problem_definition ( input_probdef, mesh, problem )

!   build transformation matrix

    call build_transformation_matrix ( mesh, problem )

!   create system vectors for the up problem

    call create_sysvector ( problem, sol, rhsd )

!   create system matrix for the stokes problem

    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric=.true. )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )
    call create_sysmatrix_data ( sysmatrix )

!   create the structure oldvectors and vector

    call create_oldvectors ( oldvectors, nsysvec=1, nvec=3 )

    call create_vector ( problem, strain, vec=4 )
    call create_vector ( problem, persist, vec=3 )
    call create_vector ( problem, fclass, vec=3 )

    oldvectors%s(1)%p => sol
    oldvectors%v(3)%p => strain

  end subroutine define_problems_create_vectors


! delete old problem

  subroutine delete_old_problems

    call delete ( problem )
    call delete ( sysmatrix )
    call delete ( input_probdef )
    call delete ( coefficients )
    call delete ( mesh )
    call delete ( sol, rhsd )
    call delete ( oldvectors )
    call delete ( persist, strain, fclass )

  end subroutine delete_old_problems


! build and solve the stokes equation

  subroutine build_solve_stokes

!   fill solution vector with essential boundary conditions

    sol%u = 0

!   inner cylinder

    call fill_sysvector ( mesh, problem, sol, curve1=3, &
      physq=physqvel, vfunc=vfunc, vfuncnr=1 )

!   outer cylinder

    call fill_sysvector ( mesh, problem, sol, curve1=4, &
      physq=physqvel, vfunc=vfunc, vfuncnr=2 )

!   point pressure

    call fill_sysvector ( mesh, problem, sol, point=1, &
      physq=physqpress, value=0._dp )

!   stokes velocity/pressure

    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients )

!   periodical condition on velocities

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients, transform=.false. )

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

!   solve sysmatrix

    solver_options%real_storage=rs_up
    solver_options%integer_storage=is_up

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options )

    call transform_to_global ( problem, sol )

  end subroutine build_solve_stokes


! function for velocity on wall

  function vfunc ( n, nr, x )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    real(dp) :: normvec(2)
    real(dp) :: r

!   compute distance from origin

    r = norm(x)

!   get vector perpendicular to position vector

    normvec = [-x(2),x(1)]

!   normalize

    normvec = normvec / norm(normvec)

    select case(nr)
      case(1)
        vfunc = Ui * r * normvec
      case(2)
        vfunc = Uo * r * normvec
      case default
        write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
        stop
    end select

  end function vfunc


! find angle between two curves
! NOTE: curves should be straight lines pointing to origin

  subroutine angle_curves ( mesh, curve1, curve2, theta )

    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: curve1, curve2
    real(dp), intent(out) :: theta

    real(dp) :: p1(2), p2(2), x0(2)
    real(dp) :: c1vec(2), c2vec(2)
    real(dp) :: cross

!   set origin coordinate

    x0 = 0._dp

!   coordinate of first node on the curves

    p1(:) = mesh%coor(mesh%curves(curve1)%nodes(1),:)
    p2(:) = mesh%coor(mesh%curves(curve2)%nodes(1),:)

!   vectors pointing from origin to p1 and p2

    c1vec(:) = p1(:) - x0(:)
    c2vec(:) = p2(:) - x0(:)

!   find direction of angle

    cross = c1vec(1) * c2vec(2) - c1vec(2) * c2vec(1)

!   compute angle theta

    theta = &
      sign ( 1._dp, cross ) * &
      acos ( dot_product(c1vec, c2vec) / &
           ( norm2(c1vec) * norm2(c2vec) ) )

  end subroutine angle_curves

end program taylor_couette2
