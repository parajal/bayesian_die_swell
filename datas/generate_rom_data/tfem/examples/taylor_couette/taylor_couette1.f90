! Stokes problem:
! Taylor-Couette flow with a spherical rigid particle subjected to a force.
! Only a segment of angle theta of the flow is considered and periodical
! boundary conditions at the cuts are assumed.
! The velocities are imposed on the inner and outer cylinder.
! At the upper and lower surfaces the vertical velocity is set to zero and
! the tangential tractions are assumed to be zero.

program taylor_couette1

  use tfem_m
  use math_defs_m
  use hsl_ma57_m
  use stokes_elements_m
  use io_utils_m
  use subs_couette_m
  use limits_m
  use timer_m

  implicit none


! constants

  integer, parameter :: &
    uintpl = 6,       & ! P2 velocities
    pintpl = 2,       & ! P1 pressures
    gauss = 24,       & ! 24-point Gauss integration of tets
    gaussb = 3,       & ! 3-point integration of boundary elements
    physqvel = 1,     & ! physical quantity nr of the velocities
    physqpress = 2      ! physical quantity nr of the pressures

! variables

  real(dp) :: &
    mu = 1._dp,      & ! viscosity
    Uo = 1._dp,      & ! magnitude of the outer wall velocity
    Ui = -0.5_dp,    & ! magnitude of the inner wall velocity
    theta = 45._dp * 2 * pi / 360, & ! angle of Couette slice
    Ro = 2._dp,         & ! outer radius of Couette device
    Ri = 1._dp,         & ! inner radius of Couette device
    height = 2._dp,     & ! height of the Couette slice
    radius = 0.1_dp,    & ! radius of the particle
    yposp = 1.5_dp,    &  ! y-position of the particle
    Fz = -2._dp,        & ! force on the particle in z-direction
    rs_gup = 2.4_dp,    & ! real_storage gradient-velocity-pressure LU (HSL)
    is_gup = 2.5_dp       ! integer_storage gradient-velocity-pressure LU (HSL)

  integer :: &
    n = 6,          & ! number of elements along gap
    np = 40,        & ! number of elements along equator particle
    n_ref = 10,     & ! number of elements on inner cylinder
    n_th = 5          ! number of elements in angle theta

! definitions for up problem

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(coefficients_t) :: coefficients
  type(subscript_t) :: velx, vely, velz
  type(solver_options_ma57_t) :: solver_options_ma57


! variables

  real(dp), dimension(3) :: pp, normal
  real(dp), dimension(6) :: up
  real(dp), dimension(3,3) :: Rmat

  real(dp) :: dx_box, dx_part, dx_refined


! set some parameters

! rotation matrix (counterclockwise, fixed coordinate system)

  Rmat=reshape([ cos(theta),sin(theta),0._dp, &
                -sin(theta),cos(theta),0._dp, &
                         0._dp, 0._dp, 1._dp],[3,3])


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i = 0
  coefficients%i(1:11) = &
     [ uintpl,   pintpl,     0,     0,        0, &
       physqvel, physqpress, 0,     0,    gauss, &
       gaussb ]

  coefficients%r = 0
  coefficients%r(1) = mu

! particle info

  allocate( xp(3), force(3) )

  rp = radius        ! particle radius
  xp(1) = 0          ! x-position of the particle
  xp(2) = yposp      ! y-position of the particle
  xp(3) = 0          ! z-position of the particle

  set_force = .true.
  force = [ 0._dp, 0._dp, Fz ]


! Generate and read mesh

  dx_box = (Ro-Ri)/n   ! element spacing on the external boundaries
  dx_part = 2*pi*rp/np ! element spacing on the particle boundary
  dx_refined = theta*Ri/n_ref ! element spacing on inner cylinder

  pp = xp

  call tic

  call generate_read_mesh

  call toc ( 'generate and read mesh')

! write vtks of the geometries
! NOTE: to visualize points and vectors in Paraview, use the "Glyph" filter

  call write_geometries_vtk ( mesh, write_normals=.true. )

  call toc ( 'write_geometries_vtk' )

! get the (outwardly directed) normal on surface 4

  normal = [-Ro*sin(theta/2),Ro*cos(theta/2),0._dp] - &
           [-Ri*sin(theta/2),Ri*cos(theta/2),0._dp]
  normal = [-normal(2),normal(1),0._dp]
  normal = normal / sqrt(dot_product(normal,normal))


! define the up problem, create (sys)vectors, subscripts, oldvectors
! and initialize sysvectors

  call define_problems_create_vectors

  call toc ( 'define_problems_create_vectors' )

! solve up problem

  call build_solve_up

! extract info on particle

  call get_sysvector_constraint ( mesh, problem, sol, constraint=1, &
    addunknowns=.true., u=up )

  write(*,*) ' up =     ', up(1:3)
  write(*,*) ' omegap = ', up(4:6)

! write VTK files

  call postprocessing ( mesh, sol )

  call toc ( 'postprocessing' )


! delete all data including all allocated memory

  call delete_old_problems

  deallocate( xp, force )

  call toc ( 'delete' )

contains


! generate and read mesh

  subroutine generate_read_mesh

    real(dp) :: rcoor

!   radial coordinate of the particle
    rcoor = sqrt(dot_product(pp(1:2),pp(1:2)))

    open ( unit=25, file='mesh.geo' )

    write ( 25, '(1X,A,F18.14,A)' ) 'x[1] = ', -Ro*sin(theta/2), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'y[1] = ',  Ro*cos(theta/2), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'x[2] = ', -Ri*sin(theta/2), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'y[2] = ',  Ri*cos(theta/2), ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'z[1] = ', -height/2, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'z[2] = ',  height/2, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'theta = ', theta, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_box = ', dx_box, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_refined = ', dx_refined, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx_part = ', dx_part, ';'
    write ( 25, '(1X,A,I0,A)' ) 'n_th = ', n_th, ';'
    write ( 25, '(1X,A,I0,A)' ) 'nobj = ', 1, ';'

    write ( 25, '(1X,A,I0,A,F18.14,A)' ) 'xp[', 1, '] = ', 0._dp, ';'
    write ( 25, '(1X,A,I0,A,F18.14,A)' ) 'yp[', 1, '] = ', rcoor, ';'
    write ( 25, '(1X,A,I0,A,F18.14,A)' ) 'zp[', 1, '] = ', 0._dp, ';'
    write ( 25, '(1X,A,I0,A,F18.14,A)' ) 'rp[', 1, '] = ', rp, ';'

    write ( 25, '(/1x,a)' ) 'Include "particle_in_couette_part3d.igo";'

    close ( 25 )

    call execute_command_line ( 'gmsh -3 -order 2 -algo front2d -o mesh.msh &
                  &mesh.geo > outputmesh.out' )

!   read mesh generated by gmsh
    call read_mesh_gmsh ( mesh, filename='mesh.msh', ndim=3, &
      physgeom=.true., sortphys=.true. )

!   match the numbering of the rotated surfaces
    call add_to_mesh ( mesh, matchingsurface=[4,2], replace=4, &
      rotation_angle=[0._dp,0._dp,theta] )

!   match the numbering of the rotated curves
    call add_to_mesh ( mesh, matchingcurve=[2,3], replace=2, &
      rotation_angle=[0._dp,0._dp,theta] )
    call add_to_mesh ( mesh, matchingcurve=[4,5], replace=4, &
      rotation_angle=[0._dp,0._dp,theta] )

    call fill_mesh_parts ( mesh )

    call write_mesh_vtk ( mesh, filename='mesh.vtk' )

  end subroutine generate_read_mesh


! start the problems

  subroutine define_problems_create_vectors

!   problem definition of velocity/pressure

    call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

    input_probdef%vec_elementdof(1)%a = &
          reshape ( [ 3,3,3,3,3,3,3,3,3,3,        & ! velocity
                      1,0,1,0,1,0,0,0,0,1,        & ! pressure
                      1,1,1,1,1,1,1,1,1,1,        & ! scalar
                      6,6,6,6,6,6,6,6,6,6 ],      & ! tensor
                      [10,4] )

    input_probdef%physq = [physqvel,physqpress]

    input_probdef%probnr = 1

!   cell boundaries
    call define_essential ( mesh, input_probdef, surface1=1, &
      physq=physqvel )
    call define_essential ( mesh, input_probdef, surface1=3, &
      physq=physqvel )
    call define_essential ( mesh, input_probdef, surface1=5, &
      physq=physqvel, degfd=[0,0,1] )
    call define_essential ( mesh, input_probdef, surface1=6, &
      physq=physqvel, degfd=[0,0,1] )

    call define_essential ( mesh, input_probdef, point=1, physq=physqpress )

!   constraint on particle boundary (freely floating)
    call define_constraint ( mesh, input_probdef, surface1=7, &
      physq=physqvel, discretization='collocation', nodedof=3, naddunknowns=6 )

!   periodic velocities
    call define_constraint ( mesh, input_probdef, &
      physq=physqvel, surface1=4, surface2=2, discretization='collocation', &
      excludesurfaces=[1,3,5,6] )
    call define_constraint ( mesh, input_probdef, &
      physq=physqvel, curve1=2, curve2=3, discretization='collocation', &
      excludesurfaces=[1,3], nodedof=2 )
    call define_constraint ( mesh, input_probdef, &
      physq=physqvel, curve1=4, curve2=5, discretization='collocation', &
      excludesurfaces=[1,3], nodedof=2 )

    call problem_definition ( input_probdef, mesh, problem )


!   create system vectors for the up problem
    call create_sysvector ( problem, sol, rhsd )

!   create subscripts for the velocity
    call create_subscript ( mesh, problem, velx, physqarr=[physqvel], &
      degfd=1, fillnodes=.true. )
    call create_subscript ( mesh, problem, vely, physqarr=[physqvel], &
      degfd=2, fillnodes=.true. )
    call create_subscript ( mesh, problem, velz, physqarr=[physqvel], &
      degfd=3, fillnodes=.true. )

!   create system matrix for the upg problem
    call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
      symmetric=.true. )
    call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
    call finalize_sysmatrix_structure ( sysmatrix )
    call create_sysmatrix_data ( sysmatrix )

  end subroutine define_problems_create_vectors


! delete old problem

  subroutine delete_old_problems

    call delete ( problem )
    call delete ( sysmatrix )
    call delete ( input_probdef )
    call delete ( mesh )
    call delete ( sol, rhsd )

  end subroutine delete_old_problems


! build and solve the up system

  subroutine build_solve_up

!   fill solution vector with essential boundary conditions
    sol%u = 0

!   boundary conditions for the velocity
!   inner cylinder
    call fill_sysvector ( mesh, problem, sol, surface1=1, &
      physq=physqvel, vfunc=vfunc, vfuncnr=1 )
!   outer cylinder
    call fill_sysvector ( mesh, problem, sol, surface1=3, &
      physq=physqvel, vfunc=vfunc, vfuncnr=2 )

!   point pressure
    call fill_sysvector ( mesh, problem, sol, point=1, &
      physq=physqpress, value=0._dp )

!   stokes velocity/pressure
    call build_system ( mesh, problem, sysmatrix, rhsd, &
      elemsub=stokes_elem, coefficients=coefficients )

!   constraints on surfaces
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=elementc_3D, addmatvec=.true. )

!   periodical condition on velocities
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=constr_node_rotated_vec, addmatvec=.true., &
      coefficients=coefficients )
    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=3, constraint2=4, elemsub=constr_node_rotated_vec2, &
      addmatvec=.true., coefficients=coefficients )

    call check ( sysmatrix )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    call toc ( ' build ' )

!   solve up system
    solver_options_ma57%real_storage=rs_gup
    solver_options_ma57%integer_storage=is_gup

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options=solver_options_ma57 )

    call toc ( ' solve MA57 ' )

  end subroutine build_solve_up


! write the data to .vtk files

  subroutine postprocessing ( mesh, sol )

    type(sysvector_t), intent(in), target :: sol
    type(mesh_t), intent(in) :: mesh
    type(oldvectors_t) :: oldvectors_dve
    type(vector_t) :: pressure, eff_shear

    call create_oldvectors ( oldvectors_dve, nsysvec=1 )
    oldvectors_dve%s(1)%p => sol

    call create_vector ( problem, pressure, vec=3 )

!   derive the pressure in all nodes
    call derive_vector ( mesh, problem, pressure, &
      elemsub=stokes_pressure, coefficients=coefficients, &
      oldvectors=oldvectors_dve )

    call create_vector ( problem, eff_shear, vec=3 )

    coefficients%i(13) = 11

!   derive the effective shear rate in all nodes
    call derive_vector ( mesh, problem, eff_shear, &
      elemsub=stokes_deriv, coefficients=coefficients, &
      oldvectors=oldvectors_dve )

    call write_scalar_vtk ( mesh, problem, vector=pressure, &
      dataname='pressure',  filename='flow.vtk' )

    call write_vector_vtk ( mesh, problem, filename='flow.vtk', &
      dataname='velocity', sysvector=sol, physq=physqvel, &
      append=.true. )

    call write_scalar_vtk ( mesh, problem, vector=eff_shear, &
      dataname='eff_shear',  filename='flow.vtk', append=.true. )

   call delete(pressure)
   call delete(oldvectors_dve)

  end subroutine postprocessing


! function for velocity on wall

  function vfunc ( n, nr, x )

    integer, intent(in) :: n, nr
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(n) :: vfunc

    real(dp) :: normvec(3)

!   get vector perpendicular to position vector
    normvec = [x(2),-x(1),0._dp]

!   normalize
    normvec = normvec / sqrt(dot_product(normvec,normvec))

    select case(nr)
      case(1)
        vfunc = Ui*normvec
      case(2)
        vfunc = Uo*normvec
      case default
        write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
        stop
    end select

  end function vfunc


! Element for the rotated constraints (connection through collocation)

  subroutine constr_node_rotated_vec ( mesh, problem, constr, elem, &
    node, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: i

!   connection through collocation

    elemmat  = 0
    elemvec  = 0
    do i = 1, size(elemmat,1)
      elemmat(i,i) = 1
    end do
    elemmat2 = -Rmat

  end subroutine constr_node_rotated_vec


! Element for the rotated constraints (connection through collocation)

  subroutine constr_node_rotated_vec2 ( mesh, problem, constr, elem, &
    node, matrix, vector, first, last, coefficients, oldvectors, elemmat, &
    elemmat2, elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: i

!   connection through collocation

    elemmat  = 0
    elemvec  = 0
    do i = 1, size(elemmat,1)
      elemmat(i,i) = 1
    end do
    elemmat2 = -Rmat(1:2,:)
    elemmat2(1:2,3) = 0

  end subroutine constr_node_rotated_vec2

end program taylor_couette1
