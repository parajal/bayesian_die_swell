! Stokes problem on a square domain with a disk at the center.
! Shear flow boundary conditions.
! Freely rotating conditions on the disk boundaries.
! Gmsh mesh.
! L2-projection of velocity solution on a second mesh

module subs_m

use tfem_elem_m

  implicit none

contains

  subroutine elementc ( mesh, problem, constr, elem, node, &
    matrix, vector, first, last, coefficients, oldvectors, elemmat, elemmat2, &
    elemmatadd, elemvec, elemvecadd )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: constr, elem, node
    logical, intent(in) :: matrix, vector, first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:) :: elemmat, elemmat2, elemmatadd
    real(dp), intent(out), dimension(:) :: elemvec, elemvecadd

    integer :: ncurveconstr
    real(dp) :: xr(1,2)

    ncurveconstr = problem%constraints(constr)%geometry1

    xr(1,:) = mesh%coor(mesh%curves(ncurveconstr)%nodes(node),:)

    if ( vector ) then

      elemvec = 0
      elemvecadd = 0

    end if

!   connection through collocation

    elemmat(1,1) = 1
    elemmat(1,2) = 0

    elemmat(2,1) = 0
    elemmat(2,2) = 1

    elemmatadd(1,:) = [ -xr(1,2) ]
    elemmatadd(2,:) = [ xr(1,1) ]

  end subroutine elementc

end module subs_m

program projection2

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use projection_elements_m
  use figplot_m
  use io_utils_m
  use subs_m
  use functions_m

  implicit none

  integer, parameter :: &
    uintpl = 6,         & ! P2 velocities
    pintpl = 2,         & ! P1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 6,          & ! 6 point integration of triangles
    gaussb = 3,         & ! 3 point integration of boundary elements
    gaussp = 6,         & ! 6 point integration of triangles for projection
    ncompv = 2            ! number of components input vector for projection

! definitions

  type(mesh_t), target :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd
  type(vector_t), target :: vec_p
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients


  type(mesh_t) :: mesh2
  type(input_probdef_t) :: input_probdef2
  type(problem_t) :: problem2
  type(sysmatrix_t) :: sysmatrix2
  type(sysvector_t) :: sol2(ncompv)
  type(sysvector_t) :: rhsd2(ncompv)
  type(oldvectors_t) :: oldvectors2
  type(coefficients_t) :: coefficients2
  type(lu_ma57_t) :: lu2

! variables

  integer :: ip, i, nb
  real(dp), parameter :: eta = 1._dp     ! viscosity

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gaussb ]
  coefficients%i(12:) = 0
  coefficients%r(1) = eta
  coefficients%r(2:) = 0

! read mesh from gmsh output file

  call read_mesh_gmsh ( mesh, filename='disk1.msh', ndim=2, physgeom=.true. )

  nb = nint( real(mesh%nelem) ** 0.25 ); print *, 'nb=', nb

  call add_to_mesh ( mesh, blocks=[nb,nb] )

  call fill_mesh_parts ( mesh )

  call printinfo ( mesh, printlevel=4 )

! plot mesh

  plot_options%fontsize = 6
  call plot_points_curves ( plot_options, mesh, 'curves.fig' )

  plot_options%fontsize = 10

  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=4, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,    &  ! pressure
                 1,1,1,1,1,1,    &  ! scalar
                 3,3,3,3,3,3 ], &  ! tensor
                 [6,4] )

  input_probdef%physq = [1,2]

! define essential boundaries

! cell boundaries
  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )

! define constraints on the disk
  call define_constraint ( mesh, input_probdef, curve1=5, physq=1, &
    discretization='collocation', nodedof=2, naddunknowns=1 )

! pressure level
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd )

! fill boundary conditions

! cell boundaries
  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, physq=1, degfd=1, func=func, funcnr=1 )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, curve2=4, physq=1, degfd=2, func=func, funcnr=2 )

! pressure level
  call fill_sysvector ( mesh, problem, sol, point=1, physq=2, value=0._dp )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1, nvec=2 )

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      elemsub=elementc, coefficients=coefficients, addmatvec=.true. )

  call check_filled_sysmatrix ( sysmatrix )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call delete ( sysmatrix )

  ip = problem%constraints(1)%addnumdegfd(1)
  print *,'Rotation rate of disk: ', sol%u( problem%degfdperm(ip+1,2) )

  oldvectors%s(1)%p => sol

! write velocity
  call write_scalar_vtk ( mesh, problem, filename='velocityxorig.vtk', &
    dataname='velocity', sysvector=sol, physq=physqvel, degfd=1 )
  call write_scalar_vtk ( mesh, problem, filename='velocityyorig.vtk', &
    dataname='velocity', sysvector=sol, physq=physqvel, degfd=2 )

! create vector and copy data for input of projection problem

  call create ( problem, vec_p, vec=1 )

  call transfer_data ( mesh, problem, sysvector1=sol, vector2=vec_p, &
    physq1=[1] )


!****************************************************************
! PROJECTION PROBLEM
!****************************************************************

! read mesh from gmsh output file

  call read_mesh_gmsh ( mesh2, filename='disk2.msh', ndim=2, physgeom=.true. )

  call fill_mesh_parts ( mesh2 )

  call printinfo ( mesh2, printlevel=4 )

! problem definition for projection

  call create_input_probdef ( mesh2, input_probdef2 )

  input_probdef2%elementdof(1)%a = 1

  call problem_definition ( input_probdef2, mesh2, problem2 )

! create system vectors (solution and right-hand side)

  call create ( problem2, sol2 )
  call create ( problem2, rhsd2 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix2, mesh2, problem2, &
    symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix2 )

! fill coefficients

  call create_coefficients ( coefficients2, ncoefi=100, ncoefr=50 )

  coefficients2%i = 0
  coefficients2%i(1) = uintpl
  coefficients2%i(2) = ncompv
  coefficients2%i(3) = uintpl
  coefficients2%i(4) = uintpl
  coefficients2%i(10) = gaussp

  coefficients2%r = 0

! oldvectors

  call create_oldvectors ( oldvectors2, nvec=1, nprob=1, nmesh=1 )

  oldvectors2%v(1)%p => vec_p
  oldvectors2%p(1)%p => problem
  oldvectors2%m(1)%p => mesh

! build (assemble) matrix and vector from elements

  call build_system ( mesh2, problem2, sysmatrix2, msysvector=rhsd2, &
    elemsub=projection_elem, coefficients=coefficients2, &
    oldvectors=oldvectors2 )

  do i = 1, ncompv

    call add_effect_of_essential_to_rhs ( problem2, sysmatrix2, sol2(i), &
      rhsd2(i) )

    call solve_system_ma57 ( sysmatrix2, rhsd2(i), sol2(i), lu=lu2 )

  end do

  call delete ( lu2 )  ! remove LU decomposition and rebuild next time step

! write velocity
  call write_scalar_vtk ( mesh2, problem2, filename='velocityx.vtk', &
    dataname='velocity', sysvector=sol2(1) )
  call write_scalar_vtk ( mesh2, problem2, filename='velocityy.vtk', &
    dataname='velocity', sysvector=sol2(2) )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( vec_p )
  call delete ( coefficients )
  call delete ( oldvectors )

  call delete ( mesh2 )
  call delete ( problem2 )
  call delete ( input_probdef2 )
  call delete ( sol2, rhsd2 )
  call delete ( sysmatrix2 )
  call delete ( oldvectors2 )

end program projection2
