! This program does not solve any problem but it
! computes persistence-of-straining parameter
! for a given velocity field on a cubic domain.

program persist1

  use tfem_m
  use stokes_elements_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    nx = 10,            & ! number of elements in x-direction
    ny = 10,             & ! number of elements in y-direction
    nz = 10,            & ! number of elements in z-direction
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3             ! 3x3x3 integration of hexahedra

  real(dp), parameter :: &
    lx = 1._dp,          & ! size in x-direction
    ly = 1._dp,          & ! size in y-direction
    lz = 1._dp             ! size in z-direction

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysvector_t), target :: sol
  type(vector_t) :: fclass
  type(vector_t), target :: strain
  type(meshgen_options_t) :: meshgen_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  call set_mesh_options ( meshgen_options, nx=nx, ny=ny, nz=nz, lx=lx, ly=ly, &
    lz=lz, elshape=14, regionshape=4 )

  call hexahedron ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=1 )

  input_probdef%vec_elementdof(1)%a(:,1) = 3  ! velocity
  input_probdef%vec_elementdof(1)%a(:,2) = 1  ! scalar
  input_probdef%vec_elementdof(1)%a(:,3) = 6  ! strain-rate tensor

  input_probdef%physq = [1]

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution)

  call create_sysvector ( problem, sol )

! fill solution vector

  call fill_sysvector ( mesh, problem, sol, &
    node1=1, node2=mesh%nnodes, physq=1, vfunc=vfunc, vfuncnr=2 )

! create the structure oldvectors

  call create_oldvectors ( oldvectors, nsysvec=1, nvec=3 )

! post-processing

  call create_vector ( problem, fclass, vec=2 )
  call create_vector ( problem, strain, vec=3 )

  oldvectors%s(1)%p => sol
  oldvectors%v(3)%p => strain

  call derive_vector ( mesh, problem, strain, elemsub=stokes_D_tensor, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, fclass, elemsub=stokes_P_scalar, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a vtk file for further post-processing

  call write_vector_vtk ( mesh, problem, filename='cavity.vtk', &
    dataname='velocity_vector', sysvector=sol )

  call write_tensor_vtk ( mesh, problem, filename='cavity.vtk', &
    dataname='D_tensor', vector=strain, append=.true. )

  call write_scalar_vtk ( mesh, problem, filename='cavity.vtk', &
    dataname='R', vector=fclass, append=.true. )


! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol )
  call delete ( strain, fclass )
  call delete ( coefficients )
  call delete ( oldvectors )

  contains

!   velocity functions for shear and shearfree flows
!   from
!   Dynamics of Polymeric Liquid, Volume 1: Fluid Mechanics
!   2nd Edition, p. 100-101 (1987)

    function vfunc ( n, nr, x )

      integer, intent(in) :: n, nr
      real(dp), intent(in), dimension(:) :: x
      real(dp), dimension(n) :: vfunc

      real(dp) :: rate, etype

!     deformation rate

      rate = 1._dp

!     extensional type ( 0 <= etype <= 1 )

      etype = 0._dp

      select case(nr)

        case(1)

!         shear flow

          vfunc(1) = rate * x(2)
          vfunc(2) = 0
          vfunc(3) = 0

        case(2)

!         shearfree flow

          vfunc(1) = -0.5_dp * rate * ( 1._dp + etype ) * x(1)
          vfunc(2) = -0.5_dp * rate * ( 1._dp - etype ) * x(2)
          vfunc(3) = rate * x(3)

        case default

          write(*,'(/a,i0/)') 'Error vfunc: wrong function number: ', nr
           stop

     end select

    end function vfunc

end program persist1
