! In this program, a scalar function is used to refine a mesh. The value of
! the scalar function indicates the mesh element size at that location
program adaptive_mesh_pos

  use tfem_m
  use io_utils_m

  implicit none

  real(dp), parameter :: &
   ox = 0.0_dp,      & ! x-coordinate of origin of the mesh
   oy = 0.0_dp,      & ! y-coordinate of origin of the mesh
   lx = 1.0_dp,      & ! size in x-direction of the mesh
   ly = 1.0_dp,      & ! size in y-direction of the mesh
   dx = 0.1_dp         ! elementsize of the initial mesh

  type(mesh_t), target :: mesh,  mesh_refined
  type(input_probdef_t) :: input_probdef
  type(problem_t), target :: problem
  type(vector_t), target :: element_size

  logical :: bg_mesh

! generate the initial mesh on which the scalar function will be defined
  bg_mesh=.false.
  call generate_read_mesh ( mesh, ox, oy, lx, ly, dx )

  call fill_mesh_parts ( mesh )

  call write_mesh_vtk ( mesh, 'mesh_init.vtk' )

! problem definition
  call create_input_probdef ( mesh, input_probdef, nvec=1 )

  input_probdef%elementdof(1)%a = 1
  input_probdef%vec_elementdof(1)%a = reshape ( [1,1,1,1,1,1], [6,1] )

  call problem_definition ( input_probdef, mesh, problem )

! create a scalar vector for the element size
  call create_vector (problem, element_size, vec=1 )

! determine the target element size
  element_size%u = (mesh%coor(:,1) + 1.0_dp)/40.0_dp

  call write_scalar_vtk ( mesh, problem, vector=element_size, &
    dataname='element_size', filename='element_size.vtk' )

! write the function values on the initial mesh to a .pos file
  call write_scalar_gmsh_parsed ( mesh, problem, vector=element_size, &
    dataname='element_size', filename='bg.pos' )

! generate the mesh using the background mesh for refinement
  bg_mesh = .true.
  call generate_read_mesh ( mesh_refined, ox, oy, lx, ly, dx )
  call fill_mesh_parts ( mesh_refined )

! write the refined mesh
  call write_mesh_vtk ( mesh_refined, filename='mesh_refined.vtk' )

! delete definitions
  call delete ( mesh, mesh_refined )
  call delete ( input_probdef )
  call delete ( problem )

contains

! generate and read mesh

  subroutine generate_read_mesh ( mesh, ox, oy, lx, ly, dx )

    type(mesh_t), intent(inout) :: mesh
    real(dp), intent(in) :: ox, oy, lx, ly, dx

    open ( unit=25, file='mesh.geo' )
    write ( 25, '(1X,A,F18.14,A)' ) 'ox = ', ox, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'oy = ', oy, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'lx = ', lx, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'ly = ', ly, ';'
    write ( 25, '(1X,A,F18.14,A)' ) 'dx = ', dx, ';'
    write ( 25, '(/1x,a)' ) 'Include "mesh_quad.igo";'
    close ( 25 )

    if ( bg_mesh ) then
      call execute_command_line ( 'gmsh -bgm bg.pos -2 -order 2 -o mesh.msh &
        & mesh.geo > outputmesh.out' )
    else
      call execute_command_line ( 'gmsh -2 -order 2 -o mesh.msh &
        & mesh.geo > outputmesh.out' )
    end if

    call read_mesh_gmsh ( mesh, filename='mesh.msh', ndim=2 )

  end subroutine generate_read_mesh

end program adaptive_mesh_pos
