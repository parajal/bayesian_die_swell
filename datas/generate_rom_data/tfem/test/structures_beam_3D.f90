
! Post processing routines for 3D beam structures

module structures_post_beam3_m

  use tfem_m
  use structures_elements_m
  use io_utils_m

  implicit none

contains

! generic routine for post-processing all data for a 3D beam structure

  subroutine postprocessing_beam3 ( mesh, problem, mshfilename, sol, load, &
    reacf, coefficients, mcoefficients )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the gmesh postprocessing output filename for writing the data to
    character (len=*), intent(in) :: mshfilename

!   sysvector containing the solution
    type(sysvector_t), target, intent(in) :: sol

!   sysvectors containing the load and the reaction forces
    type(sysvector_t), intent(in), optional :: load, reacf

!   if coefficients is present, it is passed through to the element subroutine
    type(coefficients_t), intent(inout), optional :: coefficients

!   array of coefficients (length=number of element groups)
!   If present only one of them is passed through to the element subroutine
!   based on the element group number
    type(coefficients_t), dimension(:), intent(inout), optional :: mcoefficients


!   write solution to ASCII file

 !   call printtofile ( mesh, problem, filename='sol.out', sysvector=sol )

!   direct write of sysvectors load and reacf

    if ( present(load) ) call write_sysvector_beam3 ( mesh, problem, &
      mshfilename, dataname='load', sysvector=load )

    if ( present(reacf) ) call write_sysvector_beam3 ( mesh, problem, &
      mshfilename, dataname='reaction', sysvector=reacf, append=present(load) )

    if ( maxval(mesh%elnumnod) > 2 ) then

!     additional nodes present: write displacements via deriv

      call derive_write_vector_beam3 ( mesh, problem, mshfilename, &
        dataname='displacement', sysvector=sol, &
        coefficients=coefficients, mcoefficients=mcoefficients, &
        append=present(load).or.present(reacf) )

    else

!     only end nodes present: directly write from sysvector

      print *, 'maxval abs(sol)', maxval(abs(sol%u))
 !     call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
 !       dataname='displacement', sysvector=sol, degfd=[1,2,3], &
 !       append=present(load).or.present(reacf) )

    end if

!   write of distributed load via deriv

    call derive_write_vector_beam3 ( mesh, problem, mshfilename, &
      dataname='distributed force', sysvector=sol, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      append=.true. )

!   write of beam scalar via deriv

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      dataname='transverse displacement v', sysvector=sol, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='slope v''', sysvector=sol, append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='second derivative v''''', sysvector=sol, append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='bending moment Mz', sysvector=sol, append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='shear force V', sysvector=sol, append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='axial displacement', sysvector=sol, append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='axial strain', sysvector=sol, append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='axial force', sysvector=sol, append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='axial (torsional) rotation theta', sysvector=sol, &
      append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='torsional strain dtheta/ds', sysvector=sol, append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='torsional moment Mx', sysvector=sol, append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      dataname='transverse displacement w', sysvector=sol, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='slope w''', sysvector=sol, append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='second derivative w''''', sysvector=sol, append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='bending moment My', sysvector=sol, append=.true. )

    call derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='shear force W', sysvector=sol, append=.true. )

  end subroutine postprocessing_beam3


! direct write of system vectors (load and readcf) for a 3D beam structure

  subroutine write_sysvector_beam3 ( mesh, problem, mshfilename, dataname, &
    sysvector, append )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the gmesh postprocessing output filename for writing the data to
    character (len=*), intent(in) :: mshfilename

!   the data to be written
    character (len=*), intent(in) :: dataname

!   sysvector to be written
    type(sysvector_t), intent(in) :: sysvector

!   if .true.: append data to existing file in order to avoid writing
!   coordinate and topology data more than ones.
!   default=.false.
    logical, intent(in), optional :: append

!    character (len=10) :: filename
    character (len=20) :: dataname1, dataname2
    type(vector_t) :: allnodes

!   write

    select case ( dataname )
    case('load')
      dataname1 = 'external forces'
      dataname2 = 'external moments'
    case('reaction')
      dataname1 = 'reaction forces'
      dataname2 = 'reaction moments'
      case default
        write(*,'(/a/)') 'Error: wrong function number: '
        stop
    end select

    if ( maxval(mesh%elnumnod) > 2 ) then

!     additional nodes present: fill only end nodes with data

      call transfer_data_endnodes

      print *, dataname1, dataname2, &
                  'maxval abs(allnodes)', maxval(abs(allnodes%u))
!!     forces
!!      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
!        dataname=dataname1, vector=allnodes, degfd=[1,2,3], &
!        append=append )
!!     moments
!      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
!        dataname=dataname2, vector=allnodes, degfd=[4,5,6], &
!        append=.true. )

      call delete ( allnodes )

    else

!     only end nodes present: directly write from sysvector

      print *, dataname1, dataname2, &
               'maxval abs(sysvector)', maxval(abs(sysvector%u))

!!     forces
!      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
!        dataname=dataname1, sysvector=sysvector, degfd=[1,2,3], &
!        append=append )
!!     moments
!      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
!        dataname=dataname2, sysvector=sysvector, degfd=[4,5,6], &
!        append=.true. )

    end if

!   print

!    select case ( dataname )
!    case('load')
!      filename = 'load.out'
!    case('reaction')
!      filename = 'reacf.out'
!    end select

!    call printtofile ( mesh, problem, filename=filename, sysvector=sysvector )

  contains

!   transfer data values of end nodes only and set internal nodes to zero
!   requires vec=3 and vec=4 to be defined.

    subroutine transfer_data_endnodes

      type(vector_t) :: endnodes

 !    create vector (vec=3) that contains all nodes for plotting
      call create ( problem, allnodes, vec=3 )

 !    create vector (vec=4) that contains only end nodes, e.g. [6,0,6]
      call create ( problem, endnodes, vec=4 )

      call transfer_data ( mesh, problem, sysvector1=sysvector, &
        vector2=endnodes, degfd1=[1,2,3,4,5,6], degfd2=[1,2,3,4,5,6] )

      allnodes%u = 0 ! set all values to zero, including internal nodes

      call transfer_data ( mesh, problem, vector1=endnodes, &
        vector2=allnodes, degfd1=[1,2,3,4,5,6], degfd2=[1,2,3,4,5,6] )

      call delete ( endnodes )

    end subroutine transfer_data_endnodes

  end subroutine write_sysvector_beam3


! post-processing scalar data using derive for a 3D beam structure

  subroutine derive_write_scalar_beam3 ( mesh, problem, mshfilename, &
    dataname, sysvector, coefficients, mcoefficients, append )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the gmesh postprocessing output filename for writing the data to
    character (len=*), intent(in) :: mshfilename

!   the data to be written
    character (len=*), intent(in) :: dataname

!   sysvector containing the solution
    type(sysvector_t), target, intent(in) :: sysvector

!   if coefficients is present, it is passed through to the element subroutine
    type(coefficients_t), intent(inout), optional :: coefficients

!   array of coefficients (length=number of element groups)
!   If present only one of them is passed through to the element subroutine
!   based on the element group number
    type(coefficients_t), dimension(:), intent(inout), optional :: mcoefficients

!   if .true.: append data to existing file in order to avoid writing
!   coordinate and topology data more than ones.
!   default=.false.
    logical, intent(in), optional :: append

!    character (len=12) :: filename
    integer :: i
    type(vector_t) :: dvec
    type(oldvectors_t) :: oldvectors

    call create ( oldvectors, nsysvec=1 )
    oldvectors%s(1)%p => sysvector

    call create ( problem, dvec, vec=1, elementwise=.true. )

    if ( present(mcoefficients) ) then

      do i = 1, size(mcoefficients)

        select case ( dataname )
        case('transverse displacement v')
          mcoefficients(i)%i(1) = 1 ! transverse displacement v
        case('slope v''')
          mcoefficients(i)%i(1) = 2 ! slope v'
        case('second derivative v''''')
          mcoefficients(i)%i(1) = 3 ! second derivative v''
        case('bending moment Mz')
          mcoefficients(i)%i(1) = 4 ! bending moment Mz
        case('shear force V')
          mcoefficients(i)%i(1) = 5 ! shear force V
        case('axial displacement')
          mcoefficients(i)%i(1) = 6 ! axial displacement
        case('axial strain')
          mcoefficients(i)%i(1) = 7 ! axial strain
        case('axial force')
          mcoefficients(i)%i(1) = 8 ! axial force
        case('axial (torsional) rotation theta')
          mcoefficients(i)%i(1) = 9 ! axial (torsional) rotation theta
        case('torsional strain dtheta/ds')
          mcoefficients(i)%i(1) = 10 ! torsional strain dtheta/ds
        case('torsional moment Mx')
          mcoefficients(i)%i(1) = 11 ! torsional moment Mx
        case('transverse displacement w')
          mcoefficients(i)%i(1) = 12 ! transverse displacement w
        case('slope w''')
          mcoefficients(i)%i(1) = 13 ! slope w'
        case('second derivative w''''')
          mcoefficients(i)%i(1) = 14 ! second derivative w''
        case('bending moment My')
          mcoefficients(i)%i(1) = 15 ! bending moment My
        case('shear force W')
          mcoefficients(i)%i(1) = 16 ! shear force W
        case default
          write(*,'(/a/)') 'Error: wrong function number: '
          stop
        end select

      end do

    else

      select case ( dataname )
      case('transverse displacement v')
        coefficients%i(1) = 1 ! transverse displacement v
      case('slope v''')
        coefficients%i(1) = 2 ! slope v'
      case('second derivative v''''')
        coefficients%i(1) = 3 ! second derivative v''
      case('bending moment Mz')
        coefficients%i(1) = 4 ! bending moment Mz
      case('shear force V')
        coefficients%i(1) = 5 ! shear force V
      case('axial displacement')
        coefficients%i(1) = 6 ! axial displacement
      case('axial strain')
        coefficients%i(1) = 7 ! axial strain
      case('axial force')
        coefficients%i(1) = 8 ! axial force
      case('axial (torsional) rotation theta')
        coefficients%i(1) = 9 ! axial (torsional) rotation theta
      case('torsional strain dtheta/ds')
        coefficients%i(1) = 10 ! torsional strain dtheta/ds
      case('torsional moment Mx')
        coefficients%i(1) = 11 ! torsional moment Mx
      case('transverse displacement w')
        coefficients%i(1) = 12 ! transverse displacement w
      case('slope w''')
        coefficients%i(1) = 13 ! slope w'
      case('second derivative w''''')
        coefficients%i(1) = 14 ! second derivative w''
      case('bending moment My')
        coefficients%i(1) = 15 ! bending moment My
      case('shear force W')
        coefficients%i(1) = 16 ! shear force W
      case default
        write(*,'(/a/)') 'Error: wrong function number: '
        stop
      end select

    end if

!   derive

    call derive_vector ( mesh, problem, dvec, elemsub=beam_deriv3, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      oldvectors=oldvectors )

!   write

      print *, dataname, 'maxval abs(dvec)', maxval(abs(dvec%u))
!    call write_scalar_gmsh ( mesh, problem, filename=mshfilename, &
!      dataname=dataname, vector=dvec, append=append )

!    select case ( dataname )
!    case('transverse displacement v')
!      filename = 'vy.out'
!    case('slope v''')
!      filename = 'dvdx.out'
!    case('second derivative v''''')
!      filename = 'd2vdx2.out'
!    case('bending moment Mz')
!      filename = 'Mz.out'
!    case('shear force V')
!      filename = 'V.out'
!    case('axial displacement')
!      filename = 'u.out'
!    case('axial strain')
!      filename = 'eps.out'
!    case('axial force')
!      filename = 'N.out'
!    case('axial (torsional) rotation theta')
!      filename = 'theta.out'
!    case('torsional strain dtheta/ds')
!      filename = 'dthetads.out'
!    case('torsional moment Mx')
!      filename = 'Mx.out'
!    case('transverse displacement w')
!      filename = 'wz.out'
!    case('slope w''')
!      filename = 'dwdx.out'
!    case('second derivative w''''')
!      filename = 'd2wdx2.out'
!    case('bending moment My')
!      filename = 'My.out'
!    case('shear force W')
!      filename = 'W.out'
!    end select
!
!!   print
!
!    call printtofile ( mesh, problem, filename=filename, vector=dvec )

    call delete ( dvec )
    call delete ( oldvectors )

  end subroutine derive_write_scalar_beam3


! post-processing vector data using derive for a 3D beam structure

  subroutine derive_write_vector_beam3 ( mesh, problem, mshfilename, &
    dataname, sysvector, coefficients, mcoefficients, append )

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the gmesh postprocessing output filename for writing the data to
    character (len=*), intent(in) :: mshfilename

!   the data to be written
    character (len=*), intent(in) :: dataname

!   sysvector containing the solution
    type(sysvector_t), target, intent(in) :: sysvector

!   if coefficients is present, it is passed through to the element subroutine
    type(coefficients_t), intent(inout), optional :: coefficients

!   array of coefficients (length=number of element groups)
!   If present only one of them is passed through to the element subroutine
!   based on the element group number
    type(coefficients_t), dimension(:), intent(inout), optional :: mcoefficients

!   if .true.: append data to existing file in order to avoid writing
!   coordinate and topology data more than ones.
!   default=.false.
    logical, intent(in), optional :: append

!    character (len=10) :: filename
    type(vector_t) :: dvec
    type(oldvectors_t) :: oldvectors

    call create ( oldvectors, nsysvec=1 )
    oldvectors%s(1)%p => sysvector

    select case ( dataname )
    case('displacement')

      call create ( problem, dvec, vec=2 )

!     derive
      call derive_vector ( mesh, problem, dvec, elemsub=beam_displacement3, &
        coefficients=coefficients, mcoefficients=mcoefficients, &
        oldvectors=oldvectors, order='ND' )

    case('distributed force')

      call create ( problem, dvec, vec=2, elementwise=.true. )

!     derive
      call derive_vector ( mesh, problem, dvec, elemsub=beam_q3, &
        coefficients=coefficients, mcoefficients=mcoefficients, &
        oldvectors=oldvectors, order='ND' )

      case default
        write(*,'(/a/)') 'Error: wrong function number: '
        stop
    end select

!   write

      print *, dataname, 'maxval abs(dvec)', maxval(abs(dvec%u))
!    call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
!      dataname=dataname, vector=dvec, append=append )

!    select case ( dataname )
!!    case('displacement')
!      filename = 'disp.out'
!    case('distributed force')
!      filename = 'q.out'
!    end select

!   print

!    call printtofile ( mesh, problem, filename=filename, vector=dvec )

    call delete ( dvec )
    call delete ( oldvectors )

  end subroutine derive_write_vector_beam3

end module structures_post_beam3_m

! A two-beam structure problem. Gravity load. 3D beam elements.
! Forces and moments in point 2 (tip of triangle).
! Gmsh mesh.
! P1/P3/P1 interpolation for axial displacement, transverse displacement
! and axial rotation, respectively.
! Extension of structures11 to a 3D beam structure.

program structures16

  use tfem_m
  use hsl_ma57_m
  use structures_elements_m
  use structures_post_beam3_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    ninti = 3 ! number of Gauss points

  real(dp), parameter :: &
    rho = 0._dp,  & ! density
    g = 1._dp,    & ! gravitational acceleration
    Emod = 1._dp, &  ! E modulus
    Gmod = 1._dp, &  ! G modulus
    Ac = 1._dp,   &  ! cross-sectional area
    Izz = 1._dp,  &  ! second moment of area Izz
    Iyy = 2._dp,  &  ! second moment of area Iyy
    Jrr = 1._dp,  &  ! polar moment of area J
!   external forces/moments in point 2
    F(3) = [ 1._dp, 2._dp, 10._dp ], &  ! external force vector
    M(3) = [ 10._dp, 0._dp, 3._dp ]    ! external moment

! definitions

  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: load, rhsd, reacf
  type(coefficients_t) :: coefficients


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i = 0
  coefficients%i(2) = ninti
  coefficients%i(3) = 0 ! P3 Hermite interpolation for v
  coefficients%i(4) = 0 ! constant distributed force
  coefficients%i(12) = 2 ! P1 interpolation for u
  coefficients%i(13) = 2 ! P1 interpolation for theta

  coefficients%r = 0
  coefficients%r(1) = Ac
  coefficients%r(2) = Emod
  coefficients%r(3) = Izz
  coefficients%r(5:7) = [ 0._dp, -rho*g*Ac, 0._dp ]
  coefficients%r(8:10) = [ 0._dp, 0._dp, 1._dp ] ! local z-axis direction
  coefficients%r(11) = Iyy
  coefficients%r(12) = 0   ! Izy = 0
  coefficients%r(13) = Jrr
  coefficients%r(14) = Gmod


! read mesh

  call read_mesh_gmsh ( mesh, filename='triangle.msh', ndim=3 )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  input_probdef%elementdof(1)%a = 6  ! (u,v,w,phix,phiy,phiz) in two end nodes
  input_probdef%vec_elementdof(1)%a = 1      ! scalar quantity
  input_probdef%vec_elementdof(1)%a(:,2) = 3 ! vector quantity

! builtin at points 1 and 3

  call define_essential ( mesh, input_probdef, point=1 )
  call define_essential ( mesh, input_probdef, point=3 )

  call problem_definition ( input_probdef, mesh, problem )

  call printinfo ( mesh, printlevel=4 )
  call printinfo ( problem, printlevel=4 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, load, rhsd, reacf )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, point=3, value=0._dp )

! fill rhs vector with discrete forces/monents at point 2

  load%u = 0

  call fill_sysvector ( mesh, problem, load, point=2, degfd=1, value=F(1) )
  call fill_sysvector ( mesh, problem, load, point=2, degfd=2, value=F(2) )
  call fill_sysvector ( mesh, problem, load, point=2, degfd=3, value=F(3) )
  call fill_sysvector ( mesh, problem, load, point=2, degfd=4, value=M(1) )
  call fill_sysvector ( mesh, problem, load, point=2, degfd=5, value=M(2) )
  call fill_sysvector ( mesh, problem, load, point=2, degfd=6, value=M(3) )

  call copy ( load, rhsd )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=beam_elem3, &
    order='ND', coefficients=coefficients, addvec=.true. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

! post processing

  call postprocessing_beam3 ( mesh, problem, coefficients=coefficients, &
    mshfilename='structures16.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures16

