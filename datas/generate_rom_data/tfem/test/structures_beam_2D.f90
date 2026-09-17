
! Post processing routines for 2D beam structures

module structures_post_beam2_m

  use tfem_m
  use structures_elements_m
  use io_utils_m

  implicit none

contains

! generic routine for post-processing all data for a 2D beam structure

  subroutine postprocessing_beam2 ( mesh, problem, mshfilename, sol, load, &
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

!    call printtofile ( mesh, problem, filename='sol.out', sysvector=sol )

!   direct write of sysvectors load and reacf

    if ( present(load) ) call write_sysvector_beam2 ( mesh, problem, &
      mshfilename, dataname='load', sysvector=load )

    if ( present(reacf) ) call write_sysvector_beam2 ( mesh, problem, &
      mshfilename, dataname='reaction', sysvector=reacf, append=present(load) )

    if ( maxval(mesh%elnumnod) > 2 ) then

!     additional nodes present: write displacements via deriv

      call derive_write_vector_beam2 ( mesh, problem, mshfilename, &
        dataname='displacement', sysvector=sol, &
        coefficients=coefficients, mcoefficients=mcoefficients, &
        append=present(load).or.present(reacf) )

    else

!     only end nodes present: directly write from sysvector

      print *, 'maxval abs(sol)', maxval(abs(sol%u))
!      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
!        dataname='displacement', sysvector=sol, degfd=[1,2,0], &
!        append=present(load).or.present(reacf) )

    end if

!   write of distributed load via deriv

    call derive_write_vector_beam2 ( mesh, problem, mshfilename, &
      dataname='distributed force', sysvector=sol, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      append=.true. )

!   write of beam scalar via deriv

    call derive_write_scalar_beam2 ( mesh, problem, mshfilename, &
      dataname='transverse displacement', sysvector=sol, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      append=.true. )

    call derive_write_scalar_beam2 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='slope', sysvector=sol, append=.true. )

    call derive_write_scalar_beam2 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='second derivative', sysvector=sol, append=.true. )

    call derive_write_scalar_beam2 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='bending moment', sysvector=sol, append=.true. )

    call derive_write_scalar_beam2 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='shear force', sysvector=sol, append=.true. )

    call derive_write_scalar_beam2 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='axial displacement', sysvector=sol, append=.true. )

    call derive_write_scalar_beam2 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='axial strain', sysvector=sol, append=.true. )

    call derive_write_scalar_beam2 ( mesh, problem, mshfilename, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      dataname='axial force', sysvector=sol, append=.true. )

  end subroutine postprocessing_beam2


! direct write of system vectors (load and readcf) for a 2D beam structure

  subroutine write_sysvector_beam2 ( mesh, problem, mshfilename, dataname, &
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
!      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
!        dataname=dataname1, vector=allnodes, degfd=[1,2,0], &
!        append=append )
!!     moments
!      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
!        dataname=dataname2, vector=allnodes, degfd=[0,0,3], &
!        append=.true. )

      call delete ( allnodes )

    else

!     only end nodes present: directly write from sysvector

      print *, dataname1, dataname2, &
               'maxval abs(sysvector)', maxval(abs(sysvector%u))
!     forces
!      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
!        dataname=dataname1, sysvector=sysvector, degfd=[1,2,0], &
!        append=append )
!     moments
!      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
!        dataname=dataname2, sysvector=sysvector, degfd=[0,0,3], &
!        append=.true. )

    end if

!   print

!    select case ( dataname )
!    case('load')
!      filename = 'load.out'
!    case('reaction')
!      filename = 'reacf.out'
!    end select
!
!    call printtofile ( mesh, problem, filename=filename, sysvector=sysvector )

  contains

!   transfer data values of end nodes only and set internal nodes to zero
!   requires vec=3 and vec=4 to be defined.

    subroutine transfer_data_endnodes

      type(vector_t) :: endnodes

 !    create vector (vec=3) that contains all nodes for plotting
      call create ( problem, allnodes, vec=3 )

 !    create vector (vec=4) that contains only end nodes, e.g. [3,0,3]
      call create ( problem, endnodes, vec=4 )

      call transfer_data ( mesh, problem, sysvector1=sysvector, &
        vector2=endnodes, degfd1=[1,2,3], degfd2=[1,2,3] )

      allnodes%u = 0 ! set all values to zero, including internal nodes

      call transfer_data ( mesh, problem, vector1=endnodes, &
        vector2=allnodes, degfd1=[1,2,3], degfd2=[1,2,3] )

      call delete ( endnodes )

    end subroutine transfer_data_endnodes

  end subroutine write_sysvector_beam2


! post-processing scalar data using derive for a 2D beam structure

  subroutine derive_write_scalar_beam2 ( mesh, problem, mshfilename, &
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
    integer :: i
    type(vector_t) :: dvec
    type(oldvectors_t) :: oldvectors

    call create ( oldvectors, nsysvec=1 )
    oldvectors%s(1)%p => sysvector

    call create ( problem, dvec, vec=1, elementwise=.true. )

    if ( present(mcoefficients) ) then

      do i = 1, size(mcoefficients)

        select case ( dataname )
        case('transverse displacement')
          mcoefficients(i)%i(1) = 1 ! transverse displacement
        case('slope')
          mcoefficients(i)%i(1) = 2 ! slope
        case('second derivative')
          mcoefficients(i)%i(1) = 3 ! second derivative
        case('bending moment')
          mcoefficients(i)%i(1) = 4 ! bending moment
        case('shear force')
          mcoefficients(i)%i(1) = 5 ! shear force
        case('axial displacement')
          mcoefficients(i)%i(1) = 6 ! axial displacement
        case('axial strain')
          mcoefficients(i)%i(1) = 7 ! axial strain
        case('axial force')
          mcoefficients(i)%i(1) = 8 ! axial force
        case default
          write(*,'(/a/)') 'Error: wrong function number: '
          stop
        end select

      end do

    else

      select case ( dataname )
      case('transverse displacement')
        coefficients%i(1) = 1 ! transverse displacement
      case('slope')
        coefficients%i(1) = 2 ! slope
      case('second derivative')
        coefficients%i(1) = 3 ! second derivative
      case('bending moment')
        coefficients%i(1) = 4 ! bending moment
      case('shear force')
        coefficients%i(1) = 5 ! shear force
      case('axial displacement')
        coefficients%i(1) = 6 ! axial displacement
      case('axial strain')
        coefficients%i(1) = 7 ! axial strain
      case('axial force')
        coefficients%i(1) = 8 ! axial force
      case default
        write(*,'(/a/)') 'Error: wrong function number: '
        stop
      end select

    end if

!   derive

    call derive_vector ( mesh, problem, dvec, elemsub=beam_deriv2, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      oldvectors=oldvectors )

!   write

      print *, dataname, 'maxval abs(dvec)', maxval(abs(dvec%u))
!    call write_scalar_gmsh ( mesh, problem, filename=mshfilename, &
!      dataname=dataname, vector=dvec, append=append )

!    select case ( dataname )
!    case('transverse displacement')
!      filename = 'vy.out'
!    case('slope')
!      filename = 'dvdx.out'
!    case('second derivative')
!      filename = 'd2vdx2.out'
!    case('bending moment')
!      filename = 'M.out'
!    case('shear force')
!      filename = 'V.out'
!    case('axial displacement')
!      filename = 'u.out'
!    case('axial strain')
!      filename = 'eps.out'
!    case('axial force')
!      filename = 'N.out'
!    end select

!   print

!    call printtofile ( mesh, problem, filename=filename, vector=dvec )

    call delete ( dvec )
    call delete ( oldvectors )

  end subroutine derive_write_scalar_beam2


! post-processing vector data using derive for a 2D beam structure

  subroutine derive_write_vector_beam2 ( mesh, problem, mshfilename, &
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
      call derive_vector ( mesh, problem, dvec, elemsub=beam_displacement2, &
        coefficients=coefficients, mcoefficients=mcoefficients, &
        oldvectors=oldvectors, order='ND' )

    case('distributed force')

      call create ( problem, dvec, vec=2, elementwise=.true. )

!     derive
      call derive_vector ( mesh, problem, dvec, elemsub=beam_q2, &
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
!    case('displacement')
!      filename = 'disp.out'
!    case('distributed force')
!      filename = 'q.out'
!    end select

!   print

!    call printtofile ( mesh, problem, filename=filename, vector=dvec )

    call delete ( dvec )
    call delete ( oldvectors )

  end subroutine derive_write_vector_beam2

end module structures_post_beam2_m

! A two-beam structure problem. Gravity load.
! Forces and moment in point 2 (tip of triangle).
! P1/P3 interpolation. 2D elements.
! Gmsh mesh.

program structures11

  use tfem_m
  use hsl_ma57_m
  use structures_elements_m
  use structures_post_beam2_m
  use io_utils_m

  implicit none


! constants

  integer, parameter :: &
    ninti = 3 ! number of Gauss points

  real(dp), parameter :: &
    rho = 1._dp,  & ! density
    g = 1._dp,    & ! gravitational acceleration
    Emod = 1._dp, &  ! E modulus
    Ac = 1._dp,   &  ! cross-sectional area
    Iz = 1._dp,   &  ! second moment of area
!   external forces/moments in point 2
    F(2) = [ 1._dp,  2._dp ], &  ! external force vector
    M = 3._dp     ! external moment

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

  coefficients%r = 0
  coefficients%r(1) = Ac
  coefficients%r(2) = Emod
  coefficients%r(3) = Iz
  coefficients%r(5:6) = [ 0._dp, -rho*g*Ac ]


! read mesh

  call read_mesh_gmsh ( mesh, filename='triangle.msh', ndim=2 )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  input_probdef%elementdof(1)%a = 3          ! (u,v,phi) in two end nodes
  input_probdef%vec_elementdof(1)%a = 1      ! scalar quantity
  input_probdef%vec_elementdof(1)%a(:,2) = 2 ! vector quantity

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
  call fill_sysvector ( mesh, problem, load, point=2, degfd=3, value=M )

  call copy ( load, rhsd )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=beam_elem2, &
    order='ND', coefficients=coefficients, addvec=.true. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

! post processing

  call postprocessing_beam2 ( mesh, problem, coefficients=coefficients, &
    mshfilename='structures11.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures11

