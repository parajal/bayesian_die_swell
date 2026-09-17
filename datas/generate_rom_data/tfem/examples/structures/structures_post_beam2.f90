
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

    call printtofile ( mesh, problem, filename='sol.out', sysvector=sol )

!   direct write of sysvectors load and reacf

    if ( present(load) ) call write_sysvector_beam2 ( mesh, problem, &
      mshfilename, dataname='load', sysvector=load )

    if ( present(reacf) ) call write_sysvector_beam2 ( mesh, problem, &
      mshfilename, dataname='reaction', sysvector=reacf, append=present(load) )

    if ( maxval(mesh%elnumnod) > 2 ) then

!     additional nodes present: write displacements and rotations via deriv

      call derive_write_vector_beam2 ( mesh, problem, mshfilename, &
        dataname='displacement', sysvector=sol, &
        coefficients=coefficients, mcoefficients=mcoefficients, &
        append=present(load).or.present(reacf) )

      call derive_write_vector_beam2 ( mesh, problem, mshfilename, &
        dataname='rotation', sysvector=sol, &
        coefficients=coefficients, mcoefficients=mcoefficients, &
        append=.true. )

    else

!     only end nodes present: directly write from sysvector

      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
        dataname='displacement', sysvector=sol, degfd=[1,2,0], &
        append=present(load).or.present(reacf) )

      call write_scalar_gmsh ( mesh, problem, filename=mshfilename, &
        dataname='rotation', sysvector=sol, degfd=3, &
        append=.true. )

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

    character (len=10) :: filename
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
      write(*,'(/a,a/)') 'Error: wrong value dataname: ', dataname
      stop
    end select

    if ( maxval(mesh%elnumnod) > 2 ) then

!     additional nodes present: fill only end nodes with data

      call transfer_data_endnodes

!     forces
      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
        dataname=dataname1, vector=allnodes, degfd=[1,2,0], &
        append=append )
!     moments
      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
        dataname=dataname2, vector=allnodes, degfd=[0,0,3], &
        append=.true. )

      call delete ( allnodes )

    else

!     only end nodes present: directly write from sysvector

!     forces
      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
        dataname=dataname1, sysvector=sysvector, degfd=[1,2,0], &
        append=append )
!     moments
      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
        dataname=dataname2, sysvector=sysvector, degfd=[0,0,3], &
        append=.true. )

    end if

!   print

    select case ( dataname )
    case('load')
      filename = 'load.out'
    case('reaction')
      filename = 'reacf.out'
    case default
      write(*,'(/a,a/)') 'Error: wrong value dataname: ', dataname
      stop
    end select

    call printtofile ( mesh, problem, filename=filename, sysvector=sysvector )

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

    character (len=10) :: filename
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
          write(*,'(/a,a/)') 'Error: wrong value dataname: ', dataname
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
        write(*,'(/a,a/)') 'Error: wrong value dataname: ', dataname
        stop
      end select

    end if

!   derive

    call derive_vector ( mesh, problem, dvec, elemsub=beam_deriv2, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      oldvectors=oldvectors )

!   write

    call write_scalar_gmsh ( mesh, problem, filename=mshfilename, &
      dataname=dataname, vector=dvec, append=append )

    select case ( dataname )
    case('transverse displacement')
      filename = 'vy.out'
    case('slope')
      filename = 'dvdx.out'
    case('second derivative')
      filename = 'd2vdx2.out'
    case('bending moment')
      filename = 'M.out'
    case('shear force')
      filename = 'V.out'
    case('axial displacement')
      filename = 'u.out'
    case('axial strain')
      filename = 'eps.out'
    case('axial force')
      filename = 'N.out'
    case default
      write(*,'(/a,a/)') 'Error: wrong value dataname: ', dataname
      stop
    end select

!   print

    call printtofile ( mesh, problem, filename=filename, vector=dvec )

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

    logical :: write_scalar
    character (len=10) :: filename
    type(vector_t) :: dvec
    type(oldvectors_t) :: oldvectors

    call create ( oldvectors, nsysvec=1 )
    oldvectors%s(1)%p => sysvector

    write_scalar = .false.

    select case ( dataname )
    case('displacement')

      call create ( problem, dvec, vec=2 )

!     derive
      call derive_vector ( mesh, problem, dvec, elemsub=beam_displacement2, &
        coefficients=coefficients, mcoefficients=mcoefficients, &
        oldvectors=oldvectors, order='ND' )

    case('rotation')

      call create ( problem, dvec, vec=1 )

!     derive
      call derive_vector ( mesh, problem, dvec, elemsub=beam_rotation2, &
        coefficients=coefficients, mcoefficients=mcoefficients, &
        oldvectors=oldvectors, order='ND' )

      write_scalar = .true.

    case('distributed force')

      call create ( problem, dvec, vec=2, elementwise=.true. )

!     derive
      call derive_vector ( mesh, problem, dvec, elemsub=beam_q2, &
        coefficients=coefficients, mcoefficients=mcoefficients, &
        oldvectors=oldvectors, order='ND' )

    case default

      write(*,'(/a,a/)') 'Error: wrong value dataname: ', dataname
      stop

    end select

!   write

    if ( write_scalar ) then
      call write_scalar_gmsh ( mesh, problem, filename=mshfilename, &
        dataname=dataname, vector=dvec, append=append )
    else
      call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
        dataname=dataname, vector=dvec, append=append )
    end if

    select case ( dataname )
    case('displacement')
      filename = 'disp.out'
    case('rotation')
      filename = 'rota.out'
    case('distributed force')
      filename = 'q.out'
    case default
      write(*,'(/a,a/)') 'Error: wrong value dataname: ', dataname
      stop
    end select

!   print

    call printtofile ( mesh, problem, filename=filename, vector=dvec )

    call delete ( dvec )
    call delete ( oldvectors )

  end subroutine derive_write_vector_beam2

end module structures_post_beam2_m

