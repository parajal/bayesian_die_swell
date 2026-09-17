
! Post processing routines for truss structures

module structures_post_truss_m

  use tfem_m
  use structures_elements_m
  use io_utils_m

  implicit none

contains

! generic routine for post-processing all data for a truss

  subroutine postprocessing_truss ( mesh, problem, mshfilename, sol, load, &
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


!   direct write of sysvectors load, displacements, reacf

    if ( present(load) ) call write_sysvector_truss ( mesh, problem, &
      mshfilename, dataname='load', sysvector=load )

    call write_sysvector_truss ( mesh, problem, mshfilename, &
      dataname='displacement', sysvector=sol, append=present(load) )

    if ( present(reacf) ) call write_sysvector_truss ( mesh, problem, &
      mshfilename, dataname='reaction forces', sysvector=reacf, append=.true. )

!   write of distributed load via deriv

    call derive_write_vector_truss ( mesh, problem, mshfilename, &
      dataname='distributed force', sysvector=sol, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      append=.true. )

!   write of truss scalars

    call derive_write_scalar_truss ( mesh, problem, mshfilename, &
      dataname='strain', sysvector=sol, coefficients=coefficients, &
      mcoefficients=mcoefficients, append=.true. )

    call derive_write_scalar_truss ( mesh, problem, mshfilename, &
      dataname='stress', sysvector=sol, coefficients=coefficients, &
      mcoefficients=mcoefficients, append=.true. )

    call derive_write_scalar_truss ( mesh, problem, mshfilename, &
      dataname='force', sysvector=sol,coefficients=coefficients, &
      mcoefficients=mcoefficients,  append=.true. )

  end subroutine postprocessing_truss


! generic routine for post-processing system vectors for a truss

  subroutine write_sysvector_truss ( mesh, problem, mshfilename, dataname, &
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

!   write

    call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
      dataname=dataname, sysvector=sysvector, append=append )

!   print

    select case ( dataname )
    case('load')
      filename = 'load.out'
    case('displacement')
      filename = 'disp.out'
    case('reaction forces')
      filename = 'reacf.out'
    case default
      write(*,'(/a,a/)') 'Error: wrong value dataname: ', dataname
      stop
    end select

    call printtofile ( mesh, problem, filename=filename, sysvector=sysvector )

  end subroutine write_sysvector_truss


! generic routine for post-processing scalar data using derive for a truss

  subroutine derive_write_scalar_truss ( mesh, problem, mshfilename, &
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

!   derive

    if ( present(mcoefficients) ) then

      do i = 1, size(mcoefficients)
        select case ( dataname )
        case('strain')
          mcoefficients(i)%i(1) = 1 ! strain
        case('stress')
          mcoefficients(i)%i(1) = 2 ! stress
        case('force')
          mcoefficients(i)%i(1) = 3 ! force
        case default
          write(*,'(/a,a/)') 'Error: wrong value dataname: ', dataname
          stop
        end select
      end do

    else

      select case ( dataname )
      case('strain')
        coefficients%i(1) = 1 ! strain
      case('stress')
        coefficients%i(1) = 2 ! stress
      case('force')
        coefficients%i(1) = 3 ! force
      case default
        write(*,'(/a,a/)') 'Error: wrong value dataname: ', dataname
        stop
      end select

    end if

    call derive_vector ( mesh, problem, dvec, elemsub=truss_deriv1, &
      coefficients=coefficients, mcoefficients=mcoefficients, &
      oldvectors=oldvectors )

!   write

    call write_scalar_gmsh ( mesh, problem, filename=mshfilename, &
      dataname=dataname, vector=dvec, append=append )

!   print

    select case ( dataname )
    case('strain')
      filename = 'strain.out'
    case('stress')
      filename = 'stress.out'
    case('force')
      filename = 'force.out'
    case default
      write(*,'(/a,a/)') 'Error: wrong value dataname: ', dataname
      stop
    end select

    call printtofile ( mesh, problem, filename=filename, vector=dvec )

    call delete ( dvec )
    call delete ( oldvectors )

  end subroutine derive_write_scalar_truss


! post-processing vector data using derive for a truss

  subroutine derive_write_vector_truss ( mesh, problem, mshfilename, &
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
    type(vector_t) :: dvec
    type(oldvectors_t) :: oldvectors

    call create ( oldvectors, nsysvec=1 )
    oldvectors%s(1)%p => sysvector

    select case ( dataname )
    case('distributed force')

      call create ( problem, dvec, vec=2, elementwise=.true. )

!     derive
      call derive_vector ( mesh, problem, dvec, elemsub=truss_q1, &
        coefficients=coefficients, mcoefficients=mcoefficients, &
        oldvectors=oldvectors, order='ND' )

    case default

      write(*,'(/a,a/)') 'Error: wrong value dataname: ', dataname
      stop

    end select

!   write

    call write_vector_gmsh ( mesh, problem, filename=mshfilename, &
      dataname=dataname, vector=dvec, append=append )

    select case ( dataname )
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

  end subroutine derive_write_vector_truss

end module structures_post_truss_m

