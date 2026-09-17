
! renumbering using Metis version 5.

module metis5_m

  use kind_defs_m

  implicit none

  interface

    integer function METIS_SetDefaultOptions (options) &
      bind(C, name="METIS_SetDefaultOptions")

      !use iso_c_binding

      implicit none

      integer, dimension(0:39) :: options

    end function METIS_SetDefaultOptions

    integer function METIS_NodeND (nvtxs, xadj, adjncy, vwgt, options,  &
      perm, iperm) bind(C, name="METIS_NodeND")

      use iso_c_binding

      implicit none

      integer, intent(in) :: nvtxs
      integer, dimension(*) :: xadj, adjncy, perm, iperm
      integer, dimension(0:39) :: options
      type(c_ptr), value :: vwgt

    end function METIS_NodeND

  end interface

! interface for generic renumber_metis subroutine

  interface renumber_metis
    module procedure renumber_metis_mesh, renumber_metis_sysmatrix
  end interface renumber_metis

contains


! Renumber the nodes of the mesh using Metis.

  subroutine renumber_metis_mesh ( mesh )

    use mesh_m
    use iso_c_binding

    implicit none

    type(mesh_t), intent(inout) :: mesh

    integer :: metis_status
    integer, dimension(0:39) :: options
    integer, dimension(mesh%nnodes) :: iperm

    type(c_ptr) :: vwgt

    vwgt = C_NULL_PTR

!   fill_mesh_parts must have been called before

    call check_mesh ( mesh, 'renumber_metis_mesh' )

!   test for renumbering

    if ( mesh%renumber ) then
      write(*,'(/a/)') &
        'Warning in renumber_metis: mesh already contains renumbering.', &
        'Old renumbering will be removed.'
    else
      allocate ( mesh%nodperm(mesh%nnodes) )
    end if

!   set options

    metis_status = METIS_SetDefaultOptions (options)

    if ( metis_status /= 1 ) then
      write(*,'(/a/)') &
        'Error in METIS_SetDefaultOptions; call status = ', metis_status
      stop
    end if

!   options(METIS_OPTION_NUMBERING=17) = 1 for fortran numbering

    options(17) = 1

!   call metis

    metis_status = METIS_NodeND ( mesh%nnodes, mesh%nodnumnod + 1, &
      mesh%nodnod, vwgt, options, mesh%nodperm, iperm )

    if ( metis_status /= 1 ) then
      write(*,'(/a/)') &
        'Error in METIS_NodeND; call status = ', metis_status
      stop
    end if

    mesh%renumber = .true.

  end subroutine renumber_metis_mesh


! Renumber the unknown degrees of freedom of the system matrix using Metis.

  subroutine renumber_metis_sysmatrix ( sysmatrix, check_matrix )

    use set_optional_m
    use system_defs_m
    use system_matrix_m, only: check
    use iso_c_binding

    implicit none

    type(sysmatrix_t), intent(inout) :: sysmatrix

!   call check(sysmatrix) first? Default=.true.
    logical, intent(in), optional :: check_matrix

    integer :: metis_status
    integer, dimension(0:39) :: options
    integer, allocatable, dimension(:) :: xadj, adjncy

    logical :: lcheck_matrix
    integer :: n, i, j, k, m

    type(c_ptr) :: vwgt

    vwgt = C_NULL_PTR

!   check matrix

    lcheck_matrix = set_optional ( variable=check_matrix, default=.true. )

    if ( lcheck_matrix ) call check ( sysmatrix )

    if ( sysmatrix%symmetric ) then
      write(*,'(/a/a/)') &
        'Error in renumber_metis_sysmatrix: ', &
        '  symmetric system matrix structure not allowed'
      stop
    end if

!   number of degrees of freedom

    n = sysmatrix%Suu%n

!   test for renumbering

    if ( sysmatrix%renumber ) then
      write(*,'(/a/a/)') &
        'Warning in renumber_metis_sysmatrix: system matrix already contains', &
        ' renumbering. Old renumbering will be removed.'
    else
      deallocate ( sysmatrix%permu, sysmatrix%ipermu )
      allocate ( sysmatrix%permu(n), sysmatrix%ipermu(n) )
    end if

!   set options

    metis_status = METIS_SetDefaultOptions (options)

    if ( metis_status /= 1 ) then
      write(*,'(/a/)') &
        'Error in METIS_SetDefaultOptions; call status = ', metis_status
      stop
    end if

!   options(METIS_OPTION_NUMBERING=17) = 1 for fortran numbering

    options(17) = 1

!   create arrays xadj, adjncy (CSR format, diagonal terms excluded)

    allocate ( xadj(n+1), adjncy(sysmatrix%Suu%nnz) )

    m = 0
    xadj(1) = 1

    do i = 1, n
      do k = sysmatrix%Suu%ia(i), sysmatrix%Suu%ia(i+1) - 1
        j = sysmatrix%Suu%ja(k)
        if ( i == j ) cycle  ! exclude diagonal
        m = m + 1
        adjncy(m) = j
      end do
      xadj(i+1) = m + 1    ! start new row
    end do

!   call metis

    metis_status = METIS_NodeND ( n, xadj, adjncy(1:m), vwgt, options, &
                                 sysmatrix%permu, sysmatrix%ipermu )

    if ( metis_status /= 1 ) then
      write(*,'(/a/)') &
        'Error in METIS_NodeND; call status = ', metis_status
      stop
    end if

    deallocate ( xadj, adjncy )

    sysmatrix%renumber = .true.

  end subroutine renumber_metis_sysmatrix

end module metis5_m
