
! Copyright (C) 2006-2018 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Types, routines, for solving the system of equations using MUMPS
!

module mumps_m

  use kind_defs_m
  use sparse_m
  use system_defs_m

  implicit none

  include 'dmumps_struc.h'
  include 'mpif.h'


! type definition for solver options

  type solver_options_mumps_t

    integer :: printlevel = 0 ! 0 no printing
                              ! 1 print info, rinfo
                              ! 2 print a lot of information from the solver

    integer :: ordering = 7  ! ICNTL(7): Ordering scheme used (see manual MUMPS)
                             ! ordering=7 means automatic choice by MUMPS.

    integer :: storage = 20  ! ICNTL(14): percentage increase in the estimated
                             ! working space.

    integer :: niter_ref = 0 ! ICNTL(10): maximum number of iterative
                             ! refinement steps

  end type solver_options_mumps_t


! interface for generic create subroutine

  interface create
    module procedure create_mumps_par
  end interface create


! interface for generic delete subroutine

  interface delete
    module procedure delete_mumps_par
  end interface delete

! interface for generic set_solver_options subroutine

  interface set_solver_options
    module procedure set_solver_options_mumps
  end interface set_solver_options



contains


! Helper routine for setting the solver options

  subroutine set_solver_options_mumps ( solver_options, keep, printlevel, &
    ordering, storage, niter_ref )

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_mumps_t), intent(inout) :: solver_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see solver_options_mumps_t
    integer, intent(in), optional :: printlevel, ordering, storage, niter_ref


    type(solver_options_mumps_t) :: solveropt


    if ( present(keep) ) then
      if ( .not. keep ) solver_options = solveropt
    else
      solver_options = solveropt
    end if

    if ( present(printlevel) ) solver_options%printlevel = printlevel
    if ( present(ordering) ) solver_options%ordering = ordering
    if ( present(storage) ) solver_options%storage = storage
    if ( present(niter_ref) ) solver_options%niter_ref = niter_ref

  end subroutine set_solver_options_mumps




! create an instance of mumps_par

  subroutine create_mumps_par ( mumps_par, symmetric, positive_definite )

!   The structure for MUMPS.
    type(dmumps_struc), intent(inout) :: mumps_par

!   Indicates whether matrix is symmetric. Default = .false. (not symmetric).
    logical, intent(in), optional :: symmetric

!   Indicates whether matrix is positive_definite. Default = .false. (not
!   positive_definite).
    logical, intent(in), optional :: positive_definite


    logical :: sym, pdef


    if ( present(symmetric) ) then
      sym = symmetric
      if ( present(positive_definite) ) then
        pdef = positive_definite
      else
        pdef = .false.
      end if
    else
      sym  = .false.
      pdef = .false.
    end if

!   define a communicator

    mumps_par%comm = MPI_COMM_WORLD

    mumps_par%job = -1  ! initialize

    if ( sym ) then
!     symmetric
      if ( pdef ) then
!       symmetric positive definite
        mumps_par%sym = 1
      else
!       general symmetric
        mumps_par%sym = 2
      end if
    else
!     unsymmetric
      mumps_par%sym = 0
    end if

    mumps_par%par = 1 ! host does factorization/solve

!   initialize an instance of MUMPS

    call dmumps ( mumps_par )

  end subroutine create_mumps_par


! Solve the system of equations using MUMPS (a parallel direct sparse solver)

  subroutine solve_system_mumps ( sysmatrix, rhsd, sol, mumps_par, &
     solver_options )

    type(sysmatrix_t), intent(in) :: sysmatrix

!   Right-hand side vector. The part corresponding to the unknowns is used as
!   the right-hand side in the system:
!
!   The effect of essential boundary conditions must already have been taken
!   into account in the rhsd vector
    type(sysvector_t), intent(in) :: rhsd

!   The solution vector. The unknown part sol_u is filled with the solution of
!   the system (*), given above. The prescribed part sol_p is not modified.
    type(sysvector_t), intent(inout) :: sol

!   The structure for MUMPS. It is used to store the LU factorization and
!   all other storage for the solver. It is the reponsibility of the user
!   to create mumps_par and delete mumps_par from memory.
!   The next time the solve routine is called with a mumps_par that still
!   contains a LU decomposition only backsubstition is done. This happens even
!   when the matrix has changed. If mumps_par has been cleared using delete,
!   a full factorization is performed.
    type(dmumps_struc), intent(inout) :: mumps_par

!   The options/parameters for the solver
!   See the type definition for the defaults.
    type(solver_options_mumps_t), intent(in), optional :: solver_options


    integer :: numundegfd
    type(solver_options_mumps_t) :: lsolveropt ! local options


!   initialize local solver options if solver_options in heading
    if ( present(solver_options) ) lsolveropt = solver_options


    if ( mumps_par%myid == 0 ) then

!     HOST

!     some testing first

      if ( .not. sysmatrix%initialized_structure ) then
        write(*,'(/a/)') &
          'Error in solve_system_mumps: no system matrix structure.'
        stop
      end if

      if ( .not. sysmatrix%finalized ) then
        write(*,'(/a/)') &
          'Error in solve_system_mumps: system matrix has not been finalized.'
        stop
      end if

      if ( .not. sysmatrix%allocated_data ) then
        write(*,'(/a/)') &
          'Error in solve_system_mumps: data in system matrix not allocated.'
        stop
      end if

      if ( sysmatrix%symmetric ) then
        if ( mumps_par%sym == 0 ) then
          write(*,'(/3(a/))') &
            'Error in solve_system_mumps: ', &
            ' for a symmetric sysmatrix the structure mumps_par needs ', &
            ' to be created with symmetric=.true.'
          stop
        end if
      else if ( .not. sysmatrix%symmetric ) then
        if ( mumps_par%sym /= 0 ) then
          write(*,'(/3(a/))') &
            'Error in solve_system_mumps: ', &
            ' a structure mumps_par created with symmetric=.true.', &
            ' requires the system matrix to be created with symmetric=.true.'
          stop
        end if
      end if

      if ( .not. rhsd%created ) then
        write(*,'(/a/)') &
          'Error in solve_system_mumps: rhsd not created.'
        stop
      end if

      if ( .not. sol%created ) then
        write(*,'(/a/)') &
          'Error in solve_system_mumps: sol not created.'
        stop
      end if

!     number of unknowns

      numundegfd = sysmatrix%Suu%n

!     copy first part of rhsd to sol

      sol%u(1:numundegfd) = rhsd%u(1:numundegfd)

!     solve, solution returned in sol(1:numundegfd)

      call solve_mumps1 ( sysmatrix%Suu, sol%u(1:numundegfd), mumps_par, &
        lsolveropt )

    else

!     SLAVE

      call solve_mumps2 ( mumps_par, lsolveropt )

    end if

  end subroutine solve_system_mumps


! solve linear system A u = b (host)

  subroutine solve_mumps1 ( A, b, mumps_par, so )

!   The matrix A in CSR format
    type(sparsematrix_t), intent(in), target :: A

!   On entry: the right-hand side vector b.
!   On exit: the solution u = A^{-1} b.
    real(dp), intent(inout), dimension(:) :: b

!   The structure for MUMPS. It is used to store the LU factorization. The next
!   time the solve routine is called with a mumps_par that contains a LU
!   decomposition only backsubstition is done. This happens even when the
!   matrix has changed. If mumps_par has been cleared using delete, a full
!   factorization is performed.
    type(dmumps_struc), intent(inout) :: mumps_par

!   The options/parameters for the solver
    type(solver_options_mumps_t), intent(in) :: so


    logical :: allsteps
    integer :: n, ne, row


!   initialize local parameters

    allsteps = .true.

!   set dimensions of problem

    n  = A%n
    ne = A%nnz

!   do all steps (analysis, factorization, solution) or only solution?

    if ( associated(mumps_par%s) ) then
!     decomposition is already done
      allsteps = .false.
    end if

!   choose job

    if ( allsteps ) then

!     create user arrays

      mumps_par%n = n
      mumps_par%nz = ne
      mumps_par%a => A%a
      mumps_par%jcn => A%ja
      allocate( mumps_par%irn ( ne ) )
      allocate( mumps_par%rhs ( n ) )

!     fill row number array to transform CSR -> coordinate format

      do row = 1, n
        mumps_par%irn( A%ia(row):A%ia(row+1)-1 ) = row
      end do

!     default values for the solver

      mumps_par%icntl(7) = so%ordering  ! Ordering

!     maximum number of iterative refinement steps

      mumps_par%icntl(10) = so%niter_ref

!     give user control over workspace

      mumps_par%icntl(14) = so%storage

      if ( so%printlevel == 2 ) then
!       do a lot of printing
        mumps_par%icntl(2) = 6
        mumps_par%icntl(3) = 6
        mumps_par%icntl(4) = 3
      else
        mumps_par%icntl(3) = 0  ! no information printing
      end if

!     analysis

      mumps_par%job = 1

      call dmumps ( mumps_par )

      call prininfo ( 'dmumps analysis', 16, 3, mumps_par%info, &
        mumps_par%rinfo, proc=0, printlevel=so%printlevel )

      mumps_par%job  = 5

    else

!     only solution step

      mumps_par%job  = 3

    end if

!   solve system

    mumps_par%rhs = b

    call dmumps ( mumps_par )

    b = mumps_par%rhs

    call prininfo ( 'dmumps solve', 16, 3, mumps_par%info, mumps_par%rinfo, &
      proc=0, printlevel=so%printlevel )

    call prininfo ( 'infog, rinfog', 25, 3, mumps_par%infog, mumps_par%rinfog, &
      proc=0, printlevel=so%printlevel )

  end subroutine solve_mumps1


! solve linear system A u = b (slaves)

  subroutine solve_mumps2 ( mumps_par, so )

!   The structure for MUMPS. It is used to store the LU factorization. The next
!   time the solve routine is called with a mumps_par that contains a LU
!   decomposition only backsubstition is done. This happens even when the
!   matrix has changed. If mumps_par has been cleared using delete, a full
!   factorization is performed.
    type(dmumps_struc), intent(inout) :: mumps_par

!   The options/parameters for the solver
    type(solver_options_mumps_t), intent(in) :: so


    logical :: allsteps


!   initialize local parameters

    allsteps = .true.

!   do all steps (analysis, factorization, solution) or only solution?

    if ( associated(mumps_par%s) ) then
!     decomposition is already done
      allsteps = .false.
    end if

!   choose job

    if ( allsteps ) then

      if ( so%printlevel == 2 ) then
!       do a lot of printing
        mumps_par%icntl(2) = 6
        mumps_par%icntl(3) = 6
        mumps_par%icntl(4) = 3
      else
        mumps_par%icntl(3) = 0  ! no information printing
      end if

!     maximum number of iterative refinement steps

      mumps_par%icntl(10) = so%niter_ref

!     analysis

      mumps_par%job = 1

      call dmumps ( mumps_par )

      call prininfo ( 'dmumps analysis', 16, 3, mumps_par%info, &
        mumps_par%rinfo, proc=mumps_par%myid, printlevel=so%printlevel )

      mumps_par%job  = 5

    else

!     only solution step

      mumps_par%job  = 3

    end if

!   solve system

    call dmumps ( mumps_par )

    call prininfo ( 'dmumps solve', 16, 3, mumps_par%info, mumps_par%rinfo, &
      proc=mumps_par%myid, printlevel=so%printlevel )

  end subroutine solve_mumps2


! prininfo

  subroutine prininfo ( sub, ni, nr, info, rinfo, proc, printlevel )

    character (len=*), intent(in) :: sub
    integer, intent(in) ::  ni, nr, printlevel
    integer, intent(in), dimension(:) :: info
    real(dp), intent(in), dimension(:) :: rinfo
    integer, intent(in) ::  proc

    integer :: i

    if ( info(1) < 0 ) then

!     Error in sub

      write(*,'(/a,i0,a/)') 'Output for processor ', proc, ':'

      write(*,'(/3a,i0/)') 'Error in ', sub, ' INFO(1) = ', info(1)

      call prinfo

      stop

    else if ( printlevel == 1 ) then

!     Print info, rinfo

      write(*,'(/a,i0,a/)') 'Output for processor ', proc, ':'

      write(*,'(/3a,i0/)') 'Info after subroutine ', sub

      call prinfo

    end if

  contains

    subroutine prinfo

      write(*,'(/a)') 'INFO:'

      do i = 1, ni
         write(*,'(i5,1x,i0)') i, info(i)
      end do

      write(*,'(/a)') 'RINFO:'

      do i = 1, nr
         write(*,*) i, rinfo(i)
      end do

      write(*,'()')

    end subroutine prinfo

  end subroutine prininfo


! delete an instance of mumps_par

  subroutine delete_mumps_par ( mumps_par )

!   The structure for MUMPS.
    type(dmumps_struc), intent(inout) :: mumps_par

!   delete user arrays

    mumps_par%n = 0
    mumps_par%nz = 0
    mumps_par%a => null()
    mumps_par%jcn => null()
    if ( associated(mumps_par%irn) ) deallocate( mumps_par%irn )
    if ( associated(mumps_par%rhs) ) deallocate( mumps_par%rhs )

!   delete an instance of MUMPS

    mumps_par%job = -2

    call dmumps ( mumps_par )

  end subroutine delete_mumps_par

end module mumps_m
