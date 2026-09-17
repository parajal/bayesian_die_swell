! A left built-in 2D beam problem with optional supports.
! Different distributed load on parts of the beam using three groups of
! elements with difference coefficients.
! Alternatively a function can be defined for q.
! Forces and/or moments in all points.
! P1/P3 interpolation. 2D elements.
! Similar to structures9 with the following extensions:
!  - axial loading
!  - optional angle with respect to the x-axis

program structures10

  use tfem_m
  use math_defs_m
  use hsl_ma57_m
  use structures_elements_m
  use structures_post_m
  use io_utils_m

  implicit none


! constants

  logical, parameter :: &
    builtin = .true., & ! builtin at x=0
    rotate = .true., &  ! rotate beam
!   horizontal and vertical supports at end of parts
    support_h(3) = [ .false., .false., .false. ], &
    support_v(3) = [ .false., .true., .false. ]

  integer, parameter :: &
    ninti = 3,  & ! number of Gauss points
    ne(3) = [ 10, 10, 10 ]  ! number of elements in each group

  real(dp), parameter :: &
    Emod = 1._dp, &  ! E modulus
    Ac = 1._dp,   &  ! cross-sectional area
    Iz = 1._dp,   &  ! second moment of area
    L(3) = [ 1._dp, 1._dp, 1._dp ], &  ! length of the beam parts
    angle = pi/4,   &  ! angle of rotation (in radians)
!   external forces/moments at the end of the parts (note: vertical positive up)
    F_h(3) = [ 0._dp, 1._dp, 1._dp ], &   ! horizontal external forces
    F_v(3) = [ 0._dp, 0._dp, -1._dp ], &  ! vertical external forces
    M(3) = [ 0._dp, 0._dp, 0._dp ], &     ! external moments
!   distributed loading on the beam parts ((note: vertical positive up)
    q_h(3) = [ 0._dp, 0._dp, 1._dp ], &  ! horizontal distributed loading
    q_v(3) = [ 0._dp, -1._dp, 0._dp ]    ! vertical distributed loading

! definitions

  type(mesh_t) :: mesh1, mesh2, mesh3, mesh4, mesh5, mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol
  type(sysvector_t) :: rhsd, reacf, load
  type(meshgen_options_t) :: mesh_options
  type(coefficients_t) :: coefficients(3)

  integer :: i, degfd(3), k
  real(dp) :: dir(2)


! fill coefficients

  call create_coefficients ( coefficients(1), ncoefi=100, ncoefr=50 )

  coefficients(1)%i = 0
  coefficients(1)%i(2) = ninti
  coefficients(1)%i(3) = 0 ! P3 Hermite interpolation for v
  coefficients(1)%i(4) = 0 ! constant distributed force
  coefficients(1)%i(12) = 2 ! P1 interpolation for u

  coefficients(1)%r = 0
  coefficients(1)%r(1) = Ac
  coefficients(1)%r(2) = Emod
  coefficients(1)%r(3) = Iz
  coefficients(1)%r(5:6) = [ q_h(1), q_v(1) ]

  coefficients(2) = coefficients(1)
  coefficients(2)%r(5:6) = [ q_h(2), q_v(2) ]

  coefficients(3) = coefficients(1)
  coefficients(3)%r(5:6) = [ q_h(3), q_v(3) ]


! create mesh

  mesh_options%elshape = 1 ! two-node line elements

  mesh_options%nx = ne(1)  ! number of elements, equidistant
  mesh_options%lx = L(1)   ! length of interval

  call line1d ( mesh1, mesh_options )

  mesh_options%nx = ne(2)  ! number of elements, equidistant
  mesh_options%lx = L(2)   ! length of interval
  mesh_options%ox = L(1)   ! origin of interval

  call line1d ( mesh2, mesh_options )

  call mesh_merge ( mesh1, mesh2, mesh3, point1=2, point2=1, &
    nogroupmerge=.true. )

  mesh_options%nx = ne(3)  ! number of elements, equidistant
  mesh_options%lx = L(3)   ! length of interval
  mesh_options%ox = L(1)+L(2) ! origin of interval

  call line1d ( mesh4, mesh_options )

  call mesh_merge ( mesh3, mesh4, mesh5, point1=3, point2=1, &
    nogroupmerge=.true. )

  call mesh_convert ( mesh5, mesh, coordinates=[1,0] )

  if ( rotate ) then
    dir = [ cos(angle), sin(angle) ]
    do k = 1, mesh%nnodes
      mesh%coor(k,:) = dir * mesh%coor(k,1)
    end do
  end if

  call fill_mesh_parts ( mesh )

  call delete ( mesh1, mesh2, mesh3, mesh4, mesh5 )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=2 )

  do i = 1, 3
    input_probdef%elementdof(i)%a = 3          ! (u,v,phi) in two end nodes
    input_probdef%vec_elementdof(i)%a = 1      ! scalar quantity
    input_probdef%vec_elementdof(i)%a(:,2) = 2 ! vector quantity
  end do

! at x=0

  if ( builtin ) then
    call define_essential ( mesh, input_probdef, point=1 )
  else
    call define_essential ( mesh, input_probdef, point=1, degfd=[1,1,0] )
  end if

! at end points of beam sections:

  do i = 1, 3
    if ( support_h(i) .or. support_v(i) ) then
      degfd = 0
      if ( support_h(i) ) degfd(1) = 1
      if ( support_v(i) ) degfd(2) = 1
      call define_essential ( mesh, input_probdef, point=i+1, degfd=degfd )
    end if
  end do

  call problem_definition ( input_probdef, mesh, problem )

  call printinfo ( mesh, printlevel=4 )
  call printinfo ( problem, printlevel=4 )


! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )
  call create_sysvector ( problem, rhsd, reacf, load )

! fill solution vector with essential boundary conditions

  if ( builtin ) then
    call fill_sysvector ( mesh, problem, sol, point=1, value=0._dp )
  else
    call fill_sysvector ( mesh, problem, sol, point=1, degfd=1, value=0._dp )
    call fill_sysvector ( mesh, problem, sol, point=1, degfd=2, value=0._dp )
  end if
  do i = 1, 3
    if ( support_h(i) .or. support_v(i) ) then
      if ( support_h(i) ) call fill_sysvector ( mesh, problem, sol, &
                                 point=i+1, degfd=1, value=0._dp )
      if ( support_v(i) ) call fill_sysvector ( mesh, problem, sol, &
                                 point=i+1, degfd=2, value=0._dp )
    end if
  end do

! fill rhs vector with discrete forces/monents at end points of beam sections:

  load%u = 0

  do i = 1, 3
    if ( .not. support_h(i) ) then
      call fill_sysvector ( mesh, problem, load, point=i+1, degfd=1, &
         value=F_h(i) )
    end if
    if ( .not. support_v(i) ) then
      call fill_sysvector ( mesh, problem, load, point=i+1, degfd=2, &
         value=F_v(i) )
    end if
    call fill_sysvector ( mesh, problem, load, point=i+1, degfd=3, value=M(i) )
  end do

  call copy ( load, rhsd )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=beam_elem2, &
    order='ND', mcoefficients=coefficients, addvec=.true. )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol )

  call reaction_forces ( problem, sysmatrix, sol, rhsd, reacf )

! post processing

  call postprocessing_beam2 ( mesh, problem, mcoefficients=coefficients, &
    mshfilename='structures10.mso', sol=sol, load=load, reacf=reacf )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd, reacf, load )
  call delete ( sysmatrix )
  call delete ( coefficients )

end program structures10

