#if PARTRACK

! Stokes problem on a unit square with Dirichlet boundary conditions.
! Periodic cavity flow.
! Tracking of an initially circular interface.

module subs_m

  use tfem_m
  use math_defs_m
  use stokes_elements_m, only: stokes_sample_velocity

  implicit none

  type(mesh_t), save :: mesh
  type(problem_t), save :: problem
  type(sample_t), save :: sample
  type(oldvectors_t), save :: oldvectors_sample
  type(coefficients_t), save :: coefficients_sample

  real(dp) :: rcirc = 1, xcen(2) = 0

contains

! objectscoor defines the coordinates of the objects

  subroutine objectscoor ( objectnr, coor )
    integer, intent(in) :: objectnr
    real(dp), dimension(:,:), intent(inout) :: coor

    integer :: np, i
    real(dp) :: p(size(coor,1))

    np = size(coor,1)
    p = [(2*pi/np*(i-1),i=1,np)]

    coor(:,1) = rcirc * sin(p) + xcen(1)
    coor(:,2) = rcirc * cos(p) + xcen(2)

  end subroutine objectscoor

  subroutine derivs ( x, y, dydx )

    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in) :: y
    real(dp), dimension(:), intent(out) :: dydx

    real(dp) :: r(2)

    mesh%objects(2)%coor(1,:) = y

    call find_refcoor_objects ( mesh, object1=2 )

!   get the velocity in the object coordinates
    call fill_sample ( mesh, problem, sample, object=2, &
      elemsub=stokes_sample_velocity, coefficients=coefficients_sample, &
      oldvectors=oldvectors_sample )

!   set dydx
    if ( sample%nnodes /= 1 ) then
!     sample is empty
      write(*,'(/a/a/a,2e20.6/a/)') 'Warning derivs: no data in sample.', &
        'Reason is probably that the point is outside the domain', &
        'Coordinates: ', y
!     give a dydx (velocity) directed towards the center of the cavity
!     very inaccurate: the step will be rejected and time step reduced
      r = y - 0.5_dp
      dydx = - r / sqrt( r(1)**2 + r(2)** 2 )
    else
      dydx = sample%u(1,:)
    end if

  end subroutine derivs

end module subs_m

program particle_tracking2

  use tfem_m
  use hsl_ma57_m
  use stokes_elements_m
  use io_utils_m
  use subs_m
  use particle_tracking_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    nx=20,              & ! number of elements in x
    ny=20                 ! number of elements in y

  integer, parameter :: &
    npnts = 20,       & ! initial number of points in the object
    npntsmax = 20000, & ! maximum number of points in the object
    nperiod = 20,     & ! number of periods
    nsubsteps = 5,    & ! number of substeps per period
    everystep = 1       ! plot every ... steps

  real(dp), parameter :: &
    eta = 1._dp,         & ! viscosity
    halfperiod = 1._dp,  & ! half of the period
    deltatime = halfperiod / nsubsteps ! time of a substep

  real(dp), parameter :: & ! odeint parameters
    eps = 1e-4_dp,   &
    hmax = 0.1_dp,   &
    hcmax = 0.05_dp, &
    alphac = 0.4_dp, &
    deltat1 = 0.05_dp

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(input_probdef_t) :: input_probdef
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol1, sol2
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity
  type(lu_ma57_t) :: lu
  type(coefficients_t) :: coefficients
  type(tracking_options_t) :: tracking_options
  type(pcurve_t) :: pcurve

  logical, parameter :: sblock_search = .true.

  integer :: i, j, k, nod, ip, periodnr, step, img, nbx, nby, node
  real(dp) :: xs(1,2), time1, time2
  real(dp) :: coor_init(npnts,2)
  character(len=20) :: filename

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=150, ncoefr=100 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
      physqvel, physqpress, 0,     0,     gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1) = eta
  coefficients%r(2:) = 0

  call write_coefficients ( coefficients, filename='coefficients.out' )

  coefficients_sample = coefficients

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

! a fixed circular object

  rcirc = 0.1_dp ! radius of the circular object
  xcen(:) = [0.5_dp,0.85_dp]  ! centre of the circular object
  call objectscoor ( 1, coor_init )

  call add_to_mesh ( mesh, object='coordinates', nnodes=npntsmax, &
    coor=coor_init )

! one object for sampling the velocity in a single point
! initial value is not used, since it is refilled in derivs

  xs(1,:) = [0.5_dp,0.85_dp]

  call add_to_mesh ( mesh, object='coordinates', nnodes=1, coor=xs )

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

! blocks/sblocks

  if ( sblock_search ) then

!   use structured blocks

    call add_to_mesh ( mesh, sblocks=[nx,ny] )

  else

!   optimize number of blocks

    nbx = nint(sqrt(real(nx)))
    nby = nint(sqrt(real(ny)))

    call add_to_mesh ( mesh, blocks=[nbx,nby] )

  end if

  call fill_mesh_parts ( mesh )


! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, curve2=4, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

  call problem_definition ( input_probdef, mesh, problem )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol1, sol2 )
  call create_sysvector ( problem, rhsd )

! fill solution vector sol1 with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol1, &
    curve1=1, curve2=4, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol1, &
    curve1=3, physq=1, degfd=1, exclude=3, value=1._dp )
  call fill_sysvector ( mesh, problem, sol1, &
    point=1, physq=2, value=0._dp )

! fill solution vector sol2 with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol2, &
    curve1=1, curve2=4, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol2, &
    curve1=1, physq=1, degfd=1, exclude=3, value=-1._dp )
  call fill_sysvector ( mesh, problem, sol2, &
    point=1, physq=2, value=0._dp )

 ! sample in one point
  call create_sample ( mesh, sample, ndegfd=2, object=2 )

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem, symmetric=.true. )

  call create_sysmatrix_data ( sysmatrix )

! build (assemble) matrix and vector from elements

  call build_system ( mesh, problem, sysmatrix, rhsd, elemsub=stokes_elem, &
    coefficients=coefficients )

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol1, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol1, lu )

! for this problem the right-hand side is zero, no need to call
! build_system again. Therefore set rshd=0.

  rhsd%u=0

  call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol2, rhsd )

  call solve_system_ma57 ( sysmatrix, rhsd, sol2, lu )

  call delete ( lu )

! create oldvectors
  call create_oldvectors ( oldvectors_sample, nsysvec=1 )

! post-processing

  call create_pcurve ( pcurve, max_nnodes=npntsmax, closed=.true., &
    coor=coor_init )


  open ( unit=11, file='length.txt', recl=300 )

  write ( 11, * ) 0, area(pcurve), length(pcurve), pcurve%nnodes

  tracking_options%eps = eps
  tracking_options%hmax = hmax
  tracking_options%hcmax = hcmax
  tracking_options%alphac = alphac
  tracking_options%reorder = .true.

  time1 = 0
  time2 = 0
  img = 1

  open(unit=10,file='pcurvecoor001.txt')
  do node = 1,pcurve%nnodes
    write(10,'(es16.8,es16.8)') pcurve%coor(node,1),pcurve%coor(node,2)
  end do
  close(10)


  do periodnr = 1, nperiod

    oldvectors_sample%s(1)%p => sol1

    do step = 1, nsubsteps

      time1 = time2
      time2 = time1 + deltatime

      call step_points ( tracking_options, pcurve, time1, time2, deltat1, &
        derivs=derivs )

      call split_elements ( tracking_options, pcurve, deltat1, derivs=derivs )

      pcurve%prevcoor(1:pcurve%nnodes,:) = pcurve%coor(1:pcurve%nnodes,:)
      pcurve%prevtime = time2

      mesh%objects(1)%coor(1:pcurve%nnodes,:) = pcurve%coor(1:pcurve%nnodes,:)
      mesh%objects(1)%nnodes = pcurve%nnodes

      img = img + 1

      if ( mod(img,everystep) == 0 ) then

        write(filename,'(a,i3.3,a)') 'pcurvecoor', img/everystep, '.txt'
        open(unit=10,file=filename)
        do node = 1,pcurve%nnodes
          write(10,'(es16.8,es16.8)') pcurve%coor(node,1),pcurve%coor(node,2)
        end do
        close(10)

      end if

      write ( 11, * ) time2, area(pcurve), length(pcurve), pcurve%nnodes

      write(*,'(a,i0,a,i0)') 'period = ', periodnr, ' step = ', step

    end do

    oldvectors_sample%s(1)%p => sol2

    do step = nsubsteps+1, 2*nsubsteps

      time1 = time2
      time2 = time1 + deltatime

      call step_points ( tracking_options, pcurve, time1, time2, deltat1, &
        derivs=derivs )

      call split_elements ( tracking_options, pcurve, deltat1, derivs=derivs )

      pcurve%prevcoor(1:pcurve%nnodes,:) = pcurve%coor(1:pcurve%nnodes,:)
      pcurve%prevtime = time2

      mesh%objects(1)%coor(1:pcurve%nnodes,:) = pcurve%coor(1:pcurve%nnodes,:)
      mesh%objects(1)%nnodes = pcurve%nnodes

      img = img + 1

      if ( mod(img,everystep) == 0 ) then

        write(filename,'(a,i3.3,a)') 'pcurvecoor', img/everystep, '.txt'
        open(unit=10,file=filename)
        do node = 1,pcurve%nnodes
          write(10,'(es16.8,es16.8)') pcurve%coor(node,1),pcurve%coor(node,2)
        end do
        close(10)

      end if

      write ( 11, * ) time2, area(pcurve), length(pcurve), pcurve%nnodes

      write(*,'(a,i0,a,i0)') 'period = ', periodnr, ' step = ', step

    end do

  end do

  close ( unit=11 )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol1, sol2 )
  call delete ( rhsd )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors_sample )

end program particle_tracking2

#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on particle_tracking', &
    ' - add the nr library for linking', &
    ' - set preprocessing macro PARTRACK in Mdefs.mk'
end
#endif

