! Elastic problem on a unit square (2D). Mooney-Rivlin solid.

program elastic1

  use tfem_m
  use hsl_ma41_m
  use elastic_elements_m
  use figplot_m
  use io_utils_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 displacements
    pintpl = 4,         & ! Q1 pressures
    physqdisp = 1,      & ! physical quantity nr of the displacements
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    model = 2,          & ! Mooney-Rivlin model
    nx=20,              & ! number of elements in x
    ny=10,              & ! number of elements in y
    numinc=3,           & ! number of increments
    itermax=10            ! maximum interations in an increment

  real(dp), parameter :: &
    C1 = 1.0, &
    C2 = 0.2

  integer :: i

  real(dp), parameter :: &
    ! incremental displacements
    xinc(numinc)=[(0.5_dp,i=1,numinc)], &
    ! stop criterium for iterative displacement differences
    epsdd(numinc)=[(1e-5_dp,i=1,numinc-1),1e-13_dp], &
    ! stop criterium for iterative pressure differences
    epspd(numinc)=[(1e-4_dp,i=1,numinc-1),1e-12_dp]

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t) :: sol, rhsd
  type(sysvector_t), target :: soln
  type(vector_t) :: displacement, pressure, energy, volumech
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(plot_options_t) :: plot_options
  type(subscript_t) :: disp, pres, ess

! variables

  integer :: iter, inc
  real(dp) :: dd, pd


! fill coefficients

  call create_coefficients ( coefficients, ncoefi=100, ncoefr=50 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,      model, 0,     0,  &
      physqdisp, physqpress,     0, 0, gauss,  &
      gauss ]
  coefficients%i(12:) = 0

  coefficients%r(1:2) = [ C1, C2 ]
  coefficients%r(3:) = 0

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! displacement
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar
                 [9,3] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=2, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=4, physq=1 )

  call problem_definition ( input_probdef, mesh, problem )

! define vector subscripts for direct manipulation of sysvector data

! all displacements
  call create_subscript ( mesh, problem, disp, physqarr=[1] )
! essential degrees
  call create_subscript ( mesh, problem, ess, physqarr=[1], &
                          unknownpart=.false. )
! pressures
  call create_subscript ( mesh, problem, pres, physqarr=[2] )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol )  ! Newton-Raphson iteration increment
  call create_sysvector ( problem, rhsd ) ! Right-hand side
  call create_sysvector ( problem, soln ) ! current solution

  sol%u = 0
  soln%u = 0

  call create ( oldvectors, nsysvec=1 )

  oldvectors%s(1)%p => soln ! used in element module

! create system matrix

  call create_sysmatrix_structure ( sysmatrix, mesh, problem )

  call create_sysmatrix_data ( sysmatrix )

! start increments

  increm: do inc = 1, numinc

    print *
    print *, 'increment ', inc
    print *

!   start Newton-Raphson iteration within increments

    iter = 0

    iterate: do

      iter = iter + 1

!     fill solution vector with essential boundary conditions
      if ( iter == 1 ) then
        call fill_sysvector ( mesh, problem, sol, &  ! increment
          curve1=2, physq=1, degfd=1, value=xinc(inc) )
      end if
      if ( iter == 2 ) sol%u(ess%s) = 0  ! set back essential bc to zero

!     build (assemble) matrix and vector from elements

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=elastic_elem, coefficients=coefficients, oldvectors=oldvectors )

      call check ( sysmatrix )

      call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

      call solve_system_ma41 ( sysmatrix, rhsd, sol )

      dd = maxval(abs(sol%u(disp%s))) ! displacement difference
      pd = maxval(abs(sol%u(pres%s)-soln%u(pres%s)))  ! pressure difference

      print *, iter, dd , pd

      soln%u(disp%s) = soln%u(disp%s) + sol%u(disp%s)  ! update solution
      soln%u(pres%s) = sol%u(pres%s) ! copy pressure

      if ( dd < epsdd(inc) .and. pd < epspd(inc) ) exit iterate

      if ( iter >= itermax ) then
        write(*,'(a,i0)') ' too many iterations: ', itermax
        stop
      end if

    end do iterate

  end do increm

! postprocessing, create and obtain vectors

  call create_vector ( problem, energy, vec=3 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, volumech, vec=3 )
  call create_vector ( problem, displacement, physq=1 )

  call derive_vector ( mesh, problem, pressure, elemsub=elastic_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13) = 1
  call derive_vector ( mesh, problem, energy, elemsub=elastic_derive, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13) = 2
  call derive_vector ( mesh, problem, volumech, elemsub=elastic_derive, &
    coefficients=coefficients, oldvectors=oldvectors )

! plot mesh before deformation

  plot_options%figsize=150/(1+sum(xinc)) ! scale original mesh back
  plot_options%meshcolor=4
  call plot_mesh ( plot_options, mesh, 'mesh.fig' )

! update coordinates of mesh

  call extract_physvector ( mesh, problem, soln, displacement )

  mesh%coor = mesh%coor + &
           transpose ( reshape (displacement%u,[2,mesh%nnodes]) )

! plot mesh after deformation

  plot_options%figsize=150
  plot_options%meshcolor=0
  call plot_mesh ( plot_options, mesh, 'mesh.fig', append=.true. )

! write to a fig file for plotting of pressure and energy

  plot_options%shiftbary=40
  call plot_color_fill ( plot_options, mesh, problem, 'pressure_color.fig', &
    vector=pressure )
  call plot_color_fill ( plot_options, mesh, problem, 'energy_color.fig', &
    vector=energy )
  call plot_color_fill ( plot_options, mesh, problem, 'J.fig', &
    vector=volumech )

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, soln, rhsd )
  call delete ( sysmatrix )
  call delete ( displacement, pressure, energy, volumech )
  call delete ( disp, pres, ess )
  call delete ( coefficients )
  call delete ( oldvectors )

end program elastic1
