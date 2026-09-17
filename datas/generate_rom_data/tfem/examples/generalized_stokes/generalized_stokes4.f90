! Generalized Stokes problem on a unit square. Periodical boundary conditions.
! Poiseuille flow. Collocation of nodes.
! Varying visosity eta using a separate element to fill the Gauss
! point values of eta according to a regularized Bingham model.
! Picard iteration.

module subs4_m

  use generalized_stokes_elements_m

  implicit none

contains

! Internal element routine to fill the position dependent coefficient eta
! in the Gauss/nodal points. This element should be used together with the routine
! loop_over_elements.

  subroutine fill_eta ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use generalized_stokes_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors

    integer :: ip, j
    real(dp) :: gammadot, tauy, eta, r


!   set globals, gauss, shapefunctions, ...

    call set_stokes_elem ( mesh, problem, elgrp, elem, first, last, &
      coefficients, oldvectors )

!   allocate more arrays

    if ( first .and. coefficients%i(255) == 1 ) then

!     first element in this group and nodal points

      nump = nodalp

      call deallocate_arrays_nodalp

      call allocate_arrays_nodalp

      call refcoor_nodal_points ( mesh, elgrp, xig ) ! use xig for xrnod

      call set_shape_function ( shapefunc, xig, phi, dphi )

    else if ( first .and. coefficients%i(37) == 0 ) then

!     first element in this group and the same Gauss rule for all elements

      nump = ninti

      call allocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      nump = ninti

      call allocate_arrays

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_deformation ( x(1:ndf,:), dphi, F, Finv, detF )

    call shape_derivative ( dphi, Finv, dphidx )

    call isoparametric_coordinates ( x, phi, xg )

!   get velocity vector

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, elem, u, &
      physq=[physqvel], layer=layer )

    uvector = reshape ( u, [ndf,ndim] )

    ugvector =  matmul ( phi, uvector )

    do j = 1, ndim
      gradu(:,:,j) = matmul ( dphidx(:,:,j), uvector )
    end do

!   determine viscosity coefficient eta

    gammadot = 0

    do ip = 1, nump

      gammadot = sum ( ( gradu(ip,:,:) + transpose(gradu(ip,:,:)) ) ** 2 ) / 2

      if ( coorsys == 1 ) then
!       axisymmetric
        if ( xg(ip,2) < 1e-10_dp ) then
!         use du_r / dr (r=0) for gradu_tt
          gammadot = gammadot + 2 * ( gradu(ip,2,2) ) ** 2
        else
!         use u_r / r for gradu_tt
          gammadot = gammadot + 2 * ( ugvector(ip,2) / xg(ip,2) ) ** 2
        end if
      end if

      gammadot = sqrt ( gammadot )

!     viscosity function (regularized Bingham model)

      eta = coefficients%r(1)
      tauy = coefficients%r(209)
      r = coefficients%r(210)

      etag(ip) = eta + tauy * ( 1._dp - exp( - r * gammadot ) ) / gammadot

    end do

    call put_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, r1=etag )

!   unset globals, gauss, shapefunctions, ...

    call unset_stokes_elem ( last, coefficients )

!   deallocate more memory

    if ( last .and. coefficients%i(255) == 1 ) then

!     last element in this group and nodal points

      call deallocate_arrays

    else if ( last .and. coefficients%i(37) == 0 ) then

!     last element in this group and the same Gauss rule for all elements

      call deallocate_arrays

    else if ( coefficients%i(37) == 1 .or. coefficients%i(37) == 2 ) then

!     different Gauss rule for each element

      call deallocate_arrays

    end if

  contains

    subroutine allocate_arrays_nodalp

      allocate ( wg(nodalp), xig(nodalp,ndim), phi(nodalp,ndf), detF(nodalp) )
      allocate ( dphi(nodalp,ndf,ndim), F(nodalp,ndim,ndim) )
      allocate ( Finv(nodalp,ndim,ndim), dphidx(nodalp,ndf,ndim) )
      allocate ( xg(nodalp,ndim) )

      allocate ( u(ndim*ndf), uvector(ndf,ndim), gradu(nodalp,ndim,ndim) )
      allocate ( etag(nodalp), ugvector(nodalp,ndim) )

    end subroutine allocate_arrays_nodalp

    subroutine deallocate_arrays_nodalp

      deallocate ( wg, xig, phi, detF )
      deallocate ( dphi, F, Finv, dphidx )
      deallocate ( xg )

    end subroutine deallocate_arrays_nodalp

    subroutine allocate_arrays

      allocate ( u(ndim*ndf), uvector(ndf,ndim), gradu(ninti,ndim,ndim) )
      allocate ( etag(ninti), ugvector(ninti,ndim) )

    end subroutine allocate_arrays

    subroutine deallocate_arrays

      deallocate ( u, uvector, gradu )
      deallocate ( etag, ugvector )

    end subroutine deallocate_arrays

  end subroutine fill_eta

end module subs4_m

program generalized_stokes4

  use tfem_m
  use hsl_ma57_m
  use generalized_stokes_elements_m
  use io_utils_m
  use figplot_m
  use subs4_m

  implicit none

! constants

  integer, parameter :: &
    uintpl = 8,         & ! Q2 velocities
    pintpl = 4,         & ! Q1 pressures
    physqvel = 1,       & ! physical quantity nr of the velocities
    physqpress = 2,     & ! physical quantity nr of the pressures
    gauss = 3,          & ! 3x3 integration of quads
    itermax=100,        & ! maximum interations
    nx=20,              & ! number of elements in x
    ny=20                 ! number of elements in y

  real(dp), parameter :: &
    eta = 1._dp,        & ! viscosity of yielded material
    r = 100._dp,        & ! exponent parameter
    tauy = 10._dp,      & ! yield stress
    epsvd = 1.e-3_dp,   & ! maximum velocity difference for iteration
    epspd = 1.e-2_dp,   & ! maximum pressure difference for iteration
    flowrate = 2._dp

! definitions

  type(meshgen_options_t) :: meshgen_options
  type(mesh_t) :: mesh
  type(input_probdef_t) :: input_probdef
  type(problem_t) :: problem
  type(sysmatrix_t) :: sysmatrix
  type(sysvector_t), target :: sol, soln
  type(sysvector_t) :: rhsd
  type(vector_t) :: velocity, pressure, vorticity, stress, viscosity
  type(elvector_t), target :: elvector
  type(plot_options_t) :: plot_options
  type(oldvectors_t) :: oldvectors
  type(coefficients_t) :: coefficients
  type(subscript_t) :: vel, pres
  type(solver_options_ma57_t) :: solver_options

  integer :: iter
  real(dp) :: vd, pd

! fill coefficients

  call create_coefficients ( coefficients, ncoefi=300, ncoefr=250 )

  coefficients%i(1:11) = &
    [ uintpl,   pintpl,     0,     0,         0,  &
       physqvel, physqpress, 0,     0,     gauss,  &
       gauss ]
  coefficients%i(12:) = 0

  coefficients%i(251) = 3  ! compute eta using Gauss values in separate routine

  coefficients%r = 0
  coefficients%r(1) = eta     ! viscosity of yielded material
  coefficients%r(6) = flowrate

  coefficients%r(209) = tauy  ! yield stress
  coefficients%r(210) = r     ! exponent parameter

  call write_coefficients ( coefficients, filename='coefficients.out' )

! create mesh

  meshgen_options%elshape = 6
  meshgen_options%nx = nx
  meshgen_options%ny = ny

  call quadrilateral2d ( mesh, meshgen_options )

  call add_to_mesh ( mesh, curve=[-4] )     ! curve 5

! write mesh (read by streamfunction computation)

  call write_mesh ( mesh, filename='mesh.out' )

  call fill_mesh_parts ( mesh )

! problem definition

  call create_input_probdef ( mesh, input_probdef, nvec=3, nphysq=2 )

  input_probdef%vec_elementdof(1)%a =   &
      reshape ( [2,2,2,2,2,2,2,2,2,    &  ! velocity
                 1,0,1,0,1,0,1,0,0,    &  ! pressure
                 1,1,1,1,1,1,1,1,1 ], &  ! scalar, such as vorticity
                 [9,3] )

  input_probdef%physq = [1,2]

  call define_essential ( mesh, input_probdef, curve1=1, physq=1 )
  call define_essential ( mesh, input_probdef, curve1=3, physq=1 )
  call define_essential ( mesh, input_probdef, point=1, physq=2 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, nglobalc=1 )

  call define_constraint ( mesh, input_probdef, &
    physq=1, curve1=2, curve2=5, discretization='collocation', exclude=3 )

  call problem_definition ( input_probdef, mesh, problem )

! define vector subscripts for direct manipulation of sysvector data

! all velocities
  call create_subscript ( mesh, problem, vel, physqarr=[physqvel] )
! pressures
  call create_subscript ( mesh, problem, pres, physqarr=[physqpress] )

! create system vectors (solution and right-hand side)

  call create_sysvector ( problem, sol, soln )
  call create_sysvector ( problem, rhsd )

! fill solution vector with essential boundary conditions

  call fill_sysvector ( mesh, problem, sol, &
    curve1=1, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    curve1=3, physq=1, value=0._dp )
  call fill_sysvector ( mesh, problem, sol, &
    point=1, physq=2, value=0._dp )

  soln%u = 0 ! initialize old solution to zero

! Create vector defined per element and store Gauss values of eta

  call create_elvector ( mesh, elvector, nreal1d=1 )

  call create ( oldvectors, nsysvec=1, nelvec=1 )

  oldvectors%s(1)%p => soln
  oldvectors%e(1)%p => elvector

! create system matrix

  call create_sysmatrix_structure_base ( sysmatrix, mesh, problem, &
    symmetric=.true. )
  call create_sysmatrix_structure_constraint ( sysmatrix, mesh, problem )
  call finalize_sysmatrix_structure ( sysmatrix )

  call create_sysmatrix_data ( sysmatrix )

! Start Picard iteration

  iter = 0

  iterate: do

    iter = iter + 1

!   build (assemble) matrix and vector from elements

    if ( iter == 1 ) then

!     use Stokes for first iteration

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=stokes_elem, coefficients=coefficients, &
        oldvectors=oldvectors )

    else

!     store Gauss values of eta

      call loop_over_elements ( mesh, problem, elemsub=fill_eta, &
        coefficients=coefficients, oldvectors=oldvectors )

      call build_system ( mesh, problem, sysmatrix, rhsd, &
        elemsub=generalized_stokes_elem, coefficients=coefficients, &
        oldvectors=oldvectors )

    end if

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=1, elemsub=stokes_constr_flowr, addmatvec=.true., &
      coefficients=coefficients )

    call build_system_constraint ( mesh, problem, sysmatrix, rhsd, &
      constraint1=2, elemsub=stokes_constr_node_conn, addmatvec=.true., &
      coefficients=coefficients  )

    call add_effect_of_essential_to_rhs ( problem, sysmatrix, sol, rhsd )

    solver_options%real_storage=1.5

    call solve_system_ma57 ( sysmatrix, rhsd, sol, &
      solver_options = solver_options )

    vd = maxval(abs(sol%u(vel%s)-soln%u(vel%s)))  ! velocity difference
    pd = maxval(abs(sol%u(pres%s)-soln%u(pres%s)))  ! pressure difference

    print *, iter, vd , pd

    call copy ( sol, soln )

    if ( vd < epsvd .and. pd < epspd ) exit iterate

    if ( iter >= itermax ) then
      write(*,'(a,i0)') ' too many iterations: ', itermax
      stop
    end if

  end do iterate

! post-processing

  call create_vector ( problem, velocity, physq=1 )
  call create_vector ( problem, pressure, vec=3 )
  call create_vector ( problem, stress, vec=3 )
  call create_vector ( problem, vorticity, vec=3, elementwise=.true. )
  call create_vector ( problem, viscosity, vec=3 )

  call extract_physvector ( mesh, problem, sol, velocity )

! write profile to a file for plotting with gnuplot

  call printtofile ( mesh, problem, 'vprofile.out', curve=2, &
    vector=velocity )

! derive vectors

  oldvectors%s(1)%p => sol

  call derive_vector ( mesh, problem, pressure, elemsub=stokes_pressure, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=5

  call derive_vector ( mesh, problem, vorticity, elemsub=stokes_deriv, &
    coefficients=coefficients, oldvectors=oldvectors )

! store nodal values of eta in a vector defined per element

  coefficients%i(255)=1

  call loop_over_elements ( mesh, problem, elemsub=fill_eta, &
    coefficients=coefficients, oldvectors=oldvectors )

  call derive_vector ( mesh, problem, viscosity, &
    elemsub=generalized_stokes_viscosity, &
    coefficients=coefficients, oldvectors=oldvectors )

  coefficients%i(13)=2

  call derive_vector ( mesh, problem, stress, &
    elemsub=generalized_stokes_stress, &
    coefficients=coefficients, oldvectors=oldvectors )

! write to a fig file for plotting

  plot_options%scalevector=0.05
  call plot_vector ( plot_options, mesh, problem, 'velocity.fig', &
    vector=velocity )

  plot_options%numlevels=9
  call plot_color_contour ( plot_options, mesh, problem, &
    'pressure_contour.fig', vector=pressure )

  call plot_color_contour ( plot_options, mesh, problem, &
    'vorticity_contour.fig', vector=vorticity )

  call plot_color_contour ( plot_options, mesh, problem, &
    'viscosity_contour.fig', vector=viscosity )

! write profile to a file for plotting with gnuplot

  call printtofile ( mesh, problem, 'stressprofile.out', curve=2, &
    vector=stress )

! write binary file for reading by streamfunction

  open(unit=10,form='unformatted',file='velocity_bin.out')

  write(10) velocity%u

  close(unit=10)

! delete all data including all allocated memory

  call delete ( problem )
  call delete ( input_probdef )
  call delete ( mesh )
  call delete ( sol, rhsd )
  call delete ( velocity, pressure, vorticity, viscosity )
  call delete ( sysmatrix )
  call delete ( coefficients )
  call delete ( oldvectors )
  call delete ( vel, pres )
  call delete ( elvector )

end program generalized_stokes4
