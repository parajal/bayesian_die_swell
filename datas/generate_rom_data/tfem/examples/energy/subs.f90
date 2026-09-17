module subs_m

  use tfem_elem_m

  implicit none

  real(dp) :: leta0, lt0, lt1

  real(dp), allocatable, dimension(:) :: lG0, llam0, modg, lambdag

contains


! compute temperature-dependent viscosity

  function compute_viscosity ( t )

    real(dp), intent(in) :: t(:)
    real(dp) :: compute_viscosity(size(t))

    compute_viscosity = leta0 * exp(-(t-lt0))

  end function compute_viscosity


! compute temperature-dependent modulus

  function compute_modulus ( t, mode )

    real(dp), intent(in) :: t(:)
    integer, intent(in) :: mode
    real(dp) :: compute_modulus(size(t))

    compute_modulus = lG0(mode) * exp(-(t-lt0))

  end function compute_modulus


! compute temperature-dependent modulus

  function compute_relaxation_time ( t, mode )

    real(dp), intent(in) :: t(:)
    integer, intent(in) :: mode
    real(dp) :: compute_relaxation_time(size(t))

    compute_relaxation_time = llam0(mode) * exp(-(t-lt0))

  end function compute_relaxation_time


! Internal element routine to fill the temperature dependent coefficient eta
! in the Gauss points. This element should be used together  with the routine
! loop_over_elements.

  subroutine fill_viscosity_gauss ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use generalized_stokes_elements_m
    use generalized_stokes_globals_m
    use convection_diffusion_supg_elements_m
    use convection_diffusion_supg_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors

    if ( first ) then

!     first element in this group
!     set globals

      call set_stokes_elem ( mesh, problem, elgrp, elem, &
        first, last, coefficients, oldvectors )

      call set_conv_diff_supg_elem ( mesh, problem, elgrp, elem, &
        first, last, coefficients, oldvectors )

      allocate ( q_n(ndfq), q_ng(ninti), etag(ninti) )


    end if

    call get_coordinates ( mesh, elgrp, elem, x )

!   determine coefficient by a function of t

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(5)%p, elgrp, &
      elem, q_n, physq=[1] )

    q_ng = matmul ( chi, q_n )

!   function to determine the local viscosity

    etag = compute_viscosity ( q_ng )

    call put_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, r1=etag )

    if ( last ) then

!     last element in this group

      call unset_stokes_elem ( last, coefficients )
      call unset_conv_diff_supg_elem ( last, coefficients )

      deallocate ( q_n, q_ng, etag )

    end if

  end subroutine fill_viscosity_gauss


! Internal element routine to fill the temperature dependent modulus
! and relaxation time in the Gauss points. This element should be used
! together  with the routine loop_over_elements.

  subroutine fill_modlam_gauss ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use stokes_elements_m
    use stokes_globals_m
    use convection_diffusion_supg_elements_m
    use convection_diffusion_supg_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors

    integer :: nmodes, m

    if ( first ) then

!     first element in this group
!     set globals

      call set_stokes_elem ( mesh, problem, elgrp, elem, &
        first, last, coefficients, oldvectors )

      call set_conv_diff_supg_elem ( mesh, problem, elgrp, elem, &
        first, last, coefficients, oldvectors )

      allocate ( q_n(ndfq), q_ng(ninti) )
      allocate ( modg(ninti), lambdag(ninti) )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

!   determine coefficient by a function of t

    call get_sysvector ( mesh, oldvectors%p(2)%p, oldvectors%s(2)%p, elgrp, &
      elem, q_n, physq=[1] )

    q_ng = matmul ( chi, q_n )

!   functions to determine the local modulus and relaxation time

    nmodes = coefficients%i(19)
    do m = 1, nmodes
      modg = compute_modulus ( q_ng, m )
      call put_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, r1=modg, nr=m )
      lambdag = compute_relaxation_time ( q_ng, m )
      call put_elvector ( mesh, oldvectors%e(2)%p, elgrp, elem, r1=lambdag, &
        nr=m )
    end do

    if ( last ) then

!     last element in this group

      call unset_stokes_elem ( last, coefficients )
      call unset_conv_diff_supg_elem ( last, coefficients )

      deallocate ( q_n, q_ng )
      deallocate ( modg, lambdag )

    end if

  end subroutine fill_modlam_gauss


! Internal element routine to fill the temperature dependent modulus
! and in the nodes. The result is stored in an array per element (elvector).
! This element should be used together with the routine loop_over_elements.

  subroutine fill_modulus_nodes ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use stokes_elements_m
    use stokes_globals_m
    use convection_diffusion_supg_elements_m
    use convection_diffusion_supg_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors

    integer :: nmodes, m

    if ( first ) then

!     first element in this group

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )

      call set_globals_conv_diff_supg ( coefficients )

      allocate ( chi(nodalp,ndfq), x(nodalp,ndim) )
      allocate ( xrnod(nodalp,ndim) )
      allocate ( q_n(ndfq), q_ng(nodalp) )
      allocate ( modg(nodalp) )

      call refcoor_nodal_points ( mesh, elgrp, xrnod )

      call set_shape_function ( shapefuncq, xrnod, chi )

    end if

!   get temperature solution

    call get_sysvector ( mesh, problem, oldvectors%s(1)%p, elgrp, &
      elem, q_n, physq=[1] )

    q_ng = matmul ( chi, q_n )

!   functions to determine the local modulus and relaxation time

    nmodes = coefficients%i(19)
    do m = 1, nmodes
      modg = compute_modulus ( q_ng, m )
      call put_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, r1=modg, nr=m )
    end do

    if ( last ) then

!     last element in this group

      deallocate ( chi, x )
      deallocate ( xrnod )
      deallocate ( q_n, q_ng )
      deallocate ( modg )

    end if

  end subroutine fill_modulus_nodes

end module subs_m
