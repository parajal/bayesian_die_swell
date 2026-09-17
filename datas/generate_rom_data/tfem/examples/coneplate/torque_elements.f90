module torque_elements_m

    use tfem_elem_m
    use stokes_set_globals_m
    use viscoelastic_elements_generic_m

    implicit none

contains

! oldvectors:
!  v(1) = "viscous stress_tensor" derived from subroutine "stokes_stress_tensor"
!  v(2) = "pressure" derived from subroutine "stokes_pressure"

  subroutine torque_from_stress ( mesh, problem, curve, elem, first, last, &
    coefficients, oldvectors, elemvec )

    use stokes_globals_m
    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:) :: elemvec

    real(dp), allocatable, dimension(:,:) :: t_theta
    real(dp), allocatable, dimension(:,:,:) :: tauten_tot, tauten_visc, &
                                             tauten_ve
    real(dp), dimension(:), allocatable :: radius, radius_g
    integer :: ip
    integer, save :: maxstress


    if ( first ) then

!     first element on this curve

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, curve=curve, &
        maxvel3D=1 )

!     number of stress tensor components

      if ( vel3D == 1 ) then
        maxstress = 6
      else
        maxstress = 3 + coorsys
      end if

!     allocate arrays

      allocate ( wg(ninti), curvel(ninti) )
      allocate ( normal(ninti,ndim) )
      allocate ( xig(ninti,1), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim) )
      allocate ( dphi(ninti,ndf,1), dxdxi(ninti,ndim) )
      allocate ( st(ndf*maxstress), pr(ndf) )
      allocate ( work2(ndf,maxstress), work(ninti), work4(ninti,ncompu) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    allocate ( t_theta(ninti,1))
    allocate ( tauten_tot(ninti,ncompu,ncompu) )
    allocate ( tauten_visc(ninti,ncompu,ncompu) )
    allocate ( tauten_ve(ninti,ncompu,ncompu) )
    allocate ( radius(ndf), radius_g(ninti) )

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )

    call isoparametric_deformation_curve ( x, dphi(:,:,1), dxdxi, curvel, &
        normal )

    call isoparametric_coordinates ( x, phi, xg )

    if ( coorsys == 1 ) then
      curvel = 2 * pi * xg(:,2) * curvel
    end if

!   get reaction forces

    call get_vector_geometry ( mesh, oldvectors%p(1)%p, oldvectors%v(1)%p, &
      elem, st, curve=curve, layer=layer )

    work2 = reshape ( st, [ndf,maxstress] )

    ! initialize
    tauten_visc = 0._dp

    if ( vel3D == 1 ) then
      tauten_visc(:,1,1) = matmul ( phi, work2(:,1) )
      tauten_visc(:,1,2) = matmul ( phi, work2(:,2) )
      tauten_visc(:,1,3) = matmul ( phi, work2(:,3) )
      tauten_visc(:,2,1) = tauten_visc(:,1,2)
      tauten_visc(:,2,2) = matmul ( phi, work2(:,4) )
      tauten_visc(:,2,3) = matmul ( phi, work2(:,5) )
      tauten_visc(:,3,1) = tauten_visc(:,1,3)
      tauten_visc(:,3,2) = tauten_visc(:,2,3)
      tauten_visc(:,3,3) = matmul ( phi, work2(:,6) )
    else
      tauten_visc(:,1,1) = matmul ( phi, work2(:,1) )
      tauten_visc(:,1,2) = matmul ( phi, work2(:,2) )
      tauten_visc(:,2,1) = tauten_visc(:,1,2)
      tauten_visc(:,2,2) = matmul ( phi, work2(:,3) )
    end if

!   pressure

    call get_vector_geometry ( mesh, oldvectors%p(1)%p, oldvectors%v(2)%p, &
      elem, pr, curve=curve, layer=layer )

    work = matmul ( phi, pr )

!   Add pressure to stress tensor

    tauten_visc(:,1,1) = tauten_visc(:,1,1) - work
    tauten_visc(:,2,2) = tauten_visc(:,2,2) - work

    if ( vel3D == 1 ) tauten_visc(:,3,3) = tauten_visc(:,3,3) - work

    if ( coefficients%i(599) == 1 ) then

      call  calc_ve_stress ( mesh, oldvectors%p(2)%p, curve, elem, first, &
        last, coefficients, oldvectors, tauten_ve )

      ! initialize
      tauten_tot = 0._dp

!     calculate total stress tensor

      do ip = 1, ninti
        tauten_tot(ip,:,:) = tauten_visc(ip,:,:) + tauten_ve(ip,:,:)
      end do

    else

!     calculate total stress tensor

      do ip = 1, ninti
        tauten_tot(ip,:,:) = tauten_visc(ip,:,:)
      end do

    end if

!   traction force in each integration point: {t}=[tau]{n}

    do ip = 1, ninti
      work4(ip,:) = matmul( tauten_tot(ip,:,1:ndim), normal(ip,:) )
    end do

!   get radius
    call get_vector_geometry ( mesh, problem, oldvectors%v(3)%p, elem, radius, &
      curve=curve, layer=layer )

    radius_g = matmul ( phi, radius )

!   t_theta*r
    do ip = 1,ninti
       t_theta(ip,1) = work4(ip,3)*abs(radius_g(ip))
    end do

!   integrate traction to force

    elemvec = sum ( t_theta(:,1) * curvel * wg  )

    deallocate ( t_theta )
    deallocate ( tauten_tot, tauten_visc, tauten_ve )
    deallocate ( radius, radius_g )

    if ( last ) then

!     last element on this curve

      deallocate ( wg, curvel )
      deallocate ( normal )
      deallocate ( xig, phi, x )
      deallocate ( xg )
      deallocate ( dphi, dxdxi )
      deallocate ( st, pr )
      deallocate ( work2, work, work4 )

    end if

  end subroutine torque_from_stress

! oldvectors:
!  v1(1)%p(1:nmodes) = conformation tensors derived from subroutine
!                "derive_conformation_tensor_std" for _both_ the standard
!                and the log conformation. Thus in the log formulation
!                they contain the nodal values of logc.

  subroutine calc_ve_stress ( mesh, problem, curve, elem, first, &
    last, coefficients, oldvectors, elemtens )

    use viscoelastic_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: curve, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(in) :: oldvectors
    real(dp), intent(out), dimension(:,:,:) :: elemtens
    real(dp), allocatable, dimension(:,:,:) :: tauten_ve
    real(dp), allocatable, dimension(:) :: ct

    integer :: ip, m


    if ( first ) then

!     first element on this curve/surface

!     set globals

      call set_globals_stokes_vp_boun ( mesh, coefficients, curve=curve, &
        maxvel3D=1 )
      call set_viscoelastic_model ( coefficients )
      call set_globals_viscoelastic_boun ( coefficients )

!     allocate arrays

      allocate ( work8(ndf*ncomp), c(ndf,ncomp), cg(ninti,ncomp,nmodes) )
      allocate ( tauten(ninti,ncompu,ncompu), tauvec(ninti,ncompt) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi, dphi )

    end if

    allocate ( tauten_ve(ninti,ncompu, ncompu) )
    allocate ( ct(ndf*ncomp) )

!   geometry

    call get_coordinates_geometry ( mesh, elem, x, curve=curve )

    call isoparametric_deformation_curve ( x, dphi(:,:,1), dxdxi, curvel, &
      normal )

    if ( coorsys == 1 ) then
      call isoparametric_coordinates ( x, phi, xg )
      curvel = 2 * pi * xg(:,2) * curvel
    end if
!   conformation tensor

    do m = mode1, mode2

      call get_vector_geometry ( mesh, oldvectors%p(2)%p, &
        oldvectors%v1(1)%p(m), elem, work8, curve=curve, layer=layer )

      c = reshape ( work8, [ndf,ncomp] )

      cg(:,:,m) = matmul ( phi, c )

    end do

!   stress tensor: use of global variables; result in tauten

    call stress_tensor_viscoelastic ( cg )

!   traction force in each integration point: {t}=[tau]{n}

    do ip=1,ninti
      elemtens(ip,:,:) = tauten(ip,:,:)
    end do

    deallocate ( ct )

    if ( last ) then

!     last element on this curve

      call delete ( vemodel )

      deallocate ( work8, c, cg )
      deallocate ( tauten, tauvec )

    end if

  end subroutine calc_ve_stress

end module torque_elements_m
