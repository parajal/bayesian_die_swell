module subs_m

  use kind_defs_m
  use tfem_elem_m

  implicit none

  real(dp) :: eta1l, eta2l

contains


! Internal element routine to fill the viscosity in the Gauss points as a
! function of the composition c

  subroutine fill_eta_gauss_c ( mesh, problem, elgrp, elem, first, last, &
    coefficients, oldvectors )

    use generalized_stokes_elements_m
    use generalized_stokes_globals_m
    use diffuse_interface_globals_m

    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem
    integer, intent(in) :: elgrp, elem
    logical, intent(in) :: first, last
    type(coefficients_t), intent(in) :: coefficients
    type(oldvectors_t), intent(inout) :: oldvectors


    if ( first ) then

!     first element in this group

!     set globals

      call set_globals_stokes_vp ( mesh, coefficients, elgrp )
      call set_globals_generalized_stokes ( mesh, coefficients, elgrp )

      allocate ( wg(ninti), xig(ninti,ndim), phi(ninti,ndf), x(nodalp,ndim) )
      allocate ( xg(ninti,ndim), etag(ninti) )
      allocate ( c_n(ndf), c_ng(ninti) )

!     set Gauss integration and shape function

      call set_Gauss_integration ( gauss, xig, wg )

      call set_shape_function ( shapefunc, xig, phi )

    end if

    call get_coordinates ( mesh, elgrp, elem, x )

    call isoparametric_coordinates ( x, phi, xg )

!   determine coefficient by a function of c

    call get_vector ( mesh, oldvectors%p(2)%p, oldvectors%v(1)%p, elgrp, &
      elem, c_n, layer=layer )

    c_ng = matmul ( phi, c_n )

!   function to determine the local viscosity

    etag = eta1l * ((c_ng+1.0_dp)/2.0_dp) - eta2l * ((c_ng-1.0_dp)/2.0_dp)

    call put_elvector ( mesh, oldvectors%e(1)%p, elgrp, elem, r1=etag )

    if ( last ) then

!     last element in this group

      deallocate ( wg, xig, phi, x )
      deallocate ( xg, etag )
      deallocate ( c_n, c_ng )

    end if

  end subroutine fill_eta_gauss_c


! Subroutine to find c=0 on a line by sampling c at the beginning of the line,
! halfway the line and at the end of the line. If c changes sign from the first
! to the second point, it is assumed that c=0 is located in the first segment,
! else it is assumed that c=0 is located in the second segment. Next, the chosen
! segment is used as a new line and the process is repeated until the segment
! has decreased in length below some tolerance ( default 1.e-12 ).
! NOTE: this routine assumes that the points c(x1) and c(x2) have opposite signs
! and there is only one position on the line where c changes sign

  subroutine find_c0 ( mesh, problemdi, soldi, x1, x2, x_c0, tol )

    type(sysvector_t), intent(in) :: soldi
    type(problem_t), intent(in) :: problemdi
    type(mesh_t), intent(in) :: mesh
    real(dp), dimension(:), intent(in) :: x1, x2
    real(dp), dimension(:), intent(out) :: x_c0
    real(dp), intent(in), optional :: tol

    real(dp), dimension(3,size(x1)) :: x
    real(dp), dimension(1,size(x1)) :: x_temp
    real(dp), dimension(3,size(x1)) :: xref
    real(dp), dimension(1,size(x1)) :: xref_temp
    integer, dimension(3,2) :: grpelm
    integer, dimension(1,2) :: grpelm_temp

    real(dp) :: c(3), ltol
    character(len=13) :: globalshape

    real(dp), allocatable :: cnod(:), phi(:,:), phi_temp(:,:)

    integer :: i,iter

    if ( .not. present(tol) ) ltol=1.e-12
    if ( present(tol) ) ltol=tol

    ! default elgrp is 1
    allocate ( cnod(mesh%element(1)%numnod) )
    allocate ( phi(3,mesh%element(1)%numnod) )
    allocate ( phi_temp(1,mesh%element(1)%numnod) )
    globalshape=mesh%element(1)%globalshape

!   first iteration
    x(1,:) = x1
    x(2,:) = (x1+x2)/2
    x(3,:) = x2

    iter=0

    finding_c0: do

      iter = iter+1

!     sample c in the three points
      if ( iter == 1 ) then

        call find_refcoor_points ( mesh, coor=x, grpelm=grpelm, refcoor=xref )

        if ( any(grpelm(:,1)==0 ) ) then ! sample point outside mesh
          write(*,'(/a/)') 'Error find_c0: sample point outside mesh.'
          stop
        end if

        call set_shape_function_simple ( globalshape, xr=xref, phi=phi )

        do i = 1,size(x,1)

          call get_sysvector ( mesh, problemdi, soldi, &
            elgrp=grpelm(i,1), elem=grpelm(i,2), u=cnod, physq=[1] )

          c(i) = dot_product ( phi(i,:), cnod )

        end do

      else  ! only update the middle point

        x_temp(1,:) = x(2,:)

        call find_refcoor_points ( mesh, coor=x_temp, grpelm=grpelm_temp, &
          refcoor=xref_temp )

        if ( any(grpelm(:,1)==0 ) ) then ! sample point outside mesh
          write(*,'(/a/)') 'Error find_c0: sample point outside mesh.'
          stop
        end if

        call set_shape_function_simple ( globalshape, xr=xref_temp, &
          phi=phi_temp )

        call get_sysvector ( mesh, problemdi, soldi, &
          elgrp=grpelm_temp(1,1), elem=grpelm_temp(1,2), u=cnod, physq=[1] )

        c(2) = dot_product ( phi_temp(1,:), cnod )

      end if


!     adapt x1 or x2 depending where c changes sign
      if ( c(1)*c(2) < 0.0_dp ) then
!       c changes sign in left segment
        x(3,:) = x(2,:)
        c(3) = c(2)
        x(2,:) = (x(1,:)+x(3,:))/2
      else
!       c changes sign in right segment
        x(1,:) = x(2,:)
        c(1) = c(2)
        x(2,:) = (x(1,:)+x(3,:))/2
      end if

!     check if search segment has decreased in length below the tolerance
      if ( iter > 1 ) then
        if ( sqrt(dot_product(x(3,:)-x(1,:),x(3,:)-x(1,:))) < ltol ) then
          x_c0 = x(2,:)
          exit finding_c0
        end if
      end if

      if ( iter > 100 ) then
        print *,'Warning: location of c=0 not found within the given &
          &tolerance in 100 iterations'
        x_c0 = x(2,:)
        exit finding_c0
      end if

    end do finding_c0

    deallocate ( cnod )
    deallocate ( phi )
    deallocate ( phi_temp )

  end subroutine find_c0

end module subs_m
