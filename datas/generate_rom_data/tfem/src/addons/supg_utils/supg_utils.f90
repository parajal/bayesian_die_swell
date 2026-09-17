
! Copyright (C) 2007-2023 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Routines for the SUPG method.
!

module supg_utils_m

  use tfem_elem_m
  use math_defs_m, only: norm, sumsq
  use set_optional_m

  implicit none

! interface statement to support older interface

  interface supg_hoverU
    module procedure supg_hU
  end interface supg_hoverU

contains


! Routine for computing the h/U and, optionally, the h*U parameter in SUPG and
! the Jacobian of h/U wrt to the velocity vector components.

  subroutine supg_hU ( ndim, globalshape, hoverU, htype, hlocation, &
    Uscaling, Uglobal, uvec, uvecc, Finv, Finvc, x, esize, htimesU, &
    checkzero, dhoverU )

!   space dimension
    integer, intent(in) :: ndim

!   Global shape of the element:
!     quadrilateral
!     triangle
!     hexahedron
!     tetrahedron
!     prism
!     pyramid
    character (len=*), intent(in) :: globalshape

!   The h/U parameter. This must be an array of the size of the number of
!   integration points.
    real(dp), intent(out), dimension(:) :: hoverU

!   The h*U parameter. This must be an array of the size of the number of
!   integration points. Optional.
    real(dp), intent(out), dimension(:), optional :: htimesU

!   The Jacobian of the h/U parameter with respect ot the velocity vector
!   components. Array dhoverU(i,j) where i is the integration point number and
!   j the velocity vector component. Optional.
    real(dp), intent(out), dimension(:,:), optional :: dhoverU

!   these parameters determine the way h/U is computed:
!     htype, the type of h computation:
!        1: the standard velocity projection
!        2: mapping of the reference element
!        3: determined by the size of the element:
!            2D: h = (esize)**0.5
!            3D: h = (esize)**1/3
!            hence: h is a constant per element and hlocation is irrelevant.
!     hlocation, the location where h is computed:
!        1: in the center of the element only
!        2: in all integration points separately
!     Uscaling, how is the U scaling computed:
!        1: |u| in the center of the element only
!        2: 1/N sum_i |u|_i , i=1,..,N where N is the number of integration p.
!        3: |u| in the integration points separately
!        4: global scaling Uglobal
!   defaults: htype=2, hlocation=2, Uscaling=3
    integer, intent(in), optional :: htype, hlocation, Uscaling

!   the global scaling velocity for U; only needed for Uscaling=4
!   default: Uglobal=1.0
    real(dp), intent(in), optional :: Uglobal

!   the velocity vector in the integration points
!   uvec(ip,dim) is velocity vector in point in direction dim
    real(dp), intent(in), dimension(:,:), optional :: uvec

!   the velocity vector in the center of the element
!   uvecc(dim) is velocity vector in direction dim
    real(dp), intent(in), dimension(:), optional :: uvecc

!   Deformation gradient matrix between the reference element and the
!   actual element.
!   F(i,j,m) means that in point i the transformation is
!
!                         F(j,m) = d x  / d xi
!                                     j       m
!
!   The inverse of F is Finv
    real(dp), intent(in), dimension(:,:,:), optional :: Finv

!   The inverse of F (Finv) in the center of the element
    real(dp), intent(in), dimension(:,:), optional :: Finvc

!   coordinates of the nodes of the element
!   x(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:), optional :: x

!   The size (area in 2D, volume in 3D) of the element
!   Only needed for htype=3
    real(dp), intent(in), optional :: esize

!   Check for zero velocity and use htype=3 and Uscaling=4 instead.
!   If a zero velocity is detected in ANY of the points, the method is
!   switched to htype=3 and Uscaling=4 for ALL points.
!   NOTE, this is a brute force approach to avoid divide by zero, when
!   the velocity is exactly zero. It should be used only for that case
!   and normally be off (checkzero=.false.).
!   Default=.false.
    logical, intent(in), optional :: checkzero


    logical :: doFref, lcheckzero, zerovec
    integer :: ninti, ip, lhtype, lhlocation, lUscaling
    real(dp) :: h0, h, Finvu(ndim), wk, lUglobal
    real(dp), dimension(size(hoverU)) :: work
!   Fref is a matrix for (pre-) deforming the reference triangle/tetrahedron to
!   a regular one with the length of all the edges equal to 1.
!   FFinv = Fref*Finv if doFref=.true. otherwise FFinv = Finv
    real(dp), dimension(ndim,ndim) :: Fref, FFinv
    real(dp), parameter :: eps = 10 * tiny(1._dp)


    ninti = size ( hoverU )

!   set optional parameters

    lhtype     = set_optional ( variable=htype, default=2 )
    lhlocation = set_optional ( variable=hlocation, default=2 )
    lUscaling  = set_optional ( variable=Uscaling, default=3 )
    lUglobal   = set_optional ( variable=Uglobal, default=1._dp )
    lcheckzero = set_optional ( variable=checkzero, default=.false. )

!   check for zero velocity

    zerovec = .false.

    if ( lcheckzero ) then

      if ( present(uvec) ) then
        if ( any( sum(abs(uvec),dim=2) < eps ) ) zerovec = .true.
      end if

      if ( present(uvecc) ) then
        if ( sum(abs(uvecc)) < eps ) zerovec = .true.
      end if

      if ( zerovec ) then
        lhtype = 3
        lUscaling = 4
      end if

    end if

!   test presence of other optional parameters

    if ( lhlocation == 2 .or. lUscaling == 2 .or. lUscaling == 3 ) then
      if ( .not. present(uvec) ) then
        write(*,'(/a/a/)') 'Error in supg_hU:', &
          ' uvec not present in the heading'
        stop
      end if
    end if

    if ( lhlocation == 1 .or. lUscaling == 1 ) then
      if ( .not. present(uvecc) ) then
        write(*,'(/a/a/)') 'Error in supg_hU:', &
          ' uvecc not present in the heading'
        stop
      end if
    end if

    if ( lhtype == 1 ) then
      if ( .not. present(x) ) then
        write(*,'(/a/a/)') 'Error in supg_hU:', &
          ' x not present in the heading'
        stop
      end if
    else if ( lhtype == 2 .and. lhlocation == 1 ) then
      if ( .not. present(Finvc) ) then
        write(*,'(/a/a/)') 'Error in supg_hU:', &
          ' Finvc not present in the heading'
        stop
      end if
    else if ( lhtype == 2 .and. lhlocation == 2 ) then
      if ( .not. present(Finv) ) then
        write(*,'(/a/a/)') 'Error in supg_hU:', &
          ' Finv not present in the heading'
        stop
      end if
    else if ( lhtype == 3 ) then
      if ( .not. present(esize) ) then
        write(*,'(/a/a/)') 'Error in supg_hU:', &
          ' esize not present in the heading'
        stop
      end if
    end if

    if ( lUscaling == 4 ) then
      if ( .not. present(Uglobal) ) then
        write(*,'(/a/a/)') 'Error in supg_hU:', &
          ' Uglobal not present in the heading'
        stop
      else if ( Uglobal < eps ) then
        write(*,'(/a/a/)') 'Error in supg_hU:', &
          ' Uglobal negative or near zero '
        stop
      end if
    end if

    if ( present(dhoverU) ) then
      if ( zerovec ) then
        dhoverU = 0
      else if ( any( [ lhtype, lhlocation, lUscaling] /= [2,2,3] ) ) then
        write(*,'(/a/a/)') 'Error in supg_hU:', &
          ' dhoverU only available for [lhtype,lhlocation,lUscaling] = [2,2,3]'
        stop
      end if
    end if

!   choose various options

    if ( lhtype == 1 ) then

!     velocity projection

      if ( lhlocation == 1 ) then

!       h in center point

        select case ( lUscaling )

        case(1)

!         U: u in center point

          hoverU = projection_hU ( uvecc )
          if ( present(htimesU) ) htimesU = hoverU ! h*U
          hoverU = hoverU / sumsq ( uvecc ) ! h/U

        case(2)

!         U: u average of |u| in integration points

          h =  projection_hU ( uvecc ) / norm ( uvecc ) ! h
          work = norm ( uvec, dim=2 )
          wk = sum(work) / ninti
          hoverU = h / wk
          if ( present(htimesU) ) htimesU = h * wk  ! h*U

        case(3)

!         U: u in integration point

          h = projection_hU ( uvecc ) / norm ( uvecc ) ! h
          work = norm ( uvec, dim=2 )
          hoverU = h / work
          if ( present(htimesU) ) htimesU = h * work

        case(4)

!         U: global velocity scale

          h =  projection_hU ( uvecc ) / norm ( uvecc ) ! h
          hoverU = h / lUglobal
          if ( present(htimesU) ) htimesU = h * lUglobal

        case default

          write(*,'(/a/a,i0/)') 'Error in supg_hU:', &
            ' incorrect Uscaling parameter = ', lUscaling
          stop

        end select

      else if ( lhlocation == 2 ) then

!       h in integration points

        select case ( lUscaling )

        case(1)

!         U: u in center point

          do ip = 1, ninti
            hoverU(ip) = projection_hU ( uvec(ip,:) ) / norm ( uvec(ip,:) ) ! h
          end do
          wk = norm ( uvecc )
          if ( present(htimesU) ) htimesU = hoverU * wk  ! h*U
          hoverU = hoverU / wk  ! h/U

        case(2)

!         U: u average of |u| in integration points

          work = norm ( uvec, dim=2 )
          do ip = 1, ninti
            hoverU(ip) = projection_hU ( uvec(ip,:) ) / work(ip)  ! h
          end do
          wk = sum(work) / ninti
          if ( present(htimesU) ) htimesU = hoverU * wk  ! h*U
          hoverU = hoverU / wk  ! h/U

        case(3)

!         U: u in integration point

          do ip = 1, ninti
            hoverU(ip) = projection_hU ( uvec(ip,:) )
          end do
          if ( present(htimesU) ) htimesU = hoverU   ! h*U
          hoverU = hoverU / sumsq ( uvec, dim=2 )  ! h/U

        case(4)

!         U: global velocity scale

          do ip = 1, ninti
            hoverU(ip) = projection_hU ( uvec(ip,:) ) / norm ( uvec(ip,:) )
          end do
          if ( present(htimesU) ) htimesU = hoverU * lUglobal  ! h*U
          hoverU = hoverU / lUglobal

        case default

          write(*,'(/a/a,i0/)') 'Error in supg_hU:', &
            ' incorrect Uscaling parameter = ', lUscaling
          stop

        end select

      else

        write(*,'(/a/a,i0/)') 'Error in supg_hU:', &
          ' incorrect hlocation parameter = ', lhlocation
        stop

      end if

    else if ( lhtype == 2 ) then

!     element mapping

!     set h0 and Fref
      select case ( globalshape )
        case ( 'quadrilateral', 'hexahedron' )
          h0 = 2   ! reference size of the element
          doFref = .false.
        case ( 'triangle' )
          h0 = 3._dp**0.25_dp/2   ! reference size of the element
          Fref = transpose ( reshape ( [ 1._dp, 0.5_dp, &
                                         0._dp, sqrt(3._dp)/2 ], [2,2] ) )
          doFref = .true.
        case ( 'tetrahedron' )
          h0 = ( sqrt(2._dp)/12 ) ** ( 1._dp / 3 )  ! ref size of the element
          Fref = transpose ( reshape ( [ 1._dp, 0.5_dp, 0.5_dp, &
                                         0._dp, sqrt(3._dp)/2, sqrt(3._dp)/6, &
                                         0._dp, 0._dp, sqrt(6._dp)/3 ], &
                                       [3,3] ) )
          doFref = .true.
        case ( 'prism' )
          h0 = ( sqrt(3._dp)/2 ) ** ( 1._dp / 3 )  ! ref size of the element
          Fref = transpose ( reshape ( [ 1._dp, 0.5_dp, 0._dp, &
                                         0._dp, sqrt(3._dp)/2, 0._dp, &
                                         0._dp, 0._dp, 1._dp ], [3,3] ) )
          doFref = .true.
        case ( 'pyramid' )
          h0 = 4._dp ** ( 1._dp / 3 )  ! reference size of the element
          doFref = .false.
        case default
          write(*,'(/a/2a/)') 'Error in supg_hU:', &
            ' invalid element shape: ', globalshape
          stop
      end select

      if ( lhlocation == 1 ) then

!       h in center point

        select case ( lUscaling )

        case(1)

!         U: u in center point

          Finvu = matmul ( Finvc, uvecc )   !  F^-1 u
          if ( doFref ) Finvu = matmul ( Fref, Finvu )   !  Fref * F^-1 u
          hoverU = h0 / norm ( Finvu )   ! h/U
          if ( present(htimesU) ) htimesU = hoverU * sumsq(uvecc)  ! h*U

        case(2)

!         U: u average of |u| in integration points

          Finvu = matmul ( Finvc, uvecc )   !  F^-1 u
          if ( doFref ) Finvu = matmul ( Fref, Finvu )   !  Fref * F^-1 u
          h = h0 * norm ( uvecc ) / norm ( Finvu )   ! h
          work = norm ( uvec, dim=2 )
          wk = sum(work) / ninti
          hoverU = h / wk  ! h/U
          if ( present(htimesU) ) htimesU = h * wk  ! h*U

        case(3)

!         U: u in integration point

          Finvu = matmul ( Finvc, uvecc )   !  F^-1 u
          if ( doFref ) Finvu = matmul ( Fref, Finvu )   !  Fref * F^-1 u
          h = h0 * norm ( uvecc ) / norm ( Finvu )   ! h
          work = norm ( uvec, dim=2 )
          hoverU = h / work
          if ( present(htimesU) ) htimesU = h * work   ! h*U

        case(4)

!         U: global velocity scale

          Finvu = matmul ( Finvc, uvecc )   !  F^-1 u
          if ( doFref ) Finvu = matmul ( Fref, Finvu )   !  Fref * F^-1 u
          h = h0 * norm ( uvecc ) / norm ( Finvu )   ! h
          hoverU = h / lUglobal
          if ( present(htimesU) ) htimesU = h * lUglobal  ! h*U

        case default

          write(*,'(/a/a,i0/)') 'Error in supg_hU:', &
            ' incorrect Uscaling parameter = ', lUscaling
          stop

        end select

      else if ( lhlocation == 2 ) then

!       h in integration points

        select case ( lUscaling )

        case(1)

!         U: u in center point

          do ip = 1, ninti
            Finvu = matmul ( Finv(ip,:,:), uvec(ip,:) )   !  F^-1 u
            if ( doFref ) Finvu = matmul ( Fref, Finvu )  !  Fref * F^-1 u
            hoverU(ip) = h0 * norm ( uvec(ip,:) ) / norm ( Finvu )   ! h
          end do
          wk = norm ( uvecc )
          if ( present(htimesU) ) htimesU = hoverU * wk  ! h*U
          hoverU = hoverU / wk  ! h/U

        case(2)

!         U: u average of |u| in integration points

          work = norm ( uvec, dim=2 )
          do ip = 1, ninti
            Finvu = matmul ( Finv(ip,:,:), uvec(ip,:) )   !  F^-1 u
            if ( doFref ) Finvu = matmul ( Fref, Finvu )  !  Fref * F^-1 u
            hoverU(ip) = h0 * work(ip) / norm ( Finvu )   ! h
          end do
          wk = sum(work) / ninti
          if ( present(htimesU) ) htimesU = hoverU * wk  ! h*U
          hoverU = hoverU / wk  ! h/U

        case(3)

!         U: u in integration point

          do ip = 1, ninti
            if ( doFref ) then
              FFinv = matmul ( Fref, Finv(ip,:,:) )  !  Fref * F^-1
            else
              FFinv = Finv(ip,:,:)  ! F^-1
            end if
            Finvu = matmul ( FFinv, uvec(ip,:) )   !  F^-1 u
            hoverU(ip) = h0 / norm ( Finvu )   ! h/U
            if ( present(dhoverU) ) then
              dhoverU(ip,:) = &
                   - h0 / norm ( Finvu ) ** 3 * matmul ( Finvu, FFinv )
            end if
          end do
          if ( present(htimesU) ) htimesU = hoverU * sumsq(uvec,dim=2)  ! h*U

        case(4)

!         U: global velocity scale

          do ip = 1, ninti
            Finvu = matmul ( Finv(ip,:,:), uvec(ip,:) )   !  F^-1 u
            if ( doFref ) Finvu = matmul ( Fref, Finvu )   !  Fref * F^-1 u
            hoverU(ip) = h0 * norm ( uvec(ip,:) ) / norm ( Finvu ) ! h
          end do
          if ( present(htimesU) ) htimesU = hoverU * lUglobal  ! h*U
          hoverU = hoverU / lUglobal   ! h/U

        case default

          write(*,'(/a/a,i0/)') 'Error in supg_hU:', &
            ' incorrect Uscaling parameter = ', lUscaling
          stop

        end select

      else

        write(*,'(/a/a,i0/)') 'Error in supg_hU:', &
          ' incorrect hlocation parameter = ', lhlocation
        stop

      end if

    else if ( lhtype == 3 ) then

!     constant h given by size of the element

!     set h
      select case ( globalshape )
        case ( 'quadrilateral' )
          h = sqrt(esize)
        case ( 'triangle' )
          h = sqrt(esize)
        case ( 'hexahedron', 'tetrahedron', 'prism', 'pyramid' )
          h = esize ** ( 1._dp / 3 )
        case default
          write(*,'(/a/2a/)') 'Error in supg_hU:', &
            ' invalid element shape: ', globalshape
          stop
      end select

      select case ( lUscaling )

      case(1)

!       U: u in center point

        wk = norm ( uvecc )
        hoverU = h / wk   ! h/U
        if ( present(htimesU) ) htimesU = h * wk   ! h*U

      case(2)

!       U: u average of |u| in integration points

        work = norm ( uvec, dim=2 )
        wk = sum(work) / ninti
        hoverU = h / wk  ! h/U
        if ( present(htimesU) ) htimesU = h * wk   ! h*U

      case(3)

!       U: u in integration point

        work = norm ( uvec, dim=2 )
        hoverU = h / work  ! h/U
        if ( present(htimesU) ) htimesU = h * work  ! h*U

      case(4)

!       U: global velocity scale

        hoverU = h / lUglobal  ! h/U
        if ( present(htimesU) ) htimesU = h * lUglobal   ! h*U

      case default

        write(*,'(/a/a,i0/)') 'Error in supg_hU:', &
          ' incorrect Uscaling parameter = ', lUscaling
        stop

      end select

    else

      write(*,'(/a/a,i0/)') 'Error in supg_hU:', &
        ' incorrect htype parameter = ', lhtype
      stop

    end if

  contains


!   inline routine for h1, h2, h3 computation

    subroutine h1h2h3 ( h1, h2, h3 )

      real(dp), intent(out), dimension(:) :: h1, h2, h3

      integer :: nodalp

      nodalp = size(x,1)

      if ( globalshape == 'quadrilateral' .and. nodalp == 9 ) then

        h1 = ( x(3,:) + x(5,:) - x(1,:) - x(7,:) ) / 2
        h2 = ( x(7,:) + x(5,:) - x(1,:) - x(3,:) ) / 2

      else if ( globalshape == 'hexahedron' .and. nodalp == 27 ) then

        h1 = (   x(3,:) + x(9,:) + x(27,:) + x(21,:) &
               - x(1,:) - x(7,:) - x(25,:) - x(19,:) ) / 4
        h2 = (   x(1,:) + x(3,:) + x(21,:) + x(19,:) &
               - x(7,:) - x(9,:) - x(27,:) - x(25,:) ) / 4
        h3 = (   x(19,:) + x(21,:) + x(27,:) + x(25,:) &
               - x(1,:) - x(3,:) - x(9,:) - x(7,:) ) / 4

      else

        write(*,'(/a/3a,i0,a/)') 'Error in routine h1h2h3 (in supg_hU):', &
          ' invalid element shape: ', globalshape, &
          ' with ', nodalp, ' nodal points.'
        stop

      end if

    end subroutine h1h2h3


!   inline routine for projection h computation

    function projection_hU ( uvector )

      real(dp), intent(in), dimension(:) :: uvector
      real(dp) :: projection_hU

      real(dp) :: h1(ndim), h2(ndim), h3(ndim)

      call h1h2h3 ( h1, h2, h3 )

      if ( ndim == 2 ) then

        projection_hU =  abs ( dot_product ( uvector, h1 ) ) +   &
                         abs ( dot_product ( uvector, h2 ) )

      else if ( ndim == 3 ) then

        projection_hU = abs ( dot_product ( uvector, h1 ) ) +   &
                        abs ( dot_product ( uvector, h2 ) ) +   &
                        abs ( dot_product ( uvector, h3 ) )

      end if

    end function projection_hU

  end subroutine supg_hU


! function routine for parameter beta

  function supg_beta ( choice, Peh, betaconst )

!   choice for the beta parameter
    integer, intent(in) :: choice

!   mesh Peclet number, Peh > 0, size gives size of output array
    real(dp), dimension(:), intent(in) :: Peh

!   constant value for beta (choice=5)
    real(dp), intent(in), optional :: betaconst

    real(dp), dimension(size(Peh)) :: supg_beta


    if ( any ( choice == [2,3,4] ) ) then
      if ( any ( Peh < 0 ) ) then
        write(*,'(/a/)') &
          'Error supg_beta: mesh Peclet number has a negative value'
        stop
      end if
    else if ( choice == 5 .and. .not. present(betaconst) )  then
      write(*,'(/a/)') &
        'Error supg_beta: betaconst not present'
      stop
    end if

!   choose beta parameter

    select case ( choice )
    case(0)
      supg_beta = 0
    case(1)
      supg_beta = 1
    case(2)
      where ( Peh < 1 )
        supg_beta = 0
      else where
        supg_beta = 1 - 1/Peh
      end where
    case(3)
      where ( Peh > 0 )
        supg_beta = 1/tanh(Peh) - 1/Peh
      else where
        supg_beta = 0
      end where
    case(4)
      where ( Peh < 3 )
        supg_beta = Peh/3
      else where
        supg_beta = 1
      end where
    case(5)
      supg_beta = betaconst
    case default
      call errormsg_case_default ( 'supg_beta', 'choice', int_value=choice )
    end select

  end function supg_beta

end module supg_utils_m

