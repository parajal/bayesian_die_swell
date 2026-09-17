
! Copyright (C) 2004-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Defines shape functions and associated routines

module shapefunc_m

  use glob_defs_m
  use shapefunc_standard_m
  use shapefunc_standard_global_m
  use shapefunc_macro_m
  use shapefunc_spectral_m
  use shapefunc_highorder_m
  use mesh_m, only: element_t, mesh_t
  use set_optional_m

  implicit none


! interface to support older interfaces for shapefunction chooser

  interface set_shape_function
    module procedure set_shape_function, set_shape_function_simple
  end interface set_shape_function

! interface to support older interfaces

  interface isoparametric_deformation_surface
    module procedure isoparametric_deformation_curved
  end interface isoparametric_deformation_surface

  interface isoparametric_deformation_geometry
    module procedure isoparametric_deformation_curved
  end interface isoparametric_deformation_geometry

! type definition of a shape function

  type shapefunc_t

!   Global shape of the element:
!     line
!     quadrilateral
!     triangle
!     hexahedron
!     tetrahedron
!     prism
!     pyramid
    character (len=13) :: globalshape = ''

!   Interpolation:
!      1:P0   2:P1   3:P1+   4:Q1   5:Q1+   6:P2   7:P2+   8:Q2
!      9:P1isoP2    10:Q1isoQ2     11:serendipity (quadratic)
!     12:P2+ on tetrahedron with face bubbles only (14 nodes)
!     13:Qp spectral, polynomial order=p
!     14:Pp spectral, Legrendre polymials, polynomial order=p
!     15:P1Q1 on prism with 6 nodes
!     16:P2Q2 on prism with 18 nodes
!     17:Q1P1r on pyramid with 5 nodes
!     18:Q2P2r on pyramid with 14 nodes
!     19:Pp high order, polynomial order=p
!     20:Qp high order, polynomial order=p
!     21:Pp high order, Legrendre polymials, polynomial order=p
!     22:P1isoPp macro high order, polynomial order=p
!     23:Q1isoQp macro high order, polynomial order=p
    integer :: interpolation = 0

!   polynomial order for spectral and high-order elements
    integer :: p = 0

!   intrule: integration rule for spectral elements of interpolation=14 and
!            spec_eval='gauss'
!   Needed for:
!     Gauss-Legendre-Lobatto on line/quadrilateral/hexahedron:
!       intrule = number of Gauss points in one dimension
    integer :: intrule = 0

!   Numbering of the degrees of freedom in the nodes according to
!     'standard'  (=nodal numbering according to the Sepran standard)
!     'regular'   (=regular numbering)
!   Default: standard.
!   Currently this only makes a difference for the Q1 on a hexahedron:
!     - On a iso-parametric Q1-element with 8 nodes, you need to use 'standard'.
!     - On an element with 27 nodes with only Q1 interpolation, you need
!       to use 'regular'. For example on a Q2/Q1 velocity/pressure element
!       the pressure interpolation needs `regular' numbering.
    character (len=8) :: numbering = 'standard'

!   The shape functions and derivatives for spectral elements need
!   to be evaluated as follows:
!    'general' : in general points
!    'gauss'   : in the gauss points corresponding to the polynomial order
!   Default: gauss.
    character (len=8) :: spec_eval = 'gauss'

  end type shapefunc_t


contains


! set number of degrees of freedom in an element

  subroutine set_ndf ( shapefunc, name_of_routine, ndf, ndfb, ndfl, intpoll )

    type(shapefunc_t), intent(in) :: shapefunc

!   name of calling routine
    character (len=*), intent(in), optional :: name_of_routine

!   number of degrees of freedom
    integer, intent(out), optional :: ndf

!   number of degrees of freedom on a single edge of an element
    integer, intent(out), optional :: ndfb

!   number of Lagrange multiplier degrees of freedom for constraints
!   taken one order lower than the primary unknown
    integer, intent(out), optional :: ndfl

!   interpolation of the Lagrange multiplier
!   taken one order lower than the primary unknown
    integer, intent(out), optional :: intpoll


    character (len=:), allocatable :: lname_of_routine
    integer :: l_ndf, l_ndfb, l_ndfl, l_intpoll

!   set name of routine

    lname_of_routine = set_optional ( variable=name_of_routine, &
                                      default='set_ndf' )

!   find ndf

    if ( shapefunc%globalshape == 'line' ) then

      select case ( shapefunc%interpolation )
        case(1) ! constant
          l_ndf = 1
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
        case(2,3,4,15) ! linear
          l_ndf = 2
          l_ndfb = 1
          l_ndfl = 1
          l_intpoll = 1
        case(6,7,8,9,10,11,16) ! quadratic
          l_ndf = 3
          l_ndfb = 1
          l_ndfl = 2
          l_intpoll = 2
        case(13,14) ! spectral
          l_ndf = shapefunc%p + 1
          l_ndfb = 1
          l_ndfl = 0   ! Lagrange multiplier undefined for spectral elements
          l_intpoll = 0
        case(19,20,22,23) ! high-order
          l_ndf = shapefunc%p + 1
          l_ndfb = 1
          l_ndfl = 0   ! Lagrange multiplier undefined for high-order elements
          l_intpoll = 0
        case default ! not set
          l_ndf = 0
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
      end select

    else if ( shapefunc%globalshape == 'triangle' ) then

      select case ( shapefunc%interpolation )
        case(1) ! constant
          l_ndf = 1
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
        case(2,15) ! linear triangle
          l_ndf = 3
          l_ndfb = 2
          l_ndfl = 1
          l_intpoll = 1
        case(3) ! extended linear triangle
          l_ndf = 4
          l_ndfb = 2
          l_ndfl = 1
          l_intpoll = 1
        case(6,16) ! quadratic triangle
          l_ndf = 6
          l_ndfb = 3
          l_ndfl = 3
          l_intpoll = 2
        case(7) ! extended quadratic triangle
          l_ndf = 7
          l_ndfb = 3
          l_ndfl = 3
          l_intpoll = 2
        case(9) ! macro-element P1isoP2 (4 linear subtriangles)
          l_ndf = 6
          l_ndfb = 3
          l_ndfl = 3
          l_intpoll = 2
        case(19,21) ! high-order
          l_ndf = (shapefunc%p+1)*(shapefunc%p+2)/2
          l_ndfb = shapefunc%p+1
          l_ndfl = 0   ! Lagrange multiplier undefined for high-order elements
          l_intpoll = 0
        case default ! not set
          l_ndf = 0
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
      end select

    else if ( shapefunc%globalshape == 'quadrilateral' ) then

      select case ( shapefunc%interpolation )
        case(1) ! constant
          l_ndf = 1
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
        case(2) ! linear
          l_ndf = 3
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
        case(4,15) ! bilinear quad
          l_ndf = 4
          l_ndfb = 2
          l_ndfl = 1
          l_intpoll = 1
        case(5) ! extended bilinear quad
          l_ndf = 5
          l_ndfb = 2
          l_ndfl = 1
          l_intpoll = 1
        case(8,16) ! biquadratic quad
          l_ndf = 9
          l_ndfb = 3
          l_ndfl = 4
          l_intpoll = 4
        case(10) ! macroelement Q1isoQ2: four bilinear subquads
          l_ndf = 9
          l_ndfb = 3
          l_ndfl = 4
          l_intpoll = 4
        case(11) ! serendipity quadratic
          l_ndf = 8
          l_ndfb = 3
          l_ndfl = 4
          l_intpoll = 4
        case(13) ! spectral
          l_ndf = ( shapefunc%p + 1 ) ** 2
          l_ndfb = shapefunc%p + 1
          l_ndfl = 0   ! Lagrange multiplier undefined for spectral elements
          l_intpoll = 0
        case(14) ! spectral Pp with Legendre polynomials
          l_ndf = ( shapefunc%p + 1 ) * ( shapefunc%p + 2 ) / 2
          l_ndfb = shapefunc%p + 1
          l_ndfl = 0   ! Lagrange multiplier undefined for spectral elements
          l_intpoll = 0
        case(20) ! high-order
          l_ndf = ( shapefunc%p + 1 ) ** 2
          l_ndfb = shapefunc%p + 1
          l_ndfl = 0   ! Lagrange multiplier undefined for high-order elements
          l_intpoll = 0
        case(21) ! high-order Pp with Legendre polynomials
          l_ndf = ( shapefunc%p + 1 ) * ( shapefunc%p + 2 ) / 2
          l_ndfb = shapefunc%p + 1
          l_ndfl = 0   ! Lagrange multiplier undefined for high-order elements
          l_intpoll = 0
        case(23) ! macro-element high-order
          l_ndf = ( shapefunc%p + 1 ) ** 2
          l_ndfb = shapefunc%p + 1
          l_ndfl = 0   ! Lagrange multiplier undefined for high-order elements
          l_intpoll = 0
        case default ! not set
          l_ndf = 0
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
      end select

    else if ( shapefunc%globalshape == 'hexahedron' ) then

      select case ( shapefunc%interpolation )
        case(1) ! constant
          l_ndf = 1
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
        case(2) ! linear
          l_ndf = 4
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
        case(4) ! trilinear
          l_ndf = 8
          l_ndfb = 4
          l_ndfl = 1
          l_intpoll = 1
        case(8) ! triquadratic
          l_ndf = 27
          l_ndfb = 9
          l_ndfl = 8
          l_intpoll = 4
        case(10) ! macroelement Q1isoQ2: eight trilinear subhexa
          l_ndf = 27
          l_ndfb = 9
          l_ndfl = 8
          l_intpoll = 4
        case(11) ! serendipity quadratic
          l_ndf = 20
          l_ndfb = 8
          l_ndfl = 8
          l_intpoll = 4
        case(20) ! high-order
          l_ndf = ( shapefunc%p + 1 ) ** 3
          l_ndfb = shapefunc%p + 1
          l_ndfl = 0   ! Lagrange multiplier undefined for high-order elements
          l_intpoll = 0
        case(23) ! macro-element high-order
          l_ndf = ( shapefunc%p + 1 ) ** 3
          l_ndfb = shapefunc%p + 1
          l_ndfl = 0   ! Lagrange multiplier undefined for high-order elements
          l_intpoll = 0
        case default ! not set
          l_ndf = 0
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
      end select

    else if ( shapefunc%globalshape == 'tetrahedron' ) then

      select case ( shapefunc%interpolation )
        case(1) ! constant
          l_ndf = 1
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
        case(2) ! linear
          l_ndf = 4
          l_ndfb = 3
          l_ndfl = 1
          l_intpoll = 1
        case(6) ! quadratic
          l_ndf = 10
          l_ndfb = 6
          l_ndfl = 4
          l_intpoll = 2
        case(7) ! extended quadratic (with 4 face bubbles and 1 volume bubble)
          l_ndf = 15
          l_ndfb = 7
          l_ndfl = 4
          l_intpoll = 2
        case(9) ! macro-element P1isoP2 (8 linear subtets)
          l_ndf = 10
          l_ndfb = 6
          l_ndfl = 4
          l_intpoll = 2
        case(12) ! extended quadratic (with 4 face bubbles)
          l_ndf = 14
          l_ndfb = 7
          l_ndfl = 4
          l_intpoll = 2
        case default ! not set
          l_ndf = 0
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
      end select

    else if ( shapefunc%globalshape == 'prism' ) then

      select case ( shapefunc%interpolation )
        case(1) ! constant
          l_ndf = 1
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
        case(2) ! linear
          l_ndf = 4
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
        case(15) ! linear-bilinear
          l_ndf = 6
          l_ndfb = 0
          l_ndfl = 1
          l_intpoll = 1
        case(16) ! quadratic-biquadratic
          l_ndf = 18
          l_ndfb = 0
          l_ndfl = 6
          l_intpoll = 15
        case default ! not set
          l_ndf = 0
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
      end select

    else if ( shapefunc%globalshape == 'pyramid' ) then

      select case ( shapefunc%interpolation )
        case(1) ! constant
          l_ndf = 1
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
        case(2) ! linear
          l_ndf = 4
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
        case(17) ! bilinear-linear (rational)
          l_ndf = 5
          l_ndfb = 0
          l_ndfl = 1
          l_intpoll = 1
        case(18) ! biquadratic-quadratic (rational)
          l_ndf = 14
          l_ndfb = 0
          l_ndfl = 5
          l_intpoll = 17
        case default ! not set
          l_ndf = 0
          l_ndfb = 0
          l_ndfl = 0
          l_intpoll = 0
      end select

    end if

    if ( l_ndf == 0 ) then
      write(*,'(/3a/a,i0,2a/)') 'Error in ', lname_of_routine, ':', &
        ' incorrect interpolation = ', shapefunc%interpolation, &
        ' for globalshape = ', shapefunc%globalshape
      stop
    end if

    if ( present(ndf) ) ndf = l_ndf

    if ( present(ndfb) ) ndfb = l_ndfb

    if ( present(ndfl) ) then
      if ( l_ndfl == 0 ) then
        write(*,'(/3a/a,i0,2a/)') 'Error in ', lname_of_routine, ':', &
          ' for interpolation = ', shapefunc%interpolation, &
          ' on globalshape = ', shapefunc%globalshape, &
          ' Lagrange multiplier not defined '
        stop
      end if
      ndfl = l_ndfl
      intpoll = l_intpoll
    end if

  end subroutine set_ndf


! shape function chooser

  subroutine set_shape_function ( shapefunc, xr, phi, dphi )

    type(shapefunc_t), intent(in) :: shapefunc

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi


!   set shape function

    logical :: found
    integer :: ndf


    found = .true.

!   test number of degrees of freedom

    call set_ndf ( shapefunc, ndf=ndf )

    if ( size(phi,2) /= ndf ) then
      write(*,'(/a/a,i0,2a/a,i0,a/a,i0/)') &
        'Error in set_shape_function: ', &
        ' for interpolation = ', shapefunc%interpolation, &
        ' on globalshape = ', shapefunc%globalshape, &
        ' leads to ndf = ', ndf, ' degrees of freedom, wheras the ', &
        ' second dimension of phi is ', size(phi,2)
      stop
    end if

!   Select on global shape of the element

    select case ( shapefunc%globalshape )

    case ( 'line' )

      if ( present(dphi) ) then

        select case (shapefunc%interpolation)
          case(1)
            phi  = 1
            dphi = 0
          case(2:5,15,17)
            call shape_line_P1 ( xr(:,1), phi, dphi(:,:,1) )
          case(6:8,11,16,18)
            call shape_line_P2 ( xr(:,1), phi, dphi(:,:,1) )
          case(9:10)
            call shape_line_P1isoP2 ( xr(:,1), phi, dphi(:,:,1) )
          case(13)
            if ( shapefunc%spec_eval=='gauss' ) then
              call shape_line_at_GLL ( xr(:,1), phi, dphi(:,:,1), shapefunc%p )
            else
              call shape_line_GLL ( xr(:,1), shapefunc%p, phi, dphi(:,:,1) )
            end if
          case(19,20)
            call shape_line_Pp ( xr(:,1), shapefunc%p, phi, dphi(:,:,1) )
          case(22,23)
            call shape_line_P1isoPp ( xr(:,1), shapefunc%p, phi, dphi(:,:,1) )
          case default
            found = .false.
        end select

      else

        select case (shapefunc%interpolation)
          case(1)
            phi  = 1
          case(2:5,15,17)
            call shape_line_P1 ( xr(:,1), phi )
          case(6:8,11,16,18)
            call shape_line_P2 ( xr(:,1), phi )
          case(9:10)
            call shape_line_P1isoP2 ( xr(:,1), phi )
          case(13)
            if ( shapefunc%spec_eval=='gauss' ) then
              call shape_line_at_GLL ( xr(:,1), phi, p=shapefunc%p )
            else
              call shape_line_GLL ( xr(:,1), shapefunc%p, phi )
            end if
          case(19,20)
            call shape_line_Pp ( xr(:,1), shapefunc%p, phi )
          case(22,23)
            call shape_line_P1isoPp ( xr(:,1), shapefunc%p, phi )
          case default
            found = .false.
        end select

      end if

    case ( 'triangle' )

      select case (shapefunc%interpolation)
        case(1)
          phi  = 1
          if ( present(dphi) ) dphi = 0
        case(2,15,17)
          call shape_triangle_P1 ( xr, phi, dphi )
        case(3)
          call shape_triangle_P1plus ( xr, phi, dphi )
        case(6,16,18)
          call shape_triangle_P2 ( xr, phi, dphi )
        case(7)
          call shape_triangle_P2plus ( xr, phi, dphi )
        case(9)
          call shape_triangle_P1isoP2 ( xr, phi, dphi )
        case(14)
          call shape_quad_Pp ( xr, shapefunc%p, phi, dphi )
        case(19)
          call shape_triangle_Pp ( xr, shapefunc%p, phi, dphi )
        case(21)
          call shape_quad_Pp ( xr, shapefunc%p, phi, dphi )
        case default
          found = .false.
      end select

    case ( 'quadrilateral' )

      select case (shapefunc%interpolation)
        case(1)
          phi  = 1
          if ( present(dphi) ) dphi = 0
        case(2)
          call shape_quad_P1 ( xr, phi, dphi )
        case(4,15,17)
          call shape_quad_Q1 ( xr, phi, dphi )
        case(5)
          call shape_quad_Q1plus ( xr, phi, dphi )
        case(8,16,18)
          call shape_quad_Q2 ( xr, phi, dphi )
        case(10)
          call shape_quad_Q1isoQ2 ( xr, phi, dphi )
        case(11)
          call shape_quad_serendipity2 ( xr, phi, dphi )
        case(13)
          if ( shapefunc%spec_eval=='gauss' ) then
            call shape_quad_at_GLL ( xr, shapefunc%p, phi, dphi )
          else
            call shape_quad_GLL ( xr, shapefunc%p, phi, dphi )
          end if
        case(14)
          if ( shapefunc%spec_eval=='gauss' .and. shapefunc%intrule /= 0 ) then
            call shape_quad_Pp_at_GLL ( xr, shapefunc%p, phi, dphi, &
              shapefunc%intrule )
          else
            call shape_quad_Pp ( xr, shapefunc%p, phi, dphi )
          end if
        case(20)
          call shape_quad_Qp ( xr, shapefunc%p, phi, dphi )
        case(21)
          call shape_quad_Pp ( xr, shapefunc%p, phi, dphi )
        case(23)
          call shape_quad_Q1isoQp ( xr, shapefunc%p, phi, dphi )
        case default
          found = .false.
      end select

    case ( 'hexahedron' )

      select case (shapefunc%interpolation)
        case(1)
          phi  = 1
          if ( present(dphi) ) dphi = 0
        case(2)
          call shape_hexa_P1 ( xr, phi, dphi )
        case(4)
          if ( shapefunc%numbering == 'standard' ) then
            call shape_hexa_Q1 ( xr, phi, dphi )
          else
            call shape_hexa_Q1_reg ( xr, phi, dphi )
          end if
        case(8)
          call shape_hexa_Q2 ( xr, phi, dphi )
        case(10)
          call shape_hexa_Q1isoQ2 ( xr, phi, dphi )
        case(11)
          call shape_hexa_serendipity2 ( xr, phi, dphi )
        case(20)
          call shape_hexa_Qp ( xr, shapefunc%p, phi, dphi )
        case(23)
          call shape_hexa_Q1isoQp ( xr, shapefunc%p, phi, dphi )
        case default
          found = .false.
      end select

    case ( 'tetrahedron' )

      select case (shapefunc%interpolation)
        case(1)
          phi  = 1
          if ( present(dphi) ) dphi = 0
        case(2)
          call shape_tetra_P1 ( xr, phi, dphi )
        case(6)
          call shape_tetra_P2 ( xr, phi, dphi )
        case(7)
          call shape_tetra_P2plus15 ( xr, phi, dphi )
        case(9)
          call shape_tetra_P1isoP2 ( xr, phi, dphi )
        case(12)
          call shape_tetra_P2plus14 ( xr, phi, dphi )
        case default
          found = .false.
      end select

    case ( 'prism' )

      select case (shapefunc%interpolation)
        case(1)
          phi  = 1
          if ( present(dphi) ) dphi = 0
        case(2)
          call shape_prism_P1 ( xr, phi, dphi )
        case(15)
          call shape_prism_P1Q1 ( xr, phi, dphi )
        case(16)
          call shape_prism_P2Q2 ( xr, phi, dphi )
        case default
          found = .false.
      end select

    case ( 'pyramid' )

      select case (shapefunc%interpolation)
        case(1)
          phi  = 1
          if ( present(dphi) ) dphi = 0
        case(2)
          call shape_pyramid_P1 ( xr, phi, dphi )
        case(17)
          call shape_pyramid_Q1P1r ( xr, phi, dphi )
        case(18)
          call shape_pyramid_Q2P2r ( xr, phi, dphi )
        case default
          found = .false.
      end select

    case default

      call errormsg_case_default ( 'set_shape_function', &
        'shapefunc%globalshape', char_value=shapefunc%globalshape )

    end select

    if ( .not. found ) then
      write(*,'(/a/a,i0,2a/)') 'Error in set_shape_function:', &
        ' incorrect interpolation = ', shapefunc%interpolation, &
        ', for globalshape = ', shapefunc%globalshape
      stop
    end if

  end subroutine set_shape_function


! shape function chooser based on shape of the element and number of degrees
! Only standard shape functions can be used.

  subroutine set_shape_function_simple ( globalshape, xr, phi, dphi, numbering )

!   Global shape of the element:
!     line
!     quadrilateral
!     triangle
!     hexahedron
!     tetrahedron
!     prism
!     pyramid
    character (len=*), intent(in) :: globalshape

!   Reference coordinates where phi and dphi must be computed:
!   xr(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: xr

!   shape function phi(i,j), with i the point in space and j the unknown
!   Also ndf=size(phi,2), the number of degrees of freedom, is used to
!   determine the type of shape function.
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function: dphi(i,j,k), with i the point in space
!   j the unknown and k the direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphi

!   Numbering of the degrees of freedom in the nodes according to
!     'standard'  (=nodal numbering according to the Sepran standard)
!     'regular'   (=regular numbering)
!   Default: standard.
!   Currently this only makes a difference for the Q1 on a hexahedron:
!     - On a iso-parametric Q1-element with 8 nodes, you need to use 'standard'.
!     - On an element with 27 nodes with only Q1 interpolation, you need
!       to use 'regular'. For example on a Q2/Q1 velocity/pressure element
!       the pressure interpolation needs `regular' numbering.
    character (len=*), intent(in), optional :: numbering


!   set shape function

    logical :: found
    integer :: ndf, interpolation
    character (len=8) :: numberingl
    type(shapefunc_t) :: shapefunc


!   test numbering

    if ( present(numbering) ) then
      if ( all ( numbering /= [ 'standard', 'regular ' ] ) ) then
        write(*,'(/a/2a/)') &
          'Error in set_shape_function_simple: ', &
          ' heading parameter numbering must be either ''standard'', ', &
          ' or ''regular''.'
        stop
      end if
      numberingl = numbering
    else
      numberingl = 'standard'
    end if

    found = .true.
    ndf = size(phi,2)

!   Select on global shape of the element

    select case ( globalshape )

    case ( 'line' )

      select case (ndf)
        case(1)
          interpolation = 1
        case(2)
          interpolation = 2
        case(3)
          interpolation = 6
        case default
          found = .false.
      end select

    case ( 'triangle' )

      select case (ndf)
        case(1)
          interpolation = 1
        case(3)
          interpolation = 2
        case(4)
          interpolation = 3
        case(6)
          interpolation = 6
        case(7)
          interpolation = 7
        case default
          found = .false.
      end select

    case ( 'quadrilateral' )

      select case (ndf)
        case(1)
          interpolation = 1
        case(3)
          interpolation = 2
        case(4)
          interpolation = 4
        case(5)
          interpolation = 5
        case(8)
          interpolation = 11
        case(9)
          interpolation = 8
        case default
          found = .false.
      end select

    case ( 'hexahedron' )

      select case (ndf)
        case(1)
          interpolation = 1
        case(4)
          interpolation = 2
        case(8)
          interpolation = 4
        case(20)
          interpolation = 11
        case(27)
          interpolation = 8
        case default
          found = .false.
      end select

    case ( 'tetrahedron' )

      select case (ndf)
        case(1)
          interpolation = 1
        case(4)
          interpolation = 2
        case(10)
          interpolation = 6
        case(14)
          interpolation = 12
        case(15)
          interpolation = 7
        case default
          found = .false.
      end select

    case ( 'prism' )

      select case (ndf)
        case(1)
          interpolation = 1
        case(4)
          interpolation = 2
        case(6)
          interpolation = 15
        case(18)
          interpolation = 16
        case default
          found = .false.
      end select

    case ( 'pyramid' )

      select case (ndf)
        case(1)
          interpolation = 1
        case(4)
          interpolation = 2
        case(5)
          interpolation = 17
        case(14)
          interpolation = 18
        case default
          found = .false.
      end select

    case default

      call errormsg_case_default ( 'set_shape_function_simple', &
        'globalshape', char_value=globalshape )

    end select

    if ( .not. found ) then
      write(*,'(/a/a,i0,2a/)') 'Error in set_shape_function_simple:', &
        ' incorrect interpolation, ndf = size(phi,2) = ', ndf, &
        ', globalshape = ', globalshape
      stop
    end if

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = interpolation
    shapefunc%numbering = numberingl

    call set_shape_function ( shapefunc, xr, phi, dphi )

  end subroutine set_shape_function_simple


! shape function chooser (based on global coordinates)

  subroutine set_shape_function_global ( shapefunc, x, phi, dphidx )

    type(shapefunc_t), intent(in) :: shapefunc

!   Global coordinates where phi and dphidx must be computed:
!   x(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: x

!   shape function phi(i,j), with i the point in space and j the unknown
    real(dp), intent(out), dimension(:,:) :: phi

!   derivative of the shape function with respect to the global coordinates:
!   dphidx(i,j,k), with i the point in space j the unknown and k the
!   direction in space.
    real(dp), intent(out), dimension(:,:,:), optional :: dphidx


!   set shape function

    integer :: ndf

!   test number of degrees of freedom

    call set_ndf ( shapefunc, ndf=ndf )

    if ( size(phi,2) /= ndf ) then
      write(*,'(/a/a,i0,2a/a,i0,a/a,i0/)') &
        'Error in set_shape_function_global: ', &
        ' for interpolation = ', shapefunc%interpolation, &
        ' leads to ndf = ', ndf, ' degrees of freedom, wheras the ', &
        ' second dimension of phi is ', size(phi,2)
      stop
    end if

!   Select on interpolation

    select case (shapefunc%interpolation)
      case(2)
        call shape_global_P1 ( x, phi, dphidx )
      case default
        write(*,'(/a/a,i0,a/)') 'Error in set_shape_function_global:', &
          ' incorrect interpolation = ', shapefunc%interpolation
        stop
    end select

  end subroutine set_shape_function_global


! Spread reference points to a single side, based on a pattern in 1D on [-1,1]
! Triangle version

  subroutine spread_to_side_triangle ( side, xigb, xigs )

!   side number
    integer, intent(in) :: side

!   reference coordinates in 1D
    real(dp), intent(in), dimension(:) :: xigb

!   reference coordinates in 2D on a single side
    real(dp), intent(out), dimension(:,:) :: xigs

    integer :: n

    n = size(xigb)

    select case ( side )
      case(1)
        xigs(:,1) = ( xigb + 1 ) / 2
        xigs(:,2) = 0
      case(2)
        xigs(:,1) = ( xigb(n:1:-1) + 1 ) / 2
        xigs(:,2) = ( xigb + 1 ) / 2
      case(3)
        xigs(:,1) = 0
        xigs(:,2) = ( xigb(n:1:-1) +1 ) / 2
      case default
        write(*,'(/a/a,i0/)') 'Error spread_to_side_triangle: ', &
          ' wrong side number for triangle: side = ', side
        stop
    end select

  end subroutine spread_to_side_triangle


! Spread reference points to a single side, based on a pattern in 1D on [-1,1]
! Quad version

  subroutine spread_to_side_quad ( side, xigb, xigs )

!   side number
    integer, intent(in) :: side

!   reference coordinates in 1D
    real(dp), intent(in), dimension(:) :: xigb

!   reference coordinates in 2D on a single side
    real(dp), intent(out), dimension(:,:) :: xigs

    integer :: n

    n = size(xigb)

    select case ( side )
      case(1)
        xigs(:,1) = xigb
        xigs(:,2) = -1
      case(2)
        xigs(:,1) = 1
        xigs(:,2) = xigb
      case(3)
        xigs(:,1) = xigb(n:1:-1)
        xigs(:,2) = 1
      case(4)
        xigs(:,1) = -1
        xigs(:,2) = xigb(n:1:-1)
      case default
        write(*,'(/a/a,i0/)') 'Error spread_to_side_quad: ', &
          ' wrong side number for quad: side = ', side
        stop
    end select

  end subroutine spread_to_side_quad


! Spread reference points to the sides (2D), based on a pattern in 1D

  subroutine spread_to_sides_2D ( xigb, xigs )

!   reference coordinates in 1D
    real(dp), intent(in), dimension(:) :: xigb

!   reference coordinates in 2D on the sides
    real(dp), intent(out), dimension(:,:,:) :: xigs

    integer :: nsides, side

    nsides = size(xigs,3)

    select case ( nsides )
      case(3)
        do side = 1, nsides
          call spread_to_side_triangle ( side, xigb, xigs(:,:,side) )
        end do
      case(4)
        do side = 1, nsides
          call spread_to_side_quad ( side, xigb, xigs(:,:,side) )
        end do
      case default
        write(*,'(/a/a,i0/)') 'Error spread_to_sides_2D: ', &
          ' element not available: nsides = ', nsides
        stop
    end select

  end subroutine spread_to_sides_2D


! Reference coordinates of the nodal points

  subroutine refcoor_nodal_points ( mesh, elgrp, xrnod, blend )

    type(mesh_t), intent(in) :: mesh

!   element group elgrp
    integer, intent(in) :: elgrp

!   Reference coordinates of the nodal points:
!   xr(i,j) with i nodal point i and j the direction in space
    real(dp), intent(out), dimension(:,:) :: xrnod

!   use the element from the mesh blend
!   default=0 (main mesh)
    integer, intent(in), optional :: blend


    integer :: i, j, k, m, p, q
    real(dp) :: dx, dy, y
    type(element_t) :: element

!   choose element

    if ( present(blend) ) then
      if ( blend > 0 ) then
        element = mesh%element_blend(elgrp,blend)
      else
        element = mesh%element(elgrp)
      end if
    else
      element = mesh%element(elgrp)
    end if

!   fill xrnod for elshape

    select case ( element%elshape )

      case(1) ! two-node line element

        xrnod(1:2,1) = [ -1, 1 ]

      case(2) ! three-node line element

        xrnod(1:3,1) = [ -1, 0, 1 ]

      case(3) ! three-node triangular element

        xrnod(1:3,1:2) = &
          reshape ( [ 0, 1, 0, &
                      0, 0, 1 ], [3,2] )

      case(4,33) ! six-node triangular element

        xrnod(1:6,1:2) = &
          reshape ( [ 0._dp, 0.5_dp, 1._dp, 0.5_dp, 0._dp, 0._dp, &
                      0._dp, 0._dp,  0._dp, 0.5_dp, 1._dp, 0.5_dp ], [6,2] )

      case(5) ! four-node quadrilateral

        xrnod(1:4,1:2) = &
          reshape ( [ -1,  1, 1, -1, &
                      -1, -1, 1,  1  ], [4,2] )

      case(6,34) ! nine-node quadrilateral

        xrnod(1:9,1:2) = &
          reshape ( [ -1,  0,  1, 1, 1, 0, -1, -1, 0, &
                      -1, -1, -1, 0, 1, 1,  1,  0, 0 ], [9,2] )

      case(7) ! seven-node triangular element

        xrnod(1:7,1:2) = &
          reshape ( [ 0._dp, 0.5_dp, 1._dp, 0.5_dp, 0._dp,  0._dp, 1._dp/3, &
                      0._dp,  0._dp, 0._dp, 0.5_dp, 1._dp, 0.5_dp, 1._dp/3 ],&
                      [7,2] )

      case(9) ! five-node quadrilateral

        xrnod(1:5,1:2) = &
          reshape ( [ -1,  1, 1, -1, 0, &
                      -1, -1, 1,  1, 0 ], [5,2] )

      case(10) ! four-node triangular element

        xrnod(1:4,1:2) = &
          reshape ( [ 0._dp, 1._dp, 0._dp, 1._dp/3, &
                      0._dp, 0._dp, 1._dp, 1._dp/3 ], [4,2] )

      case(11) ! four-node tetrahedron element

        xrnod(1:4,1:3) = &
          reshape ( [ 0, 1, 0, 0, &
                      0, 0, 1, 0, &
                      0, 0, 0, 1 ], [4,3] )

      case(12,35) ! ten-node tetrahedron element

        xrnod(1:6,1:3) = &
          reshape ( [ 0._dp, 0.5_dp, 1._dp, 0.5_dp, 0._dp,  0._dp, &
                      0._dp, 0._dp,  0._dp, 0.5_dp, 1._dp, 0.5_dp, &
                      0._dp, 0._dp,  0._dp,  0._dp, 0._dp,  0._dp ], [6,3] )
        xrnod(7:10,1:3) = &
          reshape ( [ 0._dp,  0.5_dp, 0._dp,  0._dp, &
                      0._dp,  0._dp,  0.5_dp, 0._dp, &
                      0.5_dp, 0.5_dp, 0.5_dp, 1._dp  ], [4,3] )

      case(13) ! eight-node hexahedron element

        xrnod(1:8,1:3) = &
          reshape ( [ -1,  1,  1, -1, -1,  1, 1, -1, &
                      -1, -1,  1,  1, -1, -1, 1,  1, &
                      -1, -1, -1, -1,  1,  1, 1,  1 ], [8,3] )

      case(14,36) ! 27-node hexahedron element

        do i = 1, 3
          do j = 1, 3
            do k = 1, 3
              m = i + 3*(j-1) + 9*(k-1)
              xrnod(m,1:3) = [ i-2, j-2, k-2 ]
            end do
          end do
        end do

      case(15,16) ! 14- or 15-node tetrahedron element

        xrnod(1:6,1:3) = &
          reshape ( [ 0._dp, 0.5_dp, 1._dp, 0.5_dp, 0._dp,  0._dp, &
                      0._dp, 0._dp,  0._dp, 0.5_dp, 1._dp, 0.5_dp, &
                      0._dp, 0._dp,  0._dp,  0._dp, 0._dp,  0._dp ], [6,3] )
        xrnod(7:10,1:3) = &
          reshape ( [ 0._dp,  0.5_dp, 0._dp,  0._dp, &
                      0._dp,  0._dp,  0.5_dp, 0._dp, &
                      0.5_dp, 0.5_dp, 0.5_dp, 1._dp  ], [4,3] )
        xrnod(11:14,1:3) = &
          reshape ( [ 1._dp/3, 1._dp/3, 1._dp/3, 0._dp, &
                      1._dp/3, 0._dp,   1._dp/3, 1._dp/3, &
                      0._dp,   1._dp/3, 1._dp/3, 1._dp/3 ], [4,3] )

        if ( element%elshape == 16 ) then
          xrnod(15,1:3) = [ 1._dp/4, 1._dp/4, 1._dp/4 ]
        end if

      case(30) ! eight-node quadrilateral element

        xrnod(1:8,1:2) = &
          reshape ( [ -1,  0,  1, 1, 1, 0, -1, -1, &
                      -1, -1, -1, 0, 1, 1,  1,  0 ], [8,2] )

      case(31) ! twenty-node hexahedron element

        xrnod(1:8,1:3) = &
          reshape ( [ -1,  0,  1, -1,  1, -1,  0,  1, &
                      -1, -1, -1,  0,  0,  1,  1,  1, &
                      -1, -1, -1, -1, -1, -1, -1, -1 ], [8,3] )
        xrnod(9:12,1:3) = &
          reshape ( [ -1,  1, -1,  1, &
                      -1, -1,  1,  1, &
                       0,  0,  0,  0  ], [4,3] )
        xrnod(13:20,1:3) = &
          reshape ( [ -1,  0,  1, -1,  1, -1,  0,  1, &
                      -1, -1, -1,  0,  0,  1,  1,  1, &
                       1,  1,  1,  1,  1,  1,  1,  1 ], [8,3] )

      case(41) ! six-node prism element

        xrnod(1:3,1:2) = &
          reshape ( [ 0, 1, 0, &
                      0, 0, 1 ], [3,2] )
        xrnod(1:3,3) = -1
        xrnod(4:6,1:2) = xrnod(1:3,1:2)
        xrnod(4:6,3) = 1

      case(42) ! fifteen-node prism element

        xrnod(1:6,1:2) = &
          reshape ( [ 0._dp, 0.5_dp, 1._dp, 0.5_dp, 0._dp, 0._dp, &
                      0._dp, 0._dp,  0._dp, 0.5_dp, 1._dp, 0.5_dp ], [6,2] )
        xrnod(1:6,3) = -1
        xrnod(7:9,1:2) = xrnod(1:6:2,1:2)
        xrnod(7:9,3) = 0
        xrnod(10:15,1:2) = xrnod(1:6,1:2)
        xrnod(10:15,3) = 1

      case(43) ! eighteen-node prism element

        xrnod(1:6,1:2) = &
          reshape ( [ 0._dp, 0.5_dp, 1._dp, 0.5_dp, 0._dp, 0._dp, &
                      0._dp, 0._dp,  0._dp, 0.5_dp, 1._dp, 0.5_dp ], [6,2] )
        xrnod(1:6,3) = -1
        xrnod(7:12,1:2) = xrnod(1:6,1:2)
        xrnod(7:12,3) = 0
        xrnod(13:18,1:2) = xrnod(1:6,1:2)
        xrnod(13:18,3) = 1

      case(51) ! five-node pyramid element

        xrnod(1:4,1:2) = &
          reshape ( [ -1,  1, 1, -1, &
                      -1, -1, 1,  1  ], [4,2] )
        xrnod(1:4,3) = 0
        xrnod(5,:) = [ 0, 0, 1 ]

      case(52) ! thirteen-node pyramid element

        xrnod(1:8,1:2) = &
          reshape ( [ -1,  0,  1, 1, 1, 0, -1, -1, &
                      -1, -1, -1, 0, 1, 1,  1,  0 ], [8,2] )
        xrnod(1:8,3) = 0
        xrnod(9:12,1:2) = xrnod(1:8:2,1:2)
        xrnod(9:12,3) = 0.5_dp
        xrnod(13,:) = [ 0, 0, 1 ]

      case(53) ! fourteen-node pyramid element

        xrnod(1:9,1:2) = &
          reshape ( [ -1,  0,  1, 1, 1, 0, -1, -1, 0, &
                      -1, -1, -1, 0, 1, 1,  1,  0, 0 ], [9,2] )
        xrnod(1:9,3) = 0
        xrnod(10:13,1:2) = xrnod(1:8:2,1:2) / 2
        xrnod(10:13,3) = 0.5_dp
        xrnod(14,:) = [ 0, 0, 1 ]

      case(101,105) ! line high-order element (including macro variant)

        if ( element%p(1,2) == 1 ) then
          call Gauss_Legendre_Lobatto_line ( element%p(1,1)+1, &
            xrnod(:,1) )
        else
          p = element%p(1,1)
          dx = 2._dp / p
          xrnod(1,1) = -1
          xrnod(2:p,1) = [ ( -1.0_dp + dx * i, i = 1, p-1 ) ]
          xrnod(p+1,1) = 1
        end if

      case(102,106) ! quad high-order element (including macro variant)

        if ( all( element%p(1:2,2) == 1 ) ) then

!         spectral layout

          if ( element%p(1,1) /= &
               element%p(2,1) ) then
            write(*,'(/a/a/)') 'Error refcoor_nodal_points:', &
             ' Unequal order for spectral quad not available yet.'
            stop
          end if

          call Gauss_Legendre_Lobatto_quad ( element%p(1,1)+1, &
            xrnod )

        else if ( all( element%p(1:2,2) == 0 ) ) then

!         equidistant layout

          p = element%p(1,1)
          q = element%p(2,1)
          dx = 2._dp / p
          dy = 2._dp / q
          do j = 1, q+1
            m = (j-1)*(p+1)
            y = -1.0_dp + dy * ( j - 1 )
            do i = 1, p+1
              xrnod(i+m,1:2) = [ -1.0_dp + dx * ( i - 1 ), y ]
            end do
          end do

        else

          write(*,'(/a/a/)') 'Error refcoor_nodal_points:', &
           ' Unequal layout for high-order quad not available yet.'
          stop

        end if

      case(104) ! triangle high-order element

        if ( element%p(1,1) /= &
             element%p(2,1) ) then
          write(*,'(/a/a/)') 'Error refcoor_nodal_points:', &
           ' Unequal order for high-order triangle not available.'
          stop
        end if

        if ( all( element%p(1:2,2) == 0 ) ) then
          p = element%p(1,1)
          dx = 1._dp / p
          do j = 1, p+1
            m = (j-1)*(2*p+4-j)/2
            y = dx * ( j - 1 )
            do i = 1, p+2-j
              xrnod(i+m,1:2) = [ dx * ( i - 1 ), y ]
            end do
          end do
        else
          write(*,'(/a/a/)') 'Error refcoor_nodal_points:', &
           ' Unequal layout not available for high-order triangle.'
          stop
        end if

      case default

        write(*,'(/a/a,i0,a/)') 'Error refcoor_nodal_points:', &
         ' shape ', element%elshape, ' unknown.'
        stop

    end select

  end subroutine refcoor_nodal_points


! isoparametric coordinates

  subroutine isoparametric_coordinates ( x, phi, xphi )

!   coordinates of the nodes = position of the unknowns for phi and dphi
!   x(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: x

!   shape function
!   phi(i,j), with i the point in space j the unknown
    real(dp), intent(in), dimension(:,:) :: phi

!   real coordinates of the points given by phi
!   xphi(i,j) with i the point in space and j the direction in space
    real(dp), intent(out), dimension(:,:) :: xphi

    xphi = matmul ( phi, x )

  end subroutine isoparametric_coordinates


! isoparametric deformation

  subroutine isoparametric_deformation ( x, dphi, F, Finv, detF )

    use limits_m, only: CHECK_POSITIVE_DETF, SET_ABSOLUTE_DETF

!   coordinates of the nodes = position of the unknowns for phi and dphi
!   x(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: x

!   derivatives of the shape function with respect to the reference coordinates
!   dphi(i,j,k), with i the point in space j the unknown and k the direction in
!   space.
    real(dp), intent(in), dimension(:,:,:) :: dphi

!   Deformation gradient matrix between the reference element and the
!   actual element.
!   F(i,j,m) means that in point i the transformation is
!
!                         F(j,m) = d x  / d xi
!                                     j       m
!
!   The inverse of F is Finv
    real(dp), intent(out), dimension(:,:,:) :: F, Finv

!   The determinant of F (Jacobian)
!   detF(i) gives the determinant in point i
!   NOTE: the default is to output the absolute value, i.e. |detF|.
!         Set SET_ABSOLUTE_DETF=.false. in the limits module if you want to
!         include the sign, i.e. output of detF.
    real(dp), intent(out), dimension(:) :: detF


    integer :: npts, ndim, ip, i, j


    npts = size(dphi,1)
    ndim = size(x,2)

    do ip = 1, npts

!     compute F for point ip

      do i = 1, ndim
        do j = 1, ndim
          F(ip,i,j) = dot_product ( x(:,i), dphi(ip,:,j) )
        end do
      end do

      if ( ndim == 1 ) then

!       1D

        detF(ip) = F(ip,1,1)
        Finv(ip,1,1) = 1._dp / F(ip,1,1)

      else if ( ndim == 2 ) then

!       2D

        detF(ip) = F(ip,1,1) * F(ip,2,2) - F(ip,1,2) * F(ip,2,1)
        Finv(ip,1,1) = F(ip,2,2) / detF(ip)
        Finv(ip,1,2) = -F(ip,1,2) / detF(ip)
        Finv(ip,2,1) = -F(ip,2,1) / detF(ip)
        Finv(ip,2,2) = F(ip,1,1) / detF(ip)

      else if ( ndim == 3 ) then

!       3D

        detF(ip) =  F(ip,1,1)*F(ip,2,2)*F(ip,3,3) &
                   -F(ip,1,1)*F(ip,2,3)*F(ip,3,2) &
                   -F(ip,2,1)*F(ip,1,2)*F(ip,3,3) &
                   +F(ip,2,1)*F(ip,1,3)*F(ip,3,2) &
                   +F(ip,3,1)*F(ip,1,2)*F(ip,2,3) &
                   -F(ip,3,1)*F(ip,1,3)*F(ip,2,2)
        Finv(ip,1,1) = F(ip,2,2)*F(ip,3,3)-F(ip,2,3)*F(ip,3,2)
        Finv(ip,1,2) = -F(ip,1,2)*F(ip,3,3)+F(ip,1,3)*F(ip,3,2)
        Finv(ip,1,3) = F(ip,1,2)*F(ip,2,3)-F(ip,1,3)*F(ip,2,2)
        Finv(ip,2,1) = -F(ip,2,1)*F(ip,3,3)+F(ip,2,3)*F(ip,3,1)
        Finv(ip,2,2) = F(ip,1,1)*F(ip,3,3)-F(ip,1,3)*F(ip,3,1)
        Finv(ip,2,3) = -F(ip,1,1)*F(ip,2,3)+F(ip,1,3)*F(ip,2,1)
        Finv(ip,3,1) = F(ip,2,1)*F(ip,3,2)-F(ip,2,2)*F(ip,3,1)
        Finv(ip,3,2) = -F(ip,1,1)*F(ip,3,2)+F(ip,1,2)*F(ip,3,1)
        Finv(ip,3,3) = F(ip,1,1)*F(ip,2,2)-F(ip,1,2)*F(ip,2,1)

        Finv(ip,:,:) = Finv(ip,:,:) / detF(ip)

      end if

    end do

    if ( SET_ABSOLUTE_DETF ) then

!     take the absolute value always

      detF = abs(detF)

    else if ( CHECK_POSITIVE_DETF ) then

!     check if detF is positive ?

      if ( any ( detF <= 0._dp ) ) then

        write(*,'(/a/a/3es16.8/)') 'Error isoparametric_deformation: ', &
          ' detF < 0 in the element with coordinates of nodal point 1 = ', &
          x(1,:)
        stop

      end if

    end if

  end subroutine isoparametric_deformation


! Derivatives of shape function in real coordinates

  subroutine shape_derivative ( dphi, Finv, dphidx )

!   derivatives of the shape function with respect to the reference coordinates
!   dphi(i,j,k), with i the point in space j the unknown and k the direction in
!   space.
    real(dp), intent(in), dimension(:,:,:) :: dphi

!   Deformation gradient matrix between the reference element and the
!   actual element.
!   F(i,j,m) means that in point i the transformation is
!
!                         F(j,m) = d x  / d xi
!                                     j       m
!
!   The inverse of F is Finv
    real(dp), intent(in), dimension(:,:,:) :: Finv

!   derivatives of the shape function with respect to the real coordinates
!   dphidx(i,j,k), with i the point in space j the unknown and k the direction
!   in space.
    real(dp), intent(out), dimension(:,:,:) :: dphidx


!   The derivative of a function f(x_i), where x_i=x_k(eta_j) is a mapping of
!   the coordinates eta_j (reference coordinates) to the real coordinates, can
!   be found by
!
!      df/deta_i   = df/dx_j dx_j/dx_i = dx_j/dx_i df/xj
!
!   or in matrix notation
!
!                 T
!      df/deta = F  df/dx
!
!   or by taking the inverse
!
!               -T
!      df/dx = F   df/deta


    integer :: ip


    do ip = 1, size(Finv,1)

      dphidx(ip,:,:) = matmul ( dphi(ip,:,:), Finv(ip,:,:) )

    end do

  end subroutine shape_derivative


! Isoparametric deformation for curved line elements.
! NOTE: this is a routine limited to curved line elements (reference
! dimension = 1). Use the more general isoparametric_deformation_curved
! for new code, which works for both curved line and surface elements.

  subroutine isoparametric_deformation_curve ( x, dphi, dxdxi, curvel, &
    normal, g11_down, g11_up, g1_up )

!   coordinates of the nodes = position of the unknowns for phi and dphi
!   x(i,j) with i the point in space and j the direction in space
    real(dp), intent(in), dimension(:,:) :: x

!   derivative of the shape function with respect to the reference coordinate
!   dphi(i,j), with i the point in space j the unknown
    real(dp), intent(in), dimension(:,:) :: dphi

!   Derivative of the coordinates x with respect to the reference coordinate
!   (tangential vector).
!   dxdxi(i,j) means that in point i the derivative is
!
!                       dxdxi(j) = d x  / d xi
!                                     j
!
!   NOTE: this could alternatively be called g1_down forming the (tangential)
!   basis vector g  of a curvilinear coordinate system x(xi)
!                -1                                    -
    real(dp), intent(out), dimension(:,:) :: dxdxi

!   The length of the vector dxdxi
!   curvel(i) gives the |dxdxi| in point i
    real(dp), intent(out), dimension(:) :: curvel

!   Unit normal vector to the curve.
!   normal(i,j) means in point i the j-th component of n (n_j)
    real(dp), intent(out), dimension(:,:), optional :: normal

!   The metric "tensor" (scalar) g   of the curvilinear coordinate system x(xi)
!                                 11                                      -
!   in the tangential direction.
!
!   NOTE that:
!               g  = g .g
!                11  -1 -1
!
!   g11_down(i) means that in point i
!
!                       g11_down = g
!                                   11

    real(dp), intent(out), dimension(:), optional :: g11_down

!                                      11
!   The dual metric "tensor" (scalar) g   of the curvilinear coordinate
!   system x(xi) in the tangential direction.
!          -
!   NOTE that:   11                    11    1  1
!               g  = 1 / g   and also g   = g .g
!                         11                -  -
!
!   g11_up(i) means that in point i

!                                 11
!                       g11_up = g

    real(dp), intent(out), dimension(:), optional :: g11_up

!                          1
!   The dual basis vector g  of the curvilinear coordinate system x(xi) in the
!                        -                                       -
!   tangential direction.
!
!   NOTE that:   1    11
!               g  = g  g
!               -       -1
!
!   g1_up(i,j) means that in point i
!                                     1
!                     g1_up(j) = e . g      j=1,2,3
!                                -j  -
!
!   where e , j=1,2,3 are the Cartesian orthonormal basis vectors.
!         -j
!                                       1 d
!   NOTE: the curve gradient nabla_s = g  --
!                                      -  dxi

    real(dp), intent(out), dimension(:,:), optional :: g1_up


    integer :: npts, ndim, ip, i
    real(dp) :: curvel2(size(curvel))


    npts = size(dphi,1)
    ndim = size(x,2)

    do ip = 1, npts

!     compute dxdxi for point ip

      do i = 1, ndim
        dxdxi(ip,i) = dot_product ( x(:,i), dphi(ip,:) )
      end do

      curvel2(ip) = dot_product( dxdxi(ip,:), dxdxi(ip,:) )
      curvel(ip) = sqrt ( curvel2(ip) )

    end do

    if ( present(normal) ) then

!     compute normal vectors

      if ( ndim == 3 ) then
        write(*,'(/a/a/)') 'Error isoparametric_deformation_curve: ', &
          ' normal vectors for 3D curves not available.'
        stop
      end if

      normal(:,1) =   dxdxi(:,2) / curvel
      normal(:,2) = - dxdxi(:,1) / curvel

    end if

    if ( present(g11_down) ) then

!     metric tensor

      g11_down = curvel2

    end if

    if ( present(g11_up) ) then

!     dual metric tensor

      g11_up = 1 / curvel2

    end if

    if ( present(g1_up) ) then

!     dual basis vector

      do i = 1, ndim
        g1_up(:,i) = dxdxi(:,i) / curvel2
      end do

    end if

  end subroutine isoparametric_deformation_curve


! isoparametric deformation for a curved line or surface element.

  subroutine isoparametric_deformation_curved ( x, dphi, dxdxis, surfl, &
    normal, gij_down, gij_up, gi_up )

!   coordinates of the nodes = position of the unknowns for phi and dphi
!   x(i,j) with i the point in space and j the direction in space (j=1,..,ndim).
    real(dp), intent(in), dimension(:,:) :: x

!   derivatives of the shape function with respect to the reference coordinates
!   dphi(i,j,k), with i the point in space j the unknown and k the direction in
!   surface space (k=1,ndimr), where ndimr is reference space dimension.
    real(dp), intent(in), dimension(:,:,:) :: dphi


!   Derivative of the coordinates x with respect to the reference coordinates
!   (tangential vectors).
!   dxdxis(i,j,k) means that in point i the derivative is
!                                        k
!               dxdxis(j,k) = d x  / d xi     j=1,..,ndim; k=1,..,ndimr
!                                j
!
!   NOTE: this could alternatively be called gi_down forming the (tangential)
!                                                          k
!   basis vectors g  of a curvilinear coordinate system x(xi ), k=1,..,ndimr
!                 -i                                    -
    real(dp), intent(out), dimension(:,:,:) :: dxdxis

!   ndimr = 1 (curve):
!     The length of the vector dxdxis
!     surfl(i) gives the |dxdxi| in point i
!   ndimr = 2 (surface):
!     The length of the vector dxdxis_1 x dxdxis_2
!     surfl(i) gives the |dxdxi_1 x dxdxi_2| in point i
    real(dp), intent(out), dimension(:) :: surfl

!   Unit normal vector to the surface
!   normal(i,j) means in point i the j-th component of n (n_j)
    real(dp), intent(out), dimension(:,:), optional :: normal
!                                                                  k
!   The metric tensor g   of the curvilinear coordinate system x(xi ),
!                      ij                                      -
!   k=1,..,ndimr in the tangential direction.
!   NOTE that:
!               g  = g .g
!                ij  -i -j
!
!   gij_down(i,j,k) means that in point i
!
!               gij_down(j,k) = g      j=1,..,ndimr; k=1,..,ndimr
!                                jk

    real(dp), intent(out), dimension(:,:,:), optional :: gij_down

!                           ij                                          k
!   The dual metric tensor g   of the curvilinear coordinate system x(xi ),
!                                                                   -
!   k=1,..,ndimr in the tangential direction.
!
!   NOTE that:   ij       -1           ij    i  j
!               g  = [g  ]   and also g   = g .g
!                      ij                   -  -
!
!   gij_up(i,j,k) means that in point i

!                              jk
!               gij_up(j,k) = g      j=1,..,ndimr; k=1,..,ndimr

    real(dp), intent(out), dimension(:,:,:), optional :: gij_up

!                           i                                          k
!   The dual basis vectors g  of the curvilinear coordinate system x(xi ),
!                         -                                       -
!   k=1,..,ndimr in the tangential direction.
!
!   NOTE that:   i    ij
!               g  = g  g     (summation over j)
!               -       -j
!
!   gi_up(i,j,k) means that in point i
!                                         k
!                       gi_up(j,k) = e . g      j=1,..,ndim; k=1,..,ndimr
!                                    -j  -
!
!   where e , j=1,..,ndim are the Cartesian orthonormal basis vectors.
!         -j
!                                         i d
!   NOTE: the surface gradient nabla_s = g  -- .
!                                        -  dxi

    real(dp), intent(out), dimension(:,:,:), optional :: gi_up


    integer :: npts, ndim, ndimr
    real(dp) :: tl1(3), tl2(3), nrml(3), g11, g12, g22, g11d, g12d, g22d, detg
    real(dp) :: surfl2(size(surfl))


    npts = size(dphi,1)
    ndimr = size(dphi,3)
    ndim = size(x,2)

!   choose line_element or surface

    select case ( ndimr )
    case(1)
      call line_element
    case(2)
      call surface_element
    case default
      call errormsg_case_default ( 'isoparametric_deformation_curved', &
        'ndimr', int_value=ndimr )
    end select


  contains


!   compute geometric properties of a curved line element

    subroutine line_element

      integer :: ip, i

      do ip = 1, npts

!       compute dxdxi for point ip

        do i = 1, ndim
          dxdxis(ip,i,1) = dot_product ( x(:,i), dphi(ip,:,1) )
        end do

        surfl2(ip) = dot_product( dxdxis(ip,:,1), dxdxis(ip,:,1) )
        surfl(ip) = sqrt ( surfl2(ip) )

      end do

      if ( present(normal) ) then

!       compute normal vectors

        if ( ndim == 3 ) then
          write(*,'(/a/a/)') 'Error isoparametric_deformation_curved: ', &
            ' normal vectors for 3D curves not available.'
          stop
        end if

        normal(:,1) =   dxdxis(:,2,1) / surfl
        normal(:,2) = - dxdxis(:,1,1) / surfl

      end if

      if ( present(gij_down) ) then

!       metric tensor

        gij_down(:,1,1) = surfl2

      end if

      if ( present(gij_up) ) then

!       dual metric tensor

        gij_up(:,1,1) = 1 / surfl2

      end if

      if ( present(gi_up) ) then

!       dual basis vector

        do i = 1, ndim
          gi_up(:,i,1) = dxdxis(:,i,1) / surfl2
        end do

      end if

    end subroutine line_element


!   compute geometric properties of a curved surface element

    subroutine surface_element

      integer :: ip, i

      do ip = 1, npts

!       compute dxdxi for point ip

        do i = 1, ndim
          dxdxis(ip,i,1) = dot_product ( x(:,i), dphi(ip,:,1) )
          dxdxis(ip,i,2) = dot_product ( x(:,i), dphi(ip,:,2) )
        end do

        tl1(:ndim) = dxdxis(ip,:,1)
        tl2(:ndim) = dxdxis(ip,:,2)
        if ( ndim == 2 ) then
          tl1(3) = 0
          tl2(3) = 0
        end if

        nrml(1) = tl1(2) * tl2(3) - tl1(3) * tl2(2)
        nrml(2) = tl1(3) * tl2(1) - tl1(1) * tl2(3)
        nrml(3) = tl1(1) * tl2(2) - tl1(2) * tl2(1)

        surfl(ip) = sqrt ( dot_product ( nrml, nrml ) )

        if ( present(normal) ) then

!         compute normal vectors

          normal(ip,:) = nrml / surfl(ip)

        end if

      end do


      if ( present(gi_up) .or. present(gij_down) .or. present(gij_up) ) then

!       surface geometry

        do ip = 1, npts

          g11 = dot_product ( dxdxis(ip,:,1), dxdxis(ip,:,1) )
          g12 = dot_product ( dxdxis(ip,:,1), dxdxis(ip,:,2) )
          g22 = dot_product ( dxdxis(ip,:,2), dxdxis(ip,:,2) )

          detg = g11 * g22 - g12 ** 2

          g11d = g11 / detg
          g12d = g12 / detg
          g22d = g22 / detg

          if ( present(gij_down) ) then

            gij_down(ip,1,1) = g11
            gij_down(ip,1,2) = g12
            gij_down(ip,2,1) = g12
            gij_down(ip,2,2) = g22

          end if

          if ( present(gij_up) ) then

            gij_up(ip,1,1) = g22d
            gij_up(ip,1,2) = - g12d
            gij_up(ip,2,1) = - g12d
            gij_up(ip,2,2) = g11d

          end if

          if ( present(gi_up) ) then

            gi_up(ip,:,1) =   g22d * dxdxis(ip,:,1) - g12d * dxdxis(ip,:,2)
            gi_up(ip,:,2) = - g12d * dxdxis(ip,:,1) + g11d * dxdxis(ip,:,2)

          end if

        end do

      end if

    end subroutine surface_element

  end subroutine isoparametric_deformation_curved

end module shapefunc_m
