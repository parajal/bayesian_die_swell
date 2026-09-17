
! Copyright (C) 2007-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Gauss integration abscissas and weights

module gauss_m

  use glob_defs_m
  use gauss_defs_m
  use gauss_standard_m
  use gauss_macro_m
  use gauss_spectral_m
  use gauss_mapped_m
  use gauss_standard_numeric_m

  implicit none


! interface to support older interfaces for setting Gauss integration

  interface set_Gauss_integration
    module procedure set_Gauss_integration, set_Gauss_integration_simple
  end interface set_Gauss_integration

contains


! set intrule for specified order of integration

  function set_intrule ( globalshape, inttype, order )

    use limits_m, only: GAUSS_RULE

    character(len=*), intent(in) :: globalshape
    integer, intent(in) :: inttype, order
    integer :: set_intrule
    integer, parameter :: ntri(0:6) = [ 1, 1, 3, 6, 6, 7, 12 ]
    integer, parameter :: ntet(0:6) = [ 1, 1, 4, 8, 15, 15, 24 ]

    set_intrule = -1 ! for testing if no rule available

    select case ( globalshape )

    case ( 'line', 'quadrilateral' )

      if ( inttype == 1 ) then
        set_intrule = nint ( ( order + 3 ) / 2._dp )
      else if ( inttype == 0 .or. inttype == 3 ) then
        set_intrule = nint ( ( order + 1 ) / 2._dp )
      end if

    case ( 'triangle', 'prism' )

      if ( inttype == 0 ) then
        if ( order >= 0 .and. order <= ubound(ntri,1) ) then
          set_intrule = ntri ( order )
        end if
      else if ( inttype == 2 ) then
        set_intrule = nint ( ( order + 1 ) / 2._dp )
      else if ( inttype == 3 ) then
        if ( order == 0 ) then
          set_intrule = 1
        else
          set_intrule = order
        end if
      end if

    case ( 'tetrahedron' )

      if ( inttype == 0 ) then
        if ( order >= 0 .and. order <= ubound(ntet,1) ) then
          set_intrule = ntet ( order )
        end if
      else if ( inttype == 3 ) then
        if ( order == 0 ) then
          set_intrule = 1
        else
          set_intrule = order
        end if
      end if

    case ( 'hexahedron' )

      if ( inttype == 0 ) then
        set_intrule = nint ( ( order + 1 ) / 2._dp )
      else if ( inttype == 1 ) then
        set_intrule = nint ( ( order + 3 ) / 2._dp )
      else if ( inttype == 3 ) then
        if ( GAUSS_RULE == 2 ) then
          if ( order == 0 ) then
            set_intrule = 1
          else
            set_intrule = order
          end if
        else
          set_intrule = nint ( ( order + 1 ) / 2._dp )
        end if
      end if

    case ( 'pyramid' )

      if ( inttype == 2 ) then
        set_intrule = nint ( ( order + 1 ) / 2._dp )
      else if ( inttype == 3 ) then
        if ( order == 0 ) then
          set_intrule = 1
        else
          set_intrule = order
        end if
      end if

    case default

      call errormsg_case_default ( 'set_intrule', 'globalshape', &
        char_value=globalshape )

    end select

!   test whether set_intrule is set

    if ( set_intrule == -1 ) then

      write(*,'(/a/a/2x,a,2(a,i0)/)') &
        'Error in set_intrule: ', &
        '  Gauss integration not available for: ', &
           trim(globalshape),', inttype = ', inttype, ', order = ', order
      stop

    end if

 end function set_intrule


! set intrule2 for specified order of integration

  function set_intrule2 ( globalshape, inttype, order )

    character(len=*), intent(in) :: globalshape
    integer, intent(in) :: inttype, order
    integer :: set_intrule2

    set_intrule2 = 0 ! value if no rule available

    select case ( globalshape )

    case ( 'prism' )

      if ( inttype == 0 .or. inttype == 3 ) then
        set_intrule2 = nint ( ( order + 1 ) / 2._dp )
      end if

    case ( 'pyramid' )

      if ( inttype == 2 ) then
        set_intrule2 = nint ( ( order + 1 ) / 2._dp )
      end if

    case default

      call errormsg_case_default ( 'set_intrule2', 'globalshape', &
        char_value=globalshape )

    end select

  end function set_intrule2


! Set number of Gauss integration points for an element

  subroutine set_ninti ( gauss, ninti, single )

    use set_optional_m

    type(gauss_t), intent(in) :: gauss

!   number of integration points
    integer, intent(out) :: ninti

!   if .true. force the integration rule for a single subdomain
!   default=.false.
    logical, optional, intent(in) :: single

    logical :: lsingle
    integer :: nintis

    lsingle = set_optional ( variable=single, default=.false. )

!   determine nintis for a single domain first

    if ( gauss%inttype == 0 ) then

!     standard Gauss (analytical expressions)

      call set_ninti_standard ( gauss, nintis )

    else if ( gauss%inttype == 1 ) then

!     Lobatto

      call set_ninti_spectral ( gauss, nintis )

    else if ( gauss%inttype == 2 ) then

!     mapped Gauss

      call set_ninti_mapped ( gauss, nintis )

    else if ( gauss%inttype == 3 ) then

!     standard Gauss-Legrendre rules (numerical values)

      call set_ninti_standard_numeric ( gauss, nintis )

    end if


!   set ninti for single of multiple subdomain (macro)

    if ( gauss%nsubint == 1 .or. lsingle ) then

      ninti = nintis

    else

!     multiple subdomains

      call set_ninti_macro ( gauss, nintis, ninti )

    end if

  end subroutine set_ninti


! Gauss rule chooser

  subroutine set_Gauss_integration ( gauss, x, w )

    use limits_m, only: GAUSS_RULE

    type(gauss_t), intent(in) :: gauss

!   the reference coordinates and weights of the integration points
!   x(ninti,ndim), w(ninti)
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w

    real(dp), dimension(:,:), allocatable :: xs
    real(dp), dimension(:), allocatable :: ws

    integer :: m, intrule, intrule2, nintis


!   set some locals

    m = gauss%nsubint
    intrule = gauss%intrule
    intrule2 = gauss%intrule2

    if ( m == 1 ) then

!     single domain

      call set_Gauss_integration_single ( x, w )

    else if ( m > 1 ) then

!     subdomains

      call set_ninti ( gauss, nintis, single=.true. )

!     Select on global shape of the element

      select case ( gauss%globalshape )

      case ( 'line' )

        allocate ( xs(nintis,1), ws(nintis) )

        call set_Gauss_integration_single ( xs, ws )

        call Gauss_Legendre_line_macro ( m, xs(:,1), ws, x(:,1), w )

      case ( 'triangle' )

        allocate ( xs(nintis,2), ws(nintis) )

        call set_Gauss_integration_single ( xs, ws )

        call Gauss_Legendre_triangle_macro ( m, xs, ws, x, w )

      case ( 'quadrilateral' )

        allocate ( xs(nintis,2), ws(nintis) )

        call set_Gauss_integration_single ( xs, ws )

        call Gauss_Legendre_quad_macro ( m, xs, ws, x, w )

      case ( 'hexahedron' )

        allocate ( xs(nintis,3), ws(nintis) )

        call set_Gauss_integration_single ( xs, ws )

        call Gauss_Legendre_hexahedron_macro ( m, xs, ws, x, w )

      case ( 'tetrahedron' )

        allocate ( xs(nintis,3), ws(nintis) )

        call set_Gauss_integration_single ( xs, ws )

        call Gauss_Legendre_tetrahedron_macro ( m, xs, ws, x, w )

      case default

        write(*,'(/a/2a/)') &
          'Error in set_Gauss_integration: ', &
          '  composite Gauss integration (subdomains) not available for: ', &
          gauss%globalshape
        stop

      end select

    else

      write(*,'(/a/a,i0/)') &
        'Error in set_Gauss_integration: ', &
        '  invalid nsubint = ', gauss%nsubint
      stop

    end if

  contains

!   Gauss for a single domain

    subroutine set_Gauss_integration_single ( x, w )

      real(dp), intent(out), dimension(:,:) :: x
      real(dp), intent(out), dimension(:) :: w

      if ( gauss%inttype == 0 ) then

!       standard Gauss

        select case ( gauss%globalshape )

        case ( 'line' )

          call Gauss_Legendre_line ( intrule, x(:,1), w )

        case ( 'triangle' )

          call Gauss_Legendre_triangle ( intrule, x, w )

        case ( 'quadrilateral' )

          call Gauss_Legendre_quad ( intrule, x, w )

        case ( 'hexahedron' )

          call Gauss_Legendre_hexahedron ( intrule, x, w )

        case ( 'tetrahedron' )

          call Gauss_Legendre_tetrahedron ( intrule, x, w )

        case ( 'prism' )

          call Gauss_Legendre_prism ( intrule, intrule2, x, w )

        case default

          write(*,'(/a/2a/)') &
            'Error in set_Gauss_integration: ', &
            '  standard Gauss integration not available for: ', &
            gauss%globalshape
          stop

        end select

      else if ( gauss%inttype == 1 ) then

!       Gauss-Legrendre-Lobatto

        select case ( gauss%globalshape )

        case ( 'line' )

          call Gauss_Legendre_Lobatto_line ( intrule, x(:,1), w )

        case ( 'quadrilateral' )

          call Gauss_Legendre_Lobatto_quad ( intrule, x, w )

        case ( 'hexahedron' )

          call Gauss_Legendre_Lobatto_hexahedron ( intrule, x, w )

        case default

          write(*,'(/a/2a/)') &
            'Error in set_Gauss_integration: ', &
            '  GLL integration not yet implemented for: ', gauss%globalshape
          stop

        end select

      else if ( gauss%inttype == 2 ) then

!       mapped Gauss

        select case ( gauss%globalshape )

        case ( 'triangle' )

          call Gauss_Legendre_mapped_triangle ( intrule, x, w )

        case ( 'pyramid' )

          call Gauss_Legendre_mapped_pyramid ( intrule, intrule2, x, w )

        case default

          write(*,'(/a/2a/)') &
            'Error in set_Gauss_integration: ', &
            '  mapped Gauss integration not available for: ', gauss%globalshape
          stop

        end select

      else if ( gauss%inttype == 3 ) then

!       standard Gauss (numerical values)

        select case ( gauss%globalshape )

        case ( 'line' )

          call Gauss_Legendre_numeric_line ( intrule, x(:,1), w )

        case ( 'triangle' )

          call Gauss_Legendre_numeric_triangle ( intrule, x, w )

        case ( 'quadrilateral' )

          call Gauss_Legendre_numeric_quad ( intrule, x, w )

        case ( 'hexahedron' )

          if ( GAUSS_RULE == 2 ) then

            call Gauss_Legendre_numeric_hexahedron2 ( intrule, x, w )

          else

            call Gauss_Legendre_numeric_hexahedron ( intrule, x, w )

          end if

        case ( 'tetrahedron' )

          call Gauss_Legendre_numeric_tetrahedron ( intrule, x, w )

        case ( 'prism' )

          call Gauss_Legendre_numeric_prism ( intrule, intrule2, x, w )

        case ( 'pyramid' )

          call Gauss_Legendre_numeric_pyramid2 ( intrule, x, w )

        case default

          write(*,'(/a/2a/)') &
            'Error in set_Gauss_integration: ', &
            '  Gauss integration (numerical values) not available for: ', &
            gauss%globalshape
          stop

        end select

      else

        write(*,'(/a/a,i0/)') &
          'Error in set_Gauss_integration: ', &
          '  invalid inttype = ', gauss%inttype
        stop

      end if

    end subroutine set_Gauss_integration_single

  end subroutine set_Gauss_integration


! Gauss rule chooser based on shape of element and number of integration points
! Can only be used for standard Gauss rules.

  subroutine set_Gauss_integration_simple ( globalshape, x, w )

!   Global shape of the element:
!     line
!     quadrilateral
!     triangle
!     hexahedron
!     tetrahedron
    character (len=*), intent(in) :: globalshape

!   the reference coordinates and weights of the integration points
!   x(ninti,ndim), w(ninti)
    real(dp), intent(out), dimension(:,:) :: x
    real(dp), intent(out), dimension(:) :: w


    integer :: ninti, nint1D


    ninti = size(x,1)

!   Select on global shape of the element

    select case ( globalshape )

    case ( 'line' )

      call Gauss_Legendre_line ( ninti, x(:,1), w )

    case ( 'triangle' )

      call Gauss_Legendre_triangle ( ninti, x, w )

    case ( 'quadrilateral' )

      nint1D = nint ( sqrt ( real ( ninti, kind=dp ) ) )

      if ( ninti == nint1D ** 2 ) then
        call Gauss_Legendre_quad ( nint1D, x, w )
      else
        write(*,'(/a/a/)') &
        'Error in set_Gauss_integration_simple: ', &
        ' size(x,1): should be an integer**2 '
        stop
      end if

    case ( 'hexahedron' )

      nint1D = nint ( real ( ninti, kind=dp ) ** (1._dp/3) )

      if ( ninti == nint1D ** 3 ) then
        call Gauss_Legendre_hexahedron ( nint1D, x, w )
      else
        write(*,'(/a/a/)') &
        'Error in set_Gauss_integration_simple: ', &
        ' size(x,1): should be an integer**2 '
        stop
      end if

    case ( 'tetrahedron' )

      call Gauss_Legendre_tetrahedron ( ninti, x, w )

    case default

      write(*,'(/a/2a/)') &
        'Error in set_Gauss_integration_simple: ', &
        '  simple Gauss integration chooser not available for: ', &
        globalshape
      stop

    end select

  end subroutine set_Gauss_integration_simple

end module gauss_m
