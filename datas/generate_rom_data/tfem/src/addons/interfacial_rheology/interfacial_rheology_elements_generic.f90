!
! Copyright (C) 2004-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system (e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.

! Element routines for the interfacial boundary condition:
!              __
!   t| + t|  = \/ . tau
!     d    m     s     s
!
! using interfacial rheological models for tau.
!
! This module contains the generic stuff (works for 2D and 3D).

module interfacial_rheology_elements_generic_m

  use tfem_elem_m

  implicit none

contains

! set global parameters (boundary element)

  subroutine set_globals_interfacial_bc ( mesh, coefficients, curve, &
    surface, volume, ndimr, geometry )

    use interfacial_rheology_globals_m
    use set_optional_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in), optional :: curve, surface, volume
    integer, intent(in), optional :: ndimr, geometry

    integer :: lcurve, lsurface, lvolume

!   check size of coefficients

    call check ( coefficients, 'set_globals_interfacial_bc', ncoefi=100, &
      ncoefr=50, indexarray=[8,9,12], minimum=[0,0,0], &
      maximum=[1,1,1] )

!   traditional interface (legacy)

    lcurve = set_optional ( variable=curve, default=0 )
    lsurface = set_optional ( variable=surface, default=0 )
    lvolume = set_optional ( variable=volume, default=0 )

!   interface with dimension of reference space of geometry

    if ( present(ndimr) .and. present(geometry) ) then
      select case (ndimr)
        case(1); lcurve = geometry
        case(2); lsurface = geometry
        case(3); lvolume = geometry
      case default
        call errormsg_case_default ( 'set_globals_interfacial_bc', &
          'ndimr', int_value=ndimr )
      end select
    end if

    if ( lcurve > 0 ) then
      ndim = mesh%curves(lcurve)%ndim
      if ( ndim == 3 ) then
        coorsys = 2
        vel3D = 0
      else
        coorsys = coefficients%i(8)
        vel3D = coefficients%i(9)
      end if
      globalshape = mesh%curves(lcurve)%element%globalshape
      nodalp = mesh%curves(lcurve)%element%numnod
    else if ( lsurface > 0 ) then
      ndim = mesh%surfaces(lsurface)%ndim
      if ( ndim == 3 ) then
        coorsys = 2
        vel3D = 0
      else
        coorsys = coefficients%i(8)
        vel3D = coefficients%i(9)
      end if
      globalshape = mesh%surfaces(lsurface)%element%globalshape
      nodalp = mesh%surfaces(lsurface)%element%numnod
    else if ( lvolume > 0 ) then
      ndim = mesh%volumes(lvolume)%ndim
      coorsys = 2 ! always 3D
      globalshape = mesh%volumes(lvolume)%element%globalshape
      nodalp = mesh%volumes(lvolume)%element%numnod
    end if

!   set number of degrees of freedom velocity

    intpol = coefficients%i(1)

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(2)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, 'set_globals_interfacial_bc', ndf=ndf )

    if ( ndf /= nodalp ) then
      write(*,'(/5(a/))') 'Error in set_globals_interfacial_bc:', &
        ' Number of degrees of freedom of the shape function (ndf) ', &
        ' is different from the number of nodal points in the element', &
        ' (nodalp) ', &
        ' This possibility (ndf /= nodalp) is not available.'
      stop
    end if

!   set integration

    inttype = coefficients%i(5)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_interfacial_bc:', &
        ' For spectral elements Gauss-Legendre-Lobatto integration ', &
        ' needs to be specified. '
      stop
    end if

    if ( coefficients%i(7) == 1 .or. SET_GAUSS_BY_ORDER ) then
      intrule = set_intrule ( globalshape, inttype, order=coefficients%i(3) )
      intrule2 = set_intrule2 ( globalshape, inttype, order=coefficients%i(6) )
    else
      intrule = coefficients%i(3)
      intrule2 = coefficients%i(6)
    end if
    nsubint = get_coefficient ( coefficients, index=4, default=1 )

    if ( globalshape == 'prism' .and. intrule2 == 0 ) then
      write(*,'(/a/3a/)') 'Error in set_globals_interfacial_bc:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape
      stop
    end if

    if ( globalshape == 'pyramid' .and. inttype == 2 .and. intrule2 == 0 ) then
      write(*,'(/a/a/3a,i0/)') 'Error in set_globals_interfacial_bc:', &
        ' Secondary integration rule (intrule2) needs to be set for ', &
        ' globalshape = ', globalshape, ' and inttype = ', inttype
      stop
    end if

    gauss%globalshape = globalshape
    gauss%intrule = intrule
    gauss%intrule2 = intrule2
    gauss%nsubint = nsubint
    gauss%inttype = inttype

    call set_ninti ( gauss, ninti )

    if ( intpol == 13 .and. ninti /= ndf ) then
      write(*,'(/5(a/))') 'Error in set_globals_interfacial_bc:', &
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_interfacial_bc


! set global parameters

  subroutine set_globals_interfacial_stress ( mesh, coefficients, elgrp )

    use interfacial_rheology_globals_m
    use limits_m, only: SET_GAUSS_BY_ORDER

    type(mesh_t), intent(in) :: mesh
    type(coefficients_t), intent(in) :: coefficients
    integer, intent(in) :: elgrp

    integer :: ndimr

!   check size of coefficients

    call check ( coefficients, 'set_globals_interfacial_stress', ncoefi=100, &
      ncoefr=50, indexarray=[8,9,12], minimum=[0,0,0], maximum=[1,1,1] )

    ndim = mesh%element(elgrp)%ndim
    nodalp = mesh%element(elgrp)%numnod
    globalshape = mesh%element(elgrp)%globalshape

!   set number of degrees of freedom for position

    intpol = coefficients%i(1)

    ndimr = mesh%element(elgrp)%ndimr
    if ( intpol == 13 .and. any(mesh%element(elgrp)%p(:ndimr,2) /= 1 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_interfacial_stress:', &
        ' For spectral elements GLL nodal distribution is required. '
      stop
    else if ( intpol == 20 .and. &
                            any(mesh%element(elgrp)%p(:ndimr,2) /= 0 ) ) then
      write(*,'(/2(a/))') 'Error in set_globals_interfacial_stress:', &
        ' For high-order elements equidistant nodal distribution is required. '
      stop
    end if

    shapefunc%globalshape = globalshape
    shapefunc%interpolation = intpol
    shapefunc%numbering = 'standard'
    shapefunc%p = coefficients%i(2)
    shapefunc%spec_eval = 'gauss'

    call set_ndf ( shapefunc, 'set_globals_interfacial_stress', ndf=ndf )

!   set integration

    inttype = coefficients%i(5)

    if ( intpol == 13 .and. inttype /= 1 ) then
      write(*,'(/3(a/))') 'Error in set_globals_interfacial_stress:', &
        ' For spectral elements Gauss-Legendre-Lobatto integration ', &
        ' needs to be specified. '
      stop
    end if

    if ( coefficients%i(7) == 1 .or. SET_GAUSS_BY_ORDER ) then
      intrule = set_intrule ( globalshape, inttype, order=coefficients%i(3) )
    else
      intrule = coefficients%i(3)
    end if
    nsubint = get_coefficient ( coefficients, index=4, default=1 )

    gauss%globalshape = globalshape
    gauss%intrule = intrule
    gauss%nsubint = nsubint
    gauss%inttype = inttype

    call set_ninti ( gauss, ninti )

    if ( intpol == 13 .and. ninti /= ndf ) then
      write(*,'(/5(a/))') 'Error in set_globals_interfacial_stress:', &
        ' Number of integration points (ninti) is different from ', &
        ' the number of degrees of freedom (ndf) ', &
        ' This possibility (ninti /= ndf) is not available if ', &
        ' spectral interpolation is used'
      stop
    end if

  end subroutine set_globals_interfacial_stress

end module interfacial_rheology_elements_generic_m
