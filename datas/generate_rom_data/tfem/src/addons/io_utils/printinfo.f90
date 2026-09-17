! Copyright (C) 2007-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Print some info of a structure to standard output

module printinfo_m

  use kind_defs_m
  use mesh_m
  use problem_defs_m
  use set_optional_m

  implicit none

! interface to support generic interface

  interface printinfo
    module procedure printinfo_mesh, printinfo_problem
  end interface printinfo

contains


! write mesh data to standard output

  subroutine printinfo_mesh ( mesh, printlevel )

    type(mesh_t), intent(in) :: mesh

!   printlevel:
!    0: no printing
!    1: minimal printing on standard output
!    2: some more printing
!    3: like 2 but include element shape numbers
!    4: like 3 but include element shape numbers curves and surfaces
!    5: like 4 but include element order
!    6: like 5 but include element order for curves and surfaces
!   default=2
    integer, intent(in), optional :: printlevel


    integer :: lprintlevel, m, i, curve, surface, volume

!   testing

    if ( .not. mesh%meshgen ) then
      write(*,'(/a/)') &
        'Error printinfo_mesh: mesh has not been generated '
      stop
    end if

    lprintlevel = set_optional ( variable=printlevel, default=2 )

    if ( lprintlevel >= 1 ) then

!     basic data

      write(unit=*,fmt='(a/)') 'Mesh info:'

      write(unit=*,fmt='(a,i0)') &
        ' Space dimension (ndim)                = ', mesh%ndim, &
        ' Number of nodes (nnodes)              = ', mesh%nnodes, &
        ' Number of elements (nelem)            = ', mesh%nelem, &
        ' Number of element groups (nelgrp)     = ', mesh%nelgrp, &
        ' Number of blend meshes (nblend)       = ', mesh%nblend

      if ( mesh%nblend > 0 ) then
        write(unit=*,fmt='(a)') ' Number of nodes in each blend mesh:'
        write(unit=*,fmt='((2x,10(i0,1x)))') &
          ( mesh%nnodes_blend(m+2)-mesh%nnodes_blend(m+1), m = 0, mesh%nblend )
      end if

    end if

    if ( lprintlevel >= 2 ) then

!     more data

      write(unit=*,fmt='(a,i0)') &
        ' Number of points (npoints)            = ', mesh%npoints, &
        ' Number of curves (ncurves)            = ', mesh%ncurves, &
        ' Number of surfaces (nsurfaces)        = ', mesh%nsurfaces, &
        ' Number of volumes (nvolumes)          = ', mesh%nvolumes, &
        ' Number of objects (nobjects)          = ', mesh%nobjects, &
        ' Number of nodesets (nnodesets)        = ', mesh%nnodesets, &
        ' Number of nelementsets (nelementsets) = ', mesh%nelementsets, &
        ' Number of blocks (nblocks)            = ', mesh%nblocks, &
        ' Number of gluepoints (ngluepoints)    = ', mesh%ngluepoints

      write(unit=*,fmt='(4(a,l1))') &
        ' meshparts = ', mesh%meshparts, &
        ' renumber = ', mesh%renumber, &
        ' sblocks = ', mesh%sblocks%filled, &
        ' nodblocks = ', mesh%nodblocks%filled

      if ( mesh%meshparts ) then
        write(unit=*,fmt='(a)') ' Number of nodes per element in each group:'
        write(unit=*,fmt='((2x,10(i0,1x)))') mesh%elnumnod
        write(unit=*,fmt='(a,i0)') &
        ' maxnodnumnod = ', mesh%maxnodnumnod, &
        ' Size of mesh%nodnod = ', mesh%nodnumnod(mesh%nnodes+1)
      end if

      if ( mesh%nelgrp > 1 ) then
        write(unit=*,fmt='(a)') ' Number of elements in each element group:'
        write(unit=*,fmt='((2x,10(i0,1x)))') mesh%grpnumel
      end if

      if ( mesh%ncurves > 0 ) then
        write(unit=*,fmt='(a)') ' Number of nodes in each curve:'
        write(unit=*,fmt='((2x,10(i0,1x)))') mesh%curves(1:mesh%ncurves)%nnodes
        if ( mesh%nblend > 0 ) then
          write(unit=*,fmt='(a)') &
                        ' Number of nodes in the curves for each blend mesh:'
          do curve = 1, mesh%ncurves
            write(unit=*,fmt='((2x,a,i0,a,1x,10(i0,1x)))') &
            'curve = ', curve, ':', ( mesh%curves(curve)%nnodes_blend(m+2) &
                 - mesh%curves(curve)%nnodes_blend(m+1), m = 0, mesh%nblend )
          end do
        end if
        write(unit=*,fmt='(a)') ' Number of elements in each curve:'
        write(unit=*,fmt='((2x,10(i0,1x)))') mesh%curves(1:mesh%ncurves)%nelem
        if ( mesh%meshparts ) then
          write(unit=*,fmt='(a)') ' Number of nodes per element in each curve:'
          write(unit=*,fmt='((2x,10(i0,1x)))') &
                               mesh%curves(1:mesh%ncurves)%elnumnod
        end if
      end if

      if ( mesh%nsurfaces > 0 ) then
        write(unit=*,fmt='(a)') ' Number of nodes in each surface:'
        write(unit=*,fmt='((2x,10(i0,1x)))') &
                                mesh%surfaces(1:mesh%nsurfaces)%nnodes
        if ( mesh%nblend > 0 ) then
          write(unit=*,fmt='(a)') &
                     ' Number of nodes in the surfaces for each blend mesh:'
          do surface = 1, mesh%nsurfaces
            write(unit=*,fmt='((2x,a,i0,a,1x,10(i0,1x)))') &
            'surface = ', surface, ':', &
                ( mesh%surfaces(surface)%nnodes_blend(m+2) &
                 - mesh%surfaces(surface)%nnodes_blend(m+1), &
                          m = 0, mesh%nblend )
          end do
        end if
        write(unit=*,fmt='(a)') ' Number of elements in each surface:'
        write(unit=*,fmt='((2x,10(i0,1x)))')  &
                                mesh%surfaces(1:mesh%nsurfaces)%nelem
        if ( mesh%meshparts ) then
          write(unit=*,fmt='(a)') &
                             ' Number of nodes per element in each surface:'
          write(unit=*,fmt='((2x,10(i0,1x)))') &
                               mesh%surfaces(1:mesh%nsurfaces)%elnumnod
        end if
      end if

      if ( mesh%nvolumes > 0 ) then
        write(unit=*,fmt='(a)') ' Number of nodes in each volume:'
        write(unit=*,fmt='((2x,10(i0,1x)))') &
                                mesh%volumes(1:mesh%nvolumes)%nnodes
        if ( mesh%nblend > 0 ) then
          write(unit=*,fmt='(a)') &
                        ' Number of nodes in the volumes for each blend mesh:'
          do volume = 1, mesh%nvolumes
            write(unit=*,fmt='((2x,a,i0,a,1x,10(i0,1x)))') &
            'volume = ', volume, ':', ( mesh%volumes(volume)%nnodes_blend(m+2) &
                 - mesh%volumes(volume)%nnodes_blend(m+1), m = 0, mesh%nblend )
          end do
        end if
        write(unit=*,fmt='(a)') ' Number of elements in each volume:'
        write(unit=*,fmt='((2x,10(i0,1x)))')  &
                                mesh%volumes(1:mesh%nvolumes)%nelem
        if ( mesh%meshparts ) then
          write(unit=*,fmt='(a)') ' Number of nodes per element in each volume:'
          write(unit=*,fmt='((2x,10(i0,1x)))') &
                               mesh%volumes(1:mesh%nvolumes)%elnumnod
        end if
      end if

    end if

    if ( lprintlevel >= 3 ) then

!     print element shapes

      if ( mesh%nelgrp > 1 ) then
        write(unit=*,fmt='(a)') ' Element shape of each element group:'
        write(unit=*,fmt='((2x,25(i0,1x)))') mesh%element(:)%elshape
        do m = 1, mesh%nblend
           write(unit=*,fmt='(a,i0,a/2x,(25(i0,1x)))') ' Blend = ', m, &
                    '. Element shape of each group: ', &
                     mesh%element_blend(:,m)%elshape
        end do
      else if ( mesh%nelgrp == 1 ) then
        write(unit=*,fmt='(a,i0)') ' Element shape = ', &
                                   mesh%element(:)%elshape
        do m = 1, mesh%nblend
           write(unit=*,fmt='(a,i0,a,(25(i0,1x)))') ' Blend = ', m, &
                    '. Element shape = ', mesh%element_blend(:,m)%elshape
        end do
      end if

    end if

    if ( lprintlevel >= 4 ) then

!     print element shapes in elements of curves

      if ( mesh%ncurves > 0 ) then
        write(unit=*,fmt='(a)') ' Element shapes of elements in each curve:'
        write(unit=*,fmt='((2x,35(i0,1x)))') &
                          mesh%curves(1:mesh%ncurves)%element%elshape
        do m = 1, minval(mesh%curves(1:mesh%ncurves)%nblend)
           write(unit=*,fmt='(a,i0,a/2x,(25(i0,1x)))') ' Blend = ', m, &
                    '. Element shape in each curve: ', &
                  ( mesh%curves(i)%element_blend(m)%elshape, i=1,mesh%ncurves )
        end do
      end if

!     print element shapes in elements of surfaces

      if ( mesh%nsurfaces > 0 ) then
        write(unit=*,fmt='(a)') ' Element shapes of elements in each surface:'
        write(unit=*,fmt='((2x,35(i0,1x)))')  &
                          mesh%surfaces(1:mesh%nsurfaces)%element%elshape
        do m = 1, minval(mesh%surfaces(1:mesh%nsurfaces)%nblend)
           write(unit=*,fmt='(a,i0,a/2x,(25(i0,1x)))') ' Blend = ', m, &
                    '. Element shape in each surface: ', &
             ( mesh%surfaces(i)%element_blend(m)%elshape,i=1,mesh%nsurfaces )
        end do
      end if

!     print element shapes in elements of volumes

      if ( mesh%nvolumes > 0 ) then
        write(unit=*,fmt='(a)') ' Element shapes of elements in each volume:'
        write(unit=*,fmt='((2x,35(i0,1x)))')  &
                          mesh%volumes(1:mesh%nvolumes)%element%elshape
        do m = 1, minval(mesh%volumes(1:mesh%nvolumes)%nblend)
           write(unit=*,fmt='(a,i0,a/2x,(25(i0,1x)))') ' Blend = ', m, &
                    '. Element shape in each volume: ', &
             ( mesh%volumes(i)%element_blend(m)%elshape,i=1,mesh%nvolumes )
        end do
      end if

    end if

    if ( lprintlevel >= 5 ) then

!     print element order

      if ( mesh%nelgrp > 1 ) then
        write(unit=*,fmt='(a)') ' Element order in each element group:'
        write(unit=*,fmt='((2x,a,25(i0,1x)))') 'p ', mesh%element(:)%p(1,1)
        write(unit=*,fmt='((2x,a,25(i0,1x)))') 'q ', mesh%element(:)%p(2,1)
        write(unit=*,fmt='((2x,a,25(i0,1x)))') 'r ', mesh%element(:)%p(3,1)
        do m = 1, mesh%nblend
           write(unit=*,fmt='(a,i0,a)') ' Blend = ', m, &
                    '. Element order in each group: '
           write(unit=*,fmt='((2x,a,25(i0,1x)))') 'p ', &
                                           mesh%element_blend(:,m)%p(1,1)
           write(unit=*,fmt='((2x,a,25(i0,1x)))') 'q ', &
                                           mesh%element_blend(:,m)%p(2,1)
           write(unit=*,fmt='((2x,a,25(i0,1x)))') 'r ', &
                                           mesh%element_blend(:,m)%p(3,1)
        end do
      else if ( mesh%nelgrp == 1 ) then
        write(unit=*,fmt='(a,3(i0,1x))') ' Element order p q r = ', &
          mesh%element(:)%p(1,1), mesh%element(:)%p(2,1), &
          mesh%element(:)%p(3,1)
        do m = 1, mesh%nblend
          write(unit=*,fmt='(a,i0,a,3(i0,1x))') ' Blend = ', m, &
                    '. Element order p q r = ', &
            mesh%element_blend(:,m)%p(1,1), mesh%element_blend(:,m)%p(2,1), &
            mesh%element_blend(:,m)%p(3,1)
        end do
      end if

    end if

    if ( lprintlevel >= 6 ) then

!     print element order in elements of curves

      if ( mesh%ncurves > 0 ) then
        write(unit=*,fmt='(a)') ' Element order of elements in each curve:'
        write(unit=*,fmt='((2x,35(i0,1x)))') &
                          mesh%curves(1:mesh%ncurves)%element%p(1,1)
        do m = 1, minval(mesh%curves(1:mesh%ncurves)%nblend)
           write(unit=*,fmt='(a,i0,a/2x,(35(i0,1x)))') ' Blend = ', m, &
                    '. Element order in each curve: ', &
                  ( mesh%curves(i)%element_blend(m)%p(1,1), i=1,mesh%ncurves )
        end do
      end if

!     print element order in elements of surfaces

      if ( mesh%nsurfaces > 0 ) then
        write(unit=*,fmt='(a)') ' Element order of elements in each surface:'
        write(unit=*,fmt='((2x,a,35(i0,1x)))') 'p ', &
                               mesh%surfaces(1:mesh%nsurfaces)%element%p(1,1)
        write(unit=*,fmt='((2x,a,35(i0,1x)))') 'q ', &
                               mesh%surfaces(1:mesh%nsurfaces)%element%p(2,1)
        do m = 1, minval(mesh%surfaces(1:mesh%nsurfaces)%nblend)
           write(unit=*,fmt='(a,i0,a)') ' Blend = ', m, &
                    '. Element order in each surface: '
           write(unit=*,fmt='((2x,a,35(i0,1x)))') 'p ', &
               ( mesh%surfaces(i)%element_blend(m)%p(1,1), i=1,mesh%nsurfaces )
           write(unit=*,fmt='((2x,a,35(i0,1x)))') 'q ', &
               ( mesh%surfaces(i)%element_blend(m)%p(2,1), i=1,mesh%nsurfaces )
        end do
      end if

!     print element order in elements of volumes

      if ( mesh%nvolumes > 0 ) then
        write(unit=*,fmt='(a)') ' Element order of elements in each volume:'
        write(unit=*,fmt='((2x,a,35(i0,1x)))') 'p ', &
                               mesh%volumes(1:mesh%nvolumes)%element%p(1,1)
        write(unit=*,fmt='((2x,a,35(i0,1x)))') 'q ', &
                               mesh%volumes(1:mesh%nvolumes)%element%p(2,1)
        write(unit=*,fmt='((2x,a,35(i0,1x)))') 'r ', &
                               mesh%volumes(1:mesh%nvolumes)%element%p(3,1)
        do m = 1, minval(mesh%volumes(1:mesh%nvolumes)%nblend)
           write(unit=*,fmt='(a,i0,a)') ' Blend = ', m, &
                    '. Element order in each volume: '
           write(unit=*,fmt='((2x,a,35(i0,1x)))') 'p ', &
               ( mesh%volumes(i)%element_blend(m)%p(1,1), i=1,mesh%nvolumes )
           write(unit=*,fmt='((2x,a,35(i0,1x)))') 'q ', &
               ( mesh%volumes(i)%element_blend(m)%p(2,1), i=1,mesh%nvolumes )
           write(unit=*,fmt='((2x,a,35(i0,1x)))') 'r ', &
               ( mesh%volumes(i)%element_blend(m)%p(3,1), i=1,mesh%nvolumes )
        end do
      end if

    end if

    if ( lprintlevel >= 1 ) write(unit=*,fmt=*)

  end subroutine printinfo_mesh


! write problem data to standard output

  subroutine printinfo_problem ( problem, printlevel )

    type(problem_t), intent(in) :: problem

!   printlevel:
!    0: no printing
!    1: minimal printing on standard output
!    2: some more printing
!    3: like 2 but include number of degrees for the vectors
!    4: like 3 but include element definition of the sysvectors and vectors and
!       info on physical quantities and inactive/active groups
!   default=2
    integer, intent(in), optional :: printlevel


    integer :: lprintlevel, vec, elgrp

!   testing

    call check ( problem, 'printinfo_problem' )

    lprintlevel = set_optional ( variable=printlevel, default=2 )

    if ( lprintlevel >= 1 ) then

!     basic data

      write(unit=*,fmt='(a/)') 'Problem info:'

      if ( problem%probnr > 0 ) then
        write(unit=*,fmt='(a,18x,a,i0)') &
        ' Problem number (probnr)                ', ' = ', problem%probnr
      end if
      write(unit=*,fmt='(a,18x,a,i0)') &
        ' Number of vectors (nvec)               ', ' = ', problem%nvec, &
        ' Number of physical quantities (nphysq) ', ' = ', problem%nphysq, &
        ' Number of physqshifted (nphysqshifted) ', ' = ', &
                                                      problem%nphysqshifted, &
        ' Number of constraints (numconstraints) ', ' = ',  &
                                                     problem%numconstraints, &
        ' Number of degrees of freedom (numdegfd)', ' = ', problem%numdegfd

    end if

    if ( lprintlevel >= 2 ) then

!     more data

      write(unit=*,fmt='(a,1x,a,i0)') &
        ' Number of nodal degrees of freedom (numnodaldegfd)      ', ' = ',  &
                                                     problem%numnodaldegfd, &
        ' Number of constraint degrees of freedom (numconstrdegfd)', ' = ',  &
                                                     problem%numconstrdegfd, &
        ' Number of unknown degrees of freedom (numundegfd)       ', ' = ',  &
                                                     problem%numundegfd, &
        ' Number of essential degrees of freedom (numessdegfd)    ', ' = ',  &
                                                     problem%numessdegfd, &
        ' Number of inactive groups (numinactivegroups)           ', ' = ', &
                                                     problem%numinactivegroups
    end if

    if ( lprintlevel >= 3 ) then

!     number of degrees in vectors

      write(unit=*,fmt='(a)') ' Number of degrees for each vector (nodalwise): '
      write(unit=*,fmt='((2x,25(i0,1x)))') problem%vec_numdegfd(:,1)
      write(unit=*,fmt='(a)')&
                          ' Number of degrees for each vector (elementwise): '
      write(unit=*,fmt='((2x,25(i0,1x)))') problem%vec_numdegfd(:,2)

    end if

    if ( lprintlevel >= 4 ) then

!     element definition of degrees in sysvectors and vectors

      write(unit=*,fmt='(a)') ' Element definition for the sysvector: '
      if ( problem%nelgrp > 1 ) then
        do elgrp = 1, problem%nelgrp
          write(unit=*,fmt='(a,i0)') ' Element group ', elgrp
          write(unit=*,fmt='((2x,25(i0,1x)))') problem%elnumdegfd(elgrp)%a
        end do
      else if ( problem%nelgrp == 1 ) then
        write(unit=*,fmt='((2x,25(i0,1x)))') problem%elnumdegfd(1)%a
      end if

      write(unit=*,fmt='(a)') ' Element definition for the vectors: '
      if ( problem%nelgrp > 1 ) then
        do elgrp = 1, problem%nelgrp
          write(unit=*,fmt='(a,i0)') ' Element group ', elgrp
          do vec = 1, problem%nvec
            write(unit=*,fmt='(a,i0)') '  Vector ', vec
            write(unit=*,fmt='((2x,25(i0,1x)))') &
                                          problem%vec_elnumdegfd(elgrp)%a(:,vec)
          end do
        end do
      else if ( problem%nelgrp == 1 ) then
        do vec = 1, problem%nvec
          write(unit=*,fmt='(a,i0)') '  Vector ', vec
          write(unit=*,fmt='((2x,25(i0,1x)))') &
                                        problem%vec_elnumdegfd(1)%a(:,vec)
        end do
      end if

      if ( problem%nphysq > 0 ) then
        write(unit=*,fmt='(a)') ' Physical quantities defined by vectors: '
          write(unit=*,fmt='((2x,25(i0,1x)))') problem%physq
      end if

      if ( problem%nphysqshifted > 0 ) then
        write(unit=*,fmt='(a)') ' Physical quantities shifted: '
          write(unit=*,fmt='((2x,25(i0,1x)))') problem%physqshifted
      end if

      if ( problem%numinactivegroups > 0 ) then
        write(unit=*,fmt='(a)') ' Active groups: '
          write(unit=*,fmt='((2x,25(i0,1x)))') problem%activegroups
        write(unit=*,fmt='(a)') ' Inactive groups: '
          write(unit=*,fmt='((2x,25(i0,1x)))') problem%inactivegroups
      end if

    end if

    if ( lprintlevel >= 1 ) write(unit=*,fmt=*)

  end subroutine printinfo_problem

end module printinfo_m
