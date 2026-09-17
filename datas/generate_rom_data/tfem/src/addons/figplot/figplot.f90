
! Copyright (C) 2004-2021 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.
!
! plotting routines producing Fig format output
!

module figplot_m

  use kind_defs_m
  use system_m
  use vector_m
  use mesh_m, only: mesh_t
  use problem_defs_m
  use limits_m, only: MAXPOLYLINEPOINTS
  use meshgen_parts_m


  implicit none


  integer, parameter  :: MMFIGUNITS = 45  ! number of figunits in 1 mm

  integer :: unit_figplot = 13 ! unit for writing


! A basic assumption underlies the plot_color_fill, plot_color_contour and
! plot_vector routines:
!
!   ALL NODES MUST CONTAIN DATA OF THE QUANTITY TO BE PLOTTED
!
! If a quantity violates this, a new derived vector should be created (using the
! routine derive_vector), that defines the quantity in all nodal points.
! Examples are: pressures in velocity/pressure elements, a seven-node triangle
! where the velocities are defined in the six boundary nodes and the
! discontinuous stresses in the center node (node 7). In the latter example
! both the velocities and stresses need derived vectors for plotting.


! type definition of the options for plotting

  type plot_options_t

!   the size of the figure in mm of the data region or the mesh
    real(dp) :: figsize = 150

!   Increase the maximum and minimum coordinate values in 2D or the
!   projected plot in 3D by the amount given. The reason for including these
!   values is to be able to `zoom out' such that objects are fully plotted
!   even if the coordinates of the objects are outside the window of the mesh.
    real(dp) :: xmax_coor_add = 0._dp
    real(dp) :: xmin_coor_add = 0._dp
    real(dp) :: ymax_coor_add = 0._dp
    real(dp) :: ymin_coor_add = 0._dp

!   Maximum and minimum coordinate values in 2D or the projected plot in 3D.
!   Setting (one or more of) these values creates the possibility to
!   `zoom in' en plot only part of the region.
!   NOTE, this is not a clipping of the plot: any plot object or element
!   that has at least one of it's coordinates outside this region will be
!   fully discarded.
    real(dp) :: xmax = 1.0e+100_dp
    real(dp) :: xmin = -1.0e+100_dp
    real(dp) :: ymax = 1.0e+100_dp
    real(dp) :: ymin = -1.0e+100_dp

!   fontsize used in text that is printed. A side effect is the scaling of
!   the size of the levels block in a contour plot: by decreasing the fontsize
!   more levels fit into to plot.
    integer :: fontsize = 20

!   add labels like C1, P1 etc. to curves and points
    logical :: printlabels = .true.

!   use a spline instead of a polyline for plotting elements in curves and
!   boundaries of internal elements. Note, that this only applies to plotting
!   curves and meshes
    logical :: spline = .false.

!   scaling of vectors relative to unit of the physical coordinates x and y
    real(dp) :: scalevector = 0.1_dp

!   the size of the points of the objects in mm
    real(dp) :: objectpointsize = 0.2

!   color of vectors, mesh, boundary in mesh, objects, points in objects
!   from manual of xfig:
!      -1 = Default
!       0 = Black
!       1 = Blue
!       2 = Green
!       3 = Cyan
!       4 = Red
!       5 = Magenta
!       6 = Yellow
!       7 = White
!    8-11 = four shades of blue (dark to lighter)
!   12-14 = three shades of green (dark to lighter)
!   15-17 = three shades of cyan (dark to lighter)
!   18-20 = three shades of red (dark to lighter)
!   21-23 = three shades of magenta (dark to lighter)
!   24-26 = three shades of brown (dark to lighter)
!   27-30 = four shades of pink (dark to lighter)
!      31 = Gold
    integer :: vectorcolor = 4, meshcolor = 0
    integer :: boundarycolor = 1, objectpointcolor = 0

!   number of colors used in the filled color plot
!   note that the number of colors must be 4*i/4 + 1, i=1,2,...
!   thus only for numcolors=5, 9, 13, 17, 21, ... this number is exact,
!   otherwise it is rounded.
    integer :: numcolors = 13

!   number of levels used in the contour plots
!   note that when colored contours are used, the number of colors must
!   be 4*i/4 + 1, i=1,2,...  thus only for numlevels=5, 9, 13, 17, 21, ...
!   this number is exact, otherwise it is rounded
    integer :: numlevels = 13

!   use floating point notation instead of exponential
    logical :: usefloating = .false.

!   plot complete outside boundary in the following cases:
!     1) if only part of the element groups is plotted
!     2) in routine plot_objects
    logical :: plotboundary = .true.

!   plot outside boundary for the elements plotted in the following cases:
!     1) plot_color_fill
!     2) plot_color_contour
!     3) plot_vector
    logical :: plotboundary2 = .true.

!   plot a colorbar
    logical :: plotcolorbar = .true.

!   plot levelvalues of contours
    logical :: plotlevelvalues = .true.

!   print filename in title
    logical :: printfilename = .true.

!   print date and time in title
    logical :: printdateandtime = .true.

!   viewpoint in case of 3D
!   For example viewpoint=(/1,1,1/) means that every point will be projected
!   on the plane with a normal in the direction of (/1,1,1/) and the view
!   is _from_ (/1,1,1/).
    real(dp), dimension(3) :: viewpoint = [ 1._dp, 1._dp, 0.8_dp ]

!   ydirection in case of 3D
!   co-ordinate direction that becomes the new y co-oordinate in the projected
!   plane, 1=x, 2=y, 3=z.
    integer :: ydirection = 3


!   options for the colorbar and/or levelvalues:
!   --------------------------------------------

!   print the maximum and the minimum values
    logical :: printmaxminvalues = .true.

!   maximum and minimum value to be used in the colors
    real(dp) :: maxvalue = 1.0e+100_dp
    real(dp) :: minvalue = -1.0e+100_dp

!   if forcemaxmin = .true. the maximum/minimum values above are always used
!   in the colors independent of the actual values. In this way the colors are
!   static and can be used in movies, for example. Don't forgot to specify
!   proper values for maxvalue and minvalue.
!   if forcemaxmin = .false. the maximum/minimum values above are used in the
!   colors only when the actual maximum (minimum) is larger (smaller).
!   The colors are dynamic and all the colors are distributed between the
!   actual maximum and minimum when maxvalue of minvalue is not crossed.
    logical :: forcemaxmin = .false.

!   shift the colorbar/levelvalues in x and y direction. values are in mm in
!   the actual print of the plot.
    real(dp) :: shiftbarx = 0.0_dp, shiftbary = 0.0_dp

!   number of numerical value intervals. (colorbar only)
!   number of labels is numvalueintervals + 1
    integer :: numvalueintervals = 5

!   rotate the colorbar 90 degree (colorbar only)
    logical :: rotatecolorbar = .false.

  end type plot_options_t


contains


! Helper routine for setting the plot options

  subroutine set_plot_options ( plot_options, keep, fontsize, vectorcolor, &
    meshcolor, boundarycolor, objectpointcolor, numcolors, numlevels, &
    numvalueintervals, figsize, scalevector, objectpointsize, maxvalue, &
    minvalue, shiftbarx, shiftbary, usefloating, plotcolorbar, &
    plotboundary, plotlevelvalues, printfilename, printdateandtime, &
    printmaxminvalues, forcemaxmin, rotatecolorbar, viewpoint, ydirection, &
    printlabels, xmin, xmax, ymin, ymax, xmin_coor_add, xmax_coor_add, &
    ymin_coor_add, ymax_coor_add, plotboundary2, spline )

!   Options for plotting
!   See type definition for possibilities and defaults.
    type(plot_options_t), intent(inout) :: plot_options

!   if .true. input values are kept. default=.false.
    logical, intent(in), optional :: keep

!   parameters: see plot_options_t

    integer, intent(in), optional :: fontsize, vectorcolor, meshcolor, &
      boundarycolor, objectpointcolor, numcolors, numlevels, &
      numvalueintervals, ydirection

    real(dp), intent(in), optional :: figsize, scalevector, objectpointsize, &
      maxvalue, minvalue, shiftbarx, shiftbary, viewpoint(3), xmin, xmax, &
      ymin, ymax, xmin_coor_add, xmax_coor_add, ymin_coor_add, ymax_coor_add

    logical, intent(in), optional :: usefloating, plotcolorbar, plotboundary, &
      plotlevelvalues, printfilename, printdateandtime, printmaxminvalues, &
      forcemaxmin, rotatecolorbar, printlabels, plotboundary2, spline

    type(plot_options_t) :: plotopt


    if ( present(keep) ) then
      if ( .not. keep ) plot_options = plotopt
    else
      plot_options = plotopt
    end if

    if ( present(fontsize) ) plot_options%fontsize = fontsize
    if ( present(vectorcolor) ) plot_options%vectorcolor = vectorcolor
    if ( present(meshcolor) ) plot_options%meshcolor = meshcolor
    if ( present(boundarycolor) ) plot_options%boundarycolor = boundarycolor
    if ( present(objectpointcolor) ) &
      plot_options%objectpointcolor = objectpointcolor
    if ( present(numcolors) ) plot_options%numcolors = numcolors
    if ( present(numlevels) ) plot_options%numlevels = numlevels
    if ( present(numvalueintervals) ) &
      plot_options%numvalueintervals = numvalueintervals
    if ( present(figsize) ) plot_options%figsize = figsize
    if ( present(scalevector) ) plot_options%scalevector = scalevector
    if ( present(objectpointsize) ) &
      plot_options%objectpointsize = objectpointsize
    if ( present(maxvalue) ) plot_options%maxvalue = maxvalue
    if ( present(minvalue) ) plot_options%minvalue = minvalue
    if ( present(shiftbarx) ) plot_options%shiftbarx = shiftbarx
    if ( present(shiftbary) ) plot_options%shiftbary = shiftbary
    if ( present(usefloating) ) plot_options%usefloating = usefloating
    if ( present(plotcolorbar) ) plot_options%plotcolorbar = plotcolorbar
    if ( present(plotboundary) ) plot_options%plotboundary = plotboundary
    if ( present(plotboundary2) ) plot_options%plotboundary2 = plotboundary2
    if ( present(plotlevelvalues) ) &
      plot_options%plotlevelvalues = plotlevelvalues
    if ( present(printfilename) ) plot_options%printfilename = printfilename
    if ( present(printdateandtime) )  &
     plot_options%printdateandtime = printdateandtime
    if ( present(printmaxminvalues) )  &
     plot_options%printmaxminvalues = printmaxminvalues
    if ( present(forcemaxmin) ) plot_options%forcemaxmin = forcemaxmin
    if ( present(rotatecolorbar) ) plot_options%rotatecolorbar = rotatecolorbar
    if ( present(viewpoint) ) plot_options%viewpoint = viewpoint
    if ( present(ydirection) ) plot_options%ydirection = ydirection
    if ( present(printlabels) ) plot_options%printlabels = printlabels
    if ( present(spline) ) plot_options%spline = spline
    if ( present(xmin) ) plot_options%xmin = xmin
    if ( present(xmax) ) plot_options%xmax = xmax
    if ( present(ymin) ) plot_options%ymin = ymin
    if ( present(ymax) ) plot_options%ymax = ymax
    if ( present(xmin_coor_add) ) plot_options%xmin_coor_add = xmin_coor_add
    if ( present(xmax_coor_add) ) plot_options%xmax_coor_add = xmax_coor_add
    if ( present(ymin_coor_add) ) plot_options%ymin_coor_add = ymin_coor_add
    if ( present(ymax_coor_add) ) plot_options%ymax_coor_add = ymax_coor_add

  end subroutine set_plot_options


! plot points and curves (and surface labels when present)

  subroutine plot_points_curves ( plot_options, mesh, filename, append, blend )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   when .true.: append plot to existing file. Default=.false.
    logical, intent(in), optional :: append

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
!   default=0
    integer, intent(in), optional :: blend


    real(dp), dimension(3) :: xunitvec, yunitvec
    type(mesh_t) :: lmesh


!   test dimension of mesh

    if ( mesh%ndim == 3 ) then

!     3D mesh: project coordinates onto a 2D plane.

!     compute projection vectors

      call projection_vectors ( &
        plot_options%viewpoint, plot_options%ydirection, xunitvec, yunitvec )

!     copy mesh to local mesh

      call copy ( mesh, lmesh )

!     project coordinates in new mesh

      call project_coordinates ( xunitvec, yunitvec, lmesh%coor )

      call plot_points_curves_2D ( plot_options, lmesh, filename, append, &
        blend )

      call delete ( lmesh )

    else

      call plot_points_curves_2D ( plot_options, mesh, filename, append, &
        blend )

    end if

  end subroutine plot_points_curves


! plot points and curves (and surface labels when present) (2D)

  subroutine plot_points_curves_2D ( plot_options, mesh, filename, append, &
    blend )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   when .true.: append plot to existing file. Default=.false.
    logical, intent(in), optional :: append

!   the mesh blend is used for the coordinates and elements:
!     blend = 0 is the main mesh with element shapes defined by mesh%element
!     blend > 0 is the blend mesh with element shapes defined by
!               mesh%element_blend
!   default=0
    integer, intent(in), optional :: blend


    logical :: fileexists, fileappend
    character (len=4) :: ch
    integer :: npoints, curve, point, surface, elem, lblend, i, m
    real(dp) :: x(2,MAXPOLYLINEPOINTS)
    real(dp) :: xc, yc
    real(dp) :: scalefac, xmin, xmax, ymin, ymax
    real(dp) :: xminc, xmaxc, yminc, ymaxc
    real(dp) :: xmins, xmaxs, ymins, ymaxs
    type(element_t), dimension(mesh%ncurves) :: elementc
    integer, dimension(mesh%ncurves,mesh%nblend+2) :: numnodtop


    inquire ( file=filename, exist=fileexists )

    if ( present(append) ) then
      fileappend = append
    else
      fileappend = .false.
    end if

    if ( fileexists .and. fileappend ) then

      open ( unit=unit_figplot, file=filename, status='old', position='append' )

    else

      open ( unit=unit_figplot, file=filename )
      call fig_header
      call print_head ( plot_options, filename )

    end if

    if ( present(blend) ) then
      if ( blend < 0 .or. blend > mesh%nblend ) then
        write(*,'(/a/a,i0/)') &
          'Error plot_points_curves: invalid argument ', &
          '  blend < 0 or larger than ', mesh%nblend
        stop
      end if
      lblend = blend
    else
      lblend = 0
    end if

!   choose element

    if ( lblend == 0 ) then
      do curve = 1, mesh%ncurves
        elementc(curve) = mesh%curves(curve)%element
      end do
    else
      do curve = 1, mesh%ncurves
        elementc(curve) = mesh%curves(curve)%element_blend(lblend)
      end do
    end if

!   fill numnodtop (works also for meshparts=.false.)

    do curve = 1, mesh%ncurves
      numnodtop(curve,1) = 0
      numnodtop(curve,2) = mesh%curves(curve)%element%numnod
      do m = 1, mesh%curves(curve)%nblend
        numnodtop(curve,m+2) = numnodtop(curve,m+1) &
                               + mesh%curves(curve)%element_blend(m)%numnod
      end do
    end do

!   fit plot into figsize x figsize

    xmin = max ( minval ( mesh%coor(:,1) ) + plot_options%xmin_coor_add, &
                 plot_options%xmin )
    xmax = min ( maxval ( mesh%coor(:,1) ) + plot_options%xmax_coor_add, &
                 plot_options%xmax )
    ymin = max ( minval ( mesh%coor(:,2) ) + plot_options%ymin_coor_add, &
                 plot_options%ymin )
    ymax = min ( maxval ( mesh%coor(:,2) ) + plot_options%ymax_coor_add, &
                 plot_options%ymax )

    scalefac = plot_options%figsize / max ( xmax - xmin , ymax - ymin )

!   plot points

    if ( plot_options%printlabels ) then

      do point = 1, mesh%npoints

        if ( mesh%points(point) == 0 ) cycle  ! sepran isolated points

        if ( mesh%points(point) <= mesh%nnodes_blend(lblend+1) .or. &
             mesh%points(point)  > mesh%nnodes_blend(lblend+2) ) cycle

!       plot point

        x(1,1) = mesh%coor(mesh%points(point),1)
        x(2,1) = mesh%coor(mesh%points(point),2)

        if ( any ( x(:,1) < [xmin,ymin] ) .or. &
             any ( x(:,1) > [xmax,ymax] ) )  cycle

        x(1,1) = scalefac * ( x(1,1) - xmin )
        x(2,1) = scalefac * ( ymax - x(2,1) )

        write ( ch, '(a,i0)' ) 'P', point

        call fig_text ( plot_options%fontsize, x(1,1), x(2,1), ch )

      end do

    end if

!   plot curves

    do curve = 1, mesh%ncurves

      if ( mesh%curves(curve)%nnodes == 0 ) cycle  ! skip empty curve (sepran?)

      npoints = elementc(curve)%numnod

      if ( npoints > MAXPOLYLINEPOINTS ) then
        write(*,'(/a/a,i0/)') &
        ' Error in plot_points_curves_2D: increase MAXPOLYLINEPOINTS ', &
        '  to at least ', npoints
        stop
      end if

      do elem = 1, mesh%curves(curve)%nelem

!       plot element on curve

        x(1,1:npoints) = &
            mesh%coor(mesh%curves(curve)%topology([(i,i=1,npoints)]+ &
                                          numnodtop(curve,lblend+1),elem,2),1)
        x(2,1:npoints) = &
            mesh%coor(mesh%curves(curve)%topology([(i,i=1,npoints)]+ &
                                          numnodtop(curve,lblend+1),elem,2),2)

        if ( any( x(1,:npoints) < xmin ) .or. any( x(1,:npoints) > xmax ) .or. &
           any( x(2,:npoints) < ymin ) .or. any( x(2,:npoints) > ymax ) ) cycle

        x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
        x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

        if ( plot_options%spline ) then
          call fig_spline ( 1, npoints, x(:,1:npoints) )
        else
          call fig_polyline ( 1, npoints, x(:,1:npoints) )
        end if

      end do

      if ( plot_options%printlabels ) then

        write ( ch, '(a,i0)' ) 'C', curve

        xminc = mesh%coor(mesh%curves(curve)%topology(&
                  1+numnodtop(curve,lblend+1),1,2),1)
        xmaxc = mesh%coor(mesh%curves(curve)%topology(&
                 elementc(curve)%numnod+numnodtop(curve,lblend+1),&
                   mesh%curves(curve)%nelem,2),1)
        yminc = mesh%coor(mesh%curves(curve)%topology(&
                  1+numnodtop(curve,lblend+1),1,2),2)
        ymaxc = mesh%coor(mesh%curves(curve)%topology(&
                 elementc(curve)%numnod+numnodtop(curve,lblend+1),&
                   mesh%curves(curve)%nelem,2),2)

        xc = 0.65_dp*xminc + 0.35_dp*xmaxc
        yc = 0.65_dp*yminc + 0.35_dp*ymaxc

        xc = scalefac * ( xc - xmin )
        yc = scalefac * ( ymax - yc )

        call fig_text ( plot_options%fontsize, xc, yc, ch )

      end if

    end do

!   plot surface labels

    if ( plot_options%printlabels ) then

      do surface = 1, mesh%nsurfaces

        xmins = minval ( mesh%coor(mesh%surfaces(surface)%nodes,1) )
        xmaxs = maxval ( mesh%coor(mesh%surfaces(surface)%nodes,1) )
        ymins = minval ( mesh%coor(mesh%surfaces(surface)%nodes,2) )
        ymaxs = maxval ( mesh%coor(mesh%surfaces(surface)%nodes,2) )

        xc = ( xmaxs + xmins ) / 2
        yc = ( ymaxs + ymins ) / 2

        if ( any ( [xc,yc] < [xmin,ymin] ) .or. &
             any ( [xc,yc] > [xmax,ymax] ) )  cycle

        xc = scalefac * ( xc - xmin )
        yc = scalefac * ( ymax - yc )

        write ( ch, '(a,i0)' ) 'S', surface

        call fig_text ( plot_options%fontsize, xc, yc, ch )

      end do

    end if

    close(unit=unit_figplot)

  end subroutine plot_points_curves_2D


! plot objects

  subroutine plot_objects ( plot_options, mesh, filename, append, object1, &
    object2 )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   when .true.: append plot to existing file in order to overlay plots.
!   default=.false.
    logical, intent(in), optional :: append

!   if these are present only objects object1,...,object2 are plotted.
!   If only object1 is present one object is plotted.
    integer, intent(in), optional :: object1, object2


    integer :: object
    real(dp), dimension(3) :: xunitvec, yunitvec
    type(mesh_t) :: lmesh


    call check ( mesh, 'plot_objects' )
    call check_blend ( mesh, 'plot_objects' )

!   test dimension of mesh

    if ( mesh%ndim == 3 ) then

!     3D mesh: project coordinates onto a 2D plane.

!     compute projection vectors

      call projection_vectors ( &
        plot_options%viewpoint, plot_options%ydirection, xunitvec, yunitvec )

!     copy mesh to local mesh

      call copy ( mesh, lmesh )
      call fill_mesh_parts ( lmesh )

!     project coordinates in new mesh

      call project_coordinates ( xunitvec, yunitvec, lmesh%coor )

      do object = 1, mesh%nobjects
        call project_coordinates ( xunitvec, yunitvec, &
          lmesh%objects(object)%coor )
      end do

      call plot_objects_2D ( plot_options, lmesh, filename, append=append, &
        object1=object1, object2=object2 )

      call delete ( lmesh )

    else

      call plot_objects_2D ( plot_options, mesh, filename, append=append, &
        object1=object1, object2=object2 )

    end if

  end subroutine plot_objects


! plot objects (2D)

  subroutine plot_objects_2D ( plot_options, mesh, filename, append, object1, &
    object2 )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   if .true.: append plot to existing file in order to overlay plots.
!   default=.false.
    logical, intent(in), optional :: append

!   if these are present only objects object1,...,object2 are plotted.
!   If only object1 is present one object is plotted.
    integer, intent(in), optional :: object1, object2


    logical :: fileexists, fileappend
    integer :: npoints
    integer :: node, object
    integer :: obj1, obj2
    real(dp) :: x(2,MAXPOLYLINEPOINTS)
    real(dp) :: scalefac, xmin, xmax, ymin, ymax


    inquire ( file=filename, exist=fileexists )

    if ( present(append) ) then
      fileappend = append
    else
      fileappend = .false.
    end if

!   which objects?

    if ( present(object1) .and. present(object2) ) then
!     specified range of objects only
      obj1 = object1
      obj2 = object2
    else if ( present(object1) ) then
!     one object only
      obj1 = object1
      obj2 = object1
    else
!     all objects
      obj1 = 1
      obj2 = mesh%nobjects
    end if

    if ( obj1 < 1 .or. obj1 > mesh%nobjects .or.  &
         obj2 < 1 .or. obj2 > mesh%nobjects ) then

      write(*,'(/2a/2(a,i0),/a,i0/)') &
        'Error: object1 and/or object2 in the heading of ', &
        'plot_objects_2D is ', &
        'out of range. object1 is ', obj1, '; object2 is ', obj2, &
        'whereas the number of objects is ', mesh%nobjects
      stop

    end if

!   open file

    if ( fileexists .and. fileappend ) then

      open ( unit=unit_figplot, file=filename, status='old', position='append' )

    else

      open ( unit=unit_figplot, file=filename )
      call fig_header
      call print_head ( plot_options, filename )

    end if

!   fit plot into figsize x figsize

    xmin = max ( minval ( mesh%coor(:,1) ) + plot_options%xmin_coor_add, &
                 plot_options%xmin )
    xmax = min ( maxval ( mesh%coor(:,1) ) + plot_options%xmax_coor_add, &
                 plot_options%xmax )
    ymin = max ( minval ( mesh%coor(:,2) ) + plot_options%ymin_coor_add, &
                 plot_options%ymin )
    ymax = min ( maxval ( mesh%coor(:,2) ) + plot_options%ymax_coor_add, &
                 plot_options%ymax )

    scalefac = plot_options%figsize / max ( xmax - xmin , ymax - ymin )

!   boundary

    if ( plot_options%plotboundary ) then

      call plot_boundary_internal ( plot_options, mesh, scalefac, &
        xmin, xmax, ymin, ymax )

    end if

!   plot coordinates of object

    npoints = 2

    do object = obj1, obj2

      do node = 1, mesh%objects(object)%nnodes

        x(:,1) = mesh%objects(object)%coor(node,1:2)

        if ( any ( x(:,1) < [xmin,ymin] ) .or. &
             any ( x(:,1) > [xmax,ymax] ) )  cycle

        x(1,1) = scalefac * ( x(1,1) - xmin )
        x(2,1) = scalefac * ( ymax - x(2,1) )

        x(:,2) = x(:,1) + [plot_options%objectpointsize,0._dp]

        call fig_filledcircle ( pen_color = plot_options%objectpointcolor, &
          x=x(:,1:npoints) )

      end do

    end do

    close(unit=unit_figplot)

  end subroutine plot_objects_2D


! plot mesh

  subroutine plot_mesh ( plot_options, mesh, filename, append, surfaces, &
    groups )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   if .true.: append plot to existing file in order to overlay plots.
!   default=.false.
    logical, intent(in), optional :: append

!   if present the mesh will be plotted on the given surfaces only
    integer, dimension(:), intent(in), optional :: surfaces

!   if present: the element groups to be plotted
!   For example groups=(/2,4/) will plot mesh for element groups 2 and 4.
!   Default: all groups are plotted
    integer, dimension(:), intent(in), optional :: groups

    real(dp), dimension(3) :: xunitvec, yunitvec
    type(mesh_t) :: lmesh


    call check ( mesh, 'plot_mesh' )
    call check_blend ( mesh, 'plot_mesh' )

    if ( present(surfaces) ) then
      if ( any(surfaces < 1) .or. any(surfaces > mesh%nsurfaces) ) then
        write(*,'(/2a/(a,i0)/a,20(i0,1x)/)') &
          'Error: surfaces in the heading of plot_mesh are ', &
          'out of range.', &
          'the number of surfaces is ', mesh%nsurfaces, &
          'whereas surfaces are ', surfaces
        stop
      end if
    end if

!   test dimension of mesh

    if ( mesh%ndim == 3 ) then

!     3D mesh: project coordinates onto a 2D plane.

!     compute projection vectors

      call projection_vectors ( &
        plot_options%viewpoint, plot_options%ydirection, xunitvec, yunitvec )

      if ( present(surfaces) ) then

        if ( present(groups) ) then
          write(*,'(/a/a,i0/)') &
            ' Error in plot_mesh: the heading parameter groups cannot be ', &
            ' used together with the heading parameter surfaces.'
          stop
        end if

!       create new mesh from surface topology

        call generate_surface_mesh ( mesh, lmesh, surfaces )

      else

        write( *, '(/5(a/))') &
          'Warning: plotting a full mesh in 3D may lead to ', &
          ' unexpected results, such as: ', &
          ' - wrong color for boundary edges',  &
          ' - edges to center nodes ', &
          ' use surfaces= argument to plot on surfaces only '

!       copy mesh to local mesh

        call copy ( mesh, lmesh )
        call fill_mesh_parts ( lmesh )

      end if

!     project coordinates in new mesh

!     allocate new data to coor array and copy mesh%coor

      call project_coordinates ( xunitvec, yunitvec, lmesh%coor )

      call plot_mesh_2D ( plot_options, lmesh, filename, append=append, &
        groups=groups )

      call delete ( lmesh )

    else

      call plot_mesh_2D ( plot_options, mesh, filename, append=append, &
        groups=groups )

    end if

  end subroutine plot_mesh


! plot mesh (2D)

  subroutine plot_mesh_2D ( plot_options, mesh, filename, append, groups )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   if .true.: append plot to existing file in order to overlay plots.
!   default=.false.
    logical, intent(in), optional :: append

!   if present: the element groups to be plotted
!   For example groups=(/2,4/) will plot mesh for element groups 2 and 4.
!   Default: all groups are plotted
    integer, dimension(:), intent(in), optional :: groups


    logical :: fileexists, fileappend
    integer :: i
    integer :: npoints, grp, elgrp, elem, side, nodes(MAXPOLYLINEPOINTS)
    integer, allocatable, dimension(:) :: lgroups
    real(dp) :: x(2,MAXPOLYLINEPOINTS)
    real(dp) :: scalefac, xmin, xmax, ymin, ymax


!   check groups

    if ( present(groups) ) then
      if ( any ( groups < 1 ) .or. any ( groups > mesh%nelgrp ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: parameter groups in the heading of plot_mesh_2D is ', &
          'out of range: some groups are < 1 or larger than ', mesh%nelgrp
        stop
      end if
      allocate(lgroups(size(groups)))
      lgroups = groups
    else
      allocate(lgroups(mesh%nelgrp))
      lgroups = [ (i,i=1,mesh%nelgrp) ]
    end if

    inquire ( file=filename, exist=fileexists )

    if ( present(append) ) then
      fileappend = append
    else
      fileappend = .false.
    end if

    if ( fileexists .and. fileappend ) then

      open ( unit=unit_figplot, file=filename, status='old', position='append' )

    else

      open ( unit=unit_figplot, file=filename )
      call fig_header
      call print_head ( plot_options, filename )

    end if

!   fit plot into figsize x figsize

    xmin = max ( minval ( mesh%coor(:,1) ) + plot_options%xmin_coor_add, &
                 plot_options%xmin )
    xmax = min ( maxval ( mesh%coor(:,1) ) + plot_options%xmax_coor_add, &
                 plot_options%xmax )
    ymin = max ( minval ( mesh%coor(:,2) ) + plot_options%ymin_coor_add, &
                 plot_options%ymin )
    ymax = min ( maxval ( mesh%coor(:,2) ) + plot_options%ymax_coor_add, &
                 plot_options%ymax )

    scalefac = plot_options%figsize / max ( xmax - xmin , ymax - ymin )

!   plot boundary if not all groups are present

    if ( plot_options%plotboundary .and. present(groups) ) then

      call plot_boundary_internal ( plot_options, mesh, scalefac, xmin, xmax, &
        ymin, ymax )

    end if

!   plot sides of elements

    do grp = 1, size(lgroups)

      elgrp = lgroups(grp)

!     choose between lines or plane/volume elements

      if ( any ( mesh%element(elgrp)%elshape == [ 1, 2, 32, 101 ] ) ) then

!       line elements

        npoints= mesh%element(elgrp)%numnod

        do elem = 1, mesh%grpnumel(elgrp)

          nodes(1:npoints) = mesh%topology(elgrp)%a(1:npoints,elem)

          x(1,1:npoints) = mesh%coor(nodes(1:npoints),1)
          x(2,1:npoints) = mesh%coor(nodes(1:npoints),2)

          if ( any( x(1,:npoints) < xmin ) .or. &
               any( x(1,:npoints) > xmax ) .or. &
               any( x(2,:npoints) < ymin ) .or. &
               any( x(2,:npoints) > ymax ) ) cycle

          x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
          x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

          if ( plot_options%spline ) then
            call fig_spline ( plot_options%meshcolor, npoints, &
              x(:,1:npoints) )
          else
            call fig_polyline ( plot_options%meshcolor, npoints, &
              x(:,1:npoints) )
          end if

        end do

      else

!       plane or volume elements

        npoints= mesh%element(elgrp)%sidnumnod

        if ( npoints == 0 ) then
          write(*,'(/a,i0/)') &
            'Error plot_mesh: element shape cannot be plotted = ', &
             mesh%element(elgrp)%elshape
          stop
        end if

        do elem = 1, mesh%grpnumel(elgrp)
          do side = 1, mesh%element(elgrp)%numsides

            nodes(1:npoints) = mesh%topology(elgrp)%a(&
                   &mesh%element(elgrp)%sidnod(1:npoints,side),elem)

!           plot sides

            if ( ( nodes(npoints) > nodes(1) .and. &
                   elgrp == mesh%sidelem(elgrp)%a(side,elem,1) ) .or. &
                   elgrp /= mesh%sidelem(elgrp)%a(side,elem,1) .or. &
                   mesh%sidelem(elgrp)%a(side,elem,4) <= 0 ) then

!             plot sides if:
!              - between elements inside a single group (plot only once)
!              - between elements from a different group (plot twice to
!                avoid plotting none in 3D surfaces.
!              - sides are near a boundary.


              x(1,1:npoints) = mesh%coor(nodes(1:npoints),1)
              x(2,1:npoints) = mesh%coor(nodes(1:npoints),2)

              if ( any( x(1,:npoints) < xmin ) .or. &
                   any( x(1,:npoints) > xmax ) .or. &
                   any( x(2,:npoints) < ymin ) .or. &
                   any( x(2,:npoints) > ymax ) ) cycle

              x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
              x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

              if ( elgrp /= mesh%sidelem(elgrp)%a(side,elem,1) .or. &
                   mesh%sidelem(elgrp)%a(side,elem,4) <= 0 ) then
!               interface between element groups or a boundary
                if ( plot_options%spline ) then
                  call fig_spline ( plot_options%boundarycolor, npoints, &
                    x(:,1:npoints) )
                else
                  call fig_polyline ( plot_options%boundarycolor, npoints, &
                    x(:,1:npoints) )
                end if
              else
!               internal elements
                if ( plot_options%spline ) then
                  call fig_spline ( plot_options%meshcolor, npoints, &
                    x(:,1:npoints) )
                else
                  call fig_polyline ( plot_options%meshcolor, npoints, &
                    x(:,1:npoints) )
                end if
              end if

            end if

          end do
        end do

      end if

    end do

    deallocate(lgroups)

    close(unit=unit_figplot)

  end subroutine plot_mesh_2D


! color plot

  subroutine plot_color_fill ( plot_options, mesh, problem, filename, physq, &
    layer, degfd, surfaces, groups, sysvector, vector )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   the physical quantity to be plotted
    integer, intent(in), optional :: physq

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered. Only elements that have all nodes within the layer
!   will be used.
    integer, intent(in), optional :: layer

!   the degree of freedom to be plotted (within the physical quantity when
!   physq is also present)
    integer, intent(in), optional :: degfd

!   when present the data will be plotted on the given surfaces only
    integer, dimension(:), intent(in), optional :: surfaces

!   if present: the element groups to be plotted
!   For example groups=(/2,4/), will plot for element groups 2 and 4 only
!   Default: all groups are plotted
    integer, dimension(:), intent(in), optional :: groups

!   sysvector to be plotted
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be plotted
    type(vector_t), intent(in), optional :: vector


!   This routine plots colors between contours. It is based on the filling of
!   triangles where the value in each of the three vertices is known. The
!   quantity is interpolated linearly within the triangle.
!   See subroutine plot_color_2D for more information.


    real(dp), dimension(3) :: xunitvec, yunitvec
    type(mesh_t) :: lmesh


    call check ( mesh, 'plot_color_fill' )
    call check_blend ( mesh, 'plot_color_fill' )
    call check ( problem, 'plot_color_fill', mesh )

    if ( present(surfaces) ) then
      if ( any(surfaces < 1) .or. any(surfaces > mesh%nsurfaces) ) then
        write(*,'(/2a/(a,i0)/a,20i0/)') &
          'Error: surfaces in the heading of plot_color_fill are ', &
          'out of range.', &
          'the number of surfaces is ', mesh%nsurfaces, &
          'whereas surfaces are ', surfaces
        stop
      end if
    end if

!   test dimension of mesh

    if ( mesh%ndim == 3 ) then

!     3D mesh: project coordinates onto a 2D plane.

      if ( present(vector) ) then
        if ( vector%elementwise ) then
          write( *, '(/2(a/))') &
            'Error plot_color_fill: ', &
            ' in 3D vectors stored elementwise cannot be plotted'
          stop
        end if
      end if

!     compute projection vectors

      call projection_vectors ( &
        plot_options%viewpoint, plot_options%ydirection, xunitvec, yunitvec )

      if ( present(surfaces) ) then

!       create new mesh from surface topology

        if ( present(groups) ) then
          write(*,'(/a/a,i0/)') &
            ' Error in plot_color_fill: the heading parameter groups cannot ', &
            ' be used together with the heading parameter surfaces.'
          stop
        end if

        if ( problem%numinactivegroups > 0 ) then
          write(*,'(/a/a,i0/)') &
            ' Warning in  plot_color_fill: inactive groups cannot be ', &
            ' used together with the heading parameter surfaces.'
        end if

        call generate_surface_mesh ( mesh, lmesh, surfaces )

      else

        write( *, '(/3(a/))') &
          'Error plot_color_fill: ', &
          ' in 3D only plotting on surfaces is supported', &
          ' use surface=argument to plot on surface only '
        stop

      end if

!     project coordinates in new mesh

      call project_coordinates ( xunitvec, yunitvec, lmesh%coor )

      call plot_color_fill_2D ( plot_options, lmesh, problem, filename, &
        physq, layer, degfd, sysvector=sysvector, vector=vector, &
        projection=.true. )

      call delete ( lmesh )

    else

!     2D mesh

      call plot_color_fill_2D ( plot_options, mesh, problem, filename, &
        physq, layer, degfd, groups=groups, sysvector=sysvector,  &
        vector=vector )

    end if

  end subroutine plot_color_fill


! color plot

  subroutine plot_color_fill_2D ( plot_options, mesh, problem, filename, &
    physq, layer, degfd, groups, sysvector, vector, projection )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   the physical quantity to be plotted
    integer, intent(in), optional :: physq

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered. Only elements that have all nodes within the layer
!   will be used.
    integer, intent(in), optional :: layer

!   the degree of freedom to be plotted (within the physical quantity if
!   physq is also present)
!   default=1
    integer, intent(in), optional :: degfd

!   the element groups to be plotted
!   default: all groups
    integer, dimension(:), intent(in), optional :: groups

!   sysvector to be plotted
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be plotted
    type(vector_t), intent(in), optional :: vector

!   logical indicating whether routine is used for 3D projection (.true.) or
!   not (.false.). The default is .false.
    logical, intent(in), optional :: projection


!   This routine plots colors between contours. It is based on the filling of
!   triangles where the value in each of the three vertices is known. The
!   quantity is interpolated linearly within the triangle.
!   The filling is based on the following model:
!
!                        med
!                      /    \                      !
!                    /      | \                    !
!                  / |      |   \             min <= med <= max
!                /   |      |   | \                !
!              / |   |      |   |   \              !
!            /   |   |      |   |   | \            !
!          / |   |   |      |   |   |  |\          !
!        /   |   |   |      |   |   |  |  \        !
!      /     |   |   |      |   |   |  |    \      !
!    min --------------------------------- max
!
!   having
!      - two filled triangles near the min and max value
!      - a number of filled quadrilaterals (can be zero)
!      - one filled five points polygon near the med value
!
!   there exist degenerate cases when there are no level lines between min and
!   med, med and max or both. In the latter case there is just one color in the
!   triangle. In the two other cases the five points polygon degenerates into
!   an additional quadrilateral, for example:
!
!                        med
!                      /    \                      !
!                    /        \                    !
!                  / |          \                  !
!                /   |            \                !
!              /     |              \              !
!            /       |                \            !
!          / |       |                  \          !
!        /   |       |                    \        !
!      /     |       |                      \      !
!    min --------------------------------- max
!


    logical :: onecolor, full, first, proj
    logical, dimension(mesh%nnodes) :: wkl
    integer :: elgrp, elem, side, nodes(MAXPOLYLINEPOINTS)
    integer :: deg, numnod, nodenr, npoints
    integer :: pos(problem%maxnoddegfd), dof, ncolors, bp, subelem
    integer :: posv(problem%maxvecnoddegfd)
    integer :: nodenrs(3), nintervals
    integer :: intervalmax, intervalmed, intervalmin
    integer :: nodmax, nodmed, nodmin, numquad, quad, corr
    integer :: interval, colbarstart, colbarend, level, ntext, fontsize
    integer :: ndofu, elnumnod, elnodes(mesh%maxelnumnod), grp, i
    integer, dimension(:), allocatable :: lgroups

    real(dp) :: x(2,MAXPOLYLINEPOINTS), xvert(2,3), values(3)
    real(dp) :: elcoor(2,mesh%maxelnumnod)
    real(dp) :: scalefac, xmin, xmax, ymin, ymax
    real(dp) :: umax, umin, delta, umx, umn
    real(dp) :: valmax, valmed, valmin, f, xc(2,6), xtmp(5)
    real(dp) :: offsetx, offsety, widthbar, heightcbox, xs, heightbar
    real(dp) :: shiftmincbox, heightmaxmincbox, deltatxt, heighttxtbox
    real(dp) :: shiftmaxcbox, dx, shiftrotate, xt, yt, chunit
    real(dp), dimension(:), allocatable :: u
    character(len=20) :: ch
    character(len=8) :: frmt


    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in plot_color_fill_2D: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(sysvector) ) then
      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in plot_color_fill_2D: sysvector has not been created '
        stop
      end if
      allocate(u(problem%maxnoddegfd*mesh%maxelnumnod))
    else if ( present(vector) ) then
      if ( .not. vector%created ) then
        write(*,'(/a/)') &
          'Error in plot_color_fill_2D: vector has not been created '
        stop
      end if
      allocate(u(problem%maxvecnoddegfd*mesh%maxelnumnod))
    end if

    if ( present(projection) ) then
      proj = projection
    else
      proj = .false.
    end if

    if ( present(degfd) ) then
      deg = degfd
    else
      deg = 1
    end if

    if ( plot_options%usefloating ) then
      frmt = '(f9.2)'
    else
      frmt = '(es9.2)'
    end if

!   check groups
    if ( present(groups) ) then
      if ( any ( groups < 1 ) .or. any ( groups > mesh%nelgrp ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: parameter groups in the heading of plot_color_fill_2D is ', &
          'out of range: some groups are < 1 or larger than ', mesh%nelgrp
        stop
      end if
      allocate(lgroups(size(groups)))
      lgroups = groups
    else
      allocate(lgroups(mesh%nelgrp))
      lgroups = [ (i,i=1,mesh%nelgrp) ]
    end if

!   check element

    if ( any ( mesh%element(lgroups)%sidnumnod == 0 ) ) then
      write(*,'(/2a,i0/)') &
        'Error plot_color_fill: some of these element shapes ', &
        'cannot be plotted = ', mesh%element(:)%elshape
      stop
    end if

!   layer

    if ( present(layer) ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in plot_color_fill: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a),i0/)') &
          'Error in plot_color_fill: ',&
          ' layer < 1 or layer > number of layers = ', &
          problem%numlayers
        stop
      end if
      wkl = btest(problem%nodlayers,layer-1)
    else if ( problem%numlayers > 0 ) then
      write(*,'(3(/a)/)') &
        'Error in plot_color_fill: ', &
        ' layers have been defined and the layer keyword is not present', &
        ' layers not implemented for this routine without the layer keyword'
      stop
    end if

    open(unit=unit_figplot,file=filename)

!   write header

    call fig_header

!   write colors

    ncolors = plot_options%numcolors

    call fig_colors ( ncolors )

    nintervals = ncolors

!   print filename in title

    call print_head ( plot_options, filename )

!   fit plot into figsize x figsize

    xmin = max ( minval ( mesh%coor(:,1) ) + plot_options%xmin_coor_add, &
                 plot_options%xmin )
    xmax = min ( maxval ( mesh%coor(:,1) ) + plot_options%xmax_coor_add, &
                 plot_options%xmax )
    ymin = max ( minval ( mesh%coor(:,2) ) + plot_options%ymin_coor_add, &
                 plot_options%ymin )
    ymax = min ( maxval ( mesh%coor(:,2) ) + plot_options%ymax_coor_add, &
                 plot_options%ymax )

    scalefac = plot_options%figsize / max ( xmax - xmin , ymax - ymin )

!   plot boundary if not all groups are present

    if ( plot_options%plotboundary .and. ( present(groups) .or. &
           ( problem%numinactivegroups > 0 .and. .not. proj ) ) ) then

      call plot_boundary_internal ( plot_options, mesh, scalefac, xmin, xmax, &
        ymin, ymax )

    end if

!   find maximum and minimum value

    if ( present(sysvector) ) then

!     sysvector

      umax = 0
      umin = 0

      first = .true.

      do nodenr = 1, mesh%nnodes

        if ( .not. node_in_groups ( mesh, lgroups, nodenr ) ) cycle
        if ( .not. node_in_groups ( mesh, problem%activegroups, nodenr ) &
               .and. .not. proj ) cycle
        if ( present(layer) ) then
          if ( .not. wkl(nodenr) ) cycle
        end if
        if ( present(physq) ) then
          call pos_array_node ( problem, nodenr, dof, pos, [physq], layer )
        else
          call pos_array_node ( problem, nodenr, dof, pos, layer=layer )
        end if
        if ( first ) then
          first = .false.
          umax=sysvector%u(pos(deg))
          umin=sysvector%u(pos(deg))
        else
          umax=max(umax,sysvector%u(pos(deg)))
          umin=min(umin,sysvector%u(pos(deg)))
        end if

      end do

    else if ( present(vector) ) then

      if ( vector%elementwise ) then

!       vector elementwise

        umax = 0
        umin = 0

        first = .true.

        do elgrp = 1, mesh%nelgrp

          if ( .not. any ( lgroups == elgrp ) ) cycle

!         skip inactive groups
          if ( any( elgrp == problem%inactivegroups ) .and. .not. proj ) cycle

          numnod = mesh%element(elgrp)%numnod

          do elem = 1, mesh%grpnumel(elgrp)

!           get element vector
            call get_vector ( mesh, problem, vector, elgrp, elem, u, &
              sloppy=.true. )
            bp = ( deg - 1 ) * numnod
            if ( first ) then
              first = .false.
              umax=maxval(u(bp+1:bp+numnod))
              umin=minval(u(bp+1:bp+numnod))
            else
              umax=max(umax,maxval(u(bp+1:bp+numnod)))
              umin=min(umin,minval(u(bp+1:bp+numnod)))
            end if

          end do

        end do

      else

!       vector in nodes

        umax = 0
        umin = 0

        first = .true.

        do nodenr = 1, mesh%nnodes
          if ( .not. node_in_groups ( mesh, lgroups, nodenr ) ) cycle
          if ( .not. node_in_groups ( mesh, problem%activegroups, nodenr ) &
                 .and. .not. proj ) cycle
          if ( present(layer) ) then
            if ( .not. wkl(nodenr) ) cycle
          end if
          call pos_array_vec_node ( problem, nodenr, dof, posv, vector%vec, &
            layer )
          if ( first ) then
            first = .false.
            umax=vector%u(posv(deg))
            umin=vector%u(posv(deg))
          else
            umax=max(umax,vector%u(posv(deg)))
            umin=min(umin,vector%u(posv(deg)))
          end if

        end do

      end if

    end if

!   determine interval between levels

    if ( ( umax > plot_options%maxvalue .and. umin < plot_options%minvalue ) &
            .or. plot_options%forcemaxmin ) then
!     one extra interval at the top and bottom
      umx = plot_options%maxvalue
      umn = plot_options%minvalue
      delta = ( umx - umn ) / ( nintervals - 2 )
      corr = 1   ! shift 1 to make all 1 <= interval <= nintervals
      colbarstart = 2
      colbarend = nintervals - 1
    else if ( umax > plot_options%maxvalue ) then
!     one extra interval at the top
      umx = plot_options%maxvalue
      umn = umin
      delta = ( umx - umn ) / ( nintervals - 1 )
      corr = 0
      colbarstart = 1
      colbarend = nintervals - 1
    else if ( umin < plot_options%minvalue ) then
!     one extra interval at the bottom
      umx = umax
      umn = plot_options%minvalue
      delta = ( umx - umn ) / ( nintervals - 1 )
      corr = 1   ! shift 1 to make all 1 <= interval <= nintervals
      colbarstart = 2
      colbarend = nintervals
    else
!     all intervals between max and min
      umx = umax
      umn = umin
      delta = ( umax - umin ) / nintervals
      corr = 0
      colbarstart = 1
      colbarend = nintervals
    end if

!   color filled values and boundary

    do grp = 1, size(lgroups)

      elgrp = lgroups(grp)

!     skip inactive groups
      if ( any( elgrp == problem%inactivegroups ) .and. .not. proj ) cycle

      do elem = 1, mesh%grpnumel(elgrp)

!       skip elements not in layer
        if ( present(layer) ) then
          if ( .not. all ( wkl(mesh%topology(elgrp)%a(:,elem)) ) ) cycle
        end if

!       check whether element in plot window: otherwise skip element

        elnumnod = mesh%element(elgrp)%numnod
        elnodes(1:elnumnod)  = mesh%topology(elgrp)%a(:,elem)

        elcoor(1,1:elnumnod) = mesh%coor(elnodes(1:elnumnod),1)
        elcoor(2,1:elnumnod) = mesh%coor(elnodes(1:elnumnod),2)

        if ( any( elcoor(1,:elnumnod) < xmin ) .or. &
             any( elcoor(1,:elnumnod) > xmax ) .or. &
             any( elcoor(2,:elnumnod) < ymin ) .or. &
             any( elcoor(2,:elnumnod) > ymax ) ) cycle

!       color plot per element

        if ( present(sysvector) ) then

!         get element vector
          if ( present(physq) ) then
            call get_sysvector ( mesh, problem, sysvector, elgrp, elem, u, &
              [physq], ndofu=ndofu, sloppy=.true., layer=layer )
          else
            call get_sysvector ( mesh, problem, sysvector, elgrp, elem, u, &
              ndofu=ndofu, sloppy=.true., layer=layer )
          end if

        else if ( present(vector) ) then

!         get element vector
          call get_vector ( mesh, problem, vector, elgrp, elem, u, &
            ndofu=ndofu, sloppy=.true., layer=layer )

        end if

        numnod = mesh%element(elgrp)%numnod

        if ( ndofu < deg * numnod ) then
          write(*,'(/a/a/)') &
            'Error in plot_color_fill_2D:', &
            ' not enough unknowns in nodes for plotting'
          stop
        else if ( elem == 1 .and. mod(ndofu,numnod) /= 0 ) then
          write(*,'(/6(a/),a,i0/)') &
            'Warning in plot_color_fill_2D:', &
            ' number of unknowns in the sysvector or vector ', &
            ' is not a multiple of the number of nodes.', &
            ' Plot_color_fill assumes that the vector/sysvector is defined ', &
            ' in all nodes having the same number of degrees of freedom.', &
            ' Know what you are doing, I warned you!', &
            ' Element group ', elgrp
        end if

        do subelem = 1, mesh%element(elgrp)%numsubdiv

!         process each subtriangle individually

!         get nodes and coordinates of vertices

          nodes(1:3) = mesh%element(elgrp)%subtopol(:,subelem)

          nodenrs = mesh%topology(elgrp)%a(nodes(1:3),elem)

          xvert(1,1:3) = mesh%coor(nodenrs,1)
          xvert(2,1:3) = mesh%coor(nodenrs,2)

!         values in the vertices

          values = u ( ( deg - 1 ) * numnod + nodes(1:3) )

          valmax = maxval(values)
          valmin = minval(values)

          if ( valmax - valmin <= 1.e-14_dp * (umax - umin) ) then

!           value is constant

            onecolor = .true.
            nodmin = 0
            nodmax = 0
            nodmed = 0
            valmed = valmin

          else

!           other cases

            onecolor = .false.
            nodmin = maxval(minloc(values))
            nodmax = maxval(maxloc(values))
            nodmed = 6 - nodmin - nodmax
            valmed = values(nodmed)

          end if

!         determine interval between contour levels for the vertices

          intervalmax = int ( ( valmax - umn ) / delta + 1 + corr )
            intervalmax = min ( intervalmax, nintervals )
            intervalmax = max ( intervalmax, 1 )
          intervalmin = int ( ( valmin - umn ) / delta + 1 + corr )
            intervalmin = min ( intervalmin, nintervals )
            intervalmin = max ( intervalmin, 1 )
          intervalmed = int ( ( valmed - umn ) / delta + 1 + corr )
            intervalmed = min ( intervalmed, nintervals )
            intervalmed = max ( intervalmed, 1 )

          if ( intervalmax < intervalmed .or. &
               intervalmax < intervalmin .or. &
               intervalmed < intervalmin ) then

            write(*,*) 'error intervalmax, intervalmed, intervalmin'
            write(*,*) intervalmax, intervalmed, intervalmin
            write(*,*) 'valmax, valmed, valmin'
            write(*,*) valmax, valmed, valmin
            write(*,*) 'umin, umax, umn, umx, delta'
            write(*,*) umin, umax, umn, umx, delta
          end if

          if ( intervalmax == intervalmin .or. onecolor ) then

!           one filled triangle

            npoints = 4

            x(:,1) = xvert(:,1)
            x(:,2) = xvert(:,2)
            x(:,3) = xvert(:,3)
            x(:,4) = xvert(:,1)

            x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
            x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

            call fig_filled_polygon ( fill_color=intervalmax+31, &
              x=x(:,1:npoints), npoints=npoints )

            full = .true.

          end if

          if ( intervalmed > intervalmin .and. .not. onecolor ) then

!           one filled triangle from min point:
!
!                        med
!                      .     .
!                     .       .
!                    .         .
!                   2           .
!                 /   \          .
!                /     \          .
!               min --- 1  .....   max
!

            npoints = 4

            f = ( ( intervalmin - corr ) * delta + umn - valmin ) &
                  / ( valmax - valmin )
            x(:,1) = (1-f)*xvert(:,nodmin)+f*xvert(:,nodmax)
            f = ( ( intervalmin - corr ) * delta + umn - valmin ) &
                  / ( valmed - valmin )
            x(:,2) = (1-f)*xvert(:,nodmin)+f*xvert(:,nodmed)
            x(:,3) = xvert(:,nodmin)
            x(:,4) = x(:,1)

            if ( intervalmax /= intervalmed ) then
!             fill coordinates for remaining polygon
              xc(:,1) = xvert(:,nodmed)
              xc(:,2) = x(:,2)
              xc(:,3) = x(:,1)
              xc(:,6) = xc(:,1)
            end if

            x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
            x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

            call fig_filled_polygon ( fill_color=intervalmin+31, &
              x=x(:,1:npoints), npoints=npoints, test=.true. )

!           quadrilaterals

            if ( intervalmax == intervalmed ) then
              numquad = intervalmed - intervalmin
              full = .true.
            else
              numquad = intervalmed - intervalmin - 1
              full = .false.
            end if

            do quad = 1, numquad

!             filled quadrilaterals from min
!
!                       med
!                      .    .
!                     3      .
!                    / \      .
!                   2   \      .
!                 .  \   \      .
!                .    \   \      .
!               min ... 1--4 ... max
!

              npoints = 5

              f = ( ( intervalmin - corr + quad - 1 ) * delta + umn - valmin ) &
                   / ( valmax - valmin )
              x(:,1) = (1-f)*xvert(:,nodmin)+f*xvert(:,nodmax)
              f = ( ( intervalmin - corr + quad - 1 ) * delta + umn - valmin ) &
                   / ( valmed - valmin )
              x(:,2) = (1-f)*xvert(:,nodmin)+f*xvert(:,nodmed)
              if ( quad == numquad .and. full ) then
!               last
                x(:,3) = xvert(:,nodmed)
                x(:,4) = xvert(:,nodmax)
              else
                f = ( ( intervalmin - corr + quad ) * delta + umn - valmin ) &
                     / ( valmed - valmin )
                x(:,3) = (1-f)*xvert(:,nodmin)+f*xvert(:,nodmed)
                f = ( ( intervalmin - corr + quad ) * delta + umn - valmin ) &
                     / ( valmax - valmin )
                x(:,4) = (1-f)*xvert(:,nodmin)+f*xvert(:,nodmax)
              end if
              x(:,5) = x(:,1)

              if ( quad == numquad .and. .not. full ) then
!               fill coordinates for remaining polygon
                xc(:,1) = xvert(:,nodmed)
                xc(:,2) = x(:,3)
                xc(:,3) = x(:,4)
                xc(:,6) = xc(:,1)
              end if

              x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
              x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

              call fig_filled_polygon ( fill_color=intervalmin+quad+31, &
                 x=x(:,1:npoints), npoints=npoints )

            end do

          end if

          if ( intervalmed < intervalmax .and. .not. onecolor ) then

!           one filled triangle from max point:
!
!                       med
!                     .     .
!                    .       .
!                   .         2
!                  .         / \        !
!                 .         /   \       !
!                .         /     \      !
!               min ..... 1 ---- max
!

            npoints = 4

            f = ( ( intervalmax-corr-1) * delta + umn - valmin ) &
                    / ( valmax - valmin )
            x(:,1) = (1-f)*xvert(:,nodmin)+f*xvert(:,nodmax)
            f = ( ( intervalmax-corr-1) * delta + umn - valmed ) &
                    / ( valmax - valmed )
            x(:,2) = (1-f)*xvert(:,nodmed)+f*xvert(:,nodmax)
            x(:,3) = xvert(:,nodmax)
            x(:,4) = x(:,1)

            if ( intervalmin /= intervalmed ) then
!             fill coordinates for remaining polygon
              xc(:,4) = x(:,1)
              xc(:,5) = x(:,2)
            end if

            x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
            x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

            call fig_filled_polygon ( fill_color=intervalmax+31, &
              x=x(:,1:npoints), npoints=npoints, test=.true. )

!           quadrilaterals

            if ( intervalmin == intervalmed ) then
              numquad = intervalmax - intervalmed
              full = .true.
            else
              numquad = intervalmax - intervalmed - 1
              full = .false.
            end if

            do quad = 1, numquad

!             filled quadrilaterals from max
!
!                       med
!                     .     .
!                    .       3
!                   .       / \       !
!                  .       /   2
!                 .       /   / .
!                .       /   /   .
!               min ... 4---1 ... max
!

              npoints = 5

              f = ( ( intervalmax - corr - quad ) * delta + umn - valmin ) &
                   / ( valmax - valmin )
              x(:,1) = (1-f)*xvert(:,nodmin)+f*xvert(:,nodmax)
              f = ( ( intervalmax - corr - quad ) * delta + umn - valmed ) &
                   / ( valmax - valmed )
              x(:,2) = (1-f)*xvert(:,nodmed)+f*xvert(:,nodmax)
              if ( quad == numquad .and. full ) then
!               last
                x(:,3) = xvert(:,nodmed)
                x(:,4) = xvert(:,nodmin)
              else
                f = ( ( intervalmax-corr - quad -1 ) * delta + umn - valmed ) &
                     / ( valmax - valmed )
                x(:,3) = (1-f)*xvert(:,nodmed)+f*xvert(:,nodmax)
                f = ( ( intervalmax-corr - quad -1 ) * delta + umn - valmin ) &
                     / ( valmax - valmin )
                x(:,4) = (1-f)*xvert(:,nodmin)+f*xvert(:,nodmax)
              end if
              x(:,5) = x(:,1)

              if ( quad == numquad .and. .not. full ) then
!               fill coordinates for remaining polygon
                xc(:,4) = x(:,4)
                xc(:,5) = x(:,3)
              end if

              x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
              x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

              call fig_filled_polygon ( fill_color=intervalmax-quad+31, &
                 x=x(:,1:npoints), npoints=npoints )

            end do

          end if

          if ( .not. full ) then

!           remaining filled polygon
!
!                       med
!                     /     \         !
!                    /       \        !
!                   2         \       !
!                  . \         5
!                 .   \       / .
!                .     \     /   .
!               min ... 3---4 ... max
!

            npoints = 6

            xc(1,1:npoints) = scalefac * ( xc(1,1:npoints) - xmin )
            xc(2,1:npoints) = scalefac * ( ymax - xc(2,1:npoints) )

            call fig_filled_polygon ( fill_color=intervalmed+31, &
               x=xc(:,1:npoints), npoints=npoints )

          end if

        end do

        if ( plot_options%plotboundary2 ) then

!         boundary

          npoints= mesh%element(elgrp)%sidnumnod

          do side = 1, mesh%element(elgrp)%numsides

!           plot side on boundary or between element groups (useful for surfaces
!           in 3D).

            if ( mesh%sidelem(elgrp)%a(side,elem,4) <= 0 .or. &
                 elgrp /= mesh%sidelem(elgrp)%a(side,elem,1) ) then

              nodes(1:npoints) = mesh%topology(elgrp)%a(&
                   &mesh%element(elgrp)%sidnod(1:npoints,side),elem)

              x(1,1:npoints) = mesh%coor(nodes(1:npoints),1)
              x(2,1:npoints) = mesh%coor(nodes(1:npoints),2)

              x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
              x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

              call fig_polyline ( 0, npoints, x(:,1:npoints) )

            end if

          end do

        end if

      end do
    end do

!   colorbar

    if ( plot_options%plotcolorbar ) then

!     plot all colors

      dx = plot_options%figsize / 150  ! scale factor (=1 for figsize=150mm)

      shiftmaxcbox     = 5   * dx  ! the maximum separate box shift upward
      shiftmincbox     = 10  * dx  ! the minimum separate box shift downward
      heightmaxmincbox = 5   * dx  ! height of the separate max of min box
      offsetx          = 10  * dx  ! initial offset of the colorbar in x direc
      offsety          = 40  * dx  ! initial offset of the colorbar in y direc
      widthbar         = 10  * dx  ! width of the colorbar
      heightbar        = 100 * dx  ! height of the colorbar
      shiftrotate      = offsetx + widthbar + 10 * dx  ! extra upward shift when
                                                       ! rotating the colorbar

!     height of a single colorbox
      heightcbox = heightbar / ( colbarend - colbarstart + 1 )
!     left x position of the colorbar
      xs = plot_options%figsize + offsetx

      npoints = 5

      if ( colbarstart == 2 ) then

!       add separate color block for small values u < umn

        x(:,1) = &
          [ xs, offsety - shiftmincbox ]
        x(:,2) = &
          [ xs + widthbar, offsety - shiftmincbox ]
        x(:,3) = x(:,2) + [0._dp, heightmaxmincbox]
        x(:,4) = x(:,1) + [0._dp, heightmaxmincbox]
        x(:,5) = x(:,1)

        x(2,1:npoints) = plot_options%figsize - x(2,1:npoints)

        if ( plot_options%rotatecolorbar ) then
          xtmp = x(1,1:npoints)
          x(1,1:npoints) = x(2,1:npoints)
          x(2,1:npoints) = xtmp - shiftrotate
        end if

        x(1,1:npoints) = x(1,1:npoints) + plot_options%shiftbarx
        x(2,1:npoints) = x(2,1:npoints) - plot_options%shiftbary

        call fig_filled_polygon ( fill_color=32, &
          x=x(:,1:npoints), npoints=npoints )

      end if

      do interval = colbarstart, colbarend

!       color block for values umn < u < umx

        x(:,1) = &
          [ xs, offsety + ( interval - colbarstart ) * heightcbox ]
        x(:,2) = &
          [ xs + widthbar, offsety + ( interval - colbarstart ) * heightcbox ]
        x(:,3) = x(:,2) + [0._dp, heightcbox]
        x(:,4) = x(:,1) + [0._dp, heightcbox]
        x(:,5) = x(:,1)

        x(2,1:npoints) = plot_options%figsize - x(2,1:npoints)

        if ( plot_options%rotatecolorbar ) then
          xtmp = x(1,1:npoints)
          x(1,1:npoints) = x(2,1:npoints)
          x(2,1:npoints) = xtmp - shiftrotate
        end if

        x(1,1:npoints) = x(1,1:npoints) + plot_options%shiftbarx
        x(2,1:npoints) = x(2,1:npoints) - plot_options%shiftbary

        call fig_filled_polygon ( fill_color=interval+31, &
          x=x(:,1:npoints), npoints=npoints )

      end do

      if ( colbarend == nintervals - 1 ) then

!       add separate color block large values u > umx

        x(:,1) = &
          [ xs, offsety + heightbar + shiftmaxcbox ]
        x(:,2) = &
          [ xs + widthbar, offsety + heightbar + shiftmaxcbox ]
        x(:,3) = x(:,2) + [0._dp, heightmaxmincbox]
        x(:,4) = x(:,1) + [0._dp, heightmaxmincbox]
        x(:,5) = x(:,1)

        x(2,1:npoints) = plot_options%figsize - x(2,1:npoints)

        if ( plot_options%rotatecolorbar ) then
          xtmp = x(1,1:npoints)
          x(1,1:npoints) = x(2,1:npoints)
          x(2,1:npoints) = xtmp - shiftrotate
        end if

        x(1,1:npoints) = x(1,1:npoints) + plot_options%shiftbarx
        x(2,1:npoints) = x(2,1:npoints) - plot_options%shiftbary

        call fig_filled_polygon ( fill_color=31+nintervals, &
          x=x(:,1:npoints), npoints=npoints )

      end if

!     text of numbers near the bar for umn <= u <= umx

      ntext = plot_options%numvalueintervals  ! number of intervals for text
      deltatxt = ( umx - umn ) / ntext      ! interval between values
      heighttxtbox = heightbar / ntext      ! distance between text
      fontsize = plot_options%fontsize      ! fontsize in points
      chunit = real(fontsize,dp)/MMFIGUNITS ! unit scaling with the fontsize

      do level = 1, ntext + 1

!       print level value

        write( ch, frmt ) umn + deltatxt * ( level - 1 )

        xt = xs + widthbar + 5*chunit
        yt = offsety + ( level - 1 ) * heighttxtbox &
                - 5*chunit
        yt = plot_options%figsize - yt
        if ( plot_options%rotatecolorbar ) then
          xtmp(1) = xt
          xt = yt - 20*chunit
          yt = xtmp(1) - shiftrotate + 14*chunit
        end if
        xt = xt + plot_options%shiftbarx
        yt = yt - plot_options%shiftbary

        call fig_text ( fontsize=fontsize, x=xt, y=yt, text=ch(1:9) )

      end do

!     print max/min values

      if ( plot_options%printmaxminvalues ) then

        write( ch, '(a,'//frmt(2:) ) 'max:', umax

        xt = xs
        yt = offsety - shiftmincbox - 20*chunit
        yt = plot_options%figsize - yt
        if ( plot_options%rotatecolorbar ) then
          xtmp(1) = xt
          xt = yt
          yt = xtmp(1) - shiftrotate + 10*chunit
        end if
        xt = xt + plot_options%shiftbarx
        yt = yt - plot_options%shiftbary

        call fig_text ( fontsize=fontsize, x=xt, y=yt, text=ch(1:13) )

        write( ch, '(a,'//frmt(2:) ) 'min:', umin

        yt = yt + 20*chunit

        call fig_text ( fontsize=fontsize, x=xt, y=yt, text=ch(1:13) )

      end if

    end if

!   deallocate memory

    if ( allocated(u) ) deallocate(u)
    deallocate(lgroups)

    close(unit=unit_figplot)

  end subroutine plot_color_fill_2D


! color plot

  subroutine plot_color_contour ( plot_options, mesh, problem, filename, &
    physq, layer, degfd, surfaces, groups, sysvector, vector )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   the physical quantity to be plotted
    integer, intent(in), optional :: physq

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered. Only elements that have all nodes within the layer
!   will be used.
    integer, intent(in), optional :: layer

!   the degree of freedom to be plotted (within the physical quantity when
!   physq is also present)
    integer, intent(in), optional :: degfd

!   when present the data will be plotted on the given surfaces only
    integer, dimension(:), intent(in), optional :: surfaces

!   if present: the element groups to be plotted
!   For example groups=(/2,4/), will plot for element groups 2 and 4 only
!   Default: all groups are plotted
    integer, dimension(:), intent(in), optional :: groups

!   sysvector to be plotted
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be plotted
    type(vector_t), intent(in), optional :: vector


!   This routine plots color contours. It is based on triangles where the value
!   in each of the three vertices is known. The quantity is interpolated
!   linearly within the triangle. The plotting is elementwise. This means that
!   for discontinuous quantities (vector%elemenwise=.true.) contour lines can
!   be discontinuous as well.


    real(dp), dimension(3) :: xunitvec, yunitvec
    type(mesh_t) :: lmesh


    call check ( mesh, 'plot_color_contour' )
    call check_blend ( mesh, 'plot_color_contour' )
    call check ( problem, 'plot_color_contour', mesh )

    if ( present(surfaces) ) then
      if ( any(surfaces < 1) .or. any(surfaces > mesh%nsurfaces) ) then
        write(*,'(/2a/(a,i0)/a,20i0/)') &
          'Error: surfaces in the heading of plot_color_contour are ', &
          'out of range.', &
          'the number of surfaces is ', mesh%nsurfaces, &
          'whereas surfaces are ', surfaces
        stop
      end if
    end if

!   test dimension of mesh

    if ( mesh%ndim == 3 ) then

!     3D mesh: project coordinates onto a 2D plane.

      if ( present(vector) ) then
        if ( vector%elementwise ) then
          write( *, '(/2(a/))') &
            'Error plot_color_contour: ', &
            ' in 3D vectors stored elementwise cannot be plotted'
          stop
        end if
      end if

!     compute projection vectors

      call projection_vectors ( &
        plot_options%viewpoint, plot_options%ydirection, xunitvec, yunitvec )

      if ( present(surfaces) ) then

!       create new mesh from surface topology

        if ( present(groups) ) then
          write(*,'(/a/a,i0/)') &
            ' Error in plot_color_contour: the heading parameter groups ', &
            ' cannot be used together with the heading parameter surfaces.'
          stop
        end if

        if ( problem%numinactivegroups > 0 ) then
          write(*,'(/a/a,i0/)') &
            ' Warning in  plot_color_contour: inactive groups cannot be ', &
            ' used together with the heading parameter surfaces.'
        end if

        call generate_surface_mesh ( mesh, lmesh, surfaces )

      else

        write( *, '(/3(a/))') &
          'Error plot_color_contour: ', &
          ' in 3D only plotting on surfaces is supported', &
          ' use surface= argument to plot on surface only '
        stop

      end if

!     project coordinates in new mesh

      call project_coordinates ( xunitvec, yunitvec, lmesh%coor )

      call plot_color_contour_2D ( plot_options, lmesh, problem, filename, &
        physq, layer, degfd, sysvector=sysvector, vector=vector, &
        projection=.true. )

      call delete ( lmesh )

    else

!     2D mesh

      call plot_color_contour_2D ( plot_options, mesh, problem, filename, &
        physq, layer, degfd, groups=groups, sysvector=sysvector,  &
        vector=vector )

    end if

  end subroutine plot_color_contour


! color plot

  subroutine plot_color_contour_2D ( plot_options, mesh, problem, filename, &
    physq, layer, degfd, groups, sysvector, vector, projection )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   the physical quantity to be plotted
    integer, intent(in), optional :: physq

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered. Only elements that have all nodes within the layer
!   will be used.
    integer, intent(in), optional :: layer

!   the degree of freedom to be plotted (within the physical quantity if
!   physq is also present)
!   default=1
    integer, intent(in), optional :: degfd

!   the element groups to be plotted
!   default: all groups
    integer, dimension(:), intent(in), optional :: groups

!   sysvector to be plotted
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be plotted
    type(vector_t), intent(in), optional :: vector

!   logical indicating whether routine is used for 3D projection (.true.) or
!   not (.false.). The default is .false.
    logical, intent(in), optional :: projection


!   This routine plots color contours. It is based on triangles where the value
!   in each of the three vertices is known. The quantity is interpolated
!   linearly within the triangle. The plotting is elementwise. This means that
!   for discontinuous quantities (vector%elemenwise=.true.) contour lines can
!   be discontinuous as well.

    logical :: constant, first, proj
    logical, dimension(mesh%nnodes) :: wkl
    integer :: elgrp, elem, side, nodes(MAXPOLYLINEPOINTS)
    integer :: deg, numnod, nodenr, npoints
    integer :: pos(problem%maxnoddegfd), dof, ncolors, bp, subelem
    integer :: posv(problem%maxvecnoddegfd)
    integer :: nodenrs(3), nintervals
    integer :: intervalmax, intervalmed, intervalmin
    integer :: nodmax, nodmed, nodmin, numcont, cont
    integer :: interval, fontsize, corr, ndofu, grp
    integer :: elnumnod, elnodes(mesh%maxelnumnod), i
    integer, dimension(:), allocatable :: lgroups

    real(dp) :: x(2,MAXPOLYLINEPOINTS), xvert(2,3), values(3)
    real(dp) :: elcoor(2,mesh%maxelnumnod)
    real(dp) :: scalefac, xmin, xmax, ymin, ymax
    real(dp) :: umax, umin, delta, umx, umn
    real(dp) :: valmax, valmed, valmin, f
    real(dp) :: dx, offsetx, offsety, widthcolor
    real(dp) :: heightcbox, xs, chunit, xt, yt
    real(dp), dimension(:), allocatable :: u
    character(len=20) :: ch
    character(len=8) :: frmt


    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in plot_color_contour_2D: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(sysvector) ) then
      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in plot_color_contour_2D: sysvector has not been created '
        stop
      end if
      allocate(u(problem%maxnoddegfd*mesh%maxelnumnod))
    else if ( present(vector) ) then
      if ( .not. vector%created ) then
        write(*,'(/a/)') &
          'Error in plot_color_contour_2D: vector has not been created '
        stop
      end if
      allocate(u(problem%maxvecnoddegfd*mesh%maxelnumnod))
    end if

    if ( present(projection) ) then
      proj = projection
    else
      proj = .false.
    end if

    if ( present(degfd) ) then
      deg = degfd
    else
      deg = 1
    end if

    if ( plot_options%usefloating ) then
      frmt = '(f9.2)'
    else
      frmt = '(es9.2)'
    end if

!   check groups
    if ( present(groups) ) then
      if ( any ( groups < 1 ) .or. any ( groups > mesh%nelgrp ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: parameter groups in the heading of plot_color_contour_2D ', &
          'is out of range: some groups are < 1 or larger than ', mesh%nelgrp
        stop
      end if
      allocate(lgroups(size(groups)))
      lgroups = groups
    else
      allocate(lgroups(mesh%nelgrp))
      lgroups = [ (i,i=1,mesh%nelgrp) ]
    end if

!   check element

    if ( any ( mesh%element(lgroups)%sidnumnod == 0 ) ) then
      write(*,'(/2a,i0/)') &
        'Error plot_color_contour: some of these element shapes ', &
        ' cannot be plotted = ', mesh%element(:)%elshape
      stop
    end if

!   layer

    if ( present(layer) ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in plot_color_contour: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a),i0/)') &
          'Error in plot_color_contour: ',&
          ' layer < 1 or layer > number of layers = ', &
          problem%numlayers
        stop
      end if
      wkl = btest(problem%nodlayers,layer-1)
    else if ( problem%numlayers > 0 ) then
      write(*,'(3(/a)/)') &
        'Error in plot_color_contour: ', &
        ' layers have been defined and the layer keyword is not present', &
        ' layers not implemented for this routine without the layer keyword'
      stop
    end if

    open(unit=unit_figplot,file=filename)

!   write header

    call fig_header

!   write colors

    ncolors = plot_options%numlevels

    call fig_colors ( ncolors )

    nintervals = ncolors + 1

!   fit plot into figsize x figsize

    xmin = max ( minval ( mesh%coor(:,1) ) + plot_options%xmin_coor_add, &
                 plot_options%xmin )
    xmax = min ( maxval ( mesh%coor(:,1) ) + plot_options%xmax_coor_add, &
                 plot_options%xmax )
    ymin = max ( minval ( mesh%coor(:,2) ) + plot_options%ymin_coor_add, &
                 plot_options%ymin )
    ymax = min ( maxval ( mesh%coor(:,2) ) + plot_options%ymax_coor_add, &
                 plot_options%ymax )

    scalefac = plot_options%figsize / max ( xmax - xmin , ymax - ymin )

!   plot boundary if not all groups are present

    if ( plot_options%plotboundary .and. ( present(groups) .or. &
           ( problem%numinactivegroups > 0 .and. .not. proj ) ) ) then

      call plot_boundary_internal ( plot_options, mesh, scalefac, xmin, xmax, &
        ymin, ymax )

    end if

!   filename in title

    call print_head ( plot_options, filename )

!   find maximum and minimum value

    if ( present(sysvector) ) then

!     sysvector

      umax = 0
      umin = 0

      first = .true.

      do nodenr = 1, mesh%nnodes
        if ( .not. node_in_groups ( mesh, lgroups, nodenr ) ) cycle
        if ( .not. node_in_groups ( mesh, problem%activegroups, nodenr ) &
               .and. .not. proj ) cycle
        if ( present(layer) ) then
          if ( .not. wkl(nodenr) ) cycle
        end if
        if ( present(physq) ) then
          call pos_array_node ( problem, nodenr, dof, pos, [physq], layer )
        else
          call pos_array_node ( problem, nodenr, dof, pos, layer=layer )
        end if
        if ( first ) then
          first = .false.
          umax=sysvector%u(pos(deg))
          umin=sysvector%u(pos(deg))
        else
          umax=max(umax,sysvector%u(pos(deg)))
          umin=min(umin,sysvector%u(pos(deg)))
        end if
      end do

    else if ( present(vector) ) then

      if ( vector%elementwise ) then

!       vector elementwise

        umax = 0
        umin = 0

        first = .true.

        do elgrp = 1, mesh%nelgrp

          if ( .not. any ( lgroups == elgrp ) ) cycle

!         skip inactive groups
          if ( any( elgrp == problem%inactivegroups ) .and. .not. proj ) cycle

          numnod = mesh%element(elgrp)%numnod

          do elem = 1, mesh%grpnumel(elgrp)

!           get element vector
            call get_vector ( mesh, problem, vector, elgrp, elem, u, &
              sloppy=.true. )
            bp = ( deg - 1 ) * numnod
            if ( first ) then
              first = .false.
              umax=maxval(u(bp+1:bp+numnod))
              umin=minval(u(bp+1:bp+numnod))
            else
              umax=max(umax,maxval(u(bp+1:bp+numnod)))
              umin=min(umin,minval(u(bp+1:bp+numnod)))
            end if

          end do

        end do

      else

!       vector in nodes

        umax = 0
        umin = 0

        first = .true.

        do nodenr = 1, mesh%nnodes
          if ( .not. node_in_groups ( mesh, lgroups, nodenr ) ) cycle
          if ( .not. node_in_groups ( mesh, problem%activegroups, nodenr ) &
                 .and. .not. proj ) cycle
          if ( present(layer) ) then
            if ( .not. wkl(nodenr) ) cycle
          end if
          call pos_array_vec_node ( problem, nodenr, dof, posv, vector%vec, &
            layer )
          if ( first ) then
            first = .false.
            umax=vector%u(posv(deg))
            umin=vector%u(posv(deg))
          else
            umax=max(umax,vector%u(posv(deg)))
            umin=min(umin,vector%u(posv(deg)))
          end if
        end do

      end if

    end if

!   determine interval between levels

    if ( ( umax > plot_options%maxvalue .and. umin < plot_options%minvalue ) &
            .or. plot_options%forcemaxmin ) then
!     one extra interval at the top and bottom
      umx = plot_options%maxvalue
      umn = plot_options%minvalue
      delta = ( umx - umn ) / ( nintervals - 2 )
      corr = 1   ! shift 1 to make all 1 <= interval <= nintervals
    else if ( umax > plot_options%maxvalue ) then
!     one extra interval at the top
      umx = plot_options%maxvalue
      umn = umin
      delta = ( umx - umn ) / ( nintervals - 1 )
      corr = 0
    else if ( umin < plot_options%minvalue ) then
!     one extra interval at the bottom
      umx = umax
      umn = plot_options%minvalue
      delta = ( umx - umn ) / ( nintervals - 1 )
      corr = 1   ! shift 1 to make all 1 <= interval <= nintervals
    else
!     all intervals between max and min
      umx = umax
      umn = umin
      delta = ( umax - umin ) / nintervals
      corr = 0
    end if

!   contours and boundary

    do grp = 1, size(lgroups)

      elgrp = lgroups(grp)

!     skip inactive groups
      if ( any( elgrp == problem%inactivegroups ) .and. .not. proj ) cycle

      do elem = 1, mesh%grpnumel(elgrp)

!       skip elements not in layer
        if ( present(layer) ) then
          if ( .not. all ( wkl(mesh%topology(elgrp)%a(:,elem)) ) ) cycle
        end if

!       check whether element in plot window: otherwise skip element

        elnumnod = mesh%element(elgrp)%numnod
        elnodes(1:elnumnod)  = mesh%topology(elgrp)%a(:,elem)

        elcoor(1,1:elnumnod) = mesh%coor(elnodes(1:elnumnod),1)
        elcoor(2,1:elnumnod) = mesh%coor(elnodes(1:elnumnod),2)

        if ( any( elcoor(1,:elnumnod) < xmin ) .or. &
             any( elcoor(1,:elnumnod) > xmax ) .or. &
             any( elcoor(2,:elnumnod) < ymin ) .or. &
             any( elcoor(2,:elnumnod) > ymax ) ) cycle

!       contour plot per element

        if ( present(sysvector) ) then

!         get element vector
          if ( present(physq) ) then
            call get_sysvector ( mesh, problem, sysvector, elgrp, elem, u, &
              [physq], ndofu=ndofu, sloppy=.true., layer=layer )
          else
            call get_sysvector ( mesh, problem, sysvector, elgrp, elem, u, &
              ndofu=ndofu, sloppy=.true., layer=layer )
          end if

        else if ( present(vector) ) then

!         get element vector
          call get_vector ( mesh, problem, vector, elgrp, elem, u, &
            ndofu=ndofu, sloppy=.true., layer=layer )

        end if

        numnod = mesh%element(elgrp)%numnod

        if ( ndofu < deg * numnod ) then
          write(*,'(/a/a/)') &
            'Error in plot_color_contour_2D:', &
            ' not enough unknowns in nodes for plotting'
          stop
        else if ( elem == 1 .and. mod(ndofu,numnod) /= 0 ) then
          write(*,'(/6(a/),a,i0/)') &
            'Warning in plot_color_contour_2D:', &
            ' number of unknowns in the sysvector or vector ', &
            ' is not a multiple of the number of nodes. Plot_color_contour', &
            ' assumes that the vector/sysvector is defined ', &
            ' in all nodes having the same number of degrees of freedom.', &
            ' Know what you are doing, I warned you!', &
            ' Element group ', elgrp
        end if

        do subelem = 1, mesh%element(elgrp)%numsubdiv

!         process each subtriangle individually

!         get nodes and coordinates of vertices

          nodes(1:3) = mesh%element(elgrp)%subtopol(:,subelem)

          nodenrs = mesh%topology(elgrp)%a(nodes(1:3),elem)

          xvert(1,1:3) = mesh%coor(nodenrs,1)
          xvert(2,1:3) = mesh%coor(nodenrs,2)

!         values in the vertices

          values = u ( ( deg - 1 ) * numnod + nodes(1:3) )

          valmax = maxval(values)
          valmin = minval(values)

          if ( valmax - valmin <= 1.e-14_dp * (umax - umin) ) then

!           value is constant

            constant = .true.
            nodmin = 0
            nodmax = 0
            nodmed = 0
            valmed = valmin

          else

!           other cases

            constant = .false.
            nodmin = maxval(minloc(values))
            nodmax = maxval(maxloc(values))
            nodmed = 6 - nodmin - nodmax
            valmed = values(nodmed)

          end if

!         determine interval between contour levels for the vertices

          intervalmax = int ( ( valmax - umn ) / delta + 1 + corr )
            intervalmax = min ( intervalmax, nintervals )
            intervalmax = max ( intervalmax, 1 )
          intervalmin = int ( ( valmin - umn ) / delta + 1 + corr )
            intervalmin = min ( intervalmin, nintervals )
            intervalmin = max ( intervalmin, 1 )
          intervalmed = int ( ( valmed - umn ) / delta + 1 + corr )
            intervalmed = min ( intervalmed, nintervals )
            intervalmed = max ( intervalmed, 1 )

          if ( intervalmax == intervalmin .or. constant ) then

!           no contour

            cycle

          end if

          if ( intervalmed > intervalmin .and. .not. constant ) then

!           contours near min point
!
!                       med
!                      .    .
!                     .      .
!                    . \      .
!                   .   \      .
!                 .  \   \      .
!                .    \   \      .
!               min ... ...  ... max
!

            npoints = 2

            numcont = intervalmed - intervalmin

            do cont = 1, numcont

              f = ( ( intervalmin - corr + cont - 1 ) * delta + umn - valmin ) &
                   / ( valmax - valmin )
              x(:,1) = (1-f)*xvert(:,nodmin)+f*xvert(:,nodmax)
              f = ( ( intervalmin - corr + cont - 1 ) * delta + umn - valmin ) &
                   / ( valmed - valmin )
              x(:,2) = (1-f)*xvert(:,nodmin)+f*xvert(:,nodmed)

              x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
              x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

              call fig_polyline ( pen_color=intervalmin+cont-1+31, &
                 x=x(:,1:npoints), npoints=npoints, test=.true. )

            end do

          end if

          if ( intervalmed < intervalmax .and. .not. constant ) then

!           contours near max point
!
!                       med
!                     .     .
!                    .       .
!                   .       / .
!                  .       /   .
!                 .       /   / .
!                .       /   /   .
!               min ... ..... ... max
!
            npoints = 2

            numcont = intervalmax - intervalmed

            do cont = 1, numcont

              f = ( ( intervalmax - corr - cont ) * delta + umn - valmin ) &
                   / ( valmax - valmin )
              x(:,1) = (1-f)*xvert(:,nodmin)+f*xvert(:,nodmax)
              f = ( ( intervalmax - corr - cont ) * delta + umn - valmed ) &
                   / ( valmax - valmed )
              x(:,2) = (1-f)*xvert(:,nodmed)+f*xvert(:,nodmax)

              x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
              x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

              call fig_polyline ( pen_color=intervalmax-cont+31, &
                 x=x(:,1:npoints), npoints=npoints, test=.true. )

            end do

          end if

        end do

        if ( plot_options%plotboundary2 ) then

!         boundary

          do side = 1, mesh%element(elgrp)%numsides

!           plot side on boundary or between element groups (useful for surfaces
!           in 3D).

            if ( mesh%sidelem(elgrp)%a(side,elem,4) <= 0 .or. &
                 elgrp /= mesh%sidelem(elgrp)%a(side,elem,1) ) then

              npoints= mesh%element(elgrp)%sidnumnod

              nodes(1:npoints) = mesh%topology(elgrp)%a(&
                   &mesh%element(elgrp)%sidnod(1:npoints,side),elem)

              x(1,1:npoints) = mesh%coor(nodes(1:npoints),1)
              x(2,1:npoints) = mesh%coor(nodes(1:npoints),2)

              x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
              x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

              call fig_polyline ( 0, npoints, x(:,1:npoints) )

            end if

          end do

        end if

      end do
    end do

!   contour values

    if ( plot_options%plotlevelvalues ) then

!     plot colors and values of levels

      dx = plot_options%figsize / 150  ! scale factor (=1 for figsize=150mm)

      offsetx          = 10  * dx  ! initial offset in x direc
      offsety          = 150 * dx  ! initial offset in y direc
      widthcolor       = 10  * dx  ! width of the color lines

!     left x position of the colorlevels
      xs = plot_options%figsize + offsetx

      npoints = 2

      fontsize = plot_options%fontsize      ! fontsize in points
      chunit = real(fontsize,dp)/MMFIGUNITS ! unit scaling with the fontsize

!     height of a single colorlevelbox
      heightcbox = 20 * chunit

      do interval = 1, nintervals - 1

!       color line for levels

        x(:,1) = [ xs, offsety - ( nintervals - interval ) * heightcbox ]
        x(:,2) = x(:,1) + [ widthcolor, 0._dp ]

        x(2,1:npoints) = plot_options%figsize - x(2,1:npoints)

        x(1,1:npoints) = x(1,1:npoints) + plot_options%shiftbarx
        x(2,1:npoints) = x(2,1:npoints) - plot_options%shiftbary

        call fig_polyline ( pen_color=interval+31, &
          x=x(:,1:npoints), npoints=npoints )

!       print level value

        write( ch, frmt ) umn + delta * ( interval - corr )

        xt = x(1,2) + 5*chunit
        yt = x(2,2) + 5*chunit

        call fig_text ( fontsize=fontsize, x=xt, y=yt, text=ch(1:9) )

      end do

!     print max/min values

      if ( plot_options%printmaxminvalues ) then

        write( ch, '(a,'//frmt(2:) ) 'max:', umax

        xt = xs
        yt = offsety - ( nintervals + 1 ) * heightcbox
        yt = plot_options%figsize - yt
        xt = xt + plot_options%shiftbarx
        yt = yt - plot_options%shiftbary

        call fig_text ( fontsize=fontsize, x=xt, y=yt, text=ch(1:13) )

        write( ch, '(a,'//frmt(2:) ) 'min:', umin

        yt = yt + 20*chunit

        call fig_text ( fontsize=fontsize, x=xt, y=yt, text=ch(1:13) )

      end if

    end if

!   deallocate memory

    if ( allocated(u) ) deallocate(u)
    deallocate(lgroups)

    close(unit=unit_figplot)

  end subroutine plot_color_contour_2D


! vector plot

  subroutine plot_vector ( plot_options, mesh, problem, filename, physq, &
    layer, degfd, surfaces, groups, sysvector, vector, append )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   the physical quantity to be plotted
    integer, intent(in), optional :: physq

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered.
    integer, intent(in), optional :: layer

!   the degrees of freedom to be plotted as a vector (within the physical
!   quantity when physq is also present).
!   the dimension of degfd must be equal to the space dimension.
!   in 3D: setting one of the components to zero means this component is
!   not included in the projection to 2D of the vector (as if being zero).
    integer, intent(in), dimension(:), optional :: degfd

!   when present the vectors will be plotted on the given surfaces only
    integer, dimension(:), intent(in), optional :: surfaces

!   if present: the element groups to be plotted
!   For example groups=(/2,4/), will plot for element groups 2 and 4 only
!   Default: all groups are plotted
    integer, dimension(:), intent(in), optional :: groups


!   sysvector to be plotted
    type(sysvector_t), intent(in), optional :: sysvector

!   vector to be plotted
    type(vector_t), intent(in), optional :: vector

!   if .true.: append plot to existing file in order to overlay plots.
!   default=.false.
    logical, intent(in), optional :: append


    real(dp), dimension(3) :: xunitvec, yunitvec
    type(mesh_t) :: lmesh


    call check ( mesh, 'plot_vector' )
    call check_blend ( mesh, 'plot_vector' )
    call check ( problem, 'plot_vector', mesh )

    if ( present(degfd) ) then
      if ( size(degfd) /= mesh%ndim ) then
        write(*,'(/a,i0/)') &
          'Error in plot_vector: degfd must have dimension ', mesh%ndim
        stop
      end if
    end if

    if ( present(surfaces) ) then
      if ( any(surfaces < 1) .or. any(surfaces > mesh%nsurfaces) ) then
        write(*,'(/2a/(a,i0)/a,20i0/)') &
          'Error: surfaces in the heading of plot_vector are ', &
          'out of range.', &
          'the number of surfaces is ', mesh%nsurfaces, &
          'whereas surfaces are ', surfaces
        stop
      end if
    end if

!   test dimension of mesh

    if ( mesh%ndim == 3 ) then

!     3D mesh: project coordinates onto a 2D plane.

      if ( present(vector) ) then
        if ( vector%elementwise ) then
          write( *, '(/2(a/))') &
            'Error plot_vector: ', &
            ' in 3D vectors stored elementwise cannot be plotted'
          stop
        end if
      end if

!     compute projection vectors

      call projection_vectors ( &
        plot_options%viewpoint, plot_options%ydirection, xunitvec, yunitvec )

      if ( present(surfaces) ) then

!       create new mesh from surface topology

        if ( present(groups) ) then
          write(*,'(/a/a,i0/)') &
            ' Error in plot_vector: the heading parameter groups ', &
            ' cannot be used together with the heading parameter surfaces.'
          stop
        end if

        if ( problem%numinactivegroups > 0 ) then
          write(*,'(/a/a,i0/)') &
            ' Warning in plot_vector: inactive groups cannot be ', &
            ' used together with the heading parameter surfaces.'
        end if

        call generate_surface_mesh ( mesh, lmesh, surfaces )

      else

        write( *, '(/3(a/))') &
          'Error plot_vector: ', &
          ' in 3D only plotting on surfaces is supported', &
          ' use surface= argument to plot on surface only '
        stop

      end if

!     project coordinates in new mesh

      call project_coordinates ( xunitvec, yunitvec, lmesh%coor )

      call plot_vector_2D ( plot_options, lmesh, problem, filename, &
        physq, layer, degfd, append=append, sysvector=sysvector, &
        vector=vector, xunitvec=xunitvec, yunitvec=yunitvec,  &
        projection=.true.  )

      call delete ( lmesh )

    else

!     2D mesh

      call plot_vector_2D ( plot_options, mesh, problem, filename, &
        physq, layer, degfd, groups=groups, append=append, &
        sysvector=sysvector, vector=vector, xunitvec=xunitvec, &
        yunitvec=yunitvec )

    end if

  end subroutine plot_vector


! vector plot

  subroutine plot_vector_2D ( plot_options, mesh, problem, filename, physq, &
    layer, degfd, groups, append, sysvector, vector, xunitvec, yunitvec, &
    projection )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh
    type(problem_t), intent(in) :: problem

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   the physical quantity to be plotted
    integer, intent(in), optional :: physq

!   if layer is present only degrees of freedom with respect to the specified
!   layer are considered.
    integer, intent(in), optional :: layer

!   the degree of freedom to be plotted (within the physical quantity if
!   physq is also present)
    integer, intent(in), dimension(:), optional :: degfd

!   the element groups to be plotted
    integer, dimension(:), intent(in), optional :: groups

!   if .true.: append plot to existing file in order to overlay plots.
!   default is .false.
    logical, intent(in), optional :: append

    type(sysvector_t), intent(in), optional :: sysvector
    type(vector_t), intent(in), optional :: vector

!   unit vectors of the new x and y directions in the projected plane
    real(dp), dimension(3), intent(in), optional :: xunitvec, yunitvec

!   logical indicating whether routine is used for 3D projection (.true.) or
!   not (.false.). The default is .false.
    logical, intent(in), optional :: projection


    logical :: fileexists, proj, fileappend
    logical :: elementwise, firstskip
    logical, dimension(mesh%nnodes) :: wkl
    integer :: elgrp, elem, side, nodes(MAXPOLYLINEPOINTS)
    integer :: node, numnod, nodenr, npoints, tmp(mesh%nelgrp)
    integer :: pos(problem%maxnoddegfd), dof, ndofu, grp
    integer :: posv(problem%maxvecnoddegfd)
    integer :: work(mesh%nnodes), i, ldeg(3)
    integer, dimension(:), allocatable :: lgroups, ldegfd

    real(dp) :: x(2,MAXPOLYLINEPOINTS)
    real(dp) :: scalefac, xmin, xmax, ymin, ymax
    real(dp), dimension(:), allocatable :: u


    if ( present(projection) ) then
      proj = projection
    else
      proj = .false.
    end if

    if ( present(degfd) ) then
      allocate(ldegfd(size(degfd)))
      ldegfd = degfd
    else if ( proj ) then
      allocate(ldegfd(3))
      ldegfd = [1,2,3]
    else
      allocate(ldegfd(2))
      ldegfd = [1,2]
    end if

    if ( size(ldegfd) == 3 ) then
      if ( .not. present(xunitvec) .or. .not. present(yunitvec) ) then
        write( *, '(/2(a/))') &
          'Error plot_vector_2D: ', &
          ' for 3D vectors xunitvec and yunitvec must be present'
        stop
      end if
!     in 3D plot only vectors in points that are in the topology are used
!     this allows for a fake mesh that consists of only a part of the
!     original mesh
      work = 0
      do elgrp = 1, mesh%nelgrp
        do elem = 1, mesh%grpnumel(elgrp)
          work ( mesh%topology(elgrp)%a(:,elem) ) = 1
        end do
      end do
!     avoid zero index in position vector
      where ( ldegfd == 0 )
        ldeg = 1
      else where
        ldeg = ldegfd
      end where
    end if

    if ( present(vector) .and. present(sysvector) ) then
      write(*,'(/a/a/)') &
        'Error in plot_vector_2D: only one of vector or sysvector ', &
        ' can be present in the heading.'
      stop
    end if

    if ( present(sysvector) ) then
      if ( .not. sysvector%created ) then
        write(*,'(/a/)') &
          'Error in plot_vector_2D: sysvector has not been created '
        stop
      end if
    else if ( present(vector) ) then
      if ( .not. vector%created ) then
        write(*,'(/a/)') &
          'Error in plot_vector_2D: vector has not been created '
        stop
      end if
    end if

    elementwise = .false.

    if ( present(vector) ) then
      if ( vector%elementwise ) then
        elementwise = .true.
        do elgrp = 1, mesh%nelgrp
          tmp(elgrp) = sum ( problem%vec_elnumdegfd(elgrp)%a(:,vector%vec) )
        end do
        allocate(u(maxval(tmp)))
      end if
    end if

!   check groups
    if ( present(groups) ) then
      if ( any ( groups < 1 ) .or. any ( groups > mesh%nelgrp ) ) then
        write(*,'(/a/a,i0/)') &
          'Error: parameter groups in the heading of plot_vector_2D is ', &
          'out of range: some groups are < 1 or larger than ', mesh%nelgrp
        stop
      end if
      allocate(lgroups(size(groups)))
      lgroups = groups
    else
      allocate(lgroups(mesh%nelgrp))
      lgroups = [ (i,i=1,mesh%nelgrp) ]
    end if

!   check element

    if ( any ( mesh%element(lgroups)%sidnumnod == 0 ) ) then
      write(*,'(/a/a,i0/)') &
        'Error plot_vector: some of these element shapes cannot be plotted.', &
        ' elshape = ', mesh%element(:)%elshape
      stop
    end if

!   layer

    if ( present(layer) ) then
      if ( problem%numlayers == 0 ) then
        write(*,'(3(/a)/)') &
          'Error in plot_vector: ',&
          ' specification of a layer ', &
          ' is only available if layers have been defined.'
        stop
      end if
      if ( layer < 1 .or. layer > problem%numlayers ) then
        write(*,'(2(/a),i0/)') &
          'Error in plot_vector: ',&
          ' layer < 1 or layer > number of layers = ', &
          problem%numlayers
        stop
      end if
      wkl = btest(problem%nodlayers,layer-1)
    else if ( problem%numlayers > 0 ) then
      write(*,'(3(/a)/)') &
        'Error in plot_vector: ', &
        ' layers have been defined and the layer keyword is not present', &
        ' layers not implemented for this routine without the layer keyword'
      stop
    end if

    inquire ( file=filename, exist=fileexists )

    if ( present(append) ) then
      fileappend = append
    else
      fileappend = .false.
    end if

    if ( fileexists .and. fileappend ) then
      open ( unit=unit_figplot, file=filename, status='old', position='append' )
    else
      open ( unit=unit_figplot, file=filename )
      call fig_header
      call print_head ( plot_options, filename )
    end if

!   fit plot into figsize x figsize

    xmin = max ( minval ( mesh%coor(:,1) ) + plot_options%xmin_coor_add, &
                 plot_options%xmin )
    xmax = min ( maxval ( mesh%coor(:,1) ) + plot_options%xmax_coor_add, &
                 plot_options%xmax )
    ymin = max ( minval ( mesh%coor(:,2) ) + plot_options%ymin_coor_add, &
                 plot_options%ymin )
    ymax = min ( maxval ( mesh%coor(:,2) ) + plot_options%ymax_coor_add, &
                 plot_options%ymax )

    scalefac = plot_options%figsize / max ( xmax - xmin , ymax - ymin )


!   plot boundary if not all groups are present

    if ( plot_options%plotboundary .and. ( present(groups) .or. &
           ( problem%numinactivegroups > 0 .and. .not. proj ) ) ) then

      call plot_boundary_internal ( plot_options, mesh, scalefac, xmin, xmax, &
        ymin, ymax )

    end if


!   plot boundary and elementwise defined vector

    do grp = 1, size(lgroups)

      elgrp = lgroups(grp)

!     skip inactive groups
      if ( any( elgrp == problem%inactivegroups ) .and. .not. proj ) cycle

      npoints= mesh%element(elgrp)%sidnumnod

      do elem = 1, mesh%grpnumel(elgrp)

        if ( plot_options%plotboundary2 ) then

          do side = 1, mesh%element(elgrp)%numsides

!           plot side on boundary or between element groups (useful for surfaces
!           in 3D).

            if ( mesh%sidelem(elgrp)%a(side,elem,4) <= 0 .or. &
                 elgrp /= mesh%sidelem(elgrp)%a(side,elem,1) ) then

              nodes(1:npoints) = mesh%topology(elgrp)%a(&
                   &mesh%element(elgrp)%sidnod(1:npoints,side),elem)

              x(1,1:npoints) = mesh%coor(nodes(1:npoints),1)
              x(2,1:npoints) = mesh%coor(nodes(1:npoints),2)

              if ( any( x(1,:npoints) < xmin ) .or. &
                   any( x(1,:npoints) > xmax ) .or. &
                   any( x(2,:npoints) < ymin ) .or. &
                   any( x(2,:npoints) > ymax ) ) cycle

              x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
              x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

              call fig_polyline ( 0, npoints, x(:,1:npoints) )

            end if

          end do

        end if

!       vector type defined elementwise must be plotted element wise

        if ( present(vector) .and. elementwise ) then

!         get element vector
          call get_vector ( mesh, problem, vector, elgrp, elem, u, &
            ndofu=ndofu, sloppy=.true. )

          numnod = mesh%element(elgrp)%numnod

          if ( elem == 1 .and. mod(ndofu,numnod) /= 0 ) then
            write(*,'(/6(a/),2(a,i0/))') &
              'Warning in plot_vector_2D:', &
              ' number of unknowns in the vector ', &
              ' is not a multiple of the number of nodes.', &
              ' Plot_vector assumes that the vector is defined ', &
              ' in all nodes having the same number of degrees of freedom.', &
              ' Know what you are doing, I warned you!', &
              ' Element group ', elgrp, &
              ' Vector number vector%vec = ', vector%vec
          end if

          do node = 1, numnod

            nodenr = mesh%topology(elgrp)%a(node,elem)

            x(:,1) = mesh%coor(nodenr,:)
            x(:,2) = x(:,1) + &
               plot_options%scalevector * u ( ( ldegfd - 1 ) * numnod + node )

            if ( any ( x(:,1) < [xmin,ymin] ) .or. &
                 any ( x(:,1) > [xmax,ymax] ) )  cycle

            x(1,1:2) = scalefac * ( x(1,1:2) - xmin )
            x(2,1:2) = scalefac * ( ymax - x(2,1:2) )

            call fig_vector ( pen_color=plot_options%vectorcolor, x=x(:,1:2) )

          end do

        end if

      end do
    end do

    if ( present(sysvector) ) then

!     sysvector in nodes

      firstskip = .true.

      do nodenr = 1, mesh%nnodes

!       check group

        if ( .not. node_in_groups ( mesh, lgroups, nodenr ) ) cycle
        if ( .not. node_in_groups ( mesh, problem%activegroups, nodenr ) &
               .and. .not. proj ) cycle

!       layer

        if ( present(layer) ) then
          if ( .not. wkl(nodenr) ) cycle
        end if

!       test whether node is in plot window

        if ( any ( mesh%coor(nodenr,1:2) < [xmin,ymin] ) .or. &
             any ( mesh%coor(nodenr,1:2) > [xmax,ymax] ) )  cycle

        if ( size(ldegfd) == 3 ) then
          if ( work(nodenr) == 0 ) cycle  ! only nodes that are in the topology
        end if

        if ( present(physq) ) then
          call pos_array_node ( problem, nodenr, dof, pos, &
            physqarr=[physq], layer=layer )
        else
          call pos_array_node ( problem, nodenr, dof, pos, layer=layer )
        end if

        if ( any(ldegfd > dof) ) then
          if ( firstskip ) then
            write(*,'(/3(a/))') &
              'Warning in plot_vector_2D:', &
              ' there are nodes where the required unknowns in the ', &
              ' sysvector are undefined. These nodes are skipped. '
            firstskip = .false.
          end if
          cycle
        end if

        x(:,1) = mesh%coor(nodenr,1:2)

        if ( size(ldegfd) == 2 ) then
!         2D vectors
          x(:,2) = &
                 x(:,1) + plot_options%scalevector * sysvector%u( pos(ldegfd) )
        else if ( size(ldegfd) == 3 ) then
!         3D vectors
          x(1,2) = x(1,1) + plot_options%scalevector * &
            sum ( sysvector%u( pos(ldeg) ) * xunitvec, mask=( ldegfd > 0 ) )
          x(2,2) = x(2,1) + plot_options%scalevector * &
            sum ( sysvector%u( pos(ldeg) ) * yunitvec, mask=( ldegfd > 0 ) )
        end if

        x(1,1:2) = scalefac * ( x(1,1:2) - xmin )
        x(2,1:2) = scalefac * ( ymax - x(2,1:2) )

        call fig_vector ( pen_color=plot_options%vectorcolor, x=x(:,1:2) )

      end do

    else if ( present(vector) .and. .not. elementwise ) then

!     vector in nodes

      firstskip = .true.

      do nodenr = 1, mesh%nnodes

!       check group

        if ( .not. node_in_groups ( mesh, lgroups, nodenr ) ) cycle
        if ( .not. node_in_groups ( mesh, problem%activegroups, nodenr ) &
               .and. .not. proj ) cycle

!       layer

        if ( present(layer) ) then
          if ( .not. wkl(nodenr) ) cycle
        end if

!       test whether node is in plot window

        if ( any ( mesh%coor(nodenr,1:2) < [xmin,ymin] ) .or. &
             any ( mesh%coor(nodenr,1:2) > [xmax,ymax] ) )  cycle

        if ( size(ldegfd) == 3 ) then
          if ( work(nodenr) == 0 ) cycle  ! only nodes that are in the topology
        end if

        call pos_array_vec_node ( problem, nodenr, dof, posv, vector%vec, &
          layer )

        if ( any(ldegfd > dof) ) then
          if ( firstskip ) then
            write(*,'(/3(a/))') &
              'Warning in plot_vector_2D:', &
              ' there are nodes where the required unknowns in the ', &
              ' vector are undefined. These nodes are skipped. '
            firstskip = .false.
          end if
          cycle
        end if

        x(:,1) = mesh%coor(nodenr,1:2)

        if ( size(ldegfd) == 2 ) then
!         2D vectors
          x(:,2) = x(:,1) + plot_options%scalevector * &
                 vector%u( posv(ldegfd) )
        else if ( size(ldegfd) == 3 ) then
!         3D vectors
          x(1,2) = x(1,1) + plot_options%scalevector * sum ( &
            vector%u( posv(ldeg) ) * xunitvec, mask=( ldegfd > 0 ) )
          x(2,2) = x(2,1) + plot_options%scalevector * sum ( &
            vector%u( posv(ldeg) ) * yunitvec, mask=( ldegfd > 0 ) )
        end if

        x(1,1:2) = scalefac * ( x(1,1:2) - xmin )
        x(2,1:2) = scalefac * ( ymax - x(2,1:2) )

        call fig_vector ( pen_color=plot_options%vectorcolor, x=x(:,1:2) )

      end do

    end if

!   deallocate memory

    if ( allocated(u) ) deallocate(u)
    deallocate(lgroups,ldegfd)

    close(unit=unit_figplot)

  end subroutine plot_vector_2D


! plot boundary (external use)

  subroutine plot_boundary ( plot_options, mesh, filename, append )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

!   if .true.: append plot to existing file in order to overlay plots.
!   default is .false.
    logical, intent(in), optional :: append


    integer :: zerogroups(0)


!   Trick: plot mesh without groups

    call plot_mesh ( plot_options, mesh, filename, groups=zerogroups, &
      append=append )

  end subroutine plot_boundary


! plot boundary (internal use only)

  subroutine plot_boundary_internal ( plot_options, mesh, scalefac, xmin, &
    xmax, ymin, ymax )

    type(plot_options_t), intent(in) :: plot_options
    type(mesh_t), intent(in) :: mesh
    real(dp), intent(in) :: scalefac, xmin, xmax, ymin, ymax


    integer :: npoints, elgrp, elem, side, nodes(MAXPOLYLINEPOINTS)
    real(dp) :: x(2,MAXPOLYLINEPOINTS)

!   check element

    if ( any ( mesh%element(:)%sidnumnod == 0 ) ) then
      write(*,'(/2a,i0/)') &
        'Error plot_boundary_internal: some of these element shapes ', &
        ' cannot be plotted = ', mesh%element(:)%elshape
      stop
    end if

    do elgrp = 1, mesh%nelgrp
      do elem = 1, mesh%grpnumel(elgrp)

        do side = 1, mesh%element(elgrp)%numsides

!         plot side on boundary

          if ( mesh%sidelem(elgrp)%a(side,elem,4) <= 0 ) then

            npoints= mesh%element(elgrp)%sidnumnod

            nodes(1:npoints) = mesh%topology(elgrp)%a(&
                 &mesh%element(elgrp)%sidnod(1:npoints,side),elem)

            x(1,1:npoints) = mesh%coor(nodes(1:npoints),1)
            x(2,1:npoints) = mesh%coor(nodes(1:npoints),2)

            if ( any( x(1,:npoints) < xmin ) .or. &
                 any( x(1,:npoints) > xmax ) .or. &
                 any( x(2,:npoints) < ymin ) .or. &
                 any( x(2,:npoints) > ymax ) ) cycle

            x(1,1:npoints) = scalefac * ( x(1,1:npoints) - xmin )
            x(2,1:npoints) = scalefac * ( ymax - x(2,1:npoints) )

            call fig_polyline ( plot_options%boundarycolor, npoints, &
              x(:,1:npoints) )

          end if

        end do

      end do
    end do

  end subroutine plot_boundary_internal


! print header

  subroutine print_head ( plot_options, filename )

    type(plot_options_t), intent(in) :: plot_options

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename

    call print_title ( plot_options, filename )

    call print_date_and_time ( plot_options )

  end subroutine print_head


! print a title

  subroutine print_title ( plot_options, filename )

    type(plot_options_t), intent(in) :: plot_options

!   the filename for writing the plot to
    character (len=*), intent(in) :: filename


    real(dp) :: fontsize, chunit


    if ( plot_options%printfilename ) then

!     print filename in the title

      fontsize = plot_options%fontsize      ! fontsize in points
      chunit = real(fontsize,dp)/MMFIGUNITS ! unit scaling with the fontsize

      call fig_text ( fontsize=plot_options%fontsize, &
        x=0._dp, y=-15*chunit, text=filename, justification=0 )

    end if

  end subroutine print_title


! print date and time

  subroutine print_date_and_time ( plot_options )

    type(plot_options_t), intent(in) :: plot_options


    character (len=21) :: line
    integer :: values(8)
    real(dp) :: fontsize, chunit


    if ( plot_options%printdateandtime ) then

!     print date and time in the footer

      fontsize = plot_options%fontsize      ! fontsize in points
      chunit = real(fontsize,dp)/MMFIGUNITS ! unit scaling with the fontsize

      call date_and_time ( values=values )

      write ( line, '(2(i2.2,a),i4,a,3(i2.2,a))' ) &
        values(3), '-', values(2), '-', values(1), ',  ', &
        values(5), ':', values(6), ':', values(7)

      call fig_text ( fontsize=plot_options%fontsize, &
        x=plot_options%figsize, y=-15*chunit, text=line, &
        justification=2 )

    end if

  end subroutine print_date_and_time


! compute projection vectors

  subroutine projection_vectors ( viewpoint, ydirection, xunitvec, yunitvec )

!   direction of the viewpoint. For example viewpoint=(/1,1,1/) means that
!   every point will be projected on the plane with a normal in the direction
!   of (/1,1,1/) and the view is _from_ (/1,1,1/).
    real(dp), dimension(3), intent(in) :: viewpoint

!   co-ordinate direction that becomes the new y co-oordinate in the projected
!   plane, 1=x, 2=y, 3=z.
    integer, intent(in) :: ydirection

!   unit vectors of the new x and y directions in the projected plane
    real(dp), dimension(3), intent(out) :: xunitvec, yunitvec


    real(dp), dimension(3) :: n
    real(dp) :: length


    n = viewpoint / sqrt ( dot_product ( viewpoint, viewpoint ) )

    select case ( ydirection )
    case(1) ! x -> new y
      length = sqrt ( n(2)**2 + n(3)**2 )
      xunitvec = [ 0._dp, - n(3), n(2) ] / length
      yunitvec = [ 1._dp - n(1)**2, - n(1)*n(2), - n(1)*n(3) ] / length
    case(2) ! y -> new y
      length = sqrt ( n(1)**2 + n(3)**2 )
      xunitvec = [ n(3), 0._dp, - n(2) ] / length
      yunitvec = [ - n(1)*n(2), 1._dp - n(2)**2, - n(2)*n(3) ] / length
    case(3) ! z -> new y
      length = sqrt ( n(1)**2 + n(2)**2 )
      xunitvec = [ - n(2), n(1), 0._dp ] / length
      yunitvec = [ - n(1)*n(3), - n(2)*n(3), 1._dp - n(3)**2 ] / length
    case default
      call errormsg_case_default ( 'projection_vectors', &
        'ydirection', int_value=ydirection )
    end select

  end subroutine projection_vectors


! project coordinates in the projection plane

  subroutine project_coordinates ( xunitvec, yunitvec, coor )

!   unit vectors of the new x and y directions in the projected plane
    real(dp), dimension(3), intent(in) :: xunitvec, yunitvec

!   coordinates coor(nnodes,ndim)
    real(dp), dimension(:,:), intent(inout) :: coor


    integer :: ndim, node
    real(dp) :: xnew, ynew


    ndim = size(coor,2)

    do node = 1, size(coor,1)
      xnew = dot_product ( coor(node,:), xunitvec(:ndim) )
      ynew = dot_product ( coor(node,:), yunitvec(:ndim) )
      coor(node,1) = xnew
      coor(node,2) = ynew
    end do

  end subroutine project_coordinates


! generate surface mesh

  subroutine generate_surface_mesh ( mesh, lmesh, surfaces )

!   original mesh
    type(mesh_t), intent(in) :: mesh

!   surface mesh
    type(mesh_t), intent(inout) :: lmesh

!   surface to be used
    integer, dimension(:), intent(in) :: surfaces


    integer :: surface, nelgrp, surf, elnumnod, nelem

!   use multiple groups

    nelgrp = size ( surfaces )

!   layout like meshgen_basic

    lmesh%ndim = 2
    lmesh%nnodes = mesh%nnodes
    lmesh%nelgrp = nelgrp
    lmesh%nelem = sum ( mesh%surfaces(surfaces)%nelem )

    allocate( lmesh%element(nelgrp) )
    allocate( lmesh%element_blend(nelgrp,0) )
    allocate( lmesh%nnodes_blend(2) ); lmesh%nnodes_blend = [0,lmesh%nnodes]
    allocate( lmesh%grpnumel(nelgrp) )
    allocate( lmesh%topology(nelgrp) )

    do surf = 1, nelgrp

      surface = surfaces(surf)
      elnumnod = mesh%surfaces(surface)%element%numnod
      nelem = mesh%surfaces(surface)%nelem

      allocate( lmesh%topology(surf)%a(elnumnod,nelem) )

!     element type

      lmesh%element(surf)%elshape = mesh%surfaces(surface)%element%elshape
      lmesh%element(surf)%numnod = elnumnod
      lmesh%element(surf)%ndim = 2

!     topology

      lmesh%grpnumel(surf) = nelem
      lmesh%topology(surf)%a = mesh%surfaces(surface)%topology(:,:,2)

    end do

!   coordinates

    allocate ( lmesh%coor(mesh%nnodes,mesh%ndim) )

    lmesh%coor = mesh%coor

!   geometries etc.

    lmesh%npoints = 0
    allocate( lmesh%points(0) )
    lmesh%ncurves = 0
    allocate( lmesh%curves(0) )
    lmesh%nsurfaces = 0
    allocate( lmesh%surfaces(0) )
    lmesh%nvolumes = 0
    allocate( lmesh%volumes(0) )
    lmesh%nnodesets = 0
    allocate ( lmesh%nodesets(0) )
    lmesh%nelementsets = 0
    allocate ( lmesh%elementsets(0) )
    lmesh%nobjects = 0
    allocate ( lmesh%objects(0) )
    lmesh%nblocks = -1  ! trick to avoid building blocks
    allocate ( lmesh%blocks(0) )
    lmesh%meshgen = .true.

    call fill_mesh_parts ( lmesh )

  end subroutine generate_surface_mesh

! write a header

  subroutine fig_header

!   write header

    write(unit=unit_figplot,fmt='(8(a/),a)')   &
      '#FIG 3.2', &
      'Landscape', &
      'Center', &
      'Metric', &
      'A4', &
      '100.00', &
      'Single', &
      '-2', &
      '1200 2'

  end subroutine fig_header


! plot a polyline

  subroutine fig_polyline ( pen_color, npoints, x, test )

    integer, intent(in) :: pen_color, npoints

!   x: coordinates in mm
    real(dp), intent(in) :: x(2,npoints)

!   test equal coordinates
    logical, optional, intent(in) :: test

    if ( present(test) ) then
      if ( all( nint(MMFIGUNITS*x(:,1)) == nint(MMFIGUNITS*x(:,2)) ) ) return
    end if

    write(unit=unit_figplot,fmt='(a,i0,a,i0)') '2 1 0 1 ', pen_color, &
      ' 7 50 -1 -1 0.000 0 0 -1 0 0 ', npoints
    write(unit=unit_figplot,fmt='(12(1x,i0))') nint(MMFIGUNITS*x)

  end subroutine fig_polyline


! plot spline

  subroutine fig_spline ( pen_color, npoints, x, test )

    integer, intent(in) :: pen_color, npoints

!   x: coordinates in mm
    real(dp), intent(in) :: x(2,npoints)

!   test equal coordinates
    logical, optional, intent(in) :: test

    if ( present(test) ) then
      if ( all( nint(MMFIGUNITS*x(:,1)) == nint(MMFIGUNITS*x(:,2)) ) ) return
    end if

    write(unit=unit_figplot,fmt='(a,i0,a,i0)') '3 2 0 1 ', pen_color, &
      ' 7 50 -1 -1 0.000 0 0 0 ', npoints
    write(unit=unit_figplot,fmt='(12(1x,i0))') nint(MMFIGUNITS*x)
    write(unit=unit_figplot,fmt='(a)') '0.000 -1.000 0.000'

  end subroutine fig_spline


! plot a filled polygon

  subroutine fig_filled_polygon ( fill_color, npoints, x, test )

    integer, intent(in) :: fill_color, npoints

!   x: coordinates in mm
    real(dp), intent(in) :: x(2,npoints)

!   test equal coordinates
    logical, optional, intent(in) :: test

    if ( present(test) ) then
      if ( all( nint(MMFIGUNITS*x(:,1)) == nint(MMFIGUNITS*x(:,2)) ) ) return
    end if

    write(unit=unit_figplot,fmt='(a,i0,a,i0)') '2 3 0 0 0 ', fill_color, &
      ' 50 -1 20 0.000 0 0 -1 0 0 ', npoints
    write(unit=unit_figplot,fmt='(12(1x,i0))') nint(MMFIGUNITS*x)

  end subroutine fig_filled_polygon


! plot a vector

  subroutine fig_vector ( pen_color, x )

    integer, intent(in) :: pen_color

!   x: coordinates in mm
    real(dp), intent(in) :: x(2,2)

    real(dp) :: length

    length = sqrt( sum((x(:,2)-x(:,1))**2) )

    if ( length < 1 ) then
!     no arrow
      write(unit=unit_figplot,fmt='(a,i0,a)') '2 1 0 1 ', pen_color, &
        ' 7 50 -1 -1 0.000 0 0 -1 0 0 2'
      write(unit=unit_figplot,fmt='(4(1x,i0))') nint(MMFIGUNITS*x)
    else
      write(unit=unit_figplot,fmt='(a,i0,a)') '2 1 0 1 ', pen_color, &
        ' 7 50 -1 -1 0.000 0 0 -1 1 0 2'
      write(unit=unit_figplot,fmt='(a,f6.2,1x,f6.2)') '1 1 1.00 ', &
        MMFIGUNITS*length*0.08_dp, MMFIGUNITS*length*0.16_dp
      write(unit=unit_figplot,fmt='(4(1x,i0))') nint(MMFIGUNITS*x)
    end if

  end subroutine fig_vector


! plot a filled circle

  subroutine fig_filledcircle ( pen_color, x )

    integer, intent(in) :: pen_color

!   x: coordinates in mm
!   x(:,1) is center
!   x(:,2) is on circle
    real(dp), intent(in) :: x(2,2)

    real(dp) :: radius

    radius = sqrt( sum((x(:,2)-x(:,1))**2) )

    write ( unit=unit_figplot, fmt='(a,2(1x,i0),a,8(1x,i0))') '1 3 0 1 ', &
      pen_color, pen_color, ' 49 -1 20 0.000 1 0.0000 ', &
      nint(MMFIGUNITS*x(:,1)), nint(MMFIGUNITS*radius), &
      nint(MMFIGUNITS*radius), nint(MMFIGUNITS*x)

  end subroutine fig_filledcircle

! write a text

  subroutine fig_text ( fontsize, x, y, text, justification )

!   font size in points
    integer, intent(in) :: fontsize

!   x, y coordinates in mm
    real(dp), intent(in) :: x, y

!   the text to be plotted
    character (len=*), intent(in) :: text

!   justification:
!     0  Left justified (default)
!     1  Center justified
!     2  Right justified
    integer, optional, intent(in) :: justification


    integer :: height, length, sub_type

    if ( present(justification) ) then
      sub_type = justification
    else
      sub_type = 0
    end if

!   these are very rough estimates of the height and length.
!   read and save the resulting file with xfig to get better values.

    height = 13 * fontsize
    length = fontsize * ( 7 * (len(text)-1) + 10 )

    write(unit=unit_figplot,fmt='(2(a,i0),a,4(i0,1x),2a)') &
      '4 ', sub_type, ' 0 50 -1 0 ', fontsize, ' 0.0000 0 ', &
      height, length, nint(MMFIGUNITS*x), nint(MMFIGUNITS*y), text, '\001'

  end subroutine fig_text


! write colors

  subroutine fig_colors ( ncolors )

!   number of colors that are defined in the fig file
!   note that ncolors on input it is just an estimate.
!   output = 4*input/4 + 1
!   thus only for ncolors=5, 9, 13, 17, 21, ... the output is not changed.
    integer, intent(inout) :: ncolors


    integer :: nc, i, red, green, blue
    character (len=7) :: ch


!   split the range into four subranges

    nc = ncolors / 4

!   first nc colors: from blue to cyan

    do i = 1, nc
      red = 0
      green = 255 * ( i - 1 ) / nc
      blue = 255
      write (ch,'(a,3z2.2)') '#', red, green, blue
      write(unit=unit_figplot,fmt='(a,i0,1x,a,i0)') '0 ', i+31, ch
    end do

!   second nc colors: from cyan to green

    do i = 1, nc
      red = 0
      green = 255
      blue = 255 - 255 * ( i - 1 ) / nc
      write (ch,'(a,3z2.2)') '#', red, green, blue
      write(unit=unit_figplot,fmt='(a,i0,1x,a,i0)') '0 ', i+nc+31, ch
    end do

!   second nc colors: from green to yellow

    do i = 1, nc
      red = 255 * ( i - 1 ) / nc
      green = 255
      blue = 0
      write (ch,'(a,3z2.2)') '#', red, green, blue
      write(unit=unit_figplot,fmt='(a,i0,1x,a,i0)') '0 ', i+2*nc+31, ch
    end do

!   last nc+1 colors: from yellow to red

    do i = 1, nc + 1
      red = 255
      green = 255 - 255 * ( i - 1 ) / nc
      blue = 0
      write (ch,'(a,3z2.2)') '#', red, green, blue
      write(unit=unit_figplot,fmt='(a,i0,1x,a,i0)') '0 ', i+3*nc+31, ch
    end do

    ncolors = 4 * nc + 1

  end subroutine fig_colors

end module figplot_m
