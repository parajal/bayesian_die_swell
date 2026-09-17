
! Copyright (C) 2007-2025 Martien A. Hulsen
! Permission to copy or distribute this software or documentation
! in hard copy or soft copy granted only by written license
! obtained from Martien A. Hulsen.
! All rights reserved. No part of this publication may be reproduced,
! stored in a retrieval system ( e.g., in memory, disk, or core)
! or be transmitted by any means, electronic, mechanical, photocopy,
! recording, or otherwise, without written permission from the
! publisher.


! Define user changeable limits and defaults set for tfem

! Changing the limits/defaults that are constants (parameters) requires a
! recompile of tfem core + all addons used.


module limits_m

  use kind_defs_m

  implicit none

!---------------------------------------------------------

! limits build in for generating or reading a mesh

  integer :: &
    MAXPOINTS   = 1000, &
    MAXCURVES   = 100,  &
    MAXOBJECTS  = 200,  &
    MAXSURFACES = 100,  &
    MAXVOLUMES  = 10,   &
    MAXNODESETS = 100,  &
    MAXELEMENTSETS = 100,  &
    MAXBLOCKS   = 1000

!---------------------------------------------------------

! maximum number of iterations for finding the factor in a distribution of
! elements on an interval (routine compute_factor).
  integer :: MAXNUMITER  = 100

!---------------------------------------------------------

! limits/accuracies build in for finding intersections of objects with the mesh

! maximum number of Newton iterations for finding the reference coordinates
! in the intersected element of a point in an object
  integer :: MAXITER = 100

! stop on 'no convergence' (MAXITER or xi > MAXVALXI) or continue searching.
  logical :: STOP_ON_NOCONV = .false.

! give a warning and continue searching on 'no convergence' (MAXITER or xi >
! MAXVALXI)
  logical :: WARN_ON_NOCONV = .false.

! perform a two-pass search. First pass always with a value of EPSELEMENT.
! If the reference coordinates are still not found and
! SEARCH_BLOCKS_PASSTWO=.true., a second pass is performed found with a
! value of EPSELEMENT2. This can be useful if points are only
! slightly outside the domain due to curved boundaries.
  logical :: SEARCH_BLOCKS_PASSTWO = .false.

! accuracy for the Newton iterations for finding the reference coordinates
! in the intersected element of a point in an object
  real(dp) :: EPSITER = 1e-10_dp

! accuracy for determining intersection of an element. The reference element
! is made slightly larger by EPSELEMENT to have a slight overlap between
! elements
  real(dp) :: EPSELEMENT = 1e-10_dp

! accuracy for determining intersection of an element. The reference element
! is made slightly larger by EPSELEMENT2 to have a slight overlap between
! elements. This is for pass two (see above).
  real(dp) :: EPSELEMENT2 = 1e-10_dp

! maximum value for xir allowed for in the iteration process
  real(dp), parameter :: MAXVALXI = 100._dp

!---------------------------------------------------------

! limits/accuracies build in for blocks

! default number of blocks (in each space direction) in fill_mesh_parts

  integer :: NUMBLOCKSX = 2

! make box around elements 1+EPSBOUND larger

  real(dp) :: EPSBOUND = 0.1_dp

!---------------------------------------------------------

! limits/accuracies build in for nodblocks

! create a box around a node of size 2*EPSNODENOX*(xmax-xmin), with
! xmax, xmin the upper right and the lower left corner of the full domain.
! The box extends EPSNODENOX*(xmax-xmin) to the positive and negative direction
! around a node, with the node at the center.

  real(dp) :: EPSNODEBOX = 1e-10_dp

!---------------------------------------------------------

! limits/accuracies build in for
!  1) removal of double nodes
!  2) check whether nodes are the same (renumber_nodes_to_match_coor)

! nodes overlap if the distance between the nodes is smaller than
! EPSNODEOVERLAP*(xmax-xmin), with xmax, xmin the upper right and the
! lower left corner of the full domain.

  real(dp) :: EPSNODEOVERLAP = 1e-10_dp

!---------------------------------------------------------

! limits build in for the coefficients

  integer :: &
    MAXFUNCTIONS = 10  ! default maximum number of pointers to functions

!---------------------------------------------------------

! limits build in for the problem definition

  integer ::      &
    MAXNUMESSPOINTS       = 50,  &
    MAXNUMESSCURVES       = 50,  &
    MAXNUMESSSURFACES     = 50,  &
    MAXNUMESSELEMENTS     = 10,  &
    MAXNUMESSGROUPS       = 5,   &
    MAXNUMESSNODESETS     = 50,  &
    MAXNUMCONSTRAINTS     = 100, &
    MAXNUMCONNECTIONS     = 10,  &
    MAXNUMTRANSFORMATIONS = 10,  &
    MAXNUMDEPENDENCIES    = 10,  &
    MAXNUMWARNINGS        = 5

!---------------------------------------------------------

! limits/defaults in add-on figplot

  integer :: MAXPOLYLINEPOINTS = 10

!---------------------------------------------------------

! limits/defaults in add-on io-utils

  integer ::  &
    MAXGROUPS  = 100,    & ! used in gbt_utils and marc_utils
    MAXWORK    = 100000, & ! used in marc_utils
    MAXDEFINES = 100       ! used in marc_utils

! give a warning if there is no data to output in some nodes
  logical :: WARN_ON_NO_DATA_IN_NODES = .true.

!---------------------------------------------------------

! limits/defaults in add-on subdivide.

! cut-off for zerolevelset, so if |phi|<=ZEROLEVELSET it is set to phi=0.
  real(dp) :: ZEROLEVELSET = 0

! Use six tetrahedrons instead of five for subdividing a hexahedron.
! It leads to a continuous zero levelset surface (on eltree level or for
! structured meshes), but is more expensive.
  logical :: SUBDIVIDE_HEX6 = .false.

!---------------------------------------------------------

! limits/defaults in fill_mesh_parts_sidelem

! Do not fill the sidelem array if
!  1) line elements for ndim=2,3 or
!  2) triangular or quadrilateral elements for ndim=3,
! because of the possibility to have multiple connected elements.
! NOTE: this also has the side effect that plotting with figplot does not
!       work correctly, such as plotting the full mesh in data plots.
  logical :: TEST_FOR_MULTIPLE_SIDELEM = .true.

! Give a warning if the sideelem array is indeed not filled.
  logical :: WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = .true.

! Give a warning for filling the sidelem array in the mesh and the
! meshnumelem/meshelem arrays in the surfaces, if prism or pyramid shaped
! elements are present in the mesh. These shapes are not supported for
! filling these arrays.
  logical :: WARN_ON_TEST_FOR_VARYING_SIDES = .true.

!---------------------------------------------------------

! limits/defaults in fill_element_vtk

! Use original vtk elements where possible
! This affects the following shapes:
!                         .true.                   .false.
!    elshape=7      seven-node triangle       use six-node triangle
!    elshape=6,34   nine-node quadrilateral   use four four-node quads
!    elshape=14,36  27-node hexahedron        use eight eight-node hexas
  logical :: USE_VTK_ELEMENTS = .true.

!---------------------------------------------------------

! limits/defaults in fill_mesh_parts_geometries_edges

! Curves are defined as aliases for element edges and elements on the curve
! are expected to be connected to nodes on an element edge. If not, a warning
! is given. In addition to the warning, the curve elements that violate this
! definition can be printed by setting the variable below to .true.
  logical :: PRINT_CURVE_ELEMENTS_NOT_ON_EDGE = .false.

! The warning can be completely shut off by setting the variable below to
! .false. Use only if you know what you are doing!
  logical :: WARN_ON_EDGES_AND_FACES = .true.

!---------------------------------------------------------

! limits/defaults in fill_mesh_parts_geometries_faces

! Surfaces are defined as aliases for element faces and elements on the surface
! are expected to be connected to nodes on an element face. If not, a warning
! is given. In addition to the warning, the surface elements that violate this
! definition can be printed by setting the variable below to .true.
  logical :: PRINT_SURFACE_ELEMENTS_NOT_ON_FACE = .false.

! The warning can be completely shut off by setting the variable
! WARN_ON_EDGES_AND_FACES to .false. (see above).
!---------------------------------------------------------

! limits/defaults in read_mesh_gmsh (module gmsh_utils_m)

! Nodes should be connected to internal elements. If not, a warning is given.
! In addition to the warning, the nodes involved can be printed by setting
! the variable below to .true.
  logical :: PRINT_NODES_NOT_CONNECTED = .true.

! The warning can be completely shut off by setting the variable below to
! .false. Use only if you know what you are doing!
  logical :: WARN_ON_NODES_NOT_CONNECTED = .true.

!---------------------------------------------------------

! limits/defaults in module vtk_utils_m

! Write double data format vtk files.
  logical :: WRITE_DOUBLE_VTK = .false.

! Write binary vtk files.
  logical :: WRITE_BINARY_VTK = .true.

!---------------------------------------------------------

! limits/defaults in routine mesh_change_geometries_remove_double_nodes of the
! module meshgen_extra_m

! remove double points (=points refering to the same nodes of the mesh)
  logical :: REMOVE_DOUBLE_POINTS = .true.

!---------------------------------------------------------

! limits in shape functions

! set detF = abs(detF)
  logical :: SET_ABSOLUTE_DETF = .true.

! check if detF is positive and stop if not.
! Note: only activated if SET_ABSOLUTE_DETF = .false.
  logical :: CHECK_POSITIVE_DETF = .false.

! cut-off value of |1-zeta| in the rational shape function of a pyramid
  real(dp) :: EPS_PYRAMID = 1e-30_dp

!---------------------------------------------------------

! limits/defaults in specifying the Gauss integration rule

! Enforce that the Gauss integration rule (intrule,intrule2) in the standard
! elements (elements defined in the addons) must be specified by the polynomial
! order of integration.
  logical :: SET_GAUSS_BY_ORDER = .false.

! Gauss rule chooser. Within a given type of integration (inttype), this
! can be used to choose between multiple available rules.
! Note, for GAUSS_RULE=1 the default rule is always chosen.
! Available:
!  inttype=3
!    hexahedron, GAUSS_RULE=2: rule from gauss_standard_numeric2_m, with
!                intrule = order of integration
  integer :: GAUSS_RULE = 1

!---------------------------------------------------------

! limits/defaults in viscoelastic models

! Use the conformation model routines and and transform to the b-formulation
! when using contravariant deformation models
  logical :: USE_CONFORMATION_MODEL = .false.

! Use the EGP model where the stress tensor is used in the evolution equation
! instead of the conformation only equation where the stress in substituted.
  logical :: USE_EGP_STRESS_TENSOR_FORM = .false.

! Use the adapted lambda in the determinant stabilization of the EGP. Otherwise
! if set to .false. the stabilization is performed with the original lambda
! constant.
  logical :: USE_ADAP_LAMBDA_FOR_EGP_DET_STAB = .true.

!---------------------------------------------------------

! limits/defaults in using the time integration schemes

! This applies to using higher-order (>=2) methods with evaluations of the
! convection term at different times within a single time step in combination
! with a moving mesh (ALE). In that case, it is required that the change in
! coordinates needs to be taken into account for the nabla operator. If not,
! the time integration becomes a first-order scheme in theory. By setting
!    BLOCK_IMPROPER_ALE = .false.
! you can bypass the blockage and still use the scheme as implemented, i.e.
! assuming a constant nabla operator within the time step.
  logical :: BLOCK_IMPROPER_ALE = .true.

! Implementation of the Jacobian for viscoelastic models that depend on J
! combined with Crank-Nicolson or theta time-integration is incomplete.
! Set to .false. to suppress the warning message.
  logical :: WARN_ON_INCOMPLETE_JACOBIAN = .true.

end module limits_m
