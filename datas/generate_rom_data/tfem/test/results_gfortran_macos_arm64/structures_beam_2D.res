
Warning fill_mesh_parts_sidelem: mesh contains 
  1) line elements for ndim=2,3 or 
  2) triangular or quadrilateral elements for ndim=3 
 and therefore mesh%sidelem is not computed because 
 multiple connected elements are possible. This can have 
 some unexpected side effects, such as mesh plots in 
 data plots with figplot. Set (use limits_m module) 
   WARN_ON_TEST_FOR_MULTIPLE_SIDELEM = .false. 
 to stop this warning from appearing or 
   TEST_FOR_MULTIPLE_SIDELEM = .false. 
 to skip the test (and compute mesh%sidelem anyway), 
 if your mesh does not have multiple connected elements.

Mesh info:

 Space dimension (ndim)                = 2
 Number of nodes (nnodes)              = 26
 Number of elements (nelem)            = 25
 Number of element groups (nelgrp)     = 1
 Number of blend meshes (nblend)       = 0
 Number of points (npoints)            = 3
 Number of curves (ncurves)            = 0
 Number of surfaces (nsurfaces)        = 0
 Number of volumes (nvolumes)          = 0
 Number of objects (nobjects)          = 0
 Number of nodesets (nnodesets)        = 0
 Number of nelementsets (nelementsets) = 0
 Number of blocks (nblocks)            = 4
 Number of gluepoints (ngluepoints)    = 0
 meshparts = T renumber = F sblocks = F nodblocks = F
 Number of nodes per element in each group:
  2
 maxnodnumnod = 2
 Size of mesh%nodnod = 50
 Element shape = 1

Problem info:

 Number of vectors (nvec)                                  = 2
 Number of physical quantities (nphysq)                    = 0
 Number of physqshifted (nphysqshifted)                    = 0
 Number of constraints (numconstraints)                    = 0
 Number of degrees of freedom (numdegfd)                   = 78
 Number of nodal degrees of freedom (numnodaldegfd)        = 78
 Number of constraint degrees of freedom (numconstrdegfd)  = 0
 Number of unknown degrees of freedom (numundegfd)         = 72
 Number of essential degrees of freedom (numessdegfd)      = 6
 Number of inactive groups (numinactivegroups)             = 0
 Number of degrees for each vector (nodalwise): 
  26 52
 Number of degrees for each vector (elementwise): 
  50 100
 Element definition for the sysvector: 
  3 3
 Element definition for the vectors: 
  Vector 1
  1 1
  Vector 2
  2 2

 external forces     external moments    maxval abs(sysvector)   3.0000000000000000     
 reaction forces     reaction moments    maxval abs(sysvector)   2.3843989011290652     
 maxval abs(sol)   2.1054426394829338     
 distributed forcemaxval abs(dvec)   1.0000000000000000     
 transverse displacementmaxval abs(dvec)   1.4990198840277880     
 slopemaxval abs(dvec)   2.1054426394829338     
 second derivativemaxval abs(dvec)   2.3852322344623955     
 bending momentmaxval abs(dvec)   2.3852322344623955     
 shear forcemaxval abs(dvec)   1.1745791899612643     
 axial displacementmaxval abs(dvec)   1.0068313313543245     
 axial strainmaxval abs(dvec)   1.0068313313543265     
 axial forcemaxval abs(dvec)   1.0068313313543265     
