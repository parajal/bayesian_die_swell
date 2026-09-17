
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

 Space dimension (ndim)                = 3
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
 Number of blocks (nblocks)            = 8
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
 Number of degrees of freedom (numdegfd)                   = 156
 Number of nodal degrees of freedom (numnodaldegfd)        = 156
 Number of constraint degrees of freedom (numconstrdegfd)  = 0
 Number of unknown degrees of freedom (numundegfd)         = 144
 Number of essential degrees of freedom (numessdegfd)      = 12
 Number of inactive groups (numinactivegroups)             = 0
 Number of degrees for each vector (nodalwise): 
  26 78
 Number of degrees for each vector (elementwise): 
  50 150
 Element definition for the sysvector: 
  6 6
 Element definition for the vectors: 
  Vector 1
  1 1
  Vector 2
  3 3

 external forces     external moments    maxval abs(sysvector)  10.0000000000000000
 reaction forces     reaction moments    maxval abs(sysvector)  10.7137551066969099
 maxval abs(sol)   2.5847358881968940
 distributed forcemaxval abs(dvec)   0.0000000000000000
 transverse displacement vmaxval abs(dvec)   1.6681963542389293
 slope v'maxval abs(dvec)   2.3314516580774525
 second derivative v''maxval abs(dvec)   3.2526553631791164
 bending moment Mzmaxval abs(dvec)   3.2526553631791164
 shear force Vmaxval abs(dvec)   1.8424074102049628
 axial displacementmaxval abs(dvec)   1.0399261289770430
 axial strainmaxval abs(dvec)   1.0399261289770447
 axial forcemaxval abs(dvec)   1.0399261289770447
 axial (torsional) rotation thetamaxval abs(dvec)   3.1939639795702037
 torsional strain dtheta/dsmaxval abs(dvec)   2.5847358881968994
 torsional moment Mxmaxval abs(dvec)   2.5847358881968994
 transverse displacement wmaxval abs(dvec)   0.9623306159764959
 slope w'maxval abs(dvec)   1.9322112894424310
 second derivative w''maxval abs(dvec)   4.1141467433354748
 bending moment Mymaxval abs(dvec)   8.2282934866709496
 shear force Wmaxval abs(dvec)  10.7137551066982954
