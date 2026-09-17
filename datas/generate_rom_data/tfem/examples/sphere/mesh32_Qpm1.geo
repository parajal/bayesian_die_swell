// Run with gmsh -2 mesh32_Qpm1.geo
// See tube_with_sphere1.igo for more info.

// constants

L1 = 24;
L2 = 6;
R1 = 1;
R2 = 2;

numelc1 = 15;
numelc2 = 6;
numelc3 = 6;
numelc4 = 6;
numelc8 = 3;

factor1 = 8;
factor2 = 4;
factor3 = 1.5;
factor4 = 4;
factor5 = 4;
factor6 = 10;

quad = 1; // quad = 0 (triangles) or 1 (quads)

Mesh.ElementOrder = 3; // Third-order elements

Include "tube_with_sphere1.igo";
