// Run with gmsh -2 mesh1.geo
// See tube_with_sphere1.igo for more info.

// constants

L1 = 24;
L2 = 6;
R1 = 1;
R2 = 2;

numelc1 = 30;
numelc2 = 12;
numelc3 = 12;
numelc4 = 12;
numelc8 = 5;

factor1 = 8;
factor2 = 4;
factor3 = 1.5;
factor4 = 8;
factor5 = 4;
factor6 = 10;

quad = 1; // quad = 0 (triangles) or 1 (quads)

Mesh.ElementOrder = 2; // Second-order elements

Include "tube_with_sphere1.igo";
