// Run with gmsh -3 mesh8.geo

 ox =  -2.00000000000000;
 oy =  -2.00000000000000;
 oz =  -2.00000000000000;
 lx =   4.00000000000000;
 ly =   4.00000000000000;
 lz =   4.00000000000000;
 dx =   0.80000000000000;

Mesh.ElementOrder = 1; // First-order elements
Mesh.Optimize = 1; // Optimize to improve the quality of tetrahedral elements

 Include "mesh_hexa.igo";
