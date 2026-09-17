// Run with gmsh -2 mesh9.geo

ox =  -2.00000000000000;
oy =  -2.00000000000000;
lx =   4.00000000000000;
ly =   4.00000000000000;
dx =   0.10000000000000;

// Mesh.Algorithm = 6; // Frontal 
Mesh.ElementOrder = 1; // First-order elements

Include "mesh_quad.igo";
