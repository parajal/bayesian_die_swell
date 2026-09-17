// Rectangular domain with a circle in the center (symmetric version: half of
// the domain only).
// Run with gmsh -2 mesh44.geo

 ox =  -200.00000000000000;
 oy =   0.00000000000000;
 lx =   400.00000000000000;
 ly =   200.0000000000000;
 xp =  0.0000000000000;
 rp =   1.00000000000000;
 dx_box =   10.000000000000;
 Nelem_part = 20;

// Mesh.Algorithm = 6; // Frontal 
Mesh.ElementOrder = 2; // Second-order elements

Include "rectangle_with_circle_sym.igo";
