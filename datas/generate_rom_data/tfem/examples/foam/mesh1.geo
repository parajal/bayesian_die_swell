// Square domain with a circle in the center
// Run with gmsh -2 mesh1.geo

ox = 0.00000000000000;
oy = 0.00000000000000;
lx = 4.00000000000000;
ly = 4.00000000000000;
dx_box = 0.3;
xp =   2.50000000000000;
yp =   2.50000000000000;
rp =   1.00000000000000;
Nelem_part = 80;

// Mesh.Algorithm = 6; // Frontal 
Mesh.ElementOrder = 2; // Second-order elements

Include "rectangle_with_circle.igo";
