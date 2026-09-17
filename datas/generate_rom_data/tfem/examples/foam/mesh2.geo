// Cubic domain with a sphere in the center
// Run with gmsh -3 mesh2.geo

ox = 0.0;
oy = 0.0;
oz = 0.0;
lx = 4.0;
ly = 4.0;
lz = 4.0;
dx_box = 0.3;
xp = 2.0;
yp = 2.0;
zp = 2.0;
rp = 1.0;
//Nelem_part = 80;
Nelem_part = 40;

//Mesh.Algorithm = 2; // Automatic (default)
//Mesh.Algorithm3D = 2; // Delauney (default)
Mesh.ElementOrder = 2; // Second-order elements
Mesh.Optimize = 1; // Optimize to improve the quality of tetrahedral elements
//Mesh.HighOrderOptimize = 1; // Higher-order optimize
//Mesh.OptimizeNetgen = 1; // Optimize using Netgen

Include "cube_with_sphere.igo";
