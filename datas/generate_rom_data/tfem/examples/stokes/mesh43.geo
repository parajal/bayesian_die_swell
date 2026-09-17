// Cubic domain with a sphere in the center
// Run with gmsh -3 mesh43.geo

ox = -150.0;
oy = -150.0;
oz = -150.0;
lx = 300.0;
ly = 300.0;
lz = 300.0;
dx_box = 30.0;
xp = 0.0;
yp = 0.0;
zp = 0.0;
rp = 1.0;
Nelem_part = 20;

//Mesh.Algorithm = 2; // Automatic (default)
//Mesh.Algorithm3D = 2; // Delauney (default)
Mesh.ElementOrder = 2; // Second-order elements
Mesh.Optimize = 1; // Optimize to improve the quality of tetrahedral elements
//Mesh.HighOrderOptimize = 1; // Higher-order optimize
//Mesh.OptimizeNetgen = 1; // Optimize using Netgen

Include "cube_with_sphere.igo";
