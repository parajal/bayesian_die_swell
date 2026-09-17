// Single sphere surface (quadratic elements)
// Run with gmsh -2 mesh2q.geo

nobj = 1;
N = 25; // number of elements on the equator (approximately)
xp[1] = 0.00000000000000;
yp[1] = 0.00000000000000;
zp[1] = 0.00000000000000;
rp[1] = 1.00000000000000;
dx_part = 2*Pi*rp[1]/N;

Mesh.Algorithm = 6; // Frontal for surfaces
Mesh.ElementOrder = 2; // Second-order elements

Include "particles_3D.igo";
